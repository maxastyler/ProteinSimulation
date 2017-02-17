#!/usr/bin/python2
#For ranaspumin:
#   don't create dihedral bonds for residues <= 16 (this is the floppy tail)
#   Add a disulfide bond between the cysteines C68-C86
#For CST3 :
#   Maybe make residue 67 - 75 floppy?
#   Add disulfide bond to C62-C70 and C84-C104

import sys
import math
import numpy as NP

# define atoms class
class Atom:
    pass
class Bond:
    pass
class Angle:
    pass
class Dihedral:
    pass

def splitp(l):
    f,i=float,int
    return [l[0:6],i(l[6:11]),l[12:16],l[16],l[17:20],l[21],i(l[22:26]),l[26],f(l[30:38]),f(l[38:46]),f(l[46:54]),f(l[54:60]),f(l[60:66])]
def writeatom(a):
    return '%6s%5d %4s %3s  %4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n'%('ATOM  ',a.id,a.type,a.res,a.rid,a.pos[0],a.pos[1],a.pos[2],1.0,1.0)
def writeatom2(a):
    return '%6s%5d %4s %3s  %4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n'%('ATOM  ',a.id,a.type,a.res,a.rid,a.pos[0],a.pos[1],a.pos[2],1.0,a.sigma)
    '''
Splits a line from a pdb file in the correct fields, and convert Angstrom to nm.
 1 -  6        Record name     "ATOM  "                                            
 7 - 11        Integer         Atom serial number.                   
13 - 16        Atom            Atom name.                            
17             Character       Alternate location indicator.         
18 - 20        Residue name    Residue name.                         
22             Character       Chain identifier.                     
23 - 26        Integer         Residue sequence number.              
27             AChar           Code for insertion of residues.       
31 - 38        Real(8.3)       Orthogonal coordinates for X in Angstroms.                       
39 - 46        Real(8.3)       Orthogonal coordinates for Y in Angstroms.                            
47 - 54        Real(8.3)       Orthogonal coordinates for Z in Angstroms.                            
55 - 60        Real(6.2)       Occupancy.                            
61 - 66        Real(6.2)       Temperature factor (Default = 0.0).                   
73 - 76        LString(4)      Segment identifier, left-justified.   
77 - 78        LString(2)      Element symbol, right-justified.      
79 - 80        LString(2)      Charge on the atom.
    '''

# define element masses
global mass
mass={}
mass['C']=12.0107
mass['F']=18.9984
mass['H']=1.00794
mass['N']=14.0067
mass['O']=15.9994
mass['S']=32.065

# center of mass
def com(atoms):
  dim=len(atoms[0].pos)
  cm=NP.array([0. for i in range(dim)])
  M=0.
  for a in atoms:
    cm+=mass[a.element]*a.pos
    M+=mass[a.element]
  return cm/M
def com_nomass(atoms):
  dim=len(atoms[0].pos)
  cm=NP.array([0. for i in range(dim)])
  M=0.
  for a in atoms:
    cm+=1.0*a.pos
    M+=1.0
  return cm/M

# distance
def norm(v):
  return math.sqrt(NP.dot(v,v))

# angle
def angle(p1,p2,p3):
  r21=p1-p2
  r23=p3-p2
  return 180.0*math.acos(NP.dot(r21,r23)/(norm(r21)*norm(r23)))/math.pi

# dihedral
def get_orth(vec,vref):
  temp = vec-vref*NP.dot(vec,vref)
  return temp/norm(temp)
def dihedral(a,b,c,d):
  bc = (c-b)/norm(c-b)
  bao = get_orth(a-b,bc)
  cdo = get_orth(d-c,bc)
  psi = math.acos(NP.dot(bao,cdo))*180.0/math.pi
  sign = NP.dot(NP.cross(bao,cdo),bc)
  if(sign<0): psi=-psi
  return NP.round(psi+180.0)

catype='CA'
sctype='SC'
def coarsegrain(atoms,aacontacts):
  cgatoms=[]
  nativepairs={} # native non-bonded interactions for the cgatoms
  ncontacts={}   # number of contacts in each residue-residue native pair
  i=0
  ncg=0
  while(i<len(atoms)):
    rid=atoms[i].rid
    res=atoms[i].res
    caflag=0
    scflag=0
    caatoms=[]
    scatoms=[]
    while(atoms[i].rid==rid):
      type=atoms[i].type
      if( type=='CA' or type=='N' or type=='C' or type=='O' ):
        caatoms.append(atoms[i])
        if(type=='CA'): caflag=1; capos=atoms[i].pos
      else:
        scatoms.append(atoms[i])
        scflag=1
      i+=1
      if(i>=len(atoms)): break
    # make ca bead
    if(caflag==0): print '# error: CA bead not found'; sys.exit()
    ncg+=1
    atom=Atom()
    atom.id=ncg
    atom.type=catype
    atom.res=res
    atom.rid=rid
    atom.pos=capos
    atom.nheavy=4
    atom.atoms=caatoms
    cgatoms.append(atom)
    if(scflag==1): # make side-chain bead
      scpos=com(scatoms)
      ncg+=1
      atom=Atom()
      atom.id=ncg
      atom.type=sctype
      atom.res=res
      atom.rid=rid
      atom.pos=scpos
      atom.nheavy=sum([1 for a in scatoms])
      atom.atoms=scatoms
      cgatoms.append(atom)
  for a in cgatoms: a.sigma=4.0*(a.nheavy/4.)**(1./3.)
  ncg=len(cgatoms)
  for i in range(ncg):
    ai=cgatoms[i]
    for j in range(i+1,ncg):
      aj=cgatoms[j]
      if(ai.rid>aj.rid): print '# error: cgatoms not well ordered'; sys.exit()
      for aii in ai.atoms:
        for ajj in aj.atoms:
          if( (aii.id,ajj.id) in aacontacts ):
            if((aj.rid-ai.rid)<=3): #print '# error: (aj.rid-ai.rid)<=3'; sys.exit()
                print "# Warning: (aj.rid-ai.rid)<=3"
                print "At aj.rid: " + str(aj.rid) + " and ai.rid: " + str(ai.rid)
            if(ai.id>=aj.id): print 'error: ai.id>=aj.id'; sys.exit()
            if( (ai.id,aj.id) not in nativepairs ):
              pair=Bond()
              pair.a1=ai
              pair.a2=aj
              pair.name='%s%d-%s%d'%(ai.type,ai.rid,aj.type,aj.rid)
              pair.ncontacts=1
              pair.rids=(ai.rid,aj.rid)
              nativepairs[(ai.id,aj.id)]=pair
            else: nativepairs[(ai.id,aj.id)].ncontacts+=1
            if( (ai.rid,aj.rid) not in ncontacts ): ncontacts[(ai.rid,aj.rid)]=1
            else: ncontacts[(ai.rid,aj.rid)]+=1
  # normalisation
  for i in nativepairs:
    p=nativepairs[i]
    p.efactor=1.0*p.ncontacts/ncontacts[p.rids]
  return cgatoms,nativepairs

# Open and read pdb file
try: infile = open(sys.argv[1],'r').readlines(); pdbfile=sys.argv[1]
except: print '# usage: ./readpdb.py example.pdb aa_contacts.dat'; sys.exit()
print '# reading ',pdbfile,'..'
atoms=[]
for l in infile:
  if(l[:6]=='CRYST1'): box=[float(l[6:15]),float(l[15:24]),float(l[24:33])]
  elif(l[:4]=='ATOM'):
    tokens=splitp(l)
    element=l[13]
    if(element not in mass): print '# element not recognised'; sys.exit()
    if(element=='H'): print '# error: remove hydrogens'; sys.exit()
    atom=Atom()
    atom.element=element
    atom.id=tokens[1]
    atom.type=tokens[2].split()[0]
    atom.res=tokens[4].split()[0]
    atom.rid=tokens[6]
    atom.pos=NP.array(tokens[8:11])
    atoms.append( atom )
print '# %d atoms found'%len(atoms)

# read native contacts
try: infile = open(sys.argv[2],'r').readlines()
except: print '# usage: ./readpdb.py example.pdb aa_contacts.dat'; sys.exit()
print '# reading ',sys.argv[2],'..'
aacontacts={}
for l in infile:
  pair=l.split()
  ida=int(pair[0])
  idb=int(pair[1])
  aacontacts[(ida,idb)]=1
  aacontacts[(idb,ida)]=1

# coarse grain atomic positions
print '# coarse graining..'
cgatoms,nativepairs=coarsegrain(atoms,aacontacts)
distances=NP.array([ [ norm(a2.pos-a1.pos) for a1 in cgatoms] for a2 in cgatoms ])
outpdb=open('cg_%s'%(pdbfile),'w')
for a in cgatoms: outpdb.write(writeatom2(a))
ncg=len(cgatoms)
nca=0; nsc=0
for a in cgatoms:
  if(a.type==catype): nca+=1
  elif(a.type==sctype): nsc+=1
  else: print '# atom type not recognised'; sys.exit()
print '# %d cg beads created'%ncg
print '# %d c_alpha beads'%nca
print '# %d side-chain beads'%nsc
print '# %d native contacts'%len(nativepairs)
# set center of mass to zero
center=com_nomass(cgatoms)
for a in cgatoms: a.pos=a.pos-center

# define all pairs
pairs=[]
for i in range(ncg):
  ai=cgatoms[i]
  for j in range(i,ncg):
    aj=cgatoms[j]
    if( (ai.id,aj.id) in nativepairs ):
      pairs.append(nativepairs[(ai.id,aj.id)])
    else:
      b=Bond()
      b.a1=ai
      b.a2=aj
      b.name='%s%d-%s%d'%(b.a1.type,b.a1.rid,b.a2.type,b.a2.rid)
      pairs.append(b)
# print out coarse grained contacts
outcontacts=open('cg_contacts.dat','w')
tmp=nativepairs.items()
tmp=sorted(tmp, key=lambda x:x[0][1])
tmp=sorted(tmp, key=lambda x:x[0][0])
outcontacts.write('# type1 resid1 id1 type2 resid2 id2 distance\n')
for pair in tmp:
  p=pair[1]
  p.eq=norm(p.a2.pos-p.a1.pos)
  outcontacts.write('%s\t%d\t%d\t%s\t%d\t%d\t%f\t%f\t%d\n'%(p.a1.type,p.a1.rid,p.a1.id,p.a2.type,p.a2.rid,p.a2.id,p.eq,p.efactor,p.ncontacts))

# define bonds
bonds=[]
angles=[]
dihedrals=[]
# make bonds
nbonds=0
bondmatrix=NP.array([ [ 0 for i in range(ncg)] for j in range(ncg) ])
for i in range(ncg):
  ai=cgatoms[i]
  typei=ai.type
  ridi=ai.rid
  for j in range(i+1,ncg):
    aj=cgatoms[j]
    typej=aj.type
    ridj=aj.rid
    # CA-CA
    if((ridj-ridi)==1 and typei==catype and typej==catype):
      bondmatrix[i,j]=1
      bondmatrix[j,i]=1
      b=Bond()
      nbonds+=1
      b.id=nbonds
      b.type=b.id
      b.a1=ai
      b.a2=aj
      b.name='%s%d-%s%d'%(b.a1.type,b.a1.rid,b.a2.type,b.a2.rid)
      bonds.append(b)
    # CA-SC
    if(ridi==ridj and typei==catype and typej==sctype):
      bondmatrix[i,j]=1
      bondmatrix[j,i]=1
      b=Bond()
      nbonds+=1
      b.id=nbonds
      b.type=b.id
      b.a1=ai
      b.a2=aj
      b.name='%s%d-%s%d'%(b.a1.type,b.a1.rid,b.a2.type,b.a2.rid)
      bonds.append(b)
# check
if(len(bonds)!=(nca-1+nsc)): print '# nbonds!=nca-1+nsc'; sys.exit()
else: print '#',len(bonds),'bonds created'
# angles
# CA-CA-CA
# CA-CA-SC
# SC-CA-CA
# search for all angles with atom i as centre
nangles=0
for i in range(ncg): # a2
  for j in range(ncg): # a1
    if(j==i): continue
    if(not bondmatrix[i,j]): continue
    for k in range(j+1,ncg): # a3
      if(k==i): continue
      if(not bondmatrix[i,k]): continue
      a=Angle()
      nangles+=1
      a.id=nangles
      a.type=a.id
      a.a1=cgatoms[j]
      a.a2=cgatoms[i]
      a.a3=cgatoms[k]
      a.name='%s%d-%s%d-%s%d'%(a.a1.type,a.a1.rid,a.a2.type,a.a2.rid,a.a3.type,a.a3.rid)
      angles.append(a)
nanglescheck=nca-2+2*nsc
if(cgatoms[1].type==sctype): nanglescheck-=1
if(cgatoms[ncg-1].type==sctype): nanglescheck-=1
if(len(angles)!=nanglescheck): print '# %d nangles incorrect'%nangles; sys.exit()
else: print '#',len(angles),'angles created'
# dihedrals
# CA-CA-CA-CA
# CA-CA-CA-SC
# SC-CA-CA-CA
# SC-CA-CA-SC
# search for dihedrals starting from the bond i-j
ndihedrals=0
#ridflex=16 # included
ridflex = (67,75) #inclusive
def inbetween(x, min_max): # inclusive
    if (x>=min_max[0] and x<=min_max[1]): return True
    else: return False
inflex = lambda x: inbetween(x, ridflex)
for b in bonds:
  i=b.a1.id-1
  j=b.a2.id-1
  ridi=b.a1.rid
  ridj=b.a2.rid
  for k in range(ncg):
    for l in range(k+1,ncg):
      ridk=cgatoms[k].rid
      ridl=cgatoms[l].rid
      if(l==i or l==j or k==i or k==j): continue
      #if(ridi<=ridflex or ridj<=ridflex or ridk<=ridflex or ridl<=ridflex): continue # this dihedral will be flexible
      if (inflex(ridi) or inflex(ridj) or inflex(ridk) or inflex(ridl)): continue
      # bond k-i-j-l
      if(bondmatrix[k,i] and bondmatrix[l,j]):
        # n=1
        d=Dihedral()
        ndihedrals+=1
        d.id=ndihedrals
        d.type=d.id
        d.a1=cgatoms[k]
        d.a2=cgatoms[i]
        d.a3=cgatoms[j]
        d.a4=cgatoms[l]
        d.n=1
        d.name='%s%d-%s%d-%s%d-%s%d-n%d'%(d.a1.type,d.a1.rid,d.a2.type,d.a2.rid,d.a3.type,d.a3.rid,d.a4.type,d.a4.rid,d.n)
        dihedrals.append(d)
        # n=3
        d=Dihedral()
        ndihedrals+=1
        d.id=ndihedrals
        d.type=d.id
        d.a1=cgatoms[k]
        d.a2=cgatoms[i]
        d.a3=cgatoms[j]
        d.a4=cgatoms[l]
        d.n=3
        d.name='%s%d-%s%d-%s%d-%s%d-n%d'%(d.a1.type,d.a1.rid,d.a2.type,d.a2.rid,d.a3.type,d.a3.rid,d.a4.type,d.a4.rid,d.n)
        dihedrals.append(d)
      # bond l-i-j-k
      if(bondmatrix[l,i] and bondmatrix[k,j]):
        # n=1
        d=Dihedral()
        ndihedrals+=1
        d.id=ndihedrals
        d.type=d.id
        d.a1=cgatoms[l]
        d.a2=cgatoms[i]
        d.a3=cgatoms[j]
        d.a4=cgatoms[k]
        d.n=1
        d.name='%s%d-%s%d-%s%d-%s%d-n%d'%(d.a1.type,d.a1.rid,d.a2.type,d.a2.rid,d.a3.type,d.a3.rid,d.a4.type,d.a4.rid,d.n)
        dihedrals.append(d)
        # n=3
        d=Dihedral()
        ndihedrals+=1
        d.id=ndihedrals
        d.type=d.id
        d.a1=cgatoms[l]
        d.a2=cgatoms[i]
        d.a3=cgatoms[j]
        d.a4=cgatoms[k]
        d.n=3
        d.name='%s%d-%s%d-%s%d-%s%d-n%d'%(d.a1.type,d.a1.rid,d.a2.type,d.a2.rid,d.a3.type,d.a3.rid,d.a4.type,d.a4.rid,d.n)
        dihedrals.append(d)
print '#',len(dihedrals),'dihedrals created'
# create disulfide bridge, after angles and dihedrals
#Create two disulfide bridges 
idcys1=62
idcys2=70
b=Bond()
nbonds+=1
b.id=nbonds
b.type=b.id
b.a1=cgatoms[idcys1-1]
b.a2=cgatoms[idcys2-1]
b.name='%s%d-%s%d'%(b.a1.type,b.a1.rid,b.a2.type,b.a2.rid)
bonds.append(b)

idcys1=84
idcys2=104
b=Bond()
nbonds+=1
b.id=nbonds
b.type=b.id
b.a1=cgatoms[idcys1-1]
b.a2=cgatoms[idcys2-1]
b.name='%s%d-%s%d'%(b.a1.type,b.a1.rid,b.a2.type,b.a2.rid)
bonds.append(b)

# define hamiltonian
eps=1.0*0.239 # energy given in kcal/mol
outparam=open('settings.param','w')
econtact=0.0
edihedral=0.0
for p in pairs:
  p.c10=0.0
  sigma=0.0
  type=' '
  epair=0.0
  if( (p.a1.id,p.a2.id) in nativepairs ):
    sigma=p.eq
    epair=p.efactor*eps
    p.c6=-2.0*epair*sigma**6
    p.c12=+1.0*epair*sigma**12
    econtact+=epair
    type='native'
  else:
    sigma=0.5*(p.a1.sigma+p.a2.sigma)
    epair=0.01*eps
    p.c6=0.0
    p.c12=epair*sigma**12
    type='non-native'
  outparam.write('pair_coeff %d %d %e %e %e # %s type=%s eps=%f sigma=%f\n'%(p.a1.id,p.a2.id,p.c6,p.c10,p.c12,p.name,type,epair,sigma))
for b in bonds:
  b.k=100.0*eps
  b.eq=norm(b.a2.pos-b.a1.pos)
  outparam.write('bond_coeff %d %f %f # %s\n'%(b.id,b.k,b.eq,b.name))
for a in angles:
  a.k=20.0*eps
  a.eq=angle(a.a1.pos,a.a2.pos,a.a3.pos)
  outparam.write('angle_coeff %d %f %f # %s\n'%(a.id,a.k,a.eq,a.name))
for d in dihedrals:
  k=0.0
  k=0.5*eps
  if(d.n==1):   d.k=1.0*k
  elif(d.n==3): d.k=0.5*k
  else: print '# n in dihedral invalid'; sys.exit()
  if(d.n==1): edihedral+=d.k
  d.eq=d.n*dihedral(d.a1.pos,d.a2.pos,d.a3.pos,d.a4.pos)
  outparam.write('dihedral_coeff %d %f %d %d %f # %s\n'%(d.id,d.k,d.n,d.eq,0.0,d.name))
print '# E_contact =',  econtact/eps, 'eps'
print '# E_dihedral =', edihedral/eps, 'eps'
print '# E_contact/E_dihedral =', econtact/edihedral

# print data file
for a in cgatoms: a.mol=1; a.ix=0; a.iy=0; a.iz=0; a.type=a.id
print '# printing data file..'
outdata=open('data.input','w')
outdata.write('LAMMPS data file\n\n')
outdata.write('%d atoms\n'%len(cgatoms))
if(len(bonds)>0): outdata.write('%d bonds\n'%len(bonds))
if(len(angles)>0): outdata.write('%d angles\n'%len(angles))
if(len(dihedrals)>0): outdata.write('%d dihedrals\n'%len(dihedrals))
outdata.write('\n')
outdata.write('%d atom types\n'%len(set(a.type for a in cgatoms)))
if(len(bonds)>0): outdata.write('%d bond types\n'%len(set(b.type for b in bonds)))
if(len(angles)>0): outdata.write('%d angle types\n'%len(set(a.type for a in angles)))
if(len(dihedrals)>0): outdata.write('%d dihedral types\n'%len(set(d.type for d in dihedrals)))
outdata.write('\n')
halfbox=200.0
outdata.write('%f %f xlo xhi\n'%(-halfbox, halfbox))
outdata.write('%f %f ylo yhi\n'%(-halfbox, halfbox))
outdata.write('%f %f zlo zhi\n'%(-halfbox, halfbox))
# print masses
outdata.write('\nMasses\n\n')
for type in sorted(set(a.type for a in cgatoms)): outdata.write('%d 1.0\n'%type)
# print atoms
outdata.write('\nAtoms\n\n')
for a in cgatoms: outdata.write('%d %d %d %f %f %f %d %d %d\n'%(a.id,a.mol,a.type,a.pos[0],a.pos[1],a.pos[2],a.ix,a.iy,a.iz))
# print bonds
if(len(bonds)>0):
  outdata.write('\nBonds\n\n')
  for b in bonds: outdata.write('%d %d %d %d # %s\n'%(b.id,b.type,b.a1.id,b.a2.id,b.name))
# print angles
if(len(angles)>0):
  outdata.write('\nAngles\n\n')
  for a in angles: outdata.write('%d %d %d %d %d # %s\n'%(a.id,a.type,a.a1.id,a.a2.id,a.a3.id,a.name))
# print dihedrals
if(len(dihedrals)>0):
  outdata.write('\nDihedrals\n\n')
  for d in dihedrals: outdata.write('%d %d %d %d %d %d # %s\n'%(d.id,d.type,d.a1.id,d.a2.id,d.a3.id,d.a4.id,d.name))

print '# done'
