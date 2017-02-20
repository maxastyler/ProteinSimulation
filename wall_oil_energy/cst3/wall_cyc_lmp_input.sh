pdbfile=$1
ewall=$2
awk '
  BEGIN{
   # free energies of partitioning to cyclohexane in kcal/mol from:
   # Wolfenden et. al.,
   # Comparing the polarities of the amino acids: 
   # side-chain distribution coefficients between the vapor phase, cyclohexane, 1-octanol, and neutral aqueous solution
   # Biochemistry, 1988
   H["ILE"]= 4.92 
   H["PHE"]= 2.98
   H["VAL"]= 4.04
   H["LEU"]= 4.92
   H["TRP"]= 2.33
   H["MET"]= 2.35
   H["ALA"]= 1.81
   H["GLY"]= 0.94
   H["CYS"]= 1.28
   H["TYR"]=-0.14
   H["PRO"]=H["VAL"]	# approximation
   H["THR"]=-2.57
   H["SER"]=-3.40
   H["HIS"]=-4.66
   H["GLU"]=-6.81
   H["ASN"]=-6.64
   H["GLN"]=-5.54
   H["ASP"]=-8.72
   H["LYS"]=-5.55
   H["ARG"]=-14.92
   Tfactor='$ewall' # ~ Tsimulation/Treal (e.g. 118.0/300.0)
   # for attraction:
   # V(r) = c9/r^9 - c3/r^3
   #      = eps [ 1/2 (sigma/r)^9 - 3/2 (sigma/r)^3 ]
   # F(r) = 9/2 eps/sigma [ (sigma/r)^10 - (sigma/r)^4 ]
   # c3 = 3/2 eps sigma^3
   # c9 = 1/2 eps sigma^9
   # for repulsion
   # V(r) = c9/r^9
   #      = eps (sigma/r)^9   # which means that at r=sigma the energy is eps
   # c3 = 0.0
   # c9 = eps sigma^9
   etot=0.0
 }
 {
  id=$2
  type=$3
  name=$4
  if( type=="SC" || name=="GLY" ) {
   n[name]++
   group[name,n[name]]=id
   if(H[name]>0.0) etot+=H[name]*Tfactor
  }
  else {
   name="BB"
   n[name]++
   group[name,n[name]]=id
  }
 }
 END{
  ntot=0
  for(name in n) {
   ntot+=n[name]
   printf("group "name" id ")
   for(i=1;i<=n[name];i++) {
    printf(group[name,i]" ")
   }
   print ""
  }
  print "# ntot =",ntot
  pos=60.0
  sigma=8.0
  for(name in n) {
   if(name=="BB")     { cut=sigma ; eps=0.239            ; c3=1.5*eps*sigma**3 ; c9=0.5*eps*sigma**9 } # attraction with cutoff
   else if(H[name]>0) { cut=30.0  ; eps=+Tfactor*H[name] ; c3=1.5*eps*sigma**3 ; c9=0.5*eps*sigma**9 } # attraction
   else               { cut=30.0  ; eps=-Tfactor*H[name] ; c3=0.0              ; c9=1.0*eps*sigma**9 } # repulsion
   printf("fix wall"name" "name" wall/lj93new zhi "pos" "c3" "c9" "cut" units box\n")
  }
  for(name in n) printf("fix_modify wall"name" energy yes\n")
  print "# etot = ", etot, "kcal/mol"
  print "# etot/efolding = ", etot/(0.239*224)
 }
' $pdbfile # cg_2WGO_0.pdb
