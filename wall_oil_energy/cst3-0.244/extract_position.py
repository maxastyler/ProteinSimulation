#!/usr/bin/env python3
import sys
import mdtraj as md
top=sys.argv[1]
traj=sys.argv[2]
frame_no=int(sys.argv[3])

prot=md.load(traj, top=top)

#Print out each particle
for i in range(len(prot[frame_no]._xyz[0][:-1])):
    part = 10*prot[frame_no]._xyz[0][i] #10* to convert from gromacs to lammps units
    print("{0} 1 {0} {1:.3f} {2:.3f} {3:.3f} 0 0 0".format(i+1, part[0], part[1], part[2]))

#Print out the wall particle
part=10*prot[frame_no]._xyz[0][-1]
print("{0} 0 {0} {1:.3f} {2:.3f} {3:.3f} 0 0 0".format(len(prot[frame_no]._xyz[0]), part[0], part[1], part[2]))
