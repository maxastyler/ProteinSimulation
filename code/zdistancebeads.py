#!/usr/bin/env python3

#Script to sort out the bead ids that interact with the wall - usage:
# ./zdistancebeads.py wall_cyc_lmp_input_0.248.dat > plumed_temp_meta.dat
import sys
file_to_open=sys.argv[1]

res_id={}
res_energy={}

with open(file_to_open) as open_file:
    for line in open_file:
        sp_line = line.split()
        if sp_line[0]=="group" and sp_line[1]!="BB":
            res_id[sp_line[1]]=sp_line[3:]
        elif sp_line[0]=="fix" and sp_line[2]!="BB":
            res_energy[sp_line[2]]=sp_line[6]

interaction_ids=[]

for key in res_energy:
    if res_energy[key]!=0:
        interaction_ids += res_id[key]

interaction_nums=[int(i) for i in interaction_ids]

interaction_nums.sort()

my_string=""
for i in interaction_nums:
    my_string+=","+str(i)
print(my_string[1:])
