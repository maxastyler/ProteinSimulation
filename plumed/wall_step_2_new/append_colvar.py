o_lines=[]
i_lines=[]
with open("colvar.dat") as o_file:
    with open("colvar_124-meta.dat") as i_file:
        for line in o_file:
            if line[:2]!='#!':
                o_lines.append(line[:-1])
        for line in i_file:
            if line[:2]!='#!':
                i_lines.append(line.split()[1])
for i in range(len(o_lines)):
    o_lines[i]+=" " + i_lines[i] + "\n"
print(o_lines[0:10])
with open("new_colvar.dat", "w") as outfile:
    for line in o_lines:
        outfile.write(line)

