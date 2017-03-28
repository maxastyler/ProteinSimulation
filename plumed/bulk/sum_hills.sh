#plumed sum_hills --hills hills_137.7.dat --outfile fes_137.7_ --kt 1.031 --stride 1000 --mintozero
plumed sum_hills --histo colvar.dat --kt 1.031 --idw alpha_beta,all --sigma 2,2 --outhisto histo.dat
