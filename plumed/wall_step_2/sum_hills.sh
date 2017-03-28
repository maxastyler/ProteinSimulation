#plumed sum_hills --hills hills_124.dat --outfile fes_124_ --kt 1.031 --mintozero --stride 500
plumed sum_hills --histo colvar.dat --kt 1.031 --idw all,alpha_beta --sigma 4,4 --outhisto histo.dat
