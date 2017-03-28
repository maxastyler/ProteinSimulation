#plumed sum_hills --hills hills_124.dat --outfile fes_124_ --kt 1.031 --stride 500 --mintozero
plumed sum_hills --histo colvar.dat --kt 1.031 --idw all,alpha_beta --sigma 2,2 --outhisto histo_all_core.dat
