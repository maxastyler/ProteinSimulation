#plumed sum_hills --hills hills_124.dat --outfile fes_124_ --kt 1.031 --mintozero --stride 500
plumed sum_hills --histo new_colvar.dat --kt 1.031 --idw alpha_beta,dcore.z --sigma 4,0.2 --outhisto histo.dat
