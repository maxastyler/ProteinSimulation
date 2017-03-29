#plumed sum_hills --hills hills_124.dat --outfile fes_124_ --idw all --kt 1.031 --mintozero --stride 500
plumed sum_hills --histo colvar.dat --kt 1.031 --idw beta_beta,alpha_beta --sigma 1,1 --outhisto histo_beta_core.dat

