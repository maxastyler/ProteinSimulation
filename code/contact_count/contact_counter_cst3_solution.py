#!/bin/env python3

#This script takes a load of .xtc trajectory files at different temperatures and works out the fraction of native contacts that is present, then graphs this.
import numpy as np
from scipy.optimize import curve_fit 
import mdtraj as md
from itertools import combinations
import matplotlib.pyplot as plt

#List of temperatures that have been done by the script
#temps=[2, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 68, 72, 76, 80, 84, 88, 92, 96, 100]
#temps=[100, 105, 110, 115, 120, 125, 130, 135, 140, 145, 150, 155, 160, 165, 170, 175, 180, 185]
temp=124

#sigmoid function of the form f(x)=A/(1+B*exp(-C*x)+D
cg_path='../../protein/cystatin/3GAX/7-coarsegrained/'
sim_path='../../lammps_scripts/cst3_solution/'

def sigmoid(x, A, B, C, D):
    return A/(1+np.exp(-B*(x-C))) + D

def best_hummer_q(traj, native):
    """Compute the fraction of native contacts according the definition from
    Best, Hummer and Eaton [1]
    ----------
    Parameters
    ----------
    traj : md.Trajectory
        The trajectory to do the computation for
    native : md.Trajectory
        The 'native state'. This can be an entire trajecory, or just a single frame.
        Only the first conformation is used
    -------
    Returns
    -------
    q : np.array, shape=(len(traj),)
        The fraction of native contacts in each frame of `traj`
    ----------
    References
    ----------
    ..[1] Best, Hummer, and Eaton, "Native contacts determine protein folding
          mechanisms in atomistic simulations" PNAS (2013)
    """
    BETA_CONST = 50  # 1/nm
    LAMBDA_CONST = 1.8
    NATIVE_CUTOFF = 0.5  # nanometers
    # get the indices of all of the heavy atoms
    heavy = native.topology.select_atom_indices('heavy')
    # get the pairs of heavy atoms which are farther than 3
    # residues apart
    heavy_pairs = np.array(
        [(i,j) for (i,j) in combinations(heavy, 2)
            if abs(native.topology.atom(i).residue.index - \
                   native.topology.atom(j).residue.index) > 3])
    # compute the distances between these pairs in the native state
    heavy_pairs_distances = md.compute_distances(native[0], heavy_pairs)[0]
    # and get the pairs s.t. the distance is less than NATIVE_CUTOFF
    native_contacts = heavy_pairs[heavy_pairs_distances < NATIVE_CUTOFF]
    print("Number of native contacts", len(native_contacts))
    # now compute these distances for the whole trajectory
    r = md.compute_distances(traj, native_contacts)
    # and recompute them for just the native state
    r0 = md.compute_distances(native[0], native_contacts)
    q = np.mean(1.0 / (1 + np.exp(BETA_CONST * (r - LAMBDA_CONST * r0))), axis=1)
    return q

native = md.load(cg_path+'cg_3gax_filtered.pdb')
traj = md.load(sim_path+'run.{0}.xtc'.format(temp), top=cg_path+'cg_3gax_filtered.pdb')
contacts = best_hummer_q(traj, native)
#traj=md.load(sim_path+'run.60.xtc', top=cg_path+'cg_3gax_filtered.pdb')
if __name__=='__main__':

    #Average the contacts, remove the first few values before the protein gets to equilibrium
    ys = np.average(contacts)
    
    plt.plot(contacts)
    plt.title('Fraction of Native Contacts for Coarsegrained CST3 vs Temperature')
    plt.xlabel('Temperature (LAMMPS units)')
    plt.ylabel('Native Contact Fraction')
    #fitted_curve=curve_fit(sigmoid, temps, ys, p0=[-0.8, 1, 150, 1])
    #print(fitted_curve[0])
    #gen_ys=[sigmoid(i, fitted_curve[0][0], fitted_curve[0][1], fitted_curve[0][2], fitted_curve[0][3]) for i in temps]
    #residuals=[gen_ys[i]-ys[i] for i in range(len(ys))]
    #plt.plot(temps, gen_ys)
    #plt.plot(temps, residuals)
    plt.show()
