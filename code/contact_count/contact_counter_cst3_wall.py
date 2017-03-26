#!/bin/env python3

#This script takes a load of .xtc trajectory files at different temperatures and works out the fraction of native contacts that is present, then graphs this.
import numpy as np
from scipy.optimize import curve_fit 
import mdtraj as md
from itertools import combinations
import matplotlib.pyplot as plt

chunked = 100
energy='0.244'
cg_path='../../wall_oil_energy/cst3-'+energy+'/'
sim_path=cg_path

def sigmoid(x, A, B, C, D):
    return A/(1+np.exp(-B*(x-C))) + D

def best_hummer_q(traj, native, atom_indices=None):
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
    #NATIVE_CUTOFF = 0.6
    # get the indices of all of the heavy atoms
    if atom_indices is None: atom_indices = native.topology.select_atom_indices('heavy')
    # get the pairs of heavy atoms which are farther than 3
    # residues apart
    index_pairs = np.array(
        [(i,j) for (i,j) in combinations(atom_indices, 2)
            if abs(native.topology.atom(i).residue.index - \
                   native.topology.atom(j).residue.index) > 3])
    # compute the distances between these pairs in the native state
    pairs_distances = md.compute_distances(native[0], index_pairs)[0]
    # and get the pairs s.t. the distance is less than NATIVE_CUTOFF
    native_contacts = index_pairs[pairs_distances < NATIVE_CUTOFF]
    print("Number of native contacts", len(native_contacts))
    # now compute these distances for the whole trajectory
    r = md.compute_distances(traj, native_contacts)
    # and recompute them for just the native state
    r0 = md.compute_distances(native[0], native_contacts)
    q = np.mean(1.0 / (1 + np.exp(BETA_CONST * (r - LAMBDA_CONST * r0))), axis=1)
    return q

def average_array(xs, ys, n):
    remainder = len(xs)%n
    added = int(remainder!=0)
    new_len=added+int(len(xs)/n)
    x_chunked=[xs[i:i+n] for i in range(0, len(xs), n)]
    y_chunked=[ys[i:i+n] for i in range(0, len(ys), n)]
    x_avg = [sum(i)/len(i) for i in x_chunked]
    y_avg = [sum(i)/len(i) for i in y_chunked]
    return x_avg, y_avg


native = md.load(cg_path+'cg_3gax_wall.pdb')
traj=md.load(sim_path+'wall{0}-124.xtc'.format(energy), top=cg_path+'cg_3gax_wall.pdb')
print(traj.top.select("resSeq 50 to 70"))
alpha_contacts=best_hummer_q(traj, native, traj.top.select("resSeq 10 to 26"))
beta_contacts=best_hummer_q(traj, native, traj.top.select("(resSeq 82 to 106) or (resSeq 31 to 62)"))
all_contacts=best_hummer_q(traj, native)
print(native[0])
if __name__=='__main__':

    plt.rc('text', usetex=True)
    #Average the contacts, remove the first few values before the protein gets to equilibrium
    all_x, all_y = average_array(range(len(all_contacts)), all_contacts, chunked)
    alpha_x, alpha_y = average_array(range(len(alpha_contacts)), alpha_contacts, chunked)
    beta_x, beta_y = average_array(range(len(beta_contacts)), beta_contacts, chunked)
    plt.plot(all_x, all_y)
    plt.plot(alpha_x, alpha_y)
    plt.plot(beta_x, beta_y)
    plt.title(r'Fraction of Native Contacts Over Simulation Time for CST3 (E_{wall}=0.248)')
    plt.xlabel('Timestep')
    plt.ylabel('Native Contact Fraction')
    plt.show()
