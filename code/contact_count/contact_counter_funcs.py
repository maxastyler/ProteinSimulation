import numpy as np
from scipy.optimize import curve_fit 
import mdtraj as md
from itertools import combinations
import matplotlib.pyplot as plt
import matplotlib

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
