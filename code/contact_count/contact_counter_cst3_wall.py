#!/bin/env python3

#This script takes a load of .xtc trajectory files at different temperatures and works out the fraction of native contacts that is present, then graphs this.
import numpy as np
from scipy.optimize import curve_fit 
import mdtraj as md
from itertools import combinations
import matplotlib.pyplot as plt
import matplotlib
from contact_counter_funcs import *

chunked = 100
energy='0.244'
cg_path='../../wall_oil_energy/cst3-'+energy+'/'
sim_path=cg_path

def sigmoid(x, A, B, C, D):
    return A/(1+np.exp(-B*(x-C))) + D

native = md.load(cg_path+'cg_3gax_wall.pdb')
traj=md.load(sim_path+'wall{0}-124.xtc'.format(energy), top=cg_path+'cg_3gax_wall.pdb')
print(traj.top.select("resSeq 50 to 70"))
alpha_contacts=best_hummer_q(traj, native, traj.top.select("resSeq 10 to 26"))
beta_contacts=best_hummer_q(traj, native, traj.top.select("(resSeq 82 to 106) or (resSeq 31 to 62)"))
all_contacts=best_hummer_q(traj, native)
print(native[0])
if __name__=='__main__':

    plt.rc('text', usetex=True)
    matplotlib.rcParams.update({'font.size': 16})

    #Average the contacts, remove the first few values before the protein gets to equilibrium
    all_x, all_y = average_array(range(len(all_contacts)), all_contacts, chunked)
    alpha_x, alpha_y = average_array(range(len(alpha_contacts)), alpha_contacts, chunked)
    beta_x, beta_y = average_array(range(len(beta_contacts)), beta_contacts, chunked)
    plt.plot(all_x, all_y, label="All Contacts")
    plt.plot(alpha_x, alpha_y, label=r"$\alpha$ Contacts")
    plt.plot(beta_x, beta_y, label=r"$\beta$ Contacts")
    #plt.title(r'Fraction of Native Contacts Over Simulation Time for CST3 (E_{wall}=0.244k_BT)')
    plt.xlabel('Timestep')
    plt.ylabel('Native Contact Fraction')
    plt.legend()
    plt.show()
