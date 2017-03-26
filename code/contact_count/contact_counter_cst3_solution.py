#!/bin/env python3

#This script takes a load of .xtc trajectory files at different temperatures and works out the fraction of native contacts that is present, then graphs this.
import numpy as np
from scipy.optimize import curve_fit 
import mdtraj as md
from itertools import combinations
import matplotlib.pyplot as plt
from contact_counter_funcs import *

#List of temperatures that have been done by the script
#temps=[2, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 68, 72, 76, 80, 84, 88, 92, 96, 100]
#temps=[100, 105, 110, 115, 120, 125, 130, 135, 140, 145, 150, 155, 160, 165, 170, 175, 180, 185]
temp=124
chunked=700

#sigmoid function of the form f(x)=A/(1+B*exp(-C*x)+D
cg_path='../../protein/cystatin/3GAX/7-coarsegrained/'
sim_path='../../lammps_scripts/cst3_solution/'

native = md.load(cg_path+'cg_3gax_filtered.pdb')
traj = md.load(sim_path+'run.{0}.xtc'.format(temp), top=cg_path+'cg_3gax_filtered.pdb')
alpha_contacts=best_hummer_q(traj, native, traj.top.select("resSeq 10 to 26"))
beta_contacts=best_hummer_q(traj, native, traj.top.select("(resSeq 82 to 106) or (resSeq 31 to 62)"))
all_contacts=best_hummer_q(traj, native)
#traj=md.load(sim_path+'run.60.xtc', top=cg_path+'cg_3gax_filtered.pdb')
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
