import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib
from matplotlib.ticker import LinearLocator, FormatStrFormatter
my_files = ["fes_124_{}.dat".format(i) for i in range(6)]

rc('text', usetex=True)

matplotlib.rcParams.update({'font.size': 16})

def extract_1d_file_data(f_path):
    x_data=[]
    y_data=[]
    with open(f_path) as my_file:
        for line in my_file:
            ls = line.split()
            if ls[0]!='#!':
                x_data.append(float(ls[0]))
                y_data.append(float(ls[1]))

    return x_data, y_data

def extract_2d_file_data(f_path):
    sections=[[]]
    with open(f_path) as my_file:
        for line in my_file:
            sp=line.split()
            if sp==[]: sections.append([])
            else:
                if sp[0]!='#!':
                    sections[-1].append((float(sp[0]), float(sp[1]), float(sp[2])))
    xs=[]
    ys=[]
    for line in sections[0]:
        xs.append(line[0])
    for sec in sections:
        ys.append(sec[0][1])
    xs = np.array(xs)
    ys = np.array(ys)
    xs, ys = np.meshgrid(xs, ys)
    zs=[]
    for sec in sections:
        zs.append([])
        for line in sec:
            zs[-1].append(line[2])
    zs = np.array(zs)
    return xs, ys, zs

def show_1d():
    xs=[]
    ys=[]

    for path in my_files:
        x,y=extract_1d_file_data(path)
        xs.append(x)
        ys.append(y)

        for i in range(len(xs)):
            plt.plot(xs[i], ys[i])
    plt.xlabel("distance from interface (nm)")
    plt.ylabel("energy (kT)")

    plt.show()

def show_contour():
    fig=plt.figure()
    xs, ys, zs = extract_2d_file_data('histo_core_beta.dat')
    xs/=133
    ys/=46
    newzs=np.copy(zs)
    newzs[np.isnan(newzs)]=0
    newzs[newzs==np.inf]=0
    newzs[newzs==-np.inf]=0
    zs[np.isnan(zs)]=np.max(newzs)
    zs[zs==np.inf]=np.max(newzs)
    levels=np.arange(np.min(zs), 0, 0.2)
    csf=plt.contourf(xs, ys, zs, levels=levels)
    cs=plt.contour(xs, ys, zs, levels=levels, colors=('k'))
    cbar=fig.colorbar(csf)
    cbar.set_label(r'Free Energy ($k_BT$)')
    plt.xlabel(r"$\beta$-sheet native contact fraction")
    plt.ylabel(r"$\alpha\leftrightarrow\beta$ native contact fraction")
    plt.show()

def show_2d():
    fig=plt.figure()
    ax = fig.gca(projection='3d')
    xs, ys, zs = extract_2d_file_data('fes_124_111.dat')
    surf = ax.plot_surface(xs, ys, zs, cmap=cm.coolwarm, linewidth=0, antialiased=False)

    plt.show()

show_contour()
