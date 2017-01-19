"""
Use this to visualize the output from
the unit test.
"""
import numpy as np
import matplotlib.pyplot as plt

#Read in MF data data
data_path = "./output/mass_function/final_mass_function/final_mass_function.txt"
data = np.genfromtxt(data_path)
M_bins = data[:,:2]
lM = np.log10(np.mean(M_bins,1))
N = data[:,2]
N_err = data[:,3]

def plot_MF(lM,N,N_err):
    fig,ax = plt.subplots(1,1)
    ax.errorbar(lM,N,yerr=N_err)
    ax.set_yscale('log')
    ax.set_xlabel(r"$\log_{10}M\ [{\rm M_\odot}/h]$",fontsize=24)
    ylims = ax.get_ylim()
    ax.set_ylabel(r"${\rmNumber\ of\ halos}$",fontsize=24)
    plt.subplots_adjust(bottom=0.15,left=0.15,hspace=0.001)
    plt.title(r"${\rm Halo\ mass\ function\ at\ }z=0.5$",fontsize=18)
    plt.gcf().savefig("plots/MF_example.png")
    plt.show()
    plt.close()

plot_MF(lM,N,N_err)

#Read in the HH_CF data
data_path = "./output/halohalo_correlation_function/final_hhcf/final_hhcf.txt"
data = np.genfromtxt(data_path)
R = data[:,0]
xi = data[:,1]
xi_err = data[:,2]
def plot_HHCF(R,xi,xi_err):
    fig,ax = plt.subplots(1,1)
    ax.errorbar(R,xi,yerr=xi_err)
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlabel(r"$R\ [{\rm Mpc}/h]$",fontsize=24)
    ylims = ax.get_ylim()
    ax.set_ylabel(r"$\xi_{\rm hh}$",fontsize=24)
    plt.subplots_adjust(bottom=0.15,left=0.15,hspace=0.001)
    plt.title(r"${\rm Halo-halo\ correlation\ function\ at\ }z=0.5$",fontsize=18)
    plt.gcf().savefig("plots/HHCF_example.png")
    plt.show()
    plt.close()
    return

plot_HHCF(R,xi,xi_err)

#Read in the HM_CF data
data_path = "./output/halomatter_correlation_function/full_hmcf/full_hmcf.txt"
data = np.genfromtxt(data_path)
R = data[:,0]
xi = data[:,3]
xi_err = data[:,4]
def plot_HMCF(R,xi,xi_err):
    fig,ax = plt.subplots(1,1)
    ax.errorbar(R,xi,yerr=xi_err)
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlabel(r"$R\ [{\rm Mpc}/h]$",fontsize=24)
    ylims = ax.get_ylim()
    ax.set_ylabel(r"$\xi_{\rm hm}$",fontsize=24)
    plt.subplots_adjust(bottom=0.15,left=0.15,hspace=0.001)
    plt.title(r"${\rm Halo-matter\ correlation\ function\ at\ }z=0.5$",fontsize=18)
    plt.gcf().savefig("plots/HMCF_example.png")
    plt.show()
    plt.close()
    return
plot_HMCF(R,xi,xi_err)
