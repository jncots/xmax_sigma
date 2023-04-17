# is it OK if I give you this file? 
# Anton can then bin the energies :slightly_smiling_face:
# the keys of the HDF5 file are PDG IDs (±12, ±13, ±14), 
# as well as 'num_primaries' which tells you the number 
# of simulated protons. for each PDG ID key there are 4 
# observation levels (1, 2, 3, 4), which correspond to 
# 143, 638, 1167, and 1195 g/cm2. then for each PDG key 
# and observation level there are two columns — energy and theta.
# so for example, if you load this file under the name 
# corsika_data, all muon energies at 143 g/cm2 can be 
# extracted as np.array(corsika_data['13']['1']['energy [GeV]']) . 
# and then num_primaries can be used to normalize it


from pathlib import Path
import h5py
import numpy as np
from particle import Particle


base_dir = Path("/hetghome/antonpr/xmax_sigma/data_comparison")
h5file = base_dir/"corsika_leptons_trial.h5"


all_particles_dict = {p.pdgid: p for p in Particle.findall()}

xdepth_list = np.array([143, 638, 1167, 1195], dtype=np.float32)
pdg_list = [-12, 12, -13, 13, -14, 14]
pdg_dict = { pdg: all_particles_dict[pdg] for pdg in pdg_list}

def corsika_hist(en_bins):
    energy_hist = dict()
    with h5py.File(h5file, "r") as corsika_data:
        num_primaries = corsika_data["num_primaries"]
        for pdg in pdg_dict:
            xd_dict = []
            for i, xdepth in enumerate(xdepth_list):
                mcdata = np.array(corsika_data[str(pdg)][str(i+1)]["energy [GeV]"])
                hist, bin_edges = np.histogram(mcdata, bins = en_bins)
                xd_dict.append((hist/num_primaries, bin_edges, xdepth))         
            energy_hist[pdg] = (xd_dict, f"${pdg_dict[pdg].latex_name}$")
    return energy_hist


def combined_data(energy_hist, pdgs, ixdepth):
    
    hist = None
    label = ""
    for pdg in pdgs:
        
        if hist is None:
            hist = energy_hist[pdg][0][ixdepth][0]
        else:
            hist = energy_hist[pdg][0][ixdepth][0] + hist    
        
        label += f"+{energy_hist[pdg][1]}"
    
    label_p = f"{label[1:]}"
    label_x = f"X={energy_hist[pdg][0][ixdepth][2]}"
    
    return hist, energy_hist[pdg][0][ixdepth][1], label_p, label_x


def corsika_en_theta_2dhist(en_bins, theta_bins):
    energy_hist = dict()
    with h5py.File(h5file, "r") as corsika_data:
        num_primaries = corsika_data["num_primaries"]
        for pdg in pdg_dict:
            xd_dict = []
            for i, xdepth in enumerate(xdepth_list):
                mc_energy = np.array(corsika_data[str(pdg)][str(i+1)]["energy [GeV]"])
                mc_theta = np.array(corsika_data[str(pdg)][str(i+1)]["theta [rad]"])
                
                # If Omega is used for binning
                # mc_omega = 2*np.pi*(1 - np.cos(mc_theta))
                # mc_theta = mc_omega
                
                hist, en_edges, theta_edges = np.histogram2d(mc_energy, mc_theta, 
                                                             bins = [en_bins, theta_bins])
                
                # If Omega is used for binning
                # omega_edges = np.arccos(1 - theta_edges/(2*np.pi))
                # theta_edges = omega_edges
                
                # if pdg in [-13, 13]:
                #     print(f"Hist = {hist}")
                #     print(f"Hist = {np.sqrt(hist)/hist}")
                #     print(f"Xdepth = {xdepth}, pdg = ${pdg_dict[pdg].latex_name}$")
                #     # print(f"Hist = {theta_edges}")
                
                xd_dict.append((hist/num_primaries, en_edges, theta_edges, xdepth, np.sqrt(hist)/num_primaries))         
            energy_hist[pdg] = (xd_dict, f"${pdg_dict[pdg].latex_name}$")
    return energy_hist


def combined_ang_data(energy_hist, pdgs):    
    dist_xdepth = []
    for ixdepth in range(4):
        
        dist_en = []
        for ind_energy in range(energy_hist[13][0][ixdepth][1].size - 1):
            
            hist = None
            hist_error = None
            label = ""
            for pdg in pdgs:    
                ang_dist = energy_hist[pdg][0][ixdepth][0][ind_energy]
                ang_dist_error = energy_hist[pdg][0][ixdepth][4][ind_energy]
                en_bins = energy_hist[pdg][0][ixdepth][1]
                cur_en_label = f"[{en_bins[ind_energy]}, {en_bins[ind_energy + 1]}]"
                cur_en = np.array([en_bins[ind_energy], en_bins[ind_energy + 1]])
                ang_bins = energy_hist[pdg][0][ixdepth][2]
                xdepth = energy_hist[pdg][0][ixdepth][3]
        
                if hist is None:
                    hist = ang_dist
                else:
                    hist = ang_dist + hist   

                if hist_error is None:
                    hist_error = ang_dist_error**2
                else:
                    hist_error = ang_dist_error**2 + hist_error   
                
                # print(f"pdg = {pdg}, sum = {np.sum(hist)}, size = {hist.size}")


                
                label += f"+{energy_hist[pdg][1]}"
    
                label_p = f"{label[1:]}"
                label_x = f"{xdepth}"

            hist_error = np.sqrt(hist_error)
            # print(f"hist_error = {hist_error/hist}")
            # print(f"ang_dist_error = {ang_dist_error/ang_dist}")

            dist_en.append((hist, ang_bins, cur_en_label, label_x, label_p, cur_en, hist_error))
        dist_xdepth.append(dist_en)    
    return dist_xdepth
  
    
if __name__ == "__main__":
    import nexusformat.nexus as nx
    # #Print a structure of db
    # f = nx.nxload(h5file)
    # print(f.tree)
    
        