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

# 1       2.00000000E+06       5.69000815E+01
#           2       1.50000000E+06       1.24139149E+02
#           3       5.00000000E+05       5.52958799E+02
#           4       0.00000000E+00       1.03610000E+03

# h=0, X=1036.0992336839001
# h=5, X=552.9465437966635
# h=15, X=124.14014337836426
# h=20, X=56.90063537140265

# h=0, X=1195.9290875457918
# h=5, X=638.1217239628874
# h=15, X=143.18394585718303
# h=20, X=65.61184736576809


all_particles_dict = {p.pdgid: p for p in Particle.findall()}

# xdepth_list = np.array([124, 552, 1036], dtype=np.float32)
# xdepth_list = np.array([1, 10, 100, 300, 700, 
#                                1000, 3000, 5000, 7000, 10000], dtype=np.float32)
# xdepth_list = np.array([910, 950, 1000, 1033], dtype=np.float32)

# xdepth_list = np.array([10, 50, 100, 133], dtype=np.float32)
xdepth_list = np.array([1, 10, 50, 100, 300, 700, 1000, 1065], dtype=np.float32)
pdg_list = [-12, 12, -13, 13, -14, 14]
pdg_dict = { pdg: all_particles_dict[pdg] for pdg in pdg_list}

def corsika_hist_en(en_bins, h5file = "corsika_leptons_trial.h5"):
    # base_dir = Path("/hetghome/antonpr/xmax_sigma/data_comparison")
    # h5file = base_dir / h5file
    # h5file = 
    energy_hist = dict()
    print(h5file)
    with h5py.File(h5file, "r") as corsika_data:
        
        num_primaries = corsika_data["num_primaries"]
        print(f"{str(h5file.name)}: Number of primaries = {np.array(num_primaries)*1.0:e}")
        
        for pdg in pdg_dict:
            xd_dict = []
            for i, xdepth in enumerate(xdepth_list):
                # rrr = np.array(corsika_data[str(pdg)][str(i)]["energy [GeV]"])
                # print(f"res = {rrr}")
                mcdata = np.array(corsika_data[str(pdg)][str(i)]["energy [GeV]"])
                hist, bin_edges = np.histogram(mcdata, bins = en_bins)
                num_tot = np.sum(hist)
                xd_dict.append((hist, bin_edges, xdepth, np.array(num_primaries)*1.0, num_tot))        
            energy_hist[pdg] = (xd_dict, f"${pdg_dict[pdg].latex_name}$")
    return energy_hist


def combined_data_en(energy_hist, pdgs, ixdepth):
    
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
    
    return (hist, energy_hist[pdg][0][ixdepth][1], 
            label_p, label_x, 
            energy_hist[pdg][0][ixdepth][3], # num_primaries
            )


def corsika_en_theta_2dhist(en_bins, theta_bins, h5file):
    
    base_dir = Path("/hetghome/antonpr/xmax_sigma/data_comparison")
    h5file = base_dir / h5file
    
    energy_hist = dict()
    with h5py.File(h5file, "r") as corsika_data:
        num_primaries = corsika_data["num_primaries"]
        for pdg in pdg_dict:
            print(f"pdg = {pdg}")
            ix = [int(key) for key in corsika_data[str(pdg)].keys()]
            xd_dict = []
            for i, xdepth in enumerate(xdepth_list):
                mc_energy = np.array(corsika_data[str(pdg)][str(ix[i])]["energy [GeV]"])
                mc_theta = np.array(corsika_data[str(pdg)][str(ix[i])]["theta [rad]"])
                print(f"xdepth={xdepth}, number={len(mc_energy)*1e0:.4e}")
                
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
    for ixdepth in range(3):
        
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

def print_h5file_structure(h5file):
    """Print a structure of db
    Args:
        h5file (pathlib.Path): path to *.h5 file
    """
    import nexusformat.nexus as nx
    f = nx.nxload(h5file)
    print(f.tree)
    
    # with h5py.File(h5file, "r") as corsika_data:
    #     num_primaries = corsika_data["num_primaries"]
    #     print(np.array(num_primaries))
      
    
if __name__ == "__main__":
    import nexusformat.nexus as nx
    base_dir = Path("/hetghome/antonpr/projects/mceq_vs/data/corsika")
    h5file = base_dir/"03_single_muon_shower/01_single_muon_run/01_single_muon.h5"
    print_h5file_structure(h5file)
    
        