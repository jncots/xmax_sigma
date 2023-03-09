import matplotlib.pyplot as plt
import numpy as np
from xdepth_on_table import XdepthOnTable
from xdepth_conversion import XdepthConversion
from MCEq.geometry.density_profiles import CorsikaAtmosphere
from cascade_driver import CascadeDriver
from pympler import asizeof
import particle
from pdg_pid_map import PdgLists


from copy import copy


class CascadeAnalysis:
    def __init__(self, cascade_driver):
        
        self.all_pdgs = {
            int(p.pdgid): f"${p.latex_name}$" for p in particle.Particle.findall()
        }
        
        self.all_pdgs_mass = {
            int(p.pdgid): p.mass*1e-3 if p.mass is not None else 0e0 for p in particle.Particle.findall()
        }
                
        self.cascade_driver = cascade_driver        
        self.final_particles = cascade_driver.get_final_particles().valid()
        self.pack_data()

    def print_stats(self):
        if self.cascade_driver.runs_number > 1:
            print(f"Number of runs = {self.cascade_driver.runs_number}")
            
        print(f"Initial state:")
        print(f"  {self.all_pdgs[self.cascade_driver.initial_pdg]}({self.cascade_driver.initial_pdg})"
              f" with energy = {self.cascade_driver.initial_energy:.3e}")                
        
        print(f"\nFinal state:")
        print(f"  Number of final particles = {len(self.cascade_driver.final_stack)}")
        print(f"  Number of interactions = {self.cascade_driver.number_of_interactions}")
        print(f"  Number of decays = {self.cascade_driver.number_of_decays}")
        print(f"  Max number of generations = {np.max(self.cascade_driver.final_stack.valid().generation_num[:])}")
        print(f"\n  Max xdepth = {self.cascade_driver.xdepth_getter.max_xdepth}")
        print(f"  Exectution time = {self.cascade_driver.loop_execution_time:.2f} s")
        if self.cascade_driver.runs_number > 1:
            exec_time = self.cascade_driver.loop_execution_time/self.cascade_driver.runs_number
            print(f"  Exectution time per run = {exec_time:.2f} s")
        print(f"  Size of cascade_driver object = {asizeof.asizeof(self.cascade_driver)/(1024**2):.2f} Mb")
        self.energy_conservation()
    
    
    def energy_conservation(self):     
        etot = np.sum(self.energy_data)/self.cascade_driver.runs_number
        init_energy = self.cascade_driver.initial_energy
        conservation = (init_energy - etot)/init_energy
        
        print("\nEnergy conservation in cascade:")
        print(f"  Initial energy = {init_energy:.5e} GeV")
        print(f"  Energy in final particles = {etot:.5e} GeV")
        print(f"  Relative loss(+)/gain(-) {conservation:.3e}")
    
    
    # def particle_flight(self):
        
    #     print(f"xdepth_decay"
    #         np.where(self.final_particles.xdepth_decay > 
    #                  self.cascade_driver.xdepth_getter.max_xdepth)[0])
        
        
            
    
    def pack_data(self):
        self.pid_data = self.final_particles.pid
        self.energy_data = self.final_particles.energy
        self.xdepth_data = self.final_particles.xdepth
        self.xdepth_stop = self.final_particles.xdepth_stop
        self.height_data = self.cascade_driver.xdepth_getter.xdepth_on_table.convert_x2h(self.final_particles.xdepth)
        self.height_data = self.height_data * 1e-5
        
        self.all_pids = list(set(self.pid_data))
        self.all_pids.sort(key=lambda x: abs(x))
        
        
    def digitize(self):
        
        xdepth_grid_bins = np.linspace(0.1, self.cascade_driver.xdepth_getter.max_xdepth + 1, 101)
        xdepth_grid_centers = (xdepth_grid_bins[0:-1] + xdepth_grid_bins[1:])/2
        xdepth_grid_steps = xdepth_grid_bins[1:] - xdepth_grid_bins[0:-1]
                
        energy_grid_bins = np.geomspace(1e-1, 1e11, 121)
        energy_grid_centers = (energy_grid_bins[0:-1] + energy_grid_bins[1:])/2
        energy_grid_steps = energy_grid_bins[1:] - energy_grid_bins[0:-1]
        
        
        hist_dict = {}
        for pid in self.all_pids:            
            en_vec = self.energy_data[np.where(self.pid_data == pid)[0]]
            xd_vec = self.xdepth_data[np.where(self.pid_data == pid)[0]]
            
            
            npart = np.histogram2d(xd_vec, en_vec, 
                           bins=(xdepth_grid_bins, energy_grid_bins))
            
            hist_dict[pid] = npart
        
        self.egrid = energy_grid_centers
        self.xgrid = xdepth_grid_centers
        self.hist_dict = hist_dict   
            
        
    def plot_pid(self, from_=None, to_=None):

        pid_dist = dict()

        for pid in self.final_particles.pid:
            pid_dist[pid] = pid_dist.get(pid, 0) + 1

        pid_dist = dict(
            sorted(pid_dist.items(), key=lambda item: item[1], reverse=True)
        )
        print(pid_dist)
        
        mceq_particles = PdgLists().mceq_particles
        for pdg in pid_dist.keys():
            if pdg not in mceq_particles:
                print(f"pdg = {pdg} is not among mceq_particles")
                

        ptypes = [self.all_pdgs[i] for i in pid_dist.keys()]
        pnum = list(pid_dist.values())
        plt.bar(ptypes[from_:to_], pnum[from_:to_])
        plt.title("Particle type distribution")

    def plot_energy(self, pid=None):
        
        energy_data = []

        if pid:
            for i, ppid in enumerate(self.pid_data):
                if ppid == pid:
                    energy_data.append(self.energy_data[i])
        else:
            energy_data = self.final_particles.energy

        gr, cnt = np.histogram(np.log10(energy_data), 100)
        plt.semilogx()
        plt.step(10 ** cnt[:-1] * 1e9, gr)
        plt.title("Energy distribution")
        plt.xlabel("Energy, eV")
        print(
            f"Min = {np.min(energy_data)*1e9:0.2e} eV, Max = {np.max(energy_data)*1e9:0.2e} eV"
        )
        
    def plot_height(self, pid=None):
        
        if not self._is_height:
            self.calc_height()
            
        height_data = []    

        if pid:
            for i, ppid in enumerate(self.pid_data):
                if ppid == pid:
                    height_data.append(self.height_data[i])
        else:
            height_data = self.height_data

        print(f"Min = {np.min(height_data):.2f} km, Max = {np.max(height_data):.2f} km")

        plt.stairs(*np.histogram(height_data, 100))
        plt.title("Height distribution")
        plt.xlabel("Height, km") 
        
    def plot_xdepth(self, pid=None):
        xdepth_data = []    

        if pid:
            for i, ppid in enumerate(self.pid_data):
                if ppid == pid:
                    xdepth_data.append(self.xdepth_data[i])
        else:
            xdepth_data = self.xdepth_data

        print(f"Min = {np.min(xdepth_data):.2f}, Max = {np.max(xdepth_data):.2f}")

        plt.stairs(*np.histogram(xdepth_data, 100))
        plt.title("Height distribution")
        plt.xlabel("Xdepth, gm/cm2")         

    def plot_energy_list(self, pids=None, all_pids = None, nbins = 100, xrange = None):
        
        energy_data = []

        if pids:
            for pid in pids:
                energy_data_pid = []
                for i, ppid in enumerate(self.pid_data):
                    if ppid == pid:
                        energy_data_pid.append(self.energy_data[i])
                energy_data.append(energy_data_pid)        
            
        if xrange:
            xrange = (np.log10(xrange[0]*1e-9), np.log10(xrange[1]*1e-9))
                    

        plt.title("Energy distribution")
        plt.xlabel("Energy, eV") 
        plt.semilogx()
            
        for i, pid in enumerate(pids):
            gr, cnt = np.histogram(np.log10(energy_data[i]), bins = nbins, range = xrange)
            plt.semilogx()
            plt.step(10 ** cnt[:-1] * 1e9, gr, label = f"{self.all_pdgs[pid]}")
            
        if all_pids:
            gr, cnt = np.histogram(np.log10(self.energy_data), bins = nbins, range = xrange)
            plt.semilogx()
            plt.step(10 ** cnt[:-1] * 1e9, gr, label = f"all")
                
            
        plt.legend()
        
    def plot_energy_kin_dist(self, pids=None, all_pids = None, nbins = 100, xrange = None, 
                             per_run = False):
        
        print(f"Muon mass = {self.all_pdgs_mass[-13]}, {self.all_pdgs_mass[13]}")
        energy_data = []

        if pids:
            for pid in pids:
                energy_data_pid = []
                for i, ppid in enumerate(self.pid_data):
                    if ppid in pid:
                        energy_data_pid.append(self.energy_data[i] - self.all_pdgs_mass[ppid])
                energy_data.append(energy_data_pid)        
            
        if xrange:
            xrange = (np.log10(xrange[0]), np.log10(xrange[1]))
                    

        plt.title(r"Kinetic Energy distribution, $\frac{dN}{dE_{kin}}E_{kin}$")
        plt.xlabel("Energy, GeV") 
        plt.semilogx()
        
        if per_run:
            runs_number = self.cascade_driver.runs_number
        else:
            runs_number = 1 
            
        for i, pid in enumerate(pids):
            gr, cnt = np.histogram(np.log10(energy_data[i]), bins = nbins, range = xrange)
            dbin = 10 **cnt[1:] - 10 **cnt[0:-1]
            cbin = (cnt[1:] + cnt[0:-1])/2
            gr = ((gr/dbin)/runs_number)*(10**cbin)
            
            ss = ""
            for pp in pid:
                ss += f"+{self.all_pdgs[pp]}"
            
            plt.semilogx()
            plt.plot(10 ** cbin, gr, label = ss[1:])
            # plt.step(10 ** cbin, gr, where='mid', label = ss[1:])
            
        if all_pids:
            gr, cnt = np.histogram(np.log10(self.energy_data), bins = nbins, range = xrange)
            dbin = 10 **cnt[1:] - 10 **cnt[0:-1]
            cbin = (cnt[1:] + cnt[0:-1])/2
            gr = ((gr/dbin)/runs_number)*(10**cbin)
            plt.semilogx()
            plt.step(10 ** cbin, gr,  where='mid', label = f"all")
                
            
        plt.legend()
        plt.grid()
        plt.grid(b=True, which='minor', linestyle='--')    
        
    def plot_height_list(self, pids=None, 
                         all_pids = None, 
                         nbins = 100, 
                         xrange = None,
                         energy_range = None):
        
        
        
        height_data = []
        
        if energy_range:
            energy_range = (energy_range[0]*1e-9, energy_range[1]*1e-9)
            
        
        height_data_e = []
        pid_data_e = []
             
        if energy_range:
            for i, hh in enumerate(self.height_data):
                if energy_range[0] <= self.energy_data[i] <= energy_range[1]:
                    height_data_e.append(hh)
                    pid_data_e.append(self.pid_data[i])
        else:
            height_data_e = self.height_data
            pid_data_e = self.pid_data
            
        if pids:
            for pid in pids:
                height_data_pid = []
                for i, ppid in enumerate(pid_data_e):
                    if ppid == pid:
                        height_data_pid.append(height_data_e[i])
                height_data.append(height_data_pid)
        
        
        height_data_tot = height_data_e                     

        plt.title("Height distribution")
        plt.xlabel("Height, km") 
            
        for i, pid in enumerate(pids):
            gr, cnt = np.histogram(height_data[i], bins = nbins, range = xrange)
            plt.step(cnt[:-1], gr, label = f"{self.all_pdgs[pid]}")
            
        if all_pids:
            gr, cnt = np.histogram(height_data_tot, bins = nbins, range = xrange)
            plt.step(cnt[:-1], gr, label = f"all")
                
            
        plt.legend()
        
    def plot_xdepth_list(self, pids=None, 
                         all_pids = None, 
                         nbins = 100, 
                         xrange = None,
                         energy_range = None,
                         per_run = None):
        
        
        
        if pids is None:
            pids = self.all_pids
        
        print(pids)
        
        xdepth_data = []
        
        if energy_range:
            energy_range = (energy_range[0]*1e-9, energy_range[1]*1e-9)
            
        
        xdepth_data_e = []
        pid_data_e = []
             
        if energy_range:
            for i, hh in enumerate(self.xdepth_data):
                if energy_range[0] <= self.energy_data[i] <= energy_range[1]:
                    xdepth_data_e.append(hh)
                    pid_data_e.append(self.pid_data[i])
        else:
            xdepth_data_e = self.xdepth_data
            pid_data_e = self.pid_data
            
        if pids:
            for pid in pids:
                xdepth_data_pid = []
                for i, ppid in enumerate(pid_data_e):
                    if ppid == pid:
                        xdepth_data_pid.append(xdepth_data_e[i])
                xdepth_data.append(xdepth_data_pid)
        
        
        xdepth_data_tot = xdepth_data_e                    

        
        print(f"Min = {np.min(self.xdepth_data):.2f} g/cm2, Max = {np.max(self.xdepth_data):.2f} g/cm2")
        plt.title("Xdepth distribution")
        plt.xlabel("Xdepth, g/cm2") 
        
        if per_run:
            runs_number = self.cascade_driver.runs_number
        else:
            runs_number = 1    
        
        if pids:    
            for i, pid in enumerate(pids):
                gr, cnt = np.histogram(xdepth_data[i], bins = nbins, range = xrange)
                grsum = np.cumsum(gr)
                plt.step(cnt[:-1], gr/runs_number, label = f"{self.all_pdgs[pid]}")
                # plt.step(cnt[:-1], grsum/runs_number, label = f"cumul_{self.all_pdgs[pid]}")
            
        if all_pids:
            gr, cnt = np.histogram(xdepth_data_tot, bins = nbins, range = xrange)
            grsum = np.cumsum(gr)
            plt.step(cnt[:-1], gr/runs_number, label = f"all")
            # plt.step(cnt[:-1], grsum/runs_number, label = f"cumul_all")
                
            
        plt.legend()          
        
        
    def plot_xdepth_stop(self, pids=None, 
                         all_pids = None, 
                         nbins = 100, 
                         xrange = None,
                         energy_range = None,
                         per_run = None):
        
        
        
        if pids is None:
            pids = self.all_pids
        
        print(pids)
        
        xdepth_data = []
        
        if energy_range:
            energy_range = (energy_range[0]*1e-9, energy_range[1]*1e-9)
            
        
        xdepth_data_e = []
        pid_data_e = []
             
        if energy_range:
            for i, hh in enumerate(self.xdepth_stop):
                if energy_range[0] <= self.energy_data[i] <= energy_range[1]:
                    xdepth_data_e.append(hh)
                    pid_data_e.append(self.pid_data[i])
        else:
            xdepth_data_e = self.xdepth_stop
            pid_data_e = self.pid_data
            
        if pids:
            for pid in pids:
                xdepth_data_pid = []
                for i, ppid in enumerate(pid_data_e):
                    if ppid == pid:
                        xdepth_data_pid.append(xdepth_data_e[i])
                xdepth_data.append(xdepth_data_pid)
        
        
        xdepth_data_tot = xdepth_data_e                    

        
        print(f"Min = {np.min(self.xdepth_data):.2f} g/cm2, Max = {np.max(self.xdepth_data):.2f} g/cm2")
        plt.title("Xdepth stop distribution")
        plt.xlabel("Xdepth, g/cm2") 
        
        if per_run:
            runs_number = self.cascade_driver.runs_number
        else:
            runs_number = 1    
        
        if pids:    
            for i, pid in enumerate(pids):
                gr, cnt = np.histogram(xdepth_data[i], bins = nbins, range = xrange)
                grsum = np.cumsum(gr)
                plt.step(cnt[:-1], gr/runs_number, label = f"{self.all_pdgs[pid]}")
                # plt.step(cnt[:-1], grsum/runs_number, label = f"cumul_{self.all_pdgs[pid]}")
            
        if all_pids:
            gr, cnt = np.histogram(xdepth_data_tot, bins = nbins, range = xrange)
            grsum = np.cumsum(gr)
            plt.step(cnt[:-1], gr/runs_number, label = f"all")
            # plt.step(cnt[:-1], grsum/runs_number, label = f"cumul_all")
                
            
        plt.legend()                       