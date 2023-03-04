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
                
        self.cascade_driver = cascade_driver        
        self.final_particles = cascade_driver.get_final_particles().valid()
        self.pack_data()

    def print_stats(self):
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
        print(f"  Size of cascade_driver object = {asizeof.asizeof(self.cascade_driver)/(1024**2):.2f} Mb")
        self.energy_conservation()
    
    
    def energy_conservation(self):     
        etot = np.sum(self.energy_data)
        init_energy = self.cascade_driver.initial_energy
        conservation = (init_energy - etot)/init_energy
        
        print("\nEnergy conservation in cascade:")
        print(f"  Initial energy = {init_energy:.5e} GeV")
        print(f"  Energy in final particles = {etot:.5e} GeV")
        print(f"  Relative loss(+)/gain(-) {conservation:.3e}")
        
    
    def pack_data(self):
        self.pid_data = self.final_particles.pid
        self.energy_data = self.final_particles.energy
        self.xdepth_data = self.final_particles.xdepth

        # for particle in self.particles:
        #     self.pid_data.append(particle.pid)
        #     self.energy_data.append(particle.energy)
        #     self.xdepth_data.append(particle.xdepth)
        
    def calc_height(self):
        xconv = XdepthConversion()
        xconv.set_length_unit("km")
        
        self.height_data = []
        for i, ppid in enumerate(self.pid_data):
            self.height_data.append(xconv.convert_x2h(self.xdepth_data[i]))
            
        self._is_height = True

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
        
    def plot_height_list(self, pids=None, 
                         all_pids = None, 
                         nbins = 100, 
                         xrange = None,
                         energy_range = None):
        
        if not self._is_height:
            self.calc_height()
        
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