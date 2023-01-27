import matplotlib.pyplot as plt
import numpy as np
from cascade.xdepth_conversion import XdepthConversion
from copy import copy


class FilterParticleData:
    property_list = ["pid", "energy", "xdepth", "production_mode", "generation_number"]

    def __init__(self, property_):
        if property_ not in self.property_list:
            raise ValueError(f'Property "{property_}" is unknown')
        self.property = property_
        self.conversion = None
        self.condition = None
        self.filtered_data = None
        self._set_filter()

    def set_property(self, property_):
        if property_ not in self.property_list:
            raise ValueError(f'Property "{property_}" is unknown')
        self.property = property_

    def set_condition(self, condition_func):
        self.condition = condition_func
        self._set_filter()
        return self

    def set_conversion(self, conversion_func):
        self.conversion = conversion_func
        self._set_filter()
        return self

    def _set_filter(self):
        if not self.condition and not self.conversion:
            self.ffilter = self._filter0
            return

        if self.condition and not self.conversion:
            self.ffilter = self._filter1
            return

        if not self.condition and self.conversion:
            self.ffilter = self._filter2
            return

        if self.condition and self.conversion:
            self.ffilter = self._filter3
            return

    def _filter0(self, data):
        return [getattr(particle, self.property) for particle in data]

    def _filter1(self, data):
        return [
            getattr(particle, self.property)
            for particle in data
            if self.condition(particle)
        ]

    def _filter2(self, data):
        return [self.conversion(getattr(particle, self.property)) for particle in data]

    def _filter3(self, data):
        return [
            self.conversion(getattr(particle, self.property))
            for particle in data
            if self.condition(particle)
        ]

    def filter(self, data):
        self.filtered_data = self.ffilter(data)
        return self

    def result_data(self):
        return self.filtered_data

    def plot_data(self, bins_number=100):
        plt.stairs(*np.histogram(self.filtered_data, bins_number))
        plt.title(f"{self.property} distribution")
        plt.xlabel(f"{self.property}")

    def print_stat(self):
        print(
            f"Min = {np.min(self.filtered_data)*1e9:0.2e} eV, Max = {np.max(en_data)*1e9:0.2e} eV"
        )


class CascadeAnalysis:
    def __init__(self, cascade_data):
        import particle
        self.data = cascade_data
        self.particles = cascade_data.get_particles()
        self.pack_data()
        self._is_height = False
        
        self.all_pdgs = {
            int(p.pdgid): f"${p.latex_name}$" for p in particle.Particle.findall()
        }

    def print_stats(self):
        print(f"Number of final particles = {self.data.num_final_particles}")
        print(f"Number of events = {self.data.num_events}:")

        print(
            f" interactions = {self.data.num_interactions},"
            f" decays = {self.data.num_decays}"
        )

    def pack_data(self):
        self.pid_data = []
        self.energy_data = []
        self.xdepth_data = []

        for particle in self.particles:
            self.pid_data.append(particle.pid)
            self.energy_data.append(particle.energy)
            self.xdepth_data.append(particle.xdepth)
        
    def calc_height(self):
        xconv = XdepthConversion()
        xconv.set_length_unit("km")
        
        self.height_data = []
        for i, ppid in enumerate(self.pid_data):
            self.height_data.append(xconv.convert_x2h(self.xdepth_data[i]))
            
        self._is_height = True

    def plot_pid(self, from_=None, to_=None):

        pid_dist = dict()

        for pid in self.pid_data:
            pid_dist[pid] = pid_dist.get(pid, 0) + 1

        pid_dist = dict(
            sorted(pid_dist.items(), key=lambda item: item[1], reverse=True)
        )
        print(pid_dist)

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
            energy_data = self.energy_data

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