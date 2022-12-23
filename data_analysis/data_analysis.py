import matplotlib.pyplot as plt
import numpy as np
from cascade.xdepth_conversion import XdepthConversion
from copy import copy


class FilterParticleData:
    property_list = ['pid', 'energy', 'xdepth', 'production_mode', 'generation_number']
    
    def __init__(self, property_):
        if property_ not in self.property_list:
            raise ValueError(f"Property \"{property_}\" is unknown")
        self.property = property_
        self.conversion = None
        self.condition = None
        self.filtered_data = None
        self._set_filter()
        
    
    def set_property(self, property_):
        if property_ not in self.property_list:
            raise ValueError(f"Property \"{property_}\" is unknown")
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
        return [getattr(particle, self.property) for particle in data if self.condition(particle)]
    
    def _filter2(self, data):
        return [self.conversion(getattr(particle, self.property)) for particle in data]
    
    def _filter3(self, data):
        return [self.conversion(getattr(particle, self.property)) for particle in data if self.condition(particle)]
    
    def filter(self, data):
        self.filtered_data = self.ffilter(data)
        return self
    
    def result_data(self):
        return self.filtered_data
    
    def plot_data(self, bins_number = 100):
        plt.stairs(*np.histogram(self.filtered_data, bins_number))
        plt.title(f"{self.property} distribution")
        plt.xlabel(f"{self.property}")
        
        
    def print_stat(self):
        print(f"Min = {np.min(self.filtered_data)*1e9:0.2e} eV, Max = {np.max(en_data)*1e9:0.2e} eV")
                

            
class CascadeAnalysis:
    def __init__(self, cascade_data):
        self.data = cascade_data
        self.particles = cascade_data.get_particles()

    def print_stats(self):
        print(f"Number of final particles = {self.data.num_final_particles}")
        print(f"Number of events = {self.data.num_events}:")

        print(
            f" interactions = {self.data.num_interactions},"
            f" decays = {self.data.num_decays}"
        )
    
    
    def plot_xdepth(self):
        pfilter = FilterParticleData("xdepth")
        pfilter.filter(self.particles)
        pfilter.plot_data()
    
    def plot_energy(self):
        
        # def to_log(x):
        #     np.log10(x)
        
        def pid_cond(x):
            return x.pid == 22
        
        pfilter = FilterParticleData("energy")
        en_data = pfilter.filter(self.particles).set_condition(pid_cond)
        en_data = pfilter.result_data()
        
        gr, cnt = np.histogram(np.log10(en_data), 100)
        plt.semilogx()
        plt.step(10 ** cnt[:-1] * 1e9, gr)
        plt.title("Energy distribution")
        plt.xlabel("Energy, eV")
        # plt.loglog(basex = 10)
        # plt.stairs(gr,cnt)    
  
    # def plot_energy(self):
        
        
    #     def to_log(x):
    #         np.log10(x)
            
    #     pfilter = FilterParticleData("energy")
    #     pfilter.filter(self.particles)

    #     etot = 0
    #     for ee in pfilter.result_data():
    #         etot += ee
        
    #     print(f"Energy conservation = {abs(etot - 1e9)/1e9}")    

    #     pfilter.set_conversion(to_log).filter(self.particles)
    #     pfilter.plot_data(100)
        
        
    #     print(f"Energy conservation = {abs(etot - 1e9)/1e9}")

    #     gr, cnt = np.histogram(np.log10(en_data), 100)
    #     plt.semilogx()
    #     plt.step(10 ** cnt[:-1] * 1e9, gr)
    #     plt.title("Energy distribution")
    #     plt.xlabel("Energy, eV")
    #     # plt.loglog(basex = 10)
    #     # plt.stairs(gr,cnt)
    #     # print(10**cnt)
    #     np.min(en_data) * 1e9 / 1e6
    #     print(f"Min = {np.min(en_data)*1e9:0.2e} eV, Max = {np.max(en_data)*1e9:0.2e} eV")


    # def plot_xdepth(self):

    #     xdepth_values = []
    #     for particle in self.particles:
    #         xdepth_values.append(particle.xdepth)
    #         # if xdepth_values[i] > 800:
    #         #     xdepth_values[i] = 800
    #         # xm_data.append(xdepth)
    #         # heightx.append(xconv.convert_x2h(xdepth))

    #     print(f"Min = {np.min(xdepth_values):.2f}, Max = {np.max(xdepth_values):.2f}")

    #     plt.stairs(*np.histogram(xdepth_values, 100))
    #     plt.title("Xdepth distribution")
    #     plt.xlabel("Xdepth")


    # def plot_height_distribution(self):
    #     xconv = XdepthConversion()
    #     xconv.set_length_unit("km")
        
    #     for i in :
    #         xdepth_values[i] = copy(final_prt[i].xdepth)
    #         height_values[i] = xconv.convert_x2h(xdepth_values[i])
    

    # print(f"Min = {np.min(height_values):.2f} km, Max = {np.max(height_values):.2f} km")

    # plt.stairs(*np.histogram(height_values, 100))
    # plt.title("Height distribution")
    # plt.xlabel("Height, km")           