import numpy as np
from MCEq.geometry.density_profiles import CorsikaAtmosphere

class XdepthConversion:
    length_units = {"cm": 1, "m": 1e2, "km": 1e5}

    def __init__(self):
        self.cka_obj = CorsikaAtmosphere("SouthPole", "December")
        self.cka_obj.set_theta(0.0)
        self.length_unit = self.length_units["cm"]

    def set_theta(self, theta):
        self.cka_obj.set_theta(theta)

    def set_length_unit(self, unit):
        self.length_unit = self.length_units[unit]

    def convert_x2h(self, x):
        """Convert xdepth to height
        """
        return self.cka_obj.X2h(x) / self.length_unit

    def get_length(self, x1, x2):
        """Return a length between atmospheric depth x1 and x2
        x1 < x2. Because X2h coverts only to height, zenith angle 
        is also taken into account.

        Args:
            x1 (float): _description_
            x2 (float): _description_

        Returns:
            float: length in a set length units
        """
        return (self.cka_obj.X2h(x2) - self.cka_obj.X2h(x1)) / (
            np.cos(self.cka_obj.thrad) * self.length_unit
        )
        
    def get_delta_xdepth(self, x1, length):
        """Returns delta x = x2 - x1 for the path starting at point
        x1 and run `length` cm
        """
        h1 = self.cka_obj.X2h(x1)
        h2 = h1 - length * np.cos(self.cka_obj.thrad)
        x2 = self.cka_obj.h2X(h2)
        return x2 - x1