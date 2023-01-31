import numpy as np
from MCEq.geometry.density_profiles import CorsikaAtmosphere


class XdepthConversion:
    length_units = {"cm": 1, "m": 1e2, "km": 1e5}
    _min_xdepth = 1e-7

    def __init__(self, theta = 0.0):
        self.cka_obj = CorsikaAtmosphere("SouthPole", "December")
        self.cka_obj.set_theta(theta)
        self.length_unit = self.length_units["cm"]

    def set_theta(self, theta):
        self.cka_obj.set_theta(theta)

    def set_length_unit(self, unit):
        self.length_unit = self.length_units[unit]

    def convert_x2h(self, x):
        """Convert xdepth to height"""
        if x < self._min_xdepth:
            x = self._min_xdepth
        return self.cka_obj.X2h(x) / self.length_unit

    def convert_h2x(self, x):
        """Convert xdepth to height"""
        return self.cka_obj.h2X(x * self.length_unit)

    def get_max_xdepth(self):
        return self.convert_h2x(0)
    
    def get_max_height(self):
        return self.convert_x2h(0)

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
        return (self.convert_x2h(x1) - self.convert_x2h(x2)) / (
            np.cos(self.cka_obj.thrad)
        )

    def get_delta_xdepth(self, x1, length):
        """Returns delta x = x2 - x1 for the path starting at point
        x1 and run `length` cm
        """
        h1 = self.convert_x2h(x1)
        h2 = h1 - length * np.cos(self.cka_obj.thrad)
        if h2 < 0:
            return None
        x2 = self.convert_h2x(h2)
        return x2 - x1



if __name__ == "__main__":
    xconv = XdepthConversion(0)
    xconv.set_length_unit("km")
    print(f"Max height = {xconv.get_max_height()}")
    print(f"Max depth = {xconv.get_max_xdepth()}")
    print(f"Height for xdepth = 500: {xconv.convert_x2h(500)}")
    print(xconv.get_delta_xdepth(500, 5.2))

# xconv = XdepthConversion()
# xconv.set_length_unit("km")
# xconv.set_theta(60)
# print(xconv.convert_x2h(0))
# print(xconv.get_max_height())
# print(xconv.get_delta_xdepth(0, 300))
# xconv.set_theta(60)
# # print(xconv.convert_x2h(100))
# print(f"x = 100 = {xconv.convert_x2h(100)} km")
# print(f"x = 700 = {xconv.convert_x2h(700)} km")

# print(f"h = 15km = {xconv.convert_h2x(xconv.convert_x2h(100))}")


# print(xconv.get_length(500, 700))
# print(xconv.get_delta_xdepth(500, xconv.get_length(500, 700)))
# print(xconv.get_delta_xdepth(300, 0.5))
# print(xconv.convert_x2h(300))
# print(xconv.convert_x2h(300 + 311))
