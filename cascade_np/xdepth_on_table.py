from xdepth_conversion import XdepthConversion
import numpy as np


class XdepthOnTable:
    def __init__(self, *, xconv = XdepthConversion(), theta_deg=0, npoints=1000):
        xconv.set_theta(theta_deg)
        # Make height and length vector
        self.height = np.linspace(0, xconv.get_max_height(), npoints, dtype="float64")
        self.length = self.height / np.cos(theta_deg * np.pi / 180)

        # Calculate xdepth vector
        xdepth_fun = np.frompyfunc(lambda x: np.float64(xconv.convert_h2x(x)), 1, 1)
        self.xdepth = xdepth_fun(self.height).astype("float64")

        # Revert vectors for interpolation
        self.rev_height = np.copy(self.height[::-1])
        self.rev_length = np.copy(self.length[::-1])
        self.rev_xdepth = np.copy(self.xdepth[::-1])

    def convert_x2h(self, xdepth_vec):
        return np.interp(xdepth_vec, self.rev_xdepth, self.rev_height)

    def convert_h2x(self, height_vec):
        return np.interp(height_vec, self.height, self.xdepth)

    def add_len2x(self, xdepth_vec, length_vec):
        """Returns final xdepth for initial xdepth (g/cm^2) and
        length (in cm)

        Args:
            xdepth_vec (np.array): initial xdepth
            length_vec (np.array): delta length
        """
        print("xdepth_vec", xdepth_vec)
        print("length0", np.interp(xdepth_vec, self.rev_xdepth, self.rev_length))
        print("length_vec", length_vec)
        
        length = np.interp(xdepth_vec, self.rev_xdepth, self.rev_length) - length_vec
        print("length", length)
        print("length", np.interp(length, self.length, self.xdepth))
        return np.interp(length, self.length, self.xdepth)


if __name__ == "__main__":
    xconv = XdepthOnTable()

    nn = 1000000

    xdepth = np.array(np.zeros(nn), dtype="float64")
    dlen = np.random.rand(nn) * 1e7
    # print(dlen)
    import time

    start = time.process_time()
    xconv.add_len2x(xdepth, dlen)
    print((time.process_time() - start) / nn)
# np.zeros(10000)

# np.random.rand()*1e7

# xdepth = np.array(np.zeros(10000),dtype='float64')
# dlen = np.array([1, 2, 3, 4],dtype='float64') * 1e8
# print(dlen)

# print(xconv.add_len2x(xdepth, dlen))

# print(np.random.rand(5)*1e7)


# from xdepth_conversion import XdepthConversion
# import matplotlib.pyplot as plt
# import numpy as np


# theta = 80
# xconv = XdepthConversion(theta)
# xconv.set_length_unit("km")

# height = np.linspace(0, xconv.get_max_height(), 1000, dtype='float64')
# length = height/np.cos(80 * np.pi/180)
# xdepth_fun = np.frompyfunc(lambda x : np.float64(xconv.convert_h2x(x)), 1, 1)
# xdepth = xdepth_fun(height).astype('float64')
# plt.loglog(height, xdepth)
# plt.grid()
# print(xconv.convert_h2x(0))

# # np.interp(np.array([1, 2, 3],dtype='float64'), height, xdepth)
# print(np.interp(np.array([1, 2, 3, 4],dtype='float64'), xdepth[::-1], height[::-1]))
# print(np.interp(np.array([0, 0, 0, 0],dtype='float64'), xdepth[::-1], length[::-1]))

# llen = np.interp(np.array([0, 0, 0, 0],dtype='float64'), xdepth[::-1], length[::-1])
# dlen = np.array([100, 200, 300, 400], dtype = np.float64)
# print(np.interp(llen - dlen, length,  xdepth))
# # print(llen - dlen)