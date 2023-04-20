import numpy as np

def merge_bins(hist, bins, ntimes = 1):
    """Merge adjacent bins"""
    for _ in range(ntimes):
        new_hist = (hist[:-1] + hist[1:])[::2]
        new_bins = bins[::2]
        if (hist.size % 2) > 0:
            new_hist = np.append(new_hist, hist[-1])
            new_bins = np.append(new_bins, bins[-1])
        
        hist = new_hist
        bins = new_bins
        
    return hist, bins