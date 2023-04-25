from parsing import BlockDataFile
import numpy as np
from merge_bins import merge_bins


def fluka_data(data_file):
    energy_ranges = {"E1": 0, "E2": 1, "E3": 2}
    xdepth_ranges = {"143": 0, "647": 1, "1033": 2}
    block_file = BlockDataFile(data_file)

    data_dict = {}
    for i in range(block_file.num_blocks()):
        # print(block_file.block_header(i).split()[4].split("-"))
        erange, xdepth = block_file.block_header(i).split()[4].split("-")

        ang1 = block_file.block_data(i, col=0)
        ang2 = block_file.block_data(i, col=1)

        ang_bins = np.concatenate([ang1, [ang2[-1]]])
        val = block_file.block_data(i, col=2)
        error = block_file.block_data(i, col=3)

        ixdepth = xdepth_ranges[xdepth]
        ienergy = energy_ranges[erange]
        if data_dict.get(ixdepth, None) is None:
            data_dict[ixdepth] = {}

        data_dict[ixdepth][ienergy] = [ang_bins, val, error]

    return data_dict


def hist_one_value(ang_bins, dist, error, bin_merge_level):
    omega_bins = 2 * np.pi * (1 - np.cos(ang_bins))
    omega_widths = omega_bins[1:] - omega_bins[:-1]

    hist_dist, new_ang_bins = merge_bins(dist * omega_widths, ang_bins, bin_merge_level)
    hist_error, _ = merge_bins(error * omega_widths, ang_bins, bin_merge_level)

    return [new_ang_bins, hist_dist, hist_error]


def fluka_files(label):
    from pathlib import Path

    base_dir = Path("/hetghome/antonpr/xmax_sigma/flincpy/scripts")
    data_dir = base_dir / "fluka_comparison" / "new_data"

    if label == "develop":
        data_file = data_dir / "atmop100gev_yld_tab.lis"
    elif label == "current":
        data_file = data_dir / "atmop100gev20212_yld_tab.lis"
    else:
        raise ValueError(
            f"No such data: {label}. Possible values: 'develop' or 'current'"
        )

    return data_file


def fluka_histogram(label, bin_merge_level=0):
    fl_data = fluka_data(fluka_files(label))

    data_dict = {}
    en_labels = ["1.0-1.3", "2.0-2.5", "4.0-5.0"]
    xdepth_labels = ["143", "647", "1033"]

    for ixdepth, xd_dict in fl_data.items():
        if data_dict.get(ixdepth, None) is None:
            data_dict[ixdepth] = {}

        for ienergy, val_list in xd_dict.items():
            vals = hist_one_value(*val_list, bin_merge_level)
            # xdepth label
            vals.append(xdepth_labels[ixdepth])
            # energy range label
            vals.append(en_labels[ienergy])
            data_dict[ixdepth][ienergy] = vals

    return data_dict


def fluka_files_en(label):
    from pathlib import Path

    base_dir = Path("/hetghome/antonpr/xmax_sigma/flincpy/scripts")
    data_dir = base_dir / "fluka_comparison" / "new_data"
    if label == "develop":
        data_file = data_dir / "atmop100gev_lgyld_tab.lis"
    elif label == "current":
        data_file = data_dir / "atmop100gev20212_lgyld_tab.lis"
    else:
        raise ValueError(
            f"No such data: {label}. Possible values: 'develop' or 'current'"
        )

    return data_file


def fluka_histogram_en(label, bin_merge_level = 0):

    data_file = fluka_files_en(label)
    block_file = BlockDataFile(data_file)
    
    xdepth_labels = ["143", "647", "1033"]
    hist_data = {}
    for i in range(block_file.num_blocks()):
        
        en1 = block_file.block_data(i, col = 0)
        en2 = block_file.block_data(i, col = 1)
        
        energy_bins = np.concatenate([en1, [en2[-1]]])
        flux = block_file.block_data(i, col=2)
        error = block_file.block_data(i, col=3)
        
        energy_width = energy_bins[1:] - energy_bins[:-1]
        
        hist_flux, new_energy_bins = merge_bins(flux * energy_width, energy_bins, bin_merge_level)
        hist_error, _ = merge_bins(error * energy_width, energy_bins, bin_merge_level)
        
        
        hist_data[i] = [new_energy_bins, hist_flux, hist_error, xdepth_labels[i]]

    return hist_data


if __name__ == "__main__":
    # res = fluka_histogram("develop", bin_merge_level=4)
    fluka_histogram_en("current")
