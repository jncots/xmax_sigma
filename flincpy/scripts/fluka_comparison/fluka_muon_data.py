from pathlib import Path
import json
import numpy as np



def parse_data_file(data_file):
    """Returns a list of tuples (header, list_of_rows_with_numbers)
    """
    
    prev_block_type = "start"    
    parsed_data = []
    
    with open(data_file) as file:
        for line in file:
            line = line.strip()
            
            # Determine type of data block
            if line:
                if line.startswith("#"):
                    block_type = "header"
                else:
                    block_type = "data"
            else:
                block_type = "separator"  
            
            # If header starts or continues
            if block_type == "header":
                if prev_block_type == "header":
                    header += f"\n{line}"
                else:
                    if prev_block_type == "separator":
                        parsed_data.append((header, data))
                    header = line
                    data = []     
                
            # If block with numerical data
            if block_type == "data":
                line_record = []
                for number in line.split():
                    line_record.append(float(number))
                data.append(line_record)
                        
            prev_block_type = block_type
        
        if prev_block_type == "data":
            parsed_data.append((header, data))
    
    return parsed_data


def parse_xdepth(header):
    splitted = header.split()
    useful_info = header.split()[4]
    xdepth = float(useful_info.split("-")[1])
    return xdepth
    

def fluka_en_data(file_name = "develop_version"):
    if file_name == "develop_version":
        data_file = Path(__file__).parent/"new_data"/"atmop100gev_lgyld_tab.lis"
        # data_file = Path(__file__).parent/"new_data"/"atmop100gev_yld_tab.lis"
    elif file_name == "current_version":
        data_file = Path(__file__).parent/"new_data"/"atmop100gev20212_lgyld_tab.lis"
    else:
        raise ValueError("This should not happen")

    parsed_data = parse_data_file(data_file)
    
    converted_data = []
    for data_block in parsed_data:
        xdepth = parse_xdepth(data_block[0])
        data = np.array(data_block[1])
        en1 = data[:,0]
        en2 = data[:,1]
        flux = data[:,2]
        error = data[:,3]
        
        energy_bins = np.append(en1, en2[-1])
        
        delta_en = energy_bins[1:] - energy_bins[:-1]
        flux1 = flux * delta_en
        error1 = error * delta_en
        converted_data.append((energy_bins, flux1, error1, xdepth))
    
    return converted_data


def fluka_data(file_name = "old"):
    
    if file_name == "old":
        data_file = Path(__file__).parent/"atmop100gev_yld_tab.lis"
    elif file_name == "devel":
        data_file = Path(__file__).parent/"new_data"/"atmop100gev_yld_tab.lis"
    elif file_name == "current":
        data_file = Path(__file__).parent/"new_data"/"atmop100gev20212_yld_tab.lis"
    else:
        raise ValueError("This should not happen")

    new_block = False

    data_dict = dict()
    dinfo =""
    data = []

    record = False
    with open(data_file) as file:
        for line in file:
            
            if record:
                # print(f"Dinfo = {dinfo.split()}, {dinfo.split()[4]}")
                data_dict[dinfo.split()[4]] = (dinfo, data)
                dinfo =""
                data = []
                record = False
                
            if not line.strip():
                if dinfo.strip():
                    record = True
                continue
            
            if line.startswith(" #"):
                dinfo += line
                new_block = True
            else:
                line_data = []
                for s in line.split():
                    line_data.append(float(s))
                data.append(line_data)    

        else:
            # print(f"Dinfo = {dinfo.split()}, {dinfo.split()[4]}")
            data_dict[dinfo.split()[4]] = (dinfo, data)
            dinfo =""
            data = []
            record = False
            
    data_dict1 = dict()
    for key, item in data_dict.items():
        om1 = []
        om2 = []
        val = []
        err = []
        for line_vals in item[1]:
            om1.append(line_vals[0])
            om2.append(line_vals[1])
            val.append(line_vals[2])
            err.append(line_vals[3])
        
        correction = 1
        if key.startswith("E1"):
            erange = [1, 1.3]
            if file_name == "old":
                correction = 0.3**2
        elif key.startswith("E2"):
            erange = [2, 2.5]
            if file_name == "old":
                correction = 0.5**2
        elif key.startswith("E3"):
            erange = [4, 5]
            if file_name == "old":
                correction = 1
        else:
            raise ValueError("Should not happen")    
                  
        #1-1.3, 2-2.5, 4-5 GeV,
        data_dict1[key] = (item[0], np.array(om1), np.array(om2), 
                           np.array(val)*correction, np.array(err), np.array(erange, dtype=np.float64))
    
    return data_dict1
   

def fluka_dists():
    data_dist = fluka_data()
    ang_dist = dict()
    en_dist = dict()
    
    for key in data_dist:
        xdepth = key.split("-")[1]
        data1 = data_dist[key]
        
        # Omega bins
        omega_bins = np.append(data1[1], data1[2][-1])
        omega_bins = 2 * np.pi * (1 - np.cos(omega_bins))
        d_omega = omega_bins[1:] - omega_bins[:-1]
        
        # Energy bins
        energy_bins = data1[5]
        d_energy = energy_bins[1:] - energy_bins[:-1]
        
        # Distribution
        dN_dOmdE = data1[3]
        
        # hist_dOmdE = (dN_dOmdE * d_omega) * d_energy[0]
        
        hist_dOmdE = dN_dOmdE
        
        fxdepth = float(xdepth)
        ang_dist.setdefault(xdepth, []).append(
            (hist_dOmdE, omega_bins, energy_bins, fxdepth))
        
        en_dist.setdefault(xdepth, []).append(
            (np.sum(hist_dOmdE), energy_bins, fxdepth))
    
    ang_dist1 = []   
    for xdepth in ang_dist:
        ang_dist1.append(ang_dist[xdepth])
    ang_dist = ang_dist1    
    
    en_dist1 = []    
    for xdepth in en_dist:
        en_dist1.append(en_dist[xdepth])
    en_dist = en_dist1   
        
    en_hists = []
    for en_dist_x in en_dist:
        ehist = np.array([en_dist_x[0][0], 0, en_dist_x[1][0], 0, en_dist_x[2][0]])
        en_bins = np.array([en_dist_x[0][1][0], en_dist_x[0][1][1], 
                            en_dist_x[1][1][0], en_dist_x[1][1][1],
                            en_dist_x[2][1][0], en_dist_x[2][1][1]])
        en_hists.append((ehist, en_bins, en_dist_x[0][2]))
    
    return en_hists, en_dist, ang_dist


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
    
    
        
    

def fluka_original_dists(data_dist = fluka_data(), binmerging_level = 3):
    ang_dist = dict()
    en_dist = dict()
    
    for key in data_dist:
        xdepth = key.split("-")[1]
        data1 = data_dist[key]
        
        # Omega bins
        theta_bins = np.append(data1[1], data1[2][-1])
        omega_bins = 2 * np.pi * (1 - np.cos(theta_bins))
        d_omega = omega_bins[1:] - omega_bins[:-1]
        
        # Energy bins
        energy_bins = data1[5]
        d_energy = energy_bins[1:] - energy_bins[:-1]
        
        # Distribution
        dN_dOmdE = data1[3]*d_omega
        
        hist_dOmdE, theta_bins = merge_bins(dN_dOmdE, theta_bins, binmerging_level)
        # hist_dOmdE = hist_dOmdE/(omega_bins[1:] - omega_bins[:-1])
        
        omega_bins = 2 * np.pi * (1 - np.cos(theta_bins))
        d_omega =  omega_bins[1:] - omega_bins[:-1]
        # hist_dOmdE = (hist_dOmdE * d_omega) * d_energy[0]
        
        fxdepth = float(xdepth)
        ang_dist.setdefault(xdepth, []).append(
            (hist_dOmdE, theta_bins, energy_bins, fxdepth))
        
        en_dist.setdefault(xdepth, []).append(
            (np.sum(hist_dOmdE), energy_bins, fxdepth))
    
    ang_dist1 = []   
    for xdepth in ang_dist:
        ang_dist1.append(ang_dist[xdepth])
    ang_dist = ang_dist1    
    
    en_dist1 = []    
    for xdepth in en_dist:
        en_dist1.append(en_dist[xdepth])
    en_dist = en_dist1   
        
    en_hists = []
    for en_dist_x in en_dist:
        ehist = np.array([en_dist_x[0][0], 0, en_dist_x[1][0], 0, en_dist_x[2][0]])
        en_bins = np.array([en_dist_x[0][1][0], en_dist_x[0][1][1], 
                            en_dist_x[1][1][0], en_dist_x[1][1][1],
                            en_dist_x[2][1][0], en_dist_x[2][1][1]])
        en_hists.append((ehist, en_bins, en_dist_x[0][2]))
    
    return en_hists, en_dist, ang_dist      


if __name__ == "__main__":
    en_hist, en_dist, ang_dist = fluka_dists()
    print(en_hist)
    
       
# ang_data = return_data()

# for key in ang_data:
#     print(ang_data[key][1])
       