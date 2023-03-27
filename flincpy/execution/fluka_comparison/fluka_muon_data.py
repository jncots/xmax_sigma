from pathlib import Path
import json
import numpy as np 


def fluka_data():
    data_file = Path(__file__).parent/"atmop100gev_yld_tab.lis"


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
        
        if key.startswith("E1"):
            erange = [1, 1.3]
        elif key.startswith("E2"):
            erange = [2, 2.5]
        elif key.startswith("E3"):
            erange = [4, 5]
        else:
            raise ValueError("Should not happen")    
                  
        #1-1.3, 2-2.5, 4-5 GeV,
        data_dict1[key] = (item[0], np.array(om1), np.array(om2), 
                           np.array(val), np.array(err), np.array(erange, dtype=np.float64))
    
    return data_dict1
   

def fluka_en_dist():
    data_dist = fluka_data()
    en_dist = dict()
    
    for i in data_dist:
        ang_dist = data_dist[i]
        num_dist = ang_dist[3]*(ang_dist[2]-ang_dist[1])*(ang_dist[5][1] - ang_dist[5][0])
        xdepth = i.split("-")[1]
        en_dist.setdefault(xdepth, []).append([np.sum(num_dist), ang_dist[5][0], ang_dist[5][1]])
        
    num_en_dist = dict()
    for key, value in en_dist.items():
        num = np.array([value[0][0], 0, value[1][0], 0, value[2][0]])
        en = np.array([value[0][1], value[0][2], value[1][1], value[1][2], value[2][1], value[2][2]])
        num_en_dist[key] = (en, num)
    
    return num_en_dist       
   
# ang_data = return_data()

# for key in ang_data:
#     print(ang_data[key][1])
       