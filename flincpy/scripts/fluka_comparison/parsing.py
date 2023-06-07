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


class BlockDataFile:
    def __init__(self, filename):
        self.data = parse_data_file(filename)
        
        np_data = []
        for block in self.data:
            np_data.append(np.array(block[1]))
            
        self.np_data = np_data
    
    def num_blocks(self):
        return len(self.data)
        
    def block(self, index):
        return self.data[index]
    
    def block_header(self, index):
        return self.data[index][0]
    
    def block_shape(self, index):
        return self.np_data[index].shape
    
    def block_data(self, index, col = None, row = None):
        if row is not None:
            if col is not None:
                return self.np_data[index][row, col]
            else:
                return self.np_data[index][row, :]
        else:
            if col is not None:
                return self.np_data[index][:, col]
            else:
                return self.np_data[index]
            



if __name__ == "__main__":
    
    from pathlib import Path
    base_dir = Path("/hetghome/antonpr/xmax_sigma/flincpy/scripts")
    data_dir = base_dir / "fluka_comparison" / "alfredo_muons"
    # data_file = data_dir / "atmop100gev_yld_tab.lis"
    data_file = data_dir / "mutest_anueyld_tab.lis"
    
    block_file = BlockDataFile(data_file)
    
    print(block_file.num_blocks())
    for i in range(block_file.num_blocks()):
        print("\n--------")
        print(block_file.block_header(i))
    
    print(block_file.block_header(2))
    print(block_file.block_shape(2))
    print(block_file.block_data(2, col = 2)[0:5])
    