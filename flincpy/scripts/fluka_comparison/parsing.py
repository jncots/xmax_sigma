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
        
        for block in self.data:
            block[1] = np.array(block[1])
    
    def num_blocks(self):
        return len(self.data)
        
    def block(self, index):
        return self.data[index]
    
    def block_header(self, index):
        return self.data[index][0]
    
    def block_data(self, index):
        return self.data[index][1]
            


if __name__ == "__main__":
    
    from pathlib import Path
    
    base_dir = Path("/hetghome/antonpr/xmax_sigma/flincpy/scripts")
    data_dir = base_dir / "fluka_comparison"
    
    data_file = data_dir / "atmop100gev_yld_tab.lis"
    
    block_file = BlockDataFile(data_file)
    
    print(block_file.num_blocks())
    print(block_file.block_data(0)[0, 1])
    