import numpy as np
from enum import Enum


class FilterCode(Enum):
    # xdepth_decay is already calculated
    XD_DECAY_ON = 21
    # xdepth_decay is empty
    XD_DECAY_OFF = 22

class ParticleArray:
    """ParticleArray is data structure to keep particles and their properties

        filter_code: is code to filter entries
        {1: interacting, 2: decaying, 3: final}
    """

    data_attributes = ["pid",
                       "energy",
                       "xdepth",
                       "generation_num",
                       "xdepth_decay",
                       "xdepth_inter",
                       "production_code",
                       "final_code",
                       "filter_code",
                       "valid_code"]

    def __init__(self, size=None):

        # Allocate memory
        if size is not None:
            self._allocate(size)
        else:
            # Use is as view
            self.pid = None
            self.energy = None
            self.xdepth = None
            self.generation_num = None
            self.xdepth_decay = None
            self.xdepth_inter = None
            self.production_code = None
            self.final_code = None
            self.filter_code = None
            self.valid_code = None
            self._len = None
            self.data = None

    def _allocate(self, size):
        self.pid = np.empty(size, dtype=np.int64)
        self.energy = np.empty(size)
        self.xdepth = np.empty(size)
        self.generation_num = np.empty(size, dtype=np.int64)
        self.xdepth_decay = np.empty(size)
        self.xdepth_inter = np.empty(size)
        self.production_code = np.empty(size, dtype=np.int64)
        self.final_code = np.empty(size, dtype=np.int64)
        self.filter_code = np.empty(size, dtype=np.int64)
        self.valid_code = np.zeros(size, dtype=np.int64)
        self.data = self
        self._len = 0

    def __len__(self):
        return self._len
    
    
    def __getitem__(self, key):
        view_stack = ParticleArray()
        for attr in self.data_attributes:
            value = getattr(self, attr)[key]
            setattr(view_stack, attr, value)

        view_stack._len = view_stack.pid.size
        view_stack.data = self
        return view_stack
        
    def __setitem__(self, slice_, other):

        if not isinstance(other, ParticleArray):
            raise ValueError("Argument is not an object of ParticleArray class")
        
        for attr in self.data_attributes:
            other_value = getattr(other, attr)
            self_value = getattr(self, attr)
            self_value[slice_] = other_value
        

    def reserved_size(self):
        return len(self.pid)

    def push(self, **kwargs):

        pid = kwargs.get("pid")
        if pid is None:
            return
        
        # If scalar (one) value is given
        if isinstance(pid, int):
            return self.push_one(**kwargs)
        
        src_slice = kwargs.get("src_slice")
        if src_slice is None:
            src_slice = slice(0, None)
        
        dst_start = self._len
        dst_end = self._len + len(pid[src_slice])
        dst_slice = slice(dst_start, dst_end)
        
        for name, value in kwargs.items():
            data_attr = getattr(self, name, None)
            if data_attr is not None:
                # Here should be an Exception is thrown
                # If you met a problem here, it means
                # that you gave len(pid) > len(other_data array)
                # Arrays should be equal
                data_attr[dst_slice] = np.copy(value[src_slice])
        self._len = dst_end
        self.valid_code[dst_slice].fill(1)
        return dst_slice
    
    def push_one(self, **kwargs):
        pid = kwargs.get("pid")
        if pid is None:
            return
        
        dst_start = self._len
        dst_end = self._len + 1
        dst_slice = slice(dst_start, dst_end)
        
        for name, value in kwargs.items():
            data_attr = getattr(self, name, None)
            if data_attr is not None:
                data_attr[dst_slice] = np.copy(value)
        self._len = dst_end
        self.valid_code[dst_slice] = 1
        return dst_slice


    def copy(self, *, src_slice=None, dst_slice=None, size=None):
        if size is None:
            copy_stack = ParticleArray(self.reserved_size())
        else:
            copy_stack = ParticleArray(size)

        if src_slice is None:
            src_slice = slice(0, None)   
            
        copy_stack._len = 0
        if dst_slice is None:
            copy_stack._len = len(self.pid[src_slice])
            dst_slice = slice(0, copy_stack._len)
        
        for attr in self.data_attributes:
            src_value = getattr(self, attr)[src_slice]
            getattr(copy_stack, attr)[dst_slice] = np.copy(src_value)
        
        copy_stack.valid_code[dst_slice] = np.copy(self.valid_code[src_slice])
        return copy_stack

    def clear(self, size=None):
        if size is None:
            self._len = 0
        else:
            self._len -= size
            if self._len < 0:
                self._len = 0
        return

    def pop(self, size=None):
        if size is None or self._len < size:
            pop_slice = slice(0, None)
        else:
            pop_slice = slice(self._len - size, self._len)
                      
        popped_stack = self.copy(src_slice=pop_slice)
        self.clear(size)
        return popped_stack
    
    def append(self, other):
        
        if not isinstance(other, ParticleArray):
            raise ValueError("argument is not a ParticleArray object")
        
        other_slice =  slice(0, len(other))
        self_slice = slice(len(self), len(self) + len(other))
        
        for attr in other.data_attributes:
            val_other = getattr(other, attr)
            val_self = getattr(self, attr)
            val_self[self_slice] = val_other[other_slice]
            
        self._len = len(self) + len(other)
        return self    
            
    def valid(self):
        return self[0:self._len]  

if __name__ == "__main__":
    # Initialize stack having size 10 reserved
    pstack = ParticleArray(10)

    # Push single values (pid is required)
    print(f"Add single element {pstack.push(pid = 2212)}")
    print(f"Add single element with energy {pstack.push(pid = 2212, energy = 789, xdepth = 0)}")
    # Push arrays
    print("Push array with 3 elements "
          f"{pstack.push(pid = np.array([888, 342, 777]), energy = np.array([20, 20, 20]))}")


    print(pstack.valid().energy)
    
    pstack1 = ParticleArray(100)
    pstack1.append(pstack).append(pstack).append(pstack)
    print(np.where(np.in1d(pstack1.valid().pid, np.array([888])))[0])
    # print(pstack1.valid()[np.where(pstack1.valid().pid in [888, 777])].pid)
    

    # # Check len of stack (number of filled elements)
    # print(f"Size of stack = {len(pstack)}")
    # print(f"Print pid array = {pstack.pid}")
    # print(f"Print data_slice array = {pstack.data_slice}")
    
    # # Pop (remove from stack and return a copy)
    # top_stack = pstack.pop(2)
    # print(f"Size of stack after pop(2) = {len(pstack)}")

    # print(f"Print pid array of pop(2) = {top_stack.pid}")
    # print(f"Print data_slice array of pop(2) = {top_stack.data_slice}")
    # print(f"Size of top_stack = {len(top_stack)}")

# print(len(pstack), pstack.reserved_size())

# print(pstack.pid)
# print(pstack.energy)

# view_p = pstack.view(slice(1, 5))
# print(view_p.pid)

# pstack.pid[4] = 888
# print(view_p.pid, len(view_p))
# print(pstack.pid[np.where(pstack.data_slice>0)])


# print(view_p.data_slice)
# print(view_p.data.view([1, 1, 1]).energy)

# pcopy = pstack.copy(src_slice=np.where((pstack.energy > 1) & (pstack.energy < 1000)),
# dst_slice = slice(1, 5))
# print(pstack.energy)

# print(pcopy.data_slice)
# print(len(pcopy))
# print(pcopy.energy)
