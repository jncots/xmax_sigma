import numpy as np


class ParticleArray:
    
    data_attributes = ["pid", 
                       "energy", 
                       "xdepth",
                       "generation_num",
                       "dxdepth_decay",
                       "dxdepth_inter",
                       "production_code",
                       "final_code"]
    
    def __init__(self, size = None):
        
        # Allocate memory
        if size is not None:
            self._allocate(size)
        else:
        # Use is as view    
            self.pid = None
            self.energy = None
            self.xdepth = None
            self.generation_num = None
            self.dxdepth_decay = None
            self.dxdepth_inter = None
            self.production_code = None
            self.final_code = None
            self.data_slice = None
            self._len = None
            self.data = None
            
        
    
    def _allocate(self, size):
        self.pid = np.empty(size, dtype = np.int64)
        self.energy = np.empty(size)
        self.xdepth = np.empty(size)
        self.generation_num = np.empty(size, dtype = np.int64)
        self.dxdepth_decay = np.empty(size)
        self.dxdepth_inter = np.empty(size)
        self.production_code = np.empty(size, dtype = np.int64)
        self.final_code = np.empty(size, dtype = np.int64)
        self.data_slice = np.empty(size, dtype = np.int64)
        self.data_slice.fill(0)
        self.data = self
        self._len = 0
            
    
    def __len__(self):
        return self._len
    
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
        self.data_slice[dst_slice].fill(1)
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
        self.data_slice[dst_slice] = 1
        return dst_slice  


    def view(self, view_slice = None):
        view_stack = ParticleArray()
        
        if view_slice is None:
            view_slice = slice(0, None)
        
        for attr in self.data_attributes:
            value = getattr(self, attr)
            setattr(view_stack, attr, value[view_slice])
        
        view_stack._len = len(view_stack.pid)
        view_stack.data_slice = np.copy(self.data_slice)
        view_stack.data_slice[np.where(view_stack.data_slice == 0)] = 333
        view_stack.data = self
        return view_stack
        
    
    def copy(self, *, src_slice = None, dst_slice = None, size = None):
        
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
        
        copy_stack.data_slice[dst_slice] = np.copy(self.data_slice[src_slice])
        return copy_stack

    def clear(self, size = None):     
        if size is None:
            self._len = 0
        else:
            self._len -= size
            if self._len < 0:
                self._len = 0
        return

    def pop(self, size = None):
        
        if size is None or self._len < size:
            pop_slice = slice(0, None)
        else:
            pop_slice = slice(self._len - size, self._len)    
            
        popped_stack = self.copy(src_slice=pop_slice)
        self.clear(size)
        return popped_stack
        
                    

    # def pop(self, size = None):
    #     if size is None:
    #         self._len = 0
    #     else:
                

    # def __getitem__(self, arg):
    #     pass    
    
    # def append_one(self, pid, energy, xdepth):
    #     dst_slice = self._len
    #     self.pid[dst_slice] = pid
    #     self.energy[dst_slice] = energy
    #     self.xdepth[dst_slice] = xdepth
    #     self._len += 1
        
        
    # def append(self, pid, energy, xdepth, src_slice = None):
    #     dst_start = self._len
    #     dst_end = self._len + len(pid[src_slice])
    #     dst_slice = slice(dst_start, dst_end)
        
    #     self.pid[dst_slice] = np.copy(pid[src_slice])
    #     self.energy[dst_slice] = np.copy(energy[src_slice])
    #     self.xdepth[dst_slice] = np.copy(xdepth[src_slice])
    #     self._len = dst_end
                

pstack = ParticleArray(10)

print(pstack.push(pid = 2212))
print(pstack.push(pid = 2212, energy = 789, xdepth = 0))
print(pstack.push(pid = np.array([2212, 342, 777]), energy = np.array([20, 20, 20])))

print(pstack.pid)

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

# pcopy = pstack.copy(src_slice=np.where((pstack.energy > 1) & (pstack.energy < 1000)), dst_slice = slice(1, 5))
# print(pstack.energy)

# print(pcopy.data_slice)
# print(len(pcopy))
# print(pcopy.energy)
# a = np.empty(10, dtype = np.int64)
# b = np.array([0, 0, 4, 5, 0, 0], dtype = np.int64)

# llen = 8
# sl = slice(llen, None)
# a[sl] = np.copy(b[2:4])
# # a[-3] = np.array(87)

# # print(a)
# # print(len(a[-3:-1]))

# # a[11] = 5

# print(a[None])


        