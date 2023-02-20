


class EventDispatcher:
    def __init__(self):
        pass
    
    
    def distribute_pstack(self, pstack):
        
        # There are one pre-interaction filter:
        # Check if there are particles with pdgs from final_pdg
        # Put them in the finals
        
        # Check if there are particles with pdgs from decay_pdg
        # Put them in the decaying
        
        # Check if there are particles with energies < threshold_energy
        # Put them in the finals
        
        # Calculate next interaction and decay xdepth for
        
        
        # Check if there are particles with next event (decay or interaction)
        # on the ground (or below)
        # Put them in the finals
        
        # Compare xdepth_decay and xdepth_inter
        
        
        
        # From current stack:
        # Append to finals
        # Append to decaying
        # Put to interacting