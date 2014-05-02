# --- Inherent python/system level imports

# --- OpenMDAO imports
from openmdao.main.api import Component
from openmdao.lib.datatypes.api import Int, Str, Event
     
class DOECounterComp(Component): 
       
    # --- Initialize counter variable ---
    config = Str('', iotype = 'out', desc = 'current design configuration number')
    
    def __init__(self, *args, **kwargs):
        # ---------------------------------------------
        # --- Constructor for the counter component ---
        # ---------------------------------------------
        super(DOECounterComp, self).__init__(*args, **kwargs)
        
        self.case_list = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        
        start_case = 2
        end_case = 4
        
        # --- Edit this line to define DOE case number
        self.case = 2

        #[self.case_list.append(case_number) for case_number in range(start_case, end_case)]
        
        self.force_execute = True    
        
    def execute(self): 
        
        self.config = str(self.case_list[self.case])
        self.case += 1
            
class SubCounterComp(Component): 
       
    # --- Initialize counter variable ---
    config = Str('', iotype = 'out', desc ='full ID string')    
    case = Str('', iotype = 'in', desc = 'current DOE case number')
    power_level = Int(100, iotype = 'in', desc = '% thrust setting')
    
    reset_iteration = Event()

    def _reset_iteration_fired(self):
        self._iteration = 0
        
    def __init__(self, *args, **kwargs):
        # ---------------------------------------------
        # --- Constructor for the counter component ---
        # ---------------------------------------------
        super(SubCounterComp, self).__init__(*args, **kwargs)

        # --- Edit this line to override default numbering
        # --- If restarting, should equal the subiteration number of the last successful case + 1
        # --- Needs to be 0 if starting a new case
        self._iteration = 7
        self.force_execute = True
        
    def execute(self): 

        self.config = '%s_%s_%s'%(self.case, str(self._iteration), str(self.power_level))
        
        self._iteration += 1

        print '--------------------------------------------------'         
        print 'Starting Simulation ', self.config
        print '--------------------------------------------------'
        
if __name__ == "__main__":
    
    # --- Default test case ---      
   c = SubCounterComp()
