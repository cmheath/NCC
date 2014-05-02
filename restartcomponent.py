''' OpenMDAO Restart Component '''
''' DOE Restart Component reads in a previously created DB file and sets design parameters to match case values ''' 

# --- OpenMDAO imports
from openmdao.main.api import Component
from openmdao.lib.datatypes.api import Float, Int, Str
from openmdao.lib.casehandlers.api import DBCaseRecorder, case_db_to_dict

DOE_OUT_DB = 'DOE_Output.db'

class RestartComp(Component):
    ''' OpenMDAO component for restarting cases from a DOE '''

    # -----------------------------------------------------
    # --- Initialize Input Design Parameters and Ranges ---
    # -----------------------------------------------------
    venturi_angle = Float(30.0, low = 20.0, high = 35.0, iotype ='out', desc = 'angle of converging/diverging venturi section', units = 'deg')
    vane_height = Float(0.455205, low = 0.317293, high = 0.540691, iotype ='out', desc = 'height of each helical vane blade')
    pilot_recession = Float(0.25, low = 0.0, high = 0.75, iotype ='out', desc = 'pilot dome recession depth', units = 'inch')      
    config = Str('0', iotype = 'in', desc = 'current design configuration number')
    
    def __init__(self):
        super(RestartComp, self).__init__()
        
        self.force_execute = True
    
    def execute(self):
    
        # --- Original Case
        # venturi_angle_list = [30.0]
        # vane_height_list = [0.455205]
        # pilot_recession_list = [0.25]

        pilot_recession_list =  [0.4453125, 0.1640625, 0.4921875, 0.1171875, 0.000000, 0.6796875, 0.7265625, 0.3046875, 0.3984375, 0.2109375, 0.2578125, 0.5859375, 0.3515625, 0.0703125, 0.5390625, 0.6328125]
        venturi_angle_list = [27.03125, 24.21875, 29.84375, 23.28125, 22.34375, 28.90625, 30.78125, 26.09375, 34.53125, 33.59375, 27.96875, 25.15625, 20.46875, 31.71875, 32.65625, 21.40625]
        vane_height_list = [0.5337098125, 0.4220108125, 0.4918226875, 0.3521989375, 0.4499355625, 0.3382365625, 0.3940860625, 0.47786031250000005, 0.4638979375, 0.38012368750000003, 0.4359731875, 0.3661613125, 0.5197474375, 0.5057850625, 0.3242741875, 0.4080484375]
        
        # --- Note: vane angles = [40.4255455865, 47.1151139333, 42.7138936666, 52.2504676998, 45.2665298066, 53.3656344731, 49.0811594199, 43.5354070399, 44.3862857533, 50.1082301732, 46.1761391999, 51.1646662665, 41.1589629399, 41.9217456332, 54.5101665864, 48.0834540066]      
        
        self.venturi_angle = venturi_angle_list[int(self.config)]
        self.vane_height = vane_height_list[int(self.config)]
        self.pilot_recession = pilot_recession_list[int(self.config)]       
        
        print '--------------------------------------------------'
        print 'Config Number    ', self.config
        print '--------------------------------------------------'
        
if __name__ == "__main__":
    
    # -------------------------
    # --- Default Test Case ---
    # ------------------------- 
    Restart_Comp = RestartComp()
    Restart_Comp.run()