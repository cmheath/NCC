''' OpenMDAO geometry component '''
''' Updates geometry via the SolidWorks API and exports as a STEP file''' 
''' Loads STEP file into cubit and applies a geometry adaptive tetrahedral mesh'''

# --- Inherent python/system level imports
import os
import sys
import shutil

# --- External python library imports (i.e. matplotlib, numpy, scipy)
from tempfile import TemporaryFile
from math import pi

# --- OpenMDAO imports
from openmdao.main.api import FileMetadata
from openmdao.lib.datatypes.api import Float, Int, Str
from openmdao.lib.components.api import ExternalCode
from pathname import Path


class GeometryComp(ExternalCode):
    ''' OpenMDAO component for Modifying Geometry '''
    # -----------------------------------------------------
    # --- File Wrapper for SolidWorks via Visual Basic  ---
    # -----------------------------------------------------
    config = Str('1', iotype ='in', desc ='current design configuration number')
    
    # -----------------------------------------------------
    # --- Initialize Input Design Parameters and Ranges ---
    # -----------------------------------------------------
    venturi_angle = Float(30.0, low = 20.0, high = 35.0, iotype ='in', desc = 'angle of converging/diverging venturi section', units = 'deg')
    vane_height = Float(0.455205, low = 0.317293, high = 0.540691, iotype ='in', desc = 'height of each helical vane blade')
    pilot_recession = Float(0.25, low = 0.0, high = 0.75, iotype ='in', desc = 'pilot dome recession depth', units = 'inch')
    
    # ---------------------------------    
    # --- Fixed Geometry Dimensions ---
    # ---------------------------------      
    injector_dia = Float(0.0, iotype ='out', desc = 'individual injector diameter', units = 'inch')    
    sector_area = Float(0.0, iotype ='out', desc = 'cross-sectional area of 15-deg combustor sector', units = 'inch**2')
    venturi_exit_area = Float(0.0, iotype ='out', desc = 'total exit area of venturis', units = 'inch**2')    
    venturi_throat_area = Float(0.0, iotype ='out', desc = 'total throat area of venturis', units = 'inch**2')
        
    def __init__(self, *args, **kwargs):
        
        # ----------------------------------------------
        # --- Constructor for the Geometry Component ---
        # ----------------------------------------------
        super(GeometryComp, self).__init__(*args, **kwargs)
        
        # --------------------------------------
        # --- External Code Public Variables ---
        # --------------------------------------
        self.stdout = Path('OpenMDAO') + r"\inj_loc.txt"
        self.stderr = Path('OpenMDAO') + r"\stderr.log"
        
        self.external_files = [FileMetadata(path = self.stdout), FileMetadata(path = self.stderr)]  
        
        try: 
            os.remove(Path('STEP') + '\Config' + self.config + '.step')
        except OSError:
            pass        
        
        self.force_execute = True

        venturi_throat_dia = 0.55 # --- inch
        sector_angle = 12.0 # --- deg
        inner_dia = 15.5 # --- inch
        outer_dia = 22.0 # --- inch
        
        self.sector_area = (pi*(outer_dia/2)**2 - pi*(inner_dia/2)**2) * sector_angle/360.0
        self.venturi_exit_area = pi*(0.85/2)**2 * 4 + pi*(1.15/2)**2
        self.venturi_throat_area = pi*(venturi_throat_dia/2)**2 * 5
        self.injector_dia = 0.85
     
    def execute(self):
        
        # -----------------------------------
        # --- Diagnostic Print Statements --- 
        # ----------------------------------- 
        print '--------------------------------------------------'        
        print 'Generating Geometry for Configuration ' + self.config + ' ...'
        print '{0:20} = {1:10.3f}'.format('Venturi Angle', self.venturi_angle)
        print '{0:20} = {1:10.3f}'.format('Vane Height', self.vane_height)
        print '{0:20} = {1:10.3f}'.format('Pilot Recession Depth', self.pilot_recession)         
        print '--------------------------------------------------' 
        
        # --------------------------------------
        # --- Execute file-wrapped component --- 
        # --------------------------------------

        self.command = [Path('VB') + '\ExportSTEP.exe', 'Config' + self.config, str(self.pilot_recession), str(self.vane_height), str(self.venturi_angle)] 

        # --------------------------------------
        # --- Execute the Geometry Component ---
        # --------------------------------------
        super(GeometryComp, self).execute()
        
        print 'Geometry Generation Complete for Configuration ', self.config  
        
if __name__ == "__main__":
    
    # -------------------------
    # --- Default Test Case ---
    # ------------------------- 
    Geom_Comp = GeometryComp()
    
    Geom_Comp.venturi_angle = 40.0
    Geom_Comp.vane_height = 0.315 #60.0 degrees
    Geom_Comp.pilot_recession = 0.175
    
    Geom_Comp.run()