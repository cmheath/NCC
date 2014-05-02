''' OpenMDAO component wrapper for Cubit mesh component'''
# --- Inherent python/system level imports
import logging 
import shutil
import os

# --- OpenMDAO main and library imports
from openmdao.main.api import Component, FileMetadata, enable_console
from openmdao.lib.datatypes.api import Float, Int, Str
from openmdao.lib.components.api import ExternalCode
from openmdao.main.resource import ResourceAllocationManager as RAM
from nas_access import NAS_Allocator

# --- OpenMDAO component imports
from pathname import Path

class MeshComp(ExternalCode):

    config = Str('1', iotype ='in', desc ='current design configuration number')     
    pilot_recession = Float(0.0, iotype = 'in', desc = 'axial location of pilot dome')
    venturi_angle = Float(30.0, iotype ='in', desc = 'angle of converging/diverging venturi section', units = 'deg')
    vane_height = Float(0.317293, iotype ='in', desc = 'height of each helical vane blade')
 
    
    def __init__(self, *args, **kwargs):
        '''
         Constructor for the NCC component ---
        '''
        super(MeshComp, self).__init__(*args, **kwargs)                
        self.force_execute = True            

    def execute(self):
        '''
         Executes file-wrapped component --- 
        '''               
                
        print 'Updating Local Directory Structure ...'        
        shutil.rmtree(Path('Patran') + '\Config' + self.config, True)
        os.mkdir(Path('Patran') + '\Config' + self.config)
        shutil.copy2(Path('STEP') + '\Config' + self.config + '.step', Path('OpenMDAO'))      
        shutil.copy2(Path('Scripts') + '/run_cubit.py', Path('OpenMDAO'))
        
        print 'Generating Mesh Remotely for Configuration ' + self.config + ' ...' 

        self.command = ['python', 'run_cubit.py', self.config, str(self.pilot_recession), str(self.venturi_angle), str(self.vane_height)]        
  
        # ------------------------------       
        # --- Job Submission Command ---
        # ------------------------------
        print 'Starting Job on NAS ...'             
        print Path('STEP') + '\Config' + self.config + '.step'
        self.external_files = [FileMetadata('run_cubit.py', input=True),       
               FileMetadata('Config' + self.config + '.step', input=True),               
               FileMetadata('patran.out', output=True, binary=True),
               FileMetadata('Config' + self.config + '.cub', output=True, binary=True)]       
                          
        self.resources = {'localhost': False}
            
        print 'Job Complete ...'   

        # -----------------------------
        # --- Execute the component ---
        # -----------------------------
        super(MeshComp, self).execute() 

        # --- Copy files to local directories
        shutil.copy2(Path('OpenMDAO') + '/patran.out', Path('Patran') + '\Config' + self.config)        
        shutil.copy2(Path('OpenMDAO') + '/Config' + self.config + '.cub', Path('Cubit'))

if __name__ == "__main__":
 
    enable_console() 
    logging.getLogger().setLevel(logging.DEBUG) 

    # --- Create 4 allocators (redundancy for reliability) adding each to the manager ---
    allocator1 = NAS_Allocator(name='PFE20_DMZ1',
                              dmz_host='dmzfs1.nas.nasa.gov',
                              server_host='pfe20')                             
                              
    RAM.add_allocator(allocator1)   

    Remote_Mesh_Comp = MeshComp()
    Remote_Mesh_Comp.run()
    