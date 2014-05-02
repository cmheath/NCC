''' OpenMDAO tecplot component '''

# --- Inherent python/system level imports
import os
import shutil
from glob import iglob

# --- OpenMDAO imports
from openmdao.lib.datatypes.api import Str, Float, Bool
from openmdao.main.api import Component
from pathname import Path
    
def copy_files(src_glob, dst_folder):
    for fname in iglob(src_glob):
        shutil.copy2(fname, dst_folder)

def remove_files(filepath):
    for fname in iglob(filepath):
        os.remove(fname)        
    
class TecplotComp(Component):
    ''' OpenMDAO component for Processing Flow Solution Data in Tecplot '''

    config = Str('1_5', iotype ='in', desc ='current design configuration number')   
    nonreacting = Bool(False, iotype = 'in', desc = 'flag specifier to run ncc nonreacting')
    reacting = Bool(False, iotype = 'in', desc = 'flag specifier to run ncc with liquid spray')    
    reacting2 = Bool(False, iotype = 'in', desc = 'flag specifier to run ncc with liquid spray')
    dPqPBase = Float(0.00, iotype ='out', desc = 'Static pressure drop through the combustor')     
    dTqTBase = Float(0.00, iotype ='out', desc = 'Static temperature drop at combustor exit between successive iterations') 
        
    def __init__(self, *args, **kwargs):
        
        # ---------------------------------------------
        # --- Constructor for the Tecplot Component ---
        # ---------------------------------------------
        super(TecplotComp, self).__init__(*args, **kwargs)
         
        self._config = '' 
        self.T4_prev = 0.0

    def execute(self):   
    
        if self.nonreacting:
            self._config = self.config + '_NR'   
        else:
            self._config = self.config + '_R'
                           
        print 'Copying Output to ' + os.path.join(Path('Tecplot'), 'Sim-' + self._config)
        
        try: 
            shutil.rmtree(os.path.join(Path('Tecplot'), 'Sim-' + self._config), True)
            os.remove('patran.out')
        except OSError:
            pass        
        
        os.mkdir(os.path.join(Path('Tecplot'), 'Sim-' + self._config))
        
        indir = os.path.join(Path('Tecplot'), 'Sim-' + self._config, 'Input')
        outdir = os.path.join(Path('Tecplot'), 'Sim-' + self._config, 'Output')
        imdir = os.path.join(Path('Tecplot'), 'Sim-' + self._config, 'Images')
       
        os.mkdir(indir) 
        os.mkdir(outdir)
        os.mkdir(imdir)
             
        try: 
            copy_files(os.path.join(Path('OpenMDAO'), 'tecplot.dat'), outdir)
            copy_files(os.path.join(Path('OpenMDAO'), 'open_ncc_resd.out'), outdir)
            copy_files(os.path.join(Path('OpenMDAO'), 'tecplot_unsteady_droplets.dat'), outdir)
            copy_files(os.path.join(Path('OpenMDAO'), '*.in*'), indir)  
            copy_files(os.path.join(Path('OpenMDAO'), '*.out*'), outdir)
            copy_files(os.path.join(Path('OpenMDAO'), '*.png*'), imdir)
            copy_files(os.path.join(Path('OpenMDAO'), 'Ts.txt'), outdir)
            copy_files(os.path.join(Path('OpenMDAO'), 'NO.txt'), outdir)
            copy_files(os.path.join(Path('OpenMDAO'), 'TKE.txt'), outdir)
            copy_files(os.path.join(Path('OpenMDAO'), 'Unmixedness.txt'), outdir)            
        except:
            pass     
            
        print '--------------------------------------------------'      
        print 'Reading Post-Processed Solution for Config ' + self._config + ' ...'            
        print '--------------------------------------------------'
        
        # --------------------------------------
        # --- Read Tecplot Component Results ---
        # --------------------------------------
        with open(os.path.join(outdir, 'post-1d.out'), 'r') as f:

            line = f.readline()
            while 'BCflag =    1' not in line:
                line = f.readline()
            while 'Ptot_avg {Pa}' not in line:
                line = f.readline()        
            P1 = float(line.split()[-1])
            
            f.seek(0)
            line = f.readline()
            while 'BCflag =    2' not in line:
                line = f.readline()
            while 'Ptot_avg {Pa}' not in line:
                line = f.readline()        
            P2 = float(line.split()[-1])
            
            if self.reacting or self.reacting2:
                f.seek(0)
                line = f.readline()
                while 'BCflag =    2' not in line:
                    line = f.readline()
                while 'T_avg {K}' not in line:
                    line = f.readline()        
                T4_avg = float(line.split()[-1])   
                print 'Average Exit Ts: ', T4_avg
                self.dTqTBase = abs(T4_avg - self.T4_prev) / T4_avg
                self.T4_prev = T4_avg
            
            self.dPqPBase = (P1 - P2) / P1
        
        print 'Config %s : dPqPBase = %.10f, P1 = %.10f, P2 = %.10f' %(self._config, self.dPqPBase, P1, P2)        
        print 'Config %s : dTqTBase = %.10f' %(self._config, self.dTqTBase)               
              
        print 'Processing Complete for Configuration ' + self._config
        print '--------------------------------------------------'        
        print 'Cleaning up Current Working Directory ...'  
        
        try:
            remove_files(os.path.join(os.getcwd(), '*.png'))
            remove_files(os.path.join(os.getcwd(), '*.out.*'))
            remove_files(os.path.join(os.getcwd(), 'Ts.txt'))
            remove_files(os.path.join(os.getcwd(), 'NO.txt'))
            remove_files(os.path.join(os.getcwd(), 'TKE.txt'))
            remove_files(os.path.join(os.getcwd(), 'Unmixedness.txt'))        
        except:
            pass
        
    
if __name__ == "__main__":
    
    # -------------------------
    # --- Default Test Case ---
    # ------------------------- 
    Tecplot_Comp = TecplotComp()
    Tecplot_Comp.reacting2 = True
    Tecplot_Comp.run()