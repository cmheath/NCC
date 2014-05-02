import os
def Path(key):
    
    pathname = {}
    pathname['SolidWorks_Multi'] = 'Z:\NASA\RTM\Research\LDI2\SolidWorks'  
    pathname['OpenMDAO'] = os.getcwd()
    pathname['NPSS'] = os.path.join(os.getcwd(), 'NPSS') 
    pathname['Patran'] = 'Z:\NASA\RTM\Research\LDI2\Patran'
    pathname['VB'] = 'Z:\NASA\RTM\Research\LDI2\VB'
    pathname['STEP'] = 'Z:\NASA\RTM\Research\LDI2\STEP'
    pathname['Cubit'] = 'Z:\NASA\RTM\Research\LDI2\Cubit'
    pathname['Tecplot'] = 'Z:\NASA\RTM\Research\LDI2\Tecplot'
    pathname['Scripts'] = 'Z:\NASA\RTM\Research\LDI2\Scripts'    
    return pathname[key]

    
