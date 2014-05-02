def init_dP():

    # --------------------------------------
    # --- Read Tecplot Component Results ---
    # --------------------------------------
    try:
        with open('post-1d.out', 'r') as f:

            line = f.readline()
            while 'BCflag =    1' not in line:
                line = f.readline()
                if line is '':
                    raise IOError
            while 'Ptot_avg {Pa}' not in line:
                line = f.readline()
                if line is '':
                    raise IOError      
            P1 = float(line.split()[-1])
            
            f.seek(0)
            line = f.readline()
            while 'BCflag =    2' not in line:
                line = f.readline()
                if line is '':
                    raise IOError                
            while 'Ptot_avg {Pa}' not in line:
                line = f.readline()
                if line is '':
                    raise IOError                  
            P2 = float(line.split()[-1])
                    
            dPqPBase = (P1 - P2) / P1   

        return dPqPBase

    except:
        return 0.03 # --- Default from init
        pass
