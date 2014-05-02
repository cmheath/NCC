import subprocess
from sys import argv

if __name__ == "__main__":
    
    # argv1, argv2, argv3
    #'/run_nonreacting', '400', '3:15:00', 'normal', Sim-1'
    
    reacting_state = argv[1]
    num_procs = argv[2]
    walltime = argv[3]
    queue = argv[4]
    config = argv[5]

    case = config.split('_')

    restart_directory = ('%s_%s_%s_' % (case[0], int(case[1])-1, case[2]))

    # --- Execute preprocessing commands
    cmd = ['csh', './preprocess', reacting_state, num_procs, restart_directory]
    subprocess.call(cmd)
    
    # --- Run NCC 
    proc = subprocess.Popen(['pwd'], stdout=subprocess.PIPE, stderr=subprocess.PIPE) 
    working_directory, err = proc.communicate()
    working_directory = working_directory[0:-1]
    
    cmd = ['chmod', '777', working_directory + reacting_state]
    print cmd
    subprocess.call(cmd)

    cmd = ['qsub', '-N', config, '-q', queue, '-l select=' + str(int(num_procs)/20) + ':ncpus=20:mpiprocs=20:model=ivy', '-l walltime=' + walltime, working_directory + reacting_state]
    print cmd
    subprocess.call(cmd)
    
    # --- Execute post-processing commands and run Tecplot
    cmd = ['csh', './postprocess', reacting_state, num_procs, config]
    print cmd
    subprocess.call(cmd)