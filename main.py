import math
import matplotlib.pyplot as plt
import matplotlib.patches as ptch
from utils import *
from algorithms import *
import time
from datetime import datetime

def main(algorithm, save_res, filename):

    #Define the number of available UAVs
    number_uavs = 5
    
    #Define the energy available to each UAV
    energy_uavs = [2000, 2000, 2000, 2000, 2000]

    #Define the decay factor used in the objective function (not used in UninfAtt and Att)
    decay_factor = 0.01
    
    #sensor data from a DJI phantom 4 pro
    hfov = 84 #field of view in degrees
    h_sensor = 50 #meters
    p_o = 50 #percentage of overlap
    pod = 0.63
    
    #grid size calculation
    d = 2*(1-p_o/100)*h_sensor*math.tan(hfov/2) #as per eq 3.1 in PAer report
    
    #cm/pixel
    #cm_p_pixel = 2*h_sensor*math.tan(hfov/2)/num_pixels*100
    
    #mission = Mission(1)
    mission = Mission(2,seed=34)
    
    start_grid = time.time()
    
    optimize_grid = False
    if not optimize_grid:
        cell_centers, sx, sy, theta = mission.generate_grid()
    else:
        cell_centers, sx, sy, theta = mission.generate_optimal_grid_bayesian_opt(100)
    
    end_grid = time.time()
    
    print("time to get grid = ", end_grid-start_grid)
    
    CG_grid, m = define_initial_grid(mission.AOI_vertices, mission.d, True)
    x_lims = (CG_grid[0] - m*mission.d/2, CG_grid[0] + m*mission.d/2)
    y_lims = (CG_grid[1] - m*mission.d/2, CG_grid[1] + m*mission.d/2)
        
    start = time.time()
    start_process = time.process_time()
    

    if algorithm == 'UninfAtt':
        uninf_cell_centers = copy.deepcopy(cell_centers)
        for i in range(len(uninf_cell_centers)):
            if uninf_cell_centers[i].status == 1:
                uninf_cell_centers[i].poc = 1
        Att = Attraction_algorithm(uninf_cell_centers,d, energy_uavs, number_uavs, pod)
        Att.generate_uav_path()
        paths = Att.paths
    elif algorithm == 'Att':
        Att = Attraction_algorithm(cell_centers,d, energy_uavs, number_uavs, pod)
        Att.generate_uav_path()
        paths = Att.paths
    elif algorithm == 'SA':
        type_init = 0
        sync = 0
        n_threads = 15
        SA = SimulatedAnnealing(cell_centers, d, type_init, n_threads, sync, energy_uavs, number_uavs, pod)
        SA.factor = decay_factor
        SA.generate_uav_path(save_res, filename)
        paths = SA.paths
    elif algorithm == 'ACO':
        ACO = AntColonyOpt(cell_centers, d, energy_uavs, number_uavs, pod)
        ACO.factor = decay_factor
        paths = ACO.generate_uav_path(save_res, filename)
    elif algorithm == 'MCTS':
        max_runtime = 200
        C_exploration = 1e-12
        MCTS = MonteCarloTreeSearch(cell_centers, d, energy_uavs, number_uavs, pod)
        MCTS.factor = decay_factor
        MCTS.c = C_exploration
        paths = MCTS.generate_uav_path(max_runtime, save_res, filename)
    
    end = time.time()
    end_process = time.process_time()
    #Do the plots
    fig,ax = plt.subplots()

    """
    #Define the AOI and NFZ as polygons
    AOI_patch = ptch.Polygon(AOI.vertices, alpha=0.2)#, label = 'AOI')
    NFZ_patch = ptch.Polygon(NFZ.vertices, alpha=0.2, facecolor = 'r')#, label = 'NFZ')
    ax.add_patch(NFZ_patch)
    ax.add_patch(AOI_patch)
    """
    #"""
    AOI_patch = ptch.Polygon(mission.AOI_vertices, alpha=0.2)#, label = 'AOI')
    for vert in mission.NFZ_vertices:
        NFZ_patch = ptch.Polygon(vert, alpha=0.2, facecolor = 'r')
        ax.add_patch(NFZ_patch)
    ax.add_patch(AOI_patch)
    blue_patch = ptch.Patch(alpha=0.2, label='AOI')
    red_patch = ptch.Patch(color='r', alpha=0.2, label='NFZ')
    #"""
    #Overlay grid centers in plot
    plot_grid(cell_centers, d, ax, sx, sy, theta)
    plt.xlabel("x (m)")
    plt.ylabel("y (m)")
    plt.xlim(x_lims)
    plt.ylim(y_lims)
    fig.set_size_inches(4.5,4.5)#(4.5,4.5)#(6, 6)
    #plt.title("Example of a mission environment")
    #fig.set_size_inches(10, 8)
    
    # print textstr
    #textstr = 'altitude = %.1f m\nresolution = %.2f cm/pixel\noverlap = %.0f %%\ngrid size = %.2f m\n'%(h_sensor, cm_p_pixel, p_o, d)
    #plt.text(-2.2*CG_grid[0], 2000, textstr, fontsize=12)
    #plt.subplots_adjust(left=0.32)
    
    #path = initial_path
    J, D, ADS = evaluate_paths_multi(paths, cell_centers, pod, save_res, filename)
    #J, D, ADS = evaluate_path(paths[0], cell_centers, save_fig, fig_name)
    #ax2 = ax.twinx()
    #plot the path
    colors = ['orange', 'blue', 'deeppink']
    cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']
    cycle = cycle + cycle + cycle
    cycle[0] = 'blue'
    for uav in range(number_uavs):
        x_coords = [cell.x for cell in paths[uav]]
        y_coords = [cell.y for cell in paths[uav]]
        #plt.scatter(x_coords[0], y_coords[0], marker='o', color=cycle[uav], label = 'start of path of UAV '+str(uav+1))
        plt.plot(x_coords[0], y_coords[0], marker='o', color=cycle[uav])
        plt.plot(x_coords, y_coords, '-', color=cycle[uav], label = 'path of UAV '+str(uav+1))
    #ax2.set_xlim(x_lims)
    #ax2.set_ylim(y_lims)
    
    plt.title( algorithm +" algorithm result, Comp. time = "+"{:.2f}".format(end-start)+"s\nJ = " + "{:.4f}".format(J) + r", $\mathrm{POS}_\mathrm{c}$ = "+"{:.2f}".format(D*100)+"%, EDT = "+"{:.2f}".format(ADS))
    
    if save_res == 1:
        plt.savefig('./results/' + filename + '.pdf', format = 'pdf', bbox_inches='tight')
        f = open('./results/' + filename+'.txt', "a")
        f.write("Bayesian opt iter for grid = 1000\n")
        f.write("time to get grid (s) = " + str(end_grid-start_grid)+"\n")
        f.write("Wall clock time (s)" + str(end-start) + "\n")
        f.write("User+Sys time main process (s)" + str(end_process-start_process) + "\n")
        #f.write("Paths:" + str(SA.paths_idx) + "\n")
        """
        f.write("Wall clock time child processes (s)" + str(list(SA.threads_time[0:SA.number_threads])) + "\n")
        f.write("User+Syst time child processes (s)" + str(list(SA.threads_process_time)) + "\n")
        f.write("Total User+Sys processing time (s)" + str(end_process-start_process + sum(list(SA.threads_process_time[0:SA.number_threads]))) + "\n")
        f.write("Real processing time (s)" + str(end_process-start_process + max(list(SA.threads_process_time[0:SA.number_threads]))) + "\n")
        """
        f.close()
    #for cell in path:
    #    if cell.not_covered == 0:
    #        plt.plot(cell.x, cell.y, 'o', color='k')
    
    print("Wall clock time (s)", end-start)
    print("User+Sys time main process (s)", end_process-start_process)
    """
    print("Wall clock time child processes (s)", list(SA.threads_time[0:SA.number_threads]))
    print("User+Syst time child processes (s)", list(SA.threads_process_time[0:SA.number_threads]))
    print("Total User+Sys processing time (s)", end_process-start_process + sum(list(SA.threads_process_time[0:SA.number_threads])))
    print("Real processing time (s)", end_process-start_process + max(list(SA.threads_process_time[0:SA.number_threads])))
    """
    #for cell in path:
    #    print(cell.poc)
    #plt.show()
    #exit()
    fig2,ax2 = plt.subplots()
    #ax2.add_patch(AOI_patch2)
    #ax2.add_patch(NFZ_patch2)
    #Overlay grid centers in plot
    #plot_grid(cell_centers, d, ax2)
    fig2.set_size_inches(5, 4.3)#(3.5, 3)#(8.17, 7)#(5, 4.3) #(7,6)
    plt.xlabel("x (m)")
    plt.ylabel("y (m)")
    #plt.title("Grid opt. w/ BO 100 iter:\n"r"$s_x = -15.4$m, $s_y = 33.4$m, $\theta = -6.52$deg"+"\n"+r"$\text{POC}_\text{total}$ = 67.32%, Grid comp. time = 580.35s")
    #plt.title(r"No grid optimization: $s_x = 0$, $s_y = 0$, $\theta = 0$"+"\n"+r"$\text{POC}_\text{total}$ = 66.14%, Grid comp. time = 6.52s")

    plt.title("Total POC = "+"{:.2f}".format(100*sum(cell.poc for cell in cell_centers if cell.status == 1))+" % Grid comp. time = "+"{:.2f}".format(end_grid-start_grid)+"s\n"+r"$s_x = 0$, $s_y = 0$, $\theta = 0$")
    #plt.title("Example of a POC map")
    plt.xlim(x_lims)
    plt.ylim(y_lims)
    plot_poc(cell_centers, d, ax2, sx, sy, theta)
    #plt.savefig('./other_figures/34_poc_opt_100.pdf', format = "pdf", bbox_inches='tight')
    #size = math.floor(math.sqrt(len(aux_poc)))
    #aux_poc = np.array(aux_poc)
    #plt.imshow(aux_poc.reshape(size,size), interpolation='bilinear')
    plt.show()
    
    
    

if __name__=='__main__':
    if len(sys.argv) < 3 or len(sys.argv) > 5:
        print("The required structure of command line arguments is the following:")
        print("     python main.py <host> <algorithm> [save result] [file name]")
        print("The following options are available:")
        print("     host (required argument): local or remote")
        print("     algorithm (required argument): UninfAtt, Att, SA, ACO or MCTS")
        print("     save result (optional argument, default 0): 0 or 1")
        print("     file name (optional argument, default is date-time): whatever_name_you_want")
        print("Example command:")
        print("     python main.py local Att 1 john_file_doe")
        exit()
    
    if sys.argv[1] not in ['local', 'remote']:
        print("The required structure of command line arguments is the following:")
        print("     python main.py <host> <algorithm> [save result] [file name]")
        print("The host argument must be: local or remote")
        exit()
    
    if sys.argv[2] not in ['UninfAtt', 'Att', 'SA', 'ACO', 'MCTS']:
        print("The required structure of command line arguments is the following:")
        print("     python main.py <host> <algorithm> [save result] [file name]")
        print("An invalid algorithm was chosed. The available options are: UninfAtt, Att, SA, ACO or MCTS.")
        exit()
    
    if len(sys.argv) > 3 and sys.argv[3] not in ['0', '1']:
        print("The required structure of command line arguments is the following:")
        print("     python main.py <host> <algorithm> [save result] [file name]")
        print("The host argument must be: 0 or 1")
        exit()
    
    #Check if the results should be saved
    if len(sys.argv) > 3:
        save_res = int(sys.argv[3])
    else:
        save_res = 0
        

    if sys.argv[1] == 'local':
        if len(sys.argv) == 5:
            filename = sys.argv[4]
        else:
            now = datetime.now()
            filename = now.strftime("%d_%m_%Y_%Hh_%Mm_%Ss")
        main(sys.argv[2], save_res, filename)
    
    if sys.argv[1] == 'remote':
        #Used for graph search
        import paramiko
        #from remote_plot import plt
        print("Running algorithm in server")
        # Connect to remote host
        client = paramiko.SSHClient()
        client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        client.connect(hostname = '', username='', password='')
        
        # Setup sftp connection and transmit this script
        sftp = client.open_sftp()
        try:
            sftp.chdir('./UAV-path-planning')  # Test if remote_path exists
        except IOError:
            sftp.mkdir('./UAV-path-planning')  # Create remote_path
            sftp.chdir('./UAV-path-planning')

        try:
            sftp.chdir('./results')  # Test if remote_path exists
            sftp.chdir('..')
        except IOError:
            sftp.mkdir('./results')  # Create remote_path
        
        #sftp.put(__file__, '/tmp/code/main.py')
        sftp.put(__file__, './main.py')
        sftp.put('./algorithms.py', './algorithms.py')
        sftp.put('./utils.py', './utils.py')
        sftp.put('./search.py', './search.py')
        sftp.put('./utils2.py', './utils2.py')
        sftp.close()
        
        if len(sys.argv) == 5:
            filename = sys.argv[4]
        else:
            now = datetime.now()
            filename = now.strftime("%d_%m_%Y_%Hh_%Mm_%Ss")
        
        stdin, stdout, stderr = client.exec_command('cd ./UAV-path-planning; python3 ./main.py local ' + sys.argv[2] +' '+ str(save_res) +' '+ filename)
        
        stdin.close()
        for line in stdout:
            # Process each line in the remote output
            print(line, end="")
        
        for line in stderr:
            print(line)
        
        if save_res == 1:
            sftp = client.open_sftp()
            sftp.get("./UAV-path-planning/results/" + filename + '.txt', "./results/"+filename + '.txt')
            sftp.get("./UAV-path-planning/results/" + filename + '.pdf', "./results/"+filename + '.pdf')
            sftp.close()
        
        client.close()
        sys.exit(0)
    
    
    
