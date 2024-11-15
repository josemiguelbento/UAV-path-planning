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

    #sensor data from a DJI phantom 4 pro
    hfov = 84 #field of view in degrees
    h_sensor = 50 #meters
    p_o = 50 #percentage of overlap
    pod = 0.63
    
    #grid size calculation
    d = 2*(1-p_o/100)*h_sensor*math.tan(hfov/2) #as per eq 3.1 in PAer report
    
    #cm/pixel
    #cm_p_pixel = 2*h_sensor*math.tan(hfov/2)/num_pixels*100
    
    #mission = Mission('simple_square')
    #mission = Mission(1)
    mission = Mission(2,seed=34)
    
    start_grid = time.time()
    #cell_centers = mission.generate_grid()
    cell_centers, sx, sy, theta = mission.generate_optimal_grid()
    #cell_centers, sx, sy, theta = mission.generate_optimal_grid_bayesian_opt(100)
    end_grid = time.time()
    
    print("time to get grid = ", end_grid-start_grid)
    
    CG_grid, m = define_initial_grid(mission.AOI_vertices, mission.d, True)
    x_lims = (CG_grid[0] - m*mission.d/2, CG_grid[0] + m*mission.d/2)
    y_lims = (CG_grid[1] - m*mission.d/2, CG_grid[1] + m*mission.d/2)
        
    start = time.time()
    start_process = time.process_time()
    

    if algorithm == 'UninfAtt':
        pass
    elif algorithm == 'Att':
        att = Attraction_algorithm(cell_centers,d, energy_uavs, number_uavs, pod)
        att.generate_uav_path()
        paths = att.paths
    elif algorithm == 'SA':
        pass
    elif algorithm == 'ACO':
        pass
    elif algorithm == 'MCTS':
        pass

    
    """
    #initial_path = []
    #cell_centers = calculate_poc(cell_centers, 1500, 1500, 1000, 1000, True)
    SA = SimulatedAnnealing_v2(cell_centers, d, 0, 15, 0, energy_uavs, number_uavs, initial_positions, pod, False) #False)
    #paths = SA.initial_guess_list[0]
    SA.L_k = 2653
    SA.initial_temp = 1.83e-4 #0.0004#0.0085#0.002526 #5e-5
    SA.factor = 0.01
    SA.cooling_factor = 0.954#0.99
    SA.final_temp = 2.11e-5 #2.755e-06 #2.755e-7
    random.seed(time.time())
    SA.generate_uav_path(save_fig, fig_name)
    paths = SA.paths
    #print(SA.paths_idx)
    """

    """
    MCS = MonteCarloSearch(cell_centers, d, energy_uavs, number_uavs, pod)
    MCS.c = 1e-12 #1.399e-5
    MCS.generate_path(200, save_fig, fig_name)
    paths = MCS.paths
    """
    
    """
    aco = AntColonyOpt_v2(cell_centers, d, energy_uavs, number_uavs, pod)
    #aco.max_iterations = 200
    #aco.alpha = 2.4826382170953
    #aco.beta = 0.4309240189521
    #aco.evaporation = 0.0423631296383
    #aco.Q = 0.0267611642549
    #aco.min_heuristic = 8.6134e-07
    aco.max_iterations = 10000
    aco.alpha = 2.73
    aco.beta = 19.96
    aco.evaporation =  0.690
    aco.Q = 9.56
    aco.min_heuristic = 9.31e-7

    paths = aco.run_single_ACO(save_fig, fig_name)
    """
    """
    valid_cell_list = [cell for cell in cell_centers if cell.status == 1]
    paths_idx = [[810, 811, 837, 836, 835, 859, 877, 878, 879, 861, 837, 862, 838, 812, 813, 789, 763, 734, 733, 762, 761, 760, 786, 787, 788, 812, 811, 810, 836, 860, 861, 837, 838, 862, 863, 864, 840, 814, 839, 863, 839, 813, 814, 790, 814, 840, 864, 863, 862, 861, 879, 878, 877, 859, 860, 859, 858, 834, 808, 784, 808, 809, 810, 811, 787, 786, 785, 809, 835, 859, 878, 860, 836, 835, 809, 785, 810, 836, 837, 838, 812, 788, 762, 734, 735, 736, 764, 790, 789, 788, 787, 786, 785, 759, 760, 761, 762, 761, 762, 733, 761, 787, 811, 838, 837, 811, 787, 788, 789, 790, 789, 812, 838, 837, 861, 862, 863, 839, 813, 812, 811, 810, 835, 836, 810, 811, 812, 838, 861, 862, 838, 812, 811, 837, 861, 879], [571, 604, 605, 604, 603, 602, 603, 570, 537, 570, 569, 568, 567, 600, 626, 627, 601, 568, 601, 628, 656, 683, 684, 657, 629, 657, 656, 655, 654, 626, 627, 628, 629, 630, 658, 631, 605, 572, 539, 538, 537, 536, 535, 536, 504, 472, 504, 505, 506, 474, 473, 472, 471, 470, 502, 501, 533, 566, 567, 600, 626, 653, 625, 599, 567, 534, 535, 536, 569, 602, 628, 629, 630, 631, 605, 572, 539, 538, 571, 604, 603, 602, 569, 568, 601, 627, 601, 568, 535, 503, 504, 505, 506, 538, 570, 603, 570, 536, 504, 505, 537, 570, 603, 629, 630, 631, 605, 572, 571, 570, 537, 505, 536, 569, 602, 601, 602, 603, 570, 571, 572, 605, 604, 603, 602, 601, 568, 569, 570, 569, 602, 628, 656, 627, 601, 568, 535, 503], [457, 458, 459, 460, 459, 458, 457, 489, 490, 491, 523, 556, 523, 491, 490, 489, 521, 522, 523, 491, 459, 430, 429, 458, 491, 492, 493, 492, 491, 490, 522, 555, 556, 524, 492, 460, 461, 460, 459, 430, 399, 428, 427, 456, 488, 520, 521, 522, 490, 489, 521, 554, 587, 554, 553, 520, 488, 456, 427, 428, 429, 430, 400, 430, 429, 428, 457, 489, 521, 554, 555, 556, 555, 522, 490, 491, 523, 522, 521, 489, 457, 458, 459, 430, 400, 368, 399, 398, 427, 456, 488, 520, 521, 522, 523, 524, 492, 460, 461, 460, 459, 491, 490, 457, 428, 457, 489, 488, 456, 427, 428, 457, 489, 457, 428, 398, 399, 429, 458, 490, 458, 429, 458, 490, 458, 429, 399, 429, 458, 490, 489, 488, 489, 490, 491, 459, 458, 490], [750, 722, 692, 693, 723, 692, 722, 721, 720, 690, 663, 664, 665, 638, 612, 638, 637, 610, 636, 663, 690, 691, 664, 665, 666, 693, 666, 667, 694, 724, 752, 751, 723, 724, 725, 695, 668, 694, 723, 722, 721, 720, 719, 747, 748, 749, 750, 722, 721, 720, 748, 774, 773, 747, 719, 747, 774, 775, 776, 775, 774, 773, 747, 748, 749, 750, 751, 752, 753, 752, 751, 776, 777, 778, 752, 724, 694, 667, 640, 639, 638, 611, 610, 636, 663, 690, 720, 690, 691, 664, 665, 692, 722, 750, 776, 749, 721, 691, 690, 720, 748, 749, 721, 691, 692, 693, 723, 751, 750, 722, 692, 691, 690, 663, 664, 665, 692, 722, 750, 749, 721, 691, 664, 637, 664, 691, 721, 691, 664, 637, 611, 612, 638, 665, 691, 721, 749, 775], [242, 274, 273, 272, 240, 241, 209, 208, 207, 238, 239, 208, 172, 138, 105, 104, 137, 171, 207, 239, 271, 272, 273, 274, 275, 276, 309, 308, 307, 337, 307, 275, 243, 211, 176, 212, 176, 175, 174, 173, 139, 140, 174, 210, 242, 241, 240, 208, 172, 171, 170, 206, 207, 208, 209, 210, 242, 241, 240, 239, 271, 304, 305, 306, 307, 308, 275, 243, 211, 175, 174, 140, 174, 210, 242, 241, 209, 173, 172, 171, 172, 173, 210, 242, 274, 275, 243, 211, 175, 141, 140, 139, 138, 139, 106, 107, 79, 108, 141, 107, 106, 105, 138, 172, 171, 137, 138, 139, 173, 208, 240, 273, 306, 305, 272, 271, 270, 238, 206, 207, 208, 209, 210, 209, 208, 207, 206, 238, 239, 240, 241, 209, 208, 240, 272, 273, 241, 209]]
    #paths_idx = [[503, 502, 533, 534, 567, 568, 535, 536, 569, 570, 571, 537, 504, 471, 437, 438, 472, 439, 473, 440, 474, 505, 506, 538, 539, 507, 475, 476, 508, 540, 572, 605, 604, 630, 631, 658, 657, 629, 603, 602, 601, 600, 599, 626, 627, 628, 656, 683, 684, 685, 711, 710, 682, 709, 740, 768, 739, 708, 681, 655, 654, 625, 653, 680, 652, 679, 706, 705, 678, 677, 676, 649, 621, 596, 597, 623, 651, 624, 598, 565, 566, 600, 627, 656, 684, 712, 742, 741, 769, 795, 770, 796, 820, 819, 845, 818, 794, 817, 793, 792, 816, 842, 843, 844, 843, 842, 841, 840, 839, 813, 790, 764, 789, 788, 787, 786, 759, 785, 810, 837, 838, 862, 861, 860, 836, 811, 812], [751, 777, 778, 802, 801, 823, 822, 821, 798, 774, 799, 800, 776, 775, 749, 721, 722, 750, 723, 693, 667, 666, 692, 691, 720, 748, 773, 797, 772, 771, 745, 744, 716, 717, 746, 718, 747, 719, 690, 663, 636, 610, 637, 611, 582, 581, 580, 548, 549, 550, 518, 519, 551, 583, 612, 638, 664, 665, 639, 640, 641, 668, 694, 695, 724, 725, 752, 753, 779, 780, 804, 826, 803, 825, 824, 847, 848, 866, 867, 849, 850, 868, 869, 870, 871, 853, 830, 831, 855, 832, 807, 783, 757, 729, 730, 731, 732, 733, 703, 677, 650, 622, 595, 563, 530, 497, 529, 562, 561, 594, 593, 592, 591, 559, 560, 528, 496, 495, 527, 494, 462, 463, 464, 465, 431, 432, 403], [879, 861, 862, 837, 810, 811, 786, 787, 788, 813, 839, 812, 838, 863, 864, 840, 865, 841, 815, 791, 814, 790, 789, 762, 763, 764, 765, 766, 767, 738, 707, 737, 736, 704, 735, 734, 733, 761, 760, 759, 758, 785, 784, 808, 809, 835, 836, 860, 859, 878, 877, 876, 858, 857, 834, 833, 856, 875, 887, 886, 874, 873, 885, 884, 872, 854, 829, 852, 851, 828, 806, 827, 805, 782, 756, 755, 781, 754, 726, 727, 728, 698, 671, 699, 672, 645, 644, 643, 670, 697, 696, 669, 642, 613, 614, 615, 584, 616, 617, 585, 552, 586, 587, 618, 646, 673, 700, 701, 702, 674, 675, 648, 647, 619, 620, 590, 557, 524, 491, 490, 489, 457, 428, 458, 459, 460], [211, 174, 210, 209, 173, 139, 138, 172, 208, 207, 171, 137, 105, 104, 136, 169, 170, 205, 236, 269, 237, 206, 238, 270, 239, 240, 273, 305, 336, 306, 274, 307, 308, 339, 372, 373, 374, 340, 310, 276, 277, 244, 212, 175, 176, 143, 177, 178, 145, 112, 111, 144, 110, 82, 81, 80, 107, 108, 109, 142, 141, 140, 106, 78, 79, 55, 54, 77, 76, 53, 34, 33, 51, 74, 101, 102, 134, 133, 167, 166, 203, 235, 204, 168, 135, 103, 75, 52, 51, 50, 73, 72, 49, 32, 31, 30, 48, 71, 100, 99, 98, 131, 165, 201, 202, 234, 266, 233, 265, 298, 297, 296, 331, 332, 333, 334, 299, 267, 268, 237, 206, 207, 208, 209, 241], [491, 459, 458, 490, 489, 522, 555, 521, 488, 520, 553, 554, 588, 589, 557, 558, 526, 493, 525, 524, 556, 523, 492, 461, 460, 430, 429, 400, 399, 398, 397, 427, 428, 457, 456, 455, 426, 396, 365, 366, 367, 368, 335, 300, 301, 302, 303, 304, 271, 272, 241, 242, 243, 275, 309, 338, 371, 337, 370, 369, 401, 402, 403, 433, 467, 466, 498, 499, 531, 564, 532, 500, 501, 468, 434, 404, 405, 435, 469, 470, 436, 407, 406, 375, 408, 376, 341, 342, 311, 343, 312, 279, 278, 245, 246, 280, 313, 314, 344, 345, 346, 381, 380, 412, 413, 443, 442, 441, 411, 379, 378, 377, 409, 410, 441, 475, 507, 506, 505, 537, 570, 603, 602, 601, 567, 568, 569, 570, 571]]
    paths = [[valid_cell_list[cell_idx] for cell_idx in path] for path in paths_idx]

    #cell_centers, sx, sy, theta = mission.generate_optimal_grid()
    end = start + 496.17#516.3015468120575
    
    """
    
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
    #plt.title("Attraction algorithm for 4000kJ available energy\n" + "J = " + "{:.4f}".format(J) + " D = "+"{:.2f}".format(D*100)+"% ADS = "+"{:.2f}".format(ADS) + " Computational time = "+"{:.2f}".format(end-start)+"s")
    #plt.title("Uninformed baseline for 4000kJ available energy\n" + "J = " + "{:.4f}".format(J) + " D = "+"{:.2f}".format(D*100)+"% ADS = "+"{:.2f}".format(ADS) + " Computational time = "+"{:.2f}".format(end-start)+"s")
    
    #plt.title("SA algorithm (Att initial guess)\n" + "J = " + "{:.4f}".format(J) + " D = "+"{:.2f}".format(D*100)+"% ADS = "+"{:.2f}".format(ADS) + " Comp. time = "+"{:.2f}".format(end-start)+"s")
    
    #plt.title("ACO algorithm\n" + "J = " + "{:.4f}".format(J) + " D = "+"{:.2f}".format(D*100)+"% ADS = "+"{:.2f}".format(ADS) + " Comp. time = "+"{:.2f}".format(end-start)+"s")
    #plt.title("MCTS algorithm result, Comp. time="+"{:.2f}".format(end-start)+"s\nJ = " + "{:.4f}".format(J) + r", $\mathrm{POS}_\mathrm{c}$ = "+"{:.2f}".format(D*100)+"%, EDT = "+"{:.2f}".format(ADS))
    
    #plt.title("MCTS algorithm result\nJ = " + "{:.4f}".format(J) + r", $\mathrm{POS}_\mathrm{c}$ = "+"{:.2f}".format(D*100)+"%, EDT = "+"{:.2f}".format(ADS)+"\nComp. time = "+"{:.2f}".format(end-start)+"s", fontsize=14)
    #plt.title("SA algorithm result\nJ = " + "{:.4f}".format(J) + r", $\mathrm{POS}_\mathrm{c}$ = "+"{:.2f}".format(D*100)+"%, EDT = "+"{:.2f}".format(ADS)+"\nComp. time = "+"{:.2f}".format(end-start)+"s", fontsize=14)
    #plt.savefig('./other_figures/SA_paths_pod10_small_12_10_2024_10h_54m_02s.pdf', format = "pdf", bbox_inches='tight')
    #plt.savefig('./other_figures/SA_paths_pod90_small_12_10_2024_10h_51m_06s.pdf', format = "pdf", bbox_inches='tight')
    
    plt.title("SA algorithm result, Comp. time = "+"{:.2f}".format(end-start)+"s\nJ = " + "{:.4f}".format(J) + r", $\mathrm{POS}_\mathrm{c}$ = "+"{:.2f}".format(D*100)+"%, EDT = "+"{:.2f}".format(ADS))
    #plt.title("Att algorithm result, Comp. time = "+"{:.2f}".format(end-start)+"s\nJ = " + "{:.4f}".format(J) + r", $\mathrm{POS}_\mathrm{c}$ = "+"{:.2f}".format(D*100)+"%, EDT = "+"{:.2f}".format(ADS))
    #plt.title("Uninf. Att algorithm result, Comp. time = 1.07s\nJ = " + "{:.4f}".format(J) + r", $\mathrm{POS}_\mathrm{c}$ = "+"{:.2f}".format(D*100)+"%, EDT = "+"{:.2f}".format(ADS))
    #plt.title("ACO algorithm result, Comp. time = "+"{:.2f}".format(end-start)+"s\nJ = " + "{:.4f}".format(J) + r", $\mathrm{POS}_\mathrm{c}$ = "+"{:.2f}".format(D*100)+"%, EDT = "+"{:.2f}".format(ADS))
    
    
    #handles, labels = ax.get_legend_handles_labels()
    #handles.reverse()
    #handles.append(blue_patch)
    #handles.append(red_patch)
    #handles.reverse()
    #lgnd = fig.legend(handles=handles, bbox_to_anchor=(0.59, 0.89), framealpha=0.8, ncol=2, columnspacing=0.5)#, bbox_to_anchor=(0.9, 0.1), framealpha=1)
    #lgnd = fig.legend(loc = "upper right", bbox_to_anchor=(0.9, 0.87), framealpha=0.8)
    #lgnd.legend_handles[2]._sizes = [4]
    #lgnd.legend_handles[3]._sizes = [4]
    
    #ax2.legend(bbox_to_anchor=(0.56, 0.72))
    #plt.savefig('./other_figures/example_aoi.pdf', format = "pdf")
    #"""
    
    #plt.savefig('./other_figures/example_mission_env.pdf', format = "pdf", bbox_inches='tight')
    if save_res == 1:
        #print(os.getcwd())
        #plt.savefig('./tese/code_multi_uavs/result_images/' + fig_name + '.pdf', format = 'pdf', bbox_inches='tight')
        plt.savefig('./tese/code_var_poc/result_images/' + filename + '.pdf', format = 'pdf', bbox_inches='tight')
        f = open('./tese/code_var_poc/result_images/' + filename+'.txt', "a")
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
    plt.title("Grid opt. w/ BO 100 iter:\n"r"$s_x = -15.4$m, $s_y = 33.4$m, $\theta = -6.52$deg"+"\n"+r"$\text{POC}_\text{total}$ = 67.32%, Grid comp. time = 580.35s")
    #plt.title(r"No grid optimization: $s_x = 0$, $s_y = 0$, $\theta = 0$"+"\n"+r"$\text{POC}_\text{total}$ = 66.14%, Grid comp. time = 6.52s")

    #plt.title("Total POC = "+"{:.2f}".format(100*sum(cell.poc for cell in att.valid_cell_list))+" % Grid comp. time = 6.52 s")
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
        #sftp.put(__file__, '/tmp/code/main.py')
        sftp.put(__file__, './main.py')
        sftp.put('./algorithms.py', './algorithms.py')
        sftp.put('./utils.py', './utils.py')
        sftp.put('./search.py', './search.py')
        sftp.put('./utils2.py', './utils2.py')
        #sftp.put('./algorithms_v2.py', './tese/code_multi_uavs/algorithms_v2.py')
        sftp.close()
        
        # Run the transmitted script remotely without args and show its output.
        # SSHClient.exec_command() returns the tuple (stdin,stdout,stderr)
        
        if len(sys.argv) == 5:
            filename = sys.argv[4]
        else:
            now = datetime.now()
            filename = now.strftime("%d_%m_%Y_%Hh_%Mm_%Ss")
        
        stdin, stdout, stderr = client.exec_command('python3 ./main.py local ' + sys.argv[2] +' '+ str(save_res) +' '+ filename)
        
        stdin.close()
        for line in stdout:
            # Process each line in the remote output
            print(line, end="")
        
        for line in stderr:
            print(line)
        
        #sftp = client.open_sftp()
        #sftp.get("./tese/code/result_images/" + img_name + '.png', "./result_images/"+img_name + '.png')
        #sftp.get("./tese/code/result_images/" + img_name + '.txt', "./result_images/"+img_name + '.txt')
        #sftp.close()
        
        client.close()
        sys.exit(0)
    
    
    
