import math
import shapely.geometry as geom
import sys
import matplotlib.pyplot as plt
import matplotlib.patches as ptch
import numpy as np
import pylab as pl
import matplotlib.colorbar as cbar
import search
import copy
import random
import heapq
from typing import List, Tuple
from ConfigSpace import Configuration, ConfigurationSpace, Float
from smac import HyperparameterOptimizationFacade, Scenario
from smac.runhistory.dataclasses import TrialValue, TrialInfo
from datetime import datetime
import time
from scipy import integrate

class Mission():
    """Class that defines a mission instance."""
    def __init__(self, mission_selection = 0, seed = 0, hfov = 84, h_sensor = 50, p_o = 50):
        if mission_selection == 0: #Pre-defined mission with a small AOI and a POC map with a single gaussian distribution.
            self.AOI_vertices = [[50, 50], [1525, 500], [2400, 1250], [1200, 2350],  [50, 1500]]
            self.AOI_polygon = geom.polygon.Polygon(self.AOI_vertices)
            self.NFZ_vertices = [[[550,550], [550, 950], [950, 950], [950, 550]]]
            self.NFZ_polygons = [geom.polygon.Polygon(vert) for vert in self.NFZ_vertices]
            self.POC_weights = [1]
            self.POC_means = [[1500,1500]]
            self.POC_variance_matrices = [[[1000**2, 0],[0, 1000**2]]]
            
        elif mission_selection == 1: #Pre-defined mission with a medium-sized AOI and a POC map with a weighted sum of five different gaussian distributions.
            self.AOI_vertices = [[50, 1000], [1500, 50], [4800, 1000], [4000, 3500], [2000, 4000], [50, 3000]]
            self.AOI_polygon = geom.polygon.Polygon(self.AOI_vertices)
            self.NFZ_vertices = [[[1050,950], [1050, 1550], [1450, 1450], [1450, 1050]],
                                  [[1350,2090], [1650, 1800], [2000, 2200], [1450, 2600]],
                                  [[2550,2550], [2550, 2950], [2950, 2950], [2950, 2550]],
                                  [[550,2550], [550, 2950], [950, 2950], [950, 2550]],
                                  [[2550,1550], [2550, 1950], [2950, 1950], [2950, 1550]],
                                  [[3550,1050], [4050, 1150], [3750, 1550]],
                                  [[2050,3050], [2050, 3450], [2450, 3450], [2450, 3050]]]
            self.NFZ_polygons = [geom.polygon.Polygon(vert) for vert in self.NFZ_vertices]
            self.POC_weights = [1,1,9/25,1,1]
            self.POC_means = [[1200,3000], [3000,3500], [2300,2100], [4000,2400], [3100,1100]]
            self.POC_variance_matrices = [[[500**2, 0],[0, 500**2]], [[500**2, 0],[0, 500**2]], [[300**2, 0],[0, 300**2]], [[500**2, 0],[0, 500**2]], [[500**2, 0],[0, 500**2]]]

        else: #randomized mission scenario and POC map.
            random.seed(seed)
            #Randomly generate an AOI
            size_scale = random.uniform(0.2, 1.7)
            
            self.AOI_vertices = generate_random_polygon(center=(2500*size_scale, 2500*size_scale), avg_radius=2200*size_scale, irregularity=0.35, spikiness=0.2, num_vertices=8)
            self.AOI_polygon = geom.polygon.Polygon(self.AOI_vertices)
            self.NFZ_vertices = randomize_NFZs(self.AOI_polygon, size_scale)
            self.NFZ_polygons = [geom.polygon.Polygon(vert) for vert in self.NFZ_vertices]
            
            num_expected_pos = random.randint(3,8)            
            self.POC_weights = [random.randint(80, 120)/100 for _ in range(num_expected_pos)]
            self.POC_means = [[2500*size_scale + random.uniform(-1500*size_scale,1500*size_scale), 2500*size_scale + random.uniform(-1500*size_scale,1500*size_scale)] for _ in range(num_expected_pos)]
            
            correlation_list = [random.uniform(-0.5,0.5) for _ in range(num_expected_pos)]
            POC_sigma_x_list = [random.randint(int(300*size_scale), int(800*size_scale)) for _ in range(num_expected_pos)]
            POC_sigma_y_list = [random.randint(int(300*size_scale), int(800*size_scale)) for _ in range(num_expected_pos)]
            self.POC_variance_matrices = [[[POC_sigma_x**2, correlation*POC_sigma_x*POC_sigma_y],[correlation*POC_sigma_x*POC_sigma_y, POC_sigma_y**2]] for (correlation, POC_sigma_x, POC_sigma_y) in zip(correlation_list, POC_sigma_x_list, POC_sigma_y_list)]
            
        #sensor data from a DJI phantom 4 pro
        self.hfov = hfov #field of view in degrees
        self.h_sensor = h_sensor #meters
        self.p_o = p_o #percentage of overlap
        
    
    def generate_grid(self, theta = 0, sx = 0, sy = 0):
        #grid size calculation
        self.d = 2*(1-self.p_o/100)*self.h_sensor*math.tan(self.hfov/2) #as per eq 3.1 in PAer report
        
        #Calculate the CG of the grid and its size
        CG_grid, m = define_initial_grid(self.AOI_vertices, self.d)
        
        #Calculate the centers of the grid cells
        self.cell_centers = calculate_grid_cells_centers_shifted_and_rotated(CG_grid, m, self.d, sx, sy, theta)
        
        #self.AOI_polygon = geom.polygon.Polygon(self.AOI_vertices)
        #self.NFZ_polygons = [geom.polygon.Polygon(vert) for vert in self.NFZ_vertices]
        
        #Check if the points are inside the AOI
        for cell in self.cell_centers:
            point = geom.Point(cell.x, cell.y)
            cell.status = self.AOI_polygon.contains(point) #A cell is valid if its center is inside the AOI
        
        #Calculate the POC map
        self.cell_centers = calculate_poc(self.d, self.cell_centers, self.POC_weights, self.POC_means, self.POC_variance_matrices, theta)
        
        r = math.sqrt(2)*self.d/2
        #Check if the points are inside the NFZ
        for cell in self.cell_centers:
            square = [[cell.x+r*math.cos(45/180*math.pi + theta), cell.y+r*math.sin(45/180*math.pi + theta)], [cell.x+r*math.cos(135/180*math.pi + theta), cell.y+r*math.sin(135/180*math.pi + theta)], [cell.x+r*math.cos(225/180*math.pi + theta), cell.y+r*math.sin(225/180*math.pi + theta)], [cell.x+r*math.cos(315/180*math.pi + theta), cell.y+r*math.sin(315/180*math.pi + theta)]]
            cell_polygon = geom.polygon.Polygon(square)
            for poly in self.NFZ_polygons:
                if cell_polygon.intersects(poly):
                    cell.status = 0
                    break
                
        #Calculate the POC map
        #self.cell_centers = calculate_poc_v2(self.cell_centers, self.POC_x0_list, self.POC_y0_list, self.POC_sigma_x_list, self.POC_sigma_y_list)
        
        #total_poc = sum([cell.poc for cell in self.cell_centers if cell.status == 1])
        #print("Total poc:", total_poc)
        
        return copy.deepcopy(self.cell_centers), sx, sy, theta

    def generate_optimal_grid_bayesian_opt(self, n_trials = 100):
        #grid size calculation
        self.d = 2*(1-self.p_o/100)*self.h_sensor*math.tan(self.hfov/2) #as per eq 3.1 in PAer report
        
        #Calculate the CG of the grid and its size
        CG_grid, m = define_initial_grid(self.AOI_vertices, self.d)
        
        now = datetime.now()
        fig_name = now.strftime("%d_%m_%Y_%Hh_%Mm_%Ss")
        
        start = time.time()
        model = black_box_function_AOI(self.d, CG_grid, m, self.AOI_polygon, self.NFZ_polygons, self.POC_weights, self.POC_means, self.POC_variance_matrices)#, self.POC_x0_list, self.POC_y0_list, self.POC_sigma_x_list, self.POC_sigma_y_list)
        
        
        model.fig_name = fig_name
        scenario = Scenario(model.configspace, deterministic=True, n_trials = n_trials)#, n_trials = 100, n_workers = 3)
        # Use SMAC to find the best configuration/hyperparameters
        smac = HyperparameterOptimizationFacade(
            scenario,
            model.function,  #model.function,  # We pass the target function here
            overwrite=True,  # Overrides any previous results that are found that are inconsistent with the meta-data
            logging_level = False,
        )
        
        info = TrialInfo(config=model.configspace.get_default_configuration(), seed = 0)
        aux_time_start = time.time()
        cost = model.function(model.configspace.get_default_configuration(), seed = 0)
        aux_time_end = time.time()
        value = TrialValue(cost=cost, time=aux_time_end-aux_time_start)
        smac.tell(info, value)
        
        incumbent = smac.optimize()
        
        end = time.time()
        
        default_cost = smac.validate(model.configspace.get_default_configuration(), seed = 0)
        config_history = smac.runhistory.get_configs()
        cost_history = []
        for config in config_history:
            cost_history.append(smac.runhistory.get_cost(config))
            
        
        if 0:
            f = open('./bayesian_logs/AOI_' + fig_name+'.txt', "a")
            #f = open('./bayesian_logs/SA_' + fig_name+'.txt', "a")
            f.write("Wall clock time (s): " + str(end-start) + "\n")
            f.write(f"Default config: {model.configspace.get_default_configuration()}\n")
            f.write(f"Default cost: {default_cost}\n")
            f.write(f"Incumbent config: {incumbent}\n")
            f.write(f"Incumbent cost: {smac.runhistory.get_cost(incumbent)}\n")
            f.write("\nRun History\n")
            for (config, cost) in zip(config_history, cost_history):
                f.write(str(config)+" Cost: "+str(cost)+"\n")
            f.close()
        
        
        #Calculate the centers of the grid cells
        #self.cell_centers = calculate_grid_cells_centers(CG_grid, m, self.d)
        sx = incumbent["sx"]
        sy = incumbent["sy"]
        theta = incumbent["theta"]
        self.cell_centers = calculate_grid_cells_centers_shifted_and_rotated(CG_grid, m, self.d, sx, sy, theta)
        
        #self.AOI_polygon = geom.polygon.Polygon(self.AOI_vertices)
        #self.NFZ_polygons = [geom.polygon.Polygon(vert) for vert in self.NFZ_vertices]
        
        #Check if the points are inside the AOI
        for cell in self.cell_centers:
            point = geom.Point(cell.x, cell.y)
            cell.status = self.AOI_polygon.contains(point) #A cell is valid if its center is inside the AOI
        
        #Calculate the POC map
        self.cell_centers = calculate_poc(self.d, self.cell_centers, self.POC_weights, self.POC_means, self.POC_variance_matrices, theta)
        
        r = math.sqrt(2)*self.d/2
        #Check if the points are inside the NFZ
        for cell in self.cell_centers:
            square = [[cell.x+r*math.cos(45/180*math.pi + theta), cell.y+r*math.sin(45/180*math.pi + theta)], [cell.x+r*math.cos(135/180*math.pi + theta), cell.y+r*math.sin(135/180*math.pi + theta)], [cell.x+r*math.cos(225/180*math.pi + theta), cell.y+r*math.sin(225/180*math.pi + theta)], [cell.x+r*math.cos(315/180*math.pi + theta), cell.y+r*math.sin(315/180*math.pi + theta)]]
            cell_polygon = geom.polygon.Polygon(square)
            for poly in self.NFZ_polygons:
                if cell_polygon.intersects(poly):
                    cell.status = 0
                    break
        
        total_poc = sum([cell.poc for cell in self.cell_centers if cell.status == 1])
        #print("Total poc:", total_poc, "  sx: ", sx, "m  sy: ", sy, "m  theta: ", theta/math.pi*180,"deg")
        
        return copy.deepcopy(self.cell_centers), sx, sy, theta

class black_box_function_AOI:
    def __init__ (self, d, CG_grid, m, AOI_polygon, NFZ_polygons, POC_weights, POC_means, POC_variance_matrices):
        self.d = d
        self.CG_grid = CG_grid
        self.m = m
        self.AOI_polygon = AOI_polygon
        self.NFZ_polygons = NFZ_polygons
        self.POC_weights = POC_weights
        self.POC_means = POC_means
        self.POC_variance_matrices = POC_variance_matrices
        
        
    @property
    def configspace(self) -> ConfigurationSpace:
        cs = ConfigurationSpace(seed=0)
        #'initial_temp': (0.0001, 0.1), 'final_temp': (1e-8, 0.0001), 'cooling_factor': (0.8, 0.96)
        #'initial_temp': 0.0004, 'final_temp': 2.755e-06, 'cooling_factor': 0.96
        #x = Float("x", (-5, 5), default=-5)
        #cs.add([x])
        sx = Float("sx", (-self.d/2, self.d/2), default=0)
        cs.add([sx])
        sy = Float("sy", (-self.d/2, self.d/2), default=0)
        cs.add([sy])
        theta = Float("theta", (-math.pi/4, math.pi/4), default=0)
        cs.add([theta])
        return cs

    def function(self, config: Configuration, seed) -> float:
        
        #Calculate the centers of the grid cells
        sx = config["sx"]
        sy = config["sy"]
        theta = config["theta"]
        cell_centers = []
        cell_centers = calculate_grid_cells_centers_shifted_and_rotated(self.CG_grid, self.m, self.d, sx, sy, theta)
        
        #Check if the points are inside the AOI
        for cell in cell_centers:
            point = geom.Point(cell.x, cell.y)
            cell.status = self.AOI_polygon.contains(point) #A cell is valid if its center is inside the AOI
        
        r = math.sqrt(2)*self.d/2
        #Check if the points are inside the NFZ
        for cell in cell_centers:
            square = [[cell.x+r*math.cos(45/180*math.pi + theta), cell.y+r*math.sin(45/180*math.pi + theta)], [cell.x+r*math.cos(135/180*math.pi + theta), cell.y+r*math.sin(135/180*math.pi + theta)], [cell.x+r*math.cos(225/180*math.pi + theta), cell.y+r*math.sin(225/180*math.pi + theta)], [cell.x+r*math.cos(315/180*math.pi + theta), cell.y+r*math.sin(315/180*math.pi + theta)]]
            cell_polygon = geom.polygon.Polygon(square)
            for poly in self.NFZ_polygons:
                if cell_polygon.intersects(poly):
                    cell.status = 0
                    break
        
        #Calculate the POC map
        cell_centers = calculate_poc(self.d, cell_centers, self.POC_weights, self.POC_means, self.POC_variance_matrices, theta)
        
        total_poc = sum([cell.poc for cell in cell_centers if cell.status == 1])
        
        return -total_poc

def generate_random_polygon(center: Tuple[float, float], avg_radius: float,
                     irregularity: float, spikiness: float,
                     num_vertices: int) -> List[Tuple[float, float]]:
    """
    Start with the center of the polygon at center, then creates the
    polygon by sampling points on a circle around the center.
    Random noise is added by varying the angular spacing between
    sequential points, and by varying the radial distance of each
    point from the centre.

    Args:
        center (Tuple[float, float]):
            a pair representing the center of the circumference used
            to generate the polygon.
        avg_radius (float):
            the average radius (distance of each generated vertex to
            the center of the circumference) used to generate points
            with a normal distribution.
        irregularity (float):
            variance of the spacing of the angles between consecutive
            vertices.
        spikiness (float):
            variance of the distance of each vertex to the center of
            the circumference.
        num_vertices (int):
            the number of vertices of the polygon.
    Returns:
        List[Tuple[float, float]]: list of vertices, in CCW order.
    """
    # Parameter check
    if irregularity < 0 or irregularity > 1:
        raise ValueError("Irregularity must be between 0 and 1.")
    if spikiness < 0 or spikiness > 1:
        raise ValueError("Spikiness must be between 0 and 1.")

    irregularity *= 2 * math.pi / num_vertices
    spikiness *= avg_radius
    angle_steps = random_angle_steps(num_vertices, irregularity)

    # now generate the points
    points = []
    angle = random.uniform(0, 2 * math.pi)
    for i in range(num_vertices):
        radius = clip(random.gauss(avg_radius, spikiness), 0, 2 * avg_radius)
        point = (center[0] + radius * math.cos(angle),
                 center[1] + radius * math.sin(angle))
        points.append(point)
        angle += angle_steps[i]

    return points

def random_angle_steps(steps: int, irregularity: float) -> List[float]:
    """Generates the division of a circumference in random angles.

    Args:
        steps (int):
            the number of angles to generate.
        irregularity (float):
            variance of the spacing of the angles between consecutive vertices.
    Returns:
        List[float]: the list of the random angles.
    """
    # generate n angle steps
    angles = []
    lower = (2 * math.pi / steps) - irregularity
    upper = (2 * math.pi / steps) + irregularity
    cumsum = 0
    for i in range(steps):
        angle = random.uniform(lower, upper)
        angles.append(angle)
        cumsum += angle

    # normalize the steps so that point 0 and point n+1 are the same
    cumsum /= (2 * math.pi)
    for i in range(steps):
        angles[i] /= cumsum
    return angles

def clip(value, lower, upper):
    """
    Given an interval, values outside the interval are clipped to the interval
    edges.
    """
    return min(upper, max(value, lower))

def randomize_NFZs(AOI_polygon, size_scale):
    NFZ_vertices = []
    NFZ_polygons = []
    number_NFZs = random.randint(3, 6)
    
    while len(NFZ_vertices) < number_NFZs:
        vertices=[]
        center = (2500*size_scale + random.uniform(-1500*size_scale,1500*size_scale), 2500*size_scale + random.uniform(-1500*size_scale,1500*size_scale))
        #center = (random.uniform(100/size_scale,4900/size_scale), random.uniform(100/size_scale,4900/size_scale))
        avg_radius = random.uniform(200*size_scale,500*size_scale)
        num_vertices = random.randint(3,6)
        vertices = generate_random_polygon(center=center, avg_radius=avg_radius, irregularity=0.35, spikiness=0.2, num_vertices=num_vertices)
        candidate_NFZ = geom.polygon.Polygon(vertices)
        if candidate_NFZ.intersects(AOI_polygon) and all(not candidate_NFZ.intersects(nfz_pol) for nfz_pol in NFZ_polygons):
            NFZ_vertices.append(copy.deepcopy(vertices))
            NFZ_polygons.append(geom.polygon.Polygon(vertices))
    
    return NFZ_vertices

class InputArea():
    def __init__(self, type):
        self.type = type # 0 - AOI, 1 - NFZ
        
        #defined in meters in the local NED ref
        if type == 0:
            self.vertices = [[50, 50], [1525, 500], [2400, 1250], [1200, 2350],  [50, 1500]] 
        elif type == 1:
            self.vertices = [[550,550], [550, 950], [950, 950], [950, 550]]
        
        self.polygon = geom.polygon.Polygon(self.vertices)
        
    def is_point_inside(self, x_pos, y_pos):
        """Checks if the point x, y is inside the area. Returns 0 if no and 1 if yes"""
        point = geom.Point(x_pos ,y_pos)
        return self.polygon.contains(point)

class InputArea_v2():
    def __init__(self, type):
        self.type = type # 0 - AOI, 1 - NFZ
        
        #defined in meters in the local NED ref
        if type == 0:
            #self.vertices = [[50, 50], [3000, 1000], [4800, 2500], [2400, 4700],  [50, 3000]]
            self.vertices = [[50, 1000], [1500, 50], [4800, 1000], [4000, 3500], [2000, 4000], [50, 3000]]
            self.polygon = [geom.polygon.Polygon(self.vertices)]
        elif type == 1:
            self.vertices_list = [[[1050,950], [1050, 1550], [1450, 1450], [1450, 1050]],
                                  [[1350,2090], [1650, 1800], [2000, 2200], [1450, 2600]],
                                  [[2550,2550], [2550, 2950], [2950, 2950], [2950, 2550]],
                                  [[550,2550], [550, 2950], [950, 2950], [950, 2550]],
                                  [[2550,1550], [2550, 1950], [2950, 1950], [2950, 1550]],
                                  [[3550,1050], [4050, 1150], [3750, 1550]],
                                  [[2050,3050], [2050, 3450], [2450, 3450], [2450, 3050]]]
            self.polygon = [geom.polygon.Polygon(vert) for vert in self.vertices_list]
        
        
    def is_point_inside(self, x_pos, y_pos):
        """Checks if the point x, y is inside the area. Returns 0 if no and 1 if yes"""
        point = geom.Point(x_pos ,y_pos)
        for poly in self.polygon:
            if poly.contains(point):
                return True
        return False
        #return self.polygon.contains(point)

class Cell():
    """Defines a cell, its position, coordinates and status"""
    def __init__(self, i_index, j_index, x_coord, y_coord, status):
        self.i = i_index
        self.j = j_index
        self.x = x_coord
        self.y = y_coord
        self.status = status
        self.poc = 0
        #self.not_covered = 1
        self.times_covered_prev = 0
        self.prev_action = 0
    
    def __eq__(self, cell):
        """ Returns True if all the attributes in both states are the same """
        return (self.i, self.j, self.times_covered_prev) == (cell.i, cell.j, cell.times_covered_prev)
    
    def __hash__(self):
        """ Defines a state object as a hashable object """
        return hash((self.i,self.j))
    
    def __lt__(self, cell):
        """ Defines the comparison between state objects.
        The "better" state is the one with mode requests
        already dropped off """
        return self.i > cell.i

def search_for_cell_ij(cell_list, i,j):
    for count, cell in enumerate(cell_list):
        if cell.i == i and cell.j == j:
            return count, cell
    return None, None

def define_initial_grid(list_of_coords, d, for_plots = False):
    x_min = min(vertex[0] for vertex in list_of_coords)
    y_min = min(vertex[1] for vertex in list_of_coords)
    x_max = max(vertex[0] for vertex in list_of_coords)
    y_max = max(vertex[1] for vertex in list_of_coords)
    
    if not for_plots:
        #eq 3.4 in PEAer report
        m = math.ceil(math.sqrt((x_max-x_min)**2+(y_max-y_min)**2)/d) + 1
    else:
        m = math.ceil(max([(x_max-x_min), (y_max-y_min)])/d) + 1
    
    #eq 3.5 in PEAer report
    CG_x = (x_max+x_min)/2
    CG_y = (y_max+y_min)/2
    
    return [CG_x,CG_y], m

def calculate_grid_cells_centers(CG_grid, m, d):
    #without considering shift and rotation
    #using eq 3.6 and 3.7
    #note that i is for rows, j is for columns and (0,0) is in the bottom left corner
    cell_centers = []
    for i in range(0,m):
        for j in range(0,m):
            x = CG_grid[0] + (j-(m-1)/2)*d
            y = CG_grid[1] + (i-(m-1)/2)*d
            
            cell_centers.append(Cell(i,j,x,y,0))
    return cell_centers

def calculate_grid_cells_centers_shifted_and_rotated(CG_grid, m, d, sx, sy, theta):
    #without considering shift and rotation
    #using eq 3.10 and 3.11
    #note that i is for rows, j is for columns and (0,0) is in the bottom left corner
    cell_centers = []
    for i in range(0,m):
        for j in range(0,m):
            x = CG_grid[0] + sx + (j-(m-1)/2)*d*math.cos(theta) - (i-(m-1)/2)*d*math.sin(theta)
            y = CG_grid[1] + sy + (j-(m-1)/2)*d*math.sin(theta) + (i-(m-1)/2)*d*math.cos(theta)
            
            cell_centers.append(Cell(i,j,x,y,0))
    return cell_centers

def classify_centers(cell_centers, area, d):
    """Classifies all points in the dictionary as to their location relative to the defined areas"""
    if area.type == 0: #means the area is an AOI
        for cell in cell_centers:
            #treats the outside of the AOI as a NFZ
            cell.status = area.is_point_inside(cell.x, cell.y)
            
    if area.type == 1: #means the area is a NFZ
        for cell in cell_centers:
            for i in [-1, 1]:
                for j in [-1, 1]:
                    val = area.is_point_inside(cell.x+d/2*i, cell.y+d/2*j)
                    val = 1 - val #flips a 0 to 1 and a 1 to 0
                    cell.status = cell.status*val
                    
    return cell_centers

def check_if_cells_adjacent(cell1, cell2):
    """checks for the 4-connected neighbourhood"""
    if abs((cell1.i - cell2.i)) + abs((cell1.j - cell2.j)) <= 1:
        return True
    else:
        return False

def check_if_cells_adjacent_8_conn(cell1, cell2):
    """checks for the 8-connected neighbourhood"""
    if abs((cell1.i - cell2.i)) <= 1 and abs((cell1.j - cell2.j)) <= 1:
        return True
    else:
        return False

def plot_grid(cell_centers, d, ax, sx, sy, theta):
    i = 0
    j = 0
    #size = 1
    size = 0.05#0.1
    r = math.sqrt(2)*d/2 
    for cell in cell_centers:
        if cell.status == 0:
            if i == 0:
                plt.scatter(cell.x, cell.y, marker='o', s=size, color='r', label = 'invalid cell')
                i += 1
            else:
                plt.scatter(cell.x, cell.y, marker='o', s=size, color='r')
            
        elif cell.status == 1:
            if j == 0:
                plt.scatter(cell.x, cell.y, marker='o', s=size, color='g', label = 'valid cell')
                j += 1
            else:
                plt.scatter(cell.x, cell.y, marker='o', s=size, color='g')
        
        square = [[cell.x+r*math.cos(45/180*math.pi + theta), cell.y+r*math.sin(45/180*math.pi + theta)], [cell.x+r*math.cos(135/180*math.pi + theta), cell.y+r*math.sin(135/180*math.pi + theta)], [cell.x+r*math.cos(225/180*math.pi + theta), cell.y+r*math.sin(225/180*math.pi + theta)], [cell.x+r*math.cos(315/180*math.pi + theta), cell.y+r*math.sin(315/180*math.pi + theta)]]
        #square = [(cell.x-d/2, cell.y+d/2),(cell.x+d/2, cell.y+d/2),
        #              (cell.x+d/2, cell.y-d/2),(cell.x-d/2, cell.y-d/2),
        #              (cell.x-d/2, cell.y+d/2)]
        square_patch = ptch.Polygon(square, alpha=0.05, facecolor = 'w',edgecolor = 'k')
        ax.add_patch(square_patch)

def gaussian_func(weight, mean, variance_matrix):
    xy_mean = np.array([mean])
    var = np.array(variance_matrix)
    return lambda y, x: weight * 1/(2*math.pi*math.sqrt(np.linalg.det(variance_matrix))) * (math.exp(-0.5 * (np.array([[x, y]])-xy_mean) @ np.linalg.inv(var) @ (np.array([[x, y]]-xy_mean).T)))

def calculate_poc(d, cell_centers, POC_weights, POC_means, POC_variance_matrices, theta):
    """Considers that the target might be inside a NFZ, and is definitely inside the AOI"""
    #Takes 4 seconds to compute the integral
    gaussian_funcs = []
    for (weight, mean, variance_matrix) in zip(POC_weights, POC_means, POC_variance_matrices):
        #xy_mean = np.array([mean])
        #var = np.array(variance_matrix)
        #gauss_f = lambda x, y: weight * 1/(2*math.pi*math.sqrt(np.linalg.det(variance_matrix))) * (math.exp(-0.5 * (np.array([[x, y]])-xy_mean) @ np.linalg.inv(var) @ (np.array([[x, y]]-xy_mean).T)))
        gauss_f = gaussian_func(weight, mean, variance_matrix)
        gaussian_funcs.append(gauss_f)
    
    for cell in cell_centers:
        if cell.status == 1:
            cell.poc = 0
            for gauss_func in gaussian_funcs:
                x0=cell.x
                y0=cell.y
                def rotated_gaussian_func(y, x):
                    qx = x0 + math.cos(theta) * (x - x0) - math.sin(theta) * (y - y0)
                    qy = y0 + math.sin(theta) * (x - x0) + math.cos(theta) * (y - y0)
                    return gauss_func(qy,qx)
                #(result, err) = integrate.dblquad(rotated_gaussian_func, cell.x-d/2, cell.x+d/2, cell.y-d/2, cell.y+d/2)
                
                x = np.linspace(cell.x-d/2, cell.x+d/2, 5)
                y = np.linspace(cell.y-d/2, cell.y+d/2, 5)
                X, Y = np.meshgrid(x, y)
                vfunc = np.vectorize(rotated_gaussian_func)
                list1=vfunc(Y, X)
                result = np.trapz(np.trapz(list1, y, axis=0), x, axis=0)
                
                cell.poc += result
            #Scaling for the sum of the weights of the elements in the list
            cell.poc = cell.poc/sum(POC_weights)
        else:
            cell.poc = 0
    
    return cell_centers



def plot_poc(cell_centers, d, ax, sx, sy, theta):
    #cmap = ListedColormap(plt.get_cmap('jet')(np.linspace(0.1, 1, 256)))
    aux = [100*cell.poc for cell in cell_centers]
    normal = pl.Normalize(min(aux), max(aux))
    colors = pl.cm.jet(normal(aux))
    
    r = math.sqrt(2)*d/2 
    for count, cell in enumerate(cell_centers):
        
        square = [[cell.x+r*math.cos(45/180*math.pi + theta), cell.y+r*math.sin(45/180*math.pi + theta)], [cell.x+r*math.cos(135/180*math.pi + theta), cell.y+r*math.sin(135/180*math.pi + theta)], [cell.x+r*math.cos(225/180*math.pi + theta), cell.y+r*math.sin(225/180*math.pi + theta)], [cell.x+r*math.cos(315/180*math.pi + theta), cell.y+r*math.sin(315/180*math.pi + theta)]]
        #square = [(cell.x-d/2, cell.y+d/2),(cell.x+d/2, cell.y+d/2),
        #              (cell.x+d/2, cell.y-d/2),(cell.x-d/2, cell.y-d/2),
        #              (cell.x-d/2, cell.y+d/2)]
        if cell.status == 0:
            plt.plot(cell.x, cell.y, marker='o', markersize=0.5, color='r')
            square_patch = ptch.Polygon(square, alpha=0.05, facecolor = 'w',edgecolor = 'k')
            
        elif cell.status == 1:
            plt.plot(cell.x, cell.y, marker='o', markersize=0.5, color='g')
            #square_patch = ptch.Polygon(square, alpha=1, facecolor = my_cmap(cell.poc), edgecolor = 'k')
            square_patch = ptch.Polygon(square, alpha=1, facecolor = colors[count])#, edgecolor = 'k')
            square_patch_edge = ptch.Polygon(square, alpha=0.05, edgecolor = 'k', fill=False, closed=True)
            ax.add_patch(square_patch_edge)
            #polyEdge = Polygon(k*(np.array(pts)- cent) + cent,fill=False,closed=True, alpha=0.8,edgecolor=color,linewidth=7)
        
        ax.add_patch(square_patch)
    
    cax, _ = cbar.make_axes(ax) 
    cb2 = cbar.ColorbarBase(cax, cmap=pl.cm.jet,norm=normal, label = "Probability of Containment (%)")
    #cbar.label('# of contacts')

def calculate_turn(cell1, cell2, cell3):
        """Calculate the degrees involved in turing from cell1, 2 and 3"""
        vec1 = np.array([cell2.x-cell1.x, cell2.y-cell1.y])
        vec2 = np.array([cell3.x-cell2.x, cell3.y-cell2.y])
        cos_angle = np.dot(vec1, vec2)/(np.linalg.norm(vec1)*np.linalg.norm(vec2))
        cos_angle = round(cos_angle,4)
        if (np.linalg.norm(vec1)*np.linalg.norm(vec2)) == 0:
            return 0
        else:
            angle = np.arccos(cos_angle)/np.pi*180
        return angle #in degrees

def calculate_turn_from_actions(action_prev, action):
        """Calculate the degrees involved in turning from cell1, 2 and 3"""
        if action_prev[0] == 'FIRST':
            return 0
        vec1 = np.array([action[1], action[2]])
        vec2 = np.array([action_prev[1], action_prev[2]])
        cos_angle = np.dot(vec1, vec2)/(np.linalg.norm(vec1)*np.linalg.norm(vec2))
        cos_angle = round(cos_angle,4)
        angle = np.arccos(cos_angle)/np.pi*180
        return angle #in degrees

def calculate_info_gain(path):
    order = 0
    decay_factor = 0.01
    objective = 0
    not_prev_visited = 1
    path_cost = 0
    for cell in path:
        if order != 0:
            idx, aux = search_for_cell_ij(path[0:order], cell.i,cell.j)
            if idx == None:
                not_prev_visited = 1
            else:
                not_prev_visited = 0
        #objective += decay_factor*math.exp(-decay_factor*order)*cell.poc*not_prev_visited
        #path_cost += decay_factor*(1-math.exp(-decay_factor*order))*cell.poc*not_prev_visited
        objective += math.exp(-decay_factor*order)*cell.poc*not_prev_visited
        path_cost += (1-math.exp(-decay_factor*order))*cell.poc*not_prev_visited
        order += 1
    
    return objective


def calculate_missed_detections(valid_cell_list, path):
    order = 0
    detection = 0
    not_prev_visited = 1
    for cell_idx in path:
        if order != 0:
            if cell_idx in path[0:order]:
                not_prev_visited = 0
            else:
                not_prev_visited = 1
        detection += valid_cell_list[cell_idx].poc*not_prev_visited
        order += 1
    
    sum_poc = sum([cell.poc for cell in valid_cell_list])
    detection = detection/sum_poc
    return 1-detection

def calculate_detections_multi(valid_cell_list, paths, pod):
    max_path_len = max([len(path) for path in paths])
    times_visited_dict = {i:0 for i in range(len(valid_cell_list))} #initialize the dict as 0s
    total_pos_cumulative = 0
    for step in range(max_path_len):
        for path in paths:
            if step < len(path):
                pos_cumulative = (1-pod)**times_visited_dict[path[step]] * pod * valid_cell_list[path[step]].poc
                total_pos_cumulative += pos_cumulative
                times_visited_dict[path[step]]+=1
    return total_pos_cumulative

def calculate_avg_detection_step(valid_cell_list, path):
    order = 0
    avg_detection_step = 0
    detection = 0
    not_prev_visited = 1
    for cell_idx in path:
        if order != 0:
            if cell_idx in path[0:order]:
                not_prev_visited = 0
            else:
                not_prev_visited = 1
        avg_detection_step += order*valid_cell_list[cell_idx].poc*not_prev_visited
        detection += valid_cell_list[cell_idx].poc*not_prev_visited
        order += 1
    
    avg_detection_step = avg_detection_step/detection
    return avg_detection_step
    
def calculate_avg_detection_step_multi(valid_cell_list, paths, pod):
    max_path_len = max([len(path) for path in paths])
    times_visited_dict = {i:0 for i in range(len(valid_cell_list))} #initialize the dict as 0s
    detect_step_cumulative = 0
    total_pos_cumulative = 0
    for step in range(max_path_len):
        for path in paths:
            if step < len(path):
                pos_cumulative = (1-pod)**times_visited_dict[path[step]] * pod * valid_cell_list[path[step]].poc
                detect_step_cumulative += step*pos_cumulative
                total_pos_cumulative += pos_cumulative
                times_visited_dict[path[step]]+=1
    return detect_step_cumulative/total_pos_cumulative

def calculate_objective_multi(valid_cell_list, paths, pod, decay_factor = 0.01): #paths is a list of indexes
    max_path_len = max([len(path) for path in paths])
    times_visited_dict = {i:0 for i in range(len(valid_cell_list))} #initialize the dict as 0s
    obj = 0
    for step in range(max_path_len):
        for path in paths:
            if step < len(path):
                pos_cumulative = (1-pod)**times_visited_dict[path[step]] * pod * valid_cell_list[path[step]].poc
                obj += math.exp(-decay_factor*step)*pos_cumulative
                times_visited_dict[path[step]]+=1
    return obj
    
def evaluate_paths_multi(paths, cell_centers, pod, save_fig, fig_name):
    """Calculates the overall distance covered, as well as the turns performed"""
    distances = [0 for path in paths]
    turns = [0 for path in paths]
    for uav, path in enumerate(paths):
        for idx in range(len(path)-1):
            distances[uav] = distances[uav] + math.sqrt((path[idx+1].x - path[idx].x)**2 + (path[idx+1].y - path[idx].y)**2)
        
        for idx in range(len(path)-3):
            turns[uav] = turns[uav] + calculate_turn(path[idx], path[idx+1], path[idx+2]) #in degrees7

    turn_weight = 0.0173
    distance_weight = 0.1164
    energies_expended = [turn_weight*turn+distance_weight*distance for turn, distance in zip(turns, distances)]
    print("Paths distance (m): ", distances)
    print("Sum of the turn angles in the path (deg): ", turns)
    print("Energy required for the paths (kJ): ", energies_expended)
    
    valid_cell_list = [cell for cell in cell_centers if cell.status == 1]
    paths_idx = []
    for path in paths:
        aux_list = []
        for cell in path:
            idx_valid_list, aux = search_for_cell_ij(valid_cell_list, cell.i,cell.j)
            aux_list.append(idx_valid_list)
        paths_idx.append(aux_list)
    
    total_poc = sum(cell.poc for cell in valid_cell_list)
    print("Total mission POC: ", total_poc*100, "%")
    
    objective = calculate_objective_multi(valid_cell_list, paths_idx, pod)
    print("Objective function J = ", objective)
    
    detect_perc = calculate_detections_multi(valid_cell_list, paths_idx, pod)
    print("Cumulative probability of success POS_c = ", detect_perc*100, "%")
    
    avg_detect_step = calculate_avg_detection_step_multi(valid_cell_list, paths_idx, pod)
    print("Expected detection time EDT (in UAV steps) = ", avg_detect_step)
    
    
    if save_fig == 1:
        f = open('./results/' + fig_name+'.txt', "a")
        f.write("Paths distance (m): " + str(distances)+'\n')
        f.write("Sum of the turn angles in the path (deg): " + str(turns)+'\n')
        f.write("Energy required for the paths (kJ): " + str(energies_expended)+'\n')
        f.write("Objective function J = " + str(objective)+'\n')
        f.write("Cumulative probability of success POS_c = " + str(detect_perc*100)+'\n')
        f.write("Expected detection time EDT (in UAV steps) = " + str(avg_detect_step)+'\n')
    return objective, detect_perc, avg_detect_step

        
class Graph():
    """From https://www.geeksforgeeks.org/prims-minimum-spanning-tree-mst-greedy-algo-5/"""
    def __init__(self, vertices):
        self.V = vertices
        self.graph = [[0 for column in range(vertices)]
                      for row in range(vertices)]
 
    # A utility function to print 
    # the constructed MST stored in parent[]
    def printMST(self, parent):
        print("Edge \tWeight")
        for i in range(1, self.V):
            print(parent[i], "-", i, "\t", self.graph[i][parent[i]])
 
    # A utility function to find the vertex with
    # minimum distance value, from the set of vertices
    # not yet included in shortest path tree
    def minKey(self, key, mstSet):
 
        # Initialize min value
        min = sys.maxsize
 
        for v in range(self.V):
            if key[v] < min and mstSet[v] == False:
                min = key[v]
                min_index = v
 
        return min_index
 
    # Function to construct and print MST for a graph
    # represented using adjacency matrix representation
    def primMST(self):
 
        # Key values used to pick minimum weight edge in cut
        key = [sys.maxsize] * self.V
        parent = [None] * self.V  # Array to store constructed MST
        # Make key 0 so that this vertex is picked as first vertex
        key[0] = 0
        mstSet = [False] * self.V
 
        parent[0] = -1  # First node is always the root of
 
        for cout in range(self.V):
 
            # Pick the minimum distance vertex from
            # the set of vertices not yet processed.
            # u is always equal to src in first iteration
            u = self.minKey(key, mstSet)
 
            # Put the minimum distance vertex in
            # the shortest path tree
            mstSet[u] = True
 
            # Update dist value of the adjacent vertices
            # of the picked vertex only if the current
            # distance is greater than new distance and
            # the vertex in not in the shortest path tree
            for v in range(self.V):
 
                # graph[u][v] is non zero only for adjacent vertices of m
                # mstSet[v] is false for vertices not yet included in MST
                # Update the key only if graph[u][v] is smaller than key[v]
                if self.graph[u][v] > 0 and mstSet[v] == False \
                and key[v] > self.graph[u][v]:
                    key[v] = self.graph[u][v]
                    parent[v] = u
 
        #self.printMST(parent)
        return parent

class TrajectoryProblem (search.Problem):

    def __init__(self, valid_cell_list, ini_cell, goal_cell):
        """ Problem definition """
        self.valid_cell_list = valid_cell_list
        self.initial = ini_cell
        self.goal = goal_cell
    
    
    #Calculate the cost associated with the solution presented in the plan file
    def cost(self, sol):
        """ Calculates the cost of a solution """
        return 

    def result(self,state,action) :
        """ Return the state that results from executing
        the given action in the given state. """
        #We assume any action here taken is valid
        idx, next_cell = search_for_cell_ij(self.valid_cell_list, state.i+action[1],state.j+action[2])
        next_cell.prev_action = action
        return next_cell
    

    def actions(self,state):
        """ Return the actions that can be executed in
        the given state. """
        actions = []
        for cell in self.valid_cell_list:
            if cell.i == state.i-1 and cell.j == state.j-1:
                actions.append(('SW', -1, -1))
            elif cell.i == state.i-1 and cell.j == state.j:
                actions.append(('S', -1, 0))
            elif cell.i == state.i-1 and cell.j == state.j+1:
                actions.append(('SE', -1, 1))
            elif cell.i == state.i and cell.j == state.j+1:
                actions.append(('E', 0, 1))
            elif cell.i == state.i+1 and cell.j == state.j+1:
                actions.append(('NE', 1, 1))
            elif cell.i == state.i+1 and cell.j == state.j:
                actions.append(('N', 1, 0))
            elif cell.i == state.i and cell.j == state.j-1:
                actions.append(('W', 0, -1))
        return actions
    
    def goal_test(self,state):
        """ Return True if the state is a goal. """
        return state == self.goal
    
    def calculate_turn(self, action_prev, action):
        """Calculate the degrees involved in turing from cell1, 2 and 3"""
        if action_prev[0] == 'FIRST':
            return 0
        vec1 = np.array([action[1], action[2]])
        vec2 = np.array([action_prev[1], action_prev[2]])
        cos_angle = np.dot(vec1, vec2)/(np.linalg.norm(vec1)*np.linalg.norm(vec2))
        cos_angle = round(cos_angle,4)
        angle = np.arccos(cos_angle)/np.pi*180
        return angle #in degrees
    
    def path_cost(self, c, state1, action, state2):
        """Return the cost of a solution path that arrives at state2 from
        state1 via action, assuming cost c to get up to state1"""
        #We can use for the cost the euclidean distance between cells
        distance_cost = math.sqrt((state2.x - state1.x)**2 + (state2.y - state1.y)**2)
        #We can also consider the turns between cells
        #turn_cost = self.calculate_turn(state1.prev_action, action)
        #In order to approximate the path cost to the energy expenditure, using the paper:
        #Modares, et al. UB-ANC planner: Energy efficient coverage path planning with multiple drones
        turn_weight = 0.0173
        distance_weight = 0.1164 #the test was done for 10 m/s, which is quite low
        #turn_weight = 0
        #distance_weight = 0
        #include a term for information gain
        #info_cost = 1-state2.poc*state2.not_covered #the greater the POC of the cells visited, the lower the cost
        #info_weight = 100 #tunable parameter
        #info_weight = 0
        total_cost = c + distance_weight*distance_cost# + turn_weight*turn_cost# + info_weight*info_cost #Results in the total cost in KJ
        return total_cost
    
    
    def h(self, state):
        """ Return the heuristic value for the given state """
        #We can use for an heuristic the euclidean distance between cells
        return 0.1164*math.sqrt((state.state.x - self.goal.x)**2 + (state.state.y - self.goal.y)**2)


    def solve(self):
        """ Calls the uninformed search algorithm
        chosen. Returns a solution using the specified format. """
        # The best method found was astar search
        #last_node = search.depth_first_graph_search(self)
        #last_node = search.uniform_cost_search(self)
        last_node = search.astar_search(self)
        #self.test_solution(last_node)
        return last_node.solution(), last_node.path()
    
    def test_solution(self, last_node):
        """ Useful function for debugging - prints the path until the given node """
        path = last_node.path()
        for node in path:
            print("Action: ", node.action)
            print(node.state.disp())
            print("path cost: ", node.path_cost, "h: ", self.h(node), "f(n): ", node.path_cost + self.h(node))
            print('\n')
            
class State():
    def __init__(self, valid_cell_list, curr_cell, initial_energy):
        self.valid_cell_list = tuple(valid_cell_list) #has the original POC map calculated
        self.curr_cell = curr_cell
        self.steps = 0
        self.energy_remaining = initial_energy #kJ
        #idx, cell = search_for_cell_ij(self.valid_cell_list, curr_cell.i,curr_cell.j)
        #self.curr_cell_idx = idx
    
    def __eq__(self, state):
        """ Returns True if all the attributes in both states are the same """
        
    
        if self.valid_cell_list == state.valid_cell_list and self.curr_cell.i == state.curr_cell.i and self.curr_cell.j == state.curr_cell.j and self.steps == state.steps and self.energy_remaining == state.energy_remaining:
            return True
        else:
            return False
        #return all((cell1.not_covered == cell2.not_covered) for cell1 in self.valid_cell_list for cell2 in state.valid_cell_list if cell1==cell2)
    
    def __hash__(self):
        """ Defines a state object as a hashable object """
        return hash((self.valid_cell_list, self.curr_cell))
    
    def __lt__(self, state):
        """ Defines the comparison between state objects."""
        return len([1 for cell in self.valid_cell_list if cell.not_covered == 0]) > len([1 for cell in state.valid_cell_list if cell.not_covered == 0])
    
    def display(self):
        for cell in self.valid_cell_list:
            print("i = ", cell.i, "j = ", cell.j, "not_covered = ", cell.not_covered)
        print("Current_Cell: ", "i = ", self.curr_cell.i, "j = ", self.curr_cell.j, "not_covered = ", self.curr_cell.not_covered)
        
class FullTrajectoryProblem (search.Problem):

    def __init__(self, valid_cell_list, initial_cell, d, initial_energy):
        """ Problem definition """
        self.initial = State(valid_cell_list, initial_cell, initial_energy)
        self.d = d
        #goal = copy.deepcopy(valid_cell_list)
        #for idx in range(len(goal)):
        #    goal[idx].not_covered = 0
        
        #self.goal = State(goal, initial_cell) #the initial cell is not relevant here
    
    
    #Calculate the cost associated with the solution presented in the plan file
    def cost(self, sol):
        """ Calculates the cost of a solution """
        return 

    def copy_state_to_editable(self, state):
        """ Make a copy of the state into a version that can be edited """
        next_state = copy.deepcopy(state) # we use copy otherwise they keep the same pointer
        next_state.valid_cell_list = [x for x in next_state.valid_cell_list]
        return next_state
    
    def convert_state_to_hashable(self, state):
        """ Convert the state into a hashable version """
        state.valid_cell_list = tuple(state.valid_cell_list)
        return state
    
    def result(self,state,action) :
        """ Return the state that results from executing
        the given action in the given state. """
        #We assume any action here taken is valid
        idx, next_cell = search_for_cell_ij(state.valid_cell_list, state.curr_cell.i+action[1],state.curr_cell.j+action[2])
        
        next_state = self.copy_state_to_editable(state)
        next_state.curr_cell = copy.deepcopy(next_cell)
        next_state.curr_cell.prev_action = action
        next_state.valid_cell_list[idx].not_covered = 0
        
        #If using the information heuristic
        next_state.steps += 1
        
        #If considering the energy remaining
        distance_cost = math.sqrt((next_cell.x - state.curr_cell.x)**2 + (next_cell.y - state.curr_cell.y)**2)
        turn_cost = self.calculate_turn(state.curr_cell.prev_action, action)
        turn_weight = 0.0173
        distance_weight = 0.1164
        
        next_state.energy_remaining = next_state.energy_remaining - distance_weight*distance_cost - turn_weight*turn_cost
        
        next_state = self.convert_state_to_hashable(next_state)
        
        if next_state.energy_remaining <= 0:
            aux_state = self.copy_state_to_editable(state)
            aux_state.energy_remaining = next_state.energy_remaining
            aux_state = self.convert_state_to_hashable(aux_state)
            return aux_state
        
        return next_state
    

    def actions(self,state):
        """ Return the actions that can be executed in
        the given state. """
        curr_cell = state.curr_cell
        actions = []
        for cell in state.valid_cell_list:
            if cell.i == curr_cell.i-1 and cell.j == curr_cell.j-1:
                actions.append(('SW', -1, -1))
            elif cell.i == curr_cell.i-1 and cell.j == curr_cell.j:
                actions.append(('S', -1, 0))
            elif cell.i == curr_cell.i-1 and cell.j == curr_cell.j+1:
                actions.append(('SE', -1, 1))
            elif cell.i == curr_cell.i and cell.j == curr_cell.j+1:
                actions.append(('E', 0, 1))
            elif cell.i == curr_cell.i+1 and cell.j == curr_cell.j+1:
                actions.append(('NE', 1, 1))
            elif cell.i == curr_cell.i+1 and cell.j == curr_cell.j:
                actions.append(('N', 1, 0))
            elif cell.i == curr_cell.i+1 and cell.j == curr_cell.j-1:
                actions.append(('NW', 1, -1))    
            elif cell.i == curr_cell.i and cell.j == curr_cell.j-1:
                actions.append(('W', 0, -1))
        return actions
    
    def goal_test(self,state):
        """ Return True if the state is a goal. """
        if state.energy_remaining <= 0:
            return True
        
        for cell in state.valid_cell_list:
            if cell.not_covered==1:
                return False
        return True
    
        count = 0
        for cell in state.valid_cell_list:
            if cell.not_covered==0:
                count = count+1
            if count >=8:
                return True
        return False
        #return all(cell.not_covered==0 for cell in state.valid_cell_list)
    
    def calculate_turn(self, action_prev, action):
        """Calculate the degrees involved in turning from cell1, 2 and 3"""
        if action_prev[0] == 'FIRST':
            return 0
        vec1 = np.array([action[1], action[2]])
        vec2 = np.array([action_prev[1], action_prev[2]])
        cos_angle = np.dot(vec1, vec2)/(np.linalg.norm(vec1)*np.linalg.norm(vec2))
        cos_angle = round(cos_angle,4)
        angle = np.arccos(cos_angle)/np.pi*180
        return angle #in degrees
    
    #Using for the path cost the energy cost (without turns) to traverse the map
    def path_cost1(self, c, state1, action, state2):
        """Return the cost of a solution path that arrives at state2 from
        state1 via action, assuming cost c to get up to state1"""
        #We can use for the cost the euclidean distance between cells
        distance_cost = math.sqrt((state2.curr_cell.x - state1.curr_cell.x)**2 + (state2.curr_cell.y - state1.curr_cell.y)**2)
        #We can also consider the turns between cells
        turn_cost = self.calculate_turn(state1.curr_cell.prev_action, action)
        #In order to approximate the path cost to the energy expenditure, using the paper:
        #Modares, et al. UB-ANC planner: Energy efficient coverage path planning with multiple drones
        #turn_weight = 0.0173
        turn_weight = 0
        distance_weight = 0.1164 #the test was done for 10 m/s, which is quite low
        total_cost = c + distance_weight*distance_cost + turn_weight*turn_cost #Results in the total cost in KJ
        return total_cost
    
    def find_closest_cell_dist(self, group1, group2):
        min_dist = None
        for cell1 in group1:
            for cell2 in group2:
                distance = math.sqrt((cell2.x - cell1.x)**2 + (cell2.y - cell1.y)**2)
                if min_dist == None or distance<min_dist:
                    min_dist = distance
        return min_dist

    def find_furthest_cell_dist(self, group1, group2):
        max_dist = None
        for cell1 in group1:
            for cell2 in group2:
                distance = math.sqrt((cell2.x - cell1.x)**2 + (cell2.y - cell1.y)**2)
                if max_dist == None or distance>max_dist:
                    max_dist = distance
        return max_dist
    
    def get_isolated_group(self, not_visited_cell_list, cell):
        group = []
        idx, curr_cell = search_for_cell_ij(not_visited_cell_list, cell.i,cell.j)
        group.append(curr_cell)
        i=0
        while(i<len(not_visited_cell_list)):
            cell2 = not_visited_cell_list[i]
            i = i+1
            if cell2 in group:
                pass
            else:
                for cell1 in group:
                    if check_if_cells_adjacent_8_conn(cell1, cell2):
                        group.append(cell2)
                        i = 0
                        break
        return group
                    
        
    def estimate_dist(self, valid_cell_list, cell):
        min_dist = None
        for cell2 in valid_cell_list:
            if cell.i != cell2.i and cell.j != cell2.j and cell2.not_covered == 1:
                distance = math.sqrt((cell2.x - cell.x)**2 + (cell2.y - cell.y)**2)
                if min_dist == None or distance<min_dist:
                    min_dist = distance
                    #we can check if it already the same as d to speed it up
                    if abs(cell.i-cell2.i)+abs(cell.j-cell2.j) == 1:
                        break
        if min_dist == None:
            return 0
        return min_dist
    
    #Using for the heuristic the expected remaining energy cost to traverse the map    
    def h1(self, node):
        """ Return the heuristic value for the given state """
        remaining_distance_cost = 0
        distance_weight = 0.1164
        
        for cell in node.state.valid_cell_list:
            if cell.not_covered == 1:
                #Calculate distance to closest not covered cell since optimally, they will not go over already covered cells
                remaining_distance_cost += distance_weight*self.estimate_dist(node.state.valid_cell_list, cell)
        #Find the closest unvisited cell to the current one
        remaining_distance_cost += distance_weight*self.estimate_dist(node.state.valid_cell_list, node.state.curr_cell)
        
        #Account for spread out groups
        not_visited_cell_list = [cell for cell in node.state.valid_cell_list if cell.not_covered == 1]
        
        if len(not_visited_cell_list)==0:
            return 0
        
        while(1):
            #Given a random cell, calculate the ones adjacent to it
            initial_cell = random.choice(not_visited_cell_list)
            isolated_group = self.get_isolated_group(not_visited_cell_list, initial_cell)
            
            #Calculate the distance to closest cell outside this group
            if len(isolated_group) < len(not_visited_cell_list):
                #remove the elements of the isolated group from the not visited list
                not_visited_cell_list = [cell for cell in not_visited_cell_list if cell not in isolated_group]
                distance_between_groups = self.find_closest_cell_dist(isolated_group, not_visited_cell_list)
                remaining_distance_cost += distance_weight*distance_between_groups
            else:
                break #There are no isolated cells
        
        return remaining_distance_cost
    
    def info_penalty(self, valid_cell_list):
        factor = 0.1
        num_prev_visited = 0
        for cell in valid_cell_list:
            if cell.not_covered == 0:
                num_prev_visited += 1
        return math.exp(-factor*num_prev_visited)
        
    def calculate_dist(self, cell1, cell2):
        return math.sqrt((cell2.x - cell1.x)**2 + (cell2.y - cell1.y)**2)
    
    #Using for the path cost the information not gained
    def path_cost(self, c, state1, action, state2):
        factor = 0.01
        #info_loss = (1-self.info_penalty(state1.valid_cell_list))*state2.curr_cell.poc
        if state2.energy_remaining <= 0:
            info_loss = 0
            for cell in state1.valid_cell_list:
                if cell.not_covered == 1:
                    info_loss += factor*1*cell.poc
        else:
            if state2.curr_cell.not_covered == 1:
                info_loss = factor*(1-math.exp(-factor*state1.steps))*state2.curr_cell.poc
            else:
                info_loss = 0
        return c + info_loss

    #Using for the path cost the estimate of the information not gained
    def h2(self, node):
        """ Return the heuristic value for the given state """
        #Get the list of cells that havent been visited
        not_visited_cell_list = [cell for cell in node.state.valid_cell_list if cell.not_covered == 1]
        #num_prev_visited_og = node.state.steps + len(not_visited_cell_list)
        num_prev_visited = node.state.steps
        
        if len(not_visited_cell_list)==0:
            return 0
        
        cost = 0
        factor = 0.01
        #sorted_list = [cell for cell in not_visited_cell_list]
        #sorted_list.sort(reverse = True, key = lambda p: p.poc)
        
        energy_rem = node.state.energy_remaining
        distance_weight = 0.1164
        
        close_sorted_list = [cell for cell in not_visited_cell_list if math.sqrt((cell.x - node.state.curr_cell.x)**2 + (cell.y - node.state.curr_cell.y)**2) <= energy_rem/distance_weight]
        close_sorted_list.sort(reverse = True, key = lambda p: p.poc)
        
        for cell in close_sorted_list:
            if energy_rem > 0:
                cost += factor*(1-math.exp(-factor*num_prev_visited))*cell.poc
                num_prev_visited += 1
                energy_rem = energy_rem - distance_weight*self.d
            else:
                cost += factor*1*cell.poc
                
        for cell in not_visited_cell_list:
            if cell not in close_sorted_list:
                cost += factor*1*cell.poc
        """
        remaining_not_visited_cell_list = [cell for cell in not_visited_cell_list]
        
        while(1):
            #Given a random cell, calculate the ones adjacent to it
            initial_cell = random.choice(remaining_not_visited_cell_list)
            isolated_group = self.get_isolated_group(remaining_not_visited_cell_list, initial_cell)
            
            
            #Calculate the distance to closest cell outside this group
            if len(isolated_group) < len(remaining_not_visited_cell_list):
                #remove the elements of the isolated group from the not visited list
                aux_list = [cell for cell in not_visited_cell_list if cell not in isolated_group]
                distance_between_groups = math.floor(self.find_closest_cell_dist(isolated_group, aux_list)/self.d)
                #num_prev_visited += distance_between_groups
                for cell in isolated_group:
                    cost += factor*(1-math.exp(-factor*(num_prev_visited_og + sorted_list.index(cell) + distance_between_groups)))*cell.poc
                    
                #Update for the remaining cells to check
                remaining_not_visited_cell_list = [cell for cell in remaining_not_visited_cell_list if cell not in isolated_group]
            else:
                for cell in isolated_group:
                    cost += factor*(1-math.exp(-factor*(num_prev_visited_og + sorted_list.index(cell))))*cell.poc
                break #There are no isolated cells
        """
        
        return cost

    def attraction(self, cell1, cell2):
        """Quantifies the attraction generated by cell2 on cell1"""
        distance = math.sqrt((cell2.x - cell1.x)**2 + (cell2.y - cell1.y)**2)
        #attraction = cell2.poc*cell2.not_covered/math.exp(distance/self.d)
        attraction = cell2.poc*cell2.not_covered/(distance/self.d)**(1/4)
        return attraction

    def calc_movement_dir(self, cell1, valid_cell_list):
        # Finds the cell that produces the biggest attraction
        max_att = 0
        for idx, cell2 in enumerate(valid_cell_list):
            if cell1 != cell2:
                att = self.attraction(cell1,cell2)
                if att> max_att:
                    max_att = att
                    max_cell = cell2
                    max_idx = idx
        #If all the cells have been visited, the max attraction is 0
        if max_att == 0:
            return 0
        # Performs graph search to find the best path to it
        traj_prob = TrajectoryProblem(valid_cell_list, cell1, max_cell)
        solution, path = traj_prob.solve()
        return solution[0]
    
    def get_next_cell(self, cell, valid_cell_list):
        solution = self.calc_movement_dir(cell, valid_cell_list)
        if solution == 0: #all the cells have been visited
            return None, None
        idx, next_cell = search_for_cell_ij(valid_cell_list, cell.i+solution[1],cell.j+solution[2])
        next_cell.prev_action = solution
        return idx, next_cell
    
    def generate_estimated_path(self, state):
        idx, initial_cell = search_for_cell_ij(state.valid_cell_list, state.curr_cell.i,state.curr_cell.j)
        cell = initial_cell
        path = []
        cell.prev_action = state.curr_cell.prev_action
        path.append(cell)
        valid_cell_list = state.valid_cell_list
        valid_cell_list[idx].not_covered = 0
        cont = 0
        #print("Current cell: ", cell.i, cell.j)
        energy_expended = 0
        action_prev = state.curr_cell.prev_action
        while(1):
            idx, next_cell = self.get_next_cell(cell, valid_cell_list)
            action = next_cell.prev_action
            cont = cont + 1
            
            #If considering the energy remaining
            distance_cost = math.sqrt((next_cell.x - state.curr_cell.x)**2 + (next_cell.y - state.curr_cell.y)**2)
            turn_cost = self.calculate_turn(action_prev, action)
            turn_weight = 0.0173
            distance_weight = 0.1164
            energy_expended += distance_weight*distance_cost + turn_weight*turn_cost
            
            if idx == None or energy_expended >= state.energy_remaining: #all the cells have been visited
                break

            action_prev = action
            path.append(copy.deepcopy(next_cell))
            valid_cell_list[idx].not_covered = 0
            cell = next_cell
            #print("Current cell: ", cell.i, cell.j)
        return path, valid_cell_list
    
    #Using a generation of the path using the attraction algorithm to calculate the heuristic
    def h(self, node):
        """ Return the heuristic value for the given state """        
        #Get the list of cells that havent been visited
        not_visited_cell_list = [cell for cell in node.state.valid_cell_list if cell.not_covered == 1]
        if len(not_visited_cell_list)==0 or node.state.energy_remaining <= 0:
            return 0
        num_prev_visited_og = node.state.steps
        
        cost = 0
        #sorted_list = [cell for cell in not_visited_cell_list]
        #sorted_list.sort(reverse = True, key = lambda p: self.attraction(node.state.curr_cell, p))
        estim_path, updated_valid_cell_list = self.generate_estimated_path(copy.deepcopy(node.state))
        
        factor = 0.01
        
        for count, cell in enumerate(estim_path[1:]):
            #info_loss = (1-self.info_penalty(state1.valid_cell_list))*state2.curr_cell.poc
            if cell.not_covered == 1:
                info_loss = factor*(1-math.exp(-factor*(num_prev_visited_og + count + 1)))*cell.poc
                #info_loss = factor*(1-math.exp(-factor*(num_prev_visited_og)))*cell.poc
            else:
                info_loss = 0
            cost = cost + info_loss
        
        for cell in updated_valid_cell_list:
            if cell.not_covered == 1:
                cost += factor*1*cell.poc
        #print(node.path_cost+cost)
        return cost

    def solve(self):
        """ Calls the uninformed search algorithm
        chosen. Returns a solution using the specified format. """
        # The best method found was astar search
        #last_node = search.depth_first_graph_search(self)
        #last_node = search.uniform_cost_search(self)
        last_node = search.astar_search(self)
        #last_node = search.weighted_astar_search(self, epsilon = 0.5)
        #self.test_solution(last_node)
        return last_node.solution(), last_node.path()
    
    def test_solution(self, last_node):
        """ Useful function for debugging - prints the path until the given node """
        path = last_node.path()
        for node in path:
            print("Action: ", node.action)
            print(node.state.disp())
            print("path cost: ", node.path_cost, "h: ", self.h(node), "f(n): ", node.path_cost + self.h(node))
            print('\n')
    
class State2():
    def __init__(self, curr_cell, curr_cell_idx, initial_energy):
        self.curr_cell = curr_cell #Contains the properties of the current cell being occupied
        self.curr_cell_idx = curr_cell_idx #Is the index of the current cell in the valid cell list
        self.energy_remaining = initial_energy #Energy remaining in kJ

    
    def __eq__(self, state):
        """ Returns True if all the attributes in both states are the same """
        
    
        if self.curr_cell.i == state.curr_cell.i and self.curr_cell.j == state.curr_cell.j and self.energy_remaining == state.energy_remaining:
            return True
        else:
            return False
        #return all((cell1.not_covered == cell2.not_covered) for cell1 in self.valid_cell_list for cell2 in state.valid_cell_list if cell1==cell2)
    
    def __hash__(self):
        """ Defines a state object as a hashable object """
        return hash((self.curr_cell_idx, self.energy_remaining))
    
    def __lt__(self, state):
        """ Defines the comparison between state objects."""
        return self.energy_remaining > state.energy_remaining
    
    def display(self):
        print("Current_Cell: ", "i = ", self.curr_cell.i, "j = ", self.curr_cell.j)

class Node:
    """A node in a search tree. Contains a pointer to the parent (the node
    that this is a successor of) and to the actual state for this node. Note
    that if a state is arrived at by two paths, then there are two nodes with
    the same state. Also includes the action that got us to this state, and
    the total path_cost (also known as g) to reach the node. Other functions
    may add an f and h value; see best_first_graph_search and astar_search for
    an explanation of how the f and h values are handled. You will not need to
    subclass this class."""

    def __init__(self, state, parent=None, action=None, path_cost=0, h = 0):
        """Create a search tree Node, derived from a parent by an action."""
        self.state = state
        self.parent = parent
        self.action = action
        self.path_cost = path_cost
        self.h = h
        self.depth = 0
        if parent:
            self.depth = parent.depth + 1

    def __repr__(self):
        return "<Node {}>".format(self.state)

    def __lt__(self, node):
        return self.state < node.state

    #def expand(self, problem):
    #    """List the nodes reachable in one step from this node."""
    #    return [self.child_node(problem, action)
    #            for action in problem.actions(self.state)]

    #def child_node(self, problem, action):
    #    """[Figure 3.10]"""
    #    next_state = problem.result(self.state, action)
    #    next_node = Node(next_state, self, action, problem.path_cost(self.path_cost, self.state, action, next_state))
    #    return next_node

    def solution(self):
        """Return the sequence of actions to go from the root to this node."""
        return [node.action for node in self.path()[1:]]

    def path(self):
        """Return a list of nodes forming the path from the root to this node."""
        node, path_back = self, []
        while node:
            path_back.append(node)
            node = node.parent
        return list(reversed(path_back))

    # We want for a queue of nodes in breadth_first_graph_search or
    # astar_search to have no duplicated states, so we treat nodes
    # with the same state as equal. [Problem: this may not be what you
    # want in other contexts.]

    def __eq__(self, other):
        return isinstance(other, Node) and self.state == other.state

    def __hash__(self):
        # We use the hash value of the state
        # stored in the node instead of the node
        # object itself to quickly search a node
        # with the same state in a Hash Table
        return hash(self.state)

class PriorityQueue:
    """A Queue in which the minimum (or maximum) element (as determined by f and
    order) is returned first.
    If order is 'min', the item with minimum f(x) is
    returned first; if order is 'max', then it is the item with maximum f(x).
    Also supports dict-like lookup."""

    def __init__(self, order='min', f=lambda x: x):
        self.heap = []
        if order == 'min':
            self.f = f
        elif order == 'max':  # now item with max f(x)
            self.f = lambda x: -f(x)  # will be popped first
        else:
            raise ValueError("Order must be either 'min' or 'max'.")

    def append(self, item):
        """Insert item at its correct position."""
        heapq.heappush(self.heap, (self.f(item), item))

    def extend(self, items):
        """Insert each item in items at its correct position."""
        for item in items:
            self.append(item)

    def pop(self):
        """Pop and return the item (with min or max f(x) value)
        depending on the order."""
        if self.heap:
            return heapq.heappop(self.heap)[1]
        else:
            raise Exception('Trying to pop from empty PriorityQueue.')

    def __len__(self):
        """Return current capacity of PriorityQueue."""
        return len(self.heap)

    def __contains__(self, key):
        """Return True if the key is in PriorityQueue."""
        return any([item == key for _, item in self.heap])

    def __getitem__(self, key):
        """Returns the first value associated with key in PriorityQueue.
        Raises KeyError if key is not present."""
        for value, item in self.heap:
            if item == key:
                return value
        raise KeyError(str(key) + " is not in the priority queue")

    def __delitem__(self, key):
        """Delete the first occurrence of key."""
        try:
            del self.heap[[item == key for _, item in self.heap].index(True)]
        except ValueError:
            raise KeyError(str(key) + " is not in the priority queue")
        heapq.heapify(self.heap)


class MCS_Node:
    def __init__(self, cells_idx, uavs_energy_rem, last_uav_idx, parent = None):
        self.last_uav_idx = last_uav_idx
        self.cells_idx = cells_idx #The index of the current cells of the UAVs
        self.uavs_energy_rem = uavs_energy_rem
        self.parent = parent
        self.children = []
        self.Q = 0
        self.N = 0
        self.Q_max = 0
        self.Q_max_max = 0
        self.fully_explored = 0 #All of its child nodes have been fully explored
        self.fully_expanded = 0
    
    def get_nodes_in_path(self):
        aux_node = self
        path = []
        while(aux_node is not None):
            path.append(aux_node)
            aux_node = aux_node.parent
        path.reverse()
        return path
    
    def print(self):
        print("N = ", self.N, " Q = ", self.Q)
        
        
class ACO_Node:
    """Node structure used in the ACO tree."""
    def __init__(self, cells_idx, uavs_energy_rem, last_uav_idx, curr_iter, evaporation, parent = None):
        self.last_uav_idx = last_uav_idx
        self.cells_idx = cells_idx #The index of the current cells of the UAVs
        self.uavs_energy_rem = uavs_energy_rem
        self.parent = parent
        self.children = []
        self.pheromones = (1-evaporation)**curr_iter
        self.n_times_visited_prev = self.number_times_visited_prev()
        self.last_pheromone_update = curr_iter
        self.heuristic = 0
    
    def get_nodes_in_path(self):
        """Returns the node path until the current node."""
        aux_node = self
        path = []
        while(aux_node is not None):
            path.append(aux_node)
            aux_node = aux_node.parent
        path.reverse()
        return path
    
    def number_times_visited_prev(self):
        """Returns the number of times that the cell visited by the UAV that moved in the current node, has been covered in the past by the any UAV."""
        n_times = 0
        aux_node = self
        aux_node = aux_node.parent
        while(aux_node is not None):
            if self.cells_idx[self.last_uav_idx] == aux_node.cells_idx[aux_node.last_uav_idx]:
                n_times += 1
            aux_node = aux_node.parent
        return n_times
    
    def update_and_get_pheromones(self, curr_iter, evaporation, deposit = 0):
        """Returns and updates the value of the pheromones in the current node."""
        self.pheromones = (1-evaporation)**(curr_iter-self.last_pheromone_update) *self.pheromones + deposit
        self.last_pheromone_update = curr_iter
        return self.pheromones