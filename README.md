# UAV-path-planning
This repository containts the code for the implementation of the algorithms proposed in "Information-Oriented and Energy-Aware Path Planning for Multiple Air-Dropped Small Unmanned Aerial Vehicles"

## Abstract
Small unmanned aerial vehicles (UAVs) have revolutionized fields such as search and rescue (SAR), surveillance and monitoring. However, their limited endurance presents a challenge in efficiently covering extensive areas during time-sensitive missions. This thesis proposes a novel approach to maximize target detection probability while minimizing search time, by using an innovative formulation of the Coverage Path Planning (CPP) problem combined with a global optimization algorithm. The algorithms presented are information-oriented, leveraging prior knowledge of the target position; and energy-aware, considering the limited endurance of the UAVs. By prioritizing areas with high likelihood of target presence, they increase the probability of target detection and reduce the detection time, thus increasing the chances of survival in SAR missions. This novel approach leverages the simultaneous search capabilities of a multi-UAV system, whose performance can be improved by the use of a larger aircraft capable of rapidly air-dropping the UAVs across large, inaccessible regions. This allows for the optimization of the UAVs' deployment position, contributing to a better mission performance. Four algorithms are proposed to generate the paths of a swarm of UAVs: Attraction (Att), Simulated Annealing (SA), Ant Colony Optimization (ACO) and Monte Carlo Tree Search (MCTS). A parameter optimization of these algorithms is carried out using Bayesian Optimization (BO) and grid search, and their performance is compared for different mission scenarios. Tests showed that SA and MCTS resulted in the best-performing paths, with all of the proposed algorithms out-performing a baseline method which is typically employed in SAR missions.

## Installation
This instruction is for the installation on a Linux system, for Windows or Mac.
The code was developed for Python 3.10.12, and the necessary dependencies can be found in the file requirements.txt. These can be installed using the following command:

```
pip install -r ./requirements.txt
```

This code was developted in order to easily be run in a remote server. If a server is available, fill out the hostname, username and password in line ??. In order to use this functionality, the following library is also required:
```
pip install paramiko
```

## Generating a random mission
The mission environment is represented by an instance of the class "Mission", defined in "utils.py".

## Generating the paths of the UAVs using the available algorithms





