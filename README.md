# UAV-path-planning
This repository containts the code for the implementation of the algorithms proposed in "Information-Oriented and Energy-Aware Path Planning for Multiple Air-Dropped Small Unmanned Aerial Vehicles"

## Abstract
Small unmanned aerial vehicles (UAVs) have revolutionized fields such as search and rescue (SAR), surveillance and monitoring. However, their limited endurance presents a challenge in efficiently covering extensive areas during time-sensitive missions. This thesis proposes a novel approach to maximize target detection probability while minimizing search time, by using an innovative formulation of the Coverage Path Planning (CPP) problem combined with a global optimization algorithm. The algorithms presented are information-oriented, leveraging prior knowledge of the target position; and energy-aware, considering the limited endurance of the UAVs. By prioritizing areas with high likelihood of target presence, they increase the probability of target detection and reduce the detection time, thus increasing the chances of survival in SAR missions. This novel approach leverages the simultaneous search capabilities of a multi-UAV system, whose performance can be improved by the use of a larger aircraft capable of rapidly air-dropping the UAVs across large, inaccessible regions. This allows for the optimization of the UAVs' deployment position, contributing to a better mission performance. Four algorithms are proposed to generate the paths of a swarm of UAVs: Attraction (Att), Simulated Annealing (SA), Ant Colony Optimization (ACO) and Monte Carlo Tree Search (MCTS). A parameter optimization of these algorithms is carried out using Bayesian Optimization (BO) and grid search, and their performance is compared for different mission scenarios. Tests showed that SA and MCTS resulted in the best-performing paths, with all of the proposed algorithms out-performing a baseline method which is typically employed in SAR missions.

## Installation
These instructions are for the installation on a Linux system, for Windows or Mac.
The code was developed for Python 3.10.12, and the necessary dependencies can be found in the file requirements.txt. These can be installed using the following command:

```
pip install -r ./requirements.txt
```

This code was developted in order to easily be run in a remote server. If a server is available, fill out the hostname, username and password in line 254 of "main.py". In order to use this functionality, the following library is also required:
```
pip install paramiko
```

## Generating a mission environment
The mission environment is represented by an instance of the class "Mission", defined in "utils.py". In order to change the desired mission environment for the mission, the mission initialization should be edited in the code in "main.py". If a custom mission needs to be created, it should be done according to the structure defined in the class "Mission".

A mission instance created using Mission(0) corresponds to simple mission environment with only 1 NFZ, and with the POC map represented by a single gaussian probability distribution. This map is typical for missions where only one report of the target's last known position is available.

A mission instance created using Mission(1) corresponds to complex mission environment with several NFZs, and with the POC map represented by a weighted sum of 5 different gaussian probability distributions. This map represents a mission where multiple and possibly conflicting target position reports are available.

In order to randomize the mission creation process, a method for generating random AOIs, NFZs and POC maps was implemented. These can be created using Mission(2), and the seed used for the generation can also be passed as a parameter. It should be noted that this randomization process sometimes results in missions where the valid cells are not fully-connected, and therefore, the UAVs might not be able to visit all of them. The following list of mission seeds will result in a set of missions that are representative of real-world SAR missions: [72, 28, 59, 76, 71, 41, 88, 63, 34, 29, 50, 37, 95, 22, 78, 23, 62].

## Generating the paths of the UAVs using the available algorithms
The paths of the swarm of UAVs can be generated using any of the proposed algorithms for a given mission, by running a command in the Linux terminal.

The required structure of command line arguments is the following:
```
python main.py <host> <algorithm> [save result] [file name]
```

The following options are available:

\begin{itemize}
    \item host (required argument): local or remote
    \item algorithm (required argument): UninfAtt, Att, SA, ACO or MCTS
    \item save result (optional argument, default 0): 0 or 1
    \item file name (optional argument, default is date-time): whatever_name_you_want
\end{itemize}

Example command:
```
python main.py local Att 1 john_file_doe
```
