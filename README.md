## On-Line Learning Control of Redundant Robots with Uncertain Dynamics
 
### Master Thesis of Marco Ficorilli 

This repository contains the software that has been developed to implement the method proposed in the Master Thesis "On-Line Learning Control of Redundant Robots with Uncertain Dynamics"

The "software" folder contains the MATLAB code, containing all the methods, the main files that have been used. Also, the results have been stored both as workspaces that as .txt files. The generated figures, instead, have been saved both as .eps and .fig files.

Main files:
- main_KUKA_linGP_sparse.m 
    - Run to simulate the 7R KUKA LWR 4+ robot with a large uncertainty in the link masses;
- main_KUKA_friction.m 
    - Run to simulate the 7R KUKA LWR 4+ robot in a more realistic case, where both 10% mass uncertainty and the presence of complex unmodeled frictions have been taken into account. 

The "videos" folder, instead, contains the videos of the several simulations performed for different trajectory tasks. 
