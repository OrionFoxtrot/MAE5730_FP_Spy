# Cornell MAE 5730-Intermediate Dynamics Final Project 

## Shortform setup
My project is to model the dynamics of a human taking a zipline through a 3 section piecewise continuous 2D course. The human, or “spy” in the context of this document consists of a 3-link body connected to a point mass. The point mass is representative of the zipline trolley, and the 3 links represent the head-neck connection, neck-waist connection and then waist to toe connection. This can be altered for an alternative representation where the connection centers of gravity go Power Sledge -> Arms -> Torso -> Legs. 

## Run Instructions
Run the second code block in each of the EL_Approach.m or DAE_Approach.m sections to generate the equations of motions. They are then stored in /EL_Eqs and /DAE_Matrix. These equations are used in the third codeblock ("Simulation") of each of the respective .m files. The first code block can be used to regenerate different course sections individually.

The output of the simulation codeblock is a matrix "X" which stores the state variables. These can be used for plotting by running the "EL_Plotting_CG_and_Linkage.m" and "DAE_Plotting_CG_and_Linkage.m" files which simulate the motion.

Trajectory plots may also be called via writing Trajj(X,n) where X is the state variable matrix and n is this method. Input "1" for the EL approach and "2" for the DAE Approach. For example, after simulating the DAE approach and wanting to inspect the trajectory, write Trajj(X,1) into the command line.

For more information such as detailed setup, constraints, derivations and equations of motion, check out /Report/Report.pdf.

