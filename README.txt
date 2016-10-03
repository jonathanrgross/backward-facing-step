INSTRUCTIONS FOR RUNNING JONATHAN GROSS' AME 535a FINAL PROJECT
===============================================================


The code can be run from the command terminal using the following format:

	
	gross_final(Re, nx, ny, IC, tolerance)




Re	  = Reynolds number

nx 	  = number of mesh points in the x direction

ny 	  = number of mesh points in the y direction

IC 	  = select initial conditions

		IC = 1 is uniform flow from inlet across domain

		IC = 2 is parabolic flow from inlet across domain

		IC = 3 is parabolic flow interpolated from inlet to full height of outlet

		IC = 4 is zero velocity at all interior points

tolerance = this allows the user to define a threshold for when to stop the program.  It is

	    based on the sum of the differences between the velocity at two successive time

	    steps.  The program will stop running automatically after 2000 iterations.  
	    
	    Results presented in the report were found with a tolerance of 0.01.





If you wish to copy and paste a command into the command window, you can use the following:



	gross_final(40, 61, 41, 3, 0.01)
