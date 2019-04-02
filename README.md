# LapSim
Vehicle Lap Simulator using three scripts: LapSim, Engine and Track. 
Track:
The track script reads the track file coordinates and creates an x and y arrays of the coordinates in reverse.

Engine:
The Engine script is split into two sections: 
One reads the engine data file coordinates (RPM and Power), the gear ratio file, and takes wheel diameter and final drive ratio.
The next iterates over the gears to determine the power and torque at a given velocity for each gear, including the RPM which the gear
change occurs. A five column array is output. 

LapSim:
The LapSim script takes the engine data, track data and various vehicle values to calculate the time taken to complete the lap.
The LapSim only considers tyre friction, centripetal forces, aerodynamic drag, down forces, mass and efficiency to calculate the lap time.
The calculation iterates for each track point given from the track data. 

To use this, the LapSim script should be run, and will call the other two scripts when required.
