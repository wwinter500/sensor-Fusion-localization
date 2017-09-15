# CarND-Controls-MPC
Self-Driving Car Engineer Nanodegree Program

-------

## The Model

kinematickinematic model used for the vehicle in this project.
The states variables are: `[x,y,psi,v]` which stand for the position in x and y, the yaw angle, and the speed.
The actuators are `[delta,a]` which stand for the steering angle and the acceleration. 

## Parameters

uses `N=15` and `dt=0.05` as time step parameters. 
The weights for the costs are: 1, 10, 1 for the reference states and  10, 1, for the actuators and 600, 1 for sequential actuation. We use a large weight on the cost relates to the change of delta to ensure a smooth change of steering angle. We use `dt=0.05` so that a the delay of 100ms can be incorporate easily. `N0=15` is to ensure that the trajectory extend long enough into the future but not to long to cost a numerical instability. The target speed is set for 60mph. The parameters tried before are all be 1 expect the change of delta to be 500, which also works.

## Way points

The way points are used to fit to a third order polynomial. The fitting process is in the car coordinates, so that the way points are transformed from the map coordinates to the car coordinates. One thing to notice is that the map coordinates are left-handed while the car coordinates are right-handed so the y values from in the map coordinates are reversed before the coordinate transformation.

## Latency

There is a 100ms delay in the code to simulate the latency of the actuations. In order to account the latency, we use the solution from the solver at a t=0.1s delay. A better way is to transform the vehicle location using the motion model.
