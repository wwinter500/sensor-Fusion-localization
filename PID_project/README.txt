# CarND-Controls-PID
Self-Driving Car Engineer Nanodegree Program

# PID parameters

In this implementation, we use the following parameters as initial:

Kp = 0.12
Ki = 0.001
Kd = 0.7

note: These parameters are chosen from a manual tuning.
------------------------------------------------------

In order to test the effect of the P, I and D components, we introduced a 0.2 bias on the driving angle.

1. With P component only, the car oscillates heavily in the road
2. With PD components, the car can drive on the road, but has a constant offset from the middle of the road
3. With both the PID components, the car now drives much smoothly in the middle of the road