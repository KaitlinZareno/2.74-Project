# 2.74-Project

This Github contains code and images for the Biped walking robot Sllipsis. Sllipsis was our team's final project for MIT's 2.74: Bio-Inspired Robotics class, taught in Fall 2021. 

Motivation:
We wanted to explore the effect of the ankle during bipdal walking; more specifically, we wanted to understand the effect different ankle stiffnesses have on the cost of transport. In addition, we wanted to determine if there would be an optimal walking gait given a specific ankle stiffness. 

Model:
We ended up using a no-slip model defined by 8 generalized coordinates: 2 DoF (x,y translation) for the robotic body, as well as 3 DoF for each leg modelling the hip, knee, and ankle joints. 

Simulation Methods:
- Cartesian Space Control (PD)
- We tracked ellipses at the foot to model different walking gaits for the prior defined optimal ankle stiffness

Simuation Optimization:
- used Fmincon to minimize mechanical cost of transport to determine the optimal gait given a specific ankle stiffness.


