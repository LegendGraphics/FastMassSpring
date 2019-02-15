# Fast Mass Spring System
This is a modified c++ implementation of fast mass spring system. The main goal of this project is to build a framework for a light-weight version of "projective dynamics"

Projective dynamics is a method for nonlinear optimization. It separates linear part and non-linear part of the problem into global step and local step. It first look for a projection for all constraints (local step, solve many small non-linear subproblems). Then it combines all the projections together and build a large sparse linear system to move one step. These two steps are repeated until convergence.


[1] Tiantian Liu, Adam W. Bargteil, James F. O'Brien, Ladislav Kavan. "Fast Simulation of Mass-Spring Systems" ACM Transaction on Graphics 32(6) [Proceedings of SIGGRAPH Asia], 2013.

[2] Sofien Bouaziz, Mario Deuss, Yuliy Schwartzburg, Thibaut Weise, Mark Pauly "Shape-Up: Shaping Discrete Geometry with Projections" [Computer Graphics Forum (Proceedings of SGP)], 2012.
