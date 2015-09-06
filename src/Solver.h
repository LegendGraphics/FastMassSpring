/*=========================================================================

* @file
* @author  Lin Ma <majcjc@gmail.com>
* @version 1.0
*
* @section LICENSE
*
* This program is free software; you can redistribute it and/or
* modify it under the terms of the GNU General Public License as
* published by the Free Software Foundation; either version 2 of
* the License, or (at your option) any later version.
*
* This program is distributed in the hope that it will be useful, but
* WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
* General Public License for more details at
* http://www.gnu.org/copyleft/gpl.html
*
* @section DESCRIPTION
*
* The ultimate goal of these files is to build a framework for a light-weight
* version of "projective dynamics"
*
* Projective dynamics is a method for nonlinear optimization. It separates
* linear part and non-linear part of the problem into global step and local
* step. It first look for a projection for all constraints (local step, solve
* many small non-linear subproblems). Then it combines all the projections
* together and build a large sparse linear system to move one step. These two
* steps are repeated until convergence.

=========================================================================*/

#ifndef Solver_H
#define Solver_H

#include "BasicHeader.h"

class Solver
{
public:
  Solver();
  ~Solver();

  Eigen::VectorXf P_Opt;

private:
  Solver(const Solver&); // not implemented
  void operator = (const Solver&); // not implemented
};

#endif