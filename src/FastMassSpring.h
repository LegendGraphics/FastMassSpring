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
*
* This class provides an implementation of fast mass spring system based on
* Tiantian Liu, Adam W. Bargteil, James F. O'Brien, Ladislav Kavan. 
* "Fast Simulation of Mass-Spring Systems" ACM Transaction on Graphics 32(6)
* [Proceedings of SIGGRAPH Asia], 2013.

=========================================================================*/

#ifndef FastMassSpring_H
#define FastMassSpring_H

#include "Constraint.h"
#include "BasicHeader.h" // to make life easier

class Solver;

class FastMassSpring : public Constraint
{
public:
  FastMassSpring();
  virtual ~FastMassSpring();

  void init();
  void initEdgeGraph(
    FACELIST& face_list,
    VERTEXLIST& vertex_list,
    ADJLIST& vertex_share_faces);
  void buildMatrix();
  void fillLMatrix(
    std::vector<Eigen::Triplet<float> >& triplets,
    std::vector<EDGE>& edges, float k);
  void fillJMatrix(
    std::vector<Eigen::Triplet<float> >& triplets,
    std::vector<EDGE>& edges, float k, int edge_counts);
  void computedVector();
  virtual void update();
  virtual void projection();

private:
  bool findShareVertex(
    int& cross_pi, int& cross_pj,
    int pi, int pj,
    ADJLIST& vertex_share_faces, FACELIST& face_list);

  std::vector<EDGE> strech_edges;
  std::vector<EDGE> bending_edges;
  std::vector<float> strech_r_length;
  std::vector<float> bending_r_length;

  Eigen::SparseMatrix<float> L_matrix;
  Eigen::SparseMatrix<float> J_matrix;
  Eigen::VectorXf d_vector;
  Eigen::VectorXf right_hand;

  Solver* solver;

  float k_strech;
  float k_bending;
  int P_Num;

private:
  FastMassSpring(const FastMassSpring&);
  void operator = (const FastMassSpring&);
};

#endif