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
    FaceList& face_list,
    VertexList& vertex_list,
    AdjList& vertex_share_faces);
  void buildMatrix();
  void fillLMatrix(TripletList& triplets, Edges& edges);
  void fillJMatrix(TripletList& triplets, Edges& edges, int edge_counts);
  void computedVector();
  inline void setkStrech(float k) { this->k_strech = k; };
  inline void setkBending(float k) { this->k_bending = k; };


  virtual void update();
  virtual void projection();
  virtual void getRightHand(VectorX& right_hand);
  virtual void getLinearSys(SparseMatrix& linear_sys);
  virtual void setSolver(Solver* solver);

private:
  bool findShareVertex(
    int& cross_pi, int& cross_pj,
    int pi, int pj,
    AdjList& vertex_share_faces, FaceList& face_list);

  Edges strech_edges;
  Edges bending_edges;
  std::vector<float> strech_r_length;
  std::vector<float> bending_r_length;

  SparseMatrix L_strech_matrix;
  SparseMatrix L_bending_matrix;
  SparseMatrix J_strech_matrix;
  SparseMatrix J_bending_matrix;
  VectorX d_vector;
  VectorX right_hand;

  Solver* solver;

  float k_strech;
  float k_bending;
  int P_Num;

private:
  FastMassSpring(const FastMassSpring&);
  void operator = (const FastMassSpring&);
};

#endif