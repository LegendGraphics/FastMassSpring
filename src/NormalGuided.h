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
* Constraint for face normal. This class give a type of constraint about
* desired new face normal to the mesh. With this constraint added into
* the solver, the deformed model will have very close face normal towards
* the given new normal.

=========================================================================*/

#ifndef NormalGuided_H
#define NormalGuided_H

#include "Constraint.h"
#include "BasicHeader.h"

class Solver;

class NormalGuided : public Constraint
{
public:
  NormalGuided();
  virtual ~NormalGuided();

  void initMatrix(
    FaceList& face_list,
    VertexList& vertex_list,
    AdjList& vertex_share_faces,
    NormalList& normal_list,
    NormalList& new_face_normal);
  void fillNormalMatrix(
    FaceList& face_list,
    AdjList& vertex_share_faces,
    NormalList& new_face_normal);
  void fillVMoveMatrix(
    VertexList& vertex_list,
    NormalList& normal_list);
  inline void setLamdNormal(float lamd) { this->lamd_normal = lamd; };
  inline void setLamdVMove(float lamd) { this->lamd_vertical_move = lamd; };

  virtual void init();
  virtual void update();
  virtual void projection();
  virtual void getRightHand(VectorXf& right_hand);
  virtual void getLinearSys(SparseMatrix& linear_sys);
  virtual void setSolver(Solver* solver);

private:
  void getConnectedPtID(int i_pt, int points_in_face[3], int connect_pt[2]);

  float lamd_normal;
  float lamd_vertical_move;
  int P_Num;

  SparseMatrix normal_matrix;
  SparseMatrix vertical_move_matrix;
  VectorXf vertical_move;

private:
  NormalGuided(const NormalGuided&);
  void operator = (const NormalGuided&);
};

#endif