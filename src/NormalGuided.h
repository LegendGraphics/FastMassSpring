#ifndef NormalGuided_H
#define NormalGuided_H

#include "Constraint.h"
#include "BasicHeader.h"

class Solver;

class NormalGuided : public Constraint
{
public:
  NormalGuided();
  ~NormalGuided();

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

  virtual void update();
  virtual void projection();
  virtual void getRightHand(VectorX& right_hand);
  virtual void getLinearSys(SparseMatrix& linear_sys);
  virtual void setSolver(Solver* solver);

private:
  void getConnectedPtID(int i_pt, int points_in_face[3], int connect_pt[2]);

  float lamd_normal;
  float lamd_vertical_move;
  int P_Num;

  SparseMatrix normal_matrix;
  SparseMatrix vertical_move_matrix;
  VectorX vertical_move;

private:
  NormalGuided(const NormalGuided&);
  void operator = (const NormalGuided&);
};

#endif