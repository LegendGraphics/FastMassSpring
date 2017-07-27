#ifndef SimCloth_H
#define SimCloth_H

#include "Constraint.h"
#include "BasicHeader.h"

class Solver;

class SimCloth : public Constraint
{
public:
  SimCloth();
  virtual ~SimCloth();

  void setTimeStep(float h) { this->h = h; };
  void buildMassMatrix(int P_Num, float nodeMass, VertexList& P_n, VertexList& P_n_1);
  virtual void init();
  virtual void update();
  virtual void projection();
  virtual void getRightHand(VectorXf& right_hand);
  virtual void getLinearSys(SparseMatrix& linear_sys);
  virtual void setSolver(Solver* solver);

private:
  SparseMatrix mass_matrix;
  VectorXf right_hand;
  VectorXf ext_force;
  VectorXf y;
  float h;

  Solver* solver;
  int P_Num;

private:
  SimCloth(const SimCloth&);
  void operator = (const SimCloth&);
};

#endif