#ifndef ProjConstraint_H
#define ProjConstraint_H

#include "Constraint.h"
#include "BasicHeader.h"

class Solver;

class ProjConstraint : public Constraint
{
public:
  ProjConstraint();
  virtual ~ProjConstraint();

  void initMatrix(
    STLVectorf& ray_list,
    STLVectori& ray_id_list);
  inline void setLamdProj(float lamd) { this->lamd_proj = lamd; };

  virtual void init();
  virtual void update();
  virtual void projection();
  virtual void getRightHand(VectorXf& right_hand);
  virtual void getLinearSys(SparseMatrix& linear_sys);
  virtual void setSolver(Solver* solver);

private:
  float lamd_proj;
  int P_Num;

  SparseMatrix constraint_matrix;

  Solver* solver;

private:
  ProjConstraint(const ProjConstraint&);
  void operator = (const ProjConstraint&);
};

#endif