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

  // ray_list stores the ray direction the origin of the
  // ray_list if the camera origin in world coordinate
  void initMatrix(
    STLVectorf& ray_list,
    STLVectori& ray_id_list,
    Vector3f cam_ori);
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
  VectorXf right_hand;

  Solver* solver;

private:
  ProjConstraint(const ProjConstraint&);
  void operator = (const ProjConstraint&);
};

#endif