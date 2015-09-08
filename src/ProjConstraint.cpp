#include "ProjConstraint.h"
#include "Solver.h"

ProjConstraint::ProjConstraint()
{
  this->init();
  std::cout << "create a set of projection constraint.\n";
}

ProjConstraint::~ProjConstraint()
{

}

void ProjConstraint::init()
{
  this->lamd_proj = 0.0;
  this->P_Num = 0;
  this->solver = NULL;
}

void ProjConstraint::initMatrix(STLVectorf& ray_list, STLVectori& ray_id_list)
{
  this->P_Num = this->solver->problem_size / 3;
  TripletList proj_triplets;

  for (decltype(ray_id_list.size()) i = 0; i < ray_id_list.size(); ++i)
  {
    int v_id = ray_id_list[i];
    Vector3f ray(ray_list[3 * i + 0], ray_list[3 * i + 1], ray_list[3 * i + 2]);
    Matrix3f differential = Matrix3f::Identity() - ray * ray.transpose();

    for (int k = 0; k < 3; ++k)
    {
      Vector3f diff_temp = differential.col(k).transpose() * differential;
      for (int l = 0; l < 3; ++l)
      {
        proj_triplets.push_back(Triplet(3 * v_id + k, 3 * v_id + l, diff_temp[l]));
      }
    }
  }

  this->constraint_matrix.resize(3 * this->P_Num, 3 * this->P_Num);
  this->constraint_matrix.setFromTriplets(proj_triplets.begin(), proj_triplets.end());
}

void ProjConstraint::update()
{
  // no need to update
}

void ProjConstraint::projection()
{
  // no need to do projection
}

void ProjConstraint::getRightHand(VectorXf& right_hand)
{
  right_hand = VectorXf::Zero(3 * this->P_Num);
}

void ProjConstraint::getLinearSys(SparseMatrix& linear_sys)
{
  linear_sys = this->lamd_proj * this->constraint_matrix;
}

void ProjConstraint::setSolver(Solver* solver)
{
  this->solver = solver;
}