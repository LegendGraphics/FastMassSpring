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

void ProjConstraint::initMatrix(
  STLVectorf& ray_list,
  STLVectori& ray_id_list,
  Vector3f cam_ori)
{
  // actually this term is similar to the vertical move term in
  // normal guided constraints
  this->P_Num = this->solver->problem_size / 3;
  TripletList proj_triplets;
  this->right_hand = VectorXf::Zero(3 * this->P_Num);

  for (decltype(ray_id_list.size()) i = 0; i < ray_id_list.size(); ++i)
  {
    //int v_id = ray_id_list[i];
    //Vector3f ray(ray_list[3 * i + 0], ray_list[3 * i + 1], ray_list[3 * i + 2]);
    //Matrix3f differential = Matrix3f::Identity()
    //                      - (ray * ray.transpose())
    //                      / (ray.transpose() * ray);

    //for (int k = 0; k < 3; ++k)
    //{
    //  Vector3f diff_temp = differential.col(k).transpose() * differential;
    //  for (int l = 0; l < 3; ++l)
    //  {
    //    proj_triplets.push_back(Triplet(3 * v_id + k, 3 * v_id + l, diff_temp[l]));
    //  }
    //}

    Vector3f ray(ray_list[3 * i + 0], ray_list[3 * i + 1], ray_list[3 * i + 2]);
    ray.normalize();

    proj_triplets.push_back(Triplet(3 * i + 0, 3 * i + 0, 1 - ray[0] * ray[0]));
    proj_triplets.push_back(Triplet(3 * i + 0, 3 * i + 1, -ray[0] * ray[1]));
    proj_triplets.push_back(Triplet(3 * i + 0, 3 * i + 2, -ray[0] * ray[2]));
    right_hand(3 * i + 0) = (1 - ray[0] * ray[0]) * cam_ori[0]
                          - ray[0] * ray[1] * cam_ori[1]
                          - ray[0] * ray[2] * cam_ori[2];

    proj_triplets.push_back(Triplet(3 * i + 1, 3 * i + 0, -ray[1] * ray[0]));
    proj_triplets.push_back(Triplet(3 * i + 1, 3 * i + 1, 1 - ray[1] * ray[1]));
    proj_triplets.push_back(Triplet(3 * i + 1, 3 * i + 2, -ray[1] * ray[2]));
    right_hand(3 * i + 1) = -ray[1] * ray[0] * cam_ori[0]
                          + (1 - ray[1] * ray[1]) * cam_ori[1]
                          - ray[1] * ray[2] * cam_ori[2];

    proj_triplets.push_back(Triplet(3 * i + 2, 3 * i + 0, -ray[2] * ray[0]));
    proj_triplets.push_back(Triplet(3 * i + 2, 3 * i + 1, -ray[2] * ray[1]));
    proj_triplets.push_back(Triplet(3 * i + 2, 3 * i + 2, 1 - ray[2] * ray[2]));
    right_hand(3 * i + 2) = -ray[2] * ray[0] * cam_ori[0]
                          - ray[2] * ray[1] * cam_ori[1]
                          + (1 - ray[2] * ray[2]) * cam_ori[2];
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