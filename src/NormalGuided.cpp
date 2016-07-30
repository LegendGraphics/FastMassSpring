#include "NormalGuided.h"
#include "Solver.h"

NormalGuided::NormalGuided()
{
  this->init();
  std::cout << "create a set of normal constraint.\n";
}

NormalGuided::~NormalGuided()
{

}

void NormalGuided::init()
{
  this->lamd_normal = 0.0;
  this->lamd_vertical_move = 0.0;
  this->P_Num = 0;
}

void NormalGuided::initMatrix(
  FaceList& face_list,
  VertexList& vertex_list,
  AdjList& vertex_share_faces,
  NormalList& normal_list,
  NormalList& new_face_normal)
{
  // normal_list is face normal
  this->P_Num = vertex_list.size() / 3;

  this->fillNormalMatrix(face_list, vertex_share_faces, new_face_normal);
  this->fillVMoveMatrix(vertex_list, normal_list);
}

void NormalGuided::fillNormalMatrix(
  FaceList& face_list,
  AdjList& vertex_share_faces,
  NormalList& new_face_normal)
{
  TripletList normal_triplets;

  for (int i = 0; i < P_Num; ++i)
  {
    std::vector<int>& one_ring_faces = vertex_share_faces[i];

    for (decltype(one_ring_faces.size()) j = 0; j < one_ring_faces.size(); ++j)
    {
      float n[3] = { new_face_normal[3 * one_ring_faces[j] + 0],
                     new_face_normal[3 * one_ring_faces[j] + 1],
                     new_face_normal[3 * one_ring_faces[j] + 2] };

      int points_in_face[3] = { face_list[3 * one_ring_faces[j] + 0],
        face_list[3 * one_ring_faces[j] + 1],
        face_list[3 * one_ring_faces[j] + 2] };
      int connect[2];
      this->getConnectedPtID(i, points_in_face, connect);

      for (int k = 0; k < 3; ++k)
      {
        for (int l = 0; l < 3; ++l)
        {
          // 2 edges first dimension
          normal_triplets.push_back(Triplet(3 * i + k, 3 * i + l, n[k] * n[l] * 2));
          normal_triplets.push_back(Triplet(3 * i + k, 3 * connect[0] + l, -n[k] * n[l]));
          normal_triplets.push_back(Triplet(3 * i + k, 3 * connect[1] + l, -n[k] * n[l]));
        }
      }
    }
  }

  normal_matrix.resize(3 * P_Num, 3 * P_Num);
  normal_matrix.setFromTriplets(normal_triplets.begin(), normal_triplets.end());
}

void NormalGuided::fillVMoveMatrix(
  VertexList& vertex_list,
  NormalList& normal_list)
{
  TripletList vertical_move_triplets;
  vertical_move = VectorXf::Zero(3 * P_Num);

  for (int i = 0; i < P_Num; ++i)
  {
    float n[3] = { normal_list[3 * i + 0],
      normal_list[3 * i + 1],
      normal_list[3 * i + 2] };

    float pt[3] = { vertex_list[3 * i + 0],
      vertex_list[3 * i + 1],
      vertex_list[3 * i + 2] };

    vertical_move_triplets.push_back(Triplet(3 * i + 0, 3 * i + 0, 1 - n[0] * n[0]));
    vertical_move_triplets.push_back(Triplet(3 * i + 0, 3 * i + 1, -n[0] * n[1]));
    vertical_move_triplets.push_back(Triplet(3 * i + 0, 3 * i + 2, -n[0] * n[2]));
    vertical_move(3 * i + 0) = (1 - n[0] * n[0]) * pt[0] - n[0] * n[1] * pt[1] - n[0] * n[2] * pt[2];

    vertical_move_triplets.push_back(Triplet(3 * i + 1, 3 * i + 0, -n[1] * n[0]));
    vertical_move_triplets.push_back(Triplet(3 * i + 1, 3 * i + 1, 1 - n[1] * n[1]));
    vertical_move_triplets.push_back(Triplet(3 * i + 1, 3 * i + 2, -n[1] * n[2]));
    vertical_move(3 * i + 1) = -n[1] * n[0] * pt[0] + (1 - n[1] * n[1]) * pt[1] - n[1] * n[2] * pt[2];

    vertical_move_triplets.push_back(Triplet(3 * i + 2, 3 * i + 0, -n[2] * n[0]));
    vertical_move_triplets.push_back(Triplet(3 * i + 2, 3 * i + 1, -n[2] * n[1]));
    vertical_move_triplets.push_back(Triplet(3 * i + 2, 3 * i + 2, 1 - n[2] * n[2]));
    vertical_move(3 * i + 2) = -n[2] * n[0] * pt[0] - n[2] * n[1] * pt[1] + (1 - n[2] * n[2]) * pt[2];
  }

  vertical_move_matrix.resize(3 * P_Num, 3 * P_Num);
  vertical_move_matrix.setFromTriplets(vertical_move_triplets.begin(), vertical_move_triplets.end());
}

void NormalGuided::getConnectedPtID(int i_pt, int points_in_face[3], int connect_pt[2])
{
  if (points_in_face[0] == i_pt) {
    connect_pt[0] = points_in_face[1];
    connect_pt[1] = points_in_face[2];
  }
  else if (points_in_face[1] == i_pt) {
    connect_pt[0] = points_in_face[0];
    connect_pt[1] = points_in_face[2];
  }
  else if (points_in_face[2] == i_pt)
  {
    connect_pt[0] = points_in_face[0];
    connect_pt[1] = points_in_face[1];
  }
  else
  {
    //std::cout << "Error: can not find point in the given face!\n";
    // get the other point's id in this cell
  }
}

void NormalGuided::projection()
{
  // no need to do projection
}

void NormalGuided::update()
{
  // no need to update
}

void NormalGuided::getRightHand(VectorXf& right_hand)
{
  right_hand = this->lamd_vertical_move * this->vertical_move;
}

void NormalGuided::getLinearSys(SparseMatrix& linear_sys)
{
  linear_sys = this->lamd_normal * this->normal_matrix
             + this->lamd_vertical_move * this->vertical_move_matrix;
}

void NormalGuided::setSolver(Solver* solver)
{
  // do nothing, this constraints doesn't update through solver
}