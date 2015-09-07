#include "FastMassSpring.h"
#include "Solver.h"

#include <iostream>

FastMassSpring::FastMassSpring()
{
  this->init();
  std::cout << "create a fast mass spring sys.\n";
}

FastMassSpring::~FastMassSpring()
{

}

void FastMassSpring::init()
{
  // init necessary variable here
  this->k_strech = 1.0;
  this->k_bending = 1.0;
  this->P_Num = 0;
  this->solver = NULL;
}

void FastMassSpring::initEdgeGraph(
  FaceList& face_list,
  VertexList& vertex_list,
  AdjList& vertex_share_faces)
{
  this->P_Num = vertex_list.size() / 3;

  // building edge graph
  this->strech_edges.clear();
  {
    int ptid[3] = {0, 0, 0};
    for (decltype(face_list.size()) i = 0; i < face_list.size() / 3; ++i)
    {
      ptid[0] = face_list[3 * i + 0];
      ptid[1] = face_list[3 * i + 1];
      ptid[2] = face_list[3 * i + 2];
      // the order of start point and end point doesn't matter
      this->strech_edges.push_back(
        ptid[0] < ptid[1] ? Edge(ptid[0], ptid[1]) : Edge(ptid[1], ptid[0]));
      this->strech_edges.push_back(
        ptid[1] < ptid[2] ? Edge(ptid[1], ptid[2]) : Edge(ptid[2], ptid[1]));
      this->strech_edges.push_back(
        ptid[2] < ptid[0] ? Edge(ptid[2], ptid[0]) : Edge(ptid[0], ptid[2]));
    }
    std::sort(this->strech_edges.begin(), this->strech_edges.end());
    std::vector<Edge>::iterator iter = 
      std::unique(this->strech_edges.begin(), this->strech_edges.end());
    this->strech_edges.erase(iter, strech_edges.end());
    this->strech_edges.shrink_to_fit();

    // store rest length
    float cur_r = 0.0;
    this->strech_r_length.clear();
    for (auto& i : this->strech_edges)
    {
      cur_r = 0.0;
      for (int j = 0; j < 3; ++j)
      {
        cur_r += pow(vertex_list[3 * i.first + j] - vertex_list[ 3 * i.second + j], 2);
      }
      this->strech_r_length.push_back(sqrt(cur_r));
    }
  }

  this->bending_edges.clear();
  {
    for (auto& i : this->strech_edges)
    {
      int cross_pi = -1;
      int cross_pj = -1;
      if (findShareVertex(
        cross_pi, cross_pj, i.first, i.second, vertex_share_faces, face_list))
      {
        this->bending_edges.push_back(
          cross_pi < cross_pj ? 
            Edge(cross_pi, cross_pj) : Edge(cross_pj, cross_pi));
      }
    }

    // store rest length
    float cur_r = 0.0;
    this->bending_r_length.clear();
    for (auto& i : this->bending_edges)
    {
      cur_r = 0.0;
      for (int j = 0; j < 3; ++j)
      {
        cur_r += pow(vertex_list[3 * i.first + j] - vertex_list[ 3 * i.second + j], 2);
      }
      this->bending_r_length.push_back(sqrt(cur_r));
    }
  }

  // init d vector
  this->computedVector();
}

void FastMassSpring::buildMatrix()
{
  TripletList L_triplets;
  this->L_strech_matrix.resize(3 * this->P_Num, 3 * this->P_Num);
  this->fillLMatrix(L_triplets, this->strech_edges);
  this->L_strech_matrix.setFromTriplets(L_triplets.begin(), L_triplets.end());

  L_triplets.clear();
  this->L_bending_matrix.resize(3 * this->P_Num, 3 * this->P_Num);
  this->fillLMatrix(L_triplets, this->bending_edges);
  this->L_bending_matrix.setFromTriplets(L_triplets.begin(), L_triplets.end());

  TripletList J_triplets;
  this->J_strech_matrix.resize(
    3 * this->P_Num,
    3 * (this->strech_edges.size() + this->bending_edges.size()));
  this->fillJMatrix(J_triplets, this->strech_edges, 0);
  this->J_strech_matrix.setFromTriplets(J_triplets.begin(), J_triplets.end());

  J_triplets.clear();
  this->J_bending_matrix.resize(
    3 * this->P_Num,
    3 * (this->strech_edges.size() + this->bending_edges.size()));
  this->fillJMatrix(J_triplets, this->bending_edges, this->strech_edges.size());
  this->J_bending_matrix.setFromTriplets(J_triplets.begin(), J_triplets.end());
}

void FastMassSpring::fillLMatrix(TripletList& triplets, Edges& edges)
{
  for (auto& i : edges)
  {
    for (int j = 0; j < 3; ++j)
    {
      triplets.push_back(Triplet(3 * i.first + j, 3 * i.first + j, 1));
      triplets.push_back(Triplet(3 * i.first + j, 3 * i.second + j, -1));
      triplets.push_back(Triplet(3 * i.second + j, 3 * i.first + j, -1));
      triplets.push_back(Triplet(3 * i.second + j, 3 * i.second + j, 1));
    }
  }
}

void FastMassSpring::fillJMatrix(TripletList& triplets, Edges& edges, int edge_counts)
{
  for (decltype(edges.size()) i = 0; i != edges.size(); ++i)
  {
    for (int j = 0; j < 3; ++j)
    {
      triplets.push_back(Triplet(3 * edges[i].first + j, 3 * (i + edge_counts) + j, 1));
      triplets.push_back(Triplet(3 * edges[i].second + j, 3 * (i + edge_counts) + j, -1));
    }
  }
}

void FastMassSpring::computedVector()
{
  this->d_vector = VectorX::Zero(
    3 * (this->strech_edges.size() + this->bending_edges.size()));

  size_t id_edge = 0;
  for (auto& i : this->strech_edges)
  {
    Vector3 p12;
    for (int j = 0; j < 3; ++j)
    {
      p12[j] =
        this->solver->P_Opt[3 * i.first + j]
      - this->solver->P_Opt[3 * i.second + j];
    }
    p12.normalize();
    this->d_vector[3 * id_edge + 0] = this->strech_r_length[id_edge] * p12[0];
    this->d_vector[3 * id_edge + 1] = this->strech_r_length[id_edge] * p12[1];
    this->d_vector[3 * id_edge + 2] = this->strech_r_length[id_edge] * p12[2];
    ++id_edge;
  }
  size_t id_edge_offset = 0;
  for (auto& i : this->bending_edges)
  {
    Vector3 p12;
    for (int j = 0; j < 3; ++j)
    {
      p12[j] =
        this->solver->P_Opt[3 * i.first + j]
      - this->solver->P_Opt[3 * i.second + j];
    }
    p12.normalize();
    this->d_vector[3 * id_edge + 0] = this->bending_r_length[id_edge_offset] * p12[0];
    this->d_vector[3 * id_edge + 1] = this->bending_r_length[id_edge_offset] * p12[1];
    this->d_vector[3 * id_edge + 2] = this->bending_r_length[id_edge_offset] * p12[2];
    ++id_edge;
    ++id_edge_offset;
  }
}

bool FastMassSpring::findShareVertex(
  int& cross_pi, int& cross_pj,
  int pi, int pj,
  AdjList& vertex_share_faces, FaceList& face_list)
{
  std::vector<int> share_face;
  std::set_intersection(
    vertex_share_faces[pi].begin(), vertex_share_faces[pi].end(),
    vertex_share_faces[pj].begin(), vertex_share_faces[pj].end(),
    std::back_inserter(share_face));

  if (share_face.size() == 2)
  {
    // set edge points
    std::vector<int> edgePoints;
    edgePoints.clear();
    edgePoints.push_back(pi);
    edgePoints.push_back(pj);

    // set all points in the two faces
    int f0_0 = face_list[3 * share_face[0] + 0];
    int f0_1 = face_list[3 * share_face[0] + 1];
    int f0_2 = face_list[3 * share_face[0] + 2];

    int f1_0 = face_list[3 * share_face[1] + 0];
    int f1_1 = face_list[3 * share_face[1] + 1];
    int f1_2 = face_list[3 * share_face[1] + 2];

    if (f0_0 != pi && f0_0 != pj) cross_pi = f0_0;
    else if (f0_1 != pi && f0_1 != pj) cross_pi = f0_1;
    else if (f0_2 != pi && f0_2 != pj) cross_pi = f0_2;

    if (f1_0 != pi && f1_0 != pj) cross_pj = f1_0;
    else if (f1_1 != pi && f1_1 != pj) cross_pj = f1_1;
    else if (f1_2 != pi && f1_2 != pj) cross_pj = f1_2;

    return true;
  }
  else
  {
    //std::cout << "Error: number of shared face is more than 2...\n";
    return false;
  }
}

void FastMassSpring::projection()
{
  this->computedVector();
}

void FastMassSpring::update()
{
  this->right_hand =
    (this->k_strech * this->J_strech_matrix
   + this->k_bending * this->J_bending_matrix)
   * this->d_vector;
}

void FastMassSpring::getRightHand(VectorX& right_hand)
{
  right_hand = this->right_hand;
}

void FastMassSpring::getLinearSys(SparseMatrix& linear_sys)
{
  linear_sys = this->k_strech * this->L_strech_matrix
             + this->k_bending * this->L_bending_matrix;
}

void FastMassSpring::setSolver(Solver* solver)
{
  this->solver = solver;
}