#include "SimCloth.h"
#include "Solver.h"

SimCloth::SimCloth()
{
    init();
}

SimCloth::~SimCloth()
{

}

void SimCloth::buildMassMatrix(int P_Num, float nodeMass, VertexList& P_n, VertexList& P_n_1)
{
    this->P_Num = P_Num;

    TripletList massTriplets;
    this->mass_matrix.resize(3 * this->P_Num, 3 * this->P_Num);
    for (int i = 0; i < 3 * this->P_Num; ++i)
    {
        massTriplets.push_back(Triplet(i, i, nodeMass));
    }
    this->mass_matrix.setFromTriplets(massTriplets.begin(), massTriplets.end());

    ext_force = VectorXf::Zero(3 * this->P_Num);
    y = 2*Eigen::Map<VectorXf>(&(P_n)[0], (P_n).size()) - Eigen::Map<VectorXf>(&(P_n_1)[0], (P_n_1).size());
}

void SimCloth::init()
{
    this->solver = NULL;
    this->P_Num = 0;
}

void SimCloth::update()
{
    this->right_hand = 1 * (mass_matrix * y - ext_force) / (h*h);
}

void SimCloth::projection()
{
    // do nothing
}

void SimCloth::getRightHand(VectorXf& right_hand)
{
    right_hand = this->right_hand;
}

void SimCloth::getLinearSys(SparseMatrix& linear_sys)
{
    linear_sys = this->mass_matrix / (h*h);
}

void SimCloth::setSolver(Solver* solver)
{
    this->solver = solver;
}