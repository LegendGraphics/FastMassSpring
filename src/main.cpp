#include "Solver.h"
#include "FastMassSpring.h"
#include "SimCloth.h"

int main(int argc, char *argv[])
{
    FastMassSpring* fms = new FastMassSpring;
    SimCloth* sim = new SimCloth;
    Solver* solver = new Solver;

    FaceList faceList{0,1,2,1,3,2};
    VertexList vertexList{0,1,0,5,1,0,0,1,5,5,1,5};

    solver->addConstraint(fms);
    solver->addConstraint(sim);
    solver->problem_size = vertexList.size();
    solver->P_Opt = Eigen::Map<VectorXf>(&(vertexList)[0], (vertexList).size());

    fms->setSolver(solver);
    fms->setkBending(0);
    fms->initEdgeGraph(faceList, vertexList);
    fms->buildMatrix();

    sim->setSolver(solver);
    sim->setTimeStep(1.0/30.0);
    sim->buildMassMatrix(vertexList.size() / 3, 1.0, vertexList, vertexList);// q_n and q_n_1 changes in every time step

    solver->initCholesky();
    solver->max_iter = 1;
    std::cout << solver->P_Opt.transpose() << std::endl;
    solver->solve();
    std::cout << solver->P_Opt.transpose() << std::endl;

    delete solver;
    delete sim;
    delete fms;
		return 0;
}
