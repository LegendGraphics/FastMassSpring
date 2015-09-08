#include "Solver.h"
#include "Constraint.h"

Solver::Solver()
{

}

Solver::~Solver()
{
}

void Solver::addConstraint(Constraint* constraint)
{
  this->constraints.push_back(constraint);
}

void Solver::initCholesky()
{
  if (this->constraints.empty())
  {
    std::cout << "no constraints now.\n";
    return;
  }
  
  std::cout << "prepare system matrix.\n";
  this->setSystemMatrix();

  // pre factorize the system matrix
  std::cout << "Pre factorize matrix.\n";
  chol.analyzePattern(this->system_matrix);
  chol.factorize(this->system_matrix);
}

void Solver::solve()
{
  for (int i = 0; i < max_iter; ++i)
  {
    this->runLocalStep();
    this->runGlobalStep();
  }
}

void Solver::runOneStep()
{
  this->runLocalStep();
  this->runGlobalStep();
}

void Solver::runGlobalStep()
{
  P_Opt = chol.solve(this->right_hand);
}

void Solver::runLocalStep()
{
  for (decltype(this->constraints.size()) i = 0; i < this->constraints.size(); ++i)
  {
    this->constraints[i]->projection();
    this->constraints[i]->update();
  }

  this->setRightHand();
}

void Solver::setSystemMatrix()
{
  if (this->constraints.empty())
  {
    std::cout << "no constraints now.\n";
    return;
  }

  this->constraints[0]->getLinearSys(this->system_matrix);
  for (decltype(this->constraints.size()) i = 1; i < this->constraints.size(); ++i)
  {
    SparseMatrix temp_matrix;
    this->constraints[i]->getLinearSys(temp_matrix);
    this->system_matrix += temp_matrix;
  }
}

void Solver::setRightHand()
{
  this->right_hand = VectorXf::Zero(problem_size);
  for (decltype(this->constraints.size()) i = 0; i < this->constraints.size(); ++i)
  {
    VectorXf temp_right_hand;
    this->constraints[i]->getRightHand(temp_right_hand);
    this->right_hand += temp_right_hand;
  }
}