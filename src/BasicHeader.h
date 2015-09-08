#ifndef BasicHeader_H
#define BasicHeader_H

#include <Eigen\Eigen>
#include <Eigen\Sparse>
#include <vector>
#include <iostream>

typedef std::pair<int, int> Edge;
typedef std::vector<std::pair<int, int> > Edges;
typedef std::vector<unsigned int> FaceList;
typedef std::vector<float> VertexList;
typedef std::vector<float> NormalList;
typedef std::vector<float> STLVectorf;
typedef std::vector<int> STLVectori;
typedef std::vector<std::vector<int> > AdjList;

typedef Eigen::SparseMatrix<float> SparseMatrix;
typedef Eigen::Triplet<float> Triplet;
typedef std::vector<Eigen::Triplet<float> > TripletList;
typedef Eigen::VectorXf VectorXf;
typedef Eigen::Vector3f Vector3f;
typedef Eigen::Matrix3f Matrix3f;
typedef Eigen::MatrixXf MatrixXf;
typedef Eigen::SimplicialCholesky<Eigen::SparseMatrix<float> > SimplicialCholesky;

#endif