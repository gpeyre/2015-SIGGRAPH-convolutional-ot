#ifndef TYPES_H
#define TYPES_H

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>

typedef Eigen::Vector2d Vec2;
typedef Eigen::Vector3d Vec3;
typedef Eigen::ArrayXd  ArrayXd;
typedef Eigen::VectorXd VectorXd;
typedef Eigen::VectorXi VectorXi;

typedef Eigen::MatrixXd Matrix;
typedef Eigen::Triplet<double> Triplet;
typedef Eigen::SparseMatrix<double> SparseMatrix;

#endif // TYPES_H
