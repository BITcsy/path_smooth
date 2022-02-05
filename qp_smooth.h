#pragma once
#ifndef QP_SMOOTH_H
#define QP_SMOOTH_H
// osqp-eigen
#include "OsqpEigen/OsqpEigen.h"
// eigen
#include <Eigen/Dense>
#include "basic_data_def.h"

using namespace std;

class QpSmooth { // convert question to QP, obtain Matrix
public:
	// knot: including init and final state, seg: segment between knots
	void SetMatSize(unit knotNum, unit segNum); 
	bool GenerateHessian(double dIn);
	bool GenerateGradiant();
	bool GenerateLinearMat(const vector<WayPoint>& WayPts, double dIn);
	bool GenerateBounds(const vector<WayPoint>& WayPts);
	bool SetQpSolver(OsqpEigen::Solver& solver);
	void DrawSmoothPath(const Eigen::VectorXd& QPSolution, const vector<WayPoint>& sceWayPt, const double reso,
		vector<double>& lPath);
private:
	Eigen::SparseMatrix<double> hessian_;
	Eigen::VectorXd gradient_;
	Eigen::SparseMatrix<double> linearMatrix_;
	Eigen::VectorXd lowerBound_;
	Eigen::VectorXd upperBound_;
	unit knotNum_;
	unit segNum_;
	unit numOfVar_;
	unit numOfCons_;
	double d_ = 0.5 * 5; // distance between knots
	double d2_;
	double d3_;
	double d4_;
	double d5_;
	double d6_;
	double d7_;
	double d8_;
	double d9_;
	double w1_ = 1.0;
	double w2_ = 1.0;
	double ds_;
};

#endif // !QP_SMOOTH_H