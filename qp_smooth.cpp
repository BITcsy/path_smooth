#include <iostream>
#include "matplotlibcpp.h"
#include "qp_smooth.h"
#include "scenarios.h"

#define PI 3.1415926

namespace plt = matplotlibcpp;
using namespace std;

void QpSmooth::SetMatSize(unit knotNum, unit segNum)
{
    knotNum_ = knotNum;
    segNum_ = segNum;
    // 5th polynominal between knots
    numOfVar_ = (knotNum - 1) * 6; 
    // 2rd continuous is considered. knotNum_ * 3 equal constrain. ((knotNum_ - 1) * segNum_ - 1) boundary constrain.
    numOfCons_ = knotNum_ * 3 + ((knotNum_ - 1) * segNum_ - 1);
    
    // Mat size
    hessian_.resize(numOfVar_, numOfVar_);
    gradient_.resize(numOfVar_);
    linearMatrix_.resize(numOfCons_, numOfVar_);
    lowerBound_.resize(numOfCons_);
    upperBound_.resize(numOfCons_);
    // debug
    cout << "----Set Mat Size-----" << endl;
    cout << "numOfVar" << numOfVar_ << endl;
    cout << "numOfCons" << numOfCons_ << endl;
    return;
}

bool QpSmooth::GenerateHessian(double dIn)
{
    cout << "----Generate Hessian-----" << endl;
    d_ = dIn;
    d2_ = d_ * d_;
    d3_ = d2_ * d_;
    d4_ = d3_ * d_;
    d5_ = d4_ * d_;
    d6_ = d5_ * d_;
    d7_ = d6_ * d_;
    d8_ = d7_ * d_;
    d9_ = d8_ * d_;

    /*double hessian11 = w1_ * d_ * 2.0;
    double hessian12 = w1_ * d2_ * 2.0;
    double hessian13 = w1_ * d3_ * 2.0;
    double hessian14 = w1_ * d4_ * 2.0;
    double hessian15 = w1_ * d5_ * 2.0;
    double hessian21 = w1_ * d2_ * 2.0;
    double hessian22 = (w1_ * 4.0 / 3.0 * d3_ + w2_ * 4.0 * d_) * 2.0;
    double hessian23 = (w1_ * 6.0 / 4.0 * d4_ + w2_ * 6.0 * d2_) * 2.0;
    double hessian24 = (w1_ * 8.0 / 5.0 * d5_ + w2_ * 8.0 * d3_) * 2.0;
    double hessian25 = (w1_ * 10.0 / 6.0 * d6_ + w2_ * 10 * d4_) * 2.0;
    double hessian31 = w1_ * d3_ * 2.0;
    double hessian32 = (w1_ * 6.0 / 4.0 * d4_ + w2_ * 6.0 * d2_) * 2.0;
    double hessian33 = (w1_ * 9.0 / 5.0 * d5_ + w2_ * 12.0 * d3_) * 2.0;
    double hessian34 = (w1_ * 2.0 * d6_ + w2_ * 18.0 * d4_) * 2.0;
    double hessian35 = (w1_ * 15.0 / 7.0 * d7_ + w2_ * 24.0 * d5_) * 2.0;
    double hessian41 = w1_ * d4_ * 2.0;
    double hessian42 = (w1_ * 8.0 / 5.0 * d5_ + w2_ * 8.0 * d3_) * 2.0;
    double hessian43 = (w1_ * 2.0 * d6_ + w2_ * 18.0 * d4_) * 2.0;
    double hessian44 = (w1_ * 16.0 / 7.0 * d7_ + w2_ * 144.0 / 5.0 * d5_) * 2.0;
    double hessian45 = (w1_ * 20.0 / 8.0 * d8_ + w2_ * 40.0 * d6_) * 2.0;
    double hessian51 = w1_ * d5_ * 2.0;
    double hessian52 = (w1_ * 10.0 / 6.0 * d6_ + w2_ * 10 * d4_) * 2.0;
    double hessian53 = (w1_ * 15.0 / 7.0 * d7_ + w2_ * 24.0 * d5_) * 2.0;
    double hessian54 = (w1_ * 20.0 / 8.0 * d8_ + w2_ * 40.0 * d6_) * 2.0;
    double hessian55 = (w1_ * 25.0 / 9.0 * d9_ + w2_ + 400.0 / 7.0 * d7_) * 2.0;*/

    for (uint i = 0; i < knotNum_ - 1; i++) {
        hessian_.insert(6 * i + 1, 6 * i + 1) = w1_ * d_ * 2.0;
        hessian_.insert(6 * i + 1, 6 * i + 2) = w1_ * d2_ * 2.0;
        hessian_.insert(6 * i + 1, 6 * i + 3) = w1_ * d3_ * 2.0;
        hessian_.insert(6 * i + 1, 6 * i + 4) = w1_ * d4_ * 2.0;
        hessian_.insert(6 * i + 1, 6 * i + 5) = w1_ * d5_ * 2.0;
        hessian_.insert(6 * i + 2, 6 * i + 1) = w1_ * d2_ * 2.0;
        hessian_.insert(6 * i + 2, 6 * i + 2) = (w1_ * 4.0 / 3.0 * d3_ + w2_ * 4.0 * d_) * 2.0;
        hessian_.insert(6 * i + 2, 6 * i + 3) = (w1_ * 6.0 / 4.0 * d4_ + w2_ * 6.0 * d2_) * 2.0;
        hessian_.insert(6 * i + 2, 6 * i + 4) = (w1_ * 8.0 / 5.0 * d5_ + w2_ * 8.0 * d3_) * 2.0;
        hessian_.insert(6 * i + 2, 6 * i + 5) = (w1_ * 10.0 / 6.0 * d6_ + w2_ * 10 * d4_) * 2.0;
        hessian_.insert(6 * i + 3, 6 * i + 1) = w1_ * d3_ * 2.0;
        hessian_.insert(6 * i + 3, 6 * i + 2) = (w1_ * 6.0 / 4.0 * d4_ + w2_ * 6.0 * d2_) * 2.0;
        hessian_.insert(6 * i + 3, 6 * i + 3) = (w1_ * 9.0 / 5.0 * d5_ + w2_ * 12.0 * d3_) * 2.0;
        hessian_.insert(6 * i + 3, 6 * i + 4) = (w1_ * 2.0 * d6_ + w2_ * 18.0 * d4_) * 2.0;
        hessian_.insert(6 * i + 3, 6 * i + 5) = (w1_ * 15.0 / 7.0 * d7_ + w2_ * 24.0 * d5_) * 2.0;
        hessian_.insert(6 * i + 4, 6 * i + 1) = w1_ * d4_ * 2.0;
        hessian_.insert(6 * i + 4, 6 * i + 2) = (w1_ * 8.0 / 5.0 * d5_ + w2_ * 8.0 * d3_) * 2.0;
        hessian_.insert(6 * i + 4, 6 * i + 3) = (w1_ * 2.0 * d6_ + w2_ * 18.0 * d4_) * 2.0;
        hessian_.insert(6 * i + 4, 6 * i + 4) = (w1_ * 16.0 / 7.0 * d7_ + w2_ * 144.0 / 5.0 * d5_) * 2.0;
        hessian_.insert(6 * i + 4, 6 * i + 5) = (w1_ * 20.0 / 8.0 * d8_ + w2_ * 40.0 * d6_) * 2.0;
        hessian_.insert(6 * i + 5, 6 * i + 1) = w1_ * d5_ * 2.0;
        hessian_.insert(6 * i + 5, 6 * i + 2) = (w1_ * 10.0 / 6.0 * d6_ + w2_ * 10 * d4_) * 2.0;
        hessian_.insert(6 * i + 5, 6 * i + 3) = (w1_ * 15.0 / 7.0 * d7_ + w2_ * 24.0 * d5_) * 2.0;
        hessian_.insert(6 * i + 5, 6 * i + 4) = (w1_ * 20.0 / 8.0 * d8_ + w2_ * 40.0 * d6_) * 2.0;
        hessian_.insert(6 * i + 5, 6 * i + 5) = (w1_ * 25.0 / 9.0 * d9_ + w2_ + 400.0 / 7.0 * d7_) * 2.0;
    }
    cout << "Hessian = " << hessian_ << endl;
    cout << "Hessian Generated Successfully" << endl;
    return true;
}

bool QpSmooth::GenerateGradiant()
{
    cout << "----Generate Gradiant-----" << endl;
    for (uint i = 0; i < gradient_.size(); i++) {
        gradient_[i] = 0.0;
    }
    cout << "Gradiant = " << gradient_ << endl;
    cout << "Gradiant Generated Successfully" << endl;
    return true;
}

bool QpSmooth::GenerateLinearMat(const vector<WayPoint>& WayPts, double reso)
{
    cout << "----Generate Linear Matrix-----" << endl;
    if (WayPts.size() <= 2) return false;
    
    // init state
    linearMatrix_.insert(0, 0) = 1.0;
    linearMatrix_.insert(1, 1) = 1.0;
    linearMatrix_.insert(2, 2) = 2.0;
    // final state
    linearMatrix_.insert(3, (knotNum_ - 2) * 6) = 1.0;
    linearMatrix_.insert(3, (knotNum_ - 2) * 6 + 1) = d_;
    linearMatrix_.insert(3, (knotNum_ - 2) * 6 + 2) = d2_;
    linearMatrix_.insert(3, (knotNum_ - 2) * 6 + 3) = d3_;
    linearMatrix_.insert(3, (knotNum_ - 2) * 6 + 4) = d4_;
    linearMatrix_.insert(3, (knotNum_ - 2) * 6 + 5) = d5_;
    linearMatrix_.insert(4, (knotNum_ - 2) * 6 + 1) = 1.0;
    linearMatrix_.insert(4, (knotNum_ - 2) * 6 + 2) = 2.0 * d_;
    linearMatrix_.insert(4, (knotNum_ - 2) * 6 + 3) = 3.0 * d2_;
    linearMatrix_.insert(4, (knotNum_ - 2) * 6 + 4) = 4.0 * d3_;
    linearMatrix_.insert(4, (knotNum_ - 2) * 6 + 5) = 5.0 * d4_;
    linearMatrix_.insert(5, (knotNum_ - 2) * 6 + 2) = 2.0;
    linearMatrix_.insert(5, (knotNum_ - 2) * 6 + 3) = 6.0 * d_;
    linearMatrix_.insert(5, (knotNum_ - 2) * 6 + 4) = 12.0 * d2_;
    linearMatrix_.insert(5, (knotNum_ - 2) * 6 + 5) = 20.0 * d3_;
    // connect state
    for (uint i = 0; i < knotNum_ - 2; i++) { // knotNum_ is larger than 2
        // function value continuous 
        // first
        linearMatrix_.insert(6 + 3 * i, 6 * i) = 1;
        linearMatrix_.insert(6 + 3 * i, 6 * i + 1) = d_;
        linearMatrix_.insert(6 + 3 * i, 6 * i + 2) = d2_;
        linearMatrix_.insert(6 + 3 * i, 6 * i + 3) = d3_;
        linearMatrix_.insert(6 + 3 * i, 6 * i + 4) = d4_;
        linearMatrix_.insert(6 + 3 * i, 6 * i + 5) = d5_;
        // second
        linearMatrix_.insert(6 + 3 * i, 6 * i + 6) = -1.0;

        // function 1st continuous 
        // first
        linearMatrix_.insert(6 + 3 * i + 1, 6 * i + 1) = 1;
        linearMatrix_.insert(6 + 3 * i + 1, 6 * i + 2) = 2.0 * d_;
        linearMatrix_.insert(6 + 3 * i + 1, 6 * i + 3) = 3.0 * d2_;
        linearMatrix_.insert(6 + 3 * i + 1, 6 * i + 4) = 4.0 * d3_;
        linearMatrix_.insert(6 + 3 * i + 1, 6 * i + 5) = 5.0 * d4_;
        // second
        linearMatrix_.insert(6 + 3 * i + 1, 6 * i + 7) = -1.0;

        // function 2rd continuous
        // first
        linearMatrix_.insert(6 + 3 * i + 2, 6 * i + 2) = 2.0;
        linearMatrix_.insert(6 + 3 * i + 2, 6 * i + 3) = 6.0 * d_;
        linearMatrix_.insert(6 + 3 * i + 2, 6 * i + 4) = 12.0 * d2_;
        linearMatrix_.insert(6 + 3 * i + 2, 6 * i + 5) = 20.0 * d3_;
        // second
        linearMatrix_.insert(6 + 3 * i + 2, 6 * i + 8) = -2.0;
        // cout << "connect pt index i = " << 6 + 3 * i << "index j = " << 6 * i + 8 << endl;
    }
    // upper & lower bound
    for (uint i = 1; i < WayPts.size() - 1; i++) {
        uint knotId = i / segNum_;
        uint segId = i % segNum_;
        // cout << "bound i = " << i << " knotId = " << knotId << endl;
        double dTemp = 0.0;
        double ds = reso * segId;
        double ds2 = ds * ds;
        double ds3 = ds2 * ds;
        double ds4 = ds3 * ds;
        double ds5 = ds4 * ds;
        linearMatrix_.insert(knotNum_ * 3 + i - 1, 6 * knotId) = 1;
        linearMatrix_.insert(knotNum_ * 3 + i - 1, 6 * knotId + 1) = ds;
        linearMatrix_.insert(knotNum_ * 3 + i - 1, 6 * knotId + 2) = ds2;
        linearMatrix_.insert(knotNum_ * 3 + i - 1, 6 * knotId + 3) = ds3;
        linearMatrix_.insert(knotNum_ * 3 + i - 1, 6 * knotId + 4) = ds4;
        linearMatrix_.insert(knotNum_ * 3 + i - 1, 6 * knotId + 5) = ds5;
        // cout << "bound pt index i = " << knotNum_ * 3 + i - 1 << " index j = " << 6 * knotId<< endl;
    }
    cout << "Linear Matrix" << linearMatrix_ << endl;
    cout << "Linear Matrix Generated Successfully" << endl;
    return true;
}
bool QpSmooth::GenerateBounds(const vector<WayPoint>& WayPts)
{
    cout << "----Generate Bounds-----" << endl;
    // concide with Linear Matrix
    // init state
    lowerBound_[0] = -1.0; // init a position offset
    upperBound_[0] = -1.0;
    lowerBound_[1] = 0.0;
    upperBound_[1] = 0.0;
    lowerBound_[2] = 0.0;
    upperBound_[2] = 0.0;

    // final state
    lowerBound_[3] = 0.0;
    upperBound_[3] = 0.0;
    lowerBound_[4] = 0.0;
    upperBound_[4] = 0.0;
    lowerBound_[5] = 0.0;
    upperBound_[5] = 0.0;

    // connect state
    for (uint i = 0; i < knotNum_ - 2; i++) { // knotNum_ is larger than 2
        lowerBound_[6 + 3 * i] = 0.0;
        upperBound_[6 + 3 * i] = 0.0;
        lowerBound_[6 + 3 * i + 1] = 0.0;
        upperBound_[6 + 3 * i + 1] = 0.0;
        lowerBound_[6 + 3 * i + 2] = 0.0;
        upperBound_[6 + 3 * i + 2] = 0.0;
    }
    // uppper & lower bound
    for (uint i = 1; i < WayPts.size() - 1; i++) {
        lowerBound_[knotNum_ * 3 + i - 1] = WayPts[i].leftBound;
        upperBound_[knotNum_ * 3 + i - 1] = WayPts[i].rightBound;
    }
    cout << "lowerBound = " << lowerBound_ << endl;
    cout << "upperBound = " << upperBound_ << endl;
    cout << "Bounds Generated Successfully" << endl;
    return true;
}

bool QpSmooth::SetQpSolver(OsqpEigen::Solver& solver)
{
    cout << "----Set Qp Solver-----" << endl;
    solver.data()->setNumberOfVariables(numOfVar_);
    solver.data()->setNumberOfConstraints(numOfCons_);
    if (!solver.data()->setHessianMatrix(hessian_)) return false;
    if (!solver.data()->setGradient(gradient_)) return false;
    if (!solver.data()->setLinearConstraintsMatrix(linearMatrix_)) return false;
    if (!solver.data()->setLowerBound(lowerBound_)) return false;
    if (!solver.data()->setUpperBound(upperBound_)) return false;
    cout << "----Set Qp Solver Successfully-----" << endl;
    return true;
}

void QpSmooth::DrawSmoothPath(const Eigen::VectorXd& QPSolution, const vector<WayPoint>& sceWayPt, const double reso,
    vector<double>& lPath)
{
    lPath.resize(sceWayPt.size());
    for (uint i = 0; i < sceWayPt.size(); i++) {
        uint knotId = i / segNum_;
        uint segId = i % segNum_;
        if (i < sceWayPt.size() - 1) {
            double ds = reso * segId;
            double ds2 = ds * ds;
            double ds3 = ds2 * ds;
            double ds4 = ds3 * ds;
            double ds5 = ds4 * ds;
            lPath[i] = QPSolution[knotId * 6] + QPSolution[knotId * 6 + 1] * ds + QPSolution[knotId * 6 + 2] * ds2 +
                QPSolution[knotId * 6 + 3] * ds3 + QPSolution[knotId * 6 + 4] * ds4 + QPSolution[knotId * 6 + 5] * ds5;
        }
        else {
            double ds = reso * (segId + segNum_);
            double ds2 = ds * ds;
            double ds3 = ds2 * ds;
            double ds4 = ds3 * ds;
            double ds5 = ds4 * ds;
            lPath[i] = QPSolution[(knotId - 1) * 6] + QPSolution[(knotId - 1) * 6 + 1] * ds + QPSolution[(knotId - 1) * 6 + 2] * ds2 +
                QPSolution[(knotId - 1) * 6 + 3] * ds3 + QPSolution[(knotId - 1) * 6 + 4] * ds4 + QPSolution[(knotId - 1) * 6 + 5] * ds5;
        }
    }
}

int main() {
    uint knotNum = 4;
    uint segNum = 10;
    double ptReso = 0.25;
    
    // Init Senarios 
    Scenarios straightSce;
    vector<WayPoint> sceWayPt;
    vector<double> s, l, leftBound, rightBound;
    straightSce.InitStraightNudgeSce(sceWayPt);
    straightSce.DrawSceInFrenet(sceWayPt, s, l, leftBound, rightBound);

    // Smoother
    QpSmooth smoother;
    smoother.SetMatSize(knotNum, segNum);
    smoother.GenerateHessian(segNum * ptReso);
    smoother.GenerateGradiant();
    smoother.GenerateLinearMat(sceWayPt, ptReso); // 0.5 resolution (each segment).
    smoother.GenerateBounds(sceWayPt);
    
    // set the solver
    OsqpEigen::Solver solver;
    solver.settings()->setWarmStart(true);

    // instantiate the solver
    if (smoother.SetQpSolver(solver)) {
        if (!solver.initSolver()) return 1;
    }
    else {
        cout << "initilize QP solver failed" << endl;
        return 1;
    }

    // solve
    solver.solve();
    Eigen::VectorXd QPSolution;
    QPSolution = solver.getSolution();
    cout << QPSolution << endl;
    vector<double> lPath;
    smoother.DrawSmoothPath(QPSolution, sceWayPt, ptReso, lPath);

    // cout << "s = " << s.size() << " l = " << l.size() << " lPath = " << lPath.size() << " left =" << leftBound.size()
    //    << " right =" << rightBound.size();
    plt::plot(s, l);
    plt::plot(s, lPath, "r");
    plt::plot(s, leftBound);
    plt::plot(s, rightBound);
    plt::save("straightSce.png", 300);
}
