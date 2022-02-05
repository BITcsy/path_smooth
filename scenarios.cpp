#include "scenarios.h"

void Scenarios::InitStraightNudgeSce(vector<WayPoint>& straightSce)
{
	uint knotNum = 4;
	uint segNum = 10;
	double reso = 0.25;
	straightSce.clear();
	uint wayPointNum = segNum * (knotNum - 1) + 1;
	straightSce.resize(wayPointNum);
	for (uint i = 0; i < wayPointNum; i++) {
		straightSce[i].x = reso * i;
		straightSce[i].y = 0.0;
		straightSce[i].theta = 0.0;
		straightSce[i].s = reso * i;
		straightSce[i].l = 0.0;
		straightSce[i].leftBound = -2.5;
		straightSce[i].rightBound = 2.5;
	}

	straightSce[10].leftBound = 0.3;
	straightSce[11].leftBound = 0.3;
	straightSce[19].rightBound = -0.3;
	straightSce[20].rightBound = -0.3;
	cout << "InitStraightNudgeSce Finished" << endl;
	return;
}

void Scenarios::DrawSceInFrenet(const vector<WayPoint>& Sce, vector<double>& s, 
	vector<double>& l, vector<double>& leftBound, vector<double>& rightBound)
{
	uint wayPtNum = Sce.size();
	s.resize(wayPtNum);
	l.resize(wayPtNum);
	leftBound.resize(wayPtNum);
	rightBound.resize(wayPtNum);
	for (uint i = 0; i < wayPtNum; i++) {
		s[i] = Sce[i].s;
		l[i] = Sce[i].l;
		leftBound[i] = Sce[i].leftBound;
		rightBound[i] = Sce[i].rightBound;
	}
	cout << "DrawSceInFrenet Convert Finished" << endl;
	/*plt::plot(s, l);
	plt::plot(s, leftBound);
	plt::plot(s, rightBound);
	plt::save("straightSce.png");*/
	return;
}