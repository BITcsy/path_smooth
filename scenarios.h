/* oxalate-c 20220205
Genearte Scenarios, including maps, objs etc.

*/
#include "basic_data_def.h"
#include "matplotlibcpp.h"
#include <vector>
#ifndef SCENARIOS_H
#define SCENARIOS_H
using namespace std;
namespace plt = matplotlibcpp;

class Scenarios {
public:
	// Straight path with some static objs 
	void InitStraightNudgeSce(vector<WayPoint>& straightSce);
	void DrawSceInFrenet(const vector<WayPoint>& Sce, vector<double>& s,
		vector<double>& l, vector<double>& leftBound, vector<double>& rightBound);
private:

};

#endif // !SCENARIOS_H
