/* oxalate-c 20220205
Basic data struct definition
*/
#ifndef BASIC_DATA_DEF_H
#define BASIC_DATA_DEF_H
#include <vector>
typedef unsigned int unit;

using namespace std;

struct WayPoint {
	// Cartasian
	double x;
	double y;
	double theta;
	// Frenet
	double s;
	double l;
	double leftBound;
	double rightBound;
	WayPoint() {};
	WayPoint(double _x, double _y, double _theta, double _s, double _l, double _leftBound, double _rightBound) {
		x = _x;
		y = _y;
		theta = _theta;
		s = _s;
		l = _l;
		leftBound = _leftBound;
		rightBound = _rightBound;
	}
};

#endif // !BASIC_DATA_DEF_H