#ifndef GREINER_H
#define GREINER_H

#include <iostream>
#include <vector>
#include <limits>
#include "common/segment.h"
#include "common/polygon.h"
#include "common/point.h"

using namespace std;

enum BoolOpType { INTERSECTION, UNION, DIFF, XOR };

struct Vertex {
	double x, y;
	Vertex *prev, *next;
	Segment s; // Segment which first point is the vertex
	bool intersect; // is the vertex an intersection point?
	bool entry; // is the intersection point an entry point on the other polygon?
	bool processed; // for step 3
	Vertex *neighbor;
	double alpha;
	Vertex (double xp, double yp, Segment sp, bool i, double a = numeric_limits<double>::max());
	Vertex () {}
};
ostream& operator<< (ostream& o, const Vertex& v);

class GreinerContour {
public:
	GreinerContour (Contour& c);
	//~GreinerContour () { deleteIntersections (); }
	~GreinerContour() {  }
	Vertex *firstVertex () { return &v[0]; }
	/** @brief Insert vertex (intersection) v before the vertex pointed by vp */
	Vertex *insert (const Vertex& v, Vertex *vp);
	void deleteIntersections ();
	bool intersectBoundingbox (GreinerContour& gc) const;
private:
	/** @brief It holds the original vertices of the polygon */
	vector<Vertex> v;
	/** @brief Number of intersection points */
	int nint;
	/** @ brief bounding box */
	Point minbox;
	Point maxbox;
};
ostream& operator<< (ostream& o, GreinerContour& gc);

class GreinerHormann {
public:
	GreinerHormann (Polygon& p1, Polygon& p2);
	~GreinerHormann ();
	int boolop (BoolOpType op, Polygon& result);
private:
	vector<GreinerContour*> gp1;
	vector<GreinerContour*> gp2;
	Polygon& subject;
	Polygon& clipping;
	int boolop (BoolOpType op, GreinerContour& gc1, GreinerContour& gc2, Polygon& result);
};

#endif
