#pragma once
#include "types.h"

// point distance
#define DIST(a, b)  sqrt(((a).x - (b).x) * ((a).x - (b).x) + \
    ((a).y - (b).y) * ((a).y - (b).y))

#define FZERO  -1e-8
#define CV_PI   3.1415926535897932384626433832795
#define CV_LOG2 0.69314718055994530941723212145818

#ifndef CV_SWAP
#define CV_SWAP(a,b,t) ((t) = (a), (a) = (b), (b) = (t))
#endif

#ifndef MIN
#  define MIN(a,b)  ((a) > (b) ? (b) : (a))
#endif

#ifndef MAX
#  define MAX(a,b)  ((a) < (b) ? (b) : (a))
#endif

#define ROT_CC(px, py, sinv, cosv, pt)  	\
	pt.x = px * cosv + py * sinv; \
	pt.y = -px * sinv + py * cosv

#define ROT_CW(px, py, sinv, cosv, pt)   \
	pt.x = px * cosv - py * sinv; \
	pt.y = px * sinv + py * cosv

//#define USE_ROT_BOX 

// Operates on Edges
// Edge i  is pi -> pi+1
// Edge -1 is pN-1 -> p0

typedef struct PolySeg {
	int type; // 0: simple line strip, 1: segment
	//std::vector<int> points;
	int start, end;
	// bound point, left top, right bottom
	Point tl, br;
	// bound
	Point rbox[4];
#ifdef USE_ROT_BOX
	// local coord
	Point pmin, pmax;
	FLOAT cosv, sinv;
#endif
	// LOD
	PolySeg *left, *right;
	PolySeg() {
		type = 0;
		this->start = 0;
		this->end = 0;
		left = right = NULL;
	}
	PolySeg(int start, int end) {
		type = 0;
		this->start = start;
		this->end = end;
		left = right = NULL;
	}
	void setRBOX(Point vertices[4]) {
		for (int i = 0; i < 4; i++) {
			rbox[i].x = vertices[i].x;
			rbox[i].y = vertices[i].y;
		}
	}
	inline int pointInRect(float x, float y) const {
		return tl.x <= x and x <= br.x and
			tl.y <= y and y <= br.y;
	}
	inline int interectRect(float x1, float y1, float x2, float y2) const {
		float t;
		if (x1 > x2) {
			CV_SWAP(x1, x2, t);
		}
		float W = MIN(x2, br.x) - MAX(x1, tl.x);
		if (W < FZERO) {
			return 0;
		}
		if (y1 > y2) {
			CV_SWAP(y1, y2, t);
		}
		float H = MIN(y2, br.y) - MAX(y1, tl.y);
		return H > FZERO;
	}
#ifndef USE_ROT_BOX
	inline int pointInRBOX(float x, float y) const {
		const Point& A = rbox[0];
		const Point& B = rbox[1];
		const Point& C = rbox[2];
		const Point& D = rbox[3];
		float a = (B.x - A.x)*(y - A.y) - (B.y - A.y)*(x - A.x);
		float b = (C.x - B.x)*(y - B.y) - (C.y - B.y)*(x - B.x);
		float c = (D.x - C.x)*(y - C.y) - (D.y - C.y)*(x - C.x);
		float d = (A.x - D.x)*(y - D.y) - (A.y - D.y)*(x - D.x);
		if ((a > 0 && b > 0 && c > 0 && d > 0) || (a < 0 && b < 0 && c < 0 && d < 0)) {
			return 1;
		}
		return 0;
	}
#else
	inline int pointInRBOX(float x, float y) const {
		Point pt;
		ROT_CC(x, y, sinv, cosv, pt);
		int ret= pmin.x <= pt.x and pt.x <= pmax.x and
			pmin.y <= pt.y and pt.y <= pmax.y;
		if (ret == 0) {
			return ret;
		}
		return ret;
	}
#endif 
	inline int interectRBOXEdge(Segment& edge, long long* cnt = NULL) const {
		int numPts = 0;
		Point pi0, pi1;
		for (int j = 0; j < 4; j++)
		{
			int p2Idx = (j + 1) % 4;
			Segment sub(Point(rbox[j].x, rbox[j].y),
				Point(rbox[p2Idx].x, rbox[p2Idx].y));
			//numPts = findIntersection(sub, edge, pi0, pi1);
			numPts = checkIntersection(sub, edge);
			if (cnt) *cnt += 1;
			if (numPts)
			{
				return 1;
			}
		}
		return 0;
	}
	inline int interectRBOX(Segment& edge, long long* cnt=NULL) const {
		if (interectRBOXEdge(edge, cnt)) {
			return 1;
		}
		if (pointInRBOX(edge.p1.x, edge.p1.y)) {
			return 1;
		}
		if (pointInRBOX(edge.p2.x, edge.p2.y)) {
			return 1;
		}
		return 0;
	}

	inline int interectRBOX(const PolySeg& clip, long long* cnt = NULL) const {
		int numPts = 0;
		// clip x sub or 
		// clip in sub
		for (int j = 0; j < 4; j++)
		{
			int p2Idx = (j + 1) % 4;
			Segment edge(Point(clip.rbox[j].x, clip.rbox[j].y),
				Point(clip.rbox[p2Idx].x, clip.rbox[p2Idx].y));
			numPts = interectRBOXEdge(edge, cnt);
			if (numPts) return 1;
			if (pointInRBOX(clip.rbox[j].x, clip.rbox[j].y)) {
				return 1;
			}
		}
		// sub in clip
		for (int j = 0; j < 4; j++)
		{
			if (clip.pointInRBOX(rbox[j].x, rbox[j].y)) {
				return 1;
			}
		}
		return 0;
	}
	inline int numEdges()
	{
		return end - start + 1;
	}
}PolySeg;

typedef struct CrossEdge {
	int sIdx, cIdx;     // subject edge index & clipper edge index
	int numPts;         // num of cross points
	Point p0, p1;
	CrossEdge() {
		sIdx = cIdx = 0;
		numPts = 0;
	}
}CrossEdge;

class PolyStruct{
public:
	const Poly& mPoly;
	std::vector<PolySeg> mSegs;
	const Poly* mClip;
	// param
	int N;
	int K;
	int MinSegEdges;
	int MaxDepth;
	std::vector< FLOAT > mLens;
	std::vector< std::vector<PolySeg>* > mLod;
	std::vector< CrossEdge > mCrossEdges;
	// debug
	long long nInterCounter;
	PolyStruct(const Poly& poly, int K = 16, int minSegEdges = 4,
		int maxDepth=4) : mPoly(poly)
	{
		this->K = K;
		N = poly.size();
		MaxDepth = maxDepth;
		MinSegEdges = minSegEdges;
	}
	~PolyStruct() {
		_Reset();
	}
	void BuildLod();
	void BuildLod(const Poly& poly, std::vector<PolySeg>& segs, std::vector< std::vector<PolySeg>* >& lod);
	void FindCrossEdges(const Poly& clip);
	void FindCrossEdgesLOD(const Poly& clip);
	void Viz(const char* savePath, Polys& extra);
protected:
	void _Reset();
	void _FindCrossEdges(Segment& edge, int cIdx, std::vector< CrossEdge >& crossEdges);
	void _FindCrossEdgesTree(const PolySeg& seg, Segment& edge, int cIdx, std::vector< CrossEdge >& crossEdges);
	void _FindCrossSegs(const PolySeg& subSeg, const Poly& clipPoly, const PolySeg& clipSeg, std::vector< CrossEdge >& crossEdges);
};
