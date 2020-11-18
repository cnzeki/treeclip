#include "PolyStruct.h"
#include "svg.h"
#include "common/utilities.h"
#include "common/ticker.h"


#define   SegFitting  SegFittingSimple

#define USE_OPENCV  0

void PrintPoints(int n, Point* pts) {
	for (int i = 0; i < n; i++)
	{
		printf("%.2f, %.2f\n", pts[i].x, pts[i].y);
	}
}
#if USE_OPENCV
#include <opencv2/core.hpp>
#include <opencv2/imgproc.hpp>

#define HAVE_OPENCV 1

//#ifndef CVLIBVER
//#define CVLIBVER "346"
//#endif
//
//#ifdef _DEBUG 
//#define CVLIBEXT "d.lib"
//#else
//#define CVLIBEXT ".lib"
//#endif
//
//#define CVLIBNAME(name,ver,ext) "opencv_" name ver ext

//#pragma comment(lib,CVLIBNAME("core",CVLIBVER,CVLIBEXT))
//#pragma comment(lib,CVLIBNAME("imgproc",CVLIBVER,CVLIBEXT))

void SegFittingOpenCV(const Poly& poly, PolySeg& seg, int minSegEdges,
	int maxDepth, int depth) {
	int numEdges = seg.numEdges();
	if (numEdges < minSegEdges) {
		return;
	}
	std::vector<cv::Point2f> pts(numEdges+1);
	cv::Point2f* ptsPtr = &pts[0];
	FLOAT minx, maxx, miny, maxy;
	minx = maxx = poly[seg.start].x;
	miny = maxy = poly[seg.start].y;
	for (int i = 0; i < numEdges + 1; i++) {
		const Point& pt = poly[(i+seg.start) % poly.size()];
		ptsPtr[i].x = pt.x;
		ptsPtr[i].y = pt.y;
		minx = MIN(minx, pt.x);
		maxx = MAX(maxx, pt.x);
		miny = MIN(miny, pt.y);
		maxy = MAX(maxy, pt.y);
	}
	// bound
	seg.tl = Point(minx, miny);
	seg.br = Point(maxx, maxy);
	// rbox
	cv::Point2f vertices[4];
	cv::RotatedRect rbox; 
	rbox = cv::minAreaRect(pts);
	rbox.points(vertices);
	for (int i = 0; i < 4; i++) {
		seg.rbox[i].x = vertices[i].x;
		seg.rbox[i].y = vertices[i].y;
	}
	//PrintPoints(4, seg.rbox);
	seg.type = 1;
}
#endif


void SimpleFitting(int n, Point* ptsPtr, PolySeg& seg, FLOAT v[2]) {
	FLOAT modv = 1.0f / sqrt(v[0] * v[0] + v[1] * v[1]);
	FLOAT sinv = v[1] * modv;
	FLOAT cosv = v[0] * modv;
	Point* vertices = seg.rbox;

	// rotate counter-clockwise
	FLOAT minx, maxx, miny, maxy;
	minx = miny = FLT_MAX;
	maxx = maxy = -FLT_MAX;
	Point pt;
	for (int i = 0; i < n; i++) {
		ROT_CC(ptsPtr[i].x, ptsPtr[i].y, sinv, cosv, pt);
		minx = MIN(minx, pt.x);
		maxx = MAX(maxx, pt.x);
		miny = MIN(miny, pt.y);
		maxy = MAX(maxy, pt.y);
	}
	// vertices
	ROT_CW(minx, miny, sinv, cosv, vertices[0]);
	ROT_CW(maxx, miny, sinv, cosv, vertices[1]);
	ROT_CW(maxx, maxy, sinv, cosv, vertices[2]);
	ROT_CW(minx, maxy, sinv, cosv, vertices[3]);
	// bound
#ifdef USE_ROT_BOX
	seg.pmin.x = minx; seg.pmin.y = miny;
	seg.pmax.x = maxx; seg.pmax.y = maxy;
	seg.sinv = sinv;
	seg.cosv = cosv;
#endif
}


void SimpleFitting(int n, Point* ptsPtr, Point vertices[4]) {
	Point& start = ptsPtr[0];
	Point& end = ptsPtr[n - 1];
	FLOAT v[2] = { end.x - start.x, end.y - start.y };
	FLOAT modv = 1.0f / sqrt(v[0] * v[0] + v[1] * v[1]);
	FLOAT sinv = v[1] * modv;
	FLOAT cosv = v[0] * modv;
	// rotate counter-clockwise
	FLOAT minx, maxx, miny, maxy;
	minx = miny = FLT_MAX;
	maxx = maxy = -FLT_MAX;
	Point pt;
	for (int i = 0; i < n; i++) {
		ROT_CC(ptsPtr[i].x, ptsPtr[i].y, sinv, cosv, pt);
		minx = MIN(minx, pt.x);
		maxx = MAX(maxx, pt.x);
		miny = MIN(miny, pt.y);
		maxy = MAX(maxy, pt.y);
	}
	// 
	ROT_CW(minx, miny, sinv, cosv, vertices[0]);
	ROT_CW(maxx, miny, sinv, cosv, vertices[1]);
	ROT_CW(maxx, maxy, sinv, cosv, vertices[2]);
	ROT_CW(minx, maxy, sinv, cosv, vertices[3]);
}

void SegFittingSimple(const Poly& poly, PolySeg& seg, int minSegEdges,
	int maxDepth, int depth) {
	int numEdges = seg.numEdges();
	if (numEdges < minSegEdges) {
		//return;
	}
	//printf("Point size: %d\n", sizeof(Point2D));
	std::vector<Point> pts(numEdges + 1);
	Point* ptsPtr = &pts[0];
	FLOAT minx, maxx, miny, maxy;
	minx = maxx = poly[seg.start].x;
	miny = maxy = poly[seg.start].y;
	for (int i = 0; i < numEdges + 1; i++) {
		const Point& pt = poly[(i + seg.start) % poly.size()];
		ptsPtr[i].x = pt.x;
		ptsPtr[i].y = pt.y;
		minx = MIN(minx, pt.x);
		maxx = MAX(maxx, pt.x);
		miny = MIN(miny, pt.y);
		maxy = MAX(maxy, pt.y);
	}
	// bound
	seg.tl = Point(minx, miny);
	seg.br = Point(maxx, maxy);
	// rbox
	//SimpleFitting(pts.size(), &pts[0], seg.rbox);
	Point& start = ptsPtr[0];
	Point& end = ptsPtr[pts.size() - 1];
	Point* vertices = seg.rbox;
	FLOAT v[2] = { end.x - start.x, end.y - start.y };
	SimpleFitting(pts.size(), &pts[0], seg, v);
	//PrintPoints(4, seg.rbox);
	seg.type = 1;
}

std::vector<PolySeg>* GroupSegments(const Poly& poly, std::vector<PolySeg>& segs, int K) {
	
	if (segs.size() <= K) {
		return &segs;
	}
	int n = segs.size()/2;
	int nceil = (segs.size() + 1) / 2;
	std::vector<PolySeg>* pOutputs = new std::vector<PolySeg>(nceil);
	std::vector<PolySeg>& outputs = *pOutputs;
	PolySeg* segPtr = &segs[0];
	PolySeg* outPtr = &outputs[0];
	// local vars
	std::vector<Point> pts2(8);
#if USE_OPENCV
	std::vector<cv::Point2f> pts(8);
	cv::Point2f vertices[4];
	cv::RotatedRect rbox;
#endif	
	for (int i = 0; i < n; i++) {
		// merge
		PolySeg& left = segPtr[i * 2];
		PolySeg& right = segPtr[i * 2+1];
		// 
		FLOAT minx, maxx, miny, maxy;
		minx = MIN(left.tl.x, right.tl.x);
		maxx = MAX(left.br.x, right.br.x);
		miny = MIN(left.tl.y, right.tl.y);
		maxy = MAX(left.br.y, right.br.y);
		
		for (int j = 0; j < 4; j++) {
#if USE_OPENCV
			cv::Point2f pt(left.rbox[j].x, left.rbox[j].y);
			pts[j] = pt;
#endif
			pts2[j] = left.rbox[j];
		}
		for (int j = 0; j < 4; j++) {
#if USE_OPENCV
			cv::Point2f pt(right.rbox[j].x, right.rbox[j].y);
			pts[4+j] = pt;
#endif
			pts2[4 + j] = right.rbox[j];
		}
		
		// bound
		PolySeg seg(left.start, right.end);
		seg.tl = Point(minx, miny);
		seg.br = Point(maxx, maxy);
		// rbox
#if USE_OPENCV		
		rbox = cv::minAreaRect(pts);
		rbox.points(vertices);
		for (int j = 0; j < 4; j++) {
			seg.rbox[j].x = vertices[j].x;
			seg.rbox[j].y = vertices[j].y;
		}
#else
		//SimpleFitting(pts2.size(), &pts2[0], seg.rbox);
		const Point& start = poly[left.start];
		const Point& end = poly[(right.end + 1) % poly.size()];
		Point* vertices = seg.rbox;
		FLOAT v[2] = { end.x - start.x, end.y - start.y };
		SimpleFitting(pts2.size(), &pts2[0], seg, v);
#endif
		seg.type = 2;
		seg.left = &left;
		seg.right = &right;
		outPtr[i] = seg;
	}
	if (segs.size() % 2) {
		outputs[nceil-1] = segs[segs.size()-1];
	}
	return pOutputs;
}

void PolyStruct::_Reset() {
	for (int i = 0; i < mLod.size(); i++) {
		//delete mLod[i];
		mLod[i]->clear();
	}
	mLod.clear();
}

void PolyStruct::BuildLod(const Poly& poly, std::vector<PolySeg>& segs, std::vector< std::vector<PolySeg>* >& lod) {
	int N = poly.size();
	// Figure out seg length
	int maxSegs = 0; 
	int edgesPerSeg = 0; 
	if (MaxDepth > 0) {
		maxSegs = K * (2 << MaxDepth);
		edgesPerSeg = N / maxSegs;
		edgesPerSeg = max(edgesPerSeg, MinSegEdges);
	}
	else {
		maxSegs = (N + MinSegEdges - 1) / MinSegEdges;
		edgesPerSeg = MinSegEdges;
	}
	
	// find a segment
	int startIdx = 0;
	int curIdx = startIdx;
	while (curIdx < N)
	{
		// add current edge
		int numEdges = curIdx - startIdx + 1;
		if (curIdx != N - 1 && numEdges < edgesPerSeg) {
			curIdx += 1;
			continue;
		}

		// new segment
		PolySeg seg(startIdx, curIdx);
		SegFitting(poly, seg, MinSegEdges, MaxDepth, 1);
		segs.push_back(seg);
		// move next
		startIdx = curIdx + 1;
		curIdx = startIdx;
	}

	// recursive fitting
	std::vector<PolySeg>* input = &segs;
	lod.push_back(input);
	do {
		if (input->size() <= K) {
			break;
		}
		std::vector<PolySeg>* outputs = GroupSegments(poly, *input, K);
		lod.push_back(outputs);
		input = lod[lod.size() - 1];
	} while (1);

#if 0
	double elapsed = double(clock() - startTime) * 1000 / CLOCKS_PER_SEC;
	printf("\n");
	printf("Param: K=%d, MaxDepth=%d, MinSegEdges=%d\n", K, MaxDepth, MinSegEdges);
	printf("Data : N=%d, edgesPerSeg=%d\n", N, edgesPerSeg);
	printf("LOD  : %.2f ms\n", elapsed);
	printf("Fit  : %.2f ms\n", tlod);
	printf("Merge: %.2f ms\n", elapsed-tlod);
#endif 
}

void PolyStruct::BuildLod()
{
	_Reset();
	BuildLod(mPoly, mSegs, mLod);
}

void PolyStruct::_FindCrossEdgesTree(const PolySeg& seg, Segment& edge, int cIdx, std::vector< CrossEdge >& crossEdges) {
	Point pi0, pi1;
	if (seg.type < 2) {
		for (int j = seg.start; j <= seg.end; j++)
		{
			Segment sub(mPoly[j], mPoly[(j + 1) % N]);
			int numPts = findIntersection(sub, edge, pi0, pi1);
			nInterCounter += 1;
			if (numPts)
			{
				CrossEdge newEdge;
				newEdge.sIdx = j;
				newEdge.cIdx = cIdx;
				newEdge.numPts = numPts;
				newEdge.p0 = pi0;
				newEdge.p1 = pi1;

				crossEdges.push_back(newEdge);
				
				//printf("Cross %d x %d\n", j, cIdx);
			}
		}
	}
	else {
		// rect
		if (not seg.interectRect(edge.p1.x, edge.p1.y, edge.p2.x, edge.p2.y)) {
			return;
		}
		// rbox
		if (not seg.interectRBOX(edge, &nInterCounter)) {
			return;
		}
		// children
		if (seg.left) {
			_FindCrossEdgesTree(*seg.left, edge, cIdx, crossEdges);
		}
		if (seg.right) {
			_FindCrossEdgesTree(*seg.right, edge, cIdx, crossEdges);
		}
	}
}

void PolyStruct::_FindCrossEdges(Segment& edge, int cIdx, std::vector< CrossEdge >& crossEdges) {
	std::vector<PolySeg>& root = *mLod[mLod.size() - 1];
	// travel through tree root
	for (int i = 0; i < root.size();i++) {
		PolySeg& seg = root[i];
		int level = mLod.size() - 1 - i;
		_FindCrossEdgesTree(seg, edge, cIdx, crossEdges);
	}
}

void PolyStruct::_FindCrossSegs(const PolySeg& subSeg, const Poly& clipPoly,
	const PolySeg& clipSeg, std::vector< CrossEdge >& crossEdges) {
	Point pi0, pi1;
	if (clipSeg.type < 2) {
		for (int j = clipSeg.start; j <= clipSeg.end; j++)
		{
			Segment edge(clipPoly[j], clipPoly[(j + 1) % clipPoly.size()]);
			_FindCrossEdgesTree(subSeg, edge, j, crossEdges);
		}
	}
	else {
		// rect
		if (not clipSeg.interectRect(subSeg.tl.x, subSeg.tl.y, subSeg.br.x, subSeg.br.y)) {
			return;
		}
		// rbox
		if (not clipSeg.interectRBOX(subSeg, &nInterCounter)) {
			return;
		}
		// children
		if (clipSeg.left) {
			if (subSeg.left) {
				_FindCrossSegs(*subSeg.left, clipPoly, *clipSeg.left, crossEdges);
			}
			if (subSeg.right) {
				_FindCrossSegs(*subSeg.right, clipPoly, *clipSeg.left, crossEdges);
			}
			if (subSeg.type < 2) {
				_FindCrossSegs(subSeg, clipPoly, *clipSeg.left, crossEdges);
			}
		}
		if (clipSeg.right) {
			if (subSeg.left) {
				_FindCrossSegs(*subSeg.left, clipPoly, *clipSeg.right, crossEdges);
			}
			if (subSeg.right) {
				_FindCrossSegs(*subSeg.right, clipPoly, *clipSeg.right, crossEdges);
			}
			if (subSeg.type < 2) {
				_FindCrossSegs(subSeg, clipPoly, *clipSeg.right, crossEdges);
			}
		}
	}
}

void PolyStruct::FindCrossEdges(const Poly& clip) {
	mCrossEdges.clear();
	mClip = &clip;
	nInterCounter = 0;

	// loop all cliper edges
	int testEdges = clip.size();
	for (int i = 0; i < clip.size() && i < testEdges; i++)
	{
		Segment edge(clip[i], clip[(i + 1) % clip.size()]);
		_FindCrossEdges(edge, i, mCrossEdges);
	}
	//printf("\n\nSubject points: %d\n", mPoly.size());
	//printf("Clipper points: %d\n", clip.size());
	//printf("Inter    edges: %d\n", mCrossEdges.size());
	//printf("Inter    count: %d\n", nInterCounter);
	//printf("Speedup       : %.2f\n", double(mPoly.size())*testEdges / nInterCounter);
	//printf("elapsed       : %.2f\n", elapsed);
}

void PolyStruct::FindCrossEdgesLOD(const Poly& clip) {
	mCrossEdges.clear();
	mClip = &clip;
	nInterCounter = 0;
	std::vector<PolySeg> segs;
	std::vector< std::vector<PolySeg>* > lod;

	BuildLod(clip, segs, lod);

	GetTicker().record("Build Clip LOD");

	// clip root
	std::vector<PolySeg>& clipRoot = *lod[lod.size() - 1];
	// travel through tree root
	int testEdges = clip.size();
	for (int i = 0; i < clipRoot.size(); i++) {
		PolySeg& clipSeg = clipRoot[i];
		// subject
		std::vector<PolySeg>& root = *mLod[mLod.size() - 1];
		for (int j = 0; j < root.size(); j++) {
			PolySeg& subSeg = root[j];
			_FindCrossSegs(subSeg, clip, clipSeg, mCrossEdges);
		}
	}
	//GetTicker().record("Find crosses");
#if 0
	double elapsed = double(clock() - startTime) * 1000 / CLOCKS_PER_SEC;
	printf("\n\nSubject points: %d\n", mPoly.size());
	printf("Clipper points: %d\n", clip.size());
	printf("Inter    edges: %d\n", mCrossEdges.size());
	printf("Inter    count: %d\n", nInterCounter);
	printf("Speedup       : %.2f\n", double(mPoly.size())*testEdges/nInterCounter);
	printf("elapsed       : %.2f\n", elapsed);
	printf("lod time      : %.2f\n", tlod);
#endif
}


void PolyStruct::Viz(const char* savePath, Polys& extra) {
	// viz
	SvgBase svg;
	svg.style.pft = ClipperLib::pftEvenOdd;

	Polys oriPoly;
	oriPoly.push_back(mPoly);

	svg.AddPath(oriPoly, COLOR_GRAY, 0, true);
	svg.AddPoly(*mClip, COLOR_MY_GREEN, 1, COLOR_MY_GREEN, false);
	//svg.AddPath(extra, COLOR_MY_RED, 0, true);

	//Polys bounds;
	for (int l = 0; l < mLod.size(); l++) {
		break;
		if (l >= mLod.size()-5) {
			continue;
		}
		std::vector<PolySeg>& segs = *mLod[l];
		//printf("Layer: %d size: %d\n", l, segs.size());
		unsigned color = GetColor(l);
		for (size_t i = 0; i < segs.size(); i++)
		{
			PolySeg& seg = segs[i];
			if (seg.type == 0) {
				continue;
			}
			Poly poly(4);
			for (int j = 0; j < 4; j++) {
				poly[j] = Point(seg.rbox[j].x, seg.rbox[j].y);
			}
#if 0
			poly[0] = Point(seg.tl);
			poly[1] = Point(seg.br.x, seg.tl.y);
			poly[2] = Point(seg.br);
			poly[3] = Point(seg.tl.x, seg.br.y);
#endif
			svg.AddPoly(poly, color);
		}
	}

	// cross edges
	Poly points;
	for (size_t i = 0; i < mCrossEdges.size(); i++)
	{
		CrossEdge& edge = mCrossEdges[i];
		Poly poly(2);
		poly[0] = mPoly[edge.sIdx];
		poly[1] = mPoly[(edge.sIdx+1) % N];
		svg.AddPoly(poly, COLOR_MY_RED, 5, COLOR_MY_RED, false);
		points.push_back(edge.p0);
		char txt[256];
		sprintf(txt, "%d", i);
		svg.SetFont("Verdana", 18, 0xFFAAAA00);
		svg.AddText(edge.p0.x, edge.p0.y, std::string(txt));
	}
	/*
	{
		Poly poly(2);
		int sIdx = 89837;
		poly[0] = mPoly[sIdx];
		poly[1] = mPoly[(sIdx + 1) % N];
		svg.AddPoly(poly, COLOR_MY_GREEN, 5, COLOR_MY_GREEN, false);
	}
	{
		Poly poly(2);
		int sIdx = 41000;
		poly[0] = (*mClip)[sIdx];
		poly[1] = (*mClip)[(sIdx + 1) % (*mClip).size()];
		printf("%f, %f\n", poly[0].x, poly[0].y);
		printf("%f, %f\n", poly[1].x, poly[1].y);
		svg.AddPoly(poly, COLOR_MY_BLUE, 5, COLOR_MY_BLUE, false);
		svg.AddPoints(poly, COLOR_MY_BLUE, 5, COLOR_MY_BLUE, false);
	}
	*/
	// end points
	Poly endPoints;
	endPoints.push_back(mPoly[0]);
	endPoints.push_back(mPoly[N-1]);
	svg.AddPoints(endPoints, COLOR_MY_BLUE, 5, COLOR_MY_BLUE, false);

	svg.AddPoints(points, COLOR_MY_GREEN, 5, COLOR_MY_GREEN, false);
	svg.SaveToFile(savePath, 1024, 800);
}