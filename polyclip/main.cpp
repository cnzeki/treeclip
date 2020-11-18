#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif

#ifdef _DEBUG
//#include <vld.h> //leak detection 
#endif

#include <cmath>
#include <ctime>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <numeric>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cfloat>
#include <exception>
#include <direct.h>
#include <io.h>

#include "clipper/clipper.hpp"          //http://sourceforge.net/projects/polyclipping/ (clipper)
#include "greiner/greiner.h"
#include "gpc/gpc.h"
#include "treeclip/treeclip.h"
//---------------------------------------------------------------------------
#include "svg.h"
#include "common/ticker.h"

using namespace std;

const double INT_SCALE = 1000;


typedef std::vector< Point > Poly;
typedef std::vector< Poly > Polys;


//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

void SimpleSVG(const string filename, Polys& subj, Polys& clip, Polys& solution,
	int width = 1024, int height = 800)
{
	SvgBase svg;
	svg.style.pft = ClipperLib::pftEvenOdd;
	svg.AddPath(subj, 0x206666AC, 0xCCD0D0DD, true);
	svg.AddPath(clip, 0x24666600, 0xCCDDDD80, true);
	svg.AddPath(solution, 0xFF99FF99, 0x40009900, true);
	//svg.SetFont("Verdana", 16, 0xFF0000AA);
	//svg.AddText(svg.bounds.left, svg.bounds.top + 20, "Clipper Benchmarks");
	//svg.SetFont("Verdana", 12, 0xFFAA0000);
	//svg.AddText(svg.bounds.left, svg.bounds.top + 36, "&#0169; Angus Johnson 2013");
	svg.SaveToFile(filename, width, height);
}
//------------------------------------------------------------------------------

inline ClipperLib::long64 Round(double val)
{
	if ((val < 0)) return (ClipperLib::long64)(val - 0.5); else return (ClipperLib::long64)(val + 0.5);
}
//------------------------------------------------------------------------------

double Area(Poly& poly)
{
	int highI = poly.size() - 1;
	if (highI < 2) return 0;
	double a;
	a = (poly[highI].x + poly[0].x) * (poly[0].y - poly[highI].y);
	for (int i = 1; i <= highI; ++i)
		a += (poly[i - 1].x + poly[i].x) * (poly[i].y - poly[i - 1].y);
	return a / 2;
}
//------------------------------------------------------------------------------

void Ellipse(Poly& p, double cx, double cy, double rx, double ry, int steps = 0)
{
	const double pi = 3.1415926535898, tolerance = 0.125;
	double r = (rx + ry) / 2;
	if (steps <= 0) steps = int(pi / (acos(1 - tolerance / r)));
	if (steps < 3) steps = 3;
	p.resize(steps);
	double sn = sin(2 * pi / steps);
	double cs = cos(2 * pi / steps);
	Point pt = Point(1.0, 0.);
	for (int i = 0; i < steps; i++)
	{
		p[i].x = cx + pt.x * rx;
		p[i].y = cy + pt.y * ry;
		//cross & dot products avoids repeat calls to sin() & cos() 
		pt = Point(pt.x * cs - sn * pt.y, pt.x * sn + pt.y * cs);
	}
}
//------------------------------------------------------------------------------

void Star(Poly& p, double cx, double cy, double radius1, double radius2,
	int count = 5, double offset_angle = 0.0)
{
	const double pi = 3.1415926535898;
	if (count <= 5) count = 5;
	count *= 2;
	p.resize(count);
	double sn = sin(2 * pi / count);
	double cs = cos(2 * pi / count);
	Point delta(1.0, 0.0);
	if (offset_angle != 0.0)
	{
		delta.x = cos(offset_angle / count);
		delta.y = sin(offset_angle / count);
	}
	for (int i = 0; i < count; i++)
	{
		double r = (i % 2 == 0 ? radius1 : radius2);
		p[i].x = cx + delta.x * r;
		p[i].y = cy + delta.y * r;
		//cross & dot products faster than repeat calls to sin() & cos() ...
		delta = Point(delta.x * cs - sn * delta.y, delta.x * sn + delta.y * cs);
	}
}
//------------------------------------------------------------------------------

void MakeRandomPoly(Poly& poly, int width, int height, unsigned vertCnt)
{
	//stress_factor > 1 causes more frequent complex intersections with GPC crashes
	const int stress_factor = 10;
	//make vertices a multiple of stress_factor ...
	poly.resize(vertCnt);
	int w = width / stress_factor, h = height / stress_factor;
	for (unsigned i = 0; i < vertCnt; ++i)
	{
		poly[i].x = (rand() % w) * stress_factor;
		poly[i].y = (rand() % h) * stress_factor;
	}
}
//---------------------------------------------------------------------------

bool LoadFromWlrFile(char *filename, Polys &pp) {
	FILE *f = fopen(filename, "r");
	if (!f) return false;
	int polyCnt, vertCnt, i, j;
	int X, Y;
	pp.clear();
	if ((fscanf(f, "%d", &i) == 1) &&
		(fscanf(f, "%d", &polyCnt) == 1) && (i == 1) && (polyCnt > 0)) {
		pp.resize(polyCnt);
		for (i = 0; i < polyCnt; i++) {
			if (fscanf(f, "%d", &vertCnt) != 1) break;
			pp[i].resize(vertCnt);
			for (j = 0; j < vertCnt; j++) {
				fscanf(f, "%d,%d", &X, &Y);
				pp[i][j].x = X; pp[i][j].y = Y;
			}
		}
	}
	fclose(f);
	return true;
}
//---------------------------------------------------------------------------


//---------------------------------------------------------------------------

void LoadClipper(ClipperLib::Polygons &p, Polys& polys)
{
	p.resize(polys.size());
	for (size_t i = 0; i < polys.size(); i++)
	{
		p[i].resize(polys[i].size());
		for (size_t j = 0; j < polys[i].size(); j++)
		{
			p[i][j].X = Round(polys[i][j].x *INT_SCALE);
			p[i][j].Y = Round(polys[i][j].y *INT_SCALE);
		};
	}
}
//---------------------------------------------------------------------------

void LoadGPC(gpc_polygon &p, Polys& polys)
{
	gpc_clear_polygon(p);
	p.num_contours = polys.size();
	p.hole = new int[p.num_contours];
	p.contour = new gpc_vertex_list[p.num_contours];
	for (int i = 0; i < p.num_contours; i++)
	{
		p.hole[i] = (Area(polys[i]) >= 0 ? 0 : 1);
		size_t s = polys[i].size();
		p.contour[i].num_vertices = s;
		p.contour[i].vertex = new gpc_vertex[s];
		for (size_t j = 0; j < s; j++)
		{
			p.contour[i].vertex[j].x = polys[i][j].x;
			p.contour[i].vertex[j].y = polys[i][j].y;
		}
	}
}
//---------------------------------------------------------------------------

void LoadGreiner(Polygon &polygon, Polys &polys)
{
	int ncontours;
	double px, py;

	ncontours = polys.size();
	for (int i = 0; i < ncontours; i++) {
		int npoints, level;
		npoints = polys[i].size();
		level = 1;
		Contour& contour = polygon.pushbackContour();
		for (int j = 0; j < npoints; j++) {
			px = polys[i][j].x;
			py = polys[i][j].y;
			if (j > 0 && px == contour.vertex(j - 1).x && py == contour.vertex(j - 1).y)
				continue;
			if (j == npoints - 1 && px == contour.vertex(0).x && py == contour.vertex(0).y)
				continue;
			contour.add(Point(px, py));
		}
		if (contour.nvertices() < 3) {
			polygon.deletebackContour();
			//printf("delete contour %d\n", i);
			continue;
		}
	}
}
//--------------------------------------------------------------------------

void LoadTreeClip(TreeClip::Polygons &p, Polys& polys)
{
	p.resize(polys.size());
	for (size_t i = 0; i < polys.size(); i++)
	{
		p[i].resize(polys[i].size());
		for (size_t j = 0; j < polys[i].size(); j++)
		{
			p[i][j].x = (float)polys[i][j].x;
			p[i][j].y = (float)polys[i][j].y;
		};
	}
}
//--------------------------------------------------------------------------

void UnloadClipper(Polys& polys, ClipperLib::Polygons &p)
{
	polys.resize(p.size());
	for (size_t i = 0; i < p.size(); i++)
	{
		polys[i].resize(p[i].size());
		for (size_t j = 0; j < p[i].size(); j++)
		{
			polys[i][j].x = (double)p[i][j].X / INT_SCALE;
			polys[i][j].y = (double)p[i][j].Y / INT_SCALE;
		};
	}
}
//---------------------------------------------------------------------------

void UnloadGPC(Polys& polys, gpc_polygon &p)
{
	polys.resize(p.num_contours);
	for (int i = 0; i < p.num_contours; i++)
	{
		gpc_vertex_list vs = p.contour[i];
		polys[i].resize(vs.num_vertices);
		for (int j = 0; j < vs.num_vertices; j++)
		{
			polys[i][j].x = vs.vertex[j].x;
			polys[i][j].y = vs.vertex[j].y;
		};
	}
}
//---------------------------------------------------------------------------

void UnloadGreiner(Polys &polys, Polygon &polygon)
{
	size_t ncontours = polygon.ncontours();
	polys.resize(ncontours);
	for (size_t i = 0; i < ncontours; i++)
	{
		Contour &contour = polygon.contour(i);
		size_t npoints = contour.nvertices();
		polys[i].resize(npoints);
		for (size_t j = 0; j < npoints; j++)
		{
			polys[i][j] = contour.vertex(j);
		};
	}
}
//-------------------------------------------------------------------

void UnloadTreeClip(Polys& polys, TreeClip::Polygons &p)
{
	polys.resize(p.size());
	for (size_t i = 0; i < p.size(); i++)
	{
		polys[i].resize(p[i].size());
		for (size_t j = 0; j < p[i].size(); j++)
		{
			polys[i][j].x = (double)p[i][j].x;
			polys[i][j].y = (double)p[i][j].y;
		};
	}
}
//-------------------------------------------------------------------

int CountVertices(Polys& p)
{
	int cnt = 0;
	for (size_t i = 0; i < p.size(); i++)
		for (size_t j = 0; j < p[i].size(); j++)
			cnt++;
	return cnt;
}
//---------------------------------------------------------------------------

double DoClipper(Polys& subj, Polys& clip, Polys& solution, BoolType bt = Intersection)
{
	ClipperLib::Polygons clipper_subj, clipper_clip, clipper_solution;
	LoadClipper(clipper_subj, subj);
	LoadClipper(clipper_clip, clip);

	GetTicker().reset();
	GetTicker().tick("start");

	ClipperLib::ClipType op = ClipperLib::ctIntersection;
	switch (bt)
	{
	case Union: op = ClipperLib::ctUnion; break;
	case Difference: op = ClipperLib::ctDifference; break;
	case Xor: op = ClipperLib::ctXor; break;
	default: op = ClipperLib::ctIntersection; break;
	}

	ClipperLib::Clipper cp;
	cp.AddPolygons(clipper_subj, ClipperLib::ptSubject);
	cp.AddPolygons(clipper_clip, ClipperLib::ptClip);

	double elapsed = 0.;
	if (cp.Execute(op, clipper_solution, ClipperLib::pftEvenOdd, ClipperLib::pftEvenOdd))
	{
		GetTicker().record("Clipper", "start");
		elapsed = GetTicker().getTime("Clipper");
		GetTicker().report("Clipper");
	}

	UnloadClipper(solution, clipper_solution);

	return elapsed;
}
//---------------------------------------------------------------------------

double DoGPC(Polys& subj, Polys& clip, Polys& solution, BoolType bt = Intersection)
{
	gpc_polygon gpc_subj, gpc_clip, gpc_solution;
	gpc_subj.num_contours = 0;
	gpc_clip.num_contours = 0;
	gpc_solution.num_contours = 0;
	LoadGPC(gpc_subj, subj);
	LoadGPC(gpc_clip, clip);

	GetTicker().reset();
	GetTicker().tick("start");

	gpc_op op = GPC_INT;
	switch (bt)
	{
	case Union: op = GPC_UNION; break;
	case Difference: op = GPC_DIFF; break;
	case Xor: op = GPC_XOR; break;
	default: op = GPC_INT; break;
	}
	double elapsed = 0;

#ifdef _WIN32
	//This Windows specific function adds Structured Exception Handling  
	//and prevents *most* GPC crashes ...
	try
	{
		gpc_polygon_clip(op, &gpc_subj, &gpc_clip, &gpc_solution);
	}
	catch (std::exception& e)
	{
		cout << "Standard exception: " << e.what() << endl;
}
#else
	gpc_polygon_clip(GPC_INT, subj, clip, sol);
#endif

	GetTicker().record("GPC", "start");
	elapsed = GetTicker().getTime("GPC");
	GetTicker().report("GPC");

	if (elapsed != elapsed) elapsed = 0; //test for NAN

	UnloadGPC(solution, gpc_solution);
	gpc_clear_polygon(gpc_subj);
	gpc_clear_polygon(gpc_clip);
	gpc_clear_polygon(gpc_solution);
	return elapsed;
}
//---------------------------------------------------------------------------

double DoGreiner(Polys& subj, Polys& clip, Polys& solution, BoolType bt = Intersection)
{
	Polygon greiner_subj, greiner_clip, greiner_solution;
	int GreinerResult;

	LoadGreiner(greiner_subj, subj);
	LoadGreiner(greiner_clip, clip);
	
	double elapsed = 0;
	GetTicker().reset();
	GetTicker().tick("start");

	BoolOpType op = INTERSECTION;
	switch (bt)
	{
	case Union: op = UNION; break;
	case Difference: op = DIFF; break;
	case Xor: op = XOR; break;
	default: op = INTERSECTION; break;
	}

	GreinerHormann gh(greiner_subj, greiner_clip);
	GreinerResult = gh.boolop(op, greiner_solution);

	GetTicker().record("Greiner", "start");
	elapsed = GetTicker().getTime("Greiner");
	GetTicker().report("Greiner");

	UnloadGreiner(solution, greiner_solution);

	return elapsed;
}

//-----------------------------------------------------------------------------
bool sortInSubj(const TreeClip::IntersectionPoint &p1, const TreeClip::IntersectionPoint &p2)
{
	return p1.step0 < p2.step0;//ÉýÐòÅÅÁÐ  
}

bool sortInClip(const TreeClip::IntersectionPoint &p1, const TreeClip::IntersectionPoint &p2)
{
	return p1.step1 < p2.step1;//ÉýÐòÅÅÁÐ  
}

vector<size_t> sort_indexes(TreeClip::IntersectionPoints &points, bool isSubj)
{

	// initialize original index locations
	vector<size_t> idx(points.size());
	std::iota(idx.begin(), idx.end(), 0);

	// sort indexes based on comparing values in points
	// using std::stable_sort instead of std::sort
	// to avoid unnecessary index re-orderings
	// when points contains elements of equal values 
	if (isSubj)
	{
		stable_sort(idx.begin(), idx.end(),
			[&points](size_t i1, size_t i2) {return points[i1].less0(points[i2]); });
	}
	else
	{
		stable_sort(idx.begin(), idx.end(),
			[&points](size_t i1, size_t i2) {return points[i1].less1(points[i2]); });
	}

	return idx;
}

double DoTreeClip(Polys& subj, Polys& clip, Polys& solution, BoolType bt = Intersection)
{
	TreeClip::Polygons treeclip_subj, treeclip_clip;

	LoadTreeClip(treeclip_subj, subj);
	LoadTreeClip(treeclip_clip, clip);

	double elapsed = 0;

	GetTicker().reset();
	GetTicker().tick("start");

	TreeClip::Polygon &subj0 = treeclip_subj[0];
	TreeClip::Polygon &clip0 = treeclip_clip[0];
	
	PolyStruct ps(clip[0], 1, 16, -1);
	ps.BuildLod();

	GetTicker().record("Build subject LOD");

	ps.FindCrossEdgesLOD(subj[0]);

	GetTicker().record("stage1", "start");

	size_t k, i = 0, j = 0;
	TreeClip::IntersectionPoints subj_inter_points;
	for (k = 0; k < ps.mCrossEdges.size(); ++k)
	{
		CrossEdge& edge = ps.mCrossEdges[k];

		if (edge.numPts == 1)
		{
			TreeClip::IntersectionPoint interPoint(edge.p0.x, edge.p0.y, edge.cIdx, edge.sIdx);

			TreeClip::Point startPoint1(ps.mPoly[edge.sIdx].x, ps.mPoly[edge.sIdx].y), endPoint1(ps.mPoly[(edge.sIdx + 1) % ps.mPoly.size()].x, ps.mPoly[(edge.sIdx + 1) % ps.mPoly.size()].y);
			TreeClip::Point startPoint0(ps.mClip[0][edge.cIdx].x, ps.mClip[0][edge.cIdx].y), endPoint0(ps.mClip[0][(edge.cIdx + 1) % ps.mClip->size()].x, ps.mClip[0][(edge.cIdx + 1) % ps.mClip->size()].y);

			double _a1 = endPoint0.y - startPoint0.y;
			double _b1 = startPoint0.x - endPoint0.x;
			double a2 = endPoint1.y - startPoint1.y;
			double b2 = startPoint1.x - endPoint1.x;
			double determinant = _a1 * b2 - a2 * _b1;

			if (fabs(determinant) < 0.001)
				continue;

			interPoint.dist0 = interPoint.calcDist(startPoint0, endPoint0);
			interPoint.step0 = interPoint.ind0 + interPoint.dist0;
			interPoint.dist1 = interPoint.calcDist(startPoint1, endPoint1);
			interPoint.step1 = interPoint.ind1 + interPoint.dist1;

			if (fabs(startPoint0.x - interPoint.x) < 0.01 && fabs(startPoint0.y - interPoint.y) < 0.01)
			{
				interPoint.flag0 = 1;
				//printf("%d flag0 is 1\n", k);
			}
			else
			{
				if (fabs(endPoint0.x - interPoint.x) < 0.01 && fabs(endPoint0.y - interPoint.y) < 0.01)
				{
					interPoint.flag0 = 2;
				}
			}

			if (fabs(startPoint1.x - interPoint.x) < 0.01 && fabs(startPoint1.y - interPoint.y) < 0.01)
			{
				interPoint.flag1 = 1;
				//printf("%d flag1 is 1\n", k);
			}
			else
			{
				if (fabs(endPoint1.x - interPoint.x) < 0.01 && fabs(endPoint1.y - interPoint.y) < 0.01)
				{
					interPoint.flag1 = 2;
				}
			}

			if (interPoint.flag0 != 1 && interPoint.flag1 != 1)
			{
				interPoint.entering = determinant < 0;

				subj_inter_points.push_back((interPoint));
			}
		}
	}
	
	TreeClip::IntersectionPoints clip_inter_points(subj_inter_points);
	sort(clip_inter_points.begin(), clip_inter_points.end(), sortInClip);
	sort(subj_inter_points.begin(), subj_inter_points.end(), sortInSubj);

	TreeClip::WeilerAtherton weiler;
	TreeClip::Polygons result_solution = weiler.clipPoints(subj0, clip0, subj_inter_points, clip_inter_points);

	GetTicker().record("TreeClip", "start");
	elapsed = GetTicker().getTime("TreeClip");
	GetTicker().report("TreeClip");

	UnloadTreeClip(solution, result_solution);

	return elapsed;
}
//--------------------------------------------------------------------------

int MakeShrinkingEllipses(Polys& p, int count, Point& center, Point& radius, double step)
{
	p.resize(count);
	int result = 0;
	for (int i = 0; i < count; i++)
	{
		if (i*step + 1 >= radius.x || i*step + 1 >= radius.y)
		{
			p.resize(i);
			break;
		}
		Ellipse(p[i], center.x, center.y, radius.x - (i*step), radius.y - (i*step), 0);
		result += p[i].size();
		if (i % 2 != 0) reverse(p[i].begin(), p[i].end());
	}
	return result;
}
//---------------------------------------------------------------------------

int MakeShrinkingRects(Polys& p, int count, Point& center, Point& radius, double step)
{
	p.resize(count);
	for (int i = 0; i < count; i++)
	{
		if (i*step + 1 >= radius.x || i*step + 1 >= radius.y) break;
		p[i].resize(4);
		p[i][0] = Point(center.x - radius.x + (i*step), center.y - radius.y + (i*step));
		p[i][1] = Point(center.x + radius.x - (i*step), center.y - radius.y + (i*step));
		p[i][2] = Point(center.x + radius.x - (i*step), center.y + radius.y - (i*step));
		p[i][3] = Point(center.x - radius.x + (i*step), center.y + radius.y - (i*step));
		if (i % 2 != 0) reverse(p[i].begin(), p[i].end());
	}
	return count * 4;
}
//---------------------------------------------------------------------------

int MakeFanBlades(Polys& p, int blade_cnt, Point& center, Point& radius)
{
	const int inner_rad = 60;
	blade_cnt *= 2;
	if (blade_cnt < 8) blade_cnt = 8;
	if (radius.x < inner_rad + 10) radius.x = inner_rad + 10;
	if (radius.y < inner_rad + 10) radius.y = inner_rad + 10;
	p.resize(1);
	Poly &pg = p.back(), inner, outer;
	Ellipse(outer, center.x, center.y, radius.x, radius.y, blade_cnt);
	Ellipse(inner, center.x, center.y, inner_rad, inner_rad, blade_cnt);
	pg.resize(blade_cnt * 2);
	for (int i = 0; i + 1 < blade_cnt; i += 2)
	{
		pg[i * 2] = inner[i];
		pg[i * 2 + 1] = outer[i];
		pg[i * 2 + 2] = outer[i + 1];
		pg[i * 2 + 3] = inner[i + 1];
	}
	return blade_cnt * 2;
}
//---------------------------------------------------------------------------

void Do(Polys& subj, Polys& clip, Polys& sol)
{
	int i, j;

	while (subj.size() > 1)
	{
		subj.pop_back();
	}

	while (clip.size() > 1)
	{
		clip.pop_back();
	}

	cout << "No. vertices in subject & clip polygons: " << CountVertices(subj) << '*' << CountVertices(clip) << '\n';

	std::vector<std::pair<std::string, double>> elapsed_arr(5);

	int iterCount = 1;
	for (i = 0; i < iterCount; ++i)
	{
		j = 0;

		//elapsed_arr[j].first = "GPC";
		//elapsed_arr[j++].second += DoGPC(subj, clip, sol);

		elapsed_arr[j].first = "Clipper";
		elapsed_arr[j++].second += DoClipper(subj, clip, sol);

		elapsed_arr[j].first = "Greiner";
		elapsed_arr[j++].second += DoGreiner(subj, clip, sol);

		elapsed_arr[j].first = "TreeClip";
		elapsed_arr[j++].second += DoTreeClip(subj, clip, sol);
	}

	for (i = 0; i < j; ++i)
	{
		if(elapsed_arr[i].second > 0.001)
			cout << elapsed_arr[i].first << " elapsed " << elapsed_arr[i].second / iterCount << "milliseconds" << endl;
	}
}

void StarTest()
{
	Polys subj(1), clip(1), sol;

	cout << "Star Test:" << endl;

	Point center1 = Point(310, 320);
	Star(subj[0], 325, 325, 300, 150, 250, 0.0);
	Star(clip[0], 325, 325, 300, 150, 250, 0.005);

	Do(subj, clip, sol);

	//clip.clear();
	//subj.clear();
	//sol.clear();

	SimpleSVG("output/st_stars.svg", subj, clip, sol, 0, 0);
	cout << "Test finished. ('output/st_classic.svg' file created)\n\n";
}
//---------------------------------------------------------------------------

void ClassicTest()
{
	Polys subj, clip, sol;

	cout << "\nClassic Test:\n" << endl;
	if (!LoadFromWlrFile("data/s.wlr", subj) || !LoadFromWlrFile("data/c.wlr", clip))
	{
		cout << "\nUnable to find or load 's.wlr' or 'c.wlr'.\n";
		cout << "Aborting test.\n";
		return;
	}

	Do(subj, clip, sol);

	SimpleSVG("output/st_classic.svg", subj, clip, sol, 600, 600); //can do this after any of the above
	//
	cout << "Test finished. ('output/st_classic.svg' file created)\n\n";
}
//---------------------------------------------------------------------------

void EllipseAndFanTest()
{
	Polys subj, clip, sol;

	cout << "Ellipses and Fan Test:" << endl;

	Point center1 = Point(310, 320), center2 = Point(410, 350);
	MakeShrinkingEllipses(subj, 80, center1, Point(290, 320), 5);
	MakeFanBlades(clip, 64, center2, Point(340, 300));

	Do(subj, clip, sol);

	SimpleSVG("output/st_ellipse_fan.svg", subj, clip, sol, 0, 0); //can do this after any of the above
	cout << "Test finished. ('output/st_ellipse_fan.svg' file created)\n\n";
}
//---------------------------------------------------------------------------

void EllipseAndRectTest()
{
	Polys subj, clip, sol;

	cout << "Ellipses and Rectangles Test:" << endl;

	Point center1 = Point(310, 320), center2 = Point(410, 350);
	MakeShrinkingEllipses(subj, 80, center1, Point(290, 320), 5);
	MakeShrinkingRects(clip, 80, center2, Point(340, 300), 5);

	Do(subj, clip, sol);

	SimpleSVG("output/st_ellipse_rect.svg", subj, clip, sol, 0, 0); //can do this after any of the above
	cout << "Test finished. ('output/st_ellipse_rect.svg' file created)\n\n";
}
//---------------------------------------------------------------------------

int main(int argc, char* argv[])
{
	char * dir_name = "output";
	if (access(dir_name, 0) != 0)
	{
		if (mkdir(dir_name) != 0)
		{
			printf("failed to create directory.\n");
			return -1;
		}
	}
	 

	EllipseAndFanTest();
	EllipseAndRectTest();
	StarTest();

	ClassicTest();
	return 0;
}
