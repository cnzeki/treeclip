#pragma once
#include "types.h"

//#define ARGB(a, r, g, b)  ((a) << 24 | (r) << 16 | (g) << 8 | b)
//#define argb(a, r, g, b) ()
#ifndef RGB
#define RGB(r,g,b)  ((unsigned)( (b) | ((g)<<8) | ((r)<<16) | 255 << 24  ))
#endif

#define COLOR_TRANS  0x00FFFFFF
#define COLOR_GRAY   0xFF808080
#define COLOR_LIGHT  0xFFA0A0A0
#define COLOR_MY_BLUE                 RGB(  0, 0, 255)      // 蓝色    
#define COLOR_MY_BLUE2                RGB(  0, 32, 164)      // 蓝色  
#define COLOR_MY_BLUE3                RGB(  0, 64, 80)      // 蓝色  

#define COLOR_MY_ORANGE               RGB(233,  82,  16)      // 橙色
#define COLOR_MY_YELLOW               RGB(255, 255,  0)      // 黄色
#define COLOR_MY_YELLOW2                RGB( 164, 164,  32)      // 黄绿色
#define COLOR_MY_YELLOW3                RGB( 80, 80,  64)      // 黄绿色

#define COLOR_MY_RED                  RGB(232, 54,   36)      // 红色

#define COLOR_MY_GREEN                RGB( 0, 255,  0)      // 黄绿色
#define COLOR_MY_CYAN                 RGB( 60, 200,  252)     // 青色

static unsigned GetColor(int idx) {
	unsigned colors[] = {
		COLOR_MY_RED,
		COLOR_MY_GREEN,
		COLOR_MY_BLUE,
		COLOR_MY_YELLOW,
		RGB(255, 192, 203),
		COLOR_MY_CYAN,
		COLOR_MY_ORANGE
	};
	return colors[idx % (sizeof(colors) / sizeof(colors[0]))];
}

static unsigned ARGB(unsigned a, unsigned r, unsigned g, unsigned b) {
	return ((a) << 24 | (r) << 16 | (g) << 8 | b);
}

static unsigned argb(float a, float r, float g, float b) {
	return (unsigned(a*255) << 24 | unsigned(r * 255) << 16 | unsigned(g * 255) << 8 | unsigned(b * 255));
}

static unsigned RandColor(float alpha = 1.0) {
	float rgb[3];
	for (int j = 0; j < 3; j++) {
		rgb[j] = float(rand()) / RAND_MAX;
	}
	return argb(1, rgb[0], rgb[1], rgb[2]);
}

static string ColorToHtml(unsigned clr)
{
	stringstream ss;
	ss << '#' << hex << std::setfill('0') << setw(6) << (clr & 0xFFFFFF);
	return ss.str();
}
//------------------------------------------------------------------------------

static float GetAlphaAsFrac(unsigned clr)
{
	return ((float)(clr >> 24) / 255);
}
//------------------------------------------------------------------------------

//a simple class to build SVG files that displays polygons
class SvgBase
{

	struct Rect
	{
		double left;
		double top;
		double right;
		double bottom;
		Rect(double l = 0, double t = 0, double r = 0, double b = 0) :
			left(l), top(t), right(r), bottom(b) {};
	};

	class StyleInfo
	{
	public:
		ClipperLib::PolyFillType pft;
		unsigned brushClr;
		unsigned penClr;
		double penWidth;
		bool closePath;
		bool showCoords;
	};

	struct FontInfo
	{
	public:
		std::string family;
		int size;
		unsigned fillColor;
	};
	typedef std::vector<FontInfo> FontInfoList;

	class TextInfo
	{
	public:
		const std::string text;
		double  x;
		double  y;
		unsigned fontIdx;
		TextInfo(double  _x, double  _y, unsigned _fi, const std::string& _text) : x(_x), y(_y), fontIdx(_fi), text(_text) {}
	};
	typedef std::vector<TextInfo> TextInfoList;

	class PolyInfo
	{
	public:
		Polys polygons;
		StyleInfo si;

		PolyInfo(Polys& p, StyleInfo style)
		{
			this->polygons = p;
			this->si = style;
		}
	};
	typedef std::vector<PolyInfo> PolyInfoList;

	class PointInfo
	{
	public:
		Poly points;
		StyleInfo si;

		PointInfo(Poly& p, StyleInfo style)
		{
			this->points = p;
			this->si = style;
		}
	};
	typedef std::vector<PointInfo> PointInfoList;

private:
	PolyInfoList polyInfos;
	FontInfoList fontInfos;
	TextInfoList textInfos;
	PointInfoList pointInfos;

	static const std::string svg_xml_start[];
	static const std::string path_end_poly[];
	static const std::string path_end_line[];

	void CheckFonts()
	{
		//if no font has been specified create a default font...
		if (!fontInfos.empty()) return;
		FontInfo fi;
		fi.family = "Verdana";
		fi.size = 24;// 15;
		fi.fillColor = 0xFFFF00FF;
		fontInfos.push_back(fi);
	}

	void UpdateBounds(Polys& p)
	{
		Rect r = GetBounds(p);
		if (r.left < bounds.left) bounds.left = r.left;
		if (r.top < bounds.top) bounds.top = r.top;
		if (r.right > bounds.right) bounds.right = r.right;
		if (r.bottom > bounds.bottom) bounds.bottom = r.bottom;
	}
	//---------------------------------------------------------------------------

public:
	StyleInfo style;
	Rect bounds;

	SvgBase()
	{
		style.pft = ClipperLib::pftNonZero;
		style.brushClr = 0xFFFFFFCC;
		style.penClr = 0xFF000000;
		style.penWidth = 0.8;
		style.closePath = true;
		style.showCoords = false;

		bounds.left = DBL_MAX;
		bounds.top = DBL_MAX;
		bounds.right = -DBL_MAX;
		bounds.bottom = -DBL_MAX;
	}

	Rect GetBounds(Polys& p)
	{
		Rect result(DBL_MAX, DBL_MAX, -DBL_MAX, -DBL_MAX);
		for (size_t i = 0; i < p.size(); i++)
			for (size_t j = 0; j < p[i].size(); j++)
			{
				if (p[i][j].x < result.left) result.left = p[i][j].x;
				if (p[i][j].x > result.right) result.right = p[i][j].x;
				if (p[i][j].y < result.top) result.top = p[i][j].y;
				if (p[i][j].y > result.bottom) result.bottom = p[i][j].y;
			}
		return result;
	}

	void AddPath(Polys& poly, unsigned brushClr, unsigned penClr, bool closed)
	{
		if (poly.size() == 0) return;
		CheckFonts();
		
		style.brushClr = brushClr;
		style.penClr = penClr;
		style.closePath = closed;
		PolyInfo pi = PolyInfo(poly, style);
		polyInfos.push_back(pi);
		UpdateBounds(poly);
	}

	void AddPoly(const Poly& poly, unsigned strokeClr, float strokeWidth=0.8, unsigned fillClr=COLOR_TRANS, bool closed=true) {
		if (poly.size() == 0) return;
		CheckFonts();
		StyleInfo style = this->style;
		style.brushClr = fillClr;
		style.penClr = strokeClr;
		style.closePath = closed;
		style.penWidth = strokeWidth;
		//style.showCoords = true;
		Polys polys;
		polys.push_back(poly);
		PolyInfo pi = PolyInfo(polys, style);
		polyInfos.push_back(pi);
		UpdateBounds(polys);
	}

	void AddPoints(const Poly& poly, unsigned strokeClr, float strokeWidth = 0.8, unsigned fillClr = COLOR_TRANS, bool closed = true) {
		if (poly.size() == 0) return;
		StyleInfo style = this->style;
		style.brushClr = fillClr;
		style.penClr = strokeClr;
		style.closePath = closed;
		style.penWidth = strokeWidth;
		//style.showCoords = true;
		Poly p = poly;
		PointInfo pi = PointInfo(p, style);
		pointInfos.push_back(pi);
	}
	//---------------------------------------------------------------------------

	void SetFont(std::string family, int size, unsigned fillColor)
	{
		FontInfo fi;
		fi.family = family;
		fi.size = size;
		fi.fillColor = fillColor;
		fontInfos.push_back(fi);
	}
	//---------------------------------------------------------------------------

	void AddText(double x, double y, const std::string& text)
	{
		CheckFonts();
		TextInfo ti = TextInfo(x, y, fontInfos.size() - 1, text);
		textInfos.push_back(ti);
	}
	//---------------------------------------------------------------------------

	bool SaveToFile(const std::string filename, int width = 0, int height = 0, int margin = 10) const
	{
		if (margin < 0) margin = 0;
		double scale = 1.0;
		if (width > 0 && height > 0)
			scale = 1.0 / max((bounds.right - bounds.left) / width,
			(bounds.bottom - bounds.top) / height);
		ofstream file;
		file.open(filename);
		if (!file.is_open()) return false;
		file.setf(ios::fixed);
		file.precision(0);
		file << svg_xml_start[0] <<
			(int)((bounds.right - bounds.left) *scale + margin * 2) << "px" << svg_xml_start[1] <<
			(int)((bounds.bottom - bounds.top) *scale + margin * 2) << "px" << svg_xml_start[2] <<
			(int)((bounds.right - bounds.left) *scale + margin * 2) << " " <<
			(int)((bounds.bottom - bounds.top) *scale + margin * 2) << svg_xml_start[3];
		setlocale(LC_NUMERIC, "C");
		file.precision(1);

		for (PolyInfoList::size_type i = 0; i < polyInfos.size(); ++i)
		{
			file << " <path d=\"";
			for (ClipperLib::Polygons::size_type j = 0; j < polyInfos[i].polygons.size(); ++j)
			{
				if (polyInfos[i].polygons[j].size() < 2) continue;
				file << " M " << ((double)(polyInfos[i].polygons[j][0].x - bounds.left) * scale + margin) <<
					" " << ((double)(polyInfos[i].polygons[j][0].y - bounds.top) * scale + margin);
				for (ClipperLib::Polygon::size_type k = 1; k < polyInfos[i].polygons[j].size(); ++k)
				{
					const Point ip = polyInfos[i].polygons[j][k];
					double x = (ip.x - bounds.left) * scale + margin;
					double y = (ip.y - bounds.top) * scale + margin;
					file << " L " << x << " " << y;
				}
				if (polyInfos[i].si.closePath) file << " z";
			}
			if (polyInfos[i].si.closePath)
				file <<
				path_end_poly[0] << ColorToHtml(polyInfos[i].si.brushClr) <<
				path_end_poly[1] << GetAlphaAsFrac(polyInfos[i].si.brushClr) <<
				path_end_poly[2] << (polyInfos[i].si.pft == ClipperLib::pftEvenOdd ? "evenodd" : "nonzero") <<
				path_end_poly[3] << ColorToHtml(polyInfos[i].si.penClr) <<
				path_end_poly[4] << GetAlphaAsFrac(polyInfos[i].si.penClr) <<
				path_end_poly[5] << polyInfos[i].si.penWidth <<
				path_end_poly[6];
			else
				file <<
				path_end_line[0] << ColorToHtml(polyInfos[i].si.penClr) <<
				path_end_line[1] << GetAlphaAsFrac(polyInfos[i].si.penClr) <<
				path_end_line[2] << polyInfos[i].si.penWidth <<
				path_end_line[3];

		}
		bool showCoords = false;
		for (size_t i = 0; i < polyInfos.size(); i++)
			if (polyInfos[i].si.showCoords) { showCoords = true; break; }

		if (!textInfos.empty() || showCoords)
		{
			if (showCoords)
			{
				FontInfo fontInfo = fontInfos.front();
				file << "<g font-family=\"" << fontInfo.family << "\" font-size=\"" <<
					(int)ceil(scale * fontInfo.size) << "\" fill=\"" << ColorToHtml(fontInfo.fillColor) << "\">\n";
				for (size_t i = 0; i < polyInfos.size(); i++)
					if (polyInfos[i].si.showCoords)
					{
						for (ClipperLib::Polygons::size_type j = 0; j < polyInfos[i].polygons.size(); ++j)
						{
							if (polyInfos[i].polygons[j].size() < 3) continue;
							for (ClipperLib::Polygon::size_type k = 0; k < polyInfos[i].polygons[j].size(); ++k)
							{
								Point ip = polyInfos[i].polygons[j][k];
								file << "  <text x=\"" << (int)((ip.x - bounds.left) * scale + margin) <<
									"\" y=\"" << (int)((ip.y - bounds.top) * scale + margin) << "\">" <<
									ip.x << ", " << ip.y << "</text>\n";
							}
						}
					}
				if (showCoords)  file << "</g>\n";
			}

			unsigned fi = INT_MAX;
			for (size_t i = 0; i < textInfos.size(); ++i)
			{
				TextInfo ti = textInfos[i];
				//if (ti.fontIdx != fi)
				{
					if (fi != INT_MAX) file << "</g>\n";
					fi = ti.fontIdx;
					FontInfo fontInfo = fontInfos[fi];
					file << "<g font-family=\"" << fontInfo.family << "\" font-size=\"" <<
						(int)ceil(fontInfo.size) << "\" fill=\"" << ColorToHtml(fontInfo.fillColor) << "\">\n";
				}
				file << "  <text x=\"" << (int)((ti.x - bounds.left) * scale + margin) <<
					"\" y=\"" << (int)((ti.y - bounds.top) * scale + margin) << "\">" <<
					ti.text << "</text>\n";
			}
			file << "</g>\n";
		}
		// points
		for (int i = 0; i < pointInfos.size(); ++i)
		{
			for (int j = 0; j < pointInfos[i].points.size(); ++j)
			{
				double x = (double)(pointInfos[i].points[j].x - bounds.left) * scale + margin;
				double y = (double)(pointInfos[i].points[j].y - bounds.top) * scale + margin;
				file << "<circle cx = \"" << x <<"\" ";
				file << "cy = \"" << y << "\" r = \"5\" ";
				file << "stroke = \"" << ColorToHtml(pointInfos[i].si.penClr) << "\" ";
				file << "stroke-width = \"2\" fill = \"red\" /> ";
			}
		}
		file << "</svg>\n";
		file.close();
		setlocale(LC_NUMERIC, "");
		return true;
	}
	//---------------------------------------------------------------------------

};
//------------------------------------------------------------------------------
