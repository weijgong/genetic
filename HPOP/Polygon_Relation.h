//------------------------------------------------------------------------------
//
// Polygon_Clip.h
// 
// Purpose:
//
// 计算区域相交的坐标（经纬度）
//
// 2020.8.20 左义新
//
//------------------------------------------------------------------------------

#ifndef INC_Polygon_CLIP_H
#define INC_Polygon_CLIP_H
#include <vector>
#include <iostream>
#include "SAT_Const.h"
using namespace std;
typedef struct Point
{
	double x;
	double y;
}Point;
bool IsRectCross(const Point &p1, const Point &p2, const Point &q1, const Point &q2);
bool IsLineSegmentCross(const Point &pFirst1, const Point &pFirst2, const Point &pSecond1, const Point &pSecond2);
bool GetCrossPoint(const Point &p1, const Point &p2, const Point &q1, const Point &q2, double &x, double &y);
bool IsPointInpolygon(std::vector<Point> poly, Point pt);
bool PointCmp(const Point &a, const Point &b, const Point &center);
vector<Point> ClockwiseSortPoints(std::vector<Point> &vPoints);
vector<Point> PolygonClip(const vector<Point> &poly1, const vector<Point> &poly2, std::vector<Point> &interPoly);
vector<Point> PolygonUnion(const vector<Point> &poly1, const vector<Point> &poly2, std::vector<Point> &unionPoly);
// vector<vector<Point>> PolygonUnionNew(const vector<Point> &poly1, const vector<Point> &poly2, std::vector<Point> &unionPoly);

double polarTriangleArea(double &tan1, double &lng1, double &tan2, double &lng2);
double CalculateArea(std::vector<Point> &polygon);

double CalculateAreaNew(std::vector<Point> &polygon);
double SphericalPolygonAreaMeters2(std::vector<Point> &polygon);
double Angle(Point &p1, Point &p2, Point &p3);
double Bearing(Point &from, Point &to);
double PlanarPolygonAreaMeters2(std::vector<Point> &polygon);

void DFS(int v, int **x, int Flag[], int StripNum);
double StripTotalCoverage(vector<vector<Point>> StripSet, int StripNum);
#endif  // include blocker