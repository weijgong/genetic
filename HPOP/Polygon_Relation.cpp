//------------------------------------------------------------------------------
//
// Ploygon_Clip.cpp
// 
// Purpose:
//
//    计算区域相交坐标（经纬度）
//
// 2020.8.20 左义新
//
//------------------------------------------------------------------------------

#include <cmath>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <vector>
#include "Polygon_Relation.h"
#include "Graph_Var.h"
#include "GNU_iomanip.h"

using namespace std;

//////////////////////////////////////////////////////////////// 求线段交点
//排斥实验
bool IsRectCross(const Point &p1, const Point &p2, const Point &q1, const Point &q2)
{
	bool ret = min(p1.x, p2.x) <= max(q1.x, q2.x) &&
		min(q1.x, q2.x) <= max(p1.x, p2.x) &&
		min(p1.y, p2.y) <= max(q1.y, q2.y) &&
		min(q1.y, q2.y) <= max(p1.y, p2.y);
	return ret;
}

//跨立判断
bool IsLineSegmentCross(const Point &pFirst1, const Point &pFirst2, const Point &pSecond1, const Point &pSecond2)
{
	double line1, line2;
	line1 = pFirst1.x * (pSecond1.y - pFirst2.y) +
		pFirst2.x * (pFirst1.y - pSecond1.y) +
		pSecond1.x * (pFirst2.y - pFirst1.y);
	line2 = pFirst1.x * (pSecond2.y - pFirst2.y) +
		pFirst2.x * (pFirst1.y - pSecond2.y) +
		pSecond2.x * (pFirst2.y - pFirst1.y);
	if (((line1 >= 0 && line2 >= 0) || (line1 <= 0 && line2 <= 0)) && !(line1 == 0 && line2 == 0))
		return false;

	line1 = pSecond1.x * (pFirst1.y - pSecond2.y) +
		pSecond2.x * (pSecond1.y - pFirst1.y) +
		pFirst1.x * (pSecond2.y - pSecond1.y);
	line2 = pSecond1.x * (pFirst2.y - pSecond2.y) +
		pSecond2.x * (pSecond1.y - pFirst2.y) +
		pFirst2.x * (pSecond2.y - pSecond1.y);
	if (((line1 >= 0 && line2 >= 0) || (line1 <= 0 && line2 <= 0)) && !(line1 == 0 && line2 == 0)) //将异或改成同号判断
		return false;
	return true;
}

// 求线段交点
bool GetCrossPoint(const Point &p1, const Point &p2, const Point &q1, const Point &q2, double &x, double &y)
{
	if (IsRectCross(p1, p2, q1, q2))
	{
		if (IsLineSegmentCross(p1, p2, q1, q2))
		{
			//求交点
			double tmpLeft, tmpRight;
			tmpLeft = (q2.x - q1.x) * (p1.y - p2.y) - (p2.x - p1.x) * (q1.y - q2.y);
			tmpRight = (p1.y - q1.y) * (p2.x - p1.x) * (q2.x - q1.x) + q1.x * (q2.y - q1.y) * (p2.x - p1.x) - p1.x * (p2.y - p1.y) * (q2.x - q1.x);

			x = (double)tmpRight / (double)tmpLeft;

			tmpLeft = (p1.x - p2.x) * (q2.y - q1.y) - (p2.y - p1.y) * (q1.x - q2.x);
			tmpRight = p2.y * (p1.x - p2.x) * (q2.y - q1.y) + (q2.x - p2.x) * (q2.y - q1.y) * (p1.y - p2.y) - q2.y * (q1.x - q2.x) * (p2.y - p1.y);
			y = (double)tmpRight / (double)tmpLeft;
			return true;
		}
	}
	return false;
}
/////////////////////////////////////////////////求线段交点

/////////////////////////////////////////////////判断点是否在多边形内部
//  The function will return YES if the point x,y is inside the polygon, or
//  NO if it is not.  If the point is exactly on the edge of the polygon,
//  then the function may return YES or NO.
bool IsPointInpolygon(std::vector<Point> poly, Point pt)
{
	int i, j;
	bool c = false;
	for (i = 0, j = poly.size() - 1; i < poly.size(); j = i++)
	{
		if ((((poly[i].y <= pt.y) && (pt.y <= poly[j].y)) ||
			((poly[j].y <= pt.y) && (pt.y <= poly[i].y)))
			&& (pt.x <= (poly[j].x - poly[i].x) * (pt.y - poly[i].y) / (poly[j].y - poly[i].y) + poly[i].x))
		{
			c = !c;
		}
	}
	return c;
}
/////////////////////////////////////////////////判断点是否在多边形内部

////////////////////////////////////////////////多边形点集排序
//若点a大于点b,即点a在点b顺时针方向,返回true,否则返回false
bool PointCmp(const Point &a, const Point &b, const Point &center)
{
	if (a.x >= 0 && b.x < 0)
		return true;
	if (a.x == 0 && b.x == 0 && center.x > 0)
		return a.y > b.y;
	if (a.x == 0 && b.x == 0 && center.x <= 0)
		return a.y < b.y;
	//向量OA和向量OB的叉积
	double det = (a.x - center.x) * (b.y - center.y) - (b.x - center.x) * (a.y - center.y);
	if (det >0)
		return true;
	if (det < 0)
		return false;
	//向量OA和向量OB共线，以距离判断大小
	double d1 = (a.x - center.x) * (a.x - center.x) + (a.y - center.y) * (a.y - center.y);
	double d2 = (b.x - center.x) * (b.x - center.y) + (b.y - center.y) * (b.y - center.y);
	return d1 > d2;
}

vector<Point> ClockwiseSortPoints(std::vector<Point> &vPoints)
{
	//计算重心
	Point center;
	double x = 0.0, y = 0.0;
	for (int i = 0; i < vPoints.size(); i++)
	{
		x += vPoints[i].x;
		y += vPoints[i].y;
	}
	center.x = x / vPoints.size();
	center.y = y / vPoints.size();

	//冒泡排序
	for (int i = 0; i < vPoints.size() - 1; i++)
	{
		for (int j = 0; j < vPoints.size() - i - 1; j++)
		{
			if (PointCmp(vPoints[j], vPoints[j + 1], center))
			{
				Point tmp = vPoints[j];
				vPoints[j] = vPoints[j + 1];
				vPoints[j + 1] = tmp;
			}
		}
	}
	// 将相交区域的横坐标精度转换为-180-180
	for (int i = 0; i < vPoints.size(); i++)
	{
		if (vPoints[i].x > 180 && vPoints[i].x <= 360)
		{
			vPoints[i].x = fmod(vPoints[i].x, -360);
		}
	}
	return vPoints;
}
////////////////////////////////////////////////多边形点集排序

///////////////////////////////////////////////求交集算法流程
vector<Point> PolygonClip(const vector<Point> &poly1, const vector<Point> &poly2, std::vector<Point> &interPoly)
{
std:vector<Point>p(1);
	if (poly1.size() < 3 || poly2.size() < 3)
	{
		//return false;

		return p;
	}

	//判断是否存在包含关系
	int flag1 = 1;
	for (int i = 0; i < poly1.size(); i++)
	{
		if (!IsPointInpolygon(poly2, poly1[i])) //多边形1包含于多边形2，返回多边形1
		{
			flag1 = 0;
			break;
		}
	}
	if (flag1 == 1)
	{
		return poly1;
	}

	int flag2 = 1;
	for (int i = 0; i < poly2.size(); i++)
	{
		if (!IsPointInpolygon(poly1, poly2[i])) //多边形2包含于多边形1，返回多边形2
		{
			flag2 = 0;
			break;
		}
	}
	if (flag2 == 1)
	{
		return poly2;
	}

	double x, y;
	//计算多边形交点，交点是交集一部分
	for (int i = 0; i < poly1.size(); i++)
	{
		int poly1_next_idx = (i + 1) % poly1.size();
		for (int j = 0; j < poly2.size(); j++)
		{
			int poly2_next_idx = (j + 1) % poly2.size();
			if (GetCrossPoint(poly1[i], poly1[poly1_next_idx],
				poly2[j], poly2[poly2_next_idx],
				x, y))
			{
				Point point1;
				point1.x = x;
				point1.y = y;;
				interPoly.push_back(point1);
			}
		}
	}

	//计算多边形内部点，内部点是交集一部分
	for (int i = 0; i < poly1.size(); i++)
	{
		if (IsPointInpolygon(poly2, poly1[i]))
		{
			interPoly.push_back(poly1[i]);
		}
	}
	for (int i = 0; i < poly2.size(); i++)
	{
		if (IsPointInpolygon(poly1, poly2[i]))
		{
			interPoly.push_back(poly2[i]);
		}
	}

	if (interPoly.size() <= 0)
		//return false;
		return p;

	//删除重复顶点元素
	std::vector<Point> p4;
	std::vector<Point> p5;
	for (int i = 0; i < interPoly.size(); i++)
	{
		if (interPoly[i].x < 180 && interPoly[i].x > -180)
		{
			p4.push_back(interPoly[i]);
		}
	}
	p5.push_back(p4[0]);
	for (int i = 1; i < p4.size(); i++)
	{
		int flag = 1;
		for (int j = 0; j < p5.size(); j++)
		{
			if (p4[i].x == p5[j].x && p4[i].y == p5[j].y)
			{
				flag = 0;
				break;
			}
		}
		if (flag == 1)
		{
			p5.push_back(p4[i]);
		}
	}
	//点集排序 
	return ClockwiseSortPoints(p5);
}
///////////////////////////////////////////////求交集算法流程

///////////////////////////////////////////////求并集算法流程
vector<Point> PolygonUnion(const vector<Point> &poly1, const vector<Point> &poly2, std::vector<Point> &unionPoly)
{
std:vector<Point>poly;
	//两个都不是多边形
	if (poly1.size() < 3 && poly2.size() < 3)
	{
		//return false;
		Point p;
		p.x = 0.0;
		p.y = 0.0;
		poly.push_back(p);

		return poly;
	}
	//有一个不是多边形
	if (poly1.size() < 3 || poly2.size() < 3)
	{
		if (poly1.size() < 3)
		{
			return poly2;
		}
		else
		{
			return poly1;
		}
	}

	//两个多边形相等
	if (poly1.size() == poly2.size())
	{
		int num = 0;
		for (int i = 0; i < poly1.size(); i++)
		{
			for (int j = 0; j < poly2.size(); j++)
			{
				if (poly1[i].x == poly2[j].x && poly1[i].y == poly2[j].y)
				{
					num = num + 1;
					break;
				}
			}
		}
		if (num == poly1.size())
		{
			return poly1;
		}
	}

	double x, y;
	//计算多边形交点，交点是并集一部分
	for (int i = 0; i < poly1.size(); i++)
	{
		int poly1_next_idx = (i + 1) % poly1.size();
		for (int j = 0; j < poly2.size(); j++)
		{
			int poly2_next_idx = (j + 1) % poly2.size();
			if (GetCrossPoint(poly1[i], poly1[poly1_next_idx],
				poly2[j], poly2[poly2_next_idx],
				x, y))
			{
				Point point1;
				point1.x = x;
				point1.y = y;;
				unionPoly.push_back(point1);
			}
		}
	}

	if (unionPoly.size() <= 0) //两多边形不存在交集,判断是否存在内部点
	{
		//判断是否存在包含关系
		for (int i = 0; i < poly1.size(); i++)
		{
			if (IsPointInpolygon(poly2, poly1[i])) //多边形1包含于多边形2，返回多边形2
			{
				return poly2;
				break;
			}
		}
		for (int i = 0; i < poly2.size(); i++)
		{
			if (IsPointInpolygon(poly1, poly2[i])) //多边形2包含于多边形1，返回多边形1
			{
				return poly1;
				break;
			}
		}
		//不存在内部点时，则多边形相离，返回所有顶点
		for (int i = 0; i < poly1.size(); i++)
		{
			poly.push_back(poly1[i]);
		}
		Point p;
		p.x = 360.0; p.y = 360.0; //隔离两个多边形
		poly.push_back(p);
		for (int i = 0; i < poly2.size(); i++)
		{
			poly.push_back(poly2[i]);
		}
		return poly;
	}

	//多边形存在交集，则计算多边形内部点,不在内部的顶点是并集的一部分
	for (int i = 0; i < poly1.size(); i++)
	{
		if (!IsPointInpolygon(poly2, poly1[i]))
		{
			unionPoly.push_back(poly1[i]);
		}
	}
	for (int i = 0; i < poly2.size(); i++)
	{
		if (!IsPointInpolygon(poly1, poly2[i]))
		{
			unionPoly.push_back(poly2[i]);
		}
	}

	//删除重复顶点元素
	std::vector<Point> p4;
	std::vector<Point> p5;
	for (int i = 0; i < unionPoly.size(); i++)
	{
		if (unionPoly[i].x < 180 && unionPoly[i].x > -180)
		{
			p4.push_back(unionPoly[i]);
		}
	}
	p5.push_back(p4[0]);
	for (int i = 1; i < p4.size(); i++)
	{
		int flag = 1;
		for (int j = 0; j < p5.size(); j++)
		{
			if (p4[i].x == p5[j].x && p4[i].y == p5[j].y || abs(p4[i].x - p5[j].x)<0.1 && abs(p4[i].y - p5[j].y)<0.1) //加入了松弛的参数0.1
			{
				flag = 0;
				break;
			}
		}
		if (flag == 1)
		{
			p5.push_back(p4[i]);
		}
	}
	//点集排序 
	return ClockwiseSortPoints(p5);
}
///////////////////////////////////////////////求并集算法流程

////////////////////////////////////////////根据极角求解球面多边形的面积
double polarTriangleArea(double &tan1, double &lng1, double &tan2, double &lng2)
{
	double deltaLng = lng1 - lng2;
	double t = tan1 * tan2;
	return 2 * atan2(t * sin(deltaLng), 1 + t * cos(deltaLng));
}
double CalculateArea(std::vector<Point> &polygon)
{
	double radius = R_Earth;
	int length = polygon.size();
	if (length < 3)
	{
		return 0.0;
	}
	double total = 0.0;
	Point prev;
	prev = polygon[length - 1];
	double prevTanLat = tan(((pi/ 2 - prev.y * Rad) / 2));
	double prevLng = prev.x * Rad;
	for (int i = 0; i < length; i++)
	{
		double tanLat = tan((pi / 2 - (polygon[i].y) * Rad) / 2);
		double lng = (polygon[i].x) * Rad;
		total = total + polarTriangleArea(tanLat, lng, prevTanLat, prevLng);
		prevTanLat = tanLat;
		prevLng = lng;
	}
	return fabs(total * (radius * radius)) / (1e+6);
}
////////////////////////////////////////////根据极角求解球面多边形的面积

//////////////////////////////////////////////////////////////////根据经纬度计算球面多边形面积
double CalculateAreaNew(std::vector<Point> &polygon)
{
	int Count;
	Count = polygon.size();
	double *PointX = new double[Count];
	double *PointY = new double[Count];
	for (int i = 0; i < Count; i++)
	{
		PointX[i] = polygon[i].x;
		PointY[i] = polygon[i].y;
	}
	string MapUnits = "DEGREES";
	if (Count > 2)
	{
		double mtotalArea = 0;

		if (MapUnits == "DEGREES")//经纬度坐标下的球面多边形
		{
			double LowX = 0.0;
			double LowY = 0.0;
			double MiddleX = 0.0;
			double MiddleY = 0.0;
			double HighX = 0.0;
			double HighY = 0.0;
			double AM = 0.0;
			double BM = 0.0;
			double CM = 0.0;
			double AL = 0.0;
			double BL = 0.0;
			double CL = 0.0;
			double AH = 0.0;
			double BH = 0.0;
			double CH = 0.0;
			double CoefficientL = 0.0;
			double CoefficientH = 0.0;
			double ALtangent = 0.0;
			double BLtangent = 0.0;
			double CLtangent = 0.0;
			double AHtangent = 0.0;
			double BHtangent = 0.0;
			double CHtangent = 0.0;
			double ANormalLine = 0.0;
			double BNormalLine = 0.0;
			double CNormalLine = 0.0;
			double OrientationValue = 0.0;
			double AngleCos = 0.0;
			double Sum1 = 0.0;
			double Sum2 = 0.0;
			double Count2 = 0;
			double Count1 = 0;

			double Sum = 0.0;
			for (int i = 0; i < Count; i++)
			{
				if (i == 0)
				{
					LowX = (double)PointX[Count - 1] * Rad;
					LowY = (double)PointY[Count - 1] * Rad;
					MiddleX = (double)PointX[0] * Rad;
					MiddleY = (double)PointY[0] * Rad;
					HighX = (double)PointX[1] * Rad;
					HighY = (double)PointY[1] * Rad;
				}
				else if (i == Count - 1)
				{
					LowX = (double)PointX[Count - 2] * Rad;
					LowY = (double)PointY[Count - 2] * Rad;
					MiddleX = (double)PointX[Count - 1] * Rad;
					MiddleY = (double)PointY[Count - 1] * Rad;
					HighX = (double)PointX[0] * Rad;
					HighY = (double)PointY[0] * Rad;
				}
				else
				{
					LowX = (double)PointX[i - 1] * Rad;
					LowY = (double)PointY[i - 1] * Rad;
					MiddleX = (double)PointX[i] * Rad;
					MiddleY = (double)PointY[i] * Rad;
					HighX = (double)PointX[i + 1] * Rad;
					HighY = (double)PointY[i + 1] * Rad;
				}
				AM = cos(MiddleY) * cos(MiddleX);
				BM = cos(MiddleY) * sin(MiddleX);
				CM = sin(MiddleY);
				AL = cos(LowY) * cos(LowX);
				BL = cos(LowY) * sin(LowX);
				CL = sin(LowY);
				AH = cos(HighY) * cos(HighX);
				BH = cos(HighY) * sin(HighX);
				CH = sin(HighY);

				CoefficientL = (AM * AM + BM * BM + CM * CM) / (AM * AL + BM * BL + CM * CL);
				CoefficientH = (AM * AM + BM * BM + CM * CM) / (AM * AH + BM * BH + CM * CH);
				ALtangent = CoefficientL * AL - AM;
				BLtangent = CoefficientL * BL - BM;
				CLtangent = CoefficientL * CL - CM;
				AHtangent = CoefficientH * AH - AM;
				BHtangent = CoefficientH * BH - BM;
				CHtangent = CoefficientH * CH - CM;

				AngleCos = (AHtangent * ALtangent + BHtangent * BLtangent + CHtangent * CLtangent) / (sqrt(AHtangent * AHtangent + BHtangent * BHtangent + CHtangent * CHtangent) * sqrt(ALtangent * ALtangent + BLtangent * BLtangent + CLtangent * CLtangent));
				AngleCos = acos(AngleCos);
				ANormalLine = BHtangent * CLtangent - CHtangent * BLtangent;
				BNormalLine = 0 - (AHtangent * CLtangent - CHtangent * ALtangent);
				CNormalLine = AHtangent * BLtangent - BHtangent * ALtangent;
				if (AM != 0)
					OrientationValue = ANormalLine / AM;
				else if (BM != 0)
					OrientationValue = BNormalLine / BM;
				else
					OrientationValue = CNormalLine / CM;
				if (OrientationValue > 0)
				{
					Sum1 += AngleCos;
					Count1++;
				}
				else
				{
					Sum2 += AngleCos;
					Count2++;
					//Sum +=2*M_PI-AngleCos;
				}
			}
			if (Sum1 > Sum2)
			{
				Sum = Sum1 + (2 * pi * Count2 - Sum2);
			}
			else
			{
				Sum = (2 * pi * Count1 - Sum1) + Sum2;
			}
			//平方米
			mtotalArea = (Sum - (Count - 2) * pi) * R_Earth * R_Earth;
		}
		else
		{ //非经纬度坐标下的平面多边形
			int i, j;
			//double j;
			double p1x, p1y;
			double p2x, p2y;
			for (i = Count - 1, j = 0; j < Count; i = j, j++)
			{
				p1x = (double)PointX[i];
				p1y = (double)PointY[i];
				p2x = (double)PointX[j];
				p2y = (double)PointY[j];
				mtotalArea += p1x * p2y - p2x * p1y;
			}
			mtotalArea /= 2.0;
		}
		return mtotalArea / (1e+6); //输出面积平方公里
	}
	return 0;
}
//////////////////////////////////////////////////////////////////根据经纬度计算球面多边形面积

////////////////////////////////////////////深度优先搜索求解子图
void DFS(int v, int **x, int Flag[], int StripNum)
{
	int w, i, k;
	Flag[v] = 1;//把这个顶点的状态设置为1
	//cout << v << " ";
	OrderGraph[PolygonNum++] = v;
	int *temp;
	temp = new int[StripNum];
	for (i = 0; i < StripNum; i++)
		temp[i] = 0;

	k = 0;
	for (i = 0; i < StripNum; i++)//找到与当前这个顶点有连接的顶点，并存入temp数组中
	{
		if (x[v][i] == 1)
			temp[k++] = i;
	}
	i = 0;
	for (i = 0; i < k; i++)//对这个几个顶点递归遍历
	{
		w = temp[i];
		if (Flag[w] == 0)
			DFS(w, x, Flag, StripNum);
	}
}


////////////////////////////////////////////深度优先搜索求解子图




////////////////////////////////////////////根据DFS求解条带覆盖区域的总面积
double StripTotalCoverage(vector<vector<Point>> StripSet, int StripNum)
{
	//初始化N*N邻接关系矩阵
	int **Relation; //动态定义二维数组 
	Relation = new int*[StripNum];
	for (int i = 0; i < StripNum; i++)
		Relation[i] = new int[StripNum];
	//初始化二维数组
	for (int i = 0; i < StripNum; i++)
		for (int j = 0; j < StripNum; j++)
			Relation[i][j] = 0;

	//初始化节点遍历标志位N*N
	int *Flag;
	Flag = new int[StripNum];
	for (int i = 0; i < StripNum; i++)
	{
		Flag[i] = 0;
	}

	//找出区域条带的相交关系，更新Relation数组
	for (int i = 0; i < StripNum - 1; i++)
	{
		for (int j = i + 1; j < StripNum; j++)
		{
			int flag = 1; //相交标志
			vector<Point> un1;
			vector<Point> TempUnion;
			if (StripSet[i].size() >= 3 && StripSet[j].size() >= 3)
			{
				TempUnion = PolygonUnion(StripSet[i], StripSet[j], un1); //通过并集判断是否存在交集可以合并
				for (int n = 0; n < TempUnion.size(); n++)
				{
					if (TempUnion[n].x == 360.0 && TempUnion[n].y == 360.0)  //不存在交集的两个多边形在求解并集时用（360，360）将两多边形的存储隔离
					{
						flag = 0;
						break;
					}
				}
			}
			else //不满足都是多边形，则不存在交集
			{
				flag = 0;
			}
			if (flag == 1)
			{
				Relation[i][j] = 1;
				Relation[j][i] = 1;
			}
		}
	}

	//初始化子图参数
	SonGraph = new int*[StripNum];
	for (int i = 0; i < StripNum; i++)
		SonGraph[i] = new int[StripNum];
	for (int i = 0; i < StripNum; i++)
		for (int j = 0; j < StripNum; j++)
			SonGraph[i][j] = 0;

	GraphAdd = new int[StripNum];
	for (int i = 0; i < StripNum; i++)
		GraphAdd[i] = 0;

	GraphFirst = new int[StripNum];
	for (int i = 0; i < StripNum; i++)
		GraphFirst[i] = 0;

	OrderGraph = new int[StripNum];
	for (int i = 0; i < StripNum; i++)
		OrderGraph[i] = 0;

	SonpointNum = new int[StripNum];
	for (int i = 0; i < StripNum; i++)
		SonpointNum[i] = 0;

	PolygonNum = 0;
	GraphNum = 0;
	//找出并集原子团,二维数组表示
	for (int i = 0; i < StripNum; i++)//可能有多个连通图
	{
		if (Flag[i] == 0)//找到一个状态为0的顶点
		{
			GraphFirst[GraphNum++] = i;
			DFS(i, Relation, Flag, StripNum);
		}
	}
	// 将序号分成二维数组存储
	for (int i = 0; i < GraphNum; i++)
	{
		for (int j = 0; j < PolygonNum; j++)
		{
			if (OrderGraph[j] == GraphFirst[i])
			{
				GraphAdd[i] = j;
			}
		}
	}
	GraphAdd[GraphNum++] = PolygonNum; //最后一个图的结尾下标定义为PolygonNum

	//为二位数组赋值，子图的边形序号
	for (int i = 0; i < GraphNum - 1; i++)
	{
		int j = 0;
		//cout << "子图" << i << endl;
		for (int k = 0; k < PolygonNum; k++)
		{
			if (k < GraphAdd[i + 1] && k >= GraphAdd[i])
			{
				SonGraph[i][j++] = OrderGraph[k]; //子图个数
				//cout << SonGraph[i][j - 1] << endl;
			}
		}
		SonpointNum[i] = GraphAdd[i + 1] - GraphAdd[i]; //各子图多边形个数
	}

	//根据SonGraph中元素和SonpointNum数值计算多边形并集面积
	vector<Point> *TempUnion2;
	TempUnion2 = new vector<Point>[StripNum];
	double *SonGraph_Area;
	SonGraph_Area = new double[StripNum];
	for (int i = 0; i < GraphNum - 1; i++)
	{
		TempUnion2[i] = StripSet[SonGraph[i][0]];
		for (int j = 0; j < SonpointNum[i]; j++)
		{
			vector<Point> un2;
			TempUnion2[i] = PolygonUnion(TempUnion2[i], StripSet[SonGraph[i][j]], un2); //求各子图的并集
		}
		SonGraph_Area[i] = CalculateAreaNew(TempUnion2[i]); //求子图面积
	}

	//计算总的覆盖面积
	double Total_Area = 0.0;
	for (int i = 0; i < GraphNum - 1; i++)
	{
		Total_Area = Total_Area + SonGraph_Area[i];
	}
	return Total_Area;
}
