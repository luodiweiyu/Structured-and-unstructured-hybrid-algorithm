#include<iostream>
#include"include/const.h"
#include"include/functions.h"
using namespace std;
using ConstPara::pi;
double distance(mesh& a, mesh& b)
{
	double dx = a.x - b.x;
	double dy = a.y - b.y;
	return sqrt(dx * dx + dy * dy);
}
double distance(double x1, double y1, double x2, double y2)
{
	double dx = x1 - x2;
	double dy = y1 - y2;
	return sqrt(dx * dx + dy * dy);
}

double get_theta(double x1, double y1, double x2, double y2)//求直线与x轴的夹角
{
	double theta;
	if (x1 == x2)
	{
		//if (y2 < y1)
		//	theta = -pi / 2;
		//else
		//	theta = pi / 2;
		theta = pi / 2;
	}
	else if (y1 == y2)
		theta = 0;
	else
	{
		theta = atan(abs((y2 - y1) / (x2 - x1)));
		if ((x2 > x1&& y2 < y1) || (x2 < x1 && y2 > y1))
			theta = -theta;
	}
	return theta;
}
double absmax(double a, double b)
{
	if (abs(a) >= abs(b))
		return a;
	else
		return b;
}
double absmin(double a, double b)
{
	if (abs(a) <= abs(b))
		return a;
	else
		return b;

}
double max(double a, double b)
{
	if (a >= b)
		return a;
	else
		return b;
}
double min(double a, double b)
{
	if (a <= b)
		return a;
	else
		return b;

}
double area(mesh A, mesh B, mesh C, mesh D)//求任意四点构成四边形面积 
{
	return 0.5 * abs(A.x * B.y + B.x * C.y + C.x * D.y + D.x * A.y - B.x * A.y - C.x * B.y - D.x * C.y - A.x * D.y);
}
Line getLine(mesh A, mesh B)
{
	Line L;
	if (A.x == B.x)
	{
		if (A.y == B.y)
		{
			std::cout << "something wrong in getLine !" << std::endl;
			L.A = 0, L.B = 0, L.C = 0;
		}
		else
			L.A = 1, L.B = 0, L.C = -A.x;
	}
	else if (A.y == B.y)
	{
		L.A = 0, L.B = 1, L.C = -A.y;
	}
	else
		L.A = (B.y - A.y) / (B.x - A.x), L.B = -1, L.C = (B.x * A.y - A.x * B.y) / (B.x - A.x);
	return L;
}
Line getLine(double x1, double y1, double x2, double y2)
{
	Line L;
	double delta = 1e-10;
	if (abs(x1 - x2) < delta)
	{
		if (abs(y1 - y2) < delta)
		{
			std::cout << "something wrong in getLine !" << std::endl;
			L.A = 0, L.B = 0, L.C = 0;
		}
		else
			L.A = 1, L.B = 0, L.C = -x1;
	}
	else if (abs(y1 - y2) < delta)
	{
		L.A = 0, L.B = 1, L.C = -y1;
	}
	else
		L.A = (y2 - y1) / (x2 - x1), L.B = -1, L.C = (x2 * y1 - x1 * y2) / (x2 - x1);
	return L;
}
Line getLine(double theta, mesh A)
{
	double k = tan(theta);
	double b = A.y - k * A.x;
	Line L;
	L.A = k, L.B = -1, L.C = b;
	return L;
}
bool judgeFieldInOut(mesh& A, vector <mesh>& poly)
//judge whether the point is inside the polygon
//html: https://blog.csdn.net/u011722133/article/details/52813374 
{
	float maxX, maxY, minX, minY;
	maxX = minX = poly[0].x;
	maxY = minY = poly[0].y;
	int i, j, c = 0;
	for (int i = 0; i < poly.size(); i++)
	{
		maxX = max(maxX, poly[i].x);
		maxY = max(maxY, poly[i].y);
		minX = min(minX, poly[i].x);
		minY = min(minY, poly[i].y);
	}
	if (A.x<minX || A.x>maxX || A.y<minY || A.y>maxY)
		c = 0;
	else
		for (i = 0, j = poly.size() - 1; i < poly.size(); j = i++)
		{
			//if (poly[j].y == poly[i].y && (poly[i].y > A.y) != (poly[j].y > A.y))
			//	return !c;
			/*else*/ if ((poly[i].y > A.y) != (poly[j].y > A.y) && (A.x < (poly[j].x - poly[i].x) * (A.y - poly[i].y) / (poly[j].y - poly[i].y) + poly[i].x))
				c = !c;
		}
	if (c == 0)
	{
		int near_id = findNearPoint(A, poly);
		double  d1, d2, d3;
		if (near_id == 0)
		{
			d1 = distance(A, poly[poly.size() - 1]);
			d2 = distance(A, poly[0]);
			d3 = distance(A, poly[1]);
		}
		else if (near_id == poly.size() - 1)
		{
			d1 = distance(A, poly[poly.size() - 2]);
			d2 = distance(A, poly[poly.size() - 1]);
			d3 = distance(A, poly[0]);
		}
		else
		{
			d1 = distance(A, poly[near_id - 1]);
			d2 = distance(A, poly[near_id]);
			d3 = distance(A, poly[near_id + 1]);
		}
		if (d2 < 0.5 * MeshPara::dx/*d1 / d2 > 2 || d2 / d1 > 2 || d2 / d3 > 2 || d3 / d2 > 2*/)
			c = 1;
	}
	return c;
}
bool judgeFieldInOut(double x, double y, vector <mesh>& poly)
//judge whether the point is inside the polygon
//html: https://blog.csdn.net/u011722133/article/details/52813374 
{
	float maxX, maxY, minX, minY;
	maxX = minX = poly[0].x;
	maxY = minY = poly[0].y;
	int i, j, c = 0;
	for (int i = 0; i < poly.size(); i++)
	{
		maxX = max(maxX, poly[i].x);
		maxY = max(maxY, poly[i].y);
		minX = min(minX, poly[i].x);
		minY = min(minY, poly[i].y);
	}
	if (x<minX || x>maxX || y<minY || y>maxY)
		c = 0;
	else
		for (i = 0, j = poly.size() - 1; i < poly.size(); j = i++)
		{
			//if (poly[j].y == poly[i].y && (poly[i].y > y) != (poly[j].y > y))
			//	return !c;
			/*else*/ if ((poly[i].y > y) != (poly[j].y > y) && (x < (poly[j].x - poly[i].x) * (y - poly[i].y) / (poly[j].y - poly[i].y) + poly[i].x))
				c = !c;
		}
	if (c == 0)
	{
		int near_id = findNearPoint(x, y, poly);
		double d1, d2, d3;
		if (near_id == 0)
		{
			d1 = distance(x, y, poly[poly.size() - 1].x, poly[poly.size() - 1].y);
			d2 = distance(x, y, poly[0].x, poly[0].y);
			d3 = distance(x, y, poly[1].x, poly[1].y);
		}
		else if (near_id == poly.size() - 1)
		{
			d1 = distance(x, y, poly[poly.size() - 2].x, poly[poly.size() - 2].y);
			d2 = distance(x, y, poly[poly.size() - 1].x, poly[poly.size() - 1].y);
			d3 = distance(x, y, poly[0].x, poly[0].y);
		}
		else
		{
			d1 = distance(x, y, poly[near_id - 1].x, poly[near_id - 1].y);
			d2 = distance(x, y, poly[near_id].x, poly[near_id].y);
			d3 = distance(x, y, poly[near_id + 1].x, poly[near_id + 1].y);
		}
		if (d2 < 0.5 * MeshPara::dx/*d1 / d2 > 2 || d2 / d1 > 2 || d2 / d3 > 2 || d3 / d2 > 2*/)
			c =!c;
	}
	return c;
}
bool judgeFieldInOut(double x, double y)
{
	int c = 0;
	double r1 = r + 0.5 * MeshPara::dx;
	if (x <= a)
	{
		if ((x - a) * (x - a) + (y - b) * (y - b) > r1* r1)
			return c;
		else
			return !c;
	}
	else
	{
		double k1 = tan(4.6 * ConstPara::pi / 180);
		double k2 = tan(-4.6 * ConstPara::pi / 180);
		double x1 = a, y1 = b + r1;
		double x2 = a, y2 = b - r1;
		double b1 = y1 - k1 * x1;
		double b2 = y2 - k2 * x2;
		if (k1 * x + b1 < y || k2 * x + b2 > y)
			return c;
		else
			return !c;
	}
}
void polygonPoint(vector <mesh>& poly)
{
	mesh c1, c2, c3, c4, c;
	VECTOR v1, v2, v3;
	double S = 400;
	if (poly.size() != 0)
		poly.clear();
	//Blunt body problem
	else
	{
		double alpha = 4.6 * pi / 180;
		c1.x = 1.0 + r;
		c1.y = 1.5 + r;
		c2.x = 3;
		c2.y = (3.0 - 1.0 - r) * tan(alpha) + 1.5 + r;
		c3.x = 3;
		c3.y = (3.0 - 1.0 - r) * tan(-alpha) + 1.5 - r;
		c4.x = 1.0 + r;
		c4.y = 1.5 - r;
		v1.x = c2.x - c1.x;
		v1.y = c2.y - c1.y;
		v2.x = c3.x - c2.x;
		v2.y = c3.y - c2.y;
		v3.x = c4.x - c3.x;
		v3.y = c4.y - c3.y;
		poly.push_back(c1);
		c = c1;
		while (c.x < c2.x)
		{
			if (c.x + v1.x / S > c2.x || c.y + v1.y / S > c2.y)
				c.x = c2.x, c.y = c2.y;
			else
				c.x = c.x + v1.x / S, c.y = c.y + v1.y / S;
			poly.push_back(c);
		}
		while (c.y > c3.y)
		{
			if (c.y + v2.y / S < c3.y)
				c.x = c3.x, c.y = c3.y;
			else
				c.x = c.x + v2.x / S, c.y = c.y + v2.y / S;
			poly.push_back(c);
		}
		while (c.x > c4.x)
		{
			if (c.x + v3.x / S < c4.x)
				c.x = c4.x, c.y = c4.y;
			else
				c.x = c.x + v3.x / S, c.y = c.y + v3.y / S;
			poly.push_back(c);
		}
		double beta = 3 * pi / 2;
		while (beta > pi / 2)
		{
			beta -= pi / S;
			c.x = r * cos(beta) + a;
			c.y = r * sin(beta) + b;
			poly.push_back(c);
		}
	}

	//cylinder problem
	//else
	//{
	//	float beta = 0;
	//	while (beta < 2 * pi)
	//	{
	//		beta += pi / 400;
	//		c.x = (r + 0.5 * MeshPara::dx) * cos(beta) + a;
	//		c.y = (r + 0.5 * MeshPara::dx) * sin(beta) + b;
	//		poly.push_back(c);
	//	}

	//}
}
mesh getCrossPoint(double theta, double a, double b, double r)
{
	mesh M;
	if (theta == pi / 2)
	{
		M.x = a;
		M.y = b + r;
	}
	else if (theta == 3 * pi / 2)
	{
		M.x = a;
		M.y = b - r;
	}
	else
	{
		M.x = r * cos(theta) + a;
		M.y = r * sin(theta) + b;
	}
	return M;
}
mesh getCrossPoint(Line L1, Line L2)
{
	mesh L;
	L.x = (L2.B * L1.C - L1.B * L2.C) / (L2.A * L1.B - L1.A * L2.B);
	L.y = (L2.A * L1.C - L1.A * L2.C) / (L1.A * L2.B - L2.A * L1.B);
	return L;
}
mesh getCrossPoint(mesh M, double a, double b, double r)//某点和圆心的连线与圆的交点
{
	mesh P;
	Line L = getLine(M.x, M.y, a, b);
	if (L.A == 0)
	{
		if (M.x < a)
			P.x = a - r, P.y = b;
		else
			P.x = a + r, P.y = b;
	}
	else if (L.B == 0)
	{
		if (M.y < b)
			P.x = a, P.y = b - r;
		else
			P.x = a, P.y = b + r;
	}
	else
	{
		double diff = 1e-15;
		double x0 = min(M.x, a);
		double x1 = max(M.x, a);
		double y0, y1;
		double x2, y2;
		double f0, f1, f2;
		while (abs(x0 - x1) > diff)
		{
			y0 = -(L.A * x0 + L.C) / L.B;
			y1 = -(L.A * x1 + L.C) / L.B;
			x2 = (x0 + x1) / 2;
			y2 = -(L.A * x2 + L.C) / L.B;
			f0 = (x0 - a) * (x0 - a) + (y0 - b) * (y0 - b) - r * r;
			f1 = (x1 - a) * (x1 - a) + (y1 - b) * (y1 - b) - r * r;
			f2 = (x2 - a) * (x2 - a) + (y2 - b) * (y2 - b) - r * r;
			if (f0 * f2 >= 0)
				x0 = x2;
			else
				x1 = x2;
		}
		P.x = (x0 + x1) / 2;
		P.y = (y0 + y1) / 2;
	}
	return P;
}
double get_beta(mesh A, mesh B)//求出两个网格点与x轴的夹角
{
	double dy = abs(A.y - B.y);
	double dx = abs(A.x - B.x);
	double beta = atan(dy / dx);
	return beta;
}
int findNearPoint(mesh A, vector<mesh>& poly)//find the closest point of given grid point
{
	int i, n = -1;
	double minD = distance(A.x, A.y, poly[0].x, poly[0].y);
	for (int i = 0; i < poly.size(); i++)
	{
		if (minD != min(minD, distance(A.x, A.y, poly[i].x, poly[i].y)))
		{
			minD = distance(A.x, A.y, poly[i].x, poly[i].y);
			n = i;
		}
	}
	return n;
}
int findNearPoint(double x, double y, vector<mesh>& poly)//find the closest point of given grid point
{
	int i, n = -1;
	double minD = distance(x, y, poly[0].x, poly[0].y);
	for ( i = 0; i < poly.size(); i++)
	{
		if (minD != min(minD, distance(x, y, poly[i].x, poly[i].y)))
		{
			minD = distance(x, y, poly[i].x, poly[i].y);
			n = i;
		}
	}
	return n;
}