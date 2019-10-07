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
double distance(double x1,double y1,double x2,double y2)
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
bool judgeFieldInOut(mesh& A, vector <mesh> &poly)
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
	return c;
}
bool judgeFieldInOut(double x, double y, vector <mesh> &poly)
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
	return c;
}

void polygonPoint(vector <mesh>& poly)
{
	mesh c;
	if (poly.size() != 0)
		poly.clear();
	//Blunt body problem
	//else
	//{
	//	float alpha = 4.6 * pi / 180;
	//	c.x = 1.0 + r;
	//	c.y = 1.5 + r;
	//	poly.push_back(c);
	//	c.x = 3;
	//	c.y = (3.0 - 1.0 - r) * tan(alpha) + 1.5 + r;
	//	poly.push_back(c);
	//	c.y = (3.0 - 1.0 - r) * tan(-alpha) + 1.5 - r;
	//	poly.push_back(c);
	//	c.x = 1.0 + r;
	//	c.y = 1.5 - r;
	//	poly.push_back(c);
	//	float beta = 3 * pi / 2;
	//	while (beta > pi / 2)
	//	{
	//		beta -= pi / 400;
	//		c.x = r * cos(beta) + a;
	//		c.y = r * sin(beta) + b;
	//		poly.push_back(c);
	//	}
	//}

	//cylinder problem
	else
	{
		float beta = 0;
		while (beta < 2 * pi)
		{
			beta += pi / 400;
			c.x = (r + 0.5 * MeshPara::dx) * cos(beta) + a;
			c.y = (r + 0.5 * MeshPara::dx) * sin(beta) + b;
			poly.push_back(c);
		}

	}
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
int findNearPoint(mesh A, vector<Coor> &poly)//find the closest point of given grid point
{
	int i, n;
	double maxD = -1;
	for (int i = 0; i < poly.size(); i++)
	{
		if (maxD != max(maxD, distance(A.x, A.y, poly[i].x, poly[i].y)))
		{
			maxD =  distance(A.x, A.y, poly[i].x, poly[i].y);
			n = i;
		}
	}
}