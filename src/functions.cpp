#include"include/const.h"
#include"include/functions.h"
#include<omp.h>
#include"include/shockwave.h"
#include<iostream>
#include"include/Prandtl-Meyer.h"
using std::vector;
using namespace ConstPara;
using namespace MeshPara;

double distance(mesh& a, mesh& b)
{
	double dx = a.x - b.x;
	double dy = a.y - b.y;
	return sqrt(dx * dx + dy * dy);
}
void get_dt()
{
	extern vector<vector <mesh*>> A;
	double maxxi = 0, maxeta = 0;
	extern double dt;
	double t;
	int i, j, k;
	double max1, max2;
	double Sxi, Seta, c, uxi, ueta;

	dt = t_end;
	max1 = max2 = 0;
	for (i = 0; i < A.size(); i++)
	{
		for (j = 0; j < A[i].size(); j++)
		{
			if (A[i][j]->type != "IN" || A[i][j]->section == 0)
				continue;
			maxxi = maxeta = 0;
			if (A[i][j]->neibor.size() > 3)
			{
				Sxi = sqrt(A[i][j]->xix[0] * A[i][j]->xix[0] + A[i][j]->xiy[0] * A[i][j]->xiy[0]);
				Seta = sqrt(A[i][j]->etax[0] * A[i][j]->etax[0] + A[i][j]->etay[0] * A[i][j]->etay[0]);
				c = sqrt(gama * A[i][j]->p / A[i][j]->rho);
				uxi = A[i][j]->u * A[i][j]->xix[0] + A[i][j]->v * A[i][j]->xiy[0];
				ueta = A[i][j]->u * A[i][j]->etax[0] + A[i][j]->v * A[i][j]->etay[0];
				//maxxi = max(Sxi, uxi);
				//maxeta = max(Seta, ueta);
				maxxi = abs(uxi) + c * Sxi;
				maxeta = abs(uxi) + c * Seta;
				max1 = max(max1, maxxi);
				max2 = max(max2, maxeta);

			}
			else
			{
				for (k = 0; k < 3; k++)
				{
					Sxi = sqrt(A[i][j]->xix[k] * A[i][j]->xix[k] + A[i][j]->xiy[k] * A[i][j]->xiy[k]);
					Seta = sqrt(A[i][j]->etax[k] * A[i][j]->etax[k] + A[i][j]->etay[k] * A[i][j]->etay[k]);
					c = sqrt(gama * A[i][j]->p / A[i][j]->rho);
					uxi = A[i][j]->u * A[i][j]->xix[k] + A[i][j]->v * A[i][j]->xiy[k];
					ueta = A[i][j]->u * A[i][j]->etax[k] + A[i][j]->v * A[i][j]->etay[k];
					//maxxi = max(Sxi, uxi);
					//maxeta = max(Seta, ueta);
					maxxi += abs(uxi) + c * Sxi;
					maxeta += abs(uxi) + c * Seta;
				}
				maxxi = maxxi / 3;
				maxeta = maxeta / 3;
				max1 = max(max1, maxxi);
				max2 = max(max2, maxeta);
			}
		}
	}
	t = CFL / (max1 + max2);
	//t = CFL / (maxxi + maxeta);
	dt = min(dt, t);
}
//void update_AfromU()
//{
//	extern vector<vector <mesh*>> A;
//	extern vector<vector<vector <double>>> U;
//	int i, j;
//	extern int step;
//
//	for (i = 0; i < A.size(); i++)
//	{
//#pragma omp parallel
//
//		for (j = 0; j < A[i].size(); j++)
//		{
//			A[i][j]->rho = U[i][j][0];
//			A[i][j]->u = U[i][j][1] / U[i][j][0];
//			A[i][j]->v = U[i][j][2] / U[i][j][0];
//			A[i][j]->p = (gama - 1) * (U[i][j][3] - 0.5 * A[i][j]->rho * (A[i][j]->u * A[i][j]->u + A[i][j]->v * A[i][j]->v));
//			A[i][j]->step = step;
//		}
//	}
//
//}
void update_bound_uniform()
{
	extern vector<vector <mesh*>> A;
	int i, j;
	using namespace Init;

#pragma omp parallel

	for (i = 0; i < A.size(); i++)
	{
		for (j = 0; j < A[i].size(); j++)
		{
			if (A[i][j]->type == "IN" || A[i][j]->type == "N")
				continue;
			A[i][j]->rho = rho0;
			A[i][j]->u = u0;
			A[i][j]->v = v0;
			A[i][j]->p = p0;
			//if (A[i].type == "L")
			//{
			//	A[i].rho = rho1;
			//	A[i].u = u1;
			//	A[i].v = v1;
			//	A[i].p = p1;

			//}
			if (A[i][j]->type != "L" || A[i][j]->type != "R")
				A[i][j]->v = 0;
		}
	}
}
void update_bound_shockwave()
{
	extern vector <mesh> A0;
	extern vector<vector <mesh*>> A;
	int i, j;
	if (FlowType == "normal")
	{
		using namespace Normal;
		for (i = 0; i < A.size(); i++)
		{
			for (j = 0; j < A[i].size(); j++)
			{
				if (A[i][j]->type == "IN")
					continue;
				else if (A[i][j]->type == "U")
				{
					A[i][j]->rho = A[i][j - Xnum]->rho;
					A[i][j]->u = A[i][j - Xnum]->u;
					A[i][j]->v = A[i][j - Xnum]->v;
					A[i][j]->p = A[i][j - Xnum]->p;
					if (A[i][j]->neibor.size() == 2)
					{
						A[i][j]->rho = (A[i][j]->neibor[0]->rho + A[i][j]->neibor[1]->rho) / 2;
						A[i][j]->u = (A[i][j]->neibor[0]->u + A[i][j]->neibor[1]->u) / 2;
						A[i][j]->v = (A[i][j]->neibor[0]->v + A[i][j]->neibor[1]->v) / 2;
						A[i][j]->p = (A[i][j]->neibor[0]->p + A[i][j]->neibor[1]->p) / 2;
					}
				}
				else if (A[i][j]->type == "D")
				{
					A[i][j]->rho = A[i][j + Xnum]->rho;
					A[i][j]->u = A[i][j + Xnum]->u;
					A[i][j]->v = A[i][j + Xnum]->v;
					A[i][j]->p = A[i][j + Xnum]->p;
					if (A[i][j]->neibor.size() == 2)
					{
						A[i][j]->rho = (A[i][j]->neibor[0]->rho + A[i][j]->neibor[1]->rho) / 2;
						A[i][j]->u = (A[i][j]->neibor[0]->u + A[i][j]->neibor[1]->u) / 2;
						A[i][j]->v = (A[i][j]->neibor[0]->v + A[i][j]->neibor[1]->v) / 2;
						A[i][j]->p = (A[i][j]->neibor[0]->p + A[i][j]->neibor[1]->p) / 2;
					}

				}
				else if (A[i][j]->type == "L")
				{
					A[i][j]->rho = rho1;
					A[i][j]->u = u1;
					A[i][j]->v = v1;
					A[i][j]->p = p1;
					if (A[i][j]->neibor.size() == 2)
					{
						A[i][j]->rho = (A[i][j]->neibor[0]->rho + A[i][j]->neibor[1]->rho) / 2;
						A[i][j]->u = (A[i][j]->neibor[0]->u + A[i][j]->neibor[1]->u) / 2;
						A[i][j]->v = (A[i][j]->neibor[0]->v + A[i][j]->neibor[1]->v) / 2;
						A[i][j]->p = (A[i][j]->neibor[0]->p + A[i][j]->neibor[1]->p) / 2;
					}

				}
				else if (A[i][j]->type == "R")
				{
					A[i][j]->rho = A[i][j - 1]->rho;
					A[i][j]->u = A[i][j - 1]->u;
					A[i][j]->v = A[i][j - 1]->v;
					A[i][j]->p = A[i][j - 1]->p;
					if (A[i][j]->neibor.size() == 2)
					{
						A[i][j]->rho = (A[i][j]->neibor[0]->rho + A[i][j]->neibor[1]->rho) / 2;
						A[i][j]->u = (A[i][j]->neibor[0]->u + A[i][j]->neibor[1]->u) / 2;
						A[i][j]->v = (A[i][j]->neibor[0]->v + A[i][j]->neibor[1]->v) / 2;
						A[i][j]->p = (A[i][j]->neibor[0]->p + A[i][j]->neibor[1]->p) / 2;
					}

				}

			}
		}
		for (i = 0; i < A.size(); i++)
		{
			for (j = 0; j < A[i].size(); j++)
			{
				if (A[i][j]->type == "U" || A[i][j]->type == "D")
				{
					if (A[i][j]->neibor.size() == 2)
					{
						if (abs(A[i][j]->x - A[i][j]->neibor[0]->x) < 1e-10 && abs(A[i][j]->y - A[i][j]->neibor[0]->y) < 1e-10)
						{
							A[i][j]->rho = A[i][j]->neibor[0]->rho;
							A[i][j]->u = A[i][j]->neibor[0]->u;
							A[i][j]->v = A[i][j]->neibor[0]->v;
							A[i][j]->p = A[i][j]->neibor[0]->p;
						}
						else
						{
							A[i][j]->rho = A[i][j]->neibor[1]->rho;
							A[i][j]->u = A[i][j]->neibor[1]->u;
							A[i][j]->v = A[i][j]->neibor[1]->v;
							A[i][j]->p = A[i][j]->neibor[1]->p;

						}
					}
				}

			}
		}
	}
	else if (FlowType == "intersection")
	{
		using namespace ShockwaveCross;
		double theta12 = -30 / 180.0 * ConstPara::pi;
		double theta13 = 30 / 180.0 * ConstPara::pi;
		mesh A12, A13, C;
		Line L12, L13;
		A12.x = -dx / 2, A12.y = A0[Pnum - 1].y - 2 * dy * 5 * Ynum / 45;
		A13.x = -dx / 2, A13.y = 2 * dy * 5 * Ynum / 45;
		L12 = getLine(theta12, A12);
		L13 = getLine(theta13, A13);
		C = getCrossPoint(L12, L13);
		L12 = getLine(beta2, C);
		L13 = getLine(beta3, C);
		mesh A1, A2, A3;
		A1.rho = rho1, A1.p = p1, A1.u = u1, A1.v = v1;
		get_down(A1, A2, beta2);
		get_down(A1, A3, beta3);

		rho2 = A2.rho, p2 = A2.p, u2 = A2.u, v2 = A2.v;
		rho3 = A3.rho, p3 = A3.p, u3 = A3.u, v3 = A3.v;
#pragma omp parallel
		for (i = 0; i < A.size(); i++)
		{


			for (j = 0; j < A[i].size(); j++)
			{
				if (A[i][j]->type == "IN")
					continue;
				else if (A[i][j]->type == "L")
				{
					if (A[i][j]->x * L12.A + A[i][j]->y * L12.B + L12.C < 0)
					{
						A[i][j]->rho = rho2;
						A[i][j]->u = u2;
						A[i][j]->v = v2;
						A[i][j]->p = p2;
					}
					else if (A[i][j]->x * L13.A + A[i][j]->y * L13.B + L13.C > 0)
					{
						A[i][j]->rho = rho3;
						A[i][j]->u = u3;
						A[i][j]->v = v3;
						A[i][j]->p = p3;
					}
					else
					{
						A[i][j]->rho = rho1;
						A[i][j]->u = u1;
						A[i][j]->v = v1;
						A[i][j]->p = p1;
					}
				}
				else if (A[i][j]->type == "R")
				{
					A[i][j]->rho = A[i][j - 1]->rho;
					A[i][j]->u = A[i][j - 1]->u;
					A[i][j]->v = A[i][j - 1]->v;
					A[i][j]->p = A[i][j - 1]->p;

				}
				else if (A[i][j]->type == "U")
				{
					A[i][j]->rho = A[i][j - Xnum]->rho;
					A[i][j]->u = A[i][j - Xnum]->u;
					A[i][j]->v = A[i][j - Xnum]->v;
					A[i][j]->p = A[i][j - Xnum]->p;

				}
				else if (A[i][j]->type == "D")
				{
					A[i][j]->rho = A[i][j + Xnum]->rho;
					A[i][j]->u = A[i][j + Xnum]->u;
					A[i][j]->v = A[i][j + Xnum]->v;
					A[i][j]->p = A[i][j + Xnum]->p;
				}
			}
		}
	}
	else if (FlowType == "oblique")
	{
		using namespace Oblique;
		mesh M1;
		Line L;
		vector<mesh> t;
		M1.x = 0.0014, M1.y = 0;

		L = getLine(beta, M1);

		for (i = 0; i < A.size(); i++)
		{
			for (j = 0; j < A[i].size(); j++)
			{
				if (A[i][j]->type == "IN")
					continue;
				else if (A[i][j]->type == "U")
				{

					A[i][j]->rho = A[i][j - Xnum]->rho;
					A[i][j]->u = A[i][j - Xnum]->u;
					A[i][j]->v = A[i][j - Xnum]->v;
					A[i][j]->p = A[i][j - Xnum]->p;

				}
				else if (A[i][j]->type == "D")
				{
					if (A[i][j]->x > 0.0075)
					{
						A[i][j]->rho = rho2;
						A[i][j]->u = u2;
						A[i][j]->v = v2;
						A[i][j]->p = p2;

					}
					else
					{
						A[i][j]->rho = A[i][j + Xnum]->rho;
						A[i][j]->u = A[i][j + Xnum]->u;
						A[i][j]->v = A[i][j + Xnum]->v;
						A[i][j]->p = A[i][j + Xnum]->p;

					}

				}
				else if (A[i][j]->type == "L")
				{
					A[i][j]->rho = rho1;
					A[i][j]->u = u1;
					A[i][j]->v = v1;
					A[i][j]->p = p1;
				}
				else if (A[i][j]->type == "R")
				{
					A[i][j]->rho = A[i][j - 1]->rho;
					A[i][j]->u = A[i][j - 1]->u;
					A[i][j]->v = A[i][j - 1]->v;
					A[i][j]->p = A[i][j - 1]->p;
				}

			}
		}
		for (i = 0; i < A.size(); i++)
		{
			for (j = 0; j < A[i].size(); j++)
			{
				if (A[i][j]->type == "U" || A[i][j]->type == "D")
				{
					if (A[i][j]->neibor.size() == 2)
					{
						if (abs(A[i][j]->x - A[i][j]->neibor[0]->x) < 1e-10 && abs(A[i][j]->y - A[i][j]->neibor[0]->y) < 1e-10)
						{
							A[i][j]->rho = A[i][j]->neibor[0]->rho;
							A[i][j]->u = A[i][j]->neibor[0]->u;
							A[i][j]->v = A[i][j]->neibor[0]->v;
							A[i][j]->p = A[i][j]->neibor[0]->p;
						}
						else
						{
							A[i][j]->rho = A[i][j]->neibor[1]->rho;
							A[i][j]->u = A[i][j]->neibor[1]->u;
							A[i][j]->v = A[i][j]->neibor[1]->v;
							A[i][j]->p = A[i][j]->neibor[1]->p;

						}
					}
				}

			}
		}

	}
	else if (FlowType == "couette")
	{
		using namespace Couette;
		extern double yU;
		for (i = 0; i < A.size(); i++)
		{
			for (j = 0; j < A[i].size(); j++)
			{
				if (A[i][j]->type == "IN")
					continue;
				else if (A[i][j]->type == "U")
				{
					A[i][j]->rho = rho1;
					A[i][j]->u = u1;
					A[i][j]->v = v1;
					A[i][j]->p = p1;
					A[i][j]->rho = A[i][j - Xnum]->rho;
					A[i][j]->u = A[i][j - Xnum]->u;
					A[i][j]->v = A[i][j - Xnum]->v;
					A[i][j]->p = A[i][j - Xnum]->p;

				}
				else if (A[i][j]->type == "D")
				{
					A[i][j]->rho = rho1;
					A[i][j]->u = 0;
					A[i][j]->v = v1;
					A[i][j]->p = p1;
				}
				else if (A[i][j]->type == "L")
				{
					A[i][j]->rho = rho1;
					A[i][j]->u = u1 * A[i][j]->y / yU;
					A[i][j]->v = v1;
					A[i][j]->p = p1;
				}
				A[i][j]->rho = rho1;
				A[i][j]->u = u1 * A[i][j]->y / yU;
				A[i][j]->v = v1;
				A[i][j]->p = p1;

				if (A[i][j]->type == "R")
				{
					A[i][j]->rho = A[i][j - 1]->rho;
					A[i][j]->u = A[i][j - 1]->u;
					A[i][j]->v = A[i][j - 1]->v;
					A[i][j]->p = A[i][j - 1]->p;
				}

			}
		}


	}
	else if (FlowType == "cylinder")
	{
		for (i = 0; i < A0.size(); i++)
		{
			if (A0[i].type == "L")
			{
				A0[i].rho = 20;
				A0[i].u = 10;
				A0[i].v = 0;
				A0[i].p = 1;
			}
			else if (A0[i].type == "R")
			{
				A0[i].rho = A0[i - 1].rho;
				A0[i].u = A0[i - 1].u;
				A0[i].v = A0[i - 1].v;
				A0[i].p = A0[i - 1].p;
			}
			else if (A0[i].type == "U")
			{
				A0[i].rho = A0[i - Xnum].rho;
				A0[i].u = A0[i - Xnum].u;
				A0[i].v = A0[i - Xnum].v;
				A0[i].p = A0[i - Xnum].p;
			}
			else if (A0[i].type == "D")
			{
				A0[i].rho = A0[i + Xnum].rho;
				A0[i].u = A0[i + Xnum].u;
				A0[i].v = A0[i + Xnum].v;
				A0[i].p = A0[i + Xnum].p;
			}
			else if (A0[i].type == "Cy")
			{

				double nx = A0[i].x - a;
				double ny = A0[i].y - b;
				double tx, ty;
				if (ny > 0)
				{
					tx = ny;
					ty = -nx;
				}
				else
				{
					tx = -ny;
					ty = nx;
				}
				double n1 = A0[i].neibor[0]->u;
				double n2 = A0[i].neibor[0]->v;
				if ((abs(n1) < 1e-10 && abs(n2) < 1e-10) || (abs(tx) < 1e-8 || abs(ty) < 1e-8))
				{
					A0[i].rho = A0[i].neibor[0]->rho;
					A0[i].u = A0[i].neibor[0]->u;
					A0[i].v = A0[i].neibor[0]->v;
					A0[i].p = A0[i].neibor[0]->p;
				}
				else
				{
					double costheta = (tx * n1 + ty * n2) / (sqrt(tx * tx + ty * ty) * sqrt(n1 * n1 + n2 * n2));
					double ex = tx / sqrt(tx * tx + ty * ty);
					double ey = ty / sqrt(tx * tx + ty * ty);
					double u = (n1 * n1 + n2 * n2) * costheta * ex;
					double v = (n1 * n1 + n2 * n2) * costheta * ey;
					//std::cout << costheta << "  " << ex << "  " << ey << std::endl;
					//std::cout << A0[i].neibor[0]->u << "  " << A0[i].neibor[0]->v << "  " << u << "  " << v << std::endl;
					A0[i].rho = A0[i].neibor[0]->rho;
					A0[i].u = u;
					A0[i].v = v;
					A0[i].p = A0[i].neibor[0]->p;
					//A0[i].rho = A0[i].neibor[0]->rho;
					//A0[i].u = A0[i].neibor[0]->u;
					//A0[i].v = A0[i].neibor[0]->v;
					//A0[i].p = A0[i].neibor[0]->p;

				}
				//A0[i].rho = 2;
				//A0[i].u = 0;
				//A0[i].v = 0;
				//A0[i].p = A0[i].p;

			}
		}
	}
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
		if ((x2 > x1 && y2 < y1) || (x2 < x1 && y2 > y1))
			theta = -theta;
	}
	return theta;
}
void update_Vm()
{
	extern double  yU, yD;
	extern vector<vector<mesh*>> A;
	extern vector<mesh>A0;
	int i, j, k, n;
	double d1;
	extern int step;
	int s = 1, m;
	double um, vm;
	if (FlowType == "oblique" || FlowType == "normal")
	{
		double x;
		if (FlowType == "oblique")
		{
			A[0][231]->um = A[0][231]->vm = 0;
			A[1][1419]->um = A[1][1419]->vm = 0;
		}

		//else if (FlowType == "normal")
		//{
		//	A[0][585].vm = 0;
		//	A[1][5915].vm = 0;
		//}
		for (i = 0; i < A.size(); i++)
		{
			for (j = 0; j < A[i].size(); j++)
			{
				if (A[i][j]->type == "SHOCK")
					continue;
				um = 0, vm = 0, x = 0;
				for (k = 0; k < A[i][j]->moveConnct.size(); k++)
				{
					m = A[i][j]->moveConnct[k];
					um += A[i][m]->um;
					vm += A[i][m]->vm;
					x += A[i][m]->x;
				}
				um = um / A[i][j]->moveConnct.size();
				vm = vm / A[i][j]->moveConnct.size();
				x = x / A[i][j]->moveConnct.size();
				if (i == 0)
				{
					A[i][j]->um = A[i][j]->x / x * um;
					A[i][j]->vm = A[i][j]->x / x * vm;
				}
				if (i == 1)
				{
					A[i][j]->um = (A0[Xnum - 1].x - A[i][j]->x) / (A0[Xnum - 1].x - x) * um;
					A[i][j]->vm = (A0[Xnum - 1].x - A[i][j]->x) / (A0[Xnum - 1].x - x) * vm;
				}
			}
		}
	}
	else if (FlowType == "intersection")
	{
		double x, y, x1, y1, x2, y2;
		double um1, vm1, um2, vm2;
		double d1, d2;
		int n1, n2;
		for (i = 0; i < A.size(); i++)
		{
			A[i][0]->um = A[i][0]->vm = 0;
			for (j = 0; j < A[i].size(); j++)
			{
				if (A[i][j]->type == "SHOCK")
				{
					continue;
				}
				else if (A[i][j]->type == "DISCON")
				{
					continue;
				}
				else if (A[i][j]->type == "CENTER")
				{
					A[i][j]->um = A[i][j]->vm = 0;
				}
				else if (A[i][j]->y == yU || A[i][j]->y == yD)
					A[i][j]->um = A[i][j]->vm = 0;
				else
				{

					if (i == 0 || i == 3 || i == 4)
					{
						A[i][j]->vm = 0;
						for (k = 0; k < A[i][j]->moveConnct.size(); k++)
						{
							m = A[i][j]->moveConnct[k];
							if (A[i][m]->y > A[i][j]->y)
								y1 = A[i][m]->y, vm1 = A[i][m]->vm;
							else
								y2 = A[i][m]->y, vm2 = A[i][m]->vm;
						}
					}

					else if (i == 1)
					{
						y1 = yU;
						y2 = A[i][A[i][j]->moveConnct[0]]->y;
						vm1 = 0;
						vm2 = A[i][A[i][j]->moveConnct[0]]->vm;
					}
					else
					{
						y1 = A[i][A[i][j]->moveConnct[0]]->y;
						y2 = yD;
						vm1 = A[i][A[i][j]->moveConnct[0]]->vm;
						vm2 = 0;
					}
					y = A[i][j]->y;
					double con = (y1 - y) / (y - y2);
					A[i][j]->vm = (y1 + vm1 + con * (y2 + vm2) - (y + con * y)) / (con + 1);
					A[i][j]->um = 0;
				}



			}

		}


	}

}
void clear_Vm()
{
	extern vector<vector<mesh*>> A;
	int i, j;
#pragma omp parallel
	for (i = 0; i < A.size(); i++)
	{
		for (j = 0; j < A[i].size(); j++)
		{
			A[i][j]->um = 0;
			A[i][j]->vm = 0;
		}
	}
}
void update_bound()
{
	extern vector<vector<mesh*>> A;

	int i, j, k, m, n;
	if (methodType == "C")
	{
		if (FlowType == "normal")
		{
			using namespace Normal;
#pragma omp parallel

			for (i = 0; i < A.size(); i++)
			{

				for (j = 0; j < A[i].size(); j++)
				{
					if (A[i][j]->type == "IN" || A[i][j]->type == "SHOCK")
						continue;
					if (i == 0)
					{
						A[i][j]->rho = rho1;
						A[i][j]->u = u1;
						A[i][j]->v = v1;
						A[i][j]->p = p1;
					}
					if (i == 1)
					{
						A[i][j]->rho = rho2;
						A[i][j]->u = u2;
						A[i][j]->v = v2;
						A[i][j]->p = p2;
					}

				}
			}
		}
		if (FlowType == "oblique")
		{
			using namespace Oblique;
#pragma omp parallel
			for (i = 0; i < A.size(); i++)
			{


				for (j = 0; j < A[i].size(); j++)
				{
					if (A[i][j]->type == "IN" || A[i][j]->type == "SHOCK")
						continue;
					if (i == 0)
					{
						A[i][j]->rho = rho1;
						A[i][j]->u = u1;
						A[i][j]->v = v1;
						A[i][j]->p = p1;
					}
					if (i == 1)
					{
						A[i][j]->rho = rho2;
						A[i][j]->u = u2;
						A[i][j]->v = v2;
						A[i][j]->p = p2;
					}
					if (A[i][j]->type == "R")
					{
						A[i][j]->rho = A[i][j - 1]->rho;
						A[i][j]->u = A[i][j - 1]->u;
						A[i][j]->v = A[i][j - 1]->v;
						A[i][j]->p = A[i][j - 1]->p;
					}

				}
			}
		}
	}
	else if (methodType == "F")
	{
		if (FlowType == "normal")
		{
			using namespace Normal;
			for (i = 0; i < A.size(); i++)
			{

				for (j = 0; j < A[i].size(); j++)
				{
					if (A[i][j]->type == "IN" || A[i][j]->type == "SHOCK" || A[i][j]->type == "DISCON")
						continue;
					else if (A[i][j]->type == "U")
					{
						for (k = 0; k < A[i][j]->neibor.size(); k++)
						{
							m = A[i][j]->neibor[k]->id;
							if (A[i][j]->neibor[k]->type == "SHOCK")
								continue;
							else if (abs(A[i][j]->y - A[i][j]->neibor[k]->y) > 1e-10)
							{
								n = m;
								break;
							}
							else
								n = m;
						}
						A[i][j]->rho = A[i][n]->rho;
						A[i][j]->u = A[i][n]->u;
						A[i][j]->v = A[i][n]->v;
						A[i][j]->p = A[i][n]->p;
					}
					else if (A[i][j]->type == "D")
					{
						if (i == 0)
						{
							A[i][j]->rho = rho1;
							A[i][j]->u = u1;
							A[i][j]->v = v1;
							A[i][j]->p = p1;
						}
						else
						{
							A[i][j]->rho = rho2;
							A[i][j]->u = u2;
							A[i][j]->v = v2;
							A[i][j]->p = p2;

						}

					}
					else if (A[i][j]->type == "L")
					{
						A[i][j]->rho = rho1;
						A[i][j]->u = u1;
						A[i][j]->v = v1;
						A[i][j]->p = p1;
					}
					else if (A[i][j]->type == "R")
					{
						for (k = 0; k < A[i][j]->neibor.size(); k++)
						{
							if (A[i][j]->neibor[k]->type == "SHOCK")
								continue;
							else if (abs(A[i][j]->x - A[i][m]->x) > 1e-10)
							{
								n = m;
								break;
							}
							else
								n = m;
						}
						A[i][j]->rho = A[i][n]->rho;
						A[i][j]->u = A[i][n]->u;
						A[i][j]->v = A[i][n]->v;
						A[i][j]->p = A[i][n]->p;
					}
				}
			}

		}
		else if (FlowType == "oblique")
		{
			using namespace Oblique;

			for (i = 0; i < A.size(); i++)
			{

				for (j = 0; j < A[i].size(); j++)
				{
					if (A[i][j]->type == "IN" || A[i][j]->type == "SHOCK" || A[i][j]->type == "DISCON")
						continue;
					else if (A[i][j]->type == "U")
					{
						for (k = 0; k < A[i][j]->neibor.size(); k++)
						{
							m = A[i][j]->neibor[k]->id;
							if (A[i][j]->neibor[k]->type == "SHOCK")
								continue;
							else if (abs(A[i][j]->y - A[i][m]->y) > 1e-10)
							{
								n = m;
								break;
							}
							else
								n = m;
						}
						A[i][j]->rho = A[i][n]->rho;
						A[i][j]->u = A[i][n]->u;
						A[i][j]->v = A[i][n]->v;
						A[i][j]->p = A[i][n]->p;
					}
					else if (A[i][j]->type == "D")
					{
						if (i == 0)
						{
							A[i][j]->rho = rho1;
							A[i][j]->u = u1;
							A[i][j]->v = v1;
							A[i][j]->p = p1;
						}
						else
						{
							A[i][j]->rho = rho2;
							A[i][j]->u = u2;
							A[i][j]->v = v2;
							A[i][j]->p = p2;

						}

					}
					else if (A[i][j]->type == "L")
					{
						A[i][j]->rho = rho1;
						A[i][j]->u = u1;
						A[i][j]->v = v1;
						A[i][j]->p = p1;
					}
					else if (A[i][j]->type == "R")
					{
						for (k = 0; k < A[i][j]->neibor.size(); k++)
						{
							if (A[i][j]->neibor[k]->type == "SHOCK")
								continue;
							else if (abs(A[i][j]->x - A[i][m]->x) > 1e-10)
							{
								n = m;
								break;
							}
							else
								n = m;
						}
						A[i][j]->rho = A[i][n]->rho;
						A[i][j]->u = A[i][n]->u;
						A[i][j]->v = A[i][n]->v;
						A[i][j]->p = A[i][n]->p;
					}
				}
			}
		}
		else if (FlowType == "intersection")
		{
			using namespace ShockwaveCross;
			for (i = 0; i < A.size(); i++)
			{

				for (j = 0; j < A[i].size(); j++)
				{
					if (A[i][j]->type == "IN")
						continue;
					else if (A[i][j]->type == "U" || A[i][j]->type == "D")
					{
						for (k = 0; k < A[i][j]->neibor.size(); k++)
						{
							m = A[i][j]->neibor[k]->id;
							if (A[i][j]->neibor[k]->type == "SHOCK")
								continue;
							else if (abs(A[i][j]->y - A[i][m]->y) > 1e-10)
							{
								n = m;
								break;
							}
							else
								n = m;
						}
						A[i][j]->rho = A[i][n]->rho;
						A[i][j]->u = A[i][n]->u;
						A[i][j]->v = A[i][n]->v;
						A[i][j]->p = A[i][n]->p;
					}
					else if (A[i][j]->type == "L")
					{
						if (i == 0)
						{
							A[i][j]->rho = rho1;
							A[i][j]->u = u1;
							A[i][j]->v = v1;
							A[i][j]->p = p1;
						}
						else if (i == 1)
						{
							A[i][j]->rho = rho2;
							A[i][j]->u = u2;
							A[i][j]->v = v2;
							A[i][j]->p = p2;
						}
						else if (i == 2)
						{
							A[i][j]->rho = rho3;
							A[i][j]->u = u3;
							A[i][j]->v = v3;
							A[i][j]->p = p3;

						}
					}
					else if (A[i][j]->type == "R")
					{
						for (k = 0; k < A[i][j]->neibor.size(); k++)
						{
							if (A[i][j]->neibor[k]->type == "SHOCK")
								continue;
							else if (abs(A[i][j]->x - A[i][m]->x) > 1e-10)
							{
								n = m;
								break;
							}
							else
								n = m;
						}
						A[i][j]->rho = A[i][n]->rho;
						A[i][j]->u = A[i][n]->u;
						A[i][j]->v = A[i][n]->v;
						A[i][j]->p = A[i][n]->p;
					}
				}
			}

		}
	}
}
void update_bound_shockwaveCross()
{
	using namespace ShockwaveCross;
	extern vector<vector <mesh*>> A;
	int i, j;
#pragma omp parallel

	for (i = 0; i < A.size(); i++)
	{
		for (j = 0; j < A[i].size(); j++)
		{
			if (A[i][j]->type == "IN")
				continue;
			else if (A[i][j]->type == "U")
			{
				A[i][j]->rho = rho2;
				A[i][j]->u = u2;
				A[i][j]->v = v2;
				A[i][j]->p = p2;
			}
			else if (A[i][j]->type == "L")
			{
				A[i][j]->rho = rho1;
				A[i][j]->u = u1;
				A[i][j]->v = v1;
				A[i][j]->p = p1;

			}
			else if (A[i][j]->type == "D")
			{
				A[i][j]->rho = rho3;
				A[i][j]->u = u3;
				A[i][j]->v = v3;
				A[i][j]->p = p3;
			}
			else if (A[i][j]->type == "R")
			{
				A[i][j]->rho = A[i][j - 1]->rho;
				A[i][j]->u = A[i][j - 1]->u;
				A[i][j]->v = A[i][j - 1]->v;
				A[i][j]->p = A[i][j - 1]->p;
			}

		}
	}
}
void update_bound_Prandtl_Meyer()
{
	using namespace Prandtl_Meyer;
	extern vector<vector <mesh*>> A;
	int i, j;
#pragma omp parallel

	for (i = 0; i < A.size(); i++)
	{
		for (j = 0; j < A[i].size(); j++)
		{
			if (A[i][j]->type == "IN")
				continue;
			else if (A[i][j]->type == "L")
			{
				A[i][j]->rho = rho1;
				A[i][j]->u = u1;
				A[i][j]->v = v1;
				A[i][j]->p = p1;
			}
			else if (A[i][j]->type == "DL")
			{
				A[i][j]->rho = rho1;
				A[i][j]->u = u1;
				A[i][j]->v = v1;
				A[i][j]->p = p1;
			}
			else if (A[i][j]->type == "DR")
			{
				A[i][j]->rho = rho2;
				A[i][j]->u = u2;
				A[i][j]->v = v2;
				A[i][j]->p = p2;
			}
			else if (A[i][j]->type == "R")
			{
				A[i][j]->rho = A[i][j - 1]->rho;
				A[i][j]->u = A[i][j - 1]->u;
				A[i][j]->v = A[i][j - 1]->v;
				A[i][j]->p = A[i][j - 1]->p;
			}
		}
	}
}



//void get_F()
//{
//	extern vector <mesh> A;
//	extern double F[Pnum][4];
//	for (int i = 0; i < Pnum; i++)
//	{
//		F[i][0] = A[i].rho*A[i].u;
//		F[i][1] = A[i].rho*A[i].u*A[i].u + A[i].p;
//		F[i][2] = A[i].rho*A[i].u*A[i].v;
//		F[i][3] = A[i].u*(A[i].p / (gama - 1) + 0.5*A[i].rho*(A[i].u*A[i].u + A[i].v*A[i].v) + A[i].p);
//	}
//}
//void get_G()
//{
//	extern vector <mesh> A;
//	extern double G[Pnum][4];
//	for (int i = 0; i < Pnum; i++)
//	{
//		G[i][0] = A[i].rho*A[i].v; extern vector<vector<vector <double>>> Utr;
//		G[i][1] = A[i].rho*A[i].u *A[i].v;
//		G[i][2] = A[i].rho*A[i].u*A[i].v + A[i].p;
//		G[i][3] = A[i].v*(A[i].p / (gama - 1) + 0.5*A[i].rho*(A[i].u*A[i].u + A[i].v*A[i].v) + A[i].p);
//	}
//}

//void update_INtrans()
//{
//	extern vector <mesh> A0;
//	extern vector<vector<vector <double>>> U;
//	extern vector<vector<vector <double>>> Utr;
//	extern int i;
//
//#pragma omp parallel 
//
//	for (int i = 0; i < A0.size(); i++)
//	{
//		if (A0[i].type != "IN")
//			continue;
//		if (A0[i].neibor.size() != 3)
//		{
//			Utr[i][0][0] = U[0][i][0] * A0[i].J[0];
//			Utr[i][0][1] = U[0][i][1] * A0[i].J[0];
//			Utr[i][0][2] = U[0][i][2] * A0[i].J[0];
//			Utr[i][0][3] = U[0][i][3] * A0[i].J[0];
//		}
//		else
//		{
//			Utr[i][0][0] = U[0][i][0] * A0[i].J[0];
//			Utr[i][0][1] = U[0][i][1] * A0[i].J[0];
//			Utr[i][0][2] = U[0][i][2] * A0[i].J[0];
//			Utr[i][0][3] = U[0][i][3] * A0[i].J[0];
//
//			Utr[i][1][0] = U[0][i][0] * A0[i].J[0];
//			Utr[i][1][1] = U[0][i][1] * A0[i].J[0];
//			Utr[i][1][2] = U[0][i][2] * A0[i].J[0];
//			Utr[i][1][3] = U[0][i][3] * A0[i].J[0];
//
//			Utr[i][2][0] = U[0][i][0] * A0[i].J[0];
//			Utr[i][2][1] = U[0][i][1] * A0[i].J[0];
//			Utr[i][2][2] = U[0][i][2] * A0[i].J[0];
//			Utr[i][2][3] = U[0][i][3] * A0[i].J[0];
//		}
//	}
//}
double compute_res()//计算残差
{
	extern vector<vector <mesh*>> A;
	extern   vector <mesh> Ar;
	extern double dt;
	int i, j;
	double res = 0;
	int n = 0;
	for (i = 0; i < A.size(); i++)
	{
		for (j = 0; j < A[i].size(); j++)
		{
			res = max(res, abs(A[i][j]->rho - Ar[A[i][j]->id].rho) / Ar[A[i][j]->id].rho);
			n++;
		}
	}
	return res / n;
}
void record()
{
	extern vector <mesh> A0;
	extern vector <mesh> Ar;
	extern double t_sim;
	int i, j;
	vector<mesh>t;
	for (i = 0; i < A0.size(); i++)
	{
		if (t_sim == 0)
			Ar.push_back(A0[i]);
		else
			Ar[i] = A0[i];
	}

}

void update_IN()
{
	extern vector<vector <mesh*>> A;
	extern vector <mesh> Ar;
	double U[4], U1[4], U2[4];
	double temp[4] = { 0 };
	extern double dt;
	int i, j, k;
	Flux F1, F2, G1, G2;
	Flux F1_, F2_, G1_, G2_;
	int n1, n2, n3, n4;

	double d1, d2;
	Flux Rf_xieta, Rg_xieta;
	Flux Rf_phiψ, Rg_phiψ;
	Flux FL, FC, FR;
	mesh CL, CR, C;
	Flux GU, GC, GD;
	mesh CU, CD;

	Flux Fll, Flr, Fcl, Fcr, Frl, Frr;
	Flux Gdd, Gdu, Gcd, Gcu, Gud, Guu;

	Flux fll, flr, fcl, fcr, frl, frr;
	Flux gdd, gdu, gcd, gcu, gud, guu;
	double c1, c2;
	int id;
	//#pragma omp parallel for private(F1,F2,G1,G2,temp)
	extern double t_sim;
	for (i = 0; i < A.size(); i++)
	{
		for (j = 0; j < A[i].size(); j++)
		{
			if (A[i][j]->type == "L" || A[i][j]->type == "R" || A[i][j]->type == "U" || A[i][j]->type == "D" || A[i][j]->type == "Cy" || A[i][j]->section == 0 /*|| A[i][j]->change == "N"*/)
				continue;
			id = A[i][j]->id;
			U[0] = Ar[id].rho;
			U[1] = Ar[id].rho * A[i][j]->u;
			U[2] = Ar[id].rho * A[i][j]->v;
			U[3] = 0.5 * Ar[id].rho * (Ar[id].u * Ar[id].u + Ar[id].v * Ar[id].v) + Ar[id].p / (gama - 1);

			if (A[i][j]->neibor.size() == 3)
			{
				F1 = F2 = G1 = G2 = { 0 };

				for (k = 0; k < 12; k++)
				{
					n1 = method[k][0];
					n2 = method[k][1];
					n3 = method[k][2];
					n4 = method[k][3];
					n1 = A[i][j]->neibor[n1]->id;
					n2 = A[i][j]->neibor[n2]->id;
					n3 = A[i][j]->neibor[n3]->id;
					n4 = A[i][j]->neibor[n4]->id;
					id = A[i][j]->id;
					F1_ = HLLC_Χ(Ar[n3], Ar[id], Ar[id], k);
					F2_ = HLLC_Χ(Ar[id], Ar[n1], Ar[id], k);
					G1_ = HLLC_Υ(Ar[n4], Ar[id], Ar[id], k);
					G2_ = HLLC_Υ(Ar[id], Ar[n2], Ar[id], k);
					F1.f1 += F1_.f1, F1.f2 += F1_.f2, F1.f3 += F1_.f3, F1.f4 += F1_.f4;
					F2.f1 += F2_.f1, F2.f2 += F2_.f2, F2.f3 += F2_.f3, F2.f4 += F2_.f4;
					G1.f1 += G1_.f1, G1.f2 += G1_.f2, G1.f3 += G1_.f3, G1.f4 += G1_.f4;
					G2.f1 += G2_.f1, G2.f2 += G2_.f2, G2.f3 += G2_.f3, G2.f4 += G2_.f4;
				}
				F1.f1 /= k, F1.f2 /= k, F1.f3 /= k, F1.f4 /= k;
				F2.f1 /= k, F2.f2 /= k, F2.f3 /= k, F2.f4 /= k;
				G1.f1 /= k, G1.f2 /= k, G1.f3 /= k, G1.f4 /= k;
				G2.f1 /= k, G2.f2 /= k, G2.f3 /= k, G2.f4 /= k;

				U[0] = U[0] - dt * (F2.f1 - F1.f1) - dt * (G2.f1 - G1.f1);
				U[1] = U[1] - dt * (F2.f2 - F1.f2) - dt * (G2.f2 - G1.f2);
				U[2] = U[2] - dt * (F2.f3 - F1.f3) - dt * (G2.f3 - G1.f3);
				U[3] = U[3] - dt * (F2.f4 - F1.f4) - dt * (G2.f4 - G1.f4);
			}

			else if (A[i][j]->neibor.size() >= 4)
			{
				F1 = F2 = G1 = G2 = { 0 };
				n1 = A[i][j]->neibor[0]->id;
				n2 = A[i][j]->neibor[1]->id;
				n3 = A[i][j]->neibor[2]->id;
				n4 = A[i][j]->neibor[3]->id;
				//VanLeer
				Fll = VanLeerB(Ar[n3], Ar[id].xix[0], Ar[id].xiy[0], Ar[id].xit[0], Ar[id].J[0]);
				Flr = VanLeerA(Ar[n3], Ar[id].xix[0], Ar[id].xiy[0], Ar[id].xit[0], Ar[id].J[0]);
				Fcl = VanLeerB(Ar[id], Ar[id].xix[0], Ar[id].xiy[0], Ar[id].xit[0], Ar[id].J[0]);
				Fcr = VanLeerA(Ar[id], Ar[id].xix[0], Ar[id].xiy[0], Ar[id].xit[0], Ar[id].J[0]);
				Frl = VanLeerB(Ar[n1], Ar[id].xix[0], Ar[id].xiy[0], Ar[id].xit[0], Ar[id].J[0]);
				Frr = VanLeerA(Ar[n1], Ar[id].xix[0], Ar[id].xiy[0], Ar[id].xit[0], Ar[id].J[0]);

				Gdd = VanLeerB(Ar[n4], Ar[id].etax[0], Ar[id].etay[0], Ar[id].etat[0], Ar[id].J[0]);
				Gdu = VanLeerA(Ar[n4], Ar[id].etax[0], Ar[id].etay[0], Ar[id].etat[0], Ar[id].J[0]);
				Gcd = VanLeerB(Ar[id], Ar[id].etax[0], Ar[id].etay[0], Ar[id].etat[0], Ar[id].J[0]);
				Gcu = VanLeerA(Ar[id], Ar[id].etax[0], Ar[id].etay[0], Ar[id].etat[0], Ar[id].J[0]);
				Gud = VanLeerB(Ar[n2], Ar[id].etax[0], Ar[id].etay[0], Ar[id].etat[0], Ar[id].J[0]);
				Guu = VanLeerA(Ar[n2], Ar[id].etax[0], Ar[id].etay[0], Ar[id].etat[0], Ar[id].J[0]);
				U[0] = U[0] - dt * A[i][j]->J[0] * (Fcr.f1 - Flr.f1 + Frl.f1 - Fcl.f1 + Gcu.f1 - Gdu.f1 + Gud.f1 - Gcd.f1);
				U[1] = U[1] - dt * A[i][j]->J[0] * (Fcr.f2 - Flr.f2 + Frl.f2 - Fcl.f2 + Gcu.f2 - Gdu.f2 + Gud.f2 - Gcd.f2);
				U[2] = U[2] - dt * A[i][j]->J[0] * (Fcr.f3 - Flr.f3 + Frl.f3 - Fcl.f3 + Gcu.f3 - Gdu.f3 + Gud.f3 - Gcd.f3);
				U[3] = U[3] - dt * A[i][j]->J[0] * (Fcr.f4 - Flr.f4 + Frl.f4 - Fcl.f4 + Gcu.f4 - Gdu.f4 + Gud.f4 - Gcd.f4);
			}
			else if (A[i][j]->neibor.size() == 2)
			{
				if (A[i][j]->type == "SHOCK")
				{
					double m, n;

					A[i][j]->rho = (A[i][j]->neibor[0]->rho + A[i][j]->neibor[1]->rho) / 2;
					A[i][j]->u = (A[i][j]->neibor[0]->u + A[i][j]->neibor[1]->u) / 2;
					A[i][j]->v = (A[i][j]->neibor[0]->v + A[i][j]->neibor[1]->v) / 2;
					A[i][j]->p = (A[i][j]->neibor[0]->p + A[i][j]->neibor[1]->p) / 2;
				}
				else
					continue;
			}
			else
				std::cout << i << "   " << j << "   " << "somthing wrong in updat_u" << std::endl;
			A[i][j]->rho = U[0];
			A[i][j]->u = U[1] / U[0];
			A[i][j]->v = U[2] / U[0];
			A[i][j]->p = (gama - 1) * (U[3] - 0.5 * A[i][j]->rho * (A[i][j]->u * A[i][j]->u + A[i][j]->v * A[i][j]->v));
			if (i == 0 & FlowType == "cylinder")
			{
				A[i][j]->rho = A[i][j]->rho * A[i][j]->sec_num;
				A[i][j]->u = A[i][j]->u * A[i][j]->sec_num;
				A[i][j]->v = A[i][j]->v * A[i][j]->sec_num;
				A[i][j]->p = A[i][j]->p * A[i][j]->sec_num;
			}
		}
	}

}
Flux get_F(mesh N, mesh C, int method)//得到当地坐标系下的通量
{
	double xix = C.xix[method];
	double xiy = C.xiy[method];
	double xit = C.xit[method];
	double J = C.J[method];
	double Dxi = sqrt(xix * xix + xiy * xiy);
	double xi1 = xix / Dxi;
	double xi2 = xiy / Dxi;
	double xi3 = xit / Dxi;
	double ub = N.u * xi1 + N.v * xi2 + xi3;
	Flux F;
	F.f1 = C.rho;
	F.f2 = C.rho * C.u;
	F.f3 = C.rho * C.v;
	F.f4 = C.p * gama / (gama - 1) + 0.5 * C.rho * C.u * C.u + 0.5 * C.rho * C.v * C.v;
	F.f1 = F.f1 * ub;
	F.f2 = F.f2 * ub + xi1 * C.p;
	F.f3 = F.f3 * ub + xi2 * C.p;
	F.f4 = F.f4 * ub - xi3 * C.p;
	F.f1 = F.f1 * Dxi;
	F.f2 = F.f2 * Dxi;
	F.f3 = F.f3 * Dxi;
	F.f4 = F.f4 * Dxi;
	return F;

}
Flux get_G(mesh N, mesh C, int method)//得到当地坐标系下的通量
{
	double etax = C.etax[method];
	double etay = C.etay[method];
	double etat = C.etat[method];
	double J = C.J[method];
	double Deta = sqrt(etax * etax + etay * etay);
	double eta1 = etax / Deta;
	double eta2 = etay / Deta;
	double eta3 = etat / Deta;
	double ub = N.u * eta1 + N.v * eta2 + eta3;

	Flux F;
	F.f1 = C.rho;
	F.f2 = C.rho * C.u;
	F.f3 = C.rho * C.v;
	F.f4 = C.p * gama / (gama - 1) + 0.5 * C.rho * C.u * C.u + 0.5 * C.rho * C.v * C.v;
	F.f1 = F.f1 * ub;
	F.f2 = F.f2 * ub + eta1 * C.p;
	F.f3 = F.f3 * ub + eta2 * C.p;
	F.f4 = F.f4 * ub - eta3 * C.p;
	F.f1 = F.f1 * Deta;
	F.f2 = F.f2 * Deta;
	F.f3 = F.f3 * Deta;
	F.f4 = F.f4 * Deta;
	return F;

}
//void choose_U(int i)
//{
//	extern vector <mesh> A;
//	extern vector <mesh> A;
//	extern vector<vector <double>> U;
//	extern double U0[Pnum][4];
//	extern double U1[Pnum][4];
//	extern double U2[Pnum][4];
//	double rho00, u00, v00, p00;
//	double rho, u, v, p;
//	double rho0, u0, v0, p0;
//	double rho1, u1, v1, p1;
//	double rho2, u2, v2, p2;
//	rho00 = A[i].rho;
//	u00 = A[i].u;
//	v00 = A[i].v;
//	p00 = A[i].p;
//
//	rho0 = U0[i][0] / A[i].J[0];
//	u0 = U0[i][1] / U0[i][0];
//	v0 = U0[i][2] / U0[i][0];
//	p0 = (gama - 1)*(U0[i][3] / A[i].J[0] - 0.5*A[i].rho*(A[i].u*A[i].u + A[i].v*A[i].v));
//
//	rho1 = U1[i][0] / A[i].J[1];
//	u1 = U1[i][1] / U1[i][0];
//	v1 = U1[i][2] / U1[i][0];
//	p1 = (gama - 1)*(U1[i][3] / A[i].J[1] - 0.5*A[i].rho*(A[i].u*A[i].u + A[i].v*A[i].v));
//	//rho0 = rho1;
//	//u0 = u1;
//	//v0 = v1;
//	//p0 = p1;
//
//	rho2 = U2[i][0] / A[i].J[2];
//	u2 = U2[i][1] / U2[i][0];
//	v2 = U2[i][2] / U2[i][0];
//	p2 = (gama - 1)*(U2[i][3] / A[i].J[1] - 0.5*A[i].rho*(A[i].u*A[i].u + A[i].v*A[i].v));
//	rho = absmax(absmax(rho0, rho1), rho2);
//	u = absmin(absmin(u0, u1), u2);
//	v = absmax(absmax(v0, v1), v2);
//	p = absmax(absmax(p0, p1), p2);
//	U[i][0] = rho;
//	U[i][1] = rho * u;
//	U[i][2] = rho * v;
//	U[i][3] = 0.5*rho*(u*u + v * v) + p / (gama - 1);
//}
//void choose_U(int i)
//{
//	extern vector <mesh> A;
//
//	extern double U[Pnum][4];
//	extern double U0[Pnum][4];
//	extern double U1[Pnum][4];
//	extern double U2[Pnum][4];
//
//	double u0 = abs(U0[i][1] / U0[i][0]);
//	double u1 = abs(U1[i][1] / U1[i][0]);
//	double u2 = abs(U2[i][1] / U2[i][0]);
//	if (u0 < u1&&u0 < u2)
//	{
//		U[i][0] = U0[i][0] / A[i].J[0];
//		U[i][1] = U0[i][1] / A[i].J[0];
//		U[i][2] = U0[i][2] / A[i].J[0];
//		U[i][3] = U0[i][3] / A[i].J[0];
//	}
//	else if (u1 < u0&&u1 < u2)
//	{
//		U[i][0] = U1[i][0] / A[i].J[1];
//		U[i][1] = U1[i][1] / A[i].J[1];
//		U[i][2] = U1[i][2] / A[i].J[1];
//		U[i][3] = U1[i][3] / A[i].J[1];
//	}
//	else
//	{
//		U[i][0] = U2[i][0] / A[i].J[2];
//		U[i][1] = U2[i][1] / A[i].J[2];
//		U[i][2] = U2[i][2] / A[i].J[2];
//		U[i][3] = U2[i][3] / A[i].J[2];
//	}
//
//}
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
double get_beta(mesh A, mesh B)//求出两个网格点与x轴的夹角
{
	double dy = abs(A.y - B.y);
	double dx = abs(A.x - B.x);
	double beta = atan(dy / dx);
	return beta;
}
void reorder_neighbor()
//make the neighbor point store in this order
//R,U,L,D
{
	extern vector<vector <mesh*>> A;
	double maxy, miny, maxx, minx;
	mesh* t;
	for (int i = 0; i < A.size(); i++)
	{
		for (int j = 0; j < A[i].size(); j++)
		{
			if (A[i][j]->type != "IN" || A[i][j]->neibor.size() < 3)
				continue;
			if (A[i][j]->neibor.size() == 4)
			{
				mesh* n1 = A[i][j]->neibor[0];
				mesh* n2 = A[i][j]->neibor[1];
				mesh* n3 = A[i][j]->neibor[2];
				mesh* n4 = A[i][j]->neibor[3];
				double maxy = max(max(max(n1->y, n2->y), n3->y), n4->y);
				if (n1->y == maxy)
					t = n1, n1 = n2, n2 = t;
				else if (n3->y == maxy)
					t = n3, n3 = n2, n2 = t;
				else if (n4->y == maxy)
					t = n4, n4 = n2, n2 = t;
				double miny = min(min(n1->y, n3->y), n4->y);
				if (n1->y == miny)
					t = n1, n1 = n4, n4 = t;
				else if (n3->y == miny)
					t = n3, n3 = n4, n4 = t;
				maxx = max(n1->x, n3->x);
				if (n3->x == maxx)
					t = n1, n1 = n3, n3 = t;
				A[i][j]->neibor[0] = n1;
				A[i][j]->neibor[1] = n2;
				A[i][j]->neibor[2] = n3;
				A[i][j]->neibor[3] = n4;
			}
			if (A[i][j]->neibor.size() > 4)
			{
				mesh* n1 = A[i][j]->neibor[0];
				mesh* n2 = A[i][j]->neibor[1];
				mesh* n3 = A[i][j]->neibor[2];
				mesh* n4 = A[i][j]->neibor[3];
				double maxy = max(max(max(n1->y, n2->y), n3->y), n4->y);
				if (n1->y == maxy)
					t = n1, n1 = n2, n2 = t;
				else if (n3->y == maxy)
					t = n3, n3 = n2, n2 = t;
				else if (n4->y == maxy)
					t = n4, n4 = n2, n2 = t;
				double miny = min(min(n1->y, n3->y), n4->y);
				if (n1->y == miny)
					t = n1, n1 = n4, n4 = t;
				else if (n3->y == miny)
					t = n3, n3 = n4, n4 = t;
				maxx = max(n1->x, n3->x);
				if (n3->x == maxx)
					t = n1, n1 = n3, n3 = t;
				A[i][j]->neibor[0] = n1;
				A[i][j]->neibor[1] = n2;
				A[i][j]->neibor[2] = n3;
				A[i][j]->neibor[3] = n4;
			}

		}

	}
}
double area(mesh A, mesh B, mesh C, mesh D)//求任意四点构成四边形面积 
{
	return 0.5 * abs(A.x * B.y + B.x * C.y + C.x * D.y + D.x * A.y - B.x * A.y - C.x * B.y - D.x * C.y - A.x * D.y);
}

void movemesh()
{
	extern double dt;
	extern vector<vector<mesh*>> A;
	int i, j;
#pragma omp parallel

	for (i = 0; i < A.size(); i++)
	{
		for (j = 0; j < A[i].size(); j++)
		{
			A[i][j]->x = A[i][j]->x + A[i][j]->um * dt;
			A[i][j]->y = A[i][j]->y + A[i][j]->vm * dt;
		}
	}
	//extern double t_sim;
	//extern int step;
	//extern vector <mesh> A0;
	//double L0 = A0[Xnum - 1].x0 - A0[0].x0;
	//double L1 = L0 / 2;
	////double L = L1 + 0.5*L1*sin(20 * pi*t_sim);
	//double L = L1 + 0.5 * L1 * sin(pi * step / 2000);
	//for (int i = 0; i < A0.size(); i++)
	//{
	//	if (A0[i].x0 < L1)
	//		A0[i].x = A0[i].x0 * (L / L1);
	//	else
	//		A0[i].x = A0[Xnum - 1].x0 - (A0[Xnum - 1].x0 - A0[i].x0) * ((L0 - L) / L1);
	//}
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
void findNeiborSec()
{
	extern vector<vector<mesh*>> A;
	int i, j, k, m;
#pragma omp parallel

	for (i = 0; i < A.size(); i++)
	{
		if (A.size() == 1)
			continue;
		for (j = 0; j < A[i].size(); j++)
		{
			for (k = i + 1; k < A.size(); k++)
			{
				for (m = 0; m < A[k].size(); m++)
				{
					if (A[i][j]->x == A[k][m]->x && A[i][j]->y == A[k][m]->y)
					{
						A[i][j]->neiborsec = k;
						A[i][j]->neiborsec_ad = m;
						A[k][m]->neiborsec = i;
						A[k][m]->neiborsec_ad = j;
					}
				}
			}

		}
	}
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
