#include"include/const.h"
#include"include/functions.h"
#include<omp.h>
#include"include/shockwave.h"
#include<iostream>
#include"include/Prandtl-Meyer.h"
using std::vector;
using namespace ConstPara;
using namespace MeshPara;

void get_dt()
{
	extern vector<mesh>AP;
	double maxxi = 0, maxeta = 0;
	extern double dt;
	double t;
	int i, j, k;
	double max1, max2;
	double Sxi, Seta, c, uxi, ueta;

	dt = t_end;
	max1 = max2 = 0;
	for (i = 0; i < AP.size(); i++)
	{
		maxxi = maxeta = 0;
		if (AP[i].type != "IN")
			continue;
		if (AP[i].neibor.size() > 3)
		{
			Sxi = sqrt(AP[i].xix[0] * AP[i].xix[0] + AP[i].xiy[0] * AP[i].xiy[0]);
			Seta = sqrt(AP[i].etax[0] * AP[i].etax[0] + AP[i].etay[0] * AP[i].etay[0]);
			c = sqrt(gama * AP[i].p / AP[i].rho);
			uxi = AP[i].u * AP[i].xix[0] + AP[i].v * AP[i].xiy[0];
			ueta = AP[i].u * AP[i].etax[0] + AP[i].v * AP[i].etay[0];
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
				Sxi = sqrt(AP[i].xix[k] * AP[i].xix[k] + AP[i].xiy[k] * AP[i].xiy[k]);
				Seta = sqrt(AP[i].etax[k] * AP[i].etax[k] + AP[i].etay[k] * AP[i].etay[k]);
				c = sqrt(gama * AP[i].p / AP[i].rho);
				uxi = AP[i].u * AP[i].xix[k] + AP[i].v * AP[i].xiy[k];
				ueta = AP[i].u * AP[i].etax[k] + AP[i].v * AP[i].etay[k];
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
	t = CFL / (max1 + max2);
	//t = CFL / (maxxi + maxeta);
	dt = min(dt, t);
}
void update_Vm()
{
	extern vector<vector<mesh*>> A;
	extern vector<mesh>AP;
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
					A[i][j]->um = (AP[Xnum - 1].x - A[i][j]->x) / (AP[Xnum - 1].x - x) * um;
					A[i][j]->vm = (AP[Xnum - 1].x - A[i][j]->x) / (AP[Xnum - 1].x - x) * vm;
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
	extern vector <mesh> AP;
	extern vector <mesh> Ar;
	extern double t_sim;
	int i, j;
	vector<mesh>t;
	if (t_sim == 0)

		for (i = 0; i < AP.size(); i++)
		{
			Ar.push_back(AP[i]);
		}
	else

#pragma omp parallel for
		for (i = 0; i < AP.size(); i++)
		{
			Ar[i] = AP[i];
		}

}

void update_p3(mesh& p)
//unstructral grid point,3 neighbor points
{
	extern vector <mesh> Ar;
	extern double dt;
	int i, j;
	int n1, n2, n3, n4;
	double U[4], U1[4], U2[4];

	Flux Fll, Flr, Fcl, Fcr, Frl, Frr;
	Flux Gdd, Gdu, Gcd, Gcu, Gud, Guu;

	int id;
	id = p.id;
	U[0] = Ar[id].rho;
	U[1] = Ar[id].rho * p.u;
	U[2] = Ar[id].rho * p.v;
	U[3] = 0.5 * Ar[id].rho * (Ar[id].u * Ar[id].u + Ar[id].v * Ar[id].v) + Ar[id].p / (gama - 1);

	for (j = 0; j < 12; j++)
	{
		n1 = method[j][0];
		n2 = method[j][1];
		n3 = method[j][2];
		n4 = method[j][3];
		n1 = p.neibor[n1]->id;
		n2 = p.neibor[n2]->id;
		n3 = p.neibor[n3]->id;
		n4 = p.neibor[n4]->id;
		Fll = Flr = Fcl = Fcr = Frl = Frr = { 0 };
		Gdd = Gdu = Gcd = Gcu = Gud = Guu = { 0 };

		Fll = Fll + VanLeerB(Ar[n3], Ar[id].xix[j], Ar[id].xiy[j], Ar[id].xit[j], Ar[id].J[j]) * p.J[j];
		Flr = Flr + VanLeerA(Ar[n3], Ar[id].xix[j], Ar[id].xiy[j], Ar[id].xit[j], Ar[id].J[j]) * p.J[j];
		Fcl = Fcl + VanLeerB(Ar[id], Ar[id].xix[j], Ar[id].xiy[j], Ar[id].xit[j], Ar[id].J[j]) * p.J[j];
		Fcr = Fcr + VanLeerA(Ar[id], Ar[id].xix[j], Ar[id].xiy[j], Ar[id].xit[j], Ar[id].J[j]) * p.J[j];
		Frl = Frl + VanLeerB(Ar[n1], Ar[id].xix[j], Ar[id].xiy[j], Ar[id].xit[j], Ar[id].J[j]) * p.J[j];
		Frr = Frr + VanLeerA(Ar[n1], Ar[id].xix[j], Ar[id].xiy[j], Ar[id].xit[j], Ar[id].J[j]) * p.J[j];

		Gdd = Gdd + VanLeerB(Ar[n4], Ar[id].etax[j], Ar[id].etay[j], Ar[id].etat[j], Ar[id].J[j]) * p.J[j];
		Gdu = Gdu + VanLeerA(Ar[n4], Ar[id].etax[j], Ar[id].etay[j], Ar[id].etat[j], Ar[id].J[j]) * p.J[j];
		Gcd = Gcd + VanLeerB(Ar[id], Ar[id].etax[j], Ar[id].etay[j], Ar[id].etat[j], Ar[id].J[j]) * p.J[j];
		Gcu = Gcu + VanLeerA(Ar[id], Ar[id].etax[j], Ar[id].etay[j], Ar[id].etat[j], Ar[id].J[j]) * p.J[j];
		Gud = Gud + VanLeerB(Ar[n2], Ar[id].etax[j], Ar[id].etay[j], Ar[id].etat[j], Ar[id].J[j]) * p.J[j];
		Guu = Guu + VanLeerA(Ar[n2], Ar[id].etax[j], Ar[id].etay[j], Ar[id].etat[j], Ar[id].J[j]) * p.J[j];
	}
	Fll = Fll / j, Flr = Flr / j, Fcl = Fcl / j, Fcr = Fcr / j, Frl = Frl / j, Frr = Frr / j;
	Gdd = Gdd / j, Gdu = Gdu / j, Gcd = Gcd / j, Gcu = Gcu / j, Gud = Gud / j, Guu = Guu / j;

	U[0] = U[0] - dt * (Fcr.f1 - Flr.f1 + Frl.f1 - Fcl.f1 + Gcu.f1 - Gdu.f1 + Gud.f1 - Gcd.f1);
	U[1] = U[1] - dt * (Fcr.f2 - Flr.f2 + Frl.f2 - Fcl.f2 + Gcu.f2 - Gdu.f2 + Gud.f2 - Gcd.f2);
	U[2] = U[2] - dt * (Fcr.f3 - Flr.f3 + Frl.f3 - Fcl.f3 + Gcu.f3 - Gdu.f3 + Gud.f3 - Gcd.f3);
	U[3] = U[3] - dt * (Fcr.f4 - Flr.f4 + Frl.f4 - Fcl.f4 + Gcu.f4 - Gdu.f4 + Gud.f4 - Gcd.f4);
	p.rho = U[0];
	p.u = U[1] / U[0];
	p.v = U[2] / U[0];
	p.p = (gama - 1) * (U[3] - 0.5 * p.rho * (p.u * p.u + p.v * p.v));
}
void update_p4_s(mesh& p)
//structral grid point,4 neighbor points
{
	extern vector <mesh> Ar;
	extern double dt;
	int i, j;
	int n1, n2, n3, n4;
	double U[4], U1[4], U2[4];

	Flux Fll, Flr, Fcl, Fcr, Frl, Frr;
	Flux Gdd, Gdu, Gcd, Gcu, Gud, Guu;

	int id;
	id = p.id;
	U[0] = Ar[id].rho;
	U[1] = Ar[id].rho * p.u;
	U[2] = Ar[id].rho * p.v;
	U[3] = 0.5 * Ar[id].rho * (Ar[id].u * Ar[id].u + Ar[id].v * Ar[id].v) + Ar[id].p / (gama - 1);

	n1 = p.neibor[0]->id;
	n2 = p.neibor[1]->id;
	n3 = p.neibor[2]->id;
	n4 = p.neibor[3]->id;


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

	U[0] = U[0] - dt * p.J[0] * (Fcr.f1 - Flr.f1 + Frl.f1 - Fcl.f1 + Gcu.f1 - Gdu.f1 + Gud.f1 - Gcd.f1);
	U[1] = U[1] - dt * p.J[0] * (Fcr.f2 - Flr.f2 + Frl.f2 - Fcl.f2 + Gcu.f2 - Gdu.f2 + Gud.f2 - Gcd.f2);
	U[2] = U[2] - dt * p.J[0] * (Fcr.f3 - Flr.f3 + Frl.f3 - Fcl.f3 + Gcu.f3 - Gdu.f3 + Gud.f3 - Gcd.f3);
	U[3] = U[3] - dt * p.J[0] * (Fcr.f4 - Flr.f4 + Frl.f4 - Fcl.f4 + Gcu.f4 - Gdu.f4 + Gud.f4 - Gcd.f4);
	p.rho = U[0] * p.sec_num;
	p.u = U[1] / U[0] * p.sec_num;
	p.v = U[2] / U[0] * p.sec_num;
	p.p = (gama - 1) * (U[3] - 0.5 * p.rho * (p.u * p.u + p.v * p.v)) * p.sec_num;
}
void update_p4_u(mesh& p)
//unstructral grid point,4 neighbor points
{
	extern vector <mesh> Ar;
	extern double dt;
	int i, j;
	int n1, n2, n3, n4;
	double U[4], U1[4], U2[4];

	Flux Fll, Flr, Fcl, Fcr, Frl, Frr;
	Flux Gdd, Gdu, Gcd, Gcu, Gud, Guu;

	int id;
	id = p.id;
	U[0] = Ar[id].rho;
	U[1] = Ar[id].rho * p.u;
	U[2] = Ar[id].rho * p.v;
	U[3] = 0.5 * Ar[id].rho * (Ar[id].u * Ar[id].u + Ar[id].v * Ar[id].v) + Ar[id].p / (gama - 1);

	n1 = p.neibor[0]->id;
	n2 = p.neibor[1]->id;
	n3 = p.neibor[2]->id;
	n4 = p.neibor[3]->id;
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

	U[0] = U[0] - dt * p.J[0] * (Fcr.f1 - Flr.f1 + Frl.f1 - Fcl.f1 + Gcu.f1 - Gdu.f1 + Gud.f1 - Gcd.f1);
	U[1] = U[1] - dt * p.J[0] * (Fcr.f2 - Flr.f2 + Frl.f2 - Fcl.f2 + Gcu.f2 - Gdu.f2 + Gud.f2 - Gcd.f2);
	U[2] = U[2] - dt * p.J[0] * (Fcr.f3 - Flr.f3 + Frl.f3 - Fcl.f3 + Gcu.f3 - Gdu.f3 + Gud.f3 - Gcd.f3);
	U[3] = U[3] - dt * p.J[0] * (Fcr.f4 - Flr.f4 + Frl.f4 - Fcl.f4 + Gcu.f4 - Gdu.f4 + Gud.f4 - Gcd.f4);
	p.rho = U[0];
	p.u = U[1] / U[0];
	p.v = U[2] / U[0];
	p.p = (gama - 1) * (U[3] - 0.5 * p.rho * (p.u * p.u + p.v * p.v));
}
void update_bound()
{
	extern vector<mesh*> bl;//左边界
	extern vector<mesh*> br;//右边界
	extern vector<mesh*> bu;//上边界
	extern vector<mesh*> bd;//下边界
	extern vector<mesh*> bb;//物体边界
	int i;
	using namespace Init;
	extern vector<mesh>AP;
	int id;

	for (i = 0; i < bl.size(); i++)
	{
		bl[i]->rho = rho0;
		bl[i]->u = 12 * sqrt(ConstPara::gama * p0 / rho0);
		bl[i]->v = v0;
		bl[i]->p = p0;
	}

	for (i = 0; i < br.size(); i++)
	{
		id = br[i]->id;
		br[i]->rho = AP[id - 1].rho;
		br[i]->u = AP[id - 1].u;
		br[i]->v = AP[id - 1].v;
		br[i]->p = AP[id - 1].p;
	}
	for (i = 0; i < bu.size(); i++)
	{
		id = bu[i]->id;
		bu[i]->rho = AP[id - Xnum].rho;
		bu[i]->u = AP[id - Xnum].u;
		bu[i]->v = AP[id - Xnum].v;
		bu[i]->p = AP[id - Xnum].p;
	}

	for (i = 0; i < bd.size(); i++)
	{
		id = bd[i]->id;
		bd[i]->rho = AP[id + Xnum].rho;
		bd[i]->u = AP[id + Xnum].u;
		bd[i]->v = AP[id + Xnum].v;
		bd[i]->p = AP[id + Xnum].p;
	}
	for (i = 0; i < bb.size(); i++)
	{
		double DELTA = 1e-10;
		double nx, ny, tx, ty, n1, n2;
		nx = bb[i]->neibor[0]->x - bb[i]->x;
		ny = bb[i]->neibor[0]->y - bb[i]->y;
		tx = -ny, ty = nx;
		//if (bb[i]->x < a)
		//{
		//	nx = bb[i]->x - a;
		//	ny = bb[i]->y - b;
		//	if (ny > 0)
		//	{
		//		tx = ny;
		//		ty = -nx;
		//	}
		//	else
		//	{
		//		tx = -ny;
		//		ty = nx;
		//	}
		//}
		//else if (bb[i]->y > b)
		//{
		//	double k = tan(4.6 * ConstPara::pi / 180);
		//	nx = -k, ny = 1;
		//	tx = 1, ty = k;
		//}
		//else
		//{
		//	double k = tan(-4.6 * ConstPara::pi / 180);
		//	nx = -k, ny = 1;
		//	tx = 1, ty = k;
		//}
		n1 = bb[i]->neibor[0]->u;
		n2 = bb[i]->neibor[0]->v;
		if ((abs(n1) < DELTA && abs(n2) < DELTA)/* || (abs(tx) < DELTA || abs(ty) < DELTA)*/)
		{
			bb[i]->rho = bb[i]->neibor[0]->rho;
			bb[i]->u = bb[i]->neibor[0]->u;
			bb[i]->v = bb[i]->neibor[0]->v;
			bb[i]->p = bb[i]->neibor[0]->p;
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
			bb[i]->rho = bb[i]->neibor[0]->rho;
			bb[i]->u = u;
			bb[i]->v = v;
			bb[i]->p = bb[i]->neibor[0]->p;
		}
		//bb[i]->rho = bb[i]->neibor[0]->rho;
		//bb[i]->u = bb[i]->neibor[0]->u;
		//bb[i]->v = bb[i]->neibor[0]->v;
		//bb[i]->p = bb[i]->neibor[0]->p;
		//bb[i]->rho = rho0;
		//bb[i]->u = 0;
		//bb[i]->v = 0;
		//bb[i]->p = p0;

	}

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
	//extern vector <mesh> AP;
	//double L0 = AP[Xnum - 1].x0 - AP[0].x0;
	//double L1 = L0 / 2;
	////double L = L1 + 0.5*L1*sin(20 * pi*t_sim);
	//double L = L1 + 0.5 * L1 * sin(pi * step / 2000);
	//for (int i = 0; i < AP.size(); i++)
	//{
	//	if (AP[i].x0 < L1)
	//		AP[i].x = AP[i].x0 * (L / L1);
	//	else
	//		AP[i].x = AP[Xnum - 1].x0 - (AP[Xnum - 1].x0 - AP[i].x0) * ((L0 - L) / L1);
	//}
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

void sortPoint()
//put mesh point into different arrays
{
	extern vector<mesh*> ps;//结构网格节点
	extern vector<mesh*> pu;//非结构网格节点
	extern vector<mesh*> bl;//左边界
	extern vector<mesh*> br;//右边界
	extern vector<mesh*> bu;//上边界
	extern vector<mesh*> bd;//下边界
	extern vector<mesh*> bb;//物体边界
	extern vector <mesh> AP;
	int n;
	for (int i = 0; i < AP.size(); i++)
	{
		if (AP[i].type == "IN")
		{
			n = 0;
			for (int j = 0; j < AP[i].neibor.size(); j++)
				if (AP[i].neibor[j]->type=="Body")
					n++;
			if(n==0&& AP[i].neibor.size()==4)
				ps.push_back(&AP[i]);
			else
				pu.push_back(&AP[i]);
		}
			//if (AP[i].neibor.size() == 4 )
			//	for
			//	ps.push_back(&AP[i]);
			//else if (AP[i].neibor.size() == 3 || AP[i].neibor.size() == 4)
			//	pu.push_back(&AP[i]);
			//else
			//	std::cout << "not a ps or pu inner point" << std::endl;
		else if (AP[i].type == "L")
			bl.push_back(&AP[i]);
		else if (AP[i].type == "R")
			br.push_back(&AP[i]);
		else if (AP[i].type == "U")
			bu.push_back(&AP[i]);
		else if (AP[i].type == "D")
			bd.push_back(&AP[i]);
		else if (AP[i].type == "Body")
			bb.push_back(&AP[i]);
		else
			std::cout << "need to confirm the point's type of ID" << AP[i].id << std::endl;
	}
}
void polymesh()
//get polygon mesh from grid points
{
	extern vector<mesh> AP;
	extern vector<polygon_mesh>PM;
	polygon_mesh Ptemp;
	int i, j;
	double DELTA = 1e-10;
	int n1, n2, n3, n4, n5;
	int size;
	int id, id1, id2;
	for (i = 0; i < Xnum * Ynum; i++)
	{
		if (i + Xnum + 1 < Xnum * Ynum && abs(AP[i + 1].x - AP[i].x - dx) < DELTA)
		{
			n1 = i;
			n2 = i + 1;
			n3 = i + Xnum + 1;
			n4 = i + Xnum;
			//if only one point is out of body
			if ((AP[n1].section != 0 && AP[n2].section == 0 && AP[n3].section == 0 && AP[n4].section == 0) ||
				(AP[n1].section == 0 && AP[n2].section != 0 && AP[n3].section == 0 && AP[n4].section == 0) ||
				(AP[n1].section == 0 && AP[n2].section == 0 && AP[n3].section != 0 && AP[n4].section == 0) ||
				(AP[n1].section == 0 && AP[n2].section == 0 && AP[n3].section == 0 && AP[n4].section != 0) ||
				(AP[n1].section == 0 && AP[n2].section == 0 && AP[n3].section == 0 && AP[n4].section == 0))
				continue;
			//if all four points are out of body
			else if (AP[n1].section != 0 && AP[n2].section != 0 && AP[n3].section != 0 && AP[n4].section != 0)
			{
				PM.push_back(Ptemp);
				size = PM.size() - 1;
				PM[size].node.push_back(n1);
				PM[size].node.push_back(n2);
				PM[size].node.push_back(n3);
				PM[size].node.push_back(n4);
				PM[size].face_start.push_back(n1);
				PM[size].face_start.push_back(n2);
				PM[size].face_start.push_back(n3);
				PM[size].face_start.push_back(n4);
				PM[size].face_end.push_back(n2);
				PM[size].face_end.push_back(n3);
				PM[size].face_end.push_back(n4);
				PM[size].face_end.push_back(n1);
			}
			else if (AP[n1].section == 0 && AP[n2].section != 0 && AP[n3].section != 0 && AP[n4].section != 0)
			{
				id = AP[n2].connectId;
				for (j = 0; j < AP[id].neibor.size(); j++)
				{
					if (AP[id].neibor[j]->type == "Body")
						id1 = AP[id].neibor[j]->id;
				}
				id = AP[n4].connectId;
				for (j = 0; j < AP[id].neibor.size(); j++)
				{
					if (AP[id].neibor[j]->type == "Body")
						id2 = AP[id].neibor[j]->id;
				}
				if (id1 == id2)
				{
					n1 = id1;
					PM.push_back(Ptemp);
					size = PM.size() - 1;
					PM[size].node.push_back(n1);
					PM[size].node.push_back(n2);
					PM[size].node.push_back(n3);
					PM[size].node.push_back(n4);
					PM[size].face_start.push_back(n1);
					PM[size].face_start.push_back(n2);
					PM[size].face_start.push_back(n3);
					PM[size].face_start.push_back(n4);
					PM[size].face_end.push_back(n2);
					PM[size].face_end.push_back(n3);
					PM[size].face_end.push_back(n4);
					PM[size].face_end.push_back(n1);
				}
				else
				{
					n1 = id1;
					n5 = id2;
					//judge whether n1 and n5 are already be marked as neighbors
					int n = 0;
					for (j = 0; j < AP[n1].neibor.size(); j++)
					{
						if (AP[n5].id == AP[n1].neibor[j]->id)
						{
							n++;
							break;
						}
					}
					if (n == 0)
						AP[n1].neibor.push_back(&AP[n5]);
					n = 0;
					for (j = 0; j < AP[n5].neibor.size(); j++)
					{
						if (AP[n1].id == AP[n5].neibor[j]->id)
						{
							n++;
							break;
						}
					}
					if (n == 0)
						AP[n5].neibor.push_back(&AP[n1]);
					PM.push_back(Ptemp);
					size = PM.size() - 1;
					PM[size].node.push_back(n1);
					PM[size].node.push_back(n2);
					PM[size].node.push_back(n3);
					PM[size].node.push_back(n4);
					PM[size].node.push_back(n5);
					PM[size].face_start.push_back(n1);
					PM[size].face_start.push_back(n2);
					PM[size].face_start.push_back(n3);
					PM[size].face_start.push_back(n4);
					PM[size].face_start.push_back(n5);
					PM[size].face_end.push_back(n2);
					PM[size].face_end.push_back(n3);
					PM[size].face_end.push_back(n4);
					PM[size].face_end.push_back(n5);
					PM[size].face_end.push_back(n1);
				}
			}
			else if (AP[n1].section != 0 && AP[n2].section == 0 && AP[n3].section != 0 && AP[n4].section != 0)
			{
				id = AP[n3].connectId;
				for (j = 0; j < AP[id].neibor.size(); j++)
				{
					if (AP[id].neibor[j]->type == "Body")
						id1 = AP[id].neibor[j]->id;
				}
				id = AP[n1].connectId;
				for (j = 0; j < AP[id].neibor.size(); j++)
				{
					if (AP[id].neibor[j]->type == "Body")
						id2 = AP[id].neibor[j]->id;
				}
				if (id1 == id2)
				{
					n2 = id1;
					PM.push_back(Ptemp);
					size = PM.size() - 1;
					PM[size].node.push_back(n1);
					PM[size].node.push_back(n2);
					PM[size].node.push_back(n3);
					PM[size].node.push_back(n4);
					PM[size].face_start.push_back(n1);
					PM[size].face_start.push_back(n2);
					PM[size].face_start.push_back(n3);
					PM[size].face_start.push_back(n4);
					PM[size].face_end.push_back(n2);
					PM[size].face_end.push_back(n3);
					PM[size].face_end.push_back(n4);
					PM[size].face_end.push_back(n1);
				}
				else
				{
					n2 = id1;
					n5 = id2;
					//judge whether n1 and n5 are already be marked as neighbors
					int n = 0;
					for (j = 0; j < AP[n2].neibor.size(); j++)
					{
						if (AP[n5].id == AP[n2].neibor[j]->id)
						{
							n++;
							break;
						}
					}
					if (n == 0)
						AP[n2].neibor.push_back(&AP[n5]);
					n = 0;
					for (j = 0; j < AP[n5].neibor.size(); j++)
					{
						if (AP[n2].id == AP[n5].neibor[j]->id)
						{
							n++;
							break;
						}
					}
					if (n == 0)
						AP[n5].neibor.push_back(&AP[n2]);
					PM.push_back(Ptemp);
					size = PM.size() - 1;
					PM[size].node.push_back(n1);
					PM[size].node.push_back(n5);
					PM[size].node.push_back(n2);
					PM[size].node.push_back(n3);
					PM[size].node.push_back(n4);
					PM[size].face_start.push_back(n1);
					PM[size].face_start.push_back(n5);
					PM[size].face_start.push_back(n2);
					PM[size].face_start.push_back(n3);
					PM[size].face_start.push_back(n4);
					PM[size].face_end.push_back(n5);
					PM[size].face_end.push_back(n2);
					PM[size].face_end.push_back(n3);
					PM[size].face_end.push_back(n4);
					PM[size].face_end.push_back(n1);
				}
			}
			else if (AP[n1].section != 0 && AP[n2].section != 0 && AP[n3].section == 0 && AP[n4].section != 0)
			{
				id = AP[n4].connectId;
				for (j = 0; j < AP[id].neibor.size(); j++)
				{
					if (AP[id].neibor[j]->type == "Body")
						id1 = AP[id].neibor[j]->id;
				}
				id = AP[n2].connectId;
				for (j = 0; j < AP[id].neibor.size(); j++)
				{
					if (AP[id].neibor[j]->type == "Body")
						id2 = AP[id].neibor[j]->id;
				}
				if (id1 == id2)
				{
					n3 = id1;
					PM.push_back(Ptemp);
					size = PM.size() - 1;
					PM[size].node.push_back(n1);
					PM[size].node.push_back(n2);
					PM[size].node.push_back(n3);
					PM[size].node.push_back(n4);
					PM[size].face_start.push_back(n1);
					PM[size].face_start.push_back(n2);
					PM[size].face_start.push_back(n3);
					PM[size].face_start.push_back(n4);
					PM[size].face_end.push_back(n2);
					PM[size].face_end.push_back(n3);
					PM[size].face_end.push_back(n4);
					PM[size].face_end.push_back(n1);
				}
				else
				{
					n3 = id1;
					n5 = id2;
					//judge whether n1 and n5 are already be marked as neighbors
					int n = 0;
					for (j = 0; j < AP[n3].neibor.size(); j++)
					{
						if (AP[n5].id == AP[n3].neibor[j]->id)
						{
							n++;
							break;
						}
					}
					if (n == 0)
						AP[n3].neibor.push_back(&AP[n5]);
					n = 0;
					for (j = 0; j < AP[n5].neibor.size(); j++)
					{
						if (AP[n3].id == AP[n5].neibor[j]->id)
						{
							n++;
							break;
						}
					}
					if (n == 0)
						AP[n5].neibor.push_back(&AP[n3]);
					PM.push_back(Ptemp);
					size = PM.size() - 1;
					PM[size].node.push_back(n1);
					PM[size].node.push_back(n2);
					PM[size].node.push_back(n5);
					PM[size].node.push_back(n3);
					PM[size].node.push_back(n4);
					PM[size].face_start.push_back(n1);
					PM[size].face_start.push_back(n2);
					PM[size].face_start.push_back(n5);
					PM[size].face_start.push_back(n3);
					PM[size].face_start.push_back(n4);
					PM[size].face_end.push_back(n2);
					PM[size].face_end.push_back(n5);
					PM[size].face_end.push_back(n3);
					PM[size].face_end.push_back(n4);
					PM[size].face_end.push_back(n1);
				}
			}
			else if (AP[n1].section != 0 && AP[n2].section != 0 && AP[n3].section != 0 && AP[n4].section == 0)
			{
				id = AP[n1].connectId;
				for (j = 0; j < AP[id].neibor.size(); j++)
				{
					if (AP[id].neibor[j]->type == "Body")
						id1 = AP[id].neibor[j]->id;
				}
				id = AP[n3].connectId;
				for (j = 0; j < AP[id].neibor.size(); j++)
				{
					if (AP[id].neibor[j]->type == "Body")
						id2 = AP[id].neibor[j]->id;
				}
				if (id1 == id2)
				{
					n4 = id1;
					PM.push_back(Ptemp);
					size = PM.size() - 1;
					PM[size].node.push_back(n1);
					PM[size].node.push_back(n2);
					PM[size].node.push_back(n3);
					PM[size].node.push_back(n4);
					PM[size].face_start.push_back(n1);
					PM[size].face_start.push_back(n2);
					PM[size].face_start.push_back(n3);
					PM[size].face_start.push_back(n4);
					PM[size].face_end.push_back(n2);
					PM[size].face_end.push_back(n3);
					PM[size].face_end.push_back(n4);
					PM[size].face_end.push_back(n1);
				}
				else
				{
					n4 = id1;
					n5 = id2;
					//judge whether n1 and n5 are already be marked as neighbors
					int n = 0;
					for (j = 0; j < AP[n4].neibor.size(); j++)
					{
						if (AP[n5].id == AP[n4].neibor[j]->id)
						{
							n++;
							break;
						}
					}
					if (n == 0)
						AP[n4].neibor.push_back(&AP[n5]);
					n = 0;
					for (j = 0; j < AP[n5].neibor.size(); j++)
					{
						if (AP[n4].id == AP[n5].neibor[j]->id)
						{
							n++;
							break;
						}
					}
					if (n == 0)
						AP[n5].neibor.push_back(&AP[n4]);
					PM.push_back(Ptemp);
					size = PM.size() - 1;
					PM[size].node.push_back(n1);
					PM[size].node.push_back(n2);
					PM[size].node.push_back(n3);
					PM[size].node.push_back(n5);
					PM[size].node.push_back(n4);
					PM[size].face_start.push_back(n1);
					PM[size].face_start.push_back(n2);
					PM[size].face_start.push_back(n3);
					PM[size].face_start.push_back(n5);
					PM[size].face_start.push_back(n4);
					PM[size].face_end.push_back(n2);
					PM[size].face_end.push_back(n3);
					PM[size].face_end.push_back(n5);
					PM[size].face_end.push_back(n4);
					PM[size].face_end.push_back(n1);
				}
			}
			else
			{
				if (AP[n1].section == 0 && AP[n2].section == 0 && AP[n3].section != 0 && AP[n4].section != 0)
				{
					id = AP[n3].connectId;
					for (j = 0; j < AP[id].neibor.size(); j++)
					{
						if (AP[id].neibor[j]->type == "Body")
							id1 = AP[id].neibor[j]->id;
					}
					id = AP[n4].connectId;
					for (j = 0; j < AP[id].neibor.size(); j++)
					{
						if (AP[id].neibor[j]->type == "Body")
							id2 = AP[id].neibor[j]->id;
					}
					n1 = id2;
					n2 = id1;
				}
				else if (AP[n1].section != 0 && AP[n2].section == 0 && AP[n3].section == 0 && AP[n4].section != 0)
				{
					id = AP[n4].connectId;
					for (j = 0; j < AP[id].neibor.size(); j++)
					{
						if (AP[id].neibor[j]->type == "Body")
							id1 = AP[id].neibor[j]->id;
					}
					id = AP[n1].connectId;
					for (j = 0; j < AP[id].neibor.size(); j++)
					{
						if (AP[id].neibor[j]->type == "Body")
							id2 = AP[id].neibor[j]->id;
					}
					n2 = id2;
					n3 = id1;
				}
				else if (AP[n1].section != 0 && AP[n2].section != 0 && AP[n3].section == 0 && AP[n4].section == 0)
				{
					id = AP[n1].connectId;
					for (j = 0; j < AP[id].neibor.size(); j++)
					{
						if (AP[id].neibor[j]->type == "Body")
							id1 = AP[id].neibor[j]->id;
					}
					id = AP[n2].connectId;
					for (j = 0; j < AP[id].neibor.size(); j++)
					{
						if (AP[id].neibor[j]->type == "Body")
							id2 = (*AP[id].neibor[j]).id;
					}
					n3 = id2;
					n4 = id1;
				}
				else if (AP[n1].section == 0 && AP[n2].section != 0 && AP[n3].section != 0 && AP[n4].section == 0)
				{
					id = AP[n2].connectId;
					for (j = 0; j < AP[id].neibor.size(); j++)
					{
						if (AP[id].neibor[j]->type == "Body")
							id1 = AP[id].neibor[j]->id;
					}
					id = AP[n3].connectId;
					for (j = 0; j < AP[id].neibor.size(); j++)
					{
						if (AP[id].neibor[j]->type == "Body")
							id2 = AP[id].neibor[j]->id;
					}
					n4 = id2;
					n1 = id1;
				}
				int n = 0;
				for (j = 0; j < AP[id1].neibor.size(); j++)
				{
					if (AP[id2].id == AP[id1].neibor[j]->id)
					{
						n++;
						break;
					}
				}
				if (n == 0)
					AP[id1].neibor.push_back(&AP[id2]);
				n = 0;
				for (j = 0; j < AP[id2].neibor.size(); j++)
				{
					if (AP[id1].id == AP[id2].neibor[j]->id)
					{
						n++;
						break;
					}
				}
				if (n == 0)
					AP[id2].neibor.push_back(&AP[id1]);
				PM.push_back(Ptemp);
				size = PM.size() - 1;
				PM[size].node.push_back(n1);
				PM[size].node.push_back(n2);
				PM[size].node.push_back(n3);
				PM[size].node.push_back(n4);
				PM[size].face_start.push_back(n1);
				PM[size].face_start.push_back(n2);
				PM[size].face_start.push_back(n3);
				PM[size].face_start.push_back(n4);
				PM[size].face_end.push_back(n2);
				PM[size].face_end.push_back(n3);
				PM[size].face_end.push_back(n4);
				PM[size].face_end.push_back(n1);

			}
		}
	}
}