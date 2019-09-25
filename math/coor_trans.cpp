#include"include/const.h"
#include"include/functions.h"
#include<vector>

void coordinate_trans()
{
	extern std::vector<std::vector<mesh*>> A;
	int i, j,k;
	int n1, n2, n3, n4;
	extern double t_sim;
	extern double dt;
	extern  std::vector <mesh> Ar;

	for (i = 0; i < A.size(); i++)
	{
		for (j = 0; j < A[i].size(); j++)
		{

			if (A[i][j]->neibor.size() < 3)//小于3的点无法做坐标变换
				continue;
			if (t_sim == 0)
			{
				if (A[i][j]->neibor.size() == 3)
				{
					using MeshPara::method;
					for (k = 0; k < 12; k++)
					{
						n1 = method[k][0];
						n2 = method[k][1];
						n3 = method[k][2];
						n4 = method[k][3];
						A[i][j]->xxi.push_back(0.5 * (A[i][j]->neibor[n1]->x - A[i][j]->neibor[n3]->x));
						A[i][j]->yxi.push_back(0.5 * (A[i][j]->neibor[n1]->y - A[i][j]->neibor[n3]->y));
						A[i][j]->xeta.push_back(0.5 * (A[i][j]->neibor[n2]->x - A[i][j]->neibor[n4]->x));
						A[i][j]->yeta.push_back(0.5 * (A[i][j]->neibor[n2]->y - A[i][j]->neibor[n4]->y));
						A[i][j]->J.push_back(1 / (A[i][j]->xxi[k] * A[i][j]->yeta[k] - A[i][j]->xeta[k] * A[i][j]->yxi[k]));
						A[i][j]->xix.push_back(A[i][j]->yeta[k] * A[i][j]->J[k]);
						A[i][j]->xiy.push_back(-A[i][j]->xeta[k] * A[i][j]->J[k]);
						A[i][j]->etax.push_back(-A[i][j]->yxi[k] * A[i][j]->J[k]);
						A[i][j]->etay.push_back(A[i][j]->xxi[k] * A[i][j]->J[k]);
						A[i][j]->xtau.push_back(0);
						A[i][j]->ytau.push_back(0);
						A[i][j]->xit.push_back(0);
						A[i][j]->etat.push_back(0);
					}
				}
				else if (A[i][j]->neibor.size() >= 4)
				{
					A[i][j]->xxi.push_back(0.5 * (A[i][j]->neibor[0]->x - A[i][j]->neibor[2]->x));
					A[i][j]->yxi.push_back(0.5 * (A[i][j]->neibor[0]->y - A[i][j]->neibor[2]->y));
					A[i][j]->xeta.push_back(0.5 * (A[i][j]->neibor[1]->x - A[i][j]->neibor[3]->x));
					A[i][j]->yeta.push_back(0.5 * (A[i][j]->neibor[1]->y - A[i][j]->neibor[3]->y));
					A[i][j]->J.push_back(1 / (A[i][j]->xxi[0] * A[i][j]->yeta[0] - A[i][j]->xeta[0] * A[i][j]->yxi[0]));
					A[i][j]->xix.push_back(A[i][j]->yeta[0] * A[i][j]->J[0]);
					A[i][j]->xiy.push_back(-A[i][j]->xeta[0] * A[i][j]->J[0]);
					A[i][j]->etax.push_back(-A[i][j]->yxi[0] * A[i][j]->J[0]);
					A[i][j]->etay.push_back(A[i][j]->xxi[0] * A[i][j]->J[0]);
					A[i][j]->xtau.push_back(0);
					A[i][j]->ytau.push_back(0);
					A[i][j]->xit.push_back(0);
					A[i][j]->etat.push_back(0);
				}
			}
			else//t_sim不等于零时还需要坐标变换，应为动网格！
			{
				if (A[i][j]->neibor.size() >= 4)
				{
					A[i][j]->xxi[0] = 0.5 * (A[i][j]->neibor[0]->x - A[i][j]->neibor[2]->x);
					A[i][j]->yxi[0] = 0.5 * (A[i][j]->neibor[0]->y - A[i][j]->neibor[2]->y);
					A[i][j]->xeta[0] = 0.5 * (A[i][j]->neibor[1]->x - A[i][j]->neibor[3]->x);
					A[i][j]->yeta[0] = 0.5 * (A[i][j]->neibor[1]->y - A[i][j]->neibor[3]->y);
					A[i][j]->J[0] = 1 / (A[i][j]->xxi[0] * A[i][j]->yeta[0] - A[i][j]->xeta[0] * A[i][j]->yxi[0]);
					A[i][j]->xix[0] = A[i][j]->yeta[0] * A[i][j]->J[0];
					A[i][j]->xiy[0] = -A[i][j]->xeta[0] * A[i][j]->J[0];
					A[i][j]->etax[0] = -A[i][j]->yxi[0] * A[i][j]->J[0];
					A[i][j]->etay[0] = A[i][j]->xxi[0] * A[i][j]->J[0];
					A[i][j]->xtau[0] = (A[i][j]->x - Ar[A[i][j]->id].x) / dt;
					A[i][j]->ytau[0] = (A[i][j]->y - Ar[A[i][j]->id].y) / dt;
					A[i][j]->xit[0] = -A[i][j]->xtau[0] * A[i][j]->xix[0] - A[i][j]->ytau[0] * A[i][j]->xiy[0];
					A[i][j]->etat[0] = -A[i][j]->xtau[0] * A[i][j]->etax[0] - A[i][j]->ytau[0] * A[i][j]->etay[0];
				}
				if (A[i][j]->neibor.size() == 3)
				{
					A[i][j]->xxi[0] = 0.5 * (A[i][j]->neibor[0]->x - A[i][j]->neibor[2]->x);
					A[i][j]->yxi[0] = 0.5 * (A[i][j]->neibor[0]->y - A[i][j]->neibor[2]->y);
					A[i][j]->xeta[0] = 0.5 * (A[i][j]->neibor[1]->x - A[i][j]->neibor[2]->x);
					A[i][j]->yeta[0] = 0.5 * (A[i][j]->neibor[1]->y - A[i][j]->neibor[2]->y);
					A[i][j]->J[0] = 1 / (A[i][j]->xxi[0] * A[i][j]->yeta[0] - A[i][j]->xeta[0] * A[i][j]->yxi[0]);
					A[i][j]->xix[0] = A[i][j]->yeta[0] * A[i][j]->J[0];
					A[i][j]->xiy[0] = -A[i][j]->xeta[0] * A[i][j]->J[0];
					A[i][j]->etax[0] = -A[i][j]->yxi[0] * A[i][j]->J[0];
					A[i][j]->etay[0] = A[i][j]->xxi[0] * A[i][j]->J[0];
					A[i][j]->xtau[0] = (A[i][j]->x - Ar[A[i][j]->id].x) / dt;
					A[i][j]->ytau[0] = (A[i][j]->y - Ar[A[i][j]->id].y) / dt;
					A[i][j]->xit[0] = -A[i][j]->xtau[0] * A[i][j]->xix[0] - A[i][j]->ytau[0] * A[i][j]->xiy[0];
					A[i][j]->etat[0] = -A[i][j]->xtau[0] * A[i][j]->etax[0] - A[i][j]->ytau[0] * A[i][j]->etay[0];

					A[i][j]->xxi[1] = 0.5 * (A[i][j]->neibor[1]->x - A[i][j]->neibor[0]->x);
					A[i][j]->yxi[1] = 0.5 * (A[i][j]->neibor[1]->y - A[i][j]->neibor[0]->y);
					A[i][j]->xeta[1] = 0.5 * (A[i][j]->neibor[2]->x - A[i][j]->neibor[0]->x);
					A[i][j]->yeta[1] = 0.5 * (A[i][j]->neibor[2]->y - A[i][j]->neibor[0]->y);
					A[i][j]->J[1] = 1 / (A[i][j]->xxi[1] * A[i][j]->yeta[1] - A[i][j]->xeta[1] * A[i][j]->yxi[1]);
					A[i][j]->xix[1] = A[i][j]->yeta[1] * A[i][j]->J[1];
					A[i][j]->xiy[1] = -A[i][j]->xeta[1] * A[i][j]->J[1];
					A[i][j]->etax[1] = -A[i][j]->yxi[1] * A[i][j]->J[1];
					A[i][j]->etay[1] = A[i][j]->xxi[1] * A[i][j]->J[1];
					A[i][j]->xtau[1] = (A[i][j]->x - Ar[A[i][j]->id].x) / dt;
					A[i][j]->ytau[1] = (A[i][j]->y - Ar[A[i][j]->id].y) / dt;
					A[i][j]->xit[1] = -A[i][j]->xtau[1] * A[i][j]->xix[1] - A[i][j]->ytau[1] * A[i][j]->xiy[1];
					A[i][j]->etat[1] = -A[i][j]->xtau[1] * A[i][j]->etax[1] - A[i][j]->ytau[1] * A[i][j]->etay[1];

					A[i][j]->xxi[2] = 0.5 * (A[i][j]->neibor[2]->x - A[i][j]->neibor[1]->x);
					A[i][j]->yxi[2] = 0.5 * (A[i][j]->neibor[2]->y - A[i][j]->neibor[1]->y);
					A[i][j]->xeta[2] = 0.5 * (A[i][j]->neibor[0]->x - A[i][j]->neibor[1]->x);
					A[i][j]->yeta[2] = 0.5 * (A[i][j]->neibor[0]->y - A[i][j]->neibor[1]->y);
					A[i][j]->J[2] = 1 / (A[i][j]->xxi[2] * A[i][j]->yeta[2] - A[i][j]->xeta[2] * A[i][j]->yxi[2]);
					A[i][j]->xix[2] = A[i][j]->yeta[2] * A[i][j]->J[2];
					A[i][j]->xiy[2] = -A[i][j]->xeta[2] * A[i][j]->J[2];
					A[i][j]->etax[2] = -A[i][j]->yxi[2] * A[i][j]->J[2];
					A[i][j]->etay[2] = A[i][j]->xxi[2] * A[i][j]->J[2];
					A[i][j]->xtau[2] = (A[i][j]->x - Ar[A[i][j]->id].x) / dt;
					A[i][j]->ytau[2] = (A[i][j]->y - Ar[A[i][j]->id].y) / dt;
					A[i][j]->xit[2] = -A[i][j]->xtau[2] * A[i][j]->xix[2] - A[i][j]->ytau[2] * A[i][j]->xiy[2];
					A[i][j]->etat[2] = -A[i][j]->xtau[2] * A[i][j]->etax[2] - A[i][j]->ytau[2] * A[i][j]->etay[2];

				}

			}
		}
	}
}
