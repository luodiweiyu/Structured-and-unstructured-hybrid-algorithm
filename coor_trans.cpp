#include"const.h"
#include"functions.h"
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
						A[i][j]->xξ.push_back(0.5 * (A[i][j]->neibor[n1]->x - A[i][j]->neibor[n3]->x));
						A[i][j]->yξ.push_back(0.5 * (A[i][j]->neibor[n1]->y - A[i][j]->neibor[n3]->y));
						A[i][j]->xη.push_back(0.5 * (A[i][j]->neibor[n2]->x - A[i][j]->neibor[n4]->x));
						A[i][j]->yη.push_back(0.5 * (A[i][j]->neibor[n2]->y - A[i][j]->neibor[n4]->y));
						A[i][j]->J.push_back(1 / (A[i][j]->xξ[k] * A[i][j]->yη[k] - A[i][j]->xη[k] * A[i][j]->yξ[k]));
						A[i][j]->ξx.push_back(A[i][j]->yη[k] * A[i][j]->J[k]);
						A[i][j]->ξy.push_back(-A[i][j]->xη[k] * A[i][j]->J[k]);
						A[i][j]->ηx.push_back(-A[i][j]->yξ[k] * A[i][j]->J[k]);
						A[i][j]->ηy.push_back(A[i][j]->xξ[k] * A[i][j]->J[k]);
						A[i][j]->xτ.push_back(0);
						A[i][j]->yτ.push_back(0);
						A[i][j]->ξt.push_back(0);
						A[i][j]->ηt.push_back(0);
					}
					//A[i][j]->xξ.push_back(0.5 * (A[i][A[i][j]->neibor[0]].x - A[i][A[i][j]->neibor[2]].x));
					//A[i][j]->yξ.push_back(0.5 * (A[i][A[i][j]->neibor[0]].y - A[i][A[i][j]->neibor[2]].y));
					//A[i][j]->xη.push_back(0.5 * (A[i][A[i][j]->neibor[1]].x - A[i][A[i][j]->neibor[2]].x));
					//A[i][j]->yη.push_back(0.5 * (A[i][A[i][j]->neibor[1]].y - A[i][A[i][j]->neibor[2]].y));
					//A[i][j]->J.push_back(1 / (A[i][j]->xξ[0] * A[i][j]->yη[0] - A[i][j]->xη[0] * A[i][j]->yξ[0]));
					//A[i][j]->ξx.push_back(A[i][j]->yη[0] * A[i][j]->J[0]);
					//A[i][j]->ξy.push_back(-A[i][j]->xη[0] * A[i][j]->J[0]);
					//A[i][j]->ηx.push_back(-A[i][j]->yξ[0] * A[i][j]->J[0]);
					//A[i][j]->ηy.push_back(A[i][j]->xξ[0] * A[i][j]->J[0]);
					//A[i][j]->xτ.push_back(0);
					//A[i][j]->yτ.push_back(0);
					//A[i][j]->ξt.push_back(0);
					//A[i][j]->ηt.push_back(0);

					//A[i][j]->xξ.push_back(0.5 * (A[i][A[i][j]->neibor[1]].x - A[i][A[i][j]->neibor[0]].x));
					//A[i][j]->yξ.push_back(0.5 * (A[i][A[i][j]->neibor[1]].y - A[i][A[i][j]->neibor[0]].y));
					//A[i][j]->xη.push_back(0.5 * (A[i][A[i][j]->neibor[2]].x - A[i][A[i][j]->neibor[0]].x));
					//A[i][j]->yη.push_back(0.5 * (A[i][A[i][j]->neibor[2]].y - A[i][A[i][j]->neibor[0]].y));
					//A[i][j]->J.push_back(1 / (A[i][j]->xξ[1] * A[i][j]->yη[1] - A[i][j]->xη[1] * A[i][j]->yξ[1]));
					//A[i][j]->ξx.push_back(A[i][j]->yη[1] * A[i][j]->J[1]);
					//A[i][j]->ξy.push_back(-A[i][j]->xη[1] * A[i][j]->J[1]);
					//A[i][j]->ηx.push_back(-A[i][j]->yξ[1] * A[i][j]->J[1]);
					//A[i][j]->ηy.push_back(A[i][j]->xξ[1] * A[i][j]->J[1]);
					//A[i][j]->xτ.push_back(0);
					//A[i][j]->yτ.push_back(0);
					//A[i][j]->ξt.push_back(0);
					//A[i][j]->ηt.push_back(0);

					//A[i][j]->xξ.push_back(0.5 * (A[i][A[i][j]->neibor[2]].x - A[i][A[i][j]->neibor[1]].x));
					//A[i][j]->yξ.push_back(0.5 * (A[i][A[i][j]->neibor[2]].y - A[i][A[i][j]->neibor[1]].y));
					//A[i][j]->xη.push_back(0.5 * (A[i][A[i][j]->neibor[0]].x - A[i][A[i][j]->neibor[1]].x));
					//A[i][j]->yη.push_back(0.5 * (A[i][A[i][j]->neibor[0]].y - A[i][A[i][j]->neibor[1]].y));
					//A[i][j]->J.push_back(1 / (A[i][j]->xξ[2] * A[i][j]->yη[2] - A[i][j]->xη[2] * A[i][j]->yξ[2]));
					//A[i][j]->ξx.push_back(A[i][j]->yη[2] * A[i][j]->J[2]);
					//A[i][j]->ξy.push_back(-A[i][j]->xη[2] * A[i][j]->J[2]);
					//A[i][j]->ηx.push_back(-A[i][j]->yξ[2] * A[i][j]->J[2]);
					//A[i][j]->ηy.push_back(A[i][j]->xξ[2] * A[i][j]->J[2]);
					//A[i][j]->xτ.push_back(0);
					//A[i][j]->yτ.push_back(0);
					//A[i][j]->ξt.push_back(0);
					//A[i][j]->ηt.push_back(0);
				}
				else if (A[i][j]->neibor.size() == 4)
				{
					A[i][j]->xξ.push_back(0.5 * (A[i][j]->neibor[0]->x - A[i][j]->neibor[2]->x));
					A[i][j]->yξ.push_back(0.5 * (A[i][j]->neibor[0]->y - A[i][j]->neibor[2]->y));
					A[i][j]->xη.push_back(0.5 * (A[i][j]->neibor[1]->x - A[i][j]->neibor[3]->x));
					A[i][j]->yη.push_back(0.5 * (A[i][j]->neibor[1]->y - A[i][j]->neibor[3]->y));
					A[i][j]->J.push_back(1 / (A[i][j]->xξ[0] * A[i][j]->yη[0] - A[i][j]->xη[0] * A[i][j]->yξ[0]));
					A[i][j]->ξx.push_back(A[i][j]->yη[0] * A[i][j]->J[0]);
					A[i][j]->ξy.push_back(-A[i][j]->xη[0] * A[i][j]->J[0]);
					A[i][j]->ηx.push_back(-A[i][j]->yξ[0] * A[i][j]->J[0]);
					A[i][j]->ηy.push_back(A[i][j]->xξ[0] * A[i][j]->J[0]);
					A[i][j]->xτ.push_back(0);
					A[i][j]->yτ.push_back(0);
					A[i][j]->ξt.push_back(0);
					A[i][j]->ηt.push_back(0);

					A[i][j]->xξ.push_back(0.5 * (A[i][A[i][j]->neibor1[0]]->x - A[i][A[i][j]->neibor1[2]]->x));
					A[i][j]->yξ.push_back(0.5 * (A[i][A[i][j]->neibor1[0]]->y - A[i][A[i][j]->neibor1[2]]->y));
					A[i][j]->xη.push_back(0.5 * (A[i][A[i][j]->neibor1[1]]->x - A[i][A[i][j]->neibor1[3]]->x));
					A[i][j]->yη.push_back(0.5 * (A[i][A[i][j]->neibor1[1]]->y - A[i][A[i][j]->neibor1[3]]->y));
					A[i][j]->J.push_back(1 / (A[i][j]->xξ[1] * A[i][j]->yη[1] - A[i][j]->xη[1] * A[i][j]->yξ[1]));
					A[i][j]->ξx.push_back(A[i][j]->yη[1] * A[i][j]->J[1]);
					A[i][j]->ξy.push_back(-A[i][j]->xη[1] * A[i][j]->J[1]);
					A[i][j]->ηx.push_back(-A[i][j]->yξ[1] * A[i][j]->J[1]);
					A[i][j]->ηy.push_back(A[i][j]->xξ[1] * A[i][j]->J[1]);
					A[i][j]->xτ.push_back(0);
					A[i][j]->yτ.push_back(0);
					A[i][j]->ξt.push_back(0);
					A[i][j]->ηt.push_back(0);

				}
				else if (A[i][j]->neibor.size() > 4)
				{
					A[i][j]->xξ.push_back(0.5 * (A[i][A[i][j]->neibor1[0]]->x - A[i][A[i][j]->neibor1[2]]->x));
					A[i][j]->yξ.push_back(0.5 * (A[i][A[i][j]->neibor1[0]]->y - A[i][A[i][j]->neibor1[2]]->y));
					A[i][j]->xη.push_back(0.5 * (A[i][A[i][j]->neibor1[1]]->x - A[i][A[i][j]->neibor1[3]]->x));
					A[i][j]->yη.push_back(0.5 * (A[i][A[i][j]->neibor1[1]]->y - A[i][A[i][j]->neibor1[3]]->y));
					A[i][j]->J.push_back(1 / (A[i][j]->xξ[0] * A[i][j]->yη[0] - A[i][j]->xη[0] * A[i][j]->yξ[0]));
					A[i][j]->ξx.push_back(A[i][j]->yη[0] * A[i][j]->J[0]);
					A[i][j]->ξy.push_back(-A[i][j]->xη[0] * A[i][j]->J[0]);
					A[i][j]->ηx.push_back(-A[i][j]->yξ[0] * A[i][j]->J[0]);
					A[i][j]->ηy.push_back(A[i][j]->xξ[0] * A[i][j]->J[0]);
					A[i][j]->xτ.push_back(0);
					A[i][j]->yτ.push_back(0);
					A[i][j]->ξt.push_back(0);
					A[i][j]->ηt.push_back(0);

				}

			}
			else//t_sim不等于零时还需要坐标变换，应为动网格！
			{
				if (A[i][j]->neibor.size() == 4)
				{
					A[i][j]->xξ[0] = 0.5 * (A[i][j]->neibor[0]->x - A[i][j]->neibor[2]->x);
					A[i][j]->yξ[0] = 0.5 * (A[i][j]->neibor[0]->y - A[i][j]->neibor[2]->y);
					A[i][j]->xη[0] = 0.5 * (A[i][j]->neibor[1]->x - A[i][j]->neibor[3]->x);
					A[i][j]->yη[0] = 0.5 * (A[i][j]->neibor[1]->y - A[i][j]->neibor[3]->y);
					A[i][j]->J[0] = 1 / (A[i][j]->xξ[0] * A[i][j]->yη[0] - A[i][j]->xη[0] * A[i][j]->yξ[0]);
					A[i][j]->ξx[0] = A[i][j]->yη[0] * A[i][j]->J[0];
					A[i][j]->ξy[0] = -A[i][j]->xη[0] * A[i][j]->J[0];
					A[i][j]->ηx[0] = -A[i][j]->yξ[0] * A[i][j]->J[0];
					A[i][j]->ηy[0] = A[i][j]->xξ[0] * A[i][j]->J[0];
					A[i][j]->xτ[0] = (A[i][j]->x - Ar[A[i][j]->id].x) / dt;
					A[i][j]->yτ[0] = (A[i][j]->y - Ar[A[i][j]->id].y) / dt;
					A[i][j]->ξt[0] = -A[i][j]->xτ[0] * A[i][j]->ξx[0] - A[i][j]->yτ[0] * A[i][j]->ξy[0];
					A[i][j]->ηt[0] = -A[i][j]->xτ[0] * A[i][j]->ηx[0] - A[i][j]->yτ[0] * A[i][j]->ηy[0];
				}
				if (A[i][j]->neibor.size() > 4)
				{
					A[i][j]->xξ[0] = 0.5 * (A[i][A[i][j]->neibor1[0]]->x - A[i][A[i][j]->neibor1[2]]->x);
					A[i][j]->yξ[0] = 0.5 * (A[i][A[i][j]->neibor1[0]]->y - A[i][A[i][j]->neibor1[2]]->y);
					A[i][j]->xη[0] = 0.5 * (A[i][A[i][j]->neibor1[1]]->x - A[i][A[i][j]->neibor1[3]]->x);
					A[i][j]->yη[0] = 0.5 * (A[i][A[i][j]->neibor1[1]]->y - A[i][A[i][j]->neibor1[3]]->y);
					A[i][j]->J[0] = 1 / (A[i][j]->xξ[0] * A[i][j]->yη[0] - A[i][j]->xη[0] * A[i][j]->yξ[0]);
					A[i][j]->ξx[0] = A[i][j]->yη[0] * A[i][j]->J[0];
					A[i][j]->ξy[0] = -A[i][j]->xη[0] * A[i][j]->J[0];
					A[i][j]->ηx[0] = -A[i][j]->yξ[0] * A[i][j]->J[0];
					A[i][j]->ηy[0] = A[i][j]->xξ[0] * A[i][j]->J[0];
					A[i][j]->xτ[0] = (A[i][j]->x - Ar[A[i][j]->id].x) / dt;
					A[i][j]->yτ[0] = (A[i][j]->y - Ar[A[i][j]->id].y) / dt;
					A[i][j]->ξt[0] = -A[i][j]->xτ[0] * A[i][j]->ξx[0] - A[i][j]->yτ[0] * A[i][j]->ξy[0];
					A[i][j]->ηt[0] = -A[i][j]->xτ[0] * A[i][j]->ηx[0] - A[i][j]->yτ[0] * A[i][j]->ηy[0];

				}
				if (A[i][j]->neibor.size() == 3)
				{
					A[i][j]->xξ[0] = 0.5 * (A[i][j]->neibor[0]->x - A[i][j]->neibor[2]->x);
					A[i][j]->yξ[0] = 0.5 * (A[i][j]->neibor[0]->y - A[i][j]->neibor[2]->y);
					A[i][j]->xη[0] = 0.5 * (A[i][j]->neibor[1]->x - A[i][j]->neibor[2]->x);
					A[i][j]->yη[0] = 0.5 * (A[i][j]->neibor[1]->y - A[i][j]->neibor[2]->y);
					A[i][j]->J[0] = 1 / (A[i][j]->xξ[0] * A[i][j]->yη[0] - A[i][j]->xη[0] * A[i][j]->yξ[0]);
					A[i][j]->ξx[0] = A[i][j]->yη[0] * A[i][j]->J[0];
					A[i][j]->ξy[0] = -A[i][j]->xη[0] * A[i][j]->J[0];
					A[i][j]->ηx[0] = -A[i][j]->yξ[0] * A[i][j]->J[0];
					A[i][j]->ηy[0] = A[i][j]->xξ[0] * A[i][j]->J[0];
					A[i][j]->xτ[0] = (A[i][j]->x - Ar[A[i][j]->id].x) / dt;
					A[i][j]->yτ[0] = (A[i][j]->y - Ar[A[i][j]->id].y) / dt;
					A[i][j]->ξt[0] = -A[i][j]->xτ[0] * A[i][j]->ξx[0] - A[i][j]->yτ[0] * A[i][j]->ξy[0];
					A[i][j]->ηt[0] = -A[i][j]->xτ[0] * A[i][j]->ηx[0] - A[i][j]->yτ[0] * A[i][j]->ηy[0];

					A[i][j]->xξ[1] = 0.5 * (A[i][j]->neibor[1]->x - A[i][j]->neibor[0]->x);
					A[i][j]->yξ[1] = 0.5 * (A[i][j]->neibor[1]->y - A[i][j]->neibor[0]->y);
					A[i][j]->xη[1] = 0.5 * (A[i][j]->neibor[2]->x - A[i][j]->neibor[0]->x);
					A[i][j]->yη[1] = 0.5 * (A[i][j]->neibor[2]->y - A[i][j]->neibor[0]->y);
					A[i][j]->J[1] = 1 / (A[i][j]->xξ[1] * A[i][j]->yη[1] - A[i][j]->xη[1] * A[i][j]->yξ[1]);
					A[i][j]->ξx[1] = A[i][j]->yη[1] * A[i][j]->J[1];
					A[i][j]->ξy[1] = -A[i][j]->xη[1] * A[i][j]->J[1];
					A[i][j]->ηx[1] = -A[i][j]->yξ[1] * A[i][j]->J[1];
					A[i][j]->ηy[1] = A[i][j]->xξ[1] * A[i][j]->J[1];
					A[i][j]->xτ[1] = (A[i][j]->x - Ar[A[i][j]->id].x) / dt;
					A[i][j]->yτ[1] = (A[i][j]->y - Ar[A[i][j]->id].y) / dt;
					A[i][j]->ξt[1] = -A[i][j]->xτ[1] * A[i][j]->ξx[1] - A[i][j]->yτ[1] * A[i][j]->ξy[1];
					A[i][j]->ηt[1] = -A[i][j]->xτ[1] * A[i][j]->ηx[1] - A[i][j]->yτ[1] * A[i][j]->ηy[1];

					A[i][j]->xξ[2] = 0.5 * (A[i][j]->neibor[2]->x - A[i][j]->neibor[1]->x);
					A[i][j]->yξ[2] = 0.5 * (A[i][j]->neibor[2]->y - A[i][j]->neibor[1]->y);
					A[i][j]->xη[2] = 0.5 * (A[i][j]->neibor[0]->x - A[i][j]->neibor[1]->x);
					A[i][j]->yη[2] = 0.5 * (A[i][j]->neibor[0]->y - A[i][j]->neibor[1]->y);
					A[i][j]->J[2] = 1 / (A[i][j]->xξ[2] * A[i][j]->yη[2] - A[i][j]->xη[2] * A[i][j]->yξ[2]);
					A[i][j]->ξx[2] = A[i][j]->yη[2] * A[i][j]->J[2];
					A[i][j]->ξy[2] = -A[i][j]->xη[2] * A[i][j]->J[2];
					A[i][j]->ηx[2] = -A[i][j]->yξ[2] * A[i][j]->J[2];
					A[i][j]->ηy[2] = A[i][j]->xξ[2] * A[i][j]->J[2];
					A[i][j]->xτ[2] = (A[i][j]->x - Ar[A[i][j]->id].x) / dt;
					A[i][j]->yτ[2] = (A[i][j]->y - Ar[A[i][j]->id].y) / dt;
					A[i][j]->ξt[2] = -A[i][j]->xτ[2] * A[i][j]->ξx[2] - A[i][j]->yτ[2] * A[i][j]->ξy[2];
					A[i][j]->ηt[2] = -A[i][j]->xτ[2] * A[i][j]->ηx[2] - A[i][j]->yτ[2] * A[i][j]->ηy[2];

				}

			}
		}
	}
}
