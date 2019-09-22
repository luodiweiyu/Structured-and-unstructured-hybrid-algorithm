#pragma once
void init_mesh();
void init_mesh1();
void getType();
void initFlow();
void init_flow_uniform();
void init_flow_shockwave();
void init_flow_shockwaveCross();
void init_polygon_mesh();
void init_flow_cylinder();//圆柱绕流

//void init_Ar();
void findConnectPoint();
void reorderMesh();
void init_flow_normal();
void remesh_bound();//将边界x或y坐标移动到和内部点相同，改善边界条件

//void init_U();
