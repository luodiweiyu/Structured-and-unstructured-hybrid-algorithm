# Structured and unstructured hybrid algorithm
 design new CFD compute method for HPC
 first version 2019.09.19
 achive simple hybrid grid compute
 2019.10.08
 此版本实现了钝头体绕流算例，但是由于物体表面的点是由非结构点连接法向而形成，物体表面的网格在实时移动，不能实现动网格。
 将采用密集散点布在物体表面，寻找最近点的算法进行尝试