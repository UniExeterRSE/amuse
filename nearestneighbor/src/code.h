#ifndef CODE_H
#define CODE_H

double distance_between_points(double x0, double y0, double z0, double x1, double y1, double z1);

class DistanceAndIndex{
    
public:
    double r;
    int index;
    
    DistanceAndIndex();
    DistanceAndIndex(DistanceAndIndex & original);
};

int find_nearest_neighbors(int npoints, double * x, double * y, double * z, int * n0, int * n1, int * n2);

#endif
