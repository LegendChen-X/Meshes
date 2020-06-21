#include "per_corner_normals.h"
#include "triangle_area_normal.h"
// Hint:
#include "vertex_triangle_adjacency.h"
#include <iostream>

void per_corner_normals(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const double corner_threshold,
  Eigen::MatrixXd & N)
{
    N = Eigen::MatrixXd::Zero(F.rows()*3,3);
    std::vector<std::vector<int>> VF;
    vertex_triangle_adjacency(F,V.rows(),VF);
    double x = cos(corner_threshold*M_PI/180.0);
    
    int index = 0;
    for(int i=0;i<F.rows();++i)
    {
        for(int j=0;j<F.cols();++j)
        {
            Eigen::RowVector3d alpha(0,0,0);
            double beta = 0;
            Eigen::RowVector3d buff_f = triangle_area_normal(V.row(F(i,0)),V.row(F(i,1)),V.row(F(i,2)));
            for(int l=0;i<VF[F(i,j)].size();++l)
            {
                Eigen::RowVector3d buff_g = triangle_area_normal(V.row(F(VF[F(i,j)][l],0)),V.row(F(VF[F(i,j)][l],1)),V.row(F(VF[F(i,j)][l],2)));
                if((buff_g.dot(buff_f) / (buff_g.norm()*buff_f.norm())) > x)
                {
                    alpha += buff_g;
                    beta += buff_g.norm();
                }
            }
            N.row(index) = (alpha/beta).normalized();
        }
    }
}
