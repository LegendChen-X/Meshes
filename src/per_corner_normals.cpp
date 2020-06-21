#include "per_corner_normals.h"
#include "triangle_area_normal.h"
// Hint:
#include "vertex_triangle_adjacency.h"
#include <iostream>
#define DEBUG 0

void per_corner_normals(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const double corner_threshold,
  Eigen::MatrixXd & N)
{
    N = Eigen::MatrixXd::Zero(F.rows()*3,3);
    std::vector<std::vector<int>> VF;
    vertex_triangle_adjacency(F,V.rows(),VF);
    
#if(DEBUG)
    printf("vertex_triangle_adjacency check\n");
#endif
    
    int index = 0;
    for(int i=0;i<F.rows();++i)
    {
        for(int j=0;j<F.cols();++j)
        {
            Eigen::RowVector3d alpha(0,0,0);
            double beta = 0;
            Eigen::RowVector3d buff_f = triangle_area_normal(V.row(F(i,0)),V.row(F(i,1)),V.row(F(i,2)));
            for(auto l:VF[F(i,j)])
            {
#if(DEBUG)
                printf("VF[F(i,j)].size(): %lu\n",VF[F(i,j)].size());
#endif
                
#if(DEBUG)
                printf("triangle_area_normal start\n");
#endif
                
                Eigen::RowVector3d buff_g = triangle_area_normal(V.row(F(l,0)),V.row(F(l,1)),V.row(F(l,2)));
                
#if(DEBUG)
                printf("triangle_area_normal check\n");
#endif
                
#if(DEBUG)
                printf("value: %f\n",(buff_g.dot(buff_f) / (buff_g.norm()*buff_f.norm())));
#endif
                
                if((buff_g.dot(buff_f) / (buff_g.norm()*buff_f.norm())) > cos(corner_threshold*M_PI/180.0))
                {
                    alpha += buff_g;
                    beta += buff_g.norm();
                }
            }
            N.row(index) = (alpha/beta).normalized();
            index++;
        }
    }
}
