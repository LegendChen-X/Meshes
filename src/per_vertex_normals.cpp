#include "per_vertex_normals.h"
#include "triangle_area_normal.h"

void per_vertex_normals(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & N)
{
    N = Eigen::MatrixXd::Zero(V.rows(),3);
    double y = 0;
    for(int i=0;i<V.rows;++i)
    {
        for(int j=0;j<F.rows();++j)
        {
            Eigen::RowVector3d x(0,0,0);
            for(int l=0;l<F.cols();++l)
            {
                if(F[j][l] == i)
                {
                    Eigen::RowVector3d buff = triangle_area_normal(V.row(F(j,0)),V.row(F(j,1)),V.row(F(j,2)));
                    x += buff;
                    y += (double)buff.norm();
                    N.row(i) = (x/y).normalized();
                }
            }
            y = 0;
        }
    }
}
