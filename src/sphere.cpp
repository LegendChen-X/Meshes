#include "sphere.h"
#include <iostream>

void sphere(
  const int num_faces_u,
  const int num_faces_v,
  Eigen::MatrixXd & V,
  Eigen::MatrixXi & F,
  Eigen::MatrixXd & UV,
  Eigen::MatrixXi & UF,
  Eigen::MatrixXd & NV,
  Eigen::MatrixXi & NF)
{
    int total_f = num_faces_u * num_faces_v;
    int total_v = (num_faces_u + 1) * (num_faces_v + 1);
    V.resize(total_v,3);
    F.resize(total_f,4);
    UV.resize(total_v,2);
    UF.resize(total_f,4);
    NV.resize(total_v,3);
    Nf.resize(total_f,4);
    
    double alpha = 2.0 * M_PI / num_faces_u;
    double beta = M_PI / num_faces_v;
    
    int index = 0;
    
    for(int i=0;i<num_faces_u+1;++i)
    {
        for(int j=0;j<num_faces_v+1;++j)
        {
            double x = cos(alpha*i) * sin(beta*j);
            double y = sin(alpha*i) * sin(beta*j);
            double z = cos(beta*j);
            
            V.row(index) = Eigen::Vector3d(x,y,z);
            
            double p_x = (double) i / (double) (num_faces_u+1);
            double p_y = (double) j / (double) (num_faces_v+1);
            UV.row(index) = Eigen::Vector2d(p_x,p_y);
            
            NV.row(index) = Eigen::Vector3d(x,y,z);
            
            index += 1;
        }
    }
    
    index = 0;
    for(int i=0;i<num_faces_u;++i)
    {
        for(int j=0;j<num_faces_v;++j)
        {
            int v_1 = i*(num_faces_v+1)+j;
            int v_2 = v_1+1;
            int v_3 = (i+1)*(num_faces_v+1)+j;
            int v_4 = v_3+1;

            F.row(index) = Eigen::RowVector4i(v_1,v_2,v_3,v_4);
            UF.row(index) = Eigen::RowVector4i(v_1,v_2,v_3,v_4);
            NF.row(index) = Eigen::RowVector4i(v_1,v_2,v_3,v_4);
            
            index += 1;
        }
    }
}
