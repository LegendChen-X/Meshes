#include "write_obj.h"
#include <fstream>
#include <cassert>
#include <iostream>

bool write_obj(
  const std::string & filename,
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::MatrixXd & UV,
  const Eigen::MatrixXi & UF,
  const Eigen::MatrixXd & NV,
  const Eigen::MatrixXi & NF)
{
    assert((F.size() == 0 || F.cols() == 3 || F.cols() == 4) && "F must have 3 or 4 columns");
    std::ofstream fp;
    fp.open(filename,std::ios::binary);
    
    if(!fp) return false;
    
    for(int i=0;i<V.rows();++i)
    {
        fp<<"v";
        for(int j=0;j<3;++j)
            fp<<" "<<V(i,j);
        fp<<"\n";
    }
    
    for(int i=0;i<UV.rows();++i)
    {
        fp<<"vt";
        for(int j=0;j<2;++j)
            fp<<" "<<UV(i,j);
        fp<<"\n";
    }
    
    for(int i=0;i<NV.rows();++i)
    {
        fp<<"vn";
        for(int j=0;j<3;++j)
            fp<<" "<<NV(i,j);
        fp<<"\n";
    }
    
    for(int i=0;i<F.rows();++i)
    {
        fp<<"f";
        for(int j=0;j<F.cols();++j)
            fp<<" "<<F(i,j)+1<<"/"<<UF(i,j)+1<<"/"<<NF(i,j)+1<<" ";
        fp<<"\n";
    }
    fp.close();
    
    return true;
}
