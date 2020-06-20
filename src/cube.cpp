#include "cube.h"

void cube(
  Eigen::MatrixXd & V,
  Eigen::MatrixXi & F,
  Eigen::MatrixXd & UV,
  Eigen::MatrixXi & UF,
  Eigen::MatrixXd & NV,
  Eigen::MatrixXi & NF)
{
    V.resize(8,3);
    F.resize(6,4);
    UV.resize(14,2);
    UF.resize(6,4);
    NV.resize(6,3);
    NF.resize(6,4);
    
    V<<0,0,0;
    V<<1,0,0;
    V<<1,0,1;
    V<<0,0,1;
    V<<1,1,0;
    V<<1,1,1;
    V<<0,1,0;
    V<<0,1,1;
    
    F<<0,3,2,1;
    F<<3,7,5,2;
    F<<7,6,4,5;
    F<<6,0,1,4;
    F<<3,0,6,7;
    F<<5,4,1,2;
    
    UV<<1,4;//0 0
    UV<<1,3;//3 1
    UV<<2,3;//2 2
    UV<<2,4;//1 3
    UV<<1,2;//7 4
    UV<<2,2;//5 5
    UV<<1,1;//6 6
    UV<<2,1;//4 7
    UV<<1,0;//0 8
    UV<<2,0;//1 9
    UV<<0,2;//3 10
    UV<<0,1;//0 11
    UV<<3,1;//1 12
    UV<<3,2;//2 13
    
    UV /= 4.0;
    
    UF<<0,1,2,3;
    UF<<1,4,5,2;
    UF<<4,6,7,5;
    UF<<6,8,9,7;
    UF<<10,11,6,4;
    UF<<5,7,9,13;
    
    NV<<0,-1,0;
    NV<<0,0,1;
    NV<<0,1,0;
    NV<<0,0,-1;
    NV<<-1,0,0;
    NV<<1,0,0;
    
    NF<<0,0,0,0;
    NF<<1,1,1,1;
    NF<<2,2,2,2;
    NF<<3,3,3,3;
    NF<<4,4,4,4;
    NF<<5,5,5,5;
}
