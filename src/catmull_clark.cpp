#include "catmull_clark.h"
#include <unordered_map>
#include <utility>
#include <functional>
#include <string>
#include <vector>
#include <iostream>

//https://en.wikipedia.org/wiki/Catmull%E2%80%93Clark_subdivision_surface

void remove_duplicate(std::vector<int> & v)
{
    std::vector<int>::iterator end = v.end();
    for(std::vector<int>::iterator i = v.begin();i!=end;++i)
        end = std::remove(i+1,end,*i);
    v.erase(end,v.end());
}

void catmull_clark(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const int num_iters,
  Eigen::MatrixXd & SV,
  Eigen::MatrixXi & SF)
{
    if(!num_iters) return;
    
    std::unordered_map<int, Eigen::RowVector3d> point_map;
    for(int i=0;i<F.rows();++i)
        point_map[i] = (V.row(F(i,0))+V.row(F(i,1))+V.row(F(i,2))+V.row(F(i,3)))/4.0;
    
    std::unordered_map<int,std::vector<int>> face_map;
    for(int i=0;i<F.rows();++i)
        for(int j=0;j<F.cols();++j)
            face_map[F(i,j)].push_back(i);
    
    std::unordered_map<int,std::vector<int>> vertice_map;
    for(int i=0;i<F.rows();++i)
    {
        for(int j=0;j<F.rows();++j)
        {
            vertice_map[F(i,j)].push_back(F(i,(j-1 + F.cols())%F.cols()));
            vertice_map[F(i,j)].push_back(F(i, (j+1)%F.cols()));
            remove_duplicate(vertice_map[F(i,j)]);
        }
    }
    
    std::unordered_map<std::string,std::vector<int>> edge_map;
    for(int i=0;i<F.rows();++i)
    {
        for(int j=0;j<F.cols();++j)
        {
            std::string key = std::to_string(F(i,j)) +  " " + std::to_string(F(i,(j+1)%F.cols()));
            std::string s_key = std::to_string(F(i,(j+1)%F.cols())) +  " " + std::to_string(F(i,j));
            edge_map[key].push_back(i);
            edge_map[s_key].push_back(i);
            
            remove_duplicate(edge_map[key]);
            remove_duplicate(edge_map[s_key]);
        }
    }
    
    SV = V;
    SF = F;
    
    for(int i=0;i<F.rows();++i)
    {
        for(int j=0;j<F.cols();++j)
        {
            std::vector<Eigen::RowVector3d> new_vertice;
            
            Eigen::RowVector3d P = V.row(F(i,j));
            
            Eigen::RowVector3d F_buff(0,0,0);
            for(int l=0;i<face_map[F(i,j)].size();++l)
                F_buff += point_map[face_map[F(i,j)][l]];
            F_buff = F_buff / face_map[F(i,j)].size();
            
            Eigen::RowVector3d R(0,0,0);
            for(int l=0;l<vertice_map[F(i,j)].size();++l)
               R += (P + V.row(vertice_map[F(i,j)][l])) / 2.0;
            R = R / vertice_map[F(i,j)].size();
            
            double n = face_map[F(i,j)].size();
            Eigen::RowVector3d barycenter = (F_buff + 2.0 * R + (n-3) * P) / n;
            
            new_vertice.push_back(barycenter);
            
            Eigen::RowVector3d new_point(0,0,0);
            std::string key = std::to_string(F(i,j)) +  " " + std::to_string(F(i,(j+1)%F.cols()));
            for(int l=0;i<edge_map[key].size();++l)
                new_point += point_map[edge_map[key][l]];
            new_point = ((new_point + V.row(F(i,j)) + V.row(F(i,(j+1)%F.cols()))) / 4.0);
            
            new_vertice.push_back(new_point);
            new_vertice.push_back(point_map[i]);
            
            Eigen::RowVector3d new_point_1(0,0,0);
            key = std::to_string(F(i,j)) +  " " + std::to_string(F(i,((j-1)+F.cols())%F.cols()));
            for(int l=0;i<edge_map[key].size();++l)
                new_point_1 += point_map[edge_map[key][l]];
            new_point_1 = (new_point_1 + V.row(F(i,j)) + V.row(F(i,((j-1)+F.cols())%F.cols()))) / 4.0;
            
            new_vertice.push_back(new_point_1);
            
            Eigen::RowVector4i newF(-1,-1,-1,-1);
            int index = 0;
            for(int i=0;i<new_vertice.size();++i)
            {
                for(int j=0;j<SV.rows();++j)
                    if((new_vertice.at(i)).isApprox(SV.row(j)))
                        newF(index) = j;
                if(newF[index]==-1)
                {
                    Eigen::MatrixXd newSV = Eigen::MatrixXd::Zero(SV.rows()+1,3);
                    newSV.topRows(SV.rows()) = SV;
                    newSV.bottomRows(1) = new_vertice.at(i);
                    SV = newSV;
                    newF(index) = SV.rows()-1;
                }
                index += 1;
            }
            
            Eigen::MatrixXi newSF = Eigen::MatrixXi::Zero(SF.rows()+1,4);
            newSF.topRows(SF.rows()) = SF;
            newSF.bottomRows(1) = newF;
            SF = newSF;
        }
    }
    catmull_clark(Eigen::MatrixXd(SV),Eigen::MatrixXi(SF),num_iters-1,SV,SF);
}
