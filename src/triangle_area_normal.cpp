#include "triangle_area_normal.h"
#include <Eigen/Geometry>

Eigen::RowVector3d triangle_area_normal(
  const Eigen::RowVector3d & a, 
  const Eigen::RowVector3d & b, 
  const Eigen::RowVector3d & c)
{
    Eigen::RowVector3d norm = (0.5 * (b-a).cross(c-a).norm()) * (b-a).cross(c-a).normalized();
    return norm;
}
