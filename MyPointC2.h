#ifndef MY_POINTC2_H
#define MY_POINTC2_H


#include <CGAL/Origin.h>
#include <CGAL/Bbox_2.h>


class MyPointC2 {
  
  private:
  double vec[2];//位置ベクトル
  double vel[2];//速度ベクトル
  int col;//細胞タイプ
  double rad;//半径
  
  public:
  
  MyPointC2()
  : col(0), rad(0.25)
  {
    *vec = 0;
    *(vec+1) = 0;
    *vel = 0;
    *(vel+1) = 0;
  }
  
  
  MyPointC2(const double x, const double y, int c=0, double r=0.25, double vx=0.0,  double vy=0.0)
  : col(c), rad(r)
  {
    *vec = x;
    *(vec+1) = y;
    *vel=(vx);
    *(vel+1)=(vy);
    
  }
  
  const double& x() const  { return *vec; }
  const double& y() const { return *(vec+1); }
  double & x() { return *vec; }
  double& y() { return *(vec+1); }
  
  const double& vx() const  { return *vel; }
  const double& vy() const { return *(vel+1); }
  double & vx() { return *vel; }
  double& vy() { return *(vel+1); }
  
  int color() const { return col; }
  int& color() { return col; }
  
  double radius() const { return rad; }
  double& radius() { return rad; }
  
  
  
  
  //値として等しいか出はなく、オブジェクトして等しいかの判定
  bool operator==(const MyPointC2 &p) const
  {
    return (( *vec == *(p.vec) )  &&  
    *(vec+1) == *(p.vec + 1) 
    && (col == (p.col)) 
    && (rad==(p.rad)) 
    && (*vel == *(p.vel)) 
    && (*(vel+1) == *(p.vel+1)));
  }
  
  bool operator!=(const MyPointC2 &p) const
  {
    return !(*this == p);
  }
  
};

#endif // MY_POINTC2_H
