#ifndef MY_POINTC2_H
#define MY_POINTC2_H


#include <CGAL/Origin.h>
#include <CGAL/Bbox_2.h>


class MyPointC2 {
  
  private:
  double vec[2];//位置ベクトル
  double vel[2];//速度ベクトル
  int col;//細胞タイプ
  int col1;//細胞タイプ
  double rad;//半径
  double pie;//角度
  
  
  public:
  
  MyPointC2()
  :col(0),col1(0), rad(0.25),pie(0.0)
  {
    
     *vec = 0;
     *(vec+1) = 0;
     *vel = 0;
     *(vel+1) = 0;
  }
  
  
  //MyPointC2(const double x, const double y, const int c=0, const double r=0.25, const double vx=0.0, const double vy=0.0)
  //{
    //col=c;
    //rad=r;
     //*vec = x;
     //*(vec+1) = y;
     //*vel=(vx);
     //*(vel+1)=(vy);
    
  //}
     MyPointC2(const double x, const double y, const int c=0, const double r=0.25, const double vx=0.0, const double vy=0.0,const int c1=0,  const double p=.0)
    //MyPointC2(const double x, const double y, const int c=0, const int c1=0,  const double p=.0, const double r=.25,  const double vx=.0, const  double vy=.0)
    {
      col=c;
      col1=c1;
      pie=p;
      rad=r;
      *vec = x;
      *(vec+1) = y;
      *vel=(vx);
      *(vel+1)=(vy);
      
    }
  
  const double& x() const  { return *vec; }
  double & x() { return *vec; }
  const double& y() const { return *(vec+1); }  
  double& y() { return *(vec+1); }
  
  const double& vx() const  { return *vel; }
  const double& vy() const { return *(vel+1); }
  double & vx() { return *vel; }
  double& vy() { return *(vel+1); }
  
  int color() const { return col; }
  int& color() { return col; }
  
  int color1() const { return col1; }
  int& color1() { return col1; }
  
  double radius() const { return rad; }
  double& radius() { return rad; }
  
  double pi() const { return pie; }
  double& pi() { return pie; }
  
  
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
