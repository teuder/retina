#ifndef SIMULATION_H
#define SIMULATION_H





#include <sys/stat.h>
#include <sstream>

using namespace std;

#include "main.h"
#include "utility.h"
#include "parameter.h"


#include <CGAL/Object.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/squared_distance_2.h>

#include "MyKernel.h"
#include "MyPointC2_iostream.h"
typedef MyKernel<double>                   MK;
typedef CGAL::Filtered_kernel_adaptor<MK>  K;
typedef CGAL::Delaunay_triangulation_2<K>  DT;

typedef K::Point_2         Point;
typedef K::Segment_2       Segment;
typedef K::Ray_2           Ray;
typedef K::Line_2          Line;

typedef DT::Vertex_iterator Vertex_iterator;
typedef DT::Vertex_handle   Vertex_handle;


typedef vector<Point>           vector_Point;
typedef vector<vector<Point>>   vector_Point_2d;
typedef vector<Vertex_handle>   vector_Vertex_handle;
typedef vector<vector_Vertex_handle>   vector_Vertex_handle_2d;
typedef vector<string>           vector_string;
typedef vector<vector_string>           vector_string_2d;






class Simulation 
{
  
  public:
  
  
  void Run(int& argc,char** &argv);
  void RunMakeInitFile();
  void RunGrowingRetina();
  

  Simulation():
    param(),
    csv_data(),
    voronoi(),
    delauny()
  {  };
  
  private:
  Parameter param;

  
  int t;//時間
  DT dt;//Delaunay三角形分割
   
  //ドローネ三角形分割上での各点へのハンドル
  vector_Vertex_handle_2d vtx;
  vector_Vertex_handle_2d vtx_R;
  vector_Vertex_handle_2d vtx_L;
  vector_Vertex_handle_2d vtx_U;
  vector_Vertex_handle_2d vtx_D;
  vector_Vertex_handle_2d vtx_RU;
  vector_Vertex_handle_2d vtx_RD;
  vector_Vertex_handle_2d vtx_LU;
  vector_Vertex_handle_2d vtx_LD;
     
  ////ファイル出力のためのバッファ
  vector<string> str_t;//時間
  vector_string_2d str_points;//点の位置
  vector<string> str_delauny;//ドローネ辺
  vector<string> str_voronoi;//ボロノイ辺
  stringstream csv_data;
  stringstream voronoi;
  stringstream delauny;
  
  
  void Init();
  void InitFirstCell();
  
  void InitPoints();
  void InitPoints(vector_Vertex_handle& vtx,const int& N,const double& L, const double& U, int t, double r);
  void InitPointsCircle(int i);
  void InitPointsForMakeInitFile();
  void InitPointsRegular();
  void InitPointsRetina();
  void InitPointsCheckered();
  
  void DevidePointsCircle();
  void DevidePointsStripe();
  void DevidePointsRandom();
  void DevidePointsRandomRegular();
  
  void InitPeriodicPoints();
  
  void CalcVelocity();
  void CalcVelocityByMembrane();
  void UpdatePosition();
  
  void Apoptosis();
  void Apoptosis(int type, double Pd);
  void Ablation();
  
  void CacheVertices();
  void CacheEdges();
  void FoutData();
  void FoutEdges();
  
  void ShowVerticesAll();
  void ShowVerticesAllNum();
  
  
};


extern string OssFlags();
extern void OutFlags(string filepath=string("gflag.txt"));










#endif