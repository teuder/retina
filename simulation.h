#ifndef SIMULATION_H
#define SIMULATION_H





#include <sys/stat.h>
#include <sstream>

using namespace std;

#include "main.h"
#include "cell.h"
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
  
  
  void Run();
  void MakeInitFile();
  
  

  Simulation():
    param(),
    csv_data(),
    voronoi(),
    delauny()
  {};
  
  private:
  Parameter param;

  stringstream csv_data;
  stringstream voronoi;
  stringstream delauny;
  
  int t;//時間
  DT dt;//Delaunay三角形分割
  
  //各点の位置
  vector_Point_2d points;
  vector_Point_2d points_L;
  vector_Point_2d points_R;
  vector_Point_2d points_U;
  vector_Point_2d points_D;
  vector_Point_2d points_LU;
  vector_Point_2d points_LD;
  vector_Point_2d points_RU;
  vector_Point_2d points_RD;
  vector_Point_2d points_all;//全てのオリジナル点と周期点を含む  
  
  //ドローネ三角形分割上での各点へのハンドル（ポインタ）
  vector_Vertex_handle_2d vtx;
  vector_Vertex_handle_2d vtx_R;
  vector_Vertex_handle_2d vtx_L;
  vector_Vertex_handle_2d vtx_U;
  vector_Vertex_handle_2d vtx_D;
  vector_Vertex_handle_2d vtx_RU;
  vector_Vertex_handle_2d vtx_RD;
  vector_Vertex_handle_2d vtx_LU;
  vector_Vertex_handle_2d vtx_LD;
    
  
  
  list<double> Speed_A1;
  list<double> Speed_A2;
  
  list<double> Vx_A1_inner;
  list<double> Vx_A1_bound;
  list<double> Vx_A1_all;
  
  list<double> Vx_A2_inner;  
  list<double> Vx_A2_bound;
  list<double> Vx_A2_all;
  
  ////ファイル出力のためのバッファ
  vector<string> str_t;//時間
  vector_string_2d str_points;//点の位置
  vector<string> str_delauny;//ドローネ辺
  vector<string> str_voronoi;//ボロノイ辺
  
  //vector<string> str_t2;//こちらのデータをサンプルした時間
  //vector<string> str_analyzed_data;//色々なデータ
  //vector<string> str_distance_inner;//距離と内積
  //vector<string> str_distance_inner_adjacent;//隣接距離と内積
  //vector<string> str_edges_inner_adjacent;//隣接エッジと内積
  //vector<string> str_diff_vx;//a8aの平均速度に対して、A9paの速度のx成分
  
  
  
  
  
  
  
  void Init();
  void InitForMakeInitFile();
  void InitPoints();
  void InitPoints(vector<Point>& points,const int& N,const double& L, const double& U, int t, double r);
  void InitPointsCircle(int i);
  void InitPointsForMakeInitFile();
  void DevidePointsCircle();
  void DevidePointsStripe();
  void InitPointsRegular();
  
  void InitPointsRetina();
  
  
  void InitPeriodicPoints();
  void InitTrianglation();
  
  void CalcVelocity();

  void UpdatePosition();
  void Apoptosis();
  void Apoptosis(int type, double Pd);
  
  void Ablation();
  
  void CacheVertices();
  void FoutVertices();
  void FoutCSVdata();
  void SaveVerticesToString();
  void SaveEdgesToString();
  
  
  void OutVertices();
  void SaveVx();
  
  
  
  void ShowPointsAll();
  void ShowPointsAllNum();
  void ShowVerticesAll();
  void ShowVerticesAllNum();
  
  
  void AnalyzeMovement();
  void AnalyzeMovement2();
  string OssMovement();
  void   OutMovement();
  
};


extern string OssFlags();
extern void OutFlags(string filepath=string("gflag.txt"));










#endif