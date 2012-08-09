#include "simulation.h"
#include "utility.h"
#include "random.h"

#include <GL/glut.h>
#include "glut.hpp"



void display(void)
{
  glClear(GL_COLOR_BUFFER_BIT);
  DrawPrimitive(GL_TRIANGLE_FAN) {
    glColor3d(1.0, 0.0, 0.0); /* 赤 */
    //glColor3d(0.0, 1.0, 0.0); /* 緑 */
    //glColor3d(0.0, 0.0, 1.0); /* 青 */
    //glColor3d(1.0, 1.0, 0.0); /* 黄 */    
    glVertex2d(0.0, 0.0);
    glVertex2d(-0.5, -0.5);
    glVertex2d(0.9, -0.9);
    glVertex2d(0.3, 0.9);
    glVertex2d(-0.9, 0.9);
  }
  
  glFlush();
}






void Simulation::RunGrowingRetina(){
  BGN;
  InitFirstCell();
  
  if(FLAGS_v) ShowVerticesAllNum();
  int tmax = static_cast<int>(FLAGS_Tm/FLAGS_H);
  bool ADD=true;
  
  //時間発展
  for(t=0;t<=tmax;++t){
    
    const double T=t*FLAGS_H;
    if(FLAGS_v){
      if(t==tmax) {
        cout << "\r";
        cout << "T =       " << flush;
        cout << "\r";
        cout << "T = " << T << endl;
      }
      else{      
        cout << "\r";
        cout << "T = " << T << flush;
      }
    }
    
    
    
    //時間を進めつつ、
    //ある時間間隔でCellを追加してゆく
    if(ADD && t%int(FLAGS_Fs/FLAGS_H)==0){
      //
      int success=0;
      int fail=0;
      while(1){
        
        int i=Prob(param.P);
        
        double Rmax=FLAGS_gr*T+FLAGS_gb;
        double Rmin=FLAGS_gr*T;
        
        const double r=NMrand(Rmin,Rmax);
        const double pi=2.0*myrand::PI*Zp1rand();
        
        const double X=r*cos(pi);
        const double Y=r*sin(pi);
        
        
        
        const MyPointC2 newpoint(X,Y,i,param.r[i],0.0,0.0);
        const MyPointC2 nearest(dt.nearest_vertex(newpoint)->point());
        const double sqr_dist=
        (nearest.x()-newpoint.x())*(nearest.x()-newpoint.x())+
        (nearest.y()-newpoint.y())*(nearest.y()-newpoint.y());
        
        
        if(sqr_dist>=FLAGS_inispace*(nearest.radius()+newpoint.radius())*(nearest.radius()+newpoint.radius())){
        
          ++success;
          vtx[i].push_back(dt.insert(newpoint));
          
          if(FLAGS_periodic>0){
           const  MyPointC2 newpoint_R(newpoint.x()+FLAGS_Xp,newpoint.y(),newpoint.color(),newpoint.radius(),0.0,0.0);
           const  MyPointC2 newpoint_L(newpoint.x()-FLAGS_Xp,newpoint.y(),newpoint.color(),newpoint.radius(),0.0,0.0);

            vtx_R[i].push_back(dt.insert(newpoint_R));
            vtx_L[i].push_back(dt.insert(newpoint_L));
            
            if(FLAGS_periodic>1){
              const MyPointC2 newpoint_U(newpoint.x(),newpoint.y()+FLAGS_Yp,newpoint.color(),newpoint.radius(),0.0,0.0);
              const MyPointC2 newpoint_D(newpoint.x(),newpoint.y()-FLAGS_Yp,newpoint.color(),newpoint.radius(),0.0,0.0);
             const  MyPointC2 newpoint_RU(newpoint.x()+FLAGS_Xp,newpoint.y()+FLAGS_Yp,newpoint.color(),newpoint.radius(),0.0,0.0);
             const  MyPointC2 newpoint_RD(newpoint.x()+FLAGS_Xp,newpoint.y()-FLAGS_Yp,newpoint.color(),newpoint.radius(),0.0,0.0);
             const  MyPointC2 newpoint_LU(newpoint.x()-FLAGS_Xp,newpoint.y()+FLAGS_Yp,newpoint.color(),newpoint.radius(),0.0,0.0);
              const MyPointC2 newpoint_LD(newpoint.x()-FLAGS_Xp,newpoint.y()-FLAGS_Yp,newpoint.color(),newpoint.radius(),0.0,0.0);
              
              vtx_U[i].push_back(dt.insert(newpoint_U));
              vtx_D[i].push_back(dt.insert(newpoint_D));
              vtx_RU[i].push_back(dt.insert(newpoint_RU));
              vtx_RD[i].push_back(dt.insert(newpoint_RD));
              vtx_LU[i].push_back(dt.insert(newpoint_LU));
              vtx_LD[i].push_back(dt.insert(newpoint_LD));
            }
          }
          
        }
        else{
          ++fail; 
        }
        if(success > 100) {break;}
        if(fail > 200) {ADD=true;break;} 
      }
    }
    
    CalcVelocity();
    
    //細胞の位置や速度を記録
    if((T >= FLAGS_Ts) && (t%int(FLAGS_Fs/FLAGS_H)==0 || t==tmax)){
      CacheVertices();
    }
    
    UpdatePosition();
  }
  
  //DevidePointsRandom();
  //DevidePointsRandomRegular();
  if(FLAGS_v) ShowVerticesAllNum();
  CacheVertices();
  FoutData();
  
  END;
}

void Simulation::Run(int& argc,char** &argv){
  BGN;
  
  Init();//普通の初期化
  
  if(FLAGS_display){
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGBA);
    glutCreateWindow(argv[0]);
    glClearColor(1.0, 1.0, 1.0, 1.0);
    glutDisplayFunc(display);
    glutMainLoop();
  }
  
  
  if(FLAGS_v){
    ShowVerticesAllNum();
  }
  
  //#########################################################
  //点の位置を更新する
  int tmax = static_cast<int>(FLAGS_Tm/FLAGS_H);
  for(t=0;t<=tmax;++t){
    
    const double T=t*FLAGS_H;
    
    if(FLAGS_v){
      if(t==tmax) {
        cout << "\r";
        cout << "T =       " << flush;
        cout << "\r";
        cout << "T = " << T << endl;
      }
      else{      
        cout << "\r";
        cout << "T = " << T << flush;
      }
    }
  
    CalcVelocity();
    
    //細胞の位置や速度を記録
    if((T >= FLAGS_Ts) && (t%int(FLAGS_Fs/FLAGS_H)==0 || t==tmax)){
      if(FLAGS_anime) CacheVertices();
      if(FLAGS_edges) CacheEdges();
    }
    
    UpdatePosition();
    

  }
  
  
  FoutData();
  if(FLAGS_edges) FoutEdges();
  END;
}





void Simulation::RunMakeInitFile(){
  BGN;
  
  
  //FLAGS_H=0.001;
  //FLAGS_Tm=10;
  //#########################################################
  //点の位置を更新する
  //const double Tm=FLAGS_Tm;w
  //const double H=FLAGS_H;
  
  //最初の細胞
  vtx.resize(1);
  vtx[0].push_back(dt.insert_first(MyPointC2(0,0,0,FLAGS_r,0.0,0.0)));
  
  if(FLAGS_v) ShowVerticesAllNum();
  int tmax = static_cast<int>(FLAGS_Tm/FLAGS_H);
  bool ADD=true;
  
  //時間発展
  for(t=0;t<=tmax;++t){
    
    const double T=t*FLAGS_H;
    if(FLAGS_v){
      if(t==tmax) {
        cout << "\r";
        cout << "T =       " << flush;
        cout << "\r";
        cout << "T = " << T << endl;
      }
      else{      
        cout << "\r";
        cout << "T = " << T << flush;
      }
    }
    
    
    
    //時間を進めつつ、
    //ある時間間隔でCellを追加してゆく
    if(ADD && t%int(FLAGS_Fs/FLAGS_H)==0){
      //
      int success=0;
      int fail=0;
      while(1){
        const double r=FLAGS_Xp*zP1rand();
        const double pi=2.0*myrand::PI*Zp1rand();
        
        const double X=r*cos(pi);
        const double Y=r*sin(pi);
        
        
        //const double X=dNMrand(-0.5*FLAGS_Xp,0.5*FLAGS_Xp);
        //const double Y=dNMrand(FLAGS_L,FLAGS_U);
        MyPointC2 newpoint(X,Y,0,FLAGS_r,0.0,0.0);
        MyPointC2 nearest(dt.nearest_vertex(newpoint)->point());
        
        if(CGAL::squared_distance(Point(nearest),Point(newpoint))>(nearest.radius()+newpoint.radius())*(nearest.radius()+newpoint.radius())){
          ++success;
          dt.insert(newpoint);
          vtx[0].push_back(dt.nearest_vertex(newpoint));
          
          if(FLAGS_periodic>0){
            MyPointC2 newpoint_R(newpoint.x()+FLAGS_Xp,newpoint.y(),newpoint.color(),newpoint.radius(),0.0,0.0);
            MyPointC2 newpoint_L(newpoint.x()-FLAGS_Xp,newpoint.y(),newpoint.color(),newpoint.radius(),0.0,0.0);
            dt.insert(newpoint_R);
            dt.insert(newpoint_L);
            vtx_R[0].push_back(dt.nearest_vertex(newpoint_R));
            vtx_L[0].push_back(dt.nearest_vertex(newpoint_L));
            
            if(FLAGS_periodic>1){
              MyPointC2 newpoint_U(newpoint.x(),newpoint.y()+FLAGS_Yp,newpoint.color(),newpoint.radius(),0.0,0.0);
              MyPointC2 newpoint_D(newpoint.x(),newpoint.y()-FLAGS_Yp,newpoint.color(),newpoint.radius(),0.0,0.0);
              MyPointC2 newpoint_RU(newpoint.x()+FLAGS_Xp,newpoint.y()+FLAGS_Yp,newpoint.color(),newpoint.radius(),0.0,0.0);
              MyPointC2 newpoint_RD(newpoint.x()+FLAGS_Xp,newpoint.y()-FLAGS_Yp,newpoint.color(),newpoint.radius(),0.0,0.0);
              MyPointC2 newpoint_LU(newpoint.x()-FLAGS_Xp,newpoint.y()+FLAGS_Yp,newpoint.color(),newpoint.radius(),0.0,0.0);
              MyPointC2 newpoint_LD(newpoint.x()-FLAGS_Xp,newpoint.y()-FLAGS_Yp,newpoint.color(),newpoint.radius(),0.0,0.0);
              
              dt.insert(newpoint_U);
              dt.insert(newpoint_D);
              dt.insert(newpoint_RU);
              dt.insert(newpoint_RD);
              dt.insert(newpoint_LU);
              dt.insert(newpoint_LD);
              
              vtx_U[0].push_back(dt.nearest_vertex(newpoint_U));
              vtx_D[0].push_back(dt.nearest_vertex(newpoint_D));
              vtx_RU[0].push_back(dt.nearest_vertex(newpoint_RU));
              vtx_RD[0].push_back(dt.nearest_vertex(newpoint_RD));
              vtx_LU[0].push_back(dt.nearest_vertex(newpoint_LU));
              vtx_LD[0].push_back(dt.nearest_vertex(newpoint_LD));
            }
          }
          
        }
        else{
          ++fail; 
        }
        if(success > 100) {break;}
        if(fail > 200) {ADD=true;break;} 
      }
    }
    
    CalcVelocity();
    
    //細胞の位置や速度を記録
    if((T >= FLAGS_Ts) && (t%int(FLAGS_Fs/FLAGS_H)==0 || t==tmax)){
      CacheVertices();
      //SaveVx();
    }
    
    UpdatePosition();
  }
  
  //DevidePointsRandom();
  DevidePointsRandomRegular();
  if(FLAGS_v) ShowVerticesAllNum();
  CacheVertices();
  FoutData();
  END;
}


























void Simulation::CacheVertices(){
  BGN;
  
  for(auto vv : vtx){
    for(auto v : vv){
      const auto &p=v->point();
      //cout << v->point() << endl;
      //fout << "t,col,x,y,vx,vy,r";
      csv_data 
      << t*FLAGS_H << "," 
      << p.color() << ","
      << p.x() << ","
      << p.y() << ","
      << p.vx() << ","
      << p.vy() << ","
      << p.radius() << endl;
      
      //DT::Vertex_circulator vi,v0;
      //vi=v0=dt.incident_vertices(v);//vと隣接する頂点をめぐるサーキュレータ
      //int counter=0;
      //do{
      //if(!dt.is_infinite(vi)){//viが無限遠にある頂点でないなら
      
      //const auto &pi=vi->point();
      //const auto X=(pi.x()-p.x());
      //const auto Y=(pi.y()-p.y());
      
      //csv_data 
      //<< sqrt(X*X+Y*Y) << ","
      //<< CalcRadian(X,Y) << ","
      //<< pi.color() << ",";
      //++counter;
      //}
      //if(counter==8) break;//隣接細胞が多すぎるのは例外的なので打ち切り
      //}while(++vi!=v0);
      
      //int n=counter;
      //for(;counter<8;++counter){
      //csv_data << ",,,";
      //}
      //csv_data << n << ",";
      //csv_data
      //<< FLAGS_Pd1 << ","
      //<< FLAGS_Pd2 << ","
      //<< FLAGS_x1 << ","
      //<< FLAGS_y1 << ","
      //<< FLAGS_d1 << ","
      //<< FLAGS_repli
      //<< endl;
      
    }
  }
  
  //シミュレーションパラメターの出力
  
  
  
  
  
  
  END;
}



void Simulation::CalcVelocity()
{
  
  BGN;
  if(dt.number_of_vertices()>1){
  
  for(auto vv : vtx){    
    for(auto v : vv){
      
      double vx=0.0;
      double vy=0.0;
      const double x=(v)->point().x();
      const double y=(v)->point().y();
      const int c=(v)->point().color();
      
      if(!param.fix[c]){
        
        DT::Vertex_circulator vi,v0;
        
        vi=v0=dt.incident_vertices(v);//vと隣接する頂点をめぐるサーキュレータ
        
        do{
          
          if(!dt.is_infinite(vi)){//viが無限遠にある頂点でないなら
            
            const double xi=(vi->point().x());
            const double yi=(vi->point().y());
            const int    ci=(vi->point().color());
            const double X=(xi-x);
            const double Y=(yi-y);
            const double r2=X*X+Y*Y;
            const double r3=r2*sqrt(r2);
            const double r5=r2*r3;
            const double r_opt=((v)->point().radius()+vi->point().radius());//最適距離
            const double sigma=r_opt*(1.0/sqrt(3.0));
            const double sigma3=3.0*sigma*sigma*sigma;
            const double ADHESION=param.A[c][ci]*((sigma/r3)-(sigma3/r5));
            const double ROTATION=param.R[c][ci]*(r_opt/r2);
            
            vx+=ADHESION*X;
            vy+=ADHESION*Y;
            if(t*FLAGS_H>FLAGS_Ta){
              vx+=ROTATION*Y;
              vy+=ROTATION*(-X);
            }
            
          }
        }while(++vi!=v0);
      }
      //cout << "x=" << x << " y=" << y << "vx=" << vx << " vy=" << vy << endl;
      
      //(v)->point().vx()=vx;
      //(v)->point().vy()=vy; 
        
      if(FLAGS_noise>0.0 &&  t*FLAGS_H>FLAGS_Tnoise){
        
        (v)->point().vx()=vx+NormRand(0,FLAGS_noise);
        (v)->point().vy()=vy+NormRand(0,FLAGS_noise);
      }
      else{
        
        (v)->point().vx()=vx;
        (v)->point().vy()=vy;
      }
    }
    }
  }
  END;
}


//細胞を決まった比率で分割、順番はランダム
void Simulation::DevidePointsRandomRegular(){
  BGN;
  
  //cout << endl;
  //ShowVerticesAllNum();
  vector_Vertex_handle_2d vtx_orig;
  swap(vtx,vtx_orig);
  
  const size_t n_cell_type=4;
  vtx.resize(n_cell_type);
  
  int Ntotal=vtx_orig[0].size();
  //cout << "Ntotal " << Ntotal << endl;
  
  vector<int> N(n_cell_type);//各細胞タイプをいくつ作るか
  int sum=0;
  for(auto i : seq(n_cell_type)){
    int hoge=floor(Ntotal*param.P[i]);
    N[i]=hoge;
    sum+=hoge;
    //cout << "hoge " << hoge << endl; 
    //sum+=N[i]=floor(Ntotal*param.P[i]);
  }
  //cout << "sum " << sum << endl; 
  //cout << "Ntotal-sum " << Ntotal-sum << endl;
  int Nres=(Ntotal-sum);
  while(Nres-- != 0){
    N[Prob(param.P)]+=1;
  }
  
  
  vector<int> colors;
  colors.reserve(Ntotal);
  
  
  for(auto i : seq(n_cell_type)){
    while(N[i]-- !=0) colors.push_back(i);
  }
  
  random_shuffle(colors.begin(),colors.end(),myrandom);
  
  auto c=colors.begin();
  auto v=vtx_orig[0].begin();
  while(c!=colors.end()){
    vtx[*c].push_back(*v);
    (*v)->point().color()=*c;
    ++c;++v;
  }
    
  END;
}


void Simulation::DevidePointsRandom(){
  BGN;
  
  
  cout << endl;
  ShowVerticesAllNum();
  vector_Vertex_handle_2d vtx_orig;
  swap(vtx,vtx_orig);
  
  const size_t n_cell_type=4;
  vtx.resize(n_cell_type);
  
  
  cout << endl;
  cout << "BEFORE" << endl;
  ShowVerticesAllNum();
  for(auto v:vtx_orig[0]){
    int i=Prob(param.P);
    v->point().color()=i;
    vtx[i].push_back(v);
    
    
    
  }
  
  cout << "AFTER" << endl;
  ShowVerticesAllNum();
  
  END;
}







void Simulation::FoutData(){
  BGN;
  
  //if(FLAGS_makeinit) mkdir(FLAGS_I.c_str(), 0777);
  //else mkdir(FLAGS_o.c_str(), 0777);
  if(FLAGS_makeinit) mkdir(FLAGS_I.c_str(), 0777);
  mkdir(FLAGS_o.c_str(), 0777);
  
  
  string filename=FLAGS_o+string("/output.csv");
  ofstream fout(filename.c_str(),ios::trunc);
  fout << "t,col,x,y,vx,vy,r";
  fout << endl;
  fout << csv_data.str();
  fout.close();
  
  
  size_t n=vtx.size();
  for(size_t i=0;i<n;++i){
    ostringstream oss; 
    if(FLAGS_makeinit) oss << FLAGS_I << "/I" << i;
    else               oss << FLAGS_o << "/I" << i;
    
    fout.open(oss.str().c_str());
    for(auto v:vtx[i]) fout << MyPointC2(v->point()) << endl;
    fout.close();
  }
  
  //fout << "t,col,x,y,vx,vy,r,r1,th1,c1,r2,th2,c2,r3,th3,c3,r4,th4,c4,r5,th5,c5,r6,th6,c6,r7,th7,c7,r8,th8,c8,adj";
  END;
}






bool out_range(Point& p){
  return (
    p.x() >  0.5*FLAGS_Xp ||
    p.x() <= (-0.5*FLAGS_Xp + 0.4*FLAGS_r) ||
    p.y() < 0 ||
    p.y() >= (FLAGS_Yp - 0.4*FLAGS_r)
    );
};











void Simulation::Init(){
  BGN;
  
  vtx.resize(FLAGS_celltype);
  vtx_L.resize(FLAGS_celltype);
  vtx_R.resize(FLAGS_celltype);
  vtx_U.resize(FLAGS_celltype);
  vtx_D.resize(FLAGS_celltype);
  vtx_LU.resize(FLAGS_celltype);
  vtx_LD.resize(FLAGS_celltype);
  vtx_RU.resize(FLAGS_celltype);
  vtx_RD.resize(FLAGS_celltype);
  
  str_t.reserve((FLAGS_Tm-FLAGS_Ts)*FLAGS_H + 2);
  if(FLAGS_edges){
    str_delauny.reserve((FLAGS_Tm-FLAGS_Ts)*FLAGS_H + 2);
    str_voronoi.reserve((FLAGS_Tm-FLAGS_Ts)*FLAGS_H + 2);
  }
  for(auto p : str_points){
    p.reserve(FLAGS_Tm*FLAGS_H + 2);
  }
  
  
  InitPoints();
  InitPeriodicPoints();
  
  
  END;
}






void Simulation::InitPeriodicPoints(){
  BGN;
  
  //if(FLAGS_periodic>0){
  
  ////for(size_t i=0;i<points.size();++i){  
  //for(vector_Point::iterator p=points[i].begin(),pe=points[i].end();p!=pe;++p){
  //const size_t p_size=points[i].size();
  
  //points_R[i].reserve(p_size);
  //points_L[i].reserve(p_size);
  //points_R[i].push_back(MyPointC2(p->x()+FLAGS_Xp,p->y(),p->color(),p->radius(),0.0,0.0));
  //points_L[i].push_back(MyPointC2(p->x()-FLAGS_Xp,p->y(),p->color(),p->radius(),0.0,0.0));
  
  //if(FLAGS_periodic>1){
  //points_U[i].reserve(p_size);
  //points_D[i].reserve(p_size);
  //points_RU[i].reserve(p_size);
  //points_LD[i].reserve(p_size);
  //points_RU[i].reserve(p_size);
  //points_LD[i].reserve(p_size);
  //points_U[i].push_back(MyPointC2(p->x(),p->y()+FLAGS_Yp,p->color(),p->radius(),0.0,0.0));
  //points_D[i].push_back(MyPointC2(p->x(),p->y()-FLAGS_Yp,p->color(),p->radius(),0.0,0.0));
  //points_RU[i].push_back(MyPointC2(p->x()+FLAGS_Xp,p->y()+FLAGS_Yp,p->color(),p->radius(),0.0,0.0));
  //points_RD[i].push_back(MyPointC2(p->x()+FLAGS_Xp,p->y()-FLAGS_Yp,p->color(),p->radius(),0.0,0.0));
  //points_LU[i].push_back(MyPointC2(p->x()-FLAGS_Xp,p->y()+FLAGS_Yp,p->color(),p->radius(),0.0,0.0));
  //points_LD[i].push_back(MyPointC2(p->x()-FLAGS_Xp,p->y()-FLAGS_Yp,p->color(),p->radius(),0.0,0.0));
  //}
  
  //}
  //points_all[i].insert(points_all[i].end(),points_R[i].begin(),points_R[i].end());
  //points_all[i].insert(points_all[i].end(),points_L[i].begin(),points_L[i].end());
  //points_all[i].insert(points_all[i].end(),points_U[i].begin(),points_U[i].end());
  //points_all[i].insert(points_all[i].end(),points_D[i].begin(),points_D[i].end());
  //points_all[i].insert(points_all[i].end(),points_RU[i].begin(),points_RU[i].end());
  //points_all[i].insert(points_all[i].end(),points_RD[i].begin(),points_RD[i].end());
  //points_all[i].insert(points_all[i].end(),points_LU[i].begin(),points_LU[i].end());
  //points_all[i].insert(points_all[i].end(),points_LD[i].begin(),points_LD[i].end());
  
  //}
  //}
  
  END;
}






void Simulation::InitPoints(){
  BGN;
  
  for(int i=0; i<FLAGS_celltype;++i){
    if(FLAGS_I!=string("")&&param.N[i]>0){
      stringstream ss;
      ss<<FLAGS_I<<"/I"<<i;
      ifstream fin(ss.str().c_str());
      MyPointC2 p;
      while(fin >> p) {
        vtx[i].push_back(dt.insert(p));
        
      }
    }
    else if(FLAGS_circle) InitPointsCircle(i);
    else InitPoints(vtx[i],param.N[i],param.L[i],param.U[i],i,FLAGS_r);
    
  }
  
  //if(FLAGS_mode=="1"){
  ////out_rangeな座標の要素を削除
  //for(vector_Point_2d::iterator pp=points.begin(),ppe=points.end();pp!=ppe;++pp){
  //vector_Point::iterator end_it = remove_if( pp->begin(), pp->end(), out_range);
  //pp->erase(end_it,pp->end());
  //}
  //}
  //points_all=points;
  END;
}









void Simulation::InitPoints(vector_Vertex_handle& vtx,const int& N,const double& L, const double& U, int t, double r)
{
  BGN;
  
  
  
  static bool firstcell=true;
  
  if(N>0){
    int n=0;//成功したセルの数
    
    
    if(firstcell){//全ての細胞型を通した、最初の細胞
      firstcell=false;
      const double X=dNMrand(-0.5*FLAGS_Xp + FLAGS_r,0.5*FLAGS_Xp);
      const double Y=dNMrand(FLAGS_L,FLAGS_U);
      MyPointC2 newpoint(X,Y,0,FLAGS_r,0.0,0.0);
      vtx.push_back(dt.insert(newpoint));
      ++n;
    }
    else{
      
      int n=0;//成功したセルの数
      int x=0;//失敗回数
      while(n<N){
        
        const double X=dNMrand(-0.5*FLAGS_Xp + r*FLAGS_inispace ,0.5*FLAGS_Xp);
        const double Y=dNMrand(L,U);
        const MyPointC2 new_point(X,Y,t,r,0.0,0.0);
        const MyPointC2 nearest(dt.nearest_vertex(new_point)->point());
        bool flag(true);//球が重ならないか示すフラグ
        if(CGAL::squared_distance(Point(nearest),Point(new_point))
          <FLAGS_inispace*(nearest.radius()+new_point.radius())*(nearest.radius()+new_point.radius()))
        {flag=false;++x;}
        if(flag){
          dt.insert(new_point);
          x=0;
          ++n;
        }
        
        if(x > 1000) {//点の挿入に1000回以上連続で失敗したら。
          if(FLAGS_fill) break;//FLAGS_fill=true なら詰められるだけ詰める
          else cerr << "Initialization failed in Simulation::InitPoints()" << endl; abort();
          
        }
        
      }
    }
    
    END;
  }
}



//Cellを円盤状に配置する
void Simulation::InitPointsCircle(int i){
  
  static bool firstcell=true;
  if(param.N[i]>0){
    
    int n=0;//成功したセルの数
    int x=0;//失敗回数
    
    if(firstcell){//全ての細胞型を通した、最初の細胞
      firstcell=false;
      const double r=param.U[i]*zP1rand();
      const double pi=2.0*myrand::PI*Zp1rand();
      const double X=r*cos(pi);
      const double Y=r*sin(pi);
      MyPointC2 newpoint(X,Y,0,FLAGS_r,0.0,0.0);
      vtx[i].push_back(dt.insert(newpoint));
      ++n;
      
    }
    else{
      
      while(n<param.N[i]){
        
        const double r=param.U[i]*zP1rand();
        const double pi=2.0*myrand::PI*Zp1rand();
        const double X=r*cos(pi);
        const double Y=r*sin(pi);
        const MyPointC2 new_point(X,Y,i,FLAGS_r,0.0,0.0);
        const MyPointC2 nearest(dt.nearest_vertex(new_point)->point());
        bool flag(true);//球が重ならないか示すフラグ
        if(CGAL::squared_distance(Point(nearest),Point(new_point))<FLAGS_inispace*(nearest.radius()+new_point.radius())*(nearest.radius()+new_point.radius())){flag=false;++x;}
        
        
        if(flag){
          vtx[i].push_back(dt.insert(new_point));
          x=0;
          ++n;
        }
        
        if(x > 100) {//点の挿入に1000回以上連続で失敗したら。
          if(FLAGS_fill) break;//FLAGS_fill=true なら詰められるだけ詰める
          else cerr << "Initialization failed in Simulation::InitPoints()" << endl; abort();
        }
        
      }
    }
  }
}












//点を１つだけ入れる
void Simulation::InitFirstCell(){
  BGN;
  
  const int number_of_celltype=4;
  vtx.resize(number_of_celltype);
  

  int i=Prob(param.P);
  const double r0=param.r[i]*zP1rand();
  const double pi=2.0*myrand::PI*Zp1rand();
  
  const double X=r0*cos(pi);
  const double Y=r0*sin(pi);
  
  vtx[i].push_back(dt.insert_first(MyPointC2(X,Y,i,param.r[i],.0,.0)));
  
  END;
}
























void Simulation::FoutEdges(){
  BGN;
  string filename;
  ofstream fout;
  
  
  if(FLAGS_edges){
    filename = FLAGS_o + string("/delaunay.txt");
    fout.open(filename.c_str(),ios::trunc);
    fout << "t x1 y1 x2 y2" << endl;
    fout << delauny.str();
    fout.close();
    
    filename = FLAGS_o + string("/voronoi.txt");
    fout.open(filename.c_str(),ios::trunc);
    fout << "t x1 y1 x2 y2" << endl;
    fout << voronoi.str();
    fout.close();
  }
  
  
 
  
  
  END;
}












void Simulation::CacheEdges(){
  BGN;
  
  str_t.push_back(tos(t*FLAGS_H,5,'0'));
  
  
  //ドロネー、ボロノイ辺を記録する場合
  
  for(DT::Edge_iterator e=dt.finite_edges_begin();e!=dt.finite_edges_end();++e){
    //delaunay_edges
    double T=t*FLAGS_H;
    delauny << T << " " << dt.segment(e) << endl; 
    
    //voronoi edges
    CGAL::Object obj;
    obj=dt.dual(e);
    K::Segment_2 segment;
    K::Ray_2 ray;
    K::Line_2 line;
    if (assign(segment, obj)) {
      voronoi << T << " " << segment << endl;
    } else if (assign(ray, obj)) {
      voronoi << T << " " << ray << endl;
    }else if(assign(line, obj)){
      voronoi << T << " " << line << endl;
    }
  }
  
  
  
  
  END;
}














void Simulation::ShowVerticesAll(){
  BGN;
  
  
  cout << "###############" << endl;
  size_t n=vtx.size();
  for(size_t i=0;i<n;++i){
    cout << "I"<<i<<" Cells"<< endl;
    for(auto v : vtx[i]) cout << MyPointC2((v)->point()) << endl;
  }
  cout << "###############" << endl;
  END;
}







void Simulation::ShowVerticesAllNum(){
  BGN;
  cout << endl;
  cout << "###############" << endl;
  size_t n=vtx.size();
  for(size_t i=0;i<n;++i){
    cout << "I"<<i<<" Cells "<< vtx[i].size() << endl;
  }
  cout << "###############" << endl;
  END;
}









void Simulation::UpdatePosition(){
  BGN;
  
  for(auto vv : vtx){
    for(auto v : vv){
      
      
      double X=(v)->point().x() + FLAGS_H*(v)->point().vx();
      double Y=(v)->point().y() + FLAGS_H*(v)->point().vy();
      
      if(FLAGS_periodic>0.0){
        if     (0.5*FLAGS_Xp <= X) X = X - FLAGS_Xp;
        else if(X < -0.5*FLAGS_Xp) X = X + FLAGS_Xp;
      }
      
      if(FLAGS_periodic>1.0){
        if(FLAGS_Yp<=Y) Y = Y - FLAGS_Yp;
        else if(Y<0.0)  Y = Y + FLAGS_Yp;
      }
      
      dt.move_if_no_collision(v,Point(MyPointC2(X,Y,(v)->point().color(),(v)->point().radius(),0.0,0.0)));
      
      //dt.move(*v,Point(MyPointC2(X,Y,(*v)->point().color(),(*v)->point().radius(),(*v)->point().vx(),(*v)->point().vy())));
    }
  }
  
  
  //周期点を動かす
  for(size_t i=0;i<vtx.size();++i){
    
    vector_Vertex_handle::iterator v=vtx[i].begin();
    vector_Vertex_handle::iterator ve=vtx[i].end();
    vector_Vertex_handle::iterator v_L;
    vector_Vertex_handle::iterator v_R;
    vector_Vertex_handle::iterator v_U;
    vector_Vertex_handle::iterator v_D;
    vector_Vertex_handle::iterator v_LU;
    vector_Vertex_handle::iterator v_RU;
    vector_Vertex_handle::iterator v_LD;
    vector_Vertex_handle::iterator v_RD;
    
    if(FLAGS_periodic>0){
      v_L=vtx_L[i].begin();
      v_R=vtx_R[i].begin();
    }
    else if(FLAGS_periodic>1){
      v_U=vtx_U[i].begin();
      v_D=vtx_D[i].begin();
      v_LU=vtx_LU[i].begin();
      v_RU=vtx_RU[i].begin();
      v_LD=vtx_LD[i].begin();
      v_RD=vtx_RD[i].begin();
    }
    
    while(v!=ve){
      const double x=(*v)->point().x();
      const double y=(*v)->point().y();
      const int    c=(*v)->point().color();
      const double r=(*v)->point().radius();
      const double vx=(*v)->point().vx();
      const double vy=(*v)->point().vy();
      const double R=x+FLAGS_Xp;
      const double L=x-FLAGS_Xp;
      const double U=y+FLAGS_Yp;
      const double D=y-FLAGS_Yp;
      if(FLAGS_periodic>0){
        
        dt.move_if_no_collision(*(v_L++),Point(MyPointC2(L,y,c,r,vx,vy)));
        dt.move_if_no_collision(*(v_R++),Point(MyPointC2(R,y,c,r,vx,vy)));
      }
      else if(FLAGS_periodic>1){
        
        dt.move_if_no_collision(*(v_U++),Point(MyPointC2(x,U,c,r,vx,vy)));
        dt.move_if_no_collision(*(v_D++),Point(MyPointC2(x,D,c,r,vx,vy)));
        dt.move_if_no_collision(*(v_LU++),Point(MyPointC2(L,U,c,r,vx,vy)));
        dt.move_if_no_collision(*(v_LD++),Point(MyPointC2(L,D,c,r,vx,vy)));
        dt.move_if_no_collision(*(v_RU++),Point(MyPointC2(R,U,c,r,vx,vy)));
        dt.move_if_no_collision(*(v_RD++),Point(MyPointC2(R,D,c,r,vx,vy)));
      }
      
      ++v;
    }
  }
  
  
  END;
}













void Simulation::Ablation(){
  BGN;
  
  {
    
    for(size_t type=0;type<vtx.size();++type){
      
      vector<int> e;
      e.reserve(vtx[type].size());
      
      //削除する頂点の番号eを取得
      for(vector<DT::Vertex_handle>::iterator v=vtx[type].begin(),ve=vtx[type].end();v!=ve;++v){
        
        
        //if(CGAL::squared_distance((*v)->point(),DT::Point(FLAGS_x1,FLAGS_y1)) <=  FLAGS_d1*FLAGS_d1)
        //e.push_back(distance(vtx[type].begin(), v));
        if(
          FLAGS_x1-FLAGS_d1*0.5 < (*v)->point().x() && (*v)->point().x() < FLAGS_x1+FLAGS_d1*0.5 &&
          FLAGS_y1-FLAGS_d1*0.5 < (*v)->point().y() && (*v)->point().y() < FLAGS_y1+FLAGS_d1*0.5)
        { 
          e.push_back(distance(vtx[type].begin(), v));
        }
      }
      
      //eを逆順に巡り要素を削除してゆく
      //vtxはvectorなので前から削除すると
      //それより後ろのイテレータが無効になってしまう
      for(vector<int>::reverse_iterator i=e.rbegin(),ie=e.rend();i!=ie;++i){
        dt.remove(vtx[type][*i]);
        vtx[type].erase(vtx[type].begin()+*i);
        if(FLAGS_periodic>0){
          dt.remove(vtx_R[type][*i]);
          dt.remove(vtx_L[type][*i]);
          vtx_R[type].erase(vtx_R[type].begin()+*i);
          vtx_L[type].erase(vtx_L[type].begin()+*i);
        }
        if(FLAGS_periodic>1){
          dt.remove(vtx_U[type][*i]);
          dt.remove(vtx_D[type][*i]);
          dt.remove(vtx_RU[type][*i]);
          dt.remove(vtx_LU[type][*i]);
          dt.remove(vtx_RD[type][*i]);
          dt.remove(vtx_LD[type][*i]);
          vtx_U[type].erase(vtx_U[type].begin()+*i);
          vtx_D[type].erase(vtx_D[type].begin()+*i);
          vtx_RU[type].erase(vtx_RU[type].begin()+*i);
          vtx_LU[type].erase(vtx_LU[type].begin()+*i);
          vtx_RD[type].erase(vtx_RD[type].begin()+*i);
          vtx_LD[type].erase(vtx_LD[type].begin()+*i);
          
        }
      }
      
      
    }
    
  }   
  
  END;
}


void Simulation::Apoptosis(){
  BGN;
  
  size_t n=param.Pd.size();
  for(size_t i=0;i<n;++i) if(param.Pd[i] > 0.0) Apoptosis(i,param.Pd[i]);
  END;
}


void Simulation::Apoptosis(int type, double Pd){
  BGN;
  
  //削除する細胞を決める
  vector<int> e;
  for(size_t i=0;i<vtx[type].size();++i){
    if(Prob1(FLAGS_H*Pd)) e.push_back(i);
  }
  
  for(vector<int>::reverse_iterator i=e.rbegin(),ie=e.rend();i!=ie;++i){
    dt.remove(vtx[type][*i]);
    vtx[type].erase(vtx[type].begin()+*i);
    if(FLAGS_periodic>0){
      dt.remove(vtx_R[type][*i]);
      dt.remove(vtx_L[type][*i]);
      vtx_R[type].erase(vtx_R[type].begin()+*i);
      vtx_L[type].erase(vtx_L[type].begin()+*i);
    }
    if(FLAGS_periodic>1){
      dt.remove(vtx_U[type][*i]);
      dt.remove(vtx_D[type][*i]);
      vtx_U[type].erase(vtx_U[type].begin()+*i);
      vtx_D[type].erase(vtx_D[type].begin()+*i);
      
      dt.remove(vtx_RU[type][*i]);
      dt.remove(vtx_RD[type][*i]);
      vtx_RU[type].erase(vtx_RU[type].begin()+*i);
      vtx_RD[type].erase(vtx_RD[type].begin()+*i);
      
      dt.remove(vtx_LU[type][*i]);
      dt.remove(vtx_LD[type][*i]);
      vtx_LU[type].erase(vtx_LU[type].begin()+*i);
      vtx_LD[type].erase(vtx_LD[type].begin()+*i);
    }
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
  }
  
  
}
