#include "simulation.h"
#include "utility.h"
#include "random.h"
//ほげ


void Simulation::MakeInitFile(){
  BGN;
  Init();
  
  
  //InitForMakeInitFile();
  
  
  

  //FLAGS_H=0.001;
  //FLAGS_Tm=10;
  //#########################################################
  //点の位置を更新する
  //const double Tm=FLAGS_Tm;
  //const double H=FLAGS_H;
  //さぶい
  
  //ほげほげ
  
  
  if(FLAGS_v) ShowVerticesAllNum();
  int tmax = static_cast<int>(FLAGS_Tm/FLAGS_H);
  bool ADD=true;
  
  //時間発展
  for(t=0;t<=tmax;++t){
    
    const double T=t*FLAGS_H;
    if(FLAGS_v){
      cout << "\r";
      cout << "T = " << T << flush;
    }
    
    
    
    //時間を進めつつ、
    //ある時間間隔でCellを追加してゆく
    if(ADD && t%int(FLAGS_Fs/FLAGS_H)==0){
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
        Vertex_handle v=dt.nearest_vertex(newpoint);
        if(CGAL::squared_distance(v->point(),Point(newpoint))>4*FLAGS_r*FLAGS_r){
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
  
  DevidePointsStripe();

  FoutVertices();
  FoutCSVdata();
  END;
}

void Simulation::InitForMakeInitFile(){
  BGN;
  const int n_cell_type=7;
  points.resize(n_cell_type);
  points_L.resize(n_cell_type);
  points_R.resize(n_cell_type);
  points_U.resize(n_cell_type);
  points_D.resize(n_cell_type);
  points_LU.resize(n_cell_type);
  points_LD.resize(n_cell_type);
  points_RU.resize(n_cell_type);
  points_RD.resize(n_cell_type);
  points_all.resize(n_cell_type);
  vtx.resize(n_cell_type);
  vtx_L.resize(n_cell_type);
  vtx_R.resize(n_cell_type);
  vtx_U.resize(n_cell_type);
  vtx_D.resize(n_cell_type);
  vtx_LU.resize(n_cell_type);
  vtx_LD.resize(n_cell_type);
  vtx_RU.resize(n_cell_type);
  vtx_RD.resize(n_cell_type);
  //str_points.resize(n_cell_type);
  //str_t.reserve(FLAGS_Tm*FLAGS_H + 2);
  //if(FLAGS_edges){
  //str_delauny.reserve(FLAGS_Tm*FLAGS_H + 2);
  //str_voronoi.reserve(FLAGS_Tm*FLAGS_H + 2);
  //}
  //for(auto p : str_points){
  //p.reserve(FLAGS_Tm*FLAGS_H + 2);
  //}
  
  
  
  InitRand();
  if(!FLAGS_makeinit) InitPoints();
  else InitPointsForMakeInitFile();
  InitPeriodicPoints();
  InitTrianglation();
  END;
}








//点を１つ入れる
void Simulation::InitPointsForMakeInitFile(){
  BGN;
  
  
  //InitPoints(points[0],10000,FLAGS_L,FLAGS_U,0,FLAGS_r);

  
  const double X=dNMrand(-0.5*FLAGS_Xp + FLAGS_r,0.5*FLAGS_Xp);
  const double Y=dNMrand(FLAGS_L,FLAGS_U);
  points[0].push_back(MyPointC2(X,Y,0,FLAGS_r,0.0,0.0));

  
  points_all=points;
  END;
}



void Simulation::DevidePointsStripe(){
  BGN;
  
  vector_Vertex_handle_2d vtx_orig;
  vtx.swap(vtx_orig);
  
  const size_t n_cell_type=7;
  vtx.resize(n_cell_type);
  cout << endl;
  cout << "BEFORE" << endl;
  ShowVerticesAllNum();
  
  
  for(auto vv:vtx_orig){
    for(auto v:vv){
      

      for(size_t i=0;i<n_cell_type;++i){
        if(param.L[i] <= v->point().y() && v->point().y() < param.U[i]){
          v->point().color()=i;
          vtx[i].push_back(v);
        }
      } 
    }
  }
  cout << "AFTER" << endl;
  ShowVerticesAllNum();
  
  END;
}


void Simulation::DevidePointsCircle(){
  BGN;
  
  vector_Vertex_handle_2d vtx_orig;
  swap(vtx,vtx_orig);
  
  const size_t n_cell_type=7;
  vtx.resize(n_cell_type);
  cout << "BEFORE" << endl;
  ShowVerticesAllNum();
  
  
  for(auto vv:vtx_orig){
    for(auto v:vv){
      
      const double sqd=CGAL::squared_distance(v->point(),Point(0,0));
      for(size_t i=0;i<n_cell_type;++i){
        if(param.L[i]*param.L[i] <= sqd && sqd < param.U[i]*param.U[i]){
          v->point().color()=i;
          vtx[i].push_back(v);
        }
      } 
    }
  }
  cout << "AFTER" << endl;
  ShowVerticesAllNum();
  
  END;
}



void Simulation::InitPointsRegular(){
  BGN; 
  const double sqrt3=sqrt(3);
  double Y1= 0.0;
  double Y2= FLAGS_r*sqrt3;
  
  
  
  vector<MyPointC2> p6;
  p6.reserve(6);
  p6.push_back(MyPointC2(0.0,  Y1,1,FLAGS_r,0.0,0.0));
  p6.push_back(MyPointC2(2*FLAGS_r,    Y1,1,FLAGS_r,0.0,0.0));
  p6.push_back(MyPointC2(4*FLAGS_r,  Y1,1,FLAGS_r,0.0,0.0));
  p6.push_back(MyPointC2(FLAGS_r,Y2,1,FLAGS_r,0.0,0.0));
  p6.push_back(MyPointC2(3*FLAGS_r,Y2,1,FLAGS_r,0.0,0.0));
  p6.push_back(MyPointC2(5*FLAGS_r,Y2,1,FLAGS_r,0.0,0.0));
  
  //
  FLAGS_Xp=6.0*FLAGS_r;
  points[1].reserve(FLAGS_U1*6);
  for(int i=0;i<FLAGS_U1;++i){
    double X= -2.5*FLAGS_r;
    double Y=0.5*sqrt3*FLAGS_r + 2*sqrt3*FLAGS_r*i;
    
    for(auto p : p6){
      p.x()+=X;
      p.y()+=Y;
    }
    points[1].insert(points[1].end(),p6.begin(),p6.end());
  }
  
  points[2].reserve(points[1].size());
  for(auto p : points[1]){
    points[2].push_back(MyPointC2(p.x(),-p.y(),2,FLAGS_r,0.0,0.0));
  }
  
  
  
  
  END;
}




void Simulation::Run(){
  BGN;
  
  //細胞初期化ファイルを作る場合
  if(FLAGS_makeinit){
    MakeInitFile();
    return;
  }
  
  Init();
  
  if(FLAGS_v){
    //ShowPointsAll();
    //ShowPointsAllNum();
    ShowVerticesAllNum();
  }
  
  //if(FileExist(FLAGS_o.c_str())) remove_all(FLAGS_o.c_str());
  //mkdir(FLAGS_o.c_str(), 0777);
  
  //OutFlags(FLAGS_o+string("/gflag.txt"));
  
  //#########################################################
  //点の位置を更新する
  //ShowVerticesAll();
  int tmax = static_cast<int>(FLAGS_Tm/FLAGS_H);
  for(t=0;t<=tmax;++t){
    
    const double T=t*FLAGS_H;
    
    
    if(FLAGS_v){
      cout << "\r";
      cout << "T = " << T << flush;
    }
    
    if(T==FLAGS_Tabl) Ablation();
    if(T>FLAGS_Td) Apoptosis();
    CalcVelocity();
    
    
    //if(t==0||t==tmax) ShowVerticesAll();
    //ShowVerticesAll();
    
    //点の位置をバッファに記録
    if(FLAGS_anime && (t%int(FLAGS_Fs/FLAGS_H)==0 || t==tmax)){
      
    }
    
    //細胞の位置や速度を記録
    if((T >= FLAGS_Ts) && (t%int(FLAGS_Fs/FLAGS_H)==0 || t==tmax)){
      if(FLAGS_anime) CacheVertices();
      SaveVerticesToString();
      SaveVx();
    }
    
    
    UpdatePosition();
  }
  
  //バッファをファイルに出力
  //点の数が多い場合や時間が長い場合は
  //適当に出力しながらヤラないとメモリ足らなくなるかも
  if(FLAGS_makeinit) mkdir(FLAGS_I.c_str(), 0777);
  else mkdir(FLAGS_o.c_str(), 0777);
  AnalyzeMovement();
  if(FLAGS_outvtx)FoutVertices();
  if(FLAGS_anime) FoutCSVdata();
  if(FLAGS_edges) OutVertices();
  END;
}










//記録した細胞の移動を解析する
void Simulation::AnalyzeMovement()
{
  BGN;
  
  
  const size_t A1=1;
  const size_t A2=2;
  
  double mean_speeds_A1=0.0;
  double mean_speeds_A2=0.0;
  
  
  if(vtx[A1].size()>0) mean_speeds_A1=mean(Speed_A1);
  if(vtx[A2].size()>0) mean_speeds_A2=mean(Speed_A2);
  
  
  
  
  double mean_Vx_A1_all=0.0;
  double mean_Vx_A2_all=0.0;
  double mean_Vx_A1_inner=0.0;
  double mean_Vx_A2_inner=0.0;
  double mean_Vx_A1_bound=0.0;
  double mean_Vx_A2_bound=0.0;
  
  double median_Vx_A1_all=0.0;
  double median_Vx_A2_all=0.0;
  double median_Vx_A1_inner=0.0;
  double median_Vx_A2_inner=0.0;
  double median_Vx_A1_bound=0.0;
  double median_Vx_A2_bound=0.0;
  
  
  if(vtx[A1].size()>0 && vtx[A2].size()>0){
    mean_Vx_A1_all =mean(Vx_A1_all);
    mean_Vx_A2_all=mean(Vx_A2_all);    
    mean_Vx_A1_inner =mean(Vx_A1_inner);
    mean_Vx_A2_inner=mean(Vx_A2_inner);
    mean_Vx_A1_bound=mean(Vx_A1_bound);
    mean_Vx_A2_bound=mean(Vx_A2_bound);
    
    median_Vx_A1_all =median(Vx_A1_all);
    median_Vx_A2_all=median(Vx_A2_all);    
    median_Vx_A1_inner =median(Vx_A1_inner);
    median_Vx_A2_inner=median(Vx_A2_inner);
    median_Vx_A1_bound=median(Vx_A1_bound);
    median_Vx_A2_bound=median(Vx_A2_bound);
    
  }
  
  string filename = FLAGS_o + string("/velocity.csv");
  ofstream fout(filename.c_str(),ios::trunc);
  
  fout <<"speed_A1,speed_A2,mean_Vx_A1_all,mean_Vx_A2_all,mean_Vx_A1_inner,mean_Vx_A2_inner,mean_Vx_A1_bound,mean_Vx_A2_bound,";
  fout <<"median_Vx_A1_all,median_Vx_A2_all,median_Vx_A1_inner,median_Vx_A2_inner,median_Vx_A1_bound,median_Vx_A2_bound,";
  fout <<"Asame1,Asame2,Adiff1,Rsame1,Rsame2,Rdiff1,Pd1,Pd2,x1,y1,d1,repli"<<endl
  << mean_speeds_A1 << ","
  << mean_speeds_A2 << ","
  << mean_Vx_A1_all  << ","
  << mean_Vx_A2_all  << ","  
  << mean_Vx_A1_inner  << ","
  << mean_Vx_A2_inner  << ","
  << mean_Vx_A1_bound << ","
  << mean_Vx_A2_bound<< ","
  << median_Vx_A1_all  << ","
  << median_Vx_A2_all  << ","  
  << median_Vx_A1_inner  << ","
  << median_Vx_A2_inner  << ","
  << median_Vx_A1_bound << ","
  << median_Vx_A2_bound << ","
  << FLAGS_Pd1 << ","
  << FLAGS_Pd2 << ","
  << FLAGS_x1 << ","
  << FLAGS_y1 << ","
  << FLAGS_d1 << ","
  << FLAGS_repli << endl;
  fout.close();
  
  END;
}















void Simulation::CacheVertices(){
  BGN;
  
  for(auto vv : vtx){
    for(auto v : vv){
      const auto &p=v->point();
      csv_data 
      << t*FLAGS_H << "," 
      << p.color() << ","
      << p.x() << ","
      << p.y() << ","
      << p.vx() << ","
      << p.vy() << ","
      << p.radius() << ",";
      
      DT::Vertex_circulator vi,v0;
      vi=v0=dt.incident_vertices(v);//vと隣接する頂点をめぐるサーキュレータ
      int counter=0;
      do{
        if(!dt.is_infinite(vi)){//viが無限遠にある頂点でないなら
          
          const auto &pi=vi->point();
          const auto X=(pi.x()-p.x());
          const auto Y=(pi.y()-p.y());
          
          csv_data 
          << sqrt(X*X+Y*Y) << ","
          << CalcRadian(X,Y) << ","
          << pi.color() << ",";
          ++counter;
        }
        if(counter==8) break;//隣接細胞が多すぎるのは例外的なので打ち切り
      }while(++vi!=v0);
      
      int n=counter;
      for(;counter<8;++counter){
        csv_data << ",,,";
      }
      csv_data << n << ",";
      csv_data
      << FLAGS_mode << ","
      //<< FLAGS_Asame1 << ","
      //<< FLAGS_Asame2 << ","
      //<< FLAGS_Adiff1 << ","
      //<< FLAGS_Rsame1 << ","
      //<< FLAGS_Rsame2 << ","
      //<< FLAGS_Rdiff1 << ","
      << FLAGS_Pd1 << ","
      << FLAGS_Pd2 << ","
      << FLAGS_x1 << ","
      << FLAGS_y1 << ","
      << FLAGS_d1 << ","
      << FLAGS_repli
      << endl;
      
    }
  }
  
  //シミュレーションパラメターの出力
  
  
  
  END;
}



void Simulation::CalcVelocity()
{
  
  BGN;
  
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
  END;
}


void Simulation::FoutCSVdata(){
  BGN;
  
  string filename=FLAGS_o+string("/output.csv");
  
  ofstream fout(filename.c_str(),ios::trunc);
  fout << "t,col,x,y,vx,vy,r,r1,th1,c1,r2,th2,c2,r3,th3,c3,r4,th4,c4,r5,th5,c5,r6,th6,c6,r7,th7,c7,r8,th8,c8,adj";
  fout << ",mode,Asame1,Asame2,Adiff1,Rsame1,Rsame2,Rdiff1,Pd1,Pd2,x1,y1,d1,repli";
  fout << endl;
  fout << csv_data.str();
  fout.close();
  
  END;
}



void Simulation::FoutVertices(){
  BGN;
  //if(FLAGS_makeinit) mkdir(FLAGS_I.c_str(), 0777);
  //else mkdir(FLAGS_o.c_str(), 0777);
  
  size_t n=vtx.size();
  for(size_t i=0;i<n;++i){
    ostringstream oss; 
    if(FLAGS_makeinit) oss << FLAGS_I << "/I" << i;
    else               oss << FLAGS_o << "/I" << i;
    ofstream fout(oss.str().c_str());
    for(auto v:vtx[i]) fout << MyPointC2(v->point()) << endl;
    fout.close();
  }
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
  const int n_cell_type=7;
  points.resize(n_cell_type);
  points_L.resize(n_cell_type);
  points_R.resize(n_cell_type);
  points_U.resize(n_cell_type);
  points_D.resize(n_cell_type);
  points_LU.resize(n_cell_type);
  points_LD.resize(n_cell_type);
  points_RU.resize(n_cell_type);
  points_RD.resize(n_cell_type);
  points_all.resize(n_cell_type);
  vtx.resize(n_cell_type);
  vtx_L.resize(n_cell_type);
  vtx_R.resize(n_cell_type);
  vtx_U.resize(n_cell_type);
  vtx_D.resize(n_cell_type);
  vtx_LU.resize(n_cell_type);
  vtx_LD.resize(n_cell_type);
  vtx_RU.resize(n_cell_type);
  vtx_RD.resize(n_cell_type);
  
  //str_points.resize(n_cell_type);
  str_t.reserve((FLAGS_Tm-FLAGS_Ts)*FLAGS_H + 2);
  if(FLAGS_edges){
    str_delauny.reserve((FLAGS_Tm-FLAGS_Ts)*FLAGS_H + 2);
    str_voronoi.reserve((FLAGS_Tm-FLAGS_Ts)*FLAGS_H + 2);
  }
  for(auto p : str_points){
    p.reserve(FLAGS_Tm*FLAGS_H + 2);
  }
  
  
  
  InitRand();
  if(!FLAGS_makeinit) InitPoints();
  else InitPointsForMakeInitFile();
  InitPeriodicPoints();
  InitTrianglation();
  
  END;
}



void Simulation::InitPoints(){
  BGN;
  
  const size_t n=points.size();
  for(size_t i=0; i<n;++i){
    if(FLAGS_I!=string("")&&param.N[i]>0){
      stringstream ss;
      ss<<FLAGS_I<<"/I"<<i;
      ifstream fin(ss.str().c_str());
      MyPointC2 p;
      while(fin >> p) points[i].push_back(p);
    }
    else if(FLAGS_circle) InitPointsCircle(i);
    else InitPoints(points[i],param.N[i],param.L[i],param.U[i],i,FLAGS_r);
    
  }
  
  //if(FLAGS_mode=="1"){
  ////out_rangeな座標の要素を削除
  //for(vector_Point_2d::iterator pp=points.begin(),ppe=points.end();pp!=ppe;++pp){
  //vector_Point::iterator end_it = remove_if( pp->begin(), pp->end(), out_range);
  //pp->erase(end_it,pp->end());
  //}
  //}
  points_all=points;
  END;
}









void Simulation::InitPoints(vector<Point>& points,const int& N,const double& L, const double& U, int t, double r)
{
  BGN;
  if(N>0){
    
    int n=0;//成功したセルの数
    int x=0;//失敗回数
    while(n<N){
      
      const double X=dNMrand(-0.5*FLAGS_Xp + r*FLAGS_inispace ,0.5*FLAGS_Xp);
      const double Y=dNMrand(L,U);
      MyPointC2 new_point(X,Y,t,r,0.0,0.0);
      
      bool flag(true);//球が重ならないか示すフラグ
      for(vector<Point>::iterator p=points.begin();p!=points.end();++p){
        if(CGAL::squared_distance(*p,Point(new_point))<FLAGS_inispace*2*2*r*r){flag=false;++x;break;}
      }
      if(flag){
        points.push_back(new_point);
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



//Cellを円盤状に配置する
void Simulation::InitPointsCircle(int i){
  
  if(param.N[i]>0){
    
    int n=0;//成功したセルの数
    int x=0;//失敗回数
    while(n<param.N[i]){
      
      const double r=param.U[i]*zP1rand();
      const double pi=2.0*myrand::PI*Zp1rand();
      
      const double X=r*cos(pi);
      const double Y=r*sin(pi);
      MyPointC2 new_point(X,Y,i,FLAGS_r,0.0,0.0);
      
      bool flag(true);//球が重ならないか示すフラグ
      for(vector<Point>::iterator p=points[i].begin();p!=points[i].end();++p){
        if(CGAL::squared_distance(*p,Point(new_point))<FLAGS_inispace*2*2*FLAGS_r*FLAGS_r){flag=false;++x;break;}
      }
      if(flag){
        points[i].push_back(new_point);
        x=0;
        ++n;
      }
      
      if(x > 100) {//点の挿入に1000回以上連続で失敗したら。
        if(FLAGS_fill) break;//FLAGS_fill=true なら詰められるだけ詰める
        else cerr << "Initialization failed in Simulation::InitPoints()" << endl; abort();
      }
      
    }
    //out_rangeな座標の要素を削除
    vector_Point::iterator end_it = remove_if(points[i].begin(), points[i].end(), [&](Point p){
        const double d=CGAL::squared_distance(p,Point(0,0));
        if( d > param.U[i]*param.U[i] || d < param.L[i]*param.L[i]) return true;
        else return false;
      });
    points[i].erase(end_it,points[i].end());
    
  }
  
}







void Simulation::InitPeriodicPoints(){
  BGN;
  
  if(FLAGS_periodic>0){
    
    for(size_t i=0;i<points.size();++i){  
      for(vector_Point::iterator p=points[i].begin(),pe=points[i].end();p!=pe;++p){
        const size_t p_size=points[i].size();
        
        points_R[i].reserve(p_size);
        points_L[i].reserve(p_size);
        points_R[i].push_back(MyPointC2(p->x()+FLAGS_Xp,p->y(),p->color(),p->radius(),0.0,0.0));
        points_L[i].push_back(MyPointC2(p->x()-FLAGS_Xp,p->y(),p->color(),p->radius(),0.0,0.0));
        
        if(FLAGS_periodic>1){
          points_U[i].reserve(p_size);
          points_D[i].reserve(p_size);
          points_RU[i].reserve(p_size);
          points_LD[i].reserve(p_size);
          points_RU[i].reserve(p_size);
          points_LD[i].reserve(p_size);
          points_U[i].push_back(MyPointC2(p->x(),p->y()+FLAGS_Yp,p->color(),p->radius(),0.0,0.0));
          points_D[i].push_back(MyPointC2(p->x(),p->y()-FLAGS_Yp,p->color(),p->radius(),0.0,0.0));
          points_RU[i].push_back(MyPointC2(p->x()+FLAGS_Xp,p->y()+FLAGS_Yp,p->color(),p->radius(),0.0,0.0));
          points_RD[i].push_back(MyPointC2(p->x()+FLAGS_Xp,p->y()-FLAGS_Yp,p->color(),p->radius(),0.0,0.0));
          points_LU[i].push_back(MyPointC2(p->x()-FLAGS_Xp,p->y()+FLAGS_Yp,p->color(),p->radius(),0.0,0.0));
          points_LD[i].push_back(MyPointC2(p->x()-FLAGS_Xp,p->y()-FLAGS_Yp,p->color(),p->radius(),0.0,0.0));
        }
        
      }
      points_all[i].insert(points_all[i].end(),points_R[i].begin(),points_R[i].end());
      points_all[i].insert(points_all[i].end(),points_L[i].begin(),points_L[i].end());
      points_all[i].insert(points_all[i].end(),points_U[i].begin(),points_U[i].end());
      points_all[i].insert(points_all[i].end(),points_D[i].begin(),points_D[i].end());
      points_all[i].insert(points_all[i].end(),points_RU[i].begin(),points_RU[i].end());
      points_all[i].insert(points_all[i].end(),points_RD[i].begin(),points_RD[i].end());
      points_all[i].insert(points_all[i].end(),points_LU[i].begin(),points_LU[i].end());
      points_all[i].insert(points_all[i].end(),points_LD[i].begin(),points_LD[i].end());
      
    }
  }
  
  END;
}




void Simulation::InitTrianglation()
{
  BGN;
  
  //点を三角形分割する
  //全ての点、及び、全ての周期点を三角形分割する
  for(size_t i=0;i<points_all.size();++i){
    //cout << "Cell " << i << " " << points_all[i].size() << endl;
    dt.insert(points_all[i].begin(), points_all[i].end());
  }
  
  //各細胞のハンドルを保持しておく
  Vertex_handle h;
  
  
  for(size_t i=0;i<points.size();++i){
    
    size_t p_size=points[i].size();
    vtx[i].reserve(p_size);
    if(FLAGS_periodic>0){
      vtx_L[i].reserve(p_size);
      vtx_R[i].reserve(p_size);
      if(FLAGS_periodic>1){
        vtx_U[i].reserve(p_size);
        vtx_D[i].reserve(p_size);
        vtx_LU[i].reserve(p_size);
        vtx_LD[i].reserve(p_size);
        vtx_RU[i].reserve(p_size);
        vtx_RD[i].reserve(p_size);
      }
    }
    
    
    for(vector_Point::iterator p=points[i].begin(),pe=points[i].end();p!=pe;++p){
      
      h=dt.nearest_vertex(*p);
      vtx[i].push_back(h);
      
      if(FLAGS_periodic>0){
        for(vector<Point>::iterator p=points_L[i].begin();p!=points_L[i].end();++p){
          h=dt.nearest_vertex(*p);
          vtx_L[i].push_back(h);
        }
        
        for(vector<Point>::iterator p=points_R[i].begin();p!=points_R[i].end();++p){
          h=dt.nearest_vertex(*p);
          vtx_R[i].push_back(h);
        }
        
      }
      if(FLAGS_periodic>1){
        
        for(vector<Point>::iterator p=points_U[i].begin();p!=points_U[i].end();++p){
          h=dt.nearest_vertex(*p);
          vtx_U[i].push_back(h);
        }
        
        
        for(vector<Point>::iterator p=points_D[i].begin();p!=points_D[i].end();++p){
          h=dt.nearest_vertex(*p);
          vtx_D[i].push_back(h);
        }
        
        
        for(vector<Point>::iterator p=points_LU[i].begin();p!=points_LU[i].end();++p){
          h=dt.nearest_vertex(*p);
          vtx_LU[i].push_back(h);
        }
        
        
        for(vector<Point>::iterator p=points_LD[i].begin();p!=points_LD[i].end();++p){
          h=dt.nearest_vertex(*p);
          vtx_LD[i].push_back(h);
        }
        
        
        for(vector<Point>::iterator p=points_RU[i].begin();p!=points_RU[i].end();++p){
          h=dt.nearest_vertex(*p);
          vtx_RU[i].push_back(h);
        }
        
        
        for(vector<Point>::iterator p=points_RD[i].begin();p!=points_RD[i].end();++p){
          h=dt.nearest_vertex(*p);
          vtx_RD[i].push_back(h);
        }
        
      }
    }
  }
  
  MemClear(points);
  MemClear(points_L);  
  MemClear(points_R);  
  MemClear(points_U);  
  MemClear(points_D);  
  MemClear(points_RU);  
  MemClear(points_RD);  
  MemClear(points_LU);  
  MemClear(points_LD);  
  
  END;
}




















void Simulation::OutVertices(){
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
  
  
  //vector<string>::iterator t,te;
  //vector<string>::iterator a0,a7,a8a,a8p,a9,d,v;
  //vector<string>::iterator d,v;
  
  
  
  //a0=str_points[0].begin();
  //a7=str_points[1].begin();
  //a8a=str_points[2].begin();
  //a8p=str_points[3].begin();
  //a9=str_points[4].begin();
  //d=str_delauny.begin();
  //v=str_voronoi.begin();
  
  //size_t n=str_points.size();
  
  //t=str_t.begin();
  //te=str_t.end();
  //while(t!=te){
  
  //for(size_t i=0;i<n;++i){
  //if(param.N[i]>0){
  //filename = FLAGS_o + string("/I") + i string("_") + *(t) + string(".txt");
  //fout.open(filename.c_str(),ios::trunc);
  //fout << *(a0++);
  //fout.close();
  //}
  //}
  
  //if(FLAGS_edges){
  //filename = FLAGS_o + string("/delaunay_") + *(t) + string(".txt");
  //fout.open(filename.c_str(),ios::trunc);
  //fout << *(d++);
  //fout.close();
  
  //filename = FLAGS_o + string("/voronoi_") + *(t) + string(".txt");
  //fout.open(filename.c_str(),ios::trunc);
  //fout << *(v++);
  //fout.close();
  //}
  
  //++t;
  //}
  
  
  END;
}





void Simulation::SaveVerticesToString(){
  BGN;
  
  str_t.push_back(tos(t*FLAGS_H,5,'0'));
  
  //この時点における点の位置を記録
  //for(size_t i=0;i<vtx.size();++i){
  //stringstream ss;
  //vector<Vertex_handle>::iterator v,ve;
  //for(v=vtx[i].begin(), ve=vtx[i].end(); v!=ve; ++v){
  //ss << MyPointC2((*v)->point()) << "\n";
  //}
  //str_points[i].push_back(ss.str());
  //}
  
  //ドロネー、ボロノイ辺を記録する場合
  if(FLAGS_edges){
    
    //stringstream delauny;
    //stringstream voronoi;
    
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
    
    //str_delauny.push_back(delauny.str());
    //str_voronoi.push_back(voronoi.str());
    
    
  }
  
  END;
}



void Simulation::SaveVx(){
  BGN;
  
  const int A1=1;
  const size_t A2=2;
  
  
  for(auto v : vtx[A1])
  Speed_A1.push_back(sqrt(v->point().vx()*v->point().vx()+v->point().vy()*v->point().vy()));
  
  for(auto v : vtx[A2])
  Speed_A2.push_back(sqrt(v->point().vx()*v->point().vx()+v->point().vy()*v->point().vy()));
  
  
  
  set<DT::Vertex_handle> all_A1,boundary_A1;
  
  if(vtx[A1].size()>0 && vtx[A2].size()>0){
    //A2を巡回し、A1と隣接する細胞を同定する
    
    
    for(vector_Vertex_handle::iterator i=vtx[A2].begin(),ie=vtx[A2].end();i!=ie;++i){
      
      bool This_A2_is_bound=false;//境界にあるのか
      const double vx_A2=(*i)->point().vx();
      Vx_A2_all.push_back(vx_A2);
      
      DT::Vertex_circulator j,je;//A2[i]の隣接頂点jを巡る
      j=je=dt.incident_vertices(*i);
      do{
        if(j->point().color()==A1){//隣接頂点jがA1であるなら
          This_A2_is_bound=true;
          boundary_A1.insert(j);//A8と隣接するA1を保持しておく
        }
      }while(++j!=je);
      
      if(This_A2_is_bound)  Vx_A2_bound.push_back(vx_A2);
      else                     Vx_A2_inner.push_back(vx_A2);
    } 
  }
  
  //境界部のA1の速度を記録
  for(set<DT::Vertex_handle>::iterator i=boundary_A1.begin(),ie=boundary_A1.end();i!=ie;++i){
    Vx_A1_bound.push_back((*i)->point().vx());
  }
  
  //全てのA1の速度を記録
  for(vector<DT::Vertex_handle>::iterator i=vtx[A1].begin(),ie=vtx[A1].end();i!=ie;++i){
    Vx_A1_all.push_back((*i)->point().vx());
    all_A1.insert(*i);
  }
  
  //中心部のA1の速度を記録
  list<Vertex_handle> center_A1;
  back_insert_iterator<list<Vertex_handle> > v(center_A1);
  set_difference(
    all_A1.begin(),all_A1.end(),
    boundary_A1.begin(),boundary_A1.end(),
    v);
  
  
  for(list<Vertex_handle>::iterator i=center_A1.begin(),ie=center_A1.end();i!=ie;++i){
    Vx_A1_inner.push_back((*i)->point().vx()); 
  }
  
  
  
  
  END;
}






void Simulation::ShowPointsAll(){
  BGN;
  
  size_t n=points.size();
  for(size_t i=0;i<n;++i){
    cout << "I" << i << " Cells" << endl;
    copy(points[i].begin(), points[i].end(),
      std::ostream_iterator<MyPointC2>( std::cout, "\n"));
    cout << endl;
  }
  
  END;
}


void Simulation::ShowPointsAllNum(){
  BGN;
  size_t n=points.size();
  for(size_t i=0;i<n;++i){
    cout << "I"<<i<<"   Cells " << points[i].size() << endl;
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
  
  for(vector_Vertex_handle_2d::iterator vv=vtx.begin(),vve=vtx.end();vv!=vve;++vv){
    for(vector_Vertex_handle::iterator v=vv->begin(),ve=vv->end();v!=ve;++v){
      
      
      double X=(*v)->point().x() + FLAGS_H*(*v)->point().vx();
      double Y=(*v)->point().y() + FLAGS_H*(*v)->point().vy();
      
      if(FLAGS_periodic>0.0){
        if     (0.5*FLAGS_Xp <= X) X = X - FLAGS_Xp;
        else if(X < -0.5*FLAGS_Xp) X = X + FLAGS_Xp;
      }
      
      if(FLAGS_periodic>1.0){
        if(FLAGS_Yp<=Y) Y = Y - FLAGS_Yp;
        else if(Y<0.0)  Y = Y + FLAGS_Yp;
      }
      dt.move_if_no_collision(*v,Point(MyPointC2(X,Y,(*v)->point().color(),(*v)->point().radius(),0.0,0.0)));
      //dt.move(*v,Point(MyPointC2(X,Y,(*v)->point().color(),(*v)->point().radius(),(*v)->point().vx(),(*v)->point().vy())));
    }
  }
  
  //周期点を動かす
  for(size_t i=0;i<vtx.size();++i){
    vector_Vertex_handle::iterator v=vtx[i].begin();
    vector_Vertex_handle::iterator ve=vtx[i].end();
    vector_Vertex_handle::iterator v_L=vtx_L[i].begin();
    vector_Vertex_handle::iterator v_R=vtx_R[i].begin();
    vector_Vertex_handle::iterator v_U=vtx_U[i].begin();
    vector_Vertex_handle::iterator v_D=vtx_D[i].begin();
    vector_Vertex_handle::iterator v_LU=vtx_LU[i].begin();
    vector_Vertex_handle::iterator v_RU=vtx_RU[i].begin();
    vector_Vertex_handle::iterator v_LD=vtx_LD[i].begin();
    vector_Vertex_handle::iterator v_RD=vtx_RD[i].begin();
    
    if(FLAGS_periodic>0){
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
        
        //Point(MyPointC2(X,Y,(*v)->point().color(),(*v)->point().radius(),(*v)->point().vx(),(*v)->point().vy())))
        
        dt.move_if_no_collision(*(v_L++),Point(MyPointC2(L,y,c,r,vx,vy)));
        dt.move_if_no_collision(*(v_R++),Point(MyPointC2(R,y,c,r,vx,vy)));
        if(FLAGS_periodic>1){
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
