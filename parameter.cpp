#include "parameter.h"
#include "utility.h"

using namespace std;

//★パラメータ
Parameter::Parameter()
{
  
  FLAGS_A10=FLAGS_A01;
  FLAGS_A20=FLAGS_A02;
  FLAGS_A30=FLAGS_A03;
  FLAGS_A21=FLAGS_A12;
  FLAGS_A31=FLAGS_A13;
  FLAGS_A32=FLAGS_A23;
  
  P={FLAGS_P0,FLAGS_P1,FLAGS_P2,FLAGS_P3};
  N={FLAGS_N0,FLAGS_N1,FLAGS_N2,FLAGS_N3};
  U={FLAGS_U0,FLAGS_U1,FLAGS_U2,FLAGS_U3};
  L={FLAGS_L0,FLAGS_L1,FLAGS_L2,FLAGS_L3};
  Pd={FLAGS_Pd0,FLAGS_Pd1,FLAGS_Pd2,FLAGS_Pd3};
  A={
    {FLAGS_A00,FLAGS_A01,FLAGS_A02,FLAGS_A03},
    {FLAGS_A10,FLAGS_A11,FLAGS_A12,FLAGS_A13},
    {FLAGS_A20,FLAGS_A21,FLAGS_A22,FLAGS_A23},
    {FLAGS_A30,FLAGS_A31,FLAGS_A32,FLAGS_A33}
  };
  R={
    {FLAGS_R00,FLAGS_R01,FLAGS_R02,FLAGS_R03},
    {FLAGS_R10,FLAGS_R11,FLAGS_R12,FLAGS_R13},
    {FLAGS_R20,FLAGS_R21,FLAGS_R22,FLAGS_R23},
    {FLAGS_R30,FLAGS_R31,FLAGS_R32,FLAGS_R33}
  };
  fix={FLAGS_fix0,FLAGS_fix1,FLAGS_fix2,FLAGS_fix3};
  I={FLAGS_I0,FLAGS_I1,FLAGS_I2,FLAGS_I3}; 
  r={FLAGS_r0,FLAGS_r1,FLAGS_r2,FLAGS_r3};
};




std::string Parameter::OssFlags(){
  BGN;
  ostringstream oss;
  
  //計算時間
  oss << "-Tm=" << FLAGS_Tm << endl;
  //データ出力の間隔
  oss << "-Fs=" << FLAGS_Fs << endl;
  
  //回転の開始時間
  oss << "-Ta=" << FLAGS_Ta << endl;
  //細胞死の開始時間
  oss << "-Td=" << FLAGS_Td << endl;
  oss << "-H=" << FLAGS_H << endl;
  oss << "-edges=" << FLAGS_edges << endl;
  
  
  //周期境界条件
  oss << "-periodic=" << FLAGS_periodic << endl;
  oss << "-Xp=" << FLAGS_Xp << endl;
  oss << "-Yp=" << FLAGS_Yp << endl;
  //最適細胞間距離
  oss << "-r=" << FLAGS_r << endl;
  
  //最大細胞数
  oss << "-fill=" << FLAGS_fill << endl;
  
  const size_t n=N.size();
  for(size_t i=0;i<n;++i){
    oss << "-N" << i << "=" << N[i] << endl;
    oss << "-I" << i << "=" << I[i] << endl;
    oss << "-U" << i << "=" << U[i] << endl;
    oss << "-L" << i << "=" << L[i] << endl;
    oss << "-Pd" << i << "=" << Pd[i] << endl;
    for(size_t j=0;j<n;++j){
      oss << "-A" << i << j << "=" << A[i][j] << endl;
      oss << "-R" << i << j << "=" << R[i][j] << endl;
    }
    
  }
  
  
  
  END;
  return(oss.str());
  
}

void Parameter::OutFlags(std::string filepath){
  BGN;
  ofstream fout(filepath.c_str(),ios::trunc);
  fout << OssFlags() << endl;
  fout.close();
  END;
}










//計算時間
DEFINE_double(Tm,10,"Max calculation time");
DEFINE_double(H,0.001,"Step size of numerical calculation");//数値計算の刻み

//データ出力の間隔
DEFINE_double(Ts,8,"Starting time of data sampling ");
DEFINE_double(Fs,0.1,"Data sampling frequency");

//回転の開始時間
DEFINE_double(Ta,0.0,"Time beginning rotation");
//細胞死の開始時間
DEFINE_double(Td,10000.0,"Starting time of cell death");

//細胞除去を起こす時間
DEFINE_double(Tabl,10000.0,"Timing of ablation");
//細胞除去をする範囲
DEFINE_double(d1,0.0,"distance of ablation 1");
DEFINE_double(x1,0.0,"center x of ablation 1");
DEFINE_double(y1,0.0,"center y of ablation 1");


DEFINE_string(o,"out","Output directory name");
DEFINE_bool(edges,false,"Output delaunay & voronoi edges");
DEFINE_bool(velocity,false,"Output velocity");
DEFINE_bool(anime,false,"Output vertex position every time step for animation");
DEFINE_bool(outvtx,false,"Output vertices position");

DEFINE_double(noise,0.0,"Noise in velocity");
DEFINE_double(Tnoise,0.0,"Time to begin Noise");

//周期境界条件
DEFINE_int32(periodic,0,"priodic boundary");
DEFINE_double(Xp,6.283185,"Period of X in periodic boundary condition");
DEFINE_double(Yp,6.283185,"Period of Y in periodic boundary condition");


//最適細胞間距離
DEFINE_double(r,0.25,"Radius of a cell");
DEFINE_string(mode,"","Mode of simulation");
DEFINE_int32(repli,1,"Replicate ID number");
DEFINE_double(inispace,0.8,"adjust initial distance");
DEFINE_bool(circle,false,"Set initial position of cells circular.");
DEFINE_bool(makeinit,false,"Make init file");


//最大細胞数
DEFINE_bool(fill,true,"Fill up field with cells");
DEFINE_int32(N0,0,"Number of A0");
DEFINE_int32(N1,0,"Number of A1");
DEFINE_int32(N2,0,"Number of A2");
DEFINE_int32(N3,0,"Number of A3");


//初期位置の入力ファイル
DEFINE_string(I,"in","Input directory");
DEFINE_string(I0,"","Input file for I0");
DEFINE_string(I1,"","Input file for I1");
DEFINE_string(I2,"","Input file for I2");
DEFINE_string(I3,"","Input file for I3");




//初期位置の範囲
DEFINE_double(U , 0.0,"Upper boundary of Cells");
DEFINE_double(L , 0.0,"Lower boundary of Cells");


DEFINE_double(U0 , 0.0,"Upper boundary of I0");
DEFINE_double(U1 , 0.0,"Upper boundary of I1");
DEFINE_double(U2, 0.0,"Upper boundary of I2");
DEFINE_double(U3, 0.0,"Upper boundary of I3");


DEFINE_double(L0, 0.0,"Lower boundary of I0");
DEFINE_double(L1, 0.0,"Lower boundary of I1");
DEFINE_double(L2, 0.0,"Lower boundary of I2");
DEFINE_double(L3, 0.0,"Lower boundary of I3");


//細胞死率
DEFINE_double(Pd,0.0,"Probability of cell death per unit time for A0");
DEFINE_double(Pd0,0.0,"Probability of cell death per unit time for A0");
DEFINE_double(Pd1,0.0,"Probability of cell death per unit time for A1");
DEFINE_double(Pd2,0.0,"Probability of cell death per unit time for A2");
DEFINE_double(Pd3,0.0,"Probability of cell death per unit time for A3");


//位置の更新をするかどうか
DEFINE_bool(fix0,false,"Forcing the position of A0 fixed to be anchor");
DEFINE_bool(fix1,false,"Forcing the position of A1 fixed to be anchor");
DEFINE_bool(fix2,false,"Forcing the position of A2 fixed to be anchor");
DEFINE_bool(fix3,false,"Forcing the position of A3 fixed to be anchor");






DEFINE_double(A00,0.0,"");
DEFINE_double(A01,0.0,"");
DEFINE_double(A02,0.0,"");
DEFINE_double(A03,0.0,"");

  
DEFINE_double(A10,0.0,"");
DEFINE_double(A11,0.0,"");
DEFINE_double(A12,0.0,"");
DEFINE_double(A13,0.0,"");

  
DEFINE_double(A20,0.0,"");
DEFINE_double(A21,0.0,"");
DEFINE_double(A22,0.0,"");
DEFINE_double(A23,0.0,"");

  
DEFINE_double(A30,0.0,"");
DEFINE_double(A31,0.0,"");
DEFINE_double(A32,0.0,"");
DEFINE_double(A33,0.0,"");

  


DEFINE_double(R00,0.0,"");
DEFINE_double(R01,0.0,"");
DEFINE_double(R02,0.0,"");
DEFINE_double(R03,0.0,"");

  
DEFINE_double(R10,0.0,"");
DEFINE_double(R11,0.0,"");
DEFINE_double(R12,0.0,"");
DEFINE_double(R13,0.0,"");

  
DEFINE_double(R20,0.0,"");
DEFINE_double(R21,0.0,"");
DEFINE_double(R22,0.0,"");
DEFINE_double(R23,0.0,"");

  
DEFINE_double(R30,0.0,"");
DEFINE_double(R31,0.0,"");
DEFINE_double(R32,0.0,"");
DEFINE_double(R33,0.0,"");

  


////細胞の存在率
DEFINE_double(P0,0.0,"");
DEFINE_double(P1,0.0,"");
DEFINE_double(P2,0.0,"");
DEFINE_double(P3,0.0,"");

DEFINE_int32(celltype,4,"Number of cell type");

DEFINE_double(ir,1.0,"Initial size of the retinal radius");
DEFINE_double(gr,1.0,"Growth rate of the retinal radius");
DEFINE_double(gb,1.0,"Growth band: Width of Celll Growth Region");

DEFINE_bool(display,false,"Display in window");

DEFINE_double(r0,0.5,"intrinsic radius of cell type 0");
DEFINE_double(r1,0.5,"intrinsic radius of cell type 1");
DEFINE_double(r2,0.5,"intrinsic radius of cell type 2");
DEFINE_double(r3,0.5,"intrinsic radius of cell type 3");


