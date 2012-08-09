#include "main.h"
#include "simulation.h"
#include "utility.h"
#include "random.h"








int main(int argc, char** argv)
{
  
  
  


  //init();
  //glutDisplayFunc(display);

  
  //コマンドライン引数の処理
  google::SetUsageMessage("Simulating Multicellular Morphogenesis");
  google::ParseCommandLineFlags(&argc, &argv, false);
  
  InitRand();
  
  //boost::timer t;
  Simulation sim;
  sim.RunGrowingRetina();
  //if(FLAGS_gr>0.0)    
  //else if(!FLAGS_makeinit)      sim.Run(argc, argv);         //普通のシミュレーション
  //else if(FLAGS_makeinit)  sim.RunMakeInitFile(); //細胞初期化ファイルを作る場合

  
  //if(FLAGS_A00==0.0) cout << FLAGS_A00 << endl;
  //cout << t.elapsed() << "s" << endl;
  //double theta=1;
  //cout << CalcRadian(cos(theta),sin(theta)) << endl;
  
  
  
  return 0; 
}














