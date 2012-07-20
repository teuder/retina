#include "main.h"
#include "simulation.h"
#include "utility.h"







int main(int argc, char** argv)
{
  
  //コマンドライン引数の処理
  google::SetUsageMessage("Simulating Multicellular Morphogenesis");
  google::ParseCommandLineFlags(&argc, &argv, false);
  
      

    
  
  

  //boost::timer t;
  Simulation sim;
  sim.Run();
  //if(FLAGS_A00==0.0) cout << FLAGS_A00 << endl;
  //cout << t.elapsed() << "s" << endl;
  //double theta=1;
  //cout << CalcRadian(cos(theta),sin(theta)) << endl;
  
  
  
  return 0; 
}














