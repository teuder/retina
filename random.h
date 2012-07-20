#ifndef RANDOM_H
#define RANDOM_H

#include <ctime>
#include <cmath>
#include <unistd.h>
#include <iostream>
#include <fstream>


#include "./SFMT/SFMT.h"





using namespace std;


namespace myrand {
  
  //円周率
  const double PI     = 3.14159265359;
  //extern const double PI;
  //uint_32 genrand()の最大値(2^32-1)
  const double MAX = 4294967295.0;
}



/*#########いろんな乱数とか確率の関数の宣言##########*/

//
inline unsigned int dev_urandom(){
  unsigned int x;
  try{
    ifstream fin("/dev/urandom", ios::binary | ios::in);
    fin.exceptions(ios::failbit | ios::badbit);
    fin.read( reinterpret_cast<char*>(&x), sizeof(x) );
  }
  catch(ios::failure& e){ throw ios::failure("/dev/urandom"); }
  return x;
}

//乱数の初期化
inline unsigned int InitRand(unsigned int seed=0){	
  if(seed==0) seed=dev_urandom();
  init_gen_rand(seed);
  return seed;
}


//基本となる整数乱数[0,MAX]
inline static unsigned int genrand(void)
{
  return gen_rand32();
}




/* These real versions are due to Isaku Wada */

//[0,1]-real-interval
inline static double ZP1rand(void)
{
  return to_real1(gen_rand32());
}

//[0,1)-real-interval
inline static double Zp1rand(void)
{
  return to_real2(gen_rand32());
}


//(0,1)-real-interval
inline static double zp1rand(void)
{
  return to_real3(gen_rand32());
}

inline static double to_real4(uint32_t v)
{
  return (v+0.5) * (1.0/4294967295.5);  
}

//(0,1]の実数乱数
inline static double zP1rand(void)
{
  return to_real4(gen_rand32());
}


////-1以上1以下の乱数を生成する
//extern double MP1rand(void);



////任意の確率でtrueを返す関数
extern bool Prob1(double x);
extern bool Prob0(double x);


////Zは0を含む
////zは0を含まない

////0以上n以下の整数乱数生成
//extern int Znirand(int n);


////0以上n以下の整数乱数生成
//extern int ZNirand(int n);

////0以上n以下の実数乱数生成
//extern double ZNdrand(int n);
//extern double ZNdrand(double n);

////0より大、n以下の実数乱数
//extern double zNdrand(int n);
//extern double zNdrand(double n);

//2分の1の確率でtrueまたはfalseを返す関数
extern bool TFrand();

//2分の1の確率で+1または-1を返す関数
extern double P1orM1rand();


//-n から n の間の整数の乱数 
int iMPNrand(int n);

//-n から n の間の乱数[-n,n], -nとnを含む
template<typename T>
double dMPNrand(T n)
{
  return(n*(2*(genrand()/myrand::MAX)-1));
}


//[n,m]の実数乱数
template<typename T>
double dNMrand(T n,T m)
{
  return((m-n)*ZP1rand()+n);
}


//[n,m]の整数乱数
inline int iNMrand(int n,int m)
{
  return(genrand()%(m-n)+n);
}





//標準正規乱数
extern double NormRand();

//平均値mu、標準偏差sigmaの正規分布
template<typename T1,typename T2>
double NormRand(T1 mu, T2 sigma)
{
  return(sigma*NormRand()+mu);
}


//平均値lambdaのポワソン乱数（整数）を返す関数
extern int p_rand(double lambda);



//0以上n以下の整数乱数生成
inline int ZNirand(int n)
{
  return (int(genrand()%(n+1)));
};


//任意の確率[0,1]でtrueを返す関数
inline bool Prob(double x)
{
  if(Zp1rand()<x) return true;
  else return false;
}

////任意の確率でtrueを返す関数
inline bool Prob1(double x)
{
  if(Zp1rand()<x) return true;
  else return false;
}
////任意の確率でfalseを返す関数
inline bool Prob0(double x)
{
  if(Zp1rand()<x) return false;
  else return true;
}





#endif
