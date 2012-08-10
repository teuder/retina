#ifndef RANDOM_H
#define RANDOM_H

#include <ctime>
#include <cmath>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <map>

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


//-1以上1以下の乱数を生成する
inline double MP1rand(void)
{
	return (2*(genrand()*(1.0/myrand::MAX))-1.0);
}



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
double NMrand(T n,T m)
{
  return((m-n)*ZP1rand()+n);
}

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
//標準正規乱数
inline double NormRand()
{
	double t=sqrt(-2.0 * log(1-Zp1rand()));
	double u=myrand::PI*2.0*Zp1rand();
	return(t*cos(u));
}


//一般的な平均値mu、標準偏差sigmaの正規分布をする乱数
//これは標準正規乱数を利用して作る
//一般正規分布f(x)があったとき
//xを変数変換し
//z = (x - mu)/sigma
//として得られる分布g(z)は標準正規分布となる
//つまり、これを逆変換すると
//x = sigma*x + mu
//であり
//標準正規分布を一般正規分布に変換できる
inline double NormRand(double mu, double sigma)
{
	return(sigma*NormRand()+mu);
}

//平均値mu、標準偏差sigmaの正規分布
//template<typename T1,typename T2>
//double NormRand(T1 mu, T2 sigma)
//{
  //return(sigma*NormRand()+mu);
//}




//平均値lambdaのポワソン乱数（整数）を返す関数
//平均値lambdaのポワソン乱数（整数）を返す関数
inline int p_rand(double lambda)
{
	lambda = exp(lambda) * ZP1rand();
	int k = 0;
	while (lambda > 1) {
		lambda *= ZP1rand();
		++k;
	}
	return k;
}



//0以上n以下の整数乱数生成
inline unsigned int ZNrand(unsigned int n)
{
  return (genrand()%(n+1));
};

inline double ZNrand(double n)
{
  return (n*ZP1rand());
};

inline double Znrand(double n)
{
  return (n*Zp1rand());
};

//0以上n未満の整数乱数生成
inline int Znirand(int n)
{
	return (int(genrand()%(n)));
}

inline int Znrand(int n)
{
	return (int(genrand()%(n)));
}


inline std::ptrdiff_t myrandom (std::ptrdiff_t i) { return genrand()%i;}



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




//与えられたシーケンスの比率（確率）でtrueとなる要素の番号を返す
template<typename Seq>
size_t Prob(Seq seq)
{
  
  //ルーレットを作る
  double sum=0.0;
  size_t i=0;
  multimap<double,size_t> border_value;
  for(auto s : seq) {
    sum+=s;
    border_value.insert(make_pair(sum, i++));
  }
  
  //ダーツを投げる
  const double dart=ZNrand(sum);
  

  //あたった場所を調べる
  const auto b=border_value.lower_bound(dart);
  return(b->second);
  
}



#endif
