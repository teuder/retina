#include "random.h"


using namespace std;



////uint_32 genrand()の最大値(2^32-1)
//const double MAX = 4294967295.0;






//標準正規乱数
double NormRand()
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
double NormRand(double mu, double sigma)
{
	return(sigma*NormRand()+mu);
}


//平均値lambdaのポワソン乱数（整数）を返す関数
int p_rand(double lambda)
{
	lambda = exp(lambda) * ZP1rand();
	int k = 0;
	while (lambda > 1) {
		lambda *= ZP1rand();
		++k;
	}
	return k;
}




//-1以上1以下の乱数を生成する
double MP1rand(void)
{
	return (2*(genrand()*(1.0/myrand::MAX))-1.0);
}





//任意の確率[0,1]でtrueを返す関数
//bool Prob1(double x)
//{
  //if(Zp1rand()<x) return true;
  //else return false;
  ////if(x > genrand()*(1.0/(MAX+1))) return true;
//}

//bool Prob0(double x)
//{
  //if(Zp1rand()<x) return false;
  //else return true;
//}


//0以上n未満の整数乱数生成
int Znirand(int n)
{
	return (int(genrand()%(n)));
}




////0以上n以下の実数乱数生成
//double ZNdrand(int n)
//{
	//return((genrand()*(1.0/MAX))*double(n));
//}

//double ZNdrand(double n)
//{
	//return((genrand()*(1.0/MAX))*n);
//}

////0より大、n以下の実数乱数
//double zNdrand(int n)
//{
	//return ((double(genrand()) + 0.5)*(1.0/4294967296.5)*double(n)); 
//}

//double zNdrand(double n)
//{
	//return ((double(genrand()) + 0.5)*(1.0/4294967296.5)*n);
//}


////★0より大、nより小の実数乱数
//double zndrand(double n)
//{
  
  
//}


////2分の1の確率でtrueまたはfalseを返す関数
//bool TFrand()
//{
	//if(0.5 > genrand()*(1.0/MAX)) return true;
	//else return false;
//}

////2分の1の確率で+1または-1を返す関数
//extern double P1orM1rand(){
	//if(0.5 > genrand()*(1.0/MAX)) return 1.0;
	//else return -1.0;	
//}


////-n から n の間の整数の乱数
//int iMPNrand(int n)
//{
	//int a;
	
	//a = (genrand()%(2*n+1))-n;
	
	//return a;
	
//}

////-n から n の間の乱数
//double dMPNrand(int n)
//{
	//double RAND;
	
	//RAND = genrand()/MAX;
	
	//return(n*(2*RAND-1));
//}

////-n から n の間の乱数
//double dMPNrand(double n)
//{
	//double RAND;
	
	//RAND = genrand()/MAX;
	
	//return(n*(2*RAND-1));
	
	
//}








