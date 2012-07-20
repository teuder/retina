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
  
  //�~����
  const double PI     = 3.14159265359;
  //extern const double PI;
  //uint_32 genrand()�̍ő�l(2^32-1)
  const double MAX = 4294967295.0;
}



/*#########�����ȗ����Ƃ��m���̊֐��̐錾##########*/

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

//�����̏�����
inline unsigned int InitRand(unsigned int seed=0){	
  if(seed==0) seed=dev_urandom();
  init_gen_rand(seed);
  return seed;
}


//��{�ƂȂ鐮������[0,MAX]
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

//(0,1]�̎�������
inline static double zP1rand(void)
{
  return to_real4(gen_rand32());
}


////-1�ȏ�1�ȉ��̗����𐶐�����
//extern double MP1rand(void);



////�C�ӂ̊m����true��Ԃ��֐�
extern bool Prob1(double x);
extern bool Prob0(double x);


////Z��0���܂�
////z��0���܂܂Ȃ�

////0�ȏ�n�ȉ��̐�����������
//extern int Znirand(int n);


////0�ȏ�n�ȉ��̐�����������
//extern int ZNirand(int n);

////0�ȏ�n�ȉ��̎�����������
//extern double ZNdrand(int n);
//extern double ZNdrand(double n);

////0����An�ȉ��̎�������
//extern double zNdrand(int n);
//extern double zNdrand(double n);

//2����1�̊m����true�܂���false��Ԃ��֐�
extern bool TFrand();

//2����1�̊m����+1�܂���-1��Ԃ��֐�
extern double P1orM1rand();


//-n ���� n �̊Ԃ̐����̗��� 
int iMPNrand(int n);

//-n ���� n �̊Ԃ̗���[-n,n], -n��n���܂�
template<typename T>
double dMPNrand(T n)
{
  return(n*(2*(genrand()/myrand::MAX)-1));
}


//[n,m]�̎�������
template<typename T>
double dNMrand(T n,T m)
{
  return((m-n)*ZP1rand()+n);
}


//[n,m]�̐�������
inline int iNMrand(int n,int m)
{
  return(genrand()%(m-n)+n);
}





//�W�����K����
extern double NormRand();

//���ϒlmu�A�W���΍�sigma�̐��K���z
template<typename T1,typename T2>
double NormRand(T1 mu, T2 sigma)
{
  return(sigma*NormRand()+mu);
}


//���ϒllambda�̃|���\�������i�����j��Ԃ��֐�
extern int p_rand(double lambda);



//0�ȏ�n�ȉ��̐�����������
inline int ZNirand(int n)
{
  return (int(genrand()%(n+1)));
};


//�C�ӂ̊m��[0,1]��true��Ԃ��֐�
inline bool Prob(double x)
{
  if(Zp1rand()<x) return true;
  else return false;
}

////�C�ӂ̊m����true��Ԃ��֐�
inline bool Prob1(double x)
{
  if(Zp1rand()<x) return true;
  else return false;
}
////�C�ӂ̊m����false��Ԃ��֐�
inline bool Prob0(double x)
{
  if(Zp1rand()<x) return false;
  else return true;
}





#endif
