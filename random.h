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


//-1�ȏ�1�ȉ��̗����𐶐�����
inline double MP1rand(void)
{
	return (2*(genrand()*(1.0/myrand::MAX))-1.0);
}



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
double NMrand(T n,T m)
{
  return((m-n)*ZP1rand()+n);
}

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
//�W�����K����
inline double NormRand()
{
	double t=sqrt(-2.0 * log(1-Zp1rand()));
	double u=myrand::PI*2.0*Zp1rand();
	return(t*cos(u));
}


//��ʓI�ȕ��ϒlmu�A�W���΍�sigma�̐��K���z�����闐��
//����͕W�����K�����𗘗p���č��
//��ʐ��K���zf(x)���������Ƃ�
//x��ϐ��ϊ���
//z = (x - mu)/sigma
//�Ƃ��ē����镪�zg(z)�͕W�����K���z�ƂȂ�
//�܂�A������t�ϊ������
//x = sigma*x + mu
//�ł���
//�W�����K���z����ʐ��K���z�ɕϊ��ł���
inline double NormRand(double mu, double sigma)
{
	return(sigma*NormRand()+mu);
}

//���ϒlmu�A�W���΍�sigma�̐��K���z
//template<typename T1,typename T2>
//double NormRand(T1 mu, T2 sigma)
//{
  //return(sigma*NormRand()+mu);
//}




//���ϒllambda�̃|���\�������i�����j��Ԃ��֐�
//���ϒllambda�̃|���\�������i�����j��Ԃ��֐�
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



//0�ȏ�n�ȉ��̐�����������
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

//0�ȏ�n�����̐�����������
inline int Znirand(int n)
{
	return (int(genrand()%(n)));
}

inline int Znrand(int n)
{
	return (int(genrand()%(n)));
}


inline std::ptrdiff_t myrandom (std::ptrdiff_t i) { return genrand()%i;}



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




//�^����ꂽ�V�[�P���X�̔䗦�i�m���j��true�ƂȂ�v�f�̔ԍ���Ԃ�
template<typename Seq>
size_t Prob(Seq seq)
{
  
  //���[���b�g�����
  double sum=0.0;
  size_t i=0;
  multimap<double,size_t> border_value;
  for(auto s : seq) {
    sum+=s;
    border_value.insert(make_pair(sum, i++));
  }
  
  //�_�[�c�𓊂���
  const double dart=ZNrand(sum);
  

  //���������ꏊ�𒲂ׂ�
  const auto b=border_value.lower_bound(dart);
  return(b->second);
  
}



#endif
