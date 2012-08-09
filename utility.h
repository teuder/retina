#ifndef UTILITY_H
#define UTILITY_H

#include<iostream>
#include<fstream>
#include <sstream> 
#include<cmath>
#include <iterator>
#include <vector>
#include <list>
#include <set>
#include <map>
#include <numeric>
#include <algorithm>

//#include <boost/lambda.hpp>
#include <boost/timer.hpp>

#include <google/gflags.h>



//using namespace boost::lambda;
using namespace std;

//*Command line flags
DECLARE_bool(debug);
DECLARE_bool(v);


//* Global variables
namespace utility{
  
  //円周率
  const double PI = 3.14159265359;

}



//*Typedefs
typedef unsigned int    uint;
typedef unsigned long   ulong;
typedef vector<int>     vector_int;
typedef vector<unsigned int>  vector_uint;
typedef vector<long>          vector_long;
typedef vector<unsigned long> vector_ulong;
typedef vector<double>        vector_double;
typedef vector< vector<double> >  vector_double_2D;
typedef vector< vector<int> >     vector_int_2D;
typedef vector< vector<uint> > 				vector_uint_2D;
typedef list<int> 							list_int;
typedef list<double> 						list_double;
typedef vector<list <double> > 				vector_list_double;
typedef vector<bool>  						vector_bool;
typedef map<int,double>						map_int_double;
typedef map<uint,double>						map_uint_double;
typedef map<int, vector<double> >			map_int_vector_double;
typedef map<uint, vector<double> >			map_uint_vector_double;
typedef map<string,double> 					map_string_double;
typedef set<int>							set_int;
typedef set<uint>							set_uint;
typedef set<double>							set_double;
typedef pair<unsigned long,double>			pair_ulong_double;
typedef pair<unsigned long,int>				pair_ulong_int;


//★関数の最初と最後に「BGN;」「END;」と記述することで、デバッグ用の文章を挿入する。
#ifdef __GNUC__
#ifdef DEBUG
#define BGN if(FLAGS_debug) std::cout << __PRETTY_FUNCTION__ << " BGN" << std::endl
#define END if(FLAGS_debug) std::cout << __PRETTY_FUNCTION__ << " END" << std::endl
#define HERE std::cout << __FILE__ << " " << __LINE__ << std::endl

#else
#define BGN 
#define END  
#define HERE
#endif /* DEBUG */  
#endif /* __GNUC__ */



//*Geometory

//★2点a,b間のユークリッド距離を返す
//double 	CalcDistance(vector_double a, vector_double b);

template<typename T>
double CalcDistance(T& a, T& b){
  
  if(a.size()!=b.size()){
    cerr << "Error: in fuction CalcDistance(): unequal vector size\n";
    abort();
  }
  
  double sum =0.0;
  typename T::iterator i=a.begin();
  typename T::iterator j=b.begin();
  while(i!=a.end()){
    sum+=(*i - *j)*(*i - *j);
    ++i;
    ++j;
  }
  
  return(std::sqrt(sum));
}

template<typename T>
double CalcSquareDistance(T& a, T& b){
  
  if(a.size()!=b.size()){
    cerr << "Error: in fuction CalcDistance(): unequal vector size\n";
    abort();
  }
  
  double sum =0.0;
  typename T::iterator i=a.begin();
  typename T::iterator j=b.begin();
  while(i!=a.end()){
    sum+=(*i - *j)*(*i - *j);
    ++i;
    ++j;
  }
  
  return(sum);
}

//★Mapで与えられた2点a,b間のユークリッド距離を返す
template<typename KEY_T,typename VAL_T>
double CalcDistanceMap(const map<KEY_T,VAL_T>& A,const map<KEY_T,VAL_T>& B){
  
  double sum(0.0);
  
  typename map<KEY_T,VAL_T>::const_iterator a=A.begin();
  while(a!=A.end()){
    typename map<KEY_T,VAL_T>::const_iterator b=B.find(a->first);
    if(b!=B.end()){
      sum+=(a->second - b->second)*(a->second - b->second);
    }
    else {
      sum+=(a->second - 0.0)*(a->second - 0.0);
    }
    ++a;
  }
  return(sqrt(sum));
}

template<typename KEY_T,typename VAL_T>
double CalcSquareDistanceMap(const map<KEY_T,VAL_T>& A,const map<KEY_T,VAL_T>& B){
  
  double sum(0.0);
  
  typename map<KEY_T,VAL_T>::const_iterator a=A.begin();
  while(a!=A.end()){
    typename map<KEY_T,VAL_T>::const_iterator b=B.find(a->first);
    if(b!=B.end()){
      sum+=(a->second - b->second)*(a->second - b->second);
    }
    else {
      sum+=(a->second - 0.0)*(a->second - 0.0);
    }
    ++a;
  }
  return(sqrt(sum));
}

//★与えられた点（X,Y）がX軸となす角を返す
template <typename T>
double CalcRadian(T x, T y){
  double theta;
  double r=1.0/sqrt(x*x+y*y);
  double X=x*r;
  double Y=y*r;
  //HERE;
  
  
  if(X>=0.0 && Y>=0.0){
    theta=acos(X);
  }else if(X<0.0 && Y>=0.0){
    theta=acos(X);
  }else if(X<0.0 && Y<0.0){
    theta=utility::PI-asin(Y);
  }else{
    theta=(2.0*utility::PI)-acos(X);
  }
  
  //cout << "theta="<< theta << " x=" << x << " y=" << y <<  " X=" << X << " Y=" << Y  << " r="<< r << endl;

  return(theta);
}


//★コンテナの要素を表示する関数テンプレート

template <typename T>
void Show(const T& container, string sep=" "){
  std::copy( container.begin(), container.end(), std::ostream_iterator< int >( std::cout, sep.c_str()) );
  cout << endl;
}

template <typename T>
void ShowSeq(const T& container)
{
  for(typename T::const_iterator i=container.begin(); i!=container.end(); ++i){
    cout << *i << "\t";
  }
  cout << endl;
}


template <typename T>
void ShowSeq2D(const T& container)
{
  for(typename T::const_iterator i=container.begin(); i!=container.end(); ++i){
    ShowSeq(*i);
  }
  cout << endl;
}



template <typename KEY,typename VAL>
void ShowMap(const map<KEY,VAL>& container)
{
  if(container.empty()){
    cout << "empty" << endl;
  }
  else{
    cout << "KEY\tVAL" << endl;
    for(typename map<KEY,VAL>::const_iterator i=container.begin(); i!=container.end(); ++i){
      cout << i->first << "\t" << i->second << endl;
    }
  }
}

template <typename KEY,typename VAL>
void ShowMap(const multimap<KEY,VAL>& container)
{
  if(container.empty()){
    cout << "empty" << endl;
  }
  else{
    cout << "KEY\tVAL" << endl;
    for(typename multimap<KEY,VAL>::const_iterator i=container.begin(); i!=container.end(); ++i){
      cout << i->first << "\t" << i->second << endl;
    }
  }
}


template <typename T>
void ShowMapValSeq(const T& container)
{
  cout << "KEY\tVAL\n";
  for(typename T::const_iterator i=container.begin(); i!=container.end(); ++i){
    cout << i->first << "\t";
    ShowSeq(i->second);
  }
}


//★vector< map<int,double> > の表示
template<typename KEY,typename VAL>
void ShowPhenotype(vector< map<KEY,VAL> > P)
{
  BGN;
  int i(0);
  for(typename vector< map<KEY,VAL> >::iterator p=P.begin(),pe=P.end();p!=pe;++p,++i){
    cout << "[" << i << "]\n";
    ShowMap(*p);
    cout << endl;
  }
  END;
}




//★Mapのキーを取得する
template <typename T> 
vector<typename T::key_type> GetMapKeyVector(const T& container)
{
  vector<typename T::key_type> result;
  result.reserve(container.size());
  for(typename T::const_iterator i=container.begin(); i!=container.end(); ++i){
    result.push_back(i->first);
  }
  return result;
}

template <typename T> 
void GetMapKeyVector(const T& container, vector<typename T::key_type>& result)
{
  result.reserve(container.size());
  for(typename T::const_iterator i=container.begin(); i!=container.end(); ++i){
    result.push_back(i->first);
  }
}

template <typename T>
set<typename T::key_type> GetMapKeySet(const T& container)
{
  set<typename T::key_type> result;
  for(typename T::const_iterator i=container.begin(); i!=container.end(); ++i){
    result.insert(i->first);
  }
  return result;
}

template <typename T>
void GetMapKeySet(const T& container, set<typename T::key_type>& result)
{
  for(typename T::const_iterator i=container.begin(); i!=container.end(); ++i){
    result.insert(i->first);
  }
}

//★Mapの値を取得する
template <typename T> inline
vector<typename T::mapped_type> GetMapValVector(const T& container)
{
  vector<typename T::mapped_type> result;
  result.reserve(container.size());
  for(typename T::const_iterator i=container.begin(); i!=container.end(); ++i){
    result.push_back(i->second);
  }
  return result;
}





//★マップの値同士の足し算
template<typename T1, typename T2>
map<T1,T2> operator+(map<T1,T2> A,const map<T1,T2> &B){
  typename map<T1,T2>::const_iterator b;
  for(b=B.begin();b!=B.end();++b){
    A[b->first] += b->second;
  }
  return A;
}


template<typename T1, typename T2>
map<T1,T2> operator-(map<T1,T2> A,const map<T1,T2> &B){
  typename map<T1,T2>::const_iterator b;
  for(b=B.begin();b!=B.end();++b){
    A[b->first] -= b->second;
  }
  return A;
}



//★ベクタの要素同士の加算と減算
template<typename T1>
vector<T1> operator+(vector<T1> A,const vector<T1>& B){
  typename vector<T1>::iterator a,ae;
  typename vector<T1>::iterator b;
  for(a=A.begin(),b=B.begin(),ae=A.end();a!=ae;++a,++b){
    *a += *b; 
  }
  return(A);
}

template<typename T1,typename T2>
vector<T1> operator+(vector<T1> A,const vector<T2>& B){
  typename vector<T1>::iterator a,ae;
  typename vector<T1>::iterator b;
  for(a=A.begin(),b=B.begin(),ae=A.end();a!=ae;++a,++b){
    *a += *b; 
  }
  return(A);
}


template<typename T1>
vector<T1> operator-(vector<T1> A,const vector<T1>& B){
  typename vector<T1>::iterator a,ae;
  typename vector<T1>::iterator b;
  for(a=A.begin(),b=B.begin(),ae=A.end();a!=ae;++a,++b){
    *a -= *b; 
  }
  return(A);
}

template<typename T1,typename T2>
vector<T1> operator-(vector<T1> A,const vector<T2>& B){
  typename vector<T1>::iterator a,ae;
  typename vector<T1>::iterator b;
  for(a=A.begin(),b=B.begin(),ae=A.end();a!=ae;++a,++b){
    *a -= *b; 
  }
  return(A);
}





//★ベクタのスカラー倍
template<typename S,typename T>
vector<T> operator*(const vector<T>& A,const S B){
  vector<T> res(A);
  typename vector<T>::iterator a;
  for(a=res.begin();a!=res.end();++a){
    (*a)*=B;
  }
  return(res);
}

template<typename S,typename T>
vector<T> operator*(const S B,const vector<T>& A){
  vector<T> res(A);
  typename vector<T>::iterator a;
  for(a=res.begin();a!=res.end();++a){
    (*a)*=B;
  }
  return(res);
}

//★マップのスカラー倍
template<typename S,typename T1,typename T2>
map<T1,T2> operator*(map<T1,T2> A,const S B){
  typename map<T1,T2>::iterator a;
  for(a=A.begin();a!=A.end();++a){
    a->second = B*a->second;
  }
  return(A);
}

template<typename S,typename T1>
map<T1,T1> operator*(map<T1,T1> A, const S B){
  typename map<T1,T1>::iterator a;
  for(a=A.begin();a!=A.end();++a){
    a.second *= B;
  }
  return(A);
}


template<typename S,typename T1,typename T2>
map<T1,T2> operator*(const S B, map<T1,T2> A){
  typename map<T1,T2>::iterator a;
  for(a=A.begin();a!=A.end();++a){
    a->second = B*a->second;
  }
  return(A);
}

template<typename S,typename T1>
map<T1,T1> operator*(const S B, map<T1,T1> A){
  typename map<T1,T1>::iterator a;
  for(a=A.begin();a!=A.end();++a){
    a.second *= B;
  }
  return(A);
}


//★マップの値とベクタの要素と足し算
template<typename KEY, typename VAL>
map<KEY,VAL> operator+(map<KEY,VAL> M, const vector<VAL>& V){
  typename map<KEY,VAL>::iterator m(M.begin());
  typename map<KEY,VAL>::iterator me(M.end());
  typename vector<VAL>::const_iterator v(V.begin());
  while(m!=me){
    (m++)->second += *(v++);
  }
  return(M);
}

template<typename KEY, typename VAL>
map<KEY,VAL> operator+(const vector<VAL>& V, map<KEY,VAL> M){
  typename map<KEY,VAL>::const_iterator m(M.begin());
  typename map<KEY,VAL>::const_iterator me(M.end());
  typename vector<VAL>::const_iterator v(V.begin());
  while(m!=me){
    (m++)->second += *(v++);
  }
  return(M);
}





//★Convert number to string
template<typename T> inline 
string itos ( T number )
{
	ostringstream oss;
	oss<< number;
	return oss.str();
}

template<typename T> inline
string tos ( T number )
{
	ostringstream oss;
	oss << number;
	return oss.str();
}

template<typename T> inline
string tos (T number, int width, char fill)
{
  ostringstream oss;
  if(width) oss.width(width);
  oss.fill(fill);
  oss.precision(1);
  oss<< fixed << number;
  return oss.str();
}







//★Check existance of the specified file or directry
inline
bool FileExist(string filename){
	fstream fst(filename.c_str(), ios::in); 
	if(!fst) return(false);
	else return(true);
}



//★シーケンシャル・コンテナの要素に対する処理

template <typename Container_t>
const typename Container_t::value_type max ( const Container_t& A) {
  typename Container_t::value_type res= *(A.begin());
  
  for(typename Container_t::const_iterator a=A.begin();a!=A.end();++a){
    res=max(res,*a);
  }
  return res;     // or: return comp(a,b)?b:a; for the comp version
}

template <typename Container_t>
const typename Container_t::value_type min ( const Container_t& A) {
  typename Container_t::value_type res= *(A.begin());
  for(typename Container_t::const_iterator a=A.begin();a!=A.end();++a){
    res=min(res,*a);
  }
  return res;     // or: return comp(a,b)?b:a; for the comp version
}

template <typename Container_t>
const Container_t range ( const  Container_t& A) {

  vector<typename Container_t::value_type> res(2);
  res[0] = res[1] = *(A.begin());
  
  for(typename Container_t::const_iterator a=A.begin();a!=A.end();++a){
    res[0]=min(res[0],*a);
    res[1]=min(res[1],*a);
  }
  
  return (Container_t(res.begin(),res.end()));
}


template<typename T> inline
double sum(T& a){
  return ((double) accumulate(a.begin(), a.end(), 0.0));
}

template<typename T> inline
double mean(T& a){
  BGN;
  double sum = 0.0;
  int c=0;
  typename T::iterator i=a.begin();
  while(isfinite(*i) && i!=a.end()){
    sum+=(*i)*(*i);
    ++i;
    ++c;
  }
  
  

  return(sum/c);
}

template<typename T> inline
double sumsq(T& a){
  double sum = 0.0;
  typename T::iterator i=a.begin();
  while(i!=a.end()){
    sum+=(*i)*(*i);
    ++i;
  }
  return sum;
}

template<typename T> inline
double variance(T& a){
  double meana=mean(a);
  
  double sum = 0.0;
  typename T::iterator i=a.begin();
  while(i!=a.end()){
    sum+=(*i - meana)*(*i - meana);
    ++i;
  }
  return sum/(a.size());
}

template<typename T> inline
double stddev(T& a){
  return(sqrt(variance(a)));
}

//中央値
template<typename T> inline
double median(T& a){
  BGN;
  //if(a.empty()) return(0/0);
  if(a.empty()) return 0;
  HERE;
  vector<double> res(a.begin(),a.end());
  HERE;
  sort(res.begin(),res.end());
  HERE;
  if(res.size()%2==1){
    HERE;
    END;
    return(res[res.size()/2]);
  }
  else{
    HERE;
    size_t i=(res.size()+1)/2;
    HERE;
    END;
    return((res[i]+res[i-1])*0.5);
  }
  
}


//★シーケンスに0,1,2,...を格納する
template<typename T> inline
void sequence(T& seq){
  typename T::iterator s=seq.begin();
  typename T::iterator se=seq.end();
  int i=0;
  while(s!=se){
    *(s++)=(i++);
  }
}

template<typename T> inline
void sequence(T s,T se){
  int i=0;
  while(s!=se){
    *(s++)=(i++);
  }
}




//★[0,N] の整数列を返す
template<typename T> 
inline vector<int> seqZN(int N) {
  
  vector<int> S(N+1);
  auto s=S.begin();
  int i=0;
  while(i!=N+1){
    *(s++)=(i++);
  }
  return(S);
}

//★[0,N-1] の整数列を返す
template<typename T> 
inline vector<int> seqZN_1(int N) {
  
  vector<int> S(N);
  auto s=S.begin();
  int i=0;
  while(i!=N){
    *(s++)=(i++);
  }
  return(S);
}


inline vector<int> seq(int N) {
  
  vector<int> S(N);
  auto s=S.begin();
  int i=0;
  while(i!=N){
    *(s++)=(i++);
  }
  return(S);
}




//★イテレータで渡された、範囲のシーケンスコンテナに対する処理

//総和
template<typename T> inline
double sum(T begin, T end){
  double sum(0.0);
  while(begin!=end){
    sum += *(begin++);
  }
  return sum;
}

//総和平均
template<typename T> inline
double mean(T begin, T end){
  BGN;
  size_t size(distance(begin,end));
  double sum(0.0);
  while(begin!=end){
    sum += *(begin++);
  }
  END;
  return(sum/size);
}

//平方和
template<typename T> inline
double sumsq(T begin, T end){
  double sum(0.0);
  while(begin!=end){
      sum += (*begin)*(*begin);
      ++begin;
  }
  return sum;
}

//分散
template<typename T> inline
double variance(T begin, T end){
  size_t size(distance(begin,end));
  double mean=mean(begin,end);
  double sum=0.0;
  while(begin!=end){
    sum += (*begin - mean)*(*begin - mean);
    ++begin;
  }
  return(sum/size);
}

//標準偏差
template<typename T> inline
double stddev(T begin, T end){
  return(sqrt(variance(begin,end)));
}

//総乗
template<typename T> inline
double product(T begin, T end){
  double sum(1.0);
  while(begin!=end){
    sum *= *(begin++);
  }
  return sum;
}


//相乗平均（幾何平均）
template<typename T> inline
double geomean(T begin, T end){
  double sum(1.0);
  size_t size(distance(begin,end));
  while(begin!=end){
    sum *= *(begin++);
  }
  return pow(sum,1.0/size);
}

//中央値
template<typename T> inline
double median(T begin, T end)
{
  vector<double> res(begin,end);
  sort(res.begin(),res.end());
  if(res.size()%2==1){
    return(res[res.size()/2]);
  }
  else{
    size_t i=(res.size()+1)/2;
    return((res[i]+res[i-1])*0.5);
  }
}





//★n次元の点を回転移動させた点に変換する
template<typename T> inline
void rotation_transform(T& point,const T& theta){
  typename T::iterator x=point.begin();
  typename T::iterator y=point.begin();++y;
  typename T::const_iterator t=theta.begin();
  while(y!=point.end()){
    double x_old=*x;
    double y_old=*y;
    *x = x_old * cos(*t) - y_old * sin(*t);
    *y = x_old * sin(*t) + y_old * cos(*t);
    ++x;++y;++t;
  }
}


//★vectorをクリアしてメモリを解放する
template<typename T>
void MemClear(vector<T>& vec){
  vec.clear();
  vector<T>(vec).swap(vec);
}



// Below 2 lines return same results.
// for_each(vec.begin(), vec.end(), function);
// for_each(vec, function);
template <typename T_container, typename T_function>
T_function for_each(T_container& rcontainer, T_function function) {
  return for_each(rcontainer.begin(), rcontainer.end(), function);
}



#endif


