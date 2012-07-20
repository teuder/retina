
#ifndef PARAMETERH
#define PARAMETERH
#include <google/gflags.h>
#include <vector>
#include <string>


class Parameter
{
  public:
  std::vector<int> N;
  std::vector<double> U,L,Pd;
  std::vector<std::vector<double>> A,R;
  std::vector<bool> fix;
  std::vector<std::string> I;
  std::string OssFlags();
  void OutFlags(std::string filepath);
  
  Parameter();
  
};

DECLARE_double(Tm);
DECLARE_double(H);
DECLARE_double(D);
DECLARE_double(Fs);
DECLARE_double(Ts);
DECLARE_double(Pd);
DECLARE_double(Td);
DECLARE_int32(periodic);
DECLARE_double(Xp);
DECLARE_double(Yp);
DECLARE_double(Ta);
DECLARE_string(o);


DECLARE_int32(N0);
DECLARE_int32(N1);
DECLARE_int32(N2);
DECLARE_int32(N3);
DECLARE_int32(N4);
DECLARE_int32(N5);
DECLARE_int32(N6);

DECLARE_double(U);
DECLARE_double(L);

DECLARE_double(U0);
DECLARE_double(U1);
DECLARE_double(U2);
DECLARE_double(U3);
DECLARE_double(U4);
DECLARE_double(U5);
DECLARE_double(U6);


DECLARE_double(L0);
DECLARE_double(L1);
DECLARE_double(L2);
DECLARE_double(L3);
DECLARE_double(L4);
DECLARE_double(L5);
DECLARE_double(L6);

DECLARE_string(I);
DECLARE_string(I0);
DECLARE_string(I1);
DECLARE_string(I2);
DECLARE_string(I3);
DECLARE_string(I4);
DECLARE_string(I5);
DECLARE_string(I6);

DECLARE_double(Pd);
DECLARE_double(Pd0);
DECLARE_double(Pd1);
DECLARE_double(Pd2);
DECLARE_double(Pd3);
DECLARE_double(Pd4);
DECLARE_double(Pd5);
DECLARE_double(Pd6);


DECLARE_double(r);
DECLARE_bool(fill);
DECLARE_bool(edges);
DECLARE_bool(velocity);
DECLARE_string(mode);
DECLARE_bool(anime);
DECLARE_int32(repli);
DECLARE_double(Tnoise);
DECLARE_double(noise);
DECLARE_double(inispace);
DECLARE_bool(outvtx);
DECLARE_bool(circle);
DECLARE_bool(makeinit);



//DECLARE_double(Asame0);
//DECLARE_double(Asame1);
//DECLARE_double(Asame2);
//DECLARE_double(Asame3);

//DECLARE_double(Adiff0);
//DECLARE_double(Adiff1);
//DECLARE_double(Adiff2);
//DECLARE_double(Adiff3);


//DECLARE_double(Rsame0);
//DECLARE_double(Rsame1);
//DECLARE_double(Rsame2);
//DECLARE_double(Rsame3);

//DECLARE_double(Rdiff0);
//DECLARE_double(Rdiff1);
//DECLARE_double(Rdiff2);
//DECLARE_double(Rdiff3);




DECLARE_double(Tabl);
DECLARE_double(d1);
DECLARE_double(x1);
DECLARE_double(y1);

DECLARE_bool(fix0);
DECLARE_bool(fix1);
DECLARE_bool(fix2);
DECLARE_bool(fix3);
DECLARE_bool(fix4);
DECLARE_bool(fix5);
DECLARE_bool(fix6);

DECLARE_double(A00);
DECLARE_double(A01);
DECLARE_double(A02);
DECLARE_double(A03);
DECLARE_double(A04);
DECLARE_double(A05);
DECLARE_double(A06);
DECLARE_double(A10);
DECLARE_double(A11);
DECLARE_double(A12);
DECLARE_double(A13);
DECLARE_double(A14);
DECLARE_double(A15);
DECLARE_double(A16);
DECLARE_double(A20);
DECLARE_double(A21);
DECLARE_double(A22);
DECLARE_double(A23);
DECLARE_double(A24);
DECLARE_double(A25);
DECLARE_double(A26);
DECLARE_double(A30);
DECLARE_double(A31);
DECLARE_double(A32);
DECLARE_double(A33);
DECLARE_double(A34);
DECLARE_double(A35);
DECLARE_double(A36);
DECLARE_double(A40);
DECLARE_double(A41);
DECLARE_double(A42);
DECLARE_double(A43);
DECLARE_double(A44);
DECLARE_double(A45);
DECLARE_double(A46);
DECLARE_double(A50);
DECLARE_double(A51);
DECLARE_double(A52);
DECLARE_double(A53);
DECLARE_double(A54);
DECLARE_double(A55);
DECLARE_double(A56);
DECLARE_double(A60);
DECLARE_double(A61);
DECLARE_double(A62);
DECLARE_double(A63);
DECLARE_double(A64);
DECLARE_double(A65);
DECLARE_double(A66);
DECLARE_double(A00);
DECLARE_double(A01);
DECLARE_double(A02);
DECLARE_double(A03);
DECLARE_double(A04);
DECLARE_double(A05);
DECLARE_double(A06);



DECLARE_double(R00);
DECLARE_double(R01);
DECLARE_double(R02);
DECLARE_double(R03);
DECLARE_double(R04);
DECLARE_double(R05);
DECLARE_double(R06);
DECLARE_double(R10);
DECLARE_double(R11);
DECLARE_double(R12);
DECLARE_double(R13);
DECLARE_double(R14);
DECLARE_double(R15);
DECLARE_double(R16);
DECLARE_double(R20);
DECLARE_double(R21);
DECLARE_double(R22);
DECLARE_double(R23);
DECLARE_double(R24);
DECLARE_double(R25);
DECLARE_double(R26);
DECLARE_double(R30);
DECLARE_double(R31);
DECLARE_double(R32);
DECLARE_double(R33);
DECLARE_double(R34);
DECLARE_double(R35);
DECLARE_double(R36);
DECLARE_double(R40);
DECLARE_double(R41);
DECLARE_double(R42);
DECLARE_double(R43);
DECLARE_double(R44);
DECLARE_double(R45);
DECLARE_double(R46);
DECLARE_double(R50);
DECLARE_double(R51);
DECLARE_double(R52);
DECLARE_double(R53);
DECLARE_double(R54);
DECLARE_double(R55);
DECLARE_double(R56);
DECLARE_double(R60);
DECLARE_double(R61);
DECLARE_double(R62);
DECLARE_double(R63);
DECLARE_double(R64);
DECLARE_double(R65);
DECLARE_double(R66);



              
              
              
#endif
