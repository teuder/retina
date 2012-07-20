#ifndef MYPOINTC2_IOSTREAM_H
#define MYPOINTC2_IOSTREAM_H

#include "MyPointC2.h"
std::ostream &
 operator<<(std::ostream &os, const MyPointC2 &p);




 std::istream &
 operator>>(std::istream &is, MyPointC2 &p);
 
 
 

 
#endif //MYPOINTC2_IOSTREAM_H
