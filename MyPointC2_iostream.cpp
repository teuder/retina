#include"MyPointC2_iostream.h"

std::ostream &
operator<<(std::ostream &os, const MyPointC2 &p)
{
  switch(os.iword(CGAL::IO::mode)) {
    case CGAL::IO::ASCII :
    return os << p.x() << ' ' << p.y() << ' ' << p.vx() << ' ' << p.vy() << ' ' << p.color() << ' ' << p.radius();
    case CGAL::IO::BINARY :
    CGAL::write(os, p.x());
    CGAL::write(os, p.y());
    CGAL::write(os, p.vx());
    CGAL::write(os, p.vy());
    CGAL::write(os, p.color());
    CGAL::write(os, p.radius());
    return os;
    default:
    return os << "MyPointC2(" << p.x() << ", " << p.y() << ") v(" << p.vx() << ", " << p.vy() << ") type(" << p.color() << ") r(" << p.radius() << ")" ;
  }
}



std::istream &
operator>>(std::istream &is, MyPointC2 &p)
{
  double x, y, vx, vy, r;
  int c;
  switch(is.iword(CGAL::IO::mode)) {
    case CGAL::IO::ASCII :
    is >> x >> y >> vx >> vy >> c  >> r;
    break;
    case CGAL::IO::BINARY :
    CGAL::read(is, x);
    CGAL::read(is, y);
    CGAL::read(is, vx);
    CGAL::read(is, vy);
    CGAL::read(is, c);
    CGAL::read(is, r);
    break;
    default:
    std::cerr << "" << std::endl;
    std::cerr << "Stream must be in ascii or binary mode" << std::endl;
    break;
  }
  if (is) {
    p = MyPointC2(x, y, c, r, vx,vy);
  }
  return is;
}



