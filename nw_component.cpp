#include <Rcpp.h>
using namespace Rcpp;

class Point {
public:
  Point( double x_, double y_) : x(x_), y(y_){}

  double x, y ;
} ;

double square( double x) {
  return x*x ;
}
double distance( const Point& p1, const Point& p2 ){
  return sqrt( square( p1.x - p2.x) + square( p1.y - p2.y ) ) ;
}

RCPP_MODULE(play){

  class_<Point>("Point")
  .constructor<double,double>()
  .field( "x", &Point::x)
  .field( "y", &Point::y);
};


/*** R
origin <- new( Point, 0, 0 )
# origin$square()
origin$distance(origin, origin)
*/
