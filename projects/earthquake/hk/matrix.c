#include <math.h>
#include "knt.h"

/******************define some operations of matrix****************/

matrix plus(matrix a, matrix b) {
  a.pp = cplus(a.pp, b.pp);
  a.ps = cplus(a.ps, b.ps);
  a.sp = cplus(a.sp, b.sp);
  a.ss = cplus(a.ss, b.ss);
  a.sh = cplus(a.sh, b.sh);
  return(a);
}

matrix mltp(matrix a, matrix b) {
  matrix c;
  c.pp = cplus(cmltp(a.pp, b.pp), cmltp(a.ps, b.sp));
  c.ps = cplus(cmltp(a.pp, b.ps), cmltp(a.ps, b.ss));
  c.sp = cplus(cmltp(a.sp, b.pp), cmltp(a.ss, b.sp));
  c.ss = cplus(cmltp(a.sp, b.ps), cmltp(a.ss, b.ss));
  c.sh = cmltp(a.sh, b.sh);
  return(c);
}

matrix ngtv(matrix a) {
  a.pp = cngtv(a.pp);
  a.ps = cngtv(a.ps);
  a.sp = cngtv(a.sp);
  a.ss = cngtv(a.ss);
  a.sh = cngtv(a.sh);
  return(a);
}

matrix invs(matrix a) {
  complex det;
  matrix  c;
  det = cplus(cmltp(a.pp, a.ss), cngtv(cmltp(a.ps, a.sp)));
  det = cinvs(det);
  c.pp = cmltp(a.ss, det);
  c.ss = cmltp(a.pp, det);
  c.ps = cngtv(cmltp(a.ps, det));
  c.sp = cngtv(cmltp(a.sp, det));
  c.sh = cinvs(a.sh);
  return(c);
}

matrix trns(matrix a) {
  complex t;
  t = a.sp;
  a.sp = a.ps;
  a.ps = t;
  return(a);
}

matrix imlt(matrix a) {
  a.pp = cmplx(-a.pp.y, a.pp.x);
  a.ps = cmplx(-a.ps.y, a.ps.x);
  a.sp = cmplx(-a.sp.y, a.sp.x);
  a.ss = cmplx(-a.ss.y, a.ss.x);
  a.sh = cmplx(-a.sh.y, a.sh.x);
  return(a);
}
