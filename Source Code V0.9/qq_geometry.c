//Geometry methods

#ifndef PARILIB
#define PARILIB
#include <pari/pari.h>
#endif

#ifndef METHDECL
#define METHDECL
#include "qquadraticdecl.h"
#endif

//The length (lg, so technically length+1) of a circle/line and arc/segment, and a normalized boundary
#define CIRCLEN 4
#define ARCLEN 9
#define NORMBOUND 9

//STATIC DECLARATIONS
static GEN mobius_arcseg(GEN M, GEN c, int isarc, GEN tol, long prec);
static GEN mobius_circle(GEN M, GEN c, GEN tol, long prec);
static GEN mobius_line(GEN M, GEN l, GEN tol, long prec);
static GEN normalizedbasis_shiftpoint(GEN c, GEN r, int initial, long prec);

//Circle is stored as [centre, radius, 0], where 0 means a bona fide cicle (and not a line)
//Circle arc is [centre, radius, start pt, end pt, start angle, end angle, dir, 0]. It is the arc counterclockwise from startpt to endpt, and dir=1 means oriented counterclockwise, and dir=-1 means oriented clockwise. This can also be left uninitialized if arc is not oriented. The final 0 represents that we have an arc and not a segment.
//Line is stored as [slope, intercept, 1], where the line is y=slope*x+intercept unless slope=oo, where it is x=intercept instead, and the 1 just means a bona fide line (not circle).
//Line segment is stored as [slope, intercept, startpt, endpt, 0, ooendptor, dir, 1] (the extra 0 is to format it the same as a circle arc). The final 1 is to signal a segment. dir=1 means we go from startpt to endpt in the upper half plane, and -1 means we go through infinity (only when neither endpoint is infinite). If one endpoint is oo, ooendptor stores which way we get to it. If ooendptor=1, this means the segment travels either vertically up or right, and -1 means the arc is vertically down or left.

//GEN tol -> The tolerance.


//BASIC LINE, CIRCLE, AND POINT OPERATIONS

//Creates the arc on circle c going from p1 to p2 counterclockwise. if dir=1, the arc is oriented counterclockwise, else it is clockwise (i.e. clockwise from p2 to p1). If dir=0, we take it to be unoriented
GEN arc_init(GEN c, GEN p1, GEN p2, int dir, long prec){
  pari_sp top=avma;
  GEN ang2=radialangle(c, p2, gen_0, prec);//No tolerance need
  GEN arc=cgetg(ARCLEN, t_VEC);
  gel(arc, 1)=gcopy(gel(c, 1));
  gel(arc, 2)=gcopy(gel(c, 2));
  gel(arc, 3)=gcopy(p1);
  gel(arc, 4)=gcopy(p2);
  gel(arc, 5)=radialangle(c, p1, gen_0, prec);//start angle
  gel(arc, 6)=shiftangle(ang2, gel(arc, 5), gen_0, prec);//end angle; no need for tolerance.
  if(dir==1) gel(arc, 7)=gen_1;
  else if(dir==0) gel(arc, 7)=gen_0;
  else gel(arc, 7)=gen_m1;
  gel(arc, 8)=gen_0;
  return gerepileupto(top, arc);
}

//arc_init with checking of c
GEN arc_init_tc(GEN c, GEN p1, GEN p2, int dir, long prec){
  int i1=geom_check(c);
  if(i1!=0 && i1!=2) pari_err_TYPE("Please input a circle", c);
  return arc_init(c, p1, p2, dir, prec);
}

//Returns the midpoint of the arc between p1 and p2 (counterclockwise) on c.
GEN arc_midpoint(GEN c, GEN p1, GEN p2, GEN tol, long prec){
  pari_sp top=avma;
  GEN pts=circleline_int(c, perpbis(p1, p2), tol, prec);
  GEN a1=radialangle(c, p1, gen_0, prec);
  GEN a2=shiftangle(radialangle(c, p2, gen_0, prec), a1, gen_0, prec);//No tolerance concerns
  GEN angint1=shiftangle(radialangle(c, gel(pts, 1), gen_0, prec), a1, gen_0, prec);//the angle formed by pts[1] to c with base a1. Again, no tolerance need.
  if(gcmp(a2, angint1)==1) return gerepilecopy(top, gel(pts, 1));//No need for tolerance, as if this is an issue our points p1 and p2 would be equal up to tolerance.
  return gerepilecopy(top, gel(pts, 2));
}

//arc_midpoint with checking c and presetting tol
GEN arc_midpoint_tc(GEN c, GEN p1, GEN p2, long prec){
  int i1=geom_check(c);
  if(i1!=0 && i1!=2) pari_err_TYPE("Please input a circle", c);
  pari_sp top=avma;
  GEN tol=deftol(prec);
  return gerepileupto(top, arc_midpoint(c, p1, p2, tol, prec));
}

//Returns the angle of intersection between circles c1 and c2 which intersect at p, with the angle formed by rotating the tangent to c1 at p counterclockwise to the tangent to c2 at p.
GEN circle_angle(GEN c1, GEN c2, GEN p, GEN tol, long prec){
  pari_sp top=avma;
  GEN s1=circle_tangentslope(c1, p, prec);
  GEN s2=circle_tangentslope(c2, p, prec);
  GEN ang=anglediff(atanoo(s2, prec), atanoo(s1, prec), tol, prec), pi=mppi(prec);//Difference in angles in [0, 2*Pi]
  int topi=tolcmp(ang, pi, tol, prec);//We want to be in the range [0, Pi), so potentially subtract off Pi
  if(topi==-1) return gerepileupto(top, ang);
  else if(topi==0){avma=top;return gen_0;}//Same angle
  return gerepileupto(top, gsub(ang, pi));//Subtract off pi
}

//Circle with centre cent passing through p
GEN circle_fromcp(GEN cent, GEN p, long prec){
  pari_sp top=avma;
  GEN pmcent=gsub(p, cent);
  GEN ret=cgetg(4, t_VEC);
  gel(ret, 1)=gcopy(cent);
  gel(ret, 2)=gabs(pmcent, prec);
  gel(ret, 3)=gen_0;
  return gerepileupto(top, ret);
}

//Circle through 3 points (with one allowed to be oo, making a line instead)
GEN circle_fromppp(GEN p1, GEN p2, GEN p3, GEN tol, long prec){
  if(typ(p1)==t_INFINITY) return line_frompp(p2, p3);
  if(typ(p2)==t_INFINITY) return line_frompp(p1, p3);
  if(typ(p3)==t_INFINITY) return line_frompp(p1, p2);//Lines
  pari_sp top=avma;
  GEN l1=perpbis(p1, p2), l2=perpbis(p1, p3);
  GEN centre=line_int(l1, l2, tol, prec);//centre is intersection of perp bisectors.
  if(typ(centre)==t_INFINITY) return gerepileupto(top, line_frompp(p1, p2));//p1, p2, p3 collinear
  return gerepileupto(top, circle_fromcp(centre, p1, prec));//The circle!
}

//circle_fromppp with presetting of tolerance.
GEN circle_fromppp_tc(GEN p1, GEN p2, GEN p3, long prec){
  pari_sp top=avma;
  GEN tol=deftol(prec);
  return gerepileupto(top, circle_fromppp(p1, p2, p3, tol, prec));
}

//Returns the slope of the tangent to c at p
GEN circle_tangentslope(GEN c, GEN p, long prec){
  pari_sp top=avma;
  GEN c1mp=gsub(gel(c, 1), p);//c[1]-p
  GEN c1mpr=real_i(c1mp);
  GEN c1mpi=imag_i(c1mp);
  return gerepileupto(top, divoo(c1mpr, gneg(c1mpi)));//divoo(real(c[1])-real(p),imag(p)-imag(c[1])));
}

//Checks c is a circle or circle arc and calls circle_tangentslope
GEN circle_tangentslope_tc(GEN c, GEN p, long prec){
  int i1=geom_check(c);
  if(i1!=0 && i1!=2) pari_err_TYPE("Please input a circle", c);
  return circle_tangentslope(c, p, prec);
}

//Angle between l1, l2, formed by rotating l1 counterclockwise (result is in [0, Pi)
GEN line_angle(GEN l1, GEN l2, long prec){
  pari_sp top=avma;
  return gerepileupto(top, gmod(gsub(atanoo(gel(l2, 1), prec), atanoo(gel(l1, 1), prec)), mppi(prec)));
}

//The line through p with slope s
GEN line_fromsp(GEN s, GEN p){
  if(typ(s)==t_INFINITY){//oo slope
    GEN ret=cgetg(4, t_VEC);
    gel(ret, 1)=mkoo();
    gel(ret, 2)=greal(p);//x-intercept
    gel(ret, 3)=gen_1;
    return ret;
  }
  pari_sp top=avma;
  GEN strealp=gmul(s, real_i(p)), imagp=imag_i(p);
  GEN ret=cgetg(4, t_VEC);
  gel(ret, 1)=gcopy(s);
  gel(ret, 2)=gsub(imagp, strealp);//y=intercept
  gel(ret, 3)=gen_1;
  return gerepileupto(top, ret);
}

//Line through two points
GEN line_frompp(GEN p1, GEN p2){
  pari_sp top=avma;
  return gerepileupto(top, line_fromsp(slope(p1, p2), p1));
}

//Midpoint of p1 and p2
GEN midpoint(GEN p1, GEN p2){
  pari_sp top=avma;
  return gerepileupto(top, gdivgs(gadd(p1, p2), 2));
}

//Returns M(c), for c a circle/line/arc/segment
GEN mobius(GEN M, GEN c, GEN tol, long prec){
  int i=geom_check(c);
  switch(i){
    case 0: return mobius_circle(M, c, tol, prec);
    case 1: return mobius_line(M, c, tol, prec);
    case 2: return mobius_arcseg(M, c, 1, tol, prec);
    case 3: return mobius_arcseg(M, c, 0, tol, prec);
  }
  pari_err_TYPE("Please input a circle/line/arc/segment", c);
  return gen_0;//Never reach this far
}

//mobius, where we check the type of M and set the default tolerance
GEN mobius_tc(GEN M, GEN c, long prec){
  if(typ(M)!=t_MAT) pari_err_TYPE("M should be a 2x2 matrix", M);
  if(lg(M)!=3 || lg(gel(M,1))!=3) pari_err(e_DIM,"M should be a 2x2 matrix");
  pari_sp top=avma;
  GEN tol=deftol(prec);
  return gerepileupto(top, mobius(M, c, tol, prec));
}

//Mobius map acting on an arc or segment
static GEN mobius_arcseg(GEN M, GEN c, int isarc, GEN tol, long prec){
  pari_sp top=avma;//We start by finding the new circle/line(ignoring the arc/segment bits)
  GEN endpt1, endpt2, extpt;//Endpoints and an extra point (used to have 3 points to translate)
  endpt1=mat_eval(M, gel(c, 3));
  endpt2=mat_eval(M, gel(c, 4));
  //Now we must find an extra point extpt
  if(isarc==1) extpt=mat_eval(M, arc_midpoint(c, gel(c, 3), gel(c, 4), tol, prec));//arc
  else{//segment
    if(typ(gel(c, 3))==t_INFINITY){//Start point infinity
      GEN u;
      if(gequal(gel(c, 6), gen_1)) u=gen_m2;//segment goes vertically up or right
      else u=gen_2;//segment goes vertically down or left
      if(typ(gel(c, 1))==t_INFINITY) extpt=mat_eval(M, gadd(gel(c, 4), gmul(u, gen_I())));//Vertical line, M(c[4]+U*I)
      else extpt=mat_eval(M, gadd(gel(c, 4), gmul(u, gaddsg(1, gmul(gel(c, 1), gen_I())))));//non-vertical line, M(c[4]+u+u*c[1]*I)
    }
    else if(typ(gel(c, 4))==t_INFINITY){//End point infinity
      GEN u;
      if(gequal(gel(c, 6), gen_1)) u=gen_2;//segment goes vertically up or right
      else u=gen_m2;//segment goes vertically down or left
      if(typ(gel(c, 1))==t_INFINITY) extpt=mat_eval(M, gadd(gel(c, 3), gmul(u, gen_I())));//Vertical line, M(c[3]+U*I)
      else extpt=mat_eval(M, gadd(gel(c, 3), gmul(u, gaddsg(1, gmul(gel(c, 1), gen_I())))));//non-vertical line, M(c[3]+u+u*c[1]*I)
    }
    else{//Start/end points in the plane
      if(gequal(gel(c, 7), gen_1)) extpt=mat_eval(M, midpoint(gel(c, 3), gel(c, 4)));//Does not go through oo, can take midpoint
      else extpt=mat_eval(M, mkoo());//Use oo, since the line goes through there
    }
  }
  //Now we have the 3 new points used to define the new arc/segment. Let's finish the process.
  GEN ret=cgetg(ARCLEN, t_VEC);//The returned arc/segment
  GEN newcirc=circle_fromppp(endpt1, endpt2, extpt, tol, prec);//The new circle/line
  gel(ret, 1)=gel(newcirc, 1);//Slope/centre
  gel(ret, 2)=gel(newcirc, 2);//x or y intercept/radius
  gel(ret, 3)=endpt1;//Start point
  gel(ret, 4)=endpt2;//end point. These may be in the wrong order, so will fix later if so.
  if(gequal0(gel(newcirc, 3))){//Circle
    gel(ret, 8)=gen_0;//arc
    gel(ret, 5)=radialangle(newcirc, endpt1, tol, prec);//angle 1
    gel(ret, 6)=shiftangle(radialangle(newcirc, endpt2, tol, prec), gel(ret, 5), tol, prec);//angle 2
    if(isarc) gel(ret, 7)=gel(c, 7);//Temporary; match the old to the new direction
    else gel(ret, 7)=gen_1;//Temporary; since for segments we have a bona fide start and end, we start with the direction being 1.
    if(!onarc(ret, extpt, tol, prec)){//Must swap start/endpoints, angles and the direction
      gel(ret, 3)=endpt2;
      gel(ret, 4)=endpt1;
      GEN tempang=gel(ret, 5);//angle to endpt1
      gel(ret, 5)=shiftangle(gel(ret, 6), gen_0, tol, prec);//The angle to endpt2 shifted to [0, 2*Pi)
      gel(ret, 6)=shiftangle(tempang, gel(ret, 5), tol, prec);//Angle to endpt1 shifted with a base of the angle to endpt2
      gel(ret, 7)=gneg(gel(ret, 7));//We now run backwards
    }
  }
  else{//Line
    if(isarc==1 && gequal(gel(c, 7), gen_m1)){//We need to reverse the order of the points, because the arc ran backwards.
      gel(ret, 3)=endpt2;
      gel(ret, 4)=endpt1;
    }
    gel(ret, 8)=gen_1; //segment
    gel(ret, 5)=gen_0;//Unused
    if(typ(endpt1)==t_INFINITY || typ(endpt2)==t_INFINITY){//oo endpoint
      gel(ret, 6)=gen_1;//Temporary, assume we go up/right
      gel(ret, 7)=gen_0;//Not used
      if(!onseg(ret, extpt, tol, prec)) gel(ret, 6)=gen_m1;//We were wrong, and go down/left.
    }
    else{//Both endpoints finite
      gel(ret, 6)=gen_0;
      gel(ret, 7)=gen_1;//Temporary, assume we go through the plane only
      if(!onseg(ret, extpt, tol, prec)) gel(ret, 7)=gen_m1;//We were wrong, and go through oo.
    }
  }
  return gerepilecopy(top, ret);
}

//Mobius map acting on circle
static GEN mobius_circle(GEN M, GEN c, GEN tol, long prec){
  pari_sp top=avma;
  GEN p1=mat_eval(M, gadd(gel(c, 1), gel(c, 2)));//M(c[1]+c[2])
  GEN p2=mat_eval(M, gadd(gel(c, 1), gmul(gel(c, 2), gen_I())));//M(c[1]+c[2]*I)
  GEN p3=mat_eval(M, gsub(gel(c, 1), gel(c, 2)));//M(c[1]-c[2])
  return gerepileupto(top, circle_fromppp(p1, p2, p3, tol, prec));
}

//Mobius map acting on line
static GEN mobius_line(GEN M, GEN l, GEN tol, long prec){
  pari_sp top=avma;
  GEN p1, p2, p3, I=gen_I();
  if(typ(gel(l, 1))==t_INFINITY){//Vertical line
    p1=mat_eval(M, gel(l, 2));//M(x-intercept)
    p2=mat_eval(M, gadd(gel(l, 2), I));//M(x-intercept+I)
    p3=mat_eval(M, gsub(gel(l, 2), I));//M(x-intercept-I)
  }
  else{//Non-vertical line
    GEN slopeIp1=gaddgs(gmul(gel(l, 1), I), 1);//1+Slope*I
    GEN p1base=gmul(gel(l, 2), I);//y-intercept
    GEN p2base=gadd(p1base, slopeIp1);//y-intercept+1+slope*I
    GEN p3base=gadd(p2base, slopeIp1);//y-intercept+2+2*slope*I
    p1=mat_eval(M, p1base);p2=mat_eval(M, p2base);p3=mat_eval(M, p3base);
  }
  return gerepileupto(top, circle_fromppp(p1, p2, p3, tol, prec));
}

//Perpendicular bisector of distinct points
GEN perpbis(GEN p1, GEN p2){
  pari_sp top=avma;
  return gerepileupto(top, line_fromsp(divoo(gen_m1, slope(p1, p2)), midpoint(p1, p2)));
}

//Angle between p and the centre of c, in the range [0, 2*Pi)
GEN radialangle(GEN c, GEN p, GEN tol, long prec){
  pari_sp top=avma;
  return gerepileupto(top, shiftangle(garg(gsub(p, gel(c, 1)), prec), gen_0, tol, prec));
}

//Radialangle with default tolerance.
GEN radialangle_tc(GEN c, GEN p, long prec){
  int i1=geom_check(c);
  if(i1!=0 && i1!=2) pari_err_TYPE("Please input a circle/arc", c);
  pari_sp top=avma;
  GEN tol=deftol(prec);
  return gerepileupto(top, radialangle(c, p, tol, prec));
}

/*
slope, intercept, startpt, endpt, 0, ooendptor, dir, 1

GEN seg_init(GEN p1, GEN p2, int thruoo, long prec){
  
  GEN l=line_frompp(p1, p2);
}*/

//The slope of the line through p1, p2
GEN slope(GEN p1, GEN p2){
  pari_sp top=avma;
  return gerepileupto(top, divoo(imag_i(gsub(p2, p1)), real_i(gsub(p2, p1))));
}


//INTERSECTION OF LINES/CIRCLES

//Returns the intersection points of two arcs
GEN arc_int(GEN c1, GEN c2, GEN tol, long prec){
  pari_sp top=avma;
  GEN ipts=circle_int(c1, c2, tol, prec);
  if(lg(ipts)==1){avma=top;return cgetg(1, t_VEC);}//No intersection
  if(lg(ipts)==2){//One intersection point (tangent circles)
    if(!onarc(c1, gel(ipts, 1), tol, prec)){avma=top;return cgetg(1, t_VEC);}//Not on arc 1
    if(!onarc(c2, gel(ipts, 1), tol, prec)){avma=top;return cgetg(1, t_VEC);}//Not on arc 2
    return gerepilecopy(top, ipts);//On arc
  }
  //Two intersections
  int i1=onarc(c1, gel(ipts, 1), tol, prec);
  if(i1==1) i1=onarc(c2, gel(ipts, 1), tol, prec);//Now i1==1 iff the ipts[1] is on both c1 and c2
  int i2=onarc(c1, gel(ipts, 2), tol, prec);
  if(i2==1) i2=onarc(c2, gel(ipts, 2), tol, prec);//Now i2==1 iff the ipts[2] is on both c1 and c2
  if(i1==1){
    if(i2==1) return gerepilecopy(top, ipts);//Both pts on the arcs
    GEN ret=cgetg(2, t_VEC);//Just point 1
    gel(ret, 1)=gcopy(gel(ipts, 1));
    return gerepileupto(top, ret);
  }
  //Now i1=0
  if(i2==0){avma=top;return cgetg(1, t_VEC);}//Not on either arc
  GEN ret=cgetg(2, t_VEC);//Just point 2
  gel(ret, 1)=gcopy(gel(ipts, 2));
  return gerepileupto(top, ret);
}

//arc_int with typecheck
GEN arc_int_tc(GEN c1, GEN c2, long prec){
  int i1=geom_check(c1), i2=geom_check(c2);
  if(i1!=2 || i2!=2) pari_err_TYPE("Please input circle arcs", c1);
  pari_sp top=avma;
  GEN tol=deftol(prec);
  return gerepileupto(top, arc_int(c1, c2, tol, prec));
}

//Returns the intersection points of an arc and a segment
GEN arcseg_int(GEN c, GEN l, GEN tol, long prec){
  pari_sp top=avma;
  GEN ipts=circleline_int(c, l, tol, prec);
  if(lg(ipts)==1){avma=top;return cgetg(1, t_VEC);}//No intersection
  if(lg(ipts)==2){//One intersection point (tangent circles)
    if(!onarc(c, gel(ipts, 1), tol, prec)){avma=top;return cgetg(1, t_VEC);}//Not on arc
    if(!onseg(l, gel(ipts, 1), tol, prec)){avma=top;return cgetg(1, t_VEC);}//Not on segment
    return gerepilecopy(top, ipts);//On both
  }
  //Two intersections
  int i1=onarc(c, gel(ipts, 1), tol, prec);
  if(i1==1) i1=onseg(l, gel(ipts, 1), tol, prec);//Now i1==1 iff the ipts[1] is on both c and l
  int i2=onarc(c, gel(ipts, 2), tol, prec);
  if(i2==1) i2=onseg(l, gel(ipts, 2), tol, prec);//Now i2==1 iff the ipts[2] is on both c and l
  if(i1==1){
    if(i2==1) return gerepilecopy(top, ipts);//Both pts on both
    GEN ret=cgetg(2, t_VEC);//Just point 1
    gel(ret, 1)=gcopy(gel(ipts, 1));
    return gerepileupto(top, ret);
  }
  //Now i1=0
  if(i2==0){avma=top;return cgetg(1, t_VEC);}//Not on either
  GEN ret=cgetg(2, t_VEC);//Just point 2
  gel(ret, 1)=gcopy(gel(ipts, 2));
  return gerepileupto(top, ret);
}

//arc_int with typecheck
GEN arcseg_int_tc(GEN c, GEN l, long prec){
  int i1=geom_check(c), i2=geom_check(l);
  if(i1!=2 || i2!=3) pari_err_TYPE("Please input a circle arc and line segment", c);
  pari_sp top=avma;
  GEN tol=deftol(prec);
  return gerepileupto(top, arcseg_int(c, l, tol, prec));
}

//Returns the set of points in the intersection of circles c1, c2
GEN circle_int(GEN c1, GEN c2, GEN tol, long prec){
  pari_sp top=avma;
  GEN a1=real_i(gel(c1, 1)), b1=imag_i(gel(c1, 1)), r1=gel(c1, 2);//x, y coords and radius of c1
  GEN a2=real_i(gel(c2, 1)), b2=imag_i(gel(c2, 1)), r2=gel(c2, 2);//x, y coords and radius of c2
  GEN a1ma2=gsub(a1, a2), b1mb2=gsub(b1, b2), x1, x2, y1, y2;
  int oneint=0;
  if(gcmp(gabs(a1ma2, prec), gabs(b1mb2, prec))>=0){//We want to divide by the larger of the two quantities to maximize precision and avoid errors when the centres are on the same line.
    if(toleq(a1ma2, gen_0, tol, prec)==1){avma=top;return cgetg(1, t_VEC);}//Same centre, cannot intersect.
    //u=(r1^2-r2^2+b2^2-b1^2+a2^2-a1^2)/(2*a2-2*a1)-a1;
    GEN u=gsub(gdiv(gsub(gadd(gsub(gadd(gsub(gsqr(r1), gsqr(r2)), gsqr(b2)), gsqr(b1)), gsqr(a2)), gsqr(a1)), gmulgs(a1ma2, -2)), a1);
    GEN v=gneg(gdiv(b1mb2, a1ma2));//v=(b1-b2)/(a2-a1), and x=a1+u+vy
    GEN uvmb1=gsub(gmul(u, v), b1);//uv-b1
    GEN vsqrp1=gaddgs(gsqr(v), 1);//v^2+1
    GEN rtpart=gsub(gsqr(uvmb1), gmul(vsqrp1, gadd(gsqr(b1), gsub(gsqr(u), gsqr(r1)))));//(u*v-b1)^2-(v^2+1)*(b1^2+u^2-r1^2)
    oneint=tolcmp(rtpart, gen_0, tol, prec);//Comparing rtpart to 0
    if(oneint==-1){avma=top;return cgetg(1, t_VEC);}//rtpart must be square rooted, so if it's negative the circles do not intersect
    if(oneint==0){//One intersection, so we take rtpart=0. This is CRUCIAL, as taking the square root kills our precision if we don't do this here.
      y1=gdiv(gneg(uvmb1), vsqrp1);//y1=(b1-u*v)/(1*v^2+1)
      x1=gadd(gadd(a1, u), gmul(v, y1));//x1=a1+u+v*y1
    }
    else{
      y1=gdiv(gadd(gneg(uvmb1), gsqrt(rtpart, prec)), vsqrp1);//y1=(b1-u*v+sqrt((u*v-b1)^2-*(v^2+1)*(b1^2+u^2-r1^2)))/(1*v^2+1)
      y2=gadd(gneg(y1), gdiv(gmulgs(uvmb1, -2), vsqrp1));//y2=-y1+(2*b1-2*u*v)/(v^2+1)
      GEN a1pu=gadd(a1, u);
      x1=gadd(a1pu, gmul(v, y1));//x1=a1+u+v*y1
      x2=gadd(a1pu, gmul(v, y2));//x1=a1+u+v*y2
    }
  }
  else{
    if(toleq(b1mb2, gen_0, tol, prec)==1){avma=top;return cgetg(1, t_VEC);}//Same centre, cannot intersect.
    //u=(r1^2-r2^2+b2^2-b1^2+a2^2-a1^2)/(2*b2-2*b1)-b1;
    GEN u=gsub(gdiv(gsub(gadd(gsub(gadd(gsub(gsqr(r1), gsqr(r2)), gsqr(b2)), gsqr(b1)), gsqr(a2)), gsqr(a1)), gmulgs(b1mb2, -2)), b1);
    GEN v=gneg(gdiv(a1ma2, b1mb2));//v=(a1-a2)/(b2-b1), and y=b1+u+vx
    GEN uvma1=gsub(gmul(u, v), a1);//uv-a1
    GEN vsqrp1=gaddgs(gsqr(v), 1);//v^2+1
    GEN rtpart=gsub(gsqr(uvma1), gmul(vsqrp1, gadd(gsqr(a1), gsub(gsqr(u), gsqr(r1)))));//(u*v-a1)^2-(v^2+1)*(a1^2+u^2-r1^2))
    oneint=tolcmp(rtpart, gen_0, tol, prec);//Comparing rtpart to 0
    if(oneint==-1){avma=top;return cgetg(1, t_VEC);}//rtpart must be square rooted, so if it's negative the circles do not intersect
    if(oneint==0){//One intersection, so we take rtpart=0. This is CRUCIAL, as taking the square root kills our precision if we don't do this here.
      x1=gdiv(gneg(uvma1), vsqrp1);//x1=(a1-u*v)/(v^2+1);
      y1=gadd(gadd(b1, u), gmul(v, x1));//y1=b1+u+v*x1;
    }
    else{
      x1=gdiv(gadd(gneg(uvma1), gsqrt(rtpart, prec)), vsqrp1);//x1=(a1-u*v+sqrt((u*v-a1)^2-(v^2+1)*(a1^2+u^2-r1^2)))/(v^2+1);
      x2=gadd(gneg(x1), gdiv(gmulgs(uvma1, -2), vsqrp1));//x2=-x1+(2*a1-2*u*v)/(v^2+1)
      GEN b1pu=gadd(b1, u);
      y1=gadd(b1pu, gmul(v, x1));//y1=b1+u+v*x1;
      y2=gadd(b1pu, gmul(v, x2));//y2=b1+u+v*x2;
    }
  }
  if(oneint==0){//One point of intersection (0 pts of intersection was already dealt with and returned)
    GEN y1I=gmul(gen_I(), y1);
    GEN ret=cgetg(2, t_VEC);
    gel(ret, 1)=gadd(x1, y1I);
    return gerepileupto(top, ret);
  }
  GEN y1I=gmul(gen_I(), y1), y2I=gmul(gen_I(), y2);
  GEN ret=cgetg(3, t_VEC);
  gel(ret, 1)=gadd(x1, y1I);
  gel(ret, 2)=gadd(x2, y2I);
  return gerepileupto(top, ret);
}

//Checks that c1, c2 are circles and calls circle_int with the default tolerance.
GEN circle_int_tc(GEN c1, GEN c2, long prec){
  int i1=geom_check(c1), i2=geom_check(c2);
  if(i1!=0 || i2!=0) pari_err_TYPE("Please input circles", c1);
  pari_sp top=avma;
  GEN tol=deftol(prec);
  return gerepileupto(top, circle_int(c1, c2, tol, prec));
}

//Returns the intersection points of c and l
GEN circleline_int(GEN c, GEN l, GEN tol, long prec){
  pari_sp top=avma;
  if(typ(gel(l, 1))==t_INFINITY){
    GEN x1=gel(l, 2);
    GEN rtpart=gsub(gsqr(gel(c, 2)), gsqr(gsub(x1, real_i(gel(c, 1)))));//c[2]^2-(x1-real(c[1]))^2
    if(gsigne(rtpart)==-1){avma=top;return cgetg(1, t_VEC);}//No intersections.
    GEN y1=gadd(imag_i(gel(c, 1)), gsqrt(rtpart, prec));//y1=imag(c[1])+sqrt(c[2]^2-(x1-real(c[1]))^2)
    if(toleq(rtpart, gen_0, tol, prec)){//Only one intersection point
      GEN ret=cgetg(2, t_VEC);
      gel(ret, 1)=cgetg(3, t_COMPLEX);
      gel(gel(ret, 1), 1)=gcopy(x1);
      gel(gel(ret, 1), 2)=gcopy(y1);
      return gerepileupto(top, ret);
    }
    //Two intersection points
    GEN y1py2=gmulgs(imag_i(gel(c, 1)), 2);//2*imag(c[1])
    GEN ret=cgetg(3, t_VEC);
    gel(ret, 1)=cgetg(3, t_COMPLEX);gel(ret, 2)=cgetg(3, t_COMPLEX);
    gel(gel(ret, 1), 1)=gcopy(x1);
    gel(gel(ret, 1), 2)=gcopy(y1);
    gel(gel(ret, 2), 1)=gcopy(x1);
    gel(gel(ret, 2), 2)=gsub(y1py2, y1);
    return gerepileupto(top, ret);
  }
  //Now y=mx+b with m finite
  GEN A=gaddgs(gsqr(gel(l, 1)), 1);//l[1]^2+1
  GEN l2mic1=gsub(gel(l, 2), imag_i(gel(c, 1)));//l[2]-imag(c[1])
  GEN B=gadd(gmulgs(real_i(gel(c, 1)), -2), gmulsg(2, gmul(gel(l, 1), l2mic1)));//-2*real(c[1])+2*l[1]*(l[2]-imag(c[1]))
  GEN C=gadd(gsqr(real_i(gel(c, 1))), gsub(gsqr(l2mic1), gsqr(gel(c, 2))));//real(c[1])^2+(l[2]-imag(c[1]))^2-c[2]^2
  GEN rtpart=gsub(gsqr(B), gmulsg(4, gmul(A, C)));
  int rtpartsig=tolcmp(rtpart, gen_0, tol, prec);
  if(rtpartsig==-1){avma=top;return cgetg(1, t_VEC);}//No intersection
  if(rtpartsig==0){//One root, and rtpart=0
    GEN x1=gdiv(B, gmulgs(A, -2));//-B/(2A)
    GEN y1part=gmul(gel(l, 1), x1);//l[1]*x1
    GEN ret=cgetg(2, t_VEC);
    gel(ret, 1)=cgetg(3, t_COMPLEX);
    gel(gel(ret, 1), 1)=gcopy(x1);
    gel(gel(ret, 1), 2)=gadd(y1part, gel(l, 2));//y1=l[1]*x1+l[2];
    return gerepileupto(top, ret);
  }
  //Two roots
  GEN x1=gdiv(gsub(gsqrt(rtpart, prec), B), gmulgs(A, 2));//x1=(-B+sqrt(B^2-4*A*C))/(2*A);
  GEN x2=gsub(gneg(gdiv(B, A)), x1);//-B/A-x1
  GEN y1part=gmul(gel(l, 1), x1);//l[1]*x1
  GEN y2part=gmul(gel(l, 1), x2);//l[1]*x2
  GEN ret=cgetg(3, t_VEC);
  gel(ret, 1)=cgetg(3, t_COMPLEX);gel(ret, 2)=cgetg(3, t_COMPLEX);
  gel(gel(ret, 1), 1)=gcopy(x1);
  gel(gel(ret, 1), 2)=gadd(y1part, gel(l, 2));//l[1]*x1+l[2];
  gel(gel(ret, 2), 1)=gcopy(x2);
  gel(gel(ret, 2), 2)=gadd(y2part, gel(l, 2));//l[1]*x2+l[2];
  return gerepileupto(top, ret);
}

//Checks that c is a circle and l is a line and calls circleline_int with the default tolerance.
GEN circleline_int_tc(GEN c, GEN l, long prec){
  int i1=geom_check(c), i2=geom_check(l);
  if(i1!=0 || i2!=1) pari_err_TYPE("Please input a line and circle", c);
  pari_sp top=avma;
  GEN tol=deftol(prec);
  return gerepileupto(top, circleline_int(c, l, tol, prec));
}

//Returns the intersection of two circles/lines/arcs/segments (in any combination)
GEN genseg_int(GEN s1, GEN s2, GEN tol, long prec){
  int t1=geom_check(s1), t2=geom_check(s2);
  if(t1==-1 || t2==-1)  pari_err_TYPE("Please input two circles/lines/arcs/segments", s1);
  switch(t1){
    case 0:
      switch(t2){
        case 0: return circle_int(s1, s2, tol, prec);//Circle circle
        case 1: return circleline_int(s1, s2, tol, prec);//Circle line
        case 2: return arc_int(s1, s2, tol, prec);//Circle arc
        case 3: return arcseg_int(s1, s2, tol, prec);//Circle segment
      }
    case 1:
      switch(t2){
        case 0: return circleline_int(s2, s1, tol, prec);//Line circle
        case 1: return line_int(s1, s2, tol, prec);//Line line
        case 2: return arcseg_int(s2, s1, tol, prec);//Line arc
        case 3: return seg_int(s1, s2, tol, prec);//Line segment
      }
    case 2:
      switch(t2){
        case 0: return arc_int(s1, s2, tol, prec);//Arc circle
        case 1: return arcseg_int(s1, s2, tol, prec);//Arc line
        case 2: return arc_int(s1, s2, tol, prec);//Arc arc
        case 3: return arcseg_int(s1, s2, tol, prec);//Arc segment
      }
    case 3:
      switch(t2){
        case 0: return arcseg_int(s2, s1, tol, prec);//segment circle
        case 1: return seg_int(s1, s2, tol, prec);//segment line
        case 2: return arcseg_int(s2, s1, tol, prec);//segment arc
        case 3: return seg_int(s1, s2, tol, prec);//segment seg
      }
  }
  pari_err_TYPE("ERROR: UNEXPECTED INPUT TYPES", mkvec2(s1, s2));
  return gen_0;
}

//genseg_int with default tolerance.
GEN genseg_int_tc(GEN s1, GEN s2, long prec){
  pari_sp top=avma;
  GEN tol=deftol(prec);
  return gerepileupto(top, genseg_int(s1, s2, tol, prec));
}

//The intersection of two lines
GEN line_int(GEN l1, GEN l2, GEN tol, long prec){
  GEN s1=gel(l1, 1), s2=gel(l2, 1);//Slopes
  if(toleq(s1, s2, tol, prec)) return mkoo();//Parallel or equal
  pari_sp top=avma;
  if(typ(s1)==t_INFINITY){//l1 vertical
    GEN ypart=gmul(s2, gel(l1, 2));//s2*l1[2]
    GEN ipt=cgetg(3, t_COMPLEX);
    gel(ipt, 1)=gcopy(gel(l1, 2));
    gel(ipt, 2)=gadd(ypart, gel(l2, 2));//s2*l1[2]+l2[2]
    return gerepileupto(top, ipt);
  }
  if(typ(s2)==t_INFINITY){//l2 vertical
    GEN ypart=gmul(s1, gel(l2, 2));//s1*l2[2]
    GEN ipt=cgetg(3, t_COMPLEX);
    gel(ipt, 1)=gcopy(gel(l2, 2));
    gel(ipt, 2)=gadd(ypart, gel(l1, 2));//s1*l2[2]+l1[2]
    return gerepileupto(top, ipt);
  }
  GEN x=gdiv(gsub(gel(l2, 2), gel(l1, 2)), gsub(s1, s2));//(l2[2]-l1[2])/(s1-s2)
  GEN ypart=gmul(s1, x);//s1*x
  GEN ipt=cgetg(3, t_COMPLEX);
  gel(ipt, 1)=gcopy(x);
  gel(ipt, 2)=gadd(ypart, gel(l1, 2));//s1*x+l1[2]
  return gerepileupto(top, ipt);
}

//Checks that l1, l2 are lines and calls line_int with the default tolerance.
GEN line_int_tc(GEN l1, GEN l2, long prec){
  int i1=geom_check(l1), i2=geom_check(l2);
  if(i1!=1 || i2!=1) pari_err_TYPE("Please input lines", l1);
  pari_sp top=avma;
  GEN tol=deftol(prec);
  return gerepileupto(top, line_int(l1, l2, tol, prec));
}

//p is assumed to be on the circle defined by c; this checks if it is actually on the arc (running counterclockwise from c[3] to c[4]).
int onarc(GEN c, GEN p, GEN tol, long prec){
  if(lg(c)==CIRCLEN) return 1;//Allow input of just a circle, so the return is trivially 1
  if(toleq(gel(c, 3), p, tol, prec)) return 1;//p=the start point. We have this done seperately in case rounding errors take the angle to <c[5], as this will cause issues with the shifting angle.
  pari_sp top=avma;
  GEN angle=shiftangle(radialangle(c, p, tol, prec), gel(c, 5), tol, prec);//Getting the angle in the range [c[5], c[5]+2*Pi)
  if(tolcmp(angle, gel(c, 6), tol, prec)<=0){avma=top;return 1;}//On the arc
  avma=top;
  return 0;//Beyond the arc.
}

//onarc with checking c and passing the default tolerance.
int onarc_tc(GEN c, GEN p, long prec){
  int i1=geom_check(c);
  if(i1!=0 && i1!=2) pari_err_TYPE("Please input a circle arc and point", c);
  pari_sp top=avma;
  GEN tol=deftol(prec);
  int oa=onarc(c, p, tol, prec);
  avma=top;
  return oa;
}

//p is assumed to be on the line defined by l; this checks if it is actually on the segment l
int onseg(GEN l, GEN p, GEN tol, long prec){
  if(lg(l)==CIRCLEN) return 1;//Allow input of a line, so return is trivially 1
  if(typ(p)==t_INFINITY){//p is the point at oo
    if(!gequal0(gel(l, 6)) || gequal(gel(l, 7), gen_m1)) return 1;//oo is an endpoint OR the seg passes through oo
    return 0;//If not, does not pass through oo
  }
  pari_sp top=avma;
  //Okay, now p is not oo and l is a line segment
  if(typ(gel(l, 1))==t_INFINITY){//Vertical line
    if(typ(gel(l, 3))==t_INFINITY){//Start point in oo
      if(equali1(gel(l, 6))){//p must lie BELOW l[4]
          if(tolcmp(imag_i(p), imag_i(gel(l, 4)), tol, prec)<=0){avma=top;return 1;}//Lies below l[4]
          avma=top;return 0;//lies above l[4]
      }
      //p must lie ABOVE l[4]
      if(tolcmp(imag_i(p), imag_i(gel(l, 4)), tol, prec)>=0){avma=top;return 1;}//Lies above l[4]
      avma=top;return 0;//lies below l[4]
    }
    if(typ(gel(l, 4))==t_INFINITY){//End point is oo
      if(equali1(gel(l, 6))){//p must lie ABOVE l[3]
          if(tolcmp(imag_i(p), imag_i(gel(l, 3)), tol, prec)>=0){avma=top;return 1;}//Lies above l[3]
          avma=top;return 0;//lies below l[3]
      }
      //p must lie BELOW l[3]
      if(tolcmp(imag_i(p), imag_i(gel(l, 3)), tol, prec)<=0){avma=top;return 1;}//Lies below l[3]
      avma=top;return 0;//lies above l[3]
    }
    //Start and end points are finite
    int i1=tolcmp(imag_i(gsub(p, gel(l, 3))), gen_0, tol, prec);//sign of imag(p)-imag(l[3])
    int i2=tolcmp(imag_i(gsub(p, gel(l, 4))), gen_0, tol, prec);//sign of imag(p)-imag(l[4])
    avma=top;
    if(i1==0 || i2==0) return 1;//endpoint
    if(i1==i2){//p on the same side of l[3] and l[4], so return 1 iff l passes through oo
        if(gequal(gel(l, 7), gen_1)) return 0;//Not through oo
        return 1;//through oo
    }
    //p is between l[3] and l[4], so return 1 iff l does not pass through oo
    if(gequal(gel(l, 7), gen_1)) return 1;//not through oo
    return 0;//through oo
  }
  //Non-vertical line
  if(typ(gel(l, 3))==t_INFINITY){//Start point in oo
    if(equali1(gel(l, 6))){//p must lie LEFT OF l[4]
      if(tolcmp(real_i(p), real_i(gel(l, 4)), tol, prec)<=0){avma=top;return 1;}//Lies left of l[4]
      avma=top;return 0;//lies right of  l[4]
    }
    //p must lie RIGHT OF l[4]
    if(tolcmp(real_i(p), real_i(gel(l, 4)), tol, prec)>=0){avma=top;return 1;}//Lies right of l[4]
    avma=top;return 0;//lies left of l[4]
  }
  if(typ(gel(l, 4))==t_INFINITY){//End point is oo
    if(equali1(gel(l, 6))){//p must lie RIGHT OF l[3]
      if(tolcmp(real_i(p), real_i(gel(l, 3)), tol, prec)>=0){avma=top;return 1;}//Lies right of l[3]
      avma=top;return 0;//lies below l[3]
    }
    //p must lie LEFT OF l[3]
    if(tolcmp(real_i(p), real_i(gel(l, 3)), tol, prec)<=0){avma=top;return 1;}//Lies left of l[3]
    avma=top;return 0;//lies above l[3]
  }
  //Start and end points are finite
  int i1=tolcmp(real_i(gsub(p, gel(l, 3))), gen_0, tol, prec);//sign of real(p)-real(l[3])
  int i2=tolcmp(real_i(gsub(p, gel(l, 4))), gen_0, tol, prec);//sign of real(p)-real(l[4])
  avma=top;
  if(i1==0 || i2==0) return 1;//endpoint
  if(i1==i2){//p on the same side of l[3] and l[4], so return 1 iff l passes through oo
    if(gequal(gel(l, 7), gen_1)) return 0;//Not through oo
    return 1;//through oo
  }
  //p is between l[3] and l[4], so return 1 iff l does not pass through oo
  if(gequal(gel(l, 7), gen_1)) return 1;//not through oo
  return 0;//through oo
}

//onseg with checking l and passing the default tolerance.
int onseg_tc(GEN l, GEN p, long prec){
  int i1=geom_check(l);
  if(i1!=1 && i1!=3) pari_err_TYPE("Please input a line segment and point", l);
  pari_sp top=avma;
  GEN tol=deftol(prec);
  int os=onseg(l, p, tol, prec);
  avma=top;
  return os;
}

//Returns the intersection points of l1 and l2 (length 0 or 1 vector)
GEN seg_int(GEN l1, GEN l2, GEN tol, long prec){
  pari_sp top=avma;
  GEN in=line_int(l1, l2, tol, prec);//Line intersection
  if(onseg(l1, in, tol, prec) && onseg(l2, in, tol, prec)){//On the segments
    GEN ret=cgetg(2, t_VEC);
    gel(ret, 1)=gcopy(in);
    return gerepileupto(top, ret);
  }
  avma=top;
  return cgetg(1, t_VEC);//Not on at least one segment.
}

//Checks that l1, l2 are line segments and calls seg_int with the default tolerance.
GEN seg_int_tc(GEN l1, GEN l2, long prec){
  int i1=geom_check(l1), i2=geom_check(l2);
  if(i1!=3 || i2!=3) pari_err_TYPE("Please input line segments", l1);
  pari_sp top=avma;
  GEN tol=deftol(prec);
  return gerepileupto(top, seg_int(l1, l2, tol, prec));
}


//DISTANCES

//z1 and z2 are complex numbers, this computes the hyperbolic distance between them.
GEN hdist(GEN z1, GEN z2, long prec){
  pari_sp top=avma;
  GEN x1=gel(z1,1);
  GEN y1=gel(z1,2);
  GEN x2=gel(z2,1);
  GEN y2=gel(z2,2);
  GEN x=gaddsg(1,gdiv(gadd(gsqr(gsub(x2,x1)),gsqr(gsub(y2,y1))),gmul(gmulsg(2,y1),y2)));
  GEN expd=gadd(x,gsqrt(gsubgs(gsqr(x), 1), prec));
  return gerepileupto(top,glog(expd, prec));
}

//hdist with typechecking
GEN hdist_tc(GEN z1, GEN z2, long prec){
  if(typ(z1)!=t_COMPLEX || signe(gel(z1,2))!=1) pari_err_TYPE("Please input complex numbers in the upper half plane", z1);
  if(typ(z2)!=t_COMPLEX || signe(gel(z2,2))!=1) pari_err_TYPE("Please input complex numbers in the upper half plane", z2);
  return hdist(z1, z2, prec);
}

//The hyperbolic distance between z1 and z2 in the unit disc model
GEN hdist_ud(GEN z1, GEN z2, long prec){
  pari_sp top=avma;
  GEN a = gabs(gsubsg(1, gmul(z1, conj_i(z2))), prec);//|1-z1*conj(z2)|
  GEN b = gabs(gsub(z1, z2), prec);//|z1-z2|
  GEN num=gadd(a, b);
  GEN denom=gsub(a, b);
  GEN ret;
  pari_CATCH(e_INV){
    avma=top;
    return mkoo();
  }
  pari_TRY{
    ret=gerepileupto(top, glog(gdiv(num, denom), prec));//log((a+b)/(a-b))
  }
  pari_ENDCATCH
  return ret;
}

//Given the isometric circles (in order) and the vertices (so that circles[i] intersect circles[i+1] is vertices[i]), returns the area of the convex hyperbolic polygon. If one of the sides is oo (an entry of circles is 0), the answer will be oo. If there are n vertices with angles a1,...,an, then the area is is (n-2)Pi-sum(ai)
GEN hpolygon_area(GEN circles, GEN vertices, GEN tol, long prec){
  pari_sp top=avma;
  long blen=lg(circles);
  if(blen==1 || gequal0(gel(circles, 1))) return mkoo();//No cicles or the first is 0, i.e. an infinite side
  GEN ang, area=gmulsg(blen-3, mppi(prec));//We subtract off from area.
  for(long i=1;i<blen-1;i++){
    if(gequal0(gel(circles, i+1))){avma=top;return mkoo();}//The next side is infinite
    ang=circle_angle(gel(circles, i+1), gel(circles, i), gel(vertices, i), tol, prec);//Do in opposite order to get correct angle
    area=gsub(area, ang);
  }
  ang=circle_angle(gel(circles, 1), gel(circles, blen-1), gel(vertices, blen-1), tol, prec);//The last side wraps around.
  area=gsub(area, ang);
  return gerepileupto(top, area);
}

//Normalized boundary area with typechecking
GEN hpolygon_area_tc(GEN circles, GEN vertices, long prec){
  pari_sp top=avma;
  if(typ(circles)!=t_VEC || typ(vertices)!=t_VEC) pari_err_TYPE("Please enter the edges and vertices as vectors", circles);
  if(lg(circles)!=lg(vertices)) pari_err_TYPE("There needs to be the same number of edges as vertices", vertices);
  GEN tol=deftol(prec);
  return gerepileupto(top, hpolygon_area(circles, vertices, tol, prec));
}


//FUNDAMENTAL DOMAIN COMPUTATION

/*Let Gamma be a Fuschian group; the following methods allow us to compute a fundamental domain for Gamma, following the paper of John Voight. We are assuming that we have an "exact" way to deal with the elements of Gamma. In otherwords, we need separate methods to:
i)  Multiply elements of Gamma: format as GEN eltmul(GEN *data, GEN x, GEN y), with data the extra data you need.
ii) Find the image of g in PSL(1, 1): format as GEN gamtopsl(GEN *data, GEN g, long prec), with data the extra data you need (e.g. the quaterniona algebra in the case of Shimura curves).
iii) Identify if g is trivial in Gamma: format as int istriv(GEN *data, GEN g). Need to take care as we care about being trivial in PSL, whereas most representations of Gamma are for SL (so need to check with +/-id).
iv) Pass in the identity element of Gamma and find the area of the fundamental domain. These methods are not passed in; just the values.
We pass references to these methods into the methods here.
We do computations mostly in PSU, and shift from PSL to PSU via phi(z):=(z-p)/(z-conj(p)) for some given upper half plane point p.
*/

//Some of the methods may not be that useful (i.e. assuming we have p only and not the matrices).

//Returns the edge pairing, as VECSMALL v where v[i]=j means i is paired with j. If not all edges can be paired, instead returns [v1, v2, ...] where vi is a VECSMALL that is either [gind, v1ind] or [gind, v1ind, v2ind]. gind is the index of the unpaired side, and the viind are the corresponding unpaired vertices (1 or 2). If rboth=1, returns [paired, unpaired], where this time paired is a vector of vecsmalls [i, j] meaning i is paired to j (differing to the output format when rboth=0 and everything is paired)
GEN edgepairing(GEN U, GEN tol, int rboth, long prec){
  pari_sp top=avma;
  GEN vangles=gel(U, 4);//Vertex angles
  GEN baseangle=gel(vangles, 1);//Base angle
  GEN toldata=cgetg(3, t_VEC);//Stores the necessary info for searching with tolerance (tolcmp_sort)
  gel(toldata, 1)=tol;
  gel(toldata, 2)=cgetg(2, t_VECSMALL);
  gel(toldata, 2)[1]=prec;
  long lU=lg(gel(U, 1));
  GEN unpair=vectrunc_init(lU+1), vim, vimang, pair=vectrunc_init(lU);//Unpair stores the unpaired edges, pair stores the paired edges
  long ind1, ind2, i1, i2;
  for(long i=1;i<lU;i++){
    if(gequal0(gel(gel(U, 5), i))){vectrunc_append(pair, mkvecsmall2(i, i));continue;}//oo side, go next (we say it is paired with itself)
    ind1=i;
    vim=mat_eval(gel(gel(U, 5), i), gel(gel(U, 3), ind1));//The new vertex
    vimang=shiftangle(garg(vim, prec), baseangle, tol, prec);//The new angle
    i1=gen_search(vangles, vimang, 0, &toldata, &tolcmp_sort);
    if(i1!=0) if(!toleq(vim, gel(gel(U, 3), i1), tol, prec)) i1=0;//Just because the angles are equal, the points don't have to be (though this occurence is expected to be extremely rare).
    if(i==1) ind2=lU-1;
    else ind2=i-1;//The two vertices of the side, this is the second one
    vim=mat_eval(gel(gel(U, 5), i), gel(gel(U, 3), ind2));//The second new vertex
    if(i1!=0){//If ind2 is paired, it MUST be paired to vertex i1+1
      i2=i1+1;
      if(i2==lU) i2=1;//Wrap around
      if(!toleq(vim, gel(gel(U, 3), i2), tol, prec)) i2=0;
      if(i2!=0){
        if(i<=i2) vectrunc_append(pair, mkvecsmall2(i, i2));//i<=i2 is put so that we are not pairing things twice.
      }
      else vectrunc_append(unpair, mkvecsmall2(i, ind2));//Second vertex not paired
    }
    else{//ind1 not paired
      vimang=shiftangle(garg(vim, prec), baseangle, tol, prec);//The second new angle
      i2=gen_search(vangles, vimang, 0, &toldata, &tolcmp_sort);
      if(i2!=0) if(!toleq(vim, gel(gel(U, 3), i2), tol, prec)) i2=0;//Just because the angles are equal, the points don't have to be (though this occurence is expected to be extremely rare).
      if(i2!=0) vectrunc_append(unpair, mkvecsmall2(i, ind1));//First vtx not paired
       else vectrunc_append(unpair, mkvecsmall3(i, ind1, ind2));//Neither vtx paired
    }
  }
  if(rboth) return gerepilecopy(top, mkvec2(pair, unpair));
  if(lg(unpair)==1){
    GEN pairvs=cgetg(lU, t_VECSMALL);
    for(long i=1;i<lg(pair);i++){
      pairvs[gel(pair, i)[1]]=gel(pair, i)[2];
      pairvs[gel(pair, i)[2]]=gel(pair, i)[1];
    }
    return gerepilecopy(top, pairvs);//No unpaired vertices
  }
  return gerepilecopy(top, unpair);//Unpaired vertices!
}

//Edgepairing with typecheck
GEN edgepairing_tc(GEN U, long prec){
  pari_sp top=avma;
  if(typ(U)!=t_VEC || lg(U)!=NORMBOUND) pari_err_TYPE("Please enter a normalized boundary", U);
  GEN tol=deftol(prec);
  return gerepileupto(top, edgepairing(U, tol, 1, prec));
}

//Computes the isometric circle for g, returning [g, image in PSU(1, 1), circle]. Must pass in mats (psltopsu_transmats(p)), and a method that translates g to an element of PSL(2, R).
GEN isometriccircle_mats(GEN g, GEN mats, GEN *data, GEN (*gamtopsl)(GEN *, GEN, long), GEN tol, long prec){
  pari_sp top=avma;
  GEN ginpsl=gamtopsl(data, g, prec);
  GEN ret=cgetg(4, t_VEC);
  gel(ret, 1)=gcopy(g);
  gel(ret, 2)=psltopsu_mats(ginpsl, mats);
  gel(ret, 3)=isometriccircle_psu(gel(ret, 2), tol, prec);
  if(!gequal(gel(ret, 3), gen_m1)) return gerepileupto(top, ret);//Everything A-OK
  //If we reach here, we did not have enough precision
  long newprec=prec;
  pari_CATCH(CATCH_ALL){
    avma=top;
    pari_CATCH_reset();
    pari_err(e_MISC,"Could not increase precision enough. Please increase precision/memory");
    return gen_0;
  }
  pari_TRY{
    do{
      if(newprec-prec==5) pari_err(e_MISC,"Throw");
      avma=top;
      newprec++;//Increase precision
      tol=deftol(newprec);
      if(precision(gel(mats, 3))>0){//p is inexact
        GEN p=gtofp(gel(mats, 3), newprec);
        mats=psltopsu_transmats(p);//Updating mats to more precision
      }
      ginpsl=gamtopsl(data, g, newprec);
      ret=cgetg(4, t_VEC);
      gel(ret, 1)=gcopy(g);
      gel(ret, 2)=psltopsu_mats(ginpsl, mats);
      gel(ret, 3)=isometriccircle_psu(gel(ret, 2), tol, newprec);
    }
    while(gequal(gel(ret, 3), gen_m1));
    pari_CATCH_reset();
    return gerepileupto(top, ret);
  }
  pari_ENDCATCH
}

//Not sure if useful
//Returns the isometric circle of an element of PSL
GEN isometriccircle_psl(GEN g, GEN p, GEN tol, long prec){
  pari_sp top=avma;
  GEN mats=psltopsu_transmats(p);
  return gerepileupto(top, isometriccircle_psl_mats(g, mats, tol, prec));
}

//Not sure if useful
//Given the mats to translate to PSU, this returns the isometric circle of g in the unit disc model
GEN isometriccircle_psl_mats(GEN g, GEN mats, GEN tol, long prec){
  pari_sp top=avma;
  GEN gpsu=psltopsu_mats(g, mats);//Image in PSU
  return gerepileupto(top, isometriccircle_psu(gpsu, tol, prec));
}

//Given an element g of PSU(1, 1), this returns the isometric circle associated to it.
GEN isometriccircle_psu(GEN g, GEN tol, long prec){
  pari_sp top=avma;
  if(toleq(gen_0, gcoeff(g, 2, 1), tol, prec)) return gen_0;//Isometric circle is everything, don't want to call it here.
  GEN geod=zerovec(8);
  gel(geod, 2)=gdivsg(1, gcoeff(g, 2, 1));//Need to take absolute value
  gel(geod, 1)=gneg(gmul(gcoeff(g, 2, 2), gel(geod, 2)));//-g[2,2]/g[2,1], the centre of the circle
  gel(geod, 2)=gabs(gel(geod, 2), prec);//We do things in this order to save a division.
  gel(geod, 7)=gen_1;//Always oriented counterclockwise
  pari_CATCH(CATCH_ALL){
    avma=top;
    pari_CATCH_reset();
    return gen_m1;//We increase precision in g and retry.
  }
  pari_TRY{
    GEN ipts=circle_int(geod, mkvec3s(0, 1, 0), tol, prec);//Intersect with x^2+y^2=1
    GEN ang=anglediff(garg(gsub(gel(ipts, 2), gel(geod, 1)), prec), garg(gsub(gel(ipts, 1), gel(geod, 1)), prec), tol, prec);
    if(gcmp(ang, mppi(prec))==1){//Properly orienting the start and endpoints
      gel(geod, 3)=gel(ipts, 2);
    gel(geod, 4)=gel(ipts, 1);
    }
    else{
      gel(geod, 3)=gel(ipts, 1);
      gel(geod, 4)=gel(ipts, 2);
    }
    gel(geod, 5)=radialangle(geod, gel(geod, 3), tol, prec);//Start angle
    gel(geod, 6)=shiftangle(radialangle(geod, gel(geod, 4), tol, prec), gel(geod, 5), tol, prec);//End angle
    pari_CATCH_reset();
    return gerepilecopy(top, geod);
  }
  pari_ENDCATCH
}

//Returns the normalized basis of G. Follows Algorithm 4.7 of Voight. Can pass in Ubase as a normalized boundary to append to, or Ubase=0 means we just start with G.
GEN normalizedbasis(GEN G, GEN Ubase, GEN mats, GEN gamid, GEN *data, GEN (*gamtopsl)(GEN *, GEN, long), GEN (*eltmul)(GEN *, GEN, GEN), GEN (*eltinv)(GEN *, GEN), int (*istriv)(GEN *, GEN), GEN tol, long prec){
  pari_sp top=avma, mid;
  long w=lg(G);
  GEN Gwithinv=cgetg(2*w-1, t_VEC);
  for(long i=1;i<lg(G);i++){
    gel(Gwithinv, i)=gel(G, i);
    gel(Gwithinv, w)=eltinv(data, gel(G, i));
    w++;
  }
  GEN U;
  if(gequal0(Ubase)) U=normalizedboundary(Gwithinv, mats, gamid, data, gamtopsl, tol, prec);//Step 2 when Ubase=0
  else U=normalizedboundary_givenU(Ubase, Gwithinv, mats, gamid, data, gamtopsl, tol, prec);//Step 2
  if(lg(gel(U, 1))==1) return gerepileupto(top, U);//No iso circles, go back now.
  GEN Gadd=Gwithinv, gbardat, gbar, unpair, v, g, scale=gdivgs(stoi(9), 10), Gaddnew, Uold=gel(U, 1);//Uold will track if we need to go back to step 3 (from step 4) or increase the scale (in step 5)
  long lunp, gind, vind;
  for(;;){
    if(gc_needed(top, 2)){//Garbage collection
      mid=avma;
      Uold=gcopy(Uold);
      U=gcopy(U);
      scale=gcopy(scale);
      Gadd=gcopy(Gadd);
      gerepileallsp(top, mid, 4, &Uold, &U, &scale, &Gadd);
    }
    Gaddnew=vectrunc_init(2*lg(Gadd));//The reductions to add to G
    for(long i=1;i<lg(Gadd);i++){//Doing step 3
      gbardat=reduceelt_givennormbound(U, gel(Gadd, i), gen_0, gamid, data, gamtopsl, eltmul, tol, prec);//Finding gbar=red_U(g). We actually only need to do this for g^(-1) not in U[1]; maybe add this optimization later.
      gbar=gel(gbardat, 1);
      if(!istriv(data, gbar)){
        vectrunc_append(Gaddnew, gbar);//If not trivial, add it.
        vectrunc_append(Gaddnew, eltinv(data, gbar));//Need to add the inverse as well.
      }
    }
    if(lg(Gaddnew)!=1){//We add Gaddnew to G, compute the normalized boundary, and go back to step 3 if U is changed
      Gadd=Gaddnew;
      U=normalizedboundary_givenU(U, Gadd, mats, gamid, data, gamtopsl, tol, prec);//Adding Gaddnew to U, recomputing normalized boundary
      if(lg(Uold)==lg(gel(U, 1))){
        if(!gequal(Uold, gel(U, 1))){Uold=gel(U, 1);continue;}//Go back to step 3.
      }
    }
    //Step 5. Find unpaired vertices, add them, go back!
    unpair=edgepairing(U, tol, 0, prec);//We only want the unpaired edges if they exist, and the paired only if there are no unpaired
    lunp=lg(unpair);
    if(typ(unpair)==t_VECSMALL){gel(U, 7)=unpair;break;}//Done!
    Gaddnew=vectrunc_init(2*lunp-1);//We add up to 2 new elements per unpaired vertex (there are lunp-1 such)
    for(long i=1;i<lunp;i++){
      gind=gel(unpair, i)[1];
      for(long k=2;k<lg(gel(unpair, i));k++){//The 1 or 2 vertices
        vind=gel(unpair, i)[k];
        v=gel(gel(U, 3), vind);
        if(toleq(gabs(v, prec), gen_1, tol, prec)){//Infinite vertex
          if(gind==vind) v=normalizedbasis_shiftpoint(gel(gel(U, 2), gind), scale, 1, prec);
          else v=normalizedbasis_shiftpoint(gel(gel(U, 2), gind), scale, 0, prec);
        }
        g=gel(gel(U, 1), gind);
        gbardat=reduceelt_givennormbound(U, g, v, gamid, data, gamtopsl, eltmul, tol, prec);//Reduce with respect to v (and not 0 like in step 3).
        vectrunc_append(Gaddnew, gel(gbardat, 1));
      } 
    }
    Gadd=Gaddnew;
    U=normalizedboundary_givenU(U, Gadd, mats, gamid, data, gamtopsl, tol, prec);//Adding Gadd to U, recomputing normalized boundary.
    if(lg(gel(U, 1))==lg(Uold)){//They might be equal
      if(gequal(gel(U, 1), Uold)){//They are equal. There must be an error where we have an oo vertex which gets shifted to being in the INTERIOR, i.e. we scaled down too much. To fix, we increase scale
        scale=gdivgs(gaddgs(scale, 9), 10);
      }
    }
    Uold=gel(U, 1);
  }
  return gerepilecopy(top, U);
}

//Given a circle arc c and either the initial(=1) or terminal(=0) point, this returns a nearby point on the arc (the new angle is r of the way between them; want r to be close to 1).
static GEN normalizedbasis_shiftpoint(GEN c, GEN r, int initial, long prec){
  pari_sp top=avma;
  GEN newang;
  if(initial==0) newang=gadd(gel(c, 5), gmul(gsub(gel(c, 6), gel(c, 5)), r));//a+(b-a)*r with a the initial and b the terminal angle
  else newang=gsub(gel(c, 6), gmul(gsub(gel(c, 6), gel(c, 5)), r));//b-(b-a)*r with a the initial and b the terminal angle
  return gerepileupto(top, gadd(gel(c, 1), gmul(expIr(gtofp(newang, prec)), gel(c, 2))));//c[1]+c[2]*e^(I*newang), the shifted point.
}

//Given a normalized boundary U, this appends the elements of G to it (G is the set of [elt, image in PSU(1, 1), isometric circle].
GEN normalizedboundary_append(GEN Ubase, GEN G, GEN mats, GEN id, GEN tol, long prec){
  pari_sp top=avma, mid;
  long nGp1=lg(G), nG=nGp1-1;//nG=number of elements in G
  long nUp1=lg(gel(Ubase, 1)), nU=nUp1-1;//nUp1=number of elts in Ubase+1
  if(nU==0) return normalizedboundary_givencircles(G, mats, id, tol, prec);//Ubase is trivial.
  long maxsides=2*nG+nUp1+4;//maxsides is 5+the maximal number of sides we can have (every side exists and there are two sides added for each element of G (the side and an infinite side). The 5 is for added security (we loop around)
  GEN U=cgetg(maxsides, t_VECSMALL);//Stores the indices in U in order. The initial ordering may not have v_1 correct (may need to do a cyclic shift at the end).
  GEN vertices=cgetg(maxsides, t_VEC);//Stores the vertices in order; an entry is [vertex, radial angle]
  GEN pi=mppi(prec);//Pi
  GEN inter, ang, ang1, ang2, sidecirc, sidecirctermang, Ltermang, ten=stoi(10);
  GEN L=gel(gel(Ubase, 2), 1);//The current segment we are looking for intersections with.
  GEN baseang=garg(gel(L, 4), prec);//Angle to the terminal point of the first element of Ubase
  GEN Utermangles=cgetg(nUp1, t_VEC);//The angles to the terminal points in Ubase
  for(long i=1;i<nUp1;i++){
    L=gel(gel(Ubase, 2), i);
    if(gequal0(L)) gel(Utermangles, i)=gen_0;//infinite side, ignoring
    else gel(Utermangles, i)=shiftangle(garg(gel(L, 4), prec), baseang, tol, prec);//Angle to origin in [baseind, baseind+2*Pi)
  }
  GEN Gtermangles=cgetg(nGp1, t_VEC);
  for(long i=1;i<nGp1;i++){
    if(gequal0(gel(gel(G, i), 3))) gel(Gtermangles, i)=gen_m2;//Does not give rise to a circle, want to ignore. -2 will ALWAYS be less than the start angle.
    else gel(Gtermangles, i)=shiftangle(garg(gel(gel(gel(G, i), 3), 4), prec), baseang, tol, prec);//Angle to origin in [baseind, baseind+2*Pi)
  }
  GEN Gord=indexsort(Gtermangles);//The order in which we want to look at the elements of G.
  long Gordind=0;
  for(long i=1;i<nGp1;i++){
    if(!gequal(gel(Gtermangles, Gord[i]), gen_m2)){Gordind=i;break;}
  }//Moving past the -2's, i.e. elements of G giving no circle. These occur first as the other angles are >-2.
  if(Gordind==0) return gerepilecopy(top, Ubase);//No new circles.
  U[1]=1;
  gel(vertices, 1)=mkvec2(gel(gel(gel(Ubase, 2), 1), 4), baseang);//The first vertex is initially set to be the start angle of the first side.
  long ulen=1;//The current length of U and vertices. If we have to delete some vertices, this can decrease.
  long side, sidem1;//We basically re-insert the first side back at the end
  GEN gang=gel(Gtermangles, Gord[Gordind]);//The angle to the terminal side of the current element of G.
  int lastsidenew=0, newsidefromG, finalstretch=0, leftoverGs=0;
  long startpt=1;//Stores what will be the start point.
  for(long sid=2;sid<=2*nU;sid++){//We move through Ubase, stopping to insert the elements of G if need be.
    if(sid>nU){
      if(!gequal(gang, ten)) leftoverGs=1;//There are still G's to do, and we have looped around
      else leftoverGs=0;
      side=sid-nU;
      if(side==1) sidem1=nU;
      else sidem1=side-1;
      finalstretch=1;//We stop once the (looped around) side gets re-inserted in.
    }
    else{
      side=sid;
      sidem1=sid-1;
    }
    mid=avma;
	if(finalstretch && !leftoverGs){//We go back and hit the start of U
	  if(U[side]==0){startpt++;continue;}
	  else if(U[side]>0){
		sidecirc=gel(gel(Ubase, 2), U[side]);
		sidecirctermang=gel(Utermangles, U[side]);
		newsidefromG=0;
	  }
	  else{
		sidecirc=gel(gel(G, -U[side]), 3);
		sidecirctermang=gel(Gtermangles, -U[side]);
		newsidefromG=1;
	  }
	  if(U[ulen]>0){
		L=gel(gel(Ubase, 2), U[ulen]);
		Ltermang=gel(Utermangles, U[ulen]);
	  }
	  else{
		L=gel(gel(G, -U[ulen]), 3);
		Ltermang=gel(Gtermangles, -U[ulen]);
	  }
	}
	else if(gequal0(gel(gel(Ubase, 2), side))) continue;//We skip past infinite sides of Ubase.
    else if(lastsidenew==0){//Working on a non-infinite side of Ubase, and the last side was also a side of Ubase
      if(tolcmp(gang, gel(Utermangles, side), tol, prec)==1 && !leftoverGs){//Consecutive old sides
	    if(gequal0(gel(gel(Ubase, 2), sidem1))){//The last side was infinite and skipped over; must re-insert it.
          ulen++;
          U[ulen]=0;
          gel(vertices, ulen)=mkvec2(gel(gel(Ubase, 3), sidem1-1), gel(gel(Ubase, 4), sidem1-1));//The infinite side
		  ulen++;
          U[ulen]=side;
          gel(vertices, ulen)=mkvec2(gel(gel(Ubase, 3), sidem1), gel(gel(Ubase, 4), sidem1));
          continue;//Go on
        }
		//Now we check if the new vertex is beyond the old one
	    ang=gel(gel(Ubase, 4), sidem1);//Angle to the intersection point
		ang1=anglediff(ang, gel(gel(vertices, ulen), 2), tol, prec);//ang-angle to the previous vertex.
        if(gequal0(ang1) || gcmp(ang1, pi)>=0){//Delete last side and go backwards. The previous side MUST be a new side.
		  avma=mid;
		  ulen--;
		  lastsidenew=1;
		  sid--;
		  continue;
		}
        ulen++;
        U[ulen]=side;
        gel(vertices, ulen)=mkvec2(gel(gel(Ubase, 3), sidem1), gel(gel(Ubase, 4), sidem1));
        continue;//Go on
      }
      //Here, this means that we need to try to insert G[Gord[Gordind]] BEFORE the current side.
      sidecirc=gel(gel(G, Gord[Gordind]), 3);//The new side we are trying to insert
      sidecirctermang=gang;//Terminal angle
      L=gel(gel(Ubase, 2), sidem1);//The previous side
      if(gequal0(L)){L=gel(gel(Ubase, 2), sidem1-1);Ltermang=gel(Utermangles, sidem1-1);}//If the previous side was oo, we go back one further, as this side is not oo.
      else Ltermang=gel(Utermangles, sidem1);
      newsidefromG=1;
    }
    else{//Now we check if we are inserting a new side or an old side
      if(tolcmp(gang, gel(Utermangles, side), tol, prec)==1 && !leftoverGs){//Inserting old side
        sidecirc=gel(gel(Ubase, 2), side);//The old side we are trying to insert
        sidecirctermang=gel(Utermangles, side);//Terminal angle
        L=gel(gel(G, -U[ulen]), 3);//The newly inserted side
        Ltermang=gel(Gtermangles, -U[ulen]);
        newsidefromG=0;
      }
      else{
        sidecirc=gel(gel(G, Gord[Gordind]), 3);//The new side we are trying to insert
        sidecirctermang=gang;//Terminal angle
        L=gel(gel(G, -U[ulen]), 3);//The newly inserted side
        Ltermang=gel(Gtermangles, -U[ulen]);
        newsidefromG=1;
      }
    }
    
    //At this point, we are either trying to insert a new on an old, or an old on a new side, or a new on a new
    inter=arc_int(L, sidecirc, tol, prec);//Intersection of L and the next side.
    if(lg(inter)==1){//It did NOT intersect.
      ang2=garg(gel(L, 3), prec);
      ang1=anglediff(ang2, Ltermang, tol, prec);//The angle between the terminal and initial points of L
      ang=anglediff(garg(gel(sidecirc, 3), prec), Ltermang, tol, prec);//Angle to the initial point of sidecirc from the terminal angle of L.
      if(tolcmp(ang, ang1, tol, prec)<=0){//sidecirc is contained within L.
        avma=mid;
		if(finalstretch && !leftoverGs) startpt++;//We add 1 to startpt to signal that we must start at a different point.
        else if(newsidefromG){//Failed to insert since it did not help.
          Gordind++;
          if(Gordind==nGp1){gang=ten;}//We are done with G, so we just want to start appending old sides.
          else gang=gel(Gtermangles, Gord[Gordind]);
          sid--;//We need to try again with the current side since we "jumped the line" with the element of G.
        }//We also want to leave lastsidenew unchanged, as we did not insert
        continue;
      }
      //We have two new sides: a side at infinity, and this side.
      ulen++;
      U[ulen]=0;
      gel(vertices, ulen)=mkvec2(gel(L, 3), ang2);
      ulen++;
	  gel(vertices, ulen)=mkvec2(gel(sidecirc, 4), sidecirctermang);
	  if(finalstretch && !leftoverGs){
	    U[ulen]=U[side];
		break;//Done!
	  }
      else if(newsidefromG){
        U[ulen]=-Gord[Gordind];
        Gordind++;
        if(Gordind==nGp1){gang=ten;}//We are done with G, so we just want to start appending old sides.
        else gang=gel(Gtermangles, Gord[Gordind]);
        sid--;//Redo this side
        lastsidenew=1;
      }
      else{
        U[ulen]=side;
        lastsidenew=0;
      }
    }
    else{//It DID intersect
      //It may be that this new side actually comes into U[1] from the bottom. Then we need a side at oo
      ang=anglediff(sidecirctermang, Ltermang, tol, prec);//Angle from the terminal angle of the last side to the terminal angle of the new side.
      if(gequal0(ang)){
        if(tolcmp(gel(sidecirc, 2), gel(L, 2), tol, prec)<=0){//This side lies inside the previous one, continue on (compared radii).
          avma=mid;
		  if(finalstretch && !leftoverGs) startpt++;
          else if(newsidefromG){//Failed to insert since it did not help.
            Gordind++;
            if(Gordind==nGp1){gang=ten;}//We are done with G, so we just want to start appending old sides.
            else gang=gel(Gtermangles, Gord[Gordind]);
            sid--;//We need to try again with the current side since we "jumped the line" with the element of G.
          }//We also want to leave lastsidenew unchanged, as we did not insert
          continue;
        }
		//Only need to update U[ulen]; the intersection point is the same.
		if(newsidefromG){U[ulen]=-Gord[Gordind];lastsidenew=1;}
		else{U[ulen]=side;lastsidenew=0;}
		continue;
      }
      else if(gcmp(ang, pi)==1){//We DID come in from below
        ulen++;
        U[ulen]=0;//Side at oo
        gel(vertices, ulen)=mkvec2(gel(L, 3), garg(gel(L, 3), prec));//Side at oo
        ulen++;
		gel(vertices, ulen)=mkvec2(gel(sidecirc, 4), sidecirctermang);
		if(finalstretch && !leftoverGs){//Done!
		  U[ulen]=U[side];
		  break;
		}
        if(newsidefromG){
          U[ulen]=-Gord[Gordind];
          Gordind++;
          if(Gordind==nGp1){gang=ten;}//We are done with G, so we just want to start appending old sides.
          else gang=gel(Gtermangles, Gord[Gordind]);
          sid--;//Redo this side
          lastsidenew=1;
        }
        else{
          U[ulen]=side;
          lastsidenew=0;
        }
        continue;
      }
      //Now we are sure did not come in from below.
      inter=gel(inter, 1);//The point
      if(toleq(inter, gel(L, 3), tol, prec)){//The side lies entirely in the previous side OR touches it at the end
	    if(toleq(inter, gel(sidecirc, 3), tol, prec)){//Lies inside
          avma=mid;
		  if(finalstretch && !leftoverGs) startpt++;
          else if(newsidefromG){//Failed to insert since it did not help.
            Gordind++;
            if(Gordind==nGp1){gang=ten;}//We are done with G, so we just want to start appending old sides.
            else gang=gel(Gtermangles, Gord[Gordind]);
            sid--;//We need to try again with the current side since we "jumped the line" with the element of G.
          }//We also want to leave lastsidenew unchanged, as we did not insert
          continue;
		}
      }
      //Now we have a proper "normal" intersection
      ang1=garg(inter, prec);
      ang=anglediff(ang1, gel(gel(vertices, ulen), 2), tol, prec);//Angle to the new vtx from the previous as a bases
      if(gcmp(ang, pi)!=-1 || toleq(ang, gen_0, tol, prec)){//We must go backwards!
        while(ulen>1){
          ulen--;
          if(U[ulen]>0) L=gel(gel(Ubase, 2), U[ulen]);
          else L=gel(gel(G, -U[ulen]), 3);//U[ulen]=0 is impossible; it is guarenteed that we do not backtrack past an infinite side in this algorithm.
          inter=gel(arc_int(L, sidecirc, tol, prec), 1);//They MUST intersect
          ang1=garg(inter, prec);
          ang=anglediff(ang1, gel(gel(vertices, ulen), 2), tol, prec);
          if(gcmp(ang, pi)!=-1 || toleq(ang, gen_0, tol, prec)) continue;//Keep going back
          break;//At this point we have reached where we need to insert the new side.
        }
      }
      //Now we are ready to insert it.
	  if(finalstretch && !leftoverGs){//Done if we intersected before the original intersection, and must continue on if not.
		if(gcmp(pi, anglediff(ang1, gel(gel(vertices, side+1), 2), tol, prec))==1){startpt++;continue;}//We have superseeded the previous side, increase start point.
		ulen++;
		U[ulen]=U[side];
	    gel(vertices, ulen)=mkvec2(inter, ang1);
		break;
	  }
	  ulen++;
	  gel(vertices, ulen)=mkvec2(inter, ang1);
      if(newsidefromG){
        U[ulen]=-Gord[Gordind];
        Gordind++;
        if(Gordind==nGp1){gang=ten;}//We are done with G, so we just want to start appending old sides.
        else gang=gel(Gtermangles, Gord[Gordind]);
        sid--;//Redo this side
        lastsidenew=1;
      }
      else{
        U[ulen]=side;
        lastsidenew=0;
      }
    }
  }
  //We potentially adjust the start position.
  long best=ulen;//How many sides we have (plus 1).
  mid=avma;
  long k=ulen-1;
  GEN L0=zerovec(8);gel(L0, 4)=gen_1;gel(L0, 7)=gen_1;gel(L0, 8)=gen_1;//L0=line segment [0, 1].
  GEN mininter=gen_2;
  ang1=anglediff(gel(Utermangles, 1), gen_0, tol, prec);
  while(U[k]<=0){//U[1] in the start position can only be superseeded by a new edge.
    if(U[k]==0){k--;continue;}//We may have added in a 0 edge
    L=gel(gel(G, -U[k]), 3);
    inter=arcseg_int(L, L0, tol, prec);
    if(lg(inter)!=1){//Intersection
      inter=real_i(gel(inter, 1));//Intersected [0, 1]
      if(tolcmp(inter, mininter, tol, prec)==-1){best=k;mininter=inter;}
	  else break;//We cannot do better than this.
    }
    else{//No intersection
      if(gequal(mininter, gen_2)) best=k;//Might have to take the side with smallest angle
      else break;//Done, as we had intersections, and now there are none.
      ang=anglediff(gel(Gtermangles, -U[k]), gen_0, tol, prec);
      if(tolcmp(ang, ang1, tol, prec)!=-1){best++;break;}//Done, and the best was actually one ago.
    }
    k--;
  }
  if(U[startpt]>0) L=gel(gel(Ubase, 2), U[startpt]);
  else L=gel(gel(G, -U[startpt]), 3);
  inter=arcseg_int(L, L0, tol, prec);
  if(lg(inter)>1){//Intersect. Now we need to go forward, as it is possible that the best intersection with [0, 1] was added right after the first side.
    inter=real_i(gel(inter, 1));
    if(tolcmp(mininter, inter, tol, prec)!=-1){best=ulen;mininter=inter;}//The first arc was better
	k=2;
	while(U[k]<0){
	  L=gel(gel(G, -U[k]), 3);
      inter=arcseg_int(L, L0, tol, prec);
	  if(lg(inter)==1) break;//No, intersection, done
      inter=real_i(gel(inter, 1));//Intersected [0, 1]
      if(tolcmp(inter, mininter, tol, prec)!=1){best=k;mininter=inter;}
	  else break;//We cannot do better than this.
	  k++;
	}
  }
  else{
	if(gequal(mininter, gen_2) && U[k]<0) best++;//There was NO intersection with [0, 1], so we are ending up on an infinite side! This is not right, so we must increment it by one. If best=ulen then we did not boop it backward, so don't need to increment if forward
  }
  
  avma=mid;
  long np1=ulen-startpt+1;
  GEN firstang=gel(gel(vertices, startpt+1), 2);//The angle to the first vertex
  //By wrapping back around, we have ulen-startpt sides: the last side is the same as the first.
  GEN ret=cgetg(NORMBOUND, t_VEC);
  for(long i=1;i<=5;i++) gel(ret, i)=cgetg(np1, t_VEC);//elements, icircs, vertices, matrices, term angles. Places 6, 7 will store the area and side pairing (0 for now)
  long j=1, h;
  for(long i=best;i<ulen;i++){
    if(U[i]==0){//Side at oo
      gel(gel(ret, 1), j)=gcopy(id);//Element
      gel(gel(ret, 2), j)=gen_0;//No circle
      gel(gel(ret, 3), j)=gcopy(gel(gel(vertices, i+1), 1));//Vertex
      gel(gel(ret, 4), j)=shiftangle(gel(gel(vertices, i+1), 2), firstang, tol, prec);//Vertex angle
      gel(gel(ret, 5), j)=gen_0;//No matrix
    }//Now we have a real side
    else if(U[i]<0){
      h=-U[i];
      gel(gel(ret, 1), j)=gcopy(gel(gel(G, h), 1));//Element
      gel(gel(ret, 2), j)=gcopy(gel(gel(G, h), 3));//Circle
      gel(gel(ret, 3), j)=gcopy(gel(gel(vertices, i+1), 1));//Vertex
      gel(gel(ret, 4), j)=shiftangle(gel(gel(vertices, i+1), 2), firstang, tol, prec);//Vertex angle
      gel(gel(ret, 5), j)=gcopy(gel(gel(G, h), 2));//Matrix
    }
    else{
      gel(gel(ret, 1), j)=gcopy(gel(gel(Ubase, 1), U[i]));//Element
      gel(gel(ret, 2), j)=gcopy(gel(gel(Ubase, 2), U[i]));//Circle
      gel(gel(ret, 3), j)=gcopy(gel(gel(vertices, i+1), 1));//Vertex
      gel(gel(ret, 4), j)=shiftangle(gel(gel(vertices, i+1), 2), firstang, tol, prec);//Vertex angle
      gel(gel(ret, 5), j)=gcopy(gel(gel(Ubase, 5), U[i]));//Matrix
    }
    j++;
  }
  for(long i=startpt;i<best;i++){
    if(U[i]==0){//Side at oo
      gel(gel(ret, 1), j)=gcopy(id);//Element
      gel(gel(ret, 2), j)=gen_0;//No circle
      gel(gel(ret, 3), j)=gcopy(gel(gel(vertices, i+1), 1));//Vertex
      gel(gel(ret, 4), j)=shiftangle(gel(gel(vertices, i+1), 2), firstang, tol, prec);//Vertex angle
      gel(gel(ret, 5), j)=gen_0;//No matrix
    }//Now we have a real side
    else if(U[i]<0){
      h=-U[i];
      gel(gel(ret, 1), j)=gcopy(gel(gel(G, h), 1));//Element
      gel(gel(ret, 2), j)=gcopy(gel(gel(G, h), 3));//Circle
      gel(gel(ret, 3), j)=gcopy(gel(gel(vertices, i+1), 1));//Vertex
      gel(gel(ret, 4), j)=shiftangle(gel(gel(vertices, i+1), 2), firstang, tol, prec);//Vertex angle
      gel(gel(ret, 5), j)=gcopy(gel(gel(G, h), 2));//Matrix
    }
    else{
      gel(gel(ret, 1), j)=gcopy(gel(gel(Ubase, 1), U[i]));//Element
      gel(gel(ret, 2), j)=gcopy(gel(gel(Ubase, 2), U[i]));//Circle
      gel(gel(ret, 3), j)=gcopy(gel(gel(vertices, i+1), 1));//Vertex
      gel(gel(ret, 4), j)=shiftangle(gel(gel(vertices, i+1), 2), firstang, tol, prec);//Vertex angle
      gel(gel(ret, 5), j)=gcopy(gel(gel(Ubase, 5), U[i]));//Matrix
    }
    j++;
  }
  gel(ret, 6)=hpolygon_area(gel(ret, 2), gel(ret, 3), tol, prec);//Area
  gel(ret, 7)=gen_0;
  gel(ret, 8)=gcopy(mats);
  return gerepileupto(top, ret);
}

//Initializes the inputs for normalizedboundary_append. Works BEST if p is given as an EXACT number.
GEN normalizedboundary_givenU(GEN Ubase, GEN G, GEN mats, GEN id, GEN *data, GEN (*gamtopsl)(GEN *, GEN, long), GEN tol, long prec){
  pari_sp top=avma;
  long lx;
  GEN Gnew=cgetg_copy(G, &lx);
  int skipped=0;
  for(long i=1;i<lx;i++){
    pari_CATCH(CATCH_ALL){
      gel(Gnew, i)=mkvec3(gen_0, gen_0, gen_0);//The third entry being 0 means this gets ignored
      skipped++;
    }
    pari_TRY{
      gel(Gnew, i)=isometriccircle_mats(gel(G, i), mats, data, gamtopsl, tol, prec);
    }
    pari_ENDCATCH
  }
  if(skipped>0){
    char *warningmsg=pari_sprintf("%d isometric circles were skipped due to insufficient precision", skipped);
    pari_warn(warner, warningmsg);
    pari_free(warningmsg);
  }
  return gerepileupto(top, normalizedboundary_append(Ubase, Gnew, mats, id, tol, prec));
}

//G is the set of [elt, image in PSU(1, 1), isometric circle]. This returns the normalized boundary of the exterior domain. The output is [elements, icircs, vertices, vertex angles, matrices, area, 0, mats]. The circle corresponding to elements[i] is icircs[i], and the vertices are vertices[i-1] and vertices[i]. matrices[i] is the image in PSU(1,1) of elements[i]. The element 1 corresponds to a section on the unit circle, which also corresponds to an angle and circle of -1. Vertex angles stores the radial angle to the ith vertex (with base angle being the first one). The area is the area, and the 0 stores the side pairing when we have a fundamental domain (so a priori stores nothing).
GEN normalizedboundary_givencircles(GEN G, GEN mats, GEN id, GEN tol, long prec){
  pari_sp top=avma, mid;
  long np1=lg(G), n=np1-1, twonp3=2*n+3;
  GEN U=cgetg(twonp3, t_VECSMALL);//Stores the indices in U in order. The initial ordering may not have v_1 correct (may need to do a cyclic shift at the end). We double since there are at most n sides coming from G and n sides coming from the unit circle.
  GEN vertices=cgetg(twonp3, t_VEC);//Stores the vertices in order; an entry is [vertex, radial angle]
  GEN pi=mppi(prec);//Pi
  GEN inter, ang, sidecirc, ang1;
  //Finding the first element.
  mid=avma;
  GEN L0=zerovec(8);gel(L0, 4)=gen_1;gel(L0, 7)=gen_1;gel(L0, 8)=gen_1;//L0=line segment [0, 1].
  long hminind=0;
  int isnew;
  GEN Hmin=gen_1;
  for(long i=1;i<np1;i++){//We start by finding intersections with L.
    sidecirc=gel(gel(G, i), 3);
    if(gequal0(sidecirc)) continue;//Ignore elts of G giving no circle.
    inter=arcseg_int(sidecirc, L0, tol, prec);
    if(lg(inter)==1) continue;//No intersection
    inter=real_i(gel(inter, 1));//Take the real part since it is actually real, and need to get the type correct.
    isnew=tolcmp(inter, Hmin, tol, prec);
    if(isnew==-1){hminind=i;Hmin=inter;continue;}//New min
    else if(isnew==1) continue;//Not new
    //Now we = the min, so need to see which initial point is longer
    ang=anglediff(garg(gel(gel(gel(G, i), 3), 3), prec), garg(gel(gel(gel(G, hminind), 3), 3), prec), tol, prec);
    if(gcmp(ang, pi)==-1) hminind=i;//Better min.
  }
  avma=mid;//Don't need these calcs other than hminind
  GEN baseang;
  if(hminind>0) baseang=garg(gel(gel(gel(G, hminind), 3), 4), prec);//The base angle to the terminal point.
  else baseang=gen_0;//We start from the 0 angle instead.
  GEN termangles=cgetg(np1, t_VEC);
  for(long i=1;i<np1;i++){
    if(gequal0(gel(gel(G, i), 3))) gel(termangles, i)=gen_m2;//Does not give rise to a circle, want to ignore. -2 will ALWAYS be less than the start angle.
    else gel(termangles, i)=shiftangle(garg(gel(gel(gel(G, i), 3), 4), prec), baseang, tol, prec);//Angle to origin in [baseind, baseind+2*Pi)
  }
  //We now order G by the angle to the terminal points.
  GEN ordering=indexsort(termangles);//The order in which we want to look at the elements of G.
  long startind=0;
  for(long i=1;i<np1;i++){
    if(!gequal(gel(termangles, ordering[i]), gen_m2)){startind=i;break;}
  }//Moving past the -2's, i.e. elements of G giving no circle. These occur first as the other angles are >=0. If hind!=0, then ordering[startind]=hind.
  //Now we start at G[ordering[ind]]. For the first element, we ONLY need to store the index.
  if(startind==0){//NO valid isometric circles inputted.
    avma=top;
    GEN retempty=cgetg(NORMBOUND, t_VEC);
    for(long i=1;i<=5;i++) gel(retempty, i)=cgetg(1, t_VEC);//elements, icircs, vertices, matrices, area, sidepairing
    gel(retempty, 6)=mkoo();
    gel(retempty, 7)=gen_0;
    return retempty;
  }
  if(hminind>0 && ordering[startind]!=hminind){//The first side has multiple circles coming out of its terminal point. I'm pretty sure this happens if and only if (well, in the Shimura curve case) Q is unramified everywhere.
    for(long i=startind+1;i<np1;i++){
      if(ordering[i]==hminind){ordering[i]=ordering[startind];ordering[startind]=hminind;break;}//Fix it
    }
  }
  U[1]=ordering[startind];
  gel(vertices, 1)=mkvec2(gel(gel(gel(G, ordering[startind]), 3), 4), gel(termangles, ordering[startind]));//The first vertex is initially set to be the start angle of the first side.
  long ulen=1;//The current length of U and vertices. If we have to delete some vertices, this can decrease.
  GEN L=gel(gel(G, ordering[startind]), 3), ang2;//The current segment we are looking for intersections with.
  long side;//We basically re-insert the first side back at the end
  for(long sid=startind+1;sid<=np1;sid++){//Doing things for the current side.
    if(sid==np1) side=startind;
    else side=sid;
    mid=avma;
    sidecirc=gel(gel(G, ordering[side]), 3);//The circle arc of the next side
    inter=arc_int(L, sidecirc, tol, prec);//Intersection of L and the next side.
    if(lg(inter)==1){//It did NOT intersect.
      ang1=garg(gel(L, 3), prec);
      ang=anglediff(ang1, gel(termangles, ordering[side]), tol, prec);//Angle to the initial point of L from the terminal point of the new side.
      if(gcmp(ang, pi)==-1 && ordering[side]!=U[ulen]){//the new side is contained entirely in L, OR comes in from below (the last check guarentees that it isn't ths same side; applicable when the final normalized boundary has 1 iso circle.
        ang2=garg(gel(L, 4), prec);
        ang=anglediff(gel(termangles, ordering[side]), ang2, tol, prec);//Angle to the terminal point of the new side from the terminal point of L
        if(gcmp(ang, pi)==-1){
          avma=mid;//the new side is contained entirely in L, discard and continue on (may as well reset avma).
          continue;
        }
        //Now we are actually coming in from below, so can continue on as normal.
      }
      //We have two new sides: a side at infinity, and this side.
      ulen++;
      U[ulen]=-1;
      gel(vertices, ulen)=mkvec2(gel(L, 3), ang1);
      ulen++;
      U[ulen]=ordering[side];
      gel(vertices, ulen)=mkvec2(gel(gel(gel(G, ordering[side]), 3), 4), gel(termangles, ordering[side]));//The terminal point of sidecirc is a vertex.
      L=gel(gel(G, ordering[side]), 3);//Setting L
    }
    else{//It DID intersect
      //It may be that this new side actually comes into U[1] from the bottom. Then we need a side at oo
      ang=anglediff(gel(termangles, U[ulen]), gel(termangles, ordering[side]), tol, prec);//Angle from the terminal angle of the last side to the terminal angle of the new side.
      if(gequal0(ang)){
        if(tolcmp(gel(gel(gel(G, ordering[side]), 3), 2), gel(gel(gel(G, U[ulen]), 3), 2), tol, prec)<=0){avma=mid;continue;}//This side lies inside the previous one, continue on (compared radii).
      }
      else if(gcmp(ang, pi)==-1){//We DID come in from below
        ulen++;
        U[ulen]=-1;//Side at oo
        gel(vertices, ulen)=mkvec2(gel(L, 3), garg(gel(L, 3), prec));//Side at oo
        ulen++;
        U[ulen]=ordering[side];
        gel(vertices, ulen)=mkvec2(gel(gel(gel(G, ordering[side]), 3), 4), gel(termangles, ordering[side]));
        L=gel(gel(G, ordering[side]), 3);
        continue;
      }
      //Now we are sure did not come in from below.
      inter=gel(inter, 1);//The point
      if(toleq(inter, gel(gel(gel(G, ordering[side]), 3), 3), tol, prec)){avma=mid;continue;}//The side lies entirely in the previous side
      ang1=garg(inter, prec);
      ang=anglediff(ang1, gel(gel(vertices, ulen), 2), tol, prec);//Angle to the new vtx from the previous as a bases
      if(gcmp(ang, pi)!=-1 || toleq(ang, gen_0, tol, prec)){//We must go backwards!
        while(ulen>1){
          ulen--;
          L=gel(gel(G, U[ulen]), 3);//The old side we look at
          inter=gel(arc_int(L, sidecirc, tol, prec), 1);//They MUST intersect
          ang1=garg(inter, prec);
          ang=anglediff(ang1, gel(gel(vertices, ulen), 2), tol, prec);
          if(gcmp(ang, pi)!=-1 || toleq(ang, gen_0, tol, prec)) continue;//Keep going back
          break;//At this point we have reached where we need to insert the new side.
        }
      }
      //Now we are ready to insert it.
      ulen++;
      L=sidecirc;
      U[ulen]=ordering[side];
      gel(vertices, ulen)=mkvec2(inter, ang1);
    }
  }
  GEN firstang=gel(gel(vertices, 2), 2);//The angle to the first vertex
  //By wrapping back around, we have ulen-1 sides: the last side is the same as the first.
  GEN ret=cgetg(NORMBOUND, t_VEC);
  for(long i=1;i<=5;i++) gel(ret, i)=cgetg(ulen, t_VEC);//elements, icircs, vertices, matrices, term angles. Places 6, 7 will store the area and side pairing (0 for now)
  for(long i=1;i<ulen;i++){
    if(U[i]==-1){//Side at oo
      gel(gel(ret, 1), i)=gcopy(id);//Element
      gel(gel(ret, 2), i)=gen_0;//No circle
      gel(gel(ret, 3), i)=gcopy(gel(gel(vertices, i+1), 1));//Vertex
      gel(gel(ret, 4), i)=shiftangle(gel(gel(vertices, i+1), 2), firstang, tol, prec);//Vertex angle
      gel(gel(ret, 5), i)=gen_0;//No matrix
      continue;
    }//Now we have a real side
    gel(gel(ret, 1), i)=gcopy(gel(gel(G, U[i]), 1));//Element
    gel(gel(ret, 2), i)=gcopy(gel(gel(G, U[i]), 3));//Circle
    gel(gel(ret, 3), i)=gcopy(gel(gel(vertices, i+1), 1));//Vertex
    gel(gel(ret, 4), i)=shiftangle(gel(gel(vertices, i+1), 2), firstang, tol, prec);//Vertex angle
    gel(gel(ret, 5), i)=gcopy(gel(gel(G, U[i]), 2));//Matrix
  }
  gel(ret, 6)=hpolygon_area(gel(ret, 2), gel(ret, 3), tol, prec);//Area
  gel(ret, 7)=gen_0;
  gel(ret, 8)=gcopy(mats);
  return gerepileupto(top, ret);
}

//Initializes the inputs for normalizedboundary_givencircles. Works BEST if p is given as an EXACT number.
GEN normalizedboundary(GEN G, GEN mats, GEN id, GEN *data, GEN (*gamtopsl)(GEN *, GEN, long), GEN tol, long prec){
  pari_sp top=avma;
  long lx;
  GEN Gnew=cgetg_copy(G, &lx);
  int skipped=0;
  for(long i=1;i<lx;i++){
    pari_CATCH(CATCH_ALL){
      gel(Gnew, i)=mkvec3(gen_0, gen_0, gen_0);//The third entry being 0 means this gets ignored
      skipped++;
    }
    pari_TRY{
      gel(Gnew, i)=isometriccircle_mats(gel(G, i), mats, data, gamtopsl, tol, prec);
    }
    pari_ENDCATCH
  }
  if(skipped>0){
    char *warningmsg=pari_sprintf("%d isometric circles were skipped due to insufficient precision", skipped);
    pari_warn(warner, warningmsg);
    pari_free(warningmsg);
  }
  return gerepileupto(top, normalizedboundary_givencircles(Gnew, mats, id, tol, prec));
}

//This returns [v, vecsmall(ind)], where the first side of U that c intersects has index ind, and the point of intersection is v. If start=1, we are searching near the start point of c, else we seach near the end point of c (c is either an arc or a segment with start/end points on the unit circle. If we intersect an oo side, returns [v, vecsmall(-1)] instead.
GEN normalizedboundary_sideint(GEN U, GEN c, int start, GEN tol, long prec){
  pari_sp top=avma;
  GEN v, ret, inter, d1, d2;
  long ind;
  if(gequal(gel(c, 8), gen_1)){//Line segment; the index found from normalizedboundary_outside is correct guarenteed.
    if(start==1) v=gel(c, 3);
    else v=gel(c, 4);
    ind=normalizedboundary_outside(U, v, tol, prec);//This is the index
    if(ind==-1){//Infinite side
      ret=cgetg(3, t_VEC);
      gel(ret, 1)=gcopy(v);
      gel(ret, 2)=mkvecsmall(-1);
      return gerepileupto(top, ret); 
    }
    inter=arcseg_int(gel(gel(U, 2), ind), c, tol, prec);
    if(lg(inter)==2) v=gel(inter, 1);//1 intersection point, the most common
    else if(lg(inter)==3){
      d1=gabs(gsub(gel(inter, 1), v), prec);
      d2=gabs(gsub(gel(inter, 2), v), prec);
      if(gcmp(d1, d2)==-1) v=gel(inter, 1);//Finding which intersection is closer to v. Tolerance not required.
      else v=gel(inter, 2);
    }
    else pari_err_TYPE("Should have been intersections, but there were none", inter);//This should never happen, but here in case.
    ret=cgetg(3, t_VEC);
    gel(ret, 1)=gcopy(v);
    gel(ret, 2)=mkvecsmall(ind);
    return gerepileupto(top, ret);
  }
  GEN v2;//The other point
  if(start==1){
    if(gequal(gel(c, 7), gen_m1)){v=gel(c, 4);v2=gel(c, 3);}//If directed backwards, need the terminal point as the start
    else {v=gel(c, 3);v2=gel(c, 4);}//If undirected or directed forwards, the initial point is the start
  }
  else{
    if(gequal(gel(c, 7), gen_m1)){v=gel(c, 3);v2=gel(c, 4);}//If directed backwards, need the initial point.
    else{v=gel(c, 4);v2=gel(c, 3);}//Undirected/forwards, the terminal point is the start.
  }
  ind=normalizedboundary_outside(U, v, tol, prec);
  long lU=lg(gel(U, 1));//This is the index
  if(ind==-1){//Infinite side
    ret=cgetg(3, t_VEC);
    gel(ret, 1)=gcopy(v);
    gel(ret, 2)=mkvecsmall(-1);
    return gerepileupto(top, ret); 
  }
  //Circle arc (which should be the much more common case).
  GEN vang, ang, pi=mppi(prec);//Base angle+Pi
  int where, ind2;
  for(;;){//We intersect, see if it is on the sement or left or right of it. If left go to next index, if right go to previous
    inter=arc_int(gel(gel(U, 2), ind), c, tol, prec);
    if(lg(inter)==2) v=gel(inter, 1);//Only 1 intersection point, the most common case
    else if(lg(inter)==3){//Two intersection points
      d1=gabs(gsub(gel(inter, 1), v), prec);
      d2=gabs(gsub(gel(inter, 1), v2), prec);
      if(gcmp(d1, d2)==-1) v=gel(inter, 1);//Intersection[1] is closer to v than v2, so this is right.
      else v=gel(inter, 2);
    }
    else pari_err_TYPE("Should have been intersections, but there were none", inter);//This should never happen, but here in case.
    vang=garg(v, prec);//Angle to v
    ang=anglediff(gel(gel(U, 4), ind), vang, tol, prec);//Angle to the vertex from v.
    where=tolcmp(ang, pi, tol, prec);
    if(where==1){//Left of the vertex, continue
      ind++;
      if(ind==lU) ind=1;
      continue;
    }
    if(ind==1) ind2=lU-1;
    else ind2=ind-1;
    ang=anglediff(vang, gel(gel(U, 4), ind2), tol, prec);//Angle to v with reference to vertex 2
    where=tolcmp(pi, ang, tol, prec);//No need for cases this time
    if(where>=0) break;//Correct side!
    ind=ind2;//We must go backwards to the side ind2
  }
  ret=cgetg(3, t_VEC);
  gel(ret, 1)=gcopy(v);
  gel(ret, 2)=mkvecsmall(ind);
  return gerepileupto(top, ret);
}

//Returns -1 if z is in the interior of the normalized boundary or on the edge, and ind if in the exterior (i.e. between the boundary and the unit circle), where ind is the index of the side it is on (projecting from the origin to z). Does not behave well if z is a vertex of U that is on the unit disc.
long normalizedboundary_outside(GEN U, GEN z, GEN tol, long prec){
  pari_sp top=avma;
  int outside;
  pari_CATCH(CATCH_ALL){//Catching if U is trivial OR z=0
    avma=top;
    pari_CATCH_reset();
    return -1;
  }
  pari_TRY{
    GEN ang=shiftangle(garg(z, prec), gel(gel(U, 4), 1), tol, prec);//Shifting to base angle
    long ind=gen_search(gel(U, 4), ang, 1, NULL, &gcmp_strict);//Index to put z. We ONLY need to search for this cicle.
    if(ind==lg(gel(U, 1))) ind=1;//Insert at the end means the first circle.
    GEN circle=gel(gel(U, 2), ind);
    if(gequal0(circle)){pari_CATCH_reset();avma=top;return -1;}//Intersects with the edge of the unit disc.
    outside=tolcmp(gel(circle, 2), gabs(gsub(z, gel(circle, 1)), prec), tol, prec);//Are we outside?
    if(outside==0) outside=-1;
    else if(outside==1) outside=ind;
  }
  pari_ENDCATCH
  avma=top;
  return outside;//There is no tolerance issues with our search for ind; they are taken care of by the tolerance check with inside (the only possible issues occur if z and v[ind] or v[ind-1] are equal up to tolerance, but of course that is solved by the one tolcmp).
}

//normalizedboundary_inside with typechecking
long normalizedboundary_outside_tc(GEN U, GEN z, long prec){
  pari_sp top=avma;
  if(typ(U)!=t_VEC || lg(U)!=NORMBOUND) pari_err_TYPE("Please enter a normalized boundary", U);
  GEN tol=deftol(prec);
  int ins=normalizedboundary_outside(U, z, tol, prec);
  avma=top;
  return ins;
}

//Given g in PSL(2, R) and p in the upper half plane, this returns the image of g in PSU(1, 1) via phi(z)=(z-p)/(z-conj(p)).
GEN psltopsu(GEN g, GEN p){
  pari_sp top=avma;
  GEN M=psltopsu_transmats(p);
  return gerepileupto(top, psltopsu_mats(g, M));
}

//Given g in PSL(2, R), M=[m1, m2, p] with m1=1/(p-conj(p))[1,-p;1,-conj(p)], m2=[conj(p), -p;1, -1], this returns m1*g*m2.
GEN psltopsu_mats(GEN g, GEN M){
  pari_sp top=avma;
  return gerepileupto(top, gmul(gel(M, 1), gmul(g, gel(M, 2))));
}

//Returns [m1, m2, p]: matrices so that if g is in PSL(2, R), then the image of g in PSU(1, 1) corresponding to phi(z)=(z-p)/(z-conj(p)) is m1*g*m2.
GEN psltopsu_transmats(GEN p){
  pari_sp top=avma;
  GEN pneg=gneg(p), pconj=conj_i(p);
  GEN m1=gdiv(mkmat22(gen_1, pneg, gen_1, gneg(pconj)), gsub(p, pconj));//m1=1/(p-conj(p))[1,-p;1,-conj(p)]
  GEN m2=mkmat22(pconj, pneg, gen_1, gen_m1);//m2=[conj(p), -p;1, -1]
  GEN ret=cgetg(4, t_VEC);
  gel(ret, 1)=gcopy(m1);
  gel(ret, 2)=gcopy(m2);
  gel(ret, 3)=gcopy(p);
  return gerepileupto(top, ret);
}

//Returns the roots of the hyperbolic matrix M in PSL(2, R) in order.
GEN psl_roots(GEN M, GEN tol, long prec){
  pari_sp top = avma;
  GEN trace=gadd(gcoeff(M, 1, 1), gcoeff(M, 2, 2));
  int sgn=tolcmp(trace, gen_0, tol, prec);
  if(sgn==0) pari_err_TYPE("Please enter a hyperbolic matrix", M);
  GEN a, b, c;
  if(sgn==1){//Positive trace, correct order
    a=gcoeff(M, 2, 1);
    b=gsub(gcoeff(M, 2, 2), gcoeff(M, 1, 1));
    c=gneg(gcoeff(M, 1, 2));
  }
  else{
    a=gneg(gcoeff(M, 2, 1));
    b=gsub(gcoeff(M, 1, 1), gcoeff(M, 2, 2));
    c=gcoeff(M, 1, 2);
  }//[a',b';c',d'] -> c'x^2+(d'-a')x-b'=0, but for the roots to be in proper order, we need the trace to be positive.
  if(toleq(a, gen_0, tol, prec)){//a=0, roots are oo and -c/b (b!=0 else M=[+/-1, x;0;+/-1], not hyperbolic.
    GEN rnum=gneg(c);
    int bsgn=tolcmp(b, gen_0, tol, prec);
    GEN ret=cgetg(3, t_VEC);
    if(bsgn==1){//b>0, first root is finite
      gel(ret, 1)=gdiv(rnum, b);
      gel(ret, 2)=mkoo();
    }
    else if(bsgn==-1){//b<0, first root is oo
      gel(ret, 1)=mkoo();
      gel(ret, 2)=gdiv(rnum, b);
    }
    else pari_err_TYPE("Please enter a hyperbolic matrix", M);
    return gerepileupto(top, ret);
  }
  //Now both roots are finite.
  GEN twoa=gmulsg(2, a);//2a
  GEN rtD=gsqrt(gsub(gsqr(b), gmulsg(2, gmul(twoa, c))), prec);//b^2-4ac
  GEN mb=gneg(b);
  GEN r1num=gadd(mb, rtD);//first root is (-b+sqrt(D))/2a
  GEN r2num=gsub(mb, rtD);//Second root is (-b-sqrt(D))/2a
  GEN ret=cgetg(3, t_VEC);
  gel(ret, 1)=gdiv(r1num, twoa);
  gel(ret, 2)=gdiv(r2num, twoa);
  return gerepileupto(top, ret);
}

//Returns a random point z in the unit disc, uniform inside the ball of radius R. See page 19 of Page (before section 2.5).
GEN randompoint_ud(GEN R, long prec){
  pari_sp top=avma;
  GEN arg=gmul(randomr(prec), Pi2n(1, prec));//Random angle
  GEN zbound=expIr(arg);//The boundary point. Now we need to scale by a random hyperbolic distance in [0, R]
  //a(r)=Area of hyperbolic disc radius r=4*Pi*sinh^2(r/2).
  GEN dist=gmul(gsqr(gsinh(gdivgs(R, 2), prec)), randomr(prec));//A random element in [0, a(R)/4Pi].
  GEN r=gmulsg(2, gasinh(gsqrt(dist, prec), prec));//The radius
  GEN e2r=gexp(r, prec);
  return gerepileupto(top, gmul(zbound, gdiv(gsubgs(e2r, 1), gaddgs(e2r, 1))));
}

//Reduces g with respect to z as in reduceelt_givenpsu, but does so much more efficiently using the normalized boundary provided.
GEN reduceelt_givennormbound(GEN U, GEN g, GEN z, GEN gamid, GEN *data, GEN (*gamtopsl)(GEN *, GEN, long), GEN (*eltmul)(GEN *, GEN, GEN), GEN tol, long prec){
  pari_sp top=avma;
  GEN delta=gamid;
  llist *decomp=NULL;
  GEN gmat=psltopsu_mats(gamtopsl(data, g, prec), gel(U, 8));
  z=mat_eval(gmat, z);//gz, the real start point
  long count=0, outside;
  for(;;){
    outside=normalizedboundary_outside(U, z, tol, prec);
    if(outside==-1) break;//Done:Either reached inside or the boundary.
    z=mat_eval(gel(gel(U, 5), outside), z);//Update z
    delta=eltmul(data, gel(gel(U, 1), outside), delta);//update delta
    llist_putstart(&decomp, outside);//add outside to the list
    count++;
  }
  GEN ret=cgetg(4, t_VEC);
  gel(ret, 1)=eltmul(data, delta, g);//gbar=delta*g
  gel(ret, 2)=gcopy(delta);
  gel(ret, 3)=llist_tovecsmall(decomp, count, 1);
  return gerepileupto(top, ret);
}

//Algorithm 4.3 of Voight. Inputs G, a finite subset of Gamma, corresponding to Gmats in PSU(1, 1), g (->gmat) an element of Gamma, z in the unit disc. This G-reduces g, i.e. translating gz to the exterior domain of G. Returns [gbar, delta, decomp], where gbar=delta*g. gbar is (G, z)-reduced, delta is in <G>, and delta=G[i1]*G[i2]*...*G[in] with decomp=[i1, i2, ..., in] (vecsmall).
GEN reduceelt_givenpsu(GEN G, GEN Gmats, GEN g, GEN gmat, GEN z, GEN gamid, GEN *data, GEN (*eltmul)(GEN *, GEN, GEN), GEN tol, long prec){
  pari_sp top=avma;
  GEN delta=gamid;
  llist *decomp=NULL;
  z=mat_eval(gmat, z);//gz, the real start point
  GEN mindist, curdist=hdist_ud(z, gen_0, prec), znew, zmin, dist;
  long ind, n=lg(G), count=0;
  for(;;){
    zmin=mat_eval(gel(Gmats, 1), z);
    mindist=hdist_ud(zmin, gen_0, prec);
    if(typ(mindist)==t_COMPLEX) mindist=mkoo();//Rounding error forced point outside the unit disc.
    ind=1;
    for(long i=2;i<n;i++){
      if(gequal0(gel(Gmats, i))) continue;//oo side, ignore
      znew=mat_eval(gel(Gmats, i), z);
      dist=hdist_ud(znew, gen_0, prec);//Distance to g_i*z
      if(typ(dist)==t_COMPLEX) continue;//Rounding error forced us onto point is too close to being outside the unit disc.
      if(tolcmp(dist, mindist, tol, prec)==-1){//Strictly smaller distance
        zmin=znew;
        mindist=dist;
        ind=i;
      }
    }
    if(tolcmp(mindist, curdist, tol, prec)!=-1) break;//Done
    count++;
    z=zmin;
    curdist=mindist;
    delta=eltmul(data, gel(G, ind), delta);
    llist_putstart(&decomp, ind);
  }
  GEN ret=cgetg(4, t_VEC);
  gel(ret, 1)=eltmul(data, delta, g);//gbar=delta*g
  gel(ret, 2)=gcopy(delta);
  gel(ret, 3)=llist_tovecsmall(decomp, count, 1);
  return gerepileupto(top, ret);
}

//Reduces z to the interior of U (Almost identical to reduceelt_givennormbound). Returns [g, z'], where g is the transition element and z' is the new point.
GEN reducepoint(GEN U, GEN z, GEN gamid, GEN *data, GEN (*eltmul)(GEN *, GEN, GEN), GEN tol, long prec){
  pari_sp top=avma;
  GEN g=gamid;
  long outside;
  for(;;){
    outside=normalizedboundary_outside(U, z, tol, prec);
    if(outside==-1) break;//Done:Either reached inside or the boundary.
    z=mat_eval(gel(gel(U, 5), outside), z);//Update z
    g=eltmul(data, gel(gel(U, 1), outside), g);//update g
  }
  GEN ret=cgetg(3, t_VEC);
  gel(ret, 1)=gcopy(g);
  gel(ret, 2)=gcopy(z);
  return gerepileupto(top, ret);
}

//Finds the image of the root geodesic of g in the fundamental domain specified by U.
GEN rootgeodesic_fd(GEN U, GEN g, GEN gamid, GEN *data, GEN (*gamtopsl)(GEN *, GEN, long), GEN (*eltmul)(GEN *, GEN, GEN), GEN (*eltinv)(GEN *, GEN), GEN tol, long prec){
  pari_sp top=avma;
  GEN gpsl=gamtopsl(data, g, prec), mats=gel(U, 8);//The image in PSL
  GEN geod=rootgeodesic_ud(gpsl, mats, tol, prec);
  GEN z;
  if(gequal0(gel(geod, 8))) z=arc_midpoint(geod, gel(geod, 3), gel(geod, 4), tol, prec);//First move the midpoint to the fundamental domain;
  else z=gdivgs(gadd(gel(geod, 3), gel(geod, 4)), 2);//When geod is a line segment, do this instead (Q ram at 11, 13 and g=[11/2, 3/2, 0, 0] e.g.)
  GEN red=reducepoint(U, z, gamid, data, eltmul, tol, prec);
  g=eltmul(data, eltmul(data, gel(red, 1), g), eltinv(data, gel(red, 1)));//Conjugating g by gel(red, 1);
  gpsl=gamtopsl(data, g, prec);//The image in PSL of the new g
  geod=rootgeodesic_ud(gpsl, mats, tol, prec);//The new geodesic, which necessarily goes through z, and hence the interior.
  //Now we need to find the start vertex, and run along.
  GEN vbaseinfo=normalizedboundary_sideint(U, geod, 1, tol, prec);//Find the vertex nearest the start of the geodesic.
  GEN vbase=gel(vbaseinfo, 1);
  GEN vstart=vbase, vend, startcentre=gel(geod, 1), starttype=gel(geod, 8);
  glist *Gs=NULL;//Tracking the g's corresponding to the arcs
  glist *circs=NULL;//Tracking the circle arcs.
  llist *sides=NULL;//Tracking the sides hit.
  llist *othersides=NULL;//Tracking the sides we are leaving from.
  llist_putstart(&othersides, gel(vbaseinfo, 2)[1]);
  long count=0;
  for(;;){
    vend=normalizedboundary_sideint(U, geod, 0, tol, prec);
    glist_putstart(&Gs, g);
    if(gequal0(gel(geod, 8))){//Arc
      if(gequal(gel(geod, 7), gen_m1)) glist_putstart(&circs, arc_init(geod, gel(vend, 1), vstart, -1, prec));//geod travelling backwards
      else glist_putstart(&circs, arc_init(geod, vstart, gel(vend, 1), 1, prec));//geod travelling normally
    }
    else{//Segment
      gel(geod, 3)=vstart;//Start
      gel(geod, 4)=gel(vend, 1);//End
      glist_putstart(&circs, geod);
    }
    llist_putstart(&sides, gel(vend, 2)[1]);//Adding the side hit.
    count++;
    vstart=mat_eval(gel(gel(U, 5), gel(vend, 2)[1]), gel(vend, 1));//The new start vertex
    g=eltmul(data, eltmul(data, gel(gel(U, 1), gel(vend, 2)[1]), g), eltinv(data, gel(gel(U, 1), gel(vend, 2)[1])));//Conjugating g by the side hit
    gpsl=gamtopsl(data, g, prec);//The image in PSL of the new g
    geod=rootgeodesic_ud(gpsl, mats, tol, prec);//The new geodesic
    if(toleq(vbase, vstart, tol, prec)){
      if(gequal(starttype, gel(geod, 8)) && toleq(startcentre, gel(geod, 1), tol, prec)) break;//Done! We need this second check since an embedding CAN have a self-intersection on the boundary (e.g. [Q, order]=qa_init_2primes(2, 3), p=I/2, g=[3674890, -1623699, 463914, -1391742]). This is sufficient because the two arcs/segments share a point AND the centre/slope => unique.
    }
    llist_putstart(&othersides, gel(U, 7)[gel(vend, 2)[1]]);//Hit it with the side pairing.
  }
  GEN ret=cgetg(5, t_VEC);
  gel(ret, 1)=glist_togvec(Gs, count, -1);
  gel(ret, 2)=glist_togvec(circs, count, -1);
  gel(ret, 3)=llist_tovecsmall(sides, count, -1);
  gel(ret, 4)=llist_tovecsmall(othersides, count, -1);
  return gerepileupto(top, ret);
}

//Returns the root geodesic in the unit disc corresponding to M in PSL(2, R) and p the reference point to mapping the upper half plane to the unit disc (mats=[m1, m2, p], with m1 being the mapping).
GEN rootgeodesic_ud(GEN M, GEN mats, GEN tol, long prec){
  pari_sp top=avma;
  GEN geod=rootgeodesic_uhp(M, tol, prec);
  return gerepileupto(top, mobius(gel(mats, 1), geod, tol, prec));
}

//Returns the upper half plane root geodesic of the hyperbolic element M in PSL(2, R)
GEN rootgeodesic_uhp(GEN M, GEN tol, long prec){
  pari_sp top=avma;
  GEN rts=psl_roots(M, tol, prec), arc;
  if(typ(gel(rts, 1))==t_INFINITY){//First root infinite; we have a vertical segment from rts[2] to oo
    arc=cgetg(ARCLEN, t_VEC);
    gel(arc, 1)=mkoo();
    gel(arc, 2)=gcopy(gel(rts, 2));
    gel(arc, 3)=gcopy(gel(rts, 2));//Starts at finite (second) root
    gel(arc, 4)=mkoo();//Ends at oo
    gel(arc, 5)=gen_0;
    gel(arc, 6)=gen_1;//Vertically up
    gel(arc, 7)=gen_0;//No dir since endpt is oo
    gel(arc, 8)=gen_1;//Segment
    return gerepileupto(top, arc);
  }
  else if(typ(gel(rts, 2))==t_INFINITY){//Second root infinite; we have a vertical segment fr
    arc=cgetg(ARCLEN, t_VEC);
    gel(arc, 1)=mkoo();
    gel(arc, 2)=gcopy(gel(rts, 1));
    gel(arc, 3)=mkoo();//Starts at oo
    gel(arc, 4)=gcopy(gel(rts, 1));//Ends at finite (second) root
    gel(arc, 5)=gen_0;
    gel(arc, 6)=gen_m1;//Vertically down
    gel(arc, 7)=gen_0;//No dir since endpt is oo
    gel(arc, 8)=gen_1;//Segment
    return gerepileupto(top, arc);
  }
  //Now both roots are finite, so we have a nice circle.
  GEN centre=gdivgs(gadd(gel(rts, 1), gel(rts, 2)), 2);
  int firstbigger=gcmp(gel(rts, 1), gel(rts, 2));//=1 if the first root is bigger, and -1 if not (no tolerance needed).
  if(firstbigger==1){
    GEN radius=gsub(gel(rts, 1), centre);
    GEN c=mkvec3(centre, radius, gen_0);//The full circle
    arc=arc_init(c, gel(rts, 1), gel(rts, 2), -1, prec);//Arc points from second to first root.
  }
  else{
    GEN radius=gabs(gsub(gel(rts, 2), centre), prec);
    GEN c=mkvec3(centre, radius, gen_0);//The full circle
    arc=arc_init(c, gel(rts, 2), gel(rts, 1), 1, prec);//Arc points from second to first root.
  }
  return gerepileupto(top, arc);
}

//rootgeodesic_uhp with typechecking
GEN rootgeodesic_uhp_tc(GEN M, long prec){
  pari_sp top=avma;
  if(typ(M)!=t_MAT || lg(M)!=3 || lg(gel(M, 1))!=3) pari_err_TYPE("Please enter a hyperbolic 2x2 matrix in PSL(2, R)", M);
  GEN tol=deftol(prec);
  return gerepileupto(top, rootgeodesic_uhp(M, tol, prec));
}




//PRINTING TO PLOTVIEWER



//Writes the set of arcs in the given colour to filename, so that it is ready to be executed by python.
void python_printarcs(GEN arcs, char *filename, int view, char *extrainput, long prec){
  pari_sp top=avma;
  if(!pari_is_dir("fdoms")){//Checking the directory
    int s=system("mkdir -p fdoms");
    if(s==-1) pari_err(e_MISC, "ERROR CREATING DIRECTORY fdoms");
  }
  char *fullfile=pari_sprintf("fdoms/%s.dat", filename);
  FILE *f=fopen(fullfile, "w");
  pari_free(fullfile);//Now we have created the output file f.
  GEN arc, fact=gdiv(stoi(180), mppi(prec));//fact=180/Pi
  for(long i=1;i<lg(arcs);i++){
    arc=gel(arcs, i);
    if(gequal0(arc)) continue;//Not a circle
    if(gequal0(gel(arc, 8))){//Arc
      pari_fprintf(f, "0 %lf %lf %lf %lf %lf %d\n", rtodbl(gtofp(real_i(gel(arc, 1)), prec)), rtodbl(gtofp(imag_i(gel(arc, 1)), prec)), rtodbl(gtofp(gel(arc, 2), prec)), rtodbl(gtofp(gmul(gel(arc, 5), fact), prec)), rtodbl(gtofp(gmul(gel(arc, 6), fact), prec)), itos(gel(arc, 7)));
    }
    else{//Segment
      pari_fprintf(f, "1 %lf %lf %lf %lf\n", rtodbl(gtofp(real_i(gel(arc, 3)), prec)), rtodbl(gtofp(imag_i(gel(arc, 3)), prec)), rtodbl(gtofp(real_i(gel(arc, 4)), prec)), rtodbl(gtofp(imag_i(gel(arc, 4)), prec)));
    }
  }
  fclose(f);
  if(view==1){
    char *line;
    if(extrainput==NULL) line=pari_sprintf("%s", filename);
    else line=pari_sprintf("%s %s", extrainput, filename);
    python_plotviewer(line);
    pari_free(line);
  }
  avma=top;
}

//Launches the plotviewer with the given inputs.
void python_plotviewer(char *input){
  char *command;
  command=pari_sprintf("cmd.exe /C start py fdviewer.py %s", input);
  int s=system(command);
  if(s==-1) pari_err(e_MISC, "ERROR EXECUTING COMMAND");
  pari_free(command);
}

//Writes the fundamental domain corresponding to U
void python_printfdom(GEN U, char *filename, long prec){
  pari_sp top=avma;
  if(!pari_is_dir("fdoms")){//Checking the directory
    int s=system("mkdir -p fdoms");
    if(s==-1) pari_err(e_MISC, "ERROR CREATING DIRECTORY fdoms");
  }
  char *fullfile=pari_sprintf("fdoms/%s.dat", filename);
  FILE *f=fopen(fullfile, "w");
  pari_free(fullfile);//Now we have created the output file f.
  GEN pair=gel(U, 7);
  pari_fprintf(f, "%d", pair[1]);
  for(long i=2;i<lg(pair);i++) pari_fprintf(f, " %d", pair[i]);//Print side pairing.
  GEN vangs=gel(U, 4);
  pari_fprintf(f, "\n%lf", rtodbl(gtofp(gel(vangs, 1), prec)));//First angle
  for(long i=2;i<lg(pair);i++) pari_fprintf(f, " %lf", rtodbl(gtofp(gel(vangs, i), prec)));//Print angles.
  pari_fprintf(f, "\n");
  GEN arcs=gel(U, 2);
  GEN arc, fact=gdiv(stoi(180), mppi(prec));//fact=180/Pi
  for(long i=1;i<lg(arcs);i++){
    arc=gel(arcs, i);
    if(gequal0(arc)) continue;//Not a circle
    pari_fprintf(f, "%lf %lf %lf %lf %lf\n", rtodbl(gtofp(real_i(gel(arc, 1)), prec)), rtodbl(gtofp(imag_i(gel(arc, 1)), prec)), rtodbl(gtofp(gel(arc, 2), prec)), rtodbl(gtofp(gmul(gel(arc, 5), fact), prec)), rtodbl(gtofp(gmul(gel(arc, 6), fact), prec)));
  }
  fclose(f);
  avma=top;
}


//HELPER METHODS



//Returns the angle ang-bot in the range [0, 2*Pi)
GEN anglediff(GEN ang, GEN bot, GEN tol, long prec){
  pari_sp top=avma;
  GEN twopi=Pi2n(1, prec);
  GEN angdiff=gmod(gsub(ang, bot), twopi);
  if(toleq(angdiff, twopi, tol, prec) || toleq(angdiff, gen_0, tol, prec)){avma=top;return gen_0;}
  return gerepileupto(top, angdiff);
}

//Arctan of x, where oo is allowed as input (returning Pi/2)
GEN atanoo(GEN x, long prec){
  pari_sp top=avma;
  if(typ(x)==t_INFINITY) return gerepileupto(top, gdivgs(mppi(prec), 2));
  return gatan(x, prec);
}

//Returns gcmp(x, y), except if x==y, returns -1. Useful for gen_search when you ALWAYS want to return the index where it should be inserted/is
int gcmp_strict(void *data, GEN x, GEN y){
  int c=gcmp(x, y);
  if(c==0) c=-1;
  return c;
}

//Returns 0 if c is a circle, 1 if c is a line, 2 if c is a circle arc, 3 if c is line segment, and -1 if none of the above
int geom_check(GEN c){
  if(typ(c)!=t_VEC) return -1;
  int len=lg(c);
  if(len==CIRCLEN){//Circle or line
    if(gequal0(gel(c, 3))) return 0;//Circle
    if(gequal(gel(c, 3), gen_1)) return 1;//Line
    return -1;//Neither
  }
  else if(len==ARCLEN){
    if(gequal0(gel(c, 8))) return 2;//Circle arc
    if(gequal(gel(c, 8), gen_1)) return 3;//Line segment
    return -1;//Neither
  }
  return -1;//Not a circle/line
}

//Returns the default tolerance given the precision.
GEN deftol(long prec){
  pari_sp top=avma;
  return gerepileupto(top, gtofp(powis(gen_2, 32*(2-prec)), prec));
}

//Shifts the given angle ang by multiples of 2*Pi into the range [bot, bot+2*Pi).
GEN shiftangle(GEN ang, GEN bot, GEN tol, long prec){
  pari_sp top=avma;
  GEN diff=anglediff(ang, bot, tol, prec);
  if(gequal0(diff)){avma=top;return gcopy(bot);}
  return gerepileupto(top, gadd(bot, diff));
}

//Returns -1 if x<y, 0 if x==y, 1 if x>y (x, y are t_REAL). Accounts for the tolerance, so will deem x==y if they are equal up to tol AND at least one is inexact
long tolcmp(GEN x, GEN y, GEN tol, long prec){
  if(typ(x)==t_INFINITY || typ(y)==t_INFINITY) return gcmp(x, y);//No precision concerns
  pari_sp top=avma;
  GEN d=gsub(x, y);
  if(precision(d)==0){
    long ret=gsigne(d);
    avma=top;
    return ret;
  }//Exact objects
  if(gcmp(gabs(d, prec), tol)==-1){avma=top;return 0;}//Within tolerance
  long ret=gsigne(d);
  avma=top;return ret;
}

//Data points to [tol, vecsmall(prec)]. Used to sort/search a list with tolerance.
int tolcmp_sort(void *data, GEN x, GEN y){
  return tolcmp(x, y, gel(*(GEN*)data, 1), gel(*(GEN*)data, 2)[1]);
}

//DELETE???
//Data points to [tol, vecsmall(prec)]. Same as tolcmp_sort, but if x==y, returns -1 (so never returns 0). Used in reduceelt_givennormbound, as we ALWAYS want to search for the index to insert the element, and not to return 0 if it exists. In particular, this only accepts finite real values.
int tolcmp_strict(void *data, GEN x, GEN y){
  pari_sp top=avma;
  GEN tol=gel(*(GEN*)data, 1);
  long prec=gel(*(GEN*)data, 2)[1];
  GEN d=gsub(x, y);
  if(precision(d)==0){
    long ret=gsigne(d);
    if(ret==0) ret=-1;
    avma=top;
    return ret;
  }//Exact objects
  if(gcmp(gabs(d, prec), tol)==-1){avma=top;return -1;}//Within tolerance
  long ret=gsigne(d);
  avma=top;return ret;
}

//Returns 1 if x==y, and 0 if x!=y. If x or y is not a precise objects (e.g. t_REAL), will return 1 if they are equal up to the tolerance tol.
int toleq(GEN x, GEN y, GEN tol, long prec){
  if(typ(x)==t_INFINITY || typ(y)==t_INFINITY || gequal0(tol)) return gequal(x, y);//No precision concerns
  pari_sp top=avma;
  GEN d=gsub(x, y);
  if(gequal0(d)){avma=top;return 1;}//Deemed equal already
  if(precision(d)==0){avma=top;return 0;}//Exact objects
  if(gcmp(gabs(d, prec), tol)==-1){avma=top;return 1;}//Within tolerance
  avma=top;return 0;
}

