//Quaternionic methods over Q

#ifndef PARILIB
#define PARILIB
#include <pari/pari.h>
#endif

#ifndef METHDECL
#define METHDECL
#include "qquadraticdecl.h"
#endif

//The length (lg, so technically length+1) of a quaternion algebra entry and an initialized quaternion order.
#define QALEN 5
#define QAORDLEN 8

//STATIC DECLARATIONS
static GEN qa_ord_type(GEN Q, GEN ord, GEN level);
static GEN qa_init_m2z(void);
static GEN qa_ord_init_trace0basis(GEN Q, GEN ord, GEN maxds);
static GEN qa_conjbasis_orient(GEN Q, GEN ord, GEN v1, GEN v2, GEN e2);
static int qa_embed_isnewoptimal(GEN Q, GEN ord, GEN ordinv, GEN D, GEN Dmod2, GEN dfacs, GEN emb, GEN gcdf1g1h1, GEN embs, long pos, long prec);
static int qa_embedor_compare(void *data, GEN pair1, GEN pair2);

//A quaternion algebra Q is stored as the vector [nf, pset, [a,b,-ab], product of finite primes in pset] where it is an algebra over the field nf, pset is the set of ramified primes, and the algebra can be represented as (a,b/nf), so that i^2=a, j^2=b, k^2=-ab. Over Q, we use nf=0, and by convention will always have gcd(a, b)=1.

//An order is stored as order=[ord, type, [d1, d2, d3, d4], level, prime factorization of the level, ord^(-1), [basis of ord intersect trace 0 elements]]. The column vectors of ord form a Z-basis of the order. type=-1 means a general order, =0 is maximal, and =1 is Eichler. d_i is the maximal denominator appearing in coefficients of 1,i,j,k. The prime factorization is stored as [[p_1, e_1],...,[p_n, e_n]].

//Orders can sometimes be passed in as the 4x4 matrix of column vectors, OR containing the full suite of computed data. The typecheck methods will in general allow either. Other methods will take the first type if the variable is named "ord", and the second type if it is named "order".




//BASIC OPERATIONS ON ELEMENTS IN QUATERNION ALGEBRAS 



//Returns the conjugate of the quaternion element x. Note that the quaternion algebra is not required as an input.
GEN qa_conj(GEN x){
  long lx;
  GEN ret=cgetg_copy(x, &lx);
  gel(ret, 1)=gcopy(gel(x, 1));
  for(int i=2;i<5;i++) gel(ret, i)=gneg(gel(x, i));
  return ret;//No garbage
}

//qa_conj with typechecking
GEN qa_conj_typecheck(GEN x){
  qa_eltcheck(x);
  return qa_conj(x);
}

//Returns yxy^(-1)
GEN qa_conjby(GEN Q, GEN x, GEN y){
  pari_sp top=avma;
  GEN yinv=qa_inv(Q, y);
  GEN yx=qa_mul(Q, y, x);
  return gerepileupto(top, qa_mul(Q, yx, yinv));
}

//qa_conjby with typecheck
GEN qa_conjby_typecheck(GEN Q, GEN x, GEN y){
  qa_check(Q);qa_eltcheck(x);qa_eltcheck(y);
  return qa_conjby(Q, x, y);
}

//Returns the inverse of the quaternion element x.
GEN qa_inv(GEN Q, GEN x){
  pari_sp top=avma;
  GEN qconj=qa_conj(x);
  GEN n=qa_norm(Q, x);
  if(gequal0(n)) pari_err_INV("Cannot invert norm 0 element", x);
  return gerepileupto(top, gdiv(qconj,n));
}

//qa_inv with typechecking
GEN qa_inv_typecheck(GEN Q, GEN x){
  qa_check(Q);qa_eltcheck(x);
  return qa_inv(Q, x);
}

//Given an indefinite quaternion algebra with a>0, this outputs the image of x under the standard embedding into M_2(R). This is given by 1->1, i->[sqrt(a),0;0,-sqrt(a)], j->[0,b;1,0], k->[0,b sqrt(a);-sqrt(a),0].
GEN qa_m2rembed(GEN Q, GEN x){
  pari_sp top=avma;
  GEN rta;
  if(lg(qa_getpram(Q))==1) rta=gen_1;//quadgen only takes non-squares
  else rta=quadroot(qa_geta(Q));
  GEN x2rt=gmul(gel(x, 2), rta), x4rt=gmul(gel(x, 4), rta);
  GEN x3px4rt=gadd(gel(x, 3), x4rt);
  GEN emb=cgetg(3, t_MAT);
  gel(emb, 1)=cgetg(3, t_COL), gel(emb, 2)=cgetg(3, t_COL);
  gcoeff(emb, 1, 1)=gadd(gel(x, 1), x2rt);
  gcoeff(emb, 1, 2)=gmul(qa_getb(Q), x3px4rt);
  gcoeff(emb, 2, 1)=gsub(gel(x, 3), x4rt);
  gcoeff(emb, 2, 2)=gsub(gel(x, 1), x2rt);
  return gerepileupto(top, emb);
}

//qa_m2r_embed with typecheck
GEN qa_m2rembed_typecheck(GEN Q, GEN x){
  qa_indefcheck(Q);
  qa_eltcheck(x);
  return qa_m2rembed(Q, x);
}

//Returns the min poly of x, i.e [1, b, c] for x^2+bx+c or [1, b] for x+b
GEN qa_minpoly(GEN Q, GEN x){
  if(gequal0(gel(x, 2)) && gequal0(gel(x, 3)) && gequal0(gel(x, 4))){//Rational
	GEN ret=cgetg(3, t_VEC);gel(ret, 1)=gen_1;gel(ret, 2)=gneg(gel(x, 1));
	return ret;
  }
  GEN ret=cgetg(4, t_VEC);
  gel(ret, 1)=gen_1;
  gel(ret, 2)=gmulgs(gel(x, 1), -2);//- the trace
  gel(ret, 3)=qa_norm(Q, x);
  return ret;
}

//qa_minpoly wiht typechecking
GEN qa_minpoly_typecheck(GEN Q, GEN x){
  qa_check(Q);qa_eltcheck(x);
  return qa_minpoly(Q, x);
}

//Multiplies x by y in Q. Note that x, y CAN be column vectors, though the result is still a row vector.
GEN qa_mul(GEN Q, GEN x, GEN y){
  pari_sp top=avma;
  GEN a=qa_geta(Q), b=qa_getb(Q), mab=qa_getmab(Q);//The a and b and -ab for the q-alg
  GEN A1A2=gmul(gel(x, 1), gel(y, 1)), aB1B2=gmul(a, gmul(gel(x, 2), gel(y, 2))), bC1C2=gmul(b, gmul(gel(x, 3), gel(y, 3))), mabD1D2=gmul(mab, gmul(gel(x, 4), gel(y, 4)));
  GEN A1A2paB1B2=gadd(A1A2, aB1B2), bC1C2mabD1D2=gadd(bC1C2, mabD1D2);//The two parts of the 1 coefficient
  GEN A1B2=gmul(gel(x, 1), gel(y, 2)), A2B1=gmul(gel(x, 2), gel(y, 1)), C1D2=gmul(gel(x, 3), gel(y, 4)), C2D1=gmul(gel(x, 4), gel(y, 3));
  GEN A1B2pA2B1=gadd(A1B2, A2B1), mbC1D2pC2D1=gmul(b, gsub(C2D1, C1D2));//The two parts of the i coefficient
  GEN A1C2=gmul(gel(x, 1), gel(y, 3)), A2C1=gmul(gel(x, 3), gel(y, 1)), B1D2=gmul(gel(x, 2), gel(y, 4)), B2D1=gmul(gel(x, 4), gel(y, 2));
  GEN A1C2pA2C1=gadd(A1C2, A2C1), aB1D2maB2D1=gmul(a, gsub(B1D2, B2D1));//The two parts of the j coefficient
  GEN A1D2=gmul(gel(x, 1), gel(y, 4)), A2D1=gmul(gel(x, 4), gel(y, 1)), B1C2=gmul(gel(x, 2), gel(y, 3)), B2C1=gmul(gel(x, 3), gel(y, 2));
  GEN A1D2pA2D1=gadd(A1D2, A2D1), B1C2mB2C1=gsub(B1C2, B2C1);//The two parts of the k coefficient
  long lx;
  GEN ret=cgetg_copy(x, &lx);
  gel(ret, 1)=gadd(A1A2paB1B2, bC1C2mabD1D2);
  gel(ret, 2)=gadd(A1B2pA2B1, mbC1D2pC2D1);
  gel(ret, 3)=gadd(A1C2pA2C1, aB1D2maB2D1);
  gel(ret, 4)=gadd(A1D2pA2D1, B1C2mB2C1);
  return gerepileupto(top, ret);
}

//qa_mul with typechecking
GEN qa_mul_typecheck(GEN Q, GEN x, GEN y){
  qa_check(Q);qa_eltcheck(x);qa_eltcheck(y);
  return qa_mul(Q, x, y);
}

//Returns the norm of x
GEN qa_norm(GEN Q, GEN x){
  pari_sp top=avma;
  GEN a=qa_geta(Q), b=qa_getb(Q);
  GEN t1=gsub(gsqr(gel(x, 1)), gmul(gsqr(gel(x, 2)), a));
  GEN t2=gmul(b, gsub(gmul(a, gsqr(gel(x, 4))), gsqr(gel(x, 3))));
  return gerepileupto(top, gadd(t1, t2));
}

//qa_norm with typechecking
GEN qa_norm_typecheck(GEN Q, GEN x){
  qa_check(Q);qa_eltcheck(x);
  return qa_norm(Q, x);
}

//Powers x to the n
GEN qa_pow(GEN Q, GEN x, GEN n){
  pari_sp top=avma;
  if(gequal0(n)) return mkvec4(gen_1, gen_0, gen_0, gen_0);//1
  if(signe(n)==-1){x=qa_inv(Q, x);n=negi(n);}
  GEN dig=binary_zv(n);
  GEN qpow=gcopy(x);
  for(long i=1;i<lg(dig);i++){
	qpow=qa_square(Q, qpow);
	if(dig[i]==1) qpow=qa_mul(Q, qpow, x);
  }
  return gerepileupto(top, qpow);
}

//qa_pow with typechecking
GEN qa_pow_typecheck(GEN Q, GEN x, GEN n){
  qa_check(Q);qa_eltcheck(x);
  if(typ(n)!=t_INT) pari_err_TYPE("Please enter an integer power", n);
  return qa_pow(Q, x, n);
}

//Given an x=[e,f,g,h] in a quaternion algebra Q, this returns the roots of [e,f,g,h] under the standard embedding into SL(2,R), with the "first root" coming first. We must have a>0 for this to be good.
GEN qa_roots(GEN Q, GEN x, long prec){
  pari_sp top=avma;
  if(gsigne(gel(x, 1))==-1) x=gneg(x);//Make x have positive trace
  GEN abvec=qa_getabvec(Q);
  GEN n=qa_norm(Q, x), x1sq=gsqr(gel(x, 1));
  if(gsigne(n)!=1 || gcmp(x1sq, n)==-1) pari_err_TYPE("Not a positive norm hyperbolic element", x);
  GEN roota=gsqrt(gel(abvec, 1), prec);
  GEN den=gsub(gel(x, 3), gmul(gel(x, 4), roota));
  if(gequal0(gel(x, 3)) && gequal0(gel(x, 4))){//x[3]=x[4]=0
	int s=gsigne(gel(x, 2));
	GEN ret=cgetg(3, t_VEC);
	if(s==1){gel(ret, 1)=mkoo();gel(ret, 2)=gen_0;}
	else if(s==-1){gel(ret, 1)=gen_0;gel(ret, 2)=mkoo();}
	else pari_err_TYPE("x must not be rational", x);
	return gerepileupto(top, ret);
  }
  else if(gequal0(den)){//We must be in the M_2(Q) case, presumably with a=1
    GEN rtpart=gdiv(gmul(gadd(gel(x, 3), gel(x, 4)), gel(abvec, 2)), gmulgs(gel(x, 2), 2));//b(x[3]+x[4])/(2x[2])
	GEN ret=cgetg(3, t_VEC);
	if(gsigne(gel(x, 2))==1){
	  gel(ret, 1)=mkoo();
	  gel(ret, 2)=gneg(rtpart);
	}
	else{
	  gel(ret, 1)=gneg(rtpart);
	  gel(ret, 2)=mkoo();
	}
	return gerepileupto(top, ret);
  }//Now x[3]-sqrt(a)*x[4]=den is nonzero
  GEN x2rt1=gmul(gel(x, 2), roota), rootpart=gsqrt(gsub(x1sq, n), prec);
  GEN r1num=gadd(x2rt1, rootpart);
  GEN r2num=gsub(x2rt1, rootpart);
  GEN ret=cgetg(3, t_VEC);
  gel(ret, 1)=gdiv(r1num, den);
  gel(ret, 2)=gdiv(r2num, den);
  return gerepileupto(top, ret);
}

//qa_roots with typecheck
GEN qa_roots_typecheck(GEN Q, GEN x, long prec){
  qa_indefcheck(Q);
  qa_eltcheck(x);
  return qa_roots(Q, x, prec);
}

//Returns the square of x in Q
GEN qa_square(GEN Q, GEN x){
  pari_sp top=avma;
  GEN a=qa_geta(Q), b=qa_getb(Q);
  GEN A2paB2=gadd(gsqr(gel(x, 1)), gmul(a, gsqr(gel(x, 2)))), bC2mabD2=gmul(b, gsub(gsqr(gel(x, 3)), gmul(a, gsqr(gel(x, 4))))), twoA=gmulgs(gel(x, 1), 2);
  long lx;
  GEN ret=cgetg_copy(x, &lx);
  gel(ret, 1)=gadd(A2paB2, bC2mabD2);
  gel(ret, 2)=gmul(twoA, gel(x, 2));
  gel(ret, 3)=gmul(twoA, gel(x, 3));
  gel(ret, 4)=gmul(twoA, gel(x, 4));
  return gerepileupto(top, ret);
}

//qa_square with typechecking
GEN qa_square_typecheck(GEN Q, GEN x){
  qa_check(Q);qa_eltcheck(x);
  return qa_square(Q, x);
}

//Returns the trace of x, an element of a quaternion algebra.
GEN qa_trace(GEN x){return gmulgs(gel(x, 1), 2);}

//qa_trace with typechecking
GEN qa_trace_typecheck(GEN x){qa_eltcheck(x);return qa_trace(x);}



//BASIC OPERATIONS ON ORDERS/LATTICES IN QUATERNION ALGEBRAS 



//Checks if x is in the order specificed by ordinv^(-1).
int qa_isinorder(GEN Q, GEN ordinv, GEN x){
  pari_sp top=avma;
  GEN xcol=gtocol(x);
  GEN v=gmul(ordinv, xcol);
  if(typ(Q_content(v))==t_INT){avma=top;return 1;}//No denominators, in order.
  avma=top;return 0;//Denominators, not in order
}

//qa_isinorder with typechecking, and taking in the order or initialized order (and not its inverse).
int qa_isinorder_typecheck(GEN Q, GEN ord, GEN x){
  pari_sp top=avma;
  qa_check(Q);qa_eltcheck(x);
  if(typ(ord)==t_VEC && lg(ord)==QAORDLEN) return qa_isinorder(Q, qa_getordinv(ord), x);
  QM_check(ord);
  int isord=qa_isinorder(Q, ginv(ord), x);
  avma=top;return isord;
}

//Conjugates the order ord by the (invertible element) c.
GEN qa_ord_conj(GEN Q, GEN ord, GEN c){
  pari_sp top=avma;
  GEN neword=cgetg(5, t_MAT);
  gel(neword, 1)=qa_conjby(Q, gel(ord, 1), c);settyp(gel(neword, 1), t_COL);//Making the columns and converting from t_VEC to t_COL
  gel(neword, 2)=qa_conjby(Q, gel(ord, 2), c);settyp(gel(neword, 2), t_COL);
  gel(neword, 3)=qa_conjby(Q, gel(ord, 3), c);settyp(gel(neword, 3), t_COL);
  gel(neword, 4)=qa_conjby(Q, gel(ord, 4), c);settyp(gel(neword, 4), t_COL);
  return gerepileupto(top, QM_hnf(neword));
}

//qa_ord_conj with typechecking
GEN qa_ord_conj_typecheck(GEN Q, GEN ord, GEN c){
  qa_check(Q);qa_eltcheck(c);
  GEN actualord=qa_ordcheck(ord);
  return qa_ord_conj(Q, actualord, c);
}

//Returns the discriminant of the order ord.
GEN qa_ord_disc(GEN Q, GEN ord){
  pari_sp top=avma;
  GEN M=cgetg(5, t_MAT);
  for(long i=1;i<5;i++) gel(M, i)=cgetg(5, t_COL);
  for(long i=1;i<5;i++){
	for(long j=1;j<5;j++){
      gcoeff(M, i, j)=qa_trace(qa_mul(Q, gel(ord, i), gel(ord, j)));
	}
  }
  GEN d=absi(ZM_det(M));
  return gerepileupto(top, sqrti(d));
}

//qa_ord with typechecking (we do NOT check that ord does give an order however).
GEN qa_ord_disc_typecheck(GEN Q, GEN ord){
  qa_check(Q);
  GEN actualord=qa_ordcheck(ord);
  return qa_ord_disc(Q, actualord);
}

//IMPLEMENT EICHLER THIS IS MISSING

//Returns the type of the order, 0 if maximal, 1 if Eichler, -1 else.
static GEN qa_ord_type(GEN Q, GEN ord, GEN level){
  if(equali1(level)) return gen_0;
  pari_warn(warner, "Checking for Eichler has not yet been implemented.");
  return gen_m1;
}



//INITIALIZATION METHODS



//Initialized a quaternion order with relevant data.
GEN qa_ord_init(GEN Q, GEN ord){
  pari_sp top=avma;
  GEN disc=qa_ord_disc(Q, ord);
  GEN level=diviiexact(disc, qa_getpramprod(Q));
  GEN lfac=Z_factor(level);
  long nprimes=lgcols(lfac);//nprimes-1=number of distinct prime divisors of level.
  GEN r1=row(ord, 1), r2=row(ord, 2), r3=row(ord, 3), r4=row(ord, 4);
  GEN ret=cgetg(QAORDLEN, t_VEC);
  gel(ret, 1)=QM_hnf(ord);
  gel(ret, 2)=qa_ord_type(Q, ord, level);
  gel(ret, 3)=cgetg(5, t_VEC);
  gel(gel(ret, 3), 1)=Q_denom(r1);gel(gel(ret, 3), 2)=Q_denom(r2);gel(gel(ret, 3), 3)=Q_denom(r3);gel(gel(ret, 3), 4)=Q_denom(r4);
  gel(ret, 4)=icopy(level);
  gel(ret, 5)=cgetg(nprimes, t_VEC);
  for(long i=1;i<nprimes;i++){gel(gel(ret, 5), i)=mkvec2copy(gcoeff(lfac, i, 1), gcoeff(lfac, i, 2));}
  gel(ret, 6)=ginv(gel(ret, 1));
  gel(ret, 7)=qa_ord_init_trace0basis(Q, gel(ret, 1), gel(ret, 3));
  return gerepileupto(top, ret);
}

//qa_ord_init with typecheck
GEN qa_ord_init_typecheck(GEN Q, GEN ord){
  qa_check(Q);QM_check(ord);
  return qa_ord_init(Q, ord);
}

//Primes must be sorted, and the oo prime is present. If type=1 assumes indefinite, and type=-1 is definite.
GEN qa_init_primes(GEN pset, int type){
  pari_sp top=avma;
  if(lg(pset)==1) return qa_init_m2z();
  GEN prodp=gen_1, relations, extracong, u;//We initiate a prime search with the correct search conditions.
  if(type==-1){//Definite
    for(long i=1;i<lg(pset)-1;i++) prodp=mulii(prodp, gel(pset, i));//Product of all the finite primes
	u=gen_m1;
	if(equalii(gel(pset, 1), gen_2)){//2 ramifies (it is the smallest prime, so must be first as pset is assumed to be sorted).
	  relations=cgetg(lg(pset)-2, t_VEC);
	  for(long i=1;i<lg(pset)-2;i++) gel(relations, i)=mkvec2(gel(pset, i+1), stoi(-kronecker(gen_m1, gel(pset, i+1))));//u*q is  a non-residue for each odd prime ramifying
	  extracong=mkvec2s(8, 3);
	}
	else{//2 does not ramify
	  relations=cgetg(lg(pset)-1, t_VEC);
	  for(long i=1;i<lg(pset)-1;i++) gel(relations, i)=mkvec2(gel(pset, i), stoi(-kronecker(gen_m1, gel(pset, i))));//u*q is  a non-residue for each odd prime ramifying
	  extracong=mkvec2s(8, 7);
	}
  }
  else{//Indefinite
    for(long i=1;i<lg(pset);i++) prodp=mulii(prodp, gel(pset, i));//Product of all the finite primes
	u=gen_1;
	if(equalii(gel(pset, 1), gen_2)){//2 ramifies (it is the smallest prime, so must be first as pset is assumed to be sorted).
	  relations=cgetg(lg(pset)-1, t_VEC);
	  for(long i=1;i<lg(pset)-1;i++) gel(relations, i)=mkvec2(gel(pset, i+1), gen_m1);//u*q is  a non-residue for each odd prime ramifying
	  extracong=mkvec2s(8, 5);
	}
	else{//2 does not ramify
	  long lx;
	  relations=cgetg_copy(pset, &lx);
	  for(long i=1;i<lg(pset);i++) gel(relations, i)=mkvec2(gel(pset, i), gen_m1);//u*q is  a non-residue for each odd prime ramifying
	  extracong=mkvec2s(8, 1);
	}
  }
  GEN p=prime_ksearch(relations, extracong);
  GEN alg=cgetg(5, t_VEC);
  gel(alg, 1)=gen_0;
  gel(alg, 2)=gcopy(pset);
  gel(alg, 3)=cgetg(4, t_VEC);
  gel(gel(alg, 3), 1)=mulii(prodp, u);
  gel(gel(alg, 3), 2)=mulii(u, p);
  gel(gel(alg, 3), 3)=mulii(prodp, p);
  togglesign_safe(&gel(gel(alg, 3), 3));//Fixing to -ab
  gel(alg, 4)=icopy(prodp);
  return gerepileupto(top, alg);
}

//qa_init_primes with sorting of pset, adding in oo if missing. Checks for positive integers but NOT distinct primes.
GEN qa_init_primes_typecheck(GEN pset){
  pari_sp top=avma;
  if(typ(pset)!=t_VEC) pari_err_TYPE("Please enter a vector of primes", pset);
  GEN psort=sort(pset);
  for(long i=1;i<lg(psort);i++){
	if(typ(gel(psort, i))!=t_INT || signe(gel(psort, i))==-1){
	  if(typ(gel(psort, i))!=t_INFINITY) pari_err_TYPE("Please enter a vector of primes", pset);
	}
  }
  long l=lg(psort);
  if(l%2==1){
	if(typ(gel(psort, l-1))==t_INFINITY) return gerepileupto(top, qa_init_primes(psort, -1));//Definite
	return gerepileupto(top, qa_init_primes(psort, 1));//Indefinite
  }
  if(typ(gel(psort, l-1))==t_INFINITY) pari_err_TYPE("Odd length list with oo included not allowed", pset);
  GEN psortoo=cgetg(l+1, t_VEC);
  for(long i=1;i<l;i++) gel(psortoo, i)=gel(psort, i);
  gel(psortoo, l)=mkoo();//Adding in oo prime
  return gerepileupto(top, qa_init_primes(psortoo, -1));//Definitie
}

//Initializes the quaternion algebra for M_2(Z)
static GEN qa_init_m2z(void){
  pari_sp top=avma;
  GEN ord=zeromatcopy(4, 4);
  gcoeff(ord, 1, 1)=gen_1;gcoeff(ord, 1, 2)=ghalf;
  gcoeff(ord, 2, 2)=gneg(ghalf);
  gcoeff(ord, 3, 3)=gen_1;gcoeff(ord, 3, 4)=ghalf;
  gcoeff(ord, 4, 4)=gneg(ghalf);
  GEN ret=cgetg(3, t_VEC);
  GEN alg=cgetg(5, t_VEC);
  gel(alg, 1)=gen_0;
  gel(alg, 2)=cgetg(1, t_VEC);
  gel(alg, 3)=mkvec3(gen_1, gen_1, gen_m1);
  gel(alg, 4)=gen_1;
  gel(ret, 1)=alg;
  gel(ret, 2)=qa_ord_init(alg, ord);
  return gerepileupto(top, ret);
}

//Initializes a quaternion algebra over Q with two primes and a corresponding maximal order
GEN qa_init_2primes(GEN p, GEN q){
  pari_sp top=avma;
  if(cmpii(p,q)==-1){GEN r=q;q=p;p=r;}//WLOG p>q
  GEN a, b, ord=zeromatcopy(4, 4);
  if(equalii(q,gen_2)){//Case of q=2
    long m8=smodis(p, 8);
	GEN r;
	switch(m8){
	  case 1://p==1(8) and q=2; (2p, -r)
	    r=prime_ksearch(mkvec(mkvec2(p, gen_m1)), mkvec2(stoi(8), stoi(3)));
		a=shifti(p, 1);
		b=negi(r);
		GEN x=Zp_sqrt(shifti(p, 1), r, 1);//sqrt(2p) modulo r
		gcoeff(ord, 1, 1)=gen_1;gcoeff(ord, 1, 2)=ghalf;
		gcoeff(ord, 2, 3)=ghalf;
		gcoeff(ord, 3, 2)=ghalf;gcoeff(ord, 3, 4)=Qdivii(x, r);
		gcoeff(ord, 4, 3)=ghalf;gcoeff(ord, 4, 4)=Qdivii(gen_1, r);
		break;
	  case 5://p==5(8) and q=2; (p, 2)
	    a=p;
		b=gen_2;
		gcoeff(ord, 1, 1)=gen_1;gcoeff(ord, 1, 2)=ghalf;
		gcoeff(ord, 2, 2)=ghalf;
		gcoeff(ord, 3, 3)=gen_1;gcoeff(ord, 3, 4)=ghalf;
		gcoeff(ord, 4, 4)=ghalf;
		break;
	  default: //p==3(4); (p, -1)
	    a=p;
		b=gen_m1; 
		gcoeff(ord, 1, 1)=gen_1;gcoeff(ord, 1, 4)=ghalf;
		gcoeff(ord, 2, 2)=gen_1;gcoeff(ord, 2, 4)=ghalf;
		gcoeff(ord, 3, 3)=gen_1;gcoeff(ord, 3, 4)=ghalf;
		gcoeff(ord, 4, 4)=ghalf;
	}
  }
  else{//Now q>2
    long p1=smodis(p, 4), q1=smodis(q, 4);
    if(p1==3 && q1==3){//p==q==3 (4); (pq, -1)
	  a=mulii(p, q);
      b=gen_m1;
	  gcoeff(ord, 1, 1)=gen_1;gcoeff(ord, 1, 2)=ghalf;
	  gcoeff(ord, 2, 2)=ghalf;
	  gcoeff(ord, 3, 3)=gen_1;gcoeff(ord, 3, 4)=ghalf;
	  gcoeff(ord, 4, 4)=ghalf;
    }
	else{//p==1 (4) or q==1 (4)
	  if(kronecker(p, q)==-1){//(p/q)=-1
		a=p;
		b=q;
		gcoeff(ord, 1, 1)=gen_1;gcoeff(ord, 1, 2)=ghalf;gcoeff(ord, 1, 4)=ghalf;
		gcoeff(ord, 2, 4)=ghalf;
		gcoeff(ord, 3, 4)=ghalf;
		gcoeff(ord, 4, 4)=ghalf;
		if(p1==1){
		  gcoeff(ord, 2, 2)=ghalf;gcoeff(ord, 3, 3)=gen_1;
		}
		else{
		  gcoeff(ord, 2, 3)=gen_1;gcoeff(ord, 3, 2)=ghalf;
		}
	  }
	  else{//(p/q)=1
	    GEN s1, s2;
		if(p1==1) s1=gen_m1;
		else s1=gen_1;
		if(q1==1) s2=gen_m1;
		else s2=gen_1;//Getting the signs of -p%4 and -q%4
		GEN r=prime_ksearch(mkvec2(mkvec2(p, s1), mkvec2(q, s2)), mkvec2(stoi(4), stoi(3)));
		a=mulii(p, q);
		b=negi(r);
		GEN x=Zp_sqrt(mulii(p, q), r, 1);//sqrt(pq) modulo r
		gcoeff(ord, 1, 1)=gen_1;gcoeff(ord, 1, 2)=ghalf;
		gcoeff(ord, 2, 3)=ghalf;
		gcoeff(ord, 3, 2)=ghalf;gcoeff(ord, 3, 4)=Qdivii(x, r);
		gcoeff(ord, 4, 3)=ghalf;gcoeff(ord, 4, 4)=Qdivii(gen_1, r);
	  }
	}
  }
  GEN ret=cgetg(3, t_VEC);
  GEN alg=cgetg(5, t_VEC);
  gel(alg, 1)=gen_0;
  gel(alg, 2)=mkvec2copy(q, p);
  gel(alg, 3)=cgetg(4, t_VEC);
  gel(gel(alg, 3), 1)=icopy(a);
  gel(gel(alg, 3), 2)=icopy(b);
  gel(gel(alg, 3), 3)=mulii(a,b);
  togglesign_safe(&gel(gel(alg, 3), 3));//Fixing to -ab
  gel(alg, 4)=mulii(p, q);
  gel(ret, 1)=alg;
  gel(ret, 2)=qa_ord_init(alg, ord);
  return gerepileupto(top, ret);
}

//Initializes a quaternion algebra over Q with two primes and checks the inputs are valid.
GEN qa_init_2primes_typecheck(GEN p, GEN q){
  if(typ(p)!=t_INT || typ(q)!=t_INT || equalii(p, q) || !isprime(p) || !isprime(q)) pari_err_TYPE("Please input two distinct prime numbers", mkvec2(p, q));
  return qa_init_2primes(p, q);
}

//Returns generators for the basis of the elements of ord with trace zero
static GEN qa_ord_init_trace0basis(GEN Q, GEN ord, GEN maxds){
  pari_sp top=avma;
  GEN tracezero=zeromatcopy(4, 3);
  gcoeff(tracezero, 2, 1)=gdivsg(1, gel(maxds, 2));
  gcoeff(tracezero, 3, 2)=gdivsg(1, gel(maxds, 3));
  gcoeff(tracezero, 4, 3)=gdivsg(1, gel(maxds, 4));
  GEN space=module_intersect(ord, tracezero);
  GEN ret=cgetg(4, t_VEC);
  gel(ret, 1)=gtovec(gel(space, 1));
  gel(ret, 2)=gtovec(gel(space, 2));
  gel(ret, 3)=gtovec(gel(space, 3));
  return gerepileupto(top, ret);
}

//Given an a and a b, returns the set of primes ramifying in the quaternion algebra ramified at a, b
GEN qa_ram_fromab(GEN a, GEN b){
  pari_sp top=avma;
  glist *S=NULL;
  long nram=0;
  if(signe(a)==-1 && signe(b)==-1){glist_putstart(&S, mkoo());nram++;}
  if(hilbertii(a, b, gen_2)==-1){glist_putstart(&S, gen_2);nram++;}
  GEN ashift, bshift;
  Z_pvalrem(a, gen_2, &ashift);
  Z_pvalrem(b, gen_2, &bshift);
  if(signe(ashift)==-1) ashift=negi(ashift);
  if(signe(bshift)==-1) bshift=negi(bshift);//Now ashift, bshift are oddd and positive, so we factorize
  GEN fact=Z_factor(ashift), p;
  for(long i=1;i<lg(gel(fact, 1));i++){
	p=gcoeff(fact, i, 1);
	if(hilbertii(a, b, p)==-1){glist_putstart(&S, p);nram++;}//-1, so ramifies!
  }
  fact=Z_factor(bshift);
  for(long i=1;i<lg(gel(fact, 1));i++){
	p=gcoeff(fact, i, 1);
	if(hilbertii(a, b, p)==-1){glist_putstart(&S, p);nram++;}//-1, so ramifies!
  }
  GEN pset=glist_togvec(S, nram, -1);
  return gerepileupto(top, gen_sort_uniq(pset, NULL, &cmp_data));
}

//qa_ram_fromab with typechecking
GEN qa_ram_fromab_typecheck(GEN a, GEN b){
  if(typ(a)!=t_INT || gequal0(a)) pari_err_TYPE("a, b must be non-zero integers", a);
  if(typ(b)!=t_INT || gequal0(b)) pari_err_TYPE("a, b must be non-zero integers", b);
  return qa_ram_fromab(a, b);
}



//CONJUGATION OF ELEMENTS IN A GIVEN ORDER



//Returns the basis for the 2-dimensional Z-module of the set of x for which x*e1=e2*x. Returns 0 if e1, e2 are not conjugate or rational.
GEN qa_conjbasis(GEN Q, GEN ord, GEN ordinv, GEN e1, GEN e2, int orient){
  pari_sp top=avma;
  if(!gequal(gel(e1, 1), gel(e2, 1))) return gen_0;//Traces must be equal.
  if((gequal0(gel(e1, 2)) && gequal0(gel(e1, 3)) && gequal0(gel(e1, 4))) || (gequal0(gel(e2, 2)) && gequal0(gel(e2, 3)) && gequal0(gel(e2, 4)))) return gen_0;//Don't allow rational.
  GEN a=qa_geta(Q), b=qa_getb(Q), mab=qa_getmab(Q);
  //Solving v*e1=e2*v yields that v must be in the kernel of W, described as follows:
  //[m,n,p,q]=e1, [m,r,x,y]=e2, W=[0,a*(n-r),b*(p-x),a*b*(y-q);n-r,0,-b*(y+q),b*(p+x);p-x,a*(q+y),0,-a*(n+r);y-q,-x-p,r+n,0]
  GEN W=zeromatcopy(4, 4);
  GEN nmr=gsub(gel(e1, 2), gel(e2, 2)), pmx=gsub(gel(e1, 3), gel(e2, 3)), ymq=gsub(gel(e2, 4), gel(e1, 4));
  GEN npr=gadd(gel(e1, 2), gel(e2, 2)), ppx=gadd(gel(e1, 3), gel(e2, 3)), ypq=gadd(gel(e1, 4), gel(e2, 4));
  gcoeff(W, 1, 2)=gmul(a, nmr);gcoeff(W, 1, 3)=gmul(b, pmx);gcoeff(W, 1, 4)=gneg(gmul(mab, ymq));
  gcoeff(W, 2, 1)=nmr;gcoeff(W, 2, 3)=gneg(gmul(b, ypq));gcoeff(W, 2, 4)=gmul(b, ppx);
  gcoeff(W, 3, 1)=pmx;gcoeff(W, 3, 2)=gmul(a, ypq);gcoeff(W, 3, 4)=gneg(gmul(a, npr));
  gcoeff(W, 4, 1)=ymq;gcoeff(W, 4, 2)=gneg(ppx);gcoeff(W, 4, 3)=npr;
  GEN ker=matkerint0(W, 0);//Integer kernel. It is of rank 0 (if not conjugate) and 2 if conjugate.
  if(lg(ker)==1){avma=top;return gen_0;}//Done, no solution.
  //Now ker is comprised of two primitive ZC's v1, v2. We need to solve Av1+Bv2=ker*[A,B] is in the order ord.
  GEN ordker=shallowtrans(QM_mul(ordinv, ker));//Now [A, B]*ordker must be integral.
  GEN H=QM_hnf(ordker);//The corret space is now H^(-1)*ordker, and we multuply by ord to get back the coefficients of i, j, k
  GEN space=QM_mul(ord, shallowtrans(QM_mul(ginv(H), ordker)));
  if(orient==0){
	GEN ret=cgetg(3, t_VEC);
	gel(ret, 1)=gtovec(gel(space, 1));
	gel(ret, 2)=gtovec(gel(space, 2));
	return gerepileupto(top, ret);
  }
  GEN v1=gtovec(gel(space, 1));
  GEN v2=gtovec(gel(space, 2));
  return gerepileupto(top, qa_conjbasis_orient(Q, ord, v1, v2, e2));
}

//qa_conjbasis with initializing ordinv and typecheck
GEN qa_conjbasis_typecheck(GEN Q, GEN ord, GEN e1, GEN e2, int orient){
  pari_sp top=avma;
  qa_check(Q);qa_eltcheck(e1);qa_eltcheck(e2);
  GEN ordinv;
  if(typ(ord)==t_VEC && lg(ord)==QAORDLEN){ordinv=qa_getordinv(ord);ord=qa_getord(ord);}
  else{QM_check(ord);ordinv=ginv(ord);}
  return gerepileupto(top, qa_conjbasis(Q, ord, ordinv, e1, e2, orient));
}

//Orients the basis [v1, v2] as either that or [v1, -v2], so that v2*conj(v1)-v1*conj(v) has a positive ratio to phi_2(sqrt(D_2)), i.e. e2-trace(e2)/2.
static GEN qa_conjbasis_orient(GEN Q, GEN ord, GEN v1, GEN v2, GEN e2){
  pari_sp top=avma;
  GEN v1conj=qa_conj(v1), v2conj=qa_conj(v2);
  GEN x=gsub(qa_mul(Q, v2, v1conj), qa_mul(Q, v1, v2conj));//Needs to be a positive multiple of e2
  GEN e2shift=zerovec(4);
  for(int i=2;i<=4;i++) gel(e2shift, i)=gel(e2, i);
  GEN rat=vecratio(x, e2shift);
  if(gsigne(rat)==-1) v2=gneg(v2);//Must negate v2
  return gerepileupto(top, mkvec2copy(v1, v2));
}

//Returns [qf, v1, v2], where [v1, v2] is the oriented basis of qa_conjbasis and qf=nrd('x*v1+'y*v2).
GEN qa_conjqf(GEN Q, GEN ord, GEN ordinv, GEN e1, GEN e2){
  pari_sp top=avma;
  GEN v=qa_conjbasis(Q, ord, ordinv, e1, e2, 0);//We orient here to save computations
  if(gequal0(v)){avma=top;return gen_0;}
  GEN middle=qa_mul(Q, gel(v, 2), qa_conj(gel(v, 1)));//v2*conj(v1)
  int toneg=0;
  for(int i=2;i<=4;i++){//Testing whether we need to negate v2 or not.
	int j=gsigne(gel(e2, i));
	if(j==0) continue;
	if(j!=gsigne(gel(middle, i))) toneg=1;
	break;
  }
  long lx;
  GEN ret=cgetg(4, t_VEC);
  gel(ret, 1)=cgetg_copy(ret, &lx);
  gel(gel(ret, 1), 1)=qa_norm(Q, gel(v, 1));
  gel(gel(ret, 1), 2)=gmulgs(gel(middle, 1), 2);//v2*conj(v1)+v1*conj(v2)
  gel(gel(ret, 1), 3)=qa_norm(Q, gel(v, 2));
  gel(ret, 2)=gcopy(gel(v, 1));
  if(toneg==0) gel(ret, 3)=gcopy(gel(v, 2));
  else{
	togglesign_safe(&gel(gel(ret, 1), 2));//must negate this too
	gel(ret, 3)=gneg(gel(v, 2));
  }
  return gerepileupto(top, ret);
}

//qa_conjqf with typechecking and converting ord to an ord if it is an order.
GEN qa_conjqf_typecheck(GEN Q, GEN ord, GEN e1, GEN e2){
  pari_sp top=avma;
  qa_check(Q);qa_eltcheck(e1);qa_eltcheck(e2);
  GEN ordinv;
  if(typ(ord)==t_VEC && lg(ord)==QAORDLEN){ordinv=qa_getordinv(ord);ord=qa_getord(ord);}
  else{QM_check(ord);ordinv=ginv(ord);}
  return gerepileupto(top, qa_conjqf(Q, ord, ordinv, e1, e2));
}

//Checks if there is an element x of ord of norm n for which xe1x^(-1)=e2. Returns 1 if so, 0 else, and returns a possible such element if retconelt=1
GEN qa_conjnorm(GEN Q, GEN ord, GEN ordinv, GEN e1, GEN e2, GEN n, int retconelt, long prec){
  pari_sp top=avma;
  GEN data=qa_conjqf(Q, ord, ordinv, e1, e2);
  if(gequal0(data)){avma=top;return gen_0;}//Not conjugate
  GEN sols=bqf_reps(gel(data, 1), n, 0, 1, prec);
  if(gequal0(sols)){avma=top;return gen_0;}//No solutions
  if(retconelt==0){avma=top;return gen_1;}//Solution
  return gerepileupto(top, gadd(gmul(gel(data, 2), gel(gel(sols, 2), 1)), gmul(gel(data, 3), gel(gel(sols, 2), 2))));//Second entry is always a solution
}

//qa_conjnorm with typechecking
GEN qa_conjnorm_typecheck(GEN Q, GEN ord, GEN e1, GEN e2, GEN n, int retconelt, long prec){
  pari_sp top=avma;
  qa_check(Q);qa_eltcheck(e1);qa_eltcheck(e2);
  GEN ordinv;
  if(typ(ord)==t_VEC && lg(ord)==QAORDLEN){ordinv=qa_getordinv(ord);ord=qa_getord(ord);}
  else{QM_check(ord);ordinv=ginv(ord);}
  if(typ(n)!=t_INT) pari_err_TYPE("n must be an integer", n);
  return gerepileupto(top, qa_conjnorm(Q, ord, ordinv, e1, e2, n, retconelt, prec));
}

//If (e1, e2) and (f1, f2) have the same pairs of min polys and trace(e1e2)=trace(f1f2), the pairs are conjugate, with the conjugation space (union 0) being 1-dimensional (assuming e1, e2, and e1e2 are all not rational and e1!=e2). Intersecting with ord gives a Z-module, and this returns a generator.
GEN qa_simulconj(GEN Q, GEN ord, GEN ordinv, GEN e1, GEN e2, GEN f1, GEN f2, long prec){
  pari_sp top=avma;
  GEN lat1=qa_conjbasis(Q, ord, ordinv, e1, f1, 1);
  if(gequal0(lat1)){avma=top;return gen_0;}//Not actually conjugate
  GEN lat2=qa_conjbasis(Q, ord, ordinv, e2, f2, 1);
  if(gequal0(lat2)){avma=top;return gen_0;}//Not actually conjugate
  GEN mat1=cgetg(3, t_MAT), mat2=cgetg(3, t_MAT);//Will represent the conjugation spaces
  for(int i=1;i<=2;i++){gel(mat1, i)=gtocol(gel(lat1, i));gel(mat2, i)=gtocol(gel(lat2, i));}
  GEN lat=module_intersect(mat1, mat2);
  if(gequal0(lat)){avma=top;return gen_0;}//Not simultaneously conjugate (trace condition must have failed)
  return gerepileupto(top, gtovec(gel(lat, 1)));//it's a matrix with one column
}

//qa_conjnorm with typechecking
GEN qa_simulconj_typecheck(GEN Q, GEN ord, GEN e1, GEN e2, GEN f1, GEN f2, long prec){
  pari_sp top=avma;
  qa_check(Q);qa_eltcheck(e1);qa_eltcheck(e2);qa_eltcheck(f1);qa_eltcheck(f2);
  GEN ordinv;
  if(typ(ord)==t_VEC && lg(ord)==QAORDLEN){ordinv=qa_getordinv(ord);ord=qa_getord(ord);}
  else{QM_check(ord);ordinv=ginv(ord);}
  return gerepileupto(top, qa_simulconj(Q, ord, ordinv, e1, e2, f1, f2, prec));
}



//EMBEDDING QUADRATIC ORDERS INTO EICHLER ORDERS



//Given a quaternion algebra Q, order order, and embedding of O_D given by emb=image of (p_D+sqrt(D))/2, this outputs the pair [em', D'] which is the embedding associated to phi and ord. It is the unique optimal embedding of an order O_D' which agrees with em on the fraction field. D can be passed in as 0.
GEN qa_associatedemb(GEN Q, GEN order, GEN emb, GEN D){
  pari_sp top=avma;
  GEN ord=qa_getord(order), ordinv=qa_getordinv(order);
  GEN rootD=cgetg(5, t_COL);
  gel(rootD, 1)=gen_0;
  for(int i=2;i<=4;i++) gel(rootD, i)=gmulgs(gel(emb, i), 2);//rootD represent the image of sqrt(D) as a column vector
  if(gequal0(D)) D=gel(qa_square(Q, rootD), 1);
  GEN v=gmul(ordinv, rootD);//sqrt(D) in terms of order basis
  GEN l=Q_denom(v);
  v=gmul(v, l);
  D=gmul(D, sqri(l));//Updating v and D. Now we land inside the order guarenteed.
  l=Q_content(v);
  v=ZV_Z_divexact(v, l);
  D=diviiexact(D, sqri(l));//Updating v and D to the smallest D landing in the order. We now need to check that adding D%2 and dividing by 2 is okay.
  GEN newemb=gdivgs(gtovec(gmul(ord, v)), 2);
  if(smodis(D, 2)==1) gel(newemb, 1)=ghalf;//The embedding in terms of i, j, k
  if(qa_isinorder(Q, ordinv, newemb)) return gerepileupto(top, mkvec2copy(newemb, D));//Good!
  GEN rvec=cgetg(3, t_VEC);//Now we must double the embedding and quadruple D
  gel(rvec, 1)=cgetg(5, t_VEC);
  gel(gel(rvec, 1), 1)=gen_0;
  for(int i=2;i<=4;i++) gel(gel(rvec, 1), i)=gmulgs(gel(newemb, i), 2);
  gel(rvec, 2)=shifti(D, 2);//4D
  return gerepileupto(top, rvec);
}

//qa_associatedemb with typecheck and presetting D
GEN qa_associatedemb_typecheck(GEN Q, GEN order, GEN emb, GEN D){
  pari_sp top=avma;
  qa_check(Q);
  qa_ordcheck(order);
  qa_eltcheck(emb);
  if(typ(D)!=t_INT) pari_err_TYPE("D is either 0 or the discriminant of emb", D);
  return gerepileupto(top, qa_associatedemb(Q, order, emb, D));
}

//Finds nembeds optimal embedding of the quadratic order into the quaternion order. Does NOT check that we can have this many distinct embeddings, so method will not terminate if this is so.
GEN qa_embed(GEN Q, GEN order, GEN D, GEN nembeds, GEN rpell, long prec){
  pari_sp top=avma;
  long pos=1, longnembeds=itos(nembeds);//Current position of solution, number of embeddings in long format
  GEN embs=cgetg(longnembeds+1, t_VEC);
  GEN dfacs=discprimeindex(D, gen_0);//Prime divisors of D for which D/p^2 is a discriminant. Used to check embedding is optimal
  GEN Dmod2=gmodgs(D, 2);//D modulo 2; we now try to embed (Dmod2+sqrt(D))/2.
  //We search for aF^2+bG^2-abH^2=D and (Dmod2+Fi+Gj+Hk)/2 is optimal in order. F, G, H don't need to be integers, but they do lie in a lattice, so we start by correcting for this (translating them to integers).
  GEN ab=qa_getabvec(Q), maxds=qa_getordmaxd(order), ord=qa_getord(order), ordinv=qa_getordinv(order);//[a, b, -ab], Maximal denominators, order, inverse order
  GEN di=Qdivii(gel(ab, 1), sqri(gel(maxds, 2))), dj=Qdivii(gel(ab, 2), sqri(gel(maxds, 3))), dk=Qdivii(gel(ab, 3), sqri(gel(maxds, 4)));
  GEN gc=lcmii(Q_denom(di), lcmii(Q_denom(dj), Q_denom(dk)));//By convention, gcd(a, b)=1 hence one of them is odd, hence 4 divides gc necessarily (as 2 divides maxds).
  GEN E=Q_muli_to_int(di, gc), F=Q_muli_to_int(dj, gc), G=Q_muli_to_int(dk, gc), res=shifti(mulii(gc, D), -2);//Now Ef1^2+Fg1^2+Gh1^2=res, where E, F, G, f1, g1, h1 are integers and gcd(E, F, G)=1. Note that we are currently ignoring optimality.
  if(signe(E)==-1){togglesign_safe(&E);togglesign_safe(&F);togglesign_safe(&G);togglesign_safe(&res);}//Making sure we work with positive definite forms.
  GEN testvec=mkvec4(Qdivii(Dmod2, gen_2), gen_0, gen_0, gen_0);//The image of (D mod 2+sqrt(D))/2.
  GEN sols, f, g, h, mf, mg, mh, gcdf1g1h1;
  if(signe(gel(ab, 1))==signe(gel(ab, 2))){//a, b have the same sign. We fix h1, and solve for the definite BQF in f1, g1, which is definite
    GEN h1=gen_0, twoh1p1G=G, twoG=shifti(G, 1);
	GEN form=mkvec3(E, gen_0, F);//Since E|a and F|b, this is a primitive form.
	GEN formred=dbqf_red_tmat(form), formdisc=bqf_disc(form);//Reduction of the form and discriminant.
	for(;;){//Solving Ef1^2+Fg1^2=res=gc*D/4-G*h1^2, starting with h1=0 and going to oo.
	  sols=dbqf_reps(formred, formdisc, res, 0, 1);//The solutions to the DBQF, and we only pick up half of them (up to +/- sign)
	  if(!gequal0(sols)){//Solutions!
		for(long i=2;i<lg(sols);i++){//For each solution, we have 4 sign cases to test: (s1*f1, s1*g1, s2*h1) with s1, s2=-1 or 1.
		f=Qdivii(gel(gel(sols, i), 1), gel(maxds, 2));mf=gneg(f);
		g=Qdivii(gel(gel(sols, i), 2), gel(maxds, 3));mg=gneg(g);
		h=Qdivii(h1, gel(maxds, 4));mh=gneg(h);
		gcdf1g1h1=gcdii(gcdii(gel(gel(sols, i), 1), gel(gel(sols, i), 2)), h1);//gcd(f1, g1, h1)
		for(int caseno=1;caseno<=4;caseno++){//Doing the signs
		  switch(caseno){
		    case 1://++
			  gel(testvec, 2)=f;gel(testvec, 3)=g;gel(testvec, 4)=h;break;
			case 2://-+
			  gel(testvec, 2)=mf;gel(testvec, 3)=mg;break;
			case 3://--
			  gel(testvec, 4)=mh;break;
			case 4://+-
			  gel(testvec, 2)=f;gel(testvec, 3)=g;break;
		    }
		  if(!qa_embed_isnewoptimal(Q, ord, ordinv, D, Dmod2, dfacs, testvec, gcdf1g1h1, embs, pos, prec)) continue;//Not new
		  gel(embs, pos)=gcopy(testvec);
		  pos++;
		  if(pos>longnembeds) goto FOUNDSOLUTIONS;//Done
		  }
		}
	  }
	  res=subii(res, twoh1p1G);//Updating res to res-G(2h_1+1)
	  h1=addis(h1, 1);//Updating h_1
	  twoh1p1G=addii(twoh1p1G, twoG);//Updating G(2h_1+1)
	}
  }
  else{//a, b have opposite sign. In this case, we iterate over g1, and solve for f1, h1.
    GEN EGgcd=gcdii(E, G);
	E=diviiexact(E, EGgcd);
	G=diviiexact(G, EGgcd);
	res=Qdivii(res, EGgcd);//res may be rational now.
	GEN g1=gen_0, twog1p1FdEGgcd=Qdivii(F, EGgcd), twoFdEGgcd=Qdivii(shifti(F, 1), EGgcd);
	GEN form=mkvec3(E, gen_0, G);
	GEN formred=dbqf_red_tmat(form), formdisc=bqf_disc(form);//Reduction of the form and discriminant
	for(;;){//Solving Ef1^2+Gh1^2=res=(gc*D/4-F*g1^2)/EGgcd, starting with g1=0 and going to oo.
	  if(typ(res)!=t_INT){
		res=gsub(res, twog1p1FdEGgcd);//Updating res to res-F(2g_1+1)/EGgcd
	    g1=addis(g1, 1);//Updating g_1
	    twog1p1FdEGgcd=gadd(twog1p1FdEGgcd, twoFdEGgcd);//Updating F(2g_1+1)/EGgcd
	    continue;
	  }
	  sols=dbqf_reps(formred, formdisc, res, 0, 1);//The solutions to the DBQF, and we only pick up half of them (up to +/- sign)
	  if(!gequal0(sols)){//Solutions!
		for(long i=2;i<lg(sols);i++){//For each solution, we have 4 sign cases to test: (s1*f1, s2*g1, s1*h1) with s1, s2=-1 or 1.
		f=Qdivii(gel(gel(sols, i), 1), gel(maxds, 2));mf=gneg(f);
		g=Qdivii(g1, gel(maxds, 3));mg=gneg(g);
		h=Qdivii(gel(gel(sols, i), 2), gel(maxds, 4));mh=gneg(h);
		gcdf1g1h1=gcdii(gcdii(gel(gel(sols, i), 1), gel(gel(sols, i), 2)), g1);//gcd(f1, g1, h1)
		for(int caseno=1;caseno<=4;caseno++){//Doing the signs
		  switch(caseno){
		    case 1://++
			  gel(testvec, 2)=f;gel(testvec, 3)=g;gel(testvec, 4)=h;break;
			case 2://-+
			  gel(testvec, 2)=mf;gel(testvec, 4)=mh;break;
			case 3://--
			  gel(testvec, 3)=mg;break;
			case 4://+-
			  gel(testvec, 2)=f;gel(testvec, 4)=h;break;
		    }
		  if(!qa_embed_isnewoptimal(Q, ord, ordinv, D, Dmod2, dfacs, testvec, gcdf1g1h1, embs, pos, prec)) continue;//Not new
		  gel(embs, pos)=gcopy(testvec);
		  pos++;
		  if(pos>longnembeds) goto FOUNDSOLUTIONS;//Done
		  }
		}
	  }
	  res=gsub(res, twog1p1FdEGgcd);//Updating res to res-F(2g_1+1)/EGgcd
	  g1=addis(g1, 1);//Updating g_1
	  twog1p1FdEGgcd=gadd(twog1p1FdEGgcd, twoFdEGgcd);//Updating F(2g_1+1)/EGgcd
	}
  }
  FOUNDSOLUTIONS:
  if(gequal0(rpell)) return gerepilecopy(top, embs);//Just return the embeddings
  //Now, rpell is assumed to  be [T, U]=pell(D).
  GEN addv;
  if(gequal0(Dmod2)) addv=mkvec4(shifti(gel(rpell, 1), -1), gen_0, gen_0, gen_0);//D even, embed sqrt(D)/2. addv=[T/2, 0, 0, 0]
  else addv=mkvec4(shifti(subii(gel(rpell, 1), gel(rpell, 2)), -1), gen_0, gen_0, gen_0);//D odd, embed (1+sqrt(D))/2. addv=[(T-U)/2, 0, 0, 0]
  GEN Uvec=cgetg_copy(embs, &pos);
  for(long i=1;i<pos;i++) gel(Uvec, i)=gmul(gel(embs, i), gel(rpell, 2));//U*embeddings
  GEN retvec=cgetg_copy(embs, &pos);
  for(long i=1;i<pos;i++) gel(retvec, i)=gadd(gel(Uvec, i), addv);
  return gerepileupto(top, retvec);
}

//qa_embed with typecheck and setting of nembeds if required. If nembeds=0, we require Q to be indefinite and order to be Eichler.
GEN qa_embed_typecheck(GEN Q, GEN order, GEN D, GEN nembeds, int rpell, long prec){
  pari_sp top=avma;
  if(!isdisc(D)) return gen_0;
  if(gequal0(nembeds)){
    qa_indefcheck(Q);qa_ordeichlercheck(order);
	nembeds=gel(qa_numemb(Q, order, D, gel(bqf_ncgp(D, prec), 1)), 1);
	if(gequal0(nembeds)){avma=top;return gen_0;}
  }
  else{qa_check(Q);qa_ordcheck(order);}
  GEN pel;
  if(rpell==1) pel=pell(D);
  else pel=gen_0;
  return gerepileupto(top, qa_embed(Q, order, D, nembeds, pel, prec));
}

//Checks if we have a new and optimal embedding.
static int qa_embed_isnewoptimal(GEN Q, GEN ord, GEN ordinv, GEN D, GEN Dmod2, GEN dfacs, GEN emb, GEN gcdf1g1h1, GEN embs, long pos, long prec){
  pari_sp top=avma;
  if(!qa_isinorder(Q, ordinv, emb)) return 0;//Not in order
  GEN pvec=zerovec(4), p;
  for(long i=1;i<lg(dfacs);i++){
	p=gel(dfacs, i);
	if(!dvdii(gcdf1g1h1, p)) continue;//Cannot be continued to D/p^2
	if(equalii(p, gen_2) && smodis(D, 8)==4) gel(pvec, 1)=ghalf;
    else gel(pvec, 1)=Qdivii(Dmod2, gen_2);//Setting first coeff.
	gel(pvec, 2)=gdiv(gel(emb, 2), p);
	gel(pvec, 3)=gdiv(gel(emb, 3), p);
	gel(pvec, 4)=gdiv(gel(emb, 4), p);//Now pvec is the embedding of D/p^2.
	if(qa_isinorder(Q, ordinv, pvec)){avma=top;return 0;}//Not optimal wrt p
  }//Now we are optimal, so we check if we are also new.
  for(long i=1;i<pos;i++) if(!gequal0(qa_conjnorm(Q, ord, ordinv, emb, gel(embs, i), gen_1, 0, prec))){avma=top;return 0;}//Testing if new
  avma=top;
  return 1;//Passed!
}

//Outputs the list of discriminants which have optimal embeddings into an Eichler order. Note that we call qa_numemb, so if you want that info calling this first is redundant.
GEN qa_embeddablediscs(GEN Q, GEN order, GEN d1, GEN d2, int fund, GEN cop){
  pari_sp top=avma;
  GEN baselist=disclist(d1, d2, fund, cop);
  long lblist=lg(baselist);
  GEN l=vectrunc_init(lblist), N;//The list
  for(long i=1;i<lblist;i++){
	N=qa_numemb(Q, order, gel(baselist, i), gen_1);
	if(gequal0(N)){cgiv(N);continue;}
	cgiv(N);
	vectrunc_append(l, icopy(gel(baselist, i)));
  }
  return gerepileupto(top, l);
}

//qa_embeddablediscs with typecheck
GEN qa_embeddablediscs_typecheck(GEN Q, GEN order, GEN d1, GEN d2, int fund, GEN cop){
  qa_indefcheck(Q);qa_ordeichlercheck(order);
  return qa_embeddablediscs(Q, order, d1, d2, fund, cop);
}

//Returns the data associated to the number of optimal embeddigns of O_D into order. The format is: [m, n, v1, v2, v3]: m is the total number of opt embs, n=h^+(D) is the number of a fixed orientation (NOTE: we just copy what is given in narclno, so if this number is wrong, m and n will be wrong. However, this can be useful if you just want to know if it's non-zero or not, so passing in 1 is OK for those purposes). v1=[x] with x being the number of embeddings at oo (1 if D>0 and 2 if D<0). v2=[x_1,...,x_r] with Q[1]=[p_1,...,p_r] and x_i=number of local embeddings at p_i for i<=r (1 if p_i divides D and 2 else). v3=[y_1,...,y_s]. where s distinct primes q_1,...,q_s divide the level of the order and y_i is the number of local embedding classes. This is 2 as long as q_i does not divide D.
GEN qa_numemb(GEN Q, GEN order, GEN D, GEN narclno){
  pari_sp top=avma;
  GEN dfund=fdisc(D), dratio=diviiexact(D, dfund), m=narclno;//Ratio of D to D_fund
  GEN qapram=qa_getpram(Q), ordpram=qa_getordlevelpfac(order), p;
  long lrams, llevelrams, kron;
  GEN rams=cgetg_copy(qapram, &lrams);//Local orientations at ramifying primes
  for(long i=1;i<lrams;i++){//We need kronecker(dfund, p)!=1 for all primes ramifying AND p does not divide dratio.
    p=gel(qapram, i);
	if(dvdii(dratio, p)){avma=top;return zerovec(4);}//Failed optimality
	kron=kronecker(dfund, p);
	switch(kron){
      case 1:
	    avma=top;return zerovec(4);//kronecker=1, no embedding
	  case -1:
	    gel(rams, i)=gen_2;m=shifti(m, 1);break;
	  case 0:
	    gel(rams, i)=gen_1;break;
	}
  }//On to the level ramifications
  GEN levelrams=cgetg_copy(ordpram, &llevelrams), Dover4=gdivgs(D, 4), infor;
  for(long i=1;i<llevelrams;i++){//Need kronecker(D, p)!=-1. The exact count is a little more complicated: if p does not divide D it is 2, otherwise it is the sum of solutions to x^2==D/4 modulo p^e PLUS the number of orbits modulo p^e of solutions modulo p^{e+1} (where p^e||level). See Voight 30.6.12.
    p=gel(gel(ordpram, i), 1);//The prime
	kron=kronecker(D, p);
	if(kron==-1){avma=top;return zerovec(4);}//Failed
	else if(kron==1){gel(levelrams, i)=gen_2;m=shifti(m, 1);continue;}//kronecker=1, so 2 orientations
	GEN e=gel(gel(ordpram, i), 2);//The exponent
	GEN p2e=powii(p, e);
	GEN fact= gtomat(gel(ordpram, i));
	GEN roots1=sqmod(Dover4, p2e, fact);//Solutions modulo p^e
	if(gequal0(roots1)){avma=top;return zerovec(4);}//Can't do it
	gcoeff(fact, 1, 2)=addis(gcoeff(fact, 1, 2), 1);
	GEN roots2=sqmod(Dover4, mulii(p2e, p), fact);//Solutions modulo p^{e+1}
	GEN nors=mulii(stoi(lg(gel(roots1, 1))-1), diviiexact(p2e, gel(roots1, 2)));//The number of solutions modulo p^e
	if(!gequal0(roots2)){//More to do!
	  if(cmpii(gel(roots2, 2), p2e)==1){
		GEN v=FpV_red(gel(roots2, 1), p2e);//Reduction modulo p^e
		v=ZV_sort_uniq(v);//Removing equal elements
		nors=addis(nors, lg(v)-1);//Adding in
	  }
	  else nors=addii(nors, mulsi(lg(gel(roots2, 1))-1, diviiexact(p2e, gel(roots2, 2))));//Solutions are already distinct modulo p^e, so need to just multiply to boost up.
	}
	gel(levelrams, i)=nors;
	m=mulii(m, nors);
  }
  if(signe(D)==1) infor=mkvecs(1);//Infinite orientation
  else{infor=mkvecs(2);m=shifti(m, 1);}
  GEN rvec=cgetg(6, t_VEC);
  gel(rvec, 1)=icopy(m);
  gel(rvec, 2)=icopy(narclno);
  gel(rvec, 3)=ZV_copy(infor);
  gel(rvec, 4)=ZV_copy(rams);
  gel(rvec, 5)=ZV_copy(levelrams);
  return gerepileupto(top, rvec);
}

//qa_numemb with typecheck and presetting of narclno if passed as 0.
GEN qa_numemb_typecheck(GEN Q, GEN order, GEN D, GEN narclno, long prec){
  pari_sp top=avma;
  if(!isdisc(D)) return gen_0;
  qa_indefcheck(Q);qa_ordeichlercheck(order);
  if(typ(narclno)!=t_INT) pari_err_TYPE("Please pass in narclno as 0 or the narrow class number of discriminant D", narclno);
  if(gequal0(narclno)) narclno=gel(bqf_ncgp(D, prec), 1);
  return gerepileupto(top, qa_numemb(Q, order, D, narclno));
}

//Computes the set of primes where the orientations of the embeddings e1 and e2 differ.
GEN qa_ordiffer(GEN Q, GEN order, GEN e1, GEN e2, GEN D){
  pari_sp top=avma;
  GEN ord=qa_getord(order), ordinv=qa_getordinv(order), qaram=qa_getpram(Q), levelram=qa_getordlevelpfac(order);
  GEN form=gel(qa_conjqf(Q, ord, ordinv, e1, e2), 1);
  GEN g=Z_content(form);
  if(g==NULL){//Same orientation at the finite places!
    if(signe(D)==1 || signe(gel(form, 1))==1){avma=top;return cgetg(1, t_VEC);}//D>0 or form is pos def
	avma=top;//D<0 and form is neg def
	GEN rvec=cgetg(2, t_VEC);
	gel(rvec, 1)=mkoo();
	return rvec;
  }
  GEN l=vectrunc_init(lg(qaram)+lg(levelram)), r;
  long e;
  for(long i=1;i<lg(qaram);i++){//Primes ramifying in Q
	e=Z_pvalrem(g, gel(qaram, i), &r);
	if(e==0) continue;//Did not divide
	g=icopy(r);//Must copy as we use its address in the next call
	vectrunc_append(l, gel(qaram, i));//We sort at the end so it's OK to append the original place
	if(equali1(g)) break;
  }
  if(!equali1(g)){
	for(long i=1;i<lg(levelram);i++){//Primes dividing the level
	  e=Z_pvalrem(g, gel(gel(levelram, i), 1), &r);
	  if(e==0) continue;//Did not divide
	  g=icopy(r);//Must copy as we use its address in the next call
	  vectrunc_append(l, gel(gel(levelram, i), 1));//We sort at the end so it's OK to append the original place
	  if(equali1(g)) break;
	}
  }
  if(signe(D)==-1 && signe(gel(form, 1))==-1) vectrunc_append(l, mkoo());//Neg definite
  return gerepileupto(top, sort(l));
}

//qa_ordiffer with typecheck
GEN qa_ordiffer_typecheck(GEN Q, GEN order, GEN e1, GEN e2, GEN D){
  pari_sp top=avma;
  qa_indefcheck(Q);qa_ordeichlercheck(order);
  qa_eltcheck(e1);qa_eltcheck(e2);
  if(gequal0(D)){
	GEN v=gmulgs(e1, 2);
	gel(v, 1)=gen_0;
	D=gel(qa_square(Q, v), 1);
  }
  else if(!isdisc(D)) pari_err_TYPE("D is not a discriminant", D);
  return gerepileupto(top, qa_ordiffer(Q, order, e1, e2, D));
}

//Returns the orientation at oo
GEN qa_orinfinite(GEN Q, GEN emb, GEN D, long prec){
  pari_sp top=avma;
  if(signe(D)==1) return gen_1;
  GEN a=qa_geta(Q);
  if(signe(a)==-1) pari_err_TYPE("Sorry, at the moment this is not initialized for a<0.", a);
  if(gcmp(gel(emb, 3), gmul(gsqrt(a, prec), gel(emb, 4)))==1){avma=top;return gen_1;}
  avma=top;
  return gen_m1;
}

//qa_orinfinite with typecheck
GEN qa_orinfinite_typecheck(GEN Q, GEN emb, GEN D, long prec){
  pari_sp top=avma;
  qa_indefcheck(Q);
  qa_eltcheck(emb);
  if(gequal0(D)){
	GEN v=gmulgs(emb, 2);
	gel(v, 1)=gen_0;
	D=gel(qa_square(Q, v), 1);
  }
  else if(!isdisc(D)) pari_err_TYPE("D is not a discriminant", D);
  return gerepileupto(top, qa_orinfinite(Q, emb, D, prec));
}

//This method finds the optimal embeddings of O_D into order, and sorts them into a matrix. A row represents an orientation, with the first entry being the orientation, and the other h^+(D) entries being the embeddings. They are organized so as to respect the group action, and the rows are organized consistently with each other (mostly related by "Atkin-Lehner involutions")
GEN qa_sortedembed(GEN Q, GEN order, GEN D, GEN rpell, GEN ncgp, long prec){
  pari_sp top=avma, mid;
  GEN nemb=qa_numemb(Q, order, D, gel(ncgp, 1));//The number of embeddings
  GEN ord=qa_getord(order), ordinv=qa_getordinv(order);//basic info
  if(gequal0(gel(nemb, 1))){avma=top;return gen_0;}//No embeddings
  GEN embs=qa_embed(Q, order, D, gel(nemb, 1), rpell, prec);//The embeddings
  long lenembs;
  mid=avma;
  GEN emborpairs=cgetg_copy(embs, &lenembs);
  for(long i=1;i<lenembs;i++) gel(emborpairs, i)=mkvec2(gel(embs, i), qa_ordiffer(Q, order, gel(embs, 1), gel(embs, i), D));//Pairing the embeddings with orientations
  emborpairs=gerepileupto(mid, gen_sort(emborpairs, NULL, &qa_embedor_compare));//Initial sorting
  long nor=itos(diviiexact(gel(nemb, 1), gel(nemb, 2)));//The number of orientations
  GEN orvecs=cgetg(nor+1, t_VEC), embvecs=cgetg(nor+1, t_VEC);//For storing the orientations and embeddings
  long baseorind=1, endorind, hplusD=itos(gel(ncgp, 1)), ind;//The base or at the start of each section
  GEN tempform, g, rootD=gsqrt(D, prec), temp;
  for(long or=1;or<=nor;or++){
	endorind=baseorind;
	do{endorind++;}while(endorind<lenembs && gequal(gel(gel(emborpairs, baseorind), 2), gel(gel(emborpairs, endorind), 2)));
	endorind--;//Now the block from baseorind to endorind inclusive constitutes a single orientation difference. This will be one row of the matrix, EXCEPT if there are primes dividing the level that also divide D, in which case it may comprise multiple.
	if(endorind-baseorind+1==hplusD){//Single, so all one line
	  gel(orvecs, or)=gel(gel(emborpairs, baseorind), 2);//The orientation of the block
	  gel(embvecs, or)=cgetg(hplusD+1, t_VEC);
	  for(long i=baseorind;i<=endorind;i++){
		tempform=gel(qa_conjqf(Q, ord, ordinv, gel(gel(emborpairs, 1), 1), gel(gel(emborpairs, i), 1)), 1);
		g=Z_content(tempform);//This is a little redundant as we did this on ordiffer, but it's what we do for now
		if(g==NULL) g=gen_1;
		tempform=ZV_Z_divexact(tempform, g);
		if(signe(D)==-1 && signe(gel(tempform, 1))==-1){gel(tempform, 1)=gneg(gel(tempform, 1));gel(tempform, 3)=gneg(gel(tempform, 3));}//Making pos def and primitive. We have to not negate the middle coefficient (<==> taking the inverse) since this is what is required to make the action correct (we are computing from the base of the WHOLE chain, not from each row, and at oo the conjugation by -1 introduces a "shift" which is exactly this.)
		ind=itos(bqf_isequiv_set(tempform, gel(ncgp, 3), rootD, signe(D), 0));
		gel(gel(embvecs, or), ind)=gel(gel(emborpairs, i), 1);
	  }
	}
	else{
	  while(baseorind<=endorind){
		gel(orvecs, or)=gel(gel(emborpairs, baseorind), 2);
		gel(embvecs, or)=cgetg(hplusD+1, t_VEC);
		long foundforms=1;
		gel(gel(embvecs, or), 1)=gel(gel(emborpairs, baseorind), 1);
		for(long i=baseorind+1;i<=endorind;i++){
		  tempform=gel(qa_conjqf(Q, ord, ordinv, gel(gel(emborpairs, baseorind), 1), gel(gel(emborpairs, i), 1)), 1);
		  g=Z_content(tempform);
		  if(g==NULL){//Found one!
		    ind=itos(bqf_isequiv_set(tempform, gel(ncgp, 3), rootD, signe(D), 0));
		    gel(gel(embvecs, or), ind)=gel(gel(emborpairs, i), 1);
			temp=gel(emborpairs, baseorind+foundforms);
			gel(emborpairs, baseorind+foundforms)=gel(emborpairs, i);
			gel(emborpairs, i)=temp;//Swaperoo
			foundforms++;
			if(foundforms==hplusD) break;//Done the block
		  }
		}
		baseorind=baseorind+hplusD;
		or++;
	  }
      or--;//or++ will occur when the for loop hits again
	}
	baseorind=endorind+1;
  }
  GEN ret=cgetg(3, t_MAT);//Putting it all together
  gel(ret, 1)=gtocol(orvecs);
  gel(ret, 2)=gtocol(embvecs);
  return gerepileupto(top, ret);
}

//qa_sortedembed with typecheck
GEN qa_sortedembed_typecheck(GEN Q, GEN order, GEN D, int rpell, GEN ncgp, long prec){
  pari_sp top=avma;
  if(!isdisc(D)) return gen_0;
  qa_indefcheck(Q);qa_ordeichlercheck(order);
  if(gequal0(ncgp)) ncgp=bqf_ncgp_lexic(D, prec);  
  GEN pel;
  if(rpell==1) pel=pell(D);
  else pel=gen_0;
  return gerepileupto(top, qa_sortedembed(Q, order, D, pel, ncgp, prec));
}

//Given two pairs [emb, or], this sorts by or (smallest to longest length first, then lexicographically
static int qa_embedor_compare(void *data, GEN pair1, GEN pair2){
  long l1=lg(gel(pair1, 2)), l2=lg(gel(pair2, 2));
  if(l1<l2) return -1;
  else if(l1>l2) return 1;
  return lexcmp(gel(pair1, 2), gel(pair2, 2));
}



//CHECKING METHODS



//Checks the format of a quaternion algebra is correct
void qa_check(GEN Q){
  if(typ(Q)!=t_VEC || lg(Q)!=QALEN) pari_err_TYPE("Not a quaternion algebra", Q);
  GEN abvec=qa_getabvec(Q);
  if(typ(abvec)!=t_VEC || lg(abvec)!=4) pari_err_TYPE("Not a quaternion algebra", Q);
}

//Checks that Q is an indefinite quaternion algebra
void qa_indefcheck(GEN Q){
  if(typ(Q)!=t_VEC || lg(Q)!=QALEN) pari_err_TYPE("Not a quaternion algebra", Q);
  GEN abvec=qa_getabvec(Q);
  if(typ(abvec)!=t_VEC || lg(abvec)!=4) pari_err_TYPE("Not a quaternion algebra", Q);
  if(signe(gel(abvec, 1))==-1 && signe(gel(abvec, 2))==-1) pari_err_TYPE("Not an indefinite quaternion algebra", Q);
}

//Checks that x is of a valid form for a quaternion algebra element
void qa_eltcheck(GEN x){
  if(typ(x)!=t_VEC || lg(x)!=5) pari_err_TYPE("Not an element of a quaternion algebra", x);
}

//Checks ord is an order or an initialized order, and returns the matrix whose column vectors generate it. NOT STACK CLEAN, returns a reference to the original object.
GEN qa_ordcheck(GEN ord){
  GEN actualord;
  if(typ(ord)==t_VEC && lg(ord)==QAORDLEN) actualord=qa_getord(ord);
  else actualord=ord;
  QM_check(actualord);
  return actualord;
}

//Checks that order in an initialized Eichler order
void qa_ordeichlercheck(GEN order){
  if(typ(order)!=t_VEC || lg(order)!=QAORDLEN) pari_err_TYPE("Not an initialized order", order);
  if(signe(qa_getordtype(order))==-1) pari_err_TYPE("Not an Eichler order", order);
}

//Checks if M is a QM
void QM_check(GEN M){
  if(typ(M)!=t_MAT) pari_err_TYPE("Please input a QM matrix", M);
  long t;
  for(long i=1;i<lg(M);i++){
	for(long j=1;j<lg(gel(M, 1));j++){
	  t=typ(gcoeff(M, j, i));
	  if(t!=t_FRAC && t!=t_INT) pari_err_TYPE("Please input a QM matrix", M);
	}
  }
}



//PROPERTY RETRIEVAL



//These methods are generally not gerepile safe, as they retrieve a member.

GEN qa_getnf(GEN Q){return gel(Q, 1);}//From Q, retrieve the number field
GEN qa_getpram(GEN Q){return gel(Q, 2);}//From Q, retrieve the ramifying primes
GEN qa_getabvec(GEN Q){return gel(Q, 3);}//From Q, retrieve the vector [a, b, -ab]
GEN qa_geta(GEN Q){return gel(gel(Q, 3), 1);}//From Q, retrieve a
GEN qa_getb(GEN Q){return gel(gel(Q, 3), 2);}//From Q, retrieve b
GEN qa_getmab(GEN Q){return gel(gel(Q, 3), 3);}//From Q, retrieve -ab
GEN qa_getpramprod(GEN Q){return gel(Q, 4);}//From Q, retrieve the product of the ramifying primes

GEN qa_getord(GEN order){return gel(order, 1);}//From order, get ord
GEN qa_getordtype(GEN order){return gel(order, 2);}//From the order, get the type
GEN qa_getordmaxd(GEN order){return gel(order, 3);}//From the order, get the max denominators of the coefficients
GEN qa_getordlevel(GEN order){return gel(order, 4);}//From order, get the level
GEN qa_getordlevelpfac(GEN order){return gel(order, 5);}//From order, get the prime factorization of the level
GEN qa_getordinv(GEN order){return gel(order, 6);}//From order, get the inverse of ord.
GEN qa_getordtrace0basis(GEN order){return gel(order, 7);}//From order, get the basis of trace 0 elements



//SUPPORTING METHODS



//Useful for having gcmp in gen_sort methods (e.g. for gen_sort_uniq).
int cmp_data(void *data, GEN x, GEN y){return gcmp(x, y);}

//Given Z-modules spanned by the colums of A, B (QM matrices), this finds their intersection
GEN module_intersect(GEN A, GEN B){
  pari_sp top=avma;
  GEN combine=shallowconcat(A, B);//Shallowly concatenating A, B
  GEN inter=matkerint0(combine, 0);
  if(gequal0(inter)){avma=top;return gen_0;}//No intersection
  GEN M=QM_mul(A, matslice(inter, 1, lg(A)-1, 1, lg(inter)-1));
  return gerepileupto(top, QM_hnf(M));
}

//module_intersect with typechecking
GEN module_intersect_typecheck(GEN A, GEN B){
  QM_check(A);QM_check(B);return module_intersect(A, B);
}

//relations=[[p_1,s_1],...,[p_k,s_k]] with the p_i distinct integers and s_i=-1 or 1 for each i, and extra=0 or [n, c]. This searches for a prime p such that p==c mod n for all i AND kronecker(p,p_i)=s_i for all i. If there is no solution, this will not terminate.
GEN prime_ksearch(GEN relations, GEN extra){
  pari_sp top=avma;
  forprime_t T;
  GEN p;
  forprime_init(&T, gen_1, NULL);//Initialize prime search.
  long rellen=lg(relations);
  int isgood;
  if(gequal0(extra)){
	for(;;){
	  p=forprime_next(&T);
	  isgood=1;
	  for(long i=1;i<rellen;i++){
		if(equalsi(kronecker(p, gel(gel(relations, i),1)), gel(gel(relations, i), 2))) continue;
		isgood=0;//If we get here, they are not equal.
		break;
	  }
	  if(isgood) return gerepilecopy(top, p);
	  if(gc_needed(top, 2)) pari_err(e_MISC, "Computation exceeded 2/3 of available memory. Either re-run with more memory, or the givens are inconsistent (which is far more likely)");
	}
  }//If not, we have a congruence to also check.
  for(;;){
	p=forprime_next(&T);
	if(!equalii(modii(p, gel(extra, 1)), gel(extra, 2))) continue;//Must get the modulus constraint in first.
	isgood=1;
	for(long i=1;i<rellen;i++){
      if(equalsi(kronecker(p, gel(gel(relations, i),1)), gel(gel(relations, i), 2))) continue;
	  isgood=0;//If we get here, they are not equal.
	  break;
	}
	if(isgood) return gerepilecopy(top, p);
	if(gc_needed(top, 2)) pari_err(e_MISC, "Computation exceeded 2/3 of available memory. Either re-run with more memory, or the givens are inconsistent (which is far more likely)");
  }
}

//prime_ksearch with typechecking and a warning if the method may not terminate.
GEN prime_ksearch_typecheck(GEN relations, GEN extra){
  pari_sp top=avma;
  if(typ(relations)!=t_VEC || lg(relations)==1) pari_err_TYPE("Please enter a vector of relations of the form [p, +/-1]", relations);
  for(long i=1;i<lg(relations);i++){
	if(typ(gel(relations, i))!=t_VEC || lg(gel(relations, i))!=3) pari_err_TYPE("Each relation must be of the form [p, +/-1]", gel(relations, i));
	if(typ(gel(gel(relations, i), 1))!=t_INT || (gequal(gel(gel(relations, i), 2), gen_m1)==0 && gequal(gel(gel(relations, i), 2), gen_1)==0)) pari_err_TYPE("Each relation must be of the form [p, +/-1]", gel(relations, i));
	if(Z_issquare(gel(gel(relations, i), 1)) && equalii(gel(gel(relations, i), 2), gen_m1)) pari_err_TYPE("kronecker(p, n^2)!=-1, do not enter a square!", gel(relations, i));
  }
  if(!gequal0(extra)){
	if(typ(extra)!=t_VEC || lg(extra)!=3 || typ(gel(extra, 1))!=t_INT || typ(gel(extra, 2))!=t_INT || gequal0(gel(extra, 1))) pari_err_TYPE("Please enter 0 or [n, c] with n>0 and 0<=c<n", extra);
	GEN newextra=cgetg(3, t_VEC);
	gel(newextra, 1)=absi(gel(extra, 1));
	gel(newextra, 2)=modii(gel(extra, 2), gel(newextra, 1));
	if(!equali1(gcdii(gel(newextra, 1), gel(newextra, 2)))) pari_err_TYPE("Relation must be 0 or a coprime pair", extra);
	return gerepileupto(top, prime_ksearch(relations, newextra));
  }
  return gerepileupto(top, prime_ksearch(relations, extra));
}

//Given a rational matrix M, this returns the Hermite normal form of M with respect to columns.
GEN QM_hnf(GEN M){
  pari_sp top=avma;
  GEN l, N=Q_primitive_part(M, &l);//N is integral and primitive, N=M/l
  GEN hN=ZM_hnf(N);
  if(l==NULL) return gerepileupto(top, hN);//l=NULL happens when l=1 actually.
  return gerepileupto(top, ZM_Q_mul(hN, l));
}

//QM_hnf with typechecking
GEN QM_hnf_typecheck(GEN M){
  QM_check(M);
  return QM_hnf(M);
}

//Z_issquareall but working on Q not Z
int Q_issquareall(GEN x, GEN *sqrtx){
  pari_sp top=avma;
  if(typ(x)==t_INT) return Z_issquareall(x, sqrtx);//Integer, already done!
  GEN num=gel(x, 1), den=gel(x, 2), numrt, denrt;//numerator and denominator
  if(!Z_issquareall(num, &numrt)){avma=top;return 0;}
  if(!Z_issquareall(den, &denrt)){avma=top;return 0;}
  *sqrtx=gdiv(numrt, denrt);
  avma=top;
  return 1;
}

//ADD TYPECHECK

//Finds all subsets of L
GEN powerset(GEN L){
  pari_sp top=avma;
  long len=lg(L);
  if(len==1){return cgetg(1, t_VEC);}
  else if(len==2){
	GEN rvec=cgetg(3, t_VEC);
	gel(rvec, 1)=cgetg(1, t_VEC);
	gel(rvec, 2)=cgetg(2, t_VEC);
	gel(gel(rvec, 2), 1)=gcopy(gel(L, 1));
	return rvec;
  }
  GEN Lcop=cgetg(len-1, t_VEC);
  for(long i=1;i<len-1;i++) gel(Lcop, i)=gel(L, i);
  GEN halfs=powerset(Lcop);
  long lhalfs=lg(halfs), j=lhalfs-1, templ, newlen=lhalfs*2-1;
  len--;
  GEN full=cgetg(newlen, t_VEC);
  for(long i=1;i<lhalfs;i++){
    j++;
	gel(full, i)=gel(halfs, i);
	templ=lg(gel(halfs, i));
	gel(full, j)=cgetg(templ+1, t_VEC);
	for(long k=1;k<templ;k++) gel(gel(full, j), k)=gel(gel(halfs, i), k);
	gel(gel(full, j), templ)=gel(L, len);
  }
  return gerepileupto(top, gtoset(full));
}

//Assuming v1 and v2 are vectors in the same 1-dim linear subspace, this returns v1/v2 (which can be 0 if v1=0 and oo if v2=0).
GEN vecratio(GEN v1, GEN v2){
  for(long i=1;i<lg(v2);i++){
	if(gequal0(gel(v2, i))) continue;//Don't want to divide by the 0
	return gdiv(gel(v1, i), gel(v2, i));
  }
  return mkoo();
}

//I HAVENT INSTALLED THIS YET

//vecratio with typecheck
GEN vecratio_typecheck(GEN v1, GEN v2){
  if(typ(v1)!=t_VEC || typ(v2)!=t_VEC || lg(v1)!=lg(v2)) pari_err_TYPE("Please enter vectors of the same length that lie in the same linear subspace.", mkvec2(v1, v2));
  return vecratio(v1, v2);
}

