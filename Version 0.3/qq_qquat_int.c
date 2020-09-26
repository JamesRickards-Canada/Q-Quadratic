//Intersections in quaternion algebras over Q

#ifndef PARILIB
#define PARILIB
#include <pari/pari.h>
#endif

#ifndef METHDECL
#define METHDECL
#include "qquadraticdecl.h"
#endif


//STATIC DECLARATIONS
static GEN qa_inum_roots_f1bds(GEN r1, GEN r2, GEN n2, GEN coef);
static GEN qa_inum_roots_ghsearch(GEN Q, GEN ordinv, GEN rta, GEN maxd3d4, GEN D1, GEN rtD1, GEN e, GEN f, GEN res, GEN r1, GEN r2, GEN n2, int orient);



//INTERSECTION NUMBER BASED OFF OF ROOTS BOUND AREA




//Must have ei being the image if (p_Di+sqrt(Di))/2.
GEN qa_inum_roots(GEN Q, GEN order, GEN e1, GEN e2, GEN D1, GEN D2, int data, long prec){
  pari_sp top=avma, mid;
  GEN abvec=qa_getabvec(Q), pel1=pell(D1), pel2=pell(D2);
  int isswapped=0;
  if(cmpii(gel(pel1, 1), gel(pel2, 1))==-1){//Swapping so pell(D2)<pell(D1); this makes it MUCH faster when there is a large discrepancy
	isswapped=1;
	GEN temp=e1;e1=e2;e2=temp;
	temp=D1;D1=D2;D2=temp;
	temp=pel2;pel2=pel1;pel1=temp;
  }
  GEN e2pell=cgetg(5, t_VEC);
  gel(e2pell, 1)=gdivgs(gel(pel2, 1), 2);//T/2
  for(int i=2;i<=4;i++) gel(e2pell, i)=gmul(gel(pel2, 2), gel(e2, i));//U*e2[i]. Now e2pell is the image of epsilon_D2 under e2
  GEN e2roots=qa_roots(Q, e2pell, prec);
  GEN e2sl2r=qa_m2rembed(Q, e2pell);
  //We now need to pick an n1 in R union oo, and find all root geodesics of equivalent embeddings to e1 that intersect ell_{e2roots} (rt geod for e2pell) and ell_{n1, n2=e2sl2r(n1)}. If [r1, r2]=[0, oo], we choose [n1, n2]=[1, bigger]. Otherwise, assume r1<r2 and we find the n1 so that n1<r1<r2<n2 and (r1-n1)+(n2-r2) is maximized (this gives us the most efficient search later on).
  GEN r1, r2, n1, n2, mat, matinv;
  if(gcmp(gel(e2roots, 1), gel(e2roots, 2))==-1){r1=gel(e2roots, 1);r2=gel(e2roots, 2);mat=ginv(e2sl2r);matinv=e2sl2r;}
  else{r2=gel(e2roots, 1);r1=gel(e2roots, 2);mat=e2sl2r;matinv=ginv(e2sl2r);}//r1<r2 are the roots of e2pell, mat is either e2sl2r or its inverse. We take n1<r1, and then mat(n1)>r2.
  if(typ(r2)==t_INFINITY){//r2=oo, so r1=0 necessarily.
	n1=gen_1;
	n2=mat_eval(mat, n1);
  }
  else{//Finite
	n1=gsub(mat_eval(matinv, mkoo()), gdivsg(1, gcoeff(mat, 2, 1)));//gcoeff(mat, 2, 1)>0 necessarily
	n2=mat_eval(mat, n1);
  }
  //Great, now n1<r1<r2<n2 or [r1,r2]=[0,oo] and [n1,n2]=[1,something bigger]
  //We are now searching for W=[e,f,g,h] in Gamma, conjugate to e1, for which the roots of W intersect ell_(r1,r2) and ell_(n1,n2). First, e=e1[1] necessarily, and it turns out there are finitely many possibilities with just this restriction, as well as assuming [e,f,g,h] is in Z/d[1]+Z/d[2]i+Z/d[3]j+Z/d[4]k. We will find this set, check which are actually in Gamma, and check for conjugation with e1 afterwards.
  GEN maxdens=qa_getordmaxd(order), ord=qa_getord(order), ordinv=qa_getordinv(order);
  GEN maxd34=lcmii(gel(maxdens, 3), gel(maxdens, 4));
  GEN minf1, maxf1, rtD1=gsqrt(D1, prec), rta=gsqrt(gel(abvec, 1), prec), f, f1, poss, coef1, coef2, coef3, d3d4sq=Qdivii(sqri(maxd34), gel(abvec, 2));
  coef1=gmul(d3d4sq, gdivgs(D1, 4));
  coef2=gmul(d3d4sq, gneg(gel(abvec, 1)));//res=coef1+coef2*f^2=maxd3d4^2(D1/4-af^2)/b
  coef3=gmul(gel(maxdens, 2), gsqrt(gdiv(D1, shifti(gel(abvec, 1), 2)), prec));//coef3=maxdens[2]*sqrt(D1/4a); passed into f1bounds method.
  glist *S=NULL;
  GEN partials=cgetg(1, t_VEC), partialstemp;
  long nint=0;
  if(typ(r2)==t_INFINITY){
	minf1=gceil(gneg(coef3));
	maxf1=gneg(minf1);//Roots separate ell_{0, oo}
	mid=avma;
	for(GEN f1=minf1;cmpii(f1, maxf1)!=1;f1=addis(f1, 1)){
	  f=gdiv(f1, gel(maxdens, 2));
	  poss=qa_inum_roots_ghsearch(Q, ordinv, rta, maxd34, D1, rtD1, gel(e1, 1), f, gadd(coef1, gmul(coef2, gsqr(f))), n1, r2, n2, 1);
	  for(long i=1;i<lg(poss);i++){
		if(gequal0(qa_conjnorm(Q, ord, ordinv, gel(poss, i), e1, gen_1, 0, prec))) continue;//Nope, not conjugate
		glist_putstart(&S, gel(poss, i));
		nint++;
	  }
	  if(gc_needed(mid, 3)){
		long j=lg(partials);
		partialstemp=cgetg(j+nint, t_VEC);
		for(long i=1;i<j;i++) gel(partialstemp, i)=gcopy(gel(partials, i));
		partials=glist_togvec(S, nint, -1);
		for(long i=1;i<=nint;i++){gel(partialstemp, j)=gel(partials, i);j++;}
		partials=partialstemp;
		gerepileall(mid, 2, &f1, &partials);
		S=NULL;//Reset S
		nint=0;//Reset nint
	  }
	}
  }
  else{
	//Case 1: bigger root is >n2
	GEN fbds=qa_inum_roots_f1bds(r1, r2, n2, coef3);
	minf1=gel(fbds, 1);maxf1=gel(fbds, 2);
	mid=avma;
	for(GEN f1=minf1;cmpii(f1, maxf1)!=1;f1=addis(f1, 1)){
	  f=gdiv(f1, gel(maxdens, 2));
	  poss=qa_inum_roots_ghsearch(Q, ordinv, rta, maxd34, D1, rtD1, gel(e1, 1), f, gadd(coef1, gmul(coef2, gsqr(f))), r1, r2, n2, 1);
	  for(long i=1;i<lg(poss);i++){
		if(gequal0(qa_conjnorm(Q, ord, ordinv, gel(poss, i), e1, gen_1, 0, prec))) continue;//Nope, not conjugate
		glist_putstart(&S, gel(poss, i));
		nint++;
	  }
	  if(gc_needed(mid, 3)){
		long j=lg(partials);
		partialstemp=cgetg(j+nint, t_VEC);
		for(long i=1;i<j;i++) gel(partialstemp, i)=gcopy(gel(partials, i));
		partials=glist_togvec(S, nint, -1);
		for(long i=1;i<=nint;i++){gel(partialstemp, j)=gel(partials, i);j++;}
		partials=partialstemp;
		gerepileall(mid, 2, &f1, &partials);
		S=NULL;//Reset S
		nint=0;//Reset nint
	  }
	}
	//Case 2: smaller root is <n1. We copy the previous code, replacing f with -f (also exchanging minf1 and maxf1), negating n1 and r1, r2
	fbds=qa_inum_roots_f1bds(gneg(r2), gneg(r1), gneg(n1), coef3);
	maxf1=gneg(gel(fbds, 1));minf1=gneg(gel(fbds, 2));
	mid=avma;
	for(GEN f1=minf1;cmpii(f1, maxf1)!=1;f1=addis(f1, 1)){
	  f=gdiv(f1, gel(maxdens, 2));
	  poss=qa_inum_roots_ghsearch(Q, ordinv, rta, maxd34, D1, rtD1, gel(e1, 1), gneg(f), gadd(coef1, gmul(coef2, gsqr(f))), gneg(r2), gneg(r1), gneg(n1), -1);
	  for(long i=1;i<lg(poss);i++){
		if(gequal0(qa_conjnorm(Q, ord, ordinv, gel(poss, i), e1, gen_1, 0, prec))) continue;//Nope, not conjugate
		glist_putstart(&S, gel(poss, i));
		nint++;
	  }
	  if(gc_needed(mid, 3)){
		long j=lg(partials);
		partialstemp=cgetg(j+nint, t_VEC);
		for(long i=1;i<j;i++) gel(partialstemp, i)=gcopy(gel(partials, i));
		partials=glist_togvec(S, nint, -1);
		for(long i=1;i<=nint;i++){gel(partialstemp, j)=gel(partials, i);j++;}
		partials=partialstemp;
		gerepileall(mid, 2, &f1, &partials);
		S=NULL;//Reset S
		nint=0;//Reset nint
	  }
	}
	//Case 3: One root is infinite, the other 0. Then e^2-af^2=nm, and we use the f1 variable to mean f. 
	if(gsigne(r1)==-1 && gsigne(r2)==1){
	  f1=Qdivii(D1, shifti(gel(abvec, 1), 2));//D1/4a
	  GEN f1root;
	  if(Q_issquareall(f1, &f1root)){
		GEN emb=mkvec4(gel(e1, 1), f1root, gen_0, gen_0);
		if(!gequal0(qa_conjnorm(Q, ord, ordinv, emb, e1, gen_1, 0, prec))){glist_putstart(&S, emb);nint++;}
		emb=qa_conj(emb);
		if(!gequal0(qa_conjnorm(Q, ord, ordinv, emb, e1, gen_1, 0, prec))){glist_putstart(&S, emb);nint++;}
	  }
	}
  }
  GEN e1set=glist_togvec(S, nint, -1);
  long j=lg(partials);
  GEN pairs=cgetg(nint+j, t_VEC);
  if(!isswapped){
	for(long i=1;i<j;i++){
	  gel(pairs, i)=cgetg(3, t_VEC);
	  gel(gel(pairs, i), 1)=gcopy(gel(partials, i));
	  gel(gel(pairs, i), 2)=gcopy(e2);
	}
    for(long i=1;i<=nint;i++){
	  gel(pairs, j)=cgetg(3, t_VEC);
	  gel(gel(pairs, j), 1)=gcopy(gel(e1set, i));
	  gel(gel(pairs, j), 2)=gcopy(e2);
	  j++;
	}
  }
  else{
	for(long i=1;i<j;i++){
	  gel(pairs, i)=cgetg(3, t_VEC);
	  gel(gel(pairs, i), 1)=gcopy(e2);
	  gel(gel(pairs, i), 2)=gcopy(gel(partials, i));
	}
	for(long i=1;i<=nint;i++){
	  gel(pairs, j)=cgetg(3, t_VEC);
	  gel(gel(pairs, j), 1)=gcopy(e2);
	  gel(gel(pairs, j), 2)=gcopy(gel(e1set, i));
	  j++;
	}
  }
  return gerepileupto(top, pairs);
}

//qa_inum_roots with typechecking
GEN qa_inum_roots_typecheck(GEN Q, GEN order, GEN e1, GEN e2, int data, long prec){
  pari_sp top=avma;
  qa_indefcheck(Q);
  qa_ordeichlercheck(order);
  qa_eltcheck(e1);qa_eltcheck(e2);
  GEN v1=qa_associatedemb(Q, order, e1, gen_0), v2=qa_associatedemb(Q, order, e2, gen_0);
  return gerepileupto(top, qa_inum_roots(Q, order, gel(v1, 1), gel(v2, 1), gel(v1, 2), gel(v2, 2), 0, prec));
}

//coef=maxdens[2]*sqrt(D1/4a)
//This generates the bounds for f1 based off of n1<r1< R1 <r2<n2 < R2, where the roots we search for are R1 and R2. Note that n1 is not needed in the input.
static GEN qa_inum_roots_f1bds(GEN r1, GEN r2, GEN n2, GEN coef){
  pari_sp top=avma;
  GEN minf1, maxf1;
  if(gsigne(n2)==1){
	if(gsigne(r2)==1){
	  maxf1=gfloor(gmul(coef, gdiv(gadd(r2, n2), gsub(n2, r2))));//Combining R2>n2 and R1<r2
	  if(gsigne(r1)==1) minf1=gceil(coef);//Use R1>r1>0
	  else minf1=gceil(gneg(coef));//Use R2>n2>0
	}
	else{
	  maxf1=gfloor(coef);//R1<r2<0
	  minf1=gceil(gneg(coef));//R2>n2>0
	}
  }
  else{
    maxf1=gfloor(coef);//R1<r2<0
	minf1=gceil(gmul(coef, gdiv(gadd(n2, r2), gsub(n2, r2))));//R1<r2 and R2>n2
  }
  GEN ret=cgetg(3, t_VEC);
  gel(ret, 1)=gcopy(minf1);
  gel(ret, 2)=gcopy(maxf1);
  return gerepileupto(top, ret);
}


//res=maxd3d4^2(D1/4-af^2)/b
//Starting with e, f, D1, r1, r2, n1, n2, we derive bounds on g, h and find the finite set of possible pairs(g, h) making embeddings [e, f, g, h] of discriminant D1 that intersect both root geodesics ell_{r1, r2} and ell_{n1, n2}. maxd3d4 is the lcm of qa_getordmaxd[3], [4]. Note that n1 is not inputted, it is not needed. Through the various cases we use this, sometimes we use r1 as n1, permute the inputs, etc. orient=1 if one root is in (r1, r2) and the other is >n2, and is 0 if the other is <n1. res=maxd3d4^2(D1/4-af^2)/b. Note that we are only finding finite rooted solutions here, must search for g=h=0 elsewhere.
static GEN qa_inum_roots_ghsearch(GEN Q, GEN ordinv, GEN rta, GEN maxd3d4, GEN D1, GEN rtD1, GEN e, GEN f, GEN res, GEN r1, GEN r2, GEN n2, int orient){
  pari_sp top=avma;
  if(gequal0(res) || typ(res)!=t_INT) return cgetg(1, t_VEC);
  GEN forient, a=qa_geta(Q);
  if(orient==1) forient=f;//At the end, if orient=-1 we are checking with -f not f
  else forient=gneg(f);
  GEN gmhmax, gmhmin;//We use the root bounds with the standard embedding to get bounds on g-hsqrt(a)
  if(typ(r2)==t_INFINITY){
	GEN ghnumer=gadd(gmul(f, rta), gdivgs(rtD1, 2));
	gmhmax=gdiv(ghnumer, r1);
	gmhmin=gdiv(ghnumer, n2);
  }
  else{
	GEN ghnumer=gsub(gmul(f, rta), gdivgs(rtD1, 2));
	gmhmin=gen_0, gmhmax=mkoo();
	if(gsigne(n2)==1) gmhmax=gdiv(gadd(ghnumer, rtD1), n2);
	else gmhmin=gmax_shallow(gdiv(gadd(ghnumer, rtD1), n2), gen_0);
	if(gsigne(r1)==1) gmhmax=gmin_shallow(gdiv(ghnumer, r1), gmhmax);
	else gmhmin=gmax_shallow(gdiv(ghnumer, r1), gmhmin);
	if(gsigne(r2)==-1) gmhmax=gmin_shallow(gdiv(ghnumer, r2), gmhmax);
	else gmhmin=gmax_shallow(gdiv(ghnumer, r2), gmhmin);
  }
  gmhmin=gmul(maxd3d4, gmhmin);//Scaling so that we get bounds on integers (maxd3d4*[g, h] is an integer vector necessarily.
  gmhmax=gmul(maxd3d4, gmhmax);
  if(gcmp(gmhmin, gmhmax)==1){avma=top;return cgetg(1, t_VEC);}//min>max
  //The geodesic condition is <==> gmhmin<g1-h1*sqrt(a)gphmax (for g-hsqrt(a)>0).
  //From (g-hsqrt(a))(g+hsqrt(a))=res, we get bounds for g+hsqrt(a)
  GEN gphmin, gphmax;//would be more accurate to name g1ph1min/max since the bounds we get everywhere are for g1 and h1 not h and h
  if(gsigne(res)==1){//res>0
    gphmin=gdiv(res, gmhmax);
	gphmax=gdiv(res, gmhmin);
  }
  else{//res<0
    gphmin=gdiv(res, gmhmin);
	gphmax=gdiv(res, gmhmax);
  }
  GEN g1min=gceil(gdivgs(gadd(gmhmin, gphmin), 2));//Bounds for g1
  GEN g1max=gfloor(gdivgs(gadd(gmhmax, gphmax), 2));
  GEN h1sq, h1abs, g, h, gmh, emb;
  glist *S=NULL;
  long nfound=0;
  for(GEN g1=g1min;cmpii(g1, g1max)!=1;g1=addis(g1, 1)){
	h1sq=Qdivii(subii(sqri(g1), res), a);
	if(typ(h1sq)!=t_INT) continue;//Not integer
	if(!Z_issquareall(h1sq, &h1abs)) continue;//Not square
	g=Qdivii(g1, maxd3d4);
	h=Qdivii(h1abs, maxd3d4);
	gmh=gsub(g1, gmul(h1abs, rta));//remember, bounds are for g1 and h1 not g and h
	if(gcmp(gmhmin, gmh)==-1 && gcmp(gmh, gmhmax)==-1){//Yes!
	  emb=mkvec4(e, forient, g, h);
	  if(qa_isinorder(Q, ordinv, emb)){glist_putstart(&S, emb);glist_putstart(&S, qa_conj(emb));nfound=nfound+2;}//Must add conjugate due to assuming g-hsqrt(a)>0
	}
	if(gequal0(h)) continue;//Now we do (g, -h), but of course only if h!=0
	h=gneg(h);
	gmh=gadd(g1, gmul(h1abs, rta));
	if(gcmp(gmhmin, gmh)==-1 && gcmp(gmh, gmhmax)==-1){//Yes!
	  emb=mkvec4(e, forient, g, h);
	  if(qa_isinorder(Q, ordinv, emb)){glist_putstart(&S, emb);glist_putstart(&S, qa_conj(emb));nfound=nfound+2;}//Must add conjugate due to assuming g-hsqrt(a)>0
	}
  }
  return gerepileupto(top, glist_togvec(S, nfound, -1));
}



//INTERSECTIONS BASED ON X-LINKAGE



//Intersection number based on x-linking.
GEN qa_inum_x(GEN Q, GEN order, GEN e1, GEN e2, GEN D1, GEN D2, int data, long prec){
  pari_sp top=avma;
  GEN xposs=qa_xposs(qa_getpram(Q), qa_getpramprod(Q), D1, D2, gen_0, gen_0);//The possible x's
  long lgxposs=lg(xposs), nx=2*lgxposs-1;//Possible x's; double since we account for -x
  if(lgxposs==1){avma=top;return cgetg(1, t_VEC);}//No solutions
  GEN pairs=cgetg(nx, t_VEC);
  long i=1, j=1, npairs=1;
  if(gequal0(gel(xposs, 1))){
	gel(pairs, 1)=qa_xlink(Q, order, e1, e2, D1, D2, gen_0, prec);
	i=2;j=2;
	npairs=npairs+lg(gel(pairs, 1))-1;
  }
  for(;i<lgxposs;i++){
    gel(pairs, j)=qa_xlink(Q, order, e1, e2, D1, D2, gel(xposs, i), prec);//x
	npairs=npairs+lg(gel(pairs, j))-1;j++;
	gel(pairs, j)=qa_xlink(Q, order, e1, e2, D1, D2, negi(gel(xposs, i)), prec);//-x
	npairs=npairs+lg(gel(pairs, j))-1;j++;
  }
  GEN allpairs=cgetg(npairs, t_VEC);
  long pos=1;
  for(i=1;i<j;i++){
	for(long k=1;k<lg(gel(pairs, i));k++){gel(allpairs, pos)=gcopy(gel(gel(pairs, i), k));pos++;}
  }
  return gerepileupto(top, allpairs);
}

//qa_inum_x with typecheck
GEN qa_inum_x_typecheck(GEN Q, GEN order, GEN e1, GEN e2, int data, long prec){
  pari_sp top=avma;
  qa_indefcheck(Q);
  qa_ordeichlercheck(order);
  qa_eltcheck(e1);qa_eltcheck(e2);
  GEN v1=qa_associatedemb(Q, order, e1, gen_0), v2=qa_associatedemb(Q, order, e2, gen_0);
  return gerepileupto(top, qa_inum_x(Q, order, gel(v1, 1), gel(v2, 1), gel(v1, 2), gel(v2, 2), 0, prec));
}


// The output is [[pairs],[signed levels]], with pairs corresponding to respective levels."); TO DO?

//Finds all x-linkage of e1 with e2.
GEN qa_xlink(GEN Q, GEN order, GEN e1, GEN e2, GEN D1, GEN D2, GEN x, long prec){
  pari_sp top=avma;
  GEN zb=qa_getordtrace0basis(order), ord=qa_getord(order), ordinv=qa_getordinv(order);
  GEN v1=gmulgs(e1, 2);gel(v1, 1)=gen_0;//v1 represents sqrt(D1) under phi1
  GEN quadf=cgetg(7, t_VEC);//The embedding of D2 is given by emb=A*zb[1]+B*zb[2]+C*zb[3]. solving emb^2=D2 gives us the quadratic form quadf.
  gel(quadf, 1)=negi(qa_norm(Q, gel(zb, 1)));//A^2 coeff
  gel(quadf, 2)=negi(qa_norm(Q, gel(zb, 2)));//B^2 coeff
  gel(quadf, 3)=negi(qa_norm(Q, gel(zb, 3)));//C^2 coeff
  gel(quadf, 4)=qa_trace(qa_mul(Q, gel(zb, 1), gel(zb, 2)));//AB coeff
  gel(quadf, 5)=qa_trace(qa_mul(Q, gel(zb, 1), gel(zb, 3)));//AC coeff
  gel(quadf, 6)=qa_trace(qa_mul(Q, gel(zb, 2), gel(zb, 3)));//BC coeff
  GEN lin=cgetg(4, t_VEC);//The x-linking gives the equation (v1*emb)[1]=x, and solving this gives the linear form lin.
  gel(lin, 1)=gel(qa_mul(Q, v1, gel(zb, 1)), 1);//A coeff
  gel(lin, 2)=gel(qa_mul(Q, v1, gel(zb, 2)), 1);//B coeff
  gel(lin, 3)=gel(qa_mul(Q, v1, gel(zb, 3)), 1);//C coeff
  GEN ABCs=bqf_linearsolve(quadf, D2, lin, x, prec);//Solving for A, B, C
  if(gequal0(ABCs)){avma=top;return cgetg(1, t_VEC);}//No solutions
  if(!equali1(gel(gel(ABCs, 1), 1))) pari_err_TYPE("Solutions are of an unexpected shape, please report!", ABCs);//The return should be type 1, unless D_1D_2-x^2=0, and we do NOT consider that case.
  //Now ABCs=[[1, M, [s1, s2, s3]], [a1, b1, c1], [a2, b2, c2], ...], where a general solution is [A;B;C]=M^k*[ai;bi;ci]+[s1;s2;s3]
  GEN scol=gtocol(gel(gel(ABCs, 1), 3)), M=gel(gel(ABCs, 1), 2);//s in column format, M
  GEN sshift=gsub(scol, gmul(M, scol));//scol-M*scol; to go from k->k+1, multiply the solution by M and add sshift
  GEN base, embofD2, baseemb=zerovec(4), emb=zerovec(4);
  if(smodis(D2, 2)==1){gel(baseemb, 1)=ghalf;gel(emb, 1)=ghalf;}//Emb represents the embedding, and baseemb the base embedding
  glist *S=NULL;//List of the embeddings
  long count=0;
  for(long i=2;i<lg(ABCs);i++){
	base=gadd(gtocol(gel(ABCs, i)), scol);//The first solution
	embofD2=gadd(gadd(gmul(gel(zb, 1), gel(base, 1)), gmul(gel(zb, 2), gel(base, 2))), gmul(gel(zb, 3), gel(base, 3)));//The embedding of sqrt(D2)
	for(int j=2;j<=4;j++) gel(baseemb, j)=gdivgs(gel(embofD2, j), 2);//The embedding of (p_{D2}+sqrt(D2))/2
	if(equali1(qa_conjnorm(Q, ord, ordinv, e2, baseemb, gen_1, 0, prec))){glist_putstart(&S, gcopy(baseemb));count++;}//Good!
    for(;;){
	  base=gadd(gmul(M, base), sshift);
	  embofD2=gadd(gadd(gmul(gel(zb, 1), gel(base, 1)), gmul(gel(zb, 2), gel(base, 2))), gmul(gel(zb, 3), gel(base, 3)));//The embedding of sqrt(D2)
	  for(int j=2;j<=4;j++) gel(emb, j)=gdivgs(gel(embofD2, j), 2);//The embedding of (p_{D2}+sqrt(D2))/2
	  if(equali1(qa_norm(Q, qa_simulconj(Q, ord, ordinv, e1, baseemb, e1, emb, prec)))) break;//Back to equivalent embeddings.
	  if(equali1(qa_conjnorm(Q, ord, ordinv, e2, emb, gen_1, 0, prec))){glist_putstart(&S, gcopy(emb));count++;}//Good!
	}
  }
  GEN allposs=gerepileupto(top, glist_togvec(S, count, -1));//All the possibilities for emb2. We do want to make sure we did not double count, so we go through and check them.
  if(gequal0(allposs)){avma=top;return cgetg(1, t_VEC);}//No solutions
  pari_sp mid;
  long pos=1, lgallpos=lg(allposs);
  GEN final=vectrunc_init(lgallpos), tempelt, n;
  vectrunc_append(final, cgetg(3, t_VEC));
  gel(gel(final, 1), 1)=gcopy(e1);
  gel(gel(final, 1), 2)=gcopy(gel(allposs, 1));
  for(long i=2;i<lgallpos;i++){
	mid=avma;
    tempelt=qa_simulconj(Q, ord, ordinv, e1, gel(gel(final, 1), 2), e1, gel(allposs, i), prec);
	n=qa_norm(Q, tempelt);
	if(equali1(n)){avma=mid;continue;}//Not a new solution.
	avma=mid;
	vectrunc_append(final, cgetg(3, t_VEC));
	pos++;
	gel(gel(final, pos), 1)=gcopy(e1);
    gel(gel(final, pos), 2)=gcopy(gel(allposs, i));
  }
  return gerepileupto(top, final);
}

//qa_xlink with typechecking
GEN qa_xlink_typecheck(GEN Q, GEN order, GEN e1, GEN e2, GEN x, long prec){
  pari_sp top=avma;
  qa_indefcheck(Q);
  qa_ordeichlercheck(order);
  qa_eltcheck(e1);qa_eltcheck(e2);
  GEN v1=qa_associatedemb(Q, order, e1, gen_0), v2=qa_associatedemb(Q, order, e2, gen_0);
  return gerepileupto(top, qa_xlink(Q, order, gel(v1, 1), gel(v2, 1), gel(v1, 2), gel(v2, 2), x, prec));
}

//Finds the set of x's bewteen xmin and xmax (inclusive) which exhibit x-linking in the indefinite quaternion algebra ramified at pset. 
GEN qa_xposs(GEN pset, GEN Psetprod, GEN D1, GEN D2, GEN xmin, GEN xmax){
  pari_sp top=avma;
  GEN D=mulii(D1, D2), r;
  GEN rtD=sqrtremi(D, &r);//D=rtD^2+r
  if(gequal0(xmax) && gequal0(xmin)){//Auto setting max
    if(gequal0(r)) xmax=subis(rtD, 1);//Must remove rtD from the range.
	else xmax=rtD;
  }
  if(smodis(D, 2)!=smodis(xmin, 2)) xmin=addis(xmin, 1);//making xmin==D (2)
  int isbad=0;
  if(gequal0(r) && cmpii(xmin, rtD)!=1 && cmpii(rtD, xmax)!=1) isbad=1;//D is a square, must avoid rootD
  GEN N=shifti(subii(sqri(xmin), D), -2), set2, x;//N=(x^2-D1D2)/4
  glist *S=NULL;//Stores the set of x's
  long count=0;
  if(!isbad){
    for(x=xmin;cmpii(x, xmax)!=1;x=addis(x, 2)){//For loop over all x's
	  if(!gequal0(modii(N, Psetprod))){
	    N=addii(N, addis(x, 1));//N->N+x+1, i.e. from (x^2-D1D2)/4 to ((x+2)^2-D1D2)/4
	    continue;//Need PsetProd to divide N
	  }
	  set2=qa_ram_fromab(D1, N);
	  if(gequal(set2, pset)){glist_putstart(&S, x);count++;}
	  N=addii(N, addis(x, 1));//N->N+x+1, i.e. from (x^2-D1D2)/4 to ((x+2)^2-D1D2)/4
    }
  }
  else{
	for(x=xmin;cmpii(x, rtD)==-1;x=addis(x, 2)){//For loop over all x's up to rtD
	  if(!gequal0(modii(N, Psetprod))){
	    N=addii(N, addis(x, 1));//N->N+x+1, i.e. from (x^2-D1D2)/4 to ((x+2)^2-D1D2)/4
	    continue;//Need PsetProd to divide N
	  }
	  set2=qa_ram_fromab(D1, N);
	  if(gequal(set2, pset)){glist_putstart(&S, x);count++;}
	  N=addii(N, addis(x, 1));//N->N+x+1, i.e. from (x^2-D1D2)/4 to ((x+2)^2-D1D2)/4
    }
	N=addii(N, addis(x, 1));//At this point, x=rtD, so we increment by 2
	for(GEN x=addis(rtD, 2);cmpii(x, xmax)!=1;x=addis(x, 2)){//For loop over all x's from rtD to xmax
	  if(!gequal0(modii(N, Psetprod))){
	    N=addii(N, addis(x, 1));//N->N+x+1, i.e. from (x^2-D1D2)/4 to ((x+2)^2-D1D2)/4
	    continue;//Need PsetProd to divide N
	  }
	  set2=qa_ram_fromab(D1, N);
	  if(gequal(set2, pset)){glist_putstart(&S, x);count++;}
	  N=addii(N, addis(x, 1));//N->N+x+1, i.e. from (x^2-D1D2)/4 to ((x+2)^2-D1D2)/4
    }
  }
  return gerepileupto(top, glist_togvec(S, count, -1));
}

//qa_xposs with typecheck and allowing either Q or pset to be inputted
GEN qa_xposs_typecheck(GEN Qorpset, GEN D1, GEN D2, GEN xmin, GEN xmax){
  pari_sp top=avma;
  if(typ(Qorpset)!=t_VEC) pari_err_TYPE("Qorpset is either a quaternion algebra or even sized set of primes", Qorpset);
  if(typ(xmin)!=t_INT || typ(xmax)!=t_INT) pari_err_TYPE("xmin and xmax must be integers", xmin);
  if(!isdisc(D1) || !isdisc(D2)) pari_err_TYPE("D1 and D2 must be discriminants", D1);
  if(lg(Qorpset)==1) return qa_xposs(Qorpset, gen_1, D1, D2, xmin, xmax);//Ramified nowhere
  if(lg(Qorpset)==2 || typ(gel(Qorpset, 2))==t_INT){//Qorpset is the prime set
	GEN prod=gen_1;
	for(long i=1;i<lg(Qorpset);i++) prod=mulii(prod, gel(Qorpset, i));
	return gerepileupto(top, qa_xposs(Qorpset, prod, D1, D2, xmin, xmax));
  }
  qa_check(Qorpset);
  return qa_xposs(qa_getpram(Qorpset), qa_getpramprod(Qorpset), D1, D2, xmin, xmax);
}



