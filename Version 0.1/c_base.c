//This is a collection of miscellaneous methods that may be useful in a variety of settings, and not just for the programs they were originally created for.

#ifndef PARILIB
#define PARILIB
#include <pari/pari.h>
#endif

#ifndef METHDECL
#define METHDECL
#include "qquadraticdecl.h"
#endif

#ifndef TIME
#define TIME
#include <time.h>
#endif

//COMPLEX GEOMETRY

/*
GP;install("crossratio","GGGG","crossratio","./c_base.so");
GP;addhelp(crossratio, "Inputs a,b,c,d complex numbers or oo with at most one being oo.\n Outputs the crossratio of [a,b;c,d].");
GP;install("mat_eval_typecheck","GG","mat_eval","./c_base.so");
GP;addhelp(mat_eval, "Inputs M,x; M a matrix, and x number.\n Outputs Mx with M acting via Mobius transformation. x=+/-oo is allowed.");
*/

//INFINITY

/*
GP;install("addoo","GG","addoo","./c_base.so");
GP;addhelp(addoo, "Inputs a,b real numbers or oo.\n Outputs a+b, where if a or b is +/- oo, returns that back. Note that this will make oo+-oo=oo and -oo+oo=-oo");
GP;install("divoo","GG","divoo","./c_base.so");
GP;addhelp(divoo, "Inputs a,b, real numbers or oo.\n Outputs a/b, where oo is output if a=oo and b>=0 or a=-oo and b<0 or b=0 and a>=0 (outputs -oo under analogous assumptions).");
*/

//LINEAR EQUATIONS AND MATRICES

/*
GP;install("lin_intsolve_typecheck","GGG","lin_intsolve","./c_base.so");
GP;addhelp(lin_intsolve, "Inputs A,B,n integers.\n Outputs the general integral solutions to Ax+By=n. The format is [[s1,s2],[x0,y0]], where the general solution is x=s1*t+x0, y=s2*t+y0 for t an integer. The output is also reduced, i.e. gcd(s1,s2)=1. If A=B=0 or there are no integer solutions, returns 0.");
GP;install("mat3_complete_typecheck", "GGG", "mat3_complete", "./c_base.so");
GP;addhelp(mat3_complete, "Inputs A,B,C integers with gcd 1.\n Outputs a 3x3 integer matrix with determinant 1 whose top row is [A, B, C].");
*/

//SOLVING EQUATIONS MOD N

/*
GP;install("sqmod_typecheck","GG","sqmod","./c_base.so");
GP;addhelp(sqmod, "Inputs: x,n: rational number x, integer n with gcd(n,denom(x))=1.\n Solves y^2==x mod n, and outputs the solutions. The output format is 0 if no solutions, and [S,m] otherwise, where the solutions are y==S[i] modulo m.");
*/

//TIME

/*
GP;install("printtime","v","printtime","./c_base.so");
GP;addhelp(printtime, "Prints current time");
*/

//GENERAL HELP

/*
GP;addhelp(base,"This package is a collection of miscellaneous methods that may be useful in a variety of settings, and not just for the programs they were originally created for \n Subtopics: \n Complex plane geometry (cpg) \n Infinity (inf) \n Integer linear equations and matrices (ilm) \n Modulo n methods (modn) \n Time (time)");
GP;addhelp(cpg,"crossratio, mat_eval.");
GP;addhelp(inf,"addoo, divoo.");
GP;addhelp(ilm,"lin_intsolve, mat3_complete.");
GP;addhelp(modn,"sqmod.");
GP;addhelp(time,"printtime.");
*/

//STATIC METHOD DECLARATIONS
static GEN sqmod_ppower(GEN x, GEN p, long n, GEN p2n, int iscoprime);


//COMPLEX GEOMETRY
  
  
//Crossratio, allows one entry to be infinite.
GEN crossratio(GEN a, GEN b, GEN c, GEN d){
  pari_sp top=avma;
  if (typ(a)==t_INFINITY) return gerepileupto(top,divoo(gsub(d, b), gsub(c, b)));
  if (typ(b)==t_INFINITY) return gerepileupto(top,divoo(gsub(c, a), gsub(d, a)));
  if (typ(c)==t_INFINITY) return gerepileupto(top,divoo(gsub(d, b), gsub(d, a)));
  if (typ(d)==t_INFINITY) return gerepileupto(top,divoo(gsub(c, a), gsub(c, b)));
  return gerepileupto(top,divoo(gmul(gsub(c, a), gsub(d, b)), gmul(gsub(c, b), gsub(d, a))));
}

//Mx, where M is a 2x2 matrix and x is complex or infinite.
GEN mat_eval(GEN M, GEN x){
  pari_sp top = avma;
  if(typ(x)==t_INFINITY) return gerepileupto(top,divoo(gcoeff(M, 1, 1), gcoeff(M, 2, 1)));
  return gerepileupto(top,divoo(gadd(gmul(gcoeff(M, 1, 1), x), gcoeff(M, 1, 2)), gadd(gmul(gcoeff(M, 2, 1), x), gcoeff(M, 2, 2))));
}

//This checks that the inputs are valid.
GEN mat_eval_typecheck(GEN M, GEN x){
  if(typ(M)!=t_MAT) pari_err(e_TYPE,"mat_eval. M should be a 2x2 matrix",M);
  if(lg(M)!=3 || lg(gel(M,1))!=3) pari_err(e_DIM,"mat_eval. M should be a 2x2 matrix");
  return mat_eval(M,x);
}


//INFINITY 


//Adds a,b, and allows for oo
GEN addoo(GEN a, GEN b){//No need to do garbage collection
  if(typ(a)==t_INFINITY){
    if(inf_get_sign(a)==1) return mkoo();
	else return mkmoo();
  }
  if(typ(b)==t_INFINITY){
    if(inf_get_sign(b)==1) return mkoo();
	else return mkmoo();
  }
  return gadd(a,b);
}

//Divides a and b, and allows for oo and division by 0. Returns oo for 0/0.
GEN divoo(GEN a, GEN b){//No garbage collection necessary
  if(gequal0(b)){//b=0
    if(gcmpgs(a,0)>=0) return mkoo();
    return mkmoo();
  }
  if(typ(a)==t_INFINITY){//a=+/-oo
	if(gsigne(a)==gsigne(b)) return mkoo();
	return mkmoo();
  }
  if(typ(b)==t_INFINITY) return gen_0;
  return gdiv(a,b);
}


//LINEAR EQUATIONS AND MATRICES


//Solves Ax+By=n
GEN lin_intsolve(GEN A, GEN B, GEN n){
  pari_sp top = avma;
  if(gequal0(A) && gequal0(B)) return gen_0;
  GEN g,X,Y;
  g=gbezout(A,B,&X,&Y);//XA+YB=g now
  GEN scale=Qdivii(n,g);
  if(typ(scale)!=t_INT){//g does not divide n
	avma=top;
	return gen_0;
  }
  GEN rvec=cgetg(3,t_VEC);//Creating the return vector
  gel(rvec,1)=cgetg(3,t_VEC);
  gel(rvec,2)=cgetg(3,t_VEC);
  gel(gel(rvec,2),1)=mulii(X,scale);
  gel(gel(rvec,2),2)=mulii(Y,scale);//Second part initialized
  if(gequal0(A)){
    gel(gel(rvec,1),1)=gen_1;
	gel(gel(rvec,1),2)=gen_0;
  }
  else if(gequal0(B)){
	gel(gel(rvec,1),1)=gen_0;
	gel(gel(rvec,1),2)=gen_1;
  }
  else{
	gel(gel(rvec,1),1)=diviiexact(negi(B),g);
	gel(gel(rvec,1),2)=diviiexact(A,g);
  }
  return(gerepilecopy(top, rvec));
}

//lin_intsolve with typecheck
GEN lin_intsolve_typecheck(GEN A, GEN B, GEN n){
  if (typ(A)!=t_INT) pari_err_TYPE("lin_intsolve. A should be an integer",A);
  if (typ(B)!=t_INT) pari_err_TYPE("lin_intsolve. B should be an integer",B);
  if (typ(n)!=t_INT) pari_err_TYPE("lin_intsolve. n should be an integer",n);
  return lin_intsolve(A,B,n);
}

//Returns a 3x3 ZM with top row A, B, C. Assumes gcd(A,B,C)=1
GEN mat3_complete(GEN A, GEN B, GEN C){
  pari_sp top=avma;
  GEN u, v, w, x;
  GEN g=bezout(A, B, &u, &v);
  bezout(g,C, &x, &w);
  u=mulii(u,x);
  v=mulii(v,x);
  togglesign_safe(&v);//Now uA-vB+wC=g2=1 since gcd(A,B,C)=1. Write the bottom two rows as [d,f,g;g,h,i]
  GEN m;//We want to now solve ei-fh=u; di-fg=v; dh-eg=w.
  if(gequal0(u) && gequal0(v)){
	m=zeromatcopy(3,3);
	gcoeff(m,2,1)=gen_1;gcoeff(m,2,2)=gen_0;gcoeff(m,2,3)=gen_0;
    gcoeff(m,3,1)=gen_0;gcoeff(m,3,3)=gen_0;
	if(equali1(w)) gcoeff(m,3,2)=gen_1;
	else gcoeff(m,3,2)=gen_m1;//Making det 1
  }
  else{
	GEN p, q;
    g=bezout(u, v, &q, &p);
    GEN H=diviiexact(u, g);togglesign_safe(&H);//We use captials to represent the matrix elements now
    GEN G=diviiexact(v, g);togglesign_safe(&G);
    GEN E=mulii(p, w);
    GEN D=mulii(q, w);togglesign_safe(&D);//DH-EG=w, F=g and I=0
    m=zeromatcopy(3,3);
	gcoeff(m,2,1)=icopy(D);gcoeff(m,2,2)=icopy(E);gcoeff(m,2,3)=icopy(g);
    gcoeff(m,3,1)=icopy(G);gcoeff(m,3,2)=icopy(H);gcoeff(m,3,3)=gen_0;
  }
  gcoeff(m,1,1)=icopy(A);gcoeff(m,1,2)=icopy(B);gcoeff(m,1,3)=icopy(C);
  return gerepileupto(top, m);
}

//Returns a 3x3 ZM with top row A, B, C. Assumes gcd(A,B,C)=1
GEN mat3_complete_typecheck(GEN A, GEN B, GEN C){
  pari_sp top=avma;
  if(typ(A)!=t_INT) pari_err_TYPE("Please enter three integers with gcd 1", A);
  if(typ(B)!=t_INT) pari_err_TYPE("Please enter three integers with gcd 1", B);
  if(typ(C)!=t_INT) pari_err_TYPE("Please enter three integers with gcd 1", C);
  GEN g=gcdii(A, B);
  GEN g1=gcdii(g,C);
  if(!equali1(g1)) pari_err_TYPE("GCD is not equal to 1", mkvec3(A, B, C));
  avma=top;
  return mat3_complete(A, B, C);
}


//SOLVING EQUATIONS MOD N


//i should update this so that we use something like forvec or w/e
//Solves y^2==x mod n. n must be positive integer, x is rational with gcd(denominator(x),n)=1, fact is the factorization of x (and is optional)
GEN sqmod(GEN x, GEN n, GEN fact){
  pari_sp top = avma;
  if(equali1(n)){//n=1
	GEN rvec=cgetg(3,t_VEC);//The return vector
	gel(rvec,1)=mkvec(gen_0);//First entry is [0]
	gel(rvec,2)=gen_1;//Second entry is [1]
    return(gerepileupto(top,rvec));
  }
  if(gequal0(fact)) fact = Z_factor(n);//Factorization matrix, which can be passed in. It is non-trivial as n>1 necessarily now.
  if(typ(x)==t_INT) x=Fp_red(x, n);//Reducing x.
  else x=Fp_div(gel(x,1), gel(x,2), n);//x is t_FRAC, so gel(x,1)=numerator and gel(x,2)=denominator
  long nprimes=lg(gel(fact,1))-1, nsols=1, lx;
  GEN residues=cgetg(nprimes+1,t_VEC);
  GEN moduli=cgetg_copy(residues, &lx), temp;
  //GEN viter=cgetg_copy(residues, &lx);
  for(long i=1;i<=nprimes;i++){
    temp=sqmod_ppower(x,gcoeff(fact,i,1),itos(gcoeff(fact,i,2)),powii(gcoeff(fact,i,1),gcoeff(fact,i,2)),0);
	if(typ(temp)==t_INT){avma = top;return gen_0;}//Cannot do it
	gel(residues,i)=gel(temp,1);
	gel(moduli,i)=gel(temp,2);
	if(lg(gel(temp,1))==3){
	  nsols=2*nsols;//Double sols if there are two solutions, else nothing required to do.
	  //gel(viter,i)=mkvec2(gen_1,gen_2);
	}
	//else gel(viter,i)=mkvec2(gen_1,gen_1);
  }//At this point, we have all the congruences solved.
  GEN allres=cgetg(nsols+1,t_VEC);//The residues
  gel(allres,1)=gel(gel(residues,1),1);//First one
  GEN mastermod=gel(moduli,1), oldmastermod=mastermod;
  long curmaxpos;
  if(lg(gel(residues,1))==3){
    gel(allres,2)=gel(gel(residues,1),2);
	curmaxpos=2;
  }
  else curmaxpos=1;
  for(long i=2;i<=nprimes;++i){//Add the new congruences
    mastermod=mulii(mastermod,gel(moduli,i));
	if(lg(gel(residues,i))==3){//Double
	  for(long j=1;j<=curmaxpos;j++){
		  gel(allres,curmaxpos+j)=Z_chinese_coprime(gel(allres,j), gel(gel(residues,i),1), oldmastermod, gel(moduli,i), mastermod);
		  gel(allres,j)=Z_chinese_coprime(gel(allres,j), gel(gel(residues,i),2), oldmastermod, gel(moduli,i), mastermod);
	  }
	  curmaxpos=curmaxpos*2;
	}
	else{
	  for(long j=1;j<=curmaxpos;++j) gel(allres,j)=Z_chinese_coprime(gel(allres,j), gel(gel(residues,i),1), oldmastermod, gel(moduli,i), mastermod);
	}
	oldmastermod=mastermod;
  }
  GEN rvec=cgetg(3,t_VEC);
  gel(rvec,1)=ZV_sort(allres);
  gel(rvec,2)=icopy(mastermod);
  return gerepileupto(top,rvec);
}

//sqmod, but we check the inputs.
GEN sqmod_typecheck(GEN x, GEN n){
  pari_sp top = avma;
  n=absi(n);
  if(typ(x)!=t_INT && typ(x)!=t_FRAC) pari_err_TYPE("sqmod. x should be rational",x);
  if(typ(n)!=t_INT) pari_err_TYPE("sqmod. n should be integral",n);
  return gerepileupto(top,sqmod(x,n,gen_0));
}

//Solves y^2==x mod p^n=p2n, returns 0 if no solutions and [[residues], modulus] if there are solutions. The modulus will be p^n if p is odd and x is coprime to p; otherwise it may be a smaller power of p.
static GEN sqmod_ppower(GEN x, GEN p, long n, GEN p2n, int iscoprime){
  pari_sp top = avma;
  x=Fp_red(x, p2n);
  if(iscoprime == 0){//Dealing with gcd(x,p)>1 first.
    if(gequal0(x)){
	  GEN tosqrt;
	  if(n%2==0) tosqrt=p2n;
	  else tosqrt=mulii(p,p2n);//The return modulus is sqrt(tosqrt)
      GEN rvec=cgetg(3,t_VEC);//The return vector
	  gel(rvec,1)=mkvec(gen_0);
	  gel(rvec,2)=sqrti(tosqrt);
	  return gerepileupto(top,rvec);
    }
    GEN xnew;
    long v=Z_pvalrem(x, p,&xnew);//x=p^v*xnew, and gcd(xnew,p)=1. We must have v<n-1, else x==0(p^n) which has been treated already.
    if(v%2==1){avma = top;return gen_0;}//Odd power, cannot do it.
    if(v>0){//If v==0, we can pass to the second half.
	  GEN ptov=divii(x,xnew);//Represents p^v
	  GEN copvec=sqmod_ppower(xnew,p,n-v,divii(p2n,ptov),1);//y^2==p^v*xnew mod p^n <==> (y/p^(v/2))^2==xnew mod p^(n-v).
	  if(typ(copvec)==t_INT){avma = top;return gen_0;}//No solution here.
	  //Else, must scale everything by p^(v/2)
	  GEN ptovo2=sqrti(ptov);//p^(v/2);
	  long lx;
	  GEN rvec=cgetg_copy(copvec,&lx);
	  gel(rvec,2)=mulii(gel(copvec,2),ptovo2);//New modulus
	  gel(rvec,1)=cgetg_copy(gel(copvec,1),&lx);
	  for(long i=1;i<lx;i++) gel(gel(rvec,1),i)=mulii(gel(gel(copvec,1),i),ptovo2);
	  return gerepileupto(top, rvec);
    }
  }//Now we have y^2==x mod p^n and x is coprime to p
  GEN root=Zp_sqrt(x,p,n);//Built in method!
  if(root==NULL){avma = top;return gen_0;}//No solution
  //The solutions are now y==+/-root mod p^n, except for p=2 where there are some further considerations.
  GEN rvec;
  if(equalis(p,2)){//Case p=2
    if(n<=3){//Since there IS a solution, these must all boil down to y==1(2).
	  rvec=cgetg(3,t_VEC);
	  gel(rvec,1)=mkvec(gen_1);
	  gel(rvec,2)=gen_2;
	}
	else{//y==+/-root mod 2^(n-1)
	  rvec=cgetg(3,t_VEC);
	  gel(rvec,2)=shifti(p2n,-1);//2^(n-1)
	  gel(rvec,1)=cgetg(3,t_VEC);
	  gel(gel(rvec,1),1)=Fp_red(root,gel(rvec,2));
	  togglesign_safe(&root);
	  gel(gel(rvec,1),2)=Fp_red(root,gel(rvec,2));
	}
  }
  else{//p odd
    rvec=cgetg(3,t_VEC);
	gel(rvec,2)=icopy(p2n);
	gel(rvec,1)=cgetg(3,t_VEC);
	gel(gel(rvec,1),1)=Fp_red(root,p2n);
	togglesign_safe(&root);
	gel(gel(rvec,1),2)=Fp_red(root,p2n);
  }
  return gerepileupto(top, rvec);
}


//INTEGER VECTORS


//Copies an integer vector
GEN ZV_copy(GEN v){
  long len=lg(v);
  GEN rvec=cgetg(len,t_VEC);
  for(long i=1;i<len;++i) gel(rvec,i)=icopy(gel(v,i));
  return rvec;
}

//Checks if the integral vectors v1, v2 are equal, returns 1 if so and 0 else. Segmentation faults will occur if the entries are not integer.
int ZV_equal(GEN v1, GEN v2){
  long len=lg(v1);
  if(lg(v2)!=len) return 0;//Different length
  for(long i=1;i<len;++i){if(!equalii(gel(v1,i),gel(v2,i))) return 0;}
  return 1;
}

//v is a Z-vector, divides v by y, assuming y divides each term exactly
GEN ZV_Z_divexact(GEN v, GEN y){
  long len=lg(v);
  GEN rvec=cgetg(len,t_VEC);
  for(long i=1;i<len;++i) gel(rvec,i)=diviiexact(gel(v,i),y);
  return rvec;
  
}

//v is a Z-vector, multiplies v by the integer x
GEN ZV_Z_mul(GEN v, GEN x){
  long len=lg(v);
  GEN rvec=cgetg(len,t_VEC);
  for(long i=1;i<len;++i) gel(rvec,i)=mulii(gel(v,i),x);
  return rvec; 
}



//TIME


//Returns current time
char *returntime(void){
  time_t rawtime;
  time(&rawtime);
  return asctime(localtime(&rawtime));
}

//Prints current time
void printtime(void){
  pari_printf("%s",returntime());
}


//LISTS


//Circular list of GENs

//Frees the memory pari_malloc'ed by clist
void clist_free(clist *l, long length){
  clist *temp=l;
  long i=1;
  while(i<length){
	temp=l;
	l=l->next;
	pari_free(temp);
	i++;
  }
  pari_free(l);
}

//Put an element before *head_ref, and update *head_ref to point there
void clist_putbefore(clist **head_ref, GEN new_data){
  clist *new_elt = (clist*)pari_malloc(sizeof(clist)); 
  new_elt->data = new_data;
  if(*head_ref!=NULL){
    new_elt->next = *head_ref; 
    new_elt->prev = (*head_ref)->prev;
    (*head_ref)->prev = new_elt;
	(new_elt->prev)->next=new_elt;
  }
  else{
    new_elt->next = new_elt; 
    new_elt->prev = new_elt;
  }
  *head_ref = new_elt;
}

//Put an element after *head_ref, and update *head_ref to point there
void clist_putafter(clist **head_ref, GEN new_data){
  clist *new_elt = (clist*)pari_malloc(sizeof(clist)); 
  new_elt->data = new_data;
  if(*head_ref!=NULL){
    new_elt->prev = *head_ref; 
    new_elt->next = (*head_ref)->next;
    (*head_ref)->next = new_elt;
	(new_elt->next)->prev=new_elt;
  }
  else{
    new_elt->next = new_elt; 
    new_elt->prev = new_elt;
  }
  *head_ref = new_elt;
}

//dir=1 means forward, dir=-1 means backwards. Returns the list as a vector, and makes a clean copy. This also frees the list, but we also need to clean up the list data at the list creation location. The passed in pointer to l should NOT be used as it no longer points to a valid address.
GEN clist_togvec(clist *l, long length, int dir){
  if(l==NULL){//Empty list, return the empty vector.
    GEN rvec=cgetg(1,t_VEC);
    return(rvec);	  
  }
  GEN rvec=cgetg(length+1, t_VEC);
  long lind=1;
  if(dir==1){
    while(lind<=length){
	  gel(rvec,lind)=gcopy(l->data);
	  l=l->next;
	  lind++;
    }
  }
  else{
    while(lind<=length){
	  gel(rvec,lind)=gcopy(l->data);
	  l=l->prev;
	  lind++;
    }
  }
  clist_free(l,length);
  return rvec;
}

//List of GENs

//Frees the memory pari_malloc'ed by glist
void glist_free(glist *l){
  glist *temp=l;
  while(l!=NULL){
	temp=l;
	l=l->next;
	pari_free(temp);
  }
}

//Removes the last element of the glist and returns it without copying
GEN glist_pop(glist **head_ref){
  glist *temp=*head_ref;
  GEN x=temp->data;
  *head_ref=temp->next;
  pari_free(temp);
  return x;
}

//Put an element at the start of the glist
void glist_putstart(glist **head_ref, GEN new_data){
  glist *new_elt = (glist*)pari_malloc(sizeof(glist)); 
  new_elt->data = new_data; 
  new_elt->next = *head_ref; 
  *head_ref = new_elt; 

}

//dir=1 means forward, dir=-1 means backwards. Returns the list as a vector, makes a clean copy. This also frees the list, but we also need to clean up the list data at the list creation location. The passed in pointer to l should NOT be used as it no longer points to a valid address.
GEN glist_togvec(glist *l, long length, int dir){
	glist *lcopy=l;
	GEN rvec=cgetg(length+1, t_VEC);
	if(dir==1){
	  long lind=1;
	  while(l!=NULL && lind<=length){
		  gel(rvec,lind)=gcopy(l->data);
		  l=l->next;
		  lind++;
	  }
	  if(lind<=length){//Couldn't finish.
	    pari_err(e_MISC,"List length is too long");
	  }
	}
	else{
      long lind=length;
	  while(l!=NULL && lind>0){
		gel(rvec,lind)=gcopy(l->data);
		l=l->next;
		lind--;
	  }
	  if(lind>0){//Couldn't finish.
	    pari_err(e_MISC,"List length is too long");
	  }
	}
	glist_free(lcopy);
	return rvec;
}


//List of longs

//Frees the memory pari_malloc'ed by llist
void llist_free(llist *l){
  llist *temp=l;
  while(l!=NULL){
	temp=l;
	l=l->next;
	pari_free(temp);
  }
}

//Removes the last element of the llist and returns it
long llist_pop(llist **head_ref){
  llist *temp=*head_ref;
  long x=temp->data;
  *head_ref=temp->next;
  pari_free(temp);
  return x;
}

//Put an element at the start of the llist
void llist_putstart(llist **head_ref, long new_data){
  llist *new_elt = (llist*)pari_malloc(sizeof(llist)); 
  new_elt->data = new_data; 
  new_elt->next = *head_ref; 
  *head_ref = new_elt; 
}

//dir=1 means forward, dir=-1 means backwards. Returns the list as a vector. This also frees the list. The passed in pointer to l should NOT be used as it no longer points to a valid address.
GEN llist_togvec(llist *l, long length, int dir){//No garbage collection necessary with longs!
  llist *lcopy=l;
  GEN rvec=cgetg(length+1, t_VEC);
  if(dir==1){
    long lind=1;
    while(l!=NULL && lind<=length){
	  gel(rvec,lind)=stoi(l->data);
	  l=l->next;
	  lind++;
    }
    if(lind<=length){//Couldn't finish.
      pari_err(e_MISC,"List length is too long");
    }
  }
  else{
    long lind=length;
    while(l!=NULL && lind>0){
      gel(rvec,lind)=stoi(l->data);
	  l=l->next;
	  lind--;
    }
    if(lind>0){//Couldn't finish.
      pari_err(e_MISC,"List length is too long");
    }
  }
  llist_free(lcopy);
  return(rvec);
}

//dir=1 means forward, dir=-1 means backwards. Returns the list as a VECSMALL. This also frees the list. The passed in pointer to l should NOT be used as it no longer points to a valid address.
GEN llist_tovecsmall(llist *l, long length, int dir){//No garbage collection necessary with longs!
  llist *lcopy=l;
  GEN rvec=cgetg(length+1, t_VECSMALL);
  if(dir==1){
    long lind=1;
    while(l!=NULL && lind<=length){
	  rvec[lind]=l->data;
	  l=l->next;
	  lind++;
    }
    if(lind<=length){//Couldn't finish.
      pari_err(e_MISC,"List length is too long");
    }
  }
  else{
    long lind=length;
    while(l!=NULL && lind>0){
      rvec[lind]=l->data;
	  l=l->next;
	  lind--;
    }
    if(lind>0){//Couldn't finish.
      pari_err(e_MISC,"List length is too long");
    }
  }
  llist_free(lcopy);
  return(rvec);
}
