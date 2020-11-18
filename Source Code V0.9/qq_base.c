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

//STATIC METHOD DECLARATIONS
static GEN sqmod_ppower(GEN x, GEN p, long n, GEN p2n, int iscoprime);
static int opp_gcmp(void *data, GEN x, GEN y);
static GEN quadraticinteger(GEN A, GEN B, GEN C);

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
GEN mat_eval_tc(GEN M, GEN x){
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
GEN lin_intsolve_tc(GEN A, GEN B, GEN n){
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
GEN mat3_complete_tc(GEN A, GEN B, GEN C){
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



//LINEAR ALGEBRA MODULO N

//Returns the vector of [eigenvalue, [eigenvectors]] of an FpM, where p is PRIME
GEN FpM_eigenvecs(GEN M, GEN p){
  pari_sp top=avma;
  GEN pol=FpM_charpoly(M, p);
  GEN rts=FpX_roots(pol, p);
  long nrts=lg(rts);
  if(nrts==1){avma=top;return cgetg(1, t_VEC);}//No roots
  GEN shiftM=cgetg(nrts, t_VEC), id=matid(4);
  for(long i=1;i<nrts;i++) gel(shiftM, i)=FpM_sub(M, FpM_Fp_mul(id, gel(rts, i), p), p);//Stores M-eval*Id
  GEN ret=cgetg(nrts, t_VEC);
  for(long i=1;i<nrts;i++){
	gel(ret, i)=cgetg(3, t_VEC);
	gel(gel(ret, i), 1)=icopy(gel(rts, i));//Eigenvalue
	gel(gel(ret, i), 2)=FpM_ker(gel(shiftM, i), p);//Eigenvectors
  }
  return gerepileupto(top, ret);
}


//SOLVING EQUATIONS MOD N



//I should update this so that we use something like forvec or w/e
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
GEN sqmod_tc(GEN x, GEN n){
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
  long lx;
  GEN rvec=cgetg_copy(v, &lx);
  for(long i=1;i<lx;i++) gel(rvec, i)=diviiexact(gel(v, i), y);
  return rvec;
  
}

//v is a Z-vector, multiplies v by the integer x
GEN ZV_Z_mul(GEN v, GEN x){
  long lx;
  GEN rvec=cgetg_copy(v, &lx);
  for(long i=1;i<lx;i++) gel(rvec, i)=mulii(gel(v, i), x);
  return rvec; 
}



//RANDOM



//Returns a random element from the vector v. NOT STACK CLEAN
GEN rand_elt(GEN v){
  if(typ(v)!=t_VEC) pari_err_TYPE("Not a vector", v);
  long i=rand_l(lg(v)-1);
  return gel(v,i);
}

//Returns random long from 1 to len
long rand_l(long len){
  pari_sp top=avma;
  GEN ind=randomi(stoi(len));
  long ret=itos(ind)+1;
  avma=top;
  return ret;
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

//Appends l to the end of v, returning a clean copy. dir=-1 means forward, dir=-1 backward. This also frees the list, but we also need to clean up the list data at the list creation location. The passed in pointer to l should NOT be used as it no longer points to a valid address.
GEN glist_togvec_append(glist *l, GEN v, long length, int dir){
	glist *lcopy=l;
	long vlen=lg(v), rveclen=length+vlen;
	GEN rvec=cgetg(rveclen, t_VEC);
	for(long i=1;i<vlen;i++) gel(rvec, i)=gcopy(gel(v, i));//Copying v
	if(dir==1){
	  long lind=vlen;
	  while(l!=NULL && lind<rveclen){
		  gel(rvec,lind)=gcopy(l->data);
		  l=l->next;
		  lind++;
	  }
	  if(lind<rveclen){//Couldn't finish.
	    pari_err(e_MISC,"List length is too long");
	  }
	}
	else{
      long lind=rveclen-1;
	  while(l!=NULL && lind>=vlen){
		gel(rvec,lind)=gcopy(l->data);
		l=l->next;
		lind--;
	  }
	  if(lind>=vlen){//Couldn't finish.
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



//SHORT VECTORS IN LATTICES

//We follow the article "Improved Methods for Calculating Vectors of Short Length in a Lattice, Including Complexity Analysis" by Fincke and Pohst (Mathematics of Computation, Vol. 44, No. 170 (Apr., 1985), pp 463-471


//Follows Algorithm 2.12 in Fincke-Pohst. If rdataonly=1, returns the data required to run the algorithm (useful if you will use multiple times only changing C).
GEN lat_smallvectors(GEN A, GEN C1, GEN C2, GEN condition, int onesign, int isintegral, int rdataonly, long prec){
  pari_sp top=avma, mid;
  long np1=lg(A);//n+1
  GEN R=mat_choleskydecomp(A, 0, prec), Rinv=ginv(R);//Step 1
  GEN Rinvtrans=shallowtrans(Rinv);//We reduce R with LLL, but this works on the columns and we care about the rows.
  GEN Uinvtrans=lllfp(Rinvtrans, 0.99, 0);//LLL reduce Rinvtrans
  GEN Uinv=shallowtrans(Uinvtrans);
  mid=avma;
  GEN Sinv=gmul(Uinv, Rinv);
  GEN U=ginv(Uinv);
  GEN S=gmul(R, U);
  gerepileallsp(top, mid, 3, &Sinv, &U, &S);//Cleaning up. We are at the end of Step 3, and know S, S^-1 and U (R is forgotten but not needed anymore).
  mid=avma;
  GEN Sprimesizes=cgetg(np1, t_VEC);//The sizes of the rows of S^-1 (will store ||s'_i||, which is the ith column of S^(-1)^T, i.e. the ith row of S^-1.
  for(long i=1;i<np1;i++){
	gel(Sprimesizes, i)=gen_0;
	for(long j=1;j<np1;j++) gel(Sprimesizes, i)=gadd(gel(Sprimesizes, i), gsqr(gcoeff(Sinv, i, j)));//Adding up
  }
  GEN perm=gerepileupto(mid, gen_indexsort(Sprimesizes, NULL, &opp_gcmp));//Sort Sprimesizes from large to small, returns VECSMALL. This is pi in Step 3.
  GEN perminv=perm_inv(perm);//pi^(-1)
  mid=avma;
  GEN Snew=vecpermute(S, perminv);//S with the columns permuted.
  GEN SnewtransSnew=gmul(shallowtrans(Snew), Snew);//S^(Tr)*S
  if(isintegral==1) SnewtransSnew=gerepileupto(mid, gdivgs(ground(gmulgs(SnewtransSnew, 2)), 2));//The entries of SnewtransSnew are half integers, so we double it, round (should generally have enough precision that this is correct), and divide by 2 again to make the entries exact.
  GEN chol=mat_choleskydecomp(SnewtransSnew, 1, prec), newcond;//Cholesky on Snew
  if(!gequal0(condition)){//Must update the condition.
    mid=avma;
	GEN newcondmat=gmul(shallowtrans(U), gmul(gel(condition, 1), U));//M->U^T*M*U
	newcondmat=vecpermute(newcondmat, perm);//Permute the columns by perm, ie *W on the right (W is the permutation matrix producing the final vector); we are now solving for y, and x=U*W*y.
	newcond=cgetg(3, t_VEC);//Remaking so as to not disrupt previous memory
	gel(newcond, 1)=rowpermute(newcondmat, perm);//Permute the rows by perm, i.e. *W^T on the left
	gel(newcond, 2)=gcopy(gel(condition, 2));
  }
  else newcond=gen_0;
  if(rdataonly==0) return gerepileupto(top, lat_smallvectors_givendata(chol, U, perminv, C1, C2, newcond, onesign, prec));
  GEN ret=cgetg(5, t_VEC);//The crucial data is chol, U, perminv, condition.
  gel(ret, 1)=gcopy(chol);
  gel(ret, 2)=gcopy(U);
  gel(ret, 3)=gcopy(perminv);
  gel(ret, 4)=gcopy(condition);
  return gerepileupto(top, ret);
}

//Given the data, finds the small vectors. Use lat_smallvectors to generate [chol, U, perminv, condition], and feed that into here along with C1, C2, onesign, and prec to get the small vectors.
GEN lat_smallvectors_givendata(GEN chol, GEN U, GEN perminv, GEN C1, GEN C2, GEN condition, int onesign, long prec){
  pari_sp top=avma;
  GEN yvals=lat_smallvectors_cholesky(chol, C1, C2, condition, onesign, prec);
  for(long i=1;i<lg(yvals);i++) gel(gel(yvals, i), 1)=vecpermute(gel(gel(yvals, i), 1), perminv);//Permuting the entries of y, now we just need to apply U to get back to x
  long lx;
  GEN ret=cgetg_copy(yvals, &lx);//The return. Must shift back.
  for(long i=1;i<lx;i++){
	gel(ret, i)=cgetg(3, t_VEC);
	gel(gel(ret, i), 1)=gmul(U, gel(gel(yvals, i), 1));//U*y
	gel(gel(ret, i), 2)=gcopy(gel(gel(yvals, i), 2));//The value
  }
  return gerepileupto(top, ret); 
}

//lat_smallvectors with typechecking We do not allow for a condition to be passed in, and if C2=0 we replace (C1, C2] by (0, C1].
GEN lat_smallvectors_tc(GEN A, GEN C1, GEN C2, int onesign, int isintegral, long prec){
  if(typ(A)!=t_MAT) pari_err_TYPE("Please input a nxn symmetric MATRIX", A);
  if(gequal0(C2)) return lat_smallvectors(A, gen_0, C1, gen_0, onesign, isintegral, 0, prec);
  return lat_smallvectors(A, C1, C2, gen_0, onesign, isintegral, 0, prec);
}

//Q is the Cholesky decomposition of a matrix A, this computes all vectors x such that C_1<x^T*A*x<=C2 (i.e. C_1<Q(x)<=C_2 where Q(x)=sum(i=1..n)q_ii(x_i+sum(j=i+1..n)q_ijxj)^2. Algorithm 2.8 of Fincke Pohst. We can also pass in a condition, which is either 0 (no condition) or [M, n] where x^T*M*x=n must also be satisfied. M and n must be integral. (the point: M gives an indefinite norm condition on the vector, and A combines this norm with other info to make a positive definite form. We use the condition when finding small norm 1 elements of a quaternion algebra.)
GEN lat_smallvectors_cholesky(GEN Q, GEN C1, GEN C2, GEN condition, int onesign, long prec){
  pari_sp top=avma, mid;
  long np1=lg(Q), n=np1-1;//Number of variables+1 and n
  GEN T=zerovec(n);//Stores the ''tail'' of x, as we work from the back to the front (x_m to x_1)
  GEN U=zerovec(n);//U[i] will store the sum of j=i+1 to n of q_{ij}x_j
  GEN UB=zerovec(n);//UB represents the upper bould for x_i
  GEN x=zerocol(n), Z, K, cond, temp;//x represents the solution
  long i=np1-1, count=0;//i represents the current index, initially set to n. count=the number of solutions
  gel(T, n)=gcopy(C2);//initialize T[n]=0
  gel(U, n)=gen_0;//Clearly U[n]=0
  int step=2;//Represents the current step of algorithm 2.8
  int C1isnot0=1-gequal0(C1), condisnot0=1-gequal0(condition), xpass0=0;
  long uvar=0;//In keeping track of the condition, we use a variable.
  GEN u=gen_0, x1sols;
  if(condisnot0){
	uvar=fetch_var();//Creating temporary variable
    u=pol_x(uvar);//The monomial
  }
  glist *S=NULL;//Pointer to the list start
  GEN v=cgetg(1, t_VEC);//The list, is used for garbage collection partway through
  while(step>0){
	if(gc_needed(top, 1)){
	  mid=avma;
	  if(onesign==0) count=2*count;//We found double the solutions in this case.
	  v=glist_togvec_append(S, v, count, 1);
	  count=0;
	  S=NULL;
	  T=gcopy(T);
	  U=gcopy(U);
	  UB=gcopy(UB);
	  x=gcopy(x);
	  gerepileallsp(top, mid, 5, &v, &T, &U, &UB, &x);
	  if(condisnot0) u=pol_x(uvar);//Resetting u
	}
	if(step==2){
	  Z=gsqrt(gabs(gdiv(gel(T, i), gcoeff(Q, i, i)), prec), prec);//The inner square root should be positive always, but could run into issue if T=0 and rounding puts it <0. Z=sqrt(T[i]/Q[i,i])
	  gel(UB, i)=gfloor(gsub(Z, gel(U, i)));//UB[i]=floor(Z-U[i]);
	  gel(x, i)=gsubgs(gceil(gneg(gadd(Z, gel(U, i)))), 1);//x[i]=ceil(-Z-U[i])-1;
	  step=3;
	}
    if(step==3){
	  gel(x, i)=gaddgs(gel(x, i), 1);//x[i]=x[i]+1
	  if(gcmp(gel(x, i), gel(UB, i))<=0) step=5;//If x[i]<=UB[i], goto step 5
	  else step=4; //If x[i]>UB[i], goto step 4
	}
	if(step==4){
	  i=i+1;
	  step=3;
	  continue;//May as well go back to start
	}
	if(step==5){
	  if(i==1) step=6;
	  else{
		i=i-1;
		gel(U, i)=gen_0;
		for(long j=i+1;j<np1;j++) gel(U, i)=gadd(gel(U, i), gmul(gcoeff(Q, i, j), gel(x, j)));//U[i]=sum(j=i+1,n,q[i,j]*x[j]);
		gel(T, i)=gsub(gel(T, i+1), gmul(gcoeff(Q, i+1, i+1), gsqr(gadd(gel(x, i+1), gel(U, i+1)))));//T[i]=T[i+1]-q[i+1,i+1]*(x[i+1]+U[i+1])^2;
		if(condisnot0 && i==1){step=7;}//We have a condition to deal with!
		else{//Go back now
		  step=2;
		  continue;
		}
	  }
	}
	if(step==6){//Solution found
	  if(gequal0(x)){step=0;continue;}//Terminate
	  K=gadd(gsub(C2, gel(T, 1)), gmul(gcoeff(Q, 1, 1), gsqr(gadd(gel(x, 1), gel(U, 1)))));//K=C-T[1]+q[1,1]*(x[1]+U[1])^2;
	  if(C1isnot0){
		if(gcmp(C1, K)!=-1){//Result was <C1. Since K is symmetric in x[1] about U[1], we can swap x to the other side.
		  if(gcmp(gel(x, 1), gel(U, 1))==-1){//-2U[1]-x[1]; i.e. going from -U[1]-V to -U[1]+V.
			temp=gfloor(gsub(gmulgs(gel(U, 1), -2), gel(x, 1)));//We include the if statement to account for rounding errors; don't want an infinite loop if we keep swapping back to the same value of gel(x, 1)!
			if(gsigne(temp)>=0){
			  gel(x, 1)=gen_0;
			  if(gequal0(x)){step=0;continue;}//If flipping over 0, we must check that we don't pass the stop condition!
			}
			gel(x, 1)=temp;
		  }
		  step=3;
		  continue;
		}
	  }
	  glist_putstart(&S, mkvec2copy(x, K));
	  count++;
	  if(onesign==0) glist_putstart(&S, mkvec2copy(gneg(x), K));
	  step=3;
	}
	if(step==7){//Dealing with extra condtions
	  gel(x, 1)=u;
	  cond=gmul(gmul(shallowtrans(x), gel(condition, 1)), x);
	  x1sols=quadraticinteger(polcoef_i(cond, 2, uvar), polcoef_i(cond, 1, uvar), subii(polcoef_i(cond, 0, uvar), gel(condition, 2)));
	  gel(x, 1)=gen_0;
	  if(gequal0(x)) xpass0=1;//This is the last check
	  for(long j=1;j<lg(x1sols);j++){
		K=gadd(gsub(C2, gel(T, 1)), gmul(gcoeff(Q, 1, 1), gsqr(gadd(gel(x1sols, j), gel(U, 1)))));
		if(gcmp(K, C2)==1) continue;//>C2
		if(C1isnot0){if(gcmp(C1, K)!=-1) continue;}//<=C1
		if(xpass0 && signe(gel(x1sols, j))!=-1) continue;//x is 0 (except the first coefficient), so the first coefficent has to be negative.
		gel(x, 1)=gel(x1sols, j);//Now we are good, all checks out.
		glist_putstart(&S, mkvec2copy(x, K));
	    count++;
	    if(onesign==0) glist_putstart(&S, mkvec2copy(gneg(x), K));
	  }
	  if(xpass0){step=0;continue;}//Game over, we are done!
	  i=2;
	  step=3;
	  continue;
	}
  }
  if(condisnot0) delete_var();//Delete the created variable.
  if(onesign==0) count=2*count;//We found double the solutions in this case.
  return gerepileupto(top, glist_togvec_append(S, v, count, -1));
}

//Solves Ax^2+Bx+C=0 in the integers and returns the solution (A, B, C need to be integers, with at least one non-zero).
static GEN quadraticinteger(GEN A, GEN B, GEN C){
  pari_sp top=avma;
  if(gequal0(A)){//Actually a linear. This case occurs when dealing with small vectors in the quaternion algebra ramified nowhere.
	if(gequal0(B)) return cgetg(1, t_VEC);//We say there are no soln's when A=B=0
	GEN x=Qdivii(C, B);
	if(typ(x)!=t_INT){avma=top;return cgetg(1, t_VEC);}//No solution!
	GEN ret=cgetg(2, t_VEC);
	gel(ret, 1)=negi(x);
	return gerepileupto(top, ret);
  }
  GEN disc=subii(sqri(B), mulii(A, shifti(C, 2)));//B^2-4AC
  GEN rt;
  long israt=Z_issquareall(disc, &rt);
  if(!israt){avma=top;return cgetg(1, t_VEC);}//Disc is not a square, no solutions.
  GEN x1, x2, sols;
  if(gequal0(rt)){
	x1=Qdivii(gneg(B), shifti(A, 1));//-B/2A
	if(typ(x1)!=t_INT){avma=top;return cgetg(1, t_VEC);}//No solution
	sols=cgetg(2, t_VEC);
	gel(sols, 1)=icopy(x1);
	return gerepileupto(top, sols);
  }
  //Now there could be 2 roots
  GEN denom=shifti(A, 1);//2A;
  x1=Qdivii(subii(rt, B), denom);//(sqrt(D)-B)/2A
  x2=Qdivii(subii(negi(rt), B), denom);//(-sqrt(D)-B)/2A
  if(typ(x1)==t_INT){
	if(typ(x2)==t_INT){
	  sols=cgetg(3, t_VEC);
	  gel(sols, 1)=icopy(x1);
	  gel(sols, 2)=icopy(x2);
	}
	else{
	  sols=cgetg(2, t_VEC);
	  gel(sols, 1)=icopy(x1);
	}
  }
  else{
	if(typ(x2)==t_INT){
	  sols=cgetg(2, t_VEC);
	  gel(sols, 1)=icopy(x2);
	}
	else{
	  avma=top;
	  return cgetg(1, t_VEC);
	}
  }
  return gerepileupto(top, sols);
}

//Returns -gcmp(x, y), and used for sorting backwards.
static int opp_gcmp(void *data, GEN x, GEN y){return -gcmp(x, y);}

//Computes the Cholesky decomposition of A (nxn symmetric matrix). In otherwords, finds R such that R^T*R=A. If rcoefs=0, returns R, otherwise returns the nxn matrix B so that x^TAx is expressible as sum(i=1..n)b_ii(x_i+sum(j=i+1..n)b_ijxj)^2.
GEN mat_choleskydecomp(GEN A, int rcoefs, long prec){
  pari_sp top=avma;
  long n=lg(A)-1;//A is nxn
  GEN M=gcopy(A);//Will be manipulating the entries, so need to copy A.
  for(long i=1;i<n;i++){
	for(long j=i+1;j<=n;j++){
	  gcoeff(M, j, i)=gcopy(gcoeff(M, i, j));//M[j,i]=M[i,j]
	  gcoeff(M, i, j)=gdiv(gcoeff(M, i, j), gcoeff(M, i, i));//M[i,j]=M[i,j]/M[i,i]
	}
	for(long j=i+1;j<=n;j++){
	  for(long k=j;k<=n;k++) gcoeff(M, j, k)=gsub(gcoeff(M, j, k), gmul(gcoeff(M, j, i), gcoeff(M, i, k)));//M[j,k]=M[j,k]-M[j,i]*M[i,k];
	}
  }
  if(rcoefs==1){
	GEN ret=cgetg_copy(M, &n);//M stores the coeff, but we should delete the lower diagonal
	for(long i=1;i<n;i++){//Column i
      gel(ret, i)=cgetg(n, t_COL);
	  for(long j=1;j<=i;j++) gcoeff(ret, j, i)=gcopy(gcoeff(M, j, i));
	  for(long j=i+1;j<n;j++) gcoeff(ret, j, i)=gen_0;
	}
	return gerepileupto(top, ret);
  }
  for(long i=1;i<=n;i++){
	for(long j=1;j<i;j++) gcoeff(M, i, j)=gen_0;
	gcoeff(M, i, i)=gsqrt(gcoeff(M, i, i), prec);
	for(long j=i+1;j<=n;j++) gcoeff(M, i, j)=gmul(gcoeff(M, i, j), gcoeff(M, i, i));
  }
  return gerepilecopy(top, M);
}

//choleskydecomp with typechecking
GEN mat_choleskydecomp_tc(GEN A, int rcoefs, long prec){
  if(typ(A)!=t_MAT) pari_err_TYPE("Please input a nxn symmetric MATRIX", A);
  if(lg(A)!=lg(gel(A, 1))) pari_err_TYPE("Please input a nxn symmetric matrix", A);
  return mat_choleskydecomp(A, rcoefs, prec);
}

//Finds unimodular S so that M'=S*M, with |M'[i,j]|<=1/2*M'[j,j] for j>i. Returns [M', S, S^(-1)].
GEN mat_uptriag_rowred(GEN M){
  pari_sp top=avma;
  long n=lg(M)-1;
  GEN A=gcopy(M), S=matid(n), Sinv=matid(n), u;//Copy M, S=Sinv=Id
  for(long i=2;i<=n;i++){
	for(long j=1;j<i;j++){
	  u=ground(gdiv(gneg(gcoeff(A, j, i)), gcoeff(A, i, i)));//u=round(-A[j,i]/A[i,i])
	  for(long k=1;k<=n;k++){
		gcoeff(A, j, k)=gadd(gcoeff(A, j, k), gmul(u, gcoeff(A, i, k)));//A[j,]=A[j,]+u*A[i,];
		gcoeff(Sinv, k, i)=gsub(gcoeff(Sinv, k, i), gmul(u, gcoeff(Sinv, k, j)));//Sinv[,i]=Sinv[,i]-u*Sinv[,j];
		gcoeff(S, j, k)=gadd(gcoeff(S, j, k), gmul(u, gcoeff(S, i, k)));//S[j,]=S[j,]+u*S[i,];
	  }
	}
  }
  GEN ret=cgetg(4, t_VEC);
  gel(ret, 1)=gcopy(A);
  gel(ret, 2)=gcopy(S);
  gel(ret, 3)=gcopy(Sinv);
  return gerepileupto(top, ret);
}

//mat_uptriag_rowred with typechecking
GEN mat_uptriag_rowred_tc(GEN M){
  if(typ(M)!=t_MAT) pari_err_TYPE("Please input a nxn upper triangular MATRIX", M);
  if(lg(M)!=lg(gel(M, 1))) pari_err_TYPE("Please input a nxn  upper triangular matrix", M);
  return mat_uptriag_rowred(M);
}






