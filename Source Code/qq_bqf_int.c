#ifndef PARILIB
#define PARILIB
#include <pari/pari.h>
#endif

#ifndef METHDECL
#define METHDECL
#include "qquadraticdecl.h"
#endif

//STATIC DECLARATIONS
static GEN bqf_iform(GEN q1, GEN q2);
static GEN ibqf_int_reverseriver(GEN r);
static GEN ibqf_intRS_splitindices(GEN river, GEN ind);
static void ibqf_intformsRS_byriver_indices(GEN r1, GEN r2, llist **inds1, llist **inds2, llist **loverlap, long *inum, int data);



//INTERSECTION DATA



//Returns B_delta(q1,q2)
GEN bqf_bdelta(GEN q1, GEN q2){
  pari_sp top=avma;
  GEN B1B2=mulii(gel(q1,2),gel(q2,2));
  GEN tA1C3=shifti(mulii(gel(q1,1),gel(q2,3)),1);
  GEN tA3C1=shifti(mulii(gel(q1,3),gel(q2,1)),1);
  return gerepileupto(top,subii(B1B2,addii(tA1C3,tA3C1)));
}

//bqf_bdelta with typechecking
GEN bqf_bdelta_tc(GEN q1, GEN q2){
  bqf_check(q1);
  bqf_check(q2);
  return bqf_bdelta(q1,q2);
}

//Returns intersection level of q1, q2
GEN bqf_intlevel(GEN q1, GEN q2){
  pari_sp top=avma;
  GEN mA1B2=mulii(gel(q1,1),gel(q2,2));togglesign_safe(&mA1B2);
  GEN A2B1=mulii(gel(q2,1),gel(q1,2));
  GEN A=addii(mA1B2,A2B1);
  GEN mtA1C2=shifti(mulii(gel(q1,1),gel(q2,3)),1);togglesign_safe(&mtA1C2);
  GEN tA2C1=shifti(mulii(gel(q2,1),gel(q1,3)),1);
  GEN B=addii(mtA1C2,tA2C1);
  GEN mB1C2=mulii(gel(q1,2),gel(q2,3));togglesign_safe(&mB1C2);
  GEN B2C1=mulii(gel(q2,2),gel(q1,3));
  GEN C=addii(mB1C2,B2C1);
  //Now [A,B,C]=[-A1B2 + A2B1,-2A1C2 + 2A2C1,-B1C2 + B2C1], so need to find gcd and multiply by the sign of the first coefficient.
  int s=signe(A);
  GEN g=gcdii(A,B);
  g=gerepileupto(top,gcdii(g,C));
  if(s==-1) togglesign_safe(&g);
  return g;
}

//bqf_intlevel with typechecking
GEN bqf_intlevel_tc(GEN q1, GEN q2){
  bqf_check(q1);
  bqf_check(q2);
  return bqf_intlevel(q1,q2);
}

//Each pair in pairs is assumed to have the first qf similar to q. This translates them to [q,q'].
GEN ibqf_intpairs_transtoq(GEN pairs, GEN q, GEN rootD){
  pari_sp top=avma;
  long lx;
  GEN tmats=cgetg_copy(pairs,&lx);
  for(long i=1;i<lg(pairs);i++) gel(tmats,i)=ibqf_isequiv_tmat(gel(gel(pairs,i),1), q, rootD);//Transition matrices
  GEN rvec=cgetg_copy(pairs,&lx);
  for(long i=1;i<lg(pairs);i++){
	gel(rvec,i)=cgetg_copy(gel(pairs,1),&lx);
	gel(gel(rvec,i),1)=ZV_copy(q);
	gel(gel(rvec,i),2)=bqf_trans(gel(gel(pairs,i),2),gel(tmats,i));//Translating the pairs
  }
  return gerepileupto(top,rvec);
}

//Outputs the upper half plane intersection point of q1, q2. If location=0, it is the intersection of q1, q2; if location=1, we translate it to the fundamental domain of PSL(2,Z); if imag(location)!=0, then location is assumed to be a point on \ell_q1. We translate the intersection point to the geodesic between location and gamma_q1(location). If the invariant automorph is large, then we need to increase the precision to ensure accurate results.
GEN ibqf_intpoint(GEN q1, GEN q2, GEN location, GEN autom){
  pari_sp top=avma;
  GEN form=bqf_iform(q1,q2);
  GEN g=ZV_content(form);
  form=ZV_Z_divexact(form,g);
  if(signe(gel(form,1))==-1) ZV_togglesign(form);
  GEN D=bqf_disc(form);
  if(signe(D)!=-1) pari_err_TYPE("q1 and q2 do not intersect!",form);
  GEN w=quadroot(D);
  if(gequal0(location)) return gerepilecopy(top, gel(bqf_roots(form, D, w),1));//Intersection as is
  if(gequal(location,gen_1)) return gerepilecopy(top, gel(bqf_roots(dbqf_red(form), D, w),1));//Fundamental domain
  GEN autominv;
  if(signe(gel(q1,1))==-1){
    autominv=autom;
	autom=ZM_inv(autominv,NULL);
  }
  else autominv=ZM_inv(autom,NULL);//We want applying autom to approach the larger root. The automorph approaches the first root, so we must swap it up if the first root is smaller, i.e. gel(q1,1)<0.
  GEN z=gel(bqf_roots(form, D, w),1);
  GEN rloc=greal(location);
  while(gcmp(greal(z),rloc)>=0){
	z=mat_eval(autominv,z);
	form=bqf_trans(form,autom);
  }
  while(gcmp(greal(z),rloc)<0){
	z=mat_eval(autom,z);
	form=bqf_trans(form,autominv);
  }
  return gerepilecopy(top,gel(bqf_roots(form, D, w),1));//We update the form to increase precision of final result.
}

//ibqf_intpoint with typechecking
GEN ibqf_intpoint_tc(GEN q1, GEN q2, GEN location){
  pari_sp top=avma;
  bqf_check(q1);
  bqf_check(q2);
  if(gequal0(location)) return ibqf_intpoint(q1, q2, gen_0, gen_0);
  if(gequal(location,gen_1)) return ibqf_intpoint(q1, q2, gen_1, gen_0);
  if(typ(location)!=t_COMPLEX || gsigne(gimag(location))!=1) pari_err_TYPE("please supply either 0, 1, or a complex number with positive imaginary part",location);
  GEN aut=ibqf_automorph_D(q1,bqf_disc(q1));
  return gerepileupto(top, ibqf_intpoint(q1, q2, location, aut));
}

// For qi=[Ai,Bi,Ci], this returns [-A1B2 + A2B1,-2A1C2 + 2A2C1,-B1C2 + B2C1]
static GEN bqf_iform(GEN q1, GEN q2){
  pari_sp top=avma;
  GEN mA1B2=mulii(gel(q1,1),gel(q2,2));togglesign_safe(&mA1B2);
  GEN A2B1=mulii(gel(q2,1),gel(q1,2));
  GEN mtA1C2=shifti(mulii(gel(q1,1),gel(q2,3)),1);togglesign_safe(&mtA1C2);
  GEN tA2C1=shifti(mulii(gel(q2,1),gel(q1,3)),1);
  GEN mB1C2=mulii(gel(q1,2),gel(q2,3));togglesign_safe(&mB1C2);
  GEN B2C1=mulii(gel(q2,2),gel(q1,3));
  long lx;
  GEN rvec=cgetg_copy(q1,&lx);
  gel(rvec,1)=addii(mA1B2,A2B1);
  gel(rvec,2)=addii(mtA1C2,tA2C1);
  gel(rvec,3)=addii(mB1C2,B2C1);
  return gerepileupto(top, rvec);
}



//INTERSECTION NUMBER COMPUTATION



//Given the output of ibqf_river_positions(qi), this computes the full intersection of q1 and q2 2(RS+RO)
GEN ibqf_int(GEN r1, GEN r2){
  pari_sp top=avma;
  GEN inumRS=ibqf_intRS_byriver(r1,r2);
  GEN r2rev=ibqf_int_reverseriver(r2);
  GEN inumRO=ibqf_intRS_byriver(r1,r2rev);
  return gerepileupto(top,shifti(addii(inumRS,inumRO),1));
}

//Given the output r of ibqf_river_positions(q), this reverses the river, i.e. computes the river corresponding to -q.
static GEN ibqf_int_reverseriver(GEN r){
  long lx;
  GEN rvec=cgetg(4,t_VEC);
  gel(rvec,1)=cgetg_copy(gel(r,2),&lx);
  gel(rvec,2)=cgetg_copy(gel(r,1),&lx);
  gel(rvec,3)=cgetg_copy(gel(r,3),&lx);
  long i1=1, i2=1, ibot=1;
  for(long itop=lg(gel(r,3))-1;itop>0;itop--){
	if(gel(r,3)[itop]==0){
	  gel(rvec,3)[ibot]=1;
	  gel(rvec,2)[i2]=ibot;
	  i2++;
	}
	else{
	  gel(rvec,3)[ibot]=0;
	  gel(rvec,1)[i1]=ibot;
	  i1++;
	}
	ibot++;
  }
  return rvec;
}

//ibqf_int with typecheck
GEN ibqf_int_tc(GEN q1, GEN q2, long prec){
  pari_sp top=avma;
  GEN D1=bqf_checkdisc(q1);
  if(signe(D1)!=1) pari_err_TYPE("please supply an indefinite binary quadratic form",q1);
  GEN D2=bqf_checkdisc(q2);
  if(signe(D2)!=1) pari_err_TYPE("please supply an indefinite binary quadratic form",q2);
  GEN r1data=ibqf_river_positions(q1, gsqrt(D1,prec));
  GEN r2data=ibqf_river_positions(q2, gsqrt(D2,prec));
  return gerepileupto(top,ibqf_int(r1data, r2data));
}

//Given the output of ibqf_river_positions(qi), this computes the RS intersection of q1 and q2.
GEN ibqf_intRS_byriver(GEN r1, GEN r2){
  pari_sp top=avma;
  GEN split;//For splitting the indices
  llist *iseq=NULL;//Stores the overlap in sequences
  glist *ind1LR=NULL;//Stores the indices of r1 which give 0/1 and preceeded by iseq
  glist *ind2LR=NULL;//Stores the indices of r2 which give 0/1 and preceeded by iseq
  GEN *r1inds, *r2inds;//Pointers to the current set of indices we are looking at
  split=ibqf_intRS_splitindices(gel(r1,3),gel(r1,1));//r1 starts with 0
  r1inds=&gel(split,1);
  split=ibqf_intRS_splitindices(gel(r2,3),gel(r2,2));//r2 starts with 1
  r2inds=&gel(split,2);
  int endoftheline=0;//Indicates if we can go on or not
  GEN inum=gen_0;
  long n;
  for(;;){//Each loop starts with non-empty r1inds and r2inds and we hit it with L, or endoftheline=1 and we backtrack to the last L and hit it with R
	if(!endoftheline) llist_putstart(&iseq,0);//Push onward L
	else{//Now we backtrack to the last zero
	  endoftheline=-1;
	  while(iseq!=NULL){
	    n=llist_pop(&iseq);
	    if(n==0){//Stop!
		  r1inds=&gel(ind1LR->data,2);
		  r2inds=&gel(ind2LR->data,2);
		  if(lg(*r1inds)>1 && lg(*r2inds)>1){endoftheline=0;break;}//No longer the end of the line.
	    }
	    glist_pop(&ind1LR);
	    glist_pop(&ind2LR);
	  }
	  if(endoftheline==-1) return gerepilecopy(top,inum);//We had to backtrack all the way, hence we are done.
	  llist_putstart(&iseq,1);//Now we add on a 1 and a 0
	  llist_putstart(&iseq,0);
	}//Update onto the next loop
	split=ibqf_intRS_splitindices(gel(r1,3),*r1inds);
	glist_putstart(&ind1LR,split);
	r1inds=&gel(split,1);
	split=ibqf_intRS_splitindices(gel(r2,3),*r2inds);
	glist_putstart(&ind2LR,split);
	r2inds=&gel(split,1);
	inum=addis(inum,(lg(gel(ind1LR->data,2))-1)*(lg(gel(ind2LR->data,1))-1));//Updating inum via ending with R for r1 and L for r2
	if(lg(*r1inds)==1 || lg(*r2inds)==1) endoftheline=1;//Must backtrack now
  }
}

//ibqf_intRS with typecheck
GEN ibqf_intRS_tc(GEN q1, GEN q2, long prec){
  pari_sp top=avma;
  GEN D1=bqf_checkdisc(q1);
  if(signe(D1)!=1) pari_err_TYPE("please supply an indefinite binary quadratic form",q1);
  GEN D2=bqf_checkdisc(q2);
  if(signe(D2)!=1) pari_err_TYPE("please supply an indefinite binary quadratic form",q2);
  GEN r1data=ibqf_river_positions(q1, gsqrt(D1,prec));
  GEN r2data=ibqf_river_positions(q2, gsqrt(D2,prec));
  return gerepileupto(top,ibqf_intRS_byriver(r1data, r2data));
}

//Supporting method to ibqf_intRS_byriver: Given the river river and a sorted t_VECSMALL of indices, this splits the set into those which are 0 and those which are 1, adds one to the indices, and returns them
static GEN ibqf_intRS_splitindices(GEN river, GEN ind){
  llist *L=NULL, *R=NULL;
  long len=lg(ind)-1, llen=0;
  if(ind[len]==lg(river)-1){//Sorting out if we pass over the end or not
	if(river[ind[len]]==0){llist_putstart(&L,1);llen++;}
	else llist_putstart(&R,1);
	len--;
  }
  for(long i=1;i<=len;i++){
	if(river[ind[i]]==0){llist_putstart(&L,ind[i]+1);llen++;}
	else llist_putstart(&R,ind[i]+1);
  }
  GEN rvec=cgetg(3,t_VEC);
  gel(rvec,1)=llist_tovecsmall(L,llen,-1);
  gel(rvec,2)=llist_tovecsmall(R,lg(ind)-1-llen,-1);
  return rvec;
}

//Finds the intersecting forms of all types, also outputs intersection data if data=1. If split=1, splits it into 4 vectors (RS, RO, LS, LO), else makes it one vector.
GEN ibqf_intforms_byriver(GEN r1, GEN r2, int data, int split){
  pari_sp top=avma;
  if(split){
    GEN rvec=cgetg(5,t_VEC);
	gel(rvec,1)=ibqf_intformsRS_byriver(r1,r2,data);
	gel(rvec,2)=ibqf_intformsRO_byriver(r1,r2,data);
	gel(rvec,3)=ibqf_intformsLS_byriver(r1,r2,data);
	gel(rvec,4)=ibqf_intformsLO_byriver(r1,r2,data);
	return rvec;
  }
  llist *indsRS1=NULL, *indsRS2=NULL, *indsRO1=NULL, *indsRO2=NULL, *indsLS1=NULL, *indsLS2=NULL, *indsLO1=NULL, *indsLO2=NULL;
  long inumRS=0, inumRO=0, inumLS=0, inumLO=0, r2len=lg(gel(r2,4))+1;
  GEN r3=ibqf_int_reverseriver(r2);//Reverse river of r2
  GEN rvec;
  if(!data){//NO DATA
    ibqf_intformsRS_byriver_indices(r1, r2, &indsRS1, &indsRS2, NULL, &inumRS, 0);//RS
    ibqf_intformsRS_byriver_indices(r3, r1, &indsRO2, &indsRO1, NULL, &inumRO, 0);//RO
    ibqf_intformsRS_byriver_indices(r2, r1, &indsLS2, &indsLS1, NULL, &inumLS, 0);//LS
    ibqf_intformsRS_byriver_indices(r1, r3, &indsLO1, &indsLO2, NULL, &inumLO, 0);//LO
    rvec=cgetg(inumRS+inumRO+inumLS+inumLO+1,t_VEC);
	
	long i=1, imax=inumRS;
    while(i<=imax){//RS
	  gel(rvec,i)=cgetg(3,t_VEC);
	  gel(gel(rvec,i),1)=ZV_copy(gel(gel(r1,4),indsRS1->data));
	  gel(gel(rvec,i),2)=ZV_copy(gel(gel(r2,4),indsRS2->data));
	  indsRS1=indsRS1->next;
	  indsRS2=indsRS2->next;
	  i++;
    }//RS end
	
	imax=imax+inumRO;
    while(i<=imax){//RO
	  gel(rvec,i)=cgetg(3,t_VEC);
	  gel(gel(rvec,i),1)=ZV_copy(gel(gel(r1,4),indsRO1->data));
	  if(indsRO2->data==1) gel(gel(rvec,i),2)=bqf_transS(gel(gel(r2,4),1));//1->1 and i->len(r2)+2-i for i!=1
	  else gel(gel(rvec,i),2)=bqf_transS(gel(gel(r2,4),r2len-(indsRO2->data)));
	  indsRO1=indsRO1->next;
	  indsRO2=indsRO2->next;
	  i++;
    }//RO end
	
	imax=imax+inumLS;
    while(i<=imax){//LS
	  gel(rvec,i)=cgetg(3,t_VEC);
	  gel(gel(rvec,i),1)=ZV_copy(gel(gel(r1,4),indsLS1->data));
	  gel(gel(rvec,i),2)=ZV_copy(gel(gel(r2,4),indsLS2->data));
	  indsLS1=indsLS1->next;
	  indsLS2=indsLS2->next;
	  i++;
    }//LS end
	
    imax=imax+inumLO;
    while(i<=imax){//LO
	  gel(rvec,i)=cgetg(3,t_VEC);
	  gel(gel(rvec,i),1)=ZV_copy(gel(gel(r1,4),indsLO1->data));
	  if(indsLO2->data==1) gel(gel(rvec,i),2)=bqf_transS(gel(gel(r2,4),1));//1->1 and i->len(r2)+2-i for i!=1
	  else gel(gel(rvec,i),2)=bqf_transS(gel(gel(r2,4),r2len-(indsLO2->data)));
	  indsLO1=indsLO1->next;
	  indsLO2=indsLO2->next;
	  i++;
    }//LO end
  }
  else{
    llist *loverlapRS=NULL, *loverlapRO=NULL, *loverlapLS=NULL, *loverlapLO=NULL;
    ibqf_intformsRS_byriver_indices(r1, r2, &indsRS1, &indsRS2, &loverlapRS, &inumRS, 1);//RS
    ibqf_intformsRS_byriver_indices(r3, r1, &indsRO2, &indsRO1, &loverlapRO, &inumRO, 1);//RO
    ibqf_intformsRS_byriver_indices(r2, r1, &indsLS2, &indsLS1, &loverlapLS, &inumLS, 1);//LS
    ibqf_intformsRS_byriver_indices(r1, r3, &indsLO1, &indsLO2, &loverlapLO, &inumLO, 1);//LO
    rvec=cgetg(inumRS+inumRO+inumLS+inumLO+1,t_VEC);
	
	long i=1, imax=inumRS;
    while(i<=imax){//RS
	  gel(rvec,i)=cgetg(6,t_VEC);
	  gel(gel(rvec,i),4)=ZV_copy(gel(gel(r1,4),indsRS1->data));
	  gel(gel(rvec,i),5)=ZV_copy(gel(gel(r2,4),indsRS2->data));
	  gel(gel(rvec,i),1)=bqf_bdelta(gel(gel(rvec,i),4),gel(gel(rvec,i),5));
	  gel(gel(rvec,i),2)=bqf_intlevel(gel(gel(rvec,i),4),gel(gel(rvec,i),5));
	  gel(gel(rvec,i),3)=stoi(loverlapRS->data);
	  indsRS1=indsRS1->next;
	  indsRS2=indsRS2->next;
	  loverlapRS=loverlapRS->next;
	  i++;
    }
	pari_free(loverlapRS);//RS end
	
	imax=imax+inumRO;
    while(i<=imax){//RO
	  gel(rvec,i)=cgetg(6,t_VEC);
	  gel(gel(rvec,i),4)=ZV_copy(gel(gel(r1,4),indsRO1->data));
	  if(indsRO2->data==1) gel(gel(rvec,i),5)=bqf_transS(gel(gel(r2,4),1));//1->1 and i->len(r2)+2-i for i!=1
	  else gel(gel(rvec,i),5)=bqf_transS(gel(gel(r2,4),r2len-(indsRO2->data)));
	  gel(gel(rvec,i),1)=bqf_bdelta(gel(gel(rvec,i),4),gel(gel(rvec,i),5));
	  gel(gel(rvec,i),2)=bqf_intlevel(gel(gel(rvec,i),4),gel(gel(rvec,i),5));
	  gel(gel(rvec,i),3)=stoi(loverlapRO->data);
	  indsRO1=indsRO1->next;
	  indsRO2=indsRO2->next;
	  loverlapRO=loverlapRO->next;
	  i++;
    }
	pari_free(loverlapRO);//RO end
	
	imax=imax+inumLS;
    while(i<=imax){//LS
	  gel(rvec,i)=cgetg(6,t_VEC);
	  gel(gel(rvec,i),4)=ZV_copy(gel(gel(r1,4),indsLS1->data));
	  gel(gel(rvec,i),5)=ZV_copy(gel(gel(r2,4),indsLS2->data));
	  gel(gel(rvec,i),1)=bqf_bdelta(gel(gel(rvec,i),4),gel(gel(rvec,i),5));
	  gel(gel(rvec,i),2)=bqf_intlevel(gel(gel(rvec,i),4),gel(gel(rvec,i),5));
	  gel(gel(rvec,i),3)=stoi(loverlapLS->data);
	  indsLS1=indsLS1->next;
	  indsLS2=indsLS2->next;
	  loverlapLS=loverlapLS->next;
	  i++;
    }
	pari_free(loverlapLS);//LS end
	
    imax=imax+inumLO;
    while(i<=imax){//LO
	  gel(rvec,i)=cgetg(6,t_VEC);
	  gel(gel(rvec,i),4)=ZV_copy(gel(gel(r1,4),indsLO1->data));
	  if(indsLO2->data==1) gel(gel(rvec,i),5)=bqf_transS(gel(gel(r2,4),1));//1->1 and i->len(r2)+2-i for i!=1
	  else gel(gel(rvec,i),5)=bqf_transS(gel(gel(r2,4),r2len-(indsLO2->data)));
	  gel(gel(rvec,i),1)=bqf_bdelta(gel(gel(rvec,i),4),gel(gel(rvec,i),5));
	  gel(gel(rvec,i),2)=bqf_intlevel(gel(gel(rvec,i),4),gel(gel(rvec,i),5));
	  gel(gel(rvec,i),3)=stoi(loverlapLO->data);
	  indsLO1=indsLO1->next;
	  indsLO2=indsLO2->next;
	  loverlapLO=loverlapLO->next;
	  i++;
    }
	pari_free(loverlapLO);//LO end
  }
  pari_free(indsRS1);
  pari_free(indsRS2);
  pari_free(indsRO1);
  pari_free(indsRO2);
  pari_free(indsLS1);
  pari_free(indsLS2);
  pari_free(indsLO1);
  pari_free(indsLO2);
  return gerepileupto(top,rvec);
}

//ibqf_intforms_byriver with typechecking and conversion to rivers
GEN ibqf_intforms_tc(GEN q1, GEN q2, int data, int split, long prec){
  pari_sp top=avma;
  GEN D1=bqf_checkdisc(q1);
  if(signe(D1)!=1) pari_err_TYPE("please supply an indefinite binary quadratic form",q1);
  GEN D2=bqf_checkdisc(q2);
  if(signe(D2)!=1) pari_err_TYPE("please supply an indefinite binary quadratic form",q2);
  GEN r1data=ibqf_river_positions_forms(q1, gsqrt(D1,prec));
  GEN r2data=ibqf_river_positions_forms(q2, gsqrt(D2,prec));
  return gerepileupto(top,ibqf_intforms_byriver(r1data, r2data, data, split));
}

//ibqf_intformsRS with the river inputted. If data=1, each output is [B_Delta(f1,f2), signed level, river overlap length, f1, f2].
GEN ibqf_intformsRS_byriver(GEN r1, GEN r2, int data){
  pari_sp top=avma;
  llist *inds1=NULL, *inds2=NULL;
  long inum=0;
  GEN rvec;
  if(data==1){
	llist *loverlap=NULL;
    ibqf_intformsRS_byriver_indices(r1, r2, &inds1, &inds2, &loverlap, &inum, 1);
    rvec=cgetg(inum+1,t_VEC);
    for(long i=1;i<=inum;i++){
	  gel(rvec,i)=cgetg(6,t_VEC);
	  gel(gel(rvec,i),4)=ZV_copy(gel(gel(r1,4),inds1->data));
	  gel(gel(rvec,i),5)=ZV_copy(gel(gel(r2,4),inds2->data));
	  gel(gel(rvec,i),1)=bqf_bdelta(gel(gel(rvec,i),4),gel(gel(rvec,i),5));
	  gel(gel(rvec,i),2)=bqf_intlevel(gel(gel(rvec,i),4),gel(gel(rvec,i),5));
	  gel(gel(rvec,i),3)=stoi(loverlap->data);
	  inds1=inds1->next;
	  inds2=inds2->next;
	  loverlap=loverlap->next;
    }
	pari_free(loverlap);
  }
  else{
	ibqf_intformsRS_byriver_indices(r1, r2, &inds1, &inds2, NULL, &inum, 0);//No data
    rvec=cgetg(inum+1,t_VEC);
    for(long i=1;i<=inum;i++){
	  gel(rvec,i)=cgetg(3,t_VEC);
	  gel(gel(rvec,i),1)=ZV_copy(gel(gel(r1,4),inds1->data));
	  gel(gel(rvec,i),2)=ZV_copy(gel(gel(r2,4),inds2->data));
	  inds1=inds1->next;
	  inds2=inds2->next;
    } 
  }
  pari_free(inds1);
  pari_free(inds2);
  return gerepileupto(top,rvec);
}

//ibqf_intformsRS_byriver with typechecking
GEN ibqf_intformsRS_tc(GEN q1, GEN q2, int data, long prec){
  pari_sp top=avma;
  GEN D1=bqf_checkdisc(q1);
  if(signe(D1)!=1) pari_err_TYPE("please supply an indefinite binary quadratic form",q1);
  GEN D2=bqf_checkdisc(q2);
  if(signe(D2)!=1) pari_err_TYPE("please supply an indefinite binary quadratic form",q2);
  GEN r1data=ibqf_river_positions_forms(q1, gsqrt(D1,prec));
  GEN r2data=ibqf_river_positions_forms(q2, gsqrt(D2,prec));
  return gerepileupto(top,ibqf_intformsRS_byriver(r1data, r2data, data));
}

//ibqf_intformsRO with the river inputted. If data=1, each output is [B_Delta(f1,f2), signed level, river overlap length, f1, f2].
//We use Int_RO(f,g)=Int_RS(-g,f)
GEN ibqf_intformsRO_byriver(GEN r1, GEN r2, int data){
  pari_sp top=avma;
  llist *inds1=NULL, *inds2=NULL;
  long inum=0, r2len=lg(gel(r2,4))+1;
  GEN r3=ibqf_int_reverseriver(r2);//Reverse river of r2
  GEN rvec;
  if(data==1){
	llist *loverlap=NULL;
    ibqf_intformsRS_byriver_indices(r3, r1, &inds2, &inds1, &loverlap, &inum, 1);
    rvec=cgetg(inum+1,t_VEC);
    for(long i=1;i<=inum;i++){
	  gel(rvec,i)=cgetg(6,t_VEC);
	  gel(gel(rvec,i),4)=ZV_copy(gel(gel(r1,4),inds1->data));
	  if(inds2->data==1) gel(gel(rvec,i),5)=bqf_transS(gel(gel(r2,4),1));//1->1 and i->len(r2)+2-i for i!=1
	  else gel(gel(rvec,i),5)=bqf_transS(gel(gel(r2,4),r2len-(inds2->data)));
	  gel(gel(rvec,i),1)=bqf_bdelta(gel(gel(rvec,i),4),gel(gel(rvec,i),5));
	  gel(gel(rvec,i),2)=bqf_intlevel(gel(gel(rvec,i),4),gel(gel(rvec,i),5));
	  gel(gel(rvec,i),3)=stoi(loverlap->data);
	  inds1=inds1->next;
	  inds2=inds2->next;
	  loverlap=loverlap->next;
    }
	pari_free(loverlap);
  }
  else{
	ibqf_intformsRS_byriver_indices(r3, r1, &inds2, &inds1, NULL, &inum, 0);//No data
    rvec=cgetg(inum+1,t_VEC);
    for(long i=1;i<=inum;i++){
	  gel(rvec,i)=cgetg(3,t_VEC);
	  gel(gel(rvec,i),1)=ZV_copy(gel(gel(r1,4),inds1->data));
	  if(inds2->data==1) gel(gel(rvec,i),2)=bqf_transS(gel(gel(r2,4),1));//1->1 and i->len(r2)+2-i for i!=1
	  else gel(gel(rvec,i),2)=bqf_transS(gel(gel(r2,4),r2len-(inds2->data)));
	  inds1=inds1->next;
	  inds2=inds2->next;
    } 
  }
  pari_free(inds1);
  pari_free(inds2);
  return gerepileupto(top,rvec);
}

//ibqf_intformsRO_byriver with typechecking
GEN ibqf_intformsRO_tc(GEN q1, GEN q2, int data, long prec){
  pari_sp top=avma;
  GEN D1=bqf_checkdisc(q1);
  if(signe(D1)!=1) pari_err_TYPE("please supply an indefinite binary quadratic form",q1);
  GEN D2=bqf_checkdisc(q2);
  if(signe(D2)!=1) pari_err_TYPE("please supply an indefinite binary quadratic form",q2);
  GEN r1data=ibqf_river_positions_forms(q1, gsqrt(D1,prec));
  GEN r2data=ibqf_river_positions_forms(q2, gsqrt(D2,prec));
  return gerepileupto(top,ibqf_intformsRO_byriver(r1data, r2data, data));
}

//ibqf_intformsLS with the river inputted. If data=1, each output is [B_Delta(f1,f2), signed level, river overlap length, f1, f2].
//We use Int_LS(f,g)=Int_RS(g,f)
GEN ibqf_intformsLS_byriver(GEN r1, GEN r2, int data){
  pari_sp top=avma;
  llist *inds1=NULL, *inds2=NULL;
  long inum=0;
  GEN rvec;
  if(data==1){
	llist *loverlap=NULL;
    ibqf_intformsRS_byriver_indices(r2, r1, &inds2, &inds1, &loverlap, &inum, 1);
    rvec=cgetg(inum+1,t_VEC);
    for(long i=1;i<=inum;i++){
	  gel(rvec,i)=cgetg(6,t_VEC);
	  gel(gel(rvec,i),4)=ZV_copy(gel(gel(r1,4),inds1->data));
	  gel(gel(rvec,i),5)=ZV_copy(gel(gel(r2,4),inds2->data));
	  gel(gel(rvec,i),1)=bqf_bdelta(gel(gel(rvec,i),4),gel(gel(rvec,i),5));
	  gel(gel(rvec,i),2)=bqf_intlevel(gel(gel(rvec,i),4),gel(gel(rvec,i),5));
	  gel(gel(rvec,i),3)=stoi(loverlap->data);
	  inds1=inds1->next;
	  inds2=inds2->next;
	  loverlap=loverlap->next;
    }
	pari_free(loverlap);
  }
  else{
	ibqf_intformsRS_byriver_indices(r2, r1, &inds2, &inds1, NULL, &inum, 0);//No data
    rvec=cgetg(inum+1,t_VEC);
    for(long i=1;i<=inum;i++){
	  gel(rvec,i)=cgetg(3,t_VEC);
	  gel(gel(rvec,i),1)=ZV_copy(gel(gel(r1,4),inds1->data));
	  gel(gel(rvec,i),2)=ZV_copy(gel(gel(r2,4),inds2->data));
	  inds1=inds1->next;
	  inds2=inds2->next;
    } 
  }
  pari_free(inds1);
  pari_free(inds2);
  return gerepileupto(top,rvec);
}

//ibqf_intformsLS_byriver with typechecking
GEN ibqf_intformsLS_tc(GEN q1, GEN q2, int data, long prec){
  pari_sp top=avma;
  GEN D1=bqf_checkdisc(q1);
  if(signe(D1)!=1) pari_err_TYPE("please supply an indefinite binary quadratic form",q1);
  GEN D2=bqf_checkdisc(q2);
  if(signe(D2)!=1) pari_err_TYPE("please supply an indefinite binary quadratic form",q2);
  GEN r1data=ibqf_river_positions_forms(q1, gsqrt(D1,prec));
  GEN r2data=ibqf_river_positions_forms(q2, gsqrt(D2,prec));
  return gerepileupto(top,ibqf_intformsLS_byriver(r1data, r2data, data));
}

//ibqf_intformsRS with the river inputted. If data=1, each output is [B_Delta(f1,f2), signed level, river overlap length, f1, f2].
//We use Int_LO(f,g)=Int_RS(f,-g)
GEN ibqf_intformsLO_byriver(GEN r1, GEN r2, int data){
  pari_sp top=avma;
  llist *inds1=NULL, *inds2=NULL;
  long inum=0, r2len=lg(gel(r2,4))+1;
  GEN r3=ibqf_int_reverseriver(r2);//Reverse river of r2
  GEN rvec;
  if(data==1){
	llist *loverlap=NULL;
    ibqf_intformsRS_byriver_indices(r1, r3, &inds1, &inds2, &loverlap, &inum, 1);
    rvec=cgetg(inum+1,t_VEC);
    for(long i=1;i<=inum;i++){
	  gel(rvec,i)=cgetg(6,t_VEC);
	  gel(gel(rvec,i),4)=ZV_copy(gel(gel(r1,4),inds1->data));
	  if(inds2->data==1) gel(gel(rvec,i),5)=bqf_transS(gel(gel(r2,4),1));//1->1 and i->len(r2)+2-i for i!=1
	  else gel(gel(rvec,i),5)=bqf_transS(gel(gel(r2,4),r2len-(inds2->data)));
	  gel(gel(rvec,i),1)=bqf_bdelta(gel(gel(rvec,i),4),gel(gel(rvec,i),5));
	  gel(gel(rvec,i),2)=bqf_intlevel(gel(gel(rvec,i),4),gel(gel(rvec,i),5));
	  gel(gel(rvec,i),3)=stoi(loverlap->data);
	  inds1=inds1->next;
	  inds2=inds2->next;
	  loverlap=loverlap->next;
    }
	pari_free(loverlap);
  }
  else{
	ibqf_intformsRS_byriver_indices(r1, r3, &inds1, &inds2, NULL, &inum, 0);//No data
    rvec=cgetg(inum+1,t_VEC);
    for(long i=1;i<=inum;i++){
	  gel(rvec,i)=cgetg(3,t_VEC);
	  gel(gel(rvec,i),1)=ZV_copy(gel(gel(r1,4),inds1->data));
	  if(inds2->data==1) gel(gel(rvec,i),2)=bqf_transS(gel(gel(r2,4),1));//1->1 and i->len(r2)+2-i for i!=1
	  else gel(gel(rvec,i),2)=bqf_transS(gel(gel(r2,4),r2len-(inds2->data)));
	  inds1=inds1->next;
	  inds2=inds2->next;
    } 
  }
  pari_free(inds1);
  pari_free(inds2);
  return gerepileupto(top,rvec);
}

//ibqf_intformsRS_byriver with typechecking
GEN ibqf_intformsLO_tc(GEN q1, GEN q2, int data, long prec){
  pari_sp top=avma;
  GEN D1=bqf_checkdisc(q1);
  if(signe(D1)!=1) pari_err_TYPE("please supply an indefinite binary quadratic form",q1);
  GEN D2=bqf_checkdisc(q2);
  if(signe(D2)!=1) pari_err_TYPE("please supply an indefinite binary quadratic form",q2);
  GEN r1data=ibqf_river_positions_forms(q1, gsqrt(D1,prec));
  GEN r2data=ibqf_river_positions_forms(q2, gsqrt(D2,prec));
  return gerepileupto(top,ibqf_intformsLO_byriver(r1data, r2data, data));
}

//ibqf_intformsRS_byriver_indices where we just update the index list and river overlap length list
static void ibqf_intformsRS_byriver_indices(GEN r1, GEN r2, llist **inds1, llist **inds2, llist **loverlap, long *inum, int data){
  pari_sp top=avma;
  GEN split;//For splitting the indices
  llist *iseq=NULL;//Stores the overlap in sequences
  long overlap=0;//length of overlap
  glist *ind1LR=NULL;//Stores the indices of r1 which give 0/1 and preceeded by iseq
  glist *ind2LR=NULL;//Stores the indices of r2 which give 0/1 and preceeded by iseq
  //glist *inumind=NULL;//Stores the pairs of indices
  GEN *r1inds, *r2inds;//Pointers to the current set of indices we are looking at
  split=ibqf_intRS_splitindices(gel(r1,3),gel(r1,1));//r1 starts with 0
  r1inds=&gel(split,1);
  split=ibqf_intRS_splitindices(gel(r2,3),gel(r2,2));//r2 starts with 1
  r2inds=&gel(split,2);
  int endoftheline=0;//Indicates if we can go on or not
  long n, l1, l2, istart, jstart, toadd, r1len=lg(gel(r1,3))-1, r2len=lg(gel(r2,3))-1;
  for(;;){//Each loop starts with non-empty r1inds and r2inds and we hit it with L, or endoftheline=1 and we backtrack to the last L and hit it with R
	if(!endoftheline){llist_putstart(&iseq,0);overlap++;}//Push onward L
	else{//Now we backtrack to the last zero
	  endoftheline=-1;
	  while(iseq!=NULL){
	    n=llist_pop(&iseq);
		overlap--;
	    if(n==0){//Stop!
		  r1inds=&gel(ind1LR->data,2);
		  r2inds=&gel(ind2LR->data,2);
		  if(lg(*r1inds)>1 && lg(*r2inds)>1){endoftheline=0;break;}//No longer the end of the line.
	    }
		cgiv(glist_pop(&ind2LR));//These are the last GENS on the stack, and we can give them back.
	    cgiv(glist_pop(&ind1LR));
	  }
	  if(endoftheline==-1){avma=top;return;}//We had to backtrack all the way, hence we are done.
	  llist_putstart(&iseq,1);//Now we add on a 1 and a 0
	  llist_putstart(&iseq,0);
	  overlap=overlap+2;
	}//Update onto the next loop
	split=ibqf_intRS_splitindices(gel(r1,3),*r1inds);
	glist_putstart(&ind1LR,split);
	r1inds=&gel(split,1);
	split=ibqf_intRS_splitindices(gel(r2,3),*r2inds);
	glist_putstart(&ind2LR,split);
	r2inds=&gel(split,1);
	l1=lg(gel(ind1LR->data,2))-1;
	l2=lg(gel(ind2LR->data,1))-1;
	toadd=l1*l2;
	*inum=*inum+toadd;//Updating inum via ending with R for r1 and L for r2
	if(data){for(long i=1;i<=toadd;i++) llist_putstart(loverlap,overlap);}//Updating overlap if needed
	if(l1>0 && gel(ind1LR->data,2)[1]==1){//Have to shift indices back by one
	  istart=2;
	  if(l2>0 && gel(ind2LR->data,1)[1]==1){
		  jstart=2;
		  llist_putstart(inds1, r1len);llist_putstart(inds2, r2len);
		  for(long j=2;j<=l2;j++){llist_putstart(inds1, r1len);llist_putstart(inds2, gel(ind2LR->data,1)[j]-1);}
		  for(long i=2;i<=l1;i++){llist_putstart(inds1, gel(ind1LR->data,2)[i]-1);llist_putstart(inds2, r2len);}
	  }
	  else{
		jstart=1;
		for(long j=1;j<=l2;j++){llist_putstart(inds1, r1len);llist_putstart(inds2, gel(ind2LR->data,1)[j]-1);}
	  }
	}
	else{
	  istart=1;
	  if(l2>0 && gel(ind2LR->data,1)[1]==1){
		  jstart=2;
		  for(long i=1;i<=l1;i++){llist_putstart(inds1, gel(ind1LR->data,2)[i]-1);llist_putstart(inds2, r2len);}
	  }
	  else jstart=1;
	}
	for(long i=istart;i<=l1;i++){
	  for(long j=jstart;j<=l2;j++){//Updating the index list
	    llist_putstart(inds1, gel(ind1LR->data,2)[i]-1);
		llist_putstart(inds2, gel(ind2LR->data,1)[j]-1);
	  }
	}
	if(lg(*r1inds)==1 || lg(*r2inds)==1) endoftheline=1;//Must backtrack now
  }
}

