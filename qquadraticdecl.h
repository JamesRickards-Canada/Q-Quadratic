//c_base.c

//STRUCTURES

typedef struct listtype1{//A circular list of GENs, stores data, next, and previous terms
  GEN data; 
  struct listtype1 *next;
  struct listtype1 *prev;
}clist;

typedef struct listtype2{//A generic linked list of GENs, stores data and next term
  GEN data; 
  struct listtype2 *next;
}glist;

typedef struct listtype3{//A generic linked list of longs, stores data and next term
  long data; 
  struct listtype3 *next;
}llist;

//METHODS

//COMPLEX GEOMETRY
GEN crossratio(GEN a, GEN b, GEN c, GEN d);
GEN mat_eval(GEN M, GEN x);
GEN mat_eval_typecheck(GEN M, GEN x);

//INFINITY 
GEN addoo(GEN a, GEN b);
GEN divoo(GEN a, GEN b);

//LINEAR EQUATIONS AND MATRICES
GEN lin_intsolve(GEN A, GEN B, GEN n);
GEN lin_intsolve_typecheck(GEN A, GEN B, GEN n);
GEN mat3_complete(GEN A, GEN B, GEN C);
GEN mat3_complete_typecheck(GEN A, GEN B, GEN C);

//SOLVING EQUATIONS MOD N
GEN sqmod(GEN x, GEN n, GEN fact);
GEN sqmod_typecheck(GEN x, GEN n);

//INTEGER VECTORS
GEN ZV_copy(GEN v);
int ZV_equal(GEN v1, GEN v2);
GEN ZV_Z_divexact(GEN v, GEN y);
GEN ZV_Z_mul(GEN v, GEN x);

//TIME
char *returntime(void);
void printtime(void);

//LISTS
void clist_free(clist *l, long length);
void clist_putbefore(clist **head_ref, GEN new_data);
void clist_putafter(clist **head_ref, GEN new_data);
GEN clist_togvec(clist *l, long length, int dir);
void glist_free(glist *l);
GEN glist_pop(glist **head_ref);
void glist_putstart(glist **head_ref, GEN new_data);
GEN glist_togvec(glist *l, long length, int dir);
void llist_free(llist *l);
long llist_pop(llist **head_ref);
void llist_putstart(llist **head_ref, long new_data);
GEN llist_togvec(llist *l, long length, int dir);
GEN llist_tovecsmall(llist *l, long length, int dir);


//c_bqf.c

//METHODS

//DISCRIMINANT METHODS
GEN disclist(GEN D1, GEN D2, int fund, GEN cop);
GEN discprimeindex(GEN D, GEN facs);
GEN discprimeindex_typecheck(GEN D);
GEN fdisc(GEN D);
GEN fdisc_typecheck(GEN D);
int isdisc(GEN D);
GEN pell(GEN D);
GEN pell_typecheck(GEN D);
GEN posreg(GEN D, long prec);
GEN posreg_typecheck(GEN D, long prec);
GEN quadroot(GEN D);
GEN quadroot_typecheck(GEN D);

//BASIC OPERATIONS ON BINARY QUADRATIC FORMS
GEN bqf_automorph_typecheck(GEN q);
int bqf_compare(void *data, GEN q1, GEN q2);
int bqf_compare_tmat(void *data, GEN d1, GEN d2);
GEN bqf_disc(GEN q);
GEN bqf_disc_typecheck(GEN q);
int bqf_isreduced(GEN q, int Dsign);
int bqf_isreduced_typecheck(GEN q);
GEN bqf_random(GEN maxc, int type, int primitive);
GEN bqf_random_D(GEN maxc, GEN D);
GEN bqf_red(GEN q, GEN rootD, int Dsign, int tmat);
GEN bqf_red_typecheck(GEN q, int tmat, long prec);
GEN bqf_roots(GEN q, GEN D, GEN w);
GEN bqf_roots_typecheck(GEN q);
GEN bqf_trans(GEN q, GEN mtx);
GEN bqf_trans_typecheck(GEN q, GEN mtx);
GEN bqf_transL(GEN q, GEN n);
GEN bqf_transR(GEN q, GEN n);
GEN bqf_transS(GEN q);
GEN bqf_trans_coprime(GEN q, GEN n);
GEN bqf_trans_coprime_typecheck(GEN q, GEN n);

//BASIC METHODS FOR NEGATIVE DISCRIMINANTS
GEN dbqf_automorph(GEN q, GEN D);
GEN dbqf_red(GEN q);
GEN dbqf_red_tmat(GEN q);

//BASIC OPERATIONS SPECIFIC TO INDEFINITE FORMS/POSITIVE DISCRIMINANTS
GEN ibqf_automorph_D(GEN q, GEN D);
GEN ibqf_automorph_pell(GEN q, GEN qpell);
int ibqf_isrecip(GEN q, GEN rootD);
int ibqf_isrecip_typecheck(GEN q, long prec);
GEN ibqf_leftnbr(GEN q, GEN rootD);
GEN ibqf_leftnbr_tmat(GEN q, GEN rootD);
GEN ibqf_leftnbr_typecheck(GEN q, int tmat, long prec);
GEN ibqf_leftnbr_update(GEN qvec, GEN rootD);
GEN ibqf_red(GEN q, GEN rootD);
GEN ibqf_red_tmat(GEN q, GEN rootD);
GEN ibqf_red_pos(GEN q, GEN rootD);
GEN ibqf_red_pos_tmat(GEN q, GEN rootD);
GEN ibqf_redorbit(GEN q, GEN rootD);
GEN ibqf_redorbit_tmat(GEN q, GEN rootD);
GEN ibqf_redorbit_posonly(GEN q, GEN rootD);
GEN ibqf_redorbit_posonly_tmat(GEN q, GEN rootD);
GEN ibqf_redorbit_typecheck(GEN q, int tmat, int posonly, long prec);
GEN ibqf_rightnbr(GEN q, GEN rootD);
GEN ibqf_rightnbr_tmat(GEN q, GEN rootD);
GEN ibqf_rightnbr_update(GEN qvec, GEN rootD);
GEN ibqf_rightnbr_typecheck(GEN q, int tmat, long prec);
GEN ibqf_river(GEN q, GEN rootD);
GEN ibqf_river_positions(GEN q, GEN rootD);
GEN ibqf_river_positions_forms(GEN q, GEN rootD);
GEN ibqf_river_typecheck(GEN q, long prec);
GEN ibqf_riverforms(GEN q, GEN rootD);
GEN ibqf_riverforms_typecheck(GEN q, long prec);
GEN ibqf_symmetricarc(GEN q, GEN D, GEN rootD, GEN qpell, long prec);
GEN ibqf_symmetricarc_typecheck(GEN q, long prec);
GEN ibqf_toriver(GEN q, GEN rootD);
GEN ibqf_toriver_tmat(GEN q, GEN rootD);
GEN mat_toibqf(GEN mtx);
GEN mat_toibqf_typecheck(GEN mtx);

//EQUIVALENCE OF BQFs
GEN bqf_isequiv(GEN q1, GEN q2, GEN rootD, int Dsign, int tmat);
GEN bqf_isequiv_set(GEN q, GEN S, GEN rootD, int Dsign, int tmat);
GEN bqf_isequiv_typecheck(GEN q1, GEN q2, int tmat, long prec);
GEN dbqf_isequiv(GEN q1, GEN q2);
GEN dbqf_isequiv_tmat(GEN q1, GEN q2);
long dbqf_isequiv_set(GEN q, GEN S);
GEN dbqf_isequiv_set_tmat(GEN q, GEN S);
GEN ibqf_isequiv(GEN q1, GEN q2, GEN rootD);
GEN ibqf_isequiv_tmat(GEN q1, GEN q2, GEN rootD);
long ibqf_isequiv_set_byq(GEN q, GEN S, GEN rootD);
long ibqf_isequiv_set_byq_presorted(GEN qredsorted, GEN S, GEN rootD);
GEN ibqf_isequiv_set_byq_tmat(GEN q, GEN S, GEN rootD);
GEN ibqf_isequiv_set_byq_tmat_presorted(GEN qredsorted, GEN S, GEN rootD);
long ibqf_isequiv_set_byS(GEN q, GEN S, GEN rootD);
long ibqf_isequiv_set_byS_presorted(GEN q, GEN Sreds, GEN perm, GEN rootD);
GEN ibqf_isequiv_set_byS_tmat(GEN q, GEN S, GEN rootD);
GEN ibqf_isequiv_set_byS_tmat_presorted(GEN q, GEN Sreds, GEN perm, GEN rootD);


//CLASS GROUPS AND COMPOSITION OF FORMS
GEN bqf_comp(GEN q1, GEN q2);
GEN bqf_comp_red(GEN q1, GEN q2, GEN rootD, int Dsign);
GEN bqf_comp_typecheck(GEN q1, GEN q2, int tored, long prec);
GEN bqf_idelt(GEN D);
GEN bqf_ncgp(GEN D, long prec);
GEN bqf_ncgp_lexic(GEN D, long prec);
GEN bqf_pow(GEN q, GEN n);
GEN bqf_pow_red(GEN q, GEN n, GEN rootD, int Dsign);
GEN bqf_pow_typecheck(GEN q, GEN n, int tored, long prec);
GEN bqf_square(GEN q);
GEN bqf_square_red(GEN q, GEN rootD, int Dsign);
GEN bqf_square_typecheck(GEN q, int tored, long prec);
GEN ideal_tobqf(GEN numf, GEN ideal);

//REPRESENTATION OF NUMBERS BY BQFs
GEN bqf_reps(GEN q, GEN n, int proper, int half, long prec);
GEN bqf_reps_typecheck(GEN q, GEN n, int proper, int half, long prec);
GEN dbqf_reps(GEN qred, GEN D, GEN n, int proper, int half);
GEN ibqf_reps(GEN qorb, GEN qautom, GEN D, GEN rootD, GEN n, int proper, int half);
GEN sbqf_reps(GEN q, GEN D, GEN rootD, GEN n, int half);
GEN zbqf_reps(GEN A, GEN B, GEN n, int half);

//MORE REPRESENTATION OF NUMBERS
GEN bqf_bigreps(GEN q, GEN n, long prec);
GEN bqf_bigreps_typecheck(GEN q, GEN n, long prec);
GEN bqf_linearsolve(GEN q, GEN n1, GEN lin, GEN n2, long prec);
GEN bqf_linearsolve_typecheck(GEN q, GEN n1, GEN lin, GEN n2, long prec);

//GENERAL CHECKING METHODS
void bqf_check(GEN q);
GEN bqf_checkdisc(GEN q);
void intmatrix_check(GEN mtx);
