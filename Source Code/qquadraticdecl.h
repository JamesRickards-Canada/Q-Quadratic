//qq_base.c

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

//RANDOM
GEN rand_elt(GEN v);
long rand_l(long len);

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


//qq_bqf.c

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


//qq_bqf_int.c

//METHODS

//INTERSECTION DATA
GEN bqf_bdelta(GEN q1, GEN q2);
GEN bqf_bdelta_typecheck(GEN q1, GEN q2);
GEN bqf_intlevel(GEN q1, GEN q2);
GEN bqf_intlevel_typecheck(GEN q1, GEN q2);
GEN ibqf_intpairs_transtoq(GEN pairs, GEN q, GEN rootD);
GEN ibqf_intpoint(GEN q1, GEN q2, GEN location, GEN autom);
GEN ibqf_intpoint_typecheck(GEN q1, GEN q2, GEN location);
GEN hdist(GEN z1, GEN z2, long prec);
GEN hdist_typecheck(GEN z1, GEN z2, long prec);

//INTERSECTION NUMBER COMPUTATION
GEN ibqf_int(GEN r1, GEN r2);
GEN ibqf_int_typecheck(GEN q1, GEN q2, long prec);
GEN ibqf_intRS_byriver(GEN r1, GEN r2);
GEN ibqf_intRS_typecheck(GEN q1, GEN q2, long prec);
GEN ibqf_intforms_byriver(GEN r1, GEN r2, int data, int split);
GEN ibqf_intforms_typecheck(GEN q1, GEN q2, int data, int split, long prec);
GEN ibqf_intformsRS_byriver(GEN r1, GEN r2, int data);
GEN ibqf_intformsRS_typecheck(GEN q1, GEN q2, int data, long prec);
GEN ibqf_intformsRO_byriver(GEN r1, GEN r2, int data);
GEN ibqf_intformsRO_typecheck(GEN q1, GEN q2, int data, long prec);
GEN ibqf_intformsLS_byriver(GEN r1, GEN r2, int data);
GEN ibqf_intformsLS_typecheck(GEN q1, GEN q2, int data, long prec);
GEN ibqf_intformsLO_byriver(GEN r1, GEN r2, int data);
GEN ibqf_intformsLO_typecheck(GEN q1, GEN q2, int data, long prec);


//qq_bqf_int_scripts.c

//METHODS

//OLD INTERSECTION ALGORITHMS
GEN ibqf_intRS_byriverOLDBAD(GEN r1, GEN r2);
GEN ibqf_intRSOLD(GEN q1, GEN q2, long prec);

//INUM DATA
GEN ibqf_int_dists(GEN q1, GEN q2, GEN rootD1, GEN q1auto, GEN z, int ispairs, long prec);
GEN ibqf_int_dists_typecheck(GEN q1, GEN q2, long prec);
GEN ibqf_intangle(GEN q1, GEN q2, GEN D1D2, long prec);
GEN ibqf_int_angles(GEN q1, GEN q2, GEN D1, GEN D2, int ispairs, long prec);
GEN ibqf_intsets(GEN S1, GEN S2, long prec);

//COMPILING DATA
void ibqf_int_timings(GEN D1range, GEN D2range, GEN filename, long trials, long prec);
GEN ibqf_int_dists_histogram(GEN q1, GEN D2, char *filename, char *autocompile, long prec);
GEN ibqf_int_angles_histogram(GEN q1, GEN D2, char *filename, char *autocompile, long prec);
GEN ibqf_CD1D2data(GEN D1range, GEN D2range, long trials, int scaled, char *filename, char *autocompile, long prec);

//HISTOGRAMS
//void hist_make(GEN sorteddata, GEN minx, GEN maxx, GEN shift, int scaled, GEN filename, GEN autocompile, long prec);
//void hist_rebin(GEN data, GEN nbins, int scale, long prec);

//PLOTTING INTERSECTION NUMBERS
GEN ibqf_plotint_arc(GEN q1, GEN q2, GEN filename, long prec);
GEN ibqf_plotint_bydist(GEN q1, GEN q2, GEN filename, long prec);
GEN ibqf_plotint_fd(GEN q1, GEN q2, GEN filename, long prec);


//qq_qquat.c

//METHODS

//BASIC OPERATIONS ON ELEMENTS IN QUATERNION ALGEBRAS
GEN qa_conj(GEN x);
GEN qa_conj_typecheck(GEN x);
GEN qa_conjby(GEN Q, GEN x, GEN y);
GEN qa_conjby_typecheck(GEN Q, GEN x, GEN y);
GEN qa_inv(GEN Q, GEN x);
GEN qa_inv_typecheck(GEN Q, GEN x);
GEN qa_m2rembed(GEN Q, GEN x);
GEN qa_m2rembed_typecheck(GEN Q, GEN x);
GEN qa_minpoly(GEN Q, GEN x);
GEN qa_minpoly_typecheck(GEN Q, GEN x);
GEN qa_mul(GEN Q, GEN x, GEN y);
GEN qa_mul_typecheck(GEN Q, GEN x, GEN y);
GEN qa_norm(GEN Q, GEN x);
GEN qa_norm_typecheck(GEN Q, GEN x);
GEN qa_pow(GEN Q, GEN x, GEN n);
GEN qa_pow_typecheck(GEN Q, GEN x, GEN n);
GEN qa_roots(GEN Q, GEN x, long prec);
GEN qa_roots_typecheck(GEN Q, GEN x, long prec);
GEN qa_square(GEN Q, GEN x);
GEN qa_square_typecheck(GEN Q, GEN x);
GEN qa_trace(GEN x);
GEN qa_trace_typecheck(GEN x);

//BASIC OPERATIONS ON ORDERS/LATTICES IN QUATERNION ALGEBRAS
int qa_isinorder(GEN Q, GEN ordinv, GEN x);
int qa_isinorder_typecheck(GEN Q, GEN ord, GEN x);
GEN qa_ord_conj(GEN Q, GEN ord, GEN c);
GEN qa_ord_conj_typecheck(GEN Q, GEN ord, GEN c);
GEN qa_ord_disc(GEN Q, GEN ord);
GEN qa_ord_disc_typecheck(GEN Q, GEN ord);

//INITIALIZATION METHODS
GEN qa_ord_init(GEN Q, GEN ord);
GEN qa_ord_init_typecheck(GEN Q, GEN ord);
GEN qa_init_primes(GEN pset, int type);
GEN qa_init_primes_typecheck(GEN pset);
GEN qa_init_2primes(GEN p, GEN q);
GEN qa_init_2primes_typecheck(GEN p, GEN q);
GEN qa_ram_fromab(GEN a, GEN b);
GEN qa_ram_fromab_typecheck(GEN a, GEN b);

//CONJUGATION OF ELEMENTS IN A GIVEN ORDER
GEN qa_conjbasis(GEN Q, GEN ord, GEN ordinv, GEN e1, GEN e2, int orient);
GEN qa_conjbasis_typecheck(GEN Q, GEN ord, GEN e1, GEN e2, int orient);
GEN qa_conjqf(GEN Q, GEN ord, GEN ordinv, GEN e1, GEN e2);
GEN qa_conjqf_typecheck(GEN Q, GEN ord, GEN e1, GEN e2);
GEN qa_conjnorm(GEN Q, GEN ord, GEN ordinv, GEN e1, GEN e2, GEN n, int retconelt, long prec);
GEN qa_conjnorm_typecheck(GEN Q, GEN ord, GEN e1, GEN e2, GEN n, int retconelt, long prec);
GEN qa_simulconj(GEN Q, GEN ord, GEN ordinv, GEN e1, GEN e2, GEN f1, GEN f2, long prec);
GEN qa_simulconj_typecheck(GEN Q, GEN ord, GEN e1, GEN e2, GEN f1, GEN f2, long prec);

//EMBEDDING QUADRATIC ORDERS INTO EICHLER ORDERS
GEN qa_associatedemb(GEN Q, GEN order, GEN emb, GEN D);
GEN qa_associatedemb_typecheck(GEN Q, GEN order, GEN emb, GEN D);
GEN qa_embed(GEN Q, GEN order, GEN D, GEN nembeds, GEN rpell, long prec);
GEN qa_embed_typecheck(GEN Q, GEN order, GEN D, GEN nembeds, int rpell, long prec);
GEN qa_embeddablediscs(GEN Q, GEN order, GEN d1, GEN d2, int fund, GEN cop);
GEN qa_embeddablediscs_typecheck(GEN Q, GEN order, GEN d1, GEN d2, int fund, GEN cop);
GEN qa_numemb(GEN Q, GEN order, GEN D, GEN narclno);
GEN qa_numemb_typecheck(GEN Q, GEN order, GEN D, GEN narclno, long prec);
GEN qa_ordiffer(GEN Q, GEN order, GEN e1, GEN e2, GEN D);
GEN qa_ordiffer_typecheck(GEN Q, GEN order, GEN e1, GEN e2, GEN D);
GEN qa_orinfinite(GEN Q, GEN emb, GEN D, long prec);
GEN qa_orinfinite_typecheck(GEN Q, GEN emb, GEN D, long prec);
GEN qa_sortedembed(GEN Q, GEN order, GEN D, GEN rpell, GEN ncgp, long prec);
GEN qa_sortedembed_typecheck(GEN Q, GEN order, GEN D, int rpell, GEN ncgp, long prec);

//CHECKING METHODS
void qa_check(GEN Q);
void qa_indefcheck(GEN Q);
void qa_eltcheck(GEN x);
GEN qa_ordcheck(GEN ord);
void qa_ordeichlercheck(GEN order);
void QM_check(GEN M);

//PROPERTY RETRIEVAL
GEN qa_getnf(GEN Q);
GEN qa_getpram(GEN Q);
GEN qa_getabvec(GEN Q);
GEN qa_geta(GEN Q);
GEN qa_getb(GEN Q);
GEN qa_getmab(GEN Q);
GEN qa_getpramprod(GEN Q);
GEN qa_getord(GEN order);
GEN qa_getordtype(GEN order);
GEN qa_getordmaxd(GEN order);
GEN qa_getordlevel(GEN order);
GEN qa_getordlevelpfac(GEN order);
GEN qa_getordinv(GEN order);
GEN qa_getordtrace0basis(GEN order);

//SUPPORTING METHODS
int cmp_data(void *data, GEN x, GEN y);
GEN module_intersect(GEN A, GEN B);
GEN module_intersect_typecheck(GEN A, GEN B);
GEN prime_ksearch(GEN relations, GEN extra);
GEN prime_ksearch_typecheck(GEN relations, GEN extra);
GEN QM_hnf(GEN M);
GEN QM_hnf_typecheck(GEN M);
int Q_issquareall(GEN x, GEN *sqrtx);
GEN powerset(GEN L);
GEN vecratio(GEN v1, GEN v2);
GEN vecratio_typecheck(GEN v1, GEN v2);

//qq_qquat_int.c

//METHODS

//INTERSECTION NUMBER BASED OFF OF ROOTS BOUND AREA
GEN qa_inum_roots(GEN Q, GEN order, GEN e1, GEN e2, GEN D1, GEN D2, int data, long prec);
GEN qa_inum_roots_typecheck(GEN Q, GEN order, GEN e1, GEN e2, int data, long prec);

//INTERSECTIONS BASED ON X-LINKAGE
GEN qa_inum_x(GEN Q, GEN order, GEN e1, GEN e2, GEN D1, GEN D2, int data, long prec);
GEN qa_inum_x_typecheck(GEN Q, GEN order, GEN e1, GEN e2, int data, long prec);
GEN qa_xlink(GEN Q, GEN order, GEN e1, GEN e2, GEN D1, GEN D2, GEN x, long prec);
GEN qa_xlink_typecheck(GEN Q, GEN order, GEN e1, GEN e2, GEN x, long prec);
GEN qa_xposs(GEN pset, GEN Psetprod, GEN D1, GEN D2, GEN xmin, GEN xmax);
GEN qa_xposs_typecheck(GEN Qorpset, GEN D1, GEN D2, GEN xmin, GEN xmax);

//qq_visual.c

//METHODS

//HISTOGRAMS
void hist_autocompile(GEN minx, GEN maxx, char *imagename, char *autofile, char *plotoptions, int open);
void hist_compile(char *imagename, char *autoname, int open);
GEN hist_make(GEN data, char *imagename, char *autofile, int compilenew, char *plotoptions, int open, long prec);
GEN hist_tobins(GEN data, GEN minx, GEN maxx, GEN nbins, int toscale, int compilenew, char *imagename, char *autofile, char *plotoptions, int open, long prec);
GEN hist_tobins_defaultbins(GEN data, GEN minx, GEN maxx, int toscale, int compilenew, char *imagename, char *autofile, char *plotoptions, int open, long prec);
GEN hist_rebin(GEN data, GEN histdata, GEN nbins, long prec);
void hist_recompile(GEN histdata);
GEN hist_rerange(GEN data, GEN histdata, GEN minx, GEN maxx, long prec);
GEN hist_rescale(GEN data, GEN histdata, int scale, long prec);

