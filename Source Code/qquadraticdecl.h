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

//INFINITY 
GEN addoo(GEN a, GEN b);
GEN divoo(GEN a, GEN b);

//LINEAR ALGEBRA
GEN FpM_eigenvecs(GEN M, GEN p);
GEN lin_intsolve(GEN A, GEN B, GEN n);
GEN lin_intsolve_tc(GEN A, GEN B, GEN n);
GEN mat3_complete(GEN A, GEN B, GEN C);
GEN mat3_complete_tc(GEN A, GEN B, GEN C);

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
GEN glist_togvec_append(glist *l, GEN v, long length, int dir);
void llist_free(llist *l);
long llist_pop(llist **head_ref);
void llist_putstart(llist **head_ref, long new_data);
GEN llist_togvec(llist *l, long length, int dir);
GEN llist_tovecsmall(llist *l, long length, int dir);

//SHORT VECTORS IN LATTICES
GEN lat_smallvectors(GEN A, GEN C1, GEN C2, GEN condition, int onesign, int isintegral, int rdataonly, long prec);
GEN lat_smallvectors_givendata(GEN chol, GEN U, GEN perminv, GEN C1, GEN C2, GEN condition, int onesign, long prec);
GEN lat_smallvectors_tc(GEN A, GEN C1, GEN C2, int onesign, int isintegral, long prec);
GEN lat_smallvectors_cholesky(GEN Q, GEN C1, GEN C2, GEN condition, int onesign, long prec);
GEN mat_choleskydecomp(GEN A, int rcoefs, long prec);
GEN mat_choleskydecomp_tc(GEN A, int rcoefs, long prec);
GEN mat_uptriag_rowred(GEN M);
GEN mat_uptriag_rowred_tc(GEN M);

//qq_bqf.c

//METHODS

//DISCRIMINANT METHODS
GEN disclist(GEN D1, GEN D2, int fund, GEN cop);
GEN discprimeindex(GEN D, GEN facs);
GEN discprimeindex_tc(GEN D);
int isdisc(GEN D);
GEN pell(GEN D);
GEN pell_tc(GEN D);
GEN posreg(GEN D, long prec);
GEN posreg_tc(GEN D, long prec);
GEN quadroot(GEN D);
GEN quadroot_tc(GEN D);

//BASIC OPERATIONS ON BINARY QUADRATIC FORMS
GEN bqf_automorph_tc(GEN q);
int bqf_compare(void *data, GEN q1, GEN q2);
int bqf_compare_tmat(void *data, GEN d1, GEN d2);
GEN bqf_disc(GEN q);
GEN bqf_disc_tc(GEN q);
int bqf_isreduced(GEN q, int Dsign);
int bqf_isreduced_tc(GEN q);
GEN bqf_random(GEN maxc, int type, int primitive);
GEN bqf_random_D(GEN maxc, GEN D);
GEN bqf_red(GEN q, GEN rootD, int Dsign, int tmat);
GEN bqf_red_tc(GEN q, int tmat, long prec);
GEN bqf_roots(GEN q, GEN D, GEN w);
GEN bqf_roots_tc(GEN q);
GEN bqf_trans(GEN q, GEN mtx);
GEN bqf_trans_tc(GEN q, GEN mtx);
GEN bqf_transL(GEN q, GEN n);
GEN bqf_transR(GEN q, GEN n);
GEN bqf_transS(GEN q);
GEN bqf_trans_coprime(GEN q, GEN n);
GEN bqf_trans_coprime_tc(GEN q, GEN n);

//BASIC METHODS FOR NEGATIVE DISCRIMINANTS
GEN dbqf_automorph(GEN q, GEN D);
GEN dbqf_red(GEN q);
GEN dbqf_red_tmat(GEN q);

//BASIC OPERATIONS SPECIFIC TO INDEFINITE FORMS/POSITIVE DISCRIMINANTS
GEN ibqf_automorph_D(GEN q, GEN D);
GEN ibqf_automorph_pell(GEN q, GEN qpell);
int ibqf_isrecip(GEN q, GEN rootD);
int ibqf_isrecip_tc(GEN q, long prec);
GEN ibqf_leftnbr(GEN q, GEN rootD);
GEN ibqf_leftnbr_tmat(GEN q, GEN rootD);
GEN ibqf_leftnbr_tc(GEN q, int tmat, long prec);
GEN ibqf_leftnbr_update(GEN qvec, GEN rootD);
GEN ibqf_red(GEN q, GEN rootD);
GEN ibqf_red_tmat(GEN q, GEN rootD);
GEN ibqf_red_pos(GEN q, GEN rootD);
GEN ibqf_red_pos_tmat(GEN q, GEN rootD);
GEN ibqf_redorbit(GEN q, GEN rootD);
GEN ibqf_redorbit_tmat(GEN q, GEN rootD);
GEN ibqf_redorbit_posonly(GEN q, GEN rootD);
GEN ibqf_redorbit_posonly_tmat(GEN q, GEN rootD);
GEN ibqf_redorbit_tc(GEN q, int tmat, int posonly, long prec);
GEN ibqf_rightnbr(GEN q, GEN rootD);
GEN ibqf_rightnbr_tmat(GEN q, GEN rootD);
GEN ibqf_rightnbr_update(GEN qvec, GEN rootD);
GEN ibqf_rightnbr_tc(GEN q, int tmat, long prec);
GEN ibqf_river(GEN q, GEN rootD);
GEN ibqf_river_positions(GEN q, GEN rootD);
GEN ibqf_river_positions_forms(GEN q, GEN rootD);
GEN ibqf_river_tc(GEN q, long prec);
GEN ibqf_riverforms(GEN q, GEN rootD);
GEN ibqf_riverforms_tc(GEN q, long prec);
GEN ibqf_symmetricarc(GEN q, GEN D, GEN rootD, GEN qpell, long prec);
GEN ibqf_symmetricarc_tc(GEN q, long prec);
GEN ibqf_toriver(GEN q, GEN rootD);
GEN ibqf_toriver_tmat(GEN q, GEN rootD);
GEN mat_toibqf(GEN mtx);
GEN mat_toibqf_tc(GEN mtx);

//EQUIVALENCE OF BQFs
GEN bqf_isequiv(GEN q1, GEN q2, GEN rootD, int Dsign, int tmat);
GEN bqf_isequiv_set(GEN q, GEN S, GEN rootD, int Dsign, int tmat);
GEN bqf_isequiv_tc(GEN q1, GEN q2, int tmat, long prec);
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
GEN bqf_comp_tc(GEN q1, GEN q2, int tored, long prec);
GEN bqf_idelt(GEN D);
GEN bqf_ncgp(GEN D, long prec);
GEN bqf_ncgp_lexic(GEN D, long prec);
GEN bqf_pow(GEN q, GEN n);
GEN bqf_pow_red(GEN q, GEN n, GEN rootD, int Dsign);
GEN bqf_pow_tc(GEN q, GEN n, int tored, long prec);
GEN bqf_square(GEN q);
GEN bqf_square_red(GEN q, GEN rootD, int Dsign);
GEN bqf_square_tc(GEN q, int tored, long prec);
GEN ideal_tobqf(GEN numf, GEN ideal);

//REPRESENTATION OF NUMBERS BY BQFs
GEN bqf_reps(GEN q, GEN n, int proper, int half, long prec);
GEN bqf_reps_tc(GEN q, GEN n, int proper, int half, long prec);
GEN dbqf_reps(GEN qred, GEN D, GEN n, int proper, int half);
GEN ibqf_reps(GEN qorb, GEN qautom, GEN D, GEN rootD, GEN n, int proper, int half);
GEN sbqf_reps(GEN q, GEN D, GEN rootD, GEN n, int half);
GEN zbqf_reps(GEN A, GEN B, GEN n, int half);

//MORE REPRESENTATION OF NUMBERS
GEN bqf_bigreps(GEN q, GEN n, long prec);
GEN bqf_bigreps_tc(GEN q, GEN n, long prec);
GEN bqf_linearsolve(GEN q, GEN n1, GEN lin, GEN n2, long prec);
GEN bqf_linearsolve_tc(GEN q, GEN n1, GEN lin, GEN n2, long prec);

//GENERAL CHECKING METHODS
void bqf_check(GEN q);
GEN bqf_checkdisc(GEN q);
void intmatrix_check(GEN mtx);


//qq_bqf_int.c

//METHODS

//INTERSECTION DATA
GEN bqf_bdelta(GEN q1, GEN q2);
GEN bqf_bdelta_tc(GEN q1, GEN q2);
GEN bqf_intlevel(GEN q1, GEN q2);
GEN bqf_intlevel_tc(GEN q1, GEN q2);
GEN ibqf_intpairs_transtoq(GEN pairs, GEN q, GEN rootD);
GEN ibqf_intpoint(GEN q1, GEN q2, GEN location, GEN autom);
GEN ibqf_intpoint_tc(GEN q1, GEN q2, GEN location);

//INTERSECTION NUMBER COMPUTATION
GEN ibqf_int(GEN r1, GEN r2);
GEN ibqf_int_tc(GEN q1, GEN q2, long prec);
GEN ibqf_intRS_byriver(GEN r1, GEN r2);
GEN ibqf_intRS_tc(GEN q1, GEN q2, long prec);
GEN ibqf_intforms_byriver(GEN r1, GEN r2, int data, int split);
GEN ibqf_intforms_tc(GEN q1, GEN q2, int data, int split, long prec);
GEN ibqf_intformsRS_byriver(GEN r1, GEN r2, int data);
GEN ibqf_intformsRS_tc(GEN q1, GEN q2, int data, long prec);
GEN ibqf_intformsRO_byriver(GEN r1, GEN r2, int data);
GEN ibqf_intformsRO_tc(GEN q1, GEN q2, int data, long prec);
GEN ibqf_intformsLS_byriver(GEN r1, GEN r2, int data);
GEN ibqf_intformsLS_tc(GEN q1, GEN q2, int data, long prec);
GEN ibqf_intformsLO_byriver(GEN r1, GEN r2, int data);
GEN ibqf_intformsLO_tc(GEN q1, GEN q2, int data, long prec);


//qq_bqf_int_scripts.c

//METHODS

//OLD INTERSECTION ALGORITHMS
GEN ibqf_intRS_byriverOLDBAD(GEN r1, GEN r2);
GEN ibqf_intRSOLD(GEN q1, GEN q2, long prec);

//INUM DATA
GEN ibqf_int_dists(GEN q1, GEN q2, GEN rootD1, GEN q1auto, GEN z, int ispairs, long prec);
GEN ibqf_int_dists_tc(GEN q1, GEN q2, long prec);
GEN ibqf_intangle(GEN q1, GEN q2, GEN D1D2, long prec);
GEN ibqf_int_angles(GEN q1, GEN q2, GEN D1, GEN D2, int ispairs, long prec);
GEN ibqf_intsets(GEN S1, GEN S2, long prec);

//COMPILING DATA
void ibqf_int_timings(GEN D1range, GEN D2range, GEN filename, long trials, long prec);
GEN ibqf_int_dists_histogram(GEN q1, GEN D2, char *filename, char *autocompile, long prec);
GEN ibqf_int_angles_histogram(GEN q1, GEN D2, char *filename, char *autocompile, long prec);
GEN ibqf_CD1D2data(GEN D1range, GEN D2range, long trials, int scaled, char *filename, char *autocompile, long prec);

//PLOTTING INTERSECTION NUMBERS
GEN ibqf_plotint_arc(GEN q1, GEN q2, GEN filename, long prec);
GEN ibqf_plotint_bydist(GEN q1, GEN q2, GEN filename, long prec);
GEN ibqf_plotint_fd(GEN q1, GEN q2, GEN filename, long prec);


//qq_geometry

//METHODS

//BASIC LINE, CIRCLE, AND POINT OPERATIONS
GEN arc_init(GEN c, GEN p1, GEN p2, int dir, long prec);
GEN arc_init_tc(GEN c, GEN p1, GEN p2, int dir, long prec);
GEN arc_midpoint(GEN c, GEN p1, GEN p2, GEN tol, long prec);
GEN arc_midpoint_tc(GEN c, GEN p1, GEN p2, long prec);
GEN circle_angle(GEN c1, GEN c2, GEN p, GEN tol, long prec);
GEN circle_angle_tc(GEN c1, GEN c2, GEN p, long prec);
GEN circle_fromcp(GEN cent, GEN p, long prec);
GEN circle_fromppp(GEN p1, GEN p2, GEN p3, GEN tol, long prec);
GEN circle_fromppp_tc(GEN p1, GEN p2, GEN p3, long prec);
GEN circle_tangentslope(GEN c, GEN p, long prec);
GEN circle_tangentslope_tc(GEN c, GEN p, long prec);
GEN crossratio(GEN a, GEN b, GEN c, GEN d);
GEN line_angle(GEN l1, GEN l2, long prec);
GEN line_fromsp(GEN s, GEN p);
GEN line_frompp(GEN p1, GEN p2);
GEN mat_eval(GEN M, GEN x);
GEN mat_eval_tc(GEN M, GEN x);
GEN midpoint(GEN p1, GEN p2);
GEN mobius(GEN M, GEN c, GEN tol, long prec);
GEN mobius_tc(GEN M, GEN c, long prec);
GEN perpbis(GEN p1, GEN p2);
GEN radialangle(GEN c, GEN p, GEN tol, long prec);
GEN radialangle_tc(GEN c, GEN p, long prec);
GEN slope(GEN p1, GEN p2);

//INTERSECTION OF LINES/CIRCLES
GEN arc_int(GEN c1, GEN c2, GEN tol, long prec);
GEN arc_int_tc(GEN c1, GEN c2, long prec);
GEN arcseg_int(GEN c, GEN l, GEN tol, long prec);
GEN arcseg_int_tc(GEN c, GEN l, long prec);
GEN circle_int(GEN c1, GEN c2, GEN tol, long prec);
GEN circle_int_tc(GEN c1, GEN c2, long prec);
GEN circleline_int(GEN c, GEN l, GEN tol, long prec);
GEN circleline_int_tc(GEN c, GEN l, long prec);
GEN genseg_int(GEN s1, GEN s2, GEN tol, long prec);
GEN genseg_int_tc(GEN s1, GEN s2, long prec);
GEN line_int(GEN l1, GEN l2, GEN tol, long prec);
GEN line_int_tc(GEN l1, GEN l2, long prec);
int onarc(GEN c, GEN p, GEN tol, long prec);
int onarc_tc(GEN c, GEN p, long prec);
int onseg(GEN l, GEN p, GEN tol, long prec);
int onseg_tc(GEN l, GEN p, long prec);
GEN seg_int(GEN l1, GEN l2, GEN tol, long prec);
GEN seg_int_tc(GEN l1, GEN l2, long prec);

//DISTANCES
GEN hdist(GEN z1, GEN z2, long prec);
GEN hdist_tc(GEN z1, GEN z2, long prec);
GEN hdist_ud(GEN z1, GEN z2, long prec);
GEN hpolygon_area(GEN circles, GEN vertices, GEN tol, long prec);
GEN hpolygon_area_tc(GEN circles, GEN vertices, long prec);

//FUNDAMENTAL DOMAIN COMPUTATION
GEN edgepairing(GEN U, GEN tol, int rboth, long prec);
GEN edgepairing_tc(GEN U, long prec);
GEN isometriccircle_mats(GEN g, GEN mats, GEN *data, GEN (*gamtopsl)(GEN *, GEN, long), GEN tol, long prec);
GEN isometriccircle_psl(GEN g, GEN p, GEN tol, long prec);
GEN isometriccircle_psl_mats(GEN g, GEN mats, GEN tol, long prec);
GEN isometriccircle_psu(GEN g, GEN tol, long prec);
GEN normalizedbasis(GEN G, GEN U, GEN mats, GEN gamid, GEN *data, GEN (*gamtopsl)(GEN *, GEN, long), GEN (*eltmul)(GEN *, GEN, GEN), GEN (*eltinv)(GEN *, GEN), int (*istriv)(GEN *, GEN), GEN tol, long prec);
GEN normalizedboundary_append(GEN Ubase, GEN G, GEN mats, GEN id, GEN tol, long prec);
GEN normalizedboundary_givenU(GEN Ubase, GEN G, GEN mats, GEN id, GEN *data, GEN (*gamtopsl)(GEN *, GEN, long), GEN tol, long prec);
GEN normalizedboundary_givencircles(GEN G, GEN mats, GEN id, GEN tol, long prec);
GEN normalizedboundary(GEN G, GEN mats, GEN id, GEN *data, GEN (*gamtopsl)(GEN *, GEN, long), GEN tol, long prec);
long normalizedboundary_outside(GEN U, GEN z, GEN tol, long prec);
long normalizedboundary_outside_tc(GEN U, GEN z, long prec);
GEN normalizedboundary_sideint(GEN U, GEN c, int start, GEN tol, long prec);
GEN psltopsu(GEN g, GEN p);
GEN psltopsu_mats(GEN g, GEN M);
GEN psltopsu_transmats(GEN p);
GEN psl_roots(GEN M, GEN tol, long prec);
GEN randompoint_ud(GEN R, long prec);
GEN reduceelt_givennormbound(GEN U, GEN g, GEN z, GEN gamid, GEN *data, GEN (*gamtopsl)(GEN *, GEN, long), GEN (*eltmul)(GEN *, GEN, GEN), GEN tol, long prec);
GEN reduceelt_givenpsu(GEN G, GEN Gmats, GEN g, GEN gmat, GEN z, GEN gamid, GEN *data, GEN (*eltmul)(GEN *, GEN, GEN), GEN tol, long prec);
GEN reducepoint(GEN U, GEN z, GEN gamid, GEN *data, GEN (*eltmul)(GEN *, GEN, GEN), GEN tol, long prec);
GEN rootgeodesic_fd(GEN U, GEN g, GEN gamid, GEN *data, GEN (*gamtopsl)(GEN *, GEN, long), GEN (*eltmul)(GEN *, GEN, GEN), GEN (*eltinv)(GEN *, GEN), GEN tol, long prec);
GEN rootgeodesic_ud(GEN M, GEN mats, GEN tol, long prec);
GEN rootgeodesic_uhp(GEN M, GEN tol, long prec);
GEN rootgeodesic_uhp_tc(GEN M, long prec);

//PRINTING TO PLOTVIEWER
void python_printarcs(GEN arcs, char *filename, int view, char *extrainput, long prec);
void python_plotviewer(char *input);
void python_printfdom(GEN U, char *filename, long prec);

//HELPER METHODS
GEN anglediff(GEN ang, GEN bot, GEN tol, long prec);
GEN atanoo(GEN x, long prec);
GEN deftol(long prec);
int gcmp_strict(void *data, GEN x, GEN y);
int geom_check(GEN c);
GEN shiftangle(GEN ang, GEN bot, GEN tol, long prec);
long tolcmp(GEN x, GEN y, GEN tol, long prec);
int tolcmp_sort(void *data, GEN x, GEN y);
int toleq(GEN x, GEN y, GEN tol, long prec);


//qq_quat.c

//METHODS

//BASIC OPERATIONS ON ELEMENTS IN QUATERNION ALGEBRAS
GEN qa_conj(GEN x);
GEN qa_conj_tc(GEN x);
GEN qa_conjby(GEN Q, GEN x, GEN y);
GEN qa_conjby_tc(GEN Q, GEN x, GEN y);
GEN qa_inv(GEN Q, GEN x);
GEN qa_inv_tc(GEN Q, GEN x);
GEN qa_m2rembed(GEN Q, GEN x);
GEN qa_m2rembed_tc(GEN Q, GEN x);
GEN qa_minpoly(GEN Q, GEN x);
GEN qa_minpoly_tc(GEN Q, GEN x);
GEN qa_mul(GEN Q, GEN x, GEN y);
GEN qa_mul_tc(GEN Q, GEN x, GEN y);
GEN qa_mulvec(GEN Q, GEN L);
GEN qa_mulvec_tc(GEN Q, GEN L);
GEN qa_mulvecindices(GEN Q, GEN L, GEN indices);
GEN qa_mulvecindices_tc(GEN Q, GEN L, GEN indices);
GEN qa_norm(GEN Q, GEN x);
GEN qa_norm_tc(GEN Q, GEN x);
GEN qa_pow(GEN Q, GEN x, GEN n);
GEN qa_pow_tc(GEN Q, GEN x, GEN n);
GEN qa_roots(GEN Q, GEN x, long prec);
GEN qa_roots_tc(GEN Q, GEN x, long prec);
GEN qa_square(GEN Q, GEN x);
GEN qa_square_tc(GEN Q, GEN x);
GEN qa_trace(GEN x);
GEN qa_trace_tc(GEN x);

//BASIC OPERATIONS ON ORDERS/LATTICES IN QUATERNION ALGEBRAS
int qa_isinorder(GEN Q, GEN ordinv, GEN x);
int qa_isinorder_tc(GEN Q, GEN ord, GEN x);
int qa_isorder(GEN Q, GEN ord, GEN ordinv);
int qa_isorder_tc(GEN Q, GEN ord);
GEN qa_leftorder(GEN Q, GEN L, GEN Linv);
GEN qa_leftorder_tc(GEN Q, GEN L);
GEN qa_rightorder(GEN Q, GEN L, GEN Linv);
GEN qa_rightorder_tc(GEN Q, GEN L);
GEN qa_ord_conj(GEN Q, GEN ord, GEN c);
GEN qa_ord_conj_tc(GEN Q, GEN ord, GEN c);
GEN qa_ord_disc(GEN Q, GEN ord);
GEN qa_ord_disc_tc(GEN Q, GEN ord);
GEN qa_ord_normform(GEN Q, GEN ord);
GEN qa_superorders(GEN Q, GEN ord, GEN n);
GEN qa_superorders_prime(GEN Q, GEN ord, GEN ordinv, GEN n);
GEN qa_superorders_tc(GEN Q, GEN ord, GEN n);

//INITIALIZATION METHODS
GEN qa_eichlerorder(GEN Q, GEN l, GEN maxord);
GEN qa_eichlerorder_tc(GEN Q, GEN l, GEN maxord);
GEN qa_ord_init(GEN Q, GEN ord);
GEN qa_ord_init_tc(GEN Q, GEN ord);
GEN qa_init_ab(GEN a, GEN b);
GEN qa_init_ab_tc(GEN a, GEN b);
GEN qa_init_primes(GEN pset, int type);
GEN qa_init_primes_tc(GEN pset);
GEN qa_init_2primes(GEN p, GEN q);
GEN qa_init_2primes_tc(GEN p, GEN q);
GEN qa_maximalorder(GEN Q, GEN baseord);
GEN qa_maximalorder_tc(GEN Q, GEN baseord);
GEN qa_ram_fromab(GEN a, GEN b);
GEN qa_ram_fromab_tc(GEN a, GEN b);

//CONJUGATION OF ELEMENTS IN A GIVEN ORDER
GEN qa_conjbasis(GEN Q, GEN ord, GEN ordinv, GEN e1, GEN e2, int orient);
GEN qa_conjbasis_tc(GEN Q, GEN ord, GEN e1, GEN e2, int orient);
GEN qa_conjqf(GEN Q, GEN ord, GEN ordinv, GEN e1, GEN e2);
GEN qa_conjqf_tc(GEN Q, GEN ord, GEN e1, GEN e2);
GEN qa_conjnorm(GEN Q, GEN ord, GEN ordinv, GEN e1, GEN e2, GEN n, int retconelt, long prec);
GEN qa_conjnorm_tc(GEN Q, GEN ord, GEN e1, GEN e2, GEN n, int retconelt, long prec);
GEN qa_simulconj(GEN Q, GEN ord, GEN ordinv, GEN e1, GEN e2, GEN f1, GEN f2, long prec);
GEN qa_simulconj_tc(GEN Q, GEN ord, GEN e1, GEN e2, GEN f1, GEN f2, long prec);

//EMBEDDING QUADRATIC ORDERS INTO EICHLER ORDERS
GEN qa_associatedemb(GEN Q, GEN order, GEN emb, GEN D);
GEN qa_associatedemb_tc(GEN Q, GEN order, GEN emb, GEN D);
GEN qa_embed(GEN Q, GEN order, GEN D, GEN nembeds, GEN rpell, long prec);
GEN qa_embed_tc(GEN Q, GEN order, GEN D, GEN nembeds, int rpell, long prec);
GEN qa_embeddablediscs(GEN Q, GEN order, GEN d1, GEN d2, int fund, GEN cop);
GEN qa_embeddablediscs_tc(GEN Q, GEN order, GEN d1, GEN d2, int fund, GEN cop);
GEN qa_numemb(GEN Q, GEN order, GEN D, GEN narclno);
GEN qa_numemb_tc(GEN Q, GEN order, GEN D, GEN narclno, long prec);
GEN qa_ordiffer(GEN Q, GEN order, GEN e1, GEN e2, GEN D);
GEN qa_ordiffer_tc(GEN Q, GEN order, GEN e1, GEN e2, GEN D);
GEN qa_orinfinite(GEN Q, GEN emb, GEN D, long prec);
GEN qa_orinfinite_tc(GEN Q, GEN emb, GEN D, long prec);
GEN qa_sortedembed(GEN Q, GEN order, GEN D, GEN rpell, GEN ncgp, long prec);
GEN qa_sortedembed_tc(GEN Q, GEN order, GEN D, int rpell, GEN ncgp, long prec);

//FUNDAMENTAL DOMAIN METHODS
GEN qa_fundamentaldomain(GEN Q, GEN order, GEN p, int dispprogress, GEN ANRdata, GEN tol, long prec);
GEN qa_fundamentaldomain_tc(GEN Q, GEN order, GEN p, int dispprogress, GEN ANRdata, long prec);
GEN qa_invradqf(GEN Q, GEN order, GEN mats, GEN z, long prec);
GEN qa_isometriccircle(GEN Q, GEN x, GEN mats, GEN tol, long prec);
GEN qa_isometriccircle_tc(GEN Q, GEN x, GEN p, long prec);
GEN qa_fdarea(GEN Q, GEN order, long prec);
GEN qa_fdarea_tc(GEN Q, GEN order, long prec);
GEN qa_fdm2rembed(GEN *data, GEN x, long prec);
GEN qa_fdinv(GEN *data, GEN x);
GEN qa_fdmul(GEN *data, GEN x, GEN y);
int qa_istriv(GEN *data, GEN x);
GEN qa_normalizedbasis(GEN Q, GEN G, GEN mats, GEN U, GEN tol, long prec);
GEN qa_normalizedbasis_tc(GEN Q, GEN G, GEN p, long prec);
GEN qa_normalizedboundary(GEN Q, GEN G, GEN mats, GEN tol, long prec);
GEN qa_normalizedboundary_tc(GEN Q, GEN G, GEN p, long prec);
void qa_printisometriccircles(GEN Q, GEN L, char *filename, GEN mats, GEN tol, long prec);
void qa_printisometriccircles_tc(GEN Q, GEN L, GEN p, char *filename, int view, long prec);
GEN qa_reduceelt(GEN Q, GEN G, GEN x, GEN z, GEN p, GEN tol, long prec);
GEN qa_reduceelt_normbound(GEN Q, GEN U, GEN x, GEN z, GEN tol, long prec);
GEN qa_reduceelt_tc(GEN Q, GEN G, GEN x, GEN z, GEN p, long prec);
GEN qa_rootgeodesic_fd(GEN Q, GEN U, GEN g, GEN tol, long prec);
GEN qa_rootgeodesic_fd_tc(GEN Q, GEN U, GEN g, long prec);
GEN qa_smallnorm1elts_invrad(GEN Q, GEN order, GEN C1, GEN C2, GEN invrad, long prec);
GEN qa_smallnorm1elts_tc(GEN Q, GEN order, GEN p, GEN z, GEN C1, GEN C2, long prec);
GEN qa_topsu(GEN Q, GEN g, GEN p, long prec);
GEN qa_topsu_mat(GEN Q, GEN g, GEN mats, long prec);
GEN qa_topsu_tc(GEN Q, GEN g, GEN p, long prec);

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
GEN module_intersect_tc(GEN A, GEN B);
GEN prime_ksearch(GEN relations, GEN extra);
GEN prime_ksearch_tc(GEN relations, GEN extra);
GEN QM_hnf(GEN M);
GEN QM_hnf_tc(GEN M);
int Q_issquareall(GEN x, GEN *sqrtx);
GEN powerset(GEN L);
GEN powerset_tc(GEN L);
GEN vecratio(GEN v1, GEN v2);
GEN vecratio_tc(GEN v1, GEN v2);


//qq_quat_int.c

//METHODS

//INTERSECTION NUMBER BASED OFF OF ROOTS BOUND AREA
GEN qa_inum_roots(GEN Q, GEN order, GEN e1, GEN e2, GEN D1, GEN D2, int data, long prec);
GEN qa_inum_roots_tc(GEN Q, GEN order, GEN e1, GEN e2, int data, long prec);

//INTERSECTIONS BASED ON X-LINKAGE
GEN qa_inum_x(GEN Q, GEN order, GEN e1, GEN e2, GEN D1, GEN D2, int data, long prec);
GEN qa_inum_x_tc(GEN Q, GEN order, GEN e1, GEN e2, int data, long prec);
GEN qa_xlink(GEN Q, GEN order, GEN e1, GEN e2, GEN D1, GEN D2, GEN x, long prec);
GEN qa_xlink_tc(GEN Q, GEN order, GEN e1, GEN e2, GEN x, long prec);
GEN qa_xposs(GEN pset, GEN Psetprod, GEN D1, GEN D2, GEN xmin, GEN xmax);
GEN qa_xposs_tc(GEN Qorpset, GEN D1, GEN D2, GEN xmin, GEN xmax);

//INTERSECTION NUMBER BASED ON FUNDAMENTAL DOMAIN
GEN qa_inum_fd_givengeod(GEN Q, GEN order, GEN U, GEN geod1, GEN geod2, GEN pell1, GEN pell2, GEN D1D2, int data, GEN tol, long prec);
GEN qa_inum_fd_tc(GEN Q, GEN order, GEN U, GEN e1, GEN e2, int data, long prec);

//INTERSECTION DATA
GEN qa_intlevel(GEN Q, GEN order, GEN e1, GEN e2, GEN D1D2, long prec);
GEN qa_intlevel_tc(GEN Q, GEN order, GEN e1, GEN e2, GEN D1, GEN D2, long prec);


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

