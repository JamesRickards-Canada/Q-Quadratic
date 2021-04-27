print("\n\nType '?qq' for help.\n\n");
addhelp(qq, "For each package P, call ?P to access a basic description and list of methods. Installed packages: \n base \n bqf \n bqf_int \n geometry \n quat \n quat_int \n visual");
default(help, "gphelp -detex");

\\qq_base

	\\INFINITY
		install("addoo","GG","addoo","./libqquadratic.so");
		addhelp(addoo, "Inputs a,b real numbers or oo.\n Outputs a+b, where if a or b is +/- oo, returns that back. Note that this will make oo+-oo=oo and -oo+oo=-oo");
		install("divoo","GG","divoo","./libqquadratic.so");
		addhelp(divoo, "Inputs a,b, real numbers or oo.\n Outputs a/b, where oo is output if a=oo and b>=0 or a=-oo and b<0 or b=0 and a>=0 (outputs -oo under analogous assumptions).");

	\\LINEAR ALGEBRA
		install("lin_intsolve_tc","GGG","lin_intsolve","./libqquadratic.so");
		addhelp(lin_intsolve, "Inputs A,B,n integers.\n Outputs the general integral solutions to Ax+By=n. The format is [[s1,s2],[x0,y0]], where the general solution is x=s1*t+x0, y=s2*t+y0 for t an integer. The output is also reduced, i.e. gcd(s1,s2)=1. If A=B=0 or there are no integer solutions, returns 0.");
		install("mat3_complete_tc", "GGG", "mat3_complete", "./libqquadratic.so");
		addhelp(mat3_complete, "Inputs A,B,C integers with gcd 1.\n Outputs a 3x3 integer matrix with determinant 1 whose top row is [A, B, C].");

	\\TIME
		install("printtime","v","printtime","./libqquadratic.so");
		addhelp(printtime, "Prints current time");

	\\SHORT VECTORS IN LATTICES
		install("lat_smallvectors_tc","GGD0,G,D1,L,D0,L,p","lat_smallvectors","./libqquadratic.so");
		addhelp(lat_smallvectors,"Inputs A, C1, {C2=0}, {onesign=1}, {isintegral=0}: positive definite symmetric matrix A, nonnegative constants C1 and C2, onesign=0,1, isintegral=0,1.\n Outputs all non-zero x such that C1<x^TAx<=C. If C2=0, instead outputs non-zero x such that x^TAx<=C1. If onesign=1, only output one of x,-x for each solution x. If the norms are always integral, (entries of A are half integers and integers on the diagonal), then pass isintegral as 1 and the values of the output are fixed to being exact integers. It is advised to pass C+0.1 or similar into the method to account for rounding errors with q(x)=C (correct solutions may be discarded otherwise).");
		install("mat_choleskydecomp_tc","GD1,L,p","mat_choleskydecomp","./libqquadratic.so");
		addhelp(mat_choleskydecomp,"Input A, {rcoefs=1}: an nxn symmetric matrix A, and rcoefs=0,1.\n If rcoefs=1, outputs B, an n x n matrix, such that the quadratic form x^TAx is expressible as sum(i=1..n)b_ii(x_i+sum(j=i+1..n)b_ijxj)^2. Otherwise, returns the corresponding matrix R, i.e. so that R^TR=A. Note that if A is not positive definite, this will return a complex matrix.");

	\\GENERAL HELP
		addhelp(base,"This package is a collection of miscellaneous methods that may be useful in a variety of settings, and not just for the programs they were originally created for \n Subtopics: \n Infinity (inf) \n Linear algebra (la) \n Random (rand) \n Square roots (root) \n Time (time)");
		addhelp(inf,"addoo, divoo.");
		addhelp(la,"lat_smallvectors, lin_intsolve, mat_choleskydecomp, mat3_complete.");
		addhelp(rand,"rand_elt.");
		addhelp(time,"printtime.");


\\qq_bqf

	\\DISCRIMINANT METHODS
		install("disclist","GGD0,L,D0,G,","disclist","./libqquadratic.so");
		addhelp(disclist, "Inputs d1, d2, {fund=0}, {cop=0}: d1 and d2 integers with d1<=d2, fund=0,1, cop integer.\n Returns the set of proper discriminants between d1 and d2, inclusive. If fund=1, only returns fundamental discriminants. If cop!=0, only returns discriminants coprime to cop.");
		install("discprimeindex_tc","G","discprimeindex","./libqquadratic.so");
		addhelp(discprimeindex, "Inputs: D, a proper discriminant.\n Returns all prime divisors p of D for which D/p^2 is a proper discriminant.");
		install("isdisc","iG","isdisc","./libqquadratic.so");
		addhelp(isdisc, "Inputs: D a real number.\n Returns 1 if D is a proper discriminant, and 0 otherwise.");
		install("pell_tc","G","pell","./libqquadratic.so");
		addhelp(pell, "Inputs: D a positive discriminant.\n Returns [T, U], which is the smallest positive integer solution to T^2-DU^2=4 (and so (T+Usqrt(D))/2 is the fundamental unit in O_D).");
		install("posreg_tc","Gp","posreg","./libqquadratic.so");
		addhelp(posreg, "Inputs: D a positive discriminant.\n Returns the positive regulator of O_D, i.e. the logarithm of the fundamental totally positive unit.");
		install("quadroot_tc","G","quadroot","./libqquadratic.so");
		addhelp(quadroot, "Input D, a non-square integer.\n Returns sqrt(D) of type t_QUAD.");

	\\BASIC OPERATIONS ON BINARY QUADRATIC FORMS
		install("bqf_automorph_tc","G","bqf_automorph","./libqquadratic.so");
		addhelp(bqf_automorph, "Inputs: q a BQF.\n Returns a generator of the automorph group of q in PSL(2,Z).");
		install("bqf_disc_tc","G","bqf_disc","./libqquadratic.so");
		addhelp(bqf_disc, "Inputs: q, quadratic form.\n Returns the discriminant of q.");
		install("bqf_isequiv_tc","GGD0,L,p","bqf_isequiv","./libqquadratic.so");
		addhelp(bqf_isequiv, "Inputs: q, S, {tmat=0}: q a BQF, S either a BQF or a set of BQFs, tmat=0,1.\n This method tests if q is PSL(2,Z) equivalent to S or any form in S. If S is a form, this returns 1 if equivalent and 0 if not (if tmat!=0, returns a possible transition matrix).\n If S is a set of forms, this returns 0 if not equivalent and an index i such that q is equivalent to S[i] otherwise. If tmat!=0, this returns [index, transition matrix].");
		install("bqf_isreduced_tc","iG","bqf_isreduced","./libqquadratic.so");
		addhelp(bqf_isreduced,"Input q, quadratic form.\n Returns 1 if q is reduced, and 0 if not.");
		install("bqf_random","GD0,L,D1,L,","bqf_random","./libqquadratic.so");
		addhelp(bqf_random,"Inputs maxc, {type=0}, {primitive=1}; maxc a positive integer, type=-1,0,1, and primitive=0,1.\n Returns a random BQF with coefficients bounded by maxc. If type=-1 it is positive definite, =1 is indefinite, and =0 means either. If primitive=1 the form is primitive, else it doesn't have to be.");
		install("bqf_random_D","GG","bqf_random_D","./libqquadratic.so");
		addhelp(bqf_random_D,"Inputs maxc, D: maxc a positive integer, and D a discriminant.\n Returns a random primitive form (positive definite if D<0) of discriminant D whose B coefficient is bounded by maxc.");
		install("bqf_red_tc","GD0,L,p","bqf_red","./libqquadratic.so");
		addhelp(bqf_red, "Inputs: q, {tmat=0}: BQF q, (tmat=0,1).\n Returns a reduced form equivalent to q, and if tmat!=0, we return [q_red, transition matrix].");
		install("bqf_roots_tc","G","bqf_roots","./libqquadratic.so");
		addhelp(bqf_roots, "Inputs q: quadratic form q.\n Returns the roots of q with the first root first.");
		install("bqf_trans_tc","GG","bqf_trans","./libqquadratic.so");
		addhelp(bqf_trans, "Inputs q, mtx: integral quadratic form q, 2x2 integral matrix mtx.\n Returns the form acquired by replacing (x,y)^T with m(x,y)^T.");
		install("bqf_trans_coprime_tc", "GG", "bqf_trans_coprime", "./libqquadratic.so");
		addhelp(bqf_trans_coprime,"Inputs q, n: q a primitive integral BQF, and n an integer.\n Returns a form similar to q whose first coefficient is coprime to n.");
		install("ideal_tobqf","GG","ideal_tobqf","./libqquadratic.so");
		addhelp(ideal_tobqf,"Inputs nf, ideal: a quadratic number field nf with ideal ideal.\n Returns the corresponding binary quadratic form.");

	\\BASIC OPERATIONS SPECIFIC TO INDEFINITE FORMS
		install("ibqf_isrecip_tc","iGp","ibqf_isrecip","./libqquadratic.so");
		addhelp(ibqf_isrecip,"Inputs: q, a PIBQF.\n Returns 1 if q is q is reciprocal, and 0 otherwise.");
		install("ibqf_leftnbr_tc","GD0,L,p","ibqf_leftnbr","./libqquadratic.so");
		addhelp(ibqf_leftnbr, "Inputs q, {tmat=0}: an indefinite binary quadratic form on the river and tmat=0,1. \n Returns q' or [q',mat] (if tmat=0,1 respectively), where q' is the left neighbour of q and mat is the bqf_transition matrix from q to q'. The left neighbour is the previous form along the flow of the river that is reduced (AC<0 and B>|A+C|, and occurs when the branches swap from being below to above or vice versa).");
		install("ibqf_redorbit_tc","GD0,L,D0,L,p","ibqf_redorbit","./libqquadratic.so");
		addhelp(ibqf_redorbit, "Inputs q, (tmat), (posonly): q a PIBQF, tmat and posonly=0,1.\n Returns the reduced orbit of q. If tmat=1, also returns the corresponding transition matrices, and if posonly=1 only returns the reduced forms with A>0.");
		install("ibqf_rightnbr_tc","GD0,L,p","ibqf_rightnbr","./libqquadratic.so");
		addhelp(ibqf_rightnbr, "Inputs q, (tmat): an indefinite binary quadratic form on the river and tmat=0,1. \n Returns q' or [q',mat] (if tmat=0,1 respectively), where q' is the right neighbour of q and mat is the bqf_transition matrix from q to q'. The right neighbour is the next form along the flow of the river that is reduced (AC<0 and B>|A+C|, and occurs when the branches swap from being below to above or vice versa).");
		install("ibqf_river_tc","Gp","ibqf_river","./libqquadratic.so");
		addhelp(ibqf_river, "Input: q an indefinite quadratic form.\n Returns the river sequence corresponding to q, where a 1 corresponds to going right and 0 corresponds to going left.");
		install("ibqf_riverforms_tc","Gp","ibqf_riverforms","./libqquadratic.so");
		addhelp(ibqf_riverforms, "Input q a PIBQF.\n This calculates all forms on the river of q, and returns those with A>0, in the order that they appear on the river.");
		install("ibqf_symmetricarc_tc","Gp","ibqf_symmetricarc","./libqquadratic.so");
		addhelp(ibqf_symmetricarc,"Input q, a PIBQF.\n Returns [z,gamma_q(z)] on the root geodesic corresponding to q so that q,gamma_q are symmetric about the arc.");
		install("mat_toibqf_tc","G","mat_toibqf","./libqquadratic.so");
		addhelp(mat_toibqf, "Inputs: mtx, a hyperbolic matrix in SL(2,Z).\n This returns the PIBQF for which it is the ibqf_automorph, namely [c,d-a,-b]/gcd(c,d-a,b) if mtx=[a,b;c,d].");

	\\CLASS GROUPS AND COMPOSITION OF FORMS
		install("bqf_comp_tc","GGD1,L,p","bqf_comp","./libqquadratic.so");
		addhelp(bqf_comp,"Inputs q1, q2, {tored=1}: BQFs q1, q2 of the same discriminant, tored=0, 1.\n Returns the composition of q1 and q2, reduced if tored=1.");
		install("bqf_ncgp","Gp","bqf_ncgp","./libqquadratic.so");
		addhelp(bqf_ncgp, "Input D, a proper discriminant.\n Returns the narrow class group, in the format [n,[n_1,...,n_r],[g_1,...,g_r]], where it has size n, is isomorphic to c_{n_1} x ... x c_{n_r} with n_1 | n_2 | ... | n_r, and g_i is a generator of the corresponding cyclic group of order n_i.");
		install("bqf_ncgp_lexic","Gp","bqf_ncgp_lexic","./libqquadratic.so");
		addhelp(bqf_ncgp_lexic, "Input D, a proper discriminant D.\n This returns [n,[n_1,...,n_r],[f1,f2,...,fl]], where n is the narrow class number of D, the narrow class group is c_{n_1} x ... x c_{n_r} with n_1 | n_2 | ... | n_r, and representative BQFs are the f1,f2,... written in lexicographic order: starting with the identity element, and the component with the highest order moves first.");
		install("bqf_pow_tc","GGD1,L,p","bqf_pow","./libqquadratic.so");
		addhelp(bqf_pow,"Inputs q, n, {tored=1}: BQF q, integer n, tored=0, 1.\n Returns q^n, reduced if tored=1.");
		install("bqf_square_tc","GD1,L,p","bqf_square","./libqquadratic.so");
		addhelp(bqf_square,"Inputs q, {tored=1}: BQF q, tored=0, 1.\n Returns q^2, reduced if tored=1.");

	\\REPRESENTATIONS OF NUMBERS BY BQFs
		install("bqf_reps_tc","GGD0,L,D1,L,p","bqf_reps","./libqquadratic.so");
		addhelp(bqf_reps,"Inputs q, n, {proper=0}, {half=1}: BQF q, integer n, (proper=0,1), (half=0,1).\n This solves the equation q(x,y)=n over the integers. We will get a finite set of (families) of solutions, and if half=0 we only return one of the families corresponding to (x,y) and (-x,-y). If we want only coprime solutions when disc!=square, pass in proper=1. If Q has discriminant D, the return is:\n\n --------If no solutions, returns 0;\n -D>0 nonsquare and n!=0, [[1,M],[x1,y1],[x2,y2],...,[xk,yk]] where the +/-(xi,yi) are representatives of the distinct classes of solutions and M is the invariant automorph;\n ----D>0 square and n!=0, [[0],[x1,y1],[x2,y2],...,[xk,yk]] where the +/-(xi,yi) are the solutions;\n --------------------D<0, [[0],[x1,y1],[x2,y2],...,[xk,yk]] where the +/-(xi,yi) are the solutions;\n -----------D=0 and n!=0, [[2],[[s1,s2],[x1,y1]]] where the solutions are (up to +/-) x=x1+Us1, y=y1+Us2;\n ---------n=0, D!=square, [[0],[0,0]];\n -------n=0, D=square!=0, [[2],[[s1,s2],[0,0]],[[s3,s4],[0,0]]], solutions are (x,y)=(s1k,s2k),(s3k,s4k) for integer k;\n ------------n=0 and D=0, if Q!=0 as above, if Q=0 then [[-1]] (everything is a solution).\n\n In general, -1=all, 0=finite, 1=positive, 2=linear");

	\\MORE REPRESENTATION OF NUMBERS
		install("bqf_bigreps_tc","GGp","bqf_bigreps","./libqquadratic.so");
		addhelp(bqf_bigreps,"Inputs: Q, n, Q=[A,B,C,D,E]=Ax^2+Bxy+Cy^2+Dx+Ey and an integer n.\n Returns the solutions to Q(x,y)=n over the integers. If D=bqf_disc(Q)=B^2-4AC, then the output is: \n\n ----------If no solutions, returns 0; \n ------------D>0 nonsquare, [[1,M,[s1,s2]],[x1,y1],[x2,y2],...,[xk,yk]] where the general solution is [x;y]=M^k*[xi;yi]+[s1;s2];\n --------D>0 square EITHER, [[0],[x1,y1],[x2,y2],...,[xk,yk]] where the (xi,yi) are the solutions;\n -----------------------OR, [[2],[[s1,t1],[x1,y1]],...,[[sk,tk],[xk,yk]]] where the solutions are x=xi+Usi and y=yi+Uti for U integer;\n ----------------------D<0, [[0],[x1,y1],[x2,y2],...,[xk,yk]] where the (xi,yi) are the solutions;\n ----------------------Q=0, [[-1]] if n=0 and 0 else;\n -D=0 and A=B=C=0 or D=E=0, [[2],[[s1,t1],[x1,y1]],...,[[sk,tk],[xk,yk]]] where the solutions are x=xi+Usi and y=yi+Uti for U integer (si=sj and ti=tj for all i,j in fact);\n ------------D=0 otherwise, [[-2],[[a1,b1,c1],[e1,f1,g1]],...,[[ak,bk,ck],[ek,fk,gk]]] where the solutions are x=ai*U^2+bi*U+ci, y=ei*U^2+fi*U+gi for U integer.\n\n In general, -2=quadratic, -1=all, 0=finite, 1=positive, 2=linear");
		install("bqf_linearsolve_tc","GGGGp","bqf_linearsolve","./libqquadratic.so");
		addhelp(bqf_linearsolve,"Inputs qf, n1, lin, n2: qf a six term integer vector representing the form Ax^2+By^2+Cz^2+Dxy+Exz+Fyz, n1 an integer, lin a three term integer vector representing Ax+By+Cz, and n2 an integer.\n Solves qf(x, y, z)=n1 and lin(x, y, z)=n2 simultaneously. If there are no solutions, this returns 0. Otherwise it returns a vector v. Let v[1][1]=t, and then the format of v is:\n\n -t=-2, v=[[-2], [[a1, a2, a3], [b1, b2, b3], [c1, c2, c3]], ...], where each general solution is x=a1U^2+a2U+a3, y=b1U^2+b2U+b3, z=c1U^2+c2U+c3 for any U integral;\n -t=-1, v=[[-1], [[a1, a2, a3], [b1, b2, b3], [c1, c2, c3]]], where the solution is x=a1U+b1V+c1, y=a2U+b2V+c2, z=a3U+v3V+c3 for any U, V integral;\n --t=0, v=[[0], [a1, b1, c1], ...], where the finite set of solutions are (x,y,z)=(ai, bi, ci);\n --t=1, v=[[1, M, [s1, s2, s3]], [a1, b1, c1], ...], where the general solution is [x;y;z]=M^k[ai;bi;ci]+[s1;s2;s3] for k integral. Note that ai and si need not be integral, though M is.\n --t=2, v=[[2], [[a1, a2, a3], [b1, b2, b3]], ...], where each general solution is x=a1U+b1, y=a2U+b2, z=a3U+b3 for any U integral;\n\n In general, -2=quadratic, -1=plane, 0=finite, 1=positive, 2=linear.");

	\\GENERAL HELP
		addhelp(bqf, "This package deals with binary quadratic forms with integer coefficients. A homogeneous binary quadratic form Ax^2+Bxy+Cy^2 is stored as [A,B,C]. A proper discriminant is an integer that is equivalent to 0 or 1 modulo 4 and is not a square. \n Subtopics:\n Discriminants (disc)\n Basic operations (bqfbasic)\n Indefinite forms (ibqf)\n Class group and composition (bqfclass)\n Representation of numbers (bqfsolve)");
		addhelp(disc,"disclist, discprimeindex, isdisc, pell, posreg, quadroot.");
		addhelp(bqfbasic,"bqf_automorph, bqf_disc, bqf_isequiv, bqf_isreduced, bqf_random, bqf_random_D, bqf_red, bqf_roots, bqf_trans, bqf_trans_coprime, ideal_tobqf.");
		addhelp(ibqf,"ibqf_isrecip, ibqf_leftnbr, ibqf_redorbit, ibqf_rightnbr, ibqf_river, ibqf_riverforms, ibqf_symmetricarc, mat_toibqf.");
		addhelp(bqfclass,"bqf_comp, bqf_ncgp, bqf_ncgp_lexic, bqf_pow, bqf_square.");
		addhelp(bqfsolve,"bqf_bigreps, bqf_linearsolve, bqf_reps.");


\\qq_bqf_int

	\\INTERSECTION DATA
		install("bqf_bdelta_tc","GG","bqf_bdelta","./libqquadratic.so");
		addhelp(bqf_bdelta, "Inputs q1,q2, integral binary quadratic forms.\n Returns B_{Delta}(q1,q2)=B1B2-2A1C2-2A2C1, where qi=[Ai,Bi,Ci].");
		install("bqf_intlevel_tc","GG","bqf_intlevel","./libqquadratic.so");
		addhelp(bqf_intlevel, "Inputs q1, q2, integral binary quadratic forms.\n Returns the signed intersection level of q1, q2, i.e. if qi=[Ai, Bi, Ci], then this is the gcd of [-A1B2 + A2B1, -2A1C2 + 2A2C1, -B1C2 + B2C1] times the sign of -A1B2+A2B1.");
		install("ibqf_intpoint_tc", "GGD0,G,","ibqf_intpoint","./libqquadratic.so");
		addhelp(ibqf_intpoint, "Inputs q1, q2, {location=0}: PIBQFs q1, q2 whose root geodesics intersect, complex number location, automorph of q1.\n Returns the upper half plane intersection point of q1, q2. If location=0, it is the intersection of q1, q2; if location=1, we translate it to the fundamental domain of PSL(2,Z); if imag(location)!=0, then location is assumed to be a point on \ell_q1. We translate the intersection point to the geodesic between location and gamma_q1(location). If the invariant automorph is large, then we need to increase the precision to ensure accurate results.");

	\\INTERSECTION NUMBER COMPUTATION
		install("ibqf_int_tc","GGp","ibqf_int","./libqquadratic.so");
		addhelp(ibqf_int, "Inputs q1, q2 PIBQFs.\n Returns the total intersection of q1 and q2, i.e. 2(Int_RS(q1,q2)+Int_RS(q1,-q2)).");
		install("ibqf_intRS_tc","GGp","ibqf_intRS","./libqquadratic.so");
		addhelp(ibqf_intRS, "Inputs q1, q2 PIBQFs.\n Returns the RS intersection of q1 and q2.");
		install("ibqf_intforms_tc","GGD0,L,D0,L,p","ibqf_intforms","./libqquadratic.so");
		addhelp(ibqf_intforms, "Inputs q1, q2, {data=0}, {split=0}: q1, q2 PIBQFs, data and split=0,1.\n Returns the intersecting forms of q1, q2 of all types. If data=1, each entry of the return is [B_delta(f1,f2), level of int, length of river overlap, f1, f2]; otherwise it is just the pair [f1, f2]. If split=0 returns a single vector of the return data, if split=1 it splits the return into [[RS], [RO], [LS], [LO]] intersection.");
		install("ibqf_intformsRS_tc","GGD0,L,p","ibqf_intformsRS","./libqquadratic.so");
		addhelp(ibqf_intformsRS,"Inputs q1, q2, {data=0}: q1, q2 PIBQFs, data=0,1.\n returns the RS intersection of q1 and q2 as a vector of non-simultaneously equivalent intersecting forms. If data=1, each return entry is instead [B_delta(f1,f2), level of int, length of river overlap, f1, f2]");
		install("ibqf_intformsRO_tc","GGD0,L,p","ibqf_intformsRO","./libqquadratic.so");
		addhelp(ibqf_intformsRO,"Inputs q1, q2, {data=0}: q1, q2 PIBQFs, data=0,1.\n Returns the RO intersection of q1 and q2 as a vector of non-simultaneously equivalent intersecting forms. If data=1, each return entry is instead [B_delta(f1,f2), level of int, length of river overlap, f1, f2]");
		install("ibqf_intformsLS_tc","GGD0,L,p","ibqf_intformsLS","./libqquadratic.so");
		addhelp(ibqf_intformsLS,"Inputs q1, q2, {data=0}: q1, q2 PIBQFs, data=0,1.\n Returns the LS intersection of q1 and q2 as a vector of non-simultaneously equivalent intersecting forms. If data=1, each return entry is instead [B_delta(f1,f2), level of int, length of river overlap, f1, f2]");
		install("ibqf_intformsLO_tc","GGD0,L,p","ibqf_intformsLO","./libqquadratic.so");
		addhelp(ibqf_intformsLO,"Inputs q1, q2, {data=0}: q1, q2 PIBQFs, data=0,1.\n Returns the LO intersection of q1 and q2 as a vector of non-simultaneously equivalent intersecting forms. If data=1, each return entry is instead [B_delta(f1,f2), level of int, length of river overlap, f1, f2]");

	\\GENERAL HELP
		addhelp(bqf_int, "This package deals with intersections of binary quadratic forms with integer coefficients. \n Subtopics:\n Computing intersection numbers (intcomp)\n Properties of an intersection (intprop)");
		addhelp(intcomp,"ibqf_int, ibqf_intRS, ibqf_intforms, ibqf_intformsRS, ibqf_intformsRO, ibqf_intformsLS, ibqf_intformsLO.");
		addhelp(intprop,"bqf_bdelta, bqf_intlevel, ibqf_intpoint.");


\\qq_geometry

	\\BASIC LINE, CIRCLE, AND POINT OPERATIONS
		install("arc_init_tc","GGGD0,L,p","arc_init","./libqquadratic.so");
		addhelp(arc_init,"Inputs c, p1, p2, {dir=0}: circle/arc c, p1 and p2 points on the circle, dir=-1, 0, 1.\n Returns the arc of c going counterclockwise from p1 to p2, where dir=1 means that is also the orientation, and dir=-1 means the orientation is backwards (and dir=0 means unoriented).");
		install("arc_midpoint_tc","GGGp","arc_midpoint","./libqquadratic.so");
		addhelp(arc_midpoint,"Inputs c, p1, p2: circle c, points p1, p2 on c.\n Returns the midpoint of the counterclockwise arc on c from p1 to p2.");
		install("circle_angle_tc","GGGp","circle_angle","./libqquadratic.so");
		addhelp(circle_angle,"Inputs c1, c2, p: circles/arcs c1, c2 that intersect in the point p.\n Returns the angle of intersection between c1 and c2 at p, with the angle formed by rotating the tangent to c1 at p counterclockwise to the tangent to c2 at p.");
		install("circle_fromcp","GGp","circle_fromcp","./libqquadratic.so");
		addhelp(circle_fromcp,"Inputs cent, p: centre cent, and point p.\n Returns the circle with centre cent and passing through point p.");
		install("circle_fromppp_tc","GGGp","circle_fromppp","./libqquadratic.so");
		addhelp(circle_fromppp,"Inputs p1, p2, p3: distinct complex points (oo is allowed).\n Returns the circle that passes through p1, p2, p3. If they are collinear (including if one is oo), this will return the line going through them instead.");
		install("circle_tangentslope_tc","GGp","circle_tangentslope","./libqquadratic.so");
		addhelp(circle_tangentslope,"Inputs c, p: circle c and point p on the circle.\n Returns the slope of the tangent to the circle at p.");
		install("crossratio","GGGG","crossratio","./libqquadratic.so");
		addhelp(crossratio, "Inputs a,b,c,d complex numbers or oo with at most one being oo.\n Returns the crossratio of [a,b;c,d].");
		install("line_angle","GGp","line_angle","./libqquadratic.so");
		addhelp(line_angle,"Inputs l1, l2: lines l1, l2.\n Returns the angle formed between l1 and l2, formed by rotating l1 counterclockwise to l2 (in the range [0, Pi)).");
		install("line_fromsp","GG","line_fromsp","./libqquadratic.so");
		addhelp(line_fromsp,"Inputs s, p: slope and point p.\n Returns the line through p with slope s.");
		install("line_frompp","GG","line_frompp","./libqquadratic.so");
		addhelp(line_frompp,"Inputs p1, p2: distinct complex points.\n Returns the line formed by p1 and p2.");
		install("mat_eval_tc","GG","mat_eval","./libqquadratic.so");
		addhelp(mat_eval, "Inputs M,x; M a matrix, and x number.\n Returns Mx with M acting via Mobius transformation. x=+/-oo is allowed.");
		install("midpoint","GG","midpoint","./libqquadratic.so");
		addhelp(midpoint,"Inputs: p1, p2 complex points.\n Returns the midpoint of the line segment from p1 to p2.");
		install("mobius_tc","GGp","mobius","./libqquadratic.so");
		addhelp(mobius,"Inputs M, c: matrix M, c a line/circle/segment/arc.\n Returns Mc, the mobius map acting on c.");
		install("perpbis","GG","perpbis","./libqquadratic.so");
		addhelp(perpbis,"Input p1, p2, two distinct complex numbers.\n Returns the perpendicular bisector.");
		install("radialangle_tc","GGp","radialangle","./libqquadratic.so");
		addhelp(radialangle,"Inputs c, p: c circle/circle arc, p point on c.\n Returns the angle formed between the horizontal radius facing right to p, and in the range [0,2*Pi).");
		install("slope","GG","slope","./libqquadratic.so");
		addhelp(slope,"Inputs p1, p2: distinct complex points.\n Returns the slope of the line through p1,p2");
	
	\\INTERSECTION OF LINES/CIRCLES
		install("arc_int_tc","GGp","arc_int","./libqquadratic.so");
		addhelp(arc_int,"Inputs c1, c2: circle arcs of circles that are not concentric.\n Returns the intersection points of c1, c2 (of size 0,1,2).");
		install("arcseg_int_tc","GGp","arcseg_int","./libqquadratic.so");
		addhelp(arcseg_int,"Inputs c, l: arc c and segment l.\n Returns the intersection points of c with l.");
		install("circle_int_tc","GGp","circle_int","./libqquadratic.so");
		addhelp(circle_int,"Input c1, c2: circles c1, c2 that are not concentric.\n Returns the set of points in their intersection.");
		install("circleline_int_tc","GGp","circleline_int","./libqquadratic.so");
		addhelp(circleline_int,"Inputs c, l: circle and line.\n Returns the intersection points of c and l.");
		install("genseg_int_tc","GGp","genseg_int","./libqquadratic.so");
		addhelp(genseg_int,"Inputs s1, s2: circles/arcs/lines/segmentss.\n Returns the intersection points of s1 and s2.");
		install("line_int_tc","GGp","line_int","./libqquadratic.so");
		addhelp(line_int,"Inputs l1, l2: lines.\n Returns the intersection point of the two lines. If the lines are parallel, returns oo.");
		install("onarc_tc","lGGp","onarc","./libqquadratic.so");
		addhelp(onarc,"Inputs c, p: circle arc c, point p assumed to be on the circle defined by c.\n Returns 1 if p is on the arc c (running counterclockwise from c[3] to c[4]). A circle can be inputted, and the function will trivially return 1.");
		install("onseg_tc","lGGp","onseg","./libqquadratic.so");
		addhelp(onseg,"Inputs l, p: line segment and point on it.\n Returns 1 if p is on the line segment l and 0 else.");
		install("seg_int_tc","GGp","seg_int","./libqquadratic.so");
		addhelp(seg_int,"Inputs l1, l2: line segmentss.\n Returns the intersection points of the line segments (length 0 or 1 vector).");

	\\DISTANCES
		install("hdist_tc","GGp","hdist","./libqquadratic.so");
		addhelp(hdist,"Inputs z1, z2 complex numbers in the upper half plane.\n Returns the hyperbolic distance between z1 and z2.");
		install("hdist_ud","GGp","hdist_ud","./libqquadratic.so");
		addhelp(hdist_ud,"Inputs z1, z2 complex numbers inside the unit disc.\n Returns the hyperbolic distance between z1 and z2 in the unit disc model.");
		install("hpolygon_area_tc","GGp","hpolygon_area","./libqquadratic.so");
		addhelp(hpolygon_area,"Input circles, vertices: the vectors of circles and vertices forming the edges of the polygon in the unit circle model (circles[i], circles[i+1] intersect at vertices[i]).\n Returns the area of the corresponding hyperbolic polygon. If there are edges on the unit circle (corresponding to circles[i]=0), the return is oo.");

	\\FUNDAMENTAL DOMAIN COMPUTATION
		install("edgepairing_tc","Gp","edgepairing","./libqquadratic.so");
		addhelp(edgepairing,"Input U, a normalized boundary.\n Returns [vp, vunp], where vp is a vector containing vecsmalls of paired edges, and vunp is a vector containing vecsmalls of [gind, i1ind] or [gind, i1ind, i2ind], where gind is the index of an unpaired side, and i1ind (and i2ind if there) are the unpaired vertices of the side.");
		install("normalizedboundary_outside_tc","lGGp","normalizedboundary_outside","./libqquadratic.so");
		addhelp(normalizedboundary_outside,"Inputs U, z: normalized boundary U, point z in the unit disc (not including the boundary).\n Returns -1 if z is inside U or on the boundary, and ind if z is outside U (i.e. between U and the edge of the unit disc), where ind is the index of the side between 0 and z.");
		install("randompoint_ud","Gp","randompoint_ud","./libqquadratic.so");
		addhelp(randompoint_ud,"Input R, positive real number.\n Returns a (uniform) random point in the unit disc of hyperbolic distance at most R from 0.");
		install("rootgeodesic_uhp_tc","Gp","rootgeodesic_uhp","./libqquadratic.so");
		addhelp(rootgeodesic_uhp,"Input M, a 2x2 hyperbolic matrix in PSL(2, R).\n Returns the upper half plane root geodesic of M.");

	\\PRINTING TO PLOTVIEWER
		install("python_printarcs","vGrD0,L,Drp","python_printarcs","./libqquadratic.so");
		addhelp(python_printarcs,"Inputs arcs, filename, {view=0}, {extrainput=NULL}: a set of arcs arcs, string filename, view=0, 1, and extrainput=NULL or a string.\n Prints the arcs specified by arcs to the file fdoms/filename.dat, ready for plotviewer. If view=1, calls plotviewer with the additional input of extrainput if you want to include other arcs/fundamental domains.");
		install("python_plotviewer","vr","python_plotviewer","./libqquadratic.so");
		addhelp(python_plotviewer,"Input S: string denoting the file names of data fundamental domains/geodesics.\n Launches the python file fdviewer.py to view the domain/geodesics. Enter the files separated by spaces (they must be stored in the sub-folder 'fdoms').");
		install("python_printfdom","vGrp","python_printfdom","./libqquadratic.so");
		addhelp(python_printfdom,"Input U, filename: fundamental domain U, string filename.\n Prints the fundamental domain U to the file fdoms/filename.dat, ready for the plot viewer. The filename must start with 'fd' to work properly.")

	\\HELPER METHODS
		install("atanoo","Gp","atanoo","./libqquadratic.so");
		addhelp(atanoo,"Inputs x, a real number or oo.\n Returns atan(x), with the convention of atan(oo)=Pi/2.");
		install("shiftangle","GGp","shiftangle","./libqquadratic.so");
		addhelp(shiftangle,"Inputs ang, bot: angles ang and bot.\n Returns ang+2*Pi*N for the unique integer N such that bot<=ang+2*Pi*N<bot+2*Pi.");
		
	\\GENERAL HELP
		addhelp(geometry,"This package deals with geometry in the complex plane. Points are complex numbers, with oo being the point at infinity. Circles/lines are represented by length 3 vectors, and circe arcs/line segments are length 8 vectors.\n Subtopics:\n geometryfunctions \n circles \n lines\n hyperbolic distances (hdist)\n fundamental domains (fd)");

		addhelp(geometryfunctions,"arc_init, arc_int, arc_midpoint, arcseg_int, atanoo, circle_angle, circle_fromcp, circle_fromppp, circle_int, circle_tangentslope, circleline_int, crossratio, genseg_int, hdist, hdist_ud, hpolygon_area, line_angle, line_fromsp, line_frompp, line_int, mat_eval, midpoint, mobius, onarc, onseg, perpbis, radialangle, seg_int, shiftangle, slope.");

		addhelp(circles,"A circle is stored as a vector [c, r, 0], where c is the centre, and r is the radius. The third entry is to distinguish a circle from a line.\n A circle arc is stored as a vector [c, r, p1, p2, p1ang, p2ang, dir, 0], where c,r are the centre and radius, p1 and p2 are the start and endpoints of the arc. We take the arc going counterclockwise from p1 to p2, and dir=1 means it is oriented from p1 to p2, and dir=-1 means it is oriented from p2 to p1. p1ang and p2ang are the radial angles to p1 and p2 from c, so that 0<p2-p1<2*Pi. The final 0 is to distinguish from a line segment.\n Relevant functions are:\n arc_init, arc_int, arc_midpoint, arcseg_int, circle_angle, circle_fromcp, circle_fromppp, circle_int, circle_tangentslope, circleline_int, genseg_int, mobius, onarc, radialangle.");

		addhelp(lines,"A line is stored as a vector [s, b, 1], where s is the slope, and b is the y-intercept if s!=oo, and the x-intercept if s=oo. The third entry is to distinguish a circle from a line.\n A line segment is stored as a vector [s, b, p1, p2, 0, ooendptor, dir, 1], where s, b are as before, and p1 and p2 are the start and endpoints of the segment. If neither endpoint is oo, then dir=1 means the segment from p1 to p2 in the upper half plane, and dir=-1 means the segment passing through oo. When one endpoint is oo, ooendptor=1 means the ray from p1 to p2 that is either pointing straight up or to the right, and ooendptor=-1 is pointing straight down or left. The final 1 is to distinguish from a circle arc.\n Relevant functions are:\n arcseg_int, circleline_int, genseg_int, line_angle, line_fromsp, line_frompp, line_int, midpoint, mobius, onseg, perpbis, seg_int, slope.");

		addhelp(hdist,"Hyperbolic distance/area functions: hdist, hdist_ud, hpolygon_area.");

		addhelp(fd,"Most of the fundamental domain methods are PARI-accessible only. The few available here are: edgepairing, normalizedboundary_outside, python_printarcs, python_plotviewer, python_printfdom, randompoint_ud, rootgeodesic_uhp.");

\\qq_quat


	\\BASIC OPERATIONS ON ELEMENTS IN QUATERNION ALGEBRAS 
		install("qa_conj_tc","G","qa_conj","./libqquadratic.so");
		addhelp(qa_conj,"Input x, element of a quaternion algebra.\n Returns the conjugate of x.");
		install("qa_conjby_tc","GGG","qa_conjby","./libqquadratic.so");
		addhelp(qa_conjby,"Inputs Q, x, y: quaternion algebra Q, elements x, y with y invertible.\n Returns yxy^(-1).");
		install("qa_inv_tc","GG","qa_inv","./libqquadratic.so");
		addhelp(qa_inv,"Inputs Q, x; quaternion algebra Q, element x. Returns x^(-1), and produces an error if x is not invertible");
		install("qa_m2rembed_tc","GG","qa_m2rembed","./libqquadratic.so");
		addhelp(qa_m2rembed,"Inputs Q, x; indefinite quaternion algebra Q, element x.\n Returns image of x under standard embedding of Q into M(2,R) (requires a>0).");
		install("qa_minpoly_tc","GG","qa_minpoly","./libqquadratic.so");
		addhelp(qa_minpoly,"Inputs Q, x; Q a quaternion algebra, x in Q.\n Returns the minimal polynomial of x. The format is [1,b,c] for x^2+bx+c, and [1,b] for x+b");
		install("qa_mul_tc","GGG","qa_mul","./libqquadratic.so");
		addhelp(qa_mul,"Inputs Q, x, y; quaternion algebra Q, elements x, y.\n Returns x*y.");
		install("qa_mulvec_tc","GG","qa_mulvec","./libqquadratic.so");
		addhelp(qa_mulvec,"Inputs Q, L: quaternion algbera Q, vector L of elements.\n Returns L[1]*L[2]*...*L[n].");
		install("qa_mulvecindices_tc","GGG","qa_mulvecindices","./libqquadratic.so");
		addhelp(qa_mulvecindices,"Inputs Q, L, indices: quaternion algebra Q, vector of elements of Q, vector (or vecsmall) of indices.\n Returns L[indices[1]]*L[indices[2]]*...*L[indices[n]].");
		install("qa_norm_tc", "GG", "qa_norm","./libqquadratic.so");
		addhelp(qa_norm,"Inputs Q, x; Q a quaternion algebra, x in Q.\n Returns the norm of x.");
		install("qa_pow_tc","GGG","qa_pow","./libqquadratic.so");
		addhelp(qa_pow,"Inputs Q, x, n; quaternion algebra Q, element x, integer n.\n Returns x^n.");
		install("qa_roots_tc","GGp","qa_roots","./libqquadratic.so");
		addhelp(qa_roots,"Inputs Q, x; indefinite quaternion algebra Q, element x with positive norm and hyperbolic (trace^2>4*norm).\n Returns the vector [rt1, rt2] of first root, second root of x under the standard embedding into SL(2,R).");
		install("qa_square_tc","GG","qa_square","./libqquadratic.so");
		addhelp(qa_square,"Inputs Q, x: quaternion algebra Q, element x.\n Returns x^2.");
		install("qa_trace_tc","G","qa_trace","./libqquadratic.so");
		addhelp(qa_trace,"Input x, element of a quaternion algebra.\n Returns the trace of x");

	\\BASIC OPERATIONS ON ORDERS/LATTICES IN QUATERNION ALGEBRAS 
		install("qa_isinorder_tc","lGGG","qa_isinorder","./libqquadratic.so");
		addhelp(qa_isinorder,"Inputs Q, ord, v; quaternion algebra Q, an order ord, and an element v.\n Returns 1 if v is in ord, and 0 if not. Can pass in an initialized order for ord.");
		install("qa_isorder_tc","lGG","qa_isorder","./libqquadratic.so");
		addhelp(qa_isorder,"Inputs Q, ord; Q a quaternion algebra, ord an 4x4 matrix whose columns span a supposed order in Q.\n Checks if ord is indeed an order, i.e. multiplication lands inside ord.");
		install("qa_leftorder_tc","GG","qa_leftorder","./libqquadratic.so");
		addhelp(qa_leftorder,"Inputs Q, L: quaternion algebra Q, lattice L.\n Returns the left order associated to L.");
		install("qa_rightorder_tc","GG","qa_rightorder","./libqquadratic.so");
		addhelp(qa_rightorder,"Inputs Q, L: quaternion algebra Q, lattice L.\n Returns the right order associated to L.");
		install("qa_ord_conj_tc","GGG","qa_ord_conj","./libqquadratic.so");
		addhelp(qa_ord_conj,"Inputs Q, ord, c; quaternion algebra Q, lattice ord, invertible element c.\n Returns the lattice c*ord*c^(-1) in Hermite normal form."); 
		install("qa_ord_disc_tc","GG","qa_ord_disc","./libqquadratic.so");
		addhelp(qa_ord_disc,"Inputs Q, ord; quaternion algebra Q, order ord.\n Ouputs the discriminant of the order ord.");
		install("qa_superorders_tc","GGG","qa_superorders","./libqquadratic.so");
		addhelp(qa_superorders,"Inputs Q, ord, n: quaternion algebra Q, (initialized) order ord, positive integer n.\n Returns the super orders of ord such that ord has index n in them.");
		
	\\INITIALIZATION METHODS
		install("qa_eichlerorder_tc","GGD0,G,","qa_eichlerorder","./libqquadratic.so");
		addhelp(qa_eichlerorder,"Inputs Q, l, {maxord=0}: quaternion algebra Q, positive integer l, maximal order of Q or 0 maxord.\n Returns an Eichler order of Q of level l, which is contained in maxord if set to being non-zero.");
		install("qa_ord_init_tc","GG","qa_ord_init","./libqquadratic.so");
		addhelp(qa_ord_init,"Inputs: Q, ord: a quaternion algebra and an order.\n Returns the initialization of the order: [ord, type, [d1, d2, d3, d4], level, prime factorization of the level]. ord is the hnf of the inputted ord, and d_n is the maximal denominator of the nth coefficient (i.e. 1, i, j, k). type=-1 means general order, type=0 means Maximal, type=1 means Eichler.");
		install("qa_init_ab_tc","GG","qa_init_ab","./libqquadratic.so");
		addhelp(qa_init_ab,"Inputs a, b: non-zero integers.\n Returns the quaternion algebra B=(a, b/Q).");
		install("qa_init_primes_tc","G","qa_init_primes","./libqquadratic.so");
		addhelp(qa_init_primes,"Input pset, a list of primes (or oo).\n Returns a quaternion algebra Q with ramification at the specified set of primes. If the list has odd length, then automatically adds oo to the ramification set.");
		install("qa_init_2primes_tc","GG","qa_init_2primes","./libqquadratic.so");
		addhelp(qa_init_2primes,"Inputs p, q distinct primes.\n Returns a quaternion algebra over Q with ramification at p,q only.");
		install("qa_maximalorder_tc","GD0,G,","qa_maximalorder","./libqquadratic.so");
		addhelp(qa_maximalorder,"Inputs Q, {baseord=0}: a quaternion algebra and an order.\n Returns a maximal order of Q, which contains baseord if baseord is non-zero.");
		install("qa_ram_fromab_tc","GG","qa_ram_fromab","./libqquadratic.so");
		addhelp(qa_ram_fromab,"Inputs a, b, non-zero integers.\n Returns the set of primes ramifying in the quaternion algebra (a, b / Q).");

	\\CONJUGATION OF ELEMENTS IN A GIVEN ORDER
		install("qa_conjbasis_tc","GGGGD0,L,","qa_conjbasis","./libqquadratic.so");
		addhelp(qa_conjbasis,"Inputs: Q, ord, e1, e2, {orient=0}; Q is the quaternion algebra, ord is an (initialized) order, e1 and e2 are the elements of Q, and orient=0,1.\n We test if e1 and e2 are conjugate in Q. If e1 or e2 is rational we automatically return 0. Else, if they are conjugate, the conjugation space is of rank 2. We return a Z basis for this space intersected with ord. If orient=0 we orient it correctly, otherwise we don't care.");
		install("qa_conjqf_tc","GGGG","qa_conjqf","./libqquadratic.so");
		addhelp(qa_conjqf,"Inputs: Q, ord, e1, e2; quaternion algebra Q, (initialized) order ord, non rational elements e1, e2 of Q with the same minimal polynomial.\n Returns [nrd(Xv_1+Yv_2), v1, v2] where Zv_1+Zv_2 is the Z-module for which ve1=e2v and v is in ord. We take it so v_2*conj(v_1)-v_1*conj(v_2) corresponds to the same embedding as e2.");
		install("qa_conjnorm_tc","GGGGGD0,L,p","qa_conjnorm","./libqquadratic.so");
		addhelp(qa_conjnorm,"Inputs: Q, ord, e1, e2, n, {retconj=0}; quaternion algebra Q, order ord, elements e1, e2 of Q, integer n, retconj=0,1.\n The return is 1 if there is an element v of ord of norm n conjugating e1 to e2 (v*e1*v^(-1)=e2), or 0 if no such element exists. If retconj=1 we also return a possible v.");
		install("qa_simulconj_tc","GGGGGGp","qa_simulconj","./libqquadratic.so");
		addhelp(qa_simulconj,"Inputs: Q, ord, e1, e2, f1, f2; quaternion algebra Q, (initialized) order ord, elements e1, e2, f1, f2 such that e1, f1 are conjugate, e2, f2 are conjugate, and trace(e1e2)=trace(f1f2). \n Then (e1 ,e2) and (f1, f2) are simultaneously conjugate, and the conjugation space has dimension 1. This returns a generator for this space intersected with ord (isomorphic to Z, so the generator is unique up to +/-).");

	\\EMBEDDING QUADRATIC ORDERS INTO EICHLER ORDERS
		install("qa_associatedemb_tc","GGGD0,G,","qa_associatedemb","./libqquadratic.so");
		addhelp(qa_associatedemb,"Inputs Q, order, emb,{D=0}: Q a quaternion algebra, order an initialized order, emb the image of image of (A+sqrt(D))/2, (D a discriminant which does not need to be passed).\n This method computes the pair [emb', D'] consisting of emb', the optimal embedding associated to emb and ord, and its corresponding discriminant D'.");
		install("qa_embed_tc","GGGD0,G,D0,L,p","qa_embed","./libqquadratic.so");
		addhelp(qa_embed,"Inputs Q, order, D, {nembs=0}, {retpell=0}; indefinite quaternion algebra Q, Eichler order order, discriminant D, nonnegative integer nembs, retpell=0,1.\n Returns the optimal embeddings of O_D into order up to conjugation by order_1^x. If retpell=1 we return the images of the fundamental unit in O_D (D>0 only), and otherwise we return the images of (D%2+sqrt(D))/2. Can pass the number of embeddings nembs if we have already counted them, or if we want a smaller number of embeddings (e.g. 1) can also pass that. If we pass a larger number than the number of embeddings the method will not end.");
		install("qa_embeddablediscs_tc","GGGGD0,L,D0,G,","qa_embeddablediscs","./libqquadratic.so");
		addhelp(qa_embeddablediscs,"Inputs Q, order, d1, d2, {fund=0}, {cop=0}: Q an indefinite quaternion algebra, initialized Eichler order order, d1<d2 integers, fund=0,1, cop a nonnegative integer.\n Returns the discriminants D with d1<=D<=d2 for which there exists an optimal embedding of discriminant D into Q[3]. If fund=1 only returns fundamental discriminants, and if cop!=0 only returns discriminants coprime to cop.");
		install("qa_numemb_tc","GGGD0,G,p","qa_numemb","./libqquadratic.so");
		addhelp(qa_numemb,"Inputs Q, order, D, {narclno=0}; Q an indefinite quaternion algebra, order an initialized Eichler order, D a discriminant, narclno the narrow class number of discriminant D.\n Returns [m, n, v1, v2, v3]. Here, m is the number of optimal embeddings of O_D into order, and n=h^+(D) is the number of a fixed orientation.\n v1=[x] with x the number of embeddings at oo (1 if D>0 and 2 if D<0).\n v2=[x_1, ..., x_r] with Q[1]=[p_1, ..., p_r] and x_i=number of local embeddings at p_i for i<=r (1 if p_i divides D and 2 else).\n v3=[y_1, ..., y_s], where s distinct primes (see order[5]) q_1,...,q_s divide the level of the order and y_i is the number of local embedding classes. This is 2 as long as q_i does not divide D.");
		install("qa_ordiffer_tc","GGGGD0,G,","qa_ordiffer","./libqquadratic.so");
		addhelp(qa_ordiffer,"Inputs Q, order, e1, e2, {D=0}: indefinite quaternion algebra Q, initialized Eichler order order, e1 and e2 optimal embeddings of O_D into order.\n Returns the set of finite primes for which the orientation of the embeddings at that prime differs.");
		install("qa_orinfinite_tc","GGD0,G,p","qa_orinfinite","./libqquadratic.so");
		addhelp(qa_orinfinite,"Inputs Q, emb, {D=0}: indefinite quaternion algebra Q, embedding emb of discriminant D.\n Returns the orientation of the embedding emb with respect to the standard embedding into M(2, R) (0 if D>0, and +/-1 if D<0).");
		install("qa_sortedembed_tc","GGGD0,L,D0,G,p","qa_sortedembed","./libqquadratic.so");
		addhelp(qa_sortedembed,"Inputs Q, order, D, {rpell=0}, {ncgp=0}; indefinite quaternion algebra Q, initialized Eichler order order, discriminant D, rpell=0,1, ncgp=0 or the output of bqf_ncgp_lexic(D).\n Returns the optimal embeddings of O_D into order sorted by embedding and group action. Rows correspond to a fixed orientation, with the first entry being the primes at which the orientation differs from the first line. The subsequent h^+(D) entries are the optimal embeddings, arranged so that lexicorder(D,1)[i] acts on form in position j takes it to the form in position i+j. Conjugating a row of embeddings by an Atkin Lehner involution takes the row to the row with the corresponding orientation.");
	
	\\ELEMENTS OF NORM N IN AN EICHLER ORDER
		install("qa_orbitreps_tc","GGGp","qa_orbitreps","./libqquadratic_testing.so");
		addhelp(qa_orbitreps,"Inputs Q, order, n: indefinite quaternion algebra Q, initialized Eichler order order, positive integer n.\n Returns a set of elements S, all of whom have reduced norm n, and are a complete set of representatives for O_{N=1}/O_{N=n}. If n is not coprime to the level of the order, returns 0.");
		install("qa_orbitrepsrange_tc","GGGp","qa_orbitrepsrange","./libqquadratic_testing.so");
		addhelp(qa_orbitrepsrange,"Inputs Q, order, n: indefinite quaternion algebra Q, initialized Eichler order order, positive integer n.\n Returns the vector of length n whose ith element is qa_orbitreps(Q, order, i).");
		install("qa_hecke_tc","GGGGp","qa_hecke","./libqquadratic_testing.so");
		addhelp(qa_hecke,"Inputs Q, order, n, emb: indefinite quaternion algebra Q, initialized Eichler order order, positive integer n, optimal embedding emb.\n Returns the action of T_n on emb. The format is a vector of [m, emb'], where emb' has multiplicity m.");

	\\FUNDAMENTAL DOMAIN METHODS
		install("qa_fundamentaldomain_tc","GGD0,G,D0,L,D0,G,p","qa_fundamentaldomain","./libqquadratic.so");
		addhelp(qa_fundamentaldomain,"Inputs Q, order, {p=0}, {dispprogress=0}, {ANRdata=0}: indefinite quaternion algebra Q, Eichler order order, {upper half plane point p (or 0)}, {dispprogress=0, 1}, {ANRdata=0 or length 5 vector}.\n Returns the fundamental domain of order, using Algorithm 4.8 of Voight with (probabilistic) enumeration via Algorithm 11 of Page. If p is passed in as 0, sets it to I/2 (which is never fixed if there is ramification). If dispprogress=1, returns partial results to the screen. The method relies on the constants [A, N, R, 1+nu, epsilon], which can be auto-set if ANRdata=0 or passed in by the user. Passing in a vector with some zero entries will automatically set those entries, and take the user input for the others (e.g. [10, 0, 0, 0, 1] sets A=10, epsilon=1, and auto-sets N, R, 1+nu).");
		install("qa_fdarea_tc","GGp","qa_fdarea","./libqquadratic.so");
		addhelp(qa_fdarea,"Inputs Q, order: indefinite quaternion algebra Q, initialized order order.\n Returns the hyperbolic area of Gamma\H, where Gamma is the image of the units of norm 1 in the order (a discrete subgroup of SL(2,R)). Only works for maximal/Eichler orders.");
		install("qa_isometriccircle_tc","GGGp","qa_isometriccircle","./libqquadratic.so");
		addhelp(qa_isometriccircle,"Inputs Q, x, p: indefinite quaternion algebra Q, norm 1 element x, upper half plane point p.\n Returns the isometric circle attached to x and p, i.e. mapping x to PSL(2, R), then to PSU(1, 1) (via the upper half plane map sending p ->0 and the upper half plane to the unit disc).");
		install("qa_normalizedbasis_tc","GGGp","qa_normalizedbasis","./libqquadratic.so");
		addhelp(qa_normalizedbasis,"Inputs Q, G, p: indefinite quaternion algebra Q, vector of norm 1 elements, upper half plane point p OR normalized boundary.\n Returns the normalized basis of G with respect to p; if p is a normalized boundary, we add G to it and compute the normalized basis. This is the normalized boundary of <G> (or <G union p[1]>), hence the fundamental domain if G (or <G union p[1]>) spans the group.");
		install("qa_normalizedboundary_tc","GGGp","qa_normalizedboundary","./libqquadratic.so");
		addhelp(qa_normalizedboundary,"Inputs Q, G, p: indefinite quaternion algebra Q, vector of norm 1 elements, upper half plane point p.\n Returns the normalized boundary of G, i.e. the union of the exteriors of the isometric circles of G. Returns [elements, icircs, vertices, vertex angles, matrices, area, 0]. The circle corresponding to elements[i] is icircs[i], and the vertices are vertices[i-1] and vertices[i]. angle[i] is the angle between elements[i] and elements[i+1], and matrices[i] is the image in PSU(1,1) of elements[i]. The element 1 corresponds to a section on the unit circle, which also corresponds to an angle and circle of -1. Vertex angles stores the angles to each of the vertices, the area is the area, and the side pairing stores the side pairing for fundamental domains (this function returns it as 0 however, as we assume that we are not a fundamental domain).");
		install("qa_printisometriccircles_tc","vGGGrD0,L,p","qa_printisometriccircles","./libqquadratic.so");
		addhelp(qa_printisometriccircles, "Inputs Q, L, p, filename, {view=0}: indefinite quaternion algebra Q, vector L of norm 1 elements, upper half plane point p, file name filename, {view=0, 1}.\n Prints the isometric circles corresponding to L and p to fdoms/filename.dat, ready for the plotviewer. If view=1, calls the plotviewer to view the geodesics.");
		install("qa_reduceelt_tc","GGGD0,G,D0,G,p","qa_reduceelt","./libqquadratic.so");
		addhelp(qa_reduceelt,"Inputs Q, G, g, {z=0}, {p=0}: indefinite quaternion algebra Q, set G of elements of norm 1 OR a normalized boundary, element g of norm 1 in Q, point z in unit disc, point p in upper half plane.\n Returns the triple [gammabar, delta, decomp], where gammabar=delta*g is (G,z)-reduced (i.e. distance between gammabar*z and 0 is less than or equal to the distance between g'*gammabar*z for all g' in G), and decomp is the vecsmall [i1, i2, ..., in] with delta=G[i1]*G[i2]*...*G[in]. If G is a normalized boundary, this method is a fair bit faster, and the value p does not need to be inputted (z can be inputted as non-zero if desired).");
		install("qa_rootgeodesic_fd_tc","GGGp","qa_rootgeodesic_fd","./libqquadratic.so");
		addhelp(qa_rootgeodesic_fd,"Inputs Q, U, g: indefinite quaternion algebra Q, fundamental domain of an Eichler order in Q U, non-rational element g.\n Returns the root geodesic of g in the fundamental domain. The format is [g's, circle arcs, vecsmall(sides hit), vecsmall(sides emenating from)].");
		install("qa_smallnorm1elts_tc","GGGGGD0,G,p","qa_smallnorm1elts","./libqquadratic.so");
		addhelp(qa_smallnorm1elts,"Inputs Q, order, p, z, C1, {C2=0}: indefinite quaternion algebra Q, initialized Eichler order order, upper half plane point p, point in the unit disc z, non-negative real numbers C1, C2.\n Finds ''small'' norm 1 elements, i.e. (if z=0) the set G(C2)-G(C1) from ''Computing fundamental domains'' by Voight. If C2=0, returns G(C1). If z!=0, shifts the basepoint (ala Page) to z, so elements with isometric circles near z are weighted more.");
		install("qa_topsu_tc","GGGp","qa_topsu","./libqquadratic.so");
		addhelp(qa_topsu,"Inputs Q, g, p: indefinite quaternion algebra Q, norm 1 element g, upper half plane point p.\n This returns the image of g in PSU(1, 1), given by mapping g to PSL(2, R) and then to acting on the unit disc via p->0.");

	\\SUPPORTING METHODS
		install("module_intersect_tc","GG","module_intersect","./libqquadratic.so");
		addhelp(module_intersect,"Inputs A, B: Z-modules spanned by the columns.\n Returns the intersection of the Z-modules, as a Z-module spanned by the colunmns.");
		install("prime_ksearch_tc","GD0,G,","prime_ksearch","./libqquadratic.so");
		addhelp(prime_ksearch, "Inputs rels, {extra=0}; rels=[[p_1,s_1],...,[p_k,s_k]], with the p_i distinct integers and s_i=-1 or 1 for each i, and extra=0 or [n, c].\n Returns the smallest prime r such that kronecker(r,p_i)=s_i for all i, and r== c_i mod n_i for all i. If the inputs are not consistent, the program will issue an error before the memory has been exhausted.");
		install("QM_hnf_tc","G","QM_hnf","./libqquadratic.so");
		addhelp(QM_hnf,"Input M a rational matrix.\n Returns the Hermite normal form of M with resepct to columns.")
		install("powerset_tc","G","powerset","./libqquadratic.so");
		addhelp(powerset,"Input L, a vector.\n Returns all subsets of L.");
		install("vecratio_tc","GG","vecratio","./libqquadratic.so");
		addhelp(vecratio,"Input v1, v2: vectors in the same 1-dim linear subspace.\n This returns v1/v2. If v1=0, returns 0, and if v2=0, returns oo.");

	\\GENERAL HELP
		addhelp(quat, "This package deals with quaternion algebras.\n Subtopics:\n Basic operations on elements (qabasic)\n Basic operations on orders (qaordbasic)\n Initialization methods (qainit)\n Conjugation of elements in a given order (qaconj)\n Embedding quadratic orders (qaemb)\n Elements of norm n (qanormn)\n Fundamental domain (qafdom)\n Supporting methods (qasupport)");
		addhelp(qabasic, "qa_conj, qa_conjby, qa_inv, qa_m2rembed, qa_minpoly, qa_mul, qa_mulvec, qa_mulvecindices, qa_norm, qa_pow, qa_roots, qa_square, qa_trace.");
		addhelp(qaordbasic, "qa_isinorder, qa_isorder, qa_leftorder, qa_rightorder, qa_ord_conj, qa_ord_disc, qa_superorders.");
		addhelp(qainit, "qa_eichlerorder, qa_ord_init, qa_init_ab, qa_init_primes, qa_init_2primes, qa_maximalorder, qa_ram_fromab");
		addhelp(qaconj, "qa_conjbasis, qa_conjqf, qa_conjnorm, qa_simulconj")
		addhelp(qaemb, "qa_associatedemb, qa_embed, qa_embeddablediscs, qa_numemb, qa_ordiffer, qa_orinfinite, qa_sortedembed.");
		addhelp(qanormn, "qa_orbitreps, qa_orbitrepsrange, qa_hecke");
		addhelp(qafdom, "qa_fundamentaldomain, qa_fdarea, qa_isometriccircle, qa_normalizedbasis, qa_normalizedboundary, qa_printisometriccircles, qa_reduceelt, qa_rootgeodesic_fd, qa_smallnorm1elts, qa_topsu.");
		addhelp(qasupport, "module_intersect, powerset, prime_ksearch, QM_hnf, vecratio.");

\\qq_quat_int

	\\INTERSECTION NUMBER BASED OFF OF ROOTS BOUND AREA
		install("qa_inum_roots_tc","GGGGD1,L,p","qa_inum_roots","./libqquadratic.so");
		addhelp(qa_inum_roots,"Inputs Q, order, e1, e2, {data=1}: indefinite quaternion algebra Q, initialized Eichler order order, e1, e2 not rational and lying in real quadratic subfields, data=0, 1.\n This calculates the intersections of e1 and e2 with respect to order using the root geodesic method. If data=0, only returns the pairs, otherwise returns [[intersections],[data]], where intersections[i] is a pair of optimal embeddings that intersect, and data[i] is the pair [signed level,x] corresponding to the embedding pair.");

	\\INTERSECTIONS BASED ON X-LINKAGE
		install("qa_inum_x_tc","GGGGD1,L,p","qa_inum_x","./libqquadratic.so");
		addhelp(qa_inum_x,"Inputs Q, order, e1, e2, {data=1}: indefinite quaternion algebra Q, initialized Eichler order order, e1, e2 optimal embeddings into order, data=0, 1.\n This calculates the intersections number of e1 and e2 with respect to order, via x-linking. If data=0, only returns the pairs, otherwise returns [[intersections],[data]], where intersections[i] is a pair of optimal embeddings that intersect, and data[i] is the pair [signed level,x] corresponding to the embedding pair.");
		install("qa_xlink_tc","GGGGGp","qa_xlink","./libqquadratic.so");
		addhelp(qa_xlink,"Inputs Q, ord, e1, e2, x: indefinite quaternion algebra Q, initialized Eichler order ord, optimal embeddings e1, e2, integer x.\n Returns the set of pairs of x-linked optimal embeddings conjugate to (e1, e2) taken up to simultaneous conjugation.");
		install("qa_xposs_tc","GGGD0,G,D0,G,","qa_xposs","./libqquadratic.so");
		addhelp(qa_xposs,"Inputs P, D1, D2, {xmin=0}, {xmax=0}: even sized set of primes OR quaternion algebra P, discriminants D1, D2, 0<=xmin<=xmax.\n Returns the set of x's between xmin and xmax (inclusive) which exhbit x-linking in the indefinite quaternion algebra ramified at P (or P itself when it is a quaternion algebra). If xmin=xmax=0, returns x's between 0 and sqrt(D_1D_2).");
	
	\\INTERSECTIONS BASED ON FUNDAMENTAL DOMAIN
		install("qa_inum_fd_tc","GGGGGD1,L,p","qa_inum_fd","./libqquadratic.so");
		addhelp(qa_inum_fd,"Inputs Q, order, U, e1, e2, {data=1}: indefinite quaternion algebra Q with order order, fundamental domain of order U, embeddings e1 and e2, data=0, 1.\n Returns the intersection number of e1 and e2 based on computing their root geodesics in the fundamental domain. If data=0, only returns the pairs, otherwise returns [[intersections],[data]], where intersections[i] is a pair of optimal embeddings that intersect, and data[i] is the pair [signed level,x] corresponding to the embedding pair.");
	
	\\INTERSECTION SERIES
		install("qa_inumseries_tc","GGGGGGD1,L,p","qa_inumseries","./libqquadratic_testing.so");
		addhelp(qa_inumseries,"Inputs Q, order, U, e1, e2, N, {type=1}: indefinite quaternion algebra Q, Eichler order order, fundamental domain U, optimal embeddings e1, e2, positive integer N, type=0, 1, prime.\n Returns the intersection series associated to e1, e2 and type, i.e. sum_{n=1}^N Int_{type}(e1, T_n e2)q^n, where type=0 means unsigned, =1 means signed, and >1 means type-weighted, where type should be a prime ramifying in Q or dividing the level of the order order (as the intersection pairing is Hecke-adjoint for those primes).");

	\\INTERSECTION DATA
		install("qa_intlevel_tc","GGGGD0,G,D0,G,p","qa_intlevel","./libqquadratic.so");
		addhelp(qa_intlevel,"Inputs Q, order, e1, e2, {D1=0}, {D2=0}: indefinite quaternion algebra Q, initialized Eichler order order, embeddings e1, e2 of discriminants D1, D2 (which do not need to be passed in).\n Returns the pair [signed level,x] corresponding to the embedding pair.");
	
	\\GENERAL HELP
		addhelp(quat_int, "This package deals with intersection numbers of optimal embeddings.\n Methods:\n qa_intlevel, qa_inum_fd, qa_inum_roots, qa_inum_x, qa_xlink, qa_xposs.");

\\qq_visual

	\\HISTOGRAMS
		install("hist_make","GrrD0,L,DrD0,L,p","hist_make","./libqquadratic.so");
		addhelp(hist_make,"Inputs data, imagename, filename, {compilenew=0}, {plotoptions=NULL}, {open=0}: sorted list of real numbers data, name of the tikz picture, name of the LaTeX file (without .tex), {compilenew=0, 1}, {plotoptions=NULL, string}, {open=0, 1}.\n Automatically bins the data, and creates a pdf of the histogram using tikz and externalize. The output is in the folder /images, with the build file (named filename) being in the folder /images/build. If compilenew=0 assumes the LaTeX document to compile the plot is pre-made, and otherwise this method automatically writes it. If additionally, plotoptions!=NULL, this string is inserted in between \\begin{axis} and \\end{axis} in the LaTeX document (allowing one to customize how the histogram looks). If open=1, the pdf is automatically opened (only works with Linux subsystem for Windows). The returned value is used to modify the histogram, e.g. changing the bins, scaling it, and changing the range.");
		install("hist_rebin","GGGp","hist_rebin","./libqquadratic.so");
		addhelp(hist_rebin,"Inputs data, histdata, nbins: the sorted data, the length 8 vector output of a hist_ method, the number of bins.\n Rebins the data according to the new number of bins, and updates histdata.");
		install("hist_recompile","vG","hist_recompile","./libqquadratic.so");
		addhelp(hist_recompile,"Input histdata, the length 8 vector output of a hist_ method.\n Recompiles the histogram, returning nothing. This is used when you edit the LaTeX document by hand.");
		install("hist_rerange","GGGGp","hist_rerange","./libqquadratic.so");
		addhelp(hist_rerange,"Inputs data, histdata, minx, maxx: the sorted data, the length 8 vector output of a hist_ method, minimum x-value, maximum x-value.\n Rebins the data according to the new minimum and maximum value. This is useful when there are outliers that skew the look of the graph. Returns the updated histdata.");
		install("hist_rescale","GGLp","hist_rescale","./libqquadratic.so");
		addhelp(hist_rescale,"Inputs data, histdata, scale: the sorted data, the length 8 vector output of a hist_ method, and scale=0, 1.\n If scale=1 scales the data so the total area is 1, and if scale=0 uses the absolute count for the y-axis. Returns the updated histdata.");

	\\GENERAL HELP
		addhelp(visual,"This package deals with visualizing data. Subtopics:\n Histograms (hist)");
		
		addhelp(hist,"hist_make, hist_rebin, hist_recompile, hist_rerange, hist_rescale.");

