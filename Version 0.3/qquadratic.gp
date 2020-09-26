\\qq_base

	\\COMPLEX GEOMETRY
		install("crossratio","GGGG","crossratio","./libqquadratic.so");
		addhelp(crossratio, "Inputs a,b,c,d complex numbers or oo with at most one being oo.\\n Outputs the crossratio of [a,b;c,d].");
		install("mat_eval_typecheck","GG","mat_eval","./libqquadratic.so");
		addhelp(mat_eval, "Inputs M,x; M a matrix, and x number.\\n Outputs Mx with M acting via Mobius transformation. x=+/-oo is allowed.");

	\\INFINITY
		install("addoo","GG","addoo","./libqquadratic.so");
		addhelp(addoo, "Inputs a,b real numbers or oo.\\n Outputs a+b, where if a or b is +/- oo, returns that back. Note that this will make oo+-oo=oo and -oo+oo=-oo");
		install("divoo","GG","divoo","./libqquadratic.so");
		addhelp(divoo, "Inputs a,b, real numbers or oo.\\n Outputs a/b, where oo is output if a=oo and b>=0 or a=-oo and b<0 or b=0 and a>=0 (outputs -oo under analogous assumptions).");

	\\LINEAR EQUATIONS AND MATRICES
		install("lin_intsolve_typecheck","GGG","lin_intsolve","./libqquadratic.so");
		addhelp(lin_intsolve, "Inputs A,B,n integers.\\n Outputs the general integral solutions to Ax+By=n. The format is [[s1,s2],[x0,y0]], where the general solution is x=s1*t+x0, y=s2*t+y0 for t an integer. The output is also reduced, i.e. gcd(s1,s2)=1. If A=B=0 or there are no integer solutions, returns 0.");
		install("mat3_complete_typecheck", "GGG", "mat3_complete", "./libqquadratic.so");
		addhelp(mat3_complete, "Inputs A,B,C integers with gcd 1.\\n Outputs a 3x3 integer matrix with determinant 1 whose top row is [A, B, C].");

	\\RANDOM
		install("rand_elt","G","rand_elt","./libqquadratic.so");
		addhelp(rand_elt, "Inputs v, a vector.\\n Outputs a random component of the vector.");

	\\SOLVING EQUATIONS MOD N
		install("sqmod_typecheck","GG","sqmod","./libqquadratic.so");
		addhelp(sqmod, "Inputs: x,n: rational number x, integer n with gcd(n,denom(x))=1.\\n Solves y^2==x mod n, and outputs the solutions. The output format is 0 if no solutions, and [S,m] otherwise, where the solutions are y==S[i] modulo m.");

	\\TIME
		install("printtime","v","printtime","./libqquadratic.so");
		addhelp(printtime, "Prints current time");

	\\GENERAL HELP
		addhelp(base,"This package is a collection of miscellaneous methods that may be useful in a variety of settings, and not just for the programs they were originally created for \\n Subtopics: \\n Complex plane geometry (cpg) \\n Infinity (inf) \\n Integer linear equations and matrices (ilm) \\n Modulo n methods (modn) \\n Random (rand)\\n Time (time)");
		addhelp(cpg,"crossratio, mat_eval.");
		addhelp(inf,"addoo, divoo.");
		addhelp(ilm,"lin_intsolve, mat3_complete.");
		addhelp(modn,"sqmod.");
		addhelp(rand,"rand_elt.");
		addhelp(time,"printtime.");


\\qq_bqf

	\\DISCRIMINANT METHODS
		install("disclist","GGD0,L,D0,G,","disclist","./libqquadratic.so");
		addhelp(disclist, "Inputs d1, d2, {fund=0}, {cop=0}: d1 and d2 integers with d1<=d2, fund=0,1, cop integer.\n Outputs the set of proper discriminants between d1 and d2, inclusive. If fund=1, only outputs fundamental discriminants. If cop!=0, only outputs discriminants coprime to cop.");
		install("discprimeindex_typecheck","G","discprimeindex","./libqquadratic.so");
		addhelp(discprimeindex, "Inputs: D, a proper discriminant.\n Outputs all prime divisors p of D for which D/p^2 is a proper discriminant.");
		install("fdisc_typecheck","G","fdisc","./libqquadratic.so");
		addhelp(fdisc, "Inputs: D, a proper discriminant.\n Outputs the fundamental discriminant associated to D, and 0 if D is not a proper discriminant.");
		install("isdisc","iG","isdisc","./libqquadratic.so");
		addhelp(isdisc, "Inputs: D a real number.\n Outputs 1 if D is a proper discriminant, and 0 otherwise.");
		install("pell_typecheck","G","pell","./libqquadratic.so");
		addhelp(pell, "Inputs: D a positive discriminant.\n Outputs [T, U], which is the smallest positive integer solution to T^2-DU^2=4 (and so (T+Usqrt(D))/2 is the fundamental unit in O_D).");
		install("posreg_typecheck","Gp","posreg","./libqquadratic.so");
		addhelp(posreg, "Inputs: D a positive discriminant.\n Outputs the positive regulator of O_D, i.e. the logarithm of the fundamental totally positive unit.");
		install("quadroot_typecheck","G","quadroot","./libqquadratic.so");
		addhelp(quadroot, "Input D, a non-square integer.\n Outputs sqrt(D) of type t_QUAD.");

	\\BASIC OPERATIONS ON BINARY QUADRATIC FORMS
		install("bqf_automorph_typecheck","G","bqf_automorph","./libqquadratic.so");
		addhelp(bqf_automorph, "Inputs: q a BQF.\n Outputs a generator of the automorph group of q in PSL(2,Z).");
		install("bqf_disc_typecheck","G","bqf_disc","./libqquadratic.so");
		addhelp(bqf_disc, "Inputs: q, quadratic form.\n Outputs the discriminant of q.");
		install("bqf_isequiv_typecheck","GGD0,L,p","bqf_isequiv","./libqquadratic.so");
		addhelp(bqf_isequiv, "Inputs: q, S, {tmat=0}: q a BQF, S either a BQF or a set of BQFs, tmat=0,1.\n This method tests if q is PSL(2,Z) equivalent to S or any form in S. If S is a form, this returns 1 if equivalent and 0 if not (if tmat!=0, returns a possible transition matrix).\n If S is a set of forms, this returns 0 if not equivalent and an index i such that q is equivalent to S[i] otherwise. If tmat!=0, this returns [index, transition matrix].");
		install("bqf_isreduced_typecheck","iG","bqf_isreduced","./libqquadratic.so");
		addhelp(bqf_isreduced,"Input q, quadratic form.\n Outputs 1 if q is reduced, and 0 if not.");
		install("bqf_random","GD0,L,D1,L,","bqf_random","./libqquadratic.so");
		addhelp(bqf_random,"Inputs maxc, {type=0}, {primitive=1}; maxc a positive integer, type=-1,0,1, and primitive=0,1.\n Outputs a random BQF with coefficients bounded by maxc. If type=-1 it is positive definite, =1 is indefinite, and =0 means either. If primitive=1 the form is primitive, else it doesn't have to be.");
		install("bqf_random_D","GG","bqf_random_D","./libqquadratic.so");
		addhelp(bqf_random_D,"Inputs maxc, D: maxc a positive integer, and D a discriminant.\n Outputs a random primitive form (positive definite if D<0) of discriminant D whose B coefficient is bounded by maxc.");
		install("bqf_red_typecheck","GD0,L,p","bqf_red","./libqquadratic.so");
		addhelp(bqf_red, "Inputs: q, {tmat=0}: BQF q, (tmat=0,1).\n Outputs a reduced form equivalent to q, and if tmat!=0, we output [q_red, transition matrix].");
		install("bqf_roots_typecheck","G","bqf_roots","./libqquadratic.so");
		addhelp(bqf_roots, "Inputs q: quadratic form q.\n Outputs the roots of q with the first root first.");
		install("bqf_trans_typecheck","GG","bqf_trans","./libqquadratic.so");
		addhelp(bqf_trans, "Inputs q, mtx: integral quadratic form q, 2x2 integral matrix mtx.\n Outputs the form acquired by replacing (x,y)^T with m(x,y)^T.");
		install("bqf_trans_coprime_typecheck", "GG", "bqf_trans_coprime", "./libqquadratic.so");
		addhelp(bqf_trans_coprime,"Inputs q, n: q a primitive integral BQF, and n an integer.\n Outputs a form similar to q whose first coefficient is coprime to n.");
		install("ideal_tobqf","GG","ideal_tobqf","./libqquadratic.so");
		addhelp(ideal_tobqf,"Inputs nf, ideal: a quadratic number field nf with ideal ideal.\n Outputs the corresponding binary quadratic form.");

	\\BASIC OPERATIONS SPECIFIC TO INDEFINITE FORMS
		install("ibqf_isrecip_typecheck","iGp","ibqf_isrecip","./libqquadratic.so");
		addhelp(ibqf_isrecip,"Inputs: q, a PIBQF.\n Outputs 1 if q is q is reciprocal, and 0 otherwise.");
		install("ibqf_leftnbr_typecheck","GD0,L,p","ibqf_leftnbr","./libqquadratic.so");
		addhelp(ibqf_leftnbr, "Inputs q, {tmat=0}: an indefinite binary quadratic form on the river and tmat=0,1. \n Outputs q' or [q',mat] (if tmat=0,1 respectively), where q' is the left neighbour of q and mat is the bqf_transition matrix from q to q'. The left neighbour is the previous form along the flow of the river that is reduced (AC<0 and B>|A+C|, and occurs when the branches swap from being below to above or vice versa).");
		install("ibqf_redorbit_typecheck","GD0,L,D0,L,p","ibqf_redorbit","./libqquadratic.so");
		addhelp(ibqf_redorbit, "Inputs q, (tmat), (posonly): q a PIBQF, tmat and posonly=0,1.\n Returns the reduced orbit of q. If tmat=1, also returns the corresponding transition matrices, and if posonly=1 only outputs the reduced forms with A>0.");
		install("ibqf_rightnbr_typecheck","GD0,L,p","ibqf_rightnbr","./libqquadratic.so");
		addhelp(ibqf_rightnbr, "Inputs q, (tmat): an indefinite binary quadratic form on the river and tmat=0,1. \n Outputs q' or [q',mat] (if tmat=0,1 respectively), where q' is the right neighbour of q and mat is the bqf_transition matrix from q to q'. The right neighbour is the next form along the flow of the river that is reduced (AC<0 and B>|A+C|, and occurs when the branches swap from being below to above or vice versa).");
		install("ibqf_river_typecheck","Gp","ibqf_river","./libqquadratic.so");
		addhelp(ibqf_river, "Input: q an indefinite quadratic form.\n Outputs the river sequence corresponding to q, where a 1 corresponds to going right and 0 corresponds to going left.");
		install("ibqf_riverforms_typecheck","Gp","ibqf_riverforms","./libqquadratic.so");
		addhelp(ibqf_riverforms, "Input q a PIBQF.\n This calculates all forms on the river of q, and outputs those with A>0, in the order that they appear on the river.");
		install("ibqf_symmetricarc_typecheck","Gp","ibqf_symmetricarc","./libqquadratic.so");
		addhelp(ibqf_symmetricarc,"Input q, a PIBQF.\n Outputs [z,gamma_q(z)] on the root geodesic corresponding to q so that q,gamma_q are symmetric about the arc.");
		install("mat_toibqf_typecheck","G","mat_toibqf","./libqquadratic.so");
		addhelp(mat_toibqf, "Inputs: mtx, a hyperbolic matrix in SL(2,Z).\n This outputs the PIBQF for which it is the ibqf_automorph, namely [c,d-a,-b]/gcd(c,d-a,b) if mtx=[a,b;c,d].");

	\\CLASS GROUPS AND COMPOSITION OF FORMS
		install("bqf_comp_typecheck","GGD1,L,p","bqf_comp","./libqquadratic.so");
		addhelp(bqf_comp,"Inputs q1, q2, {tored=1}: BQFs q1, q2 of the same discriminant, tored=0, 1.\n Outputs the composition of q1 and q2, reduced if tored=1.");
		install("bqf_ncgp","Gp","bqf_ncgp","./libqquadratic.so");
		addhelp(bqf_ncgp, "Input D, a proper discriminant.\n Outputs the narrow class group, in the format [n,[n_1,...,n_r],[g_1,...,g_r]], where it has size n, is isomorphic to c_{n_1} x ... x c_{n_r} with n_1 | n_2 | ... | n_r, and g_i is a generator of the corresponding cyclic group of order n_i.");
		install("bqf_ncgp_lexic","Gp","bqf_ncgp_lexic","./libqquadratic.so");
		addhelp(bqf_ncgp_lexic, "Input D, a proper discriminant D.\n This outputs [n,[n_1,...,n_r],[f1,f2,...,fl]], where n is the narrow class number of D, the narrow class group is c_{n_1} x ... x c_{n_r} with n_1 | n_2 | ... | n_r, and representative BQFs are the f1,f2,... written in lexicographic order: starting with the identity element, and the component with the highest order moves first.");
		install("bqf_pow_typecheck","GGD1,L,p","bqf_pow","./libqquadratic.so");
		addhelp(bqf_pow,"Inputs q, n, {tored=1}: BQF q, integer n, tored=0, 1.\n Outputs q^n, reduced if tored=1.");
		install("bqf_square_typecheck","GD1,L,p","bqf_square","./libqquadratic.so");
		addhelp(bqf_square,"Inputs q, {tored=1}: BQF q, tored=0, 1.\n Outputs q^2, reduced if tored=1.");

	\\REPRESENTATIONS OF NUMBERS BY BQFs
		install("bqf_reps_typecheck","GGD0,L,D1,L,p","bqf_reps","./libqquadratic.so");
		addhelp(bqf_reps,"Inputs q, n, {proper=0}, {half=1}: BQF q, integer n, (proper=0,1), (half=0,1).\n This solves the equation q(x,y)=n over the integers. We will get a finite set of (families) of solutions, and if half=0 we only output one of the families corresponding to (x,y) and (-x,-y). If we want only coprime solutions when disc!=square, pass in proper=1. If Q has discriminant D, the output is:\n\n --------If no solutions, returns 0;\n -D>0 nonsquare and n!=0, [[1,M],[x1,y1],[x2,y2],...,[xk,yk]] where the +/-(xi,yi) are representatives of the distinct classes of solutions and M is the invariant automorph;\n ----D>0 square and n!=0, [[0],[x1,y1],[x2,y2],...,[xk,yk]] where the +/-(xi,yi) are the solutions;\n --------------------D<0, [[0],[x1,y1],[x2,y2],...,[xk,yk]] where the +/-(xi,yi) are the solutions;\n -----------D=0 and n!=0, [[2],[[s1,s2],[x1,y1]]] where the solutions are (up to +/-) x=x1+Us1, y=y1+Us2;\n ---------n=0, D!=square, [[0],[0,0]];\n -------n=0, D=square!=0, [[2],[[s1,s2],[0,0]],[[s3,s4],[0,0]]], solutions are (x,y)=(s1k,s2k),(s3k,s4k) for integer k;\n ------------n=0 and D=0, if Q!=0 as above, if Q=0 then [[-1]] (everything is a solution).\n\n In general, -1=all, 0=finite, 1=positive, 2=linear");

	\\MORE REPRESENTATION OF NUMBERS
		install("bqf_bigreps_typecheck","GGp","bqf_bigreps","./libqquadratic.so");
		addhelp(bqf_bigreps,"Inputs: Q, n, Q=[A,B,C,D,E]=Ax^2+Bxy+Cy^2+Dx+Ey and an integer n.\n Outputs the solutions to Q(x,y)=n over the integers. If D=bqf_disc(Q)=B^2-4AC, then the output is: \n\n ----------If no solutions, returns 0; \n ------------D>0 nonsquare, [[1,M,[s1,s2]],[x1,y1],[x2,y2],...,[xk,yk]] where the general solution is [x;y]=M^k*[xi;yi]+[s1;s2];\n --------D>0 square EITHER, [[0],[x1,y1],[x2,y2],...,[xk,yk]] where the (xi,yi) are the solutions;\n -----------------------OR, [[2],[[s1,t1],[x1,y1]],...,[[sk,tk],[xk,yk]]] where the solutions are x=xi+Usi and y=yi+Uti for U integer;\n ----------------------D<0, [[0],[x1,y1],[x2,y2],...,[xk,yk]] where the (xi,yi) are the solutions;\n ----------------------Q=0, [[-1]] if n=0 and 0 else;\n -D=0 and A=B=C=0 or D=E=0, [[2],[[s1,t1],[x1,y1]],...,[[sk,tk],[xk,yk]]] where the solutions are x=xi+Usi and y=yi+Uti for U integer (si=sj and ti=tj for all i,j in fact);\n ------------D=0 otherwise, [[-2],[[a1,b1,c1],[e1,f1,g1]],...,[[ak,bk,ck],[ek,fk,gk]]] where the solutions are x=ai*U^2+bi*U+ci, y=ei*U^2+fi*U+gi for U integer.\n\n In general, -2=quadratic, -1=all, 0=finite, 1=positive, 2=linear");
		install("bqf_linearsolve_typecheck","GGGGp","bqf_linearsolve","./libqquadratic.so");
		addhelp(bqf_linearsolve,"Inputs qf, n1, lin, n2: qf a six term integer vector representing the form Ax^2+By^2+Cz^2+Dxy+Exz+Fyz, n1 an integer, lin a three term integer vector representing Ax+By+Cz, and n2 an integer.\n Solves qf(x, y, z)=n1 and lin(x, y, z)=n2 simultaneously. If there are no solutions, this returns 0. Otherwise it returns a vector v. Let v[1][1]=t, and then the format of v is:\n\n -t=-2, v=[[-2], [[a1, a2, a3], [b1, b2, b3], [c1, c2, c3]], ...], where each general solution is x=a1U^2+a2U+a3, y=b1U^2+b2U+b3, z=c1U^2+c2U+c3 for any U integral;\n -t=-1, v=[[-1], [[a1, a2, a3], [b1, b2, b3], [c1, c2, c3]]], where the solution is x=a1U+b1V+c1, y=a2U+b2V+c2, z=a3U+v3V+c3 for any U, V integral;\n --t=0, v=[[0], [a1, b1, c1], ...], where the finite set of solutions are (x,y,z)=(ai, bi, ci);\n --t=1, v=[[1, M, [s1, s2, s3]], [a1, b1, c1], ...], where the general solution is [x;y;z]=M^k[ai;bi;ci]+[s1;s2;s3] for k integral. Note that ai and si need not be integral, though M is.\n --t=2, v=[[2], [[a1, a2, a3], [b1, b2, b3]], ...], where each general solution is x=a1U+b1, y=a2U+b2, z=a3U+b3 for any U integral;\n\n In general, -2=quadratic, -1=plane, 0=finite, 1=positive, 2=linear.");

	\\GENERAL HELP
		addhelp(bqf, "This package deals with binary quadratic forms with integer coefficients. A homogeneous binary quadratic form Ax^2+Bxy+Cy^2 is stored as [A,B,C]. A proper discriminant is an integer that is equivalent to 0 or 1 modulo 4 and is not a square. \n Subtopics:\n Discriminants (disc)\n Basic operations (bqfbasic)\n Indefinite forms (ibqf)\n Class group and composition (bqfclass)\n Representation of numbers (bqfsolve)");
		addhelp(disc,"disclist, discprimeindex, fdisc, isdisc, pell, posreg, quadroot.");
		addhelp(bqfbasic,"bqf_automorph, bqf_disc, bqf_isequiv, bqf_isreduced, bqf_random, bqf_random_D, bqf_red, bqf_roots, bqf_trans, bqf_trans_coprime, ideal_tobqf.");
		addhelp(ibqf,"ibqf_isrecip, ibqf_leftnbr, ibqf_redorbit, ibqf_rightnbr, ibqf_river, ibqf_riverforms, ibqf_symmetricarc, mat_toibqf.");
		addhelp(bqfclass,"bqf_comp, bqf_ncgp, bqf_ncgp_lexic, bqf_pow, bqf_square.");
		addhelp(bqfsolve,"bqf_bigreps, bqf_linearsolve, bqf_reps.");


\\qq_bqf_int

	\\INTERSECTION DATA
		install("bqf_bdelta_typecheck","GG","bqf_bdelta","./libqquadratic.so");
		addhelp(bqf_bdelta, "Inputs q1,q2, integral binary quadratic forms.\n Outputs B_{Delta}(q1,q2)=B1B2-2A1C2-2A2C1, where qi=[Ai,Bi,Ci].");
		install("bqf_intlevel_typecheck","GG","bqf_intlevel","./libqquadratic.so");
		addhelp(bqf_intlevel, "Inputs q1, q2, integral binary quadratic forms.\n Outputs the signed intersection level of q1, q2, i.e. if qi=[Ai, Bi, Ci], then this is the gcd of [-A1B2 + A2B1, -2A1C2 + 2A2C1, -B1C2 + B2C1] times the sign of -A1B2+A2B1.");
		install("ibqf_intpoint_typecheck", "GGD0,G,","ibqf_intpoint","./libqquadratic.so");
		addhelp(ibqf_intpoint, "Inputs q1, q2, {location=0}: PIBQFs q1, q2 whose root geodesics intersect, complex number location, automorph of q1.\n Outputs the upper half plane intersection point of q1, q2. If location=0, it is the intersection of q1, q2; if location=1, we translate it to the fundamental domain of PSL(2,Z); if imag(location)!=0, then location is assumed to be a point on \ell_q1. We translate the intersection point to the geodesic between location and gamma_q1(location). If the invariant automorph is large, then we need to increase the precision to ensure accurate results.");
		install("hdist_typecheck","GGp","hdist","./libqquadratic.so");
		addhelp(hdist,"Inputs z1, z2 complex numbers.\n Outputs the hyperbolic distance between z1 and z2.");

	\\INTERSECTION NUMBER COMPUTATION
		install("ibqf_int_typecheck","GGp","ibqf_int","./libqquadratic.so");
		addhelp(ibqf_int, "Inputs q1, q2 PIBQFs.\n Outputs the total intersection of q1 and q2, i.e. 2(Int_RS(q1,q2)+Int_RS(q1,-q2)).");
		install("ibqf_intRS_typecheck","GGp","ibqf_intRS","./libqquadratic.so");
		addhelp(ibqf_intRS, "Inputs q1, q2 PIBQFs.\n Outputs RS intersection of q1 and q2.");
		install("ibqf_intforms_typecheck","GGD0,L,D0,L,p","ibqf_intforms","./libqquadratic.so");
		addhelp(ibqf_intforms, "Inputs q1, q2, {data=0}, {split=0}: q1, q2 PIBQFs, data and split=0,1.\n Outputs the intersecting forms of q1, q2 of all types. If data=1, each entry of the output is [B_delta(f1,f2), level of int, length of river overlap, f1, f2]; otherwise it is just the pair [f1, f2]. If split=0 outputs a single vector of the return data, if split=1 it splits the output into [[RS], [RO], [LS], [LO]] intersection.");
		install("ibqf_intformsRS_typecheck","GGD0,L,p","ibqf_intformsRS","./libqquadratic.so");
		addhelp(ibqf_intformsRS,"Inputs q1, q2, {data=0}: q1, q2 PIBQFs, data=0,1.\n Outputs the RS intersection of q1 and q2 as a vector of non-simultaneously equivalent intersecting forms. If data=1, each output entry is instead [B_delta(f1,f2), level of int, length of river overlap, f1, f2]");
		install("ibqf_intformsRO_typecheck","GGD0,L,p","ibqf_intformsRO","./libqquadratic.so");
		addhelp(ibqf_intformsRO,"Inputs q1, q2, {data=0}: q1, q2 PIBQFs, data=0,1.\n Outputs the RO intersection of q1 and q2 as a vector of non-simultaneously equivalent intersecting forms. If data=1, each output entry is instead [B_delta(f1,f2), level of int, length of river overlap, f1, f2]");
		install("ibqf_intformsLS_typecheck","GGD0,L,p","ibqf_intformsLS","./libqquadratic.so");
		addhelp(ibqf_intformsLS,"Inputs q1, q2, {data=0}: q1, q2 PIBQFs, data=0,1.\n Outputs the LS intersection of q1 and q2 as a vector of non-simultaneously equivalent intersecting forms. If data=1, each output entry is instead [B_delta(f1,f2), level of int, length of river overlap, f1, f2]");
		install("ibqf_intformsLO_typecheck","GGD0,L,p","ibqf_intformsLO","./libqquadratic.so");
		addhelp(ibqf_intformsLO,"Inputs q1, q2, {data=0}: q1, q2 PIBQFs, data=0,1.\n Outputs the LO intersection of q1 and q2 as a vector of non-simultaneously equivalent intersecting forms. If data=1, each output entry is instead [B_delta(f1,f2), level of int, length of river overlap, f1, f2]");

	\\GENERAL HELP
		addhelp(bqf_int, "This package deals with intersections of binary quadratic forms with integer coefficients. \n Subtopics:\n Computing intersection numbers (intcomp)\n Properties of an intersection (intprop)");
		addhelp(intcomp,"ibqf_int, ibqf_intRS, ibqf_intforms, ibqf_intformsRS, ibqf_intformsRO, ibqf_intformsLS, ibqf_intformsLO.");
		addhelp(intprop,"bqf_bdelta, bqf_intlevel, ibqf_intpoint, hdist.");


\\qq_qquat

	\\BASIC OPERATIONS ON ELEMENTS IN QUATERNION ALGEBRAS 
		install("qa_conj_typecheck","G","qa_conj","./libqquadratic.so");
		addhelp(qa_conj,"Input x, element of a quaternion algebra.\n Outputs the conjugate of x.");
		install("qa_conjby_typecheck","GGG","qa_conjby","./libqquadratic.so");
		addhelp(qa_conjby,"Inputs Q, x, y: quaternion algebra Q, elements x, y with y invertible.\n Outputs yxy^(-1).");
		install("qa_inv_typecheck","GG","qa_inv","./libqquadratic.so");
		addhelp(qa_inv,"Inputs Q, x; quaternion algebra Q, element x. Returns x^(-1), and produces an error if x is not invertible");
		install("qa_m2rembed_typecheck","GG","qa_m2rembed","./libqquadratic.so");
		addhelp(qa_m2rembed,"Inputs Q, x; indefinite quaternion algebra Q, element x.\n Outputs image of x under standard embedding of Q into M(2,R) (requires a>0).");
		install("qa_minpoly_typecheck","GG","qa_minpoly","./libqquadratic.so");
		addhelp(qa_minpoly,"Inputs Q, x; Q a quaternion algebra, x in Q.\n Outputs the minimal polynomial of x. The format is [1,b,c] for x^2+bx+c, and [1,b] for x+b");
		install("qa_mul_typecheck","GGG","qa_mul","./libqquadratic.so");
		addhelp(qa_mul,"Inputs Q, x, y; quaternion algebra Q, elements x, y.\n Returns x*y.");
		install("qa_norm_typecheck", "GG", "qa_norm","./libqquadratic.so");
		addhelp(qa_norm,"Inputs Q, x; Q a quaternion algebra, x in Q.\n Outputs the norm of x.");
		install("qa_pow_typecheck","GGG","qa_pow","./libqquadratic.so");
		addhelp(qa_pow,"Inputs Q, x, n; quaternion algebra Q, element x, integer n.\n Outputs x^n.");
		install("qa_roots_typecheck","GGp","qa_roots","./libqquadratic.so");
		addhelp(qa_roots,"Inputs Q, x; indefinite quaternion algebra Q, element x with positive norm and hyperbolic (trace^2>4*norm).\n Outputs the vector [rt1, rt2] of first root, second root of x under the standard embedding into SL(2,R).");
		install("qa_square_typecheck","GG","qa_square","./libqquadratic.so");
		addhelp(qa_square,"Inputs Q, x: quaternion algebra Q, element x.\n Outputs x^2.");
		install("qa_trace_typecheck","G","qa_trace","./libqquadratic.so");
		addhelp(qa_trace,"Input x, element of a quaternion algebra.\n Outputs the trace of x");

	\\BASIC OPERATIONS ON ORDERS/LATTICES IN QUATERNION ALGEBRAS 
		install("qa_isinorder_typecheck","lGGG","qa_isinorder","./libqquadratic.so");
		addhelp(qa_isinorder,"Inputs Q, ord, v; quaternion algebra Q, an order ord, and an element v.\n Outputs 1 if v is in ord, and 0 if not. Can pass in an initialized order for ord.");
		install("qa_ord_conj_typecheck","GGG","qa_ord_conj","./libqquadratic.so");
		addhelp(qa_ord_conj,"Inputs Q, ord, c; quaternion algebra Q, lattice ord, invertible element c.\n Outputs the lattice c*ord*c^(-1) in Hermite normal form."); 
		install("qa_ord_disc_typecheck","GG","qa_ord_disc","./libqquadratic.so");
		addhelp(qa_ord_disc,"Inputs Q, ord; quaternion algebra Q, order ord.\n Ouputs the discriminant of the order ord.");

	\\INITIALIZATION METHODS
		install("qa_ord_init_typecheck","GG","qa_ord_init","./libqquadratic.so");
		addhelp(qa_ord_init,"Inputs: Q, ord: a quaternion algebra and an order.\n Outputs the initialization of the order: [ord, type, [d1, d2, d3, d4], level, prime factorization of the level]. ord is the hnf of the inputted ord, and d_n is the maximal denominator of the nth coefficient (i.e. 1, i, j, k). type=-1 means general order, type=0 means Maximal, type=1 means Eichler.");
		install("qa_init_primes_typecheck","G","qa_init_primes","./libqquadratic.so");
		addhelp(qa_init_primes,"Input pset, a list of primes (or oo).\n Outputs a quaternion algebra Q with ramification at the specified set of primes. If the list has odd length, then automatically adds oo to the ramification set.");
		install("qa_init_2primes_typecheck","GG","qa_init_2primes","./libqquadratic.so");
		addhelp(qa_init_2primes,"Inputs p, q distinct primes.\n Outputs a quaternion algebra over Q with ramification at p,q only.");
		install("qa_ram_fromab_typecheck","GG","qa_ram_fromab","./libqquadratic.so");
		addhelp(qa_ram_fromab,"Inputs a, b, non-zero integers.\n Outputs the set of primes ramifying in the quaternion algebra (a, b / Q).");

	\\CONJUGATION OF ELEMENTS IN A GIVEN ORDER
		install("qa_conjbasis_typecheck","GGGGD0,L,","qa_conjbasis","./libqquadratic.so");
		addhelp(qa_conjbasis,"Inputs: Q, ord, e1, e2, {orient=0}; Q is the quaternion algebra, ord is an (initialized) order, e1 and e2 are the elements of Q, and orient=0,1.\n We test if e1 and e2 are conjugate in Q. If e1 or e2 is rational we automatically return 0. Else, if they are conjugate, the conjugation space is of rank 2. We return a Z basis for this space intersected with ord. If orient=0 we orient it correctly, otherwise we don't care.");
		install("qa_conjqf_typecheck","GGGG","qa_conjqf","./libqquadratic.so");
		addhelp(qa_conjqf,"Inputs: Q, ord, e1, e2; quaternion algebra Q, (initialized) order ord, non rational elements e1, e2 of Q with the same minimal polynomial.\n Outputs [nrd(Xv_1+Yv_2), v1, v2] where Zv_1+Zv_2 is the Z-module for which ve1=e2v and v is in ord. We take it so v_2*conj(v_1)-v_1*conj(v_2) corresponds to the same embedding as e2.");
		install("qa_conjnorm_typecheck","GGGGGD0,L,p","qa_conjnorm","./libqquadratic.so");
		addhelp(qa_conjnorm,"Inputs: Q, ord, e1, e2, n, {retconj=0}; quaternion algebra Q, order ord, elements e1, e2 of Q, integer n, retconj=0,1.\n Output is 1 if there is an element v of ord of norm n conjugating e1 to e2 (v*e1*v^(-1)=e2), or 0 if no such element exists. If retconj=1 we also output a possible v.");
		install("qa_simulconj_typecheck","GGGGGGp","qa_simulconj","./libqquadratic.so");
		addhelp(qa_simulconj,"Inputs: Q, ord, e1, e2, f1, f2; quaternion algebra Q, (initialized) order ord, elements e1, e2, f1, f2 such that e1, f1 are conjugate, e2, f2 are conjugate, and trace(e1e2)=trace(f1f2). \n Then (e1 ,e2) and (f1, f2) are simultaneously conjugate, and the conjugation space has dimension 1. This returns a generator for this space intersected with ord (isomorphic to Z, so the generator is unique up to +/-).");

	\\EMBEDDING QUADRATIC ORDERS INTO EICHLER ORDERS
		install("qa_associatedemb_typecheck","GGGD0,G,","qa_associatedemb","./libqquadratic.so");
		addhelp(qa_associatedemb,"Inputs Q, order, emb,{D=0}: Q a quaternion algebra, order an initialized order, emb the image of image of (A+sqrt(D))/2, (D a discriminant which does not need to be passed).\n This method computes the pair [emb', D'] consisting of emb', the optimal embedding associated to emb and ord, and its corresponding discriminant D'.");
		install("qa_embed_typecheck","GGGD0,G,D0,L,p","qa_embed","./libqquadratic.so");
		addhelp(qa_embed,"Inputs Q, order, D, {nembs=0}, {retpell=0}; indefinite quaternion algebra Q, Eichler order order, discriminant D, nonnegative integer nembs, retpell=0,1.\n Outputs the optimal embeddings of O_D into order up to conjugation by order_1^x. If retpell=1 we return the images of the fundamental unit in O_D (D>0 only), and otherwise we return the images of (D%2+sqrt(D))/2. Can pass the number of embeddings nembs if we have already counted them, or if we want a smaller number of embeddings (e.g. 1) can also pass that. If we pass a larger number than the number of embeddings the method will not end.");
		install("qa_embeddablediscs_typecheck","GGGGD0,L,D0,G,","qa_embeddablediscs","./libqquadratic.so");
		addhelp(qa_embeddablediscs,"Inputs Q, order, d1, d2, {fund=0}, {cop=0}: Q an indefinite quaternion algebra, initialized Eichler order order, d1<d2 integers, fund=0,1, cop a nonnegative integer.\n Outputs the discriminants D with d1<=D<=d2 for which there exists an optimal embedding of discriminant D into Q[3]. If fund=1 only outputs fundamental discriminants, and if cop!=0 only outputs discriminants coprime to cop.");
		install("qa_numemb_typecheck","GGGD0,G,p","qa_numemb","./libqquadratic.so");
		addhelp(qa_numemb,"Inputs Q, order, D, {narclno=0}; Q an indefinite quaternion algebra, order an initialized Eichler order, D a discriminant, narclno the narrow class number of discriminant D.\n Outputs [m, n, v1, v2, v3]. Here, m is the number of optimal embeddings of O_D into order, and n=h^+(D) is the number of a fixed orientation.\n v1=[x] with x the number of embeddings at oo (1 if D>0 and 2 if D<0).\n v2=[x_1, ..., x_r] with Q[1]=[p_1, ..., p_r] and x_i=number of local embeddings at p_i for i<=r (1 if p_i divides D and 2 else).\n v3=[y_1, ..., y_s], where s distinct primes (see order[5]) q_1,...,q_s divide the level of the order and y_i is the number of local embedding classes. This is 2 as long as q_i does not divide D.");
		install("qa_ordiffer_typecheck","GGGGD0,G,","qa_ordiffer","./libqquadratic.so");
		addhelp(qa_ordiffer,"Inputs Q, order, e1, e2, {D=0}: indefinite quaternion algebra Q, initialized Eichler order order, e1 and e2 optimal embeddings of O_D into order.\n Outputs the set of finite primes for which the orientation of the embeddings at that prime differs.");
		install("qa_orinfinite_typecheck","GGD0,G,p","qa_orinfinite","./libqquadratic.so");
		addhelp(qa_orinfinite,"Inputs Q, emb, {D=0}: indefinite quaternion algebra Q, embedding emb of discriminant D.\n Outputs the orientation of the embedding emb with respect to the standard embedding into M(2, R) (0 if D>0, and +/-1 if D<0).");
		install("qa_sortedembed_typecheck","GGGD0,L,D0,G,p","qa_sortedembed","./libqquadratic.so");
		addhelp(qa_sortedembed,"Inputs Q, order, D, {rpell=0}, {ncgp=0}; indefinite quaternion algebra Q, initialized Eichler order order, discriminant D, rpell=0,1, ncgp=0 or the output of bqf_ncgp_lexic(D).\n Outputs the optimal embeddings of O_D into order sorted by embedding and group action. Rows correspond to a fixed orientation, with the first entry being the primes at which the orientation differs from the first line. The subsequent h^+(D) entries are the optimal embeddings, arranged so that lexicorder(D,1)[i] acts on form in position j takes it to the form in position i+j. Conjugating a row of embeddings by an Atkin Lehner involution takes the row to the row with the corresponding orientation.");

	\\SUPPORTING METHODS
		install("module_intersect_typecheck","GG","module_intersect","./libqquadratic.so");
		addhelp(module_intersect,"Inputs A, B: Z-modules spanned by the columns.\n Outputs the intersection of the Z-modules, as a Z-module spanned by the colunmns.");
		install("prime_ksearch_typecheck","GD0,G,","prime_ksearch","./libqquadratic.so");
		addhelp(prime_ksearch, "Inputs rels, {extra=0}; rels=[[p_1,s_1],...,[p_k,s_k]], with the p_i distinct integers and s_i=-1 or 1 for each i, and extra=0 or [n, c].\n Outputs the smallest prime r such that kronecker(r,p_i)=s_i for all i, and r== c_i mod n_i for all i. If the inputs are not consistent, the program will issue an error before the memory has been exhausted.");
		install("QM_hnf_typecheck","G","QM_hnf","./libqquadratic.so");
		addhelp(QM_hnf,"Input M a rational matrix.\n Outputs the Hermite normal form of M with resepct to columns.")
		install("powerset","G","powerset","./libqquadratic.so");
		addhelp(powerset,"Input L, a vector.\n Outputs all subsets of L.");

	\\GENERAL HELP
		\\TO DO!!!!!!!!!!!!!!!!!!!!!!


\\qq_qquat_int

	\\INTERSECTION NUMBER BASED OFF OF ROOTS BOUND AREA
		install("qa_inum_roots_typecheck","GGGGD1,L,p","qa_inum_roots","./libqquadratic.so");
		addhelp(qa_inum_roots,"Inputs Q, order, e1, e2, {data=1}: indefinite quaternion algebra Q, initialized Eichler order order, e1, e2 not rational and lying in real quadratic subfields, data=0, 1.\n This calculates the intersections of e1 and e2 with respect to order using the root geodesic method. If data=0, only outputs the pairs, otherwise the output is [[intersections],[data]], where intersections[i] is a pair of optimal embeddings that intersect, and data[i] is the pair [signed level,x] corresponding to the embedding pair.");

	\\INTERSECTIONS BASED ON X-LINKAGE
		install("qa_inum_x_typecheck","GGGGD0,L,p","qa_inum_x","./libqquadratic.so");
		addhelp(qa_inum_x,"Inputs Q, order, e1, e2, {data=1}: indefinite quaternion algebra Q, initialized Eichler order order, e1, e2 optimal embeddings into order, data=0, 1.\n This calculates the intersections number of e1 and e2 with respect to order, via x-linking.");
		install("qa_xlink_typecheck","GGGGGp","qa_xlink","./libqquadratic.so");
		addhelp(qa_xlink,"Inputs Q, ord, e1, e2, x: indefinite quaternion algebra Q, initialized Eichler order ord, optimal embeddings e1, e2, integer x.\n Outputs the set of pairs of x-linked optimal embeddings conjugate to (e1,e2) taken up to simultaneous conjugation.");
		install("qa_xposs_typecheck","GGGD0,G,D0,G,","qa_xposs","./libqquadratic.so");
		addhelp(qa_xposs,"Inputs P, D1, D2, {xmin=0}, {xmax=0}: even sized set of primes OR quaternion algebra P, discriminants D1, D2, 0<=xmin<=xmax.\n Outputs the set of x's between xmin and xmax (inclusive) which exhbit x-linking in the indefinite quaternion algebra ramified at P (or P itself when it is a quaternion algebra). If xmin=xmax=0, outputs x's between 0 and sqrt(D_1D_2).");
	
	\\GENERAL HELP
		\\TO DO!!!!
		

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
		addhelp(hist,"hist_make, hist_rebin, hist_recompile, hist_rerange, hist_rescale");

