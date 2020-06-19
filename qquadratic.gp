\\c_base

install("crossratio","GGGG","crossratio","./libqquadratic.so");
addhelp(crossratio, "Inputs a,b,c,d complex numbers or oo with at most one being oo.\n Outputs the crossratio of [a,b;c,d].");
install("mat_eval_typecheck","GG","mat_eval","./libqquadratic.so");
addhelp(mat_eval, "Inputs M,x; M a matrix, and x number.\n Outputs Mx with M acting via Mobius transformation. x=+/-oo is allowed.");
install("addoo","GG","addoo","./libqquadratic.so");
addhelp(addoo, "Inputs a,b real numbers or oo.\n Outputs a+b, where if a or b is +/- oo, returns that back. Note that this will make oo+-oo=oo and -oo+oo=-oo");
install("divoo","GG","divoo","./libqquadratic.so");
addhelp(divoo, "Inputs a,b, real numbers or oo.\n Outputs a/b, where oo is output if a=oo and b>=0 or a=-oo and b<0 or b=0 and a>=0 (outputs -oo under analogous assumptions).");
install("lin_intsolve_typecheck","GGG","lin_intsolve","./libqquadratic.so");
addhelp(lin_intsolve, "Inputs A,B,n integers.\n Outputs the general integral solutions to Ax+By=n. The format is [[s1,s2],[x0,y0]], where the general solution is x=s1*t+x0, y=s2*t+y0 for t an integer. The output is also reduced, i.e. gcd(s1,s2)=1. If A=B=0 or there are no integer solutions, returns 0.");
install("mat3_complete_typecheck", "GGG", "mat3_complete", "./libqquadratic.so");
addhelp(mat3_complete, "Inputs A,B,C integers with gcd 1.\n Outputs a 3x3 integer matrix with determinant 1 whose top row is [A, B, C].");
install("sqmod_typecheck","GG","sqmod","./libqquadratic.so");
addhelp(sqmod, "Inputs: x,n: rational number x, integer n with gcd(n,denom(x))=1.\n Solves y^2==x mod n, and outputs the solutions. The output format is 0 if no solutions, and [S,m] otherwise, where the solutions are y==S[i] modulo m.");
install("printtime","v","printtime","./libqquadratic.so");
addhelp(printtime, "Prints current time");
addhelp(base,"This package is a collection of miscellaneous methods that may be useful in a variety of settings, and not just for the programs they were originally created for \n Subtopics: \n Complex plane geometry (cpg) \n Infinity (inf) \n Integer linear equations and matrices (ilm) \n Modulo n methods (modn) \n Time (time)");
addhelp(cpg,"crossratio, mat_eval.");
addhelp(inf,"addoo, divoo.");
addhelp(ilm,"lin_intsolve, mat3_complete.");
addhelp(modn,"sqmod.");
addhelp(time,"printtime.");

\\c_bqf

install("disclist","GGD0,L,D0,G,","disclist","./libqquadratic.so");
addhelp(disclist, "Inputs d1, d2,{fund=0}, {cop=0}: d1 and d2 integers with d1<=d2, fund=0,1, cop integer.\n Outputs the set of proper discriminants between d1 and d2, inclusive. If fund=1, only outputs fundamental discriminants. If cop!=0, only outputs discriminants coprime to cop.");
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
addhelp(quadroot, "Input D, a proper discriminant.\n Outputs sqrt(D) of type t_QUAD.");
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
install("bqf_reps_typecheck","GGD0,L,D1,L,p","bqf_reps","./libqquadratic.so");
addhelp(bqf_reps,"Inputs q, n, {proper=0}, {half=1}: BQF q, integer n, (proper=0,1), (half=0,1).\n This solves the equation q(x,y)=n over the integers. We will get a finite set of (families) of solutions, and if half=0 we only output one of the families corresponding to (x,y) and (-x,-y). If we want only coprime solutions when disc!=square, pass in proper=1. If Q has discriminant D, the output is:\n\n --------If no solutions, returns 0;\n -D>0 nonsquare and n!=0, [[1,M],[x1,y1],[x2,y2],...,[xk,yk]] where the +/-(xi,yi) are representatives of the distinct classes of solutions and M is the invariant automorph;\n ----D>0 square and n!=0, [[0],[x1,y1],[x2,y2],...,[xk,yk]] where the +/-(xi,yi) are the solutions;\n --------------------D<0, [[0],[x1,y1],[x2,y2],...,[xk,yk]] where the +/-(xi,yi) are the solutions;\n -----------D=0 and n!=0, [[2],[[s1,s2],[x1,y1]]] where the solutions are (up to +/-) x=x1+Us1, y=y1+Us2;\n ---------n=0, D!=square, [[0],[0,0]];\n -------n=0, D=square!=0, [[2],[[s1,s2],[0,0]],[[s3,s4],[0,0]]], solutions are (x,y)=(s1k,s2k),(s3k,s4k) for integer k;\n ------------n=0 and D=0, if Q!=0 as above, if Q=0 then [[-1]] (everything is a solution).\n\n In general, -1=all, 0=finite, 1=positive, 2=linear");
install("bqf_bigreps_typecheck","GGp","bqf_bigreps","./libqquadratic.so");
addhelp(bqf_bigreps,"Inputs: Q, n, Q=[A,B,C,D,E]=Ax^2+Bxy+Cy^2+Dx+Ey and an integer n.\n Outputs the solutions to Q(x,y)=n over the integers. If D=bqf_disc(Q)=B^2-4AC, then the output is: \n\n ----------If no solutions, returns 0; \n ------------D>0 nonsquare, [[1,M,[s1,s2]],[x1,y1],[x2,y2],...,[xk,yk]] where the general solution is [x;y]=M^k*[xi;yi]+[s1;s2];\n --------D>0 square EITHER, [[0],[x1,y1],[x2,y2],...,[xk,yk]] where the (xi,yi) are the solutions;\n -----------------------OR, [[2],[[s1,t1],[x1,y1]],...,[[sk,tk],[xk,yk]]] where the solutions are x=xi+Usi and y=yi+Uti for U integer;\n ----------------------D<0, [[0],[x1,y1],[x2,y2],...,[xk,yk]] where the (xi,yi) are the solutions;\n ----------------------Q=0, [[-1]] if n=0 and 0 else;\n -D=0 and A=B=C=0 or D=E=0, [[2],[[s1,t1],[x1,y1]],...,[[sk,tk],[xk,yk]]] where the solutions are x=xi+Usi and y=yi+Uti for U integer (si=sj and ti=tj for all i,j in fact);\n ------------D=0 otherwise, [[-2],[[a1,b1,c1],[e1,f1,g1]],...,[[ak,bk,ck],[ek,fk,gk]]] where the solutions are x=ai*U^2+bi*U+ci, y=ei*U^2+fi*U+gi for U integer.\n\n In general, -2=quadratic, -1=all, 0=finite, 1=positive, 2=linear");
install("bqf_linearsolve_typecheck","GGGGp","bqf_linearsolve","./libqquadratic.so");
addhelp(bqf_linearsolve,"Inputs qf, n1, lin, n2: qf a six term integer vector representing the form Ax^2+By^2+Cz^2+Dxy+Exz+Fyz, n1 an integer, lin a three term integer vector representing Ax+By+Cz, and n2 an integer.\n Solves qf(x, y, z)=n1 and lin(x, y, z)=n2 simultaneously. If there are no solutions, this returns 0. Otherwise it returns a vector v. Let v[1][1]=t, and then the format of v is:\n\n -t=-2, v=[[-2], [[a1, a2, a3], [b1, b2, b3], [c1, c2, c3]], ...], where each general solution is x=a1U^2+a2U+a3, y=b1U^2+b2U+b3, z=c1U^2+c2U+c3 for any U integral;\n -t=-1, v=[[-1], [[a1, a2, a3], [b1, b2, b3], [c1, c2, c3]]], where the solution is x=a1U+b1V+c1, y=a2U+b2V+c2, z=a3U+v3V+c3 for any U, V integral;\n --t=0, v=[[0], [a1, b1, c1], ...], where the finite set of solutions are (x,y,z)=(ai, bi, ci);\n --t=1, v=[[1, M, [s1, s2, s3]], [a1, b1, c1], ...], where the general solution is [x;y;z]=M^k[ai;bi;ci]+[s1;s2;s3] for k integral. Note that ai and si need not be integral, though M is.\n --t=2, v=[[2], [[a1, a2, a3], [b1, b2, b3]], ...], where each general solution is x=a1U+b1, y=a2U+b2, z=a3U+b3 for any U integral;\n\n In general, -2=quadratic, -1=plane, 0=finite, 1=positive, 2=linear.");
addhelp(bqf, "This package deals with binary quadratic forms with integer coefficients. A homogeneous binary quadratic form Ax^2+Bxy+Cy^2 is stored as [A,B,C]. A proper discriminant is an integer that is equivalent to 0 or 1 modulo 4 and is not a square. \n Subtopics:\n Discriminants (disc)\n Basic operations (bqfbasic)\n Indefinite forms (ibqf)\n Class group and composition (bqfclass)\n Representation of numbers (bqfsolve)");
addhelp(disc,"disclist, discprimeindex, fdisc, isdisc, pell, posreg, quadroot.");
addhelp(bqfbasic,"bqf_automorph, bqf_disc, bqf_isequiv, bqf_isreduced, bqf_random, bqf_random_D, bqf_red, bqf_roots, bqf_trans, bqf_trans_coprime, ideal_tobqf.");
addhelp(ibqf,"ibqf_isrecip, ibqf_leftnbr, ibqf_redorbit, ibqf_rightnbr, ibqf_river, ibqf_riverforms, ibqf_symmetricarc, mat_toibqf.");
addhelp(bqfclass,"bqf_comp, bqf_ncgp, bqf_ncgp_lexic, bqf_pow, bqf_square.");
addhelp(bqfsolve,"bqf_bigreps, bqf_linearsolve, bqf_reps.");
