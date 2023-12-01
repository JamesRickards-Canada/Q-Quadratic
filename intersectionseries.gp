/*
This file details the computations in "HECKE OPERATORS ACTING ON OPTIMAL EMBEDDINGS IN INDEFINITE QUATERNION ALGEBRAS"

Load both qquadratic and this file into a GP session, and call any of the 4 is_ex() methods here to test that method.
Call "makeseries" on any output to make a LaTeX version of the output.
*/

is_ex11(N=50, dooutput=1)={\\N coefficients
  my(Q, order, U, e1, e2, S);
  [Q, order]=qa_init_2primes(5, 7);\\Algebra and order
  if(dooutput,
    printf("The algebra is B=(%d, %d/Q).\nThe order is spanned by the columns of :", Q[3][1], Q[3][2]);
    printp(order[1]);
  );
  U=qa_fundamentaldomain(Q, order);\\Fundamental domain
  e1=qa_embed(Q, order, 5)[1];\\First embedding
  e2=qa_embed(Q, order, 12)[2];\\Second embedding
  if(dooutput,
    printf("The two embeddings are e1=%Ps and e2=%Ps,\n", e1, e2);
  );
  S=qa_inumseries(Q, order, U, e1, e2, N);
  if(dooutput,
    print("The non-zero coefficients are:");
    for(i=1,N,
	  if(S[i]!=0,printf("%d q^%d\n", S[i], i));
	);
  );
  return(S);
}

is_ex12(N=50, dooutput=1)={\\N coefficients
  my(Q, order, U, e1, e2, S);
  [Q, order]=qa_init_2primes(5, 7);\\Algebra and order
  if(dooutput,
    printf("The algebra is B=(%d, %d/Q).\nThe order is spanned by the columns of :", Q[3][1], Q[3][2]);
    printp(order[1]);
  );
  U=qa_fundamentaldomain(Q, order);\\Fundamental domain
  e1=qa_embed(Q, order, 12)[2];\\First embedding
  e2=qa_embed(Q, order, 173)[3];\\Second embedding
  if(dooutput,
    printf("The two embeddings are e1=%Ps and e2=%Ps,\n", e1, e2);
  );
  S=qa_inumseries(Q, order, U, e1, e2, N);
  if(dooutput,
    print("The non-zero coefficients are:");
    for(i=1,N,
	  if(S[i]!=0,printf("%d q^%d\n", S[i], i));
	);
  );
  return(S);
}

is_ex21(N=100, dooutput=1)={\\N coefficients
  my(Q, order, U, e1, e2, S);
  [Q, order]=qa_init_2primes(2, 7);\\Algebra and order
  order=qa_eichlerorder(Q, 3, order[1]);\\Eichler order level 3
  if(dooutput,
    printf("The algebra is B=(%d, %d/Q).\nThe order is spanned by the columns of :", Q[3][1], Q[3][2]);
    printp(order[1]);
  );
  U=qa_fundamentaldomain(Q, order);\\Fundamental domain
  e1=qa_embed(Q, order, 13)[1];\\First embedding
  e2=qa_embed(Q, order, 24)[2];\\Second embedding
  if(dooutput,
    printf("The two embeddings are e1=%Ps and e2=%Ps,\n", e1, e2);
  );
  S=qa_inumseries(Q, order, U, e1, e2, N);
  if(dooutput,
    print("The non-zero coefficients are:");
    for(i=1,N,
	  if(S[i]!=0,printf("%d q^%d\n", S[i], i));
	);
  );
  return(S);
}

is_ex22(N=100, dooutput=1)={\\N coefficients
  my(Q, order, U, e1, e2, S);
  [Q, order]=qa_init_2primes(2, 7);\\Algebra and order
  order=qa_eichlerorder(Q, 3, order[1]);\\Eichler order level 3
  if(dooutput,
    printf("The algebra is B=(%d, %d/Q).\nThe order is spanned by the columns of :", Q[3][1], Q[3][2]);
    printp(order[1]);
  );
  U=qa_fundamentaldomain(Q, order);\\Fundamental domain
  e1=qa_embed(Q, order, 24)[2];\\First embedding
  e2=qa_embed(Q, order, 45)[2];\\Second embedding
  if(dooutput,
    printf("The two embeddings are e1=%Ps and e2=%Ps,\n", e1, e2);
  );
  S=qa_inumseries(Q, order, U, e1, e2, N);
  if(dooutput,
    print("The non-zero coefficients are:");
    for(i=1,N,
	  if(S[i]!=0,printf("%d q^%d\n", S[i], i));
	);
  );
  return(S);
}

makeseries(v)={
  my(S, fnonz=1);
  S="$";
  for(i=1,length(v),
    if(v[i]>0,
	  if(fnonz,\\First time non-zero
	    fnonz=0;
		if(v[i]>1,
	      S=strprintf("%s%dq^{%d}", S, v[i], i);
		,
		  S=strprintf("%sq^{%d}", S, i);
		);
	  ,
	    if(v[i]>1,
	      S=strprintf("%s+%dq^{%d}", S, v[i], i)
	    ,
	      S=strprintf("%s+q^{%d}", S, i)
	    );
      );
	,
	  if(v[i]<0,
	    if(fnonz,fnonz=0);
	    if(v[i]<-1,
	      S=strprintf("%s%dq^{%d}", S, v[i], i)
	    ,
	      S=strprintf("%s-q^{%d}", S, i)
	    );
	  );
	);
  );
  S=strprintf("%s$", S);
  return(S);
}
