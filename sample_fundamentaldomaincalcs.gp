Q=qa_init_primes([3, 11]);\\Quaternion algebra ramified at 3, 11
order=qa_eichlerorder(Q, 5);\\Eichler order of level 5
S=qa_embeddablediscs(Q, order, 2000, 2100, 1, 3);\\Discriminants for which there exist optimal embeddings into order that are between 2000 and 2100, fundamental, and coprime to 3
E=qa_embed(Q, order, S[1], 0, 1);\\All non-conjugate embeddings of discriminant S[1], returning the images of the fundamental units
U=qa_fundamentaldomain(Q, order, I/2, 1);\\The fundamental domain with centre I/2
g1=qa_rootgeodesic_fd(Q, U, E[1]);\\The root geodesic of E[1]
g2=qa_rootgeodesic_fd(Q, U, E[3]);\\The root geodesic of E[3]
python_printfdom(U, "fd_33_5");\\Saves the fundamental domain to "fdoms/fd_407_23". It is important to start the filename with "fd".
python_printarcs(g1[2], "geod1", 0);\\Saves g1. It is important to start the filename with anything other than "fd".
python_printarcs(g2[2], "geod2", 0);\\Saves g2.
python_plotviewer("fd_33_5 geod1 geod2");\\Displays the fundamental domain and the two geodesics, but ONLY if you are using Windows subsystem for Linux. Otherwise, you need to call "py fdviewer.py fd_33_5 geod1 geod2" from the command line.
