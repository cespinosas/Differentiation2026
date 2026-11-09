nupre_2=3;
numnompremot_2=2;
LM_2= new GraphB[nupre_2];
for (i=0; i<nupre_2; i++) {
	LM_2[i].start_rng(est);
	LM_2[i].make_nw(2, false);
}
nompremot_2= new int*[nupre_2];
for (i=0; i<nupre_2; i++) {
	nompremot_2[i]= new int[numnompremot_2];
}
socont_2= new int[nupre_2];
basic.fillv0(socont_2, nupre_2);
basic.fillmat0(nompremot_2, nupre_2, numnompremot_2);
socont_2[0] = 1;
nompremot_2[0][0] = 2;
nompremot_2[0][1] = 0;
LM_2[0].force_interaction_undir(0, 1, true);
socont_2[1] = 2;
nompremot_2[1][0] = 1;
nompremot_2[1][1] = 1;
LM_2[1].force_interaction_undir(0, 0, true);
LM_2[1].force_interaction_undir(0, 1, true);
socont_2[2] = 1;
nompremot_2[2][0] = 0;
nompremot_2[2][1] = 2;
LM_2[2].force_interaction_undir(0, 0, true);
LM_2[2].force_interaction_undir(0, 1, true);
LM_2[2].force_interaction_undir(1, 1, true);
