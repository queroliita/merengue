#include "merengue.c"
#include <stdio.h>

/*
gcc testsuite.c testsuite
./testsuite {ID word} {ID bit} {OD word} {OD bit}
*/

/* Recommended arguments:
   ID = {7,31} OD = {1,31}
   ID = {8,31} OD = {6,31}
   ID = {7,31} OD = {1,31} CIV1 = {7,19,12}
   ID = {8,31} OD = {6,31} CIV1 = {8,19,12}
   ID = {7,31} OD = {1,31} CIV1 = {7,19,12} CIV2 = {6,25,10}
 */

int ID[2], OD[2];
PNB bpnb, fpnb;
CIV civ;

int main(int argc, char *argv[])
{ 
  if ( argc < 5 ) exit(1);

  ID[0] = atoi(argv[1]);
  ID[1] = atoi(argv[2]);
  OD[0] = atoi(argv[3]);
  OD[1] = atoi(argv[4]);
  
  if ( argc >= 8 ){
    civ.word[civ.len] = atoi(argv[5]); 
    civ.ini[civ.len] = atoi(argv[6]);
    civ.end[civ.len] = atoi(argv[7]);
    civ.len += 1;
  }
  if ( argc == 11 ){
    civ.word[civ.len] = atoi(argv[8]); 
    civ.ini[civ.len] = atoi(argv[9]);
    civ.end[civ.len] = atoi(argv[10]);
    civ.len += 1;
  }

  clock_t t;
  t = clock(); 

  /* TEST 1: get bias Ed */
  /*
  printf("%f = Ed with ID = [%d,%d]\n", getEd(ID,OD,NULL),ID[0],ID[1]);
  if (argc >= 8) {
    printf("%f = conditioned Ed with ID = [%d,%d]\n", getEd(ID,OD,&civ),ID[0],ID[1]);
  }
  */
  /* 6.5 seconds */
  
  /* TEST 2: compute PNBs */
  /*
  getBPNBs(&bpnb,0.12,ID,OD,&civ);
  for (int i = 0; i < bpnb.len; i++) {
    if ( iskey(bpnb.word[i]) ) printf("BPNKB (%d,%d) ",bpnb.word[i],bpnb.bit[i]);
    else printf("BPNB (%d,%d) ",bpnb.word[i],bpnb.bit[i]);
    printf("with neutrality measure = %f\n",bpnb.bias[i]);
    }
  } printf("-> Found %d BPNBs and %d key\n",bpnb.len,bpnb.n);
  */
  /* 730 seconds */

  /* TEST 3: get bias Ea after fixing PNBs */
  /*PNB_list PNB_1_14 = { 63, 36,
    {  1, 1, 1, 1, 1, 1, 3, 3, 4, 4, 4,11,
      12,12,12,12,12,12,12,12,12,12,12,12,
      12,13,13,13,14,14,14,14,14,14,14,14 },
    { 26,27,28,29,30,31, 7, 8,24,25,26,20,
       5, 6, 7, 8, 9,10,11,12,13,14,15,16,
      17,18,19,20, 0, 1,18,19,20,21,22,23 }
  };*/
  /*PNB_list PNB_6_14 = { 63, 30,
    { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
      1, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3,
      4, 4,11,11,11,12},
    { 5, 6, 7, 8, 9,10,11,12,13,14,15,16,
     17,18,19,20, 0, 1,18,19,20,21,22,23
      7, 8,13,14,31,20}
  };*/
      /*
  printf("Bias |Ea*| = %f\n",getEa(&PNB,ID,OD));*/
  /* 3.5 seconds */

  /* TEST : get forwards PNBs */
  getFPNBs(&fpnb,0.12,ID,OD,NULL);
  for (int i = 0; i < fpnb.len; i++) {
    printf("FPNB (%d,%d) ",fpnb.word[i],fpnb.bit[i]);
    printf("with neutrality measure = %f\n",fpnb.bias[i]);
  } printf("-> Found %d FPNBs and %d ivs\n",fpnb.len,fpnb.n);
  /**/

  t = clock() - t;
  double time_taken = ((double)t)/CLOCKS_PER_SEC;
  printf("\n%f seconds\n",time_taken);
}

