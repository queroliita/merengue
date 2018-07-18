#include "merengue.c"
#include <time.h>
#include <sys/time.h>

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

int ID[2], OD[2], r;
PNB bpnb, fpnb, fpsb;
CIV civ;

PNB pnb_1_14 = { 63, 36,
    {  1, 1, 1, 1, 1, 1, 3, 3, 4, 4, 4,11,
      12,12,12,12,12,12,12,12,12,12,12,12,
      12,13,13,13,14,14,14,14,14,14,14,14 },
    { 26,27,28,29,30,31, 7, 8,24,25,26,20,
       5, 6, 7, 8, 9,10,11,12,13,14,15,16,
      17,18,19,20, 0, 1,18,19,20,21,22,23 } };
PNB pnb_6_14 = { 63, 30,
  { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 
    2, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4,11,11,11,12},
  { 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,
    20, 0, 1,18,19,20,21,22,23,7, 8,13,14,31,20} };
PNB pnbciv7 = { 68, 37,
  {  1, 1, 1, 1, 1, 1, 3, 3, 4, 4, 4,11,
    12,12,12,12,12,12,12,12,12,12,12,12,12,
    13,13,13,14,14,14,14,14,14,14,14,14 },
  { 26,27,28,29,30,31, 7, 8,24,25,26,20,
    5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,
    18,19,20, 0, 1, 2,18,19,20,21,22,23 } };
PNB fpnb_std = { 134, 47,
  {  6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
     7, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
     8, 8, 8, 8, 9, 9, 9, 9, 9, 9, 9, 9,
     9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9 },
  {  0,10,22,23,24,25,26,27,28,29,30,31,
    18,19, 1, 2, 3,11,12,23,24,25,26,27,
    28,29,30,31, 5, 6, 7, 8, 9,10,11,13,
    14,15,16,17,24,25,27,28,29,30,31 } };


static clock_t timeit(clock_t t){
  t = clock() - t;
  double time_taken = ((double)t)/CLOCKS_PER_SEC;
  printf("\n%f seconds\n",time_taken);
  return t;
}


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

  // Seed random generator
  struct timeval t1;
  gettimeofday(&t1, NULL);
  srand(t1.tv_usec * t1.tv_sec); //srand(time(NULL));

  clock_t t;
  t = clock(); 

  r = 4;

  /* TEST 1: get bias Ed */
  /*
  if (argc == 5) {
    printf("%f = Ed with ID = [%d,%d]\n", getEd(r,ID,OD,NULL,NULL),ID[0],ID[1]);
  } else {
    printf("%f = conditioned Ed with ID = [%d,%d]\n", getEd(r,ID,OD,&civ,NULL),ID[0],ID[1]);
  }
  */
    
  /* TEST 2: compute BPNBs */
  /*
  if (argc<=5) getBPNBs(&bpnb,0.12,ID,OD,NULL);
  else getBPNBs(&bpnb,0.12,ID,OD,&civ);
  for (int i = 0; i < bpnb.len; i++) {
    if ( iskey(bpnb.word[i]) ) printf("BPNKB (%d,%d) ",bpnb.word[i],bpnb.bit[i]);
    else printf("BPNB (%d,%d) ",bpnb.word[i],bpnb.bit[i]);
    printf("with neutrality measure = %f\n",bpnb.bias[i]);
  } printf("-> Found %d BPNBs and %d key\n",bpnb.len,bpnb.n);
  */

  /* TEST : get bias Ea after fixing PNBs */
  /*
  ID[0] = 7; ID[1] = 31; OD[0] = 1; OD[1] = 14; 
  printf("Bias |Ea*| = %f\n",getEa(&pnb_1_14,ID,OD,NULL));
  */
  /*
  ID[0] = 8; ID[1] = 31; OD[0] = 6; OD[1] = 14;
  printf("Bias |Ea*| = %f\n",getEa(&pnb_6_14,ID,OD,NULL));
  */

  /* TEST: Ea conditioned on X7 -> 7 31 1 14 7 9 12 */
  /*
  ID[0] = 7; ID[1] = 31; OD[0] = 1; OD[1] = 14; 
  civ.word[0] = 7; civ.ini[0] = 19; civ.end[0] = 12; civ.len = 1;
  printf("Bias |Ea*| = %f\n",getEa(&pnbciv7,ID,OD,&civ));
  */

  /* TEST : get forwards PNBs */
  if (argc==5) getFPNBs(&fpnb,0.12,r,ID,OD,NULL);
  else getFPNBs(&fpnb,0.12,r,ID,OD,&civ);
  for (int i = 0; i < fpnb.len; i++) {
    printf("FPNB (%d,%d) ",fpnb.word[i],fpnb.bit[i]);
    printf("with neutrality measure = %f\n",fpnb.bias[i]);
  } printf("-> Found %d FPNBs and %d ivs\n",fpnb.len,fpnb.n);
  

  /* TEST: get bias E */
  //printf("Bias |E*| = %f\n",getE(&pnb_1_14,ID,OD,NULL));
  /*printf("Bias |E*| = %f conditioned \n",getE(&pnbciv7,ID,OD,&civ));
  */

  /* TEST : get bias Ed fixing non FPNB*/
  /*
  getFPSBs(&fpnb, &fpsb);
  for ( int i = 0; i < fpsb.n; ++i) printf("(%d,%d) is FPSB and IV\n",fpsb.word[i],fpsb.bit[i]);
  printf("%f = Ed with FPSBs with ID = [%d,%d]\n", getEd(ID,OD,NULL,&fpnb),ID[0],ID[1]);

  if (argc<=5)
    printf("%f = Ed FPNBs with ID = [%d,%d]\n", getEd(ID,OD,NULL,&fpnb),ID[0],ID[1]);
  else
    printf("%f = Ed FPNBs CIV x7 with ID = [%d,%d]\n", getEd(ID,OD,&civ,&fpnb),ID[0],ID[1]);
  
  if (argc==5) getBPNBs(&bpnb,0.12,ID,OD,NULL,&fpnb);
  else getBPNBs(&bpnb,0.12,ID,OD,&civ,&fpnb);
  for (int i = 0; i < bpnb.len; i++) {
    if ( iskey(bpnb.word[i]) ) printf("BPNKB (%d,%d) ",bpnb.word[i],bpnb.bit[i]);
    else printf("BPNB (%d,%d) ",bpnb.word[i],bpnb.bit[i]);
    printf("with neutrality measure = %f\n",bpnb.bias[i]);
  } printf("-> Found %d BPNBs and %d key using 2 CIV and FPNBs \n",bpnb.len,bpnb.n);
  */
  timeit(t);
  
}

