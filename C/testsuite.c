#include "merengue.c"
#include <time.h>
#include <sys/time.h>

/*
gcc testsuite.c -o testsuite -lm
./testsuite {ID word} {ID bit} {OD word} {OD bit}
*/

/* Recommended arguments:
   ID = {7,31} OD = {1,31}
   ID = {8,31} OD = {6,31}
   ID = {7,31} OD = {1,31} CIV1 = {7,19,12}
   ID = {8,31} OD = {6,31} CIV1 = {8,19,12}
   ID = {7,31} OD = {1,31} CIV1 = {7,19,12} CIV2 = {6,25,10}
 */

int ID[2], OD[2], R, r, eNk, eNiv;
unsigned int Nk, Niv;
float th;
PNB bpnb, fpnb, fpsb;
CIV civ;

CIV civx7 = { 1, {7}, {19}, {12} };
CIV civx8 = { 1, {8}, {19}, {12} };
CIV civx7x6 = { 2, {7, 6}, {19, 25}, {12, 10} };

PNB bpnb_1_14 = { 63, 36,
    {  1, 1, 1, 1, 1, 1, 3, 3, 4, 4, 4,11,
      12,12,12,12,12,12,12,12,12,12,12,12,
      12,13,13,13,14,14,14,14,14,14,14,14 },
    { 26,27,28,29,30,31, 7, 8,24,25,26,20,
       5, 6, 7, 8, 9,10,11,12,13,14,15,16,
      17,18,19,20, 0, 1,18,19,20,21,22,23 } };
PNB bpnb_6_14 = { 63, 30,
  { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 
    2, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4,11,11,11,12},
  { 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,
    20, 0, 1,18,19,20,21,22,23,7, 8,13,14,31,20} };
PNB bpnbx7 = { 68, 37,
  {  1, 1, 1, 1, 1, 1, 3, 3, 4, 4, 4,11,
    12,12,12,12,12,12,12,12,12,12,12,12,12,
    13,13,13,14,14,14,14,14,14,14,14,14 },
  { 26,27,28,29,30,31, 7, 8,24,25,26,20,
    5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,
    18,19,20, 0, 1, 2,18,19,20,21,22,23 } };
PNB bpnbx7x6 = { 69, 38,
  {  1, 1, 1, 1, 1, 1, 3, 3, 4, 4, 4, 4,11,
    12,12,12,12,12,12,12,12,12,12,12,12,12,
    13,13,13,14,14,14,14,14,14,14,14,14 },
  { 26,27,28,29,30,31, 7, 8,12,24,25,26,20,
    5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,
    18,19,20, 0, 1, 2,18,19,20,21,22,23 } };
PNB fpnb_ = { 134, 47,
  {  6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
     7, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
     8, 8, 8, 8, 9, 9, 9, 9, 9, 9, 9, 9,
     9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9 },
  {  0,10,22,23,24,25,26,27,28,29,30,31,
    18,19, 1, 2, 3,11,12,23,24,25,26,27,
    28,29,30,31, 5, 6, 7, 8, 9,10,11,13,
    14,15,16,17,24,25,27,28,29,30,31 } };
PNB fpnbx7 = { 139, 52,
  {  6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
     7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
     8, 8, 8, 8, 9, 9, 9, 9, 9, 9, 9, 9,
     9, 9, 9, 9, 9,9, 9, 9, 9, 9, 9, 9 },
  {  0,10,11,22,23,24,25,26,27,28,29,30,31,
     1,18,19,30, 1, 2, 3,10,11,12,23,24,25,26,27,
    28,29,30,31, 5, 6, 7, 8, 9,10,11,13,
    14,15,16,17,18,24,25,27,28,29,30,31 } };
PNB fpnbx7x6 = { 141, 54,
  {  6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
     7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
     8, 8, 8, 8, 9, 9, 9, 9, 9, 9, 9, 9,
     9, 9, 9, 9, 9,9, 9, 9, 9, 9, 9, 9 },
  {  0,10,11,22,23,24,25,26,27,28,29,30,31,
     1,18,19,20,30, 1, 2, 3,10,11,12,13,23,24,25,26,27,
    28,29,30,31, 5, 6, 7, 8, 9,10,11,13,
    14,15,16,17,18,24,25,27,28,29,30,31 } };


static clock_t timeit(clock_t t){
  t = clock() - t;
  double time_taken = ((double)t)/CLOCKS_PER_SEC;
  printf("\n%f seconds\n",time_taken);
  return t;
}

int main(int argc, char *argv[])
{ 
  /*
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
  */

  // Seed random generator
  struct timeval t1;
  gettimeofday(&t1, NULL);
  srand(t1.tv_usec * t1.tv_sec); //srand(time(NULL));
  
  clock_t t;
  t = clock(); 

  /* TEST : Salsa and Alsa working */
  /*
  ECRYPT_ctx X, Z, Y;
  bitstring(k,32);
  ECRYPT_keysetup(&X,k,256,64);
  bitstring(v,16);
  ECRYPT_ivsetup(&X,v);
  salsa(4,Z.state,X.state,0);
  prettyprint(&Z);
  salsa(8,Z.state,X.state,0);
  aslas(8,4,Y.state,Z.state);
  prettyprint(&Y);
  salsa(5,Z.state,X.state,0);
  prettyprint(&Z);
  salsa(9,Z.state,X.state,0);
  aslas(9,5,Y.state,Z.state);
  prettyprint(&Y);
  */

  /* TEST : get Ed -> 277 seconds */
  /*
  r = 4; Nk = pow(2,12); Niv = pow(2,12);
  ID[0] = 7; ID[1] = 31; OD[0] = 1; OD[1] = 14;
  printf("|Ed*| = %f : ID = [%d,%d]\n", 
    getEf(r,ID,OD,NULL,NULL,Nk,Niv),ID[0],ID[1]);
  printf("|Ed*| = %f : ID = [%d,%d] : CIV X7_12\n", 
    getEf(r,ID,OD,&civx7,NULL,Nk,Niv),ID[0],ID[1]);
  printf("|Ed*| = %f : ID = [%d,%d] : CIV X7_12 X6_10\n", 
    getEf(r,ID,OD,&civx7x6,NULL,Nk,Niv),ID[0],ID[1]);

  printf("|Ed*| = %f : ID = [%d,%d] :                 : FPNB 47\n", 
    getEf(r,ID,OD,NULL,&fpnb_,Nk,Niv),ID[0],ID[1]);
  printf("|Ed*| = %f : ID = [%d,%d] : CIV X7_12       : FPNB 52\n", 
    getEf(r,ID,OD,&civx7,&fpnbx7,Nk,Niv),ID[0],ID[1]);
  printf("|Ed*| = %f : ID = [%d,%d] : CIV X7_12 X6_10 : FPNB 54\n", 
    getEf(r,ID,OD,&civx7x6,&fpnbx7x6,Nk,Niv),ID[0],ID[1]);

  ID[0] = 8; ID[1] = 31; OD[0] = 6; OD[1] = 14;
  printf("|Ed*| = %f : ID = [%d,%d]\n", 
    getEf(r,ID,OD,NULL,NULL,Nk,Niv),ID[0],ID[1]);
  printf("|Ed*| = %f : ID = [%d,%d] : CIV X7_12\n", 
    getEf(r,ID,OD,&civx8,NULL,Nk,Niv),ID[0],ID[1]);
  */

  /* TEST : get FPNBs */
  /*
  r = 4; th = 0.12; Nk = pow(2,10); Niv = pow(2,14);
  ID[0] = 7; ID[1] = 31; OD[0] = 1; OD[1] = 14;

  getFPNBs(&fpnb,0.12,r,ID,OD,NULL,Nk,Niv);             // 30m
  for (int i = 0; i < fpnb.len; i++) {
    printf("FPNB (%d,%d) with neutrality measure = %f\n",
      fpnb.word[i],fpnb.bit[i],fpnb.bias[i]);
  } printf("-> %d FPNBs (%d IVs) : ID=[%d,%d]\n",
    fpnb.len,fpnb.n,ID[0],ID[1]);

  getFPNBs(&fpnb,0.12,r,ID,OD,&civx7,Nk,Niv);           // 1h
  for (int i = 0; i < fpnb.len; i++) {
    printf("FPNB (%d,%d) with neutrality measure = %f\n",
      fpnb.word[i],fpnb.bit[i],fpnb.bias[i]);
  } printf("-> %d FPNBs (%d IVs) : ID=[%d,%d] : CIV X7_12\n",
    fpnb.len,fpnb.n,ID[0],ID[1]);

  getFPNBs(&fpnb,0.12,r,ID,OD,&civx7x6,Nk,Niv);         // 2h
  for (int i = 0; i < fpnb.len; i++) {
    printf("FPNB (%d,%d) with neutrality measure = %f\n",
      fpnb.word[i],fpnb.bit[i],fpnb.bias[i]);
  } printf("-> %d FPNBs (%d IVs) : ID=[%d,%d] : CIV X7_12 X6_10\n",
    fpnb.len,fpnb.n,ID[0],ID[1]);
  */
  


  /* TEST 2: get BPNBs */
  
  R = 9; r = 5; th = 0.12; eNk = 10 ; eNiv = 14; Nk = pow(2,eNk); Niv = pow(2,eNiv);
  ID[0] = 7; ID[1] = 31; OD[0] = 1; OD[1] = 14;

  getBPNBs(&bpnb,R,r,th,ID,OD,NULL,NULL,Nk,Niv);                      // 1h 30m
  for (int i = 0; i < bpnb.len; i++) {
    if ( iskey(bpnb.word[i]) ) printf("BPNKB (%d,%d) ",bpnb.word[i],bpnb.bit[i]);
    else printf("BPNB (%d,%d) ",bpnb.word[i],bpnb.bit[i]);
    printf("with neutrality measure = %f\n",bpnb.bias[i]);
  } printf("-> %d BPNBs (%d key) : ID=[%d,%d] : Nk 2^%d Niv 2^%d\n",
    bpnb.len,bpnb.n,ID[0],ID[1],eNk,eNiv);
  /*
  getBPNBs(&bpnb,0.12,ID,OD,&civx7);                    // 3h 30m
  for (int i = 0; i < bpnb.len; i++) {
    if ( iskey(bpnb.word[i]) ) printf("BPNKB (%d,%d) ",bpnb.word[i],bpnb.bit[i]);
    else printf("BPNB (%d,%d) with neutrality measure = %f\n",
      bpnb.word[i],bpnb.bit[i],bpnb.bias[i]);
  } printf("-> Found %d BPNBs (%d key) : ID=[%d,%d] : CIV X7_12\n",
    bpnb.len,bpnb.n,ID[0],ID[1]);
  
  getBPNBs(&bpnb,R,r,th,ID,OD,&civx7x6,NULL,Nk,Niv);                  // 7h
  for (int i = 0; i < bpnb.len; i++) {
    if ( iskey(bpnb.word[i]) ) printf("BPNKB (%d,%d) ",bpnb.word[i],bpnb.bit[i]);
    else printf("BPNB (%d,%d) ",bpnb.word[i],bpnb.bit[i]);
    printf("with neutrality measure = %f\n",bpnb.bias[i]);
  } printf("-> Found %d BPNBs (%d key) : ID=[%d,%d] : R=%d r=%d : Nk 2^%d Niv^%d : CIV X7_12 X6_10\n",
    bpnb.len,bpnb.n,ID[0],ID[1],R,r,eNk,eNiv);
  */



  /* TEST : get Ea */

  r = 4; eNk = 8. ; eNiv = 25. ; Nk = pow(2,eNk); Niv = pow(2,eNiv);
  ID[0] = 7; ID[1] = 31; OD[0] = 1; OD[1] = 14; 
  
  /*                                                             
  printf("|Ea*| = %f : Nk 2^%d Niv 2^%d : BPNB 36\n",             // 1d 5h  2^10 2^26? -> 0.002266 :: 1h 2^10 2^21 -> 0.002357
    getEg(&bpnb_1_14,ID,OD,NULL,NULL,Nk,Niv),eNk,eNiv);           // 1h 15m -> 0.0015625 2^20 2^11 :: // 7 h -> 0.006348 Nk 20 Niv 14
  
  printf("|Ea*| = %f : Nk 2^%d Niv 2^%d : BPNB 37 : CIV X7_12\n",     // 2h 20m -> 0.002007 2^10 2^21
    getEg(&bpnbx7,ID,OD,&civx7,NULL,Nk,Niv),eNk,eNiv);                // 2d 21h 2^36 -> 0.001496
  
  printf("|Ea*| = %f : Nk 2^%d Niv 2^%d : BPNB 38 : CIV X7_12 X6_10\n",  // 4h 30m 2^20 2^11 -> 0.031250 :: 
    getEg(&bpnbx7x6,ID,OD,&civx7x6,NULL,Nk,Niv),eNk,eNiv);               // 7h     2^10 2^21 -> 0.001873 :: 5h 2^10 2^21 -> 0.001811
  
  printf("|Ea*| = %f : Nk 2^%d Niv 2^%d : BPNB 36 :                 : FPNB 47\n",  // 3h 30m 2^10 2^21 -> 0.0025
    getEg(&bpnb_1_14,ID,OD,NULL,&fpnb_,Nk,Niv),eNk,eNiv);                         
  
  printf("|Ea*| = %f : Nk 2^%d Niv 2^%d : BPNB 37 : CIV X7_12       : FPNB 52\n",  // 5h 2^10 2^21 -> 0.001969
    getEg(&bpnbx7,ID,OD,&civx7,&fpnbx7,Nk,Niv),eNk,eNiv);
  
  printf("|Ea*| = %f : Nk 2^%d Niv 2^%d : BPNB 38 : CIV X7_12 X6_10 : FPNB 54\n",       
    getEg(&bpnbx7x6,ID,OD,&civx7x6,&fpnbx7x6,Nk,Niv),eNk,eNiv);                          // 5 h -> 0.031250 2^20 2^11
  */


  /* TEST: get bias E */
  
  r = 4; eNk = 4 ; eNiv = 21 ; Nk = pow(2,eNk); Niv = pow(2,eNiv);
  ID[0] = 7; ID[1] = 31; OD[0] = 1; OD[1] = 14; 
  /*
  printf("|E*| = %f : BPNB 36 :\n",                           // 1k 2^31 -> 0.005714 1h
    getE(&bpnb_1_14,ID,OD,NULL,NULL,Nk,Niv));                 // 2^31 -> 0.015 40m

  printf("|E*| = %f : BPNB 37 : CIV X7_12\n",
    getE(&bpnbx7,ID,OD,&civx7,NULL,Nk,Niv));                    // 2d 21h
  printf("|E*| = %f : BPNB 38 : CIV X7_12 X6_10\n",
    getE(&bpnbx7,ID,OD,&civx7,NULL,Nk,Niv));
  printf("|E*| = %f : BPNB 36 :                 : FPNB 47\n",
    getE(&bpnb_1_14,ID,OD,NULL,&fpnb_,Nk,Niv));
  printf("|E*| = %f : BPNB 37 : CIV X7_12       : FPNB 52\n",
    getE(&bpnbx7,ID,OD,&civx7,&fpnbx7,Nk,Niv));
  
  printf("|E*| = %f : BPNB 38 : CIV X7_12 X6_10 : FPNB 54\n",
    getE(&bpnbx7,ID,OD,&civx7,&fpnbx7x6,Nk,Niv));
  */

  t = timeit(t);
  
}

