#include "merengue.c"
#include "parameters.h"
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

static clock_t timeit(clock_t t){
  t = clock() - t;
  double time_taken = ((double)t)/CLOCKS_PER_SEC;
  printf("\n%f seconds\n",time_taken);
  return t;
}

static int max(int a, int b) { return a*(a>=b)+b*(b>a); }

static int same(char *str, char *txt){
  if (strlen(str)!=strlen(txt)) return 0;
  int ret = 1;
  int maxlen = max(strlen(str),strlen(txt));
  for (int c = 0 ; c < maxlen ; c++)
    ret = ret * str[c]==txt[c];
  return ret;
}

int main(int argc, char *argv[])
{ 
  char stflg[8], stret[100], stciv[4];
  PNB pnb = { 0, 0, {0}, {0}, {0.0}};
  int eNk, eNiv, flg, nfpnb, nbpnb;
  float ret;

  bpnb = NULL; fpnb = NULL; civ = NULL; 
  
  // Seed random generator
  struct timeval seed;
  gettimeofday(&seed, NULL);
  srand(seed.tv_usec * seed.tv_sec);
  
  // Disable stdout buffering
  setbuf(stdout, NULL);

  printf("Compute: "); scanf("%s", stflg);

  if      ( same(stflg,"Ef") )   { flg = flgEf; } 
  else if ( same(stflg,"Eg") )   { flg = flgEg; }
  else if ( same(stflg,"E") )    { flg = flgE; }
  else if ( same(stflg,"FPNB") ) { flg = flgF; }
  else if ( same(stflg,"BPNB") ) { flg = flgB; }
  else { printf("Ef | Eg | E | FPNB | BPNB\n"); exit(1); }

  printf("Nk = 2^");  scanf("%d", &eNk ); Nk  = pow(2,eNk );
  printf("Niv = 2^"); scanf("%d", &eNiv); Niv = pow(2,eNiv);
  
  printf("ID = "); scanf("%d %d", &ID[0], &ID[1]);
  if ( ID[0]<0 || ID[0]>15 || ID[1]<0 || ID[1]>31 ) { exit(1); }

  printf("#ODs: "); scanf("%d",&nODs); if ( nODs < 1 || nODs > 256 ) { exit(1);}
  for (int o = 0; o < nODs; o++) {
    printf("OD[%d] = ",o); scanf("%d %d", &OD[o][0], &OD[o][1]);
    if ( OD[o][0]<0 || OD[o][0]>15 || OD[o][1]<0 || OD[o][1]>31 ) { exit(1); }
  }
  
  printf("r = "); scanf("%d", &r);
  if ( flg!=flgF && flg!=flgEf ) {
    printf("R = "); scanf("%d", &R);
    if ( R < r ) { exit(1); }
  } else { R = r; }

  printf("CIV = "); scanf("%s",stciv);
  if      ( same(stciv,"7") ) { civ = &civx7; }
  else if ( same(stciv,"8") ) { civ = &civx8; }
  else if ( same(stciv,"6") ) { civ = &civx7x6; }
  else if ( !same(stciv,"n") ) { exit(1);}

  if (flg!=flgF) {
    printf("FPNB = "); scanf("%d",&nfpnb);
    if      ( nfpnb == 11 ) { fpnb = &fpnb11; }
    else if ( civ == NULL && nfpnb == 47 ) { fpnb = &fpnb47; }
    else if ( civ == &civx7 && nfpnb == 52 ) { fpnb = &fpnb52; }
    else if ( civ == &civx7x6 && nfpnb == 54 ) { fpnb = &fpnb54; }
    else { printf("11 | 47 | 52 | 54\n"); exit(1); }
  }

  if ( flg==flgF || flg==flgB ) {
    printf("th = "); scanf("%f",&th);
  }

  if ( flg==flgEg || flg==flgE ) {
    printf("BPNB = "); scanf("%d",&nbpnb);
    if ( R == 8 && r == 4 ) {
      if      ( nbpnb == 48 ) { bpnb = &bpnb48; }
      else if ( civ==NULL && OD[0][0]==1 && OD[0][1]==14 && nbpnb==36) { bpnb = &bpnb36; }
      else if ( civ==NULL && OD[0][0]==6 && OD[0][1]==14 && nbpnb==30) { bpnb = &bpnb30; }
      else if ( civ==&civx7 && nbpnb==37 ) { bpnb = &bpnb37; }
      else if ( civ==&civx7 && nbpnb==33 ) { bpnb = &bpnb33; }
      else if ( civ==&civx8 && nbpnb==31 ) { bpnb = &bpnb31; }
      else if ( civ==&civx7x6 && nbpnb==38 ) { bpnb = &bpnb38; }
    } else if ( R==9 && r==5 && OD[0][0]==9 && OD[0][1]==21 && nbpnb==25) { bpnb = &bpnb25;
    } else if ( R==8 && r==5 && nODs==3 && nbpnb==36) { bpnb = &bpnb_1_13; }
  }

  clock_t t;
  t = clock(); 

  if ( flg == flgEf || flg == flgEg || flg == flgE ) {ret = getbias(flg);}

  if    ( flg == flgEf ) {
    printf("|Ef*| = %f = 2^(%.1f)\n",ret,log2(ret));
  } 
  else if ( flg == flgEg ) {
    printf("|Eg*| = %f = 2^(%.1f)\n",ret,log2(ret));
  } 
  else if ( flg == flgE ) {
    printf("|E*| = %f = 2^(%.1f)\n",ret,log2(ret));
  } 
  else if ( flg == flgF ) {
    getPNBs(flg,&pnb);
    for (int i = 0; i < pnb.len; ++i) {
      printf("FPNB (%d,%d) ",pnb.word[i],pnb.bit[i]);
      printf("with neutrality measure = %f\n",pnb.bias[i]);
    } printf("-> %d FPNBs (%d iv)\n",pnb.len,pnb.n);   
  } 
  else if ( flg == flgB ) {
    getPNBs(flg,&pnb);
    for (int i = 0; i < pnb.len; ++i) {
      printf("BPNB (%d,%d) ",pnb.word[i],pnb.bit[i]);
      printf("with neutrality measure = %f\n",pnb.bias[i]);
    } printf("-> %d BPNBs (%d key)\n",pnb.len,pnb.n);
  }

  t = timeit(t);

  return 0;

}

