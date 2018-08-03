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

int main(int argc, char *argv[])
{ 
  char stget[8], stret[100], stciv[8], stfpnb[8];
  PNB pnb = { 0, 0, {0}, {0}, {0.0}};
  int eNk, eNiv;

  bpnb = NULL; fpnb = NULL; civ = NULL; 
  
  // Seed random generator
  struct timeval seed;
  gettimeofday(&seed, NULL);
  srand(seed.tv_usec * seed.tv_sec);

  printf("Compute: "); scanf("%s", stget);

  if ( !same(stget,"Ef")   
    && !same(stget,"Eg") 
    && !same(stget,"E")  
    && !same(stget,"FPNB") 
    && !same(stget,"BPNB") ) {
    printf("Ef | Eg | E | FPNB | BPNB\n");
    exit(1);
  }

  printf("Nk = 2^");  scanf("%d", &eNk ); Nk  = pow(2,eNk );
  printf("Niv = 2^"); scanf("%d", &eNiv); Niv = pow(2,eNiv);
  printf("ID = "); scanf("%d %d", &ID[0], &ID[1]);
  printf("OD = "); scanf("%d %d", &OD[0], &OD[1]);

  if ( ID[0]<0  || ID[1]<0   || OD[0]<0  || OD[0]<0
    || ID[0]>15 || ID[1]>31  || OD[0]>15 || OD[0]>31 ) { exit(1); }
  
  printf("r = "); scanf("%d", &r);

  if ( !same(stget,"FPNB") && !same(stget,"Ef") ) {
    printf("R = "); scanf("%d", &R);
    if ( R < r ) { exit(1); }
  } else { R = r; }

  printf("CIV = "); scanf("%s",stciv);  

  if      ( same(stciv,"7") ) { civ = &civx7; }
  else if ( same(stciv,"8") ) { civ = &civx8; }
  else if ( same(stciv,"6") ) { civ = &civx7x6; }

  if (!same(stget,"FPNB")) {
    printf("FPNB = "); scanf("%s",stfpnb);
    if (stfpnb[0]!='n'){
      if      ( civ == NULL ) { fpnb = &fpnb_; }
      else if ( civ == &civx7 ) { fpnb = &fpnbx7; }
      else if ( civ == &civx7x6 ) { fpnb = &fpnbx7x6; }
    }
  }

  if (same(stget,"FPNB") || same(stget,"BPNB") ) {
    printf("th = "); scanf("%f",&th);
  }

  if (same(stget,"Eg") || same(stget,"E") ) {
    if ( R == 8 && r == 4 ) {
      if      ( civ==NULL && OD[0]==1 && OD[1]==14 ) { bpnb = &bpnb_1_14; }
      else if ( civ==NULL && OD[0]==6 && OD[1]==14 ) { bpnb = &bpnb_6_14; }
      else if ( civ==&civx7 ) { bpnb = &bpnbx7; }
      else if ( civ==&civx8 ) { bpnb = &bpnbx8; }
      else if ( civ==&civx7x6 ) { bpnb = &bpnbx7x6; }
    } else if ( R==9 && r==5 && OD[0]==9 && OD[1]==21 ) { bpnb = &bpnb_9_21; }
  }

  clock_t t;
  t = clock(); 

  if (same(stget,"Ef")) {
    printf("|Ef*|=%f",getEf(stget));
  } 
  else if (same(stget,"Eg")) {
    printf("|Eg*|=%f",getEg());
  } 
  else if (same(stget,"E")) {
    printf("|E*|=%f",getE());
  } 
  else if (same(stget,"FPNB")) {
    getPNBs('f',&pnb);
    for (int i = 0; i < pnb.len; ++i) {
      printf("FPNB (%d,%d) ",pnb.word[i],pnb.bit[i]);
      printf("with neutrality measure = %f\n",pnb.bias[i]);
    } printf("-> %d FPNBs (%d iv)\n",pnb.len,pnb.n);   
  } 
  else if (same(stget,"BPNB")) {
    getPNBs('b',&pnb);
    for (int i = 0; i < pnb.len; ++i) {
      printf("BPNB (%d,%d) ",pnb.word[i],pnb.bit[i]);
      printf("with neutrality measure = %f\n",pnb.bias[i]);
    } printf("-> %d BPNBs (%d key)\n",pnb.len,pnb.n);
  }

  t = timeit(t);

  return 0;

}

