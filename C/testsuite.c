#include "merengue.c"
#include "parameters.h"
#include <time.h>
#include <sys/time.h>
#include <string.h>

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
  int ID[2], OD[2], R, r, eNk, eNiv;
  unsigned int Nk;
  unsigned long int Niv;
  float th;
  PNB *bpnb, *fpnb;
  CIV *civ;
  
  // Seed random generator
  struct timeval t1;
  gettimeofday(&t1, NULL);
  srand(t1.tv_usec * t1.tv_sec);

  printf("Compute: "); scanf("%s", stget);
  printf("Nk = 2^"); scanf("%d", &eNk);
  printf("Niv = 2^"); scanf("%d", &eNiv);
  printf("ID = "); scanf("%d %d", &ID[0], &ID[1]);
  printf("OD = "); scanf("%d %d", &OD[0], &OD[1]);
  printf("CIV = "); scanf("%s",stciv);
  printf("FPNB = "); scanf("%s",stfpnb);
  printf("r = "); scanf("%d", &r);

  if (strncmp(stget,"FPNB",8)!=0 && strncmp(stget,"Ef",8)!=0 ) {
    printf("R = "); scanf("%d", &R);
  } 
  if (strncmp(stget,"FPNB",8)==0 || strncmp(stget,"BPNB",8)==0 ) {
    printf("th = "); scanf("%f",&th);
  }

  if      (strncmp(stciv,"NULL",8)==0) { civ = NULL; }
  else if (strncmp(stciv,"civx7",8)==0) { civ = &civx7; }
  else if (strncmp(stciv,"civx8",8)==0) { civ = &civx8; }
  else if (strncmp(stciv,"civx7x6",8)==0) { civ = &civx7x6; }

  if      (strncmp(stfpnb,"NULL",8)==0) { fpnb = NULL; }
  else if (strncmp(stfpnb,"fpnb_",8)==0) { fpnb = &fpnb_; }
  else if (strncmp(stfpnb,"fpnbx7",8)==0) { fpnb = &fpnbx7; }
  else if (strncmp(stfpnb,"fpnbx7x6",8)==0) { fpnb = &fpnbx7x6; } 

  if (strncmp(stget,"Eg",8)==0 || strncmp(stget,"E",8)==0) {
    if      (civ==NULL && R==8 && r==4 ){
      if      ( OD[0]==1 && OD[1]==14 ) { bpnb = &bpnb_1_14; }
      else if ( OD[0]==6 && OD[1]==14 ) { bpnb = &bpnb_6_14; }}
    else if (civ==&civx7) { bpnb = &bpnbx7; }
    else if (civ==&civx8) { bpnb = &bpnbx8; }
    else if (civ==&civx7x6) { bpnb = &bpnbx7x6; }
  }

  Nk = pow(2,eNk); Niv = pow(2,eNiv);    

  clock_t t;
  t = clock(); 

  if (strncmp(stget,"Ef",8)==0) {
    printf("|Ef*|=%f",getEf(r,ID,OD,civ,fpnb,Nk,Niv));
  } 
  else if (strncmp(stget,"Eg",8)==0) {
    printf("|Eg*|=%f",getEg(bpnb,ID,OD,civ,fpnb,Nk,Niv,R,r));
  } 
  else if (strncmp(stget,"E",8)==0) {
    printf("|E*|=%f",getE(bpnb,ID,OD,civ,fpnb,Nk,Niv,R,r));
  } 
  else if (strncmp(stget,"FPNB",8)==0) {
    getFPNBs(fpnb,th,r,ID,OD,civ,Nk,Niv);
    for (int i = 0; i < fpnb->len; ++i) {
      printf("FPNB (%d,%d) ",fpnb->word[i],fpnb->bit[i]);
      printf("with neutrality measure = %f\n",fpnb->bias[i]);
    } printf("-> %d FPNBs (%d iv)\n",fpnb->len,fpnb->n);   
  } 
  else if (strncmp(stget,"BPNB",8)==0) {
    getBPNBs(bpnb,R,r,th,ID,OD,civ,fpnb,Nk,Niv);
    for (int i = 0; i < bpnb->len; ++i) {
      if ( iskey(bpnb->word[i]) ) printf("BPNKB (%d,%d) ",bpnb->word[i],bpnb->bit[i]);
      else printf("BPNB (%d,%d) ",bpnb->word[i],bpnb->bit[i]);
      printf("with neutrality measure = %f\n",bpnb->bias[i]);
    } printf("-> %d BPNBs (%d key)\n",bpnb->len,bpnb->n);
  }
      
  t = timeit(t);

  return 0;

}

