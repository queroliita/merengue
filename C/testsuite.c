#include "merengue.c"
#include <stdio.h>
#include <time.h>

/*
gcc testsuite.c testsuite
./testsuite
*/

int ID[2], OD[2], CIV[3];
PNB_list PNB;

int main()
{  
  clock_t t;
  t = clock();
  
  /* TEST 1: get bias Ed with Aumasson's ( /\4_1,14 | /\0_7,31 ) *//*
  ID[0] = 7; ID[1] = 31;
  OD[0] = 1; OD[1] = 14;
  CIV[0] = 7; CIV[1] = 19; CIV[2] = 12;
  printf("%f = Ed with ID = [%d,%d]\n", getEd(ID,OD,NULL),ID[0],ID[1]);
  printf("%f = conditioned Ed with ID = [%d,%d]\n", getEd(ID,OD,CIV),ID[0],ID[1]);
  *//* 6.5 seconds */
  
  /* TEST 2: compute PNBs with Aumasson's differential *//*
  getPNBs(&PNB,0.12,ID,OD,NULL);
  for (int i = 0; i < PNB.n; i++){
    printf("PNB (%d,%d) ",PNB.word[i],PNB.bit[i]);
    printf("with neutrality measure = %f\n",PNB.bias[i]);
  } printf("-> Found %d PNBs\n",PNB.n);
  *//* 375 seconds */

  /* TEST 3: get bias Ed with our differential ( /\6_1,14 | /\0_8,31 ) *//*
  ID[0] = 8; ID[1] = 31;
  OD[0] = 6; OD[1] = 14;
  CIV[0] = 8; CIV[1] = 19; CIV[2] = 12;
  printf("%f = Ed with ID = [%d,%d]\n", getEd(ID,OD,NULL),ID[0],ID[1]);
  printf("%f = conditioned Ed with ID = [%d,%d]\n", getEd(ID,OD,CIV),ID[0],ID[1]);
  *//* 6.5 seconds */
  
  /* TEST 4: compute PNBs with our differential */
  getPNBs(&PNB,0.12,ID,OD,NULL);
  for (int i = 0; i < PNB.n; i++){
    printf("PNB (%d,%d) ",PNB.word[i],PNB.bit[i]);
    printf("with neutrality measure = %f\n",PNB.bias[i]);
  } printf("-> Found %d PNBs\n",PNB.n);
  /* 375 seconds */

  t = clock() - t;
  double time_taken = ((double)t)/CLOCKS_PER_SEC;
  printf("\n%f seconds\n",time_taken);

}

