#include "ecrypt-sync.h"
#include "merengue.c"
#include <stdio.h>
#include <time.h>

/*
gcc -c ecrypt.c -o ecrypt.o
gcc -c testsuite.c -o testsuite.o
gcc -c testsuite.o -o test -lcrypto
./test
*/

u8 m[64], c[64], d[64];
int ID[2], OD[2], CIV[3];
int Nk, Niv;

int main()
{  
  clock_t t;
  t = clock();
  Nk = pow(2,10);
  Niv = pow(2,12);
  
  /* TEST 1: get bias Ed with Aumasson's ( /\4_1,14 | /\0_7,31 ) */
  ID[0] = 7; ID[1] = 31;
  OD[0] = 1; OD[1] = 14;
  CIV[0] = 7; CIV[1] = 19; CIV[2] = 12;
  printf("%f = Ed with ID = [%d,%d]\n", Ed(ID,OD,Nk,Niv,NULL),ID[0],ID[1]);
  printf("%f = conditioned Ed with ID = [%d,%d]\n", Ed(ID,OD,Nk,Niv,CIV),ID[0],ID[1]);

  /* TEST 2: get bias Ed with our differential ( /\6_1,14 | /\0_8,31 ) */
  ID[0] = 8; ID[1] = 31;
  OD[0] = 6; OD[1] = 14;
  CIV[0] = 8; CIV[1] = 19; CIV[2] = 12;
  printf("%f = Ed with ID = [%d,%d]\n", Ed(ID,OD,Nk,Niv,NULL),ID[0],ID[1]);
  printf("%f = conditioned Ed with ID = [%d,%d]\n", Ed(ID,OD,Nk,Niv,CIV),ID[0],ID[1]);
  
  t = clock() - t;
  double time_taken = ((double)t)/CLOCKS_PER_SEC;
  printf("\n%f seconds\n",time_taken);
}

