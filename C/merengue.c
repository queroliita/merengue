#include "ecrypt-sync.h"
#include "ecrypt.c"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

u8 k[32];
u8 v[16];

/* Modified Salsa20/R function */
static void salsa(int R, u32 output[16], const u32 input[16]) {
  u32 x[16];
  int i;
  for (i = 0;i < 16;++i) x[i] = input[i];
  for (i = R;i > 0;i -= 2) {
    x[ 4] = XOR(x[ 4],ROTATE(PLUS(x[ 0],x[12]), 7));
    x[ 8] = XOR(x[ 8],ROTATE(PLUS(x[ 4],x[ 0]), 9));
    x[12] = XOR(x[12],ROTATE(PLUS(x[ 8],x[ 4]),13));
    x[ 0] = XOR(x[ 0],ROTATE(PLUS(x[12],x[ 8]),18));
    x[ 9] = XOR(x[ 9],ROTATE(PLUS(x[ 5],x[ 1]), 7));
    x[13] = XOR(x[13],ROTATE(PLUS(x[ 9],x[ 5]), 9));
    x[ 1] = XOR(x[ 1],ROTATE(PLUS(x[13],x[ 9]),13));
    x[ 5] = XOR(x[ 5],ROTATE(PLUS(x[ 1],x[13]),18));
    x[14] = XOR(x[14],ROTATE(PLUS(x[10],x[ 6]), 7));
    x[ 2] = XOR(x[ 2],ROTATE(PLUS(x[14],x[10]), 9));
    x[ 6] = XOR(x[ 6],ROTATE(PLUS(x[ 2],x[14]),13));
    x[10] = XOR(x[10],ROTATE(PLUS(x[ 6],x[ 2]),18));
    x[ 3] = XOR(x[ 3],ROTATE(PLUS(x[15],x[11]), 7));
    x[ 7] = XOR(x[ 7],ROTATE(PLUS(x[ 3],x[15]), 9));
    x[11] = XOR(x[11],ROTATE(PLUS(x[ 7],x[ 3]),13));
    x[15] = XOR(x[15],ROTATE(PLUS(x[11],x[ 7]),18));
    if  ((i%2)==0){
      x[ 1] = XOR(x[ 1],ROTATE(PLUS(x[ 0],x[ 3]), 7));
      x[ 2] = XOR(x[ 2],ROTATE(PLUS(x[ 1],x[ 0]), 9));
      x[ 3] = XOR(x[ 3],ROTATE(PLUS(x[ 2],x[ 1]),13));
      x[ 0] = XOR(x[ 0],ROTATE(PLUS(x[ 3],x[ 2]),18));
      x[ 6] = XOR(x[ 6],ROTATE(PLUS(x[ 5],x[ 4]), 7));
      x[ 7] = XOR(x[ 7],ROTATE(PLUS(x[ 6],x[ 5]), 9));
      x[ 4] = XOR(x[ 4],ROTATE(PLUS(x[ 7],x[ 6]),13));
      x[ 5] = XOR(x[ 5],ROTATE(PLUS(x[ 4],x[ 7]),18));
      x[11] = XOR(x[11],ROTATE(PLUS(x[10],x[ 9]), 7));
      x[ 8] = XOR(x[ 8],ROTATE(PLUS(x[11],x[10]), 9));
      x[ 9] = XOR(x[ 9],ROTATE(PLUS(x[ 8],x[11]),13));
      x[10] = XOR(x[10],ROTATE(PLUS(x[ 9],x[ 8]),18));
      x[12] = XOR(x[12],ROTATE(PLUS(x[15],x[14]), 7));
      x[13] = XOR(x[13],ROTATE(PLUS(x[12],x[15]), 9));
      x[14] = XOR(x[14],ROTATE(PLUS(x[13],x[12]),13));
      x[15] = XOR(x[15],ROTATE(PLUS(x[14],x[13]),18));
    }
  }
  for (i = 0;i < 16;i++) output[i] = x[i];
}

/* Store random bitstring of given byte length */
void bitstring(u8 *array, int bytes) {
  for (int i = 0; i < bytes; i++) array[i] = rand();
}

/* Compare two values */
int compare(const void* a, const void* b) {
     int int_a = * ( (int*) a );
     int int_b = * ( (int*) b );
     if ( int_a == int_b ) return 0;
     else if ( int_a < int_b ) return -1;
     else return 1;
}

/* Set X[word][bit] = 0 */
void set0(ECRYPT_ctx *X, int word, int bit) {
  X->state[word] = X->state[word] & ~ ( 1 << bit );
}

/* Set X[word][bit] = 1 */
void set1(ECRYPT_ctx *X, int word, int bit){
  X->state[word] = (X->state[word] & ~(1<<bit)) | (1<<bit);
}

/* Return bit value at X[word][bit] */
int getbit(ECRYPT_ctx *X, int word, int bit){
  return ( X->state[word] >> bit ) & 1;
}

/* Count number of bit differences between two states */
int numdif(ECRYPT_ctx Z0, ECRYPT_ctx Z1){
  int notequal = 0;
  for (int i = 0; i < 16; i++)
    for (int j = 0; j < 32; j++)
      notequal += ((Z0.state[i] >> j) & 1) ^ ((Z1.state[i] >> j) & 1);
  return notequal;
}

/* Print matrix state in hexadecimal format */
void prettyprint(ECRYPT_ctx *X) {
  for (int i = 0; i < 16; i++) {
    if ( i % 4 == 0 ) printf("\n");
    printf("%02x ",X->state[i]);
  } printf("\n");
}

/* Set input difference at initial states */
void setID(ECRYPT_ctx *X0, ECRYPT_ctx *X1, int *ID) {
  set0(X0,ID[0],ID[1]); 
  set1(X1,ID[0],ID[1]); 
}

/* Fix conditioned bits of IV */
void fix(ECRYPT_ctx *X0, ECRYPT_ctx *X1, ECRYPT_ctx *aux, int *CIV){
  for (int j = 0; j < 32; j++)
    if (j >= CIV[1] || j < CIV[2]) {
      if ( getbit(aux,CIV[0],j) == 1 ) {
        set1(X0,CIV[0],j);
        set1(X1,CIV[0],j);
      } else {
        set0(X0,CIV[0],j);
        set0(X1,CIV[0],j);
      }
    }
}

int delta(ECRYPT_ctx *Z0, ECRYPT_ctx *Z1, int *OD){
  return ((Z0->state[OD[0]] >> OD[1]) & 1) ^ ((Z1->state[OD[0]] >> OD[1]) & 1);
}

/* Estimate forward bias, conditioned optional */
float Ed( int *ID, int *OD, int Nk, int Niv, int *CIV) {
  ECRYPT_ctx X0, X1, Z0, Z1, aux;
  int ones0, ones1;
  float bias0, bias1;
  float bias[Nk];
  for (int keys = 0; keys < Nk; keys++){
    bitstring(k,32); // Fix key
    ECRYPT_keysetup(&X0,k,256,64);
    ECRYPT_keysetup(&X1,k,256,64);
    ones0 = 0; 
    ones1 = 0;
    for (int ivs = 0; ivs < Niv; ivs++){
      bitstring(v,16); // Random IVs for each key
      ECRYPT_ivsetup(&X0,v);
      ECRYPT_ivsetup(&X1,v);
      if ( CIV != NULL ) {
        if (ivs == 0) ECRYPT_ivsetup(&aux,v);
        else fix(&X0,&X1,&aux,CIV);
      }
      setID(&X0,&X1,ID);
      if ( CIV != NULL ) {
        set0(&X0,CIV[0],CIV[2]); 
        set0(&X1,CIV[0],CIV[2]);
        salsa(4,Z0.state,X0.state);
        salsa(4,Z1.state,X1.state);
        ones0 += delta(&Z0,&Z1,OD);
        set1(&X0,CIV[0],CIV[2]); 
        set1(&X1,CIV[0],CIV[2]);
      }
      salsa(4,Z0.state,X0.state);
      salsa(4,Z1.state,X1.state);
      ones1 += delta(&Z0,&Z1,OD);
    }
    bias[keys] = fabs(2.0*ones1/Niv - 1);
    if ( CIV != NULL ) {
      bias0 = fabs(2.0*ones0/Niv - 1);
      bias1 = fabs(2.0*ones1/Niv - 1);
      bias[keys] = (bias0>=bias1)*bias0+(bias1>bias0)*bias1;
    }
  } 
  qsort(bias, Nk, sizeof(float), compare );
  return bias[Nk/2-1];
}
