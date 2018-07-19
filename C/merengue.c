#include "ecrypt-sync.h"
#include "ecrypt.c"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

u8 k[32]; u8 v[16];
const int new = 1, old = 0;

typedef struct
{
  int n;
  int len;
  int word[256];
  int bit[256];
  float bias[256];
} PNB;

typedef struct
{
  int len;
  int word[4];
  int ini[4];
  int end[4];
} CIV;

#define MINUS(v,w) (U32V((v) - (w)))

/* Modified Salsa20/R function */
static void salsa(int R, u32 output[16], const u32 input[16], int feedforward) {
  u32 x[16];
  int i;
  for (i = 0; i < 16; ++i) x[i] = input[i];
  for (i = R; i > 0 ;i -= 2) {
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
    if ( i != 1 ){
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
  if ( feedforward ) {
    for (i = 0;i < 16;++i) output[i] = PLUS(x[i],input[i]);
  } else for (i = 0;i < 16;++i) output[i] = x[i];
}

/* Inverse Salsa function */
static void aslas(int R, u32 output[16], const u32 input[16]) {
  u32 x[16];
  int i;
  for (i = 0; i < 16; ++i) x[i] = input[i];
  for (i = 0; i < R; i += 2) {
    if ( i > 0 || (R % 2) == 0 ) {
      x[15] = XOR(x[15],ROTATE(PLUS(x[14],x[13]),18));
      x[14] = XOR(x[14],ROTATE(PLUS(x[13],x[12]),13));
      x[13] = XOR(x[13],ROTATE(PLUS(x[12],x[15]), 9));
      x[12] = XOR(x[12],ROTATE(PLUS(x[15],x[14]), 7));
      x[10] = XOR(x[10],ROTATE(PLUS(x[ 9],x[ 8]),18));
      x[ 9] = XOR(x[ 9],ROTATE(PLUS(x[ 8],x[11]),13));
      x[ 8] = XOR(x[ 8],ROTATE(PLUS(x[11],x[10]), 9));
      x[11] = XOR(x[11],ROTATE(PLUS(x[10],x[ 9]), 7));
      x[ 5] = XOR(x[ 5],ROTATE(PLUS(x[ 4],x[ 7]),18));
      x[ 4] = XOR(x[ 4],ROTATE(PLUS(x[ 7],x[ 6]),13));
      x[ 7] = XOR(x[ 7],ROTATE(PLUS(x[ 6],x[ 5]), 9));
      x[ 6] = XOR(x[ 6],ROTATE(PLUS(x[ 5],x[ 4]), 7));
      x[ 0] = XOR(x[ 0],ROTATE(PLUS(x[ 3],x[ 2]),18));
      x[ 3] = XOR(x[ 3],ROTATE(PLUS(x[ 2],x[ 1]),13));
      x[ 2] = XOR(x[ 2],ROTATE(PLUS(x[ 1],x[ 0]), 9));
      x[ 1] = XOR(x[ 1],ROTATE(PLUS(x[ 0],x[ 3]), 7));
    }
    x[15] = XOR(x[15],ROTATE(PLUS(x[11],x[ 7]),18));
    x[11] = XOR(x[11],ROTATE(PLUS(x[ 7],x[ 3]),13));
    x[ 7] = XOR(x[ 7],ROTATE(PLUS(x[ 3],x[15]), 9));
    x[ 3] = XOR(x[ 3],ROTATE(PLUS(x[15],x[11]), 7));
    x[10] = XOR(x[10],ROTATE(PLUS(x[ 6],x[ 2]),18));
    x[ 6] = XOR(x[ 6],ROTATE(PLUS(x[ 2],x[14]),13));
    x[ 2] = XOR(x[ 2],ROTATE(PLUS(x[14],x[10]), 9));
    x[14] = XOR(x[14],ROTATE(PLUS(x[10],x[ 6]), 7));
    x[ 5] = XOR(x[ 5],ROTATE(PLUS(x[ 1],x[13]),18));
    x[ 1] = XOR(x[ 1],ROTATE(PLUS(x[13],x[ 9]),13));
    x[13] = XOR(x[13],ROTATE(PLUS(x[ 9],x[ 5]), 9));
    x[ 9] = XOR(x[ 9],ROTATE(PLUS(x[ 5],x[ 1]), 7));
    x[ 0] = XOR(x[ 0],ROTATE(PLUS(x[12],x[ 8]),18));
    x[12] = XOR(x[12],ROTATE(PLUS(x[ 8],x[ 4]),13));
    x[ 8] = XOR(x[ 8],ROTATE(PLUS(x[ 4],x[ 0]), 9));
    x[ 4] = XOR(x[ 4],ROTATE(PLUS(x[ 0],x[12]), 7));
  }
  for (i = 0;i < 16;++i) output[i] = x[i];
}

/* Store random bitstring of given byte length */
static void bitstring(u8 *array, int bytes) { 
  for (int i = 0; i < bytes; ++i) array[i] = rand() % 256;
}

static void keysetup(int flag, u8 k[32], ECRYPT_ctx *X0, ECRYPT_ctx *X1){
  if (flag == new) bitstring(k,32);
  ECRYPT_keysetup(X0,k,256,64);
  ECRYPT_keysetup(X1,k,256,64);
}

static void ivsetup(int flag, u8 iv[16], ECRYPT_ctx *X0, ECRYPT_ctx *X1){
  if (flag == new) bitstring(iv,16);
  ECRYPT_ivsetup(X0,iv);
  if (X1!=NULL) ECRYPT_ivsetup(X1,iv);
}

/* Answers whether the word position is a key */
static int iskey(int word) {
  return (word > 0 && word < 5) || (word > 10 && word < 15);
}

/* Answers whether the word position is a iv */
static int isiv(int word) {
  return (word > 5 && word < 10);
}

/* Compare two values */
static int compare(const void* a, const void* b) {
     int int_a = * ( (int*) a );
     int int_b = * ( (int*) b );
     if ( int_a == int_b ) return 0;
     else if ( int_a < int_b ) return -1;
     else return 1;
}

/* Set X[word][bit] = 0 */
static void set0(ECRYPT_ctx *X, ECRYPT_ctx *Y, int word, int bit) {
  X->state[word] = X->state[word] & ~ ( 1 << bit );
  if (Y!=NULL) Y->state[word] = Y->state[word] & ~ ( 1 << bit );
}

/* Set X[word][bit] = 1 */
static void set1(ECRYPT_ctx *X, ECRYPT_ctx *Y, int word, int bit){
  X->state[word] = (X->state[word] & ~(1<<bit)) | (1<<bit);
  if (Y!=NULL) Y->state[word] = (Y->state[word] & ~(1<<bit)) | (1<<bit);
}

/* Return bit value at X[word][bit] */
static int getbit(ECRYPT_ctx *X, int word, int bit){
  return ( X->state[word] >> bit ) & 1;
}

/* Flip bit value */
static void flip(ECRYPT_ctx *X, int word, int bit) {
  if ( getbit(X,word,bit) ) set0(X,NULL,word,bit);
  else set1(X,NULL,word,bit);
}

/* Count number of bit differences between two words */
static int numdifword(ECRYPT_ctx *Z0, ECRYPT_ctx *Z1, int word){
  int notequal = 0;
  for (int j = 0; j < 32; ++j)
    notequal += ((Z0->state[word] >> j) & 1) ^ ((Z1->state[word] >> j) & 1);
  return notequal;
}

/* Count number of bit differences between two states */
static int numdif(ECRYPT_ctx *Z0, ECRYPT_ctx *Z1){
  int notequal = 0;
  for (int i = 0; i < 16; ++i)
    notequal += numdifword(Z0,Z1,i);
  return notequal;
}

/* Print matrix state in hexadecimal format */
static void prettyprint(ECRYPT_ctx *X) {
  for (int i = 0; i < 16; ++i) {
    if ( i % 4 == 0 ) printf("\n");
    printf("%02x ",X->state[i]);
  } printf("\n");
}

/* Set input difference at initial states */
static void setID(ECRYPT_ctx *X0, ECRYPT_ctx *X1, int *ID) {
  set0(X0,NULL,ID[0],ID[1]); 
  set1(X1,NULL,ID[0],ID[1]); 
}

/* Copy bit to other states */
static void copybit(ECRYPT_ctx *X0, ECRYPT_ctx *X1, ECRYPT_ctx *aux, int word, int bit){
  if ( getbit(aux,word,bit) == 1 ) set1(X0,X1,word,bit);
  else set0(X0,X1,word,bit);
}

/* Copy fragment of bitstring to other states */
static void copyrange(ECRYPT_ctx *X0, ECRYPT_ctx *X1, ECRYPT_ctx *aux, int word, int ini, int end){
  for (int j = 0; j < 32; ++j)
    if (j >= ini || j <= end) {
      copybit(X0,X1,aux,word,j);
    }
}

/* Fix conditioned bits of IV */
static void fixCIV(ECRYPT_ctx *X0, ECRYPT_ctx *X1, ECRYPT_ctx *aux, CIV *civ, int ivs){
  if ( civ != NULL ) { // Fix conditioned bits
    if (ivs == 0) ivsetup(old,v,aux,NULL);
    else { // Need to fix whole first word if many conditions
      if ( civ->len > 1 ) copyrange(X0,X1,aux,civ->word[0],0,31);
      for (int c = 0; c < civ->len; ++c){
        copyrange(X0,X1,aux,civ->word[c],civ->ini[c],civ->end[c]);
      } 
    }
  }
}

static void fixFPSBs(ECRYPT_ctx *X0, ECRYPT_ctx *X1, ECRYPT_ctx *aux, PNB *fpnb, int ivs){
  if (fpnb!=NULL) {
    if (ivs == 0) ivsetup(old,v,aux,NULL);
    else { // Fix FPNBs of the IV
      for ( int i = 0 ; i < fpnb->n; ++i ) {
        if (isiv(fpnb->word[i])) copybit(X0,X1,aux,fpnb->word[i],fpnb->bit[i]);
      }
    }
  }
}

/* Set random bit */
static void setrandom (ECRYPT_ctx *X0, ECRYPT_ctx *X1, int word, int bit) {
  int r = rand() % 2 ;
  if ( r ) set1(X0,X1,word,bit);
  else set0(X0,X1,word,bit);
}

/* Fix PNKB to zero to build function g */
static void fixPNKBs(int type, PNB *pnb, ECRYPT_ctx *X0g, ECRYPT_ctx *X1g) {
  for ( int i = 0 ; i < pnb->len ; ++i ) {
    if ( iskey(pnb->word[i]) ) {
      if ( type == 1 ) set1(X0g,X1g,pnb->word[i],pnb->bit[i]);
      else if ( type == 0 ) set0(X0g,X1g,pnb->word[i],pnb->bit[i]);
      else setrandom(X0g,X1g,pnb->word[i],pnb->bit[i]);
    }
  }
}

/* Compute XOR difference of two output states at one position */
static int delta(ECRYPT_ctx *Z0, ECRYPT_ctx *Z1, int *OD){
  return ((Z0->state[OD[0]] >> OD[1]) & 1) ^ ((Z1->state[OD[0]] >> OD[1]) & 1);
}

/* Compute bias */
static float biasformula(int ones, int Niv, int absolute) {
  if (absolute) return fabs(2.0*ones/Niv -1);
  else return (2.0*ones)/Niv-1;
}

/* Returns highest value */
static float highest(float bias0, float bias1){
  return (bias0>=bias1)*bias0+(bias1>bias0)*bias1;
}

/* Return best bias between two threads */
static float bestbias(int ones0, int ones1, int Niv, int absolute) {
  float bias0,bias1;
  bias0 = biasformula(ones0,Niv,1);
  bias1 = biasformula(ones1,Niv,1);
  if (absolute) return highest(bias0,bias1);
  else return (bias0>=bias1)*(2.0*ones0/Niv-1)+(bias1>bias0)*(2.0*ones1/Niv-1);
}

static int salsadelta(int r, ECRYPT_ctx *X0, ECRYPT_ctx *X1, int *OD){
  ECRYPT_ctx Z0,Z1;
  salsa(r,Z0.state,X0->state,0);
  salsa(r,Z1.state,X1->state,0);
  return delta(&Z0,&Z1,OD);
}

/* Median value in array */
static float median(float *array, int num){
  qsort(array, num, sizeof(float), compare );
  return array[num/2-1];
}

static int isFPNB(PNB *fpnb, int word, int bit){
  for (int i = 0; i < fpnb->n; ++i) {
    if (word == fpnb->word[i] && bit == fpnb->bit[i]) return 1; 
  }
  return 0;
}

void getFPSBs(PNB *fpnb, PNB *fpsb) {
  if (fpnb!=NULL){
    fpsb->n = 0;
    for (int i = 6; i < 10; ++i) {
      for (int j = 0; j < 32; ++j) {
        if (!isFPNB(fpnb,i,j)) {
          fpsb->word[fpsb->n] = i;
          fpsb->bit[fpsb->n] = j;
          fpsb->n += 1;
        }
      }
    }
  }
}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////


/* Estimate forward bias, conditioned optional */
float getEd( int r, int *ID, int *OD, CIV *civ, PNB *fpnb, int Nk, int Niv) {
  ECRYPT_ctx X0, X1, aux;
  PNB fpsb;
  int ones00, ones01, ones10, ones11;
  float Ed[Nk];
  getFPSBs(fpnb,&fpsb);
  for (int keys = 0; keys < Nk; ++keys){
    keysetup(new,k,&X0,&X1);
    ones00 = 0; ones01 = 0, ones10 = 0; ones11 = 0;
    for (int ivs = 0; ivs < Niv; ++ivs){
      ivsetup(new,v,&X0,&X1);
      fixFPSBs(&X0,&X1,&aux,&fpsb,ivs);
      fixCIV(&X0,&X1,&aux,civ,ivs);
      setID(&X0,&X1,ID);
      if ( civ != NULL ) {
        if ( civ->len > 1 ) {
          set0(&X0,&X1,civ->word[1],civ->end[1]);
          set0(&X0,&X1,civ->word[0],civ->end[0]); 
          ones00 += salsadelta(r,&X0,&X1,OD);
          set1(&X0,&X1,civ->word[0],civ->end[0]); 
          ones01 += salsadelta(r,&X0,&X1,OD);
          set1(&X0,&X1,civ->word[1],civ->end[1]);
        }
        set1(&X0,&X1,civ->word[0],civ->end[0]);
        ones11 += salsadelta(r,&X0,&X1,OD);
        set0(&X0,&X1,civ->word[0],civ->end[0]); 
      }
      ones10 += salsadelta(r,&X0,&X1,OD);
    }
    if ( civ == NULL ) Ed[keys] = biasformula(ones10,Niv,1);
    else if ( civ->len == 1) Ed[keys] = bestbias(ones11,ones10,Niv,1);
    else Ed[keys] = highest(bestbias(ones00,ones01,Niv,1),bestbias(ones10,ones11,Niv,1)); 
  }
  return median(Ed,Nk);
}

/* Auxiliar thread function for neutrality */
static int backflip(ECRYPT_ctx *X0, ECRYPT_ctx *X1, int i, int j, int *OD) {
  int D, G;
  ECRYPT_ctx Y0, Y1, Z0, Z1;
  D = salsadelta(4,X0,X1,OD);
  salsa(8,Z0.state,X0->state,1);
  salsa(8,Z1.state,X1->state,1);
  flip(X0,i,j); 
  flip(X1,i,j);
  for ( int w = 0; w < 16; ++w ) {
    Z0.state[w] = MINUS(Z0.state[w],X0->state[w]);
    Z1.state[w] = MINUS(Z1.state[w],X1->state[w]);
  }
  aslas(4,Y0.state,Z0.state);
  aslas(4,Y1.state,Z1.state);
  G = delta(&Y0,&Y1,OD);
  return 1 - (D ^ G);
}

/* Compute neutrality measure of a bit */
static float neutrality(int i, int j, int *ID, int *OD, CIV *civ, PNB *fpnb, int Nk, int Niv) {
  ECRYPT_ctx X0, X1, Y0, Y1, Z0, Z1, aux;
  int equal00, equal01, equal10, equal11;
  PNB fpsb;
  float bias[Nk], bias0, bias1;
  getFPSBs(fpnb,&fpsb);
  for (int keys = 0; keys < Nk; ++keys){
    keysetup(new,k,&X0,&X1);
    equal00 = 0; equal01 = 0; equal10 = 0; equal11 = 0;
    for (int ivs = 0; ivs < Niv; ++ivs){
      ivsetup(new,v,&X0,&X1);
      fixFPSBs(&X0,&X1,&aux,&fpsb,ivs);
      fixCIV(&X0,&X1,&aux,civ,ivs);
      setID(&X0,&X1,ID);
      if ( civ != NULL ) {
        if ( civ->len == 2){
          set0(&X0,&X1,civ->word[1],civ->end[1]); 
          set0(&X0,&X1,civ->word[0],civ->end[0]); 
          equal00 += backflip(&X0,&X1,i,j,OD);
          set1(&X0,&X1,civ->word[0],civ->end[0]); 
          equal01 += backflip(&X0,&X1,i,j,OD);
          set1(&X0,&X1,civ->word[1],civ->end[1]); 
        }
        set1(&X0,&X1,civ->word[0],civ->end[0]); 
        equal11 += backflip(&X0,&X1,i,j,OD);
        set0(&X0,&X1,civ->word[0],civ->end[0]); 
      }  
      equal10 += backflip(&X0,&X1,i,j,OD);
    }
    if ( civ == NULL ) bias[keys] = biasformula(equal10,Niv,0);
    else if ( civ->len == 1) bias[keys] = bestbias(equal11,equal10,Niv,0);
    else {
      bias0 = bestbias(equal00,equal01,Niv,0);
      bias1 = bestbias(equal10,equal11,Niv,0);
      if ( fabs(bias0) > fabs(bias1) ) bias[keys] = bias0;
      else bias[keys] = bias1; 
    }
  }
  return median(bias,Nk);
}

/* Update list of PNBs depending on flag is forwards or backwards */
static void updatePNB(PNB *pnb, float bias, float gama, int i, int j, int flag){
  if (bias > gama) {
    printf("is PNB\n");
    if (iskey(i) && flag ) pnb->n += 1;
    if (isiv(i) && !flag ) pnb->n += 1;
    pnb->bias[pnb->len] = bias;
    pnb->word[pnb->len] = i;
    pnb->bit[pnb->len] = j;
    pnb->len += 1;
  }
}

/* Compute backwards probabilistic neutral bits */
void getBPNBs( PNB *bpnb, float gama, int *ID, int *OD, CIV *civ, PNB *fpnb, int Nk, int Niv) {
  float bias;
  for ( int i = 0 ; i < 16 ; ++i ){
    if ( iskey(i)){
    for ( int j = 0 ; j < 32 ; ++j ) {
      bias = neutrality(i,j,ID,OD,civ,fpnb,Nk,Niv);
      printf("(%d,%d) bias %f\n",i,j,bias);
      updatePNB(bpnb,bias,gama,i,j,1);
    }
    } 
  }
}

/* Auxiliar function for forwards() */
static int forflip(int r, ECRYPT_ctx *X0, ECRYPT_ctx *X1, int i, int j, int *OD){
  int f0, f1;
  set0(X0,X1,i,j);
  f0 = salsadelta(r,X0,X1,OD);
  set1(X0,X1,i,j);
  f1 = salsadelta(r,X0,X1,OD);
  return 1 - (f0 ^ f1);
}

/* Compute bias of flipping one initial bit */
static float forwards(int r, int i, int j, int *ID, int *OD, CIV *civ, int Nk, int Niv){
  ECRYPT_ctx X0,X1,aux;
  int equal00, equal01, equal10,equal11;
  float bias[Nk], bias0, bias1;
  for (int keys = 0; keys < Nk; keys++){
    keysetup(new,k,&X0,&X1);
    equal00 = 0; equal01 = 0; equal10 = 0; equal11 = 0;
    for (int ivs = 0; ivs < Niv; ivs++){
      ivsetup(new,v,&X0,&X1);
      fixCIV(&X0,&X1,&aux,civ,ivs);
      setID(&X0,&X1,ID);
      if ( civ != NULL ){
        if ( civ->len > 1 ){
          set0(&X0,&X1,civ->word[1],civ->end[1]);
          set0(&X0,&X1,civ->word[0],civ->end[0]);
          equal00 += forflip(r,&X0,&X1,i,j,OD);
          set1(&X0,&X1,civ->word[0],civ->end[0]);
          equal01 += forflip(r,&X0,&X1,i,j,OD);
          set1(&X0,&X1,civ->word[1],civ->end[1]);
        }
        set1(&X0,&X1,civ->word[0],civ->end[0]);
        equal11 += forflip(r,&X0,&X1,i,j,OD);
        set0(&X0,&X1,civ->word[0],civ->end[0]); 
      }
      equal10 += forflip(r,&X0,&X1,i,j,OD);
    }
    if ( civ == NULL ) bias[keys] = biasformula(equal10,Niv,0);
    else if ( civ->len == 1) bias[keys] = bestbias(equal11,equal10,Niv,0);
    else {
      bias0 = bestbias(equal00,equal01,Niv,0);
      bias1 = bestbias(equal10,equal11,Niv,0);
      if ( fabs(bias0) > fabs(bias1) ) bias[keys] = bias0;
      else bias[keys] = bias1; 
    }
  }
  return median(bias,Nk);
}

/* Compute forwards probabilistic neutral bits */
void getFPNBs(PNB *fpnb, float gama, int r, int *ID, int *OD, CIV *civ, int Nk, int Niv) {
  float bias;
  for ( int i = 0 ; i < 16 ; ++i ){
    if ( isiv(i) ){
      for ( int j = 0 ; j < 32 ; ++j ) {
        bias = forwards(r,i,j,ID,OD,civ,Nk,Niv);
        printf("(%d,%d) bias %f\n",i,j,bias);
        updatePNB(fpnb,bias,gama,i,j,0);
      }
    }
  }
}

/* Compute backwards equality f = g */
static int backwards(ECRYPT_ctx *X0f, ECRYPT_ctx *X1f, ECRYPT_ctx *X0g, ECRYPT_ctx *X1g, int *OD) {
  ECRYPT_ctx Y0, Y1, Z0, Z1;
  salsa(8,Z0.state,X0f->state,1);
  salsa(8,Z1.state,X1f->state,1);
  for ( int w = 0; w < 16; ++w ) {
    Z0.state[w] = MINUS(Z0.state[w],X0g->state[w]);
    Z1.state[w] = MINUS(Z1.state[w],X1g->state[w]);
  }
  aslas(4,Y0.state,Z0.state);
  aslas(4,Y1.state,Z1.state);
  return delta(&Y0,&Y1,OD);      
}

static int fequalg(ECRYPT_ctx *X0f, ECRYPT_ctx *X1f, ECRYPT_ctx *X0g, ECRYPT_ctx *X1g, int *OD){
  int f = salsadelta(4,X0f,X1f,OD);
  int g = backwards(X0f,X1f,X0g,X1g,OD);
  return 1 - ( f ^ g );
}

/* Estimate backwards bias */
float getEa(PNB *bpnb, int *ID, int *OD, CIV *civ, PNB *fpnb, int Nk, int Niv) {
  ECRYPT_ctx X0f, X1f, X0g, X1g, aux;
  int equal00, equal01, equal10, equal11;
  PNB fpsb;
  float Ea[Nk];
  getFPSBs(fpnb,&fpsb);
  for ( int keys = 0 ; keys < Nk ; ++keys ) {
    keysetup(new,k,&X0f,&X1f);
    keysetup(old,k,&X0g,&X1g);
    fixPNKBs(-1, bpnb,&X0g,&X1g);
    equal00 = 0; equal01 = 0; equal10 = 0; equal11 = 0;
    printf("key %d\n",keys);
    for ( int ivs = 0 ; ivs < Niv ; ++ivs ){
      ivsetup(new,v,&X0f,&X1f);
      ivsetup(old,v,&X0g,&X1g);
      fixFPSBs(&X0f,&X1f,&aux,&fpsb,ivs);
      fixFPSBs(&X0g,&X1g,&aux,&fpsb,ivs);
      fixCIV(&X0f,&X1f,&aux,civ,ivs);
      fixCIV(&X0g,&X1g,&aux,civ,ivs);
      setID(&X0f,&X1f,ID);
      setID(&X0g,&X1g,ID);
      if ( civ != NULL ) {
        if (civ->len == 2 ){
          set0(&X0f,&X1f,civ->word[1],civ->end[1]); 
          set0(&X0g,&X1g,civ->word[1],civ->end[1]);
          set0(&X0f,&X1f,civ->word[0],civ->end[0]); 
          set0(&X0g,&X1g,civ->word[0],civ->end[0]);
          equal00 += fequalg(&X0f,&X1f,&X0g,&X1g,OD);
          set1(&X0f,&X1f,civ->word[0],civ->end[0]); 
          set1(&X0g,&X1g,civ->word[0],civ->end[0]);
          equal01 += fequalg(&X0f,&X1f,&X0g,&X1g,OD);
          set1(&X0f,&X1f,civ->word[1],civ->end[1]); 
          set1(&X0g,&X1g,civ->word[1],civ->end[1]);
        }
        set1(&X0f,&X1f,civ->word[0],civ->end[0]); 
        set1(&X0g,&X1g,civ->word[0],civ->end[0]);
        equal11 += fequalg(&X0f,&X1f,&X0g,&X1g,OD);
        set0(&X0f,&X1f,civ->word[0],civ->end[0]); 
        set0(&X0g,&X1g,civ->word[0],civ->end[0]);
      }  
      equal10 += fequalg(&X0f,&X1f,&X0g,&X1g,OD);
    }
    if ( civ == NULL ) Ea[keys] = biasformula(equal10,Niv,1);
    else if ( civ->len == 1) Ea[keys] = bestbias(equal11,equal10,Niv,1);
    else Ea[keys] = highest(bestbias(equal00,equal01,Niv,1),bestbias(equal10,equal11,Niv,1)); 
  }
  return median(Ea,Nk);
}

float getE(PNB *bpnb, int *ID, int *OD, CIV *civ, PNB *fpnb, int Nk, int Niv) {
  ECRYPT_ctx X0f, X1f, X0g, X1g, aux;
  int ones00, ones01, ones10, ones11;
  PNB fpsb;
  float E[Nk];
  getFPSBs(fpnb,&fpsb);
  for ( int keys = 0 ; keys < Nk ; ++keys ) {
    printf("key %d\n",keys);
    keysetup(new,k,&X0f,&X1f);
    keysetup(old,k,&X0g,&X1g);
    fixPNKBs(-1,bpnb,&X0g,&X1g);
    ones00 = 0; ones01 = 0; ones10 = 0; ones11 = 0;
    for ( int ivs = 0 ; ivs < Niv ; ++ivs ){
      ivsetup(new,v,&X0f,&X1f);
      ivsetup(old,v,&X0g,&X1g);
      fixFPSBs(&X0f,&X1f,&aux,&fpsb,ivs);
      fixFPSBs(&X0g,&X1g,&aux,&fpsb,ivs);
      fixCIV(&X0f,&X1f,&aux,civ,ivs);
      fixCIV(&X0g,&X1g,&aux,civ,ivs);
      setID(&X0g,&X1g,ID);
      if ( civ != NULL ) {
        if (civ->len == 2 ){
          set0(&X0f,&X1f,civ->word[1],civ->end[1]); 
          set0(&X0g,&X1g,civ->word[1],civ->end[1]);
          set0(&X0f,&X1f,civ->word[0],civ->end[0]); 
          set0(&X0g,&X1g,civ->word[0],civ->end[0]);
          ones00 += backwards(&X0f,&X1f,&X0g,&X1g,OD);
          set1(&X0f,&X1f,civ->word[0],civ->end[0]); 
          set1(&X0g,&X1g,civ->word[0],civ->end[0]);
          ones01 += backwards(&X0f,&X1f,&X0g,&X1g,OD);
          set1(&X0f,&X1f,civ->word[1],civ->end[1]); 
          set1(&X0g,&X1g,civ->word[1],civ->end[1]);
        }
        set1(&X0f,&X1f,civ->word[0],civ->end[0]); 
        set1(&X0g,&X1g,civ->word[0],civ->end[0]);
        ones11 += backwards(&X0f,&X1f,&X0g,&X1g,OD);
        set0(&X0f,&X1f,civ->word[0],civ->end[0]); 
        set0(&X0g,&X1g,civ->word[0],civ->end[0]);
      }  
      ones10 += backwards(&X0f,&X1f,&X0g,&X1g,OD);
    }
    if ( civ == NULL ) E[keys] = biasformula(ones10,Niv,1);
    else if ( civ->len == 1) E[keys] = bestbias(ones11,ones10,Niv,1);
    else E[keys] = highest(bestbias(ones00,ones01,Niv,1),bestbias(ones10,ones11,Niv,1)); 
  }  
  return median(E,Nk);
}