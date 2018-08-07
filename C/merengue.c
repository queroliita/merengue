#include "ecrypt-sync.h"
#include "ecrypt.c"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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

u8 k[32]; u8 v[16];
const int new = 1, old = 0;

CIV *civ;
PNB *fpnb, *bpnb;
int ID[2], OD[256][2], R, r, nODs, nIDs;
float th;
unsigned int Nk;
unsigned long int Niv;

static void colround (u32 x[16]){
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
}

static void locround (u32 x[16]){
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

static void rowround (u32 x[16]) {
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

static void worround (u32 x[16]) {
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

// Modified Salsa20/R function
static void salsa(int rounds, u32 output[16], const u32 input[16], int feedforward) {
  u32 x[16];
  int i;
  for (i = 0; i < 16; ++i) x[i] = input[i];
  for (i = 0; i < rounds ;i += 2) {
    colround(x);
    if ( i != rounds - 1 ) rowround(x);
  }
  if ( feedforward ) {
    for (i = 0;i < 16;++i) output[i] = PLUS(x[i],input[i]);
  } else for (i = 0;i < 16;++i) output[i] = x[i];
}

// Inverse Salsa function 
static void aslas(u32 output[16], const u32 input[16]) {
  u32 x[16];
  int i, Rounds = R;
  for (i = 0; i < 16; ++i) x[i] = input[i];
  if ( R % 2 == 1 ) { locround(x); Rounds = R - 1; }
  for (i = Rounds; i > r; i -= 2) {
    worround(x);
    if ( i != r+1 ) locround(x);
  }
  for (i = 0;i < 16;++i) output[i] = x[i];
}

// Store random bitstring of given byte length
static void bitstring(u8 *array, int bytes) { 
  for (int i = 0; i < bytes; ++i) array[i] = rand() % 256;
}

static void keysetup(int flag, ECRYPT_ctx *X0, ECRYPT_ctx *X1){
  if (flag == new) bitstring(k,32);
  if (!X0->null) ECRYPT_keysetup(X0,k,256,64);
  if (!X1->null) ECRYPT_keysetup(X1,k,256,64);
}

static void ivsetup(int flag, ECRYPT_ctx *X0, ECRYPT_ctx *X1){
  if (flag == new) bitstring(v,16);
  if (!X0->null) ECRYPT_ivsetup(X0,v);
  if (!X1->null) ECRYPT_ivsetup(X1,v);
}

// Answers whether the word position is a key
static int iskey(int word) {
  return (word > 0 && word < 5) || (word > 10 && word < 15);
}

// Answers whether the word position is a iv
static int isiv(int word) {
  return (word > 5 && word < 10);
}

static int isFPNB(int word, int bit){
  for (int i = 0; i < fpnb->n; ++i) {
    if (word == fpnb->word[i] && bit == fpnb->bit[i]) return 1; 
  }
  return 0;
}

// Return bit value at X[word][bit] 
static int getbit(ECRYPT_ctx *X, int word, int bit){
  return ( X->state[word] >> bit ) & 1;
}

// Set X[word][bit] = 0 
static void set0(ECRYPT_ctx *X, int word, int bit) {
  if (!X->null) X->state[word] = X->state[word] & ~ ( 1 << bit );
}

// Set X[word][bit] = 1
static void set1(ECRYPT_ctx *X, int word, int bit){
  if (!X->null) X->state[word] = (X->state[word] & ~(1<<bit)) | (1<<bit);
}

// Flip bit value 
static void flip(ECRYPT_ctx *X, int word, int bit) {
  if ( getbit(X,word,bit) ) set0(X,word,bit);
  else set1(X,word,bit);
}

// Set input difference at initial states 
static void setID(ECRYPT_ctx *X0, ECRYPT_ctx *X1) {
  set0(X0,ID[0],ID[1]); 
  set1(X1,ID[0],ID[1]);
}

// Set random bits
static void setrandom (ECRYPT_ctx *X0, ECRYPT_ctx *X1, int word, int bit) {
  if ( rand() % 2 ) { set1(X0,word,bit); set1(X1,word,bit); }
  else              { set0(X0,word,bit); set0(X1,word,bit); }
}

// Count number of bit differences between two words
static int numdifword(ECRYPT_ctx *Z0, ECRYPT_ctx *Z1, int word){
  int notequal = 0;
  for (int j = 0; j < 32; ++j)
    notequal += ((Z0->state[word] >> j) & 1) ^ ((Z1->state[word] >> j) & 1);
  return notequal;
}

// Count number of bit differences between two states
static int numdif(ECRYPT_ctx *Z0, ECRYPT_ctx *Z1){
  int notequal = 0;
  for (int i = 0; i < 16; ++i)
    notequal += numdifword(Z0,Z1,i);
  return notequal;
}

// Copy bit to other states 
static void copybit(ECRYPT_ctx *X0, ECRYPT_ctx *X1, ECRYPT_ctx *aux, int word, int bit){
  if ( getbit(aux,word,bit) ) { set1(X0,word,bit); set1(X1,word,bit); }
  else                        { set0(X0,word,bit); set0(X1,word,bit); }
}

// Copy fragment of bitstring to other states
static void copyrange(ECRYPT_ctx *X0, ECRYPT_ctx *X1, ECRYPT_ctx *aux, int word, int ini, int end){
  for (int j = 0; j < 32; ++j)
    if (j >= ini || j <= end) {
      copybit(X0,X1,aux,word,j);
    }
}

// Fix conditioned bits of IV 
static void fixCIV(ECRYPT_ctx *X0, ECRYPT_ctx *X1, ECRYPT_ctx *aux, int ivs){
  if ( civ != NULL ) { // Fix conditioned bits
    if (ivs == 0) ECRYPT_ivsetup(aux,v);
    else { // Need to fix whole first word if many conditions
      if ( civ->len > 1 ) copyrange(X0,X1,aux,civ->word[0],0,31);
      for (int c = 0; c < civ->len; ++c){
        copyrange(X0,X1,aux,civ->word[c],civ->ini[c],civ->end[c]);
      } 
    }
  }
}

static void fixFPSBs(ECRYPT_ctx *X0, ECRYPT_ctx *X1, ECRYPT_ctx *aux, PNB *fpsb, int ivs){
  if (fpnb!=NULL) {
    if (ivs == 0) ECRYPT_ivsetup(aux,v);
    else { // Fix FPNBs of the IV
      for ( int i = 0 ; i < fpnb->n; ++i ) {
        if (isiv(fpsb->word[i])) copybit(X0,X1,aux,fpsb->word[i],fpsb->bit[i]);
      }
    }
  }
}

// Fix PNKB to zero to build function g 
static void fixPNKBs(int type, ECRYPT_ctx *X0g, ECRYPT_ctx *X1g) {
  for ( int i = 0 ; i < bpnb->len ; ++i ) {
    if ( iskey(bpnb->word[i]) ) {
      if      ( type == 1 ) { set1(X0g,bpnb->word[i],bpnb->bit[i]); set1(X1g,bpnb->word[i],bpnb->bit[i]); }
      else if ( type == 0 ) { set0(X0g,bpnb->word[i],bpnb->bit[i]); set1(X1g,bpnb->word[i],bpnb->bit[i]); }
      else setrandom(X0g,X1g,bpnb->word[i],bpnb->bit[i]);
    }
  }
}

// Compare two values 
static int compare(const void* a, const void* b) {
     int int_a = * ( (int*) a );
     int int_b = * ( (int*) b );
     if ( int_a == int_b ) return 0;
     else if ( int_a < int_b ) return -1;
     else return 1;
}

static int same(char *str, char *txt){
  return strncmp(str,txt,8) == 0;
}

// Print matrix state in hexadecimal format 
static void prettyprint(ECRYPT_ctx *X) {
  for (int i = 0; i < 16; ++i) {
    if ( i % 4 == 0 ) printf("\n");
    printf("%02x ",X->state[i]);
  } printf("\n");
}

// Compute bias
static float biasformula(unsigned long int ones, int absolute) {
  if (absolute) return fabs(2.0*ones/Niv -1);
  else return (2.0*ones)/Niv-1;
}

// Returns highest value 
static float highest(float bias0, float bias1){
  return (bias0>=bias1)*bias0+(bias1>bias0)*bias1;
}

// Return best bias between two threads
static float bestbias(unsigned long int ones0, unsigned long int ones1, int absolute) {
  float bias0,bias1;
  bias0 = biasformula(ones0,1);
  bias1 = biasformula(ones1,1);
  if (absolute) return highest(bias0,bias1);
  else return (bias0>=bias1)*(2.0*ones0/Niv-1)+(bias1>bias0)*(2.0*ones1/Niv-1);
}

static float maximum(float *array, int len) {
  float max = array[0];
  for (int i = 0; i < len; i++)
    max = array[i]*(array[i]>max);
  return max;
}

// Average value in array
static float avg(float array[Nk]) {
  float sum = 0.0;
  for (unsigned int i = 0; i < Nk; i++) { sum += array[i]; }
  return sum / Nk ; 
}
// Median value in array
static float median(float *array, unsigned int num){
  qsort(array, num, sizeof(float), compare );
  if ( num % 2 ) return array[(num-1)/2];       // Odd
  else return (array[num/2-1]+array[num/2])/2;  // Even
}

// Compute XOR difference of two output states at one position
static int delta(ECRYPT_ctx *Z0, ECRYPT_ctx *Z1){
  int xor = 0;
  for (int o = 0; o < nODs ; o++) 
    xor = xor ^ ((Z0->state[OD[o][0]] >> OD[o][1]) & 1) ^ ((Z1->state[OD[o][0]] >> OD[o][1]) & 1);
  return xor;
}

static int salsadelta(int rounds, ECRYPT_ctx *X0, ECRYPT_ctx *X1){
  ECRYPT_ctx Z0,Z1;
  salsa(rounds,Z0.state,X0->state,0);
  salsa(rounds,Z1.state,X1->state,0);
  return delta(&Z0,&Z1);
}

void getFPSBs(PNB *fpsb) {
  if (fpnb!=NULL){
    fpsb->n = 0;
    for (int i = 6; i < 10; ++i) {
      for (int j = 0; j < 32; ++j) {
        if (!isFPNB(i,j)) {
          fpsb->word[fpsb->n] = i;
          fpsb->bit[fpsb->n] = j;
          fpsb->n += 1;
        }
      }
    }
  }
}

// Auxiliar thread function for neutrality
static int backflip(ECRYPT_ctx *X0, ECRYPT_ctx *X1, int i, int j) {
  int D, G;
  ECRYPT_ctx Y0, Y1, Z0, Z1;
  D = salsadelta(r,X0,X1);
  salsa(R,Z0.state,X0->state,1);
  salsa(R,Z1.state,X1->state,1);
  flip(X0,i,j); 
  flip(X1,i,j);
  for ( int w = 0; w < 16; ++w ) {
    Z0.state[w] = MINUS(Z0.state[w],X0->state[w]);
    Z1.state[w] = MINUS(Z1.state[w],X1->state[w]);
  }
  aslas(Y0.state,Z0.state);
  aslas(Y1.state,Z1.state);
  G = delta(&Y0,&Y1);
  return 1 - (D ^ G);
}

// Auxiliar function for forwards()
static int forflip(ECRYPT_ctx *X0, ECRYPT_ctx *X1, int i, int j){
  int f0, f1;
  set0(X0,i,j); set0(X1,i,j);
  f0 = salsadelta(r,X0,X1);
  set1(X0,i,j); set1(X1,i,j);
  f1 = salsadelta(r,X0,X1);
  return 1 - (f0 ^ f1);
}

// Compute backwards equality f = g 
static int backwards(ECRYPT_ctx *X0f, ECRYPT_ctx *X1f, ECRYPT_ctx *X0g, ECRYPT_ctx *X1g) {
  ECRYPT_ctx Y0, Y1, Z0, Z1;
  salsa(R,Z0.state,X0f->state,1);
  salsa(R,Z1.state,X1f->state,1);
  for ( int w = 0; w < 16; ++w ) {
    Z0.state[w] = MINUS(Z0.state[w],X0g->state[w]);
    Z1.state[w] = MINUS(Z1.state[w],X1g->state[w]);
  }
  aslas(Y0.state,Z0.state);
  aslas(Y1.state,Z1.state);
  return delta(&Y0,&Y1);      
}

static int fequalg(ECRYPT_ctx *X0f, ECRYPT_ctx *X1f, ECRYPT_ctx *X0g, ECRYPT_ctx *X1g){
  int f = salsadelta(r,X0f,X1f);
  int g = backwards(X0f,X1f,X0g,X1g);
  return 1 - ( f ^ g );
}

///////////////////////////////////////////////////////////////////////////////

static float updatebias(unsigned long int count[4], int abso){
  if      ( civ == NULL )   return biasformula(count[2],abso);
  else if ( civ->len == 1 ) return bestbias(count[2],count[3],abso);
  else if ( abso ) 
    return highest(bestbias(count[0],count[1],abso),bestbias(count[2],count[3],abso)); 
  else {
    float bias0 = bestbias(count[0],count[1],0);
    float bias1 = bestbias(count[2],count[3],0);
    if ( fabs(bias0) > fabs(bias1) ) return bias0;
    else return bias1; 
  }
}

static int measurement(char *flag, ECRYPT_ctx *X0f, ECRYPT_ctx *X1f, ECRYPT_ctx *X0g, ECRYPT_ctx *X1g,
  int i, int j, int rounds) {
  if      ( same(flag,"Ef") )   return salsadelta(rounds,X0f,X1f);
  else if ( same(flag,"Eg") )   return fequalg(X0f,X1f,X0g,X1g);
  else if ( same(flag,"E") )    return backwards(X0f,X1f,X0g,X1g);
  else if ( same(flag,"back") ) return backflip(X0f,X1f,i,j);
  else if ( same(flag,"for") )  return forflip(X0f,X1f,i,j);
  return 0;
}

static void set(int f, int g, ECRYPT_ctx Xf[2], ECRYPT_ctx Xg[2], int word, int bit) {
  if      ( f == 0 ) { set0(&Xf[0],word,bit); set0(&Xf[1],word,bit); }
  else if ( f == 1 ) { set1(&Xf[0],word,bit); set1(&Xf[1],word,bit); }
  if      ( g == 0 ) { set0(&Xg[0],word,bit); set0(&Xg[1],word,bit); }
  else if ( g == 1 ) { set1(&Xg[0],word,bit); set1(&Xg[1],word,bit); }
}

// Compute neutrality measure of a bit (backwards or forwards)
static float neutrality(char *flag, int i, int j) {
  ECRYPT_ctx Xf[2], Xg[2], aux;
  PNB fpsb;
  unsigned long int count[4];
  int absol = flag[0]=='E';
  float bias[Nk];
  Xf[0].null = 0; Xf[1].null = 0;
  if ( !same(flag,"Eg") && !same(flag,"E") ) { Xg[0].null = 1; Xg[1].null = 1; }
  getFPSBs(&fpsb);
  for (unsigned int key = 0; key < Nk; ++key) {
    keysetup(new,&Xf[0],&Xf[1]);
    keysetup(new,&Xg[0],&Xg[1]);
    for (int n = 0; n < 4; n++ ) { count[n] = 0; }
    for (unsigned long int ivs = 0; ivs < Niv; ++ivs) {
      ivsetup(new,&Xf[0],&Xf[1]);
      ivsetup(old,&Xg[0],&Xg[1]);
      fixFPSBs(&Xf[0],&Xf[1],&aux,&fpsb,ivs);
      fixFPSBs(&Xg[0],&Xg[1],&aux,&fpsb,ivs);
      fixCIV(&Xf[0],&Xf[1],&aux,ivs);
      fixCIV(&Xg[0],&Xg[1],&aux,ivs);
      setID(&Xf[0],&Xf[1]);
      setID(&Xg[0],&Xg[1]);
      if ( civ != NULL ) {
        if ( civ->len == 2){
          set(0,0,Xf,Xg,civ->word[1],civ->end[1]); // Second condition
          set(0,0,Xf,Xg,civ->word[0],civ->end[0]); // First condition
          count[0] += measurement(flag,&Xf[0],&Xf[1],&Xg[0],&Xg[1],i,j,r);
          set(1,1,Xf,Xg,civ->word[0],civ->end[0]);
          count[1] += measurement(flag,&Xf[0],&Xf[1],&Xg[0],&Xg[1],i,j,r);
          set(1,1,Xf,Xg,civ->word[1],civ->end[1]);
        }
        set(1,1,Xf,Xg,civ->word[0],civ->end[0]);
        count[3] += measurement(flag,&Xf[0],&Xf[1],&Xg[0],&Xg[1],i,j,r);
        set(0,0,Xf,Xg,civ->word[0],civ->end[0]);
      }  
      count[2] += measurement(flag,&Xf[0],&Xf[1],&Xg[0],&Xg[1],i,j,r);
    }
    bias[key] = updatebias(count,absol);
    printf("%u key with bias %f\n",key,bias[key]);
  }
  printf("avg(bias)=%f\n",avg(bias));
  return median(bias,Nk);
}

// Update list of PNBs depending on flag is forwards or backwards
static void updatePNB(char *flag, PNB *pnb, float bias, int i, int j){
  if (bias > th) {
    printf("is PNB\n");
    if (iskey(i) && same(flag,"back") ) pnb->n += 1;
    if (isiv(i) && same(flag,"for") ) pnb->n += 1;
    pnb->bias[pnb->len] = bias;
    pnb->word[pnb->len] = i;
    pnb->bit[pnb->len] = j;
    pnb->len += 1;
  }
}

// Compute forwards or backwards probabilistic neutral bits
void getPNBs(char *flag, PNB *pnb) {
  float bias;
  for ( int i = 0 ; i < 16 ; ++i ) {
    if ( (iskey(i) && same(flag,"back") ) || (isiv(i) && same(flag,"for")) ) {
      for ( int j = 0 ; j < 32 ; ++j ) {
        bias = neutrality(flag,i,j);
        printf("(%d,%d) bias %f\n",i,j,bias);
        updatePNB(flag,pnb,bias,i,j);
      }
    } 
  }
}

// Get bias Ef, Eg or E
float getbias(char *flag) {
  return neutrality(flag,-1,-1); 
}