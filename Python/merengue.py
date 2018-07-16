# Inria SECRET 2018
# anais.querol-cruz@inria.fr
# 
# Cryptanalysis tool for Bernstein's Salsa20


##########################
# BASIC CIPHER FUNCTIONS # 
##########################

def transpose (X):
    for row in range(4):
        for col in range(4):
            if row < col:               # For all cells in upper right triangle
                X = swap(row, col, X)   # Swap to lower left triangle accordingly
    return X


def swap (row, col, X):
    X[row*4+col], X[col*4+row] = X[col*4+row], X[row*4+col]
    return X                       # (i,j) <-> (j,i)


def QuarterRound(X, mode):
    if mode == "inv":               # Inverse
        X = transpose(X)            # Transpose before
    X = column(X,4,mode)            # First column starts from X4
    X = column(X,9,mode)            # Second column starts from X9
    X = column_s(X,14,mode)         # Third column starts from X14
    X = column_s(X,3,mode)          # Fourth column starts from X3
    if mode != "inv":               # Encrypt, Carrify, Linearize, Propagate, Differential,
        X = transpose(X)            # Transpose after operations on columns
    return X


def column(X, k, mode):
    r = [7, 9, 13, 18]              # First reverse diagonal element
    if mode == "inv" :              # Inverse
        for n in range(4)[::-1]:    # Apply ARX to this cell, using two above
            X[(k+n*4)%16] = ARX(X[(k+(n-2)*4)%16],X[(k+(n-1)*4)%16],r[n],X[(k+n*4)%16])
    else:                           # Encryption, propagate, carrify, linearize, Differential
        for n in range(4):          # Apply ARX to this cell, using two above
            X[(k+n*4)%16] = ARX(X[(k+(n-2)*4)%16],X[(k+(n-1)*4)%16],r[n],X[(k+n*4)%16],mode)
    return X


def ARX(a1, a2, r, x, mode):
    if mode == "pda":
        return XOR(x, ROT(r, ADP(a1,a2,mode),mode),mode)
    else: return XOR(x, ROT(r, ADD(a1,a2,mode),mode),mode)    # b (+)= [ (a + d) <<< r ]


def ADD(A, B, mode):
    if mode == "enc" :  # Encrypt
        added = ( bit2int(A) + bit2int(B) ) % 2**32   # Sum of two words
        return int2bit(added)
    elif mode == "pro": # Propagate
        return [ A[n] + B[n] for n in range(32) ]     # Dependencies in same bit gain importance
    elif mode == "car": # Carrify
        C = [ 0.0 for _ in range(32) ]   
        return [ FA(A[n],B[n],C,n) for n in range(32)[::-1] ][::-1]
    elif mode == "lin": # Linear
        added = [0]*32
        added[31] = str(int(A[31]) ^ int(B[31]))
        for n in range(31)[::-1]:
            added[n] = str( int(A[n]) ^ int(B[n]) ^ int(A[n+1]) )
        return added


def ROT(r, o, mode):
    if mode != "pro":           # Encrypt, Carrify, Linear
        return o[r:] + o[:r]    # Rotate bits to the left
    elif mode == "pro":         # Propagate
        return o[-r:] + o[:-r]  # Dependencies are rotated to the left


def XOR(A, B, mode):
    if mode == "enc" or mode == "lin":  # Encrypt, Linear
        xored = bit2int(A) ^ bit2int(B)
        return int2bit(xored)           # Exclusive OR
    elif mode == "pro":                 # Propagate
        return ADD(A, B, mode)          # XORing dependencies acts like adding
    elif mode == "car" or mode == "pda" # Carrify, Differential
        if isinstance(A,list) and isinstance(B,list):
            return [ A[i]+B[i]-2*A[i]*B[i] for i in range(min(len(A),len(B),32))]
        else:
            return A+B-2*A*B


####################
# INITIALIZE STATE #
####################

from random import SystemRandom

C = [0]*4
C[0] = bin(int("61707865", 16))[2:].zfill(32) # X[0]
C[1] = bin(int("3320646e", 16))[2:].zfill(32) # X[5]
C[2] = bin(int("79622d32", 16))[2:].zfill(32) # X[10]
C[3] = bin(int("6b206574", 16))[2:].zfill(32) # X[15]


def bitstring(bits):
    """
    Random bitstring of certain size in bits.
    Input:  bits    bitstring length in bits, must be a multiple of 32
    Output  bs      random bitstring in list of words of 32 bits
    """
    if mode == "enc" or mode == "lin" :  # Encrypt, Linearized
        cryptogen = SystemRandom()
        return [ [ str(cryptogen.randrange(2)) for i in range(32)] for j in range(bits/32)]
    elif mode == "car": # Carrify
        bs = [ [random.choice([0.0,1.0]) for _ in range(32)] for n in range(bits/32)]
        return bs


def diagonal(X, mode = None):
    if mode == "enc":   # Encrypt
        X[0]  = hex2bit("61707865") # c0
        X[5]  = hex2bit("3320646e") # c1
        X[10] = hex2bit("79622d32") # c2
        X[15] = hex2bit("6b206574") # c3
    elif mode == "car"  # Carrify 
        X[0]  = hex2prb("61707865") # c0
        X[5]  = hex2prb("3320646e") # c1
        X[10] = hex2prb("79622d32") # c2
        X[15] = hex2prb("6b206574") # c3
    return


def inivalue(X, iv, ID=None, bit=None):
    X[6:10] = copy(iv)
    if ID is not None and bit is not None:
        X[ID[0]][31-ID[1]] = bit
    return


def keywords(X,key):
    X[1:5] = copy(key[0:4])
    X[11:15] = copy(key[4:])
    return


################
# Salsa Cipher #
################

def differential(R, key, ID, ff = None):
    iv = bitstring(128,"enc")
    X0, iv0, Z0 = Salsa(R,key,iv,ff,ID,'0')    
    X1, iv1, Z1 = Salsa(R,key,iv,ff,ID,'1')
    return X0, X1, iv, Z0, Z1


def Salsa (R, key, iv = None, ff = None, ID = None, bit = None):
    """
    Compute two Salsa keystreams in R rounds with ID in IV
    Input:  R   number of rounds
            key list of 8 words for 256 bit secret key
            ID  tuple with position of input difference ([6..9],[0..31])
    Output: X0  512bit keystream with word of ID value 0
            X1  512bit keystream with word of ID value 1
    """
    mode = "enc"
    X = [ [] for n in range(16)]    # Initialize "array"
    diagonal(X,mode)
    if iv is None:
        iv = bitstring(128,mode)
    inivalue(X, iv, ID, bit)
    keywords(X, key)
    Z = copy(X)
    for n in range(R):
        QuarterRound(Z,mode)        # Apply QuarterRound R times
    if R % 2 == 1:                  # Odd rounds over rows instead
        transpose(Z)                # Need to transpose again before output
    if ff is not None:              # Apply feedforward after all rounds
        for word in range(16):
            Z[word] = ADD(X[word], Z[word], mode)
    return X, iv, Z


#################
# Inverse Salsa #
#################

def Aslas(R,X):
    mode = "inv"
    Z = copy(X)
    if R % 2:
        Z = transpose(Z)
    for n in range(R):
        Z = QuarterRound(Z,mode)
    return Z


def feedbackwards(Z, iv, key=None):
    Z_X = [['0']*32 for n in range(16)]
    Z_X[0] = MINUS(Z[0],hex2bit("61707865"))
    Z_X[5] = MINUS(Z[5],hex2bit("3320646e"))
    Z_X[10] = MINUS(Z[10],hex2bit("79622d32"))
    Z_X[15] = MINUS(Z[15],hex2bit("6b206574"))
    Z_X[6] = MINUS(Z[6],iv[0])
    Z_X[7] = MINUS(Z[7],iv[1])
    Z_X[8] = MINUS(Z[8],iv[2])
    Z_X[9] = MINUS(Z[9],iv[3])
    if key is not None:
        Z_X[1] = MINUS(Z[1],key[0])
        Z_X[2] = MINUS(Z[2],key[1])
        Z_X[3] = MINUS(Z[3],key[2])
        Z_X[4] = MINUS(Z[4],key[3])
        Z_X[11] = MINUS(Z[11],key[4])
        Z_X[12] = MINUS(Z[12],key[5])
        Z_X[13] = MINUS(Z[13],key[6])
        Z_X[14] = MINUS(Z[14],key[7])
    return Z_X


def MINUS(z, x):
    return int2bit((bit2int(z) - bit2int(x))%2**32)

def subtract(Z,Y):
    X = [['0']*32 for n in range(16)]
    for n in range(16):
        X[n] = MINUS(Z[n],Y[n])
    return X

def addition(Z,Y):
    A = [['0']*32 for n in range(16)]
    for n in range(16):
        A[n] = int2bit( (bit2int(Z[n]) + bit2int(Y[n]))%2**32 )
    return A


################################
# Linear propagation functions #
################################

def propagate (R, i, j):
    """
    Compute dependencies of differential (i,j) over R rounds of linear Salsa.
    Input:  R   number of rounds
            i   index of word [0..15]
            j   index of bit [0..31]
    Output: X   array of size 16x32 
    """
    mode = "pro"
    X = [ [0]*32 for n in range(16)]    # Initialize "array"
    X[i][j] = 1                         # Only bit difference in Xi,j
    for n in range(R):
        X = QuarterRound(X,mode)        # Apply QuarterRound R times
    if R % 2 == 1:                      # Odd rounds over rows instead
        X = transpose(X)                # Need to transpose again before output
    return X


#################
# Salsa carries #
#################

def carrify(R, ff=None, key=None, iv=None, rand=None):
    """
    Probabilistic outputs
    """
    mode = "car"
    X = [ [] for n in range(R+1)]
    if key is not None and iv is not None:
        X[0] = [ [0.0]*32 for n in range(16)]
        keywords(X[0], key)
        inivalue(X[0], iv)
    elif rand is not None:
        X[0] = [ [0.0]*32 for n in range(16)]
        keywords(X[0], bitstring(256,mode))
        inivalue(X[0], bitstring(128,mode))
    else:
        X[0] = [ [0.5]*32 for n in range(16)]    # Same probability 0/1
    diagonal(X[0],0)
    for n in range(R):
        X[n+1] = copy(X[n])
        QuarterRound(X[n+1],mode)
    for n in range(1,R+1)[::2]:
        transpose(X[n])
    if ff is not None:
        for word in range(16):
            X[R][word] = ADD(X[0][word], X[R][word],mode)
    return X


def FA(A,B,C,n):
    """1 bit Full Adder"""
    AxB, AyB = HA(A, B)             # HA(A,B) = A x B , A & B
    AxBxC, AxByC = HA(AxB, C[n])    # HA(XOR(A,B),C) = (A x B) x C, (A x B) & C
    C[(n-1)%32] = OR(AxByC, AyB)    # Update output carry
    return AxBxC                    # Return sum value

def HA(A, B): return XOR(A,B), AND(A,B) # 1 bit Half Adder

def AND(A,B): return A * B

def OR(A,B): return A + B - A * B


####################
# Linearized Salsa #
####################

def linearized(R,key,iv,ff=True):
    mode = "lin"
    X = [ [] for n in range(16)]    # Initialize "array"
    diagonal(X,mode)
    inivalue(X, iv)
    keywords(X, key)
    Z = copy(X)
    for n in range(R):
        QuarterRound(Z,mode)        # Apply QuarterRound R times
    if R % 2 == 1:                  # Odd rounds over rows instead
        transpose(Z)                # Need to transpose again before output
    if ff is True:                  # Apply feedforward after all rounds
        for word in range(16):
            Z[word] = ADD(X[word], Z[word], mode)
    return Z


##########################
# Unroll Salsa functions #
##########################

def rowof(p): return p / 4

def colof(p): return p % 4

def trans4x4(p):
    i = rowof(p)
    j = colof(p)
    return j * 4 + i

def mix(X1,X2):
    X = [ [ X1[i][j] + bit for j, bit in enumerate(word)] for i, word in enumerate(X2) ]
    return X

def constant(X):
    global C
    non = nnZ(X)
    if sum(non) == 1:
        word = non.index(1)
        if word % 5 == 0:       # Diagonal word
            bit = X[word].index(1)
            cnt = C[word/5][31-bit]
            return (word, bit, int(cnt))
    return False

def unroll (R, p, q, mode):
    """
    Recursive inverse Salsa function.
    Input   R   number of rounds
            p   Salsa output word index
            q   Salsa output bit index
            mode  "ful" or "red"
    Output  X   Updated dependency array
            S   Updated dependency array
    """
    if not R:
        Z = [[0]*32 for _ in range(16)] # Initialize array
        Z[p][q] += 1                    # Arrived to initial cell
        T = "X[%s][%s]" % (p,q)         # Stringfy operation
        return Z, T
    if R % 2 == 0:                      # Rows
        p = trans4x4(p)                 # Transpose word index
    r = [7, 9, 13, 18]                  # Rotation constants
    b = [4, 9, 14, 3]                   # Below diagonal words
    col = colof(p)                      # Index of column used
    row = rowof(p)                      # Index of row used
    p1 = (p-4) % 16                     # Index of first word argument
    p2 = (p-2*4) % 16                   # Index of second word argument    
    k = colof(row-rowof(b[col]))        # Order of present cell
    qr = (q - r[k]) % 32                # Index of rotated bit
    if R % 2 == 0:                      # Rows
        p = trans(p)                    # Transpose word index
        p1 = trans4x4(p1)               # Transpose first word index
        p2 = trans4x4(p2)               # Transpose second word index
    if k == 0:                          # Below diagonal
        Z, T = xor( unroll(R-1,p,q), add(unroll(R-1,p1,qr), unroll(R-1,p2,qr),mode), mode )
    elif k == 1:                        # Below-below diagonal
        Z, T = xor( unroll(R-1,p,q), add(unroll(R,p1,qr), unroll(R-1,p2,qr),mode), mode )
    else:                               # Above diagonal and diagonal
        Z, T = xor( unroll(R-1,p,q), add(unroll(R,p1,qr), unroll(R,p2,qr)mode), mode )  
    return Z, T                         # Return array and string


def xor( (X1,S1), (X2,S2), mode ):
    X = mix(X1,X2)
    if mode == "red":
        c1 = constant(X1)
        c2 = constant(X2)
        if c1 and c1[2]:    return X, neg(S2)
        elif c2 and c2[2]:  return X, neg(S1)
        elif c1:            return X, S2
        elif c2:            return X, S1
    elif mode == "ful":     return X, "xor(%s , %s)" % (S1, S2)

def add( (X1,S1), (X2,S2), mode ):
    X = mix(X1,X2)
     if mode == "red":
        c1 = constant(X1)           # (word1, bit1, value1)
        c2 = constant(X2)           # (word2, bit2, value2)
        if c1 and c2:       return X, add2ct(c1[2], c2[2])  
        elif c1 and not c2: return X, add1ct(c1[2], S2)
        elif not c1 and c2: return X, add1ct(c2[2], S1)
    elif mode == "ful"      return X, "add(%s , %s)" % (S1, S2)

def add1ct(x, S):
    if x:                       # 1 + expression
        return "(%s, carry(0.5))" % neg(S)
    return S                    # 0 + expression

def add2ct(x1,x2):
    carry = x1 and x2           # Carry only if 1 + 1
    if carry:
        return "(0, carry(1))"
    return "%s" % (int(x1) ^ x2)     # ADD without carry behaves as XOR
    
def neg(S):
    if S[0:4] == "not(":
        return S[4:-1]          # not(not( * )) = *
    return "not(%s)" % S        # not(S)

def xorof(Z):
    return [ [ b for b, n in enumerate(z) if n % 2 ] for z in Z]

def xordiag(X):
    x = 0
    for n in range(4):
        idx = X[n*5]
        for dep in idx:
            x = x ^ int(C[n][31-dep])
    return x 


##############
# Unroll 3/4 #
##############

def unrollin (R, p, q):
    """
    Recursive inverse Salsa function. Counts dependencies of linerized 3/4 probability
    Input   R   number of rounds
            p   Salsa output word index
            q   Salsa output bit index
    Output  X   Updated dependency array
            S   Updated dependency array
    """
    if not R:
        Z = [[0]*32 for _ in range(16)] # Initialize array
        Z[p][q] += 1                    # Arrived to initial cell
        return Z
    if R % 2 == 0:                  # Rows
        p = trans(p)                # Transpose word index
    r = [7, 9, 13, 18]              # Rotation constants
    b = [4, 9, 14, 3]               # Below diagonal words
    col = colof(p)                  # Index of column used
    row = rowof(p)                  # Index of row used
    p1 = (p-4) % 16                 # Index of first word argument
    p2 = (p-2*4) % 16               # Index of second word argument    
    k = colof(row-rowof(b[col]))    # Order of present cell
    qr = (q - r[k]) % 32            # Index of rotated bit
    if R % 2 == 0:                  # Rows
        p = trans(p)                # Transpose word index
        p1 = trans(p1)              # Transpose first word index
        p2 = trans(p2)              # Transpose second word index
    lsb = qr == 0;                  # Whether in LSB or not. No need to linearize, LSB already behaves like XOR
    if k == 0:                      # Below diagonal
        Z = xorlin( unrollin(R-1,p,q), addlin(unrollin(R-1,p1,qr), unrollin(R-1,p2,qr), unrollin(R-1,p1,qr-1),lsb ) )
    elif k == 1:                    # Below-below diagonal
        Z = xorlin( unrollin(R-1,p,q), addlin(unrollin(R,p1,qr), unrollin(R-1,p2,qr), unrollin(R-1,p1,qr-1)),lsb )
    else:                           # Above diagonal and diagonal
        Z = xorlin( unrollin(R-1,p,q), addlin(unrollin(R,p1,qr), unrollin(R,p2,qr), unrollin(R-1,p1,qr-1)),lsb )   
    return Z                        # Return array and string

def xorlin( X1, X2 ):
    X = mix(X1,X2)
    return X

def addlin( A, B, A_1, lsb):
    X = mix(A,B)
    if not lsb:
        X = mix(X, A_1)
    return X


##################################
# Probabilistic Differential ARX #
##################################

def diffMachine (R, Delta):
    """ Probabilistic Differentials of ARX """ 
    mode = "pda"
    for r in range(R):
        Delta[r+1] = copy(Delta[r])
        QuarterRound(Delta[r+1],mode)
    for r in range(1,R+1)[::2]:
        transpose(Delta[r])
    return Delta


def ADP(alfa,beta):
    carry = 0.0
    gama = [0.0]*len(alfa)
    A = [ [0.0, 0.0], [0.0, 0.0] ]
    P = [ [0.0, 0.0], [0.0, 0.0] ]
    for i in range(len(alfa))[::-1]:
        A[0][0] = [carry, carry/2.0]
        A[0][1] = [1-carry, 0.5]   
        A[1][0] = [1-carry, 0.5]    
        A[1][1] = [carry, (1+carry)/2.0]
        P[0][0] = (1-alfa[i][0])*(1-beta[i])
        P[0][1] = (1-alfa[i][0])*beta[i]
        P[1][0] = alfa[i][0]*(1-beta[i])
        P[1][1] = alfa[i][0]*beta[i] 
        carry = 0.0
        for n in range(2):
            for m in range(2):
                gama[i] += P[n][m]*A[n][m][0] # value
                carry += P[n][m]*A[n][m][1]   # carry
    return gama


def DP(alfa,beta):
    lon = len(alfa)     # Exhaustive
    diff = [ [0]*lon for _ in range(2**(2*lon)) ]
    for n1 in range(2**lon):
        for n2 in range(2**lon):
            dx1 = int2bit(n1)
            dx2 = int2bit(n2)
            for i in range(lon):
                dx1[31-i] = str((int(dx1[31-i]) + alfa[lon-1-i])%2)
                dx2[31-i] = str((int(dx2[31-i]) + beta[lon-1-i])%2)
            sx = (n1 + n2)%2**lon
            sdx = (bit2int(dx1) + bit2int(dx2)) % 2**lon
            sx = int2bit(sx)
            sdx = int2bit(sdx)
            for j in range(lon):
                diff[n1*lon+n2][lon-1-j] += sx[31-j] != sdx[31-j]
    gama = [0]*lon
    for i in range(len(diff)):
        for j in range(lon):
            gama[j] += diff[i][j]/float(len(diff))
    return gama


######################################
###### TENTATIVE CRYPTANALYSIS #######
######################################

def compare(S,Z,i=None,j=None):
    eq = 0
    if i is None and j is None:
        for i in range(16):
            for j in range(32):
                eq += S[i][j] == Z[i][j]
    else:
        eq = S[i][31-j] == Z[i][31-j]
    return eq


def multibit(R,i,out):
    """
    Compute minimum effect of ID over all 16 words and all 32 bits after many rounds
    Input   R   number of rounds for Salsa
            i   word position of ID -> [0..15]
    Output  T   list of lists with tuples of quasineutral bits
    """
    T = [ neutral(R,i,j) for j in range(32) ]
    prettytuple(T,i,out)
    return T


def neutral(R, i, j):
    """
    Compute minimum effect of ID over all 16 words after many rounds
    Input   R   number of rounds for Salsa
            i   word position of ID -> [0..15]
            j   bit  position of ID -> [0..31]
    Output  T   list of 16 tuples containing (minimum effect, bit position)
    """
    Z = propagate(R,i,j)
    T = [ ( min(z), indexof(z, min(z)) )  for z in Z ]
    return T


def indexof(X, n):
    """
    List of bit positions where X[j] == n
    Input   X   list of values
            n   value inside X
    Output  idx index(es) of value n in X
    """
    idx = [j for j,x in enumerate(X) if x == n]
    return idx


def PNBs (R,p,q):
    """
    Compute truly neutral bits in whole block for a given differential
    Input   R   number of rounds for Salsa
            p   index of word for OD
            q   index of bit  for OD
    Output  pnb list of PNBs found in whole block
    """
    pnb = [ [] for _ in range(16)]
    for i in range(16):
        for j in range(32):
            Z = propagate(R,i,j)
            if not Z[p][q]:
                pnb[i].append(j)
    return pnb


def PNBkey (pnb):
    """
    Select neutral bits corresponding to key positions
    Input   pnb list of PNBs found for whole block
    Output  n   number of PNBs found (the larger the better)
            key list of PNBs in key positions
    """
    key = pnb[1:5]+pnb[11:15] # Key is X1,X2,X3,X4,X11,X12,X13,X14
    n = sum([ len(k) for k in key if k ])
    return n, key


def PNBest (R):
    """
    Compute best bits for OD that give the largest number of neutral bits
    Input   R   number of rounds
    Output  n   largest number of neutral bits
            p   word index with largest number of neutral bits
            q   bit  index with largest number of neutral bits
    """
    n = 0       # Initialize. Invariant n >= 0
    p = -1      # Initialize. Invariant p : [0..15]
    q = -1      # Initialize. Invariant q = [0..31]
    for P in range(16):
        for Q in range(32):
            N, key = PNBkey(PNBs(R,P,Q))
            if n < N:   # If larger number of PNBs
                n = N   # Update value
                p = P   # Update word coordinate
                q = Q   # Update bit  coordinate
    return n, p, q


def struct(R):
    """
    Creates large structure with all bits dependencies
    Input   R   number of rounds
    Output  A   large structure should be used A[i][j][p][q]
    """
    A = [ [ [] for _ in range(32) ] for _ in range(16)]
    for i in range(16):
        for j in range(32):
            for p in range(16):
                for q in range(32):
                    A[i][j] = propagate(R,i,j)
    return A


def nonegl():
    """
    Compute number of rounds until all bits depend on all bits
    Input:  -
    Output: r   maximum number of rounds for non-negligible differential
    """
    maxR = [0]*16 
    for idx in range(16):           # For all words
        R = 0                       # Initialize zero rounds
        d = 0                       # Initialize zero dependencies
        while d < 32:               # Until all bits 
            R = R + 1               # Apply one more round
            d = bestdif(R,idx,0)    # Number of dependencies of best differential
        maxR[idx] = R - 1           # Maximum of rounds with non-negl is one before
    return min(maxR)                # Rounds is minimum of all words (some may be slower)


def nnZ (Z):
    return [ len(nnz(z)) for z in Z]     # Number of nonzero dependencies in each word


def nnz (l):
    return [ p for p,e in enumerate(l) if e ]   # Word positions with nonzero elements in list


def check (R):
    """
    Check smallest number of dependencies with R rounds
    Input:  R   number of rounds for Salsa
    Output: D   minimum number of dependencies after R rounds
            i   index of input word with minimum d
            j   index of input bit with minimum d
            p   index of output word that causes minimum dependencies
    """
    D = 33                          # Initialize. Invariant D <= 32
    d = 33                          # Initialize. Invariant d <= 32
    i = -1                          # Initialize. Invariant i = [0..15]
    j = -1                          # Initialize. Invariant j = [0..31]
    for idx in range(16):           # For all words for ID
        for jdx in range(32):       # For all bits in word for ID
            Z = Salsa(R,idx,jdx)    # Execute Salsa with 5 rounds
            n = nnZ(Z)              # Number of dependencies for this ID
            d = min(n)              # Get minimum number
            if D > d:               # If it is smaller
                D = d               # Update minimum value
                i = idx             # Update location of word
                j = jdx             # Update location of bit
                p = n.index(d)      # Index of output word
    return D, i, j, p


def bestdif (R, i, out):
    """ 
    Compute best single bit differential for a word. ( [r/_\p]q | [0/_\i]j )
    Input:  R       number of rounds
            i       index of word [0..15]
            out     flag for print or no print
    Output: mind    minimum number of bits depending on ID 
            bitj    bit in Xi of the best input differential
            word    index of word of the best output differential
            bitq    bit in Zp of the best output differential
            maxe    maximum number of occurences of bit Xij in the OD word
    """
    mind = 33       # Initialize. Invariant mind <= 32 (#bits)
    maxe = -1       # Initialize. Invariant maxe >= 0
    bitj = -1       # Initialize. Invariant bitj : [0..31]
    bitq = -1       # Initialize. Invariant bitq : [0..31]
    word = -1       # Initialize. Invariant word : [0..15]
    for j in range(32):             # For all possible bits for ID
        Z = propagate(R,i,j)        # Compute Salsa dependencies with ID in Xij
        m, p = mindep(Z)            # Compute slower words with smallest number of dependencies
        if m <= mind:               # If this value is smaller so far, need to update 
            for z in p:             # For each slower word (may be a few with same number)
                M,B = maxeff(Z[z])  # Find bits more affected by ID
                if M >= maxe:       # If this value is larger so far, need to update
                    maxe = M        # Update maximum number of effects
                    word = z        # Update index of new best word for OD
                    bitq = B[0]     # If many, choose any of them
            mind = m        # Update minimum number of dependencies
            bitj = j        # Update new best bit for Xi for ID
    if out:
        print "Best differential: ( [ "+str(R)+"/_\ "+str(word)+"]"+str(bitq)+" | [ "+str(0)+"/_\ "+str(i)+"]"+str(bitj)+" ) -> #depend =",mind, ", #effect =", maxe
        print "     Bit in word X"+str(i)+" bit "+str(bitj)+" in round 0"
        print "     correlated to word Z"+str(word)+" bit "+str(bitq)+" after round "+str(R)
    else: return mind #, maxe, bitj, word, bitq


def fullword (R, i):
    for j in range(32):             # For all bits in the word
        Z = propagate(R,i,j)        # Compute Salsa dependencies of the differential
        m, p = mindep(Z)            # Compute words with fewer dependencies -> more time to propagate
        print 'X'+str(i)+','+str(j), '-> mindep =', m, 'in: '
        for z in p:                 # For each slow word in output
            M, B = maxeff(Z[z])     # Compute most dependent bit
            print '    -','z'+str(z), '-> maxeff = ', M, 'bit', [b for b in B]
    return


def mindep (Z):
    m = 33                      # Initialize. Invariant m <= 32
    p = [-1]                    # Initialize list of word indexes : [0..15]
    for k, z in enumerate(Z):   # List of words depending on ID
        nz = nnz(z)             # Count how many bits in each word depend on ID
        if m == len(nz):        # If word with same value
            p.append(k)         # Include index of this slow word
        elif m > len(nz):       # The smaller the better
            m = len(nz)         # Update new best value
            p = [k]             # Update list only with this new best word
    return m, p


def maxeff (z):
    M = -1                      # Initialize. Invariant M >= 0
    B = [-1]                    # Initialize list of bit positions : [0..31]
    for b, e in enumerate(z):   # List of times each bit is affected by ID
        if M == e:              # If bit with same value
            B.append(b)         # Include also this bit index
        elif M < e:             # The larger the better
            M = e               # Update new best value
            B = [b]             # Update list only with this new best bit
    return M, B


####################
# Format functions #
####################

def copy(old):
    return [ list(word) for word in old]

def printall(stream):
    st = [ ''.join( [str(b) for b in stream[n]] ) for n in range(len(stream)) ]
    return st

def joinbits(bs):
    bst = ''.join(bs)
    return bst

def hexall(X):
    return [ hex(bit2int(x))[2:].zfill(8) for x in X]

def bit2int(word):
    return int(joinbits(word),2)

def prb2hex(prb):
    return hex(int(joinbits(prb2bit(prb)), 2))

def int2bit(word):
    word = list(bin(int(word))[2:].zfill(32))
    return word

def hex2bit(hexa):
    word = list(bin(int(hexa, 16))[2:].zfill(32))
    return word

def hex2prb(hexa):
    word = list(bin(int(hexa, 16))[2:].zfill(32))
    word = [float(bit) for bit in word]
    return word

def prb2bit(prbs):
    if isinstance(prbs[0], list):
        return [ [ '0' if bit==0.0 else '1' for bit in word] for word in prbs ]
    else:
        return [ '0' if bit==0.0 else '1' for bit in prbs]

######################
# Printing functions #
######################

def beautify (Z):
    for p,z in enumerate(Z):    # For each word in Salsa output array
        nz = nnz(z)
        # word : |bits depend on ID| , indexes of such bits 
        print "Z"+str(p), ": dep = ", len(nz), " bit =", prettybit(nz)

def prettybit (x):
    s = '[ '
    for b in x:
        s += str(b)+' '     # Compact print list of numbers without ',' in between
    s += ']'
    return s

def prettytuple(TT,i,out):
    Xi = "X"+str(i)
    print "Minimum effects of ID "+Xi+" :"
    for j, T in enumerate(TT):  # 32 times
        if out or j==31:
            print " -------"
            print "|",Xi+","+str(j), "|"
            print " -------"
            for p, t in enumerate(T):
                print "      Z"+str(p)+" : only "+ str(t[0]) + " effects in bits " + str(t[1])

def question(D):
    return [ [ '?' if bit!=0.0 and bit!=1.0 else bit for bit in word ] for word in D]
    
def abbrev(D):
    return [ [ "%.3f"%bit for bit in word ] for word in D]

def simply(S):
    global C
    T = S + ""
    for n in range(4):
        for b in range(32):
            T = T.replace("X[%s][%s]"%(n*5,b), C[n][31-b])
    return T
