#ifndef _PROTO_H_
#define _PROTO_H_

#define DIM 3

#ifndef NMK_DOUBLE_PREC
#define NMK_SINGLE_PREC // Using single precision as default, same in most other codes.
#endif 

#ifdef NMK_SINGLE_PREC // NMK_SINGLE_PREC
	#define SQRT sqrtf
	#define INVSQRT invsqrtf
	#define EXP expf
	#define ERFC (Real) 1.0 - erff
    typedef float Real;
#else // NMK_SINGLE_PREC
	#define SQRT sqrt
	#define INVSQRT invsqrt
	#define EXP exp
	#define ERFC (Real) 1.0 - erf
    typedef double Real;
#endif // NMK_SINGLE_PREC


typedef Real			Vect3[DIM];
typedef unsigned int	Short;
typedef unsigned int	Int;
typedef unsigned long	Long;
typedef unsigned long	Peano;
typedef unsigned int	Morton;

typedef struct {
    Long    id;      // constant particle ID
    Peano   key;     // peano-hilbert key
    Int     tag;     // 
    Int     group;   // comes from cuboid
	Morton	mkey;    // morton key in cuboid.
    Vect3   pos;     // position of particle
    Vect3   vel;     // velocity of particle
    Vect3   acc;     // acceration of particle
	Vect3   acc_pm;
    Real    mass;
    int     nstep;
} Body;

typedef struct
{
    Long    id;      // constant particle ID
    Int     tag;     // 
    Vect3   pos;     // position of particle
    Vect3   vel;     // velocity of particle
    Real	mass;
} TransBody;  // For transition between the processes.

typedef struct
{
    Vect3	pos;     // position of particle
    Vect3   acc;     // acceration of particle
    Real	mass;
} CalcBody;  // Body of N-Body, You knew it. */

typedef struct {
    Int     nPart;
    Int     firstPart;
    Int     firstChild;
    Int     nChildren; // childnum
    int     level;
    Morton  mortonKey;
    Real    width;
    Real    mass;
    Vect3   massCenter;
	Real	rmax;
} Node;

#endif /* _PROTO_H_ */

