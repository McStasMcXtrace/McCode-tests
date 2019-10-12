/* Automatically generated file. Do not edit. 
 * Format:     ANSI C source code
 * Creator:    McStas <http://www.mcstas.org>
 * Instrument: HZB_NEAT.instr (HZB_NEAT)
 * Date:       Sat Oct 12 09:05:00 2019
 * File:       HZB_NEAT.c
 * CFLAGS=
 */

#define MCCODE_STRING "McStas 3.0-dev - Oct. 11, 2019"
#define FLAVOR        "mcstas"
#define FLAVOR_UPPER  "MCSTAS"

#define MC_USE_DEFAULT_MAIN
#define MC_TRACE_ENABLED
#define PI 3.14159265358979323846
#ifndef M_PI
#define M_PI PI
#endif
#ifndef M_2_PI
#define M_2_PI 0.63661977236758134308
#endif
#ifndef M_PI_2
#define M_PI_2 1.57079632679489661923  /* pi/2 */
#endif
#ifndef M_SQRT2
#define M_SQRT2 1.4142135623730951  /* sqrt(2) */
#endif

#ifdef USE_PGI
#include <openacc_curand.h>
#endif

struct _struct_particle {
  double x,y,z; /* position [m] */
  double vx,vy,vz; /* velocity [m/s] */
  double sx,sy,sz; /* spin [0-1] */
#ifdef USE_PGI
  curandState_t MCRANDstate; /* CUDA RNG state */
#endif
  double t, p;    /* time, event weight */
  long long _uid;  /* event ID */
  long _index;     /* component index where to send this event */
  long _absorbed;  /* flag set to TRUE when this event is to be removed/ignored */
  long _scattered; /* flag set to TRUE when this event has interacted with the last component instance */
  long _restore;   /* set to true if neutron event must be restored */
};
typedef struct _struct_particle _class_particle;

#define MC_EMBEDDED_RUNTIME
/* embedding file "mccode-r.h" */

/*******************************************************************************
*
* McCode, neutron/xray ray-tracing package
*         Copyright (C) 1997-2009, All rights reserved
*         Risoe National Laboratory, Roskilde, Denmark
*         Institut Laue Langevin, Grenoble, France
*
* Runtime: share/mccode-r.h
*
* %Identification
* Written by: KN
* Date:    Aug 29, 1997
* Release: McStas 3.0-dev
* Version: $Revision$
*
* Runtime system header for McStas/McXtrace.
*
* In order to use this library as an external library, the following variables
* and macros must be declared (see details in the code)
*
*   struct mcinputtable_struct mcinputtable[];
*   int numipar;
*   char instrument_name[], instrument_source[];
*   int traceenabled, defaultmain;
*   extern MCNUM  mccomp_storein[];
*   extern MCNUM  mcAbsorbProp[];
*   extern MCNUM  mcScattered;
*   #define MCCODE_STRING "the McStas/McXtrace version"
*
* Usage: Automatically embbeded in the c code.
*
* $Id$
*
*******************************************************************************/

#ifndef MCCODE_R_H
#define MCCODE_R_H "$Revision$"

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <limits.h>
#include <errno.h>
#include <time.h>
#include <float.h>
#include <inttypes.h>

/* If the runtime is embedded in the simulation program, some definitions can
   be made static. */

#ifdef MC_EMBEDDED_RUNTIME
#define mcstatic
#else
#define mcstatic
#endif

#ifdef __dest_os
#if (__dest_os == __mac_os)
#define MAC
#endif
#endif

#ifdef __FreeBSD__
#define NEED_STAT_H
#endif

#if defined(__APPLE__) && defined(__GNUC__)
#define NEED_STAT_H
#endif

#ifdef NEED_STAT_H
#include <sys/stat.h>
#endif

#ifndef MC_PATHSEP_C
#ifdef WIN32
#define MC_PATHSEP_C '\\'
#define MC_PATHSEP_S "\\"
#else  /* !WIN32 */
#define MC_PATHSEP_C '/'
#define MC_PATHSEP_S "/"
#endif /* !WIN32 */
#endif /* MC_PATHSEP_C */



/* the version string is replaced when building distribution with mkdist */
#ifndef MCCODE_STRING
#define MCCODE_STRING "McStas 3.0-dev - Oct. 11, 2019"
#endif

#ifndef MCCODE_DATE
#define MCCODE_DATE "Oct. 11, 2019"
#endif

#ifndef MCCODE_VERSION
#define MCCODE_VERSION "3.0-dev"
#endif

#ifndef MCCODE_NAME
#define MCCODE_NAME "McStas"
#endif

#ifndef MCCODE_PARTICLE
#define MCCODE_PARTICLE "neutron"
#endif

#ifndef MCCODE_PARTICLE_CODE
#define MCCODE_PARTICLE_CODE 2112
#endif

#ifndef MCCODE_LIBENV
#define MCCODE_LIBENV "MCSTAS"
#endif

#ifndef FLAVOR_UPPER
#define FLAVOR_UPPER MCCODE_NAME
#endif

#ifdef MC_PORTABLE
#ifndef NOSIGNALS
#define NOSIGNALS 1
#endif
#endif

#ifdef MAC
#ifndef NOSIGNALS
#define NOSIGNALS 1
#endif
#endif

#if (USE_MPI == 0)
#undef USE_MPI
#endif

#ifdef USE_MPI  /* default is to disable signals with MPI, as MPICH uses them to communicate */
#ifndef NOSIGNALS
#define NOSIGNALS 1
#endif
#endif

#if (NOSIGNALS == 0)
#undef NOSIGNALS
#endif

/* Note: the enum instr_formal_types definition MUST be kept
   synchronized with the one in mccode.h and with the
   instr_formal_type_names array in cogen.c. */
enum instr_formal_types
  {
    instr_type_int,
    instr_type_string, instr_type_char,
    instr_type_vector, instr_type_double
  };
struct mcinputtable_struct { /* defines instrument parameters */
  char *name; /* name of parameter */
  void *par;  /* pointer to instrument parameter (variable) */
  enum instr_formal_types type;
  char *val;  /* default value */
};

typedef double MCNUM;
typedef struct {MCNUM x, y, z;} Coords;
typedef MCNUM Rotation[3][3];

/* the following variables are defined in the McStas generated C code
   but should be defined externally in case of independent library usage */
#ifndef DANSE
extern struct mcinputtable_struct mcinputtable[];         /* list of instrument parameters */
extern int    numipar;                                    /* number of instrument parameters */
extern char   instrument_name[], instrument_source[]; /* instrument name and filename */
extern char  *instrument_exe;                           /* executable path = argv[0] or NULL */
extern char   instrument_code[];                        /* contains the initial 'instr' file */

#ifndef MC_ANCIENT_COMPATIBILITY
extern int traceenabled, defaultmain;
#endif
#endif


/* Useful macros ============================================================ */


/* SECTION: Dynamic Arrays */
typedef double* DArray1d;
DArray1d create_darr1d(int n);
void destroy_darr1d(DArray1d a);

typedef double** DArray2d;
DArray2d create_darr2d(int nx, int ny);
void destroy_darr2d(DArray2d a);

typedef double*** DArray3d;
DArray3d create_darr3d(int nx, int ny, int nz);
void destroy_darr3d(DArray3d a);


/* MPI stuff */
#ifdef USE_MPI
#include "mpi.h"

#ifdef OMPI_MPI_H  /* openmpi does not use signals: we may install our sighandler */
#undef NOSIGNALS
#endif

/*
 * MPI_MASTER(i):
 * execution of i only on master node
 */
#define MPI_MASTER(statement) { \
  if(mpi_node_rank == mpi_node_root)\
  { statement; } \
}

#ifndef MPI_REDUCE_BLOCKSIZE
#define MPI_REDUCE_BLOCKSIZE 1000
#endif

int mc_MPI_Sum(double* buf, long count);
int mc_MPI_Send(void *sbuf, long count, MPI_Datatype dtype, int dest);
int mc_MPI_Recv(void *rbuf, long count, MPI_Datatype dtype, int source);

/* MPI_Finalize exits gracefully and should be preferred to MPI_Abort */
#define exit(code) do {                                   \
    MPI_Finalize();                                       \
    exit(code);                                           \
  } while(0)

#else /* !USE_MPI */
#define MPI_MASTER(instr) instr
#endif /* USE_MPI */

#ifdef USE_MPI
static int mpi_node_count;
#endif

#ifdef USE_THREADS  /* user want threads */
#error Threading (USE_THREADS) support has been removed for very poor efficiency. Use MPI/SSH grid instead.
#endif


void   mcset_ncount(unsigned long long count);    /* wrapper to get mcncount */
unsigned long long int mcget_ncount(void);            /* wrapper to set mcncount */
unsigned long long mcget_run_num(void);           /* wrapper to get mcrun_num=0:mcncount */


/* Following part is only embedded when not redundant with mccode.h ========= */

#ifndef MCCODE_H

#ifndef NOSIGNALS
#include <signal.h>
char  *mcsig_message;
#define SIG_MESSAGE(msg) mcsig_message=(char *)(msg);
#else
#define SIG_MESSAGE
#endif /* !NOSIGNALS */


/* Useful macros and constants ============================================== */

#ifndef FLT_MAX
#define FLT_MAX         3.40282347E+38F /* max decimal value of a "float" */
#endif

#ifndef MIN
#define MIN(a, b)  (((a) < (b)) ? (a) : (b))
#endif
#ifndef MAX
#define MAX(a, b)  (((a) > (b)) ? (a) : (b))
#endif
#ifndef SQR
#define SQR(x) ( (x) * (x) )
#endif
#ifndef SIGN
#define SIGN(x) (((x)>0.0)?(1):(-1))
#endif

#ifndef PI
# ifdef M_PI
#  define PI M_PI
# else
/* When using c99 in the CFLAGS, some of these consts
   are lost... Perhaps we should in fact include everything from
   https://www.gnu.org/software/libc/manual/html_node/Mathematical-Constants.html
*/
#  define PI 3.14159265358979323846
#  define M_PI PI
#  define M_PI_2 M_PI/2.0
#  define M_PI_4 M_PI/4.0
#  define M_1_PI 1.0/M_PI
#  define M_2_PI 2*M_1_PI
# endif
#endif

#define RAD2MIN  ((180*60)/PI)
#define MIN2RAD  (PI/(180*60))
#define DEG2RAD  (PI/180)
#define RAD2DEG  (180/PI)
#define FWHM2RMS 0.424660900144    /* Convert between full-width-half-max and */
#define RMS2FWHM 2.35482004503     /* root-mean-square (standard deviation) */
#define HBAR     1.05457168e-34    /* [Js] h bar Planck constant CODATA 2002 */
#define MNEUTRON 1.67492728e-27    /* [kg] mass of neutron CODATA 2002 */
#define GRAVITY  9.81              /* [m/s^2] gravitational acceleration */
#define NA       6.02214179e23     /* [#atoms/g .mole] Avogadro's number*/


/* wrapper to get absolute and relative position of comp */
/* mccomp_posa and mccomp_posr are defined in McStas generated C code */
#define POS_A_COMP_INDEX(index) (instrument->_position_absolute[index])
#define POS_R_COMP_INDEX(index) (instrument->_position_relative[index])

/* setting parameters based MC_GETPAR */
#define MC_GETPAR3(type, compname, par) \
    &( ((_class_ ## type ##_parameters *) _getvar_parameters(compname))->par )
/* the body of this function depends on component instances, and is cogen'd */
void* _getvar_parameters(char* compname);

#define INSTRUMENT_GETPAR(par) (instrument->_parameters._ ## par)

/* Current component name, index, position and orientation */
/* These macros work because, using class-based functions, "comp" is usually
*  the local variable of the active/current component. */
#define INDEX_CURRENT_COMP (_comp->_index)
#define NAME_CURRENT_COMP (_comp->_name)
#define TYPE_CURRENT_COMP (_comp->_type)
#define POS_A_CURRENT_COMP (_comp->_position_absolute)
#define POS_R_CURRENT_COMP (_comp->_position_relative)
#define ROT_A_CURRENT_COMP (_comp->_rotation_absolute)
#define ROT_R_CURRENT_COMP (_comp->_rotation_relative)

#define NAME_INSTRUMENT (instrument->_name)


/* MCDISPLAY/trace and debugging message sent to stdout */
#ifdef MC_TRACE_ENABLED
#define DEBUG
#endif

#ifdef DEBUG
#define DEBUG_INSTR() if(!mcdotrace); else { printf("INSTRUMENT:\n"); printf("Instrument '%s' (%s)\n", instrument_name, instrument_source); }
#define DEBUG_COMPONENT(name,c,t) if(!mcdotrace); else {\
  printf("COMPONENT: \"%s\"\n" \
         "POS: %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g\n", \
         name, c.x, c.y, c.z, t[0][0], t[0][1], t[0][2], \
         t[1][0], t[1][1], t[1][2], t[2][0], t[2][1], t[2][2]); \
  printf("Component %30s AT (%g,%g,%g)\n", name, c.x, c.y, c.z); \
  }
#define DEBUG_INSTR_END() if(!mcdotrace); else printf("INSTRUMENT END:\n");
#define DEBUG_ENTER() if(!mcdotrace); else printf("ENTER:\n");
#define DEBUG_COMP(c) if(!mcdotrace); else printf("COMP: \"%s\"\n", c);
#define DEBUG_LEAVE() if(!mcdotrace); else printf("LEAVE:\n");
#define DEBUG_ABSORB() if(!mcdotrace); else printf("ABSORB:\n");
#else
#define DEBUG_INSTR()
#define DEBUG_COMPONENT(name,c,t)
#define DEBUG_INSTR_END()
#define DEBUG_ENTER()
#define DEBUG_COMP(c)
#define DEBUG_LEAVE()
#define DEBUG_ABSORB()
#endif

// mcDEBUG_STATE and mcDEBUG_SCATTER are defined by mcstas-r.h and mcxtrace-r.h



#ifdef TEST
#define test_printf printf
#else
#define test_printf while(0) printf
#endif

/* send MCDISPLAY message to stdout to show gemoetry */
void mcdis_magnify(char *what);
void mcdis_line(double x1, double y1, double z1,
                double x2, double y2, double z2);
void mcdis_dashed_line(double x1, double y1, double z1,
		       double x2, double y2, double z2, int n);
void mcdis_multiline(int count, ...);
void mcdis_rectangle(char* plane, double x, double y, double z,
		     double width, double height);
void mcdis_box(double x, double y, double z,
	       double width, double height, double length);
void mcdis_circle(char *plane, double x, double y, double z, double r);
void mcdis_Circle(double x, double y, double z, double r, double nx, double ny, double nz);
void mcdis_cylinder( double x, double y, double z,
        double r, double height, int N, double nx, double ny, double nz);
void mcdis_sphere(double x, double y, double z, double r, int N);

#  define MC_RAND_MAX ((unsigned long)0xffffffff)
/* selection of random number generator. default is MT */
#ifndef USE_PGI
/* "Mersenne Twister", by Makoto Matsumoto and Takuji Nishimura. */
#include <math.h>
#  define random mt_random
#  define srandom mt_srandom

#  define rand01 rand01_cpu
#  define randpm1 randpm1_cpu
#  define rand0max rand0max_cpu
#  define randminmax randminmax_cpu
#  define randnorm randnorm_cpu
#  define randtriangle randtriangle_cpu

#else
/* Use CUDA rand algo - we should check for availability
   of CUDA libs - and otherwise fail */
#  include <openacc_curand.h>
#  include <accelmath.h>
#  define random() /* TODO: set this */
#  define srandom mt_srandom

#  define rand01() curand_uniform(&_particle->MCRANDstate)
#  define randpm1() randpm1_gpu(&_particle->MCRANDstate)
#  define rand0max(p1) rand0max_gpu(p1, &_particle->MCRANDstate)
#  define randminmax(p1, p2) randminmax_gpu(p1, p2, &_particle->MCRANDstate)
#  define randnorm() curand_normal(&_particle->MCRANDstate)
#  define randtriangle() randtriangle_gpu(&_particle->MCRANDstate)

#  define printf_sys printf
#  define fprintf_sys fprintf
#endif


/*
* Random number generation.
*
* On the CPU, random alg. Marsenne Twister iteration state is handled
* internally. On the GPU however, every ray needs to keep track of its iteration
* state. The state variable is tied to the particle truct instance.
*
* This means that GPU random number functions need access to the particle
* context. Thus,separate functions exist for the CPU/GPU, and universally
* used symbols like "randnorm" and "rand01" are set as defines.
*/

// random number generation (MT on the CPU)
typedef int mc_int32_t;
mc_int32_t mc_random(void);
void mc_srandom (unsigned int x);
unsigned long mt_random(void);
void mt_srandom (unsigned long x);

// GPUs don't have i/o
/*int printf_GPU(const char *format, ...);
  int fprintf_GPU(FILE *stream, const char *format, ...);*/

// random number genration CPU
double rand01_cpu();
double randpm1_cpu();
double rand0max_cpu(double max);
double randminmax_cpu(double min, double max);
double randnorm_cpu(void);
double randtriangle_cpu(void);

// random number generation GPU
#ifdef USE_PGI
double randpm1_gpu(curandState_t* state);
double rand0max_gpu(double max, curandState_t* state);
double randminmax_gpu(double min, double max, curandState_t* state);
double randtriangle_gpu(curandState_t* state);
#endif

#ifdef USE_OPENCL
#include "opencl-lib.h"
#include "opencl-lib.c"
#endif

#ifndef DANSE
int init(void);
int raytrace(_class_particle*);
int save(FILE *);
int finally(void);
int display(void);
#endif

/* simple vector algebra ==================================================== */
#define vec_prod(x, y, z, x1, y1, z1, x2, y2, z2) \
	vec_prod_func(&x, &y, &z, x1, y1, z1, x2, y2, z2)
mcstatic void vec_prod_func(double *x, double *y, double *z,
		double x1, double y1, double z1, double x2, double y2, double z2);

mcstatic double scalar_prod(
		double x1, double y1, double z1, double x2, double y2, double z2);

mcstatic void norm_func(double *x, double *y, double *z);
#define NORM(x,y,z)	norm_func(&x, &y, &z)


void normal_vec(double *nx, double *ny, double *nz,
    double x, double y, double z);

/**
 * Rotate the vector vx,vy,vz psi radians around the vector ax,ay,az
 * and put the result in x,y,z.
 */
#define rotate(x, y, z, vx, vy, vz, phi, ax, ay, az) \
  do { \
    double mcrt_tmpx = (ax), mcrt_tmpy = (ay), mcrt_tmpz = (az); \
    double mcrt_vp, mcrt_vpx, mcrt_vpy, mcrt_vpz; \
    double mcrt_vnx, mcrt_vny, mcrt_vnz, mcrt_vn1x, mcrt_vn1y, mcrt_vn1z; \
    double mcrt_bx, mcrt_by, mcrt_bz; \
    double mcrt_cos, mcrt_sin; \
    NORM(mcrt_tmpx, mcrt_tmpy, mcrt_tmpz); \
    mcrt_vp = scalar_prod((vx), (vy), (vz), mcrt_tmpx, mcrt_tmpy, mcrt_tmpz); \
    mcrt_vpx = mcrt_vp*mcrt_tmpx; \
    mcrt_vpy = mcrt_vp*mcrt_tmpy; \
    mcrt_vpz = mcrt_vp*mcrt_tmpz; \
    mcrt_vnx = (vx) - mcrt_vpx; \
    mcrt_vny = (vy) - mcrt_vpy; \
    mcrt_vnz = (vz) - mcrt_vpz; \
    vec_prod(mcrt_bx, mcrt_by, mcrt_bz, \
             mcrt_tmpx, mcrt_tmpy, mcrt_tmpz, mcrt_vnx, mcrt_vny, mcrt_vnz); \
    mcrt_cos = cos((phi)); mcrt_sin = sin((phi)); \
    mcrt_vn1x = mcrt_vnx*mcrt_cos + mcrt_bx*mcrt_sin; \
    mcrt_vn1y = mcrt_vny*mcrt_cos + mcrt_by*mcrt_sin; \
    mcrt_vn1z = mcrt_vnz*mcrt_cos + mcrt_bz*mcrt_sin; \
    (x) = mcrt_vpx + mcrt_vn1x; \
    (y) = mcrt_vpy + mcrt_vn1y; \
    (z) = mcrt_vpz + mcrt_vn1z; \
  } while(0)

/**
 * Mirror (xyz) in the plane given by the point (rx,ry,rz) and normal (nx,ny,nz)
 *
 * TODO: This define is seemingly never used...
 */
#define mirror(x,y,z,rx,ry,rz,nx,ny,nz) \
  do { \
    double mcrt_tmpx= (nx), mcrt_tmpy = (ny), mcrt_tmpz = (nz); \
    double mcrt_tmpt; \
    NORM(mcrt_tmpx, mcrt_tmpy, mcrt_tmpz); \
    mcrt_tmpt=scalar_prod((rx),(ry),(rz),mcrt_tmpx,mcrt_tmpy,mcrt_tmpz); \
    (x) = rx -2 * mcrt_tmpt*mcrt_rmpx; \
    (y) = ry -2 * mcrt_tmpt*mcrt_rmpy; \
    (z) = rz -2 * mcrt_tmpt*mcrt_rmpz; \
  } while (0)

Coords coords_set(MCNUM x, MCNUM y, MCNUM z);
Coords coords_get(Coords a, MCNUM *x, MCNUM *y, MCNUM *z);
Coords coords_add(Coords a, Coords b);
Coords coords_sub(Coords a, Coords b);
Coords coords_neg(Coords a);
Coords coords_scale(Coords b, double scale);
double coords_sp(Coords a, Coords b);
Coords coords_xp(Coords b, Coords c);
double coords_len(Coords a);
void   coords_print(Coords a);
mcstatic void coords_norm(Coords* c);

void rot_set_rotation(Rotation t, double phx, double phy, double phz);
int  rot_test_identity(Rotation t);
void rot_mul(Rotation t1, Rotation t2, Rotation t3);
void rot_copy(Rotation dest, Rotation src);
void rot_transpose(Rotation src, Rotation dst);
Coords rot_apply(Rotation t, Coords a);

void mccoordschange(Coords a, Rotation t, _class_particle *particle);
void mccoordschange_polarisation(Rotation t, double *sx, double *sy, double *sz);

double mcestimate_error(double N, double p1, double p2);
void mcreadparams(void);

/* this is now in mcstas-r.h and mcxtrace-r.h as the number of state parameters is no longer equal*/
_class_particle mcgenstate(void);

/* trajectory/shape intersection routines */
int inside_rectangle(double, double, double, double);
int box_intersect(double *dt_in, double *dt_out, double x, double y, double z,
      double vx, double vy, double vz, double dx, double dy, double dz);
int cylinder_intersect(double *t0, double *t1, double x, double y, double z,
      double vx, double vy, double vz, double r, double h);
int sphere_intersect(double *t0, double *t1, double x, double y, double z,
      double vx, double vy, double vz, double r);
/* second order equation roots */
int solve_2nd_order(double *t1, double *t2,
      double A,  double B,  double C);


/* random vector generation to shape */

/* jg-20191003: Gpu versions need to apppend _particle for gpu rng functions,
 * which depend on the particle struct MCRANDstate variable, but this has to be
 * transparent to component use of randvec shape functions. */

// interface (shared)
#  define randvec_target_sphere randvec_target_circle
#  define randvec_target_rect(p0,p1,p2,p3,p4,p5,p6,p7,p8,p9) \
        randvec_target_rect_real(p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,0,0,0,1)
#ifndef USE_PGI
// interface (cpu)
#  define randvec_target_circle       randvec_target_circle_cpu
#  define randvec_target_rect_angular randvec_target_rect_angular_cpu
#  define randvec_target_rect_real    randvec_target_rect_real_cpu
// headers (cpu)
void randvec_target_circle_cpu(double *xo, double *yo, double *zo,
      double *solid_angle, double xi, double yi, double zi, double radius);
void randvec_target_rect_angular_cpu(double *xo, double *yo, double *zo,
      double *solid_angle,
      double xi, double yi, double zi, double height, double width, Rotation A);
void randvec_target_rect_real_cpu(double *xo, double *yo, double *zo, double *solid_angle,
      double xi, double yi, double zi, double height, double width, Rotation A,
      double lx, double ly, double lz, int order);
#else
// interface (gpu, appending _particle)
#  define randvec_target_circle(xo, yo, zo, solid_angle, xi, yi, zi, radius) \
        randvec_target_circle_gpu(xo, yo, zo, solid_angle, xi, yi, zi, radius, _particle)
#  define randvec_target_rect_angular(xo, yo, zo, solid_angle, xi, yi, zi, height, width, Rotation A) \
        randvec_target_rect_angular_gpu(xo, yo, zo, solid_angle, xi, yi, zi, height, width, Rotation A, _particle)
#  define randvec_target_rect_real(xo, yo, zo, solid_angle, xi, yi, zi, height, width, A, lx, ly, lz, order) \
        randvec_target_rect_real_gpu(xo, yo, zo, solid_angle, xi, yi, zi, height, width, A, lx, ly, lz, order, _particle)
// heaaders (gpu)
void randvec_target_circle_gpu(double *xo, double *yo, double *zo,
      double *solid_angle, double xi, double yi, double zi, double radius,
      _class_particle* _particle);
void randvec_target_rect_angular_gpu(double *xo, double *yo, double *zo,
      double *solid_angle, double xi, double yi, double zi, double height,
      double width, Rotation A, _class_particle* _particle);
void randvec_target_rect_real_gpu(double *xo, double *yo, double *zo, double *solid_angle,
      double xi, double yi, double zi, double height, double width, Rotation A,
      double lx, double ly, double lz, int order, _class_particle* _particle);
#endif

/* this is the main() */
int mccode_main(int argc, char *argv[]);


#endif /* !MCCODE_H */

#ifndef MCCODE_R_IO_H
#define MCCODE_R_IO_H "$Revision$"

#if (USE_NEXUS == 0)
#undef USE_NEXUS
#endif

#ifndef CHAR_BUF_LENGTH
#define CHAR_BUF_LENGTH 1024
#endif


/* I/O section part ========================================================= */

/* ========================================================================== */

/*                               MCCODE_R_IO_C                                */

/* ========================================================================== */


/* main DETECTOR structure which stores most information to write to data files */
struct mcdetector_struct {
  char   filename[CHAR_BUF_LENGTH];   /* file name of monitor */
  char   position[CHAR_BUF_LENGTH];   /* position of detector component */
  char   component[CHAR_BUF_LENGTH];  /* component instance name */
  char   instrument[CHAR_BUF_LENGTH]; /* instrument name */
  char   type[CHAR_BUF_LENGTH];       /* data type, e.g. 0d, 1d, 2d, 3d */
  char   user[CHAR_BUF_LENGTH];       /* user name, e.g. HOME */
  char   date[CHAR_BUF_LENGTH];       /* date of simulation end/write time */
  char   title[CHAR_BUF_LENGTH];      /* title of detector */
  char   xlabel[CHAR_BUF_LENGTH];     /* X axis label */
  char   ylabel[CHAR_BUF_LENGTH];     /* Y axis label */
  char   zlabel[CHAR_BUF_LENGTH];     /* Z axis label */
  char   xvar[CHAR_BUF_LENGTH];       /* X variable name */
  char   yvar[CHAR_BUF_LENGTH];       /* Y variable name */
  char   zvar[CHAR_BUF_LENGTH];       /* Z variable name */
  char   ncount[CHAR_BUF_LENGTH];     /* number of events initially generated */
  char   limits[CHAR_BUF_LENGTH];     /* X Y Z limits, e.g. [xmin xmax ymin ymax zmin zmax] */
  char   variables[CHAR_BUF_LENGTH];  /* variables written into data block */
  char   statistics[CHAR_BUF_LENGTH]; /* center, mean and half width along axis */
  char   signal[CHAR_BUF_LENGTH];     /* min max and mean of signal (data block) */
  char   values[CHAR_BUF_LENGTH];     /* integrated values e.g. [I I_err N] */
  double xmin,xmax;                   /* min max of axes */
  double ymin,ymax;
  double zmin,zmax;
  double intensity;                   /* integrated values for data block */
  double error;
  double events;
  double min;                         /* statistics for data block */
  double max;
  double mean;
  double centerX;                     /* statistics for axes */
  double halfwidthX;
  double centerY;
  double halfwidthY;
  int    rank;                        /* dimensionaly of monitor, e.g. 0 1 2 3 */
  char   istransposed;                /* flag to transpose matrix for some formats */

  long   m,n,p;                       /* dimensions of data block and along axes */
  long   date_l;                      /* same as date, but in sec since 1970 */

  double *p0, *p1, *p2;               /* pointers to saved data, NULL when freed */
  char   format[CHAR_BUF_LENGTH];    /* format for file generation */
};

typedef struct mcdetector_struct MCDETECTOR;

static   char *dirname             = NULL;      /* name of output directory */
static   char *siminfo_name        = "mccode";  /* default output sim file name */
char    *mcformat                    = NULL;      /* NULL (default) or a specific format */

/* file I/O definitions and function prototypes */

#ifndef MC_EMBEDDED_RUNTIME /* the mcstatic variables (from mccode-r.c) */
extern FILE * siminfo_file;     /* handle to the output siminfo file */
extern int    mcgravitation;      /* flag to enable gravitation */
#pragma acc declare create ( mcgravitation )
extern int    mcdotrace;          /* flag to print MCDISPLAY messages */
#else
mcstatic FILE *siminfo_file        = NULL;
#endif

/* I/O function prototypes ================================================== */

// from msysgit: https://code.google.com/p/msysgit/source/browse/compat/strcasestr.c
char *strcasestr(const char *haystack, const char *needle);

/* output functions */
MCDETECTOR mcdetector_out_0D(char *t, double p0, double p1, double p2, char *c, Coords pos);
MCDETECTOR mcdetector_out_1D(char *t, char *xl, char *yl,
                  char *xvar, double x1, double x2, long n,
                  double *p0, double *p1, double *p2, char *f, char *c, Coords pos);
MCDETECTOR mcdetector_out_2D(char *t, char *xl, char *yl,
                  double x1, double x2, double y1, double y2, long m,
                  long n, double *p0, double *p1, double *p2, char *f,
                  char *c, Coords pos);
MCDETECTOR mcdetector_out_list(char *t, char *xl, char *yl,
                  long m, long n,
                  double *p1, char *f,
                  char *c, Coords posa);

/* wrappers to output functions, that automatically set NAME and POSITION */
#define DETECTOR_OUT(p0,p1,p2) mcdetector_out_0D(NAME_CURRENT_COMP,p0,p1,p2,NAME_CURRENT_COMP,POS_A_CURRENT_COMP)
#define DETECTOR_OUT_0D(t,p0,p1,p2) mcdetector_out_0D(t,p0,p1,p2,NAME_CURRENT_COMP,POS_A_CURRENT_COMP)
#define DETECTOR_OUT_1D(t,xl,yl,xvar,x1,x2,n,p0,p1,p2,f) \
     mcdetector_out_1D(t,xl,yl,xvar,x1,x2,n,p0,p1,p2,f,NAME_CURRENT_COMP,POS_A_CURRENT_COMP)
#define DETECTOR_OUT_2D(t,xl,yl,x1,x2,y1,y2,m,n,p0,p1,p2,f) \
     mcdetector_out_2D(t,xl,yl,x1,x2,y1,y2,m,n,p0,p1,p2,f,NAME_CURRENT_COMP,POS_A_CURRENT_COMP)

#ifdef USE_NEXUS
#include "napi.h"
NXhandle nxhandle;
#endif

#endif /* ndef MCCODE_R_IO_H */

#endif /* MCCODE_R_H */
/* End of file "mccode-r.h". */

/* embedding file "mcstas-r.h" */

/*******************************************************************************
*
* McStas, neutron ray-tracing package
*         Copyright (C) 1997-2009, All rights reserved
*         Risoe National Laboratory, Roskilde, Denmark
*         Institut Laue Langevin, Grenoble, France
*
* Runtime: share/mcstas-r.h
*
* %Identification
* Written by: KN
* Date:    Aug 29, 1997
* Release: McStas X.Y
* Version: $Revision$
*
* Runtime system header for McStas.
*
* In order to use this library as an external library, the following variables
* and macros must be declared (see details in the code)
*
*   struct mcinputtable_struct mcinputtable[];
*   int mcnumipar;
*   char instrument_name[], instrument_source[];
*   int traceenabled, defaultmain;
*   extern MCNUM  mccomp_storein[];
*   extern MCNUM  instrument.counter_AbsorbProp[];
*   extern MCNUM  mcScattered;
*   #define MCCODE_STRING "the McStas version"
*
* Usage: Automatically embbeded in the c code.
*
* $Id$
*
*******************************************************************************/

#ifndef MCSTAS_R_H
#define MCSTAS_R_H "$Revision$"

/* Following part is only embedded when not redundent with mcstas.h ========= */

#ifndef MCCODE_H

#define AA2MS    629.622368        /* Convert k[1/AA] to v[m/s] */
#define MS2AA    1.58825361e-3     /* Convert v[m/s] to k[1/AA] */
#define K2V      AA2MS
#define V2K      MS2AA
#define Q2V      AA2MS
#define V2Q      MS2AA
#define SE2V     437.393377        /* Convert sqrt(E)[meV] to v[m/s] */
#define VS2E     5.22703725e-6     /* Convert (v[m/s])**2 to E[meV] */

#define SCATTER0 do {DEBUG_SCATTER(); SCATTERED++;} while(0)
#define SCATTER SCATTER0

#define JUMPTOCOMP(comp) mcneutron->_index = INDEX_COMP(comp);

#define MAGNET_ON \
  do { \
    mcMagnet = 1; \
  } while(0)

#define MAGNET_OFF \
  do { \
    mcMagnet = 0; \
  } while(0)

#define ALLOW_BACKPROP \
  do { \
    mcallowbackprop = 1; \
  } while(0)

#define DISALLOW_BACKPROP \
  do { \
    mcallowbackprop = 0; \
  } while(0)

#define PROP_MAGNET(dt) \
  do { \
  } while (0)
    /* change coordinates from local system to magnet system */
/*    Rotation rotLM, rotTemp; \
      Coords   posLM = coords_sub(POS_A_CURRENT_COMP, mcMagnetPos); \
      rot_transpose(ROT_A_CURRENT_COMP, rotTemp); \
      rot_mul(rotTemp, mcMagnetRot, rotLM); \
      mcMagnetPrecession(x, y, z, t, vx, vy, vz, \
               &sx, &sy, &sz, dt, posLM, rotLM); \
      } while(0)
*/

#define mcPROP_DT(dt) \
  do { \
    if (mcMagnet && dt > 0) PROP_MAGNET(dt);\
    x += vx*(dt); \
    y += vy*(dt); \
    z += vz*(dt); \
    t += (dt); \
    if (isnan(p) || isinf(p)) { instrument->counter_AbsorbProp[INDEX_CURRENT_COMP]++; ABSORB; }\
  } while(0)

/* ADD: E. Farhi, Aug 6th, 2001 PROP_GRAV_DT propagation with acceleration */
#define PROP_GRAV_DT(dt, Ax, Ay, Az) \
  do { \
    if(dt < 0 && mcallowbackprop == 0) { instrument->counter_AbsorbProp[INDEX_CURRENT_COMP]++; ABSORB; }\
    if (mcMagnet) /*printf("Spin precession gravity\n")*/; \
    x  += vx*(dt) + (Ax)*(dt)*(dt)/2; \
    y  += vy*(dt) + (Ay)*(dt)*(dt)/2; \
    z  += vz*(dt) + (Az)*(dt)*(dt)/2; \
    vx += (Ax)*(dt); \
    vy += (Ay)*(dt); \
    vz += (Az)*(dt); \
    t  += (dt); \
    DISALLOW_BACKPROP;\
  } while(0)


#define PROP_DT(dt) \
  do { \
    if(dt < 0) { RESTORE=1; ABSORB; }; \
    if (mcgravitation) { Coords mcLocG; double mc_gx, mc_gy, mc_gz; \
    mcLocG = rot_apply(ROT_A_CURRENT_COMP, coords_set(0,-GRAVITY,0)); \
    coords_get(mcLocG, &mc_gx, &mc_gy, &mc_gz); \
    PROP_GRAV_DT(dt, mc_gx, mc_gy, mc_gz); } \
    else mcPROP_DT(dt); \
    DISALLOW_BACKPROP;\
  } while(0)


#define PROP_Z0 \
  do { \
    if (mcgravitation) { Coords mcLocG; int mc_ret; \
    double mc_dt, mc_gx, mc_gy, mc_gz; \
    mcLocG = rot_apply(ROT_A_CURRENT_COMP, coords_set(0,-GRAVITY,0)); \
    coords_get(mcLocG, &mc_gx, &mc_gy, &mc_gz); \
    mc_ret = solve_2nd_order(&mc_dt, NULL, -mc_gz/2, -vz, -z); \
    if (mc_ret && mc_dt>=0) PROP_GRAV_DT(mc_dt, mc_gx, mc_gy, mc_gz); \
    else { if (mcallowbackprop ==0) {instrument->counter_AbsorbProp[INDEX_CURRENT_COMP]++; ABSORB; }}; }\
    else mcPROP_Z0; \
    DISALLOW_BACKPROP;\
  } while(0)

#define mcPROP_Z0 \
  do { \
    double mc_dt; \
    if(vz == 0) { instrument->counter_AbsorbProp[INDEX_CURRENT_COMP]++; ABSORB; }; \
    mc_dt = -z/vz; \
    if(mc_dt < 0 && mcallowbackprop == 0) { instrument->counter_AbsorbProp[INDEX_CURRENT_COMP]++; ABSORB; }; \
    mcPROP_DT(mc_dt); \
    z = 0; \
    DISALLOW_BACKPROP;\
  } while(0)

#define PROP_X0 \
  do { \
    if (mcgravitation) { Coords mcLocG; int mc_ret; \
    double mc_dt, mc_gx, mc_gy, mc_gz; \
    mcLocG = rot_apply(ROT_A_CURRENT_COMP, coords_set(0,-GRAVITY,0)); \
    coords_get(mcLocG, &mc_gx, &mc_gy, &mc_gz); \
    mc_ret = solve_2nd_order(&mc_dt, NULL, -mc_gx/2, -vx, -x); \
    if (mc_ret && mc_dt>=0) PROP_GRAV_DT(mc_dt, mc_gx, mc_gy, mc_gz); \
    else { if (mcallowbackprop ==0) {instrument->counter_AbsorbProp[INDEX_CURRENT_COMP]++; ABSORB; }}; }\
    else mcPROP_X0; \
    DISALLOW_BACKPROP;\
  } while(0)

#define mcPROP_X0 \
  do { \
    double mc_dt; \
    if(vx == 0) { instrument->counter_AbsorbProp[INDEX_CURRENT_COMP]++; ABSORB; }; \
    mc_dt = -x/vx; \
    if(mc_dt < 0 && mcallowbackprop == 0) { instrument->counter_AbsorbProp[INDEX_CURRENT_COMP]++; ABSORB; }; \
    mcPROP_DT(mc_dt); \
    x = 0; \
    DISALLOW_BACKPROP;\
  } while(0)

#define PROP_Y0 \
  do { \
    if (mcgravitation) { Coords mcLocG; int mc_ret; \
    double mc_dt, mc_gx, mc_gy, mc_gz; \
    mcLocG = rot_apply(ROT_A_CURRENT_COMP, coords_set(0,-GRAVITY,0)); \
    coords_get(mcLocG, &mc_gx, &mc_gy, &mc_gz); \
    mc_ret = solve_2nd_order(&mc_dt, NULL, -mc_gy/2, -vy, -y); \
    if (mc_ret && mc_dt>=0) PROP_GRAV_DT(mc_dt, mc_gx, mc_gy, mc_gz); \
    else { if (mcallowbackprop ==0) {instrument->counter_AbsorbProp[INDEX_CURRENT_COMP]++; ABSORB; }}; }\
    else mcPROP_Y0; \
    DISALLOW_BACKPROP;\
  } while(0)


#define mcPROP_Y0 \
  do { \
    double mc_dt; \
    if(vy == 0) { instrument->counter_AbsorbProp[INDEX_CURRENT_COMP]++; ABSORB; }; \
    mc_dt = -y/vy; \
    if(mc_dt < 0 && mcallowbackprop == 0) { instrument->counter_AbsorbProp[INDEX_CURRENT_COMP]++; ABSORB; }; \
    mcPROP_DT(mc_dt); \
    y = 0; \
    DISALLOW_BACKPROP; \
  } while(0)

#pragma acc routine seq
_class_particle mcsetstate(double x, double y, double z, double vx, double vy, double vz,
                double t, double sx, double sy, double sz, double p);

#ifdef DEBUG

#define DEBUG_STATE() if(!mcdotrace); else \
  printf("STATE: %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g\n", \
         x,y,z,vx,vy,vz,t,sx,sy,sz,p);
#define DEBUG_SCATTER() if(!mcdotrace); else \
  printf("SCATTER: %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g\n", \
         x,y,z,vx,vy,vz,t,sx,sy,sz,p);

#else

#define DEBUG_STATE()
#define DEBUG_SCATTER()

#endif

#endif /* !MCCODE_H */

#endif /* MCSTAS_R_H */
/* End of file "mcstas-r.h". */

/* embedding file "mccode-r.c" */

/*******************************************************************************
*
* McCode, neutron/xray ray-tracing package
*         Copyright (C) 1997-2009, All rights reserved
*         Risoe National Laboratory, Roskilde, Denmark
*         Institut Laue Langevin, Grenoble, France
*
* Runtime: share/mccode-r.c
*
* %Identification
* Written by: KN
* Date:    Aug 29, 1997
* Release: McStas X.Y/McXtrace X.Y
* Version: $Revision$
*
* Runtime system for McStas and McXtrace.
* Embedded within instrument in runtime mode.
* Contains SECTIONS:
*   MPI handling (sum, send, recv)
*   format definitions
*   I/O
*   mcdisplay support
*   random numbers
*   coordinates handling
*   vectors math (solve 2nd order, normals, randvec...)
*   parameter handling
*   signal and main handlers
*
* Usage: Automatically embbeded in the c code whenever required.
*
* $Id$
*
*******************************************************************************/

/*******************************************************************************
* The I/O format definitions and functions
*******************************************************************************/


/** Include header files to avoid implicit declarations (not allowed on LLVM) */
#include <ctype.h>
#include <sys/types.h>

// UNIX specific headers (non-Windows)
#if defined(__unix__) || defined(__APPLE__)
#include <unistd.h>
#include <sys/stat.h>
#endif


#ifndef DANSE
#ifdef MC_ANCIENT_COMPATIBILITY
int traceenabled = 0;
int defaultmain  = 0;
#endif
/* else defined directly in the McCode generated C code */

static   long mcseed                 = 0; /* seed for random generator */
#pragma acc declare create ( mcseed )
static   long mcstartdate            = 0; /* start simulation time */
static   int  mcdisable_output_files = 0; /* --no-output-files */
mcstatic int  mcgravitation          = 0; /* use gravitation flag, for PROP macros */
#pragma acc declare create ( mcgravitation )
int      mcMagnet                    = 0; /* magnet stack flag */
#pragma acc declare create ( mcMagnet )
mcstatic int  mcdotrace              = 0; /* flag for --trace and messages for DISPLAY */
int      mcallowbackprop             = 0;         /* flag to enable negative/backprop */
#pragma acc declare create ( mcallowbackprop )

/* Number of particle histories to simulate. */
#ifdef NEUTRONICS
mcstatic unsigned long long int mcncount             = 1;
mcstatic unsigned long long int mcrun_num            = 0;
#else
#ifdef MCDEFAULT_NCOUNT
mcstatic unsigned long long int mcncount             = MCDEFAULT_NCOUNT;
#else
mcstatic unsigned long long int mcncount             = 1000000;
#endif
#pragma acc declare create ( mcncount )
mcstatic unsigned long long int mcrun_num            = 0;
#endif /* NEUTRONICS */

#else
#include "mcstas-globals.h"
#endif /* !DANSE */


/* SECTION: Dynamic Arrays ================================================== */

#pragma acc routine seq
DArray1d create_darr1d(int n){
  DArray1d arr2d;
  //arr2d = calloc(n, sizeof(double));
  return arr2d;
}
#pragma acc routine seq
void destroy_darr1d(DArray1d a){
  //free(a);
}

#pragma acc routine seq
DArray2d create_darr2d(int nx, int ny){
  DArray2d arr2d;
  /*arr2d = calloc(nx, sizeof(double *));

  double *p1;
  p1 = calloc(nx*ny, sizeof(double));

  int i;
  for (i=0; i<nx; i++){
    arr2d[i] = &(p1[i*ny]);
  }*/
  return arr2d;
}
#pragma acc routine seq
void destroy_darr2d(DArray2d a){
  //free(a[0]);
  //free(a);
}

#pragma acc routine seq
DArray3d create_darr3d(int nx, int ny, int nz){
  DArray3d arr3d;
  /*
  int i, j;

  // 1d
  arr3d = calloc(nx, sizeof(double **));

  // d2
  double **p1;
  p1 = calloc(nx*ny, sizeof(double *));

  for (i=0; i<nx; i++){
    arr3d[i] = &(p1[i*ny]);
  }

  // 3d
  double *p2;
  p2 = calloc(nx*ny*nz, sizeof(double));
  for (i=0; i<nx; i++){
    for (j=0; j<ny; j++){
      arr3d[i][j] = &(p2[(i*ny+j)*nz]);
    }
  }*/
  return arr3d;
}
#pragma acc routine seq
void destroy_darr3d(DArray3d a){
  //free(a[0][0]);
  //free(a[0]);
  //free(a);
}


/* SECTION: MPI handling ==================================================== */

#ifdef USE_MPI
/* MPI rank */
static int mpi_node_rank;
static int mpi_node_root = 0;


/*******************************************************************************
* mc_MPI_Reduce: Gathers arrays from MPI nodes using Reduce function.
*******************************************************************************/
int mc_MPI_Sum(double *sbuf, long count)
{
  if (!sbuf || count <= 0) return(MPI_SUCCESS); /* nothing to reduce */
  else {
    /* we must cut the buffer into blocks not exceeding the MPI max buffer size of 32000 */
    long   offset=0;
    double *rbuf=NULL;
    int    length=MPI_REDUCE_BLOCKSIZE; /* defined in mccode-r.h */
    int    i=0;
    rbuf = calloc(count, sizeof(double));
    if (!rbuf)
      exit(-fprintf(stderr, "Error: Out of memory %li (mc_MPI_Sum)\n", count*sizeof(double)));
    while (offset < count) {
      if (!length || offset+length > count-1) length=count-offset;
      else length=MPI_REDUCE_BLOCKSIZE;
      if (MPI_Allreduce((double*)(sbuf+offset), (double*)(rbuf+offset),
              length, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD) != MPI_SUCCESS)
        return MPI_ERR_COUNT;
      offset += length;
    }

    for (i=0; i<count; i++) sbuf[i] = rbuf[i];
    free(rbuf);
  }
  return MPI_SUCCESS;
} /* mc_MPI_Sum */

/*******************************************************************************
* mc_MPI_Send: Send array to MPI node by blocks to avoid buffer limit
*******************************************************************************/
int mc_MPI_Send(void *sbuf,
                  long count, MPI_Datatype dtype,
                  int dest)
{
  int dsize;
  long offset=0;
  int  tag=1;
  int  length=MPI_REDUCE_BLOCKSIZE; /* defined in mccode-r.h */

  if (!sbuf || count <= 0) return(MPI_SUCCESS); /* nothing to send */
  MPI_Type_size(dtype, &dsize);

  while (offset < count) {
    if (offset+length > count-1) length=count-offset;
    else length=MPI_REDUCE_BLOCKSIZE;
    if (MPI_Send((void*)(sbuf+offset*dsize), length, dtype, dest, tag++, MPI_COMM_WORLD) != MPI_SUCCESS)
      return MPI_ERR_COUNT;
    offset += length;
  }

  return MPI_SUCCESS;
} /* mc_MPI_Send */

/*******************************************************************************
* mc_MPI_Recv: Receives arrays from MPI nodes by blocks to avoid buffer limit
*             the buffer must have been allocated previously.
*******************************************************************************/
int mc_MPI_Recv(void *sbuf,
                  long count, MPI_Datatype dtype,
                  int source)
{
  int dsize;
  long offset=0;
  int  tag=1;
  int  length=MPI_REDUCE_BLOCKSIZE; /* defined in mccode-r.h */

  if (!sbuf || count <= 0) return(MPI_SUCCESS); /* nothing to recv */
  MPI_Type_size(dtype, &dsize);

  while (offset < count) {
    if (offset+length > count-1) length=count-offset;
    else length=MPI_REDUCE_BLOCKSIZE;
    if (MPI_Recv((void*)(sbuf+offset*dsize), length, dtype, source, tag++,
            MPI_COMM_WORLD, MPI_STATUS_IGNORE) != MPI_SUCCESS)
      return MPI_ERR_COUNT;
    offset += length;
  }

  return MPI_SUCCESS;
} /* mc_MPI_Recv */

#endif /* USE_MPI */

/* SECTION: parameters handling ============================================= */

/* Instrument input parameter type handling. */
/*******************************************************************************
* mcparm_double: extract double value from 's' into 'vptr'
*******************************************************************************/
static int
mcparm_double(char *s, void *vptr)
{
  char *p;
  double *v = (double *)vptr;

  if (!s) { *v = 0; return(1); }
  *v = strtod(s, &p);
  if(*s == '\0' || (p != NULL && *p != '\0') || errno == ERANGE)
    return 0;                        /* Failed */
  else
    return 1;                        /* Success */
}

/*******************************************************************************
* mcparminfo_double: display parameter type double
*******************************************************************************/
static char *
mcparminfo_double(char *parmname)
{
  return "double";
}

/*******************************************************************************
* mcparmerror_double: display error message when failed extract double
*******************************************************************************/
static void
mcparmerror_double(char *parm, char *val)
{
  fprintf(stderr, "Error: Invalid value '%s' for floating point parameter %s (mcparmerror_double)\n",
          val, parm);
}

/*******************************************************************************
* mcparmprinter_double: convert double to string
*******************************************************************************/
static void
mcparmprinter_double(char *f, void *vptr)
{
  double *v = (double *)vptr;
  sprintf(f, "%g", *v);
}

/*******************************************************************************
* mcparm_int: extract int value from 's' into 'vptr'
*******************************************************************************/
static int
mcparm_int(char *s, void *vptr)
{
  char *p;
  int *v = (int *)vptr;
  long x;

  if (!s) { *v = 0; return(1); }
  *v = 0;
  x = strtol(s, &p, 10);
  if(x < INT_MIN || x > INT_MAX)
    return 0;                        /* Under/overflow */
  *v = x;
  if(*s == '\0' || (p != NULL && *p != '\0') || errno == ERANGE)
    return 0;                        /* Failed */
  else
    return 1;                        /* Success */
}

/*******************************************************************************
* mcparminfo_int: display parameter type int
*******************************************************************************/
static char *
mcparminfo_int(char *parmname)
{
  return "int";
}

/*******************************************************************************
* mcparmerror_int: display error message when failed extract int
*******************************************************************************/
static void
mcparmerror_int(char *parm, char *val)
{
  fprintf(stderr, "Error: Invalid value '%s' for integer parameter %s (mcparmerror_int)\n",
          val, parm);
}

/*******************************************************************************
* mcparmprinter_int: convert int to string
*******************************************************************************/
static void
mcparmprinter_int(char *f, void *vptr)
{
  int *v = (int *)vptr;
  sprintf(f, "%d", *v);
}

/*******************************************************************************
* mcparm_string: extract char* value from 's' into 'vptr' (copy)
*******************************************************************************/
static int
mcparm_string(char *s, void *vptr)
{
  char **v = (char **)vptr;
  if (!s) { *v = NULL; return(1); }
  *v = (char *)malloc(strlen(s) + 1);
  if(*v == NULL)
  {
    exit(-fprintf(stderr, "Error: Out of memory %li (mcparm_string).\n", (long)strlen(s) + 1));
  }
  strcpy(*v, s);
  return 1;                        /* Success */
}

/*******************************************************************************
* mcparminfo_string: display parameter type string
*******************************************************************************/
static char *
mcparminfo_string(char *parmname)
{
  return "string";
}

/*******************************************************************************
* mcparmerror_string: display error message when failed extract string
*******************************************************************************/
static void
mcparmerror_string(char *parm, char *val)
{
  fprintf(stderr, "Error: Invalid value '%s' for string parameter %s (mcparmerror_string)\n",
          val, parm);
}

/*******************************************************************************
* mcparmprinter_string: convert string to string (including esc chars)
*******************************************************************************/
static void
mcparmprinter_string(char *f, void *vptr)
{
  char **v = (char **)vptr;
  char *p;

  if (!*v) { *f='\0'; return; }
  strcpy(f, "");
  for(p = *v; *p != '\0'; p++)
  {
    switch(*p)
    {
      case '\n':
        strcat(f, "\\n");
        break;
      case '\r':
        strcat(f, "\\r");
        break;
      case '"':
        strcat(f, "\\\"");
        break;
      case '\\':
        strcat(f, "\\\\");
        break;
      default:
        strncat(f, p, 1);
    }
  }
  /* strcat(f, "\""); */
} /* mcparmprinter_string */

/* now we may define the parameter structure, using previous functions */
static struct
  {
    int (*getparm)(char *, void *);
    char * (*parminfo)(char *);
    void (*error)(char *, char *);
    void (*printer)(char *, void *);
} mcinputtypes[] = {
  {
    mcparm_int, mcparminfo_int, mcparmerror_int,
    mcparmprinter_int
  }, {
    mcparm_string, mcparminfo_string, mcparmerror_string,
    mcparmprinter_string
  }, {
    mcparm_string, mcparminfo_string, mcparmerror_string,
    mcparmprinter_string
  }, {
    mcparm_double, mcparminfo_double, mcparmerror_double,
    mcparmprinter_double
  }, {
    mcparm_double, mcparminfo_double, mcparmerror_double,
    mcparmprinter_double
  }
};

/*******************************************************************************
* mcestimate_error: compute sigma from N,p,p2 in Gaussian large numbers approx
*******************************************************************************/
double mcestimate_error(double N, double p1, double p2)
{
  double pmean, n1;
  if(N <= 1)
    return p1;
  pmean = p1 / N;
  n1 = N - 1;
  /* Note: underflow may cause p2 to become zero; the fabs() below guards
     against this. */
  return sqrt((N/n1)*fabs(p2 - pmean*pmean));
}

double (*mcestimate_error_p)
  (double V2, double psum, double p2sum)=mcestimate_error;

/* ========================================================================== */

/*                               MCCODE_R_IO_C                                */

/* ========================================================================== */

#ifndef MCCODE_R_IO_C
#define MCCODE_R_IO_C "$Revision$"

/* SECTION: file i/o handling ================================================ */

#ifndef HAVE_STRCASESTR
// from msysgit: https://code.google.com/p/msysgit/source/browse/compat/strcasestr.c
char *strcasestr(const char *haystack, const char *needle)
{
  int nlen = strlen(needle);
  int hlen = strlen(haystack) - nlen + 1;
  int i;

  for (i = 0; i < hlen; i++) {
    int j;
    for (j = 0; j < nlen; j++) {
            unsigned char c1 = haystack[i+j];
            unsigned char c2 = needle[j];
            if (toupper(c1) != toupper(c2))
                    goto next;
    }
    return (char *) haystack + i;
  next:
    ;
  }
  return NULL;
}


#endif
#ifndef HAVE_STRCASECMP
int strcasecmp( const char *s1, const char *s2 )
{
  int c1, c2;
  do {
    c1 = tolower( (unsigned char) *s1++ );
    c2 = tolower( (unsigned char) *s2++ );
  } while (c1 == c2 && c1 != 0);
  return c2 > c1 ? -1 : c1 > c2;
}
#endif

#ifndef STRACPY
/* this is a replacement to strncpy, but ensures that the copy ends with NULL */
/* http://stracpy.blogspot.fr/2011/04/stracpy-strncpy-replacement.html */
#define STRACPY
char *stracpy(char *destination, const char *source, size_t amount)
{
        if (!destination || !source || !amount) return(NULL);
        while(amount--)
          if((*destination++ = *source++) == '\0') break;
        *destination = '\0';
        return destination;
}
#endif

/*******************************************************************************
* mcfull_file: allocates a full file name=dirname+file. Catenate extension if missing.
*******************************************************************************/
char *mcfull_file(char *name, char *ext)
{
  int   dirlen=0;
  char *mem   =NULL;

  dirlen = dirname ? strlen(dirname) : 0;
  mem = (char*)malloc(dirlen + strlen(name) + CHAR_BUF_LENGTH);
  if(!mem) {
    exit(-fprintf(stderr, "Error: Out of memory %li (mcfull_file)\n", (long)(dirlen + strlen(name) + 256)));
  }
  strcpy(mem, "");

  /* prepend directory name to path if name does not contain a path */
  if (dirlen > 0 && !strchr(name, MC_PATHSEP_C)) {
    strcat(mem, dirname);
    strcat(mem, MC_PATHSEP_S);
  } /* dirlen */

  strcat(mem, name);
  if (!strchr(name, '.') && ext && strlen(ext))
  { /* add extension if not in file name already */
    strcat(mem, ".");
    strcat(mem, ext);
  }
  return(mem);
} /* mcfull_file */

/*******************************************************************************
* mcnew_file: opens a new file within dirname if non NULL
*             the file is opened in "a" (append, create if does not exist)
*             the extension 'ext' is added if the file name does not include one.
*             the last argument is set to 0 if file did not exist, else to 1.
*******************************************************************************/
FILE *mcnew_file(char *name, char *ext, int *exists)
{
  char *mem;
  FILE *file=NULL;

  if (!name || strlen(name) == 0 || mcdisable_output_files) return(NULL);

  mem  = mcfull_file(name, ext); /* create dirname/name.ext */

  /* check for existence */
  file = fopen(mem, "r"); /* for reading -> fails if does not exist */
  if (file) {
    fclose(file);
    *exists=1;
  } else
    *exists=0;

  /* open the file for writing/appending */
#ifdef USE_NEXUS
  if (mcformat && strcasestr(mcformat, "NeXus")) {
    /* NXhandle nxhandle is defined in the .h with USE_NEXUS */
    NXaccess mode = (*exists ? NXACC_CREATE5 | NXACC_RDWR : NXACC_CREATE5);

    if (NXopen(mem, mode, &nxhandle) != NX_OK)
      file = NULL;
    else
      file = (FILE*)&nxhandle; /* to make it non NULL */
  } else
#endif
    file = fopen(mem, "a+");

  if(!file)
    fprintf(stderr, "Warning: could not open output file '%s' for %s (mcnew_file)\n",
      mem, *exists ? "append" : "create");
  free(mem);

  return file;
} /* mcnew_file */

/*******************************************************************************
* mcdetector_statistics: compute detector statistics, error bars, [x I I_err N] 1D
* RETURN:            updated detector structure
* Used by: detector_import
*******************************************************************************/
MCDETECTOR mcdetector_statistics(
  MCDETECTOR detector)
{

  if (!detector.p1 || !detector.m || !detector.filename)
    return(detector);

  /* compute statistics and update MCDETECTOR structure ===================== */
  double sum_z  = 0, min_z  = 0, max_z  = 0;
  double fmon_x =0,  smon_x = 0, fmon_y =0, smon_y=0, mean_z=0;
  double Nsum=0, P2sum=0;

  double sum_xz = 0, sum_yz = 0, sum_x = 0, sum_y = 0, sum_x2z = 0, sum_y2z = 0;
  int    i,j;
  char   hasnan=0, hasinf=0;
  char   israw = ((char*)strcasestr(detector.format,"raw") != NULL);
  double *this_p1=NULL; /* new 1D McCode array [x I E N]. Freed after writing data */

  /* if McCode/PGPLOT and rank==1 we create a new m*4 data block=[x I E N] */
  if (detector.rank == 1 && strcasestr(detector.format,"McCode")) {
    this_p1 = (double *)calloc(detector.m*detector.n*detector.p*4, sizeof(double));
    if (!this_p1)
      exit(-fprintf(stderr, "Error: Out of memory creating %li 1D " MCCODE_STRING " data set for file '%s' (detector_import)\n",
        detector.m*detector.n*detector.p*4*sizeof(double*), detector.filename));
  }

  max_z = min_z = detector.p1[0];

  /* compute sum and moments (not for lists) */
  if (!strcasestr(detector.format,"list") && detector.m)
  for(j = 0; j < detector.n*detector.p; j++)
  {
    for(i = 0; i < detector.m; i++)
    {
      double x,y,z;
      double N, E;
      long   index= !detector.istransposed ? i*detector.n*detector.p + j : i+j*detector.m;
      char   hasnaninf=0;

      if (detector.m)
        x = detector.xmin + (i + 0.5)/detector.m*(detector.xmax - detector.xmin);
      else x = 0;
      if (detector.n && detector.p)
        y = detector.ymin + (j + 0.5)/detector.n/detector.p*(detector.ymax - detector.ymin);
      else y = 0;
      z = detector.p1[index];
      N = detector.p0 ? detector.p0[index] : 1;
      E = detector.p2 ? detector.p2[index] : 0;
      if (detector.p2 && !israw)
        detector.p2[index] = (*mcestimate_error_p)(detector.p0[index],detector.p1[index],detector.p2[index]); /* set sigma */

      if (detector.rank == 1 && this_p1 && strcasestr(detector.format,"McCode")) {
        /* fill-in 1D McCode array [x I E N] */
        this_p1[index*4]   = x;
        this_p1[index*4+1] = z;
        this_p1[index*4+2] = detector.p2 ? detector.p2[index] : 0;
        this_p1[index*4+3] = N;
      }

      if (isnan(z) || isnan(E) || isnan(N)) hasnaninf=hasnan=1;
      if (isinf(z) || isinf(E) || isinf(N)) hasnaninf=hasinf=1;

      /* compute stats integrals */
      if (!hasnaninf) {
        sum_xz += x*z;
        sum_yz += y*z;
        sum_x  += x;
        sum_y  += y;
        sum_z  += z;
        sum_x2z += x*x*z;
        sum_y2z += y*y*z;
        if (z > max_z) max_z = z;
        if (z < min_z) min_z = z;

        Nsum += N;
        P2sum += E;
      }

    }
  } /* for j */

  /* compute 1st and 2nd moments. For lists, sum_z=0 so this is skipped. */
  if (sum_z && detector.n*detector.m*detector.p)
  {
    fmon_x = sum_xz/sum_z;
    fmon_y = sum_yz/sum_z;
    smon_x = sum_x2z/sum_z-fmon_x*fmon_x; smon_x = smon_x > 0 ? sqrt(smon_x) : 0;
    smon_y = sum_y2z/sum_z-fmon_y*fmon_y; smon_y = smon_y > 0 ? sqrt(smon_y) : 0;
    mean_z = sum_z/detector.n/detector.m/detector.p;
  }
  /* store statistics into detector */
  detector.intensity = sum_z;
  detector.error     = Nsum ? (*mcestimate_error_p)(Nsum, sum_z, P2sum) : 0;
  detector.events    = Nsum;
  detector.min       = min_z;
  detector.max       = max_z;
  detector.mean      = mean_z;
  detector.centerX   = fmon_x;
  detector.halfwidthX= smon_x;
  detector.centerY   = fmon_y;
  detector.halfwidthY= smon_y;

  /* if McCode/PGPLOT and rank==1 replace p1 with new m*4 1D McCode and clear others */
  if (detector.rank == 1 && this_p1 && strcasestr(detector.format,"McCode")) {

    detector.p1 = this_p1;
    detector.n  = detector.m; detector.m  = 4;
    detector.p0 = detector.p2 = NULL;
    detector.istransposed = 1;
  }

  if (detector.n*detector.m*detector.p > 1)
    snprintf(detector.signal, CHAR_BUF_LENGTH,
      "Min=%g; Max=%g; Mean=%g;", detector.min, detector.max, detector.mean);
  else
    strcpy(detector.signal, "None");
  snprintf(detector.values, CHAR_BUF_LENGTH,
    "%g %g %g", detector.intensity, detector.error, detector.events);

  switch (detector.rank) {
    case 1:  snprintf(detector.statistics, CHAR_BUF_LENGTH, "X0=%g; dX=%g;",
      detector.centerX, detector.halfwidthX); break;
    case 2:
    case 3:  snprintf(detector.statistics, CHAR_BUF_LENGTH, "X0=%g; dX=%g; Y0=%g; dY=%g;",
      detector.centerX, detector.halfwidthX, detector.centerY, detector.halfwidthY);
      break;
    default: strcpy(detector.statistics, "None");
  }

  if (hasnan)
    printf("WARNING: Nan detected in component/file %s %s\n",
      detector.component, strlen(detector.filename) ? detector.filename : "");
  if (hasinf)
    printf("WARNING: Inf detected in component/file %s %s\n",
      detector.component, strlen(detector.filename) ? detector.filename : "");

  return(detector);

} /* mcdetector_statistics */

/*******************************************************************************
* detector_import: build detector structure, merge non-lists from MPI
*                    compute basic stat, write "Detector:" line
* RETURN:            detector structure. Invalid data if detector.p1 == NULL
*                    Invalid detector sets m=0 and filename=""
*                    Simulation data  sets m=0 and filename=siminfo_name
* This function is equivalent to the old 'mcdetector_out', returning a structure
*******************************************************************************/
MCDETECTOR detector_import(
  char *format,
  char *component, char *title,
  long m, long n,  long p,
  char *xlabel, char *ylabel, char *zlabel,
  char *xvar, char *yvar, char *zvar,
  double x1, double x2, double y1, double y2, double z1, double z2,
  char *filename,
  double *p0, double *p1, double *p2,
  Coords position)
{
  time_t t;       /* for detector.date */
  long   date_l;  /* date as a long number */
  char   istransposed=0;
  char   c[CHAR_BUF_LENGTH]; /* temp var for signal label */

  MCDETECTOR detector;

  /* build MCDETECTOR structure ============================================= */
  /* make sure we do not have NULL for char fields */

  /* these also apply to simfile */
  strncpy (detector.filename,  filename ? filename : "",        CHAR_BUF_LENGTH);
  strncpy (detector.format,    format   ? format   : "McCode" , CHAR_BUF_LENGTH);
  /* add extension if missing */
  if (strlen(detector.filename) && !strchr(detector.filename, '.'))
  { /* add extension if not in file name already */
    strcat(detector.filename, ".dat");
  }
  strncpy (detector.component, component ? component : MCCODE_STRING " component", CHAR_BUF_LENGTH);

  snprintf(detector.instrument, CHAR_BUF_LENGTH, "%s (%s)", instrument_name, instrument_source);
  snprintf(detector.user, CHAR_BUF_LENGTH,      "%s on %s",
        getenv("USER") ? getenv("USER") : MCCODE_NAME,
        getenv("HOST") ? getenv("HOST") : "localhost");
  time(&t);         /* get current write time */
  date_l = (long)t; /* same but as a long */
  snprintf(detector.date, CHAR_BUF_LENGTH, "%s", ctime(&t));
  if (strlen(detector.date))   detector.date[strlen(detector.date)-1] = '\0'; /* remove last \n in date */
  detector.date_l = date_l;

  if (!mcget_run_num() || mcget_run_num() >= mcget_ncount())
    snprintf(detector.ncount, CHAR_BUF_LENGTH, "%llu", mcget_ncount()
#ifdef USE_MPI
*mpi_node_count
#endif
  );
  else
    snprintf(detector.ncount, CHAR_BUF_LENGTH, "%g/%g", (double)mcget_run_num(), (double)mcget_ncount());

  detector.p0         = p0;
  detector.p1         = p1;
  detector.p2         = p2;

  /* handle transposition (not for NeXus) */
  if (!strcasestr(detector.format, "NeXus")) {
    if (m<0 || n<0 || p<0)             istransposed = !istransposed;
    if (strcasestr(detector.format, "transpose")) istransposed = !istransposed;
    if (istransposed) { /* do the swap once for all */
      long i=m; m=n; n=i;
    }
  }

  m=abs(m); n=abs(n); p=abs(p); /* make sure dimensions are positive */
  detector.istransposed = istransposed;

  /* determine detector rank (dimensionality) */
  if (!m || !n || !p || !p1) detector.rank = 4; /* invalid: exit with m=0 filename="" */
  else if (m*n*p == 1)       detector.rank = 0; /* 0D */
  else if (n == 1 || m == 1) detector.rank = 1; /* 1D */
  else if (p == 1)           detector.rank = 2; /* 2D */
  else                       detector.rank = 3; /* 3D */

  /* from rank, set type */
  switch (detector.rank) {
    case 0:  strcpy(detector.type,  "array_0d"); m=n=p=1; break;
    case 1:  snprintf(detector.type, CHAR_BUF_LENGTH, "array_1d(%ld)", m*n*p); m *= n*p; n=p=1; break;
    case 2:  snprintf(detector.type, CHAR_BUF_LENGTH, "array_2d(%ld, %ld)", m, n*p); n *= p; p=1; break;
    case 3:  snprintf(detector.type, CHAR_BUF_LENGTH, "array_3d(%ld, %ld, %ld)", m, n, p); break;
    default: m=0; strcpy(detector.type, ""); strcpy(detector.filename, "");/* invalid */
  }

  detector.m    = m;
  detector.n    = n;
  detector.p    = p;

  /* these only apply to detector files ===================================== */

  snprintf(detector.position, CHAR_BUF_LENGTH, "%g %g %g", position.x, position.y, position.z);
  /* may also store actual detector orientation in the future */

  strncpy(detector.title,      title && strlen(title) ? title : component,       CHAR_BUF_LENGTH);
  strncpy(detector.xlabel,     xlabel && strlen(xlabel) ? xlabel : "X", CHAR_BUF_LENGTH); /* axis labels */
  strncpy(detector.ylabel,     ylabel && strlen(ylabel) ? ylabel : "Y", CHAR_BUF_LENGTH);
  strncpy(detector.zlabel,     zlabel && strlen(zlabel) ? zlabel : "Z", CHAR_BUF_LENGTH);
  strncpy(detector.xvar,       xvar && strlen(xvar) ? xvar :       "x", CHAR_BUF_LENGTH); /* axis variables */
  strncpy(detector.yvar,       yvar && strlen(yvar) ? yvar :       detector.xvar, CHAR_BUF_LENGTH);
  strncpy(detector.zvar,       zvar && strlen(zvar) ? zvar :       detector.yvar, CHAR_BUF_LENGTH);

  /* set "variables" as e.g. "I I_err N" */
  strcpy(c, "I ");
  if (strlen(detector.zvar))      strncpy(c, detector.zvar,32);
  else if (strlen(detector.yvar)) strncpy(c, detector.yvar,32);
  else if (strlen(detector.xvar)) strncpy(c, detector.xvar,32);

  if (detector.rank == 1)
    snprintf(detector.variables, CHAR_BUF_LENGTH, "%s %s %s_err N", detector.xvar, c, c);
  else
    snprintf(detector.variables, CHAR_BUF_LENGTH, "%s %s_err N", c, c);

  /* limits */
  detector.xmin = x1;
  detector.xmax = x2;
  detector.ymin = y1;
  detector.ymax = y2;
  detector.zmin = z1;
  detector.zmax = z2;
  if (abs(detector.rank) == 1)
    snprintf(detector.limits, CHAR_BUF_LENGTH, "%g %g", x1, x2);
  else if (detector.rank == 2)
    snprintf(detector.limits, CHAR_BUF_LENGTH, "%g %g %g %g", x1, x2, y1, y2);
  else
    snprintf(detector.limits, CHAR_BUF_LENGTH, "%g %g %g %g %g %g", x1, x2, y1, y2, z1, z2);

  /* if MPI and nodes_nb > 1: reduce data sets when using MPI =============== */
#ifdef USE_MPI
  if (!strcasestr(detector.format,"list") && mpi_node_count > 1 && m) {
    /* we save additive data: reduce everything into mpi_node_root */
    if (p0) mc_MPI_Sum(p0, m*n*p);
    if (p1) mc_MPI_Sum(p1, m*n*p);
    if (p2) mc_MPI_Sum(p2, m*n*p);
    if (!p0) {  /* additive signal must be then divided by the number of nodes */
      int i;
      for (i=0; i<m*n*p; i++) {
        p1[i] /= mpi_node_count;
        if (p2) p2[i] /= mpi_node_count;
      }
    }
  }
#endif /* USE_MPI */

  /* compute statistics, Nsum, intensity, Error bars */
  detector = mcdetector_statistics(detector);

#ifdef USE_MPI
  /* slaves are done */
  if(mpi_node_rank != mpi_node_root) {
    return detector;
  }
#endif

  /* output "Detector:" line ================================================ */
  /* when this is a detector written by a component (not the SAVE from instrument),
     not an event lists */
  if (!m) return(detector);
  if (!strcasestr(detector.format,"list")) {
    if (!strcmp(detector.component, instrument_name)) {
      if (strlen(detector.filename))  /* we name it from its filename, or from its title */
        strncpy(c, detector.filename, CHAR_BUF_LENGTH);
      else
        snprintf(c, CHAR_BUF_LENGTH, "%s", instrument_name);
    } else
      strncpy(c, detector.component, CHAR_BUF_LENGTH);  /* usual detectors written by components */

    printf("Detector: %s_I=%g %s_ERR=%g %s_N=%g",
           c, detector.intensity,
           c, detector.error,
           c, detector.events);
    printf(" \"%s\"\n", strlen(detector.filename) ? detector.filename : detector.component);
  }


  return(detector);
} /* detector_import */

/* end MCDETECTOR import section ============================================ */

















/* ========================================================================== */

/*                               ASCII output                                 */
/*     The SIM file is YAML based, the data files have '#' headers            */

/* ========================================================================== */


/*******************************************************************************
* mcinfo_out: output instrument tags/info (only in SIM)
* Used in: siminfo_init (ascii), mcinfo(stdout)
*******************************************************************************/
static void mcinfo_out(char *pre, FILE *f)
{
  char Parameters[CHAR_BUF_LENGTH] = "";
  int  i;

  if (!f || mcdisable_output_files) return;

  /* create parameter string ================================================ */
  for(i = 0; i < numipar; i++)
  {
    char ThisParam[CHAR_BUF_LENGTH];
    if (strlen(mcinputtable[i].name) > CHAR_BUF_LENGTH) break;
    snprintf(ThisParam, CHAR_BUF_LENGTH, " %s(%s)", mcinputtable[i].name,
            (*mcinputtypes[mcinputtable[i].type].parminfo)
                (mcinputtable[i].name));
    strcat(Parameters, ThisParam);
    if (strlen(Parameters) >= CHAR_BUF_LENGTH-64) break;
  }

  /* output data ============================================================ */
  if (f != stdout)
    fprintf(f, "%sFile: %s%c%s\n",    pre, dirname, MC_PATHSEP_C, siminfo_name);
  else
    fprintf(f, "%sCreator: %s\n",     pre, MCCODE_STRING);

  fprintf(f, "%sSource: %s\n",   pre, instrument_source);
  fprintf(f, "%sParameters: %s\n",    pre, Parameters);

  fprintf(f, "%sTrace_enabled: %s\n", pre, traceenabled ? "yes" : "no");
  fprintf(f, "%sDefault_main: %s\n",  pre, defaultmain ?  "yes" : "no");
  fprintf(f, "%sEmbedded_runtime: %s\n", pre,
#ifdef MC_EMBEDDED_RUNTIME
         "yes"
#else
         "no"
#endif
         );

  fflush(f);
} /* mcinfo_out */

/*******************************************************************************
* mcruninfo_out: output simulation tags/info (both in SIM and data files)
* Used in: siminfo_init (ascii case), mcdetector_out_xD_ascii
*******************************************************************************/
static void mcruninfo_out(char *pre, FILE *f)
{
  int i;
  char Parameters[CHAR_BUF_LENGTH];

  if (!f || mcdisable_output_files) return;

  fprintf(f, "%sFormat: %s%s\n",      pre,
    mcformat && strlen(mcformat) ? mcformat : MCCODE_NAME,
    mcformat && strcasestr(mcformat,"McCode") ? " with text headers" : "");
  fprintf(f, "%sURL: %s\n",         pre, "http://www.mccode.org");
  fprintf(f, "%sCreator: %s\n",     pre, MCCODE_STRING);
  fprintf(f, "%sInstrument: %s\n", pre, instrument_source);
  fprintf(f, "%sNcount: %llu\n",        pre, mcget_ncount());
  fprintf(f, "%sTrace: %s\n",       pre, mcdotrace ? "yes" : "no");
  fprintf(f, "%sGravitation: %s\n", pre, mcgravitation ? "yes" : "no");
  snprintf(Parameters, CHAR_BUF_LENGTH, "%ld", mcseed);
  fprintf(f, "%sSeed: %s\n",        pre, Parameters);
  fprintf(f, "%sDirectory: %s\n",        pre, dirname ? dirname : ".");
#ifdef USE_MPI
  if (mpi_node_count > 1)
    fprintf(f, "%sNodes: %i\n",        pre, mpi_node_count);
#endif

  /* output parameter string ================================================ */
  for(i = 0; i < numipar; i++) {
      if (mcinputtable[i].par){
	/* Parameters with a default value */
	if(mcinputtable[i].val && strlen(mcinputtable[i].val)){
	  (*mcinputtypes[mcinputtable[i].type].printer)(Parameters, mcinputtable[i].par);
	  fprintf(f, "%sParam: %s=%s\n", pre, mcinputtable[i].name, Parameters);
        /* ... and those without */
	}else{
	  fprintf(f, "%sParam: %s=NULL\n", pre, mcinputtable[i].name);
	}
      }
  }
  fflush(f);
} /* mcruninfo_out */

/*******************************************************************************
* siminfo_out:    wrapper to fprintf(siminfo_file)
*******************************************************************************/
void siminfo_out(char *format, ...)
{
  va_list ap;

  if(siminfo_file && !mcdisable_output_files)
  {
    va_start(ap, format);
    vfprintf(siminfo_file, format, ap);
    va_end(ap);
  }
} /* siminfo_out */


/*******************************************************************************
* mcdatainfo_out: output detector header
*   mcdatainfo_out(prefix, file_handle, detector) writes info to data file
*******************************************************************************/
static void
mcdatainfo_out(char *pre, FILE *f, MCDETECTOR detector)
{
  if (!f || !detector.m || mcdisable_output_files) return;

  /* output data ============================================================ */
  fprintf(f, "%sDate: %s (%li)\n",       pre, detector.date, detector.date_l);
  fprintf(f, "%stype: %s\n",       pre, detector.type);
  fprintf(f, "%sSource: %s\n",     pre, detector.instrument);
  fprintf(f, "%scomponent: %s\n",  pre, detector.component);
  fprintf(f, "%sposition: %s\n",   pre, detector.position);

  fprintf(f, "%stitle: %s\n",      pre, detector.title);
  fprintf(f, !mcget_run_num() || mcget_run_num() >= mcget_ncount() ?
             "%sNcount: %s\n" :
             "%sratio: %s\n",  pre, detector.ncount);

  if (strlen(detector.filename)) {
    fprintf(f, "%sfilename: %s\n", pre, detector.filename);
  }

  fprintf(f, "%sstatistics: %s\n", pre, detector.statistics);
  fprintf(f, "%ssignal: %s\n",     pre, detector.signal);
  fprintf(f, "%svalues: %s\n",     pre, detector.values);

  if (detector.rank >= 1)
  {
    fprintf(f, "%sxvar: %s\n",     pre, detector.xvar);
    fprintf(f, "%syvar: %s\n",     pre, detector.yvar);
    fprintf(f, "%sxlabel: %s\n",   pre, detector.xlabel);
    fprintf(f, "%sylabel: %s\n",   pre, detector.ylabel);
    if (detector.rank > 1) {
      fprintf(f, "%szvar: %s\n",   pre, detector.zvar);
      fprintf(f, "%szlabel: %s\n", pre, detector.zlabel);
    }
  }

  fprintf(f,
    abs(detector.rank)==1 ?
             "%sxlimits: %s\n" :
             "%sxylimits: %s\n", pre, detector.limits);
  fprintf(f, "%svariables: %s\n", pre,
    strcasestr(detector.format, "list") ? detector.ylabel : detector.variables);

  fflush(f);

} /* mcdatainfo_out */

/* mcdetector_out_array_ascii: output a single array to a file
 *   m: columns
 *   n: rows
 *   p: array
 *   f: file handle (already opened)
 */
static void mcdetector_out_array_ascii(long m, long n, double *p, FILE *f, char istransposed)
{
  if(f)
  {
    int i,j;
    for(j = 0; j < n; j++)
    {
      for(i = 0; i < m; i++)
      {
          fprintf(f, "%.10g ", p[!istransposed ? i*n + j : j*m+i]);
      }
      fprintf(f,"\n");
    }
  }
} /* mcdetector_out_array_ascii */

/*******************************************************************************
* mcdetector_out_0D_ascii: called by mcdetector_out_0D for ascii output
*******************************************************************************/
MCDETECTOR mcdetector_out_0D_ascii(MCDETECTOR detector)
{
  int exists=0;
  FILE *outfile = NULL;

  /* Write data set information to simulation description file. */
  MPI_MASTER(
    siminfo_out("\nbegin data\n"); // detector.component
    mcdatainfo_out("  ", siminfo_file, detector);
    siminfo_out("end data\n");
    /* Don't write if filename is NULL: mcnew_file handles this (return NULL) */
    outfile = mcnew_file(detector.component, "dat", &exists);
    if(outfile)
    {
      /* write data file header and entry in simulation description file */
      mcruninfo_out( "# ", outfile);
      mcdatainfo_out("# ", outfile, detector);
      /* write I I_err N */
      fprintf(outfile, "%g %g %g\n",
        detector.intensity, detector.error, detector.events);
      fclose(outfile);
    }
  ); /* MPI_MASTER */
  return(detector);
} /* mcdetector_out_0D_ascii */

/*******************************************************************************
* mcdetector_out_1D_ascii: called by mcdetector_out_1D for ascii output
*******************************************************************************/
MCDETECTOR mcdetector_out_1D_ascii(MCDETECTOR detector)
{
  int exists=0;
  FILE *outfile = NULL;

  MPI_MASTER(
    /* Write data set information to simulation description file. */
    siminfo_out("\nbegin data\n"); // detector.filename
    mcdatainfo_out("  ", siminfo_file, detector);
    siminfo_out("end data\n");
    /* Loop over array elements, writing to file. */
    /* Don't write if filename is NULL: mcnew_file handles this (return NULL) */
    outfile = mcnew_file(detector.filename, "dat", &exists);
    if(outfile)
    {
      /* write data file header and entry in simulation description file */
      mcruninfo_out( "# ", outfile);
      mcdatainfo_out("# ", outfile, detector);
      /* output the 1D array columns */
      mcdetector_out_array_ascii(detector.m, detector.n, detector.p1, outfile, detector.istransposed);

      fclose(outfile);
    }
  ); /* MPI_MASTER */
  return(detector);

}  /* mcdetector_out_1D_ascii */

/*******************************************************************************
* mcdetector_out_2D_ascii: called by mcdetector_out_2D for ascii output
*******************************************************************************/
MCDETECTOR mcdetector_out_2D_ascii(MCDETECTOR detector)
{
  int exists=0;
  FILE *outfile = NULL;

  MPI_MASTER(
    /* Loop over array elements, writing to file. */
    /* Don't write if filename is NULL: mcnew_file handles this (return NULL) */
    outfile = mcnew_file(detector.filename, "dat", &exists);
    if(outfile)
    {
      /* write header only if file has just been created (not appending) */
      if (!exists) {
        /* Write data set information to simulation description file. */
        siminfo_out("\nbegin data\n"); // detector.filename
        mcdatainfo_out("  ", siminfo_file, detector);
        siminfo_out("end data\n");

        mcruninfo_out( "# ", outfile);
        mcdatainfo_out("# ", outfile,   detector);
        fprintf(outfile, "# Data [%s/%s] %s:\n", detector.component, detector.filename, detector.zvar);
      }
      mcdetector_out_array_ascii(detector.m, detector.n*detector.p, detector.p1,
        outfile, detector.istransposed);
      if (detector.p2) {
        fprintf(outfile, "# Errors [%s/%s] %s_err:\n", detector.component, detector.filename, detector.zvar);
        mcdetector_out_array_ascii(detector.m, detector.n*detector.p, detector.p2,
          outfile, detector.istransposed);
      }
      if (detector.p0) {
        fprintf(outfile, "# Events [%s/%s] N:\n", detector.component, detector.filename);
        mcdetector_out_array_ascii(detector.m, detector.n*detector.p, detector.p0,
          outfile, detector.istransposed);
      }
      fclose(outfile);

      if (!exists) {
        if (strcasestr(detector.format, "list"))
          printf("Events:   \"%s\"\n",
            strlen(detector.filename) ? detector.filename : detector.component);
      }
    } /* if outfile */
  ); /* MPI_MASTER */
#ifdef USE_MPI
  if (strcasestr(detector.format, "list") && mpi_node_count > 1) {
    int node_i=0;
    /* loop along MPI nodes to write sequentially */
    for(node_i=0; node_i<mpi_node_count; node_i++) {
      /* MPI: slaves wait for the master to write its block, then append theirs */
      MPI_Barrier(MPI_COMM_WORLD);
      if (node_i != mpi_node_root && node_i == mpi_node_rank) {
        if(strlen(detector.filename) && !mcdisable_output_files)	/* Don't write if filename is NULL */
          outfile = mcnew_file(detector.filename, "dat", &exists);
        if (!exists)
          fprintf(stderr, "Warning: [MPI node %i] file '%s' does not exist yet, "
                          "MASTER should have opened it before.\n",
            mpi_node_rank, detector.filename);
        if(outfile) {
          mcdetector_out_array_ascii(detector.m, detector.n*detector.p, detector.p1,
            outfile, detector.istransposed);
          fclose(outfile);
        }
      }
    }
  } /* if strcasestr list */
#endif
  return(detector);
} /* mcdetector_out_2D_ascii */

/*******************************************************************************
* strcpy_valid: makes a valid string for variable names.
*   copy 'original' into 'valid', replacing invalid characters by '_'
*   char arrays must be pre-allocated
*******************************************************************************/
static char *strcpy_valid(char *valid, char *original)
{
  long i;
  int  n=32; /* max length of valid names */

  if (original == NULL || !strlen(original)) return(NULL);

  if (n > strlen(original)) n = strlen(original);
  else original += strlen(original)-n;
  strncpy(valid, original, n);

  for (i=0; i < n; i++)
  {
    if ( (valid[i] > 122)
      || (valid[i] < 32)
      || (strchr("!\"#$%&'()*+,-.:;<=>?@[\\]^`/ \n\r\t", valid[i]) != NULL) )
    {
      if (i) valid[i] = '_'; else valid[i] = 'm';
    }
  }
  valid[i] = '\0';

  return(valid);
} /* strcpy_valid */

/* end ascii output section ================================================= */







#ifdef USE_NEXUS

/* ========================================================================== */

/*                               NeXus output                                 */

/* ========================================================================== */

#define nxprintf(...)    nxstr('d', __VA_ARGS__)
#define nxprintattr(...) nxstr('a', __VA_ARGS__)

/*******************************************************************************
* nxstr: output a tag=value data set (char) in NeXus/current group
*   when 'format' is larger that 1024 chars it is used as value for the 'tag'
*   else the value is assembled with format and following arguments.
*   type='d' -> data set
*        'a' -> attribute for current data set
*******************************************************************************/
static int nxstr(char type, NXhandle *f, char *tag, char *format, ...)
{
  va_list ap;
  char value[CHAR_BUF_LENGTH];
  int  i;
  int  ret=NX_OK;

  if (!tag || !format || !strlen(tag) || !strlen(format)) return(NX_OK);

  /* assemble the value string */
  if (strlen(format) < CHAR_BUF_LENGTH) {
    va_start(ap, format);
    ret = vsnprintf(value, CHAR_BUF_LENGTH, format, ap);
    va_end(ap);

    i = strlen(value);
  } else {
    i = strlen(format);
  }

  if (type == 'd') {
    /* open/put/close data set */
    if (NXmakedata (f, tag, NX_CHAR, 1, &i) != NX_OK) return(NX_ERROR);
    NXopendata (f, tag);
    if (strlen(format) < CHAR_BUF_LENGTH)
      ret = NXputdata  (f, value);
    else
      ret = NXputdata  (f, format);
    NXclosedata(f);
  } else {
    if (strlen(format) < CHAR_BUF_LENGTH)
      ret = NXputattr  (f, tag, value, strlen(value), NX_CHAR);
    else
      ret = NXputattr  (f, tag, format, strlen(format), NX_CHAR);
  }

  return(ret);

} /* nxstr */

/*******************************************************************************
* mcinfo_readfile: read a full file into a string buffer which is allocated
*   Think to free the buffer after use.
* Used in: mcinfo_out_nexus (nexus)
*******************************************************************************/
char *mcinfo_readfile(char *filename)
{
  FILE *f = fopen(filename, "rb");
  if (!f) return(NULL);
  fseek(f, 0, SEEK_END);
  long fsize = ftell(f);
  rewind(f);
  char *string = malloc(fsize + 1);
  if (string) {
    int n = fread(string, fsize, 1, f);
    fclose(f);

    string[fsize] = 0;
  }
  return(string);
}

/*******************************************************************************
* mcinfo_out: output instrument/simulation groups in NeXus file
* Used in: siminfo_init (nexus)
*******************************************************************************/
static void mcinfo_out_nexus(NXhandle f)
{
  FILE  *fid;     /* for intrument source code/C/IDF */
  char  *buffer=NULL;
  time_t t     =time(NULL); /* for date */
  char   entry0[CHAR_BUF_LENGTH];
  int    count=0;
  char   name[CHAR_BUF_LENGTH];
  char   class[CHAR_BUF_LENGTH];

  if (!f || mcdisable_output_files) return;

  /* write NeXus NXroot attributes */
  /* automatically added: file_name, HDF5_Version, file_time, NeXus_version */
  nxprintattr(f, "creator",   "%s generated with " MCCODE_STRING, instrument_name);

  /* count the number of existing NXentry and create the next one */
  NXgetgroupinfo(f, &count, name, class);
  sprintf(entry0, "entry%i", count+1);

  /* create the main NXentry (mandatory in NeXus) */
  if (NXmakegroup(f, entry0, "NXentry") == NX_OK)
  if (NXopengroup(f, entry0, "NXentry") == NX_OK) {

    nxprintf(nxhandle, "program_name", MCCODE_STRING);
    nxprintf(f, "start_time", ctime(&t));
    nxprintf(f, "title", "%s%s%s simulation generated by instrument %s",
      dirname && strlen(dirname) ? dirname : ".", MC_PATHSEP_S, siminfo_name,
      instrument_name);
    nxprintattr(f, "program_name", MCCODE_STRING);
    nxprintattr(f, "instrument",   instrument_name);
    nxprintattr(f, "simulation",   "%s%s%s",
        dirname && strlen(dirname) ? dirname : ".", MC_PATHSEP_S, siminfo_name);

    /* write NeXus instrument group */
    if (NXmakegroup(f, "instrument", "NXinstrument") == NX_OK)
    if (NXopengroup(f, "instrument", "NXinstrument") == NX_OK) {
      int   i;
      char *string=NULL;

      /* write NeXus parameters(types) data =================================== */
      string = (char*)malloc(CHAR_BUF_LENGTH);
      if (string) {
        strcpy(string, "");
        for(i = 0; i < numipar; i++)
        {
          char ThisParam[CHAR_BUF_LENGTH];
          snprintf(ThisParam, CHAR_BUF_LENGTH, " %s(%s)", mcinputtable[i].name,
                  (*mcinputtypes[mcinputtable[i].type].parminfo)
                      (mcinputtable[i].name));
          if (strlen(string) + strlen(ThisParam) < CHAR_BUF_LENGTH)
            strcat(string, ThisParam);
        }
        nxprintattr(f, "Parameters",    string);
        free(string);
      }

      nxprintattr(f, "name",          instrument_name);
      nxprintf   (f, "name",          instrument_name);
      nxprintattr(f, "Source",        instrument_source);

      nxprintattr(f, "Trace_enabled", traceenabled ? "yes" : "no");
      nxprintattr(f, "Default_main",  defaultmain ?  "yes" : "no");
      nxprintattr(f, "Embedded_runtime",
  #ifdef MC_EMBEDDED_RUNTIME
           "yes"
  #else
           "no"
  #endif
           );

      /* add instrument source code when available */
      buffer = mcinfo_readfile(instrument_source);
      if (buffer && strlen(buffer)) {
        long length=strlen(buffer);
        nxprintf (f, "description", buffer);
        NXopendata(f,"description");
        nxprintattr(f, "file_name", instrument_source);
        nxprintattr(f, "file_size", "%li", length);
        nxprintattr(f, "MCCODE_STRING", MCCODE_STRING);
        NXclosedata(f);
        nxprintf (f,"instrument_source", "%s " MCCODE_NAME " " MCCODE_PARTICLE " Monte Carlo simulation", instrument_name);
        free(buffer);
      } else
        nxprintf (f, "description", "File %s not found (instrument description %s is missing)",
          instrument_source, instrument_name);

      /* add Mantid/IDF.xml when available */
      char *IDFfile=NULL;
      IDFfile = (char*)malloc(CHAR_BUF_LENGTH);
      sprintf(IDFfile,"%s%s",instrument_source,".xml");
      buffer = mcinfo_readfile(IDFfile);
      if (buffer && strlen(buffer)) {
        NXmakegroup (nxhandle, "instrument_xml", "NXnote");
        NXopengroup (nxhandle, "instrument_xml", "NXnote");
        nxprintf(f, "data", buffer);
        nxprintf(f, "description", "IDF.xml file found with instrument %s", instrument_source);
        nxprintf(f, "type", "text/xml");
        NXclosegroup(f); /* instrument_xml */
        free(buffer);
      }
      free(IDFfile);
      NXclosegroup(f); /* instrument */
    } /* NXinstrument */

    /* write NeXus simulation group */
    if (NXmakegroup(f, "simulation", "NXnote") == NX_OK)
    if (NXopengroup(f, "simulation", "NXnote") == NX_OK) {

      nxprintattr(f, "name",   "%s%s%s",
        dirname && strlen(dirname) ? dirname : ".", MC_PATHSEP_S, siminfo_name);

      nxprintf   (f, "name",      "%s",     siminfo_name);
      nxprintattr(f, "Format",    mcformat && strlen(mcformat) ? mcformat : MCCODE_NAME);
      nxprintattr(f, "URL",       "http://www.mccode.org");
      nxprintattr(f, "program",   MCCODE_STRING);
      nxprintattr(f, "Instrument",instrument_source);
      nxprintattr(f, "Trace",     mcdotrace ?     "yes" : "no");
      nxprintattr(f, "Gravitation",mcgravitation ? "yes" : "no");
      nxprintattr(f, "Seed",      "%li", mcseed);
      nxprintattr(f, "Directory", dirname);
    #ifdef USE_MPI
      if (mpi_node_count > 1)
        nxprintf(f, "Nodes", "%i",        mpi_node_count);
    #endif

      /* output parameter string ================================================ */
      if (NXmakegroup(f, "Param", "NXparameters") == NX_OK)
      if (NXopengroup(f, "Param", "NXparameters") == NX_OK) {
        int i;
        char string[CHAR_BUF_LENGTH];
        for(i = 0; i < numipar; i++) {
          if (mcget_run_num() || (mcinputtable[i].val && strlen(mcinputtable[i].val))) {
            if (mcinputtable[i].par == NULL)
              strncpy(string, (mcinputtable[i].val ? mcinputtable[i].val : ""), CHAR_BUF_LENGTH);
            else
              (*mcinputtypes[mcinputtable[i].type].printer)(string, mcinputtable[i].par);

            nxprintf(f,  mcinputtable[i].name, "%s", string);
            nxprintattr(f, mcinputtable[i].name, string);
          }
        }
        NXclosegroup(f); /* Param */
      } /* NXparameters */

      NXclosegroup(f); /* simulation */
    } /* NXsimulation */

    /* create a group to hold all monitors */
    NXmakegroup(f, "data", "NXdetector");

    /* leave the NXentry opened (closed at exit) */
  } /* NXentry */
} /* mcinfo_out_nexus */

/*******************************************************************************
* mcdatainfo_out_nexus: output detector header
*   mcdatainfo_out_nexus(detector) create group and write info to NeXus data file
*   open data:NXdetector then filename:NXdata and write headers/attributes
*   requires: NXentry to be opened
*******************************************************************************/
static void
mcdatainfo_out_nexus(NXhandle f, MCDETECTOR detector)
{
  char data_name[32];
  if (!f || !detector.m || mcdisable_output_files) return;

  strcpy_valid(data_name,
    detector.filename && strlen(detector.filename) ?
      detector.filename : detector.component);

  /* the NXdetector group has been created in mcinfo_out_nexus (siminfo_init) */
  if (NXopengroup(f, "data", "NXdetector") == NX_OK) {

    /* create and open the data group */
    /* this may fail when appending to list -> ignore/skip */
    NXMDisableErrorReporting(); /* unactivate NeXus error messages, as creation may fail */

    if (NXmakegroup(f, data_name, "NXdata") == NX_OK)
    if (NXopengroup(f, data_name, "NXdata") == NX_OK) {

      /* output metadata (as attributes) ======================================== */
      nxprintattr(f, "Date",       detector.date);
      nxprintattr(f, "type",       detector.type);
      nxprintattr(f, "Source",     detector.instrument);
      nxprintattr(f, "component",  detector.component);
      nxprintattr(f, "position",   detector.position);

      nxprintattr(f, "title",      detector.title);
      nxprintattr(f, !mcget_run_num() || mcget_run_num() >= mcget_ncount() ?
                 "Ncount" :
                 "ratio",  detector.ncount);

      if (strlen(detector.filename)) {
        nxprintattr(f, "filename", detector.filename);
      }

      nxprintattr(f, "statistics", detector.statistics);
      nxprintattr(f, "signal",     detector.signal);
      nxprintattr(f, "values",     detector.values);

      if (detector.rank >= 1)
      {
        nxprintattr(f, "xvar",     detector.xvar);
        nxprintattr(f, "yvar",     detector.yvar);
        nxprintattr(f, "xlabel",   detector.xlabel);
        nxprintattr(f, "ylabel",   detector.ylabel);
        if (detector.rank > 1) {
          nxprintattr(f, "zvar",   detector.zvar);
          nxprintattr(f, "zlabel", detector.zlabel);
        }
      }

      nxprintattr(f, abs(detector.rank)==1 ?
                 "xlimits" :
                 "xylimits", detector.limits);
      nxprintattr(f, "variables",
        strcasestr(detector.format, "list") ? detector.ylabel : detector.variables);
      nxprintf(f, "distance", detector.position);
      nxprintf(f, "acquisition_mode",
        strcasestr(detector.format, "list") ? "event" : "summed");

      NXclosegroup(f);
    } /* NXdata (filename) */
    NXMEnableErrorReporting();  /* re-enable NeXus error messages */
    NXclosegroup(f);
  } /* NXdetector (data) */

} /* mcdatainfo_out_nexus */

/*******************************************************************************
* mcdetector_out_axis_nexus: write detector axis into current NXdata
*   requires: NXdata to be opened
*******************************************************************************/
int mcdetector_out_axis_nexus(NXhandle f, char *label, char *var, int rank, long length, double min, double max)
{
  if (!f || length <= 1 || mcdisable_output_files || max == min) return(NX_OK);
  else {
    double axis[length];
    char valid[32];
    int dim=(int)length;
    int i;
    int nprimary=1;
    /* create an axis from [min:max] */
    for(i = 0; i < length; i++)
      axis[i] = min+(max-min)*(i+0.5)/length;
    /* create the data set */
    strcpy_valid(valid, label);
    NXcompmakedata(f, valid, NX_FLOAT64, 1, &dim, NX_COMP_LZW, &dim);
    /* open it */
    if (NXopendata(f, valid) != NX_OK) {
      fprintf(stderr, "Warning: could not open axis rank %i '%s' (NeXus)\n",
        rank, valid);
      return(NX_ERROR);
    }
    /* put the axis and its attributes */
    NXputdata  (f, axis);
    nxprintattr(f, "long_name",  label);
    nxprintattr(f, "short_name", var);
    NXputattr  (f, "axis",       &rank,     1, NX_INT32);
    nxprintattr(f, "units",      var);
    NXputattr  (f, "primary",    &nprimary, 1, NX_INT32);
    NXclosedata(f);

    return(NX_OK);
  }
} /* mcdetector_out_axis_nexus */

/*******************************************************************************
* mcdetector_out_array_nexus: write detector array into current NXdata (1D,2D)
*   requires: NXdata to be opened
*******************************************************************************/
int mcdetector_out_array_nexus(NXhandle f, char *part, double *data, MCDETECTOR detector)
{

  int dims[3]={detector.m,detector.n,detector.p};  /* number of elements to write */
  int signal=1;
  int exists=0;
  int current_dims[3]={0,0,0};
  int ret=NX_OK;

  if (!f || !data || !detector.m || mcdisable_output_files) return(NX_OK);

  /* when this is a list, we set 1st dimension to NX_UNLIMITED for creation */
  if (strcasestr(detector.format, "list")) dims[0] = NX_UNLIMITED;

  /* create the data set in NXdata group */
  NXMDisableErrorReporting(); /* unactivate NeXus error messages, as creation may fail */
  /* NXcompmakedata fails with NX_UNLIMITED */
  if (strcasestr(detector.format, "list"))
    ret = NXmakedata(    f, part, NX_FLOAT64, detector.rank, dims);
  else
    ret = NXcompmakedata(f, part, NX_FLOAT64, detector.rank, dims, NX_COMP_LZW, dims);
  if (ret != NX_OK) {
    /* failed: data set already exists */
    int datatype=0;
    int rank=0;
    exists=1;
    /* inquire current size of data set (nb of events stored) */
    NXopendata(f, part);
    NXgetinfo(f, &rank, current_dims, &datatype);
    NXclosedata(f);
  }
  NXMEnableErrorReporting();  /* re-enable NeXus error messages */
  dims[0] = detector.m; /* restore actual dimension from data writing */

  /* open the data set */
  if (NXopendata(f, part) == NX_ERROR) {
    fprintf(stderr, "Warning: could not open DataSet %s '%s' (NeXus)\n",
      part, detector.title);
    return(NX_ERROR);
  }
  if (strcasestr(detector.format, "list")) {
    current_dims[1] = current_dims[2] = 0; /* set starting location for writing slab */
    NXputslab(f, data, current_dims, dims);
    if (!exists)
      printf("Events:   \"%s\"\n",
        strlen(detector.filename) ? detector.filename : detector.component);
  } else {
    NXputdata (f, data);
  }

  if (strstr(part,"data") || strstr(part, "events")) {
    NXputattr(f, "signal", &signal, 1, NX_INT32);
    nxprintattr(f, "short_name", detector.filename && strlen(detector.filename) ?
      detector.filename : detector.component);
  }
  nxprintattr(f, "long_name", "%s '%s'", part, detector.title);
  NXclosedata(f);

  return(NX_OK);
} /* mcdetector_out_array_nexus */

/*******************************************************************************
* mcdetector_out_data_nexus: write detector axes+data into current NXdata
*   The data:NXdetector is opened, then filename:NXdata
*   requires: NXentry to be opened
*******************************************************************************/
int mcdetector_out_data_nexus(NXhandle f, MCDETECTOR detector)
{
  char data_name[32];

  if (!f || !detector.m || mcdisable_output_files) return(NX_OK);

  strcpy_valid(data_name,
    detector.filename && strlen(detector.filename) ?
      detector.filename : detector.component);

  /* the NXdetector group has been created in mcinfo_out_nexus (siminfo_init) */
  if (NXopengroup(f, "data", "NXdetector") == NX_OK) {

    /* the NXdata group has been created in mcdatainfo_out_nexus */
    if (NXopengroup(f, data_name, "NXdata") == NX_OK) {

      /* write axes, for histogram data sets, not for lists */
      if (!strcasestr(detector.format, "list")) {
        mcdetector_out_axis_nexus(f, detector.xlabel, detector.xvar,
          1, detector.m, detector.xmin, detector.xmax);

        mcdetector_out_axis_nexus(f, detector.ylabel, detector.yvar,
          2, detector.n, detector.ymin, detector.ymax);

        mcdetector_out_axis_nexus(f, detector.zlabel, detector.zvar,
          3, detector.p, detector.zmin, detector.zmax);

      } /* !list */

      /* write the actual data (appended if already exists) */
      if (!strcasestr(detector.format, "list")) {
        mcdetector_out_array_nexus(f, "data", detector.p1, detector);
        mcdetector_out_array_nexus(f, "errors", detector.p2, detector);
        mcdetector_out_array_nexus(f, "ncount", detector.p0, detector);
      } else
        mcdetector_out_array_nexus(  f, "events", detector.p1, detector);

      NXclosegroup(f);
    } /* NXdata */
    NXclosegroup(f);
  } /* NXdetector */

  return(NX_OK);
} /* mcdetector_out_array_nexus */

#ifdef USE_MPI
/*******************************************************************************
* mcdetector_out_list_slaves: slaves send their list data to master which writes
*   requires: NXentry to be opened
* WARNING: this method has a flaw: it requires all nodes to flush the lists
*   the same number of times. In case one node is just below the buffer size
*   when finishing (e.g. monitor_nd), it may not trigger save but others may.
*   Then the number of recv/send is not constant along nodes, and simulation stalls.
*******************************************************************************/
MCDETECTOR mcdetector_out_list_slaves(MCDETECTOR detector)
{
  int     node_i=0;
  MPI_MASTER(
	     printf("\n** MPI master gathering slave node list data ** \n");
  );

  if (mpi_node_rank != mpi_node_root) {
    /* MPI slave: slaves send their data to master: 2 MPI_Send calls */
    /* m, n, p must be sent first, since all slaves do not have the same number of events */
    int mnp[3]={detector.m,detector.n,detector.p};

    if (mc_MPI_Send(mnp, 3, MPI_INT, mpi_node_root)!= MPI_SUCCESS)
      fprintf(stderr, "Warning: proc %i to master: MPI_Send mnp list error (mcdetector_out_list_slaves)\n", mpi_node_rank);
    if (!detector.p1
     || mc_MPI_Send(detector.p1, mnp[0]*mnp[1]*mnp[2], MPI_DOUBLE, mpi_node_root) != MPI_SUCCESS)
      fprintf(stderr, "Warning: proc %i to master: MPI_Send p1 list error: mnp=%i (mcdetector_out_list_slaves)\n", mpi_node_rank, abs(mnp[0]*mnp[1]*mnp[2]));
    /* slaves are done: sent mnp and p1 */
    return (detector);
  } /* end slaves */

  /* MPI master: receive data from slaves sequentially: 2 MPI_Recv calls */

  if (mpi_node_rank == mpi_node_root) {
    for(node_i=0; node_i<mpi_node_count; node_i++) {
      double *this_p1=NULL;                               /* buffer to hold the list from slaves */
      int     mnp[3]={0,0,0};  /* size of this buffer */
      if (node_i != mpi_node_root) { /* get data from slaves */
	if (mc_MPI_Recv(mnp, 3, MPI_INT, node_i) != MPI_SUCCESS)
	  fprintf(stderr, "Warning: master from proc %i: "
		  "MPI_Recv mnp list error (mcdetector_write_data)\n", node_i);
	if (mnp[0]*mnp[1]*mnp[2]) {
	  this_p1 = (double *)calloc(mnp[0]*mnp[1]*mnp[2], sizeof(double));
	  if (!this_p1 || mc_MPI_Recv(this_p1, abs(mnp[0]*mnp[1]*mnp[2]), MPI_DOUBLE, node_i)!= MPI_SUCCESS)
	    fprintf(stderr, "Warning: master from proc %i: "
		    "MPI_Recv p1 list error: mnp=%i (mcdetector_write_data)\n", node_i, mnp[0]*mnp[1]*mnp[2]);
	  else {
	    printf(". MPI master writing data for slave node %i\n",node_i);
	    detector.p1 = this_p1;
	    detector.m  = mnp[0]; detector.n  = mnp[1]; detector.p  = mnp[2];

	    mcdetector_out_data_nexus(nxhandle, detector);
	  }
	}
      } /* if not master */
    } /* for */
  MPI_MASTER(
	     printf("\n** Done ** \n");
  );
  }
}
#endif

MCDETECTOR mcdetector_out_0D_nexus(MCDETECTOR detector)
{
  /* Write data set information to NeXus file. */
  MPI_MASTER(
    mcdatainfo_out_nexus(nxhandle, detector);
  );

  return(detector);
} /* mcdetector_out_0D_ascii */

MCDETECTOR mcdetector_out_1D_nexus(MCDETECTOR detector)
{
  MPI_MASTER(
  mcdatainfo_out_nexus(nxhandle, detector);
  mcdetector_out_data_nexus(nxhandle, detector);
  );
  return(detector);
} /* mcdetector_out_1D_ascii */

MCDETECTOR mcdetector_out_2D_nexus(MCDETECTOR detector)
{
  MPI_MASTER(
  mcdatainfo_out_nexus(nxhandle, detector);
  mcdetector_out_data_nexus(nxhandle, detector);
  );

#ifdef USE_MPI // and USE_NEXUS
  /* NeXus: slave nodes have master write their lists */
  if (strcasestr(detector.format, "list") && mpi_node_count > 1) {
    mcdetector_out_list_slaves(detector);
  }
#endif /* USE_MPI */

  return(detector);
} /* mcdetector_out_2D_nexus */

#endif /* USE_NEXUS*/








/* ========================================================================== */

/*                            Main input functions                            */
/*            DETECTOR_OUT_xD function calls -> ascii or NeXus                */

/* ========================================================================== */

/*******************************************************************************
* siminfo_init:   open SIM and write header
*******************************************************************************/
FILE *siminfo_init(FILE *f)
{
  int exists=0;
  int index;

  /* check format */
  if (!mcformat || !strlen(mcformat)
   || !strcasecmp(mcformat, "MCSTAS") || !strcasecmp(mcformat, "MCXTRACE")
   || !strcasecmp(mcformat, "PGPLOT") || !strcasecmp(mcformat, "GNUPLOT") || !strcasecmp(mcformat, "MCCODE")
   || !strcasecmp(mcformat, "MATLAB")) {
    mcformat="McCode";
#ifdef USE_NEXUS
  } else if (strcasestr(mcformat, "NeXus")) {
    /* Do nothing */
#endif
  } else {
    fprintf(stderr,
	    "Warning: You have requested the output format %s which is unsupported by this binary. Resetting to standard %s format.\n",mcformat ,"McCode");
    mcformat="McCode";
  }

  /* open the SIM file if not defined yet */
  if (siminfo_file || mcdisable_output_files)
    return (siminfo_file);

#ifdef USE_NEXUS
  /* only master writes NeXus header: calls NXopen(nxhandle) */
  if (mcformat && strcasestr(mcformat, "NeXus")) {
	  MPI_MASTER(
	  siminfo_file = mcnew_file(siminfo_name, "h5", &exists);
    if(!siminfo_file)
      fprintf(stderr,
	      "Warning: could not open simulation description file '%s'\n",
	      siminfo_name);
	  else
	    mcinfo_out_nexus(nxhandle);
	  );
    return(siminfo_file); /* points to nxhandle */
  }
#endif

  /* write main description file (only MASTER) */
  MPI_MASTER(

  siminfo_file = mcnew_file(siminfo_name, "sim", &exists);
  if(!siminfo_file)
    fprintf(stderr,
	    "Warning: could not open simulation description file '%s'\n",
	    siminfo_name);
  else
  {
    /* write SIM header */
    time_t t=time(NULL);
    siminfo_out("%s simulation description file for %s.\n",
      MCCODE_NAME, instrument_name);
    siminfo_out("Date:    %s", ctime(&t)); /* includes \n */
    siminfo_out("Program: %s\n\n", MCCODE_STRING);

    siminfo_out("begin instrument: %s\n", instrument_name);
    mcinfo_out(   "  ", siminfo_file);
    siminfo_out("end instrument\n");

    siminfo_out("\nbegin simulation: %s\n", dirname);
    mcruninfo_out("  ", siminfo_file);
    siminfo_out("end simulation\n");

  }
  return (siminfo_file);

  ); /* MPI_MASTER */

} /* siminfo_init */

/*******************************************************************************
*   siminfo_close:  close SIM
*******************************************************************************/
void siminfo_close()
{
#ifdef USE_MPI
  if(mpi_node_rank == mpi_node_root) {
#endif
  if(siminfo_file && !mcdisable_output_files) {
#ifdef USE_NEXUS
    if (mcformat && strcasestr(mcformat, "NeXus")) {
      time_t t=time(NULL);
      nxprintf(nxhandle, "end_time", ctime(&t));
      nxprintf(nxhandle, "duration", "%li", (long)t-mcstartdate);
      NXclosegroup(nxhandle); /* NXentry */
      NXclose(&nxhandle);
    } else {
#endif
      fclose(siminfo_file);
#ifdef USE_NEXUS
    }
#endif
#ifdef USE_MPI
  }
#endif
    siminfo_file = NULL;
  }
} /* siminfo_close */

/*******************************************************************************
* mcdetector_out_0D: wrapper for 0D (single value).
*   Output single detector/monitor data (p0, p1, p2).
*   Title is t, component name is c.
*******************************************************************************/
MCDETECTOR mcdetector_out_0D(char *t, double p0, double p1, double p2,
                         char *c, Coords posa)
{
  /* import and perform basic detector analysis (and handle MPI reduce) */
  MCDETECTOR detector = detector_import(mcformat,
    c, (t ? t : MCCODE_STRING " data"),
    1, 1, 1,
    "I", "", "",
    "I", "", "",
    0, 0, 0, 0, 0, 0, "",
    &p0, &p1, &p2, posa); /* write Detector: line */

#ifdef USE_NEXUS
  if (strcasestr(detector.format, "NeXus"))
    return(mcdetector_out_0D_nexus(detector));
  else
#endif
    return(mcdetector_out_0D_ascii(detector));

} /* mcdetector_out_0D */



/*******************************************************************************
* mcdetector_out_1D: wrapper for 1D.
*   Output 1d detector data (p0, p1, p2) for n bins linearly
*   distributed across the range x1..x2 (x1 is lower limit of first
*   bin, x2 is upper limit of last bin). Title is t, axis labels are xl
*   and yl. File name is f, component name is c.
*******************************************************************************/
MCDETECTOR mcdetector_out_1D(char *t, char *xl, char *yl,
        char *xvar, double x1, double x2,
        long n,
        double *p0, double *p1, double *p2, char *f,
        char *c, Coords posa)
{
  /* import and perform basic detector analysis (and handle MPI_Reduce) */
  MCDETECTOR detector = detector_import(mcformat,
    c, (t ? t : MCCODE_STRING " 1D data"),
    n, 1, 1,
    xl, yl, (n > 1 ? "Signal per bin" : " Signal"),
    xvar, "(I,I_err)", "I",
    x1, x2, 0, 0, 0, 0, f,
    p0, p1, p2, posa); /* write Detector: line */
  if (!detector.p1 || !detector.m) return(detector);

#ifdef USE_NEXUS
  if (strcasestr(detector.format, "NeXus"))
    return(mcdetector_out_1D_nexus(detector));
  else
#endif
    return(mcdetector_out_1D_ascii(detector));

} /* mcdetector_out_1D */

/*******************************************************************************
* mcdetector_out_2D: wrapper for 2D.
*   special case for list: master creates file first, then slaves append their blocks without header
*******************************************************************************/
MCDETECTOR mcdetector_out_2D(char *t, char *xl, char *yl,
                  double x1, double x2, double y1, double y2,
                  long m, long n,
                  double *p0, double *p1, double *p2, char *f,
                  char *c, Coords posa)
{
  char xvar[CHAR_BUF_LENGTH];
  char yvar[CHAR_BUF_LENGTH];

  /* create short axes labels */
  if (xl && strlen(xl)) { strncpy(xvar, xl, CHAR_BUF_LENGTH); xvar[2]='\0'; }
  else strcpy(xvar, "x");
  if (yl && strlen(yl)) { strncpy(yvar, yl, CHAR_BUF_LENGTH); yvar[2]='\0'; }
  else strcpy(yvar, "y");

  MCDETECTOR detector;

  /* import and perform basic detector analysis (and handle MPI_Reduce) */
  if (abs(m) == 1) {/* n>1 on Y, m==1 on X: 1D, no X axis*/
    detector = detector_import(mcformat,
      c, (t ? t : MCCODE_STRING " 1D data"),
      n, 1, 1,
      yl, "", "Signal per bin",
      yvar, "(I,Ierr)", "I",
      y1, y2, x1, x2, 0, 0, f,
      p0, p1, p2, posa); /* write Detector: line */
  } else if (abs(n)==1) {/* m>1 on X, n==1 on Y: 1D, no Y axis*/
    detector = detector_import(mcformat,
      c, (t ? t : MCCODE_STRING " 1D data"),
      m, 1, 1,
      xl, "", "Signal per bin",
      xvar, "(I,Ierr)", "I",
      x1, x2, y1, y2, 0, 0, f,
      p0, p1, p2, posa); /* write Detector: line */
  }else {
    detector = detector_import(mcformat,
      c, (t ? t : MCCODE_STRING " 2D data"),
      m, n, 1,
      xl, yl, "Signal per bin",
      xvar, yvar, "I",
      x1, x2, y1, y2, 0, 0, f,
      p0, p1, p2, posa); /* write Detector: line */
  }

  if (!detector.p1 || !detector.m) return(detector);

#ifdef USE_NEXUS
  if (strcasestr(detector.format, "NeXus"))
    return(mcdetector_out_2D_nexus(detector));
  else
#endif
    return(mcdetector_out_2D_ascii(detector));

} /* mcdetector_out_2D */

/*******************************************************************************
* mcdetector_out_list: wrapper for list output (calls out_2D with mcformat+"list").
*   m=number of events, n=size of each event
*******************************************************************************/
MCDETECTOR mcdetector_out_list(char *t, char *xl, char *yl,
                  long m, long n,
                  double *p1, char *f,
                  char *c, Coords posa)
{
  char       format_new[CHAR_BUF_LENGTH];
  char      *format_org;
  MCDETECTOR detector;

  format_org = mcformat;
  strcpy(format_new, mcformat);
  strcat(format_new, " list");
  mcformat = format_new;

  detector = mcdetector_out_2D(t, xl, yl,
                  1,abs(m),1,abs(n),
                  m,n,
                  NULL, p1, NULL, f,
                  c, posa);

  mcformat = format_org;
  return(detector);
}

/*******************************************************************************
 * mcuse_dir: set data/sim storage directory and create it,
 * or exit with error if exists
 ******************************************************************************/
static void
mcuse_dir(char *dir)
{
  if (!dir || !strlen(dir)) return;
#ifdef MC_PORTABLE
  fprintf(stderr, "Error: "
          "Directory output cannot be used with portable simulation (mcuse_dir)\n");
  exit(1);
#else  /* !MC_PORTABLE */
  /* handle file://directory URL type */
  if (strncmp(dir, "file://", strlen("file://")))
    dirname = dir;
  else
    dirname = dir+strlen("file://");


#ifdef USE_MPI
  if(mpi_node_rank == mpi_node_root) {
#endif
    if(mkdir(dirname, 0777)) {
#ifndef DANSE
      fprintf(stderr, "Error: unable to create directory '%s' (mcuse_dir)\n", dir);
      fprintf(stderr, "(Maybe the directory already exists?)\n");
#endif
#ifdef USE_MPI
      MPI_Abort(MPI_COMM_WORLD, -1);
#endif
    exit(-1);
    }
#ifdef USE_MPI
    }
#endif

  /* remove trailing PATHSEP (if any) */
  while (strlen(dirname) && dirname[strlen(dirname) - 1] == MC_PATHSEP_C)
    dirname[strlen(dirname) - 1]='\0';
#endif /* !MC_PORTABLE */
} /* mcuse_dir */

/*******************************************************************************
* mcinfo: display instrument simulation info to stdout and exit
*******************************************************************************/
static void
mcinfo(void)
{
  fprintf(stdout, "begin instrument: %s\n", instrument_name);
  mcinfo_out("  ", stdout);
  fprintf(stdout, "end instrument\n");
  fprintf(stdout, "begin simulation: %s\n", dirname ? dirname : ".");
  mcruninfo_out("  ", stdout);
  fprintf(stdout, "end simulation\n");
  exit(0); /* includes MPI_Finalize in MPI mode */
} /* mcinfo */

#endif /* ndef MCCODE_R_IO_C */

/* end of the I/O section =================================================== */







/*******************************************************************************
* mcset_ncount: set total number of rays to generate
*******************************************************************************/
void mcset_ncount(unsigned long long int count)
{
  mcncount = count;
}

/* mcget_ncount: get total number of rays to generate */
#pragma acc routine seq
unsigned long long int mcget_ncount(void)
{
  return mcncount;
}

/* mcget_run_num: get curent number of rays in TRACE */
#pragma acc routine seq
unsigned long long int mcget_run_num(void)
{
  //FIXME!!
  return 100000;//mcrun_num;
}

/* mcsetn_arg: get ncount from a string argument */
static void
mcsetn_arg(char *arg)
{
  mcset_ncount((long long int) strtod(arg, NULL));
}

/* mcsetseed: set the random generator seed from a string argument */
static void
mcsetseed(char *arg)
{
  mcseed = atol(arg);
  if(mcseed) {
    srandom(mcseed);
  } else {
    fprintf(stderr, "Error: seed must not be zero (mcsetseed)\n");
    exit(1);
  }
}

/* Following part is only embedded when not redundent with mccode-r.h ========= */

#ifndef MCCODE_H

/* SECTION: MCDISPLAY support. =============================================== */

/*******************************************************************************
* Just output MCDISPLAY keywords to be caught by an external plotter client.
*******************************************************************************/

void mcdis_magnify(char *what){
  // Do nothing here, better use interactive zoom from the tools
}

void mcdis_line(double x1, double y1, double z1,
                double x2, double y2, double z2){
  printf("MCDISPLAY: multiline(2,%g,%g,%g,%g,%g,%g)\n",
         x1,y1,z1,x2,y2,z2);
}

void mcdis_dashed_line(double x1, double y1, double z1,
		       double x2, double y2, double z2, int n){
  int i;
  const double dx = (x2-x1)/(2*n+1);
  const double dy = (y2-y1)/(2*n+1);
  const double dz = (z2-z1)/(2*n+1);

  for(i = 0; i < n+1; i++)
    mcdis_line(x1 + 2*i*dx,     y1 + 2*i*dy,     z1 + 2*i*dz,
	       x1 + (2*i+1)*dx, y1 + (2*i+1)*dy, z1 + (2*i+1)*dz);
}

void mcdis_multiline(int count, ...){
  va_list ap;
  double x,y,z;

  printf("MCDISPLAY: multiline(%d", count);
  va_start(ap, count);
  while(count--)
    {
    x = va_arg(ap, double);
    y = va_arg(ap, double);
    z = va_arg(ap, double);
    printf(",%g,%g,%g", x, y, z);
    }
  va_end(ap);
  printf(")\n");
}

void mcdis_rectangle(char* plane, double x, double y, double z,
		     double width, double height){
  /* draws a rectangle in the plane           */
  /* x is ALWAYS width and y is ALWAYS height */
  if (strcmp("xy", plane)==0) {
    mcdis_multiline(5,
		    x - width/2, y - height/2, z,
		    x + width/2, y - height/2, z,
		    x + width/2, y + height/2, z,
		    x - width/2, y + height/2, z,
		    x - width/2, y - height/2, z);
  } else if (strcmp("xz", plane)==0) {
    mcdis_multiline(5,
		    x - width/2, y, z - height/2,
		    x + width/2, y, z - height/2,
		    x + width/2, y, z + height/2,
		    x - width/2, y, z + height/2,
		    x - width/2, y, z - height/2);
  } else if (strcmp("yz", plane)==0) {
    mcdis_multiline(5,
		    x, y - height/2, z - width/2,
		    x, y - height/2, z + width/2,
		    x, y + height/2, z + width/2,
		    x, y + height/2, z - width/2,
		    x, y - height/2, z - width/2);
  } else {

    fprintf(stderr, "Error: Definition of plane %s unknown\n", plane);
    exit(1);
  }
}

/*  draws a box with center at (x, y, z) and
    width (deltax), height (deltay), length (deltaz) */
void mcdis_box(double x, double y, double z,
	       double width, double height, double length){

  mcdis_rectangle("xy", x, y, z-length/2, width, height);
  mcdis_rectangle("xy", x, y, z+length/2, width, height);
  mcdis_line(x-width/2, y-height/2, z-length/2,
	     x-width/2, y-height/2, z+length/2);
  mcdis_line(x-width/2, y+height/2, z-length/2,
	     x-width/2, y+height/2, z+length/2);
  mcdis_line(x+width/2, y-height/2, z-length/2,
	     x+width/2, y-height/2, z+length/2);
  mcdis_line(x+width/2, y+height/2, z-length/2,
	     x+width/2, y+height/2, z+length/2);
}

void mcdis_circle(char *plane, double x, double y, double z, double r){
  printf("MCDISPLAY: circle('%s',%g,%g,%g,%g)\n", plane, x, y, z, r);
}

/* Draws a circle with center (x,y,z), radius (r), and in the plane
 * with normal (nx,ny,nz)*/
void mcdis_Circle(double x, double y, double z, double r, double nx, double ny, double nz){
    int i;
    if(nx==0 && ny && nz==0){
        for (i=0;i<24; i++){
            mcdis_line(x+r*sin(i*2*PI/24),y,z+r*cos(i*2*PI/24),
                    x+r*sin((i+1)*2*PI/24),y,z+r*cos((i+1)*2*PI/24));
        }
    }else{
        double mx,my,mz;
        /*generate perpendicular vector using (nx,ny,nz) and (0,1,0)*/
        vec_prod(mx,my,mz, 0,1,0, nx,ny,nz);
        NORM(mx,my,mz);
        /*draw circle*/
        for (i=0;i<24; i++){
            double ux,uy,uz;
            double wx,wy,wz;
            rotate(ux,uy,uz, mx,my,mz, i*2*PI/24, nx,ny,nz);
            rotate(wx,wy,wz, mx,my,mz, (i+1)*2*PI/24, nx,ny,nz);
            mcdis_line(x+ux*r,y+uy*r,z+uz*r,
                    x+wx*r,y+wy*r,z+wz*r);
        }
    }
}

/* Draws a cylinder with center at (x,y,z) with extent (r,height).
 * The cylinder axis is along the vector nx,ny,nz.
 * N determines how many vertical lines are drawn.*/
void mcdis_cylinder( double x, double y, double z,
        double r, double height, int N, double nx, double ny, double nz){
    int i;
    /*no lines make little sense - so trigger the default*/
    if(N<=0) N=5;

    NORM(nx,ny,nz);
    double h_2=height/2.0;
    mcdis_Circle(x+nx*h_2,y+ny*h_2,z+nz*h_2,r,nx,ny,nz);
    mcdis_Circle(x-nx*h_2,y-ny*h_2,z-nz*h_2,r,nx,ny,nz);

    double mx,my,mz;
    /*generate perpendicular vector using (nx,ny,nz) and (0,1,0)*/
    if(nx==0 && ny && nz==0){
        mx=my=0;mz=1;
    }else{
        vec_prod(mx,my,mz, 0,1,0, nx,ny,nz);
        NORM(mx,my,mz);
    }
    /*draw circle*/
    for (i=0; i<24; i++){
        double ux,uy,uz;
        rotate(ux,uy,uz, mx,my,mz, i*2*PI/24, nx,ny,nz);
        mcdis_line(x+nx*h_2+ux*r, y+ny*h_2+uy*r, z+nz*h_2+uz*r,
                 x-nx*h_2+ux*r, y-ny*h_2+uy*r, z-nz*h_2+uz*r);
    }
}

/* draws a sphere with center at (x,y,z) with extent (r)
 * The sphere is drawn using N longitudes and N latitudes.*/
void mcdis_sphere(double x, double y, double z, double r, int N){
    double nx,ny,nz;
    int i;
    /*no lines make little sense - so trigger the default*/
    if(N<=0) N=5;

    nx=0;ny=0;nz=1;
    mcdis_Circle(x,y,z,r,nx,ny,nz);
    for (i=1;i<N;i++){
        rotate(nx,ny,nz, nx,ny,nz, PI/N, 0,1,0);
        mcdis_Circle(x,y,z,r,nx,ny,nz);
    }
    /*lastly draw a great circle perpendicular to all N circles*/
    //mcdis_Circle(x,y,z,radius,1,0,0);

    for (i=1;i<=N;i++){
        double yy=-r+ 2*r*((double)i/(N+1));
        mcdis_Circle(x,y+yy ,z,  sqrt(r*r-yy*yy) ,0,1,0);
    }
}

/* SECTION: coordinates handling ============================================ */

/*******************************************************************************
* Since we use a lot of geometric calculations using Cartesian coordinates,
* we collect some useful routines here. However, it is also permissible to
* work directly on the underlying struct coords whenever that is most
* convenient (that is, the type Coords is not abstract).
*
* Coordinates are also used to store rotation angles around x/y/z axis.
*
* Since coordinates are used much like a basic type (such as double), the
* structure itself is passed and returned, rather than a pointer.
*
* At compile-time, the values of the coordinates may be unknown (for example
* a motor position). Hence coordinates are general expressions and not simple
* numbers. For this we used the type Coords_exp which has three CExp
* fields. For runtime (or calculations possible at compile time), we use
* Coords which contains three double fields.
*******************************************************************************/

/* coords_set: Assign coordinates. */
#pragma acc routine seq
Coords coords_set(MCNUM x, MCNUM y, MCNUM z)
{
  Coords a;

  a.x = x;
  a.y = y;
  a.z = z;
  return a;
}

/* coords_get: get coordinates. Required when 'x','y','z' are #defined as ray pars */
#pragma acc routine seq
Coords coords_get(Coords a, MCNUM *x, MCNUM *y, MCNUM *z)
{
  *x = a.x;
  *y = a.y;
  *z = a.z;
  return a;
}

/* coords_add: Add two coordinates. */
#pragma acc routine seq
Coords coords_add(Coords a, Coords b)
{
  Coords c;

  c.x = a.x + b.x;
  c.y = a.y + b.y;
  c.z = a.z + b.z;
  if (fabs(c.z) < 1e-14) c.z=0.0;
  return c;
}

/* coords_sub: Subtract two coordinates. */
#pragma acc routine seq
Coords coords_sub(Coords a, Coords b)
{
  Coords c;

  c.x = a.x - b.x;
  c.y = a.y - b.y;
  c.z = a.z - b.z;
  if (fabs(c.z) < 1e-14) c.z=0.0;
  return c;
}

/* coords_neg: Negate coordinates. */
#pragma acc routine seq
Coords coords_neg(Coords a)
{
  Coords b;

  b.x = -a.x;
  b.y = -a.y;
  b.z = -a.z;
  return b;
}

/* coords_scale: Scale a vector. */
#pragma acc routine seq
Coords coords_scale(Coords b, double scale) {
  Coords a;

  a.x = b.x*scale;
  a.y = b.y*scale;
  a.z = b.z*scale;
  return a;
}

/* coords_sp: Scalar product: a . b */
#pragma acc routine seq
double coords_sp(Coords a, Coords b) {
  double value;

  value = a.x*b.x + a.y*b.y + a.z*b.z;
  return value;
}

/* coords_xp: Cross product: a = b x c. */
#pragma acc routine seq
Coords coords_xp(Coords b, Coords c) {
  Coords a;

  a.x = b.y*c.z - c.y*b.z;
  a.y = b.z*c.x - c.z*b.x;
  a.z = b.x*c.y - c.x*b.y;
  return a;
}

/* coords_len: Gives length of coords set. */
#pragma acc routine seq
double coords_len(Coords a) {
  return sqrt(a.x*a.x + a.y*a.y + a.z*a.z);
}

/* coords_mirror: Mirror a in plane (through the origin) defined by normal n*/
#pragma acc routine seq
Coords coords_mirror(Coords a, Coords n) {
  double t = scalar_prod(n.x, n.y, n.z, n.x, n.y, n.z);
  Coords b;
  if (t!=1) {
    t = sqrt(t);
    n.x /= t;
    n.y /= t;
    n.z /= t;
  }
  t=scalar_prod(a.x, a.y, a.z, n.x, n.y, n.z);
  b.x = a.x-2*t*n.x;
  b.y = a.y-2*t*n.y;
  b.z = a.z-2*t*n.z;
  return b;
}

/* coords_print: Print out vector values. */
#pragma acc routine seq
void coords_print(Coords a) {

  //  fprintf(stdout, "(%f, %f, %f)\n", a.x, a.y, a.z);
  return;
}

mcstatic void coords_norm(Coords* c) {
	double temp = coords_sp(*c,*c);

	// Skip if we will end dividing by zero
	if (temp == 0) return;

	temp = sqrt(temp);

	c->x /= temp;
	c->y /= temp;
	c->z /= temp;
}

/*******************************************************************************
* The Rotation type implements a rotation transformation of a coordinate
* system in the form of a double[3][3] matrix.
*
* Contrary to the Coords type in coords.c, rotations are passed by
* reference. Functions that yield new rotations do so by writing to an
* explicit result parameter; rotations are not returned from functions. The
* reason for this is that arrays cannot by returned from functions (though
* structures can; thus an alternative would have been to wrap the
* double[3][3] array up in a struct). Such are the ways of C programming.
*
* A rotation represents the tranformation of the coordinates of a vector when
* changing between coordinate systems that are rotated with respect to each
* other. For example, suppose that coordinate system Q is rotated 45 degrees
* around the Z axis with respect to coordinate system P. Let T be the
* rotation transformation representing a 45 degree rotation around Z. Then to
* get the coordinates of a vector r in system Q, apply T to the coordinates
* of r in P. If r=(1,0,0) in P, it will be (sqrt(1/2),-sqrt(1/2),0) in
* Q. Thus we should be careful when interpreting the sign of rotation angles:
* they represent the rotation of the coordinate systems, not of the
* coordinates (which has opposite sign).
*******************************************************************************/

/*******************************************************************************
* rot_set_rotation: Get transformation for rotation first phx around x axis,
* then phy around y, then phz around z.
*******************************************************************************/
#pragma acc routine seq
void rot_set_rotation(Rotation t, double phx, double phy, double phz)
{
  if ((phx == 0) && (phy == 0) && (phz == 0)) {
    t[0][0] = 1.0;
    t[0][1] = 0.0;
    t[0][2] = 0.0;
    t[1][0] = 0.0;
    t[1][1] = 1.0;
    t[1][2] = 0.0;
    t[2][0] = 0.0;
    t[2][1] = 0.0;
    t[2][2] = 1.0;
  } else {
    double cx = cos(phx);
    double sx = sin(phx);
    double cy = cos(phy);
    double sy = sin(phy);
    double cz = cos(phz);
    double sz = sin(phz);

    t[0][0] = cy*cz;
    t[0][1] = sx*sy*cz + cx*sz;
    t[0][2] = sx*sz - cx*sy*cz;
    t[1][0] = -cy*sz;
    t[1][1] = cx*cz - sx*sy*sz;
    t[1][2] = sx*cz + cx*sy*sz;
    t[2][0] = sy;
    t[2][1] = -sx*cy;
    t[2][2] = cx*cy;
  }
}

/*******************************************************************************
* rot_test_identity: Test if rotation is identity
*******************************************************************************/
#pragma acc routine seq
int rot_test_identity(Rotation t)
{
  return (t[0][0] + t[1][1] + t[2][2] == 3);
}

/*******************************************************************************
* rot_mul: Matrix multiplication of transformations (this corresponds to
* combining transformations). After rot_mul(T1, T2, T3), doing T3 is
* equal to doing first T2, then T1.
* Note that T3 must not alias (use the same array as) T1 or T2.
*******************************************************************************/
#pragma acc routine seq
void rot_mul(Rotation t1, Rotation t2, Rotation t3)
{
  if (rot_test_identity(t1)) {
    rot_copy(t3, t2);
  } else if (rot_test_identity(t2)) {
    rot_copy(t3, t1);
  } else {
    int i,j;
    for(i = 0; i < 3; i++)
      for(j = 0; j < 3; j++)
	t3[i][j] = t1[i][0]*t2[0][j] + t1[i][1]*t2[1][j] + t1[i][2]*t2[2][j];
  }
}

/*******************************************************************************
* rot_copy: Copy a rotation transformation (arrays cannot be assigned in C).
*******************************************************************************/
#pragma acc routine seq
void rot_copy(Rotation dest, Rotation src)
{
  int i,j;
  for(i = 0; i < 3; i++)
    for(j = 0; j < 3; j++)
      dest[i][j] = src[i][j];
}

/*******************************************************************************
* rot_transpose: Matrix transposition, which is inversion for Rotation matrices
*******************************************************************************/
#pragma acc routine seq
void rot_transpose(Rotation src, Rotation dst)
{
  dst[0][0] = src[0][0];
  dst[0][1] = src[1][0];
  dst[0][2] = src[2][0];
  dst[1][0] = src[0][1];
  dst[1][1] = src[1][1];
  dst[1][2] = src[2][1];
  dst[2][0] = src[0][2];
  dst[2][1] = src[1][2];
  dst[2][2] = src[2][2];
}

/*******************************************************************************
* rot_apply: returns t*a
*******************************************************************************/
#pragma acc routine seq
Coords rot_apply(Rotation t, Coords a)
{
  Coords b;
  if (rot_test_identity(t)) {
    return a;
  } else {
    b.x = t[0][0]*a.x + t[0][1]*a.y + t[0][2]*a.z;
    b.y = t[1][0]*a.x + t[1][1]*a.y + t[1][2]*a.z;
    b.z = t[2][0]*a.x + t[2][1]*a.y + t[2][2]*a.z;
    return b;
  }
}

/**
 * Pretty-printing of rotation matrices.
 */
void rot_print(Rotation rot) {
	printf("[ %4.2f %4.2f %4.2f ]\n",
			rot[0][0], rot[0][1], rot[0][2]);
	printf("[ %4.2f %4.2f %4.2f ]\n",
			rot[1][0], rot[1][1], rot[1][2]);
	printf("[ %4.2f %4.2f %4.2f ]\n\n",
			rot[2][0], rot[2][1], rot[2][2]);
}

/**
 * Vector product: used by vec_prod (mccode-r.h). Use coords_xp for Coords.
 */
#pragma acc routine seq
void vec_prod_func(double *x, double *y, double *z,
		double x1, double y1, double z1,
		double x2, double y2, double z2) {
    *x = (y1)*(z2) - (y2)*(z1);
    *y = (z1)*(x2) - (z2)*(x1);
    *z = (x1)*(y2) - (x2)*(y1);
}

/**
 * Scalar product: use coords_sp for Coords.
 */
#pragma acc routine seq
double scalar_prod(
		double x1, double y1, double z1,
		double x2, double y2, double z2) {
	return ((x1 * x2) + (y1 * y2) + (z1 * z2));
}

#pragma acc routine seq
mcstatic void norm_func(double *x, double *y, double *z) {
	double temp = (*x * *x) + (*y * *y) + (*z * *z);
	if (temp != 0) {
		temp = sqrt(temp);
		*x /= temp;
		*y /= temp;
		*z /= temp;
	}
}

/*******************************************************************************
* mccoordschange: applies rotation to (x y z) and (vx vy vz) and Spin (sx,sy,sz)
*******************************************************************************/
#pragma acc routine seq
void mccoordschange(Coords a, Rotation t, _class_particle *particle)
{
  Coords b, c;

  b.x = particle->x;
  b.y = particle->y;
  b.z = particle->z;
  c = rot_apply(t, b);
  b = coords_add(c, a);
  particle->x = b.x;
  particle->y = b.y;
  particle->z = b.z;

#if MCCODE_PARTICLE_CODE == 2112
    if (particle->vz != 0.0 || particle->vx != 0.0 || particle->vy != 0.0)
      mccoordschange_polarisation(t, &(particle->vx), &(particle->vy), &(particle->vz));

    if (particle->sz != 0.0 || particle->sx != 0.0 || particle->sy != 0.0)
      mccoordschange_polarisation(t, &(particle->sx), &(particle->sy), &(particle->sz));
#elif MCCODE_PARTICLE_CODE == 22
    if (particle->kz != 0.0 || particle->kx != 0.0 || particle->ky != 0.0)
      mccoordschange_polarisation(t, &(particle->kx), &(particle->ky), &(particle->kz));

    if (particle->Ez != 0.0 || particle->Ex != 0.0 || particle->Ey != 0.0)
      mccoordschange_polarisation(t, &(particle->Ex), &(particle->Ey), &(particle->Ez));
#endif
}

/*******************************************************************************
* mccoordschange_polarisation: applies rotation to vector (sx sy sz)
*******************************************************************************/
#pragma acc routine seq
void mccoordschange_polarisation(Rotation t, double *sx, double *sy, double *sz)
{
  Coords b, c;

  b.x = *sx;
  b.y = *sy;
  b.z = *sz;
  c = rot_apply(t, b);
  *sx = c.x;
  *sy = c.y;
  *sz = c.z;
}

/* SECTION: vector math  ==================================================== */

/* normal_vec_func: Compute normal vector to (x,y,z). */
#pragma acc routine seq
void normal_vec(double *nx, double *ny, double *nz,
                double x, double y, double z)
{
  double ax = fabs(x);
  double ay = fabs(y);
  double az = fabs(z);
  double l;
  if(x == 0 && y == 0 && z == 0)
  {
    *nx = 0;
    *ny = 0;
    *nz = 0;
    return;
  }
  if(ax < ay)
  {
    if(ax < az)
    {                           /* Use X axis */
      l = sqrt(z*z + y*y);
      *nx = 0;
      *ny = z/l;
      *nz = -y/l;
      return;
    }
  }
  else
  {
    if(ay < az)
    {                           /* Use Y axis */
      l = sqrt(z*z + x*x);
      *nx = z/l;
      *ny = 0;
      *nz = -x/l;
      return;
    }
  }
  /* Use Z axis */
  l = sqrt(y*y + x*x);
  *nx = y/l;
  *ny = -x/l;
  *nz = 0;
} /* normal_vec */

/*******************************************************************************
 * solve_2nd_order: second order equation solve: A*t^2 + B*t + C = 0
 * solve_2nd_order(&t1, NULL, A,B,C)
 *   returns 0 if no solution was found, or set 't1' to the smallest positive
 *   solution.
 * solve_2nd_order(&t1, &t2, A,B,C)
 *   same as with &t2=NULL, but also returns the second solution.
 * EXAMPLE usage for intersection of a trajectory with a plane in gravitation
 * field (gx,gy,gz):
 * The neutron starts at point r=(x,y,z) with velocityv=(vx vy vz). The plane
 * has a normal vector n=(nx,ny,nz) and contains the point W=(wx,wy,wz).
 * The problem consists in solving the 2nd order equation:
 *      1/2.n.g.t^2 + n.v.t + n.(r-W) = 0
 * so that A = 0.5 n.g; B = n.v; C = n.(r-W);
 * Without acceleration, t=-n.(r-W)/n.v
 ******************************************************************************/
#pragma acc routine seq
int solve_2nd_order(double *t1, double *t2,
                  double A,  double B,  double C)
{
  int ret=0;

  if (!t1) return 0;
  *t1 = 0;
  if (t2) *t2=0;

  if (fabs(A) < 1E-10) /* approximate to linear equation: A ~ 0 */
  {
    if (B) {  *t1 = -C/B; ret=1; if (t2) *t2=*t1; }
    /* else no intersection: A=B=0 ret=0 */
  }
  else
  {
    double D;
    D = B*B - 4*A*C;
    if (D >= 0) /* Delta > 0: two solutions */
    {
      double sD, dt1, dt2;
      sD = sqrt(D);
      dt1 = (-B + sD)/2/A;
      dt2 = (-B - sD)/2/A;
      /* we identify very small values with zero */
      if (fabs(dt1) < 1e-10) dt1=0.0;
      if (fabs(dt2) < 1e-10) dt2=0.0;

      /* now we choose the smallest positive solution */
      if      (dt1<=0.0 && dt2>0.0) ret=2; /* dt2 positive */
      else if (dt2<=0.0 && dt1>0.0) ret=1; /* dt1 positive */
      else if (dt1> 0.0 && dt2>0.0)
      {  if (dt1 < dt2) ret=1; else ret=2; } /* all positive: min(dt1,dt2) */
      /* else two solutions are negative. ret=-1 */
      if (ret==1) { *t1 = dt1;  if (t2) *t2=dt2; }
      else        { *t1 = dt2;  if (t2) *t2=dt1; }
      ret=2;  /* found 2 solutions and t1 is the positive one */
    } /* else Delta <0: no intersection. ret=0 */
  }
  return(ret);
} /* solve_2nd_order */

/*******************************************************************************
 * randvec_target_circle: Choose random direction towards target at (x,y,z)
 * with given radius.
 * If radius is zero, choose random direction in full 4PI, no target.
 ******************************************************************************/
#ifndef USE_PGI
void randvec_target_circle_cpu(double *xo, double *yo, double *zo, double *solid_angle,
        double xi, double yi, double zi, double radius)
{
  double l2, phi, theta, nx, ny, nz, xt, yt, zt, xu, yu, zu;

  if(radius == 0.0)
  {
    /* No target, choose uniformly a direction in full 4PI solid angle. */
    theta = acos (1 - rand0max(2));
    phi = rand0max(2 * PI);
    if(solid_angle)
      *solid_angle = 4*PI;
    nx = 1;
    ny = 0;
    nz = 0;
    yi = sqrt(xi*xi+yi*yi+zi*zi);
    zi = 0;
    xi = 0;
  }
  else
  {
    double costheta0;
    l2 = xi*xi + yi*yi + zi*zi; /* sqr Distance to target. */
    costheta0 = sqrt(l2/(radius*radius+l2));
    if (radius < 0) costheta0 *= -1;
    if(solid_angle)
    {
      /* Compute solid angle of target as seen from origin. */
        *solid_angle = 2*PI*(1 - costheta0);
    }

    /* Now choose point uniformly on circle surface within angle theta0 */
    theta = acos (1 - rand0max(1 - costheta0)); /* radius on circle */
    phi = rand0max(2 * PI); /* rotation on circle at given radius */
    /* Now, to obtain the desired vector rotate (xi,yi,zi) angle theta around a
       perpendicular axis u=i x n and then angle phi around i. */
    if(xi == 0 && zi == 0)
    {
      nx = 1;
      ny = 0;
      nz = 0;
    }
    else
    {
      nx = -zi;
      nz = xi;
      ny = 0;
    }
  }

  /* [xyz]u = [xyz]i x n[xyz] (usually vertical) */
  vec_prod(xu,  yu,  zu, xi, yi, zi,        nx, ny, nz);
  /* [xyz]t = [xyz]i rotated theta around [xyz]u */
  rotate  (xt,  yt,  zt, xi, yi, zi, theta, xu, yu, zu);
  /* [xyz]o = [xyz]t rotated phi around n[xyz] */
  rotate (*xo, *yo, *zo, xt, yt, zt, phi, xi, yi, zi);
}
#endif
#pragma acc routine seq nohost
void randvec_target_circle_gpu(double *xo, double *yo, double *zo, double *solid_angle,
        double xi, double yi, double zi, double radius,
        _class_particle* _particle)
{
  double l2, phi, theta, nx, ny, nz, xt, yt, zt, xu, yu, zu;

  if(radius == 0.0)
  {
    /* No target, choose uniformly a direction in full 4PI solid angle. */
    theta = acos (1 - rand0max(2));
    phi = rand0max(2 * PI);
    if(solid_angle)
      *solid_angle = 4*PI;
    nx = 1;
    ny = 0;
    nz = 0;
    yi = sqrt(xi*xi+yi*yi+zi*zi);
    zi = 0;
    xi = 0;
  }
  else
  {
    double costheta0;
    l2 = xi*xi + yi*yi + zi*zi; /* sqr Distance to target. */
    costheta0 = sqrt(l2/(radius*radius+l2));
    if (radius < 0) costheta0 *= -1;
    if(solid_angle)
    {
      /* Compute solid angle of target as seen from origin. */
        *solid_angle = 2*PI*(1 - costheta0);
    }

    /* Now choose point uniformly on circle surface within angle theta0 */
    theta = acos (1 - rand0max(1 - costheta0)); /* radius on circle */
    phi = rand0max(2 * PI); /* rotation on circle at given radius */
    /* Now, to obtain the desired vector rotate (xi,yi,zi) angle theta around a
       perpendicular axis u=i x n and then angle phi around i. */
    if(xi == 0 && zi == 0)
    {
      nx = 1;
      ny = 0;
      nz = 0;
    }
    else
    {
      nx = -zi;
      nz = xi;
      ny = 0;
    }
  }

  /* [xyz]u = [xyz]i x n[xyz] (usually vertical) */
  vec_prod(xu,  yu,  zu, xi, yi, zi,        nx, ny, nz);
  /* [xyz]t = [xyz]i rotated theta around [xyz]u */
  rotate  (xt,  yt,  zt, xi, yi, zi, theta, xu, yu, zu);
  /* [xyz]o = [xyz]t rotated phi around n[xyz] */
  rotate (*xo, *yo, *zo, xt, yt, zt, phi, xi, yi, zi);
}
/* randvec_target_circle */

/*******************************************************************************
 * randvec_target_rect_angular: Choose random direction towards target at
 * (xi,yi,zi) with given ANGULAR dimension height x width. height=phi_x=[0,PI],
 * width=phi_y=[0,2*PI] (radians)
 * If height or width is zero, choose random direction in full 4PI, no target.
 *******************************************************************************/
#ifndef USE_PGI
void randvec_target_rect_angular_cpu(double *xo, double *yo, double *zo, double *solid_angle,
               double xi, double yi, double zi, double width, double height, Rotation A)
{
  double theta, phi, nx, ny, nz, xt, yt, zt, xu, yu, zu;
  Coords tmp;
  Rotation Ainverse;

  rot_transpose(A, Ainverse);

  if(height == 0.0 || width == 0.0)
  {
    randvec_target_circle(xo, yo, zo, solid_angle,
               xi, yi, zi, 0);
    return;
  }
  else
  {
    if(solid_angle)
    {
      /* Compute solid angle of target as seen from origin. */
      *solid_angle = 2*fabs(width*sin(height/2));
    }

    /* Go to global coordinate system */

    tmp = coords_set(xi, yi, zi);
    tmp = rot_apply(Ainverse, tmp);
    coords_get(tmp, &xi, &yi, &zi);

    /* Now choose point uniformly on the unit sphere segment with angle theta/phi */
    phi   = width*randpm1()/2.0;
    theta = asin(randpm1()*sin(height/2.0));
    /* Now, to obtain the desired vector rotate (xi,yi,zi) angle theta around
       n, and then phi around u. */
    if(xi == 0 && zi == 0)
    {
      nx = 1;
      ny = 0;
      nz = 0;
    }
    else
    {
      nx = -zi;
      nz = xi;
      ny = 0;
    }
  }

  /* [xyz]u = [xyz]i x n[xyz] (usually vertical) */
  vec_prod(xu,  yu,  zu, xi, yi, zi,        nx, ny, nz);
  /* [xyz]t = [xyz]i rotated theta around [xyz]u */
  rotate  (xt,  yt,  zt, xi, yi, zi, theta, nx, ny, nz);
  /* [xyz]o = [xyz]t rotated phi around n[xyz] */
  rotate (*xo, *yo, *zo, xt, yt, zt, phi, xu,  yu,  zu);

  /* Go back to local coordinate system */
  tmp = coords_set(*xo, *yo, *zo);
  tmp = rot_apply(A, tmp);
  coords_get(tmp, &*xo, &*yo, &*zo);

}
#endif
#pragma acc routine seq nohost
void randvec_target_rect_angular_gpu(double *xo, double *yo, double *zo, double *solid_angle,
        double xi, double yi, double zi, double width, double height, Rotation A,
        _class_particle* _particle)
{
  double theta, phi, nx, ny, nz, xt, yt, zt, xu, yu, zu;
  Coords tmp;
  Rotation Ainverse;

  rot_transpose(A, Ainverse);

  if(height == 0.0 || width == 0.0)
  {
    randvec_target_circle(xo, yo, zo, solid_angle,
               xi, yi, zi, 0);
    return;
  }
  else
  {
    if(solid_angle)
    {
      /* Compute solid angle of target as seen from origin. */
      *solid_angle = 2*fabs(width*sin(height/2));
    }

    /* Go to global coordinate system */

    tmp = coords_set(xi, yi, zi);
    tmp = rot_apply(Ainverse, tmp);
    coords_get(tmp, &xi, &yi, &zi);

    /* Now choose point uniformly on the unit sphere segment with angle theta/phi */
    phi   = width*randpm1()/2.0;
    theta = asin(randpm1()*sin(height/2.0));
    /* Now, to obtain the desired vector rotate (xi,yi,zi) angle theta around
       n, and then phi around u. */
    if(xi == 0 && zi == 0)
    {
      nx = 1;
      ny = 0;
      nz = 0;
    }
    else
    {
      nx = -zi;
      nz = xi;
      ny = 0;
    }
  }

  /* [xyz]u = [xyz]i x n[xyz] (usually vertical) */
  vec_prod(xu,  yu,  zu, xi, yi, zi,        nx, ny, nz);
  /* [xyz]t = [xyz]i rotated theta around [xyz]u */
  rotate  (xt,  yt,  zt, xi, yi, zi, theta, nx, ny, nz);
  /* [xyz]o = [xyz]t rotated phi around n[xyz] */
  rotate (*xo, *yo, *zo, xt, yt, zt, phi, xu,  yu,  zu);

  /* Go back to local coordinate system */
  tmp = coords_set(*xo, *yo, *zo);
  tmp = rot_apply(A, tmp);
  coords_get(tmp, &*xo, &*yo, &*zo);

}
/* randvec_target_rect_angular */

/*******************************************************************************
 * randvec_target_rect_real: Choose random direction towards target at (xi,yi,zi)
 * with given dimension height x width (in meters !).
 *
 * Local emission coordinate is taken into account and corrected for 'order' times.
 * (See remarks posted to mcstas-users by George Apostolopoulus <gapost@ipta.demokritos.gr>)
 *
 * If height or width is zero, choose random direction in full 4PI, no target.
 *
 * Traditionally, this routine had the name randvec_target_rect - this is now a
 * a define (see mcstas-r.h) pointing here. If you use the old rouine, you are NOT
 * taking the local emmission coordinate into account.
*******************************************************************************/
#ifndef USE_PGI
void randvec_target_rect_real_cpu(double *xo, double *yo, double *zo, double *solid_angle,
        double xi, double yi, double zi,
        double width, double height, Rotation A,
        double lx, double ly, double lz, int order)
{
  double dx, dy, dist, dist_p, nx, ny, nz, mx, my, mz, n_norm, m_norm;
  double cos_theta;
  Coords tmp;
  Rotation Ainverse;

  rot_transpose(A, Ainverse);

  if(height == 0.0 || width == 0.0)
  {
    randvec_target_circle(xo, yo, zo, solid_angle,
               xi, yi, zi, 0);
    return;
  }
  else
  {
    /* Now choose point uniformly on rectangle within width x height */
    dx = width*randpm1()/2.0;
    dy = height*randpm1()/2.0;

    /* Determine distance to target plane*/
    dist = sqrt(xi*xi + yi*yi + zi*zi);
    /* Go to global coordinate system */

    tmp = coords_set(xi, yi, zi);
    tmp = rot_apply(Ainverse, tmp);
    coords_get(tmp, &xi, &yi, &zi);

    /* Determine vector normal to trajectory axis (z) and gravity [0 1 0] */
    vec_prod(nx, ny, nz, xi, yi, zi, 0, 1, 0);

    /* This now defines the x-axis, normalize: */
    n_norm=sqrt(nx*nx + ny*ny + nz*nz);
    nx = nx/n_norm;
    ny = ny/n_norm;
    nz = nz/n_norm;

    /* Now, determine our y-axis (vertical in many cases...) */
    vec_prod(mx, my, mz, xi, yi, zi, nx, ny, nz);
    m_norm=sqrt(mx*mx + my*my + mz*mz);
    mx = mx/m_norm;
    my = my/m_norm;
    mz = mz/m_norm;

    /* Our output, random vector can now be defined by linear combination: */

    *xo = xi + dx * nx + dy * mx;
    *yo = yi + dx * ny + dy * my;
    *zo = zi + dx * nz + dy * mz;

    /* Go back to local coordinate system */
    tmp = coords_set(*xo, *yo, *zo);
    tmp = rot_apply(A, tmp);
    coords_get(tmp, &*xo, &*yo, &*zo);

    /* Go back to local coordinate system */
    tmp = coords_set(xi, yi, zi);
    tmp = rot_apply(A, tmp);
    coords_get(tmp, &xi, &yi, &zi);

    if (solid_angle) {
      /* Calculate vector from local point to remote random point */
      lx = *xo - lx;
      ly = *yo - ly;
      lz = *zo - lz;
      dist_p = sqrt(lx*lx + ly*ly + lz*lz);

      /* Adjust the 'solid angle' */
      /* 1/r^2 to the chosen point times cos(\theta) between the normal */
      /* vector of the target rectangle and direction vector of the chosen point. */
      cos_theta = (xi * lx + yi * ly + zi * lz) / (dist * dist_p);
      *solid_angle = width * height / (dist_p * dist_p);
      int counter;
      for (counter = 0; counter < order; counter++) {
        *solid_angle = *solid_angle * cos_theta;
      }
    }
  }
}
#endif
#pragma acc routine seq nohost
void randvec_target_rect_real_gpu(double *xo, double *yo, double *zo, double *solid_angle,
        double xi, double yi, double zi,
        double width, double height, Rotation A,
        double lx, double ly, double lz, int order,
        _class_particle* _particle)
{
  double dx, dy, dist, dist_p, nx, ny, nz, mx, my, mz, n_norm, m_norm;
  double cos_theta;
  Coords tmp;
  Rotation Ainverse;

  rot_transpose(A, Ainverse);

  if(height == 0.0 || width == 0.0)
  {
    randvec_target_circle(xo, yo, zo, solid_angle,
               xi, yi, zi, 0);
    return;
  }
  else
  {
    /* Now choose point uniformly on rectangle within width x height */
    dx = width*randpm1()/2.0;
    dy = height*randpm1()/2.0;

    /* Determine distance to target plane*/
    dist = sqrt(xi*xi + yi*yi + zi*zi);
    /* Go to global coordinate system */

    tmp = coords_set(xi, yi, zi);
    tmp = rot_apply(Ainverse, tmp);
    coords_get(tmp, &xi, &yi, &zi);

    /* Determine vector normal to trajectory axis (z) and gravity [0 1 0] */
    vec_prod(nx, ny, nz, xi, yi, zi, 0, 1, 0);

    /* This now defines the x-axis, normalize: */
    n_norm=sqrt(nx*nx + ny*ny + nz*nz);
    nx = nx/n_norm;
    ny = ny/n_norm;
    nz = nz/n_norm;

    /* Now, determine our y-axis (vertical in many cases...) */
    vec_prod(mx, my, mz, xi, yi, zi, nx, ny, nz);
    m_norm=sqrt(mx*mx + my*my + mz*mz);
    mx = mx/m_norm;
    my = my/m_norm;
    mz = mz/m_norm;

    /* Our output, random vector can now be defined by linear combination: */

    *xo = xi + dx * nx + dy * mx;
    *yo = yi + dx * ny + dy * my;
    *zo = zi + dx * nz + dy * mz;

    /* Go back to local coordinate system */
    tmp = coords_set(*xo, *yo, *zo);
    tmp = rot_apply(A, tmp);
    coords_get(tmp, &*xo, &*yo, &*zo);

    /* Go back to local coordinate system */
    tmp = coords_set(xi, yi, zi);
    tmp = rot_apply(A, tmp);
    coords_get(tmp, &xi, &yi, &zi);

    if (solid_angle) {
      /* Calculate vector from local point to remote random point */
      lx = *xo - lx;
      ly = *yo - ly;
      lz = *zo - lz;
      dist_p = sqrt(lx*lx + ly*ly + lz*lz);

      /* Adjust the 'solid angle' */
      /* 1/r^2 to the chosen point times cos(\theta) between the normal */
      /* vector of the target rectangle and direction vector of the chosen point. */
      cos_theta = (xi * lx + yi * ly + zi * lz) / (dist * dist_p);
      *solid_angle = width * height / (dist_p * dist_p);
      int counter;
      for (counter = 0; counter < order; counter++) {
        *solid_angle = *solid_angle * cos_theta;
      }
    }
  }
}
/* randvec_target_rect_real */


/* SECTION: random numbers ================================================== */


/* "Mersenne Twister", by Makoto Matsumoto and Takuji Nishimura. */
/* See http://www.math.keio.ac.jp/~matumoto/emt.html for original source. */
/*
   A C-program for MT19937, with initialization improved 2002/1/26.
   Coded by Takuji Nishimura and Makoto Matsumoto.

   Before using, initialize the state by using mt_srandom(seed)
   or init_by_array(init_key, key_length).

   Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
   All rights reserved.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

     1. Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.

     2. Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.

     3. The names of its contributors may not be used to endorse or promote
        products derived from this software without specific prior written
        permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


   Any feedback is very welcome.
   http://www.math.keio.ac.jp/matumoto/emt.html
   email: matumoto@math.keio.ac.jp
*/
#include <stdio.h>
/* Period parameters */
#define N 624
#define M 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define LOWER_MASK 0x7fffffffUL /* least significant r bits */

unsigned long mt[N]; /* the array for the state vector  */
int mti=N+1; /* mti==N+1 means mt[N] is not initialized */

// initializes mt[N] with a seed
void mt_srandom(unsigned long s)
{
    mt[0]= s & 0xffffffffUL;
    for (mti=1; mti<N; mti++) {
        mt[mti] =
            (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti);
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array mt[].                        */
        /* 2002/01/09 modified by Makoto Matsumoto             */
        mt[mti] &= 0xffffffffUL;
        /* for >32 bit machines */
    }
}
/* Initialize by an array with array-length.
   Init_key is the array for initializing keys.
   key_length is its length. */
void init_by_array(unsigned long init_key[], unsigned long key_length)
{
    int i, j, k;
    mt_srandom(19650218UL);
    i=1; j=0;
    k = (N>key_length ? N : key_length);
    for (; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1664525UL))
          + init_key[j] + j; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++; j++;
        if (i>=N) { mt[0] = mt[N-1]; i=1; }
        if (j>=key_length) j=0;
    }
    for (k=N-1; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1566083941UL))
          - i; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++;
        if (i>=N) { mt[0] = mt[N-1]; i=1; }
    }

    mt[0] = 0x80000000UL; /* MSB is 1; assuring non-zero initial array */
}
/* generates a random number on [0,0xffffffff]-interval */
unsigned long mt_random(void)
{
    unsigned long y;
    unsigned long mag01[2]={0x0UL, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    if (mti >= N) { /* generate N words at one time */
        int kk;

        if (mti == N+1)   /* if mt_srandom() has not been called, */
            mt_srandom(5489UL); /* a default initial seed is used */

        for (kk=0;kk<N-M;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for (;kk<N-1;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

        mti = 0;
    }

    y = mt[mti++];

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
}
#undef N
#undef M
#undef MATRIX_A
#undef UPPER_MASK
#undef LOWER_MASK
/* End of "Mersenne Twister". */

// randnorm: generate a random number from normal law
double randnorm_cpu(void)
{
  double v1, v2, s;
  int phase = 0;
  double X, u1, u2;

  if(phase == 0)
  {
    do
    {
      u1 = rand01_cpu();
      u2 = rand01_cpu();
      v1 = 2*u1 - 1;
      v2 = 2*u2 - 1;
      s = v1*v1 + v2*v2;
    } while(s >= 1 || s == 0);

    X = v1*sqrt(-2*log(s)/s);
  }
  else
  {
    X = v2*sqrt(-2*log(s)/s);
  }

  phase = 1 - phase;
  return X;
}
// Generate a random number from -1 to 1 with triangle distribution
double randtriangle_cpu(void) {
	double randnum = rand01_cpu();
	if (randnum>0.5) return(1-sqrt(2*(randnum-0.5)));
	else return(sqrt(2*randnum)-1);
}
// Random number between 0.0 and 1.0 (including?)
double rand01_cpu() {
	double randnum;
	randnum = (double) mt_random();
	randnum /= (double) MC_RAND_MAX + 1;
	return randnum;
}
// Return a random number between 1 and -1
double randpm1_cpu() {
	double randnum;
	randnum = (double) mt_random();
	randnum /= ((double) MC_RAND_MAX + 1) / 2;
	randnum -= 1;
	return randnum;
}
// Return a random number between 0 and max.
double rand0max_cpu(double max) {
	double randnum;
	randnum = (double) mt_random();
	randnum /= ((double) MC_RAND_MAX + 1) / max;
	return randnum;
}
// Return a random number between min and max.
double randminmax_cpu(double min, double max) {
	return rand0max_cpu(max - min) + max;
}


/*
RNG for GPU specific routines below, all based on the native cuda rand01
equivalent, curand_uniform (see mccode-r.h).

The footprint of the _gpu variants need to include a state variable, since each
gpu thread needs to keep track of rng iteration independently, and this must
be handled explicitly.
*/
#ifdef USE_PGI
#pragma acc routine seq nohost
double randpm1_gpu(curandState_t* state) {
	return ((double) curand_uniform(state) - 0.5) * 2;
}
#pragma acc routine seq nohost
double rand0max_gpu(double max, curandState_t* state) {
	return ((double) curand_uniform(state)) * max;
}
#pragma acc routine seq nohost
double randminmax_gpu(double min, double max, curandState_t* state) {
  return ((double) curand_uniform(state)) * (max - min) + max;
}
#pragma acc routine seq nohost
double randtriangle_gpu(curandState_t* state) {
	double randnum = curand_uniform(state);
	if (randnum>0.5) return(1-sqrt(2*(randnum-0.5)));
	else return(sqrt(2*randnum)-1);
}
#endif


/**
 * Wrapper functions for fprintf / printf on GPU
 */

/* These would ideally be #pragma acc routine seq */
/*int printf_GPU(const char *format,...) {
  return 0;
  }*/

/* These would ideally be #pragma acc routine seq */
/*int fprintf_GPU(FILE *stream, const char *format, ...) {
  return 0;
  }*/


/* SECTION: main and signal handlers ======================================== */

/*******************************************************************************
* mchelp: displays instrument executable help with possible options
*******************************************************************************/
static void
mchelp(char *pgmname)
{
  int i;

  fprintf(stderr, "%s (%s) instrument simulation, generated with " MCCODE_STRING " (" MCCODE_DATE ")\n", instrument_name, instrument_source);
  fprintf(stderr, "Usage: %s [options] [parm=value ...]\n", pgmname);
  fprintf(stderr,
"Options are:\n"
"  -s SEED   --seed=SEED      Set random seed (must be != 0)\n"
"  -n COUNT  --ncount=COUNT   Set number of particles to simulate.\n"
"  -d DIR    --dir=DIR        Put all data files in directory DIR.\n"
"  -t        --trace          Enable trace of " MCCODE_PARTICLE "s through instrument.\n"
"  -g        --gravitation    Enable gravitation for all trajectories.\n"
"  --no-output-files          Do not write any data files.\n"
"  -h        --help           Show this help message.\n"
"  -i        --info           Detailed instrument information.\n"
"  --source                   Show the instrument code which was compiled.\n"
"  --format=FORMAT            Output data files using FORMAT="
   FLAVOR_UPPER
#ifdef USE_NEXUS
   " NEXUS"
#endif
"\n\n"
);
#ifdef USE_MPI
  fprintf(stderr,
  "This instrument has been compiled with MPI support.\n  Use 'mpirun %s [options] [parm=value ...]'.\n", pgmname);
#endif
  if(numipar > 0)
  {
    fprintf(stderr, "Instrument parameters are:\n");
    for(i = 0; i < numipar; i++)
      if (mcinputtable[i].val && strlen(mcinputtable[i].val))
        fprintf(stderr, "  %-16s(%s) [default='%s']\n", mcinputtable[i].name,
        (*mcinputtypes[mcinputtable[i].type].parminfo)(mcinputtable[i].name),
        mcinputtable[i].val);
      else
        fprintf(stderr, "  %-16s(%s)\n", mcinputtable[i].name,
        (*mcinputtypes[mcinputtable[i].type].parminfo)(mcinputtable[i].name));
  }

#ifndef NOSIGNALS
  fprintf(stderr, "Known signals are: "
#ifdef SIGUSR1
  "USR1 (status) "
#endif
#ifdef SIGUSR2
  "USR2 (save) "
#endif
#ifdef SIGBREAK
  "BREAK (save) "
#endif
#ifdef SIGTERM
  "TERM (save and exit)"
#endif
  "\n");
#endif /* !NOSIGNALS */
} /* mchelp */


/* mcshowhelp: show help and exit with 0 */
static void
mcshowhelp(char *pgmname)
{
  mchelp(pgmname);
  exit(0);
}

/* mcusage: display usage when error in input arguments and exit with 1 */
static void
mcusage(char *pgmname)
{
  fprintf(stderr, "Error: incorrect command line arguments\n");
  mchelp(pgmname);
  exit(1);
}

/* mcenabletrace: enable trace/mcdisplay or error if requires recompile */
static void
mcenabletrace(void)
{
 if(traceenabled)
  mcdotrace = 1;
 else
 {
   fprintf(stderr,
           "Error: trace not enabled (mcenabletrace)\n"
           "Please re-run the " MCCODE_NAME " compiler "
                   "with the --trace option, or rerun the\n"
           "C compiler with the MC_TRACE_ENABLED macro defined.\n");
   exit(1);
 }
}

/*******************************************************************************
* mcreadparams: request parameters from the prompt (or use default)
*******************************************************************************/
void
mcreadparams(void)
{
  int i,j,status;
  static char buf[CHAR_BUF_LENGTH];
  char *p;
  int len;

  MPI_MASTER(printf("Instrument parameters for %s (%s)\n",
                    instrument_name, instrument_source));

  for(i = 0; mcinputtable[i].name != 0; i++)
  {
    do
    {
      MPI_MASTER(
                 if (mcinputtable[i].val && strlen(mcinputtable[i].val))
                   printf("Set value of instrument parameter %s (%s) [default='%s']:\n",
                          mcinputtable[i].name,
                          (*mcinputtypes[mcinputtable[i].type].parminfo)
                          (mcinputtable[i].name), mcinputtable[i].val);
                 else
                   printf("Set value of instrument parameter %s (%s):\n",
                          mcinputtable[i].name,
                          (*mcinputtypes[mcinputtable[i].type].parminfo)
                          (mcinputtable[i].name));
                 fflush(stdout);
                 );
#ifdef USE_MPI
      if(mpi_node_rank == mpi_node_root)
        {
          p = fgets(buf, CHAR_BUF_LENGTH, stdin);
          if(p == NULL)
            {
              fprintf(stderr, "Error: empty input for paramater %s (mcreadparams)\n", mcinputtable[i].name);
              exit(1);
            }
        }
      else
        p = buf;
      MPI_Bcast(buf, CHAR_BUF_LENGTH, MPI_CHAR, mpi_node_root, MPI_COMM_WORLD);
#else /* !USE_MPI */
      p = fgets(buf, CHAR_BUF_LENGTH, stdin);
      if(p == NULL)
        {
          fprintf(stderr, "Error: empty input for paramater %s (mcreadparams)\n", mcinputtable[i].name);
          exit(1);
        }
#endif /* USE_MPI */
      len = strlen(buf);
      if (!len || (len == 1 && (buf[0] == '\n' || buf[0] == '\r')))
      {
        if (mcinputtable[i].val && strlen(mcinputtable[i].val)) {
          strncpy(buf, mcinputtable[i].val, CHAR_BUF_LENGTH);  /* use default value */
          len = strlen(buf);
        }
      }
      for(j = 0; j < 2; j++)
      {
        if(len > 0 && (buf[len - 1] == '\n' || buf[len - 1] == '\r'))
        {
          len--;
          buf[len] = '\0';
        }
      }

      status = (*mcinputtypes[mcinputtable[i].type].getparm)
                   (buf, mcinputtable[i].par);
      if(!status)
      {
        (*mcinputtypes[mcinputtable[i].type].error)(mcinputtable[i].name, buf);
        if (!mcinputtable[i].val || strlen(mcinputtable[i].val)) {
          fprintf(stderr, "       Change %s default value in instrument definition.\n", mcinputtable[i].name);
          exit(1);
        }
      }
    } while(!status);
  }
} /* mcreadparams */

/*******************************************************************************
* mcparseoptions: parse command line arguments (options, parameters)
*******************************************************************************/
void
mcparseoptions(int argc, char *argv[])
{
  int i, j;
  char *p;
  int paramset = 0, *paramsetarray;
  char *usedir=NULL;

  /* Add one to numipar to avoid allocating zero size memory block. */
  paramsetarray = (int*)malloc((numipar + 1)*sizeof(*paramsetarray));
  if(paramsetarray == NULL)
  {
    fprintf(stderr, "Error: insufficient memory (mcparseoptions)\n");
    exit(1);
  }
  for(j = 0; j < numipar; j++)
    {
      paramsetarray[j] = 0;
      if (mcinputtable[j].val != NULL && strlen(mcinputtable[j].val))
      {
        int  status;
        char buf[CHAR_BUF_LENGTH];
        strncpy(buf, mcinputtable[j].val, CHAR_BUF_LENGTH);
        status = (*mcinputtypes[mcinputtable[j].type].getparm)
                   (buf, mcinputtable[j].par);
        if(!status) fprintf(stderr, "Invalid '%s' default value %s in instrument definition (mcparseoptions)\n", mcinputtable[j].name, buf);
        else paramsetarray[j] = 1;
      } else {
        (*mcinputtypes[mcinputtable[j].type].getparm)
          (NULL, mcinputtable[j].par);
        paramsetarray[j] = 0;
      }
    }
  for(i = 1; i < argc; i++)
  {
    if(!strcmp("-s", argv[i]) && (i + 1) < argc)
      mcsetseed(argv[++i]);
    else if(!strncmp("-s", argv[i], 2))
      mcsetseed(&argv[i][2]);
    else if(!strcmp("--seed", argv[i]) && (i + 1) < argc)
      mcsetseed(argv[++i]);
    else if(!strncmp("--seed=", argv[i], 7))
      mcsetseed(&argv[i][7]);
    else if(!strcmp("-n", argv[i]) && (i + 1) < argc)
      mcsetn_arg(argv[++i]);
    else if(!strncmp("-n", argv[i], 2))
      mcsetn_arg(&argv[i][2]);
    else if(!strcmp("--ncount", argv[i]) && (i + 1) < argc)
      mcsetn_arg(argv[++i]);
    else if(!strncmp("--ncount=", argv[i], 9))
      mcsetn_arg(&argv[i][9]);
    else if(!strcmp("-d", argv[i]) && (i + 1) < argc)
      usedir=argv[++i];  /* will create directory after parsing all arguments (end of this function) */
    else if(!strncmp("-d", argv[i], 2))
      usedir=&argv[i][2];
    else if(!strcmp("--dir", argv[i]) && (i + 1) < argc)
      usedir=argv[++i];
    else if(!strncmp("--dir=", argv[i], 6))
      usedir=&argv[i][6];
    else if(!strcmp("-h", argv[i]))
      mcshowhelp(argv[0]);
    else if(!strcmp("--help", argv[i]) || !strcmp("--version", argv[i]))
      mcshowhelp(argv[0]);
    else if(!strcmp("-i", argv[i])) {
      mcformat=FLAVOR_UPPER;
      mcinfo();
    }
    else if(!strcmp("--info", argv[i]))
      mcinfo();
    else if(!strcmp("-t", argv[i]))
      mcenabletrace();
    else if(!strcmp("--trace", argv[i]) || !strcmp("--verbose", argv[i]))
      mcenabletrace();
    else if(!strcmp("--gravitation", argv[i]))
      mcgravitation = 1;
    else if(!strcmp("-g", argv[i]))
      mcgravitation = 1;
    else if(!strncmp("--format=", argv[i], 9)) {
      mcformat=&argv[i][9];
    }
    else if(!strcmp("--format", argv[i]) && (i + 1) < argc) {
      mcformat=argv[++i];
    }
    else if(!strcmp("--no-output-files", argv[i]))
      mcdisable_output_files = 1;
    else if(!strcmp("--source", argv[i])) {
      printf("/* Source code %s from %s: */\n"
        "/******************************************************************************/\n"
        "%s\n"
        "/******************************************************************************/\n"
        "/* End of source code %s from %s */\n",
        instrument_name, instrument_source, instrument_code,
        instrument_name, instrument_source);
      exit(1);
    }
    else if(argv[i][0] != '-' && (p = strchr(argv[i], '=')) != NULL)
    {
      *p++ = '\0';

      for(j = 0; j < numipar; j++)
        if(!strcmp(mcinputtable[j].name, argv[i]))
        {
          int status;
          status = (*mcinputtypes[mcinputtable[j].type].getparm)(p,
                        mcinputtable[j].par);
          if(!status || !strlen(p))
          {
            (*mcinputtypes[mcinputtable[j].type].error)
              (mcinputtable[j].name, p);
            exit(1);
          }
          paramsetarray[j] = 1;
          paramset = 1;
          break;
        }
      if(j == numipar)
      {                                /* Unrecognized parameter name */
        fprintf(stderr, "Error: unrecognized parameter %s (mcparseoptions)\n", argv[i]);
        exit(1);
      }
    }
    else if(argv[i][0] == '-') {
      fprintf(stderr, "Error: unrecognized option argument %s (mcparseoptions). Ignored.\n", argv[i++]);
    }
    else {
      fprintf(stderr, "Error: unrecognized argument %s (mcparseoptions). Aborting.\n", argv[i]);
      mcusage(argv[0]);
    }
  }
  if(!paramset)
    mcreadparams();                /* Prompt for parameters if not specified. */
  else
  {
    for(j = 0; j < numipar; j++)
      if(!paramsetarray[j])
      {
        fprintf(stderr, "Error: Instrument parameter %s left unset (mcparseoptions)\n",
                mcinputtable[j].name);
        exit(1);
      }
  }
  free(paramsetarray);
#ifdef USE_MPI
  if (mcdotrace) mpi_node_count=1; /* disable threading when in trace mode */
#endif
  if (usedir && strlen(usedir) && !mcdisable_output_files) mcuse_dir(usedir);
} /* mcparseoptions */

#ifndef NOSIGNALS
/*******************************************************************************
* sighandler: signal handler that makes simulation stop, and save results
*******************************************************************************/
void sighandler(int sig)
{
  /* MOD: E. Farhi, Sep 20th 2001: give more info */
  time_t t1, t0;
#define SIG_SAVE 0
#define SIG_TERM 1
#define SIG_STAT 2
#define SIG_ABRT 3

  printf("\n# " MCCODE_STRING ": [pid %i] Signal %i detected", getpid(), sig);
#ifdef USE_MPI
  printf(" [proc %i]", mpi_node_rank);
#endif
#if defined(SIGUSR1) && defined(SIGUSR2) && defined(SIGKILL)
  if (!strcmp(mcsig_message, "sighandler") && (sig != SIGUSR1) && (sig != SIGUSR2))
  {
    printf("\n# Fatal : unrecoverable loop ! Suicide (naughty boy).\n");
    kill(0, SIGKILL); /* kill myself if error occurs within sighandler: loops */
  }
#endif
  switch (sig) {
#ifdef SIGINT
    case SIGINT : printf(" SIGINT (interrupt from terminal, Ctrl-C)"); sig = SIG_TERM; break;
#endif
#ifdef SIGILL
    case SIGILL  : printf(" SIGILL (Illegal instruction)"); sig = SIG_ABRT; break;
#endif
#ifdef SIGFPE
    case SIGFPE  : printf(" SIGFPE (Math Error)"); sig = SIG_ABRT; break;
#endif
#ifdef SIGSEGV
    case SIGSEGV : printf(" SIGSEGV (Mem Error)"); sig = SIG_ABRT; break;
#endif
#ifdef SIGTERM
    case SIGTERM : printf(" SIGTERM (Termination)"); sig = SIG_TERM; break;
#endif
#ifdef SIGABRT
    case SIGABRT : printf(" SIGABRT (Abort)"); sig = SIG_ABRT; break;
#endif
#ifdef SIGQUIT
    case SIGQUIT : printf(" SIGQUIT (Quit from terminal)"); sig = SIG_TERM; break;
#endif
#ifdef SIGTRAP
    case SIGTRAP : printf(" SIGTRAP (Trace trap)"); sig = SIG_ABRT; break;
#endif
#ifdef SIGPIPE
    case SIGPIPE : printf(" SIGPIPE (Broken pipe)"); sig = SIG_ABRT; break;
#endif
#ifdef SIGUSR1
    case SIGUSR1 : printf(" SIGUSR1 (Display info)"); sig = SIG_STAT; break;
#endif
#ifdef SIGUSR2
    case SIGUSR2 : printf(" SIGUSR2 (Save simulation)"); sig = SIG_SAVE; break;
#endif
#ifdef SIGHUP
    case SIGHUP  : printf(" SIGHUP (Hangup/update)"); sig = SIG_SAVE; break;
#endif
#ifdef SIGBUS
    case SIGBUS  : printf(" SIGBUS (Bus error)"); sig = SIG_ABRT; break;
#endif
#ifdef SIGURG
    case SIGURG  : printf(" SIGURG (Urgent socket condition)"); sig = SIG_ABRT; break;
#endif
#ifdef SIGBREAK
    case SIGBREAK: printf(" SIGBREAK (Break signal, Ctrl-Break)"); sig = SIG_SAVE; break;
#endif
    default : printf(" (look at signal list for signification)"); sig = SIG_ABRT; break;
  }
  printf("\n");
  printf("# Simulation: %s (%s) \n", instrument_name, instrument_source);
  printf("# Breakpoint: %s ", mcsig_message);
  if (strstr(mcsig_message, "Save") && (sig == SIG_SAVE))
    sig = SIG_STAT;
  SIG_MESSAGE("sighandler");
  if (mcget_ncount() == 0)
    printf("(0 %%)\n" );
  else
  {
    printf("%.2f %% (%10.1f/%10.1f)\n", 100.0*mcget_run_num()/mcget_ncount(), 1.0*mcget_run_num(), 1.0*mcget_ncount());
  }
  t0 = (time_t)mcstartdate;
  t1 = time(NULL);
  printf("# Date:      %s", ctime(&t1));
  printf("# Started:   %s", ctime(&t0));

  if (sig == SIG_STAT)
  {
    printf("# " MCCODE_STRING ": Resuming simulation (continue)\n");
    fflush(stdout);
    return;
  }
  else
  if (sig == SIG_SAVE)
  {
    printf("# " MCCODE_STRING ": Saving data and resume simulation (continue)\n");
    save(NULL);
    fflush(stdout);
    return;
  }
  else
  if (sig == SIG_TERM)
  {
    printf("# " MCCODE_STRING ": Finishing simulation (save results and exit)\n");
    finally();
    exit(0);
  }
  else
  {
    fflush(stdout);
    perror("# Last I/O Error");
    printf("# " MCCODE_STRING ": Simulation stop (abort).\n");
// This portion of the signal handling only works on UNIX
#if defined(__unix__) || defined(__APPLE__)
    signal(sig, SIG_DFL); /* force to use default sighandler now */
    kill(getpid(), sig);  /* and trigger it with the current signal */
#endif
    exit(-1);
  }
#undef SIG_SAVE
#undef SIG_TERM
#undef SIG_STAT
#undef SIG_ABRT

} /* sighandler */
#endif /* !NOSIGNALS */

#ifdef NEUTRONICS
/*Main neutronics function steers the McStas calls, initializes parameters etc */
/* Only called in case NEUTRONICS = TRUE */
void neutronics_main_(float *inx, float *iny, float *inz, float *invx, float *invy, float *invz, float *intime, float *insx, float *insy, float *insz, float *inw, float *outx, float *outy, float *outz, float *outvx, float *outvy, float *outvz, float *outtime, float *outsx, float *outsy, float *outsz, float *outwgt)
{

  extern double mcnx, mcny, mcnz, mcnvx, mcnvy, mcnvz;
  extern double mcnt, mcnsx, mcnsy, mcnsz, mcnp;

  /* External code governs iteration - McStas is iterated once per call to neutronics_main. I.e. below counter must be initiancated for each call to neutronics_main*/
  mcrun_num=0;

  time_t t;
  t = (time_t)mcstartdate;
  mcstartdate = t;  /* set start date before parsing options and creating sim file */
  init();

  /* *** parse options *** */
  SIG_MESSAGE("[" __FILE__ "] main START");
  mcformat=getenv(FLAVOR_UPPER "_FORMAT") ?
           getenv(FLAVOR_UPPER "_FORMAT") : FLAVOR_UPPER;

  /* Set neutron state based on input from neutronics code */
  mcsetstate(*inx,*iny,*inz,*invx,*invy,*invz,*intime,*insx,*insy,*insz,*inw);

  /* main neutron event loop - runs only one iteration */

  //mcstas_raytrace(&mcncount); /* prior to McStas 1.12 */

  mcallowbackprop = 1; //avoid absorbtion from negative dt
  int argc=1;
  char *argv[0];
  int dummy = mccode_main(argc, argv);

  *outx =  mcnx;
  *outy =  mcny;
  *outz =  mcnz;
  *outvx =  mcnvx;
  *outvy =  mcnvy;
  *outvz =  mcnvz;
  *outtime =  mcnt;
  *outsx =  mcnsx;
  *outsy =  mcnsy;
  *outsz =  mcnsz;
  *outwgt =  mcnp;

  return;
} /* neutronics_main */

#endif /*NEUTRONICS*/

#endif /* !MCCODE_H */
/* End of file "mccode-r.c". */
/* End of file "mccode-r.c". */

/* embedding file "mcstas-r.c" */

/*******************************************************************************
*
* McStas, neutron ray-tracing package
*         Copyright (C) 1997-2009, All rights reserved
*         Risoe National Laboratory, Roskilde, Denmark
*         Institut Laue Langevin, Grenoble, France
*
* Runtime: share/mcstas-r.c
*
* %Identification
* Written by: KN
* Date:    Aug 29, 1997
* Release: McStas X.Y
* Version: $Revision$
*
* Runtime system for McStas.
* Embedded within instrument in runtime mode.
*
* Usage: Automatically embbeded in the c code whenever required.
*
* $Id$
*
*******************************************************************************/

#ifndef MCSTAS_R_H
#include "mcstas-r.h"
#endif
#ifdef DANSE
#include "mcstas-globals.h"
#endif

/*******************************************************************************
* The I/O format definitions and functions
*******************************************************************************/

/*the magnet stack*/
#ifdef MC_POL_COMPAT
void (*mcMagnetPrecession) (double, double, double, double, double, double,
    double, double*, double*, double*, double, Coords, Rotation)=NULL;
Coords   mcMagnetPos;
Rotation mcMagnetRot;
double*  mcMagnetData                = NULL;
/* mcMagneticField(x, y, z, t, Bx, By, Bz) */
int (*mcMagneticField) (double, double, double, double,
    double*, double*, double*, void *) = NULL;
#endif

#ifndef MCSTAS_H

/*******************************************************************************
* mcsetstate: transfer parameters into global McStas variables
*******************************************************************************/
#pragma acc routine seq
_class_particle mcsetstate(double x, double y, double z, double vx, double vy, double vz,
           double t, double sx, double sy, double sz, double p)
{
  _class_particle mcneutron;

  mcneutron.x  = x;
  mcneutron.y  = y;
  mcneutron.z  = z;
  mcneutron.vx = vx;
  mcneutron.vy = vy;
  mcneutron.vz = vz;
  mcneutron.t  = t;
  mcneutron.sx = sx;
  mcneutron.sy = sy;
  mcneutron.sz = sz;
  mcneutron.p  = p;
  mcneutron._uid       = 0;
  mcneutron._index     = 1;
  mcneutron._absorbed  = 0;
  mcneutron._restore   = 0;
  mcneutron._scattered = 0;

  return(mcneutron);
} /* mcsetstate */

/*******************************************************************************
* mcgetstate: get neutron parameters from particle structure
*******************************************************************************/
#pragma acc routine seq
_class_particle mcgetstate(_class_particle mcneutron, double *x, double *y, double *z,
               double *vx, double *vy, double *vz, double *t,
               double *sx, double *sy, double *sz, double *p)
{
  *x  =  mcneutron.x;
  *y  =  mcneutron.y;
  *z  =  mcneutron.z;
  *vx =  mcneutron.vx;
  *vy =  mcneutron.vy;
  *vz =  mcneutron.vz;
  *t  =  mcneutron.t;
  *sx =  mcneutron.sx;
  *sy =  mcneutron.sy;
  *sz =  mcneutron.sz;
  *p  =  mcneutron.p;

  return(mcneutron);
} /* mcgetstate */


/*******************************************************************************
* mcgenstate: set default neutron parameters
*******************************************************************************/
#pragma acc routine seq
_class_particle mcgenstate(void)
{
  return(mcsetstate(0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1));
  /* old initialisation: mcsetstate(0, 0, 0, 0, 0, 1, 0, sx=0, sy=1, sz=0, 1); */
}

/*******************************************************************************
* mccoordschanges: old style rotation routine rot -> (x y z) ,(vx vy vz),(sx,sy,sz)
*******************************************************************************/
void
mccoordschanges(Coords a, Rotation t, double *x, double *y, double *z,
               double *vx, double *vy, double *vz, double *sx, double *sy, double *sz)
{
  Coords b, c;

  b.x = *x;
  b.y = *y;
  b.z = *z;
  c = rot_apply(t, b);
  b = coords_add(c, a);
  *x = b.x;
  *y = b.y;
  *z = b.z;

  if ( (vz && vy  && vx) && (*vz != 0.0 || *vx != 0.0 || *vy != 0.0) )
    mccoordschange_polarisation(t, vx, vy, vz);

  if ( (sz && sy  && sx) && (*sz != 0.0 || *sx != 0.0 || *sy != 0.0) )
    mccoordschange_polarisation(t, sx, sy, sz);

}

/*******************************************************************************
* mcrestore_neutron: restores neutron coodinates from global array
*******************************************************************************/
#pragma acc routine seq
_class_particle
mcrestore_neutron(MCNUM *s, int index, double *x, double *y, double *z,
               double *vx, double *vy, double *vz, double *t,
               double *sx, double *sy, double *sz, double *p)
{
    double *dptr = &s[11*index];
    *x  =  *dptr++;
    *y  =  *dptr++;
    *z  =  *dptr++;
    *vx =  *dptr++;
    *vy =  *dptr++;
    *vz =  *dptr++;
    *t  =  *dptr++;
    *sx =  *dptr++;
    *sy =  *dptr++;
    *sz =  *dptr++;
    *p  =  *dptr;

    return mcsetstate(*x, *y, *z, *vx, *vy, *vz, *t, *sx, *sy, *sz, *p);
} /* mcrestore_neutron */

/* intersection routines ==================================================== */

/*******************************************************************************
* inside_rectangle: Check if (x,y) is inside rectangle (xwidth, yheight)
* return 0 if outside and 1 if inside
*******************************************************************************/
int inside_rectangle(double x, double y, double xwidth, double yheight)
{
  if (x>-xwidth/2 && x<xwidth/2 && y>-yheight/2 && y<yheight/2)
    return 1;
  else
    return 0;
}

/*******************************************************************************
 * box_intersect: compute time intersection with a box
 * returns 0 when no intersection is found
 *      or 1 in case of intersection with resulting times dt_in and dt_out
 * This function written by Stine Nyborg, 1999.
 *******************************************************************************/
#pragma acc routine sequential
int box_intersect(double *dt_in, double *dt_out,
                  double x, double y, double z,
                  double vx, double vy, double vz,
                  double dx, double dy, double dz)
{
  double x_in, y_in, z_in, tt, t[6], a, b;
  int i, count, s;

      /* Calculate intersection time for each of the six box surface planes
       *  If the box surface plane is not hit, the result is zero.*/

  if(vx != 0)
   {
    tt = -(dx/2 + x)/vx;
    y_in = y + tt*vy;
    z_in = z + tt*vz;
    if( y_in > -dy/2 && y_in < dy/2 && z_in > -dz/2 && z_in < dz/2)
      t[0] = tt;
    else
      t[0] = 0;

    tt = (dx/2 - x)/vx;
    y_in = y + tt*vy;
    z_in = z + tt*vz;
    if( y_in > -dy/2 && y_in < dy/2 && z_in > -dz/2 && z_in < dz/2)
      t[1] = tt;
    else
      t[1] = 0;
   }
  else
    t[0] = t[1] = 0;

  if(vy != 0)
   {
    tt = -(dy/2 + y)/vy;
    x_in = x + tt*vx;
    z_in = z + tt*vz;
    if( x_in > -dx/2 && x_in < dx/2 && z_in > -dz/2 && z_in < dz/2)
      t[2] = tt;
    else
      t[2] = 0;

    tt = (dy/2 - y)/vy;
    x_in = x + tt*vx;
    z_in = z + tt*vz;
    if( x_in > -dx/2 && x_in < dx/2 && z_in > -dz/2 && z_in < dz/2)
      t[3] = tt;
    else
      t[3] = 0;
   }
  else
    t[2] = t[3] = 0;

  if(vz != 0)
   {
    tt = -(dz/2 + z)/vz;
    x_in = x + tt*vx;
    y_in = y + tt*vy;
    if( x_in > -dx/2 && x_in < dx/2 && y_in > -dy/2 && y_in < dy/2)
      t[4] = tt;
    else
      t[4] = 0;

    tt = (dz/2 - z)/vz;
    x_in = x + tt*vx;
    y_in = y + tt*vy;
    if( x_in > -dx/2 && x_in < dx/2 && y_in > -dy/2 && y_in < dy/2)
      t[5] = tt;
    else
      t[5] = 0;
   }
  else
    t[4] = t[5] = 0;

  /* The intersection is evaluated and *dt_in and *dt_out are assigned */

  a = b = s = 0;
  count = 0;

  for( i = 0; i < 6; i = i + 1 )
    if( t[i] == 0 )
      s = s+1;
    else if( count == 0 )
    {
      a = t[i];
      count = 1;
    }
    else
    {
      b = t[i];
      count = 2;
    }

  if ( a == 0 && b == 0 )
    return 0;
  else if( a < b )
  {
    *dt_in = a;
    *dt_out = b;
    return 1;
  }
  else
  {
    *dt_in = b;
    *dt_out = a;
    return 1;
  }

} /* box_intersect */

/*******************************************************************************
 * cylinder_intersect: compute intersection with a cylinder
 * returns 0 when no intersection is found
 *      or 2/4/8/16 bits depending on intersection,
 *     and resulting times t0 and t1
 * Written by: EM,NB,ABA 4.2.98
  *******************************************************************************/
#pragma acc routine sequential
int cylinder_intersect(double *t0, double *t1, double x, double y, double z,
                   double vx, double vy, double vz, double r, double h)
{
  double D, t_in, t_out, y_in, y_out;
  int ret=1;

  D = (2*vx*x + 2*vz*z)*(2*vx*x + 2*vz*z)
    - 4*(vx*vx + vz*vz)*(x*x + z*z - r*r);

  if (D>=0)
  {
    if (vz*vz + vx*vx) {
      t_in  = (-(2*vz*z + 2*vx*x) - sqrt(D))/(2*(vz*vz + vx*vx));
      t_out = (-(2*vz*z + 2*vx*x) + sqrt(D))/(2*(vz*vz + vx*vx));
    } else if (vy) { /* trajectory parallel to cylinder axis */
      t_in = (-h/2-y)/vy;
      t_out = (h/2-y)/vy;
      if (t_in>t_out){
        double tmp=t_in;
        t_in=t_out;t_out=tmp;
      }
    } else return 0;
    y_in = vy*t_in + y;
    y_out =vy*t_out + y;

    if ( (y_in > h/2 && y_out > h/2) || (y_in < -h/2 && y_out < -h/2) )
      return 0;
    else
    {
      if (y_in > h/2)
        { t_in = ((h/2)-y)/vy; ret += 2; }
      else if (y_in < -h/2)
        { t_in = ((-h/2)-y)/vy; ret += 4; }
      if (y_out > h/2)
        { t_out = ((h/2)-y)/vy; ret += 8; }
      else if (y_out < -h/2)
        { t_out = ((-h/2)-y)/vy; ret += 16; }
    }
    *t0 = t_in;
    *t1 = t_out;
    return ret;
  }
  else
  {
    *t0 = *t1 = 0;
    return 0;
  }
} /* cylinder_intersect */


/*******************************************************************************
 * sphere_intersect: Calculate intersection between a line and a sphere.
 * returns 0 when no intersection is found
 *      or 1 in case of intersection with resulting times t0 and t1
 *******************************************************************************/
#pragma acc routine sequential
int sphere_intersect(double *t0, double *t1, double x, double y, double z,
                 double vx, double vy, double vz, double r)
{
  double A, B, C, D, v;

  v = sqrt(vx*vx + vy*vy + vz*vz);
  A = v*v;
  B = 2*(x*vx + y*vy + z*vz);
  C = x*x + y*y + z*z - r*r;
  D = B*B - 4*A*C;
  if(D < 0)
    return 0;
  D = sqrt(D);
  *t0 = (-B - D) / (2*A);
  *t1 = (-B + D) / (2*A);
  return 1;
} /* sphere_intersect */

/*******************************************************************************
 * plane_intersect: Calculate intersection between a plane and a line.
 * returns 0 when no intersection is found (i.e. line is parallel to the plane)
 * returns 1 or -1 when intersection time is positive and negative respectively
 *******************************************************************************/
#pragma acc routine sequential
int plane_intersect(double *t, double x, double y, double z,
                 double vx, double vy, double vz, double nx, double ny, double nz, double wx, double wy, double wz)
{
  double s;
  if (fabs(s=scalar_prod(nx,ny,nz,vx,vy,vz))<FLT_EPSILON) return 0;
  *t = - scalar_prod(nx,ny,nz,x-wx,y-wy,z-wz)/s;
  if (*t<0) return -1;
  else return 1;
} /* plane_intersect */

#endif /* !MCSTAS_H */
/* End of file "mcstas-r.c". */


/* *****************************************************************************
* Start of instrument 'HZB_NEAT' generated code
***************************************************************************** */

#ifdef MC_TRACE_ENABLED
int traceenabled = 1;
#else
int traceenabled = 0;
#endif
#define MCSTAS "/usr/share/mcstas/3.0-dev/"
int   defaultmain         = 1;
char  instrument_name[]   = "HZB_NEAT";
char  instrument_source[] = "HZB_NEAT.instr";
char *instrument_exe      = NULL; /* will be set to argv[0] in main */
char  instrument_code[]   = "Instrument HZB_NEAT source code HZB_NEAT.instr is not embedded in this executable.\n  Use --source option when running McStas.\n";

int main(int argc, char *argv[]){return mccode_main(argc, argv);}

/* *****************************************************************************
* instrument 'HZB_NEAT' and components DECLARE
***************************************************************************** */

/* Instrument parameters: structure and a table for the initialisation
   (Used in e.g. inputparse and I/O function (e.g. detector_out) */

struct _struct_instrument_parameters {
  MCNUM _lambda;
  MCNUM _dlambda;
  MCNUM _rpm;
  char* _coh;
  char* _inc;
};
typedef struct _struct_instrument_parameters _class_instrument_parameters;

/* instrument SPLIT, GROUP and JUMP control logic */
struct instrument_logic_struct {
  long Split_Sample; /* this is the SPLIT counter decremented down to 0 */
  _class_particle Split_Sample_particle; /* this is the particle to duplicate */
};

struct _instrument_struct {
  char   _name[256]; /* the name of this instrument e.g. 'HZB_NEAT' */
/* Counters per component instance */
  double counter_AbsorbProp[38]; /* absorbed events in PROP routines */
  double counter_N[38], counter_P[38], counter_P2[38]; /* event counters after each component instance */
  _class_particle _trajectory[38]; /* current trajectory for STORE/RESTORE */
/* Components position table (absolute and relative coords) */
  Coords _position_relative[38]; /* positions of all components */
  Coords _position_absolute[38];
  _class_instrument_parameters _parameters; /* instrument parameters */
  struct instrument_logic_struct logic; /* instrument logic */
} _instrument_var;
struct _instrument_struct *instrument = & _instrument_var;
#pragma acc declare create ( _instrument_var )
#pragma acc declare create ( instrument )

int numipar = 5;
struct mcinputtable_struct mcinputtable[] = {
  "lambda", &(_instrument_var._parameters._lambda), instr_type_double, "6", 
  "dlambda", &(_instrument_var._parameters._dlambda), instr_type_double, "0.05", 
  "rpm", &(_instrument_var._parameters._rpm), instr_type_double, "10000", 
  "coh", &(_instrument_var._parameters._coh), instr_type_string, "Rb_liq_coh.sqw", 
  "inc", &(_instrument_var._parameters._inc), instr_type_string, "Rb_liq_inc.sqw", 
  NULL, NULL, instr_type_double, ""
};


/* ************************************************************************** */
/*             SHARE user declarations for all components                     */
/* ************************************************************************** */

/* Shared user declarations for all components types 'Source_gen'. */
/*******************************************************************************
*
* McStas, neutron ray-tracing package
*         Copyright 1997-2002, All rights reserved
*         Risoe National Laboratory, Roskilde, Denmark
*         Institut Laue Langevin, Grenoble, France
*
* Library: share/read_table-lib.h
*
* %Identification
* Written by: EF
* Date: Aug 28, 2002
* Origin: ILL
* Release: McStas 1.6
* Version: $Revision$
*
* This file is to be imported by components that may read data from table files
* It handles some shared functions.
*
* This library may be used directly as an external library. It has no dependency
*
* Usage: within SHARE
* %include "read_table-lib"
*
*******************************************************************************/

#ifndef READ_TABLE_LIB_H
#define READ_TABLE_LIB_H "$Revision$"

#define READ_TABLE_STEPTOL  0.04 /* tolerancy for constant step approx */

#ifndef MC_PATHSEP_C
#ifdef WIN32
#define MC_PATHSEP_C '\\'
#define MC_PATHSEP_S "\\"
#else  /* !WIN32 */
#ifdef MAC
#define MC_PATHSEP_C ':'
#define MC_PATHSEP_S ":"
#else  /* !MAC */
#define MC_PATHSEP_C '/'
#define MC_PATHSEP_S "/"
#endif /* !MAC */
#endif /* !WIN32 */
#endif /* !MC_PATHSEP_C */

#ifndef MCSTAS
#ifdef WIN32
#define MCSTAS "C:\\mcstas\\lib"
#else  /* !WIN32 */
#ifdef MAC
#define MCSTAS ":mcstas:lib" /* ToDo: What to put here? */
#else  /* !MAC */
#define MCSTAS "/usr/local/lib/mcstas"
#endif /* !MAC */
#endif /* !WIN32 */
#endif /* !MCSTAS */

#include <sys/stat.h>
#include <stdio.h>
#include <stdlib.h>

  typedef struct struct_table
  {
    char    filename[1024];
    long    filesize;
    char   *header;  /* text header, e.g. comments */
    double *data;    /* vector { x[0], y[0], ... x[n-1], y[n-1]... } */
    double  min_x;   /* min value of first column */
    double  max_x;   /* max value of first column */
    double  step_x;  /* minimal step value of first column */
    long    rows;    /* number of rows in matrix block */
    long    columns; /* number of columns in matrix block */

    long    begin;   /* start fseek index of block */
    long    end;     /* stop  fseek index of block */
    long    block_number;  /* block index. 0 is catenation of all */
    long    array_length;  /* number of elements in the t_Table array */
    char    monotonic;     /* true when 1st column/vector data is monotonic */
    char    constantstep;  /* true when 1st column/vector data has constant step */
    char    method[32];    /* interpolation method: nearest, linear */
  } t_Table;

typedef struct t_Read_table_file_item {
    int ref_count;
    t_Table *table_ref;
} t_Read_table_file_item;

typedef enum enum_Read_table_file_actions {STORE,FIND,GC}  t_Read_table_file_actions;

/* read_table-lib function prototypes */
/* ========================================================================= */

/* 'public' functions */
long     Table_Read              (t_Table *Table, char *File, long block_number);
long     Table_Read_Offset       (t_Table *Table, char *File, long block_number,
                                  long *offset, long max_lines);
long     Table_Read_Offset_Binary(t_Table *Table, char *File, char *Type,
                                  long *Offset, long Rows, long Columns);
long     Table_Rebin(t_Table *Table); /* rebin table with regular 1st column and interpolate all columns 2:end */
long     Table_Info (t_Table Table);
double   Table_Index(t_Table Table,   long i, long j); /* get indexed value */
double   Table_Value(t_Table Table, double X, long j); /* search X in 1st column and return interpolated value in j-column */
t_Table *Table_Read_Array(char *File, long *blocks);
void     Table_Free_Array(t_Table *Table);
long     Table_Info_Array(t_Table *Table);
int      Table_SetElement(t_Table *Table, long i, long j, double value);
long     Table_Init(t_Table *Table, long rows, long columns); /* create a Table */
double   Table_Value2d(t_Table Table, double X, double Y);    /* same as Table_Index with non-integer indices and 2d interpolation */
MCDETECTOR Table_Write(t_Table Table, char*file, char*xl, char*yl, 
           double x1, double x2, double y1, double y2); /* write Table to disk */
void * Table_File_List_Handler(t_Read_table_file_actions action, void *item, void *item_modifier);
t_Table *Table_File_List_find(char *name, int block, int offset);
int Table_File_List_gc(t_Table *tab);
void *Table_File_List_store(t_Table *tab);

#define Table_ParseHeader(header, ...) \
  Table_ParseHeader_backend(header,__VA_ARGS__,NULL);

char **Table_ParseHeader_backend(char *header, ...);

/* private functions */
void Table_Free(t_Table *Table);
long Table_Read_Handle(t_Table *Table, FILE *fid, long block_number, long max_lines, char *name);
static void Table_Stat(t_Table *Table);
double Table_Interp1d(double x, double x1, double y1, double x2, double y2);
double Table_Interp1d_nearest(double x, double x1, double y1, double x2, double y2);
double Table_Interp2d(double x, double y, double x1, double y1, double x2, double y2,
double z11, double z12, double z21, double z22);

#endif

/* end of read_table-lib.h */
/*******************************************************************************
*
* McStas, neutron ray-tracing package
*         Copyright (C) 1997-2009, All rights reserved
*         Risoe National Laboratory, Roskilde, Denmark
*         Institut Laue Langevin, Grenoble, France
*
* Library: share/read_table-lib.c
*
* %Identification
* Written by: EF
* Date: Aug 28, 2002
* Origin: ILL
* Release: McStas CVS_090504
* Version: $Revision$
*
* This file is to be imported by components that may read data from table files
* It handles some shared functions. Embedded within instrument in runtime mode.
*
* Usage: within SHARE
* %include "read_table-lib"
*
*******************************************************************************/

#ifndef READ_TABLE_LIB_H
#include "read_table-lib.h"
#endif


/*******************************************************************************
 * void *Table_File_List_Handler(action, item, item_modifier)
 *   ACTION: handle file entries in the read_table-lib file list. If a file is read - it is supposed to be
 *   stored in a list such that we can avoid reading the same file many times.
 *   input  action: FIND, STORE, GC. check if file exists in the list, store an item in the list, or check if it can be garbage collected.
 *   input item: depends on the action.
 *    FIND)  item is a filename, and item_modifier is the block number
 *    STORE) item is the Table to store - item_modifier is ignored
 *    GC)    item is the Table to check. If it has a ref_count >1 then this is simply decremented.
 *   return  depends on the action
 *    FIND)  return a reference to a table+ref_count item if found - NULL otherwise. I.e. NULL means the file has not been read before and must be read again.
 *    STORE) return NULL always
 *    GC)    return NULL if no garbage collection is needed, return an adress to the t_Table which should be garbage collected. 0x1 is returned if
 *           the item is not found in the list
*******************************************************************************/
void * Table_File_List_Handler(t_Read_table_file_actions action, void *item, void *item_modifier){

    /* logic here is Read_Table should include a call to FIND. If found the return value shoud just be used as
     * if the table had been read. If not found then read the table and STORE.
     * Table_Free should include a call to GC. If this returns non-NULL then we shoudl proceed with freeing the memory
     * associated with the table item - otherwise do nothing since there are more references that may need it.*/

    static t_Read_table_file_item read_table_file_list[1024];
    static int read_table_file_count=0;

    t_Read_table_file_item *tr;
    switch(action){
        case FIND:
            /*interpret data item as a filename, if it is found return a pointer to the table and increment refcount.
             * if not found return the item itself*/
            tr=read_table_file_list;
            while ( tr->table_ref!=NULL ){
                int i=*((int*) item_modifier);
                int j=*( ((int*) item_modifier)+1);
                if ( !strcmp(tr->table_ref->filename,(char *) item) &&
                        tr->table_ref->block_number==i && tr->table_ref->begin==j ){
                    tr->ref_count++;
                    return (void *) tr;
                }
                tr++;
            }
            return NULL;
        case STORE:
            /*find an available slot and store references to table there*/
            tr=&(read_table_file_list[read_table_file_count++]);
            tr->table_ref=(t_Table *)calloc(1,sizeof(t_Table));
            /*copy the contents of the table handle*/
            *(tr->table_ref)= *((t_Table *) item);
            tr->ref_count++;
            return NULL;
        case GC:
            /* Should this item be garbage collected (freed) - if so scratch the entry and return the address of the item -
             * else decrement ref_count and return NULL.
             * A non-NULL return expects the item to actually be freed afterwards.*/
            tr=read_table_file_list;
            while ( tr->table_ref!=NULL ){
                if ( tr->table_ref->data ==((t_Table *)item)->data &&
                        tr->table_ref->block_number == ((t_Table *)item)->block_number){
                    /*matching item found*/
                    if (tr->ref_count>1){
                        /*the item is found - no garbage collection needed*/
                        tr->ref_count--;
                        return NULL;
                    }else{
                        /* The item is found - move remaining list items up one slot,
                         * and return the table for garbage collection by caller*/
                        while (tr->table_ref!=NULL){
                            *tr=*(tr+1);
                            tr++;
                        }
                        read_table_file_count--;
                        return (t_Table *) item;
                    }
                }
                tr++;
            }
            return (void *)0x1 ;/*item not found*/
    }

}

/* Access functions to the handler*/

/********************************************
 * t_Table *Table_File_List_find(char *name, int block, int offset)
 * input name: filename to search for in the file list
 * input block: data block in the file as each file may contain more than 1 data block.
 * return a ref. to a table if it is found (you may use this pointer and skip reading the file), NULL otherwise (i.e. go ahead and read the file)
*********************************************/
t_Table *Table_File_List_find(char *name, int block, int offset){
    int vars[2]={block,offset};
    t_Read_table_file_item *item = Table_File_List_Handler(FIND,name, vars);
    if (item == NULL){
        return NULL;
    }else{
        return item->table_ref;
    }
}
/********************************************
 * int Table_File_List_gc(t_Table *tab)
 * input tab: the table to check for references.
 * return 0: no garbage collection needed
 *        1: Table's data and header (at least) should be freed.
*********************************************/
int Table_File_List_gc(t_Table *tab){
    void *rval=Table_File_List_Handler(GC,tab,0);
    if (rval==NULL) return 0;
    else return 1;
}


/*****************************************************************************
 * void *Table_File_List_store(t_Table *tab)
 * input tab: pointer to table to store.
 * return None.
*******************************************************************************/
void *Table_File_List_store(t_Table *tab){
    Table_File_List_Handler(STORE,tab,0);
}


/*******************************************************************************
* FILE *Open_File(char *name, char *Mode, char *path)
*   ACTION: search for a file and open it. Optionally return the opened path.
*   input   name:  file name from which table should be extracted
*           mode: "r", "w", "a" or any valid fopen mode
*           path:  NULL or a pointer to at least 1024 allocated chars
*   return  initialized file handle or NULL in case of error
*******************************************************************************/

  FILE *Open_File(char *File, const char *Mode, char *Path)
  {
    char path[1024];
    FILE *hfile = NULL;

    if (!File || File[0]=='\0')                     return(NULL);
    if (!strcmp(File,"NULL") || !strcmp(File,"0"))  return(NULL);

    /* search in current or full path */
    strncpy(path, File, 1024);
    hfile = fopen(path, Mode);
    if(!hfile)
    {
      char dir[1024];

      if (!hfile && instrument_source && strlen(instrument_source)) /* search in instrument source location */
      {
        char *path_pos   = NULL;
        /* extract path: searches for last file separator */
        path_pos    = strrchr(instrument_source, MC_PATHSEP_C);  /* last PATHSEP */
        if (path_pos) {
          long path_length = path_pos +1 - instrument_source;  /* from start to path+sep */
          if (path_length) {
            strncpy(dir, instrument_source, path_length);
            dir[path_length] = '\0';
            snprintf(path, 1024, "%s%c%s", dir, MC_PATHSEP_C, File);
            hfile = fopen(path, Mode);
          }
        }
      }
      if (!hfile && instrument_exe && strlen(instrument_exe)) /* search in PWD instrument executable location */
      {
        char *path_pos   = NULL;
        /* extract path: searches for last file separator */
        path_pos    = strrchr(instrument_exe, MC_PATHSEP_C);  /* last PATHSEP */
        if (path_pos) {
          long path_length = path_pos +1 - instrument_exe;  /* from start to path+sep */
          if (path_length) {
            strncpy(dir, instrument_exe, path_length);
            dir[path_length] = '\0';
            snprintf(path, 1024, "%s%c%s", dir, MC_PATHSEP_C, File);
            hfile = fopen(path, Mode);
          }
        }
      }
      if (!hfile) /* search in HOME or . */
      {
        strcpy(dir, getenv("HOME") ? getenv("HOME") : ".");
        snprintf(path, 1024, "%s%c%s", dir, MC_PATHSEP_C, File);
        hfile = fopen(path, Mode);
      }
      if (!hfile) /* search in MCSTAS/data */
      {
        strcpy(dir, getenv(FLAVOR_UPPER) ? getenv(FLAVOR_UPPER) : MCSTAS);
        snprintf(path, 1024, "%s%c%s%c%s", dir, MC_PATHSEP_C, "data", MC_PATHSEP_C, File);
        hfile = fopen(path, Mode);
      }
      if (!hfile) /* search in MVCSTAS/contrib */
      {
        strcpy(dir, getenv(FLAVOR_UPPER) ? getenv(FLAVOR_UPPER) : MCSTAS);
        snprintf(path, 1024, "%s%c%s%c%s", dir, MC_PATHSEP_C, "contrib", MC_PATHSEP_C, File);
        hfile = fopen(path, Mode);
      }
      if(!hfile)
      {
        fprintf(stderr, "Error: Could not open input file '%s' (Open_File)\n", File);
        return (NULL);
      }
    }
    if (Path) strncpy(Path, path, 1024);
    return(hfile);
  } /* end Open_File */

/*******************************************************************************
* long Read_Table(t_Table *Table, char *name, int block_number)
*   ACTION: read a single Table from a text file
*   input   Table: pointer to a t_Table structure
*           name:  file name from which table should be extracted
*           block_number: if the file does contain more than one
*                 data block, then indicates which one to get (from index 1)
*                 a 0 value means append/catenate all
*   return  initialized single Table t_Table structure containing data, header, ...
*           number of read elements (-1: error, 0:header only)
* The routine stores any line starting with '#', '%' and ';' into the header
* File is opened, read and closed
* Other lines are interpreted as numerical data, and stored.
* Data block should be a rectangular matrix or vector.
* Data block may be rebinned with Table_Rebin (also sort in ascending order)
*******************************************************************************/
  long Table_Read(t_Table *Table, char *File, long block_number)
  { /* reads all or a single data block from 'file' and returns a Table structure  */
    return(Table_Read_Offset(Table, File, block_number, NULL, 0));
  } /* end Table_Read */

/*******************************************************************************
* long Table_Read_Offset(t_Table *Table, char *name, int block_number, long *offset
*                        long max_rows)
*   ACTION: read a single Table from a text file, starting at offset
*     Same as Table_Read(..) except:
*   input   offset:    pointer to an offset (*offset should be 0 at start)
*           max_rows: max number of data rows to read from file (0 means all)
*   return  initialized single Table t_Table structure containing data, header, ...
*           number of read elements (-1: error, 0:header only)
*           updated *offset position (where end of reading occured)
*******************************************************************************/
  long Table_Read_Offset(t_Table *Table, char *File,
                         long block_number, long *offset,
                         long max_rows)
  { /* reads all/a data block in 'file' and returns a Table structure  */
    FILE *hfile;
    long  nelements=0;
    long  begin=0;
    long  filesize=0;
    char  name[1024];
    char  path[1024];
    struct stat stfile;

    /*Need to be able to store the pointer*/
    if (!Table) return(-1);

    //if (offset && *offset) snprintf(name, 1024, "%s@%li", File, *offset);
    //else
    strncpy(name, File, 1024);
    if(offset && *offset){
        begin=*offset;
    }
    /* Check if the table has already been read from file.
     * If so just reuse the table, if not (this is flagged by returning NULL
     * set up a new table and read the data into it */
    t_Table *tab_p= Table_File_List_find(name,block_number,begin);
    if ( tab_p!=NULL ){
        /*table was found in the Table_File_List*/
        // printf("Reusing input file '%s' (Table_Read_Offset)\n", name);
        *Table=*tab_p;
        return Table->rows*Table->columns;
    }

    /* open the file */
    hfile = Open_File(File, "r", path);
    if (!hfile) return(-1);
    else {
      MPI_MASTER(
      printf("Opening input file '%s' (Table_Read_Offset)\n", path);
      );
    }

    /* read file state */
    stat(path,&stfile); filesize = stfile.st_size;
    if (offset && *offset) fseek(hfile, *offset, SEEK_SET);
    begin     = ftell(hfile);

    Table_Init(Table, 0, 0);

    /* read file content and set the Table */
    nelements = Table_Read_Handle(Table, hfile, block_number, max_rows, name);
    Table->begin = begin;
    Table->end   = ftell(hfile);
    Table->filesize = (filesize>0 ? filesize : 0);
    Table_Stat(Table);

    Table_File_List_store(Table);

    if (offset) *offset=Table->end;
    fclose(hfile);
    return(nelements);

  } /* end Table_Read_Offset */

/*******************************************************************************
* long Table_Read_Offset_Binary(t_Table *Table, char *File, char *type,
*                               long *offset, long rows, long columns)
*   ACTION: read a single Table from a binary file, starting at offset
*     Same as Table_Read_Offset(..) except that it handles binary files.
*   input   type: may be "float"/NULL or "double"
*           offset: pointer to an offset (*offset should be 0 at start)
*           rows   : number of rows (0 means read all)
*           columns: number of columns
*   return  initialized single Table t_Table structure containing data, header, ...
*           number of read elements (-1: error, 0:header only)
*           updated *offset position (where end of reading occured)
*******************************************************************************/
  long Table_Read_Offset_Binary(t_Table *Table, char *File, char *type,
                                long *offset, long rows, long columns)
  { /* reads all/a data block in binary 'file' and returns a Table structure  */
    long    nelements, sizeofelement;
    long    filesize;
    FILE   *hfile;
    char    path[1024];
    struct stat stfile;
    double *data;
    long    i;
    long    begin;

    if (!Table) return(-1);

    Table_Init(Table, 0, 0);

    /* open the file */
    hfile = Open_File(File, "r", path);
    if (!hfile) return(-1);
    else {
      MPI_MASTER(
      printf("Opening input file '%s' (Table_Read, Binary)\n", path);
      );
    }

    /* read file state */
    stat(File,&stfile);
    filesize = stfile.st_size;
    Table->filesize=filesize;

    /* read file content */
    if (type && !strcmp(type,"double")) sizeofelement = sizeof(double);
    else  sizeofelement = sizeof(float);
    if (offset && *offset) fseek(hfile, *offset, SEEK_SET);
    begin     = ftell(hfile);
    if (rows && filesize > sizeofelement*columns*rows)
      nelements = columns*rows;
    else nelements = (long)(filesize/sizeofelement);
    if (!nelements || filesize <= *offset) return(0);
    data    = (double*)malloc(nelements*sizeofelement);
    if (!data) {
      fprintf(stderr,"Error: allocating %ld elements for %s file '%s'. Too big (Table_Read_Offset_Binary).\n", nelements, type, File);
      exit(-1);
    }
    nelements = fread(data, sizeofelement, nelements, hfile);

    if (!data || !nelements)
    {
      fprintf(stderr,"Error: reading %ld elements from %s file '%s' (Table_Read_Offset_Binary)\n", nelements, type, File);
      exit(-1);
    }
    Table->begin   = begin;
    Table->end     = ftell(hfile);
    if (offset) *offset=Table->end;
    fclose(hfile);
    data = (double*)realloc(data, (double)nelements*sizeofelement);
    /* copy file data into Table */
    if (type && !strcmp(type,"double")) Table->data = data;
    else {
      float  *s;
      double *dataf;
      s     = (float*)data;
      dataf = (double*)malloc(sizeof(double)*nelements);
      for (i=0; i<nelements; i++)
        dataf[i]=s[i];
      free(data);
      Table->data = dataf;
    }
    strncpy(Table->filename, File, 1024);
    Table->rows    = nelements/columns;
    Table->columns = columns;
    Table->array_length = 1;
    Table->block_number = 1;

    Table_Stat(Table);

    return(nelements);
  } /* end Table_Read_Offset_Binary */

/*******************************************************************************
* long Table_Read_Handle(t_Table *Table, FILE *fid, int block_number, long max_rows, char *name)
*   ACTION: read a single Table from a text file handle (private)
*   input   Table:pointer to a t_Table structure
*           fid:  pointer to FILE handle
*           block_number: if the file does contain more than one
*                 data block, then indicates which one to get (from index 1)
*                 a 0 value means append/catenate all
*           max_rows: if non 0, only reads that number of lines
*   return  initialized single Table t_Table structure containing data, header, ...
*           modified Table t_Table structure containing data, header, ...
*           number of read elements (-1: error, 0:header only)
* The routine stores any line starting with '#', '%' and ';' into the header
* Other lines are interpreted as numerical data, and stored.
* Data block should be a rectangular matrix or vector.
* Data block may be rebined with Table_Rebin (also sort in ascending order)
*******************************************************************************/
  long Table_Read_Handle(t_Table *Table, FILE *hfile,
                         long block_number, long max_rows, char *name)
  { /* reads all/a data block from 'file' handle and returns a Table structure  */
    double *Data;
    char *Header              = NULL;
    long  malloc_size         = CHAR_BUF_LENGTH;
    long  malloc_size_h       = 4096;
    long  Rows = 0,   Columns = 0;
    long  count_in_array      = 0;
    long  count_in_header     = 0;
    long  block_Current_index = 0;
    char  flag_End_row_loop   = 0;

    if (!Table) return(-1);
    Table_Init(Table, 0, 0);
    if (name && name[0]!='\0') strncpy(Table->filename, name, 1024);

    if(!hfile) {
       fprintf(stderr, "Error: File handle is NULL (Table_Read_Handle).\n");
       return (-1);
    }
    Header = (char*)  calloc(malloc_size_h, sizeof(char));
    Data   = (double*)calloc(malloc_size,   sizeof(double));
    if ((Header == NULL) || (Data == NULL)) {
       fprintf(stderr, "Error: Could not allocate Table and Header (Table_Read_Handle).\n");
       return (-1);
    }

    int flag_In_array = 0;
    do { /* while (!flag_End_row_loop) */
      char  line[1024*CHAR_BUF_LENGTH];
      long  back_pos=0;   /* ftell start of line */

      back_pos = ftell(hfile);
      if (fgets(line, 1024*CHAR_BUF_LENGTH, hfile) != NULL) { /* analyse line */
        /* first skip blank and tabulation characters */
        int i = strspn(line, " \t");

        /* handle comments: stored in header */
        if (NULL != strchr("#%;/", line[i]))
        { /* line is a comment */
          count_in_header += strlen(line);
          if (count_in_header >= malloc_size_h) {
            /* if succeed and in array : add (and realloc if necessary) */
            malloc_size_h = count_in_header+4096;
            Header        = (char*)realloc(Header, malloc_size_h*sizeof(char));
          }
          strncat(Header, line, 4096);
          flag_In_array=0;
          /* exit line and file if passed desired block */
          if (block_number > 0 && block_number == block_Current_index) {
            flag_End_row_loop = 1;
          }

          /* Continue with next line */
          continue;
        }

        /* get the number of columns splitting line with strtok */
        char  *lexeme;
        char  flag_End_Line = 0;
        long  block_Num_Columns = 0;
        const char seps[] = " ,;\t\n\r";

        lexeme = strtok(line, seps);
        while (!flag_End_Line) {
          if ((lexeme != NULL) && (lexeme[0] != '\0')) {
            /* reading line: the token is not empty */
            double X;
            int    count=1;
            /* test if we have 'NaN','Inf' */
            if (!strncasecmp(lexeme,"NaN",3))
              X = 0;
            else if (!strncasecmp(lexeme,"Inf",3) || !strncasecmp(lexeme,"+Inf",4))
              X = FLT_MAX;
            else if (!strncasecmp(lexeme,"-Inf",4))
              X = -FLT_MAX;
            else
              count = sscanf(lexeme,"%lg",&X);
            if (count == 1) {
              /* reading line: the token is a number in the line */
              if (!flag_In_array) {
                /* reading num: not already in a block: starts a new data block */
                block_Current_index++;
                flag_In_array    = 1;
                block_Num_Columns= 0;
                if (block_number > 0) {
                  /* initialise a new data block */
                  Rows = 0;
                  count_in_array = 0;
                } /* else append */
              }
              /* reading num: all blocks or selected block */
              if (flag_In_array && (block_number == 0 ||
                  block_number == block_Current_index)) {
                /* starting block: already the desired number of rows ? */
                if (block_Num_Columns == 0 &&
                    max_rows > 0 && Rows >= max_rows) {
                  flag_End_Line      = 1;
                  flag_End_row_loop  = 1;
                  flag_In_array      = 0;
                  /* reposition to begining of line (ignore line) */
                  fseek(hfile, back_pos, SEEK_SET);
                } else { /* store into data array */
                  if (count_in_array >= malloc_size) {
                    /* realloc data buffer if necessary */
                    malloc_size = count_in_array+CHAR_BUF_LENGTH;
                    Data = (double*) realloc(Data, malloc_size*sizeof(double));
                    if (Data == NULL) {
                      fprintf(stderr, "Error: Can not re-allocate memory %li (Table_Read_Handle).\n",
                              malloc_size*sizeof(double));
                      return (-1);
                    }
                  }
                  if (0 == block_Num_Columns) Rows++;
                  Data[count_in_array] = X;
                  count_in_array++;
                  block_Num_Columns++;
                }
              } /* reading num: end if flag_In_array */
            } /* end reading num: end if sscanf lexeme -> numerical */
            else {
              /* reading line: the token is not numerical in that line. end block */
              if (block_Current_index == block_number) {
                flag_End_Line = 1;
                flag_End_row_loop = 1;
              } else {
                flag_In_array = 0;
                flag_End_Line = 1;
              }
            }
          }
          else {
            /* no more tokens in line */
            flag_End_Line = 1;
            if (block_Num_Columns > 0) Columns = block_Num_Columns;
          }

          // parse next token
          lexeme = strtok(NULL, seps);

        } /* while (!flag_End_Line) */
      } /* end: if fgets */
      else flag_End_row_loop = 1; /* else fgets : end of file */

    } while (!flag_End_row_loop); /* end while flag_End_row_loop */

    Table->block_number = block_number;
    Table->array_length = 1;

    // shrink header to actual size (plus terminating 0-byte)
    if (count_in_header) {
      Header = (char*)realloc(Header, count_in_header*sizeof(char) + 1);
    }
    Table->header = Header;

    if (count_in_array*Rows*Columns == 0)
    {
      Table->rows         = 0;
      Table->columns      = 0;
      free(Data);
      return (0);
    }
    if (Rows * Columns != count_in_array)
    {
      fprintf(stderr, "Warning: Read_Table :%s %s Data has %li values that should be %li x %li\n",
        (Table->filename ? Table->filename : ""),
        (!block_number ? " catenated" : ""),
        count_in_array, Rows, Columns);
      Columns = count_in_array; Rows = 1;
    }
    Data     = (double*)realloc(Data, count_in_array*sizeof(double));
    Table->data         = Data;
    Table->rows         = Rows;
    Table->columns      = Columns;

    return (count_in_array);

  } /* end Table_Read_Handle */

/*******************************************************************************
* long Table_Rebin(t_Table *Table)
*   ACTION: rebin a single Table, sorting 1st column in ascending order
*   input   Table: single table containing data.
*                  The data block is reallocated in this process
*   return  updated Table with increasing, evenly spaced first column (index 0)
*           number of data elements (-1: error, 0:empty data)
*******************************************************************************/
  long Table_Rebin(t_Table *Table)
  {
    double new_step=0;
    long   i;
    /* performs linear interpolation on X axis (0-th column) */

    if (!Table) return(-1);
    if (!Table->data
    || Table->rows*Table->columns == 0 || !Table->step_x)
      return(0);
    Table_Stat(Table); /* recompute statitstics and minimal step */
    new_step = Table->step_x; /* minimal step in 1st column */

    if (!(Table->constantstep)) /* not already evenly spaced */
    {
      long Length_Table;
      double *New_Table;

      Length_Table = ceil(fabs(Table->max_x - Table->min_x)/new_step)+1;
      New_Table    = (double*)malloc(Length_Table*Table->columns*sizeof(double));

      for (i=0; i < Length_Table; i++)
      {
        long   j;
        double X;
        X = Table->min_x + i*new_step;
        New_Table[i*Table->columns] = X;
        for (j=1; j < Table->columns; j++)
          New_Table[i*Table->columns+j]
                = Table_Value(*Table, X, j);
      } /* end for i */

      Table->rows = Length_Table;
      Table->step_x = new_step;
      Table->max_x = Table->min_x + (Length_Table-1)*new_step;
      /*max might not be the same anymore
       * Use Length_Table -1 since the first and laset rows are the limits of the defined interval.*/
      free(Table->data);
      Table->data = New_Table;
      Table->constantstep=1;
    } /* end else (!constantstep) */
    return (Table->rows*Table->columns);
  } /* end Table_Rebin */

/*******************************************************************************
* double Table_Index(t_Table Table, long i, long j)
*   ACTION: read an element [i,j] of a single Table
*   input   Table: table containing data
*           i : index of row      (0:Rows-1)
*           j : index of column   (0:Columns-1)
*   return  Value = data[i][j]
* Returns Value from the i-th row, j-th column of Table
* Tests are performed on indexes i,j to avoid errors
*******************************************************************************/

#ifndef MIN
#define MIN(a, b)  (((a) < (b)) ? (a) : (b))
#endif
#ifndef MAX
#define MAX(a, b)  (((a) > (b)) ? (a) : (b))
#endif

double Table_Index(t_Table Table, long i, long j)
{
  long AbsIndex;

  if (Table.rows == 1 || Table.columns == 1) {
    /* vector */
    j = MIN(MAX(0, i+j), Table.columns*Table.rows - 1);
    i = 0;
  } else {
    /* matrix */
    i = MIN(MAX(0, i), Table.rows - 1);
    j = MIN(MAX(0, j), Table.columns - 1);
  }

  /* handle vectors specifically */
  AbsIndex = i*(Table.columns)+j;

  if (Table.data != NULL)
    return (Table.data[AbsIndex]);
  else
    return 0;
} /* end Table_Index */

/*******************************************************************************
* void Table_SetElement(t_Table *Table, long i, long j, double value)
*   ACTION: set an element [i,j] of a single Table
*   input   Table: table containing data
*           i : index of row      (0:Rows-1)
*           j : index of column   (0:Columns-1)
*           value = data[i][j]
* Returns 0 in case of error
* Tests are performed on indexes i,j to avoid errors
*******************************************************************************/
int Table_SetElement(t_Table *Table, long i, long j,
                     double value)
{
  long AbsIndex;

  if (Table->rows == 1 || Table->columns == 1) {
    /* vector */
    j = MIN(MAX(0, i+j), Table->columns*Table->rows - 1); i=0;
  } else {
    /* matrix */
    i = MIN(MAX(0, i), Table->rows - 1);
    j = MIN(MAX(0, j), Table->columns - 1);
  }

  AbsIndex = i*(Table->columns)+j;
  if (Table->data != NULL) {
    Table->data[AbsIndex] = value;
    return 1;
  }

  return 0;
} /* end Table_SetElement */

/*******************************************************************************
* double Table_Value(t_Table Table, double X, long j)
*   ACTION: read column [j] of a single Table at row which 1st column is X
*   input   Table: table containing data.
*           X : data value in the first column (index 0)
*           j : index of column from which is extracted the Value (0:Columns-1)
*   return  Value = data[index for X][j] with linear interpolation
* Returns Value from the j-th column of Table corresponding to the
* X value for the 1st column (index 0)
* Tests are performed (within Table_Index) on indexes i,j to avoid errors
* NOTE: data should rather be monotonic, and evenly sampled.
*******************************************************************************/
double Table_Value(t_Table Table, double X, long j)
{
  long   Index = -1;
  double X1=0, Y1=0, X2=0, Y2=0;
  double ret=0;

  if (X > Table.max_x) return Table_Index(Table,Table.rows-1  ,j);
  if (X < Table.min_x) return Table_Index(Table,0  ,j);

  // Use constant-time lookup when possible
  if(Table.constantstep) {
    Index = (long)floor(
              (X - Table.min_x) / (Table.max_x - Table.min_x) * (Table.rows-1));
    X1 = Table_Index(Table,Index  ,0);
    X2 = Table_Index(Table,Index+1,0);
  }
  // Use binary search on large, monotonic tables
  else if(Table.monotonic && Table.rows > 100) {
    long left = Table.min_x;
    long right = Table.max_x;

    while (!((X1 <= X) && (X < X2)) && (right - left > 1)) {
      Index = (left + right) / 2;

      X1 = Table_Index(Table, Index-1, 0);
      X2 = Table_Index(Table, Index,   0);

      if (X < X1) {
        right = Index;
      } else {
        left  = Index;
      }
    }
  }

  // Fall back to linear search, if no-one else has set X1, X2 correctly
  if (!((X1 <= X) && (X < X2))) {
    /* look for index surrounding X in the table -> Index */
    for (Index=1; Index <= Table.rows-1; Index++) {
        X1 = Table_Index(Table, Index-1,0);
        X2 = Table_Index(Table, Index  ,0);
        if ((X1 <= X) && (X < X2)) break;
      } /* end for Index */
  }

  Y1 = Table_Index(Table,Index-1,j);
  Y2 = Table_Index(Table,Index  ,j);

  if (!strcmp(Table.method,"linear")) {
    ret = Table_Interp1d(X, X1,Y1, X2,Y2);
  }
  else if (!strcmp(Table.method,"nearest")) {
    ret = Table_Interp1d_nearest(X, X1,Y1, X2,Y2);
  }

  return ret;
} /* end Table_Value */

/*******************************************************************************
* double Table_Value2d(t_Table Table, double X, double Y)
*   ACTION: read element [X,Y] of a matrix Table
*   input   Table: table containing data.
*           X : row index, may be non integer
*           Y : column index, may be non integer
*   return  Value = data[index X][index Y] with bi-linear interpolation
* Returns Value for the indices [X,Y]
* Tests are performed (within Table_Index) on indexes i,j to avoid errors
* NOTE: data should rather be monotonic, and evenly sampled.
*******************************************************************************/
  double Table_Value2d(t_Table Table, double X, double Y)
  {
    long   x1,x2,y1,y2;
    double z11,z12,z21,z22;
    double ret=0;

    x1 = (long)floor(X);
    y1 = (long)floor(Y);

    if (x1 > Table.rows-1 || x1 < 0) {
      x2 = x1;
    } else {
      x2 = x1 + 1;
    }

    if (y1 > Table.columns-1 || y1 < 0) {
      y2 = y1;
    } else {
      y2 = y1 + 1;
    }

    z11 = Table_Index(Table, x1, y1);

    if (y2 != y1) z12=Table_Index(Table, x1, y2); else z12 = z11;
    if (x2 != x1) z21=Table_Index(Table, x2, y1); else z21 = z11;
    if (y2 != y1) z22=Table_Index(Table, x2, y2); else z22 = z21;

    if (!strcmp(Table.method,"linear"))
      ret = Table_Interp2d(X,Y, x1,y1,x2,y2, z11,z12,z21,z22);
    else {
      if (fabs(X-x1) < fabs(X-x2)) {
        if (fabs(Y-y1) < fabs(Y-y2)) ret = z11; else ret = z12;
      } else {
        if (fabs(Y-y1) < fabs(Y-y2)) ret = z21; else ret = z22;
      }
    }
    return ret;
  } /* end Table_Value2d */


/*******************************************************************************
* void Table_Free(t_Table *Table)
*   ACTION: free a single Table
*   return: empty Table
*******************************************************************************/
  void Table_Free(t_Table *Table)
  {
    if( !Table_File_List_gc(Table) ){
       return;
    }
    if (!Table) return;
    if (Table->data   != NULL) free(Table->data);
    if (Table->header != NULL) free(Table->header);
    Table->data   = NULL;
    Table->header = NULL;
  } /* end Table_Free */

/******************************************************************************
* void Table_Info(t_Table Table)
*    ACTION: print informations about a single Table
*******************************************************************************/
  long Table_Info(t_Table Table)
  {
    char buffer[256];
    long ret=0;

    if (!Table.block_number) strcpy(buffer, "catenated");
    else sprintf(buffer, "block %li", Table.block_number);
    printf("Table from file '%s' (%s)",
      Table.filename ? Table.filename : "", buffer);
    if ((Table.data != NULL) && (Table.rows*Table.columns))
    {
      printf(" is %li x %li ", Table.rows, Table.columns);
      if (Table.rows*Table.columns > 1)
        printf("(x=%g:%g)", Table.min_x, Table.max_x);
      else printf("(x=%g) ", Table.min_x);
      ret = Table.rows*Table.columns;
      if (Table.monotonic)    printf(", monotonic");
      if (Table.constantstep) printf(", constant step");
      printf(". interpolation: %s\n", Table.method);
    }
    else printf(" is empty.\n");

    if (Table.header && strlen(Table.header)) {
      char *header;
      int  i;
      header = malloc(80);
      if (!header) return(ret);
      for (i=0; i<80; header[i++]=0);
      strncpy(header, Table.header, 75);
      if (strlen(Table.header) > 75) {
        strcat( header, " ...");
      }
      for (i=0; i<strlen(header); i++)
        if (header[i] == '\n' || header[i] == '\r') header[i] = ';';
      printf("  '%s'\n", header);
      free(header);
    }

    return(ret);
  } /* end Table_Info */

/******************************************************************************
* long Table_Init(t_Table *Table, m, n)
*   ACTION: initialise a Table to empty m by n table
*   return: empty Table
******************************************************************************/
long Table_Init(t_Table *Table, long rows, long columns)
{
  double *data=NULL;
  long   i;

  if (!Table) return(0);

  Table->header  = NULL;
  Table->filename[0]= '\0';
  Table->filesize= 0;
  Table->min_x   = 0;
  Table->max_x   = 0;
  Table->step_x  = 0;
  Table->block_number = 0;
  Table->array_length = 0;
  Table->monotonic    = 0;
  Table->constantstep = 0;
  Table->begin   = 0;
  Table->end     = 0;
  strcpy(Table->method,"linear");

  if (rows*columns >= 1) {
    data    = (double*)malloc(rows*columns*sizeof(double));
    if (data) for (i=0; i < rows*columns; data[i++]=0);
    else {
      fprintf(stderr,"Error: allocating %ld double elements."
                     "Too big (Table_Init).\n", rows*columns);
      rows = columns = 0;
    }
  }
  Table->rows    = (rows >= 1 ? rows : 0);
  Table->columns = (columns >= 1 ? columns : 0);
  Table->data    = data;
  return(Table->rows*Table->columns);
} /* end Table_Init */

/******************************************************************************
* long Table_Write(t_Table Table, char *file, x1,x2, y1,y2)
*   ACTION: write a Table to disk (ascii).
*     when x1=x2=0 or y1=y2=0, the table default limits are used.
*   return: 0=all is fine, non-0: error
*******************************************************************************/
MCDETECTOR Table_Write(t_Table Table, char *file, char *xl, char *yl,
  double x1, double x2, double y1, double y2)
{
  long    i =0;
  MCDETECTOR detector;

  if ((Table.data == NULL) && (Table.rows*Table.columns)) {
    detector.m = 0;
    return(detector); /* Table is empty - nothing to do */
  }
  if (!x1 && !x2) {
    x1 = Table.min_x;
    x2 = Table.max_x;
  }
  if (!y1 && !y2) {
    y1 = 1;
    y2 = Table.columns;
  }

  /* transfer content of the Table into a 2D detector */
  Coords coords = { 0, 0, 0};

  if (Table.rows == 1 || Table.columns == 1) {
    detector = mcdetector_out_1D(Table.filename,
                      xl ? xl : "", yl ? yl : "",
                      "x", x1, x2,
                      Table.rows * Table.columns,
                      NULL, Table.data, NULL,
                      file, file, coords);
  } else {
    detector = mcdetector_out_2D(Table.filename,
                      xl ? xl : "", yl ? yl : "",
                      x1, x2, y1, y2,
                      Table.rows, Table.columns,
                      NULL, Table.data, NULL,
                      file, file, coords);
  }
  return(detector);
}

/******************************************************************************
* void Table_Stat(t_Table *Table)
*   ACTION: computes min/max/mean step of 1st column for a single table (private)
*   return: updated Table
*******************************************************************************/
  static void Table_Stat(t_Table *Table)
  {
    long   i;
    double max_x, min_x;
    double row=1;
    char   monotonic=1;
    char   constantstep=1;
    double step=0;
    long n;

    if (!Table) return;
    if (!Table->rows || !Table->columns) return;
    if (Table->rows == 1) row=0; // single row
    max_x = -FLT_MAX;
    min_x =  FLT_MAX;
    n     = (row ? Table->rows : Table->columns);
    /* get min and max of first column/vector */
    for (i=0; i < n; i++)
    {
      double X;
      X = (row ? Table_Index(*Table,i  ,0)
                               : Table_Index(*Table,0, i));
      if (X < min_x) min_x = X;
      if (X > max_x) max_x = X;
    } /* for */

    /* test for monotonicity and constant step if the table is an XY or single vector */
    if (n > 1) {
      /* mean step */
      step = (max_x - min_x)/(n-1);
      /* now test if table is monotonic on first column, and get minimal step size */
      for (i=0; i < n-1; i++) {
        double X, diff;;
        X    = (row ? Table_Index(*Table,i  ,0)
                    : Table_Index(*Table,0,  i));
        diff = (row ? Table_Index(*Table,i+1,0)
                    : Table_Index(*Table,0,  i+1)) - X;
        if (diff && fabs(diff) < fabs(step)) step = diff;
        /* change sign ? */
        if ((max_x - min_x)*diff < 0 && monotonic)
          monotonic = 0;
      } /* end for */

      /* now test if steps are constant within READ_TABLE_STEPTOL */
      if(!step){
        /*means there's a disconitnuity -> not constantstep*/
        constantstep=0;
      }else if (monotonic) {
        for (i=0; i < n-1; i++) {
          double X, diff;
          X    = (row ? Table_Index(*Table,i  ,0)
              : Table_Index(*Table,0,  i));
          diff = (row ? Table_Index(*Table,i+1,0)
              : Table_Index(*Table,0,  i+1)) - X;
          if ( fabs(step)*(1+READ_TABLE_STEPTOL) < fabs(diff) ||
                fabs(diff) < fabs(step)*(1-READ_TABLE_STEPTOL) )
          { constantstep = 0; break; }
        }
      }

    }
    Table->step_x= step;
    Table->max_x = max_x;
    Table->min_x = min_x;
    Table->monotonic = monotonic;
    Table->constantstep = constantstep;
  } /* end Table_Stat */

/******************************************************************************
* t_Table *Table_Read_Array(char *File, long *blocks)
*   ACTION: read as many data blocks as available, iteratively from file
*   return: initialized t_Table array, last element is an empty Table.
*           the number of extracted blocks in non NULL pointer *blocks
*******************************************************************************/
  t_Table *Table_Read_Array(char *File, long *blocks)
  {
    t_Table *Table_Array=NULL;
    long offset=0;
    long block_number=0;
    long allocated=256;
    long nelements=1;

    /* fisrt allocate an initial empty t_Table array */
    Table_Array = (t_Table *)malloc(allocated*sizeof(t_Table));
    if (!Table_Array) {
      fprintf(stderr, "Error: Can not allocate memory %li (Table_Read_Array).\n",
         allocated*sizeof(t_Table));
      *blocks = 0;
      return (NULL);
    }

    while (nelements > 0)
    {
      t_Table Table;

      /* if ok, set t_Table block number else exit loop */
      block_number++;
      Table.block_number = block_number;

      /* access file at offset and get following block. Block number is from the set offset
       * hence the hardcoded 1 - i.e. the next block counted from offset.*/
      nelements = Table_Read_Offset(&Table, File, 1, &offset,0);
      /* if t_Table array is not long enough, expand and realocate */
      if (block_number >= allocated-1) {
        allocated += 256;
        Table_Array = (t_Table *)realloc(Table_Array,
           allocated*sizeof(t_Table));
        if (!Table_Array) {
          fprintf(stderr, "Error: Can not re-allocate memory %li (Table_Read_Array).\n",
              allocated*sizeof(t_Table));
          *blocks = 0;
          return (NULL);
        }
      }
      /* store it into t_Table array */
      //snprintf(Table.filename, 1024, "%s#%li", File, block_number-1);
      Table_Array[block_number-1] = Table;
      /* continues until we find an empty block */
    }
    /* send back number of extracted blocks */
    if (blocks) *blocks = block_number-1;

    /* now store total number of elements in Table array */
    for (offset=0; offset < block_number;
      Table_Array[offset++].array_length = block_number-1);

    return(Table_Array);
  } /* end Table_Read_Array */
/*******************************************************************************
* void Table_Free_Array(t_Table *Table)
*   ACTION: free a Table array
*******************************************************************************/
  void Table_Free_Array(t_Table *Table)
  {
    long index=0;
    if (!Table) return;
    while (Table[index].data || Table[index].header){
            Table_Free(&Table[index]);
            index++;
    }
    free(Table);
  } /* end Table_Free_Array */

/******************************************************************************
* long Table_Info_Array(t_Table *Table)
*    ACTION: print informations about a Table array
*    return: number of elements in the Table array
*******************************************************************************/
  long Table_Info_Array(t_Table *Table)
  {
    long index=0;

    if (!Table) return(-1);
    while (index < Table[index].array_length
       && (Table[index].data || Table[index].header)
       && (Table[index].rows*Table[index].columns) ) {
      Table_Info(Table[index]);
      index++;
    }
    printf("This Table array contains %li elements\n", index);
    return(index);
  } /* end Table_Info_Array */

/******************************************************************************
* char **Table_ParseHeader(char *header, symbol1, symbol2, ..., NULL)
*    ACTION: search for char* symbols in header and return their value or NULL
*            the search is not case sensitive.
*            Last argument MUST be NULL
*    return: array of char* with line following each symbol, or NULL if not found
*******************************************************************************/
#ifndef MyNL_ARGMAX
#define MyNL_ARGMAX 50
#endif

char **Table_ParseHeader_backend(char *header, ...){
  va_list ap;
  char exit_flag=0;
  int counter   =0;
  char **ret    =NULL;
  if (!header || header[0]=='\0') return(NULL);

  ret = (char**)calloc(MyNL_ARGMAX, sizeof(char*));
  if (!ret) {
    printf("Table_ParseHeader: Cannot allocate %i values array for Parser (Table_ParseHeader).\n",
      MyNL_ARGMAX);
    return(NULL);
  }
  for (counter=0; counter < MyNL_ARGMAX; ret[counter++] = NULL);
  counter=0;

  va_start(ap, header);
  while(!exit_flag && counter < MyNL_ARGMAX-1)
  {
    char *arg_char=NULL;
    char *pos     =NULL;
    /* get variable argument value as a char */
    arg_char = va_arg(ap, char *);
    if (!arg_char || arg_char[0]=='\0'){
      exit_flag = 1; break;
    }
    /* search for the symbol in the header */
    pos = (char*)strcasestr(header, arg_char);
    if (pos) {
      char *eol_pos;
      eol_pos = strchr(pos+strlen(arg_char), '\n');
      if (!eol_pos)
        eol_pos = strchr(pos+strlen(arg_char), '\r');
      if (!eol_pos)
        eol_pos = pos+strlen(pos)-1;
      ret[counter] = (char*)malloc(eol_pos - pos);
      if (!ret[counter]) {
        printf("Table_ParseHeader: Cannot allocate value[%i] array for Parser searching for %s (Table_ParseHeader).\n",
          counter, arg_char);
        exit_flag = 1; break;
      }
      strncpy(ret[counter], pos+strlen(arg_char), eol_pos - pos - strlen(arg_char));
      ret[counter][eol_pos - pos - strlen(arg_char)]='\0';
    }
    counter++;
  }
  va_end(ap);
  return(ret);
} /* Table_ParseHeader */

/******************************************************************************
* double Table_Interp1d(x, x1, y1, x2, y2)
*    ACTION: interpolates linearly at x between y1=f(x1) and y2=f(x2)
*    return: y=f(x) value
*******************************************************************************/
double Table_Interp1d(double x,
  double x1, double y1,
  double x2, double y2)
{
  double slope;
  if (x2 == x1) return (y1+y2)/2;
  if (y1 == y2) return  y1;
  slope = (y2 - y1)/(x2 - x1);
  return y1+slope*(x - x1);
} /* Table_Interp1d */

/******************************************************************************
* double Table_Interp1d_nearest(x, x1, y1, x2, y2)
*    ACTION: table lookup with nearest method at x between y1=f(x1) and y2=f(x2)
*    return: y=f(x) value
*******************************************************************************/
double Table_Interp1d_nearest(double x,
  double x1, double y1,
  double x2, double y2)
{
  if (fabs(x-x1) < fabs(x-x2)) return (y1);
  else return(y2);
} /* Table_Interp1d_nearest */

/******************************************************************************
* double Table_Interp2d(x,y, x1,y1, x2,y2, z11,z12,z21,z22)
*    ACTION: interpolates bi-linearly at (x,y) between z1=f(x1,y1) and z2=f(x2,y2)
*    return: z=f(x,y) value
*    x,y |   x1   x2
*    ----------------
*     y1 |   z11  z21
*     y2 |   z12  z22
*******************************************************************************/
double Table_Interp2d(double x, double y,
  double x1, double y1,
  double x2, double y2,
  double z11, double z12, double z21, double z22)
{
  double ratio_x, ratio_y;
  if (x2 == x1) return Table_Interp1d(y, y1,z11, y2,z12);
  if (y1 == y2) return Table_Interp1d(x, x1,z11, x2,z21);

  ratio_y = (y - y1)/(y2 - y1);
  ratio_x = (x - x1)/(x2 - x1);
  return (1-ratio_x)*(1-ratio_y)*z11 + ratio_x*(1-ratio_y)*z21
    + ratio_x*ratio_y*z22         + (1-ratio_x)*ratio_y*z12;
} /* Table_Interp2d */

/* end of read_table-lib.c */


#ifndef SOURCE_GEN_DEF
#define SOURCE_GEN_DEF
/*******************************************************************************
* str_dup_numeric: replaces non 'valid name' chars with spaces
*******************************************************************************/
char *str_dup_numeric(char *orig)
  {
    long i;

    if (!orig || !strlen(orig)) return(NULL);

    for (i=0; i < strlen(orig); i++)
    {
      if ( (orig[i] > 122)
        || (orig[i] < 32)
        || (strchr("!\"#$%&'()*,:;<=>?@[\\]^`/ ", orig[i]) != NULL) )
      {
        orig[i] = ' ';
      }
    }
    orig[i] = '\0';
    /* now skip spaces */
    for (i=0; i < strlen(orig); i++) {
      if (*orig == ' ') orig++;
      else break;
    }

    return(orig);
  } /* str_dup_numeric */

  /* A normalised Maxwellian distribution : Integral over all l = 1 */
  double SG_Maxwell(double l, double temp)
  {
    double a=949.0/temp;
    return 2*a*a*exp(-a/(l*l))/(l*l*l*l*l);
  }
#endif


/* Shared user declarations for all components types 'Guide_gravity'. */
/*****************************************************************************
*
* McStas, neutron ray-tracing package
*         Copyright 1997-2006, All rights reserved
*         Risoe National Laboratory, Roskilde, Denmark
*         Institut Laue Langevin, Grenoble, France
*
* Library: share/ref-lib.h
*
* %Identification
* Written by: Peter Christiansen
* Date: August, 2006
* Origin: RISOE
* Release: McStas 1.10
* Version: $Revision$
*
* Commonly used reflection functions are declared in this file which
* are used by some guide and mirror components.
*
* Depends on read_table-lib
*
* Usage: within SHARE
* %include "ref-lib"
*
****************************************************************************/


#ifndef REF_LIB_H
#define REF_LIB_H "$Revision$"

void StdReflecFunc(double, double*, double*);
void TableReflecFunc(double, t_Table*, double*);

#endif

/* end of ref-lib.h */
/****************************************************************************
*
* McStas, neutron ray-tracing package
*         Copyright 1997-2006, All rights reserved
*         Risoe National Laboratory, Roskilde, Denmark
*         Institut Laue Langevin, Grenoble, France
*
* Library: share/ref-lib.c
*
* %Identification
* Written by: Peter Christiansen
* Date: August, 2006
* Origin: RISOE
* Release: McStas 1.10
* Version: $Revision$
*
* Commonly used reflection functions are declared in this file which
* are used by some guide and mirror components.
*
* Variable names have prefix 'mc_ref_' for 'McStas Reflection' 
* to avoid conflicts
*
* Usage: within SHARE
* %include "ref-lib"
*
****************************************************************************/

#ifndef REF_LIB_H
#include "ref-lib.h"
#endif

#ifndef READ_TABLE_LIB_H
#include "read_table-lib.h"
#include "read_table-lib.c"
#endif

/****************************************************************************
* void StdReflecFunc(double q, double *par, double *r)
* 
* The McStas standard analytic parametrization of the reflectivity.
* The parameters are:
* R0:      [1]    Low-angle reflectivity
* Qc:      [AA-1] Critical scattering vector
* alpha:   [AA]   Slope of reflectivity
* m:       [1]    m-value of material. Zero means completely absorbing.
* W:       [AA-1] Width of supermirror cut-off
*****************************************************************************/
void StdReflecFunc(double mc_pol_q, double *mc_pol_par, double *mc_pol_r) {
    double R0    = mc_pol_par[0];
    double Qc    = mc_pol_par[1];
    double alpha = mc_pol_par[2];
    double m     = mc_pol_par[3];
    double W     = mc_pol_par[4];
    double beta  = 0;
    mc_pol_q     = fabs(mc_pol_q);
    double arg;
        
    /* Simpler parametrization from Henrik Jacobsen uses these values that depend on m only.
       double m_value=m*0.9853+0.1978;
       double W=-0.0002*m_value+0.0022;
       double alpha=0.2304*m_value+5.0944;
       double beta=-7.6251*m_value+68.1137; 
       If W and alpha are set to 0, use Henrik's approach for estimating these parameters
       and apply the formulation:
       arg = R0*0.5*(1-tanh(arg))*(1-alpha*(q-Qc)+beta*(q-Qc)*(q-Qc));
    */  
    if (W==0 && alpha==0) {
      m=m*0.9853+0.1978;
      W=-0.0002*m+0.0022;
      alpha=0.2304*m+5.0944;
      beta=-7.6251*m+68.1137;
      if (m<=3) {
	alpha=m;
	beta=0;
      }
    }
    
    arg = W > 0 ? (mc_pol_q - m*Qc)/W : 11;

    if (arg > 10 || m <= 0 || Qc <=0 || R0 <= 0) {
      *mc_pol_r = 0;
      return;
    }
    
    if (m < 1) { Qc *= m; m=1; }
    
    if(mc_pol_q <= Qc) {      
      *mc_pol_r = R0;
      return;
    }
    
    
    *mc_pol_r = R0*0.5*(1 - tanh(arg))*(1 - alpha*(mc_pol_q - Qc) + beta*(mc_pol_q - Qc)*(mc_pol_q - Qc));
    
    return;
  }

/****************************************************************************
* void TableReflecFunc(double q, t_Table *par, double *r) {
* 
* Looks up the reflectivity in a table using the routines in read_table-lib.
*****************************************************************************/
void TableReflecFunc(double mc_pol_q, t_Table *mc_pol_par, double *mc_pol_r) {
    
  *mc_pol_r = Table_Value(*mc_pol_par, mc_pol_q, 1);
  if(*mc_pol_r>1)
    *mc_pol_r = 1;
  return;
}

/* end of ref-lib.c */

#ifndef Gravity_guide_Version
#define Gravity_guide_Version "$Revision$"

#ifndef PROP_GRAV_DT
#error McStas : You need PROP_GRAV_DT (McStas >= 1.4.3) to run this component
#endif

/*
* G:       (m/s^2) Gravitation acceleration along y axis [-9.81]
* Gx:      (m/s^2) Gravitation acceleration along x axis [0]
* Gy:      (m/s^2) Gravitation acceleration along y axis [-9.81]
* Gz:      (m/s^2) Gravitation acceleration along z axis [0]
* mh:      (1)    m-value of material for left/right vert. mirrors
* mv:      (1)    m-value of material for top/bottom horz. mirrors
* mx:      (1)    m-value of material for left/right vert. mirrors
* my:      (1)    m-value of material for top/bottom horz. mirrors
*/

  typedef struct Gravity_guide_Vars
  {
    double gx;
    double gy;
    double gz;
    double nx[6], ny[6], nz[6];
    double wx[6], wy[6], wz[6];
    double A[6], norm_n2[6], norm_n[6];
    long   N_reflection[7];
    double w1c, h1c;
    double w2c, h2c;
    double M[5];
    double Alpha[5];
    double nzC[5], norm_n2xy[5], Axy[5];
    double wav_lr, wav_tb, wav_z;
    double chamfer_z, chamfer_lr, chamfer_tb;
    char   compcurname[256];
    double fc_freq, fc_phase;
    double warnings;
  } Gravity_guide_Vars_type;

  void Gravity_guide_Init(Gravity_guide_Vars_type *aVars,
    MCNUM a_w1, MCNUM a_h1, MCNUM a_w2, MCNUM a_h2, MCNUM a_l, MCNUM a_R0,
    MCNUM a_Qc, MCNUM a_alpha, MCNUM a_m, MCNUM a_W, MCNUM a_nslit, MCNUM a_d,
    MCNUM a_Gx, MCNUM a_Gy, MCNUM a_Gz,
    MCNUM a_mleft, MCNUM a_mright, MCNUM a_mtop, MCNUM a_mbottom, MCNUM a_nhslit,
    MCNUM a_wavy_lr, MCNUM a_wavy_tb, MCNUM a_wavy_z, MCNUM a_wavy,
    MCNUM a_chamfers_z, MCNUM a_chamfers_lr, MCNUM a_chamfers_tb, MCNUM a_chamfers,
    MCNUM a_nu, MCNUM a_phase, MCNUM a_aleft, MCNUM a_aright, MCNUM a_atop, MCNUM a_abottom)
  {
    int i;

    for (i=0; i<7; aVars->N_reflection[i++] = 0);
    for (i=0; i<5; aVars->M[i++] = 0);
    for (i=0; i<5; aVars->Alpha[i++] = 0);

    aVars->gx = a_Gx; /* The gravitation vector in the current component axis system */
    aVars->gy = a_Gy;
    aVars->gz = a_Gz;
    aVars->warnings=0;

    if (a_nslit <= 0 || a_nhslit <= 0) { fprintf(stderr,"%s: Fatal: no channel in this guide (nhslit or nslit=0).\n", aVars->compcurname); exit(-1); }
    if (a_d < 0) { fprintf(stderr,"%s: Fatal: subdividing walls have negative thickness in this guide (d<0).\n", aVars->compcurname); exit(-1); }
    aVars->w1c = (a_w1 - (a_nslit-1) *a_d)/(double)a_nslit;
    aVars->w2c = (a_w2 - (a_nslit-1) *a_d)/(double)a_nslit;
    aVars->h1c = (a_h1 - (a_nhslit-1)*a_d)/(double)a_nhslit;
    aVars->h2c = (a_h2 - (a_nhslit-1)*a_d)/(double)a_nhslit;

    for (i=0; i <= 4;   aVars->M[i++]=a_m);
    for (i=0; i <= 4;   aVars->Alpha[i++]=a_alpha);
    if (a_mleft   >= 0) aVars->M[1] =a_mleft  ;
    if (a_mright  >= 0) aVars->M[2] =a_mright ;
    if (a_mtop    >= 0) aVars->M[3] =a_mtop   ;
    if (a_mbottom >= 0) aVars->M[4] =a_mbottom;
    if (a_aleft   >= 0) aVars->Alpha[1] =a_aleft  ;
    if (a_aright  >= 0) aVars->Alpha[2] =a_aright ;
    if (a_atop    >= 0) aVars->Alpha[3] =a_atop   ;
    if (a_abottom >= 0) aVars->Alpha[4] =a_abottom;

    /* n: normal vectors to surfaces */
    aVars->nx[1] =  a_l; aVars->ny[1] =  0;   aVars->nz[1] =  0.5*(aVars->w2c-aVars->w1c);  /* 1:+X left       */
    aVars->nx[2] = -a_l; aVars->ny[2] =  0;   aVars->nz[2] = -aVars->nz[1];             /* 2:-X right      */
    aVars->nx[3] =  0;   aVars->ny[3] =  a_l; aVars->nz[3] =  0.5*(aVars->h2c-aVars->h1c);  /* 3:+Y top        */
    aVars->nx[4] =  0;   aVars->ny[4] = -a_l; aVars->nz[4] = -aVars->nz[3];             /* 4:-Y bottom     */
    aVars->nx[5] =  0;   aVars->ny[5] =  0;   aVars->nz[5] =  a_l;                      /* 5:+Z exit       */
    aVars->nx[0] =  0;   aVars->ny[0] =  0;   aVars->nz[0] = -a_l;                      /* 0:Z0 input      */
    /* w: a point on these surfaces */
    aVars->wx[1] = +(aVars->w1c)/2; aVars->wy[1] =  0;              aVars->wz[1] = 0;   /* 1:+X left       */
    aVars->wx[2] = -(aVars->w1c)/2; aVars->wy[2] =  0;              aVars->wz[2] = 0;   /* 2:-X right      */
    aVars->wx[3] =  0;              aVars->wy[3] = +(aVars->h1c)/2; aVars->wz[3] = 0;   /* 3:+Y top        */
    aVars->wx[4] =  0;              aVars->wy[4] = -(aVars->h1c)/2; aVars->wz[4] = 0;   /* 4:-Y bottom     */
    aVars->wx[5] =  0;              aVars->wy[5] =  0;              aVars->wz[5] = a_l; /* 5:+Z exit       */
    aVars->wx[0] =  0;              aVars->wy[0] =  0;              aVars->wz[0] = 0;   /* 0:Z0 input      */

    for (i=0; i <= 5; i++)
    {
      aVars->A[i] = scalar_prod(aVars->nx[i], aVars->ny[i], aVars->nz[i], aVars->gx, aVars->gy, aVars->gz)/2;
      aVars->norm_n2[i] = aVars->nx[i]*aVars->nx[i] + aVars->ny[i]*aVars->ny[i] + aVars->nz[i]*aVars->nz[i];
      if (aVars->norm_n2[i] <= 0)
        { fprintf(stderr,"%s: Fatal: normal vector norm %i is null/negative ! check guide dimensions.\n", aVars->compcurname, i); exit(-1); } /* should never occur */
      else
        aVars->norm_n[i] = sqrt(aVars->norm_n2[i]);
    }
    /* partial computations for l/r/t/b sides, to save computing time */
    for (i=1; i <= 4; i++)
    { /* stores nz that changes in case non box element (focus/defocus) */
      aVars->nzC[i]      =  aVars->nz[i]; /* partial xy terms */
      aVars->norm_n2xy[i]=  aVars->nx[i]*aVars->nx[i] + aVars->ny[i]*aVars->ny[i];
      aVars->Axy[i]      = (aVars->nx[i]*aVars->gx    + aVars->ny[i]*aVars->gy)/2;
    }
    /* handle waviness init */
    if (a_wavy && (!a_wavy_tb && !a_wavy_lr && !a_wavy_z))
    { aVars->wav_tb=aVars->wav_lr=aVars->wav_z=a_wavy; }
    else
    { aVars->wav_tb=a_wavy_tb; aVars->wav_lr=a_wavy_lr; aVars->wav_z=a_wavy_z; }
    aVars->wav_tb *= DEG2RAD/(sqrt(8*log(2)));   /* Convert from deg FWHM to rad Gaussian sigma */
    aVars->wav_lr *= DEG2RAD/(sqrt(8*log(2)));
    aVars->wav_z  *= DEG2RAD/(sqrt(8*log(2)));
    /* handle chamfers init */
    if (a_chamfers && (!a_chamfers_z && !a_chamfers_lr && !a_chamfers_tb))
    { aVars->chamfer_z=aVars->chamfer_lr=aVars->chamfer_tb=a_chamfers; }
    else
    {
      aVars->chamfer_z=a_chamfers_z;
      aVars->chamfer_lr=a_chamfers_lr;
      aVars->chamfer_tb=a_chamfers_tb;
    }

    aVars->fc_freq  = a_nu;
    aVars->fc_phase = a_phase;
  }

  int Gravity_guide_Trace(double *dt,
        Gravity_guide_Vars_type *aVars,
        double cx, double cy, double cz,
        double cvx, double cvy, double cvz,
        double cxnum, double cxk, double cynum, double cyk,
        double *cnx, double *cny,double *cnz)
  {
    double B, C;
    int    ret=0;
    int    side=0;
    double n1;
    double dt0, dt_min=0;
    int    i;
    double loc_num, loc_nslit;
    int    i_slope=3;

    /* look if there is a previous intersection with guide sides */
    /* A = 0.5 n.g; B = n.v; C = n.(r-W); */
    /* 5=+Z side: n=(0, 0, -l) ; W = (0, 0, l) (at z=l, guide exit)*/
    B = aVars->nz[5]*cvz; C = aVars->nz[5]*(cz - aVars->wz[5]);
    ret = solve_2nd_order(&dt0, NULL, aVars->A[5], B, C);
    if (ret && dt0>1e-10) { dt_min = dt0; side=5; }

    loc_num = cynum; loc_nslit = cyk;
    for (i=4; i>0; i--)
    {
      if (i == 2) { i_slope=1; loc_num = cxnum; loc_nslit = cxk; }

      if (aVars->nzC[i_slope] != 0) {
        n1 = loc_nslit - 2*(loc_num);  /* slope of l/r/u/d sides depends on the channel ! */
        loc_num++; /* use partial computations to alter nz and A */
        aVars->nz[i]= aVars->nzC[i]*n1;
        aVars->A[i] = aVars->Axy[i] + aVars->nz[i]*aVars->gz/2;
      }
      if (i < 3)
      {      B = aVars->nx[i]*cvx + aVars->nz[i]*cvz; C = aVars->nx[i]*(cx-aVars->wx[i]) + aVars->nz[i]*cz; }
      else { B = aVars->ny[i]*cvy + aVars->nz[i]*cvz; C = aVars->ny[i]*(cy-aVars->wy[i]) + aVars->nz[i]*cz; }
      ret = solve_2nd_order(&dt0, NULL, aVars->A[i], B, C);
      if (ret && dt0>1e-10 && (dt0<dt_min || !dt_min))
      { dt_min = dt0; side=i;
        if (aVars->nzC[i] != 0)
        { aVars->norm_n2[i] = aVars->norm_n2xy[i] + aVars->nz[i]*aVars->nz[i];
          aVars->norm_n[i]  = sqrt(aVars->norm_n2[i]); }
      }
     }

    *dt = dt_min;
    /* handles waviness: rotate n vector */
    if (side > 0 && side < 5 && (aVars->wav_z || aVars->wav_lr || aVars->wav_tb))
    {
      double nt_x, nt_y, nt_z;  /* transverse vector */
      double nn_x, nn_y, nn_z;  /* normal vector (tmp) */
      double phi;
      /* normal vector n_z = [ 0,0,1], n_t = n x n_z; */
      vec_prod(nt_x,nt_y,nt_z, aVars->nx[side],aVars->ny[side],aVars->nz[side], 0,0,1);
      /* rotate n with angle wavy_z around n_t -> nn */
      if (aVars->wav_z) {
        phi = aVars->wav_z;
        rotate(nn_x,nn_y,nn_z, aVars->nx[side],aVars->ny[side],aVars->nz[side], aVars->wav_z*randnorm(), nt_x,nt_y,nt_z);
      } else { nn_x=aVars->nx[side]; nn_y=aVars->ny[side]; nn_z=aVars->nz[side]; }
      /* rotate n with angle wavy_{x|y} around n_z -> nt */
      phi = (side <=2) ? aVars->wav_lr : aVars->wav_tb;
      if (phi) {
        rotate(nt_x,nt_y,nt_z, nn_x,nn_y,nn_z, phi*randnorm(), 0,0,1);
      } else { nt_x=nn_x; nt_y=nn_y; nt_z=nn_z; }
      *cnx=nt_x; *cny=nt_y; *cnz=nt_z;
    } else
    { *cnx=aVars->nx[side]; *cny=aVars->ny[side]; *cnz=aVars->nz[side]; }
    return (side);
  }



#endif

/* Shared user declarations for all components types 'Monitor_nD'. */
/*******************************************************************************
*
* McStas, neutron ray-tracing package
*         Copyright 1997-2002, All rights reserved
*         Risoe National Laboratory, Roskilde, Denmark
*         Institut Laue Langevin, Grenoble, France
*
* Library: share/monitor_nd-lib.h
*
* %Identification
* Written by: EF
* Date: Aug 28, 2002
* Origin: ILL
* Release: McStas 1.6
* Version: $Revision$
*
* This file is to be imported by the monitor_nd related components
* It handles some shared functions.
*
* Usage: within SHARE
* %include "monitor_nd-lib"
*
*******************************************************************************/

#ifndef MONITOR_ND_LIB_H

#define MONITOR_ND_LIB_H "$Revision$"
#define MONnD_COORD_NMAX  30  /* max number of variables to record */

  typedef struct MonitornD_Defines
  {
    int COORD_NONE  ;
    int COORD_X     ;
    int COORD_Y     ;
    int COORD_Z     ;
    int COORD_RADIUS; 
    int COORD_VX    ;
    int COORD_VY    ;
    int COORD_VZ    ;
    int COORD_V     ;
    int COORD_T     ;
    int COORD_P     ;
    int COORD_SX    ;
    int COORD_SY    ;
    int COORD_SZ    ;
    int COORD_KX    ;
    int COORD_KY    ;
    int COORD_KZ    ;
    int COORD_K     ;
    int COORD_ENERGY;
    int COORD_LAMBDA;
    int COORD_KXY   ;
    int COORD_KYZ   ;
    int COORD_KXZ   ;
    int COORD_VXY   ;
    int COORD_VYZ   ;
    int COORD_VXZ   ;
    int COORD_HDIV  ;
    int COORD_VDIV  ;
    int COORD_ANGLE ;
    int COORD_NCOUNT;
    int COORD_THETA ;
    int COORD_PHI   ;
    int COORD_USER1 ;
    int COORD_USER2 ;
    int COORD_USER3 ;
    int COORD_XY    ;
    int COORD_XZ    ;
    int COORD_YZ    ;
    int COORD_PIXELID;

    /* token modifiers */
    int COORD_VAR   ; /* next token should be a variable or normal option */
    int COORD_MIN   ; /* next token is a min value */
    int COORD_MAX   ; /* next token is a max value */
    int COORD_DIM   ; /* next token is a bin value */
    int COORD_FIL   ; /* next token is a filename */
    int COORD_EVNT  ; /* next token is a buffer size value */
    int COORD_3HE   ; /* next token is a 3He pressure value */
    int COORD_LOG   ; /* next variable will be in log scale */
    int COORD_ABS   ; /* next variable will be in abs scale */
    int COORD_SIGNAL; /* next variable will be the signal var */
    int COORD_AUTO  ; /* set auto limits */

    char TOKEN_DEL[32]; /* token separators */

    char SHAPE_SQUARE; /* shape of the monitor */
    char SHAPE_DISK  ;
    char SHAPE_SPHERE;
    char SHAPE_CYLIND;
    char SHAPE_BANANA; /* cylinder without top/bottom, on restricted angular area */
    char SHAPE_BOX   ;
    char SHAPE_PREVIOUS;
    char SHAPE_OFF;

  } MonitornD_Defines_type;

  typedef struct MonitornD_Variables
  {
    double area;
    double Sphere_Radius     ;
    double Cylinder_Height   ;
    char   Flag_With_Borders ;   /* 2 means xy borders too */
    char   Flag_List         ;   /* 1 store 1 buffer, 2 is list all, 3 list all+append */
    char   Flag_Multiple     ;   /* 1 when n1D, 0 for 2D */
    char   Flag_Verbose      ;
    int    Flag_Shape        ;
    char   Flag_Auto_Limits  ;   /* get limits from first Buffer */
    char   Flag_Absorb       ;   /* monitor is also a slit */
    char   Flag_per_cm2      ;   /* flux is per cm2 */
    char   Flag_log          ;   /* log10 of the flux */
    char   Flag_parallel     ;   /* set neutron state back after detection (parallel components) */
    char   Flag_Binary_List  ;
    char   Flag_capture      ;   /* lambda monitor with lambda/lambda(2200m/s = 1.7985 Angs) weightening */
    int    Flag_signal       ;   /* 0:monitor p, else monitor a mean value */
    int    Flag_mantid       ;   /* 0:normal monitor, else do mantid-event specifics */
    int    Flag_OFF          ;   /* Flag to indicate external geometry from OFF file */
    unsigned long OFF_polyidx;   /* When intersection is done externally by off_intersect, this gives the 
				    polygon number, i.e. pixel index */

    unsigned long Coord_Number      ;   /* total number of variables to monitor, plus intensity (0) */
    unsigned long Coord_NumberNoPixel;  /* same but without counting PixelID */
    unsigned long Buffer_Block      ;   /* Buffer size for list or auto limits */
    unsigned long Neutron_Counter   ;   /* event counter, simulation total counts is mcget_ncount() */
    unsigned long Buffer_Counter    ;   /* index in Buffer size (for realloc) */
    unsigned long Buffer_Size       ;
    int    Coord_Type[MONnD_COORD_NMAX];      /* type of variable */
    char   Coord_Label[MONnD_COORD_NMAX][30]; /* label of variable */
    char   Coord_Var[MONnD_COORD_NMAX][30];   /* short id of variable */
    long   Coord_Bin[MONnD_COORD_NMAX];       /* bins of variable array */
    long   Coord_BinProd[MONnD_COORD_NMAX];   /* product of bins of variable array */
    double Coord_Min[MONnD_COORD_NMAX];
    double Coord_Max[MONnD_COORD_NMAX];
    char   Monitor_Label[MONnD_COORD_NMAX*30];/* Label for monitor */
    char   Mon_File[128];                     /* output file name */

    double cx,cy,cz;
    double cvx, cvy, cvz;
    double ckx, cky, ckz;
    double csx, csy, csz;
    double cEx, cEy, cEz;
    double cs1, cs2, ct, cphi, cp;
    double He3_pressure;
    char   Flag_UsePreMonitor    ;   /* use a previously stored neutron parameter set */
    char   UserName1[128];
    char   UserName2[128];
    char   UserName3[128];
    double UserVariable1;
    double UserVariable2;
    double UserVariable3;
    char   option[CHAR_BUF_LENGTH];

    long long int Nsum;
    double psum, p2sum;
    double **Mon2D_N;
    double **Mon2D_p;
    double **Mon2D_p2;
    double *Mon2D_Buffer;
    unsigned long PixelID;

    double mxmin,mxmax,mymin,mymax,mzmin,mzmax;
    double mean_dx, mean_dy, min_x, min_y, max_x, max_y, mean_p;

    char   compcurname[128];
    Coords compcurpos;

  } MonitornD_Variables_type;

/* monitor_nd-lib function prototypes */
/* ========================================================================= */

void Monitor_nD_Init(MonitornD_Defines_type *, MonitornD_Variables_type *, MCNUM, MCNUM, MCNUM, MCNUM, MCNUM, MCNUM, MCNUM, MCNUM, MCNUM, int);
int Monitor_nD_Trace(MonitornD_Defines_type *, MonitornD_Variables_type *);
MCDETECTOR Monitor_nD_Save(MonitornD_Defines_type *, MonitornD_Variables_type *);
void Monitor_nD_Finally(MonitornD_Defines_type *, MonitornD_Variables_type *);
void Monitor_nD_McDisplay(MonitornD_Defines_type *,
 MonitornD_Variables_type *);

#define MONND_DECLARE(monname) \
  struct MonitornD_Variables *mcmonnd ## monname;
#define MONND_USER_TITLE(monname, num, title) \
  { mcmonnd ## monname = &(MC_GETPAR(monname, Vars)); \
    strcpy(mcmonnd ## monname->UserName ## num, title); }
#define MONND_USER_VALUE(monname, num, value) \
  { mcmonnd ## monname = &(MC_GETPAR(monname, Vars)); \
    mcmonnd ## monname->UserVariable ## num = (value); }

#endif

/* end of monitor_nd-lib.h */
/*******************************************************************************
*
* McStas, neutron ray-tracing package
*         Copyright 1997-2002, All rights reserved
*         Risoe National Laboratory, Roskilde, Denmark
*         Institut Laue Langevin, Grenoble, France
*
* Library: share/monitor_nd-lib.c
*
* %Identification
* Written by: EF
* Date: Aug 28, 2002
* Origin: ILL
* Release: McStas 1.6
* Version: $Revision$
*
* This file is to be imported by the monitor_nd related components
* It handles some shared functions. Embedded within instrument in runtime mode.
*
* Usage: within SHARE
* %include "monitor_nd-lib"
*
*******************************************************************************/

#ifndef MONITOR_ND_LIB_H
#error McStas : please import this library with %include "monitor_nd-lib"
#endif

/* ========================================================================= */
/* Monitor_nD_Init: this routine is used to parse options                    */
/* ========================================================================= */

void Monitor_nD_Init(MonitornD_Defines_type *DEFS,
  MonitornD_Variables_type *Vars,
  MCNUM xwidth,
  MCNUM yheight,
  MCNUM zdepth,
  MCNUM xmin,
  MCNUM xmax,
  MCNUM ymin,
  MCNUM ymax,
  MCNUM zmin,
  MCNUM zmax,
  int offflag)
  {
    long carg = 1;
    char *option_copy, *token;
    char Flag_New_token = 1;
    char Flag_End       = 1;
    char Flag_All       = 0;
    char Flag_No        = 0;
    char Flag_abs       = 0;
    int  Flag_auto      = 0;  /* -1: all, 1: the current variable */
    int  Set_Vars_Coord_Type;
    char Set_Vars_Coord_Label[64];
    char Set_Vars_Coord_Var[64];
    char Short_Label[MONnD_COORD_NMAX][64];
    int  Set_Coord_Mode;
    long i=0, j=0;
    double lmin, lmax, XY=0;
    long t;


    t = (long)time(NULL);

/* initialize DEFS */
/* Variables to monitor */
    DEFS->COORD_NONE   =0;
    DEFS->COORD_X      =1;
    DEFS->COORD_Y      =2;
    DEFS->COORD_Z      =3;
    DEFS->COORD_RADIUS =19;
    DEFS->COORD_VX     =4;
    DEFS->COORD_VY     =5;
    DEFS->COORD_VZ     =6;
    DEFS->COORD_V      =16;
    DEFS->COORD_T      =7;
    DEFS->COORD_P      =8;
    DEFS->COORD_SX     =9;
    DEFS->COORD_SY     =10;
    DEFS->COORD_SZ     =11;
    DEFS->COORD_KX     =12;
    DEFS->COORD_KY     =13;
    DEFS->COORD_KZ     =14;
    DEFS->COORD_K      =15;
    DEFS->COORD_ENERGY =17;
    DEFS->COORD_LAMBDA =18;
    DEFS->COORD_HDIV   =20;
    DEFS->COORD_VDIV   =21;
    DEFS->COORD_ANGLE  =22;
    DEFS->COORD_NCOUNT =23;
    DEFS->COORD_THETA  =24;
    DEFS->COORD_PHI    =25;
    DEFS->COORD_USER1  =26;
    DEFS->COORD_USER2  =27;
    DEFS->COORD_USER3  =28;
    DEFS->COORD_XY     =37;
    DEFS->COORD_YZ     =31;
    DEFS->COORD_XZ     =32;
    DEFS->COORD_VXY    =30;
    DEFS->COORD_VYZ    =34;
    DEFS->COORD_VXZ    =36;
    DEFS->COORD_KXY    =29;
    DEFS->COORD_KYZ    =33;
    DEFS->COORD_KXZ    =35;
    DEFS->COORD_PIXELID=38;

/* token modifiers */
    DEFS->COORD_VAR    =0;    /* next token should be a variable or normal option */
    DEFS->COORD_MIN    =1;    /* next token is a min value */
    DEFS->COORD_MAX    =2;    /* next token is a max value */
    DEFS->COORD_DIM    =3;    /* next token is a bin value */
    DEFS->COORD_FIL    =4;    /* next token is a filename */
    DEFS->COORD_EVNT   =5;    /* next token is a buffer size value */
    DEFS->COORD_3HE    =6;    /* next token is a 3He pressure value */
    DEFS->COORD_LOG    =64;   /* next variable will be in log scale */
    DEFS->COORD_ABS    =128;  /* next variable will be in abs scale */
    DEFS->COORD_SIGNAL =256;  /* next variable will be the signal var */
    DEFS->COORD_AUTO   =512;  /* set auto limits */

    strcpy(DEFS->TOKEN_DEL, " =,;[](){}:");  /* token separators */

    DEFS->SHAPE_SQUARE =0;    /* shape of the monitor */
    DEFS->SHAPE_DISK   =1;
    DEFS->SHAPE_SPHERE =2;
    DEFS->SHAPE_CYLIND =3;
    DEFS->SHAPE_BANANA =4;
    DEFS->SHAPE_BOX    =5;
    DEFS->SHAPE_PREVIOUS=6;
    DEFS->SHAPE_OFF=7;

    Vars->Sphere_Radius     = 0;
    Vars->Cylinder_Height   = 0;
    Vars->Flag_With_Borders = 0;   /* 2 means xy borders too */
    Vars->Flag_List         = 0;   /* 1=store 1 buffer, 2=list all, 3=re-use buffer */
    Vars->Flag_Multiple     = 0;   /* 1 when n1D, 0 for 2D */
    Vars->Flag_Verbose      = 0;
    Vars->Flag_Shape        = DEFS->SHAPE_SQUARE;
    Vars->Flag_Auto_Limits  = 0;   /* get limits from first Buffer */
    Vars->Flag_Absorb       = 0;   /* monitor is also a slit */
    Vars->Flag_per_cm2      = 0;   /* flux is per cm2 */
    Vars->Flag_log          = 0;   /* log10 of the flux */
    Vars->Flag_parallel     = 0;   /* set neutron state back after detection (parallel components) */
    Vars->Flag_Binary_List  = 0;   /* save list as a binary file (smaller) */
    Vars->Coord_Number      = 0;   /* total number of variables to monitor, plus intensity (0) */
    Vars->Coord_NumberNoPixel=0;   /* same but without counting PixelID */

/* Allow to specify size of Monitor_nD buffer via a define*/
#ifndef MONND_BUFSIZ
    Vars->Buffer_Block      = 100000;     /* Buffer size for list or auto limits */
#else
	Vars->Buffer_Block      = MONND_BUFSIZ;     /* Buffer size for list or auto limits */	
#endif
    Vars->Neutron_Counter   = 0;   /* event counter, simulation total counts is mcget_ncount() */
    Vars->Buffer_Counter    = 0;   /* index in Buffer size (for realloc) */
    Vars->Buffer_Size       = 0;
    Vars->UserVariable1     = 0;
    Vars->UserVariable2     = 0;
    Vars->He3_pressure      = 0;
    Vars->Flag_capture      = 0;
    Vars->Flag_signal       = DEFS->COORD_P;
    Vars->Flag_mantid       = 0;
    Vars->Flag_OFF          = offflag;
    Vars->OFF_polyidx       = -1;
    Vars->mean_dx=Vars->mean_dy=0;
    Vars->min_x = Vars->max_x  =0;
    Vars->min_y = Vars->max_y  =0;

    Set_Vars_Coord_Type = DEFS->COORD_NONE;
    Set_Coord_Mode = DEFS->COORD_VAR;

    /* handle size parameters */
    /* normal use is with xwidth, yheight, zdepth */
    /* if xmin,xmax,ymin,ymax,zmin,zmax are non 0, use them */
    if (fabs(xmin-xmax) == 0)
      { Vars->mxmin = -fabs(xwidth)/2; Vars->mxmax = fabs(xwidth)/2; }
    else
      { if (xmin < xmax) {Vars->mxmin = xmin; Vars->mxmax = xmax;}
        else {Vars->mxmin = xmax; Vars->mxmax = xmin;}
      }
    if (fabs(ymin-ymax) == 0)
      { Vars->mymin = -fabs(yheight)/2; Vars->mymax = fabs(yheight)/2; }
    else
      { if (ymin < ymax) {Vars->mymin = ymin; Vars->mymax = ymax;}
        else {Vars->mymin = ymax; Vars->mymax = ymin;}
      }
    if (fabs(zmin-zmax) == 0)
      { Vars->mzmin = -fabs(zdepth)/2; Vars->mzmax = fabs(zdepth)/2; }
    else
      { if (zmin < zmax) {Vars->mzmin = zmin; Vars->mzmax = zmax; }
        else {Vars->mzmin = zmax; Vars->mzmax = zmin; }
      }

    if (fabs(Vars->mzmax-Vars->mzmin) == 0)
      Vars->Flag_Shape        = DEFS->SHAPE_SQUARE;
    else
      Vars->Flag_Shape        = DEFS->SHAPE_BOX;

    if (Vars->Flag_OFF) {
      Vars->Flag_Shape        = DEFS->SHAPE_OFF;
    }
    
    /* parse option string */

    option_copy = (char*)malloc(strlen(Vars->option)+1);
    if (option_copy == NULL)
    {
      fprintf(stderr,"Monitor_nD: %s cannot allocate 'options' copy (%li). Fatal.\n", Vars->compcurname, (long)strlen(Vars->option));
      exit(-1);
    }

    if (strlen(Vars->option))
    {
      Flag_End = 0;
      strcpy(option_copy, Vars->option);
    }

    if (strstr(Vars->option, "cm2") || strstr(Vars->option, "cm^2")) Vars->Flag_per_cm2 = 1;

    if (strstr(Vars->option, "binary") || strstr(Vars->option, "float"))
      Vars->Flag_Binary_List  = 1;
    if (strstr(Vars->option, "double"))
      Vars->Flag_Binary_List  = 2;

    strcpy(Vars->Coord_Label[0],"Intensity");
    strncpy(Vars->Coord_Var[0],"p",30);
    Vars->Coord_Type[0] = DEFS->COORD_P;
    Vars->Coord_Bin[0] = 1;
    Vars->Coord_Min[0] = 0;
    Vars->Coord_Max[0] = FLT_MAX;

    /* default file name is comp_name+dateID */
    sprintf(Vars->Mon_File, "%s_%li", Vars->compcurname, t);

    carg = 1;
    while((Flag_End == 0) && (carg < 128))
    {

      if (Flag_New_token) /* retain previous token or get a new one */
      {
        if (carg == 1) token=(char *)strtok(option_copy,DEFS->TOKEN_DEL);
        else token=(char *)strtok(NULL,DEFS->TOKEN_DEL);
        if (token == NULL) Flag_End=1;
      }
      Flag_New_token = 1;
      if ((token != NULL) && (strlen(token) != 0))
      {
        char iskeyword=0; /* left at 0 when variables are processed, 1 for modifiers */
        int  old_Mode;
        /* change token to lower case */
        for (i=0; i<strlen(token); i++) token[i]=tolower(token[i]);
        /* first handle option values from preceeding keyword token detected */
        old_Mode = Set_Coord_Mode;
        if (Set_Coord_Mode == DEFS->COORD_MAX)  /* max=%i */
        {
          if (!Flag_All)
            Vars->Coord_Max[Vars->Coord_Number] = atof(token);
          else
            for (i = 0; i <= Vars->Coord_Number; Vars->Coord_Max[i++] = atof(token));
          Set_Coord_Mode = DEFS->COORD_VAR; Flag_All = 0;
        }
        if (Set_Coord_Mode == DEFS->COORD_MIN)  /* min=%i */
        {
          if (!Flag_All)
            Vars->Coord_Min[Vars->Coord_Number] = atof(token);
          else
            for (i = 0; i <= Vars->Coord_Number; Vars->Coord_Min[i++] = atof(token));
          Set_Coord_Mode = DEFS->COORD_MAX;
        }
        if (Set_Coord_Mode == DEFS->COORD_DIM)  /* bins=%i */
        {
          if (!Flag_All)
            Vars->Coord_Bin[Vars->Coord_Number] = atoi(token);
          else
            for (i = 0; i <= Vars->Coord_Number; Vars->Coord_Bin[i++] = atoi(token));
          Set_Coord_Mode = DEFS->COORD_VAR; Flag_All = 0;
        }
        if (Set_Coord_Mode == DEFS->COORD_FIL)  /* file=%s */
        {
          if (!Flag_No) strncpy(Vars->Mon_File,token,128);
          else { strcpy(Vars->Mon_File,""); Vars->Coord_Number = 0; Flag_End = 1;}
          Set_Coord_Mode = DEFS->COORD_VAR;
        }
        if (Set_Coord_Mode == DEFS->COORD_EVNT) /* list=%i */
        {
          if (!strcmp(token, "all") || Flag_All) Vars->Flag_List = 2;
          else { i = (long)ceil(atof(token)); if (i) Vars->Buffer_Block = i;
            Vars->Flag_List = 1; }
          Set_Coord_Mode = DEFS->COORD_VAR; Flag_All = 0;
        }
        if (Set_Coord_Mode == DEFS->COORD_3HE)  /* pressure=%g */
        {
            Vars->He3_pressure = atof(token);
            Set_Coord_Mode = DEFS->COORD_VAR; Flag_All = 0;
        }

        /* now look for general option keywords */
        if (!strcmp(token, "borders"))  {Vars->Flag_With_Borders = 1; iskeyword=1; }
        if (!strcmp(token, "verbose"))  {Vars->Flag_Verbose      = 1; iskeyword=1; }
        if (!strcmp(token, "log"))      {Vars->Flag_log          = 1; iskeyword=1; }
        if (!strcmp(token, "abs"))      {Flag_abs                = 1; iskeyword=1; }
        if (!strcmp(token, "multiple")) {Vars->Flag_Multiple     = 1; iskeyword=1; }
        if (!strcmp(token, "list") || !strcmp(token, "events")) {
          Vars->Flag_List = 1; Set_Coord_Mode = DEFS->COORD_EVNT;  }
        if (!strcmp(token, "limits") || !strcmp(token, "min"))
          Set_Coord_Mode = DEFS->COORD_MIN;
        if (!strcmp(token, "slit") || !strcmp(token, "absorb")) {
          Vars->Flag_Absorb = 1;  iskeyword=1; }
        if (!strcmp(token, "max"))  Set_Coord_Mode = DEFS->COORD_MAX;
        if (!strcmp(token, "bins") || !strcmp(token, "dim")) Set_Coord_Mode = DEFS->COORD_DIM;
        if (!strcmp(token, "file") || !strcmp(token, "filename")) {
          Set_Coord_Mode = DEFS->COORD_FIL;
          if (Flag_No) { strcpy(Vars->Mon_File,""); Vars->Coord_Number = 0; Flag_End = 1; }
        }
        if (!strcmp(token, "unactivate")) {
          Flag_End = 1; Vars->Coord_Number = 0; iskeyword=1; }
        if (!strcmp(token, "all"))    { Flag_All = 1;  iskeyword=1; }
        if (!strcmp(token, "sphere")) { Vars->Flag_Shape = DEFS->SHAPE_SPHERE; iskeyword=1; }
        if (!strcmp(token, "cylinder")) { Vars->Flag_Shape = DEFS->SHAPE_CYLIND; iskeyword=1; }
        if (!strcmp(token, "banana")) { Vars->Flag_Shape = DEFS->SHAPE_BANANA; iskeyword=1; }
        if (!strcmp(token, "square")) { Vars->Flag_Shape = DEFS->SHAPE_SQUARE; iskeyword=1; }
        if (!strcmp(token, "disk"))   { Vars->Flag_Shape = DEFS->SHAPE_DISK; iskeyword=1; }
        if (!strcmp(token, "box"))     { Vars->Flag_Shape = DEFS->SHAPE_BOX; iskeyword=1; }
        if (!strcmp(token, "previous")) { Vars->Flag_Shape = DEFS->SHAPE_PREVIOUS; iskeyword=1; }
        if (!strcmp(token, "parallel")){ Vars->Flag_parallel = 1; iskeyword=1; }
        if (!strcmp(token, "capture")) { Vars->Flag_capture = 1; iskeyword=1; }
        if (!strcmp(token, "auto") && (Flag_auto != -1)) {
          Vars->Flag_Auto_Limits = 1;
          if (Flag_All) Flag_auto = -1;
          else          Flag_auto = 1;
          iskeyword=1; Flag_All=0; }
        if (!strcmp(token, "premonitor")) {
          Vars->Flag_UsePreMonitor = 1; iskeyword=1; }
        if (!strcmp(token, "3He_pressure") || !strcmp(token, "pressure")) {
          Vars->He3_pressure = 3; iskeyword=1; }
        if (!strcmp(token, "no") || !strcmp(token, "not")) { Flag_No = 1;  iskeyword=1; }
        if (!strcmp(token, "signal")) Set_Coord_Mode = DEFS->COORD_SIGNAL;
        if (!strcmp(token, "mantid")) { Vars->Flag_mantid = 1; iskeyword=1; }

        /* Mode has changed: this was a keyword or value  ? */
        if (Set_Coord_Mode != old_Mode) iskeyword=1;

        /* now look for variable names to monitor */
        Set_Vars_Coord_Type = DEFS->COORD_NONE; lmin = 0; lmax = 0;

        if (!strcmp(token, "x"))
          { Set_Vars_Coord_Type = DEFS->COORD_X; strcpy(Set_Vars_Coord_Label,"x [m]"); strcpy(Set_Vars_Coord_Var,"x");
          lmin = Vars->mxmin; lmax = Vars->mxmax;
          Vars->Coord_Min[Vars->Coord_Number+1] = Vars->mxmin;
          Vars->Coord_Max[Vars->Coord_Number+1] = Vars->mxmax;}
        if (!strcmp(token, "y"))
          { Set_Vars_Coord_Type = DEFS->COORD_Y; strcpy(Set_Vars_Coord_Label,"y [m]"); strcpy(Set_Vars_Coord_Var,"y");
          lmin = Vars->mymin; lmax = Vars->mymax;
          Vars->Coord_Min[Vars->Coord_Number+1] = Vars->mymin;
          Vars->Coord_Max[Vars->Coord_Number+1] = Vars->mymax;}
        if (!strcmp(token, "z"))
          { Set_Vars_Coord_Type = DEFS->COORD_Z; strcpy(Set_Vars_Coord_Label,"z [m]"); strcpy(Set_Vars_Coord_Var,"z"); lmin = Vars->mzmin; lmax = Vars->mzmax; }
        if (!strcmp(token, "k") || !strcmp(token, "wavevector"))
          { Set_Vars_Coord_Type = DEFS->COORD_K; strcpy(Set_Vars_Coord_Label,"|k| [Angs-1]"); strcpy(Set_Vars_Coord_Var,"k"); lmin = 0; lmax = 10; }
        if (!strcmp(token, "v"))
          { Set_Vars_Coord_Type = DEFS->COORD_V; strcpy(Set_Vars_Coord_Label,"Velocity [m/s]"); strcpy(Set_Vars_Coord_Var,"v"); lmin = 0; lmax = 10000; }
        if (!strcmp(token, "t") || !strcmp(token, "time") || !strcmp(token, "tof"))
          { Set_Vars_Coord_Type = DEFS->COORD_T; strcpy(Set_Vars_Coord_Label,"TOF [s]"); strcpy(Set_Vars_Coord_Var,"t"); lmin = 0; lmax = .1; }
        if ((!strcmp(token, "p") || !strcmp(token, "i") || !strcmp(token, "intensity") || !strcmp(token, "flux")))
          { Set_Vars_Coord_Type = DEFS->COORD_P;
            strcpy(Set_Vars_Coord_Label,"Intensity");
            strncat(Set_Vars_Coord_Label, " [n/s", 30);
            if (Vars->Flag_per_cm2) strncat(Set_Vars_Coord_Label, "/cm2", 30);
            if (XY > 1 && Vars->Coord_Number)
              strncat(Set_Vars_Coord_Label, "/bin", 30);
            strncat(Set_Vars_Coord_Label, "]", 30);
            strcpy(Set_Vars_Coord_Var,"I");
            lmin = 0; lmax = FLT_MAX;
            if (Flag_auto>0) Flag_auto=0;
          }

        if (!strcmp(token, "vx"))
          { Set_Vars_Coord_Type = DEFS->COORD_VX; strcpy(Set_Vars_Coord_Label,"vx [m/s]"); strcpy(Set_Vars_Coord_Var,"vx"); lmin = -1000; lmax = 1000; }
        if (!strcmp(token, "vy"))
          { Set_Vars_Coord_Type = DEFS->COORD_VY; strcpy(Set_Vars_Coord_Label,"vy [m/s]"); strcpy(Set_Vars_Coord_Var,"vy"); lmin = -1000; lmax = 1000; }
        if (!strcmp(token, "vz"))
          { Set_Vars_Coord_Type = DEFS->COORD_VZ; strcpy(Set_Vars_Coord_Label,"vz [m/s]"); strcpy(Set_Vars_Coord_Var,"vz"); lmin = -10000; lmax = 10000; }
        if (!strcmp(token, "kx"))
          { Set_Vars_Coord_Type = DEFS->COORD_KX; strcpy(Set_Vars_Coord_Label,"kx [Angs-1]"); strcpy(Set_Vars_Coord_Var,"kx"); lmin = -1; lmax = 1; }
        if (!strcmp(token, "ky"))
          { Set_Vars_Coord_Type = DEFS->COORD_KY; strcpy(Set_Vars_Coord_Label,"ky [Angs-1]"); strcpy(Set_Vars_Coord_Var,"ky"); lmin = -1; lmax = 1; }
        if (!strcmp(token, "kz"))
          { Set_Vars_Coord_Type = DEFS->COORD_KZ; strcpy(Set_Vars_Coord_Label,"kz [Angs-1]"); strcpy(Set_Vars_Coord_Var,"kz"); lmin = -10; lmax = 10; }
        if (!strcmp(token, "sx"))
          { Set_Vars_Coord_Type = DEFS->COORD_SX; strcpy(Set_Vars_Coord_Label,"sx [1]"); strcpy(Set_Vars_Coord_Var,"sx"); lmin = -1; lmax = 1; }
        if (!strcmp(token, "sy"))
          { Set_Vars_Coord_Type = DEFS->COORD_SY; strcpy(Set_Vars_Coord_Label,"sy [1]"); strcpy(Set_Vars_Coord_Var,"sy"); lmin = -1; lmax = 1; }
        if (!strcmp(token, "sz"))
          { Set_Vars_Coord_Type = DEFS->COORD_SZ; strcpy(Set_Vars_Coord_Label,"sz [1]"); strcpy(Set_Vars_Coord_Var,"sz"); lmin = -1; lmax = 1; }

        if (!strcmp(token, "energy") || !strcmp(token, "omega") || !strcmp(token, "e"))
          { Set_Vars_Coord_Type = DEFS->COORD_ENERGY; strcpy(Set_Vars_Coord_Label,"Energy [meV]"); strcpy(Set_Vars_Coord_Var,"E"); lmin = 0; lmax = 100; }
        if (!strcmp(token, "lambda") || !strcmp(token, "wavelength") || !strcmp(token, "l"))
          { Set_Vars_Coord_Type = DEFS->COORD_LAMBDA; strcpy(Set_Vars_Coord_Label,"Wavelength [Angs]"); strcpy(Set_Vars_Coord_Var,"L"); lmin = 0; lmax = 100; }
        if (!strcmp(token, "radius") || !strcmp(token, "r"))
          { Set_Vars_Coord_Type = DEFS->COORD_RADIUS; strcpy(Set_Vars_Coord_Label,"Radius [m]"); strcpy(Set_Vars_Coord_Var,"xy"); lmin = 0; lmax = xmax; }
        if (!strcmp(token, "xy"))
          { Set_Vars_Coord_Type = DEFS->COORD_XY; strcpy(Set_Vars_Coord_Label,"Radius (xy) [m]"); strcpy(Set_Vars_Coord_Var,"xy"); lmin = 0; lmax = xmax; }
        if (!strcmp(token, "yz"))
          { Set_Vars_Coord_Type = DEFS->COORD_YZ; strcpy(Set_Vars_Coord_Label,"Radius (yz) [m]"); strcpy(Set_Vars_Coord_Var,"yz"); lmin = 0; lmax = xmax; }
        if (!strcmp(token, "xz"))
          { Set_Vars_Coord_Type = DEFS->COORD_XZ; strcpy(Set_Vars_Coord_Label,"Radius (xz) [m]"); strcpy(Set_Vars_Coord_Var,"xz"); lmin = 0; lmax = xmax; }
        if (!strcmp(token, "vxy"))
          { Set_Vars_Coord_Type = DEFS->COORD_VXY; strcpy(Set_Vars_Coord_Label,"Radial Velocity (xy) [m]"); strcpy(Set_Vars_Coord_Var,"Vxy"); lmin = 0; lmax = 2000; }
        if (!strcmp(token, "kxy"))
          { Set_Vars_Coord_Type = DEFS->COORD_KXY; strcpy(Set_Vars_Coord_Label,"Radial Wavevector (xy) [Angs-1]"); strcpy(Set_Vars_Coord_Var,"Kxy"); lmin = 0; lmax = 2; }
        if (!strcmp(token, "vyz"))
          { Set_Vars_Coord_Type = DEFS->COORD_VYZ; strcpy(Set_Vars_Coord_Label,"Radial Velocity (yz) [m]"); strcpy(Set_Vars_Coord_Var,"Vyz"); lmin = 0; lmax = 2000; }
        if (!strcmp(token, "kyz"))
          { Set_Vars_Coord_Type = DEFS->COORD_KYZ; strcpy(Set_Vars_Coord_Label,"Radial Wavevector (yz) [Angs-1]"); strcpy(Set_Vars_Coord_Var,"Kyz"); lmin = 0; lmax = 2; }
        if (!strcmp(token, "vxz"))
          { Set_Vars_Coord_Type = DEFS->COORD_VXZ; strcpy(Set_Vars_Coord_Label,"Radial Velocity (xz) [m]"); strcpy(Set_Vars_Coord_Var,"Vxz"); lmin = 0; lmax = 2000; }
        if (!strcmp(token, "kxz"))
          { Set_Vars_Coord_Type = DEFS->COORD_KXZ; strcpy(Set_Vars_Coord_Label,"Radial Wavevector (xz) [Angs-1]"); strcpy(Set_Vars_Coord_Var,"Kxz"); lmin = 0; lmax = 2; }
        if (!strcmp(token, "angle") || !strcmp(token, "a"))
          { Set_Vars_Coord_Type = DEFS->COORD_ANGLE; strcpy(Set_Vars_Coord_Label,"Angle [deg]"); strcpy(Set_Vars_Coord_Var,"A"); lmin = -50; lmax = 50; }
        if (!strcmp(token, "hdiv")|| !strcmp(token, "divergence") || !strcmp(token, "xdiv") || !strcmp(token, "hd") || !strcmp(token, "dx"))
          { Set_Vars_Coord_Type = DEFS->COORD_HDIV; strcpy(Set_Vars_Coord_Label,"Hor. Divergence [deg]"); strcpy(Set_Vars_Coord_Var,"hd"); lmin = -5; lmax = 5; }
        if (!strcmp(token, "vdiv") || !strcmp(token, "ydiv") || !strcmp(token, "vd") || !strcmp(token, "dy"))
          { Set_Vars_Coord_Type = DEFS->COORD_VDIV; strcpy(Set_Vars_Coord_Label,"Vert. Divergence [deg]"); strcpy(Set_Vars_Coord_Var,"vd"); lmin = -5; lmax = 5; }
        if (!strcmp(token, "theta") || !strcmp(token, "longitude") || !strcmp(token, "th"))
          { Set_Vars_Coord_Type = DEFS->COORD_THETA; strcpy(Set_Vars_Coord_Label,"Longitude [deg]"); strcpy(Set_Vars_Coord_Var,"th"); lmin = -180; lmax = 180; }
        if (!strcmp(token, "phi") || !strcmp(token, "lattitude") || !strcmp(token, "ph"))
          { Set_Vars_Coord_Type = DEFS->COORD_PHI; strcpy(Set_Vars_Coord_Label,"Lattitude [deg]"); strcpy(Set_Vars_Coord_Var,"ph"); lmin = -180; lmax = 180; }
        if (!strcmp(token, "ncounts") || !strcmp(token, "n") || !strcmp(token, "neutron"))
          { Set_Vars_Coord_Type = DEFS->COORD_NCOUNT; strcpy(Set_Vars_Coord_Label,"Neutron ID [1]"); strcpy(Set_Vars_Coord_Var,"n"); lmin = 0; lmax = mcget_ncount(); if (Flag_auto>0) Flag_auto=0; }
        if (!strcmp(token, "id") || !strcmp(token, "pixel"))
          { Set_Vars_Coord_Type = DEFS->COORD_PIXELID; 
            strcpy(Set_Vars_Coord_Label,"Pixel ID [1]"); 
            strcpy(Set_Vars_Coord_Var,"id"); lmin = 0; lmax = FLT_MAX; 
            if (Flag_auto>0) Flag_auto=0;
            Vars->Flag_List = 1; }
        if (!strcmp(token, "user") || !strcmp(token, "user1") || !strcmp(token, "u1"))
          { Set_Vars_Coord_Type = DEFS->COORD_USER1; strncpy(Set_Vars_Coord_Label,Vars->UserName1,30); strcpy(Set_Vars_Coord_Var,"U1"); lmin = -1e10; lmax = 1e10; }
        if (!strcmp(token, "user2") || !strcmp(token, "u2"))
          { Set_Vars_Coord_Type = DEFS->COORD_USER2; strncpy(Set_Vars_Coord_Label,Vars->UserName2,30); strcpy(Set_Vars_Coord_Var,"U2"); lmin = -1e10; lmax = 1e10; }
        if (!strcmp(token, "user3") || !strcmp(token, "u3"))
          { Set_Vars_Coord_Type = DEFS->COORD_USER3; strncpy(Set_Vars_Coord_Label,Vars->UserName3,30); strcpy(Set_Vars_Coord_Var,"U3"); lmin = -1e10; lmax = 1e10; }

        /* now stores variable keywords detected, if any */
        if (Set_Vars_Coord_Type != DEFS->COORD_NONE)
        {
          int Coord_Number = Vars->Coord_Number;
          if (Vars->Flag_log) { Set_Vars_Coord_Type |= DEFS->COORD_LOG; Vars->Flag_log = 0; }
          if (Flag_abs) { Set_Vars_Coord_Type |= DEFS->COORD_ABS; Flag_abs = 0; }
          if (Flag_auto != 0) { Set_Vars_Coord_Type |= DEFS->COORD_AUTO; 
            if (Flag_auto > 0) Flag_auto = 0; }
          if (Set_Coord_Mode == DEFS->COORD_SIGNAL)
          {
            Coord_Number = 0;
            Vars->Flag_signal = Set_Vars_Coord_Type;
          }
          else
          {
            if (Coord_Number < MONnD_COORD_NMAX)
            { Coord_Number++;
              Vars->Coord_Number = Coord_Number; 
              if (Set_Vars_Coord_Type != DEFS->COORD_PIXELID)
                Vars->Coord_NumberNoPixel++;
            }
            else if (Vars->Flag_Verbose) printf("Monitor_nD: %s reached max number of variables (%i).\n", Vars->compcurname, MONnD_COORD_NMAX);
          }
          Vars->Coord_Type[Coord_Number] = Set_Vars_Coord_Type;
          strncpy(Vars->Coord_Label[Coord_Number], Set_Vars_Coord_Label,30);
          strncpy(Vars->Coord_Var[Coord_Number], Set_Vars_Coord_Var,30);
          if (lmin > lmax) { XY = lmin; lmin=lmax; lmax = XY; }
          Vars->Coord_Min[Coord_Number] = lmin;
          Vars->Coord_Max[Coord_Number] = lmax;
          if (Set_Vars_Coord_Type == DEFS->COORD_NCOUNT || Set_Vars_Coord_Type == DEFS->COORD_PIXELID || Set_Vars_Coord_Type == DEFS->COORD_SIGNAL)
            Vars->Coord_Bin[Coord_Number] = 1;
          else
            Vars->Coord_Bin[Coord_Number] = 20;
          Set_Coord_Mode = DEFS->COORD_VAR;
          Flag_All = 0;
          Flag_No  = 0;
        } else {
          /* no variable name could be read from options */
          if (!iskeyword) {
            if (strcmp(token, "cm2") && strcmp(token, "incoming")
             && strcmp(token, "outgoing") && strcmp(token, "cm2")
             && strcmp(token, "cm^2") && strcmp(token, "float")
             && strcmp(token, "double") && strcmp(token, "binary")
             && strcmp(token, "steradian") && Vars->Flag_Verbose)
              printf("Monitor_nD: %s: unknown '%s' keyword in 'options'. Ignoring.\n", Vars->compcurname, token);
          }
        }
      carg++;
      } /* end if token */
    } /* end while carg */
    free(option_copy);
    if (carg == 128) printf("Monitor_nD: %s reached max number of tokens (%i). Skipping.\n", Vars->compcurname, 128);

    if ((Vars->Flag_Shape == DEFS->SHAPE_BOX) && (fabs(Vars->mzmax - Vars->mzmin) == 0)) Vars->Flag_Shape = DEFS->SHAPE_SQUARE;

    if (Vars->Flag_log == 1) Vars->Coord_Type[0] |= DEFS->COORD_LOG;
    if (Vars->Coord_Number == 0)
    { Vars->Flag_Auto_Limits=0; Vars->Flag_Multiple=0; Vars->Flag_List=0; }

    /* now setting Monitor Name from variable labels */
    strcpy(Vars->Monitor_Label,"");
    XY = 1; /* will contain total bin number */
    for (i = 0; i <= Vars->Coord_Number; i++)
    {
      if (Flag_auto != 0) Vars->Coord_Type[i] |= DEFS->COORD_AUTO;
      Set_Vars_Coord_Type = (Vars->Coord_Type[i] & (DEFS->COORD_LOG-1));
      if ((Set_Vars_Coord_Type == DEFS->COORD_X)
       || (Set_Vars_Coord_Type == DEFS->COORD_Y)
       || (Set_Vars_Coord_Type == DEFS->COORD_Z))
       strcpy(Short_Label[i],"Position");
      else
      if ((Set_Vars_Coord_Type == DEFS->COORD_THETA)
       || (Set_Vars_Coord_Type == DEFS->COORD_PHI)
       || (Set_Vars_Coord_Type == DEFS->COORD_ANGLE))
       strcpy(Short_Label[i],"Angle");
      else
      if ((Set_Vars_Coord_Type == DEFS->COORD_XY)
       || (Set_Vars_Coord_Type == DEFS->COORD_XZ)
       || (Set_Vars_Coord_Type == DEFS->COORD_YZ)
       || (Set_Vars_Coord_Type == DEFS->COORD_RADIUS))
       strcpy(Short_Label[i],"Radius");
      else
      if ((Set_Vars_Coord_Type == DEFS->COORD_VX)
       || (Set_Vars_Coord_Type == DEFS->COORD_VY)
       || (Set_Vars_Coord_Type == DEFS->COORD_VZ)
       || (Set_Vars_Coord_Type == DEFS->COORD_V)
       || (Set_Vars_Coord_Type == DEFS->COORD_VXY)
       || (Set_Vars_Coord_Type == DEFS->COORD_VYZ)
       || (Set_Vars_Coord_Type == DEFS->COORD_VXZ))
       strcpy(Short_Label[i],"Velocity");
      else
      if ((Set_Vars_Coord_Type == DEFS->COORD_KX)
       || (Set_Vars_Coord_Type == DEFS->COORD_KY)
       || (Set_Vars_Coord_Type == DEFS->COORD_KZ)
       || (Set_Vars_Coord_Type == DEFS->COORD_KXY)
       || (Set_Vars_Coord_Type == DEFS->COORD_KYZ)
       || (Set_Vars_Coord_Type == DEFS->COORD_KXZ)
       || (Set_Vars_Coord_Type == DEFS->COORD_K))
       strcpy(Short_Label[i],"Wavevector");
      else
      if ((Set_Vars_Coord_Type == DEFS->COORD_SX)
       || (Set_Vars_Coord_Type == DEFS->COORD_SY)
       || (Set_Vars_Coord_Type == DEFS->COORD_SZ))
       strcpy(Short_Label[i],"Spin");
      else
      if ((Set_Vars_Coord_Type == DEFS->COORD_HDIV)
       || (Set_Vars_Coord_Type == DEFS->COORD_VDIV))
       strcpy(Short_Label[i],"Divergence");
      else
      if (Set_Vars_Coord_Type == DEFS->COORD_ENERGY)
       strcpy(Short_Label[i],"Energy");
      else
      if (Set_Vars_Coord_Type == DEFS->COORD_LAMBDA)
       strcpy(Short_Label[i],"Wavelength");
      else
      if (Set_Vars_Coord_Type == DEFS->COORD_NCOUNT)
       strcpy(Short_Label[i],"Neutron_ID");
      else
      if (Set_Vars_Coord_Type == DEFS->COORD_PIXELID)
       strcpy(Short_Label[i],"Pixel_ID");
      else
      if (Set_Vars_Coord_Type == DEFS->COORD_T)
          strcpy(Short_Label[i],"Time_Of_Flight");
      else
      if (Set_Vars_Coord_Type == DEFS->COORD_P)
          strcpy(Short_Label[i],"Intensity");
      else
      if (Set_Vars_Coord_Type == DEFS->COORD_USER1)
          strncpy(Short_Label[i],Vars->UserName1,30);
      else
      if (Set_Vars_Coord_Type == DEFS->COORD_USER2)
          strncpy(Short_Label[i],Vars->UserName2,30);
      else
      if (Set_Vars_Coord_Type == DEFS->COORD_USER3)
          strncpy(Short_Label[i],Vars->UserName3,30);
      else
          strcpy(Short_Label[i],"Unknown");

      if (Vars->Coord_Type[i] & DEFS->COORD_ABS)
      { strcat(Vars->Coord_Label[i]," (abs)"); }

      if (Vars->Coord_Type[i] & DEFS->COORD_LOG)
      { strcat(Vars->Coord_Label[i]," (log)"); }

      strcat(Vars->Monitor_Label, " ");
      strcat(Vars->Monitor_Label, Short_Label[i]);
      XY *= Vars->Coord_Bin[i];

    } /* end for Short_Label */

    if ((Vars->Coord_Type[0] & (DEFS->COORD_LOG-1)) == DEFS->COORD_P) {
      strncat(Vars->Coord_Label[0], " [n/s", 30);
      if (Vars->Flag_per_cm2) strncat(Vars->Coord_Label[0], "/cm2", 30);

      if (XY > 1 && Vars->Coord_Number)
        strncat(Vars->Coord_Label[0], "/bin", 30);
      strncat(Vars->Coord_Label[0], "]", 30);
    }

    /* update label 'signal per bin' if more than 1 bin */
    if (XY > 1 && Vars->Coord_Number) {
      if (Vars->Flag_capture)
        printf("Monitor_nD: %s: Using capture flux weightening on %ld bins.\n"
               "WARNING     Use binned data with caution, and prefer monitor integral value (I,Ierr).\n", Vars->compcurname, (long)XY);
    }

    strcat(Vars->Monitor_Label, " Monitor");
    if (Vars->Flag_Shape == DEFS->SHAPE_SQUARE) strcat(Vars->Monitor_Label, " (Square)");
    if (Vars->Flag_Shape == DEFS->SHAPE_DISK)   strcat(Vars->Monitor_Label, " (Disk)");
    if (Vars->Flag_Shape == DEFS->SHAPE_SPHERE) strcat(Vars->Monitor_Label, " (Sphere)");
    if (Vars->Flag_Shape == DEFS->SHAPE_CYLIND) strcat(Vars->Monitor_Label, " (Cylinder)");
    if (Vars->Flag_Shape == DEFS->SHAPE_BANANA) strcat(Vars->Monitor_Label, " (Banana)");
    if (Vars->Flag_Shape == DEFS->SHAPE_BOX)    strcat(Vars->Monitor_Label, " (Box)");
    if (Vars->Flag_Shape == DEFS->SHAPE_PREVIOUS) strcat(Vars->Monitor_Label, " (on PREVIOUS)");
    if (Vars->Flag_Shape == DEFS->SHAPE_OFF) strcat(Vars->Monitor_Label, " (OFF geometry)");
    if ((Vars->Flag_Shape == DEFS->SHAPE_CYLIND) || (Vars->Flag_Shape == DEFS->SHAPE_BANANA) || (Vars->Flag_Shape == DEFS->SHAPE_SPHERE) || (Vars->Flag_Shape == DEFS->SHAPE_BOX))
    {
      if (strstr(Vars->option, "incoming"))
      {
        Vars->Flag_Shape = abs(Vars->Flag_Shape);
        strcat(Vars->Monitor_Label, " [in]");
      }
      else /* if strstr(Vars->option, "outgoing")) */
      {
        Vars->Flag_Shape = -abs(Vars->Flag_Shape);
        strcat(Vars->Monitor_Label, " [out]");
      }
    }
    if (Vars->Flag_UsePreMonitor == 1)
    {
        strcat(Vars->Monitor_Label, " at ");
        strncat(Vars->Monitor_Label, Vars->UserName1,30);
    }
    if (Vars->Flag_log == 1) strcat(Vars->Monitor_Label, " [log] ");

    /* now allocate memory to store variables in TRACE */

    /* Vars->Coord_Number  0   : intensity or signal
     * Vars->Coord_Number  1:n : detector variables */

    if ((Vars->Coord_NumberNoPixel != 2) && !Vars->Flag_Multiple && !Vars->Flag_List)
    { Vars->Flag_Multiple = 1; /* default is n1D */
      if (Vars->Coord_Number != Vars->Coord_NumberNoPixel) Vars->Flag_List = 1; }

    /* list and auto limits case : Vars->Flag_List or Vars->Flag_Auto_Limits
     * -> Buffer to flush and suppress after Vars->Flag_Auto_Limits
     */
    if ((Vars->Flag_Auto_Limits || Vars->Flag_List) && Vars->Coord_Number)
    { /* Dim : (Vars->Coord_Number+1)*Vars->Buffer_Block matrix (for p, dp) */
      Vars->Mon2D_Buffer = (double *)malloc((Vars->Coord_Number+1)*Vars->Buffer_Block*sizeof(double));
      if (Vars->Mon2D_Buffer == NULL)
      { printf("Monitor_nD: %s cannot allocate Vars->Mon2D_Buffer (%li). No list and auto limits.\n", Vars->compcurname, Vars->Buffer_Block*(Vars->Coord_Number+1)*sizeof(double)); Vars->Flag_List = 0; Vars->Flag_Auto_Limits = 0; }
      else
      {
        for (i=0; i < (Vars->Coord_Number+1)*Vars->Buffer_Block; Vars->Mon2D_Buffer[i++] = (double)0);
      }
      Vars->Buffer_Size = Vars->Buffer_Block;
    }

    /* 1D and n1D case : Vars->Flag_Multiple */
    if (Vars->Flag_Multiple && Vars->Coord_NumberNoPixel)
    { /* Dim : Vars->Coord_Number*Vars->Coord_Bin[i] vectors */
      Vars->Mon2D_N  = (double **)malloc((Vars->Coord_Number)*sizeof(double *));
      Vars->Mon2D_p  = (double **)malloc((Vars->Coord_Number)*sizeof(double *));
      Vars->Mon2D_p2 = (double **)malloc((Vars->Coord_Number)*sizeof(double *));
      if ((Vars->Mon2D_N == NULL) || (Vars->Mon2D_p == NULL) || (Vars->Mon2D_p2 == NULL))
      { fprintf(stderr,"Monitor_nD: %s n1D cannot allocate Vars->Mon2D_N/p/p2 (%li). Fatal.\n", Vars->compcurname, (Vars->Coord_Number)*sizeof(double *)); exit(-1); }
      for (i= 1; i <= Vars->Coord_Number; i++)
      {
        Vars->Mon2D_N[i-1]  = (double *)malloc(Vars->Coord_Bin[i]*sizeof(double));
        Vars->Mon2D_p[i-1]  = (double *)malloc(Vars->Coord_Bin[i]*sizeof(double));
        Vars->Mon2D_p2[i-1] = (double *)malloc(Vars->Coord_Bin[i]*sizeof(double));
        if ((Vars->Mon2D_N == NULL) || (Vars->Mon2D_p == NULL) || (Vars->Mon2D_p2 == NULL))
        { fprintf(stderr,"Monitor_nD: %s n1D cannot allocate %s Vars->Mon2D_N/p/p2[%li] (%li). Fatal.\n", Vars->compcurname, Vars->Coord_Var[i], i, (Vars->Coord_Bin[i])*sizeof(double *)); exit(-1); }
        else
        {
          for (j=0; j < Vars->Coord_Bin[i]; j++ )
          { Vars->Mon2D_N[i-1][j] = (double)0; Vars->Mon2D_p[i-1][j] = (double)0; Vars->Mon2D_p2[i-1][j] = (double)0; }
        }
      }
    }
    else /* 2D case : Vars->Coord_Number==2 and !Vars->Flag_Multiple and !Vars->Flag_List */
    if ((Vars->Coord_NumberNoPixel == 2) && !Vars->Flag_Multiple)
    { /* Dim : Vars->Coord_Bin[1]*Vars->Coord_Bin[2] matrix */
      Vars->Mon2D_N  = (double **)malloc((Vars->Coord_Bin[1])*sizeof(double *));
      Vars->Mon2D_p  = (double **)malloc((Vars->Coord_Bin[1])*sizeof(double *));
      Vars->Mon2D_p2 = (double **)malloc((Vars->Coord_Bin[1])*sizeof(double *));
      if ((Vars->Mon2D_N == NULL) || (Vars->Mon2D_p == NULL) || (Vars->Mon2D_p2 == NULL))
      { fprintf(stderr,"Monitor_nD: %s 2D cannot allocate %s Vars->Mon2D_N/p/p2 (%li). Fatal.\n", Vars->compcurname, Vars->Coord_Var[1], (Vars->Coord_Bin[1])*sizeof(double *)); exit(-1); }
      for (i= 0; i < Vars->Coord_Bin[1]; i++)
      {
        Vars->Mon2D_N[i]  = (double *)malloc(Vars->Coord_Bin[2]*sizeof(double));
        Vars->Mon2D_p[i]  = (double *)malloc(Vars->Coord_Bin[2]*sizeof(double));
        Vars->Mon2D_p2[i] = (double *)malloc(Vars->Coord_Bin[2]*sizeof(double));
        if ((Vars->Mon2D_N == NULL) || (Vars->Mon2D_p == NULL) || (Vars->Mon2D_p2 == NULL))
        { fprintf(stderr,"Monitor_nD: %s 2D cannot allocate %s Vars->Mon2D_N/p/p2[%li] (%li). Fatal.\n", Vars->compcurname, Vars->Coord_Var[1], i, (Vars->Coord_Bin[2])*sizeof(double *)); exit(-1); }
        else
        {
          for (j=0; j < Vars->Coord_Bin[2]; j++ )
          { Vars->Mon2D_N[i][j] = (double)0; Vars->Mon2D_p[i][j] = (double)0; Vars->Mon2D_p2[i][j] = (double)0; }
        }
      }
    }
    else {
      Vars->Mon2D_N = Vars->Mon2D_p = Vars->Mon2D_p2 = NULL;
    }
      /* no Mon2D allocated for
       * (Vars->Coord_Number != 2) && !Vars->Flag_Multiple && Vars->Flag_List */

    Vars->psum  = 0;
    Vars->p2sum = 0;
    Vars->Nsum  = 0;

    Vars->area  = fabs(Vars->mxmax - Vars->mxmin)*fabs(Vars->mymax - Vars->mymin)*1E4; /* in cm**2 for square and box shapes */
    Vars->Sphere_Radius = fabs(Vars->mxmax - Vars->mxmin)/2;
    if ((abs(Vars->Flag_Shape) == DEFS->SHAPE_DISK) || (abs(Vars->Flag_Shape) == DEFS->SHAPE_SPHERE))
    {
      Vars->area = PI*Vars->Sphere_Radius*Vars->Sphere_Radius*1E4; /* disk shapes */
    }


    if (Vars->area == 0 && abs(Vars->Flag_Shape) != DEFS->SHAPE_PREVIOUS ) {
      if (abs(Vars->Flag_Shape) != DEFS->SHAPE_OFF) {  
	Vars->Coord_Number = 0;
      }
    }
    if (Vars->Coord_Number == 0 && Vars->Flag_Verbose)
      printf("Monitor_nD: %s is unactivated (0D)\n", Vars->compcurname);
    Vars->Cylinder_Height = fabs(Vars->mymax - Vars->mymin);

    if (Vars->Flag_Verbose)
    {
      printf("Monitor_nD: %s is a %s.\n", Vars->compcurname, Vars->Monitor_Label);
      printf("Monitor_nD: version %s with options=%s\n", MONITOR_ND_LIB_H, Vars->option);
    }
    
    /* compute the product of bin dimensions for PixelID */
    Vars->Coord_BinProd[0]=1;
    for (i = 1; i <= Vars->Coord_Number; i++)
      Vars->Coord_BinProd[i]=Vars->Coord_Bin[i]*Vars->Coord_BinProd[i-1];
  } /* end Monitor_nD_Init */

/* ========================================================================= */
/* Monitor_nD_Trace: this routine is used to monitor one propagating neutron */
/* return values: 0=neutron was absorbed, -1=neutron was outside bounds, 1=neutron was measured*/
/* ========================================================================= */

int Monitor_nD_Trace(MonitornD_Defines_type *DEFS, MonitornD_Variables_type *Vars)
{

  double  XY=0, pp=0;
  int     retval;
  long    i =0, j =0;
  double  Coord[MONnD_COORD_NMAX];
  long    Coord_Index[MONnD_COORD_NMAX];
  char    While_End   =0;
  long    While_Buffer=0;
  char    Set_Vars_Coord_Type = DEFS->COORD_NONE;
  
  /* the logic below depends mainly on:
       Flag_List:        1=store 1 buffer, 2=list all, 3=re-use buffer 
       Flag_Auto_Limits: 0 (no auto limits/list), 1 (store events into Buffer), 2 (re-emit store events)
   */

  /* Vars->Flag_Auto_Limits=1: buffer full, we read the Buffer, and determine min and max bounds */
  if ((Vars->Buffer_Counter >= Vars->Buffer_Block) && (Vars->Flag_Auto_Limits == 1) && (Vars->Coord_Number > 0))
  {
    /* auto limits case : get limits in Buffer for each variable */
          /* Dim : (Vars->Coord_Number+1)*Vars->Buffer_Block matrix (for p, dp) */
    if (Vars->Flag_Verbose) printf("Monitor_nD: %s getting %li Auto Limits from List (%li events) in TRACE.\n", Vars->compcurname, Vars->Coord_Number, Vars->Buffer_Counter);
    for (i = 1; i <= Vars->Coord_Number; i++)
    {
      if (Vars->Coord_Type[i] & DEFS->COORD_AUTO)
      {
        Vars->Coord_Min[i] =  FLT_MAX;
        Vars->Coord_Max[i] = -FLT_MAX;
        for (j = 0; j < Vars->Buffer_Counter; j++)
        {
          XY = Vars->Mon2D_Buffer[i+j*(Vars->Coord_Number+1)];  /* scanning variables in Buffer */
          if (XY < Vars->Coord_Min[i]) Vars->Coord_Min[i] = XY;
          if (XY > Vars->Coord_Max[i]) Vars->Coord_Max[i] = XY;
        }
        if  (Vars->Flag_Verbose)  
          printf("  %s: min=%g max=%g\n", Vars->Coord_Var[i], Vars->Coord_Min[i], Vars->Coord_Max[i]);
      }
    }
    Vars->Flag_Auto_Limits = 2;  /* pass to 2nd auto limits step (read Buffer and generate new events to store in histograms) */
  } /* end if Flag_Auto_Limits == 1 */

  /* manage realloc for 'list all' if Buffer size exceeded: flush Buffer to file */
  if ((Vars->Buffer_Counter >= Vars->Buffer_Block) && (Vars->Flag_List >= 2))
  {
    if (Vars->Buffer_Size >= 1000000 || Vars->Flag_List == 3)
    { /* save current (possibly append) and re-use Buffer */
      Monitor_nD_Save(DEFS, Vars);
      Vars->Flag_List = 3;
      Vars->Buffer_Block = Vars->Buffer_Size;
      Vars->Buffer_Counter  = 0;
      Vars->Neutron_Counter = 0;
    }
    else
    {
      Vars->Mon2D_Buffer  = (double *)realloc(Vars->Mon2D_Buffer, (Vars->Coord_Number+1)*(Vars->Neutron_Counter+Vars->Buffer_Block)*sizeof(double));
      if (Vars->Mon2D_Buffer == NULL)
            { printf("Monitor_nD: %s cannot reallocate Vars->Mon2D_Buffer[%li] (%li). Skipping.\n", Vars->compcurname, i, (Vars->Neutron_Counter+Vars->Buffer_Block)*sizeof(double)); Vars->Flag_List = 1; }
      else { Vars->Buffer_Counter = 0; Vars->Buffer_Size = Vars->Neutron_Counter+Vars->Buffer_Block; }
    }
  } /* end if Buffer realloc */

  char    outsidebounds=0;
  while (!While_End)
  { /* we generate Coord[] and Coord_index[] from Buffer (auto limits) or passing neutron */
    if ((Vars->Flag_Auto_Limits == 2) && (Vars->Coord_Number > 0))
    { /* Vars->Flag_Auto_Limits == 2: read back from Buffer (Buffer is filled or auto limits have been computed) */
      if (While_Buffer < Vars->Buffer_Block)
      {
        /* first while loop (While_Buffer) */
        /* auto limits case : scan Buffer within limits and store in Mon2D */
        Coord[0] = pp = Vars->Mon2D_Buffer[While_Buffer*(Vars->Coord_Number+1)];

        for (i = 1; i <= Vars->Coord_Number; i++)
        {
          /* scanning variables in Buffer */
          if (Vars->Coord_Bin[i] <= 1) continue;
          XY = (Vars->Coord_Max[i]-Vars->Coord_Min[i]);

          Coord[i] = Vars->Mon2D_Buffer[i+While_Buffer*(Vars->Coord_Number+1)];
          if (XY > 0) Coord_Index[i] = floor((Coord[i]-Vars->Coord_Min[i])*Vars->Coord_Bin[i]/XY);
          else        Coord_Index[i] = 0;
          if (Vars->Flag_With_Borders)
          {
            if (Coord_Index[i] < 0)                   Coord_Index[i] = 0;
            if (Coord_Index[i] >= Vars->Coord_Bin[i]) Coord_Index[i] = Vars->Coord_Bin[i] - 1;
          }
        } /* end for */
        
        /* update the PixelID, we compute it from the previous variables index */
        if (Vars->Coord_NumberNoPixel < Vars->Coord_Number) /* there is a Pixel variable */
        for (i = 1; i <= Vars->Coord_Number; i++) {
          char Set_Vars_Coord_Type = (Vars->Coord_Type[i] & (DEFS->COORD_LOG-1));
          if (Set_Vars_Coord_Type == DEFS->COORD_PIXELID) {
            char flag_outside=0;
            Coord_Index[i] = Coord[i] = 0;
            for (j= 1; j < i; j++) {
              /* not for 1D variables with Bin=1 such as PixelID, NCOUNT, Intensity */
              if (Vars->Coord_Bin[j] == 1) continue; 
              if (0 > Coord_Index[j] || Coord_Index[j] >= Vars->Coord_Bin[j]) {
                flag_outside=1;
                Coord[i] = 0;
                break;
              }
              Coord[i] += Coord_Index[j]*Vars->Coord_BinProd[j-1];
            }
            if (!flag_outside) {
              Vars->Mon2D_Buffer[i+While_Buffer*(Vars->Coord_Number+1)] = Coord[i];
            }
          } /* end if PixelID */
        }
        While_Buffer++;
      } /* end if in Buffer */
      else /* (While_Buffer >= Vars->Buffer_Block) && (Vars->Flag_Auto_Limits == 2) */
      {
        Vars->Flag_Auto_Limits = 0;
        if (!Vars->Flag_List) /* free Buffer not needed anymore (no list to output) */
        { /* Dim : (Vars->Coord_Number+1)*Vars->Buffer_Block matrix (for p, p2) */
          free(Vars->Mon2D_Buffer); Vars->Mon2D_Buffer = NULL;
        }
        if (Vars->Flag_Verbose) printf("Monitor_nD: %s flushed %li Auto Limits from List (%li) in TRACE.\n", Vars->compcurname, Vars->Coord_Number, Vars->Buffer_Counter);
      }
    } /* if Vars->Flag_Auto_Limits == 2 */
    
    if (Vars->Flag_Auto_Limits != 2 || !Vars->Coord_Number) /* Vars->Flag_Auto_Limits == 0 (no auto limits/list) or 1 (store events into Buffer) */
    {
      /* automatically compute area and steradian solid angle when in AUTO mode */
      /* compute the steradian solid angle incoming on the monitor */
      double v;
      v=sqrt(Vars->cvx*Vars->cvx
            +Vars->cvy*Vars->cvy
            +Vars->cvz*Vars->cvz);
      if (Vars->min_x > Vars->cx) Vars->min_x = Vars->cx;
      if (Vars->max_x < Vars->cx) Vars->max_x = Vars->cx;
      if (Vars->min_y > Vars->cy) Vars->min_y = Vars->cy;
      if (Vars->max_y < Vars->cy) Vars->max_y = Vars->cy;
      Vars->mean_p  += Vars->cp;
      if (v) {
        Vars->mean_dx += Vars->cp*fabs(Vars->cvx/v);
        Vars->mean_dy += Vars->cp*fabs(Vars->cvy/v);
      }

      for (i = 0; i <= Vars->Coord_Number; i++)
      { /* handle current neutron : last while */
        XY = 0;
        Set_Vars_Coord_Type = (Vars->Coord_Type[i] & (DEFS->COORD_LOG-1));
        /* get values for variables to monitor */
        if (Set_Vars_Coord_Type == DEFS->COORD_X) XY = Vars->cx;
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_Y) XY = Vars->cy;
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_Z) XY = Vars->cz;
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_VX) XY = Vars->cvx;
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_VY) XY = Vars->cvy;
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_VZ) XY = Vars->cvz;
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_KX) XY = V2K*Vars->cvx;
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_KY) XY = V2K*Vars->cvy;
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_KZ) XY = V2K*Vars->cvz;
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_SX) XY = Vars->csx;
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_SY) XY = Vars->csy;
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_SZ) XY = Vars->csz;
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_T) XY = Vars->ct;
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_P) XY = Vars->cp;
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_HDIV) XY = RAD2DEG*atan2(Vars->cvx,Vars->cvz);
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_VDIV) XY = RAD2DEG*atan2(Vars->cvy,Vars->cvz);
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_V) XY = sqrt(Vars->cvx*Vars->cvx+Vars->cvy*Vars->cvy+Vars->cvz*Vars->cvz);
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_RADIUS)
          XY = sqrt(Vars->cx*Vars->cx+Vars->cy*Vars->cy+Vars->cz*Vars->cz);
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_XY)
          XY = sqrt(Vars->cx*Vars->cx+Vars->cy*Vars->cy)*(Vars->cx > 0 ? 1 : -1);
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_YZ) XY = sqrt(Vars->cy*Vars->cy+Vars->cz*Vars->cz);
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_XZ)
          XY = sqrt(Vars->cx*Vars->cx+Vars->cz*Vars->cz);
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_VXY) XY = sqrt(Vars->cvx*Vars->cvx+Vars->cvy*Vars->cvy);
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_VXZ) XY = sqrt(Vars->cvx*Vars->cvx+Vars->cvz*Vars->cvz);
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_VYZ) XY = sqrt(Vars->cvy*Vars->cvy+Vars->cvz*Vars->cvz);
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_K) { XY = sqrt(Vars->cvx*Vars->cvx+Vars->cvy*Vars->cvy+Vars->cvz*Vars->cvz);  XY *= V2K; }
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_KXY) { XY = sqrt(Vars->cvx*Vars->cvx+Vars->cvy*Vars->cvy);  XY *= V2K; }
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_KXZ) { XY = sqrt(Vars->cvx*Vars->cvx+Vars->cvz*Vars->cvz);  XY *= V2K; }
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_KYZ) { XY = sqrt(Vars->cvy*Vars->cvy+Vars->cvz*Vars->cvz);  XY *= V2K; }
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_ENERGY) { XY = Vars->cvx*Vars->cvx+Vars->cvy*Vars->cvy+Vars->cvz*Vars->cvz;  XY *= VS2E; }
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_LAMBDA) { XY = sqrt(Vars->cvx*Vars->cvx+Vars->cvy*Vars->cvy+Vars->cvz*Vars->cvz);  XY *= V2K; if (XY != 0) XY = 2*PI/XY; }
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_NCOUNT) XY = Vars->Neutron_Counter;
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_ANGLE)
        {  XY = sqrt(Vars->cvx*Vars->cvx+Vars->cvy*Vars->cvy);
           if (Vars->cvz != 0)
                XY = RAD2DEG*atan2(XY,Vars->cvz)*(Vars->cx > 0 ? 1 : -1);
           else XY = 0;
        }
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_THETA)  { if (Vars->cz != 0) XY = RAD2DEG*atan2(Vars->cx,Vars->cz); }
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_PHI) { if (Vars->cz != 0) XY = RAD2DEG*asin(Vars->cy/Vars->cz); }
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_USER1) XY = Vars->UserVariable1;
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_USER2) XY = Vars->UserVariable2;
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_USER3) XY = Vars->UserVariable3;
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_PIXELID && !Vars->Flag_Auto_Limits) {
          /* compute the PixelID from previous coordinates 
             the PixelID is the product of Coord_Index[i] in the detector geometry 
             pixelID = sum( Coord_Index[j]*prod(Vars->Coord_Bin[1:(j-1)]) )
             
             this does not apply when we store events in the buffer as Coord_Index
             is not set. Then the pixelID will be re-computed during SAVE.
          */
          char flag_outside=0;
          for (j= 1; j < i; j++) {
            /* not for 1D variables with Bin=1 such as PixelID, NCOUNT, Intensity */
            if (Vars->Coord_Bin[j] <= 1) continue; 
            if (0 > Coord_Index[j] || Coord_Index[j] >= Vars->Coord_Bin[j]) { 
              flag_outside=1; XY=0; break;
            }
            XY += Coord_Index[j]*Vars->Coord_BinProd[j-1];
          }
	  if (Vars->Flag_mantid && Vars->Flag_OFF && Vars->OFF_polyidx >=0) XY=Vars->OFF_polyidx;
          if (!flag_outside) XY += Vars->Coord_Min[i];
        }
        
        /* handle 'abs' and 'log' keywords */
        if (Vars->Coord_Type[i] & DEFS->COORD_ABS) XY=fabs(XY);

        if (Vars->Coord_Type[i] & DEFS->COORD_LOG) /* compute log of variable if requested */
        {  if (XY > 0) XY = log(XY)/log(10);
           else        XY = -100; }

        Coord[i] = XY; Coord_Index[i] = 0;
        if (i == 0) { pp = XY; Coord_Index[i] = 0; }
        else {
        /* check bounds for variables which have no automatic limits */
          if ((!Vars->Flag_Auto_Limits || !(Vars->Coord_Type[i] & DEFS->COORD_AUTO)) && Vars->Coord_Bin[i]>1)
          { /* compute index in histograms for each variable to monitor */
            XY = (Vars->Coord_Max[i]-Vars->Coord_Min[i]);
            if (XY > 0) Coord_Index[i] = floor((Coord[i]-Vars->Coord_Min[i])*Vars->Coord_Bin[i]/XY);
            if (Vars->Flag_With_Borders)
            {
              if (Coord_Index[i] >= Vars->Coord_Bin[i]) Coord_Index[i] = Vars->Coord_Bin[i] - 1;
              if (Coord_Index[i] < 0) Coord_Index[i] = 0;
            }
            //if (0 > Coord_Index[i] || Coord_Index[i] >= Vars->Coord_Bin[i])
            //  outsidebounds=1;
          } /* else will get Index later from Buffer when Flag_Auto_Limits == 2 */
        }
        
      } /* end for i */
      While_End = 1;
    }/* end else if Vars->Flag_Auto_Limits == 2 */
    
    /* ====================================================================== */
    /* store n1d/2d neutron from Buffer (Auto_Limits == 2) or current neutron in while */
    if (Vars->Flag_Auto_Limits != 1) /* not when storing auto limits Buffer */
    {
      /* apply per cm2 */
      if (Vars->Flag_per_cm2 && Vars->area != 0)
        pp /= Vars->area;

      /* 2D case : Vars->Coord_Number==2 and !Vars->Flag_Multiple and !Vars->Flag_List */
      if ( Vars->Coord_NumberNoPixel == 2 && !Vars->Flag_Multiple)
      { /* Dim : Vars->Coord_Bin[1]*Vars->Coord_Bin[2] matrix */
        
        i = Coord_Index[1];
        j = Coord_Index[2];
        if (i >= 0 && i < Vars->Coord_Bin[1] && j >= 0 && j < Vars->Coord_Bin[2])
        {
          if (Vars->Mon2D_N) { 
            Vars->Mon2D_N[i][j]++;
            Vars->Mon2D_p[i][j] += pp;
            Vars->Mon2D_p2[i][j] += pp*pp;
          }
        } else {
          outsidebounds=1; 
        }
      } else {
        /* 1D and n1D case : Vars->Flag_Multiple */
        /* Dim : Vars->Coord_Number*Vars->Coord_Bin[i] vectors (intensity is not included) */
          
        for (i= 1; i <= Vars->Coord_Number; i++) {
          j = Coord_Index[i];
          if (j >= 0 && j < Vars->Coord_Bin[i]) {
            if  (Vars->Flag_Multiple && Vars->Mon2D_N) {
              Vars->Mon2D_N[i-1][j]++;
              Vars->Mon2D_p[i-1][j]  += pp;
              Vars->Mon2D_p2[i-1][j] += pp*pp;
            }
          } else { 
            outsidebounds=1;
            break;
          }
        }
      }
    } /* end (Vars->Flag_Auto_Limits != 1) */
    
    if (Vars->Flag_Auto_Limits != 2 && !outsidebounds) /* not when reading auto limits Buffer */
    { /* now store Coord into Buffer (no index needed) if necessary (list or auto limits) */
      if ((Vars->Buffer_Counter < Vars->Buffer_Block) && ((Vars->Flag_List) || (Vars->Flag_Auto_Limits == 1)))
      {
          
        for (i = 0; i <= Vars->Coord_Number; i++)
        {
          Vars->Mon2D_Buffer[i + Vars->Neutron_Counter*(Vars->Coord_Number+1)] = Coord[i];
        }
        Vars->Buffer_Counter++;
        if (Vars->Flag_Verbose && (Vars->Buffer_Counter >= Vars->Buffer_Block) && (Vars->Flag_List == 1)) 
          printf("Monitor_nD: %s %li neutrons stored in List.\n", Vars->compcurname, Vars->Buffer_Counter);
      }
      Vars->Neutron_Counter++;
    } /* end (Vars->Flag_Auto_Limits != 2) */
    
  } /* end while */
  Vars->Nsum++;
  Vars->psum  += pp;
  Vars->p2sum += pp*pp;

  /*determine return value: 1:neutron was in bounds and measured, -1: outside bounds, 0: outside bounds, should be absorbed.*/
  if(outsidebounds){
      if(Vars->Flag_Absorb){
          return 0;
      }else{
          return -1;
      }
  }
  return 1;
} /* end Monitor_nD_Trace */

/* ========================================================================= */
/* Monitor_nD_Save: this routine is used to save data files                  */
/* ========================================================================= */

MCDETECTOR Monitor_nD_Save(MonitornD_Defines_type *DEFS, MonitornD_Variables_type *Vars)
  {
    char   *fname;
    long    i,j;
    double *p0m = NULL;
    double *p1m = NULL;
    double *p2m = NULL;
    char    Coord_X_Label[CHAR_BUF_LENGTH];
    double  min1d, max1d;
    double  min2d, max2d;
    long    bin1d, bin2d;
    char    While_End = 0;
    long    While_Buffer = 0;
    double  XY=0, pp=0;
    double  Coord[MONnD_COORD_NMAX];
    long    Coord_Index[MONnD_COORD_NMAX];
    char    label[CHAR_BUF_LENGTH];
    double  ratio;

    MCDETECTOR detector;

    ratio = 100.0*mcget_run_num()/mcget_ncount();
    if (Vars->Flag_Verbose && Vars->Flag_per_cm2) {
      printf("Monitor_nD: %s: active flat detector area is %g [cm^2], total area is %g [cm^2]\n",
        Vars->compcurname, (Vars->max_x-Vars->min_x)
                          *(Vars->max_y-Vars->min_y)*1E4, Vars->area);
      printf("Monitor_nD: %s: beam solid angle is %g [st] (%g x %g [deg^2])\n",
        Vars->compcurname,
        2*fabs(2*atan(Vars->mean_dx/Vars->mean_p)
         *sin(2*atan(Vars->mean_dy/Vars->mean_p)/2)),
        atan(Vars->mean_dx/Vars->mean_p)*RAD2DEG,
        atan(Vars->mean_dy/Vars->mean_p)*RAD2DEG);
    }

    /* check Buffer flush when end of simulation reached */
    if ((Vars->Buffer_Counter <= Vars->Buffer_Block) && Vars->Flag_Auto_Limits && Vars->Mon2D_Buffer && Vars->Buffer_Counter)
    {
      /* Get Auto Limits */
      if (Vars->Flag_Verbose) printf("Monitor_nD: %s getting %li Auto Limits from List (%li events).\n", Vars->compcurname, Vars->Coord_Number, Vars->Buffer_Counter);

      for (i = 1; i <= Vars->Coord_Number; i++)
      {
        if ((Vars->Coord_Type[i] & DEFS->COORD_AUTO) && Vars->Coord_Bin[i] > 1)
        {
          Vars->Coord_Min[i] = FLT_MAX;
          Vars->Coord_Max[i] = -FLT_MAX;
          for (j = 0; j < Vars->Buffer_Counter; j++)
          {
            XY = Vars->Mon2D_Buffer[i+j*(Vars->Coord_Number+1)];  /* scanning variables in Buffer */
            if (XY < Vars->Coord_Min[i]) Vars->Coord_Min[i] = XY;
            if (XY > Vars->Coord_Max[i]) Vars->Coord_Max[i] = XY;
          }
          if  (Vars->Flag_Verbose)  
            printf("  %s: min=%g max=%g in %li bins\n", Vars->Coord_Var[i], Vars->Coord_Min[i], Vars->Coord_Max[i], Vars->Coord_Bin[i]);
        }
      }
      Vars->Flag_Auto_Limits = 2;  /* pass to 2nd auto limits step */
      Vars->Buffer_Block = Vars->Buffer_Counter;

      while (!While_End)
      { /* we generate Coord[] and Coord_index[] from Buffer (auto limits) */
        /* simulation ended before Buffer was filled. Limits have to be computed, and stored events must be sent into histograms */
        
        if (While_Buffer < Vars->Buffer_Block)
        {
          /* first while loops (While_Buffer) */
          Coord[0] = Vars->Mon2D_Buffer[While_Buffer*(Vars->Coord_Number+1)];

          /* auto limits case : scan Buffer within limits and store in Mon2D */
          for (i = 1; i <= Vars->Coord_Number; i++)
          {
            /* scanning variables in Buffer */
            if (Vars->Coord_Bin[i] <= 1) Coord_Index[i] = 0;
            else {
              XY = (Vars->Coord_Max[i]-Vars->Coord_Min[i]);
              Coord[i] = Vars->Mon2D_Buffer[i+While_Buffer*(Vars->Coord_Number+1)];
              if (XY > 0) Coord_Index[i] = floor((Coord[i]-Vars->Coord_Min[i])*Vars->Coord_Bin[i]/XY);
              else Coord_Index[i] = 0;
              if (Vars->Flag_With_Borders)
              {
                if (Coord_Index[i] < 0) Coord_Index[i] = 0;
                if (Coord_Index[i] >= Vars->Coord_Bin[i]) Coord_Index[i] = Vars->Coord_Bin[i] - 1;
              }
            }
          } /* end for */

          /* update the PixelID, we compute it from the previous variables index */
          for (i = 1; i <= Vars->Coord_Number; i++) {
            char Set_Vars_Coord_Type = (Vars->Coord_Type[i] & (DEFS->COORD_LOG-1));
            if (Set_Vars_Coord_Type == DEFS->COORD_PIXELID) {
              char outsidebounds=0;
              Coord_Index[i] = Coord[i] = 0;
              for (j= 1; j < i; j++) {
                /* not for 1D variables with Bin=1 such as PixelID, NCOUNT, Intensity */
                if (Vars->Coord_Bin[j] == 1) continue; 
                if (0 > Coord_Index[j] || Coord_Index[j] >= Vars->Coord_Bin[j]) {
                  outsidebounds=1;
                  Coord[i] = 0;
                  break;
                }
                Coord[i] += Coord_Index[j]*Vars->Coord_BinProd[j-1];
              }
              if (!outsidebounds) {
                Vars->Mon2D_Buffer[i+While_Buffer*(Vars->Coord_Number+1)] = Coord[i];
              }
            } /* end if PixelID */
          }
          While_Buffer++;
        } /* end if in Buffer */
        else /* (While_Buffer >= Vars->Buffer_Block) && (Vars->Flag_Auto_Limits == 2) */
        {
          Vars->Flag_Auto_Limits = 0;
          While_End = 1;
          if (Vars->Flag_Verbose) printf("Monitor_nD: %s flushed %li Auto Limits from List (%li).\n", Vars->compcurname, Vars->Coord_Number, Vars->Buffer_Counter);
        }

        /* store n1d/2d section from Buffer */

        pp = Coord[0];
        /* apply per cm2 or per st */
        if (Vars->Flag_per_cm2 && Vars->area      != 0)
          pp /= Vars->area;
        
        /* 2D case : Vars->Coord_Number==2 and !Vars->Flag_Multiple and !Vars->Flag_List */
        if (!Vars->Flag_Multiple && Vars->Coord_NumberNoPixel == 2)
        { /* Dim : Vars->Coord_Bin[1]*Vars->Coord_Bin[2] matrix */
          i = Coord_Index[1];
          j = Coord_Index[2];
          if (i >= 0 && i < Vars->Coord_Bin[1] && j >= 0 && j < Vars->Coord_Bin[2])
          {
            if (Vars->Mon2D_N) {
              Vars->Mon2D_N[i][j]++;
              Vars->Mon2D_p[i][j] += pp;
              Vars->Mon2D_p2[i][j] += pp*pp;
            }
          } else if (Vars->Flag_Absorb) pp=0;
        }
        else
        /* 1D and n1D case : Vars->Flag_Multiple */
        { /* Dim : Vars->Coord_Number*Vars->Coord_Bin[i] vectors (intensity is not included) */
          for (i= 1; i <= Vars->Coord_Number; i++)
          {
            j = Coord_Index[i];
            if (j >= 0 && j < Vars->Coord_Bin[i])
            {
              if (Vars->Flag_Multiple && Vars->Mon2D_N) {
                Vars->Mon2D_N[i-1][j]++;
                Vars->Mon2D_p[i-1][j] += pp;
                Vars->Mon2D_p2[i-1][j] += pp*pp;
              }
            } else if (Vars->Flag_Absorb) {
              pp=0; break;
            }
          }
        } /* end store 2D/1D */
        
      } /* end while */
    } /* end Force Get Limits */

    /* write output files (sent to file as p[i*n + j] vectors) */
    if (Vars->Coord_Number == 0)
    {
      double Nsum;
      double psum, p2sum;
      Nsum = Vars->Nsum;
      psum = Vars->psum;
      p2sum= Vars->p2sum;
      if (Vars->Flag_signal != DEFS->COORD_P && Nsum > 0)
      { psum /=Nsum; p2sum /= Nsum*Nsum; }
      /* DETECTOR_OUT_0D(Vars->Monitor_Label, Vars->Nsum, Vars->psum, Vars->p2sum); */
      detector = mcdetector_out_0D(Vars->Monitor_Label, Nsum, psum, p2sum, Vars->compcurname, Vars->compcurpos);
    }
    else
    if (strlen(Vars->Mon_File) > 0)
    {
      fname = (char*)malloc(strlen(Vars->Mon_File)+10*Vars->Coord_Number);
      if (Vars->Flag_List && Vars->Mon2D_Buffer) /* List: DETECTOR_OUT_2D */
      {
       
        if (Vars->Flag_List >= 2) Vars->Buffer_Size = Vars->Neutron_Counter;
        if (Vars->Buffer_Size >= Vars->Neutron_Counter)
          Vars->Buffer_Size = Vars->Neutron_Counter;
        strcpy(fname,Vars->Mon_File);
        if (strchr(Vars->Mon_File,'.') == NULL) strcat(fname, "_list");

        strcpy(Coord_X_Label,"");
        for (i= 0; i <= Vars->Coord_Number; i++)
        {
          strcat(Coord_X_Label, Vars->Coord_Var[i]);
          strcat(Coord_X_Label, " ");
          if (strchr(Vars->Mon_File,'.') == NULL)
          { strcat(fname, "."); strcat(fname, Vars->Coord_Var[i]); }
        }
        if (Vars->Flag_Verbose) printf("Monitor_nD: %s write monitor file %s List (%lix%li).\n", Vars->compcurname, fname,bin2d,bin1d);

        /* handle the type of list output */
        strcpy(label, Vars->Monitor_Label);
        
        detector = mcdetector_out_list(
              label, "List of neutron events", Coord_X_Label,
              -Vars->Buffer_Size, Vars->Coord_Number+1,
              Vars->Mon2D_Buffer,
              fname, Vars->compcurname, Vars->compcurpos);
      }
      if (Vars->Flag_Multiple) /* n1D: DETECTOR_OUT_1D */
      {
        for (i= 0; i < Vars->Coord_Number; i++)
        {

          strcpy(fname,Vars->Mon_File);
          if (strchr(Vars->Mon_File,'.') == NULL)
          { strcat(fname, "."); strcat(fname, Vars->Coord_Var[i+1]); }
          sprintf(Coord_X_Label, "%s monitor", Vars->Coord_Label[i+1]);
          strcpy(label, Coord_X_Label);
          if (Vars->Coord_Bin[i+1] > 0) { /* 1D monitor */
            if (Vars->Flag_Verbose) printf("Monitor_nD: %s write monitor file %s 1D (%li).\n", Vars->compcurname, fname, Vars->Coord_Bin[i+1]);
            min1d = Vars->Coord_Min[i+1];
            max1d = Vars->Coord_Max[i+1];
            if (min1d == max1d) max1d = min1d+1e-6;
            p1m = (double *)malloc(Vars->Coord_Bin[i+1]*sizeof(double));
            p2m = (double *)malloc(Vars->Coord_Bin[i+1]*sizeof(double));
            if (p2m == NULL) /* use Raw Buffer line output */
            {
              if (Vars->Flag_Verbose) printf("Monitor_nD: %s cannot allocate memory for output. Using raw data.\n", Vars->compcurname);
              if (p1m != NULL) free(p1m);
              detector = mcdetector_out_1D(
              label,
              Vars->Coord_Label[i+1],
              Vars->Coord_Label[0],
              Vars->Coord_Var[i+1],
              min1d, max1d,
              Vars->Coord_Bin[i+1],
              Vars->Mon2D_N[i],Vars->Mon2D_p[i],Vars->Mon2D_p2[i],
              fname, Vars->compcurname, Vars->compcurpos);
            } /* if (p2m == NULL) */
            else
            {
              if (Vars->Flag_log != 0)
              {
                XY = FLT_MAX;
                for (j=0; j < Vars->Coord_Bin[i+1]; j++) /* search min of signal */
                  if ((XY > Vars->Mon2D_p[i][j]) && (Vars->Mon2D_p[i][j] > 0)) XY = Vars->Mon2D_p[i][j];
                if (XY <= 0) XY = -log(FLT_MAX)/log(10); else XY = log(XY)/log(10)-1;
              } /* if */

              for (j=0; j < Vars->Coord_Bin[i+1]; j++)
              {
                p1m[j] = Vars->Mon2D_p[i][j];
                p2m[j] = Vars->Mon2D_p2[i][j];
                if (Vars->Flag_signal != DEFS->COORD_P && Vars->Mon2D_N[i][j] > 0)
                { /* normalize mean signal to the number of events */
                  p1m[j] /= Vars->Mon2D_N[i][j];
                  p2m[j] /= Vars->Mon2D_N[i][j]*Vars->Mon2D_N[i][j];
                }
                if (Vars->Flag_log != 0)
                {
                  if ((p1m[j] > 0) && (p2m[j] > 0))
                  {
                    p2m[j] /= p1m[j]*p1m[j];
                    p1m[j] = log(p1m[j])/log(10);
                  }
                  else
                  {
                    p1m[j] = XY;
                    p2m[j] = 0;
                  }
                }
              } /* for */
              detector = mcdetector_out_1D(
                label,
                Vars->Coord_Label[i+1],
                Vars->Coord_Label[0],
                Vars->Coord_Var[i+1],
                min1d, max1d,
                Vars->Coord_Bin[i+1],
                Vars->Mon2D_N[i],p1m,p2m,
                fname, Vars->compcurname, Vars->compcurpos);

            } /* else */
            /* comment out 'free memory' lines to avoid loosing arrays if
               'detector' structure is used by other instrument parts
            if (p1m != NULL) free(p1m); p1m=NULL;
            if (p2m != NULL) free(p2m); p2m=NULL;
            */
          } else { /* 0d monitor */
            detector = mcdetector_out_0D(label, Vars->Mon2D_p[i][0], Vars->Mon2D_p2[i][0], Vars->Mon2D_N[i][0], Vars->compcurname, Vars->compcurpos);
          }


        } /* for */
      } /* if 1D */
      else
      if (Vars->Coord_NumberNoPixel == 2)  /* 2D: DETECTOR_OUT_2D */
      {
        strcpy(fname,Vars->Mon_File);

        p0m = (double *)malloc(Vars->Coord_Bin[1]*Vars->Coord_Bin[2]*sizeof(double));
        p1m = (double *)malloc(Vars->Coord_Bin[1]*Vars->Coord_Bin[2]*sizeof(double));
        p2m = (double *)malloc(Vars->Coord_Bin[1]*Vars->Coord_Bin[2]*sizeof(double));
        if (p2m == NULL)
        {
          if (Vars->Flag_Verbose) printf("Monitor_nD: %s cannot allocate memory for 2D array (%li). Skipping.\n", Vars->compcurname, 3*Vars->Coord_Bin[1]*Vars->Coord_Bin[2]*sizeof(double));
          /* comment out 'free memory' lines to avoid loosing arrays if
               'detector' structure is used by other instrument parts
          if (p0m != NULL) free(p0m);
          if (p1m != NULL) free(p1m);
          */
        }
        else
        {
          if (Vars->Flag_log != 0)
          {
            XY = FLT_MAX;
            for (i= 0; i < Vars->Coord_Bin[1]; i++)
              for (j= 0; j < Vars->Coord_Bin[2]; j++) /* search min of signal */
                if ((XY > Vars->Mon2D_p[i][j]) && (Vars->Mon2D_p[i][j]>0)) XY = Vars->Mon2D_p[i][j];
            if (XY <= 0) XY = -log(FLT_MAX)/log(10); else XY = log(XY)/log(10)-1;
          }
          for (i= 0; i < Vars->Coord_Bin[1]; i++)
          {
            for (j= 0; j < Vars->Coord_Bin[2]; j++)
            {
              long index;
              index = j + i*Vars->Coord_Bin[2];
              p0m[index] = Vars->Mon2D_N[i][j];
              p1m[index] = Vars->Mon2D_p[i][j];
              p2m[index] = Vars->Mon2D_p2[i][j];
              if (Vars->Flag_signal != DEFS->COORD_P && p0m[index] > 0)
              {
                  p1m[index] /= p0m[index];
                  p2m[index] /= p0m[index]*p0m[index];
              }

              if (Vars->Flag_log != 0)
              {
                if ((p1m[index] > 0) && (p2m[index] > 0))
                {
                  p2m[index] /= (p1m[index]*p1m[index]);
                  p1m[index] = log(p1m[index])/log(10);

                }
                else
                {
                  p1m[index] = XY;
                  p2m[index] = 0;
                }
              }
            }
          }
          if (strchr(Vars->Mon_File,'.') == NULL)
          { strcat(fname, "."); strcat(fname, Vars->Coord_Var[1]);
              strcat(fname, "_"); strcat(fname, Vars->Coord_Var[2]); }
          if (Vars->Flag_Verbose) printf("Monitor_nD: %s write monitor file %s 2D (%lix%li).\n", Vars->compcurname, fname, Vars->Coord_Bin[1], Vars->Coord_Bin[2]);

          min1d = Vars->Coord_Min[1];
          max1d = Vars->Coord_Max[1];
          if (min1d == max1d) max1d = min1d+1e-6;
          min2d = Vars->Coord_Min[2];
          max2d = Vars->Coord_Max[2];
          if (min2d == max2d) max2d = min2d+1e-6;
          strcpy(label, Vars->Monitor_Label);
          if (Vars->Coord_Bin[1]*Vars->Coord_Bin[2] > 1
           && Vars->Flag_signal == DEFS->COORD_P)
            strcat(label, " per bin");

          detector = mcdetector_out_2D(
            label,
            Vars->Coord_Label[1],
            Vars->Coord_Label[2],
            min1d, max1d,
            min2d, max2d,
            Vars->Coord_Bin[1],
            Vars->Coord_Bin[2],
            p0m,p1m,p2m,
            fname, Vars->compcurname, Vars->compcurpos);

          /* comment out 'free memory' lines to avoid loosing arrays if
               'detector' structure is used by other instrument parts
          if (p0m != NULL) free(p0m);
          if (p1m != NULL) free(p1m);
          if (p2m != NULL) free(p2m);
          */
        }
      }
      free(fname);
    }
    return(detector);
  } /* end Monitor_nD_Save */

/* ========================================================================= */
/* Monitor_nD_Finally: this routine is used to free memory                   */
/* ========================================================================= */

void Monitor_nD_Finally(MonitornD_Defines_type *DEFS,
  MonitornD_Variables_type *Vars)
  {
    int i;

    /* Now Free memory Mon2D.. */
    if ((Vars->Flag_Auto_Limits || Vars->Flag_List) && Vars->Coord_Number)
    { /* Dim : (Vars->Coord_Number+1)*Vars->Buffer_Block matrix (for p, dp) */
      if (Vars->Mon2D_Buffer != NULL) free(Vars->Mon2D_Buffer);
    }

    /* 1D and n1D case : Vars->Flag_Multiple */
    if (Vars->Flag_Multiple && Vars->Coord_Number)
    { /* Dim : Vars->Coord_Number*Vars->Coord_Bin[i] vectors */
      for (i= 0; i < Vars->Coord_Number; i++)
      {
        free(Vars->Mon2D_N[i]);
        free(Vars->Mon2D_p[i]);
        free(Vars->Mon2D_p2[i]);
      }
      free(Vars->Mon2D_N);
      free(Vars->Mon2D_p);
      free(Vars->Mon2D_p2);
    }


    /* 2D case : Vars->Coord_Number==2 and !Vars->Flag_Multiple and !Vars->Flag_List */
    if ((Vars->Coord_NumberNoPixel == 2) && !Vars->Flag_Multiple)
    { /* Dim : Vars->Coord_Bin[1]*Vars->Coord_Bin[2] matrix */
      for (i= 0; i < Vars->Coord_Bin[1]; i++)
      {
        free(Vars->Mon2D_N[i]);
        free(Vars->Mon2D_p[i]);
        free(Vars->Mon2D_p2[i]);
      }
      free(Vars->Mon2D_N);
      free(Vars->Mon2D_p);
      free(Vars->Mon2D_p2);
    }
  } /* end Monitor_nD_Finally */

/* ========================================================================= */
/* Monitor_nD_McDisplay: this routine is used to display component           */
/* ========================================================================= */

void Monitor_nD_McDisplay(MonitornD_Defines_type *DEFS,
  MonitornD_Variables_type *Vars)
  {
    double radius, h;
    double xmin;
    double xmax;
    double ymin;
    double ymax;
    double zmin;
    double zmax;
    int    i;
    double hdiv_min=-180, hdiv_max=180, vdiv_min=-90, vdiv_max=90;
    char   restricted = 0;

    radius = Vars->Sphere_Radius;
    h = Vars->Cylinder_Height;
    xmin = Vars->mxmin;
    xmax = Vars->mxmax;
    ymin = Vars->mymin;
    ymax = Vars->mymax;
    zmin = Vars->mzmin;
    zmax = Vars->mzmax;

    /* determine if there are angular limits set at start (no auto) in coord_types
     * cylinder/banana: look for hdiv
     * sphere: look for angle, radius (->atan2(val,radius)), hdiv, vdiv
     * this activates a 'restricted' flag, to draw a region as blades on cylinder/sphere
     */
    for (i= 0; i <= Vars->Coord_Number; i++)
    {
      int Set_Vars_Coord_Type;
      Set_Vars_Coord_Type = (Vars->Coord_Type[i] & (DEFS->COORD_LOG-1));
      if (Set_Vars_Coord_Type == DEFS->COORD_HDIV || Set_Vars_Coord_Type == DEFS->COORD_THETA)
      { hdiv_min = Vars->Coord_Min[i]; hdiv_max = Vars->Coord_Max[i]; restricted = 1; }
      else if (Set_Vars_Coord_Type == DEFS->COORD_VDIV || Set_Vars_Coord_Type == DEFS->COORD_PHI)
      { vdiv_min = Vars->Coord_Min[i]; vdiv_max = Vars->Coord_Max[i];restricted = 1;  }
      else if (Set_Vars_Coord_Type == DEFS->COORD_ANGLE)
      { hdiv_min = vdiv_min = Vars->Coord_Min[i];
        hdiv_max = vdiv_max = Vars->Coord_Max[i];
        restricted = 1; }
      else if (Set_Vars_Coord_Type == DEFS->COORD_RADIUS)
      { double angle;
        angle = RAD2DEG*atan2(Vars->Coord_Max[i], radius);
        hdiv_min = vdiv_min = angle;
        hdiv_max = vdiv_max = angle;
        restricted = 1; }
      else if (Set_Vars_Coord_Type == DEFS->COORD_Y && abs(Vars->Flag_Shape) == DEFS->SHAPE_SPHERE)
      {
        vdiv_min = atan2(ymin,radius)*RAD2DEG;
        vdiv_max = atan2(ymax,radius)*RAD2DEG;
        restricted = 1;
      }
    }
    /* full sphere */
    if ((!restricted && (abs(Vars->Flag_Shape) == DEFS->SHAPE_SPHERE))
    || abs(Vars->Flag_Shape) == DEFS->SHAPE_PREVIOUS)
    {
      mcdis_magnify("");
      mcdis_circle("xy",0,0,0,radius);
      mcdis_circle("xz",0,0,0,radius);
      mcdis_circle("yz",0,0,0,radius);
    }
    /* banana/cylinder/sphere portion */
    else
    if (restricted && ((abs(Vars->Flag_Shape) == DEFS->SHAPE_CYLIND)
                    || (abs(Vars->Flag_Shape) == DEFS->SHAPE_BANANA)
                    || (abs(Vars->Flag_Shape) == DEFS->SHAPE_SPHERE)))
    {
      int NH=24, NV=24;
      int ih, iv;
      double width, height;
      int issphere;
      issphere = (abs(Vars->Flag_Shape) == DEFS->SHAPE_SPHERE);
      width = (hdiv_max-hdiv_min)/NH;
      if (!issphere) NV=1; /* cylinder has vertical axis */
      else height= (vdiv_max-vdiv_min)/NV;
      
      /* check width and height of elements (sphere) to make sure the nb
         of plates remains limited */
      if (width < 10  && NH > 1) { width = 10;  NH=(hdiv_max-hdiv_min)/width; width=(hdiv_max-hdiv_min)/NH; }
      if (height < 10 && NV > 1) { height = 10; NV=(vdiv_max-vdiv_min)/height; height= (vdiv_max-vdiv_min)/NV; }
      
      mcdis_magnify("xyz");
      for(ih = 0; ih < NH; ih++)
        for(iv = 0; iv < NV; iv++)
        {
          double theta0, phi0, theta1, phi1;          /* angles in spherical coordinates */
          double x0,y0,z0,x1,y1,z1,x2,y2,z2,x3,y3,z3; /* vertices at plate edges */
          phi0 = (hdiv_min+ width*ih-90)*DEG2RAD;        /* in xz plane */
          phi1 = (hdiv_min+ width*(ih+1)-90)*DEG2RAD;
          if (issphere)
          {
            theta0= (vdiv_min+height* iv + 90)   *DEG2RAD; /* in vertical plane */
            theta1= (vdiv_min+height*(iv+1) + 90)*DEG2RAD;
            if (y0 < ymin) y0=ymin; 
            if (y0 > ymax) y0=ymax;
            if (y1 < ymin) y1=ymin; 
            if (y1 > ymax) y1=ymax;
            
            y0 = -radius*cos(theta0);            /* z with Z vertical */
            y1 = -radius*cos(theta1);
            if (y0 < ymin) y0=ymin;
            if (y0 > ymax) y0=ymax;
            if (y1 < ymin) y1=ymin;
            if (y1 > ymax) y1=ymax;
          } else {
            y0 = ymin;
            y1 = ymax;
            theta0=theta1=90*DEG2RAD;
          }

          x0 = radius*sin(theta0)*cos(phi0); /* x with Z vertical */
          z0 =-radius*sin(theta0)*sin(phi0); /* y with Z vertical */
          x1 = radius*sin(theta1)*cos(phi0); 
          z1 =-radius*sin(theta1)*sin(phi0);
          x2 = radius*sin(theta1)*cos(phi1); 
          z2 =-radius*sin(theta1)*sin(phi1);
          x3 = radius*sin(theta0)*cos(phi1); 
          z3 =-radius*sin(theta0)*sin(phi1);
          y2 = y1; y3 = y0;

          mcdis_multiline(5,
            x0,y0,z0,
            x1,y1,z1,
            x2,y2,z2,
            x3,y3,z3,
            x0,y0,z0);
        }
      if (Vars->Flag_mantid) {
	/* First define the base pixel type */
	double dt, dy;
	dt = (Vars->Coord_Max[1]-Vars->Coord_Min[1])/Vars->Coord_Bin[1];
	dy = (Vars->Coord_Max[2]-Vars->Coord_Min[2])/Vars->Coord_Bin[2];
	printf("MANTID_BANANA_DET:  %g, %g, %g, %g, %g, %li, %li, %g\n", radius, 
	       Vars->Coord_Min[1],Vars->Coord_Max[1], Vars->Coord_Min[2],Vars->Coord_Max[2], Vars->Coord_Bin[1], Vars->Coord_Bin[2], Vars->Coord_Min[4]); 
      }
    }
    /* disk (circle) */
    else
    if (abs(Vars->Flag_Shape) == DEFS->SHAPE_DISK)
    {
      mcdis_magnify("");
      mcdis_circle("xy",0,0,0,radius);
    }
    /* rectangle (square) */
    else
    if (abs(Vars->Flag_Shape) == DEFS->SHAPE_SQUARE)
    {
      mcdis_magnify("xy");
      mcdis_multiline(5, (double)xmin, (double)ymin, 0.0,
             (double)xmax, (double)ymin, 0.0,
             (double)xmax, (double)ymax, 0.0,
             (double)xmin, (double)ymax, 0.0,
             (double)xmin, (double)ymin, 0.0);
      
      if (Vars->Flag_mantid) {
	/* First define the base pixel type */
	double dx, dy;
	dx = (Vars->Coord_Max[1]-Vars->Coord_Min[1])/Vars->Coord_Bin[1];
	dy = (Vars->Coord_Max[2]-Vars->Coord_Min[2])/Vars->Coord_Bin[2];
	printf("MANTID_RECTANGULAR_DET:  %g, %g, %g, %g, %li, %li, %g\n", 
	       Vars->Coord_Min[1],Vars->Coord_Max[1], Vars->Coord_Min[2],Vars->Coord_Max[2], Vars->Coord_Bin[1], Vars->Coord_Bin[2], Vars->Coord_Min[4]);
      }
    }
    /* full cylinder/banana */
    else
    if (!restricted && ((abs(Vars->Flag_Shape) == DEFS->SHAPE_CYLIND) || (abs(Vars->Flag_Shape) == DEFS->SHAPE_BANANA)))
    {
      mcdis_magnify("xyz");
      mcdis_circle("xz", 0,  h/2.0, 0, radius);
      mcdis_circle("xz", 0, -h/2.0, 0, radius);
      mcdis_line(-radius, -h/2.0, 0, -radius, +h/2.0, 0);
      mcdis_line(+radius, -h/2.0, 0, +radius, +h/2.0, 0);
      mcdis_line(0, -h/2.0, -radius, 0, +h/2.0, -radius);
      mcdis_line(0, -h/2.0, +radius, 0, +h/2.0, +radius);
    }
    else
    /* box */
    if (abs(Vars->Flag_Shape) == DEFS->SHAPE_BOX)
    {
      mcdis_magnify("xyz");
      mcdis_multiline(5, xmin, ymin, zmin,
                   xmax, ymin, zmin,
                   xmax, ymax, zmin,
                   xmin, ymax, zmin,
                   xmin, ymin, zmin);
      mcdis_multiline(5, xmin, ymin, zmax,
                   xmax, ymin, zmax,
                   xmax, ymax, zmax,
                   xmin, ymax, zmax,
                   xmin, ymin, zmax);
      mcdis_line(xmin, ymin, zmin, xmin, ymin, zmax);
      mcdis_line(xmax, ymin, zmin, xmax, ymin, zmax);
      mcdis_line(xmin, ymax, zmin, xmin, ymax, zmax);
      mcdis_line(xmax, ymax, zmin, xmax, ymax, zmax);
    }
  } /* end Monitor_nD_McDisplay */

/* end of monitor_nd-lib.c */


/*******************************************************************************
*
* McStas, neutron ray-tracing package
*         Copyright (C) 1997-2008, All rights reserved
*         Risoe National Laboratory, Roskilde, Denmark
*         Institut Laue Langevin, Grenoble, France
*
* Runtime: share/interoff.h
*
* %Identification
* Written by: Reynald Arnerin
* Date:    Jun 12, 2008
* Release: 
* Version: 
*
* Object File Format intersection header for McStas. Requires the qsort function.
*
* Such files may be obtained with e.g.
*   qhull < points.xyz Qx Qv Tv o > points.off
* where points.xyz has format:
*   3
*   <nb_points>
*   <x> <y> <z>
*   ...
* The resulting file should have its first line being changed from '3' into 'OFF'.
* It can then be displayed with geomview.
* A similar, but somewhat older solution is to use 'powercrust' with e.g.
*   powercrust -i points.xyz
* which will generate a 'pc.off' file to be renamed as suited.
*
*******************************************************************************/

#ifndef INTEROFF_LIB_H
#define INTEROFF_LIB_H "$Revision$"

#ifndef EPSILON
#define EPSILON 1e-13
#endif

#define OFF_INTERSECT_MAX 100

//#include <float.h>

#define N_VERTEX_DISPLAYED    200000

typedef struct intersection {
	MCNUM time;  	  //time of the intersection
	Coords v;	      //intersection point
	Coords normal;  //normal vector of the surface intersected
	short in_out;	  //1 if the ray enters the volume, -1 otherwise
	short edge;	    //1 if the intersection is on the boundary of the polygon, and error is possible
	unsigned long index; // index of the face
} intersection;

typedef struct polygon {
  MCNUM* p;       //vertices of the polygon in adjacent order, this way : x1 | y1 | z1 | x2 | y2 | z2 ...
  int npol;       //number of vertices
  Coords normal;
} polygon;

typedef struct off_struct {
    long vtxSize;
    long polySize;
    long faceSize;
    Coords* vtxArray;
    Coords* normalArray;
    unsigned long* faceArray;
    char *filename;
    int mantidflag;
    long mantidoffset;
    intersection intersects[OFF_INTERSECT_MAX]; // After a call to off_intersect_all contains the list of intersections.
    int nextintersect;                 // 'Next' intersection (first t>0) solution after call to off_intersect_all
    int numintersect;               // Number of intersections after call to off_intersect_all
} off_struct;

/*******************************************************************************
* long off_init(  char *offfile, double xwidth, double yheight, double zdepth, off_struct* data)
* ACTION: read an OFF file, optionally center object and rescale, initialize OFF data structure
* INPUT: 'offfile' OFF file to read
*        'xwidth,yheight,zdepth' if given as non-zero, apply bounding box. 
*           Specifying only one of these will also use the same ratio on all axes
*        'notcenter' center the object to the (0,0,0) position in local frame when set to zero
* RETURN: number of polyhedra and 'data' OFF structure 
*******************************************************************************/
long off_init(  char *offfile, double xwidth, double yheight, double zdepth, 
                int notcenter, off_struct* data);

/*******************************************************************************
* int off_intersect_all(double* t0, double* t3, 
     Coords *n0, Coords *n3,
     double x, double y, double z, 
     double vx, double vy, double vz, 
     off_struct *data )
* ACTION: computes intersection of neutron trajectory with an object. 
* INPUT:  x,y,z and vx,vy,vz are the position and velocity of the neutron
*         data points to the OFF data structure
* RETURN: the number of polyhedra which trajectory intersects
*         t0 and t3 are the smallest incoming and outgoing intersection times
*         n0 and n3 are the corresponding normal vectors to the surface
*         data is the full OFF structure, including a list intersection type
*******************************************************************************/
int off_intersect_all(double* t0, double* t3, 
     Coords *n0, Coords *n3,
     double x, double y, double z, 
     double vx, double vy, double vz, 
     off_struct *data );

/*******************************************************************************
* int off_intersect(double* t0, double* t3, 
     Coords *n0, Coords *n3,
     double x, double y, double z, 
     double vx, double vy, double vz, 
     off_struct data )
* ACTION: computes intersection of neutron trajectory with an object. 
* INPUT:  x,y,z and vx,vy,vz are the position and velocity of the neutron
*         data points to the OFF data structure
* RETURN: the number of polyhedra which trajectory intersects
*         t0 and t3 are the smallest incoming and outgoing intersection times
*         n0 and n3 are the corresponding normal vectors to the surface
*******************************************************************************/
int off_intersect(double* t0, double* t3, 
     Coords *n0, Coords *n3,
     double x, double y, double z, 
     double vx, double vy, double vz, 
     off_struct data );

/*****************************************************************************
* int off_intersectx(double* l0, double* l3, 
     Coords *n0, Coords *n3,
     double x, double y, double z, 
     double kx, double ky, double kz, 
     off_struct data )
* ACTION: computes intersection of an xray trajectory with an object.
* INPUT:  x,y,z and kx,ky,kz, are spatial coordinates and wavevector of the x-ray
*         respectively. data points to the OFF data structure.
* RETURN: the number of polyhedra the trajectory intersects
*         l0 and l3 are the smallest incoming and outgoing intersection lengths
*         n0 and n3 are the corresponding normal vectors to the surface
*******************************************************************************/
int off_x_intersect(double *l0,double *l3,
     Coords *n0, Coords *n3,
     double x,  double y,  double z, 
     double kx, double ky, double kz, 
     off_struct data );

/*******************************************************************************
* void off_display(off_struct data)
* ACTION: display up to N_VERTEX_DISPLAYED points from the object
*******************************************************************************/
void off_display(off_struct);

#endif

/* end of interoff-lib.h */
/*******************************************************************************
*
* McStas, neutron ray-tracing package
*         Copyright (C) 1997-2008, All rights reserved
*         Risoe National Laboratory, Roskilde, Denmark
*         Institut Laue Langevin, Grenoble, France
*
* Runtime: share/interoff-lib.c
*
* %Identification
* Written by: Reynald Arnerin
* Date:    Jun 12, 2008
* Origin: ILL
* Release: $Revision$
* Version: McStas X.Y
*
* Object File Format intersection library for McStas. Requires the qsort function.
*
* Such files may be obtained with e.g.
*   qhull < points.xyz Qx Qv Tv o > points.off
* where points.xyz has format (it supports comments):
*   3
*   <nb_points>
*   <x> <y> <z>
*   ...
* The resulting file should have its first line being changed from '3' into 'OFF'.
* It can then be displayed with geomview.
* A similar, but somewhat older solution is to use 'powercrust' with e.g.
*   powercrust -i points.xyz
* which will generate a 'pc.off' file to be renamed as suited.
*
*******************************************************************************/

#ifndef INTEROFF_LIB_H
#include "interoff-lib.h"
#endif

double off_F(double x, double y,double z,double A,double B,double C,double D) {
  return ( A*x + B*y + C*z + D );
}

char off_sign(double a) {
  if (a<0)       return(-1);
  else if (a==0) return(0);
  else           return(1);
}

// off_normal ******************************************************************
//gives the normal vector of p
void off_normal(Coords* n, polygon p)
{
  //using Newell method
  int i=0,j=0;
  n->x=0;n->y=0;n->z=0;
  for (i = 0, j = p.npol-1; i < p.npol; j = i++)
  {
    MCNUM x1=p.p[3*i],
          y1=p.p[3*i+1],
          z1=p.p[3*i+2];
    MCNUM x2=p.p[3*j],
          y2=p.p[3*j+1],
          z2=p.p[3*j+2];
    // n is the cross product of v1*v2
    n->x += (y1 - y2) * (z1 + z2);
    n->y += (z1 - z2) * (x1 + x2);
    n->z += (x1 - x2) * (y1 + y2);
  }
} /* off_normal */

// off_pnpoly ******************************************************************
//based on http://www.ecse.rpi.edu/Homepages/wrf/Research/Short_Notes/pnpoly.html
//return 0 if the vertex is out
//    1 if it is in
//   -1 if on the boundary
int off_pnpoly(polygon p, Coords v)
{
  int i=0, c = 0;
  MCNUM minx=FLT_MAX,maxx=-FLT_MAX,miny=FLT_MAX,maxy=-FLT_MAX,minz=FLT_MAX,maxz=-FLT_MAX;
  MCNUM rangex=0,rangey=0,rangez=0;

  int pol2dx=0,pol2dy=1;          //2d restriction of the poly
  MCNUM x=v.x,y=v.y;


  //take the most relevant 2D projection (prevent from instability)
  for (i=0; i<p.npol; ++i)
  {
    if (p.p[3*i]<minx)   minx=p.p[3*i];
    if (p.p[3*i]>maxx)   maxx=p.p[3*i];
    if (p.p[3*i+1]<miny) miny=p.p[3*i+1];
    if (p.p[3*i+1]>maxy) maxy=p.p[3*i+1];
    if (p.p[3*i+2]<minz) minz=p.p[3*i+2];
    if (p.p[3*i+2]>maxz) maxz=p.p[3*i+2];
  }
  rangex=maxx-minx;
  rangey=maxy-miny;
  rangez=maxz-minz;

  if (rangex<rangez)
  {
    if (rangex<rangey) {
      pol2dx=2;
      x=v.z;
    } else {
      pol2dy=2;
      y=v.z;
    }
  }
  else if (rangey<rangez) {
    pol2dy=2;
    y=v.z;
  }

  //trace rays and test number of intersection
  int j;
  for (i = 0, j = p.npol-1; i < p.npol; j = i++) {
    if (((((p.p[3*i+pol2dy])<=y) && (y<(p.p[3*j+pol2dy]))) ||
         (((p.p[3*j+pol2dy])<=y) && (y<(p.p[3*i+pol2dy])))) &&
        (x < ( (p.p[3*j+pol2dx] - p.p[3*i+pol2dx]) * (y - p.p[3*i+pol2dy])
             / (p.p[3*j+pol2dy] - p.p[3*i+pol2dy]) + p.p[3*i+pol2dx]) ))
      c = !c;

    if (((fabs(p.p[3*i+pol2dy]-y)<=EPSILON) || ((fabs(p.p[3*j+pol2dy]-y)<=EPSILON))) &&
        fabs(x -((p.p[3*j+pol2dx] - p.p[3*i+pol2dx]) * (y - p.p[3*i+pol2dy])
          / (p.p[3*j+pol2dy] - p.p[3*i+pol2dy]) + p.p[3*i+pol2dx])) < EPSILON)
    {
      //the point lies on the edge
      c=-1;
      break;
    }
  }

  return c;
} /* off_pnpoly */

// off_intersectPoly ***********************************************************
//gives the intersection vertex between ray [a,b) and polygon p and its parametric value on (a b)
//based on http://geometryalgorithms.com/Archive/algorithm_0105/algorithm_0105.htm
int off_intersectPoly(intersection *inter, Coords a, Coords b, polygon p)
{
  //direction vector of [a,b]
  Coords dir = {b.x-a.x, b.y-a.y, b.z-a.z};

  //the normal vector to the polygon
  Coords normale=p.normal;
  //off_normal(&normale, p); done at the init stage

  //direction vector from a to a vertex of the polygon
  Coords w0 = {a.x-p.p[0], a.y-p.p[1], a.z-p.p[2]};

  //scalar product
  MCNUM nw0  =-scalar_prod(normale.x,normale.y,normale.z,w0.x,w0.y,w0.z);
  MCNUM ndir = scalar_prod(normale.x,normale.y,normale.z,dir.x,dir.y,dir.z);
  inter->time = inter->edge = inter->in_out=0;
  inter->v = inter->normal = coords_set(0,0,1);

  if (fabs(ndir) < EPSILON)    // ray is parallel to polygon plane
  {
    if (nw0 == 0)              // ray lies in polygon plane (infinite number of solution)
      return 0;
    else return 0;             // ray disjoint from plane (no solution)
  }

  // get intersect point of ray with polygon plane
  inter->time = nw0 / ndir;            //parametric value the point on line (a,b)

  inter->v = coords_set(a.x + inter->time * dir.x,// intersect point of ray and plane
    a.y + inter->time * dir.y,
    a.z + inter->time * dir.z);

  int res=off_pnpoly(p,inter->v);

  inter->edge=(res==-1);
  if (ndir<0)
    inter->in_out=1;  //the negative dot product means we enter the surface
  else
    inter->in_out=-1;

  inter->normal=p.normal;

  return res;         //true if the intersection point lies inside the poly
} /* off_intersectPoly */


// off_getBlocksIndex **********************************************************
/*reads the indexes at the beginning of the off file as this :
line 1  OFF
line 2  nbVertex nbFaces nbEdges
*/
FILE *off_getBlocksIndex(char* filename, long* vtxSize, long* polySize )
{
  FILE* f = Open_File(filename,"r", NULL); /* from read_table-lib: FILE *Open_File(char *name, char *Mode, char *path) */
  if (!f) return (f);
  
  char line[CHAR_BUF_LENGTH];
  char *ret=0;
  *vtxSize = *polySize = 0;

  /* **************** start to read the file header */
  /* OFF file:
     'OFF' or '3'
   */

  ret=fgets(line,CHAR_BUF_LENGTH , f);// line 1 = "OFF"
  if (ret == NULL)
  {
    fprintf(stderr, "Error: Can not read 1st line in file %s (interoff/off_getBlocksIndex)\n", filename);
    exit(1);
  }
  if (strlen(line)>5)
  {
      fprintf(stderr,"Error: First line in %s is too long (=%lu). Possibly the line is not terminated by '\\n'.\n" 
              "       The first line is required to be exactly 'OFF', '3' or 'ply'.\n",
              filename,(long unsigned)strlen(line));
      fclose(f);
      return(NULL);
  }

  if (strncmp(line,"OFF",3) && strncmp(line,"3",1) && strncmp(line,"ply",1))
  {
    fprintf(stderr, "Error: %s is probably not an OFF, NOFF or PLY file (interoff/off_getBlocksIndex).\n"
                    "       Requires first line to be 'OFF', '3' or 'ply'.\n",filename);
    fclose(f);
    return(NULL);
  }

  if (!strncmp(line,"OFF",3) || !strncmp(line,"3",1)) {
    do  /* OFF file: skip # comments which may be there */
    {
      ret=fgets(line,CHAR_BUF_LENGTH , f);
      if (ret == NULL)
      {
        fprintf(stderr, "Error: Can not read line in file %s (interoff/off_getBlocksIndex)\n", filename);
        exit(1);
      }
    } while (line[0]=='#');
    //line = nblines of vertex,faces and edges arrays
    sscanf(line,"%lu %lu",vtxSize,polySize);
  } else {
    do  /* PLY file: read all lines until find 'end_header'
           and locate 'element faces' and 'element vertex' */
    {
      ret=fgets(line,CHAR_BUF_LENGTH , f);
      if (ret == NULL)
      {
        fprintf(stderr, "Error: Can not read line in file %s (interoff/off_getBlocksIndex)\n", filename);
        exit(1);
      }
      if (!strncmp(line,"element face",12))
        sscanf(line,"element face %lu",polySize);
      else if (!strncmp(line,"element vertex",14))
        sscanf(line,"element vertex %lu",vtxSize);
      else if (!strncmp(line,"format binary",13))
        exit(fprintf(stderr,
          "Error: Can not read binary PLY file %s, only 'format ascii' (interoff/off_getBlocksIndex)\n%s\n",
          filename, line));
    } while (strncmp(line,"end_header",10));
  }
  
  /* The FILE is left opened ready to read 'vtxSize' vertices (vtxSize *3 numbers)
     and then polySize polygons (rows) */

  return(f);
} /* off_getBlocksIndex */

// off_init_planes *************************************************************
//gives the equations of 2 perpandicular planes of [ab]
void off_init_planes(Coords a, Coords b,
  MCNUM* A1, MCNUM* C1, MCNUM* D1, MCNUM *A2, MCNUM* B2, MCNUM* C2, MCNUM* D2)
{
  //direction vector of [a b]
  Coords dir={b.x-a.x, b.y-a.y, b.z-a.z};

  //the plane parallel to the 'y' is computed with the normal vector of the projection of [ab] on plane 'xz'
  *A1= dir.z;
  *C1=-dir.x;
  if(*A1!=0 || *C1!=0)
    *D1=-(a.x)*(*A1)-(a.z)*(*C1);
  else
  {
    //the plane does not support the vector, take the one parallel to 'z''
    *A1=1;
    //B1=dir.x=0
    *D1=-(a.x);
  }
  //the plane parallel to the 'x' is computed with the normal vector of the projection of [ab] on plane 'yz'
  *B2= dir.z;
  *C2=-dir.y;
  *A2= 0;
  if (*B2==0 && *C2==0)
  {
    //the plane does not support the vector, take the one parallel to 'z'
    *B2=1;
    //B1=dir.x=0
    *D2=-(a.y);
  }
  else {
    if (dir.z==0)
    {
      //the planes are the same, take the one parallel to 'z'
      *A2= dir.y;
      *B2=-dir.x;
      *D2=-(a.x)*(*A2)-(a.y)*(*B2);
    }
    else
      *D2=-(a.y)**B2-(a.z)**C2;
  }
} /* off_init_planes */

// off_clip_3D_mod *************************************************************
int off_clip_3D_mod(intersection* t, Coords a, Coords b,
  Coords* vtxArray, unsigned long vtxSize, unsigned long* faceArray,
  unsigned long faceSize, Coords* normalArray)
{
  MCNUM A1=0, C1=0, D1=0, A2=0, B2=0, C2=0, D2=0;      //perpendicular plane equations to [a,b]
  off_init_planes(a, b, &A1, &C1, &D1, &A2, &B2, &C2, &D2);

  int t_size=0;
  //unsigned long vtxSize=vtxTable.rows, faceSize=faceTable.columns;  //Size of the corresponding tables
  char sg[vtxSize];  //array telling if vertex is left or right of the plane
  MCNUM popol[3*CHAR_BUF_LENGTH];
  unsigned long i=0,indPoly=0;
  for (i=0; i < vtxSize; ++i)
  {
    sg[i]=off_sign(off_F(vtxArray[i].x,vtxArray[i].y,vtxArray[i].z,A1,0,C1,D1));
  }

  //exploring the polygons :
  i=indPoly=0;
  while (i<faceSize)
  {
    polygon pol;
    pol.npol  = faceArray[i];                //nb vertex of polygon
    pol.p     = popol;
    pol.normal= coords_set(0,0,1);
    unsigned long indVertP1=faceArray[++i];  //polygon's first vertex index in vtxTable
    int j=1;
    while (j<pol.npol)
    {
      //polygon's j-th vertex index in vtxTable
      if (sg[indVertP1]!=sg[faceArray[i+j]]) //if the plane intersect the polygon
        break;

      ++j;
    }

    if (j<pol.npol)          //ok, let's test with the second plane
    {
      char sg1=off_sign(off_F(vtxArray[indVertP1].x,vtxArray[indVertP1].y,vtxArray[indVertP1].z,A2,B2,C2,D2));//tells if vertex is left or right of the plane

      j=1;
      while (j<pol.npol)
      {
        //unsigned long indVertPi=faceArray[i+j];  //polyg's j-th vertex index in vtxTable
        Coords vertPi=vtxArray[faceArray[i+j]];
        if (sg1!=off_sign(off_F(vertPi.x,vertPi.y,vertPi.z,A2,B2,C2,D2)))//if the plane intersect the polygon
          break;
        ++j;
      }
      if (j<pol.npol)
      {
        if (t_size>CHAR_BUF_LENGTH)
        {
          fprintf(stderr, "Warning: number of intersection exceeded (%d) (interoff-lib/off_clip_3D_mod)\n", CHAR_BUF_LENGTH);
            return (t_size);
        }
        //both planes intersect the polygon, let's find the intersection point
        //our polygon :
        int k;
        for (k=0; k<pol.npol; ++k)
        {
          Coords vertPk=vtxArray[faceArray[i+k]];
          pol.p[3*k]  =vertPk.x;
          pol.p[3*k+1]=vertPk.y;
          pol.p[3*k+2]=vertPk.z;
        }
        pol.normal=normalArray[indPoly];
        intersection x;
        if (off_intersectPoly(&x, a, b, pol))
        {
          x.index = indPoly;
          t[t_size++]=x;
        }
      } /* if (j<pol.npol) */
    } /* if (j<pol.npol) */
    i += pol.npol;
    indPoly++;
  } /* while i<faceSize */
  return t_size;
} /* off_clip_3D_mod */


// off_compare *****************************************************************
int off_compare (void const *a, void const *b)
{
   intersection const *pa = a;
   intersection const *pb = b;

   return off_sign(pa->time - pb->time);
} /* off_compare */

// off_cleanDouble *************************************************************
//given an array of intersections throw those which appear several times
//returns 1 if there is a possibility of error
int off_cleanDouble(intersection* t, int* t_size)
{
  int i=1;
  intersection prev=t[0];
  while (i<*t_size)
  {
    int j=i;
    //for each intersection with the same time
    while (j<*t_size && fabs(prev.time-t[j].time)<EPSILON)
    {
      //if the intersection is the exact same erase it
      if (prev.in_out==t[j].in_out)
      {
        int k;
        for (k=j+1; k<*t_size; ++k)
        {
          t[k-1]=t[k];
        }
        *t_size-=1;
      }
      else
        ++j;
    }
    prev=t[i];
    ++i;

  }
  return 1;
} /* off_cleanDouble */

// off_cleanInOut **************************************************************
//given an array of intesections throw those which enter and exit in the same time
//Meaning the ray passes very close to the volume
//returns 1 if there is a possibility of error
int off_cleanInOut(intersection* t, int* t_size)
{
  int i=1;
  intersection prev=t[0];
  while (i<*t_size)
  {
    //if two intersection have the same time but one enters and the other exits erase both
    //(such intersections must be adjacent in the array : run off_cleanDouble before)
    if (fabs(prev.time-t[i].time)<EPSILON && prev.in_out!=t[i].in_out)
    {
      int j=0;
      for (j=i+1; j<*t_size; ++j)
      {
        t[j-2]=t[j];
      }
      *t_size-=2;
      prev=t[i-1];
    }
    else
    {
      prev=t[i];
      ++i;
    }
  }
  return (*t_size);
} /* off_cleanInOut */

/* PUBLIC functions ******************************************************** */

/*******************************************************************************
* long off_init(  char *offfile, double xwidth, double yheight, double zdepth, off_struct* data)
* ACTION: read an OFF file, optionally center object and rescale, initialize OFF data structure
* INPUT: 'offfile' OFF file to read
*        'xwidth,yheight,zdepth' if given as non-zero, apply bounding box.
*           Specifying only one of these will also use the same ratio on all axes
*        'notcenter' center the object to the (0,0,0) position in local frame when set to zero
* RETURN: number of polyhedra and 'data' OFF structure
*******************************************************************************/
long off_init(  char *offfile, double xwidth, double yheight, double zdepth,
                int notcenter, off_struct* data)
{
  // data to be initialized
  long    vtxSize =0, polySize=0, i=0, ret=0, faceSize=0;
  Coords* vtxArray        =NULL;
  Coords* normalArray     =NULL;
  unsigned long* faceArray=NULL;
  FILE*   f               =NULL; /* the FILE with vertices and polygons */
  double minx=FLT_MAX,maxx=-FLT_MAX,miny=FLT_MAX,maxy=-FLT_MAX,minz=FLT_MAX,maxz=-FLT_MAX;

  // get the indexes
  if (!data) return(0);
  
  MPI_MASTER(
  printf("Loading geometry file (OFF/PLY): %s\n", offfile);
  );
  
  f=off_getBlocksIndex(offfile,&vtxSize,&polySize);
  if (!f) return(0);
  
  // read vertex table = [x y z | x y z | ...] =================================
  // now we read the vertices as 'vtxSize*3' numbers and store it in vtxArray 
  MPI_MASTER(
  printf("  Number of vertices: %ld\n", vtxSize);
  );
  vtxArray   = malloc(vtxSize*sizeof(Coords));
  if (!vtxArray) return(0);
  i=0;
  while (i<vtxSize && ~feof(f))
  {
    double x,y,z;
    ret=fscanf(f, "%lg%lg%lg", &x,&y,&z);
    if (!ret) { 
      // invalid line: we skip it (probably a comment)
      char line[CHAR_BUF_LENGTH];
      char *s=fgets(line, CHAR_BUF_LENGTH, f);
      continue; 
    }
    if (ret != 3) {
      fprintf(stderr, "Error: can not read [xyz] coordinates for vertex %li in file %s (interoff/off_init). Read %li values.\n", 
        i, offfile, ret);
      exit(2);
    }
    vtxArray[i].x=x;
    vtxArray[i].y=y;
    vtxArray[i].z=z;

    //bounding box
    if (vtxArray[i].x<minx) minx=vtxArray[i].x;
    if (vtxArray[i].x>maxx) maxx=vtxArray[i].x;
    if (vtxArray[i].y<miny) miny=vtxArray[i].y;
    if (vtxArray[i].y>maxy) maxy=vtxArray[i].y;
    if (vtxArray[i].z<minz) minz=vtxArray[i].z;
    if (vtxArray[i].z>maxz) maxz=vtxArray[i].z;
    i++; // inquire next vertex
  }

  // resizing and repositioning params
  double centerx=0, centery=0, centerz=0;
  if (!notcenter) {
    centerx=(minx+maxx)*0.5;
    centery=(miny+maxy)*0.5;
    centerz=(minz+maxz)*0.5;
  }

  double rangex=-minx+maxx,
         rangey=-miny+maxy,
         rangez=-minz+maxz;

  double ratiox=1,ratioy=1,ratioz=1;

  if (xwidth && rangex)
  {
    ratiox=xwidth/rangex;
    ratioy=ratiox;
    ratioz=ratiox;
  }

  if (yheight && rangey)
  {
    ratioy=yheight/rangey;
    if(!xwidth)  ratiox=ratioy;
    ratioz=ratioy;
  }

  if (zdepth && rangez)
  {
    ratioz=zdepth/rangez;
    if(!xwidth)  ratiox=ratioz;
    if(!yheight) ratioy=ratioz;
  }

  rangex *= ratiox;
  rangey *= ratioy;
  rangez *= ratioz;

  //center and resize the object
  for (i=0; i<vtxSize; ++i)
  {
    vtxArray[i].x=(vtxArray[i].x-centerx)*ratiox+(!notcenter ? 0 : centerx);
    vtxArray[i].y=(vtxArray[i].y-centery)*ratioy+(!notcenter ? 0 : centery);
    vtxArray[i].z=(vtxArray[i].z-centerz)*ratioz+(!notcenter ? 0 : centerz);
  }
  
  // read face table = [nbvertex v1 v2 vn | nbvertex v1 v2 vn ...] =============
  MPI_MASTER(
  printf("  Number of polygons: %ld\n", polySize);
  );
  normalArray= malloc(polySize*sizeof(Coords));
  faceArray  = malloc(polySize*10*sizeof(unsigned long)); // we assume polygons have less than 9 vertices
  if (!normalArray || !faceArray) return(0);
  
  // fill faces
  faceSize=0;
  i=0;
  while (i<polySize && ~feof(f)) {
    int  nbVertex=0, j=0;
    // read the length of this polygon
    ret=fscanf(f, "%d", &nbVertex);
    if (!ret) { 
      // invalid line: we skip it (probably a comment)
      char line[CHAR_BUF_LENGTH];
      char *s=fgets(line, CHAR_BUF_LENGTH, f);
      continue; 
    }
    if (ret != 1) {
      fprintf(stderr, "Error: can not read polygon %li length in file %s (interoff/off_init)\n", 
        i, offfile);
      exit(3);
    }
    if (faceSize > polySize*10) {
      fprintf(stderr, "Error: %li exceeded allocated polygon array[%li] in file %s (interoff/off_init)\n", 
        faceSize, polySize*10, offfile);
    }
    faceArray[faceSize++] = nbVertex; // length of the polygon/face
    // then read the vertex ID's
    for (j=0; j<nbVertex; j++) {
      double vtx=0;
      ret=fscanf(f, "%lg", &vtx);
      faceArray[faceSize++] = vtx;   // add vertices index after length of polygon
    }
    i++;
  }

  // precomputes normals
  long indNormal=0;//index in polyArray
  i=0;    //index in faceArray
  while (i<faceSize)
  {
    int    nbVertex=faceArray[i];//nb of vertices of this polygon
    double vertices[3*nbVertex];
    int j;

    for (j=0; j<nbVertex; ++j)
    {
      unsigned long indVertPj=faceArray[i+j+1];
      vertices[3*j]  =vtxArray[indVertPj].x;
      vertices[3*j+1]=vtxArray[indVertPj].y;
      vertices[3*j+2]=vtxArray[indVertPj].z;
    }

    polygon p;
    p.p   =vertices;
    p.npol=nbVertex;
    off_normal(&(p.normal),p);

    normalArray[indNormal]=p.normal;

    i += nbVertex+1;
    indNormal++;

  }
  
  MPI_MASTER(
  if (ratiox!=ratioy || ratiox!=ratioz || ratioy!=ratioz)
    printf("Warning: Aspect ratio of the geometry %s was modified.\n"
           "         If you want to keep the original proportions, specifiy only one of the dimensions.\n",
           offfile);
  if ( xwidth==0 && yheight==0 && zdepth==0 ) {
    printf("Warning: Neither xwidth, yheight or zdepth are defined.\n"
	   "           The file-defined (non-scaled) geometry the OFF geometry %s will be applied!\n", 
           offfile);
  }
  printf("  Bounding box dimensions for geometry %s:\n", offfile);
  printf("    Length=%f (%.3f%%)\n", rangex, ratiox*100);
  printf("    Width= %f (%.3f%%)\n", rangey, ratioy*100);
  printf("    Depth= %f (%.3f%%)\n", rangez, ratioz*100);
  );

  data->vtxArray   = vtxArray;
  data->normalArray= normalArray;
  data->faceArray  = faceArray;
  data->vtxSize    = vtxSize;
  data->polySize   = polySize;
  data->faceSize   = faceSize;
  data->filename   = offfile;
  return(polySize);
} /* off_init */

/*******************************************************************************
* int off_intersect_all(double* t0, double* t3,
     Coords *n0, Coords *n3,
     double x, double y, double z,
     double vx, double vy, double vz,
     off_struct *data )
* ACTION: computes intersection of neutron trajectory with an object.
* INPUT:  x,y,z and vx,vy,vz are the position and velocity of the neutron
*         data points to the OFF data structure
* RETURN: the number of polyhedra which trajectory intersects
*         t0 and t3 are the smallest incoming and outgoing intersection times
*         n0 and n3 are the corresponding normal vectors to the surface
*         data is the full OFF structure, including a list intersection type
*******************************************************************************/
int off_intersect_all(double* t0, double* t3,
     Coords *n0, Coords *n3,
     double x,  double y,  double z,
     double vx, double vy, double vz,
     off_struct *data )
{
    Coords A={x, y, z};
    Coords B={x+vx, y+vy, z+vz};
    int t_size=off_clip_3D_mod(data->intersects, A, B,
      data->vtxArray, data->vtxSize, data->faceArray, data->faceSize, data->normalArray );
    qsort(data->intersects, t_size, sizeof(intersection),  off_compare);
    off_cleanDouble(data->intersects, &t_size);
    off_cleanInOut(data->intersects,  &t_size);

    /*find intersections "closest" to 0 (favouring positive ones)*/
    if(t_size>0){
      int i=0;
      if(t_size>1) {
        for (i=1; i < t_size-1; i++){
          if (data->intersects[i-1].time > 0 && data->intersects[i].time > 0)
            break;
        }
	
	data->nextintersect=i-1;
	data->numintersect=t_size;

        if (t0) *t0 = data->intersects[i-1].time;
        if (n0) *n0 = data->intersects[i-1].normal;
        if (t3) *t3 = data->intersects[i].time;
        if (n3) *n3 = data->intersects[i].normal;
      } else {
        if (t0) *t0 = data->intersects[0].time; 	 
	      if (n0) *n0 = data->intersects[0].normal;
      }
      /* should also return t[0].index and t[i].index as polygon ID */
      return t_size;
    }
    return 0;
} /* off_intersect */

/*******************************************************************************
* int off_intersect(double* t0, double* t3,
     Coords *n0, Coords *n3,
     double x, double y, double z,
     double vx, double vy, double vz,
     off_struct data )
* ACTION: computes intersection of neutron trajectory with an object.
* INPUT:  x,y,z and vx,vy,vz are the position and velocity of the neutron
*         data points to the OFF data structure
* RETURN: the number of polyhedra which trajectory intersects
*         t0 and t3 are the smallest incoming and outgoing intersection times
*         n0 and n3 are the corresponding normal vectors to the surface
*******************************************************************************/
int off_intersect(double* t0, double* t3,
     Coords *n0, Coords *n3,
     double x,  double y,  double z,
     double vx, double vy, double vz,
     off_struct data )
{
  return off_intersect_all(t0, t3, n0, n3, x, y, z, vx, vy, vz, &data );
} /* off_intersect */

/*****************************************************************************
* int off_x_intersect(double* l0, double* l3,
     Coords *n0, Coords *n3,
     double x, double y, double z,
     double kx, double ky, double kz,
     off_struct data )
* ACTION: computes intersection of an xray trajectory with an object.
* INPUT:  x,y,z and kx,ky,kz, are spatial coordinates and wavevector of the x-ray
*         respectively. data points to the OFF data structure.
* RETURN: the number of polyhedra the trajectory intersects
*         l0 and l3 are the smallest incoming and outgoing intersection lengths
*         n0 and n3 are the corresponding normal vectors to the surface
*******************************************************************************/
int off_x_intersect(double *l0,double *l3,
     Coords *n0, Coords *n3,
     double x,  double y,  double z,
     double kx, double ky, double kz,
     off_struct data )
{
  /*This function simply reformats and calls off_intersect (as for neutrons)
   *by normalizing the wavevector - this will yield the intersection lengths
   *in m*/
  double jx,jy,jz,invk;
  int n;
  invk=1/sqrt(scalar_prod(kx,ky,kz,kx,ky,kz));
  jx=kx*invk;jy=ky*invk;jz=kz*invk;
  n=off_intersect(l0,l3,n0,n3,x,y,z,jx,jy,jz,data);
  return n;
}


/*******************************************************************************
* void off_display(off_struct data)
* ACTION: display up to N_VERTEX_DISPLAYED polygons from the object
*******************************************************************************/
void off_display(off_struct data)
{
  unsigned int i;
  double ratio=(double)(N_VERTEX_DISPLAYED)/(double)data.faceSize;
  unsigned int pixel=0;
  for (i=0; i<data.faceSize-1; i++) {
    int j;
    int nbVertex = data.faceArray[i];
    double x0,y0,z0;
    x0 = data.vtxArray[data.faceArray[i+1]].x;
    y0 = data.vtxArray[data.faceArray[i+1]].y;
    z0 = data.vtxArray[data.faceArray[i+1]].z;
    double x1=x0,y1=y0,z1=z0;
    double cmx=0,cmy=0,cmz=0;
    
    int drawthis = rand01() < ratio;
    // First pass, calculate center of mass location...
    for (j=1; j<=nbVertex; j++) {
      cmx = cmx+data.vtxArray[data.faceArray[i+j]].x;
      cmy = cmy+data.vtxArray[data.faceArray[i+j]].y;
      cmz = cmz+data.vtxArray[data.faceArray[i+j]].z;
    }
    cmx /= nbVertex;
    cmy /= nbVertex;
    cmz /= nbVertex;
    
    char pixelinfo[1024];    
    sprintf(pixelinfo, "%li,%li,%li,%i,%g,%g,%g,%g,%g,%g", data.mantidoffset+pixel, data.mantidoffset, data.mantidoffset+data.polySize-1, nbVertex, cmx, cmy, cmz, x1-cmx, y1-cmy, z1-cmz);
    for (j=2; j<=nbVertex; j++) {
      double x2,y2,z2;
      x2 = data.vtxArray[data.faceArray[i+j]].x;
      y2 = data.vtxArray[data.faceArray[i+j]].y;
      z2 = data.vtxArray[data.faceArray[i+j]].z;
      sprintf(pixelinfo, "%s,%g,%g,%g", pixelinfo, x2-cmx, y2-cmy, z2-cmz); 
      if (ratio > 1 || drawthis) {
	mcdis_line(x1,y1,z1,x2,y2,z2);
      }
      x1 = x2; y1 = y2; z1 = z2;
    }
    if (ratio > 1 || drawthis) {
	mcdis_line(x1,y1,z1,x0,y0,z0);
      }
    if (data.mantidflag) {
      printf("MANTID_PIXEL: %s\n", pixelinfo);
      pixel++;
    }
    i += nbVertex;
  }
} /* off_display */

/* end of interoff-lib.c */


/* Shared user declarations for all components types 'Isotropic_Sqw'. */

#ifndef ISOTROPIC_SQW
#define ISOTROPIC_SQW $Revision$

/* {j d F2 DW Dd inv2d q F} + { Sq if j == -1}*/

#ifndef Crystallographica
#define Crystallographica { 4,5,7,0,0,0,0, 0,0 }
#define Fullprof          { 4,0,8,0,0,5,0, 0,0 }
#define Undefined         { 0,0,0,0,0,0,0, 0,0 }
#define Lazy              {17,6,0,0,0,0,0,13,0 }
#endif
/* special case for [q,Sq] table */
#define qSq               {-1,0,0,0,0,0,1, 0,0 }




/* For the density of states S(w) */
struct Sqw_W_struct
{
  double omega;        /* omega value for the data block */
  double cumul_proba;  /* cumulated intensity (between 0 and 1) */
};

/* For the S(q|w) probabilities */
struct Sqw_Q_struct
{
   double Q;           /* omega value for the data block */
   double cumul_proba; /* normalized cumulated probability */
};

struct Sqw_Data_struct /* contains normalized Sqw data for probabilities, coh and inc */
{
  struct Sqw_W_struct  *SW;     /* P(w)  ~ density of states */
  struct Sqw_Q_struct **SQW;    /* P(Q|w)= probability of each Q with w */

  long  *SW_lookup;
  long **QW_lookup;
  t_Table Sqw; /* S(q,w) rebin from file in range -w_max:w_max and 0:q_max, with exp(-hw/kT) weight */
  t_Table iqSq;         /* sigma(Ei) = sigma/2/Ki^2 * \int q S(q,w) dq dw up to 2*Ki_max */
  long   q_bins;
  long   w_bins;        /* length of q and w vectors/axes from file */
  double q_max, q_step; /* min=0      */
  double w_max, w_step; /* min=-w_max */
  long   lookup_length;
  char   filename[80];
  double intensity;
  double Ei_max;        /* max neutron incoming energy for Sigma=iqSq table */
  long   iqSq_length;
  char   type;
  double q_min_file;
};

struct Sqw_sample_struct { /* global parameters gathered as a structure */
  char   compname[256];

  struct Sqw_Data_struct Data_inc;
  struct Sqw_Data_struct Data_coh;

  double s_abs, s_coh, s_inc; /* material constants */
  double my_s;
  double my_a_v;
  double mat_rho;
  double mat_weight;
  double mat_density;
  double Temperature; /* temperature from the data file */
  int    shape;       /* 0:cylinder, 1:box, 2:sphere 3:any shape*/

  double sqw_threshold;       /* options to tune S(q,w) */
  double sqw_classical;
  double sqw_norm;

  double barns;               /* for powders */
  double Dd, DWfactor;

  double T2E;                 /* constants */
  char   Q_correction[256];
  double sqSE2K;

  int    maxloop;             /* flags to monitor caught warnings */
  int    minevents;
  long   neutron_removed;
  long   neutron_enter;
  long   neutron_pmult;
  long   neutron_exit;
  char   verbose_output;
  int    column_order[9];     /* column signification */
  long   lookup_length;

  double dq, dw; /* q/w transfer */
  char   type;   /* interaction type: c(coherent),             i(incoherent),
                                      V(isotropic incoherent), t(transmitted) */
  /* store information from the last event */
  double ki_x,ki_y,ki_z,kf_x,kf_y,kf_z;
  double ti, tf;
  double vi, vf;
  double ki, kf;
  double theta;

  double mean_scatt;      /* stat to show at the end */
  double mean_abs;
  double psum_scatt;
  double single_coh;
  double single_inc;
  double multi;

  double rw, rq;
};

#include <stdio.h>
#include <math.h>

/* sets a Data S(q,w) to 'NULL' */
void Sqw_Data_init(struct Sqw_Data_struct *Sqw_Data)
{
  Sqw_Data->q_bins       =0;
  Sqw_Data->w_bins       =0;
  Sqw_Data->q_max        =0;
  Sqw_Data->q_step       =1;
  Sqw_Data->w_max        =0;
  Sqw_Data->w_step       =1;
  Sqw_Data->Ei_max       = 0;
  Sqw_Data->lookup_length=100; /* length of lookup tables */
  Sqw_Data->intensity    =0;
  strcpy(Sqw_Data->filename, "");
  Sqw_Data->SW           =NULL;
  Sqw_Data->SQW          =NULL;
  Sqw_Data->SW_lookup    =NULL;
  Sqw_Data->QW_lookup    =NULL;
  Sqw_Data->iqSq_length  =100;
  Sqw_Data->type         = ' ';
  Sqw_Data->q_min_file   = 0;
}

off_struct offdata;

/* gaussian distribution to appply around Bragg peaks in a powder */
double Sqw_powder_gauss(double x, double mean, double rms) {
  return (exp(-(x-mean)*(x-mean)/(2*rms*rms))/(sqrt(2*PI)*rms));
}

/* Sqw_quantum_correction
*
* Return the 'quantum correction factor Q so that:
*
*   S(q, w) = Q(w) S*(q,w)
*   S(q,-w) = exp(-hw/kT) S(q,w)
*   S(q, w) = exp( hw/kT) S(q,-w)
*
* with S*=classical limit and Q(w) defined below. For omega > 0, S(q,w) > S(q,-w)
*
* input:
*   w: energy      [meV]
*   T: temperature [K]
*   type: 'Schofield' or 'Boltzmann'        Q = exp(hw/kT/2)
*         'harmonic'  or 'Bader'            Q = hw/kT./(1-exp(-hw/kT))
*         'standard'  or 'Frommhold'        Q = 2./(1+exp(-hw/kT)) [recommended]
*
* References:
*  B. Hehr, http://www.lib.ncsu.edu/resolver/1840.16/7422 PhD manuscript (2010).
*  S. A. Egorov, K. F. Everitt and J. L. Skinner. J. Phys. Chem., 103, 9494 (1999).
*  P. Schofield. Phys. Rev. Lett., 4, 239 (1960).
*  J. S. Bader and B. J. Berne. J. Chem. Phys., 100, 8359 (1994).
*  T. D. Hone and G. A. Voth. J. Chem. Phys., 121, 6412 (2004).
*  L. Frommhold. Collision-induced absorption in gases, 1 st ed., Cambridge
*    Monographs on Atomic, Molecular, and Chemical Physics, Vol. 2,
*    Cambridge Univ. Press: London (1993).

 */
double Sqw_quantum_correction(double hw, double T, char *type) {
  double Q   = 1;
  double kT  = T/11.605;  /* [K] -> [meV = 1000*KB/e] */
  if (!hw || !T) return 1;
  if (type == NULL || !strcmp(type, "standard")
                   || !strcmp(type, "Frommhold") || !strcmp(type, "default"))
    Q = 2/(1+exp(-hw/kT));
  if (!strcmp(type, "Schofield") || !strcmp(type, "Boltzmann"))
    Q = exp(hw/kT/2);
  if (!strcmp(type, "harmonic") || !strcmp(type, "Bader"))
    Q = hw/kT/(1-exp(-hw/kT));

  return Q;
}

/*****************************************************************************
* Sqw_read_PowderN: Read PowderN data files
*   Returns t_Table array or NULL in case of error
* Used in : Sqw_readfile (1)
*****************************************************************************/
t_Table *Sqw_read_PowderN(struct Sqw_sample_struct *Sqw, t_Table sqwTable)
{
  struct line_data
  {
    double F2;                  /* Value of structure factor */
    double q;                   /* Q vector */
    int j;                      /* Multiplicity */
    double DWfactor;            /* Debye-Waller factor */
    double w;                   /* Intrinsic line width */
  };
  struct line_data *list = NULL;
  double q_count=0, j_count=0, F2_count=0;
  int    mult_count  =0;
  double q_step      =FLT_MAX;
  long   size        =0;
  int    i, index;
  double q_min=0, q_max=0;
  char   flag=0;
  int    list_count=0;
  double q_step_cur;
  char   flag_qSq = 0;

  t_Table *retTable;

  flag_qSq = (Sqw->column_order[8]>0 && Sqw->column_order[6]>0);

  MPI_MASTER(
  if (Sqw->column_order[0] == 4 && Sqw->barns !=0)
    printf("Isotropic_sqw: %s: Powder file probably of type Crystallographica/Fullprof (lau)\n"
           "WARNING:       but F2 unit is set to powder_barns=1 (barns). Intensity might be 100 times too high.\n",
           Sqw->compname);
  if (Sqw->column_order[0] == 17 && Sqw->barns == 0)
    printf("Isotropic_sqw: %s: Powder file probably of type Lazy Pulver (laz)\n"
           "WARNING:       but F2 unit is set to powder_barns=0 (fm^2). Intensity might be 100 times too low.\n",
           Sqw->compname);
  );
  size = sqwTable.rows;
  MPI_MASTER(
  if (Sqw->verbose_output > 0) {
    printf("Isotropic_sqw: Converting %ld %s from %s into S(q,w) data\n",
        size, flag_qSq ? "S(q)" : "powder lines", sqwTable.filename);
  }
  );
  /* allocate line_data array */
  list = (struct line_data*)malloc(size*sizeof(struct line_data));

  for (i=0; i<size; i++)
    {
      double j=0, d=0, w=0, DWfactor=0, F2=0, Sq=-1, q=0;
      int index;

      if (Sqw->Dd >= 0)      w         = Sqw->Dd;
      if (Sqw->DWfactor > 0) DWfactor  = Sqw->DWfactor;

      /* get data from table using columns {j d F2 DW Dd inv2d q} + { Sq }*/
      /* column indexes start at 1, thus need to substract 1 */
      if (Sqw->column_order[0]>0)
        j = Table_Index(sqwTable, i, Sqw->column_order[0]-1);
      if (Sqw->column_order[1]>0)
        d = Table_Index(sqwTable, i, Sqw->column_order[1]-1);
      if (Sqw->column_order[2]>0)
        F2 = Table_Index(sqwTable, i, Sqw->column_order[2]-1);
      if (Sqw->column_order[3]>0)
        DWfactor = Table_Index(sqwTable, i, Sqw->column_order[3]-1);
      if (Sqw->column_order[4]>0)
        w = Table_Index(sqwTable, i, Sqw->column_order[4]-1);
      if (Sqw->column_order[5]>0)  {
        d = Table_Index(sqwTable, i, Sqw->column_order[5]-1); if (d) d = 1/d/2; }
      if (Sqw->column_order[6]>0)
        q = Table_Index(sqwTable, i, Sqw->column_order[6]-1);
      if (Sqw->column_order[7]>0 && !F2)
        {F2= Table_Index(sqwTable, i, Sqw->column_order[7]-1); F2 *= F2;}

      if (Sqw->column_order[8]>0)
        Sq= Table_Index(sqwTable, i, Sqw->column_order[8]-1);

      if (q > 0 && Sq >= 0) F2 = Sq;
      if (d > 0 && q <= 0)  q = 2*PI/d;

      /* assign and check values */
      j = (j > 0 ? j : 0);
      if (flag_qSq) j=1;
      DWfactor = (DWfactor > 0 ? DWfactor : 1);
      w = (w>0 ? w : 0);
      F2 = (F2 >= 0 ? F2 : 0);
      d = (q > 0 ? 2*PI/d : 0);
      if (j == 0 || d == 0 || q == 0) {
        MPI_MASTER(
        printf("Isotropic_sqw: %s: Warning: line %i has invalid definition\n"
               "         (mult=0 or q=0 or d=0)\n", Sqw->compname, i);
        );
        continue;
      }
      list[list_count].j = j;
      list[list_count].q = q;
      list[list_count].DWfactor = DWfactor;
      list[list_count].w = w;
      list[list_count].F2= F2; /* or S(q) if flag_qSq */

      if (q_max < d) q_max = q;
      if (q_min > d) q_min = q;
      if (list_count > 1) {
        q_step_cur = fabs(list[list_count].q - list[list_count-1].q);
        if (q_step_cur > 1e-5 && (!q_step || q_step_cur < q_step))
         q_step = q_step_cur;
      }

      /* adjust multiplicity if j-column + multiple d-spacing lines */
      /* if  d = previous d, increase line duplication index */
      if (!q_count)     q_count = q;
      if (!j_count)     j_count = j;
      if (!F2_count)    F2_count= F2;
      if (fabs(q_count-q) < 0.0001*fabs(q)
       && fabs(F2_count-F2) < 0.0001*fabs(F2) && j_count == j) {
       mult_count++; flag=0; }
      else flag=1;
      if (i == size-1) flag=1;
      /* else if d != previous d : just passed equivalent lines */
      if (flag) {
        if (i == size-1) list_count++;
      /*   if duplication index == previous multiplicity */
      /*      set back multiplicity of previous lines to 1 */
        if (Sqw->verbose_output > 2 && (mult_count == list[list_count-1].j
        || (mult_count == list[list_count].j && i == size-1))) {
          MPI_MASTER(
          printf("Isotropic_Sqw: %s: Setting multiplicity to 1 for lines [%i:%i]\n"
                  "         (d-spacing %g is duplicated %i times)\n",
            Sqw->compname, list_count-mult_count, list_count-1, list[list_count-1].q, mult_count);
          );
          for (index=list_count-mult_count; index<list_count; list[index++].j = 1);
          mult_count   = 1;
          q_count = q;
          j_count = j;
          F2_count= F2;
        }
        if (i == size-1) list_count--;
        flag=0;
      }
      list_count++;
    } /* end for */

  /* now builds new Table_Array to continue with Sqw_readfile */
  if (q_max == q_min || !q_step) return(NULL);
  if (!flag_qSq)
    size = 3*q_max/q_step; /* set a default of 3 q values per line */
  else size = list_count;
  /* update the value of q_step */
  q_step = q_max/size;
  MPI_MASTER(
  if (Sqw->verbose_output > 0)
    printf("Isotropic_sqw: q range [%g:%g], creating %li elements vector\n",
        q_min, q_max, size);
  );

  retTable = (t_Table*)calloc(4, sizeof(t_Table));
  if (!retTable) printf("Isotropic_Sqw: ERROR: Cannot allocate PowderN->Sqw table.\n");
  else {
    char *header;
    if (!Table_Init(&retTable[0], size, 1))
      { printf("Isotropic_Sqw: ERROR Cannot allocate q-axis [%li] from Powder lines.\n", size); return(NULL); }
    if (!Table_Init(&retTable[1], 1, 1))
      { printf("Isotropic_Sqw: ERROR Cannot allocate w-axis from Powder lines.\n"); return(NULL); }
    if (!Table_Init(&retTable[2], size, 1))
      { printf("Isotropic_Sqw: ERROR Cannot allocate Sqw [%li] from Powder lines.\n", size); return(NULL); }
    Table_Init(&retTable[3], 0,0);

    header = malloc(64); if (header)
    { retTable[0].header = header; strcpy(retTable[0].header, "q"); }
    header = malloc(64); if (header)
    { retTable[1].header = header; strcpy(retTable[1].header, "w"); }
    header = malloc(64); if (header)
    { retTable[2].header = header; strcpy(retTable[2].header, "Sqw"); }
    for (i=0; i < 4; i++) {
      retTable[i].array_length = 3;
      retTable[i].block_number = i+1;
    }
    if (!flag_qSq)
      for (i=0; i<size; i++)
        retTable[0].data[i]  = i*q_max/size;
    for (i=0; i<list_count; i++) { /* loop on each Bragg peak */
      double peak_qmin, peak_qmax,factor,q;
      if (list[i].w > 0 && !flag_qSq) {
        peak_qmin = list[i].q*(1 - list[i].w*3);
        peak_qmax = list[i].q*(1 + list[i].w*3);
      } else { /* Dirac peak, no width */
        peak_qmin = peak_qmax = list[i].q;
      }
      /* S(q) intensity is here */
      factor = list[i].j*(list[i].DWfactor ? list[i].DWfactor : 1)
               *Sqw->mat_rho*PI/2
               /(Sqw->type == 'c' ? Sqw->s_coh : Sqw->s_inc)*list[i].F2/list[i].q/list[i].q;
      if (Sqw->barns) factor *= 100;
      for (q=peak_qmin; q <= peak_qmax; q += q_step) {
        index = (long)floor(size*q/q_max);
        if (index < 0) index=0;
        else if (index >= size) index = size-1;
        if (flag_qSq) {
          retTable[2].data[index] += list[i].F2;
          retTable[0].data[index]  = list[i].q;
        } else {
          if (list[i].w <=0 || list[i].w*q < q_step) /* step function */
            retTable[2].data[index] += factor/q_step;
          else /* gaussian */
            retTable[2].data[index] += factor
                  * Sqw_powder_gauss(q, list[i].q, list[i].w*list[i].q);
        }
      }
    } /* end for i */
    Table_Stat(&retTable[0]); Table_Stat(&retTable[1]); Table_Stat(&retTable[2]);
    Sqw->sqw_norm = 0; /* F2 are normalized already */
  }

  return(retTable);
} /* Sqw_read_PowderN */

/*****************************************************************************
*  Sqw_search_SW: For a given random number 'randnum', search for the bin
*   containing  the corresponding Sqw->SW
*  Choose an energy in the projected S(w) distribution
* Used in : TRACE (1)
*****************************************************************************/
int Sqw_search_SW(struct Sqw_Data_struct Sqw, double randnum)
{
  int index_w=0;

  if (randnum <0) randnum=0;
  if (randnum >1) randnum=1;

  if (Sqw.w_bins == 1) return(0);
  /* benefit from fast lookup table if exists */
  if (Sqw.SW_lookup) {
    index_w = Sqw.SW_lookup[(long)floor(randnum*Sqw.lookup_length)]-1;
    if (index_w<0) index_w=0;
  }

  while (index_w < Sqw.w_bins && (&(Sqw.SW[index_w]) != NULL) && (randnum > Sqw.SW[index_w].cumul_proba))
      index_w++;
  if (index_w >= Sqw.w_bins) index_w = Sqw.w_bins-1;

  if (&(Sqw.SW[index_w]) == NULL)
  {
      printf("Isotropic_Sqw: Warning: No corresponding value in the SW. randnum too big.\n");
      printf("  index_w=%i ; randnum=%f ; Sqw.SW[index_w-1].cumul_proba=%f (Sqw_search_SW)\n",
            index_w, randnum, Sqw.SW[index_w-1].cumul_proba);
      return index_w-1;
  }
  else
      return (index_w);
}

/*****************************************************************************
*  Sqw_search_Q_proba_per_w: For a given random number randnum, search for
*   the bin containing the corresponding Sqw.SW in the Q probablility grid
*  Choose a momentum in the S(q|w) distribution
*  index is given by Sqw_search_SW
* Used in : TRACE (1)
*****************************************************************************/
int Sqw_search_Q_proba_per_w(struct Sqw_Data_struct Sqw, double randnum, int index_w)
{
  int index_q=0;

  if (randnum <0) randnum=0;
  if (randnum >1) randnum=1;

  /* benefit from fast lookup table if exists */
  if (Sqw.QW_lookup && Sqw.QW_lookup[index_w]) {
    index_q = Sqw.QW_lookup[index_w][(long)floor(randnum*Sqw.lookup_length)]-1;
    if (index_q<0) index_q=0;
  }

  while (index_q < Sqw.q_bins && (&(Sqw.SQW[index_w][index_q]) != NULL)
    && (randnum > Sqw.SQW[index_w][index_q].cumul_proba)) {
      index_q++;
  }
  if (index_q >= Sqw.q_bins) index_q = Sqw.q_bins-1;

  if (&(Sqw.SQW[index_w][index_q]) == NULL)
    return -1;
  else
    return (index_q);
}

/*****************************************************************************
* compute the effective total cross section \int q S(q,w) dw dq
* for incoming neutron energy 0 < Ei < 2*w_max, and
* integration range w=-w_max:Ei and q=Q0:Q1 with
*   Q0 = SE2Q*(sqrt(Ei)-sqrt(Ei-w))=|Ki-Kf|
*   Q1 = SE2Q*(sqrt(Ei)+sqrt(Ei-w))=|Ki+Kf|
* The data to use is Sqw_Data->Sqw, and the limits are Sqw_Data->w_max Sqw_Data->q_max
*   Returns the integral value
* Used in: Sqw_readfile (1)
*****************************************************************************/
double Sqw_integrate_iqSq(struct Sqw_Data_struct *Sqw_Data, double Ei)
{
  long   index_w;
  double iqSq = 0;
  /* w=Ei-Ef q=ki-kf w>0 neutron looses energy, Stokes, Ef = Ei-w > 0, Kf =|Ki-q| > 0 */
  for (index_w=0; index_w < Sqw_Data->w_bins; index_w++) {
    long   index_q;
    double w = -Sqw_Data->w_max + index_w * Sqw_Data->w_step; /* in the Sqw table */
    if (w <= Ei) {       /* integration range w=-w_max:Ei, Ef = Ei-w > 0 */
      double sq=0, Q0=0, Q1=0;
      sq = sqrt(Ei-w);  /* always real as test was true before */
      Q0 = SE2V*V2K*(sqrt(Ei)-sq);
      Q1 = SE2V*V2K*(sqrt(Ei)+sq);

      for (index_q=0; index_q < Sqw_Data->q_bins; index_q++) {
        double q=(double)index_q * Sqw_Data->q_step;
        /* add 'pixel' = q S(q,w) */
        if (Q0 <= q && q <= Q1) iqSq += q*Table_Index(Sqw_Data->Sqw, index_q, index_w);
      }
    }
  }
  /* multiply by 'pixel' size = dq dw */
  return(iqSq * Sqw_Data->q_step * Sqw_Data->w_step);
} /* Sqw_integrate_iqSq */


/*****************************************************************************
* Sqw_diagnosis: Computes Sqw_classical, moments and physical quantities
*                make consistency checks, and output some data files
*   Return: output files and information displayed
* Used in: Sqw_init (2) only by MASTER node with MPI
*****************************************************************************/
void Sqw_diagnosis(struct Sqw_sample_struct *Sqw, struct Sqw_Data_struct *Sqw_Data)
{

  t_Table  Sqw_cl;             /* the Sqw symmetric/classical version (T-> Inf) */
  t_Table  Gqw;                /* the generalized density of states as of Carpenter and Price, J Non Cryst Sol 92 (1987) 153 */
  t_Table  Sqw_moments[7];     /* M0=S(q) M1=E_r M3 w_c w_l M0_cl=S_cl(q) G(w) */
  t_Table  w_c, w_l;
  long     index_q, index_w;
  char     c[CHAR_BUF_LENGTH]; /* temporary variable */
  long     q_min_index = 0;

  char     do_coh=0, do_inc=0;
  double   q_min =0;
  double   u2    =0, S0=1;
  long     u2_count=0;

  if (!Sqw_Data || !Sqw_Data->intensity) return; /* nothing to do with empty S(q,w) */

  if (Sqw_Data->type=='c') do_coh = 1;
  if (Sqw_Data->type=='i') do_inc = 1;

  q_min = Sqw_Data->q_min_file;
  if (q_min <= 0) q_min = Sqw_Data->q_step;

  /* test if there is only one S(q,w) available */
  if (!((Sqw->Data_inc).intensity) || !((Sqw->Data_coh).intensity))
    do_coh = do_inc = 1; /* do both if only one file given */

  if (Sqw->Temperature > 0) {
    if (!Table_Init(&Sqw_cl, Sqw_Data->q_bins, Sqw_Data->w_bins)) {
      printf("Isotropic_Sqw: %s: Cannot allocate S_cl(q,w) Table (%lix%i).\n"
             "WARNING          Skipping S(q,w) diagnosis.\n",
             Sqw->compname, Sqw_Data->q_bins, 1);
      return;
    }
    sprintf(Sqw_cl.filename,
      "S(q,w)_cl from %s (dynamic structure factor, classical)",
      Sqw_Data->filename);
    Sqw_cl.block_number = 1;
    Sqw_cl.min_x        = 0;
    Sqw_cl.max_x        = Sqw_Data->q_max;
    Sqw_cl.step_x       = Sqw_Data->q_step;
  }

  /* initialize moments and 1D stuff */
  for (index_q=0; index_q < 6; index_q++) {
    if (!Table_Init(&Sqw_moments[index_q], Sqw_Data->q_bins, 1)) {
      printf("Isotropic_Sqw: %s: Cannot allocate S(q,w) moment %ld Table (%lix%i).\n"
           "WARNING          Skipping S(q,w) diagnosis.\n",
           Sqw->compname, index_q, Sqw_Data->q_bins, 1);
      Table_Free(&Sqw_cl);
      return;
    }
    Sqw_moments[index_q].block_number = 1;
    Sqw_moments[index_q].min_x  = 0;
    Sqw_moments[index_q].max_x  = Sqw_Data->q_max;
    Sqw_moments[index_q].step_x = Sqw_Data->q_step;
  }
  index_q=6;
  Table_Init(&Sqw_moments[index_q], Sqw_Data->w_bins, 1);
  Sqw_moments[index_q].block_number = 1;
  Sqw_moments[index_q].min_x  = -Sqw_Data->w_max;
  Sqw_moments[index_q].max_x  =  Sqw_Data->w_max;
  Sqw_moments[index_q].step_x =  Sqw_Data->w_step;

  /* set Table titles */
  sprintf(Sqw_moments[0].filename,
    "S(q)=M0(q) from %s [int S(q,w) dw]",
    Sqw_Data->filename);
  sprintf(Sqw_moments[1].filename,
    "M1(q) 1-st moment from %s [int w S(q,w) dw] = HBAR^2*q^2/2/m (f-sum rule, recoil, Lovesey T1 Eq 3.63 p72, Egelstaff p196)",
    Sqw_Data->filename);
  sprintf(Sqw_moments[2].filename,
    "M3(q) 3-rd moment from %s [int w^3 S(q,w) dw] = M1(q)*w_l^2(q)",
    Sqw_Data->filename);
  sprintf(Sqw_moments[3].filename,
    "w_c(q) = sqrt(M1(q)/M0(q)*2kT) collective excitation from %s (Lovesey T1 Eq 5.38 p180, p211 Eq 5.204). Gaussian half-width of the S(q,w) classical",
    Sqw_Data->filename);
  sprintf(Sqw_moments[4].filename,
    "w_l(q) = sqrt(M3(q)/M1(q)) harmonic frequency from %s (Lovesey T1 5.39 p 180)",
    Sqw_Data->filename);
  sprintf(Sqw_moments[5].filename,
    "S_cl(q)=M0_cl(q) from %s [int S_cl(q,w) dw]",
    Sqw_Data->filename);
  sprintf(Sqw_moments[6].filename,
    "G(w) generalized effective density of states from %s (Carpenter J Non Cryst Sol 92 (1987) 153)",
    Sqw_Data->filename);

  for   (index_q=0; index_q < Sqw_Data->q_bins; index_q++) {
    double q           = index_q*Sqw_Data->q_step; /* q value in Sqw_full ; q_min = 0 */
    double sq          = 0;              /* S(q) = w0 = 0-th moment */
    double w1          = 0;              /* first  moment      \int w     Sqw dw */
    double w3          = 0;              /* third  moment      \int w^3   Sqw dw */
    double sq_cl       = 0;              /* S(q) = M0 = 0-th moment classical */
    double w_c         = 0;
    double w_l         = 0;

    for (index_w=0; index_w < Sqw_Data->w_bins; index_w++) {

      double w = -Sqw_Data->w_max + index_w*Sqw_Data->w_step; /* w value in Sqw_full */
      double sqw_cl      =0;
      double sqw_full    =0;

      sqw_full = Table_Index(Sqw_Data->Sqw, index_q, index_w);

      /* Sqw moments */
      if (w && Sqw_Data->w_bins) {
        double tmp;
        tmp  = sqw_full*Sqw_Data->w_step;
        tmp *= w;   w1 += tmp;
        tmp *= w*w; w3 += tmp;
      }

      /* compute classical Sqw and S(q)_cl */
      if (Sqw->Temperature > 0) {
        double n;
        sqw_cl = sqw_full * Sqw_quantum_correction(-w,Sqw->Temperature,Sqw->Q_correction);
        if (!Table_SetElement(&Sqw_cl, index_q, index_w, sqw_cl))
          printf("Isotropic_Sqw: %s: "
                 "Error when setting Sqw_cl[%li q=%g,%li w=%g]=%g from file %s\n",
                 Sqw->compname, index_q, q, index_w, w, sqw_cl, Sqw_Data->filename);
        sq_cl += sqw_cl;
      }
      sq    += sqw_full;
    } /* for index_w */

    sq    *= Sqw_Data->w_step;         /* S(q) = \int S(q,w) dw = structure factor */
    sq_cl *= Sqw_Data->w_step;
    /* find minimal reliable q value (not interpolated) */
    if (q >= q_min && !q_min_index && sq) {
      q_min_index = index_q;
      q_min       = q;
      if (0.9 < sq)
        S0          = sq; /* minimum reliable S(q) */
      else S0 = 1;
    }
    /* compute <u^2> = <3 * ln(S(q)) / q^2> */
    if (q_min_index && q && S0 && sq) {
      u2      += 3 * log(sq/S0) /q/q;
      u2_count++;
    }

    /* store moment values (q) as M0=S(q) M1=E_r M3 w_c w_l M0_cl=S_cl(q) */
    Table_SetElement(&Sqw_moments[0],    index_q, 0, sq);
    Table_SetElement(&Sqw_moments[1],    index_q, 0, w1);
    Table_SetElement(&Sqw_moments[2],    index_q, 0, w3);
    if (w1 > 0 && sq && Sqw->Temperature > 0) {
      double w_c = sqrt(w1/sq*2*Sqw->Temperature*Sqw->T2E);  /* HBAR^2 q^2 kT /m/ S(q) */
      Table_SetElement(&Sqw_moments[3],    index_q, 0, w_c); /* collective dispersion */
    }
    if (w1 && w3*w1 > 0) {
      double w_l = sqrt(w3/w1);
      Table_SetElement(&Sqw_moments[4],    index_q, 0, w_l); /* harmonic dispersion */
    }
    if (Sqw->Temperature > 0)
      Table_SetElement(&Sqw_moments[5],    index_q, 0, sq_cl);

  } /* for index_q */



  /* display some usefull information, only once in MPI mode (MASTER) */
  if (Sqw->Temperature > 0) {
    double Da         = 1.660538921e-27;  /* [kg] unified atomic mass unit = Dalton = 1 g/mol */
  #ifndef KB
    double KB         = 1.3806503e-23;    /* [J/K] */
    /* HBAR   = 1.05457168e-34 */
  #endif
    double CELE       = 1.602176487e-19;  /* [C] Elementary charge CODATA 2006 'e' */
    double meV2Hz     = CELE/HBAR/1000/2/PI; /* 1 meV = 241.80e9 GHz */
    double gqw_sum    = 0;

    /* classical Sqw */
    sprintf(c, "%s_%s_cl.sqw", Sqw->compname, Sqw_Data->type=='c' ? "coh" : "inc");
    Table_Write(Sqw_cl, c, "Momentum [Angs-1]", "'S(q,w)*exp(hw/2kT) classical limit' Energy [meV]",
        0,Sqw_Data->q_max,-Sqw_Data->w_max,Sqw_Data->w_max);
    Table_Free(&Sqw_cl);

    if (u2_count) u2 /= u2_count;

    MPI_MASTER(
    if (do_coh || do_inc)
      printf("Isotropic_Sqw: %s: "
             "Physical constants from the S(q,w) %s for T=%g [K]. Values are estimates.\n",
             Sqw->compname, Sqw_Data->filename, Sqw->Temperature);
    if (do_coh) {
      if (Sqw->mat_weight) {
        double LAMBDA     = HBAR*2*PI/sqrt(2*PI*Sqw->mat_weight*Da*KB*Sqw->Temperature)*1e10;   /* in [Angs] */
        double z          = Sqw->mat_rho * LAMBDA*LAMBDA*LAMBDA;  /* fugacity , rho=N/V in [Angs-3]*/
        double mu         = KB*Sqw->Temperature*log(z);       /* perfect gas chemical potential */
        printf("# De Broglie wavelength     LAMBDA=%g [Angs]\n", LAMBDA);
        printf("# Fugacity                       z=%g (from Egelstaff p32 Eq 2.31)\n", z);
        printf("# Chemical potential            mu=%g [eV] (eq. perfect gas)\n", mu/CELE);
      }

      /* compute isothermal sound velocity and compressibility */
      /* get the S(q_min) value and the corresponding w_c */

      if (q_min_index > 0 && q_min && q_min < 0.6) {
        double w_c = Table_Index(Sqw_moments[3], q_min_index, 0); /* meV */
        /* HBAR = [J*s] */
        double c_T = 2*PI*w_c*meV2Hz/q_min/1e10;                  /* meV*Angs -> m/s */
        double ChiT= S0/(KB*Sqw->Temperature*Sqw->mat_rho*1e30);
        printf("# Isothermal compressibility Chi_T=%g [Pa-1] (Egelstaff  p201 Eq 10.21) at q=%g [Angs-1]\n",
          ChiT, q_min);
        printf("# Isothermal sound velocity    c_T=%g [m/s]  (Lovesey T1 p210 Eq 5.197) at q=%g [Angs-1]\n",
          c_T, q_min);

        /* Computation if C11 is rather tricky as it is obtained from w_l, which is usually quite noisy
         * This means that the obtained values are not reliable from C = rho c_l^2 (Egelstaff Eq 14.10b p284)
         * C44 = rho c_c^2 ~ C11/3
         */
        double w_l = Table_Index(Sqw_moments[4], q_min_index, 0); /* meV */
        double c_l = 2*PI*w_l*meV2Hz/q_min/1e10;                  /* meV*Angs -> m/s */
        double C11 = (Sqw->mat_weight*Da)*(Sqw->mat_rho*1e30)*c_l*c_l;
        printf("# Elastic modulus              C11=%g [GPa]  (Egelstaff Eq 14.10b p284) [rough estimate] at q=%g [Angs-1]\n",
            C11/1e9, q_min);
      }
    }
    if (do_inc) {
      /* display the mean square displacement from S(q) = exp(-<u^2>q^2/3)
           <u^2>= <3 * ln(S(q)) / q^2>
       */
      if (u2_count && u2) {
        printf("# Mean square displacement   <u^2>=%g [Angs^2] (<3 * ln(S(q)) / q^2>)\n", u2);
      }

      /* compute the mean diffusion coefficient D=w_c/q^2 */
      /* FWHM of gaussian is Gamma*RMS2FWHM, only in diffusive regime (Q < 0.2 Angs-1) */
      if (q_min_index > 0 && q_min && q_min < 0.6) {
        double w_c = Table_Index(Sqw_moments[3], q_min_index, 0);
        double D   = 2*PI*w_c*meV2Hz/q_min/q_min/1e14*RMS2FWHM/2; /* meV*Angs^2 -> mm^2/s */
        printf("# Diffusion coefficient          D=%g [mm^2/s] (Egelstaff p220)\n", D);
        if (u2_count && u2 && D)
          printf("# Jump relaxation time         tau=%g [ns] (Egelstaff Eq 11.8 p220)\n", u2*1e-2/6/D);
      }
    }
    ); /* MPI_MASTER */

    /* density of states (generalized) */
    if (!Table_Init(&Gqw, Sqw_Data->q_bins, Sqw_Data->w_bins)) {
      printf("Isotropic_Sqw: %s: Cannot allocate G(q,w) Table (%lix%i).\n"
             "WARNING          Skipping S(q,w) diagnosis.\n",
             Sqw->compname, Sqw_Data->q_bins, 1);
        return;
    }
    sprintf(Gqw.filename,
      "G(q,w) from %s (generalized density of states, Carpenter J Non Cryst Sol 92 (1987) 153)",
      Sqw_Data->filename);
    Gqw.block_number = 1;
    Gqw.min_x        = 0;
    Gqw.max_x        = Sqw_Data->q_max;
    Gqw.step_x       = Sqw_Data->q_step;

    for (index_w=0; index_w < Sqw_Data->w_bins; index_w++) {
      double w        = -Sqw_Data->w_max + index_w*Sqw_Data->w_step; /* w value in Sqw_full */
      double gw       = 0;
      for   (index_q=0; index_q < Sqw_Data->q_bins; index_q++) {
        double q        = index_q*Sqw_Data->q_step; /* q value in Sqw_full ; q_min = 0 */
        double sqw_full = Table_Index(Sqw_Data->Sqw, index_q, index_w);
        double n        = 1/(exp(w/(Sqw->Temperature*Sqw->T2E))-1); /* Bose factor */
        double DW       = q && u2 ? exp(2*u2*q*q/6) : 1;            /* Debye-Weller factor */
        double gqw      = q && n+1 ? sqw_full*DW*2*(Sqw->mat_weight*Da)*w/(n+1)/q/q : 0;
        if (!Table_SetElement(&Gqw, index_q, index_w, gqw))
          printf("Isotropic_Sqw: %s: "
                 "Error when setting Gqw[%li q=%g,%li w=%g]=%g from file %s\n",
                 Sqw->compname, index_q, q, index_w, w, gqw, Sqw_Data->filename);
        gw      += gqw;
        gqw_sum += gqw;
      }
      Table_SetElement(&Sqw_moments[6],    index_w, 0, gw);
    }

    /* normalize the density of states */
    for (index_w=0; index_w < Sqw_Data->w_bins; index_w++) {
      double gw = Table_Index(Sqw_moments[6], index_w, 0);
      Table_SetElement(&Sqw_moments[6], index_w, 0, gw / gqw_sum);
      for   (index_q=0; index_q < Sqw_Data->q_bins; index_q++) {
        double gqw = Table_Index(Gqw, index_q, index_w);
        Table_SetElement(&Gqw, index_q, index_w, gqw / gqw_sum);
      }
    }

    /* write Gqw and free memory */
    if (Sqw_Data->w_bins > 1) {
      sprintf(c, "%s_%s.gqw", Sqw->compname, Sqw_Data->type=='c' ? "coh" : "inc");
        Table_Write(Gqw, c, "Momentum [Angs-1]", "'Generalized density of states' Energy [meV]",
          0,Sqw_Data->q_max,-Sqw_Data->w_max,Sqw_Data->w_max);
      Table_Free(&Gqw);
    }
  } /* if T>0 */

  /* write all tables to disk M0=S(q) M1=E_r M3 w_c w_l M0_cl=S_cl(q) */
  if (Sqw_Data->w_bins > 1) {
    sprintf(c, "%s_%s.m1",  Sqw->compname, Sqw_Data->type=='c' ? "coh" : "inc");
    Table_Write(Sqw_moments[1], c, "Momentum [Angs-1]", "int w S(q,w) dw (recoil) q^2/2m [meV]",
      0,Sqw_Data->q_max,0,0);
    sprintf(c, "%s_%s.w_l", Sqw->compname, Sqw_Data->type=='c' ? "coh" : "inc");
    Table_Write(Sqw_moments[4], c, "Momentum [Angs-1]", "w_l(q) harmonic frequency [meV]",
      0,Sqw_Data->q_max,0,0);
    sprintf(c, "%s_%s.sqw", Sqw->compname, Sqw_Data->type=='c' ? "coh" : "inc");
    Table_Write(Sqw_Data->Sqw, c, "Momentum [Angs-1]", "'S(q,w) dynamical structure factor [meV-1]' Energy [meV]",
      0,Sqw_Data->q_max,-Sqw_Data->w_max,Sqw_Data->w_max);

    if (Sqw->Temperature > 0) {
      sprintf(c, "%s_%s.w_c",    Sqw->compname, Sqw_Data->type=='c' ? "coh" : "inc");
      Table_Write(Sqw_moments[3], c, "Momentum [Angs-1]", "w_c(q) collective excitation [meV]", 0,Sqw_Data->q_max,0,0);
      sprintf(c, "%s_%s_cl.sq",  Sqw->compname, Sqw_Data->type=='c' ? "coh" : "inc");
      Table_Write(Sqw_moments[5], c, "Momentum [Angs-1]", "int S_cl(q,w) dw",
        0,Sqw_Data->q_max,0,0);
      sprintf(c, "%s_%s.gw",  Sqw->compname, Sqw_Data->type=='c' ? "coh" : "inc");
      Table_Write(Sqw_moments[6], c, "Energy [meV]", "'Generalized effective density of states' Energy [meV]",
        -Sqw_Data->w_max,Sqw_Data->w_max,0,0);

    }
  }
  sprintf(c, "%s_%s.sq",    Sqw->compname, Sqw_Data->type=='c' ? "coh" : "inc");
  Table_Write(Sqw_moments[0], c, "Momentum [Angs-1]","S(q) = int S(q,w) dw", 0,Sqw_Data->q_max,0,0);
  sprintf(c, "%s_%s.sigma", Sqw->compname, Sqw_Data->type=='c' ? "coh" : "inc");
  Table_Write(Sqw_Data->iqSq, c, "Energy [meV]", "sigma kf/ki int q S(q,w) dw scattering cross section [barns]", 0,0,0,0);

  /* free Tables */
  for (index_q=0; index_q < 7; Table_Free(&Sqw_moments[index_q++]));

} /* Sqw_diagnosis */

/*****************************************************************************
* Sqw_readfile: Read Sqw data files
*   Returns Sqw_Data_struct or NULL in case of error
* Used in : Sqw_init (2)
*****************************************************************************/
struct Sqw_Data_struct *Sqw_readfile(
  struct Sqw_sample_struct *Sqw, char *file, struct Sqw_Data_struct *Sqw_Data)
{

  t_Table *Table_Array= NULL;
  long     nblocks    = 0;
  char     flag       = 0;

  t_Table  Sqw_full, iqSq; /* the Sqw (non symmetric) and total scattering X section */

  double   sum=0;
  double   mat_at_nb=1;
  double   iq2Sq=0;
  long    *SW_lookup=NULL;
  long   **QW_lookup=NULL;
  char   **parsing  =NULL;

  long   index_q, index_w;
  double q_min_file, q_max_file, q_step_file;
  long   q_bins_file;
  double w_min_file, w_max_file, w_step_file;
  long   w_bins_file;
  double q_max, q_step;
  long   q_bins;
  double w_max, w_step;
  long   w_bins;

  double alpha=0;

  double M1          = 0;
  double M1_cl       = 0;
  double T           = 0;
  double T_file      = 0;
  long   T_count     = 0;
  long   M1_count    = 0;
  long   M1_cl_count = 0;

  /* setup default */
  Sqw_Data_init(Sqw_Data);

  if (!file || !strlen(file) || !strcmp(file, "NULL") || !strcmp(file, "0")) return(Sqw_Data);
  /* read the Sqw file */
  Table_Array = Table_Read_Array(file, &nblocks);
  strncpy(Sqw_Data->filename, file, 80);
  if (!Table_Array) return(NULL);

  /* (1) parsing of header ================================================== */
  parsing = Table_ParseHeader(Table_Array[0].header,
    "Vc","V_0",
    "sigma_abs","sigma_a ",
    "sigma_inc","sigma_i ",
    "column_j", /* 6 */
    "column_d",
    "column_F2",
    "column_DW",
    "column_Dd",
    "column_inv2d", "column_1/2d", "column_sintheta_lambda",
    "column_q", /* 14 */
    "sigma_coh","sigma_c ",
    "Temperature",
    "column_Sq",
    "column_F ", /* 19 */
    "V_rho",
    "density",
    "weight",
    "nb_atoms","multiplicity",
    "classical",
    NULL);
  if (parsing) {
    int i;
    if (parsing[0] && !Sqw->mat_rho)      Sqw->mat_rho    =1/atof(parsing[0]);
    if (parsing[1] && !Sqw->mat_rho)      Sqw->mat_rho    =1/atof(parsing[1]);
    if (parsing[2] && !Sqw->s_abs)    Sqw->s_abs  =  atof(parsing[2]);
    if (parsing[3] && !Sqw->s_abs)    Sqw->s_abs  =  atof(parsing[3]);
    if (parsing[4] && !Sqw->s_inc)    Sqw->s_inc  =  atof(parsing[4]);
    if (parsing[5] && !Sqw->s_inc)    Sqw->s_inc  =  atof(parsing[5]);
    if (parsing[6])                   Sqw->column_order[0]=atoi(parsing[6]);
    if (parsing[7])                   Sqw->column_order[1]=atoi(parsing[7]);
    if (parsing[8])                   Sqw->column_order[2]=atoi(parsing[8]);
    if (parsing[9])                   Sqw->column_order[3]=atoi(parsing[9]);
    if (parsing[10])                  Sqw->column_order[4]=atoi(parsing[10]);
    if (parsing[11])                  Sqw->column_order[5]=atoi(parsing[11]);
    if (parsing[12])                  Sqw->column_order[5]=atoi(parsing[12]);
    if (parsing[13])                  Sqw->column_order[5]=atoi(parsing[13]);
    if (parsing[14])                  Sqw->column_order[6]=atoi(parsing[14]);
    if (parsing[15] && !Sqw->s_coh)   Sqw->s_coh=atof(parsing[15]);
    if (parsing[16] && !Sqw->s_coh)   Sqw->s_coh=atof(parsing[16]);
    if (parsing[17] && !Sqw->Temperature) Sqw->Temperature=atof(parsing[17]); /* from user or file */
    if (parsing[17] )                 T_file=atof(parsing[17]); /* from file */
    if (parsing[18])                  Sqw->column_order[8]=atoi(parsing[18]);
    if (parsing[19])                  Sqw->column_order[7]=atoi(parsing[19]);
    if (parsing[20] && !Sqw->mat_rho)     Sqw->mat_rho    =atof(parsing[20]);
    if (parsing[21] && !Sqw->mat_density) Sqw->mat_density=atof(parsing[21]);
    if (parsing[22] && !Sqw->mat_weight)  Sqw->mat_weight =atof(parsing[22]);
    if (parsing[23] )                 mat_at_nb   =atof(parsing[23]);
    if (parsing[24] )                 mat_at_nb   =atof(parsing[24]);
    if (parsing[25] )                 { /* classical is found in the header */
      char *endptr;
      double value = strtod(parsing[25], &endptr);
      if (*endptr == *parsing[25]) {
        if (Sqw->sqw_classical < 0) Sqw->sqw_classical = 1;
      } else                        Sqw->sqw_classical = value;
    }
    for (i=0; i<=25; i++) if (parsing[i]) free(parsing[i]);
    free(parsing);
  }

  /* compute the scattering unit density from material weight and density */
  /* the weight of the scattering element is the chemical formula molecular weight
   * times the nb of chemical formulae in the scattering element (nb_atoms) */
  if (!Sqw->mat_rho && Sqw->mat_density > 0 && Sqw->mat_weight > 0 && mat_at_nb > 0) {
    /* molar volume [cm^3/mol] = weight [g/mol] / density [g/cm^3] */
    /* atom density per Angs^3 = [mol/cm^3] * N_Avogadro *(1e-8)^3 */
    Sqw->mat_rho = Sqw->mat_density/(Sqw->mat_weight*mat_at_nb)/1e24*NA;
    MPI_MASTER(
    if (Sqw->verbose_output > 0)
      printf("Isotropic_Sqw: %s: Computing scattering unit density V_rho=%g [AA^-3] from density=%g [g/cm^3] weight=%g [g/mol].\n",
        Sqw->compname, Sqw->mat_rho, Sqw->mat_density, Sqw->mat_weight);
    );
  }

  /* the scattering unit cross sections are the chemical formula ones
   * times the nb of chemical formulae in the scattering element */
  if (mat_at_nb > 0) {
    Sqw->s_abs *= mat_at_nb; Sqw->s_inc *= mat_at_nb; Sqw->s_coh *= mat_at_nb;
  }

  if (nblocks) {
    if (nblocks == 1) {
      /* import Powder file */
      t_Table *newTable   = NULL;
      newTable = Sqw_read_PowderN(Sqw, Table_Array[0]);
      if (!newTable) {
        MPI_MASTER(
        printf("Isotropic_Sqw: %s: ERROR importing powder line file %s.\n"
               "               Check format definition.\n",
              Sqw->compname, file);
        );
        exit(-1);
      } else flag=0;
      Table_Free_Array(Table_Array);
      Table_Array = newTable;
    } else if (nblocks != 3) {
      MPI_MASTER(
      printf("Isotropic_Sqw: %s: ERROR "
             "File %s contains %li block%s instead of 3.\n",
              Sqw->compname, file, nblocks, (nblocks == 1 ? "" : "s"));
      );
    } else { flag=0; Sqw->barns=0; /* Sqw files do not use powder_barns */ }
  }

  /* print some info about Sqw files */
  if (flag) Sqw->verbose_output = 2;

  if (flag) {
    MPI_MASTER(
    if (nblocks) printf("ERROR          Wrong file format.\n"
      "               Disabling contribution.\n"
      "               File must contain 3 blocks for [q,w,sqw] or Powder file (1 block, laz,lau).\n");
    );
    return(Sqw_Data);
  }

  sprintf(Table_Array[0].filename, "%s#q",   file);
  sprintf(Table_Array[1].filename, "%s#w",   file);
  sprintf(Table_Array[2].filename, "%s#sqw", file);

  MPI_MASTER(
  if (nblocks && Sqw->verbose_output > 2) {
    printf("Isotropic_Sqw: %s file read, analysing...\n", file);
    Table_Info_Array(Table_Array);
  }
  );

  /* (2) compute range for full +/- w and allocate S(q,w) =================== */

  /* get the q,w extend of the table from the file */
  q_bins_file = Table_Array[0].rows*Table_Array[0].columns;
  w_bins_file = Table_Array[1].rows*Table_Array[1].columns;

  /* is there enough qw data in file to proceed ? */
  if (q_bins_file <= 1 || w_bins_file <= 0) {
    MPI_MASTER(
    printf("Isotropic_Sqw: %s: Data file %s has incomplete q or omega information (%lix%li).\n"
           "ERROR          Exiting.\n",
      Sqw->compname, file, q_bins_file, w_bins_file);
    );
    return(Sqw_Data);
  }

  q_min_file  = Table_Array[0].min_x; q_max_file = Table_Array[0].max_x;
  q_step_file = Table_Array[0].step_x ? Table_Array[0].step_x : (q_max_file - q_min_file)/(Table_Array[0].rows*Table_Array[0].columns);
  w_min_file  = Table_Array[1].min_x; w_max_file = Table_Array[1].max_x;
  w_step_file = Table_Array[1].step_x;

  /* create a regular extended q,w and Sqw tables applying the exp(-hw/kT) factor */
  q_max  = q_max_file;
  q_bins = (q_step_file ?   q_max/q_step_file : q_bins_file)+1;
  q_step = q_bins-1 > 0 ?   q_max/(q_bins-1) : 1;
  w_max  = fabs(w_max_file);
  if (fabs(w_min_file) > fabs(w_max_file)) w_max = fabs(w_min_file);
  /* w_min =-w_max */
  w_bins = (w_step_file ? (long)(2*w_max/w_step_file) : 0)+1; /* twice the initial w range */
  w_step = w_bins-1 > 0 ? 2*w_max/(w_bins-1) : 1;             /* that is +/- w_max         */

  /* create the Sqw table in full range */
  if (!Table_Init(&Sqw_full, q_bins, w_bins)) {
    printf("Isotropic_Sqw: %s: Cannot allocate Sqw_full Table (%lix%li).\n"
           "ERROR          Exiting.\n",
      Sqw->compname, q_bins, w_bins);
    return(NULL);
  }
  sprintf(Sqw_full.filename, "S(q,w) from %s (dynamic structure factor)", file);
  Sqw_full.block_number = 1;

  Sqw_Data->q_bins = q_bins; Sqw_Data->q_max = q_max; Sqw_Data->q_step= q_step;
  Sqw_Data->w_bins = w_bins; Sqw_Data->w_max = w_max; Sqw_Data->w_step= w_step;
  Sqw_Data->q_min_file = q_min_file;

  /* build an energy symmetric Sqw data set with detailed balance there-in, so
   * that we can both compute effective scattering Xsection, probability distributions
   * that is S(q) and \int q S(q).
   * We scan the new Sqw table elements with regular qw binning and search for their
   * equivalent element in the Sqw file data set. This is slower than doing the opposite.
   * We could be scanning all file elements, and fill the new table, but in the
   * process some empty spaces may appear when the initial file binning is not regular
   * in qw, leading to gaps in the new table.
   */

  /* (3) we build q and w lookup table for conversion file -> sqw_full ====== */
  MPI_MASTER(
  if (Sqw->verbose_output > 2)
    printf("Isotropic_Sqw: %s: Creating Sqw_full... (%s, %s)\n",
      Sqw->compname, file, Sqw->type=='c' ? "coh" : "inc");
  );

  double w_file2full[w_bins]; /* lookup table for fast file -> Sqw_full allocation */

  for (index_w=0; index_w < w_bins; w_file2full[index_w++]=0);

  for (index_w=0; index_w < w_bins; index_w++) {

    double w = -w_max + index_w*w_step; /* w value in Sqw_full */
    double index_w_file=0;              /* w index in Sqw file */
    char   found=0;
    for (index_w_file=0; index_w_file < w_bins_file; index_w_file++) {
      double w0=Table_Index(Table_Array[1], (long)index_w_file,  0);
      double w1=Table_Index(Table_Array[1], (long)index_w_file+1,0);
      /* test if we are in Stokes */
      if (w0 > w1) {
        double tmp=w0; w0=w1; w1=tmp;
      }
      if (w0 <= w && w < w1) {
        /* w ~ w_file exists in file, usually on w > 0 side Stokes, neutron looses energy */
        index_w_file += w1-w0 ? (w-w0)/(w1-w0) : 0; /* may correspond with a position in-betwwen two w elements */
        found=1;
        break;
      }
    }
    /* test if we are in anti-Stokes */
    if (!found)
    for (index_w_file=0; index_w_file < w_bins_file; index_w_file++) {
      double w0=Table_Index(Table_Array[1], (long)index_w_file,  0);
      double w1=Table_Index(Table_Array[1], (long)index_w_file+1,0);
      /* test if we are in anti-Stokes */
      if (w0 > w1) {
        double tmp=w0; w0=w1; w1=tmp;
      }
      if (w0 <= -w && -w < w1) {            /* w value is mirrored from the opposite side in file */
        index_w_file += w1-w0 ? (-w-w0)/(w1-w0) : 0;
        index_w_file = -index_w_file;               /* in this case, index value is set to negative */
        break;
      }
    }
    w_file2full[index_w] = index_w_file;
  }

  double q_file2full[q_bins];
  for   (index_q=0; index_q < q_bins; q_file2full[index_q++]=0);

  for (index_q=0; index_q < q_bins; index_q++) {

    double q           = index_q*q_step; /* q value in Sqw_full ; q_min = 0 */
    double index_q_file= 0;              /* q index in Sqw file */

    /* search for q value in the initial file data set */
    if      (q <= q_min_file) index_q_file=0;
    else if (q >= q_max_file) index_q_file=q_bins_file-1;
    else
    for (index_q_file=0; index_q_file < q_bins_file; index_q_file++) {
      double q0=Table_Index(Table_Array[0], (long)index_q_file,  0);
      double q1=Table_Index(Table_Array[0], (long)index_q_file+1,0);
      if (q0 <= q && q <= q1) {
        index_q_file += q1-q0 ? (q-q0)/(q1-q0) : 0; /* may correspond with a position in-betwwen two q elements */
        break;
      }
    }
    q_file2full[index_q] = index_q_file;
  }

  /* (4) now we build Sqw on full Q,W ranges, using the Q,W table lookup above -> Sqw_full */
  for (index_q=0; index_q < q_bins; index_q++) {

    double q           = index_q*q_step; /* q value in Sqw_full ; q_min = 0 */
    double index_q_file= 0;              /* q index in Sqw file */

    /* get q value in the initial file data set */
    index_q_file = q_file2full[index_q];

    /* now scan energy elements in Sqw full, and search these in file data */
    for (index_w=0; index_w < w_bins; index_w++) {
      double w = -w_max + index_w*w_step; /* w value in Sqw_full */
      double index_w_file=0;              /* w index in Sqw file */
      double sqw_file    =0;              /* Sqw(index_q, index_w) value interpolated from file */

      /* search for w value in the file data set, negative when mirrored */
      index_w_file = w_file2full[index_w];
      /* get Sqw_file element, with bi-linear interpolation from file */
      /* when the initial file does not contain the energy, the opposite element (-w) is used */
      sqw_file     = Table_Value2d(Table_Array[2], index_q_file, abs(index_w_file));
      /* apply the minimum threshold to remove noisy background in S(q,w) */
      if (sqw_file < Sqw->sqw_threshold) sqw_file = 0;
      else if (index_w_file < 0)         sqw_file = -sqw_file; /* negative == mirrored from other side */

      if (!Table_SetElement(&Sqw_full, index_q, index_w, sqw_file))
        printf("Isotropic_Sqw: %s: "
               "Error when setting Sqw[%li q=%g,%li w=%g]=%g from file %s\n",
               Sqw->compname, index_q, q, index_w, w, fabs(sqw_file), file);
    } /* for index_w */
  } /* for index_q */

  /* free memory and store limits for new full Sqw table */
  Table_Free_Array(Table_Array);

  /* if only one S(q,w) side is given, it is symmetrised by mirroring, then M1=0 */

  /* (5) test if the Sqw_full is classical or not by computing the 1st moment (=0 for classical) */
  /* also compute temperature (quantum case) from file if not set */
  for (index_q=0; index_q < q_bins; index_q++) {

    double q           = index_q*q_step; /* q value in Sqw_full ; q_min = 0 */

    for (index_w=0; index_w < w_bins; index_w++) {
      double w        = -w_max + index_w*w_step; /* w value in Sqw_full */
      double sqw_full = Table_Index(Sqw_full, index_q, index_w);
      long   index_mw = w_bins-1-index_w;        /* opposite w index in S(q,w) */
      double sqw_opp  = Table_Index(Sqw_full, index_q, index_mw);
      double T_defined= T_file ? T_file : Sqw->Temperature; /* T better from file, else from user */

      /* the analysis must be done only on values which exist on both sides */
      /* as integrals must be symmetric, and Bose factor requires both sides as well */
      if (sqw_full > 0 && sqw_opp > 0) {
        /* compute temperature from Bose factor */
        if (sqw_opp != sqw_full) {
          T      += fabs(w/log(sqw_opp/sqw_full)/Sqw->T2E);
          T_count++;
        }
        /* we first assume Sqw is quantum. M1_cl should be 0, M1 should be recoil */
        M1      += w*sqw_full*w_step;
        M1_count++;
        /* we assume it is quantum (non symmetric) and check that its symmetrized version has M1_cl=0 */
        if (T_defined > 0) {
          sqw_opp = sqw_full * Sqw_quantum_correction(-w, T_defined,Sqw->Q_correction); /* Sqw_cl */
          M1_cl      += w*sqw_opp*w_step;
          M1_cl_count++;
        } else if (Sqw->mat_weight) {
          /* T=0 ? would compute the M1_cl = M1 - recoil energy, assuming we have a quantum S(q,w) in file */
          /* the M1(quantum) = (MNEUTRON/m)*2.0725*q^2 recoil energy */
          double Da = 1.660538921e-27; /* atomic mass unit */
          double Er = (MNEUTRON/Sqw->mat_weight/Da)*2.0725*q*q; /* recoil for one scattering unit in the cell [meV] Schober JDN16 p239 */
          M1_cl      += M1 - Er;
          M1_cl_count++;
        }
      } /* both side from file */
    } /*index_w */
  } /*index_q */

  if (T_count)     T     /= T_count;     /* mean temperature from Bose ratio */
  if (M1_count)    M1    /= M1_count;
  if (M1_cl_count) M1_cl /= M1_cl_count; /* mean energy value along q range */

  /* determine if we use a classical or quantum S(q,w) */
  if (Sqw->sqw_classical < 0) {
    if (fabs(M1) < 2*w_step) {
      Sqw->sqw_classical = 1; /* the initial Sqw from file seems to be centered, thus classical */
    } else if (fabs(M1_cl) < fabs(M1)) {
      /* M1 for classical is closer to 0 than for quantum one */
      Sqw->sqw_classical = 0; /* initial data from file seems to be quantum (non classical) */
    } else { /* M1_cl > M1 > 2*w_step */
      MPI_MASTER(
      printf("Isotropic_Sqw: %s: I do not know if S(q,w) data is classical or quantum.\n"
             "WARNING        First moment M1=%g M1_cl=%g for file %s. Defaulting to classical case.\n",
                 Sqw->compname, M1, M1_cl, file);
      );
    }
  }
  if (Sqw->sqw_classical < 0) Sqw->sqw_classical=1; /* default when quantum/classical choice is not set */

  /* (5b) set temperature. Temperatures defined are:
  *   T_file:           temperature read from the .sqw file
  *   T:                temperature computed from the quantum Sqw using detailed balance
  *   Sqw->Temperature: temperature defined by user, or read from file when not set
  */


  /* display some warnings about the computed temperature */
  if (T > 3000) T=0; /* unrealistic */
  if (!T_file && T) {
    MPI_MASTER(
    if (Sqw->verbose_output > 0) {
      printf(  "Isotropic_Sqw: %s: Temperature computed from S(q,w) data from %s is T=%g [K].\n",
        Sqw->compname, file, T);
    );
    }
  }

  if (Sqw->Temperature == 0) {
    Sqw->Temperature = T_file ? T_file : T; /* 0:  not set: we use file value, else computed */
  } else if (Sqw->Temperature ==-1) {
    Sqw->Temperature = 0;                   /* -1: no detailed balance. Display message at end of INIT */
  } else if (Sqw->Temperature ==-2) {
    Sqw->Temperature = T ? T : T_file;      /* -2: use guessed value when available */
  } /* else let value as it is (e.g. >0) */

  if (Sqw->verbose_output > 0 && Sqw->Temperature) {
    MPI_MASTER(
    printf(  "Isotropic_Sqw: %s: Temperature set to T=%g [K]\n", Sqw->compname, Sqw->Temperature);
    );
  }

  MPI_MASTER(
  if (Sqw->verbose_output > 0 && w_bins > 1)
    printf("Isotropic_Sqw: %s: S(q,w) data from %s (%s) assumed to be %s.\n",
      Sqw->compname, file, Sqw->type=='c' ? "coh" : "inc",
      Sqw->sqw_classical ? "classical (symmetrised in energy)" : "non-classical (includes Bose factor, non symmetric in energy)");
  );

  /* (6) apply detailed balance on Sqw_full for classical or quantum S(q,w) */
  /* compute the \int q^2 S(q) for normalisation */

  MPI_MASTER(
  if (Sqw->sqw_classical && Sqw->verbose_output > 0 && Sqw->Temperature > 0)
    printf("Isotropic_Sqw: %s: Applying exp(hw/2kT) factor with T=%g [K] on %s file (classical/symmetric) using '%s' quantum correction\n",
      Sqw->compname, Sqw->Temperature, file, Sqw->Q_correction);
  );
  for   (index_q=0; index_q < q_bins; index_q++) {
    double sq          = 0;
    double q           = index_q*q_step;  /* q value in Sqw_full ; q_min = 0 */
    for (index_w=0; index_w < w_bins; index_w++) {
      double w = -w_max + index_w*w_step; /* w value in Sqw_full */
      double balance   = 1;               /* detailed balance factor, default is 1 */
      double sqw_full  = Table_Index(Sqw_full, index_q, index_w);

      /* do we use a symmetric S(q,w) from real G(r,t) from e.g. MD ? */

      if (Sqw->sqw_classical && Sqw->Temperature > 0) /* data is symmetric, we apply Bose factor */
        /* un-symmetrize Sqw(file) * exp(hw/kT/2) on both sides */
        balance = Sqw_quantum_correction(w, Sqw->Temperature, Sqw->Q_correction);
      else if (!Sqw->sqw_classical) {  /* data is quantum (contains Bose) */
        if (sqw_full < 0) { /* quantum but mirrored/symmetric value (was missing in file) */
          if (T)
            /* prefer to use T computed from file for mirroring */
            balance *= exp(w/(T*Sqw->T2E));                /* apply Bose where missing */
          else if (Sqw->Temperature)
            balance *= exp(w/(Sqw->Temperature*Sqw->T2E)); /* apply Bose where missing */
        }
        /* test if T computed from file matches requested T, else apply correction */
        if (T && Sqw->Temperature > 0 && Sqw->Temperature != T) {
          balance *= exp(-w/(T*Sqw->T2E)/2);                /* make it a classical data set: remove computed/read T from quantum data file */
          balance *= exp( w/(Sqw->Temperature*Sqw->T2E)/2); /* then apply Bose to requested temperature */
        }
      }

      /* update Sqw value */
      sqw_full = fabs(sqw_full)*balance;
      Table_SetElement(&Sqw_full, index_q, index_w, sqw_full);
      /* sum up the S(q) (0-th moment) = w0 */
      sq       += sqw_full;
    } /* index_w */
    sq    *= w_step;         /* S(q)  = \int S(q,w) dw    = structure factor */
    iq2Sq += q*q*sq*q_step;  /* iq2Sq = \int q^2 S(q) dq  = used for g-sum rule (normalisation) */
    sum   += sq*q_step;      /* |S|   = \int S(q,w) dq dw = total integral value in file */
  } /* index_q */

  if (!sum) {
    MPI_MASTER(
    printf("Isotropic_Sqw: %s: No valid data in the selected (Q,w) range: sum(S)=0\n"
           "ERROR          Available Sqw data is\n",
      Sqw->compname);
    printf("                 q=[%g:%g] w=[%g:%g]\n",
           q_min_file, q_max_file,
           w_min_file, w_max_file);
    );
    Table_Free(&Sqw_full);
    return(NULL);
  }

  /* norm S(q ,w) = sum(S)*q_range/q_bins on total q,w range from file */
  sum *= (q_max_file - q_min_file)/q_bins_file;

  /* (7) renormalization of S(q,w) ========================================== */

  if      ( Sqw->sqw_norm >0) alpha=Sqw->sqw_norm;
  else if (!Sqw->sqw_norm)    alpha=1;

  if (!alpha && iq2Sq) { /* compute theoretical |S| norm */
    /* Eq (2.44) from H.E. Fischer et al, Rep. Prog. Phys., 69 (2006) 233 */
    alpha =
      (q_max*q_max*q_max/3 - (Sqw->type == 'c' ? 2*PI*PI*Sqw->mat_rho : 0))
      /iq2Sq;
  }

  if (alpha < 0) {
    MPI_MASTER(
    printf("Isotropic_Sqw: %s: normalisation factor is negative. rho=%g [Angs^-3] may be too high.\n"
           "WARNING        Disabling renormalization i.e. keeping initial S(q,w).\n",
      Sqw->compname, Sqw->mat_rho);
    );
    alpha=0;
  }

  /* apply normalization on S(q,w) */
  if (alpha && alpha != 1) {
    sum *= alpha;
    for (index_q=0; index_q < q_bins ; index_q++) {
      for (index_w=0; index_w < w_bins; index_w++)
        Table_SetElement(&Sqw_full, index_q, index_w,
          Table_Index(Sqw_full, index_q, index_w) * alpha);
    }
  }

  Sqw_Data->intensity       = sum;

  Table_Stat(&Sqw_full);
  Sqw_full.min_x        = 0;
  Sqw_full.max_x        = q_max;
  Sqw_full.step_x       = q_step;

  MPI_MASTER(
  if (Sqw->verbose_output > 0) {
    printf("Isotropic_Sqw: %s: Generated %s %scoherent Sqw\n"
           "                   q=[%g:%g Angs-1] w=[%g:%g meV] |S|=%g size=[%lix%li] sigma=%g [barns]\n",
           Sqw->compname, file, (Sqw->type == 'i' ? "in" : ""),
           q_min_file, q_max_file,
           w_min_file, w_max_file, Sqw_Data->intensity,
           q_bins, Sqw_Data->w_bins,
           (Sqw->type == 'i' ? Sqw->s_inc : Sqw->s_coh));
    if (w_max < 1e-2)
      printf("               Mainly elastic scattering.\n");
    if (Sqw->sqw_norm >0 && Sqw->sqw_norm != 1)
      printf("                   normalization factor S(q,w)*%g (user)\n", alpha);
    else if (Sqw->sqw_norm<0)
      printf("                   normalization factor S(q,w)*%g (auto) \\int q^2 S(q) dq=%g\n", alpha, iq2Sq);
  }
  );

  /* (8) Compute total cross section ======================================== */

  /* now compute the effective total cross section  (Sqw_integrate_iqSq)
        sigma(Ei) = sigma/2/Ki^2 * \int q S(q,w) dw dq
   * for each incoming neutron energy 0 < Ei < 2*w_max, and
   * integration range w=-Ei:w_max and q=Q0:Q1 with
   *   Q0 = SE2Q*(sqrt(E)-sqrt(E+w))
   *   Q1 = SE2Q*(sqrt(E)+sqrt(E+w))
   */

  Sqw_Data->lookup_length = Sqw->lookup_length;
  Sqw_Data->iqSq_length   = Sqw->lookup_length;
  /* this length should be greater when w_max=0 for e.g. elastic only */
  if (w_bins <= 1) Sqw_Data->iqSq_length = q_bins;

  if (!Table_Init(&iqSq, Sqw_Data->iqSq_length, 1)) {
    /* from 0 to 2*Ki_max */
    printf("Isotropic_Sqw: %s: Cannot allocate [int q S(q,w) dq dw] array (%li bytes).\n"
           "ERROR          Exiting.\n",
      Sqw->compname, Sqw->lookup_length*sizeof(double));
    Table_Free(&Sqw_full);
    return(NULL);
  }

  /* compute the maximum incoming energy that can be handled */
  Sqw_Data->Ei_max = 2*w_max;

  /* Checked in different ways in Powder and "proper inelastic" case */
  if (w_step==1) {
    /* Powder */
    double Ei_max_q = (q_max*K2V)*(q_max*K2V)*VS2E/2;
    if (Ei_max_q > Sqw_Data->Ei_max) Sqw_Data->Ei_max = Ei_max_q;
  } else {
    /* Proper inelastic */
    /* check if the q-range will limit the integration */
    if ((q_max*K2V)*(q_max*K2V)*VS2E/2 > Sqw_Data->Ei_max) {
      /* then scan Ei until we pass q_max */
      for (index_w=0; index_w < Sqw_Data->iqSq_length; index_w++) {
	double Ei = index_w*2*w_max/Sqw_Data->iqSq_length;
	if ( (Ei > w_max && sqrt(Ei)+sqrt(Ei-w_max) >= q_max/(SE2V*V2K))
	     || sqrt(Ei)+sqrt(Ei+w_max) >= q_max/(SE2V*V2K))
	  if (Ei < Sqw_Data->Ei_max) {
	    Sqw_Data->Ei_max = Ei;
	    break;
	  }
      }
    }
  }

  MPI_MASTER(
  if (Sqw->verbose_output >= 2)
    printf("Isotropic_Sqw: %s: Creating Sigma(Ei=0:%g [meV]) with %li entries...(%s %s)\n",
      Sqw->compname, Sqw_Data->Ei_max, Sqw_Data->iqSq_length, file, Sqw->type=='c' ? "coh" : "inc");
  );
  Sqw_Data->Sqw  = Sqw_full; /* store the S(q,w) Table (matrix) for Sqw_integrate_iqSq */

  /* this loop takes time: use MPI when available */

  for (index_w=0; index_w < Sqw_Data->iqSq_length; index_w++) {

    /* Ei = energy of incoming neutron, typically 0:w_max or 0:2*q_max */
    double Ei; /* neutron energy value in Sqw_full, up to 2*w_max */
    double Ki, Vi;
    double Sigma=0;
    Ei = index_w*Sqw_Data->Ei_max/Sqw_Data->iqSq_length;
    Vi = sqrt(Ei/VS2E);
    Ki = V2K*Vi;
    /* sigma(Ei) = sigma/2/Ki^2 * \int q S(q,w) dq dw */
    /* Eq (6) from E. Farhi et al. J. Comp. Phys. 228 (2009) 5251 */
    Sigma = Ki <= 0 ? 0 : (Sqw->type=='c' ? Sqw->s_coh : Sqw->s_inc)
                          /2/Ki/Ki * Sqw_integrate_iqSq(Sqw_Data, Ei);
    Table_SetElement(&iqSq, index_w, 0, Sigma );
  }

  sprintf(iqSq.filename, "[sigma/2Ki^2 int q S(q,w) dq dw] from %s", file);
  iqSq.min_x  = 0;
  iqSq.max_x  = Sqw_Data->Ei_max;
  iqSq.step_x = Sqw_Data->Ei_max/Sqw_Data->iqSq_length;
  iqSq.block_number = 1;

  Sqw_Data->iqSq = iqSq;     /* store the sigma(Ei) = \int q S(q,w) dq dw Table (vector) */

  /* (9) Compute P(w) probability =========================================== */

  /* set up 'density of states' */
  /* uses: Sqw_full and w axes */
  Sqw_Data->SW =
    (struct Sqw_W_struct*)calloc(w_bins, sizeof(struct Sqw_W_struct));

  if (!Sqw_Data->SW) {
    printf("Isotropic_Sqw: %s: Cannot allocate SW (%li bytes).\n"
           "ERROR          Exiting.\n",
      Sqw->compname, w_bins*sizeof(struct Sqw_W_struct));
    Table_Free(&Sqw_full);
    Table_Free(&iqSq);
    return(NULL);
  }
  sum = 0;
  for (index_w=0; index_w < w_bins ; index_w++) {
    double local_val = 0;
    double w         = -w_max + index_w * w_step;
    for (index_q=0; index_q < q_bins ; index_q++) { /* integrate on all q values */
      local_val += Table_Index(Sqw_full, index_q, index_w)*q_step*index_q*q_step; /* q*S(q,w) */
    }
    Sqw_Data->SW[index_w].omega = w;
    sum                  += local_val; /* S(w)=\int S(q,w) dq */
    /* compute cumulated probability */
    Sqw_Data->SW[index_w].cumul_proba = local_val;
    if (index_w) Sqw_Data->SW[index_w].cumul_proba += Sqw_Data->SW[index_w-1].cumul_proba;
    else         Sqw_Data->SW[index_w].cumul_proba = 0;
  }

  if (!sum) {
    MPI_MASTER(
    printf("Isotropic_Sqw: %s: Total S(q,w) intensity is NULL.\n"
           "ERROR          Exiting.\n", Sqw->compname);
    );
    Table_Free(&Sqw_full);
    Table_Free(&iqSq);
    return(NULL);
  }

  /* normalize Pw distribution to 0:1 range */
  for (index_w=0; index_w < w_bins ; index_w++) {
    Sqw_Data->SW[index_w].cumul_proba /= Sqw_Data->SW[w_bins-1].cumul_proba;
  }

  if (Sqw->verbose_output > 2) {
    MPI_MASTER(
    printf("Isotropic_Sqw: %s: Generated normalized SW[%li] in range [0:%g]\n",
      Sqw->compname, w_bins, Sqw_Data->SW[w_bins-1].cumul_proba);
    );
  }

  /* (10) Compute P(Q|w) probability ======================================== */

  /* set up Q probability table per w bin */
  /* uses:  Sqw_full */
  Sqw_Data->SQW =
    (struct Sqw_Q_struct**)calloc(w_bins, sizeof(struct Sqw_Q_struct*));

  if (!Sqw_Data->SQW) {
    printf("Isotropic_Sqw: %s: Cannot allocate P(Q|w) array (%li bytes).\n"
           "ERROR          Exiting.\n",
      Sqw->compname, w_bins*sizeof(struct Sqw_Q_struct*));
    Table_Free(&Sqw_full);
    Table_Free(&iqSq);
    return(NULL);
  }
  for (index_w=0; index_w < w_bins ; index_w++) {
    Sqw_Data->SQW[index_w]=
        (struct Sqw_Q_struct*)calloc(q_bins, sizeof(struct Sqw_Q_struct));

    if (!Sqw_Data->SQW[index_w]) {
      printf("Isotropic_Sqw: %s: Cannot allocate P(Q|w)[%li] (%li bytes).\n"
             "ERROR          Exiting.\n",
        Sqw->compname, index_w, q_bins*sizeof(struct Sqw_Q_struct));
      Table_Free(&Sqw_full);
      Table_Free(&iqSq);
      return(NULL);
    }
    /* set P(Q|W) and compute total intensity */
    for (index_q=0; index_q < q_bins ; index_q++) {
      double q  = index_q * q_step;
      Sqw_Data->SQW[index_w][index_q].Q     = q;

      /* compute cumulated probability and take into account Jacobian with additional 'q' factor */
      Sqw_Data->SQW[index_w][index_q].cumul_proba = q*Table_Index(Sqw_full, index_q, index_w); /* q*S(q,w) */
      if (index_q) Sqw_Data->SQW[index_w][index_q].cumul_proba += Sqw_Data->SQW[index_w][index_q-1].cumul_proba;
      else Sqw_Data->SQW[index_w][index_q].cumul_proba = 0;
    }
    /* normalize P(q|w) distribution to 0:1 range */
    for (index_q=0; index_q < q_bins ;
    	Sqw_Data->SQW[index_w][index_q++].cumul_proba /= Sqw_Data->SQW[index_w][q_bins-1].cumul_proba
    );

  }
  if (Sqw->verbose_output > 2) {
    MPI_MASTER(
    printf("Isotropic_Sqw: %s: Generated P(Q|w)\n",
      Sqw->compname);
    );
  }

  /* (11) generate quick lookup tables for SW and SQW ======================= */

  SW_lookup = (long*)calloc(Sqw->lookup_length, sizeof(long));

  if (!SW_lookup) {
    printf("Isotropic_Sqw: %s: Cannot allocate SW_lookup (%li bytes).\n"
           "Warning        Will be slower.\n",
      Sqw->compname, Sqw->lookup_length*sizeof(long));
  } else {
    int i;
    for (i=0; i < Sqw->lookup_length; i++) {
      double w = (double)i/(double)Sqw->lookup_length; /* a random number tabulated value */
      SW_lookup[i] = Sqw_search_SW(*Sqw_Data, w);
    }
    Sqw_Data->SW_lookup = SW_lookup;
  }
  QW_lookup = (long**)calloc(w_bins, sizeof(long*));

  if (!QW_lookup) {
    printf("Isotropic_Sqw: %s: Cannot allocate QW_lookup (%li bytes).\n"
           "Warning        Will be slower.\n",
      Sqw->compname, w_bins*sizeof(long*));
  } else {
    for (index_w=0; index_w < w_bins ; index_w++) {
      QW_lookup[index_w] =
        (long*)calloc(Sqw->lookup_length, sizeof(long));
      if (!QW_lookup[index_w]) {
        printf("Isotropic_Sqw: %s: Cannot allocate QW_lookup[%li] (%li bytes).\n"
               "Warning        Will be slower.\n",
        Sqw->compname, index_w, Sqw->lookup_length*sizeof(long));
        free(QW_lookup); Sqw_Data->QW_lookup = QW_lookup = NULL; break;
      } else {
        int i;
        for (i=0; i < Sqw->lookup_length; i++) {
          double w = (double)i/(double)Sqw->lookup_length; /* a random number tabulated value */
          QW_lookup[index_w][i] = Sqw_search_Q_proba_per_w(*Sqw_Data, w, index_w);
        }
      }
    }
    Sqw_Data->QW_lookup = QW_lookup;
  }
  if ((Sqw_Data->QW_lookup || Sqw_Data->SW_lookup) && Sqw->verbose_output > 2) {
    MPI_MASTER(
    printf("Isotropic_Sqw: %s: Generated lookup tables with %li entries\n",
      Sqw->compname, Sqw->lookup_length);
    );
  }

  return(Sqw_Data);
} /* end Sqw_readfile */

/*****************************************************************************
* Sqw_init_read: Read coherent/incoherent Sqw data files
*   Returns Sqw total intensity, or 0 (error)
* Used in : INITIALIZE (1)
*****************************************************************************/
double Sqw_init(struct Sqw_sample_struct *Sqw, char *file_coh, char *file_inc)
{
  double ret=0;

  /* read files */
  struct Sqw_Data_struct *d_inc, *d_coh;
  Sqw->type = 'i';
  d_inc = Sqw_readfile(Sqw, file_inc, &(Sqw->Data_inc));
  Sqw->type = 'c';
  d_coh = Sqw_readfile(Sqw, file_coh, &(Sqw->Data_coh));

  if (d_inc && !d_inc->intensity && Sqw->s_inc>0) {
    MPI_MASTER(
    if (Sqw->verbose_output > 0)
      printf("Isotropic_Sqw: %s: Using Isotropic elastic incoherent scattering (sigma=%g [barns])\n", Sqw->compname, Sqw->s_inc);
    );
    ret=1;
  }

  if (!d_inc || !d_coh) return(0);

  d_coh->type = 'c';
  Sqw->Data_inc.type = 'i';
  MPI_MASTER(
  if (d_coh && !d_coh->intensity && Sqw->s_coh)
    printf("Isotropic_Sqw: %s: Coherent scattering Sqw intensity is null.\n"
           "Warning        Disabling coherent scattering.\n", Sqw->compname);
  );
  if (d_inc && d_coh && d_inc->intensity && d_coh->intensity) {
    char msg[80];
    strcpy(msg, "");
    /* check dimensions/limits for Q, Energy in coh and inc Tables */
    if (d_inc->q_bins  != d_coh->q_bins)
      strcpy(msg, "Q axis size");
    if (d_inc->w_bins  != d_coh->w_bins)
      strcpy(msg, "Energy axis size");
    if (d_inc->q_max != d_coh->q_max)
      strcpy(msg, "Q axis limits");
    if (d_inc->w_max != d_coh->w_max)
      strcpy(msg, "Energy axis limits");
    MPI_MASTER(
    if (strlen(msg)) {
      printf("Isotropic_Sqw: %s: Sqw data from files %s and %s do not match\n"
             "WARNING        wrong %s\n",
             Sqw->compname, file_coh, file_inc, msg);
    }
    );
  }

  if (!ret) ret=d_inc->intensity+d_coh->intensity;
  return(ret);
} /* Sqw_init */

#endif /* definied ISOTROPIC_SQW */

/* Shared user declarations for all components types 'PSD_Detector'. */



#ifndef PSD_Detector_SHARE
#define PSD_Detector_SHARE
double PSD_He_interp1value(double *array_x,double *array_y,long n_elements,double interp_x) {
    long cnt1 = 0;
    double a;
    // while (interp_x > array_x[++cnt1]);
    while (cnt1 < n_elements-1 && interp_x >= array_x[++cnt1]);

    a = (interp_x - array_x[cnt1 - 1]) / (array_x[cnt1] - array_x[cnt1 - 1]);
    return array_y[cnt1 - 1] + a * (array_y[cnt1] - array_y[cnt1 - 1]);
  }
  #endif


/* ************************************************************************** */
/*             End of SHARE user declarations for all components              */
/* ************************************************************************** */


/* ********************** component definition declarations. **************** */

/* component PG=Progress_bar() [1] DECLARE */
/* Parameter definition for component type 'Progress_bar' */
struct _struct_Progress_bar_parameters {
  /* Component type 'Progress_bar' setting parameters */
  char profile[16384];
  MCNUM percent;
  MCNUM flag_save;
  MCNUM minutes;
  /* Component type 'Progress_bar' private parameters */
  /* Component type 'Progress_bar' DECLARE code stored as structure members */
#ifndef PROGRESS_BAR
#define PROGRESS_BAR
#else
#error Only one Progress_bar component may be used in an instrument definition.
#endif

double IntermediateCnts;
time_t StartTime;
time_t EndTime;
time_t CurrentTime;
}; /* _struct_Progress_bar_parameters */
typedef struct _struct_Progress_bar_parameters _class_Progress_bar_parameters;

/* Parameters for component type 'Progress_bar' */
struct _struct_Progress_bar {
  char     _name[256]; /* e.g. PG */
  char     _type[256]; /* Progress_bar */
  long     _index; /* e.g. 2 index in TRACE list */
  Coords   _position_absolute;
  Coords   _position_relative; /* wrt PREVIOUS */
  Rotation _rotation_absolute;
  Rotation _rotation_relative; /* wrt PREVIOUS */
  int      _rotation_is_identity;
  _class_Progress_bar_parameters _parameters;
};
typedef struct _struct_Progress_bar _class_Progress_bar;
_class_Progress_bar _PG_var;
_class_Progress_bar *_PG = &_PG_var;
#pragma acc declare create ( _PG_var )
#pragma acc declare create ( _PG )

/* component Source=Source_gen() [2] DECLARE */
/* Parameter definition for component type 'Source_gen' */
struct _struct_Source_gen_parameters {
  /* Component type 'Source_gen' setting parameters */
  char flux_file[16384];
  char xdiv_file[16384];
  char ydiv_file[16384];
  MCNUM radius;
  MCNUM dist;
  MCNUM focus_xw;
  MCNUM focus_yh;
  MCNUM focus_aw;
  MCNUM focus_ah;
  MCNUM E0;
  MCNUM dE;
  MCNUM lambda0;
  MCNUM dlambda;
  MCNUM I1;
  MCNUM yheight;
  MCNUM xwidth;
  MCNUM verbose;
  MCNUM T1;
  MCNUM flux_file_perAA;
  MCNUM flux_file_log;
  MCNUM Lmin;
  MCNUM Lmax;
  MCNUM Emin;
  MCNUM Emax;
  MCNUM T2;
  MCNUM I2;
  MCNUM T3;
  MCNUM I3;
  MCNUM zdepth;
  long target_index;
  /* Component type 'Source_gen' private parameters */
  /* Component type 'Source_gen' DECLARE code stored as structure members */

  double p_in;
  double lambda1;                               
  double lambda2;                                
  double lambda3;                               
  t_Table pTable;
  t_Table pTable_x;
  t_Table pTable_y;
  double pTable_xmin;
  double pTable_xmax;
  double pTable_xsum;
  double pTable_ymin;
  double pTable_ymax;
  double pTable_ysum;
  double pTable_dxmin;
  double pTable_dxmax;
  double pTable_dymin;
  double pTable_dymax;

}; /* _struct_Source_gen_parameters */
typedef struct _struct_Source_gen_parameters _class_Source_gen_parameters;

/* Parameters for component type 'Source_gen' */
struct _struct_Source_gen {
  char     _name[256]; /* e.g. Source */
  char     _type[256]; /* Source_gen */
  long     _index; /* e.g. 2 index in TRACE list */
  Coords   _position_absolute;
  Coords   _position_relative; /* wrt PREVIOUS */
  Rotation _rotation_absolute;
  Rotation _rotation_relative; /* wrt PREVIOUS */
  int      _rotation_is_identity;
  _class_Source_gen_parameters _parameters;
};
typedef struct _struct_Source_gen _class_Source_gen;
_class_Source_gen _Source_var;
_class_Source_gen *_Source = &_Source_var;
#pragma acc declare create ( _Source_var )
#pragma acc declare create ( _Source )

/* component Guide0=Guide_gravity() [3] DECLARE */
/* Parameter definition for component type 'Guide_gravity' */
struct _struct_Guide_gravity_parameters {
  /* Component type 'Guide_gravity' setting parameters */
  MCNUM w1;
  MCNUM h1;
  MCNUM w2;
  MCNUM h2;
  MCNUM l;
  MCNUM R0;
  MCNUM Qc;
  MCNUM alpha;
  MCNUM m;
  MCNUM W;
  MCNUM nslit;
  MCNUM d;
  MCNUM mleft;
  MCNUM mright;
  MCNUM mtop;
  MCNUM mbottom;
  MCNUM nhslit;
  MCNUM G;
  MCNUM aleft;
  MCNUM aright;
  MCNUM atop;
  MCNUM abottom;
  MCNUM wavy;
  MCNUM wavy_z;
  MCNUM wavy_tb;
  MCNUM wavy_lr;
  MCNUM chamfers;
  MCNUM chamfers_z;
  MCNUM chamfers_lr;
  MCNUM chamfers_tb;
  MCNUM nelements;
  MCNUM nu;
  MCNUM phase;
  char reflect[16384];
  /* Component type 'Guide_gravity' private parameters */
  /* Component type 'Guide_gravity' DECLARE code stored as structure members */
  Gravity_guide_Vars_type GVars;
  t_Table pTable;
}; /* _struct_Guide_gravity_parameters */
typedef struct _struct_Guide_gravity_parameters _class_Guide_gravity_parameters;

/* Parameters for component type 'Guide_gravity' */
struct _struct_Guide_gravity {
  char     _name[256]; /* e.g. Guide0 */
  char     _type[256]; /* Guide_gravity */
  long     _index; /* e.g. 2 index in TRACE list */
  Coords   _position_absolute;
  Coords   _position_relative; /* wrt PREVIOUS */
  Rotation _rotation_absolute;
  Rotation _rotation_relative; /* wrt PREVIOUS */
  int      _rotation_is_identity;
  _class_Guide_gravity_parameters _parameters;
};
typedef struct _struct_Guide_gravity _class_Guide_gravity;
_class_Guide_gravity _Guide0_var;
_class_Guide_gravity *_Guide0 = &_Guide0_var;
#pragma acc declare create ( _Guide0_var )
#pragma acc declare create ( _Guide0 )

_class_Guide_gravity _Guide0_4_var;
_class_Guide_gravity *_Guide0_4 = &_Guide0_4_var;
#pragma acc declare create ( _Guide0_4_var )
#pragma acc declare create ( _Guide0_4 )

_class_Guide_gravity _Guide0_5_var;
_class_Guide_gravity *_Guide0_5 = &_Guide0_5_var;
#pragma acc declare create ( _Guide0_5_var )
#pragma acc declare create ( _Guide0_5 )

_class_Guide_gravity _Guide0_6_var;
_class_Guide_gravity *_Guide0_6 = &_Guide0_6_var;
#pragma acc declare create ( _Guide0_6_var )
#pragma acc declare create ( _Guide0_6 )

_class_Guide_gravity _Guide0_7_var;
_class_Guide_gravity *_Guide0_7 = &_Guide0_7_var;
#pragma acc declare create ( _Guide0_7_var )
#pragma acc declare create ( _Guide0_7 )

_class_Guide_gravity _Guide0_8_var;
_class_Guide_gravity *_Guide0_8 = &_Guide0_8_var;
#pragma acc declare create ( _Guide0_8_var )
#pragma acc declare create ( _Guide0_8 )

_class_Guide_gravity _Guide0_9_var;
_class_Guide_gravity *_Guide0_9 = &_Guide0_9_var;
#pragma acc declare create ( _Guide0_9_var )
#pragma acc declare create ( _Guide0_9 )

_class_Guide_gravity _Guide0_10_var;
_class_Guide_gravity *_Guide0_10 = &_Guide0_10_var;
#pragma acc declare create ( _Guide0_10_var )
#pragma acc declare create ( _Guide0_10 )

_class_Guide_gravity _Guide0_11_var;
_class_Guide_gravity *_Guide0_11 = &_Guide0_11_var;
#pragma acc declare create ( _Guide0_11_var )
#pragma acc declare create ( _Guide0_11 )

_class_Guide_gravity _Guide0_12_var;
_class_Guide_gravity *_Guide0_12 = &_Guide0_12_var;
#pragma acc declare create ( _Guide0_12_var )
#pragma acc declare create ( _Guide0_12 )

_class_Guide_gravity _Guide0_13_var;
_class_Guide_gravity *_Guide0_13 = &_Guide0_13_var;
#pragma acc declare create ( _Guide0_13_var )
#pragma acc declare create ( _Guide0_13 )

_class_Guide_gravity _Guide0_14_var;
_class_Guide_gravity *_Guide0_14 = &_Guide0_14_var;
#pragma acc declare create ( _Guide0_14_var )
#pragma acc declare create ( _Guide0_14 )

_class_Guide_gravity _Guide0_15_var;
_class_Guide_gravity *_Guide0_15 = &_Guide0_15_var;
#pragma acc declare create ( _Guide0_15_var )
#pragma acc declare create ( _Guide0_15 )

_class_Guide_gravity _Guide0_16_var;
_class_Guide_gravity *_Guide0_16 = &_Guide0_16_var;
#pragma acc declare create ( _Guide0_16_var )
#pragma acc declare create ( _Guide0_16 )

_class_Guide_gravity _Guide0_17_var;
_class_Guide_gravity *_Guide0_17 = &_Guide0_17_var;
#pragma acc declare create ( _Guide0_17_var )
#pragma acc declare create ( _Guide0_17 )

_class_Guide_gravity _Guide0_18_var;
_class_Guide_gravity *_Guide0_18 = &_Guide0_18_var;
#pragma acc declare create ( _Guide0_18_var )
#pragma acc declare create ( _Guide0_18 )

_class_Guide_gravity _Guide0_19_var;
_class_Guide_gravity *_Guide0_19 = &_Guide0_19_var;
#pragma acc declare create ( _Guide0_19_var )
#pragma acc declare create ( _Guide0_19 )

/* component Guide_PSD=Monitor_nD() [20] DECLARE */
/* Parameter definition for component type 'Monitor_nD' */
struct _struct_Monitor_nD_parameters {
  /* Component type 'Monitor_nD' setting parameters */
  double user1;
  double user2;
  double user3;
  MCNUM xwidth;
  MCNUM yheight;
  MCNUM zdepth;
  MCNUM xmin;
  MCNUM xmax;
  MCNUM ymin;
  MCNUM ymax;
  MCNUM zmin;
  MCNUM zmax;
  MCNUM bins;
  MCNUM min;
  MCNUM max;
  MCNUM restore_neutron;
  MCNUM radius;
  char options[16384];
  char filename[16384];
  char geometry[16384];
  char username1[16384];
  char username2[16384];
  char username3[16384];
  /* Component type 'Monitor_nD' private parameters */
  /* Component type 'Monitor_nD' DECLARE code stored as structure members */
  MonitornD_Defines_type DEFS;
  MonitornD_Variables_type Vars;
  MCDETECTOR detector;
  off_struct offdata;
}; /* _struct_Monitor_nD_parameters */
typedef struct _struct_Monitor_nD_parameters _class_Monitor_nD_parameters;

/* Parameters for component type 'Monitor_nD' */
struct _struct_Monitor_nD {
  char     _name[256]; /* e.g. Guide_PSD */
  char     _type[256]; /* Monitor_nD */
  long     _index; /* e.g. 2 index in TRACE list */
  Coords   _position_absolute;
  Coords   _position_relative; /* wrt PREVIOUS */
  Rotation _rotation_absolute;
  Rotation _rotation_relative; /* wrt PREVIOUS */
  int      _rotation_is_identity;
  _class_Monitor_nD_parameters _parameters;
};
typedef struct _struct_Monitor_nD _class_Monitor_nD;
_class_Monitor_nD _Guide_PSD_var;
_class_Monitor_nD *_Guide_PSD = &_Guide_PSD_var;
#pragma acc declare create ( _Guide_PSD_var )
#pragma acc declare create ( _Guide_PSD )

_class_Monitor_nD _Guide_Lambda_var;
_class_Monitor_nD *_Guide_Lambda = &_Guide_Lambda_var;
#pragma acc declare create ( _Guide_Lambda_var )
#pragma acc declare create ( _Guide_Lambda )

/* component Chopper1=DiskChopper() [22] DECLARE */
/* Parameter definition for component type 'DiskChopper' */
struct _struct_DiskChopper_parameters {
  /* Component type 'DiskChopper' setting parameters */
  MCNUM theta_0;
  MCNUM radius;
  MCNUM yheight;
  MCNUM nu;
  MCNUM nslit;
  MCNUM jitter;
  MCNUM delay;
  MCNUM isfirst;
  MCNUM n_pulse;
  MCNUM abs_out;
  MCNUM phase;
  MCNUM xwidth;
  MCNUM verbose;
  /* Component type 'DiskChopper' private parameters */
  /* Component type 'DiskChopper' DECLARE code stored as structure members */
double Tg,To,delta_y,height,omega;
}; /* _struct_DiskChopper_parameters */
typedef struct _struct_DiskChopper_parameters _class_DiskChopper_parameters;

/* Parameters for component type 'DiskChopper' */
struct _struct_DiskChopper {
  char     _name[256]; /* e.g. Chopper1 */
  char     _type[256]; /* DiskChopper */
  long     _index; /* e.g. 2 index in TRACE list */
  Coords   _position_absolute;
  Coords   _position_relative; /* wrt PREVIOUS */
  Rotation _rotation_absolute;
  Rotation _rotation_relative; /* wrt PREVIOUS */
  int      _rotation_is_identity;
  _class_DiskChopper_parameters _parameters;
};
typedef struct _struct_DiskChopper _class_DiskChopper;
_class_DiskChopper _Chopper1_var;
_class_DiskChopper *_Chopper1 = &_Chopper1_var;
#pragma acc declare create ( _Chopper1_var )
#pragma acc declare create ( _Chopper1 )

_class_Guide_gravity _Guide_SGS1_var;
_class_Guide_gravity *_Guide_SGS1 = &_Guide_SGS1_var;
#pragma acc declare create ( _Guide_SGS1_var )
#pragma acc declare create ( _Guide_SGS1 )

_class_Guide_gravity _Guide_SGS2_var;
_class_Guide_gravity *_Guide_SGS2 = &_Guide_SGS2_var;
#pragma acc declare create ( _Guide_SGS2_var )
#pragma acc declare create ( _Guide_SGS2 )

_class_Guide_gravity _Guide_SGS3_var;
_class_Guide_gravity *_Guide_SGS3 = &_Guide_SGS3_var;
#pragma acc declare create ( _Guide_SGS3_var )
#pragma acc declare create ( _Guide_SGS3 )

_class_Guide_gravity _Guide_CGS_var;
_class_Guide_gravity *_Guide_CGS = &_Guide_CGS_var;
#pragma acc declare create ( _Guide_CGS_var )
#pragma acc declare create ( _Guide_CGS )

_class_Monitor_nD _Guide2_PSD_var;
_class_Monitor_nD *_Guide2_PSD = &_Guide2_PSD_var;
#pragma acc declare create ( _Guide2_PSD_var )
#pragma acc declare create ( _Guide2_PSD )

_class_Monitor_nD _Guide2_Lambda_var;
_class_Monitor_nD *_Guide2_Lambda = &_Guide2_Lambda_var;
#pragma acc declare create ( _Guide2_Lambda_var )
#pragma acc declare create ( _Guide2_Lambda )

_class_Monitor_nD _Guide2_dXY_var;
_class_Monitor_nD *_Guide2_dXY = &_Guide2_dXY_var;
#pragma acc declare create ( _Guide2_dXY_var )
#pragma acc declare create ( _Guide2_dXY )

_class_DiskChopper _Chopper6_var;
_class_DiskChopper *_Chopper6 = &_Chopper6_var;
#pragma acc declare create ( _Chopper6_var )
#pragma acc declare create ( _Chopper6 )

_class_Monitor_nD _Guide_time_var;
_class_Monitor_nD *_Guide_time = &_Guide_time_var;
#pragma acc declare create ( _Guide_time_var )
#pragma acc declare create ( _Guide_time )

_class_Guide_gravity _Guide_DGS_var;
_class_Guide_gravity *_Guide_DGS = &_Guide_DGS_var;
#pragma acc declare create ( _Guide_DGS_var )
#pragma acc declare create ( _Guide_DGS )

/* component Sample_pos=Arm() [33] DECLARE */
/* Parameter definition for component type 'Arm' */
struct _struct_Arm_parameters {
  char Arm_has_no_parameters;
}; /* _struct_Arm_parameters */
typedef struct _struct_Arm_parameters _class_Arm_parameters;

/* Parameters for component type 'Arm' */
struct _struct_Arm {
  char     _name[256]; /* e.g. Sample_pos */
  char     _type[256]; /* Arm */
  long     _index; /* e.g. 2 index in TRACE list */
  Coords   _position_absolute;
  Coords   _position_relative; /* wrt PREVIOUS */
  Rotation _rotation_absolute;
  Rotation _rotation_relative; /* wrt PREVIOUS */
  int      _rotation_is_identity;
  _class_Arm_parameters _parameters;
};
typedef struct _struct_Arm _class_Arm;
_class_Arm _Sample_pos_var;
_class_Arm *_Sample_pos = &_Sample_pos_var;
#pragma acc declare create ( _Sample_pos_var )
#pragma acc declare create ( _Sample_pos )

/* component Sample=Isotropic_Sqw() [34] DECLARE */
/* Parameter definition for component type 'Isotropic_Sqw' */
struct _struct_Isotropic_Sqw_parameters {
  /* Component type 'Isotropic_Sqw' setting parameters */
  MCNUM* powder_format;
  char Sqw_coh[16384];
  char Sqw_inc[16384];
  char geometry[16384];
  MCNUM radius;
  MCNUM thickness;
  MCNUM xwidth;
  MCNUM yheight;
  MCNUM zdepth;
  MCNUM threshold;
  long order;
  MCNUM T;
  MCNUM verbose;
  MCNUM d_phi;
  long concentric;
  MCNUM rho;
  MCNUM sigma_abs;
  MCNUM sigma_coh;
  MCNUM sigma_inc;
  MCNUM classical;
  MCNUM powder_Dd;
  MCNUM powder_DW;
  MCNUM powder_Vc;
  MCNUM density;
  MCNUM weight;
  MCNUM p_interact;
  MCNUM norm;
  MCNUM powder_barns;
  char quantum_correction[16384];
  /* Component type 'Isotropic_Sqw' private parameters */
  /* Component type 'Isotropic_Sqw' DECLARE code stored as structure members */
  struct Sqw_sample_struct VarSqw;
  int *columns;
  off_struct offdata;
}; /* _struct_Isotropic_Sqw_parameters */
typedef struct _struct_Isotropic_Sqw_parameters _class_Isotropic_Sqw_parameters;

/* Parameters for component type 'Isotropic_Sqw' */
struct _struct_Isotropic_Sqw {
  char     _name[256]; /* e.g. Sample */
  char     _type[256]; /* Isotropic_Sqw */
  long     _index; /* e.g. 2 index in TRACE list */
  Coords   _position_absolute;
  Coords   _position_relative; /* wrt PREVIOUS */
  Rotation _rotation_absolute;
  Rotation _rotation_relative; /* wrt PREVIOUS */
  int      _rotation_is_identity;
  _class_Isotropic_Sqw_parameters _parameters;
};
typedef struct _struct_Isotropic_Sqw _class_Isotropic_Sqw;
_class_Isotropic_Sqw _Sample_var;
_class_Isotropic_Sqw *_Sample = &_Sample_var;
#pragma acc declare create ( _Sample_var )
#pragma acc declare create ( _Sample )

_class_Monitor_nD _Ideal_Det_var;
_class_Monitor_nD *_Ideal_Det = &_Ideal_Det_var;
#pragma acc declare create ( _Ideal_Det_var )
#pragma acc declare create ( _Ideal_Det )

/* component Detector=PSD_Detector() [36] DECLARE */
/* Parameter definition for component type 'PSD_Detector' */
struct _struct_PSD_Detector_parameters {
  /* Component type 'PSD_Detector' setting parameters */
  MCNUM nx;
  MCNUM ny;
  MCNUM xwidth;
  MCNUM radius;
  MCNUM awidth;
  MCNUM yheight;
  MCNUM zdepth;
  MCNUM threshold;
  MCNUM PressureConv;
  MCNUM PressureStop;
  MCNUM interpolate;
  MCNUM p_interact;
  MCNUM verbose;
  MCNUM LensOn;
  MCNUM dc;
  MCNUM borderx;
  MCNUM bordery;
  MCNUM xChDivRelSigma;
  MCNUM yChDivRelSigma;
  MCNUM bufsize;
  MCNUM restore_neutron;
  MCNUM angle;
  char type[16384];
  char filename[16384];
  char FN_Conv[16384];
  char FN_Stop[16384];
  /* Component type 'PSD_Detector' private parameters */
  /* Component type 'PSD_Detector' DECLARE code stored as structure members */
  MonitornD_Defines_type   DEFS;                       
  MonitornD_Variables_type Vars;

  DArray2d PSD_N;
  DArray2d PSD_p;
  DArray2d PSD_p2;

  double *EAP;                                                                          
  double *EAT;                                                                          
  double *M1P1;                                                      
  double *M1T1;                                                      
  double *PosAP;                                                          
  double *PosAT;                                                          
  double *PHSpectrum0;                      
  double *PHSpectrum;                                                                  
  double *PHSpectrum2;                     
  long PHSpectrum_n;
  double CrossSectionHe;                                                       

  double CountNeutrons;                                                   
  double GeomCumul;                                           
                                                                          
  double AbsCumul;                               
                                                                       
                              
  double SensVolCumul;                               
                                                                       
                                                                
                                        
  double DetCumul;                               
                                                                     
                                           
  long nH_p;                                                       
  long nH_t;                                                       
  double FullEnergyP;
  double FullEnergyT;
  long VariousErrors;                                           
                                                                 
                                     
  long DetectorType;                                                                                   
}; /* _struct_PSD_Detector_parameters */
typedef struct _struct_PSD_Detector_parameters _class_PSD_Detector_parameters;

/* Parameters for component type 'PSD_Detector' */
struct _struct_PSD_Detector {
  char     _name[256]; /* e.g. Detector */
  char     _type[256]; /* PSD_Detector */
  long     _index; /* e.g. 2 index in TRACE list */
  Coords   _position_absolute;
  Coords   _position_relative; /* wrt PREVIOUS */
  Rotation _rotation_absolute;
  Rotation _rotation_relative; /* wrt PREVIOUS */
  int      _rotation_is_identity;
  _class_PSD_Detector_parameters _parameters;
};
typedef struct _struct_PSD_Detector _class_PSD_Detector;
_class_PSD_Detector _Detector_var;
_class_PSD_Detector *_Detector = &_Detector_var;
#pragma acc declare create ( _Detector_var )
#pragma acc declare create ( _Detector )

int mcNUMCOMP = 36;

/* User declarations from instrument definition. Can define functions. */
  double time_elastic;
  char   detector_options[1024];

#undef compcurname
#undef compcurtype
#undef compcurindex
/* end of instrument 'HZB_NEAT' and components DECLARE */

/* *****************************************************************************
* instrument 'HZB_NEAT' and components INITIALISE
***************************************************************************** */

/* component PG=Progress_bar() SETTING, POSITION/ROTATION */
int _PG_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_PG_setpos] component PG=Progress_bar() SETTING [/usr/share/mcstas/3.0-dev/misc/Progress_bar.comp:57]");
  stracpy(_PG->_name, "PG", 16384);
  stracpy(_PG->_type, "Progress_bar", 16384);
  _PG->_index=1;
  if("NULL" && strlen("NULL"))
    stracpy(_PG->_parameters.profile, "NULL" ? "NULL" : "", 16384);
  else 
  _PG->_parameters.profile[0]='\0';
  #define profile (_PG->_parameters.profile)
  _PG->_parameters.percent = 10;
  #define percent (_PG->_parameters.percent)
  _PG->_parameters.flag_save = 0;
  #define flag_save (_PG->_parameters.flag_save)
  _PG->_parameters.minutes = 0;
  #define minutes (_PG->_parameters.minutes)

  #define IntermediateCnts (_PG->_parameters.IntermediateCnts)
  #define StartTime (_PG->_parameters.StartTime)
  #define EndTime (_PG->_parameters.EndTime)
  #define CurrentTime (_PG->_parameters.CurrentTime)

  #undef profile
  #undef percent
  #undef flag_save
  #undef minutes
  #undef IntermediateCnts
  #undef StartTime
  #undef EndTime
  #undef CurrentTime
  /* component PG=Progress_bar() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(_PG->_rotation_absolute,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_copy(_PG->_rotation_relative, _PG->_rotation_absolute);
    _PG->_rotation_is_identity =  rot_test_identity(_PG->_rotation_relative);
    _PG->_position_absolute = coords_set(
      0, 0, 0);
    tc1 = coords_neg(_PG->_position_absolute);
    _PG->_position_relative = rot_apply(_PG->_rotation_absolute, tc1);
  } /* PG=Progress_bar() AT ROTATED */
  DEBUG_COMPONENT("PG", _PG->_position_absolute, _PG->_rotation_absolute);
  instrument->_position_absolute[1] = _PG->_position_absolute;
  instrument->_position_relative[1] = _PG->_position_relative;
  instrument->counter_N[1]  = instrument->counter_P[1] = instrument->counter_P2[1] = 0;
  instrument->counter_AbsorbProp[1]= 0;
  return(0);
} /* _PG_setpos */

/* component Source=Source_gen() SETTING, POSITION/ROTATION */
int _Source_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Source_setpos] component Source=Source_gen() SETTING [/usr/share/mcstas/3.0-dev/sources/Source_gen.comp:206]");
  stracpy(_Source->_name, "Source", 16384);
  stracpy(_Source->_type, "Source_gen", 16384);
  _Source->_index=2;
  if("NULL" && strlen("NULL"))
    stracpy(_Source->_parameters.flux_file, "NULL" ? "NULL" : "", 16384);
  else 
  _Source->_parameters.flux_file[0]='\0';
  #define flux_file (_Source->_parameters.flux_file)
  if("NULL" && strlen("NULL"))
    stracpy(_Source->_parameters.xdiv_file, "NULL" ? "NULL" : "", 16384);
  else 
  _Source->_parameters.xdiv_file[0]='\0';
  #define xdiv_file (_Source->_parameters.xdiv_file)
  if("NULL" && strlen("NULL"))
    stracpy(_Source->_parameters.ydiv_file, "NULL" ? "NULL" : "", 16384);
  else 
  _Source->_parameters.ydiv_file[0]='\0';
  #define ydiv_file (_Source->_parameters.ydiv_file)
  _Source->_parameters.radius = .155;
  #define radius (_Source->_parameters.radius)
  _Source->_parameters.dist = 0;
  #define dist (_Source->_parameters.dist)
  _Source->_parameters.focus_xw = 0.03;
  #define focus_xw (_Source->_parameters.focus_xw)
  _Source->_parameters.focus_yh = 0.055;
  #define focus_yh (_Source->_parameters.focus_yh)
  _Source->_parameters.focus_aw = 0;
  #define focus_aw (_Source->_parameters.focus_aw)
  _Source->_parameters.focus_ah = 0;
  #define focus_ah (_Source->_parameters.focus_ah)
  _Source->_parameters.E0 = 0;
  #define E0 (_Source->_parameters.E0)
  _Source->_parameters.dE = 0;
  #define dE (_Source->_parameters.dE)
  _Source->_parameters.lambda0 = instrument->_parameters._lambda;
  #define lambda0 (_Source->_parameters.lambda0)
  _Source->_parameters.dlambda = instrument->_parameters._dlambda;
  #define dlambda (_Source->_parameters.dlambda)
  _Source->_parameters.I1 = 1.4e12;
  #define I1 (_Source->_parameters.I1)
  _Source->_parameters.yheight = 0.1;
  #define yheight (_Source->_parameters.yheight)
  _Source->_parameters.xwidth = 0.1;
  #define xwidth (_Source->_parameters.xwidth)
  _Source->_parameters.verbose = 0;
  #define verbose (_Source->_parameters.verbose)
  _Source->_parameters.T1 = 43.7;
  #define T1 (_Source->_parameters.T1)
  _Source->_parameters.flux_file_perAA = 0;
  #define flux_file_perAA (_Source->_parameters.flux_file_perAA)
  _Source->_parameters.flux_file_log = 0;
  #define flux_file_log (_Source->_parameters.flux_file_log)
  _Source->_parameters.Lmin = 0;
  #define Lmin (_Source->_parameters.Lmin)
  _Source->_parameters.Lmax = 0;
  #define Lmax (_Source->_parameters.Lmax)
  _Source->_parameters.Emin = 0;
  #define Emin (_Source->_parameters.Emin)
  _Source->_parameters.Emax = 0;
  #define Emax (_Source->_parameters.Emax)
  _Source->_parameters.T2 = 137.2;
  #define T2 (_Source->_parameters.T2)
  _Source->_parameters.I2 = 2.08e12;
  #define I2 (_Source->_parameters.I2)
  _Source->_parameters.T3 = 0;
  #define T3 (_Source->_parameters.T3)
  _Source->_parameters.I3 = 0;
  #define I3 (_Source->_parameters.I3)
  _Source->_parameters.zdepth = 0;
  #define zdepth (_Source->_parameters.zdepth)
  _Source->_parameters.target_index = + 1;
  #define target_index (_Source->_parameters.target_index)

  #define p_in (_Source->_parameters.p_in)
  #define lambda1 (_Source->_parameters.lambda1)
  #define lambda2 (_Source->_parameters.lambda2)
  #define lambda3 (_Source->_parameters.lambda3)
  #define pTable (_Source->_parameters.pTable)
  #define pTable_x (_Source->_parameters.pTable_x)
  #define pTable_y (_Source->_parameters.pTable_y)
  #define pTable_xmin (_Source->_parameters.pTable_xmin)
  #define pTable_xmax (_Source->_parameters.pTable_xmax)
  #define pTable_xsum (_Source->_parameters.pTable_xsum)
  #define pTable_ymin (_Source->_parameters.pTable_ymin)
  #define pTable_ymax (_Source->_parameters.pTable_ymax)
  #define pTable_ysum (_Source->_parameters.pTable_ysum)
  #define pTable_dxmin (_Source->_parameters.pTable_dxmin)
  #define pTable_dxmax (_Source->_parameters.pTable_dxmax)
  #define pTable_dymin (_Source->_parameters.pTable_dymin)
  #define pTable_dymax (_Source->_parameters.pTable_dymax)

  #undef flux_file
  #undef xdiv_file
  #undef ydiv_file
  #undef radius
  #undef dist
  #undef focus_xw
  #undef focus_yh
  #undef focus_aw
  #undef focus_ah
  #undef E0
  #undef dE
  #undef lambda0
  #undef dlambda
  #undef I1
  #undef yheight
  #undef xwidth
  #undef verbose
  #undef T1
  #undef flux_file_perAA
  #undef flux_file_log
  #undef Lmin
  #undef Lmax
  #undef Emin
  #undef Emax
  #undef T2
  #undef I2
  #undef T3
  #undef I3
  #undef zdepth
  #undef target_index
  #undef p_in
  #undef lambda1
  #undef lambda2
  #undef lambda3
  #undef pTable
  #undef pTable_x
  #undef pTable_y
  #undef pTable_xmin
  #undef pTable_xmax
  #undef pTable_xsum
  #undef pTable_ymin
  #undef pTable_ymax
  #undef pTable_ysum
  #undef pTable_dxmin
  #undef pTable_dxmax
  #undef pTable_dymin
  #undef pTable_dymax
  /* component Source=Source_gen() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(_Source->_rotation_absolute,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_transpose(_PG->_rotation_absolute, tr1);
    rot_mul(_Source->_rotation_absolute, tr1, _Source->_rotation_relative);
    _Source->_rotation_is_identity =  rot_test_identity(_Source->_rotation_relative);
    _Source->_position_absolute = coords_set(
      0, 0, 0);
    tc1 = coords_sub(_PG->_position_absolute, _Source->_position_absolute);
    _Source->_position_relative = rot_apply(_Source->_rotation_absolute, tc1);
  } /* Source=Source_gen() AT ROTATED */
  DEBUG_COMPONENT("Source", _Source->_position_absolute, _Source->_rotation_absolute);
  instrument->_position_absolute[2] = _Source->_position_absolute;
  instrument->_position_relative[2] = _Source->_position_relative;
  instrument->counter_N[2]  = instrument->counter_P[2] = instrument->counter_P2[2] = 0;
  instrument->counter_AbsorbProp[2]= 0;
  return(0);
} /* _Source_setpos */

/* component Guide0=Guide_gravity() SETTING, POSITION/ROTATION */
int _Guide0_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Guide0_setpos] component Guide0=Guide_gravity() SETTING [/usr/share/mcstas/3.0-dev/optics/Guide_gravity.comp:339]");
  stracpy(_Guide0->_name, "Guide0", 16384);
  stracpy(_Guide0->_type, "Guide_gravity", 16384);
  _Guide0->_index=3;
  _Guide0->_parameters.w1 = 0.03;
  #define w1 (_Guide0->_parameters.w1)
  _Guide0->_parameters.h1 = 0.055;
  #define h1 (_Guide0->_parameters.h1)
  _Guide0->_parameters.w2 = 0;
  #define w2 (_Guide0->_parameters.w2)
  _Guide0->_parameters.h2 = 0;
  #define h2 (_Guide0->_parameters.h2)
  _Guide0->_parameters.l = 2;
  #define l (_Guide0->_parameters.l)
  _Guide0->_parameters.R0 = 0.995;
  #define R0 (_Guide0->_parameters.R0)
  _Guide0->_parameters.Qc = 0.0218;
  #define Qc (_Guide0->_parameters.Qc)
  _Guide0->_parameters.alpha = 4.38;
  #define alpha (_Guide0->_parameters.alpha)
  _Guide0->_parameters.m = 1.2;
  #define m (_Guide0->_parameters.m)
  _Guide0->_parameters.W = 0.003;
  #define W (_Guide0->_parameters.W)
  _Guide0->_parameters.nslit = 1;
  #define nslit (_Guide0->_parameters.nslit)
  _Guide0->_parameters.d = 0.0005;
  #define d (_Guide0->_parameters.d)
  _Guide0->_parameters.mleft = -1;
  #define mleft (_Guide0->_parameters.mleft)
  _Guide0->_parameters.mright = -1;
  #define mright (_Guide0->_parameters.mright)
  _Guide0->_parameters.mtop = -1;
  #define mtop (_Guide0->_parameters.mtop)
  _Guide0->_parameters.mbottom = -1;
  #define mbottom (_Guide0->_parameters.mbottom)
  _Guide0->_parameters.nhslit = 1;
  #define nhslit (_Guide0->_parameters.nhslit)
  _Guide0->_parameters.G = 0;
  #define G (_Guide0->_parameters.G)
  _Guide0->_parameters.aleft = -1;
  #define aleft (_Guide0->_parameters.aleft)
  _Guide0->_parameters.aright = -1;
  #define aright (_Guide0->_parameters.aright)
  _Guide0->_parameters.atop = -1;
  #define atop (_Guide0->_parameters.atop)
  _Guide0->_parameters.abottom = -1;
  #define abottom (_Guide0->_parameters.abottom)
  _Guide0->_parameters.wavy = 0;
  #define wavy (_Guide0->_parameters.wavy)
  _Guide0->_parameters.wavy_z = 0;
  #define wavy_z (_Guide0->_parameters.wavy_z)
  _Guide0->_parameters.wavy_tb = 0;
  #define wavy_tb (_Guide0->_parameters.wavy_tb)
  _Guide0->_parameters.wavy_lr = 0;
  #define wavy_lr (_Guide0->_parameters.wavy_lr)
  _Guide0->_parameters.chamfers = 0;
  #define chamfers (_Guide0->_parameters.chamfers)
  _Guide0->_parameters.chamfers_z = 0;
  #define chamfers_z (_Guide0->_parameters.chamfers_z)
  _Guide0->_parameters.chamfers_lr = 0;
  #define chamfers_lr (_Guide0->_parameters.chamfers_lr)
  _Guide0->_parameters.chamfers_tb = 0;
  #define chamfers_tb (_Guide0->_parameters.chamfers_tb)
  _Guide0->_parameters.nelements = 1;
  #define nelements (_Guide0->_parameters.nelements)
  _Guide0->_parameters.nu = 0;
  #define nu (_Guide0->_parameters.nu)
  _Guide0->_parameters.phase = 0;
  #define phase (_Guide0->_parameters.phase)
  if("NULL" && strlen("NULL"))
    stracpy(_Guide0->_parameters.reflect, "NULL" ? "NULL" : "", 16384);
  else 
  _Guide0->_parameters.reflect[0]='\0';
  #define reflect (_Guide0->_parameters.reflect)

  #define GVars (_Guide0->_parameters.GVars)
  #define pTable (_Guide0->_parameters.pTable)

  #undef w1
  #undef h1
  #undef w2
  #undef h2
  #undef l
  #undef R0
  #undef Qc
  #undef alpha
  #undef m
  #undef W
  #undef nslit
  #undef d
  #undef mleft
  #undef mright
  #undef mtop
  #undef mbottom
  #undef nhslit
  #undef G
  #undef aleft
  #undef aright
  #undef atop
  #undef abottom
  #undef wavy
  #undef wavy_z
  #undef wavy_tb
  #undef wavy_lr
  #undef chamfers
  #undef chamfers_z
  #undef chamfers_lr
  #undef chamfers_tb
  #undef nelements
  #undef nu
  #undef phase
  #undef reflect
  #undef GVars
  #undef pTable
  /* component Guide0=Guide_gravity() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Source->_rotation_absolute, _Guide0->_rotation_absolute);
    rot_transpose(_Source->_rotation_absolute, tr1);
    rot_mul(_Guide0->_rotation_absolute, tr1, _Guide0->_rotation_relative);
    _Guide0->_rotation_is_identity =  rot_test_identity(_Guide0->_rotation_relative);
    tc1 = coords_set(
      0, 0, 2);
    rot_transpose(_Source->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Guide0->_position_absolute = coords_add(_Source->_position_absolute, tc2);
    tc1 = coords_sub(_Source->_position_absolute, _Guide0->_position_absolute);
    _Guide0->_position_relative = rot_apply(_Guide0->_rotation_absolute, tc1);
  } /* Guide0=Guide_gravity() AT ROTATED */
  DEBUG_COMPONENT("Guide0", _Guide0->_position_absolute, _Guide0->_rotation_absolute);
  instrument->_position_absolute[3] = _Guide0->_position_absolute;
  instrument->_position_relative[3] = _Guide0->_position_relative;
  instrument->counter_N[3]  = instrument->counter_P[3] = instrument->counter_P2[3] = 0;
  instrument->counter_AbsorbProp[3]= 0;
  return(0);
} /* _Guide0_setpos */

/* component Guide0_4=Guide_gravity() SETTING, POSITION/ROTATION */
int _Guide0_4_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Guide0_4_setpos] component Guide0_4=Guide_gravity() SETTING [/usr/share/mcstas/3.0-dev/optics/Guide_gravity.comp:339]");
  stracpy(_Guide0_4->_name, "Guide0_4", 16384);
  stracpy(_Guide0_4->_type, "Guide_gravity", 16384);
  _Guide0_4->_index=4;
  _Guide0_4->_parameters.w1 = 0.03;
  #define w1 (_Guide0_4->_parameters.w1)
  _Guide0_4->_parameters.h1 = 0.055;
  #define h1 (_Guide0_4->_parameters.h1)
  _Guide0_4->_parameters.w2 = 0;
  #define w2 (_Guide0_4->_parameters.w2)
  _Guide0_4->_parameters.h2 = 0;
  #define h2 (_Guide0_4->_parameters.h2)
  _Guide0_4->_parameters.l = 2;
  #define l (_Guide0_4->_parameters.l)
  _Guide0_4->_parameters.R0 = 0.995;
  #define R0 (_Guide0_4->_parameters.R0)
  _Guide0_4->_parameters.Qc = 0.0218;
  #define Qc (_Guide0_4->_parameters.Qc)
  _Guide0_4->_parameters.alpha = 4.38;
  #define alpha (_Guide0_4->_parameters.alpha)
  _Guide0_4->_parameters.m = 1.2;
  #define m (_Guide0_4->_parameters.m)
  _Guide0_4->_parameters.W = 0.003;
  #define W (_Guide0_4->_parameters.W)
  _Guide0_4->_parameters.nslit = 1;
  #define nslit (_Guide0_4->_parameters.nslit)
  _Guide0_4->_parameters.d = 0.0005;
  #define d (_Guide0_4->_parameters.d)
  _Guide0_4->_parameters.mleft = -1;
  #define mleft (_Guide0_4->_parameters.mleft)
  _Guide0_4->_parameters.mright = -1;
  #define mright (_Guide0_4->_parameters.mright)
  _Guide0_4->_parameters.mtop = -1;
  #define mtop (_Guide0_4->_parameters.mtop)
  _Guide0_4->_parameters.mbottom = -1;
  #define mbottom (_Guide0_4->_parameters.mbottom)
  _Guide0_4->_parameters.nhslit = 1;
  #define nhslit (_Guide0_4->_parameters.nhslit)
  _Guide0_4->_parameters.G = 0;
  #define G (_Guide0_4->_parameters.G)
  _Guide0_4->_parameters.aleft = -1;
  #define aleft (_Guide0_4->_parameters.aleft)
  _Guide0_4->_parameters.aright = -1;
  #define aright (_Guide0_4->_parameters.aright)
  _Guide0_4->_parameters.atop = -1;
  #define atop (_Guide0_4->_parameters.atop)
  _Guide0_4->_parameters.abottom = -1;
  #define abottom (_Guide0_4->_parameters.abottom)
  _Guide0_4->_parameters.wavy = 0;
  #define wavy (_Guide0_4->_parameters.wavy)
  _Guide0_4->_parameters.wavy_z = 0;
  #define wavy_z (_Guide0_4->_parameters.wavy_z)
  _Guide0_4->_parameters.wavy_tb = 0;
  #define wavy_tb (_Guide0_4->_parameters.wavy_tb)
  _Guide0_4->_parameters.wavy_lr = 0;
  #define wavy_lr (_Guide0_4->_parameters.wavy_lr)
  _Guide0_4->_parameters.chamfers = 0;
  #define chamfers (_Guide0_4->_parameters.chamfers)
  _Guide0_4->_parameters.chamfers_z = 0;
  #define chamfers_z (_Guide0_4->_parameters.chamfers_z)
  _Guide0_4->_parameters.chamfers_lr = 0;
  #define chamfers_lr (_Guide0_4->_parameters.chamfers_lr)
  _Guide0_4->_parameters.chamfers_tb = 0;
  #define chamfers_tb (_Guide0_4->_parameters.chamfers_tb)
  _Guide0_4->_parameters.nelements = 1;
  #define nelements (_Guide0_4->_parameters.nelements)
  _Guide0_4->_parameters.nu = 0;
  #define nu (_Guide0_4->_parameters.nu)
  _Guide0_4->_parameters.phase = 0;
  #define phase (_Guide0_4->_parameters.phase)
  if("NULL" && strlen("NULL"))
    stracpy(_Guide0_4->_parameters.reflect, "NULL" ? "NULL" : "", 16384);
  else 
  _Guide0_4->_parameters.reflect[0]='\0';
  #define reflect (_Guide0_4->_parameters.reflect)

  #define GVars (_Guide0_4->_parameters.GVars)
  #define pTable (_Guide0_4->_parameters.pTable)

  #undef w1
  #undef h1
  #undef w2
  #undef h2
  #undef l
  #undef R0
  #undef Qc
  #undef alpha
  #undef m
  #undef W
  #undef nslit
  #undef d
  #undef mleft
  #undef mright
  #undef mtop
  #undef mbottom
  #undef nhslit
  #undef G
  #undef aleft
  #undef aright
  #undef atop
  #undef abottom
  #undef wavy
  #undef wavy_z
  #undef wavy_tb
  #undef wavy_lr
  #undef chamfers
  #undef chamfers_z
  #undef chamfers_lr
  #undef chamfers_tb
  #undef nelements
  #undef nu
  #undef phase
  #undef reflect
  #undef GVars
  #undef pTable
  /* component Guide0_4=Guide_gravity() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (2.0 / 2500 * RAD2DEG)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _Guide0->_rotation_absolute, _Guide0_4->_rotation_absolute);
    rot_transpose(_Guide0->_rotation_absolute, tr1);
    rot_mul(_Guide0_4->_rotation_absolute, tr1, _Guide0_4->_rotation_relative);
    _Guide0_4->_rotation_is_identity =  rot_test_identity(_Guide0_4->_rotation_relative);
    tc1 = coords_set(
      0, 0, 2 + 1e-4);
    rot_transpose(_Guide0->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Guide0_4->_position_absolute = coords_add(_Guide0->_position_absolute, tc2);
    tc1 = coords_sub(_Guide0->_position_absolute, _Guide0_4->_position_absolute);
    _Guide0_4->_position_relative = rot_apply(_Guide0_4->_rotation_absolute, tc1);
  } /* Guide0_4=Guide_gravity() AT ROTATED */
  DEBUG_COMPONENT("Guide0_4", _Guide0_4->_position_absolute, _Guide0_4->_rotation_absolute);
  instrument->_position_absolute[4] = _Guide0_4->_position_absolute;
  instrument->_position_relative[4] = _Guide0_4->_position_relative;
  instrument->counter_N[4]  = instrument->counter_P[4] = instrument->counter_P2[4] = 0;
  instrument->counter_AbsorbProp[4]= 0;
  return(0);
} /* _Guide0_4_setpos */

/* component Guide0_5=Guide_gravity() SETTING, POSITION/ROTATION */
int _Guide0_5_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Guide0_5_setpos] component Guide0_5=Guide_gravity() SETTING [/usr/share/mcstas/3.0-dev/optics/Guide_gravity.comp:339]");
  stracpy(_Guide0_5->_name, "Guide0_5", 16384);
  stracpy(_Guide0_5->_type, "Guide_gravity", 16384);
  _Guide0_5->_index=5;
  _Guide0_5->_parameters.w1 = 0.03;
  #define w1 (_Guide0_5->_parameters.w1)
  _Guide0_5->_parameters.h1 = 0.055;
  #define h1 (_Guide0_5->_parameters.h1)
  _Guide0_5->_parameters.w2 = 0;
  #define w2 (_Guide0_5->_parameters.w2)
  _Guide0_5->_parameters.h2 = 0;
  #define h2 (_Guide0_5->_parameters.h2)
  _Guide0_5->_parameters.l = 2;
  #define l (_Guide0_5->_parameters.l)
  _Guide0_5->_parameters.R0 = 0.995;
  #define R0 (_Guide0_5->_parameters.R0)
  _Guide0_5->_parameters.Qc = 0.0218;
  #define Qc (_Guide0_5->_parameters.Qc)
  _Guide0_5->_parameters.alpha = 4.38;
  #define alpha (_Guide0_5->_parameters.alpha)
  _Guide0_5->_parameters.m = 1.2;
  #define m (_Guide0_5->_parameters.m)
  _Guide0_5->_parameters.W = 0.003;
  #define W (_Guide0_5->_parameters.W)
  _Guide0_5->_parameters.nslit = 1;
  #define nslit (_Guide0_5->_parameters.nslit)
  _Guide0_5->_parameters.d = 0.0005;
  #define d (_Guide0_5->_parameters.d)
  _Guide0_5->_parameters.mleft = -1;
  #define mleft (_Guide0_5->_parameters.mleft)
  _Guide0_5->_parameters.mright = -1;
  #define mright (_Guide0_5->_parameters.mright)
  _Guide0_5->_parameters.mtop = -1;
  #define mtop (_Guide0_5->_parameters.mtop)
  _Guide0_5->_parameters.mbottom = -1;
  #define mbottom (_Guide0_5->_parameters.mbottom)
  _Guide0_5->_parameters.nhslit = 1;
  #define nhslit (_Guide0_5->_parameters.nhslit)
  _Guide0_5->_parameters.G = 0;
  #define G (_Guide0_5->_parameters.G)
  _Guide0_5->_parameters.aleft = -1;
  #define aleft (_Guide0_5->_parameters.aleft)
  _Guide0_5->_parameters.aright = -1;
  #define aright (_Guide0_5->_parameters.aright)
  _Guide0_5->_parameters.atop = -1;
  #define atop (_Guide0_5->_parameters.atop)
  _Guide0_5->_parameters.abottom = -1;
  #define abottom (_Guide0_5->_parameters.abottom)
  _Guide0_5->_parameters.wavy = 0;
  #define wavy (_Guide0_5->_parameters.wavy)
  _Guide0_5->_parameters.wavy_z = 0;
  #define wavy_z (_Guide0_5->_parameters.wavy_z)
  _Guide0_5->_parameters.wavy_tb = 0;
  #define wavy_tb (_Guide0_5->_parameters.wavy_tb)
  _Guide0_5->_parameters.wavy_lr = 0;
  #define wavy_lr (_Guide0_5->_parameters.wavy_lr)
  _Guide0_5->_parameters.chamfers = 0;
  #define chamfers (_Guide0_5->_parameters.chamfers)
  _Guide0_5->_parameters.chamfers_z = 0;
  #define chamfers_z (_Guide0_5->_parameters.chamfers_z)
  _Guide0_5->_parameters.chamfers_lr = 0;
  #define chamfers_lr (_Guide0_5->_parameters.chamfers_lr)
  _Guide0_5->_parameters.chamfers_tb = 0;
  #define chamfers_tb (_Guide0_5->_parameters.chamfers_tb)
  _Guide0_5->_parameters.nelements = 1;
  #define nelements (_Guide0_5->_parameters.nelements)
  _Guide0_5->_parameters.nu = 0;
  #define nu (_Guide0_5->_parameters.nu)
  _Guide0_5->_parameters.phase = 0;
  #define phase (_Guide0_5->_parameters.phase)
  if("NULL" && strlen("NULL"))
    stracpy(_Guide0_5->_parameters.reflect, "NULL" ? "NULL" : "", 16384);
  else 
  _Guide0_5->_parameters.reflect[0]='\0';
  #define reflect (_Guide0_5->_parameters.reflect)

  #define GVars (_Guide0_5->_parameters.GVars)
  #define pTable (_Guide0_5->_parameters.pTable)

  #undef w1
  #undef h1
  #undef w2
  #undef h2
  #undef l
  #undef R0
  #undef Qc
  #undef alpha
  #undef m
  #undef W
  #undef nslit
  #undef d
  #undef mleft
  #undef mright
  #undef mtop
  #undef mbottom
  #undef nhslit
  #undef G
  #undef aleft
  #undef aright
  #undef atop
  #undef abottom
  #undef wavy
  #undef wavy_z
  #undef wavy_tb
  #undef wavy_lr
  #undef chamfers
  #undef chamfers_z
  #undef chamfers_lr
  #undef chamfers_tb
  #undef nelements
  #undef nu
  #undef phase
  #undef reflect
  #undef GVars
  #undef pTable
  /* component Guide0_5=Guide_gravity() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (2.0 / 2500 * RAD2DEG)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _Guide0_4->_rotation_absolute, _Guide0_5->_rotation_absolute);
    rot_transpose(_Guide0_4->_rotation_absolute, tr1);
    rot_mul(_Guide0_5->_rotation_absolute, tr1, _Guide0_5->_rotation_relative);
    _Guide0_5->_rotation_is_identity =  rot_test_identity(_Guide0_5->_rotation_relative);
    tc1 = coords_set(
      0, 0, 2 + 1e-4);
    rot_transpose(_Guide0_4->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Guide0_5->_position_absolute = coords_add(_Guide0_4->_position_absolute, tc2);
    tc1 = coords_sub(_Guide0_4->_position_absolute, _Guide0_5->_position_absolute);
    _Guide0_5->_position_relative = rot_apply(_Guide0_5->_rotation_absolute, tc1);
  } /* Guide0_5=Guide_gravity() AT ROTATED */
  DEBUG_COMPONENT("Guide0_5", _Guide0_5->_position_absolute, _Guide0_5->_rotation_absolute);
  instrument->_position_absolute[5] = _Guide0_5->_position_absolute;
  instrument->_position_relative[5] = _Guide0_5->_position_relative;
  instrument->counter_N[5]  = instrument->counter_P[5] = instrument->counter_P2[5] = 0;
  instrument->counter_AbsorbProp[5]= 0;
  return(0);
} /* _Guide0_5_setpos */

/* component Guide0_6=Guide_gravity() SETTING, POSITION/ROTATION */
int _Guide0_6_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Guide0_6_setpos] component Guide0_6=Guide_gravity() SETTING [/usr/share/mcstas/3.0-dev/optics/Guide_gravity.comp:339]");
  stracpy(_Guide0_6->_name, "Guide0_6", 16384);
  stracpy(_Guide0_6->_type, "Guide_gravity", 16384);
  _Guide0_6->_index=6;
  _Guide0_6->_parameters.w1 = 0.03;
  #define w1 (_Guide0_6->_parameters.w1)
  _Guide0_6->_parameters.h1 = 0.055;
  #define h1 (_Guide0_6->_parameters.h1)
  _Guide0_6->_parameters.w2 = 0;
  #define w2 (_Guide0_6->_parameters.w2)
  _Guide0_6->_parameters.h2 = 0;
  #define h2 (_Guide0_6->_parameters.h2)
  _Guide0_6->_parameters.l = 2;
  #define l (_Guide0_6->_parameters.l)
  _Guide0_6->_parameters.R0 = 0.995;
  #define R0 (_Guide0_6->_parameters.R0)
  _Guide0_6->_parameters.Qc = 0.0218;
  #define Qc (_Guide0_6->_parameters.Qc)
  _Guide0_6->_parameters.alpha = 4.38;
  #define alpha (_Guide0_6->_parameters.alpha)
  _Guide0_6->_parameters.m = 1.2;
  #define m (_Guide0_6->_parameters.m)
  _Guide0_6->_parameters.W = 0.003;
  #define W (_Guide0_6->_parameters.W)
  _Guide0_6->_parameters.nslit = 1;
  #define nslit (_Guide0_6->_parameters.nslit)
  _Guide0_6->_parameters.d = 0.0005;
  #define d (_Guide0_6->_parameters.d)
  _Guide0_6->_parameters.mleft = -1;
  #define mleft (_Guide0_6->_parameters.mleft)
  _Guide0_6->_parameters.mright = -1;
  #define mright (_Guide0_6->_parameters.mright)
  _Guide0_6->_parameters.mtop = -1;
  #define mtop (_Guide0_6->_parameters.mtop)
  _Guide0_6->_parameters.mbottom = -1;
  #define mbottom (_Guide0_6->_parameters.mbottom)
  _Guide0_6->_parameters.nhslit = 1;
  #define nhslit (_Guide0_6->_parameters.nhslit)
  _Guide0_6->_parameters.G = 0;
  #define G (_Guide0_6->_parameters.G)
  _Guide0_6->_parameters.aleft = -1;
  #define aleft (_Guide0_6->_parameters.aleft)
  _Guide0_6->_parameters.aright = -1;
  #define aright (_Guide0_6->_parameters.aright)
  _Guide0_6->_parameters.atop = -1;
  #define atop (_Guide0_6->_parameters.atop)
  _Guide0_6->_parameters.abottom = -1;
  #define abottom (_Guide0_6->_parameters.abottom)
  _Guide0_6->_parameters.wavy = 0;
  #define wavy (_Guide0_6->_parameters.wavy)
  _Guide0_6->_parameters.wavy_z = 0;
  #define wavy_z (_Guide0_6->_parameters.wavy_z)
  _Guide0_6->_parameters.wavy_tb = 0;
  #define wavy_tb (_Guide0_6->_parameters.wavy_tb)
  _Guide0_6->_parameters.wavy_lr = 0;
  #define wavy_lr (_Guide0_6->_parameters.wavy_lr)
  _Guide0_6->_parameters.chamfers = 0;
  #define chamfers (_Guide0_6->_parameters.chamfers)
  _Guide0_6->_parameters.chamfers_z = 0;
  #define chamfers_z (_Guide0_6->_parameters.chamfers_z)
  _Guide0_6->_parameters.chamfers_lr = 0;
  #define chamfers_lr (_Guide0_6->_parameters.chamfers_lr)
  _Guide0_6->_parameters.chamfers_tb = 0;
  #define chamfers_tb (_Guide0_6->_parameters.chamfers_tb)
  _Guide0_6->_parameters.nelements = 1;
  #define nelements (_Guide0_6->_parameters.nelements)
  _Guide0_6->_parameters.nu = 0;
  #define nu (_Guide0_6->_parameters.nu)
  _Guide0_6->_parameters.phase = 0;
  #define phase (_Guide0_6->_parameters.phase)
  if("NULL" && strlen("NULL"))
    stracpy(_Guide0_6->_parameters.reflect, "NULL" ? "NULL" : "", 16384);
  else 
  _Guide0_6->_parameters.reflect[0]='\0';
  #define reflect (_Guide0_6->_parameters.reflect)

  #define GVars (_Guide0_6->_parameters.GVars)
  #define pTable (_Guide0_6->_parameters.pTable)

  #undef w1
  #undef h1
  #undef w2
  #undef h2
  #undef l
  #undef R0
  #undef Qc
  #undef alpha
  #undef m
  #undef W
  #undef nslit
  #undef d
  #undef mleft
  #undef mright
  #undef mtop
  #undef mbottom
  #undef nhslit
  #undef G
  #undef aleft
  #undef aright
  #undef atop
  #undef abottom
  #undef wavy
  #undef wavy_z
  #undef wavy_tb
  #undef wavy_lr
  #undef chamfers
  #undef chamfers_z
  #undef chamfers_lr
  #undef chamfers_tb
  #undef nelements
  #undef nu
  #undef phase
  #undef reflect
  #undef GVars
  #undef pTable
  /* component Guide0_6=Guide_gravity() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (2.0 / 2500 * RAD2DEG)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _Guide0_5->_rotation_absolute, _Guide0_6->_rotation_absolute);
    rot_transpose(_Guide0_5->_rotation_absolute, tr1);
    rot_mul(_Guide0_6->_rotation_absolute, tr1, _Guide0_6->_rotation_relative);
    _Guide0_6->_rotation_is_identity =  rot_test_identity(_Guide0_6->_rotation_relative);
    tc1 = coords_set(
      0, 0, 2 + 1e-4);
    rot_transpose(_Guide0_5->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Guide0_6->_position_absolute = coords_add(_Guide0_5->_position_absolute, tc2);
    tc1 = coords_sub(_Guide0_5->_position_absolute, _Guide0_6->_position_absolute);
    _Guide0_6->_position_relative = rot_apply(_Guide0_6->_rotation_absolute, tc1);
  } /* Guide0_6=Guide_gravity() AT ROTATED */
  DEBUG_COMPONENT("Guide0_6", _Guide0_6->_position_absolute, _Guide0_6->_rotation_absolute);
  instrument->_position_absolute[6] = _Guide0_6->_position_absolute;
  instrument->_position_relative[6] = _Guide0_6->_position_relative;
  instrument->counter_N[6]  = instrument->counter_P[6] = instrument->counter_P2[6] = 0;
  instrument->counter_AbsorbProp[6]= 0;
  return(0);
} /* _Guide0_6_setpos */

/* component Guide0_7=Guide_gravity() SETTING, POSITION/ROTATION */
int _Guide0_7_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Guide0_7_setpos] component Guide0_7=Guide_gravity() SETTING [/usr/share/mcstas/3.0-dev/optics/Guide_gravity.comp:339]");
  stracpy(_Guide0_7->_name, "Guide0_7", 16384);
  stracpy(_Guide0_7->_type, "Guide_gravity", 16384);
  _Guide0_7->_index=7;
  _Guide0_7->_parameters.w1 = 0.03;
  #define w1 (_Guide0_7->_parameters.w1)
  _Guide0_7->_parameters.h1 = 0.055;
  #define h1 (_Guide0_7->_parameters.h1)
  _Guide0_7->_parameters.w2 = 0;
  #define w2 (_Guide0_7->_parameters.w2)
  _Guide0_7->_parameters.h2 = 0;
  #define h2 (_Guide0_7->_parameters.h2)
  _Guide0_7->_parameters.l = 2;
  #define l (_Guide0_7->_parameters.l)
  _Guide0_7->_parameters.R0 = 0.995;
  #define R0 (_Guide0_7->_parameters.R0)
  _Guide0_7->_parameters.Qc = 0.0218;
  #define Qc (_Guide0_7->_parameters.Qc)
  _Guide0_7->_parameters.alpha = 4.38;
  #define alpha (_Guide0_7->_parameters.alpha)
  _Guide0_7->_parameters.m = 1.2;
  #define m (_Guide0_7->_parameters.m)
  _Guide0_7->_parameters.W = 0.003;
  #define W (_Guide0_7->_parameters.W)
  _Guide0_7->_parameters.nslit = 1;
  #define nslit (_Guide0_7->_parameters.nslit)
  _Guide0_7->_parameters.d = 0.0005;
  #define d (_Guide0_7->_parameters.d)
  _Guide0_7->_parameters.mleft = -1;
  #define mleft (_Guide0_7->_parameters.mleft)
  _Guide0_7->_parameters.mright = -1;
  #define mright (_Guide0_7->_parameters.mright)
  _Guide0_7->_parameters.mtop = -1;
  #define mtop (_Guide0_7->_parameters.mtop)
  _Guide0_7->_parameters.mbottom = -1;
  #define mbottom (_Guide0_7->_parameters.mbottom)
  _Guide0_7->_parameters.nhslit = 1;
  #define nhslit (_Guide0_7->_parameters.nhslit)
  _Guide0_7->_parameters.G = 0;
  #define G (_Guide0_7->_parameters.G)
  _Guide0_7->_parameters.aleft = -1;
  #define aleft (_Guide0_7->_parameters.aleft)
  _Guide0_7->_parameters.aright = -1;
  #define aright (_Guide0_7->_parameters.aright)
  _Guide0_7->_parameters.atop = -1;
  #define atop (_Guide0_7->_parameters.atop)
  _Guide0_7->_parameters.abottom = -1;
  #define abottom (_Guide0_7->_parameters.abottom)
  _Guide0_7->_parameters.wavy = 0;
  #define wavy (_Guide0_7->_parameters.wavy)
  _Guide0_7->_parameters.wavy_z = 0;
  #define wavy_z (_Guide0_7->_parameters.wavy_z)
  _Guide0_7->_parameters.wavy_tb = 0;
  #define wavy_tb (_Guide0_7->_parameters.wavy_tb)
  _Guide0_7->_parameters.wavy_lr = 0;
  #define wavy_lr (_Guide0_7->_parameters.wavy_lr)
  _Guide0_7->_parameters.chamfers = 0;
  #define chamfers (_Guide0_7->_parameters.chamfers)
  _Guide0_7->_parameters.chamfers_z = 0;
  #define chamfers_z (_Guide0_7->_parameters.chamfers_z)
  _Guide0_7->_parameters.chamfers_lr = 0;
  #define chamfers_lr (_Guide0_7->_parameters.chamfers_lr)
  _Guide0_7->_parameters.chamfers_tb = 0;
  #define chamfers_tb (_Guide0_7->_parameters.chamfers_tb)
  _Guide0_7->_parameters.nelements = 1;
  #define nelements (_Guide0_7->_parameters.nelements)
  _Guide0_7->_parameters.nu = 0;
  #define nu (_Guide0_7->_parameters.nu)
  _Guide0_7->_parameters.phase = 0;
  #define phase (_Guide0_7->_parameters.phase)
  if("NULL" && strlen("NULL"))
    stracpy(_Guide0_7->_parameters.reflect, "NULL" ? "NULL" : "", 16384);
  else 
  _Guide0_7->_parameters.reflect[0]='\0';
  #define reflect (_Guide0_7->_parameters.reflect)

  #define GVars (_Guide0_7->_parameters.GVars)
  #define pTable (_Guide0_7->_parameters.pTable)

  #undef w1
  #undef h1
  #undef w2
  #undef h2
  #undef l
  #undef R0
  #undef Qc
  #undef alpha
  #undef m
  #undef W
  #undef nslit
  #undef d
  #undef mleft
  #undef mright
  #undef mtop
  #undef mbottom
  #undef nhslit
  #undef G
  #undef aleft
  #undef aright
  #undef atop
  #undef abottom
  #undef wavy
  #undef wavy_z
  #undef wavy_tb
  #undef wavy_lr
  #undef chamfers
  #undef chamfers_z
  #undef chamfers_lr
  #undef chamfers_tb
  #undef nelements
  #undef nu
  #undef phase
  #undef reflect
  #undef GVars
  #undef pTable
  /* component Guide0_7=Guide_gravity() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (2.0 / 2500 * RAD2DEG)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _Guide0_6->_rotation_absolute, _Guide0_7->_rotation_absolute);
    rot_transpose(_Guide0_6->_rotation_absolute, tr1);
    rot_mul(_Guide0_7->_rotation_absolute, tr1, _Guide0_7->_rotation_relative);
    _Guide0_7->_rotation_is_identity =  rot_test_identity(_Guide0_7->_rotation_relative);
    tc1 = coords_set(
      0, 0, 2 + 1e-4);
    rot_transpose(_Guide0_6->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Guide0_7->_position_absolute = coords_add(_Guide0_6->_position_absolute, tc2);
    tc1 = coords_sub(_Guide0_6->_position_absolute, _Guide0_7->_position_absolute);
    _Guide0_7->_position_relative = rot_apply(_Guide0_7->_rotation_absolute, tc1);
  } /* Guide0_7=Guide_gravity() AT ROTATED */
  DEBUG_COMPONENT("Guide0_7", _Guide0_7->_position_absolute, _Guide0_7->_rotation_absolute);
  instrument->_position_absolute[7] = _Guide0_7->_position_absolute;
  instrument->_position_relative[7] = _Guide0_7->_position_relative;
  instrument->counter_N[7]  = instrument->counter_P[7] = instrument->counter_P2[7] = 0;
  instrument->counter_AbsorbProp[7]= 0;
  return(0);
} /* _Guide0_7_setpos */

/* component Guide0_8=Guide_gravity() SETTING, POSITION/ROTATION */
int _Guide0_8_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Guide0_8_setpos] component Guide0_8=Guide_gravity() SETTING [/usr/share/mcstas/3.0-dev/optics/Guide_gravity.comp:339]");
  stracpy(_Guide0_8->_name, "Guide0_8", 16384);
  stracpy(_Guide0_8->_type, "Guide_gravity", 16384);
  _Guide0_8->_index=8;
  _Guide0_8->_parameters.w1 = 0.03;
  #define w1 (_Guide0_8->_parameters.w1)
  _Guide0_8->_parameters.h1 = 0.055;
  #define h1 (_Guide0_8->_parameters.h1)
  _Guide0_8->_parameters.w2 = 0;
  #define w2 (_Guide0_8->_parameters.w2)
  _Guide0_8->_parameters.h2 = 0;
  #define h2 (_Guide0_8->_parameters.h2)
  _Guide0_8->_parameters.l = 2;
  #define l (_Guide0_8->_parameters.l)
  _Guide0_8->_parameters.R0 = 0.995;
  #define R0 (_Guide0_8->_parameters.R0)
  _Guide0_8->_parameters.Qc = 0.0218;
  #define Qc (_Guide0_8->_parameters.Qc)
  _Guide0_8->_parameters.alpha = 4.38;
  #define alpha (_Guide0_8->_parameters.alpha)
  _Guide0_8->_parameters.m = 1.2;
  #define m (_Guide0_8->_parameters.m)
  _Guide0_8->_parameters.W = 0.003;
  #define W (_Guide0_8->_parameters.W)
  _Guide0_8->_parameters.nslit = 1;
  #define nslit (_Guide0_8->_parameters.nslit)
  _Guide0_8->_parameters.d = 0.0005;
  #define d (_Guide0_8->_parameters.d)
  _Guide0_8->_parameters.mleft = -1;
  #define mleft (_Guide0_8->_parameters.mleft)
  _Guide0_8->_parameters.mright = -1;
  #define mright (_Guide0_8->_parameters.mright)
  _Guide0_8->_parameters.mtop = -1;
  #define mtop (_Guide0_8->_parameters.mtop)
  _Guide0_8->_parameters.mbottom = -1;
  #define mbottom (_Guide0_8->_parameters.mbottom)
  _Guide0_8->_parameters.nhslit = 1;
  #define nhslit (_Guide0_8->_parameters.nhslit)
  _Guide0_8->_parameters.G = 0;
  #define G (_Guide0_8->_parameters.G)
  _Guide0_8->_parameters.aleft = -1;
  #define aleft (_Guide0_8->_parameters.aleft)
  _Guide0_8->_parameters.aright = -1;
  #define aright (_Guide0_8->_parameters.aright)
  _Guide0_8->_parameters.atop = -1;
  #define atop (_Guide0_8->_parameters.atop)
  _Guide0_8->_parameters.abottom = -1;
  #define abottom (_Guide0_8->_parameters.abottom)
  _Guide0_8->_parameters.wavy = 0;
  #define wavy (_Guide0_8->_parameters.wavy)
  _Guide0_8->_parameters.wavy_z = 0;
  #define wavy_z (_Guide0_8->_parameters.wavy_z)
  _Guide0_8->_parameters.wavy_tb = 0;
  #define wavy_tb (_Guide0_8->_parameters.wavy_tb)
  _Guide0_8->_parameters.wavy_lr = 0;
  #define wavy_lr (_Guide0_8->_parameters.wavy_lr)
  _Guide0_8->_parameters.chamfers = 0;
  #define chamfers (_Guide0_8->_parameters.chamfers)
  _Guide0_8->_parameters.chamfers_z = 0;
  #define chamfers_z (_Guide0_8->_parameters.chamfers_z)
  _Guide0_8->_parameters.chamfers_lr = 0;
  #define chamfers_lr (_Guide0_8->_parameters.chamfers_lr)
  _Guide0_8->_parameters.chamfers_tb = 0;
  #define chamfers_tb (_Guide0_8->_parameters.chamfers_tb)
  _Guide0_8->_parameters.nelements = 1;
  #define nelements (_Guide0_8->_parameters.nelements)
  _Guide0_8->_parameters.nu = 0;
  #define nu (_Guide0_8->_parameters.nu)
  _Guide0_8->_parameters.phase = 0;
  #define phase (_Guide0_8->_parameters.phase)
  if("NULL" && strlen("NULL"))
    stracpy(_Guide0_8->_parameters.reflect, "NULL" ? "NULL" : "", 16384);
  else 
  _Guide0_8->_parameters.reflect[0]='\0';
  #define reflect (_Guide0_8->_parameters.reflect)

  #define GVars (_Guide0_8->_parameters.GVars)
  #define pTable (_Guide0_8->_parameters.pTable)

  #undef w1
  #undef h1
  #undef w2
  #undef h2
  #undef l
  #undef R0
  #undef Qc
  #undef alpha
  #undef m
  #undef W
  #undef nslit
  #undef d
  #undef mleft
  #undef mright
  #undef mtop
  #undef mbottom
  #undef nhslit
  #undef G
  #undef aleft
  #undef aright
  #undef atop
  #undef abottom
  #undef wavy
  #undef wavy_z
  #undef wavy_tb
  #undef wavy_lr
  #undef chamfers
  #undef chamfers_z
  #undef chamfers_lr
  #undef chamfers_tb
  #undef nelements
  #undef nu
  #undef phase
  #undef reflect
  #undef GVars
  #undef pTable
  /* component Guide0_8=Guide_gravity() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (2.0 / 2500 * RAD2DEG)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _Guide0_7->_rotation_absolute, _Guide0_8->_rotation_absolute);
    rot_transpose(_Guide0_7->_rotation_absolute, tr1);
    rot_mul(_Guide0_8->_rotation_absolute, tr1, _Guide0_8->_rotation_relative);
    _Guide0_8->_rotation_is_identity =  rot_test_identity(_Guide0_8->_rotation_relative);
    tc1 = coords_set(
      0, 0, 2 + 1e-4);
    rot_transpose(_Guide0_7->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Guide0_8->_position_absolute = coords_add(_Guide0_7->_position_absolute, tc2);
    tc1 = coords_sub(_Guide0_7->_position_absolute, _Guide0_8->_position_absolute);
    _Guide0_8->_position_relative = rot_apply(_Guide0_8->_rotation_absolute, tc1);
  } /* Guide0_8=Guide_gravity() AT ROTATED */
  DEBUG_COMPONENT("Guide0_8", _Guide0_8->_position_absolute, _Guide0_8->_rotation_absolute);
  instrument->_position_absolute[8] = _Guide0_8->_position_absolute;
  instrument->_position_relative[8] = _Guide0_8->_position_relative;
  instrument->counter_N[8]  = instrument->counter_P[8] = instrument->counter_P2[8] = 0;
  instrument->counter_AbsorbProp[8]= 0;
  return(0);
} /* _Guide0_8_setpos */

/* component Guide0_9=Guide_gravity() SETTING, POSITION/ROTATION */
int _Guide0_9_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Guide0_9_setpos] component Guide0_9=Guide_gravity() SETTING [/usr/share/mcstas/3.0-dev/optics/Guide_gravity.comp:339]");
  stracpy(_Guide0_9->_name, "Guide0_9", 16384);
  stracpy(_Guide0_9->_type, "Guide_gravity", 16384);
  _Guide0_9->_index=9;
  _Guide0_9->_parameters.w1 = 0.03;
  #define w1 (_Guide0_9->_parameters.w1)
  _Guide0_9->_parameters.h1 = 0.055;
  #define h1 (_Guide0_9->_parameters.h1)
  _Guide0_9->_parameters.w2 = 0;
  #define w2 (_Guide0_9->_parameters.w2)
  _Guide0_9->_parameters.h2 = 0;
  #define h2 (_Guide0_9->_parameters.h2)
  _Guide0_9->_parameters.l = 2;
  #define l (_Guide0_9->_parameters.l)
  _Guide0_9->_parameters.R0 = 0.995;
  #define R0 (_Guide0_9->_parameters.R0)
  _Guide0_9->_parameters.Qc = 0.0218;
  #define Qc (_Guide0_9->_parameters.Qc)
  _Guide0_9->_parameters.alpha = 4.38;
  #define alpha (_Guide0_9->_parameters.alpha)
  _Guide0_9->_parameters.m = 1.2;
  #define m (_Guide0_9->_parameters.m)
  _Guide0_9->_parameters.W = 0.003;
  #define W (_Guide0_9->_parameters.W)
  _Guide0_9->_parameters.nslit = 1;
  #define nslit (_Guide0_9->_parameters.nslit)
  _Guide0_9->_parameters.d = 0.0005;
  #define d (_Guide0_9->_parameters.d)
  _Guide0_9->_parameters.mleft = -1;
  #define mleft (_Guide0_9->_parameters.mleft)
  _Guide0_9->_parameters.mright = -1;
  #define mright (_Guide0_9->_parameters.mright)
  _Guide0_9->_parameters.mtop = -1;
  #define mtop (_Guide0_9->_parameters.mtop)
  _Guide0_9->_parameters.mbottom = -1;
  #define mbottom (_Guide0_9->_parameters.mbottom)
  _Guide0_9->_parameters.nhslit = 1;
  #define nhslit (_Guide0_9->_parameters.nhslit)
  _Guide0_9->_parameters.G = 0;
  #define G (_Guide0_9->_parameters.G)
  _Guide0_9->_parameters.aleft = -1;
  #define aleft (_Guide0_9->_parameters.aleft)
  _Guide0_9->_parameters.aright = -1;
  #define aright (_Guide0_9->_parameters.aright)
  _Guide0_9->_parameters.atop = -1;
  #define atop (_Guide0_9->_parameters.atop)
  _Guide0_9->_parameters.abottom = -1;
  #define abottom (_Guide0_9->_parameters.abottom)
  _Guide0_9->_parameters.wavy = 0;
  #define wavy (_Guide0_9->_parameters.wavy)
  _Guide0_9->_parameters.wavy_z = 0;
  #define wavy_z (_Guide0_9->_parameters.wavy_z)
  _Guide0_9->_parameters.wavy_tb = 0;
  #define wavy_tb (_Guide0_9->_parameters.wavy_tb)
  _Guide0_9->_parameters.wavy_lr = 0;
  #define wavy_lr (_Guide0_9->_parameters.wavy_lr)
  _Guide0_9->_parameters.chamfers = 0;
  #define chamfers (_Guide0_9->_parameters.chamfers)
  _Guide0_9->_parameters.chamfers_z = 0;
  #define chamfers_z (_Guide0_9->_parameters.chamfers_z)
  _Guide0_9->_parameters.chamfers_lr = 0;
  #define chamfers_lr (_Guide0_9->_parameters.chamfers_lr)
  _Guide0_9->_parameters.chamfers_tb = 0;
  #define chamfers_tb (_Guide0_9->_parameters.chamfers_tb)
  _Guide0_9->_parameters.nelements = 1;
  #define nelements (_Guide0_9->_parameters.nelements)
  _Guide0_9->_parameters.nu = 0;
  #define nu (_Guide0_9->_parameters.nu)
  _Guide0_9->_parameters.phase = 0;
  #define phase (_Guide0_9->_parameters.phase)
  if("NULL" && strlen("NULL"))
    stracpy(_Guide0_9->_parameters.reflect, "NULL" ? "NULL" : "", 16384);
  else 
  _Guide0_9->_parameters.reflect[0]='\0';
  #define reflect (_Guide0_9->_parameters.reflect)

  #define GVars (_Guide0_9->_parameters.GVars)
  #define pTable (_Guide0_9->_parameters.pTable)

  #undef w1
  #undef h1
  #undef w2
  #undef h2
  #undef l
  #undef R0
  #undef Qc
  #undef alpha
  #undef m
  #undef W
  #undef nslit
  #undef d
  #undef mleft
  #undef mright
  #undef mtop
  #undef mbottom
  #undef nhslit
  #undef G
  #undef aleft
  #undef aright
  #undef atop
  #undef abottom
  #undef wavy
  #undef wavy_z
  #undef wavy_tb
  #undef wavy_lr
  #undef chamfers
  #undef chamfers_z
  #undef chamfers_lr
  #undef chamfers_tb
  #undef nelements
  #undef nu
  #undef phase
  #undef reflect
  #undef GVars
  #undef pTable
  /* component Guide0_9=Guide_gravity() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (2.0 / 2500 * RAD2DEG)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _Guide0_8->_rotation_absolute, _Guide0_9->_rotation_absolute);
    rot_transpose(_Guide0_8->_rotation_absolute, tr1);
    rot_mul(_Guide0_9->_rotation_absolute, tr1, _Guide0_9->_rotation_relative);
    _Guide0_9->_rotation_is_identity =  rot_test_identity(_Guide0_9->_rotation_relative);
    tc1 = coords_set(
      0, 0, 2 + 1e-4);
    rot_transpose(_Guide0_8->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Guide0_9->_position_absolute = coords_add(_Guide0_8->_position_absolute, tc2);
    tc1 = coords_sub(_Guide0_8->_position_absolute, _Guide0_9->_position_absolute);
    _Guide0_9->_position_relative = rot_apply(_Guide0_9->_rotation_absolute, tc1);
  } /* Guide0_9=Guide_gravity() AT ROTATED */
  DEBUG_COMPONENT("Guide0_9", _Guide0_9->_position_absolute, _Guide0_9->_rotation_absolute);
  instrument->_position_absolute[9] = _Guide0_9->_position_absolute;
  instrument->_position_relative[9] = _Guide0_9->_position_relative;
  instrument->counter_N[9]  = instrument->counter_P[9] = instrument->counter_P2[9] = 0;
  instrument->counter_AbsorbProp[9]= 0;
  return(0);
} /* _Guide0_9_setpos */

/* component Guide0_10=Guide_gravity() SETTING, POSITION/ROTATION */
int _Guide0_10_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Guide0_10_setpos] component Guide0_10=Guide_gravity() SETTING [/usr/share/mcstas/3.0-dev/optics/Guide_gravity.comp:339]");
  stracpy(_Guide0_10->_name, "Guide0_10", 16384);
  stracpy(_Guide0_10->_type, "Guide_gravity", 16384);
  _Guide0_10->_index=10;
  _Guide0_10->_parameters.w1 = 0.03;
  #define w1 (_Guide0_10->_parameters.w1)
  _Guide0_10->_parameters.h1 = 0.055;
  #define h1 (_Guide0_10->_parameters.h1)
  _Guide0_10->_parameters.w2 = 0;
  #define w2 (_Guide0_10->_parameters.w2)
  _Guide0_10->_parameters.h2 = 0;
  #define h2 (_Guide0_10->_parameters.h2)
  _Guide0_10->_parameters.l = 2;
  #define l (_Guide0_10->_parameters.l)
  _Guide0_10->_parameters.R0 = 0.995;
  #define R0 (_Guide0_10->_parameters.R0)
  _Guide0_10->_parameters.Qc = 0.0218;
  #define Qc (_Guide0_10->_parameters.Qc)
  _Guide0_10->_parameters.alpha = 4.38;
  #define alpha (_Guide0_10->_parameters.alpha)
  _Guide0_10->_parameters.m = 1.2;
  #define m (_Guide0_10->_parameters.m)
  _Guide0_10->_parameters.W = 0.003;
  #define W (_Guide0_10->_parameters.W)
  _Guide0_10->_parameters.nslit = 1;
  #define nslit (_Guide0_10->_parameters.nslit)
  _Guide0_10->_parameters.d = 0.0005;
  #define d (_Guide0_10->_parameters.d)
  _Guide0_10->_parameters.mleft = -1;
  #define mleft (_Guide0_10->_parameters.mleft)
  _Guide0_10->_parameters.mright = -1;
  #define mright (_Guide0_10->_parameters.mright)
  _Guide0_10->_parameters.mtop = -1;
  #define mtop (_Guide0_10->_parameters.mtop)
  _Guide0_10->_parameters.mbottom = -1;
  #define mbottom (_Guide0_10->_parameters.mbottom)
  _Guide0_10->_parameters.nhslit = 1;
  #define nhslit (_Guide0_10->_parameters.nhslit)
  _Guide0_10->_parameters.G = 0;
  #define G (_Guide0_10->_parameters.G)
  _Guide0_10->_parameters.aleft = -1;
  #define aleft (_Guide0_10->_parameters.aleft)
  _Guide0_10->_parameters.aright = -1;
  #define aright (_Guide0_10->_parameters.aright)
  _Guide0_10->_parameters.atop = -1;
  #define atop (_Guide0_10->_parameters.atop)
  _Guide0_10->_parameters.abottom = -1;
  #define abottom (_Guide0_10->_parameters.abottom)
  _Guide0_10->_parameters.wavy = 0;
  #define wavy (_Guide0_10->_parameters.wavy)
  _Guide0_10->_parameters.wavy_z = 0;
  #define wavy_z (_Guide0_10->_parameters.wavy_z)
  _Guide0_10->_parameters.wavy_tb = 0;
  #define wavy_tb (_Guide0_10->_parameters.wavy_tb)
  _Guide0_10->_parameters.wavy_lr = 0;
  #define wavy_lr (_Guide0_10->_parameters.wavy_lr)
  _Guide0_10->_parameters.chamfers = 0;
  #define chamfers (_Guide0_10->_parameters.chamfers)
  _Guide0_10->_parameters.chamfers_z = 0;
  #define chamfers_z (_Guide0_10->_parameters.chamfers_z)
  _Guide0_10->_parameters.chamfers_lr = 0;
  #define chamfers_lr (_Guide0_10->_parameters.chamfers_lr)
  _Guide0_10->_parameters.chamfers_tb = 0;
  #define chamfers_tb (_Guide0_10->_parameters.chamfers_tb)
  _Guide0_10->_parameters.nelements = 1;
  #define nelements (_Guide0_10->_parameters.nelements)
  _Guide0_10->_parameters.nu = 0;
  #define nu (_Guide0_10->_parameters.nu)
  _Guide0_10->_parameters.phase = 0;
  #define phase (_Guide0_10->_parameters.phase)
  if("NULL" && strlen("NULL"))
    stracpy(_Guide0_10->_parameters.reflect, "NULL" ? "NULL" : "", 16384);
  else 
  _Guide0_10->_parameters.reflect[0]='\0';
  #define reflect (_Guide0_10->_parameters.reflect)

  #define GVars (_Guide0_10->_parameters.GVars)
  #define pTable (_Guide0_10->_parameters.pTable)

  #undef w1
  #undef h1
  #undef w2
  #undef h2
  #undef l
  #undef R0
  #undef Qc
  #undef alpha
  #undef m
  #undef W
  #undef nslit
  #undef d
  #undef mleft
  #undef mright
  #undef mtop
  #undef mbottom
  #undef nhslit
  #undef G
  #undef aleft
  #undef aright
  #undef atop
  #undef abottom
  #undef wavy
  #undef wavy_z
  #undef wavy_tb
  #undef wavy_lr
  #undef chamfers
  #undef chamfers_z
  #undef chamfers_lr
  #undef chamfers_tb
  #undef nelements
  #undef nu
  #undef phase
  #undef reflect
  #undef GVars
  #undef pTable
  /* component Guide0_10=Guide_gravity() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (2.0 / 2500 * RAD2DEG)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _Guide0_9->_rotation_absolute, _Guide0_10->_rotation_absolute);
    rot_transpose(_Guide0_9->_rotation_absolute, tr1);
    rot_mul(_Guide0_10->_rotation_absolute, tr1, _Guide0_10->_rotation_relative);
    _Guide0_10->_rotation_is_identity =  rot_test_identity(_Guide0_10->_rotation_relative);
    tc1 = coords_set(
      0, 0, 2 + 1e-4);
    rot_transpose(_Guide0_9->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Guide0_10->_position_absolute = coords_add(_Guide0_9->_position_absolute, tc2);
    tc1 = coords_sub(_Guide0_9->_position_absolute, _Guide0_10->_position_absolute);
    _Guide0_10->_position_relative = rot_apply(_Guide0_10->_rotation_absolute, tc1);
  } /* Guide0_10=Guide_gravity() AT ROTATED */
  DEBUG_COMPONENT("Guide0_10", _Guide0_10->_position_absolute, _Guide0_10->_rotation_absolute);
  instrument->_position_absolute[10] = _Guide0_10->_position_absolute;
  instrument->_position_relative[10] = _Guide0_10->_position_relative;
  instrument->counter_N[10]  = instrument->counter_P[10] = instrument->counter_P2[10] = 0;
  instrument->counter_AbsorbProp[10]= 0;
  return(0);
} /* _Guide0_10_setpos */

/* component Guide0_11=Guide_gravity() SETTING, POSITION/ROTATION */
int _Guide0_11_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Guide0_11_setpos] component Guide0_11=Guide_gravity() SETTING [/usr/share/mcstas/3.0-dev/optics/Guide_gravity.comp:339]");
  stracpy(_Guide0_11->_name, "Guide0_11", 16384);
  stracpy(_Guide0_11->_type, "Guide_gravity", 16384);
  _Guide0_11->_index=11;
  _Guide0_11->_parameters.w1 = 0.03;
  #define w1 (_Guide0_11->_parameters.w1)
  _Guide0_11->_parameters.h1 = 0.055;
  #define h1 (_Guide0_11->_parameters.h1)
  _Guide0_11->_parameters.w2 = 0;
  #define w2 (_Guide0_11->_parameters.w2)
  _Guide0_11->_parameters.h2 = 0;
  #define h2 (_Guide0_11->_parameters.h2)
  _Guide0_11->_parameters.l = 2;
  #define l (_Guide0_11->_parameters.l)
  _Guide0_11->_parameters.R0 = 0.995;
  #define R0 (_Guide0_11->_parameters.R0)
  _Guide0_11->_parameters.Qc = 0.0218;
  #define Qc (_Guide0_11->_parameters.Qc)
  _Guide0_11->_parameters.alpha = 4.38;
  #define alpha (_Guide0_11->_parameters.alpha)
  _Guide0_11->_parameters.m = 1.2;
  #define m (_Guide0_11->_parameters.m)
  _Guide0_11->_parameters.W = 0.003;
  #define W (_Guide0_11->_parameters.W)
  _Guide0_11->_parameters.nslit = 1;
  #define nslit (_Guide0_11->_parameters.nslit)
  _Guide0_11->_parameters.d = 0.0005;
  #define d (_Guide0_11->_parameters.d)
  _Guide0_11->_parameters.mleft = -1;
  #define mleft (_Guide0_11->_parameters.mleft)
  _Guide0_11->_parameters.mright = -1;
  #define mright (_Guide0_11->_parameters.mright)
  _Guide0_11->_parameters.mtop = -1;
  #define mtop (_Guide0_11->_parameters.mtop)
  _Guide0_11->_parameters.mbottom = -1;
  #define mbottom (_Guide0_11->_parameters.mbottom)
  _Guide0_11->_parameters.nhslit = 1;
  #define nhslit (_Guide0_11->_parameters.nhslit)
  _Guide0_11->_parameters.G = 0;
  #define G (_Guide0_11->_parameters.G)
  _Guide0_11->_parameters.aleft = -1;
  #define aleft (_Guide0_11->_parameters.aleft)
  _Guide0_11->_parameters.aright = -1;
  #define aright (_Guide0_11->_parameters.aright)
  _Guide0_11->_parameters.atop = -1;
  #define atop (_Guide0_11->_parameters.atop)
  _Guide0_11->_parameters.abottom = -1;
  #define abottom (_Guide0_11->_parameters.abottom)
  _Guide0_11->_parameters.wavy = 0;
  #define wavy (_Guide0_11->_parameters.wavy)
  _Guide0_11->_parameters.wavy_z = 0;
  #define wavy_z (_Guide0_11->_parameters.wavy_z)
  _Guide0_11->_parameters.wavy_tb = 0;
  #define wavy_tb (_Guide0_11->_parameters.wavy_tb)
  _Guide0_11->_parameters.wavy_lr = 0;
  #define wavy_lr (_Guide0_11->_parameters.wavy_lr)
  _Guide0_11->_parameters.chamfers = 0;
  #define chamfers (_Guide0_11->_parameters.chamfers)
  _Guide0_11->_parameters.chamfers_z = 0;
  #define chamfers_z (_Guide0_11->_parameters.chamfers_z)
  _Guide0_11->_parameters.chamfers_lr = 0;
  #define chamfers_lr (_Guide0_11->_parameters.chamfers_lr)
  _Guide0_11->_parameters.chamfers_tb = 0;
  #define chamfers_tb (_Guide0_11->_parameters.chamfers_tb)
  _Guide0_11->_parameters.nelements = 1;
  #define nelements (_Guide0_11->_parameters.nelements)
  _Guide0_11->_parameters.nu = 0;
  #define nu (_Guide0_11->_parameters.nu)
  _Guide0_11->_parameters.phase = 0;
  #define phase (_Guide0_11->_parameters.phase)
  if("NULL" && strlen("NULL"))
    stracpy(_Guide0_11->_parameters.reflect, "NULL" ? "NULL" : "", 16384);
  else 
  _Guide0_11->_parameters.reflect[0]='\0';
  #define reflect (_Guide0_11->_parameters.reflect)

  #define GVars (_Guide0_11->_parameters.GVars)
  #define pTable (_Guide0_11->_parameters.pTable)

  #undef w1
  #undef h1
  #undef w2
  #undef h2
  #undef l
  #undef R0
  #undef Qc
  #undef alpha
  #undef m
  #undef W
  #undef nslit
  #undef d
  #undef mleft
  #undef mright
  #undef mtop
  #undef mbottom
  #undef nhslit
  #undef G
  #undef aleft
  #undef aright
  #undef atop
  #undef abottom
  #undef wavy
  #undef wavy_z
  #undef wavy_tb
  #undef wavy_lr
  #undef chamfers
  #undef chamfers_z
  #undef chamfers_lr
  #undef chamfers_tb
  #undef nelements
  #undef nu
  #undef phase
  #undef reflect
  #undef GVars
  #undef pTable
  /* component Guide0_11=Guide_gravity() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (2.0 / 2500 * RAD2DEG)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _Guide0_10->_rotation_absolute, _Guide0_11->_rotation_absolute);
    rot_transpose(_Guide0_10->_rotation_absolute, tr1);
    rot_mul(_Guide0_11->_rotation_absolute, tr1, _Guide0_11->_rotation_relative);
    _Guide0_11->_rotation_is_identity =  rot_test_identity(_Guide0_11->_rotation_relative);
    tc1 = coords_set(
      0, 0, 2 + 1e-4);
    rot_transpose(_Guide0_10->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Guide0_11->_position_absolute = coords_add(_Guide0_10->_position_absolute, tc2);
    tc1 = coords_sub(_Guide0_10->_position_absolute, _Guide0_11->_position_absolute);
    _Guide0_11->_position_relative = rot_apply(_Guide0_11->_rotation_absolute, tc1);
  } /* Guide0_11=Guide_gravity() AT ROTATED */
  DEBUG_COMPONENT("Guide0_11", _Guide0_11->_position_absolute, _Guide0_11->_rotation_absolute);
  instrument->_position_absolute[11] = _Guide0_11->_position_absolute;
  instrument->_position_relative[11] = _Guide0_11->_position_relative;
  instrument->counter_N[11]  = instrument->counter_P[11] = instrument->counter_P2[11] = 0;
  instrument->counter_AbsorbProp[11]= 0;
  return(0);
} /* _Guide0_11_setpos */

/* component Guide0_12=Guide_gravity() SETTING, POSITION/ROTATION */
int _Guide0_12_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Guide0_12_setpos] component Guide0_12=Guide_gravity() SETTING [/usr/share/mcstas/3.0-dev/optics/Guide_gravity.comp:339]");
  stracpy(_Guide0_12->_name, "Guide0_12", 16384);
  stracpy(_Guide0_12->_type, "Guide_gravity", 16384);
  _Guide0_12->_index=12;
  _Guide0_12->_parameters.w1 = 0.03;
  #define w1 (_Guide0_12->_parameters.w1)
  _Guide0_12->_parameters.h1 = 0.055;
  #define h1 (_Guide0_12->_parameters.h1)
  _Guide0_12->_parameters.w2 = 0;
  #define w2 (_Guide0_12->_parameters.w2)
  _Guide0_12->_parameters.h2 = 0;
  #define h2 (_Guide0_12->_parameters.h2)
  _Guide0_12->_parameters.l = 2;
  #define l (_Guide0_12->_parameters.l)
  _Guide0_12->_parameters.R0 = 0.995;
  #define R0 (_Guide0_12->_parameters.R0)
  _Guide0_12->_parameters.Qc = 0.0218;
  #define Qc (_Guide0_12->_parameters.Qc)
  _Guide0_12->_parameters.alpha = 4.38;
  #define alpha (_Guide0_12->_parameters.alpha)
  _Guide0_12->_parameters.m = 1.2;
  #define m (_Guide0_12->_parameters.m)
  _Guide0_12->_parameters.W = 0.003;
  #define W (_Guide0_12->_parameters.W)
  _Guide0_12->_parameters.nslit = 1;
  #define nslit (_Guide0_12->_parameters.nslit)
  _Guide0_12->_parameters.d = 0.0005;
  #define d (_Guide0_12->_parameters.d)
  _Guide0_12->_parameters.mleft = -1;
  #define mleft (_Guide0_12->_parameters.mleft)
  _Guide0_12->_parameters.mright = -1;
  #define mright (_Guide0_12->_parameters.mright)
  _Guide0_12->_parameters.mtop = -1;
  #define mtop (_Guide0_12->_parameters.mtop)
  _Guide0_12->_parameters.mbottom = -1;
  #define mbottom (_Guide0_12->_parameters.mbottom)
  _Guide0_12->_parameters.nhslit = 1;
  #define nhslit (_Guide0_12->_parameters.nhslit)
  _Guide0_12->_parameters.G = 0;
  #define G (_Guide0_12->_parameters.G)
  _Guide0_12->_parameters.aleft = -1;
  #define aleft (_Guide0_12->_parameters.aleft)
  _Guide0_12->_parameters.aright = -1;
  #define aright (_Guide0_12->_parameters.aright)
  _Guide0_12->_parameters.atop = -1;
  #define atop (_Guide0_12->_parameters.atop)
  _Guide0_12->_parameters.abottom = -1;
  #define abottom (_Guide0_12->_parameters.abottom)
  _Guide0_12->_parameters.wavy = 0;
  #define wavy (_Guide0_12->_parameters.wavy)
  _Guide0_12->_parameters.wavy_z = 0;
  #define wavy_z (_Guide0_12->_parameters.wavy_z)
  _Guide0_12->_parameters.wavy_tb = 0;
  #define wavy_tb (_Guide0_12->_parameters.wavy_tb)
  _Guide0_12->_parameters.wavy_lr = 0;
  #define wavy_lr (_Guide0_12->_parameters.wavy_lr)
  _Guide0_12->_parameters.chamfers = 0;
  #define chamfers (_Guide0_12->_parameters.chamfers)
  _Guide0_12->_parameters.chamfers_z = 0;
  #define chamfers_z (_Guide0_12->_parameters.chamfers_z)
  _Guide0_12->_parameters.chamfers_lr = 0;
  #define chamfers_lr (_Guide0_12->_parameters.chamfers_lr)
  _Guide0_12->_parameters.chamfers_tb = 0;
  #define chamfers_tb (_Guide0_12->_parameters.chamfers_tb)
  _Guide0_12->_parameters.nelements = 1;
  #define nelements (_Guide0_12->_parameters.nelements)
  _Guide0_12->_parameters.nu = 0;
  #define nu (_Guide0_12->_parameters.nu)
  _Guide0_12->_parameters.phase = 0;
  #define phase (_Guide0_12->_parameters.phase)
  if("NULL" && strlen("NULL"))
    stracpy(_Guide0_12->_parameters.reflect, "NULL" ? "NULL" : "", 16384);
  else 
  _Guide0_12->_parameters.reflect[0]='\0';
  #define reflect (_Guide0_12->_parameters.reflect)

  #define GVars (_Guide0_12->_parameters.GVars)
  #define pTable (_Guide0_12->_parameters.pTable)

  #undef w1
  #undef h1
  #undef w2
  #undef h2
  #undef l
  #undef R0
  #undef Qc
  #undef alpha
  #undef m
  #undef W
  #undef nslit
  #undef d
  #undef mleft
  #undef mright
  #undef mtop
  #undef mbottom
  #undef nhslit
  #undef G
  #undef aleft
  #undef aright
  #undef atop
  #undef abottom
  #undef wavy
  #undef wavy_z
  #undef wavy_tb
  #undef wavy_lr
  #undef chamfers
  #undef chamfers_z
  #undef chamfers_lr
  #undef chamfers_tb
  #undef nelements
  #undef nu
  #undef phase
  #undef reflect
  #undef GVars
  #undef pTable
  /* component Guide0_12=Guide_gravity() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (2.0 / 2500 * RAD2DEG)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _Guide0_11->_rotation_absolute, _Guide0_12->_rotation_absolute);
    rot_transpose(_Guide0_11->_rotation_absolute, tr1);
    rot_mul(_Guide0_12->_rotation_absolute, tr1, _Guide0_12->_rotation_relative);
    _Guide0_12->_rotation_is_identity =  rot_test_identity(_Guide0_12->_rotation_relative);
    tc1 = coords_set(
      0, 0, 2 + 1e-4);
    rot_transpose(_Guide0_11->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Guide0_12->_position_absolute = coords_add(_Guide0_11->_position_absolute, tc2);
    tc1 = coords_sub(_Guide0_11->_position_absolute, _Guide0_12->_position_absolute);
    _Guide0_12->_position_relative = rot_apply(_Guide0_12->_rotation_absolute, tc1);
  } /* Guide0_12=Guide_gravity() AT ROTATED */
  DEBUG_COMPONENT("Guide0_12", _Guide0_12->_position_absolute, _Guide0_12->_rotation_absolute);
  instrument->_position_absolute[12] = _Guide0_12->_position_absolute;
  instrument->_position_relative[12] = _Guide0_12->_position_relative;
  instrument->counter_N[12]  = instrument->counter_P[12] = instrument->counter_P2[12] = 0;
  instrument->counter_AbsorbProp[12]= 0;
  return(0);
} /* _Guide0_12_setpos */

/* component Guide0_13=Guide_gravity() SETTING, POSITION/ROTATION */
int _Guide0_13_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Guide0_13_setpos] component Guide0_13=Guide_gravity() SETTING [/usr/share/mcstas/3.0-dev/optics/Guide_gravity.comp:339]");
  stracpy(_Guide0_13->_name, "Guide0_13", 16384);
  stracpy(_Guide0_13->_type, "Guide_gravity", 16384);
  _Guide0_13->_index=13;
  _Guide0_13->_parameters.w1 = 0.03;
  #define w1 (_Guide0_13->_parameters.w1)
  _Guide0_13->_parameters.h1 = 0.055;
  #define h1 (_Guide0_13->_parameters.h1)
  _Guide0_13->_parameters.w2 = 0;
  #define w2 (_Guide0_13->_parameters.w2)
  _Guide0_13->_parameters.h2 = 0;
  #define h2 (_Guide0_13->_parameters.h2)
  _Guide0_13->_parameters.l = 2;
  #define l (_Guide0_13->_parameters.l)
  _Guide0_13->_parameters.R0 = 0.995;
  #define R0 (_Guide0_13->_parameters.R0)
  _Guide0_13->_parameters.Qc = 0.0218;
  #define Qc (_Guide0_13->_parameters.Qc)
  _Guide0_13->_parameters.alpha = 4.38;
  #define alpha (_Guide0_13->_parameters.alpha)
  _Guide0_13->_parameters.m = 1.2;
  #define m (_Guide0_13->_parameters.m)
  _Guide0_13->_parameters.W = 0.003;
  #define W (_Guide0_13->_parameters.W)
  _Guide0_13->_parameters.nslit = 1;
  #define nslit (_Guide0_13->_parameters.nslit)
  _Guide0_13->_parameters.d = 0.0005;
  #define d (_Guide0_13->_parameters.d)
  _Guide0_13->_parameters.mleft = -1;
  #define mleft (_Guide0_13->_parameters.mleft)
  _Guide0_13->_parameters.mright = -1;
  #define mright (_Guide0_13->_parameters.mright)
  _Guide0_13->_parameters.mtop = -1;
  #define mtop (_Guide0_13->_parameters.mtop)
  _Guide0_13->_parameters.mbottom = -1;
  #define mbottom (_Guide0_13->_parameters.mbottom)
  _Guide0_13->_parameters.nhslit = 1;
  #define nhslit (_Guide0_13->_parameters.nhslit)
  _Guide0_13->_parameters.G = 0;
  #define G (_Guide0_13->_parameters.G)
  _Guide0_13->_parameters.aleft = -1;
  #define aleft (_Guide0_13->_parameters.aleft)
  _Guide0_13->_parameters.aright = -1;
  #define aright (_Guide0_13->_parameters.aright)
  _Guide0_13->_parameters.atop = -1;
  #define atop (_Guide0_13->_parameters.atop)
  _Guide0_13->_parameters.abottom = -1;
  #define abottom (_Guide0_13->_parameters.abottom)
  _Guide0_13->_parameters.wavy = 0;
  #define wavy (_Guide0_13->_parameters.wavy)
  _Guide0_13->_parameters.wavy_z = 0;
  #define wavy_z (_Guide0_13->_parameters.wavy_z)
  _Guide0_13->_parameters.wavy_tb = 0;
  #define wavy_tb (_Guide0_13->_parameters.wavy_tb)
  _Guide0_13->_parameters.wavy_lr = 0;
  #define wavy_lr (_Guide0_13->_parameters.wavy_lr)
  _Guide0_13->_parameters.chamfers = 0;
  #define chamfers (_Guide0_13->_parameters.chamfers)
  _Guide0_13->_parameters.chamfers_z = 0;
  #define chamfers_z (_Guide0_13->_parameters.chamfers_z)
  _Guide0_13->_parameters.chamfers_lr = 0;
  #define chamfers_lr (_Guide0_13->_parameters.chamfers_lr)
  _Guide0_13->_parameters.chamfers_tb = 0;
  #define chamfers_tb (_Guide0_13->_parameters.chamfers_tb)
  _Guide0_13->_parameters.nelements = 1;
  #define nelements (_Guide0_13->_parameters.nelements)
  _Guide0_13->_parameters.nu = 0;
  #define nu (_Guide0_13->_parameters.nu)
  _Guide0_13->_parameters.phase = 0;
  #define phase (_Guide0_13->_parameters.phase)
  if("NULL" && strlen("NULL"))
    stracpy(_Guide0_13->_parameters.reflect, "NULL" ? "NULL" : "", 16384);
  else 
  _Guide0_13->_parameters.reflect[0]='\0';
  #define reflect (_Guide0_13->_parameters.reflect)

  #define GVars (_Guide0_13->_parameters.GVars)
  #define pTable (_Guide0_13->_parameters.pTable)

  #undef w1
  #undef h1
  #undef w2
  #undef h2
  #undef l
  #undef R0
  #undef Qc
  #undef alpha
  #undef m
  #undef W
  #undef nslit
  #undef d
  #undef mleft
  #undef mright
  #undef mtop
  #undef mbottom
  #undef nhslit
  #undef G
  #undef aleft
  #undef aright
  #undef atop
  #undef abottom
  #undef wavy
  #undef wavy_z
  #undef wavy_tb
  #undef wavy_lr
  #undef chamfers
  #undef chamfers_z
  #undef chamfers_lr
  #undef chamfers_tb
  #undef nelements
  #undef nu
  #undef phase
  #undef reflect
  #undef GVars
  #undef pTable
  /* component Guide0_13=Guide_gravity() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (2.0 / 2500 * RAD2DEG)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _Guide0_12->_rotation_absolute, _Guide0_13->_rotation_absolute);
    rot_transpose(_Guide0_12->_rotation_absolute, tr1);
    rot_mul(_Guide0_13->_rotation_absolute, tr1, _Guide0_13->_rotation_relative);
    _Guide0_13->_rotation_is_identity =  rot_test_identity(_Guide0_13->_rotation_relative);
    tc1 = coords_set(
      0, 0, 2 + 1e-4);
    rot_transpose(_Guide0_12->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Guide0_13->_position_absolute = coords_add(_Guide0_12->_position_absolute, tc2);
    tc1 = coords_sub(_Guide0_12->_position_absolute, _Guide0_13->_position_absolute);
    _Guide0_13->_position_relative = rot_apply(_Guide0_13->_rotation_absolute, tc1);
  } /* Guide0_13=Guide_gravity() AT ROTATED */
  DEBUG_COMPONENT("Guide0_13", _Guide0_13->_position_absolute, _Guide0_13->_rotation_absolute);
  instrument->_position_absolute[13] = _Guide0_13->_position_absolute;
  instrument->_position_relative[13] = _Guide0_13->_position_relative;
  instrument->counter_N[13]  = instrument->counter_P[13] = instrument->counter_P2[13] = 0;
  instrument->counter_AbsorbProp[13]= 0;
  return(0);
} /* _Guide0_13_setpos */

/* component Guide0_14=Guide_gravity() SETTING, POSITION/ROTATION */
int _Guide0_14_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Guide0_14_setpos] component Guide0_14=Guide_gravity() SETTING [/usr/share/mcstas/3.0-dev/optics/Guide_gravity.comp:339]");
  stracpy(_Guide0_14->_name, "Guide0_14", 16384);
  stracpy(_Guide0_14->_type, "Guide_gravity", 16384);
  _Guide0_14->_index=14;
  _Guide0_14->_parameters.w1 = 0.03;
  #define w1 (_Guide0_14->_parameters.w1)
  _Guide0_14->_parameters.h1 = 0.055;
  #define h1 (_Guide0_14->_parameters.h1)
  _Guide0_14->_parameters.w2 = 0;
  #define w2 (_Guide0_14->_parameters.w2)
  _Guide0_14->_parameters.h2 = 0;
  #define h2 (_Guide0_14->_parameters.h2)
  _Guide0_14->_parameters.l = 2;
  #define l (_Guide0_14->_parameters.l)
  _Guide0_14->_parameters.R0 = 0.995;
  #define R0 (_Guide0_14->_parameters.R0)
  _Guide0_14->_parameters.Qc = 0.0218;
  #define Qc (_Guide0_14->_parameters.Qc)
  _Guide0_14->_parameters.alpha = 4.38;
  #define alpha (_Guide0_14->_parameters.alpha)
  _Guide0_14->_parameters.m = 1.2;
  #define m (_Guide0_14->_parameters.m)
  _Guide0_14->_parameters.W = 0.003;
  #define W (_Guide0_14->_parameters.W)
  _Guide0_14->_parameters.nslit = 1;
  #define nslit (_Guide0_14->_parameters.nslit)
  _Guide0_14->_parameters.d = 0.0005;
  #define d (_Guide0_14->_parameters.d)
  _Guide0_14->_parameters.mleft = -1;
  #define mleft (_Guide0_14->_parameters.mleft)
  _Guide0_14->_parameters.mright = -1;
  #define mright (_Guide0_14->_parameters.mright)
  _Guide0_14->_parameters.mtop = -1;
  #define mtop (_Guide0_14->_parameters.mtop)
  _Guide0_14->_parameters.mbottom = -1;
  #define mbottom (_Guide0_14->_parameters.mbottom)
  _Guide0_14->_parameters.nhslit = 1;
  #define nhslit (_Guide0_14->_parameters.nhslit)
  _Guide0_14->_parameters.G = 0;
  #define G (_Guide0_14->_parameters.G)
  _Guide0_14->_parameters.aleft = -1;
  #define aleft (_Guide0_14->_parameters.aleft)
  _Guide0_14->_parameters.aright = -1;
  #define aright (_Guide0_14->_parameters.aright)
  _Guide0_14->_parameters.atop = -1;
  #define atop (_Guide0_14->_parameters.atop)
  _Guide0_14->_parameters.abottom = -1;
  #define abottom (_Guide0_14->_parameters.abottom)
  _Guide0_14->_parameters.wavy = 0;
  #define wavy (_Guide0_14->_parameters.wavy)
  _Guide0_14->_parameters.wavy_z = 0;
  #define wavy_z (_Guide0_14->_parameters.wavy_z)
  _Guide0_14->_parameters.wavy_tb = 0;
  #define wavy_tb (_Guide0_14->_parameters.wavy_tb)
  _Guide0_14->_parameters.wavy_lr = 0;
  #define wavy_lr (_Guide0_14->_parameters.wavy_lr)
  _Guide0_14->_parameters.chamfers = 0;
  #define chamfers (_Guide0_14->_parameters.chamfers)
  _Guide0_14->_parameters.chamfers_z = 0;
  #define chamfers_z (_Guide0_14->_parameters.chamfers_z)
  _Guide0_14->_parameters.chamfers_lr = 0;
  #define chamfers_lr (_Guide0_14->_parameters.chamfers_lr)
  _Guide0_14->_parameters.chamfers_tb = 0;
  #define chamfers_tb (_Guide0_14->_parameters.chamfers_tb)
  _Guide0_14->_parameters.nelements = 1;
  #define nelements (_Guide0_14->_parameters.nelements)
  _Guide0_14->_parameters.nu = 0;
  #define nu (_Guide0_14->_parameters.nu)
  _Guide0_14->_parameters.phase = 0;
  #define phase (_Guide0_14->_parameters.phase)
  if("NULL" && strlen("NULL"))
    stracpy(_Guide0_14->_parameters.reflect, "NULL" ? "NULL" : "", 16384);
  else 
  _Guide0_14->_parameters.reflect[0]='\0';
  #define reflect (_Guide0_14->_parameters.reflect)

  #define GVars (_Guide0_14->_parameters.GVars)
  #define pTable (_Guide0_14->_parameters.pTable)

  #undef w1
  #undef h1
  #undef w2
  #undef h2
  #undef l
  #undef R0
  #undef Qc
  #undef alpha
  #undef m
  #undef W
  #undef nslit
  #undef d
  #undef mleft
  #undef mright
  #undef mtop
  #undef mbottom
  #undef nhslit
  #undef G
  #undef aleft
  #undef aright
  #undef atop
  #undef abottom
  #undef wavy
  #undef wavy_z
  #undef wavy_tb
  #undef wavy_lr
  #undef chamfers
  #undef chamfers_z
  #undef chamfers_lr
  #undef chamfers_tb
  #undef nelements
  #undef nu
  #undef phase
  #undef reflect
  #undef GVars
  #undef pTable
  /* component Guide0_14=Guide_gravity() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (2.0 / 2500 * RAD2DEG)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _Guide0_13->_rotation_absolute, _Guide0_14->_rotation_absolute);
    rot_transpose(_Guide0_13->_rotation_absolute, tr1);
    rot_mul(_Guide0_14->_rotation_absolute, tr1, _Guide0_14->_rotation_relative);
    _Guide0_14->_rotation_is_identity =  rot_test_identity(_Guide0_14->_rotation_relative);
    tc1 = coords_set(
      0, 0, 2 + 1e-4);
    rot_transpose(_Guide0_13->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Guide0_14->_position_absolute = coords_add(_Guide0_13->_position_absolute, tc2);
    tc1 = coords_sub(_Guide0_13->_position_absolute, _Guide0_14->_position_absolute);
    _Guide0_14->_position_relative = rot_apply(_Guide0_14->_rotation_absolute, tc1);
  } /* Guide0_14=Guide_gravity() AT ROTATED */
  DEBUG_COMPONENT("Guide0_14", _Guide0_14->_position_absolute, _Guide0_14->_rotation_absolute);
  instrument->_position_absolute[14] = _Guide0_14->_position_absolute;
  instrument->_position_relative[14] = _Guide0_14->_position_relative;
  instrument->counter_N[14]  = instrument->counter_P[14] = instrument->counter_P2[14] = 0;
  instrument->counter_AbsorbProp[14]= 0;
  return(0);
} /* _Guide0_14_setpos */

/* component Guide0_15=Guide_gravity() SETTING, POSITION/ROTATION */
int _Guide0_15_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Guide0_15_setpos] component Guide0_15=Guide_gravity() SETTING [/usr/share/mcstas/3.0-dev/optics/Guide_gravity.comp:339]");
  stracpy(_Guide0_15->_name, "Guide0_15", 16384);
  stracpy(_Guide0_15->_type, "Guide_gravity", 16384);
  _Guide0_15->_index=15;
  _Guide0_15->_parameters.w1 = 0.03;
  #define w1 (_Guide0_15->_parameters.w1)
  _Guide0_15->_parameters.h1 = 0.055;
  #define h1 (_Guide0_15->_parameters.h1)
  _Guide0_15->_parameters.w2 = 0;
  #define w2 (_Guide0_15->_parameters.w2)
  _Guide0_15->_parameters.h2 = 0;
  #define h2 (_Guide0_15->_parameters.h2)
  _Guide0_15->_parameters.l = 2;
  #define l (_Guide0_15->_parameters.l)
  _Guide0_15->_parameters.R0 = 0.995;
  #define R0 (_Guide0_15->_parameters.R0)
  _Guide0_15->_parameters.Qc = 0.0218;
  #define Qc (_Guide0_15->_parameters.Qc)
  _Guide0_15->_parameters.alpha = 4.38;
  #define alpha (_Guide0_15->_parameters.alpha)
  _Guide0_15->_parameters.m = 1.2;
  #define m (_Guide0_15->_parameters.m)
  _Guide0_15->_parameters.W = 0.003;
  #define W (_Guide0_15->_parameters.W)
  _Guide0_15->_parameters.nslit = 1;
  #define nslit (_Guide0_15->_parameters.nslit)
  _Guide0_15->_parameters.d = 0.0005;
  #define d (_Guide0_15->_parameters.d)
  _Guide0_15->_parameters.mleft = -1;
  #define mleft (_Guide0_15->_parameters.mleft)
  _Guide0_15->_parameters.mright = -1;
  #define mright (_Guide0_15->_parameters.mright)
  _Guide0_15->_parameters.mtop = -1;
  #define mtop (_Guide0_15->_parameters.mtop)
  _Guide0_15->_parameters.mbottom = -1;
  #define mbottom (_Guide0_15->_parameters.mbottom)
  _Guide0_15->_parameters.nhslit = 1;
  #define nhslit (_Guide0_15->_parameters.nhslit)
  _Guide0_15->_parameters.G = 0;
  #define G (_Guide0_15->_parameters.G)
  _Guide0_15->_parameters.aleft = -1;
  #define aleft (_Guide0_15->_parameters.aleft)
  _Guide0_15->_parameters.aright = -1;
  #define aright (_Guide0_15->_parameters.aright)
  _Guide0_15->_parameters.atop = -1;
  #define atop (_Guide0_15->_parameters.atop)
  _Guide0_15->_parameters.abottom = -1;
  #define abottom (_Guide0_15->_parameters.abottom)
  _Guide0_15->_parameters.wavy = 0;
  #define wavy (_Guide0_15->_parameters.wavy)
  _Guide0_15->_parameters.wavy_z = 0;
  #define wavy_z (_Guide0_15->_parameters.wavy_z)
  _Guide0_15->_parameters.wavy_tb = 0;
  #define wavy_tb (_Guide0_15->_parameters.wavy_tb)
  _Guide0_15->_parameters.wavy_lr = 0;
  #define wavy_lr (_Guide0_15->_parameters.wavy_lr)
  _Guide0_15->_parameters.chamfers = 0;
  #define chamfers (_Guide0_15->_parameters.chamfers)
  _Guide0_15->_parameters.chamfers_z = 0;
  #define chamfers_z (_Guide0_15->_parameters.chamfers_z)
  _Guide0_15->_parameters.chamfers_lr = 0;
  #define chamfers_lr (_Guide0_15->_parameters.chamfers_lr)
  _Guide0_15->_parameters.chamfers_tb = 0;
  #define chamfers_tb (_Guide0_15->_parameters.chamfers_tb)
  _Guide0_15->_parameters.nelements = 1;
  #define nelements (_Guide0_15->_parameters.nelements)
  _Guide0_15->_parameters.nu = 0;
  #define nu (_Guide0_15->_parameters.nu)
  _Guide0_15->_parameters.phase = 0;
  #define phase (_Guide0_15->_parameters.phase)
  if("NULL" && strlen("NULL"))
    stracpy(_Guide0_15->_parameters.reflect, "NULL" ? "NULL" : "", 16384);
  else 
  _Guide0_15->_parameters.reflect[0]='\0';
  #define reflect (_Guide0_15->_parameters.reflect)

  #define GVars (_Guide0_15->_parameters.GVars)
  #define pTable (_Guide0_15->_parameters.pTable)

  #undef w1
  #undef h1
  #undef w2
  #undef h2
  #undef l
  #undef R0
  #undef Qc
  #undef alpha
  #undef m
  #undef W
  #undef nslit
  #undef d
  #undef mleft
  #undef mright
  #undef mtop
  #undef mbottom
  #undef nhslit
  #undef G
  #undef aleft
  #undef aright
  #undef atop
  #undef abottom
  #undef wavy
  #undef wavy_z
  #undef wavy_tb
  #undef wavy_lr
  #undef chamfers
  #undef chamfers_z
  #undef chamfers_lr
  #undef chamfers_tb
  #undef nelements
  #undef nu
  #undef phase
  #undef reflect
  #undef GVars
  #undef pTable
  /* component Guide0_15=Guide_gravity() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (2.0 / 2500 * RAD2DEG)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _Guide0_14->_rotation_absolute, _Guide0_15->_rotation_absolute);
    rot_transpose(_Guide0_14->_rotation_absolute, tr1);
    rot_mul(_Guide0_15->_rotation_absolute, tr1, _Guide0_15->_rotation_relative);
    _Guide0_15->_rotation_is_identity =  rot_test_identity(_Guide0_15->_rotation_relative);
    tc1 = coords_set(
      0, 0, 2 + 1e-4);
    rot_transpose(_Guide0_14->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Guide0_15->_position_absolute = coords_add(_Guide0_14->_position_absolute, tc2);
    tc1 = coords_sub(_Guide0_14->_position_absolute, _Guide0_15->_position_absolute);
    _Guide0_15->_position_relative = rot_apply(_Guide0_15->_rotation_absolute, tc1);
  } /* Guide0_15=Guide_gravity() AT ROTATED */
  DEBUG_COMPONENT("Guide0_15", _Guide0_15->_position_absolute, _Guide0_15->_rotation_absolute);
  instrument->_position_absolute[15] = _Guide0_15->_position_absolute;
  instrument->_position_relative[15] = _Guide0_15->_position_relative;
  instrument->counter_N[15]  = instrument->counter_P[15] = instrument->counter_P2[15] = 0;
  instrument->counter_AbsorbProp[15]= 0;
  return(0);
} /* _Guide0_15_setpos */

/* component Guide0_16=Guide_gravity() SETTING, POSITION/ROTATION */
int _Guide0_16_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Guide0_16_setpos] component Guide0_16=Guide_gravity() SETTING [/usr/share/mcstas/3.0-dev/optics/Guide_gravity.comp:339]");
  stracpy(_Guide0_16->_name, "Guide0_16", 16384);
  stracpy(_Guide0_16->_type, "Guide_gravity", 16384);
  _Guide0_16->_index=16;
  _Guide0_16->_parameters.w1 = 0.03;
  #define w1 (_Guide0_16->_parameters.w1)
  _Guide0_16->_parameters.h1 = 0.055;
  #define h1 (_Guide0_16->_parameters.h1)
  _Guide0_16->_parameters.w2 = 0;
  #define w2 (_Guide0_16->_parameters.w2)
  _Guide0_16->_parameters.h2 = 0;
  #define h2 (_Guide0_16->_parameters.h2)
  _Guide0_16->_parameters.l = 2;
  #define l (_Guide0_16->_parameters.l)
  _Guide0_16->_parameters.R0 = 0.995;
  #define R0 (_Guide0_16->_parameters.R0)
  _Guide0_16->_parameters.Qc = 0.0218;
  #define Qc (_Guide0_16->_parameters.Qc)
  _Guide0_16->_parameters.alpha = 4.38;
  #define alpha (_Guide0_16->_parameters.alpha)
  _Guide0_16->_parameters.m = 1.2;
  #define m (_Guide0_16->_parameters.m)
  _Guide0_16->_parameters.W = 0.003;
  #define W (_Guide0_16->_parameters.W)
  _Guide0_16->_parameters.nslit = 1;
  #define nslit (_Guide0_16->_parameters.nslit)
  _Guide0_16->_parameters.d = 0.0005;
  #define d (_Guide0_16->_parameters.d)
  _Guide0_16->_parameters.mleft = -1;
  #define mleft (_Guide0_16->_parameters.mleft)
  _Guide0_16->_parameters.mright = -1;
  #define mright (_Guide0_16->_parameters.mright)
  _Guide0_16->_parameters.mtop = -1;
  #define mtop (_Guide0_16->_parameters.mtop)
  _Guide0_16->_parameters.mbottom = -1;
  #define mbottom (_Guide0_16->_parameters.mbottom)
  _Guide0_16->_parameters.nhslit = 1;
  #define nhslit (_Guide0_16->_parameters.nhslit)
  _Guide0_16->_parameters.G = 0;
  #define G (_Guide0_16->_parameters.G)
  _Guide0_16->_parameters.aleft = -1;
  #define aleft (_Guide0_16->_parameters.aleft)
  _Guide0_16->_parameters.aright = -1;
  #define aright (_Guide0_16->_parameters.aright)
  _Guide0_16->_parameters.atop = -1;
  #define atop (_Guide0_16->_parameters.atop)
  _Guide0_16->_parameters.abottom = -1;
  #define abottom (_Guide0_16->_parameters.abottom)
  _Guide0_16->_parameters.wavy = 0;
  #define wavy (_Guide0_16->_parameters.wavy)
  _Guide0_16->_parameters.wavy_z = 0;
  #define wavy_z (_Guide0_16->_parameters.wavy_z)
  _Guide0_16->_parameters.wavy_tb = 0;
  #define wavy_tb (_Guide0_16->_parameters.wavy_tb)
  _Guide0_16->_parameters.wavy_lr = 0;
  #define wavy_lr (_Guide0_16->_parameters.wavy_lr)
  _Guide0_16->_parameters.chamfers = 0;
  #define chamfers (_Guide0_16->_parameters.chamfers)
  _Guide0_16->_parameters.chamfers_z = 0;
  #define chamfers_z (_Guide0_16->_parameters.chamfers_z)
  _Guide0_16->_parameters.chamfers_lr = 0;
  #define chamfers_lr (_Guide0_16->_parameters.chamfers_lr)
  _Guide0_16->_parameters.chamfers_tb = 0;
  #define chamfers_tb (_Guide0_16->_parameters.chamfers_tb)
  _Guide0_16->_parameters.nelements = 1;
  #define nelements (_Guide0_16->_parameters.nelements)
  _Guide0_16->_parameters.nu = 0;
  #define nu (_Guide0_16->_parameters.nu)
  _Guide0_16->_parameters.phase = 0;
  #define phase (_Guide0_16->_parameters.phase)
  if("NULL" && strlen("NULL"))
    stracpy(_Guide0_16->_parameters.reflect, "NULL" ? "NULL" : "", 16384);
  else 
  _Guide0_16->_parameters.reflect[0]='\0';
  #define reflect (_Guide0_16->_parameters.reflect)

  #define GVars (_Guide0_16->_parameters.GVars)
  #define pTable (_Guide0_16->_parameters.pTable)

  #undef w1
  #undef h1
  #undef w2
  #undef h2
  #undef l
  #undef R0
  #undef Qc
  #undef alpha
  #undef m
  #undef W
  #undef nslit
  #undef d
  #undef mleft
  #undef mright
  #undef mtop
  #undef mbottom
  #undef nhslit
  #undef G
  #undef aleft
  #undef aright
  #undef atop
  #undef abottom
  #undef wavy
  #undef wavy_z
  #undef wavy_tb
  #undef wavy_lr
  #undef chamfers
  #undef chamfers_z
  #undef chamfers_lr
  #undef chamfers_tb
  #undef nelements
  #undef nu
  #undef phase
  #undef reflect
  #undef GVars
  #undef pTable
  /* component Guide0_16=Guide_gravity() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (2.0 / 2500 * RAD2DEG)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _Guide0_15->_rotation_absolute, _Guide0_16->_rotation_absolute);
    rot_transpose(_Guide0_15->_rotation_absolute, tr1);
    rot_mul(_Guide0_16->_rotation_absolute, tr1, _Guide0_16->_rotation_relative);
    _Guide0_16->_rotation_is_identity =  rot_test_identity(_Guide0_16->_rotation_relative);
    tc1 = coords_set(
      0, 0, 2 + 1e-4);
    rot_transpose(_Guide0_15->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Guide0_16->_position_absolute = coords_add(_Guide0_15->_position_absolute, tc2);
    tc1 = coords_sub(_Guide0_15->_position_absolute, _Guide0_16->_position_absolute);
    _Guide0_16->_position_relative = rot_apply(_Guide0_16->_rotation_absolute, tc1);
  } /* Guide0_16=Guide_gravity() AT ROTATED */
  DEBUG_COMPONENT("Guide0_16", _Guide0_16->_position_absolute, _Guide0_16->_rotation_absolute);
  instrument->_position_absolute[16] = _Guide0_16->_position_absolute;
  instrument->_position_relative[16] = _Guide0_16->_position_relative;
  instrument->counter_N[16]  = instrument->counter_P[16] = instrument->counter_P2[16] = 0;
  instrument->counter_AbsorbProp[16]= 0;
  return(0);
} /* _Guide0_16_setpos */

/* component Guide0_17=Guide_gravity() SETTING, POSITION/ROTATION */
int _Guide0_17_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Guide0_17_setpos] component Guide0_17=Guide_gravity() SETTING [/usr/share/mcstas/3.0-dev/optics/Guide_gravity.comp:339]");
  stracpy(_Guide0_17->_name, "Guide0_17", 16384);
  stracpy(_Guide0_17->_type, "Guide_gravity", 16384);
  _Guide0_17->_index=17;
  _Guide0_17->_parameters.w1 = 0.03;
  #define w1 (_Guide0_17->_parameters.w1)
  _Guide0_17->_parameters.h1 = 0.055;
  #define h1 (_Guide0_17->_parameters.h1)
  _Guide0_17->_parameters.w2 = 0;
  #define w2 (_Guide0_17->_parameters.w2)
  _Guide0_17->_parameters.h2 = 0;
  #define h2 (_Guide0_17->_parameters.h2)
  _Guide0_17->_parameters.l = 2;
  #define l (_Guide0_17->_parameters.l)
  _Guide0_17->_parameters.R0 = 0.995;
  #define R0 (_Guide0_17->_parameters.R0)
  _Guide0_17->_parameters.Qc = 0.0218;
  #define Qc (_Guide0_17->_parameters.Qc)
  _Guide0_17->_parameters.alpha = 4.38;
  #define alpha (_Guide0_17->_parameters.alpha)
  _Guide0_17->_parameters.m = 1.2;
  #define m (_Guide0_17->_parameters.m)
  _Guide0_17->_parameters.W = 0.003;
  #define W (_Guide0_17->_parameters.W)
  _Guide0_17->_parameters.nslit = 1;
  #define nslit (_Guide0_17->_parameters.nslit)
  _Guide0_17->_parameters.d = 0.0005;
  #define d (_Guide0_17->_parameters.d)
  _Guide0_17->_parameters.mleft = -1;
  #define mleft (_Guide0_17->_parameters.mleft)
  _Guide0_17->_parameters.mright = -1;
  #define mright (_Guide0_17->_parameters.mright)
  _Guide0_17->_parameters.mtop = -1;
  #define mtop (_Guide0_17->_parameters.mtop)
  _Guide0_17->_parameters.mbottom = -1;
  #define mbottom (_Guide0_17->_parameters.mbottom)
  _Guide0_17->_parameters.nhslit = 1;
  #define nhslit (_Guide0_17->_parameters.nhslit)
  _Guide0_17->_parameters.G = 0;
  #define G (_Guide0_17->_parameters.G)
  _Guide0_17->_parameters.aleft = -1;
  #define aleft (_Guide0_17->_parameters.aleft)
  _Guide0_17->_parameters.aright = -1;
  #define aright (_Guide0_17->_parameters.aright)
  _Guide0_17->_parameters.atop = -1;
  #define atop (_Guide0_17->_parameters.atop)
  _Guide0_17->_parameters.abottom = -1;
  #define abottom (_Guide0_17->_parameters.abottom)
  _Guide0_17->_parameters.wavy = 0;
  #define wavy (_Guide0_17->_parameters.wavy)
  _Guide0_17->_parameters.wavy_z = 0;
  #define wavy_z (_Guide0_17->_parameters.wavy_z)
  _Guide0_17->_parameters.wavy_tb = 0;
  #define wavy_tb (_Guide0_17->_parameters.wavy_tb)
  _Guide0_17->_parameters.wavy_lr = 0;
  #define wavy_lr (_Guide0_17->_parameters.wavy_lr)
  _Guide0_17->_parameters.chamfers = 0;
  #define chamfers (_Guide0_17->_parameters.chamfers)
  _Guide0_17->_parameters.chamfers_z = 0;
  #define chamfers_z (_Guide0_17->_parameters.chamfers_z)
  _Guide0_17->_parameters.chamfers_lr = 0;
  #define chamfers_lr (_Guide0_17->_parameters.chamfers_lr)
  _Guide0_17->_parameters.chamfers_tb = 0;
  #define chamfers_tb (_Guide0_17->_parameters.chamfers_tb)
  _Guide0_17->_parameters.nelements = 1;
  #define nelements (_Guide0_17->_parameters.nelements)
  _Guide0_17->_parameters.nu = 0;
  #define nu (_Guide0_17->_parameters.nu)
  _Guide0_17->_parameters.phase = 0;
  #define phase (_Guide0_17->_parameters.phase)
  if("NULL" && strlen("NULL"))
    stracpy(_Guide0_17->_parameters.reflect, "NULL" ? "NULL" : "", 16384);
  else 
  _Guide0_17->_parameters.reflect[0]='\0';
  #define reflect (_Guide0_17->_parameters.reflect)

  #define GVars (_Guide0_17->_parameters.GVars)
  #define pTable (_Guide0_17->_parameters.pTable)

  #undef w1
  #undef h1
  #undef w2
  #undef h2
  #undef l
  #undef R0
  #undef Qc
  #undef alpha
  #undef m
  #undef W
  #undef nslit
  #undef d
  #undef mleft
  #undef mright
  #undef mtop
  #undef mbottom
  #undef nhslit
  #undef G
  #undef aleft
  #undef aright
  #undef atop
  #undef abottom
  #undef wavy
  #undef wavy_z
  #undef wavy_tb
  #undef wavy_lr
  #undef chamfers
  #undef chamfers_z
  #undef chamfers_lr
  #undef chamfers_tb
  #undef nelements
  #undef nu
  #undef phase
  #undef reflect
  #undef GVars
  #undef pTable
  /* component Guide0_17=Guide_gravity() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (2.0 / 2500 * RAD2DEG)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _Guide0_16->_rotation_absolute, _Guide0_17->_rotation_absolute);
    rot_transpose(_Guide0_16->_rotation_absolute, tr1);
    rot_mul(_Guide0_17->_rotation_absolute, tr1, _Guide0_17->_rotation_relative);
    _Guide0_17->_rotation_is_identity =  rot_test_identity(_Guide0_17->_rotation_relative);
    tc1 = coords_set(
      0, 0, 2 + 1e-4);
    rot_transpose(_Guide0_16->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Guide0_17->_position_absolute = coords_add(_Guide0_16->_position_absolute, tc2);
    tc1 = coords_sub(_Guide0_16->_position_absolute, _Guide0_17->_position_absolute);
    _Guide0_17->_position_relative = rot_apply(_Guide0_17->_rotation_absolute, tc1);
  } /* Guide0_17=Guide_gravity() AT ROTATED */
  DEBUG_COMPONENT("Guide0_17", _Guide0_17->_position_absolute, _Guide0_17->_rotation_absolute);
  instrument->_position_absolute[17] = _Guide0_17->_position_absolute;
  instrument->_position_relative[17] = _Guide0_17->_position_relative;
  instrument->counter_N[17]  = instrument->counter_P[17] = instrument->counter_P2[17] = 0;
  instrument->counter_AbsorbProp[17]= 0;
  return(0);
} /* _Guide0_17_setpos */

/* component Guide0_18=Guide_gravity() SETTING, POSITION/ROTATION */
int _Guide0_18_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Guide0_18_setpos] component Guide0_18=Guide_gravity() SETTING [/usr/share/mcstas/3.0-dev/optics/Guide_gravity.comp:339]");
  stracpy(_Guide0_18->_name, "Guide0_18", 16384);
  stracpy(_Guide0_18->_type, "Guide_gravity", 16384);
  _Guide0_18->_index=18;
  _Guide0_18->_parameters.w1 = 0.03;
  #define w1 (_Guide0_18->_parameters.w1)
  _Guide0_18->_parameters.h1 = 0.055;
  #define h1 (_Guide0_18->_parameters.h1)
  _Guide0_18->_parameters.w2 = 0;
  #define w2 (_Guide0_18->_parameters.w2)
  _Guide0_18->_parameters.h2 = 0;
  #define h2 (_Guide0_18->_parameters.h2)
  _Guide0_18->_parameters.l = 2;
  #define l (_Guide0_18->_parameters.l)
  _Guide0_18->_parameters.R0 = 0.995;
  #define R0 (_Guide0_18->_parameters.R0)
  _Guide0_18->_parameters.Qc = 0.0218;
  #define Qc (_Guide0_18->_parameters.Qc)
  _Guide0_18->_parameters.alpha = 4.38;
  #define alpha (_Guide0_18->_parameters.alpha)
  _Guide0_18->_parameters.m = 1.2;
  #define m (_Guide0_18->_parameters.m)
  _Guide0_18->_parameters.W = 0.003;
  #define W (_Guide0_18->_parameters.W)
  _Guide0_18->_parameters.nslit = 1;
  #define nslit (_Guide0_18->_parameters.nslit)
  _Guide0_18->_parameters.d = 0.0005;
  #define d (_Guide0_18->_parameters.d)
  _Guide0_18->_parameters.mleft = -1;
  #define mleft (_Guide0_18->_parameters.mleft)
  _Guide0_18->_parameters.mright = -1;
  #define mright (_Guide0_18->_parameters.mright)
  _Guide0_18->_parameters.mtop = -1;
  #define mtop (_Guide0_18->_parameters.mtop)
  _Guide0_18->_parameters.mbottom = -1;
  #define mbottom (_Guide0_18->_parameters.mbottom)
  _Guide0_18->_parameters.nhslit = 1;
  #define nhslit (_Guide0_18->_parameters.nhslit)
  _Guide0_18->_parameters.G = 0;
  #define G (_Guide0_18->_parameters.G)
  _Guide0_18->_parameters.aleft = -1;
  #define aleft (_Guide0_18->_parameters.aleft)
  _Guide0_18->_parameters.aright = -1;
  #define aright (_Guide0_18->_parameters.aright)
  _Guide0_18->_parameters.atop = -1;
  #define atop (_Guide0_18->_parameters.atop)
  _Guide0_18->_parameters.abottom = -1;
  #define abottom (_Guide0_18->_parameters.abottom)
  _Guide0_18->_parameters.wavy = 0;
  #define wavy (_Guide0_18->_parameters.wavy)
  _Guide0_18->_parameters.wavy_z = 0;
  #define wavy_z (_Guide0_18->_parameters.wavy_z)
  _Guide0_18->_parameters.wavy_tb = 0;
  #define wavy_tb (_Guide0_18->_parameters.wavy_tb)
  _Guide0_18->_parameters.wavy_lr = 0;
  #define wavy_lr (_Guide0_18->_parameters.wavy_lr)
  _Guide0_18->_parameters.chamfers = 0;
  #define chamfers (_Guide0_18->_parameters.chamfers)
  _Guide0_18->_parameters.chamfers_z = 0;
  #define chamfers_z (_Guide0_18->_parameters.chamfers_z)
  _Guide0_18->_parameters.chamfers_lr = 0;
  #define chamfers_lr (_Guide0_18->_parameters.chamfers_lr)
  _Guide0_18->_parameters.chamfers_tb = 0;
  #define chamfers_tb (_Guide0_18->_parameters.chamfers_tb)
  _Guide0_18->_parameters.nelements = 1;
  #define nelements (_Guide0_18->_parameters.nelements)
  _Guide0_18->_parameters.nu = 0;
  #define nu (_Guide0_18->_parameters.nu)
  _Guide0_18->_parameters.phase = 0;
  #define phase (_Guide0_18->_parameters.phase)
  if("NULL" && strlen("NULL"))
    stracpy(_Guide0_18->_parameters.reflect, "NULL" ? "NULL" : "", 16384);
  else 
  _Guide0_18->_parameters.reflect[0]='\0';
  #define reflect (_Guide0_18->_parameters.reflect)

  #define GVars (_Guide0_18->_parameters.GVars)
  #define pTable (_Guide0_18->_parameters.pTable)

  #undef w1
  #undef h1
  #undef w2
  #undef h2
  #undef l
  #undef R0
  #undef Qc
  #undef alpha
  #undef m
  #undef W
  #undef nslit
  #undef d
  #undef mleft
  #undef mright
  #undef mtop
  #undef mbottom
  #undef nhslit
  #undef G
  #undef aleft
  #undef aright
  #undef atop
  #undef abottom
  #undef wavy
  #undef wavy_z
  #undef wavy_tb
  #undef wavy_lr
  #undef chamfers
  #undef chamfers_z
  #undef chamfers_lr
  #undef chamfers_tb
  #undef nelements
  #undef nu
  #undef phase
  #undef reflect
  #undef GVars
  #undef pTable
  /* component Guide0_18=Guide_gravity() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (2.0 / 2500 * RAD2DEG)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _Guide0_17->_rotation_absolute, _Guide0_18->_rotation_absolute);
    rot_transpose(_Guide0_17->_rotation_absolute, tr1);
    rot_mul(_Guide0_18->_rotation_absolute, tr1, _Guide0_18->_rotation_relative);
    _Guide0_18->_rotation_is_identity =  rot_test_identity(_Guide0_18->_rotation_relative);
    tc1 = coords_set(
      0, 0, 2 + 1e-4);
    rot_transpose(_Guide0_17->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Guide0_18->_position_absolute = coords_add(_Guide0_17->_position_absolute, tc2);
    tc1 = coords_sub(_Guide0_17->_position_absolute, _Guide0_18->_position_absolute);
    _Guide0_18->_position_relative = rot_apply(_Guide0_18->_rotation_absolute, tc1);
  } /* Guide0_18=Guide_gravity() AT ROTATED */
  DEBUG_COMPONENT("Guide0_18", _Guide0_18->_position_absolute, _Guide0_18->_rotation_absolute);
  instrument->_position_absolute[18] = _Guide0_18->_position_absolute;
  instrument->_position_relative[18] = _Guide0_18->_position_relative;
  instrument->counter_N[18]  = instrument->counter_P[18] = instrument->counter_P2[18] = 0;
  instrument->counter_AbsorbProp[18]= 0;
  return(0);
} /* _Guide0_18_setpos */

/* component Guide0_19=Guide_gravity() SETTING, POSITION/ROTATION */
int _Guide0_19_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Guide0_19_setpos] component Guide0_19=Guide_gravity() SETTING [/usr/share/mcstas/3.0-dev/optics/Guide_gravity.comp:339]");
  stracpy(_Guide0_19->_name, "Guide0_19", 16384);
  stracpy(_Guide0_19->_type, "Guide_gravity", 16384);
  _Guide0_19->_index=19;
  _Guide0_19->_parameters.w1 = 0.03;
  #define w1 (_Guide0_19->_parameters.w1)
  _Guide0_19->_parameters.h1 = 0.055;
  #define h1 (_Guide0_19->_parameters.h1)
  _Guide0_19->_parameters.w2 = 0;
  #define w2 (_Guide0_19->_parameters.w2)
  _Guide0_19->_parameters.h2 = 0;
  #define h2 (_Guide0_19->_parameters.h2)
  _Guide0_19->_parameters.l = 2;
  #define l (_Guide0_19->_parameters.l)
  _Guide0_19->_parameters.R0 = 0.995;
  #define R0 (_Guide0_19->_parameters.R0)
  _Guide0_19->_parameters.Qc = 0.0218;
  #define Qc (_Guide0_19->_parameters.Qc)
  _Guide0_19->_parameters.alpha = 4.38;
  #define alpha (_Guide0_19->_parameters.alpha)
  _Guide0_19->_parameters.m = 1.2;
  #define m (_Guide0_19->_parameters.m)
  _Guide0_19->_parameters.W = 0.003;
  #define W (_Guide0_19->_parameters.W)
  _Guide0_19->_parameters.nslit = 1;
  #define nslit (_Guide0_19->_parameters.nslit)
  _Guide0_19->_parameters.d = 0.0005;
  #define d (_Guide0_19->_parameters.d)
  _Guide0_19->_parameters.mleft = -1;
  #define mleft (_Guide0_19->_parameters.mleft)
  _Guide0_19->_parameters.mright = -1;
  #define mright (_Guide0_19->_parameters.mright)
  _Guide0_19->_parameters.mtop = -1;
  #define mtop (_Guide0_19->_parameters.mtop)
  _Guide0_19->_parameters.mbottom = -1;
  #define mbottom (_Guide0_19->_parameters.mbottom)
  _Guide0_19->_parameters.nhslit = 1;
  #define nhslit (_Guide0_19->_parameters.nhslit)
  _Guide0_19->_parameters.G = 0;
  #define G (_Guide0_19->_parameters.G)
  _Guide0_19->_parameters.aleft = -1;
  #define aleft (_Guide0_19->_parameters.aleft)
  _Guide0_19->_parameters.aright = -1;
  #define aright (_Guide0_19->_parameters.aright)
  _Guide0_19->_parameters.atop = -1;
  #define atop (_Guide0_19->_parameters.atop)
  _Guide0_19->_parameters.abottom = -1;
  #define abottom (_Guide0_19->_parameters.abottom)
  _Guide0_19->_parameters.wavy = 0;
  #define wavy (_Guide0_19->_parameters.wavy)
  _Guide0_19->_parameters.wavy_z = 0;
  #define wavy_z (_Guide0_19->_parameters.wavy_z)
  _Guide0_19->_parameters.wavy_tb = 0;
  #define wavy_tb (_Guide0_19->_parameters.wavy_tb)
  _Guide0_19->_parameters.wavy_lr = 0;
  #define wavy_lr (_Guide0_19->_parameters.wavy_lr)
  _Guide0_19->_parameters.chamfers = 0;
  #define chamfers (_Guide0_19->_parameters.chamfers)
  _Guide0_19->_parameters.chamfers_z = 0;
  #define chamfers_z (_Guide0_19->_parameters.chamfers_z)
  _Guide0_19->_parameters.chamfers_lr = 0;
  #define chamfers_lr (_Guide0_19->_parameters.chamfers_lr)
  _Guide0_19->_parameters.chamfers_tb = 0;
  #define chamfers_tb (_Guide0_19->_parameters.chamfers_tb)
  _Guide0_19->_parameters.nelements = 1;
  #define nelements (_Guide0_19->_parameters.nelements)
  _Guide0_19->_parameters.nu = 0;
  #define nu (_Guide0_19->_parameters.nu)
  _Guide0_19->_parameters.phase = 0;
  #define phase (_Guide0_19->_parameters.phase)
  if("NULL" && strlen("NULL"))
    stracpy(_Guide0_19->_parameters.reflect, "NULL" ? "NULL" : "", 16384);
  else 
  _Guide0_19->_parameters.reflect[0]='\0';
  #define reflect (_Guide0_19->_parameters.reflect)

  #define GVars (_Guide0_19->_parameters.GVars)
  #define pTable (_Guide0_19->_parameters.pTable)

  #undef w1
  #undef h1
  #undef w2
  #undef h2
  #undef l
  #undef R0
  #undef Qc
  #undef alpha
  #undef m
  #undef W
  #undef nslit
  #undef d
  #undef mleft
  #undef mright
  #undef mtop
  #undef mbottom
  #undef nhslit
  #undef G
  #undef aleft
  #undef aright
  #undef atop
  #undef abottom
  #undef wavy
  #undef wavy_z
  #undef wavy_tb
  #undef wavy_lr
  #undef chamfers
  #undef chamfers_z
  #undef chamfers_lr
  #undef chamfers_tb
  #undef nelements
  #undef nu
  #undef phase
  #undef reflect
  #undef GVars
  #undef pTable
  /* component Guide0_19=Guide_gravity() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (2.0 / 2500 * RAD2DEG)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _Guide0_18->_rotation_absolute, _Guide0_19->_rotation_absolute);
    rot_transpose(_Guide0_18->_rotation_absolute, tr1);
    rot_mul(_Guide0_19->_rotation_absolute, tr1, _Guide0_19->_rotation_relative);
    _Guide0_19->_rotation_is_identity =  rot_test_identity(_Guide0_19->_rotation_relative);
    tc1 = coords_set(
      0, 0, 2 + 1e-4);
    rot_transpose(_Guide0_18->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Guide0_19->_position_absolute = coords_add(_Guide0_18->_position_absolute, tc2);
    tc1 = coords_sub(_Guide0_18->_position_absolute, _Guide0_19->_position_absolute);
    _Guide0_19->_position_relative = rot_apply(_Guide0_19->_rotation_absolute, tc1);
  } /* Guide0_19=Guide_gravity() AT ROTATED */
  DEBUG_COMPONENT("Guide0_19", _Guide0_19->_position_absolute, _Guide0_19->_rotation_absolute);
  instrument->_position_absolute[19] = _Guide0_19->_position_absolute;
  instrument->_position_relative[19] = _Guide0_19->_position_relative;
  instrument->counter_N[19]  = instrument->counter_P[19] = instrument->counter_P2[19] = 0;
  instrument->counter_AbsorbProp[19]= 0;
  return(0);
} /* _Guide0_19_setpos */

/* component Guide_PSD=Monitor_nD() SETTING, POSITION/ROTATION */
int _Guide_PSD_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Guide_PSD_setpos] component Guide_PSD=Monitor_nD() SETTING [/usr/share/mcstas/3.0-dev/monitors/Monitor_nD.comp:228]");
  stracpy(_Guide_PSD->_name, "Guide_PSD", 16384);
  stracpy(_Guide_PSD->_type, "Monitor_nD", 16384);
  _Guide_PSD->_index=20;
  _Guide_PSD->_parameters.user1 = '\0';
  #define user1 (_Guide_PSD->_parameters.user1)
  _Guide_PSD->_parameters.user2 = '\0';
  #define user2 (_Guide_PSD->_parameters.user2)
  _Guide_PSD->_parameters.user3 = '\0';
  #define user3 (_Guide_PSD->_parameters.user3)
  _Guide_PSD->_parameters.xwidth = 0.03;
  #define xwidth (_Guide_PSD->_parameters.xwidth)
  _Guide_PSD->_parameters.yheight = 0.055;
  #define yheight (_Guide_PSD->_parameters.yheight)
  _Guide_PSD->_parameters.zdepth = 0;
  #define zdepth (_Guide_PSD->_parameters.zdepth)
  _Guide_PSD->_parameters.xmin = 0;
  #define xmin (_Guide_PSD->_parameters.xmin)
  _Guide_PSD->_parameters.xmax = 0;
  #define xmax (_Guide_PSD->_parameters.xmax)
  _Guide_PSD->_parameters.ymin = 0;
  #define ymin (_Guide_PSD->_parameters.ymin)
  _Guide_PSD->_parameters.ymax = 0;
  #define ymax (_Guide_PSD->_parameters.ymax)
  _Guide_PSD->_parameters.zmin = 0;
  #define zmin (_Guide_PSD->_parameters.zmin)
  _Guide_PSD->_parameters.zmax = 0;
  #define zmax (_Guide_PSD->_parameters.zmax)
  _Guide_PSD->_parameters.bins = 0;
  #define bins (_Guide_PSD->_parameters.bins)
  _Guide_PSD->_parameters.min = -1e40;
  #define min (_Guide_PSD->_parameters.min)
  _Guide_PSD->_parameters.max = 1e40;
  #define max (_Guide_PSD->_parameters.max)
  _Guide_PSD->_parameters.restore_neutron = 0;
  #define restore_neutron (_Guide_PSD->_parameters.restore_neutron)
  _Guide_PSD->_parameters.radius = 0;
  #define radius (_Guide_PSD->_parameters.radius)
  if("y x" && strlen("y x"))
    stracpy(_Guide_PSD->_parameters.options, "y x" ? "y x" : "", 16384);
  else 
  _Guide_PSD->_parameters.options[0]='\0';
  #define options (_Guide_PSD->_parameters.options)
  if("NULL" && strlen("NULL"))
    stracpy(_Guide_PSD->_parameters.filename, "NULL" ? "NULL" : "", 16384);
  else 
  _Guide_PSD->_parameters.filename[0]='\0';
  #define filename (_Guide_PSD->_parameters.filename)
  if("NULL" && strlen("NULL"))
    stracpy(_Guide_PSD->_parameters.geometry, "NULL" ? "NULL" : "", 16384);
  else 
  _Guide_PSD->_parameters.geometry[0]='\0';
  #define geometry (_Guide_PSD->_parameters.geometry)
  if("NULL" && strlen("NULL"))
    stracpy(_Guide_PSD->_parameters.username1, "NULL" ? "NULL" : "", 16384);
  else 
  _Guide_PSD->_parameters.username1[0]='\0';
  #define username1 (_Guide_PSD->_parameters.username1)
  if("NULL" && strlen("NULL"))
    stracpy(_Guide_PSD->_parameters.username2, "NULL" ? "NULL" : "", 16384);
  else 
  _Guide_PSD->_parameters.username2[0]='\0';
  #define username2 (_Guide_PSD->_parameters.username2)
  if("NULL" && strlen("NULL"))
    stracpy(_Guide_PSD->_parameters.username3, "NULL" ? "NULL" : "", 16384);
  else 
  _Guide_PSD->_parameters.username3[0]='\0';
  #define username3 (_Guide_PSD->_parameters.username3)

  #define DEFS (_Guide_PSD->_parameters.DEFS)
  #define Vars (_Guide_PSD->_parameters.Vars)
  #define detector (_Guide_PSD->_parameters.detector)
  #define offdata (_Guide_PSD->_parameters.offdata)

  #undef user1
  #undef user2
  #undef user3
  #undef xwidth
  #undef yheight
  #undef zdepth
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef zmin
  #undef zmax
  #undef bins
  #undef min
  #undef max
  #undef restore_neutron
  #undef radius
  #undef options
  #undef filename
  #undef geometry
  #undef username1
  #undef username2
  #undef username3
  #undef DEFS
  #undef Vars
  #undef detector
  #undef offdata
  /* component Guide_PSD=Monitor_nD() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Guide0_19->_rotation_absolute, _Guide_PSD->_rotation_absolute);
    rot_transpose(_Guide0_19->_rotation_absolute, tr1);
    rot_mul(_Guide_PSD->_rotation_absolute, tr1, _Guide_PSD->_rotation_relative);
    _Guide_PSD->_rotation_is_identity =  rot_test_identity(_Guide_PSD->_rotation_relative);
    tc1 = coords_set(
      0, 0, 2 + 1e-4);
    rot_transpose(_Guide0_19->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Guide_PSD->_position_absolute = coords_add(_Guide0_19->_position_absolute, tc2);
    tc1 = coords_sub(_Guide0_19->_position_absolute, _Guide_PSD->_position_absolute);
    _Guide_PSD->_position_relative = rot_apply(_Guide_PSD->_rotation_absolute, tc1);
  } /* Guide_PSD=Monitor_nD() AT ROTATED */
  DEBUG_COMPONENT("Guide_PSD", _Guide_PSD->_position_absolute, _Guide_PSD->_rotation_absolute);
  instrument->_position_absolute[20] = _Guide_PSD->_position_absolute;
  instrument->_position_relative[20] = _Guide_PSD->_position_relative;
  instrument->counter_N[20]  = instrument->counter_P[20] = instrument->counter_P2[20] = 0;
  instrument->counter_AbsorbProp[20]= 0;
  return(0);
} /* _Guide_PSD_setpos */

/* component Guide_Lambda=Monitor_nD() SETTING, POSITION/ROTATION */
int _Guide_Lambda_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Guide_Lambda_setpos] component Guide_Lambda=Monitor_nD() SETTING [/usr/share/mcstas/3.0-dev/monitors/Monitor_nD.comp:228]");
  stracpy(_Guide_Lambda->_name, "Guide_Lambda", 16384);
  stracpy(_Guide_Lambda->_type, "Monitor_nD", 16384);
  _Guide_Lambda->_index=21;
  _Guide_Lambda->_parameters.user1 = '\0';
  #define user1 (_Guide_Lambda->_parameters.user1)
  _Guide_Lambda->_parameters.user2 = '\0';
  #define user2 (_Guide_Lambda->_parameters.user2)
  _Guide_Lambda->_parameters.user3 = '\0';
  #define user3 (_Guide_Lambda->_parameters.user3)
  _Guide_Lambda->_parameters.xwidth = 0.03;
  #define xwidth (_Guide_Lambda->_parameters.xwidth)
  _Guide_Lambda->_parameters.yheight = 0.055;
  #define yheight (_Guide_Lambda->_parameters.yheight)
  _Guide_Lambda->_parameters.zdepth = 0;
  #define zdepth (_Guide_Lambda->_parameters.zdepth)
  _Guide_Lambda->_parameters.xmin = 0;
  #define xmin (_Guide_Lambda->_parameters.xmin)
  _Guide_Lambda->_parameters.xmax = 0;
  #define xmax (_Guide_Lambda->_parameters.xmax)
  _Guide_Lambda->_parameters.ymin = 0;
  #define ymin (_Guide_Lambda->_parameters.ymin)
  _Guide_Lambda->_parameters.ymax = 0;
  #define ymax (_Guide_Lambda->_parameters.ymax)
  _Guide_Lambda->_parameters.zmin = 0;
  #define zmin (_Guide_Lambda->_parameters.zmin)
  _Guide_Lambda->_parameters.zmax = 0;
  #define zmax (_Guide_Lambda->_parameters.zmax)
  _Guide_Lambda->_parameters.bins = 0;
  #define bins (_Guide_Lambda->_parameters.bins)
  _Guide_Lambda->_parameters.min = -1e40;
  #define min (_Guide_Lambda->_parameters.min)
  _Guide_Lambda->_parameters.max = 1e40;
  #define max (_Guide_Lambda->_parameters.max)
  _Guide_Lambda->_parameters.restore_neutron = 0;
  #define restore_neutron (_Guide_Lambda->_parameters.restore_neutron)
  _Guide_Lambda->_parameters.radius = 0;
  #define radius (_Guide_Lambda->_parameters.radius)
  if("lambda, auto" && strlen("lambda, auto"))
    stracpy(_Guide_Lambda->_parameters.options, "lambda, auto" ? "lambda, auto" : "", 16384);
  else 
  _Guide_Lambda->_parameters.options[0]='\0';
  #define options (_Guide_Lambda->_parameters.options)
  if("NULL" && strlen("NULL"))
    stracpy(_Guide_Lambda->_parameters.filename, "NULL" ? "NULL" : "", 16384);
  else 
  _Guide_Lambda->_parameters.filename[0]='\0';
  #define filename (_Guide_Lambda->_parameters.filename)
  if("NULL" && strlen("NULL"))
    stracpy(_Guide_Lambda->_parameters.geometry, "NULL" ? "NULL" : "", 16384);
  else 
  _Guide_Lambda->_parameters.geometry[0]='\0';
  #define geometry (_Guide_Lambda->_parameters.geometry)
  if("NULL" && strlen("NULL"))
    stracpy(_Guide_Lambda->_parameters.username1, "NULL" ? "NULL" : "", 16384);
  else 
  _Guide_Lambda->_parameters.username1[0]='\0';
  #define username1 (_Guide_Lambda->_parameters.username1)
  if("NULL" && strlen("NULL"))
    stracpy(_Guide_Lambda->_parameters.username2, "NULL" ? "NULL" : "", 16384);
  else 
  _Guide_Lambda->_parameters.username2[0]='\0';
  #define username2 (_Guide_Lambda->_parameters.username2)
  if("NULL" && strlen("NULL"))
    stracpy(_Guide_Lambda->_parameters.username3, "NULL" ? "NULL" : "", 16384);
  else 
  _Guide_Lambda->_parameters.username3[0]='\0';
  #define username3 (_Guide_Lambda->_parameters.username3)

  #define DEFS (_Guide_Lambda->_parameters.DEFS)
  #define Vars (_Guide_Lambda->_parameters.Vars)
  #define detector (_Guide_Lambda->_parameters.detector)
  #define offdata (_Guide_Lambda->_parameters.offdata)

  #undef user1
  #undef user2
  #undef user3
  #undef xwidth
  #undef yheight
  #undef zdepth
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef zmin
  #undef zmax
  #undef bins
  #undef min
  #undef max
  #undef restore_neutron
  #undef radius
  #undef options
  #undef filename
  #undef geometry
  #undef username1
  #undef username2
  #undef username3
  #undef DEFS
  #undef Vars
  #undef detector
  #undef offdata
  /* component Guide_Lambda=Monitor_nD() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Guide_PSD->_rotation_absolute, _Guide_Lambda->_rotation_absolute);
    rot_transpose(_Guide_PSD->_rotation_absolute, tr1);
    rot_mul(_Guide_Lambda->_rotation_absolute, tr1, _Guide_Lambda->_rotation_relative);
    _Guide_Lambda->_rotation_is_identity =  rot_test_identity(_Guide_Lambda->_rotation_relative);
    tc1 = coords_set(
      0, 0, 1e-4);
    rot_transpose(_Guide_PSD->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Guide_Lambda->_position_absolute = coords_add(_Guide_PSD->_position_absolute, tc2);
    tc1 = coords_sub(_Guide_PSD->_position_absolute, _Guide_Lambda->_position_absolute);
    _Guide_Lambda->_position_relative = rot_apply(_Guide_Lambda->_rotation_absolute, tc1);
  } /* Guide_Lambda=Monitor_nD() AT ROTATED */
  DEBUG_COMPONENT("Guide_Lambda", _Guide_Lambda->_position_absolute, _Guide_Lambda->_rotation_absolute);
  instrument->_position_absolute[21] = _Guide_Lambda->_position_absolute;
  instrument->_position_relative[21] = _Guide_Lambda->_position_relative;
  instrument->counter_N[21]  = instrument->counter_P[21] = instrument->counter_P2[21] = 0;
  instrument->counter_AbsorbProp[21]= 0;
  return(0);
} /* _Guide_Lambda_setpos */

/* component Chopper1=DiskChopper() SETTING, POSITION/ROTATION */
int _Chopper1_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Chopper1_setpos] component Chopper1=DiskChopper() SETTING [/usr/share/mcstas/3.0-dev/optics/DiskChopper.comp:67]");
  stracpy(_Chopper1->_name, "Chopper1", 16384);
  stracpy(_Chopper1->_type, "DiskChopper", 16384);
  _Chopper1->_index=22;
  _Chopper1->_parameters.theta_0 = 0;
  #define theta_0 (_Chopper1->_parameters.theta_0)
  _Chopper1->_parameters.radius = 0.24;
  #define radius (_Chopper1->_parameters.radius)
  _Chopper1->_parameters.yheight = 0.06;
  #define yheight (_Chopper1->_parameters.yheight)
  _Chopper1->_parameters.nu = instrument->_parameters._rpm / 60;
  #define nu (_Chopper1->_parameters.nu)
  _Chopper1->_parameters.nslit = 3;
  #define nslit (_Chopper1->_parameters.nslit)
  _Chopper1->_parameters.jitter = 0;
  #define jitter (_Chopper1->_parameters.jitter)
  _Chopper1->_parameters.delay = 0;
  #define delay (_Chopper1->_parameters.delay)
  _Chopper1->_parameters.isfirst = 1;
  #define isfirst (_Chopper1->_parameters.isfirst)
  _Chopper1->_parameters.n_pulse = 1;
  #define n_pulse (_Chopper1->_parameters.n_pulse)
  _Chopper1->_parameters.abs_out = 1;
  #define abs_out (_Chopper1->_parameters.abs_out)
  _Chopper1->_parameters.phase = 0;
  #define phase (_Chopper1->_parameters.phase)
  _Chopper1->_parameters.xwidth = 0.03;
  #define xwidth (_Chopper1->_parameters.xwidth)
  _Chopper1->_parameters.verbose = 0;
  #define verbose (_Chopper1->_parameters.verbose)

  #define Tg (_Chopper1->_parameters.Tg)
  #define To (_Chopper1->_parameters.To)
  #define delta_y (_Chopper1->_parameters.delta_y)
  #define height (_Chopper1->_parameters.height)
  #define omega (_Chopper1->_parameters.omega)

  #undef theta_0
  #undef radius
  #undef yheight
  #undef nu
  #undef nslit
  #undef jitter
  #undef delay
  #undef isfirst
  #undef n_pulse
  #undef abs_out
  #undef phase
  #undef xwidth
  #undef verbose
  #undef Tg
  #undef To
  #undef delta_y
  #undef height
  #undef omega
  /* component Chopper1=DiskChopper() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Guide_Lambda->_rotation_absolute, _Chopper1->_rotation_absolute);
    rot_transpose(_Guide_Lambda->_rotation_absolute, tr1);
    rot_mul(_Chopper1->_rotation_absolute, tr1, _Chopper1->_rotation_relative);
    _Chopper1->_rotation_is_identity =  rot_test_identity(_Chopper1->_rotation_relative);
    tc1 = coords_set(
      0, 0, 1.2e-2 / 2);
    rot_transpose(_Guide_Lambda->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Chopper1->_position_absolute = coords_add(_Guide_Lambda->_position_absolute, tc2);
    tc1 = coords_sub(_Guide_Lambda->_position_absolute, _Chopper1->_position_absolute);
    _Chopper1->_position_relative = rot_apply(_Chopper1->_rotation_absolute, tc1);
  } /* Chopper1=DiskChopper() AT ROTATED */
  DEBUG_COMPONENT("Chopper1", _Chopper1->_position_absolute, _Chopper1->_rotation_absolute);
  instrument->_position_absolute[22] = _Chopper1->_position_absolute;
  instrument->_position_relative[22] = _Chopper1->_position_relative;
  instrument->counter_N[22]  = instrument->counter_P[22] = instrument->counter_P2[22] = 0;
  instrument->counter_AbsorbProp[22]= 0;
  return(0);
} /* _Chopper1_setpos */

/* component Guide_SGS1=Guide_gravity() SETTING, POSITION/ROTATION */
int _Guide_SGS1_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Guide_SGS1_setpos] component Guide_SGS1=Guide_gravity() SETTING [/usr/share/mcstas/3.0-dev/optics/Guide_gravity.comp:339]");
  stracpy(_Guide_SGS1->_name, "Guide_SGS1", 16384);
  stracpy(_Guide_SGS1->_type, "Guide_gravity", 16384);
  _Guide_SGS1->_index=23;
  _Guide_SGS1->_parameters.w1 = 0.03;
  #define w1 (_Guide_SGS1->_parameters.w1)
  _Guide_SGS1->_parameters.h1 = 0.055;
  #define h1 (_Guide_SGS1->_parameters.h1)
  _Guide_SGS1->_parameters.w2 = 0;
  #define w2 (_Guide_SGS1->_parameters.w2)
  _Guide_SGS1->_parameters.h2 = 0;
  #define h2 (_Guide_SGS1->_parameters.h2)
  _Guide_SGS1->_parameters.l = 2.41 -2 * 6e-3;
  #define l (_Guide_SGS1->_parameters.l)
  _Guide_SGS1->_parameters.R0 = 0.995;
  #define R0 (_Guide_SGS1->_parameters.R0)
  _Guide_SGS1->_parameters.Qc = 0.0218;
  #define Qc (_Guide_SGS1->_parameters.Qc)
  _Guide_SGS1->_parameters.alpha = 4.38;
  #define alpha (_Guide_SGS1->_parameters.alpha)
  _Guide_SGS1->_parameters.m = 1.2;
  #define m (_Guide_SGS1->_parameters.m)
  _Guide_SGS1->_parameters.W = 0.003;
  #define W (_Guide_SGS1->_parameters.W)
  _Guide_SGS1->_parameters.nslit = 1;
  #define nslit (_Guide_SGS1->_parameters.nslit)
  _Guide_SGS1->_parameters.d = 0.0005;
  #define d (_Guide_SGS1->_parameters.d)
  _Guide_SGS1->_parameters.mleft = -1;
  #define mleft (_Guide_SGS1->_parameters.mleft)
  _Guide_SGS1->_parameters.mright = -1;
  #define mright (_Guide_SGS1->_parameters.mright)
  _Guide_SGS1->_parameters.mtop = -1;
  #define mtop (_Guide_SGS1->_parameters.mtop)
  _Guide_SGS1->_parameters.mbottom = -1;
  #define mbottom (_Guide_SGS1->_parameters.mbottom)
  _Guide_SGS1->_parameters.nhslit = 1;
  #define nhslit (_Guide_SGS1->_parameters.nhslit)
  _Guide_SGS1->_parameters.G = 0;
  #define G (_Guide_SGS1->_parameters.G)
  _Guide_SGS1->_parameters.aleft = -1;
  #define aleft (_Guide_SGS1->_parameters.aleft)
  _Guide_SGS1->_parameters.aright = -1;
  #define aright (_Guide_SGS1->_parameters.aright)
  _Guide_SGS1->_parameters.atop = -1;
  #define atop (_Guide_SGS1->_parameters.atop)
  _Guide_SGS1->_parameters.abottom = -1;
  #define abottom (_Guide_SGS1->_parameters.abottom)
  _Guide_SGS1->_parameters.wavy = 0;
  #define wavy (_Guide_SGS1->_parameters.wavy)
  _Guide_SGS1->_parameters.wavy_z = 0;
  #define wavy_z (_Guide_SGS1->_parameters.wavy_z)
  _Guide_SGS1->_parameters.wavy_tb = 0;
  #define wavy_tb (_Guide_SGS1->_parameters.wavy_tb)
  _Guide_SGS1->_parameters.wavy_lr = 0;
  #define wavy_lr (_Guide_SGS1->_parameters.wavy_lr)
  _Guide_SGS1->_parameters.chamfers = 0;
  #define chamfers (_Guide_SGS1->_parameters.chamfers)
  _Guide_SGS1->_parameters.chamfers_z = 0;
  #define chamfers_z (_Guide_SGS1->_parameters.chamfers_z)
  _Guide_SGS1->_parameters.chamfers_lr = 0;
  #define chamfers_lr (_Guide_SGS1->_parameters.chamfers_lr)
  _Guide_SGS1->_parameters.chamfers_tb = 0;
  #define chamfers_tb (_Guide_SGS1->_parameters.chamfers_tb)
  _Guide_SGS1->_parameters.nelements = 1;
  #define nelements (_Guide_SGS1->_parameters.nelements)
  _Guide_SGS1->_parameters.nu = 0;
  #define nu (_Guide_SGS1->_parameters.nu)
  _Guide_SGS1->_parameters.phase = 0;
  #define phase (_Guide_SGS1->_parameters.phase)
  if("NULL" && strlen("NULL"))
    stracpy(_Guide_SGS1->_parameters.reflect, "NULL" ? "NULL" : "", 16384);
  else 
  _Guide_SGS1->_parameters.reflect[0]='\0';
  #define reflect (_Guide_SGS1->_parameters.reflect)

  #define GVars (_Guide_SGS1->_parameters.GVars)
  #define pTable (_Guide_SGS1->_parameters.pTable)

  #undef w1
  #undef h1
  #undef w2
  #undef h2
  #undef l
  #undef R0
  #undef Qc
  #undef alpha
  #undef m
  #undef W
  #undef nslit
  #undef d
  #undef mleft
  #undef mright
  #undef mtop
  #undef mbottom
  #undef nhslit
  #undef G
  #undef aleft
  #undef aright
  #undef atop
  #undef abottom
  #undef wavy
  #undef wavy_z
  #undef wavy_tb
  #undef wavy_lr
  #undef chamfers
  #undef chamfers_z
  #undef chamfers_lr
  #undef chamfers_tb
  #undef nelements
  #undef nu
  #undef phase
  #undef reflect
  #undef GVars
  #undef pTable
  /* component Guide_SGS1=Guide_gravity() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Guide_Lambda->_rotation_absolute, _Guide_SGS1->_rotation_absolute);
    rot_transpose(_Chopper1->_rotation_absolute, tr1);
    rot_mul(_Guide_SGS1->_rotation_absolute, tr1, _Guide_SGS1->_rotation_relative);
    _Guide_SGS1->_rotation_is_identity =  rot_test_identity(_Guide_SGS1->_rotation_relative);
    tc1 = coords_set(
      0, 0, 1.2e-2);
    rot_transpose(_Guide_Lambda->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Guide_SGS1->_position_absolute = coords_add(_Guide_Lambda->_position_absolute, tc2);
    tc1 = coords_sub(_Chopper1->_position_absolute, _Guide_SGS1->_position_absolute);
    _Guide_SGS1->_position_relative = rot_apply(_Guide_SGS1->_rotation_absolute, tc1);
  } /* Guide_SGS1=Guide_gravity() AT ROTATED */
  DEBUG_COMPONENT("Guide_SGS1", _Guide_SGS1->_position_absolute, _Guide_SGS1->_rotation_absolute);
  instrument->_position_absolute[23] = _Guide_SGS1->_position_absolute;
  instrument->_position_relative[23] = _Guide_SGS1->_position_relative;
  instrument->counter_N[23]  = instrument->counter_P[23] = instrument->counter_P2[23] = 0;
  instrument->counter_AbsorbProp[23]= 0;
  return(0);
} /* _Guide_SGS1_setpos */

/* component Guide_SGS2=Guide_gravity() SETTING, POSITION/ROTATION */
int _Guide_SGS2_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Guide_SGS2_setpos] component Guide_SGS2=Guide_gravity() SETTING [/usr/share/mcstas/3.0-dev/optics/Guide_gravity.comp:339]");
  stracpy(_Guide_SGS2->_name, "Guide_SGS2", 16384);
  stracpy(_Guide_SGS2->_type, "Guide_gravity", 16384);
  _Guide_SGS2->_index=24;
  _Guide_SGS2->_parameters.w1 = 0.03;
  #define w1 (_Guide_SGS2->_parameters.w1)
  _Guide_SGS2->_parameters.h1 = 0.055;
  #define h1 (_Guide_SGS2->_parameters.h1)
  _Guide_SGS2->_parameters.w2 = 0;
  #define w2 (_Guide_SGS2->_parameters.w2)
  _Guide_SGS2->_parameters.h2 = 0;
  #define h2 (_Guide_SGS2->_parameters.h2)
  _Guide_SGS2->_parameters.l = 5.065 -2 * 6e-3;
  #define l (_Guide_SGS2->_parameters.l)
  _Guide_SGS2->_parameters.R0 = 0.995;
  #define R0 (_Guide_SGS2->_parameters.R0)
  _Guide_SGS2->_parameters.Qc = 0.0218;
  #define Qc (_Guide_SGS2->_parameters.Qc)
  _Guide_SGS2->_parameters.alpha = 4.38;
  #define alpha (_Guide_SGS2->_parameters.alpha)
  _Guide_SGS2->_parameters.m = 1.2;
  #define m (_Guide_SGS2->_parameters.m)
  _Guide_SGS2->_parameters.W = 0.003;
  #define W (_Guide_SGS2->_parameters.W)
  _Guide_SGS2->_parameters.nslit = 1;
  #define nslit (_Guide_SGS2->_parameters.nslit)
  _Guide_SGS2->_parameters.d = 0.0005;
  #define d (_Guide_SGS2->_parameters.d)
  _Guide_SGS2->_parameters.mleft = -1;
  #define mleft (_Guide_SGS2->_parameters.mleft)
  _Guide_SGS2->_parameters.mright = -1;
  #define mright (_Guide_SGS2->_parameters.mright)
  _Guide_SGS2->_parameters.mtop = -1;
  #define mtop (_Guide_SGS2->_parameters.mtop)
  _Guide_SGS2->_parameters.mbottom = -1;
  #define mbottom (_Guide_SGS2->_parameters.mbottom)
  _Guide_SGS2->_parameters.nhslit = 1;
  #define nhslit (_Guide_SGS2->_parameters.nhslit)
  _Guide_SGS2->_parameters.G = 0;
  #define G (_Guide_SGS2->_parameters.G)
  _Guide_SGS2->_parameters.aleft = -1;
  #define aleft (_Guide_SGS2->_parameters.aleft)
  _Guide_SGS2->_parameters.aright = -1;
  #define aright (_Guide_SGS2->_parameters.aright)
  _Guide_SGS2->_parameters.atop = -1;
  #define atop (_Guide_SGS2->_parameters.atop)
  _Guide_SGS2->_parameters.abottom = -1;
  #define abottom (_Guide_SGS2->_parameters.abottom)
  _Guide_SGS2->_parameters.wavy = 0;
  #define wavy (_Guide_SGS2->_parameters.wavy)
  _Guide_SGS2->_parameters.wavy_z = 0;
  #define wavy_z (_Guide_SGS2->_parameters.wavy_z)
  _Guide_SGS2->_parameters.wavy_tb = 0;
  #define wavy_tb (_Guide_SGS2->_parameters.wavy_tb)
  _Guide_SGS2->_parameters.wavy_lr = 0;
  #define wavy_lr (_Guide_SGS2->_parameters.wavy_lr)
  _Guide_SGS2->_parameters.chamfers = 0;
  #define chamfers (_Guide_SGS2->_parameters.chamfers)
  _Guide_SGS2->_parameters.chamfers_z = 0;
  #define chamfers_z (_Guide_SGS2->_parameters.chamfers_z)
  _Guide_SGS2->_parameters.chamfers_lr = 0;
  #define chamfers_lr (_Guide_SGS2->_parameters.chamfers_lr)
  _Guide_SGS2->_parameters.chamfers_tb = 0;
  #define chamfers_tb (_Guide_SGS2->_parameters.chamfers_tb)
  _Guide_SGS2->_parameters.nelements = 1;
  #define nelements (_Guide_SGS2->_parameters.nelements)
  _Guide_SGS2->_parameters.nu = 0;
  #define nu (_Guide_SGS2->_parameters.nu)
  _Guide_SGS2->_parameters.phase = 0;
  #define phase (_Guide_SGS2->_parameters.phase)
  if("NULL" && strlen("NULL"))
    stracpy(_Guide_SGS2->_parameters.reflect, "NULL" ? "NULL" : "", 16384);
  else 
  _Guide_SGS2->_parameters.reflect[0]='\0';
  #define reflect (_Guide_SGS2->_parameters.reflect)

  #define GVars (_Guide_SGS2->_parameters.GVars)
  #define pTable (_Guide_SGS2->_parameters.pTable)

  #undef w1
  #undef h1
  #undef w2
  #undef h2
  #undef l
  #undef R0
  #undef Qc
  #undef alpha
  #undef m
  #undef W
  #undef nslit
  #undef d
  #undef mleft
  #undef mright
  #undef mtop
  #undef mbottom
  #undef nhslit
  #undef G
  #undef aleft
  #undef aright
  #undef atop
  #undef abottom
  #undef wavy
  #undef wavy_z
  #undef wavy_tb
  #undef wavy_lr
  #undef chamfers
  #undef chamfers_z
  #undef chamfers_lr
  #undef chamfers_tb
  #undef nelements
  #undef nu
  #undef phase
  #undef reflect
  #undef GVars
  #undef pTable
  /* component Guide_SGS2=Guide_gravity() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Guide_SGS1->_rotation_absolute, _Guide_SGS2->_rotation_absolute);
    rot_transpose(_Guide_SGS1->_rotation_absolute, tr1);
    rot_mul(_Guide_SGS2->_rotation_absolute, tr1, _Guide_SGS2->_rotation_relative);
    _Guide_SGS2->_rotation_is_identity =  rot_test_identity(_Guide_SGS2->_rotation_relative);
    tc1 = coords_set(
      0, 0, 2.41);
    rot_transpose(_Guide_SGS1->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Guide_SGS2->_position_absolute = coords_add(_Guide_SGS1->_position_absolute, tc2);
    tc1 = coords_sub(_Guide_SGS1->_position_absolute, _Guide_SGS2->_position_absolute);
    _Guide_SGS2->_position_relative = rot_apply(_Guide_SGS2->_rotation_absolute, tc1);
  } /* Guide_SGS2=Guide_gravity() AT ROTATED */
  DEBUG_COMPONENT("Guide_SGS2", _Guide_SGS2->_position_absolute, _Guide_SGS2->_rotation_absolute);
  instrument->_position_absolute[24] = _Guide_SGS2->_position_absolute;
  instrument->_position_relative[24] = _Guide_SGS2->_position_relative;
  instrument->counter_N[24]  = instrument->counter_P[24] = instrument->counter_P2[24] = 0;
  instrument->counter_AbsorbProp[24]= 0;
  return(0);
} /* _Guide_SGS2_setpos */

/* component Guide_SGS3=Guide_gravity() SETTING, POSITION/ROTATION */
int _Guide_SGS3_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Guide_SGS3_setpos] component Guide_SGS3=Guide_gravity() SETTING [/usr/share/mcstas/3.0-dev/optics/Guide_gravity.comp:339]");
  stracpy(_Guide_SGS3->_name, "Guide_SGS3", 16384);
  stracpy(_Guide_SGS3->_type, "Guide_gravity", 16384);
  _Guide_SGS3->_index=25;
  _Guide_SGS3->_parameters.w1 = 0.03;
  #define w1 (_Guide_SGS3->_parameters.w1)
  _Guide_SGS3->_parameters.h1 = 0.055;
  #define h1 (_Guide_SGS3->_parameters.h1)
  _Guide_SGS3->_parameters.w2 = 0;
  #define w2 (_Guide_SGS3->_parameters.w2)
  _Guide_SGS3->_parameters.h2 = 0;
  #define h2 (_Guide_SGS3->_parameters.h2)
  _Guide_SGS3->_parameters.l = 2.352 -2 * 6e-3;
  #define l (_Guide_SGS3->_parameters.l)
  _Guide_SGS3->_parameters.R0 = 0.995;
  #define R0 (_Guide_SGS3->_parameters.R0)
  _Guide_SGS3->_parameters.Qc = 0.0218;
  #define Qc (_Guide_SGS3->_parameters.Qc)
  _Guide_SGS3->_parameters.alpha = 4.38;
  #define alpha (_Guide_SGS3->_parameters.alpha)
  _Guide_SGS3->_parameters.m = 1.2;
  #define m (_Guide_SGS3->_parameters.m)
  _Guide_SGS3->_parameters.W = 0.003;
  #define W (_Guide_SGS3->_parameters.W)
  _Guide_SGS3->_parameters.nslit = 1;
  #define nslit (_Guide_SGS3->_parameters.nslit)
  _Guide_SGS3->_parameters.d = 0.0005;
  #define d (_Guide_SGS3->_parameters.d)
  _Guide_SGS3->_parameters.mleft = -1;
  #define mleft (_Guide_SGS3->_parameters.mleft)
  _Guide_SGS3->_parameters.mright = -1;
  #define mright (_Guide_SGS3->_parameters.mright)
  _Guide_SGS3->_parameters.mtop = -1;
  #define mtop (_Guide_SGS3->_parameters.mtop)
  _Guide_SGS3->_parameters.mbottom = -1;
  #define mbottom (_Guide_SGS3->_parameters.mbottom)
  _Guide_SGS3->_parameters.nhslit = 1;
  #define nhslit (_Guide_SGS3->_parameters.nhslit)
  _Guide_SGS3->_parameters.G = 0;
  #define G (_Guide_SGS3->_parameters.G)
  _Guide_SGS3->_parameters.aleft = -1;
  #define aleft (_Guide_SGS3->_parameters.aleft)
  _Guide_SGS3->_parameters.aright = -1;
  #define aright (_Guide_SGS3->_parameters.aright)
  _Guide_SGS3->_parameters.atop = -1;
  #define atop (_Guide_SGS3->_parameters.atop)
  _Guide_SGS3->_parameters.abottom = -1;
  #define abottom (_Guide_SGS3->_parameters.abottom)
  _Guide_SGS3->_parameters.wavy = 0;
  #define wavy (_Guide_SGS3->_parameters.wavy)
  _Guide_SGS3->_parameters.wavy_z = 0;
  #define wavy_z (_Guide_SGS3->_parameters.wavy_z)
  _Guide_SGS3->_parameters.wavy_tb = 0;
  #define wavy_tb (_Guide_SGS3->_parameters.wavy_tb)
  _Guide_SGS3->_parameters.wavy_lr = 0;
  #define wavy_lr (_Guide_SGS3->_parameters.wavy_lr)
  _Guide_SGS3->_parameters.chamfers = 0;
  #define chamfers (_Guide_SGS3->_parameters.chamfers)
  _Guide_SGS3->_parameters.chamfers_z = 0;
  #define chamfers_z (_Guide_SGS3->_parameters.chamfers_z)
  _Guide_SGS3->_parameters.chamfers_lr = 0;
  #define chamfers_lr (_Guide_SGS3->_parameters.chamfers_lr)
  _Guide_SGS3->_parameters.chamfers_tb = 0;
  #define chamfers_tb (_Guide_SGS3->_parameters.chamfers_tb)
  _Guide_SGS3->_parameters.nelements = 1;
  #define nelements (_Guide_SGS3->_parameters.nelements)
  _Guide_SGS3->_parameters.nu = 0;
  #define nu (_Guide_SGS3->_parameters.nu)
  _Guide_SGS3->_parameters.phase = 0;
  #define phase (_Guide_SGS3->_parameters.phase)
  if("NULL" && strlen("NULL"))
    stracpy(_Guide_SGS3->_parameters.reflect, "NULL" ? "NULL" : "", 16384);
  else 
  _Guide_SGS3->_parameters.reflect[0]='\0';
  #define reflect (_Guide_SGS3->_parameters.reflect)

  #define GVars (_Guide_SGS3->_parameters.GVars)
  #define pTable (_Guide_SGS3->_parameters.pTable)

  #undef w1
  #undef h1
  #undef w2
  #undef h2
  #undef l
  #undef R0
  #undef Qc
  #undef alpha
  #undef m
  #undef W
  #undef nslit
  #undef d
  #undef mleft
  #undef mright
  #undef mtop
  #undef mbottom
  #undef nhslit
  #undef G
  #undef aleft
  #undef aright
  #undef atop
  #undef abottom
  #undef wavy
  #undef wavy_z
  #undef wavy_tb
  #undef wavy_lr
  #undef chamfers
  #undef chamfers_z
  #undef chamfers_lr
  #undef chamfers_tb
  #undef nelements
  #undef nu
  #undef phase
  #undef reflect
  #undef GVars
  #undef pTable
  /* component Guide_SGS3=Guide_gravity() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Guide_SGS2->_rotation_absolute, _Guide_SGS3->_rotation_absolute);
    rot_transpose(_Guide_SGS2->_rotation_absolute, tr1);
    rot_mul(_Guide_SGS3->_rotation_absolute, tr1, _Guide_SGS3->_rotation_relative);
    _Guide_SGS3->_rotation_is_identity =  rot_test_identity(_Guide_SGS3->_rotation_relative);
    tc1 = coords_set(
      0, 0, 5.065);
    rot_transpose(_Guide_SGS2->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Guide_SGS3->_position_absolute = coords_add(_Guide_SGS2->_position_absolute, tc2);
    tc1 = coords_sub(_Guide_SGS2->_position_absolute, _Guide_SGS3->_position_absolute);
    _Guide_SGS3->_position_relative = rot_apply(_Guide_SGS3->_rotation_absolute, tc1);
  } /* Guide_SGS3=Guide_gravity() AT ROTATED */
  DEBUG_COMPONENT("Guide_SGS3", _Guide_SGS3->_position_absolute, _Guide_SGS3->_rotation_absolute);
  instrument->_position_absolute[25] = _Guide_SGS3->_position_absolute;
  instrument->_position_relative[25] = _Guide_SGS3->_position_relative;
  instrument->counter_N[25]  = instrument->counter_P[25] = instrument->counter_P2[25] = 0;
  instrument->counter_AbsorbProp[25]= 0;
  return(0);
} /* _Guide_SGS3_setpos */

/* component Guide_CGS=Guide_gravity() SETTING, POSITION/ROTATION */
int _Guide_CGS_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Guide_CGS_setpos] component Guide_CGS=Guide_gravity() SETTING [/usr/share/mcstas/3.0-dev/optics/Guide_gravity.comp:339]");
  stracpy(_Guide_CGS->_name, "Guide_CGS", 16384);
  stracpy(_Guide_CGS->_type, "Guide_gravity", 16384);
  _Guide_CGS->_index=26;
  _Guide_CGS->_parameters.w1 = 0.03;
  #define w1 (_Guide_CGS->_parameters.w1)
  _Guide_CGS->_parameters.h1 = 0.055;
  #define h1 (_Guide_CGS->_parameters.h1)
  _Guide_CGS->_parameters.w2 = 0.015;
  #define w2 (_Guide_CGS->_parameters.w2)
  _Guide_CGS->_parameters.h2 = 0;
  #define h2 (_Guide_CGS->_parameters.h2)
  _Guide_CGS->_parameters.l = 2.143 -2 * 6e-3;
  #define l (_Guide_CGS->_parameters.l)
  _Guide_CGS->_parameters.R0 = 0.995;
  #define R0 (_Guide_CGS->_parameters.R0)
  _Guide_CGS->_parameters.Qc = 0.0218;
  #define Qc (_Guide_CGS->_parameters.Qc)
  _Guide_CGS->_parameters.alpha = 4.38;
  #define alpha (_Guide_CGS->_parameters.alpha)
  _Guide_CGS->_parameters.m = 2.4;
  #define m (_Guide_CGS->_parameters.m)
  _Guide_CGS->_parameters.W = 0.003;
  #define W (_Guide_CGS->_parameters.W)
  _Guide_CGS->_parameters.nslit = 1;
  #define nslit (_Guide_CGS->_parameters.nslit)
  _Guide_CGS->_parameters.d = 0.0005;
  #define d (_Guide_CGS->_parameters.d)
  _Guide_CGS->_parameters.mleft = -1;
  #define mleft (_Guide_CGS->_parameters.mleft)
  _Guide_CGS->_parameters.mright = -1;
  #define mright (_Guide_CGS->_parameters.mright)
  _Guide_CGS->_parameters.mtop = -1;
  #define mtop (_Guide_CGS->_parameters.mtop)
  _Guide_CGS->_parameters.mbottom = -1;
  #define mbottom (_Guide_CGS->_parameters.mbottom)
  _Guide_CGS->_parameters.nhslit = 1;
  #define nhslit (_Guide_CGS->_parameters.nhslit)
  _Guide_CGS->_parameters.G = 0;
  #define G (_Guide_CGS->_parameters.G)
  _Guide_CGS->_parameters.aleft = -1;
  #define aleft (_Guide_CGS->_parameters.aleft)
  _Guide_CGS->_parameters.aright = -1;
  #define aright (_Guide_CGS->_parameters.aright)
  _Guide_CGS->_parameters.atop = -1;
  #define atop (_Guide_CGS->_parameters.atop)
  _Guide_CGS->_parameters.abottom = -1;
  #define abottom (_Guide_CGS->_parameters.abottom)
  _Guide_CGS->_parameters.wavy = 0;
  #define wavy (_Guide_CGS->_parameters.wavy)
  _Guide_CGS->_parameters.wavy_z = 0;
  #define wavy_z (_Guide_CGS->_parameters.wavy_z)
  _Guide_CGS->_parameters.wavy_tb = 0;
  #define wavy_tb (_Guide_CGS->_parameters.wavy_tb)
  _Guide_CGS->_parameters.wavy_lr = 0;
  #define wavy_lr (_Guide_CGS->_parameters.wavy_lr)
  _Guide_CGS->_parameters.chamfers = 0;
  #define chamfers (_Guide_CGS->_parameters.chamfers)
  _Guide_CGS->_parameters.chamfers_z = 0;
  #define chamfers_z (_Guide_CGS->_parameters.chamfers_z)
  _Guide_CGS->_parameters.chamfers_lr = 0;
  #define chamfers_lr (_Guide_CGS->_parameters.chamfers_lr)
  _Guide_CGS->_parameters.chamfers_tb = 0;
  #define chamfers_tb (_Guide_CGS->_parameters.chamfers_tb)
  _Guide_CGS->_parameters.nelements = 1;
  #define nelements (_Guide_CGS->_parameters.nelements)
  _Guide_CGS->_parameters.nu = 0;
  #define nu (_Guide_CGS->_parameters.nu)
  _Guide_CGS->_parameters.phase = 0;
  #define phase (_Guide_CGS->_parameters.phase)
  if("NULL" && strlen("NULL"))
    stracpy(_Guide_CGS->_parameters.reflect, "NULL" ? "NULL" : "", 16384);
  else 
  _Guide_CGS->_parameters.reflect[0]='\0';
  #define reflect (_Guide_CGS->_parameters.reflect)

  #define GVars (_Guide_CGS->_parameters.GVars)
  #define pTable (_Guide_CGS->_parameters.pTable)

  #undef w1
  #undef h1
  #undef w2
  #undef h2
  #undef l
  #undef R0
  #undef Qc
  #undef alpha
  #undef m
  #undef W
  #undef nslit
  #undef d
  #undef mleft
  #undef mright
  #undef mtop
  #undef mbottom
  #undef nhslit
  #undef G
  #undef aleft
  #undef aright
  #undef atop
  #undef abottom
  #undef wavy
  #undef wavy_z
  #undef wavy_tb
  #undef wavy_lr
  #undef chamfers
  #undef chamfers_z
  #undef chamfers_lr
  #undef chamfers_tb
  #undef nelements
  #undef nu
  #undef phase
  #undef reflect
  #undef GVars
  #undef pTable
  /* component Guide_CGS=Guide_gravity() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Guide_SGS3->_rotation_absolute, _Guide_CGS->_rotation_absolute);
    rot_transpose(_Guide_SGS3->_rotation_absolute, tr1);
    rot_mul(_Guide_CGS->_rotation_absolute, tr1, _Guide_CGS->_rotation_relative);
    _Guide_CGS->_rotation_is_identity =  rot_test_identity(_Guide_CGS->_rotation_relative);
    tc1 = coords_set(
      0, 0, 2.352);
    rot_transpose(_Guide_SGS3->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Guide_CGS->_position_absolute = coords_add(_Guide_SGS3->_position_absolute, tc2);
    tc1 = coords_sub(_Guide_SGS3->_position_absolute, _Guide_CGS->_position_absolute);
    _Guide_CGS->_position_relative = rot_apply(_Guide_CGS->_rotation_absolute, tc1);
  } /* Guide_CGS=Guide_gravity() AT ROTATED */
  DEBUG_COMPONENT("Guide_CGS", _Guide_CGS->_position_absolute, _Guide_CGS->_rotation_absolute);
  instrument->_position_absolute[26] = _Guide_CGS->_position_absolute;
  instrument->_position_relative[26] = _Guide_CGS->_position_relative;
  instrument->counter_N[26]  = instrument->counter_P[26] = instrument->counter_P2[26] = 0;
  instrument->counter_AbsorbProp[26]= 0;
  return(0);
} /* _Guide_CGS_setpos */

/* component Guide2_PSD=Monitor_nD() SETTING, POSITION/ROTATION */
int _Guide2_PSD_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Guide2_PSD_setpos] component Guide2_PSD=Monitor_nD() SETTING [/usr/share/mcstas/3.0-dev/monitors/Monitor_nD.comp:228]");
  stracpy(_Guide2_PSD->_name, "Guide2_PSD", 16384);
  stracpy(_Guide2_PSD->_type, "Monitor_nD", 16384);
  _Guide2_PSD->_index=27;
  _Guide2_PSD->_parameters.user1 = '\0';
  #define user1 (_Guide2_PSD->_parameters.user1)
  _Guide2_PSD->_parameters.user2 = '\0';
  #define user2 (_Guide2_PSD->_parameters.user2)
  _Guide2_PSD->_parameters.user3 = '\0';
  #define user3 (_Guide2_PSD->_parameters.user3)
  _Guide2_PSD->_parameters.xwidth = 0.03;
  #define xwidth (_Guide2_PSD->_parameters.xwidth)
  _Guide2_PSD->_parameters.yheight = 0.055;
  #define yheight (_Guide2_PSD->_parameters.yheight)
  _Guide2_PSD->_parameters.zdepth = 0;
  #define zdepth (_Guide2_PSD->_parameters.zdepth)
  _Guide2_PSD->_parameters.xmin = 0;
  #define xmin (_Guide2_PSD->_parameters.xmin)
  _Guide2_PSD->_parameters.xmax = 0;
  #define xmax (_Guide2_PSD->_parameters.xmax)
  _Guide2_PSD->_parameters.ymin = 0;
  #define ymin (_Guide2_PSD->_parameters.ymin)
  _Guide2_PSD->_parameters.ymax = 0;
  #define ymax (_Guide2_PSD->_parameters.ymax)
  _Guide2_PSD->_parameters.zmin = 0;
  #define zmin (_Guide2_PSD->_parameters.zmin)
  _Guide2_PSD->_parameters.zmax = 0;
  #define zmax (_Guide2_PSD->_parameters.zmax)
  _Guide2_PSD->_parameters.bins = 0;
  #define bins (_Guide2_PSD->_parameters.bins)
  _Guide2_PSD->_parameters.min = -1e40;
  #define min (_Guide2_PSD->_parameters.min)
  _Guide2_PSD->_parameters.max = 1e40;
  #define max (_Guide2_PSD->_parameters.max)
  _Guide2_PSD->_parameters.restore_neutron = 0;
  #define restore_neutron (_Guide2_PSD->_parameters.restore_neutron)
  _Guide2_PSD->_parameters.radius = 0;
  #define radius (_Guide2_PSD->_parameters.radius)
  if("x y" && strlen("x y"))
    stracpy(_Guide2_PSD->_parameters.options, "x y" ? "x y" : "", 16384);
  else 
  _Guide2_PSD->_parameters.options[0]='\0';
  #define options (_Guide2_PSD->_parameters.options)
  if("NULL" && strlen("NULL"))
    stracpy(_Guide2_PSD->_parameters.filename, "NULL" ? "NULL" : "", 16384);
  else 
  _Guide2_PSD->_parameters.filename[0]='\0';
  #define filename (_Guide2_PSD->_parameters.filename)
  if("NULL" && strlen("NULL"))
    stracpy(_Guide2_PSD->_parameters.geometry, "NULL" ? "NULL" : "", 16384);
  else 
  _Guide2_PSD->_parameters.geometry[0]='\0';
  #define geometry (_Guide2_PSD->_parameters.geometry)
  if("NULL" && strlen("NULL"))
    stracpy(_Guide2_PSD->_parameters.username1, "NULL" ? "NULL" : "", 16384);
  else 
  _Guide2_PSD->_parameters.username1[0]='\0';
  #define username1 (_Guide2_PSD->_parameters.username1)
  if("NULL" && strlen("NULL"))
    stracpy(_Guide2_PSD->_parameters.username2, "NULL" ? "NULL" : "", 16384);
  else 
  _Guide2_PSD->_parameters.username2[0]='\0';
  #define username2 (_Guide2_PSD->_parameters.username2)
  if("NULL" && strlen("NULL"))
    stracpy(_Guide2_PSD->_parameters.username3, "NULL" ? "NULL" : "", 16384);
  else 
  _Guide2_PSD->_parameters.username3[0]='\0';
  #define username3 (_Guide2_PSD->_parameters.username3)

  #define DEFS (_Guide2_PSD->_parameters.DEFS)
  #define Vars (_Guide2_PSD->_parameters.Vars)
  #define detector (_Guide2_PSD->_parameters.detector)
  #define offdata (_Guide2_PSD->_parameters.offdata)

  #undef user1
  #undef user2
  #undef user3
  #undef xwidth
  #undef yheight
  #undef zdepth
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef zmin
  #undef zmax
  #undef bins
  #undef min
  #undef max
  #undef restore_neutron
  #undef radius
  #undef options
  #undef filename
  #undef geometry
  #undef username1
  #undef username2
  #undef username3
  #undef DEFS
  #undef Vars
  #undef detector
  #undef offdata
  /* component Guide2_PSD=Monitor_nD() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Guide_CGS->_rotation_absolute, _Guide2_PSD->_rotation_absolute);
    rot_transpose(_Guide_CGS->_rotation_absolute, tr1);
    rot_mul(_Guide2_PSD->_rotation_absolute, tr1, _Guide2_PSD->_rotation_relative);
    _Guide2_PSD->_rotation_is_identity =  rot_test_identity(_Guide2_PSD->_rotation_relative);
    tc1 = coords_set(
      0, 0, 2.143 -2 * 6e-3 + 1e-4);
    rot_transpose(_Guide_CGS->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Guide2_PSD->_position_absolute = coords_add(_Guide_CGS->_position_absolute, tc2);
    tc1 = coords_sub(_Guide_CGS->_position_absolute, _Guide2_PSD->_position_absolute);
    _Guide2_PSD->_position_relative = rot_apply(_Guide2_PSD->_rotation_absolute, tc1);
  } /* Guide2_PSD=Monitor_nD() AT ROTATED */
  DEBUG_COMPONENT("Guide2_PSD", _Guide2_PSD->_position_absolute, _Guide2_PSD->_rotation_absolute);
  instrument->_position_absolute[27] = _Guide2_PSD->_position_absolute;
  instrument->_position_relative[27] = _Guide2_PSD->_position_relative;
  instrument->counter_N[27]  = instrument->counter_P[27] = instrument->counter_P2[27] = 0;
  instrument->counter_AbsorbProp[27]= 0;
  return(0);
} /* _Guide2_PSD_setpos */

/* component Guide2_Lambda=Monitor_nD() SETTING, POSITION/ROTATION */
int _Guide2_Lambda_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Guide2_Lambda_setpos] component Guide2_Lambda=Monitor_nD() SETTING [/usr/share/mcstas/3.0-dev/monitors/Monitor_nD.comp:228]");
  stracpy(_Guide2_Lambda->_name, "Guide2_Lambda", 16384);
  stracpy(_Guide2_Lambda->_type, "Monitor_nD", 16384);
  _Guide2_Lambda->_index=28;
  _Guide2_Lambda->_parameters.user1 = '\0';
  #define user1 (_Guide2_Lambda->_parameters.user1)
  _Guide2_Lambda->_parameters.user2 = '\0';
  #define user2 (_Guide2_Lambda->_parameters.user2)
  _Guide2_Lambda->_parameters.user3 = '\0';
  #define user3 (_Guide2_Lambda->_parameters.user3)
  _Guide2_Lambda->_parameters.xwidth = 0.03;
  #define xwidth (_Guide2_Lambda->_parameters.xwidth)
  _Guide2_Lambda->_parameters.yheight = 0.055;
  #define yheight (_Guide2_Lambda->_parameters.yheight)
  _Guide2_Lambda->_parameters.zdepth = 0;
  #define zdepth (_Guide2_Lambda->_parameters.zdepth)
  _Guide2_Lambda->_parameters.xmin = 0;
  #define xmin (_Guide2_Lambda->_parameters.xmin)
  _Guide2_Lambda->_parameters.xmax = 0;
  #define xmax (_Guide2_Lambda->_parameters.xmax)
  _Guide2_Lambda->_parameters.ymin = 0;
  #define ymin (_Guide2_Lambda->_parameters.ymin)
  _Guide2_Lambda->_parameters.ymax = 0;
  #define ymax (_Guide2_Lambda->_parameters.ymax)
  _Guide2_Lambda->_parameters.zmin = 0;
  #define zmin (_Guide2_Lambda->_parameters.zmin)
  _Guide2_Lambda->_parameters.zmax = 0;
  #define zmax (_Guide2_Lambda->_parameters.zmax)
  _Guide2_Lambda->_parameters.bins = 0;
  #define bins (_Guide2_Lambda->_parameters.bins)
  _Guide2_Lambda->_parameters.min = -1e40;
  #define min (_Guide2_Lambda->_parameters.min)
  _Guide2_Lambda->_parameters.max = 1e40;
  #define max (_Guide2_Lambda->_parameters.max)
  _Guide2_Lambda->_parameters.restore_neutron = 0;
  #define restore_neutron (_Guide2_Lambda->_parameters.restore_neutron)
  _Guide2_Lambda->_parameters.radius = 0;
  #define radius (_Guide2_Lambda->_parameters.radius)
  if("lambda, auto" && strlen("lambda, auto"))
    stracpy(_Guide2_Lambda->_parameters.options, "lambda, auto" ? "lambda, auto" : "", 16384);
  else 
  _Guide2_Lambda->_parameters.options[0]='\0';
  #define options (_Guide2_Lambda->_parameters.options)
  if("NULL" && strlen("NULL"))
    stracpy(_Guide2_Lambda->_parameters.filename, "NULL" ? "NULL" : "", 16384);
  else 
  _Guide2_Lambda->_parameters.filename[0]='\0';
  #define filename (_Guide2_Lambda->_parameters.filename)
  if("NULL" && strlen("NULL"))
    stracpy(_Guide2_Lambda->_parameters.geometry, "NULL" ? "NULL" : "", 16384);
  else 
  _Guide2_Lambda->_parameters.geometry[0]='\0';
  #define geometry (_Guide2_Lambda->_parameters.geometry)
  if("NULL" && strlen("NULL"))
    stracpy(_Guide2_Lambda->_parameters.username1, "NULL" ? "NULL" : "", 16384);
  else 
  _Guide2_Lambda->_parameters.username1[0]='\0';
  #define username1 (_Guide2_Lambda->_parameters.username1)
  if("NULL" && strlen("NULL"))
    stracpy(_Guide2_Lambda->_parameters.username2, "NULL" ? "NULL" : "", 16384);
  else 
  _Guide2_Lambda->_parameters.username2[0]='\0';
  #define username2 (_Guide2_Lambda->_parameters.username2)
  if("NULL" && strlen("NULL"))
    stracpy(_Guide2_Lambda->_parameters.username3, "NULL" ? "NULL" : "", 16384);
  else 
  _Guide2_Lambda->_parameters.username3[0]='\0';
  #define username3 (_Guide2_Lambda->_parameters.username3)

  #define DEFS (_Guide2_Lambda->_parameters.DEFS)
  #define Vars (_Guide2_Lambda->_parameters.Vars)
  #define detector (_Guide2_Lambda->_parameters.detector)
  #define offdata (_Guide2_Lambda->_parameters.offdata)

  #undef user1
  #undef user2
  #undef user3
  #undef xwidth
  #undef yheight
  #undef zdepth
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef zmin
  #undef zmax
  #undef bins
  #undef min
  #undef max
  #undef restore_neutron
  #undef radius
  #undef options
  #undef filename
  #undef geometry
  #undef username1
  #undef username2
  #undef username3
  #undef DEFS
  #undef Vars
  #undef detector
  #undef offdata
  /* component Guide2_Lambda=Monitor_nD() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Guide2_PSD->_rotation_absolute, _Guide2_Lambda->_rotation_absolute);
    rot_transpose(_Guide2_PSD->_rotation_absolute, tr1);
    rot_mul(_Guide2_Lambda->_rotation_absolute, tr1, _Guide2_Lambda->_rotation_relative);
    _Guide2_Lambda->_rotation_is_identity =  rot_test_identity(_Guide2_Lambda->_rotation_relative);
    tc1 = coords_set(
      0, 0, 1e-4);
    rot_transpose(_Guide2_PSD->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Guide2_Lambda->_position_absolute = coords_add(_Guide2_PSD->_position_absolute, tc2);
    tc1 = coords_sub(_Guide2_PSD->_position_absolute, _Guide2_Lambda->_position_absolute);
    _Guide2_Lambda->_position_relative = rot_apply(_Guide2_Lambda->_rotation_absolute, tc1);
  } /* Guide2_Lambda=Monitor_nD() AT ROTATED */
  DEBUG_COMPONENT("Guide2_Lambda", _Guide2_Lambda->_position_absolute, _Guide2_Lambda->_rotation_absolute);
  instrument->_position_absolute[28] = _Guide2_Lambda->_position_absolute;
  instrument->_position_relative[28] = _Guide2_Lambda->_position_relative;
  instrument->counter_N[28]  = instrument->counter_P[28] = instrument->counter_P2[28] = 0;
  instrument->counter_AbsorbProp[28]= 0;
  return(0);
} /* _Guide2_Lambda_setpos */

/* component Guide2_dXY=Monitor_nD() SETTING, POSITION/ROTATION */
int _Guide2_dXY_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Guide2_dXY_setpos] component Guide2_dXY=Monitor_nD() SETTING [/usr/share/mcstas/3.0-dev/monitors/Monitor_nD.comp:228]");
  stracpy(_Guide2_dXY->_name, "Guide2_dXY", 16384);
  stracpy(_Guide2_dXY->_type, "Monitor_nD", 16384);
  _Guide2_dXY->_index=29;
  _Guide2_dXY->_parameters.user1 = '\0';
  #define user1 (_Guide2_dXY->_parameters.user1)
  _Guide2_dXY->_parameters.user2 = '\0';
  #define user2 (_Guide2_dXY->_parameters.user2)
  _Guide2_dXY->_parameters.user3 = '\0';
  #define user3 (_Guide2_dXY->_parameters.user3)
  _Guide2_dXY->_parameters.xwidth = 0.03;
  #define xwidth (_Guide2_dXY->_parameters.xwidth)
  _Guide2_dXY->_parameters.yheight = 0.055;
  #define yheight (_Guide2_dXY->_parameters.yheight)
  _Guide2_dXY->_parameters.zdepth = 0;
  #define zdepth (_Guide2_dXY->_parameters.zdepth)
  _Guide2_dXY->_parameters.xmin = 0;
  #define xmin (_Guide2_dXY->_parameters.xmin)
  _Guide2_dXY->_parameters.xmax = 0;
  #define xmax (_Guide2_dXY->_parameters.xmax)
  _Guide2_dXY->_parameters.ymin = 0;
  #define ymin (_Guide2_dXY->_parameters.ymin)
  _Guide2_dXY->_parameters.ymax = 0;
  #define ymax (_Guide2_dXY->_parameters.ymax)
  _Guide2_dXY->_parameters.zmin = 0;
  #define zmin (_Guide2_dXY->_parameters.zmin)
  _Guide2_dXY->_parameters.zmax = 0;
  #define zmax (_Guide2_dXY->_parameters.zmax)
  _Guide2_dXY->_parameters.bins = 0;
  #define bins (_Guide2_dXY->_parameters.bins)
  _Guide2_dXY->_parameters.min = -1e40;
  #define min (_Guide2_dXY->_parameters.min)
  _Guide2_dXY->_parameters.max = 1e40;
  #define max (_Guide2_dXY->_parameters.max)
  _Guide2_dXY->_parameters.restore_neutron = 0;
  #define restore_neutron (_Guide2_dXY->_parameters.restore_neutron)
  _Guide2_dXY->_parameters.radius = 0;
  #define radius (_Guide2_dXY->_parameters.radius)
  if("dx dy" && strlen("dx dy"))
    stracpy(_Guide2_dXY->_parameters.options, "dx dy" ? "dx dy" : "", 16384);
  else 
  _Guide2_dXY->_parameters.options[0]='\0';
  #define options (_Guide2_dXY->_parameters.options)
  if("NULL" && strlen("NULL"))
    stracpy(_Guide2_dXY->_parameters.filename, "NULL" ? "NULL" : "", 16384);
  else 
  _Guide2_dXY->_parameters.filename[0]='\0';
  #define filename (_Guide2_dXY->_parameters.filename)
  if("NULL" && strlen("NULL"))
    stracpy(_Guide2_dXY->_parameters.geometry, "NULL" ? "NULL" : "", 16384);
  else 
  _Guide2_dXY->_parameters.geometry[0]='\0';
  #define geometry (_Guide2_dXY->_parameters.geometry)
  if("NULL" && strlen("NULL"))
    stracpy(_Guide2_dXY->_parameters.username1, "NULL" ? "NULL" : "", 16384);
  else 
  _Guide2_dXY->_parameters.username1[0]='\0';
  #define username1 (_Guide2_dXY->_parameters.username1)
  if("NULL" && strlen("NULL"))
    stracpy(_Guide2_dXY->_parameters.username2, "NULL" ? "NULL" : "", 16384);
  else 
  _Guide2_dXY->_parameters.username2[0]='\0';
  #define username2 (_Guide2_dXY->_parameters.username2)
  if("NULL" && strlen("NULL"))
    stracpy(_Guide2_dXY->_parameters.username3, "NULL" ? "NULL" : "", 16384);
  else 
  _Guide2_dXY->_parameters.username3[0]='\0';
  #define username3 (_Guide2_dXY->_parameters.username3)

  #define DEFS (_Guide2_dXY->_parameters.DEFS)
  #define Vars (_Guide2_dXY->_parameters.Vars)
  #define detector (_Guide2_dXY->_parameters.detector)
  #define offdata (_Guide2_dXY->_parameters.offdata)

  #undef user1
  #undef user2
  #undef user3
  #undef xwidth
  #undef yheight
  #undef zdepth
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef zmin
  #undef zmax
  #undef bins
  #undef min
  #undef max
  #undef restore_neutron
  #undef radius
  #undef options
  #undef filename
  #undef geometry
  #undef username1
  #undef username2
  #undef username3
  #undef DEFS
  #undef Vars
  #undef detector
  #undef offdata
  /* component Guide2_dXY=Monitor_nD() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Guide2_Lambda->_rotation_absolute, _Guide2_dXY->_rotation_absolute);
    rot_transpose(_Guide2_Lambda->_rotation_absolute, tr1);
    rot_mul(_Guide2_dXY->_rotation_absolute, tr1, _Guide2_dXY->_rotation_relative);
    _Guide2_dXY->_rotation_is_identity =  rot_test_identity(_Guide2_dXY->_rotation_relative);
    tc1 = coords_set(
      0, 0, 1e-4);
    rot_transpose(_Guide2_Lambda->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Guide2_dXY->_position_absolute = coords_add(_Guide2_Lambda->_position_absolute, tc2);
    tc1 = coords_sub(_Guide2_Lambda->_position_absolute, _Guide2_dXY->_position_absolute);
    _Guide2_dXY->_position_relative = rot_apply(_Guide2_dXY->_rotation_absolute, tc1);
  } /* Guide2_dXY=Monitor_nD() AT ROTATED */
  DEBUG_COMPONENT("Guide2_dXY", _Guide2_dXY->_position_absolute, _Guide2_dXY->_rotation_absolute);
  instrument->_position_absolute[29] = _Guide2_dXY->_position_absolute;
  instrument->_position_relative[29] = _Guide2_dXY->_position_relative;
  instrument->counter_N[29]  = instrument->counter_P[29] = instrument->counter_P2[29] = 0;
  instrument->counter_AbsorbProp[29]= 0;
  return(0);
} /* _Guide2_dXY_setpos */

/* component Chopper6=DiskChopper() SETTING, POSITION/ROTATION */
int _Chopper6_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Chopper6_setpos] component Chopper6=DiskChopper() SETTING [/usr/share/mcstas/3.0-dev/optics/DiskChopper.comp:67]");
  stracpy(_Chopper6->_name, "Chopper6", 16384);
  stracpy(_Chopper6->_type, "DiskChopper", 16384);
  _Chopper6->_index=30;
  _Chopper6->_parameters.theta_0 = 0;
  #define theta_0 (_Chopper6->_parameters.theta_0)
  _Chopper6->_parameters.radius = 0.24;
  #define radius (_Chopper6->_parameters.radius)
  _Chopper6->_parameters.yheight = 0.06;
  #define yheight (_Chopper6->_parameters.yheight)
  _Chopper6->_parameters.nu = instrument->_parameters._rpm / 60;
  #define nu (_Chopper6->_parameters.nu)
  _Chopper6->_parameters.nslit = 3;
  #define nslit (_Chopper6->_parameters.nslit)
  _Chopper6->_parameters.jitter = 0;
  #define jitter (_Chopper6->_parameters.jitter)
  _Chopper6->_parameters.delay = 252.78 * 11.97 * instrument->_parameters._lambda * 1e-6;
  #define delay (_Chopper6->_parameters.delay)
  _Chopper6->_parameters.isfirst = 0;
  #define isfirst (_Chopper6->_parameters.isfirst)
  _Chopper6->_parameters.n_pulse = 1;
  #define n_pulse (_Chopper6->_parameters.n_pulse)
  _Chopper6->_parameters.abs_out = 1;
  #define abs_out (_Chopper6->_parameters.abs_out)
  _Chopper6->_parameters.phase = 0;
  #define phase (_Chopper6->_parameters.phase)
  _Chopper6->_parameters.xwidth = 0.03;
  #define xwidth (_Chopper6->_parameters.xwidth)
  _Chopper6->_parameters.verbose = 0;
  #define verbose (_Chopper6->_parameters.verbose)

  #define Tg (_Chopper6->_parameters.Tg)
  #define To (_Chopper6->_parameters.To)
  #define delta_y (_Chopper6->_parameters.delta_y)
  #define height (_Chopper6->_parameters.height)
  #define omega (_Chopper6->_parameters.omega)

  #undef theta_0
  #undef radius
  #undef yheight
  #undef nu
  #undef nslit
  #undef jitter
  #undef delay
  #undef isfirst
  #undef n_pulse
  #undef abs_out
  #undef phase
  #undef xwidth
  #undef verbose
  #undef Tg
  #undef To
  #undef delta_y
  #undef height
  #undef omega
  /* component Chopper6=DiskChopper() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Guide2_dXY->_rotation_absolute, _Chopper6->_rotation_absolute);
    rot_transpose(_Guide2_dXY->_rotation_absolute, tr1);
    rot_mul(_Chopper6->_rotation_absolute, tr1, _Chopper6->_rotation_relative);
    _Chopper6->_rotation_is_identity =  rot_test_identity(_Chopper6->_rotation_relative);
    tc1 = coords_set(
      0, 0, 1.2e-2 / 2);
    rot_transpose(_Guide2_dXY->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Chopper6->_position_absolute = coords_add(_Guide2_dXY->_position_absolute, tc2);
    tc1 = coords_sub(_Guide2_dXY->_position_absolute, _Chopper6->_position_absolute);
    _Chopper6->_position_relative = rot_apply(_Chopper6->_rotation_absolute, tc1);
  } /* Chopper6=DiskChopper() AT ROTATED */
  DEBUG_COMPONENT("Chopper6", _Chopper6->_position_absolute, _Chopper6->_rotation_absolute);
  instrument->_position_absolute[30] = _Chopper6->_position_absolute;
  instrument->_position_relative[30] = _Chopper6->_position_relative;
  instrument->counter_N[30]  = instrument->counter_P[30] = instrument->counter_P2[30] = 0;
  instrument->counter_AbsorbProp[30]= 0;
  return(0);
} /* _Chopper6_setpos */

/* component Guide_time=Monitor_nD() SETTING, POSITION/ROTATION */
int _Guide_time_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Guide_time_setpos] component Guide_time=Monitor_nD() SETTING [/usr/share/mcstas/3.0-dev/monitors/Monitor_nD.comp:228]");
  stracpy(_Guide_time->_name, "Guide_time", 16384);
  stracpy(_Guide_time->_type, "Monitor_nD", 16384);
  _Guide_time->_index=31;
  _Guide_time->_parameters.user1 = '\0';
  #define user1 (_Guide_time->_parameters.user1)
  _Guide_time->_parameters.user2 = '\0';
  #define user2 (_Guide_time->_parameters.user2)
  _Guide_time->_parameters.user3 = '\0';
  #define user3 (_Guide_time->_parameters.user3)
  _Guide_time->_parameters.xwidth = 0.03;
  #define xwidth (_Guide_time->_parameters.xwidth)
  _Guide_time->_parameters.yheight = 0.055;
  #define yheight (_Guide_time->_parameters.yheight)
  _Guide_time->_parameters.zdepth = 0;
  #define zdepth (_Guide_time->_parameters.zdepth)
  _Guide_time->_parameters.xmin = 0;
  #define xmin (_Guide_time->_parameters.xmin)
  _Guide_time->_parameters.xmax = 0;
  #define xmax (_Guide_time->_parameters.xmax)
  _Guide_time->_parameters.ymin = 0;
  #define ymin (_Guide_time->_parameters.ymin)
  _Guide_time->_parameters.ymax = 0;
  #define ymax (_Guide_time->_parameters.ymax)
  _Guide_time->_parameters.zmin = 0;
  #define zmin (_Guide_time->_parameters.zmin)
  _Guide_time->_parameters.zmax = 0;
  #define zmax (_Guide_time->_parameters.zmax)
  _Guide_time->_parameters.bins = 0;
  #define bins (_Guide_time->_parameters.bins)
  _Guide_time->_parameters.min = -1e40;
  #define min (_Guide_time->_parameters.min)
  _Guide_time->_parameters.max = 1e40;
  #define max (_Guide_time->_parameters.max)
  _Guide_time->_parameters.restore_neutron = 0;
  #define restore_neutron (_Guide_time->_parameters.restore_neutron)
  _Guide_time->_parameters.radius = 0;
  #define radius (_Guide_time->_parameters.radius)
  if("time auto" && strlen("time auto"))
    stracpy(_Guide_time->_parameters.options, "time auto" ? "time auto" : "", 16384);
  else 
  _Guide_time->_parameters.options[0]='\0';
  #define options (_Guide_time->_parameters.options)
  if("NULL" && strlen("NULL"))
    stracpy(_Guide_time->_parameters.filename, "NULL" ? "NULL" : "", 16384);
  else 
  _Guide_time->_parameters.filename[0]='\0';
  #define filename (_Guide_time->_parameters.filename)
  if("NULL" && strlen("NULL"))
    stracpy(_Guide_time->_parameters.geometry, "NULL" ? "NULL" : "", 16384);
  else 
  _Guide_time->_parameters.geometry[0]='\0';
  #define geometry (_Guide_time->_parameters.geometry)
  if("NULL" && strlen("NULL"))
    stracpy(_Guide_time->_parameters.username1, "NULL" ? "NULL" : "", 16384);
  else 
  _Guide_time->_parameters.username1[0]='\0';
  #define username1 (_Guide_time->_parameters.username1)
  if("NULL" && strlen("NULL"))
    stracpy(_Guide_time->_parameters.username2, "NULL" ? "NULL" : "", 16384);
  else 
  _Guide_time->_parameters.username2[0]='\0';
  #define username2 (_Guide_time->_parameters.username2)
  if("NULL" && strlen("NULL"))
    stracpy(_Guide_time->_parameters.username3, "NULL" ? "NULL" : "", 16384);
  else 
  _Guide_time->_parameters.username3[0]='\0';
  #define username3 (_Guide_time->_parameters.username3)

  #define DEFS (_Guide_time->_parameters.DEFS)
  #define Vars (_Guide_time->_parameters.Vars)
  #define detector (_Guide_time->_parameters.detector)
  #define offdata (_Guide_time->_parameters.offdata)

  #undef user1
  #undef user2
  #undef user3
  #undef xwidth
  #undef yheight
  #undef zdepth
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef zmin
  #undef zmax
  #undef bins
  #undef min
  #undef max
  #undef restore_neutron
  #undef radius
  #undef options
  #undef filename
  #undef geometry
  #undef username1
  #undef username2
  #undef username3
  #undef DEFS
  #undef Vars
  #undef detector
  #undef offdata
  /* component Guide_time=Monitor_nD() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Chopper6->_rotation_absolute, _Guide_time->_rotation_absolute);
    rot_transpose(_Chopper6->_rotation_absolute, tr1);
    rot_mul(_Guide_time->_rotation_absolute, tr1, _Guide_time->_rotation_relative);
    _Guide_time->_rotation_is_identity =  rot_test_identity(_Guide_time->_rotation_relative);
    tc1 = coords_set(
      0, 0, 1.2e-2 / 2);
    rot_transpose(_Chopper6->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Guide_time->_position_absolute = coords_add(_Chopper6->_position_absolute, tc2);
    tc1 = coords_sub(_Chopper6->_position_absolute, _Guide_time->_position_absolute);
    _Guide_time->_position_relative = rot_apply(_Guide_time->_rotation_absolute, tc1);
  } /* Guide_time=Monitor_nD() AT ROTATED */
  DEBUG_COMPONENT("Guide_time", _Guide_time->_position_absolute, _Guide_time->_rotation_absolute);
  instrument->_position_absolute[31] = _Guide_time->_position_absolute;
  instrument->_position_relative[31] = _Guide_time->_position_relative;
  instrument->counter_N[31]  = instrument->counter_P[31] = instrument->counter_P2[31] = 0;
  instrument->counter_AbsorbProp[31]= 0;
  return(0);
} /* _Guide_time_setpos */

/* component Guide_DGS=Guide_gravity() SETTING, POSITION/ROTATION */
int _Guide_DGS_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Guide_DGS_setpos] component Guide_DGS=Guide_gravity() SETTING [/usr/share/mcstas/3.0-dev/optics/Guide_gravity.comp:339]");
  stracpy(_Guide_DGS->_name, "Guide_DGS", 16384);
  stracpy(_Guide_DGS->_type, "Guide_gravity", 16384);
  _Guide_DGS->_index=32;
  _Guide_DGS->_parameters.w1 = 0.015;
  #define w1 (_Guide_DGS->_parameters.w1)
  _Guide_DGS->_parameters.h1 = 0.055;
  #define h1 (_Guide_DGS->_parameters.h1)
  _Guide_DGS->_parameters.w2 = 0.023;
  #define w2 (_Guide_DGS->_parameters.w2)
  _Guide_DGS->_parameters.h2 = 0;
  #define h2 (_Guide_DGS->_parameters.h2)
  _Guide_DGS->_parameters.l = 1.0 -2 * 6e-3;
  #define l (_Guide_DGS->_parameters.l)
  _Guide_DGS->_parameters.R0 = 0.995;
  #define R0 (_Guide_DGS->_parameters.R0)
  _Guide_DGS->_parameters.Qc = 0.0218;
  #define Qc (_Guide_DGS->_parameters.Qc)
  _Guide_DGS->_parameters.alpha = 4.38;
  #define alpha (_Guide_DGS->_parameters.alpha)
  _Guide_DGS->_parameters.m = 2.4;
  #define m (_Guide_DGS->_parameters.m)
  _Guide_DGS->_parameters.W = 0.003;
  #define W (_Guide_DGS->_parameters.W)
  _Guide_DGS->_parameters.nslit = 1;
  #define nslit (_Guide_DGS->_parameters.nslit)
  _Guide_DGS->_parameters.d = 0.0005;
  #define d (_Guide_DGS->_parameters.d)
  _Guide_DGS->_parameters.mleft = -1;
  #define mleft (_Guide_DGS->_parameters.mleft)
  _Guide_DGS->_parameters.mright = -1;
  #define mright (_Guide_DGS->_parameters.mright)
  _Guide_DGS->_parameters.mtop = -1;
  #define mtop (_Guide_DGS->_parameters.mtop)
  _Guide_DGS->_parameters.mbottom = -1;
  #define mbottom (_Guide_DGS->_parameters.mbottom)
  _Guide_DGS->_parameters.nhslit = 1;
  #define nhslit (_Guide_DGS->_parameters.nhslit)
  _Guide_DGS->_parameters.G = 0;
  #define G (_Guide_DGS->_parameters.G)
  _Guide_DGS->_parameters.aleft = -1;
  #define aleft (_Guide_DGS->_parameters.aleft)
  _Guide_DGS->_parameters.aright = -1;
  #define aright (_Guide_DGS->_parameters.aright)
  _Guide_DGS->_parameters.atop = -1;
  #define atop (_Guide_DGS->_parameters.atop)
  _Guide_DGS->_parameters.abottom = -1;
  #define abottom (_Guide_DGS->_parameters.abottom)
  _Guide_DGS->_parameters.wavy = 0;
  #define wavy (_Guide_DGS->_parameters.wavy)
  _Guide_DGS->_parameters.wavy_z = 0;
  #define wavy_z (_Guide_DGS->_parameters.wavy_z)
  _Guide_DGS->_parameters.wavy_tb = 0;
  #define wavy_tb (_Guide_DGS->_parameters.wavy_tb)
  _Guide_DGS->_parameters.wavy_lr = 0;
  #define wavy_lr (_Guide_DGS->_parameters.wavy_lr)
  _Guide_DGS->_parameters.chamfers = 0;
  #define chamfers (_Guide_DGS->_parameters.chamfers)
  _Guide_DGS->_parameters.chamfers_z = 0;
  #define chamfers_z (_Guide_DGS->_parameters.chamfers_z)
  _Guide_DGS->_parameters.chamfers_lr = 0;
  #define chamfers_lr (_Guide_DGS->_parameters.chamfers_lr)
  _Guide_DGS->_parameters.chamfers_tb = 0;
  #define chamfers_tb (_Guide_DGS->_parameters.chamfers_tb)
  _Guide_DGS->_parameters.nelements = 1;
  #define nelements (_Guide_DGS->_parameters.nelements)
  _Guide_DGS->_parameters.nu = 0;
  #define nu (_Guide_DGS->_parameters.nu)
  _Guide_DGS->_parameters.phase = 0;
  #define phase (_Guide_DGS->_parameters.phase)
  if("NULL" && strlen("NULL"))
    stracpy(_Guide_DGS->_parameters.reflect, "NULL" ? "NULL" : "", 16384);
  else 
  _Guide_DGS->_parameters.reflect[0]='\0';
  #define reflect (_Guide_DGS->_parameters.reflect)

  #define GVars (_Guide_DGS->_parameters.GVars)
  #define pTable (_Guide_DGS->_parameters.pTable)

  #undef w1
  #undef h1
  #undef w2
  #undef h2
  #undef l
  #undef R0
  #undef Qc
  #undef alpha
  #undef m
  #undef W
  #undef nslit
  #undef d
  #undef mleft
  #undef mright
  #undef mtop
  #undef mbottom
  #undef nhslit
  #undef G
  #undef aleft
  #undef aright
  #undef atop
  #undef abottom
  #undef wavy
  #undef wavy_z
  #undef wavy_tb
  #undef wavy_lr
  #undef chamfers
  #undef chamfers_z
  #undef chamfers_lr
  #undef chamfers_tb
  #undef nelements
  #undef nu
  #undef phase
  #undef reflect
  #undef GVars
  #undef pTable
  /* component Guide_DGS=Guide_gravity() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Guide_CGS->_rotation_absolute, _Guide_DGS->_rotation_absolute);
    rot_transpose(_Guide_time->_rotation_absolute, tr1);
    rot_mul(_Guide_DGS->_rotation_absolute, tr1, _Guide_DGS->_rotation_relative);
    _Guide_DGS->_rotation_is_identity =  rot_test_identity(_Guide_DGS->_rotation_relative);
    tc1 = coords_set(
      0, 0, 2.143);
    rot_transpose(_Guide_CGS->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Guide_DGS->_position_absolute = coords_add(_Guide_CGS->_position_absolute, tc2);
    tc1 = coords_sub(_Guide_time->_position_absolute, _Guide_DGS->_position_absolute);
    _Guide_DGS->_position_relative = rot_apply(_Guide_DGS->_rotation_absolute, tc1);
  } /* Guide_DGS=Guide_gravity() AT ROTATED */
  DEBUG_COMPONENT("Guide_DGS", _Guide_DGS->_position_absolute, _Guide_DGS->_rotation_absolute);
  instrument->_position_absolute[32] = _Guide_DGS->_position_absolute;
  instrument->_position_relative[32] = _Guide_DGS->_position_relative;
  instrument->counter_N[32]  = instrument->counter_P[32] = instrument->counter_P2[32] = 0;
  instrument->counter_AbsorbProp[32]= 0;
  return(0);
} /* _Guide_DGS_setpos */

/* component Sample_pos=Arm() SETTING, POSITION/ROTATION */
int _Sample_pos_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Sample_pos_setpos] component Sample_pos=Arm() SETTING [Arm:0]");
  stracpy(_Sample_pos->_name, "Sample_pos", 16384);
  stracpy(_Sample_pos->_type, "Arm", 16384);
  _Sample_pos->_index=33;
  /* component Sample_pos=Arm() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Guide_DGS->_rotation_absolute, _Sample_pos->_rotation_absolute);
    rot_transpose(_Guide_DGS->_rotation_absolute, tr1);
    rot_mul(_Sample_pos->_rotation_absolute, tr1, _Sample_pos->_rotation_relative);
    _Sample_pos->_rotation_is_identity =  rot_test_identity(_Sample_pos->_rotation_relative);
    tc1 = coords_set(
      0, 0, 2.143 + 0.165);
    rot_transpose(_Guide_DGS->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Sample_pos->_position_absolute = coords_add(_Guide_DGS->_position_absolute, tc2);
    tc1 = coords_sub(_Guide_DGS->_position_absolute, _Sample_pos->_position_absolute);
    _Sample_pos->_position_relative = rot_apply(_Sample_pos->_rotation_absolute, tc1);
  } /* Sample_pos=Arm() AT ROTATED */
  DEBUG_COMPONENT("Sample_pos", _Sample_pos->_position_absolute, _Sample_pos->_rotation_absolute);
  instrument->_position_absolute[33] = _Sample_pos->_position_absolute;
  instrument->_position_relative[33] = _Sample_pos->_position_relative;
  instrument->counter_N[33]  = instrument->counter_P[33] = instrument->counter_P2[33] = 0;
  instrument->counter_AbsorbProp[33]= 0;
  return(0);
} /* _Sample_pos_setpos */

/* component Sample=Isotropic_Sqw() SETTING, POSITION/ROTATION */
int _Sample_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Sample_setpos] component Sample=Isotropic_Sqw() SETTING [/usr/share/mcstas/3.0-dev/samples/Isotropic_Sqw.comp:2015]");
  stracpy(_Sample->_name, "Sample", 16384);
  stracpy(_Sample->_type, "Isotropic_Sqw", 16384);
  _Sample->_index=34;
  _Sample->_parameters.powder_format = (double []){ 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 }; // static pointer allocation
  #define powder_format (_Sample->_parameters.powder_format)
  if(instrument->_parameters._coh && strlen(instrument->_parameters._coh))
    stracpy(_Sample->_parameters.Sqw_coh, instrument->_parameters._coh ? instrument->_parameters._coh : "", 16384);
  else 
  _Sample->_parameters.Sqw_coh[0]='\0';
  #define Sqw_coh (_Sample->_parameters.Sqw_coh)
  if(instrument->_parameters._inc && strlen(instrument->_parameters._inc))
    stracpy(_Sample->_parameters.Sqw_inc, instrument->_parameters._inc ? instrument->_parameters._inc : "", 16384);
  else 
  _Sample->_parameters.Sqw_inc[0]='\0';
  #define Sqw_inc (_Sample->_parameters.Sqw_inc)
  _Sample->_parameters.geometry[0]='\0';
  #define geometry (_Sample->_parameters.geometry)
  _Sample->_parameters.radius = 0;
  #define radius (_Sample->_parameters.radius)
  _Sample->_parameters.thickness = 0;
  #define thickness (_Sample->_parameters.thickness)
  _Sample->_parameters.xwidth = 0.05;
  #define xwidth (_Sample->_parameters.xwidth)
  _Sample->_parameters.yheight = 0.06;
  #define yheight (_Sample->_parameters.yheight)
  _Sample->_parameters.zdepth = 0.002;
  #define zdepth (_Sample->_parameters.zdepth)
  _Sample->_parameters.threshold = 1e-20;
  #define threshold (_Sample->_parameters.threshold)
  _Sample->_parameters.order = 0;
  #define order (_Sample->_parameters.order)
  _Sample->_parameters.T = 0;
  #define T (_Sample->_parameters.T)
  _Sample->_parameters.verbose = 1;
  #define verbose (_Sample->_parameters.verbose)
  _Sample->_parameters.d_phi = 0;
  #define d_phi (_Sample->_parameters.d_phi)
  _Sample->_parameters.concentric = 0;
  #define concentric (_Sample->_parameters.concentric)
  _Sample->_parameters.rho = 0;
  #define rho (_Sample->_parameters.rho)
  _Sample->_parameters.sigma_abs = 0;
  #define sigma_abs (_Sample->_parameters.sigma_abs)
  _Sample->_parameters.sigma_coh = 0;
  #define sigma_coh (_Sample->_parameters.sigma_coh)
  _Sample->_parameters.sigma_inc = 0;
  #define sigma_inc (_Sample->_parameters.sigma_inc)
  _Sample->_parameters.classical = -1;
  #define classical (_Sample->_parameters.classical)
  _Sample->_parameters.powder_Dd = 0;
  #define powder_Dd (_Sample->_parameters.powder_Dd)
  _Sample->_parameters.powder_DW = 0;
  #define powder_DW (_Sample->_parameters.powder_DW)
  _Sample->_parameters.powder_Vc = 0;
  #define powder_Vc (_Sample->_parameters.powder_Vc)
  _Sample->_parameters.density = 0;
  #define density (_Sample->_parameters.density)
  _Sample->_parameters.weight = 0;
  #define weight (_Sample->_parameters.weight)
  _Sample->_parameters.p_interact = 0.99;
  #define p_interact (_Sample->_parameters.p_interact)
  _Sample->_parameters.norm = -1;
  #define norm (_Sample->_parameters.norm)
  _Sample->_parameters.powder_barns = 1;
  #define powder_barns (_Sample->_parameters.powder_barns)
  if("Frommhold" && strlen("Frommhold"))
    stracpy(_Sample->_parameters.quantum_correction, "Frommhold" ? "Frommhold" : "", 16384);
  else 
  _Sample->_parameters.quantum_correction[0]='\0';
  #define quantum_correction (_Sample->_parameters.quantum_correction)

  #define VarSqw (_Sample->_parameters.VarSqw)
  #define columns (_Sample->_parameters.columns)
  #define offdata (_Sample->_parameters.offdata)

  #undef powder_format
  #undef Sqw_coh
  #undef Sqw_inc
  #undef geometry
  #undef radius
  #undef thickness
  #undef xwidth
  #undef yheight
  #undef zdepth
  #undef threshold
  #undef order
  #undef T
  #undef verbose
  #undef d_phi
  #undef concentric
  #undef rho
  #undef sigma_abs
  #undef sigma_coh
  #undef sigma_inc
  #undef classical
  #undef powder_Dd
  #undef powder_DW
  #undef powder_Vc
  #undef density
  #undef weight
  #undef p_interact
  #undef norm
  #undef powder_barns
  #undef quantum_correction
  #undef VarSqw
  #undef columns
  #undef offdata
  /* component Sample=Isotropic_Sqw() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (45)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _Sample_pos->_rotation_absolute, _Sample->_rotation_absolute);
    rot_transpose(_Sample_pos->_rotation_absolute, tr1);
    rot_mul(_Sample->_rotation_absolute, tr1, _Sample->_rotation_relative);
    _Sample->_rotation_is_identity =  rot_test_identity(_Sample->_rotation_relative);
    tc1 = coords_set(
      0, 0, 0);
    rot_transpose(_Sample_pos->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Sample->_position_absolute = coords_add(_Sample_pos->_position_absolute, tc2);
    tc1 = coords_sub(_Sample_pos->_position_absolute, _Sample->_position_absolute);
    _Sample->_position_relative = rot_apply(_Sample->_rotation_absolute, tc1);
  } /* Sample=Isotropic_Sqw() AT ROTATED */
  DEBUG_COMPONENT("Sample", _Sample->_position_absolute, _Sample->_rotation_absolute);
  instrument->_position_absolute[34] = _Sample->_position_absolute;
  instrument->_position_relative[34] = _Sample->_position_relative;
  instrument->counter_N[34]  = instrument->counter_P[34] = instrument->counter_P2[34] = 0;
  instrument->counter_AbsorbProp[34]= 0;
  return(0);
} /* _Sample_setpos */

/* component Ideal_Det=Monitor_nD() SETTING, POSITION/ROTATION */
int _Ideal_Det_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Ideal_Det_setpos] component Ideal_Det=Monitor_nD() SETTING [/usr/share/mcstas/3.0-dev/monitors/Monitor_nD.comp:228]");
  stracpy(_Ideal_Det->_name, "Ideal_Det", 16384);
  stracpy(_Ideal_Det->_type, "Monitor_nD", 16384);
  _Ideal_Det->_index=35;
  _Ideal_Det->_parameters.user1 = '\0';
  #define user1 (_Ideal_Det->_parameters.user1)
  _Ideal_Det->_parameters.user2 = '\0';
  #define user2 (_Ideal_Det->_parameters.user2)
  _Ideal_Det->_parameters.user3 = '\0';
  #define user3 (_Ideal_Det->_parameters.user3)
  _Ideal_Det->_parameters.xwidth = 0;
  #define xwidth (_Ideal_Det->_parameters.xwidth)
  _Ideal_Det->_parameters.yheight = 2;
  #define yheight (_Ideal_Det->_parameters.yheight)
  _Ideal_Det->_parameters.zdepth = 0;
  #define zdepth (_Ideal_Det->_parameters.zdepth)
  _Ideal_Det->_parameters.xmin = 0;
  #define xmin (_Ideal_Det->_parameters.xmin)
  _Ideal_Det->_parameters.xmax = 0;
  #define xmax (_Ideal_Det->_parameters.xmax)
  _Ideal_Det->_parameters.ymin = 0;
  #define ymin (_Ideal_Det->_parameters.ymin)
  _Ideal_Det->_parameters.ymax = 0;
  #define ymax (_Ideal_Det->_parameters.ymax)
  _Ideal_Det->_parameters.zmin = 0;
  #define zmin (_Ideal_Det->_parameters.zmin)
  _Ideal_Det->_parameters.zmax = 0;
  #define zmax (_Ideal_Det->_parameters.zmax)
  _Ideal_Det->_parameters.bins = 0;
  #define bins (_Ideal_Det->_parameters.bins)
  _Ideal_Det->_parameters.min = -1e40;
  #define min (_Ideal_Det->_parameters.min)
  _Ideal_Det->_parameters.max = 1e40;
  #define max (_Ideal_Det->_parameters.max)
  _Ideal_Det->_parameters.restore_neutron = 1;
  #define restore_neutron (_Ideal_Det->_parameters.restore_neutron)
  _Ideal_Det->_parameters.radius = 2.5;
  #define radius (_Ideal_Det->_parameters.radius)
  if(detector_options && strlen(detector_options))
    stracpy(_Ideal_Det->_parameters.options, detector_options ? detector_options : "", 16384);
  else 
  _Ideal_Det->_parameters.options[0]='\0';
  #define options (_Ideal_Det->_parameters.options)
  if("NULL" && strlen("NULL"))
    stracpy(_Ideal_Det->_parameters.filename, "NULL" ? "NULL" : "", 16384);
  else 
  _Ideal_Det->_parameters.filename[0]='\0';
  #define filename (_Ideal_Det->_parameters.filename)
  if("NULL" && strlen("NULL"))
    stracpy(_Ideal_Det->_parameters.geometry, "NULL" ? "NULL" : "", 16384);
  else 
  _Ideal_Det->_parameters.geometry[0]='\0';
  #define geometry (_Ideal_Det->_parameters.geometry)
  if("NULL" && strlen("NULL"))
    stracpy(_Ideal_Det->_parameters.username1, "NULL" ? "NULL" : "", 16384);
  else 
  _Ideal_Det->_parameters.username1[0]='\0';
  #define username1 (_Ideal_Det->_parameters.username1)
  if("NULL" && strlen("NULL"))
    stracpy(_Ideal_Det->_parameters.username2, "NULL" ? "NULL" : "", 16384);
  else 
  _Ideal_Det->_parameters.username2[0]='\0';
  #define username2 (_Ideal_Det->_parameters.username2)
  if("NULL" && strlen("NULL"))
    stracpy(_Ideal_Det->_parameters.username3, "NULL" ? "NULL" : "", 16384);
  else 
  _Ideal_Det->_parameters.username3[0]='\0';
  #define username3 (_Ideal_Det->_parameters.username3)

  #define DEFS (_Ideal_Det->_parameters.DEFS)
  #define Vars (_Ideal_Det->_parameters.Vars)
  #define detector (_Ideal_Det->_parameters.detector)
  #define offdata (_Ideal_Det->_parameters.offdata)

  #undef user1
  #undef user2
  #undef user3
  #undef xwidth
  #undef yheight
  #undef zdepth
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef zmin
  #undef zmax
  #undef bins
  #undef min
  #undef max
  #undef restore_neutron
  #undef radius
  #undef options
  #undef filename
  #undef geometry
  #undef username1
  #undef username2
  #undef username3
  #undef DEFS
  #undef Vars
  #undef detector
  #undef offdata
  /* component Ideal_Det=Monitor_nD() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Sample_pos->_rotation_absolute, _Ideal_Det->_rotation_absolute);
    rot_transpose(_Sample->_rotation_absolute, tr1);
    rot_mul(_Ideal_Det->_rotation_absolute, tr1, _Ideal_Det->_rotation_relative);
    _Ideal_Det->_rotation_is_identity =  rot_test_identity(_Ideal_Det->_rotation_relative);
    tc1 = coords_set(
      0, 0, 0);
    rot_transpose(_Sample_pos->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Ideal_Det->_position_absolute = coords_add(_Sample_pos->_position_absolute, tc2);
    tc1 = coords_sub(_Sample->_position_absolute, _Ideal_Det->_position_absolute);
    _Ideal_Det->_position_relative = rot_apply(_Ideal_Det->_rotation_absolute, tc1);
  } /* Ideal_Det=Monitor_nD() AT ROTATED */
  DEBUG_COMPONENT("Ideal_Det", _Ideal_Det->_position_absolute, _Ideal_Det->_rotation_absolute);
  instrument->_position_absolute[35] = _Ideal_Det->_position_absolute;
  instrument->_position_relative[35] = _Ideal_Det->_position_relative;
  instrument->counter_N[35]  = instrument->counter_P[35] = instrument->counter_P2[35] = 0;
  instrument->counter_AbsorbProp[35]= 0;
  return(0);
} /* _Ideal_Det_setpos */

/* component Detector=PSD_Detector() SETTING, POSITION/ROTATION */
int _Detector_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Detector_setpos] component Detector=PSD_Detector() SETTING [/usr/share/mcstas/3.0-dev/contrib/PSD_Detector.comp:258]");
  stracpy(_Detector->_name, "Detector", 16384);
  stracpy(_Detector->_type, "PSD_Detector", 16384);
  _Detector->_index=36;
  _Detector->_parameters.nx = 640;
  #define nx (_Detector->_parameters.nx)
  _Detector->_parameters.ny = 256;
  #define ny (_Detector->_parameters.ny)
  _Detector->_parameters.xwidth = 0;
  #define xwidth (_Detector->_parameters.xwidth)
  _Detector->_parameters.radius = 2.5;
  #define radius (_Detector->_parameters.radius)
  _Detector->_parameters.awidth = 0;
  #define awidth (_Detector->_parameters.awidth)
  _Detector->_parameters.yheight = 2.0;
  #define yheight (_Detector->_parameters.yheight)
  _Detector->_parameters.zdepth = 0.01;
  #define zdepth (_Detector->_parameters.zdepth)
  _Detector->_parameters.threshold = 100;
  #define threshold (_Detector->_parameters.threshold)
  _Detector->_parameters.PressureConv = 8;
  #define PressureConv (_Detector->_parameters.PressureConv)
  _Detector->_parameters.PressureStop = 1;
  #define PressureStop (_Detector->_parameters.PressureStop)
  _Detector->_parameters.interpolate = 1;
  #define interpolate (_Detector->_parameters.interpolate)
  _Detector->_parameters.p_interact = 0;
  #define p_interact (_Detector->_parameters.p_interact)
  _Detector->_parameters.verbose = 0;
  #define verbose (_Detector->_parameters.verbose)
  _Detector->_parameters.LensOn = 1;
  #define LensOn (_Detector->_parameters.LensOn)
  _Detector->_parameters.dc = 0;
  #define dc (_Detector->_parameters.dc)
  _Detector->_parameters.borderx = -1;
  #define borderx (_Detector->_parameters.borderx)
  _Detector->_parameters.bordery = -1;
  #define bordery (_Detector->_parameters.bordery)
  _Detector->_parameters.xChDivRelSigma = 0;
  #define xChDivRelSigma (_Detector->_parameters.xChDivRelSigma)
  _Detector->_parameters.yChDivRelSigma = 0.0037;
  #define yChDivRelSigma (_Detector->_parameters.yChDivRelSigma)
  _Detector->_parameters.bufsize = 0;
  #define bufsize (_Detector->_parameters.bufsize)
  _Detector->_parameters.restore_neutron = 0;
  #define restore_neutron (_Detector->_parameters.restore_neutron)
  _Detector->_parameters.angle = 136 -13;
  #define angle (_Detector->_parameters.angle)
  _Detector->_parameters.type[0]='\0';
  #define type (_Detector->_parameters.type)
  if("NEAT.psd" && strlen("NEAT.psd"))
    stracpy(_Detector->_parameters.filename, "NEAT.psd" ? "NEAT.psd" : "", 16384);
  else 
  _Detector->_parameters.filename[0]='\0';
  #define filename (_Detector->_parameters.filename)
  if("Gas_tables/He3inHe.table" && strlen("Gas_tables/He3inHe.table"))
    stracpy(_Detector->_parameters.FN_Conv, "Gas_tables/He3inHe.table" ? "Gas_tables/He3inHe.table" : "", 16384);
  else 
  _Detector->_parameters.FN_Conv[0]='\0';
  #define FN_Conv (_Detector->_parameters.FN_Conv)
  if("Gas_tables/He3inCF4.table" && strlen("Gas_tables/He3inCF4.table"))
    stracpy(_Detector->_parameters.FN_Stop, "Gas_tables/He3inCF4.table" ? "Gas_tables/He3inCF4.table" : "", 16384);
  else 
  _Detector->_parameters.FN_Stop[0]='\0';
  #define FN_Stop (_Detector->_parameters.FN_Stop)

  #define PSD_N (_Detector->_parameters.PSD_N)
  #define PSD_p (_Detector->_parameters.PSD_p)
  #define PSD_p2 (_Detector->_parameters.PSD_p2)
  #define EAP (_Detector->_parameters.EAP)
  #define M1P1 (_Detector->_parameters.M1P1)
  #define PosAP (_Detector->_parameters.PosAP)
  #define EAT (_Detector->_parameters.EAT)
  #define M1T1 (_Detector->_parameters.M1T1)
  #define PosAT (_Detector->_parameters.PosAT)
  #define CrossSectionHe (_Detector->_parameters.CrossSectionHe)
  #define CountNeutrons (_Detector->_parameters.CountNeutrons)
  #define GeomCumul (_Detector->_parameters.GeomCumul)
  #define AbsCumul (_Detector->_parameters.AbsCumul)
  #define SensVolCumul (_Detector->_parameters.SensVolCumul)
  #define DetCumul (_Detector->_parameters.DetCumul)
  #define PHSpectrum (_Detector->_parameters.PHSpectrum)
  #define PHSpectrum0 (_Detector->_parameters.PHSpectrum0)
  #define PHSpectrum2 (_Detector->_parameters.PHSpectrum2)
  #define PHSpectrum_n (_Detector->_parameters.PHSpectrum_n)
  #define nH_p (_Detector->_parameters.nH_p)
  #define nH_t (_Detector->_parameters.nH_t)
  #define FullEnergyP (_Detector->_parameters.FullEnergyP)
  #define FullEnergyT (_Detector->_parameters.FullEnergyT)
  #define VariousErrors (_Detector->_parameters.VariousErrors)
  #define DetectorType (_Detector->_parameters.DetectorType)
  #define DEFS (_Detector->_parameters.DEFS)
  #define Vars (_Detector->_parameters.Vars)

  #undef nx
  #undef ny
  #undef xwidth
  #undef radius
  #undef awidth
  #undef yheight
  #undef zdepth
  #undef threshold
  #undef PressureConv
  #undef PressureStop
  #undef interpolate
  #undef p_interact
  #undef verbose
  #undef LensOn
  #undef dc
  #undef borderx
  #undef bordery
  #undef xChDivRelSigma
  #undef yChDivRelSigma
  #undef bufsize
  #undef restore_neutron
  #undef angle
  #undef type
  #undef filename
  #undef FN_Conv
  #undef FN_Stop
  #undef PSD_N
  #undef PSD_p
  #undef PSD_p2
  #undef EAP
  #undef M1P1
  #undef PosAP
  #undef EAT
  #undef M1T1
  #undef PosAT
  #undef CrossSectionHe
  #undef CountNeutrons
  #undef GeomCumul
  #undef AbsCumul
  #undef SensVolCumul
  #undef DetCumul
  #undef PHSpectrum
  #undef PHSpectrum0
  #undef PHSpectrum2
  #undef PHSpectrum_n
  #undef nH_p
  #undef nH_t
  #undef FullEnergyP
  #undef FullEnergyT
  #undef VariousErrors
  #undef DetectorType
  #undef DEFS
  #undef Vars
  /* component Detector=PSD_Detector() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (( 136 -13 ) / 2 + 13)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _Sample_pos->_rotation_absolute, _Detector->_rotation_absolute);
    rot_transpose(_Ideal_Det->_rotation_absolute, tr1);
    rot_mul(_Detector->_rotation_absolute, tr1, _Detector->_rotation_relative);
    _Detector->_rotation_is_identity =  rot_test_identity(_Detector->_rotation_relative);
    tc1 = coords_set(
      0, 0, 0);
    rot_transpose(_Sample_pos->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Detector->_position_absolute = coords_add(_Sample_pos->_position_absolute, tc2);
    tc1 = coords_sub(_Ideal_Det->_position_absolute, _Detector->_position_absolute);
    _Detector->_position_relative = rot_apply(_Detector->_rotation_absolute, tc1);
  } /* Detector=PSD_Detector() AT ROTATED */
  DEBUG_COMPONENT("Detector", _Detector->_position_absolute, _Detector->_rotation_absolute);
  instrument->_position_absolute[36] = _Detector->_position_absolute;
  instrument->_position_relative[36] = _Detector->_position_relative;
  instrument->counter_N[36]  = instrument->counter_P[36] = instrument->counter_P2[36] = 0;
  instrument->counter_AbsorbProp[36]= 0;
  return(0);
} /* _Detector_setpos */

_class_Progress_bar *class_Progress_bar_init(_class_Progress_bar *_comp
) {
  #define profile (_comp->_parameters.profile)
  #define percent (_comp->_parameters.percent)
  #define flag_save (_comp->_parameters.flag_save)
  #define minutes (_comp->_parameters.minutes)
  #define IntermediateCnts (_comp->_parameters.IntermediateCnts)
  #define StartTime (_comp->_parameters.StartTime)
  #define EndTime (_comp->_parameters.EndTime)
  #define CurrentTime (_comp->_parameters.CurrentTime)
IntermediateCnts=0;
StartTime=0;
EndTime=0;
CurrentTime=0;

fprintf(stdout, "[%s] Initialize\n", instrument_name);
  if (percent*mcget_ncount()/100 < 1e5) {
    percent=1e5*100.0/mcget_ncount();
  }
  #undef profile
  #undef percent
  #undef flag_save
  #undef minutes
  #undef IntermediateCnts
  #undef StartTime
  #undef EndTime
  #undef CurrentTime
  return(_comp);
} /* class_Progress_bar_init */

_class_Source_gen *class_Source_gen_init(_class_Source_gen *_comp
) {
  #define flux_file (_comp->_parameters.flux_file)
  #define xdiv_file (_comp->_parameters.xdiv_file)
  #define ydiv_file (_comp->_parameters.ydiv_file)
  #define radius (_comp->_parameters.radius)
  #define dist (_comp->_parameters.dist)
  #define focus_xw (_comp->_parameters.focus_xw)
  #define focus_yh (_comp->_parameters.focus_yh)
  #define focus_aw (_comp->_parameters.focus_aw)
  #define focus_ah (_comp->_parameters.focus_ah)
  #define E0 (_comp->_parameters.E0)
  #define dE (_comp->_parameters.dE)
  #define lambda0 (_comp->_parameters.lambda0)
  #define dlambda (_comp->_parameters.dlambda)
  #define I1 (_comp->_parameters.I1)
  #define yheight (_comp->_parameters.yheight)
  #define xwidth (_comp->_parameters.xwidth)
  #define verbose (_comp->_parameters.verbose)
  #define T1 (_comp->_parameters.T1)
  #define flux_file_perAA (_comp->_parameters.flux_file_perAA)
  #define flux_file_log (_comp->_parameters.flux_file_log)
  #define Lmin (_comp->_parameters.Lmin)
  #define Lmax (_comp->_parameters.Lmax)
  #define Emin (_comp->_parameters.Emin)
  #define Emax (_comp->_parameters.Emax)
  #define T2 (_comp->_parameters.T2)
  #define I2 (_comp->_parameters.I2)
  #define T3 (_comp->_parameters.T3)
  #define I3 (_comp->_parameters.I3)
  #define zdepth (_comp->_parameters.zdepth)
  #define target_index (_comp->_parameters.target_index)
  #define p_in (_comp->_parameters.p_in)
  #define lambda1 (_comp->_parameters.lambda1)
  #define lambda2 (_comp->_parameters.lambda2)
  #define lambda3 (_comp->_parameters.lambda3)
  #define pTable (_comp->_parameters.pTable)
  #define pTable_x (_comp->_parameters.pTable_x)
  #define pTable_y (_comp->_parameters.pTable_y)
  #define pTable_xmin (_comp->_parameters.pTable_xmin)
  #define pTable_xmax (_comp->_parameters.pTable_xmax)
  #define pTable_xsum (_comp->_parameters.pTable_xsum)
  #define pTable_ymin (_comp->_parameters.pTable_ymin)
  #define pTable_ymax (_comp->_parameters.pTable_ymax)
  #define pTable_ysum (_comp->_parameters.pTable_ysum)
  #define pTable_dxmin (_comp->_parameters.pTable_dxmin)
  #define pTable_dxmax (_comp->_parameters.pTable_dxmax)
  #define pTable_dymin (_comp->_parameters.pTable_dymin)
  #define pTable_dymax (_comp->_parameters.pTable_dymax)
  pTable_xsum=0;
  pTable_ysum=0;


  double source_area, k;

  if (target_index && !dist)
  {
    Coords ToTarget;
    double tx,ty,tz;
    ToTarget = coords_sub(POS_A_COMP_INDEX(INDEX_CURRENT_COMP+target_index),POS_A_CURRENT_COMP);
    ToTarget = rot_apply(ROT_A_CURRENT_COMP, ToTarget);
    coords_get(ToTarget, &tx, &ty, &tz);
    dist=sqrt(tx*tx+ty*ty+tz*tz);
  }

  /* spectrum characteristics */
  if (flux_file && strlen(flux_file) && strcmp(flux_file,"NULL") && strcmp(flux_file,"0")) {
    if (Table_Read(&pTable, flux_file, 1) <= 0) /* read 1st block data from file into pTable */
      exit(fprintf(stderr, "Source_gen: %s: can not read flux file %s\n", NAME_CURRENT_COMP, flux_file));
    /* put table in Log scale */
    int i;
    if (pTable.columns < 2) exit(fprintf(stderr, "Source_gen: %s: Flux file %s should contain at least 2 columns [wavelength in Angs,flux].\n", NAME_CURRENT_COMP, flux_file));
    double table_lmin=FLT_MAX, table_lmax=-FLT_MAX;
    double tmin=FLT_MAX, tmax=-FLT_MAX;
    for (i=0; i<pTable.rows; i++) {
      double val = Table_Index(pTable, i,1);
      val = Table_Index(pTable, i,0); /* lambda */
      if (val > tmax) tmax=val;
      if (val < tmin) tmin=val;
    }
    for (i=0; i<pTable.rows; i++) {
      double val = Table_Index(pTable, i,1);
      if (val < 0) fprintf(stderr, "Source_gen: %s: File %s has negative flux at row %i.\n", NAME_CURRENT_COMP, flux_file, i+1);
      if (flux_file_log)
        val = log(val > 0 ? val : tmin/10);
      Table_SetElement(&pTable, i, 1, val);
      val = Table_Index(pTable, i,0); /* lambda */
      if (val > table_lmax) table_lmax=val;
      if (val < table_lmin) table_lmin=val;
    }
    if (!Lmin && !Lmax && !lambda0 && !dlambda && !E0 && !dE && !Emin && !Emax) {
      Lmin = table_lmin; Lmax = table_lmax;
    }
    if (Lmax > table_lmax) {
      if (verbose) fprintf(stderr, "Source_gen: %s: Maximum wavelength %g is beyond table range upper limit %g. Constraining.\n", NAME_CURRENT_COMP, Lmax, table_lmax);
      Lmax = table_lmax;
    }
    if (Lmin < table_lmin) {
      if (verbose) fprintf(stderr, "Source_gen: %s: Minimum wavelength %g is below table range lower limit %g. Constraining.\n", NAME_CURRENT_COMP, Lmin, table_lmin);
      Lmin = table_lmin;
    }
  }  /* end flux file */
  else
  {
    k  = 1.38066e-23; /* k_B */
    if (T1 > 0)
    {
      lambda1  = 1.0e10*sqrt(HBAR*HBAR*4.0*PI*PI/2.0/MNEUTRON/k/T1);
    }
    else
      { lambda1 = lambda0; }

    if (T2 > 0)
    {
      lambda2  = 1.0e10*sqrt(HBAR*HBAR*4.0*PI*PI/2.0/MNEUTRON/k/T2);
    }
    else
      { lambda2 = lambda0; }

    if (T3 > 0)
    {
      lambda3  = 1.0e10*sqrt(HBAR*HBAR*4.0*PI*PI/2.0/MNEUTRON/k/T3);
    }
    else
      { lambda3 = lambda0; }
  }

  /* now read position-divergence files, if any */
  if (xdiv_file && strlen(xdiv_file) && strcmp(xdiv_file,"NULL") && strcmp(xdiv_file,"0")) {
    int i,j;
    if (Table_Read(&pTable_x, xdiv_file, 1) <= 0) /* read 1st block data from file into pTable */
      exit(fprintf(stderr, "Source_gen: %s: can not read XDiv file %s\n", NAME_CURRENT_COMP, xdiv_file));
    pTable_xsum = 0;
    for (i=0; i<pTable_x.rows; i++)
      for (j=0; j<pTable_x.columns; j++)
        pTable_xsum += Table_Index(pTable_x, i,j);

    /* now extract limits */
    char **parsing;
    char xylimits[1024];
    strcpy(xylimits, "");
    parsing = Table_ParseHeader(pTable_x.header,
      "xlimits", "xylimits",
      NULL);

    if (parsing) {
      if (parsing[0])  strcpy(xylimits, str_dup_numeric(parsing[0]));
      if (parsing[1] && !strlen(xylimits))
                       strcpy(xylimits, str_dup_numeric(parsing[1]));
      for (i=0; i<=1; i++) {
        if (parsing[i]) free(parsing[i]);
      }
      free(parsing);
    }
    i = sscanf(xylimits, "%lg %lg %lg %lg",
      &(pTable_xmin),  &(pTable_xmax),
      &(pTable_dxmin), &(pTable_dxmax));
    if (i != 2 && i != 4 && verbose)
      fprintf(stderr, "Source_gen: %s: invalid xylimits '%s' from file %s. extracted %i values\n",
        NAME_CURRENT_COMP, xylimits, xdiv_file, i);

    if (!xwidth) xwidth=pTable_xmax-pTable_xmin;
    if (!focus_xw && !dist) focus_xw=fabs(pTable_dxmax-pTable_dxmin);
  } /* end xdiv file */

  if (ydiv_file && strlen(ydiv_file) && strcmp(ydiv_file,"NULL") && strcmp(ydiv_file,"0")) {
    int i,j;
    if (Table_Read(&pTable_y, ydiv_file, 1) <= 0) /* read 1st block data from file into pTable */
      exit(fprintf(stderr, "Source_gen: %s: can not read YDiv file %s\n", NAME_CURRENT_COMP, ydiv_file));
    pTable_ysum = 0;
    for (i=0; i<pTable_y.rows; i++)
      for (j=0; j<pTable_y.columns; j++)
        pTable_ysum += Table_Index(pTable_y, i,j);

    /* now extract limits */
    char **parsing;
    char xylimits[1024];
    strcpy(xylimits, "");
    parsing = Table_ParseHeader(pTable_y.header,
      "xlimits", "xylimits",
      NULL);

    if (parsing) {
      if (parsing[0])  strcpy(xylimits,str_dup_numeric(parsing[0]));
      if (parsing[1] && !strlen(xylimits))
                       strcpy(xylimits,str_dup_numeric(parsing[1]));
      for (i=0; i<=1; i++) {
        if (parsing[i]) free(parsing[i]);
      }
      free(parsing);
    }
    i = sscanf(xylimits, "%lg %lg %lg %lg",
      &(pTable_ymin),  &(pTable_ymax),
      &(pTable_dymin), &(pTable_dymax));
    if (i != 2 && i != 4 && verbose)
      fprintf(stderr, "Source_gen: %s: invalid xylimits '%s' from file %s. extracted %i values\n",
        NAME_CURRENT_COMP, xylimits, ydiv_file, i);
    if (!yheight)  yheight=pTable_ymax-pTable_ymin;
    if (!focus_yh && !dist) focus_yh=fabs(pTable_dymax-pTable_dymin);
  } /* end ydiv file */

  /* tests for parameter values */
  if (Emin < 0 || Emax < 0 || Lmin < 0 || Lmax < 0 || E0 < 0 || dE < 0 || lambda0 < 0 || dlambda < 0)
  {
    fprintf(stderr,"Source_gen: %s: Error: Negative average\n"
                   "            or range values for wavelength or energy encountered\n",
                   NAME_CURRENT_COMP);
    exit(-1);
  }
  if ((Emin == 0 && Emax > 0) || (dE > 0 && dE >= E0))
  {
    fprintf(stderr,"Source_gen: %s: Error: minimal energy cannot be less or equal zero\n",
      NAME_CURRENT_COMP);
    exit(-1);
  }
  if ((Emax >= Emin) && (Emin > 0))
  { E0 = (Emax+Emin)/2;
    dE = (Emax-Emin)/2;
  }
  if ((E0 > dE) && (dE >= 0))
  {
    Lmin = sqrt(81.81/(E0+dE)); /* Angstroem */
    Lmax = sqrt(81.81/(E0-dE));
  }
  if (Lmax > 0)
  { lambda0 = (Lmax+Lmin)/2;
    dlambda = (Lmax-Lmin)/2;
  }
  if (lambda0 <= 0 || (lambda0 < dlambda) || (dlambda < 0))
  { fprintf(stderr,"Source_gen: %s: Error: Wavelength range %.3f +/- %.3f AA calculated \n",
      NAME_CURRENT_COMP, lambda0, dlambda);
    fprintf(stderr,"- whole wavelength range must be >= 0 \n");
    fprintf(stderr,"- range must be > 0; otherwise intensity gets zero, use other sources in this case \n\n");
    exit(-1);
  }

  radius = fabs(radius); xwidth=fabs(xwidth); yheight=fabs(yheight);  I1=fabs(I1);
  lambda0=fabs(lambda0); dlambda=fabs(dlambda);
  focus_xw = fabs(focus_xw); focus_yh=fabs(focus_yh); dist=fabs(dist);

  if ((!focus_ah && !focus_aw) && (!focus_xw && !focus_yh))
  {
    fprintf(stderr,"Source_gen: %s: Error: No focusing information.\n"
                   "            Specify focus_xw, focus_yh or focus_aw, focus_ah\n",
                   NAME_CURRENT_COMP);
    exit(-1);
  }
  Lmin = lambda0 - dlambda; /* Angstroem */
  Lmax = lambda0 + dlambda;

  /* compute initial weight factor p_in to get [n/s] */
  if ((I1 > 0  && T1 >= 0)
     || (flux_file && strlen(flux_file) && strcmp(flux_file,"NULL") && strcmp(flux_file,"0")))
  { /* the I1,2,3 are usually in [n/s/cm2/st/AA] */
    if (radius)
      source_area = radius*radius*PI*1e4; /* circular cm^2 */
    else
      source_area = yheight*xwidth*1e4; /* square cm^2 */
    p_in  = source_area; /* cm2 */
    p_in *= (Lmax-Lmin); /* AA. 1 bin=AA/n */
    if (flux_file && strlen(flux_file) && strcmp(flux_file,"NULL") && strcmp(flux_file,"0")
      && !flux_file_perAA)  p_in *= pTable.rows/(Lmax-Lmin);
  }
  else
    p_in = 1.0/4/PI; /* Small angle approx. */
  p_in /= mcget_ncount();
  if (!T1 && I1) p_in *= I1;

  if (radius == 0 && yheight == 0 && xwidth == 0)
  {
    fprintf(stderr,"Source_gen: %s: Error: Please specify source geometry (radius, yheight, xwidth)\n",
      NAME_CURRENT_COMP);
    exit(-1);
  }
  if (focus_xw*focus_yh == 0)
  {
    fprintf(stderr,"Source_gen: %s: Error: Please specify source target (focus_xw, focus_yh)\n",
      NAME_CURRENT_COMP);
    exit(-1);
  }
  MPI_MASTER(
  if (verbose)
  {
    printf("Source_gen: component %s ", NAME_CURRENT_COMP);
    if ((yheight == 0) || (xwidth == 0))
      printf("(disk, radius=%g)", radius);
    else
      printf("(square %g x %g)",xwidth,yheight);
    if (dist) printf("\n            focusing distance dist=%g area=%g x %g\n", dist, focus_xw, focus_yh);
    printf("            spectra ");
    printf("%.3f to %.3f AA (%.3f to %.3f meV)", Lmin, Lmax, 81.81/Lmax/Lmax, 81.81/Lmin/Lmin);
    printf("\n");
    if (flux_file && strlen(flux_file) && strcmp(flux_file,"NULL") && strcmp(flux_file,"0"))
    { printf("  File %s for flux distribution used. Flux is dPhi/dlambda in [n/s/AA]. \n", flux_file);
      Table_Info(pTable);
    }
    else if (T1>=0 && I1)
    { if (T1 != 0)
        printf("            T1=%.1f K (%.3f AA)", T1, lambda1);
      if (T2*I2 != 0)
        printf(", T2=%.1f K (%.3f AA)", T2, lambda2);
      if (T3*I3 != 0)
        printf(", T3=%.1f K (%.3f AA)", T3, lambda3);
      if (T1) printf("\n");
      printf("  Flux is dPhi/dlambda in [n/s/cm2].\n");
    }
    else
    { printf("  Flux is Phi in [n/s].\n");
    }
    if (xdiv_file && strlen(xdiv_file) && strcmp(xdiv_file,"NULL") && strcmp(xdiv_file,"0"))
      printf("  File %s x=[%g:%g] [m] xdiv=[%g:%g] [deg] used as horizontal phase space distribution.\n", xdiv_file, pTable_xmin, pTable_xmax, pTable_dxmin, pTable_dxmax);
    if (ydiv_file && strlen(ydiv_file) && strcmp(ydiv_file,"NULL") && strcmp(ydiv_file,"0"))
      printf("  File %s y=[%g:%g] [m] ydiv=[%g:%g] [deg] used as vertical phase space distribution.\n", ydiv_file, pTable_ymin, pTable_ymax, pTable_dymin, pTable_dymax);
  }
  else
    if (verbose == -1)
      printf("Source_gen: component %s unactivated", NAME_CURRENT_COMP);
  );
  #undef flux_file
  #undef xdiv_file
  #undef ydiv_file
  #undef radius
  #undef dist
  #undef focus_xw
  #undef focus_yh
  #undef focus_aw
  #undef focus_ah
  #undef E0
  #undef dE
  #undef lambda0
  #undef dlambda
  #undef I1
  #undef yheight
  #undef xwidth
  #undef verbose
  #undef T1
  #undef flux_file_perAA
  #undef flux_file_log
  #undef Lmin
  #undef Lmax
  #undef Emin
  #undef Emax
  #undef T2
  #undef I2
  #undef T3
  #undef I3
  #undef zdepth
  #undef target_index
  #undef p_in
  #undef lambda1
  #undef lambda2
  #undef lambda3
  #undef pTable
  #undef pTable_x
  #undef pTable_y
  #undef pTable_xmin
  #undef pTable_xmax
  #undef pTable_xsum
  #undef pTable_ymin
  #undef pTable_ymax
  #undef pTable_ysum
  #undef pTable_dxmin
  #undef pTable_dxmax
  #undef pTable_dymin
  #undef pTable_dymax
  return(_comp);
} /* class_Source_gen_init */

_class_Guide_gravity *class_Guide_gravity_init(_class_Guide_gravity *_comp
) {
  #define w1 (_comp->_parameters.w1)
  #define h1 (_comp->_parameters.h1)
  #define w2 (_comp->_parameters.w2)
  #define h2 (_comp->_parameters.h2)
  #define l (_comp->_parameters.l)
  #define R0 (_comp->_parameters.R0)
  #define Qc (_comp->_parameters.Qc)
  #define alpha (_comp->_parameters.alpha)
  #define m (_comp->_parameters.m)
  #define W (_comp->_parameters.W)
  #define nslit (_comp->_parameters.nslit)
  #define d (_comp->_parameters.d)
  #define mleft (_comp->_parameters.mleft)
  #define mright (_comp->_parameters.mright)
  #define mtop (_comp->_parameters.mtop)
  #define mbottom (_comp->_parameters.mbottom)
  #define nhslit (_comp->_parameters.nhslit)
  #define G (_comp->_parameters.G)
  #define aleft (_comp->_parameters.aleft)
  #define aright (_comp->_parameters.aright)
  #define atop (_comp->_parameters.atop)
  #define abottom (_comp->_parameters.abottom)
  #define wavy (_comp->_parameters.wavy)
  #define wavy_z (_comp->_parameters.wavy_z)
  #define wavy_tb (_comp->_parameters.wavy_tb)
  #define wavy_lr (_comp->_parameters.wavy_lr)
  #define chamfers (_comp->_parameters.chamfers)
  #define chamfers_z (_comp->_parameters.chamfers_z)
  #define chamfers_lr (_comp->_parameters.chamfers_lr)
  #define chamfers_tb (_comp->_parameters.chamfers_tb)
  #define nelements (_comp->_parameters.nelements)
  #define nu (_comp->_parameters.nu)
  #define phase (_comp->_parameters.phase)
  #define reflect (_comp->_parameters.reflect)
  #define GVars (_comp->_parameters.GVars)
  #define pTable (_comp->_parameters.pTable)
  double Gx=0, Gy=-GRAVITY, Gz=0;
  Coords mcLocG;
  int i;

  if (reflect && strlen(reflect) && strcmp(reflect,"NULL") && strcmp(reflect,"0")) {
    if (Table_Read(&pTable, reflect, 1) <= 0) /* read 1st block data from file into pTable */
      exit(fprintf(stderr,"Guide_gravity: %s: can not read file %s\n", NAME_CURRENT_COMP, reflect));
  } else {
    if (W < 0 || R0 < 0 || Qc < 0)
    { fprintf(stderr,"Guide_gravity: %s: W R0 Qc must be >0.\n", NAME_CURRENT_COMP);
      exit(-1); }
  }

  if (nslit <= 0 || nhslit <= 0)
  { fprintf(stderr,"Guide_gravity: %s: nslit nhslit must be >0.\n", NAME_CURRENT_COMP);
    exit(-1); }

  if (!w1 || !h1)
  { fprintf(stderr,"Guide_gravity: %s: input window is closed (w1=h1=0).\n", NAME_CURRENT_COMP);
    exit(-1); }

  if (d*nslit > w1) exit(fprintf(stderr, "Guide_gravity: %s: absorbing walls fill input window. No space left for transmission (d*nslit > w1).\n", NAME_CURRENT_COMP));

  if (!w2) w2=w1;
  if (!h2) h2=h1;

  if (mcgravitation) G=-GRAVITY;
  mcLocG = rot_apply(ROT_A_CURRENT_COMP, coords_set(0,G,0));
  coords_get(mcLocG, &Gx, &Gy, &Gz);

  strcpy(GVars.compcurname, NAME_CURRENT_COMP);

  if (l > 0 && nelements > 0) {

    Gravity_guide_Init(&GVars,
      w1, h1, w2, h2, l, R0,
      Qc, alpha, m, W, nslit, d,
      Gx, Gy, Gz, mleft, mright, mtop,
      mbottom, nhslit, wavy_lr, wavy_tb, wavy_z, wavy,
      chamfers_z, chamfers_lr, chamfers_tb, chamfers,nu,phase,aleft,aright,atop,abottom);
    if (!G) for (i=0; i<5; GVars.A[i++] = 0);
    if (GVars.fc_freq != 0 || GVars.fc_phase != 0) {
      if (w1 != w2 || h1 != h2)
      exit(fprintf(stderr,"Guide_gravity: %s: rotating slit pack must be straight (w1=w2 and h1=h2).\n", NAME_CURRENT_COMP));
      printf("Guide_gravity: %s: Fermi Chopper mode: frequency=%g [Hz] phase=%g [deg]\n",
        NAME_CURRENT_COMP, GVars.fc_freq, GVars.fc_phase);
    }
  } else printf("Guide_gravity: %s: unactivated (l=0 or nelements=0)\n", NAME_CURRENT_COMP);

  #undef w1
  #undef h1
  #undef w2
  #undef h2
  #undef l
  #undef R0
  #undef Qc
  #undef alpha
  #undef m
  #undef W
  #undef nslit
  #undef d
  #undef mleft
  #undef mright
  #undef mtop
  #undef mbottom
  #undef nhslit
  #undef G
  #undef aleft
  #undef aright
  #undef atop
  #undef abottom
  #undef wavy
  #undef wavy_z
  #undef wavy_tb
  #undef wavy_lr
  #undef chamfers
  #undef chamfers_z
  #undef chamfers_lr
  #undef chamfers_tb
  #undef nelements
  #undef nu
  #undef phase
  #undef reflect
  #undef GVars
  #undef pTable
  return(_comp);
} /* class_Guide_gravity_init */

_class_Monitor_nD *class_Monitor_nD_init(_class_Monitor_nD *_comp
) {
  #define user1 (_comp->_parameters.user1)
  #define user2 (_comp->_parameters.user2)
  #define user3 (_comp->_parameters.user3)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define zdepth (_comp->_parameters.zdepth)
  #define xmin (_comp->_parameters.xmin)
  #define xmax (_comp->_parameters.xmax)
  #define ymin (_comp->_parameters.ymin)
  #define ymax (_comp->_parameters.ymax)
  #define zmin (_comp->_parameters.zmin)
  #define zmax (_comp->_parameters.zmax)
  #define bins (_comp->_parameters.bins)
  #define min (_comp->_parameters.min)
  #define max (_comp->_parameters.max)
  #define restore_neutron (_comp->_parameters.restore_neutron)
  #define radius (_comp->_parameters.radius)
  #define options (_comp->_parameters.options)
  #define filename (_comp->_parameters.filename)
  #define geometry (_comp->_parameters.geometry)
  #define username1 (_comp->_parameters.username1)
  #define username2 (_comp->_parameters.username2)
  #define username3 (_comp->_parameters.username3)
  #define DEFS (_comp->_parameters.DEFS)
  #define Vars (_comp->_parameters.Vars)
  #define detector (_comp->_parameters.detector)
  #define offdata (_comp->_parameters.offdata)
  char tmp[CHAR_BUF_LENGTH];
  strcpy(Vars.compcurname, NAME_CURRENT_COMP);
  if (options != NULL)
    strncpy(Vars.option, options, CHAR_BUF_LENGTH);
  else {
    strcpy(Vars.option, "x y");
    printf("Monitor_nD: %s has no option specified. Setting to PSD ('x y') monitor.\n", NAME_CURRENT_COMP);
  }
  Vars.compcurpos = POS_A_CURRENT_COMP;

  if (strstr(Vars.option, "source"))
    strcat(Vars.option, " list, x y z vx vy vz t sx sy sz ");

  if (bins) { sprintf(tmp, " all bins=%ld ", (long)bins); strcat(Vars.option, tmp); }
  if (min > -FLT_MAX && max < FLT_MAX) { sprintf(tmp, " all limits=[%g %g]", min, max); strcat(Vars.option, tmp); }
  else if (min > -FLT_MAX) { sprintf(tmp, " all min=%g", min); strcat(Vars.option, tmp); }
  else if (max <  FLT_MAX) { sprintf(tmp, " all max=%g", max); strcat(Vars.option, tmp); }

  strncpy(Vars.UserName1,
    username1 && strlen(username1) && strcmp(username1, "0") && strcmp(username1, "NULL") ?
    username1 : "", 128);
  strncpy(Vars.UserName2,
    username2 && strlen(username2) && strcmp(username2, "0") && strcmp(username2, "NULL") ?
    username2 : "", 128);
  strncpy(Vars.UserName3,
    username3 && strlen(username3) && strcmp(username3, "0") && strcmp(username3, "NULL") ?
    username3 : "", 128);
  if (radius) {
    xwidth = zdepth = 2*radius;
    if (yheight && !strstr(Vars.option, "cylinder") && !strstr(Vars.option, "banana") && !strstr(Vars.option, "sphere"))
      strcat(Vars.option, " banana");
    else if (!yheight && !strstr(Vars.option ,"sphere")) {
      strcat(Vars.option, " sphere");
      yheight=2*radius;
    }
  }
  int offflag=0;
  if (geometry && strlen(geometry) && strcmp(geometry,"0") && strcmp(geometry, "NULL"))
    if (!off_init(  geometry, xwidth, yheight, zdepth, 1, &offdata )) {
      printf("Monitor_nD: %s could not initiate the OFF geometry %s. \n"
             "            Defaulting to normal Monitor dimensions.\n",
             NAME_CURRENT_COMP, geometry);
      strcpy(geometry, "");
    } else {
      offflag=1;
    }

  if (!radius && !xwidth && !yheight && !zdepth && !xmin && !xmax && !ymin && !ymax &&
    !strstr(Vars.option, "previous") && (!geometry || !strlen(geometry)))
    exit(printf("Monitor_nD: %s has no dimension specified. Aborting (radius, xwidth, yheight, zdepth, previous, geometry).\n", NAME_CURRENT_COMP));

  Monitor_nD_Init(&DEFS, &Vars, xwidth, yheight, zdepth, xmin,xmax,ymin,ymax,zmin,zmax,offflag);

  if (Vars.Flag_OFF) {
    offdata.mantidflag=Vars.Flag_mantid;
    offdata.mantidoffset=Vars.Coord_Min[Vars.Coord_Number-1];
  }


  if (filename && strlen(filename) && strcmp(filename,"NULL") && strcmp(filename,"0"))
    strncpy(Vars.Mon_File, filename, 128);

  /* check if user given filename with ext will be used more than once */
  if ( ((Vars.Flag_Multiple && Vars.Coord_Number > 1) || Vars.Flag_List) && strchr(Vars.Mon_File,'.') )
  { char *XY; XY = strrchr(Vars.Mon_File,'.'); *XY='_'; }

  if (restore_neutron) Vars.Flag_parallel=1;
  detector.m = 0;

#ifdef USE_MPI
MPI_MASTER(
  if (strstr(Vars.option, "auto") && mpi_node_count > 1)
    printf("Monitor_nD: %s is using automatic limits option 'auto' together with MPI.\n"
           "WARNING     this may create incorrect distributions (but integrated flux will be right).\n", NAME_CURRENT_COMP);
);
#endif
  #undef user1
  #undef user2
  #undef user3
  #undef xwidth
  #undef yheight
  #undef zdepth
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef zmin
  #undef zmax
  #undef bins
  #undef min
  #undef max
  #undef restore_neutron
  #undef radius
  #undef options
  #undef filename
  #undef geometry
  #undef username1
  #undef username2
  #undef username3
  #undef DEFS
  #undef Vars
  #undef detector
  #undef offdata
  return(_comp);
} /* class_Monitor_nD_init */

_class_DiskChopper *class_DiskChopper_init(_class_DiskChopper *_comp
) {
  #define theta_0 (_comp->_parameters.theta_0)
  #define radius (_comp->_parameters.radius)
  #define yheight (_comp->_parameters.yheight)
  #define nu (_comp->_parameters.nu)
  #define nslit (_comp->_parameters.nslit)
  #define jitter (_comp->_parameters.jitter)
  #define delay (_comp->_parameters.delay)
  #define isfirst (_comp->_parameters.isfirst)
  #define n_pulse (_comp->_parameters.n_pulse)
  #define abs_out (_comp->_parameters.abs_out)
  #define phase (_comp->_parameters.phase)
  #define xwidth (_comp->_parameters.xwidth)
  #define verbose (_comp->_parameters.verbose)
  #define Tg (_comp->_parameters.Tg)
  #define To (_comp->_parameters.To)
  #define delta_y (_comp->_parameters.delta_y)
  #define height (_comp->_parameters.height)
  #define omega (_comp->_parameters.omega)
/* If slit height 'unset', assume full opening */
if (yheight == 0) {
        height=radius;
      } else {
        height=yheight;
      }
      delta_y = radius-height/2; /* radius at beam center */
      omega=2.0*PI*nu; /* rad/s */
      if (xwidth && !theta_0 && radius) theta_0 = 2*RAD2DEG*asin(xwidth/2/delta_y);

      if (nslit<=0 || theta_0 <= 0 || radius <=0)
      { fprintf(stderr,"DiskChopper: %s: nslit, theta_0 and radius must be > 0\n", NAME_CURRENT_COMP);
        exit(-1); }
      if (nslit*theta_0 >= 360)
      { fprintf(stderr,"DiskChopper: %s: nslit * theta_0 exceeds 2PI\n", NAME_CURRENT_COMP);
        exit(-1); }
      if (yheight && yheight>radius) {
        fprintf(stderr,"DiskChopper: %s: yheight must be < radius\n", NAME_CURRENT_COMP);
        exit(-1); }
      if (isfirst && n_pulse <=0)
      { fprintf(stderr,"DiskChopper: %s: wrong First chopper pulse number (n_pulse=%g)\n", NAME_CURRENT_COMP, n_pulse);
        exit(-1); }
      if (!omega) {
        fprintf(stderr,"DiskChopper: %s WARNING: chopper frequency is 0!\n", NAME_CURRENT_COMP);
        omega = 1e-15; /* We should actually use machine epsilon here... */
      }
      if (!abs_out) {
        fprintf(stderr,"DiskChopper: %s WARNING: chopper will NOT absorb neutrons outside radius %g [m]\n", NAME_CURRENT_COMP, radius);
      }

      theta_0*=DEG2RAD;


      /* Calulate delay from phase and vice versa */
      if (phase) {
        if (delay) {
          fprintf(stderr,"DiskChopper: %s WARNING: delay AND phase specified. Using phase setting\n", NAME_CURRENT_COMP);
        }
        phase*=DEG2RAD;
        /* 'Delay' should always be a delay, taking rotation direction into account: */
        delay=phase/fabs(omega);
      } else {
        phase=delay*omega;  /* rad */
      }

      /* Time from opening of slit to next opening of slit */
      Tg=2.0*PI/fabs(omega)/nslit;

      /* How long can neutrons pass the Chopper at a single point */
      To=theta_0/fabs(omega);

      if (!xwidth) xwidth=2*delta_y*sin(theta_0/2);

      if (verbose && nu) {
        printf("DiskChopper: %s: frequency=%g [Hz] %g [rpm], time frame=%g [s] phase=%g [deg]\n",
          NAME_CURRENT_COMP, nu, nu*60, Tg, phase*RAD2DEG);
        printf("             %g slits, angle=%g [deg] height=%g [m], width=%g [m] at radius=%g [m]\n",
          nslit, theta_0*RAD2DEG, height, xwidth, delta_y);
      }
  #undef theta_0
  #undef radius
  #undef yheight
  #undef nu
  #undef nslit
  #undef jitter
  #undef delay
  #undef isfirst
  #undef n_pulse
  #undef abs_out
  #undef phase
  #undef xwidth
  #undef verbose
  #undef Tg
  #undef To
  #undef delta_y
  #undef height
  #undef omega
  return(_comp);
} /* class_DiskChopper_init */

_class_Isotropic_Sqw *class_Isotropic_Sqw_init(_class_Isotropic_Sqw *_comp
) {
  #define powder_format (_comp->_parameters.powder_format)
  #define Sqw_coh (_comp->_parameters.Sqw_coh)
  #define Sqw_inc (_comp->_parameters.Sqw_inc)
  #define geometry (_comp->_parameters.geometry)
  #define radius (_comp->_parameters.radius)
  #define thickness (_comp->_parameters.thickness)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define zdepth (_comp->_parameters.zdepth)
  #define threshold (_comp->_parameters.threshold)
  #define order (_comp->_parameters.order)
  #define T (_comp->_parameters.T)
  #define verbose (_comp->_parameters.verbose)
  #define d_phi (_comp->_parameters.d_phi)
  #define concentric (_comp->_parameters.concentric)
  #define rho (_comp->_parameters.rho)
  #define sigma_abs (_comp->_parameters.sigma_abs)
  #define sigma_coh (_comp->_parameters.sigma_coh)
  #define sigma_inc (_comp->_parameters.sigma_inc)
  #define classical (_comp->_parameters.classical)
  #define powder_Dd (_comp->_parameters.powder_Dd)
  #define powder_DW (_comp->_parameters.powder_DW)
  #define powder_Vc (_comp->_parameters.powder_Vc)
  #define density (_comp->_parameters.density)
  #define weight (_comp->_parameters.weight)
  #define p_interact (_comp->_parameters.p_interact)
  #define norm (_comp->_parameters.norm)
  #define powder_barns (_comp->_parameters.powder_barns)
  #define quantum_correction (_comp->_parameters.quantum_correction)
  #define VarSqw (_comp->_parameters.VarSqw)
  #define columns (_comp->_parameters.columns)
  #define offdata (_comp->_parameters.offdata)
  int i;
  /* check for parameters */
  columns = (int*)powder_format;

  VarSqw.verbose_output= verbose;
  VarSqw.shape = -1; /* -1:no shape, 0:cyl, 1:box, 2:sphere, 3:any-shape  */
  if (geometry && strlen(geometry) && strcmp(geometry, "NULL") && strcmp(geometry, "0")) {
	  if (off_init(geometry, xwidth, yheight, zdepth, 0, &offdata)) {
      VarSqw.shape=3; thickness=0; concentric=0;
    }
  }
  else if (xwidth && yheight && zdepth)  VarSqw.shape=1; /* box */
  else if (radius > 0 && yheight)        VarSqw.shape=0; /* cylinder */
  else if (radius > 0 && !yheight)       VarSqw.shape=2; /* sphere */

  if (VarSqw.shape < 0)
    exit(fprintf(stderr,"Isotropic_Sqw: %s: sample has invalid dimensions.\n"
                        "ERROR          Please check parameter values (xwidth, yheight, zdepth, radius).\n", NAME_CURRENT_COMP));



  if (thickness) {
    if (radius && (radius < fabs(thickness) )) {
      MPI_MASTER(
      fprintf(stderr,"Isotropic_Sqw: %s: hollow sample thickness is larger than its volume (sphere/cylinder).\n"
                     "WARNING        Please check parameter values. Using bulk sample (thickness=0).\n", NAME_CURRENT_COMP);
      );
      thickness=0;
    }
    else if (!radius && (xwidth < 2*fabs(thickness) || yheight < 2*fabs(thickness) || zdepth < 2*fabs(thickness))) {
      MPI_MASTER(
      fprintf(stderr,"Isotropic_Sqw: %s: hollow sample thickness is larger than its volume (box).\n"
                     "WARNING        Please check parameter values.\n", NAME_CURRENT_COMP);
      );
    }
  }
  MPI_MASTER(
  if (VarSqw.verbose_output) {
    switch (VarSqw.shape) {
      case 0: printf("Isotropic_Sqw: %s: is a %scylinder: radius=%f thickness=%f height=%f [J Comp Phys 228 (2009) 5251]\n",
              NAME_CURRENT_COMP, (thickness ? "hollow " : ""),
              radius,fabs(thickness),yheight);
              break;
      case 1: printf("Isotropic_Sqw: %s: is a %sbox: width=%f height=%f depth=%f \n",
              NAME_CURRENT_COMP, (thickness ? "hollow " : ""), xwidth,yheight,zdepth);
              break;
      case 2: printf("Isotropic_Sqw: %s: is a %ssphere: radius=%f thickness=%f\n",
              NAME_CURRENT_COMP, (thickness ? "hollow " : ""),
              radius,fabs(thickness));
              break;
      case 3: printf("Isotropic_Sqw: %s: is a volume defined from file %s\n",
              NAME_CURRENT_COMP, geometry);
    }
  }
  );

  if (concentric && !thickness) {
    MPI_MASTER(
    printf("Isotropic_Sqw: %s:Can not use concentric mode\n"
           "WARNING        on non hollow shape. Ignoring.\n",
           NAME_CURRENT_COMP);
    );
    concentric=0;
  }

  strncpy(VarSqw.compname, NAME_CURRENT_COMP, 256);
  VarSqw.T2E       =(1/11.605);   /* Kelvin to meV = 1000*KB/e */
  VarSqw.sqSE2K    = (V2K*SE2V)*(V2K*SE2V);
  VarSqw.sqw_threshold = (threshold > 0 ? threshold : 0);
  VarSqw.s_abs     = sigma_abs;
  VarSqw.s_coh     = sigma_coh;
  VarSqw.s_inc     = sigma_inc; /* s_scatt member initialized in Sqw_init */
  VarSqw.maxloop   = 100;       /* atempts to close triangle */
  VarSqw.minevents = 100;       /* minimal # of events required to get dynamical range */
  VarSqw.neutron_removed = 0;
  VarSqw.neutron_enter   = 0;
  VarSqw.neutron_pmult   = 0;
  VarSqw.neutron_exit    = 0;
  VarSqw.mat_rho       = rho;
  VarSqw.sqw_norm  = norm;
  VarSqw.mean_scatt= 0;
  VarSqw.mean_abs  = 0;
  VarSqw.psum_scatt= 0;
  VarSqw.single_coh= 0;
  VarSqw.single_inc= 0;
  VarSqw.multi     = 0;
  VarSqw.barns     = powder_barns;
  VarSqw.sqw_classical = classical;
  VarSqw.lookup_length=100;
  VarSqw.mat_weight    = weight;
  VarSqw.mat_density   = density;
  if (quantum_correction && strlen(quantum_correction))
    strncpy(VarSqw.Q_correction, quantum_correction, 256);
  else
    strncpy(VarSqw.Q_correction, "default", 256);

  /* PowderN compatibility members */
  VarSqw.Dd        = powder_Dd;
  VarSqw.DWfactor  = powder_DW;
  VarSqw.Temperature= T;
  for (i=0; i< 9; i++) VarSqw.column_order[i] = columns[i];
  VarSqw.column_order[8] = (VarSqw.column_order[0] >= 0 ? 0 : 2);

  /* optional ways to define rho */
  if (!VarSqw.mat_rho && powder_Vc > 0)
    VarSqw.mat_rho = 1/powder_Vc;
  /* import the data files ================================================== */
  if (!Sqw_init(&VarSqw, Sqw_coh, Sqw_inc)) {
    MPI_MASTER(
    printf("Isotropic_Sqw: %s: ERROR importing data files (Sqw_init coh=%s inc=%s).\n",NAME_CURRENT_COMP, Sqw_coh, Sqw_inc);
    );
  }
  if ( VarSqw.s_coh < 0) VarSqw.s_coh=0;
  if ( VarSqw.s_inc < 0) VarSqw.s_inc=0;
  if ( VarSqw.s_abs < 0) VarSqw.s_abs=0;
  if ((VarSqw.s_coh > 0 || VarSqw.s_inc > 0) && VarSqw.mat_rho <= 0) {
    MPI_MASTER(
    printf("Isotropic_Sqw: %s: WARNING: Null density (V_rho). Unactivating component.\n",NAME_CURRENT_COMP);
    );
    VarSqw.s_coh=VarSqw.s_inc=0;
  }
  /* 100: convert from barns to fm^2 */
  VarSqw.my_a_v  =(VarSqw.mat_rho*100*VarSqw.s_abs*2200);
  VarSqw.my_s    =(VarSqw.mat_rho*100*(VarSqw.s_coh>0 ? VarSqw.s_coh : 0
                                     +VarSqw.s_inc>0 ? VarSqw.s_inc : 0));
  MPI_MASTER(
  if ((VarSqw.s_coh > 0 || VarSqw.s_inc > 0) && !VarSqw.Temperature
   && (VarSqw.Data_coh.intensity || VarSqw.Data_inc.intensity)
   && VarSqw.verbose_output)
    printf("Isotropic_Sqw: %s: Sample temperature not defined (T=0).\n"
           "Warning        Disabling detailed balance.\n", NAME_CURRENT_COMP);
  if (VarSqw.s_coh<=0 && VarSqw.s_inc<=0) {
    printf("Isotropic_Sqw: %s: Scattering cross section is zero\n"
           "ERROR          (sigma_coh, sigma_inc).\n",NAME_CURRENT_COMP);
  }
  );
  if (d_phi) d_phi = fabs(d_phi)*DEG2RAD;

  if (d_phi > PI) d_phi = 0; /* V_scatt on 4*PI */

  if (d_phi && order != 1) {
    MPI_MASTER(
    printf("Isotropic_Sqw: %s: Focusing can only apply for single\n"
           "               scattering. Setting to order=1.\n",
           NAME_CURRENT_COMP);
    );
    order = 1;
  }

  /* request statistics */
  if (VarSqw.verbose_output > 1) {
    Sqw_diagnosis(&VarSqw, &VarSqw.Data_coh);
    Sqw_diagnosis(&VarSqw, &VarSqw.Data_inc);
  }

  for (i=0; i < 2; i++) {
    struct Sqw_Data_struct Data_sqw;
    Data_sqw =  (i == 0 ? VarSqw.Data_coh : VarSqw.Data_inc);
    Table_Free(&(Data_sqw.Sqw));
  }

/* end INITIALIZE */
  #undef powder_format
  #undef Sqw_coh
  #undef Sqw_inc
  #undef geometry
  #undef radius
  #undef thickness
  #undef xwidth
  #undef yheight
  #undef zdepth
  #undef threshold
  #undef order
  #undef T
  #undef verbose
  #undef d_phi
  #undef concentric
  #undef rho
  #undef sigma_abs
  #undef sigma_coh
  #undef sigma_inc
  #undef classical
  #undef powder_Dd
  #undef powder_DW
  #undef powder_Vc
  #undef density
  #undef weight
  #undef p_interact
  #undef norm
  #undef powder_barns
  #undef quantum_correction
  #undef VarSqw
  #undef columns
  #undef offdata
  return(_comp);
} /* class_Isotropic_Sqw_init */

_class_PSD_Detector *class_PSD_Detector_init(_class_PSD_Detector *_comp
) {
  #define nx (_comp->_parameters.nx)
  #define ny (_comp->_parameters.ny)
  #define xwidth (_comp->_parameters.xwidth)
  #define radius (_comp->_parameters.radius)
  #define awidth (_comp->_parameters.awidth)
  #define yheight (_comp->_parameters.yheight)
  #define zdepth (_comp->_parameters.zdepth)
  #define threshold (_comp->_parameters.threshold)
  #define PressureConv (_comp->_parameters.PressureConv)
  #define PressureStop (_comp->_parameters.PressureStop)
  #define interpolate (_comp->_parameters.interpolate)
  #define p_interact (_comp->_parameters.p_interact)
  #define verbose (_comp->_parameters.verbose)
  #define LensOn (_comp->_parameters.LensOn)
  #define dc (_comp->_parameters.dc)
  #define borderx (_comp->_parameters.borderx)
  #define bordery (_comp->_parameters.bordery)
  #define xChDivRelSigma (_comp->_parameters.xChDivRelSigma)
  #define yChDivRelSigma (_comp->_parameters.yChDivRelSigma)
  #define bufsize (_comp->_parameters.bufsize)
  #define restore_neutron (_comp->_parameters.restore_neutron)
  #define angle (_comp->_parameters.angle)
  #define type (_comp->_parameters.type)
  #define filename (_comp->_parameters.filename)
  #define FN_Conv (_comp->_parameters.FN_Conv)
  #define FN_Stop (_comp->_parameters.FN_Stop)
  #define PSD_N (_comp->_parameters.PSD_N)
  #define PSD_p (_comp->_parameters.PSD_p)
  #define PSD_p2 (_comp->_parameters.PSD_p2)
  #define EAP (_comp->_parameters.EAP)
  #define M1P1 (_comp->_parameters.M1P1)
  #define PosAP (_comp->_parameters.PosAP)
  #define EAT (_comp->_parameters.EAT)
  #define M1T1 (_comp->_parameters.M1T1)
  #define PosAT (_comp->_parameters.PosAT)
  #define CrossSectionHe (_comp->_parameters.CrossSectionHe)
  #define CountNeutrons (_comp->_parameters.CountNeutrons)
  #define GeomCumul (_comp->_parameters.GeomCumul)
  #define AbsCumul (_comp->_parameters.AbsCumul)
  #define SensVolCumul (_comp->_parameters.SensVolCumul)
  #define DetCumul (_comp->_parameters.DetCumul)
  #define PHSpectrum (_comp->_parameters.PHSpectrum)
  #define PHSpectrum0 (_comp->_parameters.PHSpectrum0)
  #define PHSpectrum2 (_comp->_parameters.PHSpectrum2)
  #define PHSpectrum_n (_comp->_parameters.PHSpectrum_n)
  #define nH_p (_comp->_parameters.nH_p)
  #define nH_t (_comp->_parameters.nH_t)
  #define FullEnergyP (_comp->_parameters.FullEnergyP)
  #define FullEnergyT (_comp->_parameters.FullEnergyT)
  #define VariousErrors (_comp->_parameters.VariousErrors)
  #define DetectorType (_comp->_parameters.DetectorType)
  #define DEFS (_comp->_parameters.DEFS)
  #define Vars (_comp->_parameters.Vars)
  CountNeutrons=0;
  GeomCumul=0;
  AbsCumul=0;
  SensVolCumul=0;
  DetCumul=0;
  VariousErrors=0;
  DetectorType=0;

  long  i,j;

  t_Table Part_He_3_p;
  t_Table Part_Stop_p;
  long nS_p; // number of table elements for the proton in the stopping gas
  t_Table Part_He_3_t;
  t_Table Part_Stop_t;
  long nS_t; // number of table elements for the triton in the stopping gas
  t_Table Part_He_3_n,Part_Stop_n; // third entry in the file, containing cross section only

  double *E_p; // energies from the table for protons, for which stopping powers are given
  double *dEdx_p; // resulting stopping powers of the counting gas mixture


  double *DAP; // Delta Array of energies for the proton from the list
  double *MuAP; // Distances in mu traversed for each delta energy in the list

  double *TempVar34;
  double *ETP1; // Energy transferred by proton
  double *PTP1; // Positions of transfer of energy by proton

  double *E_t; // energies from the table for tritons, for which stopping powers are given
  double *dEdx_t; // resulting stopping powers of the counting gas mixture


  double *DAT; // Delta Array of energies for the triton from the list
  double *MuAT; // Distances in mu traversed for each delta energy in the list

  double *TempVar35;
  double *ETT1; // Energy transferred by triton
  double *PTT1; // Positions of transfer of energy by triton

  /* GEOMETRY stuff ********************************************************* */
  if (borderx==-2) borderx = xwidth/nx;
  else if (borderx==-1) borderx = zdepth;
  else if (borderx<0) {
    fprintf(stderr,"PSD_Detector: %s: Negative x border zone specified. Exit.\n",
      NAME_CURRENT_COMP);
    exit(0);
  }
  if (bordery==-2) bordery = yheight/ny;
  else if (bordery==-1) bordery = zdepth;
  else if (bordery<0) {
    fprintf(stderr,"PSD_Detector: %s: Negative y border zone specified. Exit.\n",
      NAME_CURRENT_COMP);
    exit(0);
  }

  if (xwidth>0) { /* panel */
    DetectorType+=1;
    if (zdepth<0) {
    fprintf(stderr,"PSD_Detector: %s: Detector box has no zdepth. Exit.\n",
                           NAME_CURRENT_COMP);
    exit(0);
    }
  }
  if (radius && angle>0)       awidth=angle/360*2*PI*radius;
  else if (radius && awidth>0) angle = awidth*360/2/PI/radius;
  if (radius>0 && !awidth) {
    DetectorType+=2; /* cylinder */
  }
  if (awidth>0 && radius) {
    DetectorType+=4; /* banana */
    if (zdepth<=0) {
    fprintf(stderr,"PSD_Detector: %s: Detector 'banana' has no zdepth. Exit.\n",
                           NAME_CURRENT_COMP);
    exit(0);
    }
    if (radius<=0) {
    fprintf(stderr,"PSD_Detector: %s: Non-positive radial distance to the sample was specified. Exit.\n",
                           NAME_CURRENT_COMP);
    exit(0);
    }
    if (!dc && zdepth && LensOn) dc = zdepth/2;
    if (((dc>zdepth) || (dc<0)) && (LensOn==1)) {
    fprintf(stderr,"PSD_Detector: %s: Electrostatic lens was turned on, but\n"
                   "              no valid zdepth of the cathode plane was specified. Exit.\n",
                                  NAME_CURRENT_COMP);
    exit(0);
    }
    if (awidth+2*borderx>2*PI*radius) {
    fprintf(stderr,"PSD_Detector: %s: Detector comes perilously close to encompassing 2pi. Exit.\n",
                           NAME_CURRENT_COMP);
    exit(0);
    }
  }
  if (yheight<=0) {
    fprintf(stderr,"PSD_Detector: %s: Detector has no height (yheight). Exit.\n",
                           NAME_CURRENT_COMP);
    exit(0);
  }
  if ( (DetectorType!=1) && (DetectorType!=2) &&(DetectorType!=4) ) {
    fprintf(stderr,"PSD_Detector: %s: Detector has conflicting size\n"
                     "              specifications, i.e. combinations of xwidth, radius\n"
                     "              or awidth, or none at all. Exit.\n",
                                    NAME_CURRENT_COMP);
    exit(0);
  }
  if (verbose) {
    printf("PSD_Detector: %s: Geometry is ", NAME_CURRENT_COMP);
    switch (DetectorType) {
    case 1: printf("box\n"); break;
    case 2: printf("tube\n"); break;
    case 4: printf("cylinder ('banana') opening angle=%g [deg] arc length=%g [m]\n", angle, awidth); break;
    default:
      printf("not defined\n");
    }
  }

/* Gas tables  ************************************************************** */

  PSD_N = create_darr2d(nx, ny);
  PSD_p = create_darr2d(nx, ny);
  PSD_p2 = create_darr2d(nx, ny);

  for (i=0; i<nx; i++)
    for (j=0; j<ny; j++)
    {
      PSD_N[i][j] = 0;
      PSD_p[i][j] = 0;
      PSD_p2[i][j] = 0;
    }

  Table_Read(&Part_He_3_p,FN_Conv,1);
  nH_p = Part_He_3_p.rows;
  Table_Read(&Part_Stop_p,FN_Stop,1);
  nS_p = Part_Stop_p.rows;
  Table_Read(&Part_He_3_t,FN_Conv,2);
  nH_t  = Part_He_3_t.rows;
  Table_Read(&Part_Stop_t,FN_Stop,2);
  nS_t = Part_Stop_t.rows;
  /* if Gas tables can not be found, try by pre-pending Gas_tables */
  char tmp[256];
  if (!nH_p) {
  	sprintf(tmp, "Gas_tables%s%s", MC_PATHSEP_S, FN_Conv);
  	Table_Read(&Part_He_3_p,tmp,1);
    nH_p = Part_He_3_p.rows;
  }
  if (!nS_p) {
  	sprintf(tmp, "Gas_tables%s%s", MC_PATHSEP_S, FN_Stop);
  	Table_Read(&Part_Stop_p,tmp,1);
    nS_p = Part_Stop_p.rows;
  }
  if (!nH_t) {
    sprintf(tmp, "Gas_tables%s%s", MC_PATHSEP_S, FN_Conv);
  	Table_Read(&Part_He_3_t,tmp,2);
    nH_t  = Part_He_3_t.rows;
  }
  if (!nS_t) {
    sprintf(tmp, "Gas_tables%s%s", MC_PATHSEP_S, FN_Stop);
  	Table_Read(&Part_Stop_t,tmp,2);
    nS_t = Part_Stop_t.rows;
  }


  if (nH_p != nS_p || nH_t != nS_t) {
    fprintf(stderr,"PSD_Detector: %s: Data files for helium %s and stopping gas %s\n"
                   "              have different number of entries. Exit.\n",
                   NAME_CURRENT_COMP, FN_Conv, FN_Stop);
    exit(0);
  }

  E_p = (double *) malloc(nH_p*sizeof(double));

  dEdx_p = (double *) malloc(nH_p*sizeof(double));
  for (i=0; i<nH_p; i++) {
    E_p[i] = Table_Index(Part_He_3_p,i,0);
    if (E_p[i] != Table_Index(Part_Stop_p,i,0)) {
        fprintf(stderr,"PSD_Detector: %s: Data files for helium %s and stopping gas %s\n"
                       "              list different energies (proton part). Exit.\n",
                       NAME_CURRENT_COMP, FN_Conv, FN_Stop);
        exit(0); }
    dEdx_p[i] = PressureConv*Table_Index(Part_He_3_p,i,1) + PressureStop*Table_Index(Part_Stop_p,i,1);
  }
  FullEnergyP = E_p[nH_p-1];
  EAP = (double *) malloc((nH_p+1)*sizeof(double));
  EAP[0] = 0;
  for (i=1; i<nH_p+1; i++) { EAP[i]       = E_p[i-1]; }
  free(E_p);
  DAP = (double *) malloc((nH_p)*sizeof(double));
  for (i=0; i<nH_p; i++)   { DAP[i]       = EAP[i+1] - EAP[i]; }

  MuAP = (double *) malloc((nH_p+1)*sizeof(double));
  MuAP[0] = 0;
  for (i=1; i<nH_p+1; i++) { MuAP[i]      = DAP[i-1] / dEdx_p[i-1]; }
  free(DAP);

  PosAP = (double *) malloc((nH_p+1)*sizeof(double));
  PosAP[0] = 0;
  for (i=1; i<nH_p+1; i++) { PosAP[i]     = PosAP[i-1] + 1e-6*MuAP[i]; }
  free(MuAP);
  // going to flip the arrays PosAP and EAP

  TempVar34 = (double *) malloc((nH_p+1)*sizeof(double));
  for (i=0; i<nH_p+1; i++) { TempVar34[i] = PosAP[nH_p] - PosAP[i]; }
  for (i=0; i<nH_p+1; i++) { PosAP[i]     = TempVar34[nH_p-i]; }
  for (i=0; i<nH_p+1; i++) { TempVar34[i] = EAP[nH_p] - EAP[i]; }
  for (i=0; i<nH_p+1; i++) { EAP[i]       = TempVar34[nH_p-i]; }
  free(TempVar34);
  // done flipping PosAP and EAP

  ETP1 = (double *) malloc((nH_p)*sizeof(double));
  for (i=0; i<nH_p; i++) { ETP1[i]        = EAP[i+1] - EAP[i]; }

  PTP1 = (double *) malloc((nH_p)*sizeof(double));
  for (i=0; i<nH_p; i++) { PTP1[i] = 0.5 * ( PosAP[i] + PosAP[i+1] ); }

  M1P1 = (double *) malloc((nH_p+1)*sizeof(double));
  M1P1[0] = 0;
  // not yet divided by EAP, see below
  for (i=1; i<nH_p+1; i++) { M1P1[i] = M1P1[i-1] + PTP1[i-1]*ETP1[i-1];  }
  for (i=1; i<nH_p+1; i++) { M1P1[i] = M1P1[i] / EAP[i]; }
  free(ETP1);
  free(PTP1);
  /* keep EAP, M1P1 and PosAP, free others */

  E_t = (double *) malloc(nH_t*sizeof(double));

  dEdx_t = (double *) malloc(nH_t*sizeof(double));
  for (i=0; i<nH_t; i++) {
    E_t[i] = Table_Index(Part_He_3_t,i,0);
    if (E_t[i] != Table_Index(Part_Stop_t,i,0)) {
        fprintf(stderr,"PSD_Detector: %s: Data files for helium %s and stopping gas %s\n"
                       "              list different energies (triton part). Exit.\n",
                       NAME_CURRENT_COMP, FN_Conv, FN_Stop);
        exit(0); }
    dEdx_t[i] = PressureConv*Table_Index(Part_He_3_t,i,1) + PressureStop*Table_Index(Part_Stop_t,i,1);
  }
  FullEnergyT = E_t[nH_t-1];

  EAT = (double *) malloc((nH_t+1)*sizeof(double));
  EAT[0] = 0;
  for (i=1; i<nH_t+1; i++) { EAT[i] = E_t[i-1]; }
  free(E_t);

  DAT = (double *) malloc((nH_t)*sizeof(double));
  for (i=0; i<nH_t; i++)   { DAT[i] = EAT[i+1] - EAT[i]; }

  MuAT = (double *) malloc((nH_t+1)*sizeof(double));
  MuAT[0] = 0;
  for (i=1; i<nH_t+1; i++) { MuAT[i] = DAT[i-1] / dEdx_t[i-1]; }
  free(DAT);

  PosAT = (double *) malloc((nH_t+1)*sizeof(double));
  PosAT[0] = 0;
  for (i=1; i<nH_t+1; i++) { PosAT[i] = PosAT[i-1] + 1e-6*MuAT[i]; }
  free(MuAT);
  // going to flip the arrays PosAT and EAT

  TempVar35 = (double *) malloc((nH_t+1)*sizeof(double));
  for (i=0; i<nH_t+1; i++) { TempVar35[i] = PosAT[nH_t] - PosAT[i]; }
  for (i=0; i<nH_t+1; i++) { PosAT[i]     = TempVar35[nH_t-i]; }
  for (i=0; i<nH_t+1; i++) { TempVar35[i] = EAT[nH_t] - EAT[i]; }
  for (i=0; i<nH_t+1; i++) { EAT[i]       = TempVar35[nH_t-i]; }
  free(TempVar35);
  // done flipping PosAT and EAT

  ETT1 = (double *) malloc((nH_t)*sizeof(double));
  for (i=0; i<nH_t; i++) { ETT1[i] = EAT[i+1] - EAT[i]; }

  PTT1 = (double *) malloc((nH_t)*sizeof(double));
  for (i=0; i<nH_t; i++) { PTT1[i] = 0.5 * ( PosAT[i] + PosAT[i+1] ); }

  M1T1 = (double *) malloc((nH_t+1)*sizeof(double));
  M1T1[0] = 0;
  // not yet divided by EAT, see below
  for (i=1; i<nH_t+1; i++) { M1T1[i] = M1T1[i-1] + PTT1[i-1]*ETT1[i-1]; }
  for (i=1; i<nH_t+1; i++) { M1T1[i] = M1T1[i] / EAT[i]; }
  free(ETT1);
  free(PTT1);

  if (verbose) {
    fprintf(stdout,"# Position (m) Energy (keV) COG position (m) PROTON\n");
    for (i=0; i<nH_p+1; i++) {
      fprintf(stdout,"  %.2e      %5.1f       %.2e\n",PosAP[i],EAP[i],M1P1[i]);
    }
    fprintf(stdout,"# Position (m) Energy (keV) COG position (m) TRITON\n");
    for (i=0; i<nH_t+1; i++) {
      fprintf(stdout,"  %.2e      %5.1f       %.2e\n",PosAT[i],EAT[i],M1T1[i]);
    }
  }

  /* keep EAT, M1T1 and PosAT, free others */

  PHSpectrum_n = (long)((FullEnergyP+FullEnergyT)*1.2);
  PHSpectrum0  = (double *) malloc(PHSpectrum_n * sizeof(double)); //
  PHSpectrum   = (double *) malloc(PHSpectrum_n * sizeof(double)); //
  PHSpectrum2  = (double *) malloc(PHSpectrum_n * sizeof(double)); //
  for (i=0; i<PHSpectrum_n; i++) {
    PHSpectrum0[i] = 0;
    PHSpectrum[i]  = 0;
    PHSpectrum2[i] = 0;
  }

  if (threshold<0) threshold=0;

  Table_Read(&Part_He_3_n,FN_Conv,3); // the absorption cross section for 0.18 nm neutrons in the converter is read from the file
  if (!Part_He_3_n.rows) {
    sprintf(tmp, "Gas_tables%s%s", MC_PATHSEP_S, FN_Conv);
  	Table_Read(&Part_He_3_n,tmp,3);
  }
  if ( (Part_He_3_n.rows!=1) || (Part_He_3_n.columns!=1) ) {
      fprintf(stderr,"PSD_Detector: %s: Problem: the third part of the converter\n"
                     "              stopping power file should contain only \n"
                     "              one value: the absorption cross section for\n"
                     "              0.18 nm neutrons in the converter.\n",
                     NAME_CURRENT_COMP);
      exit(-1);
  }
  Table_Read(&Part_Stop_n,FN_Stop,3); // the absorption cross section for 0.18 nm neutrons in the quench is read from the file
  if (!Part_Stop_n.rows) {
    sprintf(tmp, "Gas_tables%s%s", MC_PATHSEP_S, FN_Stop);
  	Table_Read(&Part_Stop_n,tmp,3);
  }
  if ( (Part_Stop_n.rows!=1) || (Part_Stop_n.columns!=1)) {
      fprintf(stderr,"PSD_Detector: %s: Problem: the third part of the stopping \n"
                     "              gas stopping power file should contain only \n"
                     "              one value: the absorption cross section for\n"
                     "              0.18 nm neutrons in the converter.\n",
                     NAME_CURRENT_COMP);
      exit(-1);
  }
  CrossSectionHe=Table_Index(Part_He_3_n,0,0);
  if (CrossSectionHe!=Table_Index(Part_Stop_n,0,0) ) {
      fprintf(stderr,"PSD_Detector: %s: Problem: the absorption cross section\n"
                     "              read from the stopping gas file does not\n"
                     "              match that read from the converter gas file.\n",
                     NAME_CURRENT_COMP);
      exit(-1);
  }

  Table_Free(&Part_He_3_p);
  Table_Free(&Part_Stop_p);
  Table_Free(&Part_He_3_t);
  Table_Free(&Part_Stop_t);
  Table_Free(&Part_He_3_n);
  Table_Free(&Part_Stop_n);

  if (p_interact>1) p_interact=1;

/* handle event files as in Virtual_output, with bufsize=0 ****************** */
  long element_size;
  if (type && strlen(type) && strcmp(type, "NULL") && strcmp(type, "0")) {
    printf("PSD_Detector: %s: saving detector signal as events\n", NAME_CURRENT_COMP);

    strcpy(Vars.compcurname, NAME_CURRENT_COMP);
    if (bufsize > 0) sprintf(Vars.option, "list=%15g", bufsize);
    else strcpy(Vars.option, "list all");
    strcat(Vars.option,", borders, x y z vx vy vz t sx sy sz");

    Monitor_nD_Init(&DEFS, &Vars, 0.1, 0.1, 0, 0,0,0,0,0,0,0); /* dims for mcdisplay */
    if (filename && strlen(filename) && strcmp(filename, "0") && strcmp(filename, "NULL"))
      strncpy(Vars.Mon_File, filename, 128);

    if (bufsize > 0)
    printf("Warning: PSD_Detector: %s: buffer size=%g not recommended\n", NAME_CURRENT_COMP, bufsize);
    if (bufsize > 0) printf(
           "PSD_Detector: %s: Beware virtual output generated file size (max %g Mo)\n"
           "WARNING         Memory required is %g Mo\n", NAME_CURRENT_COMP,
           bufsize*element_size/1e6, bufsize*sizeof(double)/1e6);
  }
  #undef nx
  #undef ny
  #undef xwidth
  #undef radius
  #undef awidth
  #undef yheight
  #undef zdepth
  #undef threshold
  #undef PressureConv
  #undef PressureStop
  #undef interpolate
  #undef p_interact
  #undef verbose
  #undef LensOn
  #undef dc
  #undef borderx
  #undef bordery
  #undef xChDivRelSigma
  #undef yChDivRelSigma
  #undef bufsize
  #undef restore_neutron
  #undef angle
  #undef type
  #undef filename
  #undef FN_Conv
  #undef FN_Stop
  #undef PSD_N
  #undef PSD_p
  #undef PSD_p2
  #undef EAP
  #undef M1P1
  #undef PosAP
  #undef EAT
  #undef M1T1
  #undef PosAT
  #undef CrossSectionHe
  #undef CountNeutrons
  #undef GeomCumul
  #undef AbsCumul
  #undef SensVolCumul
  #undef DetCumul
  #undef PHSpectrum
  #undef PHSpectrum0
  #undef PHSpectrum2
  #undef PHSpectrum_n
  #undef nH_p
  #undef nH_t
  #undef FullEnergyP
  #undef FullEnergyT
  #undef VariousErrors
  #undef DetectorType
  #undef DEFS
  #undef Vars
  return(_comp);
} /* class_PSD_Detector_init */



int init(void) { /* called by mccode_main for HZB_NEAT:INITIALISE */
  DEBUG_INSTR();

  /* code_main/parseoptions/readparams sets instrument parameters value */
  stracpy(instrument->_name, "HZB_NEAT", 256);

  /* Instrument 'HZB_NEAT' INITIALISE */
  SIG_MESSAGE("[HZB_NEAT] INITIALISE [HZB_NEAT.instr:60]");
  #define lambda (instrument->_parameters._lambda)
  #define dlambda (instrument->_parameters._dlambda)
  #define rpm (instrument->_parameters._rpm)
  #define coh (instrument->_parameters._coh)
  #define inc (instrument->_parameters._inc)
{
  double Chopper2Detector = 11.97+2.143+0.165;
  double KI, Vi, EI;
  double Emin, Emax, Vmin, Vmax, Tmin, Tmax;

  KI= 2*PI/lambda;
  Vi = K2V*fabs(KI);
  EI = VS2E*Vi*Vi;

  printf("%s: Detailed NEAT/TOF configuration\n", NAME_INSTRUMENT);
  printf("* Incoming beam: lambda=%.4g [Angs] EI=%.4g [meV]  KI=%.4g [Angs-1] Vi=%g [m/s]\n",
    lambda, EI, KI, Vi);
  time_elastic = 252.78*Chopper2Detector*lambda*1e-6; /* time from Chopper1 to reach Detector */
  printf("  Sample: coh=%s inc=%s\n", coh, inc);
  printf("          Elastic line time: %g [us]\n", time_elastic*1e6);

  /* compute a +/- 20 meV */
  Emin = EI-20; if (Emin < EI/3) Emin=EI/3;
  Emax = EI+20;
  Vmin = sqrt(Emin/VS2E);
  Vmax = sqrt(Emax/VS2E);
  Tmin = 2.5/Vmax;
  Tmax = 2.5/Vmin;

  sprintf(detector_options,
    "banana, abs angle limits=[13 136] bins=300, time limits=[%g %g] bins=256",
    time_elastic+Tmin, time_elastic+Tmax);
}
  #undef lambda
  #undef dlambda
  #undef rpm
  #undef coh
  #undef inc
  _PG_setpos(); /* type Progress_bar */
  _Source_setpos(); /* type Source_gen */
  _Guide0_setpos(); /* type Guide_gravity */
  _Guide0_4_setpos(); /* type Guide_gravity */
  _Guide0_5_setpos(); /* type Guide_gravity */
  _Guide0_6_setpos(); /* type Guide_gravity */
  _Guide0_7_setpos(); /* type Guide_gravity */
  _Guide0_8_setpos(); /* type Guide_gravity */
  _Guide0_9_setpos(); /* type Guide_gravity */
  _Guide0_10_setpos(); /* type Guide_gravity */
  _Guide0_11_setpos(); /* type Guide_gravity */
  _Guide0_12_setpos(); /* type Guide_gravity */
  _Guide0_13_setpos(); /* type Guide_gravity */
  _Guide0_14_setpos(); /* type Guide_gravity */
  _Guide0_15_setpos(); /* type Guide_gravity */
  _Guide0_16_setpos(); /* type Guide_gravity */
  _Guide0_17_setpos(); /* type Guide_gravity */
  _Guide0_18_setpos(); /* type Guide_gravity */
  _Guide0_19_setpos(); /* type Guide_gravity */
  _Guide_PSD_setpos(); /* type Monitor_nD */
  _Guide_Lambda_setpos(); /* type Monitor_nD */
  _Chopper1_setpos(); /* type DiskChopper */
  _Guide_SGS1_setpos(); /* type Guide_gravity */
  _Guide_SGS2_setpos(); /* type Guide_gravity */
  _Guide_SGS3_setpos(); /* type Guide_gravity */
  _Guide_CGS_setpos(); /* type Guide_gravity */
  _Guide2_PSD_setpos(); /* type Monitor_nD */
  _Guide2_Lambda_setpos(); /* type Monitor_nD */
  _Guide2_dXY_setpos(); /* type Monitor_nD */
  _Chopper6_setpos(); /* type DiskChopper */
  _Guide_time_setpos(); /* type Monitor_nD */
  _Guide_DGS_setpos(); /* type Guide_gravity */
  _Sample_pos_setpos(); /* type Arm */
  _Sample_setpos(); /* type Isotropic_Sqw */
  _Ideal_Det_setpos(); /* type Monitor_nD */
  _Detector_setpos(); /* type PSD_Detector */

  /* call iteratively all components INITIALISE */
  class_Progress_bar_init(_PG);

  class_Source_gen_init(_Source);

  class_Guide_gravity_init(_Guide0);

  class_Guide_gravity_init(_Guide0_4);

  class_Guide_gravity_init(_Guide0_5);

  class_Guide_gravity_init(_Guide0_6);

  class_Guide_gravity_init(_Guide0_7);

  class_Guide_gravity_init(_Guide0_8);

  class_Guide_gravity_init(_Guide0_9);

  class_Guide_gravity_init(_Guide0_10);

  class_Guide_gravity_init(_Guide0_11);

  class_Guide_gravity_init(_Guide0_12);

  class_Guide_gravity_init(_Guide0_13);

  class_Guide_gravity_init(_Guide0_14);

  class_Guide_gravity_init(_Guide0_15);

  class_Guide_gravity_init(_Guide0_16);

  class_Guide_gravity_init(_Guide0_17);

  class_Guide_gravity_init(_Guide0_18);

  class_Guide_gravity_init(_Guide0_19);

  class_Monitor_nD_init(_Guide_PSD);

  class_Monitor_nD_init(_Guide_Lambda);

  class_DiskChopper_init(_Chopper1);

  class_Guide_gravity_init(_Guide_SGS1);

  class_Guide_gravity_init(_Guide_SGS2);

  class_Guide_gravity_init(_Guide_SGS3);

  class_Guide_gravity_init(_Guide_CGS);

  class_Monitor_nD_init(_Guide2_PSD);

  class_Monitor_nD_init(_Guide2_Lambda);

  class_Monitor_nD_init(_Guide2_dXY);

  class_DiskChopper_init(_Chopper6);

  class_Monitor_nD_init(_Guide_time);

  class_Guide_gravity_init(_Guide_DGS);


  class_Isotropic_Sqw_init(_Sample);

  class_Monitor_nD_init(_Ideal_Det);

  class_PSD_Detector_init(_Detector);

  if (mcdotrace) display();
  DEBUG_INSTR_END();

  return(0);
} /* init */

/*******************************************************************************
* components TRACE
*******************************************************************************/

#define x (_particle->x)
#define y (_particle->y)
#define z (_particle->z)
#define vx (_particle->vx)
#define vy (_particle->vy)
#define vz (_particle->vz)
#define t (_particle->t)
#define sx (_particle->sx)
#define sy (_particle->sy)
#define sz (_particle->sz)
#define p (_particle->p)
// user variables:

#define SCATTERED (_particle->_scattered)
#define RESTORE (_particle->_restore)
#define RESTORE_NEUTRON(_index, ...) _particle->_restore = _index;
#define ABSORBED (_particle->_absorbed)
#define ABSORB0 do { DEBUG_STATE(); DEBUG_ABSORB(); MAGNET_OFF; ABSORBED++; return(_comp); } while(0)
#define ABSORB ABSORB0
#pragma acc routine seq
_class_Progress_bar *class_Progress_bar_trace(_class_Progress_bar *_comp
  , _class_particle *_particle) {
  ABSORBED=SCATTERED=RESTORE=0;

  #define profile (_comp->_parameters.profile)
  #define percent (_comp->_parameters.percent)
  #define flag_save (_comp->_parameters.flag_save)
  #define minutes (_comp->_parameters.minutes)
  #define IntermediateCnts (_comp->_parameters.IntermediateCnts)
  #define StartTime (_comp->_parameters.StartTime)
  #define EndTime (_comp->_parameters.EndTime)
  #define CurrentTime (_comp->_parameters.CurrentTime)
  double ncount;
  ncount = mcget_run_num();
  if (!StartTime) {
    time(&StartTime); /* compute starting time */
    IntermediateCnts = 1e3;
  }
  time_t NowTime;
  time(&NowTime);
  /* compute initial estimate of computation duration */
  if (!EndTime && ncount >= IntermediateCnts) {
    CurrentTime = NowTime;
    if (difftime(NowTime,StartTime) > 10 && ncount) { /* wait 10 sec before writing ETA */
      EndTime = StartTime + (time_t)(difftime(NowTime,StartTime)
				     *(double)mcget_ncount()/ncount);
      IntermediateCnts = 0;
      fprintf(stdout, "\nTrace ETA ");
      if (difftime(EndTime,StartTime) < 60.0)
        fprintf(stdout, "%g [s] %% ", difftime(EndTime,StartTime));
      else if (difftime(EndTime,StartTime) > 3600.0)
        fprintf(stdout, "%g [h] %% ", difftime(EndTime,StartTime)/3600.0);
      else
        fprintf(stdout, "%g [min] %% ", difftime(EndTime,StartTime)/60.0);
    } else IntermediateCnts += 1e3;
    fflush(stdout);
  }

  /* display percentage when percent or minutes have reached step */
  if (EndTime && mcget_ncount() &&
    (    (minutes && difftime(NowTime,CurrentTime) > minutes*60)
      || (percent && !minutes && ncount >= IntermediateCnts))   )
  {
    fprintf(stdout, "%d ", (int)(ncount*100.0/mcget_ncount())); fflush(stdout);
    CurrentTime = NowTime;

    IntermediateCnts = ncount + percent*mcget_ncount()/100;
    /* check that next intermediate ncount check is a multiple of the desired percentage */
    IntermediateCnts = floor(IntermediateCnts*100/percent/mcget_ncount())*percent*mcget_ncount()/100;
    /* raise flag to indicate that we did something */
    SCATTER;
    if (flag_save) save(NULL);
  }
  #undef profile
  #undef percent
  #undef flag_save
  #undef minutes
  #undef IntermediateCnts
  #undef StartTime
  #undef EndTime
  #undef CurrentTime
  return(_comp);
} /* class_Progress_bar_trace */

#pragma acc routine seq
_class_Source_gen *class_Source_gen_trace(_class_Source_gen *_comp
  , _class_particle *_particle) {
  ABSORBED=SCATTERED=RESTORE=0;

  #define flux_file (_comp->_parameters.flux_file)
  #define xdiv_file (_comp->_parameters.xdiv_file)
  #define ydiv_file (_comp->_parameters.ydiv_file)
  #define radius (_comp->_parameters.radius)
  #define dist (_comp->_parameters.dist)
  #define focus_xw (_comp->_parameters.focus_xw)
  #define focus_yh (_comp->_parameters.focus_yh)
  #define focus_aw (_comp->_parameters.focus_aw)
  #define focus_ah (_comp->_parameters.focus_ah)
  #define E0 (_comp->_parameters.E0)
  #define dE (_comp->_parameters.dE)
  #define lambda0 (_comp->_parameters.lambda0)
  #define dlambda (_comp->_parameters.dlambda)
  #define I1 (_comp->_parameters.I1)
  #define yheight (_comp->_parameters.yheight)
  #define xwidth (_comp->_parameters.xwidth)
  #define verbose (_comp->_parameters.verbose)
  #define T1 (_comp->_parameters.T1)
  #define flux_file_perAA (_comp->_parameters.flux_file_perAA)
  #define flux_file_log (_comp->_parameters.flux_file_log)
  #define Lmin (_comp->_parameters.Lmin)
  #define Lmax (_comp->_parameters.Lmax)
  #define Emin (_comp->_parameters.Emin)
  #define Emax (_comp->_parameters.Emax)
  #define T2 (_comp->_parameters.T2)
  #define I2 (_comp->_parameters.I2)
  #define T3 (_comp->_parameters.T3)
  #define I3 (_comp->_parameters.I3)
  #define zdepth (_comp->_parameters.zdepth)
  #define target_index (_comp->_parameters.target_index)
  #define p_in (_comp->_parameters.p_in)
  #define lambda1 (_comp->_parameters.lambda1)
  #define lambda2 (_comp->_parameters.lambda2)
  #define lambda3 (_comp->_parameters.lambda3)
  #define pTable (_comp->_parameters.pTable)
  #define pTable_x (_comp->_parameters.pTable_x)
  #define pTable_y (_comp->_parameters.pTable_y)
  #define pTable_xmin (_comp->_parameters.pTable_xmin)
  #define pTable_xmax (_comp->_parameters.pTable_xmax)
  #define pTable_xsum (_comp->_parameters.pTable_xsum)
  #define pTable_ymin (_comp->_parameters.pTable_ymin)
  #define pTable_ymax (_comp->_parameters.pTable_ymax)
  #define pTable_ysum (_comp->_parameters.pTable_ysum)
  #define pTable_dxmin (_comp->_parameters.pTable_dxmin)
  #define pTable_dxmax (_comp->_parameters.pTable_dxmax)
  #define pTable_dymin (_comp->_parameters.pTable_dymin)
  #define pTable_dymax (_comp->_parameters.pTable_dymax)
  double dx=0,dy=0,xf,yf,rf,pdir,chi,v,r, lambda;
  double Maxwell;

  if (verbose >= 0)
  {

    z=0;

    if (radius)
    {
      chi=2*PI*rand01();                          /* Choose point on source */
      r=sqrt(rand01())*radius;                    /* with uniform distribution. */
      x=r*cos(chi);
      y=r*sin(chi);
    }
    else
    {
      x = xwidth*randpm1()/2;   /* select point on source (uniform) */
      y = yheight*randpm1()/2;
    }
    if (zdepth != 0)
      z = zdepth*randpm1()/2;
  /* Assume linear wavelength distribution */
    lambda = lambda0+dlambda*randpm1();
    if (lambda <= 0) ABSORB;
    v = K2V*(2*PI/lambda);

    if (!focus_ah && !focus_aw) {
      randvec_target_rect_real(&xf, &yf, &rf, &pdir,
       0, 0, dist, focus_xw, focus_yh, ROT_A_CURRENT_COMP, x, y, z, 2);

      dx = xf-x;
      dy = yf-y;
      rf = sqrt(dx*dx+dy*dy+dist*dist);

      vz=v*dist/rf;
      vy=v*dy/rf;
      vx=v*dx/rf;
    } else {

      randvec_target_rect_angular(&vx, &vy, &vz, &pdir,
          0, 0, 1, focus_aw*DEG2RAD, focus_ah*DEG2RAD, ROT_A_CURRENT_COMP);
      dx = vx; dy = vy; /* from unit vector */
      vx *= v; vy *= v; vz *= v;
    }
    p = p_in*pdir;

    /* spectral dependency from files or Maxwellians */
    if (flux_file && strlen(flux_file) && strcmp(flux_file,"NULL") && strcmp(flux_file,"0"))
    {
       double binwidth=Table_Value(pTable, lambda, 1);
       if (flux_file_log) binwidth=exp(binwidth);
       p *= binwidth;
    }
    else if (T1 > 0 && I1 > 0)
    {
      Maxwell = I1 * SG_Maxwell(lambda, T1);;  /* 1/AA */

      if ((T2 > 0) && (I2 > 0))
      {
        Maxwell += I2 * SG_Maxwell(lambda, T2);
      }
      if ((T3 > 0) && (I3 > 0))
      {
        Maxwell += I3 * SG_Maxwell(lambda, T3);;
      }
      p *= Maxwell;
    }

    /* optional x-xdiv and y-ydiv weightening: position=along columns, div=along rows */
    if (xdiv_file && strlen(xdiv_file)
      && strcmp(xdiv_file,"NULL") && strcmp(xdiv_file,"0") && pTable_xsum > 0) {
      double i,j;
      j = (x-            pTable_xmin) /(pTable_xmax -pTable_xmin) *pTable_x.columns;
      i = (atan2(dx,rf)*RAD2DEG-pTable_dxmin)/(pTable_dxmax-pTable_dxmin)*pTable_x.rows;
      r = Table_Value2d(pTable_x, i,j); /* row, column */
      p *= r/pTable_xsum;
    }
    if (ydiv_file && strlen(ydiv_file)
       && strcmp(ydiv_file,"NULL") && strcmp(ydiv_file,"0") && pTable_ysum > 0) {
      double i,j;
      j = (y-            pTable_ymin) /(pTable_ymax -pTable_ymin) *pTable_y.columns;
      i = (atan2(dy,rf)*RAD2DEG-  pTable_dymin)/(pTable_dymax-pTable_dymin)*pTable_y.rows;
      r = Table_Value2d(pTable_y, i,j);
      p *= r/pTable_ysum;
    }
    SCATTER;
  }
  #undef flux_file
  #undef xdiv_file
  #undef ydiv_file
  #undef radius
  #undef dist
  #undef focus_xw
  #undef focus_yh
  #undef focus_aw
  #undef focus_ah
  #undef E0
  #undef dE
  #undef lambda0
  #undef dlambda
  #undef I1
  #undef yheight
  #undef xwidth
  #undef verbose
  #undef T1
  #undef flux_file_perAA
  #undef flux_file_log
  #undef Lmin
  #undef Lmax
  #undef Emin
  #undef Emax
  #undef T2
  #undef I2
  #undef T3
  #undef I3
  #undef zdepth
  #undef target_index
  #undef p_in
  #undef lambda1
  #undef lambda2
  #undef lambda3
  #undef pTable
  #undef pTable_x
  #undef pTable_y
  #undef pTable_xmin
  #undef pTable_xmax
  #undef pTable_xsum
  #undef pTable_ymin
  #undef pTable_ymax
  #undef pTable_ysum
  #undef pTable_dxmin
  #undef pTable_dxmax
  #undef pTable_dymin
  #undef pTable_dymax
  return(_comp);
} /* class_Source_gen_trace */

#pragma acc routine seq
_class_Guide_gravity *class_Guide_gravity_trace(_class_Guide_gravity *_comp
  , _class_particle *_particle) {
  ABSORBED=SCATTERED=RESTORE=0;

  #define w1 (_comp->_parameters.w1)
  #define h1 (_comp->_parameters.h1)
  #define w2 (_comp->_parameters.w2)
  #define h2 (_comp->_parameters.h2)
  #define l (_comp->_parameters.l)
  #define R0 (_comp->_parameters.R0)
  #define Qc (_comp->_parameters.Qc)
  #define alpha (_comp->_parameters.alpha)
  #define m (_comp->_parameters.m)
  #define W (_comp->_parameters.W)
  #define nslit (_comp->_parameters.nslit)
  #define d (_comp->_parameters.d)
  #define mleft (_comp->_parameters.mleft)
  #define mright (_comp->_parameters.mright)
  #define mtop (_comp->_parameters.mtop)
  #define mbottom (_comp->_parameters.mbottom)
  #define nhslit (_comp->_parameters.nhslit)
  #define G (_comp->_parameters.G)
  #define aleft (_comp->_parameters.aleft)
  #define aright (_comp->_parameters.aright)
  #define atop (_comp->_parameters.atop)
  #define abottom (_comp->_parameters.abottom)
  #define wavy (_comp->_parameters.wavy)
  #define wavy_z (_comp->_parameters.wavy_z)
  #define wavy_tb (_comp->_parameters.wavy_tb)
  #define wavy_lr (_comp->_parameters.wavy_lr)
  #define chamfers (_comp->_parameters.chamfers)
  #define chamfers_z (_comp->_parameters.chamfers_z)
  #define chamfers_lr (_comp->_parameters.chamfers_lr)
  #define chamfers_tb (_comp->_parameters.chamfers_tb)
  #define nelements (_comp->_parameters.nelements)
  #define nu (_comp->_parameters.nu)
  #define phase (_comp->_parameters.phase)
  #define reflect (_comp->_parameters.reflect)
  #define GVars (_comp->_parameters.GVars)
  #define pTable (_comp->_parameters.pTable)
  if (l > 0 && nelements > 0) {
    double B, C, dt;
    int    ret, bounces = 0, i=0;
    double this_width, this_height;
    double angle=0;

    if (GVars.fc_freq != 0 || GVars.fc_phase != 0) { /* rotate neutron w/r to guide element */
      /* approximation of rotating straight Fermi Chopper */
      Coords   X = coords_set(x,y,z-l/2);  /* current coordinates of neutron in centered static frame */
      Rotation R;
      double dt=(-z+l/2)/vz; /* time shift to each center of slit package */
      angle=fmod(360*GVars.fc_freq*(t+dt)+GVars.fc_phase, 360); /* in deg */
      /* modify angle so that Z0 guide side is always in front of incoming neutron */
      if (angle > 90 && angle < 270) { angle -= 180; }
      angle *= DEG2RAD;
      rot_set_rotation(R, 0, -angle, 0); /* will rotate neutron instead of comp: negative side */
      /* apply rotation to centered coordinates */
      Coords   RX = rot_apply(R, X);
      coords_get(RX, &x, &y, &z);
      z = z+l/2;
      /* rotate speed */
      X  = coords_set(vx,vy,vz);
      RX = rot_apply(R, X);
      coords_get(RX, &vx, &vy, &vz);
    }

    for (i=0; i<7; GVars.N_reflection[i++] = 0);

    /* propagate to box input (with gravitation) in comp local coords */
    /* A = 0.5 n.g; B = n.v; C = n.(r-W); */
    /* 0=Z0 side: n=(0, 0, -l) ; W = (0, 0, 0) (at z=0, guide input)*/
    B = -l*vz; C = -l*z;

    ret = solve_2nd_order(&dt, NULL, GVars.A[0], B, C);
    if (ret==0) ABSORB;

    if (dt>0.0) PROP_GRAV_DT(dt, GVars.gx, GVars.gy, GVars.gz); else if (angle) ABSORB;
    GVars.N_reflection[6]++;

    this_width  = w1;
    this_height = h1;

  /* check if we are in the box input, else absorb */
    if (fabs(x) > this_width/2 || fabs(y) > this_height/2)
      ABSORB;
    else
    {
      double w_edge, w_adj; /* Channel displacement on X */
      double h_edge, h_adj; /* Channel displacement on Y */
      double w_chnum,h_chnum; /* channel indexes */

      SCATTER;

      /* X: Shift origin to center of channel hit (absorb if hit dividing walls) */
      x += w1/2.0;
      w_chnum = floor(x/(GVars.w1c+d));  /* 0= right side, nslit+1=left side  */
      w_edge  = w_chnum*(GVars.w1c+d);
      if(x - w_edge > GVars.w1c)
      {
        x -= w1/2.0; /* Re-adjust origin */
        ABSORB;
      }
      w_adj = w_edge + (GVars.w1c)/2.0;
      x -= w_adj; w_adj -=  w1/2.0;

      /* Y: Shift origin to center of channel hit (absorb if hit dividing walls) */
      y += h1/2.0;
      h_chnum = floor(y/(GVars.h1c+d));  /* 0= lower side, nslit+1=upper side  */
      h_edge  = h_chnum*(GVars.h1c+d);
      if(y - h_edge > GVars.h1c)
      {
        y -= h1/2.0; /* Re-adjust origin */
        ABSORB;
      }
      h_adj = h_edge + (GVars.h1c)/2.0;
      y -= h_adj; h_adj -=  h1/2.0;

      /* neutron is now in the input window of the guide */
      /* do loops on reflections in the box */
      for(;;)
      {
        /* get intersections for all box sides */
        double q, nx,ny,nz;
        double this_length;
        int side=0;

        bounces++;
        /* now look for intersection with guide sides and exit */
        side = Gravity_guide_Trace(&dt, &GVars, x, y, z,
            vx, vy, vz, w_chnum, nslit, h_chnum, nhslit,
            &nx, &ny, &nz);

        /* only positive dt are valid */
        /* exit reflection loops if no intersection (neutron is after box) */
        if (side == 0 || dt <= 0)
          { if (GVars.warnings < 100)
              fprintf(stderr,"%s: warning: neutron has entered guide, but can not exit !\n", GVars.compcurname);
            GVars.warnings++;
            x += w_adj; y += h_adj; ABSORB; } /* should never occur */

        /* propagate to dt */
        PROP_GRAV_DT(dt, GVars.gx, GVars.gy, GVars.gz);

        /* do reflection on speed for l/r/u/d sides */
        if (side == 5) /* neutron reaches end of guide: end loop and exit comp */
          { GVars.N_reflection[side]++; x += w_adj; y += h_adj; SCATTER; x -= w_adj; y -= h_adj; break; }
        /* else reflection on a guide wall */
        if(GVars.M[side] == 0 || Qc == 0 || R0 == 0)  /* walls are absorbing */
          { x += w_adj; y += h_adj; ABSORB; }
        /* handle chamfers */
        this_width = w1+(w2-w1)*z/l;
        this_height= h1+(h2-h1)*z/l;
        this_length= fmod(z, l/nelements);
        /* absorb on input/output of element parts */
        if (GVars.chamfer_z && (this_length<GVars.chamfer_z || this_length>l/nelements-GVars.chamfer_z))
        { x += w_adj; y += h_adj; ABSORB; }
        /* absorb on l/r/t/b sides */
        if (GVars.chamfer_lr && (side==1 || side==2) && (fabs(y+h_adj)>this_height/2-GVars.chamfer_lr))
        { x += w_adj; y += h_adj; ABSORB; }
        if (GVars.chamfer_tb && (side==3 || side==4) && (fabs(x+w_adj)>this_width/2- GVars.chamfer_tb))
        { x += w_adj; y += h_adj; ABSORB; }
        /* change/mirror velocity: h_f = v - n.2*n.v/|n|^2 */
        GVars.N_reflection[side]++; /* GVars.norm_n2 > 0 was checked at INIT */
        /* compute n.v using current values */
        B = scalar_prod(vx,vy,vz,nx,ny,nz);
        dt = 2*B/GVars.norm_n2[side]; /* 2*n.v/|n|^2 */
        vx -= nx*dt;
        vy -= ny*dt;
        vz -= nz*dt;

        /* compute q and modify neutron weight */
        /* scattering q=|n_i-n_f| = V2Q*|vf - v| = V2Q*2*n.v/|n| */
        q = 2*V2Q*fabs(B)/GVars.norm_n[side];

        if (reflect && strlen(reflect) && strcmp(reflect,"NULL") && strcmp(reflect,"0"))
          TableReflecFunc(q, &pTable, &B);
        else {
          double par[] = {R0, Qc, GVars.Alpha[side], GVars.M[side], W};
          StdReflecFunc(q, par, &B);
        }
        if (B <= 0) { x += w_adj; y += h_adj; ABSORB; }
        else p *= B;
        x += w_adj; y += h_adj; SCATTER; x -= w_adj; y -= h_adj;
        GVars.N_reflection[0]++;
        /* go to the next reflection */
        if (bounces > 1000) ABSORB;
      } /* end for */
      x += w_adj; y += h_adj; /* Re-adjust origin after SCATTER */
    }

    if (GVars.fc_freq != 0 || GVars.fc_phase != 0) { /* rotate back neutron w/r to guide element */
      /* approximation of rotating straight Fermi Chopper */
      Coords   X = coords_set(x,y,z-l/2);  /* current coordinates of neutron in centered static frame */
      Rotation R;
      rot_set_rotation(R, 0, angle, 0); /* will rotate back neutron: positive side */
      /* apply rotation to centered coordinates */
      Coords   RX = rot_apply(R, X);
      coords_get(RX, &x, &y, &z);
      z = z+l/2;
      /* rotate speed */
      X  = coords_set(vx,vy,vz);
      RX = rot_apply(R, X);
      coords_get(RX, &vx, &vy, &vz);
    }

  } /* if l */
  #undef w1
  #undef h1
  #undef w2
  #undef h2
  #undef l
  #undef R0
  #undef Qc
  #undef alpha
  #undef m
  #undef W
  #undef nslit
  #undef d
  #undef mleft
  #undef mright
  #undef mtop
  #undef mbottom
  #undef nhslit
  #undef G
  #undef aleft
  #undef aright
  #undef atop
  #undef abottom
  #undef wavy
  #undef wavy_z
  #undef wavy_tb
  #undef wavy_lr
  #undef chamfers
  #undef chamfers_z
  #undef chamfers_lr
  #undef chamfers_tb
  #undef nelements
  #undef nu
  #undef phase
  #undef reflect
  #undef GVars
  #undef pTable
  return(_comp);
} /* class_Guide_gravity_trace */

#pragma acc routine seq
_class_Monitor_nD *class_Monitor_nD_trace(_class_Monitor_nD *_comp
  , _class_particle *_particle) {
  ABSORBED=SCATTERED=RESTORE=0;

  #define user1 (_comp->_parameters.user1)
  #define user2 (_comp->_parameters.user2)
  #define user3 (_comp->_parameters.user3)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define zdepth (_comp->_parameters.zdepth)
  #define xmin (_comp->_parameters.xmin)
  #define xmax (_comp->_parameters.xmax)
  #define ymin (_comp->_parameters.ymin)
  #define ymax (_comp->_parameters.ymax)
  #define zmin (_comp->_parameters.zmin)
  #define zmax (_comp->_parameters.zmax)
  #define bins (_comp->_parameters.bins)
  #define min (_comp->_parameters.min)
  #define max (_comp->_parameters.max)
  #define restore_neutron (_comp->_parameters.restore_neutron)
  #define radius (_comp->_parameters.radius)
  #define options (_comp->_parameters.options)
  #define filename (_comp->_parameters.filename)
  #define geometry (_comp->_parameters.geometry)
  #define username1 (_comp->_parameters.username1)
  #define username2 (_comp->_parameters.username2)
  #define username3 (_comp->_parameters.username3)
  #define DEFS (_comp->_parameters.DEFS)
  #define Vars (_comp->_parameters.Vars)
  #define detector (_comp->_parameters.detector)
  #define offdata (_comp->_parameters.offdata)
  if(!strcmp("Guide_PSD", _comp->_name)) {
    user1 = FLT_MAX;
    user2 = FLT_MAX;
    user3 = FLT_MAX;
  }
  if(!strcmp("Guide_Lambda", _comp->_name)) {
    user1 = FLT_MAX;
    user2 = FLT_MAX;
    user3 = FLT_MAX;
  }
  if(!strcmp("Guide2_PSD", _comp->_name)) {
    user1 = FLT_MAX;
    user2 = FLT_MAX;
    user3 = FLT_MAX;
  }
  if(!strcmp("Guide2_Lambda", _comp->_name)) {
    user1 = FLT_MAX;
    user2 = FLT_MAX;
    user3 = FLT_MAX;
  }
  if(!strcmp("Guide2_dXY", _comp->_name)) {
    user1 = FLT_MAX;
    user2 = FLT_MAX;
    user3 = FLT_MAX;
  }
  if(!strcmp("Guide_time", _comp->_name)) {
    user1 = FLT_MAX;
    user2 = FLT_MAX;
    user3 = FLT_MAX;
  }
  if(!strcmp("Ideal_Det", _comp->_name)) {
    user1 = FLT_MAX;
    user2 = FLT_MAX;
    user3 = FLT_MAX;
  }

  double  XY=0;
  double  t0 = 0;
  double  t1 = 0;
  double  pp;
  int     intersect   = 0;
  char    Flag_Restore = 0;

  Vars.UserVariable1 = user1;
  Vars.UserVariable2 = user2;
  Vars.UserVariable3 = user3;

  /* this is done automatically
    STORE_NEUTRON(INDEX_CURRENT_COMP, x, y, z, vx, vy, vz, t, sx, sy, sz, p);
  */

  if (geometry && strlen(geometry) && strcmp(geometry,"0") && strcmp(geometry, "NULL"))
  {
    /* determine intersections with object */
    intersect = off_intersect_all(&t0, &t1, NULL, NULL,
       x,y,z, vx, vy, vz, &offdata );
    if (Vars.Flag_mantid) {
      if(intersect) {
        Vars.OFF_polyidx=(offdata.intersects[offdata.nextintersect]).index;
      } else {
        Vars.OFF_polyidx=-1;
      }
    }
  }
  else if ( (abs(Vars.Flag_Shape) == DEFS.SHAPE_SQUARE)
            || (abs(Vars.Flag_Shape) == DEFS.SHAPE_DISK) ) /* square xy or disk xy */
  {
    // propagate to xy plane and find intersection
    // make sure the event is recoverable afterwards
    t0 = t;
    ALLOW_BACKPROP;
    PROP_Z0;
    if ( (t>=t0) && (z==0.0) ) // forward propagation to xy plane was successful
    {
      if (abs(Vars.Flag_Shape) == DEFS.SHAPE_SQUARE)
      {
        // square xy
        intersect = (x>=Vars.mxmin && x<=Vars.mxmax && y>=Vars.mymin && y<=Vars.mymax);
      }
      else
      {
        // disk xy
        intersect = (SQR(x) + SQR(y)) <= SQR(Vars.Sphere_Radius);
      }
    }
    else
    {
      intersect=0;
    }
  }
  else if (abs(Vars.Flag_Shape) == DEFS.SHAPE_SPHERE) /* sphere */
  {
    intersect = sphere_intersect(&t0, &t1, x, y, z, vx, vy, vz, Vars.Sphere_Radius);
  /*      intersect = (intersect && t0 > 0); */
  }
  else if ((abs(Vars.Flag_Shape) == DEFS.SHAPE_CYLIND) || (abs(Vars.Flag_Shape) == DEFS.SHAPE_BANANA)) /* cylinder */
  {
    intersect = cylinder_intersect(&t0, &t1, x, y, z, vx, vy, vz, Vars.Sphere_Radius, Vars.Cylinder_Height);
  }
  else if (abs(Vars.Flag_Shape) == DEFS.SHAPE_BOX) /* box */
  {
    intersect = box_intersect(&t0, &t1, x, y, z, vx, vy, vz,
                              fabs(Vars.mxmax-Vars.mxmin), fabs(Vars.mymax-Vars.mymin), fabs(Vars.mzmax-Vars.mzmin));
  }
  else if (abs(Vars.Flag_Shape) == DEFS.SHAPE_PREVIOUS) /* previous comp */
  { intersect = 1; }

  if (intersect)
  {
    if ((abs(Vars.Flag_Shape) == DEFS.SHAPE_SPHERE) || (abs(Vars.Flag_Shape) == DEFS.SHAPE_CYLIND)
     || (abs(Vars.Flag_Shape) == DEFS.SHAPE_BOX) || (abs(Vars.Flag_Shape) == DEFS.SHAPE_BANANA)
     || (geometry && strlen(geometry) && strcmp(geometry,"0") && strcmp(geometry, "NULL")) )
    {
      /* check if we have to remove the top/bottom with BANANA shape */
      if ((abs(Vars.Flag_Shape) == DEFS.SHAPE_BANANA) && (intersect != 1)) {
        double y0,y1;
        /* propagate to intersection point as temporary variable to check top/bottom */
        y0 = y+t0*vy;
        y1 = y+t1*vy;
        if (fabs(y0) >= Vars.Cylinder_Height/2*0.99) t0 = t1;
        if (fabs(y1) >= Vars.Cylinder_Height/2*0.99) t1 = t0;
      }
      if (t0 < 0 && t1 > 0)
        t0 = t;  /* neutron was already inside ! */
      if (t1 < 0 && t0 > 0) /* neutron exit before entering !! */
        t1 = t;
      /* t0 is now time of incoming intersection with the detection area */
      if ((Vars.Flag_Shape < 0) && (t1 > 0))
        PROP_DT(t1); /* t1 outgoing beam */
      else
        PROP_DT(t0); /* t0 incoming beam */
      /* Final test if we are on lid / bottom of banana/sphere */
      if (abs(Vars.Flag_Shape) == DEFS.SHAPE_BANANA || abs(Vars.Flag_Shape) == DEFS.SHAPE_SPHERE) {
        if (fabs(y) >= Vars.Cylinder_Height/2*0.99) {
          intersect=0;
          Flag_Restore=1;
        }
      }
    }
  }

  if (intersect)
  {
    /* Now get the data to monitor: current or keep from PreMonitor */
    if (Vars.Flag_UsePreMonitor != 1)
    {
      Vars.cp  = p;
      Vars.cx  = x;
      Vars.cvx = vx;
      Vars.csx = sx;
      Vars.cy  = y;
      Vars.cvy = vy;
      Vars.csy = sy;
      Vars.cz  = z;
      Vars.cvz = vz;
      Vars.csz = sz;
      Vars.ct  = t;
    }

    if ((Vars.He3_pressure > 0) && (t1 != t0) && ((abs(Vars.Flag_Shape) == DEFS.SHAPE_SPHERE) || (abs(Vars.Flag_Shape) == DEFS.SHAPE_CYLIND) || (abs(Vars.Flag_Shape) == DEFS.SHAPE_BOX)))
    {
      XY = exp(-7.417*Vars.He3_pressure*fabs(t1-t0)*2*PI*K2V);
      /* will monitor the absorbed part */
      Vars.cp *= 1-XY;
      /* and modify the neutron weight after monitor, only remains 1-p_detect */
      p *= XY;
    }

    if (Vars.Flag_capture)
    {
      XY = sqrt(Vars.cvx*Vars.cvx+Vars.cvy*Vars.cvy+Vars.cvz*Vars.cvz);
      XY *= V2K;
      if (XY != 0) XY = 2*PI/XY; /* lambda. lambda(2200 m/2) = 1.7985 Angs  */
      Vars.cp *= XY/1.7985;
    }

    pp = Monitor_nD_Trace(&DEFS, &Vars);
    if (pp==0.0)
    { ABSORB;
    }
    else if(pp==1)
    {
      SCATTER;
    }

    if (Vars.Flag_parallel) /* back to neutron state before detection */
      Flag_Restore = 1;
  } /* end if intersection */
  else {
    if (Vars.Flag_Absorb && !Vars.Flag_parallel)
    {
      // restore neutron ray before absorbing for correct mcdisplay
      RESTORE_NEUTRON(INDEX_CURRENT_COMP, x, y, z, vx, vy, vz, t, sx, sy, sz, p);
      ABSORB;
    }
    else Flag_Restore = 1;  /* no intersection, back to previous state */
  }

  if (Flag_Restore)
  {
    RESTORE_NEUTRON(INDEX_CURRENT_COMP, x, y, z, vx, vy, vz, t, sx, sy, sz, p);
  }
  #undef user1
  #undef user2
  #undef user3
  #undef xwidth
  #undef yheight
  #undef zdepth
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef zmin
  #undef zmax
  #undef bins
  #undef min
  #undef max
  #undef restore_neutron
  #undef radius
  #undef options
  #undef filename
  #undef geometry
  #undef username1
  #undef username2
  #undef username3
  #undef DEFS
  #undef Vars
  #undef detector
  #undef offdata
  return(_comp);
} /* class_Monitor_nD_trace */

#pragma acc routine seq
_class_DiskChopper *class_DiskChopper_trace(_class_DiskChopper *_comp
  , _class_particle *_particle) {
  ABSORBED=SCATTERED=RESTORE=0;

  #define theta_0 (_comp->_parameters.theta_0)
  #define radius (_comp->_parameters.radius)
  #define yheight (_comp->_parameters.yheight)
  #define nu (_comp->_parameters.nu)
  #define nslit (_comp->_parameters.nslit)
  #define jitter (_comp->_parameters.jitter)
  #define delay (_comp->_parameters.delay)
  #define isfirst (_comp->_parameters.isfirst)
  #define n_pulse (_comp->_parameters.n_pulse)
  #define abs_out (_comp->_parameters.abs_out)
  #define phase (_comp->_parameters.phase)
  #define xwidth (_comp->_parameters.xwidth)
  #define verbose (_comp->_parameters.verbose)
  #define Tg (_comp->_parameters.Tg)
  #define To (_comp->_parameters.To)
  #define delta_y (_comp->_parameters.delta_y)
  #define height (_comp->_parameters.height)
  #define omega (_comp->_parameters.omega)
    double toff;
    double yprime;
    PROP_Z0;
    yprime = y+delta_y;

    /* Is neutron outside the vertical slit range and should we absorb? */
    if (abs_out && (x*x+yprime*yprime)>radius*radius) {
      ABSORB;
    }
    /* Does neutron hit inner solid part of chopper in case of yheight!=radius? */
    if ((x*x+yprime*yprime)<(radius-height)*(radius-height)) {
      ABSORB;
    }


    if (isfirst)
      {
        /* all events are put in the transmitted time frame */
        t=atan2(x,yprime)/omega + To*randpm1()/2.0 + delay + (jitter ? jitter*randnorm():0) + (n_pulse > 1 ? floor(n_pulse*rand01())*Tg : 0);
        /* correction: chopper slits transmission opening/full disk */
        p *= nslit*theta_0/2.0/PI;
      }
    else
      {
        toff=fabs(t-atan2(x,yprime)/omega - delay - (jitter ? jitter*randnorm():0));

        /* does neutron hit outside slit? */
        if (fmod(toff+To/2.0,Tg)>To) ABSORB;
      }
    SCATTER;

  #undef theta_0
  #undef radius
  #undef yheight
  #undef nu
  #undef nslit
  #undef jitter
  #undef delay
  #undef isfirst
  #undef n_pulse
  #undef abs_out
  #undef phase
  #undef xwidth
  #undef verbose
  #undef Tg
  #undef To
  #undef delta_y
  #undef height
  #undef omega
  return(_comp);
} /* class_DiskChopper_trace */

#pragma acc routine seq
_class_Isotropic_Sqw *class_Isotropic_Sqw_trace(_class_Isotropic_Sqw *_comp
  , _class_particle *_particle) {
  ABSORBED=SCATTERED=RESTORE=0;

  #define powder_format (_comp->_parameters.powder_format)
  #define Sqw_coh (_comp->_parameters.Sqw_coh)
  #define Sqw_inc (_comp->_parameters.Sqw_inc)
  #define geometry (_comp->_parameters.geometry)
  #define radius (_comp->_parameters.radius)
  #define thickness (_comp->_parameters.thickness)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define zdepth (_comp->_parameters.zdepth)
  #define threshold (_comp->_parameters.threshold)
  #define order (_comp->_parameters.order)
  #define T (_comp->_parameters.T)
  #define verbose (_comp->_parameters.verbose)
  #define d_phi (_comp->_parameters.d_phi)
  #define concentric (_comp->_parameters.concentric)
  #define rho (_comp->_parameters.rho)
  #define sigma_abs (_comp->_parameters.sigma_abs)
  #define sigma_coh (_comp->_parameters.sigma_coh)
  #define sigma_inc (_comp->_parameters.sigma_inc)
  #define classical (_comp->_parameters.classical)
  #define powder_Dd (_comp->_parameters.powder_Dd)
  #define powder_DW (_comp->_parameters.powder_DW)
  #define powder_Vc (_comp->_parameters.powder_Vc)
  #define density (_comp->_parameters.density)
  #define weight (_comp->_parameters.weight)
  #define p_interact (_comp->_parameters.p_interact)
  #define norm (_comp->_parameters.norm)
  #define powder_barns (_comp->_parameters.powder_barns)
  #define quantum_correction (_comp->_parameters.quantum_correction)
  #define VarSqw (_comp->_parameters.VarSqw)
  #define columns (_comp->_parameters.columns)
  #define offdata (_comp->_parameters.offdata)

int    intersect=0;     /* flag to continue/stop */
double t0,  t1,  t2,  t3; /* times for intersections */
double dt0, dt1, dt2, dt; /* time intervals */
double k=0, Ei=0;
double v=0, vf=0;
double d_path;        /* total path length for straight trajectory */
double my_a;          /* absorption cross-section scaled to velocity (2200) */
double ws, p_scatt;   /* probability for scattering/absorption and for */
                      /* interaction along d_path */
double tmp_rand;      /* temporary var */
double ratio_w=0, ratio_q=0; /* variables for bilinear interpolation */
double q11, q21, q22, q12;
double omega=0;       /* energy transfer */
double q=0;           /* wavevector transfer */
long   index_w;       /* energy index for table look-up SW */
long   index_q;       /* Q index for table look-up P(Q|w) */
double theta=0, costheta=0; /* for the choice of kf direction */
double u1x,u1y,u1z;
double u2x,u2y,u2z;
double u0x,u0y,u0z;
int    index_counter;
int    flag=0;
int    flag_concentric=0;
int    flag_ishollow=0;
double solid_angle=0;
double my_t=0;
double p_mult=1;
double mc_trans, p_trans, mc_scatt;
double coh=0, inc=0;
struct Sqw_Data_struct Data_sqw;


/* Store Initial neutron state */

VarSqw.ki_x = V2K*vx;
VarSqw.ki_y = V2K*vy;
VarSqw.ki_z = V2K*vz;
VarSqw.ti   = t;
VarSqw.vi   = 0;
VarSqw.ki   = 0;
VarSqw.type = '\0';

do { /* Main interaction loop. Ends with intersect=0 */

  /* Intersection neutron trajectory / sample (sample surface) */
  if (VarSqw.s_coh > 0 || VarSqw.s_inc > 0) {
    if (thickness >= 0) {
      if (VarSqw.shape==0)
        intersect=cylinder_intersect(&t0,&t3, x,y,z,vx,vy,vz, radius,yheight);
      else if (VarSqw.shape==1)
        intersect=box_intersect     (&t0,&t3, x,y,z,vx,vy,vz, xwidth,yheight,zdepth);
      else if (VarSqw.shape==2)
        intersect=sphere_intersect  (&t0,&t3, x,y,z,vx,vy,vz, radius);
      else if (VarSqw.shape == 3)
        intersect=off_intersect(&t0, &t3, NULL, NULL, x, y, z, vx, vy, vz, offdata );
    } else {
      if (VarSqw.shape==0)
        intersect=cylinder_intersect(&t0,&t3, x,y,z,vx,vy,vz, radius-thickness,
          yheight-2*thickness > 0 ? yheight-2*thickness : yheight);
      else if (VarSqw.shape==1)
        intersect=box_intersect     (&t0,&t3, x,y,z,vx,vy,vz,
          xwidth-2*thickness > 0 ?  xwidth-2*thickness : xwidth,
          yheight-2*thickness > 0 ? yheight-2*thickness : yheight,
          zdepth-2*thickness > 0 ?  zdepth-2*thickness : zdepth);
      else if (VarSqw.shape==2)
        intersect=sphere_intersect  (&t0,&t3, x,y,z,vx,vy,vz, radius-thickness);
      else if (VarSqw.shape == 3)
        intersect=off_intersect(&t0, &t3, NULL, NULL, x, y, z, vx, vy, vz, offdata );
    }
  } else intersect=0;

  /* Computing the intermediate times */
  if (intersect) {
    flag_ishollow = 0;
    if (thickness > 0) {
      if (VarSqw.shape==0 && cylinder_intersect(&t1,&t2, x,y,z,vx,vy,vz, radius-thickness,
        yheight-2*thickness > 0 ? yheight-2*thickness : yheight))
        flag_ishollow=1;
      else if (VarSqw.shape==2 && sphere_intersect   (&t1,&t2, x,y,z,vx,vy,vz, radius-thickness))
        flag_ishollow=1;
      else if (VarSqw.shape==1 && box_intersect(&t1,&t2, x,y,z,vx,vy,vz,
        xwidth-2*thickness > 0 ? xwidth-2*thickness : xwidth,
        yheight-2*thickness > 0 ? yheight-2*thickness : yheight,
        zdepth-2*thickness > 0 ? zdepth-2*thickness : zdepth))
        flag_ishollow=1;
    } else if (thickness<0) {
      if (VarSqw.shape==0 && cylinder_intersect(&t1,&t2, x,y,z,vx,vy,vz, radius,yheight))
        flag_ishollow=1;
      else if (VarSqw.shape==2 && sphere_intersect   (&t1,&t2, x,y,z,vx,vy,vz, radius))
        flag_ishollow=1;
      else if (VarSqw.shape==1 && box_intersect(&t1,&t2, x,y,z,vx,vy,vz, xwidth, yheight, zdepth))
        flag_ishollow=1;
    }
    if (!flag_ishollow) t1 = t2 = t3; /* no empty space inside */
  } else break; /* neutron does not hit sample: transmitted  */

  if (intersect) { /* the neutron hits the sample */

    if (t0 > 0) {  /* we are before the sample */
      PROP_DT(t0); /* propagates neutron to the entry of the sample */
    } else if (t1 > 0 && t1 > t0) { /* we are inside first part of the sample */
      /* no propagation, stay inside */
    } else if (t2 > 0 && t2 > t1) { /* we are in the hole */
      PROP_DT(t2); /* propagate to inner surface of 2nd part of sample */
    } else if (t3 > 0 && t3 > t2) { /* we are in the 2nd part of sample */
      /* no propagation, stay inside */
    }

    dt0=t1-(t0 > 0 ? t0 : 0); /* Time in first part of hollow/cylinder/box */
    dt1=t2-(t1 > 0 ? t1 : 0); /* Time in hole */
    dt2=t3-(t2 > 0 ? t2 : 0); /* Time in 2nd part of hollow cylinder */

    if (dt0 < 0) dt0 = 0;
    if (dt1 < 0) dt1 = 0;
    if (dt2 < 0) dt2 = 0;

    /* initialize concentric mode */
    if (concentric && !flag_concentric && t0 >= 0
     && VarSqw.shape==0 && thickness) {
      flag_concentric=1;
    }

    if (flag_concentric == 1) {
      dt1=dt2=0; /* force exit when reaching hole/2nd part */
    }

    if (!dt0 && !dt2) {
      intersect = 0; /* the sample was passed entirely */
      break;
    }

    VarSqw.neutron_enter++;
    p_mult = 1;
    if (!v) {
      v  = vx*vx+vy*vy+vz*vz;
      v = sqrt(v);
    }
    k  = V2K*v;
    Ei = VS2E*v*v;

    if (!VarSqw.vi) VarSqw.vi = v;
    if (!VarSqw.ki) VarSqw.ki = k;

    if (v <= 0) {
      printf("Isotropic_Sqw: %s: ERROR: Null velocity !\n",NAME_CURRENT_COMP);
      VarSqw.neutron_removed++;
      ABSORB; /* should never occur */
    }

    /* check for scattering event */
    my_a   = VarSqw.my_a_v / v; /* absorption 'mu' */
    /* compute total scattering X section */
    /* \int q S(q) dq /2 /ki^2 sigma  OR  bare Xsection*/
    /* contains the 4*PI*kf/ki factor */
    coh = VarSqw.s_coh;
    inc = VarSqw.s_inc;
    if (k && VarSqw.s_coh>0 && VarSqw.Data_coh.intensity) {
      double Ei       = VS2E*v*v;
      double index_Ei = Ei / (VarSqw.Data_coh.Ei_max/VarSqw.Data_coh.iqSq_length);
      coh = Table_Value2d(VarSqw.Data_coh.iqSq, index_Ei, 0);
    }
    if (k && VarSqw.s_inc>0 && VarSqw.Data_inc.intensity) {
      double Ei       = VS2E*v*v;
      double index_Ei = Ei / (VarSqw.Data_inc.Ei_max/VarSqw.Data_inc.iqSq_length);
      inc = Table_Value2d(VarSqw.Data_inc.iqSq, index_Ei, 0);
    }
    if (coh<0) coh=0;
    if (inc<0) inc=0;
    VarSqw.my_s    =(VarSqw.mat_rho*100*(coh + inc));

    my_t = my_a + VarSqw.my_s;  /* total scattering Xsect */
    if (my_t <= 0) {
      if (VarSqw.neutron_removed<VarSqw.maxloop) printf("Isotropic_Sqw: %s: ERROR: Null total cross section %g. Removing event.\n",
        NAME_CURRENT_COMP, my_t);
      VarSqw.neutron_removed++;
      ABSORB; /* should never occur */
    } else if (VarSqw.my_s <= 0) {
      if (VarSqw.verbose_output > 1 && VarSqw.neutron_removed<VarSqw.maxloop)
        printf("Isotropic_Sqw: %s: Warning: Null scattering cross section %g. Ignoring.\n",
          NAME_CURRENT_COMP, VarSqw.my_s);
      VarSqw.my_s = 0;
    }

    /* Proba of scattering vs absorption (integrating along the whole trajectory) */
    ws = VarSqw.my_s/my_t;  /* (inc+coh)/(inc+coh+abs) */
    d_path = v*( dt0 +dt2 );    /* total path lenght in sample */
    /* Proba of transmission/interaction along length d_path */
    p_trans = exp(-my_t*d_path);
    p_scatt = 1 - p_trans; /* portion of beam which scatters */

    flag = 0; /* flag used for propagation to exit point before ending */

    /* are we next to the exit ? probably no scattering (avoid rounding errors) */
    if (VarSqw.my_s*d_path <= 4e-7) {
      flag = 1;           /* No interaction before the exit */
    }
    /* force a given fraction of the beam to scatter */
    if (p_interact>0 && p_interact<=1) {
      /* we force a portion of the beam to interact */
      /* This is used to improve statistics on single scattering (and multiple) */
      if (!SCATTERED) mc_trans = 1-p_interact;
      else            mc_trans = 1-p_interact/(4*SCATTERED+1); /* reduce effect on multi scatt */
    } else {
      mc_trans = p_trans; /* 1 - p_scatt */
    }
    mc_scatt = 1 - mc_trans; /* portion of beam to scatter (or force to) */
    if (mc_scatt <= 0 || mc_scatt>1) flag=1;
    /* MC choice: Interaction or transmission ? */
    if (!flag && mc_scatt > 0 && (mc_scatt >= 1 || rand01() < mc_scatt)) { /* Interaction neutron/sample */
      p_mult *= ws; /* Update weight ; account for absorption and retain scattered fraction */
      /* we have chosen portion mc_scatt of beam instead of p_scatt, so we compensate */
      if (!mc_scatt) ABSORB;
      p_mult *= fabs(p_scatt/mc_scatt); /* lower than 1 */
    } else {
      flag = 1; /* Transmission : no interaction neutron/sample */
      if (!VarSqw.type) VarSqw.type = 't';
      if (!mc_trans) ABSORB;
      p_mult *= fabs(p_trans/mc_trans);  /* attenuate beam by portion which is scattered (and left along) */
    }

    if (flag) { /* propagate to exit of sample and finish */
      intersect = 0;
      p *= p_mult; /* apply absorption correction */
      PROP_DT(dt0+dt2);
      break; /* exit main multi scatt while loop */
    }
  } /* end if intersect the neutron hits the sample */
  else break;

  if (intersect) { /* scattering event */
    double kf=0, kf1, kf2;
    /* mean scattering probability and absorption fraction */
    VarSqw.mean_scatt += (1-exp(-VarSqw.my_s*d_path))*p;
    VarSqw.mean_abs   += (1-ws)*p;
    VarSqw.psum_scatt += p;

    /* Decaying exponential distribution of the path length before scattering */
    /* Select a point at which to scatter the neutron, taking
         secondary extinction into account. */
    if (my_t*d_path < 1e-6)
    /* For very weak scattering, use simple uniform sampling of scattering
       point to avoid rounding errors. */
      dt = rand0max(d_path); /* length */
    else
      dt = -log(1 - rand0max((1 - exp(-my_t*d_path)))) / my_t; /* length */
    dt /= v; /* Time from present position to scattering point */

    /* If t0 is in hole, propagate to next part of the hollow cylinder */
    if (dt1 > 0 && dt0 > 0 && dt > dt0) dt += dt1;

    /* Neutron propagation to the scattering point */
    PROP_DT(dt);

    /* choice between coherent/incoherent scattering */
    tmp_rand = rand01();
    /* local description at the scattering point (scat probability for atom) */
    tmp_rand *= (coh+inc);

    flag=0;
    if (VarSqw.s_inc>0 && tmp_rand < inc) {
      /* CASE 1: incoherent case */
      if (!VarSqw.Data_inc.intensity) {
        /* CASE 1a: no incoherent Sqw from file, use isotropic V-like */
        if (d_phi && order == 1) {
          randvec_target_rect_angular(&u1x, &u1y, &u1z, &solid_angle,
              vx, vy, vz, 2*PI, d_phi, ROT_A_CURRENT_COMP);
          p_mult *= solid_angle/4/PI; /* weighted by focused range to total range */
        } else
          randvec_target_circle(&u1x, &u1y, &u1z, NULL, vx, vy, vz, 0);

        vx = u1x; vy = u1y; vz = u1z;
        vf = v; kf = k;
        if (!VarSqw.type) VarSqw.type = 'v';
        SCATTER;
      } else {
        /* CASE 1b: incoherent Sqw from file */
        if (VarSqw.Data_inc.intensity) {
          Data_sqw = VarSqw.Data_inc;
          if (!VarSqw.type) VarSqw.type = 'i';
          flag = 1;
        }
      }
    } else if (VarSqw.s_coh>0 && tmp_rand > VarSqw.s_inc) {
      if (VarSqw.Data_coh.intensity) {
        /* CASE2: coherent case */
        Data_sqw = VarSqw.Data_coh;
        if (!VarSqw.type) VarSqw.type = 'c';
        flag = 1;
      }
    }

    if (flag) { /* true when S(q,w) table exists (Data_sqw) */

      double alpha=0, alpha0;
      /* give us a limited number of tries for scattering: choose W then Q */
      for (index_counter=VarSqw.maxloop; index_counter > 0 ; index_counter--) {

        /* MC choice: energy transfer w=Ei-Ef in the S(w) = SW */
        omega = 0;
        tmp_rand = rand01();
        /* energy index for rand > cumul SW */
        index_w  = Sqw_search_SW(Data_sqw, tmp_rand);
        VarSqw.rw = (double)index_w;
        if (index_w >= 0 && &(Data_sqw.SW[index_w]) != NULL) {
          if (Data_sqw.w_bins > 1) {
            double w1, w2;
            if (index_w > 0) { /* interpolate linearly energy */
              ratio_w = (tmp_rand                         - Data_sqw.SW[index_w-1].cumul_proba)
                       /(Data_sqw.SW[index_w].cumul_proba - Data_sqw.SW[index_w-1].cumul_proba);
              /* ratio_w=0 omega[index_w-1], ratio=1 omega[index] */
              w1 = Data_sqw.SW[index_w-1].omega; w2 = Data_sqw.SW[index_w].omega;
            } else { /* index_w = 0 interpolate to 0 energy */
              /* ratio_w=0 omega=0, ratio=1 omega[index] */
              w1 = Data_sqw.SW[index_w].omega; w2= Data_sqw.SW[index_w+1].omega;
              if (!w2 && index_w+1 < Data_sqw.w_bins)
                w2= Data_sqw.SW[index_w+1].omega;
              if (Data_sqw.w_bins && Data_sqw.SW[index_w].cumul_proba) {
                ratio_w = tmp_rand/Data_sqw.SW[index_w].cumul_proba;
              } else ratio_w=0;
            }
            if (ratio_w<0) ratio_w=0; else if (ratio_w>1) ratio_w=1;
            omega = (1-ratio_w)*w1 + ratio_w*w2;
          } else {
            ratio_w = 0;
            omega = Data_sqw.SW[index_w].omega;
          }
        } else {
          if (VarSqw.verbose_output >= 3 && VarSqw.neutron_removed<VarSqw.maxloop)
            printf("Isotropic_Sqw: %s: Warning: No suitable w transfer for index_w=%li.\n",
              NAME_CURRENT_COMP, index_w);
          continue; /* no W value: try again with an other energy transfer */
        }

        /* MC choice: momentum transfer Q in P(Q|w) */
        tmp_rand = rand01();

        /* momentum index for rand > cumul SQ|W */
        index_q  = Sqw_search_Q_proba_per_w(Data_sqw, tmp_rand, index_w);
        VarSqw.rq = (double)index_q;

        if (index_q >= 0 && &(Data_sqw.SQW[index_w]) != NULL) {
          if (Data_sqw.q_bins > 1 && index_q > 0) {
            if (index_w > 0 && Data_sqw.w_bins > 1) {
              /* bilinear interpolation on - side: index_w > 0, index_q > 0 */
              ratio_q = (tmp_rand - Data_sqw.SQW[index_w][index_q-1].cumul_proba)
                       /(Data_sqw.SQW[index_w][index_q].cumul_proba
                       - Data_sqw.SQW[index_w][index_q-1].cumul_proba);
              q22 = Data_sqw.SQW[index_w]  [index_q].Q;
              q11 = Data_sqw.SQW[index_w-1][index_q-1].Q;
              q21 = Data_sqw.SQW[index_w]  [index_q-1].Q;
              q12 = Data_sqw.SQW[index_w-1][index_q].Q;
              if (ratio_q<0) ratio_q=0; else if (ratio_q>1) ratio_q=1;
              q = (1-ratio_w)*(1-ratio_q)*q11+ratio_w*(1-ratio_q)*q21
                + ratio_w*ratio_q*q22        +(1-ratio_w)*ratio_q*q12;
            } else { /* bilinear interpolation on + side: index_w=0, index_q > 0 */
              ratio_q = (tmp_rand - Data_sqw.SQW[index_w][index_q-1].cumul_proba)
                       /(Data_sqw.SQW[index_w][index_q].cumul_proba
                       - Data_sqw.SQW[index_w][index_q-1].cumul_proba);
              q11 = Data_sqw.SQW[index_w]  [index_q-1].Q;
              q12 = Data_sqw.SQW[index_w]  [index_q].Q;
              if (ratio_q<0) ratio_q=0; else if (ratio_q>1) ratio_q=1;
              if (index_w < Data_sqw.w_bins-1 && Data_sqw.w_bins > 1) {
                q22 = Data_sqw.SQW[index_w+1][index_q].Q;
                q21 = Data_sqw.SQW[index_w+1][index_q-1].Q;
                q = (1-ratio_w)*(1-ratio_q)*q11+ratio_w*(1-ratio_q)*q21
                  + ratio_w*ratio_q*q22        +(1-ratio_w)*ratio_q*q12;
              } else {
                q    = (1-ratio_q)*q11  + ratio_q*q12;
              }
            }
          } else {
            q    = Data_sqw.SQW[index_w][index_q].Q;
          }
        } else {
          if (VarSqw.verbose_output >= 3 && VarSqw.neutron_removed<VarSqw.maxloop)
            printf("Isotropic_Sqw: %s: Warning: No suitable q transfer for w=%g.\n",
              NAME_CURRENT_COMP, omega);
          VarSqw.neutron_removed++;
          continue; /* no Q value for this w choice */
        }

        /* Search for length of final wave vector kf */
        /* kf is such that : hbar*w = hbar*hbar/2/m*(k*k - kf*kf) */
        /* acceptable values for kf are kf1 and kf2 */
        if (!solve_2nd_order(&kf1, &kf2, 1, 0, -k*k + VarSqw.sqSE2K*omega)) {
          if (VarSqw.verbose_output >= 3 && VarSqw.neutron_removed<VarSqw.maxloop)
            printf("Isotropic_Sqw: %s: Warning: imaginary root for w=%g q=%g Ei=%g (triangle can not close)\n",
            NAME_CURRENT_COMP, omega, q, Ei);
          VarSqw.neutron_removed++;
          continue; /* all roots are imaginary */
        }

        /* kf1 and kf2 are opposite */
        kf = fabs(kf1);
        vf = K2V*kf;

        /* Search of the direction of kf such that : q = ki - kf */
        /* cos theta = (ki2+kf2-q2)/(2ki kf) */

        costheta= (k*k+kf*kf-q*q)/(2*kf*k); /* this is cos(theta) */

        if (-1 < costheta && costheta < 1) {
          break; /* satisfies q momentum conservation */
        }
/*      else continue; */

        /* exit for loop on success */
      } /* end for index_counter */

      if (!index_counter) { /* for loop ended: failure for scattering */
        intersect=0; /* Could not scatter: finish multiple scattering loop */
        if (VarSqw.verbose_output >= 2 && VarSqw.neutron_removed<VarSqw.maxloop)
          printf("Isotropic_Sqw: %s: Warning: No scattering [q,w] conditions\n"
               "               last try (%i): type=%c w=%g q=%g cos(theta)=%g k=%g\n",
          NAME_CURRENT_COMP, VarSqw.maxloop, (VarSqw.type ? VarSqw.type : '-'), omega, q, costheta, k);
        VarSqw.neutron_removed++;
        if (order && SCATTERED != order) ABSORB;
        break;       /* finish multiple scattering loop */
      }

      /* scattering angle from ki to DS cone */
      theta = acos(costheta);

      /* Choose point on Debye-Scherrer cone */
      if (order == 1 && d_phi)
      { /* relate height of detector to the height on DS cone */
        double cone_focus;
        cone_focus = sin(d_phi/2)/sin(theta);
        /* If full Debye-Scherrer cone is within d_phi, don't focus */
        if (cone_focus < -1 || cone_focus > 1) d_phi = 0;
        /* Otherwise, determine alpha to rotate from scattering plane
            into d_phi focusing area*/
        else alpha = 2*asin(cone_focus);
        if (d_phi) p_mult *= alpha/PI;
      }
      if (d_phi) {
        /* Focusing */
        alpha = fabs(alpha);
        /* Trick to get scattering for pos/neg theta's */
        alpha0= 2*rand01()*alpha;
        if (alpha0 > alpha) {
          alpha0=PI+(alpha0-1.5*alpha);
        } else {
          alpha0=alpha0-0.5*alpha;
        }
      }
      else
        alpha0 = PI*randpm1();

      /* now find a nearly vertical rotation axis (u1) :
	       * Either
	       *  (v along Z) x (X axis) -> nearly Y axis
	       * Or
	       *  (v along X) x (Z axis) -> nearly Y axis
	       */
	    if (fabs(scalar_prod(1,0,0,vx/v,vy/v,vz/v)) < fabs(scalar_prod(0,0,1,vx/v,vy/v,vz/v))) {
        u1x = 1; u1y = u1z = 0;
	    } else {
        u1x = u1y = 0; u1z = 1;
	    }
	    vec_prod(u2x,u2y,u2z, vx,vy,vz, u1x,u1y,u1z);

      /* handle case where v and aim are parallel */
      if (!u2x && !u2y && !u2z) { u2x=u2z=0; u2y=1; }

      /* u1 = rotate 'v' by theta around u2: DS scattering angle, nearly in horz plane */
      rotate(u1x,u1y,u1z, vx,vy,vz, theta, u2x,u2y,u2z);

      /* u0 = rotate u1 by alpha0 around 'v' (Debye-Scherrer cone) */
      rotate(u0x,u0y,u0z, u1x,u1y,u1z, alpha0, vx, vy, vz);
      NORM(u0x,u0y,u0z);
      vx = u0x*vf;
      vy = u0y*vf;
      vz = u0z*vf;

      SCATTER;

      v = vf; k = kf; /* for next iteration */

    } /* end if (flag) */

    VarSqw.neutron_exit++;
    p *= p_mult;
    if (p_mult > 1) VarSqw.neutron_pmult++;

    /* test for a given multiple order */
    if (order && SCATTERED >= order) {
      intersect=0; /* reached required number of SCATTERing */
      break;       /* finish multiple scattering loop */
    }

  } /* end if (intersect) scattering event  */

} while (intersect); /* end do (intersect) (multiple scattering loop) */

/* Store Final neutron state */
VarSqw.kf_x = V2K*vx;
VarSqw.kf_y = V2K*vy;
VarSqw.kf_z = V2K*vz;
VarSqw.tf   = t;
VarSqw.vf   = v;
VarSqw.kf   = k;
VarSqw.theta= theta;

if (SCATTERED) {



  if (SCATTERED == 1) {
    if (VarSqw.type == 'c') VarSqw.single_coh += p;
    else                    VarSqw.single_inc += p;
    VarSqw.dq = sqrt((VarSqw.kf_x-VarSqw.ki_x)*(VarSqw.kf_x-VarSqw.ki_x)
                  +(VarSqw.kf_y-VarSqw.ki_y)*(VarSqw.kf_y-VarSqw.ki_y)
                  +(VarSqw.kf_z-VarSqw.ki_z)*(VarSqw.kf_z-VarSqw.ki_z));
    VarSqw.dw = VS2E*(VarSqw.vf*VarSqw.vf - VarSqw.vi*VarSqw.vi);
  } else VarSqw.multi += p;

} else VarSqw.dq=VarSqw.dw=0;

/* end TRACE */

  // EXTEND code here
  if (!strcmp(_comp->_name, "Sample")) {
  if (!SCATTERED) ABSORB;
  }

  #undef powder_format
  #undef Sqw_coh
  #undef Sqw_inc
  #undef geometry
  #undef radius
  #undef thickness
  #undef xwidth
  #undef yheight
  #undef zdepth
  #undef threshold
  #undef order
  #undef T
  #undef verbose
  #undef d_phi
  #undef concentric
  #undef rho
  #undef sigma_abs
  #undef sigma_coh
  #undef sigma_inc
  #undef classical
  #undef powder_Dd
  #undef powder_DW
  #undef powder_Vc
  #undef density
  #undef weight
  #undef p_interact
  #undef norm
  #undef powder_barns
  #undef quantum_correction
  #undef VarSqw
  #undef columns
  #undef offdata
  return(_comp);
} /* class_Isotropic_Sqw_trace */

#pragma acc routine seq
_class_PSD_Detector *class_PSD_Detector_trace(_class_PSD_Detector *_comp
  , _class_particle *_particle) {
  ABSORBED=SCATTERED=RESTORE=0;

  #define nx (_comp->_parameters.nx)
  #define ny (_comp->_parameters.ny)
  #define xwidth (_comp->_parameters.xwidth)
  #define radius (_comp->_parameters.radius)
  #define awidth (_comp->_parameters.awidth)
  #define yheight (_comp->_parameters.yheight)
  #define zdepth (_comp->_parameters.zdepth)
  #define threshold (_comp->_parameters.threshold)
  #define PressureConv (_comp->_parameters.PressureConv)
  #define PressureStop (_comp->_parameters.PressureStop)
  #define interpolate (_comp->_parameters.interpolate)
  #define p_interact (_comp->_parameters.p_interact)
  #define verbose (_comp->_parameters.verbose)
  #define LensOn (_comp->_parameters.LensOn)
  #define dc (_comp->_parameters.dc)
  #define borderx (_comp->_parameters.borderx)
  #define bordery (_comp->_parameters.bordery)
  #define xChDivRelSigma (_comp->_parameters.xChDivRelSigma)
  #define yChDivRelSigma (_comp->_parameters.yChDivRelSigma)
  #define bufsize (_comp->_parameters.bufsize)
  #define restore_neutron (_comp->_parameters.restore_neutron)
  #define angle (_comp->_parameters.angle)
  #define type (_comp->_parameters.type)
  #define filename (_comp->_parameters.filename)
  #define FN_Conv (_comp->_parameters.FN_Conv)
  #define FN_Stop (_comp->_parameters.FN_Stop)
  #define PSD_N (_comp->_parameters.PSD_N)
  #define PSD_p (_comp->_parameters.PSD_p)
  #define PSD_p2 (_comp->_parameters.PSD_p2)
  #define EAP (_comp->_parameters.EAP)
  #define M1P1 (_comp->_parameters.M1P1)
  #define PosAP (_comp->_parameters.PosAP)
  #define EAT (_comp->_parameters.EAT)
  #define M1T1 (_comp->_parameters.M1T1)
  #define PosAT (_comp->_parameters.PosAT)
  #define CrossSectionHe (_comp->_parameters.CrossSectionHe)
  #define CountNeutrons (_comp->_parameters.CountNeutrons)
  #define GeomCumul (_comp->_parameters.GeomCumul)
  #define AbsCumul (_comp->_parameters.AbsCumul)
  #define SensVolCumul (_comp->_parameters.SensVolCumul)
  #define DetCumul (_comp->_parameters.DetCumul)
  #define PHSpectrum (_comp->_parameters.PHSpectrum)
  #define PHSpectrum0 (_comp->_parameters.PHSpectrum0)
  #define PHSpectrum2 (_comp->_parameters.PHSpectrum2)
  #define PHSpectrum_n (_comp->_parameters.PHSpectrum_n)
  #define nH_p (_comp->_parameters.nH_p)
  #define nH_t (_comp->_parameters.nH_t)
  #define FullEnergyP (_comp->_parameters.FullEnergyP)
  #define FullEnergyT (_comp->_parameters.FullEnergyT)
  #define VariousErrors (_comp->_parameters.VariousErrors)
  #define DetectorType (_comp->_parameters.DetectorType)
  #define DEFS (_comp->_parameters.DEFS)
  #define Vars (_comp->_parameters.Vars)

  long   i,j;
  double RangeProton,RangeTriton;
  double rb=radius+zdepth; // radius of the back plane of the detector
  double NAvogadro,MolVolume,v,vt,mu;
  double length,GM,ZB,pa,RandomNumber,la,ta;
  double phi,theta;   // direction of emission of the particles
  double vxp,vyp,vzp; // normalized speed of the particles
  double mlp,mlt;     // 1st moment (COG) of energy deposition of the particles
  double xCOG,yCOG,zCOG,aCOG,rCOG,Energy;
  double tin,tout,tine,toute,tinb,toutb; // times in and out of cylindrical entrance window and back plane boxes
  long   ISe,ISb; // type of intersection with a cylindrical box
  double EnergyP,EnergyT,SigmaEnergy;
  double lp,lt; // length proton- and triton tracks
  long   intersect;
  long   Impossibility; // used to check for impossible collisions
  double alimit; // borders of the sensitive volume + edges
  double aa; // alpha coordinate (along the arc of the detector) of the absorption
  double t1,t2; // time for the first and second charged particles to the next wall

  CountNeutrons += p;

  RangeProton = PosAP[nH_p];
  RangeTriton = PosAT[nH_t];
  v = sqrt(vx*vx + vy*vy + vz*vz) ; // speed of the neutron
  if (v>440000 && !(VariousErrors&1) ) {
      fprintf(stderr,"PSD_Detector: %s: Component cannot deal with\n"
                     "              high-energy neutrons. Absorbing.\n",
                                    NAME_CURRENT_COMP);
      VariousErrors+=1;
      ABSORB; }
  vt = sqrt( 2*0.025243*1.60218e-19 / ( 1.0086649 * 1.66053886e-27 ) ); // speed of 25.243 meV neutron (0.18 nm)
  NAvogadro = 6.022045e23; // Number of atoms per mol
  MolVolume = 24.7796e-3; // in m3, 24.7796 liter per mol at 1 bar at T=298 K
  mu = (NAvogadro*PressureConv*CrossSectionHe* vt ) / (MolVolume* v ); // in inverse m

  switch(DetectorType) {
    case 1:
      intersect = box_intersect(&tin,&tout,x,y,z,vx,vy,vz,xwidth+2*borderx,yheight+2*bordery,zdepth);
      if (tin==tout) intersect=0;
      break;
    case 2:
      intersect = cylinder_intersect(&tin,&tout,x,y,z,vx,vy,vz,radius,yheight+2*bordery);
      if (tin==tout) intersect=0;
      break;
    case 4:
      if (cylinder_intersect(&tine,&toute,x,y,z,vx,vy,vz,radius,yheight+2*bordery)==0 ) {
          ISe=4; } // case for no intersection with the entrance-window-cylinder
      else {
        if (tine<=0 && toute<=0) {
            ISe=1; // case for neutron outside the volume and moving away
        }
        if (tine<=0 && toute>0) {
            ISe=2; // case for neutron inside the volume and moving (obviously) out
        }
        if (tine>0 && toute>0) {
            ISe=3; // case for neutron outside the volume and moving towards it
        }
        if (tine>toute && !(VariousErrors&2) ) {
          fprintf(stderr,"PSD_Detector: %s: Something is seriously wrong.\n"
                         "              cylinder_intersect reported a later time\n"
                         "              for entering the cylinder than for leaving it.\n",
                                        NAME_CURRENT_COMP);
          VariousErrors+=2;
          ABSORB;
        }
      }
      if (cylinder_intersect(&tinb,&toutb,x,y,z,vx,vy,vz,rb,yheight+2*bordery)==0 || tinb==toutb) {
          ISb=4; } // case for no intersection with the back-plane-cylinder
      else {
        if (tinb<=0 && toutb<=0) {
            ISb=1; // case for neutron outside the volume and moving away
        }
        if (tinb<=0 && toutb>0) {
            ISb=2; // case for neutron inside the volume and moving (obviously) out
        }
        if (tinb>0 && toutb>0) {
            ISb=3; // case for neutron outside the volume and moving towards it
        }
        if (tinb>toutb && !(VariousErrors&2) ) {
          fprintf(stderr,"PSD_Detector: %s: Something is seriously wrong.\n"
                         "              cylinder_intersect reported a later time\n"
                         "              for entering the cylinder than for leaving it.\n",
                                        NAME_CURRENT_COMP);
          VariousErrors+=2;
          ABSORB;
        }
      }
    /********************************************************************
    * Schematic representations of the possibilities for the path of the neutron
    * -----------------------------------------------------------------
    *        |  ISb=1       ISb=2           ISb=3           ISb=4
    * -----------------------------------------------------------------
    * ISe=1  |  Gone        already in      Impossible      Impossible
    *        |              out thr. back
    * -----------------------------------------------------------------
    * ISe=2  |  Impossible  in thr. entr.   Impossible      Impossible
    *        |              out thr. back
    * -----------------------------------------------------------------
    * ISe=3  |  Impossible  already in      in thr. back    Impossible
    *        |              out thr. entr.  out thr. entr.
    * -----------------------------------------------------------------
    * ISe=4  |  Gone        already in      in thr. back    Gone
    *        |              out thr. back   out thr. back
    * -----------------------------------------------------------------
    * All we have to do is choose the appropriate entrance and exit times
    * from the tin and tout of both the entrance window cylinder and the
    * back plane cylinder.
    * There is a problem on the left and right sides of the detector.
    * We cannot check entrance or exit using cylinder_intersect there, so
    * we simply find interaction points in 2pi cylinder wall, and then
    * reject events that are too far from the sensitive volume.
    *********************************************************************/


      Impossibility=1; // set this to 1: if it is not subsequently set to 0, an impossible neutron track has occurred
      if ( ((ISb==1) && ((ISe==1)||(ISe==4))) || ((ISb==4)&&(ISe==4)) ) {
        Impossibility=0; // these represent the cases that the neutron misses the detector
        intersect=0;
      }
      if ( (ISe==1)&&(ISb==2) ) {
        Impossibility = 0;
        intersect=1;
        tin=0;     // neutron is already located within the sensitive volume
        tout=toutb; // neutron flies out the sensitive volume into the back plane
      }
      if ( (ISe==2)&&(ISb==2) ) {
        Impossibility = 0;
        intersect=1;
        tin=toute; // neutron comes in through the entrance window
        tout=toutb; // and flies out through the back plane
      }
      if ( (ISe==3)&&(ISb==2) ) {
        Impossibility = 0;
        intersect=1;
        tin=0;    // neutron is already located within the sensitive volume
        tout=tine; // neutron flies out the sensitive volume into the entrance window
      }
      if ( (ISe==3)&&(ISb==3) ) {
        Impossibility = 0;
        intersect=1;
        tin=tinb; // neutron comes in through the back plane
        tout=tine; // and flies out through the entrance window
      }
      if ( (ISe==4)&&(ISb==2) ) {
        Impossibility = 0;
        intersect=1;
        tin=0;     // neutron is already located within the sensitive volume
        tout=toutb; // neutron flies out the sensitive volume into the back plane
      }
      if ( (ISe==4)&&(ISb==3) ) {
        Impossibility = 0;
        intersect=1;
        tin=tinb;  // neutron comes in through the back plane
        tout=toutb; // and flies out through the back plane
      }
      if ( tin==tout ) {
        intersect=0; /* This happens for instance when the neutron flies along
        the axis of the cylinders. The two cylinders overlap there, so zero time
        is spent in the gas volume, therefore no absorption. */
      }
      if ( tin>tout && !(VariousErrors&4) ) {
        intersect=0;
        fprintf(stderr,"PSD_Detector: %s: A serious error occurred.\n"
                       "              A later time for entering the detector was calculated\n"
                       "              than for leaving it.\n",
                                      NAME_CURRENT_COMP);
        VariousErrors+=4;
      }
      if ( Impossibility==1 && !(VariousErrors&8) ) {
        fprintf(stderr,"PSD_Detector: %s: Something strange happened. A neutron\n"
                       "              followed a path deemed impossible in the algorithm.\n",
                                      NAME_CURRENT_COMP);
        intersect=0;
        VariousErrors+=8;
      }
    /* NB: at this point the neutron does not yet have to be absorbed
    within the sensitive volume of the detector - we still have to check
    polar angle to see if the event falls within the part of the cylinder
    wall that is covered by the detector. */


      break;  /* end switch case DetectorType==4 */
    default:
      exit(0);
  } // end switch DetectorType (1)

  if (intersect && tout<=0) {
    intersect=0;
    /* This is the case that the intersection along the trajectory was in the past,
    which is probably a common occurrance in a simulation, e.g. backscattering off
    the entrance window, and which should not be interpreted as a new detector event.*/
  }
  if (intersect && tin<=0 && !(VariousErrors&16) ) {
    tin=0;
    fprintf(stderr,"PSD_Detector: %s: Warning: a neutron has been found\n"
                   "              'tunneled' into the detector, rather than just entering\n"
                   "              through the entrance window or one of the other sides.\n",
                                  NAME_CURRENT_COMP);
    VariousErrors+=16;
  }

  if (intersect) {
    /* there is intersection, not 'grazing' of the detector.
    What kind of trajectory is followed is not important, since we have
    in- and out-going time. */
    length = v*(tout-tin); // length (in m) of the trajectory through the sensitive volume
    GM = mu*length ; // number of radiation lengths in the gap thickness
    ZB = exp( -GM ); // zone boundary for generating uniformly distributed random numbers [ZB,1]
    pa = 1 - ZB; // Probability of absorption somewhere along the length of the trajectory through the sensitive volume

    if ( (p_interact <= 0 && (pa >= 1  || rand01() < pa))
        || (p_interact > 0 && (p_interact >= 1 || rand01() < p_interact)) ) {
      if (p_interact > 0 && p_interact < 1) pa /= p_interact;
      if (p_interact <= 0) pa=1;
      RandomNumber = (1-ZB)*rand01() + ZB; // uniformly distributed random number between ZB and 1
      la = -log(RandomNumber) / mu; // (in m) Absorption location (along trajectory), distributed exponentially declining
      ta = tin + la/v; // Absorption time
      alimit= (0.5*awidth + borderx)/radius; /* we take the gas-filled volume a bit larger
        than the sensitive volume alone, which is realistic. Later we throw away events
        that do not end up in the sensitive volume. */
      aa=atan2(x+vx*ta,z+vz*ta);
      if ( DetectorType==1 || DetectorType==2 ||
            (DetectorType==4 && aa>-alimit && aa<alimit) ) {
        /* this is always executed for a box or a tube, but for a banana
          only when the interaction position is not too far left or right
          of the sensitive volume (in order to account for border effects) */
        GeomCumul += p/(1-ZB); /* Cumulative probability of neutrons to encounter the
            detector is incremented by the weight of the neutron. Corrected for
      the fact that this part of the algorithm is only executed with
      probability pa=1-ZB.
      By later dividing by the total sum of weights, we will find the
      'geometric efficiency'.
      This part of the algorithm cannot be executed before this point,
      because we are not sure yet if it falls within the detector. */
        AbsCumul += p; /* Cumulative probability of absorption. */
        PROP_DT(ta);
        SCATTER; // show point of creation of charges
        // select random direction in 4 PI
        theta = acos( rand01()*2-1 ); // polar angle, distributed like a sine from 0 to pi
        phi = rand01() *2*PI; // azimuth angle, uniformly distributed from 0 to 2pi

        vxp = sin(theta) * cos(phi); // unit vector, interpreted as normalized speed (1 m/s) of the emitted proton.
        vyp = sin(theta) * sin(phi); // this is used to obtain both proton and triton track end points.
        vzp = cos(theta);

        // check intersection of charge trajectory
        switch(DetectorType) {
          case 1:
            box_intersect(&tin,&tout,x,y,z,vxp,vyp,vzp,xwidth+2*borderx,yheight+2*bordery,zdepth);
            if (tin>=0) {
              t1=tin;
              t2=tout;
            } else {
              t1=tout;
              t2=tin;
            }
            break;
          case 2:
            intersect=cylinder_intersect(&tin,&tout,x,y,z,vxp,vyp,vzp,radius,yheight+2*bordery);
            if (tin>=0) {
              t1=tin;
              t2=tout;
            } else {
              t1=tout;
              t2=tin;
            }
            break;
          case 4:
            /******************************************************************
            * What we have now is an absorption point somewhere between the two
            * cylinder walls, with distance 'zdepth' of the sensitive volume, and
            * a direction of emission of the reaction products.
            * What we need is very simple: for the first reaction product (the
            * proton for absorption in He-3) we need the smallest positive time,
            * for the second reaction product (triton) we choose the smallest
            * negative time.
            ******************************************************************/
            if (cylinder_intersect(&tinb,&toutb,x,y,z,vxp,vyp,vzp,rb,yheight+2*bordery)==0 &&
                !(VariousErrors&32) ) {
              fprintf(stderr,"PSD_Detector: %s: Oops. A neutron absorbed in the detector"
                             "              failed to send its reaction products towards the detector walls.\n",
                                            NAME_CURRENT_COMP);
            VariousErrors+=32;
            ABSORB;
            }
            if (cylinder_intersect(&tine,&toute,x,y,z,vxp,vyp,vzp,radius,yheight+2*bordery) ==0) {
              tine=-HUGE_VAL; // ugly fix for cases where the inner cylinder is missed
              toute=HUGE_VAL;
            }
            t1=tine;
            if (t1<0) {
              t1=toute;
              /* We are looking for a positive number. If t1 is negative, then
              choose the next number, toute, no matter what it is. */
            } else {
              if ((toute<t1)&&(toute>=0)) {
                t1=toute;
                /* If the next number, toute, is positive yet smaller than t1,
                we want it. */
              }
            }
            if (t1<0) {
              t1=tinb;
              /* We are looking for a positive number. If t1 is negative, then
              choose the next number, tinb, no matter what it is. */
            } else {
              if ((tinb<t1)&&(tinb>=0)) {
                t1=tinb;
                /* If the next number, tinb, is positive yet smaller than t1,
                we want it. */
              }
            }
            if (t1<0) {
              t1=toutb;
              /* We are looking for a positive number. If t1 is negative, then
              choose the next number, toutb, no matter what it is. */
            } else {
              if ((toutb<t1)&&(toutb>=0)) {
                t1=toutb;
                /* If the next number, toutb, is positive yet smaller than t1,
                we want it. */
              }
            }
            t2=tine;
            if (t2>0) {
              t2=toute;
              /* We are looking for a negative number. If t2 is positive, then
              choose the next number, toute, no matter what it is. */
            }
            else {
              if ((toute>t2)&&(toute<=0)) {
                t2=toute;
                /* If the next number, toute, is negative yet larger than t2,
                we want it. */
              }
            }
            if (t2>0) {
              t2=tinb;
              /* We are looking for a negative number. If t2 is positive, then
              choose the next number, tinb, no matter what it is. */
            } else {
              if ((tinb>t2)&&(tinb<=0)) {
                t2=tinb;
                /* If the next number, tinb, is negative yet larger than t2,
                we want it. */
              }
            }
            if (t2>0) {
              t2=toutb;
              /* We are looking for a negative number. If t2 is positive, then
              choose the next number, toutb, no matter what it is. */
            } else {
              if ((toutb>t2)&&(toutb<=0)) {
                t2=toutb;
                /* If the next number, toutb, is negative yet larger than t2,
                we want it. */
              }
            }
            break;
          default:
            fprintf(stderr,"PSD_Detector: %s: Detector has conflicting size\n"
                           "              specifications, i.e. combinations of xwidth, radius\n"
                           "              or awidth, or none at all. Exit.\n",
                                          NAME_CURRENT_COMP);
            exit(0);
        } // end switch detectortype (2)
        if (t1*t2>=0 && !(VariousErrors&64) ) {
          fprintf(stderr,"PSD_Detector: %s: t1 was %g and t2 was %g.\n"
                         "              One is supposed to be negative, the other positive,\n"
                         "              neither zero.\n",
                                        NAME_CURRENT_COMP,t1,t2);
          VariousErrors+=64;
        }
        if ( t1<RangeProton) lp =  t1; // t1 is the time of flight to the next wall
        else lp = RangeProton; // with normalized speed (1 m/s) it is also the distance
        if (-t2<RangeTriton) lt = -t2;
        else lt = RangeTriton;
        if (interpolate==1) {
          // take Bragg curves into account
          EnergyP = PSD_He_interp1value(PosAP,EAP,nH_p+1,lp);
          EnergyT = PSD_He_interp1value(PosAT,EAT,nH_t+1,lt);
          // first moment of the deposited energy distribution for the proton
          mlp = PSD_He_interp1value(PosAP,M1P1,nH_p+1,lp);
          // first moment of the deposited energy distribution for the triton
          mlt = PSD_He_interp1value(PosAT,M1T1,nH_t+1,lt);
        } else {
          // do not take Bragg curve into account, model energy deposition linearly
          EnergyP = FullEnergyP * lp / RangeProton;
          EnergyT = FullEnergyT * lt / RangeTriton;
          // first moment of the deposited energy distribution for the proton
          mlp = 0.5 * lp;
          // first moment of the deposited energy distribution for the triton
          mlt = 0.5 * lt;
        }
        Energy = EnergyP + EnergyT; // in keV, deposited in the detector
        // coordinates of the Center Of Gravity of the charged particle tracks
        xCOG = x + vxp * (mlp*EnergyP-mlt*EnergyT) / (Energy);
        yCOG = y + vyp * (mlp*EnergyP-mlt*EnergyT) / (Energy);
        zCOG = z + vzp * (mlp*EnergyP-mlt*EnergyT) / (Energy);
        if (DetectorType==4) {
          aCOG = atan2(xCOG,zCOG); // alpha coordinate (along the arc of the detector) of the Center Of Gravity of the charged particle tracks
          if (LensOn==1) {
            rCOG = sqrt(xCOG*xCOG+zCOG*zCOG);
            if (rCOG<radius+dc) {
              yCOG = yCOG * exp( (radius+zdepth-rCOG)*(radius+zdepth-rCOG) / (2*zdepth*radius) );
                /* correction of the Y center-of-gravity coordinate by a lens
                optimized according to P. Van Esch, NIM A 540 (2005) pp. 361-367   */
            }
          }
        }

/* STORAGE of detected neutrons ********************************************* */
        if (DetectorType==1 || DetectorType==4) {
          xCOG=xCOG + xwidth * xChDivRelSigma * randnorm();
          /* For a box or a banana, insert charge division relative sigma to
             obtain an extra error due to the charge division. This should be
             zero if independent channels are used along the x dimension. */
        }
        yCOG=yCOG + yheight * yChDivRelSigma * randnorm();
      /* For a box, tube or banana, insert charge division relative sigma to
         obtain an extra error due to the charge division. This should be
         zero if independent channels are used along the y dimension (which
         is never the case in a tube). */
        // Insert a reasonable-sounding energy resolution for the simulated detector.
        SigmaEnergy = 0.5*sqrt(Energy);
        Energy = Energy + SigmaEnergy * randnorm();

        // find the coordinates on the grid of pixels
        switch (DetectorType) {
          case 1:
            i = (long)floor( (xCOG + 0.5*xwidth)   *nx   /   xwidth );
            break;
          case 2:
            i = (long)floor( (xCOG + radius)       *nx   / (2*radius)  );
            break;
          case 4:
            i = (long)floor( (aCOG + 0.5*awidth/radius)*nx   / (awidth/radius) );
            break;
          default:
            break;
        } // end switch detectortype (3)
        j = (long)floor( (yCOG + 0.5*yheight)*ny / yheight );
        // events beyond the sensitive volume tend to induce signals on the border pixels of the detector
        if (i == -1) i = 0; if (i == nx) i = nx-1;
        if (j == -1) j = 0; if (j == ny) j = ny-1;
        if (i >= 0 && i < nx && j >= 0 && j < ny ) {
           SensVolCumul += p*pa; /* probability that a signal was produced in the
          sensitive volume, though not necessarily above the threshold.
          Note that, AT THIS POINT, the p has not been multiplied by pa yet. */
          if (Energy>threshold) {
            // if the particles have deposited sufficient energy to be
            // recognised as a neutron event (and not a gamma event)
            double P = p*pa;
            x = xCOG; y = yCOG; z = zCOG;
            SCATTER; // show point of detection for 3D view
            PSD_N[i][j] ++;       // one neutron tallied in the appropriate pixel
            PSD_p[i][j] += P;     // incremented by neutron weight
            PSD_p2[i][j] += P*P;  // 2nd order moments

            DetCumul += P; /* probability that a signal was produced in the
              sensitive volume above the energy threshold */
            i = (long)floor(Energy);
            if (i >= 0 &&  i < PHSpectrum_n ) {
              PHSpectrum0[i] ++;     // one neutron tallied per pixel
              PHSpectrum[i]  += P;   // energy bin in pulse height spectrum incremented by neutron weight
              PHSpectrum2[i] += P*P; // 2nd order moments
            }
            /* save event file if activated */
            if (type && strlen(type) && strcmp(type, "NULL") && strcmp(type, "0")) {
              double pp;
              Vars.cp  = P;
              Vars.cx  = x;
              Vars.cvx = vx;
              Vars.csx = sx;
              Vars.cy  = y;
              Vars.cvy = vy;
              Vars.csy = sy;
              Vars.cz  = z;
              Vars.cvz = vz;
              Vars.csz = sz;
              Vars.ct  = t;

              pp = Monitor_nD_Trace(&DEFS, &Vars);

              Vars.Nsum++;
              Vars.psum  += P;
              Vars.p2sum += P*P;
            }
            // initial version ABSORB detected neutrons.
            // This was removed if user wants to analyze the behaviour of the detector
          }
        } /* storage */
        p *= 1-pa;
      } // end if absorption in case of border effects (left and right) in the banana
    } /* end if p_interact */
  } /* end if (intersect) */
  if (restore_neutron) {
    RESTORE_NEUTRON(INDEX_CURRENT_COMP, x, y, z, vx, vy, vz, t, sx, sy, sz, p);
  }

  #undef nx
  #undef ny
  #undef xwidth
  #undef radius
  #undef awidth
  #undef yheight
  #undef zdepth
  #undef threshold
  #undef PressureConv
  #undef PressureStop
  #undef interpolate
  #undef p_interact
  #undef verbose
  #undef LensOn
  #undef dc
  #undef borderx
  #undef bordery
  #undef xChDivRelSigma
  #undef yChDivRelSigma
  #undef bufsize
  #undef restore_neutron
  #undef angle
  #undef type
  #undef filename
  #undef FN_Conv
  #undef FN_Stop
  #undef PSD_N
  #undef PSD_p
  #undef PSD_p2
  #undef EAP
  #undef M1P1
  #undef PosAP
  #undef EAT
  #undef M1T1
  #undef PosAT
  #undef CrossSectionHe
  #undef CountNeutrons
  #undef GeomCumul
  #undef AbsCumul
  #undef SensVolCumul
  #undef DetCumul
  #undef PHSpectrum
  #undef PHSpectrum0
  #undef PHSpectrum2
  #undef PHSpectrum_n
  #undef nH_p
  #undef nH_t
  #undef FullEnergyP
  #undef FullEnergyT
  #undef VariousErrors
  #undef DetectorType
  #undef DEFS
  #undef Vars
  return(_comp);
} /* class_PSD_Detector_trace */

/* *****************************************************************************
* instrument 'HZB_NEAT' TRACE
***************************************************************************** */

#pragma acc routine seq
int raytrace(_class_particle* _particle) { /* called by mccode_main for HZB_NEAT:TRACE */

  /* init variables and counters for TRACE */
  #undef ABSORB0
  #undef ABSORB
  #define ABSORB0 do { DEBUG_STATE(); DEBUG_ABSORB(); MAGNET_OFF; ABSORBED++; return(ABSORBED);} while(0)
  #define ABSORB ABSORB0
  DEBUG_ENTER();
  DEBUG_STATE();
  /* the main iteration loop for one incoming neutron */
  while (!ABSORBED) { /* iterate neutron event until absorbed */
    _class_particle _particle_save;
    /* send neutron event to component instance, one after the other */
    char flag_nocoordschange=0;
    if (!ABSORBED && _particle->_index == 1) {
      /* component PG=Progress_bar() [1] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_PG->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _PG->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_PG->_position_relative, _PG->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Progress_bar_trace(_PG, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component PG [1] */
    if (!ABSORBED && _particle->_index == 2) {
      /* component Source=Source_gen() [2] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Source->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Source->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Source->_position_relative, _Source->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Source_gen_trace(_Source, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component Source [2] */
    if (!ABSORBED && _particle->_index == 3) {
      /* component Guide0=Guide_gravity() [3] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Guide0->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Guide0->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Guide0->_position_relative, _Guide0->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Guide_gravity_trace(_Guide0, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component Guide0 [3] */
    if (!ABSORBED && _particle->_index == 4) {
      /* component Guide0_4=Guide_gravity() [4] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Guide0_4->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Guide0_4->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Guide0_4->_position_relative, _Guide0_4->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Guide_gravity_trace(_Guide0_4, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component Guide0_4 [4] */
    if (!ABSORBED && _particle->_index == 5) {
      /* component Guide0_5=Guide_gravity() [5] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Guide0_5->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Guide0_5->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Guide0_5->_position_relative, _Guide0_5->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Guide_gravity_trace(_Guide0_5, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component Guide0_5 [5] */
    if (!ABSORBED && _particle->_index == 6) {
      /* component Guide0_6=Guide_gravity() [6] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Guide0_6->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Guide0_6->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Guide0_6->_position_relative, _Guide0_6->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Guide_gravity_trace(_Guide0_6, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component Guide0_6 [6] */
    if (!ABSORBED && _particle->_index == 7) {
      /* component Guide0_7=Guide_gravity() [7] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Guide0_7->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Guide0_7->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Guide0_7->_position_relative, _Guide0_7->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Guide_gravity_trace(_Guide0_7, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component Guide0_7 [7] */
    if (!ABSORBED && _particle->_index == 8) {
      /* component Guide0_8=Guide_gravity() [8] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Guide0_8->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Guide0_8->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Guide0_8->_position_relative, _Guide0_8->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Guide_gravity_trace(_Guide0_8, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component Guide0_8 [8] */
    if (!ABSORBED && _particle->_index == 9) {
      /* component Guide0_9=Guide_gravity() [9] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Guide0_9->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Guide0_9->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Guide0_9->_position_relative, _Guide0_9->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Guide_gravity_trace(_Guide0_9, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component Guide0_9 [9] */
    if (!ABSORBED && _particle->_index == 10) {
      /* component Guide0_10=Guide_gravity() [10] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Guide0_10->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Guide0_10->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Guide0_10->_position_relative, _Guide0_10->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Guide_gravity_trace(_Guide0_10, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component Guide0_10 [10] */
    if (!ABSORBED && _particle->_index == 11) {
      /* component Guide0_11=Guide_gravity() [11] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Guide0_11->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Guide0_11->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Guide0_11->_position_relative, _Guide0_11->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Guide_gravity_trace(_Guide0_11, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component Guide0_11 [11] */
    if (!ABSORBED && _particle->_index == 12) {
      /* component Guide0_12=Guide_gravity() [12] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Guide0_12->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Guide0_12->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Guide0_12->_position_relative, _Guide0_12->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Guide_gravity_trace(_Guide0_12, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component Guide0_12 [12] */
    if (!ABSORBED && _particle->_index == 13) {
      /* component Guide0_13=Guide_gravity() [13] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Guide0_13->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Guide0_13->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Guide0_13->_position_relative, _Guide0_13->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Guide_gravity_trace(_Guide0_13, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component Guide0_13 [13] */
    if (!ABSORBED && _particle->_index == 14) {
      /* component Guide0_14=Guide_gravity() [14] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Guide0_14->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Guide0_14->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Guide0_14->_position_relative, _Guide0_14->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Guide_gravity_trace(_Guide0_14, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component Guide0_14 [14] */
    if (!ABSORBED && _particle->_index == 15) {
      /* component Guide0_15=Guide_gravity() [15] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Guide0_15->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Guide0_15->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Guide0_15->_position_relative, _Guide0_15->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Guide_gravity_trace(_Guide0_15, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component Guide0_15 [15] */
    if (!ABSORBED && _particle->_index == 16) {
      /* component Guide0_16=Guide_gravity() [16] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Guide0_16->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Guide0_16->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Guide0_16->_position_relative, _Guide0_16->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Guide_gravity_trace(_Guide0_16, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component Guide0_16 [16] */
    if (!ABSORBED && _particle->_index == 17) {
      /* component Guide0_17=Guide_gravity() [17] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Guide0_17->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Guide0_17->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Guide0_17->_position_relative, _Guide0_17->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Guide_gravity_trace(_Guide0_17, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component Guide0_17 [17] */
    if (!ABSORBED && _particle->_index == 18) {
      /* component Guide0_18=Guide_gravity() [18] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Guide0_18->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Guide0_18->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Guide0_18->_position_relative, _Guide0_18->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Guide_gravity_trace(_Guide0_18, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component Guide0_18 [18] */
    if (!ABSORBED && _particle->_index == 19) {
      /* component Guide0_19=Guide_gravity() [19] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Guide0_19->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Guide0_19->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Guide0_19->_position_relative, _Guide0_19->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Guide_gravity_trace(_Guide0_19, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component Guide0_19 [19] */
    if (!ABSORBED && _particle->_index == 20) {
      /* component Guide_PSD=Monitor_nD() [20] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Guide_PSD->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Guide_PSD->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Guide_PSD->_position_relative, _Guide_PSD->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Monitor_nD_trace(_Guide_PSD, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component Guide_PSD [20] */
    if (!ABSORBED && _particle->_index == 21) {
      /* component Guide_Lambda=Monitor_nD() [21] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Guide_Lambda->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Guide_Lambda->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Guide_Lambda->_position_relative, _Guide_Lambda->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Monitor_nD_trace(_Guide_Lambda, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component Guide_Lambda [21] */
    if (!ABSORBED && _particle->_index == 22) {
      /* component Chopper1=DiskChopper() [22] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Chopper1->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Chopper1->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Chopper1->_position_relative, _Chopper1->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_DiskChopper_trace(_Chopper1, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component Chopper1 [22] */
    if (!ABSORBED && _particle->_index == 23) {
      /* component Guide_SGS1=Guide_gravity() [23] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Guide_SGS1->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Guide_SGS1->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Guide_SGS1->_position_relative, _Guide_SGS1->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Guide_gravity_trace(_Guide_SGS1, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component Guide_SGS1 [23] */
    if (!ABSORBED && _particle->_index == 24) {
      /* component Guide_SGS2=Guide_gravity() [24] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Guide_SGS2->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Guide_SGS2->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Guide_SGS2->_position_relative, _Guide_SGS2->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Guide_gravity_trace(_Guide_SGS2, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component Guide_SGS2 [24] */
    if (!ABSORBED && _particle->_index == 25) {
      /* component Guide_SGS3=Guide_gravity() [25] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Guide_SGS3->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Guide_SGS3->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Guide_SGS3->_position_relative, _Guide_SGS3->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Guide_gravity_trace(_Guide_SGS3, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component Guide_SGS3 [25] */
    if (!ABSORBED && _particle->_index == 26) {
      /* component Guide_CGS=Guide_gravity() [26] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Guide_CGS->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Guide_CGS->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Guide_CGS->_position_relative, _Guide_CGS->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Guide_gravity_trace(_Guide_CGS, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component Guide_CGS [26] */
    if (!ABSORBED && _particle->_index == 27) {
      /* component Guide2_PSD=Monitor_nD() [27] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Guide2_PSD->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Guide2_PSD->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Guide2_PSD->_position_relative, _Guide2_PSD->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Monitor_nD_trace(_Guide2_PSD, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component Guide2_PSD [27] */
    if (!ABSORBED && _particle->_index == 28) {
      /* component Guide2_Lambda=Monitor_nD() [28] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Guide2_Lambda->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Guide2_Lambda->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Guide2_Lambda->_position_relative, _Guide2_Lambda->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Monitor_nD_trace(_Guide2_Lambda, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component Guide2_Lambda [28] */
    if (!ABSORBED && _particle->_index == 29) {
      /* component Guide2_dXY=Monitor_nD() [29] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Guide2_dXY->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Guide2_dXY->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Guide2_dXY->_position_relative, _Guide2_dXY->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Monitor_nD_trace(_Guide2_dXY, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component Guide2_dXY [29] */
    if (!ABSORBED && _particle->_index == 30) {
      /* component Chopper6=DiskChopper() [30] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Chopper6->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Chopper6->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Chopper6->_position_relative, _Chopper6->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_DiskChopper_trace(_Chopper6, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component Chopper6 [30] */
    if (!ABSORBED && _particle->_index == 31) {
      /* component Guide_time=Monitor_nD() [31] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Guide_time->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Guide_time->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Guide_time->_position_relative, _Guide_time->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Monitor_nD_trace(_Guide_time, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component Guide_time [31] */
    if (!ABSORBED && _particle->_index == 32) {
      /* component Guide_DGS=Guide_gravity() [32] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Guide_DGS->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Guide_DGS->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Guide_DGS->_position_relative, _Guide_DGS->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Guide_gravity_trace(_Guide_DGS, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component Guide_DGS [32] */
    if (!ABSORBED && _particle->_index == 33) {
      /* component Sample_pos=Arm() [33] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Sample_pos->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Sample_pos->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Sample_pos->_position_relative, _Sample_pos->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      _particle->_index++;
    } /* end component Sample_pos [33] */
    /* start SPLIT at Sample */
    instrument->logic.Split_Sample_particle=_particle;
    if (!ABSORBED && _particle->_index == 34) 
    for (instrument->logic.Split_Sample = 0; instrument->logic.Split_Sample< 10; instrument->logic.Split_Sample++) {
    {
      /* component Sample=Isotropic_Sqw() [34] */
      _particle=instrument->logic.Split_Sample_particle;
      p /= 10 > 0 ? 10 : 1;
  #define powder_format (_Sample->_parameters.powder_format)
  #define Sqw_coh (_Sample->_parameters.Sqw_coh)
  #define Sqw_inc (_Sample->_parameters.Sqw_inc)
  #define geometry (_Sample->_parameters.geometry)
  #define radius (_Sample->_parameters.radius)
  #define thickness (_Sample->_parameters.thickness)
  #define xwidth (_Sample->_parameters.xwidth)
  #define yheight (_Sample->_parameters.yheight)
  #define zdepth (_Sample->_parameters.zdepth)
  #define threshold (_Sample->_parameters.threshold)
  #define order (_Sample->_parameters.order)
  #define T (_Sample->_parameters.T)
  #define verbose (_Sample->_parameters.verbose)
  #define d_phi (_Sample->_parameters.d_phi)
  #define concentric (_Sample->_parameters.concentric)
  #define rho (_Sample->_parameters.rho)
  #define sigma_abs (_Sample->_parameters.sigma_abs)
  #define sigma_coh (_Sample->_parameters.sigma_coh)
  #define sigma_inc (_Sample->_parameters.sigma_inc)
  #define classical (_Sample->_parameters.classical)
  #define powder_Dd (_Sample->_parameters.powder_Dd)
  #define powder_DW (_Sample->_parameters.powder_DW)
  #define powder_Vc (_Sample->_parameters.powder_Vc)
  #define density (_Sample->_parameters.density)
  #define weight (_Sample->_parameters.weight)
  #define p_interact (_Sample->_parameters.p_interact)
  #define norm (_Sample->_parameters.norm)
  #define powder_barns (_Sample->_parameters.powder_barns)
  #define quantum_correction (_Sample->_parameters.quantum_correction)
  #define VarSqw (_Sample->_parameters.VarSqw)
  #define columns (_Sample->_parameters.columns)
  #define offdata (_Sample->_parameters.offdata)
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Sample->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Sample->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Sample->_position_relative, _Sample->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Isotropic_Sqw_trace(_Sample, _particle); /* contains EXTEND code */
      if (_particle->_restore)
        *_particle = _particle_save;
  #undef powder_format
  #undef Sqw_coh
  #undef Sqw_inc
  #undef geometry
  #undef radius
  #undef thickness
  #undef xwidth
  #undef yheight
  #undef zdepth
  #undef threshold
  #undef order
  #undef T
  #undef verbose
  #undef d_phi
  #undef concentric
  #undef rho
  #undef sigma_abs
  #undef sigma_coh
  #undef sigma_inc
  #undef classical
  #undef powder_Dd
  #undef powder_DW
  #undef powder_Vc
  #undef density
  #undef weight
  #undef p_interact
  #undef norm
  #undef powder_barns
  #undef quantum_correction
  #undef VarSqw
  #undef columns
  #undef offdata
      _particle->_index++;
    } /* end component Sample [34] */
    if (!ABSORBED && _particle->_index == 35) {
      /* component Ideal_Det=Monitor_nD() [35] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Ideal_Det->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Ideal_Det->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Ideal_Det->_position_relative, _Ideal_Det->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Monitor_nD_trace(_Ideal_Det, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component Ideal_Det [35] */
    if (!ABSORBED && _particle->_index == 36) {
      /* component Detector=PSD_Detector() [36] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Detector->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Detector->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Detector->_position_relative, _Detector->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_PSD_Detector_trace(_Detector, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component Detector [36] */
    } /* end SPLIT at Sample */
    if (_particle->_index > 36)
      ABSORBED++; /* absorbed when passed all components */
  } /* while !ABSORBED */

  DEBUG_LEAVE()
  DEBUG_STATE()

  return(_particle->_index);
} /* raytrace */
#undef x
#undef y
#undef z
#undef vx
#undef vy
#undef vz
#undef t
#undef sx
#undef sy
#undef sz
#undef p
// user variables:

#undef SCATTERED
#undef RESTORE
#undef RESTORE_NEUTRON
#undef STORE_NEUTRON
#undef ABSORBED
#undef ABSORB
#undef ABSORB0
/* *****************************************************************************
* instrument 'HZB_NEAT' and components SAVE
***************************************************************************** */

_class_Progress_bar *class_Progress_bar_save(_class_Progress_bar *_comp
) {
  #define profile (_comp->_parameters.profile)
  #define percent (_comp->_parameters.percent)
  #define flag_save (_comp->_parameters.flag_save)
  #define minutes (_comp->_parameters.minutes)
  #define IntermediateCnts (_comp->_parameters.IntermediateCnts)
  #define StartTime (_comp->_parameters.StartTime)
  #define EndTime (_comp->_parameters.EndTime)
  #define CurrentTime (_comp->_parameters.CurrentTime)
  MPI_MASTER(fprintf(stdout, "\nSave [%s]\n", instrument_name););
  if (profile && strlen(profile) && strcmp(profile,"NULL") && strcmp(profile,"0")) {
    char filename[256];
    if (!strlen(profile) || !strcmp(profile,"NULL") || !strcmp(profile,"0")) strcpy(filename, instrument_name);
    else strcpy(filename, profile);
    DETECTOR_OUT_1D(
        "Intensity profiler",
        "Component index [1]",
        "Intensity",
        "prof", 1, mcNUMCOMP, mcNUMCOMP-1,
        &(instrument->counter_N[1]),&(instrument->counter_P[1]),&(instrument->counter_P2[1]),
        filename);

  }
  #undef profile
  #undef percent
  #undef flag_save
  #undef minutes
  #undef IntermediateCnts
  #undef StartTime
  #undef EndTime
  #undef CurrentTime
  return(_comp);
} /* class_Progress_bar_save */

_class_Monitor_nD *class_Monitor_nD_save(_class_Monitor_nD *_comp
) {
  #define user1 (_comp->_parameters.user1)
  #define user2 (_comp->_parameters.user2)
  #define user3 (_comp->_parameters.user3)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define zdepth (_comp->_parameters.zdepth)
  #define xmin (_comp->_parameters.xmin)
  #define xmax (_comp->_parameters.xmax)
  #define ymin (_comp->_parameters.ymin)
  #define ymax (_comp->_parameters.ymax)
  #define zmin (_comp->_parameters.zmin)
  #define zmax (_comp->_parameters.zmax)
  #define bins (_comp->_parameters.bins)
  #define min (_comp->_parameters.min)
  #define max (_comp->_parameters.max)
  #define restore_neutron (_comp->_parameters.restore_neutron)
  #define radius (_comp->_parameters.radius)
  #define options (_comp->_parameters.options)
  #define filename (_comp->_parameters.filename)
  #define geometry (_comp->_parameters.geometry)
  #define username1 (_comp->_parameters.username1)
  #define username2 (_comp->_parameters.username2)
  #define username3 (_comp->_parameters.username3)
  #define DEFS (_comp->_parameters.DEFS)
  #define Vars (_comp->_parameters.Vars)
  #define detector (_comp->_parameters.detector)
  #define offdata (_comp->_parameters.offdata)
  /* save results, but do not free pointers */
  detector = Monitor_nD_Save(&DEFS, &Vars);
  #undef user1
  #undef user2
  #undef user3
  #undef xwidth
  #undef yheight
  #undef zdepth
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef zmin
  #undef zmax
  #undef bins
  #undef min
  #undef max
  #undef restore_neutron
  #undef radius
  #undef options
  #undef filename
  #undef geometry
  #undef username1
  #undef username2
  #undef username3
  #undef DEFS
  #undef Vars
  #undef detector
  #undef offdata
  return(_comp);
} /* class_Monitor_nD_save */

_class_PSD_Detector *class_PSD_Detector_save(_class_PSD_Detector *_comp
) {
  #define nx (_comp->_parameters.nx)
  #define ny (_comp->_parameters.ny)
  #define xwidth (_comp->_parameters.xwidth)
  #define radius (_comp->_parameters.radius)
  #define awidth (_comp->_parameters.awidth)
  #define yheight (_comp->_parameters.yheight)
  #define zdepth (_comp->_parameters.zdepth)
  #define threshold (_comp->_parameters.threshold)
  #define PressureConv (_comp->_parameters.PressureConv)
  #define PressureStop (_comp->_parameters.PressureStop)
  #define interpolate (_comp->_parameters.interpolate)
  #define p_interact (_comp->_parameters.p_interact)
  #define verbose (_comp->_parameters.verbose)
  #define LensOn (_comp->_parameters.LensOn)
  #define dc (_comp->_parameters.dc)
  #define borderx (_comp->_parameters.borderx)
  #define bordery (_comp->_parameters.bordery)
  #define xChDivRelSigma (_comp->_parameters.xChDivRelSigma)
  #define yChDivRelSigma (_comp->_parameters.yChDivRelSigma)
  #define bufsize (_comp->_parameters.bufsize)
  #define restore_neutron (_comp->_parameters.restore_neutron)
  #define angle (_comp->_parameters.angle)
  #define type (_comp->_parameters.type)
  #define filename (_comp->_parameters.filename)
  #define FN_Conv (_comp->_parameters.FN_Conv)
  #define FN_Stop (_comp->_parameters.FN_Stop)
  #define PSD_N (_comp->_parameters.PSD_N)
  #define PSD_p (_comp->_parameters.PSD_p)
  #define PSD_p2 (_comp->_parameters.PSD_p2)
  #define EAP (_comp->_parameters.EAP)
  #define M1P1 (_comp->_parameters.M1P1)
  #define PosAP (_comp->_parameters.PosAP)
  #define EAT (_comp->_parameters.EAT)
  #define M1T1 (_comp->_parameters.M1T1)
  #define PosAT (_comp->_parameters.PosAT)
  #define CrossSectionHe (_comp->_parameters.CrossSectionHe)
  #define CountNeutrons (_comp->_parameters.CountNeutrons)
  #define GeomCumul (_comp->_parameters.GeomCumul)
  #define AbsCumul (_comp->_parameters.AbsCumul)
  #define SensVolCumul (_comp->_parameters.SensVolCumul)
  #define DetCumul (_comp->_parameters.DetCumul)
  #define PHSpectrum (_comp->_parameters.PHSpectrum)
  #define PHSpectrum0 (_comp->_parameters.PHSpectrum0)
  #define PHSpectrum2 (_comp->_parameters.PHSpectrum2)
  #define PHSpectrum_n (_comp->_parameters.PHSpectrum_n)
  #define nH_p (_comp->_parameters.nH_p)
  #define nH_t (_comp->_parameters.nH_t)
  #define FullEnergyP (_comp->_parameters.FullEnergyP)
  #define FullEnergyT (_comp->_parameters.FullEnergyT)
  #define VariousErrors (_comp->_parameters.VariousErrors)
  #define DetectorType (_comp->_parameters.DetectorType)
  #define DEFS (_comp->_parameters.DEFS)
  #define Vars (_comp->_parameters.Vars)
  char   file[128];
  if (type && strlen(type) && strcmp(type, "NULL") && strcmp(type, "0")) {
    /* event file */
    Monitor_nD_Save(&DEFS, &Vars);
    DETECTOR_OUT(Vars.Nsum, Vars.psum, Vars.p2sum);
  } else {
    double width;

    if (xwidth>0) width=xwidth;
    if (radius>0) width=2*radius;
    if (awidth>0) width=awidth;
    if (filename && strlen(filename) && strcmp(filename, "0") && strcmp(filename, "NULL"))
      strncpy(file, filename, 128);
    else
      sprintf(file, "%s.dat", NAME_CURRENT_COMP);
    if (nx > 1 && ny > 1)
    DETECTOR_OUT_2D(
          "PSD Detector",
          "X position [m]",
          "Y position [m]",
          -width/2, width/2, -yheight/2, yheight/2,
          (double)nx, (double)ny,
          &PSD_N[0][0],&PSD_p[0][0],&PSD_p2[0][0],
          file);
    else if (nx == 1)
    DETECTOR_OUT_1D(
          "PSD Detector",
          "Y position [m]","Counts","Y",
          -yheight/2, yheight/2, (double)ny,
          &PSD_N[0][0],&PSD_p[0][0],&PSD_p2[0][0],
          file);
    else if (ny == 1)
    DETECTOR_OUT_1D(
          "PSD Detector",
          "X position [m]","Counts","X",
          -width/2, width/2, (double)nx,
          &PSD_N[0][0],&PSD_p[0][0],&PSD_p2[0][0],
          file);
  }
  if (verbose) {
    char xlabelstr[64];
    snprintf(file, 128, "%s.en",
      filename && strlen(filename) && strcmp(filename,"0") && strcmp(filename,"NULL") ?
            filename : NAME_CURRENT_COMP);
    fprintf(stdout,"PSD_Detector: %s: statistics\n", NAME_CURRENT_COMP);
    fprintf(stdout,"  %g neutrons in the simulation, %g encounter the detector.\n",
                   CountNeutrons, GeomCumul );
    fprintf(stdout,"  Probability for a neutron to be absorbed in the detector is %g percent.\n",
                   100*AbsCumul/CountNeutrons );
    fprintf(stdout,"  Probability for a neutron to be detected is %g percent.\n",
                   100*DetCumul/CountNeutrons );
    fprintf(stdout,"  Fraction of neutrons not counted because their COG is outside\n"
                   "the sensitive volume is %g.\n",
                   1-SensVolCumul/AbsCumul );
    fprintf(stdout,"  Fraction of neutrons not counted because their signal is reduced below\n"
                   "%g keV due to the wall effect is %g.\n",
                   threshold,1-DetCumul/SensVolCumul );
    fprintf(stdout,"  Theoretical limit to the position resolution in this gas is %g m FWHM\n"
                   "  (of the rectangular distribution). This implies a sigma of %g m.\n",
                   2*(M1P1[nH_p]-M1T1[nH_t]),(M1P1[nH_p]-M1T1[nH_t])/sqrt(3) );

    snprintf(xlabelstr,64,"Energy [keV], threshold set to %g keV",threshold);
    DETECTOR_OUT_1D(
      "Pulse Height Spectrum",
      xlabelstr,
      "Counts [a.u]",
      "E",
      0.0, (double)(PHSpectrum_n-1), PHSpectrum_n,
      PHSpectrum0, PHSpectrum, PHSpectrum2,
      file);
  }
  #undef nx
  #undef ny
  #undef xwidth
  #undef radius
  #undef awidth
  #undef yheight
  #undef zdepth
  #undef threshold
  #undef PressureConv
  #undef PressureStop
  #undef interpolate
  #undef p_interact
  #undef verbose
  #undef LensOn
  #undef dc
  #undef borderx
  #undef bordery
  #undef xChDivRelSigma
  #undef yChDivRelSigma
  #undef bufsize
  #undef restore_neutron
  #undef angle
  #undef type
  #undef filename
  #undef FN_Conv
  #undef FN_Stop
  #undef PSD_N
  #undef PSD_p
  #undef PSD_p2
  #undef EAP
  #undef M1P1
  #undef PosAP
  #undef EAT
  #undef M1T1
  #undef PosAT
  #undef CrossSectionHe
  #undef CountNeutrons
  #undef GeomCumul
  #undef AbsCumul
  #undef SensVolCumul
  #undef DetCumul
  #undef PHSpectrum
  #undef PHSpectrum0
  #undef PHSpectrum2
  #undef PHSpectrum_n
  #undef nH_p
  #undef nH_t
  #undef FullEnergyP
  #undef FullEnergyT
  #undef VariousErrors
  #undef DetectorType
  #undef DEFS
  #undef Vars
  return(_comp);
} /* class_PSD_Detector_save */



int save(FILE *handle) { /* called by mccode_main for HZB_NEAT:SAVE */
  if (!handle) siminfo_init(NULL);

  /* call iteratively all components SAVE */
  class_Progress_bar_save(_PG);



















  class_Monitor_nD_save(_Guide_PSD);

  class_Monitor_nD_save(_Guide_Lambda);






  class_Monitor_nD_save(_Guide2_PSD);

  class_Monitor_nD_save(_Guide2_Lambda);

  class_Monitor_nD_save(_Guide2_dXY);


  class_Monitor_nD_save(_Guide_time);




  class_Monitor_nD_save(_Ideal_Det);

  class_PSD_Detector_save(_Detector);

  if (!handle) siminfo_close(); 

  return(0);
} /* save */

/* *****************************************************************************
* instrument 'HZB_NEAT' and components FINALLY
***************************************************************************** */

_class_Progress_bar *class_Progress_bar_finally(_class_Progress_bar *_comp
) {
  #define profile (_comp->_parameters.profile)
  #define percent (_comp->_parameters.percent)
  #define flag_save (_comp->_parameters.flag_save)
  #define minutes (_comp->_parameters.minutes)
  #define IntermediateCnts (_comp->_parameters.IntermediateCnts)
  #define StartTime (_comp->_parameters.StartTime)
  #define EndTime (_comp->_parameters.EndTime)
  #define CurrentTime (_comp->_parameters.CurrentTime)
  time_t NowTime;
  time(&NowTime);
  fprintf(stdout, "\nFinally [%s: %s]. Time: ", instrument_name, dirname ? dirname : ".");
  if (difftime(NowTime,StartTime) < 60.0)
    fprintf(stdout, "%g [s] ", difftime(NowTime,StartTime));
  else if (difftime(NowTime,StartTime) > 3600.0)
    fprintf(stdout, "%g [h] ", difftime(NowTime,StartTime)/3660.0);
  else
    fprintf(stdout, "%g [min] ", difftime(NowTime,StartTime)/60.0);
  fprintf(stdout, "\n");
  #undef profile
  #undef percent
  #undef flag_save
  #undef minutes
  #undef IntermediateCnts
  #undef StartTime
  #undef EndTime
  #undef CurrentTime
  return(_comp);
} /* class_Progress_bar_finally */

_class_Source_gen *class_Source_gen_finally(_class_Source_gen *_comp
) {
  #define flux_file (_comp->_parameters.flux_file)
  #define xdiv_file (_comp->_parameters.xdiv_file)
  #define ydiv_file (_comp->_parameters.ydiv_file)
  #define radius (_comp->_parameters.radius)
  #define dist (_comp->_parameters.dist)
  #define focus_xw (_comp->_parameters.focus_xw)
  #define focus_yh (_comp->_parameters.focus_yh)
  #define focus_aw (_comp->_parameters.focus_aw)
  #define focus_ah (_comp->_parameters.focus_ah)
  #define E0 (_comp->_parameters.E0)
  #define dE (_comp->_parameters.dE)
  #define lambda0 (_comp->_parameters.lambda0)
  #define dlambda (_comp->_parameters.dlambda)
  #define I1 (_comp->_parameters.I1)
  #define yheight (_comp->_parameters.yheight)
  #define xwidth (_comp->_parameters.xwidth)
  #define verbose (_comp->_parameters.verbose)
  #define T1 (_comp->_parameters.T1)
  #define flux_file_perAA (_comp->_parameters.flux_file_perAA)
  #define flux_file_log (_comp->_parameters.flux_file_log)
  #define Lmin (_comp->_parameters.Lmin)
  #define Lmax (_comp->_parameters.Lmax)
  #define Emin (_comp->_parameters.Emin)
  #define Emax (_comp->_parameters.Emax)
  #define T2 (_comp->_parameters.T2)
  #define I2 (_comp->_parameters.I2)
  #define T3 (_comp->_parameters.T3)
  #define I3 (_comp->_parameters.I3)
  #define zdepth (_comp->_parameters.zdepth)
  #define target_index (_comp->_parameters.target_index)
  #define p_in (_comp->_parameters.p_in)
  #define lambda1 (_comp->_parameters.lambda1)
  #define lambda2 (_comp->_parameters.lambda2)
  #define lambda3 (_comp->_parameters.lambda3)
  #define pTable (_comp->_parameters.pTable)
  #define pTable_x (_comp->_parameters.pTable_x)
  #define pTable_y (_comp->_parameters.pTable_y)
  #define pTable_xmin (_comp->_parameters.pTable_xmin)
  #define pTable_xmax (_comp->_parameters.pTable_xmax)
  #define pTable_xsum (_comp->_parameters.pTable_xsum)
  #define pTable_ymin (_comp->_parameters.pTable_ymin)
  #define pTable_ymax (_comp->_parameters.pTable_ymax)
  #define pTable_ysum (_comp->_parameters.pTable_ysum)
  #define pTable_dxmin (_comp->_parameters.pTable_dxmin)
  #define pTable_dxmax (_comp->_parameters.pTable_dxmax)
  #define pTable_dymin (_comp->_parameters.pTable_dymin)
  #define pTable_dymax (_comp->_parameters.pTable_dymax)
  Table_Free(&pTable);
  Table_Free(&pTable_x);
  Table_Free(&pTable_y);
  #undef flux_file
  #undef xdiv_file
  #undef ydiv_file
  #undef radius
  #undef dist
  #undef focus_xw
  #undef focus_yh
  #undef focus_aw
  #undef focus_ah
  #undef E0
  #undef dE
  #undef lambda0
  #undef dlambda
  #undef I1
  #undef yheight
  #undef xwidth
  #undef verbose
  #undef T1
  #undef flux_file_perAA
  #undef flux_file_log
  #undef Lmin
  #undef Lmax
  #undef Emin
  #undef Emax
  #undef T2
  #undef I2
  #undef T3
  #undef I3
  #undef zdepth
  #undef target_index
  #undef p_in
  #undef lambda1
  #undef lambda2
  #undef lambda3
  #undef pTable
  #undef pTable_x
  #undef pTable_y
  #undef pTable_xmin
  #undef pTable_xmax
  #undef pTable_xsum
  #undef pTable_ymin
  #undef pTable_ymax
  #undef pTable_ysum
  #undef pTable_dxmin
  #undef pTable_dxmax
  #undef pTable_dymin
  #undef pTable_dymax
  return(_comp);
} /* class_Source_gen_finally */

_class_Guide_gravity *class_Guide_gravity_finally(_class_Guide_gravity *_comp
) {
  #define w1 (_comp->_parameters.w1)
  #define h1 (_comp->_parameters.h1)
  #define w2 (_comp->_parameters.w2)
  #define h2 (_comp->_parameters.h2)
  #define l (_comp->_parameters.l)
  #define R0 (_comp->_parameters.R0)
  #define Qc (_comp->_parameters.Qc)
  #define alpha (_comp->_parameters.alpha)
  #define m (_comp->_parameters.m)
  #define W (_comp->_parameters.W)
  #define nslit (_comp->_parameters.nslit)
  #define d (_comp->_parameters.d)
  #define mleft (_comp->_parameters.mleft)
  #define mright (_comp->_parameters.mright)
  #define mtop (_comp->_parameters.mtop)
  #define mbottom (_comp->_parameters.mbottom)
  #define nhslit (_comp->_parameters.nhslit)
  #define G (_comp->_parameters.G)
  #define aleft (_comp->_parameters.aleft)
  #define aright (_comp->_parameters.aright)
  #define atop (_comp->_parameters.atop)
  #define abottom (_comp->_parameters.abottom)
  #define wavy (_comp->_parameters.wavy)
  #define wavy_z (_comp->_parameters.wavy_z)
  #define wavy_tb (_comp->_parameters.wavy_tb)
  #define wavy_lr (_comp->_parameters.wavy_lr)
  #define chamfers (_comp->_parameters.chamfers)
  #define chamfers_z (_comp->_parameters.chamfers_z)
  #define chamfers_lr (_comp->_parameters.chamfers_lr)
  #define chamfers_tb (_comp->_parameters.chamfers_tb)
  #define nelements (_comp->_parameters.nelements)
  #define nu (_comp->_parameters.nu)
  #define phase (_comp->_parameters.phase)
  #define reflect (_comp->_parameters.reflect)
  #define GVars (_comp->_parameters.GVars)
  #define pTable (_comp->_parameters.pTable)
if (GVars.warnings > 100) {
  fprintf(stderr,"%s: warning: neutron has entered guide, but can not exit !\n", GVars.compcurname);
  fprintf(stderr,"%s: warning: This message has been repeated %g times\n", GVars.compcurname, GVars.warnings);
}
  #undef w1
  #undef h1
  #undef w2
  #undef h2
  #undef l
  #undef R0
  #undef Qc
  #undef alpha
  #undef m
  #undef W
  #undef nslit
  #undef d
  #undef mleft
  #undef mright
  #undef mtop
  #undef mbottom
  #undef nhslit
  #undef G
  #undef aleft
  #undef aright
  #undef atop
  #undef abottom
  #undef wavy
  #undef wavy_z
  #undef wavy_tb
  #undef wavy_lr
  #undef chamfers
  #undef chamfers_z
  #undef chamfers_lr
  #undef chamfers_tb
  #undef nelements
  #undef nu
  #undef phase
  #undef reflect
  #undef GVars
  #undef pTable
  return(_comp);
} /* class_Guide_gravity_finally */

_class_Monitor_nD *class_Monitor_nD_finally(_class_Monitor_nD *_comp
) {
  #define user1 (_comp->_parameters.user1)
  #define user2 (_comp->_parameters.user2)
  #define user3 (_comp->_parameters.user3)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define zdepth (_comp->_parameters.zdepth)
  #define xmin (_comp->_parameters.xmin)
  #define xmax (_comp->_parameters.xmax)
  #define ymin (_comp->_parameters.ymin)
  #define ymax (_comp->_parameters.ymax)
  #define zmin (_comp->_parameters.zmin)
  #define zmax (_comp->_parameters.zmax)
  #define bins (_comp->_parameters.bins)
  #define min (_comp->_parameters.min)
  #define max (_comp->_parameters.max)
  #define restore_neutron (_comp->_parameters.restore_neutron)
  #define radius (_comp->_parameters.radius)
  #define options (_comp->_parameters.options)
  #define filename (_comp->_parameters.filename)
  #define geometry (_comp->_parameters.geometry)
  #define username1 (_comp->_parameters.username1)
  #define username2 (_comp->_parameters.username2)
  #define username3 (_comp->_parameters.username3)
  #define DEFS (_comp->_parameters.DEFS)
  #define Vars (_comp->_parameters.Vars)
  #define detector (_comp->_parameters.detector)
  #define offdata (_comp->_parameters.offdata)
  /* free pointers */
  Monitor_nD_Finally(&DEFS, &Vars);
  #undef user1
  #undef user2
  #undef user3
  #undef xwidth
  #undef yheight
  #undef zdepth
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef zmin
  #undef zmax
  #undef bins
  #undef min
  #undef max
  #undef restore_neutron
  #undef radius
  #undef options
  #undef filename
  #undef geometry
  #undef username1
  #undef username2
  #undef username3
  #undef DEFS
  #undef Vars
  #undef detector
  #undef offdata
  return(_comp);
} /* class_Monitor_nD_finally */

_class_Isotropic_Sqw *class_Isotropic_Sqw_finally(_class_Isotropic_Sqw *_comp
) {
  #define powder_format (_comp->_parameters.powder_format)
  #define Sqw_coh (_comp->_parameters.Sqw_coh)
  #define Sqw_inc (_comp->_parameters.Sqw_inc)
  #define geometry (_comp->_parameters.geometry)
  #define radius (_comp->_parameters.radius)
  #define thickness (_comp->_parameters.thickness)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define zdepth (_comp->_parameters.zdepth)
  #define threshold (_comp->_parameters.threshold)
  #define order (_comp->_parameters.order)
  #define T (_comp->_parameters.T)
  #define verbose (_comp->_parameters.verbose)
  #define d_phi (_comp->_parameters.d_phi)
  #define concentric (_comp->_parameters.concentric)
  #define rho (_comp->_parameters.rho)
  #define sigma_abs (_comp->_parameters.sigma_abs)
  #define sigma_coh (_comp->_parameters.sigma_coh)
  #define sigma_inc (_comp->_parameters.sigma_inc)
  #define classical (_comp->_parameters.classical)
  #define powder_Dd (_comp->_parameters.powder_Dd)
  #define powder_DW (_comp->_parameters.powder_DW)
  #define powder_Vc (_comp->_parameters.powder_Vc)
  #define density (_comp->_parameters.density)
  #define weight (_comp->_parameters.weight)
  #define p_interact (_comp->_parameters.p_interact)
  #define norm (_comp->_parameters.norm)
  #define powder_barns (_comp->_parameters.powder_barns)
  #define quantum_correction (_comp->_parameters.quantum_correction)
  #define VarSqw (_comp->_parameters.VarSqw)
  #define columns (_comp->_parameters.columns)
  #define offdata (_comp->_parameters.offdata)
  int  k;

  if (VarSqw.s_coh > 0 || VarSqw.s_inc > 0)
  for (k=0; k < 2; k++) {
    struct Sqw_Data_struct Data_sqw;

    Data_sqw =  (k == 0 ? VarSqw.Data_coh : VarSqw.Data_inc);
    /* Data_sqw->Sqw has already been freed at end of INIT */
    Table_Free(&(Data_sqw.iqSq));

    if (Data_sqw.SW)           free(Data_sqw.SW);
    if (Data_sqw.SQW)          free(Data_sqw.SQW);
    if (Data_sqw.SW_lookup)    free(Data_sqw.SW_lookup);
    if (Data_sqw.QW_lookup)    free(Data_sqw.QW_lookup);
  } /* end for */

#ifdef USE_MPI
  if (mpi_node_count > 1) {
    double tmp;
    tmp = (double)VarSqw.neutron_removed; mc_MPI_Sum(&tmp, 1); VarSqw.neutron_removed=(long)tmp;
    tmp = (double)VarSqw.neutron_exit;    mc_MPI_Sum(&tmp, 1); VarSqw.neutron_exit=(long)tmp;
    tmp = (double)VarSqw.neutron_pmult;   mc_MPI_Sum(&tmp, 1); VarSqw.neutron_pmult=(long)tmp;
    mc_MPI_Sum(&VarSqw.mean_scatt, 1);
    mc_MPI_Sum(&VarSqw.psum_scatt, 1);
    mc_MPI_Sum(&VarSqw.mean_abs, 1);
    mc_MPI_Sum(&VarSqw.single_coh, 1);
    mc_MPI_Sum(&VarSqw.single_inc, 1);
    mc_MPI_Sum(&VarSqw.multi, 1);
  }
#endif
  MPI_MASTER(
  if (VarSqw.neutron_removed)
    printf("Isotropic_Sqw: %s: %li neutron events (out of %li) that should have\n"
           "               scattered were transmitted because scattering conditions\n"
           "WARNING        could not be satisfied after %i tries.\n",
          NAME_CURRENT_COMP, VarSqw.neutron_removed,
          VarSqw.neutron_exit+VarSqw.neutron_removed, VarSqw.maxloop);
  if (VarSqw.neutron_pmult)
    printf("Isotropic_Sqw: %s: %li neutron events (out of %li) reached\n"
           "WARNING        unrealistic weight. The S(q,w) norm might be too high.\n",
          NAME_CURRENT_COMP, VarSqw.neutron_pmult, VarSqw.neutron_exit);

  if (VarSqw.verbose_output >= 1 && VarSqw.psum_scatt > 0) {
    printf("Isotropic_Sqw: %s: Scattering fraction=%g of incoming intensity\n"
           "               Absorption fraction           =%g\n",
           NAME_CURRENT_COMP,
           VarSqw.mean_scatt/VarSqw.psum_scatt, VarSqw.mean_abs/VarSqw.psum_scatt);
    printf("               Single   scattering intensity =%g (coh=%g inc=%g)\n"
           "               Multiple scattering intensity =%g\n",
           VarSqw.single_coh+VarSqw.single_inc, VarSqw.single_coh, VarSqw.single_inc, VarSqw.multi);
    );
  }

/* end FINALLY */
  #undef powder_format
  #undef Sqw_coh
  #undef Sqw_inc
  #undef geometry
  #undef radius
  #undef thickness
  #undef xwidth
  #undef yheight
  #undef zdepth
  #undef threshold
  #undef order
  #undef T
  #undef verbose
  #undef d_phi
  #undef concentric
  #undef rho
  #undef sigma_abs
  #undef sigma_coh
  #undef sigma_inc
  #undef classical
  #undef powder_Dd
  #undef powder_DW
  #undef powder_Vc
  #undef density
  #undef weight
  #undef p_interact
  #undef norm
  #undef powder_barns
  #undef quantum_correction
  #undef VarSqw
  #undef columns
  #undef offdata
  return(_comp);
} /* class_Isotropic_Sqw_finally */

_class_PSD_Detector *class_PSD_Detector_finally(_class_PSD_Detector *_comp
) {
  #define nx (_comp->_parameters.nx)
  #define ny (_comp->_parameters.ny)
  #define xwidth (_comp->_parameters.xwidth)
  #define radius (_comp->_parameters.radius)
  #define awidth (_comp->_parameters.awidth)
  #define yheight (_comp->_parameters.yheight)
  #define zdepth (_comp->_parameters.zdepth)
  #define threshold (_comp->_parameters.threshold)
  #define PressureConv (_comp->_parameters.PressureConv)
  #define PressureStop (_comp->_parameters.PressureStop)
  #define interpolate (_comp->_parameters.interpolate)
  #define p_interact (_comp->_parameters.p_interact)
  #define verbose (_comp->_parameters.verbose)
  #define LensOn (_comp->_parameters.LensOn)
  #define dc (_comp->_parameters.dc)
  #define borderx (_comp->_parameters.borderx)
  #define bordery (_comp->_parameters.bordery)
  #define xChDivRelSigma (_comp->_parameters.xChDivRelSigma)
  #define yChDivRelSigma (_comp->_parameters.yChDivRelSigma)
  #define bufsize (_comp->_parameters.bufsize)
  #define restore_neutron (_comp->_parameters.restore_neutron)
  #define angle (_comp->_parameters.angle)
  #define type (_comp->_parameters.type)
  #define filename (_comp->_parameters.filename)
  #define FN_Conv (_comp->_parameters.FN_Conv)
  #define FN_Stop (_comp->_parameters.FN_Stop)
  #define PSD_N (_comp->_parameters.PSD_N)
  #define PSD_p (_comp->_parameters.PSD_p)
  #define PSD_p2 (_comp->_parameters.PSD_p2)
  #define EAP (_comp->_parameters.EAP)
  #define M1P1 (_comp->_parameters.M1P1)
  #define PosAP (_comp->_parameters.PosAP)
  #define EAT (_comp->_parameters.EAT)
  #define M1T1 (_comp->_parameters.M1T1)
  #define PosAT (_comp->_parameters.PosAT)
  #define CrossSectionHe (_comp->_parameters.CrossSectionHe)
  #define CountNeutrons (_comp->_parameters.CountNeutrons)
  #define GeomCumul (_comp->_parameters.GeomCumul)
  #define AbsCumul (_comp->_parameters.AbsCumul)
  #define SensVolCumul (_comp->_parameters.SensVolCumul)
  #define DetCumul (_comp->_parameters.DetCumul)
  #define PHSpectrum (_comp->_parameters.PHSpectrum)
  #define PHSpectrum0 (_comp->_parameters.PHSpectrum0)
  #define PHSpectrum2 (_comp->_parameters.PHSpectrum2)
  #define PHSpectrum_n (_comp->_parameters.PHSpectrum_n)
  #define nH_p (_comp->_parameters.nH_p)
  #define nH_t (_comp->_parameters.nH_t)
  #define FullEnergyP (_comp->_parameters.FullEnergyP)
  #define FullEnergyT (_comp->_parameters.FullEnergyT)
  #define VariousErrors (_comp->_parameters.VariousErrors)
  #define DetectorType (_comp->_parameters.DetectorType)
  #define DEFS (_comp->_parameters.DEFS)
  #define Vars (_comp->_parameters.Vars)
  /* free pointers */
  if (type && strlen(type) && strcmp(type, "NULL") && strcmp(type, "0")) {
    Monitor_nD_Finally(&DEFS, &Vars);
    if (bufsize) {
      printf("PSD_Detector: %s: Saved %Li events (from buffer) in file %s\n",
             NAME_CURRENT_COMP, Vars.Nsum, Vars.Mon_File);
      if (bufsize < Vars.Nsum)
        printf("WARNING         When using this source, intensities must be multiplied\n"
               "                by a factor %g\n", (double)Vars.Nsum/bufsize);
    } else
      printf("PSD_Detector: %s: Saved %Li events (all) in file %s\n", NAME_CURRENT_COMP, Vars.Nsum, Vars.Mon_File);
  }

  destroy_darr2d(PSD_N);
  destroy_darr2d(PSD_p);
  destroy_darr2d(PSD_p2);
  #undef nx
  #undef ny
  #undef xwidth
  #undef radius
  #undef awidth
  #undef yheight
  #undef zdepth
  #undef threshold
  #undef PressureConv
  #undef PressureStop
  #undef interpolate
  #undef p_interact
  #undef verbose
  #undef LensOn
  #undef dc
  #undef borderx
  #undef bordery
  #undef xChDivRelSigma
  #undef yChDivRelSigma
  #undef bufsize
  #undef restore_neutron
  #undef angle
  #undef type
  #undef filename
  #undef FN_Conv
  #undef FN_Stop
  #undef PSD_N
  #undef PSD_p
  #undef PSD_p2
  #undef EAP
  #undef M1P1
  #undef PosAP
  #undef EAT
  #undef M1T1
  #undef PosAT
  #undef CrossSectionHe
  #undef CountNeutrons
  #undef GeomCumul
  #undef AbsCumul
  #undef SensVolCumul
  #undef DetCumul
  #undef PHSpectrum
  #undef PHSpectrum0
  #undef PHSpectrum2
  #undef PHSpectrum_n
  #undef nH_p
  #undef nH_t
  #undef FullEnergyP
  #undef FullEnergyT
  #undef VariousErrors
  #undef DetectorType
  #undef DEFS
  #undef Vars
  return(_comp);
} /* class_PSD_Detector_finally */



int finally(void) { /* called by mccode_main for HZB_NEAT:FINALLY */
  siminfo_init(NULL);
  save(siminfo_file); /* save data when simulation ends */

  /* call iteratively all components FINALLY */
  class_Progress_bar_finally(_PG);

  class_Source_gen_finally(_Source);

  class_Guide_gravity_finally(_Guide0);

  class_Guide_gravity_finally(_Guide0_4);

  class_Guide_gravity_finally(_Guide0_5);

  class_Guide_gravity_finally(_Guide0_6);

  class_Guide_gravity_finally(_Guide0_7);

  class_Guide_gravity_finally(_Guide0_8);

  class_Guide_gravity_finally(_Guide0_9);

  class_Guide_gravity_finally(_Guide0_10);

  class_Guide_gravity_finally(_Guide0_11);

  class_Guide_gravity_finally(_Guide0_12);

  class_Guide_gravity_finally(_Guide0_13);

  class_Guide_gravity_finally(_Guide0_14);

  class_Guide_gravity_finally(_Guide0_15);

  class_Guide_gravity_finally(_Guide0_16);

  class_Guide_gravity_finally(_Guide0_17);

  class_Guide_gravity_finally(_Guide0_18);

  class_Guide_gravity_finally(_Guide0_19);

  class_Monitor_nD_finally(_Guide_PSD);

  class_Monitor_nD_finally(_Guide_Lambda);


  class_Guide_gravity_finally(_Guide_SGS1);

  class_Guide_gravity_finally(_Guide_SGS2);

  class_Guide_gravity_finally(_Guide_SGS3);

  class_Guide_gravity_finally(_Guide_CGS);

  class_Monitor_nD_finally(_Guide2_PSD);

  class_Monitor_nD_finally(_Guide2_Lambda);

  class_Monitor_nD_finally(_Guide2_dXY);


  class_Monitor_nD_finally(_Guide_time);

  class_Guide_gravity_finally(_Guide_DGS);


  class_Isotropic_Sqw_finally(_Sample);

  class_Monitor_nD_finally(_Ideal_Det);

  class_PSD_Detector_finally(_Detector);

  siminfo_close(); 

  return(0);
} /* finally */

/* *****************************************************************************
* instrument 'HZB_NEAT' and components DISPLAY
***************************************************************************** */

  #define magnify     mcdis_magnify
  #define line        mcdis_line
  #define dashed_line mcdis_dashed_line
  #define multiline   mcdis_multiline
  #define rectangle   mcdis_rectangle
  #define box         mcdis_box
  #define circle      mcdis_circle
  #define cylinder    mcdis_cylinder
  #define sphere      mcdis_sphere
_class_Progress_bar *class_Progress_bar_display(_class_Progress_bar *_comp
) {
  #define profile (_comp->_parameters.profile)
  #define percent (_comp->_parameters.percent)
  #define flag_save (_comp->_parameters.flag_save)
  #define minutes (_comp->_parameters.minutes)
  #define IntermediateCnts (_comp->_parameters.IntermediateCnts)
  #define StartTime (_comp->_parameters.StartTime)
  #define EndTime (_comp->_parameters.EndTime)
  #define CurrentTime (_comp->_parameters.CurrentTime)

  #undef profile
  #undef percent
  #undef flag_save
  #undef minutes
  #undef IntermediateCnts
  #undef StartTime
  #undef EndTime
  #undef CurrentTime
  return(_comp);
} /* class_Progress_bar_display */

_class_Source_gen *class_Source_gen_display(_class_Source_gen *_comp
) {
  #define flux_file (_comp->_parameters.flux_file)
  #define xdiv_file (_comp->_parameters.xdiv_file)
  #define ydiv_file (_comp->_parameters.ydiv_file)
  #define radius (_comp->_parameters.radius)
  #define dist (_comp->_parameters.dist)
  #define focus_xw (_comp->_parameters.focus_xw)
  #define focus_yh (_comp->_parameters.focus_yh)
  #define focus_aw (_comp->_parameters.focus_aw)
  #define focus_ah (_comp->_parameters.focus_ah)
  #define E0 (_comp->_parameters.E0)
  #define dE (_comp->_parameters.dE)
  #define lambda0 (_comp->_parameters.lambda0)
  #define dlambda (_comp->_parameters.dlambda)
  #define I1 (_comp->_parameters.I1)
  #define yheight (_comp->_parameters.yheight)
  #define xwidth (_comp->_parameters.xwidth)
  #define verbose (_comp->_parameters.verbose)
  #define T1 (_comp->_parameters.T1)
  #define flux_file_perAA (_comp->_parameters.flux_file_perAA)
  #define flux_file_log (_comp->_parameters.flux_file_log)
  #define Lmin (_comp->_parameters.Lmin)
  #define Lmax (_comp->_parameters.Lmax)
  #define Emin (_comp->_parameters.Emin)
  #define Emax (_comp->_parameters.Emax)
  #define T2 (_comp->_parameters.T2)
  #define I2 (_comp->_parameters.I2)
  #define T3 (_comp->_parameters.T3)
  #define I3 (_comp->_parameters.I3)
  #define zdepth (_comp->_parameters.zdepth)
  #define target_index (_comp->_parameters.target_index)
  #define p_in (_comp->_parameters.p_in)
  #define lambda1 (_comp->_parameters.lambda1)
  #define lambda2 (_comp->_parameters.lambda2)
  #define lambda3 (_comp->_parameters.lambda3)
  #define pTable (_comp->_parameters.pTable)
  #define pTable_x (_comp->_parameters.pTable_x)
  #define pTable_y (_comp->_parameters.pTable_y)
  #define pTable_xmin (_comp->_parameters.pTable_xmin)
  #define pTable_xmax (_comp->_parameters.pTable_xmax)
  #define pTable_xsum (_comp->_parameters.pTable_xsum)
  #define pTable_ymin (_comp->_parameters.pTable_ymin)
  #define pTable_ymax (_comp->_parameters.pTable_ymax)
  #define pTable_ysum (_comp->_parameters.pTable_ysum)
  #define pTable_dxmin (_comp->_parameters.pTable_dxmin)
  #define pTable_dxmax (_comp->_parameters.pTable_dxmax)
  #define pTable_dymin (_comp->_parameters.pTable_dymin)
  #define pTable_dymax (_comp->_parameters.pTable_dymax)
  double xmin;
  double xmax;
  double ymin;
  double ymax;

  if (radius)
  {
    
    circle("xy",0,0,0,radius);
    if (zdepth) {
      circle("xy",0,0,-zdepth/2,radius);
      circle("xy",0,0, zdepth/2,radius);
    }
  }
  else
  {
    xmin = -xwidth/2; xmax = xwidth/2;
    ymin = -yheight/2; ymax = yheight/2;

    
    multiline(5, (double)xmin, (double)ymin, 0.0,
             (double)xmax, (double)ymin, 0.0,
             (double)xmax, (double)ymax, 0.0,
             (double)xmin, (double)ymax, 0.0,
             (double)xmin, (double)ymin, 0.0);
    if (zdepth) {
      multiline(5, (double)xmin, (double)ymin, -zdepth/2,
             (double)xmax, (double)ymin, -zdepth/2,
             (double)xmax, (double)ymax, -zdepth/2,
             (double)xmin, (double)ymax, -zdepth/2,
             (double)xmin, (double)ymin, -zdepth/2);
      multiline(5, (double)xmin, (double)ymin, zdepth/2,
             (double)xmax, (double)ymin, zdepth/2,
             (double)xmax, (double)ymax, zdepth/2,
             (double)xmin, (double)ymax, zdepth/2,
             (double)xmin, (double)ymin, zdepth/2);
    }
  }
  if (dist) {
    if (focus_aw) focus_xw=dist*tan(focus_aw*DEG2RAD);
    if (focus_ah) focus_yh=dist*tan(focus_ah*DEG2RAD);
    dashed_line(0,0,0, -focus_xw/2,-focus_yh/2,dist, 4);
    dashed_line(0,0,0,  focus_xw/2,-focus_yh/2,dist, 4);
    dashed_line(0,0,0,  focus_xw/2, focus_yh/2,dist, 4);
    dashed_line(0,0,0, -focus_xw/2, focus_yh/2,dist, 4);
  }
  #undef flux_file
  #undef xdiv_file
  #undef ydiv_file
  #undef radius
  #undef dist
  #undef focus_xw
  #undef focus_yh
  #undef focus_aw
  #undef focus_ah
  #undef E0
  #undef dE
  #undef lambda0
  #undef dlambda
  #undef I1
  #undef yheight
  #undef xwidth
  #undef verbose
  #undef T1
  #undef flux_file_perAA
  #undef flux_file_log
  #undef Lmin
  #undef Lmax
  #undef Emin
  #undef Emax
  #undef T2
  #undef I2
  #undef T3
  #undef I3
  #undef zdepth
  #undef target_index
  #undef p_in
  #undef lambda1
  #undef lambda2
  #undef lambda3
  #undef pTable
  #undef pTable_x
  #undef pTable_y
  #undef pTable_xmin
  #undef pTable_xmax
  #undef pTable_xsum
  #undef pTable_ymin
  #undef pTable_ymax
  #undef pTable_ysum
  #undef pTable_dxmin
  #undef pTable_dxmax
  #undef pTable_dymin
  #undef pTable_dymax
  return(_comp);
} /* class_Source_gen_display */

_class_Guide_gravity *class_Guide_gravity_display(_class_Guide_gravity *_comp
) {
  #define w1 (_comp->_parameters.w1)
  #define h1 (_comp->_parameters.h1)
  #define w2 (_comp->_parameters.w2)
  #define h2 (_comp->_parameters.h2)
  #define l (_comp->_parameters.l)
  #define R0 (_comp->_parameters.R0)
  #define Qc (_comp->_parameters.Qc)
  #define alpha (_comp->_parameters.alpha)
  #define m (_comp->_parameters.m)
  #define W (_comp->_parameters.W)
  #define nslit (_comp->_parameters.nslit)
  #define d (_comp->_parameters.d)
  #define mleft (_comp->_parameters.mleft)
  #define mright (_comp->_parameters.mright)
  #define mtop (_comp->_parameters.mtop)
  #define mbottom (_comp->_parameters.mbottom)
  #define nhslit (_comp->_parameters.nhslit)
  #define G (_comp->_parameters.G)
  #define aleft (_comp->_parameters.aleft)
  #define aright (_comp->_parameters.aright)
  #define atop (_comp->_parameters.atop)
  #define abottom (_comp->_parameters.abottom)
  #define wavy (_comp->_parameters.wavy)
  #define wavy_z (_comp->_parameters.wavy_z)
  #define wavy_tb (_comp->_parameters.wavy_tb)
  #define wavy_lr (_comp->_parameters.wavy_lr)
  #define chamfers (_comp->_parameters.chamfers)
  #define chamfers_z (_comp->_parameters.chamfers_z)
  #define chamfers_lr (_comp->_parameters.chamfers_lr)
  #define chamfers_tb (_comp->_parameters.chamfers_tb)
  #define nelements (_comp->_parameters.nelements)
  #define nu (_comp->_parameters.nu)
  #define phase (_comp->_parameters.phase)
  #define reflect (_comp->_parameters.reflect)
  #define GVars (_comp->_parameters.GVars)
  #define pTable (_comp->_parameters.pTable)

  if (l > 0 && nelements > 0) {
    int i,j,n;
    double x1,x2,x3,x4;
    double y1,y2,y3,y4;
    double nel = (nelements > 11 ? 11 : nelements);

    
    for (n=0; n<nel; n++)
    {
      double z0, z1;
      z0 =     n*(l/nel);
      z1 = (n+1)*(l/nel);

      for(j = 0; j < nhslit; j++)
      {
        y1 = j*(GVars.h1c+d)         - h1/2.0;
        y2 = j*(GVars.h2c+d)         - h2/2.0;
        y3 = (j+1)*(GVars.h1c+d) - d - h1/2.0;
        y4 = (j+1)*(GVars.h2c+d) - d - h2/2.0;
        for(i = 0; i < nslit; i++)
        {
          x1 = i*(GVars.w1c+d)         - w1/2.0;
          x2 = i*(GVars.w2c+d)         - w2/2.0;
          x3 = (i+1)*(GVars.w1c+d) - d - w1/2.0;
          x4 = (i+1)*(GVars.w2c+d) - d - w2/2.0;
          multiline(5,
                    x1, y1, z0,
                    x2, y2, z1,
                    x2, y4, z1,
                    x1, y3, z0,
                    x1, y1, z0);
          multiline(5,
                    x3, y1, z0,
                    x4, y2, z1,
                    x4, y4, z1,
                    x3, y3, z0,
                    x3, y1, z0);
        }
        line(-w1/2.0, y1, z0, w1/2.0, y1, z0);
        line(-w2/2.0, y2, z1, w2/2.0, y2, z1);
      }
    }

    if (nu || phase) {
      double radius = sqrt(w1*w1+l*l);
      /* cylinder top/center/bottom  */
      circle("xz", 0,-h1/2,l/2,radius);
      circle("xz", 0,0    ,l/2,radius);
      circle("xz", 0, h1/2,l/2,radius);
    }
  }
  else {
    /* A bit ugly; hard-coded dimensions. */
    
    line(0,0,0,0.2,0,0);
    line(0,0,0,0,0.2,0);
    line(0,0,0,0,0,0.2);
  }

  #undef w1
  #undef h1
  #undef w2
  #undef h2
  #undef l
  #undef R0
  #undef Qc
  #undef alpha
  #undef m
  #undef W
  #undef nslit
  #undef d
  #undef mleft
  #undef mright
  #undef mtop
  #undef mbottom
  #undef nhslit
  #undef G
  #undef aleft
  #undef aright
  #undef atop
  #undef abottom
  #undef wavy
  #undef wavy_z
  #undef wavy_tb
  #undef wavy_lr
  #undef chamfers
  #undef chamfers_z
  #undef chamfers_lr
  #undef chamfers_tb
  #undef nelements
  #undef nu
  #undef phase
  #undef reflect
  #undef GVars
  #undef pTable
  return(_comp);
} /* class_Guide_gravity_display */

_class_Monitor_nD *class_Monitor_nD_display(_class_Monitor_nD *_comp
) {
  #define user1 (_comp->_parameters.user1)
  #define user2 (_comp->_parameters.user2)
  #define user3 (_comp->_parameters.user3)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define zdepth (_comp->_parameters.zdepth)
  #define xmin (_comp->_parameters.xmin)
  #define xmax (_comp->_parameters.xmax)
  #define ymin (_comp->_parameters.ymin)
  #define ymax (_comp->_parameters.ymax)
  #define zmin (_comp->_parameters.zmin)
  #define zmax (_comp->_parameters.zmax)
  #define bins (_comp->_parameters.bins)
  #define min (_comp->_parameters.min)
  #define max (_comp->_parameters.max)
  #define restore_neutron (_comp->_parameters.restore_neutron)
  #define radius (_comp->_parameters.radius)
  #define options (_comp->_parameters.options)
  #define filename (_comp->_parameters.filename)
  #define geometry (_comp->_parameters.geometry)
  #define username1 (_comp->_parameters.username1)
  #define username2 (_comp->_parameters.username2)
  #define username3 (_comp->_parameters.username3)
  #define DEFS (_comp->_parameters.DEFS)
  #define Vars (_comp->_parameters.Vars)
  #define detector (_comp->_parameters.detector)
  #define offdata (_comp->_parameters.offdata)
  if (geometry && strlen(geometry) && strcmp(geometry,"0") && strcmp(geometry, "NULL"))
  {
    off_display(offdata);
  } else {
    Monitor_nD_McDisplay(&DEFS, &Vars);
  }
  #undef user1
  #undef user2
  #undef user3
  #undef xwidth
  #undef yheight
  #undef zdepth
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef zmin
  #undef zmax
  #undef bins
  #undef min
  #undef max
  #undef restore_neutron
  #undef radius
  #undef options
  #undef filename
  #undef geometry
  #undef username1
  #undef username2
  #undef username3
  #undef DEFS
  #undef Vars
  #undef detector
  #undef offdata
  return(_comp);
} /* class_Monitor_nD_display */

_class_DiskChopper *class_DiskChopper_display(_class_DiskChopper *_comp
) {
  #define theta_0 (_comp->_parameters.theta_0)
  #define radius (_comp->_parameters.radius)
  #define yheight (_comp->_parameters.yheight)
  #define nu (_comp->_parameters.nu)
  #define nslit (_comp->_parameters.nslit)
  #define jitter (_comp->_parameters.jitter)
  #define delay (_comp->_parameters.delay)
  #define isfirst (_comp->_parameters.isfirst)
  #define n_pulse (_comp->_parameters.n_pulse)
  #define abs_out (_comp->_parameters.abs_out)
  #define phase (_comp->_parameters.phase)
  #define xwidth (_comp->_parameters.xwidth)
  #define verbose (_comp->_parameters.verbose)
  #define Tg (_comp->_parameters.Tg)
  #define To (_comp->_parameters.To)
  #define delta_y (_comp->_parameters.delta_y)
  #define height (_comp->_parameters.height)
  #define omega (_comp->_parameters.omega)

  int j;
  /* Arrays for storing geometry of slit/beamstop */
  
  circle("xy", 0, -delta_y, 0, radius);

  /* Drawing the slit(s) */
  for (j=0; j<nslit; j++) {
    /* Angular start/end of slit */
    double tmin = j*(2.0*PI/nslit) - theta_0/2.0 + phase;
    double tmax = tmin+theta_0;
    /* Draw lines for each slit. */

    line(
      radius*sin(tmin),          radius*cos(tmin)-delta_y,          0,
      (radius-height)*sin(tmin), (radius-height)*cos(tmin)-delta_y, 0
      );
    line(
      (radius-height)*sin(tmin), (radius-height)*cos(tmin)-delta_y, 0,
      (radius-height)*sin(tmax), (radius-height)*cos(tmax)-delta_y, 0);
    line(
      (radius-height)*sin(tmax), (radius-height)*cos(tmax)-delta_y, 0,
      radius*sin(tmax),          radius*cos(tmax)-delta_y,          0);
  }
  #undef theta_0
  #undef radius
  #undef yheight
  #undef nu
  #undef nslit
  #undef jitter
  #undef delay
  #undef isfirst
  #undef n_pulse
  #undef abs_out
  #undef phase
  #undef xwidth
  #undef verbose
  #undef Tg
  #undef To
  #undef delta_y
  #undef height
  #undef omega
  return(_comp);
} /* class_DiskChopper_display */

_class_Arm *class_Arm_display(_class_Arm *_comp
) {
  /* A bit ugly; hard-coded dimensions. */
  
  line(0,0,0,0.2,0,0);
  line(0,0,0,0,0.2,0);
  line(0,0,0,0,0,0.2);
  return(_comp);
} /* class_Arm_display */

_class_Isotropic_Sqw *class_Isotropic_Sqw_display(_class_Isotropic_Sqw *_comp
) {
  #define powder_format (_comp->_parameters.powder_format)
  #define Sqw_coh (_comp->_parameters.Sqw_coh)
  #define Sqw_inc (_comp->_parameters.Sqw_inc)
  #define geometry (_comp->_parameters.geometry)
  #define radius (_comp->_parameters.radius)
  #define thickness (_comp->_parameters.thickness)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define zdepth (_comp->_parameters.zdepth)
  #define threshold (_comp->_parameters.threshold)
  #define order (_comp->_parameters.order)
  #define T (_comp->_parameters.T)
  #define verbose (_comp->_parameters.verbose)
  #define d_phi (_comp->_parameters.d_phi)
  #define concentric (_comp->_parameters.concentric)
  #define rho (_comp->_parameters.rho)
  #define sigma_abs (_comp->_parameters.sigma_abs)
  #define sigma_coh (_comp->_parameters.sigma_coh)
  #define sigma_inc (_comp->_parameters.sigma_inc)
  #define classical (_comp->_parameters.classical)
  #define powder_Dd (_comp->_parameters.powder_Dd)
  #define powder_DW (_comp->_parameters.powder_DW)
  #define powder_Vc (_comp->_parameters.powder_Vc)
  #define density (_comp->_parameters.density)
  #define weight (_comp->_parameters.weight)
  #define p_interact (_comp->_parameters.p_interact)
  #define norm (_comp->_parameters.norm)
  #define powder_barns (_comp->_parameters.powder_barns)
  #define quantum_correction (_comp->_parameters.quantum_correction)
  #define VarSqw (_comp->_parameters.VarSqw)
  #define columns (_comp->_parameters.columns)
  #define offdata (_comp->_parameters.offdata)
  if (VarSqw.s_coh > 0 || VarSqw.s_inc > 0) {

    if(VarSqw.shape==1)
    {
      double xmin = -0.5*xwidth;
      double xmax =  0.5*xwidth;
      double ymin = -0.5*yheight;
      double ymax =  0.5*yheight;
      double zmin = -0.5*zdepth;
      double zmax =  0.5*zdepth;
      multiline(5, xmin, ymin, zmin,
                  xmax, ymin, zmin,
                  xmax, ymax, zmin,
                  xmin, ymax, zmin,
                  xmin, ymin, zmin);
      multiline(5, xmin, ymin, zmax,
                  xmax, ymin, zmax,
                  xmax, ymax, zmax,
                  xmin, ymax, zmax,
                  xmin, ymin, zmax);
      line(xmin, ymin, zmin, xmin, ymin, zmax);
      line(xmax, ymin, zmin, xmax, ymin, zmax);
      line(xmin, ymax, zmin, xmin, ymax, zmax);
      line(xmax, ymax, zmin, xmax, ymax, zmax);

      if (thickness) {
        xmin = -0.5*xwidth+thickness;
        xmax = -xmin;
        ymin = -0.5*yheight+thickness;
        ymax = -ymin;
        zmin = -0.5*zdepth+thickness;
        zmax = -zmin;
        multiline(5, xmin, ymin, zmin,
                    xmax, ymin, zmin,
                    xmax, ymax, zmin,
                    xmin, ymax, zmin,
                    xmin, ymin, zmin);
        multiline(5, xmin, ymin, zmax,
                    xmax, ymin, zmax,
                    xmax, ymax, zmax,
                    xmin, ymax, zmax,
                    xmin, ymin, zmax);
        line(xmin, ymin, zmin, xmin, ymin, zmax);
        line(xmax, ymin, zmin, xmax, ymin, zmax);
        line(xmin, ymax, zmin, xmin, ymax, zmax);
        line(xmax, ymax, zmin, xmax, ymax, zmax);
      }
    }
    else if(VarSqw.shape==0)
    {
      circle("xz", 0,  yheight/2.0, 0, radius);
      circle("xz", 0, -yheight/2.0, 0, radius);
      line(-radius, -yheight/2.0, 0, -radius, +yheight/2.0, 0);
      line(+radius, -yheight/2.0, 0, +radius, +yheight/2.0, 0);
      line(0, -yheight/2.0, -radius, 0, +yheight/2.0, -radius);
      line(0, -yheight/2.0, +radius, 0, +yheight/2.0, +radius);
      if (thickness) {
        double radius_i=radius-thickness;
        circle("xz", 0,  yheight/2.0, 0, radius_i);
        circle("xz", 0, -yheight/2.0, 0, radius_i);
        line(-radius_i, -yheight/2.0, 0, -radius_i, +yheight/2.0, 0);
        line(+radius_i, -yheight/2.0, 0, +radius_i, +yheight/2.0, 0);
        line(0, -yheight/2.0, -radius_i, 0, +yheight/2.0, -radius_i);
        line(0, -yheight/2.0, +radius_i, 0, +yheight/2.0, +radius_i);
      }
    } else if(VarSqw.shape==2) {
      if (thickness) {
        double radius_i=radius-thickness;
        circle("xy",0,0,0,radius_i);
        circle("xz",0,0,0,radius_i);
        circle("yz",0,0,0,radius_i);
      }
      circle("xy",0,0,0,radius);
      circle("xz",0,0,0,radius);
      circle("yz",0,0,0,radius);
    } else if (VarSqw.shape == 3) {	/* OFF file */
      off_display(offdata);
    }
  }
/* end MCDISPLAY */
  #undef powder_format
  #undef Sqw_coh
  #undef Sqw_inc
  #undef geometry
  #undef radius
  #undef thickness
  #undef xwidth
  #undef yheight
  #undef zdepth
  #undef threshold
  #undef order
  #undef T
  #undef verbose
  #undef d_phi
  #undef concentric
  #undef rho
  #undef sigma_abs
  #undef sigma_coh
  #undef sigma_inc
  #undef classical
  #undef powder_Dd
  #undef powder_DW
  #undef powder_Vc
  #undef density
  #undef weight
  #undef p_interact
  #undef norm
  #undef powder_barns
  #undef quantum_correction
  #undef VarSqw
  #undef columns
  #undef offdata
  return(_comp);
} /* class_Isotropic_Sqw_display */

_class_PSD_Detector *class_PSD_Detector_display(_class_PSD_Detector *_comp
) {
  #define nx (_comp->_parameters.nx)
  #define ny (_comp->_parameters.ny)
  #define xwidth (_comp->_parameters.xwidth)
  #define radius (_comp->_parameters.radius)
  #define awidth (_comp->_parameters.awidth)
  #define yheight (_comp->_parameters.yheight)
  #define zdepth (_comp->_parameters.zdepth)
  #define threshold (_comp->_parameters.threshold)
  #define PressureConv (_comp->_parameters.PressureConv)
  #define PressureStop (_comp->_parameters.PressureStop)
  #define interpolate (_comp->_parameters.interpolate)
  #define p_interact (_comp->_parameters.p_interact)
  #define verbose (_comp->_parameters.verbose)
  #define LensOn (_comp->_parameters.LensOn)
  #define dc (_comp->_parameters.dc)
  #define borderx (_comp->_parameters.borderx)
  #define bordery (_comp->_parameters.bordery)
  #define xChDivRelSigma (_comp->_parameters.xChDivRelSigma)
  #define yChDivRelSigma (_comp->_parameters.yChDivRelSigma)
  #define bufsize (_comp->_parameters.bufsize)
  #define restore_neutron (_comp->_parameters.restore_neutron)
  #define angle (_comp->_parameters.angle)
  #define type (_comp->_parameters.type)
  #define filename (_comp->_parameters.filename)
  #define FN_Conv (_comp->_parameters.FN_Conv)
  #define FN_Stop (_comp->_parameters.FN_Stop)
  #define PSD_N (_comp->_parameters.PSD_N)
  #define PSD_p (_comp->_parameters.PSD_p)
  #define PSD_p2 (_comp->_parameters.PSD_p2)
  #define EAP (_comp->_parameters.EAP)
  #define M1P1 (_comp->_parameters.M1P1)
  #define PosAP (_comp->_parameters.PosAP)
  #define EAT (_comp->_parameters.EAT)
  #define M1T1 (_comp->_parameters.M1T1)
  #define PosAT (_comp->_parameters.PosAT)
  #define CrossSectionHe (_comp->_parameters.CrossSectionHe)
  #define CountNeutrons (_comp->_parameters.CountNeutrons)
  #define GeomCumul (_comp->_parameters.GeomCumul)
  #define AbsCumul (_comp->_parameters.AbsCumul)
  #define SensVolCumul (_comp->_parameters.SensVolCumul)
  #define DetCumul (_comp->_parameters.DetCumul)
  #define PHSpectrum (_comp->_parameters.PHSpectrum)
  #define PHSpectrum0 (_comp->_parameters.PHSpectrum0)
  #define PHSpectrum2 (_comp->_parameters.PHSpectrum2)
  #define PHSpectrum_n (_comp->_parameters.PHSpectrum_n)
  #define nH_p (_comp->_parameters.nH_p)
  #define nH_t (_comp->_parameters.nH_t)
  #define FullEnergyP (_comp->_parameters.FullEnergyP)
  #define FullEnergyT (_comp->_parameters.FullEnergyT)
  #define VariousErrors (_comp->_parameters.VariousErrors)
  #define DetectorType (_comp->_parameters.DetectorType)
  #define DEFS (_comp->_parameters.DEFS)
  #define Vars (_comp->_parameters.Vars)
  double h;
  h=yheight;

  if (xwidth>0) { /* box */
    double xmin = -0.5*xwidth;
    double xmax =  0.5*xwidth;
    double ymin = -0.5*yheight;
    double ymax =  0.5*yheight;
    double zmin = -0.5*zdepth;
    double zmax =  0.5*zdepth;
    multiline(5, xmin, ymin, zmin,
                 xmax, ymin, zmin,
                 xmax, ymax, zmin,
                 xmin, ymax, zmin,
                 xmin, ymin, zmin);
    multiline(5, xmin, ymin, zmax,
                 xmax, ymin, zmax,
                 xmax, ymax, zmax,
                 xmin, ymax, zmax,
                 xmin, ymin, zmax);
    line(xmin, ymin, zmin, xmin, ymin, zmax);
    line(xmax, ymin, zmin, xmax, ymin, zmax);
    line(xmin, ymax, zmin, xmin, ymax, zmax);
    line(xmax, ymax, zmin, xmax, ymax, zmax);
  }
  if (radius>0) {
		if (!awidth) { /* cylinder */
	    circle("xz", 0,  h/2.0, 0, radius);
	    circle("xz", 0, -h/2.0, 0, radius);
	    line(-radius, -h/2.0, 0, -radius, +h/2.0, 0);
	    line(+radius, -h/2.0, 0, +radius, +h/2.0, 0);
	    line(0, -h/2.0, -radius, 0, +h/2.0, -radius);
	    line(0, -h/2.0, +radius, 0, +h/2.0, +radius);
	    if (borderx>0){
	       circle("xz", 0,  h/2.0, 0, radius+borderx);
	       circle("xz", 0, -h/2.0, 0, radius+borderx);
	    }
    } else {
    	int NH=24;
      int ih;

			for(ih = 0; ih < NH; ih++) {
		    double phi0, phi1;
		    double x0,y0,z0,x1,y1,z1,x2,y2,z2,x3,y3,z3;
		    phi0 = (-angle/2+(angle/NH)*ih)    *DEG2RAD; /* in xz plane */
		    phi1 = (-angle/2+(angle/NH)*(ih+1))*DEG2RAD;

		    z0 = radius*cos(phi0);
		    x0 = radius*sin(phi0);
		    y0 = -yheight/2;
		    z1 = radius*cos(phi0);
		    x1 = radius*sin(phi0);
		    y1 = yheight/2;
		    z2 = radius*cos(phi1);
		    x2 = radius*sin(phi1);
		    y2 = y1;
		    z3 = radius*cos(phi1);
		    x3 = radius*sin(phi1);
		    y3 = y0;
		    mcdis_multiline(5,
		      x0,y0,z0,
		      x1,y1,z1,
		      x2,y2,z2,
		      x3,y3,z3,
		      x0,y0,z0);
			}
    }
  }
  #undef nx
  #undef ny
  #undef xwidth
  #undef radius
  #undef awidth
  #undef yheight
  #undef zdepth
  #undef threshold
  #undef PressureConv
  #undef PressureStop
  #undef interpolate
  #undef p_interact
  #undef verbose
  #undef LensOn
  #undef dc
  #undef borderx
  #undef bordery
  #undef xChDivRelSigma
  #undef yChDivRelSigma
  #undef bufsize
  #undef restore_neutron
  #undef angle
  #undef type
  #undef filename
  #undef FN_Conv
  #undef FN_Stop
  #undef PSD_N
  #undef PSD_p
  #undef PSD_p2
  #undef EAP
  #undef M1P1
  #undef PosAP
  #undef EAT
  #undef M1T1
  #undef PosAT
  #undef CrossSectionHe
  #undef CountNeutrons
  #undef GeomCumul
  #undef AbsCumul
  #undef SensVolCumul
  #undef DetCumul
  #undef PHSpectrum
  #undef PHSpectrum0
  #undef PHSpectrum2
  #undef PHSpectrum_n
  #undef nH_p
  #undef nH_t
  #undef FullEnergyP
  #undef FullEnergyT
  #undef VariousErrors
  #undef DetectorType
  #undef DEFS
  #undef Vars
  return(_comp);
} /* class_PSD_Detector_display */


  #undef magnify
  #undef line
  #undef dashed_line
  #undef multiline
  #undef rectangle
  #undef box
  #undef circle
  #undef cylinder
  #undef sphere

int display(void) { /* called by mccode_main for HZB_NEAT:DISPLAY */
  printf("MCDISPLAY: start\n");

  /* call iteratively all components DISPLAY */
  class_Progress_bar_display(_PG);

  class_Source_gen_display(_Source);

  class_Guide_gravity_display(_Guide0);

  class_Guide_gravity_display(_Guide0_4);

  class_Guide_gravity_display(_Guide0_5);

  class_Guide_gravity_display(_Guide0_6);

  class_Guide_gravity_display(_Guide0_7);

  class_Guide_gravity_display(_Guide0_8);

  class_Guide_gravity_display(_Guide0_9);

  class_Guide_gravity_display(_Guide0_10);

  class_Guide_gravity_display(_Guide0_11);

  class_Guide_gravity_display(_Guide0_12);

  class_Guide_gravity_display(_Guide0_13);

  class_Guide_gravity_display(_Guide0_14);

  class_Guide_gravity_display(_Guide0_15);

  class_Guide_gravity_display(_Guide0_16);

  class_Guide_gravity_display(_Guide0_17);

  class_Guide_gravity_display(_Guide0_18);

  class_Guide_gravity_display(_Guide0_19);

  class_Monitor_nD_display(_Guide_PSD);

  class_Monitor_nD_display(_Guide_Lambda);

  class_DiskChopper_display(_Chopper1);

  class_Guide_gravity_display(_Guide_SGS1);

  class_Guide_gravity_display(_Guide_SGS2);

  class_Guide_gravity_display(_Guide_SGS3);

  class_Guide_gravity_display(_Guide_CGS);

  class_Monitor_nD_display(_Guide2_PSD);

  class_Monitor_nD_display(_Guide2_Lambda);

  class_Monitor_nD_display(_Guide2_dXY);

  class_DiskChopper_display(_Chopper6);

  class_Monitor_nD_display(_Guide_time);

  class_Guide_gravity_display(_Guide_DGS);

  class_Arm_display(_Sample_pos);

  class_Isotropic_Sqw_display(_Sample);

  class_Monitor_nD_display(_Ideal_Det);

  class_PSD_Detector_display(_Detector);

  printf("MCDISPLAY: end\n");

  return(0);
} /* display */

void* _getvar_parameters(char* compname)
/* enables settings parameters based use of the GETPAR macro */
{
  if (!strcmp(compname, "PG")) return (void *) &(_PG->_parameters);
  if (!strcmp(compname, "Source")) return (void *) &(_Source->_parameters);
  if (!strcmp(compname, "Guide0")) return (void *) &(_Guide0->_parameters);
  if (!strcmp(compname, "Guide0_4")) return (void *) &(_Guide0_4->_parameters);
  if (!strcmp(compname, "Guide0_5")) return (void *) &(_Guide0_5->_parameters);
  if (!strcmp(compname, "Guide0_6")) return (void *) &(_Guide0_6->_parameters);
  if (!strcmp(compname, "Guide0_7")) return (void *) &(_Guide0_7->_parameters);
  if (!strcmp(compname, "Guide0_8")) return (void *) &(_Guide0_8->_parameters);
  if (!strcmp(compname, "Guide0_9")) return (void *) &(_Guide0_9->_parameters);
  if (!strcmp(compname, "Guide0_10")) return (void *) &(_Guide0_10->_parameters);
  if (!strcmp(compname, "Guide0_11")) return (void *) &(_Guide0_11->_parameters);
  if (!strcmp(compname, "Guide0_12")) return (void *) &(_Guide0_12->_parameters);
  if (!strcmp(compname, "Guide0_13")) return (void *) &(_Guide0_13->_parameters);
  if (!strcmp(compname, "Guide0_14")) return (void *) &(_Guide0_14->_parameters);
  if (!strcmp(compname, "Guide0_15")) return (void *) &(_Guide0_15->_parameters);
  if (!strcmp(compname, "Guide0_16")) return (void *) &(_Guide0_16->_parameters);
  if (!strcmp(compname, "Guide0_17")) return (void *) &(_Guide0_17->_parameters);
  if (!strcmp(compname, "Guide0_18")) return (void *) &(_Guide0_18->_parameters);
  if (!strcmp(compname, "Guide0_19")) return (void *) &(_Guide0_19->_parameters);
  if (!strcmp(compname, "Guide_PSD")) return (void *) &(_Guide_PSD->_parameters);
  if (!strcmp(compname, "Guide_Lambda")) return (void *) &(_Guide_Lambda->_parameters);
  if (!strcmp(compname, "Chopper1")) return (void *) &(_Chopper1->_parameters);
  if (!strcmp(compname, "Guide_SGS1")) return (void *) &(_Guide_SGS1->_parameters);
  if (!strcmp(compname, "Guide_SGS2")) return (void *) &(_Guide_SGS2->_parameters);
  if (!strcmp(compname, "Guide_SGS3")) return (void *) &(_Guide_SGS3->_parameters);
  if (!strcmp(compname, "Guide_CGS")) return (void *) &(_Guide_CGS->_parameters);
  if (!strcmp(compname, "Guide2_PSD")) return (void *) &(_Guide2_PSD->_parameters);
  if (!strcmp(compname, "Guide2_Lambda")) return (void *) &(_Guide2_Lambda->_parameters);
  if (!strcmp(compname, "Guide2_dXY")) return (void *) &(_Guide2_dXY->_parameters);
  if (!strcmp(compname, "Chopper6")) return (void *) &(_Chopper6->_parameters);
  if (!strcmp(compname, "Guide_time")) return (void *) &(_Guide_time->_parameters);
  if (!strcmp(compname, "Guide_DGS")) return (void *) &(_Guide_DGS->_parameters);
  if (!strcmp(compname, "Sample_pos")) return (void *) &(_Sample_pos->_parameters);
  if (!strcmp(compname, "Sample")) return (void *) &(_Sample->_parameters);
  if (!strcmp(compname, "Ideal_Det")) return (void *) &(_Ideal_Det->_parameters);
  if (!strcmp(compname, "Detector")) return (void *) &(_Detector->_parameters);
}

void* _get_particle_var(char *token, _class_particle *p)
/* enables setpars based use of GET_PARTICLE_DVAR macro and similar */
{
  return 0;
}

/* embedding file "mccode_main.c" */

/*******************************************************************************
* mccode_main: McCode main() function.
*******************************************************************************/
int mccode_main(int argc, char *argv[])
{
  /*  double run_num = 0; */
  time_t  t;
  clock_t ct;


#ifdef USE_MPI
  char mpi_node_name[MPI_MAX_PROCESSOR_NAME];
  int  mpi_node_name_len;
#endif /* USE_MPI */

#ifdef MAC
  argc = ccommand(&argv);
#endif

#ifdef USE_MPI
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_node_count); /* get number of nodes */
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_node_rank);
  MPI_Comm_set_name(MPI_COMM_WORLD, instrument_name);
  MPI_Get_processor_name(mpi_node_name, &mpi_node_name_len);
#endif /* USE_MPI */


  ct = clock();    /* we use clock rather than time to set the default seed */
  mcseed=(long)ct;


#ifdef USE_MPI
  /* *** print number of nodes *********************************************** */
  if (mpi_node_count > 1) {
    MPI_MASTER(
    printf("Simulation '%s' (%s): running on %i nodes (master is '%s', MPI version %i.%i).\n",
      instrument_name, instrument_source, mpi_node_count, mpi_node_name, MPI_VERSION, MPI_SUBVERSION);
    );
    /* share the same seed, then adapt random seed for each node */
    MPI_Bcast(&mcseed, 1, MPI_LONG, 0, MPI_COMM_WORLD); /* root sends its seed to slaves */
    mcseed += mpi_node_rank; /* make sure we use different seeds per noe */
  }
#endif /* USE_MPI */


  srandom(mcseed);
  mcstartdate = (long)t;  /* set start date before parsing options and creating sim file */

  /* *** parse options ******************************************************* */
  SIG_MESSAGE("[" __FILE__ "] main START");
  mcformat = getenv(FLAVOR_UPPER "_FORMAT") ?
             getenv(FLAVOR_UPPER "_FORMAT") : FLAVOR_UPPER;
  instrument_exe = argv[0]; /* store the executable path */
  /* read simulation parameters and options */
  mcparseoptions(argc, argv); /* sets output dir and format */


#ifdef USE_MPI
  if (mpi_node_count > 1) {
    /* share the same seed, then adapt random seed for each node */
    MPI_Bcast(&mcseed, 1, MPI_LONG, 0, MPI_COMM_WORLD); /* root sends its seed to slaves */
    mcseed += mpi_node_rank; /* make sure we use different seeds per node */
  }
#endif


  srandom(mcseed);


/* *** install sig handler, but only once !! after parameters parsing ******* */
#ifndef NOSIGNALS
#ifdef SIGQUIT
  if (signal( SIGQUIT ,sighandler) == SIG_IGN)
    signal( SIGQUIT,SIG_IGN);   /* quit (ASCII FS) */
#endif
#ifdef SIGABRT
  if (signal( SIGABRT ,sighandler) == SIG_IGN)
    signal( SIGABRT,SIG_IGN);   /* used by abort, replace SIGIOT in the future */
#endif
#ifdef SIGTERM
  if (signal( SIGTERM ,sighandler) == SIG_IGN)
    signal( SIGTERM,SIG_IGN);   /* software termination signal from kill */
#endif
#ifdef SIGUSR1
  if (signal( SIGUSR1 ,sighandler) == SIG_IGN)
    signal( SIGUSR1,SIG_IGN);   /* display simulation status */
#endif
#ifdef SIGUSR2
  if (signal( SIGUSR2 ,sighandler) == SIG_IGN)
    signal( SIGUSR2,SIG_IGN);
#endif
#ifdef SIGHUP
  if (signal( SIGHUP ,sighandler) == SIG_IGN)
    signal( SIGHUP,SIG_IGN);
#endif
#ifdef SIGILL
  if (signal( SIGILL ,sighandler) == SIG_IGN)
    signal( SIGILL,SIG_IGN);    /* illegal instruction (not reset when caught) */
#endif
#ifdef SIGFPE
  if (signal( SIGFPE ,sighandler) == SIG_IGN)
    signal( SIGSEGV,SIG_IGN);    /* floating point exception */
#endif
#ifdef SIGBUS
  if (signal( SIGBUS ,sighandler) == SIG_IGN)
    signal( SIGSEGV,SIG_IGN);    /* bus error */
#endif
#ifdef SIGSEGV
  if (signal( SIGSEGV ,sighandler) == SIG_IGN)
    signal( SIGSEGV,SIG_IGN);   /* segmentation violation */
#endif
#endif /* !NOSIGNALS */


  // init
  siminfo_init(NULL); /* open SIM */
  SIG_MESSAGE("[" __FILE__ "] main INITIALISE");
  init();


#ifndef NOSIGNALS
#ifdef SIGINT
  if (signal( SIGINT ,sighandler) == SIG_IGN)
    signal( SIGINT,SIG_IGN);    /* interrupt (rubout) only after INIT */
#endif
#endif /* !NOSIGNALS */

/* ================ main particle generation/propagation loop ================ */
#ifdef USE_MPI
  /* sliced Ncount on each MPI node */
  mcncount = mpi_node_count > 1 ?
    floor(mcncount / mpi_node_count) :
    mcncount; /* number of rays per node */
#endif

/* main particle event loop */
#ifdef USE_PGI
#include <openacc.h>
#include "mccode_attaches.c"
#endif


  #pragma acc parallel loop
  /* old init: mcsetstate(0, 0, 0, 0, 0, 1, 0, sx=0, sy=1, sz=0, 1); */
  for (unsigned long long Xmcrun_num=0 ; Xmcrun_num < mcncount ; Xmcrun_num++) {

    _class_particle particleN = mcgenstate(); // initial particle
    particleN._uid = Xmcrun_num;


/* CUDA */
#ifdef USE_PGI
    curandState_t MCRANDstate;
    /* To really make sense this should be long long, but for now 
       compilation for GPU seems to need longs only */
    //long long seq = Xmcrun_num;
    long seq = Xmcrun_num+mcseed;
    curand_init(seq, seq-mcseed, 0ULL, &MCRANDstate);
    particleN.MCRANDstate = MCRANDstate;
#endif


    raytrace(&particleN);
  }

  /* Likely we need an undef random here... */


#ifdef USE_MPI
 /* merge run_num from MPI nodes */
  if (mpi_node_count > 1) {
  double mcrun_num_double = (double)mcrun_num;
  mc_MPI_Sum(&mcrun_num_double, 1);
  mcrun_num = (unsigned long long)mcrun_num_double;
  }
#endif


  // save/finally executed by master node/thread
  finally();


#ifdef USE_MPI
  MPI_Finalize();
#endif /* USE_MPI */


  return 0;
} /* mccode_main */
/* End of file "mccode_main.c". */

/* end of generated C code HZB_NEAT.c */