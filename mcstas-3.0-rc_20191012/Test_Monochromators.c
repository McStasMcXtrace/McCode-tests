/* Automatically generated file. Do not edit. 
 * Format:     ANSI C source code
 * Creator:    McStas <http://www.mcstas.org>
 * Instrument: Test_Monochromators.instr (Test_Monochromators)
 * Date:       Sat Oct 12 09:06:16 2019
 * File:       Test_Monochromators.c
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
* Start of instrument 'Test_Monochromators' generated code
***************************************************************************** */

#ifdef MC_TRACE_ENABLED
int traceenabled = 1;
#else
int traceenabled = 0;
#endif
#define MCSTAS "/usr/share/mcstas/3.0-dev/"
int   defaultmain         = 1;
char  instrument_name[]   = "Test_Monochromators";
char  instrument_source[] = "Test_Monochromators.instr";
char *instrument_exe      = NULL; /* will be set to argv[0] in main */
char  instrument_code[]   = "Instrument Test_Monochromators source code Test_Monochromators.instr is not embedded in this executable.\n  Use --source option when running McStas.\n";

int main(int argc, char *argv[]){return mccode_main(argc, argv);}

/* *****************************************************************************
* instrument 'Test_Monochromators' and components DECLARE
***************************************************************************** */

/* Instrument parameters: structure and a table for the initialisation
   (Used in e.g. inputparse and I/O function (e.g. detector_out) */

struct _struct_instrument_parameters {
  long _Mono;
  MCNUM _lambda;
  MCNUM _Moz;
  MCNUM _PG;
  MCNUM _powder;
  MCNUM _order;
};
typedef struct _struct_instrument_parameters _class_instrument_parameters;

struct _instrument_struct {
  char   _name[256]; /* the name of this instrument e.g. 'Test_Monochromators' */
/* Counters per component instance */
  double counter_AbsorbProp[16]; /* absorbed events in PROP routines */
  double counter_N[16], counter_P[16], counter_P2[16]; /* event counters after each component instance */
  _class_particle _trajectory[16]; /* current trajectory for STORE/RESTORE */
/* Components position table (absolute and relative coords) */
  Coords _position_relative[16]; /* positions of all components */
  Coords _position_absolute[16];
  _class_instrument_parameters _parameters; /* instrument parameters */
} _instrument_var;
struct _instrument_struct *instrument = & _instrument_var;
#pragma acc declare create ( _instrument_var )
#pragma acc declare create ( instrument )

int numipar = 6;
struct mcinputtable_struct mcinputtable[] = {
  "Mono", &(_instrument_var._parameters._Mono), instr_type_int, "1", 
  "lambda", &(_instrument_var._parameters._lambda), instr_type_double, "2.0", 
  "Moz", &(_instrument_var._parameters._Moz), instr_type_double, "40", 
  "PG", &(_instrument_var._parameters._PG), instr_type_double, "0", 
  "powder", &(_instrument_var._parameters._powder), instr_type_double, "0", 
  "order", &(_instrument_var._parameters._order), instr_type_double, "1", 
  NULL, NULL, instr_type_double, ""
};


/* ************************************************************************** */
/*             SHARE user declarations for all components                     */
/* ************************************************************************** */

/* Shared user declarations for all components types 'Monochromator_flat'. */
#ifndef GAUSS
/* Define these arrays only once for all instances. */
/* Values for Gauss quadrature. Taken from Brice Carnahan, H. A. Luther and
James O Wilkes, "Applied numerical methods", Wiley, 1969, page 103.
This reference is available from the Copenhagen UB2 library */
double Gauss_X[] = {-0.987992518020485, -0.937273392400706, -0.848206583410427,
-0.724417731360170, -0.570972172608539, -0.394151347077563,
-0.201194093997435, 0, 0.201194093997435,
0.394151347077563, 0.570972172608539, 0.724417731360170,
0.848206583410427, 0.937273392400706, 0.987992518020485};
double Gauss_W[] = {0.030753241996117, 0.070366047488108, 0.107159220467172,
0.139570677926154, 0.166269205816994, 0.186161000115562,
0.198431485327111, 0.202578241925561, 0.198431485327111,
0.186161000115562, 0.166269205816994, 0.139570677926154,
0.107159220467172, 0.070366047488108, 0.030753241996117};


#define GAUSS(x,mean,rms) \
  (exp(-((x)-(mean))*((x)-(mean))/(2*(rms)*(rms)))/(sqrt(2*PI)*(rms)))
#endif

/* Shared user declarations for all components types 'Monochromator_pol'. */
/*****************************************************************************
*
* McStas, neutron ray-tracing package
*         Copyright 1997-2006, All rights reserved
*         Risoe National Laboratory, Roskilde, Denmark
*         Institut Laue Langevin, Grenoble, France
*
* Library: share/pol-lib.h
*
* %Identification
* Written by: Peter Christiansen
* Date: August, 2006
* Origin: RISOE
* Release: McStas 1.10
* Version: $Revision: 4382 $
*
* This file is to be imported by polarisation components.
* It handles some shared functions.
*
* This library may be used directly as an external library. 
* It has no dependency.
*
* Usage: within SHARE
* %include "pol-lib"
*
****************************************************************************/

#ifndef POL_LIB_H
#define POL_LIB_H "$Revision: 4382 $"

// Constant used 
#define mc_pol_omegaL (-2 * PI * 29.16e6) /* MHz*rad/Tesla */
#define mc_pol_mu0 (4*M_PI*1e-7)

/*example field functions should have a variable set of arguments*/
#include <stdarg.h>
#include <stddef.h>
/*macros for some stuff*/
#ifndef MCSTAS_R_H
#include <mcstas-r.h>
#endif


typedef int mcmagnet_field_func (double, double, double, double, double *, double *, double *, void *);
typedef void mcmagnet_prec_func (double, double, double, double, double, double, double, double*, double*, double*, double, Coords, Rotation);
typedef va_list mcmagnet_data;

/*here's where the mcstas magnet stack is declared*/
/*the magnet stack*/

typedef struct mcmagnet_field_info {
  mcmagnet_field_func *func;
  Rotation *rot;
  Coords *pos;
  void *data;
  int stop;
} mcmagnet_field_info;

void mc_pol_set_timestep(double);
void mc_pol_set_angular_accuracy(double);

#define mcmagnet_sizeof (sizeof(mcmagnet_field_func *)+ sizeof(Rotation *)+ sizeof(Coords *)+ sizeof(double *))
#define mcmagnet_malloc(n) malloc( (n)*sizeof(mcmagnet_field_info) );

#define mcmagnet_pack(dest,funk,rotation,position,stopbit,args) \
  do { \
    mcmagnet_field_info * mctmp_p; \
    mctmp_p=(dest); \
    mctmp_p->func=(mcmagnet_field_func *)(funk); \
    mctmp_p->rot=(rotation); \
    mctmp_p->pos=(position); \
    mctmp_p->stop=(stopbit); \
    mctmp_p->data=(args); \
  } while (0);

#define mcmagnet_reset() \
  do { \
    mcMagneticField=NULL; \
    mcMagnetData=NULL; \
    MAGNET_OFF; \
  } while (0);

#define mcmagnet_set_active(mcmagnet_new) \
  do { \
    if (mcmagnet_new!=NULL){ \
      mcMagneticField=(mcmagnet_new)->func; \
      rot_copy(mcMagnetRot, *((mcmagnet_new)->rot)); \
      mcMagnetPos=*((mcmagnet_new)->pos); \
      mcMagnetData=(double *)(mcmagnet_new)->data; \
    }else{ \
      mcmagnet_reset(); \
    } \
  } while (0);

#define mcmagnet_free(mcmagnet_desc) \
  do { \
    mcmagnet_field_info * mctmp_p=(mcmagnet_desc); \
    if (mctmp_p!=NULL) { \
      if (mctmp_p->data!=NULL) free(mctmp_p->data); \
      free(mctmp_p); \
    } \
  } while(0);

#define MCMAGNET_STOP_ARG INT_MIN

#define mcmagnet_init_par(...) \
  mcmagnet_init_par_backend(0, __VA_ARGS__, MCMAGNET_STOP_ARG);

void mcmagnet_print_active();
void mcmagnet_print_field(mcmagnet_field_info *);
void mcmagnet_print_stack();

void *mcmagnet_init_par_backend(int dummy, ...);

int mcmagnet_get_field(double x, double y, double z, double t, double *bx,double *by, double *bz, void *dummy);
void *mcmagnet_push(mcmagnet_field_func *func,  Rotation *magnet_rot, Coords *magnet_pos, int stopbit, void * prms);
void *mcmagnet_pop(void);

/*example functions for magnetic fields*/
int const_magnetic_field(double x, double y, double z, double t, double *bx, double *by, double *bz, void *data);
int rot_magnetic_field(double x, double y, double z, double t, double *bx, double *by, double *bz, void *data);
int majorana_magnetic_field(double x, double y, double z, double t, double *bx, double *by, double *bz, void *data);
int table_magnetic_field(double x, double y, double z, double t,
                         double *bx, double *by, double *bz,
                         void *data);

/* Routines used for Monochromator and guides/mirrors 
 * in the special (usual) case where
 * the up direction is parallel to the y-axis and 
 * the down direction is anti-parallel to the y-axis */
void GetMonoPolFNFM(double, double, double*, double*);
void GetMonoPolRefProb(double, double, double, double*);
void SetMonoPolRefOut(double, double, double, double*, double*, double*);
void SetMonoPolTransOut(double, double, double, double*, double*, double*);

// Routines for spin precession in magnetic fields
void SimpleNumMagnetPrecession(double, double, double, double, double, double, 
			       double, double*, double*, double*, double, 
			       Coords, Rotation);

void SimpleNumMagnetPrecession___(double, double, double, double, double, double, 
			       double, double*, double*, double*, double, 
			       Coords, Rotation);
void SeegerNumMagnetPrecession(double, double, double, double, double, double, 
			       double, double*, double*, double*, double, 
			       Coords, Rotation);


// Routines to help calculate the rquired magnetic field
double GetConstantField(double, double, double);

#endif

/* end of pol-lib.h */
/****************************************************************************
*
* McStas, neutron ray-tracing package
*         Copyright 1997-2006, All rights reserved
*         Risoe National Laboratory, Roskilde, Denmark
*         Institut Laue Langevin, Grenoble, France
*
* Library: share/pol-lib.c
*
* %Identification
* Written by: Erik Knudsen, Astrid Rømer & Peter Christiansen
* Date: Oct 08
* Origin: RISOE
* Release: McStas 1.12
* Version: $Revision: 4466 $
*
* This file is to be imported by polarisation components.
* It handles some shared functions.
* Embedded within instrument in runtime mode.
* Variable names have prefix 'mc_pol_' for 'McStas Polarisation'
* to avoid conflicts
*
* Usage: within SHARE
* %include "pol-lib"
*
****************************************************************************/

#ifndef POL_LIB_H
#include "pol-lib.h"
#endif

#include<sys/stat.h>

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

/*******************************************************************************
*
* McCode, neutron/xray ray-tracing package
*         Copyright (C) 1997-2015, All rights reserved
*         Risoe National Laboratory, Roskilde, Denmark
*         Institut Laue Langevin, Grenoble, France
*
* Runtime: share/interpolation.h
*
* %Identification
* Written by: EF
* Date:    May 5th 2015
* Release: McStas X.Y/McXtrace X.Y
* Version: $Revision: 5455 $
*
* Table interpolation routines (header)
*
* Usage: Automatically embbeded in the c code whenever required, with e.g.:
*   %include "interpolation-lib"
*
* public function:
* interpolator = interpolator_load(filename, 0, 0, NULL);
*   or
* interpolator = interpolator_load(filename, space_dim, field_dim, "regular" or "kdtree");
*
* interpolator_info(interpolator);
* 
* interpolator_interpolate(interpolator, {x,y,z...}, {bx,by,bz...});
*   or 
* interpolator_interpolate3_3(interpolator, x,y,z, &bx,&by,&bz);
* 
* interpolator_save(interpolator);
*
* Example:
*   struct interpolator_struct interpolator = 
*             interpolator_load("filename", space_dim, field_dim, NULL);
*   interpolator_info(interpolator);
*   double space[space_dim]={x,y,z};
*   double field[field_dim]; // will contain interpolated values
*   interpolator_interpolate(interpolator, space, field); 
*
* Data file format:
* file is a list of rows [x,y,z...    field_x, field_y, ... ]
*                        | space ... | field  ... |
*/

/*******************************************************************************
 * begin declaration (.h) section
 ******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

typedef struct
{
  // This is the location of this point (space).
  short  space_dimensionality;
  double *v;    // e.g. double []
  
  // These are the values for our field at this location.
  double *data; // e.g. double []

  // This is the point index in the point list.
  int    index;

} vertex;
 
/* This struct will store each node of our kdtree. */
typedef struct _treeNode {
  vertex   *point;
  int       depth;
  struct _treeNode *rChild;
  struct _treeNode *lChild;
} treeNode;


#define INTERPOLATOR_DIMENSIONS 10
  
struct interpolator_struct {
  char  method[256];
  long  space_dimensionality; // [x,y,z...]
  long  field_dimensionality; // [bx,by,bz...]
  long  points;
  char  filename[1024];
  treeNode *kdtree;    /* for k-d tree */
  double  *grid[INTERPOLATOR_DIMENSIONS];  /* each grid contains a component of the field */
  double   min[INTERPOLATOR_DIMENSIONS];
  double   max[INTERPOLATOR_DIMENSIONS];
  long     bin[INTERPOLATOR_DIMENSIONS];
  double   step[INTERPOLATOR_DIMENSIONS];
  long     constant_step[INTERPOLATOR_DIMENSIONS];
};

#undef INTERPOLATOR_DIMENSIONS

/******************************************************************************/
// interpolator_info: print information about the interpolator
void interpolator_info(struct interpolator_struct *interpolator);
 
/*******************************************************************************
 * interpolator_load: interpolation initialiser, from point cloud
 *   returns the interpolator structure
 * The input is mainly the file name, which is a column based text format.
 * The interpolator->method is set as 'kdtree' or 'regular' as set at points load
 ******************************************************************************/ 
struct interpolator_struct *interpolator_load(char *filename, 
   long space_dimensionality, long field_dimensionality,
   char *method);
     
/*******************************************************************************
 * interpolator_interpolate: main interpolation routine.
 *   returns the 'field' value (of length interpolator->field_dimensionality)
 *   at the given 'space' location (of length interpolator->space_dimensionality)
 *   The returned array 'field' MUST be pre-allocated.
 ******************************************************************************/ 
double *interpolator_interpolate(struct interpolator_struct *interpolator,
  double *space, double *field);


/*******************************************************************************
 * interpolator_interpolate3_3: main interpolation routine for 3D space
 *   returns the 'field' value (e.g. 3d)
 *   at the given 'coord' location (e.g. 3d)
 * The interpolator->method can be 'kdtree' or 'regular' as set at points load
 ******************************************************************************/ 
double *interpolator_interpolate3_3(struct interpolator_struct *interpolator,
                    double  x,  double  y,  double  z,
                    double *bx, double *by, double *bz);

/*******************************************************************************
*
* McCode, neutron/xray ray-tracing package
*         Copyright (C) 1997-2015, All rights reserved
*         Risoe National Laboratory, Roskilde, Denmark
*         Institut Laue Langevin, Grenoble, France
*
* Runtime: share/interpolation.c
*
* %Identification
* Written by: EF
* Date:    May 5th 2015
* Release: McStas X.Y/McXtrace X.Y
* Version: $Revision: 5455 $
*
* Table interpolation routines
*
* Usage: Automatically embbeded in the c code whenever required, with e.g.:
*   %include "interpolation-lib"
*
* public function:
* interpolator = interpolator_load(filename, 0, 0, NULL);
*   or
* interpolator = interpolator_load(filename, space_dim, field_dim, "regular" or "kdtree");
*
* interpolator_info(interpolator);
* 
* interpolator_interpolate(interpolator, {x,y,z...}, {bx,by,bz...});
*   or 
* interpolator_interpolate3_3(interpolator, x,y,z, &bx,&by,&bz);
* 
* interpolator_save(interpolator);
*
* Example:
*   struct interpolator_struct interpolator = 
*             interpolator_load("filename", space_dim, field_dim, NULL);
*   interpolator_info(interpolator);
*   double space[space_dim]={x,y,z};
*   double field[field_dim]; // will contain interpolated values
*   interpolator_interpolate(interpolator, space, field); 
*
* Data file format:
* file is a list of rows [x,y,z...    field_x, field_y, ... ]
*                        | space ... | field  ... |
*/

/*******************************************************************************
 * begin declaration (.h) section
 ******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/*******************************************************************************
 * begin k-D tree section
 ******************************************************************************/

#define R_SQR(x)        ((x) * (x))
#define R_SWAP(x, y, t) {t tmp; tmp=x; x=y; y=tmp;}
 
 
/******************************************************************************/

// kdtree_squaredDistance: Calculate the standard Euclidean distance between 
//   these two points in whatever dimension we are considering.
double kdtree_squaredDistance(vertex* a, vertex* b)
{
  int i;
  double sum = 0;
  if (!a || !b || a->space_dimensionality != b->space_dimensionality) return 0;
  
  for (i = 0; i < a->space_dimensionality; i++) {
    sum += R_SQR(a->v[i] - b->v[i]);
  }
  return sum;
} // kdtree_squaredDistance

/******************************************************************************/
// kdtree_borderCheck: Check to see whether or not this node provides a better 
//   nearest neighbour.
void kdtree_borderCheck(vertex *v, treeNode *thisNode,
                 vertex **currentBest, double *sDist)
{
  if (!thisNode || !v || !sDist) return;
  
  double thisDist = kdtree_squaredDistance(thisNode->point,v);
  if (thisDist < *sDist)
  {
    *sDist        = thisDist;
    *currentBest  = thisNode->point;
  }
  // Now recurse down the children, checking whether or not we should
  // go both sides of the splitting plane, or just down one side.
  int k = (thisNode->depth) % v->space_dimensionality;
  if (R_SQR(thisNode->point->v[k] - v->v[k]) <= *sDist)
  {
   // The distance to the current spliting plane is less than our current
   // estimate for the shortest distance, we are going to have to traverse
   // both sides of the splitting plane.
    kdtree_borderCheck(v, thisNode->lChild, currentBest, sDist);
    kdtree_borderCheck(v, thisNode->rChild, currentBest, sDist);
  } else {
    // We only have to consider one side of the splitting plane.
    if (thisNode->point->v[k] > (*currentBest)->v[k])
      kdtree_borderCheck(v, thisNode->lChild, currentBest, sDist);
    else
      kdtree_borderCheck(v, thisNode->rChild, currentBest, sDist);
  }
} // kdtree_borderCheck

/******************************************************************************/

// kdtree_partition: Note we slightly modify the standard partition algorithm, 
//   so that we can partition based on only one dimension of the pointset.
int kdtree_partition(vertex **points, int d, int left, int right, int pivot)
{
  double pivotValue = points[pivot]->v[d];
  int i;
  int storeIndex = left;
  
  if (!points) return 0;

  R_SWAP(points[pivot], points[right], vertex*);

  for (i = left; i < right; i++) {
    if (points[i]->v[d] < pivotValue) {
      R_SWAP(points[storeIndex], points[i], vertex*);
      storeIndex ++;
    }
  }
  R_SWAP(points[right], points[storeIndex], vertex*);

  return storeIndex;
} // kdtree_partition

/******************************************************************************/
// kdtree_splitAboutMedian: Find the median in expected linear time. - We will 
//   also pivot all the data about the found median, returning the integer giving
//   the pivot value.

int kdtree_splitAboutMedian(vertex **points, int d, int left, int right)
{
  int k = (right-left)/2 +left;
  if (!points) return 0;
  
  // This isn't a perfect uniform distribution, but it doesn't really matter
  // for this application.
  while (left < right)
  {
    int pivotIndex    = rand() % (right-left)+left;
    int pivotNewIndex = kdtree_partition(points,d,left,right,pivotIndex);
    if (k == pivotNewIndex)
      return k;
    else if (k < pivotNewIndex)
      right = pivotNewIndex-1;
    else
      left  = pivotNewIndex+1;
  }

  return left;
} // kdtree_splitAboutMedian

/******************************************************************************/
// kdtree_addToTree: create a kd-tree out of a point set
treeNode* kdtree_addToTree(vertex **points, int left, int right, int depth)
{
  // We can modify the number of dimensions in use. This is defined in the
  // header file.

  if (right < left || !points) return NULL;

  int d = depth % points[0]->space_dimensionality;

  treeNode *node = malloc(sizeof(treeNode));
  node->depth    = depth;

  int med      = kdtree_splitAboutMedian(points, d, left, right);
  node->point  = points[med];

  node->lChild = kdtree_addToTree(points, left,  med-1, depth + 1);
  node->rChild = kdtree_addToTree(points, med+1, right, depth + 1);

  return node;
} // kdtree_addToTree

/******************************************************************************/
// kdtree_nearestNeighbour_helper: helper function for kdtree_nearestNeighbour
//   used recursively until a close vertex is found
void kdtree_nearestNeighbour_helper(vertex* v, treeNode *tree,
                             vertex **bestV, double *bestDist)
{
  if (!v || !tree || !bestDist) return;
  
  int k = tree->depth % v->space_dimensionality;

  int left = tree->point->v[k] > v->v[k];

  treeNode *first  = left ? tree->lChild : tree->rChild;
  treeNode *second = left ? tree->rChild : tree->lChild;

  // investigate first child if present
  if (first != NULL) {
    kdtree_nearestNeighbour_helper(v, first, bestV, bestDist);
  }

  // update result
  double thisDist = kdtree_squaredDistance(tree->point, v);
  if ((*bestV == NULL) || (thisDist < *bestDist)) {
    *bestDist = thisDist;
    *bestV    = tree->point;
  }

  // no second child to investigate
  if (second == NULL) {
    return;
  }

  // we only investigate second child if necessary
  int treek = tree->point->v[k];

  if (R_SQR(treek - v->v[k]) <= *bestDist) {
    kdtree_borderCheck(v, second, bestV, bestDist);
  }
} // kdtree_nearestNeighbour_helper

/******************************************************************************/
// kdtree_nearestNeighbour: find closest vertex in tree to given vertex coords
vertex* kdtree_nearestNeighbour(vertex* v, treeNode *tree) {
  vertex *bestV = NULL;
  double bestDist = 0;
  if (!v || !tree) return NULL;

  kdtree_nearestNeighbour_helper(v, tree, &bestV, &bestDist);
  v->data = bestV->data;
  
  return bestV;
} // kdtree_nearestNeighbour

#undef R_SQR
#undef R_SWAP

/*******************************************************************************
 * end k-D tree section
 ******************************************************************************/


/*******************************************************************************
 * begin interpolator section
 ******************************************************************************/
 
#define INTERPOLATOR_DIMENSIONS 10


/******************************************************************************/
/* interpolator_double_vector_compare: comparator for double qsort */
int interpolator_double_vector_compare(void const *a, void const *b) {
  return ( *(double*)a - *(double*)b );
}

/******************************************************************************/
/* interpolator_init: initialise an empty interpolator structure */
struct interpolator_struct *interpolator_init(void) {
  int dim=0;
  struct interpolator_struct *interpolator = malloc(sizeof(struct interpolator_struct));
  
  if (!interpolator) return NULL;
  
  strcpy(interpolator->method,"NULL");
  strcpy(interpolator->filename,"NULL");
  interpolator->points = interpolator->space_dimensionality 
                       = interpolator->field_dimensionality = 0;
  interpolator->kdtree = NULL;
  for (dim=0; dim < INTERPOLATOR_DIMENSIONS; dim++) {
    interpolator->min[dim] = +FLT_MAX;
    interpolator->max[dim] = -FLT_MAX;
    interpolator->bin[dim] = 0;
    interpolator->step[dim]= 0;
    interpolator->constant_step[dim] = 1; /* assumes we have constant step. Check done at load. */
    interpolator->grid[dim] = NULL;
  }
  return interpolator;
} /* interpolator_init */

/******************************************************************************/
// interpolator_offset: determine element offset for an n-dimensional array
//   used in: interpolator_load and interpolator_interpolate
long interpolator_offset(int dim, long *dimInfo, long *indices) {
  
  long result;  // where the resultant offset will be stored 
  int  i;       // loop counter 
  
  /* indices check */
  for (i=0; i < dim; i++) {
    if (indices[i] < 0)           indices[i]=0;
    if (indices[i] >= dimInfo[i]) indices[i]=dimInfo[i]-1;
  }
  // Perform the general offset calculation for an n-dimensional array 
  for (i=0; i < dim; i++) {
    result = i == 0 ? indices[0]
                    : result * dimInfo[i] + indices[i];
  }
  return result; 
} // interpolator_offset

/******************************************************************************/
// interpolator_info: print information about the interpolator
void interpolator_info(struct interpolator_struct *interpolator) {
  if (!interpolator) return;
  MPI_MASTER(
    printf("interpolator: file '%s' with %ld points. Space is %ldD, Field is %ldD. Using method '%s'.\n",
      interpolator->filename, interpolator->points, 
      interpolator->space_dimensionality, interpolator->field_dimensionality,
      interpolator->method);
  );
} /* interpolator_info */
 
/*******************************************************************************
 * interpolator_load: interpolation initialiser, from point cloud
 *   returns the interpolator structure
 * The input is mainly the file name, which is a column based text format.
 * The interpolator->method is set as 'kdtree' or 'regular' as set at points load
 ******************************************************************************/ 
struct interpolator_struct *interpolator_load(char *filename, 
   long space_dimensionality, long field_dimensionality,
   char *method) {

  struct interpolator_struct *interpolator = interpolator_init();
  int dim=0;
  
  // Read the table with Read Table Lib
  t_Table table;

  if(!Table_Read(&table, filename, 0) || table.rows <= 0 || !filename || strlen(filename) > 1024) {
    // Give up!
    fprintf(stderr, "interpolator_load: ERROR: Could not open file: '%s'.\n", filename);
    Table_Free(&table);
    return NULL;
  }
  
  strcpy(interpolator->filename, filename);
  interpolator->space_dimensionality = space_dimensionality;
  interpolator->field_dimensionality = field_dimensionality;
  interpolator->points = table.rows; /* rows= [x,y,z,... field_x, field_y, ... ] */
  if (method && strlen(method) && strlen(method) < 32)
    strcpy(interpolator->method, method);
  else
    strcpy(interpolator->method, "NULL");
  
  /* get columns and determine dimensionality if not set */
  if (!interpolator->space_dimensionality) {
    if (table.columns >= 4)
      interpolator->space_dimensionality=3;
    else if (table.columns == 2)
      interpolator->space_dimensionality=1;
  }
  if (interpolator->space_dimensionality <= 0 
   || interpolator->space_dimensionality > INTERPOLATOR_DIMENSIONS) {
    fprintf(stderr, "interpolator_load: ERROR: Invalid space dimensionality "
                    "(0 < dim=%li < %i) from file '%s'.\n",
      interpolator->space_dimensionality, INTERPOLATOR_DIMENSIONS, filename);
    return NULL;
  }
  
  interpolator->field_dimensionality = table.columns - space_dimensionality;
  if (interpolator->field_dimensionality <= 0 
   || interpolator->field_dimensionality > INTERPOLATOR_DIMENSIONS) {
    fprintf(stderr, "interpolator_load: ERROR: Invalid field dimensionality "
                    "(0 < dim=%li < %i) from file '%s'.\n",
      interpolator->field_dimensionality, INTERPOLATOR_DIMENSIONS, filename);
    return NULL;
  }
  
  /* read space columns to determine if sampling is regular */
  for (dim=0; dim<interpolator->space_dimensionality; dim++) {
    double  x_prev=0;
    long    index;
    double  vector[table.rows];
    
    /* get min/max and fill vector for sorting */
    for (index=0; index<table.rows; index++) {
      double x = Table_Index(table, index, dim);
      if (x < interpolator->min[dim]) interpolator->min[dim] = x;
      if (x > interpolator->max[dim]) interpolator->max[dim] = x;
      vector[index] = x;
    }
    /* sort vector */
    qsort(vector, table.rows, sizeof(double), interpolator_double_vector_compare);
    
    /* now count the number of unique values and check constant step */
    for (index=0; index<table.rows; index++) {
      double x = vector[index];
      double this_step = 0;
      if (!index) x_prev = x;
      this_step = fabs(x - x_prev);
      if (this_step)
        interpolator->bin[dim]++; /* count unique values */
      if (interpolator->step[dim] <= 0) 
        interpolator->step[dim] = this_step;
      if (fabs(this_step - interpolator->step[dim]) > interpolator->step[dim]*READ_TABLE_STEPTOL) {
        /* difference of this step with the first one is 'large' */
        interpolator->constant_step[dim] = 0; /* not constant step -> kd-tree should be used */
        if (!strcmp(interpolator->method, "NULL") || !strcmp(interpolator->method, "0"))
          strcpy(interpolator->method, "kdtree");          
      }
      x_prev = x;
    }
    printf("interpolator_load: Axis %d: step=%g, unique values=%li, from file '%s'.\n",
        dim, interpolator->step[dim], interpolator->bin[dim], filename);

    if (interpolator->step[dim]<=0 || interpolator->bin[dim]<=1) {
      fprintf(stderr, "interpolator_load: ERROR: Invalid axis %d: step=%g, unique values=%li, from file '%s'.\n",
        dim, interpolator->step[dim], interpolator->bin[dim], filename);
      strcpy(interpolator->method,"NULL");
      return NULL;
    }
  } /* end for dim(space/axis) */

  /* check kd-tree method */
  if (!strlen(interpolator->method) || !strcmp(interpolator->method, "NULL") || !strcmp(interpolator->method, "0"))
    if (strcmp(interpolator->method, "kdtree"))  /* not kdtree ? -> use direct indexing */
      strcpy(interpolator->method, "regular");
  
  /* assign interpolation technique: 'regular' direct indexing */
  if (!strcmp(interpolator->method, "regular")) {
    interpolator->kdtree = NULL;
    /* store table values onto the grid: each field component is stored on the
     * interpolator->grid, and has size=prod(interpolator->bin)
     */
    
    long prod=1; /* the number of elements in the grid */
    for (dim=0; dim<interpolator->space_dimensionality; dim++)
      prod *= interpolator->bin[dim];
    for (dim=0; dim<interpolator->field_dimensionality; dim++) {
      double *array = (double*)calloc(prod, sizeof(double));
      printf("interpolator_load: allocating %g Gb for dim=%d\n",
        (double)prod/1073741824.0, dim); fflush(NULL);
      long index;
      if (!array) {
        fprintf(stderr, "interpolator_load: ERROR: Not enough memory for field component %i\n"
                        "  which requires %g Gb, from file '%s'. Will use kd-tree method.\n",
        dim, (double)prod/1073741824.0, filename);
        strcpy(interpolator->method,"kdtree");
        break;
      }
      for (index=0; index<table.rows; index++) {
        long indices[interpolator->space_dimensionality];
        long this_index;
        int  axis=0;

        /* compute index 'space' elements of this 'field' value */
        for (axis=0; axis < interpolator->space_dimensionality; axis++) {
          double x      = Table_Index(table, index, axis);
          indices[axis] = floor((x - interpolator->min[axis])/interpolator->step[axis]);
        }
        this_index = interpolator_offset(interpolator->space_dimensionality,
                       interpolator->bin, indices);
        // array[axis1][axis2][...] = field[dim] column after [space] elements
        array[this_index] = Table_Index(table, index, interpolator->space_dimensionality+dim);
      }
      interpolator->grid[dim] = array;
    } // end for dim(field)
  } 

  /* assign interpolation technique: kd-tree (when nearest direct indexing fails) */
  if (!strcmp(interpolator->method, "kdtree")) {
    // Allocate array of vertex pointers
    vertex **vertices = calloc(table.rows, sizeof(vertex*));
    if (!vertices) {
      fprintf(stderr, "interpolator_load: ERROR: Not enough memory when allocating field with %li vertices from file '%s'\n",
        interpolator->bin[dim], filename);
      strcpy(interpolator->method,"NULL");
      return NULL;
    }

    // Convert from table to array layout
    int i, j;
    long count=0;
    for (i=0; i < table.rows; i++)
    {
      vertex *v    = malloc(sizeof(vertex));
      double *field= calloc(interpolator->field_dimensionality, sizeof(double));
      double *coord= calloc(interpolator->space_dimensionality, sizeof(double));
      if (v && field && coord) {
        for (j = 0; j < interpolator->space_dimensionality; j++) {
          coord[j]    = Table_Index(table, i,     j);
        }
        for (j = 0; j < interpolator->field_dimensionality; j++) {
          field[j] = Table_Index(table, i, interpolator->space_dimensionality + j);
        }
        v->space_dimensionality = interpolator->space_dimensionality;
        v->v    = coord;
        v->data = field;
        v->index= i;
      }
      vertices[i] = v;
    }

    interpolator->kdtree = kdtree_addToTree(vertices, 0, table.rows-1, 0); // build treeNode
    for (i=0; i<INTERPOLATOR_DIMENSIONS; interpolator->grid[i++] = NULL);  // inactivate grid method
    free(vertices);
  } 
  else
    fprintf(stderr, "interpolator_load: ERROR: unknown interpolator method %s [file '%s'].\n",
      interpolator->method, filename);
  
  // Free table
  Table_Free(&table);
  return interpolator;
} /* end interpolator_load */
     
/*******************************************************************************
 * interpolator_interpolate: main interpolation routine.
 *   returns the 'field' value (of length interpolator->field_dimensionality)
 *   at the given 'space' location (of length interpolator->space_dimensionality)
 *   The returned array 'field' MUST be pre-allocated.
 ******************************************************************************/ 
double *interpolator_interpolate(struct interpolator_struct *interpolator,
  double *space, double *field)
{
  if (!space || !interpolator || !field) return NULL;
  
  /* k-d tree call ************************************************************/
  if (!strcmp(interpolator->method, "kdtree") && interpolator->kdtree) {
    vertex v;
    int i;
    v.v = space; 
    v.space_dimensionality=interpolator->space_dimensionality;
    vertex *w =kdtree_nearestNeighbour(&v, interpolator->kdtree);
    if (!w) return NULL;
    for (i=0; i<interpolator->field_dimensionality; field[i]=w->data[i++]);
    return (w->data);

  } else 
  
  /* nearest direct grid element call *****************************************/
  if (!strcmp(interpolator->method, "regular") && interpolator->grid[0]) {
    int axis;
    long indices[interpolator->space_dimensionality];
    for (axis=0; axis < interpolator->space_dimensionality; axis++) {
      indices[axis] = (space[axis]-interpolator->min[axis])/interpolator->step[axis];
    }
    long index = interpolator_offset(3, interpolator->bin, indices);
    for (axis=0; axis < interpolator->field_dimensionality; axis++) {
      field[axis] = interpolator->grid[axis][index];
    }
    return field;
  } else {
    fprintf(stderr, "interpolator_interpolate: ERROR: invalid interpolator method %s from file '%s'.\n",
      interpolator->method, interpolator->filename);
    exit(-1);
  }
  
} // interpolator_interpolate


/*******************************************************************************
 * interpolator_interpolate3_3: main interpolation routine for 3D space
 *   returns the 'field' value (e.g. 3d)
 *   at the given 'coord' location (e.g. 3d)
 * The interpolator->method can be 'kdtree' or 'regular' as set at points load
 ******************************************************************************/ 
double *interpolator_interpolate3_3(struct interpolator_struct *interpolator,
                    double  x,  double  y,  double  z,
                    double *bx, double *by, double *bz)
{
  double coord[3] = { x,y,z };
  double field[3] = { 0,0,0 };
  double *ret=NULL;
  if (interpolator->space_dimensionality != 3 
   || interpolator->field_dimensionality != 3) return 0;
  ret = interpolator_interpolate(interpolator, coord, field);
  *bx = field[0]; *by = field[1]; *bz = field[2];
  return(ret);
} /* interpolator_interpolate3_3 */

#undef INTERPOLATOR_DIMENSIONS



enum {MCMAGNET_STACKSIZE=12} mcmagnet_constants;

/*definition of the magnetic stack*/
static mcmagnet_field_info *stack[MCMAGNET_STACKSIZE];
/*definition of the precession function*/
#ifdef MC_POL_COMPAT
extern mcmagnet_prec_func *mcMagnetPrecession;
extern Coords   mcMagnetPos;
extern Rotation mcMagnetRot;
extern double*  mcMagnetData;
/* mcMagneticField(x, y, z, t, Bx, By, Bz) */
extern int (*mcMagneticField) (double, double, double, double,
    double*, double*, double*, void *);
#else
#ifndef POL_LIB_C
static mcmagnet_prec_func *mcMagnetPrecession=SimpleNumMagnetPrecession;
static Coords mcMagnetPos;
static Rotation mcMagnetRot;
static double*  mcMagnetData                = NULL;
static int (*mcMagneticField) (double, double, double, double,
    double*, double*, double*, void *);
/*Threshold below which two magnetic fields are considered to be
 * in the same direction.*/
static double mc_pol_angular_accuracy = 1.0*DEG2RAD; /*rad.*/
/*The maximal timestep taken by neutrons in a const field*/
static double mc_pol_initial_timestep = 1e-5;
#define POL_LIB_C 1
#endif
#endif

int mcmagnet_init(){
  mcMagnetPrecession=SimpleNumMagnetPrecession;
  return 1;
}

void mc_pol_set_angular_accuracy(double domega){
    mc_pol_angular_accuracy = domega;
}

void mc_pol_set_timestep(double dt){
    mc_pol_initial_timestep=dt;
}

#ifdef PROP_MAGNET
#undef PROP_MAGNET
#define PROP_MAGNET(dt) \
  do { \
    /* change coordinates from local system to magnet system */ \
    Rotation rotLM; \
    Coords   posLM = POS_A_CURRENT_COMP; \
    rot_transpose(ROT_A_CURRENT_COMP, rotLM); \
    mcMagnetPrecession(x, y, z, t, vx, vy, vz, \
	   	       &sx, &sy, &sz, dt, posLM, rotLM); \
  } while(0)
#endif

/*traverse the stack and return the magnetic field*/
int mcmagnet_get_field(double x, double y, double z, double t, double *bx,double *by, double *bz, void *dummy){
  mcmagnet_field_info *p=stack[0];
  Coords in,loc,b,bsum={0,0,0},zero={0,0,0};
  Rotation r;

  /*PROP_MAGNET takes care of transforming local "PROP" coordinates to lab system*/
  in.x=x;in.y=y;in.z=z;

  int i=0,stat=1;
  p=stack[i];
  *bx=0;*by=0;*bz=0;
  if (!p) return 0;
  //mcmagnet_print_stack();
  //printf("getfield_(lab):_(xyz,t)=( %g %g %g %g )\n",x,y,z,t);
  while(p){
    /*transform to the coordinate system of the particular magnetic function*/
    loc=coords_sub(rot_apply(*(p->rot),in),*(p->pos));
    stat=(p->func) (loc.x,loc.y,loc.z,t,&(b.x),&(b.y),&(b.z),p->data);
    /*check if the field function should be garbage collected*/
    //printf("getfield_(loc):_(xyz,t)=( %g %g %g %g )\n",loc.x,loc.y,loc.z,t);
    if (stat){
      /*transform to the lab system and add up. (resusing loc variable - to now contain the field in lab coords)*/
      rot_transpose(*(p->rot),r);
      loc=rot_apply(r,b);
      bsum.x+=loc.x;bsum.y+=loc.y;bsum.z+=loc.z;
      //printf("Bs=(%g %g %g), B=(%g %g %g)\n",bsum.x,bsum.y,bsum.z,loc.x,loc.y,loc.z);
    }
    if (p->stop) break;
    p=stack[++i];
  }
  /*we now have the magnetic field in lab coords in loc., transfer it back to caller*/
  *bx=bsum.x;
  *by=bsum.y;
  *bz=bsum.z;
  return 1;
}

/*void mcmagnet_init(void){
  mcmagnet_field_info *p;
  for (p=&(stack[0]);p<&(stack[MCMAGNET_STACKSIZE]);p++){
    *p = malloc (sizeof(mcmagnet_field_info));
  }
}
*/
void *mcmagnet_push(mcmagnet_field_func *func,  Rotation *magnet_rot, Coords *magnet_pos, int stopbit, void * prms){
  mcmagnet_field_info *p;
  int i;
  /*move the stack one step down start from -2 since we have 0-indexing (i.e. last item is stacksize-1) */
  for (i=MCMAGNET_STACKSIZE-2;i>=0;i--){
    stack[i+1]=stack[i];
  }
  stack[0]=(mcmagnet_field_info *)malloc(sizeof(mcmagnet_field_info));
  mcmagnet_pack(stack[0],func,magnet_rot,magnet_pos,stopbit,prms);
  mcmagnet_set_active(stack[0]);
  if(stack[0] && stack[0]->func){
    MAGNET_ON;
  }
  return (void *) stack[0];
}

void *mcmagnet_pop(void) {
  mcmagnet_field_info **p,*t;
  /*move the stack one step up*/
  int i;
  t=stack[0];
  for (i=0;i<MCMAGNET_STACKSIZE-2;i++){
    stack[i]=stack[i+1];
  }
  stack[MCMAGNET_STACKSIZE-1]=NULL;
  mcmagnet_set_active(stack[0]);
  if(stack[0] && stack[0]->func){
    MAGNET_ON;
  }else{
    MAGNET_OFF;
  }
  return (void*) t;
}

void mcmagnet_free_stack(void){
  mcmagnet_field_info **p;
  for (p=&(stack[0]);p<&(stack[MCMAGNET_STACKSIZE]);p++){
    free(*p);
  }
}

void *mcmagnet_init_par_backend(int dummy, ...){
  void * data;
  unsigned char *p=NULL;
  int q,dp=0;
  va_list arg_list;

  va_start(arg_list,dummy);
  p=(unsigned char *)arg_list;
  q=va_arg(arg_list,int);
  while (q!=MCMAGNET_STOP_ARG){
    q=va_arg(arg_list,int);
  }
  dp=(unsigned char *)arg_list-p;
  data=(void *) malloc(sizeof(int)*dp);
  memcpy(data,p,sizeof(int)*dp);
  return data;
}

void mcmagnet_print_active(){
  Rotation *p;
  printf("address of magnetic field function:%p\n",mcMagneticField);
  p=&mcMagnetRot;
  printf("rotation matrix of magnetic field:[%g %g %g; %g %g %g; %g %g %g]\n",(*p)[0][0],(*p)[0][1],(*p)[0][2],(*p)[1][0],(*p)[1][1],(*p)[1][2],(*p)[2][0],(*p)[2][1],(*p)[2][2]);
  printf("origin position of magnet (x,y,z) :[%g %g %g]\n",mcMagnetPos.x,mcMagnetPos.y,mcMagnetPos.z);
  printf("address of magnetic field parameters: %p\n",mcMagnetData);
}

void mcmagnet_print_field(mcmagnet_field_info *magnet){
  Rotation *p;
  if (magnet!=NULL){
    printf("address of magnetic field function:%p\n",magnet->func);
    p=magnet->rot;
    printf("rotation matrix of magnetic field:[%g %g %g; %g %g %g; %g %g %g]\n",(*p)[0][0],(*p)[0][1],(*p)[0][2],(*p)[1][0],(*p)[1][1],(*p)[1][2],(*p)[2][0],(*p)[2][1],(*p)[2][2]);
    printf("origin position of magnet (x,y,z) :[%g %g %g]\n",magnet->pos->x,magnet->pos->y,magnet->pos->z);
    printf("address of magnetic field parameters: %p\n",magnet->data);
  } else {
    printf("magnet is NULL\n");
  }
}

void mcmagnet_print_stack(){
  mcmagnet_field_info *p=stack[0];
  int i=0;
  p=stack[i];
  printf("magnetic stack info:\n");
  if (!p) return;
  while(p) {
    printf("magnet %d:\n",i);
    mcmagnet_print_field(p);
    if (p->stop) break;
    p=stack[++i];
  }
}


/*Example magnetic field functions*/
int const_magnetic_field(double x, double y, double z, double t,
    double *bx, double *by, double *bz, void *data) {
  int stat=1;
  if (!data) return 0;
  *bx=((double *)data)[0];
  *by=((double *)data)[1];
  *bz=((double *)data)[2];
  return stat;
}

int rot_magnetic_field(double x, double y, double z, double t,
    double *bx, double *by, double *bz, void *data) {
  /* Field of magnitude By that rotates to x in magnetLength m*/
  
  if (!data) return 0;
  double Bmagnitude=((double *)data)[0];//   = mcMagnetData[1];
  double magnetLength=((double *)data)[1];// = mcMagnetData[5];
  *bx =  Bmagnitude * sin(PI/2*z/magnetLength);
  *by =  Bmagnitude * cos(PI/2*z/magnetLength);
  *bz =  0;
  //printf("mag field at (x,y,z)=( %g %g %g ) t=%g is B=( %g %g %g )\n",x,y,z,t,*bx,*by,*bz);
  return 1;
}

int majorana_magnetic_field(double x, double y, double z, double t,
    double *bx, double *by, double *bz, void *data) {
  /* Large linearly decreasing (from +Bx to -Bx in magnetLength) component along x axis,
   * small constant component along y axis
   */
  if (!data) return 0;
  double Blarge       = ((double *)data)[0];
  double Bsmall       = ((double *)data)[1];
  double magnetLength = ((double *)data)[2];
  *bx =  Blarge -2*Blarge*z/magnetLength;
  *by =  Bsmall;
  *bz =  0;
  return 1;
}

int table_magnetic_field(double x, double y, double z, double t,
                         double *bx, double *by, double *bz,
                         void *data)
{
  if (!data) return 0;
  struct interpolator_struct *interpolator = (struct interpolator_struct*)data;
  return(interpolator_interpolate3_3(interpolator, x,y,z, bx,by,bz) != NULL);
}


/****************************************************************************
* void GetMonoPolFNFM(double Rup, double Rdown, double *FN, double *FM)
*
* ACTION: Calculate FN and FM from reflectivities Rup and Rdown
*
* For a monochromator (nuclear and magnetic scattering), the
* polarisation is done by defining the reflectivity for spin up (Rup)
* and spin down (Rdown) (which can be negative, see now!) and based on
* this the nuclear and magnetic structure factors are calculated:
* FM = sign(Rup)*sqrt(|Rup|) + sign(Rdown)*sqrt(|Rdown|)
* FN = sign(Rup)*sqrt(|Rup|) - sign(Rdown)*sqrt(|Rdown|)
*****************************************************************************/
void GetMonoPolFNFM(double mc_pol_Rup, double mc_pol_Rdown,
		    double *mc_pol_FN, double *mc_pol_FM) {
  if (mc_pol_Rup>0)
    mc_pol_Rup   = sqrt(fabs(mc_pol_Rup));
  else
    mc_pol_Rup   = -sqrt(fabs(mc_pol_Rup));

  if (mc_pol_Rdown>0)
    mc_pol_Rdown = sqrt(fabs(mc_pol_Rdown));
  else
    mc_pol_Rdown = -sqrt(fabs(mc_pol_Rdown));

  *mc_pol_FN = 0.5*(mc_pol_Rup + mc_pol_Rdown);
  *mc_pol_FM = 0.5*(mc_pol_Rup - mc_pol_Rdown);
  return;
}

/****************************************************************************
* void GetMonoPolRefProb(double FN, double FM, double sy, double *prob)
*
* ACTION: Calculate reflection probability from sy, FN and FM
*
* For a monochromator with up direction along y the reflection
* probability is given as:
* prob = FN*FN + 2*FN*FM*sy_in + FM*FM
*     (= |Rup| + |Rdown| (for sy_in=0))
* where FN and FM are calculated from Rup and Rdown by GetMonoPolFNFM
*****************************************************************************/
void GetMonoPolRefProb(double mc_pol_FN, double mc_pol_FM,
		       double mc_pol_sy, double *mc_pol_prob) {
  *mc_pol_prob = mc_pol_FN*mc_pol_FN + mc_pol_FM*mc_pol_FM
    + 2*mc_pol_FN*mc_pol_FM*mc_pol_sy;
  return;
}

/****************************************************************************
* void SetMonoPolRefOut(double FN, double FM, double refProb,
*		     double* sx, double* sy, double* sz) {
*
* ACTION: Set the outgoing polarisation vector of the reflected neutrons
* given FN, FM and the reflection probability.
*
* For a monochromator with up direction along y the outgoing polarisation
* is given as:
*	sx = (FN*FN - FM*FM)*sx_in/R0;
*	sy = ((FN*FN - FM*FM)*sy_in + 2*FN*FM + FM*FM*sy_in)/R0;
*	sz = (FN*FN - FM*FM)*sz_in/R0;
* where sx_in, sy_in, and sz_in is the incoming polarisation, and
* FN and FM are calculated from Rup and Rdown by GetMonoPolFNFM
*****************************************************************************/
void SetMonoPolRefOut(double mc_pol_FN, double mc_pol_FM,
		      double mc_pol_refProb, double* mc_pol_sx,
		      double* mc_pol_sy, double* mc_pol_sz) {
  *mc_pol_sx = (mc_pol_FN*mc_pol_FN - mc_pol_FM*mc_pol_FM)*(*mc_pol_sx)
    /mc_pol_refProb;
  *mc_pol_sy = ((mc_pol_FN*mc_pol_FN - mc_pol_FM*mc_pol_FM)*(*mc_pol_sy)
		+ 2*mc_pol_FN*mc_pol_FM + 2*mc_pol_FM*mc_pol_FM*(*mc_pol_sy))
    /mc_pol_refProb;
  *mc_pol_sz = (mc_pol_FN*mc_pol_FN - mc_pol_FM*mc_pol_FM)*(*mc_pol_sz)
    /mc_pol_refProb;
  return;
}

/****************************************************************************
* void SetMonoPolTransOut(double FN, double FM, double refProb,
*			  double* sx, double* sy, double* sz) {
*
* ACTION: Set the outgoing polarisation vector of the transmitted neutrons
* given FN, FM and the REFLECTION probability.
*
* We use that the polarization is conserved so:
* s_in = refProb*s_ref+(1-refProb)*s_trans, and then
* s_trans = (s_in-refProb*s_ref)/(1-refProb)
* where refProb is calculated using the routine GetMonoPolRefProb
* and s_ref is calculated by SetMonoPolRefOut
*****************************************************************************/
void SetMonoPolTransOut(double mc_pol_FN, double mc_pol_FM,
			double mc_pol_refProb, double* mc_pol_sx,
			double* mc_pol_sy, double* mc_pol_sz) {
  double mc_pol_sx_ref = *mc_pol_sx, mc_pol_sy_ref = *mc_pol_sy;
  double mc_pol_sz_ref = *mc_pol_sz;

  // By passing 1 as probability we get mc_pol_refProb*s_out_ref
  SetMonoPolRefOut(mc_pol_FN, mc_pol_FM, 1,
		   &mc_pol_sx_ref, &mc_pol_sy_ref, &mc_pol_sz_ref);
  *mc_pol_sx = (*mc_pol_sx - mc_pol_sx_ref)/(1 - mc_pol_refProb);
  *mc_pol_sy = (*mc_pol_sy - mc_pol_sy_ref)/(1 - mc_pol_refProb);
  *mc_pol_sz = (*mc_pol_sz - mc_pol_sz_ref)/(1 - mc_pol_refProb);
  return;
}

/****************************************************************************
* void SimpleNumMagnetPrecession(double x, double y, double z, double t,
*			         double vx, double vy, double vz,
*			         double* sx, double* sy, double* sz, double dt)
*
*****************************************************************************/
void SimpleNumMagnetPrecession(double mc_pol_x, double mc_pol_y,
			       double mc_pol_z, double mc_pol_time,
			       double mc_pol_vx, double mc_pol_vy,
			       double mc_pol_vz,
			       double* mc_pol_sx, double* mc_pol_sy,
			       double* mc_pol_sz, double mc_pol_deltaT,
			       Coords mc_pol_posLM, Rotation mc_pol_rotLM) {

  double mc_pol_Bx, mc_pol_By, mc_pol_Bz, mc_pol_phiz;
  double mc_pol_BxStart, mc_pol_ByStart, mc_pol_BzStart, mc_pol_Bstart;
  double mc_pol_BxTemp, mc_pol_ByTemp, mc_pol_BzTemp, mc_pol_Btemp;
  double mc_pol_Bstep, mc_pol_timeStep, mc_pol_sp;
  const double mc_pol_spThreshold  = cos(mc_pol_angular_accuracy);
  const double mc_pol_startTimeStep = mc_pol_initial_timestep; // s
  double dummy1, dummy2;
  Rotation mc_pol_rotBack;

  mcMagneticField=mcmagnet_get_field;

  //printf("pos_at_caller(xyz)( %g %g %g )\n", mc_pol_x,mc_pol_y,mc_pol_z);
  // change coordinates from current local system to lab system
  mccoordschanges(mc_pol_posLM, mc_pol_rotLM,
		 &mc_pol_x, &mc_pol_y, &mc_pol_z,
		 &mc_pol_vx, &mc_pol_vy, &mc_pol_vz, mc_pol_sx, mc_pol_sy, mc_pol_sz);
  //printf("pos_at_labaftertranformation(xyz)( %g %g %g )\n", mc_pol_x,mc_pol_y,mc_pol_z);

  // get initial B-field value
  mcMagneticField(mc_pol_x, mc_pol_y, mc_pol_z, mc_pol_time,
		  &mc_pol_BxTemp, &mc_pol_ByTemp, &mc_pol_BzTemp,NULL);

  do {

    mc_pol_Bx = 0; mc_pol_By = 0; mc_pol_Bz = 0; mc_pol_phiz = 0;
    mc_pol_BxStart = mc_pol_BxTemp; mc_pol_ByStart = mc_pol_ByTemp;
    mc_pol_BzStart = mc_pol_BzTemp;
    mc_pol_Bstart =
      sqrt(mc_pol_BxStart*mc_pol_BxStart + mc_pol_ByStart*mc_pol_ByStart
	   + mc_pol_BzStart*mc_pol_BzStart);
    mc_pol_timeStep = mc_pol_startTimeStep;

    if(mc_pol_deltaT<mc_pol_timeStep)
      mc_pol_timeStep = mc_pol_deltaT;

    do {

      mcMagneticField(mc_pol_x+mc_pol_vx*mc_pol_timeStep,
		      mc_pol_y+mc_pol_vy*mc_pol_timeStep,
		      mc_pol_z+mc_pol_vz*mc_pol_timeStep,
		      mc_pol_time+mc_pol_timeStep,
		      &mc_pol_BxTemp, &mc_pol_ByTemp, &mc_pol_BzTemp,NULL);
      // not so elegant, but this is how we make sure that the steps decrease
      // when the WHILE condition is not met
      mc_pol_timeStep *= 0.5;

      mc_pol_Btemp =
	sqrt(mc_pol_BxTemp*mc_pol_BxTemp + mc_pol_ByTemp*mc_pol_ByTemp
	     + mc_pol_BzTemp*mc_pol_BzTemp);

      mc_pol_sp =
	scalar_prod(mc_pol_BxStart, mc_pol_ByStart, mc_pol_BzStart,
		    mc_pol_BxTemp, mc_pol_ByTemp, mc_pol_BzTemp);
      mc_pol_sp /= mc_pol_Bstart*mc_pol_Btemp;

    } while (mc_pol_sp<mc_pol_spThreshold && mc_pol_timeStep>FLT_EPSILON);

    mc_pol_timeStep*=2;

    // update coordinate values
    mc_pol_x += mc_pol_vx*mc_pol_timeStep;
    mc_pol_y += mc_pol_vy*mc_pol_timeStep;
    mc_pol_z += mc_pol_vz*mc_pol_timeStep;
    mc_pol_time += mc_pol_timeStep;
    mc_pol_deltaT -= mc_pol_timeStep;

    mc_pol_Bx = 0.5 * (mc_pol_BxStart + mc_pol_BxTemp);
    mc_pol_By = 0.5 * (mc_pol_ByStart + mc_pol_ByTemp);
    mc_pol_Bz = 0.5 * (mc_pol_BzStart + mc_pol_BzTemp);
    mc_pol_phiz = fmod(sqrt(mc_pol_Bx*mc_pol_Bx+
			    mc_pol_By*mc_pol_By+
			    mc_pol_Bz*mc_pol_Bz)
		       *mc_pol_timeStep*mc_pol_omegaL, 2*PI);

    // Do the neutron spin precession

    if(!(mc_pol_Bx==0 && mc_pol_By==0 && mc_pol_Bz==0)) {

      double mc_pol_sx_in = *mc_pol_sx;
      double mc_pol_sy_in = *mc_pol_sy;
      double mc_pol_sz_in = *mc_pol_sz;

      rotate(*mc_pol_sx, *mc_pol_sy, *mc_pol_sz,
	     mc_pol_sx_in, mc_pol_sy_in, mc_pol_sz_in,
	     mc_pol_phiz, mc_pol_Bx, mc_pol_By, mc_pol_Bz);
    }

  } while (mc_pol_deltaT>0);

  // change back spin coordinates from lab system to local system
  rot_transpose(mc_pol_rotLM, mc_pol_rotBack);
  mccoordschange_polarisation(mc_pol_rotBack, mc_pol_sx, mc_pol_sy, mc_pol_sz);

}

/****************************************************************************
* double GetConstantField(double length, double lambda, double angle)
*
* Return the magnetic field in Tesla required to flip a neutron with
* wavelength lambda(1/velocity), angle degrees, over the specified
* length(=time*velocity).
*
*****************************************************************************/
double GetConstantField(double mc_pol_length, double mc_pol_lambda,
			double mc_pol_angle)
{
  const double mc_pol_velocity = K2V*2*PI/mc_pol_lambda;
  const double mc_pol_time = mc_pol_length/mc_pol_velocity;

  // B*omegaL*time = angle
  return mc_pol_angle*DEG2RAD/mc_pol_omegaL/mc_pol_time; // T
}

/* end of regular pol-lib.c */


/* Shared user declarations for all components types 'Monochromator_curved'. */
#ifndef GAUSS
/* Define these arrays only once for all instances. */
/* Values for Gauss quadrature. Taken from Brice Carnahan, H. A. Luther and
James O Wilkes, "Applied numerical methods", Wiley, 1969, page 103.
This reference is available from the Copenhagen UB2 library */
double Gauss_X[] = {-0.987992518020485, -0.937273392400706, -0.848206583410427,
-0.724417731360170, -0.570972172608539, -0.394151347077563,
-0.201194093997435, 0, 0.201194093997435,
0.394151347077563, 0.570972172608539, 0.724417731360170,
0.848206583410427, 0.937273392400706, 0.987992518020485};
double Gauss_W[] = {0.030753241996117, 0.070366047488108, 0.107159220467172,
0.139570677926154, 0.166269205816994, 0.186161000115562,
0.198431485327111, 0.202578241925561, 0.198431485327111,
0.186161000115562, 0.166269205816994, 0.139570677926154,
0.107159220467172, 0.070366047488108, 0.030753241996117};


#define GAUSS(x,mean,rms) \
  (exp(-((x)-(mean))*((x)-(mean))/(2*(rms)*(rms)))/(sqrt(2*PI)*(rms)))
#endif



/* Shared user declarations for all components types 'Monochromator_2foc'. */

#ifndef DIV_CUTOFF
#define DIV_CUTOFF 2            /* ~ 10^-5 cutoff. */
#endif

/* Shared user declarations for all components types 'Single_crystal'. */
/* used for reading data table from file */

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


//////////////////////////////////////////////////////////////////////////////
//OpenCL lilbrary file for McStas
//
// global variables assigned to the OpenCL mechanism
// only initialised once
//Created by Jinyan LIU
//Date july 2014
//////////////////////////////////////////////////////////////////////////////



#ifdef USE_OPENCL
#ifndef OPENCL_LIB_H
#define OPENCL_LIB_H

#include <>

// this structure holds the full GPU context initialised from a Kernel source file
// to get if use: s = oclInitKernel(filename, 1);
#define MAX_GPU_COUNT 8

struct opencl_context {
  cl_context       GPUContext;                     // OpenCL context
  cl_command_queue CommandQueue[MAX_GPU_COUNT];      // OpenCL command queue
  cl_kernel        Kernel;
};

struct opencl_context oclInitKernel(char *filename, cl_uint nDevice);

// Mersenne OCL implementation
#define   MT_RNG_COUNT 4096
#define          MT_MM 397            
#define          MT_NN 624               
#define       MT_WMASK 0xFFFFFFFFU
#define       MT_UMASK 0x80000000U      
#define       MT_LMASK 0x7fffffffU       
#define      MT_SHIFT0 12
#define      MT_SHIFTB 7
#define      MT_SHIFTC 15
#define      MT_SHIFT1 18
#define      DCMT_SEED 4172

////////////////////////////////////////////////////////////////////////////////
//OpenCL Global Variable 
//////////////////////////////////////////////////////////////////////////////// 
struct opencl_context oclContext_mt;
unsigned int        **oclContext_mt_buffer=NULL;      // buffer holding a list of random numbers, on the CPU
cl_mem               *oclContext_mt_buffer_GPU=NULL;  // buffer holding a list of random numbers, on the GPU/OpenCL
int                   oclContext_mt_counter = -2;     // index of random number read from buffer

void         mt_srandom_opencl(unsigned long s);
unsigned int mt_random_opencl(void);

#endif
#endif

//////////////////////////////////////////////////////////////////////////////
//OpenCL lilbrary file for McStas
//
// global variables assigned to the OpenCL mechanism
// only initialised once
//Created by Jinyan LIU
//Date july 2014
//////////////////////////////////////////////////////////////////////////////

#ifdef USE_OPENCL
#ifndef OPENCL_LIB_C
#define OPENCL_LIB_C

//////////////////////////////////////////////////////////////////////////////
//! fileparts: return the name of the file between '/' and '.'
//!
//! @return the name of the file if succeeded, 0 otherwise
//! @param name      filename maybe with '/' , '\' or '.'
//! 
//////////////////////////////////////////////////////////////////////////////

#if defined(WIN32) || defined(_WIN32) || defined(_WIN64)
#define PATHSEP_C '\\'
#define PATHSEP_S "\\"
#else  /* !WIN32 */
#define PATHSEP_C '/'
#define PATHSEP_S "/"
#endif /* !WIN32 */

char* fileparts(char *name)
{
  char *Name=NULL;

  if (name && strlen(name)) {
    char *dot_pos    = NULL;
    char *path_pos   = NULL;
    char *end_pos    = NULL;
    char *name_pos   = NULL;
    size_t  name_length= 0;

    end_pos = name+strlen(name);  /* end of file name */

    /* extract path: searches for last file separator */
    path_pos= strrchr(name, PATHSEP_C);  /* last PATHSEP */

    if (!path_pos) {
      path_pos   =name;
      name_pos   =name;
    } else {
      name_pos    = path_pos+1; /* from start to path+sep */
    }

    /* extract ext: now looks for the 'dot' */
    dot_pos = strrchr(name_pos, '.');           /* last dot */
    if (dot_pos <= name_pos)
      dot_pos = end_pos;

    /* extract Name (without extension) */
    name_length = dot_pos - name_pos; /* from path to dot */
    if (name_length) {
      Name = (char*)malloc(name_length);
      if (Name) strncpy(Name, name_pos, name_length);
      Name[name_length]='\0';
    }
  } /* if (name) */
  return (Name);

} /* fileparts */


//////////////////////////////////////////////////////////////////////////////
//! oclLoadProgSource: Loads a Program file and prepends the cPreamble to the code.
//!
//! @return the source string if succeeded, 0 otherwise
//! @param cFilename        program filename
//! @param cPreamble        code that is prepended to the loaded file, typically a set of #defines or a header
//! @param szFinalLength    returned length of the code string
//////////////////////////////////////////////////////////////////////////////

char* oclLoadProgSource(const char* cFilename, const char* cPreamble, size_t* szFinalLength)
{
    // locals 
    FILE* pFileStream = NULL;
    size_t szSourceLength;

    // open the OpenCL source code file
    #ifdef _WIN32   // Windows version
        if(fopen_s(&pFileStream, cFilename, "rb") != 0) 
        {       
            return NULL;
        }
    #else           // Linux version
        pFileStream = fopen(cFilename, "rb");
        if(pFileStream == 0) 
        {       
            return NULL;
        }
    #endif

    size_t szPreambleLength = strlen(cPreamble);

    // get the length of the source code
    fseek(pFileStream, 0, SEEK_END); 
    szSourceLength = ftell(pFileStream);
    fseek(pFileStream, 0, SEEK_SET); 

    // allocate a buffer for the source code string and read it in
    char* cSourceString = (char *)malloc(szSourceLength + szPreambleLength + 1); 
    memcpy(cSourceString, cPreamble, szPreambleLength);
    if (fread((cSourceString) + szPreambleLength, szSourceLength, 1, pFileStream) != 1)
    {
        fclose(pFileStream);
        free(cSourceString);
        return NULL;
    }

    // close the file and return the total length of the combined (preamble + source) string
    fclose(pFileStream);
    if(szFinalLength != 0)
    {
        *szFinalLength = szSourceLength + szPreambleLength;
    }
    cSourceString[szSourceLength + szPreambleLength] = '\0';

    return cSourceString;
} /* end oclLoadProgSource */

//////////////////////////////////////////////////////////////////////////////
//! oclGetPlatformID: Gets the platform ID for NVIDIA if available, otherwise default
//!
//! @return the id 
//! @param clSelectedPlatformID         OpenCL plateform ID
//////////////////////////////////////////////////////////////////////////////

cl_int oclGetPlatformID(cl_platform_id* clSelectedPlatformID)
{
    char chBuffer[1024];
    cl_uint num_platforms,i; 
    cl_platform_id* clPlatformIDs;
    cl_int ciErrNum;
    *clSelectedPlatformID = NULL;

    // Get OpenCL platform count
    ciErrNum = clGetPlatformIDs (0, NULL, &num_platforms);
    if (ciErrNum != CL_SUCCESS)
    {
        printf(" Error %i in clGetPlatformIDs Call !!!\n\n", ciErrNum);
        return -1000;
    }
    else 
    {
        if(num_platforms == 0)
        {
            printf("No OpenCL platform found!\n\n");
            return -2000;
        }
        else 
        {
            // if there's a platform or more, make space for ID's
            if ((clPlatformIDs = (cl_platform_id*)malloc(num_platforms * sizeof(cl_platform_id))) == NULL)
            {
                printf("Failed to allocate memory for cl_platform ID's!\n\n");
                return -3000;
            }

            // get platform info for each platform and trap the NVIDIA platform if found
            ciErrNum = clGetPlatformIDs (num_platforms, clPlatformIDs, NULL);
	    printf("We found %i CL platforms, selecting the first of them:\n",num_platforms);
	    
            ciErrNum = clGetPlatformInfo (clPlatformIDs[0], CL_PLATFORM_NAME, 1024, &chBuffer, NULL);
            if(ciErrNum == CL_SUCCESS)
            {
		printf("Platform: %s\n",chBuffer);
                *clSelectedPlatformID = clPlatformIDs[0];
            }

            free(clPlatformIDs);
        }
    }

    return CL_SUCCESS;
} /* end oclGetPlatformID */

//////////////////////////////////////////////////////////////////////////////
//! oclInit: initialize 'nDevice' GPU's
//!
//! @return error code or CL_SUCCESS
//! @param nDevice         nb of devices requested OpenCL
//
// This function should initialize cxGPUContext and cqCommandQueue
//////////////////////////////////////////////////////////////////////////////


struct opencl_context oclInitKernel(char *filename, cl_uint nDevice) {

  cl_program cpProgram;                           // OpenCL program
  
  cl_platform_id cpPlatform; 

  cl_device_id* cdDevices;                        // OpenCL device list   
  
  cl_int ciErr1, ciErr2;                          // Error code var
  
  char *basename = NULL;
  const int nPerRng = 5860;                       // # of recurrence steps, must be even if do Box-Muller transformation
  const int nRand = MT_RNG_COUNT * nPerRng;       // Output size  

  cl_uint i;
  
  struct opencl_context oclContext;
  oclContext.GPUContext = NULL;
  oclContext.Kernel     = NULL;
  
  ciErr1 = oclGetPlatformID(&cpPlatform);
  if(ciErr1 != CL_SUCCESS) 
    return(oclContext);
  ciErr1 = clGetDeviceIDs(NULL, CL_DEVICE_TYPE_GPU, 0, NULL, &nDevice);

  cdDevices = (cl_device_id *)malloc(nDevice * sizeof(cl_device_id) );
  ciErr1 =clGetDeviceIDs(NULL, CL_DEVICE_TYPE_GPU,  nDevice, cdDevices, NULL);
  
  
  oclContext.GPUContext = clCreateContext(0, nDevice, cdDevices, NULL, NULL, &ciErr1);
  if(ciErr1 != CL_SUCCESS) 
    return(oclContext);
        
  for (i = 0; i < nDevice; i++) 
  {
      oclContext.CommandQueue[i] = clCreateCommandQueue(oclContext.GPUContext, cdDevices[i], 0, &ciErr1);
  }
  
  if (!oclContext.GPUContext) return (oclContext);

  /* load mersenne twister kernel source from disk */
  size_t szKernelLength; // Byte size of kernel code
  char   *cKernel = oclLoadProgSource(filename, "// My comment\n", &szKernelLength);
  
  if(cKernel == NULL) {
    printf("%s:%u: OpenCL error: Failed to open the OpenCL program source file %s\n", __FILE__, __LINE__, filename);
    return(oclContext);
  }
  
  cpProgram = clCreateProgramWithSource(oclContext.GPUContext, 1, (const char **)&cKernel, &szKernelLength, &ciErr1);
  if(ciErr1 != CL_SUCCESS) 
    return(oclContext);  
       
  if (cpProgram== NULL) {
    printf("%s:%u: OpenCL error: Failed to create %s OpenCL program\n", __FILE__, __LINE__, filename);
    return(oclContext);
  }  
       
  ciErr1 |= clBuildProgram(cpProgram, 0, NULL, "-cl-fast-relaxed-math", NULL, NULL);
  if (ciErr1 !=CL_SUCCESS){
    printf("%s:%u: OpenCL error: Failed to build %s kernel\n", __FILE__, __LINE__, filename);
    for (i=0; i < nDevice; i++) {
      // Determine the size of the log
      size_t len;
      clGetProgramBuildInfo(cpProgram, cdDevices[i], CL_PROGRAM_BUILD_LOG, 0, NULL, &len); 
      // Allocate memory for the log
      char *build_log = (char*)malloc(len);
      // Get the log
      ciErr1 |= clGetProgramBuildInfo(cpProgram, cdDevices[i], CL_PROGRAM_BUILD_LOG, len, build_log, NULL);
      // Print the log
      printf("-------Build log-----------\n");
      printf("\nBuildInfo:%s\n", build_log);
      free(build_log);
      clFinish(oclContext.CommandQueue[i]);
    }
    return(oclContext);
  }
  
  basename = fileparts(filename);
  
  oclContext.Kernel = clCreateKernel(cpProgram, basename, &ciErr1);
  
  free(basename);

  if (oclContext.Kernel == NULL || ciErr1 !=CL_SUCCESS) {
    printf("%s:%u: OpenCL error: Failed to create %s kernel\n", __FILE__, __LINE__, filename);
    return(oclContext);
  }

  for (i = 0; i < nDevice; i++)
    {
      clFinish(oclContext.CommandQueue[i]);
    }
    
  clReleaseProgram(cpProgram); 
  free(cKernel);
  free(cdDevices);
  
  return oclContext;
}

/* ========================================================================== */
/*        OpenCL implementation of random MersenneTwister generator           */
/* ========================================================================== */



unsigned int **mt_random_opencl_fill_buffer(cl_uint nDevice, struct opencl_context oclContext)
{

  cl_int ciErr1=0, ciErr2=0;                      // Error code var
  size_t globalWorkSize[1] = {MT_RNG_COUNT};      // 1D var for Total # of work items
  size_t localWorkSize[1] = {128};                // 1D var for # of work items in the work group

  const int nPerRng =5860;                        // # of recurrence steps, must be even if do Box-Muller transformation
  const int nRand = MT_RNG_COUNT * nPerRng;       // Output size   
  char      allocate_memory=0;
  cl_uint   iDevice,i;

  printf("Filling buffer\n");

  /* this section launches the kernel to fill the buffer (when counter = -1) */
  
  //  Using GPU(s)...
  if (oclContext.Kernel == NULL || oclContext.GPUContext==NULL)
    return(NULL);

  //  Initialization: load MT parameters and init host buffers
  if (!oclContext_mt_buffer) {
    // first call to this function: allocate memory
    allocate_memory = 1;
    oclContext_mt_buffer = (unsigned int**)malloc(nDevice*sizeof(unsigned int*));
    if (!oclContext_mt_buffer) return NULL;
    for (iDevice = 0; iDevice < nDevice; oclContext_mt_buffer[iDevice++]=NULL);
  }

  //  Allocate memory
  if (allocate_memory)
    for (iDevice = 0; iDevice < nDevice; iDevice++)
    {
       if (!oclContext_mt_buffer[iDevice])
           oclContext_mt_buffer[iDevice]  = (unsigned int*)malloc(sizeof(unsigned int)*nRand); // Host buffers for GPU output
       if (!oclContext_mt_buffer[iDevice]) return(NULL);
    }

  if (allocate_memory && !oclContext_mt_buffer_GPU) {
    oclContext_mt_buffer_GPU = (cl_mem*)malloc(nDevice*sizeof(cl_mem));
    if (!oclContext_mt_buffer_GPU) return NULL;
    for (iDevice = 0; iDevice < nDevice; oclContext_mt_buffer_GPU[iDevice++]=NULL);
  }
  
  if (allocate_memory)
    for (iDevice = 0; iDevice < nDevice; iDevice++)
    {
        if (!oclContext_mt_buffer_GPU[iDevice]) {
          oclContext_mt_buffer_GPU[iDevice] = clCreateBuffer(oclContext.GPUContext, CL_MEM_READ_WRITE, sizeof(cl_uint)*nRand, NULL, &ciErr2);
        
          if (ciErr2 !=CL_SUCCESS) return(NULL);
        }
        
    } /* for */
  
  for (iDevice = 0; iDevice < nDevice; iDevice++)
  {
      clFinish(oclContext.CommandQueue[iDevice]);
  }
  
  for (iDevice = 0; iDevice < nDevice; iDevice++)
  {
    time_t t;
    int seed=(int)time(&t);
    if (allocate_memory) {
      ciErr1 |= clSetKernelArg(oclContext.Kernel, 0, sizeof(cl_mem), (void*)&oclContext_mt_buffer_GPU[iDevice]);
      if (ciErr1 !=CL_SUCCESS) return(NULL);
        
      ciErr1 |= clSetKernelArg(oclContext.Kernel, 1, sizeof(int),    (void*)&nPerRng);
      if (ciErr1 !=CL_SUCCESS) return(NULL);
    }

    // Each call to fill buffer must use a different seed
    ciErr1 |= clSetKernelArg(oclContext.Kernel, 2, sizeof(int),    (void*)&seed);
      if (ciErr1 !=CL_SUCCESS) return(NULL);
      
    // call the Kernel with the queue
    ciErr1 |= clEnqueueNDRangeKernel(oclContext.CommandQueue[iDevice], oclContext.Kernel, 1, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);

    if (ciErr1 !=CL_SUCCESS) return (NULL);
  }
    
  for (iDevice = 0; iDevice < nDevice; iDevice++)
  {
      clFinish(oclContext.CommandQueue[iDevice]);
  }

  //  Read back results
  
  for (iDevice = 0; iDevice < nDevice; iDevice++)
  {
      ciErr1 = clEnqueueReadBuffer(oclContext.CommandQueue[iDevice], oclContext_mt_buffer_GPU[iDevice], CL_TRUE, 0, 
          sizeof(cl_uint) * nRand, oclContext_mt_buffer[iDevice], 0, NULL, NULL);
  
      if (ciErr1 !=CL_SUCCESS) return(NULL);
  }  

  return (unsigned int**) oclContext_mt_buffer;
} /* end mt_random_opencl_fill_buffer */



// main Mersenne random call ---------------------------------------------------
void mt_srandom_opencl(unsigned long s) {
  fprintf(stderr, "%s: The GPU implementation does not support manual seed setting. Using seed=clock().\n", __FILE__);
}

unsigned int mt_random_opencl(void) // Should be called by others 
{
  unsigned int nDevice = 1;                              // Use 1st opencl Device, > 0
  const int    nPerRng = 5860;                           // # of recurrence steps, must be even if do Box-Muller transformation
  const int    nRand = MT_RNG_COUNT * nPerRng;           // Output size
 
  if (oclContext_mt_counter >= nRand) {    
    printf("%s:%i: INFO: reached max (oclContext_mt_counter=%i)\n", __FILE__, __LINE__, oclContext_mt_counter); fflush(NULL);
    oclContext_mt_counter = -1;
  }  
  
  if (oclContext_mt_counter <= -2) {
    /* the first time we create the OCL kernel and the buffer */
    oclContext_mt = oclInitKernel("MersenneTwister.cl", nDevice);
    oclContext_mt_counter = -1;
  }
  
  // in case the OpenCL context has noot been initialized, default to serial MT
  if (oclContext_mt.Kernel == NULL) {
    return(mt_random());
  }
    
  if (oclContext_mt_counter <= -1 /*|| oclContext_mt_counter >= nRand*/ ) {
    /* when the buffer is just created or fully read, we re-fill the buffer */
    if (mt_random_opencl_fill_buffer(nDevice, oclContext_mt) == NULL) {
      fprintf(stderr, "%s: Could not find the OpenCL MT buffer. Using serial MT.\n", __FILE__);
      return(mt_random());
    }

    oclContext_mt_counter = 0;
  } 
  
  /* return the value of the random number read from the buffer */
  printf("Rng returned %u as %u\n",oclContext_mt_counter++,oclContext_mt_buffer[nDevice-1][oclContext_mt_counter++]);
  return (oclContext_mt_buffer[nDevice-1][oclContext_mt_counter++]);
} /* mt_random_opencl */

#endif
#endif


/* Declare structures and functions only once in each instrument. */
#ifndef SINGLE_CRYSTAL_DECL
#define SINGLE_CRYSTAL_DECL

#ifndef Mosaic_AB_Undefined
#define Mosaic_AB_Undefined {0,0, 0,0,0, 0,0,0}
#endif

struct hkl_data
{
int h,k,l;                  /* Indices for this reflection */
double F2;                  /* Value of structure factor */
double tau_x, tau_y, tau_z; /* Coordinates in reciprocal space */
double tau;                 /* Length of (tau_x, tau_y, tau_z) */
      double u1x, u1y, u1z;       /* First axis of local coordinate system */
      double u2x, u2y, u2z;       /* Second axis of local coordinate system */
      double u3x, u3y, u3z;       /* Third axis of local coordinate system */
      double sig123;              /* The product sig1*sig2*sig3 = volume of spot */
      double m1, m2, m3;          /* Diagonal matrix representation of Gauss */
      double cutoff;              /* Cutoff value for Gaussian tails */
    };

  struct tau_data
    {
      int index;                  /* Index into reflection table */
      double refl;
      double xsect;
      /* The following vectors are in local koordinates. */
      double rho_x, rho_y, rho_z; /* The vector ki - tau */
      double rho;                 /* Length of rho vector */
      double ox, oy, oz;          /* Origin of Ewald sphere tangent plane */
      double b1x, b1y, b1z;       /* Spanning vectors of Ewald sphere tangent */
      double b2x, b2y, b2z;
      double l11, l12, l22;       /* Cholesky decomposition L of 2D Gauss */
      double y0x, y0y;            /* 2D Gauss center in tangent plane */
    };

  struct hkl_info_struct
    {
      struct hkl_data *list;      /* Reflection array */
      int count;                  /* Number of reflections */
      struct tau_data *tau_list;  /* Reflections close to Ewald Sphere */
      double m_delta_d_d;         /* Delta-d/d FWHM */
      double m_ax,m_ay,m_az;      /* First unit cell axis (direct space, AA) */
      double m_bx,m_by,m_bz;      /* Second unit cell axis */
      double m_cx,m_cy,m_cz;      /* Third unit cell axis */
      double asx,asy,asz;         /* First reciprocal lattice axis (1/AA) */
      double bsx,bsy,bsz;         /* Second reciprocal lattice axis */
      double csx,csy,csz;         /* Third reciprocal lattice axis */
      double m_a, m_b, m_c;       /* length of lattice parameter lengths */
      double m_aa, m_bb, m_cc;    /* lattice angles */
      double sigma_a, sigma_i;    /* abs and inc X sect */
      double rho;                 /* density */
      double at_weight;           /* atomic weight */
      double at_nb;               /* nb of atoms in a cell */
      double V0;                  /* Unit cell volume (AA**3) */
      int    column_order[5];     /* column signification [h,k,l,F,F2] */
      int    recip;               /* Flag to indicate if recip or direct cell axes given */
      int    shape;               /* 0:cylinder, 1:box, 2:sphere 3:any shape*/
      int    flag_warning;        /* number of warnings */
      char   type;                /* type of last event: t=transmit,c=coherent or i=incoherent */
      int    h,k,l;               /* last coherent scattering momentum transfer indices */
      int    tau_count;           /* Number of reflections within cutoff */
      double coh_refl, coh_xsect; /* cross section computed with last tau_list */
      double kix, kiy, kiz;       /* last incoming neutron ki */
      int    nb_reuses, nb_refl, nb_refl_count;
    };

  int SX_list_compare (void const *a, void const *b)
  {
     struct hkl_data const *pa = a;
     struct hkl_data const *pb = b;
     double s = pa->tau - pb->tau;

     if (!s) return 0;
     else    return (s < 0 ? -1 : 1);
  } /* PN_list_compare */

  /* ------------------------------------------------------------------------ */
  int
  read_hkl_data(char *SC_file, struct hkl_info_struct *info,
      double SC_mosaic, double SC_mosaic_a, double SC_mosaic_b, double SC_mosaic_c, double *SC_mosaic_AB)
  {
    struct hkl_data *list = NULL;
    int size = 0;
    t_Table sTable; /* sample data table structure from SC_file */
    int i=0;
    double tmp_x, tmp_y, tmp_z;
    char **parsing;
    char flag=0;
    double nb_atoms=1;

    if (!SC_file || !strlen(SC_file) || !strcmp(SC_file,"NULL") || !strcmp(SC_file,"0")) {
      info->count = 0;
      flag=1;
    }
    if (!flag) {
      Table_Read(&sTable, SC_file, 1); /* read 1st block data from SC_file into sTable*/
      if (sTable.columns < 4) {
        fprintf(stderr, "Single_crystal: Error: The number of columns in %s should be at least %d for [h,k,l,F2]\n", SC_file, 4);
        return(0);
      }
      if (!sTable.rows) {
        fprintf(stderr, "Single_crystal: Error: The number of rows in %s should be at least %d\n", SC_file, 1);
        return(0);
      } else size = sTable.rows;

      /* parsing of header */
      parsing = Table_ParseHeader(sTable.header,
        "sigma_abs","sigma_a ",
        "sigma_inc","sigma_i ",
        "column_h",
        "column_k",
        "column_l",
        "column_F ",
        "column_F2",
        "Delta_d/d",
        "lattice_a ",
        "lattice_b ",
        "lattice_c ",
        "lattice_aa",
        "lattice_bb",
        "lattice_cc",
        "nb_atoms","multiplicity",
        NULL);

      if (parsing) {
        if (parsing[0] && !info->sigma_a) info->sigma_a=atof(parsing[0]);
        if (parsing[1] && !info->sigma_a) info->sigma_a=atof(parsing[1]);
        if (parsing[2] && !info->sigma_i) info->sigma_i=atof(parsing[2]);
        if (parsing[3] && !info->sigma_i) info->sigma_i=atof(parsing[3]);
        if (parsing[4])                   info->column_order[0]=atoi(parsing[4]);
        if (parsing[5])                   info->column_order[1]=atoi(parsing[5]);
        if (parsing[6])                   info->column_order[2]=atoi(parsing[6]);
        if (parsing[7])                   info->column_order[3]=atoi(parsing[7]);
        if (parsing[8])                   info->column_order[4]=atoi(parsing[8]);
        if (parsing[9] && info->m_delta_d_d <0) info->m_delta_d_d=atof(parsing[9]);
        if (parsing[10] && !info->m_a)    info->m_a =atof(parsing[10]);
        if (parsing[11] && !info->m_b)    info->m_b =atof(parsing[11]);
        if (parsing[12] && !info->m_c)    info->m_c =atof(parsing[12]);
        if (parsing[13] && !info->m_aa)   info->m_aa=atof(parsing[13]);
        if (parsing[14] && !info->m_bb)   info->m_bb=atof(parsing[14]);
        if (parsing[15] && !info->m_cc)   info->m_cc=atof(parsing[15]);
        if (parsing[16])   nb_atoms=atof(parsing[16]);
        if (parsing[17])   nb_atoms=atof(parsing[17]);
        for (i=0; i<=17; i++) if (parsing[i]) free(parsing[i]);
        free(parsing);
      }
    }

    if (nb_atoms > 1) { info->sigma_a *= nb_atoms; info->sigma_i *= nb_atoms; }

    /* special cases for the structure definition */
    if (info->m_ax || info->m_ay || info->m_az) info->m_a=0; /* means we specify by hand the vectors */
    if (info->m_bx || info->m_by || info->m_bz) info->m_b=0;
    if (info->m_cx || info->m_cy || info->m_cz) info->m_c=0;

    /* compute the norm from vector a if missing */
    if (info->m_ax || info->m_ay || info->m_az) {
      double as=sqrt(info->m_ax*info->m_ax+info->m_ay*info->m_ay+info->m_az*info->m_az);
      if (!info->m_bx && !info->m_by && !info->m_bz) info->m_a=info->m_b=as;
      if (!info->m_cx && !info->m_cy && !info->m_cz) info->m_a=info->m_c=as;
    }
    if (info->m_a && !info->m_b) info->m_b=info->m_a;
    if (info->m_b && !info->m_c) info->m_c=info->m_b;

    /* compute the lattive angles if not set from data file. Not used when in vector mode. */
    if (info->m_a && !info->m_aa) info->m_aa=90;
    if (info->m_aa && !info->m_bb) info->m_bb=info->m_aa;
    if (info->m_bb && !info->m_cc) info->m_cc=info->m_bb;

    /* parameters consistency checks */
    if (!info->m_ax && !info->m_ay && !info->m_az && !info->m_a) {
      fprintf(stderr,
              "Single_crystal: Error: Wrong a lattice vector definition\n");
      return(0);
    }
    if (!info->m_bx && !info->m_by && !info->m_bz && !info->m_b) {
      fprintf(stderr,
              "Single_crystal: Error: Wrong b lattice vector definition\n");
      return(0);
    }
    if (!info->m_cx && !info->m_cy && !info->m_cz && !info->m_c) {
      fprintf(stderr,
              "Single_crystal: Error: Wrong c lattice vector definition\n");
      return(0);
    }
    if (info->m_aa && info->m_bb && info->m_cc && info->recip) {
      fprintf(stderr,
              "Single_crystal: Error: Selecting reciprocal cell and angles is unmeaningful\n");
      return(0);
    }

    /* when lengths a,b,c + angles are given (instead of vectors a,b,c) */
    if (info->m_aa && info->m_bb && info->m_cc)
    {
      double as,bs,cs;
      if (info->m_a) as = info->m_a;
      else as = sqrt(info->m_ax*info->m_ax+info->m_ay*info->m_ay+info->m_az*info->m_az);
      if (info->m_b) bs = info->m_b;
      else bs = sqrt(info->m_bx*info->m_bx+info->m_by*info->m_by+info->m_bz*info->m_bz);
      if (info->m_c) cs = info->m_c;
      else cs =  sqrt(info->m_cx*info->m_cx+info->m_cy*info->m_cy+info->m_cz*info->m_cz);

      info->m_bz = as; info->m_by = 0; info->m_bx = 0;
      info->m_az = bs*cos(info->m_cc*DEG2RAD);
      info->m_ay = bs*sin(info->m_cc*DEG2RAD);
      info->m_ax = 0;
      info->m_cz = cs*cos(info->m_bb*DEG2RAD);
      info->m_cy = cs*(cos(info->m_aa*DEG2RAD)-cos(info->m_cc*DEG2RAD)*cos(info->m_bb*DEG2RAD))
                     /sin(info->m_cc*DEG2RAD);
      info->m_cx = sqrt(cs*cs - info->m_cz*info->m_cz - info->m_cy*info->m_cy);

      printf("Single_crystal: %s structure a=%g b=%g c=%g aa=%g bb=%g cc=%g ",
        (flag ? "INC" : SC_file), as, bs, cs, info->m_aa, info->m_bb, info->m_cc);
    } else {
      if (!info->recip) {
        printf("Single_crystal: %s structure a=[%g,%g,%g] b=[%g,%g,%g] c=[%g,%g,%g] ",
	       (flag ? "INC" : SC_file), info->m_ax ,info->m_ay ,info->m_az,
	       info->m_bx ,info->m_by ,info->m_bz,
	       info->m_cx ,info->m_cy ,info->m_cz);
      } else {
        printf("Single_crystal: %s structure a*=[%g,%g,%g] b*=[%g,%g,%g] c*=[%g,%g,%g] ",
	       (flag ? "INC" : SC_file), info->m_ax ,info->m_ay ,info->m_az,
	       info->m_bx ,info->m_by ,info->m_bz,
	       info->m_cx ,info->m_cy ,info->m_cz);
      }
    }
    /* Compute reciprocal or direct lattice vectors. */
    if (!info->recip) {
      vec_prod(tmp_x, tmp_y, tmp_z,
	       info->m_bx, info->m_by, info->m_bz,
	       info->m_cx, info->m_cy, info->m_cz);
      info->V0 = fabs(scalar_prod(info->m_ax, info->m_ay, info->m_az, tmp_x, tmp_y, tmp_z));
      printf("V0=%g\n", info->V0);

      info->asx = 2*PI/info->V0*tmp_x;
      info->asy = 2*PI/info->V0*tmp_y;
      info->asz = 2*PI/info->V0*tmp_z;
      vec_prod(tmp_x, tmp_y, tmp_z, info->m_cx, info->m_cy, info->m_cz, info->m_ax, info->m_ay, info->m_az);
      info->bsx = 2*PI/info->V0*tmp_x;
      info->bsy = 2*PI/info->V0*tmp_y;
      info->bsz = 2*PI/info->V0*tmp_z;
      vec_prod(tmp_x, tmp_y, tmp_z, info->m_ax, info->m_ay, info->m_az, info->m_bx, info->m_by, info->m_bz);
      info->csx = 2*PI/info->V0*tmp_x;
      info->csy = 2*PI/info->V0*tmp_y;
      info->csz = 2*PI/info->V0*tmp_z;
    } else {
      info->asx = info->m_ax;
      info->asy = info->m_ay;
      info->asz = info->m_az;
      info->bsx = info->m_bx;
      info->bsy = info->m_by;
      info->bsz = info->m_bz;
      info->csx = info->m_cx;
      info->csy = info->m_cy;
      info->csz = info->m_cz;

      vec_prod(tmp_x, tmp_y, tmp_z,
	       info->bsx/(2*PI), info->bsy/(2*PI), info->bsz/(2*PI),
	       info->csx/(2*PI), info->csy/(2*PI), info->csz/(2*PI));
      info->V0 = 1/fabs(scalar_prod(info->asx/(2*PI), info->asy/(2*PI), info->asz/(2*PI), tmp_x, tmp_y, tmp_z));
      printf("V0=%g\n", info->V0);

      /*compute the direct cell parameters, ofr completeness*/
      info->m_ax = tmp_x*info->V0;
      info->m_ay = tmp_y*info->V0;
      info->m_az = tmp_z*info->V0;
      vec_prod(tmp_x, tmp_y, tmp_z,info->csx/(2*PI), info->csy/(2*PI), info->csz/(2*PI),info->asx/(2*PI), info->asy/(2*PI), info->asz/(2*PI));
      info->m_bx = tmp_x*info->V0;
      info->m_by = tmp_y*info->V0;
      info->m_bz = tmp_z*info->V0;
      vec_prod(tmp_x, tmp_y, tmp_z,info->asx/(2*PI), info->asy/(2*PI), info->asz/(2*PI),info->bsx/(2*PI), info->bsy/(2*PI), info->bsz/(2*PI));
      info->m_cx = tmp_x*info->V0;
      info->m_cy = tmp_y*info->V0;
      info->m_cz = tmp_z*info->V0;
    }

    if (flag) return(-1);

    if (!info->column_order[0] || !info->column_order[1] || !info->column_order[2]) {
      fprintf(stderr,
              "Single_crystal: Error: Wrong h,k,l column definition\n");
      return(0);
    }
    if (!info->column_order[3] && !info->column_order[4]) {
      fprintf(stderr,
              "Single_crystal: Error: Wrong F,F2 column definition\n");
      return(0);
    }

    /* allocate hkl_data array */
    list = (struct hkl_data*)malloc(size*sizeof(struct hkl_data));

    for (i=0; i<size; i++)
    {
      double h=0, k=0, l=0, F2=0;
      double b1[3], b2[3];
      double sig1, sig2, sig3;

      /* get data from table */
      h = Table_Index(sTable, i, info->column_order[0]-1);
      k = Table_Index(sTable, i, info->column_order[1]-1);
      l = Table_Index(sTable, i, info->column_order[2]-1);
      if (info->column_order[3])
      { F2= Table_Index(sTable, i, info->column_order[3]-1); F2 *= F2; }
      else if (info->column_order[4])
        F2= Table_Index(sTable, i, info->column_order[4]-1);

      list[i].h = h;
      list[i].k = k;
      list[i].l = l;
      list[i].F2 = F2;

      /* Precompute some values */
      list[i].tau_x = h*info->asx + k*info->bsx + l*info->csx;
      list[i].tau_y = h*info->asy + k*info->bsy + l*info->csy;
      list[i].tau_z = h*info->asz + k*info->bsz + l*info->csz;
      list[i].tau = sqrt(list[i].tau_x*list[i].tau_x +
                         list[i].tau_y*list[i].tau_y +
                         list[i].tau_z*list[i].tau_z);
      list[i].u1x = list[i].tau_x/list[i].tau;
      list[i].u1y = list[i].tau_y/list[i].tau;
      list[i].u1z = list[i].tau_z/list[i].tau;
      sig1 = FWHM2RMS*info->m_delta_d_d*list[i].tau;

      /* Find two arbitrary axes perpendicular to tau and each other. */
      normal_vec(&b1[0], &b1[1], &b1[2],
                 list[i].u1x, list[i].u1y, list[i].u1z);
      vec_prod(b2[0], b2[1], b2[2],
               list[i].u1x, list[i].u1y, list[i].u1z,
               b1[0], b1[1], b1[2]);

      /* Find the two mosaic axes perpendicular to tau. */
      if(SC_mosaic > 0) {
        /* Use isotropic mosaic. */
        list[i].u2x = b1[0];
        list[i].u2y = b1[1];
        list[i].u2z = b1[2];
        sig2 = FWHM2RMS*list[i].tau*MIN2RAD*SC_mosaic;
        list[i].u3x = b2[0];
        list[i].u3y = b2[1];
        list[i].u3z = b2[2];
        sig3 = FWHM2RMS*list[i].tau*MIN2RAD*SC_mosaic;
      } else if(SC_mosaic_a > 0 && SC_mosaic_b > 0 && SC_mosaic_c > 0) {
        /* Use anisotropic mosaic. */
        fprintf(stderr,"Single_crystal: Warning: you are using an experimental feature:\n"
          "  anistropic mosaicity. Please examine your data carefully.\n");
        /* compute the jacobian of (tau_v,tau_n) from rotations around the unit cell vectors. */
        struct hkl_data *l =&(list[i]);
        double xia_x,xia_y,xia_z,xib_x,xib_y,xib_z,xic_x,xic_y,xic_z;
        /*input parameters are in arc minutes*/
        double sig_fi_a=SC_mosaic_a*MIN2RAD;
        double sig_fi_b=SC_mosaic_b*MIN2RAD;
        double sig_fi_c=SC_mosaic_c*MIN2RAD;
        if(info->m_a==0) info->m_a=sqrt(scalar_prod( info->m_ax,info->m_ay,info->m_az,info->m_ax,info->m_ay,info->m_az));
        if(info->m_b==0) info->m_b=sqrt(scalar_prod( info->m_bx,info->m_by,info->m_bz,info->m_bx,info->m_by,info->m_bz));
        if(info->m_c==0) info->m_c=sqrt(scalar_prod( info->m_cx,info->m_cy,info->m_cz,info->m_cx,info->m_cy,info->m_cz));

        l->u2x = b1[0];
        l->u2y = b1[1];
        l->u2z = b1[2];
        l->u3x = b2[0];
        l->u3y = b2[1];
        l->u3z = b2[2];

        xia_x=l->tau_x-(M_2_PI*h/info->m_a)*info->asx;
        xia_y=l->tau_y-(M_2_PI*h/info->m_a)*info->asy;
        xia_z=l->tau_z-(M_2_PI*h/info->m_a)*info->asz;
        xib_x=l->tau_x-(M_2_PI*h/info->m_b)*info->bsx;
        xib_y=l->tau_y-(M_2_PI*h/info->m_b)*info->bsy;
        xib_z=l->tau_z-(M_2_PI*h/info->m_b)*info->bsz;
        xic_x=l->tau_x-(M_2_PI*h/info->m_c)*info->csx;
        xic_y=l->tau_y-(M_2_PI*h/info->m_c)*info->csy;
        xic_z=l->tau_z-(M_2_PI*h/info->m_c)*info->csz;

        double xia=sqrt(xia_x*xia_x + xia_y*xia_y + xia_z*xia_z);
        double xib=sqrt(xib_x*xib_x + xib_y*xib_y + xib_z*xib_z);
        double xic=sqrt(xic_x*xic_x + xic_y*xic_y + xic_z*xic_z);

        vec_prod(tmp_x,tmp_y,tmp_z,l->tau_x,l->tau_y,l->tau_z, l->u2x,l->u2y,l->u2z);
        double J_n_fia= xia/info->m_a/l->tau*scalar_prod(info->asx,info->asy,info->asz,tmp_x,tmp_y,tmp_z);
        vec_prod(tmp_x,tmp_y,tmp_z,l->tau_x,l->tau_y,l->tau_z, l->u2x,l->u2y,l->u2z);
        double J_n_fib= xib/info->m_b/l->tau*scalar_prod(info->bsx,info->bsy,info->bsz,tmp_x,tmp_y,tmp_z);
        vec_prod(tmp_x,tmp_y,tmp_z,l->tau_x,l->tau_y,l->tau_z, l->u2x,l->u2y,l->u2z);
        double J_n_fic= xic/info->m_c/l->tau*scalar_prod(info->csx,info->csy,info->csz,tmp_x,tmp_y,tmp_z);

        vec_prod(tmp_x,tmp_y,tmp_z,l->tau_x,l->tau_y,l->tau_z, l->u3x,l->u3y,l->u3z);
        double J_v_fia= xia/info->m_a/l->tau*scalar_prod(info->asx,info->asy,info->asz,tmp_x,tmp_y,tmp_z);
        vec_prod(tmp_x,tmp_y,tmp_z,l->tau_x,l->tau_y,l->tau_z, l->u3x,l->u3y,l->u3z);
        double J_v_fib= xib/info->m_b/l->tau*scalar_prod(info->bsx,info->bsy,info->bsz,tmp_x,tmp_y,tmp_z);
        vec_prod(tmp_x,tmp_y,tmp_z,l->tau_x,l->tau_y,l->tau_z, l->u3x,l->u3y,l->u3z);
        double J_v_fic= xic/info->m_c/l->tau*scalar_prod(info->csx,info->csy,info->csz,tmp_x,tmp_y,tmp_z);

        /*with the jacobian we can compute the sigmas in terms of the orthogonal vectors u2 and u3*/
        sig2=sig_fi_a*fabs(J_v_fia) + sig_fi_b*fabs(J_v_fib) + sig_fi_c*fabs(J_v_fic);
        sig3=sig_fi_a*fabs(J_n_fia) + sig_fi_b*fabs(J_n_fib) + sig_fi_c*fabs(J_n_fic);
      } else if (SC_mosaic_AB[0]!=0 && SC_mosaic_AB[1]!=0){
        if ( (SC_mosaic_AB[2]==0 && SC_mosaic_AB[3]==0 && SC_mosaic_AB[4]==0) || (SC_mosaic_AB[5]==0 && SC_mosaic_AB[6]==0 && SC_mosaic_AB[7]==0) ){
          fprintf(stderr,"Single_crystal: Error: in-plane mosaics are specified but one (or both)\n"
              "  in-plane reciprocal vector is the zero vector\n");
          return(0);
        }
        fprintf(stderr,"Single_crystal: Warning: you are using an experimental feature: \n"
              "  \"in-plane\" anistropic mosaicity. Please examine your data carefully.\n");

        /*for given reflection in list - compute linear comb of tau_a and tau_b*/
        /*check for not in plane - f.i. check if (tau_a X tau_b).tau_i)==0*/
        struct hkl_data *l =&(list[i]);
        double det,c1,c2,sig_tau_c;
        double em_x,em_y,em_z, tmp_x,tmp_y,tmp_z;
        double tau_a[3],tau_b[3];
        /*convert Miller indices to taus*/
        if(info->m_a==0) info->m_a=sqrt(scalar_prod( info->m_ax,info->m_ay,info->m_az,info->m_ax,info->m_ay,info->m_az));
        if(info->m_b==0) info->m_b=sqrt(scalar_prod( info->m_bx,info->m_by,info->m_bz,info->m_bx,info->m_by,info->m_bz));
        if(info->m_c==0) info->m_c=sqrt(scalar_prod( info->m_cx,info->m_cy,info->m_cz,info->m_cx,info->m_cy,info->m_cz));
        tau_a[0]=M_2_PI*( (SC_mosaic_AB[2]/info->m_a)*info->asx + (SC_mosaic_AB[3]/info->m_b)*info->bsx + (SC_mosaic_AB[4]/info->m_c)*info->csx );
        tau_a[1]=M_2_PI*( (SC_mosaic_AB[2]/info->m_a)*info->asy + (SC_mosaic_AB[3]/info->m_b)*info->bsy + (SC_mosaic_AB[4]/info->m_c)*info->csy );
        tau_a[2]=M_2_PI*( (SC_mosaic_AB[2]/info->m_a)*info->asz + (SC_mosaic_AB[3]/info->m_b)*info->bsz + (SC_mosaic_AB[4]/info->m_c)*info->csz );
        tau_b[0]=M_2_PI*( (SC_mosaic_AB[5]/info->m_a)*info->asx + (SC_mosaic_AB[6]/info->m_b)*info->bsx + (SC_mosaic_AB[7]/info->m_c)*info->csx );
        tau_b[1]=M_2_PI*( (SC_mosaic_AB[5]/info->m_a)*info->asy + (SC_mosaic_AB[6]/info->m_b)*info->bsy + (SC_mosaic_AB[7]/info->m_c)*info->csy );
        tau_b[2]=M_2_PI*( (SC_mosaic_AB[5]/info->m_a)*info->asz + (SC_mosaic_AB[6]/info->m_b)*info->bsz + (SC_mosaic_AB[7]/info->m_c)*info->csz );

        /*check determinants to see how we should compute the linear combination of a and b (to match c)*/
        if ((det=tau_a[0]*tau_b[1]-tau_a[1]*tau_b[0])!=0){
          c1= (l->tau_x*tau_b[1] - l->tau_y*tau_b[0])/det;
          c2= (tau_a[0]*l->tau_y - tau_a[1]*l->tau_x)/det;
        }else if ((det=tau_a[1]*tau_b[2]-tau_a[2]*tau_b[1])!=0){
          c1= (l->tau_y*tau_b[2] - l->tau_z*tau_b[1])/det;
          c2= (tau_a[1]*l->tau_z - tau_a[2]*l->tau_y)/det;
        }else if ((det=tau_a[0]*tau_b[2]-tau_a[2]*tau_b[0])!=0){
          c1= (l->tau_x*tau_b[2] - l->tau_z*tau_b[0])/det;
          c2= (tau_a[0]*l->tau_z - tau_a[2]*l->tau_x)/det;
        }
        if ((c1==0) && (c2==0)){
          fprintf(stderr,"Single_crystal: Warning: reflection tau[%i]=(%g %g %g) "
          "has no component in defined mosaic plane\n",
          i, l->tau_x,l->tau_y,l->tau_z);
        }
        /*compute linear combination => sig_tau_i = | c1*sig_tau_a + c2*sig_tau_b |  - also add in the minute to radian scaling factor*/;
        sig_tau_c = MIN2RAD*sqrt(c1*SC_mosaic_AB[0]*c1*SC_mosaic_AB[0] + c2*SC_mosaic_AB[1]*c2*SC_mosaic_AB[1]);
        l->u2x = b1[0]; l->u2y = b1[1]; l->u2z = b1[2];
        l->u3x = b2[0]; l->u3y = b2[1]; l->u3z = b2[2];

        /*so now let's compute the rotation around planenormal tau_a X tau_b*/
        /*g_bar (unit normal of rotation plane) = tau_a X tau_b / norm(tau_a X tau_b)*/
        vec_prod(tmp_x,tmp_y,tmp_z, tau_a[0],tau_a[1],tau_a[2],tau_b[0],tau_b[1],tau_b[2]);
        vec_prod(em_x,em_y,em_z, l->tau_x, l->tau_y, l->tau_z, tmp_x,tmp_y,tmp_z);
        NORM(em_x,em_y,em_z);
        sig2 = l->tau*sig_tau_c*fabs(scalar_prod(em_x,em_y,em_z, l->u2x,l->u2y,l->u2z));
        sig3 = l->tau*sig_tau_c*fabs(scalar_prod(em_x,em_y,em_z, l->u3x,l->u3y,l->u3z));
        /*protect against collapsing gaussians. These seem to be sensible values.*/
        if (sig2<1e-5) sig2=1e-5;
        if (sig3<1e-5) sig3=1e-5;
      }
      else {
        fprintf(stderr,
                "Single_crystal: Error: EITHER mosaic OR (mosaic_a, mosaic_b, mosaic_c)\n"
                "  must be given and be >0.\n");
        return(0);
      }
      list[i].sig123 = sig1*sig2*sig3;
      list[i].m1 = 1/(2*sig1*sig1);
      list[i].m2 = 1/(2*sig2*sig2);
      list[i].m3 = 1/(2*sig3*sig3);
      /* Set Gauss cutoff to 5 times the maximal sigma. */
      if(sig1 > sig2)
        if(sig1 > sig3)
          list[i].cutoff = 5*sig1;
        else
          list[i].cutoff = 5*sig3;
      else
        if(sig2 > sig3)
          list[i].cutoff = 5*sig2;
        else
          list[i].cutoff = 5*sig3;
    }
    Table_Free(&sTable);

    /* sort the list with increasing tau */
    qsort(list, i, sizeof(struct hkl_data),  SX_list_compare);

    info->list = list;
    info->count = i;
    info->tau_list = malloc(i*sizeof(*info->tau_list));
    if(!info->tau_list)
    {
      fprintf(stderr, "Single_crystal: Error: Out of memory!\n");
      return(0);
    }
    return(info->count);
  } /* read_hkl_data */

  /* ------------------------------------------------------------------------ */
  /* hkl_search
    search the HKL reflections which are on the Ewald sphere
    input:
      L,T,count,V0: constants for all calls
      kix,kiy,kiz: may be different for each call
    this function returns:
      tau_count (return), coh_refl, coh_xsect, T (updated elements in the array up to [j])
   */

#ifdef USE_OPENCL // when using OpenCL
int hkl_search_opencl(struct hkl_data *h_L, struct tau_data* h_T, int count, double V0, 
  double kix, double kiy, double kiz, double tau_max,
  double* h_coh_refl, double* h_coh_xsect, struct opencl_context oclContext,
  cl_mem d_L, cl_mem d_T, cl_mem d_tau_count, cl_mem d_coh_refl, cl_mem d_coh_xsect)
{ 
    cl_int             ciErr1;                              // Error code var
    size_t             globalWorkSize[1] = {4096};          // 1D var for Total # of work items
    size_t             localWorkSize[1] = {64};            // 1D var for # of work items in the work group	
    
    int h_tau_count=0;

    //printf("In GPU-based hkl_search");

    //  Using GPU(s)...
    if (oclContext.Kernel == NULL || oclContext.GPUContext == NULL) {
      printf("OpenCL Kernel could not be initialized\n");
      return(-10); // OpenCL Kernel could not be initialized
    }
    
    // Write data 
    int iDevice = 0;
    ciErr1 = clEnqueueWriteBuffer(oclContext.CommandQueue[iDevice], d_L, CL_TRUE, 0, sizeof(struct hkl_data)*count, h_L, 0, NULL, NULL);
    if (ciErr1 !=CL_SUCCESS) {
      return(-1);
    }
    ciErr1 = clEnqueueWriteBuffer(oclContext.CommandQueue[iDevice], d_T, CL_TRUE, 0, sizeof(struct tau_data)*count, h_T, 0, NULL, NULL);
    if (ciErr1 !=CL_SUCCESS) {
      return(-2);
    }
    ciErr1 = clEnqueueWriteBuffer(oclContext.CommandQueue[iDevice], d_tau_count, CL_TRUE, 0, sizeof(int), &h_tau_count, 0, NULL, NULL);
    if (ciErr1 !=CL_SUCCESS){
      return(-3);
    }
    ciErr1 = clEnqueueWriteBuffer(oclContext.CommandQueue[iDevice], d_coh_refl, CL_TRUE, 0, sizeof(double), h_coh_refl, 0, NULL, NULL); 
    if (ciErr1 !=CL_SUCCESS){
      return(-4);
    }
    ciErr1 = clEnqueueWriteBuffer(oclContext.CommandQueue[iDevice], d_coh_xsect, CL_TRUE, 0, sizeof(double), h_coh_xsect , 0, NULL, NULL); 
    if (ciErr1 !=CL_SUCCESS) {
      return(-5);
    }

    clFinish(oclContext.CommandQueue[iDevice]);    
	
    ciErr1  = clSetKernelArg(oclContext.Kernel, 0, sizeof(cl_mem), (void*)&d_L);
    ciErr1 |= clSetKernelArg(oclContext.Kernel, 1, sizeof(cl_mem), (void*)&d_T);
    ciErr1 |= clSetKernelArg(oclContext.Kernel, 2, sizeof(int),    (void*)&count);
    ciErr1 |= clSetKernelArg(oclContext.Kernel, 3, sizeof(double), (void*)&V0);
    ciErr1 |= clSetKernelArg(oclContext.Kernel, 4, sizeof(double), (void*)&kix);
    ciErr1 |= clSetKernelArg(oclContext.Kernel, 5, sizeof(double), (void*)&kiy);
    ciErr1 |= clSetKernelArg(oclContext.Kernel, 6, sizeof(double), (void*)&kiz);
    ciErr1 |= clSetKernelArg(oclContext.Kernel, 7, sizeof(double), (void*)&tau_max);
    ciErr1 |= clSetKernelArg(oclContext.Kernel, 8, sizeof(cl_mem), (void*)&d_tau_count);
    ciErr1 |= clSetKernelArg(oclContext.Kernel, 9, sizeof(cl_mem), (void*)&d_coh_refl);
    ciErr1 |= clSetKernelArg(oclContext.Kernel,10, sizeof(cl_mem), (void*)&d_coh_xsect);
    if (ciErr1 !=CL_SUCCESS) return(0);

    ciErr1 = clEnqueueNDRangeKernel(oclContext.CommandQueue[iDevice], oclContext.Kernel, 1, NULL, 
        globalWorkSize, localWorkSize, 0, NULL, NULL); 
    if (ciErr1 !=CL_SUCCESS) return (0);
                     	
    clFinish(oclContext.CommandQueue[iDevice]);

    //Read back results...
    ciErr1 = clEnqueueReadBuffer(oclContext.CommandQueue[iDevice], d_T, CL_TRUE, 0, sizeof(struct tau_data)*count, h_T, 0, NULL, NULL);
    if (ciErr1 !=CL_SUCCESS) {
      return(-6);
    }
    
    ciErr1 = clEnqueueReadBuffer(oclContext.CommandQueue[iDevice], d_tau_count, CL_TRUE, 0, sizeof(int), &h_tau_count, 0, NULL, NULL);
    if (ciErr1 !=CL_SUCCESS) {
      return (-7); 
    }
    
    ciErr1 = clEnqueueReadBuffer(oclContext.CommandQueue[iDevice], d_coh_refl, CL_TRUE, 0, sizeof(double), h_coh_refl, 0, NULL, NULL);
    if (ciErr1 !=CL_SUCCESS) {
      return (-8);      
    }
    
    ciErr1 = clEnqueueReadBuffer(oclContext.CommandQueue[iDevice], d_coh_xsect, CL_TRUE, 0, sizeof(double), h_coh_xsect, 0, NULL, NULL);
    if (ciErr1 !=CL_SUCCESS) {
      return (-9);
    }
    return h_tau_count;
} 
#endif 
  // non OpenCL case
  int hkl_search(struct hkl_data *L, struct tau_data *T, int count, double V0, 
    double kix, double kiy, double kiz, double tau_max,
    double *coh_refl, double *coh_xsect)
  {
    double rho, rho_x, rho_y, rho_z;
    double diff;
    int    i,j;
    double ox,oy,oz;
    double b1x,b1y,b1z, b2x,b2y,b2z, kx, ky, kz, nx, ny, nz;
    double n11, n22, n12, det_N, inv_n11, inv_n22, inv_n12, l11, l22, l12,  det_L;
    double Bt_D_O_x, Bt_D_O_y, y0x, y0y, alpha;

    double ki = sqrt(kix*kix+kiy*kiy+kiz*kiz);

    /* Common factor in coherent cross-section */
    double xsect_factor = pow(2*PI, 5.0/2.0)/(V0*ki*ki);

    for(i = j = 0; i < count; i++)
      {
    /* Assuming reflections are sorted, stop search when max tau exceeded. */
        if(L[i].tau > tau_max)
          break;
        /* Check if this reciprocal lattice point is close enough to the
           Ewald sphere to make scattering possible. */
        rho_x = kix - L[i].tau_x;
        rho_y = kiy - L[i].tau_y;
        rho_z = kiz - L[i].tau_z;
        rho = sqrt(rho_x*rho_x + rho_y*rho_y + rho_z*rho_z);
        diff = fabs(rho - ki);

        /* Check if scattering is possible (cutoff of Gaussian tails). */
        if(diff <= L[i].cutoff)
        {
          /* Store reflection. */
          T[j].index = i;
          /* Get ki vector in local coordinates. */
          kx = kix*L[i].u1x + kiy*L[i].u1y + kiz*L[i].u1z;
          ky = kix*L[i].u2x + kiy*L[i].u2y + kiz*L[i].u2z;
          kz = kix*L[i].u3x + kiy*L[i].u3y + kiz*L[i].u3z;
          T[j].rho_x = kx - L[i].tau;
          T[j].rho_y = ky;
          T[j].rho_z = kz;
          T[j].rho = rho;
          /* Compute the tangent plane of the Ewald sphere. */
          nx = T[j].rho_x/T[j].rho;
          ny = T[j].rho_y/T[j].rho;
          nz = T[j].rho_z/T[j].rho;
          ox = (ki - T[j].rho)*nx;
          oy = (ki - T[j].rho)*ny;
          oz = (ki - T[j].rho)*nz;
          T[j].ox = ox;
          T[j].oy = oy;
          T[j].oz = oz;
          /* Compute unit vectors b1 and b2 that span the tangent plane. */
          normal_vec(&b1x, &b1y, &b1z, nx, ny, nz);
          vec_prod(b2x, b2y, b2z, nx, ny, nz, b1x, b1y, b1z);
          T[j].b1x = b1x;
          T[j].b1y = b1y;
          T[j].b1z = b1z;
          T[j].b2x = b2x;
          T[j].b2y = b2y;
          T[j].b2z = b2z;
          /* Compute the 2D projection of the 3D Gauss of the reflection. */
          /* The symmetric 2x2 matrix N describing the 2D gauss. */
          n11 = L[i].m1*b1x*b1x + L[i].m2*b1y*b1y + L[i].m3*b1z*b1z;
          n12 = L[i].m1*b1x*b2x + L[i].m2*b1y*b2y + L[i].m3*b1z*b2z;
          n22 = L[i].m1*b2x*b2x + L[i].m2*b2y*b2y + L[i].m3*b2z*b2z;
          /* The (symmetric) inverse matrix of N. */
          det_N = n11*n22 - n12*n12;
          inv_n11 = n22/det_N;
          inv_n12 = -n12/det_N;
          inv_n22 = n11/det_N;
          /* The Cholesky decomposition of 1/2*inv_n (lower triangular L). */
          l11 = sqrt(inv_n11/2);
          l12 = inv_n12/(2*l11);
          l22 = sqrt(inv_n22/2 - l12*l12);
          T[j].l11 = l11;
          T[j].l12 = l12;
          T[j].l22 = l22;
          det_L = l11*l22;
          /* The product B^T D o. */
          Bt_D_O_x = b1x*L[i].m1*ox + b1y*L[i].m2*oy + b1z*L[i].m3*oz;
          Bt_D_O_y = b2x*L[i].m1*ox + b2y*L[i].m2*oy + b2z*L[i].m3*oz;
          /* Center of 2D Gauss in plane coordinates. */
          y0x = -(Bt_D_O_x*inv_n11 + Bt_D_O_y*inv_n12);
          y0y = -(Bt_D_O_x*inv_n12 + Bt_D_O_y*inv_n22);
          T[j].y0x = y0x;
          T[j].y0y = y0y;
          /* Factor alpha for the distance of the 2D Gauss from the origin. */
          alpha = L[i].m1*ox*ox + L[i].m2*oy*oy + L[i].m3*oz*oz -
                       (y0x*y0x*n11 + y0y*y0y*n22 + 2*y0x*y0y*n12);
          T[j].refl = xsect_factor*det_L*exp(-alpha)/L[i].sig123; /* intensity of that Bragg */
          *coh_refl += T[j].refl;                                  /* total scatterable intensity */
          T[j].xsect = T[j].refl*L[i].F2;
          *coh_xsect += T[j].xsect;
          j++;
        }

      } /* end for */
        return (j); // this is 'tau_count', i.e. number of reachable reflections
    } /* end hkl_search */

    int hkl_select(struct tau_data *T, int tau_count, double coh_refl, double *sum) {
      int j;
      double r = rand0max(coh_refl);
      *sum = 0;
      for(j = 0; j < tau_count; j++)
      {
        *sum += T[j].refl;
        if(*sum > r) break;
      }
      return j;
    }

    /* Functions for "reorientation", powder and PG modes */
    /* Powder, forward */
    void randrotate(double *nx, double *ny, double *nz, double a, double b) {
      double nvx, nvy, nvz;
      rotate(nvx,nvy,nvz, *nx, *ny, *nz, a, 1, 0, 0);
      rotate(*nx,*ny,*nz, nvx,nvy,nvz, b, 0, 1, 0);
    }
    /* Powder, back */
    void randderotate(double *nx, double *ny, double *nz, double a, double b) {
      double nvx, nvy, nvz;
      rotate(nvx,nvy,nvz,*nx,*ny,*nz, -b, 0, 1, 0);
      rotate(*nx, *ny, *nz, nvx,nvy,nvz, -a, 1, 0, 0);
    }
    /* PG, forward */
    void PGrotate(double *nx, double *ny, double *nz, double a, double csx, double csy, double csz) {
      /* Currently assumes c-axis along 'x', ought to be generalized... */
      double nvx, nvy, nvz;
      rotate(nvx,nvy,nvz, *nx, *ny, *nz, a, csx, csy, csz);
      *nx = nvx; *ny = nvy; *nz = nvz;
    }
    /* PG, back */
    void PGderotate(double *nx, double *ny, double *nz, double a, double csx, double csy, double csz) {
      /* Currently assumes c-axis along 'x', ought to be generalized... */
      double nvx, nvy, nvz;
      rotate(nvx,nvy,nvz, *nx, *ny, *nz, -a, csx, csy, csz);
      *nx = nvx; *ny = nvy; *nz = nvz;
    }

  #ifdef USE_OPENCL
  struct opencl_context oclContext_SX;
  cl_mem             d_T=NULL, d_tau_count=NULL, d_coh_refl=NULL, d_coh_xsect=NULL;    //OpenCL device buffers
  cl_mem             d_L=NULL;
  #endif



#endif /* !SINGLE_CRYSTAL_DECL */


/* ************************************************************************** */
/*             End of SHARE user declarations for all components              */
/* ************************************************************************** */


/* ********************** component definition declarations. **************** */

/* component Origin=Progress_bar() [1] DECLARE */
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
  char     _name[256]; /* e.g. Origin */
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
_class_Progress_bar _Origin_var;
_class_Progress_bar *_Origin = &_Origin_var;
#pragma acc declare create ( _Origin_var )
#pragma acc declare create ( _Origin )

/* component Source=Source_simple() [2] DECLARE */
/* Parameter definition for component type 'Source_simple' */
struct _struct_Source_simple_parameters {
  /* Component type 'Source_simple' setting parameters */
  MCNUM radius;
  MCNUM yheight;
  MCNUM xwidth;
  MCNUM dist;
  MCNUM focus_xw;
  MCNUM focus_yh;
  MCNUM E0;
  MCNUM dE;
  MCNUM lambda0;
  MCNUM dlambda;
  MCNUM flux;
  MCNUM gauss;
  long target_index;
  /* Component type 'Source_simple' private parameters */
  /* Component type 'Source_simple' DECLARE code stored as structure members */
double pmul, srcArea;
int square;
}; /* _struct_Source_simple_parameters */
typedef struct _struct_Source_simple_parameters _class_Source_simple_parameters;

/* Parameters for component type 'Source_simple' */
struct _struct_Source_simple {
  char     _name[256]; /* e.g. Source */
  char     _type[256]; /* Source_simple */
  long     _index; /* e.g. 2 index in TRACE list */
  Coords   _position_absolute;
  Coords   _position_relative; /* wrt PREVIOUS */
  Rotation _rotation_absolute;
  Rotation _rotation_relative; /* wrt PREVIOUS */
  int      _rotation_is_identity;
  _class_Source_simple_parameters _parameters;
};
typedef struct _struct_Source_simple _class_Source_simple;
_class_Source_simple _Source_var;
_class_Source_simple *_Source = &_Source_var;
#pragma acc declare create ( _Source_var )
#pragma acc declare create ( _Source )

/* component lamStart=L_monitor() [3] DECLARE */
/* Parameter definition for component type 'L_monitor' */
struct _struct_L_monitor_parameters {
  /* Component type 'L_monitor' setting parameters */
  MCNUM nL;
  char filename[16384];
  MCNUM xmin;
  MCNUM xmax;
  MCNUM ymin;
  MCNUM ymax;
  MCNUM xwidth;
  MCNUM yheight;
  MCNUM Lmin;
  MCNUM Lmax;
  MCNUM restore_neutron;
  /* Component type 'L_monitor' private parameters */
  /* Component type 'L_monitor' DECLARE code stored as structure members */
  DArray1d L_N;
  DArray1d L_p;
  DArray1d L_p2;
}; /* _struct_L_monitor_parameters */
typedef struct _struct_L_monitor_parameters _class_L_monitor_parameters;

/* Parameters for component type 'L_monitor' */
struct _struct_L_monitor {
  char     _name[256]; /* e.g. lamStart */
  char     _type[256]; /* L_monitor */
  long     _index; /* e.g. 2 index in TRACE list */
  Coords   _position_absolute;
  Coords   _position_relative; /* wrt PREVIOUS */
  Rotation _rotation_absolute;
  Rotation _rotation_relative; /* wrt PREVIOUS */
  int      _rotation_is_identity;
  _class_L_monitor_parameters _parameters;
};
typedef struct _struct_L_monitor _class_L_monitor;
_class_L_monitor _lamStart_var;
_class_L_monitor *_lamStart = &_lamStart_var;
#pragma acc declare create ( _lamStart_var )
#pragma acc declare create ( _lamStart )

/* component Mono_Arm=Arm() [4] DECLARE */
/* Parameter definition for component type 'Arm' */
struct _struct_Arm_parameters {
  char Arm_has_no_parameters;
}; /* _struct_Arm_parameters */
typedef struct _struct_Arm_parameters _class_Arm_parameters;

/* Parameters for component type 'Arm' */
struct _struct_Arm {
  char     _name[256]; /* e.g. Mono_Arm */
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
_class_Arm _Mono_Arm_var;
_class_Arm *_Mono_Arm = &_Mono_Arm_var;
#pragma acc declare create ( _Mono_Arm_var )
#pragma acc declare create ( _Mono_Arm )

_class_Arm _Mono_Out_var;
_class_Arm *_Mono_Out = &_Mono_Out_var;
#pragma acc declare create ( _Mono_Out_var )
#pragma acc declare create ( _Mono_Out )

/* component Mono1=Monochromator_flat() [6] DECLARE */
/* Parameter definition for component type 'Monochromator_flat' */
struct _struct_Monochromator_flat_parameters {
  /* Component type 'Monochromator_flat' setting parameters */
  MCNUM zmin;
  MCNUM zmax;
  MCNUM ymin;
  MCNUM ymax;
  MCNUM zwidth;
  MCNUM yheight;
  MCNUM mosaich;
  MCNUM mosaicv;
  MCNUM r0;
  MCNUM Q;
  MCNUM DM;
  /* Component type 'Monochromator_flat' private parameters */
  /* Component type 'Monochromator_flat' DECLARE code stored as structure members */
  double mos_rms_y;                                             
  double mos_rms_z;
  double mos_rms_max;
  double mono_Q;
}; /* _struct_Monochromator_flat_parameters */
typedef struct _struct_Monochromator_flat_parameters _class_Monochromator_flat_parameters;

/* Parameters for component type 'Monochromator_flat' */
struct _struct_Monochromator_flat {
  char     _name[256]; /* e.g. Mono1 */
  char     _type[256]; /* Monochromator_flat */
  long     _index; /* e.g. 2 index in TRACE list */
  Coords   _position_absolute;
  Coords   _position_relative; /* wrt PREVIOUS */
  Rotation _rotation_absolute;
  Rotation _rotation_relative; /* wrt PREVIOUS */
  int      _rotation_is_identity;
  _class_Monochromator_flat_parameters _parameters;
};
typedef struct _struct_Monochromator_flat _class_Monochromator_flat;
_class_Monochromator_flat _Mono1_var;
_class_Monochromator_flat *_Mono1 = &_Mono1_var;
#pragma acc declare create ( _Mono1_var )
#pragma acc declare create ( _Mono1 )

/* component Mono2=Monochromator_pol() [7] DECLARE */
/* Parameter definition for component type 'Monochromator_pol' */
struct _struct_Monochromator_pol_parameters {
  /* Component type 'Monochromator_pol' setting parameters */
  MCNUM zwidth;
  MCNUM yheight;
  MCNUM mosaic;
  MCNUM dspread;
  MCNUM Q;
  MCNUM DM;
  MCNUM pThreshold;
  MCNUM Rup;
  MCNUM Rdown;
  long debug;
  /* Component type 'Monochromator_pol' private parameters */
  /* Component type 'Monochromator_pol' DECLARE code stored as structure members */
double mos_rms;                                             
double d_rms;                                            
double mono_Q;

double FN;                                         
double FM;                                          
}; /* _struct_Monochromator_pol_parameters */
typedef struct _struct_Monochromator_pol_parameters _class_Monochromator_pol_parameters;

/* Parameters for component type 'Monochromator_pol' */
struct _struct_Monochromator_pol {
  char     _name[256]; /* e.g. Mono2 */
  char     _type[256]; /* Monochromator_pol */
  long     _index; /* e.g. 2 index in TRACE list */
  Coords   _position_absolute;
  Coords   _position_relative; /* wrt PREVIOUS */
  Rotation _rotation_absolute;
  Rotation _rotation_relative; /* wrt PREVIOUS */
  int      _rotation_is_identity;
  _class_Monochromator_pol_parameters _parameters;
};
typedef struct _struct_Monochromator_pol _class_Monochromator_pol;
_class_Monochromator_pol _Mono2_var;
_class_Monochromator_pol *_Mono2 = &_Mono2_var;
#pragma acc declare create ( _Mono2_var )
#pragma acc declare create ( _Mono2 )

_class_Monochromator_pol _Mono3_var;
_class_Monochromator_pol *_Mono3 = &_Mono3_var;
#pragma acc declare create ( _Mono3_var )
#pragma acc declare create ( _Mono3 )

/* component Mono4=Monochromator_curved() [9] DECLARE */
/* Parameter definition for component type 'Monochromator_curved' */
struct _struct_Monochromator_curved_parameters {
  /* Component type 'Monochromator_curved' setting parameters */
  char reflect[16384];
  char transmit[16384];
  MCNUM zwidth;
  MCNUM yheight;
  MCNUM gap;
  MCNUM NH;
  MCNUM NV;
  MCNUM mosaich;
  MCNUM mosaicv;
  MCNUM r0;
  MCNUM t0;
  MCNUM Q;
  MCNUM RV;
  MCNUM RH;
  MCNUM DM;
  MCNUM mosaic;
  MCNUM width;
  MCNUM height;
  MCNUM verbose;
  MCNUM order;
  /* Component type 'Monochromator_curved' private parameters */
  /* Component type 'Monochromator_curved' DECLARE code stored as structure members */
  double mos_rms_y;                                             
  double mos_rms_z;
  double mos_rms_max;
  double mono_Q;
  double SlabWidth, SlabHeight;
  t_Table rTable, tTable;
  double row,col;
  double* tiltH;
  double* tiltV;
}; /* _struct_Monochromator_curved_parameters */
typedef struct _struct_Monochromator_curved_parameters _class_Monochromator_curved_parameters;

/* Parameters for component type 'Monochromator_curved' */
struct _struct_Monochromator_curved {
  char     _name[256]; /* e.g. Mono4 */
  char     _type[256]; /* Monochromator_curved */
  long     _index; /* e.g. 2 index in TRACE list */
  Coords   _position_absolute;
  Coords   _position_relative; /* wrt PREVIOUS */
  Rotation _rotation_absolute;
  Rotation _rotation_relative; /* wrt PREVIOUS */
  int      _rotation_is_identity;
  _class_Monochromator_curved_parameters _parameters;
};
typedef struct _struct_Monochromator_curved _class_Monochromator_curved;
_class_Monochromator_curved _Mono4_var;
_class_Monochromator_curved *_Mono4 = &_Mono4_var;
#pragma acc declare create ( _Mono4_var )
#pragma acc declare create ( _Mono4 )

/* component Mono5=Monochromator_2foc() [10] DECLARE */
/* Parameter definition for component type 'Monochromator_2foc' */
struct _struct_Monochromator_2foc_parameters {
  /* Component type 'Monochromator_2foc' setting parameters */
  char reflect[16384];
  MCNUM zwidth;
  MCNUM yheight;
  MCNUM gap;
  MCNUM NH;
  MCNUM NV;
  MCNUM mosaich;
  MCNUM mosaicv;
  MCNUM r0;
  MCNUM Q;
  MCNUM RV;
  MCNUM RH;
  MCNUM DM;
  MCNUM mosaic;
  MCNUM width;
  MCNUM height;
  MCNUM verbose;
  /* Component type 'Monochromator_2foc' private parameters */
  /* Component type 'Monochromator_2foc' DECLARE code stored as structure members */
double mos_y;                         
double mos_z;
double mono_Q;
double SlabWidth;
double SlabHeight;
t_Table rTable;
}; /* _struct_Monochromator_2foc_parameters */
typedef struct _struct_Monochromator_2foc_parameters _class_Monochromator_2foc_parameters;

/* Parameters for component type 'Monochromator_2foc' */
struct _struct_Monochromator_2foc {
  char     _name[256]; /* e.g. Mono5 */
  char     _type[256]; /* Monochromator_2foc */
  long     _index; /* e.g. 2 index in TRACE list */
  Coords   _position_absolute;
  Coords   _position_relative; /* wrt PREVIOUS */
  Rotation _rotation_absolute;
  Rotation _rotation_relative; /* wrt PREVIOUS */
  int      _rotation_is_identity;
  _class_Monochromator_2foc_parameters _parameters;
};
typedef struct _struct_Monochromator_2foc _class_Monochromator_2foc;
_class_Monochromator_2foc _Mono5_var;
_class_Monochromator_2foc *_Mono5 = &_Mono5_var;
#pragma acc declare create ( _Mono5_var )
#pragma acc declare create ( _Mono5 )

/* component Mono6=Single_crystal() [11] DECLARE */
/* Parameter definition for component type 'Single_crystal' */
struct _struct_Single_crystal_parameters {
  /* Component type 'Single_crystal' setting parameters */
  char reflections[16384];
  char geometry[16384];
  MCNUM* mosaic_AB;
  MCNUM xwidth;
  MCNUM yheight;
  MCNUM zdepth;
  MCNUM radius;
  MCNUM delta_d_d;
  MCNUM mosaic;
  MCNUM mosaic_a;
  MCNUM mosaic_b;
  MCNUM mosaic_c;
  MCNUM recip_cell;
  MCNUM barns;
  MCNUM ax;
  MCNUM ay;
  MCNUM az;
  MCNUM bx;
  MCNUM by;
  MCNUM bz;
  MCNUM cx;
  MCNUM cy;
  MCNUM cz;
  MCNUM p_transmit;
  MCNUM sigma_abs;
  MCNUM sigma_inc;
  MCNUM aa;
  MCNUM bb;
  MCNUM cc;
  MCNUM order;
  MCNUM RX;
  MCNUM RY;
  MCNUM RZ;
  MCNUM powder;
  MCNUM PG;
  /* Component type 'Single_crystal' private parameters */
  /* Component type 'Single_crystal' DECLARE code stored as structure members */
  struct hkl_info_struct hkl_info;
  off_struct             offdata;
}; /* _struct_Single_crystal_parameters */
typedef struct _struct_Single_crystal_parameters _class_Single_crystal_parameters;

/* Parameters for component type 'Single_crystal' */
struct _struct_Single_crystal {
  char     _name[256]; /* e.g. Mono6 */
  char     _type[256]; /* Single_crystal */
  long     _index; /* e.g. 2 index in TRACE list */
  Coords   _position_absolute;
  Coords   _position_relative; /* wrt PREVIOUS */
  Rotation _rotation_absolute;
  Rotation _rotation_relative; /* wrt PREVIOUS */
  int      _rotation_is_identity;
  _class_Single_crystal_parameters _parameters;
};
typedef struct _struct_Single_crystal _class_Single_crystal;
_class_Single_crystal _Mono6_var;
_class_Single_crystal *_Mono6 = &_Mono6_var;
#pragma acc declare create ( _Mono6_var )
#pragma acc declare create ( _Mono6 )

/* component Sphere1=PSD_monitor_4PI() [12] DECLARE */
/* Parameter definition for component type 'PSD_monitor_4PI' */
struct _struct_PSD_monitor_4PI_parameters {
  /* Component type 'PSD_monitor_4PI' setting parameters */
  MCNUM nx;
  MCNUM ny;
  char filename[16384];
  MCNUM radius;
  MCNUM restore_neutron;
  /* Component type 'PSD_monitor_4PI' private parameters */
  /* Component type 'PSD_monitor_4PI' DECLARE code stored as structure members */
  DArray2d PSD_N;
  DArray2d PSD_p;
  DArray2d PSD_p2;
}; /* _struct_PSD_monitor_4PI_parameters */
typedef struct _struct_PSD_monitor_4PI_parameters _class_PSD_monitor_4PI_parameters;

/* Parameters for component type 'PSD_monitor_4PI' */
struct _struct_PSD_monitor_4PI {
  char     _name[256]; /* e.g. Sphere1 */
  char     _type[256]; /* PSD_monitor_4PI */
  long     _index; /* e.g. 2 index in TRACE list */
  Coords   _position_absolute;
  Coords   _position_relative; /* wrt PREVIOUS */
  Rotation _rotation_absolute;
  Rotation _rotation_relative; /* wrt PREVIOUS */
  int      _rotation_is_identity;
  _class_PSD_monitor_4PI_parameters _parameters;
};
typedef struct _struct_PSD_monitor_4PI _class_PSD_monitor_4PI;
_class_PSD_monitor_4PI _Sphere1_var;
_class_PSD_monitor_4PI *_Sphere1 = &_Sphere1_var;
#pragma acc declare create ( _Sphere1_var )
#pragma acc declare create ( _Sphere1 )

_class_L_monitor _lam1_var;
_class_L_monitor *_lam1 = &_lam1_var;
#pragma acc declare create ( _lam1_var )
#pragma acc declare create ( _lam1 )

/* component psd1=PSD_monitor() [14] DECLARE */
/* Parameter definition for component type 'PSD_monitor' */
struct _struct_PSD_monitor_parameters {
  /* Component type 'PSD_monitor' setting parameters */
  MCNUM nx;
  MCNUM ny;
  char filename[16384];
  MCNUM xmin;
  MCNUM xmax;
  MCNUM ymin;
  MCNUM ymax;
  MCNUM xwidth;
  MCNUM yheight;
  MCNUM restore_neutron;
  /* Component type 'PSD_monitor' private parameters */
  /* Component type 'PSD_monitor' DECLARE code stored as structure members */
  DArray2d PSD_N;
  DArray2d PSD_p;
  DArray2d PSD_p2;
}; /* _struct_PSD_monitor_parameters */
typedef struct _struct_PSD_monitor_parameters _class_PSD_monitor_parameters;

/* Parameters for component type 'PSD_monitor' */
struct _struct_PSD_monitor {
  char     _name[256]; /* e.g. psd1 */
  char     _type[256]; /* PSD_monitor */
  long     _index; /* e.g. 2 index in TRACE list */
  Coords   _position_absolute;
  Coords   _position_relative; /* wrt PREVIOUS */
  Rotation _rotation_absolute;
  Rotation _rotation_relative; /* wrt PREVIOUS */
  int      _rotation_is_identity;
  _class_PSD_monitor_parameters _parameters;
};
typedef struct _struct_PSD_monitor _class_PSD_monitor;
_class_PSD_monitor _psd1_var;
_class_PSD_monitor *_psd1 = &_psd1_var;
#pragma acc declare create ( _psd1_var )
#pragma acc declare create ( _psd1 )

int mcNUMCOMP = 14;

/* User declarations from instrument definition. Can define functions. */
  double DM  = 3.3539; /* Monochromator d-spacing in Angs */
                       /* PG002 Orders : 1st 3.355 2e 1.6775, 3e 1.1183 */

  /* to compute */
  double A1,A2;
  double mono_q;

  /* This variable helps us switch on and off the different setups*/
  double filterProb;

#undef compcurname
#undef compcurtype
#undef compcurindex
/* end of instrument 'Test_Monochromators' and components DECLARE */

/* *****************************************************************************
* instrument 'Test_Monochromators' and components INITIALISE
***************************************************************************** */

/* component Origin=Progress_bar() SETTING, POSITION/ROTATION */
int _Origin_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Origin_setpos] component Origin=Progress_bar() SETTING [/usr/share/mcstas/3.0-dev/misc/Progress_bar.comp:57]");
  stracpy(_Origin->_name, "Origin", 16384);
  stracpy(_Origin->_type, "Progress_bar", 16384);
  _Origin->_index=1;
  if("NULL" && strlen("NULL"))
    stracpy(_Origin->_parameters.profile, "NULL" ? "NULL" : "", 16384);
  else 
  _Origin->_parameters.profile[0]='\0';
  #define profile (_Origin->_parameters.profile)
  _Origin->_parameters.percent = 10;
  #define percent (_Origin->_parameters.percent)
  _Origin->_parameters.flag_save = 0;
  #define flag_save (_Origin->_parameters.flag_save)
  _Origin->_parameters.minutes = 0;
  #define minutes (_Origin->_parameters.minutes)

  #define IntermediateCnts (_Origin->_parameters.IntermediateCnts)
  #define StartTime (_Origin->_parameters.StartTime)
  #define EndTime (_Origin->_parameters.EndTime)
  #define CurrentTime (_Origin->_parameters.CurrentTime)

  #undef profile
  #undef percent
  #undef flag_save
  #undef minutes
  #undef IntermediateCnts
  #undef StartTime
  #undef EndTime
  #undef CurrentTime
  /* component Origin=Progress_bar() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(_Origin->_rotation_absolute,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_copy(_Origin->_rotation_relative, _Origin->_rotation_absolute);
    _Origin->_rotation_is_identity =  rot_test_identity(_Origin->_rotation_relative);
    _Origin->_position_absolute = coords_set(
      0, 0, 0);
    tc1 = coords_neg(_Origin->_position_absolute);
    _Origin->_position_relative = rot_apply(_Origin->_rotation_absolute, tc1);
  } /* Origin=Progress_bar() AT ROTATED */
  DEBUG_COMPONENT("Origin", _Origin->_position_absolute, _Origin->_rotation_absolute);
  instrument->_position_absolute[1] = _Origin->_position_absolute;
  instrument->_position_relative[1] = _Origin->_position_relative;
  instrument->counter_N[1]  = instrument->counter_P[1] = instrument->counter_P2[1] = 0;
  instrument->counter_AbsorbProp[1]= 0;
  return(0);
} /* _Origin_setpos */

/* component Source=Source_simple() SETTING, POSITION/ROTATION */
int _Source_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Source_setpos] component Source=Source_simple() SETTING [/usr/share/mcstas/3.0-dev/sources/Source_simple.comp:64]");
  stracpy(_Source->_name, "Source", 16384);
  stracpy(_Source->_type, "Source_simple", 16384);
  _Source->_index=2;
  _Source->_parameters.radius = 0.01;
  #define radius (_Source->_parameters.radius)
  _Source->_parameters.yheight = 0;
  #define yheight (_Source->_parameters.yheight)
  _Source->_parameters.xwidth = 0;
  #define xwidth (_Source->_parameters.xwidth)
  _Source->_parameters.dist = 1.0;
  #define dist (_Source->_parameters.dist)
  _Source->_parameters.focus_xw = 0.02;
  #define focus_xw (_Source->_parameters.focus_xw)
  _Source->_parameters.focus_yh = 0.02;
  #define focus_yh (_Source->_parameters.focus_yh)
  _Source->_parameters.E0 = 0;
  #define E0 (_Source->_parameters.E0)
  _Source->_parameters.dE = 0;
  #define dE (_Source->_parameters.dE)
  _Source->_parameters.lambda0 = instrument->_parameters._lambda;
  #define lambda0 (_Source->_parameters.lambda0)
  _Source->_parameters.dlambda = 0.2 * instrument->_parameters._lambda;
  #define dlambda (_Source->_parameters.dlambda)
  _Source->_parameters.flux = 1;
  #define flux (_Source->_parameters.flux)
  _Source->_parameters.gauss = 0;
  #define gauss (_Source->_parameters.gauss)
  _Source->_parameters.target_index = + 1;
  #define target_index (_Source->_parameters.target_index)

  #define pmul (_Source->_parameters.pmul)
  #define square (_Source->_parameters.square)
  #define srcArea (_Source->_parameters.srcArea)

  #undef radius
  #undef yheight
  #undef xwidth
  #undef dist
  #undef focus_xw
  #undef focus_yh
  #undef E0
  #undef dE
  #undef lambda0
  #undef dlambda
  #undef flux
  #undef gauss
  #undef target_index
  #undef pmul
  #undef square
  #undef srcArea
  /* component Source=Source_simple() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Origin->_rotation_absolute, _Source->_rotation_absolute);
    rot_transpose(_Origin->_rotation_absolute, tr1);
    rot_mul(_Source->_rotation_absolute, tr1, _Source->_rotation_relative);
    _Source->_rotation_is_identity =  rot_test_identity(_Source->_rotation_relative);
    tc1 = coords_set(
      0, 0, 0);
    rot_transpose(_Origin->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Source->_position_absolute = coords_add(_Origin->_position_absolute, tc2);
    tc1 = coords_sub(_Origin->_position_absolute, _Source->_position_absolute);
    _Source->_position_relative = rot_apply(_Source->_rotation_absolute, tc1);
  } /* Source=Source_simple() AT ROTATED */
  DEBUG_COMPONENT("Source", _Source->_position_absolute, _Source->_rotation_absolute);
  instrument->_position_absolute[2] = _Source->_position_absolute;
  instrument->_position_relative[2] = _Source->_position_relative;
  instrument->counter_N[2]  = instrument->counter_P[2] = instrument->counter_P2[2] = 0;
  instrument->counter_AbsorbProp[2]= 0;
  return(0);
} /* _Source_setpos */

/* component lamStart=L_monitor() SETTING, POSITION/ROTATION */
int _lamStart_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_lamStart_setpos] component lamStart=L_monitor() SETTING [/usr/share/mcstas/3.0-dev/monitors/L_monitor.comp:65]");
  stracpy(_lamStart->_name, "lamStart", 16384);
  stracpy(_lamStart->_type, "L_monitor", 16384);
  _lamStart->_index=3;
  _lamStart->_parameters.nL = 200;
  #define nL (_lamStart->_parameters.nL)
  if("lambdaStart.dat" && strlen("lambdaStart.dat"))
    stracpy(_lamStart->_parameters.filename, "lambdaStart.dat" ? "lambdaStart.dat" : "", 16384);
  else 
  _lamStart->_parameters.filename[0]='\0';
  #define filename (_lamStart->_parameters.filename)
  _lamStart->_parameters.xmin = -0.05;
  #define xmin (_lamStart->_parameters.xmin)
  _lamStart->_parameters.xmax = 0.05;
  #define xmax (_lamStart->_parameters.xmax)
  _lamStart->_parameters.ymin = -0.05;
  #define ymin (_lamStart->_parameters.ymin)
  _lamStart->_parameters.ymax = 0.05;
  #define ymax (_lamStart->_parameters.ymax)
  _lamStart->_parameters.xwidth = 0.10;
  #define xwidth (_lamStart->_parameters.xwidth)
  _lamStart->_parameters.yheight = 0.10;
  #define yheight (_lamStart->_parameters.yheight)
  _lamStart->_parameters.Lmin = 0.7 * instrument->_parameters._lambda;
  #define Lmin (_lamStart->_parameters.Lmin)
  _lamStart->_parameters.Lmax = 1.3 * instrument->_parameters._lambda;
  #define Lmax (_lamStart->_parameters.Lmax)
  _lamStart->_parameters.restore_neutron = 0;
  #define restore_neutron (_lamStart->_parameters.restore_neutron)

  #define L_N (_lamStart->_parameters.L_N)
  #define L_p (_lamStart->_parameters.L_p)
  #define L_p2 (_lamStart->_parameters.L_p2)

  #undef nL
  #undef filename
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef xwidth
  #undef yheight
  #undef Lmin
  #undef Lmax
  #undef restore_neutron
  #undef L_N
  #undef L_p
  #undef L_p2
  /* component lamStart=L_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Origin->_rotation_absolute, _lamStart->_rotation_absolute);
    rot_transpose(_Source->_rotation_absolute, tr1);
    rot_mul(_lamStart->_rotation_absolute, tr1, _lamStart->_rotation_relative);
    _lamStart->_rotation_is_identity =  rot_test_identity(_lamStart->_rotation_relative);
    tc1 = coords_set(
      0, 0, 0.5);
    rot_transpose(_Origin->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _lamStart->_position_absolute = coords_add(_Origin->_position_absolute, tc2);
    tc1 = coords_sub(_Source->_position_absolute, _lamStart->_position_absolute);
    _lamStart->_position_relative = rot_apply(_lamStart->_rotation_absolute, tc1);
  } /* lamStart=L_monitor() AT ROTATED */
  DEBUG_COMPONENT("lamStart", _lamStart->_position_absolute, _lamStart->_rotation_absolute);
  instrument->_position_absolute[3] = _lamStart->_position_absolute;
  instrument->_position_relative[3] = _lamStart->_position_relative;
  instrument->counter_N[3]  = instrument->counter_P[3] = instrument->counter_P2[3] = 0;
  instrument->counter_AbsorbProp[3]= 0;
  return(0);
} /* _lamStart_setpos */

/* component Mono_Arm=Arm() SETTING, POSITION/ROTATION */
int _Mono_Arm_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Mono_Arm_setpos] component Mono_Arm=Arm() SETTING [Arm:0]");
  stracpy(_Mono_Arm->_name, "Mono_Arm", 16384);
  stracpy(_Mono_Arm->_type, "Arm", 16384);
  _Mono_Arm->_index=4;
  /* component Mono_Arm=Arm() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (A1)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _Source->_rotation_absolute, _Mono_Arm->_rotation_absolute);
    rot_transpose(_lamStart->_rotation_absolute, tr1);
    rot_mul(_Mono_Arm->_rotation_absolute, tr1, _Mono_Arm->_rotation_relative);
    _Mono_Arm->_rotation_is_identity =  rot_test_identity(_Mono_Arm->_rotation_relative);
    tc1 = coords_set(
      0, 0, 1.0);
    rot_transpose(_Source->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Mono_Arm->_position_absolute = coords_add(_Source->_position_absolute, tc2);
    tc1 = coords_sub(_lamStart->_position_absolute, _Mono_Arm->_position_absolute);
    _Mono_Arm->_position_relative = rot_apply(_Mono_Arm->_rotation_absolute, tc1);
  } /* Mono_Arm=Arm() AT ROTATED */
  DEBUG_COMPONENT("Mono_Arm", _Mono_Arm->_position_absolute, _Mono_Arm->_rotation_absolute);
  instrument->_position_absolute[4] = _Mono_Arm->_position_absolute;
  instrument->_position_relative[4] = _Mono_Arm->_position_relative;
  instrument->counter_N[4]  = instrument->counter_P[4] = instrument->counter_P2[4] = 0;
  instrument->counter_AbsorbProp[4]= 0;
  return(0);
} /* _Mono_Arm_setpos */

/* component Mono_Out=Arm() SETTING, POSITION/ROTATION */
int _Mono_Out_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Mono_Out_setpos] component Mono_Out=Arm() SETTING [Arm:0]");
  stracpy(_Mono_Out->_name, "Mono_Out", 16384);
  stracpy(_Mono_Out->_type, "Arm", 16384);
  _Mono_Out->_index=5;
  /* component Mono_Out=Arm() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (A2)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _Source->_rotation_absolute, _Mono_Out->_rotation_absolute);
    rot_transpose(_Mono_Arm->_rotation_absolute, tr1);
    rot_mul(_Mono_Out->_rotation_absolute, tr1, _Mono_Out->_rotation_relative);
    _Mono_Out->_rotation_is_identity =  rot_test_identity(_Mono_Out->_rotation_relative);
    tc1 = coords_set(
      0, 0, 0);
    rot_transpose(_Mono_Arm->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Mono_Out->_position_absolute = coords_add(_Mono_Arm->_position_absolute, tc2);
    tc1 = coords_sub(_Mono_Arm->_position_absolute, _Mono_Out->_position_absolute);
    _Mono_Out->_position_relative = rot_apply(_Mono_Out->_rotation_absolute, tc1);
  } /* Mono_Out=Arm() AT ROTATED */
  DEBUG_COMPONENT("Mono_Out", _Mono_Out->_position_absolute, _Mono_Out->_rotation_absolute);
  instrument->_position_absolute[5] = _Mono_Out->_position_absolute;
  instrument->_position_relative[5] = _Mono_Out->_position_relative;
  instrument->counter_N[5]  = instrument->counter_P[5] = instrument->counter_P2[5] = 0;
  instrument->counter_AbsorbProp[5]= 0;
  return(0);
} /* _Mono_Out_setpos */

/* component Mono1=Monochromator_flat() SETTING, POSITION/ROTATION */
int _Mono1_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Mono1_setpos] component Mono1=Monochromator_flat() SETTING [/usr/share/mcstas/3.0-dev/optics/Monochromator_flat.comp:102]");
  stracpy(_Mono1->_name, "Mono1", 16384);
  stracpy(_Mono1->_type, "Monochromator_flat", 16384);
  _Mono1->_index=6;
  _Mono1->_parameters.zmin = -0.05;
  #define zmin (_Mono1->_parameters.zmin)
  _Mono1->_parameters.zmax = 0.05;
  #define zmax (_Mono1->_parameters.zmax)
  _Mono1->_parameters.ymin = -0.05;
  #define ymin (_Mono1->_parameters.ymin)
  _Mono1->_parameters.ymax = 0.05;
  #define ymax (_Mono1->_parameters.ymax)
  _Mono1->_parameters.zwidth = 0.10;
  #define zwidth (_Mono1->_parameters.zwidth)
  _Mono1->_parameters.yheight = 0.10;
  #define yheight (_Mono1->_parameters.yheight)
  _Mono1->_parameters.mosaich = instrument->_parameters._Moz;
  #define mosaich (_Mono1->_parameters.mosaich)
  _Mono1->_parameters.mosaicv = instrument->_parameters._Moz;
  #define mosaicv (_Mono1->_parameters.mosaicv)
  _Mono1->_parameters.r0 = 0.9;
  #define r0 (_Mono1->_parameters.r0)
  _Mono1->_parameters.Q = mono_q;
  #define Q (_Mono1->_parameters.Q)
  _Mono1->_parameters.DM = 0;
  #define DM (_Mono1->_parameters.DM)

  #define mos_rms_y (_Mono1->_parameters.mos_rms_y)
  #define mos_rms_z (_Mono1->_parameters.mos_rms_z)
  #define mos_rms_max (_Mono1->_parameters.mos_rms_max)
  #define mono_Q (_Mono1->_parameters.mono_Q)

  #undef zmin
  #undef zmax
  #undef ymin
  #undef ymax
  #undef zwidth
  #undef yheight
  #undef mosaich
  #undef mosaicv
  #undef r0
  #undef Q
  #undef DM
  #undef mos_rms_y
  #undef mos_rms_z
  #undef mos_rms_max
  #undef mono_Q
  /* component Mono1=Monochromator_flat() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Mono_Arm->_rotation_absolute, _Mono1->_rotation_absolute);
    rot_transpose(_Mono_Out->_rotation_absolute, tr1);
    rot_mul(_Mono1->_rotation_absolute, tr1, _Mono1->_rotation_relative);
    _Mono1->_rotation_is_identity =  rot_test_identity(_Mono1->_rotation_relative);
    tc1 = coords_set(
      0, 0, 0);
    rot_transpose(_Mono_Arm->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Mono1->_position_absolute = coords_add(_Mono_Arm->_position_absolute, tc2);
    tc1 = coords_sub(_Mono_Out->_position_absolute, _Mono1->_position_absolute);
    _Mono1->_position_relative = rot_apply(_Mono1->_rotation_absolute, tc1);
  } /* Mono1=Monochromator_flat() AT ROTATED */
  DEBUG_COMPONENT("Mono1", _Mono1->_position_absolute, _Mono1->_rotation_absolute);
  instrument->_position_absolute[6] = _Mono1->_position_absolute;
  instrument->_position_relative[6] = _Mono1->_position_relative;
  instrument->counter_N[6]  = instrument->counter_P[6] = instrument->counter_P2[6] = 0;
  instrument->counter_AbsorbProp[6]= 0;
  return(0);
} /* _Mono1_setpos */

/* component Mono2=Monochromator_pol() SETTING, POSITION/ROTATION */
int _Mono2_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Mono2_setpos] component Mono2=Monochromator_pol() SETTING [/usr/share/mcstas/3.0-dev/optics/Monochromator_pol.comp:109]");
  stracpy(_Mono2->_name, "Mono2", 16384);
  stracpy(_Mono2->_type, "Monochromator_pol", 16384);
  _Mono2->_index=7;
  _Mono2->_parameters.zwidth = 0.10;
  #define zwidth (_Mono2->_parameters.zwidth)
  _Mono2->_parameters.yheight = 0.10;
  #define yheight (_Mono2->_parameters.yheight)
  _Mono2->_parameters.mosaic = instrument->_parameters._Moz;
  #define mosaic (_Mono2->_parameters.mosaic)
  _Mono2->_parameters.dspread = 0.0;
  #define dspread (_Mono2->_parameters.dspread)
  _Mono2->_parameters.Q = mono_q;
  #define Q (_Mono2->_parameters.Q)
  _Mono2->_parameters.DM = 0;
  #define DM (_Mono2->_parameters.DM)
  _Mono2->_parameters.pThreshold = 0;
  #define pThreshold (_Mono2->_parameters.pThreshold)
  _Mono2->_parameters.Rup = 0.9;
  #define Rup (_Mono2->_parameters.Rup)
  _Mono2->_parameters.Rdown = 0.9;
  #define Rdown (_Mono2->_parameters.Rdown)
  _Mono2->_parameters.debug = 0;
  #define debug (_Mono2->_parameters.debug)

  #define mos_rms (_Mono2->_parameters.mos_rms)
  #define d_rms (_Mono2->_parameters.d_rms)
  #define mono_Q (_Mono2->_parameters.mono_Q)
  #define FN (_Mono2->_parameters.FN)
  #define FM (_Mono2->_parameters.FM)

  #undef zwidth
  #undef yheight
  #undef mosaic
  #undef dspread
  #undef Q
  #undef DM
  #undef pThreshold
  #undef Rup
  #undef Rdown
  #undef debug
  #undef mos_rms
  #undef d_rms
  #undef mono_Q
  #undef FN
  #undef FM
  /* component Mono2=Monochromator_pol() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Mono_Arm->_rotation_absolute, _Mono2->_rotation_absolute);
    rot_transpose(_Mono1->_rotation_absolute, tr1);
    rot_mul(_Mono2->_rotation_absolute, tr1, _Mono2->_rotation_relative);
    _Mono2->_rotation_is_identity =  rot_test_identity(_Mono2->_rotation_relative);
    tc1 = coords_set(
      0, 0, 0);
    rot_transpose(_Mono_Arm->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Mono2->_position_absolute = coords_add(_Mono_Arm->_position_absolute, tc2);
    tc1 = coords_sub(_Mono1->_position_absolute, _Mono2->_position_absolute);
    _Mono2->_position_relative = rot_apply(_Mono2->_rotation_absolute, tc1);
  } /* Mono2=Monochromator_pol() AT ROTATED */
  DEBUG_COMPONENT("Mono2", _Mono2->_position_absolute, _Mono2->_rotation_absolute);
  instrument->_position_absolute[7] = _Mono2->_position_absolute;
  instrument->_position_relative[7] = _Mono2->_position_relative;
  instrument->counter_N[7]  = instrument->counter_P[7] = instrument->counter_P2[7] = 0;
  instrument->counter_AbsorbProp[7]= 0;
  return(0);
} /* _Mono2_setpos */

/* component Mono3=Monochromator_pol() SETTING, POSITION/ROTATION */
int _Mono3_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Mono3_setpos] component Mono3=Monochromator_pol() SETTING [/usr/share/mcstas/3.0-dev/optics/Monochromator_pol.comp:109]");
  stracpy(_Mono3->_name, "Mono3", 16384);
  stracpy(_Mono3->_type, "Monochromator_pol", 16384);
  _Mono3->_index=8;
  _Mono3->_parameters.zwidth = 0.10;
  #define zwidth (_Mono3->_parameters.zwidth)
  _Mono3->_parameters.yheight = 0.10;
  #define yheight (_Mono3->_parameters.yheight)
  _Mono3->_parameters.mosaic = instrument->_parameters._Moz;
  #define mosaic (_Mono3->_parameters.mosaic)
  _Mono3->_parameters.dspread = 0.0;
  #define dspread (_Mono3->_parameters.dspread)
  _Mono3->_parameters.Q = mono_q;
  #define Q (_Mono3->_parameters.Q)
  _Mono3->_parameters.DM = 0;
  #define DM (_Mono3->_parameters.DM)
  _Mono3->_parameters.pThreshold = 0.01;
  #define pThreshold (_Mono3->_parameters.pThreshold)
  _Mono3->_parameters.Rup = 0.9;
  #define Rup (_Mono3->_parameters.Rup)
  _Mono3->_parameters.Rdown = 0.9;
  #define Rdown (_Mono3->_parameters.Rdown)
  _Mono3->_parameters.debug = 0;
  #define debug (_Mono3->_parameters.debug)

  #define mos_rms (_Mono3->_parameters.mos_rms)
  #define d_rms (_Mono3->_parameters.d_rms)
  #define mono_Q (_Mono3->_parameters.mono_Q)
  #define FN (_Mono3->_parameters.FN)
  #define FM (_Mono3->_parameters.FM)

  #undef zwidth
  #undef yheight
  #undef mosaic
  #undef dspread
  #undef Q
  #undef DM
  #undef pThreshold
  #undef Rup
  #undef Rdown
  #undef debug
  #undef mos_rms
  #undef d_rms
  #undef mono_Q
  #undef FN
  #undef FM
  /* component Mono3=Monochromator_pol() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Mono_Arm->_rotation_absolute, _Mono3->_rotation_absolute);
    rot_transpose(_Mono2->_rotation_absolute, tr1);
    rot_mul(_Mono3->_rotation_absolute, tr1, _Mono3->_rotation_relative);
    _Mono3->_rotation_is_identity =  rot_test_identity(_Mono3->_rotation_relative);
    tc1 = coords_set(
      0, 0, 0);
    rot_transpose(_Mono_Arm->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Mono3->_position_absolute = coords_add(_Mono_Arm->_position_absolute, tc2);
    tc1 = coords_sub(_Mono2->_position_absolute, _Mono3->_position_absolute);
    _Mono3->_position_relative = rot_apply(_Mono3->_rotation_absolute, tc1);
  } /* Mono3=Monochromator_pol() AT ROTATED */
  DEBUG_COMPONENT("Mono3", _Mono3->_position_absolute, _Mono3->_rotation_absolute);
  instrument->_position_absolute[8] = _Mono3->_position_absolute;
  instrument->_position_relative[8] = _Mono3->_position_relative;
  instrument->counter_N[8]  = instrument->counter_P[8] = instrument->counter_P2[8] = 0;
  instrument->counter_AbsorbProp[8]= 0;
  return(0);
} /* _Mono3_setpos */

/* component Mono4=Monochromator_curved() SETTING, POSITION/ROTATION */
int _Mono4_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Mono4_setpos] component Mono4=Monochromator_curved() SETTING [/usr/share/mcstas/3.0-dev/optics/Monochromator_curved.comp:148]");
  stracpy(_Mono4->_name, "Mono4", 16384);
  stracpy(_Mono4->_type, "Monochromator_curved", 16384);
  _Mono4->_index=9;
  if("NULL" && strlen("NULL"))
    stracpy(_Mono4->_parameters.reflect, "NULL" ? "NULL" : "", 16384);
  else 
  _Mono4->_parameters.reflect[0]='\0';
  #define reflect (_Mono4->_parameters.reflect)
  if("NULL" && strlen("NULL"))
    stracpy(_Mono4->_parameters.transmit, "NULL" ? "NULL" : "", 16384);
  else 
  _Mono4->_parameters.transmit[0]='\0';
  #define transmit (_Mono4->_parameters.transmit)
  _Mono4->_parameters.zwidth = 0.01;
  #define zwidth (_Mono4->_parameters.zwidth)
  _Mono4->_parameters.yheight = 0.01;
  #define yheight (_Mono4->_parameters.yheight)
  _Mono4->_parameters.gap = 0;
  #define gap (_Mono4->_parameters.gap)
  _Mono4->_parameters.NH = 11;
  #define NH (_Mono4->_parameters.NH)
  _Mono4->_parameters.NV = 11;
  #define NV (_Mono4->_parameters.NV)
  _Mono4->_parameters.mosaich = instrument->_parameters._Moz;
  #define mosaich (_Mono4->_parameters.mosaich)
  _Mono4->_parameters.mosaicv = instrument->_parameters._Moz;
  #define mosaicv (_Mono4->_parameters.mosaicv)
  _Mono4->_parameters.r0 = 0.9;
  #define r0 (_Mono4->_parameters.r0)
  _Mono4->_parameters.t0 = 1.0;
  #define t0 (_Mono4->_parameters.t0)
  _Mono4->_parameters.Q = mono_q;
  #define Q (_Mono4->_parameters.Q)
  _Mono4->_parameters.RV = 0;
  #define RV (_Mono4->_parameters.RV)
  _Mono4->_parameters.RH = 0;
  #define RH (_Mono4->_parameters.RH)
  _Mono4->_parameters.DM = 0;
  #define DM (_Mono4->_parameters.DM)
  _Mono4->_parameters.mosaic = 0;
  #define mosaic (_Mono4->_parameters.mosaic)
  _Mono4->_parameters.width = 0.10;
  #define width (_Mono4->_parameters.width)
  _Mono4->_parameters.height = 0.10;
  #define height (_Mono4->_parameters.height)
  _Mono4->_parameters.verbose = 0;
  #define verbose (_Mono4->_parameters.verbose)
  _Mono4->_parameters.order = 0;
  #define order (_Mono4->_parameters.order)

  #define mos_rms_y (_Mono4->_parameters.mos_rms_y)
  #define mos_rms_z (_Mono4->_parameters.mos_rms_z)
  #define mos_rms_max (_Mono4->_parameters.mos_rms_max)
  #define mono_Q (_Mono4->_parameters.mono_Q)
  #define SlabWidth (_Mono4->_parameters.SlabWidth)
  #define SlabHeight (_Mono4->_parameters.SlabHeight)
  #define rTable (_Mono4->_parameters.rTable)
  #define tTable (_Mono4->_parameters.tTable)
  #define row (_Mono4->_parameters.row)
  #define col (_Mono4->_parameters.col)
  #define tiltH (_Mono4->_parameters.tiltH)
  #define tiltV (_Mono4->_parameters.tiltV)

  #undef reflect
  #undef transmit
  #undef zwidth
  #undef yheight
  #undef gap
  #undef NH
  #undef NV
  #undef mosaich
  #undef mosaicv
  #undef r0
  #undef t0
  #undef Q
  #undef RV
  #undef RH
  #undef DM
  #undef mosaic
  #undef width
  #undef height
  #undef verbose
  #undef order
  #undef mos_rms_y
  #undef mos_rms_z
  #undef mos_rms_max
  #undef mono_Q
  #undef SlabWidth
  #undef SlabHeight
  #undef rTable
  #undef tTable
  #undef row
  #undef col
  #undef tiltH
  #undef tiltV
  /* component Mono4=Monochromator_curved() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Mono_Arm->_rotation_absolute, _Mono4->_rotation_absolute);
    rot_transpose(_Mono3->_rotation_absolute, tr1);
    rot_mul(_Mono4->_rotation_absolute, tr1, _Mono4->_rotation_relative);
    _Mono4->_rotation_is_identity =  rot_test_identity(_Mono4->_rotation_relative);
    tc1 = coords_set(
      0, 0, 0);
    rot_transpose(_Mono_Arm->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Mono4->_position_absolute = coords_add(_Mono_Arm->_position_absolute, tc2);
    tc1 = coords_sub(_Mono3->_position_absolute, _Mono4->_position_absolute);
    _Mono4->_position_relative = rot_apply(_Mono4->_rotation_absolute, tc1);
  } /* Mono4=Monochromator_curved() AT ROTATED */
  DEBUG_COMPONENT("Mono4", _Mono4->_position_absolute, _Mono4->_rotation_absolute);
  instrument->_position_absolute[9] = _Mono4->_position_absolute;
  instrument->_position_relative[9] = _Mono4->_position_relative;
  instrument->counter_N[9]  = instrument->counter_P[9] = instrument->counter_P2[9] = 0;
  instrument->counter_AbsorbProp[9]= 0;
  return(0);
} /* _Mono4_setpos */

/* component Mono5=Monochromator_2foc() SETTING, POSITION/ROTATION */
int _Mono5_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Mono5_setpos] component Mono5=Monochromator_2foc() SETTING [/usr/share/mcstas/3.0-dev/contrib/Monochromator_2foc.comp:100]");
  stracpy(_Mono5->_name, "Mono5", 16384);
  stracpy(_Mono5->_type, "Monochromator_2foc", 16384);
  _Mono5->_index=10;
  _Mono5->_parameters.reflect[0]='\0';
  #define reflect (_Mono5->_parameters.reflect)
  _Mono5->_parameters.zwidth = 0.01;
  #define zwidth (_Mono5->_parameters.zwidth)
  _Mono5->_parameters.yheight = 0.01;
  #define yheight (_Mono5->_parameters.yheight)
  _Mono5->_parameters.gap = 0;
  #define gap (_Mono5->_parameters.gap)
  _Mono5->_parameters.NH = 11;
  #define NH (_Mono5->_parameters.NH)
  _Mono5->_parameters.NV = 11;
  #define NV (_Mono5->_parameters.NV)
  _Mono5->_parameters.mosaich = instrument->_parameters._Moz;
  #define mosaich (_Mono5->_parameters.mosaich)
  _Mono5->_parameters.mosaicv = instrument->_parameters._Moz;
  #define mosaicv (_Mono5->_parameters.mosaicv)
  _Mono5->_parameters.r0 = 0.9;
  #define r0 (_Mono5->_parameters.r0)
  _Mono5->_parameters.Q = mono_q;
  #define Q (_Mono5->_parameters.Q)
  _Mono5->_parameters.RV = 0;
  #define RV (_Mono5->_parameters.RV)
  _Mono5->_parameters.RH = 0;
  #define RH (_Mono5->_parameters.RH)
  _Mono5->_parameters.DM = 0;
  #define DM (_Mono5->_parameters.DM)
  _Mono5->_parameters.mosaic = 0;
  #define mosaic (_Mono5->_parameters.mosaic)
  _Mono5->_parameters.width = 0.10;
  #define width (_Mono5->_parameters.width)
  _Mono5->_parameters.height = 0.10;
  #define height (_Mono5->_parameters.height)
  _Mono5->_parameters.verbose = 0;
  #define verbose (_Mono5->_parameters.verbose)

  #define mos_y (_Mono5->_parameters.mos_y)
  #define mos_z (_Mono5->_parameters.mos_z)
  #define mono_Q (_Mono5->_parameters.mono_Q)
  #define SlabWidth (_Mono5->_parameters.SlabWidth)
  #define SlabHeight (_Mono5->_parameters.SlabHeight)
  #define rTable (_Mono5->_parameters.rTable)

  #undef reflect
  #undef zwidth
  #undef yheight
  #undef gap
  #undef NH
  #undef NV
  #undef mosaich
  #undef mosaicv
  #undef r0
  #undef Q
  #undef RV
  #undef RH
  #undef DM
  #undef mosaic
  #undef width
  #undef height
  #undef verbose
  #undef mos_y
  #undef mos_z
  #undef mono_Q
  #undef SlabWidth
  #undef SlabHeight
  #undef rTable
  /* component Mono5=Monochromator_2foc() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Mono_Arm->_rotation_absolute, _Mono5->_rotation_absolute);
    rot_transpose(_Mono4->_rotation_absolute, tr1);
    rot_mul(_Mono5->_rotation_absolute, tr1, _Mono5->_rotation_relative);
    _Mono5->_rotation_is_identity =  rot_test_identity(_Mono5->_rotation_relative);
    tc1 = coords_set(
      0, 0, 0);
    rot_transpose(_Mono_Arm->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Mono5->_position_absolute = coords_add(_Mono_Arm->_position_absolute, tc2);
    tc1 = coords_sub(_Mono4->_position_absolute, _Mono5->_position_absolute);
    _Mono5->_position_relative = rot_apply(_Mono5->_rotation_absolute, tc1);
  } /* Mono5=Monochromator_2foc() AT ROTATED */
  DEBUG_COMPONENT("Mono5", _Mono5->_position_absolute, _Mono5->_rotation_absolute);
  instrument->_position_absolute[10] = _Mono5->_position_absolute;
  instrument->_position_relative[10] = _Mono5->_position_relative;
  instrument->counter_N[10]  = instrument->counter_P[10] = instrument->counter_P2[10] = 0;
  instrument->counter_AbsorbProp[10]= 0;
  return(0);
} /* _Mono5_setpos */

/* component Mono6=Single_crystal() SETTING, POSITION/ROTATION */
int _Mono6_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Mono6_setpos] component Mono6=Single_crystal() SETTING [/usr/share/mcstas/3.0-dev/samples/Single_crystal.comp:1045]");
  stracpy(_Mono6->_name, "Mono6", 16384);
  stracpy(_Mono6->_type, "Single_crystal", 16384);
  _Mono6->_index=11;
  if("C_graphite.lau" && strlen("C_graphite.lau"))
    stracpy(_Mono6->_parameters.reflections, "C_graphite.lau" ? "C_graphite.lau" : "", 16384);
  else 
  _Mono6->_parameters.reflections[0]='\0';
  #define reflections (_Mono6->_parameters.reflections)
  _Mono6->_parameters.geometry[0]='\0';
  #define geometry (_Mono6->_parameters.geometry)
  _Mono6->_parameters.mosaic_AB = (double []){ 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 }; // static pointer allocation
  #define mosaic_AB (_Mono6->_parameters.mosaic_AB)
  _Mono6->_parameters.xwidth = 0.002;
  #define xwidth (_Mono6->_parameters.xwidth)
  _Mono6->_parameters.yheight = 0.1;
  #define yheight (_Mono6->_parameters.yheight)
  _Mono6->_parameters.zdepth = 0.1;
  #define zdepth (_Mono6->_parameters.zdepth)
  _Mono6->_parameters.radius = 0;
  #define radius (_Mono6->_parameters.radius)
  _Mono6->_parameters.delta_d_d = 1e-2;
  #define delta_d_d (_Mono6->_parameters.delta_d_d)
  _Mono6->_parameters.mosaic = instrument->_parameters._Moz;
  #define mosaic (_Mono6->_parameters.mosaic)
  _Mono6->_parameters.mosaic_a = -1;
  #define mosaic_a (_Mono6->_parameters.mosaic_a)
  _Mono6->_parameters.mosaic_b = -1;
  #define mosaic_b (_Mono6->_parameters.mosaic_b)
  _Mono6->_parameters.mosaic_c = -1;
  #define mosaic_c (_Mono6->_parameters.mosaic_c)
  _Mono6->_parameters.recip_cell = 0;
  #define recip_cell (_Mono6->_parameters.recip_cell)
  _Mono6->_parameters.barns = 0;
  #define barns (_Mono6->_parameters.barns)
  _Mono6->_parameters.ax = 0;
  #define ax (_Mono6->_parameters.ax)
  _Mono6->_parameters.ay = 2.14;
  #define ay (_Mono6->_parameters.ay)
  _Mono6->_parameters.az = -1.24;
  #define az (_Mono6->_parameters.az)
  _Mono6->_parameters.bx = 0;
  #define bx (_Mono6->_parameters.bx)
  _Mono6->_parameters.by = 0;
  #define by (_Mono6->_parameters.by)
  _Mono6->_parameters.bz = 2.47;
  #define bz (_Mono6->_parameters.bz)
  _Mono6->_parameters.cx = 6.71;
  #define cx (_Mono6->_parameters.cx)
  _Mono6->_parameters.cy = 0;
  #define cy (_Mono6->_parameters.cy)
  _Mono6->_parameters.cz = 0;
  #define cz (_Mono6->_parameters.cz)
  _Mono6->_parameters.p_transmit = -1;
  #define p_transmit (_Mono6->_parameters.p_transmit)
  _Mono6->_parameters.sigma_abs = 0.014;
  #define sigma_abs (_Mono6->_parameters.sigma_abs)
  _Mono6->_parameters.sigma_inc = 0.004;
  #define sigma_inc (_Mono6->_parameters.sigma_inc)
  _Mono6->_parameters.aa = 0;
  #define aa (_Mono6->_parameters.aa)
  _Mono6->_parameters.bb = 0;
  #define bb (_Mono6->_parameters.bb)
  _Mono6->_parameters.cc = 0;
  #define cc (_Mono6->_parameters.cc)
  _Mono6->_parameters.order = instrument->_parameters._order;
  #define order (_Mono6->_parameters.order)
  _Mono6->_parameters.RX = 0;
  #define RX (_Mono6->_parameters.RX)
  _Mono6->_parameters.RY = 0;
  #define RY (_Mono6->_parameters.RY)
  _Mono6->_parameters.RZ = 0;
  #define RZ (_Mono6->_parameters.RZ)
  _Mono6->_parameters.powder = instrument->_parameters._powder;
  #define powder (_Mono6->_parameters.powder)
  _Mono6->_parameters.PG = instrument->_parameters._PG;
  #define PG (_Mono6->_parameters.PG)

  #define hkl_info (_Mono6->_parameters.hkl_info)
  #define offdata (_Mono6->_parameters.offdata)

  #undef reflections
  #undef geometry
  #undef mosaic_AB
  #undef xwidth
  #undef yheight
  #undef zdepth
  #undef radius
  #undef delta_d_d
  #undef mosaic
  #undef mosaic_a
  #undef mosaic_b
  #undef mosaic_c
  #undef recip_cell
  #undef barns
  #undef ax
  #undef ay
  #undef az
  #undef bx
  #undef by
  #undef bz
  #undef cx
  #undef cy
  #undef cz
  #undef p_transmit
  #undef sigma_abs
  #undef sigma_inc
  #undef aa
  #undef bb
  #undef cc
  #undef order
  #undef RX
  #undef RY
  #undef RZ
  #undef powder
  #undef PG
  #undef hkl_info
  #undef offdata
  /* component Mono6=Single_crystal() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Mono_Arm->_rotation_absolute, _Mono6->_rotation_absolute);
    rot_transpose(_Mono5->_rotation_absolute, tr1);
    rot_mul(_Mono6->_rotation_absolute, tr1, _Mono6->_rotation_relative);
    _Mono6->_rotation_is_identity =  rot_test_identity(_Mono6->_rotation_relative);
    tc1 = coords_set(
      0, 0, 0);
    rot_transpose(_Mono_Arm->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Mono6->_position_absolute = coords_add(_Mono_Arm->_position_absolute, tc2);
    tc1 = coords_sub(_Mono5->_position_absolute, _Mono6->_position_absolute);
    _Mono6->_position_relative = rot_apply(_Mono6->_rotation_absolute, tc1);
  } /* Mono6=Single_crystal() AT ROTATED */
  DEBUG_COMPONENT("Mono6", _Mono6->_position_absolute, _Mono6->_rotation_absolute);
  instrument->_position_absolute[11] = _Mono6->_position_absolute;
  instrument->_position_relative[11] = _Mono6->_position_relative;
  instrument->counter_N[11]  = instrument->counter_P[11] = instrument->counter_P2[11] = 0;
  instrument->counter_AbsorbProp[11]= 0;
  return(0);
} /* _Mono6_setpos */

/* component Sphere1=PSD_monitor_4PI() SETTING, POSITION/ROTATION */
int _Sphere1_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Sphere1_setpos] component Sphere1=PSD_monitor_4PI() SETTING [/usr/share/mcstas/3.0-dev/monitors/PSD_monitor_4PI.comp:62]");
  stracpy(_Sphere1->_name, "Sphere1", 16384);
  stracpy(_Sphere1->_type, "PSD_monitor_4PI", 16384);
  _Sphere1->_index=12;
  _Sphere1->_parameters.nx = 90;
  #define nx (_Sphere1->_parameters.nx)
  _Sphere1->_parameters.ny = 90;
  #define ny (_Sphere1->_parameters.ny)
  if("sphere" && strlen("sphere"))
    stracpy(_Sphere1->_parameters.filename, "sphere" ? "sphere" : "", 16384);
  else 
  _Sphere1->_parameters.filename[0]='\0';
  #define filename (_Sphere1->_parameters.filename)
  _Sphere1->_parameters.radius = 1;
  #define radius (_Sphere1->_parameters.radius)
  _Sphere1->_parameters.restore_neutron = 1;
  #define restore_neutron (_Sphere1->_parameters.restore_neutron)

  #define PSD_N (_Sphere1->_parameters.PSD_N)
  #define PSD_p (_Sphere1->_parameters.PSD_p)
  #define PSD_p2 (_Sphere1->_parameters.PSD_p2)

  #undef nx
  #undef ny
  #undef filename
  #undef radius
  #undef restore_neutron
  #undef PSD_N
  #undef PSD_p
  #undef PSD_p2
  /* component Sphere1=PSD_monitor_4PI() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Mono_Out->_rotation_absolute, _Sphere1->_rotation_absolute);
    rot_transpose(_Mono6->_rotation_absolute, tr1);
    rot_mul(_Sphere1->_rotation_absolute, tr1, _Sphere1->_rotation_relative);
    _Sphere1->_rotation_is_identity =  rot_test_identity(_Sphere1->_rotation_relative);
    tc1 = coords_set(
      0, 0, 0);
    rot_transpose(_Mono_Out->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Sphere1->_position_absolute = coords_add(_Mono_Out->_position_absolute, tc2);
    tc1 = coords_sub(_Mono6->_position_absolute, _Sphere1->_position_absolute);
    _Sphere1->_position_relative = rot_apply(_Sphere1->_rotation_absolute, tc1);
  } /* Sphere1=PSD_monitor_4PI() AT ROTATED */
  DEBUG_COMPONENT("Sphere1", _Sphere1->_position_absolute, _Sphere1->_rotation_absolute);
  instrument->_position_absolute[12] = _Sphere1->_position_absolute;
  instrument->_position_relative[12] = _Sphere1->_position_relative;
  instrument->counter_N[12]  = instrument->counter_P[12] = instrument->counter_P2[12] = 0;
  instrument->counter_AbsorbProp[12]= 0;
  return(0);
} /* _Sphere1_setpos */

/* component lam1=L_monitor() SETTING, POSITION/ROTATION */
int _lam1_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_lam1_setpos] component lam1=L_monitor() SETTING [/usr/share/mcstas/3.0-dev/monitors/L_monitor.comp:65]");
  stracpy(_lam1->_name, "lam1", 16384);
  stracpy(_lam1->_type, "L_monitor", 16384);
  _lam1->_index=13;
  _lam1->_parameters.nL = 200;
  #define nL (_lam1->_parameters.nL)
  if("lambda1.dat" && strlen("lambda1.dat"))
    stracpy(_lam1->_parameters.filename, "lambda1.dat" ? "lambda1.dat" : "", 16384);
  else 
  _lam1->_parameters.filename[0]='\0';
  #define filename (_lam1->_parameters.filename)
  _lam1->_parameters.xmin = -0.05;
  #define xmin (_lam1->_parameters.xmin)
  _lam1->_parameters.xmax = 0.05;
  #define xmax (_lam1->_parameters.xmax)
  _lam1->_parameters.ymin = -0.05;
  #define ymin (_lam1->_parameters.ymin)
  _lam1->_parameters.ymax = 0.05;
  #define ymax (_lam1->_parameters.ymax)
  _lam1->_parameters.xwidth = 0.10;
  #define xwidth (_lam1->_parameters.xwidth)
  _lam1->_parameters.yheight = 0.10;
  #define yheight (_lam1->_parameters.yheight)
  _lam1->_parameters.Lmin = 0.7 * instrument->_parameters._lambda;
  #define Lmin (_lam1->_parameters.Lmin)
  _lam1->_parameters.Lmax = 1.3 * instrument->_parameters._lambda;
  #define Lmax (_lam1->_parameters.Lmax)
  _lam1->_parameters.restore_neutron = 0;
  #define restore_neutron (_lam1->_parameters.restore_neutron)

  #define L_N (_lam1->_parameters.L_N)
  #define L_p (_lam1->_parameters.L_p)
  #define L_p2 (_lam1->_parameters.L_p2)

  #undef nL
  #undef filename
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef xwidth
  #undef yheight
  #undef Lmin
  #undef Lmax
  #undef restore_neutron
  #undef L_N
  #undef L_p
  #undef L_p2
  /* component lam1=L_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Mono_Out->_rotation_absolute, _lam1->_rotation_absolute);
    rot_transpose(_Sphere1->_rotation_absolute, tr1);
    rot_mul(_lam1->_rotation_absolute, tr1, _lam1->_rotation_relative);
    _lam1->_rotation_is_identity =  rot_test_identity(_lam1->_rotation_relative);
    tc1 = coords_set(
      0, 0, 0.25);
    rot_transpose(_Mono_Out->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _lam1->_position_absolute = coords_add(_Mono_Out->_position_absolute, tc2);
    tc1 = coords_sub(_Sphere1->_position_absolute, _lam1->_position_absolute);
    _lam1->_position_relative = rot_apply(_lam1->_rotation_absolute, tc1);
  } /* lam1=L_monitor() AT ROTATED */
  DEBUG_COMPONENT("lam1", _lam1->_position_absolute, _lam1->_rotation_absolute);
  instrument->_position_absolute[13] = _lam1->_position_absolute;
  instrument->_position_relative[13] = _lam1->_position_relative;
  instrument->counter_N[13]  = instrument->counter_P[13] = instrument->counter_P2[13] = 0;
  instrument->counter_AbsorbProp[13]= 0;
  return(0);
} /* _lam1_setpos */

/* component psd1=PSD_monitor() SETTING, POSITION/ROTATION */
int _psd1_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_psd1_setpos] component psd1=PSD_monitor() SETTING [/usr/share/mcstas/3.0-dev/monitors/PSD_monitor.comp:64]");
  stracpy(_psd1->_name, "psd1", 16384);
  stracpy(_psd1->_type, "PSD_monitor", 16384);
  _psd1->_index=14;
  _psd1->_parameters.nx = 20;
  #define nx (_psd1->_parameters.nx)
  _psd1->_parameters.ny = 20;
  #define ny (_psd1->_parameters.ny)
  if("psd1.dat" && strlen("psd1.dat"))
    stracpy(_psd1->_parameters.filename, "psd1.dat" ? "psd1.dat" : "", 16384);
  else 
  _psd1->_parameters.filename[0]='\0';
  #define filename (_psd1->_parameters.filename)
  _psd1->_parameters.xmin = -0.05;
  #define xmin (_psd1->_parameters.xmin)
  _psd1->_parameters.xmax = 0.05;
  #define xmax (_psd1->_parameters.xmax)
  _psd1->_parameters.ymin = -0.05;
  #define ymin (_psd1->_parameters.ymin)
  _psd1->_parameters.ymax = 0.05;
  #define ymax (_psd1->_parameters.ymax)
  _psd1->_parameters.xwidth = 0.10;
  #define xwidth (_psd1->_parameters.xwidth)
  _psd1->_parameters.yheight = 0.10;
  #define yheight (_psd1->_parameters.yheight)
  _psd1->_parameters.restore_neutron = 0;
  #define restore_neutron (_psd1->_parameters.restore_neutron)

  #define PSD_N (_psd1->_parameters.PSD_N)
  #define PSD_p (_psd1->_parameters.PSD_p)
  #define PSD_p2 (_psd1->_parameters.PSD_p2)

  #undef nx
  #undef ny
  #undef filename
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef xwidth
  #undef yheight
  #undef restore_neutron
  #undef PSD_N
  #undef PSD_p
  #undef PSD_p2
  /* component psd1=PSD_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Mono_Out->_rotation_absolute, _psd1->_rotation_absolute);
    rot_transpose(_lam1->_rotation_absolute, tr1);
    rot_mul(_psd1->_rotation_absolute, tr1, _psd1->_rotation_relative);
    _psd1->_rotation_is_identity =  rot_test_identity(_psd1->_rotation_relative);
    tc1 = coords_set(
      0, 0, 0.5);
    rot_transpose(_Mono_Out->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _psd1->_position_absolute = coords_add(_Mono_Out->_position_absolute, tc2);
    tc1 = coords_sub(_lam1->_position_absolute, _psd1->_position_absolute);
    _psd1->_position_relative = rot_apply(_psd1->_rotation_absolute, tc1);
  } /* psd1=PSD_monitor() AT ROTATED */
  DEBUG_COMPONENT("psd1", _psd1->_position_absolute, _psd1->_rotation_absolute);
  instrument->_position_absolute[14] = _psd1->_position_absolute;
  instrument->_position_relative[14] = _psd1->_position_relative;
  instrument->counter_N[14]  = instrument->counter_P[14] = instrument->counter_P2[14] = 0;
  instrument->counter_AbsorbProp[14]= 0;
  return(0);
} /* _psd1_setpos */

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

_class_Source_simple *class_Source_simple_init(_class_Source_simple *_comp
) {
  #define radius (_comp->_parameters.radius)
  #define yheight (_comp->_parameters.yheight)
  #define xwidth (_comp->_parameters.xwidth)
  #define dist (_comp->_parameters.dist)
  #define focus_xw (_comp->_parameters.focus_xw)
  #define focus_yh (_comp->_parameters.focus_yh)
  #define E0 (_comp->_parameters.E0)
  #define dE (_comp->_parameters.dE)
  #define lambda0 (_comp->_parameters.lambda0)
  #define dlambda (_comp->_parameters.dlambda)
  #define flux (_comp->_parameters.flux)
  #define gauss (_comp->_parameters.gauss)
  #define target_index (_comp->_parameters.target_index)
  #define pmul (_comp->_parameters.pmul)
  #define square (_comp->_parameters.square)
  #define srcArea (_comp->_parameters.srcArea)
square = 0;
/* Determine source area */
if (radius && !yheight && !xwidth ) {
    square = 0;
    srcArea = PI*radius*radius;
  } else if(yheight && xwidth) {
    square = 1;
    srcArea = xwidth * yheight;
  }

  if (flux) {
    pmul=flux*1e4*srcArea/mcget_ncount();
    if (dlambda)
      pmul *= 2*dlambda;
    else if (dE)
      pmul *= 2*dE;
  } else {
    gauss = 0;
    pmul=1.0/(mcget_ncount()*4*PI);
  }

  if (target_index && !dist)
  {
    Coords ToTarget;
    double tx,ty,tz;
    ToTarget = coords_sub(POS_A_COMP_INDEX(INDEX_CURRENT_COMP+target_index),POS_A_CURRENT_COMP);
    ToTarget = rot_apply(ROT_A_CURRENT_COMP, ToTarget);
    coords_get(ToTarget, &tx, &ty, &tz);
    dist=sqrt(tx*tx+ty*ty+tz*tz);
  }

  if (srcArea <= 0) {
    printf("Source_simple: %s: Source area is <= 0 !\n ERROR - Exiting\n",
           NAME_CURRENT_COMP);
    exit(0);
  }
  if (dist <= 0 || focus_xw <= 0 || focus_yh <= 0) {
    printf("Source_simple: %s: Target area unmeaningful! (negative dist / focus_xw / focus_yh)\n ERROR - Exiting\n",
           NAME_CURRENT_COMP);
    exit(0);
  }

  if ((!lambda0 && !E0 && !dE && !dlambda)) {
    printf("Source_simple: %s: You must specify either a wavelength or energy range!\n ERROR - Exiting\n",
           NAME_CURRENT_COMP);
    exit(0);
  }
  if ((!lambda0 && !dlambda && (E0 <= 0 || dE < 0 || E0-dE <= 0))
    || (!E0 && !dE && (lambda0 <= 0 || dlambda < 0 || lambda0-dlambda <= 0))) {
    printf("Source_simple: %s: Unmeaningful definition of wavelength or energy range!\n ERROR - Exiting\n",
           NAME_CURRENT_COMP);
      exit(0);
  }
  #undef radius
  #undef yheight
  #undef xwidth
  #undef dist
  #undef focus_xw
  #undef focus_yh
  #undef E0
  #undef dE
  #undef lambda0
  #undef dlambda
  #undef flux
  #undef gauss
  #undef target_index
  #undef pmul
  #undef square
  #undef srcArea
  return(_comp);
} /* class_Source_simple_init */

_class_L_monitor *class_L_monitor_init(_class_L_monitor *_comp
) {
  #define nL (_comp->_parameters.nL)
  #define filename (_comp->_parameters.filename)
  #define xmin (_comp->_parameters.xmin)
  #define xmax (_comp->_parameters.xmax)
  #define ymin (_comp->_parameters.ymin)
  #define ymax (_comp->_parameters.ymax)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define Lmin (_comp->_parameters.Lmin)
  #define Lmax (_comp->_parameters.Lmax)
  #define restore_neutron (_comp->_parameters.restore_neutron)
  #define L_N (_comp->_parameters.L_N)
  #define L_p (_comp->_parameters.L_p)
  #define L_p2 (_comp->_parameters.L_p2)
  if (xwidth  > 0) { xmax = xwidth/2;  xmin = -xmax; }
  if (yheight > 0) { ymax = yheight/2; ymin = -ymax; }

  if ((xmin >= xmax) || (ymin >= ymax)) {
    printf("L_monitor: %s: Null detection area !\n"
      "ERROR      (xwidth,yheight,xmin,xmax,ymin,ymax). Exiting",
      NAME_CURRENT_COMP);
    exit(0);
  }

  L_N = create_darr1d(nL);
  L_p = create_darr1d(nL);
  L_p2 = create_darr1d(nL);

  int i;
  for (i=0; i<nL; i++)
  {
    L_N[i] = 0;
    L_p[i] = 0;
    L_p2[i] = 0;
  }
  #undef nL
  #undef filename
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef xwidth
  #undef yheight
  #undef Lmin
  #undef Lmax
  #undef restore_neutron
  #undef L_N
  #undef L_p
  #undef L_p2
  return(_comp);
} /* class_L_monitor_init */

_class_Monochromator_flat *class_Monochromator_flat_init(_class_Monochromator_flat *_comp
) {
  #define zmin (_comp->_parameters.zmin)
  #define zmax (_comp->_parameters.zmax)
  #define ymin (_comp->_parameters.ymin)
  #define ymax (_comp->_parameters.ymax)
  #define zwidth (_comp->_parameters.zwidth)
  #define yheight (_comp->_parameters.yheight)
  #define mosaich (_comp->_parameters.mosaich)
  #define mosaicv (_comp->_parameters.mosaicv)
  #define r0 (_comp->_parameters.r0)
  #define Q (_comp->_parameters.Q)
  #define DM (_comp->_parameters.DM)
  #define mos_rms_y (_comp->_parameters.mos_rms_y)
  #define mos_rms_z (_comp->_parameters.mos_rms_z)
  #define mos_rms_max (_comp->_parameters.mos_rms_max)
  #define mono_Q (_comp->_parameters.mono_Q)
  mos_rms_y = MIN2RAD*mosaicv/sqrt(8*log(2));
  mos_rms_z = MIN2RAD*mosaich/sqrt(8*log(2));
  mos_rms_max = mos_rms_y > mos_rms_z ? mos_rms_y : mos_rms_z;

  mono_Q = Q;
  if (DM != 0) mono_Q = 2*PI/DM;

  if (zwidth>0)  { zmax = zwidth/2;  zmin=-zmax; }
  if (yheight>0) { ymax = yheight/2; ymin=-ymax; }

  if (zmin==zmax || ymin==ymax)
    exit(fprintf(stderr, "Monochromator_flat: %s : Surface is null (zmin,zmax,ymin,ymax)\n", NAME_CURRENT_COMP));
  #undef zmin
  #undef zmax
  #undef ymin
  #undef ymax
  #undef zwidth
  #undef yheight
  #undef mosaich
  #undef mosaicv
  #undef r0
  #undef Q
  #undef DM
  #undef mos_rms_y
  #undef mos_rms_z
  #undef mos_rms_max
  #undef mono_Q
  return(_comp);
} /* class_Monochromator_flat_init */

_class_Monochromator_pol *class_Monochromator_pol_init(_class_Monochromator_pol *_comp
) {
  #define zwidth (_comp->_parameters.zwidth)
  #define yheight (_comp->_parameters.yheight)
  #define mosaic (_comp->_parameters.mosaic)
  #define dspread (_comp->_parameters.dspread)
  #define Q (_comp->_parameters.Q)
  #define DM (_comp->_parameters.DM)
  #define pThreshold (_comp->_parameters.pThreshold)
  #define Rup (_comp->_parameters.Rup)
  #define Rdown (_comp->_parameters.Rdown)
  #define debug (_comp->_parameters.debug)
  #define mos_rms (_comp->_parameters.mos_rms)
  #define d_rms (_comp->_parameters.d_rms)
  #define mono_Q (_comp->_parameters.mono_Q)
  #define FN (_comp->_parameters.FN)
  #define FM (_comp->_parameters.FM)
mos_rms = MIN2RAD*mosaic/sqrt(8*log(2));

  mono_Q = Q;
  if (DM != 0)
    mono_Q = 2*PI/DM;

  DM = 2*PI/mono_Q;
  d_rms = dspread*DM/sqrt(8*log(2));

  // calculate the unit cell nuclear and magnetic structure factors
  if(debug > 0)
    printf("Rup: %f, Rdown: %f\n", Rup, Rdown);

  GetMonoPolFNFM(Rup, Rdown, &FN, &FM);

  if(debug > 0)
    printf("FN: %f, FM: %f\n", FN, FM);
  #undef zwidth
  #undef yheight
  #undef mosaic
  #undef dspread
  #undef Q
  #undef DM
  #undef pThreshold
  #undef Rup
  #undef Rdown
  #undef debug
  #undef mos_rms
  #undef d_rms
  #undef mono_Q
  #undef FN
  #undef FM
  return(_comp);
} /* class_Monochromator_pol_init */

_class_Monochromator_curved *class_Monochromator_curved_init(_class_Monochromator_curved *_comp
) {
  #define reflect (_comp->_parameters.reflect)
  #define transmit (_comp->_parameters.transmit)
  #define zwidth (_comp->_parameters.zwidth)
  #define yheight (_comp->_parameters.yheight)
  #define gap (_comp->_parameters.gap)
  #define NH (_comp->_parameters.NH)
  #define NV (_comp->_parameters.NV)
  #define mosaich (_comp->_parameters.mosaich)
  #define mosaicv (_comp->_parameters.mosaicv)
  #define r0 (_comp->_parameters.r0)
  #define t0 (_comp->_parameters.t0)
  #define Q (_comp->_parameters.Q)
  #define RV (_comp->_parameters.RV)
  #define RH (_comp->_parameters.RH)
  #define DM (_comp->_parameters.DM)
  #define mosaic (_comp->_parameters.mosaic)
  #define width (_comp->_parameters.width)
  #define height (_comp->_parameters.height)
  #define verbose (_comp->_parameters.verbose)
  #define order (_comp->_parameters.order)
  #define mos_rms_y (_comp->_parameters.mos_rms_y)
  #define mos_rms_z (_comp->_parameters.mos_rms_z)
  #define mos_rms_max (_comp->_parameters.mos_rms_max)
  #define mono_Q (_comp->_parameters.mono_Q)
  #define SlabWidth (_comp->_parameters.SlabWidth)
  #define SlabHeight (_comp->_parameters.SlabHeight)
  #define rTable (_comp->_parameters.rTable)
  #define tTable (_comp->_parameters.tTable)
  #define row (_comp->_parameters.row)
  #define col (_comp->_parameters.col)
  #define tiltH (_comp->_parameters.tiltH)
  #define tiltV (_comp->_parameters.tiltV)
  int i;

  if (mosaic != 0) {
    mos_rms_y = MIN2RAD*mosaic/sqrt(8*log(2));
    mos_rms_z = mos_rms_y; }
  else {
    mos_rms_y = MIN2RAD*mosaich/sqrt(8*log(2));
    mos_rms_z = MIN2RAD*mosaicv/sqrt(8*log(2)); }
  mos_rms_max = mos_rms_y > mos_rms_z ? mos_rms_y : mos_rms_z;

  mono_Q = Q;
  if (DM != 0) mono_Q = 2*PI/DM;

  if (mono_Q <= 0) { fprintf(stderr,"Monochromator_curved: %s: Error scattering vector Q = 0\n", NAME_CURRENT_COMP); exit(-1); }
  if (r0 <  0) { fprintf(stderr,"Monochromator_curved: %s: Error reflectivity r0 is negative\n", NAME_CURRENT_COMP); exit(-1); }
  if (r0 == 0) { fprintf(stderr,"Monochromator_curved: %s: Reflectivity r0 is null. Ignoring component.\n", NAME_CURRENT_COMP); }
  if (NH*NV == 0) { fprintf(stderr,"Monochromator_curved: %s: no slabs ??? (NH or NV=0)\n", NAME_CURRENT_COMP); exit(-1); }


  if (verbose && r0)
  {
    printf("Monochromator_curved: component %s Q=%.3g Angs-1 (DM=%.4g Angs)\n", NAME_CURRENT_COMP, mono_Q, 2*PI/mono_Q);
    if (NH*NV == 1) printf("            flat.\n");
    else
    { if (NH > 1)
      { printf("            horizontal: %i blades", (int)NH);
        if (RH != 0) printf(" focusing with RH=%.3g [m]", RH);
        printf("\n");
      }
      if (NV > 1)
      { printf("            vertical:   %i blades", (int)NV);
        if (RV != 0) printf(" focusing with RV=%.3g [m]", RV);
        printf("\n");
      }
    }
  }

  if (reflect != NULL && r0 && strlen(reflect) && strcmp(reflect,"NULL") && strcmp(reflect,"0"))
  {
    if (verbose) fprintf(stdout, "Monochromator_curved: %s: Reflectivity data (k, R) from %s\n", NAME_CURRENT_COMP, reflect);
    Table_Read(&rTable, reflect, 1); /* read 1st block data from file into rTable */
    Table_Rebin(&rTable);         /* rebin as evenly, increasing array */
    if (rTable.rows < 2) Table_Free(&rTable);
    if (verbose) Table_Info(rTable);
  } else rTable.data = NULL;

  if (transmit != NULL && strlen(transmit) && strcmp(transmit,"NULL") && strcmp(transmit,"0"))
  {
    if (verbose) fprintf(stdout, "Monochromator_curved: %s: Transmission data (k, T) from %s\n", NAME_CURRENT_COMP, transmit);
    Table_Read(&tTable, transmit, 1); /* read 1st block data from file into rTable */
    Table_Rebin(&tTable);         /* rebin as evenly, increasing array */
    if (tTable.rows < 2) Table_Free(&tTable);
    if (verbose) Table_Info(tTable);
  } else tTable.data = NULL;

  if (width == 0) SlabWidth = zwidth;
  else SlabWidth = (width+gap)/NH - gap;
  if (height == 0) SlabHeight = yheight;
  else SlabHeight = (height+gap)/NV - gap;

  tiltH=malloc((int)(NH+1)*sizeof(double));
  tiltV=malloc((int)(NV+1)*sizeof(double));

  if (!tiltH) printf("Monochromator_curved: %s: Warning: not enough memory to allocate tilts (NH=%g).\n", NAME_CURRENT_COMP, NH);
  else if (RH) { /* pre-compute tilts */
    for (i=0;i<=NH;i++)
    {
      tiltH[i]=asin((i-(NH+1)/2)*(SlabWidth+gap)/RH);
    }
  }
  if (!tiltV) printf("Monochromator_curved: %s: Warning: not enough memory to allocate tilts (NV=%g).\n", NAME_CURRENT_COMP, NV);
  else if (RV) {
    for (i=0;i<=NV;i++)
    {
      tiltV[i]=-asin((i-(NV+1)/2)*(SlabHeight+gap)/RV);
    }
  }

  #undef reflect
  #undef transmit
  #undef zwidth
  #undef yheight
  #undef gap
  #undef NH
  #undef NV
  #undef mosaich
  #undef mosaicv
  #undef r0
  #undef t0
  #undef Q
  #undef RV
  #undef RH
  #undef DM
  #undef mosaic
  #undef width
  #undef height
  #undef verbose
  #undef order
  #undef mos_rms_y
  #undef mos_rms_z
  #undef mos_rms_max
  #undef mono_Q
  #undef SlabWidth
  #undef SlabHeight
  #undef rTable
  #undef tTable
  #undef row
  #undef col
  #undef tiltH
  #undef tiltV
  return(_comp);
} /* class_Monochromator_curved_init */

_class_Monochromator_2foc *class_Monochromator_2foc_init(_class_Monochromator_2foc *_comp
) {
  #define reflect (_comp->_parameters.reflect)
  #define zwidth (_comp->_parameters.zwidth)
  #define yheight (_comp->_parameters.yheight)
  #define gap (_comp->_parameters.gap)
  #define NH (_comp->_parameters.NH)
  #define NV (_comp->_parameters.NV)
  #define mosaich (_comp->_parameters.mosaich)
  #define mosaicv (_comp->_parameters.mosaicv)
  #define r0 (_comp->_parameters.r0)
  #define Q (_comp->_parameters.Q)
  #define RV (_comp->_parameters.RV)
  #define RH (_comp->_parameters.RH)
  #define DM (_comp->_parameters.DM)
  #define mosaic (_comp->_parameters.mosaic)
  #define width (_comp->_parameters.width)
  #define height (_comp->_parameters.height)
  #define verbose (_comp->_parameters.verbose)
  #define mos_y (_comp->_parameters.mos_y)
  #define mos_z (_comp->_parameters.mos_z)
  #define mono_Q (_comp->_parameters.mono_Q)
  #define SlabWidth (_comp->_parameters.SlabWidth)
  #define SlabHeight (_comp->_parameters.SlabHeight)
  #define rTable (_comp->_parameters.rTable)
if (mosaic != 0) {
    mos_y = mosaic;
    mos_z = mos_y; }
  else {
    mos_y = mosaich;
    mos_z = mosaicv; }

  mono_Q = Q;
  if (DM != 0) mono_Q = 2*PI/DM;

  if (mono_Q == 0) { fprintf(stderr,"Monochromator_2foc: %s: Error scattering vector Q = 0\n", NAME_CURRENT_COMP); exit(-1); }
  if (r0 == 0) { fprintf(stderr,"Monochromator_2foc: %s: Error reflectivity r0 is null\n", NAME_CURRENT_COMP); exit(-1); }
  if (NH*NV == 0) { fprintf(stderr,"Monochromator_2foc: %s: no slabs ??? (NH or NV=0)\n", NAME_CURRENT_COMP); exit(-1); }

  if (verbose)
  {
    printf("Monochromator_2foc: component %s Q=%.3g Angs-1 (DM=%.4g Angs)\n", NAME_CURRENT_COMP, mono_Q, 2*PI/mono_Q);
    if (NH*NV == 1) printf("            flat.\n");
    else
    { if (NH > 1)
      { printf("            horizontal: %i blades", (int)NH);
        if (RH != 0) printf(" focusing with RH=%.3g [m]", RH);
        printf("\n");
      }
      if (NV > 1)
      { printf("            vertical:   %i blades", (int)NV);
        if (RV != 0) printf(" focusing with RV=%.3g [m]", RV);
        printf("\n");
      }
    }
  }

  if (reflect != NULL)
  {
    if (verbose) fprintf(stdout, "Monochromator_2foc: %s : Reflectivity data (k, R)\n", NAME_CURRENT_COMP);
    Table_Read(&rTable, reflect, 1); /* read 1st block data from file into rTable */
    Table_Rebin(&rTable);         /* rebin as evenly, increasing array */
    if (rTable.rows < 2) Table_Free(&rTable);
    Table_Info(rTable);
  } else rTable.data = NULL;

  if (width == 0) SlabWidth = zwidth;
  else SlabWidth = (width+gap)/NH - gap;
  if (height == 0) SlabHeight = yheight;
  else SlabHeight = (height+gap)/NV - gap;
  #undef reflect
  #undef zwidth
  #undef yheight
  #undef gap
  #undef NH
  #undef NV
  #undef mosaich
  #undef mosaicv
  #undef r0
  #undef Q
  #undef RV
  #undef RH
  #undef DM
  #undef mosaic
  #undef width
  #undef height
  #undef verbose
  #undef mos_y
  #undef mos_z
  #undef mono_Q
  #undef SlabWidth
  #undef SlabHeight
  #undef rTable
  return(_comp);
} /* class_Monochromator_2foc_init */

_class_Single_crystal *class_Single_crystal_init(_class_Single_crystal *_comp
) {
  #define reflections (_comp->_parameters.reflections)
  #define geometry (_comp->_parameters.geometry)
  #define mosaic_AB (_comp->_parameters.mosaic_AB)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define zdepth (_comp->_parameters.zdepth)
  #define radius (_comp->_parameters.radius)
  #define delta_d_d (_comp->_parameters.delta_d_d)
  #define mosaic (_comp->_parameters.mosaic)
  #define mosaic_a (_comp->_parameters.mosaic_a)
  #define mosaic_b (_comp->_parameters.mosaic_b)
  #define mosaic_c (_comp->_parameters.mosaic_c)
  #define recip_cell (_comp->_parameters.recip_cell)
  #define barns (_comp->_parameters.barns)
  #define ax (_comp->_parameters.ax)
  #define ay (_comp->_parameters.ay)
  #define az (_comp->_parameters.az)
  #define bx (_comp->_parameters.bx)
  #define by (_comp->_parameters.by)
  #define bz (_comp->_parameters.bz)
  #define cx (_comp->_parameters.cx)
  #define cy (_comp->_parameters.cy)
  #define cz (_comp->_parameters.cz)
  #define p_transmit (_comp->_parameters.p_transmit)
  #define sigma_abs (_comp->_parameters.sigma_abs)
  #define sigma_inc (_comp->_parameters.sigma_inc)
  #define aa (_comp->_parameters.aa)
  #define bb (_comp->_parameters.bb)
  #define cc (_comp->_parameters.cc)
  #define order (_comp->_parameters.order)
  #define RX (_comp->_parameters.RX)
  #define RY (_comp->_parameters.RY)
  #define RZ (_comp->_parameters.RZ)
  #define powder (_comp->_parameters.powder)
  #define PG (_comp->_parameters.PG)
  #define hkl_info (_comp->_parameters.hkl_info)
  #define offdata (_comp->_parameters.offdata)

  double as, bs, cs;
  int i=0;

  /* transfer input parameters */
  hkl_info.m_delta_d_d = delta_d_d;
  hkl_info.m_a  = 0;
  hkl_info.m_b  = 0;
  hkl_info.m_c  = 0;
  hkl_info.m_aa = aa;
  hkl_info.m_bb = bb;
  hkl_info.m_cc = cc;
  hkl_info.m_ax = ax;
  hkl_info.m_ay = ay;
  hkl_info.m_az = az;
  hkl_info.m_bx = bx;
  hkl_info.m_by = by;
  hkl_info.m_bz = bz;
  hkl_info.m_cx = cx;
  hkl_info.m_cy = cy;
  hkl_info.m_cz = cz;
  hkl_info.sigma_a = sigma_abs;
  hkl_info.sigma_i = sigma_inc;
  hkl_info.recip   = recip_cell;

  /* default format h,k,l,F,F2  */
  hkl_info.column_order[0]=1;
  hkl_info.column_order[1]=2;
  hkl_info.column_order[2]=3;
  hkl_info.column_order[3]=0;
  hkl_info.column_order[4]=7;
  hkl_info.kix = hkl_info.kiy = hkl_info.kiz = 0;
  hkl_info.nb_reuses = hkl_info.nb_refl = hkl_info.nb_refl_count = 0;
  hkl_info.tau_count = 0;

  /* ought to be cleaned up as mosaic_AB now is a proper vector/array and not a define */
  double* mosaic_ABin = mosaic_AB;
  /* Read in structure factors, and do some pre-calculations. */
  if (!read_hkl_data(reflections, &hkl_info, mosaic, mosaic_a, mosaic_b, mosaic_c, mosaic_ABin)) {
    printf("Single_crystal: %s: Error: Aborting.\n", NAME_CURRENT_COMP);
    exit(0);
  }

  if (hkl_info.sigma_a<0) hkl_info.sigma_a=0;
  if (hkl_info.sigma_i<0) hkl_info.sigma_i=0;

  if (hkl_info.count)
    printf("Single_crystal: %s: Read %d reflections from file '%s'\n",
      NAME_CURRENT_COMP, hkl_info.count, reflections);
  else printf("Single_crystal: %s: Using incoherent elastic scattering only sigma=%g.\n",
      NAME_CURRENT_COMP, hkl_info.sigma_i);

  hkl_info.shape=-1; /* -1:no shape, 0:cyl, 1:box, 2:sphere, 3:any-shape  */
  if (geometry && strlen(geometry) && strcmp(geometry, "NULL") && strcmp(geometry, "0")) {
	  if (off_init(geometry, xwidth, yheight, zdepth, 0, &offdata)) {
      hkl_info.shape=3;
    }
  }
  else if (xwidth && yheight && zdepth)  hkl_info.shape=1; /* box */
  else if (radius > 0 && yheight)        hkl_info.shape=0; /* cylinder */
  else if (radius > 0 && !yheight)       hkl_info.shape=2; /* sphere */

  if (hkl_info.shape < 0)
    exit(fprintf(stderr,"Single_crystal: %s: sample has invalid dimensions.\n"
                        "ERROR           Please check parameter values (xwidth, yheight, zdepth, radius).\n", NAME_CURRENT_COMP));

  printf("Single_crystal: %s: Vc=%g [Angs] sigma_abs=%g [barn] sigma_inc=%g [barn] reflections=%s\n",
      NAME_CURRENT_COMP, hkl_info.V0, hkl_info.sigma_a, hkl_info.sigma_i,
      reflections && strlen(reflections) ? reflections : "NULL");
      
  #ifdef USE_OPENCL
  oclContext_SX = oclInitKernel("hkl_search.cl", 1);
  
  if (oclContext_SX.Kernel) {
    int iDevice;
    int nDevice=1;
    for (iDevice = 0; iDevice < nDevice; iDevice++) {
      cl_int             ciErr1;                              // Error code var
      //  Initialization: load MT parameters and init host buffers    
      d_L = clCreateBuffer(oclContext_SX.GPUContext, CL_MEM_READ_WRITE, sizeof(struct hkl_data)*hkl_info.count, NULL, &ciErr1);
      if (ciErr1 !=CL_SUCCESS) { oclContext_SX.Kernel=NULL; break; }
    
      d_T = clCreateBuffer(oclContext_SX.GPUContext, CL_MEM_READ_WRITE, sizeof(struct tau_data)*hkl_info.count, NULL, &ciErr1);
      if (ciErr1 !=CL_SUCCESS) { oclContext_SX.Kernel=NULL; break; }
              
      d_tau_count = clCreateBuffer(oclContext_SX.GPUContext, CL_MEM_READ_WRITE, sizeof(int), NULL, &ciErr1);
      if (ciErr1 !=CL_SUCCESS) { oclContext_SX.Kernel=NULL; break; }
         
      d_coh_refl = clCreateBuffer(oclContext_SX.GPUContext, CL_MEM_READ_WRITE, sizeof(double), NULL, &ciErr1);
      if (ciErr1 !=CL_SUCCESS) { oclContext_SX.Kernel=NULL; break; }
          
      d_coh_xsect = clCreateBuffer(oclContext_SX.GPUContext, CL_MEM_READ_WRITE, sizeof(double), NULL, &ciErr1);
      if (ciErr1 !=CL_SUCCESS) { oclContext_SX.Kernel=NULL; break; }
      clFinish(oclContext_SX.CommandQueue[iDevice]);
      
    }
  }
  #endif

  if (PG && !(order==1)) {
    fprintf(stderr,"Single_crystal: %s: PG mode means implicit choice of no multiple scattering!\n"
	    "WARNING setting order=1\n", NAME_CURRENT_COMP);
    order=1;
  }


  #undef reflections
  #undef geometry
  #undef mosaic_AB
  #undef xwidth
  #undef yheight
  #undef zdepth
  #undef radius
  #undef delta_d_d
  #undef mosaic
  #undef mosaic_a
  #undef mosaic_b
  #undef mosaic_c
  #undef recip_cell
  #undef barns
  #undef ax
  #undef ay
  #undef az
  #undef bx
  #undef by
  #undef bz
  #undef cx
  #undef cy
  #undef cz
  #undef p_transmit
  #undef sigma_abs
  #undef sigma_inc
  #undef aa
  #undef bb
  #undef cc
  #undef order
  #undef RX
  #undef RY
  #undef RZ
  #undef powder
  #undef PG
  #undef hkl_info
  #undef offdata
  return(_comp);
} /* class_Single_crystal_init */

_class_PSD_monitor_4PI *class_PSD_monitor_4PI_init(_class_PSD_monitor_4PI *_comp
) {
  #define nx (_comp->_parameters.nx)
  #define ny (_comp->_parameters.ny)
  #define filename (_comp->_parameters.filename)
  #define radius (_comp->_parameters.radius)
  #define restore_neutron (_comp->_parameters.restore_neutron)
  #define PSD_N (_comp->_parameters.PSD_N)
  #define PSD_p (_comp->_parameters.PSD_p)
  #define PSD_p2 (_comp->_parameters.PSD_p2)
  int i, j;

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
  #undef nx
  #undef ny
  #undef filename
  #undef radius
  #undef restore_neutron
  #undef PSD_N
  #undef PSD_p
  #undef PSD_p2
  return(_comp);
} /* class_PSD_monitor_4PI_init */

_class_PSD_monitor *class_PSD_monitor_init(_class_PSD_monitor *_comp
) {
  #define nx (_comp->_parameters.nx)
  #define ny (_comp->_parameters.ny)
  #define filename (_comp->_parameters.filename)
  #define xmin (_comp->_parameters.xmin)
  #define xmax (_comp->_parameters.xmax)
  #define ymin (_comp->_parameters.ymin)
  #define ymax (_comp->_parameters.ymax)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define restore_neutron (_comp->_parameters.restore_neutron)
  #define PSD_N (_comp->_parameters.PSD_N)
  #define PSD_p (_comp->_parameters.PSD_p)
  #define PSD_p2 (_comp->_parameters.PSD_p2)
  if (xwidth  > 0) { xmax = xwidth/2;  xmin = -xmax; }
  if (yheight > 0) { ymax = yheight/2; ymin = -ymax; }

  if ((xmin >= xmax) || (ymin >= ymax)){
    printf("PSD_monitor: %s: Null detection area !\n"
           "ERROR        (xwidth,yheight,xmin,xmax,ymin,ymax). Exiting",
    NAME_CURRENT_COMP);
    exit(0);
  }

  PSD_N = create_darr2d(nx, ny);
  PSD_p = create_darr2d(nx, ny);
  PSD_p2 = create_darr2d(nx, ny);

  int i, j;
  for (i=0; i<nx; i++){
    for (j=0; j<ny; j++){
      PSD_N[i][j] = 0;
      PSD_p[i][j] = 0;
      PSD_p2[i][j] = 0;
    }
  }
  #undef nx
  #undef ny
  #undef filename
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef xwidth
  #undef yheight
  #undef restore_neutron
  #undef PSD_N
  #undef PSD_p
  #undef PSD_p2
  return(_comp);
} /* class_PSD_monitor_init */



int init(void) { /* called by mccode_main for Test_Monochromators:INITIALISE */
  DEBUG_INSTR();

  /* code_main/parseoptions/readparams sets instrument parameters value */
  stracpy(instrument->_name, "Test_Monochromators", 256);

  /* Instrument 'Test_Monochromators' INITIALISE */
  SIG_MESSAGE("[Test_Monochromators] INITIALISE [Test_Monochromators.instr:59]");
  #define Mono (instrument->_parameters._Mono)
  #define lambda (instrument->_parameters._lambda)
  #define Moz (instrument->_parameters._Moz)
  #define PG (instrument->_parameters._PG)
  #define powder (instrument->_parameters._powder)
  #define order (instrument->_parameters._order)
{
  int    ORDER = 1;
  double Ki;
  int    SM;

  /* SM : scattering at mono to the right (-1)/left(+1) */
  SM = 1;

  mono_q = 2*PI*ORDER/DM;  /* Q mono in Angs-1 */

  Ki = 2*PI/lambda;

  A2 = asin(mono_q/2/Ki)*RAD2DEG*2;

  A2 *= SM;   /* A1 : mono theta (crystal) */
  A1 = A2/2;  /* A2 : mono 2 theta (arm to sample) */

  printf("\n%s: ", NAME_INSTRUMENT);
  switch (Mono) {
  case 1:
    printf("Using Monochromator_flat\n"); break;
  case 2:
    printf("Using Monochromator_pol\n"); break;
  case 3:
    printf("Using Monochromator_pol (with 1%% events reflected)\n"); break;
  case 4:
    printf("Using Monochromator_curved (flat mode)\n"); break;
  case 5:
    printf("Using Monochromator_2foc (contrib, flat mode)\n"); break;
  case 6:
    printf("Using Single_crystal\n"); break;
  }
}
  #undef Mono
  #undef lambda
  #undef Moz
  #undef PG
  #undef powder
  #undef order
  _Origin_setpos(); /* type Progress_bar */
  _Source_setpos(); /* type Source_simple */
  _lamStart_setpos(); /* type L_monitor */
  _Mono_Arm_setpos(); /* type Arm */
  _Mono_Out_setpos(); /* type Arm */
  _Mono1_setpos(); /* type Monochromator_flat */
  _Mono2_setpos(); /* type Monochromator_pol */
  _Mono3_setpos(); /* type Monochromator_pol */
  _Mono4_setpos(); /* type Monochromator_curved */
  _Mono5_setpos(); /* type Monochromator_2foc */
  _Mono6_setpos(); /* type Single_crystal */
  _Sphere1_setpos(); /* type PSD_monitor_4PI */
  _lam1_setpos(); /* type L_monitor */
  _psd1_setpos(); /* type PSD_monitor */

  /* call iteratively all components INITIALISE */
  class_Progress_bar_init(_Origin);

  class_Source_simple_init(_Source);

  class_L_monitor_init(_lamStart);



  class_Monochromator_flat_init(_Mono1);

  class_Monochromator_pol_init(_Mono2);

  class_Monochromator_pol_init(_Mono3);

  class_Monochromator_curved_init(_Mono4);

  class_Monochromator_2foc_init(_Mono5);

  class_Single_crystal_init(_Mono6);

  class_PSD_monitor_4PI_init(_Sphere1);

  class_L_monitor_init(_lam1);

  class_PSD_monitor_init(_psd1);

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
_class_Source_simple *class_Source_simple_trace(_class_Source_simple *_comp
  , _class_particle *_particle) {
  ABSORBED=SCATTERED=RESTORE=0;

  #define radius (_comp->_parameters.radius)
  #define yheight (_comp->_parameters.yheight)
  #define xwidth (_comp->_parameters.xwidth)
  #define dist (_comp->_parameters.dist)
  #define focus_xw (_comp->_parameters.focus_xw)
  #define focus_yh (_comp->_parameters.focus_yh)
  #define E0 (_comp->_parameters.E0)
  #define dE (_comp->_parameters.dE)
  #define lambda0 (_comp->_parameters.lambda0)
  #define dlambda (_comp->_parameters.dlambda)
  #define flux (_comp->_parameters.flux)
  #define gauss (_comp->_parameters.gauss)
  #define target_index (_comp->_parameters.target_index)
  #define pmul (_comp->_parameters.pmul)
  #define square (_comp->_parameters.square)
  #define srcArea (_comp->_parameters.srcArea)
 double chi,E,lambda,v,r, xf, yf, rf, dx, dy, pdir;

 t=0;
 z=0;

 if (square == 1) {
   x = xwidth * (rand01() - 0.5);
   y = yheight * (rand01() - 0.5);
 } else {
   chi=2*PI*rand01();                          /* Choose point on source */
   r=sqrt(rand01())*radius;                    /* with uniform distribution. */
   x=r*cos(chi);
   y=r*sin(chi);
 }
 randvec_target_rect_real(&xf, &yf, &rf, &pdir,
			  0, 0, dist, focus_xw, focus_yh, ROT_A_CURRENT_COMP, x, y, z, 2);

 dx = xf-x;
 dy = yf-y;
 rf = sqrt(dx*dx+dy*dy+dist*dist);

 p = pdir*pmul;

 if(lambda0==0) {
   if (!gauss) {
     E=E0+dE*randpm1();              /*  Choose from uniform distribution */
   } else {
     E=E0+randnorm()*dE;
   }
   v=sqrt(E)*SE2V;
 } else {
   if (!gauss) {
     lambda=lambda0+dlambda*randpm1();
   } else {
     lambda=lambda0+randnorm()*dlambda;
   }
   v = K2V*(2*PI/lambda);
 }

 vz=v*dist/rf;
 vy=v*dy/rf;
 vx=v*dx/rf;
  #undef radius
  #undef yheight
  #undef xwidth
  #undef dist
  #undef focus_xw
  #undef focus_yh
  #undef E0
  #undef dE
  #undef lambda0
  #undef dlambda
  #undef flux
  #undef gauss
  #undef target_index
  #undef pmul
  #undef square
  #undef srcArea
  return(_comp);
} /* class_Source_simple_trace */

#pragma acc routine seq
_class_L_monitor *class_L_monitor_trace(_class_L_monitor *_comp
  , _class_particle *_particle) {
  ABSORBED=SCATTERED=RESTORE=0;

  #define nL (_comp->_parameters.nL)
  #define filename (_comp->_parameters.filename)
  #define xmin (_comp->_parameters.xmin)
  #define xmax (_comp->_parameters.xmax)
  #define ymin (_comp->_parameters.ymin)
  #define ymax (_comp->_parameters.ymax)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define Lmin (_comp->_parameters.Lmin)
  #define Lmax (_comp->_parameters.Lmax)
  #define restore_neutron (_comp->_parameters.restore_neutron)
  #define L_N (_comp->_parameters.L_N)
  #define L_p (_comp->_parameters.L_p)
  #define L_p2 (_comp->_parameters.L_p2)
  PROP_Z0;
  if (x>xmin && x<xmax && y>ymin && y<ymax)
  {
    double L = (2*PI/V2K)/sqrt(vx*vx + vy*vy + vz*vz);
    int i = floor((L-Lmin)*nL/(Lmax-Lmin));
    if(i >= 0 && i < nL)
    {
      L_N[i]++;
      L_p[i] += p;
      L_p2[i] += p*p;
      SCATTER;
    }
  }
  if (restore_neutron) {
    RESTORE_NEUTRON(INDEX_CURRENT_COMP, x, y, z, vx, vy, vz, t, sx, sy, sz, p);
  }
  #undef nL
  #undef filename
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef xwidth
  #undef yheight
  #undef Lmin
  #undef Lmax
  #undef restore_neutron
  #undef L_N
  #undef L_p
  #undef L_p2
  return(_comp);
} /* class_L_monitor_trace */

#pragma acc routine seq
_class_Monochromator_flat *class_Monochromator_flat_trace(_class_Monochromator_flat *_comp
  , _class_particle *_particle) {
  ABSORBED=SCATTERED=RESTORE=0;

  #define zmin (_comp->_parameters.zmin)
  #define zmax (_comp->_parameters.zmax)
  #define ymin (_comp->_parameters.ymin)
  #define ymax (_comp->_parameters.ymax)
  #define zwidth (_comp->_parameters.zwidth)
  #define yheight (_comp->_parameters.yheight)
  #define mosaich (_comp->_parameters.mosaich)
  #define mosaicv (_comp->_parameters.mosaicv)
  #define r0 (_comp->_parameters.r0)
  #define Q (_comp->_parameters.Q)
  #define DM (_comp->_parameters.DM)
  #define mos_rms_y (_comp->_parameters.mos_rms_y)
  #define mos_rms_z (_comp->_parameters.mos_rms_z)
  #define mos_rms_max (_comp->_parameters.mos_rms_max)
  #define mono_Q (_comp->_parameters.mono_Q)
  double y1,z1,t1,dt,kix,kiy,kiz,ratio,order,q0x,k,q0,theta;
  double bx,by,bz,kux,kuy,kuz,ax,ay,az,phi;
  double cos_2theta,k_sin_2theta,cos_phi,sin_phi,q_x,q_y,q_z;
  double delta,p_reflect,total,c1x,c1y,c1z,width,mos_sample;
  int i;

  if(vx != 0.0 && (dt = -x/vx) >= 0.0)
  {                             /* Moving towards crystal? */
    y1 = y + vy*dt;             /* Propagate to crystal plane */
    z1 = z + vz*dt;
    t1 = t + dt;
    if (z1>zmin && z1<zmax && y1>ymin && y1<ymax)
    {                           /* Intersect the crystal? */
      kix = V2K*vx;             /* Initial wave vector */
      kiy = V2K*vy;
      kiz = V2K*vz;
      /* Get reflection order and corresponding nominal scattering vector q0
         of correct length and direction. Only the order with the closest
         scattering vector is considered */
      ratio = -2*kix/mono_Q;
      order = floor(ratio + .5);
      if(order == 0.0)
        order = ratio < 0 ? -1 : 1;
      /* Order will be negative when the neutron enters from the back, in
         which case the direction of Q0 is flipped. */
      if(order < 0)
        order = -order;
      /* Make sure the order is small enough to allow Bragg scattering at the
         given neutron wavelength */
      k = sqrt(kix*kix + kiy*kiy + kiz*kiz);
      kux = kix/k;              /* Unit vector along ki */
      kuy = kiy/k;
      kuz = kiz/k;
      if(order > 2*k/mono_Q)
        order--;
      if(order > 0)             /* Bragg scattering possible? */
      {
        q0 = order*mono_Q;
        q0x = ratio < 0 ? -q0 : q0;
        theta = asin(q0/(2*k)); /* Actual bragg angle */
        /* Make MC choice: reflect or transmit? */
        delta = asin(fabs(kux)) - theta;
        p_reflect = r0*exp(-kiy*kiy/(kiy*kiy + kiz*kiz)*(delta*delta)/
                           (2*mos_rms_y*mos_rms_y))*
                       exp(-kiz*kiz/(kiy*kiy + kiz*kiz)*(delta*delta)/
                           (2*mos_rms_z*mos_rms_z));
        if(rand01() < p_reflect)
        {                       /* Reflect */
          cos_2theta = cos(2*theta);
          k_sin_2theta = k*sin(2*theta);
          /* Get unit normal to plane containing ki and most probable kf */
          vec_prod(bx, by, bz, kix, kiy, kiz, q0x, 0, 0);
          NORM(bx,by,bz);
          bx *= k_sin_2theta;
          by *= k_sin_2theta;
          bz *= k_sin_2theta;
          /* Get unit vector normal to ki and b */
          vec_prod(ax, ay, az, bx, by, bz, kux, kuy, kuz);
          /* Compute the total scattering probability at this ki */
          total = 0;
          /* Choose width of Gaussian distribution to sample the angle
           * phi on the Debye-Scherrer cone for the scattered neutron.
           * The radius of the Debye-Scherrer cone is smaller by a
           * factor 1/cos(theta) than the radius of the (partial) sphere
           * describing the possible orientations of Q due to mosaicity, so we
           * start with a width 1/cos(theta) greater than the largest of
           * the two mosaics. */
          mos_sample = mos_rms_max/cos(theta);
          c1x = kix*(cos_2theta-1);
          c1y = kiy*(cos_2theta-1);
          c1z = kiz*(cos_2theta-1);
          /* Loop, repeatedly reducing the sample width until it is small
           * enough to avoid sampling scattering directions with
           * ridiculously low scattering probability.
           * Use a cut-off at 5 times the gauss width for considering
           * scattering probability as well as for integration limits
           * when integrating the sampled distribution below. */
          for(i=0; i<100; i++) {
            width = 5*mos_sample;
            cos_phi = cos(width);
            sin_phi = sin(width);
            q_x = c1x + cos_phi*ax + sin_phi*bx;
            q_y = (c1y + cos_phi*ay + sin_phi*by)/mos_rms_y;
            q_z = (c1z + cos_phi*az + sin_phi*bz)/mos_rms_z;
            /* Stop when we get near a factor of 25=5^2. */
            if(q_z*q_z + q_y*q_y < (25/(2.0/3.0))*(q_x*q_x))
              break;
            mos_sample *= (2.0/3.0);
          }
          /* Now integrate the chosen sampling distribution, using a
           * cut-off at five times sigma. */
          for(i = 0; i < (sizeof(Gauss_X)/sizeof(double)); i++)
          {
            phi = width*Gauss_X[i];
            cos_phi = cos(phi);
            sin_phi = sin(phi);
            q_x = c1x + cos_phi*ax + sin_phi*bx;
            q_y = c1y + cos_phi*ay + sin_phi*by;
            q_z = c1z + cos_phi*az + sin_phi*bz;
            p_reflect = GAUSS((q_y/q_x),0,mos_rms_y)*
                        GAUSS((q_z/q_x),0,mos_rms_z);
            total += Gauss_W[i]*p_reflect;
          }
          total *= width;
          /* Choose point on Debye-Scherrer cone. Sample from a Gaussian of
           * width 1/cos(theta) greater than the mosaic and correct for any
           * error by adjusting the neutron weight later. */
          phi = mos_sample*randnorm();
          /* Compute final wave vector kf and scattering vector q = ki - kf */
          cos_phi = cos(phi);
          sin_phi = sin(phi);
          q_x = c1x + cos_phi*ax + sin_phi*bx;
          q_y = c1y + cos_phi*ay + sin_phi*by;
          q_z = c1z + cos_phi*az + sin_phi*bz;
          p_reflect = GAUSS((q_y/q_x),0,mos_rms_y)*
                      GAUSS((q_z/q_x),0,mos_rms_z);
          x = 0;
          y = y1;
          z = z1;
          t = t1;
          vx = K2V*(kix+q_x);
          vy = K2V*(kiy+q_y);
          vz = K2V*(kiz+q_z);
          p_reflect /= total*GAUSS(phi,0,mos_sample);
          if (p_reflect <= 0) ABSORB;
          if (p_reflect > 1)  p_reflect = 1;
          p *= p_reflect;
          SCATTER;
        } /* End MC choice to reflect or transmit neutron */
      } /* End bragg scattering possible */
    } /* End intersect the crystal */
  } /* End neutron moving towards crystal */
  #undef zmin
  #undef zmax
  #undef ymin
  #undef ymax
  #undef zwidth
  #undef yheight
  #undef mosaich
  #undef mosaicv
  #undef r0
  #undef Q
  #undef DM
  #undef mos_rms_y
  #undef mos_rms_z
  #undef mos_rms_max
  #undef mono_Q
  return(_comp);
} /* class_Monochromator_flat_trace */

#pragma acc routine seq
_class_Monochromator_pol *class_Monochromator_pol_trace(_class_Monochromator_pol *_comp
  , _class_particle *_particle) {
  ABSORBED=SCATTERED=RESTORE=0;

  #define zwidth (_comp->_parameters.zwidth)
  #define yheight (_comp->_parameters.yheight)
  #define mosaic (_comp->_parameters.mosaic)
  #define dspread (_comp->_parameters.dspread)
  #define Q (_comp->_parameters.Q)
  #define DM (_comp->_parameters.DM)
  #define pThreshold (_comp->_parameters.pThreshold)
  #define Rup (_comp->_parameters.Rup)
  #define Rdown (_comp->_parameters.Rdown)
  #define debug (_comp->_parameters.debug)
  #define mos_rms (_comp->_parameters.mos_rms)
  #define d_rms (_comp->_parameters.d_rms)
  #define mono_Q (_comp->_parameters.mono_Q)
  #define FN (_comp->_parameters.FN)
  #define FM (_comp->_parameters.FM)
  double y1, z1, t1, dt, vel;
  double sinTheta, lambdaBragg, lambda, dlambda2, sigmaLambda2, p_reflect;
  double R0; /* reflection probability based on FN and FM */
  double sx_in, sy_in, sz_in;
  int i;

  /* Propagate to crystal */
  PROP_X0;

  if (inside_rectangle(z, y, zwidth, yheight)) {/* Intersect the crystal? */

    // calculate sin(Bragg angle)
    vel = sqrt(vx*vx + vy*vy + vz*vz);
    sinTheta = abs(vx)/vel;

    // calculate lambdaBragg
    lambdaBragg = 2.0*DM*sinTheta;

    // calculate lambda of neutron
    lambda = 2*PI/(V2K*vel);


    // calculate deltalambda squared and sigmaLambda squared
    dlambda2 = (lambda-lambdaBragg)*(lambda-lambdaBragg);
    // The sigmaLambda is propagated by differentiating the Bragg
    // condition: lambda = 2*d*sinTheta
    sigmaLambda2 = 2.0*2.0 * sinTheta*sinTheta * d_rms*d_rms+
      2.0*2.0 * DM*DM * (1.0-sinTheta*sinTheta) * mos_rms*mos_rms;

    // calculate peak reflection probability
    GetMonoPolRefProb(FN, FM, sy, &R0);

    // calculate reflection probability
    p_reflect = R0*exp(-dlambda2/(2.0*sigmaLambda2));

    if(debug > 0) {
      printf("\n lambda: %f, Lambda_Bragg: %f\n", lambda, lambdaBragg);
      printf("sigmaLambda: %f, R0: %f, p_reflect: %f\n",
	     sqrt(sigmaLambda2), R0, p_reflect);
      printf("S_in:  (%f, %f, %f)\n", sx, sy, sz);
    }

    if((pThreshold>0 && p_reflect>pThreshold) || rand01()<p_reflect) {
      /* Reflect */

      // scale weight if neutron was accepted because of threshold
      if(pThreshold>0 && p_reflect>pThreshold)
	p*=p_reflect;

      vx = -vx;

      // Outgoing polarisation
      SetMonoPolRefOut(FN, FM, R0, &sx, &sy, &sz);

      if(debug > 0)
	printf("S_out: (%f, %f, %f)\n", sx, sy, sz);

      if(sx*sx+sy*sy+sz*sz>1)
        fprintf(stderr,"Pol_mirror: %s: Warning: polarisation |s| = %g > 1\n",
	      NAME_CURRENT_COMP, sx*sx+sy*sy+sz*sz); // check that polarisation is meaningfull

      SCATTER;
    } /* End MC choice to reflect or transmit neutron */
  } /* End intersect the crystal */

  #undef zwidth
  #undef yheight
  #undef mosaic
  #undef dspread
  #undef Q
  #undef DM
  #undef pThreshold
  #undef Rup
  #undef Rdown
  #undef debug
  #undef mos_rms
  #undef d_rms
  #undef mono_Q
  #undef FN
  #undef FM
  return(_comp);
} /* class_Monochromator_pol_trace */

#pragma acc routine seq
_class_Monochromator_curved *class_Monochromator_curved_trace(_class_Monochromator_curved *_comp
  , _class_particle *_particle) {
  ABSORBED=SCATTERED=RESTORE=0;

  #define reflect (_comp->_parameters.reflect)
  #define transmit (_comp->_parameters.transmit)
  #define zwidth (_comp->_parameters.zwidth)
  #define yheight (_comp->_parameters.yheight)
  #define gap (_comp->_parameters.gap)
  #define NH (_comp->_parameters.NH)
  #define NV (_comp->_parameters.NV)
  #define mosaich (_comp->_parameters.mosaich)
  #define mosaicv (_comp->_parameters.mosaicv)
  #define r0 (_comp->_parameters.r0)
  #define t0 (_comp->_parameters.t0)
  #define Q (_comp->_parameters.Q)
  #define RV (_comp->_parameters.RV)
  #define RH (_comp->_parameters.RH)
  #define DM (_comp->_parameters.DM)
  #define mosaic (_comp->_parameters.mosaic)
  #define width (_comp->_parameters.width)
  #define height (_comp->_parameters.height)
  #define verbose (_comp->_parameters.verbose)
  #define order (_comp->_parameters.order)
  #define mos_rms_y (_comp->_parameters.mos_rms_y)
  #define mos_rms_z (_comp->_parameters.mos_rms_z)
  #define mos_rms_max (_comp->_parameters.mos_rms_max)
  #define mono_Q (_comp->_parameters.mono_Q)
  #define SlabWidth (_comp->_parameters.SlabWidth)
  #define SlabHeight (_comp->_parameters.SlabHeight)
  #define rTable (_comp->_parameters.rTable)
  #define tTable (_comp->_parameters.tTable)
  #define row (_comp->_parameters.row)
  #define col (_comp->_parameters.col)
  #define tiltH (_comp->_parameters.tiltH)
  #define tiltV (_comp->_parameters.tiltV)
  double dt;

  if(vx != 0.0 && (dt = -x/vx) >= 0.0 && r0)
  {                             /* Moving towards crystal? */
    double zmin,zmax, ymin,ymax;

    /* Propagate to crystal plane */
    PROP_DT(dt);    /* now in the vertical plane of monochromator */

    zmax = ((NH*(SlabWidth+gap))-gap)/2;
    zmin = -zmax;
    ymax = ((NV*(SlabHeight+gap))-gap)/2;
    ymin = -ymax;

    /* hit a slab or a gap ? */

    if (z>zmin && z<zmax && y>ymin && y<ymax) { /* Intersect the crystal? */
      double tilth,tiltv;         /* used to calculate tilt angle of slab */
      double ratio, Q_order, k, kux,kuy,kuz;
      double kix,kiy,kiz;
      int    do_transmit = 0;

      col = ceil ( (z-zmin)/(SlabWidth +gap));  /* which slab hit ? */
      row = ceil ( (y-ymin)/(SlabHeight+gap));
      if (RH != 0) tilth = tiltH ? tiltH[(int)col] :  asin((col-(NH+1)/2)*(SlabWidth+gap)/RH);
      else tilth=0;
      if (RV != 0) tiltv = tiltV ? tiltV[(int)row] : -asin((row-(NV+1)/2)*(SlabHeight+gap)/RV);
      else tiltv=0;

      /* restore neutron in order to transform to slab coordinates */
      RESTORE_NEUTRON(INDEX_CURRENT_COMP, x, y, z, vx, vy, vz, t, sx, sy, sz, p);

      /* rotate with tilth (around Y) and tiltv (around Z), center on plate */
      double center_z=zmin+(col-0.5)*(SlabWidth+gap) -gap/2;
      double center_y=ymin+(row-0.5)*(SlabHeight+gap)-gap/2;
      Rotation T;
      rot_set_rotation(T, 0, tilth,    tiltv);
      /* now make the coordinate system change */
      mccoordschange_polarisation(T, &vx, &vy, &vz);
      y-=center_y;
      z-=center_z;
      coords_get(rot_apply(T,coords_set(x,y,z)),&x,&y,&z);

      /* this is where polaisation should be handled, plus further down */
      /* mccoordschange_polarisation(t, &sx, &sy, &sz); */

      /* now propagate to slab plane */
      PROP_X0;

      if (fabs(z) <= SlabWidth/2 && fabs(y) <= SlabHeight/2) { /* not in gap ? */
        kix = V2K*vx;             /* Initial wave vector */
        kiy = V2K*vy;
        kiz = V2K*vz;
        /* Get reflection order and corresponding nominal scattering vector q0
          of correct length and direction. Only the order with the closest
          scattering vector is considered */
        ratio = -2*kix/mono_Q;
        Q_order = floor(ratio + .5);
        if(Q_order == 0.0) Q_order = ratio < 0 ? -1 : 1;
        /* Order will be negative when the neutron enters from the back, in
          which case the direction of Q0 is flipped. */
        if(Q_order < 0) Q_order = -Q_order;
        /* Make sure the order is small enough to allow Bragg scattering at the
          given neutron wavelength */
        k = sqrt(kix*kix + kiy*kiy + kiz*kiz);
        kux = kix/k;              /* Unit vector along ki */
        kuy = kiy/k;
        kuz = kiz/k;
        if(Q_order > 2*k/mono_Q) Q_order--;
        if((!order && Q_order > 0) || (Q_order == fabs(order) && order)) {           /* Bragg scattering possible? */
          double q0, q0x, theta, delta, p_reflect, my_r0;

          q0 = Q_order*mono_Q;
          q0x = ratio < 0 ? -q0 : q0;
          theta = asin(q0/(2*k)); /* Actual bragg angle */
          /* Make MC choice: reflect or transmit? */
          delta = asin(fabs(kux)) - theta;

          if (rTable.data != NULL)
          {
            my_r0 = r0*Table_Value(rTable, k, 1); /* 2nd column */
          }
          else my_r0 = r0;
          if (my_r0 > 1)
          {
            if (my_r0 > 1.01 && verbose) fprintf(stdout, "Warning: Monochromator_curved : lowered reflectivity from %f to 1 (k=%f)\n", my_r0, k);
            my_r0=0.999;
          }
          if (my_r0 < 0)
          {
            if (verbose) fprintf(stdout, "Warning: Monochromator_curved : raised reflectivity from %f to 0 (k=%f)\n", my_r0, k);
            my_r0=0;
          }

          p_reflect = fabs(my_r0)*exp(-kiz*kiz/(kiy*kiy + kiz*kiz)*(delta*delta)/
                            (2*mos_rms_y*mos_rms_y))*
                        exp(-kiy*kiy/(kiy*kiy + kiz*kiz)*(delta*delta)/
                            (2*mos_rms_z*mos_rms_z));

          if(rand01() <= p_reflect) { /* Reflect */
            double bx,by,bz,ax,ay,az,phi;
            double cos_2theta,k_sin_2theta,cos_phi,sin_phi,q_x,q_y,q_z;
            double total,c1x,c1y,c1z,w,mos_sample;
            int i=0;

            cos_2theta   = cos(2*theta);
            k_sin_2theta = k*sin(2*theta);
            /* Get unit normal to plane containing ki and most probable kf */
            vec_prod(bx, by, bz, kix, kiy, kiz, q0x, 0, 0);
            NORM(bx,by,bz);
            bx *= k_sin_2theta;
            by *= k_sin_2theta;
            bz *= k_sin_2theta;
            /* Get unit vector normal to ki and b */
            vec_prod(ax, ay, az, bx, by, bz, kux, kuy, kuz);
            /* Compute the total scattering probability at this ki */
            total = 0;
            /* Choose width of Gaussian distribution to sample the angle
            * phi on the Debye-Scherrer cone for the scattered neutron.
            * The radius of the Debye-Scherrer cone is smaller by a
            * factor 1/cos(theta) than the radius of the (partial) sphere
            * describing the possible orientations of Q due to mosaicity, so we
            * start with a width 1/cos(theta) greater than the largest of
            * the two mosaics. */
            mos_sample = mos_rms_max/cos(theta);
            c1x = kix*(cos_2theta-1);
            c1y = kiy*(cos_2theta-1);
            c1z = kiz*(cos_2theta-1);
            /* Loop, repeatedly reducing the sample width until it is small
            * enough to avoid sampling scattering directions with
            * ridiculously low scattering probability.
            * Use a cut-off at 5 times the gauss width for considering
            * scattering probability as well as for integration limits
            * when integrating the sampled distribution below. */
            for(i=0; i<100; i++) {
              w = 5*mos_sample;
              cos_phi = cos(w);
              sin_phi = sin(w);
              q_x =  c1x + cos_phi*ax + sin_phi*bx;
              q_y = (c1y + cos_phi*ay + sin_phi*by)/mos_rms_z;
              q_z = (c1z + cos_phi*az + sin_phi*bz)/mos_rms_y;
              /* Stop when we get near a factor of 25=5^2. */
              if(q_z*q_z + q_y*q_y < (25/(2.0/3.0))*(q_x*q_x))
                break;
              mos_sample *= (2.0/3.0);
            }
            /* Now integrate the chosen sampling distribution, using a
            * cut-off at five times sigma. */
            for(i = 0; i < (sizeof(Gauss_X)/sizeof(double)); i++)
            {
              phi = w*Gauss_X[i];
              cos_phi = cos(phi);
              sin_phi = sin(phi);
              q_x = c1x + cos_phi*ax + sin_phi*bx;
              q_y = c1y + cos_phi*ay + sin_phi*by;
              q_z = c1z + cos_phi*az + sin_phi*bz;
              p_reflect = GAUSS((q_z/q_x),0,mos_rms_y)*
                          GAUSS((q_y/q_x),0,mos_rms_z);
              total += Gauss_W[i]*p_reflect;
            }
            total *= w;
            /* Choose point on Debye-Scherrer cone. Sample from a Gaussian of
             * width 1/cos(theta) greater than the mosaic and correct for any
             * error by adjusting the neutron weight later. */
            phi = mos_sample*randnorm();
            /* Compute final wave vector kf and scattering vector q = ki - kf */
            cos_phi = cos(phi);
            sin_phi = sin(phi);
            q_x = c1x + cos_phi*ax + sin_phi*bx;
            q_y = c1y + cos_phi*ay + sin_phi*by;
            q_z = c1z + cos_phi*az + sin_phi*bz;
            p_reflect = GAUSS((q_z/q_x),0,mos_rms_y)*
                        GAUSS((q_y/q_x),0,mos_rms_z);

            vx = K2V*(kix+q_x);
            vy = K2V*(kiy+q_y);
            vz = K2V*(kiz+q_z);
            p_reflect /= total*GAUSS(phi,0,mos_sample);
            if (p_reflect <= 0) ABSORB;
            if (p_reflect > 1)  p_reflect = 1;
            p *= p_reflect;

            SCATTER;
          } /* End MC choice to reflect or transmit neutron (if tmp<p_reflect) */
          else do_transmit = 1;
            /* else transmit neutron */
        } /* End bragg scattering possible (if order) */
        else do_transmit=1;
        if (do_transmit)
        {
          double my_t0;
          if (tTable.data != NULL)
          {
            my_t0 = t0*Table_Value(tTable, k, 1); /* 2nd column */
          }
          else my_t0 = t0;
          /* do not SCATTER, else GROUP does not work */
          if (my_t0 > 1)
          {
            if (my_t0 > 1.01 && verbose) fprintf(stdout, "Warning: Monochromator_curved : lowered transmission from %f to 1 (k=%f)\n", my_t0, k);
            my_t0=0.999;
          }
          if (my_t0 > 0) p*= my_t0;
          else ABSORB;
        }
      } /* end if not in gap */
      /* rotate back in component frame */
      Rotation TT;
      rot_transpose(T, TT);
      /* now make the coordinate system change */
      mccoordschange_polarisation(TT, &vx, &vy, &vz);
      coords_get(rot_apply(TT,coords_set(x,y,z)),&x,&y,&z);
      y+=center_y;
      z+=center_z;

      /* mccoordschange_polarisation(tt, &sx, &sy, &sz); */
    } /* End intersect the crystal (if z) */
    else {
      /* restore neutron state when no interaction */
      RESTORE_NEUTRON(INDEX_CURRENT_COMP, x, y, z, vx, vy, vz, t, sx, sy, sz, p);
    }
  } /* End neutron moving towards crystal (if vx)*/
  #undef reflect
  #undef transmit
  #undef zwidth
  #undef yheight
  #undef gap
  #undef NH
  #undef NV
  #undef mosaich
  #undef mosaicv
  #undef r0
  #undef t0
  #undef Q
  #undef RV
  #undef RH
  #undef DM
  #undef mosaic
  #undef width
  #undef height
  #undef verbose
  #undef order
  #undef mos_rms_y
  #undef mos_rms_z
  #undef mos_rms_max
  #undef mono_Q
  #undef SlabWidth
  #undef SlabHeight
  #undef rTable
  #undef tTable
  #undef row
  #undef col
  #undef tiltH
  #undef tiltV
  return(_comp);
} /* class_Monochromator_curved_trace */

#pragma acc routine seq
_class_Monochromator_2foc *class_Monochromator_2foc_trace(_class_Monochromator_2foc *_comp
  , _class_particle *_particle) {
  ABSORBED=SCATTERED=RESTORE=0;

  #define reflect (_comp->_parameters.reflect)
  #define zwidth (_comp->_parameters.zwidth)
  #define yheight (_comp->_parameters.yheight)
  #define gap (_comp->_parameters.gap)
  #define NH (_comp->_parameters.NH)
  #define NV (_comp->_parameters.NV)
  #define mosaich (_comp->_parameters.mosaich)
  #define mosaicv (_comp->_parameters.mosaicv)
  #define r0 (_comp->_parameters.r0)
  #define Q (_comp->_parameters.Q)
  #define RV (_comp->_parameters.RV)
  #define RH (_comp->_parameters.RH)
  #define DM (_comp->_parameters.DM)
  #define mosaic (_comp->_parameters.mosaic)
  #define width (_comp->_parameters.width)
  #define height (_comp->_parameters.height)
  #define verbose (_comp->_parameters.verbose)
  #define mos_y (_comp->_parameters.mos_y)
  #define mos_z (_comp->_parameters.mos_z)
  #define mono_Q (_comp->_parameters.mono_Q)
  #define SlabWidth (_comp->_parameters.SlabWidth)
  #define SlabHeight (_comp->_parameters.SlabHeight)
  #define rTable (_comp->_parameters.rTable)
  double dt;

  if(vx != 0.0 && (dt = -x/vx) >= 0.0)
  {
    double zmin,zmax, ymin,ymax, zp,yp, y1,z1,t1;

    zmax = ((NH*(SlabWidth+gap))-gap)/2;
    zmin = -1*zmax;
    ymax = ((NV*(SlabHeight+gap))-gap)/2;
    ymin = -1*ymax;
    y1 = y + vy*dt;             /* Propagate to crystal plane */
    z1 = z + vz*dt;
    t1 = t + dt;
    zp = fmod ( (z1-zmin),(SlabWidth+gap) );
    yp = fmod ( (y1-ymin),(SlabHeight+gap) );


    /* hit a slab or a gap ? */

    if (z1>zmin && z1<zmax && y1>ymin && y1<ymax && zp<SlabWidth && yp< SlabHeight)
    {
      double row,col, sna,snb,csa,csb,vxp,vyp,vzp;
      double v, theta0, theta, tmp3;
      double tilth,tiltv;         /* used to calculate tilt angle of slab */

      col = ceil ( (z1-zmin)/(SlabWidth+gap));  /* which slab hit ? */
      row = ceil ( (y1-ymin)/(SlabHeight+gap));
      if (RH != 0) tilth = asin((col-(NH+1)/2)*(SlabWidth+gap)/RH);
      else tilth=0;
      if (RV != 0) tiltv = -asin((row-(NV+1)/2)*(SlabHeight+gap)/RV);
      else tiltv=0;

      /* rotate with tilth and tiltv */

      sna = sin(tilth);
      snb = sin(tiltv);
      csa = cos(tilth);
      csb = cos(tiltv);
      vxp = vx*csa*csb+vy*snb-vz*sna*csb;
      vyp = -vx*csa*snb+vy*csb+vz*sna*snb;
      vzp = vx*sna+vz*csa;

      /* First: scattering in plane */
      /* theta0 = atan2(vx,vz);  neutron angle to slab Risoe version */

      v = sqrt(vxp*vxp+vyp*vyp+vzp*vzp);
      theta0 = asin(vxp/v);                /* correct neutron angle to slab */

      theta = asin(Q2V*mono_Q/(2.0*v));               /* Bragg's law */
      if (theta0 < 0)
              theta = -theta;
      tmp3 = (theta-theta0)/(MIN2RAD*mos_y);
      if (tmp3 < DIV_CUTOFF)
      {
        double my_r0, k;
        double dphi,tmp1,tmp2,tmp4,vratio,phi,cs,sn;

        k = V2K*v;

        if (rTable.data != NULL)
        {
          my_r0 = r0*Table_Value(rTable, k, 1); /* 2nd column */
        }
        else my_r0 = r0;
        if (my_r0 >= 1)
        {
          if (verbose) fprintf(stdout, "Warning: Monochromator_2foc: %s: lowered reflectivity from %f to 0.99 (k=%f)\n", 
            NAME_CURRENT_COMP, my_r0, k);
          my_r0=0.99;
        }
        if (my_r0 < 0)
        {
          if (verbose) fprintf(stdout, "Warning: Monochromator_2foc: %s: raised reflectivity from %f to 0 (k=%f)\n", 
          NAME_CURRENT_COMP, my_r0, k);
          my_r0=0;
        }
        x = 0.0;
        y = y1;
        z = z1;
        t = t1;
        
        /* reflectivity */
        t1 = fabs(my_r0)*exp(-tmp3*tmp3*4*log(2));
        if (t1 <= 0) ABSORB;
        if (t1 > 1)  t1 = 1;
        p *= t1; /* Use mosaics */
        
        tmp1 = 2*theta;
        cs = cos(tmp1);
        sn = sin(tmp1);
        tmp2 = cs*vxp - sn*vzp;
        vyp = vyp;
        /* vz = cs*vz + sn*vx; diese Zeile wurde durch die folgende ersetzt */
        tmp4 = vyp/vzp;  /* korrigiert den schr�en Einfall aufs Pl�tchen  */
        vzp = cs*(-vyp*sin(tmp4)+vzp*cos(tmp4)) + sn*vxp;
        vxp = tmp2;

        /* Second: scatering out of plane.
           Approximation is that Debye-Scherrer cone is a plane */

        phi = atan2(vyp,vzp);  /* out-of plane angle */
        dphi = (MIN2RAD*mos_z)/(2*sqrt(2*log(2)))*randnorm();  /* MC choice: */
        /* Vertical angle of the crystallite */
        vyp = vzp*tan(phi+2*dphi*sin(theta));
        vratio = v/sqrt(vxp*vxp+vyp*vyp+vzp*vzp);
        vzp = vzp*vratio;
        vyp = vyp*vratio;                             /* Renormalize v */
        vxp = vxp*vratio;

        /* rotate v coords back */
        vx = vxp*csb*csa-vyp*snb*csa+vzp*sna;
        vy = vxp*snb+vyp*csb;
        vz = -vxp*csb*sna+vyp*snb*sna+vzp*csa;
        /* v=sqrt(vx*vx+vy*vy+vz*vz);  */
        SCATTER;
      } /* end if Bragg ok */
    } /* End intersect the crystal (if z1) */
  } /* End neutron moving towards crystal (if vx)*/
  #undef reflect
  #undef zwidth
  #undef yheight
  #undef gap
  #undef NH
  #undef NV
  #undef mosaich
  #undef mosaicv
  #undef r0
  #undef Q
  #undef RV
  #undef RH
  #undef DM
  #undef mosaic
  #undef width
  #undef height
  #undef verbose
  #undef mos_y
  #undef mos_z
  #undef mono_Q
  #undef SlabWidth
  #undef SlabHeight
  #undef rTable
  return(_comp);
} /* class_Monochromator_2foc_trace */

#pragma acc routine seq
_class_Single_crystal *class_Single_crystal_trace(_class_Single_crystal *_comp
  , _class_particle *_particle) {
  ABSORBED=SCATTERED=RESTORE=0;

  #define reflections (_comp->_parameters.reflections)
  #define geometry (_comp->_parameters.geometry)
  #define mosaic_AB (_comp->_parameters.mosaic_AB)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define zdepth (_comp->_parameters.zdepth)
  #define radius (_comp->_parameters.radius)
  #define delta_d_d (_comp->_parameters.delta_d_d)
  #define mosaic (_comp->_parameters.mosaic)
  #define mosaic_a (_comp->_parameters.mosaic_a)
  #define mosaic_b (_comp->_parameters.mosaic_b)
  #define mosaic_c (_comp->_parameters.mosaic_c)
  #define recip_cell (_comp->_parameters.recip_cell)
  #define barns (_comp->_parameters.barns)
  #define ax (_comp->_parameters.ax)
  #define ay (_comp->_parameters.ay)
  #define az (_comp->_parameters.az)
  #define bx (_comp->_parameters.bx)
  #define by (_comp->_parameters.by)
  #define bz (_comp->_parameters.bz)
  #define cx (_comp->_parameters.cx)
  #define cy (_comp->_parameters.cy)
  #define cz (_comp->_parameters.cz)
  #define p_transmit (_comp->_parameters.p_transmit)
  #define sigma_abs (_comp->_parameters.sigma_abs)
  #define sigma_inc (_comp->_parameters.sigma_inc)
  #define aa (_comp->_parameters.aa)
  #define bb (_comp->_parameters.bb)
  #define cc (_comp->_parameters.cc)
  #define order (_comp->_parameters.order)
  #define RX (_comp->_parameters.RX)
  #define RY (_comp->_parameters.RY)
  #define RZ (_comp->_parameters.RZ)
  #define powder (_comp->_parameters.powder)
  #define PG (_comp->_parameters.PG)
  #define hkl_info (_comp->_parameters.hkl_info)
  #define offdata (_comp->_parameters.offdata)
  double t1, t2=0;                /* Entry and exit times in sample */
  struct hkl_data *L;           /* Structure factor list */
  int i;                        /* Index into structure factor list */
  struct tau_data *T;           /* List of reflections close to Ewald sphere */
  int j;                        /* Index into reflection list */
  int event_counter;            /* scattering event counter */
  double kix, kiy, kiz, ki;     /* Initial wave vector [1/AA] */
  double kfx, kfy, kfz;         /* Final wave vector */
  double v;                     /* Neutron velocity */
  double rho_x, rho_y, rho_z;   /* the vector ki - tau */
  double rho;
  double diff;                  /* Deviation from Bragg condition */
  double ox, oy, oz;            /* Origin of Ewald sphere tangent plane */
  double b1x, b1y, b1z;         /* First vector spanning tangent plane */
  double b2x, b2y, b2z;         /* Second vector spanning tangent plane */
  double n11, n12, n22;         /* 2D Gauss description matrix N */
  double det_N;                 /* Determinant of N */
  double inv_n11, inv_n12, inv_n22; /* Inverse of N */
  double l11, l12, l22;         /* Cholesky decomposition L of 1/2*inv(N) */
  double det_L;                 /* Determinant of L */
  double Bt_D_O_x, Bt_D_O_y;    /* Temporaries */
  double y0x, y0y;              /* Center of 2D Gauss in plane coordinates */
  double alpha;                 /* Offset of 2D Gauss center from 3D center */
  double V0;                    /* Volume of unit cell */
  double l_full;                /* Neutron path length for transmission */
  double l;                     /* Path length to scattering event */
  double abs_xsect, abs_xlen;   /* Absorption cross section and length */
  double inc_xsect, inc_xlen;   /* Incoherent scattering cross section and length */
  double coh_xlen;              /* Coherent cross section and length */
  double tot_xsect, tot_xlen;   /* Total cross section and length */
  double z1, z2, y1, y2;        /* Temporaries to choose kf from 2D Gauss */
  double adjust, sum;           /* Temporaries */

  double p_trans;               /* Transmission probability */
  double mc_trans, mc_interact; /* Transmission, interaction MC choices */
  int    intersect=0;
  double theta, phi;            /* rotation angles for curved lattice option */

  /* Intersection neutron trajectory / sample (sample surface) */
  if (hkl_info.shape == 0)
    intersect = cylinder_intersect(&t1, &t2, x, y, z, vx, vy, vz, radius, yheight);
  else if (hkl_info.shape == 1)
    intersect = box_intersect(&t1, &t2, x, y, z, vx, vy, vz, xwidth, yheight, zdepth);
  else if (hkl_info.shape == 2)
    intersect = sphere_intersect(&t1, &t2, x, y, z, vx, vy, vz, radius);
  else if (hkl_info.shape == 3)
    intersect = off_intersect(&t1, &t2, NULL, NULL, x, y, z, vx, vy, vz, offdata );

  if (t2 < 0) intersect=0;  /* we passed sample volume already */

  if(intersect)
  {                         /* Neutron intersects crystal */
    if(t1 > 0)
      PROP_DT(t1);          /* Move to crystal surface if not inside */
    v  = sqrt(vx*vx + vy*vy + vz*vz);
    ki = V2K*v;
    event_counter = 0;
    abs_xsect = hkl_info.sigma_a*2200/v;
    inc_xsect = hkl_info.sigma_i;
    V0= hkl_info.V0;
    abs_xlen  = abs_xsect/V0;
    inc_xlen  = inc_xsect/V0;

    /* Scalar cross sections for inc/abs are given in barns, so we need a scaling factor of 100
       to get scattering lengths in m, since V0 is assumed to be in AA*/
    abs_xlen *= 100; inc_xlen *= 100;

    L = hkl_info.list;
    T = hkl_info.tau_list;
    hkl_info.type = '\0';

    do {  /* Loop over multiple scattering events */
      /* Angles for powder randomization */
      double Alpha, Beta;

      if (hkl_info.shape == 0)
        intersect = cylinder_intersect(&t1, &t2, x, y, z, vx, vy, vz, radius, yheight);
      else if (hkl_info.shape == 1)
        intersect = box_intersect(&t1, &t2, x, y, z, vx, vy, vz, xwidth, yheight, zdepth);
      else if (hkl_info.shape == 2)
        intersect = sphere_intersect(&t1, &t2, x, y, z, vx, vy, vz, radius);
      else if (hkl_info.shape == 3)
        intersect = off_intersect(&t1, &t2, NULL, NULL, x, y, z, vx, vy, vz, offdata );
      if(!intersect || t2*v < -1e-9 || t1*v > 1e-9)
      {
        /* neutron is leaving the sample */
        if (hkl_info.flag_warning < 100)
          fprintf(stderr,
                "Single_crystal: %s: Warning: neutron has unexpectedly left the crystal!\n"
                "                t1=%g t2=%g x=%g y=%g z=%g vx=%g vy=%g vz=%g\n",
                NAME_CURRENT_COMP, t1, t2, x, y, z, vx, vy, vz);
        hkl_info.flag_warning++;
        break;
      }

      l_full = t2*v;

      /* (1). Compute incoming wave vector ki */
      if (powder) { /* orientation of crystallite is random */
	Alpha = randpm1()*PI*powder;
	Beta = randpm1()*PI/2;
	randrotate(&vx, &vy, &vz, Alpha, Beta);
      }
      if (PG) { /* orientation of crystallite is random */
	Alpha = rand01()*2*PI*PG;
	PGrotate(&vx, &vy, &vz, Alpha, hkl_info.csx, hkl_info.csy, hkl_info.csz);
      }


      /* lattice curvature option: rotate neutron velocity */
      if (RX || RY || RZ) {
        if (RX) { rotate(b1x,b1y,b1z,vx,vy,vz,    atan2(x,RX),0,0,1); }
		    else { b1x=vx; b1y=vy; b1z=vz; }
		    if (RY) { rotate(b2x,b2y,b2z,b1x,b1y,b1z, atan2(y,RY),1,0,0); }
		    else { b2x=b1x; b2y=b1y; b2z=b1z; }
		    if (RZ) { rotate(b1x,b1y,b1z,b2x,b2y,b2z, atan2(z,RZ),0,1,0); }
        else { b1x=b2x; b1y=b2y; b1z=b2z; }
		    kix = V2K*b1x;
        kiy = V2K*b1y;
        kiz = V2K*b1z;
      } else {
        kix = V2K*vx;
        kiy = V2K*vy;
        kiz = V2K*vz;
      }

      /* (2). Intersection of Ewald sphere with reciprocal lattice points */

      /* in case we use 'SPLIT' then consecutive neutrons can be identical when entering here
         and we may skip the hkl_search call */
      if ( fabs(kix - hkl_info.kix) < 1e-6
        && fabs(kiy - hkl_info.kiy) < 1e-6
        && fabs(kiz - hkl_info.kiz) < 1e-6) {
        hkl_info.nb_reuses++;
      } else {
        /* Max possible tau for this ki with 5*sigma delta-d/d cutoff. */
        double tau_max   = 2*ki/(1 - 5*hkl_info.m_delta_d_d);
        double coh_xsect = 0, coh_refl = 0;

        /* call hkl_search */
        #ifdef USE_OPENCL
        if (oclContext_SX.Kernel != NULL){ // the Kernel could be initialised
          hkl_info.tau_count = hkl_search_opencl(L, T, hkl_info.count, hkl_info.V0, 
              kix, kiy, kiz, tau_max,
              &coh_refl, &coh_xsect, oclContext_SX,
              d_L, d_T, d_tau_count, d_coh_refl, d_coh_xsect);
	  if (hkl_info.tau_count != 0)
	    printf("\nGPU hkl_info.tau_count:%i\n",hkl_info.tau_count);
        }
        else 
        #endif
	  hkl_info.tau_count = hkl_search(L, T, hkl_info.count, hkl_info.V0, 
              kix, kiy, kiz, tau_max,
              &coh_refl, &coh_xsect);
	
        /* store ki so that we can check for further SPLIT iterations */
        if (event_counter == 0 ) { /* only for incoming neutron */
          hkl_info.kix = kix;
          hkl_info.kiy = kiy;
          hkl_info.kiz = kiz;
        }

        hkl_info.coh_refl  = coh_refl;
        hkl_info.coh_xsect = coh_xsect;
        hkl_info.nb_refl += hkl_info.tau_count;
        hkl_info.nb_refl_count++;
      }

      /* (3). Probabilities of the different possible interactions. */
      /* Cross-sections are in barns = 10**-28 m**2, and unit cell volumes are
         in AA**3 = 10**-30 m**2. Hence a factor of 100 is used to convert
         scattering lengths to m**-1 */
      coh_xlen = hkl_info.coh_xsect/V0;
      if (barns) {
        coh_xlen *= 100;
      } /* else assume fm^2 */
      tot_xlen = abs_xlen + inc_xlen + coh_xlen;

      /* (5). Transmission */
      p_trans = exp(-tot_xlen*l_full);
      if(!event_counter && p_transmit >= 0 && p_transmit <= 1) {
        mc_trans = p_transmit; /* first event */
      } else {
        mc_trans = p_trans;
      }
      hkl_info.type = 't';
      mc_interact = 1 - mc_trans;
      if(mc_trans > 0 && (mc_trans >= 1 || rand01() < mc_trans))  /* Transmit */
      {
        p *= p_trans/mc_trans;
        intersect=0;
	if (powder) { /* orientation of crystallite is longer random */
	  randderotate(&vx, &vy, &vz, Alpha, Beta);
	}
	if (PG) { /* orientation of crystallite is longer random */
	  PGderotate(&vx, &vy, &vz, Alpha, hkl_info.m_cx, hkl_info.m_cz, hkl_info.m_cz);
	}
	break;
      }
      if(tot_xlen <= 0){
        ABSORB;
      }
      if(mc_interact <= 0)        /* Protect against rounding errors */
        { intersect=0;
	  if (powder) { /* orientation of crystallite is no longer random */
	    randderotate(&vx, &vy, &vz, Alpha, Beta);
	  }
	  if (PG) { /* orientation of crystallite is longer random */
	    PGderotate(&vx, &vy, &vz, Alpha, hkl_info.csx, hkl_info.csy, hkl_info.csz);
	  }
	  break;
	}
      if (!event_counter) p *= fabs(1 - p_trans)/mc_interact;
      /* Select a point at which to scatter the neutron, taking
         secondary extinction into account. */
      /* dP(l) = exp(-tot_xlen*l)dl
         P(l<l_0) = [-1/tot_xlen*exp(-tot_xlen*l)]_0^l_0
                  = (1 - exp(-tot_xlen*l0))/tot_xlen
         l = -log(1 - tot_xlen*rand0max(P(l<l_full)))/tot_xlen
       */
      if(tot_xlen*l_full < 1e-6)
        /* For very weak scattering, use simple uniform sampling of scattering
           point to avoid rounding errors. */
        l = rand0max(l_full);
      else
        l = -log(1 - rand0max((1 - exp(-tot_xlen*l_full))))/tot_xlen;
      PROP_DT(l/v);
      event_counter++;

      /* (4). Account for the probability of sigma_abs */
      p *= (coh_xlen + inc_xlen)/tot_xlen;
      /* Choose between coherent and incoherent scattering */
      if(coh_xlen == 0 || rand0max(coh_xlen + inc_xlen) <= inc_xlen)
      {
        /* (6). Incoherent scattering */
        randvec_target_circle(&kix, &kiy, &kiz, NULL, vx, vy, vz, 0);
        vx = kix; /* ki vector is used as tmp var with norm v */
        vy = kiy;
        vz = kiz; /* Go for next scattering event */
        hkl_info.type = 'i';
      } else {
        /* 7. Coherent scattering. Select reciprocal lattice point. */
        if(hkl_info.coh_refl <= 0){
          ABSORB;
        }
        sum = 0;
        j = hkl_select(T, hkl_info.tau_count, hkl_info.coh_refl, &sum);

        if(j >= hkl_info.tau_count)
        {
          if (hkl_info.flag_warning < 100)
            fprintf(stderr, "Single_crystal: Error: Illegal tau search "
              "(sum=%g, j=%i, tau_count=%i).\n", sum, j , hkl_info.tau_count);
          hkl_info.flag_warning++;
          j = hkl_info.tau_count - 1;
        }
        i = T[j].index;
        /* (8). Pick scattered wavevector kf from 2D Gauss distribution. */
        z1 = randnorm();
        z2 = randnorm();
        y1 = T[j].l11*z1 + T[j].y0x;
        y2 = T[j].l12*z1 + T[j].l22*z2 + T[j].y0y;
        kfx = T[j].rho_x + T[j].ox + T[j].b1x*y1 + T[j].b2x*y2;
        kfy = T[j].rho_y + T[j].oy + T[j].b1y*y1 + T[j].b2y*y2;
        kfz = T[j].rho_z + T[j].oz + T[j].b1z*y1 + T[j].b2z*y2;
        /* Normalize kf to length of ki, to account for planer
          approximation of the Ewald sphere. */
        adjust = ki/sqrt(kfx*kfx + kfy*kfy + kfz*kfz);
        kfx *= adjust;
        kfy *= adjust;
        kfz *= adjust;
        /* Adjust neutron weight (see manual for explanation). */
        p *= T[j].xsect*hkl_info.coh_refl/(hkl_info.coh_xsect*T[j].refl);

        vx = K2V*(L[i].u1x*kfx + L[i].u2x*kfy + L[i].u3x*kfz);
        vy = K2V*(L[i].u1y*kfx + L[i].u2y*kfy + L[i].u3y*kfz);
        vz = K2V*(L[i].u1z*kfx + L[i].u2z*kfy + L[i].u3z*kfz);
        hkl_info.type = 'c';
        hkl_info.h    = L[i].h;
        hkl_info.k    = L[i].k;
        hkl_info.l    = L[i].l;
      }
      /* lattice curvature option: rotate back neutron velocity */
      if (RZ || RY) {
        if (RX) { rotate(b1x,b1y,b1z,vx,vy,vz,    -atan2(x,RX),0,0,1); }
		    else { b1x=vx; b1y=vy; b1z=vz; }
		    if (RY) { rotate(b2x,b2y,b2z,b1x,b1y,b1z, -atan2(y,RY),1,0,0); }
		    else { b2x=b1x; b2y=b1y; b2z=b1z; }
		    if (RZ) { rotate(b1x,b1y,b1z,b2x,b2y,b2z, -atan2(z,RZ),0,1,0); }
        else { b1x=b2x; b1y=b2y; b1z=b2z; }
        vx=b1x; vy=b1y; vz=b1z;
	    }

      SCATTER;
      if (powder) { /* orientation of crystallite is no longer random */
	  randderotate(&vx, &vy, &vz, Alpha, Beta);
      }
      if (PG) { /* orientation of crystallite is longer random */
	PGderotate(&vx, &vy, &vz, Alpha, hkl_info.csx, hkl_info.csy, hkl_info.csz);
      }
      /* exit if multiple scattering order has been reached */
      if (order && event_counter >= order) { intersect=0; break; }
      /* Repeat loop for next scattering event. */
    } while (intersect); /* end do (intersect) (multiple scattering loop) */
  } /* if intersect */

  // EXTEND code here
  if (!strcmp(_comp->_name, "Mono6")) {
  if(!SCATTERED) ABSORB;
  }

  #undef reflections
  #undef geometry
  #undef mosaic_AB
  #undef xwidth
  #undef yheight
  #undef zdepth
  #undef radius
  #undef delta_d_d
  #undef mosaic
  #undef mosaic_a
  #undef mosaic_b
  #undef mosaic_c
  #undef recip_cell
  #undef barns
  #undef ax
  #undef ay
  #undef az
  #undef bx
  #undef by
  #undef bz
  #undef cx
  #undef cy
  #undef cz
  #undef p_transmit
  #undef sigma_abs
  #undef sigma_inc
  #undef aa
  #undef bb
  #undef cc
  #undef order
  #undef RX
  #undef RY
  #undef RZ
  #undef powder
  #undef PG
  #undef hkl_info
  #undef offdata
  return(_comp);
} /* class_Single_crystal_trace */

#pragma acc routine seq
_class_PSD_monitor_4PI *class_PSD_monitor_4PI_trace(_class_PSD_monitor_4PI *_comp
  , _class_particle *_particle) {
  ABSORBED=SCATTERED=RESTORE=0;

  #define nx (_comp->_parameters.nx)
  #define ny (_comp->_parameters.ny)
  #define filename (_comp->_parameters.filename)
  #define radius (_comp->_parameters.radius)
  #define restore_neutron (_comp->_parameters.restore_neutron)
  #define PSD_N (_comp->_parameters.PSD_N)
  #define PSD_p (_comp->_parameters.PSD_p)
  #define PSD_p2 (_comp->_parameters.PSD_p2)
  double t0, t1, phi, theta;
  int i,j;

  if(sphere_intersect(&t0, &t1, x, y, z, vx, vy, vz, radius) && t1 > 0)
  {
    if(t0 < 0)
      t0 = t1;
    /* t0 is now time of intersection with the sphere. */
    mcPROP_DT(t0);
    phi = atan2(x,z);
    i = floor(nx*(phi/(2*PI)+0.5));
    if(i >= nx)
      i = nx-1;                      /* Special case for phi = PI. */
    else if(i < 0)
      i = 0;
    theta=asin(y/radius);
    j = floor(ny*(theta+PI/2)/PI+0.5);
    if(j >= ny)
      j = ny-1;                      /* Special case for y = radius. */
    else if(j < 0)
      j = 0;
    PSD_N[i][j]++;
    PSD_p[i][j] += p;
    PSD_p2[i][j] += p*p;
    SCATTER;
  }
  if (restore_neutron) {
    RESTORE_NEUTRON(INDEX_CURRENT_COMP, x, y, z, vx, vy, vz, t, sx, sy, sz, p);
  }
  #undef nx
  #undef ny
  #undef filename
  #undef radius
  #undef restore_neutron
  #undef PSD_N
  #undef PSD_p
  #undef PSD_p2
  return(_comp);
} /* class_PSD_monitor_4PI_trace */

#pragma acc routine seq
_class_PSD_monitor *class_PSD_monitor_trace(_class_PSD_monitor *_comp
  , _class_particle *_particle) {
  ABSORBED=SCATTERED=RESTORE=0;

  #define nx (_comp->_parameters.nx)
  #define ny (_comp->_parameters.ny)
  #define filename (_comp->_parameters.filename)
  #define xmin (_comp->_parameters.xmin)
  #define xmax (_comp->_parameters.xmax)
  #define ymin (_comp->_parameters.ymin)
  #define ymax (_comp->_parameters.ymax)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define restore_neutron (_comp->_parameters.restore_neutron)
  #define PSD_N (_comp->_parameters.PSD_N)
  #define PSD_p (_comp->_parameters.PSD_p)
  #define PSD_p2 (_comp->_parameters.PSD_p2)
  PROP_Z0;
  if (x>xmin && x<xmax && y>ymin && y<ymax){
    int i = floor((x - xmin)*nx/(xmax - xmin));
    int j = floor((y - ymin)*ny/(ymax - ymin));
    PSD_N[i][j]++;
    PSD_p[i][j] += p;
    PSD_p2[i][j] += p*p;
    SCATTER;
  }
  if (restore_neutron) {
    RESTORE_NEUTRON(INDEX_CURRENT_COMP, x, y, z, vx, vy, vz, t, sx, sy, sz, p);
  }
  #undef nx
  #undef ny
  #undef filename
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef xwidth
  #undef yheight
  #undef restore_neutron
  #undef PSD_N
  #undef PSD_p
  #undef PSD_p2
  return(_comp);
} /* class_PSD_monitor_trace */

/* *****************************************************************************
* instrument 'Test_Monochromators' TRACE
***************************************************************************** */

#pragma acc routine seq
int raytrace(_class_particle* _particle) { /* called by mccode_main for Test_Monochromators:TRACE */

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
      /* component Origin=Progress_bar() [1] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Origin->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Origin->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Origin->_position_relative, _Origin->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Progress_bar_trace(_Origin, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component Origin [1] */
    if (!ABSORBED && _particle->_index == 2) {
      /* component Source=Source_simple() [2] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Source->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Source->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Source->_position_relative, _Source->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Source_simple_trace(_Source, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component Source [2] */
    if (!ABSORBED && _particle->_index == 3) {
      /* component lamStart=L_monitor() [3] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_lamStart->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _lamStart->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_lamStart->_position_relative, _lamStart->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_L_monitor_trace(_lamStart, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component lamStart [3] */
    if (!ABSORBED && _particle->_index == 4) {
      /* component Mono_Arm=Arm() [4] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Mono_Arm->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Mono_Arm->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Mono_Arm->_position_relative, _Mono_Arm->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      _particle->_index++;
    } /* end component Mono_Arm [4] */
    if (!ABSORBED && _particle->_index == 5) {
      /* component Mono_Out=Arm() [5] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Mono_Out->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Mono_Out->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Mono_Out->_position_relative, _Mono_Out->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      _particle->_index++;
    } /* end component Mono_Out [5] */
    if (!ABSORBED && _particle->_index == 6) {
      /* component Mono1=Monochromator_flat() [6] */
  #define zmin (_Mono1->_parameters.zmin)
  #define zmax (_Mono1->_parameters.zmax)
  #define ymin (_Mono1->_parameters.ymin)
  #define ymax (_Mono1->_parameters.ymax)
  #define zwidth (_Mono1->_parameters.zwidth)
  #define yheight (_Mono1->_parameters.yheight)
  #define mosaich (_Mono1->_parameters.mosaich)
  #define mosaicv (_Mono1->_parameters.mosaicv)
  #define r0 (_Mono1->_parameters.r0)
  #define Q (_Mono1->_parameters.Q)
  #define DM (_Mono1->_parameters.DM)
  #define mos_rms_y (_Mono1->_parameters.mos_rms_y)
  #define mos_rms_z (_Mono1->_parameters.mos_rms_z)
  #define mos_rms_max (_Mono1->_parameters.mos_rms_max)
  #define mono_Q (_Mono1->_parameters.mono_Q)
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Mono1->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Mono1->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Mono1->_position_relative, _Mono1->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      if (( instrument->_parameters._Mono == 1 )) // conditional WHEN execution
      class_Monochromator_flat_trace(_Mono1, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
  #undef zmin
  #undef zmax
  #undef ymin
  #undef ymax
  #undef zwidth
  #undef yheight
  #undef mosaich
  #undef mosaicv
  #undef r0
  #undef Q
  #undef DM
  #undef mos_rms_y
  #undef mos_rms_z
  #undef mos_rms_max
  #undef mono_Q
      _particle->_index++;
    } /* end component Mono1 [6] */
    if (!ABSORBED && _particle->_index == 7) {
      /* component Mono2=Monochromator_pol() [7] */
  #define zwidth (_Mono2->_parameters.zwidth)
  #define yheight (_Mono2->_parameters.yheight)
  #define mosaic (_Mono2->_parameters.mosaic)
  #define dspread (_Mono2->_parameters.dspread)
  #define Q (_Mono2->_parameters.Q)
  #define DM (_Mono2->_parameters.DM)
  #define pThreshold (_Mono2->_parameters.pThreshold)
  #define Rup (_Mono2->_parameters.Rup)
  #define Rdown (_Mono2->_parameters.Rdown)
  #define debug (_Mono2->_parameters.debug)
  #define mos_rms (_Mono2->_parameters.mos_rms)
  #define d_rms (_Mono2->_parameters.d_rms)
  #define mono_Q (_Mono2->_parameters.mono_Q)
  #define FN (_Mono2->_parameters.FN)
  #define FM (_Mono2->_parameters.FM)
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Mono2->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Mono2->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Mono2->_position_relative, _Mono2->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      if (( instrument->_parameters._Mono == 2 )) // conditional WHEN execution
      class_Monochromator_pol_trace(_Mono2, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
  #undef zwidth
  #undef yheight
  #undef mosaic
  #undef dspread
  #undef Q
  #undef DM
  #undef pThreshold
  #undef Rup
  #undef Rdown
  #undef debug
  #undef mos_rms
  #undef d_rms
  #undef mono_Q
  #undef FN
  #undef FM
      _particle->_index++;
    } /* end component Mono2 [7] */
    if (!ABSORBED && _particle->_index == 8) {
      /* component Mono3=Monochromator_pol() [8] */
  #define zwidth (_Mono3->_parameters.zwidth)
  #define yheight (_Mono3->_parameters.yheight)
  #define mosaic (_Mono3->_parameters.mosaic)
  #define dspread (_Mono3->_parameters.dspread)
  #define Q (_Mono3->_parameters.Q)
  #define DM (_Mono3->_parameters.DM)
  #define pThreshold (_Mono3->_parameters.pThreshold)
  #define Rup (_Mono3->_parameters.Rup)
  #define Rdown (_Mono3->_parameters.Rdown)
  #define debug (_Mono3->_parameters.debug)
  #define mos_rms (_Mono3->_parameters.mos_rms)
  #define d_rms (_Mono3->_parameters.d_rms)
  #define mono_Q (_Mono3->_parameters.mono_Q)
  #define FN (_Mono3->_parameters.FN)
  #define FM (_Mono3->_parameters.FM)
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Mono3->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Mono3->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Mono3->_position_relative, _Mono3->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      if (( instrument->_parameters._Mono == 3 )) // conditional WHEN execution
      class_Monochromator_pol_trace(_Mono3, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
  #undef zwidth
  #undef yheight
  #undef mosaic
  #undef dspread
  #undef Q
  #undef DM
  #undef pThreshold
  #undef Rup
  #undef Rdown
  #undef debug
  #undef mos_rms
  #undef d_rms
  #undef mono_Q
  #undef FN
  #undef FM
      _particle->_index++;
    } /* end component Mono3 [8] */
    if (!ABSORBED && _particle->_index == 9) {
      /* component Mono4=Monochromator_curved() [9] */
  #define reflect (_Mono4->_parameters.reflect)
  #define transmit (_Mono4->_parameters.transmit)
  #define zwidth (_Mono4->_parameters.zwidth)
  #define yheight (_Mono4->_parameters.yheight)
  #define gap (_Mono4->_parameters.gap)
  #define NH (_Mono4->_parameters.NH)
  #define NV (_Mono4->_parameters.NV)
  #define mosaich (_Mono4->_parameters.mosaich)
  #define mosaicv (_Mono4->_parameters.mosaicv)
  #define r0 (_Mono4->_parameters.r0)
  #define t0 (_Mono4->_parameters.t0)
  #define Q (_Mono4->_parameters.Q)
  #define RV (_Mono4->_parameters.RV)
  #define RH (_Mono4->_parameters.RH)
  #define DM (_Mono4->_parameters.DM)
  #define mosaic (_Mono4->_parameters.mosaic)
  #define width (_Mono4->_parameters.width)
  #define height (_Mono4->_parameters.height)
  #define verbose (_Mono4->_parameters.verbose)
  #define order (_Mono4->_parameters.order)
  #define mos_rms_y (_Mono4->_parameters.mos_rms_y)
  #define mos_rms_z (_Mono4->_parameters.mos_rms_z)
  #define mos_rms_max (_Mono4->_parameters.mos_rms_max)
  #define mono_Q (_Mono4->_parameters.mono_Q)
  #define SlabWidth (_Mono4->_parameters.SlabWidth)
  #define SlabHeight (_Mono4->_parameters.SlabHeight)
  #define rTable (_Mono4->_parameters.rTable)
  #define tTable (_Mono4->_parameters.tTable)
  #define row (_Mono4->_parameters.row)
  #define col (_Mono4->_parameters.col)
  #define tiltH (_Mono4->_parameters.tiltH)
  #define tiltV (_Mono4->_parameters.tiltV)
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Mono4->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Mono4->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Mono4->_position_relative, _Mono4->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      if (( instrument->_parameters._Mono == 4 )) // conditional WHEN execution
      class_Monochromator_curved_trace(_Mono4, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
  #undef reflect
  #undef transmit
  #undef zwidth
  #undef yheight
  #undef gap
  #undef NH
  #undef NV
  #undef mosaich
  #undef mosaicv
  #undef r0
  #undef t0
  #undef Q
  #undef RV
  #undef RH
  #undef DM
  #undef mosaic
  #undef width
  #undef height
  #undef verbose
  #undef order
  #undef mos_rms_y
  #undef mos_rms_z
  #undef mos_rms_max
  #undef mono_Q
  #undef SlabWidth
  #undef SlabHeight
  #undef rTable
  #undef tTable
  #undef row
  #undef col
  #undef tiltH
  #undef tiltV
      _particle->_index++;
    } /* end component Mono4 [9] */
    if (!ABSORBED && _particle->_index == 10) {
      /* component Mono5=Monochromator_2foc() [10] */
  #define reflect (_Mono5->_parameters.reflect)
  #define zwidth (_Mono5->_parameters.zwidth)
  #define yheight (_Mono5->_parameters.yheight)
  #define gap (_Mono5->_parameters.gap)
  #define NH (_Mono5->_parameters.NH)
  #define NV (_Mono5->_parameters.NV)
  #define mosaich (_Mono5->_parameters.mosaich)
  #define mosaicv (_Mono5->_parameters.mosaicv)
  #define r0 (_Mono5->_parameters.r0)
  #define Q (_Mono5->_parameters.Q)
  #define RV (_Mono5->_parameters.RV)
  #define RH (_Mono5->_parameters.RH)
  #define DM (_Mono5->_parameters.DM)
  #define mosaic (_Mono5->_parameters.mosaic)
  #define width (_Mono5->_parameters.width)
  #define height (_Mono5->_parameters.height)
  #define verbose (_Mono5->_parameters.verbose)
  #define mos_y (_Mono5->_parameters.mos_y)
  #define mos_z (_Mono5->_parameters.mos_z)
  #define mono_Q (_Mono5->_parameters.mono_Q)
  #define SlabWidth (_Mono5->_parameters.SlabWidth)
  #define SlabHeight (_Mono5->_parameters.SlabHeight)
  #define rTable (_Mono5->_parameters.rTable)
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Mono5->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Mono5->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Mono5->_position_relative, _Mono5->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      if (( instrument->_parameters._Mono == 5 )) // conditional WHEN execution
      class_Monochromator_2foc_trace(_Mono5, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
  #undef reflect
  #undef zwidth
  #undef yheight
  #undef gap
  #undef NH
  #undef NV
  #undef mosaich
  #undef mosaicv
  #undef r0
  #undef Q
  #undef RV
  #undef RH
  #undef DM
  #undef mosaic
  #undef width
  #undef height
  #undef verbose
  #undef mos_y
  #undef mos_z
  #undef mono_Q
  #undef SlabWidth
  #undef SlabHeight
  #undef rTable
      _particle->_index++;
    } /* end component Mono5 [10] */
    if (!ABSORBED && _particle->_index == 11) {
      /* component Mono6=Single_crystal() [11] */
  #define reflections (_Mono6->_parameters.reflections)
  #define geometry (_Mono6->_parameters.geometry)
  #define mosaic_AB (_Mono6->_parameters.mosaic_AB)
  #define xwidth (_Mono6->_parameters.xwidth)
  #define yheight (_Mono6->_parameters.yheight)
  #define zdepth (_Mono6->_parameters.zdepth)
  #define radius (_Mono6->_parameters.radius)
  #define delta_d_d (_Mono6->_parameters.delta_d_d)
  #define mosaic (_Mono6->_parameters.mosaic)
  #define mosaic_a (_Mono6->_parameters.mosaic_a)
  #define mosaic_b (_Mono6->_parameters.mosaic_b)
  #define mosaic_c (_Mono6->_parameters.mosaic_c)
  #define recip_cell (_Mono6->_parameters.recip_cell)
  #define barns (_Mono6->_parameters.barns)
  #define ax (_Mono6->_parameters.ax)
  #define ay (_Mono6->_parameters.ay)
  #define az (_Mono6->_parameters.az)
  #define bx (_Mono6->_parameters.bx)
  #define by (_Mono6->_parameters.by)
  #define bz (_Mono6->_parameters.bz)
  #define cx (_Mono6->_parameters.cx)
  #define cy (_Mono6->_parameters.cy)
  #define cz (_Mono6->_parameters.cz)
  #define p_transmit (_Mono6->_parameters.p_transmit)
  #define sigma_abs (_Mono6->_parameters.sigma_abs)
  #define sigma_inc (_Mono6->_parameters.sigma_inc)
  #define aa (_Mono6->_parameters.aa)
  #define bb (_Mono6->_parameters.bb)
  #define cc (_Mono6->_parameters.cc)
  #define order (_Mono6->_parameters.order)
  #define RX (_Mono6->_parameters.RX)
  #define RY (_Mono6->_parameters.RY)
  #define RZ (_Mono6->_parameters.RZ)
  #define powder (_Mono6->_parameters.powder)
  #define PG (_Mono6->_parameters.PG)
  #define hkl_info (_Mono6->_parameters.hkl_info)
  #define offdata (_Mono6->_parameters.offdata)
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Mono6->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Mono6->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Mono6->_position_relative, _Mono6->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      if (( instrument->_parameters._Mono == 6 )) // conditional WHEN execution
      class_Single_crystal_trace(_Mono6, _particle); /* contains EXTEND code */
      if (_particle->_restore)
        *_particle = _particle_save;
  #undef reflections
  #undef geometry
  #undef mosaic_AB
  #undef xwidth
  #undef yheight
  #undef zdepth
  #undef radius
  #undef delta_d_d
  #undef mosaic
  #undef mosaic_a
  #undef mosaic_b
  #undef mosaic_c
  #undef recip_cell
  #undef barns
  #undef ax
  #undef ay
  #undef az
  #undef bx
  #undef by
  #undef bz
  #undef cx
  #undef cy
  #undef cz
  #undef p_transmit
  #undef sigma_abs
  #undef sigma_inc
  #undef aa
  #undef bb
  #undef cc
  #undef order
  #undef RX
  #undef RY
  #undef RZ
  #undef powder
  #undef PG
  #undef hkl_info
  #undef offdata
      _particle->_index++;
    } /* end component Mono6 [11] */
    if (!ABSORBED && _particle->_index == 12) {
      /* component Sphere1=PSD_monitor_4PI() [12] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Sphere1->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Sphere1->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Sphere1->_position_relative, _Sphere1->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_PSD_monitor_4PI_trace(_Sphere1, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component Sphere1 [12] */
    if (!ABSORBED && _particle->_index == 13) {
      /* component lam1=L_monitor() [13] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_lam1->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _lam1->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_lam1->_position_relative, _lam1->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_L_monitor_trace(_lam1, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component lam1 [13] */
    if (!ABSORBED && _particle->_index == 14) {
      /* component psd1=PSD_monitor() [14] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_psd1->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _psd1->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_psd1->_position_relative, _psd1->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_PSD_monitor_trace(_psd1, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component psd1 [14] */
    if (_particle->_index > 14)
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
* instrument 'Test_Monochromators' and components SAVE
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

_class_L_monitor *class_L_monitor_save(_class_L_monitor *_comp
) {
  #define nL (_comp->_parameters.nL)
  #define filename (_comp->_parameters.filename)
  #define xmin (_comp->_parameters.xmin)
  #define xmax (_comp->_parameters.xmax)
  #define ymin (_comp->_parameters.ymin)
  #define ymax (_comp->_parameters.ymax)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define Lmin (_comp->_parameters.Lmin)
  #define Lmax (_comp->_parameters.Lmax)
  #define restore_neutron (_comp->_parameters.restore_neutron)
  #define L_N (_comp->_parameters.L_N)
  #define L_p (_comp->_parameters.L_p)
  #define L_p2 (_comp->_parameters.L_p2)
  DETECTOR_OUT_1D(
    "Wavelength monitor",
    "Wavelength [AA]",
    "Intensity",
    "L", Lmin, Lmax, nL,
    &L_N[0],&L_p[0],&L_p2[0],
    filename);
  #undef nL
  #undef filename
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef xwidth
  #undef yheight
  #undef Lmin
  #undef Lmax
  #undef restore_neutron
  #undef L_N
  #undef L_p
  #undef L_p2
  return(_comp);
} /* class_L_monitor_save */

_class_PSD_monitor_4PI *class_PSD_monitor_4PI_save(_class_PSD_monitor_4PI *_comp
) {
  #define nx (_comp->_parameters.nx)
  #define ny (_comp->_parameters.ny)
  #define filename (_comp->_parameters.filename)
  #define radius (_comp->_parameters.radius)
  #define restore_neutron (_comp->_parameters.restore_neutron)
  #define PSD_N (_comp->_parameters.PSD_N)
  #define PSD_p (_comp->_parameters.PSD_p)
  #define PSD_p2 (_comp->_parameters.PSD_p2)
  DETECTOR_OUT_2D(
    "4PI PSD monitor",
    "Longitude [deg]",
    "Lattitude [deg]",
    -180, 180, -90, 90,
    nx, ny,
    &PSD_N[0][0],&PSD_p[0][0],&PSD_p2[0][0],
    filename);
  #undef nx
  #undef ny
  #undef filename
  #undef radius
  #undef restore_neutron
  #undef PSD_N
  #undef PSD_p
  #undef PSD_p2
  return(_comp);
} /* class_PSD_monitor_4PI_save */

_class_PSD_monitor *class_PSD_monitor_save(_class_PSD_monitor *_comp
) {
  #define nx (_comp->_parameters.nx)
  #define ny (_comp->_parameters.ny)
  #define filename (_comp->_parameters.filename)
  #define xmin (_comp->_parameters.xmin)
  #define xmax (_comp->_parameters.xmax)
  #define ymin (_comp->_parameters.ymin)
  #define ymax (_comp->_parameters.ymax)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define restore_neutron (_comp->_parameters.restore_neutron)
  #define PSD_N (_comp->_parameters.PSD_N)
  #define PSD_p (_comp->_parameters.PSD_p)
  #define PSD_p2 (_comp->_parameters.PSD_p2)
  DETECTOR_OUT_2D(
    "PSD monitor",
    "X position [cm]",
    "Y position [cm]",
    xmin*100.0, xmax*100.0, ymin*100.0, ymax*100.0,
    nx, ny,
    &PSD_N[0][0],&PSD_p[0][0],&PSD_p2[0][0],
    filename);
  #undef nx
  #undef ny
  #undef filename
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef xwidth
  #undef yheight
  #undef restore_neutron
  #undef PSD_N
  #undef PSD_p
  #undef PSD_p2
  return(_comp);
} /* class_PSD_monitor_save */



int save(FILE *handle) { /* called by mccode_main for Test_Monochromators:SAVE */
  if (!handle) siminfo_init(NULL);

  /* call iteratively all components SAVE */
  class_Progress_bar_save(_Origin);


  class_L_monitor_save(_lamStart);









  class_PSD_monitor_4PI_save(_Sphere1);

  class_L_monitor_save(_lam1);

  class_PSD_monitor_save(_psd1);

  if (!handle) siminfo_close(); 

  return(0);
} /* save */

/* *****************************************************************************
* instrument 'Test_Monochromators' and components FINALLY
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

_class_L_monitor *class_L_monitor_finally(_class_L_monitor *_comp
) {
  #define nL (_comp->_parameters.nL)
  #define filename (_comp->_parameters.filename)
  #define xmin (_comp->_parameters.xmin)
  #define xmax (_comp->_parameters.xmax)
  #define ymin (_comp->_parameters.ymin)
  #define ymax (_comp->_parameters.ymax)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define Lmin (_comp->_parameters.Lmin)
  #define Lmax (_comp->_parameters.Lmax)
  #define restore_neutron (_comp->_parameters.restore_neutron)
  #define L_N (_comp->_parameters.L_N)
  #define L_p (_comp->_parameters.L_p)
  #define L_p2 (_comp->_parameters.L_p2)
  destroy_darr1d(L_N);
  destroy_darr1d(L_p);
  destroy_darr1d(L_p2);
  #undef nL
  #undef filename
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef xwidth
  #undef yheight
  #undef Lmin
  #undef Lmax
  #undef restore_neutron
  #undef L_N
  #undef L_p
  #undef L_p2
  return(_comp);
} /* class_L_monitor_finally */

_class_Monochromator_curved *class_Monochromator_curved_finally(_class_Monochromator_curved *_comp
) {
  #define reflect (_comp->_parameters.reflect)
  #define transmit (_comp->_parameters.transmit)
  #define zwidth (_comp->_parameters.zwidth)
  #define yheight (_comp->_parameters.yheight)
  #define gap (_comp->_parameters.gap)
  #define NH (_comp->_parameters.NH)
  #define NV (_comp->_parameters.NV)
  #define mosaich (_comp->_parameters.mosaich)
  #define mosaicv (_comp->_parameters.mosaicv)
  #define r0 (_comp->_parameters.r0)
  #define t0 (_comp->_parameters.t0)
  #define Q (_comp->_parameters.Q)
  #define RV (_comp->_parameters.RV)
  #define RH (_comp->_parameters.RH)
  #define DM (_comp->_parameters.DM)
  #define mosaic (_comp->_parameters.mosaic)
  #define width (_comp->_parameters.width)
  #define height (_comp->_parameters.height)
  #define verbose (_comp->_parameters.verbose)
  #define order (_comp->_parameters.order)
  #define mos_rms_y (_comp->_parameters.mos_rms_y)
  #define mos_rms_z (_comp->_parameters.mos_rms_z)
  #define mos_rms_max (_comp->_parameters.mos_rms_max)
  #define mono_Q (_comp->_parameters.mono_Q)
  #define SlabWidth (_comp->_parameters.SlabWidth)
  #define SlabHeight (_comp->_parameters.SlabHeight)
  #define rTable (_comp->_parameters.rTable)
  #define tTable (_comp->_parameters.tTable)
  #define row (_comp->_parameters.row)
  #define col (_comp->_parameters.col)
  #define tiltH (_comp->_parameters.tiltH)
  #define tiltV (_comp->_parameters.tiltV)
  Table_Free(&rTable);
  Table_Free(&tTable);
  if (tiltH) free(tiltH);
  if (tiltV) free(tiltV);
  #undef reflect
  #undef transmit
  #undef zwidth
  #undef yheight
  #undef gap
  #undef NH
  #undef NV
  #undef mosaich
  #undef mosaicv
  #undef r0
  #undef t0
  #undef Q
  #undef RV
  #undef RH
  #undef DM
  #undef mosaic
  #undef width
  #undef height
  #undef verbose
  #undef order
  #undef mos_rms_y
  #undef mos_rms_z
  #undef mos_rms_max
  #undef mono_Q
  #undef SlabWidth
  #undef SlabHeight
  #undef rTable
  #undef tTable
  #undef row
  #undef col
  #undef tiltH
  #undef tiltV
  return(_comp);
} /* class_Monochromator_curved_finally */

_class_Single_crystal *class_Single_crystal_finally(_class_Single_crystal *_comp
) {
  #define reflections (_comp->_parameters.reflections)
  #define geometry (_comp->_parameters.geometry)
  #define mosaic_AB (_comp->_parameters.mosaic_AB)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define zdepth (_comp->_parameters.zdepth)
  #define radius (_comp->_parameters.radius)
  #define delta_d_d (_comp->_parameters.delta_d_d)
  #define mosaic (_comp->_parameters.mosaic)
  #define mosaic_a (_comp->_parameters.mosaic_a)
  #define mosaic_b (_comp->_parameters.mosaic_b)
  #define mosaic_c (_comp->_parameters.mosaic_c)
  #define recip_cell (_comp->_parameters.recip_cell)
  #define barns (_comp->_parameters.barns)
  #define ax (_comp->_parameters.ax)
  #define ay (_comp->_parameters.ay)
  #define az (_comp->_parameters.az)
  #define bx (_comp->_parameters.bx)
  #define by (_comp->_parameters.by)
  #define bz (_comp->_parameters.bz)
  #define cx (_comp->_parameters.cx)
  #define cy (_comp->_parameters.cy)
  #define cz (_comp->_parameters.cz)
  #define p_transmit (_comp->_parameters.p_transmit)
  #define sigma_abs (_comp->_parameters.sigma_abs)
  #define sigma_inc (_comp->_parameters.sigma_inc)
  #define aa (_comp->_parameters.aa)
  #define bb (_comp->_parameters.bb)
  #define cc (_comp->_parameters.cc)
  #define order (_comp->_parameters.order)
  #define RX (_comp->_parameters.RX)
  #define RY (_comp->_parameters.RY)
  #define RZ (_comp->_parameters.RZ)
  #define powder (_comp->_parameters.powder)
  #define PG (_comp->_parameters.PG)
  #define hkl_info (_comp->_parameters.hkl_info)
  #define offdata (_comp->_parameters.offdata)
#ifdef USE_MPI
  if(mpi_node_rank == mpi_node_root) {
#endif

  if (hkl_info.flag_warning)
    fprintf(stderr, "Single_crystal: %s: Error message was repeated %i times with absorbed neutrons.\n",
      NAME_CURRENT_COMP, hkl_info.flag_warning);

  /* in case this instance is used in a SPLIT, we can recommend the
     optimal iteration value */
  if (hkl_info.nb_refl_count) {
    double split_iterations = (double)hkl_info.nb_reuses/hkl_info.nb_refl_count + 1;
    double split_optimal    = (double)hkl_info.nb_refl/hkl_info.nb_refl_count;
    if (split_optimal > split_iterations + 5)
      printf("Single_crystal: %s: Info: you may highly improve the computation efficiency by using\n"
        "    SPLIT %i COMPONENT %s=Single_crystal(order=1, ...)\n"
        "  in the instrument description %s.\n",
        NAME_CURRENT_COMP, (int)split_optimal, NAME_CURRENT_COMP, instrument_source);
  }
  #ifdef USE_OPENCL
  if (oclContext_SX.Kernel) {
    int iDevice=0;
    // clear OpenCL memory
   
    if (oclContext_SX.GPUContext)            clReleaseContext(oclContext_SX.GPUContext);
    if (oclContext_SX.Kernel)                clReleaseKernel(oclContext_SX.Kernel);    
    if (oclContext_SX.CommandQueue[iDevice]) clReleaseCommandQueue(oclContext_SX.CommandQueue[iDevice]);
    
    if (d_L)          clReleaseMemObject(d_L);
    if (d_T)          clReleaseMemObject(d_T);     
    if (d_tau_count)  clReleaseMemObject(d_tau_count);
    if (d_coh_refl)   clReleaseMemObject(d_coh_refl);
    if (d_coh_xsect)  clReleaseMemObject(d_coh_xsect);

  }
  #endif

#ifdef USE_MPI
  }
#endif
  #undef reflections
  #undef geometry
  #undef mosaic_AB
  #undef xwidth
  #undef yheight
  #undef zdepth
  #undef radius
  #undef delta_d_d
  #undef mosaic
  #undef mosaic_a
  #undef mosaic_b
  #undef mosaic_c
  #undef recip_cell
  #undef barns
  #undef ax
  #undef ay
  #undef az
  #undef bx
  #undef by
  #undef bz
  #undef cx
  #undef cy
  #undef cz
  #undef p_transmit
  #undef sigma_abs
  #undef sigma_inc
  #undef aa
  #undef bb
  #undef cc
  #undef order
  #undef RX
  #undef RY
  #undef RZ
  #undef powder
  #undef PG
  #undef hkl_info
  #undef offdata
  return(_comp);
} /* class_Single_crystal_finally */

_class_PSD_monitor_4PI *class_PSD_monitor_4PI_finally(_class_PSD_monitor_4PI *_comp
) {
  #define nx (_comp->_parameters.nx)
  #define ny (_comp->_parameters.ny)
  #define filename (_comp->_parameters.filename)
  #define radius (_comp->_parameters.radius)
  #define restore_neutron (_comp->_parameters.restore_neutron)
  #define PSD_N (_comp->_parameters.PSD_N)
  #define PSD_p (_comp->_parameters.PSD_p)
  #define PSD_p2 (_comp->_parameters.PSD_p2)
  destroy_darr2d(PSD_N);
  destroy_darr2d(PSD_p);
  destroy_darr2d(PSD_p2);
  #undef nx
  #undef ny
  #undef filename
  #undef radius
  #undef restore_neutron
  #undef PSD_N
  #undef PSD_p
  #undef PSD_p2
  return(_comp);
} /* class_PSD_monitor_4PI_finally */

_class_PSD_monitor *class_PSD_monitor_finally(_class_PSD_monitor *_comp
) {
  #define nx (_comp->_parameters.nx)
  #define ny (_comp->_parameters.ny)
  #define filename (_comp->_parameters.filename)
  #define xmin (_comp->_parameters.xmin)
  #define xmax (_comp->_parameters.xmax)
  #define ymin (_comp->_parameters.ymin)
  #define ymax (_comp->_parameters.ymax)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define restore_neutron (_comp->_parameters.restore_neutron)
  #define PSD_N (_comp->_parameters.PSD_N)
  #define PSD_p (_comp->_parameters.PSD_p)
  #define PSD_p2 (_comp->_parameters.PSD_p2)
  destroy_darr2d(PSD_N);
  destroy_darr2d(PSD_p);
  destroy_darr2d(PSD_p2);
  #undef nx
  #undef ny
  #undef filename
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef xwidth
  #undef yheight
  #undef restore_neutron
  #undef PSD_N
  #undef PSD_p
  #undef PSD_p2
  return(_comp);
} /* class_PSD_monitor_finally */



int finally(void) { /* called by mccode_main for Test_Monochromators:FINALLY */
  siminfo_init(NULL);
  save(siminfo_file); /* save data when simulation ends */

  /* call iteratively all components FINALLY */
  class_Progress_bar_finally(_Origin);


  class_L_monitor_finally(_lamStart);






  class_Monochromator_curved_finally(_Mono4);


  class_Single_crystal_finally(_Mono6);

  class_PSD_monitor_4PI_finally(_Sphere1);

  class_L_monitor_finally(_lam1);

  class_PSD_monitor_finally(_psd1);

  siminfo_close(); 

  return(0);
} /* finally */

/* *****************************************************************************
* instrument 'Test_Monochromators' and components DISPLAY
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

_class_Source_simple *class_Source_simple_display(_class_Source_simple *_comp
) {
  #define radius (_comp->_parameters.radius)
  #define yheight (_comp->_parameters.yheight)
  #define xwidth (_comp->_parameters.xwidth)
  #define dist (_comp->_parameters.dist)
  #define focus_xw (_comp->_parameters.focus_xw)
  #define focus_yh (_comp->_parameters.focus_yh)
  #define E0 (_comp->_parameters.E0)
  #define dE (_comp->_parameters.dE)
  #define lambda0 (_comp->_parameters.lambda0)
  #define dlambda (_comp->_parameters.dlambda)
  #define flux (_comp->_parameters.flux)
  #define gauss (_comp->_parameters.gauss)
  #define target_index (_comp->_parameters.target_index)
  #define pmul (_comp->_parameters.pmul)
  #define square (_comp->_parameters.square)
  #define srcArea (_comp->_parameters.srcArea)
  if (square == 1) {
    
    rectangle("xy",0,0,0,xwidth,yheight);
  } else {
    
    circle("xy",0,0,0,radius);
  }
  if (dist) {
    dashed_line(0,0,0, -focus_xw/2,-focus_yh/2,dist, 4);
    dashed_line(0,0,0,  focus_xw/2,-focus_yh/2,dist, 4);
    dashed_line(0,0,0,  focus_xw/2, focus_yh/2,dist, 4);
    dashed_line(0,0,0, -focus_xw/2, focus_yh/2,dist, 4);
  }
  #undef radius
  #undef yheight
  #undef xwidth
  #undef dist
  #undef focus_xw
  #undef focus_yh
  #undef E0
  #undef dE
  #undef lambda0
  #undef dlambda
  #undef flux
  #undef gauss
  #undef target_index
  #undef pmul
  #undef square
  #undef srcArea
  return(_comp);
} /* class_Source_simple_display */

_class_L_monitor *class_L_monitor_display(_class_L_monitor *_comp
) {
  #define nL (_comp->_parameters.nL)
  #define filename (_comp->_parameters.filename)
  #define xmin (_comp->_parameters.xmin)
  #define xmax (_comp->_parameters.xmax)
  #define ymin (_comp->_parameters.ymin)
  #define ymax (_comp->_parameters.ymax)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define Lmin (_comp->_parameters.Lmin)
  #define Lmax (_comp->_parameters.Lmax)
  #define restore_neutron (_comp->_parameters.restore_neutron)
  #define L_N (_comp->_parameters.L_N)
  #define L_p (_comp->_parameters.L_p)
  #define L_p2 (_comp->_parameters.L_p2)
  multiline(5, (double)xmin, (double)ymin, 0.0,
               (double)xmax, (double)ymin, 0.0,
               (double)xmax, (double)ymax, 0.0,
               (double)xmin, (double)ymax, 0.0,
               (double)xmin, (double)ymin, 0.0);
  #undef nL
  #undef filename
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef xwidth
  #undef yheight
  #undef Lmin
  #undef Lmax
  #undef restore_neutron
  #undef L_N
  #undef L_p
  #undef L_p2
  return(_comp);
} /* class_L_monitor_display */

_class_Arm *class_Arm_display(_class_Arm *_comp
) {
  /* A bit ugly; hard-coded dimensions. */
  
  line(0,0,0,0.2,0,0);
  line(0,0,0,0,0.2,0);
  line(0,0,0,0,0,0.2);
  return(_comp);
} /* class_Arm_display */

_class_Monochromator_flat *class_Monochromator_flat_display(_class_Monochromator_flat *_comp
) {
  #define zmin (_comp->_parameters.zmin)
  #define zmax (_comp->_parameters.zmax)
  #define ymin (_comp->_parameters.ymin)
  #define ymax (_comp->_parameters.ymax)
  #define zwidth (_comp->_parameters.zwidth)
  #define yheight (_comp->_parameters.yheight)
  #define mosaich (_comp->_parameters.mosaich)
  #define mosaicv (_comp->_parameters.mosaicv)
  #define r0 (_comp->_parameters.r0)
  #define Q (_comp->_parameters.Q)
  #define DM (_comp->_parameters.DM)
  #define mos_rms_y (_comp->_parameters.mos_rms_y)
  #define mos_rms_z (_comp->_parameters.mos_rms_z)
  #define mos_rms_max (_comp->_parameters.mos_rms_max)
  #define mono_Q (_comp->_parameters.mono_Q)
  
  multiline(5, 0.0, (double)ymin, (double)zmin,
               0.0, (double)ymax, (double)zmin,
               0.0, (double)ymax, (double)zmax,
               0.0, (double)ymin, (double)zmax,
               0.0, (double)ymin, (double)zmin);
  #undef zmin
  #undef zmax
  #undef ymin
  #undef ymax
  #undef zwidth
  #undef yheight
  #undef mosaich
  #undef mosaicv
  #undef r0
  #undef Q
  #undef DM
  #undef mos_rms_y
  #undef mos_rms_z
  #undef mos_rms_max
  #undef mono_Q
  return(_comp);
} /* class_Monochromator_flat_display */

_class_Monochromator_pol *class_Monochromator_pol_display(_class_Monochromator_pol *_comp
) {
  #define zwidth (_comp->_parameters.zwidth)
  #define yheight (_comp->_parameters.yheight)
  #define mosaic (_comp->_parameters.mosaic)
  #define dspread (_comp->_parameters.dspread)
  #define Q (_comp->_parameters.Q)
  #define DM (_comp->_parameters.DM)
  #define pThreshold (_comp->_parameters.pThreshold)
  #define Rup (_comp->_parameters.Rup)
  #define Rdown (_comp->_parameters.Rdown)
  #define debug (_comp->_parameters.debug)
  #define mos_rms (_comp->_parameters.mos_rms)
  #define d_rms (_comp->_parameters.d_rms)
  #define mono_Q (_comp->_parameters.mono_Q)
  #define FN (_comp->_parameters.FN)
  #define FM (_comp->_parameters.FM)
  
  rectangle("yz", 0, 0, 0, zwidth, yheight);
  #undef zwidth
  #undef yheight
  #undef mosaic
  #undef dspread
  #undef Q
  #undef DM
  #undef pThreshold
  #undef Rup
  #undef Rdown
  #undef debug
  #undef mos_rms
  #undef d_rms
  #undef mono_Q
  #undef FN
  #undef FM
  return(_comp);
} /* class_Monochromator_pol_display */

_class_Monochromator_curved *class_Monochromator_curved_display(_class_Monochromator_curved *_comp
) {
  #define reflect (_comp->_parameters.reflect)
  #define transmit (_comp->_parameters.transmit)
  #define zwidth (_comp->_parameters.zwidth)
  #define yheight (_comp->_parameters.yheight)
  #define gap (_comp->_parameters.gap)
  #define NH (_comp->_parameters.NH)
  #define NV (_comp->_parameters.NV)
  #define mosaich (_comp->_parameters.mosaich)
  #define mosaicv (_comp->_parameters.mosaicv)
  #define r0 (_comp->_parameters.r0)
  #define t0 (_comp->_parameters.t0)
  #define Q (_comp->_parameters.Q)
  #define RV (_comp->_parameters.RV)
  #define RH (_comp->_parameters.RH)
  #define DM (_comp->_parameters.DM)
  #define mosaic (_comp->_parameters.mosaic)
  #define width (_comp->_parameters.width)
  #define height (_comp->_parameters.height)
  #define verbose (_comp->_parameters.verbose)
  #define order (_comp->_parameters.order)
  #define mos_rms_y (_comp->_parameters.mos_rms_y)
  #define mos_rms_z (_comp->_parameters.mos_rms_z)
  #define mos_rms_max (_comp->_parameters.mos_rms_max)
  #define mono_Q (_comp->_parameters.mono_Q)
  #define SlabWidth (_comp->_parameters.SlabWidth)
  #define SlabHeight (_comp->_parameters.SlabHeight)
  #define rTable (_comp->_parameters.rTable)
  #define tTable (_comp->_parameters.tTable)
  #define row (_comp->_parameters.row)
  #define col (_comp->_parameters.col)
  #define tiltH (_comp->_parameters.tiltH)
  #define tiltV (_comp->_parameters.tiltV)
  int ih;

  
  for(ih = 0; ih < NH; ih++)
  {
    int iv;
    for(iv = 0; iv < NV; iv++)
    {
      double zmin,zmax,ymin,ymax;
      double xt, yt;

      zmin = (SlabWidth+gap)*(ih-NH/2.0)+gap/2;
      zmax = zmin+SlabWidth;
      ymin = (SlabHeight+gap)*(iv-NV/2.0)+gap/2;
      ymax = ymin+SlabHeight;

      if (RH) xt = -(zmax*zmax - zmin*zmin)/RH/2;
      else    xt = 0;

      if (RV) yt = -(ymax*ymax - ymin*ymin)/RV/2;
      else    yt = 0;
      multiline(5, xt+yt, (double)ymin, (double)zmin,
                   xt-yt, (double)ymax, (double)zmin,
                  -xt-yt, (double)ymax, (double)zmax,
                  -xt+yt, (double)ymin, (double)zmax,
                   xt+yt, (double)ymin, (double)zmin);
     }
   }
  #undef reflect
  #undef transmit
  #undef zwidth
  #undef yheight
  #undef gap
  #undef NH
  #undef NV
  #undef mosaich
  #undef mosaicv
  #undef r0
  #undef t0
  #undef Q
  #undef RV
  #undef RH
  #undef DM
  #undef mosaic
  #undef width
  #undef height
  #undef verbose
  #undef order
  #undef mos_rms_y
  #undef mos_rms_z
  #undef mos_rms_max
  #undef mono_Q
  #undef SlabWidth
  #undef SlabHeight
  #undef rTable
  #undef tTable
  #undef row
  #undef col
  #undef tiltH
  #undef tiltV
  return(_comp);
} /* class_Monochromator_curved_display */

_class_Monochromator_2foc *class_Monochromator_2foc_display(_class_Monochromator_2foc *_comp
) {
  #define reflect (_comp->_parameters.reflect)
  #define zwidth (_comp->_parameters.zwidth)
  #define yheight (_comp->_parameters.yheight)
  #define gap (_comp->_parameters.gap)
  #define NH (_comp->_parameters.NH)
  #define NV (_comp->_parameters.NV)
  #define mosaich (_comp->_parameters.mosaich)
  #define mosaicv (_comp->_parameters.mosaicv)
  #define r0 (_comp->_parameters.r0)
  #define Q (_comp->_parameters.Q)
  #define RV (_comp->_parameters.RV)
  #define RH (_comp->_parameters.RH)
  #define DM (_comp->_parameters.DM)
  #define mosaic (_comp->_parameters.mosaic)
  #define width (_comp->_parameters.width)
  #define height (_comp->_parameters.height)
  #define verbose (_comp->_parameters.verbose)
  #define mos_y (_comp->_parameters.mos_y)
  #define mos_z (_comp->_parameters.mos_z)
  #define mono_Q (_comp->_parameters.mono_Q)
  #define SlabWidth (_comp->_parameters.SlabWidth)
  #define SlabHeight (_comp->_parameters.SlabHeight)
  #define rTable (_comp->_parameters.rTable)
  int ih;

  
  for(ih = 0; ih < NH; ih++)
  {
    int iv;
    for(iv = 0; iv < NV; iv++)
    {
      double zmin,zmax,ymin,ymax;
      double xt, xt1, yt, yt1;

      zmin = (SlabWidth+gap)*(ih-NH/2.0)+gap/2;
      zmax = zmin+SlabWidth;
      ymin = (SlabHeight+gap)*(iv-NV/2.0)+gap/2;
      ymax = ymin+SlabHeight;

      if (RH)
      { xt = zmin*zmin/RH;
        xt1 = zmax*zmax/RH; }
      else { xt = 0; xt1 = 0; }

      if (RV)
      { yt = ymin*ymin/RV;
        yt1 = ymax*ymax/RV; }
      else { yt = 0; yt1 = 0; }
      multiline(5, xt+yt, (double)ymin, (double)zmin,
                   xt+yt1, (double)ymax, (double)zmin,
                   xt1+yt1, (double)ymax, (double)zmax,
                   xt1+yt, (double)ymin, (double)zmax,
                   xt+yt, (double)ymin, (double)zmin);
     }
   }
  #undef reflect
  #undef zwidth
  #undef yheight
  #undef gap
  #undef NH
  #undef NV
  #undef mosaich
  #undef mosaicv
  #undef r0
  #undef Q
  #undef RV
  #undef RH
  #undef DM
  #undef mosaic
  #undef width
  #undef height
  #undef verbose
  #undef mos_y
  #undef mos_z
  #undef mono_Q
  #undef SlabWidth
  #undef SlabHeight
  #undef rTable
  return(_comp);
} /* class_Monochromator_2foc_display */

_class_Single_crystal *class_Single_crystal_display(_class_Single_crystal *_comp
) {
  #define reflections (_comp->_parameters.reflections)
  #define geometry (_comp->_parameters.geometry)
  #define mosaic_AB (_comp->_parameters.mosaic_AB)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define zdepth (_comp->_parameters.zdepth)
  #define radius (_comp->_parameters.radius)
  #define delta_d_d (_comp->_parameters.delta_d_d)
  #define mosaic (_comp->_parameters.mosaic)
  #define mosaic_a (_comp->_parameters.mosaic_a)
  #define mosaic_b (_comp->_parameters.mosaic_b)
  #define mosaic_c (_comp->_parameters.mosaic_c)
  #define recip_cell (_comp->_parameters.recip_cell)
  #define barns (_comp->_parameters.barns)
  #define ax (_comp->_parameters.ax)
  #define ay (_comp->_parameters.ay)
  #define az (_comp->_parameters.az)
  #define bx (_comp->_parameters.bx)
  #define by (_comp->_parameters.by)
  #define bz (_comp->_parameters.bz)
  #define cx (_comp->_parameters.cx)
  #define cy (_comp->_parameters.cy)
  #define cz (_comp->_parameters.cz)
  #define p_transmit (_comp->_parameters.p_transmit)
  #define sigma_abs (_comp->_parameters.sigma_abs)
  #define sigma_inc (_comp->_parameters.sigma_inc)
  #define aa (_comp->_parameters.aa)
  #define bb (_comp->_parameters.bb)
  #define cc (_comp->_parameters.cc)
  #define order (_comp->_parameters.order)
  #define RX (_comp->_parameters.RX)
  #define RY (_comp->_parameters.RY)
  #define RZ (_comp->_parameters.RZ)
  #define powder (_comp->_parameters.powder)
  #define PG (_comp->_parameters.PG)
  #define hkl_info (_comp->_parameters.hkl_info)
  #define offdata (_comp->_parameters.offdata)
  
  if (hkl_info.shape == 0) {	/* cylinder */
    circle("xz", 0,  yheight/2.0, 0, radius);
    circle("xz", 0, -yheight/2.0, 0, radius);
    line(-radius, -yheight/2.0, 0, -radius, +yheight/2.0, 0);
    line(+radius, -yheight/2.0, 0, +radius, +yheight/2.0, 0);
    line(0, -yheight/2.0, -radius, 0, +yheight/2.0, -radius);
    line(0, -yheight/2.0, +radius, 0, +yheight/2.0, +radius);
  }
  else if (hkl_info.shape == 1) { 	/* box */
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
  else if (hkl_info.shape == 2) {	/* sphere */
    circle("xy", 0,  0.0, 0, radius);
    circle("xz", 0,  0.0, 0, radius);
    circle("yz", 0,  0.0, 0, radius);
  }
  else if (hkl_info.shape == 3) {	/* OFF file */
    off_display(offdata);
  }
  #undef reflections
  #undef geometry
  #undef mosaic_AB
  #undef xwidth
  #undef yheight
  #undef zdepth
  #undef radius
  #undef delta_d_d
  #undef mosaic
  #undef mosaic_a
  #undef mosaic_b
  #undef mosaic_c
  #undef recip_cell
  #undef barns
  #undef ax
  #undef ay
  #undef az
  #undef bx
  #undef by
  #undef bz
  #undef cx
  #undef cy
  #undef cz
  #undef p_transmit
  #undef sigma_abs
  #undef sigma_inc
  #undef aa
  #undef bb
  #undef cc
  #undef order
  #undef RX
  #undef RY
  #undef RZ
  #undef powder
  #undef PG
  #undef hkl_info
  #undef offdata
  return(_comp);
} /* class_Single_crystal_display */

_class_PSD_monitor_4PI *class_PSD_monitor_4PI_display(_class_PSD_monitor_4PI *_comp
) {
  #define nx (_comp->_parameters.nx)
  #define ny (_comp->_parameters.ny)
  #define filename (_comp->_parameters.filename)
  #define radius (_comp->_parameters.radius)
  #define restore_neutron (_comp->_parameters.restore_neutron)
  #define PSD_N (_comp->_parameters.PSD_N)
  #define PSD_p (_comp->_parameters.PSD_p)
  #define PSD_p2 (_comp->_parameters.PSD_p2)
  circle("xy",0,0,0,radius);
  circle("xz",0,0,0,radius);
  circle("yz",0,0,0,radius);
  #undef nx
  #undef ny
  #undef filename
  #undef radius
  #undef restore_neutron
  #undef PSD_N
  #undef PSD_p
  #undef PSD_p2
  return(_comp);
} /* class_PSD_monitor_4PI_display */

_class_PSD_monitor *class_PSD_monitor_display(_class_PSD_monitor *_comp
) {
  #define nx (_comp->_parameters.nx)
  #define ny (_comp->_parameters.ny)
  #define filename (_comp->_parameters.filename)
  #define xmin (_comp->_parameters.xmin)
  #define xmax (_comp->_parameters.xmax)
  #define ymin (_comp->_parameters.ymin)
  #define ymax (_comp->_parameters.ymax)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define restore_neutron (_comp->_parameters.restore_neutron)
  #define PSD_N (_comp->_parameters.PSD_N)
  #define PSD_p (_comp->_parameters.PSD_p)
  #define PSD_p2 (_comp->_parameters.PSD_p2)
  multiline(5,
    (double)xmin, (double)ymin, 0.0,
    (double)xmax, (double)ymin, 0.0,
    (double)xmax, (double)ymax, 0.0,
    (double)xmin, (double)ymax, 0.0,
    (double)xmin, (double)ymin, 0.0);
  #undef nx
  #undef ny
  #undef filename
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef xwidth
  #undef yheight
  #undef restore_neutron
  #undef PSD_N
  #undef PSD_p
  #undef PSD_p2
  return(_comp);
} /* class_PSD_monitor_display */


  #undef magnify
  #undef line
  #undef dashed_line
  #undef multiline
  #undef rectangle
  #undef box
  #undef circle
  #undef cylinder
  #undef sphere

int display(void) { /* called by mccode_main for Test_Monochromators:DISPLAY */
  printf("MCDISPLAY: start\n");

  /* call iteratively all components DISPLAY */
  class_Progress_bar_display(_Origin);

  class_Source_simple_display(_Source);

  class_L_monitor_display(_lamStart);

  class_Arm_display(_Mono_Arm);

  class_Arm_display(_Mono_Out);

  class_Monochromator_flat_display(_Mono1);

  class_Monochromator_pol_display(_Mono2);

  class_Monochromator_pol_display(_Mono3);

  class_Monochromator_curved_display(_Mono4);

  class_Monochromator_2foc_display(_Mono5);

  class_Single_crystal_display(_Mono6);

  class_PSD_monitor_4PI_display(_Sphere1);

  class_L_monitor_display(_lam1);

  class_PSD_monitor_display(_psd1);

  printf("MCDISPLAY: end\n");

  return(0);
} /* display */

void* _getvar_parameters(char* compname)
/* enables settings parameters based use of the GETPAR macro */
{
  if (!strcmp(compname, "Origin")) return (void *) &(_Origin->_parameters);
  if (!strcmp(compname, "Source")) return (void *) &(_Source->_parameters);
  if (!strcmp(compname, "lamStart")) return (void *) &(_lamStart->_parameters);
  if (!strcmp(compname, "Mono_Arm")) return (void *) &(_Mono_Arm->_parameters);
  if (!strcmp(compname, "Mono_Out")) return (void *) &(_Mono_Out->_parameters);
  if (!strcmp(compname, "Mono1")) return (void *) &(_Mono1->_parameters);
  if (!strcmp(compname, "Mono2")) return (void *) &(_Mono2->_parameters);
  if (!strcmp(compname, "Mono3")) return (void *) &(_Mono3->_parameters);
  if (!strcmp(compname, "Mono4")) return (void *) &(_Mono4->_parameters);
  if (!strcmp(compname, "Mono5")) return (void *) &(_Mono5->_parameters);
  if (!strcmp(compname, "Mono6")) return (void *) &(_Mono6->_parameters);
  if (!strcmp(compname, "Sphere1")) return (void *) &(_Sphere1->_parameters);
  if (!strcmp(compname, "lam1")) return (void *) &(_lam1->_parameters);
  if (!strcmp(compname, "psd1")) return (void *) &(_psd1->_parameters);
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

/* end of generated C code Test_Monochromators.c */
