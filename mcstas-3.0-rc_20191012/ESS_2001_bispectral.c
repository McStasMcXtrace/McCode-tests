/* Automatically generated file. Do not edit. 
 * Format:     ANSI C source code
 * Creator:    McStas <http://www.mcstas.org>
 * Instrument: ESS_2001_bispectral.instr (ESS_2001_bispectral)
 * Date:       Sat Oct 12 09:04:54 2019
 * File:       ESS_2001_bispectral.c
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
* Start of instrument 'ESS_2001_bispectral' generated code
***************************************************************************** */

#ifdef MC_TRACE_ENABLED
int traceenabled = 1;
#else
int traceenabled = 0;
#endif
#define MCSTAS "/usr/share/mcstas/3.0-dev/"
int   defaultmain         = 1;
char  instrument_name[]   = "ESS_2001_bispectral";
char  instrument_source[] = "ESS_2001_bispectral.instr";
char *instrument_exe      = NULL; /* will be set to argv[0] in main */
char  instrument_code[]   = "Instrument ESS_2001_bispectral source code ESS_2001_bispectral.instr is not embedded in this executable.\n  Use --source option when running McStas.\n";

int main(int argc, char *argv[]){return mccode_main(argc, argv);}

/* *****************************************************************************
* instrument 'ESS_2001_bispectral' and components DECLARE
***************************************************************************** */

/* Instrument parameters: structure and a table for the initialisation
   (Used in e.g. inputparse and I/O function (e.g. detector_out) */

struct _struct_instrument_parameters {
  long _man_lam;
  long _thermal;
  MCNUM _cold_Lam_min;
  MCNUM _cold_Lam_max;
  MCNUM _thermal_Lam_min;
  MCNUM _thermal_Lam_max;
  MCNUM _thermal_hotspot_fac;
  MCNUM _cold_hotspot_fac;
  long _use_full_guide_flag;
  MCNUM _guide_start;
  MCNUM _Length;
  MCNUM _focus_start_w;
  MCNUM _focus_end_w;
  MCNUM _smallaxis_w;
  MCNUM _focus_start_h;
  MCNUM _focus_end_h;
  MCNUM _smallaxis_h;
  MCNUM _maxdiv;
  MCNUM _coldperthermal;
  long _mirror_type;
  long _mirror_coating_type;
  MCNUM _mirror_offset;
  MCNUM _theta1;
  MCNUM _theta2;
  MCNUM _theta3;
  MCNUM _m_mirror;
  MCNUM _h_mirror;
  MCNUM _Pulse_width;
  MCNUM _frequency;
  MCNUM _gravity;
  MCNUM _substrate_thickness;
  MCNUM _coating_thickness;
  MCNUM _m1;
  MCNUM _m2;
};
typedef struct _struct_instrument_parameters _class_instrument_parameters;

/* instrument SPLIT, GROUP and JUMP control logic */
struct instrument_logic_struct {
};

struct _instrument_struct {
  char   _name[256]; /* the name of this instrument e.g. 'ESS_2001_bispectral' */
/* Counters per component instance */
  double counter_AbsorbProp[28]; /* absorbed events in PROP routines */
  double counter_N[28], counter_P[28], counter_P2[28]; /* event counters after each component instance */
  _class_particle _trajectory[28]; /* current trajectory for STORE/RESTORE */
/* Components position table (absolute and relative coords) */
  Coords _position_relative[28]; /* positions of all components */
  Coords _position_absolute[28];
  _class_instrument_parameters _parameters; /* instrument parameters */
  struct instrument_logic_struct logic; /* instrument logic */
} _instrument_var;
struct _instrument_struct *instrument = & _instrument_var;
#pragma acc declare create ( _instrument_var )
#pragma acc declare create ( instrument )

int numipar = 34;
struct mcinputtable_struct mcinputtable[] = {
  "man_lam", &(_instrument_var._parameters._man_lam), instr_type_int, "1", 
  "thermal", &(_instrument_var._parameters._thermal), instr_type_int, "2", 
  "cold_Lam_min", &(_instrument_var._parameters._cold_Lam_min), instr_type_double, "0.1", 
  "cold_Lam_max", &(_instrument_var._parameters._cold_Lam_max), instr_type_double, "10", 
  "thermal_Lam_min", &(_instrument_var._parameters._thermal_Lam_min), instr_type_double, "0.1", 
  "thermal_Lam_max", &(_instrument_var._parameters._thermal_Lam_max), instr_type_double, "10", 
  "thermal_hotspot_fac", &(_instrument_var._parameters._thermal_hotspot_fac), instr_type_double, "1.0", 
  "cold_hotspot_fac", &(_instrument_var._parameters._cold_hotspot_fac), instr_type_double, "1.0", 
  "use_full_guide_flag", &(_instrument_var._parameters._use_full_guide_flag), instr_type_int, "0", 
  "guide_start", &(_instrument_var._parameters._guide_start), instr_type_double, "2", 
  "Length", &(_instrument_var._parameters._Length), instr_type_double, "150", 
  "focus_start_w", &(_instrument_var._parameters._focus_start_w), instr_type_double, "-1.2620", 
  "focus_end_w", &(_instrument_var._parameters._focus_end_w), instr_type_double, "150.3187", 
  "smallaxis_w", &(_instrument_var._parameters._smallaxis_w), instr_type_double, "0.2550", 
  "focus_start_h", &(_instrument_var._parameters._focus_start_h), instr_type_double, "-2.1415", 
  "focus_end_h", &(_instrument_var._parameters._focus_end_h), instr_type_double, "150.0778", 
  "smallaxis_h", &(_instrument_var._parameters._smallaxis_h), instr_type_double, "0.3847", 
  "maxdiv", &(_instrument_var._parameters._maxdiv), instr_type_double, "2.0", 
  "coldperthermal", &(_instrument_var._parameters._coldperthermal), instr_type_double, "30", 
  "mirror_type", &(_instrument_var._parameters._mirror_type), instr_type_int, "2", 
  "mirror_coating_type", &(_instrument_var._parameters._mirror_coating_type), instr_type_int, "1", 
  "mirror_offset", &(_instrument_var._parameters._mirror_offset), instr_type_double, "0", 
  "theta1", &(_instrument_var._parameters._theta1), instr_type_double, "1.25", 
  "theta2", &(_instrument_var._parameters._theta2), instr_type_double, "1.25", 
  "theta3", &(_instrument_var._parameters._theta3), instr_type_double, "1.25", 
  "m_mirror", &(_instrument_var._parameters._m_mirror), instr_type_double, "5", 
  "h_mirror", &(_instrument_var._parameters._h_mirror), instr_type_double, "0.15", 
  "Pulse_width", &(_instrument_var._parameters._Pulse_width), instr_type_double, "0.00286", 
  "frequency", &(_instrument_var._parameters._frequency), instr_type_double, "14.0", 
  "gravity", &(_instrument_var._parameters._gravity), instr_type_double, "-9.81", 
  "substrate_thickness", &(_instrument_var._parameters._substrate_thickness), instr_type_double, "0.0005", 
  "coating_thickness", &(_instrument_var._parameters._coating_thickness), instr_type_double, "10e-6", 
  "m1", &(_instrument_var._parameters._m1), instr_type_double, "5", 
  "m2", &(_instrument_var._parameters._m2), instr_type_double, "4", 
  NULL, NULL, instr_type_double, ""
};


/* ************************************************************************** */
/*             SHARE user declarations for all components                     */
/* ************************************************************************** */

/* Shared user declarations for all components types 'ESS_moderator_long_2001'. */
double Mezei_M_fct(double l, double temp)
  {
    double a=949.0/temp;
    return 2*a*a*exp(-a/(l*l))/(l*l*l*l*l);
  }

  double Mezei_F_fct(double t, double tau, int n)
  {
    return (exp(-t/tau)-exp(-n*t/tau))*n/(n-1)/tau;
  }

/* Shared user declarations for all components types 'Mirror_Curved_Bispectral'. */
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


/* Shared user declarations for all components types 'Mirror_Elliptic_Bispectral'. */



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

/* component ArmForGuideRight=Arm() [2] DECLARE */
/* Parameter definition for component type 'Arm' */
struct _struct_Arm_parameters {
  char Arm_has_no_parameters;
}; /* _struct_Arm_parameters */
typedef struct _struct_Arm_parameters _class_Arm_parameters;

/* Parameters for component type 'Arm' */
struct _struct_Arm {
  char     _name[256]; /* e.g. ArmForGuideRight */
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
_class_Arm _ArmForGuideRight_var;
_class_Arm *_ArmForGuideRight = &_ArmForGuideRight_var;
#pragma acc declare create ( _ArmForGuideRight_var )
#pragma acc declare create ( _ArmForGuideRight )

_class_Arm _ArmForGuideBottom_var;
_class_Arm *_ArmForGuideBottom = &_ArmForGuideBottom_var;
#pragma acc declare create ( _ArmForGuideBottom_var )
#pragma acc declare create ( _ArmForGuideBottom )

_class_Arm _ArmForGuideTop_var;
_class_Arm *_ArmForGuideTop = &_ArmForGuideTop_var;
#pragma acc declare create ( _ArmForGuideTop_var )
#pragma acc declare create ( _ArmForGuideTop )

_class_Arm _ArmForGuideLeft_var;
_class_Arm *_ArmForGuideLeft = &_ArmForGuideLeft_var;
#pragma acc declare create ( _ArmForGuideLeft_var )
#pragma acc declare create ( _ArmForGuideLeft )

/* component cold_source=ESS_moderator_long_2001() [6] DECLARE */
/* Parameter definition for component type 'ESS_moderator_long_2001' */
struct _struct_ESS_moderator_long_2001_parameters {
  /* Component type 'ESS_moderator_long_2001' setting parameters */
  MCNUM size;
  MCNUM l_low;
  MCNUM l_high;
  MCNUM dist;
  MCNUM xw;
  MCNUM yh;
  MCNUM freq;
  MCNUM T;
  MCNUM tau;
  MCNUM tau1;
  MCNUM tau2;
  MCNUM d;
  MCNUM n;
  MCNUM n2;
  MCNUM chi2;
  MCNUM I0;
  MCNUM I2;
  MCNUM branch1;
  MCNUM branch2;
  MCNUM branch_tail;
  MCNUM twopulses;
  long target_index;
  /* Component type 'ESS_moderator_long_2001' private parameters */
  /* Component type 'ESS_moderator_long_2001' DECLARE code stored as structure members */
  double l_range, w_mult, branchframe, tx, ty, tz;
}; /* _struct_ESS_moderator_long_2001_parameters */
typedef struct _struct_ESS_moderator_long_2001_parameters _class_ESS_moderator_long_2001_parameters;

/* Parameters for component type 'ESS_moderator_long_2001' */
struct _struct_ESS_moderator_long_2001 {
  char     _name[256]; /* e.g. cold_source */
  char     _type[256]; /* ESS_moderator_long_2001 */
  long     _index; /* e.g. 2 index in TRACE list */
  Coords   _position_absolute;
  Coords   _position_relative; /* wrt PREVIOUS */
  Rotation _rotation_absolute;
  Rotation _rotation_relative; /* wrt PREVIOUS */
  int      _rotation_is_identity;
  _class_ESS_moderator_long_2001_parameters _parameters;
};
typedef struct _struct_ESS_moderator_long_2001 _class_ESS_moderator_long_2001;
_class_ESS_moderator_long_2001 _cold_source_var;
_class_ESS_moderator_long_2001 *_cold_source = &_cold_source_var;
#pragma acc declare create ( _cold_source_var )
#pragma acc declare create ( _cold_source )

_class_ESS_moderator_long_2001 _thermal_source_var;
_class_ESS_moderator_long_2001 *_thermal_source = &_thermal_source_var;
#pragma acc declare create ( _thermal_source_var )
#pragma acc declare create ( _thermal_source )

_class_Arm _ColdFocus_var;
_class_Arm *_ColdFocus = &_ColdFocus_var;
#pragma acc declare create ( _ColdFocus_var )
#pragma acc declare create ( _ColdFocus )

_class_Arm _ArmMidOne_var;
_class_Arm *_ArmMidOne = &_ArmMidOne_var;
#pragma acc declare create ( _ArmMidOne_var )
#pragma acc declare create ( _ArmMidOne )

/* component mirror_full_center=Mirror_Curved_Bispectral() [10] DECLARE */
/* Parameter definition for component type 'Mirror_Curved_Bispectral' */
struct _struct_Mirror_Curved_Bispectral_parameters {
  /* Component type 'Mirror_Curved_Bispectral' setting parameters */
  char reflect[16384];
  MCNUM focus_s;
  MCNUM focus_e;
  MCNUM mirror_start;
  MCNUM guide_start;
  MCNUM yheight;
  MCNUM smallaxis;
  MCNUM length;
  MCNUM m;
  MCNUM transmit;
  MCNUM substrate_thickness;
  MCNUM coating_thickness;
  MCNUM theta_1;
  MCNUM theta_2;
  MCNUM theta_3;
  /* Component type 'Mirror_Curved_Bispectral' private parameters */
  /* Component type 'Mirror_Curved_Bispectral' DECLARE code stored as structure members */
t_Table pTable;

}; /* _struct_Mirror_Curved_Bispectral_parameters */
typedef struct _struct_Mirror_Curved_Bispectral_parameters _class_Mirror_Curved_Bispectral_parameters;

/* Parameters for component type 'Mirror_Curved_Bispectral' */
struct _struct_Mirror_Curved_Bispectral {
  char     _name[256]; /* e.g. mirror_full_center */
  char     _type[256]; /* Mirror_Curved_Bispectral */
  long     _index; /* e.g. 2 index in TRACE list */
  Coords   _position_absolute;
  Coords   _position_relative; /* wrt PREVIOUS */
  Rotation _rotation_absolute;
  Rotation _rotation_relative; /* wrt PREVIOUS */
  int      _rotation_is_identity;
  _class_Mirror_Curved_Bispectral_parameters _parameters;
};
typedef struct _struct_Mirror_Curved_Bispectral _class_Mirror_Curved_Bispectral;
_class_Mirror_Curved_Bispectral _mirror_full_center_var;
_class_Mirror_Curved_Bispectral *_mirror_full_center = &_mirror_full_center_var;
#pragma acc declare create ( _mirror_full_center_var )
#pragma acc declare create ( _mirror_full_center )

_class_Arm _ArmForNeutronPropState_2_var;
_class_Arm *_ArmForNeutronPropState_2 = &_ArmForNeutronPropState_2_var;
#pragma acc declare create ( _ArmForNeutronPropState_2_var )
#pragma acc declare create ( _ArmForNeutronPropState_2 )

/* component guide_right=Mirror_Elliptic_Bispectral() [12] DECLARE */
/* Parameter definition for component type 'Mirror_Elliptic_Bispectral' */
struct _struct_Mirror_Elliptic_Bispectral_parameters {
  /* Component type 'Mirror_Elliptic_Bispectral' setting parameters */
  char reflect[16384];
  MCNUM focus_start_w;
  MCNUM focus_end_w;
  MCNUM focus_start_h;
  MCNUM focus_end_h;
  MCNUM mirror_start;
  MCNUM m;
  MCNUM smallaxis_w;
  MCNUM smallaxis_h;
  MCNUM length;
  MCNUM transmit;
  MCNUM substrate_thickness;
  MCNUM coating_thickness;
  /* Component type 'Mirror_Elliptic_Bispectral' private parameters */
  /* Component type 'Mirror_Elliptic_Bispectral' DECLARE code stored as structure members */
  t_Table pTable;
}; /* _struct_Mirror_Elliptic_Bispectral_parameters */
typedef struct _struct_Mirror_Elliptic_Bispectral_parameters _class_Mirror_Elliptic_Bispectral_parameters;

/* Parameters for component type 'Mirror_Elliptic_Bispectral' */
struct _struct_Mirror_Elliptic_Bispectral {
  char     _name[256]; /* e.g. guide_right */
  char     _type[256]; /* Mirror_Elliptic_Bispectral */
  long     _index; /* e.g. 2 index in TRACE list */
  Coords   _position_absolute;
  Coords   _position_relative; /* wrt PREVIOUS */
  Rotation _rotation_absolute;
  Rotation _rotation_relative; /* wrt PREVIOUS */
  int      _rotation_is_identity;
  _class_Mirror_Elliptic_Bispectral_parameters _parameters;
};
typedef struct _struct_Mirror_Elliptic_Bispectral _class_Mirror_Elliptic_Bispectral;
_class_Mirror_Elliptic_Bispectral _guide_right_var;
_class_Mirror_Elliptic_Bispectral *_guide_right = &_guide_right_var;
#pragma acc declare create ( _guide_right_var )
#pragma acc declare create ( _guide_right )

_class_Arm _ArmForNeutronPropState_4_var;
_class_Arm *_ArmForNeutronPropState_4 = &_ArmForNeutronPropState_4_var;
#pragma acc declare create ( _ArmForNeutronPropState_4_var )
#pragma acc declare create ( _ArmForNeutronPropState_4 )

_class_Mirror_Elliptic_Bispectral _guide_bottom_var;
_class_Mirror_Elliptic_Bispectral *_guide_bottom = &_guide_bottom_var;
#pragma acc declare create ( _guide_bottom_var )
#pragma acc declare create ( _guide_bottom )

_class_Arm _ArmForNeutronPropState_5_var;
_class_Arm *_ArmForNeutronPropState_5 = &_ArmForNeutronPropState_5_var;
#pragma acc declare create ( _ArmForNeutronPropState_5_var )
#pragma acc declare create ( _ArmForNeutronPropState_5 )

/* component cold_lambda_guidestart=L_monitor() [16] DECLARE */
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
  char     _name[256]; /* e.g. cold_lambda_guidestart */
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
_class_L_monitor _cold_lambda_guidestart_var;
_class_L_monitor *_cold_lambda_guidestart = &_cold_lambda_guidestart_var;
#pragma acc declare create ( _cold_lambda_guidestart_var )
#pragma acc declare create ( _cold_lambda_guidestart )

_class_L_monitor _thermal_lambda_guidestart_var;
_class_L_monitor *_thermal_lambda_guidestart = &_thermal_lambda_guidestart_var;
#pragma acc declare create ( _thermal_lambda_guidestart_var )
#pragma acc declare create ( _thermal_lambda_guidestart )

_class_L_monitor _lambda_guidestart_var;
_class_L_monitor *_lambda_guidestart = &_lambda_guidestart_var;
#pragma acc declare create ( _lambda_guidestart_var )
#pragma acc declare create ( _lambda_guidestart )

_class_Mirror_Elliptic_Bispectral _guide_top_var;
_class_Mirror_Elliptic_Bispectral *_guide_top = &_guide_top_var;
#pragma acc declare create ( _guide_top_var )
#pragma acc declare create ( _guide_top )

_class_Arm _ArmForNeutronPropState_6_var;
_class_Arm *_ArmForNeutronPropState_6 = &_ArmForNeutronPropState_6_var;
#pragma acc declare create ( _ArmForNeutronPropState_6_var )
#pragma acc declare create ( _ArmForNeutronPropState_6 )

_class_Mirror_Elliptic_Bispectral _guide_Left_var;
_class_Mirror_Elliptic_Bispectral *_guide_Left = &_guide_Left_var;
#pragma acc declare create ( _guide_Left_var )
#pragma acc declare create ( _guide_Left )

_class_Arm _ArmForNeutronPropState_7_var;
_class_Arm *_ArmForNeutronPropState_7 = &_ArmForNeutronPropState_7_var;
#pragma acc declare create ( _ArmForNeutronPropState_7_var )
#pragma acc declare create ( _ArmForNeutronPropState_7 )

_class_Arm _ArmMidTwo_var;
_class_Arm *_ArmMidTwo = &_ArmMidTwo_var;
#pragma acc declare create ( _ArmMidTwo_var )
#pragma acc declare create ( _ArmMidTwo )

_class_Arm _ArmForNeutronPropState_8_var;
_class_Arm *_ArmForNeutronPropState_8 = &_ArmForNeutronPropState_8_var;
#pragma acc declare create ( _ArmForNeutronPropState_8_var )
#pragma acc declare create ( _ArmForNeutronPropState_8 )

_class_Arm _ArmMidThree_var;
_class_Arm *_ArmMidThree = &_ArmMidThree_var;
#pragma acc declare create ( _ArmMidThree_var )
#pragma acc declare create ( _ArmMidThree )

_class_Arm _ArmExit_var;
_class_Arm *_ArmExit = &_ArmExit_var;
#pragma acc declare create ( _ArmExit_var )
#pragma acc declare create ( _ArmExit )

int mcNUMCOMP = 26;

/* User declarations from instrument definition. Can define functions. */
double cold_moderator_x=0.12;

// parameters for the elliptic guide
double distance[202], length[201], h[201], w[201], alpha[201], guide_piece_length[201];
double guide_dist, Pulse_freq;
int n_elements=50;
//%include "mdeclare.c"
double m[50];
double guide_length;

double focus_s_w;
double focus_e_w;
double focus_s_h;
double focus_e_h;

double W_par=0.003;
double R0_par=0.99;
double waviness=0.0101;


// The standard wavelengths when man_lam=0
double thermal_lam_max;
double thermal_lam_min;
double cold_lam_max;
double cold_lam_min;

// parameters to switch between sources
double dummy;
int two_sources;
int flag;
double coldperthermal;
double coldmultiplier;
double thermalmultiplier;

double x_one; //lower focus point for source
double x_two; //upper focus point for source
int r; //number of first guide element that is outside mirror
double x_mid; //center to focus on
double x_focus; //width of focusing rectangle
double y_focus;

// Parameters to find the origin of neutrons hitting the sample
double xmidlertidig;
double ymidlertidig;
double xcold;
double xthermal;
double ycold;
double ythermal;
int hit_sample_flag=0;

// Parameters to remove the too divergent neutrons from the monitors
double x_div;
double y_div;

//to remove neutrons outside the wavelength band
double lambda;

// to find where the simulated rays come from
double p_old;

// Parameters for hotspot
double size=0.12;
double thermal_hotspot_dia= 0.03;
double thermal_hotspot_factor;
double thermal_hotspot_x_center=-0.01;
double thermal_hotspot_y_center=0.0;

double cold_hotspot_dia= 0.03;
double cold_hotspot_factor;
double cold_hotspot_x_center=0.01;
double cold_hotspot_y_center=0.03;

// double cold_hotspot_xwidth=0.018;
// double cold_hotspot_yheight=0.039;




// Parametes for the mirror
//---------------------------------------------------------------
double L_moderator_mirror=3.25;   //distance to center of mirror
double mirror_full_length=2.5; //length of mirror
double extraction_start_pos=2.0; //where the mirror starts
double guide_left_start=3.970; //where the left part of the guide starts
double R0_mirror=0.99;

double W_mirror=0.003;

double alpha_mirror;

int n_elements_mirror=16;
int mirror_part_in_guide_number[51];
double mirror_rot[50];
double x_mirror[50];
double mirror_part_length[50];
double mirror_height_in_guide[51];
double y_mirror[50];
double z_mirror[50];
double mirror_rotation[50]; 

double mirror_start;
double mirror_end;
double h_mirror_part;
double L_mirror_part;


int k;
int h_index;

//Rotation of mirror is theta(z)=a*z^2+b*z+c, m is a_m*z^2+b_m*z+c_m
double x1, x2, x3;
double a, b, c;
double a_m, b_m, c_m;

//end of parameters for the mirror
//-------------------------------------------------


// for the guide made of mirrors
int guide_scatt=0;
double guide_bottom_height[202];
double guide_bottom_rotation[202];
double guide_right_height[202];
double guide_right_rotation[202];
double guide_h_pos[202];
double guide_w_pos[202];
double guide_z_pos[202];
double guide_right_rot[202];
double guide_bottom_rot[202];
double guide_top_rot[202];
double guide_left_rot[202];
int use_guidegravity_flag[202];
int use_guide_left_part[202];

double ArmExitPos;
double ArmMidOnePos;





double old_x_prop;
double old_y_prop;
double old_z_prop;

double old_vx_prop;
double old_vy_prop;
double old_vz_prop;

double old_t_prop;
double old_p_prop;

double new_x_prop;
double new_y_prop;
double new_z_prop;

double new_vx_prop;
double new_vy_prop;
double new_vz_prop;

double new_t_prop;
double new_p_prop;



double w_extractionstart;
double h_extractionstart;
double w_guide_leftstart;
double h_guide_leftstart;



#undef compcurname
#undef compcurtype
#undef compcurindex
/* end of instrument 'ESS_2001_bispectral' and components DECLARE */

/* *****************************************************************************
* instrument 'ESS_2001_bispectral' and components INITIALISE
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

/* component ArmForGuideRight=Arm() SETTING, POSITION/ROTATION */
int _ArmForGuideRight_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_ArmForGuideRight_setpos] component ArmForGuideRight=Arm() SETTING [Arm:0]");
  stracpy(_ArmForGuideRight->_name, "ArmForGuideRight", 16384);
  stracpy(_ArmForGuideRight->_type, "Arm", 16384);
  _ArmForGuideRight->_index=2;
  /* component ArmForGuideRight=Arm() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (0)*DEG2RAD, (90)*DEG2RAD);
    rot_mul(tr1, _Origin->_rotation_absolute, _ArmForGuideRight->_rotation_absolute);
    rot_transpose(_Origin->_rotation_absolute, tr1);
    rot_mul(_ArmForGuideRight->_rotation_absolute, tr1, _ArmForGuideRight->_rotation_relative);
    _ArmForGuideRight->_rotation_is_identity =  rot_test_identity(_ArmForGuideRight->_rotation_relative);
    tc1 = coords_set(
      0, 0, 0);
    rot_transpose(_Origin->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _ArmForGuideRight->_position_absolute = coords_add(_Origin->_position_absolute, tc2);
    tc1 = coords_sub(_Origin->_position_absolute, _ArmForGuideRight->_position_absolute);
    _ArmForGuideRight->_position_relative = rot_apply(_ArmForGuideRight->_rotation_absolute, tc1);
  } /* ArmForGuideRight=Arm() AT ROTATED */
  DEBUG_COMPONENT("ArmForGuideRight", _ArmForGuideRight->_position_absolute, _ArmForGuideRight->_rotation_absolute);
  instrument->_position_absolute[2] = _ArmForGuideRight->_position_absolute;
  instrument->_position_relative[2] = _ArmForGuideRight->_position_relative;
  instrument->counter_N[2]  = instrument->counter_P[2] = instrument->counter_P2[2] = 0;
  instrument->counter_AbsorbProp[2]= 0;
  return(0);
} /* _ArmForGuideRight_setpos */

/* component ArmForGuideBottom=Arm() SETTING, POSITION/ROTATION */
int _ArmForGuideBottom_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_ArmForGuideBottom_setpos] component ArmForGuideBottom=Arm() SETTING [Arm:0]");
  stracpy(_ArmForGuideBottom->_name, "ArmForGuideBottom", 16384);
  stracpy(_ArmForGuideBottom->_type, "Arm", 16384);
  _ArmForGuideBottom->_index=3;
  /* component ArmForGuideBottom=Arm() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (0)*DEG2RAD, (90)*DEG2RAD);
    rot_mul(tr1, _Origin->_rotation_absolute, _ArmForGuideBottom->_rotation_absolute);
    rot_transpose(_ArmForGuideRight->_rotation_absolute, tr1);
    rot_mul(_ArmForGuideBottom->_rotation_absolute, tr1, _ArmForGuideBottom->_rotation_relative);
    _ArmForGuideBottom->_rotation_is_identity =  rot_test_identity(_ArmForGuideBottom->_rotation_relative);
    tc1 = coords_set(
      0, 0, 0);
    rot_transpose(_Origin->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _ArmForGuideBottom->_position_absolute = coords_add(_Origin->_position_absolute, tc2);
    tc1 = coords_sub(_ArmForGuideRight->_position_absolute, _ArmForGuideBottom->_position_absolute);
    _ArmForGuideBottom->_position_relative = rot_apply(_ArmForGuideBottom->_rotation_absolute, tc1);
  } /* ArmForGuideBottom=Arm() AT ROTATED */
  DEBUG_COMPONENT("ArmForGuideBottom", _ArmForGuideBottom->_position_absolute, _ArmForGuideBottom->_rotation_absolute);
  instrument->_position_absolute[3] = _ArmForGuideBottom->_position_absolute;
  instrument->_position_relative[3] = _ArmForGuideBottom->_position_relative;
  instrument->counter_N[3]  = instrument->counter_P[3] = instrument->counter_P2[3] = 0;
  instrument->counter_AbsorbProp[3]= 0;
  return(0);
} /* _ArmForGuideBottom_setpos */

/* component ArmForGuideTop=Arm() SETTING, POSITION/ROTATION */
int _ArmForGuideTop_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_ArmForGuideTop_setpos] component ArmForGuideTop=Arm() SETTING [Arm:0]");
  stracpy(_ArmForGuideTop->_name, "ArmForGuideTop", 16384);
  stracpy(_ArmForGuideTop->_type, "Arm", 16384);
  _ArmForGuideTop->_index=4;
  /* component ArmForGuideTop=Arm() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (0)*DEG2RAD, (-90)*DEG2RAD);
    rot_mul(tr1, _Origin->_rotation_absolute, _ArmForGuideTop->_rotation_absolute);
    rot_transpose(_ArmForGuideBottom->_rotation_absolute, tr1);
    rot_mul(_ArmForGuideTop->_rotation_absolute, tr1, _ArmForGuideTop->_rotation_relative);
    _ArmForGuideTop->_rotation_is_identity =  rot_test_identity(_ArmForGuideTop->_rotation_relative);
    tc1 = coords_set(
      0, 0, 0);
    rot_transpose(_Origin->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _ArmForGuideTop->_position_absolute = coords_add(_Origin->_position_absolute, tc2);
    tc1 = coords_sub(_ArmForGuideBottom->_position_absolute, _ArmForGuideTop->_position_absolute);
    _ArmForGuideTop->_position_relative = rot_apply(_ArmForGuideTop->_rotation_absolute, tc1);
  } /* ArmForGuideTop=Arm() AT ROTATED */
  DEBUG_COMPONENT("ArmForGuideTop", _ArmForGuideTop->_position_absolute, _ArmForGuideTop->_rotation_absolute);
  instrument->_position_absolute[4] = _ArmForGuideTop->_position_absolute;
  instrument->_position_relative[4] = _ArmForGuideTop->_position_relative;
  instrument->counter_N[4]  = instrument->counter_P[4] = instrument->counter_P2[4] = 0;
  instrument->counter_AbsorbProp[4]= 0;
  return(0);
} /* _ArmForGuideTop_setpos */

/* component ArmForGuideLeft=Arm() SETTING, POSITION/ROTATION */
int _ArmForGuideLeft_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_ArmForGuideLeft_setpos] component ArmForGuideLeft=Arm() SETTING [Arm:0]");
  stracpy(_ArmForGuideLeft->_name, "ArmForGuideLeft", 16384);
  stracpy(_ArmForGuideLeft->_type, "Arm", 16384);
  _ArmForGuideLeft->_index=5;
  /* component ArmForGuideLeft=Arm() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (0)*DEG2RAD, (180)*DEG2RAD);
    rot_mul(tr1, _Origin->_rotation_absolute, _ArmForGuideLeft->_rotation_absolute);
    rot_transpose(_ArmForGuideTop->_rotation_absolute, tr1);
    rot_mul(_ArmForGuideLeft->_rotation_absolute, tr1, _ArmForGuideLeft->_rotation_relative);
    _ArmForGuideLeft->_rotation_is_identity =  rot_test_identity(_ArmForGuideLeft->_rotation_relative);
    tc1 = coords_set(
      0, 0, 0);
    rot_transpose(_Origin->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _ArmForGuideLeft->_position_absolute = coords_add(_Origin->_position_absolute, tc2);
    tc1 = coords_sub(_ArmForGuideTop->_position_absolute, _ArmForGuideLeft->_position_absolute);
    _ArmForGuideLeft->_position_relative = rot_apply(_ArmForGuideLeft->_rotation_absolute, tc1);
  } /* ArmForGuideLeft=Arm() AT ROTATED */
  DEBUG_COMPONENT("ArmForGuideLeft", _ArmForGuideLeft->_position_absolute, _ArmForGuideLeft->_rotation_absolute);
  instrument->_position_absolute[5] = _ArmForGuideLeft->_position_absolute;
  instrument->_position_relative[5] = _ArmForGuideLeft->_position_relative;
  instrument->counter_N[5]  = instrument->counter_P[5] = instrument->counter_P2[5] = 0;
  instrument->counter_AbsorbProp[5]= 0;
  return(0);
} /* _ArmForGuideLeft_setpos */

/* component cold_source=ESS_moderator_long_2001() SETTING, POSITION/ROTATION */
int _cold_source_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_cold_source_setpos] component cold_source=ESS_moderator_long_2001() SETTING [/usr/share/mcstas/3.0-dev/obsolete/ESS_moderator_long_2001.comp:114]");
  stracpy(_cold_source->_name, "cold_source", 16384);
  stracpy(_cold_source->_type, "ESS_moderator_long_2001", 16384);
  _cold_source->_index=6;
  _cold_source->_parameters.size = 0.12;
  #define size (_cold_source->_parameters.size)
  _cold_source->_parameters.l_low = instrument->_parameters._cold_Lam_min;
  #define l_low (_cold_source->_parameters.l_low)
  _cold_source->_parameters.l_high = instrument->_parameters._cold_Lam_max;
  #define l_high (_cold_source->_parameters.l_high)
  _cold_source->_parameters.dist = 0;
  #define dist (_cold_source->_parameters.dist)
  _cold_source->_parameters.xw = x_focus + 0.15;
  #define xw (_cold_source->_parameters.xw)
  _cold_source->_parameters.yh = y_focus + 0.15;
  #define yh (_cold_source->_parameters.yh)
  _cold_source->_parameters.freq = instrument->_parameters._frequency;
  #define freq (_cold_source->_parameters.freq)
  _cold_source->_parameters.T = 50;
  #define T (_cold_source->_parameters.T)
  _cold_source->_parameters.tau = 287e-6;
  #define tau (_cold_source->_parameters.tau)
  _cold_source->_parameters.tau1 = 0;
  #define tau1 (_cold_source->_parameters.tau1)
  _cold_source->_parameters.tau2 = 20e-6;
  #define tau2 (_cold_source->_parameters.tau2)
  _cold_source->_parameters.d = instrument->_parameters._Pulse_width;
  #define d (_cold_source->_parameters.d)
  _cold_source->_parameters.n = 20;
  #define n (_cold_source->_parameters.n)
  _cold_source->_parameters.n2 = 5;
  #define n2 (_cold_source->_parameters.n2)
  _cold_source->_parameters.chi2 = 0.9;
  #define chi2 (_cold_source->_parameters.chi2)
  _cold_source->_parameters.I0 = 6.9e11;
  #define I0 (_cold_source->_parameters.I0)
  _cold_source->_parameters.I2 = 27.6e10;
  #define I2 (_cold_source->_parameters.I2)
  _cold_source->_parameters.branch1 = 1.0;
  #define branch1 (_cold_source->_parameters.branch1)
  _cold_source->_parameters.branch2 = 0.5;
  #define branch2 (_cold_source->_parameters.branch2)
  _cold_source->_parameters.branch_tail = 0.1;
  #define branch_tail (_cold_source->_parameters.branch_tail)
  _cold_source->_parameters.twopulses = 0;
  #define twopulses (_cold_source->_parameters.twopulses)
  _cold_source->_parameters.target_index = 2;
  #define target_index (_cold_source->_parameters.target_index)

  #define l_range (_cold_source->_parameters.l_range)
  #define w_mult (_cold_source->_parameters.w_mult)
  #define branchframe (_cold_source->_parameters.branchframe)
  #define tx (_cold_source->_parameters.tx)
  #define ty (_cold_source->_parameters.ty)
  #define tz (_cold_source->_parameters.tz)

  #undef size
  #undef l_low
  #undef l_high
  #undef dist
  #undef xw
  #undef yh
  #undef freq
  #undef T
  #undef tau
  #undef tau1
  #undef tau2
  #undef d
  #undef n
  #undef n2
  #undef chi2
  #undef I0
  #undef I2
  #undef branch1
  #undef branch2
  #undef branch_tail
  #undef twopulses
  #undef target_index
  #undef l_range
  #undef w_mult
  #undef branchframe
  #undef tx
  #undef ty
  #undef tz
  /* component cold_source=ESS_moderator_long_2001() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Origin->_rotation_absolute, _cold_source->_rotation_absolute);
    rot_transpose(_ArmForGuideLeft->_rotation_absolute, tr1);
    rot_mul(_cold_source->_rotation_absolute, tr1, _cold_source->_rotation_relative);
    _cold_source->_rotation_is_identity =  rot_test_identity(_cold_source->_rotation_relative);
    tc1 = coords_set(
      cold_moderator_x, 0, 0);
    rot_transpose(_Origin->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _cold_source->_position_absolute = coords_add(_Origin->_position_absolute, tc2);
    tc1 = coords_sub(_ArmForGuideLeft->_position_absolute, _cold_source->_position_absolute);
    _cold_source->_position_relative = rot_apply(_cold_source->_rotation_absolute, tc1);
  } /* cold_source=ESS_moderator_long_2001() AT ROTATED */
  DEBUG_COMPONENT("cold_source", _cold_source->_position_absolute, _cold_source->_rotation_absolute);
  instrument->_position_absolute[6] = _cold_source->_position_absolute;
  instrument->_position_relative[6] = _cold_source->_position_relative;
  instrument->counter_N[6]  = instrument->counter_P[6] = instrument->counter_P2[6] = 0;
  instrument->counter_AbsorbProp[6]= 0;
  return(0);
} /* _cold_source_setpos */

/* component thermal_source=ESS_moderator_long_2001() SETTING, POSITION/ROTATION */
int _thermal_source_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_thermal_source_setpos] component thermal_source=ESS_moderator_long_2001() SETTING [/usr/share/mcstas/3.0-dev/obsolete/ESS_moderator_long_2001.comp:114]");
  stracpy(_thermal_source->_name, "thermal_source", 16384);
  stracpy(_thermal_source->_type, "ESS_moderator_long_2001", 16384);
  _thermal_source->_index=7;
  _thermal_source->_parameters.size = 0.12;
  #define size (_thermal_source->_parameters.size)
  _thermal_source->_parameters.l_low = instrument->_parameters._thermal_Lam_min;
  #define l_low (_thermal_source->_parameters.l_low)
  _thermal_source->_parameters.l_high = instrument->_parameters._thermal_Lam_max;
  #define l_high (_thermal_source->_parameters.l_high)
  _thermal_source->_parameters.dist = extraction_start_pos;
  #define dist (_thermal_source->_parameters.dist)
  _thermal_source->_parameters.xw = 2 * w_guide_leftstart + 0.1;
  #define xw (_thermal_source->_parameters.xw)
  _thermal_source->_parameters.yh = 2 * h_guide_leftstart + 0.1;
  #define yh (_thermal_source->_parameters.yh)
  _thermal_source->_parameters.freq = instrument->_parameters._frequency;
  #define freq (_thermal_source->_parameters.freq)
  _thermal_source->_parameters.T = 325;
  #define T (_thermal_source->_parameters.T)
  _thermal_source->_parameters.tau = 80e-6;
  #define tau (_thermal_source->_parameters.tau)
  _thermal_source->_parameters.tau1 = 400e-6;
  #define tau1 (_thermal_source->_parameters.tau1)
  _thermal_source->_parameters.tau2 = 12e-6;
  #define tau2 (_thermal_source->_parameters.tau2)
  _thermal_source->_parameters.d = instrument->_parameters._Pulse_width;
  #define d (_thermal_source->_parameters.d)
  _thermal_source->_parameters.n = 20;
  #define n (_thermal_source->_parameters.n)
  _thermal_source->_parameters.n2 = 5;
  #define n2 (_thermal_source->_parameters.n2)
  _thermal_source->_parameters.chi2 = 2.5;
  #define chi2 (_thermal_source->_parameters.chi2)
  _thermal_source->_parameters.I0 = 13.5e11;
  #define I0 (_thermal_source->_parameters.I0)
  _thermal_source->_parameters.I2 = 27.6e10;
  #define I2 (_thermal_source->_parameters.I2)
  _thermal_source->_parameters.branch1 = 0.5;
  #define branch1 (_thermal_source->_parameters.branch1)
  _thermal_source->_parameters.branch2 = 0.5;
  #define branch2 (_thermal_source->_parameters.branch2)
  _thermal_source->_parameters.branch_tail = 0.1;
  #define branch_tail (_thermal_source->_parameters.branch_tail)
  _thermal_source->_parameters.twopulses = 0;
  #define twopulses (_thermal_source->_parameters.twopulses)
  _thermal_source->_parameters.target_index = 0;
  #define target_index (_thermal_source->_parameters.target_index)

  #define l_range (_thermal_source->_parameters.l_range)
  #define w_mult (_thermal_source->_parameters.w_mult)
  #define branchframe (_thermal_source->_parameters.branchframe)
  #define tx (_thermal_source->_parameters.tx)
  #define ty (_thermal_source->_parameters.ty)
  #define tz (_thermal_source->_parameters.tz)

  #undef size
  #undef l_low
  #undef l_high
  #undef dist
  #undef xw
  #undef yh
  #undef freq
  #undef T
  #undef tau
  #undef tau1
  #undef tau2
  #undef d
  #undef n
  #undef n2
  #undef chi2
  #undef I0
  #undef I2
  #undef branch1
  #undef branch2
  #undef branch_tail
  #undef twopulses
  #undef target_index
  #undef l_range
  #undef w_mult
  #undef branchframe
  #undef tx
  #undef ty
  #undef tz
  /* component thermal_source=ESS_moderator_long_2001() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Origin->_rotation_absolute, _thermal_source->_rotation_absolute);
    rot_transpose(_cold_source->_rotation_absolute, tr1);
    rot_mul(_thermal_source->_rotation_absolute, tr1, _thermal_source->_rotation_relative);
    _thermal_source->_rotation_is_identity =  rot_test_identity(_thermal_source->_rotation_relative);
    tc1 = coords_set(
      0, 0, 0);
    rot_transpose(_Origin->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _thermal_source->_position_absolute = coords_add(_Origin->_position_absolute, tc2);
    tc1 = coords_sub(_cold_source->_position_absolute, _thermal_source->_position_absolute);
    _thermal_source->_position_relative = rot_apply(_thermal_source->_rotation_absolute, tc1);
  } /* thermal_source=ESS_moderator_long_2001() AT ROTATED */
  DEBUG_COMPONENT("thermal_source", _thermal_source->_position_absolute, _thermal_source->_rotation_absolute);
  instrument->_position_absolute[7] = _thermal_source->_position_absolute;
  instrument->_position_relative[7] = _thermal_source->_position_relative;
  instrument->counter_N[7]  = instrument->counter_P[7] = instrument->counter_P2[7] = 0;
  instrument->counter_AbsorbProp[7]= 0;
  return(0);
} /* _thermal_source_setpos */

/* component ColdFocus=Arm() SETTING, POSITION/ROTATION */
int _ColdFocus_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_ColdFocus_setpos] component ColdFocus=Arm() SETTING [Arm:0]");
  stracpy(_ColdFocus->_name, "ColdFocus", 16384);
  stracpy(_ColdFocus->_type, "Arm", 16384);
  _ColdFocus->_index=8;
  /* component ColdFocus=Arm() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Origin->_rotation_absolute, _ColdFocus->_rotation_absolute);
    rot_transpose(_thermal_source->_rotation_absolute, tr1);
    rot_mul(_ColdFocus->_rotation_absolute, tr1, _ColdFocus->_rotation_relative);
    _ColdFocus->_rotation_is_identity =  rot_test_identity(_ColdFocus->_rotation_relative);
    tc1 = coords_set(
      x_mid, 0, extraction_start_pos);
    rot_transpose(_Origin->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _ColdFocus->_position_absolute = coords_add(_Origin->_position_absolute, tc2);
    tc1 = coords_sub(_thermal_source->_position_absolute, _ColdFocus->_position_absolute);
    _ColdFocus->_position_relative = rot_apply(_ColdFocus->_rotation_absolute, tc1);
  } /* ColdFocus=Arm() AT ROTATED */
  DEBUG_COMPONENT("ColdFocus", _ColdFocus->_position_absolute, _ColdFocus->_rotation_absolute);
  instrument->_position_absolute[8] = _ColdFocus->_position_absolute;
  instrument->_position_relative[8] = _ColdFocus->_position_relative;
  instrument->counter_N[8]  = instrument->counter_P[8] = instrument->counter_P2[8] = 0;
  instrument->counter_AbsorbProp[8]= 0;
  return(0);
} /* _ColdFocus_setpos */

/* component ArmMidOne=Arm() SETTING, POSITION/ROTATION */
int _ArmMidOne_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_ArmMidOne_setpos] component ArmMidOne=Arm() SETTING [Arm:0]");
  stracpy(_ArmMidOne->_name, "ArmMidOne", 16384);
  stracpy(_ArmMidOne->_type, "Arm", 16384);
  _ArmMidOne->_index=9;
  /* component ArmMidOne=Arm() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Origin->_rotation_absolute, _ArmMidOne->_rotation_absolute);
    rot_transpose(_ColdFocus->_rotation_absolute, tr1);
    rot_mul(_ArmMidOne->_rotation_absolute, tr1, _ArmMidOne->_rotation_relative);
    _ArmMidOne->_rotation_is_identity =  rot_test_identity(_ArmMidOne->_rotation_relative);
    tc1 = coords_set(
      0, 0, ArmMidOnePos);
    rot_transpose(_Origin->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _ArmMidOne->_position_absolute = coords_add(_Origin->_position_absolute, tc2);
    tc1 = coords_sub(_ColdFocus->_position_absolute, _ArmMidOne->_position_absolute);
    _ArmMidOne->_position_relative = rot_apply(_ArmMidOne->_rotation_absolute, tc1);
  } /* ArmMidOne=Arm() AT ROTATED */
  DEBUG_COMPONENT("ArmMidOne", _ArmMidOne->_position_absolute, _ArmMidOne->_rotation_absolute);
  instrument->_position_absolute[9] = _ArmMidOne->_position_absolute;
  instrument->_position_relative[9] = _ArmMidOne->_position_relative;
  instrument->counter_N[9]  = instrument->counter_P[9] = instrument->counter_P2[9] = 0;
  instrument->counter_AbsorbProp[9]= 0;
  return(0);
} /* _ArmMidOne_setpos */

/* component mirror_full_center=Mirror_Curved_Bispectral() SETTING, POSITION/ROTATION */
int _mirror_full_center_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_mirror_full_center_setpos] component mirror_full_center=Mirror_Curved_Bispectral() SETTING [/usr/share/mcstas/3.0-dev/contrib/Mirror_Curved_Bispectral.comp:65]");
  stracpy(_mirror_full_center->_name, "mirror_full_center", 16384);
  stracpy(_mirror_full_center->_type, "Mirror_Curved_Bispectral", 16384);
  _mirror_full_center->_index=10;
  _mirror_full_center->_parameters.reflect[0]='\0';
  #define reflect (_mirror_full_center->_parameters.reflect)
  _mirror_full_center->_parameters.focus_s = instrument->_parameters._focus_start_h;
  #define focus_s (_mirror_full_center->_parameters.focus_s)
  _mirror_full_center->_parameters.focus_e = instrument->_parameters._focus_end_h;
  #define focus_e (_mirror_full_center->_parameters.focus_e)
  _mirror_full_center->_parameters.mirror_start = extraction_start_pos;
  #define mirror_start (_mirror_full_center->_parameters.mirror_start)
  _mirror_full_center->_parameters.guide_start = instrument->_parameters._guide_start;
  #define guide_start (_mirror_full_center->_parameters.guide_start)
  _mirror_full_center->_parameters.yheight = instrument->_parameters._h_mirror;
  #define yheight (_mirror_full_center->_parameters.yheight)
  _mirror_full_center->_parameters.smallaxis = instrument->_parameters._smallaxis_h;
  #define smallaxis (_mirror_full_center->_parameters.smallaxis)
  _mirror_full_center->_parameters.length = 2.5;
  #define length (_mirror_full_center->_parameters.length)
  _mirror_full_center->_parameters.m = instrument->_parameters._m_mirror;
  #define m (_mirror_full_center->_parameters.m)
  _mirror_full_center->_parameters.transmit = 1;
  #define transmit (_mirror_full_center->_parameters.transmit)
  _mirror_full_center->_parameters.substrate_thickness = instrument->_parameters._substrate_thickness;
  #define substrate_thickness (_mirror_full_center->_parameters.substrate_thickness)
  _mirror_full_center->_parameters.coating_thickness = instrument->_parameters._coating_thickness;
  #define coating_thickness (_mirror_full_center->_parameters.coating_thickness)
  _mirror_full_center->_parameters.theta_1 = instrument->_parameters._theta1;
  #define theta_1 (_mirror_full_center->_parameters.theta_1)
  _mirror_full_center->_parameters.theta_2 = instrument->_parameters._theta1;
  #define theta_2 (_mirror_full_center->_parameters.theta_2)
  _mirror_full_center->_parameters.theta_3 = instrument->_parameters._theta1;
  #define theta_3 (_mirror_full_center->_parameters.theta_3)

  #define pTable (_mirror_full_center->_parameters.pTable)

  #undef reflect
  #undef focus_s
  #undef focus_e
  #undef mirror_start
  #undef guide_start
  #undef yheight
  #undef smallaxis
  #undef length
  #undef m
  #undef transmit
  #undef substrate_thickness
  #undef coating_thickness
  #undef theta_1
  #undef theta_2
  #undef theta_3
  #undef pTable
  /* component mirror_full_center=Mirror_Curved_Bispectral() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (-90)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _Origin->_rotation_absolute, _mirror_full_center->_rotation_absolute);
    rot_transpose(_ArmMidOne->_rotation_absolute, tr1);
    rot_mul(_mirror_full_center->_rotation_absolute, tr1, _mirror_full_center->_rotation_relative);
    _mirror_full_center->_rotation_is_identity =  rot_test_identity(_mirror_full_center->_rotation_relative);
    tc1 = coords_set(
      instrument->_parameters._mirror_offset, 0, extraction_start_pos + mirror_full_length * 0.5);
    rot_transpose(_Origin->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _mirror_full_center->_position_absolute = coords_add(_Origin->_position_absolute, tc2);
    tc1 = coords_sub(_ArmMidOne->_position_absolute, _mirror_full_center->_position_absolute);
    _mirror_full_center->_position_relative = rot_apply(_mirror_full_center->_rotation_absolute, tc1);
  } /* mirror_full_center=Mirror_Curved_Bispectral() AT ROTATED */
  DEBUG_COMPONENT("mirror_full_center", _mirror_full_center->_position_absolute, _mirror_full_center->_rotation_absolute);
  instrument->_position_absolute[10] = _mirror_full_center->_position_absolute;
  instrument->_position_relative[10] = _mirror_full_center->_position_relative;
  instrument->counter_N[10]  = instrument->counter_P[10] = instrument->counter_P2[10] = 0;
  instrument->counter_AbsorbProp[10]= 0;
  return(0);
} /* _mirror_full_center_setpos */

/* component ArmForNeutronPropState_2=Arm() SETTING, POSITION/ROTATION */
int _ArmForNeutronPropState_2_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_ArmForNeutronPropState_2_setpos] component ArmForNeutronPropState_2=Arm() SETTING [Arm:0]");
  stracpy(_ArmForNeutronPropState_2->_name, "ArmForNeutronPropState_2", 16384);
  stracpy(_ArmForNeutronPropState_2->_type, "Arm", 16384);
  _ArmForNeutronPropState_2->_index=11;
  /* component ArmForNeutronPropState_2=Arm() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Origin->_rotation_absolute, _ArmForNeutronPropState_2->_rotation_absolute);
    rot_transpose(_mirror_full_center->_rotation_absolute, tr1);
    rot_mul(_ArmForNeutronPropState_2->_rotation_absolute, tr1, _ArmForNeutronPropState_2->_rotation_relative);
    _ArmForNeutronPropState_2->_rotation_is_identity =  rot_test_identity(_ArmForNeutronPropState_2->_rotation_relative);
    tc1 = coords_set(
      0, 0, ArmMidOnePos);
    rot_transpose(_Origin->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _ArmForNeutronPropState_2->_position_absolute = coords_add(_Origin->_position_absolute, tc2);
    tc1 = coords_sub(_mirror_full_center->_position_absolute, _ArmForNeutronPropState_2->_position_absolute);
    _ArmForNeutronPropState_2->_position_relative = rot_apply(_ArmForNeutronPropState_2->_rotation_absolute, tc1);
  } /* ArmForNeutronPropState_2=Arm() AT ROTATED */
  DEBUG_COMPONENT("ArmForNeutronPropState_2", _ArmForNeutronPropState_2->_position_absolute, _ArmForNeutronPropState_2->_rotation_absolute);
  instrument->_position_absolute[11] = _ArmForNeutronPropState_2->_position_absolute;
  instrument->_position_relative[11] = _ArmForNeutronPropState_2->_position_relative;
  instrument->counter_N[11]  = instrument->counter_P[11] = instrument->counter_P2[11] = 0;
  instrument->counter_AbsorbProp[11]= 0;
  return(0);
} /* _ArmForNeutronPropState_2_setpos */

/* component guide_right=Mirror_Elliptic_Bispectral() SETTING, POSITION/ROTATION */
int _guide_right_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_guide_right_setpos] component guide_right=Mirror_Elliptic_Bispectral() SETTING [/usr/share/mcstas/3.0-dev/contrib/Mirror_Elliptic_Bispectral.comp:67]");
  stracpy(_guide_right->_name, "guide_right", 16384);
  stracpy(_guide_right->_type, "Mirror_Elliptic_Bispectral", 16384);
  _guide_right->_index=12;
  _guide_right->_parameters.reflect[0]='\0';
  #define reflect (_guide_right->_parameters.reflect)
  _guide_right->_parameters.focus_start_w = instrument->_parameters._focus_start_w;
  #define focus_start_w (_guide_right->_parameters.focus_start_w)
  _guide_right->_parameters.focus_end_w = instrument->_parameters._focus_end_w;
  #define focus_end_w (_guide_right->_parameters.focus_end_w)
  _guide_right->_parameters.focus_start_h = instrument->_parameters._focus_start_h;
  #define focus_start_h (_guide_right->_parameters.focus_start_h)
  _guide_right->_parameters.focus_end_h = instrument->_parameters._focus_end_h;
  #define focus_end_h (_guide_right->_parameters.focus_end_h)
  _guide_right->_parameters.mirror_start = instrument->_parameters._guide_start;
  #define mirror_start (_guide_right->_parameters.mirror_start)
  _guide_right->_parameters.m = m [ 0 ];
  #define m (_guide_right->_parameters.m)
  _guide_right->_parameters.smallaxis_w = instrument->_parameters._smallaxis_w;
  #define smallaxis_w (_guide_right->_parameters.smallaxis_w)
  _guide_right->_parameters.smallaxis_h = instrument->_parameters._smallaxis_h;
  #define smallaxis_h (_guide_right->_parameters.smallaxis_h)
  _guide_right->_parameters.length = extraction_start_pos + mirror_full_length - instrument->_parameters._guide_start;
  #define length (_guide_right->_parameters.length)
  _guide_right->_parameters.transmit = 0;
  #define transmit (_guide_right->_parameters.transmit)
  _guide_right->_parameters.substrate_thickness = 0;
  #define substrate_thickness (_guide_right->_parameters.substrate_thickness)
  _guide_right->_parameters.coating_thickness = 0;
  #define coating_thickness (_guide_right->_parameters.coating_thickness)

  #define pTable (_guide_right->_parameters.pTable)

  #undef reflect
  #undef focus_start_w
  #undef focus_end_w
  #undef focus_start_h
  #undef focus_end_h
  #undef mirror_start
  #undef m
  #undef smallaxis_w
  #undef smallaxis_h
  #undef length
  #undef transmit
  #undef substrate_thickness
  #undef coating_thickness
  #undef pTable
  /* component guide_right=Mirror_Elliptic_Bispectral() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (-90)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _Origin->_rotation_absolute, _guide_right->_rotation_absolute);
    rot_transpose(_ArmForNeutronPropState_2->_rotation_absolute, tr1);
    rot_mul(_guide_right->_rotation_absolute, tr1, _guide_right->_rotation_relative);
    _guide_right->_rotation_is_identity =  rot_test_identity(_guide_right->_rotation_relative);
    tc1 = coords_set(
      0, 0, instrument->_parameters._guide_start);
    rot_transpose(_Origin->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _guide_right->_position_absolute = coords_add(_Origin->_position_absolute, tc2);
    tc1 = coords_sub(_ArmForNeutronPropState_2->_position_absolute, _guide_right->_position_absolute);
    _guide_right->_position_relative = rot_apply(_guide_right->_rotation_absolute, tc1);
  } /* guide_right=Mirror_Elliptic_Bispectral() AT ROTATED */
  DEBUG_COMPONENT("guide_right", _guide_right->_position_absolute, _guide_right->_rotation_absolute);
  instrument->_position_absolute[12] = _guide_right->_position_absolute;
  instrument->_position_relative[12] = _guide_right->_position_relative;
  instrument->counter_N[12]  = instrument->counter_P[12] = instrument->counter_P2[12] = 0;
  instrument->counter_AbsorbProp[12]= 0;
  return(0);
} /* _guide_right_setpos */

/* component ArmForNeutronPropState_4=Arm() SETTING, POSITION/ROTATION */
int _ArmForNeutronPropState_4_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_ArmForNeutronPropState_4_setpos] component ArmForNeutronPropState_4=Arm() SETTING [Arm:0]");
  stracpy(_ArmForNeutronPropState_4->_name, "ArmForNeutronPropState_4", 16384);
  stracpy(_ArmForNeutronPropState_4->_type, "Arm", 16384);
  _ArmForNeutronPropState_4->_index=13;
  /* component ArmForNeutronPropState_4=Arm() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Origin->_rotation_absolute, _ArmForNeutronPropState_4->_rotation_absolute);
    rot_transpose(_guide_right->_rotation_absolute, tr1);
    rot_mul(_ArmForNeutronPropState_4->_rotation_absolute, tr1, _ArmForNeutronPropState_4->_rotation_relative);
    _ArmForNeutronPropState_4->_rotation_is_identity =  rot_test_identity(_ArmForNeutronPropState_4->_rotation_relative);
    tc1 = coords_set(
      0, 0, ArmMidOnePos);
    rot_transpose(_Origin->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _ArmForNeutronPropState_4->_position_absolute = coords_add(_Origin->_position_absolute, tc2);
    tc1 = coords_sub(_guide_right->_position_absolute, _ArmForNeutronPropState_4->_position_absolute);
    _ArmForNeutronPropState_4->_position_relative = rot_apply(_ArmForNeutronPropState_4->_rotation_absolute, tc1);
  } /* ArmForNeutronPropState_4=Arm() AT ROTATED */
  DEBUG_COMPONENT("ArmForNeutronPropState_4", _ArmForNeutronPropState_4->_position_absolute, _ArmForNeutronPropState_4->_rotation_absolute);
  instrument->_position_absolute[13] = _ArmForNeutronPropState_4->_position_absolute;
  instrument->_position_relative[13] = _ArmForNeutronPropState_4->_position_relative;
  instrument->counter_N[13]  = instrument->counter_P[13] = instrument->counter_P2[13] = 0;
  instrument->counter_AbsorbProp[13]= 0;
  return(0);
} /* _ArmForNeutronPropState_4_setpos */

/* component guide_bottom=Mirror_Elliptic_Bispectral() SETTING, POSITION/ROTATION */
int _guide_bottom_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_guide_bottom_setpos] component guide_bottom=Mirror_Elliptic_Bispectral() SETTING [/usr/share/mcstas/3.0-dev/contrib/Mirror_Elliptic_Bispectral.comp:67]");
  stracpy(_guide_bottom->_name, "guide_bottom", 16384);
  stracpy(_guide_bottom->_type, "Mirror_Elliptic_Bispectral", 16384);
  _guide_bottom->_index=14;
  _guide_bottom->_parameters.reflect[0]='\0';
  #define reflect (_guide_bottom->_parameters.reflect)
  _guide_bottom->_parameters.focus_start_w = instrument->_parameters._focus_start_h;
  #define focus_start_w (_guide_bottom->_parameters.focus_start_w)
  _guide_bottom->_parameters.focus_end_w = instrument->_parameters._focus_end_h;
  #define focus_end_w (_guide_bottom->_parameters.focus_end_w)
  _guide_bottom->_parameters.focus_start_h = instrument->_parameters._focus_start_w;
  #define focus_start_h (_guide_bottom->_parameters.focus_start_h)
  _guide_bottom->_parameters.focus_end_h = instrument->_parameters._focus_end_w;
  #define focus_end_h (_guide_bottom->_parameters.focus_end_h)
  _guide_bottom->_parameters.mirror_start = instrument->_parameters._guide_start;
  #define mirror_start (_guide_bottom->_parameters.mirror_start)
  _guide_bottom->_parameters.m = m [ 0 ];
  #define m (_guide_bottom->_parameters.m)
  _guide_bottom->_parameters.smallaxis_w = instrument->_parameters._smallaxis_h;
  #define smallaxis_w (_guide_bottom->_parameters.smallaxis_w)
  _guide_bottom->_parameters.smallaxis_h = instrument->_parameters._smallaxis_w;
  #define smallaxis_h (_guide_bottom->_parameters.smallaxis_h)
  _guide_bottom->_parameters.length = extraction_start_pos + mirror_full_length - instrument->_parameters._guide_start;
  #define length (_guide_bottom->_parameters.length)
  _guide_bottom->_parameters.transmit = 0;
  #define transmit (_guide_bottom->_parameters.transmit)
  _guide_bottom->_parameters.substrate_thickness = 0;
  #define substrate_thickness (_guide_bottom->_parameters.substrate_thickness)
  _guide_bottom->_parameters.coating_thickness = 0;
  #define coating_thickness (_guide_bottom->_parameters.coating_thickness)

  #define pTable (_guide_bottom->_parameters.pTable)

  #undef reflect
  #undef focus_start_w
  #undef focus_end_w
  #undef focus_start_h
  #undef focus_end_h
  #undef mirror_start
  #undef m
  #undef smallaxis_w
  #undef smallaxis_h
  #undef length
  #undef transmit
  #undef substrate_thickness
  #undef coating_thickness
  #undef pTable
  /* component guide_bottom=Mirror_Elliptic_Bispectral() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (-90)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _ArmForGuideBottom->_rotation_absolute, _guide_bottom->_rotation_absolute);
    rot_transpose(_ArmForNeutronPropState_4->_rotation_absolute, tr1);
    rot_mul(_guide_bottom->_rotation_absolute, tr1, _guide_bottom->_rotation_relative);
    _guide_bottom->_rotation_is_identity =  rot_test_identity(_guide_bottom->_rotation_relative);
    tc1 = coords_set(
      0, 0, instrument->_parameters._guide_start);
    rot_transpose(_ArmForGuideBottom->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _guide_bottom->_position_absolute = coords_add(_ArmForGuideBottom->_position_absolute, tc2);
    tc1 = coords_sub(_ArmForNeutronPropState_4->_position_absolute, _guide_bottom->_position_absolute);
    _guide_bottom->_position_relative = rot_apply(_guide_bottom->_rotation_absolute, tc1);
  } /* guide_bottom=Mirror_Elliptic_Bispectral() AT ROTATED */
  DEBUG_COMPONENT("guide_bottom", _guide_bottom->_position_absolute, _guide_bottom->_rotation_absolute);
  instrument->_position_absolute[14] = _guide_bottom->_position_absolute;
  instrument->_position_relative[14] = _guide_bottom->_position_relative;
  instrument->counter_N[14]  = instrument->counter_P[14] = instrument->counter_P2[14] = 0;
  instrument->counter_AbsorbProp[14]= 0;
  return(0);
} /* _guide_bottom_setpos */

/* component ArmForNeutronPropState_5=Arm() SETTING, POSITION/ROTATION */
int _ArmForNeutronPropState_5_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_ArmForNeutronPropState_5_setpos] component ArmForNeutronPropState_5=Arm() SETTING [Arm:0]");
  stracpy(_ArmForNeutronPropState_5->_name, "ArmForNeutronPropState_5", 16384);
  stracpy(_ArmForNeutronPropState_5->_type, "Arm", 16384);
  _ArmForNeutronPropState_5->_index=15;
  /* component ArmForNeutronPropState_5=Arm() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Origin->_rotation_absolute, _ArmForNeutronPropState_5->_rotation_absolute);
    rot_transpose(_guide_bottom->_rotation_absolute, tr1);
    rot_mul(_ArmForNeutronPropState_5->_rotation_absolute, tr1, _ArmForNeutronPropState_5->_rotation_relative);
    _ArmForNeutronPropState_5->_rotation_is_identity =  rot_test_identity(_ArmForNeutronPropState_5->_rotation_relative);
    tc1 = coords_set(
      0, 0, ArmMidOnePos);
    rot_transpose(_Origin->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _ArmForNeutronPropState_5->_position_absolute = coords_add(_Origin->_position_absolute, tc2);
    tc1 = coords_sub(_guide_bottom->_position_absolute, _ArmForNeutronPropState_5->_position_absolute);
    _ArmForNeutronPropState_5->_position_relative = rot_apply(_ArmForNeutronPropState_5->_rotation_absolute, tc1);
  } /* ArmForNeutronPropState_5=Arm() AT ROTATED */
  DEBUG_COMPONENT("ArmForNeutronPropState_5", _ArmForNeutronPropState_5->_position_absolute, _ArmForNeutronPropState_5->_rotation_absolute);
  instrument->_position_absolute[15] = _ArmForNeutronPropState_5->_position_absolute;
  instrument->_position_relative[15] = _ArmForNeutronPropState_5->_position_relative;
  instrument->counter_N[15]  = instrument->counter_P[15] = instrument->counter_P2[15] = 0;
  instrument->counter_AbsorbProp[15]= 0;
  return(0);
} /* _ArmForNeutronPropState_5_setpos */

/* component cold_lambda_guidestart=L_monitor() SETTING, POSITION/ROTATION */
int _cold_lambda_guidestart_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_cold_lambda_guidestart_setpos] component cold_lambda_guidestart=L_monitor() SETTING [/usr/share/mcstas/3.0-dev/monitors/L_monitor.comp:65]");
  stracpy(_cold_lambda_guidestart->_name, "cold_lambda_guidestart", 16384);
  stracpy(_cold_lambda_guidestart->_type, "L_monitor", 16384);
  _cold_lambda_guidestart->_index=16;
  _cold_lambda_guidestart->_parameters.nL = 100;
  #define nL (_cold_lambda_guidestart->_parameters.nL)
  if("cold_lambda_guidestart" && strlen("cold_lambda_guidestart"))
    stracpy(_cold_lambda_guidestart->_parameters.filename, "cold_lambda_guidestart" ? "cold_lambda_guidestart" : "", 16384);
  else 
  _cold_lambda_guidestart->_parameters.filename[0]='\0';
  #define filename (_cold_lambda_guidestart->_parameters.filename)
  _cold_lambda_guidestart->_parameters.xmin = -0.05;
  #define xmin (_cold_lambda_guidestart->_parameters.xmin)
  _cold_lambda_guidestart->_parameters.xmax = 0.05;
  #define xmax (_cold_lambda_guidestart->_parameters.xmax)
  _cold_lambda_guidestart->_parameters.ymin = -0.05;
  #define ymin (_cold_lambda_guidestart->_parameters.ymin)
  _cold_lambda_guidestart->_parameters.ymax = 0.05;
  #define ymax (_cold_lambda_guidestart->_parameters.ymax)
  _cold_lambda_guidestart->_parameters.xwidth = w [ 0 ];
  #define xwidth (_cold_lambda_guidestart->_parameters.xwidth)
  _cold_lambda_guidestart->_parameters.yheight = h [ 0 ];
  #define yheight (_cold_lambda_guidestart->_parameters.yheight)
  _cold_lambda_guidestart->_parameters.Lmin = 0.01;
  #define Lmin (_cold_lambda_guidestart->_parameters.Lmin)
  _cold_lambda_guidestart->_parameters.Lmax = 20;
  #define Lmax (_cold_lambda_guidestart->_parameters.Lmax)
  _cold_lambda_guidestart->_parameters.restore_neutron = 1;
  #define restore_neutron (_cold_lambda_guidestart->_parameters.restore_neutron)

  #define L_N (_cold_lambda_guidestart->_parameters.L_N)
  #define L_p (_cold_lambda_guidestart->_parameters.L_p)
  #define L_p2 (_cold_lambda_guidestart->_parameters.L_p2)

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
  /* component cold_lambda_guidestart=L_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Origin->_rotation_absolute, _cold_lambda_guidestart->_rotation_absolute);
    rot_transpose(_ArmForNeutronPropState_5->_rotation_absolute, tr1);
    rot_mul(_cold_lambda_guidestart->_rotation_absolute, tr1, _cold_lambda_guidestart->_rotation_relative);
    _cold_lambda_guidestart->_rotation_is_identity =  rot_test_identity(_cold_lambda_guidestart->_rotation_relative);
    tc1 = coords_set(
      0, 0, guide_z_pos [ 0 ]);
    rot_transpose(_Origin->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _cold_lambda_guidestart->_position_absolute = coords_add(_Origin->_position_absolute, tc2);
    tc1 = coords_sub(_ArmForNeutronPropState_5->_position_absolute, _cold_lambda_guidestart->_position_absolute);
    _cold_lambda_guidestart->_position_relative = rot_apply(_cold_lambda_guidestart->_rotation_absolute, tc1);
  } /* cold_lambda_guidestart=L_monitor() AT ROTATED */
  DEBUG_COMPONENT("cold_lambda_guidestart", _cold_lambda_guidestart->_position_absolute, _cold_lambda_guidestart->_rotation_absolute);
  instrument->_position_absolute[16] = _cold_lambda_guidestart->_position_absolute;
  instrument->_position_relative[16] = _cold_lambda_guidestart->_position_relative;
  instrument->counter_N[16]  = instrument->counter_P[16] = instrument->counter_P2[16] = 0;
  instrument->counter_AbsorbProp[16]= 0;
  return(0);
} /* _cold_lambda_guidestart_setpos */

/* component thermal_lambda_guidestart=L_monitor() SETTING, POSITION/ROTATION */
int _thermal_lambda_guidestart_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_thermal_lambda_guidestart_setpos] component thermal_lambda_guidestart=L_monitor() SETTING [/usr/share/mcstas/3.0-dev/monitors/L_monitor.comp:65]");
  stracpy(_thermal_lambda_guidestart->_name, "thermal_lambda_guidestart", 16384);
  stracpy(_thermal_lambda_guidestart->_type, "L_monitor", 16384);
  _thermal_lambda_guidestart->_index=17;
  _thermal_lambda_guidestart->_parameters.nL = 100;
  #define nL (_thermal_lambda_guidestart->_parameters.nL)
  if("thermal_lambda_guidestart" && strlen("thermal_lambda_guidestart"))
    stracpy(_thermal_lambda_guidestart->_parameters.filename, "thermal_lambda_guidestart" ? "thermal_lambda_guidestart" : "", 16384);
  else 
  _thermal_lambda_guidestart->_parameters.filename[0]='\0';
  #define filename (_thermal_lambda_guidestart->_parameters.filename)
  _thermal_lambda_guidestart->_parameters.xmin = -0.05;
  #define xmin (_thermal_lambda_guidestart->_parameters.xmin)
  _thermal_lambda_guidestart->_parameters.xmax = 0.05;
  #define xmax (_thermal_lambda_guidestart->_parameters.xmax)
  _thermal_lambda_guidestart->_parameters.ymin = -0.05;
  #define ymin (_thermal_lambda_guidestart->_parameters.ymin)
  _thermal_lambda_guidestart->_parameters.ymax = 0.05;
  #define ymax (_thermal_lambda_guidestart->_parameters.ymax)
  _thermal_lambda_guidestart->_parameters.xwidth = w [ 0 ];
  #define xwidth (_thermal_lambda_guidestart->_parameters.xwidth)
  _thermal_lambda_guidestart->_parameters.yheight = h [ 0 ];
  #define yheight (_thermal_lambda_guidestart->_parameters.yheight)
  _thermal_lambda_guidestart->_parameters.Lmin = 0.01;
  #define Lmin (_thermal_lambda_guidestart->_parameters.Lmin)
  _thermal_lambda_guidestart->_parameters.Lmax = 20;
  #define Lmax (_thermal_lambda_guidestart->_parameters.Lmax)
  _thermal_lambda_guidestart->_parameters.restore_neutron = 1;
  #define restore_neutron (_thermal_lambda_guidestart->_parameters.restore_neutron)

  #define L_N (_thermal_lambda_guidestart->_parameters.L_N)
  #define L_p (_thermal_lambda_guidestart->_parameters.L_p)
  #define L_p2 (_thermal_lambda_guidestart->_parameters.L_p2)

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
  /* component thermal_lambda_guidestart=L_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Origin->_rotation_absolute, _thermal_lambda_guidestart->_rotation_absolute);
    rot_transpose(_cold_lambda_guidestart->_rotation_absolute, tr1);
    rot_mul(_thermal_lambda_guidestart->_rotation_absolute, tr1, _thermal_lambda_guidestart->_rotation_relative);
    _thermal_lambda_guidestart->_rotation_is_identity =  rot_test_identity(_thermal_lambda_guidestart->_rotation_relative);
    tc1 = coords_set(
      0, 0, guide_z_pos [ 0 ]);
    rot_transpose(_Origin->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _thermal_lambda_guidestart->_position_absolute = coords_add(_Origin->_position_absolute, tc2);
    tc1 = coords_sub(_cold_lambda_guidestart->_position_absolute, _thermal_lambda_guidestart->_position_absolute);
    _thermal_lambda_guidestart->_position_relative = rot_apply(_thermal_lambda_guidestart->_rotation_absolute, tc1);
  } /* thermal_lambda_guidestart=L_monitor() AT ROTATED */
  DEBUG_COMPONENT("thermal_lambda_guidestart", _thermal_lambda_guidestart->_position_absolute, _thermal_lambda_guidestart->_rotation_absolute);
  instrument->_position_absolute[17] = _thermal_lambda_guidestart->_position_absolute;
  instrument->_position_relative[17] = _thermal_lambda_guidestart->_position_relative;
  instrument->counter_N[17]  = instrument->counter_P[17] = instrument->counter_P2[17] = 0;
  instrument->counter_AbsorbProp[17]= 0;
  return(0);
} /* _thermal_lambda_guidestart_setpos */

/* component lambda_guidestart=L_monitor() SETTING, POSITION/ROTATION */
int _lambda_guidestart_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_lambda_guidestart_setpos] component lambda_guidestart=L_monitor() SETTING [/usr/share/mcstas/3.0-dev/monitors/L_monitor.comp:65]");
  stracpy(_lambda_guidestart->_name, "lambda_guidestart", 16384);
  stracpy(_lambda_guidestart->_type, "L_monitor", 16384);
  _lambda_guidestart->_index=18;
  _lambda_guidestart->_parameters.nL = 100;
  #define nL (_lambda_guidestart->_parameters.nL)
  if("lambda_guidestart" && strlen("lambda_guidestart"))
    stracpy(_lambda_guidestart->_parameters.filename, "lambda_guidestart" ? "lambda_guidestart" : "", 16384);
  else 
  _lambda_guidestart->_parameters.filename[0]='\0';
  #define filename (_lambda_guidestart->_parameters.filename)
  _lambda_guidestart->_parameters.xmin = -0.05;
  #define xmin (_lambda_guidestart->_parameters.xmin)
  _lambda_guidestart->_parameters.xmax = 0.05;
  #define xmax (_lambda_guidestart->_parameters.xmax)
  _lambda_guidestart->_parameters.ymin = -0.05;
  #define ymin (_lambda_guidestart->_parameters.ymin)
  _lambda_guidestart->_parameters.ymax = 0.05;
  #define ymax (_lambda_guidestart->_parameters.ymax)
  _lambda_guidestart->_parameters.xwidth = w [ 0 ];
  #define xwidth (_lambda_guidestart->_parameters.xwidth)
  _lambda_guidestart->_parameters.yheight = h [ 0 ];
  #define yheight (_lambda_guidestart->_parameters.yheight)
  _lambda_guidestart->_parameters.Lmin = 0.01;
  #define Lmin (_lambda_guidestart->_parameters.Lmin)
  _lambda_guidestart->_parameters.Lmax = 20;
  #define Lmax (_lambda_guidestart->_parameters.Lmax)
  _lambda_guidestart->_parameters.restore_neutron = 1;
  #define restore_neutron (_lambda_guidestart->_parameters.restore_neutron)

  #define L_N (_lambda_guidestart->_parameters.L_N)
  #define L_p (_lambda_guidestart->_parameters.L_p)
  #define L_p2 (_lambda_guidestart->_parameters.L_p2)

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
  /* component lambda_guidestart=L_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Origin->_rotation_absolute, _lambda_guidestart->_rotation_absolute);
    rot_transpose(_thermal_lambda_guidestart->_rotation_absolute, tr1);
    rot_mul(_lambda_guidestart->_rotation_absolute, tr1, _lambda_guidestart->_rotation_relative);
    _lambda_guidestart->_rotation_is_identity =  rot_test_identity(_lambda_guidestart->_rotation_relative);
    tc1 = coords_set(
      0, 0, guide_z_pos [ 0 ]);
    rot_transpose(_Origin->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _lambda_guidestart->_position_absolute = coords_add(_Origin->_position_absolute, tc2);
    tc1 = coords_sub(_thermal_lambda_guidestart->_position_absolute, _lambda_guidestart->_position_absolute);
    _lambda_guidestart->_position_relative = rot_apply(_lambda_guidestart->_rotation_absolute, tc1);
  } /* lambda_guidestart=L_monitor() AT ROTATED */
  DEBUG_COMPONENT("lambda_guidestart", _lambda_guidestart->_position_absolute, _lambda_guidestart->_rotation_absolute);
  instrument->_position_absolute[18] = _lambda_guidestart->_position_absolute;
  instrument->_position_relative[18] = _lambda_guidestart->_position_relative;
  instrument->counter_N[18]  = instrument->counter_P[18] = instrument->counter_P2[18] = 0;
  instrument->counter_AbsorbProp[18]= 0;
  return(0);
} /* _lambda_guidestart_setpos */

/* component guide_top=Mirror_Elliptic_Bispectral() SETTING, POSITION/ROTATION */
int _guide_top_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_guide_top_setpos] component guide_top=Mirror_Elliptic_Bispectral() SETTING [/usr/share/mcstas/3.0-dev/contrib/Mirror_Elliptic_Bispectral.comp:67]");
  stracpy(_guide_top->_name, "guide_top", 16384);
  stracpy(_guide_top->_type, "Mirror_Elliptic_Bispectral", 16384);
  _guide_top->_index=19;
  _guide_top->_parameters.reflect[0]='\0';
  #define reflect (_guide_top->_parameters.reflect)
  _guide_top->_parameters.focus_start_w = instrument->_parameters._focus_start_h;
  #define focus_start_w (_guide_top->_parameters.focus_start_w)
  _guide_top->_parameters.focus_end_w = instrument->_parameters._focus_end_h;
  #define focus_end_w (_guide_top->_parameters.focus_end_w)
  _guide_top->_parameters.focus_start_h = instrument->_parameters._focus_start_w;
  #define focus_start_h (_guide_top->_parameters.focus_start_h)
  _guide_top->_parameters.focus_end_h = instrument->_parameters._focus_end_w;
  #define focus_end_h (_guide_top->_parameters.focus_end_h)
  _guide_top->_parameters.mirror_start = instrument->_parameters._guide_start;
  #define mirror_start (_guide_top->_parameters.mirror_start)
  _guide_top->_parameters.m = m [ 0 ];
  #define m (_guide_top->_parameters.m)
  _guide_top->_parameters.smallaxis_w = instrument->_parameters._smallaxis_h;
  #define smallaxis_w (_guide_top->_parameters.smallaxis_w)
  _guide_top->_parameters.smallaxis_h = instrument->_parameters._smallaxis_w;
  #define smallaxis_h (_guide_top->_parameters.smallaxis_h)
  _guide_top->_parameters.length = extraction_start_pos + mirror_full_length - instrument->_parameters._guide_start;
  #define length (_guide_top->_parameters.length)
  _guide_top->_parameters.transmit = 0;
  #define transmit (_guide_top->_parameters.transmit)
  _guide_top->_parameters.substrate_thickness = 0;
  #define substrate_thickness (_guide_top->_parameters.substrate_thickness)
  _guide_top->_parameters.coating_thickness = 0;
  #define coating_thickness (_guide_top->_parameters.coating_thickness)

  #define pTable (_guide_top->_parameters.pTable)

  #undef reflect
  #undef focus_start_w
  #undef focus_end_w
  #undef focus_start_h
  #undef focus_end_h
  #undef mirror_start
  #undef m
  #undef smallaxis_w
  #undef smallaxis_h
  #undef length
  #undef transmit
  #undef substrate_thickness
  #undef coating_thickness
  #undef pTable
  /* component guide_top=Mirror_Elliptic_Bispectral() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (-90)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _ArmForGuideTop->_rotation_absolute, _guide_top->_rotation_absolute);
    rot_transpose(_lambda_guidestart->_rotation_absolute, tr1);
    rot_mul(_guide_top->_rotation_absolute, tr1, _guide_top->_rotation_relative);
    _guide_top->_rotation_is_identity =  rot_test_identity(_guide_top->_rotation_relative);
    tc1 = coords_set(
      0, 0, instrument->_parameters._guide_start);
    rot_transpose(_ArmForGuideTop->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _guide_top->_position_absolute = coords_add(_ArmForGuideTop->_position_absolute, tc2);
    tc1 = coords_sub(_lambda_guidestart->_position_absolute, _guide_top->_position_absolute);
    _guide_top->_position_relative = rot_apply(_guide_top->_rotation_absolute, tc1);
  } /* guide_top=Mirror_Elliptic_Bispectral() AT ROTATED */
  DEBUG_COMPONENT("guide_top", _guide_top->_position_absolute, _guide_top->_rotation_absolute);
  instrument->_position_absolute[19] = _guide_top->_position_absolute;
  instrument->_position_relative[19] = _guide_top->_position_relative;
  instrument->counter_N[19]  = instrument->counter_P[19] = instrument->counter_P2[19] = 0;
  instrument->counter_AbsorbProp[19]= 0;
  return(0);
} /* _guide_top_setpos */

/* component ArmForNeutronPropState_6=Arm() SETTING, POSITION/ROTATION */
int _ArmForNeutronPropState_6_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_ArmForNeutronPropState_6_setpos] component ArmForNeutronPropState_6=Arm() SETTING [Arm:0]");
  stracpy(_ArmForNeutronPropState_6->_name, "ArmForNeutronPropState_6", 16384);
  stracpy(_ArmForNeutronPropState_6->_type, "Arm", 16384);
  _ArmForNeutronPropState_6->_index=20;
  /* component ArmForNeutronPropState_6=Arm() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Origin->_rotation_absolute, _ArmForNeutronPropState_6->_rotation_absolute);
    rot_transpose(_guide_top->_rotation_absolute, tr1);
    rot_mul(_ArmForNeutronPropState_6->_rotation_absolute, tr1, _ArmForNeutronPropState_6->_rotation_relative);
    _ArmForNeutronPropState_6->_rotation_is_identity =  rot_test_identity(_ArmForNeutronPropState_6->_rotation_relative);
    tc1 = coords_set(
      0, 0, ArmMidOnePos);
    rot_transpose(_Origin->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _ArmForNeutronPropState_6->_position_absolute = coords_add(_Origin->_position_absolute, tc2);
    tc1 = coords_sub(_guide_top->_position_absolute, _ArmForNeutronPropState_6->_position_absolute);
    _ArmForNeutronPropState_6->_position_relative = rot_apply(_ArmForNeutronPropState_6->_rotation_absolute, tc1);
  } /* ArmForNeutronPropState_6=Arm() AT ROTATED */
  DEBUG_COMPONENT("ArmForNeutronPropState_6", _ArmForNeutronPropState_6->_position_absolute, _ArmForNeutronPropState_6->_rotation_absolute);
  instrument->_position_absolute[20] = _ArmForNeutronPropState_6->_position_absolute;
  instrument->_position_relative[20] = _ArmForNeutronPropState_6->_position_relative;
  instrument->counter_N[20]  = instrument->counter_P[20] = instrument->counter_P2[20] = 0;
  instrument->counter_AbsorbProp[20]= 0;
  return(0);
} /* _ArmForNeutronPropState_6_setpos */

/* component guide_Left=Mirror_Elliptic_Bispectral() SETTING, POSITION/ROTATION */
int _guide_Left_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_guide_Left_setpos] component guide_Left=Mirror_Elliptic_Bispectral() SETTING [/usr/share/mcstas/3.0-dev/contrib/Mirror_Elliptic_Bispectral.comp:67]");
  stracpy(_guide_Left->_name, "guide_Left", 16384);
  stracpy(_guide_Left->_type, "Mirror_Elliptic_Bispectral", 16384);
  _guide_Left->_index=21;
  _guide_Left->_parameters.reflect[0]='\0';
  #define reflect (_guide_Left->_parameters.reflect)
  _guide_Left->_parameters.focus_start_w = instrument->_parameters._focus_start_w;
  #define focus_start_w (_guide_Left->_parameters.focus_start_w)
  _guide_Left->_parameters.focus_end_w = instrument->_parameters._focus_end_w;
  #define focus_end_w (_guide_Left->_parameters.focus_end_w)
  _guide_Left->_parameters.focus_start_h = instrument->_parameters._focus_start_h;
  #define focus_start_h (_guide_Left->_parameters.focus_start_h)
  _guide_Left->_parameters.focus_end_h = instrument->_parameters._focus_end_h;
  #define focus_end_h (_guide_Left->_parameters.focus_end_h)
  _guide_Left->_parameters.mirror_start = guide_left_start;
  #define mirror_start (_guide_Left->_parameters.mirror_start)
  _guide_Left->_parameters.m = m [ 0 ];
  #define m (_guide_Left->_parameters.m)
  _guide_Left->_parameters.smallaxis_w = instrument->_parameters._smallaxis_w;
  #define smallaxis_w (_guide_Left->_parameters.smallaxis_w)
  _guide_Left->_parameters.smallaxis_h = instrument->_parameters._smallaxis_h;
  #define smallaxis_h (_guide_Left->_parameters.smallaxis_h)
  _guide_Left->_parameters.length = extraction_start_pos + mirror_full_length - guide_left_start;
  #define length (_guide_Left->_parameters.length)
  _guide_Left->_parameters.transmit = 0;
  #define transmit (_guide_Left->_parameters.transmit)
  _guide_Left->_parameters.substrate_thickness = 0;
  #define substrate_thickness (_guide_Left->_parameters.substrate_thickness)
  _guide_Left->_parameters.coating_thickness = 0;
  #define coating_thickness (_guide_Left->_parameters.coating_thickness)

  #define pTable (_guide_Left->_parameters.pTable)

  #undef reflect
  #undef focus_start_w
  #undef focus_end_w
  #undef focus_start_h
  #undef focus_end_h
  #undef mirror_start
  #undef m
  #undef smallaxis_w
  #undef smallaxis_h
  #undef length
  #undef transmit
  #undef substrate_thickness
  #undef coating_thickness
  #undef pTable
  /* component guide_Left=Mirror_Elliptic_Bispectral() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (-90)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _ArmForGuideLeft->_rotation_absolute, _guide_Left->_rotation_absolute);
    rot_transpose(_ArmForNeutronPropState_6->_rotation_absolute, tr1);
    rot_mul(_guide_Left->_rotation_absolute, tr1, _guide_Left->_rotation_relative);
    _guide_Left->_rotation_is_identity =  rot_test_identity(_guide_Left->_rotation_relative);
    tc1 = coords_set(
      0, 0, guide_left_start);
    rot_transpose(_ArmForGuideLeft->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _guide_Left->_position_absolute = coords_add(_ArmForGuideLeft->_position_absolute, tc2);
    tc1 = coords_sub(_ArmForNeutronPropState_6->_position_absolute, _guide_Left->_position_absolute);
    _guide_Left->_position_relative = rot_apply(_guide_Left->_rotation_absolute, tc1);
  } /* guide_Left=Mirror_Elliptic_Bispectral() AT ROTATED */
  DEBUG_COMPONENT("guide_Left", _guide_Left->_position_absolute, _guide_Left->_rotation_absolute);
  instrument->_position_absolute[21] = _guide_Left->_position_absolute;
  instrument->_position_relative[21] = _guide_Left->_position_relative;
  instrument->counter_N[21]  = instrument->counter_P[21] = instrument->counter_P2[21] = 0;
  instrument->counter_AbsorbProp[21]= 0;
  return(0);
} /* _guide_Left_setpos */

/* component ArmForNeutronPropState_7=Arm() SETTING, POSITION/ROTATION */
int _ArmForNeutronPropState_7_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_ArmForNeutronPropState_7_setpos] component ArmForNeutronPropState_7=Arm() SETTING [Arm:0]");
  stracpy(_ArmForNeutronPropState_7->_name, "ArmForNeutronPropState_7", 16384);
  stracpy(_ArmForNeutronPropState_7->_type, "Arm", 16384);
  _ArmForNeutronPropState_7->_index=22;
  /* component ArmForNeutronPropState_7=Arm() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Origin->_rotation_absolute, _ArmForNeutronPropState_7->_rotation_absolute);
    rot_transpose(_guide_Left->_rotation_absolute, tr1);
    rot_mul(_ArmForNeutronPropState_7->_rotation_absolute, tr1, _ArmForNeutronPropState_7->_rotation_relative);
    _ArmForNeutronPropState_7->_rotation_is_identity =  rot_test_identity(_ArmForNeutronPropState_7->_rotation_relative);
    tc1 = coords_set(
      0, 0, ArmMidOnePos);
    rot_transpose(_Origin->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _ArmForNeutronPropState_7->_position_absolute = coords_add(_Origin->_position_absolute, tc2);
    tc1 = coords_sub(_guide_Left->_position_absolute, _ArmForNeutronPropState_7->_position_absolute);
    _ArmForNeutronPropState_7->_position_relative = rot_apply(_ArmForNeutronPropState_7->_rotation_absolute, tc1);
  } /* ArmForNeutronPropState_7=Arm() AT ROTATED */
  DEBUG_COMPONENT("ArmForNeutronPropState_7", _ArmForNeutronPropState_7->_position_absolute, _ArmForNeutronPropState_7->_rotation_absolute);
  instrument->_position_absolute[22] = _ArmForNeutronPropState_7->_position_absolute;
  instrument->_position_relative[22] = _ArmForNeutronPropState_7->_position_relative;
  instrument->counter_N[22]  = instrument->counter_P[22] = instrument->counter_P2[22] = 0;
  instrument->counter_AbsorbProp[22]= 0;
  return(0);
} /* _ArmForNeutronPropState_7_setpos */

/* component ArmMidTwo=Arm() SETTING, POSITION/ROTATION */
int _ArmMidTwo_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_ArmMidTwo_setpos] component ArmMidTwo=Arm() SETTING [Arm:0]");
  stracpy(_ArmMidTwo->_name, "ArmMidTwo", 16384);
  stracpy(_ArmMidTwo->_type, "Arm", 16384);
  _ArmMidTwo->_index=23;
  /* component ArmMidTwo=Arm() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _ArmMidOne->_rotation_absolute, _ArmMidTwo->_rotation_absolute);
    rot_transpose(_ArmForNeutronPropState_7->_rotation_absolute, tr1);
    rot_mul(_ArmMidTwo->_rotation_absolute, tr1, _ArmMidTwo->_rotation_relative);
    _ArmMidTwo->_rotation_is_identity =  rot_test_identity(_ArmMidTwo->_rotation_relative);
    tc1 = coords_set(
      0, 0, 0);
    rot_transpose(_ArmMidOne->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _ArmMidTwo->_position_absolute = coords_add(_ArmMidOne->_position_absolute, tc2);
    tc1 = coords_sub(_ArmForNeutronPropState_7->_position_absolute, _ArmMidTwo->_position_absolute);
    _ArmMidTwo->_position_relative = rot_apply(_ArmMidTwo->_rotation_absolute, tc1);
  } /* ArmMidTwo=Arm() AT ROTATED */
  DEBUG_COMPONENT("ArmMidTwo", _ArmMidTwo->_position_absolute, _ArmMidTwo->_rotation_absolute);
  instrument->_position_absolute[23] = _ArmMidTwo->_position_absolute;
  instrument->_position_relative[23] = _ArmMidTwo->_position_relative;
  instrument->counter_N[23]  = instrument->counter_P[23] = instrument->counter_P2[23] = 0;
  instrument->counter_AbsorbProp[23]= 0;
  return(0);
} /* _ArmMidTwo_setpos */

/* component ArmForNeutronPropState_8=Arm() SETTING, POSITION/ROTATION */
int _ArmForNeutronPropState_8_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_ArmForNeutronPropState_8_setpos] component ArmForNeutronPropState_8=Arm() SETTING [Arm:0]");
  stracpy(_ArmForNeutronPropState_8->_name, "ArmForNeutronPropState_8", 16384);
  stracpy(_ArmForNeutronPropState_8->_type, "Arm", 16384);
  _ArmForNeutronPropState_8->_index=24;
  /* component ArmForNeutronPropState_8=Arm() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Origin->_rotation_absolute, _ArmForNeutronPropState_8->_rotation_absolute);
    rot_transpose(_ArmMidTwo->_rotation_absolute, tr1);
    rot_mul(_ArmForNeutronPropState_8->_rotation_absolute, tr1, _ArmForNeutronPropState_8->_rotation_relative);
    _ArmForNeutronPropState_8->_rotation_is_identity =  rot_test_identity(_ArmForNeutronPropState_8->_rotation_relative);
    tc1 = coords_set(
      0, 0, ArmMidOnePos);
    rot_transpose(_Origin->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _ArmForNeutronPropState_8->_position_absolute = coords_add(_Origin->_position_absolute, tc2);
    tc1 = coords_sub(_ArmMidTwo->_position_absolute, _ArmForNeutronPropState_8->_position_absolute);
    _ArmForNeutronPropState_8->_position_relative = rot_apply(_ArmForNeutronPropState_8->_rotation_absolute, tc1);
  } /* ArmForNeutronPropState_8=Arm() AT ROTATED */
  DEBUG_COMPONENT("ArmForNeutronPropState_8", _ArmForNeutronPropState_8->_position_absolute, _ArmForNeutronPropState_8->_rotation_absolute);
  instrument->_position_absolute[24] = _ArmForNeutronPropState_8->_position_absolute;
  instrument->_position_relative[24] = _ArmForNeutronPropState_8->_position_relative;
  instrument->counter_N[24]  = instrument->counter_P[24] = instrument->counter_P2[24] = 0;
  instrument->counter_AbsorbProp[24]= 0;
  return(0);
} /* _ArmForNeutronPropState_8_setpos */

/* component ArmMidThree=Arm() SETTING, POSITION/ROTATION */
int _ArmMidThree_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_ArmMidThree_setpos] component ArmMidThree=Arm() SETTING [Arm:0]");
  stracpy(_ArmMidThree->_name, "ArmMidThree", 16384);
  stracpy(_ArmMidThree->_type, "Arm", 16384);
  _ArmMidThree->_index=25;
  /* component ArmMidThree=Arm() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _ArmMidOne->_rotation_absolute, _ArmMidThree->_rotation_absolute);
    rot_transpose(_ArmForNeutronPropState_8->_rotation_absolute, tr1);
    rot_mul(_ArmMidThree->_rotation_absolute, tr1, _ArmMidThree->_rotation_relative);
    _ArmMidThree->_rotation_is_identity =  rot_test_identity(_ArmMidThree->_rotation_relative);
    tc1 = coords_set(
      0, 0, 0);
    rot_transpose(_ArmMidOne->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _ArmMidThree->_position_absolute = coords_add(_ArmMidOne->_position_absolute, tc2);
    tc1 = coords_sub(_ArmForNeutronPropState_8->_position_absolute, _ArmMidThree->_position_absolute);
    _ArmMidThree->_position_relative = rot_apply(_ArmMidThree->_rotation_absolute, tc1);
  } /* ArmMidThree=Arm() AT ROTATED */
  DEBUG_COMPONENT("ArmMidThree", _ArmMidThree->_position_absolute, _ArmMidThree->_rotation_absolute);
  instrument->_position_absolute[25] = _ArmMidThree->_position_absolute;
  instrument->_position_relative[25] = _ArmMidThree->_position_relative;
  instrument->counter_N[25]  = instrument->counter_P[25] = instrument->counter_P2[25] = 0;
  instrument->counter_AbsorbProp[25]= 0;
  return(0);
} /* _ArmMidThree_setpos */

/* component ArmExit=Arm() SETTING, POSITION/ROTATION */
int _ArmExit_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_ArmExit_setpos] component ArmExit=Arm() SETTING [Arm:0]");
  stracpy(_ArmExit->_name, "ArmExit", 16384);
  stracpy(_ArmExit->_type, "Arm", 16384);
  _ArmExit->_index=26;
  /* component ArmExit=Arm() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Origin->_rotation_absolute, _ArmExit->_rotation_absolute);
    rot_transpose(_ArmMidThree->_rotation_absolute, tr1);
    rot_mul(_ArmExit->_rotation_absolute, tr1, _ArmExit->_rotation_relative);
    _ArmExit->_rotation_is_identity =  rot_test_identity(_ArmExit->_rotation_relative);
    tc1 = coords_set(
      0, 0, ArmExitPos);
    rot_transpose(_Origin->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _ArmExit->_position_absolute = coords_add(_Origin->_position_absolute, tc2);
    tc1 = coords_sub(_ArmMidThree->_position_absolute, _ArmExit->_position_absolute);
    _ArmExit->_position_relative = rot_apply(_ArmExit->_rotation_absolute, tc1);
  } /* ArmExit=Arm() AT ROTATED */
  DEBUG_COMPONENT("ArmExit", _ArmExit->_position_absolute, _ArmExit->_rotation_absolute);
  instrument->_position_absolute[26] = _ArmExit->_position_absolute;
  instrument->_position_relative[26] = _ArmExit->_position_relative;
  instrument->counter_N[26]  = instrument->counter_P[26] = instrument->counter_P2[26] = 0;
  instrument->counter_AbsorbProp[26]= 0;
  return(0);
} /* _ArmExit_setpos */

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

_class_ESS_moderator_long_2001 *class_ESS_moderator_long_2001_init(_class_ESS_moderator_long_2001 *_comp
) {
  #define size (_comp->_parameters.size)
  #define l_low (_comp->_parameters.l_low)
  #define l_high (_comp->_parameters.l_high)
  #define dist (_comp->_parameters.dist)
  #define xw (_comp->_parameters.xw)
  #define yh (_comp->_parameters.yh)
  #define freq (_comp->_parameters.freq)
  #define T (_comp->_parameters.T)
  #define tau (_comp->_parameters.tau)
  #define tau1 (_comp->_parameters.tau1)
  #define tau2 (_comp->_parameters.tau2)
  #define d (_comp->_parameters.d)
  #define n (_comp->_parameters.n)
  #define n2 (_comp->_parameters.n2)
  #define chi2 (_comp->_parameters.chi2)
  #define I0 (_comp->_parameters.I0)
  #define I2 (_comp->_parameters.I2)
  #define branch1 (_comp->_parameters.branch1)
  #define branch2 (_comp->_parameters.branch2)
  #define branch_tail (_comp->_parameters.branch_tail)
  #define twopulses (_comp->_parameters.twopulses)
  #define target_index (_comp->_parameters.target_index)
  #define l_range (_comp->_parameters.l_range)
  #define w_mult (_comp->_parameters.w_mult)
  #define branchframe (_comp->_parameters.branchframe)
  #define tx (_comp->_parameters.tx)
  #define ty (_comp->_parameters.ty)
  #define tz (_comp->_parameters.tz)
  if (n == 1 || n2 == 1 || l_low<=0 || l_high <=0 ||
      branch2 == 0 || branch_tail == 0 || tau == 0)
  {
    printf("ESS_moderator_long_2001: %s: Check parameters (lead to Math Error).\n Avoid 0 value for {l_low l_high d tau branch1/2/tail} and 1 value for {n n2 branch1/2/tail}\n", NAME_CURRENT_COMP);
    exit(-1);
  }

  if (tau1==0 && !(branch1==1)) {
    branch1=1;
    printf("ESS_moderator_long_2001: %s: WARNING: Setting tau1 to zero implies branch 1=1.\n", NAME_CURRENT_COMP);
  }

  if (twopulses) {
    branchframe = 0.5;
    printf("ESS_moderator_long_2001: %s: INFO: Running with TWO pulses\n", NAME_CURRENT_COMP);
  } else {
    branchframe = 0;
    printf("ESS_moderator_long_2001: %s: INFO: Running with ONE pulse\n", NAME_CURRENT_COMP);
  }

  if (target_index)  {
    Coords ToTarget;
    ToTarget = coords_sub(POS_A_COMP_INDEX(INDEX_CURRENT_COMP+target_index),POS_A_CURRENT_COMP);
    ToTarget = rot_apply(ROT_A_CURRENT_COMP, ToTarget);
    coords_get(ToTarget, &tx, &ty, &tz);
    if (dist) {
      printf("ESS_moderator_long_2001: %s: WARNING: You have specified BOTH dist and non-zero target_index parameters.\n   !! Overwriting dist by value defined by the target_index. !!\n", NAME_CURRENT_COMP);
      exit(-1);
    }
    dist = sqrt(tx*tx+ty*ty+tz*tz);
  } else {
    if (!dist) {
      printf("ESS_moderator_long_2001: %s: ERROR: Target undefined, neither dist or target_index were set!\n", NAME_CURRENT_COMP);
      exit(-1);
    }
    tx = 0;  ty = 0; tz = dist;
  }

  l_range = l_high-l_low;
  w_mult = size*size*1.0e4;     /* source area correction */
  w_mult *= l_range;            /* wavelength range correction */
  w_mult *= 1.0/mcget_ncount();   /* Correct for number of rays */
  w_mult *= 50.0/3.0;           /* Correct for baseline frequency setting */
  #undef size
  #undef l_low
  #undef l_high
  #undef dist
  #undef xw
  #undef yh
  #undef freq
  #undef T
  #undef tau
  #undef tau1
  #undef tau2
  #undef d
  #undef n
  #undef n2
  #undef chi2
  #undef I0
  #undef I2
  #undef branch1
  #undef branch2
  #undef branch_tail
  #undef twopulses
  #undef target_index
  #undef l_range
  #undef w_mult
  #undef branchframe
  #undef tx
  #undef ty
  #undef tz
  return(_comp);
} /* class_ESS_moderator_long_2001_init */

_class_Mirror_Curved_Bispectral *class_Mirror_Curved_Bispectral_init(_class_Mirror_Curved_Bispectral *_comp
) {
  #define reflect (_comp->_parameters.reflect)
  #define focus_s (_comp->_parameters.focus_s)
  #define focus_e (_comp->_parameters.focus_e)
  #define mirror_start (_comp->_parameters.mirror_start)
  #define guide_start (_comp->_parameters.guide_start)
  #define yheight (_comp->_parameters.yheight)
  #define smallaxis (_comp->_parameters.smallaxis)
  #define length (_comp->_parameters.length)
  #define m (_comp->_parameters.m)
  #define transmit (_comp->_parameters.transmit)
  #define substrate_thickness (_comp->_parameters.substrate_thickness)
  #define coating_thickness (_comp->_parameters.coating_thickness)
  #define theta_1 (_comp->_parameters.theta_1)
  #define theta_2 (_comp->_parameters.theta_2)
  #define theta_3 (_comp->_parameters.theta_3)
  #define pTable (_comp->_parameters.pTable)
  if (reflect && strlen(reflect) && strcmp(reflect,"NULL") && strcmp(reflect,"0")) {
    if (Table_Read(&pTable, reflect, 1) <= 0) /* read 1st block data from file into pTable */
      exit(fprintf(stderr,"Mirror: %s: can not read file %s\n", NAME_CURRENT_COMP, reflect));
  }
  #undef reflect
  #undef focus_s
  #undef focus_e
  #undef mirror_start
  #undef guide_start
  #undef yheight
  #undef smallaxis
  #undef length
  #undef m
  #undef transmit
  #undef substrate_thickness
  #undef coating_thickness
  #undef theta_1
  #undef theta_2
  #undef theta_3
  #undef pTable
  return(_comp);
} /* class_Mirror_Curved_Bispectral_init */

_class_Mirror_Elliptic_Bispectral *class_Mirror_Elliptic_Bispectral_init(_class_Mirror_Elliptic_Bispectral *_comp
) {
  #define reflect (_comp->_parameters.reflect)
  #define focus_start_w (_comp->_parameters.focus_start_w)
  #define focus_end_w (_comp->_parameters.focus_end_w)
  #define focus_start_h (_comp->_parameters.focus_start_h)
  #define focus_end_h (_comp->_parameters.focus_end_h)
  #define mirror_start (_comp->_parameters.mirror_start)
  #define m (_comp->_parameters.m)
  #define smallaxis_w (_comp->_parameters.smallaxis_w)
  #define smallaxis_h (_comp->_parameters.smallaxis_h)
  #define length (_comp->_parameters.length)
  #define transmit (_comp->_parameters.transmit)
  #define substrate_thickness (_comp->_parameters.substrate_thickness)
  #define coating_thickness (_comp->_parameters.coating_thickness)
  #define pTable (_comp->_parameters.pTable)
  if (reflect && strlen(reflect) && strcmp(reflect,"NULL") && strcmp(reflect,"0")) {
    if (Table_Read(&pTable, reflect, 1) <= 0) /* read 1st block data from file into pTable */
      exit(fprintf(stderr,"Mirror: %s: can not read file %s\n", NAME_CURRENT_COMP, reflect));
  }
  #undef reflect
  #undef focus_start_w
  #undef focus_end_w
  #undef focus_start_h
  #undef focus_end_h
  #undef mirror_start
  #undef m
  #undef smallaxis_w
  #undef smallaxis_h
  #undef length
  #undef transmit
  #undef substrate_thickness
  #undef coating_thickness
  #undef pTable
  return(_comp);
} /* class_Mirror_Elliptic_Bispectral_init */

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



int init(void) { /* called by mccode_main for ESS_2001_bispectral:INITIALISE */
  DEBUG_INSTR();

  /* code_main/parseoptions/readparams sets instrument parameters value */
  stracpy(instrument->_name, "ESS_2001_bispectral", 256);

  /* Instrument 'ESS_2001_bispectral' INITIALISE */
  SIG_MESSAGE("[ESS_2001_bispectral] INITIALISE [ESS_2001_bispectral.instr:241]");
  #define man_lam (instrument->_parameters._man_lam)
  #define thermal (instrument->_parameters._thermal)
  #define cold_Lam_min (instrument->_parameters._cold_Lam_min)
  #define cold_Lam_max (instrument->_parameters._cold_Lam_max)
  #define thermal_Lam_min (instrument->_parameters._thermal_Lam_min)
  #define thermal_Lam_max (instrument->_parameters._thermal_Lam_max)
  #define thermal_hotspot_fac (instrument->_parameters._thermal_hotspot_fac)
  #define cold_hotspot_fac (instrument->_parameters._cold_hotspot_fac)
  #define use_full_guide_flag (instrument->_parameters._use_full_guide_flag)
  #define guide_start (instrument->_parameters._guide_start)
  #define Length (instrument->_parameters._Length)
  #define focus_start_w (instrument->_parameters._focus_start_w)
  #define focus_end_w (instrument->_parameters._focus_end_w)
  #define smallaxis_w (instrument->_parameters._smallaxis_w)
  #define focus_start_h (instrument->_parameters._focus_start_h)
  #define focus_end_h (instrument->_parameters._focus_end_h)
  #define smallaxis_h (instrument->_parameters._smallaxis_h)
  #define maxdiv (instrument->_parameters._maxdiv)
  #define coldperthermal (instrument->_parameters._coldperthermal)
  #define mirror_type (instrument->_parameters._mirror_type)
  #define mirror_coating_type (instrument->_parameters._mirror_coating_type)
  #define mirror_offset (instrument->_parameters._mirror_offset)
  #define theta1 (instrument->_parameters._theta1)
  #define theta2 (instrument->_parameters._theta2)
  #define theta3 (instrument->_parameters._theta3)
  #define m_mirror (instrument->_parameters._m_mirror)
  #define h_mirror (instrument->_parameters._h_mirror)
  #define Pulse_width (instrument->_parameters._Pulse_width)
  #define frequency (instrument->_parameters._frequency)
  #define gravity (instrument->_parameters._gravity)
  #define substrate_thickness (instrument->_parameters._substrate_thickness)
  #define coating_thickness (instrument->_parameters._coating_thickness)
  #define m1 (instrument->_parameters._m1)
  #define m2 (instrument->_parameters._m2)
{

//input focus points are in real coordinates: focus_start_w=0 means at origin, focus_start_w=-1 means 1m behind source. In elliptic guide, they are in coordinates relative to the guide
focus_s_w=focus_start_w-guide_start;
focus_s_h=focus_start_h-guide_start;
focus_e_w=focus_end_w-guide_start;
focus_e_h=focus_end_h-guide_start;


focus_s_w=focus_start_w-4.5;
focus_s_h=focus_start_h-4.5;
focus_e_w=focus_end_w-4.5;
focus_e_h=focus_end_h-4.5;


guide_length=Length-0.5-4.5; //length of guide made of Guide_gravity_NewCoating2 components

if (focus_s_w>guide_start){
printf("-------warning-------- focus_s_w is inside guide"); }

if (focus_s_h>guide_start) {
printf("---------warning -------------- focus_s_h is inside guide");}


guide_dist=1e-6; //Distance between guide elements


//Define the coating of the guide
k=0;
for (k=0; k<n_elements+1; k++){
m[k]=m1;
if (k<(n_elements/10*3))
{m[k]=m2;}
if (k>(n_elements/10*7)-1)
{m[k]=m1;}
}



//calculate the guide shape
//	%include "earray_corrected.c"

// File earray_corrected.c inserted directly in instrument for inclusion in McStas 2.0

// BEGIN earray_corrected.c - K. Klenoe
// January 2010 Generic

// temp is the file which data is (optionally) written to
 FILE *data;
 data = fopen("temp", "w");

// Variables used only in this file
double expcoeff, elength_w, elength_h, coating_price_w, coating_price_h;
int i,j, n_elements_half, n_elements_half_plus_1;

n_elements_half=n_elements/2;
n_elements_half_plus_1=n_elements_half+1;


//printf("focus_s_w, focus_e_w, smallaxis_w %f, %f, %f",focus_s_w, focus_e_w, smallaxis_w);


// Next 4 loops calculate spacing of guide elements


expcoeff=log(guide_length/2+1)/n_elements_half;

i=0;
for (i=0; i<(n_elements_half_plus_1); i++)
distance[i]=exp(expcoeff*i)-1;

j=0;
i=n_elements_half;
for (i=n_elements_half; i<(n_elements+1); i++)
{distance[n_elements-j]=guide_length-exp(expcoeff*(j))+1; j++;}

i=0;
for (i=0; i<(n_elements); i++)
length[i]=distance[i+1]-distance[i];


// Next 4 loops calculate the shape of the guide. Modify for non-elliptical geometries.


double f_h=focus_e_h-focus_s_h;
double f_w=focus_e_w-focus_s_w;

elength_h=sqrt(f_h*f_h+smallaxis_h*smallaxis_h);
elength_w=sqrt(f_w*f_w+smallaxis_w*smallaxis_w);

i=0;
for (i=0; i<(n_elements_half_plus_1); i++){
h[i]=smallaxis_h*sqrt(1-(((distance[i]-focus_s_h)-f_h/2)/(elength_h/2))*(((distance[i]-focus_s_h)-f_h/2)/(elength_h/2)));
if (h[i]<0.005) h[i]=0.005;}

i=n_elements_half;
for (i=n_elements_half; i<(n_elements+1); i++){
h[i]=smallaxis_h*sqrt(1-(((distance[i]-focus_s_h)-f_h/2)/(elength_h/2))*(((distance[i]-focus_s_h)-f_h/2)/(elength_h/2)));
if (h[i]<0.005) h[i]=0.005;}

i=0;
for (i=0; i<(n_elements_half_plus_1); i++){
w[i]=smallaxis_w*sqrt(1-(((distance[i]-focus_s_w)-f_w/2)/(elength_w/2))*(((distance[i]-focus_s_w)-f_w/2)/(elength_w/2)));
if (w[i]<0.005) w[i]=0.005;}

i=n_elements_half;
for (i=n_elements_half; i<(n_elements+1); i++){
w[i]=smallaxis_w*sqrt(1-(((distance[i]-focus_s_w)-f_w/2)/(elength_w/2))*(((distance[i]-focus_s_w)-f_w/2)/(elength_w/2)));
if (w[i]<0.005) w[i]=0.005;}


// calculating alpha from m

i=0;
for (i=0; i<(n_elements); i++)
if (m[i]>5) {alpha[i]=3.5+1.02*(m[i]-5);}
else {alpha[i]=3.5;}


// Calculating dimensions if guide was extended to source and sample
/*
double w_at_source, h_at_source, w_at_sample, h_at_sample;

h_at_source=2*sqrt(smallaxis_h/10*(1-1/elength_h/elength_h*4*(elength_h/2+focus_s_h-(-guide_start))*(elength_h/2+focus_s_h-(-guide_start))));


w_at_source=2*sqrt(smallaxis_w/10*(1-1/elength_w/elength_w*4*(elength_w/2+focus_s_w-(-guide_start))*(elength_w/2+focus_s_w-distance_s[i])));


h_at_sample=2*sqrt(smallaxis_h/10*(1-1/elength_h/elength_h*4*(elength_h/2+focus_s_h-(guide_length+sample_dist))*(elength_h/2+focus_s_h-(guide_length+sample_dist))));


w_at_sample=2*sqrt(smallaxis_w/10*(1-1/elength_w/elength_w*4*(elength_w/2+focus_s_w-(guide_length+sample_dist))*(elength_w/2+focus_s_w-(guide_length+sample_dist))));
*/
//printf("\ndistance_s[24] = %f\n",distance_s[24]);
/*printf("h_at_source = %f\n",h_at_source);
printf("w_at_source = %f\n",w_at_source);
printf("h_at_sample = %f\n",h_at_sample);
printf("w_at_sample = %f\n",w_at_sample);*/

// For printing out the calculated values (if desired)


i=0;
for (i=0; i<(n_elements+1); i++)
fprintf(data,"distance[%i] %f\n",i, distance[i]);


i=0;
for (i=0; i<(n_elements); i++)
fprintf(data,"length[%i] %f\n",i, length[i]);

i=0;
for (i=0; i<(n_elements); i++)
fprintf(data,"h[%i] %f     w[%i] %f\n",i, h[i], i, w[i]);


//fprintf(data,"\n elength_h, expcoeff, * %f, %f",elength_h, expcoeff);
//}


// END earray_corrected.c

//	%include "cost.c"		 


// Everything with sources
//---------------------------------------------------------------------------	 
// Switch settings. thermal=0: cold, thermal=1: thermal, thermal=2: bispectral
if(thermal==0){flag=1; two_sources=0;   }
if(thermal==1){flag=0; two_sources=0;   }
if(thermal==2){flag=1; two_sources=1;   }

//if man_lam=0, set emitted wavelength interval to standard - not used very often

if(Length==50){
thermal_lam_min=0.2; 
thermal_lam_max=4.7;
cold_lam_min=2.75;
cold_lam_max=7.25;
}

if(Length==70){
thermal_lam_min=0.2; 
thermal_lam_max=4.7;
cold_lam_min=2.75;
cold_lam_max=7.25;
}

if(Length==150){
thermal_lam_min=0.75; 
thermal_lam_max=2.25;
cold_lam_min=4.25;
cold_lam_max=5.75;
}

thermal_lam_min=0.75; 
thermal_lam_max=2.25;
cold_lam_min=4.25;
cold_lam_max=5.75;

if(man_lam==0){cold_Lam_min=cold_lam_min; cold_Lam_max=cold_lam_max; thermal_Lam_min=thermal_lam_min; thermal_Lam_max=thermal_lam_max;}

//calculate the propability multiplier for each source: if more cold neutrons than thermal neutrons are simulated, their weight is less
coldmultiplier=(coldperthermal+1)/coldperthermal;
thermalmultiplier=(coldperthermal+1) ;

// for hot spot
thermal_hotspot_factor=thermal_hotspot_fac;
cold_hotspot_factor=cold_hotspot_fac;

//end of calculations for sources
//-------------------------------------------




//calculate parameters for guide to speed up TRACE
k=0;
for (k=0; k<n_elements+1; k++){
guide_piece_length[k]=length[k]-guide_dist;
}

//Calculate everything needed for the mirror

 //The positions where the angles for the mirror are given: Start of mirror, center of mirror, end of mirror;
x1=0; 
x2=mirror_full_length/2;  
x3=mirror_full_length;



//calculate  alpha_mirror
	if (m_mirror>5) {alpha_mirror=3.5+1.02*(m_mirror-5);}
	else {alpha_mirror=3.5;}



//end of calculations for mirror
//------------------------------------------------------


k=0;
for (k=0; k<n_elements+1; k++)
{guide_z_pos[k]=4.5+distance[k];
}



//Position the arms for proper propagation of neutrons
ArmMidOnePos=extraction_start_pos+mirror_full_length/2;
ArmExitPos=extraction_start_pos+mirror_full_length+1e-7;

// end of calculations for guide consisting of mirrors
//---------------------------------------------------------------------


//For focusing of the source

w_extractionstart=0.5*smallaxis_w*sqrt(1-   (f_w+2*focus_start_w-2*extraction_start_pos)*(f_w+2*focus_start_w-2*extraction_start_pos)/(elength_w*elength_w));
w_guide_leftstart=0.5*smallaxis_w*sqrt(1-   (f_w+2*focus_start_w-2*guide_left_start)    *(f_w+2*focus_start_w-2*guide_left_start)/(elength_w*elength_w));
h_extractionstart=0.5*smallaxis_h*sqrt(1-   (f_h+2*focus_start_h-2*extraction_start_pos)*(f_h+2*focus_start_h-2*extraction_start_pos)/(elength_h*elength_h));
h_guide_leftstart=0.5*smallaxis_h*sqrt(1-   (f_h+2*focus_start_h-2*guide_left_start    )*(f_h+2*focus_start_h-2*guide_left_start)/(elength_h*elength_h));



if (isnan(w_extractionstart)){
w_extractionstart=0;
} 

if (isnan(h_extractionstart)){
h_extractionstart=0;
} 

x_two=(guide_left_start-extraction_start_pos)/(guide_left_start)*(0.18-w_guide_leftstart/2)+w_guide_leftstart/2;

x_mid=0.5*(x_two-w_extractionstart);
x_focus=x_two+w_extractionstart;

//printf("--------w_extractionstart=%f, h_extractionstart=%f, w_guide_leftstart=%f,h_guide_leftstart=%f x_mid=%f, x_two=%f, x_focus=%f\n\n", w_extractionstart, h_extractionstart, w_guide_leftstart, h_guide_leftstart, x_mid, x_two, x_focus);
//printf("smallaxis_w=%f, f_w=%f, focus_start_w=%f, guide_start=%f, elength_w=%f, f_w+2*focus_start_w-2*guide_start=%f\n", smallaxis_w, f_w, focus_start_w, guide_start, elength_w, f_w+2*focus_start_w-2*guide_start);

if (guide_start>extraction_start_pos) {
y_focus=h_mirror;
}else{
y_focus=h_extractionstart; // h[0] is the yheight at 4.5 m: a bit larger focus because part of the guide is missing
}

printf("h[0]=%f, w[0]=%f\n",h[0],w[0]);

}
  #undef man_lam
  #undef thermal
  #undef cold_Lam_min
  #undef cold_Lam_max
  #undef thermal_Lam_min
  #undef thermal_Lam_max
  #undef thermal_hotspot_fac
  #undef cold_hotspot_fac
  #undef use_full_guide_flag
  #undef guide_start
  #undef Length
  #undef focus_start_w
  #undef focus_end_w
  #undef smallaxis_w
  #undef focus_start_h
  #undef focus_end_h
  #undef smallaxis_h
  #undef maxdiv
  #undef coldperthermal
  #undef mirror_type
  #undef mirror_coating_type
  #undef mirror_offset
  #undef theta1
  #undef theta2
  #undef theta3
  #undef m_mirror
  #undef h_mirror
  #undef Pulse_width
  #undef frequency
  #undef gravity
  #undef substrate_thickness
  #undef coating_thickness
  #undef m1
  #undef m2
  _Origin_setpos(); /* type Progress_bar */
  _ArmForGuideRight_setpos(); /* type Arm */
  _ArmForGuideBottom_setpos(); /* type Arm */
  _ArmForGuideTop_setpos(); /* type Arm */
  _ArmForGuideLeft_setpos(); /* type Arm */
  _cold_source_setpos(); /* type ESS_moderator_long_2001 */
  _thermal_source_setpos(); /* type ESS_moderator_long_2001 */
  _ColdFocus_setpos(); /* type Arm */
  _ArmMidOne_setpos(); /* type Arm */
  _mirror_full_center_setpos(); /* type Mirror_Curved_Bispectral */
  _ArmForNeutronPropState_2_setpos(); /* type Arm */
  _guide_right_setpos(); /* type Mirror_Elliptic_Bispectral */
  _ArmForNeutronPropState_4_setpos(); /* type Arm */
  _guide_bottom_setpos(); /* type Mirror_Elliptic_Bispectral */
  _ArmForNeutronPropState_5_setpos(); /* type Arm */
  _cold_lambda_guidestart_setpos(); /* type L_monitor */
  _thermal_lambda_guidestart_setpos(); /* type L_monitor */
  _lambda_guidestart_setpos(); /* type L_monitor */
  _guide_top_setpos(); /* type Mirror_Elliptic_Bispectral */
  _ArmForNeutronPropState_6_setpos(); /* type Arm */
  _guide_Left_setpos(); /* type Mirror_Elliptic_Bispectral */
  _ArmForNeutronPropState_7_setpos(); /* type Arm */
  _ArmMidTwo_setpos(); /* type Arm */
  _ArmForNeutronPropState_8_setpos(); /* type Arm */
  _ArmMidThree_setpos(); /* type Arm */
  _ArmExit_setpos(); /* type Arm */

  /* call iteratively all components INITIALISE */
  class_Progress_bar_init(_Origin);





  class_ESS_moderator_long_2001_init(_cold_source);

  class_ESS_moderator_long_2001_init(_thermal_source);



  class_Mirror_Curved_Bispectral_init(_mirror_full_center);


  class_Mirror_Elliptic_Bispectral_init(_guide_right);


  class_Mirror_Elliptic_Bispectral_init(_guide_bottom);


  class_L_monitor_init(_cold_lambda_guidestart);

  class_L_monitor_init(_thermal_lambda_guidestart);

  class_L_monitor_init(_lambda_guidestart);

  class_Mirror_Elliptic_Bispectral_init(_guide_top);


  class_Mirror_Elliptic_Bispectral_init(_guide_Left);






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

  // EXTEND code here
  if (!strcmp(_comp->_name, "Origin")) {
	if(two_sources!=0){                //switch between sources
		dummy=thermalmultiplier*rand01();
		if (dummy>1){
		flag=1;
		}
	else {flag=0;}
	}
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
_class_Arm *class_Arm_trace(_class_Arm *_comp
  , _class_particle *_particle) {
  ABSORBED=SCATTERED=RESTORE=0;


  // EXTEND code here
  if (!strcmp(_comp->_name, "ColdFocus")) {
PROP_Z0;
  }
  if (!strcmp(_comp->_name, "ArmMidOne")) {
 guide_scatt=0; 

//save the current state of the neutron. consider including spin
old_x_prop=x;
old_y_prop=y;
old_z_prop=z;

old_vx_prop=vx;
old_vy_prop=vy;
old_vz_prop=vz;

old_t_prop=t;
old_p_prop=p;

new_x_prop=x;
new_y_prop=y;
new_z_prop=z;

new_vx_prop=vx;
new_vy_prop=vy;
new_vz_prop=vz;

new_t_prop=1e15; //any large value to make sure t<new_t_prop the first time the neutron hits any of the components near the mirror
new_p_prop=p;
SCATTER;
  }
  if (!strcmp(_comp->_name, "ArmForNeutronPropState_2")) {
//save new parameters if the time it took to reach the component is less than the time of the previous component
if (guide_scatt==5){
if (t<new_t_prop){
new_x_prop=x;
new_y_prop=y;
new_z_prop=z;

new_vx_prop=vx;
new_vy_prop=vy;
new_vz_prop=vz;

new_t_prop=t;
new_p_prop=p;
}}
//reset neutron to where it was before mirror
x=old_x_prop;
y=old_y_prop;
z=old_z_prop;

vx=old_vx_prop;
vy=old_vy_prop;
vz=old_vz_prop;

t=old_t_prop;
p=old_p_prop;
SCATTER;
//printf("mirror1_scatt=%i\,z=%f, t=%f, old_t_prop=%f, new_t_prop=%f\n",guide_scatt,z,t,old_t_prop, new_t_prop);

  }
  if (!strcmp(_comp->_name, "ArmForNeutronPropState_4")) {

//save new parameters if the time it took to reach the component is less than the time of the previous component
//printf("3guide_scatt=%i\,z=%f, t=%f, old_t_prop=%f, new_t_prop=%f\n",guide_scatt,z,t,old_t_prop, new_t_prop);
if (guide_scatt==3){

if (t<new_t_prop){
new_x_prop=x;
new_y_prop=y;
new_z_prop=z;

new_vx_prop=vx;
new_vy_prop=vy;
new_vz_prop=vz;

new_t_prop=t;
new_p_prop=p;
}}
//reset neutron to where it was before mirror
x=old_x_prop;
y=old_y_prop;
z=old_z_prop;

vx=old_vx_prop;
vy=old_vy_prop;
vz=old_vz_prop;

t=old_t_prop;
p=old_p_prop;

SCATTER;

//printf("guide_scatt=%i\,z=%f, t=%f, old_t_prop=%f, new_t_prop=%f\n",guide_scatt,z,t,old_t_prop, new_t_prop);

  }
  if (!strcmp(_comp->_name, "ArmForNeutronPropState_5")) {

if (guide_scatt==2){

//save new parameters if the time it took to reach the component is less than the time of the previous component

if (t<new_t_prop){
new_x_prop=x;
new_y_prop=y;
new_z_prop=z;

new_vx_prop=vx;
new_vy_prop=vy;
new_vz_prop=vz;

new_t_prop=t;
new_p_prop=p;
}}
//reset neutron to where it was before mirror
x=old_x_prop;
y=old_y_prop;
z=old_z_prop;

vx=old_vx_prop;
vy=old_vy_prop;
vz=old_vz_prop;

t=old_t_prop;
p=old_p_prop;

SCATTER;

  }
  if (!strcmp(_comp->_name, "ArmForNeutronPropState_6")) {

if (guide_scatt==1){

//save new parameters if the time it took to reach the component is less than the time of the previous component

if (t<new_t_prop){
new_x_prop=x;
new_y_prop=y;
new_z_prop=z;

new_vx_prop=vx;
new_vy_prop=vy;
new_vz_prop=vz;

new_t_prop=t;
new_p_prop=p;
}}
//reset neutron to where it was before mirror
x=old_x_prop;
y=old_y_prop;
z=old_z_prop;

vx=old_vx_prop;
vy=old_vy_prop;
vz=old_vz_prop;

t=old_t_prop;
p=old_p_prop;

SCATTER;


  }
  if (!strcmp(_comp->_name, "ArmForNeutronPropState_7")) {

//save new parameters if the time it took to reach the component is less than the time of the previous component
if (guide_scatt==4){

if (t<new_t_prop){
new_x_prop=x;
new_y_prop=y;
new_z_prop=z;

new_vx_prop=vx;
new_vy_prop=vy;
new_vz_prop=vz;

new_t_prop=t;
new_p_prop=p;
}}
//reset neutron to where it was before mirror
x=old_x_prop;
y=old_y_prop;
z=old_z_prop;

vx=old_vx_prop;
vy=old_vy_prop;
vz=old_vz_prop;

t=old_t_prop;
p=old_p_prop;


SCATTER;


  }
  if (!strcmp(_comp->_name, "ArmMidTwo")) {
 if (guide_scatt==0)
 {SCATTER; }
  }
  if (!strcmp(_comp->_name, "ArmForNeutronPropState_8")) {
//printf("1guide_scatt=%i\,z=%f, t=%f, old_t_prop=%f, new_t_prop=%f\n",guide_scatt,z,t*1000,old_t_prop*1000, new_t_prop*1000);

if (guide_scatt>0)
{
//let the neutrons be scattered from the component with the lowest value of t
x=new_x_prop;
y=new_y_prop;
z=new_z_prop;

vx=new_vx_prop;
vy=new_vy_prop;
vz=new_vz_prop;

//printf("2guide_scatt=%i\,z=%f, t=%f, old_t_prop=%f, new_t_prop=%f\n\n",guide_scatt,z,t*1000,old_t_prop*1000, new_t_prop*1000);

t=new_t_prop;
if (t>9e14){
t=old_t_prop;}

p=new_p_prop;
//printf("2guide_scatt=%i\,z=%f, t=%f, old_t_prop=%f, new_t_prop=%f\n\n",guide_scatt,z,t*1000,old_t_prop*1000, new_t_prop*1000);

}

  }

  return(_comp);
} /* class_Arm_trace */

#pragma acc routine seq
_class_ESS_moderator_long_2001 *class_ESS_moderator_long_2001_trace(_class_ESS_moderator_long_2001 *_comp
  , _class_particle *_particle) {
  ABSORBED=SCATTERED=RESTORE=0;

  #define size (_comp->_parameters.size)
  #define l_low (_comp->_parameters.l_low)
  #define l_high (_comp->_parameters.l_high)
  #define dist (_comp->_parameters.dist)
  #define xw (_comp->_parameters.xw)
  #define yh (_comp->_parameters.yh)
  #define freq (_comp->_parameters.freq)
  #define T (_comp->_parameters.T)
  #define tau (_comp->_parameters.tau)
  #define tau1 (_comp->_parameters.tau1)
  #define tau2 (_comp->_parameters.tau2)
  #define d (_comp->_parameters.d)
  #define n (_comp->_parameters.n)
  #define n2 (_comp->_parameters.n2)
  #define chi2 (_comp->_parameters.chi2)
  #define I0 (_comp->_parameters.I0)
  #define I2 (_comp->_parameters.I2)
  #define branch1 (_comp->_parameters.branch1)
  #define branch2 (_comp->_parameters.branch2)
  #define branch_tail (_comp->_parameters.branch_tail)
  #define twopulses (_comp->_parameters.twopulses)
  #define target_index (_comp->_parameters.target_index)
  #define l_range (_comp->_parameters.l_range)
  #define w_mult (_comp->_parameters.w_mult)
  #define branchframe (_comp->_parameters.branchframe)
  #define tx (_comp->_parameters.tx)
  #define ty (_comp->_parameters.ty)
  #define tz (_comp->_parameters.tz)
  double v,tau_l,E,lambda,k,r,xf,yf,dx,dy,w_focus,tail_flag;

  z=0;

  x = 0.5*size*randpm1();
  y = 0.5*size*randpm1();         /* Choose initial position */

  randvec_target_rect_real(&xf, &yf, &r, &w_focus,
			   tx, ty, tz, xw, yh, ROT_A_CURRENT_COMP, x, y, z, 2);

  dx = xf-x;
  dy = yf-y;
  r = sqrt(dx*dx+dy*dy+dist*dist);

  lambda = l_low+l_range*rand01();    /* Choose from uniform distribution */
  k = 2*PI/lambda;
  v = K2V*k;

  vz = v*dist/r;
  vy = v*dy/r;
  vx = v*dx/r;


/*  printf("pos0 (%g %g %g), pos1 (%g %g %g), r: %g, v (%g %g %g), v %g\n",
  x,y,z,xf,yf,dist,r,vx,vy,vz, v);
  printf("l %g, w_focus %g \n", lambda, w_focus);  */

  tail_flag = (rand01()<branch_tail);   /* Choose tail/bulk */
 if (tail_flag)
 {
  if (rand01() < branch2)
  {
    if (tau1>0)
      if (rand01() < branch1)     /* Quick and dirty non-general solution */
      {  /* FIRST CASE a */
        tau_l = tau;
        p = 1/(branch1*branch2*branch_tail); /* Correct for switching prob. */
      }
      else
      {  /* FIRST CASE b */
        tau_l = tau1;
        p = 1/((1-branch1)*branch2*branch_tail); /* Correct for switching prob. */
      }
    else
      {
        tau_l = tau;
        p = 1/(branch2*branch_tail); /* Correct for switching prob. */
      }
    t = -tau_l*log(1e-12+rand01());       /* Sample from long-time tail a */
 /* Correct for true pulse shape */
    p *= w_focus;                         /* Correct for target focusing */
    p *= tau_l/d;                         /* Correct for tail part */
    p *= I0*w_mult*Mezei_M_fct(lambda,T);           /* Calculate true intensity */
  }
  else
  {
    /* SECOND CASE */
    tau_l = tau2*lambda;
    t = -tau_l*log(1e-12+rand01());       /* Sample from long-time tail */
    p = n2/(n2-1)*((1-exp(-d/tau_l))-(1-exp(-n2*d/tau_l))*exp(-(n2-1)*t/tau_l)/n);
                                          /* Correct for true pulse shape */
    p /= (1-branch2)*branch_tail;          /* Correct for switching prob. */
    p *= tau_l/d;                         /* Correct for tail part */
    p *= w_focus;                         /* Correct for target focusing */
    p *= I2*w_mult/(1+exp(chi2*lambda-2.2))/lambda;
                                          /* Calculate true intensity */
  }
  t += d;                                 /* Add pulse length */
 }
 else
 {
  t = d*rand01();                        /* Sample from bulk pulse */
  if (rand01() < branch2)
  {
    if (rand01() < branch1)     /* Quick and dirty non-general solution */
    {  /* FIRST CASE a */
      tau_l = tau;
      p = 1/(branch1*branch2*(1-branch_tail)); /* Correct for switching prob. */
    }
    else
    {  /* FIRST CASE b */
      tau_l = tau1;
      p = 1/((1-branch1)*branch2*(1-branch_tail)); /* Correct for switching prob. */
    }
    p *= 1-n/(n-1)*(exp(-t/tau_l)-exp(-n*t/tau_l)/n); /* Correct for true pulse shape */
    p *= w_focus;                         /* Correct for target focusing */
    p *= I0*w_mult*Mezei_M_fct(lambda,T);           /* Calculate true intensity */
  }
  else
  {
    /* SECOND CASE */
    tau_l = tau2*lambda;
    p = 1-n2/(n2-1)*(exp(-t/tau_l)-exp(-n2*t/tau_l)/n2); /* Correct for true pulse shape */
    p /= (1-branch2)*(1-branch_tail);      /* Correct for switching prob. */
    p *= w_focus;                         /* Correct for target focusing */
    p *= I2*w_mult/(1+exp(chi2*lambda-2.2))/lambda;
                                          /* Calculate true intensity */
  }
 }
 if (rand01()<branchframe){
   t+=1/freq;
 }

  // EXTEND code here
  if (!strcmp(_comp->_name, "cold_source")) {
if (flag==1){
/*
	if (cold_hotspot_factor>1){
		if((x-cold_hotspot_x_center)*(x-cold_hotspot_x_center)+(y-cold_hotspot_y_center)*(y-cold_hotspot_y_center) < cold_hotspot_dia/2.0*cold_hotspot_dia/2.0){
			p=p*cold_hotspot_factor;
		}
		else {
			p=p*(size*size-cold_hotspot_dia/2.0*cold_hotspot_dia/2.0*3.1416*cold_hotspot_factor)/(size*size-cold_hotspot_dia/2.0*cold_hotspot_dia/2.0*3.1416);}
	}
*/
	if(two_sources!=0){
		p=p*coldmultiplier;   //increase intensity because not all neutrons come from this source
	}
}
  }
  if (!strcmp(_comp->_name, "thermal_source")) {
if (flag==0){
/*
	if (thermal_hotspot_factor>1){
if((x-thermal_hotspot_x_center)*(x-thermal_hotspot_x_center)+(y-thermal_hotspot_y_center)*(y-thermal_hotspot_y_center) < thermal_hotspot_dia/2.0*thermal_hotspot_dia/2.0){p=p*thermal_hotspot_factor;}
else {
p=p*(size*size-thermal_hotspot_dia/2.0*thermal_hotspot_dia/2.0*3.1416*thermal_hotspot_factor)/(size*size-thermal_hotspot_dia/2.0*thermal_hotspot_dia/2.0*3.1416);}
}*/
xthermal=x; ythermal=y;  //for origin info
if(two_sources!=0){  //increase intensity because not all neutrons come from this source
p=p*thermalmultiplier;
}
}
  }

  #undef size
  #undef l_low
  #undef l_high
  #undef dist
  #undef xw
  #undef yh
  #undef freq
  #undef T
  #undef tau
  #undef tau1
  #undef tau2
  #undef d
  #undef n
  #undef n2
  #undef chi2
  #undef I0
  #undef I2
  #undef branch1
  #undef branch2
  #undef branch_tail
  #undef twopulses
  #undef target_index
  #undef l_range
  #undef w_mult
  #undef branchframe
  #undef tx
  #undef ty
  #undef tz
  return(_comp);
} /* class_ESS_moderator_long_2001_trace */

#pragma acc routine seq
_class_Mirror_Curved_Bispectral *class_Mirror_Curved_Bispectral_trace(_class_Mirror_Curved_Bispectral *_comp
  , _class_particle *_particle) {
  ABSORBED=SCATTERED=RESTORE=0;

  #define reflect (_comp->_parameters.reflect)
  #define focus_s (_comp->_parameters.focus_s)
  #define focus_e (_comp->_parameters.focus_e)
  #define mirror_start (_comp->_parameters.mirror_start)
  #define guide_start (_comp->_parameters.guide_start)
  #define yheight (_comp->_parameters.yheight)
  #define smallaxis (_comp->_parameters.smallaxis)
  #define length (_comp->_parameters.length)
  #define m (_comp->_parameters.m)
  #define transmit (_comp->_parameters.transmit)
  #define substrate_thickness (_comp->_parameters.substrate_thickness)
  #define coating_thickness (_comp->_parameters.coating_thickness)
  #define theta_1 (_comp->_parameters.theta_1)
  #define theta_2 (_comp->_parameters.theta_2)
  #define theta_3 (_comp->_parameters.theta_3)
  #define pTable (_comp->_parameters.pTable)
double f; //half of distance between focal points
double asquared;
double a; //half of ellipse length
double b; //half of ellipse width

double xprime; //x in coordinates with center of ellipse at (xprime,zprime)=(0,0)
double ymirror; //height of the mirror


//Defining the mirror
double a1;
double b1;
double c1;

//solving the time the neutron hits the sample
double A, B, C, D, E, P, Q, R, U, V, I, J, K;

//finding rotation of mirror
double alpha1, beta1, gamma1;
double theta_m;
double sin_theta_m, cos_theta_m;

double tan_theta_1;
double tan_theta_2;
double tan_theta_3;


double v_n; //speed of neutron perpendicular to surface

double Ref; //reflectivity

double dt;
double q;
  int intersect;

double discriminant;




double dt_2;
double dt_3;
int prop_case;
double x_2;
double y_2;
double z_2;
double t_2;
double x_3;
double y_3;
double z_3;
double t_3;

int x_hit;
int x_hit_2;
int x_hit_3;
double xprime_2;
double ymirror_2;
double xprime_3;
double ymirror_3;
int intersect_2;
int intersect_3;


intersect=0;
x_hit=0;
x_hit_2=0;
x_hit_3=0;
intersect_2=0;
intersect_3=0;
prop_case=0;

//printf("\n\n\n");
   double old_x = x, old_y = y, old_z = z, old_t=t, old_vx=vx, old_vz=vz, old_vy=vy;

// printf("x=%f, y=%f, z=%f, vx=%f, vy=%f, vz=%f\n",x,y,z,vx,vy,vz);

// Check if neutron hits mirror. First find which z,x coordinates it hits.

//mirror is defined by z(x)=a1x^3+b1x^2+c1x+d1, with dz/dx|x=-length/2=tan(theta_1), dz/dx|x=0=tan(theta_2), dz/dx|x=length/2=tan(theta3), z(0)=0. (d1=0)

tan_theta_1=tan(theta_1*DEG2RAD);
tan_theta_2=tan(theta_2*DEG2RAD);
tan_theta_3=tan(theta_3*DEG2RAD);


a1=2.0/3.0*(tan_theta_1+tan_theta_3-2.0*tan_theta_2)/(length*length);
b1=(tan_theta_3-tan_theta_1)/(2.0*length);
c1=tan_theta_2;


//neutron trajectory is defined by x=x0+vx*t, z=z0+vz*t. setting z=a1*x^3+b1*x^2+c1*x gives the equation A*t^3+B*t^2+C*t+D=0, with
A=a1*vx*vx*vx;
B=3.0*a1*x*vx*vx+b1*vx*vx;
C=3.0*a1*x*x*vx+2.0*b1*x*vx+c1*vx-vz;
D=a1*x*x*x+b1*x*x+c1*x-z;

//printf("a1=%f,b1=%f,c1=%f",a1,b1,c1);

//this equation must now be solved for t;

if (A!=0){
P=1/3.0*(3.0*C/A-B*B/(A*A));
Q=1/27.0*(2.0*B*B*B/(A*A*A)-9.0*B*C/(A*A)+27.0*D/A);

E=P*P*P/27.0+Q*Q/4.0;

// printf("A=%f, B=%f, C=%f, D=%f, 1e6P=%f, 1e6Q=%f, 1e6E=%f\n", A, B, C, D, 1e6*P, 1e6*Q, 1e6*E);

prop_case=0;
if (E>=0){

U=cbrt(-Q/2.0+sqrt(E));
V=cbrt(-Q/2.0-sqrt(E));

I=U+V-B/(3.0*A);
dt=I;
dt_2=I;
dt_3=I;
// printf("I=%f\n",I);

// J=-(U+V)/2+1i*(U-V)*sqrt(3)/2-B/(3*A) //complex solution
// K=-(U+V)/2-1i*(U-V)*sqrt(3)/2-B/(3*A) //complex solution
}else{
    R=acos(-Q/(2.0*sqrt(-P*P*P/27.0)));

// printf("R=%f\n",R);

   
   I=2.0*sqrt(fabs(P)/3.0)*cos(R/3.0)-B/A/3.0;
   J=-2.0*sqrt(fabs(P)/3.0)*cos(R/3.0 + 3.1415926535/3.0)-B/A/3.0;
   K=-2.0*sqrt(fabs(P)/3.0)*cos(R/3.0 - 3.1415926535/3.0)-B/A/3.0;

// printf("2.0*sqrt(abs(P)/3.0)=%f", 2.0*sqrt(abs(P)/3.0));
// printf("cos(R/3.0)=%f, cos(R/3.0 + 3.1415926535/3.0)=%f, cos(R/3.0 - 3.1415926535/3.0)=%f, -B/A/3.0=%f\n", cos(R/3.0), cos(R/3.0 + 3.1415926535/3.0), cos(R/3.0 - 3.1415926535/3.0), -B/A/3.0);

// printf("I=%f, J=%f, K=%f, \n",I, J, K);
// printf("P=%f, R=%f, A=%f, B=%f, \n",P, R, A, B);



// Three solutions. Find the smallest positive of these.
//there are problems with the solutions....
	if (I<=0){
		if (J<=0 && K<=0){dt=-1.0;} //if all three are negative, dt<0 and nothing happens
		if (J<=0 && K>0){dt=K;}  //if only K>0, dt=K
		if (J>0 && K<=0){dt=J;} //if only J>0, dt=J

		if (J>0 && K>0){	//if both J>0 and K>0, compare
			if (J>=K){dt=K; prop_case=1; dt_2=J;}else{dt=J; dt_2=I; prop_case=2;} } //dt is the smallest value
	}else{ //end if (I<=0)
		if (J<=0 && K<=0){dt=I;} //if only I>0, dt=I;

		if (J<=0 && K>0){ //if both I>0 and K>0, compare
			if (K>=I){dt=I; dt_2=K; prop_case=3;}else{dt=K; dt_2=I; prop_case=4;} } //dt is the smallest value

		if (J>0 && K<=0){ //if both I>0 and J>0, compare
			if (J>=I){dt=I; dt_2=J; prop_case=5;}else{dt=J; dt_2=I; prop_case=6;} } //dt is the smallest value

		if (J>0 && K>0){ //if all three>0, compare
			if (J>=K){ //either K or I is smallest
				if (K>=I){dt=I; if(J>=K){ dt_2=K; dt_3=J; prop_case=9;}else{dt_2=K; dt_3=J; prop_case=15;}}else{dt=K; if (J>=I){dt_2=I; dt_3=J; prop_case=10;}else{dt_2=J; dt_3=I; prop_case=11;} } //if K is smallest, compare it to I  
			}else{
				if (J>=I){dt=I; if (K>J){dt_2=J; dt_3=K; prop_case=12;}else{dt_2=J; dt_3=J; prop_case=16;}}else{dt=J; if (K>I){dt_2=I; dt_3=K; prop_case=13;}else{{dt_2=K; dt_3=I; prop_case=14;}}  }}  //else compare J to I
			} //end if(J>0 && K>0)
				
	} //end }else{ for if(I<=0)



}    // end }else{ for if (E>=0)
  

}else{ //end if (A!=0)  
if (B!=0){

discriminant=C*C-4*B*D;

if (discriminant<0){dt=-1.0;}else{ //only complex solutions: set dt<0 to avoid interaction
I=(-C-sqrt(discriminant))/(2.0*B);
J=(-C+sqrt(discriminant))/(2.0*B);

if (I<=0 && J<=0){dt = -1.0;} //both times are negative.
if (I<=0 && J>0 ){dt = J;} //set dt to only positive value.
if (I>0  && J<=0){dt = I;} //set dt to only positive value.
if (I>0  && J>0 ){if (I>J) {dt=J; dt_2=I; prop_case=7;}else{dt=I; dt_2=J; prop_case=8;} } //set dt to smallest positive value  

} //end if (discriminant<0){}else{
}else{ //end if (B!)=0
if (C!=0) { dt = -D/C;}else{
 printf("warning: A=B=C=0. Neutron is ignored\n"); }
} //end if(B!=0){}else{
} //end if (A!=0){}else{
//now intersection time has been found.

if (dt>0) { //if time is positive, propagate neutron to where it hits mirror. This is done without gravity.
// printf("before anything: x=%f,y=%f,z=%f,vx=%f,vy=%f,vz=%f, dt=%f\n",x,y,z,vx,vy,vz,dt);

    x += vx*dt;
    y += vy*dt;
    z += vz*dt;
    t += dt;


x_hit=(x >=-length/2 && x<=length/2);


if (prop_case==0){
x_2=x;
y_2=y;
z_2=z;
t_2=t;
x_3=x;
y_3=y;
z_3=z;
t_3=t;
}

if (prop_case>0)
{
x_2=old_x+vx*dt_2;
y_2=old_y+vy*dt_2;
z_2=old_z+vz*dt_2;
t_2=old_t+dt_2;
x_hit_2=(x_2 >=-length/2 && x_2<=length/2);
}

if (prop_case>8)
{
x_3=old_x+vx*dt_3;
y_3=old_y+vy*dt_3;
z_3=old_z+vz*dt_3;
t_3=old_t+dt_3;
x_hit_3=(x_3 >=-length/2 && x_3<=length/2);
}

//printf("x_hit=%d, x_hit_2=%d, x_hit_3=%d\n",x_hit, x_hit_2, x_hit_3);
//printf("dt=%f, dt_2=%f, dt_3=%f\n",dt,dt_2,dt_3);
// printf("x=%f,y=%f,z=%f,vx=%f,vy=%f,vz=%f\n",x,y,z,vx,vy,vz);

// printf("x=%f, length/2=%f\n",x, length/2);


if (x_hit || x_hit_2 || x_hit_3){
//if (x >=-length/2 && x<=length/2){ //check if neutron is within x limits of the mirror. If so, check if it is within y limits.


//define the ellipse
b=smallaxis/2;

f=(focus_e-focus_s)*0.5;

 asquared=f*f+b*b;
 a=sqrt(asquared);

xprime=-f-focus_s+mirror_start+length/2+x; //xprime is the x-coordinate in a coordinate system centered at the center of the ellipse

//ymirror=b*sqrt(1-xprime*xprime/(f*f)); //following Kaspars convention, assuming f~=a (valid for most elliptic guides normally used)

ymirror=b*sqrt(1-xprime*xprime/asquared);



xprime_2=-f-focus_s+mirror_start+length/2+x_2; //xprime is the x-coordinate in a coordinate system centered at the center of the ellipse
ymirror_2=b*sqrt(1-xprime_2*xprime_2/asquared);

xprime_3=-f-focus_s+mirror_start+length/2+x_3; //xprime is the x-coordinate in a coordinate system centered at the center of the ellipse
ymirror_3=b*sqrt(1-xprime_3*xprime_3/asquared);

if (guide_start>mirror_start){ //If (part of the) mirror is outside the guide, the mirror can be extended
if (  x<-length/2+guide_start-mirror_start) {
ymirror=yheight/2;
}

if (  x_2<-length/2+guide_start-mirror_start) {
ymirror_2=yheight/2;
}

if (  x_3<-length/2+guide_start-mirror_start) {
ymirror_3=yheight/2;
}




}











// printf("ymirror=%f, y=%f\n",ymirror, y);
intersect = ( y>=-ymirror && y<=ymirror && x >=-length/2 && x<=length/2);

if (prop_case>0) {
intersect_2 = ( y_2>=-ymirror && y_2<=ymirror && x_2 >=-length/2 && x_2<=length/2);
}
if (prop_case>8){
intersect_3 = ( y_3>=-ymirror && y_3<=ymirror && x_3 >=-length/2 && x_3<=length/2);
}

//printf("y_2=%f, ymirror=%f\n",y_2,ymirror);

//printf("\nintersect=%d, intersect_2=%d, intersect_3=%d, prop_case=%d\n",intersect, intersect_2, intersect_3, prop_case);

//printf("x=%f,y=%f,z=%f,t=%f\n",x,y,z,t);
//printf("x_2=%f,y_2=%f,z_2=%f,t_2=%f\n",x_2,y_2,z_2,t_2);
//printf("x_3=%f,y_3=%f,z_3=%f,t_3=%f\n",x_3,y_3,z_3,t_3);

if (!intersect){
if (!intersect_2){
intersect=intersect_3;
x=x_3;
y=y_3;
z=z_3;
t=t_3;
}else{
intersect=intersect_2;
x=x_2;
y=y_2;
z=z_2;
t=t_2;
}
}

//printf("intersect=%d, intersect_2=%d, intersect_3=%d, prop_case=%d\n\n",intersect, intersect_2, intersect_3, prop_case);
//printf("x=%f,y=%f,z=%f,t=%f\n",x,y,z,t);

//printf("z=%f, zcalc=%f\n",z,a1*x*x*x+b1*x*x+c1*x);
//printf("z=%f, zcalc=%f\n",z_2,a1*x_2*x_2*x_2+b1*x_2*x_2+c1*x_2);
//printf("z=%f, zcalc=%f\n",z_3,a1*x_3*x_3*x_3+b1*x_3*x_3+c1*x_3);

    if (intersect) { //if neutron is within ylimits of the mirror handle reflection/transmission

//first find the angle of the mirror. It is given by theta(x)=alpha*x^2+beta*x+gamma1, with theta(-l/2)=theta1, theta(0)=theta2, theta(l/2)=theta3

alpha1=2*(theta_1+theta_3-2*theta_2)/(length*length);
beta1=(theta_3-theta_1)/length;
gamma1=theta_2;

theta_m=alpha1*x*x+beta1*x+gamma1; // angle of mirror.

//The vector normal to the mirror is e_n= sin(theta)*e_x-cos(theta)*e_z

//find amplitude of v in direction of e_n:

sin_theta_m=sin(theta_m*DEG2RAD);
cos_theta_m=cos(theta_m*DEG2RAD);

v_n=sin_theta_m*vx-cos_theta_m*vz;


q=fabs(2.0*v_n*V2Q);

double R0=0.99;
double Qc=0.0217;
double m_value=m*0.9853+0.1978;
double W=-0.0002*m_value+0.0022;
double alpha=0.1204*m_value+5.0944;
double beta=-7.6251*m_value+68.1137;

if (m_value<=3)
{alpha=m_value;
beta=0;}





      /* Reflectivity (see component Guide). */
      if(m == 0)
        ABSORB;
      if (reflect && strlen(reflect) && strcmp(reflect,"NULL") && strcmp(reflect,"0"))
         Ref=Table_Value(pTable, q, 1);
      else {
          Ref = R0;
          if(q > Qc)
          {
            double arg = (q-m_value*Qc)/W;
            if(arg < 10)
              Ref *= .5*(1-tanh(arg))*(1-alpha*(q-Qc)+beta*(q-Qc)*(q-Qc)); //matches data from Swiss Neutronics
            else  Ref=0;
          }
      }
      if (Ref < 0) Ref=0;
      else if (Ref > 1) Ref=1;


//Now comes actual reflection/transmission
      if (!transmit) { //all neutrons are reflected
        if (!Ref) ABSORB;
        p *= Ref;

//handle reflection: change v_n -->-v_n

vx=old_vx*(cos_theta_m*cos_theta_m-sin_theta_m*sin_theta_m)+old_vz*(2*cos_theta_m*sin_theta_m);
vz=old_vx*(2*cos_theta_m*sin_theta_m)+old_vz*(sin_theta_m*sin_theta_m-cos_theta_m*cos_theta_m);

// printf("theta_m=%f, sin_theta_m=%f, cos_theta_m=%f, v_n=%f, old_vx=%f, vx=%f, old_vz=%f, vz=%f\n\n", theta_m, sin_theta_m, cos_theta_m, v_n, old_vx, vx, old_vz, vz);


        SCATTER; 
//printf("line 471.In mirror: x=%f,y=%f,z=%f,t=%f\n",x,y,z,t);
//printf("In mirror: old_vx=%f,old_vy=%f,old_vz=%f,vx=%f,vy=%f,vz=%f,v_n=%f\n",old_vx,old_vy,old_vz,vx,vy,vz,v_n);

      } else { //if neutrons can be transmitted



//calculate absorption.
// substrate
double lambda=(2*PI/V2K)/sqrt(vx*vx + vy*vy + vz*vz);
double sin_theta=lambda*q/(4*PI);

//double substrate_path_length=substrate_thickness/sin_theta;
//double coating_path_length=coating_thickness/sin_theta;

double sin_theta_c=Qc/(4*PI);

double theta_diff;
double substrate_path_length;
double coating_path_length;

double remaining_length_through_mirror;

int hit_back_mirror;

if (v_n>0) {
hit_back_mirror=1;} else{
hit_back_mirror=0;}

remaining_length_through_mirror=length/2-x;


if (sin_theta>sin_theta_c*lambda) {
theta_diff=sqrt(sin_theta*sin_theta-sin_theta_c*sin_theta_c*lambda*lambda);
coating_path_length=coating_thickness/theta_diff;
substrate_path_length=substrate_thickness/theta_diff;

	if (coating_path_length>remaining_length_through_mirror){
coating_path_length=remaining_length_through_mirror;
substrate_path_length=0; 
}

	if (substrate_path_length>remaining_length_through_mirror){
substrate_path_length=remaining_length_through_mirror; 
}












} else{

if (hit_back_mirror==0){ //neutron comes from front of mirror
substrate_path_length=0;
coating_path_length=remaining_length_through_mirror;
}else {//neutron comes from behind mirror

substrate_path_length=remaining_length_through_mirror;
coating_path_length=0;
}
}


double mu_substrate=0.0318/lambda+0.0055*lambda-0.0050; //unit: cm^-1
mu_substrate=mu_substrate*100; //unit: m^-1;

//For nickel and titanium coating, the following formular is used:
// mu = rho/m(atom)*sigma_a,thermal*lambda/lambda_thermal

// lambda_thermal=1.798 Å

// rho_nickel=8.908g/cm^3
// m(atom)_nickel=58.6934*1.661*10^-27 kg
// sigma_a,thermal_nickel=4.49*10^-28 m^2

// rho_titanium=4.506g/cm^3
// m(atom)_titanium=47.867*1.661*10^-27 kg
// sigma_a,thermal_titanium=6.09*10^-28 m^2

double Ni_coefficient=22.8180; 
double Ti_coefficient=19.1961;

double mu_coating=(0.5*Ni_coefficient+0.5*Ti_coefficient)*lambda; //it is roughly 50% nickel and 50% titanium



        // transmit when rand > R
        if (Ref == 0 || rand01() >= Ref) { //transmit
if (substrate_thickness>0){ p=p*exp(-mu_substrate*substrate_path_length-mu_coating*coating_path_length); //reduce weight of neutrons due to attenuation in the mirror
//x+=(coating_path_length+substrate_path_length)-(coating_thickness+substrate_thickness)/sin_theta;
//printf("xshift is %f \n",(coating_path_length+substrate_path_length)-(coating_thickness+substrate_thickness)/sin_theta);
} 
// printf("line 380\n");
/*
if (v_n>0) {
printf("neutron is transmitted from back of mirror. %f\n",exp(-mu_substrate*substrate_path_length-mu_coating*coating_path_length));
}else{
printf("neutron is transmitted from front of mirror. %f\n",exp(-mu_substrate*substrate_path_length-mu_coating*coating_path_length));
}
*/
} else {//neutron is reflected
		if (v_n>0 && substrate_thickness>0) { //if neutron comes from behind the mirror
// printf("neutron is reflected from back of mirror. %f\n",Ref*exp(-2*mu_substrate*substrate_path_length-2*mu_coating*coating_path_length));
			p=p*exp(-2*mu_substrate*substrate_path_length-2*mu_coating*coating_path_length);} //else{ //reduce weight of neutrons due to attenuation in the mirror
 // printf("neutron is reflected from front of mirror. %f\n", Ref);}
//handle reflection: change v_n -->-v_n
vx=old_vx*(cos_theta_m*cos_theta_m-sin_theta_m*sin_theta_m)+old_vz*(2*cos_theta_m*sin_theta_m);
vz=old_vx*(2*cos_theta_m*sin_theta_m)+old_vz*(sin_theta_m*sin_theta_m-cos_theta_m*cos_theta_m);
// printf("line 388\n");

}

// printf("theta_m=%f, sin_theta_m=%f, cos_theta_m=%f, v_n=%f, old_vx=%f, vx=%f, old_vz=%f, vz=%f\n\n", theta_m, sin_theta_m, cos_theta_m, v_n, old_vx, vx, old_vz, vz);

//printf("vxvx+vzvz=%f, oldvxoldvx+oldvzoldvz=%f", vx*vx+vz*vz, old_vx*old_vx+old_vz*old_vz);

        SCATTER; 
//printf("line 524.In mirror: x=%f,y=%f,z=%f,t=%f\n",x,y,z,t);
//printf("old_vx=%f,old_vy=%f,old_vz=%f,vx=%f,vy=%f,vz=%f,v_n=%f\n",old_vx,old_vy,old_vz,vx,vy,vz,v_n);
//after transmission or reflection
      } //end } else { after if (!transmit) {
    } 


 


} // end if (x >=-length/2 && x<=length/2)

// printf("intersect=%d\n",intersect);

   if (!intersect) {
      /* No intersection: restore neutron position. */
      x = old_x;
      y = old_y;

      z = old_z;
      t = old_t;
// printf("line 409\n");

    }
  

} //end if (dt>0)



  // EXTEND code here
  if (!strcmp(_comp->_name, "mirror_full_center")) {
if (SCATTERED) {
//{printf("I scatter\n");
guide_scatt=5; PROP_DT(1e-9); SCATTER; }
  }

  #undef reflect
  #undef focus_s
  #undef focus_e
  #undef mirror_start
  #undef guide_start
  #undef yheight
  #undef smallaxis
  #undef length
  #undef m
  #undef transmit
  #undef substrate_thickness
  #undef coating_thickness
  #undef theta_1
  #undef theta_2
  #undef theta_3
  #undef pTable
  return(_comp);
} /* class_Mirror_Curved_Bispectral_trace */

#pragma acc routine seq
_class_Mirror_Elliptic_Bispectral *class_Mirror_Elliptic_Bispectral_trace(_class_Mirror_Elliptic_Bispectral *_comp
  , _class_particle *_particle) {
  ABSORBED=SCATTERED=RESTORE=0;

  #define reflect (_comp->_parameters.reflect)
  #define focus_start_w (_comp->_parameters.focus_start_w)
  #define focus_end_w (_comp->_parameters.focus_end_w)
  #define focus_start_h (_comp->_parameters.focus_start_h)
  #define focus_end_h (_comp->_parameters.focus_end_h)
  #define mirror_start (_comp->_parameters.mirror_start)
  #define m (_comp->_parameters.m)
  #define smallaxis_w (_comp->_parameters.smallaxis_w)
  #define smallaxis_h (_comp->_parameters.smallaxis_h)
  #define length (_comp->_parameters.length)
  #define transmit (_comp->_parameters.transmit)
  #define substrate_thickness (_comp->_parameters.substrate_thickness)
  #define coating_thickness (_comp->_parameters.coating_thickness)
  #define pTable (_comp->_parameters.pTable)
double f; //half of distance between focal points
double asquared;
double a; //half of ellipse length
double b; //half of ellipse width

double xprime; //x in coordinates with center of ellipse at (xprime,zprime)=(0,0)
double ymirror; //height of the mirror


//Defining the mirror
double a1;
double b1;
double c1;

//solving the time the neutron hits the sample
double A, B, C, D, E, P, Q, R, U, V, I, J, K;

//finding rotation of mirror
double alpha1, beta1, gamma1;
double theta_m;
double xhi, zeta;


double b_w;
double f_w;
double asquared_w;
double A_w;
double B_w;
double C_w;
double determinant_w;


double b_h;
double f_h;
double asquared_h;
double xprime_h;


double xprime_w;
double zprime_w;
double asquare_z;
double bsquare_x;


double v_n; //speed of neutron perpendicular to surface

double Ref; //reflectivity

double dt;
double q;
  int intersect;

double discriminant;

double xprime_start_w;
double zprime_start_w;

double z_test;
double x_test;
double z_prime_test;

int hit_back_flag;

int prop_case;


double x_2;
double y_2;
double z_2;
double t_2;
int x_hit;
int x_hit_2;

double xprime_h_2;
double ymirror_2;
int intersect_2;

intersect=0;
x_2=0;
y_2=0;
z_2=0;
t_2=-1;
prop_case=0;
   double old_x = x, old_y = y, old_z = z, old_t=t, old_vx=vx, old_vz=vz, old_vy=vy;

// Check if neutron hits mirror. First find which z,x coordinates it hits.



//define the ellipse
b_w=smallaxis_w/2;
f_w=(focus_end_w-focus_start_w)*0.5;
asquared_w=f_w*f_w+b_w*b_w;

//in coordinate system of mirror: xprime_w(t)=xprime_w+vx*t. z-value is zprime(t)=b*sqrt(1-x'^2/a^2)+old_z+vz*t. This gives equation for t: A_w*t^2+B_w*t+C_w=0;

xprime_start_w=old_x-f_w-focus_start_w+mirror_start;
zprime_start_w=z;

A_w=b_w*b_w*vx*vx+asquared_w*vz*vz;
B_w=2*b_w*b_w*xprime_start_w*vx+2*asquared_w*old_z*vz;
C_w=b_w*b_w*xprime_start_w*xprime_start_w+asquared_w*old_z*old_z-asquared_w*b_w*b_w;

//this equation must now be solved for t;

if (A_w!=0){
determinant_w=B_w*B_w-4.0*A_w*C_w;
if (determinant_w>=0){
I=(-B_w-sqrt(determinant_w))/(2.0*A_w);
J=(-B_w+sqrt(determinant_w))/(2.0*A_w);

if (I<=0 && J<=0){dt = -1.0;} //both times are negative.
if (I<=0 && J>0 ){dt = J;} //set dt to only positive value.
if (I>0  && J<=0){dt = I;} //set dt to only positive value.

if (I>0  && J>0 ){prop_case=1; if (I>J) {dt=J;}else{dt=I;}} //set dt to smallest positive value  


} else {dt=-1.0;} //only complex solutions: set dt negative so no scattering

}else{ //end if (A!=0)  
if (B_w!=0){
dt=-C/B_w;
}else{ //end if (B!)=0
 printf("warning: A_w=B_w=C_w=0. Neutron is ignored\n"); }
} //end if (A!=0){}else{
//now intersection time has been found.

//printf("dt=%f\n",dt);

if (dt>0) { //if time is positive, propagate neutron to where it hits mirror. This is done without gravity.

    x += vx*dt;
    y += vy*dt;
    z += vz*dt;
    t += dt;


if (prop_case>0) //also check if neutron can hit mirror at second solution - it might not be in y-range for first solution
{
    x_2=x+vx*fabs(J-I); 
    y_2=y+vy*fabs(J-I);
    z_2=z+vz*fabs(J-I); 
    t_2=t+fabs(J-I); 

}else{
x_2=x;
y_2=y;
z_2=z;
t_2=t;
}

x_hit=(x>=0 &&x<=length);
x_hit_2=(x_2>=0 &&x_2<=length);

// printf("x=%f,y=%f,z=%f\n",x,y,z);
//if (x >=0 && x<=length){ //check if neutron is within x limits of the mirror. If so, check if it is within y limits.
if (x_hit || x_hit_2){

//define the ellipse
b_h=smallaxis_h/2;

f_h=(focus_end_h-focus_start_h)*0.5;

 asquared_h=f_h*f_h+b_h*b_h;

xprime_h=-f_h-focus_start_h+mirror_start+x; //xprime is the x-coordinate in a coordinate system centered at the center of the ellipse

ymirror=b_h*sqrt(1-xprime_h*xprime_h/asquared_h);


xprime_h_2=-f_h-focus_start_h+mirror_start+x_2; //xprime is the x-coordinate in a coordinate system centered at the center of the ellipse

ymirror_2=b_h*sqrt(1-xprime_h_2*xprime_h_2/asquared_h);


intersect = ( y>=-ymirror && y<=ymirror && x>=0 && x<=length && zprime_start_w+vz*dt>=0);

intersect_2 = ( y_2>=-ymirror_2 && y_2<=ymirror_2 && x_2>=0 && x_2<=length && zprime_start_w+vz*(dt+fabs(J-I))>=0);

if (!intersect && intersect_2){ //if neutron doesn't hit mirror with smallest t, but hits with largest t, propagte to largest t
intersect=intersect_2;
x=x_2;
y=y_2;
z=z_2;
t=t_2;
dt=t_2-old_t;
//printf("x=%f,y=%f,z=%f\n",x,y,z);
}

    if (intersect) { //if neutron is within ylimits of the mirror handle reflection/transmission


//now perform reflection. 

//First find out if neutron hits front or back of mirror: propagate backwards a bit and check if neutron is outside ellipse: if so, it hits back of mirror


b_w=smallaxis_w/2.0;

f_w=(focus_end_w-focus_start_w)*0.5;
 asquared_w=f_w*f_w+b_w*b_w;

z_test=zprime_start_w+vz*(dt-1e-6);
x_test=xprime_start_w+vx*(dt-1e-6);
z_prime_test=b_w*sqrt(1-x_test*x_test/asquared_w);

//find velocity in q direction.

xprime_w=xprime_start_w+vx*dt;
zprime_w=zprime_start_w+vz*dt;
asquare_z=asquared_w*zprime_w;
bsquare_x=b_w*b_w*xprime_w;


zeta=(asquare_z)/(sqrt(asquare_z*asquare_z+bsquare_x*bsquare_x));
xhi=-(bsquare_x)/(sqrt(asquare_z*asquare_z+bsquare_x*bsquare_x));

//printf("z_test=%f, z_prime_test=%f\n",z_test,z_prime_test);

if (z_test>z_prime_test) {
hit_back_flag=1;
}
//printf("vx=%f, vz=%f, vy=%f, xhi=%f, zeta=%f, prop_case=%d\n",vx,vz,vy,xhi,zeta,prop_case);
v_n=-xhi*vx+zeta*vz;

q=fabs(2.0*v_n*V2Q);


 //Reflectivity parameters calculated from SWISS neutronics data.
double R0=0.99;
double Qc=0.0217;
double m_value=m*0.9853+0.1978;
double W=-0.0002*m_value+0.0022;
double alpha=0.1204*m_value+5.0944;
double beta=-7.6251*m_value+68.1137;

if (m_value<=3)
{alpha=m_value;
beta=0;}


      /* Reflectivity (see component Guide). */
      if(m == 0)
        ABSORB;
      if (reflect && strlen(reflect) && strcmp(reflect,"NULL") && strcmp(reflect,"0"))
         Ref=Table_Value(pTable, q, 1);
      else {
          Ref = R0;
          if(q > Qc)
          {
            double arg = (q-m_value*Qc)/W;
            if(arg < 10)
              Ref *= .5*(1-tanh(arg))*(1-alpha*(q-Qc)+beta*(q-Qc)*(q-Qc)); //matches data from Swiss Neutronics
            else  Ref=0;
          }
      }
      if (Ref < 0) Ref=0;
      else if (Ref > 1) Ref=1;


//Now comes actual reflection/transmission
      if (!transmit) { //all neutrons are reflected
//printf("v_n=%f,q=%f, Ref=%f, lambda=%f, theta=%f, 1p_before=%f",v_n,q, Ref, (2*PI/V2K)/sqrt(vx*vx + vy*vy + vz*vz),asin((2*PI/V2K)/sqrt(vx*vx + vy*vy + vz*vz)*q/(4*PI))*RAD2DEG,p);
        if (!Ref) ABSORB;
        p *= Ref;
//printf("p_after=%f\n",p);
//handle reflection: change v_n -->-v_n

vx=old_vx*(zeta*zeta-xhi*xhi)+old_vz*(2*zeta*xhi);
vz=+old_vx*(2*zeta*xhi)+old_vz*(xhi*xhi-zeta*zeta);


        SCATTER; //after transmission or reflection

      } else { //if neutrons can be transmitted



//calculate absorption.
// substrate
double lambda=(2*PI/V2K)/sqrt(vx*vx + vy*vy + vz*vz);
double sin_theta=lambda*q/(4*PI);
double substrate_path_length=substrate_thickness/sin_theta;
double mu_substrate=0.0318/lambda+0.0055*lambda-0.0050; //unit: cm^-1
mu_substrate=mu_substrate*100; //unit: m^-1;
               
//For nickel and titanium coating, the following formular is used:
// mu = rho/m(atom)*sigma_a,thermal*lambda/lambda_thermal

// lambda_thermal=1.798 Å

// rho_nickel=8.908g/cm^3
// m(atom)_nickel=58.6934*1.661*10^-27 kg
// sigma_a,thermal_nickel=4.49*10^-28 m^2

// rho_titanium=4.506g/cm^3
// m(atom)_titanium=47.867*1.661*10^-27 kg
// sigma_a,thermal_titanium=6.09*10^-28 m^2

double coating_path_length=coating_thickness/sin_theta;
double Ni_coefficient=22.8180; 
double Ti_coefficient=19.1961;

double mu_coating=(0.5*Ni_coefficient+0.5*Ti_coefficient)*lambda; //it is roughly 50% nickel and 50% titanium



        // transmit when rand > R
        if (Ref == 0 || rand01() >= Ref) { //transmit
if (substrate_thickness>0){ p=p*exp(-mu_substrate*substrate_path_length-mu_coating*coating_path_length);} //reduce weight of neutrons due to attenuation in the mirror

} else {//neutron is reflected
		if (hit_back_flag==1 && substrate_thickness>0) { //if neutron comes from behind the mirror
			p=p*exp(-2*mu_substrate*substrate_path_length-2*mu_coating*coating_path_length);} //reduce weight of neutrons due to attenuation in the mirror
//handle reflection: change v_n -->-v_n

vx=old_vx*(zeta*zeta-xhi*xhi)-old_vz*(2*zeta*xhi);
vz=-old_vx*(2*zeta*xhi)+old_vz*(xhi*xhi-zeta*zeta);
}

//printf("p_before=%f, q=%f",p,q);
        SCATTER; //after transmission or reflection
//printf("p_after=%f\n",p);
      } //end } else { after if (!transmit) {
    } 


 


} // end if (x >=-length/2 && x<=length/2)


   if (!intersect) {
      /* No intersection: restore neutron position. */
      x = old_x;
      y = old_y;

      z = old_z;
      t = old_t;


    }
  

} //end if (dt>0)



  // EXTEND code here
  if (!strcmp(_comp->_name, "guide_right")) {
if (SCATTERED){
//{printf("I scatter\n");
guide_scatt=3; PROP_DT(1e-9); SCATTER; }
  }
  if (!strcmp(_comp->_name, "guide_bottom")) {
if (SCATTERED){
//printf("I scatter\n");
guide_scatt=2; PROP_DT(1e-9); SCATTER; }
  }
  if (!strcmp(_comp->_name, "guide_top")) {
if (SCATTERED){
//{printf("I scatter\n");
guide_scatt=1; PROP_DT(1e-9); SCATTER; }
  }
  if (!strcmp(_comp->_name, "guide_Left")) {
if (SCATTERED){
//{printf("I scatter\n");
guide_scatt=4; PROP_DT(1e-9); SCATTER; }
  }

  #undef reflect
  #undef focus_start_w
  #undef focus_end_w
  #undef focus_start_h
  #undef focus_end_h
  #undef mirror_start
  #undef m
  #undef smallaxis_w
  #undef smallaxis_h
  #undef length
  #undef transmit
  #undef substrate_thickness
  #undef coating_thickness
  #undef pTable
  return(_comp);
} /* class_Mirror_Elliptic_Bispectral_trace */

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

/* *****************************************************************************
* instrument 'ESS_2001_bispectral' TRACE
***************************************************************************** */

#pragma acc routine seq
int raytrace(_class_particle* _particle) { /* called by mccode_main for ESS_2001_bispectral:TRACE */

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
      class_Progress_bar_trace(_Origin, _particle); /* contains EXTEND code */
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component Origin [1] */
    if (!ABSORBED && _particle->_index == 2) {
      /* component ArmForGuideRight=Arm() [2] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_ArmForGuideRight->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _ArmForGuideRight->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_ArmForGuideRight->_position_relative, _ArmForGuideRight->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      _particle->_index++;
    } /* end component ArmForGuideRight [2] */
    if (!ABSORBED && _particle->_index == 3) {
      /* component ArmForGuideBottom=Arm() [3] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_ArmForGuideBottom->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _ArmForGuideBottom->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_ArmForGuideBottom->_position_relative, _ArmForGuideBottom->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      _particle->_index++;
    } /* end component ArmForGuideBottom [3] */
    if (!ABSORBED && _particle->_index == 4) {
      /* component ArmForGuideTop=Arm() [4] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_ArmForGuideTop->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _ArmForGuideTop->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_ArmForGuideTop->_position_relative, _ArmForGuideTop->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      _particle->_index++;
    } /* end component ArmForGuideTop [4] */
    if (!ABSORBED && _particle->_index == 5) {
      /* component ArmForGuideLeft=Arm() [5] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_ArmForGuideLeft->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _ArmForGuideLeft->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_ArmForGuideLeft->_position_relative, _ArmForGuideLeft->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      _particle->_index++;
    } /* end component ArmForGuideLeft [5] */
    if (!ABSORBED && _particle->_index == 6) {
      /* component cold_source=ESS_moderator_long_2001() [6] */
  #define size (_cold_source->_parameters.size)
  #define l_low (_cold_source->_parameters.l_low)
  #define l_high (_cold_source->_parameters.l_high)
  #define dist (_cold_source->_parameters.dist)
  #define xw (_cold_source->_parameters.xw)
  #define yh (_cold_source->_parameters.yh)
  #define freq (_cold_source->_parameters.freq)
  #define T (_cold_source->_parameters.T)
  #define tau (_cold_source->_parameters.tau)
  #define tau1 (_cold_source->_parameters.tau1)
  #define tau2 (_cold_source->_parameters.tau2)
  #define d (_cold_source->_parameters.d)
  #define n (_cold_source->_parameters.n)
  #define n2 (_cold_source->_parameters.n2)
  #define chi2 (_cold_source->_parameters.chi2)
  #define I0 (_cold_source->_parameters.I0)
  #define I2 (_cold_source->_parameters.I2)
  #define branch1 (_cold_source->_parameters.branch1)
  #define branch2 (_cold_source->_parameters.branch2)
  #define branch_tail (_cold_source->_parameters.branch_tail)
  #define twopulses (_cold_source->_parameters.twopulses)
  #define target_index (_cold_source->_parameters.target_index)
  #define l_range (_cold_source->_parameters.l_range)
  #define w_mult (_cold_source->_parameters.w_mult)
  #define branchframe (_cold_source->_parameters.branchframe)
  #define tx (_cold_source->_parameters.tx)
  #define ty (_cold_source->_parameters.ty)
  #define tz (_cold_source->_parameters.tz)
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_cold_source->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _cold_source->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_cold_source->_position_relative, _cold_source->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      if (( flag == 1 )) // conditional WHEN execution
      class_ESS_moderator_long_2001_trace(_cold_source, _particle); /* contains EXTEND code */
      if (_particle->_restore)
        *_particle = _particle_save;
  #undef size
  #undef l_low
  #undef l_high
  #undef dist
  #undef xw
  #undef yh
  #undef freq
  #undef T
  #undef tau
  #undef tau1
  #undef tau2
  #undef d
  #undef n
  #undef n2
  #undef chi2
  #undef I0
  #undef I2
  #undef branch1
  #undef branch2
  #undef branch_tail
  #undef twopulses
  #undef target_index
  #undef l_range
  #undef w_mult
  #undef branchframe
  #undef tx
  #undef ty
  #undef tz
      _particle->_index++;
    } /* end component cold_source [6] */
    if (!ABSORBED && _particle->_index == 7) {
      /* component thermal_source=ESS_moderator_long_2001() [7] */
  #define size (_thermal_source->_parameters.size)
  #define l_low (_thermal_source->_parameters.l_low)
  #define l_high (_thermal_source->_parameters.l_high)
  #define dist (_thermal_source->_parameters.dist)
  #define xw (_thermal_source->_parameters.xw)
  #define yh (_thermal_source->_parameters.yh)
  #define freq (_thermal_source->_parameters.freq)
  #define T (_thermal_source->_parameters.T)
  #define tau (_thermal_source->_parameters.tau)
  #define tau1 (_thermal_source->_parameters.tau1)
  #define tau2 (_thermal_source->_parameters.tau2)
  #define d (_thermal_source->_parameters.d)
  #define n (_thermal_source->_parameters.n)
  #define n2 (_thermal_source->_parameters.n2)
  #define chi2 (_thermal_source->_parameters.chi2)
  #define I0 (_thermal_source->_parameters.I0)
  #define I2 (_thermal_source->_parameters.I2)
  #define branch1 (_thermal_source->_parameters.branch1)
  #define branch2 (_thermal_source->_parameters.branch2)
  #define branch_tail (_thermal_source->_parameters.branch_tail)
  #define twopulses (_thermal_source->_parameters.twopulses)
  #define target_index (_thermal_source->_parameters.target_index)
  #define l_range (_thermal_source->_parameters.l_range)
  #define w_mult (_thermal_source->_parameters.w_mult)
  #define branchframe (_thermal_source->_parameters.branchframe)
  #define tx (_thermal_source->_parameters.tx)
  #define ty (_thermal_source->_parameters.ty)
  #define tz (_thermal_source->_parameters.tz)
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_thermal_source->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _thermal_source->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_thermal_source->_position_relative, _thermal_source->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      if (( flag == 0 )) // conditional WHEN execution
      class_ESS_moderator_long_2001_trace(_thermal_source, _particle); /* contains EXTEND code */
      if (_particle->_restore)
        *_particle = _particle_save;
  #undef size
  #undef l_low
  #undef l_high
  #undef dist
  #undef xw
  #undef yh
  #undef freq
  #undef T
  #undef tau
  #undef tau1
  #undef tau2
  #undef d
  #undef n
  #undef n2
  #undef chi2
  #undef I0
  #undef I2
  #undef branch1
  #undef branch2
  #undef branch_tail
  #undef twopulses
  #undef target_index
  #undef l_range
  #undef w_mult
  #undef branchframe
  #undef tx
  #undef ty
  #undef tz
      _particle->_index++;
    } /* end component thermal_source [7] */
    if (!ABSORBED && _particle->_index == 8) {
      /* component ColdFocus=Arm() [8] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_ColdFocus->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _ColdFocus->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_ColdFocus->_position_relative, _ColdFocus->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Arm_trace(_ColdFocus, _particle); /* contains EXTEND code */
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component ColdFocus [8] */
    if (!ABSORBED && _particle->_index == 9) {
      /* component ArmMidOne=Arm() [9] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_ArmMidOne->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _ArmMidOne->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_ArmMidOne->_position_relative, _ArmMidOne->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Arm_trace(_ArmMidOne, _particle); /* contains EXTEND code */
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component ArmMidOne [9] */
    if (!ABSORBED && _particle->_index == 10) {
      /* component mirror_full_center=Mirror_Curved_Bispectral() [10] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_mirror_full_center->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _mirror_full_center->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_mirror_full_center->_position_relative, _mirror_full_center->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Mirror_Curved_Bispectral_trace(_mirror_full_center, _particle); /* contains EXTEND code */
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component mirror_full_center [10] */
    if (!ABSORBED && _particle->_index == 11) {
      /* component ArmForNeutronPropState_2=Arm() [11] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_ArmForNeutronPropState_2->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _ArmForNeutronPropState_2->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_ArmForNeutronPropState_2->_position_relative, _ArmForNeutronPropState_2->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Arm_trace(_ArmForNeutronPropState_2, _particle); /* contains EXTEND code */
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component ArmForNeutronPropState_2 [11] */
    if (!ABSORBED && _particle->_index == 12) {
      /* component guide_right=Mirror_Elliptic_Bispectral() [12] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_guide_right->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _guide_right->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_guide_right->_position_relative, _guide_right->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Mirror_Elliptic_Bispectral_trace(_guide_right, _particle); /* contains EXTEND code */
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component guide_right [12] */
    if (!ABSORBED && _particle->_index == 13) {
      /* component ArmForNeutronPropState_4=Arm() [13] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_ArmForNeutronPropState_4->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _ArmForNeutronPropState_4->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_ArmForNeutronPropState_4->_position_relative, _ArmForNeutronPropState_4->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Arm_trace(_ArmForNeutronPropState_4, _particle); /* contains EXTEND code */
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component ArmForNeutronPropState_4 [13] */
    if (!ABSORBED && _particle->_index == 14) {
      /* component guide_bottom=Mirror_Elliptic_Bispectral() [14] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_guide_bottom->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _guide_bottom->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_guide_bottom->_position_relative, _guide_bottom->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Mirror_Elliptic_Bispectral_trace(_guide_bottom, _particle); /* contains EXTEND code */
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component guide_bottom [14] */
    if (!ABSORBED && _particle->_index == 15) {
      /* component ArmForNeutronPropState_5=Arm() [15] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_ArmForNeutronPropState_5->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _ArmForNeutronPropState_5->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_ArmForNeutronPropState_5->_position_relative, _ArmForNeutronPropState_5->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Arm_trace(_ArmForNeutronPropState_5, _particle); /* contains EXTEND code */
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component ArmForNeutronPropState_5 [15] */
    if (!ABSORBED && _particle->_index == 16) {
      /* component cold_lambda_guidestart=L_monitor() [16] */
  #define nL (_cold_lambda_guidestart->_parameters.nL)
  #define filename (_cold_lambda_guidestart->_parameters.filename)
  #define xmin (_cold_lambda_guidestart->_parameters.xmin)
  #define xmax (_cold_lambda_guidestart->_parameters.xmax)
  #define ymin (_cold_lambda_guidestart->_parameters.ymin)
  #define ymax (_cold_lambda_guidestart->_parameters.ymax)
  #define xwidth (_cold_lambda_guidestart->_parameters.xwidth)
  #define yheight (_cold_lambda_guidestart->_parameters.yheight)
  #define Lmin (_cold_lambda_guidestart->_parameters.Lmin)
  #define Lmax (_cold_lambda_guidestart->_parameters.Lmax)
  #define restore_neutron (_cold_lambda_guidestart->_parameters.restore_neutron)
  #define L_N (_cold_lambda_guidestart->_parameters.L_N)
  #define L_p (_cold_lambda_guidestart->_parameters.L_p)
  #define L_p2 (_cold_lambda_guidestart->_parameters.L_p2)
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_cold_lambda_guidestart->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _cold_lambda_guidestart->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_cold_lambda_guidestart->_position_relative, _cold_lambda_guidestart->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      if (( flag == 1 )) // conditional WHEN execution
      class_L_monitor_trace(_cold_lambda_guidestart, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
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
      _particle->_index++;
    } /* end component cold_lambda_guidestart [16] */
    if (!ABSORBED && _particle->_index == 17) {
      /* component thermal_lambda_guidestart=L_monitor() [17] */
  #define nL (_thermal_lambda_guidestart->_parameters.nL)
  #define filename (_thermal_lambda_guidestart->_parameters.filename)
  #define xmin (_thermal_lambda_guidestart->_parameters.xmin)
  #define xmax (_thermal_lambda_guidestart->_parameters.xmax)
  #define ymin (_thermal_lambda_guidestart->_parameters.ymin)
  #define ymax (_thermal_lambda_guidestart->_parameters.ymax)
  #define xwidth (_thermal_lambda_guidestart->_parameters.xwidth)
  #define yheight (_thermal_lambda_guidestart->_parameters.yheight)
  #define Lmin (_thermal_lambda_guidestart->_parameters.Lmin)
  #define Lmax (_thermal_lambda_guidestart->_parameters.Lmax)
  #define restore_neutron (_thermal_lambda_guidestart->_parameters.restore_neutron)
  #define L_N (_thermal_lambda_guidestart->_parameters.L_N)
  #define L_p (_thermal_lambda_guidestart->_parameters.L_p)
  #define L_p2 (_thermal_lambda_guidestart->_parameters.L_p2)
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_thermal_lambda_guidestart->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _thermal_lambda_guidestart->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_thermal_lambda_guidestart->_position_relative, _thermal_lambda_guidestart->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      if (( flag == 0 )) // conditional WHEN execution
      class_L_monitor_trace(_thermal_lambda_guidestart, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
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
      _particle->_index++;
    } /* end component thermal_lambda_guidestart [17] */
    if (!ABSORBED && _particle->_index == 18) {
      /* component lambda_guidestart=L_monitor() [18] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_lambda_guidestart->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _lambda_guidestart->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_lambda_guidestart->_position_relative, _lambda_guidestart->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_L_monitor_trace(_lambda_guidestart, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component lambda_guidestart [18] */
    if (!ABSORBED && _particle->_index == 19) {
      /* component guide_top=Mirror_Elliptic_Bispectral() [19] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_guide_top->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _guide_top->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_guide_top->_position_relative, _guide_top->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Mirror_Elliptic_Bispectral_trace(_guide_top, _particle); /* contains EXTEND code */
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component guide_top [19] */
    if (!ABSORBED && _particle->_index == 20) {
      /* component ArmForNeutronPropState_6=Arm() [20] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_ArmForNeutronPropState_6->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _ArmForNeutronPropState_6->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_ArmForNeutronPropState_6->_position_relative, _ArmForNeutronPropState_6->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Arm_trace(_ArmForNeutronPropState_6, _particle); /* contains EXTEND code */
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component ArmForNeutronPropState_6 [20] */
    if (!ABSORBED && _particle->_index == 21) {
      /* component guide_Left=Mirror_Elliptic_Bispectral() [21] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_guide_Left->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _guide_Left->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_guide_Left->_position_relative, _guide_Left->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Mirror_Elliptic_Bispectral_trace(_guide_Left, _particle); /* contains EXTEND code */
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component guide_Left [21] */
    if (!ABSORBED && _particle->_index == 22) {
      /* component ArmForNeutronPropState_7=Arm() [22] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_ArmForNeutronPropState_7->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _ArmForNeutronPropState_7->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_ArmForNeutronPropState_7->_position_relative, _ArmForNeutronPropState_7->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Arm_trace(_ArmForNeutronPropState_7, _particle); /* contains EXTEND code */
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component ArmForNeutronPropState_7 [22] */
    if (!ABSORBED && _particle->_index == 23) {
      /* component ArmMidTwo=Arm() [23] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_ArmMidTwo->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _ArmMidTwo->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_ArmMidTwo->_position_relative, _ArmMidTwo->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Arm_trace(_ArmMidTwo, _particle); /* contains EXTEND code */
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component ArmMidTwo [23] */
    if (!ABSORBED && _particle->_index == 24) {
      /* component ArmForNeutronPropState_8=Arm() [24] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_ArmForNeutronPropState_8->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _ArmForNeutronPropState_8->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_ArmForNeutronPropState_8->_position_relative, _ArmForNeutronPropState_8->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Arm_trace(_ArmForNeutronPropState_8, _particle); /* contains EXTEND code */
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component ArmForNeutronPropState_8 [24] */
    if (!ABSORBED && _particle->_index == 25) {
      /* component ArmMidThree=Arm() [25] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_ArmMidThree->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _ArmMidThree->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_ArmMidThree->_position_relative, _ArmMidThree->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      if (( guide_scatt > 0 )) {/* conditional JUMP to ArmMidOne */
        _particle->_index=8;
        flag_nocoordschange=1; // pass corrdinate transformations when jumping
      }
      _particle->_index++;
    } /* end component ArmMidThree [25] */
    if (!ABSORBED && _particle->_index == 26) {
      /* component ArmExit=Arm() [26] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_ArmExit->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _ArmExit->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_ArmExit->_position_relative, _ArmExit->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      _particle->_index++;
    } /* end component ArmExit [26] */
    if (_particle->_index > 26)
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
* instrument 'ESS_2001_bispectral' and components SAVE
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



int save(FILE *handle) { /* called by mccode_main for ESS_2001_bispectral:SAVE */
  if (!handle) siminfo_init(NULL);

  /* call iteratively all components SAVE */
  class_Progress_bar_save(_Origin);















  class_L_monitor_save(_cold_lambda_guidestart);

  class_L_monitor_save(_thermal_lambda_guidestart);

  class_L_monitor_save(_lambda_guidestart);









  if (!handle) siminfo_close(); 

  return(0);
} /* save */

/* *****************************************************************************
* instrument 'ESS_2001_bispectral' and components FINALLY
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



int finally(void) { /* called by mccode_main for ESS_2001_bispectral:FINALLY */
  siminfo_init(NULL);
  save(siminfo_file); /* save data when simulation ends */

  /* call iteratively all components FINALLY */
  class_Progress_bar_finally(_Origin);















  class_L_monitor_finally(_cold_lambda_guidestart);

  class_L_monitor_finally(_thermal_lambda_guidestart);

  class_L_monitor_finally(_lambda_guidestart);









  siminfo_close(); 

  return(0);
} /* finally */

/* *****************************************************************************
* instrument 'ESS_2001_bispectral' and components DISPLAY
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

_class_Arm *class_Arm_display(_class_Arm *_comp
) {
  /* A bit ugly; hard-coded dimensions. */
  
  line(0,0,0,0.2,0,0);
  line(0,0,0,0,0.2,0);
  line(0,0,0,0,0,0.2);
  return(_comp);
} /* class_Arm_display */

_class_ESS_moderator_long_2001 *class_ESS_moderator_long_2001_display(_class_ESS_moderator_long_2001 *_comp
) {
  #define size (_comp->_parameters.size)
  #define l_low (_comp->_parameters.l_low)
  #define l_high (_comp->_parameters.l_high)
  #define dist (_comp->_parameters.dist)
  #define xw (_comp->_parameters.xw)
  #define yh (_comp->_parameters.yh)
  #define freq (_comp->_parameters.freq)
  #define T (_comp->_parameters.T)
  #define tau (_comp->_parameters.tau)
  #define tau1 (_comp->_parameters.tau1)
  #define tau2 (_comp->_parameters.tau2)
  #define d (_comp->_parameters.d)
  #define n (_comp->_parameters.n)
  #define n2 (_comp->_parameters.n2)
  #define chi2 (_comp->_parameters.chi2)
  #define I0 (_comp->_parameters.I0)
  #define I2 (_comp->_parameters.I2)
  #define branch1 (_comp->_parameters.branch1)
  #define branch2 (_comp->_parameters.branch2)
  #define branch_tail (_comp->_parameters.branch_tail)
  #define twopulses (_comp->_parameters.twopulses)
  #define target_index (_comp->_parameters.target_index)
  #define l_range (_comp->_parameters.l_range)
  #define w_mult (_comp->_parameters.w_mult)
  #define branchframe (_comp->_parameters.branchframe)
  #define tx (_comp->_parameters.tx)
  #define ty (_comp->_parameters.ty)
  #define tz (_comp->_parameters.tz)
  
  rectangle("xy", 0, 0, 0, size, size);
  #undef size
  #undef l_low
  #undef l_high
  #undef dist
  #undef xw
  #undef yh
  #undef freq
  #undef T
  #undef tau
  #undef tau1
  #undef tau2
  #undef d
  #undef n
  #undef n2
  #undef chi2
  #undef I0
  #undef I2
  #undef branch1
  #undef branch2
  #undef branch_tail
  #undef twopulses
  #undef target_index
  #undef l_range
  #undef w_mult
  #undef branchframe
  #undef tx
  #undef ty
  #undef tz
  return(_comp);
} /* class_ESS_moderator_long_2001_display */

_class_Mirror_Curved_Bispectral *class_Mirror_Curved_Bispectral_display(_class_Mirror_Curved_Bispectral *_comp
) {
  #define reflect (_comp->_parameters.reflect)
  #define focus_s (_comp->_parameters.focus_s)
  #define focus_e (_comp->_parameters.focus_e)
  #define mirror_start (_comp->_parameters.mirror_start)
  #define guide_start (_comp->_parameters.guide_start)
  #define yheight (_comp->_parameters.yheight)
  #define smallaxis (_comp->_parameters.smallaxis)
  #define length (_comp->_parameters.length)
  #define m (_comp->_parameters.m)
  #define transmit (_comp->_parameters.transmit)
  #define substrate_thickness (_comp->_parameters.substrate_thickness)
  #define coating_thickness (_comp->_parameters.coating_thickness)
  #define theta_1 (_comp->_parameters.theta_1)
  #define theta_2 (_comp->_parameters.theta_2)
  #define theta_3 (_comp->_parameters.theta_3)
  #define pTable (_comp->_parameters.pTable)



/*
if (xcenter==0){

xstart=0;
xend=length;
xprime_start=-a+mirror_start+xstart;
ystart=b*sqrt(1-xprime_start*xprime_start/asquared);

xprime_end=-a+mirror_start+xend;
yend=b*sqrt(1-xprime_end*xprime_end/asquared);

}
*/

/*
if (xcenter==1){
xstart=-length/2;
xend=length/2;


xprime_start=-a+mirror_start+xstart+length/2;
ystart=b*sqrt(1-xprime_start*xprime_start/asquared);

xprime_end=-a+mirror_start+xend+length/2;
yend=b*sqrt(1-xprime_end*xprime_end/asquared);
}

line(xstart,-ystart,0,xstart,ystart,0);
line(xend,-yend,0,xend,yend,0);
line(xstart,-ystart,0,xend,-yend,0);
line(xstart,ystart,0,xend,yend,0);
*/

double xstart;
double xend;
double xprime_start;
double ystart;

double xprime_end;
double yend;

double focus_2;
double focus_1;
double b;
double f;
double asquared;
double a;


int n_lines;
int j=0;
double xprimepos[51];
double ypos[51];
double x_plot[51];
double zpos[51];
double xstep;



focus_2=focus_e-mirror_start; //focus in local coordinates
focus_1=focus_s-mirror_start;

b=smallaxis/2;

f=(focus_2-focus_1)*0.5;
 asquared=f*f+b*b;
 a=sqrt(asquared);


xstart=-length/2;
xend=length/2;

n_lines=50;

xstep=length/n_lines;







double a1, b1, c1;
double tan_theta_1;
double tan_theta_2;
double tan_theta_3;



//mirror is defined by z(x)=a1x^3+b1x^2+c1x+d1, with dz/dx|x=-length/2=tan(theta_1), dz/dx|x=0=tan(theta_2), dz/dx|x=length/2=tan(theta3), z(0)=0. (d1=0)

tan_theta_1=tan(theta_1*DEG2RAD);
tan_theta_2=tan(theta_2*DEG2RAD);
tan_theta_3=tan(theta_3*DEG2RAD);


a1=2.0/3.0*(tan_theta_1+tan_theta_3-2.0*tan_theta_2)/(length*length);
b1=(tan_theta_3-tan_theta_1)/(2.0*length);
c1=tan_theta_2;



for (j=0; j<n_lines+1; j++)
{
xprimepos[j]=-f-focus_s+mirror_start+length/2+xstart+xstep*j;

ypos[j]=b*sqrt(1-xprimepos[j]*xprimepos[j]/asquared); //correct

if (guide_start>mirror_start){
if (  xstart+xstep*j<-length/2+guide_start-mirror_start) {
ypos[j]=yheight/2;
}
}



// ypos[j]=b*sqrt(1-xprimepos[j]*xprimepos[j]/(f*f)); //following convention in Kaspar's elliptic guide..
// printf("xprimepos[j]=%f,f*f=%f, ypos[j]=%f\n",xprimepos[j],f*f,ypos[j]);

x_plot[j]=xstart+xstep*j;
zpos[j]=a1*x_plot[j]*x_plot[j]*x_plot[j]+b1*x_plot[j]*x_plot[j]+c1*x_plot[j];
}

for (j=0; j<n_lines; j++)
{
line(x_plot[j], -ypos[j], zpos[j], x_plot[j+1], -ypos[j+1],zpos[j+1]);
line(x_plot[j], ypos[j], zpos[j], x_plot[j+1], ypos[j+1],zpos[j+1]);
}


line(x_plot[0],-ypos[0],zpos[0],x_plot[0],ypos[0],zpos[0]);
line(x_plot[50],-ypos[50],zpos[50],x_plot[50],ypos[50],zpos[50]);




//printf("ypos0=%f xpos0=%f ypos50=%f xpos50=%f",ypos[0], x_plot[0], ypos[50], x_plot[50]);

/*  double xmax, xmin, ymax, ymin;
  

  if (center == 0) {
    xmax= x1; xmin=0;
    ymax= yheight; ymin=0;
  } else {
    xmax= x1/2; xmin=-xmax;
    ymax= yheight/2; ymin=-ymax;
  }
  multiline(5, (double)xmin, (double)ymin, 0.0,
               (double)xmax, (double)ymin, 0.0,
               (double)xmax, (double)ymax, 0.0,
               (double)xmin, (double)ymax, 0.0,
               (double)xmin, (double)ymin, 0.0);
*/
  #undef reflect
  #undef focus_s
  #undef focus_e
  #undef mirror_start
  #undef guide_start
  #undef yheight
  #undef smallaxis
  #undef length
  #undef m
  #undef transmit
  #undef substrate_thickness
  #undef coating_thickness
  #undef theta_1
  #undef theta_2
  #undef theta_3
  #undef pTable
  return(_comp);
} /* class_Mirror_Curved_Bispectral_display */

_class_Mirror_Elliptic_Bispectral *class_Mirror_Elliptic_Bispectral_display(_class_Mirror_Elliptic_Bispectral *_comp
) {
  #define reflect (_comp->_parameters.reflect)
  #define focus_start_w (_comp->_parameters.focus_start_w)
  #define focus_end_w (_comp->_parameters.focus_end_w)
  #define focus_start_h (_comp->_parameters.focus_start_h)
  #define focus_end_h (_comp->_parameters.focus_end_h)
  #define mirror_start (_comp->_parameters.mirror_start)
  #define m (_comp->_parameters.m)
  #define smallaxis_w (_comp->_parameters.smallaxis_w)
  #define smallaxis_h (_comp->_parameters.smallaxis_h)
  #define length (_comp->_parameters.length)
  #define transmit (_comp->_parameters.transmit)
  #define substrate_thickness (_comp->_parameters.substrate_thickness)
  #define coating_thickness (_comp->_parameters.coating_thickness)
  #define pTable (_comp->_parameters.pTable)




double xstart;
double xend;
double xprime_start;
double ystart;

double xprime_end;
double yend;

double focus_2_h;
double focus_1_h;
double b_h;
double f_h;
double asquared_h;

double focus_2_w;
double focus_1_w;
double b_w;
double f_w;
double asquared_w;


int n_lines;
int j=0;
double xprimepos[51];
double ypos[51];
double x_plot[51];
double zpos[51];
double xstep;
double xprime_w[51];



focus_2_h=focus_end_h; //focus in local coordinates
focus_1_h=focus_start_h;

b_h=smallaxis_h/2.0;

f_h=(focus_2_h-focus_1_h)*0.5;
 asquared_h=f_h*f_h+b_h*b_h;



focus_2_w=focus_end_w; //focus in local coordinates
focus_1_w=focus_start_w;

b_w=smallaxis_w/2.0;

f_w=(focus_2_w-focus_1_w)*0.5;
 asquared_w=f_w*f_w+b_w*b_w;



xstart=0;
xend=length;

n_lines=50;

xstep=length/n_lines;













for (j=0; j<n_lines+1; j++)
{
xprimepos[j]=-f_h-focus_start_h+mirror_start+xstart+xstep*j;

ypos[j]=b_h*sqrt(1-xprimepos[j]*xprimepos[j]/asquared_h); //correct


x_plot[j]=xstart+xstep*j;



xprime_w[j]=-f_w-focus_start_w+mirror_start+xstart+xstep*j; //xprime is the x-coordinate in a coordinate system centered at the center of the ellipse

zpos[j]=b_w*sqrt(1-xprime_w[j]*xprime_w[j]/asquared_w);
//printf("xprime=%f, zpos=%f\n",xprime_w[j], zpos[j]);
}

for (j=0; j<n_lines; j++)
{
line(x_plot[j], -ypos[j], zpos[j], x_plot[j+1], -ypos[j+1],zpos[j+1]);
line(x_plot[j], ypos[j], zpos[j], x_plot[j+1], ypos[j+1],zpos[j+1]);
}


line(x_plot[0],-ypos[0],zpos[0],x_plot[0],ypos[0],zpos[0]);
line(x_plot[50],-ypos[50],zpos[50],x_plot[50],ypos[50],zpos[50]);




//printf("ypos0=%f xpos0=%f ypos50=%f xpos50=%f",ypos[0], x_plot[0], ypos[50], x_plot[50]);

/*  double xmax, xmin, ymax, ymin;
  

  if (center == 0) {
    xmax= x1; xmin=0;
    ymax= yheight; ymin=0;
  } else {
    xmax= x1/2; xmin=-xmax;
    ymax= yheight/2; ymin=-ymax;
  }
  multiline(5, (double)xmin, (double)ymin, 0.0,
               (double)xmax, (double)ymin, 0.0,
               (double)xmax, (double)ymax, 0.0,
               (double)xmin, (double)ymax, 0.0,
               (double)xmin, (double)ymin, 0.0);
*/
  #undef reflect
  #undef focus_start_w
  #undef focus_end_w
  #undef focus_start_h
  #undef focus_end_h
  #undef mirror_start
  #undef m
  #undef smallaxis_w
  #undef smallaxis_h
  #undef length
  #undef transmit
  #undef substrate_thickness
  #undef coating_thickness
  #undef pTable
  return(_comp);
} /* class_Mirror_Elliptic_Bispectral_display */

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


  #undef magnify
  #undef line
  #undef dashed_line
  #undef multiline
  #undef rectangle
  #undef box
  #undef circle
  #undef cylinder
  #undef sphere

int display(void) { /* called by mccode_main for ESS_2001_bispectral:DISPLAY */
  printf("MCDISPLAY: start\n");

  /* call iteratively all components DISPLAY */
  class_Progress_bar_display(_Origin);

  class_Arm_display(_ArmForGuideRight);

  class_Arm_display(_ArmForGuideBottom);

  class_Arm_display(_ArmForGuideTop);

  class_Arm_display(_ArmForGuideLeft);

  class_ESS_moderator_long_2001_display(_cold_source);

  class_ESS_moderator_long_2001_display(_thermal_source);

  class_Arm_display(_ColdFocus);

  class_Arm_display(_ArmMidOne);

  class_Mirror_Curved_Bispectral_display(_mirror_full_center);

  class_Arm_display(_ArmForNeutronPropState_2);

  class_Mirror_Elliptic_Bispectral_display(_guide_right);

  class_Arm_display(_ArmForNeutronPropState_4);

  class_Mirror_Elliptic_Bispectral_display(_guide_bottom);

  class_Arm_display(_ArmForNeutronPropState_5);

  class_L_monitor_display(_cold_lambda_guidestart);

  class_L_monitor_display(_thermal_lambda_guidestart);

  class_L_monitor_display(_lambda_guidestart);

  class_Mirror_Elliptic_Bispectral_display(_guide_top);

  class_Arm_display(_ArmForNeutronPropState_6);

  class_Mirror_Elliptic_Bispectral_display(_guide_Left);

  class_Arm_display(_ArmForNeutronPropState_7);

  class_Arm_display(_ArmMidTwo);

  class_Arm_display(_ArmForNeutronPropState_8);

  class_Arm_display(_ArmMidThree);

  class_Arm_display(_ArmExit);

  printf("MCDISPLAY: end\n");

  return(0);
} /* display */

void* _getvar_parameters(char* compname)
/* enables settings parameters based use of the GETPAR macro */
{
  if (!strcmp(compname, "Origin")) return (void *) &(_Origin->_parameters);
  if (!strcmp(compname, "ArmForGuideRight")) return (void *) &(_ArmForGuideRight->_parameters);
  if (!strcmp(compname, "ArmForGuideBottom")) return (void *) &(_ArmForGuideBottom->_parameters);
  if (!strcmp(compname, "ArmForGuideTop")) return (void *) &(_ArmForGuideTop->_parameters);
  if (!strcmp(compname, "ArmForGuideLeft")) return (void *) &(_ArmForGuideLeft->_parameters);
  if (!strcmp(compname, "cold_source")) return (void *) &(_cold_source->_parameters);
  if (!strcmp(compname, "thermal_source")) return (void *) &(_thermal_source->_parameters);
  if (!strcmp(compname, "ColdFocus")) return (void *) &(_ColdFocus->_parameters);
  if (!strcmp(compname, "ArmMidOne")) return (void *) &(_ArmMidOne->_parameters);
  if (!strcmp(compname, "mirror_full_center")) return (void *) &(_mirror_full_center->_parameters);
  if (!strcmp(compname, "ArmForNeutronPropState_2")) return (void *) &(_ArmForNeutronPropState_2->_parameters);
  if (!strcmp(compname, "guide_right")) return (void *) &(_guide_right->_parameters);
  if (!strcmp(compname, "ArmForNeutronPropState_4")) return (void *) &(_ArmForNeutronPropState_4->_parameters);
  if (!strcmp(compname, "guide_bottom")) return (void *) &(_guide_bottom->_parameters);
  if (!strcmp(compname, "ArmForNeutronPropState_5")) return (void *) &(_ArmForNeutronPropState_5->_parameters);
  if (!strcmp(compname, "cold_lambda_guidestart")) return (void *) &(_cold_lambda_guidestart->_parameters);
  if (!strcmp(compname, "thermal_lambda_guidestart")) return (void *) &(_thermal_lambda_guidestart->_parameters);
  if (!strcmp(compname, "lambda_guidestart")) return (void *) &(_lambda_guidestart->_parameters);
  if (!strcmp(compname, "guide_top")) return (void *) &(_guide_top->_parameters);
  if (!strcmp(compname, "ArmForNeutronPropState_6")) return (void *) &(_ArmForNeutronPropState_6->_parameters);
  if (!strcmp(compname, "guide_Left")) return (void *) &(_guide_Left->_parameters);
  if (!strcmp(compname, "ArmForNeutronPropState_7")) return (void *) &(_ArmForNeutronPropState_7->_parameters);
  if (!strcmp(compname, "ArmMidTwo")) return (void *) &(_ArmMidTwo->_parameters);
  if (!strcmp(compname, "ArmForNeutronPropState_8")) return (void *) &(_ArmForNeutronPropState_8->_parameters);
  if (!strcmp(compname, "ArmMidThree")) return (void *) &(_ArmMidThree->_parameters);
  if (!strcmp(compname, "ArmExit")) return (void *) &(_ArmExit->_parameters);
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

/* end of generated C code ESS_2001_bispectral.c */
