/* Automatically generated file. Do not edit. 
 * Format:     ANSI C source code
 * Creator:    McStas <http://www.mcstas.org>
 * Instrument: ISIS_CRISP.instr (ISIS_CRISP)
 * Date:       Sat Oct 12 09:04:53 2019
 * File:       ISIS_CRISP.c
 * CFLAGS= -lgsl -lgslcblas
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
* Start of instrument 'ISIS_CRISP' generated code
***************************************************************************** */

#ifdef MC_TRACE_ENABLED
int traceenabled = 1;
#else
int traceenabled = 0;
#endif
#define MCSTAS "/usr/share/mcstas/3.0-dev/"
int   defaultmain         = 1;
char  instrument_name[]   = "ISIS_CRISP";
char  instrument_source[] = "ISIS_CRISP.instr";
char *instrument_exe      = NULL; /* will be set to argv[0] in main */
char  instrument_code[]   = "Instrument ISIS_CRISP source code ISIS_CRISP.instr is not embedded in this executable.\n  Use --source option when running McStas.\n";

int main(int argc, char *argv[]){return mccode_main(argc, argv);}

/* *****************************************************************************
* instrument 'ISIS_CRISP' and components DECLARE
***************************************************************************** */

/* Instrument parameters: structure and a table for the initialisation
   (Used in e.g. inputparse and I/O function (e.g. detector_out) */

struct _struct_instrument_parameters {
  MCNUM _glen;
  MCNUM _flen;
  MCNUM _w1;
  MCNUM _vm;
  MCNUM _FRAC;
};
typedef struct _struct_instrument_parameters _class_instrument_parameters;

struct _instrument_struct {
  char   _name[256]; /* the name of this instrument e.g. 'ISIS_CRISP' */
/* Counters per component instance */
  double counter_AbsorbProp[20]; /* absorbed events in PROP routines */
  double counter_N[20], counter_P[20], counter_P2[20]; /* event counters after each component instance */
  _class_particle _trajectory[20]; /* current trajectory for STORE/RESTORE */
/* Components position table (absolute and relative coords) */
  Coords _position_relative[20]; /* positions of all components */
  Coords _position_absolute[20];
  _class_instrument_parameters _parameters; /* instrument parameters */
} _instrument_var;
struct _instrument_struct *instrument = & _instrument_var;
#pragma acc declare create ( _instrument_var )
#pragma acc declare create ( instrument )

int numipar = 5;
struct mcinputtable_struct mcinputtable[] = {
  "glen", &(_instrument_var._parameters._glen), instr_type_double, "1.4", 
  "flen", &(_instrument_var._parameters._flen), instr_type_double, "0.4", 
  "w1", &(_instrument_var._parameters._w1), instr_type_double, "0.05", 
  "vm", &(_instrument_var._parameters._vm), instr_type_double, "3.0", 
  "FRAC", &(_instrument_var._parameters._FRAC), instr_type_double, "0", 
  NULL, NULL, instr_type_double, ""
};


/* ************************************************************************** */
/*             SHARE user declarations for all components                     */
/* ************************************************************************** */

/* Shared user declarations for all components types 'ISIS_moderator'. */
typedef struct
{
int nEnergy;        ///< Number of energy bins
int nTime;          ///< number of time bins

double* TimeBin;    ///< Time bins
double* EnergyBin;  ///< Energy bins

double** Flux;       ///< Flux per bin (integrated)
    double* EInt;        ///< Integrated Energy point
    double Total;        ///< Integrated Total
  } Source;

  /* New functions */

  int cmdnumberD(char *,double*);
  int cmdnumberI(char *,int*,const int);
  double polInterp(double*,double*,int,double);
  FILE* openFile(char*);
  FILE* openFileTest(char*);
  int readHtable(FILE*,const double,const double, Source*);
  int timeStart(char*);
  int timeEnd(char*);
  int energyBin(char*,double,double,double*,double*);
  int notComment(char*);
  double strArea(double dist, double rtmodX, double rtmodY, double focus_xw, double focus_yh);


double** matrix(const int m,const int n)
 /*!
   Determine a double matrix
 */
{
  int i;
  double* pv;
  double** pd;

  if (m<1) return 0;
  if (n<1) return 0;
  pv = (double*) malloc(m*n*sizeof(double));
  pd = (double**) malloc(m*sizeof(double*));
  if (!pd)
    {
      fprintf(stderr,"No room for matrix!\n");
      exit(1);
    }
  for (i=0;i<m;i++)
    pd[i]=pv + (i*n);
  return pd;
}


double polInterp(double* X,double* Y,int Psize,double Aim)
  /*!
    returns the interpolated polynomial between Epnts
    and the integration
    \param X :: X coordinates
    \param Y :: Y coordinates
    \param Psize :: number of valid point in array to use
    \param Aim :: Aim point to intepolate result (X coordinate)
    \returns Energy point
  */
{
  double out,errOut;         /* out put variables */
  double C[Psize],D[Psize];
  double testDiff,diff;

  double w,den,ho,hp;           /* intermediate variables */
  int i,m,ns;


  ns=0;
  diff=fabs(Aim-X[0]);
  C[0]=Y[0];
  D[0]=Y[0];
  for(i=1;i<Psize;i++)
    {
      testDiff=fabs(Aim-X[i]);
      if (diff>testDiff)
	{
	  ns=i;
	  diff=testDiff;
	}
      C[i]=Y[i];
      D[i]=Y[i];
    }

  out=Y[ns];
  ns--;              /* Now can be -1 !!!! */

  for(m=1;m<Psize;m++)
    {
      for(i=0;i<Psize-m;i++)
	{
	  ho=X[i]-Aim;
	  hp=X[i+m]-Aim;
	  w=C[i+1]-D[i];
	  /*	  den=ho-hp;  -- test !=0.0 */
	  den=w/(ho-hp);
	  D[i]=hp*den;
	  C[i]=ho*den;
	}

      errOut= (2*(ns+1)<(Psize-m)) ? C[ns+1] : D[ns--];
      out+=errOut;
    }
  return out;
}

int binSearch(int Npts,double* AR,double V)
  /*!
    Object is to find the point in
    array AR, closest to the value V
    Checked for ordered array returns lower of backeting objects
  */
{
  int klo,khi,k;
  if (Npts<=0)
    return 0;
  if (V>AR[Npts-1])
    return Npts;

  if(AR[0]>0.0)AR[0]=0.0;

  if (V<AR[0])
    {
      // if(AR[0]>0.0)AR[0]=0.0;
      fprintf(stderr,"here");
      return 0;
    }
  klo=0;
  khi= Npts-1;
  while (khi-klo >1)
    {
      k=(khi+klo) >> 1;    // quick division by 2
      if (AR[k]>V)
	khi=k;
      else
	klo=k;
    }
  return khi;
}

int cmdnumberD(char *mc,double* num)
 /*!
   \returns 1 on success 0 on failure
 */
{
  int i,j;
  char* ss;
  char **endptr;
  double nmb;
  int len;

  len=strlen(mc);
  j=0;

  for(i=0;i<len && mc[i] &&
	(mc[i]=='\t' || mc[i]==' '  || mc[i]==',');i++);
  if(i==len || !mc[i]) return 0;
  ss=malloc(sizeof(char)*(len+1));

  for(;i<len && mc[i]!='\n' && mc[i]
	&& mc[i]!='\t' && mc[i]!=' ' && mc[i]!=',';i++)
    {
      ss[j]=mc[i];
      j++;
    }
  if (!j)
    {
      free(ss);
      return 0;         //This should be impossible
    }
  ss[j]=0;
  endptr=malloc(sizeof(char*));
  nmb = strtod(ss,endptr);
  if (*endptr != ss+j)
    {
      free(endptr);
      free(ss);
      return 0;
    }
  *num = (double) nmb;
  for(j=0;j<i && mc[j];j++)
    mc[j]=' ';
  free(endptr);
  free(ss);
  return 1;
}

int notComment(char* Line)
 /*!
   \returns 0 on a comment, 1 on a non-comment
 */
{
  int len,i;

  len=strlen(Line);
  for(i=0;i<len && isspace(Line[i]);i++);

  if (!Line[i] || Line[i]=='c' || Line[i]=='C' ||
      Line[i]=='!' || Line[i]=='#')
    return 0;
  return 1;
}

int timeStart(char* Line)
 /*!
   Search for a word time at the start of
   the line.
   \param Line :: Line to search
   \returns 1 on success 0 on failure
 */
{
  int len,i;

  len=strlen(Line);
  for(i=0;i<len && isspace(Line[i]);i++);
  if (len-i<4) return 0;
  return (strncmp(Line+i,"time",4)) ? 0 : 1;
}

int timeEnd(char* Line)
 /*!
   Search for a word time at the start of
   the line.
   \param Line :: Line to search
   \returns 1 on success 0 on failure
 */
{
  int len,i;

  len=strlen(Line);
  for(i=0;i<len && isspace(Line[i]);i++);
  if (len-i<5) return 0;
  return (strncmp(Line+i,"total",5)) ? 0 : 1;
}

int energyBin(char* Line,double Einit,double Eend,double* Ea,double* Eb)
     /*!
       Search for a word "energy bin:" at the start of
       the line. Then separte off the energy bin values
       \param Line :: Line to search
       \param Ea :: first energy bin [meV]
       \param Eb :: second energy bin [meV]
       \returns 1 on success 0 on failure
     */
{
  int len,i;
  double A,B;

  len=strlen(Line);
  for(i=0;i<len && isspace(Line[i]);i++);
  if (len-i<11) return 0;


  if (strncmp(Line+i,"energy bin:",11))
    return 0;

  i+=11;
  if (!cmdnumberD(Line+i,&A))
    return 0;
  // remove 'to'
  for(;i<len-1 && Line[i]!='o';i++);
  i++;
  if (!cmdnumberD(Line+i,&B))
    return 0;
  A*=1e9;
  B*=1e9;
  *Ea=A;
  *Eb=B;
  if (*Eb>Einit && *Ea<Eend)
    return 1;
  return 0;
}

double calcFraction(double EI,double EE,double Ea,double Eb)
 /*!
   Calculate the fraction of the bin between Ea -> Eb
   that is encompassed by EI->EE
 */
{
  double frac;
  double dRange;

  if (EI>Eb)
    return 0.0;
  if (EE<Ea)
    return 0.0;

  dRange=Eb-Ea;
  frac=(EI>Ea) ? (Eb-EI)/dRange : 1.0;


  frac-=(EE<Eb) ? (Eb-EE)/dRange : 0.0;

  //  if(frac != 1.0)
  //  fprintf(stderr,"frac %g, Ea %g,Eb %g, EI %g, EE %g\n",frac,Ea,Eb,EI,EE);

  return frac;
}

int readHtable(FILE* TFile,const double Einit,const double Eend, Source *TS)
     /*!
       Process a general h.o file to create an integrated
       table of results from Einit -> Eend
       \param Einit :: inital Energy
       \parma Eend  :: final energy
     */
{
  char ss[255];          /* BIG space for line */
  double Ea,Eb;
  double T,D;
  double Efrac;          // Fraction of an Energy Bin
  int Ftime;             // time Flag
  int eIndex;             // energy Index
  int tIndex;             // time Index
  double Tsum;           // Running integration
  double Efraction;      // Amount to use for an energy/time bin

  // extern Source TS;

  int DebugCnt;
  int i;
  /*!
    Status Flag::
    Ftime=1 :: [time ] Reading Time : Data : Err [Exit on Total]


    /*
    Double Read File to determine how many bins and
    memery size
  */
  if (!TFile) return(0);
  Ea=0.0;
  Eb=0.0;
  fprintf(stderr,"Energy == %g %g\n",Einit,Eend);
  eIndex= -1;
  DebugCnt=0;
  Ftime=0;
  tIndex=0;
  TS->nTime=0;
  TS->nEnergy=0;
  // Read file and get time bins
  while(fgets(ss,255,TFile) && Eend>Ea)
    {
      if (notComment(ss))
	{
	  DebugCnt++;
          if (!Ftime)
	    {
	      if (energyBin(ss,Einit,Eend,&Ea,&Eb))
		{
		  if (eIndex==0)
		    TS->nTime=tIndex;
		  eIndex++;
		}
	      else if (timeStart(ss))
		{
		  Ftime=1;
		  tIndex=0;
		}
	    }
	  else  // In the time section
	    {
	      if (timeEnd(ss))     // Found "total"
		Ftime=0;
	      else
		{
		  // Need to read the line in the case of first run
		  if (TS->nTime==0)
		    {
		      if (cmdnumberD(ss,&T) &&
			  cmdnumberD(ss,&D))
			tIndex++;
		    }
		}
	    }
	}
    }
  // Plus 2 since we have a 0 counter and we have missed the last line.
  TS->nEnergy=eIndex+2;
  if (!TS->nTime && tIndex)
    TS->nTime=tIndex;
  // printf("tIndex %d %d %d %d \n",tIndex,eIndex,TS->nEnergy,TS->nTime);

  /* SECOND TIME THROUGH:: */
  rewind(TFile);

  TS->Flux=matrix(TS->nEnergy,TS->nTime);
  TS->EInt=(double*) malloc(TS->nEnergy*sizeof(double));
  TS->TimeBin=(double*) malloc(TS->nTime*sizeof(double));
  TS->EnergyBin=(double*) malloc(TS->nEnergy*sizeof(double));

  Tsum=0.0;
  Ea=0.0;
  Eb=0.0;
  eIndex=-1;
  DebugCnt=0;
  Ftime=0;
  tIndex=0;
  TS->EInt[0]=0.0;
  // Read file and get time bins
  while(fgets(ss,255,TFile) && Eend>Ea)
    {
      if (notComment(ss))
	{
	  DebugCnt++;
          if (!Ftime)
	    {
	      if (energyBin(ss,Einit,Eend,&Ea,&Eb))
		{
		  eIndex++;
		  TS->EnergyBin[eIndex]=(Einit>Ea) ? Einit : Ea;
		  Efraction=calcFraction(Einit,Eend,Ea,Eb);
		  Ftime++;
		}
	    }
	  else if (Ftime==1)
	    {
	      if (timeStart(ss))
		{
		  Ftime=2;
		  tIndex=0;
		}
	    }

	  else           // In the time section
	    {
	      if (timeEnd(ss))     // Found "total"
		{
		  Ftime=0;
		  TS->EInt[eIndex+1]=Tsum;
		}
	      else
		{
		  // Need to read the line in the case of first run
		  if (cmdnumberD(ss,&T) &&
		      cmdnumberD(ss,&D))
		    {
		      TS->TimeBin[tIndex]=T/1e8;     // convert Time into second (from shakes)
		      Tsum+=D*Efraction;
		      TS->Flux[eIndex][tIndex]=Tsum;
		      tIndex++;
		    }
		}
	    }
	}
    }

  TS->EnergyBin[eIndex+1]=Eend;
  TS->Total=Tsum;

  //  printf("tIndex %d %d %d \n",tIndex,eIndex,TS.nTime);
  //printf("Tsum %g \n",Tsum);
  //fprintf(stderr,"ebin1 ebinN %g %g\n",TS.EnergyBin[0],TS.EnergyBin[TS.nEnergy-1]);

  return 1;
} // readHtable

void getPoint(double* TV,double* EV,double* lim1, double* lim2, Source TS)
 /*!
   Calculate the Time and Energy
   by sampling the file.
   Uses TS table to find the point
   \param TV ::
   \param EV ::
   \param lim1 ::
   \param lim2 ::
 */
{
  int i;

  // extern Source TS;
  double R0,R1,R,Rend;
  int Epnt;       ///< Points to the next higher index of the neutron integral
  int Tpnt;
  int iStart,iEnd;
  double TRange,Tspread;
  double Espread,Estart;
  double *EX;

  // So that lowPoly+highPoly==maxPoly
  const int maxPoly=6;
  const int highPoly=maxPoly/2;
  const int lowPoly=maxPoly-highPoly;

  // static int testVar=0;

  R0=rand01();
  /* if (testVar==0)
    {
    R0=1.0e-8;
    testVar=1;
    }
  */
  Rend=R=TS.Total*R0;
  // This gives Eint[Epnt-1] > R > Eint[Epnt]
  Epnt=binSearch(TS.nEnergy-1,TS.EInt,R);

  //      if (Epnt < 0)
  //   Epnt=1;
  Tpnt=binSearch(TS.nTime-1,TS.Flux[Epnt-1],R);
  //  fprintf(stderr,"TBoundaryX == %12.6e %12.6e \n",TS.TimeBin[Tpnt-1],TS.TimeBin[Tpnt]);
  //  fprintf(stderr,"TFlux == %12.6e %12.6e %12.6e \n\n",TS.Flux[Epnt-1][Tpnt-1],R,TS.Flux[Epnt-1][Tpnt]);
  //  if (Epnt == -1)
  //{
  //    Epnt=0;
  // fprintf(stderr,"\n Rvals == %g %d %d %g\n",R,Epnt,Tpnt,TS.TimeBin[0]);
  //  fprintf(stderr,"EInt == %d %12.6e %12.6e %12.6e %12.6e \n",Epnt,TS.EInt[Epnt-1],R,TS.EInt[Epnt],TS.EInt[Epnt+1]);
  // printf("EBoundary == %12.6e %12.6e \n",TS.EnergyBin[Epnt],TS.EnergyBin[Epnt+1]);

  //  fprintf(stderr,"TFlux == %12.6e %12.6e %12.6e \n\n",TS.Flux[Epnt+1][Tpnt],R,TS.Flux[Epnt+1][Tpnt+1]);
  // }

  if(R < TS.Flux[Epnt-1][Tpnt-1] || R >TS.Flux[Epnt-1][Tpnt] )
    {
      fprintf(stderr,"outside bin limits Tpnt/Epnt problem  %12.6e %12.6e %12.6e \n",TS.Flux[Epnt-1][Tpnt-1],R,TS.Flux[Epnt-1][Tpnt]);
    }

  if(Epnt == 0)
    {
      Estart=0.0;
      Espread=TS.EInt[0];
      *EV=TS.EnergyBin[1];
    }
  else
    {
      Estart=TS.EInt[Epnt-1];
      Espread=TS.EInt[Epnt]-TS.EInt[Epnt-1];
      *EV=TS.EnergyBin[Epnt+1];
    }

  if (Tpnt==0 || Epnt==0)
    {
      fprintf(stderr,"BIG ERROR WITH Tpnt: %d and Epnt: %d\n",Tpnt,Epnt);
      exit(1);
    }
  if (Tpnt==TS.nTime)
    {
      fprintf(stderr,"BIG ERROR WITH Tpnt and Epnt\n");
      exit(1);
      *TV=0.0;
      Tspread=TS.Flux[Epnt-1][0]-TS.EInt[Epnt-1];
      TRange=TS.TimeBin[0];
      R-=TS.EInt[Epnt-1];
    }
  else
    {
      *TV=TS.TimeBin[Tpnt-1];
      TRange=TS.TimeBin[Tpnt]-TS.TimeBin[Tpnt-1];
      Tspread=TS.Flux[Epnt-1][Tpnt]-TS.Flux[Epnt-1][Tpnt-1];
      R-=TS.Flux[Epnt-1][Tpnt-1];
    }
  //  printf("R == %12.6e\n",R);
  R/=Tspread;
  //  printf("R == %12.6e\n",R);
  *TV+=TRange*R;


  R1=TS.EInt[Epnt-1]+Espread*rand01();
  iStart=Epnt>lowPoly ? Epnt-lowPoly : 0;                  // max(Epnt-halfPoly,0)
  iEnd=TS.nEnergy>Epnt+highPoly ? Epnt+highPoly : TS.nEnergy-1;  // min(nEnergy-1,Epnt+highPoly

  *EV=polInterp(TS.EInt+iStart,TS.EnergyBin+iStart,1+iEnd-iStart,R1);

  //  fprintf(stderr,"Energy == %d %d %12.6e %12.6e \n",iStart,iEnd,R1,*EV);
  //  fprintf(stderr,"bins == %12.6e %12.6e %12.6e %12.6e \n",TS.EnergyBin[iStart],TS.EnergyBin[iEnd],
  //	  TS.EInt[Epnt],TS.EInt[Epnt-1]);

    if(*TV < TS.TimeBin[Tpnt-1] || *TV > TS.TimeBin[Tpnt])
    {
      fprintf(stderr,"%d Tpnt %d Tval %g Epnt %d \n",TS.nTime,Tpnt,*TV,Epnt);
      fprintf(stderr,"TBoundary == %12.6e,%g , %12.6e \n\n",TS.TimeBin[Tpnt-1],*TV,TS.TimeBin[Tpnt]);
    }


  if(*EV < *lim1 || *EV > *lim2)
    {
      fprintf(stderr,"outside boundaries\n Epnt= %d, Tpnt= %d binlo %g|%g| binhi %g \n",Epnt,Tpnt,TS.EnergyBin[Epnt-1],*EV,TS.EnergyBin[Epnt]);


      fprintf(stderr,"TS == %g %g :: %d %d \n",TS.EInt[Epnt-1],TS.EInt[Epnt],iStart,iEnd);
      fprintf(stderr,"Points (%g) == ",R1);
      for(i=0;i<iEnd-iStart;i++)
	fprintf(stderr," %g %g",*(TS.EInt+i+iStart),*(TS.EnergyBin+iStart+i));
      fprintf(stderr,"\n");

      //fprintf(stderr,"energy value %g\n",*EV);
      //  fprintf(stderr,"TFlux == %12.6e %12.6e %12.6e \n",TS.Flux[Epnt-1][Tpnt-1],Rend,TS.Flux[Epnt-1][Tpnt]);
    }
  return;
} // getPoint

int cmdnumberI(char *mc,int* num,const int len)
  /*!
    \param mc == character string to use
    \param num :: Place to put output
    \param len == length of the character string to process
    returns 1 on success and 0 on failure
  */
{
  int i,j;
  char* ss;
  char **endptr;
  double nmb;

      if (len<1)
	return 0;
      j=0;

      for(i=0;i<len && mc[i] &&
	    (mc[i]=='\t' || mc[i]==' '  || mc[i]==',');i++);
      if(i==len || !mc[i]) return 0;
      ss=malloc(sizeof(char)*(len+1));
      /*  char *ss=new char[len+1]; */
      for(;i<len && mc[i]!='\n' && mc[i]
	    && mc[i]!='\t' && mc[i]!=' ' && mc[i]!=',';i++)
	{
	  ss[j]=mc[i];
	  j++;
	}
      if (!j)
	{
	  free(ss);
	  return 0;         //This should be impossible
	}
      ss[j]=0;
      endptr=malloc(sizeof(char*));
      nmb = strtod(ss,endptr);
      if (*endptr != ss+j)
	{
	  free(endptr);
	  free(ss);
	  return 0;
	}
      *num = (double) nmb;
      for(j=0;j<i && mc[j];j++)
	mc[j]=' ';
      free(endptr);
      free(ss);
      return 1;
    }


  FILE* openFile(char* FileName)
    {
      FILE* efile=0;
      char ss[1024];
      if (!FileName) return(NULL);
      
      if (!efile && getenv("MCTABLES")) {
        /* Is MCTABLES set, files located there? */
        sprintf(ss, "%s%c%s", getenv("MCTABLES"), MC_PATHSEP_C, FileName);
        efile=fopen(ss,"r");
      }

      /* Is the file located in working dir? */
      if (!efile) {
      sprintf(ss,"%s", FileName);
        efile=fopen(FileName,"r");
      }
      if (!efile) {
        /* Try locating the file in ./ISIS_tables library */
        sprintf(ss,"%s%c%s%c%s", ".", MC_PATHSEP_C, "ISIS_tables", MC_PATHSEP_C, FileName);
        efile=fopen(ss,"r");
      }
      if (!efile) {
        /* Try locating the file in the MCSTAS contrib/ISIS_tables library */
        sprintf(ss, "%s%c%s%c%s%c%s", getenv("MCSTAS") ? getenv("MCSTAS") : MCSTAS, 
          MC_PATHSEP_C, "contrib", MC_PATHSEP_C, "ISIS_tables", MC_PATHSEP_C, FileName);
        efile=fopen(ss,"r");
      }
      if (!efile) {
        /* Try locating the file in the MCSTAS data library */
        sprintf(ss, "%s%c%s%c%s", getenv("MCSTAS") ? getenv("MCSTAS") : MCSTAS, 
          MC_PATHSEP_C, "contrib", MC_PATHSEP_C, FileName);
        efile=fopen(ss,"r");
      }
      if (!efile) {
        /* Try locating the file in the MCSTAS data/ISIS_tables library */
        sprintf(ss, "%s%c%s%c%s%c%s", getenv("MCSTAS") ? getenv("MCSTAS") : MCSTAS, 
          MC_PATHSEP_C, "data", MC_PATHSEP_C, "ISIS_tables", MC_PATHSEP_C, FileName);
        efile=fopen(ss,"r");
      }
      if (!efile) {
        /* Try locating the file in the MCSTAS data library */
        sprintf(ss, "%s%c%s%c%s", getenv("MCSTAS") ? getenv("MCSTAS") : MCSTAS, 
          MC_PATHSEP_C, "data", MC_PATHSEP_C, FileName);
        efile=fopen(ss,"r");
      }
      if (!efile) { /* Still no file - die! */
        fprintf(stderr,"ISIS_moderator: ERROR: Could not read Etable file %s.\n", FileName);
        fprintf(stderr,"                Please check your McStas installation and/or MCTABLES/ISIS_tables setting!\n");
        exit(1);
      }
      else 
        printf("Opening -- %s\n",ss);
      return efile;
    }

  double strArea(double dist, double rtmodX, double rtmodY, double focus_xw, double focus_yh)
    {
      /*
	 Returns the mean Str view of the viewport
	 This integrates over each point on the window focus_xw to focus_yh
	 View port is symmetric so use only 1/4 of the view
	 for the calcuation.
	 Control Values rtmodY rtmodX focus_xw focus_yh
      */

      double A;
      double Vx,Vy;        // view temp points
      double Mx,My;        // moderator x,y
      double D2;           // Distance ^2
      int i,j,aa,bb;       // loop variables

      D2=dist*dist;
      A=0.0;

      for(i=0;i<50;i++)              // Mod X
	{
	  Mx=i*rtmodX/100.0;
	  for(j=0;j<50;j++)         // Mod Y
	    {
	      My=j*rtmodY/100.0;
	      // Position on moderator == (Mx,My)
	      for(aa=-50;aa<51;aa++)  //view port
		for(bb=-50;bb<51;bb++)
		  {
		    Vx=aa*focus_xw/101.0;
		    Vy=bb*focus_yh/101.0;
		    A+=1.0/((Mx-Vx)*(Mx-Vx)+(My-Vy)*(My-Vy)+D2);
		  }
	    }
	}
	//change to Mx*My
      A*=(rtmodY*rtmodX)/(10201.0*2500.0);
      // Correct for the area of the viewport. (tables are per cm^2)
      A*=focus_xw*focus_yh*10000;

      fprintf(stderr,"Viewport == %g %g Moderator size == (%g * %g) m^2 \n",focus_xw,focus_yh,rtmodX,rtmodY);
      fprintf(stderr,"Dist == %g (metres) \n",dist);
      fprintf(stderr,"Viewport Solid angle == %g str\n",A/(focus_xw*focus_yh*10000));
      fprintf(stderr,"Solid angle used == %g str\n",A);
      return A;
    }

/* Shared user declarations for all components types 'Guide_tapering'. */
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


/* Shared user declarations for all components types 'Multilayer_Sample'. */
#ifndef GSL_VERSION
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_cblas.h>
#endif


/* ************************************************************************** */
/*             End of SHARE user declarations for all components              */
/* ************************************************************************** */


/* ********************** component definition declarations. **************** */

/* component a1=Arm() [1] DECLARE */
/* Parameter definition for component type 'Arm' */
struct _struct_Arm_parameters {
  char Arm_has_no_parameters;
}; /* _struct_Arm_parameters */
typedef struct _struct_Arm_parameters _class_Arm_parameters;

/* Parameters for component type 'Arm' */
struct _struct_Arm {
  char     _name[256]; /* e.g. a1 */
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
_class_Arm _a1_var;
_class_Arm *_a1 = &_a1_var;
#pragma acc declare create ( _a1_var )
#pragma acc declare create ( _a1 )

/* component isis_source=ISIS_moderator() [2] DECLARE */
/* Parameter definition for component type 'ISIS_moderator' */
struct _struct_ISIS_moderator_parameters {
  /* Component type 'ISIS_moderator' setting parameters */
  char Face[16384];
  MCNUM Emin;
  MCNUM Emax;
  MCNUM dist;
  MCNUM focus_xw;
  MCNUM focus_yh;
  MCNUM xwidth;
  MCNUM yheight;
  MCNUM CAngle;
  MCNUM SAC;
  MCNUM Lmin;
  MCNUM Lmax;
  long target_index;
  MCNUM verbose;
  /* Component type 'ISIS_moderator' private parameters */
  /* Component type 'ISIS_moderator' DECLARE code stored as structure members */
  #include <ctype.h>
                        

  double p_in;                                              
  int Tnpts;                                                    
  double scaleSize;                                                                
  double angleArea;                                  
  double Nsim;	                                                     
  int Ncount;                                                
  Source TS;

                        

  double rtE0,rtE1;                                                                                       
  double rtmodX,rtmodY;                                                                                      

  int TargetStation;
  double CurrentWeight;
}; /* _struct_ISIS_moderator_parameters */
typedef struct _struct_ISIS_moderator_parameters _class_ISIS_moderator_parameters;

/* Parameters for component type 'ISIS_moderator' */
struct _struct_ISIS_moderator {
  char     _name[256]; /* e.g. isis_source */
  char     _type[256]; /* ISIS_moderator */
  long     _index; /* e.g. 2 index in TRACE list */
  Coords   _position_absolute;
  Coords   _position_relative; /* wrt PREVIOUS */
  Rotation _rotation_absolute;
  Rotation _rotation_relative; /* wrt PREVIOUS */
  int      _rotation_is_identity;
  _class_ISIS_moderator_parameters _parameters;
};
typedef struct _struct_ISIS_moderator _class_ISIS_moderator;
_class_ISIS_moderator _isis_source_var;
_class_ISIS_moderator *_isis_source = &_isis_source_var;
#pragma acc declare create ( _isis_source_var )
#pragma acc declare create ( _isis_source )

/* component defslit1=Slit() [3] DECLARE */
/* Parameter definition for component type 'Slit' */
struct _struct_Slit_parameters {
  /* Component type 'Slit' setting parameters */
  MCNUM xmin;
  MCNUM xmax;
  MCNUM ymin;
  MCNUM ymax;
  MCNUM radius;
  MCNUM xwidth;
  MCNUM yheight;
}; /* _struct_Slit_parameters */
typedef struct _struct_Slit_parameters _class_Slit_parameters;

/* Parameters for component type 'Slit' */
struct _struct_Slit {
  char     _name[256]; /* e.g. defslit1 */
  char     _type[256]; /* Slit */
  long     _index; /* e.g. 2 index in TRACE list */
  Coords   _position_absolute;
  Coords   _position_relative; /* wrt PREVIOUS */
  Rotation _rotation_absolute;
  Rotation _rotation_relative; /* wrt PREVIOUS */
  int      _rotation_is_identity;
  _class_Slit_parameters _parameters;
};
typedef struct _struct_Slit _class_Slit;
_class_Slit _defslit1_var;
_class_Slit *_defslit1 = &_defslit1_var;
#pragma acc declare create ( _defslit1_var )
#pragma acc declare create ( _defslit1 )

_class_Slit _defslit2_var;
_class_Slit *_defslit2 = &_defslit2_var;
#pragma acc declare create ( _defslit2_var )
#pragma acc declare create ( _defslit2 )

_class_Slit _coarseslit1_var;
_class_Slit *_coarseslit1 = &_coarseslit1_var;
#pragma acc declare create ( _coarseslit1_var )
#pragma acc declare create ( _coarseslit1 )

_class_Slit _coarseslit2_var;
_class_Slit *_coarseslit2 = &_coarseslit2_var;
#pragma acc declare create ( _coarseslit2_var )
#pragma acc declare create ( _coarseslit2 )

/* component lmon1=L_monitor() [7] DECLARE */
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
  char     _name[256]; /* e.g. lmon1 */
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
_class_L_monitor _lmon1_var;
_class_L_monitor *_lmon1 = &_lmon1_var;
#pragma acc declare create ( _lmon1_var )
#pragma acc declare create ( _lmon1 )

/* component PSDmon1=PSD_monitor() [8] DECLARE */
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
  char     _name[256]; /* e.g. PSDmon1 */
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
_class_PSD_monitor _PSDmon1_var;
_class_PSD_monitor *_PSDmon1 = &_PSDmon1_var;
#pragma acc declare create ( _PSDmon1_var )
#pragma acc declare create ( _PSDmon1 )

_class_Slit _slit1_var;
_class_Slit *_slit1 = &_slit1_var;
#pragma acc declare create ( _slit1_var )
#pragma acc declare create ( _slit1 )

_class_PSD_monitor _PSDmons1_var;
_class_PSD_monitor *_PSDmons1 = &_PSDmons1_var;
#pragma acc declare create ( _PSDmons1_var )
#pragma acc declare create ( _PSDmons1 )

/* component eguide1=Guide_tapering() [11] DECLARE */
/* Parameter definition for component type 'Guide_tapering' */
struct _struct_Guide_tapering_parameters {
  /* Component type 'Guide_tapering' setting parameters */
  char option[16384];
  MCNUM w1;
  MCNUM h1;
  MCNUM l;
  MCNUM linw;
  MCNUM loutw;
  MCNUM linh;
  MCNUM louth;
  MCNUM R0;
  MCNUM Qcx;
  MCNUM Qcy;
  MCNUM alphax;
  MCNUM alphay;
  MCNUM W;
  MCNUM mx;
  MCNUM my;
  MCNUM segno;
  MCNUM curvature;
  MCNUM curvature_v;
  /* Component type 'Guide_tapering' private parameters */
  /* Component type 'Guide_tapering' DECLARE code stored as structure members */
double *w1c;
double *w2c;
double *ww, *hh;
double *whalf, *hhalf;
double *lwhalf, *lhhalf;
double *h1_in, *h2_out, *w1_in, *w2_out;
double l_seg, h12, h2, w12, w2, a_ell_q, b_ell_q, lbw, lbh;
double mxi ,u1 ,u2 ,div1, p2_para;
double test,Div1;
int seg;
char *fu;
char *pos;
char file_name[1024];
char *ep;
FILE *num;
double rotation_h, rotation_v;
}; /* _struct_Guide_tapering_parameters */
typedef struct _struct_Guide_tapering_parameters _class_Guide_tapering_parameters;

/* Parameters for component type 'Guide_tapering' */
struct _struct_Guide_tapering {
  char     _name[256]; /* e.g. eguide1 */
  char     _type[256]; /* Guide_tapering */
  long     _index; /* e.g. 2 index in TRACE list */
  Coords   _position_absolute;
  Coords   _position_relative; /* wrt PREVIOUS */
  Rotation _rotation_absolute;
  Rotation _rotation_relative; /* wrt PREVIOUS */
  int      _rotation_is_identity;
  _class_Guide_tapering_parameters _parameters;
};
typedef struct _struct_Guide_tapering _class_Guide_tapering;
_class_Guide_tapering _eguide1_var;
_class_Guide_tapering *_eguide1 = &_eguide1_var;
#pragma acc declare create ( _eguide1_var )
#pragma acc declare create ( _eguide1 )

_class_Slit _slit2_var;
_class_Slit *_slit2 = &_slit2_var;
#pragma acc declare create ( _slit2_var )
#pragma acc declare create ( _slit2 )

_class_L_monitor _lmon2_var;
_class_L_monitor *_lmon2 = &_lmon2_var;
#pragma acc declare create ( _lmon2_var )
#pragma acc declare create ( _lmon2 )

/* component samp1=Multilayer_Sample() [14] DECLARE */
/* Parameter definition for component type 'Multilayer_Sample' */
struct _struct_Multilayer_Sample_parameters {
  /* Component type 'Multilayer_Sample' setting parameters */
  MCNUM xwidth;
  MCNUM zlength;
  MCNUM nlayer;
  MCNUM* sldPar;
  MCNUM* dPar;
  MCNUM* sigmaPar;
  MCNUM frac_inc;
  MCNUM ythick;
  MCNUM mu_inc;
  long target_index;
  MCNUM focus_xw;
  MCNUM focus_yh;
  /* Component type 'Multilayer_Sample' private parameters */
  /* Component type 'Multilayer_Sample' DECLARE code stored as structure members */
double xmin, xmax, zmin, zmax;
double tx, ty, tz;
}; /* _struct_Multilayer_Sample_parameters */
typedef struct _struct_Multilayer_Sample_parameters _class_Multilayer_Sample_parameters;

/* Parameters for component type 'Multilayer_Sample' */
struct _struct_Multilayer_Sample {
  char     _name[256]; /* e.g. samp1 */
  char     _type[256]; /* Multilayer_Sample */
  long     _index; /* e.g. 2 index in TRACE list */
  Coords   _position_absolute;
  Coords   _position_relative; /* wrt PREVIOUS */
  Rotation _rotation_absolute;
  Rotation _rotation_relative; /* wrt PREVIOUS */
  int      _rotation_is_identity;
  _class_Multilayer_Sample_parameters _parameters;
};
typedef struct _struct_Multilayer_Sample _class_Multilayer_Sample;
_class_Multilayer_Sample _samp1_var;
_class_Multilayer_Sample *_samp1 = &_samp1_var;
#pragma acc declare create ( _samp1_var )
#pragma acc declare create ( _samp1 )

_class_Slit _slit3_var;
_class_Slit *_slit3 = &_slit3_var;
#pragma acc declare create ( _slit3_var )
#pragma acc declare create ( _slit3 )

_class_Guide_tapering _eguide2_var;
_class_Guide_tapering *_eguide2 = &_eguide2_var;
#pragma acc declare create ( _eguide2_var )
#pragma acc declare create ( _eguide2 )

_class_L_monitor _lmon3_var;
_class_L_monitor *_lmon3 = &_lmon3_var;
#pragma acc declare create ( _lmon3_var )
#pragma acc declare create ( _lmon3 )

_class_PSD_monitor _PSDdet1_var;
_class_PSD_monitor *_PSDdet1 = &_PSDdet1_var;
#pragma acc declare create ( _PSDdet1_var )
#pragma acc declare create ( _PSDdet1 )

int mcNUMCOMP = 18;

/* User declarations from instrument definition. Can define functions. */
  double glen,flen,spos,tend,t15,fp1;
  double linw,loutw,l,u1,u2,div1,b_ell_q,w1,w12,a_ell_q,lbw,w2;

#undef compcurname
#undef compcurtype
#undef compcurindex
/* end of instrument 'ISIS_CRISP' and components DECLARE */

/* *****************************************************************************
* instrument 'ISIS_CRISP' and components INITIALISE
***************************************************************************** */

/* component a1=Arm() SETTING, POSITION/ROTATION */
int _a1_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_a1_setpos] component a1=Arm() SETTING [Arm:0]");
  stracpy(_a1->_name, "a1", 16384);
  stracpy(_a1->_type, "Arm", 16384);
  _a1->_index=1;
  /* component a1=Arm() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(_a1->_rotation_absolute,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_copy(_a1->_rotation_relative, _a1->_rotation_absolute);
    _a1->_rotation_is_identity =  rot_test_identity(_a1->_rotation_relative);
    _a1->_position_absolute = coords_set(
      0, 0, 0);
    tc1 = coords_neg(_a1->_position_absolute);
    _a1->_position_relative = rot_apply(_a1->_rotation_absolute, tc1);
  } /* a1=Arm() AT ROTATED */
  DEBUG_COMPONENT("a1", _a1->_position_absolute, _a1->_rotation_absolute);
  instrument->_position_absolute[1] = _a1->_position_absolute;
  instrument->_position_relative[1] = _a1->_position_relative;
  instrument->counter_N[1]  = instrument->counter_P[1] = instrument->counter_P2[1] = 0;
  instrument->counter_AbsorbProp[1]= 0;
  return(0);
} /* _a1_setpos */

/* component isis_source=ISIS_moderator() SETTING, POSITION/ROTATION */
int _isis_source_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_isis_source_setpos] component isis_source=ISIS_moderator() SETTING [/usr/share/mcstas/3.0-dev/contrib/ISIS_moderator.comp:836]");
  stracpy(_isis_source->_name, "isis_source", 16384);
  stracpy(_isis_source->_type, "ISIS_moderator", 16384);
  _isis_source->_index=2;
  if("crisp" && strlen("crisp"))
    stracpy(_isis_source->_parameters.Face, "crisp" ? "crisp" : "", 16384);
  else 
  _isis_source->_parameters.Face[0]='\0';
  #define Face (_isis_source->_parameters.Face)
  _isis_source->_parameters.Emin = -6.5;
  #define Emin (_isis_source->_parameters.Emin)
  _isis_source->_parameters.Emax = -0.55;
  #define Emax (_isis_source->_parameters.Emax)
  _isis_source->_parameters.dist = 7.2695;
  #define dist (_isis_source->_parameters.dist)
  _isis_source->_parameters.focus_xw = 0.045;
  #define focus_xw (_isis_source->_parameters.focus_xw)
  _isis_source->_parameters.focus_yh = 0.0045;
  #define focus_yh (_isis_source->_parameters.focus_yh)
  _isis_source->_parameters.xwidth = -1;
  #define xwidth (_isis_source->_parameters.xwidth)
  _isis_source->_parameters.yheight = -1;
  #define yheight (_isis_source->_parameters.yheight)
  _isis_source->_parameters.CAngle = 0.0;
  #define CAngle (_isis_source->_parameters.CAngle)
  _isis_source->_parameters.SAC = 1;
  #define SAC (_isis_source->_parameters.SAC)
  _isis_source->_parameters.Lmin = 0;
  #define Lmin (_isis_source->_parameters.Lmin)
  _isis_source->_parameters.Lmax = 0;
  #define Lmax (_isis_source->_parameters.Lmax)
  _isis_source->_parameters.target_index = + 1;
  #define target_index (_isis_source->_parameters.target_index)
  _isis_source->_parameters.verbose = 0;
  #define verbose (_isis_source->_parameters.verbose)

  #define p_in (_isis_source->_parameters.p_in)
  #define Tnpts (_isis_source->_parameters.Tnpts)
  #define scaleSize (_isis_source->_parameters.scaleSize)
  #define angleArea (_isis_source->_parameters.angleArea)
  #define Nsim (_isis_source->_parameters.Nsim)
  #define Ncount (_isis_source->_parameters.Ncount)
  #define TS (_isis_source->_parameters.TS)
  #define rtE0 (_isis_source->_parameters.rtE0)
  #define rtE1 (_isis_source->_parameters.rtE1)
  #define rtmodX (_isis_source->_parameters.rtmodX)
  #define rtmodY (_isis_source->_parameters.rtmodY)
  #define TargetStation (_isis_source->_parameters.TargetStation)
  #define CurrentWeight (_isis_source->_parameters.CurrentWeight)

  #undef Face
  #undef Emin
  #undef Emax
  #undef dist
  #undef focus_xw
  #undef focus_yh
  #undef xwidth
  #undef yheight
  #undef CAngle
  #undef SAC
  #undef Lmin
  #undef Lmax
  #undef target_index
  #undef verbose
  #undef p_in
  #undef Tnpts
  #undef scaleSize
  #undef angleArea
  #undef Nsim
  #undef Ncount
  #undef TS
  #undef rtE0
  #undef rtE1
  #undef rtmodX
  #undef rtmodY
  #undef TargetStation
  #undef CurrentWeight
  /* component isis_source=ISIS_moderator() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(_isis_source->_rotation_absolute,
      (1.5)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_transpose(_a1->_rotation_absolute, tr1);
    rot_mul(_isis_source->_rotation_absolute, tr1, _isis_source->_rotation_relative);
    _isis_source->_rotation_is_identity =  rot_test_identity(_isis_source->_rotation_relative);
    tc1 = coords_set(
      0.0, 0.0, 0.00001);
    rot_transpose(_a1->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _isis_source->_position_absolute = coords_add(_a1->_position_absolute, tc2);
    tc1 = coords_sub(_a1->_position_absolute, _isis_source->_position_absolute);
    _isis_source->_position_relative = rot_apply(_isis_source->_rotation_absolute, tc1);
  } /* isis_source=ISIS_moderator() AT ROTATED */
  DEBUG_COMPONENT("isis_source", _isis_source->_position_absolute, _isis_source->_rotation_absolute);
  instrument->_position_absolute[2] = _isis_source->_position_absolute;
  instrument->_position_relative[2] = _isis_source->_position_relative;
  instrument->counter_N[2]  = instrument->counter_P[2] = instrument->counter_P2[2] = 0;
  instrument->counter_AbsorbProp[2]= 0;
  return(0);
} /* _isis_source_setpos */

/* component defslit1=Slit() SETTING, POSITION/ROTATION */
int _defslit1_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_defslit1_setpos] component defslit1=Slit() SETTING [/usr/share/mcstas/3.0-dev/optics/Slit.comp:47]");
  stracpy(_defslit1->_name, "defslit1", 16384);
  stracpy(_defslit1->_type, "Slit", 16384);
  _defslit1->_index=3;
  _defslit1->_parameters.xmin = -0.0345;
  #define xmin (_defslit1->_parameters.xmin)
  _defslit1->_parameters.xmax = 0.0345;
  #define xmax (_defslit1->_parameters.xmax)
  _defslit1->_parameters.ymin = -0.01731;
  #define ymin (_defslit1->_parameters.ymin)
  _defslit1->_parameters.ymax = 0.01731;
  #define ymax (_defslit1->_parameters.ymax)
  _defslit1->_parameters.radius = 0;
  #define radius (_defslit1->_parameters.radius)
  _defslit1->_parameters.xwidth = 0;
  #define xwidth (_defslit1->_parameters.xwidth)
  _defslit1->_parameters.yheight = 0;
  #define yheight (_defslit1->_parameters.yheight)

  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef radius
  #undef xwidth
  #undef yheight
  /* component defslit1=Slit() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(_defslit1->_rotation_absolute,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_transpose(_isis_source->_rotation_absolute, tr1);
    rot_mul(_defslit1->_rotation_absolute, tr1, _defslit1->_rotation_relative);
    _defslit1->_rotation_is_identity =  rot_test_identity(_defslit1->_rotation_relative);
    _defslit1->_position_absolute = coords_set(
      0.0, -0.09829, 3.7537);
    tc1 = coords_sub(_isis_source->_position_absolute, _defslit1->_position_absolute);
    _defslit1->_position_relative = rot_apply(_defslit1->_rotation_absolute, tc1);
  } /* defslit1=Slit() AT ROTATED */
  DEBUG_COMPONENT("defslit1", _defslit1->_position_absolute, _defslit1->_rotation_absolute);
  instrument->_position_absolute[3] = _defslit1->_position_absolute;
  instrument->_position_relative[3] = _defslit1->_position_relative;
  instrument->counter_N[3]  = instrument->counter_P[3] = instrument->counter_P2[3] = 0;
  instrument->counter_AbsorbProp[3]= 0;
  return(0);
} /* _defslit1_setpos */

/* component defslit2=Slit() SETTING, POSITION/ROTATION */
int _defslit2_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_defslit2_setpos] component defslit2=Slit() SETTING [/usr/share/mcstas/3.0-dev/optics/Slit.comp:47]");
  stracpy(_defslit2->_name, "defslit2", 16384);
  stracpy(_defslit2->_type, "Slit", 16384);
  _defslit2->_index=4;
  _defslit2->_parameters.xmin = -0.0289;
  #define xmin (_defslit2->_parameters.xmin)
  _defslit2->_parameters.xmax = 0.0289;
  #define xmax (_defslit2->_parameters.xmax)
  _defslit2->_parameters.ymin = -0.0123;
  #define ymin (_defslit2->_parameters.ymin)
  _defslit2->_parameters.ymax = 0.0123;
  #define ymax (_defslit2->_parameters.ymax)
  _defslit2->_parameters.radius = 0;
  #define radius (_defslit2->_parameters.radius)
  _defslit2->_parameters.xwidth = 0;
  #define xwidth (_defslit2->_parameters.xwidth)
  _defslit2->_parameters.yheight = 0;
  #define yheight (_defslit2->_parameters.yheight)

  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef radius
  #undef xwidth
  #undef yheight
  /* component defslit2=Slit() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(_defslit2->_rotation_absolute,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_transpose(_defslit1->_rotation_absolute, tr1);
    rot_mul(_defslit2->_rotation_absolute, tr1, _defslit2->_rotation_relative);
    _defslit2->_rotation_is_identity =  rot_test_identity(_defslit2->_rotation_relative);
    _defslit2->_position_absolute = coords_set(
      0.0, -0.15976, 6.101);
    tc1 = coords_sub(_defslit1->_position_absolute, _defslit2->_position_absolute);
    _defslit2->_position_relative = rot_apply(_defslit2->_rotation_absolute, tc1);
  } /* defslit2=Slit() AT ROTATED */
  DEBUG_COMPONENT("defslit2", _defslit2->_position_absolute, _defslit2->_rotation_absolute);
  instrument->_position_absolute[4] = _defslit2->_position_absolute;
  instrument->_position_relative[4] = _defslit2->_position_relative;
  instrument->counter_N[4]  = instrument->counter_P[4] = instrument->counter_P2[4] = 0;
  instrument->counter_AbsorbProp[4]= 0;
  return(0);
} /* _defslit2_setpos */

/* component coarseslit1=Slit() SETTING, POSITION/ROTATION */
int _coarseslit1_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_coarseslit1_setpos] component coarseslit1=Slit() SETTING [/usr/share/mcstas/3.0-dev/optics/Slit.comp:47]");
  stracpy(_coarseslit1->_name, "coarseslit1", 16384);
  stracpy(_coarseslit1->_type, "Slit", 16384);
  _coarseslit1->_index=5;
  _coarseslit1->_parameters.xmin = -0.03;
  #define xmin (_coarseslit1->_parameters.xmin)
  _coarseslit1->_parameters.xmax = 0.03;
  #define xmax (_coarseslit1->_parameters.xmax)
  _coarseslit1->_parameters.ymin = -0.1;
  #define ymin (_coarseslit1->_parameters.ymin)
  _coarseslit1->_parameters.ymax = 0.1;
  #define ymax (_coarseslit1->_parameters.ymax)
  _coarseslit1->_parameters.radius = 0;
  #define radius (_coarseslit1->_parameters.radius)
  _coarseslit1->_parameters.xwidth = 0;
  #define xwidth (_coarseslit1->_parameters.xwidth)
  _coarseslit1->_parameters.yheight = 0;
  #define yheight (_coarseslit1->_parameters.yheight)

  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef radius
  #undef xwidth
  #undef yheight
  /* component coarseslit1=Slit() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(_coarseslit1->_rotation_absolute,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_transpose(_defslit2->_rotation_absolute, tr1);
    rot_mul(_coarseslit1->_rotation_absolute, tr1, _coarseslit1->_rotation_relative);
    _coarseslit1->_rotation_is_identity =  rot_test_identity(_coarseslit1->_rotation_relative);
    _coarseslit1->_position_absolute = coords_set(
      0.0, -0.1757, 6.7096);
    tc1 = coords_sub(_defslit2->_position_absolute, _coarseslit1->_position_absolute);
    _coarseslit1->_position_relative = rot_apply(_coarseslit1->_rotation_absolute, tc1);
  } /* coarseslit1=Slit() AT ROTATED */
  DEBUG_COMPONENT("coarseslit1", _coarseslit1->_position_absolute, _coarseslit1->_rotation_absolute);
  instrument->_position_absolute[5] = _coarseslit1->_position_absolute;
  instrument->_position_relative[5] = _coarseslit1->_position_relative;
  instrument->counter_N[5]  = instrument->counter_P[5] = instrument->counter_P2[5] = 0;
  instrument->counter_AbsorbProp[5]= 0;
  return(0);
} /* _coarseslit1_setpos */

/* component coarseslit2=Slit() SETTING, POSITION/ROTATION */
int _coarseslit2_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_coarseslit2_setpos] component coarseslit2=Slit() SETTING [/usr/share/mcstas/3.0-dev/optics/Slit.comp:47]");
  stracpy(_coarseslit2->_name, "coarseslit2", 16384);
  stracpy(_coarseslit2->_type, "Slit", 16384);
  _coarseslit2->_index=6;
  _coarseslit2->_parameters.xmin = -0.1;
  #define xmin (_coarseslit2->_parameters.xmin)
  _coarseslit2->_parameters.xmax = 0.1;
  #define xmax (_coarseslit2->_parameters.xmax)
  _coarseslit2->_parameters.ymin = -0.01;
  #define ymin (_coarseslit2->_parameters.ymin)
  _coarseslit2->_parameters.ymax = 0.01;
  #define ymax (_coarseslit2->_parameters.ymax)
  _coarseslit2->_parameters.radius = 0;
  #define radius (_coarseslit2->_parameters.radius)
  _coarseslit2->_parameters.xwidth = 0;
  #define xwidth (_coarseslit2->_parameters.xwidth)
  _coarseslit2->_parameters.yheight = 0;
  #define yheight (_coarseslit2->_parameters.yheight)

  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef radius
  #undef xwidth
  #undef yheight
  /* component coarseslit2=Slit() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(_coarseslit2->_rotation_absolute,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_transpose(_coarseslit1->_rotation_absolute, tr1);
    rot_mul(_coarseslit2->_rotation_absolute, tr1, _coarseslit2->_rotation_relative);
    _coarseslit2->_rotation_is_identity =  rot_test_identity(_coarseslit2->_rotation_relative);
    _coarseslit2->_position_absolute = coords_set(
      0.0, -0.17832, 6.8096);
    tc1 = coords_sub(_coarseslit1->_position_absolute, _coarseslit2->_position_absolute);
    _coarseslit2->_position_relative = rot_apply(_coarseslit2->_rotation_absolute, tc1);
  } /* coarseslit2=Slit() AT ROTATED */
  DEBUG_COMPONENT("coarseslit2", _coarseslit2->_position_absolute, _coarseslit2->_rotation_absolute);
  instrument->_position_absolute[6] = _coarseslit2->_position_absolute;
  instrument->_position_relative[6] = _coarseslit2->_position_relative;
  instrument->counter_N[6]  = instrument->counter_P[6] = instrument->counter_P2[6] = 0;
  instrument->counter_AbsorbProp[6]= 0;
  return(0);
} /* _coarseslit2_setpos */

/* component lmon1=L_monitor() SETTING, POSITION/ROTATION */
int _lmon1_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_lmon1_setpos] component lmon1=L_monitor() SETTING [/usr/share/mcstas/3.0-dev/monitors/L_monitor.comp:65]");
  stracpy(_lmon1->_name, "lmon1", 16384);
  stracpy(_lmon1->_type, "L_monitor", 16384);
  _lmon1->_index=7;
  _lmon1->_parameters.nL = 500;
  #define nL (_lmon1->_parameters.nL)
  if("lmon1.dat" && strlen("lmon1.dat"))
    stracpy(_lmon1->_parameters.filename, "lmon1.dat" ? "lmon1.dat" : "", 16384);
  else 
  _lmon1->_parameters.filename[0]='\0';
  #define filename (_lmon1->_parameters.filename)
  _lmon1->_parameters.xmin = -0.06;
  #define xmin (_lmon1->_parameters.xmin)
  _lmon1->_parameters.xmax = 0.06;
  #define xmax (_lmon1->_parameters.xmax)
  _lmon1->_parameters.ymin = -0.04;
  #define ymin (_lmon1->_parameters.ymin)
  _lmon1->_parameters.ymax = 0.04;
  #define ymax (_lmon1->_parameters.ymax)
  _lmon1->_parameters.xwidth = 0;
  #define xwidth (_lmon1->_parameters.xwidth)
  _lmon1->_parameters.yheight = 0;
  #define yheight (_lmon1->_parameters.yheight)
  _lmon1->_parameters.Lmin = 0.0;
  #define Lmin (_lmon1->_parameters.Lmin)
  _lmon1->_parameters.Lmax = 10.0;
  #define Lmax (_lmon1->_parameters.Lmax)
  _lmon1->_parameters.restore_neutron = 0;
  #define restore_neutron (_lmon1->_parameters.restore_neutron)

  #define L_N (_lmon1->_parameters.L_N)
  #define L_p (_lmon1->_parameters.L_p)
  #define L_p2 (_lmon1->_parameters.L_p2)

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
  /* component lmon1=L_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(_lmon1->_rotation_absolute,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_transpose(_coarseslit2->_rotation_absolute, tr1);
    rot_mul(_lmon1->_rotation_absolute, tr1, _lmon1->_rotation_relative);
    _lmon1->_rotation_is_identity =  rot_test_identity(_lmon1->_rotation_relative);
    _lmon1->_position_absolute = coords_set(
      0.0, -0.1833, 7.0);
    tc1 = coords_sub(_coarseslit2->_position_absolute, _lmon1->_position_absolute);
    _lmon1->_position_relative = rot_apply(_lmon1->_rotation_absolute, tc1);
  } /* lmon1=L_monitor() AT ROTATED */
  DEBUG_COMPONENT("lmon1", _lmon1->_position_absolute, _lmon1->_rotation_absolute);
  instrument->_position_absolute[7] = _lmon1->_position_absolute;
  instrument->_position_relative[7] = _lmon1->_position_relative;
  instrument->counter_N[7]  = instrument->counter_P[7] = instrument->counter_P2[7] = 0;
  instrument->counter_AbsorbProp[7]= 0;
  return(0);
} /* _lmon1_setpos */

/* component PSDmon1=PSD_monitor() SETTING, POSITION/ROTATION */
int _PSDmon1_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_PSDmon1_setpos] component PSDmon1=PSD_monitor() SETTING [/usr/share/mcstas/3.0-dev/monitors/PSD_monitor.comp:64]");
  stracpy(_PSDmon1->_name, "PSDmon1", 16384);
  stracpy(_PSDmon1->_type, "PSD_monitor", 16384);
  _PSDmon1->_index=8;
  _PSDmon1->_parameters.nx = 100;
  #define nx (_PSDmon1->_parameters.nx)
  _PSDmon1->_parameters.ny = 100;
  #define ny (_PSDmon1->_parameters.ny)
  if("PSD1.dat" && strlen("PSD1.dat"))
    stracpy(_PSDmon1->_parameters.filename, "PSD1.dat" ? "PSD1.dat" : "", 16384);
  else 
  _PSDmon1->_parameters.filename[0]='\0';
  #define filename (_PSDmon1->_parameters.filename)
  _PSDmon1->_parameters.xmin = -0.05;
  #define xmin (_PSDmon1->_parameters.xmin)
  _PSDmon1->_parameters.xmax = 0.05;
  #define xmax (_PSDmon1->_parameters.xmax)
  _PSDmon1->_parameters.ymin = -0.01;
  #define ymin (_PSDmon1->_parameters.ymin)
  _PSDmon1->_parameters.ymax = 0.01;
  #define ymax (_PSDmon1->_parameters.ymax)
  _PSDmon1->_parameters.xwidth = 0;
  #define xwidth (_PSDmon1->_parameters.xwidth)
  _PSDmon1->_parameters.yheight = 0;
  #define yheight (_PSDmon1->_parameters.yheight)
  _PSDmon1->_parameters.restore_neutron = 0;
  #define restore_neutron (_PSDmon1->_parameters.restore_neutron)

  #define PSD_N (_PSDmon1->_parameters.PSD_N)
  #define PSD_p (_PSDmon1->_parameters.PSD_p)
  #define PSD_p2 (_PSDmon1->_parameters.PSD_p2)

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
  /* component PSDmon1=PSD_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(_PSDmon1->_rotation_absolute,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_transpose(_lmon1->_rotation_absolute, tr1);
    rot_mul(_PSDmon1->_rotation_absolute, tr1, _PSDmon1->_rotation_relative);
    _PSDmon1->_rotation_is_identity =  rot_test_identity(_PSDmon1->_rotation_relative);
    _PSDmon1->_position_absolute = coords_set(
      0.0, -0.1833, 7.01);
    tc1 = coords_sub(_lmon1->_position_absolute, _PSDmon1->_position_absolute);
    _PSDmon1->_position_relative = rot_apply(_PSDmon1->_rotation_absolute, tc1);
  } /* PSDmon1=PSD_monitor() AT ROTATED */
  DEBUG_COMPONENT("PSDmon1", _PSDmon1->_position_absolute, _PSDmon1->_rotation_absolute);
  instrument->_position_absolute[8] = _PSDmon1->_position_absolute;
  instrument->_position_relative[8] = _PSDmon1->_position_relative;
  instrument->counter_N[8]  = instrument->counter_P[8] = instrument->counter_P2[8] = 0;
  instrument->counter_AbsorbProp[8]= 0;
  return(0);
} /* _PSDmon1_setpos */

/* component slit1=Slit() SETTING, POSITION/ROTATION */
int _slit1_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_slit1_setpos] component slit1=Slit() SETTING [/usr/share/mcstas/3.0-dev/optics/Slit.comp:47]");
  stracpy(_slit1->_name, "slit1", 16384);
  stracpy(_slit1->_type, "Slit", 16384);
  _slit1->_index=9;
  _slit1->_parameters.xmin = -0.01;
  #define xmin (_slit1->_parameters.xmin)
  _slit1->_parameters.xmax = 0.01;
  #define xmax (_slit1->_parameters.xmax)
  _slit1->_parameters.ymin = -0.01;
  #define ymin (_slit1->_parameters.ymin)
  _slit1->_parameters.ymax = 0.01;
  #define ymax (_slit1->_parameters.ymax)
  _slit1->_parameters.radius = 0;
  #define radius (_slit1->_parameters.radius)
  _slit1->_parameters.xwidth = 0.04;
  #define xwidth (_slit1->_parameters.xwidth)
  _slit1->_parameters.yheight = 0.004;
  #define yheight (_slit1->_parameters.yheight)

  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef radius
  #undef xwidth
  #undef yheight
  /* component slit1=Slit() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(_slit1->_rotation_absolute,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_transpose(_PSDmon1->_rotation_absolute, tr1);
    rot_mul(_slit1->_rotation_absolute, tr1, _slit1->_rotation_relative);
    _slit1->_rotation_is_identity =  rot_test_identity(_slit1->_rotation_relative);
    _slit1->_position_absolute = coords_set(
      0.0, -0.1904, 7.2695);
    tc1 = coords_sub(_PSDmon1->_position_absolute, _slit1->_position_absolute);
    _slit1->_position_relative = rot_apply(_slit1->_rotation_absolute, tc1);
  } /* slit1=Slit() AT ROTATED */
  DEBUG_COMPONENT("slit1", _slit1->_position_absolute, _slit1->_rotation_absolute);
  instrument->_position_absolute[9] = _slit1->_position_absolute;
  instrument->_position_relative[9] = _slit1->_position_relative;
  instrument->counter_N[9]  = instrument->counter_P[9] = instrument->counter_P2[9] = 0;
  instrument->counter_AbsorbProp[9]= 0;
  return(0);
} /* _slit1_setpos */

/* component PSDmons1=PSD_monitor() SETTING, POSITION/ROTATION */
int _PSDmons1_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_PSDmons1_setpos] component PSDmons1=PSD_monitor() SETTING [/usr/share/mcstas/3.0-dev/monitors/PSD_monitor.comp:64]");
  stracpy(_PSDmons1->_name, "PSDmons1", 16384);
  stracpy(_PSDmons1->_type, "PSD_monitor", 16384);
  _PSDmons1->_index=10;
  _PSDmons1->_parameters.nx = 100;
  #define nx (_PSDmons1->_parameters.nx)
  _PSDmons1->_parameters.ny = 100;
  #define ny (_PSDmons1->_parameters.ny)
  if("PSDs1.dat" && strlen("PSDs1.dat"))
    stracpy(_PSDmons1->_parameters.filename, "PSDs1.dat" ? "PSDs1.dat" : "", 16384);
  else 
  _PSDmons1->_parameters.filename[0]='\0';
  #define filename (_PSDmons1->_parameters.filename)
  _PSDmons1->_parameters.xmin = -0.05;
  #define xmin (_PSDmons1->_parameters.xmin)
  _PSDmons1->_parameters.xmax = 0.05;
  #define xmax (_PSDmons1->_parameters.xmax)
  _PSDmons1->_parameters.ymin = -0.01;
  #define ymin (_PSDmons1->_parameters.ymin)
  _PSDmons1->_parameters.ymax = 0.01;
  #define ymax (_PSDmons1->_parameters.ymax)
  _PSDmons1->_parameters.xwidth = 0;
  #define xwidth (_PSDmons1->_parameters.xwidth)
  _PSDmons1->_parameters.yheight = 0;
  #define yheight (_PSDmons1->_parameters.yheight)
  _PSDmons1->_parameters.restore_neutron = 0;
  #define restore_neutron (_PSDmons1->_parameters.restore_neutron)

  #define PSD_N (_PSDmons1->_parameters.PSD_N)
  #define PSD_p (_PSDmons1->_parameters.PSD_p)
  #define PSD_p2 (_PSDmons1->_parameters.PSD_p2)

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
  /* component PSDmons1=PSD_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(_PSDmons1->_rotation_absolute,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_transpose(_slit1->_rotation_absolute, tr1);
    rot_mul(_PSDmons1->_rotation_absolute, tr1, _PSDmons1->_rotation_relative);
    _PSDmons1->_rotation_is_identity =  rot_test_identity(_PSDmons1->_rotation_relative);
    _PSDmons1->_position_absolute = coords_set(
      0.0, 7.27 * t15, 7.27);
    tc1 = coords_sub(_slit1->_position_absolute, _PSDmons1->_position_absolute);
    _PSDmons1->_position_relative = rot_apply(_PSDmons1->_rotation_absolute, tc1);
  } /* PSDmons1=PSD_monitor() AT ROTATED */
  DEBUG_COMPONENT("PSDmons1", _PSDmons1->_position_absolute, _PSDmons1->_rotation_absolute);
  instrument->_position_absolute[10] = _PSDmons1->_position_absolute;
  instrument->_position_relative[10] = _PSDmons1->_position_relative;
  instrument->counter_N[10]  = instrument->counter_P[10] = instrument->counter_P2[10] = 0;
  instrument->counter_AbsorbProp[10]= 0;
  return(0);
} /* _PSDmons1_setpos */

/* component eguide1=Guide_tapering() SETTING, POSITION/ROTATION */
int _eguide1_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_eguide1_setpos] component eguide1=Guide_tapering() SETTING [/usr/share/mcstas/3.0-dev/optics/Guide_tapering.comp:115]");
  stracpy(_eguide1->_name, "eguide1", 16384);
  stracpy(_eguide1->_type, "Guide_tapering", 16384);
  _eguide1->_index=11;
  if("elliptical" && strlen("elliptical"))
    stracpy(_eguide1->_parameters.option, "elliptical" ? "elliptical" : "", 16384);
  else 
  _eguide1->_parameters.option[0]='\0';
  #define option (_eguide1->_parameters.option)
  _eguide1->_parameters.w1 = instrument->_parameters._w1;
  #define w1 (_eguide1->_parameters.w1)
  _eguide1->_parameters.h1 = 0.1;
  #define h1 (_eguide1->_parameters.h1)
  _eguide1->_parameters.l = instrument->_parameters._glen;
  #define l (_eguide1->_parameters.l)
  _eguide1->_parameters.linw = 10.0;
  #define linw (_eguide1->_parameters.linw)
  _eguide1->_parameters.loutw = instrument->_parameters._flen;
  #define loutw (_eguide1->_parameters.loutw)
  _eguide1->_parameters.linh = 0;
  #define linh (_eguide1->_parameters.linh)
  _eguide1->_parameters.louth = 0;
  #define louth (_eguide1->_parameters.louth)
  _eguide1->_parameters.R0 = 0.99;
  #define R0 (_eguide1->_parameters.R0)
  _eguide1->_parameters.Qcx = 0.021;
  #define Qcx (_eguide1->_parameters.Qcx)
  _eguide1->_parameters.Qcy = 0.021;
  #define Qcy (_eguide1->_parameters.Qcy)
  _eguide1->_parameters.alphax = 6.07;
  #define alphax (_eguide1->_parameters.alphax)
  _eguide1->_parameters.alphay = 6.07;
  #define alphay (_eguide1->_parameters.alphay)
  _eguide1->_parameters.W = 0.003;
  #define W (_eguide1->_parameters.W)
  _eguide1->_parameters.mx = instrument->_parameters._vm;
  #define mx (_eguide1->_parameters.mx)
  _eguide1->_parameters.my = 0;
  #define my (_eguide1->_parameters.my)
  _eguide1->_parameters.segno = 10;
  #define segno (_eguide1->_parameters.segno)
  _eguide1->_parameters.curvature = 0;
  #define curvature (_eguide1->_parameters.curvature)
  _eguide1->_parameters.curvature_v = 0;
  #define curvature_v (_eguide1->_parameters.curvature_v)

  #define w1c (_eguide1->_parameters.w1c)
  #define w2c (_eguide1->_parameters.w2c)
  #define ww (_eguide1->_parameters.ww)
  #define hh (_eguide1->_parameters.hh)
  #define whalf (_eguide1->_parameters.whalf)
  #define hhalf (_eguide1->_parameters.hhalf)
  #define lwhalf (_eguide1->_parameters.lwhalf)
  #define lhhalf (_eguide1->_parameters.lhhalf)
  #define h1_in (_eguide1->_parameters.h1_in)
  #define h2_out (_eguide1->_parameters.h2_out)
  #define w1_in (_eguide1->_parameters.w1_in)
  #define w2_out (_eguide1->_parameters.w2_out)
  #define l_seg (_eguide1->_parameters.l_seg)
  #define seg (_eguide1->_parameters.seg)
  #define h12 (_eguide1->_parameters.h12)
  #define h2 (_eguide1->_parameters.h2)
  #define w12 (_eguide1->_parameters.w12)
  #define w2 (_eguide1->_parameters.w2)
  #define a_ell_q (_eguide1->_parameters.a_ell_q)
  #define b_ell_q (_eguide1->_parameters.b_ell_q)
  #define lbw (_eguide1->_parameters.lbw)
  #define lbh (_eguide1->_parameters.lbh)
  #define mxi (_eguide1->_parameters.mxi)
  #define u1 (_eguide1->_parameters.u1)
  #define u2 (_eguide1->_parameters.u2)
  #define div1 (_eguide1->_parameters.div1)
  #define p2_para (_eguide1->_parameters.p2_para)
  #define test (_eguide1->_parameters.test)
  #define Div1 (_eguide1->_parameters.Div1)
  #define seg (_eguide1->_parameters.seg)
  #define fu (_eguide1->_parameters.fu)
  #define pos (_eguide1->_parameters.pos)
  #define file_name (_eguide1->_parameters.file_name)
  #define ep (_eguide1->_parameters.ep)
  #define num (_eguide1->_parameters.num)
  #define rotation_h (_eguide1->_parameters.rotation_h)
  #define rotation_v (_eguide1->_parameters.rotation_v)

  #undef option
  #undef w1
  #undef h1
  #undef l
  #undef linw
  #undef loutw
  #undef linh
  #undef louth
  #undef R0
  #undef Qcx
  #undef Qcy
  #undef alphax
  #undef alphay
  #undef W
  #undef mx
  #undef my
  #undef segno
  #undef curvature
  #undef curvature_v
  #undef w1c
  #undef w2c
  #undef ww
  #undef hh
  #undef whalf
  #undef hhalf
  #undef lwhalf
  #undef lhhalf
  #undef h1_in
  #undef h2_out
  #undef w1_in
  #undef w2_out
  #undef l_seg
  #undef seg
  #undef h12
  #undef h2
  #undef w12
  #undef w2
  #undef a_ell_q
  #undef b_ell_q
  #undef lbw
  #undef lbh
  #undef mxi
  #undef u1
  #undef u2
  #undef div1
  #undef p2_para
  #undef test
  #undef Div1
  #undef seg
  #undef fu
  #undef pos
  #undef file_name
  #undef ep
  #undef num
  #undef rotation_h
  #undef rotation_v
  /* component eguide1=Guide_tapering() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(_eguide1->_rotation_absolute,
      (1.5)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_transpose(_PSDmons1->_rotation_absolute, tr1);
    rot_mul(_eguide1->_rotation_absolute, tr1, _eguide1->_rotation_relative);
    _eguide1->_rotation_is_identity =  rot_test_identity(_eguide1->_rotation_relative);
    _eguide1->_position_absolute = coords_set(
      0.0, ( spos - instrument->_parameters._flen - instrument->_parameters._glen ) * t15, spos - instrument->_parameters._flen - instrument->_parameters._glen);
    tc1 = coords_sub(_PSDmons1->_position_absolute, _eguide1->_position_absolute);
    _eguide1->_position_relative = rot_apply(_eguide1->_rotation_absolute, tc1);
  } /* eguide1=Guide_tapering() AT ROTATED */
  DEBUG_COMPONENT("eguide1", _eguide1->_position_absolute, _eguide1->_rotation_absolute);
  instrument->_position_absolute[11] = _eguide1->_position_absolute;
  instrument->_position_relative[11] = _eguide1->_position_relative;
  instrument->counter_N[11]  = instrument->counter_P[11] = instrument->counter_P2[11] = 0;
  instrument->counter_AbsorbProp[11]= 0;
  return(0);
} /* _eguide1_setpos */

/* component slit2=Slit() SETTING, POSITION/ROTATION */
int _slit2_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_slit2_setpos] component slit2=Slit() SETTING [/usr/share/mcstas/3.0-dev/optics/Slit.comp:47]");
  stracpy(_slit2->_name, "slit2", 16384);
  stracpy(_slit2->_type, "Slit", 16384);
  _slit2->_index=12;
  _slit2->_parameters.xmin = -0.01;
  #define xmin (_slit2->_parameters.xmin)
  _slit2->_parameters.xmax = 0.01;
  #define xmax (_slit2->_parameters.xmax)
  _slit2->_parameters.ymin = -0.01;
  #define ymin (_slit2->_parameters.ymin)
  _slit2->_parameters.ymax = 0.01;
  #define ymax (_slit2->_parameters.ymax)
  _slit2->_parameters.radius = 0;
  #define radius (_slit2->_parameters.radius)
  _slit2->_parameters.xwidth = tend;
  #define xwidth (_slit2->_parameters.xwidth)
  _slit2->_parameters.yheight = 0.0025;
  #define yheight (_slit2->_parameters.yheight)

  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef radius
  #undef xwidth
  #undef yheight
  /* component slit2=Slit() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(_slit2->_rotation_absolute,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_transpose(_eguide1->_rotation_absolute, tr1);
    rot_mul(_slit2->_rotation_absolute, tr1, _slit2->_rotation_relative);
    _slit2->_rotation_is_identity =  rot_test_identity(_slit2->_rotation_relative);
    _slit2->_position_absolute = coords_set(
      0.0, -0.2581, 9.8555);
    tc1 = coords_sub(_eguide1->_position_absolute, _slit2->_position_absolute);
    _slit2->_position_relative = rot_apply(_slit2->_rotation_absolute, tc1);
  } /* slit2=Slit() AT ROTATED */
  DEBUG_COMPONENT("slit2", _slit2->_position_absolute, _slit2->_rotation_absolute);
  instrument->_position_absolute[12] = _slit2->_position_absolute;
  instrument->_position_relative[12] = _slit2->_position_relative;
  instrument->counter_N[12]  = instrument->counter_P[12] = instrument->counter_P2[12] = 0;
  instrument->counter_AbsorbProp[12]= 0;
  return(0);
} /* _slit2_setpos */

/* component lmon2=L_monitor() SETTING, POSITION/ROTATION */
int _lmon2_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_lmon2_setpos] component lmon2=L_monitor() SETTING [/usr/share/mcstas/3.0-dev/monitors/L_monitor.comp:65]");
  stracpy(_lmon2->_name, "lmon2", 16384);
  stracpy(_lmon2->_type, "L_monitor", 16384);
  _lmon2->_index=13;
  _lmon2->_parameters.nL = 500;
  #define nL (_lmon2->_parameters.nL)
  if("lmon2.dat" && strlen("lmon2.dat"))
    stracpy(_lmon2->_parameters.filename, "lmon2.dat" ? "lmon2.dat" : "", 16384);
  else 
  _lmon2->_parameters.filename[0]='\0';
  #define filename (_lmon2->_parameters.filename)
  _lmon2->_parameters.xmin = -0.06;
  #define xmin (_lmon2->_parameters.xmin)
  _lmon2->_parameters.xmax = 0.06;
  #define xmax (_lmon2->_parameters.xmax)
  _lmon2->_parameters.ymin = -0.04;
  #define ymin (_lmon2->_parameters.ymin)
  _lmon2->_parameters.ymax = 0.04;
  #define ymax (_lmon2->_parameters.ymax)
  _lmon2->_parameters.xwidth = 0;
  #define xwidth (_lmon2->_parameters.xwidth)
  _lmon2->_parameters.yheight = 0;
  #define yheight (_lmon2->_parameters.yheight)
  _lmon2->_parameters.Lmin = 0.0;
  #define Lmin (_lmon2->_parameters.Lmin)
  _lmon2->_parameters.Lmax = 10.0;
  #define Lmax (_lmon2->_parameters.Lmax)
  _lmon2->_parameters.restore_neutron = 0;
  #define restore_neutron (_lmon2->_parameters.restore_neutron)

  #define L_N (_lmon2->_parameters.L_N)
  #define L_p (_lmon2->_parameters.L_p)
  #define L_p2 (_lmon2->_parameters.L_p2)

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
  /* component lmon2=L_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(_lmon2->_rotation_absolute,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_transpose(_slit2->_rotation_absolute, tr1);
    rot_mul(_lmon2->_rotation_absolute, tr1, _lmon2->_rotation_relative);
    _lmon2->_rotation_is_identity =  rot_test_identity(_lmon2->_rotation_relative);
    _lmon2->_position_absolute = coords_set(
      0.0, -0.2608, 9.96);
    tc1 = coords_sub(_slit2->_position_absolute, _lmon2->_position_absolute);
    _lmon2->_position_relative = rot_apply(_lmon2->_rotation_absolute, tc1);
  } /* lmon2=L_monitor() AT ROTATED */
  DEBUG_COMPONENT("lmon2", _lmon2->_position_absolute, _lmon2->_rotation_absolute);
  instrument->_position_absolute[13] = _lmon2->_position_absolute;
  instrument->_position_relative[13] = _lmon2->_position_relative;
  instrument->counter_N[13]  = instrument->counter_P[13] = instrument->counter_P2[13] = 0;
  instrument->counter_AbsorbProp[13]= 0;
  return(0);
} /* _lmon2_setpos */

/* component samp1=Multilayer_Sample() SETTING, POSITION/ROTATION */
int _samp1_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_samp1_setpos] component samp1=Multilayer_Sample() SETTING [/usr/share/mcstas/3.0-dev/contrib/Multilayer_Sample.comp:77]");
  stracpy(_samp1->_name, "samp1", 16384);
  stracpy(_samp1->_type, "Multilayer_Sample", 16384);
  _samp1->_index=14;
  _samp1->_parameters.xwidth = 0.05;
  #define xwidth (_samp1->_parameters.xwidth)
  _samp1->_parameters.zlength = 0.15;
  #define zlength (_samp1->_parameters.zlength)
  _samp1->_parameters.nlayer = 0;
  #define nlayer (_samp1->_parameters.nlayer)
  _samp1->_parameters.sldPar = (double []){ 0.0 , 6.35e-6 }; // static pointer allocation
  #define sldPar (_samp1->_parameters.sldPar)
  _samp1->_parameters.dPar = (double []){ 0.0 }; // static pointer allocation
  #define dPar (_samp1->_parameters.dPar)
  _samp1->_parameters.sigmaPar = (double []){ 5.0 }; // static pointer allocation
  #define sigmaPar (_samp1->_parameters.sigmaPar)
  _samp1->_parameters.frac_inc = instrument->_parameters._FRAC;
  #define frac_inc (_samp1->_parameters.frac_inc)
  _samp1->_parameters.ythick = 0.01;
  #define ythick (_samp1->_parameters.ythick)
  _samp1->_parameters.mu_inc = 0.138;
  #define mu_inc (_samp1->_parameters.mu_inc)
  _samp1->_parameters.target_index = 1;
  #define target_index (_samp1->_parameters.target_index)
  _samp1->_parameters.focus_xw = 2 * tend;
  #define focus_xw (_samp1->_parameters.focus_xw)
  _samp1->_parameters.focus_yh = 0.01;
  #define focus_yh (_samp1->_parameters.focus_yh)

  #define tx (_samp1->_parameters.tx)
  #define ty (_samp1->_parameters.ty)
  #define tz (_samp1->_parameters.tz)
  #define xmin (_samp1->_parameters.xmin)
  #define xmax (_samp1->_parameters.xmax)
  #define zmin (_samp1->_parameters.zmin)
  #define zmax (_samp1->_parameters.zmax)

  #undef xwidth
  #undef zlength
  #undef nlayer
  #undef sldPar
  #undef dPar
  #undef sigmaPar
  #undef frac_inc
  #undef ythick
  #undef mu_inc
  #undef target_index
  #undef focus_xw
  #undef focus_yh
  #undef tx
  #undef ty
  #undef tz
  #undef xmin
  #undef xmax
  #undef zmin
  #undef zmax
  /* component samp1=Multilayer_Sample() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(_samp1->_rotation_absolute,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_transpose(_lmon2->_rotation_absolute, tr1);
    rot_mul(_samp1->_rotation_absolute, tr1, _samp1->_rotation_relative);
    _samp1->_rotation_is_identity =  rot_test_identity(_samp1->_rotation_relative);
    _samp1->_position_absolute = coords_set(
      0.0, spos * t15, spos);
    tc1 = coords_sub(_lmon2->_position_absolute, _samp1->_position_absolute);
    _samp1->_position_relative = rot_apply(_samp1->_rotation_absolute, tc1);
  } /* samp1=Multilayer_Sample() AT ROTATED */
  DEBUG_COMPONENT("samp1", _samp1->_position_absolute, _samp1->_rotation_absolute);
  instrument->_position_absolute[14] = _samp1->_position_absolute;
  instrument->_position_relative[14] = _samp1->_position_relative;
  instrument->counter_N[14]  = instrument->counter_P[14] = instrument->counter_P2[14] = 0;
  instrument->counter_AbsorbProp[14]= 0;
  return(0);
} /* _samp1_setpos */

/* component slit3=Slit() SETTING, POSITION/ROTATION */
int _slit3_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_slit3_setpos] component slit3=Slit() SETTING [/usr/share/mcstas/3.0-dev/optics/Slit.comp:47]");
  stracpy(_slit3->_name, "slit3", 16384);
  stracpy(_slit3->_type, "Slit", 16384);
  _slit3->_index=15;
  _slit3->_parameters.xmin = -0.01;
  #define xmin (_slit3->_parameters.xmin)
  _slit3->_parameters.xmax = 0.01;
  #define xmax (_slit3->_parameters.xmax)
  _slit3->_parameters.ymin = -0.01;
  #define ymin (_slit3->_parameters.ymin)
  _slit3->_parameters.ymax = 0.01;
  #define ymax (_slit3->_parameters.ymax)
  _slit3->_parameters.radius = 0;
  #define radius (_slit3->_parameters.radius)
  _slit3->_parameters.xwidth = tend;
  #define xwidth (_slit3->_parameters.xwidth)
  _slit3->_parameters.yheight = 0.003;
  #define yheight (_slit3->_parameters.yheight)

  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef radius
  #undef xwidth
  #undef yheight
  /* component slit3=Slit() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(_slit3->_rotation_absolute,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_transpose(_samp1->_rotation_absolute, tr1);
    rot_mul(_slit3->_rotation_absolute, tr1, _slit3->_rotation_relative);
    _slit3->_rotation_is_identity =  rot_test_identity(_slit3->_rotation_relative);
    _slit3->_position_absolute = coords_set(
      0.0, ( spos + instrument->_parameters._flen -0.01 ) * t15 -2.0 * ( instrument->_parameters._flen -0.01 ) * t15, spos + instrument->_parameters._flen -0.01);
    tc1 = coords_sub(_samp1->_position_absolute, _slit3->_position_absolute);
    _slit3->_position_relative = rot_apply(_slit3->_rotation_absolute, tc1);
  } /* slit3=Slit() AT ROTATED */
  DEBUG_COMPONENT("slit3", _slit3->_position_absolute, _slit3->_rotation_absolute);
  instrument->_position_absolute[15] = _slit3->_position_absolute;
  instrument->_position_relative[15] = _slit3->_position_relative;
  instrument->counter_N[15]  = instrument->counter_P[15] = instrument->counter_P2[15] = 0;
  instrument->counter_AbsorbProp[15]= 0;
  return(0);
} /* _slit3_setpos */

/* component eguide2=Guide_tapering() SETTING, POSITION/ROTATION */
int _eguide2_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_eguide2_setpos] component eguide2=Guide_tapering() SETTING [/usr/share/mcstas/3.0-dev/optics/Guide_tapering.comp:115]");
  stracpy(_eguide2->_name, "eguide2", 16384);
  stracpy(_eguide2->_type, "Guide_tapering", 16384);
  _eguide2->_index=16;
  if("elliptical" && strlen("elliptical"))
    stracpy(_eguide2->_parameters.option, "elliptical" ? "elliptical" : "", 16384);
  else 
  _eguide2->_parameters.option[0]='\0';
  #define option (_eguide2->_parameters.option)
  _eguide2->_parameters.w1 = tend;
  #define w1 (_eguide2->_parameters.w1)
  _eguide2->_parameters.h1 = 0.1;
  #define h1 (_eguide2->_parameters.h1)
  _eguide2->_parameters.l = instrument->_parameters._glen;
  #define l (_eguide2->_parameters.l)
  _eguide2->_parameters.linw = instrument->_parameters._flen;
  #define linw (_eguide2->_parameters.linw)
  _eguide2->_parameters.loutw = 10.0;
  #define loutw (_eguide2->_parameters.loutw)
  _eguide2->_parameters.linh = 0;
  #define linh (_eguide2->_parameters.linh)
  _eguide2->_parameters.louth = 0;
  #define louth (_eguide2->_parameters.louth)
  _eguide2->_parameters.R0 = 0.99;
  #define R0 (_eguide2->_parameters.R0)
  _eguide2->_parameters.Qcx = 0.021;
  #define Qcx (_eguide2->_parameters.Qcx)
  _eguide2->_parameters.Qcy = 0.021;
  #define Qcy (_eguide2->_parameters.Qcy)
  _eguide2->_parameters.alphax = 6.07;
  #define alphax (_eguide2->_parameters.alphax)
  _eguide2->_parameters.alphay = 6.07;
  #define alphay (_eguide2->_parameters.alphay)
  _eguide2->_parameters.W = 0.003;
  #define W (_eguide2->_parameters.W)
  _eguide2->_parameters.mx = instrument->_parameters._vm;
  #define mx (_eguide2->_parameters.mx)
  _eguide2->_parameters.my = 0;
  #define my (_eguide2->_parameters.my)
  _eguide2->_parameters.segno = 10;
  #define segno (_eguide2->_parameters.segno)
  _eguide2->_parameters.curvature = 0;
  #define curvature (_eguide2->_parameters.curvature)
  _eguide2->_parameters.curvature_v = 0;
  #define curvature_v (_eguide2->_parameters.curvature_v)

  #define w1c (_eguide2->_parameters.w1c)
  #define w2c (_eguide2->_parameters.w2c)
  #define ww (_eguide2->_parameters.ww)
  #define hh (_eguide2->_parameters.hh)
  #define whalf (_eguide2->_parameters.whalf)
  #define hhalf (_eguide2->_parameters.hhalf)
  #define lwhalf (_eguide2->_parameters.lwhalf)
  #define lhhalf (_eguide2->_parameters.lhhalf)
  #define h1_in (_eguide2->_parameters.h1_in)
  #define h2_out (_eguide2->_parameters.h2_out)
  #define w1_in (_eguide2->_parameters.w1_in)
  #define w2_out (_eguide2->_parameters.w2_out)
  #define l_seg (_eguide2->_parameters.l_seg)
  #define seg (_eguide2->_parameters.seg)
  #define h12 (_eguide2->_parameters.h12)
  #define h2 (_eguide2->_parameters.h2)
  #define w12 (_eguide2->_parameters.w12)
  #define w2 (_eguide2->_parameters.w2)
  #define a_ell_q (_eguide2->_parameters.a_ell_q)
  #define b_ell_q (_eguide2->_parameters.b_ell_q)
  #define lbw (_eguide2->_parameters.lbw)
  #define lbh (_eguide2->_parameters.lbh)
  #define mxi (_eguide2->_parameters.mxi)
  #define u1 (_eguide2->_parameters.u1)
  #define u2 (_eguide2->_parameters.u2)
  #define div1 (_eguide2->_parameters.div1)
  #define p2_para (_eguide2->_parameters.p2_para)
  #define test (_eguide2->_parameters.test)
  #define Div1 (_eguide2->_parameters.Div1)
  #define seg (_eguide2->_parameters.seg)
  #define fu (_eguide2->_parameters.fu)
  #define pos (_eguide2->_parameters.pos)
  #define file_name (_eguide2->_parameters.file_name)
  #define ep (_eguide2->_parameters.ep)
  #define num (_eguide2->_parameters.num)
  #define rotation_h (_eguide2->_parameters.rotation_h)
  #define rotation_v (_eguide2->_parameters.rotation_v)

  #undef option
  #undef w1
  #undef h1
  #undef l
  #undef linw
  #undef loutw
  #undef linh
  #undef louth
  #undef R0
  #undef Qcx
  #undef Qcy
  #undef alphax
  #undef alphay
  #undef W
  #undef mx
  #undef my
  #undef segno
  #undef curvature
  #undef curvature_v
  #undef w1c
  #undef w2c
  #undef ww
  #undef hh
  #undef whalf
  #undef hhalf
  #undef lwhalf
  #undef lhhalf
  #undef h1_in
  #undef h2_out
  #undef w1_in
  #undef w2_out
  #undef l_seg
  #undef seg
  #undef h12
  #undef h2
  #undef w12
  #undef w2
  #undef a_ell_q
  #undef b_ell_q
  #undef lbw
  #undef lbh
  #undef mxi
  #undef u1
  #undef u2
  #undef div1
  #undef p2_para
  #undef test
  #undef Div1
  #undef seg
  #undef fu
  #undef pos
  #undef file_name
  #undef ep
  #undef num
  #undef rotation_h
  #undef rotation_v
  /* component eguide2=Guide_tapering() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(_eguide2->_rotation_absolute,
      (-1.5)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_transpose(_slit3->_rotation_absolute, tr1);
    rot_mul(_eguide2->_rotation_absolute, tr1, _eguide2->_rotation_relative);
    _eguide2->_rotation_is_identity =  rot_test_identity(_eguide2->_rotation_relative);
    _eguide2->_position_absolute = coords_set(
      0.0, ( spos + instrument->_parameters._flen ) * t15 -2.0 * instrument->_parameters._flen * t15, spos + instrument->_parameters._flen);
    tc1 = coords_sub(_slit3->_position_absolute, _eguide2->_position_absolute);
    _eguide2->_position_relative = rot_apply(_eguide2->_rotation_absolute, tc1);
  } /* eguide2=Guide_tapering() AT ROTATED */
  DEBUG_COMPONENT("eguide2", _eguide2->_position_absolute, _eguide2->_rotation_absolute);
  instrument->_position_absolute[16] = _eguide2->_position_absolute;
  instrument->_position_relative[16] = _eguide2->_position_relative;
  instrument->counter_N[16]  = instrument->counter_P[16] = instrument->counter_P2[16] = 0;
  instrument->counter_AbsorbProp[16]= 0;
  return(0);
} /* _eguide2_setpos */

/* component lmon3=L_monitor() SETTING, POSITION/ROTATION */
int _lmon3_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_lmon3_setpos] component lmon3=L_monitor() SETTING [/usr/share/mcstas/3.0-dev/monitors/L_monitor.comp:65]");
  stracpy(_lmon3->_name, "lmon3", 16384);
  stracpy(_lmon3->_type, "L_monitor", 16384);
  _lmon3->_index=17;
  _lmon3->_parameters.nL = 500;
  #define nL (_lmon3->_parameters.nL)
  if("lmon3.dat" && strlen("lmon3.dat"))
    stracpy(_lmon3->_parameters.filename, "lmon3.dat" ? "lmon3.dat" : "", 16384);
  else 
  _lmon3->_parameters.filename[0]='\0';
  #define filename (_lmon3->_parameters.filename)
  _lmon3->_parameters.xmin = -0.05;
  #define xmin (_lmon3->_parameters.xmin)
  _lmon3->_parameters.xmax = 0.05;
  #define xmax (_lmon3->_parameters.xmax)
  _lmon3->_parameters.ymin = -0.05;
  #define ymin (_lmon3->_parameters.ymin)
  _lmon3->_parameters.ymax = 0.05;
  #define ymax (_lmon3->_parameters.ymax)
  _lmon3->_parameters.xwidth = 0;
  #define xwidth (_lmon3->_parameters.xwidth)
  _lmon3->_parameters.yheight = 0;
  #define yheight (_lmon3->_parameters.yheight)
  _lmon3->_parameters.Lmin = 0.0;
  #define Lmin (_lmon3->_parameters.Lmin)
  _lmon3->_parameters.Lmax = 10.0;
  #define Lmax (_lmon3->_parameters.Lmax)
  _lmon3->_parameters.restore_neutron = 0;
  #define restore_neutron (_lmon3->_parameters.restore_neutron)

  #define L_N (_lmon3->_parameters.L_N)
  #define L_p (_lmon3->_parameters.L_p)
  #define L_p2 (_lmon3->_parameters.L_p2)

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
  /* component lmon3=L_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(_lmon3->_rotation_absolute,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_transpose(_eguide2->_rotation_absolute, tr1);
    rot_mul(_lmon3->_rotation_absolute, tr1, _lmon3->_rotation_relative);
    _lmon3->_rotation_is_identity =  rot_test_identity(_lmon3->_rotation_relative);
    _lmon3->_position_absolute = coords_set(
      0.0, 12.11 * t15 -2.0 * ( 12.11 - spos ) * t15, 12.11);
    tc1 = coords_sub(_eguide2->_position_absolute, _lmon3->_position_absolute);
    _lmon3->_position_relative = rot_apply(_lmon3->_rotation_absolute, tc1);
  } /* lmon3=L_monitor() AT ROTATED */
  DEBUG_COMPONENT("lmon3", _lmon3->_position_absolute, _lmon3->_rotation_absolute);
  instrument->_position_absolute[17] = _lmon3->_position_absolute;
  instrument->_position_relative[17] = _lmon3->_position_relative;
  instrument->counter_N[17]  = instrument->counter_P[17] = instrument->counter_P2[17] = 0;
  instrument->counter_AbsorbProp[17]= 0;
  return(0);
} /* _lmon3_setpos */

/* component PSDdet1=PSD_monitor() SETTING, POSITION/ROTATION */
int _PSDdet1_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_PSDdet1_setpos] component PSDdet1=PSD_monitor() SETTING [/usr/share/mcstas/3.0-dev/monitors/PSD_monitor.comp:64]");
  stracpy(_PSDdet1->_name, "PSDdet1", 16384);
  stracpy(_PSDdet1->_type, "PSD_monitor", 16384);
  _PSDdet1->_index=18;
  _PSDdet1->_parameters.nx = 100;
  #define nx (_PSDdet1->_parameters.nx)
  _PSDdet1->_parameters.ny = 100;
  #define ny (_PSDdet1->_parameters.ny)
  if("PSD3.dat" && strlen("PSD3.dat"))
    stracpy(_PSDdet1->_parameters.filename, "PSD3.dat" ? "PSD3.dat" : "", 16384);
  else 
  _PSDdet1->_parameters.filename[0]='\0';
  #define filename (_PSDdet1->_parameters.filename)
  _PSDdet1->_parameters.xmin = -0.05;
  #define xmin (_PSDdet1->_parameters.xmin)
  _PSDdet1->_parameters.xmax = 0.05;
  #define xmax (_PSDdet1->_parameters.xmax)
  _PSDdet1->_parameters.ymin = -0.05;
  #define ymin (_PSDdet1->_parameters.ymin)
  _PSDdet1->_parameters.ymax = 0.05;
  #define ymax (_PSDdet1->_parameters.ymax)
  _PSDdet1->_parameters.xwidth = 0;
  #define xwidth (_PSDdet1->_parameters.xwidth)
  _PSDdet1->_parameters.yheight = 0;
  #define yheight (_PSDdet1->_parameters.yheight)
  _PSDdet1->_parameters.restore_neutron = 0;
  #define restore_neutron (_PSDdet1->_parameters.restore_neutron)

  #define PSD_N (_PSDdet1->_parameters.PSD_N)
  #define PSD_p (_PSDdet1->_parameters.PSD_p)
  #define PSD_p2 (_PSDdet1->_parameters.PSD_p2)

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
  /* component PSDdet1=PSD_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(_PSDdet1->_rotation_absolute,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_transpose(_lmon3->_rotation_absolute, tr1);
    rot_mul(_PSDdet1->_rotation_absolute, tr1, _PSDdet1->_rotation_relative);
    _PSDdet1->_rotation_is_identity =  rot_test_identity(_PSDdet1->_rotation_relative);
    _PSDdet1->_position_absolute = coords_set(
      0.0, 12.111 * t15 -2.0 * ( 12.111 - spos ) * t15, 12.111);
    tc1 = coords_sub(_lmon3->_position_absolute, _PSDdet1->_position_absolute);
    _PSDdet1->_position_relative = rot_apply(_PSDdet1->_rotation_absolute, tc1);
  } /* PSDdet1=PSD_monitor() AT ROTATED */
  DEBUG_COMPONENT("PSDdet1", _PSDdet1->_position_absolute, _PSDdet1->_rotation_absolute);
  instrument->_position_absolute[18] = _PSDdet1->_position_absolute;
  instrument->_position_relative[18] = _PSDdet1->_position_relative;
  instrument->counter_N[18]  = instrument->counter_P[18] = instrument->counter_P2[18] = 0;
  instrument->counter_AbsorbProp[18]= 0;
  return(0);
} /* _PSDdet1_setpos */

_class_ISIS_moderator *class_ISIS_moderator_init(_class_ISIS_moderator *_comp
) {
  #define Face (_comp->_parameters.Face)
  #define Emin (_comp->_parameters.Emin)
  #define Emax (_comp->_parameters.Emax)
  #define dist (_comp->_parameters.dist)
  #define focus_xw (_comp->_parameters.focus_xw)
  #define focus_yh (_comp->_parameters.focus_yh)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define CAngle (_comp->_parameters.CAngle)
  #define SAC (_comp->_parameters.SAC)
  #define Lmin (_comp->_parameters.Lmin)
  #define Lmax (_comp->_parameters.Lmax)
  #define target_index (_comp->_parameters.target_index)
  #define verbose (_comp->_parameters.verbose)
  #define p_in (_comp->_parameters.p_in)
  #define Tnpts (_comp->_parameters.Tnpts)
  #define scaleSize (_comp->_parameters.scaleSize)
  #define angleArea (_comp->_parameters.angleArea)
  #define Nsim (_comp->_parameters.Nsim)
  #define Ncount (_comp->_parameters.Ncount)
  #define TS (_comp->_parameters.TS)
  #define rtE0 (_comp->_parameters.rtE0)
  #define rtE1 (_comp->_parameters.rtE1)
  #define rtmodX (_comp->_parameters.rtmodX)
  #define rtmodY (_comp->_parameters.rtmodY)
  #define TargetStation (_comp->_parameters.TargetStation)
  #define CurrentWeight (_comp->_parameters.CurrentWeight)
  /* READ IN THE ENERGY FILE */

  char fname[256];   /* Variables */
  FILE* TFile;
  double tmp;
  char lowerFace[255];
  int Bcnt;
  int i;
  struct BeamLine
    {
      char Name[50];
      double Xsize;
      double Ysize;
    } Olist[50];
    
  if (target_index && !dist)
  {
    Coords ToTarget;
    double tx,ty,tz;
    ToTarget = coords_sub(POS_A_COMP_INDEX(INDEX_CURRENT_COMP+target_index),POS_A_CURRENT_COMP);
    ToTarget = rot_apply(ROT_A_CURRENT_COMP, ToTarget);
    coords_get(ToTarget, &tx, &ty, &tz);
    dist=sqrt(tx*tx+ty*ty+tz*tz);
  }
  
  Nsim=(double)mcget_ncount();
  Bcnt=0;
  // CH4 face 1 (north)
  strcpy(Olist[Bcnt].Name,"mari"); Olist[Bcnt].Xsize=0.12; Olist[Bcnt].Ysize=0.115; Bcnt++;
  strcpy(Olist[Bcnt].Name,"gem"); Olist[Bcnt].Xsize=0.12; Olist[Bcnt].Ysize=0.115; Bcnt++;
  strcpy(Olist[Bcnt].Name,"hrpd"); Olist[Bcnt].Xsize=0.12; Olist[Bcnt].Ysize=0.115; Bcnt++;
  strcpy(Olist[Bcnt].Name,"pearl"); Olist[Bcnt].Xsize=0.12; Olist[Bcnt].Ysize=0.115; Bcnt++;
  // CH4 face 2 (south)
  strcpy(Olist[Bcnt].Name,"sandals"); Olist[Bcnt].Xsize=0.12; Olist[Bcnt].Ysize=0.115; Bcnt++;
  strcpy(Olist[Bcnt].Name,"prisma"); Olist[Bcnt].Xsize=0.12; Olist[Bcnt].Ysize=0.115; Bcnt++;

  // H2 face
  strcpy(Olist[Bcnt].Name,"surf"); Olist[Bcnt].Xsize=0.12; Olist[Bcnt].Ysize=0.11; Bcnt++;
  strcpy(Olist[Bcnt].Name,"crisp"); Olist[Bcnt].Xsize=0.12; Olist[Bcnt].Ysize=0.11; Bcnt++;
  strcpy(Olist[Bcnt].Name,"iris"); Olist[Bcnt].Xsize=0.12; Olist[Bcnt].Ysize=0.11; Bcnt++;

  // Water face 1
  strcpy(Olist[Bcnt].Name,"polaris"); Olist[Bcnt].Xsize=0.12; Olist[Bcnt].Ysize=0.115; Bcnt++;
  strcpy(Olist[Bcnt].Name,"het"); Olist[Bcnt].Xsize=0.12; Olist[Bcnt].Ysize=0.115; Bcnt++;
  strcpy(Olist[Bcnt].Name,"tosca"); Olist[Bcnt].Xsize=0.12; Olist[Bcnt].Ysize=0.115; Bcnt++;

  // Water face 2
  strcpy(Olist[Bcnt].Name,"maps"); Olist[Bcnt].Xsize=0.12; Olist[Bcnt].Ysize=0.115; Bcnt++;
  strcpy(Olist[Bcnt].Name,"evs"); Olist[Bcnt].Xsize=0.12; Olist[Bcnt].Ysize=0.115; Bcnt++;
  strcpy(Olist[Bcnt].Name,"sxd"); Olist[Bcnt].Xsize=0.12; Olist[Bcnt].Ysize=0.115; Bcnt++;

  // TS1 Generics
  strcpy(Olist[Bcnt].Name,"ch4"); Olist[Bcnt].Xsize=0.12; Olist[Bcnt].Ysize=0.115; Bcnt++;
  strcpy(Olist[Bcnt].Name,"h2"); Olist[Bcnt].Xsize=0.12; Olist[Bcnt].Ysize=0.11; Bcnt++;
  strcpy(Olist[Bcnt].Name,"water"); Olist[Bcnt].Xsize=0.12; Olist[Bcnt].Ysize=0.115; Bcnt++;

  // TS2 Generics
  strcpy(Olist[Bcnt].Name,"groove"); Olist[Bcnt].Xsize=0.08333; Olist[Bcnt].Ysize=0.03; Bcnt++;
  strcpy(Olist[Bcnt].Name,"hydrogen"); Olist[Bcnt].Xsize=0.12; Olist[Bcnt].Ysize=0.11; Bcnt++;
  strcpy(Olist[Bcnt].Name,"narrow"); Olist[Bcnt].Xsize=0.12; Olist[Bcnt].Ysize=0.12; Bcnt++;
  strcpy(Olist[Bcnt].Name,"broad"); Olist[Bcnt].Xsize=0.12; Olist[Bcnt].Ysize=0.12; Bcnt++;

  // TS2 groove
  strcpy(Olist[Bcnt].Name,"e1"); Olist[Bcnt].Xsize=0.08333; Olist[Bcnt].Ysize=0.03; Bcnt++;
  strcpy(Olist[Bcnt].Name,"e2"); Olist[Bcnt].Xsize=0.08333; Olist[Bcnt].Ysize=0.03; Bcnt++;
  strcpy(Olist[Bcnt].Name,"e3"); Olist[Bcnt].Xsize=0.08333; Olist[Bcnt].Ysize=0.03; Bcnt++;
  strcpy(Olist[Bcnt].Name,"e4"); Olist[Bcnt].Xsize=0.08333; Olist[Bcnt].Ysize=0.03; Bcnt++;
  strcpy(Olist[Bcnt].Name,"e5"); Olist[Bcnt].Xsize=0.08333; Olist[Bcnt].Ysize=0.03; Bcnt++;

  //Broad face
  strcpy(Olist[Bcnt].Name,"e6"); Olist[Bcnt].Xsize=0.12; Olist[Bcnt].Ysize=0.12; Bcnt++;
  strcpy(Olist[Bcnt].Name,"e7"); Olist[Bcnt].Xsize=0.12; Olist[Bcnt].Ysize=0.12; Bcnt++;
  strcpy(Olist[Bcnt].Name,"e8"); Olist[Bcnt].Xsize=0.12; Olist[Bcnt].Ysize=0.12; Bcnt++;
  strcpy(Olist[Bcnt].Name,"e9"); Olist[Bcnt].Xsize=0.12; Olist[Bcnt].Ysize=0.12; Bcnt++;
  // Narrow face

  strcpy(Olist[Bcnt].Name,"w1"); Olist[Bcnt].Xsize=0.12; Olist[Bcnt].Ysize=0.12; Bcnt++;
  strcpy(Olist[Bcnt].Name,"w2"); Olist[Bcnt].Xsize=0.12; Olist[Bcnt].Ysize=0.12; Bcnt++;
  strcpy(Olist[Bcnt].Name,"w3"); Olist[Bcnt].Xsize=0.12; Olist[Bcnt].Ysize=0.12; Bcnt++;
  strcpy(Olist[Bcnt].Name,"w4"); Olist[Bcnt].Xsize=0.12; Olist[Bcnt].Ysize=0.12; Bcnt++;

  //Hydrogen face
  strcpy(Olist[Bcnt].Name,"w5"); Olist[Bcnt].Xsize=0.12; Olist[Bcnt].Ysize=0.11; Bcnt++;
  strcpy(Olist[Bcnt].Name,"w6"); Olist[Bcnt].Xsize=0.12; Olist[Bcnt].Ysize=0.11; Bcnt++;
  strcpy(Olist[Bcnt].Name,"w7"); Olist[Bcnt].Xsize=0.12; Olist[Bcnt].Ysize=0.11; Bcnt++;
  strcpy(Olist[Bcnt].Name,"w8"); Olist[Bcnt].Xsize=0.12; Olist[Bcnt].Ysize=0.11; Bcnt++;
  strcpy(Olist[Bcnt].Name,"w9"); Olist[Bcnt].Xsize=0.12; Olist[Bcnt].Ysize=0.11; Bcnt++;


  /* write out version number */
  fprintf(stderr,"**********************************************************************\n");
  fprintf(stderr,"****   This is ISIS_moderator.comp version 2.0 (25/8/05)          ****\n");
  fprintf(stderr,"****   Please check to see if your files are up-to-date           ****\n");
  fprintf(stderr,"****   http://www.isis.rl.ac.uk/Computing/Software/MC/index.htm   ****\n");
  fprintf(stderr,"**********************************************************************\n\n");



  /* convert arguments to runtime variables so that they may be altered */
  rtE0=Emin;
  rtE1=Emax;
  rtmodX=xwidth;
  rtmodY=yheight;


  /* Convert NEGATIVE energy (denoting angstroms) into meV */
  if ( (rtE0<0 && Emax>0) | (rtE0>0 && Emax<0))
    {
      fprintf(stderr,"Cannot have differing signs for Emin and Emax, choose Angstroms or meV!\n");
      exit(1);
    }


  if (rtE0<0 && Emax<0)
    {
      fprintf (stderr,"converting Angstroms to meV\n");
      rtE0=81.793936/(rtE0*rtE0);
      rtE1=81.793936/(rtE1*rtE1);
    }
  if (Lmin && Lmax)
    {
      fprintf (stderr,"converting Angstroms to meV\n");
      rtE0=81.793936/(Lmin*Lmin);
      rtE1=81.793936/(Lmax*Lmax);
    }
  if (rtE0>rtE1)
    {
      tmp=rtE1;
      rtE1=rtE0;
      rtE0=tmp;
      fprintf (stderr,"%g A -> %g A =>  %g meV -> %g meV\n",-Emin,-Emax,rtE0,rtE1);
    }






  /**********************************************************************/

  Tnpts=0;
  Ncount=0;
  fprintf(stderr,"Face == %s \n",Face);

  for(i=0;Face[i] && Face[i]!=' ';i++)
    lowerFace[i]=tolower(Face[i]);
  lowerFace[i]=0;

  for(i=0;i<Bcnt;i++)
    if (strcmp(lowerFace,Olist[i].Name)==0)
      {
	if (rtmodX<=0.0)
	  {
	    rtmodX=Olist[i].Xsize;
	    fprintf(stderr,"default xwidth used %g m\n",rtmodX);
	  }
	if (rtmodY<=0.0)
	  {
	    rtmodY=Olist[i].Ysize;
	    fprintf(stderr,"default yheight used %g m\n",rtmodY);
	  }
	/* Input file naming according to "beamline" list above */
	if (i < 18) {
	  sprintf(fname,"TS1.%s",Olist[i].Name);
	  TargetStation = 1;
	} else {
	  sprintf(fname,"TS2.%s",Olist[i].Name);
	  TargetStation = 2;
	}
	scaleSize=(SAC) ? 1.0 : rtmodY*rtmodX*10000.0;
	break;
      }

  if(i==Bcnt)   /* Error condition */
    {
      fprintf(stderr,"Unknown moderator type ::%s::\n",lowerFace);
      fprintf(stderr,"Valid options == > \n");
      for(i=0;i<Bcnt;i++)
	{
	  fprintf(stderr," %s ",Olist[i].Name);
/* 	  if (!((i+1) % 4)) */
	    fprintf(stderr,"\n");
	}
      scaleSize=xwidth*yheight/0.0025;
      exit(1);
    }

  rtmodY*=cos(CAngle);

  /* READ PARAMETER FILE */

  TFile=openFile(fname);
  
  if (!readHtable(TFile,rtE0,rtE1,&TS))
    {
      fprintf(stderr,"Failed to read the Hzone from file %s\n", fname);
      exit(1);
    }
  fclose(TFile);

  fprintf(stderr,"nEnergy == %d\n",TS.nEnergy);

  /* Do solid angle correction if required */
  // if SAC=0/1 solid angle is determined
  if (SAC)
    angleArea=(dist>0.0) ? strArea(dist, rtmodX, rtmodY, focus_xw, focus_yh) : 2*3.141592654;
  else
    angleArea=1.0;
  
  /* 
  TS1: MCNPX runs were done for 60 mu-A, but the source runs at 160 mu-A, 40 Hz.
  TS2: MCNPX runs were done for 60 mu-A, but the source runs at 40-mu-A, 10 Hz.
  */
  
  if (TargetStation == 1) {
    CurrentWeight = 160.0/60.0;
  } else {
    CurrentWeight = 40.0/60.0;
  }
  
  #undef Face
  #undef Emin
  #undef Emax
  #undef dist
  #undef focus_xw
  #undef focus_yh
  #undef xwidth
  #undef yheight
  #undef CAngle
  #undef SAC
  #undef Lmin
  #undef Lmax
  #undef target_index
  #undef verbose
  #undef p_in
  #undef Tnpts
  #undef scaleSize
  #undef angleArea
  #undef Nsim
  #undef Ncount
  #undef TS
  #undef rtE0
  #undef rtE1
  #undef rtmodX
  #undef rtmodY
  #undef TargetStation
  #undef CurrentWeight
  return(_comp);
} /* class_ISIS_moderator_init */

_class_Slit *class_Slit_init(_class_Slit *_comp
) {
  #define xmin (_comp->_parameters.xmin)
  #define xmax (_comp->_parameters.xmax)
  #define ymin (_comp->_parameters.ymin)
  #define ymax (_comp->_parameters.ymax)
  #define radius (_comp->_parameters.radius)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
if (xwidth > 0)  { xmax=xwidth/2;  xmin=-xmax; }
  if (yheight > 0) { ymax=yheight/2; ymin=-ymax; }
  if (xmin == 0 && xmax == 0 && ymin == 0 && ymax == 0 && radius == 0)
    { fprintf(stderr,"Slit: %s: Error: give geometry\n", NAME_CURRENT_COMP); exit(-1); }

  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef radius
  #undef xwidth
  #undef yheight
  return(_comp);
} /* class_Slit_init */

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

_class_Guide_tapering *class_Guide_tapering_init(_class_Guide_tapering *_comp
) {
  #define option (_comp->_parameters.option)
  #define w1 (_comp->_parameters.w1)
  #define h1 (_comp->_parameters.h1)
  #define l (_comp->_parameters.l)
  #define linw (_comp->_parameters.linw)
  #define loutw (_comp->_parameters.loutw)
  #define linh (_comp->_parameters.linh)
  #define louth (_comp->_parameters.louth)
  #define R0 (_comp->_parameters.R0)
  #define Qcx (_comp->_parameters.Qcx)
  #define Qcy (_comp->_parameters.Qcy)
  #define alphax (_comp->_parameters.alphax)
  #define alphay (_comp->_parameters.alphay)
  #define W (_comp->_parameters.W)
  #define mx (_comp->_parameters.mx)
  #define my (_comp->_parameters.my)
  #define segno (_comp->_parameters.segno)
  #define curvature (_comp->_parameters.curvature)
  #define curvature_v (_comp->_parameters.curvature_v)
  #define w1c (_comp->_parameters.w1c)
  #define w2c (_comp->_parameters.w2c)
  #define ww (_comp->_parameters.ww)
  #define hh (_comp->_parameters.hh)
  #define whalf (_comp->_parameters.whalf)
  #define hhalf (_comp->_parameters.hhalf)
  #define lwhalf (_comp->_parameters.lwhalf)
  #define lhhalf (_comp->_parameters.lhhalf)
  #define h1_in (_comp->_parameters.h1_in)
  #define h2_out (_comp->_parameters.h2_out)
  #define w1_in (_comp->_parameters.w1_in)
  #define w2_out (_comp->_parameters.w2_out)
  #define l_seg (_comp->_parameters.l_seg)
  #define seg (_comp->_parameters.seg)
  #define h12 (_comp->_parameters.h12)
  #define h2 (_comp->_parameters.h2)
  #define w12 (_comp->_parameters.w12)
  #define w2 (_comp->_parameters.w2)
  #define a_ell_q (_comp->_parameters.a_ell_q)
  #define b_ell_q (_comp->_parameters.b_ell_q)
  #define lbw (_comp->_parameters.lbw)
  #define lbh (_comp->_parameters.lbh)
  #define mxi (_comp->_parameters.mxi)
  #define u1 (_comp->_parameters.u1)
  #define u2 (_comp->_parameters.u2)
  #define div1 (_comp->_parameters.div1)
  #define p2_para (_comp->_parameters.p2_para)
  #define test (_comp->_parameters.test)
  #define Div1 (_comp->_parameters.Div1)
  #define seg (_comp->_parameters.seg)
  #define fu (_comp->_parameters.fu)
  #define pos (_comp->_parameters.pos)
  #define file_name (_comp->_parameters.file_name)
  #define ep (_comp->_parameters.ep)
  #define num (_comp->_parameters.num)
  #define rotation_h (_comp->_parameters.rotation_h)
  #define rotation_v (_comp->_parameters.rotation_v)
int i,ii;

rotation_h=0;
rotation_v=0;

// dynamic memory allocation is good
w1c = (double*)malloc(sizeof(double)*segno);
  w2c = (double*)malloc(sizeof(double)*segno);
  ww = (double*)malloc(sizeof(double)*segno);
  hh = (double*)malloc(sizeof(double)*segno);
  whalf = (double*)malloc(sizeof(double)*segno);
  hhalf = (double*)malloc(sizeof(double)*segno);
  lwhalf = (double*)malloc(sizeof(double)*segno);
  lhhalf = (double*)malloc(sizeof(double)*segno);
  h1_in = (double*)malloc(sizeof(double)*(segno+1));
  h2_out = (double*)malloc(sizeof(double)*(segno+1));
  w1_in = (double*)malloc(sizeof(double)*(segno+1));
  w2_out = (double*)malloc(sizeof(double)*(segno+1));

  struct para {
    char st[128];
  } segment[800];
  
  if (W <=0)
  {
    fprintf(stderr,"Component: %s (Guide_tapering) W must \n", NAME_CURRENT_COMP);
    fprintf(stderr,"           be positive\n");
    exit(-1);
  }
  if (l <= 0)
  {
    fprintf(stderr,"Component: %s (Guide_tapering) real guide length \n",
    NAME_CURRENT_COMP);
    fprintf(stderr,"           is <= ZERO ! \n");
    exit(-1);
  }
  if (mcgravitation) fprintf(stderr,"WARNING: Guide_tapering: %s: "
    "This component produces wrong results with gravitation !\n"
    "Use Guide_gravity.\n",
    NAME_CURRENT_COMP);
  seg=segno;
  l_seg=l/(seg);
  h12 = h1/2.0;
  if (option != NULL)
  {
     fu = (char*)malloc(sizeof(char)*(strlen(option)+1));
     strcpy(fu,option);
  } else {
     exit(-1);
  }
  /* handle guide geometry ================================================== */
  if (!strcmp(fu,"elliptical"))
  {
     /* calculate parameter b of elliptical equestion - vertical mirrors */
     /* (l+linh+louth) -> distance between focal points */
     /*  printf("A1 \n"); */
     lbh = l + linh + louth;
     if (linh == 0 && louth == 0 )
     {
        /* plane mirrors (vertical) */
        b_ell_q = 0;
        h2 = h1;
     } else {
        /* elliptical mirrors */
        u1 = sqrt((linh*linh)+(h12*h12));
        u2 = sqrt((h12*h12) + ((l+louth)*(l+louth)));
        a_ell_q = ((u1 + u2)/2.0)*((u1 + u2)/2.0);
        b_ell_q = a_ell_q - ((lbh/2.0)*(lbh/2.0));
        /* calculate heigth of guide exit (h2) */
        div1 = ((lbh/2.0-louth)*(lbh/2.0-louth))/a_ell_q;
        h2 = sqrt(b_ell_q*(1.0-div1));
        h2 = h2*2.0;
     }
  } else if (!strcmp(fu,"parabolical")) {
     if ((linh > 0) && (louth > 0))
     {
       fprintf(stderr,"Component: %s (Guide_tapering) Two focal\n",NAME_CURRENT_COMP);
       fprintf(stderr,"            points lout and linh are not allowed! \n");
        free(fu);exit(-1);
     }
     if (louth == 0 && linh == 0)
     {
        /* plane mirrors (vertical) */
        h2 = h1;
     } else {
        /* parabolical mirrors */
        if (linh == 0)
        {
           Div1=((2.0*louth+2.0*l)*(2.0*louth+2.0*l))/4.0;
           p2_para=((sqrt(Div1+(h12*h12)))-(louth+l))*2.0;
           /* calculate heigth of guide exit (h2) */
           h2 = sqrt(p2_para*(louth+p2_para/4.0));
           h2 = h2*2.0;
         } else {
            /* anti-trompete */
           Div1=((2.0*linh)*(2.0*linh))/4.0;
           p2_para=((sqrt(Div1+(h12*h12)))-linh)*2.0;
           /* calculate heigth of guide exit (h2) */
           h2 = sqrt(p2_para*(l+linh+p2_para/4.0));
           h2 = h2*2.0;
         }
     }
  } else if (!strncmp(fu,"file",4)) {
     pos = strtok(fu,"=");
     while (pos=strtok(0,"="))
     {
        strcpy(file_name,pos);
     }
     if ((num=fopen(file_name,"r")) == NULL)
     {
        fprintf(stderr,"Component: %s (Guide_tapering)\n",NAME_CURRENT_COMP);
        fprintf(stderr,"           File %s not found! \n", file_name);
         free(fu);exit(-1);
     } else {
        ii = 0;
        while (!feof(num))
        {
          char *ret = fgets(segment[ii].st,128,num);
          if (ii >  799 || !ret) {
             fprintf(stderr,"%s: Number of segments is limited to 800 !! \n",NAME_CURRENT_COMP);
              free(fu);exit(-1);
          }
          ii++;
        }
        fclose(num);
        ii--;
     }
     seg = ii-3;
     l_seg=l/seg;
     for (i=3;i<ii;i++)
     {
        if (strlen(segment[i].st) < 4)
        {
          fprintf(stderr,"Component: %s (Guide_tapering)\n",NAME_CURRENT_COMP);
          fprintf(stderr,"           Data Format Error! \n");
          free(fu);exit(-1);
        }
        h1_in[i-3] = strtod(strtok(segment[i].st," "), &ep);
        h2_out[i-3] = strtod(strtok(0," "), &ep);
        w1_in[i-3] = strtod(strtok(0," "), &ep);
        w2_out[i-3] = strtod(strtok(0," "), &ep);
     }
     h1 = h1_in[0];
     h2 = h2_out[seg-1];
     w1 = w1_in[0];
     w2 = w2_out[seg-1];
     for (i=0;i<seg;i++)
     {
      fprintf(stderr,"%d: %lf %lf %lf %lf \n",i,h1_in[i],h2_out[i],w1_in[i],w2_out[i]);
     }
  } else if (!strcmp(fu,"straight")) {
    for (i=0;i<seg;i++) {
      h1_in[i] = h2_out[i] = h2 = h1;
      w1_in[i] = w2_out[i] = w2 = w1;
    }
  } else {
     fprintf(stderr,"Component: %s (Guide_tapering)\n",NAME_CURRENT_COMP);
     fprintf(stderr,"           Unknown KEYWORD: %s \n", fu);
     free(fu);exit(-1);
  }
  fprintf(stderr,"Component: %s (Guide_tapering)\n",NAME_CURRENT_COMP);
  fprintf(stderr,"           Height at the guide exit (h2): %lf \n", h2);
  if (h2 <= 0)
  {
   fprintf(stderr,"Component: %s (Guide_tapering)\n", NAME_CURRENT_COMP);
   fprintf(stderr,"           Height at the guide exit (h2) was calculated\n");
   fprintf(stderr,"           <=0; Please change the parameter h1 and/or\n");
   fprintf(stderr,"           linh and/or louth! \n");
    free(fu);exit(-1);
  }
  if (!strcmp(fu,"elliptical"))
  {
     h1_in[0] = h1;
     for (i=1;i<seg;i++)
     {
       if (b_ell_q == 0)
       {
         h1_in[i]=h1;
       } else {
         mxi = (((lbh/2.0)-linh) - (l_seg * i));
         h1_in[i] = (sqrt((1.0-((mxi*mxi)/a_ell_q))*b_ell_q))*2.0;
       }
     h2_out[i-1] = h1_in[i];
     }
     h2_out[seg-1]=h2;
  } else if (!strcmp(fu,"parabolical")) {
     h1_in[0] = h1;
     ii=seg-1;
     if (louth == 0 && linh == 0)
     {
        for (i=1;i<(seg+1);i++)
        {
           h1_in[i]=h1;
           ii=ii-1;
           h2_out[i-1] = h1_in[i];
        }
     } else {
        if ((linh == 0) && (louth > 0))
        {
           for (i=1;i<(seg+1);i++)
           {
             h1_in[i] = (sqrt((p2_para/4.0+louth+(l_seg*ii))*p2_para))*2.0;
             ii=ii-1;
             h2_out[i-1] = h1_in[i];
           }
        } else {
           for (i=1;i<(seg+1);i++)
           {
             h1_in[i] = (sqrt((p2_para/4.0+linh+(l_seg*i))*p2_para))*2.0;
             h2_out[i-1] = h1_in[i];
           }
        }
     }
  }
  /* compute each value for horizontal mirrors */
  w12 = w1/2.0;
  if (!strcmp(fu,"elliptical"))
  {
    /* calculate lbw the distance between focal points of horizontal mirrors */
    lbw = l + linw + loutw;
    /* calculate parameter b of elliptical equestion - horizontal mirrors */
    if (linw == 0 && loutw == 0 )
    {
       /* plane mirrors (horizontal) */
       b_ell_q = 0;
       w2 = w1;
    } else {
       /* elliptical mirrors */
       u1 = sqrt((linw*linw)+(w12*w12));
       u2 = sqrt((w12*w12) + ((l+loutw)*(l+loutw)));
       a_ell_q = ((u1 + u2)/2.0)*((u1 + u2)/2.0);
       b_ell_q = a_ell_q - ((lbw/2.0)*(lbw/2.0));
       /* calculate weigth of guide exit (w2) */
       div1 = ((lbw/2.0-loutw)*(lbw/2.0-loutw))/a_ell_q;
       w2 = sqrt(b_ell_q*(1.0-div1));
       w2 = w2*2.0;
     }
  } else if (!strcmp(fu,"parabolical")) {
     if ((linw > 0) && (loutw > 0))
     {
       fprintf(stderr,"Component: %s (Guide_tapering) Two focal\n",NAME_CURRENT_COMP);
       fprintf(stderr,"           points linw and loutw are not allowed! \n");
         free(fu);exit(-1);
     }
     if (loutw == 0 && linw == 0)
     {
        /* plane mirrors (horizontal) */
        w2 = w1;
     } else {
       if (linw == 0)
       {
          /* parabolical mirrors */
          Div1=((2.0*loutw+2.0*l)*(2.0*loutw+2.0*l))/4.0;
          p2_para=((sqrt(Div1+(w12*w12)))-(loutw+l))*2.0;
          /* calculate weigth of guide exit (w2) */
          w2 = sqrt(p2_para*(loutw+p2_para/4.0));
          w2 = w2*2.0;
       } else {
          /* anti-trompete */
          Div1=((2.0*linw)*(2.0*linw))/4.0;
          p2_para=((sqrt(Div1+(w12*w12)))-linw)*2.0;
          /* calculate heigth of guide exit (w2) */
          w2 = sqrt(p2_para*(l+linw+p2_para/4.0));
          w2 = w2*2.0;
       }
     }
  }
  fprintf(stderr,"Component: %s (Guide_tapering)\n",NAME_CURRENT_COMP);
  fprintf(stderr,"           Width at the guide exit (w2): %lf \n", w2);
  if (w2 <= 0)
  {
   fprintf(stderr,"Component: %s (Guide_tapering)\n", NAME_CURRENT_COMP);
   fprintf(stderr,"           Width at the guide exit (w2) was calculated\n");
   fprintf(stderr,"           <=0; Please change the parameter w1 and/or\n");
   fprintf(stderr,"           l! \n");
    free(fu);exit(-1);
  }
  if (!strcmp(fu,"elliptical"))
  {
     w1_in[0]=w1;
     for (i=1;i<seg;i++)
     {
       if (b_ell_q == 0)
       {
         w1_in[i]=w1;
       } else {
         mxi = (((lbw/2.0)-linw) - (l_seg * i));
         w1_in[i] = (sqrt((1.0-((mxi*mxi)/a_ell_q))*b_ell_q))*2.0;
       }
       w2_out[i-1] = w1_in[i];
     }
     w2_out[seg-1]=w2;
  } else if (!strcmp(fu,"parabolical")) {
     w1_in[0]=w1;
     ii=seg-1;
     if (loutw == 0 && linw == 0)
     {
        for (i=1;i<(seg+1);i++)
        {
           w1_in[i]=w1;
           ii=ii-1;
           w2_out[i-1] = w1_in[i];
        }
     } else {
        if ((linw == 0) && (loutw > 0))
        {
           for (i=1;i<(seg+1);i++)
           {
             w1_in[i] = (sqrt((p2_para/4+loutw+(l_seg*ii))*p2_para))*2;
             ii=ii-1;
             w2_out[i-1] = w1_in[i];
           }
        } else {
           for (i=1;i<(seg+1);i++)
           {
             w1_in[i] = (sqrt((p2_para/4+linw+(l_seg*i))*p2_para))*2;
             w2_out[i-1] = w1_in[i];
           }
        }
     }
  }
  free(fu);
  for (i=0;i<seg;i++)
  {
    w1c[i] = w1_in[i];
    w2c[i] = w2_out[i];
    ww[i] = .5*(w2c[i] - w1c[i]);
    hh[i] = .5*(h2_out[i] - h1_in[i]);
    whalf[i] = .5*w1c[i];
    hhalf[i] = .5*h1_in[i];
    lwhalf[i] = l_seg*whalf[i];
    lhhalf[i] = l_seg*hhalf[i];
  }
  /* guide curvature: rotation angle [rad] between each guide segment */
  if (curvature && l && segno)   rotation_h = l/curvature/segno;
  if (curvature_v && l && segno) rotation_v = l/curvature_v/segno;
  #undef option
  #undef w1
  #undef h1
  #undef l
  #undef linw
  #undef loutw
  #undef linh
  #undef louth
  #undef R0
  #undef Qcx
  #undef Qcy
  #undef alphax
  #undef alphay
  #undef W
  #undef mx
  #undef my
  #undef segno
  #undef curvature
  #undef curvature_v
  #undef w1c
  #undef w2c
  #undef ww
  #undef hh
  #undef whalf
  #undef hhalf
  #undef lwhalf
  #undef lhhalf
  #undef h1_in
  #undef h2_out
  #undef w1_in
  #undef w2_out
  #undef l_seg
  #undef seg
  #undef h12
  #undef h2
  #undef w12
  #undef w2
  #undef a_ell_q
  #undef b_ell_q
  #undef lbw
  #undef lbh
  #undef mxi
  #undef u1
  #undef u2
  #undef div1
  #undef p2_para
  #undef test
  #undef Div1
  #undef seg
  #undef fu
  #undef pos
  #undef file_name
  #undef ep
  #undef num
  #undef rotation_h
  #undef rotation_v
  return(_comp);
} /* class_Guide_tapering_init */

_class_Multilayer_Sample *class_Multilayer_Sample_init(_class_Multilayer_Sample *_comp
) {
  #define xwidth (_comp->_parameters.xwidth)
  #define zlength (_comp->_parameters.zlength)
  #define nlayer (_comp->_parameters.nlayer)
  #define sldPar (_comp->_parameters.sldPar)
  #define dPar (_comp->_parameters.dPar)
  #define sigmaPar (_comp->_parameters.sigmaPar)
  #define frac_inc (_comp->_parameters.frac_inc)
  #define ythick (_comp->_parameters.ythick)
  #define mu_inc (_comp->_parameters.mu_inc)
  #define target_index (_comp->_parameters.target_index)
  #define focus_xw (_comp->_parameters.focus_xw)
  #define focus_yh (_comp->_parameters.focus_yh)
  #define tx (_comp->_parameters.tx)
  #define ty (_comp->_parameters.ty)
  #define tz (_comp->_parameters.tz)
  #define xmin (_comp->_parameters.xmin)
  #define xmax (_comp->_parameters.xmax)
  #define zmin (_comp->_parameters.zmin)
  #define zmax (_comp->_parameters.zmax)
if (frac_inc>0) {
    if (!(ythick) || !(mu_inc)) {
      fprintf(stderr,"Multilayer: error: %s: You requested a non-meaningful combination of frac_inc, ythick, mu_inc. EXIT\n", NAME_CURRENT_COMP);
      exit(1);
    }
  }
  xmin = -xwidth/2.0;
  xmax =  xwidth/2.0;
  zmin = -zlength/2.0;
  zmax =  zlength/2.0;
  if (target_index) {
    Coords ToTarget;
    ToTarget = coords_sub(POS_A_COMP_INDEX(INDEX_CURRENT_COMP+target_index),POS_A_CURRENT_COMP);
    ToTarget = rot_apply(ROT_A_CURRENT_COMP, ToTarget);
    coords_get(ToTarget, &tx, &ty, &tz);
  } else {
    tx = 0; ty = 0; tz = 0;
  }

  #undef xwidth
  #undef zlength
  #undef nlayer
  #undef sldPar
  #undef dPar
  #undef sigmaPar
  #undef frac_inc
  #undef ythick
  #undef mu_inc
  #undef target_index
  #undef focus_xw
  #undef focus_yh
  #undef tx
  #undef ty
  #undef tz
  #undef xmin
  #undef xmax
  #undef zmin
  #undef zmax
  return(_comp);
} /* class_Multilayer_Sample_init */



int init(void) { /* called by mccode_main for ISIS_CRISP:INITIALISE */
  DEBUG_INSTR();

  /* code_main/parseoptions/readparams sets instrument parameters value */
  stracpy(instrument->_name, "ISIS_CRISP", 256);

  /* Instrument 'ISIS_CRISP' INITIALISE */
  SIG_MESSAGE("[ISIS_CRISP] INITIALISE [ISIS_CRISP.instr:50]");
  #define glen (instrument->_parameters._glen)
  #define flen (instrument->_parameters._flen)
  #define w1 (instrument->_parameters._w1)
  #define vm (instrument->_parameters._vm)
  #define FRAC (instrument->_parameters._FRAC)
{
  spos=10.2055;
  l=glen;
  loutw=flen;
  linw=10.0;
  w12=w1/2.0;

// calculate the width of the guide exit
  lbw = l + linw + loutw;
  u1 = sqrt((linw*linw)+(w12*w12));
  u2 = sqrt((w12*w12) + ((l+loutw)*(l+loutw)));
  a_ell_q = ((u1 + u2)/2)*((u1 + u2)/2);
  b_ell_q = a_ell_q - ((lbw/2)*(lbw/2));

/* calculate width of guide exit (w2) */
  div1 = ((lbw/2-loutw)*(lbw/2-loutw))/a_ell_q;
  w2 = sqrt(b_ell_q*(1-div1));
  w2 = w2*2;

  tend=w2;
  t15=tan(-1.5*DEG2RAD);
}
  #undef glen
  #undef flen
  #undef w1
  #undef vm
  #undef FRAC
  _a1_setpos(); /* type Arm */
  _isis_source_setpos(); /* type ISIS_moderator */
  _defslit1_setpos(); /* type Slit */
  _defslit2_setpos(); /* type Slit */
  _coarseslit1_setpos(); /* type Slit */
  _coarseslit2_setpos(); /* type Slit */
  _lmon1_setpos(); /* type L_monitor */
  _PSDmon1_setpos(); /* type PSD_monitor */
  _slit1_setpos(); /* type Slit */
  _PSDmons1_setpos(); /* type PSD_monitor */
  _eguide1_setpos(); /* type Guide_tapering */
  _slit2_setpos(); /* type Slit */
  _lmon2_setpos(); /* type L_monitor */
  _samp1_setpos(); /* type Multilayer_Sample */
  _slit3_setpos(); /* type Slit */
  _eguide2_setpos(); /* type Guide_tapering */
  _lmon3_setpos(); /* type L_monitor */
  _PSDdet1_setpos(); /* type PSD_monitor */

  /* call iteratively all components INITIALISE */

  class_ISIS_moderator_init(_isis_source);

  class_Slit_init(_defslit1);

  class_Slit_init(_defslit2);

  class_Slit_init(_coarseslit1);

  class_Slit_init(_coarseslit2);

  class_L_monitor_init(_lmon1);

  class_PSD_monitor_init(_PSDmon1);

  class_Slit_init(_slit1);

  class_PSD_monitor_init(_PSDmons1);

  class_Guide_tapering_init(_eguide1);

  class_Slit_init(_slit2);

  class_L_monitor_init(_lmon2);

  class_Multilayer_Sample_init(_samp1);

  class_Slit_init(_slit3);

  class_Guide_tapering_init(_eguide2);

  class_L_monitor_init(_lmon3);

  class_PSD_monitor_init(_PSDdet1);

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
_class_ISIS_moderator *class_ISIS_moderator_trace(_class_ISIS_moderator *_comp
  , _class_particle *_particle) {
  ABSORBED=SCATTERED=RESTORE=0;

  #define Face (_comp->_parameters.Face)
  #define Emin (_comp->_parameters.Emin)
  #define Emax (_comp->_parameters.Emax)
  #define dist (_comp->_parameters.dist)
  #define focus_xw (_comp->_parameters.focus_xw)
  #define focus_yh (_comp->_parameters.focus_yh)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define CAngle (_comp->_parameters.CAngle)
  #define SAC (_comp->_parameters.SAC)
  #define Lmin (_comp->_parameters.Lmin)
  #define Lmax (_comp->_parameters.Lmax)
  #define target_index (_comp->_parameters.target_index)
  #define verbose (_comp->_parameters.verbose)
  #define p_in (_comp->_parameters.p_in)
  #define Tnpts (_comp->_parameters.Tnpts)
  #define scaleSize (_comp->_parameters.scaleSize)
  #define angleArea (_comp->_parameters.angleArea)
  #define Nsim (_comp->_parameters.Nsim)
  #define Ncount (_comp->_parameters.Ncount)
  #define TS (_comp->_parameters.TS)
  #define rtE0 (_comp->_parameters.rtE0)
  #define rtE1 (_comp->_parameters.rtE1)
  #define rtmodX (_comp->_parameters.rtmodX)
  #define rtmodY (_comp->_parameters.rtmodY)
  #define TargetStation (_comp->_parameters.TargetStation)
  #define CurrentWeight (_comp->_parameters.CurrentWeight)
  double v,r,E;
  double xf,yf,dx,dy,w_focus;    /* mxp ->max var in param space */
  double Ival,Tval,Eval;
  double Ddist;   /* Temp versions of dist */


  Ncount++;

  p=p_in;

  p=1.0;         /* forcing */
  z=0;
  x = 0.5*rtmodX*randpm1();            /* Get point +/-0.5 *  */
  y = 0.5*rtmodY*randpm1();
  xf = 0.5*focus_xw*randpm1();          /* Choose focusing position uniformly */
  yf = 0.5*focus_yh*randpm1();
  dx = xf-x;
  dy = yf-y;
  if (dist>0.0)
    {
      r = sqrt(dx*dx+dy*dy+dist*dist);                 /* Actual distance to point */
      Ddist=dist;
      w_focus = (SAC) ? angleArea : scaleSize*(dist*dist)/(r*r);
    }
  else   /* Assume that we have a window 1metre infront of the moderator */
	 /*   with size area of detector and solid angle 1.0 */
    {
      r=1.0;
      w_focus=scaleSize;
      Ddist=1.0;
    }

  getPoint(&Tval,&Eval,&rtE0,&rtE1, TS);

  //fprintf(stderr,"outside %g mev\n", TS.Total );
  if(Eval>rtE1 || Eval<rtE0)
      fprintf(stderr,"outside %g mev\n", Eval );

  Ival=TS.Total*3.744905847e14*1.1879451;  /* ( of proton in 60uAmp) * (1-cos(30))*2*Pi  */


  v = SE2V*sqrt(Eval);      /* Calculate the velocity */
  vz = v*Ddist/r;
  vy = v*dy/r;
  vx = v*dx/r;

  if (Ncount==1)
    fprintf(stderr,"Totals:: %g %d %d \n",TS.Total,TS.nEnergy,TS.nTime);
  if (!(Ncount % 100000) && verbose)
    fprintf(stderr,"FF[%d]=> %g %g %g %g \n",Ncount,Eval,Tval,TS.Total,Ival);

  t=Tval;
  
  p=w_focus*Ival*CurrentWeight/Nsim;
  #undef Face
  #undef Emin
  #undef Emax
  #undef dist
  #undef focus_xw
  #undef focus_yh
  #undef xwidth
  #undef yheight
  #undef CAngle
  #undef SAC
  #undef Lmin
  #undef Lmax
  #undef target_index
  #undef verbose
  #undef p_in
  #undef Tnpts
  #undef scaleSize
  #undef angleArea
  #undef Nsim
  #undef Ncount
  #undef TS
  #undef rtE0
  #undef rtE1
  #undef rtmodX
  #undef rtmodY
  #undef TargetStation
  #undef CurrentWeight
  return(_comp);
} /* class_ISIS_moderator_trace */

#pragma acc routine seq
_class_Slit *class_Slit_trace(_class_Slit *_comp
  , _class_particle *_particle) {
  ABSORBED=SCATTERED=RESTORE=0;

  #define xmin (_comp->_parameters.xmin)
  #define xmax (_comp->_parameters.xmax)
  #define ymin (_comp->_parameters.ymin)
  #define ymax (_comp->_parameters.ymax)
  #define radius (_comp->_parameters.radius)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
    mcPROP_Z0;
    if (((radius == 0) && (x<xmin || x>xmax || y<ymin || y>ymax))
    || ((radius != 0) && (x*x + y*y > radius*radius)))
      ABSORB;
    else
        SCATTER;
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef radius
  #undef xwidth
  #undef yheight
  return(_comp);
} /* class_Slit_trace */

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

#pragma acc routine seq
_class_Guide_tapering *class_Guide_tapering_trace(_class_Guide_tapering *_comp
  , _class_particle *_particle) {
  ABSORBED=SCATTERED=RESTORE=0;

  #define option (_comp->_parameters.option)
  #define w1 (_comp->_parameters.w1)
  #define h1 (_comp->_parameters.h1)
  #define l (_comp->_parameters.l)
  #define linw (_comp->_parameters.linw)
  #define loutw (_comp->_parameters.loutw)
  #define linh (_comp->_parameters.linh)
  #define louth (_comp->_parameters.louth)
  #define R0 (_comp->_parameters.R0)
  #define Qcx (_comp->_parameters.Qcx)
  #define Qcy (_comp->_parameters.Qcy)
  #define alphax (_comp->_parameters.alphax)
  #define alphay (_comp->_parameters.alphay)
  #define W (_comp->_parameters.W)
  #define mx (_comp->_parameters.mx)
  #define my (_comp->_parameters.my)
  #define segno (_comp->_parameters.segno)
  #define curvature (_comp->_parameters.curvature)
  #define curvature_v (_comp->_parameters.curvature_v)
  #define w1c (_comp->_parameters.w1c)
  #define w2c (_comp->_parameters.w2c)
  #define ww (_comp->_parameters.ww)
  #define hh (_comp->_parameters.hh)
  #define whalf (_comp->_parameters.whalf)
  #define hhalf (_comp->_parameters.hhalf)
  #define lwhalf (_comp->_parameters.lwhalf)
  #define lhhalf (_comp->_parameters.lhhalf)
  #define h1_in (_comp->_parameters.h1_in)
  #define h2_out (_comp->_parameters.h2_out)
  #define w1_in (_comp->_parameters.w1_in)
  #define w2_out (_comp->_parameters.w2_out)
  #define l_seg (_comp->_parameters.l_seg)
  #define seg (_comp->_parameters.seg)
  #define h12 (_comp->_parameters.h12)
  #define h2 (_comp->_parameters.h2)
  #define w12 (_comp->_parameters.w12)
  #define w2 (_comp->_parameters.w2)
  #define a_ell_q (_comp->_parameters.a_ell_q)
  #define b_ell_q (_comp->_parameters.b_ell_q)
  #define lbw (_comp->_parameters.lbw)
  #define lbh (_comp->_parameters.lbh)
  #define mxi (_comp->_parameters.mxi)
  #define u1 (_comp->_parameters.u1)
  #define u2 (_comp->_parameters.u2)
  #define div1 (_comp->_parameters.div1)
  #define p2_para (_comp->_parameters.p2_para)
  #define test (_comp->_parameters.test)
  #define Div1 (_comp->_parameters.Div1)
  #define seg (_comp->_parameters.seg)
  #define fu (_comp->_parameters.fu)
  #define pos (_comp->_parameters.pos)
  #define file_name (_comp->_parameters.file_name)
  #define ep (_comp->_parameters.ep)
  #define num (_comp->_parameters.num)
  #define rotation_h (_comp->_parameters.rotation_h)
  #define rotation_v (_comp->_parameters.rotation_v)
  double t1,t2,ts,zr;                           /* Intersection times. */
  double av,ah,bv,bh,cv1,cv2,ch1,ch2,dd;        /* Intermediate values */
  double vdotn_v1,vdotn_v2,vdotn_h1,vdotn_h2;   /* Dot products. */
  int i;                                        /* Which mirror hit? */
  double q;                                     /* Q [1/AA] of reflection */
  double vlen2,nlen2;                           /* Vector lengths squared */
  double edge;
  double hadj;                                  /* Channel displacement */
  int ii;

  /* Propagate neutron to guide entrance. */
  PROP_Z0;
  for (ii=0;ii<seg;ii++)
  {
    zr=ii*l_seg;
    /* Propagate neutron to segment entrance. */
    ts=(zr-z)/vz;
    PROP_DT(ts);
    if(x <= w1_in[ii]/-2.0 || x >= w1_in[ii]/2.0 || y <= -hhalf[ii] || y >= hhalf[ii])
      ABSORB;
    /* Shift origin to center of channel hit (absorb if hit dividing walls) */
    x += w1_in[ii]/2.0;
    edge = floor(x/w1c[ii])*w1c[ii];
    if(x - edge > w1c[ii])
    {
      x -= w1_in[ii]/2.0; /* Re-adjust origin */
      ABSORB;
    }
    x -= (edge + (w1c[ii]/2.0));
    hadj = edge + (w1c[ii]/2.0) - w1_in[ii]/2.0;
    for(;;)
    {
      /* Compute the dot products of v and n for the four mirrors. */
      ts=(zr-z)/vz;
      av = l_seg*vx; bv = ww[ii]*vz;
      ah = l_seg*vy; bh = hh[ii]*vz;
      vdotn_v1 = bv + av;         /* Left vertical */
      vdotn_v2 = bv - av;         /* Right vertical */
      vdotn_h1 = bh + ah;         /* Lower horizontal */
      vdotn_h2 = bh - ah;         /* Upper horizontal */
      /* Compute the dot products of (O - r) and n as c1+c2 and c1-c2 */
      cv1 = -whalf[ii]*l_seg - (z-zr)*ww[ii]; cv2 = x*l_seg;
      ch1 = -hhalf[ii]*l_seg - (z-zr)*hh[ii]; ch2 = y*l_seg;
      /* Compute intersection times. */
      t1 = (zr + l_seg - z)/vz;
      i = 0;
      if(vdotn_v1 < 0 && (t2 = (cv1 - cv2)/vdotn_v1) < t1)
      {
        t1 = t2;
        i = 1;
      }
      if(vdotn_v2 < 0 && (t2 = (cv1 + cv2)/vdotn_v2) < t1)
      {
        t1 = t2;
        i = 2;
      }
      if(vdotn_h1 < 0 && (t2 = (ch1 - ch2)/vdotn_h1) < t1)
      {
        t1 = t2;
        i = 3;
      }
      if(vdotn_h2 < 0 && (t2 = (ch1 + ch2)/vdotn_h2) < t1)
      {
        t1 = t2;
        i = 4;
      }
      if(i == 0)
      {
        break;                    /* Neutron left guide. */
      }
      PROP_DT(t1);
      switch(i)
      {
        case 1:                   /* Left vertical mirror */
          nlen2 = l_seg*l_seg + ww[ii]*ww[ii];
          q = V2Q*(-2)*vdotn_v1/sqrt(nlen2);
          dd = 2*vdotn_v1/nlen2;
          vx = vx - dd*l_seg;
          vz = vz - dd*ww[ii];
          break;
        case 2:                   /* Right vertical mirror */
          nlen2 = l_seg*l_seg + ww[ii]*ww[ii];
          q = V2Q*(-2)*vdotn_v2/sqrt(nlen2);
          dd = 2*vdotn_v2/nlen2;
          vx = vx + dd*l_seg;
          vz = vz - dd*ww[ii];
          break;
        case 3:                   /* Lower horizontal mirror */
          nlen2 = l_seg*l_seg + hh[ii]*hh[ii];
          q = V2Q*(-2)*vdotn_h1/sqrt(nlen2);
          dd = 2*vdotn_h1/nlen2;
          vy = vy - dd*l_seg;
          vz = vz - dd*hh[ii];
          break;
        case 4:                   /* Upper horizontal mirror */
          nlen2 = l_seg*l_seg + hh[ii]*hh[ii];
          q = V2Q*(-2)*vdotn_h2/sqrt(nlen2);
          dd = 2*vdotn_h2/nlen2;
          vy = vy + dd*l_seg;
          vz = vz - dd*hh[ii];
          break;
      }
      /* Now compute reflectivity. */
      if((i <= 2 && mx == 0) || (i > 2 && my == 0))
      {
        x += hadj; /* Re-adjust origin */
        ABSORB;
      } else {
        double ref=1;
        if (i <= 2)
        {
          double m     = (mx > 0 ? mx : fabs(mx*w1/w1_in[ii]));
          double par[] = {R0, Qcx, alphax, m, W};
          StdReflecFunc(q, par, &ref);
          if (ref > 0)
            p *= ref;
          else {
            x += hadj; /* Re-adjust origin */
            ABSORB;                               /* Cutoff ~ 1E-10 */
          }
        } else {
          double m     = (my > 0 ? my : fabs(my*h1/h1_in[ii]));
          double par[] = {R0, Qcy, alphay, m, W};
          StdReflecFunc(q, par, &ref);
          if (ref > 0)
            p *= ref;
          else {
            x += hadj; /* Re-adjust origin */
            ABSORB;                               /* Cutoff ~ 1E-10 */
          }
        }
      }
      x += hadj; SCATTER; x -= hadj;
    } /* loop on reflections inside segment */
    x += hadj; /* Re-adjust origin */

    /* rotate neutron according to actual guide curvature */
    if (rotation_h) {
      double nvx, nvy, nvz;
      rotate(nvx,nvy,nvz, vx,vy,vz, -rotation_h, 0,1,0);
      vx = nvx; vy=nvy; vz=nvz;
    }
    if (rotation_v) {
      double nvx, nvy, nvz;
      rotate(nvx,nvy,nvz, vx,vy,vz, -rotation_v, 1,0,0);
      vx = nvx; vy=nvy; vz=nvz;
    }
  } /* loop on segments */

  #undef option
  #undef w1
  #undef h1
  #undef l
  #undef linw
  #undef loutw
  #undef linh
  #undef louth
  #undef R0
  #undef Qcx
  #undef Qcy
  #undef alphax
  #undef alphay
  #undef W
  #undef mx
  #undef my
  #undef segno
  #undef curvature
  #undef curvature_v
  #undef w1c
  #undef w2c
  #undef ww
  #undef hh
  #undef whalf
  #undef hhalf
  #undef lwhalf
  #undef lhhalf
  #undef h1_in
  #undef h2_out
  #undef w1_in
  #undef w2_out
  #undef l_seg
  #undef seg
  #undef h12
  #undef h2
  #undef w12
  #undef w2
  #undef a_ell_q
  #undef b_ell_q
  #undef lbw
  #undef lbh
  #undef mxi
  #undef u1
  #undef u2
  #undef div1
  #undef p2_para
  #undef test
  #undef Div1
  #undef seg
  #undef fu
  #undef pos
  #undef file_name
  #undef ep
  #undef num
  #undef rotation_h
  #undef rotation_v
  return(_comp);
} /* class_Guide_tapering_trace */

#pragma acc routine seq
_class_Multilayer_Sample *class_Multilayer_Sample_trace(_class_Multilayer_Sample *_comp
  , _class_particle *_particle) {
  ABSORBED=SCATTERED=RESTORE=0;

  #define xwidth (_comp->_parameters.xwidth)
  #define zlength (_comp->_parameters.zlength)
  #define nlayer (_comp->_parameters.nlayer)
  #define sldPar (_comp->_parameters.sldPar)
  #define dPar (_comp->_parameters.dPar)
  #define sigmaPar (_comp->_parameters.sigmaPar)
  #define frac_inc (_comp->_parameters.frac_inc)
  #define ythick (_comp->_parameters.ythick)
  #define mu_inc (_comp->_parameters.mu_inc)
  #define target_index (_comp->_parameters.target_index)
  #define focus_xw (_comp->_parameters.focus_xw)
  #define focus_yh (_comp->_parameters.focus_yh)
  #define tx (_comp->_parameters.tx)
  #define ty (_comp->_parameters.ty)
  #define tz (_comp->_parameters.tz)
  #define xmin (_comp->_parameters.xmin)
  #define xmax (_comp->_parameters.xmax)
  #define zmin (_comp->_parameters.zmin)
  #define zmax (_comp->_parameters.zmax)
  double dt,q,n1,n2,pfn,R0,lambda,theta,s0,c0,kx,ky,kz,t0,t1,t2,l_i,l_o,v, solid_angle;
  double intersect = 0;

  /* First check if neutron has the right direction. */
  /* calculate time to reach the mirror i.e. y=0*/
  if(vy != 0.0 && (dt = -y/vy) >= 0)
  {
    double old_x = x, old_y = y, old_z = z;
    //printf("x %g y %g z %g vx %g vy %g vz %g \n",x,y,z,vx,vy,vz);
    x += vx*dt;
    //y += vy*dt;
    z += vz*dt;
    y=0;
    /* Now check if neutron intersects mirror. */
    if(x >= xmin && x <= xmax && z >= zmin && z <= zmax)
    {
      // Incoherent scattering from substrate or coherent scattering from thin film?
      if (rand01()<frac_inc) {
	RESTORE_NEUTRON(INDEX_CURRENT_COMP, x, y, z, vx, vy, vz, t, sx, sy, sz, p);
	// This part is basically from V_sample

	intersect = box_intersect(&t0, &t2, x, y+ythick/2.0, z, vx, vy, vz, xwidth, ythick, zlength);
	if (intersect) {
	  if (t0 < 0) ABSORB; /* we already passed the sample; this is illegal */
	  dt = rand01()*(t2-t0);/* Time of scattering (relative to t0) */
	  PROP_DT(dt+t0);
	  SCATTER;
	  v = sqrt(vx*vx + vy*vy + vz*vz);
	  l_i = v*dt;                 /* Penetration in sample until scattering */

	  // If target comp and focus area set, work with that. Otherwise scatter in 4PI
	  if (target_index && (focus_xw>0) && (focus_yh>0)) {
	    randvec_target_rect(&vx, &vy, &vz, &solid_angle, tx, ty, tz, focus_xw, focus_yh, ROT_A_CURRENT_COMP);
	  } else {
	    if (tx == ty == tz == 0) {
	      ty = 1e-9;
	    }
	    randvec_target_circle(&vx, &vy, &vz, &solid_angle, tx, ty, tz, 0);
	  }
	  NORM(vx, vy, vz);
	  vx *= v;
	  vy *= v;
	  vz *= v;
	  intersect = box_intersect(&t0, &t2, x, y, z, vx, vy, vz, xwidth, ythick, zlength);
	  if (intersect) {
	    l_o = v*t2;
	    p *= (l_i+l_o)*(mu_inc/100.0)*exp(mu_inc*(l_i+l_o)/100.0);
	    p /= 4*PI/solid_angle;
	    p /= frac_inc;
	  } else {
	    // Kill neutron!
	    printf("Could not hit sample from inside. ABSORBED\n");
	    ABSORB;
	  }
	} else { // Otherwise simply leave alone
	  RESTORE_NEUTRON(INDEX_CURRENT_COMP, x, y, z, vx, vy, vz, t, sx, sy, sz, p);
	}

      } else {
//
// If neutron intersect mirror calculate reflectivity from a thin
// film using simple fresnel coefficients formula found in e.g. Born and Wolf
// this could be generalised to many layers using the matrix formalism
// but we'll get this bit to work first.
//
      double dbl0,dbl1,tvar;
      int arrsize=nlayer+2;
      int i;

      gsl_complex rnf,rnf1;
      gsl_complex a12t,a22t,cr,c0,ci;
      gsl_complex btm,btm1,cbtm,cbtm1;
      gsl_complex ac1,ac2,ac3,ac4;
      gsl_complex cans,tcvar1,tcvar2,tcvar3,tcvar4;

      gsl_matrix_complex * ris = gsl_matrix_complex_alloc(500,1);
      gsl_matrix_complex * pfn = gsl_matrix_complex_alloc(500,1);
      gsl_matrix_complex * betan = gsl_matrix_complex_alloc(500,1);
      gsl_matrix_complex * a1 = gsl_matrix_complex_alloc(2,2);
      gsl_matrix_complex * a2 = gsl_matrix_complex_alloc(2,2);
      gsl_matrix_complex * a3 = gsl_matrix_complex_alloc(2,2);

      t += dt;
      //q = fabs(2*vy*V2Q);
      kx = vx*V2K;
      ky = vy*V2K;
      kz = vz*V2K;
      lambda = 2*PI/sqrt(kx*kx+ky*ky+kz*kz);
      theta = atan(fabs(old_y)/fabs(z-old_z));

      double lsq=lambda*lambda;
      double tpi=2*PI;
      double tlc=8.0*PI*PI/lsq;

      double st0 = sin(theta);
      double ct0 = cos(theta);
      dbl0=0.0;
      dbl1=1.0;
      cr=gsl_complex_rect(dbl1,dbl0);
      ci=gsl_complex_rect(dbl0,dbl1);
      c0=gsl_complex_rect(dbl0,dbl0);
      a12t=c0;
      a22t=c0;

      for(i=0; i< nlayer+2; i++)
      {
        tcvar1=gsl_complex_rect(1.0 - (lsq * sldPar[i] / tpi),dbl0);
	gsl_matrix_complex_set(ris,i,0,tcvar1);
      }

      gsl_matrix_complex_set(pfn,0,0,gsl_complex_mul_real(gsl_matrix_complex_get(ris,0,0),st0));
      rnf1=gsl_complex_mul(gsl_matrix_complex_get(ris,0,0),gsl_matrix_complex_get(ris,0,0));
      if(nlayer > 0){
            for(i=1;i<nlayer+1;i++){
	            rnf=gsl_complex_mul(gsl_matrix_complex_get(ris,i,0),gsl_matrix_complex_get(ris,i,0));
		    tcvar1=gsl_complex_sub(rnf,gsl_complex_mul_real(rnf1,ct0*ct0));
	            gsl_matrix_complex_set(pfn,i,0,gsl_complex_sqrt(tcvar1));
            }
      }
      tcvar1=gsl_matrix_complex_get(ris,nlayer+1,0);
      rnf=gsl_complex_mul(tcvar1,tcvar1);
      tcvar1=gsl_complex_sub(rnf,gsl_complex_mul_real(rnf1,ct0*ct0));
      gsl_matrix_complex_set(pfn,nlayer + 1,0,gsl_complex_sqrt(tcvar1));

      if(nlayer > 0){
            for(i=0;i<nlayer;i++){
		tcvar1=gsl_matrix_complex_get(pfn,i+1,0);
		gsl_matrix_complex_set(betan,i+1,0,gsl_complex_mul_real(tcvar1,tpi*dPar[i]/lambda));
	    }
      }

      gsl_matrix_complex_set(a1,0,0,cr);
      tcvar1=gsl_matrix_complex_get(pfn,0,0);
      tcvar2=gsl_matrix_complex_get(pfn,1,0);
      if(GSL_REAL(gsl_complex_add(tcvar1,tcvar2))!=0.0 || GSL_IMAG(gsl_complex_add(tcvar1,tcvar2))!=0.0){
	a12t=gsl_complex_div(gsl_complex_sub(tcvar1,tcvar2),gsl_complex_add(tcvar1,tcvar2));
      }else{
	a12t=c0;
      }
      tcvar3=gsl_complex_mul(tcvar1,tcvar2);
      tcvar4=gsl_complex_mul_real(tcvar3,-1.0*tlc*sigmaPar[0]*sigmaPar[0]);
      gsl_matrix_complex_set(a1,0,1,gsl_complex_mul(a12t,gsl_complex_exp(tcvar4)));
      gsl_matrix_complex_set(a1,1,0,gsl_matrix_complex_get(a1,0,1));
      gsl_matrix_complex_set(a1,1,1,cr);

      if(nlayer > 0){
            for(i=1;i<nlayer+1;i++){
		btm=gsl_complex_mul(gsl_matrix_complex_get(betan,i,0),ci);
                btm1=gsl_complex_mul(gsl_complex_mul_real(gsl_matrix_complex_get(betan,i,0),-1.0),ci);
	        cbtm=gsl_complex_exp(btm);
                cbtm1=gsl_complex_exp(btm1);
                gsl_matrix_complex_set(a2,0,0,cbtm);
        	tcvar1=gsl_matrix_complex_get(pfn,i,0);
      		tcvar2=gsl_matrix_complex_get(pfn,i+1,0);
      		if(GSL_REAL(gsl_complex_add(tcvar1,tcvar2))!=0.0 || GSL_IMAG(gsl_complex_add(tcvar1,tcvar2))!=0.0){
		  a22t=gsl_complex_div(gsl_complex_sub(tcvar1,tcvar2),gsl_complex_add(tcvar1,tcvar2));
      		}else{
		  a22t=c0;
		}
		tcvar3=gsl_complex_mul(tcvar1,tcvar2);
      		tcvar4=gsl_complex_mul_real(tcvar3,-1.0*tlc*sigmaPar[i]*sigmaPar[i]);
		a22t=gsl_complex_mul(a22t,gsl_complex_exp(tcvar4));
      		gsl_matrix_complex_set(a2,0,1,gsl_complex_mul(a22t,cbtm));
      		gsl_matrix_complex_set(a2,1,0,gsl_complex_mul(a22t,cbtm1));
                gsl_matrix_complex_set(a2,1,1,cbtm1);

		tcvar1=gsl_matrix_complex_get(a1,0,0);
		tcvar2=gsl_complex_mul(tcvar1,gsl_matrix_complex_get(a2,0,0));
		tcvar3=gsl_matrix_complex_get(a1,0,1);
		tcvar4=gsl_complex_mul(tcvar3,gsl_matrix_complex_get(a2,1,0));
                gsl_matrix_complex_set(a3,0,0,gsl_complex_add(tcvar2,tcvar4));

		tcvar1=gsl_matrix_complex_get(a1,0,0);
		tcvar2=gsl_complex_mul(tcvar1,gsl_matrix_complex_get(a2,0,1));
		tcvar3=gsl_matrix_complex_get(a1,0,1);
		tcvar4=gsl_complex_mul(tcvar3,gsl_matrix_complex_get(a2,1,1));
                gsl_matrix_complex_set(a3,0,1,gsl_complex_add(tcvar2,tcvar4));

       		tcvar1=gsl_matrix_complex_get(a1,1,0);
		tcvar2=gsl_complex_mul(tcvar1,gsl_matrix_complex_get(a2,0,0));
		tcvar3=gsl_matrix_complex_get(a1,1,1);
		tcvar4=gsl_complex_mul(tcvar3,gsl_matrix_complex_get(a2,1,0));
                gsl_matrix_complex_set(a3,1,0,gsl_complex_add(tcvar2,tcvar4));

		tcvar1=gsl_matrix_complex_get(a1,1,0);
		tcvar2=gsl_complex_mul(tcvar1,gsl_matrix_complex_get(a2,0,1));
		tcvar3=gsl_matrix_complex_get(a1,1,1);
		tcvar4=gsl_complex_mul(tcvar3,gsl_matrix_complex_get(a2,1,1));
                gsl_matrix_complex_set(a3,1,1,gsl_complex_add(tcvar2,tcvar4));

		gsl_matrix_complex_set(a1,0,0,gsl_matrix_complex_get(a3,0,0));
		gsl_matrix_complex_set(a1,0,1,gsl_matrix_complex_get(a3,0,1));
		gsl_matrix_complex_set(a1,1,0,gsl_matrix_complex_get(a3,1,0));
		gsl_matrix_complex_set(a1,1,1,gsl_matrix_complex_get(a3,1,1));
            }
      }
      ac1=gsl_matrix_complex_get(a1,1,0);
      ac2=gsl_complex_conjugate(ac1);
      ac3=gsl_matrix_complex_get(a1,0,0);
      ac4=gsl_complex_conjugate(ac3);
      cans=gsl_complex_div(gsl_complex_mul(ac1,ac2),gsl_complex_mul(ac3,ac4));
      R0 = gsl_complex_abs(cans);

      //printf("Q %g pfn %g lambda %g \n",q,pfn,lambda);
      /*reflect off horizontal surface so reverse y component of velocity*/
      vy = -vy;
      /* Reflectivity (see component Guide). */
      p *= R0;
      if (frac_inc>0) {
	p /= (1-frac_inc);
      }
      SCATTER;

      gsl_matrix_complex_free(ris);
      gsl_matrix_complex_free(pfn);
      gsl_matrix_complex_free(betan);
      gsl_matrix_complex_free(a1);
      gsl_matrix_complex_free(a2);
      gsl_matrix_complex_free(a3);
      }
    }
    else
    {
      /* No intersection: restore neutron state. */
      RESTORE_NEUTRON(INDEX_CURRENT_COMP, x, y, z, vx, vy, vz, t, sx, sy, sz, p);
    }
  }
  #undef xwidth
  #undef zlength
  #undef nlayer
  #undef sldPar
  #undef dPar
  #undef sigmaPar
  #undef frac_inc
  #undef ythick
  #undef mu_inc
  #undef target_index
  #undef focus_xw
  #undef focus_yh
  #undef tx
  #undef ty
  #undef tz
  #undef xmin
  #undef xmax
  #undef zmin
  #undef zmax
  return(_comp);
} /* class_Multilayer_Sample_trace */

/* *****************************************************************************
* instrument 'ISIS_CRISP' TRACE
***************************************************************************** */

#pragma acc routine seq
int raytrace(_class_particle* _particle) { /* called by mccode_main for ISIS_CRISP:TRACE */

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
      /* component a1=Arm() [1] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_a1->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _a1->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_a1->_position_relative, _a1->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      _particle->_index++;
    } /* end component a1 [1] */
    if (!ABSORBED && _particle->_index == 2) {
      /* component isis_source=ISIS_moderator() [2] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_isis_source->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _isis_source->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_isis_source->_position_relative, _isis_source->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_ISIS_moderator_trace(_isis_source, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component isis_source [2] */
    if (!ABSORBED && _particle->_index == 3) {
      /* component defslit1=Slit() [3] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_defslit1->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _defslit1->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_defslit1->_position_relative, _defslit1->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Slit_trace(_defslit1, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component defslit1 [3] */
    if (!ABSORBED && _particle->_index == 4) {
      /* component defslit2=Slit() [4] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_defslit2->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _defslit2->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_defslit2->_position_relative, _defslit2->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Slit_trace(_defslit2, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component defslit2 [4] */
    if (!ABSORBED && _particle->_index == 5) {
      /* component coarseslit1=Slit() [5] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_coarseslit1->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _coarseslit1->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_coarseslit1->_position_relative, _coarseslit1->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Slit_trace(_coarseslit1, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component coarseslit1 [5] */
    if (!ABSORBED && _particle->_index == 6) {
      /* component coarseslit2=Slit() [6] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_coarseslit2->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _coarseslit2->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_coarseslit2->_position_relative, _coarseslit2->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Slit_trace(_coarseslit2, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component coarseslit2 [6] */
    if (!ABSORBED && _particle->_index == 7) {
      /* component lmon1=L_monitor() [7] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_lmon1->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _lmon1->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_lmon1->_position_relative, _lmon1->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_L_monitor_trace(_lmon1, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component lmon1 [7] */
    if (!ABSORBED && _particle->_index == 8) {
      /* component PSDmon1=PSD_monitor() [8] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_PSDmon1->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _PSDmon1->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_PSDmon1->_position_relative, _PSDmon1->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_PSD_monitor_trace(_PSDmon1, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component PSDmon1 [8] */
    if (!ABSORBED && _particle->_index == 9) {
      /* component slit1=Slit() [9] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_slit1->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _slit1->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_slit1->_position_relative, _slit1->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Slit_trace(_slit1, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component slit1 [9] */
    if (!ABSORBED && _particle->_index == 10) {
      /* component PSDmons1=PSD_monitor() [10] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_PSDmons1->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _PSDmons1->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_PSDmons1->_position_relative, _PSDmons1->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_PSD_monitor_trace(_PSDmons1, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component PSDmons1 [10] */
    if (!ABSORBED && _particle->_index == 11) {
      /* component eguide1=Guide_tapering() [11] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_eguide1->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _eguide1->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_eguide1->_position_relative, _eguide1->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Guide_tapering_trace(_eguide1, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component eguide1 [11] */
    if (!ABSORBED && _particle->_index == 12) {
      /* component slit2=Slit() [12] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_slit2->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _slit2->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_slit2->_position_relative, _slit2->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Slit_trace(_slit2, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component slit2 [12] */
    if (!ABSORBED && _particle->_index == 13) {
      /* component lmon2=L_monitor() [13] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_lmon2->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _lmon2->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_lmon2->_position_relative, _lmon2->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_L_monitor_trace(_lmon2, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component lmon2 [13] */
    if (!ABSORBED && _particle->_index == 14) {
      /* component samp1=Multilayer_Sample() [14] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_samp1->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _samp1->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_samp1->_position_relative, _samp1->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Multilayer_Sample_trace(_samp1, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component samp1 [14] */
    if (!ABSORBED && _particle->_index == 15) {
      /* component slit3=Slit() [15] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_slit3->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _slit3->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_slit3->_position_relative, _slit3->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Slit_trace(_slit3, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component slit3 [15] */
    if (!ABSORBED && _particle->_index == 16) {
      /* component eguide2=Guide_tapering() [16] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_eguide2->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _eguide2->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_eguide2->_position_relative, _eguide2->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Guide_tapering_trace(_eguide2, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component eguide2 [16] */
    if (!ABSORBED && _particle->_index == 17) {
      /* component lmon3=L_monitor() [17] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_lmon3->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _lmon3->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_lmon3->_position_relative, _lmon3->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_L_monitor_trace(_lmon3, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component lmon3 [17] */
    if (!ABSORBED && _particle->_index == 18) {
      /* component PSDdet1=PSD_monitor() [18] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_PSDdet1->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _PSDdet1->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_PSDdet1->_position_relative, _PSDdet1->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_PSD_monitor_trace(_PSDdet1, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component PSDdet1 [18] */
    if (_particle->_index > 18)
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
* instrument 'ISIS_CRISP' and components SAVE
***************************************************************************** */

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



int save(FILE *handle) { /* called by mccode_main for ISIS_CRISP:SAVE */
  if (!handle) siminfo_init(NULL);

  /* call iteratively all components SAVE */






  class_L_monitor_save(_lmon1);

  class_PSD_monitor_save(_PSDmon1);


  class_PSD_monitor_save(_PSDmons1);



  class_L_monitor_save(_lmon2);




  class_L_monitor_save(_lmon3);

  class_PSD_monitor_save(_PSDdet1);

  if (!handle) siminfo_close(); 

  return(0);
} /* save */

/* *****************************************************************************
* instrument 'ISIS_CRISP' and components FINALLY
***************************************************************************** */

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

_class_Guide_tapering *class_Guide_tapering_finally(_class_Guide_tapering *_comp
) {
  #define option (_comp->_parameters.option)
  #define w1 (_comp->_parameters.w1)
  #define h1 (_comp->_parameters.h1)
  #define l (_comp->_parameters.l)
  #define linw (_comp->_parameters.linw)
  #define loutw (_comp->_parameters.loutw)
  #define linh (_comp->_parameters.linh)
  #define louth (_comp->_parameters.louth)
  #define R0 (_comp->_parameters.R0)
  #define Qcx (_comp->_parameters.Qcx)
  #define Qcy (_comp->_parameters.Qcy)
  #define alphax (_comp->_parameters.alphax)
  #define alphay (_comp->_parameters.alphay)
  #define W (_comp->_parameters.W)
  #define mx (_comp->_parameters.mx)
  #define my (_comp->_parameters.my)
  #define segno (_comp->_parameters.segno)
  #define curvature (_comp->_parameters.curvature)
  #define curvature_v (_comp->_parameters.curvature_v)
  #define w1c (_comp->_parameters.w1c)
  #define w2c (_comp->_parameters.w2c)
  #define ww (_comp->_parameters.ww)
  #define hh (_comp->_parameters.hh)
  #define whalf (_comp->_parameters.whalf)
  #define hhalf (_comp->_parameters.hhalf)
  #define lwhalf (_comp->_parameters.lwhalf)
  #define lhhalf (_comp->_parameters.lhhalf)
  #define h1_in (_comp->_parameters.h1_in)
  #define h2_out (_comp->_parameters.h2_out)
  #define w1_in (_comp->_parameters.w1_in)
  #define w2_out (_comp->_parameters.w2_out)
  #define l_seg (_comp->_parameters.l_seg)
  #define seg (_comp->_parameters.seg)
  #define h12 (_comp->_parameters.h12)
  #define h2 (_comp->_parameters.h2)
  #define w12 (_comp->_parameters.w12)
  #define w2 (_comp->_parameters.w2)
  #define a_ell_q (_comp->_parameters.a_ell_q)
  #define b_ell_q (_comp->_parameters.b_ell_q)
  #define lbw (_comp->_parameters.lbw)
  #define lbh (_comp->_parameters.lbh)
  #define mxi (_comp->_parameters.mxi)
  #define u1 (_comp->_parameters.u1)
  #define u2 (_comp->_parameters.u2)
  #define div1 (_comp->_parameters.div1)
  #define p2_para (_comp->_parameters.p2_para)
  #define test (_comp->_parameters.test)
  #define Div1 (_comp->_parameters.Div1)
  #define seg (_comp->_parameters.seg)
  #define fu (_comp->_parameters.fu)
  #define pos (_comp->_parameters.pos)
  #define file_name (_comp->_parameters.file_name)
  #define ep (_comp->_parameters.ep)
  #define num (_comp->_parameters.num)
  #define rotation_h (_comp->_parameters.rotation_h)
  #define rotation_v (_comp->_parameters.rotation_v)
  free(w1c);
  free(w2c);
  free(ww);
  free(hh);
  free(whalf);
  free(hhalf);
  free(lwhalf);
  free(lhhalf);
  free(h1_in);
  free(h2_out);
  free(w1_in);
  free(w2_out);
  #undef option
  #undef w1
  #undef h1
  #undef l
  #undef linw
  #undef loutw
  #undef linh
  #undef louth
  #undef R0
  #undef Qcx
  #undef Qcy
  #undef alphax
  #undef alphay
  #undef W
  #undef mx
  #undef my
  #undef segno
  #undef curvature
  #undef curvature_v
  #undef w1c
  #undef w2c
  #undef ww
  #undef hh
  #undef whalf
  #undef hhalf
  #undef lwhalf
  #undef lhhalf
  #undef h1_in
  #undef h2_out
  #undef w1_in
  #undef w2_out
  #undef l_seg
  #undef seg
  #undef h12
  #undef h2
  #undef w12
  #undef w2
  #undef a_ell_q
  #undef b_ell_q
  #undef lbw
  #undef lbh
  #undef mxi
  #undef u1
  #undef u2
  #undef div1
  #undef p2_para
  #undef test
  #undef Div1
  #undef seg
  #undef fu
  #undef pos
  #undef file_name
  #undef ep
  #undef num
  #undef rotation_h
  #undef rotation_v
  return(_comp);
} /* class_Guide_tapering_finally */



int finally(void) { /* called by mccode_main for ISIS_CRISP:FINALLY */
  siminfo_init(NULL);
  save(siminfo_file); /* save data when simulation ends */

  /* call iteratively all components FINALLY */






  class_L_monitor_finally(_lmon1);

  class_PSD_monitor_finally(_PSDmon1);


  class_PSD_monitor_finally(_PSDmons1);

  class_Guide_tapering_finally(_eguide1);


  class_L_monitor_finally(_lmon2);



  class_Guide_tapering_finally(_eguide2);

  class_L_monitor_finally(_lmon3);

  class_PSD_monitor_finally(_PSDdet1);

  siminfo_close(); 

  return(0);
} /* finally */

/* *****************************************************************************
* instrument 'ISIS_CRISP' and components DISPLAY
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
_class_Arm *class_Arm_display(_class_Arm *_comp
) {
  /* A bit ugly; hard-coded dimensions. */
  
  line(0,0,0,0.2,0,0);
  line(0,0,0,0,0.2,0);
  line(0,0,0,0,0,0.2);
  return(_comp);
} /* class_Arm_display */

_class_ISIS_moderator *class_ISIS_moderator_display(_class_ISIS_moderator *_comp
) {
  #define Face (_comp->_parameters.Face)
  #define Emin (_comp->_parameters.Emin)
  #define Emax (_comp->_parameters.Emax)
  #define dist (_comp->_parameters.dist)
  #define focus_xw (_comp->_parameters.focus_xw)
  #define focus_yh (_comp->_parameters.focus_yh)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define CAngle (_comp->_parameters.CAngle)
  #define SAC (_comp->_parameters.SAC)
  #define Lmin (_comp->_parameters.Lmin)
  #define Lmax (_comp->_parameters.Lmax)
  #define target_index (_comp->_parameters.target_index)
  #define verbose (_comp->_parameters.verbose)
  #define p_in (_comp->_parameters.p_in)
  #define Tnpts (_comp->_parameters.Tnpts)
  #define scaleSize (_comp->_parameters.scaleSize)
  #define angleArea (_comp->_parameters.angleArea)
  #define Nsim (_comp->_parameters.Nsim)
  #define Ncount (_comp->_parameters.Ncount)
  #define TS (_comp->_parameters.TS)
  #define rtE0 (_comp->_parameters.rtE0)
  #define rtE1 (_comp->_parameters.rtE1)
  #define rtmodX (_comp->_parameters.rtmodX)
  #define rtmodY (_comp->_parameters.rtmodY)
  #define TargetStation (_comp->_parameters.TargetStation)
  #define CurrentWeight (_comp->_parameters.CurrentWeight)
  double cirp=0.0,cirq=0.3,pi=3.141592654;
  int pp=0; /* circle drawing parameter*/



  
  multiline(5,-0.5*rtmodX,-0.5*rtmodY,0.0,
	    0.5*rtmodX,-0.5*rtmodY,0.0,
	    0.5*rtmodX,0.5*rtmodY,0.0,
	    -0.5*rtmodX,0.5*rtmodY,0.0,
	    -0.5*rtmodX,-0.5*rtmodY,0.0);
  /* circle("xy",0.0,0.0,0.0,cos(cirp)); */

  /*line(0.5*sin(cirp),0.0,0.5*cos(cirp),0.5*sin(cirq),0.0,0.5*cos(cirq));*/

  /*line(-0.5,0.0,0.0,0.0,0.0,0.5);*/

  for (pp=0;pp<=20;pp=pp+2)
    {
      cirp= (pp*(pi/21.0))-(0.5*pi);
      cirq= ((pp+1)*(pi/21.0))-(0.5*pi);
      line(0.5*sin(cirp),0.0,0.5*cos(cirp),0.5*sin(cirq),0.0,0.5*cos(cirq));
    }

  #undef Face
  #undef Emin
  #undef Emax
  #undef dist
  #undef focus_xw
  #undef focus_yh
  #undef xwidth
  #undef yheight
  #undef CAngle
  #undef SAC
  #undef Lmin
  #undef Lmax
  #undef target_index
  #undef verbose
  #undef p_in
  #undef Tnpts
  #undef scaleSize
  #undef angleArea
  #undef Nsim
  #undef Ncount
  #undef TS
  #undef rtE0
  #undef rtE1
  #undef rtmodX
  #undef rtmodY
  #undef TargetStation
  #undef CurrentWeight
  return(_comp);
} /* class_ISIS_moderator_display */

_class_Slit *class_Slit_display(_class_Slit *_comp
) {
  #define xmin (_comp->_parameters.xmin)
  #define xmax (_comp->_parameters.xmax)
  #define ymin (_comp->_parameters.ymin)
  #define ymax (_comp->_parameters.ymax)
  #define radius (_comp->_parameters.radius)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  
  if (radius == 0) {
    double xw, yh;
    xw = (xmax - xmin)/2.0;
    yh = (ymax - ymin)/2.0;
    multiline(3, xmin-xw, (double)ymax, 0.0,
              (double)xmin, (double)ymax, 0.0,
              (double)xmin, ymax+yh, 0.0);
    multiline(3, xmax+xw, (double)ymax, 0.0,
              (double)xmax, (double)ymax, 0.0,
              (double)xmax, ymax+yh, 0.0);
    multiline(3, xmin-xw, (double)ymin, 0.0,
              (double)xmin, (double)ymin, 0.0,
              (double)xmin, ymin-yh, 0.0);
    multiline(3, xmax+xw, (double)ymin, 0.0,
              (double)xmax, (double)ymin, 0.0,
              (double)xmax, ymin-yh, 0.0);
  } else {
    circle("xy",0,0,0,radius);
  }
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef radius
  #undef xwidth
  #undef yheight
  return(_comp);
} /* class_Slit_display */

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

_class_Guide_tapering *class_Guide_tapering_display(_class_Guide_tapering *_comp
) {
  #define option (_comp->_parameters.option)
  #define w1 (_comp->_parameters.w1)
  #define h1 (_comp->_parameters.h1)
  #define l (_comp->_parameters.l)
  #define linw (_comp->_parameters.linw)
  #define loutw (_comp->_parameters.loutw)
  #define linh (_comp->_parameters.linh)
  #define louth (_comp->_parameters.louth)
  #define R0 (_comp->_parameters.R0)
  #define Qcx (_comp->_parameters.Qcx)
  #define Qcy (_comp->_parameters.Qcy)
  #define alphax (_comp->_parameters.alphax)
  #define alphay (_comp->_parameters.alphay)
  #define W (_comp->_parameters.W)
  #define mx (_comp->_parameters.mx)
  #define my (_comp->_parameters.my)
  #define segno (_comp->_parameters.segno)
  #define curvature (_comp->_parameters.curvature)
  #define curvature_v (_comp->_parameters.curvature_v)
  #define w1c (_comp->_parameters.w1c)
  #define w2c (_comp->_parameters.w2c)
  #define ww (_comp->_parameters.ww)
  #define hh (_comp->_parameters.hh)
  #define whalf (_comp->_parameters.whalf)
  #define hhalf (_comp->_parameters.hhalf)
  #define lwhalf (_comp->_parameters.lwhalf)
  #define lhhalf (_comp->_parameters.lhhalf)
  #define h1_in (_comp->_parameters.h1_in)
  #define h2_out (_comp->_parameters.h2_out)
  #define w1_in (_comp->_parameters.w1_in)
  #define w2_out (_comp->_parameters.w2_out)
  #define l_seg (_comp->_parameters.l_seg)
  #define seg (_comp->_parameters.seg)
  #define h12 (_comp->_parameters.h12)
  #define h2 (_comp->_parameters.h2)
  #define w12 (_comp->_parameters.w12)
  #define w2 (_comp->_parameters.w2)
  #define a_ell_q (_comp->_parameters.a_ell_q)
  #define b_ell_q (_comp->_parameters.b_ell_q)
  #define lbw (_comp->_parameters.lbw)
  #define lbh (_comp->_parameters.lbh)
  #define mxi (_comp->_parameters.mxi)
  #define u1 (_comp->_parameters.u1)
  #define u2 (_comp->_parameters.u2)
  #define div1 (_comp->_parameters.div1)
  #define p2_para (_comp->_parameters.p2_para)
  #define test (_comp->_parameters.test)
  #define Div1 (_comp->_parameters.Div1)
  #define seg (_comp->_parameters.seg)
  #define fu (_comp->_parameters.fu)
  #define pos (_comp->_parameters.pos)
  #define file_name (_comp->_parameters.file_name)
  #define ep (_comp->_parameters.ep)
  #define num (_comp->_parameters.num)
  #define rotation_h (_comp->_parameters.rotation_h)
  #define rotation_v (_comp->_parameters.rotation_v)
  int i,ii;

  for (ii=0; ii < segno; ii++)
  {
     multiline(5,
        -w1_in[ii]/2.0, -h1_in[ii]/2.0,l_seg*(double)ii,
        -w2_out[ii]/2.0, -h2_out[ii]/2.0,l_seg*((double)ii+1.0),
        -w2_out[ii]/2.0,  h2_out[ii]/2.0,l_seg*((double)ii+1.0),
        -w1_in[ii]/2.0,  h1_in[ii]/2.0,l_seg*(double)ii,
        -w1_in[ii]/2.0, -h1_in[ii]/2.0,l_seg*(double)ii);
     multiline(5,
        w1_in[ii]/2.0, -h1_in[ii]/2.0,l_seg*(double)ii,
        w2_out[ii]/2.0, -h2_out[ii]/2.0,l_seg*((double)ii+1.0),
        w2_out[ii]/2.0,  h2_out[ii]/2.0,l_seg*((double)ii+1.0),
        w1_in[ii]/2.0,  h1_in[ii]/2.0,l_seg*(double)ii,
        w1_in[ii]/2.0, -h1_in[ii]/2.0,l_seg*(double)ii);
  }
  line(-w1/2.0, -h1/2.0, 0.0, w1/2.0, -h1/2.0, 0.0);
  line(-w1/2.0, h1/2.0, 0.0, w1/2.0, h1/2.0, 0.0);
  for(i=0; i<segno;i++)
  {
     line(-w2_out[i]/2.0, -h2_out[i]/2.0, l_seg*(double)(i+1),
     w2_out[i]/2.0, -h2_out[i]/2.0, l_seg*(double)(i+1));
     line(-w2_out[i]/2.0, h2_out[i]/2.0, l_seg*(double)(i+1),
     w2_out[i]/2.0, h2_out[i]/2.0, l_seg*(double)(i+1));
  }

  #undef option
  #undef w1
  #undef h1
  #undef l
  #undef linw
  #undef loutw
  #undef linh
  #undef louth
  #undef R0
  #undef Qcx
  #undef Qcy
  #undef alphax
  #undef alphay
  #undef W
  #undef mx
  #undef my
  #undef segno
  #undef curvature
  #undef curvature_v
  #undef w1c
  #undef w2c
  #undef ww
  #undef hh
  #undef whalf
  #undef hhalf
  #undef lwhalf
  #undef lhhalf
  #undef h1_in
  #undef h2_out
  #undef w1_in
  #undef w2_out
  #undef l_seg
  #undef seg
  #undef h12
  #undef h2
  #undef w12
  #undef w2
  #undef a_ell_q
  #undef b_ell_q
  #undef lbw
  #undef lbh
  #undef mxi
  #undef u1
  #undef u2
  #undef div1
  #undef p2_para
  #undef test
  #undef Div1
  #undef seg
  #undef fu
  #undef pos
  #undef file_name
  #undef ep
  #undef num
  #undef rotation_h
  #undef rotation_v
  return(_comp);
} /* class_Guide_tapering_display */

_class_Multilayer_Sample *class_Multilayer_Sample_display(_class_Multilayer_Sample *_comp
) {
  #define xwidth (_comp->_parameters.xwidth)
  #define zlength (_comp->_parameters.zlength)
  #define nlayer (_comp->_parameters.nlayer)
  #define sldPar (_comp->_parameters.sldPar)
  #define dPar (_comp->_parameters.dPar)
  #define sigmaPar (_comp->_parameters.sigmaPar)
  #define frac_inc (_comp->_parameters.frac_inc)
  #define ythick (_comp->_parameters.ythick)
  #define mu_inc (_comp->_parameters.mu_inc)
  #define target_index (_comp->_parameters.target_index)
  #define focus_xw (_comp->_parameters.focus_xw)
  #define focus_yh (_comp->_parameters.focus_yh)
  #define tx (_comp->_parameters.tx)
  #define ty (_comp->_parameters.ty)
  #define tz (_comp->_parameters.tz)
  #define xmin (_comp->_parameters.xmin)
  #define xmax (_comp->_parameters.xmax)
  #define zmin (_comp->_parameters.zmin)
  #define zmax (_comp->_parameters.zmax)
  box(0.0, (double)-ythick/2.0, 0.0, (double)xwidth, (double)ythick, (double)zlength);
  #undef xwidth
  #undef zlength
  #undef nlayer
  #undef sldPar
  #undef dPar
  #undef sigmaPar
  #undef frac_inc
  #undef ythick
  #undef mu_inc
  #undef target_index
  #undef focus_xw
  #undef focus_yh
  #undef tx
  #undef ty
  #undef tz
  #undef xmin
  #undef xmax
  #undef zmin
  #undef zmax
  return(_comp);
} /* class_Multilayer_Sample_display */


  #undef magnify
  #undef line
  #undef dashed_line
  #undef multiline
  #undef rectangle
  #undef box
  #undef circle
  #undef cylinder
  #undef sphere

int display(void) { /* called by mccode_main for ISIS_CRISP:DISPLAY */
  printf("MCDISPLAY: start\n");

  /* call iteratively all components DISPLAY */
  class_Arm_display(_a1);

  class_ISIS_moderator_display(_isis_source);

  class_Slit_display(_defslit1);

  class_Slit_display(_defslit2);

  class_Slit_display(_coarseslit1);

  class_Slit_display(_coarseslit2);

  class_L_monitor_display(_lmon1);

  class_PSD_monitor_display(_PSDmon1);

  class_Slit_display(_slit1);

  class_PSD_monitor_display(_PSDmons1);

  class_Guide_tapering_display(_eguide1);

  class_Slit_display(_slit2);

  class_L_monitor_display(_lmon2);

  class_Multilayer_Sample_display(_samp1);

  class_Slit_display(_slit3);

  class_Guide_tapering_display(_eguide2);

  class_L_monitor_display(_lmon3);

  class_PSD_monitor_display(_PSDdet1);

  printf("MCDISPLAY: end\n");

  return(0);
} /* display */

void* _getvar_parameters(char* compname)
/* enables settings parameters based use of the GETPAR macro */
{
  if (!strcmp(compname, "a1")) return (void *) &(_a1->_parameters);
  if (!strcmp(compname, "isis_source")) return (void *) &(_isis_source->_parameters);
  if (!strcmp(compname, "defslit1")) return (void *) &(_defslit1->_parameters);
  if (!strcmp(compname, "defslit2")) return (void *) &(_defslit2->_parameters);
  if (!strcmp(compname, "coarseslit1")) return (void *) &(_coarseslit1->_parameters);
  if (!strcmp(compname, "coarseslit2")) return (void *) &(_coarseslit2->_parameters);
  if (!strcmp(compname, "lmon1")) return (void *) &(_lmon1->_parameters);
  if (!strcmp(compname, "PSDmon1")) return (void *) &(_PSDmon1->_parameters);
  if (!strcmp(compname, "slit1")) return (void *) &(_slit1->_parameters);
  if (!strcmp(compname, "PSDmons1")) return (void *) &(_PSDmons1->_parameters);
  if (!strcmp(compname, "eguide1")) return (void *) &(_eguide1->_parameters);
  if (!strcmp(compname, "slit2")) return (void *) &(_slit2->_parameters);
  if (!strcmp(compname, "lmon2")) return (void *) &(_lmon2->_parameters);
  if (!strcmp(compname, "samp1")) return (void *) &(_samp1->_parameters);
  if (!strcmp(compname, "slit3")) return (void *) &(_slit3->_parameters);
  if (!strcmp(compname, "eguide2")) return (void *) &(_eguide2->_parameters);
  if (!strcmp(compname, "lmon3")) return (void *) &(_lmon3->_parameters);
  if (!strcmp(compname, "PSDdet1")) return (void *) &(_PSDdet1->_parameters);
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

/* end of generated C code ISIS_CRISP.c */
