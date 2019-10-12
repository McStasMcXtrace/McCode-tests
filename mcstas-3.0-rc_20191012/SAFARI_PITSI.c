/* Automatically generated file. Do not edit. 
 * Format:     ANSI C source code
 * Creator:    McStas <http://www.mcstas.org>
 * Instrument: SAFARI_PITSI.instr (SAFARI_PITSI)
 * Date:       Sat Oct 12 09:04:41 2019
 * File:       SAFARI_PITSI.c
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
* Start of instrument 'SAFARI_PITSI' generated code
***************************************************************************** */

#ifdef MC_TRACE_ENABLED
int traceenabled = 1;
#else
int traceenabled = 0;
#endif
#define MCSTAS "/usr/share/mcstas/3.0-dev/"
int   defaultmain         = 1;
char  instrument_name[]   = "SAFARI_PITSI";
char  instrument_source[] = "SAFARI_PITSI.instr";
char *instrument_exe      = NULL; /* will be set to argv[0] in main */
char  instrument_code[]   = "Instrument SAFARI_PITSI source code SAFARI_PITSI.instr is not embedded in this executable.\n  Use --source option when running McStas.\n";

int main(int argc, char *argv[]){return mccode_main(argc, argv);}

/* *****************************************************************************
* instrument 'SAFARI_PITSI' and components DECLARE
***************************************************************************** */

/* Instrument parameters: structure and a table for the initialisation
   (Used in e.g. inputparse and I/O function (e.g. detector_out) */

struct _struct_instrument_parameters {
  MCNUM _source_lam_min;
  MCNUM _source_lam_max;
  MCNUM _hi_res;
  MCNUM _mono_Si_type;
  MCNUM _mono_mosh;
  MCNUM _mono_mosv;
  MCNUM _mono_dx;
  MCNUM _mono_dy;
  MCNUM _mono_dz;
  MCNUM _mono_takeoff;
  MCNUM _mono_dtilt;
  MCNUM _mono_r_h;
  MCNUM _mono_r_v;
  MCNUM _port_takeoff;
  MCNUM _inc_slit_rot;
  MCNUM _inc_slit_dx;
  MCNUM _inc_slit_to_cor;
  MCNUM _inc_slit_width;
  MCNUM _inc_slit_height;
  MCNUM _inc_slit_sep;
  MCNUM _mono_to_cor;
  MCNUM _sample_dx;
  MCNUM _sample_dy;
  MCNUM _sample_dz;
  MCNUM _sample_dom;
  MCNUM _det_takeoff;
  MCNUM _cor_to_det;
  MCNUM _dangle_interest;
  MCNUM _full_instrument;
};
typedef struct _struct_instrument_parameters _class_instrument_parameters;

/* instrument SPLIT, GROUP and JUMP control logic */
struct instrument_logic_struct {
  long Group_Monochro; /* equals index of scattering comp when in group */
  long Group_Detectors; /* equals index of scattering comp when in group */
  long Split_Sample; /* this is the SPLIT counter decremented down to 0 */
  _class_particle Split_Sample_particle; /* this is the particle to duplicate */
};

struct _instrument_struct {
  char   _name[256]; /* the name of this instrument e.g. 'SAFARI_PITSI' */
/* Counters per component instance */
  double counter_AbsorbProp[55]; /* absorbed events in PROP routines */
  double counter_N[55], counter_P[55], counter_P2[55]; /* event counters after each component instance */
  _class_particle _trajectory[55]; /* current trajectory for STORE/RESTORE */
/* Components position table (absolute and relative coords) */
  Coords _position_relative[55]; /* positions of all components */
  Coords _position_absolute[55];
  _class_instrument_parameters _parameters; /* instrument parameters */
  struct instrument_logic_struct logic; /* instrument logic */
} _instrument_var;
struct _instrument_struct *instrument = & _instrument_var;
#pragma acc declare create ( _instrument_var )
#pragma acc declare create ( instrument )

int numipar = 29;
struct mcinputtable_struct mcinputtable[] = {
  "source_lam_min", &(_instrument_var._parameters._source_lam_min), instr_type_double, "0.5", 
  "source_lam_max", &(_instrument_var._parameters._source_lam_max), instr_type_double, "2.0", 
  "hi_res", &(_instrument_var._parameters._hi_res), instr_type_double, "0", 
  "mono_Si_type", &(_instrument_var._parameters._mono_Si_type), instr_type_double, "551", 
  "mono_mosh", &(_instrument_var._parameters._mono_mosh), instr_type_double, "30", 
  "mono_mosv", &(_instrument_var._parameters._mono_mosv), instr_type_double, "30", 
  "mono_dx", &(_instrument_var._parameters._mono_dx), instr_type_double, "0", 
  "mono_dy", &(_instrument_var._parameters._mono_dy), instr_type_double, "0", 
  "mono_dz", &(_instrument_var._parameters._mono_dz), instr_type_double, "0", 
  "mono_takeoff", &(_instrument_var._parameters._mono_takeoff), instr_type_double, "90.0", 
  "mono_dtilt", &(_instrument_var._parameters._mono_dtilt), instr_type_double, "0", 
  "mono_r_h", &(_instrument_var._parameters._mono_r_h), instr_type_double, "5.0", 
  "mono_r_v", &(_instrument_var._parameters._mono_r_v), instr_type_double, "3.572", 
  "port_takeoff", &(_instrument_var._parameters._port_takeoff), instr_type_double, "90.0", 
  "inc_slit_rot", &(_instrument_var._parameters._inc_slit_rot), instr_type_double, "0", 
  "inc_slit_dx", &(_instrument_var._parameters._inc_slit_dx), instr_type_double, "0", 
  "inc_slit_to_cor", &(_instrument_var._parameters._inc_slit_to_cor), instr_type_double, "0.01", 
  "inc_slit_width", &(_instrument_var._parameters._inc_slit_width), instr_type_double, "0.006", 
  "inc_slit_height", &(_instrument_var._parameters._inc_slit_height), instr_type_double, "0.05", 
  "inc_slit_sep", &(_instrument_var._parameters._inc_slit_sep), instr_type_double, "0", 
  "mono_to_cor", &(_instrument_var._parameters._mono_to_cor), instr_type_double, "2.5", 
  "sample_dx", &(_instrument_var._parameters._sample_dx), instr_type_double, "0", 
  "sample_dy", &(_instrument_var._parameters._sample_dy), instr_type_double, "0", 
  "sample_dz", &(_instrument_var._parameters._sample_dz), instr_type_double, "0", 
  "sample_dom", &(_instrument_var._parameters._sample_dom), instr_type_double, "0", 
  "det_takeoff", &(_instrument_var._parameters._det_takeoff), instr_type_double, "-114.375", 
  "cor_to_det", &(_instrument_var._parameters._cor_to_det), instr_type_double, "1.179", 
  "dangle_interest", &(_instrument_var._parameters._dangle_interest), instr_type_double, "125", 
  "full_instrument", &(_instrument_var._parameters._full_instrument), instr_type_double, "1", 
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


/* Shared user declarations for all components types 'Filter_gen'. */
#ifndef FILTER_GEN
  #define FILTER_GEN $Revision$
  #define UNKNOWN_TABLE    0
  #define ENERGY_TABLE     1
  #define WAVEVECTOR_TABLE 2
  #define WAVELENGTH_TABLE 3
  #define FLUX_ADAPT_SET   0
  #define FLUX_ADAPT_MULT  1
  #define FLUX_ADAPT_ADD   2

  char FilterGen_Mode(char *str, char *Mode, char *Type, double *verbose)
  {
    long i;
    char *c;
    if (!str || !strlen(str)) return(0);
    c = malloc(strlen(str));
    for (i=0; i<strlen(str); i++) c[i] = tolower(str[i]);
    /* setup options */
    if (strstr(str," k ") || strstr(str," q ") || strstr(str,"wavevector"))
      *Type = WAVEVECTOR_TABLE;
    if (strstr(str,"omega") || strstr(str," e ") || strstr(str,"energy"))
      *Type = ENERGY_TABLE;
    if (strstr(str,"lambda") || strstr(str,"wavelength") || strstr(str," L "))
      *Type = WAVELENGTH_TABLE;
    if (strstr(str,"set")) *Mode  = FLUX_ADAPT_SET;
    if (strstr(str,"add")) *Mode  = FLUX_ADAPT_ADD;
    if (strstr(str,"multiply")) *Mode  = FLUX_ADAPT_MULT;
    if (strstr(str,"verbose")) *verbose = 1;

    return(*Mode);
  }

#endif

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



/* Shared user declarations for all components types 'PowderN'. */
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

/* Declare structures and functions only once in each instrument. */
#ifndef POWDERN_DECL
#define POWDERN_DECL
/* format definitions in the order {j d F2 DW Dd inv2d q F strain} */
#ifndef Crystallographica
#define Crystallographica { 4,5,7,0,0,0,0,0,0 }
#define Fullprof          { 4,0,8,0,0,5,0,0,0 }
#define Lazy              {17,6,0,0,0,0,0,13,0 }
#define Undefined         { 0,0,0,0,0,0,0,0,0 }
#endif

struct line_data
{
double F2;                  /* Value of structure factor */
double q;                   /* Qvector */
int j;                      /* Multiplicity */
double DWfactor;            /* Debye-Waller factor */
double w;                   /* Intrinsic line width */
double Epsilon;             /* Strain=delta_d_d/d shift in ppm */
};

struct line_info_struct
{
struct line_data *list;     /* Reflection array */
int  count;                  /* Number of reflections */
double Dd;
double DWfactor;
double V_0;
double rho;
double at_weight;
double at_nb;
double sigma_a;
double sigma_i;
char   compname[256];
double flag_barns;
int    shape; /* 0 cylinder, 1 box, 2 sphere, 3 OFF file */
int    column_order[9]; /* column signification */
int    flag_warning;
char   type;  /* interaction type of event t=Transmit, i=Incoherent, c=Coherent */
double dq;    /* wavevector transfer [Angs-1] */
double Epsilon; /* global strain in ppm */
double XsectionFactor;
double my_s_v2_sum;
double my_a_v;
double my_inc;
double *w_v,*q_v, *my_s_v2;
double radius_i,xwidth_i,yheight_i,zdepth_i;
double v; /* last velocity (cached) */
      double Nq;
      int    nb_reuses, nb_refl, nb_refl_count;
      double v_min, v_max;
      double xs_Nq[CHAR_BUF_LENGTH];
      double xs_sum[CHAR_BUF_LENGTH];
      double neutron_passed;
      long   xs_compute, xs_reuse, xs_calls;
    };

  off_struct offdata;

  // PN_list_compare *****************************************************************

  int PN_list_compare (void const *a, void const *b)
  {
     struct line_data const *pa = a;
     struct line_data const *pb = b;
     double s = pa->q - pb->q;

     if (!s) return 0;
     else    return (s < 0 ? -1 : 1);
  } /* PN_list_compare */

  int read_line_data(char *SC_file, struct line_info_struct *info)
  {
    struct line_data *list = NULL;
    int    size = 0;
    t_Table sTable; /* sample data table structure from SC_file */
    int    i=0;
    int    mult_count  =0;
    char   flag=0;
    double q_count=0, j_count=0, F2_count=0;
    char **parsing;
    int    list_count=0;

    if (!SC_file || !strlen(SC_file) || !strcmp(SC_file, "NULL")) {
      MPI_MASTER(
      printf("PowderN: %s: Using incoherent elastic scattering only.\n",
          info->compname);
      );
      info->count = 0;
      return(0);
    }
    Table_Read(&sTable, SC_file, 1); /* read 1st block data from SC_file into sTable*/

    /* parsing of header */
    parsing = Table_ParseHeader(sTable.header,
      "Vc","V_0",
      "sigma_abs","sigma_a ",
      "sigma_inc","sigma_i ",
      "column_j",
      "column_d",
      "column_F2",
      "column_DW",
      "column_Dd",
      "column_inv2d", "column_1/2d", "column_sintheta/lambda",
      "column_q", /* 14 */
      "DW", "Debye_Waller",
      "delta_d_d/d",
      "column_F ",
      "V_rho",
      "density",
      "weight",
      "nb_atoms","multiplicity", /* 23 */
      "column_ppm","column_strain",
      NULL);

    if (parsing) {
      if (parsing[0] && !info->V_0)     info->V_0    =atof(parsing[0]);
      if (parsing[1] && !info->V_0)     info->V_0    =atof(parsing[1]);
      if (parsing[2] && !info->sigma_a) info->sigma_a=atof(parsing[2]);
      if (parsing[3] && !info->sigma_a) info->sigma_a=atof(parsing[3]);
      if (parsing[4] && !info->sigma_i) info->sigma_i=atof(parsing[4]);
      if (parsing[5] && !info->sigma_i) info->sigma_i=atof(parsing[5]);
      if (parsing[6])                   info->column_order[0]=atoi(parsing[6]);
      if (parsing[7])                   info->column_order[1]=atoi(parsing[7]);
      if (parsing[8])                   info->column_order[2]=atoi(parsing[8]);
      if (parsing[9])                   info->column_order[3]=atoi(parsing[9]);
      if (parsing[10])                  info->column_order[4]=atoi(parsing[10]);
      if (parsing[11])                  info->column_order[5]=atoi(parsing[11]);
      if (parsing[12])                  info->column_order[5]=atoi(parsing[12]);
      if (parsing[13])                  info->column_order[5]=atoi(parsing[13]);
      if (parsing[14])                  info->column_order[6]=atoi(parsing[14]);
      if (parsing[15] && info->DWfactor<=0)    info->DWfactor=atof(parsing[15]);
      if (parsing[16] && info->DWfactor<=0)    info->DWfactor=atof(parsing[16]);
      if (parsing[17] && info->Dd <0)          info->Dd      =atof(parsing[17]);
      if (parsing[18])                  info->column_order[7]=atoi(parsing[18]);
      if (parsing[19] && !info->V_0)    info->V_0    =1/atof(parsing[19]);
      if (parsing[20] && !info->rho)    info->rho    =atof(parsing[20]);
      if (parsing[21] && !info->at_weight)     info->at_weight    =atof(parsing[21]);
      if (parsing[22] && info->at_nb <= 1)  info->at_nb    =atof(parsing[22]);
      if (parsing[23] && info->at_nb <= 1)  info->at_nb    =atof(parsing[23]);
      if (parsing[24])                  info->column_order[8]=atoi(parsing[24]);
      if (parsing[25])                  info->column_order[8]=atoi(parsing[25]);
      for (i=0; i<=25; i++) if (parsing[i]) free(parsing[i]);
      free(parsing);
    }

    if (!sTable.rows)
      exit(fprintf(stderr, "PowderN: %s: Error: The number of rows in %s "
       "should be at least %d\n", info->compname, SC_file, 1));
    else
      size = sTable.rows;

    MPI_MASTER(
    Table_Info(sTable);
    printf("PowderN: %s: Reading %d rows from %s\n",
          info->compname, size, SC_file);
    );

    if (info->column_order[0] == 4 && info->flag_barns !=0)
    MPI_MASTER(
      printf("PowderN: %s: Powder file probably of type Crystallographica/Fullprof (lau)\n"
           "WARNING: but F2 unit is set to barns=1 (barns). Intensity might be 100 times too high.\n",
           info->compname);
    );
    if (info->column_order[0] == 17 && info->flag_barns == 0)
    MPI_MASTER(
      printf("PowderN: %s: Powder file probably of type Lazy Pulver (laz)\n"
           "WARNING: but F2 unit is set to barns=0 (fm^2). Intensity might be 100 times too low.\n",
           info->compname);
    );
    /* allocate line_data array */
    list = (struct line_data*)malloc(size*sizeof(struct line_data));

    for (i=0; i<size; i++)
    {
      /*      printf("Reading in line %i\n",i);*/
      double j=0, d=0, w=0, q=0, DWfactor=0, F2=0, Epsilon=0;
      int index;

      if (info->Dd >= 0)      w         = info->Dd;
      if (info->DWfactor > 0) DWfactor  = info->DWfactor;
      if (info->Epsilon)      Epsilon   = info->Epsilon*1e-6;

      /* get data from table using columns {j d F2 DW Dd inv2d q F} */
      /* column indexes start at 1, thus need to substract 1 */
      if (info->column_order[0] >0)
        j = Table_Index(sTable, i, info->column_order[0]-1);
      if (info->column_order[1] >0)
        d = Table_Index(sTable, i, info->column_order[1]-1);
      if (info->column_order[2] >0)
        F2 = Table_Index(sTable, i, info->column_order[2]-1);
      if (info->column_order[3] >0)
        DWfactor = Table_Index(sTable, i, info->column_order[3]-1);
      if (info->column_order[4] >0)
        w = Table_Index(sTable, i, info->column_order[4]-1);
      if (info->column_order[5] >0)
        { d = Table_Index(sTable, i, info->column_order[5]-1);
          d = (d > 0? 1/d/2 : 0); }
      if (info->column_order[6] >0)
        { q = Table_Index(sTable, i, info->column_order[6]-1);
          d = (q > 0 ? 2*PI/q : 0); }
      if (info->column_order[7] >0  && !F2)
        { F2 = Table_Index(sTable, i, info->column_order[7]-1); F2 *= F2; }
      if (info->column_order[8] >0  && !Epsilon)
        { Epsilon = Table_Index(sTable, i, info->column_order[8]-1)*1e-6; }

      /* assign and check values */
      j        = (j > 0 ? j : 0);
      q        = (d > 0 ? 2*PI/d : 0); /* this is q */
      if (Epsilon && fabs(Epsilon) < 1e6) {
        q     -= Epsilon*q; /* dq/q = -delta_d_d/d = -Epsilon */
      }
      DWfactor = (DWfactor > 0 ? DWfactor : 1);
      w        = (w>0 ? w : 0); /* this is q and d relative spreading */
      F2       = (F2 >= 0 ? F2 : 0);
      if (j == 0 || q == 0) {
        MPI_MASTER(
        printf("PowderN: %s: line %i has invalid definition\n"
               "         (mult=0 or q=0 or d=0)\n", info->compname, i);
        );
        continue;
      }
      list[list_count].j = j;
      list[list_count].q = q;
      list[list_count].DWfactor = DWfactor;
      list[list_count].w = w;
      list[list_count].F2= F2;
      list[list_count].Epsilon = Epsilon;

      /* adjust multiplicity if j-column + multiple d-spacing lines */
      /* if  d = previous d, increase line duplication index */
      if (!q_count)      q_count  = q;
      if (!j_count)      j_count  = j;
      if (!F2_count)     F2_count = F2;
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
        if ((mult_count && list_count>0)
            && (mult_count == list[list_count-1].j
                || ((list_count < size) && (i == size - 1)
                    && (mult_count == list[list_count].j))) ) {
          MPI_MASTER(
          printf("PowderN: %s: Set multiplicity to 1 for lines [%i:%i]\n"
                  "         (d-spacing %g is duplicated %i times)\n",
            info->compname, list_count-mult_count, list_count-1, list[list_count-1].q, mult_count);
          );
          for (index=list_count-mult_count; index<list_count; list[index++].j = 1);
          mult_count = 1;
          q_count   = q;
          j_count   = j;
          F2_count  = F2;
        }
        if (i == size-1) list_count--;
        flag=0;
      }
      list_count++;
    } /* end for */

    Table_Free(&sTable);

    /* sort the list with increasing q */
    qsort(list, list_count, sizeof(struct line_data),  PN_list_compare);

    MPI_MASTER(
    printf("PowderN: %s: Read %i reflections from file '%s'\n",
      info->compname, list_count, SC_file);
    );

    info->list  = list;
    info->count = list_count;

    return(list_count);
  } /* read_line_data */


/* computes the number of possible reflections (return value), and the total xsection 'sum' */
/* this routine looks for a pre-computed value in the Nq and sum cache tables               */
/* when found, the earch starts from the corresponding lower element in the table           */
int calc_xsect(double v, double *qv, double *my_sv2, int count, double *sum,
  struct line_info_struct *line_info) {
  int    Nq = 0, line=0, line0=0;
  *sum=0;

  /* check if a line_info element has been recorded already */
  if (v >= line_info->v_min && v <= line_info->v_max && line_info->neutron_passed >= CHAR_BUF_LENGTH) {
    line = (int)floor(v - line_info->v_min)*CHAR_BUF_LENGTH/(line_info->v_max - line_info->v_min);
    Nq    = line_info->xs_Nq[line];
    *sum  = line_info->xs_sum[line];
    if (!Nq && *sum == 0) {
      /* not yet set: we compute the sum up to the corresponding speed in the table cache */
      double line_v = line_info->v_min + line*(line_info->v_max - line_info->v_min)/CHAR_BUF_LENGTH;
      for(line0=0; line0<count; line0++) {
        if (qv[line0] <= 2*line_v) { /* q < 2*kf: restrict structural range */
          *sum += my_sv2[line0];
          if (Nq < line0+1) Nq=line0+1; /* determine maximum line index which can scatter */
        } else break;
      }
      line_info->xs_Nq[line] = Nq;
      line_info->xs_sum[line]= *sum;
      line_info->xs_compute++;
    } else line_info->xs_reuse++;
    line0 = Nq - 1;
  }

  line_info->xs_calls++;

  for(line=line0; line<count; line++) {
    if (qv[line] <= 2*v) { /* q < 2*kf: restrict structural range */
      *sum += my_sv2[line];
      if (Nq < line+1) Nq=line+1; /* determine maximum line index which can scatter */
    } else break;
  }

  return(Nq);
} /* calc_xsect */

#endif /* !POWDERN_DECL */



/* ************************************************************************** */
/*             End of SHARE user declarations for all components              */
/* ************************************************************************** */


/* ********************** component definition declarations. **************** */

/* component Progress=Progress_bar() [1] DECLARE */
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
  char     _name[256]; /* e.g. Progress */
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
_class_Progress_bar _Progress_var;
_class_Progress_bar *_Progress = &_Progress_var;
#pragma acc declare create ( _Progress_var )
#pragma acc declare create ( _Progress )

/* component Reactorbeam=Arm() [2] DECLARE */
/* Parameter definition for component type 'Arm' */
struct _struct_Arm_parameters {
  char Arm_has_no_parameters;
}; /* _struct_Arm_parameters */
typedef struct _struct_Arm_parameters _class_Arm_parameters;

/* Parameters for component type 'Arm' */
struct _struct_Arm {
  char     _name[256]; /* e.g. Reactorbeam */
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
_class_Arm _Reactorbeam_var;
_class_Arm *_Reactorbeam = &_Reactorbeam_var;
#pragma acc declare create ( _Reactorbeam_var )
#pragma acc declare create ( _Reactorbeam )

_class_Arm _Prim_axes_var;
_class_Arm *_Prim_axes = &_Prim_axes_var;
#pragma acc declare create ( _Prim_axes_var )
#pragma acc declare create ( _Prim_axes )

/* component Source=Source_gen() [4] DECLARE */
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

/* component PSD_Source=PSD_monitor() [5] DECLARE */
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
  char     _name[256]; /* e.g. PSD_Source */
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
_class_PSD_monitor _PSD_Source_var;
_class_PSD_monitor *_PSD_Source = &_PSD_Source_var;
#pragma acc declare create ( _PSD_Source_var )
#pragma acc declare create ( _PSD_Source )

/* component LAM_Source=L_monitor() [6] DECLARE */
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
  char     _name[256]; /* e.g. LAM_Source */
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
_class_L_monitor _LAM_Source_var;
_class_L_monitor *_LAM_Source = &_LAM_Source_var;
#pragma acc declare create ( _LAM_Source_var )
#pragma acc declare create ( _LAM_Source )

/* component Window_before_filter=Slit() [7] DECLARE */
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
  char     _name[256]; /* e.g. Window_before_filter */
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
_class_Slit _Window_before_filter_var;
_class_Slit *_Window_before_filter = &_Window_before_filter_var;
#pragma acc declare create ( _Window_before_filter_var )
#pragma acc declare create ( _Window_before_filter )

/* component Sapphire_filter=Filter_gen() [8] DECLARE */
/* Parameter definition for component type 'Filter_gen' */
struct _struct_Filter_gen_parameters {
  /* Component type 'Filter_gen' setting parameters */
  char filename[16384];
  char options[16384];
  MCNUM xmin;
  MCNUM xmax;
  MCNUM ymin;
  MCNUM ymax;
  MCNUM xwidth;
  MCNUM yheight;
  MCNUM thickness;
  MCNUM scaling;
  MCNUM verbose;
  /* Component type 'Filter_gen' private parameters */
  /* Component type 'Filter_gen' DECLARE code stored as structure members */
  char Mode_Table;
  char Type_Table;
  t_Table pTable;
}; /* _struct_Filter_gen_parameters */
typedef struct _struct_Filter_gen_parameters _class_Filter_gen_parameters;

/* Parameters for component type 'Filter_gen' */
struct _struct_Filter_gen {
  char     _name[256]; /* e.g. Sapphire_filter */
  char     _type[256]; /* Filter_gen */
  long     _index; /* e.g. 2 index in TRACE list */
  Coords   _position_absolute;
  Coords   _position_relative; /* wrt PREVIOUS */
  Rotation _rotation_absolute;
  Rotation _rotation_relative; /* wrt PREVIOUS */
  int      _rotation_is_identity;
  _class_Filter_gen_parameters _parameters;
};
typedef struct _struct_Filter_gen _class_Filter_gen;
_class_Filter_gen _Sapphire_filter_var;
_class_Filter_gen *_Sapphire_filter = &_Sapphire_filter_var;
#pragma acc declare create ( _Sapphire_filter_var )
#pragma acc declare create ( _Sapphire_filter )

_class_Slit _Window_after_filter_var;
_class_Slit *_Window_after_filter = &_Window_after_filter_var;
#pragma acc declare create ( _Window_after_filter_var )
#pragma acc declare create ( _Window_after_filter )

_class_L_monitor _LAM_After_sapphire_var;
_class_L_monitor *_LAM_After_sapphire = &_LAM_After_sapphire_var;
#pragma acc declare create ( _LAM_After_sapphire_var )
#pragma acc declare create ( _LAM_After_sapphire )

_class_PSD_monitor _PSD_After_sapphire_var;
_class_PSD_monitor *_PSD_After_sapphire = &_PSD_After_sapphire_var;
#pragma acc declare create ( _PSD_After_sapphire_var )
#pragma acc declare create ( _PSD_After_sapphire )

_class_Slit _HighResOutlet_var;
_class_Slit *_HighResOutlet = &_HighResOutlet_var;
#pragma acc declare create ( _HighResOutlet_var )
#pragma acc declare create ( _HighResOutlet )

_class_Slit _HighIntensityOutlet_var;
_class_Slit *_HighIntensityOutlet = &_HighIntensityOutlet_var;
#pragma acc declare create ( _HighIntensityOutlet_var )
#pragma acc declare create ( _HighIntensityOutlet )

_class_PSD_monitor _PSD_After_Outlet_var;
_class_PSD_monitor *_PSD_After_Outlet = &_PSD_After_Outlet_var;
#pragma acc declare create ( _PSD_After_Outlet_var )
#pragma acc declare create ( _PSD_After_Outlet )

_class_Arm _Mono_axis_var;
_class_Arm *_Mono_axis = &_Mono_axis_var;
#pragma acc declare create ( _Mono_axis_var )
#pragma acc declare create ( _Mono_axis )

/* component Blade_1=Monochromator_curved() [16] DECLARE */
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
  char     _name[256]; /* e.g. Blade_1 */
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
_class_Monochromator_curved _Blade_1_var;
_class_Monochromator_curved *_Blade_1 = &_Blade_1_var;
#pragma acc declare create ( _Blade_1_var )
#pragma acc declare create ( _Blade_1 )

_class_Monochromator_curved _Blade_2_var;
_class_Monochromator_curved *_Blade_2 = &_Blade_2_var;
#pragma acc declare create ( _Blade_2_var )
#pragma acc declare create ( _Blade_2 )

_class_Monochromator_curved _Blade_3_var;
_class_Monochromator_curved *_Blade_3 = &_Blade_3_var;
#pragma acc declare create ( _Blade_3_var )
#pragma acc declare create ( _Blade_3 )

_class_Monochromator_curved _Blade_4_var;
_class_Monochromator_curved *_Blade_4 = &_Blade_4_var;
#pragma acc declare create ( _Blade_4_var )
#pragma acc declare create ( _Blade_4 )

_class_Monochromator_curved _Blade_5_var;
_class_Monochromator_curved *_Blade_5 = &_Blade_5_var;
#pragma acc declare create ( _Blade_5_var )
#pragma acc declare create ( _Blade_5 )

_class_Monochromator_curved _Blade_6_var;
_class_Monochromator_curved *_Blade_6 = &_Blade_6_var;
#pragma acc declare create ( _Blade_6_var )
#pragma acc declare create ( _Blade_6 )

_class_Monochromator_curved _Blade_7_var;
_class_Monochromator_curved *_Blade_7 = &_Blade_7_var;
#pragma acc declare create ( _Blade_7_var )
#pragma acc declare create ( _Blade_7 )

_class_Monochromator_curved _Blade_8_var;
_class_Monochromator_curved *_Blade_8 = &_Blade_8_var;
#pragma acc declare create ( _Blade_8_var )
#pragma acc declare create ( _Blade_8 )

_class_Monochromator_curved _Blade_9_var;
_class_Monochromator_curved *_Blade_9 = &_Blade_9_var;
#pragma acc declare create ( _Blade_9_var )
#pragma acc declare create ( _Blade_9 )

_class_Monochromator_curved _Blade_10_var;
_class_Monochromator_curved *_Blade_10 = &_Blade_10_var;
#pragma acc declare create ( _Blade_10_var )
#pragma acc declare create ( _Blade_10 )

_class_Monochromator_curved _Blade_11_var;
_class_Monochromator_curved *_Blade_11 = &_Blade_11_var;
#pragma acc declare create ( _Blade_11_var )
#pragma acc declare create ( _Blade_11 )

_class_Monochromator_curved _Blade_12_var;
_class_Monochromator_curved *_Blade_12 = &_Blade_12_var;
#pragma acc declare create ( _Blade_12_var )
#pragma acc declare create ( _Blade_12 )

_class_Monochromator_curved _Blade_13_var;
_class_Monochromator_curved *_Blade_13 = &_Blade_13_var;
#pragma acc declare create ( _Blade_13_var )
#pragma acc declare create ( _Blade_13 )

_class_PSD_monitor _PSD_At_sec_shutter_var;
_class_PSD_monitor *_PSD_At_sec_shutter = &_PSD_At_sec_shutter_var;
#pragma acc declare create ( _PSD_At_sec_shutter_var )
#pragma acc declare create ( _PSD_At_sec_shutter )

_class_L_monitor _LAM_At_sec_shutter_var;
_class_L_monitor *_LAM_At_sec_shutter = &_LAM_At_sec_shutter_var;
#pragma acc declare create ( _LAM_At_sec_shutter_var )
#pragma acc declare create ( _LAM_At_sec_shutter )

_class_Slit _Inside_chamber_collimator_var;
_class_Slit *_Inside_chamber_collimator = &_Inside_chamber_collimator_var;
#pragma acc declare create ( _Inside_chamber_collimator_var )
#pragma acc declare create ( _Inside_chamber_collimator )

_class_PSD_monitor _PSD_After_inside_chamber_collimator_var;
_class_PSD_monitor *_PSD_After_inside_chamber_collimator = &_PSD_After_inside_chamber_collimator_var;
#pragma acc declare create ( _PSD_After_inside_chamber_collimator_var )
#pragma acc declare create ( _PSD_After_inside_chamber_collimator )

_class_Slit _Outside_chamber_collimator_var;
_class_Slit *_Outside_chamber_collimator = &_Outside_chamber_collimator_var;
#pragma acc declare create ( _Outside_chamber_collimator_var )
#pragma acc declare create ( _Outside_chamber_collimator )

_class_PSD_monitor _PSD_Outside_chamber_collimator_1_var;
_class_PSD_monitor *_PSD_Outside_chamber_collimator_1 = &_PSD_Outside_chamber_collimator_1_var;
#pragma acc declare create ( _PSD_Outside_chamber_collimator_1_var )
#pragma acc declare create ( _PSD_Outside_chamber_collimator_1 )

_class_PSD_monitor _PSD_Outside_chamber_collimator_var;
_class_PSD_monitor *_PSD_Outside_chamber_collimator = &_PSD_Outside_chamber_collimator_var;
#pragma acc declare create ( _PSD_Outside_chamber_collimator_var )
#pragma acc declare create ( _PSD_Outside_chamber_collimator )

_class_Slit _Incident_slit_h_var;
_class_Slit *_Incident_slit_h = &_Incident_slit_h_var;
#pragma acc declare create ( _Incident_slit_h_var )
#pragma acc declare create ( _Incident_slit_h )

_class_Slit _Incident_slit_w_var;
_class_Slit *_Incident_slit_w = &_Incident_slit_w_var;
#pragma acc declare create ( _Incident_slit_w_var )
#pragma acc declare create ( _Incident_slit_w )

_class_PSD_monitor _PSD_After_Incident_slit_w_var;
_class_PSD_monitor *_PSD_After_Incident_slit_w = &_PSD_After_Incident_slit_w_var;
#pragma acc declare create ( _PSD_After_Incident_slit_w_var )
#pragma acc declare create ( _PSD_After_Incident_slit_w )

_class_Arm _Center_of_rotation_var;
_class_Arm *_Center_of_rotation = &_Center_of_rotation_var;
#pragma acc declare create ( _Center_of_rotation_var )
#pragma acc declare create ( _Center_of_rotation )

_class_Arm _Sample_rotation_var;
_class_Arm *_Sample_rotation = &_Sample_rotation_var;
#pragma acc declare create ( _Sample_rotation_var )
#pragma acc declare create ( _Sample_rotation )

_class_Arm _Sample_location_var;
_class_Arm *_Sample_location = &_Sample_location_var;
#pragma acc declare create ( _Sample_location_var )
#pragma acc declare create ( _Sample_location )

_class_PSD_monitor _PSD_Center_of_rotation_var;
_class_PSD_monitor *_PSD_Center_of_rotation = &_PSD_Center_of_rotation_var;
#pragma acc declare create ( _PSD_Center_of_rotation_var )
#pragma acc declare create ( _PSD_Center_of_rotation )

/* component DIV_Center_of_rotation=Divergence_monitor() [43] DECLARE */
/* Parameter definition for component type 'Divergence_monitor' */
struct _struct_Divergence_monitor_parameters {
  /* Component type 'Divergence_monitor' setting parameters */
  MCNUM nh;
  MCNUM nv;
  char filename[16384];
  MCNUM xmin;
  MCNUM xmax;
  MCNUM ymin;
  MCNUM ymax;
  MCNUM xwidth;
  MCNUM yheight;
  MCNUM maxdiv_h;
  MCNUM maxdiv_v;
  MCNUM restore_neutron;
  MCNUM nx;
  MCNUM ny;
  MCNUM nz;
  /* Component type 'Divergence_monitor' private parameters */
  /* Component type 'Divergence_monitor' DECLARE code stored as structure members */
  DArray2d Div_N;
  DArray2d Div_p;
  DArray2d Div_p2;
}; /* _struct_Divergence_monitor_parameters */
typedef struct _struct_Divergence_monitor_parameters _class_Divergence_monitor_parameters;

/* Parameters for component type 'Divergence_monitor' */
struct _struct_Divergence_monitor {
  char     _name[256]; /* e.g. DIV_Center_of_rotation */
  char     _type[256]; /* Divergence_monitor */
  long     _index; /* e.g. 2 index in TRACE list */
  Coords   _position_absolute;
  Coords   _position_relative; /* wrt PREVIOUS */
  Rotation _rotation_absolute;
  Rotation _rotation_relative; /* wrt PREVIOUS */
  int      _rotation_is_identity;
  _class_Divergence_monitor_parameters _parameters;
};
typedef struct _struct_Divergence_monitor _class_Divergence_monitor;
_class_Divergence_monitor _DIV_Center_of_rotation_var;
_class_Divergence_monitor *_DIV_Center_of_rotation = &_DIV_Center_of_rotation_var;
#pragma acc declare create ( _DIV_Center_of_rotation_var )
#pragma acc declare create ( _DIV_Center_of_rotation )

/* component Sample=PowderN() [44] DECLARE */
/* Parameter definition for component type 'PowderN' */
struct _struct_PowderN_parameters {
  /* Component type 'PowderN' setting parameters */
  char reflections[16384];
  char geometry[16384];
  MCNUM* format;
  MCNUM radius;
  MCNUM yheight;
  MCNUM xwidth;
  MCNUM zdepth;
  MCNUM thickness;
  MCNUM pack;
  MCNUM Vc;
  MCNUM sigma_abs;
  MCNUM sigma_inc;
  MCNUM delta_d_d;
  MCNUM p_inc;
  MCNUM p_transmit;
  MCNUM DW;
  MCNUM nb_atoms;
  MCNUM d_phi;
  MCNUM p_interact;
  MCNUM concentric;
  MCNUM density;
  MCNUM weight;
  MCNUM barns;
  MCNUM Strain;
  MCNUM focus_flip;
  /* Component type 'PowderN' private parameters */
  /* Component type 'PowderN' DECLARE code stored as structure members */
  struct line_info_struct line_info;
  double *columns;
  off_struct offdata;
}; /* _struct_PowderN_parameters */
typedef struct _struct_PowderN_parameters _class_PowderN_parameters;

/* Parameters for component type 'PowderN' */
struct _struct_PowderN {
  char     _name[256]; /* e.g. Sample */
  char     _type[256]; /* PowderN */
  long     _index; /* e.g. 2 index in TRACE list */
  Coords   _position_absolute;
  Coords   _position_relative; /* wrt PREVIOUS */
  Rotation _rotation_absolute;
  Rotation _rotation_relative; /* wrt PREVIOUS */
  int      _rotation_is_identity;
  _class_PowderN_parameters _parameters;
};
typedef struct _struct_PowderN _class_PowderN;
_class_PowderN _Sample_var;
_class_PowderN *_Sample = &_Sample_var;
#pragma acc declare create ( _Sample_var )
#pragma acc declare create ( _Sample )

_class_Arm _Det_axis_var;
_class_Arm *_Det_axis = &_Det_axis_var;
#pragma acc declare create ( _Det_axis_var )
#pragma acc declare create ( _Det_axis )

_class_Arm _Det_axis_2_var;
_class_Arm *_Det_axis_2 = &_Det_axis_2_var;
#pragma acc declare create ( _Det_axis_2_var )
#pragma acc declare create ( _Det_axis_2 )

_class_Arm _Det_axis_3_var;
_class_Arm *_Det_axis_3 = &_Det_axis_3_var;
#pragma acc declare create ( _Det_axis_3_var )
#pragma acc declare create ( _Det_axis_3 )

_class_Arm _Det_axis_4_var;
_class_Arm *_Det_axis_4 = &_Det_axis_4_var;
#pragma acc declare create ( _Det_axis_4_var )
#pragma acc declare create ( _Det_axis_4 )

/* component RadColl=Collimator_radial() [49] DECLARE */
/* Parameter definition for component type 'Collimator_radial' */
struct _struct_Collimator_radial_parameters {
  /* Component type 'Collimator_radial' setting parameters */
  MCNUM xwidth;
  MCNUM yheight;
  MCNUM length;
  MCNUM divergence;
  MCNUM transmission;
  MCNUM theta_min;
  MCNUM theta_max;
  MCNUM nchan;
  MCNUM radius;
  MCNUM nslit;
  MCNUM roc;
  MCNUM verbose;
  MCNUM approx;
  /* Component type 'Collimator_radial' private parameters */
  /* Component type 'Collimator_radial' DECLARE code stored as structure members */
double width_of_slit;
double width_of_Soller;
double slit_theta;
}; /* _struct_Collimator_radial_parameters */
typedef struct _struct_Collimator_radial_parameters _class_Collimator_radial_parameters;

/* Parameters for component type 'Collimator_radial' */
struct _struct_Collimator_radial {
  char     _name[256]; /* e.g. RadColl */
  char     _type[256]; /* Collimator_radial */
  long     _index; /* e.g. 2 index in TRACE list */
  Coords   _position_absolute;
  Coords   _position_relative; /* wrt PREVIOUS */
  Rotation _rotation_absolute;
  Rotation _rotation_relative; /* wrt PREVIOUS */
  int      _rotation_is_identity;
  _class_Collimator_radial_parameters _parameters;
};
typedef struct _struct_Collimator_radial _class_Collimator_radial;
_class_Collimator_radial _RadColl_var;
_class_Collimator_radial *_RadColl = &_RadColl_var;
#pragma acc declare create ( _RadColl_var )
#pragma acc declare create ( _RadColl )

_class_PSD_monitor _PSD_Detector_var;
_class_PSD_monitor *_PSD_Detector = &_PSD_Detector_var;
#pragma acc declare create ( _PSD_Detector_var )
#pragma acc declare create ( _PSD_Detector )

_class_PSD_monitor _PSD_Detector_2_var;
_class_PSD_monitor *_PSD_Detector_2 = &_PSD_Detector_2_var;
#pragma acc declare create ( _PSD_Detector_2_var )
#pragma acc declare create ( _PSD_Detector_2 )

_class_PSD_monitor _PSD_Detector_3_var;
_class_PSD_monitor *_PSD_Detector_3 = &_PSD_Detector_3_var;
#pragma acc declare create ( _PSD_Detector_3_var )
#pragma acc declare create ( _PSD_Detector_3 )

_class_PSD_monitor _PSD_Detector_4_var;
_class_PSD_monitor *_PSD_Detector_4 = &_PSD_Detector_4_var;
#pragma acc declare create ( _PSD_Detector_4_var )
#pragma acc declare create ( _PSD_Detector_4 )

int mcNUMCOMP = 53;

/* User declarations from instrument definition. Can define functions. */
  double hi_res, port_takeoff;

  double mono_Si_type, mono_mosh, mono_mosv;
  double mono_dx, mono_dy, mono_dz;
  double mono_takeoff, mono_dtilt, mono_r_h, mono_r_v;
  double mono_r_req_v;
  double mono_pts_v, mono_pts_h, focal_dist;
  double a,b,c,d,y;      //Gaussian fit params

  double inc_slit_rot, inc_slit_dx, inc_slit_to_cor;
  double inc_slit_width, inc_slit_height, inc_slit_sep;
  double mono_to_cor;
  double sample_dx, sample_dy, sample_dz, sample_dom;

  double det_takeoff, cor_to_det, dangle_interest, ndet, det_cover_angle;
  double det_width, det_height;

  double diff_slit_dx, diff_slit_to_cor;
  double diff_slit_width, diff_slit_height;

  double inc_slit_xmin, inc_slit_xmax, inc_slit_ymin, inc_slit_ymax;
  double inc_slit_xmin_h, inc_slit_xmax_h, inc_slit_ymin_w, inc_slit_ymax_w;
  double diff_slit_xmin, diff_slit_xmax, diff_slit_ymin, diff_slit_ymax;
  double wafer_d, start_wafer_pos, mono_turns, mono_Rh_req, mono_d, mono_q;
  double lam;
  double as;
  int msw;

  double chamber_col_start, chamber_col_length, outside_chamber_collimator_w, outside_chamber_collimator_h;

  double from_col=1;
  double full_instrument;

#undef compcurname
#undef compcurtype
#undef compcurindex
/* end of instrument 'SAFARI_PITSI' and components DECLARE */

/* *****************************************************************************
* instrument 'SAFARI_PITSI' and components INITIALISE
***************************************************************************** */

/* component Progress=Progress_bar() SETTING, POSITION/ROTATION */
int _Progress_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Progress_setpos] component Progress=Progress_bar() SETTING [/usr/share/mcstas/3.0-dev/misc/Progress_bar.comp:57]");
  stracpy(_Progress->_name, "Progress", 16384);
  stracpy(_Progress->_type, "Progress_bar", 16384);
  _Progress->_index=1;
  if("NULL" && strlen("NULL"))
    stracpy(_Progress->_parameters.profile, "NULL" ? "NULL" : "", 16384);
  else 
  _Progress->_parameters.profile[0]='\0';
  #define profile (_Progress->_parameters.profile)
  _Progress->_parameters.percent = 1;
  #define percent (_Progress->_parameters.percent)
  _Progress->_parameters.flag_save = 0;
  #define flag_save (_Progress->_parameters.flag_save)
  _Progress->_parameters.minutes = 0;
  #define minutes (_Progress->_parameters.minutes)

  #define IntermediateCnts (_Progress->_parameters.IntermediateCnts)
  #define StartTime (_Progress->_parameters.StartTime)
  #define EndTime (_Progress->_parameters.EndTime)
  #define CurrentTime (_Progress->_parameters.CurrentTime)

  #undef profile
  #undef percent
  #undef flag_save
  #undef minutes
  #undef IntermediateCnts
  #undef StartTime
  #undef EndTime
  #undef CurrentTime
  /* component Progress=Progress_bar() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(_Progress->_rotation_absolute,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_copy(_Progress->_rotation_relative, _Progress->_rotation_absolute);
    _Progress->_rotation_is_identity =  rot_test_identity(_Progress->_rotation_relative);
    _Progress->_position_absolute = coords_set(
      0, 0, 0);
    tc1 = coords_neg(_Progress->_position_absolute);
    _Progress->_position_relative = rot_apply(_Progress->_rotation_absolute, tc1);
  } /* Progress=Progress_bar() AT ROTATED */
  DEBUG_COMPONENT("Progress", _Progress->_position_absolute, _Progress->_rotation_absolute);
  instrument->_position_absolute[1] = _Progress->_position_absolute;
  instrument->_position_relative[1] = _Progress->_position_relative;
  instrument->counter_N[1]  = instrument->counter_P[1] = instrument->counter_P2[1] = 0;
  instrument->counter_AbsorbProp[1]= 0;
  return(0);
} /* _Progress_setpos */

/* component Reactorbeam=Arm() SETTING, POSITION/ROTATION */
int _Reactorbeam_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Reactorbeam_setpos] component Reactorbeam=Arm() SETTING [Arm:0]");
  stracpy(_Reactorbeam->_name, "Reactorbeam", 16384);
  stracpy(_Reactorbeam->_type, "Arm", 16384);
  _Reactorbeam->_index=2;
  /* component Reactorbeam=Arm() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(_Reactorbeam->_rotation_absolute,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_transpose(_Progress->_rotation_absolute, tr1);
    rot_mul(_Reactorbeam->_rotation_absolute, tr1, _Reactorbeam->_rotation_relative);
    _Reactorbeam->_rotation_is_identity =  rot_test_identity(_Reactorbeam->_rotation_relative);
    _Reactorbeam->_position_absolute = coords_set(
      0, 0, 0);
    tc1 = coords_sub(_Progress->_position_absolute, _Reactorbeam->_position_absolute);
    _Reactorbeam->_position_relative = rot_apply(_Reactorbeam->_rotation_absolute, tc1);
  } /* Reactorbeam=Arm() AT ROTATED */
  DEBUG_COMPONENT("Reactorbeam", _Reactorbeam->_position_absolute, _Reactorbeam->_rotation_absolute);
  instrument->_position_absolute[2] = _Reactorbeam->_position_absolute;
  instrument->_position_relative[2] = _Reactorbeam->_position_relative;
  instrument->counter_N[2]  = instrument->counter_P[2] = instrument->counter_P2[2] = 0;
  instrument->counter_AbsorbProp[2]= 0;
  return(0);
} /* _Reactorbeam_setpos */

/* component Prim_axes=Arm() SETTING, POSITION/ROTATION */
int _Prim_axes_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Prim_axes_setpos] component Prim_axes=Arm() SETTING [Arm:0]");
  stracpy(_Prim_axes->_name, "Prim_axes", 16384);
  stracpy(_Prim_axes->_type, "Arm", 16384);
  _Prim_axes->_index=3;
  /* component Prim_axes=Arm() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (instrument->_parameters._port_takeoff)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _Reactorbeam->_rotation_absolute, _Prim_axes->_rotation_absolute);
    rot_transpose(_Reactorbeam->_rotation_absolute, tr1);
    rot_mul(_Prim_axes->_rotation_absolute, tr1, _Prim_axes->_rotation_relative);
    _Prim_axes->_rotation_is_identity =  rot_test_identity(_Prim_axes->_rotation_relative);
    tc1 = coords_set(
      0, 0, 5.140);
    rot_transpose(_Reactorbeam->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Prim_axes->_position_absolute = coords_add(_Reactorbeam->_position_absolute, tc2);
    tc1 = coords_sub(_Reactorbeam->_position_absolute, _Prim_axes->_position_absolute);
    _Prim_axes->_position_relative = rot_apply(_Prim_axes->_rotation_absolute, tc1);
  } /* Prim_axes=Arm() AT ROTATED */
  DEBUG_COMPONENT("Prim_axes", _Prim_axes->_position_absolute, _Prim_axes->_rotation_absolute);
  instrument->_position_absolute[3] = _Prim_axes->_position_absolute;
  instrument->_position_relative[3] = _Prim_axes->_position_relative;
  instrument->counter_N[3]  = instrument->counter_P[3] = instrument->counter_P2[3] = 0;
  instrument->counter_AbsorbProp[3]= 0;
  return(0);
} /* _Prim_axes_setpos */

/* component Source=Source_gen() SETTING, POSITION/ROTATION */
int _Source_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Source_setpos] component Source=Source_gen() SETTING [/usr/share/mcstas/3.0-dev/sources/Source_gen.comp:206]");
  stracpy(_Source->_name, "Source", 16384);
  stracpy(_Source->_type, "Source_gen", 16384);
  _Source->_index=4;
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
  _Source->_parameters.radius = 0.0905;
  #define radius (_Source->_parameters.radius)
  _Source->_parameters.dist = 2.86805;
  #define dist (_Source->_parameters.dist)
  _Source->_parameters.focus_xw = 0.1;
  #define focus_xw (_Source->_parameters.focus_xw)
  _Source->_parameters.focus_yh = 0.05;
  #define focus_yh (_Source->_parameters.focus_yh)
  _Source->_parameters.focus_aw = 0;
  #define focus_aw (_Source->_parameters.focus_aw)
  _Source->_parameters.focus_ah = 0;
  #define focus_ah (_Source->_parameters.focus_ah)
  _Source->_parameters.E0 = 0;
  #define E0 (_Source->_parameters.E0)
  _Source->_parameters.dE = 0;
  #define dE (_Source->_parameters.dE)
  _Source->_parameters.lambda0 = 0;
  #define lambda0 (_Source->_parameters.lambda0)
  _Source->_parameters.dlambda = 0;
  #define dlambda (_Source->_parameters.dlambda)
  _Source->_parameters.I1 = 0;
  #define I1 (_Source->_parameters.I1)
  _Source->_parameters.yheight = 0.1;
  #define yheight (_Source->_parameters.yheight)
  _Source->_parameters.xwidth = 0.1;
  #define xwidth (_Source->_parameters.xwidth)
  _Source->_parameters.verbose = 0;
  #define verbose (_Source->_parameters.verbose)
  _Source->_parameters.T1 = 0;
  #define T1 (_Source->_parameters.T1)
  _Source->_parameters.flux_file_perAA = 0;
  #define flux_file_perAA (_Source->_parameters.flux_file_perAA)
  _Source->_parameters.flux_file_log = 0;
  #define flux_file_log (_Source->_parameters.flux_file_log)
  _Source->_parameters.Lmin = instrument->_parameters._source_lam_min;
  #define Lmin (_Source->_parameters.Lmin)
  _Source->_parameters.Lmax = instrument->_parameters._source_lam_max;
  #define Lmax (_Source->_parameters.Lmax)
  _Source->_parameters.Emin = 0;
  #define Emin (_Source->_parameters.Emin)
  _Source->_parameters.Emax = 0;
  #define Emax (_Source->_parameters.Emax)
  _Source->_parameters.T2 = 0;
  #define T2 (_Source->_parameters.T2)
  _Source->_parameters.I2 = 0;
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
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Reactorbeam->_rotation_absolute, _Source->_rotation_absolute);
    rot_transpose(_Prim_axes->_rotation_absolute, tr1);
    rot_mul(_Source->_rotation_absolute, tr1, _Source->_rotation_relative);
    _Source->_rotation_is_identity =  rot_test_identity(_Source->_rotation_relative);
    tc1 = coords_set(
      0, 0, 0);
    rot_transpose(_Reactorbeam->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Source->_position_absolute = coords_add(_Reactorbeam->_position_absolute, tc2);
    tc1 = coords_sub(_Prim_axes->_position_absolute, _Source->_position_absolute);
    _Source->_position_relative = rot_apply(_Source->_rotation_absolute, tc1);
  } /* Source=Source_gen() AT ROTATED */
  DEBUG_COMPONENT("Source", _Source->_position_absolute, _Source->_rotation_absolute);
  instrument->_position_absolute[4] = _Source->_position_absolute;
  instrument->_position_relative[4] = _Source->_position_relative;
  instrument->counter_N[4]  = instrument->counter_P[4] = instrument->counter_P2[4] = 0;
  instrument->counter_AbsorbProp[4]= 0;
  return(0);
} /* _Source_setpos */

/* component PSD_Source=PSD_monitor() SETTING, POSITION/ROTATION */
int _PSD_Source_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_PSD_Source_setpos] component PSD_Source=PSD_monitor() SETTING [/usr/share/mcstas/3.0-dev/monitors/PSD_monitor.comp:64]");
  stracpy(_PSD_Source->_name, "PSD_Source", 16384);
  stracpy(_PSD_Source->_type, "PSD_monitor", 16384);
  _PSD_Source->_index=5;
  _PSD_Source->_parameters.nx = 100;
  #define nx (_PSD_Source->_parameters.nx)
  _PSD_Source->_parameters.ny = 100;
  #define ny (_PSD_Source->_parameters.ny)
  if("PSD_Source" && strlen("PSD_Source"))
    stracpy(_PSD_Source->_parameters.filename, "PSD_Source" ? "PSD_Source" : "", 16384);
  else 
  _PSD_Source->_parameters.filename[0]='\0';
  #define filename (_PSD_Source->_parameters.filename)
  _PSD_Source->_parameters.xmin = -0.1;
  #define xmin (_PSD_Source->_parameters.xmin)
  _PSD_Source->_parameters.xmax = 0.1;
  #define xmax (_PSD_Source->_parameters.xmax)
  _PSD_Source->_parameters.ymin = -0.1;
  #define ymin (_PSD_Source->_parameters.ymin)
  _PSD_Source->_parameters.ymax = 0.1;
  #define ymax (_PSD_Source->_parameters.ymax)
  _PSD_Source->_parameters.xwidth = 0;
  #define xwidth (_PSD_Source->_parameters.xwidth)
  _PSD_Source->_parameters.yheight = 0;
  #define yheight (_PSD_Source->_parameters.yheight)
  _PSD_Source->_parameters.restore_neutron = 0;
  #define restore_neutron (_PSD_Source->_parameters.restore_neutron)

  #define PSD_N (_PSD_Source->_parameters.PSD_N)
  #define PSD_p (_PSD_Source->_parameters.PSD_p)
  #define PSD_p2 (_PSD_Source->_parameters.PSD_p2)

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
  /* component PSD_Source=PSD_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Reactorbeam->_rotation_absolute, _PSD_Source->_rotation_absolute);
    rot_transpose(_Source->_rotation_absolute, tr1);
    rot_mul(_PSD_Source->_rotation_absolute, tr1, _PSD_Source->_rotation_relative);
    _PSD_Source->_rotation_is_identity =  rot_test_identity(_PSD_Source->_rotation_relative);
    tc1 = coords_set(
      0, 0, 0);
    rot_transpose(_Reactorbeam->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _PSD_Source->_position_absolute = coords_add(_Reactorbeam->_position_absolute, tc2);
    tc1 = coords_sub(_Source->_position_absolute, _PSD_Source->_position_absolute);
    _PSD_Source->_position_relative = rot_apply(_PSD_Source->_rotation_absolute, tc1);
  } /* PSD_Source=PSD_monitor() AT ROTATED */
  DEBUG_COMPONENT("PSD_Source", _PSD_Source->_position_absolute, _PSD_Source->_rotation_absolute);
  instrument->_position_absolute[5] = _PSD_Source->_position_absolute;
  instrument->_position_relative[5] = _PSD_Source->_position_relative;
  instrument->counter_N[5]  = instrument->counter_P[5] = instrument->counter_P2[5] = 0;
  instrument->counter_AbsorbProp[5]= 0;
  return(0);
} /* _PSD_Source_setpos */

/* component LAM_Source=L_monitor() SETTING, POSITION/ROTATION */
int _LAM_Source_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_LAM_Source_setpos] component LAM_Source=L_monitor() SETTING [/usr/share/mcstas/3.0-dev/monitors/L_monitor.comp:65]");
  stracpy(_LAM_Source->_name, "LAM_Source", 16384);
  stracpy(_LAM_Source->_type, "L_monitor", 16384);
  _LAM_Source->_index=6;
  _LAM_Source->_parameters.nL = 100;
  #define nL (_LAM_Source->_parameters.nL)
  if("LAM_Source.out" && strlen("LAM_Source.out"))
    stracpy(_LAM_Source->_parameters.filename, "LAM_Source.out" ? "LAM_Source.out" : "", 16384);
  else 
  _LAM_Source->_parameters.filename[0]='\0';
  #define filename (_LAM_Source->_parameters.filename)
  _LAM_Source->_parameters.xmin = -0.1;
  #define xmin (_LAM_Source->_parameters.xmin)
  _LAM_Source->_parameters.xmax = 0.1;
  #define xmax (_LAM_Source->_parameters.xmax)
  _LAM_Source->_parameters.ymin = -0.1;
  #define ymin (_LAM_Source->_parameters.ymin)
  _LAM_Source->_parameters.ymax = 0.1;
  #define ymax (_LAM_Source->_parameters.ymax)
  _LAM_Source->_parameters.xwidth = 0;
  #define xwidth (_LAM_Source->_parameters.xwidth)
  _LAM_Source->_parameters.yheight = 0;
  #define yheight (_LAM_Source->_parameters.yheight)
  _LAM_Source->_parameters.Lmin = instrument->_parameters._source_lam_min;
  #define Lmin (_LAM_Source->_parameters.Lmin)
  _LAM_Source->_parameters.Lmax = instrument->_parameters._source_lam_max;
  #define Lmax (_LAM_Source->_parameters.Lmax)
  _LAM_Source->_parameters.restore_neutron = 0;
  #define restore_neutron (_LAM_Source->_parameters.restore_neutron)

  #define L_N (_LAM_Source->_parameters.L_N)
  #define L_p (_LAM_Source->_parameters.L_p)
  #define L_p2 (_LAM_Source->_parameters.L_p2)

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
  /* component LAM_Source=L_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Reactorbeam->_rotation_absolute, _LAM_Source->_rotation_absolute);
    rot_transpose(_PSD_Source->_rotation_absolute, tr1);
    rot_mul(_LAM_Source->_rotation_absolute, tr1, _LAM_Source->_rotation_relative);
    _LAM_Source->_rotation_is_identity =  rot_test_identity(_LAM_Source->_rotation_relative);
    tc1 = coords_set(
      0, 0, 0);
    rot_transpose(_Reactorbeam->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _LAM_Source->_position_absolute = coords_add(_Reactorbeam->_position_absolute, tc2);
    tc1 = coords_sub(_PSD_Source->_position_absolute, _LAM_Source->_position_absolute);
    _LAM_Source->_position_relative = rot_apply(_LAM_Source->_rotation_absolute, tc1);
  } /* LAM_Source=L_monitor() AT ROTATED */
  DEBUG_COMPONENT("LAM_Source", _LAM_Source->_position_absolute, _LAM_Source->_rotation_absolute);
  instrument->_position_absolute[6] = _LAM_Source->_position_absolute;
  instrument->_position_relative[6] = _LAM_Source->_position_relative;
  instrument->counter_N[6]  = instrument->counter_P[6] = instrument->counter_P2[6] = 0;
  instrument->counter_AbsorbProp[6]= 0;
  return(0);
} /* _LAM_Source_setpos */

/* component Window_before_filter=Slit() SETTING, POSITION/ROTATION */
int _Window_before_filter_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Window_before_filter_setpos] component Window_before_filter=Slit() SETTING [/usr/share/mcstas/3.0-dev/optics/Slit.comp:47]");
  stracpy(_Window_before_filter->_name, "Window_before_filter", 16384);
  stracpy(_Window_before_filter->_type, "Slit", 16384);
  _Window_before_filter->_index=7;
  _Window_before_filter->_parameters.xmin = -0.1;
  #define xmin (_Window_before_filter->_parameters.xmin)
  _Window_before_filter->_parameters.xmax = 0.1;
  #define xmax (_Window_before_filter->_parameters.xmax)
  _Window_before_filter->_parameters.ymin = -0.1;
  #define ymin (_Window_before_filter->_parameters.ymin)
  _Window_before_filter->_parameters.ymax = 0.1;
  #define ymax (_Window_before_filter->_parameters.ymax)
  _Window_before_filter->_parameters.radius = 0;
  #define radius (_Window_before_filter->_parameters.radius)
  _Window_before_filter->_parameters.xwidth = 0;
  #define xwidth (_Window_before_filter->_parameters.xwidth)
  _Window_before_filter->_parameters.yheight = 0;
  #define yheight (_Window_before_filter->_parameters.yheight)

  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef radius
  #undef xwidth
  #undef yheight
  /* component Window_before_filter=Slit() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Reactorbeam->_rotation_absolute, _Window_before_filter->_rotation_absolute);
    rot_transpose(_LAM_Source->_rotation_absolute, tr1);
    rot_mul(_Window_before_filter->_rotation_absolute, tr1, _Window_before_filter->_rotation_relative);
    _Window_before_filter->_rotation_is_identity =  rot_test_identity(_Window_before_filter->_rotation_relative);
    tc1 = coords_set(
      0, 0, 3.24045);
    rot_transpose(_Reactorbeam->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Window_before_filter->_position_absolute = coords_add(_Reactorbeam->_position_absolute, tc2);
    tc1 = coords_sub(_LAM_Source->_position_absolute, _Window_before_filter->_position_absolute);
    _Window_before_filter->_position_relative = rot_apply(_Window_before_filter->_rotation_absolute, tc1);
  } /* Window_before_filter=Slit() AT ROTATED */
  DEBUG_COMPONENT("Window_before_filter", _Window_before_filter->_position_absolute, _Window_before_filter->_rotation_absolute);
  instrument->_position_absolute[7] = _Window_before_filter->_position_absolute;
  instrument->_position_relative[7] = _Window_before_filter->_position_relative;
  instrument->counter_N[7]  = instrument->counter_P[7] = instrument->counter_P2[7] = 0;
  instrument->counter_AbsorbProp[7]= 0;
  return(0);
} /* _Window_before_filter_setpos */

/* component Sapphire_filter=Filter_gen() SETTING, POSITION/ROTATION */
int _Sapphire_filter_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Sapphire_filter_setpos] component Sapphire_filter=Filter_gen() SETTING [/usr/share/mcstas/3.0-dev/optics/Filter_gen.comp:126]");
  stracpy(_Sapphire_filter->_name, "Sapphire_filter", 16384);
  stracpy(_Sapphire_filter->_type, "Filter_gen", 16384);
  _Sapphire_filter->_index=8;
  if("Al2O3_sapphire.trm" && strlen("Al2O3_sapphire.trm"))
    stracpy(_Sapphire_filter->_parameters.filename, "Al2O3_sapphire.trm" ? "Al2O3_sapphire.trm" : "", 16384);
  else 
  _Sapphire_filter->_parameters.filename[0]='\0';
  #define filename (_Sapphire_filter->_parameters.filename)
  if("multiply" && strlen("multiply"))
    stracpy(_Sapphire_filter->_parameters.options, "multiply" ? "multiply" : "", 16384);
  else 
  _Sapphire_filter->_parameters.options[0]='\0';
  #define options (_Sapphire_filter->_parameters.options)
  _Sapphire_filter->_parameters.xmin = -0.053975;
  #define xmin (_Sapphire_filter->_parameters.xmin)
  _Sapphire_filter->_parameters.xmax = 0.053975;
  #define xmax (_Sapphire_filter->_parameters.xmax)
  _Sapphire_filter->_parameters.ymin = -0.0381;
  #define ymin (_Sapphire_filter->_parameters.ymin)
  _Sapphire_filter->_parameters.ymax = 0.0381;
  #define ymax (_Sapphire_filter->_parameters.ymax)
  _Sapphire_filter->_parameters.xwidth = 0;
  #define xwidth (_Sapphire_filter->_parameters.xwidth)
  _Sapphire_filter->_parameters.yheight = 0;
  #define yheight (_Sapphire_filter->_parameters.yheight)
  _Sapphire_filter->_parameters.thickness = 3.125;
  #define thickness (_Sapphire_filter->_parameters.thickness)
  _Sapphire_filter->_parameters.scaling = 1;
  #define scaling (_Sapphire_filter->_parameters.scaling)
  _Sapphire_filter->_parameters.verbose = 0;
  #define verbose (_Sapphire_filter->_parameters.verbose)

  #define pTable (_Sapphire_filter->_parameters.pTable)
  #define Mode_Table (_Sapphire_filter->_parameters.Mode_Table)
  #define Type_Table (_Sapphire_filter->_parameters.Type_Table)

  #undef filename
  #undef options
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef xwidth
  #undef yheight
  #undef thickness
  #undef scaling
  #undef verbose
  #undef pTable
  #undef Mode_Table
  #undef Type_Table
  /* component Sapphire_filter=Filter_gen() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Window_before_filter->_rotation_absolute, _Sapphire_filter->_rotation_absolute);
    rot_transpose(_Window_before_filter->_rotation_absolute, tr1);
    rot_mul(_Sapphire_filter->_rotation_absolute, tr1, _Sapphire_filter->_rotation_relative);
    _Sapphire_filter->_rotation_is_identity =  rot_test_identity(_Sapphire_filter->_rotation_relative);
    tc1 = coords_set(
      0, 0, 0);
    rot_transpose(_Window_before_filter->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Sapphire_filter->_position_absolute = coords_add(_Window_before_filter->_position_absolute, tc2);
    tc1 = coords_sub(_Window_before_filter->_position_absolute, _Sapphire_filter->_position_absolute);
    _Sapphire_filter->_position_relative = rot_apply(_Sapphire_filter->_rotation_absolute, tc1);
  } /* Sapphire_filter=Filter_gen() AT ROTATED */
  DEBUG_COMPONENT("Sapphire_filter", _Sapphire_filter->_position_absolute, _Sapphire_filter->_rotation_absolute);
  instrument->_position_absolute[8] = _Sapphire_filter->_position_absolute;
  instrument->_position_relative[8] = _Sapphire_filter->_position_relative;
  instrument->counter_N[8]  = instrument->counter_P[8] = instrument->counter_P2[8] = 0;
  instrument->counter_AbsorbProp[8]= 0;
  return(0);
} /* _Sapphire_filter_setpos */

/* component Window_after_filter=Slit() SETTING, POSITION/ROTATION */
int _Window_after_filter_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Window_after_filter_setpos] component Window_after_filter=Slit() SETTING [/usr/share/mcstas/3.0-dev/optics/Slit.comp:47]");
  stracpy(_Window_after_filter->_name, "Window_after_filter", 16384);
  stracpy(_Window_after_filter->_type, "Slit", 16384);
  _Window_after_filter->_index=9;
  _Window_after_filter->_parameters.xmin = -0.05;
  #define xmin (_Window_after_filter->_parameters.xmin)
  _Window_after_filter->_parameters.xmax = 0.05;
  #define xmax (_Window_after_filter->_parameters.xmax)
  _Window_after_filter->_parameters.ymin = -0.032005;
  #define ymin (_Window_after_filter->_parameters.ymin)
  _Window_after_filter->_parameters.ymax = 0.032005;
  #define ymax (_Window_after_filter->_parameters.ymax)
  _Window_after_filter->_parameters.radius = 0;
  #define radius (_Window_after_filter->_parameters.radius)
  _Window_after_filter->_parameters.xwidth = 0;
  #define xwidth (_Window_after_filter->_parameters.xwidth)
  _Window_after_filter->_parameters.yheight = 0;
  #define yheight (_Window_after_filter->_parameters.yheight)

  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef radius
  #undef xwidth
  #undef yheight
  /* component Window_after_filter=Slit() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Sapphire_filter->_rotation_absolute, _Window_after_filter->_rotation_absolute);
    rot_transpose(_Sapphire_filter->_rotation_absolute, tr1);
    rot_mul(_Window_after_filter->_rotation_absolute, tr1, _Window_after_filter->_rotation_relative);
    _Window_after_filter->_rotation_is_identity =  rot_test_identity(_Window_after_filter->_rotation_relative);
    tc1 = coords_set(
      0, 0, 0.15876);
    rot_transpose(_Sapphire_filter->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Window_after_filter->_position_absolute = coords_add(_Sapphire_filter->_position_absolute, tc2);
    tc1 = coords_sub(_Sapphire_filter->_position_absolute, _Window_after_filter->_position_absolute);
    _Window_after_filter->_position_relative = rot_apply(_Window_after_filter->_rotation_absolute, tc1);
  } /* Window_after_filter=Slit() AT ROTATED */
  DEBUG_COMPONENT("Window_after_filter", _Window_after_filter->_position_absolute, _Window_after_filter->_rotation_absolute);
  instrument->_position_absolute[9] = _Window_after_filter->_position_absolute;
  instrument->_position_relative[9] = _Window_after_filter->_position_relative;
  instrument->counter_N[9]  = instrument->counter_P[9] = instrument->counter_P2[9] = 0;
  instrument->counter_AbsorbProp[9]= 0;
  return(0);
} /* _Window_after_filter_setpos */

/* component LAM_After_sapphire=L_monitor() SETTING, POSITION/ROTATION */
int _LAM_After_sapphire_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_LAM_After_sapphire_setpos] component LAM_After_sapphire=L_monitor() SETTING [/usr/share/mcstas/3.0-dev/monitors/L_monitor.comp:65]");
  stracpy(_LAM_After_sapphire->_name, "LAM_After_sapphire", 16384);
  stracpy(_LAM_After_sapphire->_type, "L_monitor", 16384);
  _LAM_After_sapphire->_index=10;
  _LAM_After_sapphire->_parameters.nL = 100;
  #define nL (_LAM_After_sapphire->_parameters.nL)
  if("LAM_After_sapphire.out" && strlen("LAM_After_sapphire.out"))
    stracpy(_LAM_After_sapphire->_parameters.filename, "LAM_After_sapphire.out" ? "LAM_After_sapphire.out" : "", 16384);
  else 
  _LAM_After_sapphire->_parameters.filename[0]='\0';
  #define filename (_LAM_After_sapphire->_parameters.filename)
  _LAM_After_sapphire->_parameters.xmin = -0.1;
  #define xmin (_LAM_After_sapphire->_parameters.xmin)
  _LAM_After_sapphire->_parameters.xmax = 0.1;
  #define xmax (_LAM_After_sapphire->_parameters.xmax)
  _LAM_After_sapphire->_parameters.ymin = -0.1;
  #define ymin (_LAM_After_sapphire->_parameters.ymin)
  _LAM_After_sapphire->_parameters.ymax = 0.1;
  #define ymax (_LAM_After_sapphire->_parameters.ymax)
  _LAM_After_sapphire->_parameters.xwidth = 0;
  #define xwidth (_LAM_After_sapphire->_parameters.xwidth)
  _LAM_After_sapphire->_parameters.yheight = 0;
  #define yheight (_LAM_After_sapphire->_parameters.yheight)
  _LAM_After_sapphire->_parameters.Lmin = instrument->_parameters._source_lam_min;
  #define Lmin (_LAM_After_sapphire->_parameters.Lmin)
  _LAM_After_sapphire->_parameters.Lmax = instrument->_parameters._source_lam_max;
  #define Lmax (_LAM_After_sapphire->_parameters.Lmax)
  _LAM_After_sapphire->_parameters.restore_neutron = 0;
  #define restore_neutron (_LAM_After_sapphire->_parameters.restore_neutron)

  #define L_N (_LAM_After_sapphire->_parameters.L_N)
  #define L_p (_LAM_After_sapphire->_parameters.L_p)
  #define L_p2 (_LAM_After_sapphire->_parameters.L_p2)

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
  /* component LAM_After_sapphire=L_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Window_after_filter->_rotation_absolute, _LAM_After_sapphire->_rotation_absolute);
    rot_transpose(_Window_after_filter->_rotation_absolute, tr1);
    rot_mul(_LAM_After_sapphire->_rotation_absolute, tr1, _LAM_After_sapphire->_rotation_relative);
    _LAM_After_sapphire->_rotation_is_identity =  rot_test_identity(_LAM_After_sapphire->_rotation_relative);
    tc1 = coords_set(
      0, 0, 0);
    rot_transpose(_Window_after_filter->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _LAM_After_sapphire->_position_absolute = coords_add(_Window_after_filter->_position_absolute, tc2);
    tc1 = coords_sub(_Window_after_filter->_position_absolute, _LAM_After_sapphire->_position_absolute);
    _LAM_After_sapphire->_position_relative = rot_apply(_LAM_After_sapphire->_rotation_absolute, tc1);
  } /* LAM_After_sapphire=L_monitor() AT ROTATED */
  DEBUG_COMPONENT("LAM_After_sapphire", _LAM_After_sapphire->_position_absolute, _LAM_After_sapphire->_rotation_absolute);
  instrument->_position_absolute[10] = _LAM_After_sapphire->_position_absolute;
  instrument->_position_relative[10] = _LAM_After_sapphire->_position_relative;
  instrument->counter_N[10]  = instrument->counter_P[10] = instrument->counter_P2[10] = 0;
  instrument->counter_AbsorbProp[10]= 0;
  return(0);
} /* _LAM_After_sapphire_setpos */

/* component PSD_After_sapphire=PSD_monitor() SETTING, POSITION/ROTATION */
int _PSD_After_sapphire_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_PSD_After_sapphire_setpos] component PSD_After_sapphire=PSD_monitor() SETTING [/usr/share/mcstas/3.0-dev/monitors/PSD_monitor.comp:64]");
  stracpy(_PSD_After_sapphire->_name, "PSD_After_sapphire", 16384);
  stracpy(_PSD_After_sapphire->_type, "PSD_monitor", 16384);
  _PSD_After_sapphire->_index=11;
  _PSD_After_sapphire->_parameters.nx = 100;
  #define nx (_PSD_After_sapphire->_parameters.nx)
  _PSD_After_sapphire->_parameters.ny = 100;
  #define ny (_PSD_After_sapphire->_parameters.ny)
  if("PSD_After_sapphire" && strlen("PSD_After_sapphire"))
    stracpy(_PSD_After_sapphire->_parameters.filename, "PSD_After_sapphire" ? "PSD_After_sapphire" : "", 16384);
  else 
  _PSD_After_sapphire->_parameters.filename[0]='\0';
  #define filename (_PSD_After_sapphire->_parameters.filename)
  _PSD_After_sapphire->_parameters.xmin = -0.1;
  #define xmin (_PSD_After_sapphire->_parameters.xmin)
  _PSD_After_sapphire->_parameters.xmax = 0.1;
  #define xmax (_PSD_After_sapphire->_parameters.xmax)
  _PSD_After_sapphire->_parameters.ymin = -0.1;
  #define ymin (_PSD_After_sapphire->_parameters.ymin)
  _PSD_After_sapphire->_parameters.ymax = 0.1;
  #define ymax (_PSD_After_sapphire->_parameters.ymax)
  _PSD_After_sapphire->_parameters.xwidth = 0;
  #define xwidth (_PSD_After_sapphire->_parameters.xwidth)
  _PSD_After_sapphire->_parameters.yheight = 0;
  #define yheight (_PSD_After_sapphire->_parameters.yheight)
  _PSD_After_sapphire->_parameters.restore_neutron = 0;
  #define restore_neutron (_PSD_After_sapphire->_parameters.restore_neutron)

  #define PSD_N (_PSD_After_sapphire->_parameters.PSD_N)
  #define PSD_p (_PSD_After_sapphire->_parameters.PSD_p)
  #define PSD_p2 (_PSD_After_sapphire->_parameters.PSD_p2)

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
  /* component PSD_After_sapphire=PSD_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Window_after_filter->_rotation_absolute, _PSD_After_sapphire->_rotation_absolute);
    rot_transpose(_LAM_After_sapphire->_rotation_absolute, tr1);
    rot_mul(_PSD_After_sapphire->_rotation_absolute, tr1, _PSD_After_sapphire->_rotation_relative);
    _PSD_After_sapphire->_rotation_is_identity =  rot_test_identity(_PSD_After_sapphire->_rotation_relative);
    tc1 = coords_set(
      0, 0, 0);
    rot_transpose(_Window_after_filter->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _PSD_After_sapphire->_position_absolute = coords_add(_Window_after_filter->_position_absolute, tc2);
    tc1 = coords_sub(_LAM_After_sapphire->_position_absolute, _PSD_After_sapphire->_position_absolute);
    _PSD_After_sapphire->_position_relative = rot_apply(_PSD_After_sapphire->_rotation_absolute, tc1);
  } /* PSD_After_sapphire=PSD_monitor() AT ROTATED */
  DEBUG_COMPONENT("PSD_After_sapphire", _PSD_After_sapphire->_position_absolute, _PSD_After_sapphire->_rotation_absolute);
  instrument->_position_absolute[11] = _PSD_After_sapphire->_position_absolute;
  instrument->_position_relative[11] = _PSD_After_sapphire->_position_relative;
  instrument->counter_N[11]  = instrument->counter_P[11] = instrument->counter_P2[11] = 0;
  instrument->counter_AbsorbProp[11]= 0;
  return(0);
} /* _PSD_After_sapphire_setpos */

/* component HighResOutlet=Slit() SETTING, POSITION/ROTATION */
int _HighResOutlet_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_HighResOutlet_setpos] component HighResOutlet=Slit() SETTING [/usr/share/mcstas/3.0-dev/optics/Slit.comp:47]");
  stracpy(_HighResOutlet->_name, "HighResOutlet", 16384);
  stracpy(_HighResOutlet->_type, "Slit", 16384);
  _HighResOutlet->_index=12;
  _HighResOutlet->_parameters.xmin = -0.0325;
  #define xmin (_HighResOutlet->_parameters.xmin)
  _HighResOutlet->_parameters.xmax = 0.0325;
  #define xmax (_HighResOutlet->_parameters.xmax)
  _HighResOutlet->_parameters.ymin = -0.05;
  #define ymin (_HighResOutlet->_parameters.ymin)
  _HighResOutlet->_parameters.ymax = 0.05;
  #define ymax (_HighResOutlet->_parameters.ymax)
  _HighResOutlet->_parameters.radius = 0;
  #define radius (_HighResOutlet->_parameters.radius)
  _HighResOutlet->_parameters.xwidth = 0;
  #define xwidth (_HighResOutlet->_parameters.xwidth)
  _HighResOutlet->_parameters.yheight = 0;
  #define yheight (_HighResOutlet->_parameters.yheight)

  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef radius
  #undef xwidth
  #undef yheight
  /* component HighResOutlet=Slit() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Reactorbeam->_rotation_absolute, _HighResOutlet->_rotation_absolute);
    rot_transpose(_PSD_After_sapphire->_rotation_absolute, tr1);
    rot_mul(_HighResOutlet->_rotation_absolute, tr1, _HighResOutlet->_rotation_relative);
    _HighResOutlet->_rotation_is_identity =  rot_test_identity(_HighResOutlet->_rotation_relative);
    tc1 = coords_set(
      0, 0, 4.84055);
    rot_transpose(_Reactorbeam->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _HighResOutlet->_position_absolute = coords_add(_Reactorbeam->_position_absolute, tc2);
    tc1 = coords_sub(_PSD_After_sapphire->_position_absolute, _HighResOutlet->_position_absolute);
    _HighResOutlet->_position_relative = rot_apply(_HighResOutlet->_rotation_absolute, tc1);
  } /* HighResOutlet=Slit() AT ROTATED */
  DEBUG_COMPONENT("HighResOutlet", _HighResOutlet->_position_absolute, _HighResOutlet->_rotation_absolute);
  instrument->_position_absolute[12] = _HighResOutlet->_position_absolute;
  instrument->_position_relative[12] = _HighResOutlet->_position_relative;
  instrument->counter_N[12]  = instrument->counter_P[12] = instrument->counter_P2[12] = 0;
  instrument->counter_AbsorbProp[12]= 0;
  return(0);
} /* _HighResOutlet_setpos */

/* component HighIntensityOutlet=Slit() SETTING, POSITION/ROTATION */
int _HighIntensityOutlet_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_HighIntensityOutlet_setpos] component HighIntensityOutlet=Slit() SETTING [/usr/share/mcstas/3.0-dev/optics/Slit.comp:47]");
  stracpy(_HighIntensityOutlet->_name, "HighIntensityOutlet", 16384);
  stracpy(_HighIntensityOutlet->_type, "Slit", 16384);
  _HighIntensityOutlet->_index=13;
  _HighIntensityOutlet->_parameters.xmin = -0.05;
  #define xmin (_HighIntensityOutlet->_parameters.xmin)
  _HighIntensityOutlet->_parameters.xmax = 0.05;
  #define xmax (_HighIntensityOutlet->_parameters.xmax)
  _HighIntensityOutlet->_parameters.ymin = -0.05;
  #define ymin (_HighIntensityOutlet->_parameters.ymin)
  _HighIntensityOutlet->_parameters.ymax = 0.05;
  #define ymax (_HighIntensityOutlet->_parameters.ymax)
  _HighIntensityOutlet->_parameters.radius = 0;
  #define radius (_HighIntensityOutlet->_parameters.radius)
  _HighIntensityOutlet->_parameters.xwidth = 0;
  #define xwidth (_HighIntensityOutlet->_parameters.xwidth)
  _HighIntensityOutlet->_parameters.yheight = 0;
  #define yheight (_HighIntensityOutlet->_parameters.yheight)

  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef radius
  #undef xwidth
  #undef yheight
  /* component HighIntensityOutlet=Slit() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Reactorbeam->_rotation_absolute, _HighIntensityOutlet->_rotation_absolute);
    rot_transpose(_HighResOutlet->_rotation_absolute, tr1);
    rot_mul(_HighIntensityOutlet->_rotation_absolute, tr1, _HighIntensityOutlet->_rotation_relative);
    _HighIntensityOutlet->_rotation_is_identity =  rot_test_identity(_HighIntensityOutlet->_rotation_relative);
    tc1 = coords_set(
      0, 0, 4.84055);
    rot_transpose(_Reactorbeam->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _HighIntensityOutlet->_position_absolute = coords_add(_Reactorbeam->_position_absolute, tc2);
    tc1 = coords_sub(_HighResOutlet->_position_absolute, _HighIntensityOutlet->_position_absolute);
    _HighIntensityOutlet->_position_relative = rot_apply(_HighIntensityOutlet->_rotation_absolute, tc1);
  } /* HighIntensityOutlet=Slit() AT ROTATED */
  DEBUG_COMPONENT("HighIntensityOutlet", _HighIntensityOutlet->_position_absolute, _HighIntensityOutlet->_rotation_absolute);
  instrument->_position_absolute[13] = _HighIntensityOutlet->_position_absolute;
  instrument->_position_relative[13] = _HighIntensityOutlet->_position_relative;
  instrument->counter_N[13]  = instrument->counter_P[13] = instrument->counter_P2[13] = 0;
  instrument->counter_AbsorbProp[13]= 0;
  return(0);
} /* _HighIntensityOutlet_setpos */

/* component PSD_After_Outlet=PSD_monitor() SETTING, POSITION/ROTATION */
int _PSD_After_Outlet_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_PSD_After_Outlet_setpos] component PSD_After_Outlet=PSD_monitor() SETTING [/usr/share/mcstas/3.0-dev/monitors/PSD_monitor.comp:64]");
  stracpy(_PSD_After_Outlet->_name, "PSD_After_Outlet", 16384);
  stracpy(_PSD_After_Outlet->_type, "PSD_monitor", 16384);
  _PSD_After_Outlet->_index=14;
  _PSD_After_Outlet->_parameters.nx = 100;
  #define nx (_PSD_After_Outlet->_parameters.nx)
  _PSD_After_Outlet->_parameters.ny = 100;
  #define ny (_PSD_After_Outlet->_parameters.ny)
  if("PSD_After_Outlet" && strlen("PSD_After_Outlet"))
    stracpy(_PSD_After_Outlet->_parameters.filename, "PSD_After_Outlet" ? "PSD_After_Outlet" : "", 16384);
  else 
  _PSD_After_Outlet->_parameters.filename[0]='\0';
  #define filename (_PSD_After_Outlet->_parameters.filename)
  _PSD_After_Outlet->_parameters.xmin = -0.1;
  #define xmin (_PSD_After_Outlet->_parameters.xmin)
  _PSD_After_Outlet->_parameters.xmax = 0.1;
  #define xmax (_PSD_After_Outlet->_parameters.xmax)
  _PSD_After_Outlet->_parameters.ymin = -0.1;
  #define ymin (_PSD_After_Outlet->_parameters.ymin)
  _PSD_After_Outlet->_parameters.ymax = 0.1;
  #define ymax (_PSD_After_Outlet->_parameters.ymax)
  _PSD_After_Outlet->_parameters.xwidth = 0;
  #define xwidth (_PSD_After_Outlet->_parameters.xwidth)
  _PSD_After_Outlet->_parameters.yheight = 0;
  #define yheight (_PSD_After_Outlet->_parameters.yheight)
  _PSD_After_Outlet->_parameters.restore_neutron = 0;
  #define restore_neutron (_PSD_After_Outlet->_parameters.restore_neutron)

  #define PSD_N (_PSD_After_Outlet->_parameters.PSD_N)
  #define PSD_p (_PSD_After_Outlet->_parameters.PSD_p)
  #define PSD_p2 (_PSD_After_Outlet->_parameters.PSD_p2)

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
  /* component PSD_After_Outlet=PSD_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Reactorbeam->_rotation_absolute, _PSD_After_Outlet->_rotation_absolute);
    rot_transpose(_HighIntensityOutlet->_rotation_absolute, tr1);
    rot_mul(_PSD_After_Outlet->_rotation_absolute, tr1, _PSD_After_Outlet->_rotation_relative);
    _PSD_After_Outlet->_rotation_is_identity =  rot_test_identity(_PSD_After_Outlet->_rotation_relative);
    tc1 = coords_set(
      0, 0, 4.84055);
    rot_transpose(_Reactorbeam->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _PSD_After_Outlet->_position_absolute = coords_add(_Reactorbeam->_position_absolute, tc2);
    tc1 = coords_sub(_HighIntensityOutlet->_position_absolute, _PSD_After_Outlet->_position_absolute);
    _PSD_After_Outlet->_position_relative = rot_apply(_PSD_After_Outlet->_rotation_absolute, tc1);
  } /* PSD_After_Outlet=PSD_monitor() AT ROTATED */
  DEBUG_COMPONENT("PSD_After_Outlet", _PSD_After_Outlet->_position_absolute, _PSD_After_Outlet->_rotation_absolute);
  instrument->_position_absolute[14] = _PSD_After_Outlet->_position_absolute;
  instrument->_position_relative[14] = _PSD_After_Outlet->_position_relative;
  instrument->counter_N[14]  = instrument->counter_P[14] = instrument->counter_P2[14] = 0;
  instrument->counter_AbsorbProp[14]= 0;
  return(0);
} /* _PSD_After_Outlet_setpos */

/* component Mono_axis=Arm() SETTING, POSITION/ROTATION */
int _Mono_axis_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Mono_axis_setpos] component Mono_axis=Arm() SETTING [Arm:0]");
  stracpy(_Mono_axis->_name, "Mono_axis", 16384);
  stracpy(_Mono_axis->_type, "Arm", 16384);
  _Mono_axis->_index=15;
  /* component Mono_axis=Arm() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (instrument->_parameters._mono_takeoff / 2.0)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _Reactorbeam->_rotation_absolute, _Mono_axis->_rotation_absolute);
    rot_transpose(_PSD_After_Outlet->_rotation_absolute, tr1);
    rot_mul(_Mono_axis->_rotation_absolute, tr1, _Mono_axis->_rotation_relative);
    _Mono_axis->_rotation_is_identity =  rot_test_identity(_Mono_axis->_rotation_relative);
    tc1 = coords_set(
      0 + instrument->_parameters._mono_dx, 0 + instrument->_parameters._mono_dy, 5.140 + instrument->_parameters._mono_dz);
    rot_transpose(_Reactorbeam->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Mono_axis->_position_absolute = coords_add(_Reactorbeam->_position_absolute, tc2);
    tc1 = coords_sub(_PSD_After_Outlet->_position_absolute, _Mono_axis->_position_absolute);
    _Mono_axis->_position_relative = rot_apply(_Mono_axis->_rotation_absolute, tc1);
  } /* Mono_axis=Arm() AT ROTATED */
  DEBUG_COMPONENT("Mono_axis", _Mono_axis->_position_absolute, _Mono_axis->_rotation_absolute);
  instrument->_position_absolute[15] = _Mono_axis->_position_absolute;
  instrument->_position_relative[15] = _Mono_axis->_position_relative;
  instrument->counter_N[15]  = instrument->counter_P[15] = instrument->counter_P2[15] = 0;
  instrument->counter_AbsorbProp[15]= 0;
  return(0);
} /* _Mono_axis_setpos */

/* component Blade_1=Monochromator_curved() SETTING, POSITION/ROTATION */
int _Blade_1_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Blade_1_setpos] component Blade_1=Monochromator_curved() SETTING [/usr/share/mcstas/3.0-dev/optics/Monochromator_curved.comp:148]");
  stracpy(_Blade_1->_name, "Blade_1", 16384);
  stracpy(_Blade_1->_type, "Monochromator_curved", 16384);
  _Blade_1->_index=16;
  if("NULL" && strlen("NULL"))
    stracpy(_Blade_1->_parameters.reflect, "NULL" ? "NULL" : "", 16384);
  else 
  _Blade_1->_parameters.reflect[0]='\0';
  #define reflect (_Blade_1->_parameters.reflect)
  if("NULL" && strlen("NULL"))
    stracpy(_Blade_1->_parameters.transmit, "NULL" ? "NULL" : "", 16384);
  else 
  _Blade_1->_parameters.transmit[0]='\0';
  #define transmit (_Blade_1->_parameters.transmit)
  _Blade_1->_parameters.zwidth = 0.22 / 51.0;
  #define zwidth (_Blade_1->_parameters.zwidth)
  _Blade_1->_parameters.yheight = 0.1395 / 9.0;
  #define yheight (_Blade_1->_parameters.yheight)
  _Blade_1->_parameters.gap = 0.0;
  #define gap (_Blade_1->_parameters.gap)
  _Blade_1->_parameters.NH = 51;
  #define NH (_Blade_1->_parameters.NH)
  _Blade_1->_parameters.NV = 9;
  #define NV (_Blade_1->_parameters.NV)
  _Blade_1->_parameters.mosaich = instrument->_parameters._mono_mosh;
  #define mosaich (_Blade_1->_parameters.mosaich)
  _Blade_1->_parameters.mosaicv = instrument->_parameters._mono_mosv;
  #define mosaicv (_Blade_1->_parameters.mosaicv)
  _Blade_1->_parameters.r0 = 1.0;
  #define r0 (_Blade_1->_parameters.r0)
  _Blade_1->_parameters.t0 = 1.0;
  #define t0 (_Blade_1->_parameters.t0)
  _Blade_1->_parameters.Q = 1.8734;
  #define Q (_Blade_1->_parameters.Q)
  _Blade_1->_parameters.RV = instrument->_parameters._mono_r_v;
  #define RV (_Blade_1->_parameters.RV)
  _Blade_1->_parameters.RH = instrument->_parameters._mono_r_h;
  #define RH (_Blade_1->_parameters.RH)
  _Blade_1->_parameters.DM = mono_d;
  #define DM (_Blade_1->_parameters.DM)
  _Blade_1->_parameters.mosaic = 0;
  #define mosaic (_Blade_1->_parameters.mosaic)
  _Blade_1->_parameters.width = 0;
  #define width (_Blade_1->_parameters.width)
  _Blade_1->_parameters.height = 0;
  #define height (_Blade_1->_parameters.height)
  _Blade_1->_parameters.verbose = 0;
  #define verbose (_Blade_1->_parameters.verbose)
  _Blade_1->_parameters.order = 0;
  #define order (_Blade_1->_parameters.order)

  #define mos_rms_y (_Blade_1->_parameters.mos_rms_y)
  #define mos_rms_z (_Blade_1->_parameters.mos_rms_z)
  #define mos_rms_max (_Blade_1->_parameters.mos_rms_max)
  #define mono_Q (_Blade_1->_parameters.mono_Q)
  #define SlabWidth (_Blade_1->_parameters.SlabWidth)
  #define SlabHeight (_Blade_1->_parameters.SlabHeight)
  #define rTable (_Blade_1->_parameters.rTable)
  #define tTable (_Blade_1->_parameters.tTable)
  #define row (_Blade_1->_parameters.row)
  #define col (_Blade_1->_parameters.col)
  #define tiltH (_Blade_1->_parameters.tiltH)
  #define tiltV (_Blade_1->_parameters.tiltV)

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
  /* component Blade_1=Monochromator_curved() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Mono_axis->_rotation_absolute, _Blade_1->_rotation_absolute);
    rot_transpose(_Mono_axis->_rotation_absolute, tr1);
    rot_mul(_Blade_1->_rotation_absolute, tr1, _Blade_1->_rotation_relative);
    _Blade_1->_rotation_is_identity =  rot_test_identity(_Blade_1->_rotation_relative);
    tc1 = coords_set(
      - start_wafer_pos, 0, 0);
    rot_transpose(_Mono_axis->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Blade_1->_position_absolute = coords_add(_Mono_axis->_position_absolute, tc2);
    tc1 = coords_sub(_Mono_axis->_position_absolute, _Blade_1->_position_absolute);
    _Blade_1->_position_relative = rot_apply(_Blade_1->_rotation_absolute, tc1);
  } /* Blade_1=Monochromator_curved() AT ROTATED */
  DEBUG_COMPONENT("Blade_1", _Blade_1->_position_absolute, _Blade_1->_rotation_absolute);
  instrument->_position_absolute[16] = _Blade_1->_position_absolute;
  instrument->_position_relative[16] = _Blade_1->_position_relative;
  instrument->counter_N[16]  = instrument->counter_P[16] = instrument->counter_P2[16] = 0;
  instrument->counter_AbsorbProp[16]= 0;
  return(0);
} /* _Blade_1_setpos */

/* component Blade_2=Monochromator_curved() SETTING, POSITION/ROTATION */
int _Blade_2_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Blade_2_setpos] component Blade_2=Monochromator_curved() SETTING [/usr/share/mcstas/3.0-dev/optics/Monochromator_curved.comp:148]");
  stracpy(_Blade_2->_name, "Blade_2", 16384);
  stracpy(_Blade_2->_type, "Monochromator_curved", 16384);
  _Blade_2->_index=17;
  if("NULL" && strlen("NULL"))
    stracpy(_Blade_2->_parameters.reflect, "NULL" ? "NULL" : "", 16384);
  else 
  _Blade_2->_parameters.reflect[0]='\0';
  #define reflect (_Blade_2->_parameters.reflect)
  if("NULL" && strlen("NULL"))
    stracpy(_Blade_2->_parameters.transmit, "NULL" ? "NULL" : "", 16384);
  else 
  _Blade_2->_parameters.transmit[0]='\0';
  #define transmit (_Blade_2->_parameters.transmit)
  _Blade_2->_parameters.zwidth = 0.22 / 51.0;
  #define zwidth (_Blade_2->_parameters.zwidth)
  _Blade_2->_parameters.yheight = 0.1395 / 9.0;
  #define yheight (_Blade_2->_parameters.yheight)
  _Blade_2->_parameters.gap = 0.0;
  #define gap (_Blade_2->_parameters.gap)
  _Blade_2->_parameters.NH = 51;
  #define NH (_Blade_2->_parameters.NH)
  _Blade_2->_parameters.NV = 9;
  #define NV (_Blade_2->_parameters.NV)
  _Blade_2->_parameters.mosaich = instrument->_parameters._mono_mosh;
  #define mosaich (_Blade_2->_parameters.mosaich)
  _Blade_2->_parameters.mosaicv = instrument->_parameters._mono_mosv;
  #define mosaicv (_Blade_2->_parameters.mosaicv)
  _Blade_2->_parameters.r0 = 1.0;
  #define r0 (_Blade_2->_parameters.r0)
  _Blade_2->_parameters.t0 = 1.0;
  #define t0 (_Blade_2->_parameters.t0)
  _Blade_2->_parameters.Q = 1.8734;
  #define Q (_Blade_2->_parameters.Q)
  _Blade_2->_parameters.RV = instrument->_parameters._mono_r_v;
  #define RV (_Blade_2->_parameters.RV)
  _Blade_2->_parameters.RH = instrument->_parameters._mono_r_h;
  #define RH (_Blade_2->_parameters.RH)
  _Blade_2->_parameters.DM = mono_d;
  #define DM (_Blade_2->_parameters.DM)
  _Blade_2->_parameters.mosaic = 0;
  #define mosaic (_Blade_2->_parameters.mosaic)
  _Blade_2->_parameters.width = 0;
  #define width (_Blade_2->_parameters.width)
  _Blade_2->_parameters.height = 0;
  #define height (_Blade_2->_parameters.height)
  _Blade_2->_parameters.verbose = 0;
  #define verbose (_Blade_2->_parameters.verbose)
  _Blade_2->_parameters.order = 0;
  #define order (_Blade_2->_parameters.order)

  #define mos_rms_y (_Blade_2->_parameters.mos_rms_y)
  #define mos_rms_z (_Blade_2->_parameters.mos_rms_z)
  #define mos_rms_max (_Blade_2->_parameters.mos_rms_max)
  #define mono_Q (_Blade_2->_parameters.mono_Q)
  #define SlabWidth (_Blade_2->_parameters.SlabWidth)
  #define SlabHeight (_Blade_2->_parameters.SlabHeight)
  #define rTable (_Blade_2->_parameters.rTable)
  #define tTable (_Blade_2->_parameters.tTable)
  #define row (_Blade_2->_parameters.row)
  #define col (_Blade_2->_parameters.col)
  #define tiltH (_Blade_2->_parameters.tiltH)
  #define tiltV (_Blade_2->_parameters.tiltV)

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
  /* component Blade_2=Monochromator_curved() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Blade_1->_rotation_absolute, _Blade_2->_rotation_absolute);
    rot_transpose(_Blade_1->_rotation_absolute, tr1);
    rot_mul(_Blade_2->_rotation_absolute, tr1, _Blade_2->_rotation_relative);
    _Blade_2->_rotation_is_identity =  rot_test_identity(_Blade_2->_rotation_relative);
    tc1 = coords_set(
      wafer_d, 0, 0);
    rot_transpose(_Blade_1->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Blade_2->_position_absolute = coords_add(_Blade_1->_position_absolute, tc2);
    tc1 = coords_sub(_Blade_1->_position_absolute, _Blade_2->_position_absolute);
    _Blade_2->_position_relative = rot_apply(_Blade_2->_rotation_absolute, tc1);
  } /* Blade_2=Monochromator_curved() AT ROTATED */
  DEBUG_COMPONENT("Blade_2", _Blade_2->_position_absolute, _Blade_2->_rotation_absolute);
  instrument->_position_absolute[17] = _Blade_2->_position_absolute;
  instrument->_position_relative[17] = _Blade_2->_position_relative;
  instrument->counter_N[17]  = instrument->counter_P[17] = instrument->counter_P2[17] = 0;
  instrument->counter_AbsorbProp[17]= 0;
  return(0);
} /* _Blade_2_setpos */

/* component Blade_3=Monochromator_curved() SETTING, POSITION/ROTATION */
int _Blade_3_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Blade_3_setpos] component Blade_3=Monochromator_curved() SETTING [/usr/share/mcstas/3.0-dev/optics/Monochromator_curved.comp:148]");
  stracpy(_Blade_3->_name, "Blade_3", 16384);
  stracpy(_Blade_3->_type, "Monochromator_curved", 16384);
  _Blade_3->_index=18;
  if("NULL" && strlen("NULL"))
    stracpy(_Blade_3->_parameters.reflect, "NULL" ? "NULL" : "", 16384);
  else 
  _Blade_3->_parameters.reflect[0]='\0';
  #define reflect (_Blade_3->_parameters.reflect)
  if("NULL" && strlen("NULL"))
    stracpy(_Blade_3->_parameters.transmit, "NULL" ? "NULL" : "", 16384);
  else 
  _Blade_3->_parameters.transmit[0]='\0';
  #define transmit (_Blade_3->_parameters.transmit)
  _Blade_3->_parameters.zwidth = 0.22 / 51.0;
  #define zwidth (_Blade_3->_parameters.zwidth)
  _Blade_3->_parameters.yheight = 0.1395 / 9.0;
  #define yheight (_Blade_3->_parameters.yheight)
  _Blade_3->_parameters.gap = 0.0;
  #define gap (_Blade_3->_parameters.gap)
  _Blade_3->_parameters.NH = 51;
  #define NH (_Blade_3->_parameters.NH)
  _Blade_3->_parameters.NV = 9;
  #define NV (_Blade_3->_parameters.NV)
  _Blade_3->_parameters.mosaich = instrument->_parameters._mono_mosh;
  #define mosaich (_Blade_3->_parameters.mosaich)
  _Blade_3->_parameters.mosaicv = instrument->_parameters._mono_mosv;
  #define mosaicv (_Blade_3->_parameters.mosaicv)
  _Blade_3->_parameters.r0 = 1.0;
  #define r0 (_Blade_3->_parameters.r0)
  _Blade_3->_parameters.t0 = 1.0;
  #define t0 (_Blade_3->_parameters.t0)
  _Blade_3->_parameters.Q = 1.8734;
  #define Q (_Blade_3->_parameters.Q)
  _Blade_3->_parameters.RV = instrument->_parameters._mono_r_v;
  #define RV (_Blade_3->_parameters.RV)
  _Blade_3->_parameters.RH = instrument->_parameters._mono_r_h;
  #define RH (_Blade_3->_parameters.RH)
  _Blade_3->_parameters.DM = mono_d;
  #define DM (_Blade_3->_parameters.DM)
  _Blade_3->_parameters.mosaic = 0;
  #define mosaic (_Blade_3->_parameters.mosaic)
  _Blade_3->_parameters.width = 0;
  #define width (_Blade_3->_parameters.width)
  _Blade_3->_parameters.height = 0;
  #define height (_Blade_3->_parameters.height)
  _Blade_3->_parameters.verbose = 0;
  #define verbose (_Blade_3->_parameters.verbose)
  _Blade_3->_parameters.order = 0;
  #define order (_Blade_3->_parameters.order)

  #define mos_rms_y (_Blade_3->_parameters.mos_rms_y)
  #define mos_rms_z (_Blade_3->_parameters.mos_rms_z)
  #define mos_rms_max (_Blade_3->_parameters.mos_rms_max)
  #define mono_Q (_Blade_3->_parameters.mono_Q)
  #define SlabWidth (_Blade_3->_parameters.SlabWidth)
  #define SlabHeight (_Blade_3->_parameters.SlabHeight)
  #define rTable (_Blade_3->_parameters.rTable)
  #define tTable (_Blade_3->_parameters.tTable)
  #define row (_Blade_3->_parameters.row)
  #define col (_Blade_3->_parameters.col)
  #define tiltH (_Blade_3->_parameters.tiltH)
  #define tiltV (_Blade_3->_parameters.tiltV)

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
  /* component Blade_3=Monochromator_curved() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Blade_2->_rotation_absolute, _Blade_3->_rotation_absolute);
    rot_transpose(_Blade_2->_rotation_absolute, tr1);
    rot_mul(_Blade_3->_rotation_absolute, tr1, _Blade_3->_rotation_relative);
    _Blade_3->_rotation_is_identity =  rot_test_identity(_Blade_3->_rotation_relative);
    tc1 = coords_set(
      wafer_d, 0, 0);
    rot_transpose(_Blade_2->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Blade_3->_position_absolute = coords_add(_Blade_2->_position_absolute, tc2);
    tc1 = coords_sub(_Blade_2->_position_absolute, _Blade_3->_position_absolute);
    _Blade_3->_position_relative = rot_apply(_Blade_3->_rotation_absolute, tc1);
  } /* Blade_3=Monochromator_curved() AT ROTATED */
  DEBUG_COMPONENT("Blade_3", _Blade_3->_position_absolute, _Blade_3->_rotation_absolute);
  instrument->_position_absolute[18] = _Blade_3->_position_absolute;
  instrument->_position_relative[18] = _Blade_3->_position_relative;
  instrument->counter_N[18]  = instrument->counter_P[18] = instrument->counter_P2[18] = 0;
  instrument->counter_AbsorbProp[18]= 0;
  return(0);
} /* _Blade_3_setpos */

/* component Blade_4=Monochromator_curved() SETTING, POSITION/ROTATION */
int _Blade_4_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Blade_4_setpos] component Blade_4=Monochromator_curved() SETTING [/usr/share/mcstas/3.0-dev/optics/Monochromator_curved.comp:148]");
  stracpy(_Blade_4->_name, "Blade_4", 16384);
  stracpy(_Blade_4->_type, "Monochromator_curved", 16384);
  _Blade_4->_index=19;
  if("NULL" && strlen("NULL"))
    stracpy(_Blade_4->_parameters.reflect, "NULL" ? "NULL" : "", 16384);
  else 
  _Blade_4->_parameters.reflect[0]='\0';
  #define reflect (_Blade_4->_parameters.reflect)
  if("NULL" && strlen("NULL"))
    stracpy(_Blade_4->_parameters.transmit, "NULL" ? "NULL" : "", 16384);
  else 
  _Blade_4->_parameters.transmit[0]='\0';
  #define transmit (_Blade_4->_parameters.transmit)
  _Blade_4->_parameters.zwidth = 0.22 / 51.0;
  #define zwidth (_Blade_4->_parameters.zwidth)
  _Blade_4->_parameters.yheight = 0.1395 / 9.0;
  #define yheight (_Blade_4->_parameters.yheight)
  _Blade_4->_parameters.gap = 0.0;
  #define gap (_Blade_4->_parameters.gap)
  _Blade_4->_parameters.NH = 51;
  #define NH (_Blade_4->_parameters.NH)
  _Blade_4->_parameters.NV = 9;
  #define NV (_Blade_4->_parameters.NV)
  _Blade_4->_parameters.mosaich = instrument->_parameters._mono_mosh;
  #define mosaich (_Blade_4->_parameters.mosaich)
  _Blade_4->_parameters.mosaicv = instrument->_parameters._mono_mosv;
  #define mosaicv (_Blade_4->_parameters.mosaicv)
  _Blade_4->_parameters.r0 = 1.0;
  #define r0 (_Blade_4->_parameters.r0)
  _Blade_4->_parameters.t0 = 1.0;
  #define t0 (_Blade_4->_parameters.t0)
  _Blade_4->_parameters.Q = 1.8734;
  #define Q (_Blade_4->_parameters.Q)
  _Blade_4->_parameters.RV = instrument->_parameters._mono_r_v;
  #define RV (_Blade_4->_parameters.RV)
  _Blade_4->_parameters.RH = instrument->_parameters._mono_r_h;
  #define RH (_Blade_4->_parameters.RH)
  _Blade_4->_parameters.DM = mono_d;
  #define DM (_Blade_4->_parameters.DM)
  _Blade_4->_parameters.mosaic = 0;
  #define mosaic (_Blade_4->_parameters.mosaic)
  _Blade_4->_parameters.width = 0;
  #define width (_Blade_4->_parameters.width)
  _Blade_4->_parameters.height = 0;
  #define height (_Blade_4->_parameters.height)
  _Blade_4->_parameters.verbose = 0;
  #define verbose (_Blade_4->_parameters.verbose)
  _Blade_4->_parameters.order = 0;
  #define order (_Blade_4->_parameters.order)

  #define mos_rms_y (_Blade_4->_parameters.mos_rms_y)
  #define mos_rms_z (_Blade_4->_parameters.mos_rms_z)
  #define mos_rms_max (_Blade_4->_parameters.mos_rms_max)
  #define mono_Q (_Blade_4->_parameters.mono_Q)
  #define SlabWidth (_Blade_4->_parameters.SlabWidth)
  #define SlabHeight (_Blade_4->_parameters.SlabHeight)
  #define rTable (_Blade_4->_parameters.rTable)
  #define tTable (_Blade_4->_parameters.tTable)
  #define row (_Blade_4->_parameters.row)
  #define col (_Blade_4->_parameters.col)
  #define tiltH (_Blade_4->_parameters.tiltH)
  #define tiltV (_Blade_4->_parameters.tiltV)

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
  /* component Blade_4=Monochromator_curved() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Blade_3->_rotation_absolute, _Blade_4->_rotation_absolute);
    rot_transpose(_Blade_3->_rotation_absolute, tr1);
    rot_mul(_Blade_4->_rotation_absolute, tr1, _Blade_4->_rotation_relative);
    _Blade_4->_rotation_is_identity =  rot_test_identity(_Blade_4->_rotation_relative);
    tc1 = coords_set(
      wafer_d, 0, 0);
    rot_transpose(_Blade_3->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Blade_4->_position_absolute = coords_add(_Blade_3->_position_absolute, tc2);
    tc1 = coords_sub(_Blade_3->_position_absolute, _Blade_4->_position_absolute);
    _Blade_4->_position_relative = rot_apply(_Blade_4->_rotation_absolute, tc1);
  } /* Blade_4=Monochromator_curved() AT ROTATED */
  DEBUG_COMPONENT("Blade_4", _Blade_4->_position_absolute, _Blade_4->_rotation_absolute);
  instrument->_position_absolute[19] = _Blade_4->_position_absolute;
  instrument->_position_relative[19] = _Blade_4->_position_relative;
  instrument->counter_N[19]  = instrument->counter_P[19] = instrument->counter_P2[19] = 0;
  instrument->counter_AbsorbProp[19]= 0;
  return(0);
} /* _Blade_4_setpos */

/* component Blade_5=Monochromator_curved() SETTING, POSITION/ROTATION */
int _Blade_5_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Blade_5_setpos] component Blade_5=Monochromator_curved() SETTING [/usr/share/mcstas/3.0-dev/optics/Monochromator_curved.comp:148]");
  stracpy(_Blade_5->_name, "Blade_5", 16384);
  stracpy(_Blade_5->_type, "Monochromator_curved", 16384);
  _Blade_5->_index=20;
  if("NULL" && strlen("NULL"))
    stracpy(_Blade_5->_parameters.reflect, "NULL" ? "NULL" : "", 16384);
  else 
  _Blade_5->_parameters.reflect[0]='\0';
  #define reflect (_Blade_5->_parameters.reflect)
  if("NULL" && strlen("NULL"))
    stracpy(_Blade_5->_parameters.transmit, "NULL" ? "NULL" : "", 16384);
  else 
  _Blade_5->_parameters.transmit[0]='\0';
  #define transmit (_Blade_5->_parameters.transmit)
  _Blade_5->_parameters.zwidth = 0.22 / 51.0;
  #define zwidth (_Blade_5->_parameters.zwidth)
  _Blade_5->_parameters.yheight = 0.1395 / 9.0;
  #define yheight (_Blade_5->_parameters.yheight)
  _Blade_5->_parameters.gap = 0.0;
  #define gap (_Blade_5->_parameters.gap)
  _Blade_5->_parameters.NH = 51;
  #define NH (_Blade_5->_parameters.NH)
  _Blade_5->_parameters.NV = 9;
  #define NV (_Blade_5->_parameters.NV)
  _Blade_5->_parameters.mosaich = instrument->_parameters._mono_mosh;
  #define mosaich (_Blade_5->_parameters.mosaich)
  _Blade_5->_parameters.mosaicv = instrument->_parameters._mono_mosv;
  #define mosaicv (_Blade_5->_parameters.mosaicv)
  _Blade_5->_parameters.r0 = 1.0;
  #define r0 (_Blade_5->_parameters.r0)
  _Blade_5->_parameters.t0 = 1.0;
  #define t0 (_Blade_5->_parameters.t0)
  _Blade_5->_parameters.Q = 1.8734;
  #define Q (_Blade_5->_parameters.Q)
  _Blade_5->_parameters.RV = instrument->_parameters._mono_r_v;
  #define RV (_Blade_5->_parameters.RV)
  _Blade_5->_parameters.RH = instrument->_parameters._mono_r_h;
  #define RH (_Blade_5->_parameters.RH)
  _Blade_5->_parameters.DM = mono_d;
  #define DM (_Blade_5->_parameters.DM)
  _Blade_5->_parameters.mosaic = 0;
  #define mosaic (_Blade_5->_parameters.mosaic)
  _Blade_5->_parameters.width = 0;
  #define width (_Blade_5->_parameters.width)
  _Blade_5->_parameters.height = 0;
  #define height (_Blade_5->_parameters.height)
  _Blade_5->_parameters.verbose = 0;
  #define verbose (_Blade_5->_parameters.verbose)
  _Blade_5->_parameters.order = 0;
  #define order (_Blade_5->_parameters.order)

  #define mos_rms_y (_Blade_5->_parameters.mos_rms_y)
  #define mos_rms_z (_Blade_5->_parameters.mos_rms_z)
  #define mos_rms_max (_Blade_5->_parameters.mos_rms_max)
  #define mono_Q (_Blade_5->_parameters.mono_Q)
  #define SlabWidth (_Blade_5->_parameters.SlabWidth)
  #define SlabHeight (_Blade_5->_parameters.SlabHeight)
  #define rTable (_Blade_5->_parameters.rTable)
  #define tTable (_Blade_5->_parameters.tTable)
  #define row (_Blade_5->_parameters.row)
  #define col (_Blade_5->_parameters.col)
  #define tiltH (_Blade_5->_parameters.tiltH)
  #define tiltV (_Blade_5->_parameters.tiltV)

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
  /* component Blade_5=Monochromator_curved() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Blade_4->_rotation_absolute, _Blade_5->_rotation_absolute);
    rot_transpose(_Blade_4->_rotation_absolute, tr1);
    rot_mul(_Blade_5->_rotation_absolute, tr1, _Blade_5->_rotation_relative);
    _Blade_5->_rotation_is_identity =  rot_test_identity(_Blade_5->_rotation_relative);
    tc1 = coords_set(
      wafer_d, 0, 0);
    rot_transpose(_Blade_4->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Blade_5->_position_absolute = coords_add(_Blade_4->_position_absolute, tc2);
    tc1 = coords_sub(_Blade_4->_position_absolute, _Blade_5->_position_absolute);
    _Blade_5->_position_relative = rot_apply(_Blade_5->_rotation_absolute, tc1);
  } /* Blade_5=Monochromator_curved() AT ROTATED */
  DEBUG_COMPONENT("Blade_5", _Blade_5->_position_absolute, _Blade_5->_rotation_absolute);
  instrument->_position_absolute[20] = _Blade_5->_position_absolute;
  instrument->_position_relative[20] = _Blade_5->_position_relative;
  instrument->counter_N[20]  = instrument->counter_P[20] = instrument->counter_P2[20] = 0;
  instrument->counter_AbsorbProp[20]= 0;
  return(0);
} /* _Blade_5_setpos */

/* component Blade_6=Monochromator_curved() SETTING, POSITION/ROTATION */
int _Blade_6_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Blade_6_setpos] component Blade_6=Monochromator_curved() SETTING [/usr/share/mcstas/3.0-dev/optics/Monochromator_curved.comp:148]");
  stracpy(_Blade_6->_name, "Blade_6", 16384);
  stracpy(_Blade_6->_type, "Monochromator_curved", 16384);
  _Blade_6->_index=21;
  if("NULL" && strlen("NULL"))
    stracpy(_Blade_6->_parameters.reflect, "NULL" ? "NULL" : "", 16384);
  else 
  _Blade_6->_parameters.reflect[0]='\0';
  #define reflect (_Blade_6->_parameters.reflect)
  if("NULL" && strlen("NULL"))
    stracpy(_Blade_6->_parameters.transmit, "NULL" ? "NULL" : "", 16384);
  else 
  _Blade_6->_parameters.transmit[0]='\0';
  #define transmit (_Blade_6->_parameters.transmit)
  _Blade_6->_parameters.zwidth = 0.22 / 51.0;
  #define zwidth (_Blade_6->_parameters.zwidth)
  _Blade_6->_parameters.yheight = 0.1395 / 9.0;
  #define yheight (_Blade_6->_parameters.yheight)
  _Blade_6->_parameters.gap = 0.0;
  #define gap (_Blade_6->_parameters.gap)
  _Blade_6->_parameters.NH = 51;
  #define NH (_Blade_6->_parameters.NH)
  _Blade_6->_parameters.NV = 9;
  #define NV (_Blade_6->_parameters.NV)
  _Blade_6->_parameters.mosaich = instrument->_parameters._mono_mosh;
  #define mosaich (_Blade_6->_parameters.mosaich)
  _Blade_6->_parameters.mosaicv = instrument->_parameters._mono_mosv;
  #define mosaicv (_Blade_6->_parameters.mosaicv)
  _Blade_6->_parameters.r0 = 1.0;
  #define r0 (_Blade_6->_parameters.r0)
  _Blade_6->_parameters.t0 = 1.0;
  #define t0 (_Blade_6->_parameters.t0)
  _Blade_6->_parameters.Q = 1.8734;
  #define Q (_Blade_6->_parameters.Q)
  _Blade_6->_parameters.RV = instrument->_parameters._mono_r_v;
  #define RV (_Blade_6->_parameters.RV)
  _Blade_6->_parameters.RH = instrument->_parameters._mono_r_h;
  #define RH (_Blade_6->_parameters.RH)
  _Blade_6->_parameters.DM = mono_d;
  #define DM (_Blade_6->_parameters.DM)
  _Blade_6->_parameters.mosaic = 0;
  #define mosaic (_Blade_6->_parameters.mosaic)
  _Blade_6->_parameters.width = 0;
  #define width (_Blade_6->_parameters.width)
  _Blade_6->_parameters.height = 0;
  #define height (_Blade_6->_parameters.height)
  _Blade_6->_parameters.verbose = 0;
  #define verbose (_Blade_6->_parameters.verbose)
  _Blade_6->_parameters.order = 0;
  #define order (_Blade_6->_parameters.order)

  #define mos_rms_y (_Blade_6->_parameters.mos_rms_y)
  #define mos_rms_z (_Blade_6->_parameters.mos_rms_z)
  #define mos_rms_max (_Blade_6->_parameters.mos_rms_max)
  #define mono_Q (_Blade_6->_parameters.mono_Q)
  #define SlabWidth (_Blade_6->_parameters.SlabWidth)
  #define SlabHeight (_Blade_6->_parameters.SlabHeight)
  #define rTable (_Blade_6->_parameters.rTable)
  #define tTable (_Blade_6->_parameters.tTable)
  #define row (_Blade_6->_parameters.row)
  #define col (_Blade_6->_parameters.col)
  #define tiltH (_Blade_6->_parameters.tiltH)
  #define tiltV (_Blade_6->_parameters.tiltV)

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
  /* component Blade_6=Monochromator_curved() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Blade_5->_rotation_absolute, _Blade_6->_rotation_absolute);
    rot_transpose(_Blade_5->_rotation_absolute, tr1);
    rot_mul(_Blade_6->_rotation_absolute, tr1, _Blade_6->_rotation_relative);
    _Blade_6->_rotation_is_identity =  rot_test_identity(_Blade_6->_rotation_relative);
    tc1 = coords_set(
      wafer_d, 0, 0);
    rot_transpose(_Blade_5->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Blade_6->_position_absolute = coords_add(_Blade_5->_position_absolute, tc2);
    tc1 = coords_sub(_Blade_5->_position_absolute, _Blade_6->_position_absolute);
    _Blade_6->_position_relative = rot_apply(_Blade_6->_rotation_absolute, tc1);
  } /* Blade_6=Monochromator_curved() AT ROTATED */
  DEBUG_COMPONENT("Blade_6", _Blade_6->_position_absolute, _Blade_6->_rotation_absolute);
  instrument->_position_absolute[21] = _Blade_6->_position_absolute;
  instrument->_position_relative[21] = _Blade_6->_position_relative;
  instrument->counter_N[21]  = instrument->counter_P[21] = instrument->counter_P2[21] = 0;
  instrument->counter_AbsorbProp[21]= 0;
  return(0);
} /* _Blade_6_setpos */

/* component Blade_7=Monochromator_curved() SETTING, POSITION/ROTATION */
int _Blade_7_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Blade_7_setpos] component Blade_7=Monochromator_curved() SETTING [/usr/share/mcstas/3.0-dev/optics/Monochromator_curved.comp:148]");
  stracpy(_Blade_7->_name, "Blade_7", 16384);
  stracpy(_Blade_7->_type, "Monochromator_curved", 16384);
  _Blade_7->_index=22;
  if("NULL" && strlen("NULL"))
    stracpy(_Blade_7->_parameters.reflect, "NULL" ? "NULL" : "", 16384);
  else 
  _Blade_7->_parameters.reflect[0]='\0';
  #define reflect (_Blade_7->_parameters.reflect)
  if("NULL" && strlen("NULL"))
    stracpy(_Blade_7->_parameters.transmit, "NULL" ? "NULL" : "", 16384);
  else 
  _Blade_7->_parameters.transmit[0]='\0';
  #define transmit (_Blade_7->_parameters.transmit)
  _Blade_7->_parameters.zwidth = 0.22 / 51.0;
  #define zwidth (_Blade_7->_parameters.zwidth)
  _Blade_7->_parameters.yheight = 0.1395 / 9.0;
  #define yheight (_Blade_7->_parameters.yheight)
  _Blade_7->_parameters.gap = 0.0;
  #define gap (_Blade_7->_parameters.gap)
  _Blade_7->_parameters.NH = 51;
  #define NH (_Blade_7->_parameters.NH)
  _Blade_7->_parameters.NV = 9;
  #define NV (_Blade_7->_parameters.NV)
  _Blade_7->_parameters.mosaich = instrument->_parameters._mono_mosh;
  #define mosaich (_Blade_7->_parameters.mosaich)
  _Blade_7->_parameters.mosaicv = instrument->_parameters._mono_mosv;
  #define mosaicv (_Blade_7->_parameters.mosaicv)
  _Blade_7->_parameters.r0 = 1.0;
  #define r0 (_Blade_7->_parameters.r0)
  _Blade_7->_parameters.t0 = 1.0;
  #define t0 (_Blade_7->_parameters.t0)
  _Blade_7->_parameters.Q = 1.8734;
  #define Q (_Blade_7->_parameters.Q)
  _Blade_7->_parameters.RV = instrument->_parameters._mono_r_v;
  #define RV (_Blade_7->_parameters.RV)
  _Blade_7->_parameters.RH = instrument->_parameters._mono_r_h;
  #define RH (_Blade_7->_parameters.RH)
  _Blade_7->_parameters.DM = mono_d;
  #define DM (_Blade_7->_parameters.DM)
  _Blade_7->_parameters.mosaic = 0;
  #define mosaic (_Blade_7->_parameters.mosaic)
  _Blade_7->_parameters.width = 0;
  #define width (_Blade_7->_parameters.width)
  _Blade_7->_parameters.height = 0;
  #define height (_Blade_7->_parameters.height)
  _Blade_7->_parameters.verbose = 0;
  #define verbose (_Blade_7->_parameters.verbose)
  _Blade_7->_parameters.order = 0;
  #define order (_Blade_7->_parameters.order)

  #define mos_rms_y (_Blade_7->_parameters.mos_rms_y)
  #define mos_rms_z (_Blade_7->_parameters.mos_rms_z)
  #define mos_rms_max (_Blade_7->_parameters.mos_rms_max)
  #define mono_Q (_Blade_7->_parameters.mono_Q)
  #define SlabWidth (_Blade_7->_parameters.SlabWidth)
  #define SlabHeight (_Blade_7->_parameters.SlabHeight)
  #define rTable (_Blade_7->_parameters.rTable)
  #define tTable (_Blade_7->_parameters.tTable)
  #define row (_Blade_7->_parameters.row)
  #define col (_Blade_7->_parameters.col)
  #define tiltH (_Blade_7->_parameters.tiltH)
  #define tiltV (_Blade_7->_parameters.tiltV)

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
  /* component Blade_7=Monochromator_curved() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Blade_6->_rotation_absolute, _Blade_7->_rotation_absolute);
    rot_transpose(_Blade_6->_rotation_absolute, tr1);
    rot_mul(_Blade_7->_rotation_absolute, tr1, _Blade_7->_rotation_relative);
    _Blade_7->_rotation_is_identity =  rot_test_identity(_Blade_7->_rotation_relative);
    tc1 = coords_set(
      wafer_d, 0, 0);
    rot_transpose(_Blade_6->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Blade_7->_position_absolute = coords_add(_Blade_6->_position_absolute, tc2);
    tc1 = coords_sub(_Blade_6->_position_absolute, _Blade_7->_position_absolute);
    _Blade_7->_position_relative = rot_apply(_Blade_7->_rotation_absolute, tc1);
  } /* Blade_7=Monochromator_curved() AT ROTATED */
  DEBUG_COMPONENT("Blade_7", _Blade_7->_position_absolute, _Blade_7->_rotation_absolute);
  instrument->_position_absolute[22] = _Blade_7->_position_absolute;
  instrument->_position_relative[22] = _Blade_7->_position_relative;
  instrument->counter_N[22]  = instrument->counter_P[22] = instrument->counter_P2[22] = 0;
  instrument->counter_AbsorbProp[22]= 0;
  return(0);
} /* _Blade_7_setpos */

/* component Blade_8=Monochromator_curved() SETTING, POSITION/ROTATION */
int _Blade_8_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Blade_8_setpos] component Blade_8=Monochromator_curved() SETTING [/usr/share/mcstas/3.0-dev/optics/Monochromator_curved.comp:148]");
  stracpy(_Blade_8->_name, "Blade_8", 16384);
  stracpy(_Blade_8->_type, "Monochromator_curved", 16384);
  _Blade_8->_index=23;
  if("NULL" && strlen("NULL"))
    stracpy(_Blade_8->_parameters.reflect, "NULL" ? "NULL" : "", 16384);
  else 
  _Blade_8->_parameters.reflect[0]='\0';
  #define reflect (_Blade_8->_parameters.reflect)
  if("NULL" && strlen("NULL"))
    stracpy(_Blade_8->_parameters.transmit, "NULL" ? "NULL" : "", 16384);
  else 
  _Blade_8->_parameters.transmit[0]='\0';
  #define transmit (_Blade_8->_parameters.transmit)
  _Blade_8->_parameters.zwidth = 0.22 / 51.0;
  #define zwidth (_Blade_8->_parameters.zwidth)
  _Blade_8->_parameters.yheight = 0.1395 / 9.0;
  #define yheight (_Blade_8->_parameters.yheight)
  _Blade_8->_parameters.gap = 0.0;
  #define gap (_Blade_8->_parameters.gap)
  _Blade_8->_parameters.NH = 51;
  #define NH (_Blade_8->_parameters.NH)
  _Blade_8->_parameters.NV = 9;
  #define NV (_Blade_8->_parameters.NV)
  _Blade_8->_parameters.mosaich = instrument->_parameters._mono_mosh;
  #define mosaich (_Blade_8->_parameters.mosaich)
  _Blade_8->_parameters.mosaicv = instrument->_parameters._mono_mosv;
  #define mosaicv (_Blade_8->_parameters.mosaicv)
  _Blade_8->_parameters.r0 = 1.0;
  #define r0 (_Blade_8->_parameters.r0)
  _Blade_8->_parameters.t0 = 1.0;
  #define t0 (_Blade_8->_parameters.t0)
  _Blade_8->_parameters.Q = 1.8734;
  #define Q (_Blade_8->_parameters.Q)
  _Blade_8->_parameters.RV = instrument->_parameters._mono_r_v;
  #define RV (_Blade_8->_parameters.RV)
  _Blade_8->_parameters.RH = instrument->_parameters._mono_r_h;
  #define RH (_Blade_8->_parameters.RH)
  _Blade_8->_parameters.DM = mono_d;
  #define DM (_Blade_8->_parameters.DM)
  _Blade_8->_parameters.mosaic = 0;
  #define mosaic (_Blade_8->_parameters.mosaic)
  _Blade_8->_parameters.width = 0;
  #define width (_Blade_8->_parameters.width)
  _Blade_8->_parameters.height = 0;
  #define height (_Blade_8->_parameters.height)
  _Blade_8->_parameters.verbose = 0;
  #define verbose (_Blade_8->_parameters.verbose)
  _Blade_8->_parameters.order = 0;
  #define order (_Blade_8->_parameters.order)

  #define mos_rms_y (_Blade_8->_parameters.mos_rms_y)
  #define mos_rms_z (_Blade_8->_parameters.mos_rms_z)
  #define mos_rms_max (_Blade_8->_parameters.mos_rms_max)
  #define mono_Q (_Blade_8->_parameters.mono_Q)
  #define SlabWidth (_Blade_8->_parameters.SlabWidth)
  #define SlabHeight (_Blade_8->_parameters.SlabHeight)
  #define rTable (_Blade_8->_parameters.rTable)
  #define tTable (_Blade_8->_parameters.tTable)
  #define row (_Blade_8->_parameters.row)
  #define col (_Blade_8->_parameters.col)
  #define tiltH (_Blade_8->_parameters.tiltH)
  #define tiltV (_Blade_8->_parameters.tiltV)

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
  /* component Blade_8=Monochromator_curved() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Blade_7->_rotation_absolute, _Blade_8->_rotation_absolute);
    rot_transpose(_Blade_7->_rotation_absolute, tr1);
    rot_mul(_Blade_8->_rotation_absolute, tr1, _Blade_8->_rotation_relative);
    _Blade_8->_rotation_is_identity =  rot_test_identity(_Blade_8->_rotation_relative);
    tc1 = coords_set(
      wafer_d, 0, 0);
    rot_transpose(_Blade_7->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Blade_8->_position_absolute = coords_add(_Blade_7->_position_absolute, tc2);
    tc1 = coords_sub(_Blade_7->_position_absolute, _Blade_8->_position_absolute);
    _Blade_8->_position_relative = rot_apply(_Blade_8->_rotation_absolute, tc1);
  } /* Blade_8=Monochromator_curved() AT ROTATED */
  DEBUG_COMPONENT("Blade_8", _Blade_8->_position_absolute, _Blade_8->_rotation_absolute);
  instrument->_position_absolute[23] = _Blade_8->_position_absolute;
  instrument->_position_relative[23] = _Blade_8->_position_relative;
  instrument->counter_N[23]  = instrument->counter_P[23] = instrument->counter_P2[23] = 0;
  instrument->counter_AbsorbProp[23]= 0;
  return(0);
} /* _Blade_8_setpos */

/* component Blade_9=Monochromator_curved() SETTING, POSITION/ROTATION */
int _Blade_9_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Blade_9_setpos] component Blade_9=Monochromator_curved() SETTING [/usr/share/mcstas/3.0-dev/optics/Monochromator_curved.comp:148]");
  stracpy(_Blade_9->_name, "Blade_9", 16384);
  stracpy(_Blade_9->_type, "Monochromator_curved", 16384);
  _Blade_9->_index=24;
  if("NULL" && strlen("NULL"))
    stracpy(_Blade_9->_parameters.reflect, "NULL" ? "NULL" : "", 16384);
  else 
  _Blade_9->_parameters.reflect[0]='\0';
  #define reflect (_Blade_9->_parameters.reflect)
  if("NULL" && strlen("NULL"))
    stracpy(_Blade_9->_parameters.transmit, "NULL" ? "NULL" : "", 16384);
  else 
  _Blade_9->_parameters.transmit[0]='\0';
  #define transmit (_Blade_9->_parameters.transmit)
  _Blade_9->_parameters.zwidth = 0.22 / 51.0;
  #define zwidth (_Blade_9->_parameters.zwidth)
  _Blade_9->_parameters.yheight = 0.1395 / 9.0;
  #define yheight (_Blade_9->_parameters.yheight)
  _Blade_9->_parameters.gap = 0.0;
  #define gap (_Blade_9->_parameters.gap)
  _Blade_9->_parameters.NH = 51;
  #define NH (_Blade_9->_parameters.NH)
  _Blade_9->_parameters.NV = 9;
  #define NV (_Blade_9->_parameters.NV)
  _Blade_9->_parameters.mosaich = instrument->_parameters._mono_mosh;
  #define mosaich (_Blade_9->_parameters.mosaich)
  _Blade_9->_parameters.mosaicv = instrument->_parameters._mono_mosv;
  #define mosaicv (_Blade_9->_parameters.mosaicv)
  _Blade_9->_parameters.r0 = 1.0;
  #define r0 (_Blade_9->_parameters.r0)
  _Blade_9->_parameters.t0 = 1.0;
  #define t0 (_Blade_9->_parameters.t0)
  _Blade_9->_parameters.Q = 1.8734;
  #define Q (_Blade_9->_parameters.Q)
  _Blade_9->_parameters.RV = instrument->_parameters._mono_r_v;
  #define RV (_Blade_9->_parameters.RV)
  _Blade_9->_parameters.RH = instrument->_parameters._mono_r_h;
  #define RH (_Blade_9->_parameters.RH)
  _Blade_9->_parameters.DM = mono_d;
  #define DM (_Blade_9->_parameters.DM)
  _Blade_9->_parameters.mosaic = 0;
  #define mosaic (_Blade_9->_parameters.mosaic)
  _Blade_9->_parameters.width = 0;
  #define width (_Blade_9->_parameters.width)
  _Blade_9->_parameters.height = 0;
  #define height (_Blade_9->_parameters.height)
  _Blade_9->_parameters.verbose = 0;
  #define verbose (_Blade_9->_parameters.verbose)
  _Blade_9->_parameters.order = 0;
  #define order (_Blade_9->_parameters.order)

  #define mos_rms_y (_Blade_9->_parameters.mos_rms_y)
  #define mos_rms_z (_Blade_9->_parameters.mos_rms_z)
  #define mos_rms_max (_Blade_9->_parameters.mos_rms_max)
  #define mono_Q (_Blade_9->_parameters.mono_Q)
  #define SlabWidth (_Blade_9->_parameters.SlabWidth)
  #define SlabHeight (_Blade_9->_parameters.SlabHeight)
  #define rTable (_Blade_9->_parameters.rTable)
  #define tTable (_Blade_9->_parameters.tTable)
  #define row (_Blade_9->_parameters.row)
  #define col (_Blade_9->_parameters.col)
  #define tiltH (_Blade_9->_parameters.tiltH)
  #define tiltV (_Blade_9->_parameters.tiltV)

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
  /* component Blade_9=Monochromator_curved() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Blade_8->_rotation_absolute, _Blade_9->_rotation_absolute);
    rot_transpose(_Blade_8->_rotation_absolute, tr1);
    rot_mul(_Blade_9->_rotation_absolute, tr1, _Blade_9->_rotation_relative);
    _Blade_9->_rotation_is_identity =  rot_test_identity(_Blade_9->_rotation_relative);
    tc1 = coords_set(
      wafer_d, 0, 0);
    rot_transpose(_Blade_8->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Blade_9->_position_absolute = coords_add(_Blade_8->_position_absolute, tc2);
    tc1 = coords_sub(_Blade_8->_position_absolute, _Blade_9->_position_absolute);
    _Blade_9->_position_relative = rot_apply(_Blade_9->_rotation_absolute, tc1);
  } /* Blade_9=Monochromator_curved() AT ROTATED */
  DEBUG_COMPONENT("Blade_9", _Blade_9->_position_absolute, _Blade_9->_rotation_absolute);
  instrument->_position_absolute[24] = _Blade_9->_position_absolute;
  instrument->_position_relative[24] = _Blade_9->_position_relative;
  instrument->counter_N[24]  = instrument->counter_P[24] = instrument->counter_P2[24] = 0;
  instrument->counter_AbsorbProp[24]= 0;
  return(0);
} /* _Blade_9_setpos */

/* component Blade_10=Monochromator_curved() SETTING, POSITION/ROTATION */
int _Blade_10_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Blade_10_setpos] component Blade_10=Monochromator_curved() SETTING [/usr/share/mcstas/3.0-dev/optics/Monochromator_curved.comp:148]");
  stracpy(_Blade_10->_name, "Blade_10", 16384);
  stracpy(_Blade_10->_type, "Monochromator_curved", 16384);
  _Blade_10->_index=25;
  if("NULL" && strlen("NULL"))
    stracpy(_Blade_10->_parameters.reflect, "NULL" ? "NULL" : "", 16384);
  else 
  _Blade_10->_parameters.reflect[0]='\0';
  #define reflect (_Blade_10->_parameters.reflect)
  if("NULL" && strlen("NULL"))
    stracpy(_Blade_10->_parameters.transmit, "NULL" ? "NULL" : "", 16384);
  else 
  _Blade_10->_parameters.transmit[0]='\0';
  #define transmit (_Blade_10->_parameters.transmit)
  _Blade_10->_parameters.zwidth = 0.22 / 51.0;
  #define zwidth (_Blade_10->_parameters.zwidth)
  _Blade_10->_parameters.yheight = 0.1395 / 9.0;
  #define yheight (_Blade_10->_parameters.yheight)
  _Blade_10->_parameters.gap = 0.0;
  #define gap (_Blade_10->_parameters.gap)
  _Blade_10->_parameters.NH = 51;
  #define NH (_Blade_10->_parameters.NH)
  _Blade_10->_parameters.NV = 9;
  #define NV (_Blade_10->_parameters.NV)
  _Blade_10->_parameters.mosaich = instrument->_parameters._mono_mosh;
  #define mosaich (_Blade_10->_parameters.mosaich)
  _Blade_10->_parameters.mosaicv = instrument->_parameters._mono_mosv;
  #define mosaicv (_Blade_10->_parameters.mosaicv)
  _Blade_10->_parameters.r0 = 1.0;
  #define r0 (_Blade_10->_parameters.r0)
  _Blade_10->_parameters.t0 = 1.0;
  #define t0 (_Blade_10->_parameters.t0)
  _Blade_10->_parameters.Q = 1.8734;
  #define Q (_Blade_10->_parameters.Q)
  _Blade_10->_parameters.RV = instrument->_parameters._mono_r_v;
  #define RV (_Blade_10->_parameters.RV)
  _Blade_10->_parameters.RH = instrument->_parameters._mono_r_h;
  #define RH (_Blade_10->_parameters.RH)
  _Blade_10->_parameters.DM = mono_d;
  #define DM (_Blade_10->_parameters.DM)
  _Blade_10->_parameters.mosaic = 0;
  #define mosaic (_Blade_10->_parameters.mosaic)
  _Blade_10->_parameters.width = 0;
  #define width (_Blade_10->_parameters.width)
  _Blade_10->_parameters.height = 0;
  #define height (_Blade_10->_parameters.height)
  _Blade_10->_parameters.verbose = 0;
  #define verbose (_Blade_10->_parameters.verbose)
  _Blade_10->_parameters.order = 0;
  #define order (_Blade_10->_parameters.order)

  #define mos_rms_y (_Blade_10->_parameters.mos_rms_y)
  #define mos_rms_z (_Blade_10->_parameters.mos_rms_z)
  #define mos_rms_max (_Blade_10->_parameters.mos_rms_max)
  #define mono_Q (_Blade_10->_parameters.mono_Q)
  #define SlabWidth (_Blade_10->_parameters.SlabWidth)
  #define SlabHeight (_Blade_10->_parameters.SlabHeight)
  #define rTable (_Blade_10->_parameters.rTable)
  #define tTable (_Blade_10->_parameters.tTable)
  #define row (_Blade_10->_parameters.row)
  #define col (_Blade_10->_parameters.col)
  #define tiltH (_Blade_10->_parameters.tiltH)
  #define tiltV (_Blade_10->_parameters.tiltV)

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
  /* component Blade_10=Monochromator_curved() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Blade_9->_rotation_absolute, _Blade_10->_rotation_absolute);
    rot_transpose(_Blade_9->_rotation_absolute, tr1);
    rot_mul(_Blade_10->_rotation_absolute, tr1, _Blade_10->_rotation_relative);
    _Blade_10->_rotation_is_identity =  rot_test_identity(_Blade_10->_rotation_relative);
    tc1 = coords_set(
      wafer_d, 0, 0);
    rot_transpose(_Blade_9->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Blade_10->_position_absolute = coords_add(_Blade_9->_position_absolute, tc2);
    tc1 = coords_sub(_Blade_9->_position_absolute, _Blade_10->_position_absolute);
    _Blade_10->_position_relative = rot_apply(_Blade_10->_rotation_absolute, tc1);
  } /* Blade_10=Monochromator_curved() AT ROTATED */
  DEBUG_COMPONENT("Blade_10", _Blade_10->_position_absolute, _Blade_10->_rotation_absolute);
  instrument->_position_absolute[25] = _Blade_10->_position_absolute;
  instrument->_position_relative[25] = _Blade_10->_position_relative;
  instrument->counter_N[25]  = instrument->counter_P[25] = instrument->counter_P2[25] = 0;
  instrument->counter_AbsorbProp[25]= 0;
  return(0);
} /* _Blade_10_setpos */

/* component Blade_11=Monochromator_curved() SETTING, POSITION/ROTATION */
int _Blade_11_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Blade_11_setpos] component Blade_11=Monochromator_curved() SETTING [/usr/share/mcstas/3.0-dev/optics/Monochromator_curved.comp:148]");
  stracpy(_Blade_11->_name, "Blade_11", 16384);
  stracpy(_Blade_11->_type, "Monochromator_curved", 16384);
  _Blade_11->_index=26;
  if("NULL" && strlen("NULL"))
    stracpy(_Blade_11->_parameters.reflect, "NULL" ? "NULL" : "", 16384);
  else 
  _Blade_11->_parameters.reflect[0]='\0';
  #define reflect (_Blade_11->_parameters.reflect)
  if("NULL" && strlen("NULL"))
    stracpy(_Blade_11->_parameters.transmit, "NULL" ? "NULL" : "", 16384);
  else 
  _Blade_11->_parameters.transmit[0]='\0';
  #define transmit (_Blade_11->_parameters.transmit)
  _Blade_11->_parameters.zwidth = 0.22 / 51.0;
  #define zwidth (_Blade_11->_parameters.zwidth)
  _Blade_11->_parameters.yheight = 0.1395 / 9.0;
  #define yheight (_Blade_11->_parameters.yheight)
  _Blade_11->_parameters.gap = 0.0;
  #define gap (_Blade_11->_parameters.gap)
  _Blade_11->_parameters.NH = 51;
  #define NH (_Blade_11->_parameters.NH)
  _Blade_11->_parameters.NV = 9;
  #define NV (_Blade_11->_parameters.NV)
  _Blade_11->_parameters.mosaich = instrument->_parameters._mono_mosh;
  #define mosaich (_Blade_11->_parameters.mosaich)
  _Blade_11->_parameters.mosaicv = instrument->_parameters._mono_mosv;
  #define mosaicv (_Blade_11->_parameters.mosaicv)
  _Blade_11->_parameters.r0 = 1.0;
  #define r0 (_Blade_11->_parameters.r0)
  _Blade_11->_parameters.t0 = 1.0;
  #define t0 (_Blade_11->_parameters.t0)
  _Blade_11->_parameters.Q = 1.8734;
  #define Q (_Blade_11->_parameters.Q)
  _Blade_11->_parameters.RV = instrument->_parameters._mono_r_v;
  #define RV (_Blade_11->_parameters.RV)
  _Blade_11->_parameters.RH = instrument->_parameters._mono_r_h;
  #define RH (_Blade_11->_parameters.RH)
  _Blade_11->_parameters.DM = mono_d;
  #define DM (_Blade_11->_parameters.DM)
  _Blade_11->_parameters.mosaic = 0;
  #define mosaic (_Blade_11->_parameters.mosaic)
  _Blade_11->_parameters.width = 0;
  #define width (_Blade_11->_parameters.width)
  _Blade_11->_parameters.height = 0;
  #define height (_Blade_11->_parameters.height)
  _Blade_11->_parameters.verbose = 0;
  #define verbose (_Blade_11->_parameters.verbose)
  _Blade_11->_parameters.order = 0;
  #define order (_Blade_11->_parameters.order)

  #define mos_rms_y (_Blade_11->_parameters.mos_rms_y)
  #define mos_rms_z (_Blade_11->_parameters.mos_rms_z)
  #define mos_rms_max (_Blade_11->_parameters.mos_rms_max)
  #define mono_Q (_Blade_11->_parameters.mono_Q)
  #define SlabWidth (_Blade_11->_parameters.SlabWidth)
  #define SlabHeight (_Blade_11->_parameters.SlabHeight)
  #define rTable (_Blade_11->_parameters.rTable)
  #define tTable (_Blade_11->_parameters.tTable)
  #define row (_Blade_11->_parameters.row)
  #define col (_Blade_11->_parameters.col)
  #define tiltH (_Blade_11->_parameters.tiltH)
  #define tiltV (_Blade_11->_parameters.tiltV)

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
  /* component Blade_11=Monochromator_curved() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Blade_10->_rotation_absolute, _Blade_11->_rotation_absolute);
    rot_transpose(_Blade_10->_rotation_absolute, tr1);
    rot_mul(_Blade_11->_rotation_absolute, tr1, _Blade_11->_rotation_relative);
    _Blade_11->_rotation_is_identity =  rot_test_identity(_Blade_11->_rotation_relative);
    tc1 = coords_set(
      wafer_d, 0, 0);
    rot_transpose(_Blade_10->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Blade_11->_position_absolute = coords_add(_Blade_10->_position_absolute, tc2);
    tc1 = coords_sub(_Blade_10->_position_absolute, _Blade_11->_position_absolute);
    _Blade_11->_position_relative = rot_apply(_Blade_11->_rotation_absolute, tc1);
  } /* Blade_11=Monochromator_curved() AT ROTATED */
  DEBUG_COMPONENT("Blade_11", _Blade_11->_position_absolute, _Blade_11->_rotation_absolute);
  instrument->_position_absolute[26] = _Blade_11->_position_absolute;
  instrument->_position_relative[26] = _Blade_11->_position_relative;
  instrument->counter_N[26]  = instrument->counter_P[26] = instrument->counter_P2[26] = 0;
  instrument->counter_AbsorbProp[26]= 0;
  return(0);
} /* _Blade_11_setpos */

/* component Blade_12=Monochromator_curved() SETTING, POSITION/ROTATION */
int _Blade_12_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Blade_12_setpos] component Blade_12=Monochromator_curved() SETTING [/usr/share/mcstas/3.0-dev/optics/Monochromator_curved.comp:148]");
  stracpy(_Blade_12->_name, "Blade_12", 16384);
  stracpy(_Blade_12->_type, "Monochromator_curved", 16384);
  _Blade_12->_index=27;
  if("NULL" && strlen("NULL"))
    stracpy(_Blade_12->_parameters.reflect, "NULL" ? "NULL" : "", 16384);
  else 
  _Blade_12->_parameters.reflect[0]='\0';
  #define reflect (_Blade_12->_parameters.reflect)
  if("NULL" && strlen("NULL"))
    stracpy(_Blade_12->_parameters.transmit, "NULL" ? "NULL" : "", 16384);
  else 
  _Blade_12->_parameters.transmit[0]='\0';
  #define transmit (_Blade_12->_parameters.transmit)
  _Blade_12->_parameters.zwidth = 0.22 / 51.0;
  #define zwidth (_Blade_12->_parameters.zwidth)
  _Blade_12->_parameters.yheight = 0.1395 / 9.0;
  #define yheight (_Blade_12->_parameters.yheight)
  _Blade_12->_parameters.gap = 0.0;
  #define gap (_Blade_12->_parameters.gap)
  _Blade_12->_parameters.NH = 51;
  #define NH (_Blade_12->_parameters.NH)
  _Blade_12->_parameters.NV = 9;
  #define NV (_Blade_12->_parameters.NV)
  _Blade_12->_parameters.mosaich = instrument->_parameters._mono_mosh;
  #define mosaich (_Blade_12->_parameters.mosaich)
  _Blade_12->_parameters.mosaicv = instrument->_parameters._mono_mosv;
  #define mosaicv (_Blade_12->_parameters.mosaicv)
  _Blade_12->_parameters.r0 = 1.0;
  #define r0 (_Blade_12->_parameters.r0)
  _Blade_12->_parameters.t0 = 1.0;
  #define t0 (_Blade_12->_parameters.t0)
  _Blade_12->_parameters.Q = 1.8734;
  #define Q (_Blade_12->_parameters.Q)
  _Blade_12->_parameters.RV = instrument->_parameters._mono_r_v;
  #define RV (_Blade_12->_parameters.RV)
  _Blade_12->_parameters.RH = instrument->_parameters._mono_r_h;
  #define RH (_Blade_12->_parameters.RH)
  _Blade_12->_parameters.DM = mono_d;
  #define DM (_Blade_12->_parameters.DM)
  _Blade_12->_parameters.mosaic = 0;
  #define mosaic (_Blade_12->_parameters.mosaic)
  _Blade_12->_parameters.width = 0;
  #define width (_Blade_12->_parameters.width)
  _Blade_12->_parameters.height = 0;
  #define height (_Blade_12->_parameters.height)
  _Blade_12->_parameters.verbose = 0;
  #define verbose (_Blade_12->_parameters.verbose)
  _Blade_12->_parameters.order = 0;
  #define order (_Blade_12->_parameters.order)

  #define mos_rms_y (_Blade_12->_parameters.mos_rms_y)
  #define mos_rms_z (_Blade_12->_parameters.mos_rms_z)
  #define mos_rms_max (_Blade_12->_parameters.mos_rms_max)
  #define mono_Q (_Blade_12->_parameters.mono_Q)
  #define SlabWidth (_Blade_12->_parameters.SlabWidth)
  #define SlabHeight (_Blade_12->_parameters.SlabHeight)
  #define rTable (_Blade_12->_parameters.rTable)
  #define tTable (_Blade_12->_parameters.tTable)
  #define row (_Blade_12->_parameters.row)
  #define col (_Blade_12->_parameters.col)
  #define tiltH (_Blade_12->_parameters.tiltH)
  #define tiltV (_Blade_12->_parameters.tiltV)

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
  /* component Blade_12=Monochromator_curved() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Blade_11->_rotation_absolute, _Blade_12->_rotation_absolute);
    rot_transpose(_Blade_11->_rotation_absolute, tr1);
    rot_mul(_Blade_12->_rotation_absolute, tr1, _Blade_12->_rotation_relative);
    _Blade_12->_rotation_is_identity =  rot_test_identity(_Blade_12->_rotation_relative);
    tc1 = coords_set(
      wafer_d, 0, 0);
    rot_transpose(_Blade_11->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Blade_12->_position_absolute = coords_add(_Blade_11->_position_absolute, tc2);
    tc1 = coords_sub(_Blade_11->_position_absolute, _Blade_12->_position_absolute);
    _Blade_12->_position_relative = rot_apply(_Blade_12->_rotation_absolute, tc1);
  } /* Blade_12=Monochromator_curved() AT ROTATED */
  DEBUG_COMPONENT("Blade_12", _Blade_12->_position_absolute, _Blade_12->_rotation_absolute);
  instrument->_position_absolute[27] = _Blade_12->_position_absolute;
  instrument->_position_relative[27] = _Blade_12->_position_relative;
  instrument->counter_N[27]  = instrument->counter_P[27] = instrument->counter_P2[27] = 0;
  instrument->counter_AbsorbProp[27]= 0;
  return(0);
} /* _Blade_12_setpos */

/* component Blade_13=Monochromator_curved() SETTING, POSITION/ROTATION */
int _Blade_13_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Blade_13_setpos] component Blade_13=Monochromator_curved() SETTING [/usr/share/mcstas/3.0-dev/optics/Monochromator_curved.comp:148]");
  stracpy(_Blade_13->_name, "Blade_13", 16384);
  stracpy(_Blade_13->_type, "Monochromator_curved", 16384);
  _Blade_13->_index=28;
  if("NULL" && strlen("NULL"))
    stracpy(_Blade_13->_parameters.reflect, "NULL" ? "NULL" : "", 16384);
  else 
  _Blade_13->_parameters.reflect[0]='\0';
  #define reflect (_Blade_13->_parameters.reflect)
  if("NULL" && strlen("NULL"))
    stracpy(_Blade_13->_parameters.transmit, "NULL" ? "NULL" : "", 16384);
  else 
  _Blade_13->_parameters.transmit[0]='\0';
  #define transmit (_Blade_13->_parameters.transmit)
  _Blade_13->_parameters.zwidth = 0.22 / 51.0;
  #define zwidth (_Blade_13->_parameters.zwidth)
  _Blade_13->_parameters.yheight = 0.1395 / 9.0;
  #define yheight (_Blade_13->_parameters.yheight)
  _Blade_13->_parameters.gap = 0.0;
  #define gap (_Blade_13->_parameters.gap)
  _Blade_13->_parameters.NH = 51;
  #define NH (_Blade_13->_parameters.NH)
  _Blade_13->_parameters.NV = 9;
  #define NV (_Blade_13->_parameters.NV)
  _Blade_13->_parameters.mosaich = instrument->_parameters._mono_mosh;
  #define mosaich (_Blade_13->_parameters.mosaich)
  _Blade_13->_parameters.mosaicv = instrument->_parameters._mono_mosv;
  #define mosaicv (_Blade_13->_parameters.mosaicv)
  _Blade_13->_parameters.r0 = 1.0;
  #define r0 (_Blade_13->_parameters.r0)
  _Blade_13->_parameters.t0 = 1.0;
  #define t0 (_Blade_13->_parameters.t0)
  _Blade_13->_parameters.Q = 1.8734;
  #define Q (_Blade_13->_parameters.Q)
  _Blade_13->_parameters.RV = instrument->_parameters._mono_r_v;
  #define RV (_Blade_13->_parameters.RV)
  _Blade_13->_parameters.RH = instrument->_parameters._mono_r_h;
  #define RH (_Blade_13->_parameters.RH)
  _Blade_13->_parameters.DM = mono_d;
  #define DM (_Blade_13->_parameters.DM)
  _Blade_13->_parameters.mosaic = 0;
  #define mosaic (_Blade_13->_parameters.mosaic)
  _Blade_13->_parameters.width = 0;
  #define width (_Blade_13->_parameters.width)
  _Blade_13->_parameters.height = 0;
  #define height (_Blade_13->_parameters.height)
  _Blade_13->_parameters.verbose = 0;
  #define verbose (_Blade_13->_parameters.verbose)
  _Blade_13->_parameters.order = 0;
  #define order (_Blade_13->_parameters.order)

  #define mos_rms_y (_Blade_13->_parameters.mos_rms_y)
  #define mos_rms_z (_Blade_13->_parameters.mos_rms_z)
  #define mos_rms_max (_Blade_13->_parameters.mos_rms_max)
  #define mono_Q (_Blade_13->_parameters.mono_Q)
  #define SlabWidth (_Blade_13->_parameters.SlabWidth)
  #define SlabHeight (_Blade_13->_parameters.SlabHeight)
  #define rTable (_Blade_13->_parameters.rTable)
  #define tTable (_Blade_13->_parameters.tTable)
  #define row (_Blade_13->_parameters.row)
  #define col (_Blade_13->_parameters.col)
  #define tiltH (_Blade_13->_parameters.tiltH)
  #define tiltV (_Blade_13->_parameters.tiltV)

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
  /* component Blade_13=Monochromator_curved() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Blade_12->_rotation_absolute, _Blade_13->_rotation_absolute);
    rot_transpose(_Blade_12->_rotation_absolute, tr1);
    rot_mul(_Blade_13->_rotation_absolute, tr1, _Blade_13->_rotation_relative);
    _Blade_13->_rotation_is_identity =  rot_test_identity(_Blade_13->_rotation_relative);
    tc1 = coords_set(
      wafer_d, 0, 0);
    rot_transpose(_Blade_12->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Blade_13->_position_absolute = coords_add(_Blade_12->_position_absolute, tc2);
    tc1 = coords_sub(_Blade_12->_position_absolute, _Blade_13->_position_absolute);
    _Blade_13->_position_relative = rot_apply(_Blade_13->_rotation_absolute, tc1);
  } /* Blade_13=Monochromator_curved() AT ROTATED */
  DEBUG_COMPONENT("Blade_13", _Blade_13->_position_absolute, _Blade_13->_rotation_absolute);
  instrument->_position_absolute[28] = _Blade_13->_position_absolute;
  instrument->_position_relative[28] = _Blade_13->_position_relative;
  instrument->counter_N[28]  = instrument->counter_P[28] = instrument->counter_P2[28] = 0;
  instrument->counter_AbsorbProp[28]= 0;
  return(0);
} /* _Blade_13_setpos */

/* component PSD_At_sec_shutter=PSD_monitor() SETTING, POSITION/ROTATION */
int _PSD_At_sec_shutter_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_PSD_At_sec_shutter_setpos] component PSD_At_sec_shutter=PSD_monitor() SETTING [/usr/share/mcstas/3.0-dev/monitors/PSD_monitor.comp:64]");
  stracpy(_PSD_At_sec_shutter->_name, "PSD_At_sec_shutter", 16384);
  stracpy(_PSD_At_sec_shutter->_type, "PSD_monitor", 16384);
  _PSD_At_sec_shutter->_index=29;
  _PSD_At_sec_shutter->_parameters.nx = 100;
  #define nx (_PSD_At_sec_shutter->_parameters.nx)
  _PSD_At_sec_shutter->_parameters.ny = 100;
  #define ny (_PSD_At_sec_shutter->_parameters.ny)
  if("PSD_At_sec_shutter.out" && strlen("PSD_At_sec_shutter.out"))
    stracpy(_PSD_At_sec_shutter->_parameters.filename, "PSD_At_sec_shutter.out" ? "PSD_At_sec_shutter.out" : "", 16384);
  else 
  _PSD_At_sec_shutter->_parameters.filename[0]='\0';
  #define filename (_PSD_At_sec_shutter->_parameters.filename)
  _PSD_At_sec_shutter->_parameters.xmin = -0.1;
  #define xmin (_PSD_At_sec_shutter->_parameters.xmin)
  _PSD_At_sec_shutter->_parameters.xmax = 0.1;
  #define xmax (_PSD_At_sec_shutter->_parameters.xmax)
  _PSD_At_sec_shutter->_parameters.ymin = -0.1;
  #define ymin (_PSD_At_sec_shutter->_parameters.ymin)
  _PSD_At_sec_shutter->_parameters.ymax = 0.1;
  #define ymax (_PSD_At_sec_shutter->_parameters.ymax)
  _PSD_At_sec_shutter->_parameters.xwidth = 0;
  #define xwidth (_PSD_At_sec_shutter->_parameters.xwidth)
  _PSD_At_sec_shutter->_parameters.yheight = 0;
  #define yheight (_PSD_At_sec_shutter->_parameters.yheight)
  _PSD_At_sec_shutter->_parameters.restore_neutron = 0;
  #define restore_neutron (_PSD_At_sec_shutter->_parameters.restore_neutron)

  #define PSD_N (_PSD_At_sec_shutter->_parameters.PSD_N)
  #define PSD_p (_PSD_At_sec_shutter->_parameters.PSD_p)
  #define PSD_p2 (_PSD_At_sec_shutter->_parameters.PSD_p2)

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
  /* component PSD_At_sec_shutter=PSD_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Prim_axes->_rotation_absolute, _PSD_At_sec_shutter->_rotation_absolute);
    rot_transpose(_Blade_13->_rotation_absolute, tr1);
    rot_mul(_PSD_At_sec_shutter->_rotation_absolute, tr1, _PSD_At_sec_shutter->_rotation_relative);
    _PSD_At_sec_shutter->_rotation_is_identity =  rot_test_identity(_PSD_At_sec_shutter->_rotation_relative);
    tc1 = coords_set(
      0, 0, chamber_col_start -0.1736);
    rot_transpose(_Prim_axes->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _PSD_At_sec_shutter->_position_absolute = coords_add(_Prim_axes->_position_absolute, tc2);
    tc1 = coords_sub(_Blade_13->_position_absolute, _PSD_At_sec_shutter->_position_absolute);
    _PSD_At_sec_shutter->_position_relative = rot_apply(_PSD_At_sec_shutter->_rotation_absolute, tc1);
  } /* PSD_At_sec_shutter=PSD_monitor() AT ROTATED */
  DEBUG_COMPONENT("PSD_At_sec_shutter", _PSD_At_sec_shutter->_position_absolute, _PSD_At_sec_shutter->_rotation_absolute);
  instrument->_position_absolute[29] = _PSD_At_sec_shutter->_position_absolute;
  instrument->_position_relative[29] = _PSD_At_sec_shutter->_position_relative;
  instrument->counter_N[29]  = instrument->counter_P[29] = instrument->counter_P2[29] = 0;
  instrument->counter_AbsorbProp[29]= 0;
  return(0);
} /* _PSD_At_sec_shutter_setpos */

/* component LAM_At_sec_shutter=L_monitor() SETTING, POSITION/ROTATION */
int _LAM_At_sec_shutter_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_LAM_At_sec_shutter_setpos] component LAM_At_sec_shutter=L_monitor() SETTING [/usr/share/mcstas/3.0-dev/monitors/L_monitor.comp:65]");
  stracpy(_LAM_At_sec_shutter->_name, "LAM_At_sec_shutter", 16384);
  stracpy(_LAM_At_sec_shutter->_type, "L_monitor", 16384);
  _LAM_At_sec_shutter->_index=30;
  _LAM_At_sec_shutter->_parameters.nL = 100;
  #define nL (_LAM_At_sec_shutter->_parameters.nL)
  if("LAM_At_sec_shutter.out" && strlen("LAM_At_sec_shutter.out"))
    stracpy(_LAM_At_sec_shutter->_parameters.filename, "LAM_At_sec_shutter.out" ? "LAM_At_sec_shutter.out" : "", 16384);
  else 
  _LAM_At_sec_shutter->_parameters.filename[0]='\0';
  #define filename (_LAM_At_sec_shutter->_parameters.filename)
  _LAM_At_sec_shutter->_parameters.xmin = -0.1;
  #define xmin (_LAM_At_sec_shutter->_parameters.xmin)
  _LAM_At_sec_shutter->_parameters.xmax = 0.1;
  #define xmax (_LAM_At_sec_shutter->_parameters.xmax)
  _LAM_At_sec_shutter->_parameters.ymin = -0.1;
  #define ymin (_LAM_At_sec_shutter->_parameters.ymin)
  _LAM_At_sec_shutter->_parameters.ymax = 0.1;
  #define ymax (_LAM_At_sec_shutter->_parameters.ymax)
  _LAM_At_sec_shutter->_parameters.xwidth = 0;
  #define xwidth (_LAM_At_sec_shutter->_parameters.xwidth)
  _LAM_At_sec_shutter->_parameters.yheight = 0;
  #define yheight (_LAM_At_sec_shutter->_parameters.yheight)
  _LAM_At_sec_shutter->_parameters.Lmin = instrument->_parameters._source_lam_min;
  #define Lmin (_LAM_At_sec_shutter->_parameters.Lmin)
  _LAM_At_sec_shutter->_parameters.Lmax = instrument->_parameters._source_lam_max;
  #define Lmax (_LAM_At_sec_shutter->_parameters.Lmax)
  _LAM_At_sec_shutter->_parameters.restore_neutron = 0;
  #define restore_neutron (_LAM_At_sec_shutter->_parameters.restore_neutron)

  #define L_N (_LAM_At_sec_shutter->_parameters.L_N)
  #define L_p (_LAM_At_sec_shutter->_parameters.L_p)
  #define L_p2 (_LAM_At_sec_shutter->_parameters.L_p2)

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
  /* component LAM_At_sec_shutter=L_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _PSD_At_sec_shutter->_rotation_absolute, _LAM_At_sec_shutter->_rotation_absolute);
    rot_transpose(_PSD_At_sec_shutter->_rotation_absolute, tr1);
    rot_mul(_LAM_At_sec_shutter->_rotation_absolute, tr1, _LAM_At_sec_shutter->_rotation_relative);
    _LAM_At_sec_shutter->_rotation_is_identity =  rot_test_identity(_LAM_At_sec_shutter->_rotation_relative);
    tc1 = coords_set(
      0, 0, 0);
    rot_transpose(_PSD_At_sec_shutter->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _LAM_At_sec_shutter->_position_absolute = coords_add(_PSD_At_sec_shutter->_position_absolute, tc2);
    tc1 = coords_sub(_PSD_At_sec_shutter->_position_absolute, _LAM_At_sec_shutter->_position_absolute);
    _LAM_At_sec_shutter->_position_relative = rot_apply(_LAM_At_sec_shutter->_rotation_absolute, tc1);
  } /* LAM_At_sec_shutter=L_monitor() AT ROTATED */
  DEBUG_COMPONENT("LAM_At_sec_shutter", _LAM_At_sec_shutter->_position_absolute, _LAM_At_sec_shutter->_rotation_absolute);
  instrument->_position_absolute[30] = _LAM_At_sec_shutter->_position_absolute;
  instrument->_position_relative[30] = _LAM_At_sec_shutter->_position_relative;
  instrument->counter_N[30]  = instrument->counter_P[30] = instrument->counter_P2[30] = 0;
  instrument->counter_AbsorbProp[30]= 0;
  return(0);
} /* _LAM_At_sec_shutter_setpos */

/* component Inside_chamber_collimator=Slit() SETTING, POSITION/ROTATION */
int _Inside_chamber_collimator_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Inside_chamber_collimator_setpos] component Inside_chamber_collimator=Slit() SETTING [/usr/share/mcstas/3.0-dev/optics/Slit.comp:47]");
  stracpy(_Inside_chamber_collimator->_name, "Inside_chamber_collimator", 16384);
  stracpy(_Inside_chamber_collimator->_type, "Slit", 16384);
  _Inside_chamber_collimator->_index=31;
  _Inside_chamber_collimator->_parameters.xmin = -0.021;
  #define xmin (_Inside_chamber_collimator->_parameters.xmin)
  _Inside_chamber_collimator->_parameters.xmax = 0.021;
  #define xmax (_Inside_chamber_collimator->_parameters.xmax)
  _Inside_chamber_collimator->_parameters.ymin = -0.0425;
  #define ymin (_Inside_chamber_collimator->_parameters.ymin)
  _Inside_chamber_collimator->_parameters.ymax = 0.0425;
  #define ymax (_Inside_chamber_collimator->_parameters.ymax)
  _Inside_chamber_collimator->_parameters.radius = 0;
  #define radius (_Inside_chamber_collimator->_parameters.radius)
  _Inside_chamber_collimator->_parameters.xwidth = 0;
  #define xwidth (_Inside_chamber_collimator->_parameters.xwidth)
  _Inside_chamber_collimator->_parameters.yheight = 0;
  #define yheight (_Inside_chamber_collimator->_parameters.yheight)

  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef radius
  #undef xwidth
  #undef yheight
  /* component Inside_chamber_collimator=Slit() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Prim_axes->_rotation_absolute, _Inside_chamber_collimator->_rotation_absolute);
    rot_transpose(_LAM_At_sec_shutter->_rotation_absolute, tr1);
    rot_mul(_Inside_chamber_collimator->_rotation_absolute, tr1, _Inside_chamber_collimator->_rotation_relative);
    _Inside_chamber_collimator->_rotation_is_identity =  rot_test_identity(_Inside_chamber_collimator->_rotation_relative);
    tc1 = coords_set(
      0, 0, chamber_col_start);
    rot_transpose(_Prim_axes->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Inside_chamber_collimator->_position_absolute = coords_add(_Prim_axes->_position_absolute, tc2);
    tc1 = coords_sub(_LAM_At_sec_shutter->_position_absolute, _Inside_chamber_collimator->_position_absolute);
    _Inside_chamber_collimator->_position_relative = rot_apply(_Inside_chamber_collimator->_rotation_absolute, tc1);
  } /* Inside_chamber_collimator=Slit() AT ROTATED */
  DEBUG_COMPONENT("Inside_chamber_collimator", _Inside_chamber_collimator->_position_absolute, _Inside_chamber_collimator->_rotation_absolute);
  instrument->_position_absolute[31] = _Inside_chamber_collimator->_position_absolute;
  instrument->_position_relative[31] = _Inside_chamber_collimator->_position_relative;
  instrument->counter_N[31]  = instrument->counter_P[31] = instrument->counter_P2[31] = 0;
  instrument->counter_AbsorbProp[31]= 0;
  return(0);
} /* _Inside_chamber_collimator_setpos */

/* component PSD_After_inside_chamber_collimator=PSD_monitor() SETTING, POSITION/ROTATION */
int _PSD_After_inside_chamber_collimator_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_PSD_After_inside_chamber_collimator_setpos] component PSD_After_inside_chamber_collimator=PSD_monitor() SETTING [/usr/share/mcstas/3.0-dev/monitors/PSD_monitor.comp:64]");
  stracpy(_PSD_After_inside_chamber_collimator->_name, "PSD_After_inside_chamber_collimator", 16384);
  stracpy(_PSD_After_inside_chamber_collimator->_type, "PSD_monitor", 16384);
  _PSD_After_inside_chamber_collimator->_index=32;
  _PSD_After_inside_chamber_collimator->_parameters.nx = 100;
  #define nx (_PSD_After_inside_chamber_collimator->_parameters.nx)
  _PSD_After_inside_chamber_collimator->_parameters.ny = 100;
  #define ny (_PSD_After_inside_chamber_collimator->_parameters.ny)
  if("PSD_After_inside_chamber_collimator.out" && strlen("PSD_After_inside_chamber_collimator.out"))
    stracpy(_PSD_After_inside_chamber_collimator->_parameters.filename, "PSD_After_inside_chamber_collimator.out" ? "PSD_After_inside_chamber_collimator.out" : "", 16384);
  else 
  _PSD_After_inside_chamber_collimator->_parameters.filename[0]='\0';
  #define filename (_PSD_After_inside_chamber_collimator->_parameters.filename)
  _PSD_After_inside_chamber_collimator->_parameters.xmin = -0.1;
  #define xmin (_PSD_After_inside_chamber_collimator->_parameters.xmin)
  _PSD_After_inside_chamber_collimator->_parameters.xmax = 0.1;
  #define xmax (_PSD_After_inside_chamber_collimator->_parameters.xmax)
  _PSD_After_inside_chamber_collimator->_parameters.ymin = -0.1;
  #define ymin (_PSD_After_inside_chamber_collimator->_parameters.ymin)
  _PSD_After_inside_chamber_collimator->_parameters.ymax = 0.1;
  #define ymax (_PSD_After_inside_chamber_collimator->_parameters.ymax)
  _PSD_After_inside_chamber_collimator->_parameters.xwidth = 0;
  #define xwidth (_PSD_After_inside_chamber_collimator->_parameters.xwidth)
  _PSD_After_inside_chamber_collimator->_parameters.yheight = 0;
  #define yheight (_PSD_After_inside_chamber_collimator->_parameters.yheight)
  _PSD_After_inside_chamber_collimator->_parameters.restore_neutron = 0;
  #define restore_neutron (_PSD_After_inside_chamber_collimator->_parameters.restore_neutron)

  #define PSD_N (_PSD_After_inside_chamber_collimator->_parameters.PSD_N)
  #define PSD_p (_PSD_After_inside_chamber_collimator->_parameters.PSD_p)
  #define PSD_p2 (_PSD_After_inside_chamber_collimator->_parameters.PSD_p2)

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
  /* component PSD_After_inside_chamber_collimator=PSD_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Inside_chamber_collimator->_rotation_absolute, _PSD_After_inside_chamber_collimator->_rotation_absolute);
    rot_transpose(_Inside_chamber_collimator->_rotation_absolute, tr1);
    rot_mul(_PSD_After_inside_chamber_collimator->_rotation_absolute, tr1, _PSD_After_inside_chamber_collimator->_rotation_relative);
    _PSD_After_inside_chamber_collimator->_rotation_is_identity =  rot_test_identity(_PSD_After_inside_chamber_collimator->_rotation_relative);
    tc1 = coords_set(
      0, 0, 0);
    rot_transpose(_Inside_chamber_collimator->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _PSD_After_inside_chamber_collimator->_position_absolute = coords_add(_Inside_chamber_collimator->_position_absolute, tc2);
    tc1 = coords_sub(_Inside_chamber_collimator->_position_absolute, _PSD_After_inside_chamber_collimator->_position_absolute);
    _PSD_After_inside_chamber_collimator->_position_relative = rot_apply(_PSD_After_inside_chamber_collimator->_rotation_absolute, tc1);
  } /* PSD_After_inside_chamber_collimator=PSD_monitor() AT ROTATED */
  DEBUG_COMPONENT("PSD_After_inside_chamber_collimator", _PSD_After_inside_chamber_collimator->_position_absolute, _PSD_After_inside_chamber_collimator->_rotation_absolute);
  instrument->_position_absolute[32] = _PSD_After_inside_chamber_collimator->_position_absolute;
  instrument->_position_relative[32] = _PSD_After_inside_chamber_collimator->_position_relative;
  instrument->counter_N[32]  = instrument->counter_P[32] = instrument->counter_P2[32] = 0;
  instrument->counter_AbsorbProp[32]= 0;
  return(0);
} /* _PSD_After_inside_chamber_collimator_setpos */

/* component Outside_chamber_collimator=Slit() SETTING, POSITION/ROTATION */
int _Outside_chamber_collimator_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Outside_chamber_collimator_setpos] component Outside_chamber_collimator=Slit() SETTING [/usr/share/mcstas/3.0-dev/optics/Slit.comp:47]");
  stracpy(_Outside_chamber_collimator->_name, "Outside_chamber_collimator", 16384);
  stracpy(_Outside_chamber_collimator->_type, "Slit", 16384);
  _Outside_chamber_collimator->_index=33;
  _Outside_chamber_collimator->_parameters.xmin = -1 * outside_chamber_collimator_w / 2.0;
  #define xmin (_Outside_chamber_collimator->_parameters.xmin)
  _Outside_chamber_collimator->_parameters.xmax = outside_chamber_collimator_w / 2.0;
  #define xmax (_Outside_chamber_collimator->_parameters.xmax)
  _Outside_chamber_collimator->_parameters.ymin = -1 * outside_chamber_collimator_h / 2.0;
  #define ymin (_Outside_chamber_collimator->_parameters.ymin)
  _Outside_chamber_collimator->_parameters.ymax = outside_chamber_collimator_h / 2.0;
  #define ymax (_Outside_chamber_collimator->_parameters.ymax)
  _Outside_chamber_collimator->_parameters.radius = 0;
  #define radius (_Outside_chamber_collimator->_parameters.radius)
  _Outside_chamber_collimator->_parameters.xwidth = 0;
  #define xwidth (_Outside_chamber_collimator->_parameters.xwidth)
  _Outside_chamber_collimator->_parameters.yheight = 0;
  #define yheight (_Outside_chamber_collimator->_parameters.yheight)

  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef radius
  #undef xwidth
  #undef yheight
  /* component Outside_chamber_collimator=Slit() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Inside_chamber_collimator->_rotation_absolute, _Outside_chamber_collimator->_rotation_absolute);
    rot_transpose(_PSD_After_inside_chamber_collimator->_rotation_absolute, tr1);
    rot_mul(_Outside_chamber_collimator->_rotation_absolute, tr1, _Outside_chamber_collimator->_rotation_relative);
    _Outside_chamber_collimator->_rotation_is_identity =  rot_test_identity(_Outside_chamber_collimator->_rotation_relative);
    tc1 = coords_set(
      0, 0, chamber_col_length);
    rot_transpose(_Inside_chamber_collimator->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Outside_chamber_collimator->_position_absolute = coords_add(_Inside_chamber_collimator->_position_absolute, tc2);
    tc1 = coords_sub(_PSD_After_inside_chamber_collimator->_position_absolute, _Outside_chamber_collimator->_position_absolute);
    _Outside_chamber_collimator->_position_relative = rot_apply(_Outside_chamber_collimator->_rotation_absolute, tc1);
  } /* Outside_chamber_collimator=Slit() AT ROTATED */
  DEBUG_COMPONENT("Outside_chamber_collimator", _Outside_chamber_collimator->_position_absolute, _Outside_chamber_collimator->_rotation_absolute);
  instrument->_position_absolute[33] = _Outside_chamber_collimator->_position_absolute;
  instrument->_position_relative[33] = _Outside_chamber_collimator->_position_relative;
  instrument->counter_N[33]  = instrument->counter_P[33] = instrument->counter_P2[33] = 0;
  instrument->counter_AbsorbProp[33]= 0;
  return(0);
} /* _Outside_chamber_collimator_setpos */

/* component PSD_Outside_chamber_collimator_1=PSD_monitor() SETTING, POSITION/ROTATION */
int _PSD_Outside_chamber_collimator_1_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_PSD_Outside_chamber_collimator_1_setpos] component PSD_Outside_chamber_collimator_1=PSD_monitor() SETTING [/usr/share/mcstas/3.0-dev/monitors/PSD_monitor.comp:64]");
  stracpy(_PSD_Outside_chamber_collimator_1->_name, "PSD_Outside_chamber_collimator_1", 16384);
  stracpy(_PSD_Outside_chamber_collimator_1->_type, "PSD_monitor", 16384);
  _PSD_Outside_chamber_collimator_1->_index=34;
  _PSD_Outside_chamber_collimator_1->_parameters.nx = 100;
  #define nx (_PSD_Outside_chamber_collimator_1->_parameters.nx)
  _PSD_Outside_chamber_collimator_1->_parameters.ny = 100;
  #define ny (_PSD_Outside_chamber_collimator_1->_parameters.ny)
  if("PSD_Outside_chamber_collimator_1.out" && strlen("PSD_Outside_chamber_collimator_1.out"))
    stracpy(_PSD_Outside_chamber_collimator_1->_parameters.filename, "PSD_Outside_chamber_collimator_1.out" ? "PSD_Outside_chamber_collimator_1.out" : "", 16384);
  else 
  _PSD_Outside_chamber_collimator_1->_parameters.filename[0]='\0';
  #define filename (_PSD_Outside_chamber_collimator_1->_parameters.filename)
  _PSD_Outside_chamber_collimator_1->_parameters.xmin = -0.1;
  #define xmin (_PSD_Outside_chamber_collimator_1->_parameters.xmin)
  _PSD_Outside_chamber_collimator_1->_parameters.xmax = 0.1;
  #define xmax (_PSD_Outside_chamber_collimator_1->_parameters.xmax)
  _PSD_Outside_chamber_collimator_1->_parameters.ymin = -0.1;
  #define ymin (_PSD_Outside_chamber_collimator_1->_parameters.ymin)
  _PSD_Outside_chamber_collimator_1->_parameters.ymax = 0.1;
  #define ymax (_PSD_Outside_chamber_collimator_1->_parameters.ymax)
  _PSD_Outside_chamber_collimator_1->_parameters.xwidth = 0;
  #define xwidth (_PSD_Outside_chamber_collimator_1->_parameters.xwidth)
  _PSD_Outside_chamber_collimator_1->_parameters.yheight = 0;
  #define yheight (_PSD_Outside_chamber_collimator_1->_parameters.yheight)
  _PSD_Outside_chamber_collimator_1->_parameters.restore_neutron = 0;
  #define restore_neutron (_PSD_Outside_chamber_collimator_1->_parameters.restore_neutron)

  #define PSD_N (_PSD_Outside_chamber_collimator_1->_parameters.PSD_N)
  #define PSD_p (_PSD_Outside_chamber_collimator_1->_parameters.PSD_p)
  #define PSD_p2 (_PSD_Outside_chamber_collimator_1->_parameters.PSD_p2)

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
  /* component PSD_Outside_chamber_collimator_1=PSD_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Outside_chamber_collimator->_rotation_absolute, _PSD_Outside_chamber_collimator_1->_rotation_absolute);
    rot_transpose(_Outside_chamber_collimator->_rotation_absolute, tr1);
    rot_mul(_PSD_Outside_chamber_collimator_1->_rotation_absolute, tr1, _PSD_Outside_chamber_collimator_1->_rotation_relative);
    _PSD_Outside_chamber_collimator_1->_rotation_is_identity =  rot_test_identity(_PSD_Outside_chamber_collimator_1->_rotation_relative);
    tc1 = coords_set(
      0, 0, 0);
    rot_transpose(_Outside_chamber_collimator->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _PSD_Outside_chamber_collimator_1->_position_absolute = coords_add(_Outside_chamber_collimator->_position_absolute, tc2);
    tc1 = coords_sub(_Outside_chamber_collimator->_position_absolute, _PSD_Outside_chamber_collimator_1->_position_absolute);
    _PSD_Outside_chamber_collimator_1->_position_relative = rot_apply(_PSD_Outside_chamber_collimator_1->_rotation_absolute, tc1);
  } /* PSD_Outside_chamber_collimator_1=PSD_monitor() AT ROTATED */
  DEBUG_COMPONENT("PSD_Outside_chamber_collimator_1", _PSD_Outside_chamber_collimator_1->_position_absolute, _PSD_Outside_chamber_collimator_1->_rotation_absolute);
  instrument->_position_absolute[34] = _PSD_Outside_chamber_collimator_1->_position_absolute;
  instrument->_position_relative[34] = _PSD_Outside_chamber_collimator_1->_position_relative;
  instrument->counter_N[34]  = instrument->counter_P[34] = instrument->counter_P2[34] = 0;
  instrument->counter_AbsorbProp[34]= 0;
  return(0);
} /* _PSD_Outside_chamber_collimator_1_setpos */

/* component PSD_Outside_chamber_collimator=PSD_monitor() SETTING, POSITION/ROTATION */
int _PSD_Outside_chamber_collimator_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_PSD_Outside_chamber_collimator_setpos] component PSD_Outside_chamber_collimator=PSD_monitor() SETTING [/usr/share/mcstas/3.0-dev/monitors/PSD_monitor.comp:64]");
  stracpy(_PSD_Outside_chamber_collimator->_name, "PSD_Outside_chamber_collimator", 16384);
  stracpy(_PSD_Outside_chamber_collimator->_type, "PSD_monitor", 16384);
  _PSD_Outside_chamber_collimator->_index=35;
  _PSD_Outside_chamber_collimator->_parameters.nx = 100;
  #define nx (_PSD_Outside_chamber_collimator->_parameters.nx)
  _PSD_Outside_chamber_collimator->_parameters.ny = 100;
  #define ny (_PSD_Outside_chamber_collimator->_parameters.ny)
  if("PSD_Outside_chamber_collimator.out" && strlen("PSD_Outside_chamber_collimator.out"))
    stracpy(_PSD_Outside_chamber_collimator->_parameters.filename, "PSD_Outside_chamber_collimator.out" ? "PSD_Outside_chamber_collimator.out" : "", 16384);
  else 
  _PSD_Outside_chamber_collimator->_parameters.filename[0]='\0';
  #define filename (_PSD_Outside_chamber_collimator->_parameters.filename)
  _PSD_Outside_chamber_collimator->_parameters.xmin = -0.1;
  #define xmin (_PSD_Outside_chamber_collimator->_parameters.xmin)
  _PSD_Outside_chamber_collimator->_parameters.xmax = 0.1;
  #define xmax (_PSD_Outside_chamber_collimator->_parameters.xmax)
  _PSD_Outside_chamber_collimator->_parameters.ymin = -0.1;
  #define ymin (_PSD_Outside_chamber_collimator->_parameters.ymin)
  _PSD_Outside_chamber_collimator->_parameters.ymax = 0.1;
  #define ymax (_PSD_Outside_chamber_collimator->_parameters.ymax)
  _PSD_Outside_chamber_collimator->_parameters.xwidth = 0;
  #define xwidth (_PSD_Outside_chamber_collimator->_parameters.xwidth)
  _PSD_Outside_chamber_collimator->_parameters.yheight = 0;
  #define yheight (_PSD_Outside_chamber_collimator->_parameters.yheight)
  _PSD_Outside_chamber_collimator->_parameters.restore_neutron = 0;
  #define restore_neutron (_PSD_Outside_chamber_collimator->_parameters.restore_neutron)

  #define PSD_N (_PSD_Outside_chamber_collimator->_parameters.PSD_N)
  #define PSD_p (_PSD_Outside_chamber_collimator->_parameters.PSD_p)
  #define PSD_p2 (_PSD_Outside_chamber_collimator->_parameters.PSD_p2)

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
  /* component PSD_Outside_chamber_collimator=PSD_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Prim_axes->_rotation_absolute, _PSD_Outside_chamber_collimator->_rotation_absolute);
    rot_transpose(_PSD_Outside_chamber_collimator_1->_rotation_absolute, tr1);
    rot_mul(_PSD_Outside_chamber_collimator->_rotation_absolute, tr1, _PSD_Outside_chamber_collimator->_rotation_relative);
    _PSD_Outside_chamber_collimator->_rotation_is_identity =  rot_test_identity(_PSD_Outside_chamber_collimator->_rotation_relative);
    tc1 = coords_set(
      0, 0, chamber_col_start + chamber_col_length);
    rot_transpose(_Prim_axes->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _PSD_Outside_chamber_collimator->_position_absolute = coords_add(_Prim_axes->_position_absolute, tc2);
    tc1 = coords_sub(_PSD_Outside_chamber_collimator_1->_position_absolute, _PSD_Outside_chamber_collimator->_position_absolute);
    _PSD_Outside_chamber_collimator->_position_relative = rot_apply(_PSD_Outside_chamber_collimator->_rotation_absolute, tc1);
  } /* PSD_Outside_chamber_collimator=PSD_monitor() AT ROTATED */
  DEBUG_COMPONENT("PSD_Outside_chamber_collimator", _PSD_Outside_chamber_collimator->_position_absolute, _PSD_Outside_chamber_collimator->_rotation_absolute);
  instrument->_position_absolute[35] = _PSD_Outside_chamber_collimator->_position_absolute;
  instrument->_position_relative[35] = _PSD_Outside_chamber_collimator->_position_relative;
  instrument->counter_N[35]  = instrument->counter_P[35] = instrument->counter_P2[35] = 0;
  instrument->counter_AbsorbProp[35]= 0;
  return(0);
} /* _PSD_Outside_chamber_collimator_setpos */

/* component Incident_slit_h=Slit() SETTING, POSITION/ROTATION */
int _Incident_slit_h_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Incident_slit_h_setpos] component Incident_slit_h=Slit() SETTING [/usr/share/mcstas/3.0-dev/optics/Slit.comp:47]");
  stracpy(_Incident_slit_h->_name, "Incident_slit_h", 16384);
  stracpy(_Incident_slit_h->_type, "Slit", 16384);
  _Incident_slit_h->_index=36;
  _Incident_slit_h->_parameters.xmin = inc_slit_xmin_h;
  #define xmin (_Incident_slit_h->_parameters.xmin)
  _Incident_slit_h->_parameters.xmax = inc_slit_xmax_h;
  #define xmax (_Incident_slit_h->_parameters.xmax)
  _Incident_slit_h->_parameters.ymin = inc_slit_ymin;
  #define ymin (_Incident_slit_h->_parameters.ymin)
  _Incident_slit_h->_parameters.ymax = inc_slit_ymax;
  #define ymax (_Incident_slit_h->_parameters.ymax)
  _Incident_slit_h->_parameters.radius = 0;
  #define radius (_Incident_slit_h->_parameters.radius)
  _Incident_slit_h->_parameters.xwidth = 0;
  #define xwidth (_Incident_slit_h->_parameters.xwidth)
  _Incident_slit_h->_parameters.yheight = 0;
  #define yheight (_Incident_slit_h->_parameters.yheight)

  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef radius
  #undef xwidth
  #undef yheight
  /* component Incident_slit_h=Slit() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Prim_axes->_rotation_absolute, _Incident_slit_h->_rotation_absolute);
    rot_transpose(_PSD_Outside_chamber_collimator->_rotation_absolute, tr1);
    rot_mul(_Incident_slit_h->_rotation_absolute, tr1, _Incident_slit_h->_rotation_relative);
    _Incident_slit_h->_rotation_is_identity =  rot_test_identity(_Incident_slit_h->_rotation_relative);
    tc1 = coords_set(
      instrument->_parameters._inc_slit_dx, 0, instrument->_parameters._mono_to_cor - instrument->_parameters._inc_slit_to_cor - instrument->_parameters._inc_slit_sep);
    rot_transpose(_Prim_axes->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Incident_slit_h->_position_absolute = coords_add(_Prim_axes->_position_absolute, tc2);
    tc1 = coords_sub(_PSD_Outside_chamber_collimator->_position_absolute, _Incident_slit_h->_position_absolute);
    _Incident_slit_h->_position_relative = rot_apply(_Incident_slit_h->_rotation_absolute, tc1);
  } /* Incident_slit_h=Slit() AT ROTATED */
  DEBUG_COMPONENT("Incident_slit_h", _Incident_slit_h->_position_absolute, _Incident_slit_h->_rotation_absolute);
  instrument->_position_absolute[36] = _Incident_slit_h->_position_absolute;
  instrument->_position_relative[36] = _Incident_slit_h->_position_relative;
  instrument->counter_N[36]  = instrument->counter_P[36] = instrument->counter_P2[36] = 0;
  instrument->counter_AbsorbProp[36]= 0;
  return(0);
} /* _Incident_slit_h_setpos */

/* component Incident_slit_w=Slit() SETTING, POSITION/ROTATION */
int _Incident_slit_w_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Incident_slit_w_setpos] component Incident_slit_w=Slit() SETTING [/usr/share/mcstas/3.0-dev/optics/Slit.comp:47]");
  stracpy(_Incident_slit_w->_name, "Incident_slit_w", 16384);
  stracpy(_Incident_slit_w->_type, "Slit", 16384);
  _Incident_slit_w->_index=37;
  _Incident_slit_w->_parameters.xmin = inc_slit_xmin;
  #define xmin (_Incident_slit_w->_parameters.xmin)
  _Incident_slit_w->_parameters.xmax = inc_slit_xmax;
  #define xmax (_Incident_slit_w->_parameters.xmax)
  _Incident_slit_w->_parameters.ymin = inc_slit_ymin_w;
  #define ymin (_Incident_slit_w->_parameters.ymin)
  _Incident_slit_w->_parameters.ymax = inc_slit_ymax_w;
  #define ymax (_Incident_slit_w->_parameters.ymax)
  _Incident_slit_w->_parameters.radius = 0;
  #define radius (_Incident_slit_w->_parameters.radius)
  _Incident_slit_w->_parameters.xwidth = 0;
  #define xwidth (_Incident_slit_w->_parameters.xwidth)
  _Incident_slit_w->_parameters.yheight = 0;
  #define yheight (_Incident_slit_w->_parameters.yheight)

  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef radius
  #undef xwidth
  #undef yheight
  /* component Incident_slit_w=Slit() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Prim_axes->_rotation_absolute, _Incident_slit_w->_rotation_absolute);
    rot_transpose(_Incident_slit_h->_rotation_absolute, tr1);
    rot_mul(_Incident_slit_w->_rotation_absolute, tr1, _Incident_slit_w->_rotation_relative);
    _Incident_slit_w->_rotation_is_identity =  rot_test_identity(_Incident_slit_w->_rotation_relative);
    tc1 = coords_set(
      instrument->_parameters._inc_slit_dx, 0, instrument->_parameters._mono_to_cor - instrument->_parameters._inc_slit_to_cor);
    rot_transpose(_Prim_axes->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Incident_slit_w->_position_absolute = coords_add(_Prim_axes->_position_absolute, tc2);
    tc1 = coords_sub(_Incident_slit_h->_position_absolute, _Incident_slit_w->_position_absolute);
    _Incident_slit_w->_position_relative = rot_apply(_Incident_slit_w->_rotation_absolute, tc1);
  } /* Incident_slit_w=Slit() AT ROTATED */
  DEBUG_COMPONENT("Incident_slit_w", _Incident_slit_w->_position_absolute, _Incident_slit_w->_rotation_absolute);
  instrument->_position_absolute[37] = _Incident_slit_w->_position_absolute;
  instrument->_position_relative[37] = _Incident_slit_w->_position_relative;
  instrument->counter_N[37]  = instrument->counter_P[37] = instrument->counter_P2[37] = 0;
  instrument->counter_AbsorbProp[37]= 0;
  return(0);
} /* _Incident_slit_w_setpos */

/* component PSD_After_Incident_slit_w=PSD_monitor() SETTING, POSITION/ROTATION */
int _PSD_After_Incident_slit_w_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_PSD_After_Incident_slit_w_setpos] component PSD_After_Incident_slit_w=PSD_monitor() SETTING [/usr/share/mcstas/3.0-dev/monitors/PSD_monitor.comp:64]");
  stracpy(_PSD_After_Incident_slit_w->_name, "PSD_After_Incident_slit_w", 16384);
  stracpy(_PSD_After_Incident_slit_w->_type, "PSD_monitor", 16384);
  _PSD_After_Incident_slit_w->_index=38;
  _PSD_After_Incident_slit_w->_parameters.nx = 100;
  #define nx (_PSD_After_Incident_slit_w->_parameters.nx)
  _PSD_After_Incident_slit_w->_parameters.ny = 100;
  #define ny (_PSD_After_Incident_slit_w->_parameters.ny)
  if("PSD_After_Incident_slit_w.out" && strlen("PSD_After_Incident_slit_w.out"))
    stracpy(_PSD_After_Incident_slit_w->_parameters.filename, "PSD_After_Incident_slit_w.out" ? "PSD_After_Incident_slit_w.out" : "", 16384);
  else 
  _PSD_After_Incident_slit_w->_parameters.filename[0]='\0';
  #define filename (_PSD_After_Incident_slit_w->_parameters.filename)
  _PSD_After_Incident_slit_w->_parameters.xmin = -0.1;
  #define xmin (_PSD_After_Incident_slit_w->_parameters.xmin)
  _PSD_After_Incident_slit_w->_parameters.xmax = 0.1;
  #define xmax (_PSD_After_Incident_slit_w->_parameters.xmax)
  _PSD_After_Incident_slit_w->_parameters.ymin = -0.1;
  #define ymin (_PSD_After_Incident_slit_w->_parameters.ymin)
  _PSD_After_Incident_slit_w->_parameters.ymax = 0.1;
  #define ymax (_PSD_After_Incident_slit_w->_parameters.ymax)
  _PSD_After_Incident_slit_w->_parameters.xwidth = 0;
  #define xwidth (_PSD_After_Incident_slit_w->_parameters.xwidth)
  _PSD_After_Incident_slit_w->_parameters.yheight = 0;
  #define yheight (_PSD_After_Incident_slit_w->_parameters.yheight)
  _PSD_After_Incident_slit_w->_parameters.restore_neutron = 0;
  #define restore_neutron (_PSD_After_Incident_slit_w->_parameters.restore_neutron)

  #define PSD_N (_PSD_After_Incident_slit_w->_parameters.PSD_N)
  #define PSD_p (_PSD_After_Incident_slit_w->_parameters.PSD_p)
  #define PSD_p2 (_PSD_After_Incident_slit_w->_parameters.PSD_p2)

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
  /* component PSD_After_Incident_slit_w=PSD_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Incident_slit_w->_rotation_absolute, _PSD_After_Incident_slit_w->_rotation_absolute);
    rot_transpose(_Incident_slit_w->_rotation_absolute, tr1);
    rot_mul(_PSD_After_Incident_slit_w->_rotation_absolute, tr1, _PSD_After_Incident_slit_w->_rotation_relative);
    _PSD_After_Incident_slit_w->_rotation_is_identity =  rot_test_identity(_PSD_After_Incident_slit_w->_rotation_relative);
    tc1 = coords_set(
      0, 0, 0);
    rot_transpose(_Incident_slit_w->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _PSD_After_Incident_slit_w->_position_absolute = coords_add(_Incident_slit_w->_position_absolute, tc2);
    tc1 = coords_sub(_Incident_slit_w->_position_absolute, _PSD_After_Incident_slit_w->_position_absolute);
    _PSD_After_Incident_slit_w->_position_relative = rot_apply(_PSD_After_Incident_slit_w->_rotation_absolute, tc1);
  } /* PSD_After_Incident_slit_w=PSD_monitor() AT ROTATED */
  DEBUG_COMPONENT("PSD_After_Incident_slit_w", _PSD_After_Incident_slit_w->_position_absolute, _PSD_After_Incident_slit_w->_rotation_absolute);
  instrument->_position_absolute[38] = _PSD_After_Incident_slit_w->_position_absolute;
  instrument->_position_relative[38] = _PSD_After_Incident_slit_w->_position_relative;
  instrument->counter_N[38]  = instrument->counter_P[38] = instrument->counter_P2[38] = 0;
  instrument->counter_AbsorbProp[38]= 0;
  return(0);
} /* _PSD_After_Incident_slit_w_setpos */

/* component Center_of_rotation=Arm() SETTING, POSITION/ROTATION */
int _Center_of_rotation_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Center_of_rotation_setpos] component Center_of_rotation=Arm() SETTING [Arm:0]");
  stracpy(_Center_of_rotation->_name, "Center_of_rotation", 16384);
  stracpy(_Center_of_rotation->_type, "Arm", 16384);
  _Center_of_rotation->_index=39;
  /* component Center_of_rotation=Arm() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Prim_axes->_rotation_absolute, _Center_of_rotation->_rotation_absolute);
    rot_transpose(_PSD_After_Incident_slit_w->_rotation_absolute, tr1);
    rot_mul(_Center_of_rotation->_rotation_absolute, tr1, _Center_of_rotation->_rotation_relative);
    _Center_of_rotation->_rotation_is_identity =  rot_test_identity(_Center_of_rotation->_rotation_relative);
    tc1 = coords_set(
      0, 0, instrument->_parameters._mono_to_cor);
    rot_transpose(_Prim_axes->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Center_of_rotation->_position_absolute = coords_add(_Prim_axes->_position_absolute, tc2);
    tc1 = coords_sub(_PSD_After_Incident_slit_w->_position_absolute, _Center_of_rotation->_position_absolute);
    _Center_of_rotation->_position_relative = rot_apply(_Center_of_rotation->_rotation_absolute, tc1);
  } /* Center_of_rotation=Arm() AT ROTATED */
  DEBUG_COMPONENT("Center_of_rotation", _Center_of_rotation->_position_absolute, _Center_of_rotation->_rotation_absolute);
  instrument->_position_absolute[39] = _Center_of_rotation->_position_absolute;
  instrument->_position_relative[39] = _Center_of_rotation->_position_relative;
  instrument->counter_N[39]  = instrument->counter_P[39] = instrument->counter_P2[39] = 0;
  instrument->counter_AbsorbProp[39]= 0;
  return(0);
} /* _Center_of_rotation_setpos */

/* component Sample_rotation=Arm() SETTING, POSITION/ROTATION */
int _Sample_rotation_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Sample_rotation_setpos] component Sample_rotation=Arm() SETTING [Arm:0]");
  stracpy(_Sample_rotation->_name, "Sample_rotation", 16384);
  stracpy(_Sample_rotation->_type, "Arm", 16384);
  _Sample_rotation->_index=40;
  /* component Sample_rotation=Arm() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (instrument->_parameters._sample_dom)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _Center_of_rotation->_rotation_absolute, _Sample_rotation->_rotation_absolute);
    rot_transpose(_Center_of_rotation->_rotation_absolute, tr1);
    rot_mul(_Sample_rotation->_rotation_absolute, tr1, _Sample_rotation->_rotation_relative);
    _Sample_rotation->_rotation_is_identity =  rot_test_identity(_Sample_rotation->_rotation_relative);
    tc1 = coords_set(
      0, 0, 0);
    rot_transpose(_Center_of_rotation->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Sample_rotation->_position_absolute = coords_add(_Center_of_rotation->_position_absolute, tc2);
    tc1 = coords_sub(_Center_of_rotation->_position_absolute, _Sample_rotation->_position_absolute);
    _Sample_rotation->_position_relative = rot_apply(_Sample_rotation->_rotation_absolute, tc1);
  } /* Sample_rotation=Arm() AT ROTATED */
  DEBUG_COMPONENT("Sample_rotation", _Sample_rotation->_position_absolute, _Sample_rotation->_rotation_absolute);
  instrument->_position_absolute[40] = _Sample_rotation->_position_absolute;
  instrument->_position_relative[40] = _Sample_rotation->_position_relative;
  instrument->counter_N[40]  = instrument->counter_P[40] = instrument->counter_P2[40] = 0;
  instrument->counter_AbsorbProp[40]= 0;
  return(0);
} /* _Sample_rotation_setpos */

/* component Sample_location=Arm() SETTING, POSITION/ROTATION */
int _Sample_location_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Sample_location_setpos] component Sample_location=Arm() SETTING [Arm:0]");
  stracpy(_Sample_location->_name, "Sample_location", 16384);
  stracpy(_Sample_location->_type, "Arm", 16384);
  _Sample_location->_index=41;
  /* component Sample_location=Arm() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Sample_rotation->_rotation_absolute, _Sample_location->_rotation_absolute);
    rot_transpose(_Sample_rotation->_rotation_absolute, tr1);
    rot_mul(_Sample_location->_rotation_absolute, tr1, _Sample_location->_rotation_relative);
    _Sample_location->_rotation_is_identity =  rot_test_identity(_Sample_location->_rotation_relative);
    tc1 = coords_set(
      instrument->_parameters._sample_dx, instrument->_parameters._sample_dy, instrument->_parameters._sample_dz);
    rot_transpose(_Sample_rotation->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Sample_location->_position_absolute = coords_add(_Sample_rotation->_position_absolute, tc2);
    tc1 = coords_sub(_Sample_rotation->_position_absolute, _Sample_location->_position_absolute);
    _Sample_location->_position_relative = rot_apply(_Sample_location->_rotation_absolute, tc1);
  } /* Sample_location=Arm() AT ROTATED */
  DEBUG_COMPONENT("Sample_location", _Sample_location->_position_absolute, _Sample_location->_rotation_absolute);
  instrument->_position_absolute[41] = _Sample_location->_position_absolute;
  instrument->_position_relative[41] = _Sample_location->_position_relative;
  instrument->counter_N[41]  = instrument->counter_P[41] = instrument->counter_P2[41] = 0;
  instrument->counter_AbsorbProp[41]= 0;
  return(0);
} /* _Sample_location_setpos */

/* component PSD_Center_of_rotation=PSD_monitor() SETTING, POSITION/ROTATION */
int _PSD_Center_of_rotation_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_PSD_Center_of_rotation_setpos] component PSD_Center_of_rotation=PSD_monitor() SETTING [/usr/share/mcstas/3.0-dev/monitors/PSD_monitor.comp:64]");
  stracpy(_PSD_Center_of_rotation->_name, "PSD_Center_of_rotation", 16384);
  stracpy(_PSD_Center_of_rotation->_type, "PSD_monitor", 16384);
  _PSD_Center_of_rotation->_index=42;
  _PSD_Center_of_rotation->_parameters.nx = 100;
  #define nx (_PSD_Center_of_rotation->_parameters.nx)
  _PSD_Center_of_rotation->_parameters.ny = 100;
  #define ny (_PSD_Center_of_rotation->_parameters.ny)
  if("PSD_Center_of_rotation" && strlen("PSD_Center_of_rotation"))
    stracpy(_PSD_Center_of_rotation->_parameters.filename, "PSD_Center_of_rotation" ? "PSD_Center_of_rotation" : "", 16384);
  else 
  _PSD_Center_of_rotation->_parameters.filename[0]='\0';
  #define filename (_PSD_Center_of_rotation->_parameters.filename)
  _PSD_Center_of_rotation->_parameters.xmin = -0.035;
  #define xmin (_PSD_Center_of_rotation->_parameters.xmin)
  _PSD_Center_of_rotation->_parameters.xmax = 0.035;
  #define xmax (_PSD_Center_of_rotation->_parameters.xmax)
  _PSD_Center_of_rotation->_parameters.ymin = -0.035;
  #define ymin (_PSD_Center_of_rotation->_parameters.ymin)
  _PSD_Center_of_rotation->_parameters.ymax = 0.035;
  #define ymax (_PSD_Center_of_rotation->_parameters.ymax)
  _PSD_Center_of_rotation->_parameters.xwidth = 0;
  #define xwidth (_PSD_Center_of_rotation->_parameters.xwidth)
  _PSD_Center_of_rotation->_parameters.yheight = 0;
  #define yheight (_PSD_Center_of_rotation->_parameters.yheight)
  _PSD_Center_of_rotation->_parameters.restore_neutron = 0;
  #define restore_neutron (_PSD_Center_of_rotation->_parameters.restore_neutron)

  #define PSD_N (_PSD_Center_of_rotation->_parameters.PSD_N)
  #define PSD_p (_PSD_Center_of_rotation->_parameters.PSD_p)
  #define PSD_p2 (_PSD_Center_of_rotation->_parameters.PSD_p2)

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
  /* component PSD_Center_of_rotation=PSD_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Center_of_rotation->_rotation_absolute, _PSD_Center_of_rotation->_rotation_absolute);
    rot_transpose(_Sample_location->_rotation_absolute, tr1);
    rot_mul(_PSD_Center_of_rotation->_rotation_absolute, tr1, _PSD_Center_of_rotation->_rotation_relative);
    _PSD_Center_of_rotation->_rotation_is_identity =  rot_test_identity(_PSD_Center_of_rotation->_rotation_relative);
    tc1 = coords_set(
      0, 0, 0);
    rot_transpose(_Center_of_rotation->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _PSD_Center_of_rotation->_position_absolute = coords_add(_Center_of_rotation->_position_absolute, tc2);
    tc1 = coords_sub(_Sample_location->_position_absolute, _PSD_Center_of_rotation->_position_absolute);
    _PSD_Center_of_rotation->_position_relative = rot_apply(_PSD_Center_of_rotation->_rotation_absolute, tc1);
  } /* PSD_Center_of_rotation=PSD_monitor() AT ROTATED */
  DEBUG_COMPONENT("PSD_Center_of_rotation", _PSD_Center_of_rotation->_position_absolute, _PSD_Center_of_rotation->_rotation_absolute);
  instrument->_position_absolute[42] = _PSD_Center_of_rotation->_position_absolute;
  instrument->_position_relative[42] = _PSD_Center_of_rotation->_position_relative;
  instrument->counter_N[42]  = instrument->counter_P[42] = instrument->counter_P2[42] = 0;
  instrument->counter_AbsorbProp[42]= 0;
  return(0);
} /* _PSD_Center_of_rotation_setpos */

/* component DIV_Center_of_rotation=Divergence_monitor() SETTING, POSITION/ROTATION */
int _DIV_Center_of_rotation_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_DIV_Center_of_rotation_setpos] component DIV_Center_of_rotation=Divergence_monitor() SETTING [/usr/share/mcstas/3.0-dev/monitors/Divergence_monitor.comp:71]");
  stracpy(_DIV_Center_of_rotation->_name, "DIV_Center_of_rotation", 16384);
  stracpy(_DIV_Center_of_rotation->_type, "Divergence_monitor", 16384);
  _DIV_Center_of_rotation->_index=43;
  _DIV_Center_of_rotation->_parameters.nh = 100;
  #define nh (_DIV_Center_of_rotation->_parameters.nh)
  _DIV_Center_of_rotation->_parameters.nv = 100;
  #define nv (_DIV_Center_of_rotation->_parameters.nv)
  if("DIV_Center_of_rotation" && strlen("DIV_Center_of_rotation"))
    stracpy(_DIV_Center_of_rotation->_parameters.filename, "DIV_Center_of_rotation" ? "DIV_Center_of_rotation" : "", 16384);
  else 
  _DIV_Center_of_rotation->_parameters.filename[0]='\0';
  #define filename (_DIV_Center_of_rotation->_parameters.filename)
  _DIV_Center_of_rotation->_parameters.xmin = -0.035;
  #define xmin (_DIV_Center_of_rotation->_parameters.xmin)
  _DIV_Center_of_rotation->_parameters.xmax = 0.035;
  #define xmax (_DIV_Center_of_rotation->_parameters.xmax)
  _DIV_Center_of_rotation->_parameters.ymin = -0.035;
  #define ymin (_DIV_Center_of_rotation->_parameters.ymin)
  _DIV_Center_of_rotation->_parameters.ymax = 0.035;
  #define ymax (_DIV_Center_of_rotation->_parameters.ymax)
  _DIV_Center_of_rotation->_parameters.xwidth = 0;
  #define xwidth (_DIV_Center_of_rotation->_parameters.xwidth)
  _DIV_Center_of_rotation->_parameters.yheight = 0;
  #define yheight (_DIV_Center_of_rotation->_parameters.yheight)
  _DIV_Center_of_rotation->_parameters.maxdiv_h = 2;
  #define maxdiv_h (_DIV_Center_of_rotation->_parameters.maxdiv_h)
  _DIV_Center_of_rotation->_parameters.maxdiv_v = 2;
  #define maxdiv_v (_DIV_Center_of_rotation->_parameters.maxdiv_v)
  _DIV_Center_of_rotation->_parameters.restore_neutron = 0;
  #define restore_neutron (_DIV_Center_of_rotation->_parameters.restore_neutron)
  _DIV_Center_of_rotation->_parameters.nx = 0;
  #define nx (_DIV_Center_of_rotation->_parameters.nx)
  _DIV_Center_of_rotation->_parameters.ny = 0;
  #define ny (_DIV_Center_of_rotation->_parameters.ny)
  _DIV_Center_of_rotation->_parameters.nz = 1;
  #define nz (_DIV_Center_of_rotation->_parameters.nz)

  #define Div_N (_DIV_Center_of_rotation->_parameters.Div_N)
  #define Div_p (_DIV_Center_of_rotation->_parameters.Div_p)
  #define Div_p2 (_DIV_Center_of_rotation->_parameters.Div_p2)

  #undef nh
  #undef nv
  #undef filename
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef xwidth
  #undef yheight
  #undef maxdiv_h
  #undef maxdiv_v
  #undef restore_neutron
  #undef nx
  #undef ny
  #undef nz
  #undef Div_N
  #undef Div_p
  #undef Div_p2
  /* component DIV_Center_of_rotation=Divergence_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Center_of_rotation->_rotation_absolute, _DIV_Center_of_rotation->_rotation_absolute);
    rot_transpose(_PSD_Center_of_rotation->_rotation_absolute, tr1);
    rot_mul(_DIV_Center_of_rotation->_rotation_absolute, tr1, _DIV_Center_of_rotation->_rotation_relative);
    _DIV_Center_of_rotation->_rotation_is_identity =  rot_test_identity(_DIV_Center_of_rotation->_rotation_relative);
    tc1 = coords_set(
      0, 0, 0);
    rot_transpose(_Center_of_rotation->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _DIV_Center_of_rotation->_position_absolute = coords_add(_Center_of_rotation->_position_absolute, tc2);
    tc1 = coords_sub(_PSD_Center_of_rotation->_position_absolute, _DIV_Center_of_rotation->_position_absolute);
    _DIV_Center_of_rotation->_position_relative = rot_apply(_DIV_Center_of_rotation->_rotation_absolute, tc1);
  } /* DIV_Center_of_rotation=Divergence_monitor() AT ROTATED */
  DEBUG_COMPONENT("DIV_Center_of_rotation", _DIV_Center_of_rotation->_position_absolute, _DIV_Center_of_rotation->_rotation_absolute);
  instrument->_position_absolute[43] = _DIV_Center_of_rotation->_position_absolute;
  instrument->_position_relative[43] = _DIV_Center_of_rotation->_position_relative;
  instrument->counter_N[43]  = instrument->counter_P[43] = instrument->counter_P2[43] = 0;
  instrument->counter_AbsorbProp[43]= 0;
  return(0);
} /* _DIV_Center_of_rotation_setpos */

/* component Sample=PowderN() SETTING, POSITION/ROTATION */
int _Sample_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Sample_setpos] component Sample=PowderN() SETTING [/usr/share/mcstas/3.0-dev/samples/PowderN.comp:547]");
  stracpy(_Sample->_name, "Sample", 16384);
  stracpy(_Sample->_type, "PowderN", 16384);
  _Sample->_index=44;
  if("Fe.laz" && strlen("Fe.laz"))
    stracpy(_Sample->_parameters.reflections, "Fe.laz" ? "Fe.laz" : "", 16384);
  else 
  _Sample->_parameters.reflections[0]='\0';
  #define reflections (_Sample->_parameters.reflections)
  if("NULL" && strlen("NULL"))
    stracpy(_Sample->_parameters.geometry, "NULL" ? "NULL" : "", 16384);
  else 
  _Sample->_parameters.geometry[0]='\0';
  #define geometry (_Sample->_parameters.geometry)
  _Sample->_parameters.format = (double []){ 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 }; // static pointer allocation
  #define format (_Sample->_parameters.format)
  _Sample->_parameters.radius = 0.006;
  #define radius (_Sample->_parameters.radius)
  _Sample->_parameters.yheight = 0.05;
  #define yheight (_Sample->_parameters.yheight)
  _Sample->_parameters.xwidth = 0;
  #define xwidth (_Sample->_parameters.xwidth)
  _Sample->_parameters.zdepth = 0;
  #define zdepth (_Sample->_parameters.zdepth)
  _Sample->_parameters.thickness = 0;
  #define thickness (_Sample->_parameters.thickness)
  _Sample->_parameters.pack = 1;
  #define pack (_Sample->_parameters.pack)
  _Sample->_parameters.Vc = 0;
  #define Vc (_Sample->_parameters.Vc)
  _Sample->_parameters.sigma_abs = 0;
  #define sigma_abs (_Sample->_parameters.sigma_abs)
  _Sample->_parameters.sigma_inc = 0;
  #define sigma_inc (_Sample->_parameters.sigma_inc)
  _Sample->_parameters.delta_d_d = 0;
  #define delta_d_d (_Sample->_parameters.delta_d_d)
  _Sample->_parameters.p_inc = 0.1;
  #define p_inc (_Sample->_parameters.p_inc)
  _Sample->_parameters.p_transmit = 0.1;
  #define p_transmit (_Sample->_parameters.p_transmit)
  _Sample->_parameters.DW = 0;
  #define DW (_Sample->_parameters.DW)
  _Sample->_parameters.nb_atoms = 1;
  #define nb_atoms (_Sample->_parameters.nb_atoms)
  _Sample->_parameters.d_phi = 0;
  #define d_phi (_Sample->_parameters.d_phi)
  _Sample->_parameters.p_interact = 0;
  #define p_interact (_Sample->_parameters.p_interact)
  _Sample->_parameters.concentric = 0;
  #define concentric (_Sample->_parameters.concentric)
  _Sample->_parameters.density = 0;
  #define density (_Sample->_parameters.density)
  _Sample->_parameters.weight = 0;
  #define weight (_Sample->_parameters.weight)
  _Sample->_parameters.barns = 1;
  #define barns (_Sample->_parameters.barns)
  _Sample->_parameters.Strain = 0;
  #define Strain (_Sample->_parameters.Strain)
  _Sample->_parameters.focus_flip = 0;
  #define focus_flip (_Sample->_parameters.focus_flip)

  #define line_info (_Sample->_parameters.line_info)
  #define columns (_Sample->_parameters.columns)
  #define offdata (_Sample->_parameters.offdata)

  #undef reflections
  #undef geometry
  #undef format
  #undef radius
  #undef yheight
  #undef xwidth
  #undef zdepth
  #undef thickness
  #undef pack
  #undef Vc
  #undef sigma_abs
  #undef sigma_inc
  #undef delta_d_d
  #undef p_inc
  #undef p_transmit
  #undef DW
  #undef nb_atoms
  #undef d_phi
  #undef p_interact
  #undef concentric
  #undef density
  #undef weight
  #undef barns
  #undef Strain
  #undef focus_flip
  #undef line_info
  #undef columns
  #undef offdata
  /* component Sample=PowderN() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Sample_location->_rotation_absolute, _Sample->_rotation_absolute);
    rot_transpose(_DIV_Center_of_rotation->_rotation_absolute, tr1);
    rot_mul(_Sample->_rotation_absolute, tr1, _Sample->_rotation_relative);
    _Sample->_rotation_is_identity =  rot_test_identity(_Sample->_rotation_relative);
    tc1 = coords_set(
      0, 0, 0);
    rot_transpose(_Sample_location->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Sample->_position_absolute = coords_add(_Sample_location->_position_absolute, tc2);
    tc1 = coords_sub(_DIV_Center_of_rotation->_position_absolute, _Sample->_position_absolute);
    _Sample->_position_relative = rot_apply(_Sample->_rotation_absolute, tc1);
  } /* Sample=PowderN() AT ROTATED */
  DEBUG_COMPONENT("Sample", _Sample->_position_absolute, _Sample->_rotation_absolute);
  instrument->_position_absolute[44] = _Sample->_position_absolute;
  instrument->_position_relative[44] = _Sample->_position_relative;
  instrument->counter_N[44]  = instrument->counter_P[44] = instrument->counter_P2[44] = 0;
  instrument->counter_AbsorbProp[44]= 0;
  return(0);
} /* _Sample_setpos */

/* component Det_axis=Arm() SETTING, POSITION/ROTATION */
int _Det_axis_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Det_axis_setpos] component Det_axis=Arm() SETTING [Arm:0]");
  stracpy(_Det_axis->_name, "Det_axis", 16384);
  stracpy(_Det_axis->_type, "Arm", 16384);
  _Det_axis->_index=45;
  /* component Det_axis=Arm() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (instrument->_parameters._det_takeoff)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _Center_of_rotation->_rotation_absolute, _Det_axis->_rotation_absolute);
    rot_transpose(_Sample->_rotation_absolute, tr1);
    rot_mul(_Det_axis->_rotation_absolute, tr1, _Det_axis->_rotation_relative);
    _Det_axis->_rotation_is_identity =  rot_test_identity(_Det_axis->_rotation_relative);
    tc1 = coords_set(
      0, 0, 0);
    rot_transpose(_Center_of_rotation->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Det_axis->_position_absolute = coords_add(_Center_of_rotation->_position_absolute, tc2);
    tc1 = coords_sub(_Sample->_position_absolute, _Det_axis->_position_absolute);
    _Det_axis->_position_relative = rot_apply(_Det_axis->_rotation_absolute, tc1);
  } /* Det_axis=Arm() AT ROTATED */
  DEBUG_COMPONENT("Det_axis", _Det_axis->_position_absolute, _Det_axis->_rotation_absolute);
  instrument->_position_absolute[45] = _Det_axis->_position_absolute;
  instrument->_position_relative[45] = _Det_axis->_position_relative;
  instrument->counter_N[45]  = instrument->counter_P[45] = instrument->counter_P2[45] = 0;
  instrument->counter_AbsorbProp[45]= 0;
  return(0);
} /* _Det_axis_setpos */

/* component Det_axis_2=Arm() SETTING, POSITION/ROTATION */
int _Det_axis_2_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Det_axis_2_setpos] component Det_axis_2=Arm() SETTING [Arm:0]");
  stracpy(_Det_axis_2->_name, "Det_axis_2", 16384);
  stracpy(_Det_axis_2->_type, "Arm", 16384);
  _Det_axis_2->_index=46;
  /* component Det_axis_2=Arm() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (det_cover_angle)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _Det_axis->_rotation_absolute, _Det_axis_2->_rotation_absolute);
    rot_transpose(_Det_axis->_rotation_absolute, tr1);
    rot_mul(_Det_axis_2->_rotation_absolute, tr1, _Det_axis_2->_rotation_relative);
    _Det_axis_2->_rotation_is_identity =  rot_test_identity(_Det_axis_2->_rotation_relative);
    tc1 = coords_set(
      0, 0, 0);
    rot_transpose(_Det_axis->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Det_axis_2->_position_absolute = coords_add(_Det_axis->_position_absolute, tc2);
    tc1 = coords_sub(_Det_axis->_position_absolute, _Det_axis_2->_position_absolute);
    _Det_axis_2->_position_relative = rot_apply(_Det_axis_2->_rotation_absolute, tc1);
  } /* Det_axis_2=Arm() AT ROTATED */
  DEBUG_COMPONENT("Det_axis_2", _Det_axis_2->_position_absolute, _Det_axis_2->_rotation_absolute);
  instrument->_position_absolute[46] = _Det_axis_2->_position_absolute;
  instrument->_position_relative[46] = _Det_axis_2->_position_relative;
  instrument->counter_N[46]  = instrument->counter_P[46] = instrument->counter_P2[46] = 0;
  instrument->counter_AbsorbProp[46]= 0;
  return(0);
} /* _Det_axis_2_setpos */

/* component Det_axis_3=Arm() SETTING, POSITION/ROTATION */
int _Det_axis_3_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Det_axis_3_setpos] component Det_axis_3=Arm() SETTING [Arm:0]");
  stracpy(_Det_axis_3->_name, "Det_axis_3", 16384);
  stracpy(_Det_axis_3->_type, "Arm", 16384);
  _Det_axis_3->_index=47;
  /* component Det_axis_3=Arm() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (det_cover_angle)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _Det_axis_2->_rotation_absolute, _Det_axis_3->_rotation_absolute);
    rot_transpose(_Det_axis_2->_rotation_absolute, tr1);
    rot_mul(_Det_axis_3->_rotation_absolute, tr1, _Det_axis_3->_rotation_relative);
    _Det_axis_3->_rotation_is_identity =  rot_test_identity(_Det_axis_3->_rotation_relative);
    tc1 = coords_set(
      0, 0, 0);
    rot_transpose(_Det_axis_2->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Det_axis_3->_position_absolute = coords_add(_Det_axis_2->_position_absolute, tc2);
    tc1 = coords_sub(_Det_axis_2->_position_absolute, _Det_axis_3->_position_absolute);
    _Det_axis_3->_position_relative = rot_apply(_Det_axis_3->_rotation_absolute, tc1);
  } /* Det_axis_3=Arm() AT ROTATED */
  DEBUG_COMPONENT("Det_axis_3", _Det_axis_3->_position_absolute, _Det_axis_3->_rotation_absolute);
  instrument->_position_absolute[47] = _Det_axis_3->_position_absolute;
  instrument->_position_relative[47] = _Det_axis_3->_position_relative;
  instrument->counter_N[47]  = instrument->counter_P[47] = instrument->counter_P2[47] = 0;
  instrument->counter_AbsorbProp[47]= 0;
  return(0);
} /* _Det_axis_3_setpos */

/* component Det_axis_4=Arm() SETTING, POSITION/ROTATION */
int _Det_axis_4_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Det_axis_4_setpos] component Det_axis_4=Arm() SETTING [Arm:0]");
  stracpy(_Det_axis_4->_name, "Det_axis_4", 16384);
  stracpy(_Det_axis_4->_type, "Arm", 16384);
  _Det_axis_4->_index=48;
  /* component Det_axis_4=Arm() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (det_cover_angle)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _Det_axis_3->_rotation_absolute, _Det_axis_4->_rotation_absolute);
    rot_transpose(_Det_axis_3->_rotation_absolute, tr1);
    rot_mul(_Det_axis_4->_rotation_absolute, tr1, _Det_axis_4->_rotation_relative);
    _Det_axis_4->_rotation_is_identity =  rot_test_identity(_Det_axis_4->_rotation_relative);
    tc1 = coords_set(
      0, 0, 0);
    rot_transpose(_Det_axis_3->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Det_axis_4->_position_absolute = coords_add(_Det_axis_3->_position_absolute, tc2);
    tc1 = coords_sub(_Det_axis_3->_position_absolute, _Det_axis_4->_position_absolute);
    _Det_axis_4->_position_relative = rot_apply(_Det_axis_4->_rotation_absolute, tc1);
  } /* Det_axis_4=Arm() AT ROTATED */
  DEBUG_COMPONENT("Det_axis_4", _Det_axis_4->_position_absolute, _Det_axis_4->_rotation_absolute);
  instrument->_position_absolute[48] = _Det_axis_4->_position_absolute;
  instrument->_position_relative[48] = _Det_axis_4->_position_relative;
  instrument->counter_N[48]  = instrument->counter_P[48] = instrument->counter_P2[48] = 0;
  instrument->counter_AbsorbProp[48]= 0;
  return(0);
} /* _Det_axis_4_setpos */

/* component RadColl=Collimator_radial() SETTING, POSITION/ROTATION */
int _RadColl_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_RadColl_setpos] component RadColl=Collimator_radial() SETTING [/usr/share/mcstas/3.0-dev/optics/Collimator_radial.comp:80]");
  stracpy(_RadColl->_name, "RadColl", 16384);
  stracpy(_RadColl->_type, "Collimator_radial", 16384);
  _RadColl->_index=49;
  _RadColl->_parameters.xwidth = 0.006096;
  #define xwidth (_RadColl->_parameters.xwidth)
  _RadColl->_parameters.yheight = 0.147;
  #define yheight (_RadColl->_parameters.yheight)
  _RadColl->_parameters.length = 0.0863;
  #define length (_RadColl->_parameters.length)
  _RadColl->_parameters.divergence = 0;
  #define divergence (_RadColl->_parameters.divergence)
  _RadColl->_parameters.transmission = 1;
  #define transmission (_RadColl->_parameters.transmission)
  _RadColl->_parameters.theta_min = -15.625;
  #define theta_min (_RadColl->_parameters.theta_min)
  _RadColl->_parameters.theta_max = 15.625 * 7;
  #define theta_max (_RadColl->_parameters.theta_max)
  _RadColl->_parameters.nchan = 35 * 4;
  #define nchan (_RadColl->_parameters.nchan)
  _RadColl->_parameters.radius = 0.20066;
  #define radius (_RadColl->_parameters.radius)
  _RadColl->_parameters.nslit = 35 * 4;
  #define nslit (_RadColl->_parameters.nslit)
  _RadColl->_parameters.roc = 0;
  #define roc (_RadColl->_parameters.roc)
  _RadColl->_parameters.verbose = 0;
  #define verbose (_RadColl->_parameters.verbose)
  _RadColl->_parameters.approx = 0;
  #define approx (_RadColl->_parameters.approx)

  #define width_of_slit (_RadColl->_parameters.width_of_slit)
  #define width_of_Soller (_RadColl->_parameters.width_of_Soller)
  #define slit_theta (_RadColl->_parameters.slit_theta)

  #undef xwidth
  #undef yheight
  #undef length
  #undef divergence
  #undef transmission
  #undef theta_min
  #undef theta_max
  #undef nchan
  #undef radius
  #undef nslit
  #undef roc
  #undef verbose
  #undef approx
  #undef width_of_slit
  #undef width_of_Soller
  #undef slit_theta
  /* component RadColl=Collimator_radial() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Det_axis->_rotation_absolute, _RadColl->_rotation_absolute);
    rot_transpose(_Det_axis_4->_rotation_absolute, tr1);
    rot_mul(_RadColl->_rotation_absolute, tr1, _RadColl->_rotation_relative);
    _RadColl->_rotation_is_identity =  rot_test_identity(_RadColl->_rotation_relative);
    tc1 = coords_set(
      0, 0, 0);
    rot_transpose(_Det_axis->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _RadColl->_position_absolute = coords_add(_Det_axis->_position_absolute, tc2);
    tc1 = coords_sub(_Det_axis_4->_position_absolute, _RadColl->_position_absolute);
    _RadColl->_position_relative = rot_apply(_RadColl->_rotation_absolute, tc1);
  } /* RadColl=Collimator_radial() AT ROTATED */
  DEBUG_COMPONENT("RadColl", _RadColl->_position_absolute, _RadColl->_rotation_absolute);
  instrument->_position_absolute[49] = _RadColl->_position_absolute;
  instrument->_position_relative[49] = _RadColl->_position_relative;
  instrument->counter_N[49]  = instrument->counter_P[49] = instrument->counter_P2[49] = 0;
  instrument->counter_AbsorbProp[49]= 0;
  return(0);
} /* _RadColl_setpos */

/* component PSD_Detector=PSD_monitor() SETTING, POSITION/ROTATION */
int _PSD_Detector_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_PSD_Detector_setpos] component PSD_Detector=PSD_monitor() SETTING [/usr/share/mcstas/3.0-dev/monitors/PSD_monitor.comp:64]");
  stracpy(_PSD_Detector->_name, "PSD_Detector", 16384);
  stracpy(_PSD_Detector->_type, "PSD_monitor", 16384);
  _PSD_Detector->_index=50;
  _PSD_Detector->_parameters.nx = 330;
  #define nx (_PSD_Detector->_parameters.nx)
  _PSD_Detector->_parameters.ny = 15;
  #define ny (_PSD_Detector->_parameters.ny)
  if("PSD_Detector" && strlen("PSD_Detector"))
    stracpy(_PSD_Detector->_parameters.filename, "PSD_Detector" ? "PSD_Detector" : "", 16384);
  else 
  _PSD_Detector->_parameters.filename[0]='\0';
  #define filename (_PSD_Detector->_parameters.filename)
  _PSD_Detector->_parameters.xmin = -1 * det_width / 2.0;
  #define xmin (_PSD_Detector->_parameters.xmin)
  _PSD_Detector->_parameters.xmax = det_width / 2.0;
  #define xmax (_PSD_Detector->_parameters.xmax)
  _PSD_Detector->_parameters.ymin = -1 * det_height / 2.0;
  #define ymin (_PSD_Detector->_parameters.ymin)
  _PSD_Detector->_parameters.ymax = det_height / 2.0;
  #define ymax (_PSD_Detector->_parameters.ymax)
  _PSD_Detector->_parameters.xwidth = 0;
  #define xwidth (_PSD_Detector->_parameters.xwidth)
  _PSD_Detector->_parameters.yheight = 0;
  #define yheight (_PSD_Detector->_parameters.yheight)
  _PSD_Detector->_parameters.restore_neutron = 0;
  #define restore_neutron (_PSD_Detector->_parameters.restore_neutron)

  #define PSD_N (_PSD_Detector->_parameters.PSD_N)
  #define PSD_p (_PSD_Detector->_parameters.PSD_p)
  #define PSD_p2 (_PSD_Detector->_parameters.PSD_p2)

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
  /* component PSD_Detector=PSD_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Det_axis->_rotation_absolute, _PSD_Detector->_rotation_absolute);
    rot_transpose(_RadColl->_rotation_absolute, tr1);
    rot_mul(_PSD_Detector->_rotation_absolute, tr1, _PSD_Detector->_rotation_relative);
    _PSD_Detector->_rotation_is_identity =  rot_test_identity(_PSD_Detector->_rotation_relative);
    tc1 = coords_set(
      0, 0, instrument->_parameters._cor_to_det);
    rot_transpose(_Det_axis->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _PSD_Detector->_position_absolute = coords_add(_Det_axis->_position_absolute, tc2);
    tc1 = coords_sub(_RadColl->_position_absolute, _PSD_Detector->_position_absolute);
    _PSD_Detector->_position_relative = rot_apply(_PSD_Detector->_rotation_absolute, tc1);
  } /* PSD_Detector=PSD_monitor() AT ROTATED */
  DEBUG_COMPONENT("PSD_Detector", _PSD_Detector->_position_absolute, _PSD_Detector->_rotation_absolute);
  instrument->_position_absolute[50] = _PSD_Detector->_position_absolute;
  instrument->_position_relative[50] = _PSD_Detector->_position_relative;
  instrument->counter_N[50]  = instrument->counter_P[50] = instrument->counter_P2[50] = 0;
  instrument->counter_AbsorbProp[50]= 0;
  return(0);
} /* _PSD_Detector_setpos */

/* component PSD_Detector_2=PSD_monitor() SETTING, POSITION/ROTATION */
int _PSD_Detector_2_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_PSD_Detector_2_setpos] component PSD_Detector_2=PSD_monitor() SETTING [/usr/share/mcstas/3.0-dev/monitors/PSD_monitor.comp:64]");
  stracpy(_PSD_Detector_2->_name, "PSD_Detector_2", 16384);
  stracpy(_PSD_Detector_2->_type, "PSD_monitor", 16384);
  _PSD_Detector_2->_index=51;
  _PSD_Detector_2->_parameters.nx = 330;
  #define nx (_PSD_Detector_2->_parameters.nx)
  _PSD_Detector_2->_parameters.ny = 15;
  #define ny (_PSD_Detector_2->_parameters.ny)
  if("PSD_Detector_2" && strlen("PSD_Detector_2"))
    stracpy(_PSD_Detector_2->_parameters.filename, "PSD_Detector_2" ? "PSD_Detector_2" : "", 16384);
  else 
  _PSD_Detector_2->_parameters.filename[0]='\0';
  #define filename (_PSD_Detector_2->_parameters.filename)
  _PSD_Detector_2->_parameters.xmin = -1 * det_width / 2.0;
  #define xmin (_PSD_Detector_2->_parameters.xmin)
  _PSD_Detector_2->_parameters.xmax = det_width / 2.0;
  #define xmax (_PSD_Detector_2->_parameters.xmax)
  _PSD_Detector_2->_parameters.ymin = -1 * det_height / 2.0;
  #define ymin (_PSD_Detector_2->_parameters.ymin)
  _PSD_Detector_2->_parameters.ymax = det_height / 2.0;
  #define ymax (_PSD_Detector_2->_parameters.ymax)
  _PSD_Detector_2->_parameters.xwidth = 0;
  #define xwidth (_PSD_Detector_2->_parameters.xwidth)
  _PSD_Detector_2->_parameters.yheight = 0;
  #define yheight (_PSD_Detector_2->_parameters.yheight)
  _PSD_Detector_2->_parameters.restore_neutron = 0;
  #define restore_neutron (_PSD_Detector_2->_parameters.restore_neutron)

  #define PSD_N (_PSD_Detector_2->_parameters.PSD_N)
  #define PSD_p (_PSD_Detector_2->_parameters.PSD_p)
  #define PSD_p2 (_PSD_Detector_2->_parameters.PSD_p2)

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
  /* component PSD_Detector_2=PSD_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Det_axis_2->_rotation_absolute, _PSD_Detector_2->_rotation_absolute);
    rot_transpose(_PSD_Detector->_rotation_absolute, tr1);
    rot_mul(_PSD_Detector_2->_rotation_absolute, tr1, _PSD_Detector_2->_rotation_relative);
    _PSD_Detector_2->_rotation_is_identity =  rot_test_identity(_PSD_Detector_2->_rotation_relative);
    tc1 = coords_set(
      0, 0, instrument->_parameters._cor_to_det);
    rot_transpose(_Det_axis_2->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _PSD_Detector_2->_position_absolute = coords_add(_Det_axis_2->_position_absolute, tc2);
    tc1 = coords_sub(_PSD_Detector->_position_absolute, _PSD_Detector_2->_position_absolute);
    _PSD_Detector_2->_position_relative = rot_apply(_PSD_Detector_2->_rotation_absolute, tc1);
  } /* PSD_Detector_2=PSD_monitor() AT ROTATED */
  DEBUG_COMPONENT("PSD_Detector_2", _PSD_Detector_2->_position_absolute, _PSD_Detector_2->_rotation_absolute);
  instrument->_position_absolute[51] = _PSD_Detector_2->_position_absolute;
  instrument->_position_relative[51] = _PSD_Detector_2->_position_relative;
  instrument->counter_N[51]  = instrument->counter_P[51] = instrument->counter_P2[51] = 0;
  instrument->counter_AbsorbProp[51]= 0;
  return(0);
} /* _PSD_Detector_2_setpos */

/* component PSD_Detector_3=PSD_monitor() SETTING, POSITION/ROTATION */
int _PSD_Detector_3_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_PSD_Detector_3_setpos] component PSD_Detector_3=PSD_monitor() SETTING [/usr/share/mcstas/3.0-dev/monitors/PSD_monitor.comp:64]");
  stracpy(_PSD_Detector_3->_name, "PSD_Detector_3", 16384);
  stracpy(_PSD_Detector_3->_type, "PSD_monitor", 16384);
  _PSD_Detector_3->_index=52;
  _PSD_Detector_3->_parameters.nx = 330;
  #define nx (_PSD_Detector_3->_parameters.nx)
  _PSD_Detector_3->_parameters.ny = 15;
  #define ny (_PSD_Detector_3->_parameters.ny)
  if("PSD_Detector_3" && strlen("PSD_Detector_3"))
    stracpy(_PSD_Detector_3->_parameters.filename, "PSD_Detector_3" ? "PSD_Detector_3" : "", 16384);
  else 
  _PSD_Detector_3->_parameters.filename[0]='\0';
  #define filename (_PSD_Detector_3->_parameters.filename)
  _PSD_Detector_3->_parameters.xmin = -1 * det_width / 2.0;
  #define xmin (_PSD_Detector_3->_parameters.xmin)
  _PSD_Detector_3->_parameters.xmax = det_width / 2.0;
  #define xmax (_PSD_Detector_3->_parameters.xmax)
  _PSD_Detector_3->_parameters.ymin = -1 * det_height / 2.0;
  #define ymin (_PSD_Detector_3->_parameters.ymin)
  _PSD_Detector_3->_parameters.ymax = det_height / 2.0;
  #define ymax (_PSD_Detector_3->_parameters.ymax)
  _PSD_Detector_3->_parameters.xwidth = 0;
  #define xwidth (_PSD_Detector_3->_parameters.xwidth)
  _PSD_Detector_3->_parameters.yheight = 0;
  #define yheight (_PSD_Detector_3->_parameters.yheight)
  _PSD_Detector_3->_parameters.restore_neutron = 0;
  #define restore_neutron (_PSD_Detector_3->_parameters.restore_neutron)

  #define PSD_N (_PSD_Detector_3->_parameters.PSD_N)
  #define PSD_p (_PSD_Detector_3->_parameters.PSD_p)
  #define PSD_p2 (_PSD_Detector_3->_parameters.PSD_p2)

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
  /* component PSD_Detector_3=PSD_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Det_axis_3->_rotation_absolute, _PSD_Detector_3->_rotation_absolute);
    rot_transpose(_PSD_Detector_2->_rotation_absolute, tr1);
    rot_mul(_PSD_Detector_3->_rotation_absolute, tr1, _PSD_Detector_3->_rotation_relative);
    _PSD_Detector_3->_rotation_is_identity =  rot_test_identity(_PSD_Detector_3->_rotation_relative);
    tc1 = coords_set(
      0, 0, instrument->_parameters._cor_to_det);
    rot_transpose(_Det_axis_3->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _PSD_Detector_3->_position_absolute = coords_add(_Det_axis_3->_position_absolute, tc2);
    tc1 = coords_sub(_PSD_Detector_2->_position_absolute, _PSD_Detector_3->_position_absolute);
    _PSD_Detector_3->_position_relative = rot_apply(_PSD_Detector_3->_rotation_absolute, tc1);
  } /* PSD_Detector_3=PSD_monitor() AT ROTATED */
  DEBUG_COMPONENT("PSD_Detector_3", _PSD_Detector_3->_position_absolute, _PSD_Detector_3->_rotation_absolute);
  instrument->_position_absolute[52] = _PSD_Detector_3->_position_absolute;
  instrument->_position_relative[52] = _PSD_Detector_3->_position_relative;
  instrument->counter_N[52]  = instrument->counter_P[52] = instrument->counter_P2[52] = 0;
  instrument->counter_AbsorbProp[52]= 0;
  return(0);
} /* _PSD_Detector_3_setpos */

/* component PSD_Detector_4=PSD_monitor() SETTING, POSITION/ROTATION */
int _PSD_Detector_4_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_PSD_Detector_4_setpos] component PSD_Detector_4=PSD_monitor() SETTING [/usr/share/mcstas/3.0-dev/monitors/PSD_monitor.comp:64]");
  stracpy(_PSD_Detector_4->_name, "PSD_Detector_4", 16384);
  stracpy(_PSD_Detector_4->_type, "PSD_monitor", 16384);
  _PSD_Detector_4->_index=53;
  _PSD_Detector_4->_parameters.nx = 330;
  #define nx (_PSD_Detector_4->_parameters.nx)
  _PSD_Detector_4->_parameters.ny = 15;
  #define ny (_PSD_Detector_4->_parameters.ny)
  if("PSD_Detector_4" && strlen("PSD_Detector_4"))
    stracpy(_PSD_Detector_4->_parameters.filename, "PSD_Detector_4" ? "PSD_Detector_4" : "", 16384);
  else 
  _PSD_Detector_4->_parameters.filename[0]='\0';
  #define filename (_PSD_Detector_4->_parameters.filename)
  _PSD_Detector_4->_parameters.xmin = -1 * det_width / 2.0;
  #define xmin (_PSD_Detector_4->_parameters.xmin)
  _PSD_Detector_4->_parameters.xmax = det_width / 2.0;
  #define xmax (_PSD_Detector_4->_parameters.xmax)
  _PSD_Detector_4->_parameters.ymin = -1 * det_height / 2.0;
  #define ymin (_PSD_Detector_4->_parameters.ymin)
  _PSD_Detector_4->_parameters.ymax = det_height / 2.0;
  #define ymax (_PSD_Detector_4->_parameters.ymax)
  _PSD_Detector_4->_parameters.xwidth = 0;
  #define xwidth (_PSD_Detector_4->_parameters.xwidth)
  _PSD_Detector_4->_parameters.yheight = 0;
  #define yheight (_PSD_Detector_4->_parameters.yheight)
  _PSD_Detector_4->_parameters.restore_neutron = 0;
  #define restore_neutron (_PSD_Detector_4->_parameters.restore_neutron)

  #define PSD_N (_PSD_Detector_4->_parameters.PSD_N)
  #define PSD_p (_PSD_Detector_4->_parameters.PSD_p)
  #define PSD_p2 (_PSD_Detector_4->_parameters.PSD_p2)

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
  /* component PSD_Detector_4=PSD_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Det_axis_4->_rotation_absolute, _PSD_Detector_4->_rotation_absolute);
    rot_transpose(_PSD_Detector_3->_rotation_absolute, tr1);
    rot_mul(_PSD_Detector_4->_rotation_absolute, tr1, _PSD_Detector_4->_rotation_relative);
    _PSD_Detector_4->_rotation_is_identity =  rot_test_identity(_PSD_Detector_4->_rotation_relative);
    tc1 = coords_set(
      0, 0, instrument->_parameters._cor_to_det);
    rot_transpose(_Det_axis_4->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _PSD_Detector_4->_position_absolute = coords_add(_Det_axis_4->_position_absolute, tc2);
    tc1 = coords_sub(_PSD_Detector_3->_position_absolute, _PSD_Detector_4->_position_absolute);
    _PSD_Detector_4->_position_relative = rot_apply(_PSD_Detector_4->_rotation_absolute, tc1);
  } /* PSD_Detector_4=PSD_monitor() AT ROTATED */
  DEBUG_COMPONENT("PSD_Detector_4", _PSD_Detector_4->_position_absolute, _PSD_Detector_4->_rotation_absolute);
  instrument->_position_absolute[53] = _PSD_Detector_4->_position_absolute;
  instrument->_position_relative[53] = _PSD_Detector_4->_position_relative;
  instrument->counter_N[53]  = instrument->counter_P[53] = instrument->counter_P2[53] = 0;
  instrument->counter_AbsorbProp[53]= 0;
  return(0);
} /* _PSD_Detector_4_setpos */

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

_class_Filter_gen *class_Filter_gen_init(_class_Filter_gen *_comp
) {
  #define filename (_comp->_parameters.filename)
  #define options (_comp->_parameters.options)
  #define xmin (_comp->_parameters.xmin)
  #define xmax (_comp->_parameters.xmax)
  #define ymin (_comp->_parameters.ymin)
  #define ymax (_comp->_parameters.ymax)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define thickness (_comp->_parameters.thickness)
  #define scaling (_comp->_parameters.scaling)
  #define verbose (_comp->_parameters.verbose)
  #define pTable (_comp->_parameters.pTable)
  #define Mode_Table (_comp->_parameters.Mode_Table)
  #define Type_Table (_comp->_parameters.Type_Table)
  Mode_Table = FLUX_ADAPT_MULT;
  Type_Table = UNKNOWN_TABLE;

  if (xwidth  > 0) { xmax = xwidth/2;  xmin = -xmax; }
  if (yheight > 0) { ymax = yheight/2; ymin = -ymax; }

  FilterGen_Mode(options, &Mode_Table, &Type_Table, &verbose);

  if (filename != NULL && strlen(filename) && strcmp(filename,"NULL") && strcmp(filename,"0"))
  {
    if (Table_Read(&pTable, filename, 1) <= 0) /* read 1st block data from filename into pTable */
      exit(fprintf(stderr,"Filter_gen: %s: can not read filename %s\n", NAME_CURRENT_COMP, filename));

    Table_Rebin(&pTable);         /* rebin as evenly, increasing array */
    if (pTable.rows < 2 || !pTable.step_x) {
      Table_Free(&pTable);
    }
    if (pTable.data)
    {
      FilterGen_Mode(pTable.header, &Mode_Table, &Type_Table, &verbose);
      if (verbose)
      {
        Table_Info(pTable);
        printf("Filter_gen: %s: Filter data [", NAME_CURRENT_COMP);
        if (Type_Table == ENERGY_TABLE) printf("Energy");
        if (Type_Table == WAVEVECTOR_TABLE) printf("Wavevector");
        if (Type_Table == WAVELENGTH_TABLE) printf("Wavelength");
        if (Type_Table == UNKNOWN_TABLE) printf("UNKNOWN (not used)");
        printf(", Flux] in ");
        if (Mode_Table == FLUX_ADAPT_MULT) printf("multiply");
        else if (Mode_Table == FLUX_ADAPT_ADD) printf("add");
        else printf("set");
        printf(" mode\n");
      }
    } else fprintf(stderr,"Filter_gen: %s: file %s contains no data.\n", NAME_CURRENT_COMP, filename);

  } else pTable.data = NULL;
  #undef filename
  #undef options
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef xwidth
  #undef yheight
  #undef thickness
  #undef scaling
  #undef verbose
  #undef pTable
  #undef Mode_Table
  #undef Type_Table
  return(_comp);
} /* class_Filter_gen_init */

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

_class_Divergence_monitor *class_Divergence_monitor_init(_class_Divergence_monitor *_comp
) {
  #define nh (_comp->_parameters.nh)
  #define nv (_comp->_parameters.nv)
  #define filename (_comp->_parameters.filename)
  #define xmin (_comp->_parameters.xmin)
  #define xmax (_comp->_parameters.xmax)
  #define ymin (_comp->_parameters.ymin)
  #define ymax (_comp->_parameters.ymax)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define maxdiv_h (_comp->_parameters.maxdiv_h)
  #define maxdiv_v (_comp->_parameters.maxdiv_v)
  #define restore_neutron (_comp->_parameters.restore_neutron)
  #define nx (_comp->_parameters.nx)
  #define ny (_comp->_parameters.ny)
  #define nz (_comp->_parameters.nz)
  #define Div_N (_comp->_parameters.Div_N)
  #define Div_p (_comp->_parameters.Div_p)
  #define Div_p2 (_comp->_parameters.Div_p2)
  int i,j;

  if (xwidth  > 0) { xmax = xwidth/2;  xmin = -xmax; }
  if (yheight > 0) { ymax = yheight/2; ymin = -ymax; }

  if ((xmin >= xmax) || (ymin >= ymax)) {
    printf("Divergence_monitor: %s: Null detection area !\n"
           "ERROR               (xwidth,yheight,xmin,xmax,ymin,ymax). Exiting",
           NAME_CURRENT_COMP);
    exit(0);
  }

  Div_N = create_darr2d(nh, nv);
  Div_p = create_darr2d(nh, nv);
  Div_p2 = create_darr2d(nh, nv);

  for (i=0; i<nh; i++)
    for (j=0; j<nv; j++)
    {
      Div_N[i][j] = 0;
      Div_p[i][j] = 0;
      Div_p2[i][j] = 0;
    }
  NORM(nx,ny,nz);
  #undef nh
  #undef nv
  #undef filename
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef xwidth
  #undef yheight
  #undef maxdiv_h
  #undef maxdiv_v
  #undef restore_neutron
  #undef nx
  #undef ny
  #undef nz
  #undef Div_N
  #undef Div_p
  #undef Div_p2
  return(_comp);
} /* class_Divergence_monitor_init */

_class_PowderN *class_PowderN_init(_class_PowderN *_comp
) {
  #define reflections (_comp->_parameters.reflections)
  #define geometry (_comp->_parameters.geometry)
  #define format (_comp->_parameters.format)
  #define radius (_comp->_parameters.radius)
  #define yheight (_comp->_parameters.yheight)
  #define xwidth (_comp->_parameters.xwidth)
  #define zdepth (_comp->_parameters.zdepth)
  #define thickness (_comp->_parameters.thickness)
  #define pack (_comp->_parameters.pack)
  #define Vc (_comp->_parameters.Vc)
  #define sigma_abs (_comp->_parameters.sigma_abs)
  #define sigma_inc (_comp->_parameters.sigma_inc)
  #define delta_d_d (_comp->_parameters.delta_d_d)
  #define p_inc (_comp->_parameters.p_inc)
  #define p_transmit (_comp->_parameters.p_transmit)
  #define DW (_comp->_parameters.DW)
  #define nb_atoms (_comp->_parameters.nb_atoms)
  #define d_phi (_comp->_parameters.d_phi)
  #define p_interact (_comp->_parameters.p_interact)
  #define concentric (_comp->_parameters.concentric)
  #define density (_comp->_parameters.density)
  #define weight (_comp->_parameters.weight)
  #define barns (_comp->_parameters.barns)
  #define Strain (_comp->_parameters.Strain)
  #define focus_flip (_comp->_parameters.focus_flip)
  #define line_info (_comp->_parameters.line_info)
  #define columns (_comp->_parameters.columns)
  #define offdata (_comp->_parameters.offdata)
  /* We ought to clean up the columns variable as format is now a proper vector/array */
  columns = format;

  int i=0;
  struct line_data *L;
  line_info.Dd       = delta_d_d;
  line_info.DWfactor = DW;
  line_info.V_0      = Vc;
  line_info.rho      = density;
  line_info.at_weight= weight;
  line_info.at_nb    = nb_atoms;
  line_info.sigma_a  = sigma_abs;
  line_info.sigma_i  = sigma_inc;
  line_info.flag_barns=barns;
  line_info.shape    = 0;
  line_info.flag_warning=0;
  line_info.Epsilon  = Strain;
  line_info.radius_i =line_info.xwidth_i=line_info.yheight_i=line_info.zdepth_i=0;
  line_info.v  = 0;
  line_info.Nq = 0;
  line_info.v_min = FLT_MAX; line_info.v_max = 0;
  line_info.neutron_passed=0;
  line_info.nb_reuses = line_info.nb_refl = line_info.nb_refl_count = 0;
  line_info.xs_compute= line_info.xs_reuse= line_info.xs_calls =0;
  for (i=0; i< 9; i++) line_info.column_order[i] = (int)columns[i];
  strncpy(line_info.compname, NAME_CURRENT_COMP, 256);

  line_info.shape=-1; /* -1:no shape, 0:cyl, 1:box, 2:sphere, 3:any-shape  */
  if (geometry && strlen(geometry) && strcmp(geometry, "NULL") && strcmp(geometry, "0")) {
	  if (off_init(geometry, xwidth, yheight, zdepth, 0, &offdata)) {
      line_info.shape=3; thickness=0; concentric=0;
    }
  }
  else if (xwidth && yheight && zdepth)  line_info.shape=1; /* box */
  else if (radius > 0 && yheight)        line_info.shape=0; /* cylinder */
  else if (radius > 0 && !yheight)       line_info.shape=2; /* sphere */

  if (line_info.shape < 0)
    exit(fprintf(stderr,"PowderN: %s: sample has invalid dimensions.\n"
                        "ERROR    Please check parameter values (xwidth, yheight, zdepth, radius).\n", NAME_CURRENT_COMP));
  if (thickness) {
    if (radius && (radius < fabs(thickness))) {
      MPI_MASTER(
      printf("PowderN: %s: hollow sample thickness is larger than its volume (sphere/cylinder).\n"
                     "WARNING  Please check parameter values. Using bulk sample (thickness=0).\n", NAME_CURRENT_COMP);
      );
      thickness=0;
    }
    else if (!radius && (xwidth < 2*fabs(thickness) || yheight < 2*fabs(thickness) || zdepth < 2*fabs(thickness))) {
      MPI_MASTER(
      printf("PowderN: %s: hollow sample thickness is larger than its volume (box).\n"
                     "WARNING  Please check parameter values.\n", NAME_CURRENT_COMP);
      );
    }
  }

  if (concentric && thickness==0) {
    MPI_MASTER(
    printf("PowderN: %s:Can not use concentric mode\n"
           "WARNING     on non hollow shape. Ignoring.\n",
           NAME_CURRENT_COMP);
    );
    concentric=0;
  }

  if (thickness>0) {
    if (radius>thickness) {
      line_info.radius_i=radius-thickness;
    } else {
      if (xwidth>2*thickness)  line_info.xwidth_i =xwidth -2*thickness;
      if (yheight>2*thickness) line_info.yheight_i=yheight-2*thickness;
      if (zdepth>2*thickness)  line_info.zdepth_i =zdepth -2*thickness;
    }
  } else if (thickness<0) {
    thickness = fabs(thickness);
    if (radius) {
      line_info.radius_i=radius;
      radius=line_info.radius_i+thickness;
    } else {
      line_info.xwidth_i =xwidth;
      line_info.yheight_i=yheight;
      line_info.zdepth_i =zdepth;
      xwidth   =xwidth +2*thickness;
      yheight  =yheight+2*thickness;
      zdepth   =zdepth +2*thickness;
    }
  }

  if (!line_info.yheight_i) {
    line_info.yheight_i = yheight;
  }
  if (p_interact) {
    if (p_interact < p_inc) { double tmp=p_interact; p_interact=p_inc; p_inc=tmp; }
    p_transmit = 1-p_interact-p_inc;
  }

  if (p_inc + p_transmit > 1) {
    MPI_MASTER(
    printf("PowderN: %s: You have requested an unmeaningful choice of the 'p_inc' and 'p_transmit' parameters (sum is %g, exeeding 1). Fixing.\n",
                 NAME_CURRENT_COMP, p_inc+p_transmit);
    );
    if (p_inc > p_transmit) p_transmit=1-2*p_inc;
    else p_transmit=1-2*p_inc;
  } else if (p_inc + p_transmit == 1) {
    MPI_MASTER(
    printf("PowderN: %s: You have requested all neutrons be attenuated\n"
           "WARNING  or incoherently scattered!\n", NAME_CURRENT_COMP);
    );
  }

  if (concentric) {
    MPI_MASTER(
    printf("PowderN: %s: Concentric mode - remember to include the 'opposite' copy of this component !\n"
           "WARNING  The equivalent, 'opposite' comp should have concentric=0\n", NAME_CURRENT_COMP);
    );
    if (p_transmit == 0) {
      MPI_MASTER(
      printf("PowderN: %s: Concentric mode and p_transmit==0 !?\n"
             "WARNING  Don't you want any transmitted neutrons?\n", NAME_CURRENT_COMP);
      );
    }
  }

  if (reflections && strlen(reflections) && strcmp(reflections, "NULL") && strcmp(reflections, "0")) {
    i = read_line_data(reflections, &line_info);
    if (i == 0)
      exit(fprintf(stderr,"PowderN: %s: reflection file %s is not valid.\n"
                          "ERROR    Please check file format (laz or lau).\n", NAME_CURRENT_COMP, reflections));
  }

  /* compute the scattering unit density from material weight and density */
  /* the weight of the scattering element is the chemical formula molecular weight
   * times the nb of chemical formulae in the scattering element (nb_atoms) */
  if (!line_info.V_0 && line_info.at_nb > 0
    && line_info.at_weight > 0 && line_info.rho > 0) {
    /* molar volume [cm^3/mol] = weight [g/mol] / density [g/cm^3] */
    /* atom density per Angs^3 = [mol/cm^3] * N_Avogadro *(1e-8)^3 */
    line_info.V_0 = line_info.at_nb
      /(line_info.rho/line_info.at_weight/1e24*6.02214199e23);
  }

  /* the scattering unit cross sections are the chemical formula onces
   * times the nb of chemical formulae in the scattering element */
  if (line_info.at_nb > 0) {
    line_info.sigma_a *= line_info.at_nb; line_info.sigma_i *= line_info.at_nb;
  }

  if (line_info.sigma_a<0) line_info.sigma_a=0;
  if (line_info.sigma_i<0) line_info.sigma_i=0;

  if (line_info.V_0 <= 0)
  MPI_MASTER(
    printf("PowderN: %s: density/unit cell volume is NULL (Vc). Unactivating component.\n", NAME_CURRENT_COMP);
  );

  if (line_info.V_0 > 0 && p_inc && !line_info.sigma_i) {
  MPI_MASTER(
    printf("PowderN: %s: WARNING: You have requested statistics for incoherent scattering but not defined sigma_inc!\n", NAME_CURRENT_COMP);
  );
  }

  if (line_info.flag_barns) { /* Factor 100 to convert from barns to fm^2 */
    line_info.XsectionFactor = 100;
  } else {
    line_info.XsectionFactor = 1;
  }

  if (line_info.V_0 > 0 && i) {
    L = line_info.list;

    line_info.q_v = malloc(line_info.count*sizeof(double));
    line_info.w_v = malloc(line_info.count*sizeof(double));
    line_info.my_s_v2 = malloc(line_info.count*sizeof(double));
    if (!line_info.q_v || !line_info.w_v || !line_info.my_s_v2)
      exit(fprintf(stderr,"PowderN: %s: ERROR allocating memory (init)\n", NAME_CURRENT_COMP));
    for(i=0; i<line_info.count; i++)
    {
      line_info.my_s_v2[i] = 4*PI*PI*PI*pack*(L[i].DWfactor ? L[i].DWfactor : 1)
                 /(line_info.V_0*line_info.V_0*V2K*V2K)
                 *(L[i].j * L[i].F2 / L[i].q)*line_info.XsectionFactor;
      /* Is not yet divided by v^2 */
      /* Squires [3.103] */
      line_info.q_v[i] = L[i].q*K2V;
      line_info.w_v[i] = L[i].w;
    }
  }
  if (line_info.V_0 > 0) {
    /* Is not yet divided by v */
    line_info.my_a_v = pack*line_info.sigma_a/line_info.V_0*2200*100;   // Factor 100 to convert from barns to fm^2
    line_info.my_inc = pack*line_info.sigma_i/line_info.V_0*100;   // Factor 100 to convert from barns to fm^2
    MPI_MASTER(
    printf("PowderN: %s: Vc=%g [Angs] sigma_abs=%g [barn] sigma_inc=%g [barn] reflections=%s\n",
      NAME_CURRENT_COMP, line_info.V_0, line_info.sigma_a, line_info.sigma_i, reflections && strlen(reflections) ? reflections : "NULL");
    );
  }

  #undef reflections
  #undef geometry
  #undef format
  #undef radius
  #undef yheight
  #undef xwidth
  #undef zdepth
  #undef thickness
  #undef pack
  #undef Vc
  #undef sigma_abs
  #undef sigma_inc
  #undef delta_d_d
  #undef p_inc
  #undef p_transmit
  #undef DW
  #undef nb_atoms
  #undef d_phi
  #undef p_interact
  #undef concentric
  #undef density
  #undef weight
  #undef barns
  #undef Strain
  #undef focus_flip
  #undef line_info
  #undef columns
  #undef offdata
  return(_comp);
} /* class_PowderN_init */

_class_Collimator_radial *class_Collimator_radial_init(_class_Collimator_radial *_comp
) {
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define length (_comp->_parameters.length)
  #define divergence (_comp->_parameters.divergence)
  #define transmission (_comp->_parameters.transmission)
  #define theta_min (_comp->_parameters.theta_min)
  #define theta_max (_comp->_parameters.theta_max)
  #define nchan (_comp->_parameters.nchan)
  #define radius (_comp->_parameters.radius)
  #define nslit (_comp->_parameters.nslit)
  #define roc (_comp->_parameters.roc)
  #define verbose (_comp->_parameters.verbose)
  #define approx (_comp->_parameters.approx)
  #define width_of_slit (_comp->_parameters.width_of_slit)
  #define width_of_Soller (_comp->_parameters.width_of_Soller)
  #define slit_theta (_comp->_parameters.slit_theta)
width_of_slit=0;
width_of_Soller=0;
slit_theta=0;

if (radius <= 0)
    exit(printf("Collimator_radial: %s: incorrect radius=%g\n", NAME_CURRENT_COMP, radius));
  if (length <= 0)
    exit(printf("Collimator_radial: %s: invalid collimator length=%g\n", NAME_CURRENT_COMP, length));
  if (transmission <= 0 || transmission >1)
    exit(printf("Collimator_radial: %s: invalid transmission=%g\n", NAME_CURRENT_COMP, transmission));

  theta_max *= DEG2RAD;
  theta_min *= DEG2RAD;
  roc       *= DEG2RAD;
  divergence*= MIN2RAD;

  if (xwidth && !nchan)
    nchan  = ceil(radius*fabs(theta_max-theta_min)/xwidth);
  else if (!xwidth && nchan)
    xwidth = radius*fabs(theta_max-theta_min)/nchan;

  /* determine total width [m] of Soller channels, containing nslit in xwidth */
  if (nchan) width_of_Soller = radius*fabs(theta_max-theta_min)/nchan;
  else       width_of_Soller = 0;

  if (!nchan || !xwidth || xwidth > width_of_Soller)
    nchan=xwidth=width_of_Soller=0; /* continuous collimator */

  /* determine width [m] of slits */
  if (divergence) {
    width_of_slit = length*tan(divergence);
    if (xwidth)           /* Soller */
      nslit = ceil(xwidth/width_of_slit);
    else if (!nchan)      /* continuous collimator */
      nslit = ceil(radius*fabs(theta_max-theta_min)/width_of_slit);
  } else {
    if (!nchan && nslit)  /* continuous collimator */
      width_of_slit = radius*fabs(theta_max-theta_min)/nslit;
    else if (nchan && nslit)  /* Soller */
      width_of_slit = xwidth/nslit;
    divergence = atan2(width_of_slit,length);
  }

  if (nslit <= 0)
    printf("Collimator_radial: %s: number of channels must be positive nslit=%g.\n"
                "WARNING            Specify alternatively divergence, xwidth and nchan.\n",
                NAME_CURRENT_COMP, nslit);

  if (verbose && nslit && width_of_slit) {
    printf("Collimator_radial: %s: divergence %g [min] %s"
           ". Total opening [%g:%g] [deg]\n",
           NAME_CURRENT_COMP, divergence*RAD2MIN,
           (roc ? "oscillating" : ""),
           theta_min*RAD2DEG, theta_max*RAD2DEG);
    if (approx)
      printf("    Using triangular approximation model");
    else if (!nchan)
      printf("    Using continuous model");
    else
      printf("    Using %i Soller channels of width %g [cm]",
        (int)nchan, width_of_Soller*100);

    printf(" with %i slits of width %g [mm] pitch %g [deg].\n",
      (int)nslit, width_of_slit*1000, atan2(width_of_slit, radius)*RAD2DEG);
  }

  #undef xwidth
  #undef yheight
  #undef length
  #undef divergence
  #undef transmission
  #undef theta_min
  #undef theta_max
  #undef nchan
  #undef radius
  #undef nslit
  #undef roc
  #undef verbose
  #undef approx
  #undef width_of_slit
  #undef width_of_Soller
  #undef slit_theta
  return(_comp);
} /* class_Collimator_radial_init */



int init(void) { /* called by mccode_main for SAFARI_PITSI:INITIALISE */
  DEBUG_INSTR();

  /* code_main/parseoptions/readparams sets instrument parameters value */
  stracpy(instrument->_name, "SAFARI_PITSI", 256);

  /* Instrument 'SAFARI_PITSI' INITIALISE */
  SIG_MESSAGE("[SAFARI_PITSI] INITIALISE [SAFARI_PITSI.instr:124]");
  #define source_lam_min (instrument->_parameters._source_lam_min)
  #define source_lam_max (instrument->_parameters._source_lam_max)
  #define hi_res (instrument->_parameters._hi_res)
  #define mono_Si_type (instrument->_parameters._mono_Si_type)
  #define mono_mosh (instrument->_parameters._mono_mosh)
  #define mono_mosv (instrument->_parameters._mono_mosv)
  #define mono_dx (instrument->_parameters._mono_dx)
  #define mono_dy (instrument->_parameters._mono_dy)
  #define mono_dz (instrument->_parameters._mono_dz)
  #define mono_takeoff (instrument->_parameters._mono_takeoff)
  #define mono_dtilt (instrument->_parameters._mono_dtilt)
  #define mono_r_h (instrument->_parameters._mono_r_h)
  #define mono_r_v (instrument->_parameters._mono_r_v)
  #define port_takeoff (instrument->_parameters._port_takeoff)
  #define inc_slit_rot (instrument->_parameters._inc_slit_rot)
  #define inc_slit_dx (instrument->_parameters._inc_slit_dx)
  #define inc_slit_to_cor (instrument->_parameters._inc_slit_to_cor)
  #define inc_slit_width (instrument->_parameters._inc_slit_width)
  #define inc_slit_height (instrument->_parameters._inc_slit_height)
  #define inc_slit_sep (instrument->_parameters._inc_slit_sep)
  #define mono_to_cor (instrument->_parameters._mono_to_cor)
  #define sample_dx (instrument->_parameters._sample_dx)
  #define sample_dy (instrument->_parameters._sample_dy)
  #define sample_dz (instrument->_parameters._sample_dz)
  #define sample_dom (instrument->_parameters._sample_dom)
  #define det_takeoff (instrument->_parameters._det_takeoff)
  #define cor_to_det (instrument->_parameters._cor_to_det)
  #define dangle_interest (instrument->_parameters._dangle_interest)
  #define full_instrument (instrument->_parameters._full_instrument)
{
  printf ("\n------------------\n");

  /* Constants */
  if (port_takeoff == 70.0) {
    chamber_col_start=0.44057;   //from the monochromator diffraction center
    chamber_col_length=1.464;
    outside_chamber_collimator_w=0.01954;
    outside_chamber_collimator_h=0.04099;
  } else {        //90 degrees
    chamber_col_start=0.414;   //from the monochromator diffraction center
    chamber_col_length=1.3458;
    outside_chamber_collimator_w=0.02186;
    outside_chamber_collimator_h=0.04853;
  }

  /* Incident slit */
  inc_slit_xmin = -inc_slit_width/2.0;
  inc_slit_xmax = inc_slit_width/2.0;
  inc_slit_ymin = -inc_slit_height/2.0;
  inc_slit_ymax = inc_slit_height/2.0;
  if (inc_slit_sep < 0) {      //Emperical formula. Linear dependance derived from measurements: sep=13.18mm when width=0.13mm; sep=4.18mm when width=4.95mm
    inc_slit_sep= -1.86721992 * inc_slit_width + 0.0134227385;
  }
  printf ("Incident slit separation = %.4lf mm\n",inc_slit_sep*1000);

  if (inc_slit_sep == 0) {      //No gap, so horisonal and vertical slits lie ontop of each other and can have same dimensions
    inc_slit_xmin_h = inc_slit_xmin;
    inc_slit_xmax_h = inc_slit_xmax;
    inc_slit_ymin_w = inc_slit_ymin;
    inc_slit_ymax_w = inc_slit_ymax;
  } else {        //Perform some trig to be sure that none of the neutrons are accidentally cut of
    double col_to_slit = (mono_to_cor-inc_slit_to_cor) - (mono_to_cor-chamber_col_start+chamber_col_length);
    inc_slit_ymax_w = 1.01*(inc_slit_ymax + (0.023745 + inc_slit_ymax)*inc_slit_sep/col_to_slit);
    if (inc_slit_ymax_w < inc_slit_ymax) inc_slit_ymax_w = inc_slit_ymax;
    inc_slit_ymin_w = -inc_slit_ymax_w;
    inc_slit_xmax_h = -1.01*(((0.01077+inc_slit_xmax)*inc_slit_sep/col_to_slit)-inc_slit_xmax);
    if (inc_slit_xmax_h < inc_slit_xmax) inc_slit_xmax_h = inc_slit_xmax;
    inc_slit_xmin_h = - inc_slit_xmax_h;
  }
  printf ("Height slit: width=%.4lfmm, height=%.4lfmm\n",2*inc_slit_xmax_h*1000, 2*inc_slit_ymax*1000);
  printf ("Width  slit: width=%.4lfmm, height=%.4lfmm\n",2*inc_slit_xmax*1000, 2*inc_slit_ymax_w*1000);


  /*
  PoiConst = 0.2;
  miu = 0.00004;
  F = 4.1534*0.01;
  D = 0.000725;*/

  /* Monochromator */
  wafer_d = (6.0/13.0/1000.0);    //Wafer diameter. Total blade diameter (6mm) made up from 13? wafers
  start_wafer_pos = 6.0 * wafer_d;  //In order to have the middle waver situated at the centre of the port takeoff

  msw=(int)mono_Si_type;
  mono_d=0.0;
  switch (msw) {
    case 422: mono_d = 1.10858; break;
    case 400: mono_d = 1.35773; break;
    case 311: mono_d = 1.63748; break;
    case 511: mono_d = 1.04518; break;
    case 111: mono_d = 3.135; break;
    case 331: mono_d = 1.24594; break;
    case 551: mono_d = 0.76049; break;

  }

  mono_q = 2*PI/mono_d;
  double mono_omega = fabs(mono_takeoff/2.0);
  lam=2.0*mono_d*sin(DEG2RAD*mono_omega);
  printf ("mono_Si_type = %i, mono_d=%.4lfA, mono_omega = %.2lfdeg, lambda = %.4lfA\n",msw, mono_d, mono_omega, lam);

  /* Monochromator-sample distance Vertical*/
  mono_pts_v=55.55+(805.55-55.55)/(1.0/0.9)*(1.0/mono_r_v); //emperical formula that equates number of Control Pts to the radius. Derived from '13-09-23 Necsa Si Monochromator Data Sheet'
  as = -tan(2*DEG2RAD*mono_omega)/(tan(DEG2RAD*mono_omega));
  focal_dist = mono_r_v*(1-(0.5/as))*sin(DEG2RAD*mono_omega);
  printf ("Monochromator Vertical to focal point = %.4lfm (%.2lf Control Pts)\n",focal_dist, mono_pts_v);
  double mono_r_req_v = mono_to_cor/((1-(0.5/as))*sin(DEG2RAD*mono_omega));
  mono_pts_v=55.55+(805.55-55.55)/(1.0/0.9)*(1.0/mono_r_req_v); //emperical formula that equates number of Control Pts to the radius. Derived from '13-09-23 Necsa Si Monochromator Data Sheet'
  printf ("To Vertical focus on COR at %.3lfm, let mono_r_v = %.4lfm (%.2lf Pts)\n",mono_to_cor,mono_r_req_v, mono_pts_v);

  /* Monochromator-sample distance Horisontal*/
  a = 0.227642476966334;  //Gausian fitted parameters for Control Pts against curvature. y=a*exp((x-b)^2/(-2*c^2))+d
  b = 427.435873967018;  //
  c = 141.018775963938;
  d = 0.000855170045959246;
  y = 1.0/mono_r_h;    //curvature
  mono_pts_h = b + c * sqrt(2*log((y-d)/a));
  as = -tan(2*DEG2RAD*mono_omega)/(tan(DEG2RAD*mono_omega));
  focal_dist = mono_r_h*(1-(0.5/as))*sin(DEG2RAD*mono_omega);
  printf ("Monochromator Horizontan to focal point = %.4lfm (%.2lf Control Pts)\n",focal_dist, mono_pts_h);

  //mono_Rh_req = mono_to_cor/((1-(0.5/as))*sin(DEG2RAD*mono_omega));
  //mono_turns=(1005.5/mono_Rh_req)-13;    //emperical formula that equates the number of turns of the motor to the curvature. Taken from Multi-wafer silicon monochromator for stress machine at SAFARI, South Africa. Mihai Popovici
  //printf ("To focus on COR at %.3lfm, let mono_r_h = %.4lfm (%.2lf turns)\n",mono_to_cor,mono_Rh_req, mono_turns);

  /* Number of detectors to construct and calculates the angle offset between them */
  det_width=0.66;
  det_height=0.38;
  det_cover_angle = 2.0*RAD2DEG*atan((det_width/2.0)/cor_to_det);
  ndet = ceil(dangle_interest/det_cover_angle);
  printf ("Single detector coverage:%.2lf deg. %i detectors coverage angle: %.2f deg\n",\
    det_cover_angle, (int)ndet, ndet*det_cover_angle);


  printf ("------------------\n\n");
}
  #undef source_lam_min
  #undef source_lam_max
  #undef hi_res
  #undef mono_Si_type
  #undef mono_mosh
  #undef mono_mosv
  #undef mono_dx
  #undef mono_dy
  #undef mono_dz
  #undef mono_takeoff
  #undef mono_dtilt
  #undef mono_r_h
  #undef mono_r_v
  #undef port_takeoff
  #undef inc_slit_rot
  #undef inc_slit_dx
  #undef inc_slit_to_cor
  #undef inc_slit_width
  #undef inc_slit_height
  #undef inc_slit_sep
  #undef mono_to_cor
  #undef sample_dx
  #undef sample_dy
  #undef sample_dz
  #undef sample_dom
  #undef det_takeoff
  #undef cor_to_det
  #undef dangle_interest
  #undef full_instrument
  _Progress_setpos(); /* type Progress_bar */
  _Reactorbeam_setpos(); /* type Arm */
  _Prim_axes_setpos(); /* type Arm */
  _Source_setpos(); /* type Source_gen */
  _PSD_Source_setpos(); /* type PSD_monitor */
  _LAM_Source_setpos(); /* type L_monitor */
  _Window_before_filter_setpos(); /* type Slit */
  _Sapphire_filter_setpos(); /* type Filter_gen */
  _Window_after_filter_setpos(); /* type Slit */
  _LAM_After_sapphire_setpos(); /* type L_monitor */
  _PSD_After_sapphire_setpos(); /* type PSD_monitor */
  _HighResOutlet_setpos(); /* type Slit */
  _HighIntensityOutlet_setpos(); /* type Slit */
  _PSD_After_Outlet_setpos(); /* type PSD_monitor */
  _Mono_axis_setpos(); /* type Arm */
  _Blade_1_setpos(); /* type Monochromator_curved */
  _Blade_2_setpos(); /* type Monochromator_curved */
  _Blade_3_setpos(); /* type Monochromator_curved */
  _Blade_4_setpos(); /* type Monochromator_curved */
  _Blade_5_setpos(); /* type Monochromator_curved */
  _Blade_6_setpos(); /* type Monochromator_curved */
  _Blade_7_setpos(); /* type Monochromator_curved */
  _Blade_8_setpos(); /* type Monochromator_curved */
  _Blade_9_setpos(); /* type Monochromator_curved */
  _Blade_10_setpos(); /* type Monochromator_curved */
  _Blade_11_setpos(); /* type Monochromator_curved */
  _Blade_12_setpos(); /* type Monochromator_curved */
  _Blade_13_setpos(); /* type Monochromator_curved */
  _PSD_At_sec_shutter_setpos(); /* type PSD_monitor */
  _LAM_At_sec_shutter_setpos(); /* type L_monitor */
  _Inside_chamber_collimator_setpos(); /* type Slit */
  _PSD_After_inside_chamber_collimator_setpos(); /* type PSD_monitor */
  _Outside_chamber_collimator_setpos(); /* type Slit */
  _PSD_Outside_chamber_collimator_1_setpos(); /* type PSD_monitor */
  _PSD_Outside_chamber_collimator_setpos(); /* type PSD_monitor */
  _Incident_slit_h_setpos(); /* type Slit */
  _Incident_slit_w_setpos(); /* type Slit */
  _PSD_After_Incident_slit_w_setpos(); /* type PSD_monitor */
  _Center_of_rotation_setpos(); /* type Arm */
  _Sample_rotation_setpos(); /* type Arm */
  _Sample_location_setpos(); /* type Arm */
  _PSD_Center_of_rotation_setpos(); /* type PSD_monitor */
  _DIV_Center_of_rotation_setpos(); /* type Divergence_monitor */
  _Sample_setpos(); /* type PowderN */
  _Det_axis_setpos(); /* type Arm */
  _Det_axis_2_setpos(); /* type Arm */
  _Det_axis_3_setpos(); /* type Arm */
  _Det_axis_4_setpos(); /* type Arm */
  _RadColl_setpos(); /* type Collimator_radial */
  _PSD_Detector_setpos(); /* type PSD_monitor */
  _PSD_Detector_2_setpos(); /* type PSD_monitor */
  _PSD_Detector_3_setpos(); /* type PSD_monitor */
  _PSD_Detector_4_setpos(); /* type PSD_monitor */

  /* call iteratively all components INITIALISE */
  class_Progress_bar_init(_Progress);



  class_Source_gen_init(_Source);

  class_PSD_monitor_init(_PSD_Source);

  class_L_monitor_init(_LAM_Source);

  class_Slit_init(_Window_before_filter);

  class_Filter_gen_init(_Sapphire_filter);

  class_Slit_init(_Window_after_filter);

  class_L_monitor_init(_LAM_After_sapphire);

  class_PSD_monitor_init(_PSD_After_sapphire);

  class_Slit_init(_HighResOutlet);

  class_Slit_init(_HighIntensityOutlet);

  class_PSD_monitor_init(_PSD_After_Outlet);


  class_Monochromator_curved_init(_Blade_1);

  class_Monochromator_curved_init(_Blade_2);

  class_Monochromator_curved_init(_Blade_3);

  class_Monochromator_curved_init(_Blade_4);

  class_Monochromator_curved_init(_Blade_5);

  class_Monochromator_curved_init(_Blade_6);

  class_Monochromator_curved_init(_Blade_7);

  class_Monochromator_curved_init(_Blade_8);

  class_Monochromator_curved_init(_Blade_9);

  class_Monochromator_curved_init(_Blade_10);

  class_Monochromator_curved_init(_Blade_11);

  class_Monochromator_curved_init(_Blade_12);

  class_Monochromator_curved_init(_Blade_13);

  class_PSD_monitor_init(_PSD_At_sec_shutter);

  class_L_monitor_init(_LAM_At_sec_shutter);

  class_Slit_init(_Inside_chamber_collimator);

  class_PSD_monitor_init(_PSD_After_inside_chamber_collimator);

  class_Slit_init(_Outside_chamber_collimator);

  class_PSD_monitor_init(_PSD_Outside_chamber_collimator_1);

  class_PSD_monitor_init(_PSD_Outside_chamber_collimator);

  class_Slit_init(_Incident_slit_h);

  class_Slit_init(_Incident_slit_w);

  class_PSD_monitor_init(_PSD_After_Incident_slit_w);




  class_PSD_monitor_init(_PSD_Center_of_rotation);

  class_Divergence_monitor_init(_DIV_Center_of_rotation);

  class_PowderN_init(_Sample);





  class_Collimator_radial_init(_RadColl);

  class_PSD_monitor_init(_PSD_Detector);

  class_PSD_monitor_init(_PSD_Detector_2);

  class_PSD_monitor_init(_PSD_Detector_3);

  class_PSD_monitor_init(_PSD_Detector_4);

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
_class_Filter_gen *class_Filter_gen_trace(_class_Filter_gen *_comp
  , _class_particle *_particle) {
  ABSORBED=SCATTERED=RESTORE=0;

  #define filename (_comp->_parameters.filename)
  #define options (_comp->_parameters.options)
  #define xmin (_comp->_parameters.xmin)
  #define xmax (_comp->_parameters.xmax)
  #define ymin (_comp->_parameters.ymin)
  #define ymax (_comp->_parameters.ymax)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define thickness (_comp->_parameters.thickness)
  #define scaling (_comp->_parameters.scaling)
  #define verbose (_comp->_parameters.verbose)
  #define pTable (_comp->_parameters.pTable)
  #define Mode_Table (_comp->_parameters.Mode_Table)
  #define Type_Table (_comp->_parameters.Type_Table)
  double v2, K, L, E, X, new_p;

  PROP_Z0;
  if (Type_Table && (x>xmin && x<xmax && y>ymin && y<ymax))
  {
    v2 = (vx*vx + vy*vy + vz*vz);
    K = V2K*sqrt(v2);        /* k */
    L = (2*PI/K);        /* lambda */
    E = VS2E*v2;        /* energy */
    if (Type_Table == ENERGY_TABLE)     X=E;
    if (Type_Table == WAVEVECTOR_TABLE) X=K;
    if (Type_Table == WAVELENGTH_TABLE) X=L;
    /* table look up */
    if (pTable.data != NULL)
    {
      double y1, y2, x1;
      long   Index;
      Index = floor((X - pTable.min_x)/pTable.step_x);
      y1 = Table_Index(pTable, Index,   1); /* 2nd column */
      x1 = Table_Index(pTable, Index,   0); /* 1st column */
      y2 = Table_Index(pTable, Index+1, 1); /* 2nd column */
      new_p = scaling*(y1+(X - x1)*(y2-y1)/pTable.step_x); /* 2nd column */
      if (thickness != 1) new_p = pow(new_p, thickness);
    }
    else new_p = 1;

    if (Mode_Table == FLUX_ADAPT_MULT) p *= new_p;
    else p = new_p;
    SCATTER;
  }
  else
    if (Type_Table) ABSORB;
  #undef filename
  #undef options
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef xwidth
  #undef yheight
  #undef thickness
  #undef scaling
  #undef verbose
  #undef pTable
  #undef Mode_Table
  #undef Type_Table
  return(_comp);
} /* class_Filter_gen_trace */

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
_class_Divergence_monitor *class_Divergence_monitor_trace(_class_Divergence_monitor *_comp
  , _class_particle *_particle) {
  ABSORBED=SCATTERED=RESTORE=0;

  #define nh (_comp->_parameters.nh)
  #define nv (_comp->_parameters.nv)
  #define filename (_comp->_parameters.filename)
  #define xmin (_comp->_parameters.xmin)
  #define xmax (_comp->_parameters.xmax)
  #define ymin (_comp->_parameters.ymin)
  #define ymax (_comp->_parameters.ymax)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define maxdiv_h (_comp->_parameters.maxdiv_h)
  #define maxdiv_v (_comp->_parameters.maxdiv_v)
  #define restore_neutron (_comp->_parameters.restore_neutron)
  #define nx (_comp->_parameters.nx)
  #define ny (_comp->_parameters.ny)
  #define nz (_comp->_parameters.nz)
  #define Div_N (_comp->_parameters.Div_N)
  #define Div_p (_comp->_parameters.Div_p)
  #define Div_p2 (_comp->_parameters.Div_p2)
  int i,j;
  double h_div, v_div;
  double v, vn;

  PROP_Z0;
  if (x>xmin && x<xmax && y>ymin && y<ymax)
  {
    /* Find length of projection onto the [nx ny nz] axis */
    vn = scalar_prod(vx, vy, vz, nx, ny, nz);
    h_div = RAD2DEG*atan2(vx,vn);
    v_div = RAD2DEG*atan2(vy,vn);
    if (h_div < maxdiv_h && h_div > -maxdiv_h &&
        v_div < maxdiv_v && v_div > -maxdiv_v)
    {
      i = floor((h_div + maxdiv_h)*nh/(2.0*maxdiv_h));
      j = floor((v_div + maxdiv_v)*nv/(2.0*maxdiv_v));
      Div_N[i][j]++;
      Div_p[i][j] += p;
      Div_p2[i][j] += p*p;
      SCATTER;
    }
  }
  if (restore_neutron) {
    RESTORE_NEUTRON(INDEX_CURRENT_COMP, x, y, z, vx, vy, vz, t, sx, sy, sz, p);
  }
  #undef nh
  #undef nv
  #undef filename
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef xwidth
  #undef yheight
  #undef maxdiv_h
  #undef maxdiv_v
  #undef restore_neutron
  #undef nx
  #undef ny
  #undef nz
  #undef Div_N
  #undef Div_p
  #undef Div_p2
  return(_comp);
} /* class_Divergence_monitor_trace */

#pragma acc routine seq
_class_PowderN *class_PowderN_trace(_class_PowderN *_comp
  , _class_particle *_particle) {
  ABSORBED=SCATTERED=RESTORE=0;

  #define reflections (_comp->_parameters.reflections)
  #define geometry (_comp->_parameters.geometry)
  #define format (_comp->_parameters.format)
  #define radius (_comp->_parameters.radius)
  #define yheight (_comp->_parameters.yheight)
  #define xwidth (_comp->_parameters.xwidth)
  #define zdepth (_comp->_parameters.zdepth)
  #define thickness (_comp->_parameters.thickness)
  #define pack (_comp->_parameters.pack)
  #define Vc (_comp->_parameters.Vc)
  #define sigma_abs (_comp->_parameters.sigma_abs)
  #define sigma_inc (_comp->_parameters.sigma_inc)
  #define delta_d_d (_comp->_parameters.delta_d_d)
  #define p_inc (_comp->_parameters.p_inc)
  #define p_transmit (_comp->_parameters.p_transmit)
  #define DW (_comp->_parameters.DW)
  #define nb_atoms (_comp->_parameters.nb_atoms)
  #define d_phi (_comp->_parameters.d_phi)
  #define p_interact (_comp->_parameters.p_interact)
  #define concentric (_comp->_parameters.concentric)
  #define density (_comp->_parameters.density)
  #define weight (_comp->_parameters.weight)
  #define barns (_comp->_parameters.barns)
  #define Strain (_comp->_parameters.Strain)
  #define focus_flip (_comp->_parameters.focus_flip)
  #define line_info (_comp->_parameters.line_info)
  #define columns (_comp->_parameters.columns)
  #define offdata (_comp->_parameters.offdata)
  double t0, t1, t2, t3, v, v1,l_full, l, l_1, dt, alpha0, alpha, theta, my_s, my_s_n;
  double solid_angle, neutrontype;
  double arg, tmp_vx, tmp_vy, tmp_vz, vout_x, vout_y, vout_z, nx, ny, nz, pmul=1;
  int    line;
  char   intersect=0;
  char   intersecti=0;

  line_info.type = '\0';

  if (line_info.V_0 > 0 && (line_info.count || line_info.my_inc)) {
    if (line_info.shape == 1) {
      intersect  = box_intersect(&t0, &t3, x, y, z, vx, vy, vz, xwidth, yheight, zdepth);
      intersecti = box_intersect(&t1, &t2, x, y, z, vx, vy, vz, line_info.xwidth_i, line_info.yheight_i, line_info.zdepth_i);
    } else if (line_info.shape == 0) {
      intersect  = cylinder_intersect(&t0, &t3, x, y, z, vx, vy, vz, radius, yheight);
      intersecti = cylinder_intersect(&t1, &t2, x, y, z, vx, vy, vz, line_info.radius_i, line_info.yheight_i);
    } else if (line_info.shape == 2) {
      intersect  = sphere_intersect  (&t0, &t3, x,y,z, vx,vy,vz, radius);
      intersecti = sphere_intersect  (&t1, &t2, x,y,z, vx,vy,vz, line_info.radius_i);
    } else if (line_info.shape == 3) {
      intersect  = off_intersect  (&t0, &t3, NULL, NULL, x,y,z, vx,vy,vz, offdata);
      intersecti = 0;
    }
  }

  if(intersect && t3 >0) {

    if (concentric) {
      /* Set up for concentric case */
      /* 'Remove' the backside of this comp */
      if (!intersecti) {
        t1 = (t3 + t0) /2;
      }
      t2 = t1;
      t3 = t1;
      dt = -1.0*rand01(); /* In case of scattering we will scatter on 'forward' part of sample */
    } else {
      if (!intersecti) {
        t1 = (t3 + t0) /2;
        t2 = t1;
      }
      dt = randpm1(); /* Possibility to scatter at all points in line of sight */
    }

    /* Neutron enters at t=t0. */
    if(t0 < 0) t0=0; /* already in sample */
    if(t1 < 0) t1=0; /* already in inner hollow */
    if(t2 < 0) t2=0; /* already past inner hollow */
    v = sqrt(vx*vx + vy*vy + vz*vz);
    l_full = v * (t3 - t2 + t1 - t0);

    if (line_info.neutron_passed < CHAR_BUF_LENGTH) {
      if (v < line_info.v_min) line_info.v_min = v;
      if (v > line_info.v_max) line_info.v_max = v;
      line_info.neutron_passed++;
    }

    /* Calculate total scattering cross section at relevant velocity */
    if ( fabs(v - line_info.v) < 1e-6) {
        line_info.nb_reuses++;
      } else {
        line_info.Nq = calc_xsect(v, line_info.q_v, line_info.my_s_v2, line_info.count, &line_info.my_s_v2_sum, &line_info);
        line_info.v = v;
        line_info.nb_refl += line_info.Nq;
        line_info.nb_refl_count++;
      }

    if (t3 < 0) {
      t3=0; /* Already past sample?! */
      if (line_info.flag_warning < 100)
      printf("PowderN: %s: Warning: Neutron has already passed us? (Skipped).\n"
             "         In concentric geometry, this may be caused by a missing concentric=0 option in 2nd enclosing instance.\n", NAME_CURRENT_COMP);
      line_info.flag_warning++;
    } else {
      if (dt<0) { /* Calculate scattering point position */
        dt = fabs(dt)*(t1 - t0); /* 'Forward' part */
      } else {
        dt = dt * (t3 - t2) + (t2-t0) ; /* Possibly also 'backside' part */
      }

      my_s = line_info.my_s_v2_sum/(v*v)+line_info.my_inc;
      /* Total attenuation from scattering */

      neutrontype = rand01();
      /* How to handle this one? Transmit (1) / Incoherent (2) / Coherent (3) ? */
      if (neutrontype < p_transmit) {
        neutrontype = 1;
        l = l_full; /* Passing through, full length */
        PROP_DT(t3);
      } else if (neutrontype >= p_transmit && neutrontype < (p_transmit + p_inc)) {
        neutrontype = 2;
        l = v*dt;       /* Penetration in sample */
        PROP_DT(dt+t0); /* Point of scattering */
        SCATTER;
      } else if (neutrontype >= p_transmit + p_inc) {
        neutrontype = 3;
        l = v*dt;       /* Penetration in sample */
        PROP_DT(dt+t0); /* Point of scattering */
        SCATTER;
      } else {
        exit(fprintf(stderr,"PowderN %s: DEAD - this shouldn't happen!\n", NAME_CURRENT_COMP));
      }

      if (neutrontype == 3) { /* Make coherent scattering event */
        if (line_info.count > 0) {
          /* choose line */
          if (line_info.Nq > 1) line=floor(line_info.Nq*rand01());  /* Select between Nq powder lines */
          else line = 0;
          if (line_info.w_v[line])
            arg = line_info.q_v[line]*(1+line_info.w_v[line]*randnorm())/(2.0*v);
          else
            arg = line_info.q_v[line]/(2.0*v);
          my_s_n = line_info.my_s_v2[line]/(v*v);
          if(fabs(arg) > 1)
            ABSORB;                   /* No bragg scattering possible*/
          theta = asin(arg);          /* Bragg scattering law */

          /* Choose point on Debye-Scherrer cone */
          if (d_phi)
          { /* relate height of detector to the height on DS cone */
            arg = sin(d_phi*DEG2RAD/2)/sin(2*theta);
            /* If full Debye-Scherrer cone is within d_phi, don't focus */
            if (arg < -1 || arg > 1) d_phi = 0;
            /* Otherwise, determine alpha to rotate from scattering plane
               into d_phi focusing area*/
            else alpha = 2*asin(arg);
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
            if(focus_flip){
                alpha0+=M_PI_2;
            }
          }
          else
            alpha0 = PI*randpm1();

          /* now find a nearly vertical rotation axis:
           * Either
           *  (v along Z) x (X axis) -> nearly Y axis
           * Or
           *  (v along X) x (Z axis) -> nearly Y axis
           */
          if (fabs(scalar_prod(1,0,0,vx/v,vy/v,vz/v)) < fabs(scalar_prod(0,0,1,vx/v,vy/v,vz/v))) {
            nx = 1; ny = 0; nz = 0;
          } else {
            nx = 0; ny = 0; nz = 1;
          }
          vec_prod(tmp_vx,tmp_vy,tmp_vz, vx,vy,vz, nx,ny,nz);

          /* v_out = rotate 'v' by 2*theta around tmp_v: Bragg angle */
          rotate(vout_x,vout_y,vout_z, vx,vy,vz, 2*theta, tmp_vx,tmp_vy,tmp_vz);

          /* tmp_v = rotate v_out by alpha0 around 'v' (Debye-Scherrer cone) */
          rotate(tmp_vx,tmp_vy,tmp_vz, vout_x,vout_y,vout_z, alpha0, vx, vy, vz);
          vx = tmp_vx;
          vy = tmp_vy;
          vz = tmp_vz;

          /* Since now scattered and new direction given, calculate path to exit */
          if (line_info.shape == 1) {
            intersect  = box_intersect(&t0, &t3, x, y, z, vx, vy, vz, xwidth, yheight, zdepth);
            intersecti = box_intersect(&t1, &t2, x, y, z, vx, vy, vz, line_info.xwidth_i, line_info.yheight_i, line_info.zdepth_i);
          } else if (line_info.shape == 0) {
            intersect  = cylinder_intersect(&t0, &t3, x, y, z, vx, vy, vz, radius, yheight);
            intersecti = cylinder_intersect(&t1, &t2, x, y, z, vx, vy, vz, line_info.radius_i, line_info.yheight_i);
          } else if (line_info.shape == 2) {
            intersect  = sphere_intersect  (&t0, &t3, x,y,z, vx,vy,vz, radius);
            intersecti = sphere_intersect  (&t1, &t2, x,y,z, vx,vy,vz, line_info.radius_i);
          } else if (line_info.shape == 3) {
            intersect  = off_intersect  (&t0, &t3, NULL, NULL, x,y,z, vx,vy,vz, offdata);
            intersecti = 0;
          }

          if (!intersect) {
            /* Strange error: did not hit cylinder */
            if (line_info.flag_warning < 100)
              printf("PowderN: %s: WARNING: Did not hit sample from inside (coh). ABSORB.\n", NAME_CURRENT_COMP);
            line_info.flag_warning++;
            ABSORB;
          }

          if (!intersecti) {
            t1 = (t3 + t0) /2;
            t2 = t1;
          }

          if (concentric && intersecti) {
            /* In case of concentricity, 'remove' backward wall of sample */
            t2 = t1;
            t3 = t1;
          }

          if(t0 < 0) t0=0; /* already in sample */
          if(t1 < 0) t1=0; /* already in inner hollow */
          if(t2 < 0) t2=0; /* already past inner hollow */


          l_1 = v*(t3 - t2 + t1 - t0); /* Length to exit */

          pmul  = line_info.Nq*l_full*my_s_n*exp(-(line_info.my_a_v/v+my_s)*(l+l_1))
                                  /(1-(p_inc+p_transmit));
          /* Correction in case of d_phi focusing - BUT only when d_phi != 0 */
          if (d_phi) pmul *= alpha/PI;

          line_info.type = 'c';
          line_info.dq = line_info.q_v[line]*V2K;
        } /* else transmit <-- No powder lines in file */
      }  /* Coherent scattering event */
      else if (neutrontype == 2) {  /* Make incoherent scattering event */
        if(d_phi) {
          randvec_target_rect_angular(&vx, &vy, &vz, &solid_angle,
                                      0, 0, 1,
                                      2*PI, d_phi*DEG2RAD, ROT_A_CURRENT_COMP);
        } else {
          randvec_target_circle(&vx, &vy, &vz,
                                &solid_angle, 0, 0, 1, 0);
        }
        v1 = sqrt(vx*vx+vy*vy+vz*vz);
        vx *= v/v1;
        vy *= v/v1;
        vz *= v/v1;

        /* Since now scattered and new direction given, calculate path to exit */
        if (line_info.shape == 1) {
          intersect  = box_intersect(&t0, &t3, x, y, z, vx, vy, vz, xwidth, yheight, zdepth);
          intersecti = box_intersect(&t1, &t2, x, y, z, vx, vy, vz, line_info.xwidth_i, line_info.yheight_i, line_info.zdepth_i);
        } else if (line_info.shape == 0) {
          intersect  = cylinder_intersect(&t0, &t3, x, y, z, vx, vy, vz, radius, yheight);
          intersecti = cylinder_intersect(&t1, &t2, x, y, z, vx, vy, vz, line_info.radius_i, line_info.yheight_i);
        } else if (line_info.shape == 2) {
          intersect  = sphere_intersect  (&t0, &t3, x,y,z, vx,vy,vz, radius);
          intersecti = sphere_intersect  (&t1, &t2, x,y,z, vx,vy,vz, line_info.radius_i);
        } else if (line_info.shape == 3) {
          intersect  = off_intersect  (&t0, &t3, NULL, NULL, x,y,z, vx,vy,vz, offdata);
          intersecti = 0;
        }

        if (!intersect) {
          /* Strange error: did not hit cylinder */
          if (line_info.flag_warning < 100)
            printf("PowderN: %s: WARNING: Did not hit sample from inside (inc). ABSORB.\n", NAME_CURRENT_COMP);
          line_info.flag_warning++;
          ABSORB;
        }

        if (!intersecti) {
          t1 = (t3 + t0) /2;
          t2 = t1;
        }

        if (concentric && intersecti) {
          /* In case of concentricity, 'remove' backward wall of sample */
          t2 = t1;
          t3 = t1;
        }

        if(t0 < 0) t0=0; /* already in sample */
        if(t1 < 0) t1=0; /* already in inner hollow */
        if(t2 < 0) t2=0; /* already past inner hollow */


        l_1 = v*(t3 - t2 + t1 - t0); /* Length to exit */

        pmul = l_full*line_info.my_inc*exp(-(line_info.my_a_v/v+my_s)*(l+l_1))/(p_inc);
        pmul *= solid_angle/(4*PI);

        line_info.type = 'i';

      }  /* Incoherent scattering event */
      else if (neutrontype == 1) {
        /* Make transmitted (absorption-corrected) event */
        /* No coordinate changes here, simply change neutron weight */
        pmul = exp(-(line_info.my_a_v/v+my_s)*(l))/(p_transmit);

        line_info.type = 't';
      }
      p *= pmul;
    } /* Neutron leaving since it has passed already */
  } /* else transmit non interacting neutrons */


  // EXTEND code here
  if (!strcmp(_comp->_name, "Sample")) {
    if (!SCATTERED) ABSORB;
  }

  #undef reflections
  #undef geometry
  #undef format
  #undef radius
  #undef yheight
  #undef xwidth
  #undef zdepth
  #undef thickness
  #undef pack
  #undef Vc
  #undef sigma_abs
  #undef sigma_inc
  #undef delta_d_d
  #undef p_inc
  #undef p_transmit
  #undef DW
  #undef nb_atoms
  #undef d_phi
  #undef p_interact
  #undef concentric
  #undef density
  #undef weight
  #undef barns
  #undef Strain
  #undef focus_flip
  #undef line_info
  #undef columns
  #undef offdata
  return(_comp);
} /* class_PowderN_trace */

#pragma acc routine seq
_class_Collimator_radial *class_Collimator_radial_trace(_class_Collimator_radial *_comp
  , _class_particle *_particle) {
  ABSORBED=SCATTERED=RESTORE=0;

  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define length (_comp->_parameters.length)
  #define divergence (_comp->_parameters.divergence)
  #define transmission (_comp->_parameters.transmission)
  #define theta_min (_comp->_parameters.theta_min)
  #define theta_max (_comp->_parameters.theta_max)
  #define nchan (_comp->_parameters.nchan)
  #define radius (_comp->_parameters.radius)
  #define nslit (_comp->_parameters.nslit)
  #define roc (_comp->_parameters.roc)
  #define verbose (_comp->_parameters.verbose)
  #define approx (_comp->_parameters.approx)
  #define width_of_slit (_comp->_parameters.width_of_slit)
  #define width_of_Soller (_comp->_parameters.width_of_Soller)
  #define slit_theta (_comp->_parameters.slit_theta)
  double intersect=0;
  double t0, t1, t2, t3;

  if (width_of_slit && nslit) {
    /* determine intersection with inner and outer cylinders */
    intersect=cylinder_intersect(&t0,&t3,x,y,z,vx,vy,vz,radius,yheight);
    if (!intersect) ABSORB;
    else if (t3 > t0) t0 = t3;

    intersect=cylinder_intersect(&t1,&t2,x,y,z,vx,vy,vz,radius+length,yheight);
    if (!intersect) ABSORB;
    else if (t2 > t1) t1 = t2;

    /* propagate/determine neutron position at inner cylinder */
    if (t0 > 0 && t1 > t0) {
      double input_chan=0;
      double input_theta=0, output_theta=0;
      double roc_theta=0;
      double input_slit=0,  output_slit=0;

      PROP_DT(t0);

      /* apply ROC oscillation with a linear random distribution */
      if (roc) roc_theta = roc*randpm1()/2; else roc_theta=0;

      /* angle on inner cylinder */
      input_theta = atan2(x, z) + roc_theta;

      /* check if we are within min/max collimator input bounds */
      if (input_theta >= theta_min && input_theta <= theta_max) {
        SCATTER;
        /* input Soller channel index */
        if (width_of_Soller) {
          input_chan = radius*(input_theta-theta_min)/width_of_Soller;
          input_chan = input_chan-floor(input_chan); /* position within Soller [0:1] */

          /* check if we hit an absorbing housing (between Sollers): ABSORB
           *   total Soller aperture is width_of_Soller,
           *   containg a slit pack  of aperture xwidth
           */
          if (input_chan < (1-xwidth/width_of_Soller)/2
           || input_chan > (1+xwidth/width_of_Soller)/2) ABSORB;
        }

        /* determine input slit index */
        input_slit = floor(input_theta*radius/width_of_slit);

      } else /* neutron missed collimator input range */
        input_theta=4*PI;

      /* propagate to outer cylinder */
      PROP_DT(t1-t0);

      /* angle on outer cylinder */
      output_theta = atan2(x, z) + roc_theta;

      /* check if we are within min/max collimator output bounds */
      if (output_theta >= theta_min && output_theta <= theta_max) {
        /* check if we come from sides:   ABSORB */
        if (input_theta > 2*PI) ABSORB; /* input_theta=4*PI when missed input */

        if (approx) {
          double phi=atan2(x, z)-atan2(vx, vz); /* difference between positional slit angle and velocity */
          if (fabs(phi) > divergence)
            ABSORB; /* get outside transmission */
          else
            p *= (1.0 - phi/divergence);
        } else {
          /* check if we have changed slit: ABSORB */
          /* slits are considered radial so that their output size is:
             width_of_slit*(radius+length)/radius and it turns out this the same exp as for input
           */
          output_slit = floor(output_theta*radius/width_of_slit);
          if (input_slit != output_slit) ABSORB;
        }
        SCATTER;
        p *= transmission;

      } /* else neutron missed collimator output range */

    } /* else did not encounter collimator cylinders */
  } /* if nslit */

  #undef xwidth
  #undef yheight
  #undef length
  #undef divergence
  #undef transmission
  #undef theta_min
  #undef theta_max
  #undef nchan
  #undef radius
  #undef nslit
  #undef roc
  #undef verbose
  #undef approx
  #undef width_of_slit
  #undef width_of_Soller
  #undef slit_theta
  return(_comp);
} /* class_Collimator_radial_trace */

/* *****************************************************************************
* instrument 'SAFARI_PITSI' TRACE
***************************************************************************** */

#pragma acc routine seq
int raytrace(_class_particle* _particle) { /* called by mccode_main for SAFARI_PITSI:TRACE */

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
      /* component Progress=Progress_bar() [1] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Progress->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Progress->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Progress->_position_relative, _Progress->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Progress_bar_trace(_Progress, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component Progress [1] */
    if (!ABSORBED && _particle->_index == 2) {
      /* component Reactorbeam=Arm() [2] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Reactorbeam->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Reactorbeam->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Reactorbeam->_position_relative, _Reactorbeam->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      _particle->_index++;
    } /* end component Reactorbeam [2] */
    if (!ABSORBED && _particle->_index == 3) {
      /* component Prim_axes=Arm() [3] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Prim_axes->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Prim_axes->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Prim_axes->_position_relative, _Prim_axes->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      _particle->_index++;
    } /* end component Prim_axes [3] */
    if (!ABSORBED && _particle->_index == 4) {
      /* component Source=Source_gen() [4] */
  #define flux_file (_Source->_parameters.flux_file)
  #define xdiv_file (_Source->_parameters.xdiv_file)
  #define ydiv_file (_Source->_parameters.ydiv_file)
  #define radius (_Source->_parameters.radius)
  #define dist (_Source->_parameters.dist)
  #define focus_xw (_Source->_parameters.focus_xw)
  #define focus_yh (_Source->_parameters.focus_yh)
  #define focus_aw (_Source->_parameters.focus_aw)
  #define focus_ah (_Source->_parameters.focus_ah)
  #define E0 (_Source->_parameters.E0)
  #define dE (_Source->_parameters.dE)
  #define lambda0 (_Source->_parameters.lambda0)
  #define dlambda (_Source->_parameters.dlambda)
  #define I1 (_Source->_parameters.I1)
  #define yheight (_Source->_parameters.yheight)
  #define xwidth (_Source->_parameters.xwidth)
  #define verbose (_Source->_parameters.verbose)
  #define T1 (_Source->_parameters.T1)
  #define flux_file_perAA (_Source->_parameters.flux_file_perAA)
  #define flux_file_log (_Source->_parameters.flux_file_log)
  #define Lmin (_Source->_parameters.Lmin)
  #define Lmax (_Source->_parameters.Lmax)
  #define Emin (_Source->_parameters.Emin)
  #define Emax (_Source->_parameters.Emax)
  #define T2 (_Source->_parameters.T2)
  #define I2 (_Source->_parameters.I2)
  #define T3 (_Source->_parameters.T3)
  #define I3 (_Source->_parameters.I3)
  #define zdepth (_Source->_parameters.zdepth)
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
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Source->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Source->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Source->_position_relative, _Source->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      if (( instrument->_parameters._full_instrument == 1 )) // conditional WHEN execution
      class_Source_gen_trace(_Source, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
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
      _particle->_index++;
    } /* end component Source [4] */
    if (!ABSORBED && _particle->_index == 5) {
      /* component PSD_Source=PSD_monitor() [5] */
  #define nx (_PSD_Source->_parameters.nx)
  #define ny (_PSD_Source->_parameters.ny)
  #define filename (_PSD_Source->_parameters.filename)
  #define xmin (_PSD_Source->_parameters.xmin)
  #define xmax (_PSD_Source->_parameters.xmax)
  #define ymin (_PSD_Source->_parameters.ymin)
  #define ymax (_PSD_Source->_parameters.ymax)
  #define xwidth (_PSD_Source->_parameters.xwidth)
  #define yheight (_PSD_Source->_parameters.yheight)
  #define restore_neutron (_PSD_Source->_parameters.restore_neutron)
  #define PSD_N (_PSD_Source->_parameters.PSD_N)
  #define PSD_p (_PSD_Source->_parameters.PSD_p)
  #define PSD_p2 (_PSD_Source->_parameters.PSD_p2)
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_PSD_Source->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _PSD_Source->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_PSD_Source->_position_relative, _PSD_Source->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      if (( instrument->_parameters._full_instrument == 1 )) // conditional WHEN execution
      class_PSD_monitor_trace(_PSD_Source, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
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
      _particle->_index++;
    } /* end component PSD_Source [5] */
    if (!ABSORBED && _particle->_index == 6) {
      /* component LAM_Source=L_monitor() [6] */
  #define nL (_LAM_Source->_parameters.nL)
  #define filename (_LAM_Source->_parameters.filename)
  #define xmin (_LAM_Source->_parameters.xmin)
  #define xmax (_LAM_Source->_parameters.xmax)
  #define ymin (_LAM_Source->_parameters.ymin)
  #define ymax (_LAM_Source->_parameters.ymax)
  #define xwidth (_LAM_Source->_parameters.xwidth)
  #define yheight (_LAM_Source->_parameters.yheight)
  #define Lmin (_LAM_Source->_parameters.Lmin)
  #define Lmax (_LAM_Source->_parameters.Lmax)
  #define restore_neutron (_LAM_Source->_parameters.restore_neutron)
  #define L_N (_LAM_Source->_parameters.L_N)
  #define L_p (_LAM_Source->_parameters.L_p)
  #define L_p2 (_LAM_Source->_parameters.L_p2)
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_LAM_Source->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _LAM_Source->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_LAM_Source->_position_relative, _LAM_Source->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      if (( instrument->_parameters._full_instrument == 1 )) // conditional WHEN execution
      class_L_monitor_trace(_LAM_Source, _particle);
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
    } /* end component LAM_Source [6] */
    if (!ABSORBED && _particle->_index == 7) {
      /* component Window_before_filter=Slit() [7] */
  #define xmin (_Window_before_filter->_parameters.xmin)
  #define xmax (_Window_before_filter->_parameters.xmax)
  #define ymin (_Window_before_filter->_parameters.ymin)
  #define ymax (_Window_before_filter->_parameters.ymax)
  #define radius (_Window_before_filter->_parameters.radius)
  #define xwidth (_Window_before_filter->_parameters.xwidth)
  #define yheight (_Window_before_filter->_parameters.yheight)
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Window_before_filter->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Window_before_filter->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Window_before_filter->_position_relative, _Window_before_filter->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      if (( instrument->_parameters._full_instrument == 1 )) // conditional WHEN execution
      class_Slit_trace(_Window_before_filter, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef radius
  #undef xwidth
  #undef yheight
      _particle->_index++;
    } /* end component Window_before_filter [7] */
    if (!ABSORBED && _particle->_index == 8) {
      /* component Sapphire_filter=Filter_gen() [8] */
  #define filename (_Sapphire_filter->_parameters.filename)
  #define options (_Sapphire_filter->_parameters.options)
  #define xmin (_Sapphire_filter->_parameters.xmin)
  #define xmax (_Sapphire_filter->_parameters.xmax)
  #define ymin (_Sapphire_filter->_parameters.ymin)
  #define ymax (_Sapphire_filter->_parameters.ymax)
  #define xwidth (_Sapphire_filter->_parameters.xwidth)
  #define yheight (_Sapphire_filter->_parameters.yheight)
  #define thickness (_Sapphire_filter->_parameters.thickness)
  #define scaling (_Sapphire_filter->_parameters.scaling)
  #define verbose (_Sapphire_filter->_parameters.verbose)
  #define pTable (_Sapphire_filter->_parameters.pTable)
  #define Mode_Table (_Sapphire_filter->_parameters.Mode_Table)
  #define Type_Table (_Sapphire_filter->_parameters.Type_Table)
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Sapphire_filter->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Sapphire_filter->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Sapphire_filter->_position_relative, _Sapphire_filter->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      if (( instrument->_parameters._full_instrument == 1 )) // conditional WHEN execution
      class_Filter_gen_trace(_Sapphire_filter, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
  #undef filename
  #undef options
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef xwidth
  #undef yheight
  #undef thickness
  #undef scaling
  #undef verbose
  #undef pTable
  #undef Mode_Table
  #undef Type_Table
      _particle->_index++;
    } /* end component Sapphire_filter [8] */
    if (!ABSORBED && _particle->_index == 9) {
      /* component Window_after_filter=Slit() [9] */
  #define xmin (_Window_after_filter->_parameters.xmin)
  #define xmax (_Window_after_filter->_parameters.xmax)
  #define ymin (_Window_after_filter->_parameters.ymin)
  #define ymax (_Window_after_filter->_parameters.ymax)
  #define radius (_Window_after_filter->_parameters.radius)
  #define xwidth (_Window_after_filter->_parameters.xwidth)
  #define yheight (_Window_after_filter->_parameters.yheight)
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Window_after_filter->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Window_after_filter->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Window_after_filter->_position_relative, _Window_after_filter->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      if (( instrument->_parameters._full_instrument == 1 )) // conditional WHEN execution
      class_Slit_trace(_Window_after_filter, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef radius
  #undef xwidth
  #undef yheight
      _particle->_index++;
    } /* end component Window_after_filter [9] */
    if (!ABSORBED && _particle->_index == 10) {
      /* component LAM_After_sapphire=L_monitor() [10] */
  #define nL (_LAM_After_sapphire->_parameters.nL)
  #define filename (_LAM_After_sapphire->_parameters.filename)
  #define xmin (_LAM_After_sapphire->_parameters.xmin)
  #define xmax (_LAM_After_sapphire->_parameters.xmax)
  #define ymin (_LAM_After_sapphire->_parameters.ymin)
  #define ymax (_LAM_After_sapphire->_parameters.ymax)
  #define xwidth (_LAM_After_sapphire->_parameters.xwidth)
  #define yheight (_LAM_After_sapphire->_parameters.yheight)
  #define Lmin (_LAM_After_sapphire->_parameters.Lmin)
  #define Lmax (_LAM_After_sapphire->_parameters.Lmax)
  #define restore_neutron (_LAM_After_sapphire->_parameters.restore_neutron)
  #define L_N (_LAM_After_sapphire->_parameters.L_N)
  #define L_p (_LAM_After_sapphire->_parameters.L_p)
  #define L_p2 (_LAM_After_sapphire->_parameters.L_p2)
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_LAM_After_sapphire->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _LAM_After_sapphire->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_LAM_After_sapphire->_position_relative, _LAM_After_sapphire->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      if (( instrument->_parameters._full_instrument == 1 )) // conditional WHEN execution
      class_L_monitor_trace(_LAM_After_sapphire, _particle);
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
    } /* end component LAM_After_sapphire [10] */
    if (!ABSORBED && _particle->_index == 11) {
      /* component PSD_After_sapphire=PSD_monitor() [11] */
  #define nx (_PSD_After_sapphire->_parameters.nx)
  #define ny (_PSD_After_sapphire->_parameters.ny)
  #define filename (_PSD_After_sapphire->_parameters.filename)
  #define xmin (_PSD_After_sapphire->_parameters.xmin)
  #define xmax (_PSD_After_sapphire->_parameters.xmax)
  #define ymin (_PSD_After_sapphire->_parameters.ymin)
  #define ymax (_PSD_After_sapphire->_parameters.ymax)
  #define xwidth (_PSD_After_sapphire->_parameters.xwidth)
  #define yheight (_PSD_After_sapphire->_parameters.yheight)
  #define restore_neutron (_PSD_After_sapphire->_parameters.restore_neutron)
  #define PSD_N (_PSD_After_sapphire->_parameters.PSD_N)
  #define PSD_p (_PSD_After_sapphire->_parameters.PSD_p)
  #define PSD_p2 (_PSD_After_sapphire->_parameters.PSD_p2)
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_PSD_After_sapphire->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _PSD_After_sapphire->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_PSD_After_sapphire->_position_relative, _PSD_After_sapphire->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      if (( instrument->_parameters._full_instrument == 1 )) // conditional WHEN execution
      class_PSD_monitor_trace(_PSD_After_sapphire, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
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
      _particle->_index++;
    } /* end component PSD_After_sapphire [11] */
    if (!ABSORBED && _particle->_index == 12) {
      /* component HighResOutlet=Slit() [12] */
  #define xmin (_HighResOutlet->_parameters.xmin)
  #define xmax (_HighResOutlet->_parameters.xmax)
  #define ymin (_HighResOutlet->_parameters.ymin)
  #define ymax (_HighResOutlet->_parameters.ymax)
  #define radius (_HighResOutlet->_parameters.radius)
  #define xwidth (_HighResOutlet->_parameters.xwidth)
  #define yheight (_HighResOutlet->_parameters.yheight)
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_HighResOutlet->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _HighResOutlet->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_HighResOutlet->_position_relative, _HighResOutlet->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      if (( ( instrument->_parameters._hi_res == 1 ) && ( instrument->_parameters._full_instrument == 1 ) )) // conditional WHEN execution
      class_Slit_trace(_HighResOutlet, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef radius
  #undef xwidth
  #undef yheight
      _particle->_index++;
    } /* end component HighResOutlet [12] */
    if (!ABSORBED && _particle->_index == 13) {
      /* component HighIntensityOutlet=Slit() [13] */
  #define xmin (_HighIntensityOutlet->_parameters.xmin)
  #define xmax (_HighIntensityOutlet->_parameters.xmax)
  #define ymin (_HighIntensityOutlet->_parameters.ymin)
  #define ymax (_HighIntensityOutlet->_parameters.ymax)
  #define radius (_HighIntensityOutlet->_parameters.radius)
  #define xwidth (_HighIntensityOutlet->_parameters.xwidth)
  #define yheight (_HighIntensityOutlet->_parameters.yheight)
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_HighIntensityOutlet->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _HighIntensityOutlet->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_HighIntensityOutlet->_position_relative, _HighIntensityOutlet->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      if (( ( instrument->_parameters._hi_res == 0 ) && ( instrument->_parameters._full_instrument == 1 ) )) // conditional WHEN execution
      class_Slit_trace(_HighIntensityOutlet, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef radius
  #undef xwidth
  #undef yheight
      _particle->_index++;
    } /* end component HighIntensityOutlet [13] */
    if (!ABSORBED && _particle->_index == 14) {
      /* component PSD_After_Outlet=PSD_monitor() [14] */
  #define nx (_PSD_After_Outlet->_parameters.nx)
  #define ny (_PSD_After_Outlet->_parameters.ny)
  #define filename (_PSD_After_Outlet->_parameters.filename)
  #define xmin (_PSD_After_Outlet->_parameters.xmin)
  #define xmax (_PSD_After_Outlet->_parameters.xmax)
  #define ymin (_PSD_After_Outlet->_parameters.ymin)
  #define ymax (_PSD_After_Outlet->_parameters.ymax)
  #define xwidth (_PSD_After_Outlet->_parameters.xwidth)
  #define yheight (_PSD_After_Outlet->_parameters.yheight)
  #define restore_neutron (_PSD_After_Outlet->_parameters.restore_neutron)
  #define PSD_N (_PSD_After_Outlet->_parameters.PSD_N)
  #define PSD_p (_PSD_After_Outlet->_parameters.PSD_p)
  #define PSD_p2 (_PSD_After_Outlet->_parameters.PSD_p2)
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_PSD_After_Outlet->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _PSD_After_Outlet->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_PSD_After_Outlet->_position_relative, _PSD_After_Outlet->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      if (( instrument->_parameters._full_instrument == 1 )) // conditional WHEN execution
      class_PSD_monitor_trace(_PSD_After_Outlet, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
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
      _particle->_index++;
    } /* end component PSD_After_Outlet [14] */
    if (!ABSORBED && _particle->_index == 15) {
      /* component Mono_axis=Arm() [15] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Mono_axis->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Mono_axis->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Mono_axis->_position_relative, _Mono_axis->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      _particle->_index++;
    } /* end component Mono_axis [15] */
    if (!ABSORBED && _particle->_index == 16) {
      /* component Blade_1=Monochromator_curved() [16] */
  #define reflect (_Blade_1->_parameters.reflect)
  #define transmit (_Blade_1->_parameters.transmit)
  #define zwidth (_Blade_1->_parameters.zwidth)
  #define yheight (_Blade_1->_parameters.yheight)
  #define gap (_Blade_1->_parameters.gap)
  #define NH (_Blade_1->_parameters.NH)
  #define NV (_Blade_1->_parameters.NV)
  #define mosaich (_Blade_1->_parameters.mosaich)
  #define mosaicv (_Blade_1->_parameters.mosaicv)
  #define r0 (_Blade_1->_parameters.r0)
  #define t0 (_Blade_1->_parameters.t0)
  #define Q (_Blade_1->_parameters.Q)
  #define RV (_Blade_1->_parameters.RV)
  #define RH (_Blade_1->_parameters.RH)
  #define DM (_Blade_1->_parameters.DM)
  #define mosaic (_Blade_1->_parameters.mosaic)
  #define width (_Blade_1->_parameters.width)
  #define height (_Blade_1->_parameters.height)
  #define verbose (_Blade_1->_parameters.verbose)
  #define order (_Blade_1->_parameters.order)
  #define mos_rms_y (_Blade_1->_parameters.mos_rms_y)
  #define mos_rms_z (_Blade_1->_parameters.mos_rms_z)
  #define mos_rms_max (_Blade_1->_parameters.mos_rms_max)
  #define mono_Q (_Blade_1->_parameters.mono_Q)
  #define SlabWidth (_Blade_1->_parameters.SlabWidth)
  #define SlabHeight (_Blade_1->_parameters.SlabHeight)
  #define rTable (_Blade_1->_parameters.rTable)
  #define tTable (_Blade_1->_parameters.tTable)
  #define row (_Blade_1->_parameters.row)
  #define col (_Blade_1->_parameters.col)
  #define tiltH (_Blade_1->_parameters.tiltH)
  #define tiltV (_Blade_1->_parameters.tiltV)
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Blade_1->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Blade_1->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Blade_1->_position_relative, _Blade_1->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      if (( instrument->_parameters._full_instrument == 1 )) // conditional WHEN execution
      class_Monochromator_curved_trace(_Blade_1, _particle);
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
      // GROUP Monochro: from Blade_1 [16] to Blade_13 [28]
      if (SCATTERED) _particle->_index = 28; // when SCATTERED in GROUP: reach exit of GROUP after Blade_13
      else ABSORBED=0; // not SCATTERED within GROUP: always tries next
      _particle->_index++;
    } /* end component Blade_1 [16] */
    if (!ABSORBED && _particle->_index == 17) {
      /* component Blade_2=Monochromator_curved() [17] */
  #define reflect (_Blade_2->_parameters.reflect)
  #define transmit (_Blade_2->_parameters.transmit)
  #define zwidth (_Blade_2->_parameters.zwidth)
  #define yheight (_Blade_2->_parameters.yheight)
  #define gap (_Blade_2->_parameters.gap)
  #define NH (_Blade_2->_parameters.NH)
  #define NV (_Blade_2->_parameters.NV)
  #define mosaich (_Blade_2->_parameters.mosaich)
  #define mosaicv (_Blade_2->_parameters.mosaicv)
  #define r0 (_Blade_2->_parameters.r0)
  #define t0 (_Blade_2->_parameters.t0)
  #define Q (_Blade_2->_parameters.Q)
  #define RV (_Blade_2->_parameters.RV)
  #define RH (_Blade_2->_parameters.RH)
  #define DM (_Blade_2->_parameters.DM)
  #define mosaic (_Blade_2->_parameters.mosaic)
  #define width (_Blade_2->_parameters.width)
  #define height (_Blade_2->_parameters.height)
  #define verbose (_Blade_2->_parameters.verbose)
  #define order (_Blade_2->_parameters.order)
  #define mos_rms_y (_Blade_2->_parameters.mos_rms_y)
  #define mos_rms_z (_Blade_2->_parameters.mos_rms_z)
  #define mos_rms_max (_Blade_2->_parameters.mos_rms_max)
  #define mono_Q (_Blade_2->_parameters.mono_Q)
  #define SlabWidth (_Blade_2->_parameters.SlabWidth)
  #define SlabHeight (_Blade_2->_parameters.SlabHeight)
  #define rTable (_Blade_2->_parameters.rTable)
  #define tTable (_Blade_2->_parameters.tTable)
  #define row (_Blade_2->_parameters.row)
  #define col (_Blade_2->_parameters.col)
  #define tiltH (_Blade_2->_parameters.tiltH)
  #define tiltV (_Blade_2->_parameters.tiltV)
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Blade_2->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Blade_2->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Blade_2->_position_relative, _Blade_2->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      if (( instrument->_parameters._full_instrument == 1 )) // conditional WHEN execution
      class_Monochromator_curved_trace(_Blade_2, _particle);
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
      // GROUP Monochro: from Blade_1 [16] to Blade_13 [28]
      if (SCATTERED) _particle->_index = 28; // when SCATTERED in GROUP: reach exit of GROUP after Blade_13
      else ABSORBED=0; // not SCATTERED within GROUP: always tries next
      _particle->_index++;
    } /* end component Blade_2 [17] */
    if (!ABSORBED && _particle->_index == 18) {
      /* component Blade_3=Monochromator_curved() [18] */
  #define reflect (_Blade_3->_parameters.reflect)
  #define transmit (_Blade_3->_parameters.transmit)
  #define zwidth (_Blade_3->_parameters.zwidth)
  #define yheight (_Blade_3->_parameters.yheight)
  #define gap (_Blade_3->_parameters.gap)
  #define NH (_Blade_3->_parameters.NH)
  #define NV (_Blade_3->_parameters.NV)
  #define mosaich (_Blade_3->_parameters.mosaich)
  #define mosaicv (_Blade_3->_parameters.mosaicv)
  #define r0 (_Blade_3->_parameters.r0)
  #define t0 (_Blade_3->_parameters.t0)
  #define Q (_Blade_3->_parameters.Q)
  #define RV (_Blade_3->_parameters.RV)
  #define RH (_Blade_3->_parameters.RH)
  #define DM (_Blade_3->_parameters.DM)
  #define mosaic (_Blade_3->_parameters.mosaic)
  #define width (_Blade_3->_parameters.width)
  #define height (_Blade_3->_parameters.height)
  #define verbose (_Blade_3->_parameters.verbose)
  #define order (_Blade_3->_parameters.order)
  #define mos_rms_y (_Blade_3->_parameters.mos_rms_y)
  #define mos_rms_z (_Blade_3->_parameters.mos_rms_z)
  #define mos_rms_max (_Blade_3->_parameters.mos_rms_max)
  #define mono_Q (_Blade_3->_parameters.mono_Q)
  #define SlabWidth (_Blade_3->_parameters.SlabWidth)
  #define SlabHeight (_Blade_3->_parameters.SlabHeight)
  #define rTable (_Blade_3->_parameters.rTable)
  #define tTable (_Blade_3->_parameters.tTable)
  #define row (_Blade_3->_parameters.row)
  #define col (_Blade_3->_parameters.col)
  #define tiltH (_Blade_3->_parameters.tiltH)
  #define tiltV (_Blade_3->_parameters.tiltV)
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Blade_3->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Blade_3->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Blade_3->_position_relative, _Blade_3->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      if (( instrument->_parameters._full_instrument == 1 )) // conditional WHEN execution
      class_Monochromator_curved_trace(_Blade_3, _particle);
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
      // GROUP Monochro: from Blade_1 [16] to Blade_13 [28]
      if (SCATTERED) _particle->_index = 28; // when SCATTERED in GROUP: reach exit of GROUP after Blade_13
      else ABSORBED=0; // not SCATTERED within GROUP: always tries next
      _particle->_index++;
    } /* end component Blade_3 [18] */
    if (!ABSORBED && _particle->_index == 19) {
      /* component Blade_4=Monochromator_curved() [19] */
  #define reflect (_Blade_4->_parameters.reflect)
  #define transmit (_Blade_4->_parameters.transmit)
  #define zwidth (_Blade_4->_parameters.zwidth)
  #define yheight (_Blade_4->_parameters.yheight)
  #define gap (_Blade_4->_parameters.gap)
  #define NH (_Blade_4->_parameters.NH)
  #define NV (_Blade_4->_parameters.NV)
  #define mosaich (_Blade_4->_parameters.mosaich)
  #define mosaicv (_Blade_4->_parameters.mosaicv)
  #define r0 (_Blade_4->_parameters.r0)
  #define t0 (_Blade_4->_parameters.t0)
  #define Q (_Blade_4->_parameters.Q)
  #define RV (_Blade_4->_parameters.RV)
  #define RH (_Blade_4->_parameters.RH)
  #define DM (_Blade_4->_parameters.DM)
  #define mosaic (_Blade_4->_parameters.mosaic)
  #define width (_Blade_4->_parameters.width)
  #define height (_Blade_4->_parameters.height)
  #define verbose (_Blade_4->_parameters.verbose)
  #define order (_Blade_4->_parameters.order)
  #define mos_rms_y (_Blade_4->_parameters.mos_rms_y)
  #define mos_rms_z (_Blade_4->_parameters.mos_rms_z)
  #define mos_rms_max (_Blade_4->_parameters.mos_rms_max)
  #define mono_Q (_Blade_4->_parameters.mono_Q)
  #define SlabWidth (_Blade_4->_parameters.SlabWidth)
  #define SlabHeight (_Blade_4->_parameters.SlabHeight)
  #define rTable (_Blade_4->_parameters.rTable)
  #define tTable (_Blade_4->_parameters.tTable)
  #define row (_Blade_4->_parameters.row)
  #define col (_Blade_4->_parameters.col)
  #define tiltH (_Blade_4->_parameters.tiltH)
  #define tiltV (_Blade_4->_parameters.tiltV)
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Blade_4->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Blade_4->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Blade_4->_position_relative, _Blade_4->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      if (( instrument->_parameters._full_instrument == 1 )) // conditional WHEN execution
      class_Monochromator_curved_trace(_Blade_4, _particle);
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
      // GROUP Monochro: from Blade_1 [16] to Blade_13 [28]
      if (SCATTERED) _particle->_index = 28; // when SCATTERED in GROUP: reach exit of GROUP after Blade_13
      else ABSORBED=0; // not SCATTERED within GROUP: always tries next
      _particle->_index++;
    } /* end component Blade_4 [19] */
    if (!ABSORBED && _particle->_index == 20) {
      /* component Blade_5=Monochromator_curved() [20] */
  #define reflect (_Blade_5->_parameters.reflect)
  #define transmit (_Blade_5->_parameters.transmit)
  #define zwidth (_Blade_5->_parameters.zwidth)
  #define yheight (_Blade_5->_parameters.yheight)
  #define gap (_Blade_5->_parameters.gap)
  #define NH (_Blade_5->_parameters.NH)
  #define NV (_Blade_5->_parameters.NV)
  #define mosaich (_Blade_5->_parameters.mosaich)
  #define mosaicv (_Blade_5->_parameters.mosaicv)
  #define r0 (_Blade_5->_parameters.r0)
  #define t0 (_Blade_5->_parameters.t0)
  #define Q (_Blade_5->_parameters.Q)
  #define RV (_Blade_5->_parameters.RV)
  #define RH (_Blade_5->_parameters.RH)
  #define DM (_Blade_5->_parameters.DM)
  #define mosaic (_Blade_5->_parameters.mosaic)
  #define width (_Blade_5->_parameters.width)
  #define height (_Blade_5->_parameters.height)
  #define verbose (_Blade_5->_parameters.verbose)
  #define order (_Blade_5->_parameters.order)
  #define mos_rms_y (_Blade_5->_parameters.mos_rms_y)
  #define mos_rms_z (_Blade_5->_parameters.mos_rms_z)
  #define mos_rms_max (_Blade_5->_parameters.mos_rms_max)
  #define mono_Q (_Blade_5->_parameters.mono_Q)
  #define SlabWidth (_Blade_5->_parameters.SlabWidth)
  #define SlabHeight (_Blade_5->_parameters.SlabHeight)
  #define rTable (_Blade_5->_parameters.rTable)
  #define tTable (_Blade_5->_parameters.tTable)
  #define row (_Blade_5->_parameters.row)
  #define col (_Blade_5->_parameters.col)
  #define tiltH (_Blade_5->_parameters.tiltH)
  #define tiltV (_Blade_5->_parameters.tiltV)
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Blade_5->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Blade_5->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Blade_5->_position_relative, _Blade_5->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      if (( instrument->_parameters._full_instrument == 1 )) // conditional WHEN execution
      class_Monochromator_curved_trace(_Blade_5, _particle);
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
      // GROUP Monochro: from Blade_1 [16] to Blade_13 [28]
      if (SCATTERED) _particle->_index = 28; // when SCATTERED in GROUP: reach exit of GROUP after Blade_13
      else ABSORBED=0; // not SCATTERED within GROUP: always tries next
      _particle->_index++;
    } /* end component Blade_5 [20] */
    if (!ABSORBED && _particle->_index == 21) {
      /* component Blade_6=Monochromator_curved() [21] */
  #define reflect (_Blade_6->_parameters.reflect)
  #define transmit (_Blade_6->_parameters.transmit)
  #define zwidth (_Blade_6->_parameters.zwidth)
  #define yheight (_Blade_6->_parameters.yheight)
  #define gap (_Blade_6->_parameters.gap)
  #define NH (_Blade_6->_parameters.NH)
  #define NV (_Blade_6->_parameters.NV)
  #define mosaich (_Blade_6->_parameters.mosaich)
  #define mosaicv (_Blade_6->_parameters.mosaicv)
  #define r0 (_Blade_6->_parameters.r0)
  #define t0 (_Blade_6->_parameters.t0)
  #define Q (_Blade_6->_parameters.Q)
  #define RV (_Blade_6->_parameters.RV)
  #define RH (_Blade_6->_parameters.RH)
  #define DM (_Blade_6->_parameters.DM)
  #define mosaic (_Blade_6->_parameters.mosaic)
  #define width (_Blade_6->_parameters.width)
  #define height (_Blade_6->_parameters.height)
  #define verbose (_Blade_6->_parameters.verbose)
  #define order (_Blade_6->_parameters.order)
  #define mos_rms_y (_Blade_6->_parameters.mos_rms_y)
  #define mos_rms_z (_Blade_6->_parameters.mos_rms_z)
  #define mos_rms_max (_Blade_6->_parameters.mos_rms_max)
  #define mono_Q (_Blade_6->_parameters.mono_Q)
  #define SlabWidth (_Blade_6->_parameters.SlabWidth)
  #define SlabHeight (_Blade_6->_parameters.SlabHeight)
  #define rTable (_Blade_6->_parameters.rTable)
  #define tTable (_Blade_6->_parameters.tTable)
  #define row (_Blade_6->_parameters.row)
  #define col (_Blade_6->_parameters.col)
  #define tiltH (_Blade_6->_parameters.tiltH)
  #define tiltV (_Blade_6->_parameters.tiltV)
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Blade_6->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Blade_6->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Blade_6->_position_relative, _Blade_6->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      if (( instrument->_parameters._full_instrument == 1 )) // conditional WHEN execution
      class_Monochromator_curved_trace(_Blade_6, _particle);
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
      // GROUP Monochro: from Blade_1 [16] to Blade_13 [28]
      if (SCATTERED) _particle->_index = 28; // when SCATTERED in GROUP: reach exit of GROUP after Blade_13
      else ABSORBED=0; // not SCATTERED within GROUP: always tries next
      _particle->_index++;
    } /* end component Blade_6 [21] */
    if (!ABSORBED && _particle->_index == 22) {
      /* component Blade_7=Monochromator_curved() [22] */
  #define reflect (_Blade_7->_parameters.reflect)
  #define transmit (_Blade_7->_parameters.transmit)
  #define zwidth (_Blade_7->_parameters.zwidth)
  #define yheight (_Blade_7->_parameters.yheight)
  #define gap (_Blade_7->_parameters.gap)
  #define NH (_Blade_7->_parameters.NH)
  #define NV (_Blade_7->_parameters.NV)
  #define mosaich (_Blade_7->_parameters.mosaich)
  #define mosaicv (_Blade_7->_parameters.mosaicv)
  #define r0 (_Blade_7->_parameters.r0)
  #define t0 (_Blade_7->_parameters.t0)
  #define Q (_Blade_7->_parameters.Q)
  #define RV (_Blade_7->_parameters.RV)
  #define RH (_Blade_7->_parameters.RH)
  #define DM (_Blade_7->_parameters.DM)
  #define mosaic (_Blade_7->_parameters.mosaic)
  #define width (_Blade_7->_parameters.width)
  #define height (_Blade_7->_parameters.height)
  #define verbose (_Blade_7->_parameters.verbose)
  #define order (_Blade_7->_parameters.order)
  #define mos_rms_y (_Blade_7->_parameters.mos_rms_y)
  #define mos_rms_z (_Blade_7->_parameters.mos_rms_z)
  #define mos_rms_max (_Blade_7->_parameters.mos_rms_max)
  #define mono_Q (_Blade_7->_parameters.mono_Q)
  #define SlabWidth (_Blade_7->_parameters.SlabWidth)
  #define SlabHeight (_Blade_7->_parameters.SlabHeight)
  #define rTable (_Blade_7->_parameters.rTable)
  #define tTable (_Blade_7->_parameters.tTable)
  #define row (_Blade_7->_parameters.row)
  #define col (_Blade_7->_parameters.col)
  #define tiltH (_Blade_7->_parameters.tiltH)
  #define tiltV (_Blade_7->_parameters.tiltV)
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Blade_7->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Blade_7->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Blade_7->_position_relative, _Blade_7->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      if (( instrument->_parameters._full_instrument == 1 )) // conditional WHEN execution
      class_Monochromator_curved_trace(_Blade_7, _particle);
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
      // GROUP Monochro: from Blade_1 [16] to Blade_13 [28]
      if (SCATTERED) _particle->_index = 28; // when SCATTERED in GROUP: reach exit of GROUP after Blade_13
      else ABSORBED=0; // not SCATTERED within GROUP: always tries next
      _particle->_index++;
    } /* end component Blade_7 [22] */
    if (!ABSORBED && _particle->_index == 23) {
      /* component Blade_8=Monochromator_curved() [23] */
  #define reflect (_Blade_8->_parameters.reflect)
  #define transmit (_Blade_8->_parameters.transmit)
  #define zwidth (_Blade_8->_parameters.zwidth)
  #define yheight (_Blade_8->_parameters.yheight)
  #define gap (_Blade_8->_parameters.gap)
  #define NH (_Blade_8->_parameters.NH)
  #define NV (_Blade_8->_parameters.NV)
  #define mosaich (_Blade_8->_parameters.mosaich)
  #define mosaicv (_Blade_8->_parameters.mosaicv)
  #define r0 (_Blade_8->_parameters.r0)
  #define t0 (_Blade_8->_parameters.t0)
  #define Q (_Blade_8->_parameters.Q)
  #define RV (_Blade_8->_parameters.RV)
  #define RH (_Blade_8->_parameters.RH)
  #define DM (_Blade_8->_parameters.DM)
  #define mosaic (_Blade_8->_parameters.mosaic)
  #define width (_Blade_8->_parameters.width)
  #define height (_Blade_8->_parameters.height)
  #define verbose (_Blade_8->_parameters.verbose)
  #define order (_Blade_8->_parameters.order)
  #define mos_rms_y (_Blade_8->_parameters.mos_rms_y)
  #define mos_rms_z (_Blade_8->_parameters.mos_rms_z)
  #define mos_rms_max (_Blade_8->_parameters.mos_rms_max)
  #define mono_Q (_Blade_8->_parameters.mono_Q)
  #define SlabWidth (_Blade_8->_parameters.SlabWidth)
  #define SlabHeight (_Blade_8->_parameters.SlabHeight)
  #define rTable (_Blade_8->_parameters.rTable)
  #define tTable (_Blade_8->_parameters.tTable)
  #define row (_Blade_8->_parameters.row)
  #define col (_Blade_8->_parameters.col)
  #define tiltH (_Blade_8->_parameters.tiltH)
  #define tiltV (_Blade_8->_parameters.tiltV)
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Blade_8->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Blade_8->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Blade_8->_position_relative, _Blade_8->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      if (( instrument->_parameters._full_instrument == 1 )) // conditional WHEN execution
      class_Monochromator_curved_trace(_Blade_8, _particle);
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
      // GROUP Monochro: from Blade_1 [16] to Blade_13 [28]
      if (SCATTERED) _particle->_index = 28; // when SCATTERED in GROUP: reach exit of GROUP after Blade_13
      else ABSORBED=0; // not SCATTERED within GROUP: always tries next
      _particle->_index++;
    } /* end component Blade_8 [23] */
    if (!ABSORBED && _particle->_index == 24) {
      /* component Blade_9=Monochromator_curved() [24] */
  #define reflect (_Blade_9->_parameters.reflect)
  #define transmit (_Blade_9->_parameters.transmit)
  #define zwidth (_Blade_9->_parameters.zwidth)
  #define yheight (_Blade_9->_parameters.yheight)
  #define gap (_Blade_9->_parameters.gap)
  #define NH (_Blade_9->_parameters.NH)
  #define NV (_Blade_9->_parameters.NV)
  #define mosaich (_Blade_9->_parameters.mosaich)
  #define mosaicv (_Blade_9->_parameters.mosaicv)
  #define r0 (_Blade_9->_parameters.r0)
  #define t0 (_Blade_9->_parameters.t0)
  #define Q (_Blade_9->_parameters.Q)
  #define RV (_Blade_9->_parameters.RV)
  #define RH (_Blade_9->_parameters.RH)
  #define DM (_Blade_9->_parameters.DM)
  #define mosaic (_Blade_9->_parameters.mosaic)
  #define width (_Blade_9->_parameters.width)
  #define height (_Blade_9->_parameters.height)
  #define verbose (_Blade_9->_parameters.verbose)
  #define order (_Blade_9->_parameters.order)
  #define mos_rms_y (_Blade_9->_parameters.mos_rms_y)
  #define mos_rms_z (_Blade_9->_parameters.mos_rms_z)
  #define mos_rms_max (_Blade_9->_parameters.mos_rms_max)
  #define mono_Q (_Blade_9->_parameters.mono_Q)
  #define SlabWidth (_Blade_9->_parameters.SlabWidth)
  #define SlabHeight (_Blade_9->_parameters.SlabHeight)
  #define rTable (_Blade_9->_parameters.rTable)
  #define tTable (_Blade_9->_parameters.tTable)
  #define row (_Blade_9->_parameters.row)
  #define col (_Blade_9->_parameters.col)
  #define tiltH (_Blade_9->_parameters.tiltH)
  #define tiltV (_Blade_9->_parameters.tiltV)
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Blade_9->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Blade_9->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Blade_9->_position_relative, _Blade_9->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      if (( instrument->_parameters._full_instrument == 1 )) // conditional WHEN execution
      class_Monochromator_curved_trace(_Blade_9, _particle);
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
      // GROUP Monochro: from Blade_1 [16] to Blade_13 [28]
      if (SCATTERED) _particle->_index = 28; // when SCATTERED in GROUP: reach exit of GROUP after Blade_13
      else ABSORBED=0; // not SCATTERED within GROUP: always tries next
      _particle->_index++;
    } /* end component Blade_9 [24] */
    if (!ABSORBED && _particle->_index == 25) {
      /* component Blade_10=Monochromator_curved() [25] */
  #define reflect (_Blade_10->_parameters.reflect)
  #define transmit (_Blade_10->_parameters.transmit)
  #define zwidth (_Blade_10->_parameters.zwidth)
  #define yheight (_Blade_10->_parameters.yheight)
  #define gap (_Blade_10->_parameters.gap)
  #define NH (_Blade_10->_parameters.NH)
  #define NV (_Blade_10->_parameters.NV)
  #define mosaich (_Blade_10->_parameters.mosaich)
  #define mosaicv (_Blade_10->_parameters.mosaicv)
  #define r0 (_Blade_10->_parameters.r0)
  #define t0 (_Blade_10->_parameters.t0)
  #define Q (_Blade_10->_parameters.Q)
  #define RV (_Blade_10->_parameters.RV)
  #define RH (_Blade_10->_parameters.RH)
  #define DM (_Blade_10->_parameters.DM)
  #define mosaic (_Blade_10->_parameters.mosaic)
  #define width (_Blade_10->_parameters.width)
  #define height (_Blade_10->_parameters.height)
  #define verbose (_Blade_10->_parameters.verbose)
  #define order (_Blade_10->_parameters.order)
  #define mos_rms_y (_Blade_10->_parameters.mos_rms_y)
  #define mos_rms_z (_Blade_10->_parameters.mos_rms_z)
  #define mos_rms_max (_Blade_10->_parameters.mos_rms_max)
  #define mono_Q (_Blade_10->_parameters.mono_Q)
  #define SlabWidth (_Blade_10->_parameters.SlabWidth)
  #define SlabHeight (_Blade_10->_parameters.SlabHeight)
  #define rTable (_Blade_10->_parameters.rTable)
  #define tTable (_Blade_10->_parameters.tTable)
  #define row (_Blade_10->_parameters.row)
  #define col (_Blade_10->_parameters.col)
  #define tiltH (_Blade_10->_parameters.tiltH)
  #define tiltV (_Blade_10->_parameters.tiltV)
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Blade_10->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Blade_10->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Blade_10->_position_relative, _Blade_10->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      if (( instrument->_parameters._full_instrument == 1 )) // conditional WHEN execution
      class_Monochromator_curved_trace(_Blade_10, _particle);
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
      // GROUP Monochro: from Blade_1 [16] to Blade_13 [28]
      if (SCATTERED) _particle->_index = 28; // when SCATTERED in GROUP: reach exit of GROUP after Blade_13
      else ABSORBED=0; // not SCATTERED within GROUP: always tries next
      _particle->_index++;
    } /* end component Blade_10 [25] */
    if (!ABSORBED && _particle->_index == 26) {
      /* component Blade_11=Monochromator_curved() [26] */
  #define reflect (_Blade_11->_parameters.reflect)
  #define transmit (_Blade_11->_parameters.transmit)
  #define zwidth (_Blade_11->_parameters.zwidth)
  #define yheight (_Blade_11->_parameters.yheight)
  #define gap (_Blade_11->_parameters.gap)
  #define NH (_Blade_11->_parameters.NH)
  #define NV (_Blade_11->_parameters.NV)
  #define mosaich (_Blade_11->_parameters.mosaich)
  #define mosaicv (_Blade_11->_parameters.mosaicv)
  #define r0 (_Blade_11->_parameters.r0)
  #define t0 (_Blade_11->_parameters.t0)
  #define Q (_Blade_11->_parameters.Q)
  #define RV (_Blade_11->_parameters.RV)
  #define RH (_Blade_11->_parameters.RH)
  #define DM (_Blade_11->_parameters.DM)
  #define mosaic (_Blade_11->_parameters.mosaic)
  #define width (_Blade_11->_parameters.width)
  #define height (_Blade_11->_parameters.height)
  #define verbose (_Blade_11->_parameters.verbose)
  #define order (_Blade_11->_parameters.order)
  #define mos_rms_y (_Blade_11->_parameters.mos_rms_y)
  #define mos_rms_z (_Blade_11->_parameters.mos_rms_z)
  #define mos_rms_max (_Blade_11->_parameters.mos_rms_max)
  #define mono_Q (_Blade_11->_parameters.mono_Q)
  #define SlabWidth (_Blade_11->_parameters.SlabWidth)
  #define SlabHeight (_Blade_11->_parameters.SlabHeight)
  #define rTable (_Blade_11->_parameters.rTable)
  #define tTable (_Blade_11->_parameters.tTable)
  #define row (_Blade_11->_parameters.row)
  #define col (_Blade_11->_parameters.col)
  #define tiltH (_Blade_11->_parameters.tiltH)
  #define tiltV (_Blade_11->_parameters.tiltV)
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Blade_11->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Blade_11->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Blade_11->_position_relative, _Blade_11->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      if (( instrument->_parameters._full_instrument == 1 )) // conditional WHEN execution
      class_Monochromator_curved_trace(_Blade_11, _particle);
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
      // GROUP Monochro: from Blade_1 [16] to Blade_13 [28]
      if (SCATTERED) _particle->_index = 28; // when SCATTERED in GROUP: reach exit of GROUP after Blade_13
      else ABSORBED=0; // not SCATTERED within GROUP: always tries next
      _particle->_index++;
    } /* end component Blade_11 [26] */
    if (!ABSORBED && _particle->_index == 27) {
      /* component Blade_12=Monochromator_curved() [27] */
  #define reflect (_Blade_12->_parameters.reflect)
  #define transmit (_Blade_12->_parameters.transmit)
  #define zwidth (_Blade_12->_parameters.zwidth)
  #define yheight (_Blade_12->_parameters.yheight)
  #define gap (_Blade_12->_parameters.gap)
  #define NH (_Blade_12->_parameters.NH)
  #define NV (_Blade_12->_parameters.NV)
  #define mosaich (_Blade_12->_parameters.mosaich)
  #define mosaicv (_Blade_12->_parameters.mosaicv)
  #define r0 (_Blade_12->_parameters.r0)
  #define t0 (_Blade_12->_parameters.t0)
  #define Q (_Blade_12->_parameters.Q)
  #define RV (_Blade_12->_parameters.RV)
  #define RH (_Blade_12->_parameters.RH)
  #define DM (_Blade_12->_parameters.DM)
  #define mosaic (_Blade_12->_parameters.mosaic)
  #define width (_Blade_12->_parameters.width)
  #define height (_Blade_12->_parameters.height)
  #define verbose (_Blade_12->_parameters.verbose)
  #define order (_Blade_12->_parameters.order)
  #define mos_rms_y (_Blade_12->_parameters.mos_rms_y)
  #define mos_rms_z (_Blade_12->_parameters.mos_rms_z)
  #define mos_rms_max (_Blade_12->_parameters.mos_rms_max)
  #define mono_Q (_Blade_12->_parameters.mono_Q)
  #define SlabWidth (_Blade_12->_parameters.SlabWidth)
  #define SlabHeight (_Blade_12->_parameters.SlabHeight)
  #define rTable (_Blade_12->_parameters.rTable)
  #define tTable (_Blade_12->_parameters.tTable)
  #define row (_Blade_12->_parameters.row)
  #define col (_Blade_12->_parameters.col)
  #define tiltH (_Blade_12->_parameters.tiltH)
  #define tiltV (_Blade_12->_parameters.tiltV)
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Blade_12->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Blade_12->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Blade_12->_position_relative, _Blade_12->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      if (( instrument->_parameters._full_instrument == 1 )) // conditional WHEN execution
      class_Monochromator_curved_trace(_Blade_12, _particle);
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
      // GROUP Monochro: from Blade_1 [16] to Blade_13 [28]
      if (SCATTERED) _particle->_index = 28; // when SCATTERED in GROUP: reach exit of GROUP after Blade_13
      else ABSORBED=0; // not SCATTERED within GROUP: always tries next
      _particle->_index++;
    } /* end component Blade_12 [27] */
    if (!ABSORBED && _particle->_index == 28) {
      /* component Blade_13=Monochromator_curved() [28] */
  #define reflect (_Blade_13->_parameters.reflect)
  #define transmit (_Blade_13->_parameters.transmit)
  #define zwidth (_Blade_13->_parameters.zwidth)
  #define yheight (_Blade_13->_parameters.yheight)
  #define gap (_Blade_13->_parameters.gap)
  #define NH (_Blade_13->_parameters.NH)
  #define NV (_Blade_13->_parameters.NV)
  #define mosaich (_Blade_13->_parameters.mosaich)
  #define mosaicv (_Blade_13->_parameters.mosaicv)
  #define r0 (_Blade_13->_parameters.r0)
  #define t0 (_Blade_13->_parameters.t0)
  #define Q (_Blade_13->_parameters.Q)
  #define RV (_Blade_13->_parameters.RV)
  #define RH (_Blade_13->_parameters.RH)
  #define DM (_Blade_13->_parameters.DM)
  #define mosaic (_Blade_13->_parameters.mosaic)
  #define width (_Blade_13->_parameters.width)
  #define height (_Blade_13->_parameters.height)
  #define verbose (_Blade_13->_parameters.verbose)
  #define order (_Blade_13->_parameters.order)
  #define mos_rms_y (_Blade_13->_parameters.mos_rms_y)
  #define mos_rms_z (_Blade_13->_parameters.mos_rms_z)
  #define mos_rms_max (_Blade_13->_parameters.mos_rms_max)
  #define mono_Q (_Blade_13->_parameters.mono_Q)
  #define SlabWidth (_Blade_13->_parameters.SlabWidth)
  #define SlabHeight (_Blade_13->_parameters.SlabHeight)
  #define rTable (_Blade_13->_parameters.rTable)
  #define tTable (_Blade_13->_parameters.tTable)
  #define row (_Blade_13->_parameters.row)
  #define col (_Blade_13->_parameters.col)
  #define tiltH (_Blade_13->_parameters.tiltH)
  #define tiltV (_Blade_13->_parameters.tiltV)
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Blade_13->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Blade_13->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Blade_13->_position_relative, _Blade_13->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      if (( instrument->_parameters._full_instrument == 1 )) // conditional WHEN execution
      class_Monochromator_curved_trace(_Blade_13, _particle);
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
      // GROUP Monochro: from Blade_1 [16] to Blade_13 [28]
      if (SCATTERED) _particle->_index = 28; // when SCATTERED in GROUP: reach exit of GROUP after Blade_13
      else ABSORB;     // not SCATTERED at end of GROUP: removes left events
      _particle->_index++;
    } /* end component Blade_13 [28] */
    if (!ABSORBED && _particle->_index == 29) {
      /* component PSD_At_sec_shutter=PSD_monitor() [29] */
  #define nx (_PSD_At_sec_shutter->_parameters.nx)
  #define ny (_PSD_At_sec_shutter->_parameters.ny)
  #define filename (_PSD_At_sec_shutter->_parameters.filename)
  #define xmin (_PSD_At_sec_shutter->_parameters.xmin)
  #define xmax (_PSD_At_sec_shutter->_parameters.xmax)
  #define ymin (_PSD_At_sec_shutter->_parameters.ymin)
  #define ymax (_PSD_At_sec_shutter->_parameters.ymax)
  #define xwidth (_PSD_At_sec_shutter->_parameters.xwidth)
  #define yheight (_PSD_At_sec_shutter->_parameters.yheight)
  #define restore_neutron (_PSD_At_sec_shutter->_parameters.restore_neutron)
  #define PSD_N (_PSD_At_sec_shutter->_parameters.PSD_N)
  #define PSD_p (_PSD_At_sec_shutter->_parameters.PSD_p)
  #define PSD_p2 (_PSD_At_sec_shutter->_parameters.PSD_p2)
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_PSD_At_sec_shutter->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _PSD_At_sec_shutter->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_PSD_At_sec_shutter->_position_relative, _PSD_At_sec_shutter->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      if (( instrument->_parameters._full_instrument == 1 )) // conditional WHEN execution
      class_PSD_monitor_trace(_PSD_At_sec_shutter, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
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
      _particle->_index++;
    } /* end component PSD_At_sec_shutter [29] */
    if (!ABSORBED && _particle->_index == 30) {
      /* component LAM_At_sec_shutter=L_monitor() [30] */
  #define nL (_LAM_At_sec_shutter->_parameters.nL)
  #define filename (_LAM_At_sec_shutter->_parameters.filename)
  #define xmin (_LAM_At_sec_shutter->_parameters.xmin)
  #define xmax (_LAM_At_sec_shutter->_parameters.xmax)
  #define ymin (_LAM_At_sec_shutter->_parameters.ymin)
  #define ymax (_LAM_At_sec_shutter->_parameters.ymax)
  #define xwidth (_LAM_At_sec_shutter->_parameters.xwidth)
  #define yheight (_LAM_At_sec_shutter->_parameters.yheight)
  #define Lmin (_LAM_At_sec_shutter->_parameters.Lmin)
  #define Lmax (_LAM_At_sec_shutter->_parameters.Lmax)
  #define restore_neutron (_LAM_At_sec_shutter->_parameters.restore_neutron)
  #define L_N (_LAM_At_sec_shutter->_parameters.L_N)
  #define L_p (_LAM_At_sec_shutter->_parameters.L_p)
  #define L_p2 (_LAM_At_sec_shutter->_parameters.L_p2)
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_LAM_At_sec_shutter->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _LAM_At_sec_shutter->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_LAM_At_sec_shutter->_position_relative, _LAM_At_sec_shutter->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      if (( instrument->_parameters._full_instrument == 1 )) // conditional WHEN execution
      class_L_monitor_trace(_LAM_At_sec_shutter, _particle);
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
    } /* end component LAM_At_sec_shutter [30] */
    if (!ABSORBED && _particle->_index == 31) {
      /* component Inside_chamber_collimator=Slit() [31] */
  #define xmin (_Inside_chamber_collimator->_parameters.xmin)
  #define xmax (_Inside_chamber_collimator->_parameters.xmax)
  #define ymin (_Inside_chamber_collimator->_parameters.ymin)
  #define ymax (_Inside_chamber_collimator->_parameters.ymax)
  #define radius (_Inside_chamber_collimator->_parameters.radius)
  #define xwidth (_Inside_chamber_collimator->_parameters.xwidth)
  #define yheight (_Inside_chamber_collimator->_parameters.yheight)
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Inside_chamber_collimator->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Inside_chamber_collimator->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Inside_chamber_collimator->_position_relative, _Inside_chamber_collimator->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      if (( instrument->_parameters._full_instrument == 1 )) // conditional WHEN execution
      class_Slit_trace(_Inside_chamber_collimator, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef radius
  #undef xwidth
  #undef yheight
      _particle->_index++;
    } /* end component Inside_chamber_collimator [31] */
    if (!ABSORBED && _particle->_index == 32) {
      /* component PSD_After_inside_chamber_collimator=PSD_monitor() [32] */
  #define nx (_PSD_After_inside_chamber_collimator->_parameters.nx)
  #define ny (_PSD_After_inside_chamber_collimator->_parameters.ny)
  #define filename (_PSD_After_inside_chamber_collimator->_parameters.filename)
  #define xmin (_PSD_After_inside_chamber_collimator->_parameters.xmin)
  #define xmax (_PSD_After_inside_chamber_collimator->_parameters.xmax)
  #define ymin (_PSD_After_inside_chamber_collimator->_parameters.ymin)
  #define ymax (_PSD_After_inside_chamber_collimator->_parameters.ymax)
  #define xwidth (_PSD_After_inside_chamber_collimator->_parameters.xwidth)
  #define yheight (_PSD_After_inside_chamber_collimator->_parameters.yheight)
  #define restore_neutron (_PSD_After_inside_chamber_collimator->_parameters.restore_neutron)
  #define PSD_N (_PSD_After_inside_chamber_collimator->_parameters.PSD_N)
  #define PSD_p (_PSD_After_inside_chamber_collimator->_parameters.PSD_p)
  #define PSD_p2 (_PSD_After_inside_chamber_collimator->_parameters.PSD_p2)
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_PSD_After_inside_chamber_collimator->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _PSD_After_inside_chamber_collimator->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_PSD_After_inside_chamber_collimator->_position_relative, _PSD_After_inside_chamber_collimator->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      if (( instrument->_parameters._full_instrument == 1 )) // conditional WHEN execution
      class_PSD_monitor_trace(_PSD_After_inside_chamber_collimator, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
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
      _particle->_index++;
    } /* end component PSD_After_inside_chamber_collimator [32] */
    if (!ABSORBED && _particle->_index == 33) {
      /* component Outside_chamber_collimator=Slit() [33] */
  #define xmin (_Outside_chamber_collimator->_parameters.xmin)
  #define xmax (_Outside_chamber_collimator->_parameters.xmax)
  #define ymin (_Outside_chamber_collimator->_parameters.ymin)
  #define ymax (_Outside_chamber_collimator->_parameters.ymax)
  #define radius (_Outside_chamber_collimator->_parameters.radius)
  #define xwidth (_Outside_chamber_collimator->_parameters.xwidth)
  #define yheight (_Outside_chamber_collimator->_parameters.yheight)
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Outside_chamber_collimator->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Outside_chamber_collimator->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Outside_chamber_collimator->_position_relative, _Outside_chamber_collimator->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      if (( instrument->_parameters._full_instrument == 1 )) // conditional WHEN execution
      class_Slit_trace(_Outside_chamber_collimator, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef radius
  #undef xwidth
  #undef yheight
      _particle->_index++;
    } /* end component Outside_chamber_collimator [33] */
    if (!ABSORBED && _particle->_index == 34) {
      /* component PSD_Outside_chamber_collimator_1=PSD_monitor() [34] */
  #define nx (_PSD_Outside_chamber_collimator_1->_parameters.nx)
  #define ny (_PSD_Outside_chamber_collimator_1->_parameters.ny)
  #define filename (_PSD_Outside_chamber_collimator_1->_parameters.filename)
  #define xmin (_PSD_Outside_chamber_collimator_1->_parameters.xmin)
  #define xmax (_PSD_Outside_chamber_collimator_1->_parameters.xmax)
  #define ymin (_PSD_Outside_chamber_collimator_1->_parameters.ymin)
  #define ymax (_PSD_Outside_chamber_collimator_1->_parameters.ymax)
  #define xwidth (_PSD_Outside_chamber_collimator_1->_parameters.xwidth)
  #define yheight (_PSD_Outside_chamber_collimator_1->_parameters.yheight)
  #define restore_neutron (_PSD_Outside_chamber_collimator_1->_parameters.restore_neutron)
  #define PSD_N (_PSD_Outside_chamber_collimator_1->_parameters.PSD_N)
  #define PSD_p (_PSD_Outside_chamber_collimator_1->_parameters.PSD_p)
  #define PSD_p2 (_PSD_Outside_chamber_collimator_1->_parameters.PSD_p2)
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_PSD_Outside_chamber_collimator_1->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _PSD_Outside_chamber_collimator_1->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_PSD_Outside_chamber_collimator_1->_position_relative, _PSD_Outside_chamber_collimator_1->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      if (( instrument->_parameters._full_instrument == 1 )) // conditional WHEN execution
      class_PSD_monitor_trace(_PSD_Outside_chamber_collimator_1, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
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
      _particle->_index++;
    } /* end component PSD_Outside_chamber_collimator_1 [34] */
    if (!ABSORBED && _particle->_index == 35) {
      /* component PSD_Outside_chamber_collimator=PSD_monitor() [35] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_PSD_Outside_chamber_collimator->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _PSD_Outside_chamber_collimator->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_PSD_Outside_chamber_collimator->_position_relative, _PSD_Outside_chamber_collimator->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_PSD_monitor_trace(_PSD_Outside_chamber_collimator, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component PSD_Outside_chamber_collimator [35] */
    if (!ABSORBED && _particle->_index == 36) {
      /* component Incident_slit_h=Slit() [36] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Incident_slit_h->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Incident_slit_h->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Incident_slit_h->_position_relative, _Incident_slit_h->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Slit_trace(_Incident_slit_h, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component Incident_slit_h [36] */
    if (!ABSORBED && _particle->_index == 37) {
      /* component Incident_slit_w=Slit() [37] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Incident_slit_w->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Incident_slit_w->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Incident_slit_w->_position_relative, _Incident_slit_w->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Slit_trace(_Incident_slit_w, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component Incident_slit_w [37] */
    if (!ABSORBED && _particle->_index == 38) {
      /* component PSD_After_Incident_slit_w=PSD_monitor() [38] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_PSD_After_Incident_slit_w->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _PSD_After_Incident_slit_w->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_PSD_After_Incident_slit_w->_position_relative, _PSD_After_Incident_slit_w->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_PSD_monitor_trace(_PSD_After_Incident_slit_w, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component PSD_After_Incident_slit_w [38] */
    if (!ABSORBED && _particle->_index == 39) {
      /* component Center_of_rotation=Arm() [39] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Center_of_rotation->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Center_of_rotation->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Center_of_rotation->_position_relative, _Center_of_rotation->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      _particle->_index++;
    } /* end component Center_of_rotation [39] */
    if (!ABSORBED && _particle->_index == 40) {
      /* component Sample_rotation=Arm() [40] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Sample_rotation->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Sample_rotation->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Sample_rotation->_position_relative, _Sample_rotation->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      _particle->_index++;
    } /* end component Sample_rotation [40] */
    if (!ABSORBED && _particle->_index == 41) {
      /* component Sample_location=Arm() [41] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Sample_location->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Sample_location->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Sample_location->_position_relative, _Sample_location->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      _particle->_index++;
    } /* end component Sample_location [41] */
    if (!ABSORBED && _particle->_index == 42) {
      /* component PSD_Center_of_rotation=PSD_monitor() [42] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_PSD_Center_of_rotation->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _PSD_Center_of_rotation->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_PSD_Center_of_rotation->_position_relative, _PSD_Center_of_rotation->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_PSD_monitor_trace(_PSD_Center_of_rotation, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component PSD_Center_of_rotation [42] */
    if (!ABSORBED && _particle->_index == 43) {
      /* component DIV_Center_of_rotation=Divergence_monitor() [43] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_DIV_Center_of_rotation->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _DIV_Center_of_rotation->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_DIV_Center_of_rotation->_position_relative, _DIV_Center_of_rotation->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Divergence_monitor_trace(_DIV_Center_of_rotation, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component DIV_Center_of_rotation [43] */
    /* start SPLIT at Sample */
    instrument->logic.Split_Sample_particle=_particle;
    if (!ABSORBED && _particle->_index == 44) 
    for (instrument->logic.Split_Sample = 0; instrument->logic.Split_Sample< 10; instrument->logic.Split_Sample++) {
    {
      /* component Sample=PowderN() [44] */
      _particle=instrument->logic.Split_Sample_particle;
      p /= 10 > 0 ? 10 : 1;
  #define reflections (_Sample->_parameters.reflections)
  #define geometry (_Sample->_parameters.geometry)
  #define format (_Sample->_parameters.format)
  #define radius (_Sample->_parameters.radius)
  #define yheight (_Sample->_parameters.yheight)
  #define xwidth (_Sample->_parameters.xwidth)
  #define zdepth (_Sample->_parameters.zdepth)
  #define thickness (_Sample->_parameters.thickness)
  #define pack (_Sample->_parameters.pack)
  #define Vc (_Sample->_parameters.Vc)
  #define sigma_abs (_Sample->_parameters.sigma_abs)
  #define sigma_inc (_Sample->_parameters.sigma_inc)
  #define delta_d_d (_Sample->_parameters.delta_d_d)
  #define p_inc (_Sample->_parameters.p_inc)
  #define p_transmit (_Sample->_parameters.p_transmit)
  #define DW (_Sample->_parameters.DW)
  #define nb_atoms (_Sample->_parameters.nb_atoms)
  #define d_phi (_Sample->_parameters.d_phi)
  #define p_interact (_Sample->_parameters.p_interact)
  #define concentric (_Sample->_parameters.concentric)
  #define density (_Sample->_parameters.density)
  #define weight (_Sample->_parameters.weight)
  #define barns (_Sample->_parameters.barns)
  #define Strain (_Sample->_parameters.Strain)
  #define focus_flip (_Sample->_parameters.focus_flip)
  #define line_info (_Sample->_parameters.line_info)
  #define columns (_Sample->_parameters.columns)
  #define offdata (_Sample->_parameters.offdata)
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Sample->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Sample->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Sample->_position_relative, _Sample->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_PowderN_trace(_Sample, _particle); /* contains EXTEND code */
      if (_particle->_restore)
        *_particle = _particle_save;
  #undef reflections
  #undef geometry
  #undef format
  #undef radius
  #undef yheight
  #undef xwidth
  #undef zdepth
  #undef thickness
  #undef pack
  #undef Vc
  #undef sigma_abs
  #undef sigma_inc
  #undef delta_d_d
  #undef p_inc
  #undef p_transmit
  #undef DW
  #undef nb_atoms
  #undef d_phi
  #undef p_interact
  #undef concentric
  #undef density
  #undef weight
  #undef barns
  #undef Strain
  #undef focus_flip
  #undef line_info
  #undef columns
  #undef offdata
      _particle->_index++;
    } /* end component Sample [44] */
    if (!ABSORBED && _particle->_index == 45) {
      /* component Det_axis=Arm() [45] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Det_axis->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Det_axis->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Det_axis->_position_relative, _Det_axis->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      _particle->_index++;
    } /* end component Det_axis [45] */
    if (!ABSORBED && _particle->_index == 46) {
      /* component Det_axis_2=Arm() [46] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Det_axis_2->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Det_axis_2->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Det_axis_2->_position_relative, _Det_axis_2->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      _particle->_index++;
    } /* end component Det_axis_2 [46] */
    if (!ABSORBED && _particle->_index == 47) {
      /* component Det_axis_3=Arm() [47] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Det_axis_3->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Det_axis_3->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Det_axis_3->_position_relative, _Det_axis_3->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      _particle->_index++;
    } /* end component Det_axis_3 [47] */
    if (!ABSORBED && _particle->_index == 48) {
      /* component Det_axis_4=Arm() [48] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Det_axis_4->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Det_axis_4->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Det_axis_4->_position_relative, _Det_axis_4->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      _particle->_index++;
    } /* end component Det_axis_4 [48] */
    if (!ABSORBED && _particle->_index == 49) {
      /* component RadColl=Collimator_radial() [49] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_RadColl->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _RadColl->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_RadColl->_position_relative, _RadColl->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Collimator_radial_trace(_RadColl, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component RadColl [49] */
    if (!ABSORBED && _particle->_index == 50) {
      /* component PSD_Detector=PSD_monitor() [50] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_PSD_Detector->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _PSD_Detector->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_PSD_Detector->_position_relative, _PSD_Detector->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_PSD_monitor_trace(_PSD_Detector, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      // GROUP Detectors: from PSD_Detector [50] to PSD_Detector_4 [53]
      if (SCATTERED) _particle->_index = 53; // when SCATTERED in GROUP: reach exit of GROUP after PSD_Detector_4
      else ABSORBED=0; // not SCATTERED within GROUP: always tries next
      _particle->_index++;
    } /* end component PSD_Detector [50] */
    if (!ABSORBED && _particle->_index == 51) {
      /* component PSD_Detector_2=PSD_monitor() [51] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_PSD_Detector_2->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _PSD_Detector_2->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_PSD_Detector_2->_position_relative, _PSD_Detector_2->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_PSD_monitor_trace(_PSD_Detector_2, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      // GROUP Detectors: from PSD_Detector [50] to PSD_Detector_4 [53]
      if (SCATTERED) _particle->_index = 53; // when SCATTERED in GROUP: reach exit of GROUP after PSD_Detector_4
      else ABSORBED=0; // not SCATTERED within GROUP: always tries next
      _particle->_index++;
    } /* end component PSD_Detector_2 [51] */
    if (!ABSORBED && _particle->_index == 52) {
      /* component PSD_Detector_3=PSD_monitor() [52] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_PSD_Detector_3->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _PSD_Detector_3->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_PSD_Detector_3->_position_relative, _PSD_Detector_3->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_PSD_monitor_trace(_PSD_Detector_3, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      // GROUP Detectors: from PSD_Detector [50] to PSD_Detector_4 [53]
      if (SCATTERED) _particle->_index = 53; // when SCATTERED in GROUP: reach exit of GROUP after PSD_Detector_4
      else ABSORBED=0; // not SCATTERED within GROUP: always tries next
      _particle->_index++;
    } /* end component PSD_Detector_3 [52] */
    if (!ABSORBED && _particle->_index == 53) {
      /* component PSD_Detector_4=PSD_monitor() [53] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_PSD_Detector_4->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _PSD_Detector_4->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_PSD_Detector_4->_position_relative, _PSD_Detector_4->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_PSD_monitor_trace(_PSD_Detector_4, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      // GROUP Detectors: from PSD_Detector [50] to PSD_Detector_4 [53]
      if (SCATTERED) _particle->_index = 53; // when SCATTERED in GROUP: reach exit of GROUP after PSD_Detector_4
      else ABSORB;     // not SCATTERED at end of GROUP: removes left events
      _particle->_index++;
    } /* end component PSD_Detector_4 [53] */
    } /* end SPLIT at Sample */
    if (_particle->_index > 53)
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
* instrument 'SAFARI_PITSI' and components SAVE
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

_class_Divergence_monitor *class_Divergence_monitor_save(_class_Divergence_monitor *_comp
) {
  #define nh (_comp->_parameters.nh)
  #define nv (_comp->_parameters.nv)
  #define filename (_comp->_parameters.filename)
  #define xmin (_comp->_parameters.xmin)
  #define xmax (_comp->_parameters.xmax)
  #define ymin (_comp->_parameters.ymin)
  #define ymax (_comp->_parameters.ymax)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define maxdiv_h (_comp->_parameters.maxdiv_h)
  #define maxdiv_v (_comp->_parameters.maxdiv_v)
  #define restore_neutron (_comp->_parameters.restore_neutron)
  #define nx (_comp->_parameters.nx)
  #define ny (_comp->_parameters.ny)
  #define nz (_comp->_parameters.nz)
  #define Div_N (_comp->_parameters.Div_N)
  #define Div_p (_comp->_parameters.Div_p)
  #define Div_p2 (_comp->_parameters.Div_p2)
  DETECTOR_OUT_2D(
      "Divergence monitor",
      "X divergence [deg]",
      "Y divergence [deg]",
      -maxdiv_h, maxdiv_h, -maxdiv_v, maxdiv_v,
      nh, nv,
      &Div_N[0][0],&Div_p[0][0],&Div_p2[0][0],
      filename);
  #undef nh
  #undef nv
  #undef filename
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef xwidth
  #undef yheight
  #undef maxdiv_h
  #undef maxdiv_v
  #undef restore_neutron
  #undef nx
  #undef ny
  #undef nz
  #undef Div_N
  #undef Div_p
  #undef Div_p2
  return(_comp);
} /* class_Divergence_monitor_save */



int save(FILE *handle) { /* called by mccode_main for SAFARI_PITSI:SAVE */
  if (!handle) siminfo_init(NULL);

  /* call iteratively all components SAVE */
  class_Progress_bar_save(_Progress);




  class_PSD_monitor_save(_PSD_Source);

  class_L_monitor_save(_LAM_Source);




  class_L_monitor_save(_LAM_After_sapphire);

  class_PSD_monitor_save(_PSD_After_sapphire);



  class_PSD_monitor_save(_PSD_After_Outlet);















  class_PSD_monitor_save(_PSD_At_sec_shutter);

  class_L_monitor_save(_LAM_At_sec_shutter);


  class_PSD_monitor_save(_PSD_After_inside_chamber_collimator);


  class_PSD_monitor_save(_PSD_Outside_chamber_collimator_1);

  class_PSD_monitor_save(_PSD_Outside_chamber_collimator);



  class_PSD_monitor_save(_PSD_After_Incident_slit_w);




  class_PSD_monitor_save(_PSD_Center_of_rotation);

  class_Divergence_monitor_save(_DIV_Center_of_rotation);







  class_PSD_monitor_save(_PSD_Detector);

  class_PSD_monitor_save(_PSD_Detector_2);

  class_PSD_monitor_save(_PSD_Detector_3);

  class_PSD_monitor_save(_PSD_Detector_4);

  if (!handle) siminfo_close(); 

  return(0);
} /* save */

/* *****************************************************************************
* instrument 'SAFARI_PITSI' and components FINALLY
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

_class_Filter_gen *class_Filter_gen_finally(_class_Filter_gen *_comp
) {
  #define filename (_comp->_parameters.filename)
  #define options (_comp->_parameters.options)
  #define xmin (_comp->_parameters.xmin)
  #define xmax (_comp->_parameters.xmax)
  #define ymin (_comp->_parameters.ymin)
  #define ymax (_comp->_parameters.ymax)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define thickness (_comp->_parameters.thickness)
  #define scaling (_comp->_parameters.scaling)
  #define verbose (_comp->_parameters.verbose)
  #define pTable (_comp->_parameters.pTable)
  #define Mode_Table (_comp->_parameters.Mode_Table)
  #define Type_Table (_comp->_parameters.Type_Table)
  Table_Free(&pTable);
  #undef filename
  #undef options
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef xwidth
  #undef yheight
  #undef thickness
  #undef scaling
  #undef verbose
  #undef pTable
  #undef Mode_Table
  #undef Type_Table
  return(_comp);
} /* class_Filter_gen_finally */

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

_class_Divergence_monitor *class_Divergence_monitor_finally(_class_Divergence_monitor *_comp
) {
  #define nh (_comp->_parameters.nh)
  #define nv (_comp->_parameters.nv)
  #define filename (_comp->_parameters.filename)
  #define xmin (_comp->_parameters.xmin)
  #define xmax (_comp->_parameters.xmax)
  #define ymin (_comp->_parameters.ymin)
  #define ymax (_comp->_parameters.ymax)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define maxdiv_h (_comp->_parameters.maxdiv_h)
  #define maxdiv_v (_comp->_parameters.maxdiv_v)
  #define restore_neutron (_comp->_parameters.restore_neutron)
  #define nx (_comp->_parameters.nx)
  #define ny (_comp->_parameters.ny)
  #define nz (_comp->_parameters.nz)
  #define Div_N (_comp->_parameters.Div_N)
  #define Div_p (_comp->_parameters.Div_p)
  #define Div_p2 (_comp->_parameters.Div_p2)
  destroy_darr2d(Div_N);
  destroy_darr2d(Div_p);
  destroy_darr2d(Div_p2);
  #undef nh
  #undef nv
  #undef filename
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef xwidth
  #undef yheight
  #undef maxdiv_h
  #undef maxdiv_v
  #undef restore_neutron
  #undef nx
  #undef ny
  #undef nz
  #undef Div_N
  #undef Div_p
  #undef Div_p2
  return(_comp);
} /* class_Divergence_monitor_finally */

_class_PowderN *class_PowderN_finally(_class_PowderN *_comp
) {
  #define reflections (_comp->_parameters.reflections)
  #define geometry (_comp->_parameters.geometry)
  #define format (_comp->_parameters.format)
  #define radius (_comp->_parameters.radius)
  #define yheight (_comp->_parameters.yheight)
  #define xwidth (_comp->_parameters.xwidth)
  #define zdepth (_comp->_parameters.zdepth)
  #define thickness (_comp->_parameters.thickness)
  #define pack (_comp->_parameters.pack)
  #define Vc (_comp->_parameters.Vc)
  #define sigma_abs (_comp->_parameters.sigma_abs)
  #define sigma_inc (_comp->_parameters.sigma_inc)
  #define delta_d_d (_comp->_parameters.delta_d_d)
  #define p_inc (_comp->_parameters.p_inc)
  #define p_transmit (_comp->_parameters.p_transmit)
  #define DW (_comp->_parameters.DW)
  #define nb_atoms (_comp->_parameters.nb_atoms)
  #define d_phi (_comp->_parameters.d_phi)
  #define p_interact (_comp->_parameters.p_interact)
  #define concentric (_comp->_parameters.concentric)
  #define density (_comp->_parameters.density)
  #define weight (_comp->_parameters.weight)
  #define barns (_comp->_parameters.barns)
  #define Strain (_comp->_parameters.Strain)
  #define focus_flip (_comp->_parameters.focus_flip)
  #define line_info (_comp->_parameters.line_info)
  #define columns (_comp->_parameters.columns)
  #define offdata (_comp->_parameters.offdata)
  free(line_info.list);
  free(line_info.q_v);
  free(line_info.w_v);
  free(line_info.my_s_v2);
  MPI_MASTER(
  if (line_info.flag_warning)
    printf("PowderN: %s: Error messages were repeated %i times with absorbed neutrons.\n",
      NAME_CURRENT_COMP, line_info.flag_warning);

  /* in case this instance is used in a SPLIT, we can recommend the
     optimal iteration value */
  if (line_info.nb_refl_count) {
    double split_iterations = (double)line_info.nb_reuses/line_info.nb_refl_count + 1;
    double split_optimal    = (double)line_info.nb_refl/line_info.nb_refl_count;
    if (split_optimal > split_iterations + 5)
      printf("PowderN: %s: Info: you may highly improve the computation efficiency by using\n"
        "    SPLIT %i COMPONENT %s=PowderN(...)\n"
        "  in the instrument description %s.\n",
        NAME_CURRENT_COMP, (int)split_optimal, NAME_CURRENT_COMP, instrument_source);
  }
  );

  #undef reflections
  #undef geometry
  #undef format
  #undef radius
  #undef yheight
  #undef xwidth
  #undef zdepth
  #undef thickness
  #undef pack
  #undef Vc
  #undef sigma_abs
  #undef sigma_inc
  #undef delta_d_d
  #undef p_inc
  #undef p_transmit
  #undef DW
  #undef nb_atoms
  #undef d_phi
  #undef p_interact
  #undef concentric
  #undef density
  #undef weight
  #undef barns
  #undef Strain
  #undef focus_flip
  #undef line_info
  #undef columns
  #undef offdata
  return(_comp);
} /* class_PowderN_finally */



int finally(void) { /* called by mccode_main for SAFARI_PITSI:FINALLY */
  siminfo_init(NULL);
  save(siminfo_file); /* save data when simulation ends */

  /* call iteratively all components FINALLY */
  class_Progress_bar_finally(_Progress);



  class_Source_gen_finally(_Source);

  class_PSD_monitor_finally(_PSD_Source);

  class_L_monitor_finally(_LAM_Source);


  class_Filter_gen_finally(_Sapphire_filter);


  class_L_monitor_finally(_LAM_After_sapphire);

  class_PSD_monitor_finally(_PSD_After_sapphire);



  class_PSD_monitor_finally(_PSD_After_Outlet);


  class_Monochromator_curved_finally(_Blade_1);

  class_Monochromator_curved_finally(_Blade_2);

  class_Monochromator_curved_finally(_Blade_3);

  class_Monochromator_curved_finally(_Blade_4);

  class_Monochromator_curved_finally(_Blade_5);

  class_Monochromator_curved_finally(_Blade_6);

  class_Monochromator_curved_finally(_Blade_7);

  class_Monochromator_curved_finally(_Blade_8);

  class_Monochromator_curved_finally(_Blade_9);

  class_Monochromator_curved_finally(_Blade_10);

  class_Monochromator_curved_finally(_Blade_11);

  class_Monochromator_curved_finally(_Blade_12);

  class_Monochromator_curved_finally(_Blade_13);

  class_PSD_monitor_finally(_PSD_At_sec_shutter);

  class_L_monitor_finally(_LAM_At_sec_shutter);


  class_PSD_monitor_finally(_PSD_After_inside_chamber_collimator);


  class_PSD_monitor_finally(_PSD_Outside_chamber_collimator_1);

  class_PSD_monitor_finally(_PSD_Outside_chamber_collimator);



  class_PSD_monitor_finally(_PSD_After_Incident_slit_w);




  class_PSD_monitor_finally(_PSD_Center_of_rotation);

  class_Divergence_monitor_finally(_DIV_Center_of_rotation);

  class_PowderN_finally(_Sample);






  class_PSD_monitor_finally(_PSD_Detector);

  class_PSD_monitor_finally(_PSD_Detector_2);

  class_PSD_monitor_finally(_PSD_Detector_3);

  class_PSD_monitor_finally(_PSD_Detector_4);

  siminfo_close(); 

  return(0);
} /* finally */

/* *****************************************************************************
* instrument 'SAFARI_PITSI' and components DISPLAY
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

_class_Filter_gen *class_Filter_gen_display(_class_Filter_gen *_comp
) {
  #define filename (_comp->_parameters.filename)
  #define options (_comp->_parameters.options)
  #define xmin (_comp->_parameters.xmin)
  #define xmax (_comp->_parameters.xmax)
  #define ymin (_comp->_parameters.ymin)
  #define ymax (_comp->_parameters.ymax)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define thickness (_comp->_parameters.thickness)
  #define scaling (_comp->_parameters.scaling)
  #define verbose (_comp->_parameters.verbose)
  #define pTable (_comp->_parameters.pTable)
  #define Mode_Table (_comp->_parameters.Mode_Table)
  #define Type_Table (_comp->_parameters.Type_Table)
  
  multiline(5, (double)xmin, (double)ymin, 0.0,
               (double)xmax, (double)ymin, 0.0,
               (double)xmax, (double)ymax, 0.0,
               (double)xmin, (double)ymax, 0.0,
               (double)xmin, (double)ymin, 0.0);
  #undef filename
  #undef options
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef xwidth
  #undef yheight
  #undef thickness
  #undef scaling
  #undef verbose
  #undef pTable
  #undef Mode_Table
  #undef Type_Table
  return(_comp);
} /* class_Filter_gen_display */

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

_class_Divergence_monitor *class_Divergence_monitor_display(_class_Divergence_monitor *_comp
) {
  #define nh (_comp->_parameters.nh)
  #define nv (_comp->_parameters.nv)
  #define filename (_comp->_parameters.filename)
  #define xmin (_comp->_parameters.xmin)
  #define xmax (_comp->_parameters.xmax)
  #define ymin (_comp->_parameters.ymin)
  #define ymax (_comp->_parameters.ymax)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define maxdiv_h (_comp->_parameters.maxdiv_h)
  #define maxdiv_v (_comp->_parameters.maxdiv_v)
  #define restore_neutron (_comp->_parameters.restore_neutron)
  #define nx (_comp->_parameters.nx)
  #define ny (_comp->_parameters.ny)
  #define nz (_comp->_parameters.nz)
  #define Div_N (_comp->_parameters.Div_N)
  #define Div_p (_comp->_parameters.Div_p)
  #define Div_p2 (_comp->_parameters.Div_p2)
  multiline(5, (double)xmin, (double)ymin, 0.0,
               (double)xmax, (double)ymin, 0.0,
               (double)xmax, (double)ymax, 0.0,
               (double)xmin, (double)ymax, 0.0,
               (double)xmin, (double)ymin, 0.0);
  #undef nh
  #undef nv
  #undef filename
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef xwidth
  #undef yheight
  #undef maxdiv_h
  #undef maxdiv_v
  #undef restore_neutron
  #undef nx
  #undef ny
  #undef nz
  #undef Div_N
  #undef Div_p
  #undef Div_p2
  return(_comp);
} /* class_Divergence_monitor_display */

_class_PowderN *class_PowderN_display(_class_PowderN *_comp
) {
  #define reflections (_comp->_parameters.reflections)
  #define geometry (_comp->_parameters.geometry)
  #define format (_comp->_parameters.format)
  #define radius (_comp->_parameters.radius)
  #define yheight (_comp->_parameters.yheight)
  #define xwidth (_comp->_parameters.xwidth)
  #define zdepth (_comp->_parameters.zdepth)
  #define thickness (_comp->_parameters.thickness)
  #define pack (_comp->_parameters.pack)
  #define Vc (_comp->_parameters.Vc)
  #define sigma_abs (_comp->_parameters.sigma_abs)
  #define sigma_inc (_comp->_parameters.sigma_inc)
  #define delta_d_d (_comp->_parameters.delta_d_d)
  #define p_inc (_comp->_parameters.p_inc)
  #define p_transmit (_comp->_parameters.p_transmit)
  #define DW (_comp->_parameters.DW)
  #define nb_atoms (_comp->_parameters.nb_atoms)
  #define d_phi (_comp->_parameters.d_phi)
  #define p_interact (_comp->_parameters.p_interact)
  #define concentric (_comp->_parameters.concentric)
  #define density (_comp->_parameters.density)
  #define weight (_comp->_parameters.weight)
  #define barns (_comp->_parameters.barns)
  #define Strain (_comp->_parameters.Strain)
  #define focus_flip (_comp->_parameters.focus_flip)
  #define line_info (_comp->_parameters.line_info)
  #define columns (_comp->_parameters.columns)
  #define offdata (_comp->_parameters.offdata)
  if (line_info.V_0) {

    if (line_info.shape == 0) { /* cyl */
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
    } else if (line_info.shape == 1) {  /* box */
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
      if (line_info.zdepth_i) {
        xmin = -0.5*line_info.xwidth_i;
        xmax =  0.5*line_info.xwidth_i;
        ymin = -0.5*line_info.yheight_i;
        ymax =  0.5*line_info.yheight_i;
        zmin = -0.5*line_info.zdepth_i;
        zmax =  0.5*line_info.zdepth_i;
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
    } if (line_info.shape == 2) { /* sphere */
      if (line_info.radius_i) {
        circle("xy",0,0,0,line_info.radius_i);
        circle("xz",0,0,0,line_info.radius_i);
        circle("yz",0,0,0,line_info.radius_i);
      }
      circle("xy",0,0,0,radius);
      circle("xz",0,0,0,radius);
      circle("yz",0,0,0,radius);
    } else if (line_info.shape == 3) {	/* OFF file */
      off_display(offdata);
    }
  }
  #undef reflections
  #undef geometry
  #undef format
  #undef radius
  #undef yheight
  #undef xwidth
  #undef zdepth
  #undef thickness
  #undef pack
  #undef Vc
  #undef sigma_abs
  #undef sigma_inc
  #undef delta_d_d
  #undef p_inc
  #undef p_transmit
  #undef DW
  #undef nb_atoms
  #undef d_phi
  #undef p_interact
  #undef concentric
  #undef density
  #undef weight
  #undef barns
  #undef Strain
  #undef focus_flip
  #undef line_info
  #undef columns
  #undef offdata
  return(_comp);
} /* class_PowderN_display */

_class_Collimator_radial *class_Collimator_radial_display(_class_Collimator_radial *_comp
) {
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define length (_comp->_parameters.length)
  #define divergence (_comp->_parameters.divergence)
  #define transmission (_comp->_parameters.transmission)
  #define theta_min (_comp->_parameters.theta_min)
  #define theta_max (_comp->_parameters.theta_max)
  #define nchan (_comp->_parameters.nchan)
  #define radius (_comp->_parameters.radius)
  #define nslit (_comp->_parameters.nslit)
  #define roc (_comp->_parameters.roc)
  #define verbose (_comp->_parameters.verbose)
  #define approx (_comp->_parameters.approx)
  #define width_of_slit (_comp->_parameters.width_of_slit)
  #define width_of_Soller (_comp->_parameters.width_of_Soller)
  #define slit_theta (_comp->_parameters.slit_theta)
  double Soller_theta;
  double height=yheight/2;
  double theta1, theta2;
  double x_in_l,  z_in_l,  x_in_r,  z_in_r;
  double x_out_l, z_out_l, x_out_r, z_out_r;
  int i;

  /* display collimator radial geometry:
     in order to avoid too many lines, we shown main housing and channels
     but no slit */

  if (!nchan || nchan > 20)     nchan=20;
  if (nchan > 64)               nchan=64;
  Soller_theta=fabs(theta_max-theta_min)/nchan; /* angular width of Soller */

  /* draw all channels, which also show housing */
  
  for (i = 0; i < nchan; i++) {

    theta1 = i*Soller_theta+theta_min;
    theta2 = theta1+Soller_theta;

    z_in_l = radius*cos(theta1);
    x_in_l = radius*sin(theta1);
    z_in_r = radius*cos(theta2);
    x_in_r = radius*sin(theta2);

    z_out_l = (radius+length)*cos(theta1);
    x_out_l = (radius+length)*sin(theta1);
    z_out_r = (radius+length)*cos(theta2);
    x_out_r = (radius+length)*sin(theta2);
    /* left side */
    multiline(6,
      x_in_l, -height, z_in_l,
      x_in_l,  height, z_in_l,
      x_out_l, height, z_out_l,
      x_out_l,-height, z_out_l,
      x_in_l, -height, z_in_l,
      x_in_r, -height, z_in_r);
   /* left -> right lines */
   line(x_in_l,   height, z_in_l,  x_in_r,  height, z_in_r);
   line(x_out_l,  height, z_out_l, x_out_r, height, z_out_r);
   line(x_out_l, -height, z_out_l, x_out_r,-height, z_out_r);
  }
  /* remaining bits */
  theta1 = nchan*Soller_theta+theta_min;
  z_in_l = radius*cos(theta1);
  x_in_l = radius*sin(theta1);
  z_out_l = (radius+length)*cos(theta1);
  x_out_l = (radius+length)*sin(theta1);
  multiline(5,
      x_in_l, -height, z_in_l,
      x_in_l,  height, z_in_l,
      x_out_l, height, z_out_l,
      x_out_l,-height, z_out_l,
      x_in_l, -height, z_in_l);
  #undef xwidth
  #undef yheight
  #undef length
  #undef divergence
  #undef transmission
  #undef theta_min
  #undef theta_max
  #undef nchan
  #undef radius
  #undef nslit
  #undef roc
  #undef verbose
  #undef approx
  #undef width_of_slit
  #undef width_of_Soller
  #undef slit_theta
  return(_comp);
} /* class_Collimator_radial_display */


  #undef magnify
  #undef line
  #undef dashed_line
  #undef multiline
  #undef rectangle
  #undef box
  #undef circle
  #undef cylinder
  #undef sphere

int display(void) { /* called by mccode_main for SAFARI_PITSI:DISPLAY */
  printf("MCDISPLAY: start\n");

  /* call iteratively all components DISPLAY */
  class_Progress_bar_display(_Progress);

  class_Arm_display(_Reactorbeam);

  class_Arm_display(_Prim_axes);

  class_Source_gen_display(_Source);

  class_PSD_monitor_display(_PSD_Source);

  class_L_monitor_display(_LAM_Source);

  class_Slit_display(_Window_before_filter);

  class_Filter_gen_display(_Sapphire_filter);

  class_Slit_display(_Window_after_filter);

  class_L_monitor_display(_LAM_After_sapphire);

  class_PSD_monitor_display(_PSD_After_sapphire);

  class_Slit_display(_HighResOutlet);

  class_Slit_display(_HighIntensityOutlet);

  class_PSD_monitor_display(_PSD_After_Outlet);

  class_Arm_display(_Mono_axis);

  class_Monochromator_curved_display(_Blade_1);

  class_Monochromator_curved_display(_Blade_2);

  class_Monochromator_curved_display(_Blade_3);

  class_Monochromator_curved_display(_Blade_4);

  class_Monochromator_curved_display(_Blade_5);

  class_Monochromator_curved_display(_Blade_6);

  class_Monochromator_curved_display(_Blade_7);

  class_Monochromator_curved_display(_Blade_8);

  class_Monochromator_curved_display(_Blade_9);

  class_Monochromator_curved_display(_Blade_10);

  class_Monochromator_curved_display(_Blade_11);

  class_Monochromator_curved_display(_Blade_12);

  class_Monochromator_curved_display(_Blade_13);

  class_PSD_monitor_display(_PSD_At_sec_shutter);

  class_L_monitor_display(_LAM_At_sec_shutter);

  class_Slit_display(_Inside_chamber_collimator);

  class_PSD_monitor_display(_PSD_After_inside_chamber_collimator);

  class_Slit_display(_Outside_chamber_collimator);

  class_PSD_monitor_display(_PSD_Outside_chamber_collimator_1);

  class_PSD_monitor_display(_PSD_Outside_chamber_collimator);

  class_Slit_display(_Incident_slit_h);

  class_Slit_display(_Incident_slit_w);

  class_PSD_monitor_display(_PSD_After_Incident_slit_w);

  class_Arm_display(_Center_of_rotation);

  class_Arm_display(_Sample_rotation);

  class_Arm_display(_Sample_location);

  class_PSD_monitor_display(_PSD_Center_of_rotation);

  class_Divergence_monitor_display(_DIV_Center_of_rotation);

  class_PowderN_display(_Sample);

  class_Arm_display(_Det_axis);

  class_Arm_display(_Det_axis_2);

  class_Arm_display(_Det_axis_3);

  class_Arm_display(_Det_axis_4);

  class_Collimator_radial_display(_RadColl);

  class_PSD_monitor_display(_PSD_Detector);

  class_PSD_monitor_display(_PSD_Detector_2);

  class_PSD_monitor_display(_PSD_Detector_3);

  class_PSD_monitor_display(_PSD_Detector_4);

  printf("MCDISPLAY: end\n");

  return(0);
} /* display */

void* _getvar_parameters(char* compname)
/* enables settings parameters based use of the GETPAR macro */
{
  if (!strcmp(compname, "Progress")) return (void *) &(_Progress->_parameters);
  if (!strcmp(compname, "Reactorbeam")) return (void *) &(_Reactorbeam->_parameters);
  if (!strcmp(compname, "Prim_axes")) return (void *) &(_Prim_axes->_parameters);
  if (!strcmp(compname, "Source")) return (void *) &(_Source->_parameters);
  if (!strcmp(compname, "PSD_Source")) return (void *) &(_PSD_Source->_parameters);
  if (!strcmp(compname, "LAM_Source")) return (void *) &(_LAM_Source->_parameters);
  if (!strcmp(compname, "Window_before_filter")) return (void *) &(_Window_before_filter->_parameters);
  if (!strcmp(compname, "Sapphire_filter")) return (void *) &(_Sapphire_filter->_parameters);
  if (!strcmp(compname, "Window_after_filter")) return (void *) &(_Window_after_filter->_parameters);
  if (!strcmp(compname, "LAM_After_sapphire")) return (void *) &(_LAM_After_sapphire->_parameters);
  if (!strcmp(compname, "PSD_After_sapphire")) return (void *) &(_PSD_After_sapphire->_parameters);
  if (!strcmp(compname, "HighResOutlet")) return (void *) &(_HighResOutlet->_parameters);
  if (!strcmp(compname, "HighIntensityOutlet")) return (void *) &(_HighIntensityOutlet->_parameters);
  if (!strcmp(compname, "PSD_After_Outlet")) return (void *) &(_PSD_After_Outlet->_parameters);
  if (!strcmp(compname, "Mono_axis")) return (void *) &(_Mono_axis->_parameters);
  if (!strcmp(compname, "Blade_1")) return (void *) &(_Blade_1->_parameters);
  if (!strcmp(compname, "Blade_2")) return (void *) &(_Blade_2->_parameters);
  if (!strcmp(compname, "Blade_3")) return (void *) &(_Blade_3->_parameters);
  if (!strcmp(compname, "Blade_4")) return (void *) &(_Blade_4->_parameters);
  if (!strcmp(compname, "Blade_5")) return (void *) &(_Blade_5->_parameters);
  if (!strcmp(compname, "Blade_6")) return (void *) &(_Blade_6->_parameters);
  if (!strcmp(compname, "Blade_7")) return (void *) &(_Blade_7->_parameters);
  if (!strcmp(compname, "Blade_8")) return (void *) &(_Blade_8->_parameters);
  if (!strcmp(compname, "Blade_9")) return (void *) &(_Blade_9->_parameters);
  if (!strcmp(compname, "Blade_10")) return (void *) &(_Blade_10->_parameters);
  if (!strcmp(compname, "Blade_11")) return (void *) &(_Blade_11->_parameters);
  if (!strcmp(compname, "Blade_12")) return (void *) &(_Blade_12->_parameters);
  if (!strcmp(compname, "Blade_13")) return (void *) &(_Blade_13->_parameters);
  if (!strcmp(compname, "PSD_At_sec_shutter")) return (void *) &(_PSD_At_sec_shutter->_parameters);
  if (!strcmp(compname, "LAM_At_sec_shutter")) return (void *) &(_LAM_At_sec_shutter->_parameters);
  if (!strcmp(compname, "Inside_chamber_collimator")) return (void *) &(_Inside_chamber_collimator->_parameters);
  if (!strcmp(compname, "PSD_After_inside_chamber_collimator")) return (void *) &(_PSD_After_inside_chamber_collimator->_parameters);
  if (!strcmp(compname, "Outside_chamber_collimator")) return (void *) &(_Outside_chamber_collimator->_parameters);
  if (!strcmp(compname, "PSD_Outside_chamber_collimator_1")) return (void *) &(_PSD_Outside_chamber_collimator_1->_parameters);
  if (!strcmp(compname, "PSD_Outside_chamber_collimator")) return (void *) &(_PSD_Outside_chamber_collimator->_parameters);
  if (!strcmp(compname, "Incident_slit_h")) return (void *) &(_Incident_slit_h->_parameters);
  if (!strcmp(compname, "Incident_slit_w")) return (void *) &(_Incident_slit_w->_parameters);
  if (!strcmp(compname, "PSD_After_Incident_slit_w")) return (void *) &(_PSD_After_Incident_slit_w->_parameters);
  if (!strcmp(compname, "Center_of_rotation")) return (void *) &(_Center_of_rotation->_parameters);
  if (!strcmp(compname, "Sample_rotation")) return (void *) &(_Sample_rotation->_parameters);
  if (!strcmp(compname, "Sample_location")) return (void *) &(_Sample_location->_parameters);
  if (!strcmp(compname, "PSD_Center_of_rotation")) return (void *) &(_PSD_Center_of_rotation->_parameters);
  if (!strcmp(compname, "DIV_Center_of_rotation")) return (void *) &(_DIV_Center_of_rotation->_parameters);
  if (!strcmp(compname, "Sample")) return (void *) &(_Sample->_parameters);
  if (!strcmp(compname, "Det_axis")) return (void *) &(_Det_axis->_parameters);
  if (!strcmp(compname, "Det_axis_2")) return (void *) &(_Det_axis_2->_parameters);
  if (!strcmp(compname, "Det_axis_3")) return (void *) &(_Det_axis_3->_parameters);
  if (!strcmp(compname, "Det_axis_4")) return (void *) &(_Det_axis_4->_parameters);
  if (!strcmp(compname, "RadColl")) return (void *) &(_RadColl->_parameters);
  if (!strcmp(compname, "PSD_Detector")) return (void *) &(_PSD_Detector->_parameters);
  if (!strcmp(compname, "PSD_Detector_2")) return (void *) &(_PSD_Detector_2->_parameters);
  if (!strcmp(compname, "PSD_Detector_3")) return (void *) &(_PSD_Detector_3->_parameters);
  if (!strcmp(compname, "PSD_Detector_4")) return (void *) &(_PSD_Detector_4->_parameters);
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

/* end of generated C code SAFARI_PITSI.c */
