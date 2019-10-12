/* Automatically generated file. Do not edit. 
 * Format:     ANSI C source code
 * Creator:    McStas <http://www.mcstas.org>
 * Instrument: ESS_2015_test.instr (ESS_2015_test)
 * Date:       Sat Oct 12 09:06:54 2019
 * File:       ESS_2015_test.c
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
* Start of instrument 'ESS_2015_test' generated code
***************************************************************************** */

#ifdef MC_TRACE_ENABLED
int traceenabled = 1;
#else
int traceenabled = 0;
#endif
#define MCSTAS "/usr/share/mcstas/3.0-dev/"
int   defaultmain         = 1;
char  instrument_name[]   = "ESS_2015_test";
char  instrument_source[] = "ESS_2015_test.instr";
char *instrument_exe      = NULL; /* will be set to argv[0] in main */
char  instrument_code[]   = "Instrument ESS_2015_test source code ESS_2015_test.instr is not embedded in this executable.\n  Use --source option when running McStas.\n";

int main(int argc, char *argv[]){return mccode_main(argc, argv);}

/* *****************************************************************************
* instrument 'ESS_2015_test' and components DECLARE
***************************************************************************** */

/* Instrument parameters: structure and a table for the initialisation
   (Used in e.g. inputparse and I/O function (e.g. detector_out) */

struct _struct_instrument_parameters {
  MCNUM _frac;
  MCNUM _lmin;
  MCNUM _lmax;
  MCNUM _ANGLE;
  MCNUM _mon_shift;
};
typedef struct _struct_instrument_parameters _class_instrument_parameters;

struct _instrument_struct {
  char   _name[256]; /* the name of this instrument e.g. 'ESS_2015_test' */
/* Counters per component instance */
  double counter_AbsorbProp[12]; /* absorbed events in PROP routines */
  double counter_N[12], counter_P[12], counter_P2[12]; /* event counters after each component instance */
  _class_particle _trajectory[12]; /* current trajectory for STORE/RESTORE */
/* Components position table (absolute and relative coords) */
  Coords _position_relative[12]; /* positions of all components */
  Coords _position_absolute[12];
  _class_instrument_parameters _parameters; /* instrument parameters */
} _instrument_var;
struct _instrument_struct *instrument = & _instrument_var;
#pragma acc declare create ( _instrument_var )
#pragma acc declare create ( instrument )

int numipar = 5;
struct mcinputtable_struct mcinputtable[] = {
  "frac", &(_instrument_var._parameters._frac), instr_type_double, "0.5", 
  "lmin", &(_instrument_var._parameters._lmin), instr_type_double, "0.3", 
  "lmax", &(_instrument_var._parameters._lmax), instr_type_double, "8.3", 
  "ANGLE", &(_instrument_var._parameters._ANGLE), instr_type_double, "-55", 
  "mon_shift", &(_instrument_var._parameters._mon_shift), instr_type_double, "0.0", 
  NULL, NULL, instr_type_double, ""
};


/* ************************************************************************** */
/*             SHARE user declarations for all components                     */
/* ************************************************************************** */

/* Shared user declarations for all components types 'ESS_moderator'. */
/*******************************************************************************
*
* McStas, neutron ray-tracing package
*         Copyright 1997-2013, All rights reserved
*         DTU Physics, Lyngby, Denmark
*         Institut Laue Langevin, Grenoble, France
*
* Library: share/ess_source-lib.h
*
* %Identification
* Written by: PW
* Date: Nov 7, 2013
* Origin: DTU Physics
* Release: McStas 2.1
* Version: 0.1
*
* This file is to be imported by the ESS_moderator_long component
* It defines a set of brilliance definitions (used via function pointer) for
* easier use of the component.
*
* Usage: within SHARE
* %include "ess_source-lib"
*
*******************************************************************************/

#ifndef ESS_SOURCE_LIB_H
#define ESS_SOURCE_LIB_H 0.1
#define ESS_SOURCE_DURATION 2.857e-3
#define ESS_SOURCE_FREQUENCY 14
#define ESS_SOURCE_POWER 5


/* Struct for extra source parameters - for future geometrical adjustments */
struct ess_struct {
  double X;
  double Y;
  double Z;
  double height_t;
  double height_c;
  double Width_c;
  double Width_t;
  double tmultiplier;
  double Radius_c;
  double beamportangle;
  int Uniform;
  double extractionangle;
  int Wasleft;
};

typedef struct ess_struct ess_moderator_struct;

typedef void (*functype)(double* t , double* p, double lambda,  double tfocus_w, double tfocus_t, double tfocus_dt, ess_moderator_struct extras);

/* Maths for the curves below */
double Mezei_M(double l, double temp);
double Mezei_F(double t, double tau, int n);
double ESS_2013_Schoenfeldt_cold_spectrum(double I_SD, double alpha_SD, double lambda_SD, double alpha_l, double lambda_l, double Exponent, double I_1, double alpha_1, double I_2, double alpha_2, double lambda);
double ESS_2013_Schoenfeldt_thermal_spectrum(double I_th, double T, double I_SD, double alpha, double lambda_cf, double lambda);
double ESS_2014_Schoenfeldt_cold_spectrum(double lambda,double height);
double ESS_2014_Schoenfeldt_thermal_spectrum(double lambda, double height);
double ESS_2015_Schoenfeldt_cold_spectrum(double lambda,double theta);
double ESS_2015_Schoenfeldt_thermal_spectrum(double lambda, double theta);
double ESS_2014_Schoenfeldt_cold_Geom_120_over_60=1;
double ESS_2014_Schoenfeldt_thermal_Geom_120_over_60=1;

/* List of brilliance definitions */
void ESS_Mezei_cold(double *t, double *p, double lambda, double tfocus_w, double tfocus_t, double tfocus_dt, ess_moderator_struct extras); 
void ESS_Mezei_thermal(double *t, double *p, double lambda, double tfocus_w, double tfocus_t, double tfocus_dt, ess_moderator_struct extras); 
void ESS_Mezei_cold_2012(double *t, double *p, double lambda, double tfocus_w, double tfocus_t, double tfocus_dt, ess_moderator_struct extras); 
void ESS_2012_Lieutenant_cold(double *t, double *p, double lambda, double tfocus_w, double tfocus_t, double tfocus_dt, ess_moderator_struct extras); 
void ESS_2013_Schoenfeldt_cold(double *t, double *p, double lambda, double tfocus_w, double tfocus_t, double tfocus_dt, ess_moderator_struct extras); 
void ESS_2013_Schoenfeldt_thermal(double *t, double *p, double lambda, double tfocus_w, double tfocus_t, double tfocus_dt, ess_moderator_struct extras);
void ESS_2014_Schoenfeldt_cold(double *t, double *p, double lambda, double tfocus_w, double tfocus_t, double tfocus_dt, ess_moderator_struct extras); 
void ESS_2014_Schoenfeldt_thermal(double *t, double *p, double lambda, double tfocus_w, double tfocus_t, double tfocus_dt, ess_moderator_struct extras);
void ESS_2015_Schoenfeldt_cold(double *t, double *p, double lambda, double tfocus_w, double tfocus_t, double tfocus_dt, ess_moderator_struct extras); 
void ESS_2015_Schoenfeldt_thermal(double *t, double *p, double lambda, double tfocus_w, double tfocus_t, double tfocus_dt, ess_moderator_struct extras);

/* List of pulse-shape definitions */
double ESS_2014_Schoenfeldt_cold_timedist(double t, double lambda, double height, double pulselength);
double ESS_2014_Schoenfeldt_thermal_timedist(double t, double lambda, double height, double pulselength);
double ESS_2015_Schoenfeldt_cold_timedist(double t, double lambda, double height, double pulselength);
double ESS_2015_Schoenfeldt_thermal_timedist(double t, double lambda, double height, double pulselength);

/* List of moderator-geometry-weighting definitions */
double ESS_2014_Schoenfeldt_cold_y0(double y0,double height);
double ESS_2014_Schoenfeldt_cold_x0(double x0,double height, double width);
double ESS_2014_Schoenfeldt_thermal_y0(double y0,double height);
double ESS_2014_Schoenfeldt_thermal_x0(double x0,double height, double width);
double ESS_2014_Schoenfeldt_cold_Y(double x0,double height);
double ESS_2014_Schoenfeldt_thermal_Y(double y0,double height);
double ESS_2014_Schoenfeldt_cold_Theta120(double x0,double height);
double ESS_2014_Schoenfeldt_thermal_Theta120(double beamportangle,int isleft);
double ESS_2015_Schoenfeldt_cold_y0(double y0);
double ESS_2015_Schoenfeldt_cold_x0(double x0,double theta);
double ESS_2015_Schoenfeldt_thermal_y0(double y0);
double ESS_2015_Schoenfeldt_thermal_x0(double x0,double theta);
double ESS_2015_Schoenfeldt_cold_Y(double x0,double height);
double ESS_2015_Schoenfeldt_thermal_Y(double y0,double height);
double ESS_2015_Schoenfeldt_cold_Theta120(double x0,double height);
double ESS_2015_Schoenfeldt_thermal_Theta120(double beamportangle,int isleft);

/* List of geometry definitions - mainly for mcdisplay... */
void ESS_mcdisplay_flat(double geometry);
void ESS_mcdisplay_TDRlike(double geometry);

double TSC2015_z0_BF3cm(const double x0);
/* end of ess_source-lib.h */
#endif

/*******************************************************************************
*
* McStas, neutron ray-tracing package
*         Copyright 1997-2013, All rights reserved
*         DTU Physics, Lyngby, Denmark
*         Institut Laue Langevin, Grenoble, France
*
* Library: share/ess_source-lib.c
*
* %Identification
* Written by: PW
* Date: Nov 7, 2013
* Origin: DTU Physics
* Release: McStas 2.1
* Version: 0.1
*
* This file is to be imported by the ESS_moderator_long component
* It defines a set of brilliance definitions (used via function pointer) for
* easier use of the component.
*
* Usage: within SHARE
* %include "ess_source-lib"
*
*******************************************************************************/

#ifndef ESS_SOURCE_LIB_H
#error McStas : please import this library with %include "ess_source-lib"
#endif

/* Base maths functions used below */
double Mezei_M(double l, double temp)
    {
      double a=949.0/temp;
      return 2*a*a*exp(-a/(l*l))/(l*l*l*l*l);
    }

/* This one is a bit strange - never used in the actual comp it seems? - Should rewrite.... */
double Mezei_F(double t, double tau, int n)
    {
      return (exp(-t/tau)-exp(-n*t/tau))*n/(n-1)/tau;
    }

double ESS_2013_Schoenfeldt_cold_spectrum(double I_SD, double alpha_SD, double lambda_SD, double alpha_l, double lambda_l, double Exponent, double I_1, double alpha_1, double I_2, double alpha_2, double lambda) 
{
  return (I_1 * exp(-alpha_1 * lambda) + I_2 * exp(-alpha_2 * lambda)) * 1 / pow(1 + exp(alpha_l * (lambda - lambda_l)),-Exponent) + I_SD * (1/lambda) * 1/( 1 + exp(alpha_SD * (lambda - lambda_SD)));
}

double ESS_2013_Schoenfeldt_thermal_spectrum(double I_th, double T, double I_SD, double alpha, double lambda_cf, double lambda) 
{
  double k_Th=949;
  return I_th * exp(-k_Th/(T*lambda*lambda))*(2*k_Th*k_Th)/(T*T*pow(lambda,5)) + I_SD * (1/lambda) * 1/(1+exp(alpha*(lambda - lambda_cf)));
}

double ESS_2014_Schoenfeldt_cold_spectrum(double lambda,double height){
if(lambda<=0)return 0;
	return pow((1+exp((-7.84092e+001*exp(-1.00000e+000*height)-8.46887e+000+3)
		*(lambda-(2.48787-0.00729329*height)))),(1.17068e+000*exp(-2.52666e-001*height)-9.29703e-001))
		*((5.17542e+014*exp(-3.91646e-001*height)+7.19417e+013)
		*(exp(-0.77*(lambda))+(-5.03315e-002*exp(-5.79174e-001*height)+5.72765e-002)*exp(-0.323784*(lambda))))
		+(4.12442e+012*exp(-1.25186e-001*height)+8.53849e+011)/((1+exp(2.5*(lambda-2.2)))*lambda);
}

double ESS_2014_Schoenfeldt_thermal_spectrum(double lambda, double height){
	if(lambda<=0)return 0;
	double aOlsqr=949./(325*lambda*lambda);
	return 2*(3.46910e+013*exp(-1.65602e-001*height)+8.08542e+012)
	  *aOlsqr*aOlsqr/lambda*pow(lambda,(3.11752e-001*(1-exp(-3.45363e-001*height))+9.17072e-002))
	  *exp(-aOlsqr)

	  +

	  (4.64873e+012*exp(-1.80747e-001*height)+1.74845e+012)/((1+exp(2.5*(lambda-0.88)))*lambda);
}

double ESS_2015_Schoenfeldt_cold_spectrum(double lambda,double theta){
  if(lambda<=0)return 0;
  double par0=8.44e13/25.;
  double par1=2.5;
  double par2=2.2;
  
  double par3=-13.-.5*(theta-5);
  double par4=2.53;
  double par5=-0.0478073-0.160*exp(-0.45186*(theta-5.)/10.);
  
  double par6;
  if(theta==5)par6=5.73745e+015/25.;
  else if(theta==15)par6=5.88284e+015/25.;
  else if(theta==25)par6=6.09573e+015/25.;
  else if(theta==35)par6=6.29116e+015/25.;
  else if(theta==45)par6=6.03436e+015/25.;
  else if(theta==55)par6=6.02045e+015/25.;
  double par7=0.788956+0.00854184*(theta-5.)/10.;
  double par8=0.0461868-0.0016464*(theta-5.)/10.;
  double par9=0.325;
  
  double SD_part=par0/((1+exp(par1*(lambda-par2)))*lambda);
  double para_part=pow((1+exp(par3*(lambda-par4))),par5)*(par6*(exp(-par7*(lambda))+par8*exp(-par9*(lambda))));
  return para_part+SD_part;
  
}

double ESS_2015_Schoenfeldt_thermal_spectrum(double lambda, double theta){
    if(lambda<=0)return 0;
    double i=(theta-5.)/10.;
    double par0=4.2906e+013-9.2758e+011*i+8.02603e+011*i*i-1.29523e+011*i*i*i;
    double par2=6.24806e+012-8.84602e+010*i;
    double par3=-0.31107+0.0221138*i;
    double aOlsqr=949./(325*lambda*lambda);
    return par0*2.*aOlsqr*aOlsqr/lambda*pow(lambda,-par3)*exp(-aOlsqr)+par2/((1+exp(2.5*(lambda-0.88)))*lambda);
	  
}


/* This is the cold Mezei moderator from 2012 (updated I0 and I2) */
void ESS_Mezei_cold_2012(double *t, double *p, double lambda, double tfocus_width, double tfocus_time, double dt, ess_moderator_struct extras) 
{
  // Spectrum related constants - ESS 2001 Cold moderator
  double T=50, tau=287e-6, tau1=0, tau2=20e-6, chi2=0.9, I0=8.21e11, I2=3.29e11, branch1=1, branch2=0.5, n2=5, n=20;
  
  // Branching
  double branch_tail=tau/ESS_SOURCE_DURATION;
  
  // Other variables
  double tail_flag, tau_l;
  
  // Taken directly from the ESS_moderator.comp:
  tail_flag = (rand01()<branch_tail);   /* Choose tail/bulk */
  if (tail_flag)
    {
      if (rand01() < branch2)
	{
	  if (tau1>0)
	    if (rand01() < branch1)     /* Quick and dirty non-general solution */
	      {  /* FIRST CASE a */
		tau_l = tau;
		*p = 1/(branch1*branch2*branch_tail); /* Correct for switching prob. */
	      }
	    else
	      {  /* FIRST CASE b */
		tau_l = tau1;
		*p = 1/((1-branch1)*branch2*branch_tail); /* Correct for switching prob. */
	      }
	  else
	    {
	      tau_l = tau;
	      *p = 1/(branch2*branch_tail); /* Correct for switching prob. */
	    }
	  *t = -tau_l*log(1e-12+rand01());       /* Sample from long-time tail a */
	  /* Correct for true pulse shape */
	  //	  p *= w_focus;                         /* Correct for target focusing */
	  *p *= tau_l/ESS_SOURCE_DURATION;                         /* Correct for tail part */
	  //p *= I0*w_mult*w_geom*Mezei_M(lambda,T);           /* Calculate true intensity */
	  *p *= I0*Mezei_M(lambda,T);
	}
      else
	{
	  /* SECOND CASE */
	  tau_l = tau2*lambda;
	  *t = -tau_l*log(1e-12+rand01());       /* Sample from long-time tail */
	  *p = n2/(n2-1)*((1-exp(-ESS_SOURCE_DURATION/tau_l))-(1-exp(-n2*ESS_SOURCE_DURATION/tau_l))*exp(-(n2-1)*(*t)/tau_l)/n);
	  /* Correct for true pulse shape */
	  *p /= (1-branch2)*branch_tail;          /* Correct for switching prob. */
	  *p *= tau_l/ESS_SOURCE_DURATION;                         /* Correct for tail part */
	  // p *= w_focus;                         /* Correct for target focusing */
	  //p *= I2*w_mult*w_geom/(1+exp(chi2*lambda-2.2))/lambda;                                         /* Calculate true intensity */
	  *p *= I2/(1+exp(chi2*lambda-2.2))/lambda;                                         /* Calculate true intensity */ 
	}
      *t += ESS_SOURCE_DURATION;                                 /* Add pulse length */
    }
  else /* Tail-flag */
    {
      if (tfocus_width>0) {
	*t = tfocus_time-dt;                    /* Set time to hit time window center */
	*t += randpm1()*tfocus_width/2.0;       /* Add random time within window width */
      } else {
	*t = ESS_SOURCE_DURATION*rand01();                        /* Sample from bulk pulse */
      }
      // FLAG to KILL these on return!

      /* if (t<0) ABSORB;                       /\* Kill neutron if outside pulse duration *\/ */
      /* if (t>ESS_SOURCE_DURATION) ABSORB; */
      if (rand01() < branch2)
	{
	  if (rand01() < branch1)     /* Quick and dirty non-general solution */
	    {  /* FIRST CASE a */
	      tau_l = tau;
	      *p = 1/(branch1*branch2*(1-branch_tail)); /* Correct for switching prob. */
	    }
	  else
	    {  /* FIRST CASE b */
	      tau_l = tau1;
	      *p = 1/((1-branch1)*branch2*(1-branch_tail)); /* Correct for switching prob. */
	    }
	  *p *= 1-n/(n-1)*(exp(-*t/tau_l)-exp(-n*(*t)/tau_l)/n); /* Correct for true pulse shape */
	  //	  p *= w_focus;                         /* Correct for target focusing */
	  if (tfocus_width>0) {
	    *p *= tfocus_width/ESS_SOURCE_DURATION;    	  	  /* Correct for time focusing */
	  }
	  //p *= I0*w_mult*w_geom*M(lambda,T);       /* Calculate true intensity */
	  *p *= I0*Mezei_M(lambda,T);       /* Calculate true intensity */
	}
      else
	{
	  /* SECOND CASE */
	  tau_l = tau2*lambda;
	  *p = 1-n2/(n2-1)*(exp(-*t/tau_l)-exp(-n2*(*t)/tau_l)/n2); /* Correct for true pulse shape */
	  *p /= (1-branch2)*(1-branch_tail);   /* Correct for switching prob. */
	  //p *= w_focus;                         /* Correct for target focusing */
	  if (tfocus_width) {
	    *p *= tfocus_width/ESS_SOURCE_DURATION;    		  /* Correct for time focusing */
	  }
	  //p *= I2*w_mult*w_geom/(1+exp(chi2*lambda-2.2))/lambda;    /* Calculate true intensity */
	  *p *= I2/(1+exp(chi2*lambda-2.2))/lambda;    /* Calculate true intensity */
	}
    }
  
} /* end of ESS_Mezei_cold_2012 */

/* This is the cold Mezei moderator from 2001 (Original I0 and I2) */
void ESS_Mezei_cold(double *t, double *p, double lambda, double tfocus_width, double tfocus_time, double dt, ess_moderator_struct extras) 
{
  // Spectrum related constants - ESS 2001 Cold moderator
  double T=50, tau=287e-6, tau1=0, tau2=20e-6, chi2=0.9, I0=6.9e11, I2=27.6e10, branch1=1, branch2=0.5, n2=5, n=20;
  
  // Branching
  double branch_tail=tau/ESS_SOURCE_DURATION;
  
  // Other variables
  double tail_flag, tau_l;
  
  // Taken directly from the ESS_moderator.comp:
  tail_flag = (rand01()<branch_tail);   /* Choose tail/bulk */
  if (tail_flag)
    {
      if (rand01() < branch2)
	{
	  if (tau1>0)
	    if (rand01() < branch1)     /* Quick and dirty non-general solution */
	      {  /* FIRST CASE a */
		tau_l = tau;
		*p = 1/(branch1*branch2*branch_tail); /* Correct for switching prob. */
	      }
	    else
	      {  /* FIRST CASE b */
		tau_l = tau1;
		*p = 1/((1-branch1)*branch2*branch_tail); /* Correct for switching prob. */
	      }
	  else
	    {
	      /* This is the only active part - tau1 is set to 0! */
	      tau_l = tau;
	      *p = 1/(branch2*branch_tail); /* Correct for switching prob. */
	    }
	  *t = -tau_l*log(1e-12+rand01());       /* Sample from long-time tail a */
	  /* Correct for true pulse shape */
	  //	  p *= w_focus;                         /* Correct for target focusing */
	  *p *= tau_l/ESS_SOURCE_DURATION;                         /* Correct for tail part */
	  //p *= I0*w_mult*w_geom*Mezei_M(lambda,T);           /* Calculate true intensity */
	  *p *= I0*Mezei_M(lambda,T);
	}
      else
	{
	  /* Shine-through from spallation, i.e. not T-dependent */
	  /* SECOND CASE */
	  tau_l = tau2*lambda;
	  *t = -tau_l*log(1e-12+rand01());       /* Sample from long-time tail */
	  *p = n2/(n2-1)*((1-exp(-ESS_SOURCE_DURATION/tau_l))-(1-exp(-n2*ESS_SOURCE_DURATION/tau_l))*exp(-(n2-1)*(*t)/tau_l)/n);
	  /* Correct for true pulse shape */
	  *p /= (1-branch2)*branch_tail;          /* Correct for switching prob. */
	  *p *= tau_l/ESS_SOURCE_DURATION;                         /* Correct for tail part */
	  // p *= w_focus;                         /* Correct for target focusing */
	  //p *= I2*w_mult*w_geom/(1+exp(chi2*lambda-2.2))/lambda;                                         /* Calculate true intensity */
	  *p *= I2/(1+exp(chi2*lambda-2.2))/lambda;                                         /* Calculate true intensity */ 
	}
      *t += ESS_SOURCE_DURATION;                                 /* Add pulse length */
    }
  else /* Tail-flag */
    {
      if (tfocus_width>0) {
	*t = tfocus_time-dt;                    /* Set time to hit time window center */
	*t += randpm1()*tfocus_width/2.0;       /* Add random time within window width */
      } else {
	*t = ESS_SOURCE_DURATION*rand01();                        /* Sample from bulk pulse */
      }
      // FLAG to KILL these on return!

      /* if (t<0) ABSORB;                       /\* Kill neutron if outside pulse duration *\/ */
      /* if (t>ESS_SOURCE_DURATION) ABSORB; */
      if (rand01() < branch2)
	{
	  if (rand01() < branch1)     /* Quick and dirty non-general solution */
	    {  /* FIRST CASE a */
	      tau_l = tau;
	      *p = 1/(branch1*branch2*(1-branch_tail)); /* Correct for switching prob. */
	    }
	  else
	    {  /* FIRST CASE b */
	      tau_l = tau1;
	      *p = 1/((1-branch1)*branch2*(1-branch_tail)); /* Correct for switching prob. */
	    }
	  *p *= 1-n/(n-1)*(exp(-*t/tau_l)-exp(-n*(*t)/tau_l)/n); /* Correct for true pulse shape */
	  //	  p *= w_focus;                         /* Correct for target focusing */
	  if (tfocus_width>0) {
	    *p *= tfocus_width/ESS_SOURCE_DURATION;    	  	  /* Correct for time focusing */
	  }
	  //p *= I0*w_mult*w_geom*M(lambda,T);       /* Calculate true intensity */
	  *p *= I0*Mezei_M(lambda,T);       /* Calculate true intensity */
	}
      else
	{
	  /* SECOND CASE */
	  /* Shine-through from spallation, i.e. not T-dependent */
	  tau_l = tau2*lambda;
	  *p = 1-n2/(n2-1)*(exp(-*t/tau_l)-exp(-n2*(*t)/tau_l)/n2); /* Correct for true pulse shape */
	  *p /= (1-branch2)*(1-branch_tail);   /* Correct for switching prob. */
	  //p *= w_focus;                         /* Correct for target focusing */
	  if (tfocus_width) {
	    *p *= tfocus_width/ESS_SOURCE_DURATION;    		  /* Correct for time focusing */
	  }
	  //p *= I2*w_mult*w_geom/(1+exp(chi2*lambda-2.2))/lambda;    /* Calculate true intensity */
	  *p *= I2/(1+exp(chi2*lambda-2.2))/lambda;    /* Calculate true intensity */
	}
    }
  
} /* end of ESS_Mezei_cold */

/* This is the thermal Mezei moderator from 2001 - also used in 2012 - TDR */
void ESS_Mezei_thermal(double *t, double *p, double lambda, double tfocus_width, double tfocus_time, double dt, ess_moderator_struct extras)
{
  // Spectrum related constants - ESS 2001 Thermal moderator       
  double T=325, tau=80e-6, tau1=400e-6, tau2=12e-6, chi2=2.5, I0=13.5e11, I2=27.6e10, branch1=0.5, branch2=0.5, n2=5, n=20;

  // Branching
  double branch_tail=tau/ESS_SOURCE_DURATION;
  
  // Other variables
  double tail_flag, tau_l;
  
  // Taken directly from the ESS_moderator.comp:
  tail_flag = (rand01()<branch_tail);   /* Choose tail/bulk */
  if (tail_flag)
    {
      if (rand01() < branch2)
	{
	  if (tau1>0)
	    if (rand01() < branch1)     /* Quick and dirty non-general solution */
	      {  /* FIRST CASE a */
		tau_l = tau;
		*p = 1/(branch1*branch2*branch_tail); /* Correct for switching prob. */
	      }
	    else
	      {  /* FIRST CASE b */
		tau_l = tau1;
		*p = 1/((1-branch1)*branch2*branch_tail); /* Correct for switching prob. */
	      }
	  else
	    {
	      tau_l = tau;
	      *p = 1/(branch2*branch_tail); /* Correct for switching prob. */
	    }
	  *t = -tau_l*log(1e-12+rand01());       /* Sample from long-time tail a */
	  /* Correct for true pulse shape */
	  //	  p *= w_focus;                         /* Correct for target focusing */
	  *p *= tau_l/ESS_SOURCE_DURATION;                         /* Correct for tail part */
	  //p *= I0*w_mult*w_geom*Mezei_M(lambda,T);           /* Calculate true intensity */
	  *p *= I0*Mezei_M(lambda,T);
	}
      else
	{
	  /* SECOND CASE */
	  tau_l = tau2*lambda;
	  *t = -tau_l*log(1e-12+rand01());       /* Sample from long-time tail */
	  *p = n2/(n2-1)*((1-exp(-ESS_SOURCE_DURATION/tau_l))-(1-exp(-n2*ESS_SOURCE_DURATION/tau_l))*exp(-(n2-1)*(*t)/tau_l)/n);
	  /* Correct for true pulse shape */
	  *p /= (1-branch2)*branch_tail;          /* Correct for switching prob. */
	  *p *= tau_l/ESS_SOURCE_DURATION;                         /* Correct for tail part */
	  // p *= w_focus;                         /* Correct for target focusing */
	  //p *= I2*w_mult*w_geom/(1+exp(chi2*lambda-2.2))/lambda;                                         /* Calculate true intensity */
	  *p *= I2/(1+exp(chi2*lambda-2.2))/lambda;                                         /* Calculate true intensity */ 
	}
      *t += ESS_SOURCE_DURATION;                                 /* Add pulse length */
    }
  else /* Tail-flag */
    {
      if (tfocus_width>0) {
	*t = tfocus_time-dt;                    /* Set time to hit time window center */
	*t += randpm1()*tfocus_width/2.0;       /* Add random time within window width */
      } else {
	*t = ESS_SOURCE_DURATION*rand01();                        /* Sample from bulk pulse */
      }
      // FLAG to KILL these on return!

      /* if (t<0) ABSORB;                       /\* Kill neutron if outside pulse duration *\/ */
      /* if (t>ESS_SOURCE_DURATION) ABSORB; */
      if (rand01() < branch2)
	{
	  if (rand01() < branch1)     /* Quick and dirty non-general solution */
	    {  /* FIRST CASE a */
	      tau_l = tau;
	      *p = 1/(branch1*branch2*(1-branch_tail)); /* Correct for switching prob. */
	    }
	  else
	    {  /* FIRST CASE b */
	      tau_l = tau1;
	      *p = 1/((1-branch1)*branch2*(1-branch_tail)); /* Correct for switching prob. */
	    }
	  *p *= 1-n/(n-1)*(exp(-*t/tau_l)-exp(-n*(*t)/tau_l)/n); /* Correct for true pulse shape */
	  //	  p *= w_focus;                         /* Correct for target focusing */
	  if (tfocus_width>0) {
	    *p *= tfocus_width/ESS_SOURCE_DURATION;    	  	  /* Correct for time focusing */
	  }
	  //p *= I0*w_mult*w_geom*M(lambda,T);       /* Calculate true intensity */
	  *p *= I0*Mezei_M(lambda,T);       /* Calculate true intensity */
	}
      else
	{
	  /* SECOND CASE */
	  tau_l = tau2*lambda;
	  *p = 1-n2/(n2-1)*(exp(-*t/tau_l)-exp(-n2*(*t)/tau_l)/n2); /* Correct for true pulse shape */
	  *p /= (1-branch2)*(1-branch_tail);   /* Correct for switching prob. */
	  //p *= w_focus;                         /* Correct for target focusing */
	  if (tfocus_width) {
	    *p *= tfocus_width/ESS_SOURCE_DURATION;    		  /* Correct for time focusing */
	  }
	  //p *= I2*w_mult*w_geom/(1+exp(chi2*lambda-2.2))/lambda;    /* Calculate true intensity */
	  *p *= I2/(1+exp(chi2*lambda-2.2))/lambda;    /* Calculate true intensity */
	}
    }

} /* end of ESS_Mezei_thermal */


/* This is the Mezei moderator with a correction term from Klaus Lieutenant */
void ESS_2012_Lieutenant_cold(double *t, double *p, double lambda, double tfocus_w, double tfocus_t, double tfocus_dt, ess_moderator_struct extras)
{
  ESS_Mezei_cold_2012(t, p,  lambda,  tfocus_w,  tfocus_t, tfocus_dt, extras);
  
  double cor;
  
  /* Correction factors to converts 'predicted' spectrum from cold moderator to the one observed in MCNPX */
  if (lambda<=2.5) cor=log(1.402+0.898*lambda)*(2.0776-4.1093*lambda+4.8836*pow(lambda,2)-2.4715*pow(lambda,3)+0.4521*pow(lambda,4));
  else if (lambda <= 3.5) cor = log(1.402 + 0.898*lambda)*(4.3369 - 1.8367*lambda + 0.2524*pow(lambda,2) );
  else if (lambda  > 3.5) cor = log(1.402 + 0.898*lambda);
  
  *p *=cor;

} /* end of ESS_2012_Lieutenant_cold */


/* This is the cold moderator with 2013 updates, fits from Troels Schoenfeldt */
/* Parametrization including moderator height for the "pancake" moderator */
void ESS_2013_Schoenfeldt_cold(double *t, double *p, double lambda, double tfocus_w, double tfocus_t, double tfocus_dt, ess_moderator_struct extras)
{
    /* From the forthcoming Schoenfeldt et al.
       S_cold(\lambda) = (I_1*exp(-\alpha_1*lambda)  + I_2*exp(-\alpha_2*\lambda)) * 1/(1+exp(\alpha_l * (\lambda-lambda_l)))^(1/\gamma)
       + I_SD * (1/lambda) * 1/(1+exp(\alpha_SD*(\lambda-\lambda_SD)))
    */

  /* As function of moderator height, parameters for the brilliance expression */
  double height[7]    = {0.10, 0.05, 0.03, 0.015, 0.01, 0.005, 0.001};
  double I_SD[7]      = {4.75401e+011, 7.0319e+011,  8.36605e+011, 9.41035e+011, 9.54305e+011, 9.83515e+011, 9.54108e+011};
  double alpha_SD[7]  = {0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9};
  double lambda_SD[7] = {2.44444, 2.44444, 2.44444, 2.44444, 2.44444, 2.44444, 2.44444};
  double alpha_l[7]   = {-11.9056, -13.8444, -17.8359, -19.6643, -23.0058, -21.6241, -18.82};
  double lambda_l[7]  = {2.53562, 2.53527, 2.53956, 2.53243, 2.53375, 2.5364, 2.51714};
  double Exponent[7]  = {-0.259162, -0.215819, -0.160541, -0.140769, -0.119278, -0.124298, -0.144056};
  double I_1[7]       = {1.22098e+013, 2.57992e+013, 4.43235e+013, 8.86873e+013, 1.26172e+014, 2.02098e+014, 3.32623e+014};
  double alpha_1[7]   = {0.653579, 0.720244, 0.772538, 0.871765, 0.927905, 1.01579, 1.11621};
  double I_2[7]       = {2.97518e+011, 1.11421e+012, 1.8961e+012,  4.00852e+012, 5.05278e+012, 6.98605e+012, 7.89424e+012};
  double alpha_2[7]   = {0.261097, 0.307898, 0.317865, 0.346354, 0.354282, 0.371298, 0.38382};

  double S_a, S_b;
  double dS;
  int j, idxa, idxb;
  if ((extras.height_c <= height[0]) && (extras.height_c >= height[6])) {
    for (j=0; j<7; j++) {
      if (extras.height_c == height[j]) {
	S_a = ESS_2013_Schoenfeldt_cold_spectrum(I_SD[j], alpha_SD[j], lambda_SD[j], alpha_l[j], lambda_l[j], Exponent[j], I_1[j], alpha_1[j], I_2[j], alpha_2[j], lambda);
	  *p = S_a;
	  break;
	  
	} else if (extras.height_c < height[j] && extras.height_c > height[j+1]) {
	  dS = (height[j]-extras.height_c)/(height[j]-height[j+1]);
	  /* Linear interpolation between the two closest heights */
	  S_a = ESS_2013_Schoenfeldt_cold_spectrum(I_SD[j], alpha_SD[j], lambda_SD[j], alpha_l[j], lambda_l[j], Exponent[j], I_1[j], alpha_1[j], I_2[j], alpha_2[j], lambda);
	  S_b = ESS_2013_Schoenfeldt_cold_spectrum(I_SD[j+1], alpha_SD[j+1], lambda_SD[j+1], alpha_l[j+1], lambda_l[j+1], Exponent[j+1], I_1[j+1], alpha_1[j+1], I_2[j+1], alpha_2[j+1], lambda);
	  *p = (1-dS)*S_a + dS*S_b;
	  break;
      }
    }
  } else {
    printf("Sorry! Moderator height must be between %g and %g m\n",height[6],height[0]);
    exit(-1);
  }

  /* Next is time structure... */
  *t=0;
  *t = extras.tmultiplier*ESS_SOURCE_DURATION*(rand01());
  /* Troels Schoenfeldt function for timestructure */
  *p *= extras.tmultiplier*ESS_2014_Schoenfeldt_cold_timedist(*t, lambda, 100*extras.height_c, ESS_SOURCE_DURATION);
} /* end of ESS_2013_Schoenfeldt_cold */

/* This is the cold moderator with 2014 updates, fits from Troels Schoenfeldt */
/* Parametrization including moderator height for the "pancake" moderator */
void ESS_2014_Schoenfeldt_cold(double *t, double *p, double lambda, double tfocus_w, double tfocus_t, double tfocus_dt, ess_moderator_struct extras)
{
   if ((extras.height_c <= 0.12) && (extras.height_c >= 0.01)) {
    *p = ESS_2014_Schoenfeldt_cold_spectrum(lambda,100*extras.height_c);
  } else {
    printf("Sorry! Moderator height must be between %g and %g m\n",0.12,0.01);
    exit(-1);
  }

  /* Next is time structure... */
  *t=0;
  double pt=0;
  
  // Potentially apply russian roulette technique 
  //while (rand01()>pt) {
  *t = extras.tmultiplier*ESS_SOURCE_DURATION*(rand01());
  /* Troels Schoenfeldt function for timestructure */
  pt = extras.tmultiplier*ESS_2014_Schoenfeldt_cold_timedist(*t, lambda, 100*extras.height_c, ESS_SOURCE_DURATION);
  //}
  
  *p *= pt;
  if (extras.Uniform==0) *p *= ESS_2014_Schoenfeldt_cold_y0(100*extras.Y, 100*extras.height_c) * ESS_2014_Schoenfeldt_cold_x0(-100*extras.X, 100*extras.height_c, 100*extras.Width_c);
} /* end of ESS_2014_Schoenfeldt_cold */

/* This is ESS_2014_Schoenfeldt_cold_y0 - vertical intensity distribution for the 2014 Schoenfeldt cold moderator */
double ESS_2014_Schoenfeldt_cold_y0(double y0,double height){
  
  double one_over_integral_y0_of_height= height/((0.36434*height*height+2.53796*height-0.107774));
  if(y0 < -height/2. || y0 > height/2. )return 0;
  double cosh_ish=(exp(-7e-1/sqrt(height)*(y0-height/2.))+exp(-7e-1/20.*height+7e-1/sqrt(height)*(y0+height/2.)));
  double sinh_ish=(exp(50/sqrt(height)*(y0-height/2.))-1)*(exp(-50/sqrt(height)*(y0+height/2.))-1);
  double tmp=one_over_integral_y0_of_height*cosh_ish*sinh_ish;
  return tmp;
} /* end of ESS_2014_Schoenfeldt_cold_y0 */

/* This is ESS_2014_Schoenfeldt_thermal_y0 - vertical intensity distribution for the 2014 Schoenfeldt cold moderator */
double ESS_2014_Schoenfeldt_thermal_y0(double y0,double height){
  /* Placeholder - we assume that this distribution is flat for now */
  return 1;
} /* end of ESS_2014_Schoenfeldt_thermal_y0 */

/* This is ESS_2014_Schoenfeldt_cold_x0 - horizontal intensity distribution for the 2014 Schoenfeldt cold moderator */
double ESS_2014_Schoenfeldt_cold_x0(double x0,double height, double width){
  double normalization=1;
  if(x0<-width||x0>width)return 0;
  return normalization*(0.008*x0+1)*(exp(height/2.*(x0-width/2))-1)*(exp(-height/2.*(x0+width/2))-1);
} /* end of ESS_2014_Schoenfeldt_cold_x0 */

/* This is ESS_2014_Schoenfeldt_thermal_x0 - horizontal intensity distribution for the 2014 Schoenfeldt cold moderator */
double ESS_2014_Schoenfeldt_thermal_x0(double x0,double height, double width){
  // Kept for reference only...
  /* if(x0>-width&&x0<width)return 0; */
  /* if(x0<0)return fmax(0,2.5*(0.0524986*fabs(x0)-1.84817-0.0189762*height+(-1.49712e+002*exp(-4.06814e-001*height))*exp(-4.48657e-001*fabs(x0)))*(exp(7*(x0+width))-1)); */
  /* return fmax(0,2.5*(0.84199+0.00307022*height)*(0.0524986*fabs(x0)-1.84817-0.0189762*height+(-1.49712e+002*exp(-4.06814e-001*height))*exp(-4.48657e-001*fabs(x0)))*(exp(-7*(x0-width))-1)); */  
  if(x0>-23./2.&&x0<23./2.)return 0;
  long double cosh_ish=fmin(0.0524986*fabs(x0)-1.84817-0.0189762*height+(-1.49712e+002*exp(-4.06814e-001*height))*exp(-4.48657e-001*fabs(x0)),0);
  if(x0<0)return (-1.73518e-003*height*height+2.10277e-002*height+7.65692e-001) // intensity
	    *cosh_ish*(exp(7.*(x0+23./2.))-1); // slope 
  return (-1.73518e-003*height*height+2.10277e-002*height+7.65692e-001) // intensity
    *(0.84199+0.00307022*height) // asumetry
    *cosh_ish*(exp(-7.*(x0-23./2.))-1); // slope
} /* end of ESS_2014_Schoenfeldt_thermal_x0 */

/* This is ESS_2014_Schoenfeldt_cold_Y - vertical intensity distribution for the 2014 Schoenfeldt cold moderator */
double ESS_2014_Schoenfeldt_cold_Y(double Y,double height){
  /* Placeholder - we assume that this distribution is flat for now */
  return 1;
} /* end of ESS_2014_Schoenfeldt_cold_Y */

/* This is ESS_2014_Schoenfeldt_thermal_Y - vertical intensity distribution for the 2014 Schoenfeldt cold moderator */
double ESS_2014_Schoenfeldt_thermal_Y(double Y,double height){
  /* Placeholder - we assume that this distribution is flat for now */
  return 1;
} /* end of ESS_2014_Schoenfeldt_thermal_Y */

/* This is ESS_2014_Schoenfeldt_cold_Theta120 - vertical intensity distribution for the 2014 Schoenfeldt cold moderator */
double ESS_2014_Schoenfeldt_cold_Theta120(double Theta120,double height){
  /* Placeholder - we assume that this distribution is flat for now */
  return 1;
} /* end of ESS_2014_Schoenfeldt_cold_Theta120 */

/* This is ESS_2014_Schoenfeldt_thermal_Theta120 - vertical intensity distribution for the 2014 Schoenfeldt cold moderator */
double ESS_2014_Schoenfeldt_thermal_Theta120(double beamportangle,int isleft){
  if(!isleft)return cos((beamportangle-30)*DEG2RAD)/cos(30*DEG2RAD);
  return cos((90-beamportangle)*DEG2RAD)/cos(30*DEG2RAD);
/* Placeholder - we assume that this distribution is flat for now */
  return 1;
} /* end of ESS_2014_Schoenfeldt_thermal_Theta120 */

/* This is ESS_2014_Schoenfeldt_cold_timedist time-distribution of the 2014 Schoenfeldt cold moderator */ 
double ESS_2014_Schoenfeldt_cold_timedist(double time,double lambda,double height, double pulselength){
        if(time<0)return 0;
        double tau=3.00094e-004*(4.15681e-003*lambda*lambda+2.96212e-001*exp(-1.78408e-001*height)+7.77496e-001)*exp(-6.63537e+001*pow(fmax(1e-13,lambda+.9),-8.64455e+000));
        if(time<pulselength)return ((1-exp(-time/tau)));
        return ((1-exp(-pulselength/tau))*exp(-(time-pulselength)/tau));

} /* end of ESS_2014_Schoenfeldt_cold_timedist */


/* This is ESS_2014_Schoenfeldt_thermal_timedist time-distribution of the 2014 Schoenfeldt cold moderator */    
double ESS_2014_Schoenfeldt_thermal_timedist(double time,double lambda,double height, double pulselength){
        if(time<0)return 0;
        double tau=3.00000e-004*(1.23048e-002*lambda*lambda+1.75628e-001*exp(-1.82452e-001*height)+9.27770e-001)*exp(-3.91090e+001*pow(fmax(1e-13,lambda+9.87990e-001),-7.65675e+000));
        if(time<pulselength)return ((1-exp(-time/tau)));
        return ((1-exp(-pulselength/tau))*exp(-(time-pulselength)/tau));
} /* end of ESS_2014_Schoenfeldt_thermal_timedist */

/* HERE IT IS */


/* This is the thermal moderator with 2013 updates, fits from Troels Schoenfeldt */
void ESS_2013_Schoenfeldt_thermal(double *t, double *p, double lambda, double tfocus_w, double tfocus_t, double tfocus_dt, ess_moderator_struct extras)
{
  /*  From the forthcoming Schoenfeldt et al.
     S_Th(\lambda) = I_Th * exp(k_Th/(T*\lambda^2))*(2*k_Th^2)/(T^2*\lambda^5) + I_SD * \lambda^-1 * 1/(1+exp(\alpha*(\lambda - \lambda_cf))
   */
  
  /* As function of moderator height, parameters for the brilliance expression */
  double height[7]    = {0.10, 0.05, 0.03, 0.015, 0.01, 0.005, 0.001};
  double I_th[7]      = {2.97527e+012, 4.35192e+012, 5.18047e+012, 6.0305e+012,  6.20079e+012, 6.44927e+012, 6.55127e+012};
  double T[7]         = {303.764, 306.099, 307.497, 311.292, 310.525, 310.822, 317.56};
  double I_SD[7]      = {5.38083e+011, 7.3059e+011,  8.94408e+011, 9.89515e+011, 1.02135e+012, 1.07415e+012, 1.12157e+012};
  double alpha[7]     = {2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5};
  double lambda_cf[7] = {0.88, 0.88, 0.88, 0.88, 0.88, 0.88, 0.88};

  double S_a, S_b;
  double dS;
  int j, idxa, idxb;
  if ((extras.height_t <= height[0]) && (extras.height_t >= height[6])) {
    for (j=0; j<6; j++) {
      if (extras.height_c == height[j]) {
	S_a = ESS_2013_Schoenfeldt_thermal_spectrum(I_th[j], T[j], I_SD[j], alpha[j], lambda_cf[j], lambda);
	*p = S_a;
	break;
      } else if (extras.height_t <= height[j] && extras.height_t > height[j+1]) {
	dS = (height[j]-extras.height_t)/(height[j]-height[j+1]);
	/* Linear interpolation between the two closest heights */
	S_a = ESS_2013_Schoenfeldt_thermal_spectrum(I_th[j], T[j], I_SD[j], alpha[j], lambda_cf[j], lambda);
	S_b = ESS_2013_Schoenfeldt_thermal_spectrum(I_th[j+1], T[j+1], I_SD[j+1], alpha[j+1], lambda_cf[j+1], lambda);;
	*p = (1-dS)*S_a + dS*S_b;
	break;
      }
    }
  } else {
    printf("Sorry! Moderator height must be between %g and %g cm\n",height[6],height[0]);
    exit(-1);
  }
  /* Next is time structure... */
  *t=0;
  *t = extras.tmultiplier*ESS_SOURCE_DURATION*(rand01());
  /* Troels Schoenfeldt function for timestructure */
  *p *= extras.tmultiplier*ESS_2014_Schoenfeldt_thermal_timedist(*t, lambda, 100*extras.height_t, ESS_SOURCE_DURATION);
  /* No corrections for "geometrical anisotropy" at this point */
  
} /* end of ESS_2013_Schoenfeldt_thermal */

/* This is the thermal moderator with 2014 updates, fits from Troels Schoenfeldt */
void ESS_2014_Schoenfeldt_thermal(double *t, double *p, double lambda, double tfocus_w, double tfocus_t, double tfocus_dt, ess_moderator_struct extras)
{
  if ((extras.height_t <= 0.12) && (extras.height_t >= 0.01)) {
    *p = ESS_2014_Schoenfeldt_thermal_spectrum(lambda, 100*extras.height_t);
  } else {
    printf("Sorry! Moderator height must be between %g and %g cm\n",0.12,0.01);
    exit(-1);
  }
 
  /* Next is time structure... */
  *t=0;
  *t = extras.tmultiplier*ESS_SOURCE_DURATION*(rand01());
  /* Troels Schoenfeldt function for timestructure */
  *p *= extras.tmultiplier*ESS_2014_Schoenfeldt_thermal_timedist(*t, lambda, 100*extras.height_t, ESS_SOURCE_DURATION);   /* Using Width_c is NOT a typo - dependent on cold moderator geometry - WTF?*/
  *p *= ESS_2014_Schoenfeldt_thermal_y0(100*extras.X, 100*extras.height_t) * ESS_2014_Schoenfeldt_thermal_x0(100*extras.X, 100*extras.height_t, 100*extras.Width_t);
  *p *= ESS_2014_Schoenfeldt_thermal_Theta120(extras.beamportangle,extras.Wasleft);
  //printf("%g\n",ESS_2014_Schoenfeldt_thermal_Theta120(extras.beamportangle,extras.Wasleft));
} /* end of ESS_2014_Schoenfeldt_thermal */

/* 2015 start */

/* This is the thermal moderator with 2015 updates, fits from Troels Schoenfeldt */
void ESS_2015_Schoenfeldt_thermal(double *t, double *p, double lambda, double tfocus_w, double tfocus_t, double tfocus_dt, ess_moderator_struct extras)
{
  if ((extras.height_t == 0.03) || (extras.height_t == 0.06)) {
    *p = ESS_2015_Schoenfeldt_thermal_spectrum(lambda, extras.beamportangle);
  } else {
    printf("Sorry! Moderator height must be either %g or %g m\n",0.03,0.06);
    exit(-1);
  }

 
  /* Next is time structure... */
  *t=0;
  *t = extras.tmultiplier*ESS_SOURCE_DURATION*(rand01());
  /* Troels Schoenfeldt function for timestructure */
  *p *= extras.tmultiplier*ESS_2015_Schoenfeldt_thermal_timedist(*t, lambda, 3 /* cm height */, ESS_SOURCE_DURATION);  
  if (extras.height_c == 0.03) {
    // 3cm case
    *p *= ESS_2015_Schoenfeldt_thermal_y0(100*extras.Y) * ESS_2015_Schoenfeldt_thermal_x0(100*extras.X, extras.beamportangle);
  } else {
    // 6cm case
    // Downscale brightness by factor from 
    // "New ESS Moderator Baseline", Ken Andersen, 9/4/2015
    *p *= (6.2e14/9.0e14);
    *p *= ESS_2014_Schoenfeldt_thermal_y0(100*extras.Y, 100*extras.height_c) * ESS_2015_Schoenfeldt_thermal_x0(100*extras.X, extras.beamportangle);
  }
} /* end of ESS_2015_Schoenfeldt_thermal */


/* This is the cold moderator with 2015 updates, fits from Troels Schoenfeldt */
/* Parametrization including moderator height for the "pancake" moderator */
void ESS_2015_Schoenfeldt_cold(double *t, double *p, double lambda, double tfocus_w, double tfocus_t, double tfocus_dt, ess_moderator_struct extras)
{
   if ((extras.height_c == 0.03) || (extras.height_c == 0.06)) {
    *p = ESS_2015_Schoenfeldt_cold_spectrum(lambda,extras.beamportangle);
  } else {
    printf("Sorry! Moderator height must be either %g or %g m\n",0.03,0.06);
    exit(-1);
  }

  /* Next is time structure... */
  *t=0;
  double pt=0;
  
  // Potentially apply russian roulette technique 
  //while (rand01()>pt) {
  *t = extras.tmultiplier*ESS_SOURCE_DURATION*(rand01());
  /* Troels Schoenfeldt function for timestructure */
  pt = extras.tmultiplier*ESS_2015_Schoenfeldt_cold_timedist(*t, lambda, 3 /* cm height */, ESS_SOURCE_DURATION);
  //}
  
  *p *= pt;
  if (extras.Uniform==0) {
    if (extras.height_c == 0.03) {
      // 3cm case
      *p *= ESS_2015_Schoenfeldt_cold_y0(100*extras.Y) * ESS_2015_Schoenfeldt_cold_x0(100*extras.X, extras.beamportangle);
    } else {
      // 6cm case
      // Downscale brightness by factor from 
      // "New ESS Moderator Baseline", Ken Andersen, 9/4/2015
      *p *= (10.1e14/16.0e14);
      *p *= ESS_2014_Schoenfeldt_cold_y0(100*extras.Y, 100*extras.height_c) * ESS_2015_Schoenfeldt_cold_x0(100*extras.X, extras.beamportangle);
    }
  }
} /* end of ESS_2015_Schoenfeldt_cold */

/* This is ESS_2015_Schoenfeldt_cold_y0 - vertical intensity distribution for the 2015 Schoenfeldt cold moderator */
double ESS_2015_Schoenfeldt_cold_y0(double y0){
    double par3=30;
    double par4=.35;
    long double cosh_ish=exp(-par4*y0)+exp(par4*y0);
    long double sinh_ish=pow(1+exp(par3*(y0-3./2.)),-1)*pow(1+exp(-par3*(y0+3./2.)),-1);
    return 1./2.*(double)((long double)cosh_ish*(long double)sinh_ish);

} /* end of ESS_2015_Schoenfeldt_cold_y0 */

/* This is ESS_2015_Schoenfeldt_thermal_y0 - vertical intensity distribution for the 2015 Schoenfeldt cold moderator */
double ESS_2015_Schoenfeldt_thermal_y0(double y0){
    if(y0<-3./2.+0.105){
        return 1.005*exp(-pow((y0+3./2.-0.105)/0.372,2));
    } else if(y0>3./2.-0.105){
        return 1.005*exp(-pow((y0-3./2.+0.105)/0.372,2));
    }
    return 1.005;
} /* end of ESS_2015_Schoenfeldt_thermal_y0 */

/* This is ESS_2015_Schoenfeldt_cold_x0 - horizontal intensity distribution for the 2015 Schoenfeldt cold moderator */
double ESS_2015_Schoenfeldt_cold_x0(double x0,double theta){
  // GEOMETRY / SAMPLING SPACE
    double i=(theta-5.)/10.;
    double par0=0.0146115+0.00797729*i-0.00279541*i*i;
    double par1=0.980886;
    if(i==1)par1=0.974217;
    if(i==2)par1=0.981462;
    if(i==3)par1=1.01466;
    if(i==4)par1=1.11707;
    if(i==5)par1=1.16057;
        
    double par2=-4-.75*i;
    if(i==0)par2=-20;
    double par3=-14.9402-0.178369*i+0.0367007*i*i;
    if(i==0)par3=-14.27;
    double par4=-15;
    if(i==3)par4=-3.5;
    if(i==5)par4=-1.9;
    double par5=-7.07979+0.0835695*i-0.0546662*i*i;
    if(i==4)par5=-8.1;

    double line=par0*(x0+12)+par1;
    double CutLeftCutRight=1./((1+exp(par2*(x0-par3)))*(1+exp(-par4*(x0-par5))));

    return line*CutLeftCutRight;
} /* end of ESS_2015_Schoenfeldt_cold_x0 */

/* This is ESS_2015_Schoenfeldt_thermal_x0 - horizontal intensity distribution for the 2015 Schoenfeldt cold moderator */
double ESS_2015_Schoenfeldt_thermal_x0(double x0,double theta){
    double i=(theta-5.)/10.;
    double par0=-5.54775+0.492804*i;
    double par1=-0.265929-0.711477*i;
    if(theta==55)par1=-2.55;

    double par2=0.821885+0.00914832*i;
    double par3=1.31108-0.00698647*i;
    if(theta==55)par3=1.23;
    double par4=-.035;
    double par5=-0.0817358+0.00807125*i;
        
    double par6=-8;
    double par7=-7.15;
    if(theta==45)par7=-8.2;
    if(theta==55)par7=-7.7;

    double par8=-8;
    double par9=7.15;
    if(theta==45)par9=7.5;
    if(theta==55)par9=8.2;

    double soften1=1./(1+exp(8.*(x0-par0)));
    double soften2=1./(1+exp(8.*(x0-par1)));
    double CutLeftCutRight=1./((1+exp(par6*(x0-par7)))*(1+exp(-par8*(x0-par9))));
    double line1=par4*(x0-par0)+par2;
    double line2=(par2-par3)/(par0-par1)*(x0-par0)+par2;
    double line3=par5*(x0-par1)+par3;
    double add45degbumb=1.2*exp(-(x0+7.55)*(x0+7.55)/.35/.35);


    return CutLeftCutRight*(
        (line1+add45degbumb)*soften1
        +line2*soften2*(1-soften1)
        +line3*(1-soften2)
        );
} /* end of ESS_2015_Schoenfeldt_thermal_x0 */

/* This is ESS_2015_Schoenfeldt_cold_Y - vertical intensity distribution for the 2015 Schoenfeldt cold moderator */
double ESS_2015_Schoenfeldt_cold_Y(double Y,double height){
  /* Placeholder - we assume that this distribution is flat for now */
  return 1;
} /* end of ESS_2015_Schoenfeldt_cold_Y */

/* This is ESS_2015_Schoenfeldt_thermal_Y - vertical intensity distribution for the 2015 Schoenfeldt cold moderator */
double ESS_2015_Schoenfeldt_thermal_Y(double Y,double height){
  /* Placeholder - we assume that this distribution is flat for now */
  return 1;
} /* end of ESS_2015_Schoenfeldt_thermal_Y */

/* This is ESS_2015_Schoenfeldt_cold_Theta120 - vertical intensity distribution for the 2015 Schoenfeldt cold moderator */
double ESS_2015_Schoenfeldt_cold_Theta120(double Theta120,double height){
  /* Placeholder - we assume that this distribution is flat for now */
  return 1;
} /* end of ESS_2015_Schoenfeldt_cold_Theta120 */

/* This is ESS_2015_Schoenfeldt_thermal_Theta120 - vertical intensity distribution for the 2015 Schoenfeldt cold moderator */
double ESS_2015_Schoenfeldt_thermal_Theta120(double beamportangle,int isleft){
  if(!isleft)return cos((beamportangle-30)*DEG2RAD)/cos(30*DEG2RAD);
  return cos((90-beamportangle)*DEG2RAD)/cos(30*DEG2RAD);
/* Placeholder - we assume that this distribution is flat for now */
  return 1;
} /* end of ESS_2015_Schoenfeldt_thermal_Theta120 */


/* This is ESS_2015_Schoenfeldt_cold_timedist time-distribution of the 2014 Schoenfeldt cold moderator */ 
double ESS_2015_Schoenfeldt_cold_timedist(double time,double lambda,double height, double pulselength){
        if(time<0)return 0;
        double tau=3.00094e-004*(4.15681e-003*lambda*lambda+2.96212e-001*exp(-1.78408e-001*height)+7.77496e-001)*exp(-6.63537e+001*pow(fmax(1e-13,lambda+.9),-8.64455e+000));
        if(time<pulselength)return ((1-exp(-time/tau)));
        return ((1-exp(-pulselength/tau))*exp(-(time-pulselength)/tau));

} /* end of ESS_2015_Schoenfeldt_cold_timedist */

/* This is ESS_2015_Schoenfeldt_thermal_timedist time-distribution of the 2015 Schoenfeldt cold moderator */    
double ESS_2015_Schoenfeldt_thermal_timedist(double time,double lambda,double height, double pulselength){
        if(time<0)return 0;
        double tau=3.00000e-004*(1.23048e-002*lambda*lambda+1.75628e-001*exp(-1.82452e-001*height)+9.27770e-001)*exp(-3.91090e+001*pow(fmax(1e-13,lambda+9.87990e-001),-7.65675e+000));
        if(time<pulselength)return ((1-exp(-time/tau)));
        return ((1-exp(-pulselength/tau))*exp(-(time-pulselength)/tau));
} /* end of ESS_2015_Schoenfeldt_thermal_timedist */

/*2015 end */

double TSC2015_z0_BF3cm(const double x0){
    if(x0<-7.16)
        return (8.27-5.1)/(-7.16+14.2)*(x0+14.2)+5.1;
    return 8.27;
}

/* Display of geometry - flat and TDR-like */
void ESS_mcdisplay_flat(double geometry)
{
}/* end of ESS_mcdisplay_flat */

void ESS_mcdisplay_TDRlike(double geometry)
{
}/* end of ESS_mcdisplay_TDRlike */


/* end of ess_source-lib.c */



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

/* component Source=ESS_moderator() [2] DECLARE */
/* Parameter definition for component type 'ESS_moderator' */
struct _struct_ESS_moderator_parameters {
  /* Component type 'ESS_moderator' setting parameters */
  MCNUM isleft;
  MCNUM Lmin;
  MCNUM Lmax;
  MCNUM cold_frac;
  MCNUM dist;
  MCNUM focus_xw;
  MCNUM focus_yh;
  long target_index;
  MCNUM tmax_multiplier;
  MCNUM yheight_c;
  MCNUM yheight_t;
  long n_pulses;
  MCNUM acc_power;
  MCNUM beamport_angle;
  char sourcedef[16384];
  MCNUM xwidth_c;
  MCNUM xwidth_t;
  MCNUM extraction_opening;
  /* Component type 'ESS_moderator' private parameters */
  /* Component type 'ESS_moderator' DECLARE code stored as structure members */
                             
double l_range, w_mult, w_geom, w_geom_c, w_geom_t, w_stat;
                                                                 
int cold;
                                                                  
int wasleft;
                                                                              
int cosine;
                                
  double tx,ty,tz;
                                                                                             
  double cx,cy,cz,cyl_radius;
                                                  
  double t1x,t1y,t1z,t2x,t2y,t2z;

                                                    
  ess_moderator_struct modextras;
                                          
  functype cold_bril;
  functype thermal_bril;
                                                   
  double tfocus_width,  tfocus_time,  dt;
                                                   
  double r_empty;
  double lambda,xprime,yprime,zprime;
  double abs_beamport_angle,cos_beamport_angle,sin_beamport_angle;
                                                                 
  double cos_thermal,cos_cold, edge_thermal;
  int sgnport;
}; /* _struct_ESS_moderator_parameters */
typedef struct _struct_ESS_moderator_parameters _class_ESS_moderator_parameters;

/* Parameters for component type 'ESS_moderator' */
struct _struct_ESS_moderator {
  char     _name[256]; /* e.g. Source */
  char     _type[256]; /* ESS_moderator */
  long     _index; /* e.g. 2 index in TRACE list */
  Coords   _position_absolute;
  Coords   _position_relative; /* wrt PREVIOUS */
  Rotation _rotation_absolute;
  Rotation _rotation_relative; /* wrt PREVIOUS */
  int      _rotation_is_identity;
  _class_ESS_moderator_parameters _parameters;
};
typedef struct _struct_ESS_moderator _class_ESS_moderator;
_class_ESS_moderator _Source_var;
_class_ESS_moderator *_Source = &_Source_var;
#pragma acc declare create ( _Source_var )
#pragma acc declare create ( _Source )

/* component PortOrig=Arm() [3] DECLARE */
/* Parameter definition for component type 'Arm' */
struct _struct_Arm_parameters {
  char Arm_has_no_parameters;
}; /* _struct_Arm_parameters */
typedef struct _struct_Arm_parameters _class_Arm_parameters;

/* Parameters for component type 'Arm' */
struct _struct_Arm {
  char     _name[256]; /* e.g. PortOrig */
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
_class_Arm _PortOrig_var;
_class_Arm *_PortOrig = &_PortOrig_var;
#pragma acc declare create ( _PortOrig_var )
#pragma acc declare create ( _PortOrig )

/* component XmonT=PSDlin_monitor() [4] DECLARE */
/* Parameter definition for component type 'PSDlin_monitor' */
struct _struct_PSDlin_monitor_parameters {
  /* Component type 'PSDlin_monitor' setting parameters */
  MCNUM nx;
  char filename[16384];
  MCNUM xmin;
  MCNUM xmax;
  MCNUM ymin;
  MCNUM ymax;
  MCNUM xwidth;
  MCNUM yheight;
  MCNUM restore_neutron;
  /* Component type 'PSDlin_monitor' private parameters */
  /* Component type 'PSDlin_monitor' DECLARE code stored as structure members */
  DArray1d PSDlin_N;
  DArray1d PSDlin_p;
  DArray1d PSDlin_p2;
}; /* _struct_PSDlin_monitor_parameters */
typedef struct _struct_PSDlin_monitor_parameters _class_PSDlin_monitor_parameters;

/* Parameters for component type 'PSDlin_monitor' */
struct _struct_PSDlin_monitor {
  char     _name[256]; /* e.g. XmonT */
  char     _type[256]; /* PSDlin_monitor */
  long     _index; /* e.g. 2 index in TRACE list */
  Coords   _position_absolute;
  Coords   _position_relative; /* wrt PREVIOUS */
  Rotation _rotation_absolute;
  Rotation _rotation_relative; /* wrt PREVIOUS */
  int      _rotation_is_identity;
  _class_PSDlin_monitor_parameters _parameters;
};
typedef struct _struct_PSDlin_monitor _class_PSDlin_monitor;
_class_PSDlin_monitor _XmonT_var;
_class_PSDlin_monitor *_XmonT = &_XmonT_var;
#pragma acc declare create ( _XmonT_var )
#pragma acc declare create ( _XmonT )

_class_PSDlin_monitor _XmonC_var;
_class_PSDlin_monitor *_XmonC = &_XmonC_var;
#pragma acc declare create ( _XmonC_var )
#pragma acc declare create ( _XmonC )

/* component tofT=TOF_monitor() [6] DECLARE */
/* Parameter definition for component type 'TOF_monitor' */
struct _struct_TOF_monitor_parameters {
  /* Component type 'TOF_monitor' setting parameters */
  MCNUM nt;
  char filename[16384];
  MCNUM xmin;
  MCNUM xmax;
  MCNUM ymin;
  MCNUM ymax;
  MCNUM xwidth;
  MCNUM yheight;
  MCNUM tmin;
  MCNUM tmax;
  MCNUM dt;
  MCNUM restore_neutron;
  /* Component type 'TOF_monitor' private parameters */
  /* Component type 'TOF_monitor' DECLARE code stored as structure members */
  DArray1d TOF_N;
  DArray1d TOF_p;
  DArray1d TOF_p2;
  double t_min, t_max, delta_t;
}; /* _struct_TOF_monitor_parameters */
typedef struct _struct_TOF_monitor_parameters _class_TOF_monitor_parameters;

/* Parameters for component type 'TOF_monitor' */
struct _struct_TOF_monitor {
  char     _name[256]; /* e.g. tofT */
  char     _type[256]; /* TOF_monitor */
  long     _index; /* e.g. 2 index in TRACE list */
  Coords   _position_absolute;
  Coords   _position_relative; /* wrt PREVIOUS */
  Rotation _rotation_absolute;
  Rotation _rotation_relative; /* wrt PREVIOUS */
  int      _rotation_is_identity;
  _class_TOF_monitor_parameters _parameters;
};
typedef struct _struct_TOF_monitor _class_TOF_monitor;
_class_TOF_monitor _tofT_var;
_class_TOF_monitor *_tofT = &_tofT_var;
#pragma acc declare create ( _tofT_var )
#pragma acc declare create ( _tofT )

_class_TOF_monitor _tofC_var;
_class_TOF_monitor *_tofC = &_tofC_var;
#pragma acc declare create ( _tofC_var )
#pragma acc declare create ( _tofC )

/* component Slit1=Slit() [8] DECLARE */
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
  char     _name[256]; /* e.g. Slit1 */
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
_class_Slit _Slit1_var;
_class_Slit *_Slit1 = &_Slit1_var;
#pragma acc declare create ( _Slit1_var )
#pragma acc declare create ( _Slit1 )

/* component LamT=L_monitor() [9] DECLARE */
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
  char     _name[256]; /* e.g. LamT */
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
_class_L_monitor _LamT_var;
_class_L_monitor *_LamT = &_LamT_var;
#pragma acc declare create ( _LamT_var )
#pragma acc declare create ( _LamT )

_class_L_monitor _LamC_var;
_class_L_monitor *_LamC = &_LamC_var;
#pragma acc declare create ( _LamC_var )
#pragma acc declare create ( _LamC )

int mcNUMCOMP = 10;

/* User declarations from instrument definition. Can define functions. */
  //double lambdamin,lambdamax;
  int IsCold,np,sgn;
  char srcdef[128];
  double Yheight,power;
  double cosa,sina;
  double ctrX,ctrZ,BX0,BZ0;
  double mon_width;

#undef compcurname
#undef compcurtype
#undef compcurindex
/* end of instrument 'ESS_2015_test' and components DECLARE */

/* *****************************************************************************
* instrument 'ESS_2015_test' and components INITIALISE
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

/* component Source=ESS_moderator() SETTING, POSITION/ROTATION */
int _Source_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Source_setpos] component Source=ESS_moderator() SETTING [/usr/share/mcstas/3.0-dev/sources/ESS_moderator.comp:128]");
  stracpy(_Source->_name, "Source", 16384);
  stracpy(_Source->_type, "ESS_moderator", 16384);
  _Source->_index=2;
  _Source->_parameters.isleft = 0.9;
  #define isleft (_Source->_parameters.isleft)
  _Source->_parameters.Lmin = instrument->_parameters._lmin;
  #define Lmin (_Source->_parameters.Lmin)
  _Source->_parameters.Lmax = instrument->_parameters._lmax;
  #define Lmax (_Source->_parameters.Lmax)
  _Source->_parameters.cold_frac = instrument->_parameters._frac;
  #define cold_frac (_Source->_parameters.cold_frac)
  _Source->_parameters.dist = 0;
  #define dist (_Source->_parameters.dist)
  _Source->_parameters.focus_xw = 0.01;
  #define focus_xw (_Source->_parameters.focus_xw)
  _Source->_parameters.focus_yh = 0.01;
  #define focus_yh (_Source->_parameters.focus_yh)
  _Source->_parameters.target_index = + 6;
  #define target_index (_Source->_parameters.target_index)
  _Source->_parameters.tmax_multiplier = 3;
  #define tmax_multiplier (_Source->_parameters.tmax_multiplier)
  _Source->_parameters.yheight_c = Yheight;
  #define yheight_c (_Source->_parameters.yheight_c)
  _Source->_parameters.yheight_t = Yheight;
  #define yheight_t (_Source->_parameters.yheight_t)
  _Source->_parameters.n_pulses = np;
  #define n_pulses (_Source->_parameters.n_pulses)
  _Source->_parameters.acc_power = power;
  #define acc_power (_Source->_parameters.acc_power)
  _Source->_parameters.beamport_angle = instrument->_parameters._ANGLE;
  #define beamport_angle (_Source->_parameters.beamport_angle)
  if(srcdef && strlen(srcdef))
    stracpy(_Source->_parameters.sourcedef, srcdef ? srcdef : "", 16384);
  else 
  _Source->_parameters.sourcedef[0]='\0';
  #define sourcedef (_Source->_parameters.sourcedef)
  _Source->_parameters.xwidth_c = 0.1;
  #define xwidth_c (_Source->_parameters.xwidth_c)
  _Source->_parameters.xwidth_t = 0.18;
  #define xwidth_t (_Source->_parameters.xwidth_t)
  _Source->_parameters.extraction_opening = 120;
  #define extraction_opening (_Source->_parameters.extraction_opening)

  #define l_range (_Source->_parameters.l_range)
  #define w_mult (_Source->_parameters.w_mult)
  #define w_geom (_Source->_parameters.w_geom)
  #define w_geom_c (_Source->_parameters.w_geom_c)
  #define w_geom_t (_Source->_parameters.w_geom_t)
  #define w_stat (_Source->_parameters.w_stat)
  #define cold (_Source->_parameters.cold)
  #define wasleft (_Source->_parameters.wasleft)
  #define cosine (_Source->_parameters.cosine)
  #define tx (_Source->_parameters.tx)
  #define ty (_Source->_parameters.ty)
  #define tz (_Source->_parameters.tz)
  #define cx (_Source->_parameters.cx)
  #define cy (_Source->_parameters.cy)
  #define cz (_Source->_parameters.cz)
  #define cyl_radius (_Source->_parameters.cyl_radius)
  #define t1x (_Source->_parameters.t1x)
  #define t1y (_Source->_parameters.t1y)
  #define t1z (_Source->_parameters.t1z)
  #define t2x (_Source->_parameters.t2x)
  #define t2y (_Source->_parameters.t2y)
  #define t2z (_Source->_parameters.t2z)
  #define modextras (_Source->_parameters.modextras)
  #define cold_bril (_Source->_parameters.cold_bril)
  #define thermal_bril (_Source->_parameters.thermal_bril)
  #define tfocus_width (_Source->_parameters.tfocus_width)
  #define tfocus_time (_Source->_parameters.tfocus_time)
  #define dt (_Source->_parameters.dt)
  #define r_empty (_Source->_parameters.r_empty)
  #define lambda (_Source->_parameters.lambda)
  #define xprime (_Source->_parameters.xprime)
  #define yprime (_Source->_parameters.yprime)
  #define zprime (_Source->_parameters.zprime)
  #define sgnport (_Source->_parameters.sgnport)
  #define abs_beamport_angle (_Source->_parameters.abs_beamport_angle)
  #define cos_beamport_angle (_Source->_parameters.cos_beamport_angle)
  #define edge_thermal (_Source->_parameters.edge_thermal)
  #define cos_thermal (_Source->_parameters.cos_thermal)
  #define cos_cold (_Source->_parameters.cos_cold)
  #define sin_beamport_angle (_Source->_parameters.sin_beamport_angle)

  #undef isleft
  #undef Lmin
  #undef Lmax
  #undef cold_frac
  #undef dist
  #undef focus_xw
  #undef focus_yh
  #undef target_index
  #undef tmax_multiplier
  #undef yheight_c
  #undef yheight_t
  #undef n_pulses
  #undef acc_power
  #undef beamport_angle
  #undef sourcedef
  #undef xwidth_c
  #undef xwidth_t
  #undef extraction_opening
  #undef l_range
  #undef w_mult
  #undef w_geom
  #undef w_geom_c
  #undef w_geom_t
  #undef w_stat
  #undef cold
  #undef wasleft
  #undef cosine
  #undef tx
  #undef ty
  #undef tz
  #undef cx
  #undef cy
  #undef cz
  #undef cyl_radius
  #undef t1x
  #undef t1y
  #undef t1z
  #undef t2x
  #undef t2y
  #undef t2z
  #undef modextras
  #undef cold_bril
  #undef thermal_bril
  #undef tfocus_width
  #undef tfocus_time
  #undef dt
  #undef r_empty
  #undef lambda
  #undef xprime
  #undef yprime
  #undef zprime
  #undef sgnport
  #undef abs_beamport_angle
  #undef cos_beamport_angle
  #undef edge_thermal
  #undef cos_thermal
  #undef cos_cold
  #undef sin_beamport_angle
  /* component Source=ESS_moderator() AT ROTATED */
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
  } /* Source=ESS_moderator() AT ROTATED */
  DEBUG_COMPONENT("Source", _Source->_position_absolute, _Source->_rotation_absolute);
  instrument->_position_absolute[2] = _Source->_position_absolute;
  instrument->_position_relative[2] = _Source->_position_relative;
  instrument->counter_N[2]  = instrument->counter_P[2] = instrument->counter_P2[2] = 0;
  instrument->counter_AbsorbProp[2]= 0;
  return(0);
} /* _Source_setpos */

/* component PortOrig=Arm() SETTING, POSITION/ROTATION */
int _PortOrig_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_PortOrig_setpos] component PortOrig=Arm() SETTING [Arm:0]");
  stracpy(_PortOrig->_name, "PortOrig", 16384);
  stracpy(_PortOrig->_type, "Arm", 16384);
  _PortOrig->_index=3;
  /* component PortOrig=Arm() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (0)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _Origin->_rotation_absolute, _PortOrig->_rotation_absolute);
    rot_transpose(_Source->_rotation_absolute, tr1);
    rot_mul(_PortOrig->_rotation_absolute, tr1, _PortOrig->_rotation_relative);
    _PortOrig->_rotation_is_identity =  rot_test_identity(_PortOrig->_rotation_relative);
    tc1 = coords_set(
      BX0, 0, BZ0);
    rot_transpose(_Origin->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _PortOrig->_position_absolute = coords_add(_Origin->_position_absolute, tc2);
    tc1 = coords_sub(_Source->_position_absolute, _PortOrig->_position_absolute);
    _PortOrig->_position_relative = rot_apply(_PortOrig->_rotation_absolute, tc1);
  } /* PortOrig=Arm() AT ROTATED */
  DEBUG_COMPONENT("PortOrig", _PortOrig->_position_absolute, _PortOrig->_rotation_absolute);
  instrument->_position_absolute[3] = _PortOrig->_position_absolute;
  instrument->_position_relative[3] = _PortOrig->_position_relative;
  instrument->counter_N[3]  = instrument->counter_P[3] = instrument->counter_P2[3] = 0;
  instrument->counter_AbsorbProp[3]= 0;
  return(0);
} /* _PortOrig_setpos */

/* component XmonT=PSDlin_monitor() SETTING, POSITION/ROTATION */
int _XmonT_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_XmonT_setpos] component XmonT=PSDlin_monitor() SETTING [/usr/share/mcstas/3.0-dev/monitors/PSDlin_monitor.comp:58]");
  stracpy(_XmonT->_name, "XmonT", 16384);
  stracpy(_XmonT->_type, "PSDlin_monitor", 16384);
  _XmonT->_index=4;
  _XmonT->_parameters.nx = 250;
  #define nx (_XmonT->_parameters.nx)
  if("XmonT.dat" && strlen("XmonT.dat"))
    stracpy(_XmonT->_parameters.filename, "XmonT.dat" ? "XmonT.dat" : "", 16384);
  else 
  _XmonT->_parameters.filename[0]='\0';
  #define filename (_XmonT->_parameters.filename)
  _XmonT->_parameters.xmin = - mon_width / 2 + instrument->_parameters._mon_shift;
  #define xmin (_XmonT->_parameters.xmin)
  _XmonT->_parameters.xmax = mon_width / 2 + instrument->_parameters._mon_shift;
  #define xmax (_XmonT->_parameters.xmax)
  _XmonT->_parameters.ymin = -0.05;
  #define ymin (_XmonT->_parameters.ymin)
  _XmonT->_parameters.ymax = 0.05;
  #define ymax (_XmonT->_parameters.ymax)
  _XmonT->_parameters.xwidth = 0;
  #define xwidth (_XmonT->_parameters.xwidth)
  _XmonT->_parameters.yheight = 0.04;
  #define yheight (_XmonT->_parameters.yheight)
  _XmonT->_parameters.restore_neutron = 1;
  #define restore_neutron (_XmonT->_parameters.restore_neutron)

  #define PSDlin_N (_XmonT->_parameters.PSDlin_N)
  #define PSDlin_p (_XmonT->_parameters.PSDlin_p)
  #define PSDlin_p2 (_XmonT->_parameters.PSDlin_p2)

  #undef nx
  #undef filename
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef xwidth
  #undef yheight
  #undef restore_neutron
  #undef PSDlin_N
  #undef PSDlin_p
  #undef PSDlin_p2
  /* component XmonT=PSDlin_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _PortOrig->_rotation_absolute, _XmonT->_rotation_absolute);
    rot_transpose(_PortOrig->_rotation_absolute, tr1);
    rot_mul(_XmonT->_rotation_absolute, tr1, _XmonT->_rotation_relative);
    _XmonT->_rotation_is_identity =  rot_test_identity(_XmonT->_rotation_relative);
    tc1 = coords_set(
      0, 0, 0.3);
    rot_transpose(_PortOrig->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _XmonT->_position_absolute = coords_add(_PortOrig->_position_absolute, tc2);
    tc1 = coords_sub(_PortOrig->_position_absolute, _XmonT->_position_absolute);
    _XmonT->_position_relative = rot_apply(_XmonT->_rotation_absolute, tc1);
  } /* XmonT=PSDlin_monitor() AT ROTATED */
  DEBUG_COMPONENT("XmonT", _XmonT->_position_absolute, _XmonT->_rotation_absolute);
  instrument->_position_absolute[4] = _XmonT->_position_absolute;
  instrument->_position_relative[4] = _XmonT->_position_relative;
  instrument->counter_N[4]  = instrument->counter_P[4] = instrument->counter_P2[4] = 0;
  instrument->counter_AbsorbProp[4]= 0;
  return(0);
} /* _XmonT_setpos */

/* component XmonC=PSDlin_monitor() SETTING, POSITION/ROTATION */
int _XmonC_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_XmonC_setpos] component XmonC=PSDlin_monitor() SETTING [/usr/share/mcstas/3.0-dev/monitors/PSDlin_monitor.comp:58]");
  stracpy(_XmonC->_name, "XmonC", 16384);
  stracpy(_XmonC->_type, "PSDlin_monitor", 16384);
  _XmonC->_index=5;
  _XmonC->_parameters.nx = 250;
  #define nx (_XmonC->_parameters.nx)
  if("XmonC.dat" && strlen("XmonC.dat"))
    stracpy(_XmonC->_parameters.filename, "XmonC.dat" ? "XmonC.dat" : "", 16384);
  else 
  _XmonC->_parameters.filename[0]='\0';
  #define filename (_XmonC->_parameters.filename)
  _XmonC->_parameters.xmin = - mon_width / 2 + instrument->_parameters._mon_shift;
  #define xmin (_XmonC->_parameters.xmin)
  _XmonC->_parameters.xmax = mon_width / 2 + instrument->_parameters._mon_shift;
  #define xmax (_XmonC->_parameters.xmax)
  _XmonC->_parameters.ymin = -0.05;
  #define ymin (_XmonC->_parameters.ymin)
  _XmonC->_parameters.ymax = 0.05;
  #define ymax (_XmonC->_parameters.ymax)
  _XmonC->_parameters.xwidth = 0;
  #define xwidth (_XmonC->_parameters.xwidth)
  _XmonC->_parameters.yheight = 0.04;
  #define yheight (_XmonC->_parameters.yheight)
  _XmonC->_parameters.restore_neutron = 1;
  #define restore_neutron (_XmonC->_parameters.restore_neutron)

  #define PSDlin_N (_XmonC->_parameters.PSDlin_N)
  #define PSDlin_p (_XmonC->_parameters.PSDlin_p)
  #define PSDlin_p2 (_XmonC->_parameters.PSDlin_p2)

  #undef nx
  #undef filename
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef xwidth
  #undef yheight
  #undef restore_neutron
  #undef PSDlin_N
  #undef PSDlin_p
  #undef PSDlin_p2
  /* component XmonC=PSDlin_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _PortOrig->_rotation_absolute, _XmonC->_rotation_absolute);
    rot_transpose(_XmonT->_rotation_absolute, tr1);
    rot_mul(_XmonC->_rotation_absolute, tr1, _XmonC->_rotation_relative);
    _XmonC->_rotation_is_identity =  rot_test_identity(_XmonC->_rotation_relative);
    tc1 = coords_set(
      0, 0, 0.3);
    rot_transpose(_PortOrig->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _XmonC->_position_absolute = coords_add(_PortOrig->_position_absolute, tc2);
    tc1 = coords_sub(_XmonT->_position_absolute, _XmonC->_position_absolute);
    _XmonC->_position_relative = rot_apply(_XmonC->_rotation_absolute, tc1);
  } /* XmonC=PSDlin_monitor() AT ROTATED */
  DEBUG_COMPONENT("XmonC", _XmonC->_position_absolute, _XmonC->_rotation_absolute);
  instrument->_position_absolute[5] = _XmonC->_position_absolute;
  instrument->_position_relative[5] = _XmonC->_position_relative;
  instrument->counter_N[5]  = instrument->counter_P[5] = instrument->counter_P2[5] = 0;
  instrument->counter_AbsorbProp[5]= 0;
  return(0);
} /* _XmonC_setpos */

/* component tofT=TOF_monitor() SETTING, POSITION/ROTATION */
int _tofT_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_tofT_setpos] component tofT=TOF_monitor() SETTING [/usr/share/mcstas/3.0-dev/monitors/TOF_monitor.comp:63]");
  stracpy(_tofT->_name, "tofT", 16384);
  stracpy(_tofT->_type, "TOF_monitor", 16384);
  _tofT->_index=6;
  _tofT->_parameters.nt = 100;
  #define nt (_tofT->_parameters.nt)
  if("tofT.dat" && strlen("tofT.dat"))
    stracpy(_tofT->_parameters.filename, "tofT.dat" ? "tofT.dat" : "", 16384);
  else 
  _tofT->_parameters.filename[0]='\0';
  #define filename (_tofT->_parameters.filename)
  _tofT->_parameters.xmin = - mon_width / 2 + instrument->_parameters._mon_shift;
  #define xmin (_tofT->_parameters.xmin)
  _tofT->_parameters.xmax = mon_width / 2 + instrument->_parameters._mon_shift;
  #define xmax (_tofT->_parameters.xmax)
  _tofT->_parameters.ymin = -0.05;
  #define ymin (_tofT->_parameters.ymin)
  _tofT->_parameters.ymax = 0.05;
  #define ymax (_tofT->_parameters.ymax)
  _tofT->_parameters.xwidth = 0;
  #define xwidth (_tofT->_parameters.xwidth)
  _tofT->_parameters.yheight = 0.04;
  #define yheight (_tofT->_parameters.yheight)
  _tofT->_parameters.tmin = 0;
  #define tmin (_tofT->_parameters.tmin)
  _tofT->_parameters.tmax = 10000;
  #define tmax (_tofT->_parameters.tmax)
  _tofT->_parameters.dt = 1.0;
  #define dt (_tofT->_parameters.dt)
  _tofT->_parameters.restore_neutron = 1;
  #define restore_neutron (_tofT->_parameters.restore_neutron)

  #define TOF_N (_tofT->_parameters.TOF_N)
  #define TOF_p (_tofT->_parameters.TOF_p)
  #define TOF_p2 (_tofT->_parameters.TOF_p2)
  #define t_min (_tofT->_parameters.t_min)
  #define t_max (_tofT->_parameters.t_max)
  #define delta_t (_tofT->_parameters.delta_t)

  #undef nt
  #undef filename
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef xwidth
  #undef yheight
  #undef tmin
  #undef tmax
  #undef dt
  #undef restore_neutron
  #undef TOF_N
  #undef TOF_p
  #undef TOF_p2
  #undef t_min
  #undef t_max
  #undef delta_t
  /* component tofT=TOF_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _PortOrig->_rotation_absolute, _tofT->_rotation_absolute);
    rot_transpose(_XmonC->_rotation_absolute, tr1);
    rot_mul(_tofT->_rotation_absolute, tr1, _tofT->_rotation_relative);
    _tofT->_rotation_is_identity =  rot_test_identity(_tofT->_rotation_relative);
    tc1 = coords_set(
      0, 0, 0.3);
    rot_transpose(_PortOrig->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _tofT->_position_absolute = coords_add(_PortOrig->_position_absolute, tc2);
    tc1 = coords_sub(_XmonC->_position_absolute, _tofT->_position_absolute);
    _tofT->_position_relative = rot_apply(_tofT->_rotation_absolute, tc1);
  } /* tofT=TOF_monitor() AT ROTATED */
  DEBUG_COMPONENT("tofT", _tofT->_position_absolute, _tofT->_rotation_absolute);
  instrument->_position_absolute[6] = _tofT->_position_absolute;
  instrument->_position_relative[6] = _tofT->_position_relative;
  instrument->counter_N[6]  = instrument->counter_P[6] = instrument->counter_P2[6] = 0;
  instrument->counter_AbsorbProp[6]= 0;
  return(0);
} /* _tofT_setpos */

/* component tofC=TOF_monitor() SETTING, POSITION/ROTATION */
int _tofC_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_tofC_setpos] component tofC=TOF_monitor() SETTING [/usr/share/mcstas/3.0-dev/monitors/TOF_monitor.comp:63]");
  stracpy(_tofC->_name, "tofC", 16384);
  stracpy(_tofC->_type, "TOF_monitor", 16384);
  _tofC->_index=7;
  _tofC->_parameters.nt = 100;
  #define nt (_tofC->_parameters.nt)
  if("tofC.dat" && strlen("tofC.dat"))
    stracpy(_tofC->_parameters.filename, "tofC.dat" ? "tofC.dat" : "", 16384);
  else 
  _tofC->_parameters.filename[0]='\0';
  #define filename (_tofC->_parameters.filename)
  _tofC->_parameters.xmin = - mon_width / 2 + instrument->_parameters._mon_shift;
  #define xmin (_tofC->_parameters.xmin)
  _tofC->_parameters.xmax = mon_width / 2 + instrument->_parameters._mon_shift;
  #define xmax (_tofC->_parameters.xmax)
  _tofC->_parameters.ymin = -0.05;
  #define ymin (_tofC->_parameters.ymin)
  _tofC->_parameters.ymax = 0.05;
  #define ymax (_tofC->_parameters.ymax)
  _tofC->_parameters.xwidth = 0;
  #define xwidth (_tofC->_parameters.xwidth)
  _tofC->_parameters.yheight = 0.04;
  #define yheight (_tofC->_parameters.yheight)
  _tofC->_parameters.tmin = 0;
  #define tmin (_tofC->_parameters.tmin)
  _tofC->_parameters.tmax = 10000;
  #define tmax (_tofC->_parameters.tmax)
  _tofC->_parameters.dt = 1.0;
  #define dt (_tofC->_parameters.dt)
  _tofC->_parameters.restore_neutron = 1;
  #define restore_neutron (_tofC->_parameters.restore_neutron)

  #define TOF_N (_tofC->_parameters.TOF_N)
  #define TOF_p (_tofC->_parameters.TOF_p)
  #define TOF_p2 (_tofC->_parameters.TOF_p2)
  #define t_min (_tofC->_parameters.t_min)
  #define t_max (_tofC->_parameters.t_max)
  #define delta_t (_tofC->_parameters.delta_t)

  #undef nt
  #undef filename
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef xwidth
  #undef yheight
  #undef tmin
  #undef tmax
  #undef dt
  #undef restore_neutron
  #undef TOF_N
  #undef TOF_p
  #undef TOF_p2
  #undef t_min
  #undef t_max
  #undef delta_t
  /* component tofC=TOF_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _PortOrig->_rotation_absolute, _tofC->_rotation_absolute);
    rot_transpose(_tofT->_rotation_absolute, tr1);
    rot_mul(_tofC->_rotation_absolute, tr1, _tofC->_rotation_relative);
    _tofC->_rotation_is_identity =  rot_test_identity(_tofC->_rotation_relative);
    tc1 = coords_set(
      0, 0, 0.3);
    rot_transpose(_PortOrig->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _tofC->_position_absolute = coords_add(_PortOrig->_position_absolute, tc2);
    tc1 = coords_sub(_tofT->_position_absolute, _tofC->_position_absolute);
    _tofC->_position_relative = rot_apply(_tofC->_rotation_absolute, tc1);
  } /* tofC=TOF_monitor() AT ROTATED */
  DEBUG_COMPONENT("tofC", _tofC->_position_absolute, _tofC->_rotation_absolute);
  instrument->_position_absolute[7] = _tofC->_position_absolute;
  instrument->_position_relative[7] = _tofC->_position_relative;
  instrument->counter_N[7]  = instrument->counter_P[7] = instrument->counter_P2[7] = 0;
  instrument->counter_AbsorbProp[7]= 0;
  return(0);
} /* _tofC_setpos */

/* component Slit1=Slit() SETTING, POSITION/ROTATION */
int _Slit1_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Slit1_setpos] component Slit1=Slit() SETTING [/usr/share/mcstas/3.0-dev/optics/Slit.comp:47]");
  stracpy(_Slit1->_name, "Slit1", 16384);
  stracpy(_Slit1->_type, "Slit", 16384);
  _Slit1->_index=8;
  _Slit1->_parameters.xmin = -0.01;
  #define xmin (_Slit1->_parameters.xmin)
  _Slit1->_parameters.xmax = 0.01;
  #define xmax (_Slit1->_parameters.xmax)
  _Slit1->_parameters.ymin = -0.01;
  #define ymin (_Slit1->_parameters.ymin)
  _Slit1->_parameters.ymax = 0.01;
  #define ymax (_Slit1->_parameters.ymax)
  _Slit1->_parameters.radius = 0;
  #define radius (_Slit1->_parameters.radius)
  _Slit1->_parameters.xwidth = 0.01;
  #define xwidth (_Slit1->_parameters.xwidth)
  _Slit1->_parameters.yheight = 0.01;
  #define yheight (_Slit1->_parameters.yheight)

  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef radius
  #undef xwidth
  #undef yheight
  /* component Slit1=Slit() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _PortOrig->_rotation_absolute, _Slit1->_rotation_absolute);
    rot_transpose(_tofC->_rotation_absolute, tr1);
    rot_mul(_Slit1->_rotation_absolute, tr1, _Slit1->_rotation_relative);
    _Slit1->_rotation_is_identity =  rot_test_identity(_Slit1->_rotation_relative);
    tc1 = coords_set(
      0, 0, 100);
    rot_transpose(_PortOrig->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Slit1->_position_absolute = coords_add(_PortOrig->_position_absolute, tc2);
    tc1 = coords_sub(_tofC->_position_absolute, _Slit1->_position_absolute);
    _Slit1->_position_relative = rot_apply(_Slit1->_rotation_absolute, tc1);
  } /* Slit1=Slit() AT ROTATED */
  DEBUG_COMPONENT("Slit1", _Slit1->_position_absolute, _Slit1->_rotation_absolute);
  instrument->_position_absolute[8] = _Slit1->_position_absolute;
  instrument->_position_relative[8] = _Slit1->_position_relative;
  instrument->counter_N[8]  = instrument->counter_P[8] = instrument->counter_P2[8] = 0;
  instrument->counter_AbsorbProp[8]= 0;
  return(0);
} /* _Slit1_setpos */

/* component LamT=L_monitor() SETTING, POSITION/ROTATION */
int _LamT_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_LamT_setpos] component LamT=L_monitor() SETTING [/usr/share/mcstas/3.0-dev/monitors/L_monitor.comp:65]");
  stracpy(_LamT->_name, "LamT", 16384);
  stracpy(_LamT->_type, "L_monitor", 16384);
  _LamT->_index=9;
  _LamT->_parameters.nL = 160;
  #define nL (_LamT->_parameters.nL)
  if("LamT.dat" && strlen("LamT.dat"))
    stracpy(_LamT->_parameters.filename, "LamT.dat" ? "LamT.dat" : "", 16384);
  else 
  _LamT->_parameters.filename[0]='\0';
  #define filename (_LamT->_parameters.filename)
  _LamT->_parameters.xmin = -0.05;
  #define xmin (_LamT->_parameters.xmin)
  _LamT->_parameters.xmax = 0.05;
  #define xmax (_LamT->_parameters.xmax)
  _LamT->_parameters.ymin = -0.05;
  #define ymin (_LamT->_parameters.ymin)
  _LamT->_parameters.ymax = 0.05;
  #define ymax (_LamT->_parameters.ymax)
  _LamT->_parameters.xwidth = 0.02;
  #define xwidth (_LamT->_parameters.xwidth)
  _LamT->_parameters.yheight = 0.02;
  #define yheight (_LamT->_parameters.yheight)
  _LamT->_parameters.Lmin = instrument->_parameters._lmin;
  #define Lmin (_LamT->_parameters.Lmin)
  _LamT->_parameters.Lmax = instrument->_parameters._lmax;
  #define Lmax (_LamT->_parameters.Lmax)
  _LamT->_parameters.restore_neutron = 1;
  #define restore_neutron (_LamT->_parameters.restore_neutron)

  #define L_N (_LamT->_parameters.L_N)
  #define L_p (_LamT->_parameters.L_p)
  #define L_p2 (_LamT->_parameters.L_p2)

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
  /* component LamT=L_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Slit1->_rotation_absolute, _LamT->_rotation_absolute);
    rot_transpose(_Slit1->_rotation_absolute, tr1);
    rot_mul(_LamT->_rotation_absolute, tr1, _LamT->_rotation_relative);
    _LamT->_rotation_is_identity =  rot_test_identity(_LamT->_rotation_relative);
    tc1 = coords_set(
      0, 0, 0.0001);
    rot_transpose(_Slit1->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _LamT->_position_absolute = coords_add(_Slit1->_position_absolute, tc2);
    tc1 = coords_sub(_Slit1->_position_absolute, _LamT->_position_absolute);
    _LamT->_position_relative = rot_apply(_LamT->_rotation_absolute, tc1);
  } /* LamT=L_monitor() AT ROTATED */
  DEBUG_COMPONENT("LamT", _LamT->_position_absolute, _LamT->_rotation_absolute);
  instrument->_position_absolute[9] = _LamT->_position_absolute;
  instrument->_position_relative[9] = _LamT->_position_relative;
  instrument->counter_N[9]  = instrument->counter_P[9] = instrument->counter_P2[9] = 0;
  instrument->counter_AbsorbProp[9]= 0;
  return(0);
} /* _LamT_setpos */

/* component LamC=L_monitor() SETTING, POSITION/ROTATION */
int _LamC_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_LamC_setpos] component LamC=L_monitor() SETTING [/usr/share/mcstas/3.0-dev/monitors/L_monitor.comp:65]");
  stracpy(_LamC->_name, "LamC", 16384);
  stracpy(_LamC->_type, "L_monitor", 16384);
  _LamC->_index=10;
  _LamC->_parameters.nL = 160;
  #define nL (_LamC->_parameters.nL)
  if("LamC.dat" && strlen("LamC.dat"))
    stracpy(_LamC->_parameters.filename, "LamC.dat" ? "LamC.dat" : "", 16384);
  else 
  _LamC->_parameters.filename[0]='\0';
  #define filename (_LamC->_parameters.filename)
  _LamC->_parameters.xmin = -0.05;
  #define xmin (_LamC->_parameters.xmin)
  _LamC->_parameters.xmax = 0.05;
  #define xmax (_LamC->_parameters.xmax)
  _LamC->_parameters.ymin = -0.05;
  #define ymin (_LamC->_parameters.ymin)
  _LamC->_parameters.ymax = 0.05;
  #define ymax (_LamC->_parameters.ymax)
  _LamC->_parameters.xwidth = 0.02;
  #define xwidth (_LamC->_parameters.xwidth)
  _LamC->_parameters.yheight = 0.02;
  #define yheight (_LamC->_parameters.yheight)
  _LamC->_parameters.Lmin = instrument->_parameters._lmin;
  #define Lmin (_LamC->_parameters.Lmin)
  _LamC->_parameters.Lmax = instrument->_parameters._lmax;
  #define Lmax (_LamC->_parameters.Lmax)
  _LamC->_parameters.restore_neutron = 1;
  #define restore_neutron (_LamC->_parameters.restore_neutron)

  #define L_N (_LamC->_parameters.L_N)
  #define L_p (_LamC->_parameters.L_p)
  #define L_p2 (_LamC->_parameters.L_p2)

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
  /* component LamC=L_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Slit1->_rotation_absolute, _LamC->_rotation_absolute);
    rot_transpose(_LamT->_rotation_absolute, tr1);
    rot_mul(_LamC->_rotation_absolute, tr1, _LamC->_rotation_relative);
    _LamC->_rotation_is_identity =  rot_test_identity(_LamC->_rotation_relative);
    tc1 = coords_set(
      0, 0, 0.0001);
    rot_transpose(_Slit1->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _LamC->_position_absolute = coords_add(_Slit1->_position_absolute, tc2);
    tc1 = coords_sub(_LamT->_position_absolute, _LamC->_position_absolute);
    _LamC->_position_relative = rot_apply(_LamC->_rotation_absolute, tc1);
  } /* LamC=L_monitor() AT ROTATED */
  DEBUG_COMPONENT("LamC", _LamC->_position_absolute, _LamC->_rotation_absolute);
  instrument->_position_absolute[10] = _LamC->_position_absolute;
  instrument->_position_relative[10] = _LamC->_position_relative;
  instrument->counter_N[10]  = instrument->counter_P[10] = instrument->counter_P2[10] = 0;
  instrument->counter_AbsorbProp[10]= 0;
  return(0);
} /* _LamC_setpos */

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

_class_ESS_moderator *class_ESS_moderator_init(_class_ESS_moderator *_comp
) {
  #define isleft (_comp->_parameters.isleft)
  #define Lmin (_comp->_parameters.Lmin)
  #define Lmax (_comp->_parameters.Lmax)
  #define cold_frac (_comp->_parameters.cold_frac)
  #define dist (_comp->_parameters.dist)
  #define focus_xw (_comp->_parameters.focus_xw)
  #define focus_yh (_comp->_parameters.focus_yh)
  #define target_index (_comp->_parameters.target_index)
  #define tmax_multiplier (_comp->_parameters.tmax_multiplier)
  #define yheight_c (_comp->_parameters.yheight_c)
  #define yheight_t (_comp->_parameters.yheight_t)
  #define n_pulses (_comp->_parameters.n_pulses)
  #define acc_power (_comp->_parameters.acc_power)
  #define beamport_angle (_comp->_parameters.beamport_angle)
  #define sourcedef (_comp->_parameters.sourcedef)
  #define xwidth_c (_comp->_parameters.xwidth_c)
  #define xwidth_t (_comp->_parameters.xwidth_t)
  #define extraction_opening (_comp->_parameters.extraction_opening)
  #define l_range (_comp->_parameters.l_range)
  #define w_mult (_comp->_parameters.w_mult)
  #define w_geom (_comp->_parameters.w_geom)
  #define w_geom_c (_comp->_parameters.w_geom_c)
  #define w_geom_t (_comp->_parameters.w_geom_t)
  #define w_stat (_comp->_parameters.w_stat)
  #define cold (_comp->_parameters.cold)
  #define wasleft (_comp->_parameters.wasleft)
  #define cosine (_comp->_parameters.cosine)
  #define tx (_comp->_parameters.tx)
  #define ty (_comp->_parameters.ty)
  #define tz (_comp->_parameters.tz)
  #define cx (_comp->_parameters.cx)
  #define cy (_comp->_parameters.cy)
  #define cz (_comp->_parameters.cz)
  #define cyl_radius (_comp->_parameters.cyl_radius)
  #define t1x (_comp->_parameters.t1x)
  #define t1y (_comp->_parameters.t1y)
  #define t1z (_comp->_parameters.t1z)
  #define t2x (_comp->_parameters.t2x)
  #define t2y (_comp->_parameters.t2y)
  #define t2z (_comp->_parameters.t2z)
  #define modextras (_comp->_parameters.modextras)
  #define cold_bril (_comp->_parameters.cold_bril)
  #define thermal_bril (_comp->_parameters.thermal_bril)
  #define tfocus_width (_comp->_parameters.tfocus_width)
  #define tfocus_time (_comp->_parameters.tfocus_time)
  #define dt (_comp->_parameters.dt)
  #define r_empty (_comp->_parameters.r_empty)
  #define lambda (_comp->_parameters.lambda)
  #define xprime (_comp->_parameters.xprime)
  #define yprime (_comp->_parameters.yprime)
  #define zprime (_comp->_parameters.zprime)
  #define sgnport (_comp->_parameters.sgnport)
  #define abs_beamport_angle (_comp->_parameters.abs_beamport_angle)
  #define cos_beamport_angle (_comp->_parameters.cos_beamport_angle)
  #define edge_thermal (_comp->_parameters.edge_thermal)
  #define cos_thermal (_comp->_parameters.cos_thermal)
  #define cos_cold (_comp->_parameters.cos_cold)
  #define sin_beamport_angle (_comp->_parameters.sin_beamport_angle)

  if (Lmin>=Lmax) {
    printf("ESS_moderator: %s: Unmeaningful definition of wavelength range!\n ERROR - Exiting\n",
           NAME_CURRENT_COMP);
      exit(0);
  }

  sgnport=(beamport_angle>0 ? 1:-1);
  if (sgnport<0) {
	printf("WARNING, %s: beamport angle < 0, experimental option\n", NAME_CURRENT_COMP);
  }

/* ESS 2015 - variables needed to correct for the emission surface angle */
  abs_beamport_angle=fabs(beamport_angle);
  cos_beamport_angle=cos(beamport_angle*DEG2RAD);
  sin_beamport_angle=sin(beamport_angle*DEG2RAD);
  /* correction for projection along the beam / projection on the z=0 plane */
  cos_thermal=cos_beamport_angle;
  cos_cold=cos((abs_beamport_angle-24.24)*DEG2RAD)/cos(24.24*DEG2RAD);
  printf("input to cosine %g, output %g\n",abs_beamport_angle, cos_cold);
  /* cross-section between the z=const and tilted surfaces [m] */
  edge_thermal=-0.0716;

  r_empty=2;

  /* Weight multipliers */
  w_mult=acc_power/5;
  w_geom=1;
  w_geom_c=1;
  w_geom_t=1;
  w_stat=1.0/mcget_ncount();

  tfocus_width=0;  tfocus_time=0;  dt=0;

  n_pulses=(double)floor(n_pulses);
  if (n_pulses == 0) n_pulses=1;

  /* Calculation of location of focusing rectangle */
  if (target_index && !dist) {
    Coords ToTarget;
    ToTarget = coords_sub(POS_A_COMP_INDEX(INDEX_CURRENT_COMP+target_index),POS_A_CURRENT_COMP);
    ToTarget = rot_apply(ROT_A_CURRENT_COMP, ToTarget);
    coords_get(ToTarget, &tx, &ty, &tz);
	printf("%s: Focusing on a window centered at [%g, %g, %g]\n", NAME_CURRENT_COMP,tx,ty,tz);
    dist=sqrt(tx*tx+ty*ty+tz*tz);
  } else if (target_index && !dist) {
    printf("ESS_moderator: %s: Please choose to set either the dist parameter or specify a target_index.\nExit\n", NAME_CURRENT_COMP);
    exit(-1);
  } else {
    tx=0, ty=0, tz=dist;
  }

  if (focus_xw < 0 || focus_yh < 0) {
    printf("ESS_moderator: %s: Please specify both focus_xw and focus_yh as positive numbers.\nExit\n", NAME_CURRENT_COMP);
    exit(-1);
  }

  if (Lmin<=0 || Lmax <=0 || dist == 0)
  {
    printf("ESS_moderator: %s: Check parameters (lead to Math Error).\n Avoid 0 value for {Lmin Lmax dist} \n", NAME_CURRENT_COMP);
    exit(-1);
  }

  /* Calculate cold cylinder radius from extraction_opening */
  cyl_radius = 360*xwidth_c/(2*PI*extraction_opening);

  /* Is abs_beamport_angle set at a meaningful angle - otherwise recalc */
  if ( abs_beamport_angle<0 || abs_beamport_angle>extraction_opening){
    beamport_angle = extraction_opening/2.0;
	abs_beamport_angle=beamport_angle;
    printf("%s: WARNING: Resetting your beamport_angle to %g (%g/2) - central beamslot\n",NAME_CURRENT_COMP,beamport_angle,extraction_opening);
  }

  /* Calculate positioning of thermal slabs and cold cylinder - NOTE: overwritten in case of ESS 2014 */
  t1z = cyl_radius*cos(-DEG2RAD*abs_beamport_angle);
  t1x = cyl_radius*sin(-DEG2RAD*abs_beamport_angle);
  t1y = 0;
  /* Wing 2 (right) is at extraction_width-beamport_angle */
  t2z = cyl_radius*cos(DEG2RAD*(extraction_opening-abs_beamport_angle));
  t2x = cyl_radius*sin(DEG2RAD*(extraction_opening-abs_beamport_angle));
  t2y = 0;
  /* We want unit vectors... */
  NORM(t1x,t1y,t1z);
  NORM(t2x,t2y,t2z);
  cx = 0; cy=0; cz=0;

  /* Geometry parameters and brilliances specified via the sourcedef string */
  /* ESS-TDR definition: */
  if (strcasestr(sourcedef,"TDR")) {
    w_mult *= ESS_SOURCE_FREQUENCY;               /* Correct for frequency */
    printf("Using ESS TDR brilliance\n");
    cold_bril=ESS_2012_Lieutenant_cold;
    thermal_bril=ESS_Mezei_thermal;
 /* ESS-2001 definition: */
  } else if (strcasestr(sourcedef,"2001")) {
    w_mult *= ESS_SOURCE_FREQUENCY;               /* Correct for frequency */
    printf("Using ESS 2001 brilliance\n");
    cold_bril=ESS_Mezei_cold;
    thermal_bril=ESS_Mezei_thermal;
 /* ESS-2013 post-TDR updated definition, before ESS "pancake": */
  } else if (strcasestr(sourcedef,"2013")) {
    w_mult *= acc_power;                             /* Is already in per-MW units */
    printf("Using ESS 2013 brilliance\n");
    modextras.height_c=yheight_c;
    modextras.height_t=yheight_t;
    modextras.tmultiplier=tmax_multiplier;
    cold_bril=ESS_2013_Schoenfeldt_cold;
    thermal_bril=ESS_2013_Schoenfeldt_thermal;
  /* ESS-2014 definition: */
  } else if (strcasestr(sourcedef,"2014")) {
    if (xwidth_c!=0.23 || xwidth_t!=0.12 || extraction_opening!=120) {
      fprintf(stderr,"FATAL: sourcedef %s only allows xwidth_c=0.23, xwidth_t=0.12 and extraction_opening=120 !\n",sourcedef);
      exit(-1);
    }
    /* Specify brilliance fct.'s */
    cold_bril=ESS_2014_Schoenfeldt_cold;
    thermal_bril=ESS_2014_Schoenfeldt_thermal;
    /* Cold moderator geometry */
    /* In this mode of operation, the emission is from the x-y plane*/
    cx = 0; cy=0; cz=-0.12; cyl_radius=0.12;
    /* 120 degree beam extraction opening, defining thermal moderator directions */
    t1x=sin(DEG2RAD*(120/2-abs_beamport_angle));t1z=cos(DEG2RAD*120/2);
    t2x=-sin(DEG2RAD*120/2);t2z=cos(DEG2RAD*120/2);
    t1y=0; t2y=0;
    /* ESS-2015 definition: */
   } else if (strcasestr(sourcedef,"2015")) {
    if (xwidth_c!=0.1 || xwidth_t!=0.18 || extraction_opening!=120) {
      fprintf(stderr,"FATAL: sourcedef %s only allows xwidth_c=0.1, xwidth_t=0.18 and extraction_opening=120 !\n",sourcedef);
      exit(-1);
    }
    /* Specify brilliance fct.'s */
    cold_bril=ESS_2015_Schoenfeldt_cold;
    thermal_bril=ESS_2015_Schoenfeldt_thermal;
    /* Moderator geometry */
    cx = 0; cy=0; cz=0; cyl_radius=0;
    /* 120 degree beam extraction opening */
    t1x=sin(DEG2RAD*(120/2+abs_beamport_angle));t1z=cos(DEG2RAD*(120/2+abs_beamport_angle));
    t2x=-sin(DEG2RAD*(120/2-abs_beamport_angle));t2z=cos(DEG2RAD*(120/2-abs_beamport_angle));
    t1y=0; t2y=0;
   } else {
    fprintf(stderr,"FATAL: sourcedef %s is not defined!\n",sourcedef);
    exit(-1);
  }
  modextras.height_c=yheight_c;
  modextras.Width_c=xwidth_c;
  modextras.Width_t=xwidth_t;
  modextras.height_t=yheight_t;
  modextras.tmultiplier=tmax_multiplier;
  modextras.extractionangle=extraction_opening;

  if ((strcasestr(sourcedef,"2014")) && abs_beamport_angle != modextras.extractionangle/2.0) {
    printf("%s: WARNING: beamport_angle is not at central slot. With sourcedef %s this means application of a partially analytical brightness model for B(\\theta) \n",NAME_CURRENT_COMP,sourcedef);
  }
  if ((strcasestr(sourcedef,"2015")) && (abs_beamport_angle!=5.0 && abs_beamport_angle!=15.0 && abs_beamport_angle!=25.0 && abs_beamport_angle!=35.0 && abs_beamport_angle!=45.0 && abs_beamport_angle!=55.0)) {
    printf("%s: FATAL: beamport_angle is not at a defined value. Allowed are: (5, 15, ...) and < 60  \n",NAME_CURRENT_COMP);
    exit(-1);
  }
  modextras.beamportangle=abs_beamport_angle;

  l_range = Lmax-Lmin;
  w_geom_c  = xwidth_c*yheight_c*1.0e4;     /* source area correction */
  w_geom_t  = xwidth_t*yheight_t*1.0e4;
  w_mult *= l_range;            /* wavelength range correction */

  #undef isleft
  #undef Lmin
  #undef Lmax
  #undef cold_frac
  #undef dist
  #undef focus_xw
  #undef focus_yh
  #undef target_index
  #undef tmax_multiplier
  #undef yheight_c
  #undef yheight_t
  #undef n_pulses
  #undef acc_power
  #undef beamport_angle
  #undef sourcedef
  #undef xwidth_c
  #undef xwidth_t
  #undef extraction_opening
  #undef l_range
  #undef w_mult
  #undef w_geom
  #undef w_geom_c
  #undef w_geom_t
  #undef w_stat
  #undef cold
  #undef wasleft
  #undef cosine
  #undef tx
  #undef ty
  #undef tz
  #undef cx
  #undef cy
  #undef cz
  #undef cyl_radius
  #undef t1x
  #undef t1y
  #undef t1z
  #undef t2x
  #undef t2y
  #undef t2z
  #undef modextras
  #undef cold_bril
  #undef thermal_bril
  #undef tfocus_width
  #undef tfocus_time
  #undef dt
  #undef r_empty
  #undef lambda
  #undef xprime
  #undef yprime
  #undef zprime
  #undef sgnport
  #undef abs_beamport_angle
  #undef cos_beamport_angle
  #undef edge_thermal
  #undef cos_thermal
  #undef cos_cold
  #undef sin_beamport_angle
  return(_comp);
} /* class_ESS_moderator_init */

_class_PSDlin_monitor *class_PSDlin_monitor_init(_class_PSDlin_monitor *_comp
) {
  #define nx (_comp->_parameters.nx)
  #define filename (_comp->_parameters.filename)
  #define xmin (_comp->_parameters.xmin)
  #define xmax (_comp->_parameters.xmax)
  #define ymin (_comp->_parameters.ymin)
  #define ymax (_comp->_parameters.ymax)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define restore_neutron (_comp->_parameters.restore_neutron)
  #define PSDlin_N (_comp->_parameters.PSDlin_N)
  #define PSDlin_p (_comp->_parameters.PSDlin_p)
  #define PSDlin_p2 (_comp->_parameters.PSDlin_p2)
  int i;

  if (xwidth  > 0) { xmax = xwidth/2;  xmin = -xmax; }
  if (yheight > 0) { ymax = yheight/2; ymin = -ymax; }

  if ((xmin >= xmax) || (ymin >= ymax)) {
    printf("PSDlin_monitor: %s: Null detection area !\n"
           "ERROR           (xwidth,yheight,xmin,xmax,ymin,ymax). Exiting",
           NAME_CURRENT_COMP);
    exit(0);
  }

  PSDlin_N = create_darr1d(nx);
  PSDlin_p = create_darr1d(nx);
  PSDlin_p2 = create_darr1d(nx);

  for (i=0; i<nx; i++)
  {
    PSDlin_N[i] = 0;
    PSDlin_p[i] = 0;
    PSDlin_p2[i] = 0;
  }
  #undef nx
  #undef filename
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef xwidth
  #undef yheight
  #undef restore_neutron
  #undef PSDlin_N
  #undef PSDlin_p
  #undef PSDlin_p2
  return(_comp);
} /* class_PSDlin_monitor_init */

_class_TOF_monitor *class_TOF_monitor_init(_class_TOF_monitor *_comp
) {
  #define nt (_comp->_parameters.nt)
  #define filename (_comp->_parameters.filename)
  #define xmin (_comp->_parameters.xmin)
  #define xmax (_comp->_parameters.xmax)
  #define ymin (_comp->_parameters.ymin)
  #define ymax (_comp->_parameters.ymax)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define tmin (_comp->_parameters.tmin)
  #define tmax (_comp->_parameters.tmax)
  #define dt (_comp->_parameters.dt)
  #define restore_neutron (_comp->_parameters.restore_neutron)
  #define TOF_N (_comp->_parameters.TOF_N)
  #define TOF_p (_comp->_parameters.TOF_p)
  #define TOF_p2 (_comp->_parameters.TOF_p2)
  #define t_min (_comp->_parameters.t_min)
  #define t_max (_comp->_parameters.t_max)
  #define delta_t (_comp->_parameters.delta_t)
  int i;

  if (xwidth  > 0) { xmax = xwidth/2;  xmin = -xmax; }
  if (yheight > 0) { ymax = yheight/2; ymin = -ymax; }

  if ((xmin >= xmax) || (ymin >= ymax)) {
          printf("TOF_monitor: %s: Null detection area !\n"
                 "ERROR        (xwidth,yheight,xmin,xmax,ymin,ymax). Exiting",
         NAME_CURRENT_COMP);
    exit(0);
  }

  TOF_N = create_darr1d(nt);
  TOF_p = create_darr1d(nt);
  TOF_p2 = create_darr1d(nt);

  for (i=0; i<nt; i++)
  {
    TOF_N[i] = 0;
    TOF_p[i] = 0;
    TOF_p2[i] = 0;
  }
  if (tmax!=0)
  {
    t_max=tmax;
    t_min=tmin;
    delta_t=(t_max-t_min)/nt;
  }
  else
  {
    delta_t=dt;
    t_min=0;
    t_max=nt*dt+tmin;
  }
  #undef nt
  #undef filename
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef xwidth
  #undef yheight
  #undef tmin
  #undef tmax
  #undef dt
  #undef restore_neutron
  #undef TOF_N
  #undef TOF_p
  #undef TOF_p2
  #undef t_min
  #undef t_max
  #undef delta_t
  return(_comp);
} /* class_TOF_monitor_init */

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



int init(void) { /* called by mccode_main for ESS_2015_test:INITIALISE */
  DEBUG_INSTR();

  /* code_main/parseoptions/readparams sets instrument parameters value */
  stracpy(instrument->_name, "ESS_2015_test", 256);

  /* Instrument 'ESS_2015_test' INITIALISE */
  SIG_MESSAGE("[ESS_2015_test] INITIALISE [ESS_2015_test.instr:52]");
  #define frac (instrument->_parameters._frac)
  #define lmin (instrument->_parameters._lmin)
  #define lmax (instrument->_parameters._lmax)
  #define ANGLE (instrument->_parameters._ANGLE)
  #define mon_shift (instrument->_parameters._mon_shift)
{
  power=5;
  Yheight=0.03;
  np=1;
  mon_width=0.25;
  sgn= (ANGLE>0 ? 1:-1);
  cosa=cos(ANGLE*DEG2RAD);
  sina=sin(ANGLE*DEG2RAD);
  /* centre of the beamline coordinates w.r.t. moderator cooordinates  */
  ctrX=-sgn*0.0716;ctrZ=0.0827;
  /* rotated by beamport angle */
  BX0=ctrX*cosa + ctrZ*sina;
  BZ0=-ctrX*sina + ctrZ*cosa;
  printf("Beamline centered at [cm]: x= %g,  z=%g \n",ctrX*100,ctrZ*100);
  printf("after rotation by beamport angle [cm]: x= %g,  z=%g \n",BX0*100,BZ0*100);
  // printf("maximum extent of cold area [cm]: x= %g\n",(0.1425*cosa + 0.05077*sina - BX0)*100);

  sprintf(srcdef,"2015");
}
  #undef frac
  #undef lmin
  #undef lmax
  #undef ANGLE
  #undef mon_shift
  _Origin_setpos(); /* type Progress_bar */
  _Source_setpos(); /* type ESS_moderator */
  _PortOrig_setpos(); /* type Arm */
  _XmonT_setpos(); /* type PSDlin_monitor */
  _XmonC_setpos(); /* type PSDlin_monitor */
  _tofT_setpos(); /* type TOF_monitor */
  _tofC_setpos(); /* type TOF_monitor */
  _Slit1_setpos(); /* type Slit */
  _LamT_setpos(); /* type L_monitor */
  _LamC_setpos(); /* type L_monitor */

  /* call iteratively all components INITIALISE */
  class_Progress_bar_init(_Origin);

  class_ESS_moderator_init(_Source);


  class_PSDlin_monitor_init(_XmonT);

  class_PSDlin_monitor_init(_XmonC);

  class_TOF_monitor_init(_tofT);

  class_TOF_monitor_init(_tofC);

  class_Slit_init(_Slit1);

  class_L_monitor_init(_LamT);

  class_L_monitor_init(_LamC);

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
_class_ESS_moderator *class_ESS_moderator_trace(_class_ESS_moderator *_comp
  , _class_particle *_particle) {
  ABSORBED=SCATTERED=RESTORE=0;

  #define isleft (_comp->_parameters.isleft)
  #define Lmin (_comp->_parameters.Lmin)
  #define Lmax (_comp->_parameters.Lmax)
  #define cold_frac (_comp->_parameters.cold_frac)
  #define dist (_comp->_parameters.dist)
  #define focus_xw (_comp->_parameters.focus_xw)
  #define focus_yh (_comp->_parameters.focus_yh)
  #define target_index (_comp->_parameters.target_index)
  #define tmax_multiplier (_comp->_parameters.tmax_multiplier)
  #define yheight_c (_comp->_parameters.yheight_c)
  #define yheight_t (_comp->_parameters.yheight_t)
  #define n_pulses (_comp->_parameters.n_pulses)
  #define acc_power (_comp->_parameters.acc_power)
  #define beamport_angle (_comp->_parameters.beamport_angle)
  #define sourcedef (_comp->_parameters.sourcedef)
  #define xwidth_c (_comp->_parameters.xwidth_c)
  #define xwidth_t (_comp->_parameters.xwidth_t)
  #define extraction_opening (_comp->_parameters.extraction_opening)
  #define l_range (_comp->_parameters.l_range)
  #define w_mult (_comp->_parameters.w_mult)
  #define w_geom (_comp->_parameters.w_geom)
  #define w_geom_c (_comp->_parameters.w_geom_c)
  #define w_geom_t (_comp->_parameters.w_geom_t)
  #define w_stat (_comp->_parameters.w_stat)
  #define cold (_comp->_parameters.cold)
  #define wasleft (_comp->_parameters.wasleft)
  #define cosine (_comp->_parameters.cosine)
  #define tx (_comp->_parameters.tx)
  #define ty (_comp->_parameters.ty)
  #define tz (_comp->_parameters.tz)
  #define cx (_comp->_parameters.cx)
  #define cy (_comp->_parameters.cy)
  #define cz (_comp->_parameters.cz)
  #define cyl_radius (_comp->_parameters.cyl_radius)
  #define t1x (_comp->_parameters.t1x)
  #define t1y (_comp->_parameters.t1y)
  #define t1z (_comp->_parameters.t1z)
  #define t2x (_comp->_parameters.t2x)
  #define t2y (_comp->_parameters.t2y)
  #define t2z (_comp->_parameters.t2z)
  #define modextras (_comp->_parameters.modextras)
  #define cold_bril (_comp->_parameters.cold_bril)
  #define thermal_bril (_comp->_parameters.thermal_bril)
  #define tfocus_width (_comp->_parameters.tfocus_width)
  #define tfocus_time (_comp->_parameters.tfocus_time)
  #define dt (_comp->_parameters.dt)
  #define r_empty (_comp->_parameters.r_empty)
  #define lambda (_comp->_parameters.lambda)
  #define xprime (_comp->_parameters.xprime)
  #define yprime (_comp->_parameters.yprime)
  #define zprime (_comp->_parameters.zprime)
  #define sgnport (_comp->_parameters.sgnport)
  #define abs_beamport_angle (_comp->_parameters.abs_beamport_angle)
  #define cos_beamport_angle (_comp->_parameters.cos_beamport_angle)
  #define edge_thermal (_comp->_parameters.edge_thermal)
  #define cos_thermal (_comp->_parameters.cos_thermal)
  #define cos_cold (_comp->_parameters.cos_cold)
  #define sin_beamport_angle (_comp->_parameters.sin_beamport_angle)
  p=1;
  double v,E,k,r,xf,yf,dx,dy,dz,w_focus,tail_flag,cor;

  /* Bispectral source - choice of spectrum and initial position */
  cold = ( rand01() < cold_frac );

  /* Use standard McStas routine or -not for specifying emission directionality */
  cosine = 0;

  /* Emission geometry adapted from ESS MCNPX model, mid 2012 */
  if (strcasestr(sourcedef,"2014")) {
    if (cold) { //case: cold moderator
      x = 0.5*randpm1()*xwidth_c;
      y = 0.5*randpm1()*yheight_c;
      z = 0;
      cosine = 0;
      w_geom = w_geom_c;
    }  else  { //case: thermal moderator
      y = 0.5*randpm1()*yheight_t;
      if (rand01()<isleft) {
	wasleft=1;
	x = -xwidth_c/2.0 - rand01()*xwidth_t;
      } else {
	wasleft=0;
	x = xwidth_c/2.0 + rand01()*xwidth_t;
      }
      z = 0;
      cosine = 0;
      w_geom = w_geom_t;
    }
  } else if (strcasestr(sourcedef,"2015")) {
    if (cold) { //case: cold moderator
      xprime = -(0.06 + 0.1*rand01());
      yprime = 0.5*randpm1()*yheight_c;
      zprime = TSC2015_z0_BF3cm(100*xprime)/100.0;
      cosine = 0;
      w_geom = w_geom_c;
    }  else  { //case: thermal moderator
      xprime = -(0.09*randpm1());
      yprime = 0.5*randpm1()*yheight_t;
      zprime = TSC2015_z0_BF3cm(100*xprime)/100.0;
      cosine = 0;
      w_geom = w_geom_t;
    }
    y = yprime;
    x = sgnport*xprime*cos_beamport_angle + zprime*sin_beamport_angle;
    z = cos_beamport_angle*zprime - sgnport*xprime*sin_beamport_angle;
	/* the position was projected in the beamport direction - we must multiply by cosine*/
	/* if (xprime>edge_thermal) { */
	/* 	w_geom *=cos_thermal; */
	/* } else { */
	/* 	w_geom *=cos_cold; */
	/* } */


  } else {
    if (cold) {          //case: cold moderator
      //choose random point on cylinder surface
      double theta_tmp = (rand01()*extraction_opening - abs_beamport_angle)*DEG2RAD;
      x     = cyl_radius * sin(theta_tmp);
      y     = 0.5*randpm1()*yheight_c;
      z     = cyl_radius * cos(theta_tmp);
      cosine = 0;
      w_geom = w_geom_c;
    }  else  {                      //case: thermal moderator
      double poshorz, posvert;
      poshorz = cyl_radius+rand01()*xwidth_t;
      posvert = 0.5*randpm1()*yheight_t;

      if (rand01()<isleft) {
	wasleft=1;
	x = t1x * poshorz;
	z = t1z * poshorz;
      } else {
	wasleft=0;
	x = t2x * poshorz;
	z = t2z * poshorz;
      }
      y = posvert;
      cosine = 1;
      w_geom = w_geom_t;
    }
    // Correct for comp origin
    x = x-cx;
    y = y-cy;
    z = z-cz;
  }

  randvec_target_rect_real(&xf, &yf, &r, &w_focus,
			   tx, ty, tz, focus_xw, focus_yh, ROT_A_CURRENT_COMP, x, y, z, cosine);

  /* In case of 2014 source, solid angle correction is done internally  in the brill fcts */
  if (strcasestr(sourcedef,"2014") || strcasestr(sourcedef,"2015")) {
    w_focus=focus_xw*focus_yh/(tx*tx+ty*ty+tz*tz);
  }

  dx = xf-x;
  dy = yf-y;
  dz = r-z;
  r = sqrt(dx*dx+dy*dy+dz*dz);

  lambda = Lmin+l_range*rand01();    /* Choose from uniform distribution */

  k = 2*PI/lambda;
  v = K2V*k;

  vz = v*dz/r;
  vy = v*dy/r;
  vx = v*dx/r;

  if (strcasestr(sourcedef,"2015")) {
    modextras.X=xprime; modextras.Y=yprime; modextras.Z=zprime; //modextras.Wasleft=1;
  } else {
    modextras.X=x; modextras.Y=y; modextras.Z=z; modextras.Wasleft=wasleft;
  }

  if (cold) {          //case: cold moderator
    cold_bril( &t,  &p,  lambda,  tfocus_width,  tfocus_time,  dt, modextras);
  }  else  {                      //case: thermal moderator
    thermal_bril( &t,  &p,  lambda,  tfocus_width,  tfocus_time,  dt, modextras);
  }

  p*=w_stat*w_focus*w_geom*w_mult;
  t+=(double)floor((n_pulses)*rand01())/ESS_SOURCE_FREQUENCY;   /* Select a random pulse */

  /* Correct weight for sampling of cold vs. thermal events. */
  if (cold) {
    p /=cold_frac;
  } else {
    p/=(1-cold_frac);
    // In the 2015 case we only describe one thermal emission area
    if (!strcasestr(sourcedef,"2015")) {
      if (wasleft) {
	p/=(isleft);
      } else {
	p/=(1-isleft);
      }
    }
  }
  SCATTER;


  // EXTEND code here
  if (!strcmp(_comp->_name, "Source")) {
  IsCold=cold;
  }

  #undef isleft
  #undef Lmin
  #undef Lmax
  #undef cold_frac
  #undef dist
  #undef focus_xw
  #undef focus_yh
  #undef target_index
  #undef tmax_multiplier
  #undef yheight_c
  #undef yheight_t
  #undef n_pulses
  #undef acc_power
  #undef beamport_angle
  #undef sourcedef
  #undef xwidth_c
  #undef xwidth_t
  #undef extraction_opening
  #undef l_range
  #undef w_mult
  #undef w_geom
  #undef w_geom_c
  #undef w_geom_t
  #undef w_stat
  #undef cold
  #undef wasleft
  #undef cosine
  #undef tx
  #undef ty
  #undef tz
  #undef cx
  #undef cy
  #undef cz
  #undef cyl_radius
  #undef t1x
  #undef t1y
  #undef t1z
  #undef t2x
  #undef t2y
  #undef t2z
  #undef modextras
  #undef cold_bril
  #undef thermal_bril
  #undef tfocus_width
  #undef tfocus_time
  #undef dt
  #undef r_empty
  #undef lambda
  #undef xprime
  #undef yprime
  #undef zprime
  #undef sgnport
  #undef abs_beamport_angle
  #undef cos_beamport_angle
  #undef edge_thermal
  #undef cos_thermal
  #undef cos_cold
  #undef sin_beamport_angle
  return(_comp);
} /* class_ESS_moderator_trace */

#pragma acc routine seq
_class_PSDlin_monitor *class_PSDlin_monitor_trace(_class_PSDlin_monitor *_comp
  , _class_particle *_particle) {
  ABSORBED=SCATTERED=RESTORE=0;

  #define nx (_comp->_parameters.nx)
  #define filename (_comp->_parameters.filename)
  #define xmin (_comp->_parameters.xmin)
  #define xmax (_comp->_parameters.xmax)
  #define ymin (_comp->_parameters.ymin)
  #define ymax (_comp->_parameters.ymax)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define restore_neutron (_comp->_parameters.restore_neutron)
  #define PSDlin_N (_comp->_parameters.PSDlin_N)
  #define PSDlin_p (_comp->_parameters.PSDlin_p)
  #define PSDlin_p2 (_comp->_parameters.PSDlin_p2)
  int i;

  PROP_Z0;
  if (x>xmin && x<xmax && y>ymin && y<ymax)
  {
    i = floor(nx*(x-xmin)/(xmax-xmin));              /* Bin number */
    if((i >= nx) || (i<0))
    {
      printf("FATAL ERROR: wrong positioning in linear PSD. i= %i \n",i);
      exit(1);
    }
    PSDlin_N[i]++;
    PSDlin_p[i] += p;
    PSDlin_p2[i] += p*p;
  }
  if (restore_neutron) {
    RESTORE_NEUTRON(INDEX_CURRENT_COMP, x, y, z, vx, vy, vz, t, sx, sy, sz, p);
  }
  #undef nx
  #undef filename
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef xwidth
  #undef yheight
  #undef restore_neutron
  #undef PSDlin_N
  #undef PSDlin_p
  #undef PSDlin_p2
  return(_comp);
} /* class_PSDlin_monitor_trace */

#pragma acc routine seq
_class_TOF_monitor *class_TOF_monitor_trace(_class_TOF_monitor *_comp
  , _class_particle *_particle) {
  ABSORBED=SCATTERED=RESTORE=0;

  #define nt (_comp->_parameters.nt)
  #define filename (_comp->_parameters.filename)
  #define xmin (_comp->_parameters.xmin)
  #define xmax (_comp->_parameters.xmax)
  #define ymin (_comp->_parameters.ymin)
  #define ymax (_comp->_parameters.ymax)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define tmin (_comp->_parameters.tmin)
  #define tmax (_comp->_parameters.tmax)
  #define dt (_comp->_parameters.dt)
  #define restore_neutron (_comp->_parameters.restore_neutron)
  #define TOF_N (_comp->_parameters.TOF_N)
  #define TOF_p (_comp->_parameters.TOF_p)
  #define TOF_p2 (_comp->_parameters.TOF_p2)
  #define t_min (_comp->_parameters.t_min)
  #define t_max (_comp->_parameters.t_max)
  #define delta_t (_comp->_parameters.delta_t)
  int i;

  PROP_Z0;
  if (x>xmin && x<xmax && y>ymin && y<ymax)
  {
    i = floor((1E6*t-t_min)/delta_t);              /* Bin number */
    if(i >= 0 && i < nt) {
      TOF_N[i]++;
      TOF_p[i] += p;
      TOF_p2[i] += p*p;
      SCATTER;
    }
  }
  if (restore_neutron) {
    RESTORE_NEUTRON(INDEX_CURRENT_COMP, x, y, z, vx, vy, vz, t, sx, sy, sz, p);
  }
  #undef nt
  #undef filename
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef xwidth
  #undef yheight
  #undef tmin
  #undef tmax
  #undef dt
  #undef restore_neutron
  #undef TOF_N
  #undef TOF_p
  #undef TOF_p2
  #undef t_min
  #undef t_max
  #undef delta_t
  return(_comp);
} /* class_TOF_monitor_trace */

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

/* *****************************************************************************
* instrument 'ESS_2015_test' TRACE
***************************************************************************** */

#pragma acc routine seq
int raytrace(_class_particle* _particle) { /* called by mccode_main for ESS_2015_test:TRACE */

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
      /* component Source=ESS_moderator() [2] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Source->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Source->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Source->_position_relative, _Source->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_ESS_moderator_trace(_Source, _particle); /* contains EXTEND code */
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component Source [2] */
    if (!ABSORBED && _particle->_index == 3) {
      /* component PortOrig=Arm() [3] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_PortOrig->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _PortOrig->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_PortOrig->_position_relative, _PortOrig->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      _particle->_index++;
    } /* end component PortOrig [3] */
    if (!ABSORBED && _particle->_index == 4) {
      /* component XmonT=PSDlin_monitor() [4] */
  #define nx (_XmonT->_parameters.nx)
  #define filename (_XmonT->_parameters.filename)
  #define xmin (_XmonT->_parameters.xmin)
  #define xmax (_XmonT->_parameters.xmax)
  #define ymin (_XmonT->_parameters.ymin)
  #define ymax (_XmonT->_parameters.ymax)
  #define xwidth (_XmonT->_parameters.xwidth)
  #define yheight (_XmonT->_parameters.yheight)
  #define restore_neutron (_XmonT->_parameters.restore_neutron)
  #define PSDlin_N (_XmonT->_parameters.PSDlin_N)
  #define PSDlin_p (_XmonT->_parameters.PSDlin_p)
  #define PSDlin_p2 (_XmonT->_parameters.PSDlin_p2)
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_XmonT->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _XmonT->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_XmonT->_position_relative, _XmonT->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      if (( ! IsCold )) // conditional WHEN execution
      class_PSDlin_monitor_trace(_XmonT, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
  #undef nx
  #undef filename
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef xwidth
  #undef yheight
  #undef restore_neutron
  #undef PSDlin_N
  #undef PSDlin_p
  #undef PSDlin_p2
      _particle->_index++;
    } /* end component XmonT [4] */
    if (!ABSORBED && _particle->_index == 5) {
      /* component XmonC=PSDlin_monitor() [5] */
  #define nx (_XmonC->_parameters.nx)
  #define filename (_XmonC->_parameters.filename)
  #define xmin (_XmonC->_parameters.xmin)
  #define xmax (_XmonC->_parameters.xmax)
  #define ymin (_XmonC->_parameters.ymin)
  #define ymax (_XmonC->_parameters.ymax)
  #define xwidth (_XmonC->_parameters.xwidth)
  #define yheight (_XmonC->_parameters.yheight)
  #define restore_neutron (_XmonC->_parameters.restore_neutron)
  #define PSDlin_N (_XmonC->_parameters.PSDlin_N)
  #define PSDlin_p (_XmonC->_parameters.PSDlin_p)
  #define PSDlin_p2 (_XmonC->_parameters.PSDlin_p2)
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_XmonC->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _XmonC->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_XmonC->_position_relative, _XmonC->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      if (( IsCold )) // conditional WHEN execution
      class_PSDlin_monitor_trace(_XmonC, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
  #undef nx
  #undef filename
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef xwidth
  #undef yheight
  #undef restore_neutron
  #undef PSDlin_N
  #undef PSDlin_p
  #undef PSDlin_p2
      _particle->_index++;
    } /* end component XmonC [5] */
    if (!ABSORBED && _particle->_index == 6) {
      /* component tofT=TOF_monitor() [6] */
  #define nt (_tofT->_parameters.nt)
  #define filename (_tofT->_parameters.filename)
  #define xmin (_tofT->_parameters.xmin)
  #define xmax (_tofT->_parameters.xmax)
  #define ymin (_tofT->_parameters.ymin)
  #define ymax (_tofT->_parameters.ymax)
  #define xwidth (_tofT->_parameters.xwidth)
  #define yheight (_tofT->_parameters.yheight)
  #define tmin (_tofT->_parameters.tmin)
  #define tmax (_tofT->_parameters.tmax)
  #define dt (_tofT->_parameters.dt)
  #define restore_neutron (_tofT->_parameters.restore_neutron)
  #define TOF_N (_tofT->_parameters.TOF_N)
  #define TOF_p (_tofT->_parameters.TOF_p)
  #define TOF_p2 (_tofT->_parameters.TOF_p2)
  #define t_min (_tofT->_parameters.t_min)
  #define t_max (_tofT->_parameters.t_max)
  #define delta_t (_tofT->_parameters.delta_t)
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_tofT->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _tofT->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_tofT->_position_relative, _tofT->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      if (( ! IsCold )) // conditional WHEN execution
      class_TOF_monitor_trace(_tofT, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
  #undef nt
  #undef filename
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef xwidth
  #undef yheight
  #undef tmin
  #undef tmax
  #undef dt
  #undef restore_neutron
  #undef TOF_N
  #undef TOF_p
  #undef TOF_p2
  #undef t_min
  #undef t_max
  #undef delta_t
      _particle->_index++;
    } /* end component tofT [6] */
    if (!ABSORBED && _particle->_index == 7) {
      /* component tofC=TOF_monitor() [7] */
  #define nt (_tofC->_parameters.nt)
  #define filename (_tofC->_parameters.filename)
  #define xmin (_tofC->_parameters.xmin)
  #define xmax (_tofC->_parameters.xmax)
  #define ymin (_tofC->_parameters.ymin)
  #define ymax (_tofC->_parameters.ymax)
  #define xwidth (_tofC->_parameters.xwidth)
  #define yheight (_tofC->_parameters.yheight)
  #define tmin (_tofC->_parameters.tmin)
  #define tmax (_tofC->_parameters.tmax)
  #define dt (_tofC->_parameters.dt)
  #define restore_neutron (_tofC->_parameters.restore_neutron)
  #define TOF_N (_tofC->_parameters.TOF_N)
  #define TOF_p (_tofC->_parameters.TOF_p)
  #define TOF_p2 (_tofC->_parameters.TOF_p2)
  #define t_min (_tofC->_parameters.t_min)
  #define t_max (_tofC->_parameters.t_max)
  #define delta_t (_tofC->_parameters.delta_t)
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_tofC->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _tofC->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_tofC->_position_relative, _tofC->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      if (( IsCold )) // conditional WHEN execution
      class_TOF_monitor_trace(_tofC, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
  #undef nt
  #undef filename
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef xwidth
  #undef yheight
  #undef tmin
  #undef tmax
  #undef dt
  #undef restore_neutron
  #undef TOF_N
  #undef TOF_p
  #undef TOF_p2
  #undef t_min
  #undef t_max
  #undef delta_t
      _particle->_index++;
    } /* end component tofC [7] */
    if (!ABSORBED && _particle->_index == 8) {
      /* component Slit1=Slit() [8] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Slit1->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Slit1->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Slit1->_position_relative, _Slit1->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Slit_trace(_Slit1, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component Slit1 [8] */
    if (!ABSORBED && _particle->_index == 9) {
      /* component LamT=L_monitor() [9] */
  #define nL (_LamT->_parameters.nL)
  #define filename (_LamT->_parameters.filename)
  #define xmin (_LamT->_parameters.xmin)
  #define xmax (_LamT->_parameters.xmax)
  #define ymin (_LamT->_parameters.ymin)
  #define ymax (_LamT->_parameters.ymax)
  #define xwidth (_LamT->_parameters.xwidth)
  #define yheight (_LamT->_parameters.yheight)
  #define Lmin (_LamT->_parameters.Lmin)
  #define Lmax (_LamT->_parameters.Lmax)
  #define restore_neutron (_LamT->_parameters.restore_neutron)
  #define L_N (_LamT->_parameters.L_N)
  #define L_p (_LamT->_parameters.L_p)
  #define L_p2 (_LamT->_parameters.L_p2)
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_LamT->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _LamT->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_LamT->_position_relative, _LamT->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      if (( ! IsCold )) // conditional WHEN execution
      class_L_monitor_trace(_LamT, _particle);
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
    } /* end component LamT [9] */
    if (!ABSORBED && _particle->_index == 10) {
      /* component LamC=L_monitor() [10] */
  #define nL (_LamC->_parameters.nL)
  #define filename (_LamC->_parameters.filename)
  #define xmin (_LamC->_parameters.xmin)
  #define xmax (_LamC->_parameters.xmax)
  #define ymin (_LamC->_parameters.ymin)
  #define ymax (_LamC->_parameters.ymax)
  #define xwidth (_LamC->_parameters.xwidth)
  #define yheight (_LamC->_parameters.yheight)
  #define Lmin (_LamC->_parameters.Lmin)
  #define Lmax (_LamC->_parameters.Lmax)
  #define restore_neutron (_LamC->_parameters.restore_neutron)
  #define L_N (_LamC->_parameters.L_N)
  #define L_p (_LamC->_parameters.L_p)
  #define L_p2 (_LamC->_parameters.L_p2)
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_LamC->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _LamC->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_LamC->_position_relative, _LamC->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      if (( IsCold )) // conditional WHEN execution
      class_L_monitor_trace(_LamC, _particle);
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
    } /* end component LamC [10] */
    if (_particle->_index > 10)
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
* instrument 'ESS_2015_test' and components SAVE
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

_class_PSDlin_monitor *class_PSDlin_monitor_save(_class_PSDlin_monitor *_comp
) {
  #define nx (_comp->_parameters.nx)
  #define filename (_comp->_parameters.filename)
  #define xmin (_comp->_parameters.xmin)
  #define xmax (_comp->_parameters.xmax)
  #define ymin (_comp->_parameters.ymin)
  #define ymax (_comp->_parameters.ymax)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define restore_neutron (_comp->_parameters.restore_neutron)
  #define PSDlin_N (_comp->_parameters.PSDlin_N)
  #define PSDlin_p (_comp->_parameters.PSDlin_p)
  #define PSDlin_p2 (_comp->_parameters.PSDlin_p2)
  DETECTOR_OUT_1D(
      "Linear PSD monitor",
      "x-Position [m]",
      "Intensity",
      "x", xmin, xmax, nx,
      &PSDlin_N[0],&PSDlin_p[0],&PSDlin_p2[0],
      filename);
  #undef nx
  #undef filename
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef xwidth
  #undef yheight
  #undef restore_neutron
  #undef PSDlin_N
  #undef PSDlin_p
  #undef PSDlin_p2
  return(_comp);
} /* class_PSDlin_monitor_save */

_class_TOF_monitor *class_TOF_monitor_save(_class_TOF_monitor *_comp
) {
  #define nt (_comp->_parameters.nt)
  #define filename (_comp->_parameters.filename)
  #define xmin (_comp->_parameters.xmin)
  #define xmax (_comp->_parameters.xmax)
  #define ymin (_comp->_parameters.ymin)
  #define ymax (_comp->_parameters.ymax)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define tmin (_comp->_parameters.tmin)
  #define tmax (_comp->_parameters.tmax)
  #define dt (_comp->_parameters.dt)
  #define restore_neutron (_comp->_parameters.restore_neutron)
  #define TOF_N (_comp->_parameters.TOF_N)
  #define TOF_p (_comp->_parameters.TOF_p)
  #define TOF_p2 (_comp->_parameters.TOF_p2)
  #define t_min (_comp->_parameters.t_min)
  #define t_max (_comp->_parameters.t_max)
  #define delta_t (_comp->_parameters.delta_t)
  DETECTOR_OUT_1D(
      "Time-of-flight monitor",
      "Time-of-flight [\\gms]",
      "Intensity",
      "t", t_min, t_max, nt,
      &TOF_N[0],&TOF_p[0],&TOF_p2[0],
      filename);
  #undef nt
  #undef filename
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef xwidth
  #undef yheight
  #undef tmin
  #undef tmax
  #undef dt
  #undef restore_neutron
  #undef TOF_N
  #undef TOF_p
  #undef TOF_p2
  #undef t_min
  #undef t_max
  #undef delta_t
  return(_comp);
} /* class_TOF_monitor_save */

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



int save(FILE *handle) { /* called by mccode_main for ESS_2015_test:SAVE */
  if (!handle) siminfo_init(NULL);

  /* call iteratively all components SAVE */
  class_Progress_bar_save(_Origin);



  class_PSDlin_monitor_save(_XmonT);

  class_PSDlin_monitor_save(_XmonC);

  class_TOF_monitor_save(_tofT);

  class_TOF_monitor_save(_tofC);


  class_L_monitor_save(_LamT);

  class_L_monitor_save(_LamC);

  if (!handle) siminfo_close(); 

  return(0);
} /* save */

/* *****************************************************************************
* instrument 'ESS_2015_test' and components FINALLY
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

_class_PSDlin_monitor *class_PSDlin_monitor_finally(_class_PSDlin_monitor *_comp
) {
  #define nx (_comp->_parameters.nx)
  #define filename (_comp->_parameters.filename)
  #define xmin (_comp->_parameters.xmin)
  #define xmax (_comp->_parameters.xmax)
  #define ymin (_comp->_parameters.ymin)
  #define ymax (_comp->_parameters.ymax)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define restore_neutron (_comp->_parameters.restore_neutron)
  #define PSDlin_N (_comp->_parameters.PSDlin_N)
  #define PSDlin_p (_comp->_parameters.PSDlin_p)
  #define PSDlin_p2 (_comp->_parameters.PSDlin_p2)
  destroy_darr1d(PSDlin_N);
  destroy_darr1d(PSDlin_p);
  destroy_darr1d(PSDlin_p2);
  #undef nx
  #undef filename
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef xwidth
  #undef yheight
  #undef restore_neutron
  #undef PSDlin_N
  #undef PSDlin_p
  #undef PSDlin_p2
  return(_comp);
} /* class_PSDlin_monitor_finally */

_class_TOF_monitor *class_TOF_monitor_finally(_class_TOF_monitor *_comp
) {
  #define nt (_comp->_parameters.nt)
  #define filename (_comp->_parameters.filename)
  #define xmin (_comp->_parameters.xmin)
  #define xmax (_comp->_parameters.xmax)
  #define ymin (_comp->_parameters.ymin)
  #define ymax (_comp->_parameters.ymax)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define tmin (_comp->_parameters.tmin)
  #define tmax (_comp->_parameters.tmax)
  #define dt (_comp->_parameters.dt)
  #define restore_neutron (_comp->_parameters.restore_neutron)
  #define TOF_N (_comp->_parameters.TOF_N)
  #define TOF_p (_comp->_parameters.TOF_p)
  #define TOF_p2 (_comp->_parameters.TOF_p2)
  #define t_min (_comp->_parameters.t_min)
  #define t_max (_comp->_parameters.t_max)
  #define delta_t (_comp->_parameters.delta_t)
  destroy_darr1d(TOF_N);
  destroy_darr1d(TOF_p);
  destroy_darr1d(TOF_p2);
  #undef nt
  #undef filename
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef xwidth
  #undef yheight
  #undef tmin
  #undef tmax
  #undef dt
  #undef restore_neutron
  #undef TOF_N
  #undef TOF_p
  #undef TOF_p2
  #undef t_min
  #undef t_max
  #undef delta_t
  return(_comp);
} /* class_TOF_monitor_finally */

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



int finally(void) { /* called by mccode_main for ESS_2015_test:FINALLY */
  siminfo_init(NULL);
  save(siminfo_file); /* save data when simulation ends */

  /* call iteratively all components FINALLY */
  class_Progress_bar_finally(_Origin);



  class_PSDlin_monitor_finally(_XmonT);

  class_PSDlin_monitor_finally(_XmonC);

  class_TOF_monitor_finally(_tofT);

  class_TOF_monitor_finally(_tofC);


  class_L_monitor_finally(_LamT);

  class_L_monitor_finally(_LamC);

  siminfo_close(); 

  return(0);
} /* finally */

/* *****************************************************************************
* instrument 'ESS_2015_test' and components DISPLAY
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

_class_ESS_moderator *class_ESS_moderator_display(_class_ESS_moderator *_comp
) {
  #define isleft (_comp->_parameters.isleft)
  #define Lmin (_comp->_parameters.Lmin)
  #define Lmax (_comp->_parameters.Lmax)
  #define cold_frac (_comp->_parameters.cold_frac)
  #define dist (_comp->_parameters.dist)
  #define focus_xw (_comp->_parameters.focus_xw)
  #define focus_yh (_comp->_parameters.focus_yh)
  #define target_index (_comp->_parameters.target_index)
  #define tmax_multiplier (_comp->_parameters.tmax_multiplier)
  #define yheight_c (_comp->_parameters.yheight_c)
  #define yheight_t (_comp->_parameters.yheight_t)
  #define n_pulses (_comp->_parameters.n_pulses)
  #define acc_power (_comp->_parameters.acc_power)
  #define beamport_angle (_comp->_parameters.beamport_angle)
  #define sourcedef (_comp->_parameters.sourcedef)
  #define xwidth_c (_comp->_parameters.xwidth_c)
  #define xwidth_t (_comp->_parameters.xwidth_t)
  #define extraction_opening (_comp->_parameters.extraction_opening)
  #define l_range (_comp->_parameters.l_range)
  #define w_mult (_comp->_parameters.w_mult)
  #define w_geom (_comp->_parameters.w_geom)
  #define w_geom_c (_comp->_parameters.w_geom_c)
  #define w_geom_t (_comp->_parameters.w_geom_t)
  #define w_stat (_comp->_parameters.w_stat)
  #define cold (_comp->_parameters.cold)
  #define wasleft (_comp->_parameters.wasleft)
  #define cosine (_comp->_parameters.cosine)
  #define tx (_comp->_parameters.tx)
  #define ty (_comp->_parameters.ty)
  #define tz (_comp->_parameters.tz)
  #define cx (_comp->_parameters.cx)
  #define cy (_comp->_parameters.cy)
  #define cz (_comp->_parameters.cz)
  #define cyl_radius (_comp->_parameters.cyl_radius)
  #define t1x (_comp->_parameters.t1x)
  #define t1y (_comp->_parameters.t1y)
  #define t1z (_comp->_parameters.t1z)
  #define t2x (_comp->_parameters.t2x)
  #define t2y (_comp->_parameters.t2y)
  #define t2z (_comp->_parameters.t2z)
  #define modextras (_comp->_parameters.modextras)
  #define cold_bril (_comp->_parameters.cold_bril)
  #define thermal_bril (_comp->_parameters.thermal_bril)
  #define tfocus_width (_comp->_parameters.tfocus_width)
  #define tfocus_time (_comp->_parameters.tfocus_time)
  #define dt (_comp->_parameters.dt)
  #define r_empty (_comp->_parameters.r_empty)
  #define lambda (_comp->_parameters.lambda)
  #define xprime (_comp->_parameters.xprime)
  #define yprime (_comp->_parameters.yprime)
  #define zprime (_comp->_parameters.zprime)
  #define sgnport (_comp->_parameters.sgnport)
  #define abs_beamport_angle (_comp->_parameters.abs_beamport_angle)
  #define cos_beamport_angle (_comp->_parameters.cos_beamport_angle)
  #define edge_thermal (_comp->_parameters.edge_thermal)
  #define cos_thermal (_comp->_parameters.cos_thermal)
  #define cos_cold (_comp->_parameters.cos_cold)
  #define sin_beamport_angle (_comp->_parameters.sin_beamport_angle)
  /* if (planar==0) { */
  /* Draw cold moderator as cylinder */
  
  if (!strcasestr(sourcedef,"2015")) {
    circle("xz", cx,  cy+yheight_c/2.0, cz, cyl_radius);
    circle("xz", cx,  cy-yheight_c/2.0, cz, cyl_radius);
    line(cx, cy-yheight_c/2.0, cz+cyl_radius, cx, cy+yheight_c/2.0, cz+cyl_radius);
    line(cx, cy-yheight_c/2.0, cz-cyl_radius, cx, cy+yheight_c/2.0, cz-cyl_radius);
    line(cx+cyl_radius, cy-yheight_c/2.0, cz, cx+cyl_radius, cy+yheight_c/2.0, cz);
    line(cx-cyl_radius, cy-yheight_c/2.0, cz, cx-cyl_radius, cy+yheight_c/2.0, cz);

    /* Draw thermal moderators as a couple of squares + some lines */

    // Left
    multiline(5, cx+t1x*cyl_radius, cy-yheight_t/2.0, cz+t1z*cyl_radius,
  	      cx+t1x*cyl_radius + t1x*xwidth_t, cy-yheight_t/2.0, cz+t1z*cyl_radius + t1z*xwidth_t,
  	      cx+t1x*cyl_radius + t1x*xwidth_t, cy+yheight_t/2.0, cz+t1z*cyl_radius + t1z*xwidth_t,
  	      cx+t1x*cyl_radius, cy+yheight_t/2.0, cz+t1z*cyl_radius, cx+t1x*cyl_radius, cy-yheight_t/2.0, cz+t1z*cyl_radius);
    // Right
    multiline(5, cx+t2x*cyl_radius, cy-yheight_t/2.0, cz+t2z*cyl_radius,
	      cx+t2x*cyl_radius + t2x*xwidth_t, cy-yheight_t/2.0, cz+t2z*cyl_radius + t2z*xwidth_t,
	      cx+t2x*cyl_radius + t2x*xwidth_t, cy+yheight_t/2.0, cz+t2z*cyl_radius + t2z*xwidth_t,
	      cx+t2x*cyl_radius, cy+yheight_t/2.0, cz+t2z*cyl_radius, cx+t2x*cyl_radius, cy-yheight_t/2.0, cz+t2z*cyl_radius);

    /* Dashed lines for indicating "beam extraction" area... */
    dashed_line(cx+t1x*cyl_radius,cy-yheight_c/2.0, cz+t1z*cyl_radius, cx+t1x*r_empty, cy-yheight_c/2.0, cz+t1z*r_empty,10);
    dashed_line(cx+t1x*cyl_radius, cy+yheight_c/2.0, cz+t1z*cyl_radius, cx+t1x*r_empty, cy+yheight_c/2.0, cz+t1z*r_empty,10);
    dashed_line(cx+t2x*cyl_radius, cy-yheight_c/2.0, cz+t2z*cyl_radius, cx+t2x*r_empty, cy-yheight_c/2.0, cz+t2z*r_empty,5);
    dashed_line(cx+t2x*cyl_radius, cy+yheight_c/2.0, cz+t2z*cyl_radius, cx+t2x*r_empty, cy+yheight_c/2.0, cz+t2z*r_empty,5);


    /* Circles indicating extent of the "empty" zone where optics is not allowed */
    circle("xz", cx,  cy+yheight_c/2.0, cz, r_empty);
    circle("xz", cx,  cy-yheight_c/2.0, cz, r_empty);

    /* Circles indicating the builk shielding of the target monolith at 6 m */
    circle("xz", cx,  cy+focus_yh/2.0 , cz, 6);
    circle("xz", cx, cy-focus_yh/2.0 , cz, 6);
    circle("xz", cx,  cy+2, cz, 6);
    circle("xz", cx, cy-2, cz, 6);
  } else { // 2015
    double dx = (0.09 + 0.16)/25;
    double xp0, yp0, zp0, xp1, yp1, zp1;
    double xr0, yr0, zr0, xr1, yr1, zr1;
    double dxp0, dzp0, dxr0, dzr0, dxr1, dzr1;
    double ax,az,bbx,bbz,ccx,ccz;
    /* Drawing the "emission planes" of the cold and thermal moderators */
    xp0 = -0.16;
    yp0 = -yheight_c/2.0;
    int k,j;
    for (k=0; k<2; k++) {
      for (j=0; j<26; j++) {
	xp1 = xp0 + dx;
	yp1 = yp0;
	zp0 = TSC2015_z0_BF3cm(100*xp0)/100.0;
	zp1 = TSC2015_z0_BF3cm(100*xp1)/100.0;
	xr0 = cos(-abs_beamport_angle*DEG2RAD)*xp0 - sin(-abs_beamport_angle*DEG2RAD)*zp0;
	zr0 = cos(-abs_beamport_angle*DEG2RAD)*zp0 + sin(-abs_beamport_angle*DEG2RAD)*xp0;
	xr1 = cos(-abs_beamport_angle*DEG2RAD)*xp1 - sin(-abs_beamport_angle*DEG2RAD)*zp1;
	zr1 = cos(-abs_beamport_angle*DEG2RAD)*zp1 + sin(-abs_beamport_angle*DEG2RAD)*xp1;
	if (j==0 || j==8 || j==11) {
	  line(xr0, -yheight_c/2.0, zr0, xr0, yheight_c/2.0, zr0);
	}
	if (j==25) {
	  line(xr1, -yheight_c/2.0, zr1, xr1, yheight_c/2.0, zr1);
	}
	if (j%2==0) { // Every other element skipped to make dashed lines...
	  line(xr0, yp0, zr0, xr1, yp1, zr1);
	}
	xp0 = xp1;
      }
      /* Drawing a sketch-version of the butterfly "cross" and wing(s).... */
      /* Thermal cross: */
      xp0 = 0.0827; xp1 = -xp0;
      zp0 = TSC2015_z0_BF3cm(100*xp0)/100.0;
      zp1 = zp0;
      dxp0 = 0.01; dzp0=0.01;
      xr0 = cos(-abs_beamport_angle*DEG2RAD)*xp0 - sin(-abs_beamport_angle*DEG2RAD)*zp0;
      zr0 = cos(-abs_beamport_angle*DEG2RAD)*zp0 + sin(-abs_beamport_angle*DEG2RAD)*xp0;
      xr1 = cos(-abs_beamport_angle*DEG2RAD)*xp1 - sin(-abs_beamport_angle*DEG2RAD)*zp1;
      zr1 = cos(-abs_beamport_angle*DEG2RAD)*zp1 + sin(-abs_beamport_angle*DEG2RAD)*xp1;
      dxr0 = cos(-abs_beamport_angle*DEG2RAD)*dxp0;
      dzr0 = sin(-abs_beamport_angle*DEG2RAD)*dxp0;
      dxr1 = - sin(-abs_beamport_angle*DEG2RAD)*dzp0;
      dzr1 = cos(-abs_beamport_angle*DEG2RAD)*dzp0;

      line(xr0+dxr0, yp0, zr0+dzr0, -xr0+dxr0, yp0, -zr0+dzr0);
      line(xr1+dxr0, yp0, zr1+dzr0, -xr1+dxr0, yp0, -zr1+dzr0);
      line(xr0-dxr0, yp0, zr0-dzr0, -xr0-dxr0, yp0, -zr0-dzr0);
      line(xr1-dxr0, yp0, zr1-dzr0, -xr1-dxr0, yp0, -zr1-dzr0);
      line(xr0+dxr0, yp0, zr0+dzr0, xr0-dxr0, yp0, zr0-dzr0);
      line(xr1+dxr0, yp0, zr1+dzr0, xr1-dxr0, yp0, zr1-dzr0);
      line(-xr0+dxr0, yp0, -zr0+dzr0, -xr0-dxr0, yp0, -zr0-dzr0);
      line(-xr1+dxr0, yp0, -zr1+dzr0, -xr1-dxr0, yp0, -zr1-dzr0);

      /* Wings: */
      line(xr0+2*dxr0, yp0, zr0+2*dzr0, 2*dxr0, yp0, 2*dzr0);
      line(-xr1+2*dxr0, yp0, -zr1+2*dzr0, 2*dxr0, yp0, 2*dzr0);
      line(xr0+2*dxr0, yp0, zr0+2*dzr0,-xr1+2*dxr0, yp0, -zr1+2*dzr0);

      line(-xr0-2*dxr0, yp0, -zr0-2*dzr0, -2*dxr0, yp0, -2*dzr0);
      line(xr1-2*dxr0, yp0, zr1-2*dzr0, -2*dxr0, yp0, -2*dzr0);
      line(xr1-2*dxr0, yp0, zr1-2*dzr0,-xr0-2*dxr0, yp0, -zr0-2*dzr0);
      /* Back to drawing the emission planes... */
      yp1 = yp0;
      yp0 = yp0 + yheight_c;
      xp0 = -0.16;

    }

    /* Dashed lines for indicating "beam extraction" area... */
    dashed_line(cx+t1x*0.2,cy-yheight_c/2.0, cz+t1z*0.2, cx+t1x*r_empty, cy-yheight_c/2.0, cz+t1z*r_empty,10);
    dashed_line(cx+t1x*0.2, cy+yheight_c/2.0, cz+t1z*0.2, cx+t1x*r_empty, cy+yheight_c/2.0, cz+t1z*r_empty,10);
    dashed_line(cx+t2x*0.2, cy-yheight_c/2.0, cz+t2z*0.2, cx+t2x*r_empty, cy-yheight_c/2.0, cz+t2z*r_empty,5);
    dashed_line(cx+t2x*0.2, cy+yheight_c/2.0, cz+t2z*0.2, cx+t2x*r_empty, cy+yheight_c/2.0, cz+t2z*r_empty,5);

    /* Rectangle indicating the chosen focus rectangle - where the optics starts... */
    rectangle("xy",tx,ty,tz,focus_xw,focus_yh);

    /* Circles indicating extent of the "empty" zone where optics is not allowed */
    circle("xz", cx,  cy+yheight_c/2.0, cz, r_empty);
    circle("xz", cx,  cy-yheight_c/2.0, cz, r_empty);

    /* Circles indicating the builk shielding of the target monolith at 6 m */
    circle("xz", cx,  cy+focus_yh/2.0 , cz, 6);
    circle("xz", cx, cy-focus_yh/2.0 , cz, 6);
    circle("xz", cx,  cy+2, cz, 6);
    circle("xz", cx, cy-2, cz, 6);

    /* Arrow indicating proton beam direction */
    ax=sin(DEG2RAD*(90+abs_beamport_angle));az=cos(DEG2RAD*(90+abs_beamport_angle));
    bbx=sin(DEG2RAD*(110+abs_beamport_angle));bbz=cos(DEG2RAD*(110+abs_beamport_angle));
    ccx=sin(DEG2RAD*(70+abs_beamport_angle));ccz=cos(DEG2RAD*(70+abs_beamport_angle));
    line(-0.15*ax,0,-0.15*az,-5.9*ax,0,-5.9*az);
    line(-0.15*ax,0,-0.15*az,-0.15*ax-0.1*bbx,0,-0.15*az-0.1*bbz);
    line(-0.15*ax,0,-0.15*az,-0.15*ax-0.1*ccx,0,-0.15*az-0.1*ccz);

    /* Rectangle indicating the chosen focus rectangle - where the optics starts... */
    rectangle("xy",tx,ty,tz,focus_xw,focus_yh);
  }

  /* This last bit is relevant only in connection with 2014 source definition */
  if (strcasestr(sourcedef,"2014")) {

    line(cx-5.9,cy,cz,cx-cyl_radius-0.1,cy,cz);
    line(cx-cyl_radius-0.2*cos(10*DEG2RAD),cy,cz+0.2*sin(10*DEG2RAD),cx-cyl_radius-0.1,cy,cz);
    line(cx-cyl_radius-0.2*cos(10*DEG2RAD),cy,cz-0.2*sin(10*DEG2RAD),cx-cyl_radius-0.1,cy,cz);

    /* Rectangle indicating the chosen focus rectangle - where the optics starts... */
    rectangle("xy",tx,ty,tz,focus_xw,focus_yh);
    rectangle("xy",0,0,0,xwidth_c,yheight_c);
    rectangle("xy",-xwidth_c/2-xwidth_t/2,0,0,xwidth_t,yheight_t);
    rectangle("xy",xwidth_c/2+xwidth_t/2,0,0,xwidth_t,yheight_t);
  }

  #undef isleft
  #undef Lmin
  #undef Lmax
  #undef cold_frac
  #undef dist
  #undef focus_xw
  #undef focus_yh
  #undef target_index
  #undef tmax_multiplier
  #undef yheight_c
  #undef yheight_t
  #undef n_pulses
  #undef acc_power
  #undef beamport_angle
  #undef sourcedef
  #undef xwidth_c
  #undef xwidth_t
  #undef extraction_opening
  #undef l_range
  #undef w_mult
  #undef w_geom
  #undef w_geom_c
  #undef w_geom_t
  #undef w_stat
  #undef cold
  #undef wasleft
  #undef cosine
  #undef tx
  #undef ty
  #undef tz
  #undef cx
  #undef cy
  #undef cz
  #undef cyl_radius
  #undef t1x
  #undef t1y
  #undef t1z
  #undef t2x
  #undef t2y
  #undef t2z
  #undef modextras
  #undef cold_bril
  #undef thermal_bril
  #undef tfocus_width
  #undef tfocus_time
  #undef dt
  #undef r_empty
  #undef lambda
  #undef xprime
  #undef yprime
  #undef zprime
  #undef sgnport
  #undef abs_beamport_angle
  #undef cos_beamport_angle
  #undef edge_thermal
  #undef cos_thermal
  #undef cos_cold
  #undef sin_beamport_angle
  return(_comp);
} /* class_ESS_moderator_display */

_class_Arm *class_Arm_display(_class_Arm *_comp
) {
  /* A bit ugly; hard-coded dimensions. */
  
  line(0,0,0,0.2,0,0);
  line(0,0,0,0,0.2,0);
  line(0,0,0,0,0,0.2);
  return(_comp);
} /* class_Arm_display */

_class_PSDlin_monitor *class_PSDlin_monitor_display(_class_PSDlin_monitor *_comp
) {
  #define nx (_comp->_parameters.nx)
  #define filename (_comp->_parameters.filename)
  #define xmin (_comp->_parameters.xmin)
  #define xmax (_comp->_parameters.xmax)
  #define ymin (_comp->_parameters.ymin)
  #define ymax (_comp->_parameters.ymax)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define restore_neutron (_comp->_parameters.restore_neutron)
  #define PSDlin_N (_comp->_parameters.PSDlin_N)
  #define PSDlin_p (_comp->_parameters.PSDlin_p)
  #define PSDlin_p2 (_comp->_parameters.PSDlin_p2)
  multiline(5, (double)xmin, (double)ymin, 0.0,
               (double)xmax, (double)ymin, 0.0,
               (double)xmax, (double)ymax, 0.0,
               (double)xmin, (double)ymax, 0.0,
               (double)xmin, (double)ymin, 0.0);
  #undef nx
  #undef filename
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef xwidth
  #undef yheight
  #undef restore_neutron
  #undef PSDlin_N
  #undef PSDlin_p
  #undef PSDlin_p2
  return(_comp);
} /* class_PSDlin_monitor_display */

_class_TOF_monitor *class_TOF_monitor_display(_class_TOF_monitor *_comp
) {
  #define nt (_comp->_parameters.nt)
  #define filename (_comp->_parameters.filename)
  #define xmin (_comp->_parameters.xmin)
  #define xmax (_comp->_parameters.xmax)
  #define ymin (_comp->_parameters.ymin)
  #define ymax (_comp->_parameters.ymax)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define tmin (_comp->_parameters.tmin)
  #define tmax (_comp->_parameters.tmax)
  #define dt (_comp->_parameters.dt)
  #define restore_neutron (_comp->_parameters.restore_neutron)
  #define TOF_N (_comp->_parameters.TOF_N)
  #define TOF_p (_comp->_parameters.TOF_p)
  #define TOF_p2 (_comp->_parameters.TOF_p2)
  #define t_min (_comp->_parameters.t_min)
  #define t_max (_comp->_parameters.t_max)
  #define delta_t (_comp->_parameters.delta_t)

  multiline(5, (double)xmin, (double)ymin, 0.0,
               (double)xmax, (double)ymin, 0.0,
               (double)xmax, (double)ymax, 0.0,
               (double)xmin, (double)ymax, 0.0,
               (double)xmin, (double)ymin, 0.0);
  #undef nt
  #undef filename
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef xwidth
  #undef yheight
  #undef tmin
  #undef tmax
  #undef dt
  #undef restore_neutron
  #undef TOF_N
  #undef TOF_p
  #undef TOF_p2
  #undef t_min
  #undef t_max
  #undef delta_t
  return(_comp);
} /* class_TOF_monitor_display */

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


  #undef magnify
  #undef line
  #undef dashed_line
  #undef multiline
  #undef rectangle
  #undef box
  #undef circle
  #undef cylinder
  #undef sphere

int display(void) { /* called by mccode_main for ESS_2015_test:DISPLAY */
  printf("MCDISPLAY: start\n");

  /* call iteratively all components DISPLAY */
  class_Progress_bar_display(_Origin);

  class_ESS_moderator_display(_Source);

  class_Arm_display(_PortOrig);

  class_PSDlin_monitor_display(_XmonT);

  class_PSDlin_monitor_display(_XmonC);

  class_TOF_monitor_display(_tofT);

  class_TOF_monitor_display(_tofC);

  class_Slit_display(_Slit1);

  class_L_monitor_display(_LamT);

  class_L_monitor_display(_LamC);

  printf("MCDISPLAY: end\n");

  return(0);
} /* display */

void* _getvar_parameters(char* compname)
/* enables settings parameters based use of the GETPAR macro */
{
  if (!strcmp(compname, "Origin")) return (void *) &(_Origin->_parameters);
  if (!strcmp(compname, "Source")) return (void *) &(_Source->_parameters);
  if (!strcmp(compname, "PortOrig")) return (void *) &(_PortOrig->_parameters);
  if (!strcmp(compname, "XmonT")) return (void *) &(_XmonT->_parameters);
  if (!strcmp(compname, "XmonC")) return (void *) &(_XmonC->_parameters);
  if (!strcmp(compname, "tofT")) return (void *) &(_tofT->_parameters);
  if (!strcmp(compname, "tofC")) return (void *) &(_tofC->_parameters);
  if (!strcmp(compname, "Slit1")) return (void *) &(_Slit1->_parameters);
  if (!strcmp(compname, "LamT")) return (void *) &(_LamT->_parameters);
  if (!strcmp(compname, "LamC")) return (void *) &(_LamC->_parameters);
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

/* end of generated C code ESS_2015_test.c */
