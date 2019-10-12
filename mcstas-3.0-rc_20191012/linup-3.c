/* Automatically generated file. Do not edit. 
 * Format:     ANSI C source code
 * Creator:    McStas <http://www.mcstas.org>
 * Instrument: linup-3.instr (TAS1_Diff_Slit)
 * Date:       Sat Oct 12 09:05:06 2019
 * File:       linup-3.c
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
* Start of instrument 'TAS1_Diff_Slit' generated code
***************************************************************************** */

#ifdef MC_TRACE_ENABLED
int traceenabled = 1;
#else
int traceenabled = 0;
#endif
#define MCSTAS "/usr/share/mcstas/3.0-dev/"
int   defaultmain         = 1;
char  instrument_name[]   = "TAS1_Diff_Slit";
char  instrument_source[] = "linup-3.instr";
char *instrument_exe      = NULL; /* will be set to argv[0] in main */
char  instrument_code[]   = "Instrument TAS1_Diff_Slit source code linup-3.instr is not embedded in this executable.\n  Use --source option when running McStas.\n";

int main(int argc, char *argv[]){return mccode_main(argc, argv);}

/* *****************************************************************************
* instrument 'TAS1_Diff_Slit' and components DECLARE
***************************************************************************** */

/* Instrument parameters: structure and a table for the initialisation
   (Used in e.g. inputparse and I/O function (e.g. detector_out) */

struct _struct_instrument_parameters {
  MCNUM _PHM;
  MCNUM _TTM;
  MCNUM _TT;
  MCNUM _C1;
  MCNUM _OMC1;
  MCNUM _C2;
  MCNUM _C3;
};
typedef struct _struct_instrument_parameters _class_instrument_parameters;

struct _instrument_struct {
  char   _name[256]; /* the name of this instrument e.g. 'TAS1_Diff_Slit' */
/* Counters per component instance */
  double counter_AbsorbProp[34]; /* absorbed events in PROP routines */
  double counter_N[34], counter_P[34], counter_P2[34]; /* event counters after each component instance */
  _class_particle _trajectory[34]; /* current trajectory for STORE/RESTORE */
/* Components position table (absolute and relative coords) */
  Coords _position_relative[34]; /* positions of all components */
  Coords _position_absolute[34];
  _class_instrument_parameters _parameters; /* instrument parameters */
} _instrument_var;
struct _instrument_struct *instrument = & _instrument_var;
#pragma acc declare create ( _instrument_var )
#pragma acc declare create ( instrument )

int numipar = 7;
struct mcinputtable_struct mcinputtable[] = {
  "PHM", &(_instrument_var._parameters._PHM), instr_type_double, "-37.077", 
  "TTM", &(_instrument_var._parameters._TTM), instr_type_double, "-74", 
  "TT", &(_instrument_var._parameters._TT), instr_type_double, "33.52", 
  "C1", &(_instrument_var._parameters._C1), instr_type_double, "30", 
  "OMC1", &(_instrument_var._parameters._OMC1), instr_type_double, "5.5", 
  "C2", &(_instrument_var._parameters._C2), instr_type_double, "28", 
  "C3", &(_instrument_var._parameters._C3), instr_type_double, "67", 
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

/* component source=Source_simple() [2] DECLARE */
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
  char     _name[256]; /* e.g. source */
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
_class_Source_simple _source_var;
_class_Source_simple *_source = &_source_var;
#pragma acc declare create ( _source_var )
#pragma acc declare create ( _source )

/* component slit1=Slit() [3] DECLARE */
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
  char     _name[256]; /* e.g. slit1 */
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
_class_Slit _slit1_var;
_class_Slit *_slit1 = &_slit1_var;
#pragma acc declare create ( _slit1_var )
#pragma acc declare create ( _slit1 )

_class_Slit _slit2_var;
_class_Slit *_slit2 = &_slit2_var;
#pragma acc declare create ( _slit2_var )
#pragma acc declare create ( _slit2 )

_class_Slit _slit3_var;
_class_Slit *_slit3 = &_slit3_var;
#pragma acc declare create ( _slit3_var )
#pragma acc declare create ( _slit3 )

_class_Arm _focus_mono_var;
_class_Arm *_focus_mono = &_focus_mono_var;
#pragma acc declare create ( _focus_mono_var )
#pragma acc declare create ( _focus_mono )

/* component m0=Monochromator_flat() [7] DECLARE */
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
  char     _name[256]; /* e.g. m0 */
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
_class_Monochromator_flat _m0_var;
_class_Monochromator_flat *_m0 = &_m0_var;
#pragma acc declare create ( _m0_var )
#pragma acc declare create ( _m0 )

_class_Monochromator_flat _m1_var;
_class_Monochromator_flat *_m1 = &_m1_var;
#pragma acc declare create ( _m1_var )
#pragma acc declare create ( _m1 )

_class_Monochromator_flat _m2_var;
_class_Monochromator_flat *_m2 = &_m2_var;
#pragma acc declare create ( _m2_var )
#pragma acc declare create ( _m2 )

_class_Monochromator_flat _m3_var;
_class_Monochromator_flat *_m3 = &_m3_var;
#pragma acc declare create ( _m3_var )
#pragma acc declare create ( _m3 )

_class_Monochromator_flat _m4_var;
_class_Monochromator_flat *_m4 = &_m4_var;
#pragma acc declare create ( _m4_var )
#pragma acc declare create ( _m4 )

_class_Monochromator_flat _m5_var;
_class_Monochromator_flat *_m5 = &_m5_var;
#pragma acc declare create ( _m5_var )
#pragma acc declare create ( _m5 )

_class_Monochromator_flat _m6_var;
_class_Monochromator_flat *_m6 = &_m6_var;
#pragma acc declare create ( _m6_var )
#pragma acc declare create ( _m6 )

_class_Monochromator_flat _m7_var;
_class_Monochromator_flat *_m7 = &_m7_var;
#pragma acc declare create ( _m7_var )
#pragma acc declare create ( _m7 )

_class_Arm _a2_var;
_class_Arm *_a2 = &_a2_var;
#pragma acc declare create ( _a2_var )
#pragma acc declare create ( _a2 )

_class_Slit _slitMS1_var;
_class_Slit *_slitMS1 = &_slitMS1_var;
#pragma acc declare create ( _slitMS1_var )
#pragma acc declare create ( _slitMS1 )

_class_Slit _slitMS2_var;
_class_Slit *_slitMS2 = &_slitMS2_var;
#pragma acc declare create ( _slitMS2_var )
#pragma acc declare create ( _slitMS2 )

/* component c1=Collimator_linear() [18] DECLARE */
/* Parameter definition for component type 'Collimator_linear' */
struct _struct_Collimator_linear_parameters {
  /* Component type 'Collimator_linear' setting parameters */
  MCNUM xmin;
  MCNUM xmax;
  MCNUM ymin;
  MCNUM ymax;
  MCNUM xwidth;
  MCNUM yheight;
  MCNUM length;
  MCNUM divergence;
  MCNUM transmission;
  MCNUM divergenceV;
  /* Component type 'Collimator_linear' private parameters */
  /* Component type 'Collimator_linear' DECLARE code stored as structure members */
double slope, slopeV;
}; /* _struct_Collimator_linear_parameters */
typedef struct _struct_Collimator_linear_parameters _class_Collimator_linear_parameters;

/* Parameters for component type 'Collimator_linear' */
struct _struct_Collimator_linear {
  char     _name[256]; /* e.g. c1 */
  char     _type[256]; /* Collimator_linear */
  long     _index; /* e.g. 2 index in TRACE list */
  Coords   _position_absolute;
  Coords   _position_relative; /* wrt PREVIOUS */
  Rotation _rotation_absolute;
  Rotation _rotation_relative; /* wrt PREVIOUS */
  int      _rotation_is_identity;
  _class_Collimator_linear_parameters _parameters;
};
typedef struct _struct_Collimator_linear _class_Collimator_linear;
_class_Collimator_linear _c1_var;
_class_Collimator_linear *_c1 = &_c1_var;
#pragma acc declare create ( _c1_var )
#pragma acc declare create ( _c1 )

_class_Slit _slitMS3_var;
_class_Slit *_slitMS3 = &_slitMS3_var;
#pragma acc declare create ( _slitMS3_var )
#pragma acc declare create ( _slitMS3 )

_class_Slit _slitMS4_var;
_class_Slit *_slitMS4 = &_slitMS4_var;
#pragma acc declare create ( _slitMS4_var )
#pragma acc declare create ( _slitMS4 )

_class_Slit _slitMS5_var;
_class_Slit *_slitMS5 = &_slitMS5_var;
#pragma acc declare create ( _slitMS5_var )
#pragma acc declare create ( _slitMS5 )

/* component mon=Monitor() [22] DECLARE */
/* Parameter definition for component type 'Monitor' */
struct _struct_Monitor_parameters {
  /* Component type 'Monitor' setting parameters */
  MCNUM xmin;
  MCNUM xmax;
  MCNUM ymin;
  MCNUM ymax;
  MCNUM xwidth;
  MCNUM yheight;
  MCNUM restore_neutron;
  /* Component type 'Monitor' private parameters */
  /* Component type 'Monitor' DECLARE code stored as structure members */
double Nsum;
double psum, p2sum;
}; /* _struct_Monitor_parameters */
typedef struct _struct_Monitor_parameters _class_Monitor_parameters;

/* Parameters for component type 'Monitor' */
struct _struct_Monitor {
  char     _name[256]; /* e.g. mon */
  char     _type[256]; /* Monitor */
  long     _index; /* e.g. 2 index in TRACE list */
  Coords   _position_absolute;
  Coords   _position_relative; /* wrt PREVIOUS */
  Rotation _rotation_absolute;
  Rotation _rotation_relative; /* wrt PREVIOUS */
  int      _rotation_is_identity;
  _class_Monitor_parameters _parameters;
};
typedef struct _struct_Monitor _class_Monitor;
_class_Monitor _mon_var;
_class_Monitor *_mon = &_mon_var;
#pragma acc declare create ( _mon_var )
#pragma acc declare create ( _mon )

/* component emon1=E_monitor() [23] DECLARE */
/* Parameter definition for component type 'E_monitor' */
struct _struct_E_monitor_parameters {
  /* Component type 'E_monitor' setting parameters */
  MCNUM nE;
  char filename[16384];
  MCNUM xmin;
  MCNUM xmax;
  MCNUM ymin;
  MCNUM ymax;
  MCNUM xwidth;
  MCNUM yheight;
  MCNUM Emin;
  MCNUM Emax;
  MCNUM restore_neutron;
  /* Component type 'E_monitor' private parameters */
  /* Component type 'E_monitor' DECLARE code stored as structure members */
  DArray1d E_N;
  DArray1d E_p;
  DArray1d E_p2;
  double S_p, S_pE, S_pE2;
}; /* _struct_E_monitor_parameters */
typedef struct _struct_E_monitor_parameters _class_E_monitor_parameters;

/* Parameters for component type 'E_monitor' */
struct _struct_E_monitor {
  char     _name[256]; /* e.g. emon1 */
  char     _type[256]; /* E_monitor */
  long     _index; /* e.g. 2 index in TRACE list */
  Coords   _position_absolute;
  Coords   _position_relative; /* wrt PREVIOUS */
  Rotation _rotation_absolute;
  Rotation _rotation_relative; /* wrt PREVIOUS */
  int      _rotation_is_identity;
  _class_E_monitor_parameters _parameters;
};
typedef struct _struct_E_monitor _class_E_monitor;
_class_E_monitor _emon1_var;
_class_E_monitor *_emon1 = &_emon1_var;
#pragma acc declare create ( _emon1_var )
#pragma acc declare create ( _emon1 )

_class_Slit _sample_var;
_class_Slit *_sample = &_sample_var;
#pragma acc declare create ( _sample_var )
#pragma acc declare create ( _sample )

_class_Slit _slit1mm_var;
_class_Slit *_slit1mm = &_slit1mm_var;
#pragma acc declare create ( _slit1mm_var )
#pragma acc declare create ( _slit1mm )

_class_Arm _a3_var;
_class_Arm *_a3 = &_a3_var;
#pragma acc declare create ( _a3_var )
#pragma acc declare create ( _a3 )

_class_Collimator_linear _c2_var;
_class_Collimator_linear *_c2 = &_c2_var;
#pragma acc declare create ( _c2_var )
#pragma acc declare create ( _c2 )

_class_Arm _ana_var;
_class_Arm *_ana = &_ana_var;
#pragma acc declare create ( _ana_var )
#pragma acc declare create ( _ana )

_class_Arm _a4_var;
_class_Arm *_a4 = &_a4_var;
#pragma acc declare create ( _a4_var )
#pragma acc declare create ( _a4 )

_class_Collimator_linear _c3_var;
_class_Collimator_linear *_c3 = &_c3_var;
#pragma acc declare create ( _c3_var )
#pragma acc declare create ( _c3 )

_class_Monitor _sng_var;
_class_Monitor *_sng = &_sng_var;
#pragma acc declare create ( _sng_var )
#pragma acc declare create ( _sng )

_class_E_monitor _Emon2_var;
_class_E_monitor *_Emon2 = &_Emon2_var;
#pragma acc declare create ( _Emon2_var )
#pragma acc declare create ( _Emon2 )

int mcNUMCOMP = 32;

/* User declarations from instrument definition. Can define functions. */
/* Mosaicity used on monochromator and analysator */
double tas1_mono_mosaic = 45; /* Measurements indicate its really 45' */
/* Q vector for bragg scattering with monochromator and analysator */
double tas1_mono_q = 2*1.87325; /* Fake 2nd order scattering for 20meV */
double tas1_mono_r0 = 0.6;

/* Collimators */
double OMC1_d;

double mpos0, mpos1, mpos2, mpos3, mpos4, mpos5, mpos6, mpos7;
double mrot0, mrot1, mrot2, mrot3, mrot4, mrot5, mrot6, mrot7;

#undef compcurname
#undef compcurtype
#undef compcurindex
/* end of instrument 'TAS1_Diff_Slit' and components DECLARE */

/* *****************************************************************************
* instrument 'TAS1_Diff_Slit' and components INITIALISE
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

/* component source=Source_simple() SETTING, POSITION/ROTATION */
int _source_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_source_setpos] component source=Source_simple() SETTING [/usr/share/mcstas/3.0-dev/sources/Source_simple.comp:64]");
  stracpy(_source->_name, "source", 16384);
  stracpy(_source->_type, "Source_simple", 16384);
  _source->_index=2;
  _source->_parameters.radius = 0.060;
  #define radius (_source->_parameters.radius)
  _source->_parameters.yheight = 0;
  #define yheight (_source->_parameters.yheight)
  _source->_parameters.xwidth = 0;
  #define xwidth (_source->_parameters.xwidth)
  _source->_parameters.dist = 3.288;
  #define dist (_source->_parameters.dist)
  _source->_parameters.focus_xw = 0.042;
  #define focus_xw (_source->_parameters.focus_xw)
  _source->_parameters.focus_yh = 0.082;
  #define focus_yh (_source->_parameters.focus_yh)
  _source->_parameters.E0 = 20;
  #define E0 (_source->_parameters.E0)
  _source->_parameters.dE = 0.82;
  #define dE (_source->_parameters.dE)
  _source->_parameters.lambda0 = 0;
  #define lambda0 (_source->_parameters.lambda0)
  _source->_parameters.dlambda = 0;
  #define dlambda (_source->_parameters.dlambda)
  _source->_parameters.flux = 1;
  #define flux (_source->_parameters.flux)
  _source->_parameters.gauss = 0;
  #define gauss (_source->_parameters.gauss)
  _source->_parameters.target_index = + 1;
  #define target_index (_source->_parameters.target_index)

  #define pmul (_source->_parameters.pmul)
  #define square (_source->_parameters.square)
  #define srcArea (_source->_parameters.srcArea)

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
  /* component source=Source_simple() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (0)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _a1->_rotation_absolute, _source->_rotation_absolute);
    rot_transpose(_a1->_rotation_absolute, tr1);
    rot_mul(_source->_rotation_absolute, tr1, _source->_rotation_relative);
    _source->_rotation_is_identity =  rot_test_identity(_source->_rotation_relative);
    tc1 = coords_set(
      0, 0, 0);
    rot_transpose(_a1->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _source->_position_absolute = coords_add(_a1->_position_absolute, tc2);
    tc1 = coords_sub(_a1->_position_absolute, _source->_position_absolute);
    _source->_position_relative = rot_apply(_source->_rotation_absolute, tc1);
  } /* source=Source_simple() AT ROTATED */
  DEBUG_COMPONENT("source", _source->_position_absolute, _source->_rotation_absolute);
  instrument->_position_absolute[2] = _source->_position_absolute;
  instrument->_position_relative[2] = _source->_position_relative;
  instrument->counter_N[2]  = instrument->counter_P[2] = instrument->counter_P2[2] = 0;
  instrument->counter_AbsorbProp[2]= 0;
  return(0);
} /* _source_setpos */

/* component slit1=Slit() SETTING, POSITION/ROTATION */
int _slit1_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_slit1_setpos] component slit1=Slit() SETTING [/usr/share/mcstas/3.0-dev/optics/Slit.comp:47]");
  stracpy(_slit1->_name, "slit1", 16384);
  stracpy(_slit1->_type, "Slit", 16384);
  _slit1->_index=3;
  _slit1->_parameters.xmin = -0.020;
  #define xmin (_slit1->_parameters.xmin)
  _slit1->_parameters.xmax = 0.065;
  #define xmax (_slit1->_parameters.xmax)
  _slit1->_parameters.ymin = -0.075;
  #define ymin (_slit1->_parameters.ymin)
  _slit1->_parameters.ymax = 0.075;
  #define ymax (_slit1->_parameters.ymax)
  _slit1->_parameters.radius = 0;
  #define radius (_slit1->_parameters.radius)
  _slit1->_parameters.xwidth = 0;
  #define xwidth (_slit1->_parameters.xwidth)
  _slit1->_parameters.yheight = 0;
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
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (0)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _a1->_rotation_absolute, _slit1->_rotation_absolute);
    rot_transpose(_source->_rotation_absolute, tr1);
    rot_mul(_slit1->_rotation_absolute, tr1, _slit1->_rotation_relative);
    _slit1->_rotation_is_identity =  rot_test_identity(_slit1->_rotation_relative);
    tc1 = coords_set(
      0, 0, 1.1215);
    rot_transpose(_a1->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _slit1->_position_absolute = coords_add(_a1->_position_absolute, tc2);
    tc1 = coords_sub(_source->_position_absolute, _slit1->_position_absolute);
    _slit1->_position_relative = rot_apply(_slit1->_rotation_absolute, tc1);
  } /* slit1=Slit() AT ROTATED */
  DEBUG_COMPONENT("slit1", _slit1->_position_absolute, _slit1->_rotation_absolute);
  instrument->_position_absolute[3] = _slit1->_position_absolute;
  instrument->_position_relative[3] = _slit1->_position_relative;
  instrument->counter_N[3]  = instrument->counter_P[3] = instrument->counter_P2[3] = 0;
  instrument->counter_AbsorbProp[3]= 0;
  return(0);
} /* _slit1_setpos */

/* component slit2=Slit() SETTING, POSITION/ROTATION */
int _slit2_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_slit2_setpos] component slit2=Slit() SETTING [/usr/share/mcstas/3.0-dev/optics/Slit.comp:47]");
  stracpy(_slit2->_name, "slit2", 16384);
  stracpy(_slit2->_type, "Slit", 16384);
  _slit2->_index=4;
  _slit2->_parameters.xmin = -0.020;
  #define xmin (_slit2->_parameters.xmin)
  _slit2->_parameters.xmax = 0.020;
  #define xmax (_slit2->_parameters.xmax)
  _slit2->_parameters.ymin = -0.040;
  #define ymin (_slit2->_parameters.ymin)
  _slit2->_parameters.ymax = 0.040;
  #define ymax (_slit2->_parameters.ymax)
  _slit2->_parameters.radius = 0;
  #define radius (_slit2->_parameters.radius)
  _slit2->_parameters.xwidth = 0;
  #define xwidth (_slit2->_parameters.xwidth)
  _slit2->_parameters.yheight = 0;
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
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (0)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _a1->_rotation_absolute, _slit2->_rotation_absolute);
    rot_transpose(_slit1->_rotation_absolute, tr1);
    rot_mul(_slit2->_rotation_absolute, tr1, _slit2->_rotation_relative);
    _slit2->_rotation_is_identity =  rot_test_identity(_slit2->_rotation_relative);
    tc1 = coords_set(
      0, 0, 1.900);
    rot_transpose(_a1->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _slit2->_position_absolute = coords_add(_a1->_position_absolute, tc2);
    tc1 = coords_sub(_slit1->_position_absolute, _slit2->_position_absolute);
    _slit2->_position_relative = rot_apply(_slit2->_rotation_absolute, tc1);
  } /* slit2=Slit() AT ROTATED */
  DEBUG_COMPONENT("slit2", _slit2->_position_absolute, _slit2->_rotation_absolute);
  instrument->_position_absolute[4] = _slit2->_position_absolute;
  instrument->_position_relative[4] = _slit2->_position_relative;
  instrument->counter_N[4]  = instrument->counter_P[4] = instrument->counter_P2[4] = 0;
  instrument->counter_AbsorbProp[4]= 0;
  return(0);
} /* _slit2_setpos */

/* component slit3=Slit() SETTING, POSITION/ROTATION */
int _slit3_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_slit3_setpos] component slit3=Slit() SETTING [/usr/share/mcstas/3.0-dev/optics/Slit.comp:47]");
  stracpy(_slit3->_name, "slit3", 16384);
  stracpy(_slit3->_type, "Slit", 16384);
  _slit3->_index=5;
  _slit3->_parameters.xmin = -0.021;
  #define xmin (_slit3->_parameters.xmin)
  _slit3->_parameters.xmax = 0.021;
  #define xmax (_slit3->_parameters.xmax)
  _slit3->_parameters.ymin = -0.041;
  #define ymin (_slit3->_parameters.ymin)
  _slit3->_parameters.ymax = 0.041;
  #define ymax (_slit3->_parameters.ymax)
  _slit3->_parameters.radius = 0;
  #define radius (_slit3->_parameters.radius)
  _slit3->_parameters.xwidth = 0;
  #define xwidth (_slit3->_parameters.xwidth)
  _slit3->_parameters.yheight = 0;
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
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (0)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _a1->_rotation_absolute, _slit3->_rotation_absolute);
    rot_transpose(_slit2->_rotation_absolute, tr1);
    rot_mul(_slit3->_rotation_absolute, tr1, _slit3->_rotation_relative);
    _slit3->_rotation_is_identity =  rot_test_identity(_slit3->_rotation_relative);
    tc1 = coords_set(
      0, 0, 3.288);
    rot_transpose(_a1->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _slit3->_position_absolute = coords_add(_a1->_position_absolute, tc2);
    tc1 = coords_sub(_slit2->_position_absolute, _slit3->_position_absolute);
    _slit3->_position_relative = rot_apply(_slit3->_rotation_absolute, tc1);
  } /* slit3=Slit() AT ROTATED */
  DEBUG_COMPONENT("slit3", _slit3->_position_absolute, _slit3->_rotation_absolute);
  instrument->_position_absolute[5] = _slit3->_position_absolute;
  instrument->_position_relative[5] = _slit3->_position_relative;
  instrument->counter_N[5]  = instrument->counter_P[5] = instrument->counter_P2[5] = 0;
  instrument->counter_AbsorbProp[5]= 0;
  return(0);
} /* _slit3_setpos */

/* component focus_mono=Arm() SETTING, POSITION/ROTATION */
int _focus_mono_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_focus_mono_setpos] component focus_mono=Arm() SETTING [Arm:0]");
  stracpy(_focus_mono->_name, "focus_mono", 16384);
  stracpy(_focus_mono->_type, "Arm", 16384);
  _focus_mono->_index=6;
  /* component focus_mono=Arm() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (instrument->_parameters._PHM)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _a1->_rotation_absolute, _focus_mono->_rotation_absolute);
    rot_transpose(_slit3->_rotation_absolute, tr1);
    rot_mul(_focus_mono->_rotation_absolute, tr1, _focus_mono->_rotation_relative);
    _focus_mono->_rotation_is_identity =  rot_test_identity(_focus_mono->_rotation_relative);
    tc1 = coords_set(
      0, 0, 3.56);
    rot_transpose(_a1->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _focus_mono->_position_absolute = coords_add(_a1->_position_absolute, tc2);
    tc1 = coords_sub(_slit3->_position_absolute, _focus_mono->_position_absolute);
    _focus_mono->_position_relative = rot_apply(_focus_mono->_rotation_absolute, tc1);
  } /* focus_mono=Arm() AT ROTATED */
  DEBUG_COMPONENT("focus_mono", _focus_mono->_position_absolute, _focus_mono->_rotation_absolute);
  instrument->_position_absolute[6] = _focus_mono->_position_absolute;
  instrument->_position_relative[6] = _focus_mono->_position_relative;
  instrument->counter_N[6]  = instrument->counter_P[6] = instrument->counter_P2[6] = 0;
  instrument->counter_AbsorbProp[6]= 0;
  return(0);
} /* _focus_mono_setpos */

/* component m0=Monochromator_flat() SETTING, POSITION/ROTATION */
int _m0_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_m0_setpos] component m0=Monochromator_flat() SETTING [/usr/share/mcstas/3.0-dev/optics/Monochromator_flat.comp:102]");
  stracpy(_m0->_name, "m0", 16384);
  stracpy(_m0->_type, "Monochromator_flat", 16384);
  _m0->_index=7;
  _m0->_parameters.zmin = -0.0375;
  #define zmin (_m0->_parameters.zmin)
  _m0->_parameters.zmax = 0.0375;
  #define zmax (_m0->_parameters.zmax)
  _m0->_parameters.ymin = -0.006;
  #define ymin (_m0->_parameters.ymin)
  _m0->_parameters.ymax = 0.006;
  #define ymax (_m0->_parameters.ymax)
  _m0->_parameters.zwidth = 0;
  #define zwidth (_m0->_parameters.zwidth)
  _m0->_parameters.yheight = 0;
  #define yheight (_m0->_parameters.yheight)
  _m0->_parameters.mosaich = tas1_mono_mosaic;
  #define mosaich (_m0->_parameters.mosaich)
  _m0->_parameters.mosaicv = tas1_mono_mosaic;
  #define mosaicv (_m0->_parameters.mosaicv)
  _m0->_parameters.r0 = tas1_mono_r0;
  #define r0 (_m0->_parameters.r0)
  _m0->_parameters.Q = tas1_mono_q;
  #define Q (_m0->_parameters.Q)
  _m0->_parameters.DM = 0;
  #define DM (_m0->_parameters.DM)

  #define mos_rms_y (_m0->_parameters.mos_rms_y)
  #define mos_rms_z (_m0->_parameters.mos_rms_z)
  #define mos_rms_max (_m0->_parameters.mos_rms_max)
  #define mono_Q (_m0->_parameters.mono_Q)

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
  /* component m0=Monochromator_flat() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (0)*DEG2RAD, (mrot0)*DEG2RAD);
    rot_mul(tr1, _focus_mono->_rotation_absolute, _m0->_rotation_absolute);
    rot_transpose(_focus_mono->_rotation_absolute, tr1);
    rot_mul(_m0->_rotation_absolute, tr1, _m0->_rotation_relative);
    _m0->_rotation_is_identity =  rot_test_identity(_m0->_rotation_relative);
    tc1 = coords_set(
      0, mpos0, 0);
    rot_transpose(_focus_mono->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _m0->_position_absolute = coords_add(_focus_mono->_position_absolute, tc2);
    tc1 = coords_sub(_focus_mono->_position_absolute, _m0->_position_absolute);
    _m0->_position_relative = rot_apply(_m0->_rotation_absolute, tc1);
  } /* m0=Monochromator_flat() AT ROTATED */
  DEBUG_COMPONENT("m0", _m0->_position_absolute, _m0->_rotation_absolute);
  instrument->_position_absolute[7] = _m0->_position_absolute;
  instrument->_position_relative[7] = _m0->_position_relative;
  instrument->counter_N[7]  = instrument->counter_P[7] = instrument->counter_P2[7] = 0;
  instrument->counter_AbsorbProp[7]= 0;
  return(0);
} /* _m0_setpos */

/* component m1=Monochromator_flat() SETTING, POSITION/ROTATION */
int _m1_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_m1_setpos] component m1=Monochromator_flat() SETTING [/usr/share/mcstas/3.0-dev/optics/Monochromator_flat.comp:102]");
  stracpy(_m1->_name, "m1", 16384);
  stracpy(_m1->_type, "Monochromator_flat", 16384);
  _m1->_index=8;
  _m1->_parameters.zmin = -0.0375;
  #define zmin (_m1->_parameters.zmin)
  _m1->_parameters.zmax = 0.0375;
  #define zmax (_m1->_parameters.zmax)
  _m1->_parameters.ymin = -0.006;
  #define ymin (_m1->_parameters.ymin)
  _m1->_parameters.ymax = 0.006;
  #define ymax (_m1->_parameters.ymax)
  _m1->_parameters.zwidth = 0;
  #define zwidth (_m1->_parameters.zwidth)
  _m1->_parameters.yheight = 0;
  #define yheight (_m1->_parameters.yheight)
  _m1->_parameters.mosaich = tas1_mono_mosaic;
  #define mosaich (_m1->_parameters.mosaich)
  _m1->_parameters.mosaicv = tas1_mono_mosaic;
  #define mosaicv (_m1->_parameters.mosaicv)
  _m1->_parameters.r0 = tas1_mono_r0;
  #define r0 (_m1->_parameters.r0)
  _m1->_parameters.Q = tas1_mono_q;
  #define Q (_m1->_parameters.Q)
  _m1->_parameters.DM = 0;
  #define DM (_m1->_parameters.DM)

  #define mos_rms_y (_m1->_parameters.mos_rms_y)
  #define mos_rms_z (_m1->_parameters.mos_rms_z)
  #define mos_rms_max (_m1->_parameters.mos_rms_max)
  #define mono_Q (_m1->_parameters.mono_Q)

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
  /* component m1=Monochromator_flat() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (0)*DEG2RAD, (mrot1)*DEG2RAD);
    rot_mul(tr1, _focus_mono->_rotation_absolute, _m1->_rotation_absolute);
    rot_transpose(_m0->_rotation_absolute, tr1);
    rot_mul(_m1->_rotation_absolute, tr1, _m1->_rotation_relative);
    _m1->_rotation_is_identity =  rot_test_identity(_m1->_rotation_relative);
    tc1 = coords_set(
      0, mpos1, 0);
    rot_transpose(_focus_mono->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _m1->_position_absolute = coords_add(_focus_mono->_position_absolute, tc2);
    tc1 = coords_sub(_m0->_position_absolute, _m1->_position_absolute);
    _m1->_position_relative = rot_apply(_m1->_rotation_absolute, tc1);
  } /* m1=Monochromator_flat() AT ROTATED */
  DEBUG_COMPONENT("m1", _m1->_position_absolute, _m1->_rotation_absolute);
  instrument->_position_absolute[8] = _m1->_position_absolute;
  instrument->_position_relative[8] = _m1->_position_relative;
  instrument->counter_N[8]  = instrument->counter_P[8] = instrument->counter_P2[8] = 0;
  instrument->counter_AbsorbProp[8]= 0;
  return(0);
} /* _m1_setpos */

/* component m2=Monochromator_flat() SETTING, POSITION/ROTATION */
int _m2_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_m2_setpos] component m2=Monochromator_flat() SETTING [/usr/share/mcstas/3.0-dev/optics/Monochromator_flat.comp:102]");
  stracpy(_m2->_name, "m2", 16384);
  stracpy(_m2->_type, "Monochromator_flat", 16384);
  _m2->_index=9;
  _m2->_parameters.zmin = -0.0375;
  #define zmin (_m2->_parameters.zmin)
  _m2->_parameters.zmax = 0.0375;
  #define zmax (_m2->_parameters.zmax)
  _m2->_parameters.ymin = -0.006;
  #define ymin (_m2->_parameters.ymin)
  _m2->_parameters.ymax = 0.006;
  #define ymax (_m2->_parameters.ymax)
  _m2->_parameters.zwidth = 0;
  #define zwidth (_m2->_parameters.zwidth)
  _m2->_parameters.yheight = 0;
  #define yheight (_m2->_parameters.yheight)
  _m2->_parameters.mosaich = tas1_mono_mosaic;
  #define mosaich (_m2->_parameters.mosaich)
  _m2->_parameters.mosaicv = tas1_mono_mosaic;
  #define mosaicv (_m2->_parameters.mosaicv)
  _m2->_parameters.r0 = tas1_mono_r0;
  #define r0 (_m2->_parameters.r0)
  _m2->_parameters.Q = tas1_mono_q;
  #define Q (_m2->_parameters.Q)
  _m2->_parameters.DM = 0;
  #define DM (_m2->_parameters.DM)

  #define mos_rms_y (_m2->_parameters.mos_rms_y)
  #define mos_rms_z (_m2->_parameters.mos_rms_z)
  #define mos_rms_max (_m2->_parameters.mos_rms_max)
  #define mono_Q (_m2->_parameters.mono_Q)

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
  /* component m2=Monochromator_flat() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (0)*DEG2RAD, (mrot2)*DEG2RAD);
    rot_mul(tr1, _focus_mono->_rotation_absolute, _m2->_rotation_absolute);
    rot_transpose(_m1->_rotation_absolute, tr1);
    rot_mul(_m2->_rotation_absolute, tr1, _m2->_rotation_relative);
    _m2->_rotation_is_identity =  rot_test_identity(_m2->_rotation_relative);
    tc1 = coords_set(
      0, mpos2, 0);
    rot_transpose(_focus_mono->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _m2->_position_absolute = coords_add(_focus_mono->_position_absolute, tc2);
    tc1 = coords_sub(_m1->_position_absolute, _m2->_position_absolute);
    _m2->_position_relative = rot_apply(_m2->_rotation_absolute, tc1);
  } /* m2=Monochromator_flat() AT ROTATED */
  DEBUG_COMPONENT("m2", _m2->_position_absolute, _m2->_rotation_absolute);
  instrument->_position_absolute[9] = _m2->_position_absolute;
  instrument->_position_relative[9] = _m2->_position_relative;
  instrument->counter_N[9]  = instrument->counter_P[9] = instrument->counter_P2[9] = 0;
  instrument->counter_AbsorbProp[9]= 0;
  return(0);
} /* _m2_setpos */

/* component m3=Monochromator_flat() SETTING, POSITION/ROTATION */
int _m3_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_m3_setpos] component m3=Monochromator_flat() SETTING [/usr/share/mcstas/3.0-dev/optics/Monochromator_flat.comp:102]");
  stracpy(_m3->_name, "m3", 16384);
  stracpy(_m3->_type, "Monochromator_flat", 16384);
  _m3->_index=10;
  _m3->_parameters.zmin = -0.0375;
  #define zmin (_m3->_parameters.zmin)
  _m3->_parameters.zmax = 0.0375;
  #define zmax (_m3->_parameters.zmax)
  _m3->_parameters.ymin = -0.006;
  #define ymin (_m3->_parameters.ymin)
  _m3->_parameters.ymax = 0.006;
  #define ymax (_m3->_parameters.ymax)
  _m3->_parameters.zwidth = 0;
  #define zwidth (_m3->_parameters.zwidth)
  _m3->_parameters.yheight = 0;
  #define yheight (_m3->_parameters.yheight)
  _m3->_parameters.mosaich = tas1_mono_mosaic;
  #define mosaich (_m3->_parameters.mosaich)
  _m3->_parameters.mosaicv = tas1_mono_mosaic;
  #define mosaicv (_m3->_parameters.mosaicv)
  _m3->_parameters.r0 = tas1_mono_r0;
  #define r0 (_m3->_parameters.r0)
  _m3->_parameters.Q = tas1_mono_q;
  #define Q (_m3->_parameters.Q)
  _m3->_parameters.DM = 0;
  #define DM (_m3->_parameters.DM)

  #define mos_rms_y (_m3->_parameters.mos_rms_y)
  #define mos_rms_z (_m3->_parameters.mos_rms_z)
  #define mos_rms_max (_m3->_parameters.mos_rms_max)
  #define mono_Q (_m3->_parameters.mono_Q)

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
  /* component m3=Monochromator_flat() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (0)*DEG2RAD, (mrot3)*DEG2RAD);
    rot_mul(tr1, _focus_mono->_rotation_absolute, _m3->_rotation_absolute);
    rot_transpose(_m2->_rotation_absolute, tr1);
    rot_mul(_m3->_rotation_absolute, tr1, _m3->_rotation_relative);
    _m3->_rotation_is_identity =  rot_test_identity(_m3->_rotation_relative);
    tc1 = coords_set(
      0, mpos3, 0);
    rot_transpose(_focus_mono->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _m3->_position_absolute = coords_add(_focus_mono->_position_absolute, tc2);
    tc1 = coords_sub(_m2->_position_absolute, _m3->_position_absolute);
    _m3->_position_relative = rot_apply(_m3->_rotation_absolute, tc1);
  } /* m3=Monochromator_flat() AT ROTATED */
  DEBUG_COMPONENT("m3", _m3->_position_absolute, _m3->_rotation_absolute);
  instrument->_position_absolute[10] = _m3->_position_absolute;
  instrument->_position_relative[10] = _m3->_position_relative;
  instrument->counter_N[10]  = instrument->counter_P[10] = instrument->counter_P2[10] = 0;
  instrument->counter_AbsorbProp[10]= 0;
  return(0);
} /* _m3_setpos */

/* component m4=Monochromator_flat() SETTING, POSITION/ROTATION */
int _m4_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_m4_setpos] component m4=Monochromator_flat() SETTING [/usr/share/mcstas/3.0-dev/optics/Monochromator_flat.comp:102]");
  stracpy(_m4->_name, "m4", 16384);
  stracpy(_m4->_type, "Monochromator_flat", 16384);
  _m4->_index=11;
  _m4->_parameters.zmin = -0.0375;
  #define zmin (_m4->_parameters.zmin)
  _m4->_parameters.zmax = 0.0375;
  #define zmax (_m4->_parameters.zmax)
  _m4->_parameters.ymin = -0.006;
  #define ymin (_m4->_parameters.ymin)
  _m4->_parameters.ymax = 0.006;
  #define ymax (_m4->_parameters.ymax)
  _m4->_parameters.zwidth = 0;
  #define zwidth (_m4->_parameters.zwidth)
  _m4->_parameters.yheight = 0;
  #define yheight (_m4->_parameters.yheight)
  _m4->_parameters.mosaich = tas1_mono_mosaic;
  #define mosaich (_m4->_parameters.mosaich)
  _m4->_parameters.mosaicv = tas1_mono_mosaic;
  #define mosaicv (_m4->_parameters.mosaicv)
  _m4->_parameters.r0 = tas1_mono_r0;
  #define r0 (_m4->_parameters.r0)
  _m4->_parameters.Q = tas1_mono_q;
  #define Q (_m4->_parameters.Q)
  _m4->_parameters.DM = 0;
  #define DM (_m4->_parameters.DM)

  #define mos_rms_y (_m4->_parameters.mos_rms_y)
  #define mos_rms_z (_m4->_parameters.mos_rms_z)
  #define mos_rms_max (_m4->_parameters.mos_rms_max)
  #define mono_Q (_m4->_parameters.mono_Q)

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
  /* component m4=Monochromator_flat() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (0)*DEG2RAD, (mrot4)*DEG2RAD);
    rot_mul(tr1, _focus_mono->_rotation_absolute, _m4->_rotation_absolute);
    rot_transpose(_m3->_rotation_absolute, tr1);
    rot_mul(_m4->_rotation_absolute, tr1, _m4->_rotation_relative);
    _m4->_rotation_is_identity =  rot_test_identity(_m4->_rotation_relative);
    tc1 = coords_set(
      0, mpos4, 0);
    rot_transpose(_focus_mono->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _m4->_position_absolute = coords_add(_focus_mono->_position_absolute, tc2);
    tc1 = coords_sub(_m3->_position_absolute, _m4->_position_absolute);
    _m4->_position_relative = rot_apply(_m4->_rotation_absolute, tc1);
  } /* m4=Monochromator_flat() AT ROTATED */
  DEBUG_COMPONENT("m4", _m4->_position_absolute, _m4->_rotation_absolute);
  instrument->_position_absolute[11] = _m4->_position_absolute;
  instrument->_position_relative[11] = _m4->_position_relative;
  instrument->counter_N[11]  = instrument->counter_P[11] = instrument->counter_P2[11] = 0;
  instrument->counter_AbsorbProp[11]= 0;
  return(0);
} /* _m4_setpos */

/* component m5=Monochromator_flat() SETTING, POSITION/ROTATION */
int _m5_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_m5_setpos] component m5=Monochromator_flat() SETTING [/usr/share/mcstas/3.0-dev/optics/Monochromator_flat.comp:102]");
  stracpy(_m5->_name, "m5", 16384);
  stracpy(_m5->_type, "Monochromator_flat", 16384);
  _m5->_index=12;
  _m5->_parameters.zmin = -0.0375;
  #define zmin (_m5->_parameters.zmin)
  _m5->_parameters.zmax = 0.0375;
  #define zmax (_m5->_parameters.zmax)
  _m5->_parameters.ymin = -0.006;
  #define ymin (_m5->_parameters.ymin)
  _m5->_parameters.ymax = 0.006;
  #define ymax (_m5->_parameters.ymax)
  _m5->_parameters.zwidth = 0;
  #define zwidth (_m5->_parameters.zwidth)
  _m5->_parameters.yheight = 0;
  #define yheight (_m5->_parameters.yheight)
  _m5->_parameters.mosaich = tas1_mono_mosaic;
  #define mosaich (_m5->_parameters.mosaich)
  _m5->_parameters.mosaicv = tas1_mono_mosaic;
  #define mosaicv (_m5->_parameters.mosaicv)
  _m5->_parameters.r0 = tas1_mono_r0;
  #define r0 (_m5->_parameters.r0)
  _m5->_parameters.Q = tas1_mono_q;
  #define Q (_m5->_parameters.Q)
  _m5->_parameters.DM = 0;
  #define DM (_m5->_parameters.DM)

  #define mos_rms_y (_m5->_parameters.mos_rms_y)
  #define mos_rms_z (_m5->_parameters.mos_rms_z)
  #define mos_rms_max (_m5->_parameters.mos_rms_max)
  #define mono_Q (_m5->_parameters.mono_Q)

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
  /* component m5=Monochromator_flat() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (0)*DEG2RAD, (mrot5)*DEG2RAD);
    rot_mul(tr1, _focus_mono->_rotation_absolute, _m5->_rotation_absolute);
    rot_transpose(_m4->_rotation_absolute, tr1);
    rot_mul(_m5->_rotation_absolute, tr1, _m5->_rotation_relative);
    _m5->_rotation_is_identity =  rot_test_identity(_m5->_rotation_relative);
    tc1 = coords_set(
      0, mpos5, 0);
    rot_transpose(_focus_mono->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _m5->_position_absolute = coords_add(_focus_mono->_position_absolute, tc2);
    tc1 = coords_sub(_m4->_position_absolute, _m5->_position_absolute);
    _m5->_position_relative = rot_apply(_m5->_rotation_absolute, tc1);
  } /* m5=Monochromator_flat() AT ROTATED */
  DEBUG_COMPONENT("m5", _m5->_position_absolute, _m5->_rotation_absolute);
  instrument->_position_absolute[12] = _m5->_position_absolute;
  instrument->_position_relative[12] = _m5->_position_relative;
  instrument->counter_N[12]  = instrument->counter_P[12] = instrument->counter_P2[12] = 0;
  instrument->counter_AbsorbProp[12]= 0;
  return(0);
} /* _m5_setpos */

/* component m6=Monochromator_flat() SETTING, POSITION/ROTATION */
int _m6_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_m6_setpos] component m6=Monochromator_flat() SETTING [/usr/share/mcstas/3.0-dev/optics/Monochromator_flat.comp:102]");
  stracpy(_m6->_name, "m6", 16384);
  stracpy(_m6->_type, "Monochromator_flat", 16384);
  _m6->_index=13;
  _m6->_parameters.zmin = -0.0375;
  #define zmin (_m6->_parameters.zmin)
  _m6->_parameters.zmax = 0.0375;
  #define zmax (_m6->_parameters.zmax)
  _m6->_parameters.ymin = -0.006;
  #define ymin (_m6->_parameters.ymin)
  _m6->_parameters.ymax = 0.006;
  #define ymax (_m6->_parameters.ymax)
  _m6->_parameters.zwidth = 0;
  #define zwidth (_m6->_parameters.zwidth)
  _m6->_parameters.yheight = 0;
  #define yheight (_m6->_parameters.yheight)
  _m6->_parameters.mosaich = tas1_mono_mosaic;
  #define mosaich (_m6->_parameters.mosaich)
  _m6->_parameters.mosaicv = tas1_mono_mosaic;
  #define mosaicv (_m6->_parameters.mosaicv)
  _m6->_parameters.r0 = tas1_mono_r0;
  #define r0 (_m6->_parameters.r0)
  _m6->_parameters.Q = tas1_mono_q;
  #define Q (_m6->_parameters.Q)
  _m6->_parameters.DM = 0;
  #define DM (_m6->_parameters.DM)

  #define mos_rms_y (_m6->_parameters.mos_rms_y)
  #define mos_rms_z (_m6->_parameters.mos_rms_z)
  #define mos_rms_max (_m6->_parameters.mos_rms_max)
  #define mono_Q (_m6->_parameters.mono_Q)

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
  /* component m6=Monochromator_flat() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (0)*DEG2RAD, (mrot6)*DEG2RAD);
    rot_mul(tr1, _focus_mono->_rotation_absolute, _m6->_rotation_absolute);
    rot_transpose(_m5->_rotation_absolute, tr1);
    rot_mul(_m6->_rotation_absolute, tr1, _m6->_rotation_relative);
    _m6->_rotation_is_identity =  rot_test_identity(_m6->_rotation_relative);
    tc1 = coords_set(
      0, mpos6, 0);
    rot_transpose(_focus_mono->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _m6->_position_absolute = coords_add(_focus_mono->_position_absolute, tc2);
    tc1 = coords_sub(_m5->_position_absolute, _m6->_position_absolute);
    _m6->_position_relative = rot_apply(_m6->_rotation_absolute, tc1);
  } /* m6=Monochromator_flat() AT ROTATED */
  DEBUG_COMPONENT("m6", _m6->_position_absolute, _m6->_rotation_absolute);
  instrument->_position_absolute[13] = _m6->_position_absolute;
  instrument->_position_relative[13] = _m6->_position_relative;
  instrument->counter_N[13]  = instrument->counter_P[13] = instrument->counter_P2[13] = 0;
  instrument->counter_AbsorbProp[13]= 0;
  return(0);
} /* _m6_setpos */

/* component m7=Monochromator_flat() SETTING, POSITION/ROTATION */
int _m7_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_m7_setpos] component m7=Monochromator_flat() SETTING [/usr/share/mcstas/3.0-dev/optics/Monochromator_flat.comp:102]");
  stracpy(_m7->_name, "m7", 16384);
  stracpy(_m7->_type, "Monochromator_flat", 16384);
  _m7->_index=14;
  _m7->_parameters.zmin = -0.0375;
  #define zmin (_m7->_parameters.zmin)
  _m7->_parameters.zmax = 0.0375;
  #define zmax (_m7->_parameters.zmax)
  _m7->_parameters.ymin = -0.006;
  #define ymin (_m7->_parameters.ymin)
  _m7->_parameters.ymax = 0.006;
  #define ymax (_m7->_parameters.ymax)
  _m7->_parameters.zwidth = 0;
  #define zwidth (_m7->_parameters.zwidth)
  _m7->_parameters.yheight = 0;
  #define yheight (_m7->_parameters.yheight)
  _m7->_parameters.mosaich = tas1_mono_mosaic;
  #define mosaich (_m7->_parameters.mosaich)
  _m7->_parameters.mosaicv = tas1_mono_mosaic;
  #define mosaicv (_m7->_parameters.mosaicv)
  _m7->_parameters.r0 = tas1_mono_r0;
  #define r0 (_m7->_parameters.r0)
  _m7->_parameters.Q = tas1_mono_q;
  #define Q (_m7->_parameters.Q)
  _m7->_parameters.DM = 0;
  #define DM (_m7->_parameters.DM)

  #define mos_rms_y (_m7->_parameters.mos_rms_y)
  #define mos_rms_z (_m7->_parameters.mos_rms_z)
  #define mos_rms_max (_m7->_parameters.mos_rms_max)
  #define mono_Q (_m7->_parameters.mono_Q)

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
  /* component m7=Monochromator_flat() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (0)*DEG2RAD, (mrot7)*DEG2RAD);
    rot_mul(tr1, _focus_mono->_rotation_absolute, _m7->_rotation_absolute);
    rot_transpose(_m6->_rotation_absolute, tr1);
    rot_mul(_m7->_rotation_absolute, tr1, _m7->_rotation_relative);
    _m7->_rotation_is_identity =  rot_test_identity(_m7->_rotation_relative);
    tc1 = coords_set(
      0, mpos7, 0);
    rot_transpose(_focus_mono->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _m7->_position_absolute = coords_add(_focus_mono->_position_absolute, tc2);
    tc1 = coords_sub(_m6->_position_absolute, _m7->_position_absolute);
    _m7->_position_relative = rot_apply(_m7->_rotation_absolute, tc1);
  } /* m7=Monochromator_flat() AT ROTATED */
  DEBUG_COMPONENT("m7", _m7->_position_absolute, _m7->_rotation_absolute);
  instrument->_position_absolute[14] = _m7->_position_absolute;
  instrument->_position_relative[14] = _m7->_position_relative;
  instrument->counter_N[14]  = instrument->counter_P[14] = instrument->counter_P2[14] = 0;
  instrument->counter_AbsorbProp[14]= 0;
  return(0);
} /* _m7_setpos */

/* component a2=Arm() SETTING, POSITION/ROTATION */
int _a2_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_a2_setpos] component a2=Arm() SETTING [Arm:0]");
  stracpy(_a2->_name, "a2", 16384);
  stracpy(_a2->_type, "Arm", 16384);
  _a2->_index=15;
  /* component a2=Arm() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (instrument->_parameters._TTM)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _a1->_rotation_absolute, _a2->_rotation_absolute);
    rot_transpose(_m7->_rotation_absolute, tr1);
    rot_mul(_a2->_rotation_absolute, tr1, _a2->_rotation_relative);
    _a2->_rotation_is_identity =  rot_test_identity(_a2->_rotation_relative);
    tc1 = coords_set(
      0, 0, 0);
    rot_transpose(_focus_mono->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _a2->_position_absolute = coords_add(_focus_mono->_position_absolute, tc2);
    tc1 = coords_sub(_m7->_position_absolute, _a2->_position_absolute);
    _a2->_position_relative = rot_apply(_a2->_rotation_absolute, tc1);
  } /* a2=Arm() AT ROTATED */
  DEBUG_COMPONENT("a2", _a2->_position_absolute, _a2->_rotation_absolute);
  instrument->_position_absolute[15] = _a2->_position_absolute;
  instrument->_position_relative[15] = _a2->_position_relative;
  instrument->counter_N[15]  = instrument->counter_P[15] = instrument->counter_P2[15] = 0;
  instrument->counter_AbsorbProp[15]= 0;
  return(0);
} /* _a2_setpos */

/* component slitMS1=Slit() SETTING, POSITION/ROTATION */
int _slitMS1_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_slitMS1_setpos] component slitMS1=Slit() SETTING [/usr/share/mcstas/3.0-dev/optics/Slit.comp:47]");
  stracpy(_slitMS1->_name, "slitMS1", 16384);
  stracpy(_slitMS1->_type, "Slit", 16384);
  _slitMS1->_index=16;
  _slitMS1->_parameters.xmin = -0.0105;
  #define xmin (_slitMS1->_parameters.xmin)
  _slitMS1->_parameters.xmax = 0.0105;
  #define xmax (_slitMS1->_parameters.xmax)
  _slitMS1->_parameters.ymin = -0.035;
  #define ymin (_slitMS1->_parameters.ymin)
  _slitMS1->_parameters.ymax = 0.035;
  #define ymax (_slitMS1->_parameters.ymax)
  _slitMS1->_parameters.radius = 0;
  #define radius (_slitMS1->_parameters.radius)
  _slitMS1->_parameters.xwidth = 0;
  #define xwidth (_slitMS1->_parameters.xwidth)
  _slitMS1->_parameters.yheight = 0;
  #define yheight (_slitMS1->_parameters.yheight)

  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef radius
  #undef xwidth
  #undef yheight
  /* component slitMS1=Slit() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (0)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _a2->_rotation_absolute, _slitMS1->_rotation_absolute);
    rot_transpose(_a2->_rotation_absolute, tr1);
    rot_mul(_slitMS1->_rotation_absolute, tr1, _slitMS1->_rotation_relative);
    _slitMS1->_rotation_is_identity =  rot_test_identity(_slitMS1->_rotation_relative);
    tc1 = coords_set(
      0, 0, 0.565);
    rot_transpose(_a2->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _slitMS1->_position_absolute = coords_add(_a2->_position_absolute, tc2);
    tc1 = coords_sub(_a2->_position_absolute, _slitMS1->_position_absolute);
    _slitMS1->_position_relative = rot_apply(_slitMS1->_rotation_absolute, tc1);
  } /* slitMS1=Slit() AT ROTATED */
  DEBUG_COMPONENT("slitMS1", _slitMS1->_position_absolute, _slitMS1->_rotation_absolute);
  instrument->_position_absolute[16] = _slitMS1->_position_absolute;
  instrument->_position_relative[16] = _slitMS1->_position_relative;
  instrument->counter_N[16]  = instrument->counter_P[16] = instrument->counter_P2[16] = 0;
  instrument->counter_AbsorbProp[16]= 0;
  return(0);
} /* _slitMS1_setpos */

/* component slitMS2=Slit() SETTING, POSITION/ROTATION */
int _slitMS2_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_slitMS2_setpos] component slitMS2=Slit() SETTING [/usr/share/mcstas/3.0-dev/optics/Slit.comp:47]");
  stracpy(_slitMS2->_name, "slitMS2", 16384);
  stracpy(_slitMS2->_type, "Slit", 16384);
  _slitMS2->_index=17;
  _slitMS2->_parameters.xmin = -0.0105;
  #define xmin (_slitMS2->_parameters.xmin)
  _slitMS2->_parameters.xmax = 0.0105;
  #define xmax (_slitMS2->_parameters.xmax)
  _slitMS2->_parameters.ymin = -0.035;
  #define ymin (_slitMS2->_parameters.ymin)
  _slitMS2->_parameters.ymax = 0.035;
  #define ymax (_slitMS2->_parameters.ymax)
  _slitMS2->_parameters.radius = 0;
  #define radius (_slitMS2->_parameters.radius)
  _slitMS2->_parameters.xwidth = 0;
  #define xwidth (_slitMS2->_parameters.xwidth)
  _slitMS2->_parameters.yheight = 0;
  #define yheight (_slitMS2->_parameters.yheight)

  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef radius
  #undef xwidth
  #undef yheight
  /* component slitMS2=Slit() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (0)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _a2->_rotation_absolute, _slitMS2->_rotation_absolute);
    rot_transpose(_slitMS1->_rotation_absolute, tr1);
    rot_mul(_slitMS2->_rotation_absolute, tr1, _slitMS2->_rotation_relative);
    _slitMS2->_rotation_is_identity =  rot_test_identity(_slitMS2->_rotation_relative);
    tc1 = coords_set(
      0, 0, 0.855);
    rot_transpose(_a2->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _slitMS2->_position_absolute = coords_add(_a2->_position_absolute, tc2);
    tc1 = coords_sub(_slitMS1->_position_absolute, _slitMS2->_position_absolute);
    _slitMS2->_position_relative = rot_apply(_slitMS2->_rotation_absolute, tc1);
  } /* slitMS2=Slit() AT ROTATED */
  DEBUG_COMPONENT("slitMS2", _slitMS2->_position_absolute, _slitMS2->_rotation_absolute);
  instrument->_position_absolute[17] = _slitMS2->_position_absolute;
  instrument->_position_relative[17] = _slitMS2->_position_relative;
  instrument->counter_N[17]  = instrument->counter_P[17] = instrument->counter_P2[17] = 0;
  instrument->counter_AbsorbProp[17]= 0;
  return(0);
} /* _slitMS2_setpos */

/* component c1=Collimator_linear() SETTING, POSITION/ROTATION */
int _c1_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_c1_setpos] component c1=Collimator_linear() SETTING [/usr/share/mcstas/3.0-dev/optics/Collimator_linear.comp:54]");
  stracpy(_c1->_name, "c1", 16384);
  stracpy(_c1->_type, "Collimator_linear", 16384);
  _c1->_index=18;
  _c1->_parameters.xmin = -0.02;
  #define xmin (_c1->_parameters.xmin)
  _c1->_parameters.xmax = 0.02;
  #define xmax (_c1->_parameters.xmax)
  _c1->_parameters.ymin = -0.0375;
  #define ymin (_c1->_parameters.ymin)
  _c1->_parameters.ymax = 0.0375;
  #define ymax (_c1->_parameters.ymax)
  _c1->_parameters.xwidth = 0;
  #define xwidth (_c1->_parameters.xwidth)
  _c1->_parameters.yheight = 0;
  #define yheight (_c1->_parameters.yheight)
  _c1->_parameters.length = 0.250;
  #define length (_c1->_parameters.length)
  _c1->_parameters.divergence = instrument->_parameters._C1;
  #define divergence (_c1->_parameters.divergence)
  _c1->_parameters.transmission = 1;
  #define transmission (_c1->_parameters.transmission)
  _c1->_parameters.divergenceV = 0;
  #define divergenceV (_c1->_parameters.divergenceV)

  #define slope (_c1->_parameters.slope)
  #define slopeV (_c1->_parameters.slopeV)

  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef xwidth
  #undef yheight
  #undef length
  #undef divergence
  #undef transmission
  #undef divergenceV
  #undef slope
  #undef slopeV
  /* component c1=Collimator_linear() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (OMC1_d)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _a2->_rotation_absolute, _c1->_rotation_absolute);
    rot_transpose(_slitMS2->_rotation_absolute, tr1);
    rot_mul(_c1->_rotation_absolute, tr1, _c1->_rotation_relative);
    _c1->_rotation_is_identity =  rot_test_identity(_c1->_rotation_relative);
    tc1 = coords_set(
      0, 0, 0.87);
    rot_transpose(_a2->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _c1->_position_absolute = coords_add(_a2->_position_absolute, tc2);
    tc1 = coords_sub(_slitMS2->_position_absolute, _c1->_position_absolute);
    _c1->_position_relative = rot_apply(_c1->_rotation_absolute, tc1);
  } /* c1=Collimator_linear() AT ROTATED */
  DEBUG_COMPONENT("c1", _c1->_position_absolute, _c1->_rotation_absolute);
  instrument->_position_absolute[18] = _c1->_position_absolute;
  instrument->_position_relative[18] = _c1->_position_relative;
  instrument->counter_N[18]  = instrument->counter_P[18] = instrument->counter_P2[18] = 0;
  instrument->counter_AbsorbProp[18]= 0;
  return(0);
} /* _c1_setpos */

/* component slitMS3=Slit() SETTING, POSITION/ROTATION */
int _slitMS3_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_slitMS3_setpos] component slitMS3=Slit() SETTING [/usr/share/mcstas/3.0-dev/optics/Slit.comp:47]");
  stracpy(_slitMS3->_name, "slitMS3", 16384);
  stracpy(_slitMS3->_type, "Slit", 16384);
  _slitMS3->_index=19;
  _slitMS3->_parameters.xmin = -0.01;
  #define xmin (_slitMS3->_parameters.xmin)
  _slitMS3->_parameters.xmax = 0.01;
  #define xmax (_slitMS3->_parameters.xmax)
  _slitMS3->_parameters.ymin = -0.01;
  #define ymin (_slitMS3->_parameters.ymin)
  _slitMS3->_parameters.ymax = 0.01;
  #define ymax (_slitMS3->_parameters.ymax)
  _slitMS3->_parameters.radius = 0.025;
  #define radius (_slitMS3->_parameters.radius)
  _slitMS3->_parameters.xwidth = 0;
  #define xwidth (_slitMS3->_parameters.xwidth)
  _slitMS3->_parameters.yheight = 0;
  #define yheight (_slitMS3->_parameters.yheight)

  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef radius
  #undef xwidth
  #undef yheight
  /* component slitMS3=Slit() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (0)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _a2->_rotation_absolute, _slitMS3->_rotation_absolute);
    rot_transpose(_c1->_rotation_absolute, tr1);
    rot_mul(_slitMS3->_rotation_absolute, tr1, _slitMS3->_rotation_relative);
    _slitMS3->_rotation_is_identity =  rot_test_identity(_slitMS3->_rotation_relative);
    tc1 = coords_set(
      0, 0, 1.130);
    rot_transpose(_a2->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _slitMS3->_position_absolute = coords_add(_a2->_position_absolute, tc2);
    tc1 = coords_sub(_c1->_position_absolute, _slitMS3->_position_absolute);
    _slitMS3->_position_relative = rot_apply(_slitMS3->_rotation_absolute, tc1);
  } /* slitMS3=Slit() AT ROTATED */
  DEBUG_COMPONENT("slitMS3", _slitMS3->_position_absolute, _slitMS3->_rotation_absolute);
  instrument->_position_absolute[19] = _slitMS3->_position_absolute;
  instrument->_position_relative[19] = _slitMS3->_position_relative;
  instrument->counter_N[19]  = instrument->counter_P[19] = instrument->counter_P2[19] = 0;
  instrument->counter_AbsorbProp[19]= 0;
  return(0);
} /* _slitMS3_setpos */

/* component slitMS4=Slit() SETTING, POSITION/ROTATION */
int _slitMS4_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_slitMS4_setpos] component slitMS4=Slit() SETTING [/usr/share/mcstas/3.0-dev/optics/Slit.comp:47]");
  stracpy(_slitMS4->_name, "slitMS4", 16384);
  stracpy(_slitMS4->_type, "Slit", 16384);
  _slitMS4->_index=20;
  _slitMS4->_parameters.xmin = -0.01;
  #define xmin (_slitMS4->_parameters.xmin)
  _slitMS4->_parameters.xmax = 0.01;
  #define xmax (_slitMS4->_parameters.xmax)
  _slitMS4->_parameters.ymin = -0.01;
  #define ymin (_slitMS4->_parameters.ymin)
  _slitMS4->_parameters.ymax = 0.01;
  #define ymax (_slitMS4->_parameters.ymax)
  _slitMS4->_parameters.radius = 0.025;
  #define radius (_slitMS4->_parameters.radius)
  _slitMS4->_parameters.xwidth = 0;
  #define xwidth (_slitMS4->_parameters.xwidth)
  _slitMS4->_parameters.yheight = 0;
  #define yheight (_slitMS4->_parameters.yheight)

  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef radius
  #undef xwidth
  #undef yheight
  /* component slitMS4=Slit() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (0)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _a2->_rotation_absolute, _slitMS4->_rotation_absolute);
    rot_transpose(_slitMS3->_rotation_absolute, tr1);
    rot_mul(_slitMS4->_rotation_absolute, tr1, _slitMS4->_rotation_relative);
    _slitMS4->_rotation_is_identity =  rot_test_identity(_slitMS4->_rotation_relative);
    tc1 = coords_set(
      0, 0, 1.180);
    rot_transpose(_a2->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _slitMS4->_position_absolute = coords_add(_a2->_position_absolute, tc2);
    tc1 = coords_sub(_slitMS3->_position_absolute, _slitMS4->_position_absolute);
    _slitMS4->_position_relative = rot_apply(_slitMS4->_rotation_absolute, tc1);
  } /* slitMS4=Slit() AT ROTATED */
  DEBUG_COMPONENT("slitMS4", _slitMS4->_position_absolute, _slitMS4->_rotation_absolute);
  instrument->_position_absolute[20] = _slitMS4->_position_absolute;
  instrument->_position_relative[20] = _slitMS4->_position_relative;
  instrument->counter_N[20]  = instrument->counter_P[20] = instrument->counter_P2[20] = 0;
  instrument->counter_AbsorbProp[20]= 0;
  return(0);
} /* _slitMS4_setpos */

/* component slitMS5=Slit() SETTING, POSITION/ROTATION */
int _slitMS5_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_slitMS5_setpos] component slitMS5=Slit() SETTING [/usr/share/mcstas/3.0-dev/optics/Slit.comp:47]");
  stracpy(_slitMS5->_name, "slitMS5", 16384);
  stracpy(_slitMS5->_type, "Slit", 16384);
  _slitMS5->_index=21;
  _slitMS5->_parameters.xmin = -0.01;
  #define xmin (_slitMS5->_parameters.xmin)
  _slitMS5->_parameters.xmax = 0.01;
  #define xmax (_slitMS5->_parameters.xmax)
  _slitMS5->_parameters.ymin = -0.01;
  #define ymin (_slitMS5->_parameters.ymin)
  _slitMS5->_parameters.ymax = 0.01;
  #define ymax (_slitMS5->_parameters.ymax)
  _slitMS5->_parameters.radius = 0.0275;
  #define radius (_slitMS5->_parameters.radius)
  _slitMS5->_parameters.xwidth = 0;
  #define xwidth (_slitMS5->_parameters.xwidth)
  _slitMS5->_parameters.yheight = 0;
  #define yheight (_slitMS5->_parameters.yheight)

  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef radius
  #undef xwidth
  #undef yheight
  /* component slitMS5=Slit() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (0)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _a2->_rotation_absolute, _slitMS5->_rotation_absolute);
    rot_transpose(_slitMS4->_rotation_absolute, tr1);
    rot_mul(_slitMS5->_rotation_absolute, tr1, _slitMS5->_rotation_relative);
    _slitMS5->_rotation_is_identity =  rot_test_identity(_slitMS5->_rotation_relative);
    tc1 = coords_set(
      0, 0, 1.230);
    rot_transpose(_a2->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _slitMS5->_position_absolute = coords_add(_a2->_position_absolute, tc2);
    tc1 = coords_sub(_slitMS4->_position_absolute, _slitMS5->_position_absolute);
    _slitMS5->_position_relative = rot_apply(_slitMS5->_rotation_absolute, tc1);
  } /* slitMS5=Slit() AT ROTATED */
  DEBUG_COMPONENT("slitMS5", _slitMS5->_position_absolute, _slitMS5->_rotation_absolute);
  instrument->_position_absolute[21] = _slitMS5->_position_absolute;
  instrument->_position_relative[21] = _slitMS5->_position_relative;
  instrument->counter_N[21]  = instrument->counter_P[21] = instrument->counter_P2[21] = 0;
  instrument->counter_AbsorbProp[21]= 0;
  return(0);
} /* _slitMS5_setpos */

/* component mon=Monitor() SETTING, POSITION/ROTATION */
int _mon_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_mon_setpos] component mon=Monitor() SETTING [/usr/share/mcstas/3.0-dev/monitors/Monitor.comp:57]");
  stracpy(_mon->_name, "mon", 16384);
  stracpy(_mon->_type, "Monitor", 16384);
  _mon->_index=22;
  _mon->_parameters.xmin = -0.025;
  #define xmin (_mon->_parameters.xmin)
  _mon->_parameters.xmax = 0.025;
  #define xmax (_mon->_parameters.xmax)
  _mon->_parameters.ymin = -0.0375;
  #define ymin (_mon->_parameters.ymin)
  _mon->_parameters.ymax = 0.0375;
  #define ymax (_mon->_parameters.ymax)
  _mon->_parameters.xwidth = 0;
  #define xwidth (_mon->_parameters.xwidth)
  _mon->_parameters.yheight = 0;
  #define yheight (_mon->_parameters.yheight)
  _mon->_parameters.restore_neutron = 0;
  #define restore_neutron (_mon->_parameters.restore_neutron)

  #define Nsum (_mon->_parameters.Nsum)
  #define psum (_mon->_parameters.psum)
  #define p2sum (_mon->_parameters.p2sum)

  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef xwidth
  #undef yheight
  #undef restore_neutron
  #undef Nsum
  #undef psum
  #undef p2sum
  /* component mon=Monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (0)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _a2->_rotation_absolute, _mon->_rotation_absolute);
    rot_transpose(_slitMS5->_rotation_absolute, tr1);
    rot_mul(_mon->_rotation_absolute, tr1, _mon->_rotation_relative);
    _mon->_rotation_is_identity =  rot_test_identity(_mon->_rotation_relative);
    tc1 = coords_set(
      0, 0, 1.280);
    rot_transpose(_a2->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _mon->_position_absolute = coords_add(_a2->_position_absolute, tc2);
    tc1 = coords_sub(_slitMS5->_position_absolute, _mon->_position_absolute);
    _mon->_position_relative = rot_apply(_mon->_rotation_absolute, tc1);
  } /* mon=Monitor() AT ROTATED */
  DEBUG_COMPONENT("mon", _mon->_position_absolute, _mon->_rotation_absolute);
  instrument->_position_absolute[22] = _mon->_position_absolute;
  instrument->_position_relative[22] = _mon->_position_relative;
  instrument->counter_N[22]  = instrument->counter_P[22] = instrument->counter_P2[22] = 0;
  instrument->counter_AbsorbProp[22]= 0;
  return(0);
} /* _mon_setpos */

/* component emon1=E_monitor() SETTING, POSITION/ROTATION */
int _emon1_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_emon1_setpos] component emon1=E_monitor() SETTING [/usr/share/mcstas/3.0-dev/monitors/E_monitor.comp:66]");
  stracpy(_emon1->_name, "emon1", 16384);
  stracpy(_emon1->_type, "E_monitor", 16384);
  _emon1->_index=23;
  _emon1->_parameters.nE = 35;
  #define nE (_emon1->_parameters.nE)
  if("linup_3_1.vmon" && strlen("linup_3_1.vmon"))
    stracpy(_emon1->_parameters.filename, "linup_3_1.vmon" ? "linup_3_1.vmon" : "", 16384);
  else 
  _emon1->_parameters.filename[0]='\0';
  #define filename (_emon1->_parameters.filename)
  _emon1->_parameters.xmin = -0.01;
  #define xmin (_emon1->_parameters.xmin)
  _emon1->_parameters.xmax = 0.01;
  #define xmax (_emon1->_parameters.xmax)
  _emon1->_parameters.ymin = -0.1;
  #define ymin (_emon1->_parameters.ymin)
  _emon1->_parameters.ymax = 0.1;
  #define ymax (_emon1->_parameters.ymax)
  _emon1->_parameters.xwidth = 0;
  #define xwidth (_emon1->_parameters.xwidth)
  _emon1->_parameters.yheight = 0;
  #define yheight (_emon1->_parameters.yheight)
  _emon1->_parameters.Emin = 19.25;
  #define Emin (_emon1->_parameters.Emin)
  _emon1->_parameters.Emax = 20.75;
  #define Emax (_emon1->_parameters.Emax)
  _emon1->_parameters.restore_neutron = 0;
  #define restore_neutron (_emon1->_parameters.restore_neutron)

  #define E_N (_emon1->_parameters.E_N)
  #define E_p (_emon1->_parameters.E_p)
  #define E_p2 (_emon1->_parameters.E_p2)
  #define S_p (_emon1->_parameters.S_p)
  #define S_pE (_emon1->_parameters.S_pE)
  #define S_pE2 (_emon1->_parameters.S_pE2)

  #undef nE
  #undef filename
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef xwidth
  #undef yheight
  #undef Emin
  #undef Emax
  #undef restore_neutron
  #undef E_N
  #undef E_p
  #undef E_p2
  #undef S_p
  #undef S_pE
  #undef S_pE2
  /* component emon1=E_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (0)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _a2->_rotation_absolute, _emon1->_rotation_absolute);
    rot_transpose(_mon->_rotation_absolute, tr1);
    rot_mul(_emon1->_rotation_absolute, tr1, _emon1->_rotation_relative);
    _emon1->_rotation_is_identity =  rot_test_identity(_emon1->_rotation_relative);
    tc1 = coords_set(
      0, 0, 1.5);
    rot_transpose(_a2->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _emon1->_position_absolute = coords_add(_a2->_position_absolute, tc2);
    tc1 = coords_sub(_mon->_position_absolute, _emon1->_position_absolute);
    _emon1->_position_relative = rot_apply(_emon1->_rotation_absolute, tc1);
  } /* emon1=E_monitor() AT ROTATED */
  DEBUG_COMPONENT("emon1", _emon1->_position_absolute, _emon1->_rotation_absolute);
  instrument->_position_absolute[23] = _emon1->_position_absolute;
  instrument->_position_relative[23] = _emon1->_position_relative;
  instrument->counter_N[23]  = instrument->counter_P[23] = instrument->counter_P2[23] = 0;
  instrument->counter_AbsorbProp[23]= 0;
  return(0);
} /* _emon1_setpos */

/* component sample=Slit() SETTING, POSITION/ROTATION */
int _sample_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_sample_setpos] component sample=Slit() SETTING [/usr/share/mcstas/3.0-dev/optics/Slit.comp:47]");
  stracpy(_sample->_name, "sample", 16384);
  stracpy(_sample->_type, "Slit", 16384);
  _sample->_index=24;
  _sample->_parameters.xmin = -0.00525;
  #define xmin (_sample->_parameters.xmin)
  _sample->_parameters.xmax = 0.00525;
  #define xmax (_sample->_parameters.xmax)
  _sample->_parameters.ymin = -0.02025;
  #define ymin (_sample->_parameters.ymin)
  _sample->_parameters.ymax = 0.02025;
  #define ymax (_sample->_parameters.ymax)
  _sample->_parameters.radius = 0;
  #define radius (_sample->_parameters.radius)
  _sample->_parameters.xwidth = 0;
  #define xwidth (_sample->_parameters.xwidth)
  _sample->_parameters.yheight = 0;
  #define yheight (_sample->_parameters.yheight)

  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef radius
  #undef xwidth
  #undef yheight
  /* component sample=Slit() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (0)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _a2->_rotation_absolute, _sample->_rotation_absolute);
    rot_transpose(_emon1->_rotation_absolute, tr1);
    rot_mul(_sample->_rotation_absolute, tr1, _sample->_rotation_relative);
    _sample->_rotation_is_identity =  rot_test_identity(_sample->_rotation_relative);
    tc1 = coords_set(
      0, 0, 1.565);
    rot_transpose(_a2->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _sample->_position_absolute = coords_add(_a2->_position_absolute, tc2);
    tc1 = coords_sub(_emon1->_position_absolute, _sample->_position_absolute);
    _sample->_position_relative = rot_apply(_sample->_rotation_absolute, tc1);
  } /* sample=Slit() AT ROTATED */
  DEBUG_COMPONENT("sample", _sample->_position_absolute, _sample->_rotation_absolute);
  instrument->_position_absolute[24] = _sample->_position_absolute;
  instrument->_position_relative[24] = _sample->_position_relative;
  instrument->counter_N[24]  = instrument->counter_P[24] = instrument->counter_P2[24] = 0;
  instrument->counter_AbsorbProp[24]= 0;
  return(0);
} /* _sample_setpos */

/* component slit1mm=Slit() SETTING, POSITION/ROTATION */
int _slit1mm_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_slit1mm_setpos] component slit1mm=Slit() SETTING [/usr/share/mcstas/3.0-dev/optics/Slit.comp:47]");
  stracpy(_slit1mm->_name, "slit1mm", 16384);
  stracpy(_slit1mm->_type, "Slit", 16384);
  _slit1mm->_index=25;
  _slit1mm->_parameters.xmin = -0.0005;
  #define xmin (_slit1mm->_parameters.xmin)
  _slit1mm->_parameters.xmax = 0.0005;
  #define xmax (_slit1mm->_parameters.xmax)
  _slit1mm->_parameters.ymin = -0.02;
  #define ymin (_slit1mm->_parameters.ymin)
  _slit1mm->_parameters.ymax = 0.02;
  #define ymax (_slit1mm->_parameters.ymax)
  _slit1mm->_parameters.radius = 0;
  #define radius (_slit1mm->_parameters.radius)
  _slit1mm->_parameters.xwidth = 0;
  #define xwidth (_slit1mm->_parameters.xwidth)
  _slit1mm->_parameters.yheight = 0;
  #define yheight (_slit1mm->_parameters.yheight)

  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef radius
  #undef xwidth
  #undef yheight
  /* component slit1mm=Slit() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (0)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _sample->_rotation_absolute, _slit1mm->_rotation_absolute);
    rot_transpose(_sample->_rotation_absolute, tr1);
    rot_mul(_slit1mm->_rotation_absolute, tr1, _slit1mm->_rotation_relative);
    _slit1mm->_rotation_is_identity =  rot_test_identity(_slit1mm->_rotation_relative);
    tc1 = coords_set(
      0, 0, 0.030);
    rot_transpose(_sample->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _slit1mm->_position_absolute = coords_add(_sample->_position_absolute, tc2);
    tc1 = coords_sub(_sample->_position_absolute, _slit1mm->_position_absolute);
    _slit1mm->_position_relative = rot_apply(_slit1mm->_rotation_absolute, tc1);
  } /* slit1mm=Slit() AT ROTATED */
  DEBUG_COMPONENT("slit1mm", _slit1mm->_position_absolute, _slit1mm->_rotation_absolute);
  instrument->_position_absolute[25] = _slit1mm->_position_absolute;
  instrument->_position_relative[25] = _slit1mm->_position_relative;
  instrument->counter_N[25]  = instrument->counter_P[25] = instrument->counter_P2[25] = 0;
  instrument->counter_AbsorbProp[25]= 0;
  return(0);
} /* _slit1mm_setpos */

/* component a3=Arm() SETTING, POSITION/ROTATION */
int _a3_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_a3_setpos] component a3=Arm() SETTING [Arm:0]");
  stracpy(_a3->_name, "a3", 16384);
  stracpy(_a3->_type, "Arm", 16384);
  _a3->_index=26;
  /* component a3=Arm() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (instrument->_parameters._TT)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _a2->_rotation_absolute, _a3->_rotation_absolute);
    rot_transpose(_slit1mm->_rotation_absolute, tr1);
    rot_mul(_a3->_rotation_absolute, tr1, _a3->_rotation_relative);
    _a3->_rotation_is_identity =  rot_test_identity(_a3->_rotation_relative);
    tc1 = coords_set(
      0, 0, 0);
    rot_transpose(_sample->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _a3->_position_absolute = coords_add(_sample->_position_absolute, tc2);
    tc1 = coords_sub(_slit1mm->_position_absolute, _a3->_position_absolute);
    _a3->_position_relative = rot_apply(_a3->_rotation_absolute, tc1);
  } /* a3=Arm() AT ROTATED */
  DEBUG_COMPONENT("a3", _a3->_position_absolute, _a3->_rotation_absolute);
  instrument->_position_absolute[26] = _a3->_position_absolute;
  instrument->_position_relative[26] = _a3->_position_relative;
  instrument->counter_N[26]  = instrument->counter_P[26] = instrument->counter_P2[26] = 0;
  instrument->counter_AbsorbProp[26]= 0;
  return(0);
} /* _a3_setpos */

/* component c2=Collimator_linear() SETTING, POSITION/ROTATION */
int _c2_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_c2_setpos] component c2=Collimator_linear() SETTING [/usr/share/mcstas/3.0-dev/optics/Collimator_linear.comp:54]");
  stracpy(_c2->_name, "c2", 16384);
  stracpy(_c2->_type, "Collimator_linear", 16384);
  _c2->_index=27;
  _c2->_parameters.xmin = -0.02;
  #define xmin (_c2->_parameters.xmin)
  _c2->_parameters.xmax = 0.02;
  #define xmax (_c2->_parameters.xmax)
  _c2->_parameters.ymin = -0.0315;
  #define ymin (_c2->_parameters.ymin)
  _c2->_parameters.ymax = 0.0315;
  #define ymax (_c2->_parameters.ymax)
  _c2->_parameters.xwidth = 0;
  #define xwidth (_c2->_parameters.xwidth)
  _c2->_parameters.yheight = 0;
  #define yheight (_c2->_parameters.yheight)
  _c2->_parameters.length = 0.300;
  #define length (_c2->_parameters.length)
  _c2->_parameters.divergence = instrument->_parameters._C2;
  #define divergence (_c2->_parameters.divergence)
  _c2->_parameters.transmission = 1;
  #define transmission (_c2->_parameters.transmission)
  _c2->_parameters.divergenceV = 0;
  #define divergenceV (_c2->_parameters.divergenceV)

  #define slope (_c2->_parameters.slope)
  #define slopeV (_c2->_parameters.slopeV)

  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef xwidth
  #undef yheight
  #undef length
  #undef divergence
  #undef transmission
  #undef divergenceV
  #undef slope
  #undef slopeV
  /* component c2=Collimator_linear() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (0)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _a3->_rotation_absolute, _c2->_rotation_absolute);
    rot_transpose(_a3->_rotation_absolute, tr1);
    rot_mul(_c2->_rotation_absolute, tr1, _c2->_rotation_relative);
    _c2->_rotation_is_identity =  rot_test_identity(_c2->_rotation_relative);
    tc1 = coords_set(
      0, 0, 0.370);
    rot_transpose(_a3->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _c2->_position_absolute = coords_add(_a3->_position_absolute, tc2);
    tc1 = coords_sub(_a3->_position_absolute, _c2->_position_absolute);
    _c2->_position_relative = rot_apply(_c2->_rotation_absolute, tc1);
  } /* c2=Collimator_linear() AT ROTATED */
  DEBUG_COMPONENT("c2", _c2->_position_absolute, _c2->_rotation_absolute);
  instrument->_position_absolute[27] = _c2->_position_absolute;
  instrument->_position_relative[27] = _c2->_position_relative;
  instrument->counter_N[27]  = instrument->counter_P[27] = instrument->counter_P2[27] = 0;
  instrument->counter_AbsorbProp[27]= 0;
  return(0);
} /* _c2_setpos */

/* component ana=Arm() SETTING, POSITION/ROTATION */
int _ana_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_ana_setpos] component ana=Arm() SETTING [Arm:0]");
  stracpy(_ana->_name, "ana", 16384);
  stracpy(_ana->_type, "Arm", 16384);
  _ana->_index=28;
  /* component ana=Arm() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (0)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _a3->_rotation_absolute, _ana->_rotation_absolute);
    rot_transpose(_c2->_rotation_absolute, tr1);
    rot_mul(_ana->_rotation_absolute, tr1, _ana->_rotation_relative);
    _ana->_rotation_is_identity =  rot_test_identity(_ana->_rotation_relative);
    tc1 = coords_set(
      0, 0, 0.770);
    rot_transpose(_a3->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _ana->_position_absolute = coords_add(_a3->_position_absolute, tc2);
    tc1 = coords_sub(_c2->_position_absolute, _ana->_position_absolute);
    _ana->_position_relative = rot_apply(_ana->_rotation_absolute, tc1);
  } /* ana=Arm() AT ROTATED */
  DEBUG_COMPONENT("ana", _ana->_position_absolute, _ana->_rotation_absolute);
  instrument->_position_absolute[28] = _ana->_position_absolute;
  instrument->_position_relative[28] = _ana->_position_relative;
  instrument->counter_N[28]  = instrument->counter_P[28] = instrument->counter_P2[28] = 0;
  instrument->counter_AbsorbProp[28]= 0;
  return(0);
} /* _ana_setpos */

/* component a4=Arm() SETTING, POSITION/ROTATION */
int _a4_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_a4_setpos] component a4=Arm() SETTING [Arm:0]");
  stracpy(_a4->_name, "a4", 16384);
  stracpy(_a4->_type, "Arm", 16384);
  _a4->_index=29;
  /* component a4=Arm() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (0)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _a3->_rotation_absolute, _a4->_rotation_absolute);
    rot_transpose(_ana->_rotation_absolute, tr1);
    rot_mul(_a4->_rotation_absolute, tr1, _a4->_rotation_relative);
    _a4->_rotation_is_identity =  rot_test_identity(_a4->_rotation_relative);
    tc1 = coords_set(
      0, 0, 0);
    rot_transpose(_ana->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _a4->_position_absolute = coords_add(_ana->_position_absolute, tc2);
    tc1 = coords_sub(_ana->_position_absolute, _a4->_position_absolute);
    _a4->_position_relative = rot_apply(_a4->_rotation_absolute, tc1);
  } /* a4=Arm() AT ROTATED */
  DEBUG_COMPONENT("a4", _a4->_position_absolute, _a4->_rotation_absolute);
  instrument->_position_absolute[29] = _a4->_position_absolute;
  instrument->_position_relative[29] = _a4->_position_relative;
  instrument->counter_N[29]  = instrument->counter_P[29] = instrument->counter_P2[29] = 0;
  instrument->counter_AbsorbProp[29]= 0;
  return(0);
} /* _a4_setpos */

/* component c3=Collimator_linear() SETTING, POSITION/ROTATION */
int _c3_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_c3_setpos] component c3=Collimator_linear() SETTING [/usr/share/mcstas/3.0-dev/optics/Collimator_linear.comp:54]");
  stracpy(_c3->_name, "c3", 16384);
  stracpy(_c3->_type, "Collimator_linear", 16384);
  _c3->_index=30;
  _c3->_parameters.xmin = -0.02;
  #define xmin (_c3->_parameters.xmin)
  _c3->_parameters.xmax = 0.02;
  #define xmax (_c3->_parameters.xmax)
  _c3->_parameters.ymin = -0.05;
  #define ymin (_c3->_parameters.ymin)
  _c3->_parameters.ymax = 0.05;
  #define ymax (_c3->_parameters.ymax)
  _c3->_parameters.xwidth = 0;
  #define xwidth (_c3->_parameters.xwidth)
  _c3->_parameters.yheight = 0;
  #define yheight (_c3->_parameters.yheight)
  _c3->_parameters.length = 0.270;
  #define length (_c3->_parameters.length)
  _c3->_parameters.divergence = instrument->_parameters._C3;
  #define divergence (_c3->_parameters.divergence)
  _c3->_parameters.transmission = 1;
  #define transmission (_c3->_parameters.transmission)
  _c3->_parameters.divergenceV = 0;
  #define divergenceV (_c3->_parameters.divergenceV)

  #define slope (_c3->_parameters.slope)
  #define slopeV (_c3->_parameters.slopeV)

  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef xwidth
  #undef yheight
  #undef length
  #undef divergence
  #undef transmission
  #undef divergenceV
  #undef slope
  #undef slopeV
  /* component c3=Collimator_linear() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (0)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _a4->_rotation_absolute, _c3->_rotation_absolute);
    rot_transpose(_a4->_rotation_absolute, tr1);
    rot_mul(_c3->_rotation_absolute, tr1, _c3->_rotation_relative);
    _c3->_rotation_is_identity =  rot_test_identity(_c3->_rotation_relative);
    tc1 = coords_set(
      0, 0, 0.104);
    rot_transpose(_a4->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _c3->_position_absolute = coords_add(_a4->_position_absolute, tc2);
    tc1 = coords_sub(_a4->_position_absolute, _c3->_position_absolute);
    _c3->_position_relative = rot_apply(_c3->_rotation_absolute, tc1);
  } /* c3=Collimator_linear() AT ROTATED */
  DEBUG_COMPONENT("c3", _c3->_position_absolute, _c3->_rotation_absolute);
  instrument->_position_absolute[30] = _c3->_position_absolute;
  instrument->_position_relative[30] = _c3->_position_relative;
  instrument->counter_N[30]  = instrument->counter_P[30] = instrument->counter_P2[30] = 0;
  instrument->counter_AbsorbProp[30]= 0;
  return(0);
} /* _c3_setpos */

/* component sng=Monitor() SETTING, POSITION/ROTATION */
int _sng_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_sng_setpos] component sng=Monitor() SETTING [/usr/share/mcstas/3.0-dev/monitors/Monitor.comp:57]");
  stracpy(_sng->_name, "sng", 16384);
  stracpy(_sng->_type, "Monitor", 16384);
  _sng->_index=31;
  _sng->_parameters.xmin = -0.01;
  #define xmin (_sng->_parameters.xmin)
  _sng->_parameters.xmax = 0.01;
  #define xmax (_sng->_parameters.xmax)
  _sng->_parameters.ymin = -0.045;
  #define ymin (_sng->_parameters.ymin)
  _sng->_parameters.ymax = 0.045;
  #define ymax (_sng->_parameters.ymax)
  _sng->_parameters.xwidth = 0;
  #define xwidth (_sng->_parameters.xwidth)
  _sng->_parameters.yheight = 0;
  #define yheight (_sng->_parameters.yheight)
  _sng->_parameters.restore_neutron = 0;
  #define restore_neutron (_sng->_parameters.restore_neutron)

  #define Nsum (_sng->_parameters.Nsum)
  #define psum (_sng->_parameters.psum)
  #define p2sum (_sng->_parameters.p2sum)

  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef xwidth
  #undef yheight
  #undef restore_neutron
  #undef Nsum
  #undef psum
  #undef p2sum
  /* component sng=Monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (0)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _a4->_rotation_absolute, _sng->_rotation_absolute);
    rot_transpose(_c3->_rotation_absolute, tr1);
    rot_mul(_sng->_rotation_absolute, tr1, _sng->_rotation_relative);
    _sng->_rotation_is_identity =  rot_test_identity(_sng->_rotation_relative);
    tc1 = coords_set(
      0, 0, 0.43);
    rot_transpose(_a4->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _sng->_position_absolute = coords_add(_a4->_position_absolute, tc2);
    tc1 = coords_sub(_c3->_position_absolute, _sng->_position_absolute);
    _sng->_position_relative = rot_apply(_sng->_rotation_absolute, tc1);
  } /* sng=Monitor() AT ROTATED */
  DEBUG_COMPONENT("sng", _sng->_position_absolute, _sng->_rotation_absolute);
  instrument->_position_absolute[31] = _sng->_position_absolute;
  instrument->_position_relative[31] = _sng->_position_relative;
  instrument->counter_N[31]  = instrument->counter_P[31] = instrument->counter_P2[31] = 0;
  instrument->counter_AbsorbProp[31]= 0;
  return(0);
} /* _sng_setpos */

/* component Emon2=E_monitor() SETTING, POSITION/ROTATION */
int _Emon2_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Emon2_setpos] component Emon2=E_monitor() SETTING [/usr/share/mcstas/3.0-dev/monitors/E_monitor.comp:66]");
  stracpy(_Emon2->_name, "Emon2", 16384);
  stracpy(_Emon2->_type, "E_monitor", 16384);
  _Emon2->_index=32;
  _Emon2->_parameters.nE = 35;
  #define nE (_Emon2->_parameters.nE)
  if("linup_3_2.vmon" && strlen("linup_3_2.vmon"))
    stracpy(_Emon2->_parameters.filename, "linup_3_2.vmon" ? "linup_3_2.vmon" : "", 16384);
  else 
  _Emon2->_parameters.filename[0]='\0';
  #define filename (_Emon2->_parameters.filename)
  _Emon2->_parameters.xmin = -0.0125;
  #define xmin (_Emon2->_parameters.xmin)
  _Emon2->_parameters.xmax = 0.0125;
  #define xmax (_Emon2->_parameters.xmax)
  _Emon2->_parameters.ymin = -0.05;
  #define ymin (_Emon2->_parameters.ymin)
  _Emon2->_parameters.ymax = 0.05;
  #define ymax (_Emon2->_parameters.ymax)
  _Emon2->_parameters.xwidth = 0;
  #define xwidth (_Emon2->_parameters.xwidth)
  _Emon2->_parameters.yheight = 0;
  #define yheight (_Emon2->_parameters.yheight)
  _Emon2->_parameters.Emin = 19.25;
  #define Emin (_Emon2->_parameters.Emin)
  _Emon2->_parameters.Emax = 20.75;
  #define Emax (_Emon2->_parameters.Emax)
  _Emon2->_parameters.restore_neutron = 0;
  #define restore_neutron (_Emon2->_parameters.restore_neutron)

  #define E_N (_Emon2->_parameters.E_N)
  #define E_p (_Emon2->_parameters.E_p)
  #define E_p2 (_Emon2->_parameters.E_p2)
  #define S_p (_Emon2->_parameters.S_p)
  #define S_pE (_Emon2->_parameters.S_pE)
  #define S_pE2 (_Emon2->_parameters.S_pE2)

  #undef nE
  #undef filename
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef xwidth
  #undef yheight
  #undef Emin
  #undef Emax
  #undef restore_neutron
  #undef E_N
  #undef E_p
  #undef E_p2
  #undef S_p
  #undef S_pE
  #undef S_pE2
  /* component Emon2=E_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (0)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _a4->_rotation_absolute, _Emon2->_rotation_absolute);
    rot_transpose(_sng->_rotation_absolute, tr1);
    rot_mul(_Emon2->_rotation_absolute, tr1, _Emon2->_rotation_relative);
    _Emon2->_rotation_is_identity =  rot_test_identity(_Emon2->_rotation_relative);
    tc1 = coords_set(
      0, 0, 0.430001);
    rot_transpose(_a4->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Emon2->_position_absolute = coords_add(_a4->_position_absolute, tc2);
    tc1 = coords_sub(_sng->_position_absolute, _Emon2->_position_absolute);
    _Emon2->_position_relative = rot_apply(_Emon2->_rotation_absolute, tc1);
  } /* Emon2=E_monitor() AT ROTATED */
  DEBUG_COMPONENT("Emon2", _Emon2->_position_absolute, _Emon2->_rotation_absolute);
  instrument->_position_absolute[32] = _Emon2->_position_absolute;
  instrument->_position_relative[32] = _Emon2->_position_relative;
  instrument->counter_N[32]  = instrument->counter_P[32] = instrument->counter_P2[32] = 0;
  instrument->counter_AbsorbProp[32]= 0;
  return(0);
} /* _Emon2_setpos */

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

_class_Collimator_linear *class_Collimator_linear_init(_class_Collimator_linear *_comp
) {
  #define xmin (_comp->_parameters.xmin)
  #define xmax (_comp->_parameters.xmax)
  #define ymin (_comp->_parameters.ymin)
  #define ymax (_comp->_parameters.ymax)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define length (_comp->_parameters.length)
  #define divergence (_comp->_parameters.divergence)
  #define transmission (_comp->_parameters.transmission)
  #define divergenceV (_comp->_parameters.divergenceV)
  #define slope (_comp->_parameters.slope)
  #define slopeV (_comp->_parameters.slopeV)
slope = tan(MIN2RAD*divergence);
  slopeV= tan(MIN2RAD*divergenceV);
  if (xwidth  > 0) { xmax = xwidth/2;  xmin = -xmax; }
  if (yheight > 0) { ymax = yheight/2; ymin = -ymax; }

  if ((xmin >= xmax) || (ymin >= ymax)) {
    printf("Collimator_linear: %s: Null slit opening area !\n"
	         "ERROR              (xwidth,yheight,xmin,xmax,ymin,ymax). Exiting",
           NAME_CURRENT_COMP);
    exit(0);
  }

  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef xwidth
  #undef yheight
  #undef length
  #undef divergence
  #undef transmission
  #undef divergenceV
  #undef slope
  #undef slopeV
  return(_comp);
} /* class_Collimator_linear_init */

_class_Monitor *class_Monitor_init(_class_Monitor *_comp
) {
  #define xmin (_comp->_parameters.xmin)
  #define xmax (_comp->_parameters.xmax)
  #define ymin (_comp->_parameters.ymin)
  #define ymax (_comp->_parameters.ymax)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define restore_neutron (_comp->_parameters.restore_neutron)
  #define Nsum (_comp->_parameters.Nsum)
  #define psum (_comp->_parameters.psum)
  #define p2sum (_comp->_parameters.p2sum)
psum = 0;
p2sum = 0;
Nsum = 0;

if (xwidth  > 0) { xmax = xwidth/2;  xmin = -xmax; }
    if (yheight > 0) { ymax = yheight/2; ymin = -ymax; }

    if ((xmin >= xmax) || (ymin >= ymax)) {
            printf("Monitor: %s: Null detection area !\n"
                   "ERROR    (xwidth,yheight,xmin,xmax,ymin,ymax). Exiting",
           NAME_CURRENT_COMP);
      exit(0);
    }
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef xwidth
  #undef yheight
  #undef restore_neutron
  #undef Nsum
  #undef psum
  #undef p2sum
  return(_comp);
} /* class_Monitor_init */

_class_E_monitor *class_E_monitor_init(_class_E_monitor *_comp
) {
  #define nE (_comp->_parameters.nE)
  #define filename (_comp->_parameters.filename)
  #define xmin (_comp->_parameters.xmin)
  #define xmax (_comp->_parameters.xmax)
  #define ymin (_comp->_parameters.ymin)
  #define ymax (_comp->_parameters.ymax)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define Emin (_comp->_parameters.Emin)
  #define Emax (_comp->_parameters.Emax)
  #define restore_neutron (_comp->_parameters.restore_neutron)
  #define E_N (_comp->_parameters.E_N)
  #define E_p (_comp->_parameters.E_p)
  #define E_p2 (_comp->_parameters.E_p2)
  #define S_p (_comp->_parameters.S_p)
  #define S_pE (_comp->_parameters.S_pE)
  #define S_pE2 (_comp->_parameters.S_pE2)
  int i;

  if (xwidth  > 0) { xmax = xwidth/2;  xmin = -xmax; }
  if (yheight > 0) { ymax = yheight/2; ymin = -ymax; }

  if ((xmin >= xmax) || (ymin >= ymax)) {
    printf("E_monitor: %s: Null detection area !\n"
           "ERROR      (xwidth,yheight,xmin,xmax,ymin,ymax). Exiting",
           NAME_CURRENT_COMP);
    exit(0);
  }

  E_N = create_darr1d(nE);
  E_p = create_darr1d(nE);
  E_p2 = create_darr1d(nE);

  for (i=0; i<nE; i++)
  {
    E_N[i] = 0;
    E_p[i] = 0;
    E_p2[i] = 0;
  }
  S_p = S_pE = S_pE2 = 0;
  #undef nE
  #undef filename
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef xwidth
  #undef yheight
  #undef Emin
  #undef Emax
  #undef restore_neutron
  #undef E_N
  #undef E_p
  #undef E_p2
  #undef S_p
  #undef S_pE
  #undef S_pE2
  return(_comp);
} /* class_E_monitor_init */



int init(void) { /* called by mccode_main for TAS1_Diff_Slit:INITIALISE */
  DEBUG_INSTR();

  /* code_main/parseoptions/readparams sets instrument parameters value */
  stracpy(instrument->_name, "TAS1_Diff_Slit", 256);

  /* Instrument 'TAS1_Diff_Slit' INITIALISE */
  SIG_MESSAGE("[TAS1_Diff_Slit] INITIALISE [linup-3.instr:62]");
  #define PHM (instrument->_parameters._PHM)
  #define TTM (instrument->_parameters._TTM)
  #define TT (instrument->_parameters._TT)
  #define C1 (instrument->_parameters._C1)
  #define OMC1 (instrument->_parameters._OMC1)
  #define C2 (instrument->_parameters._C2)
  #define C3 (instrument->_parameters._C3)
{
  double d = 0.0125;    /* 12.5 mm between slab centers. */
  double phi = 0.5443;    /* Rotation between adjacent slabs. */
  mpos0 = -3.5*d; mrot0 = -3.5*phi;
  mpos1 = -2.5*d; mrot1 = -2.5*phi;
  mpos2 = -1.5*d; mrot2 = -1.5*phi;
  mpos3 = -0.5*d; mrot3 = -0.5*phi;
  mpos4 =  0.5*d; mrot4 =  0.5*phi;
  mpos5 =  1.5*d; mrot5 =  1.5*phi;
  mpos6 =  2.5*d; mrot6 =  2.5*phi;
  mpos7 =  3.5*d; mrot7 =  3.5*phi;

  OMC1_d = OMC1/60.0;
}
  #undef PHM
  #undef TTM
  #undef TT
  #undef C1
  #undef OMC1
  #undef C2
  #undef C3
  _a1_setpos(); /* type Arm */
  _source_setpos(); /* type Source_simple */
  _slit1_setpos(); /* type Slit */
  _slit2_setpos(); /* type Slit */
  _slit3_setpos(); /* type Slit */
  _focus_mono_setpos(); /* type Arm */
  _m0_setpos(); /* type Monochromator_flat */
  _m1_setpos(); /* type Monochromator_flat */
  _m2_setpos(); /* type Monochromator_flat */
  _m3_setpos(); /* type Monochromator_flat */
  _m4_setpos(); /* type Monochromator_flat */
  _m5_setpos(); /* type Monochromator_flat */
  _m6_setpos(); /* type Monochromator_flat */
  _m7_setpos(); /* type Monochromator_flat */
  _a2_setpos(); /* type Arm */
  _slitMS1_setpos(); /* type Slit */
  _slitMS2_setpos(); /* type Slit */
  _c1_setpos(); /* type Collimator_linear */
  _slitMS3_setpos(); /* type Slit */
  _slitMS4_setpos(); /* type Slit */
  _slitMS5_setpos(); /* type Slit */
  _mon_setpos(); /* type Monitor */
  _emon1_setpos(); /* type E_monitor */
  _sample_setpos(); /* type Slit */
  _slit1mm_setpos(); /* type Slit */
  _a3_setpos(); /* type Arm */
  _c2_setpos(); /* type Collimator_linear */
  _ana_setpos(); /* type Arm */
  _a4_setpos(); /* type Arm */
  _c3_setpos(); /* type Collimator_linear */
  _sng_setpos(); /* type Monitor */
  _Emon2_setpos(); /* type E_monitor */

  /* call iteratively all components INITIALISE */

  class_Source_simple_init(_source);

  class_Slit_init(_slit1);

  class_Slit_init(_slit2);

  class_Slit_init(_slit3);


  class_Monochromator_flat_init(_m0);

  class_Monochromator_flat_init(_m1);

  class_Monochromator_flat_init(_m2);

  class_Monochromator_flat_init(_m3);

  class_Monochromator_flat_init(_m4);

  class_Monochromator_flat_init(_m5);

  class_Monochromator_flat_init(_m6);

  class_Monochromator_flat_init(_m7);


  class_Slit_init(_slitMS1);

  class_Slit_init(_slitMS2);

  class_Collimator_linear_init(_c1);

  class_Slit_init(_slitMS3);

  class_Slit_init(_slitMS4);

  class_Slit_init(_slitMS5);

  class_Monitor_init(_mon);

  class_E_monitor_init(_emon1);

  class_Slit_init(_sample);

  class_Slit_init(_slit1mm);


  class_Collimator_linear_init(_c2);



  class_Collimator_linear_init(_c3);

  class_Monitor_init(_sng);

  class_E_monitor_init(_Emon2);

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
_class_Collimator_linear *class_Collimator_linear_trace(_class_Collimator_linear *_comp
  , _class_particle *_particle) {
  ABSORBED=SCATTERED=RESTORE=0;

  #define xmin (_comp->_parameters.xmin)
  #define xmax (_comp->_parameters.xmax)
  #define ymin (_comp->_parameters.ymin)
  #define ymax (_comp->_parameters.ymax)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define length (_comp->_parameters.length)
  #define divergence (_comp->_parameters.divergence)
  #define transmission (_comp->_parameters.transmission)
  #define divergenceV (_comp->_parameters.divergenceV)
  #define slope (_comp->_parameters.slope)
  #define slopeV (_comp->_parameters.slopeV)
    double phi, dt;

    PROP_Z0;
    if (x<xmin || x>xmax || y<ymin || y>ymax)
      ABSORB;
    dt = length/vz;
    PROP_DT(dt);
    if (x<xmin || x>xmax || y<ymin || y>ymax)
      ABSORB;

    if(slope > 0.0)
    {
      phi = fabs(vx/vz);
      if (phi > slope)
        ABSORB;
      else
        p *= transmission*(1.0 - phi/slope);
      SCATTER;
    }
    if (slopeV > 0) {
      phi = fabs(vy/vz);
      if (phi > slopeV)
        ABSORB;
      else
        p *= transmission*(1.0 - phi/slopeV);
      SCATTER;
    }
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef xwidth
  #undef yheight
  #undef length
  #undef divergence
  #undef transmission
  #undef divergenceV
  #undef slope
  #undef slopeV
  return(_comp);
} /* class_Collimator_linear_trace */

#pragma acc routine seq
_class_Monitor *class_Monitor_trace(_class_Monitor *_comp
  , _class_particle *_particle) {
  ABSORBED=SCATTERED=RESTORE=0;

  #define xmin (_comp->_parameters.xmin)
  #define xmax (_comp->_parameters.xmax)
  #define ymin (_comp->_parameters.ymin)
  #define ymax (_comp->_parameters.ymax)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define restore_neutron (_comp->_parameters.restore_neutron)
  #define Nsum (_comp->_parameters.Nsum)
  #define psum (_comp->_parameters.psum)
  #define p2sum (_comp->_parameters.p2sum)
    PROP_Z0;
    if (x>xmin && x<xmax && y>ymin && y<ymax)
    {
      Nsum++;
      psum += p;
      p2sum += p*p;
      SCATTER;
    }
    if (restore_neutron) {
      RESTORE_NEUTRON(INDEX_CURRENT_COMP, x, y, z, vx, vy, vz, t, sx, sy, sz, p);
    }
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef xwidth
  #undef yheight
  #undef restore_neutron
  #undef Nsum
  #undef psum
  #undef p2sum
  return(_comp);
} /* class_Monitor_trace */

#pragma acc routine seq
_class_E_monitor *class_E_monitor_trace(_class_E_monitor *_comp
  , _class_particle *_particle) {
  ABSORBED=SCATTERED=RESTORE=0;

  #define nE (_comp->_parameters.nE)
  #define filename (_comp->_parameters.filename)
  #define xmin (_comp->_parameters.xmin)
  #define xmax (_comp->_parameters.xmax)
  #define ymin (_comp->_parameters.ymin)
  #define ymax (_comp->_parameters.ymax)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define Emin (_comp->_parameters.Emin)
  #define Emax (_comp->_parameters.Emax)
  #define restore_neutron (_comp->_parameters.restore_neutron)
  #define E_N (_comp->_parameters.E_N)
  #define E_p (_comp->_parameters.E_p)
  #define E_p2 (_comp->_parameters.E_p2)
  #define S_p (_comp->_parameters.S_p)
  #define S_pE (_comp->_parameters.S_pE)
  #define S_pE2 (_comp->_parameters.S_pE2)
  int i;
  double E;

  PROP_Z0;
  if (x>xmin && x<xmax && y>ymin && y<ymax)
  {
    E = VS2E*(vx*vx + vy*vy + vz*vz);

    S_p += p;
    S_pE += p*E;
    S_pE2 += p*E*E;

    i = floor((E-Emin)*nE/(Emax-Emin));
    if(i >= 0 && i < nE)
    {
      E_N[i]++;
      E_p[i] += p;
      E_p2[i] += p*p;
      SCATTER;
    }
  }
  if (restore_neutron) {
    RESTORE_NEUTRON(INDEX_CURRENT_COMP, x, y, z, vx, vy, vz, t, sx, sy, sz, p);
  }
  #undef nE
  #undef filename
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef xwidth
  #undef yheight
  #undef Emin
  #undef Emax
  #undef restore_neutron
  #undef E_N
  #undef E_p
  #undef E_p2
  #undef S_p
  #undef S_pE
  #undef S_pE2
  return(_comp);
} /* class_E_monitor_trace */

/* *****************************************************************************
* instrument 'TAS1_Diff_Slit' TRACE
***************************************************************************** */

#pragma acc routine seq
int raytrace(_class_particle* _particle) { /* called by mccode_main for TAS1_Diff_Slit:TRACE */

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
      /* component source=Source_simple() [2] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_source->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _source->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_source->_position_relative, _source->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Source_simple_trace(_source, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component source [2] */
    if (!ABSORBED && _particle->_index == 3) {
      /* component slit1=Slit() [3] */
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
    } /* end component slit1 [3] */
    if (!ABSORBED && _particle->_index == 4) {
      /* component slit2=Slit() [4] */
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
    } /* end component slit2 [4] */
    if (!ABSORBED && _particle->_index == 5) {
      /* component slit3=Slit() [5] */
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
    } /* end component slit3 [5] */
    if (!ABSORBED && _particle->_index == 6) {
      /* component focus_mono=Arm() [6] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_focus_mono->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _focus_mono->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_focus_mono->_position_relative, _focus_mono->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      _particle->_index++;
    } /* end component focus_mono [6] */
    if (!ABSORBED && _particle->_index == 7) {
      /* component m0=Monochromator_flat() [7] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_m0->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _m0->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_m0->_position_relative, _m0->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Monochromator_flat_trace(_m0, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component m0 [7] */
    if (!ABSORBED && _particle->_index == 8) {
      /* component m1=Monochromator_flat() [8] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_m1->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _m1->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_m1->_position_relative, _m1->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Monochromator_flat_trace(_m1, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component m1 [8] */
    if (!ABSORBED && _particle->_index == 9) {
      /* component m2=Monochromator_flat() [9] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_m2->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _m2->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_m2->_position_relative, _m2->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Monochromator_flat_trace(_m2, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component m2 [9] */
    if (!ABSORBED && _particle->_index == 10) {
      /* component m3=Monochromator_flat() [10] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_m3->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _m3->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_m3->_position_relative, _m3->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Monochromator_flat_trace(_m3, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component m3 [10] */
    if (!ABSORBED && _particle->_index == 11) {
      /* component m4=Monochromator_flat() [11] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_m4->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _m4->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_m4->_position_relative, _m4->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Monochromator_flat_trace(_m4, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component m4 [11] */
    if (!ABSORBED && _particle->_index == 12) {
      /* component m5=Monochromator_flat() [12] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_m5->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _m5->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_m5->_position_relative, _m5->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Monochromator_flat_trace(_m5, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component m5 [12] */
    if (!ABSORBED && _particle->_index == 13) {
      /* component m6=Monochromator_flat() [13] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_m6->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _m6->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_m6->_position_relative, _m6->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Monochromator_flat_trace(_m6, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component m6 [13] */
    if (!ABSORBED && _particle->_index == 14) {
      /* component m7=Monochromator_flat() [14] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_m7->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _m7->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_m7->_position_relative, _m7->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Monochromator_flat_trace(_m7, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component m7 [14] */
    if (!ABSORBED && _particle->_index == 15) {
      /* component a2=Arm() [15] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_a2->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _a2->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_a2->_position_relative, _a2->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      _particle->_index++;
    } /* end component a2 [15] */
    if (!ABSORBED && _particle->_index == 16) {
      /* component slitMS1=Slit() [16] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_slitMS1->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _slitMS1->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_slitMS1->_position_relative, _slitMS1->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Slit_trace(_slitMS1, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component slitMS1 [16] */
    if (!ABSORBED && _particle->_index == 17) {
      /* component slitMS2=Slit() [17] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_slitMS2->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _slitMS2->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_slitMS2->_position_relative, _slitMS2->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Slit_trace(_slitMS2, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component slitMS2 [17] */
    if (!ABSORBED && _particle->_index == 18) {
      /* component c1=Collimator_linear() [18] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_c1->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _c1->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_c1->_position_relative, _c1->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Collimator_linear_trace(_c1, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component c1 [18] */
    if (!ABSORBED && _particle->_index == 19) {
      /* component slitMS3=Slit() [19] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_slitMS3->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _slitMS3->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_slitMS3->_position_relative, _slitMS3->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Slit_trace(_slitMS3, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component slitMS3 [19] */
    if (!ABSORBED && _particle->_index == 20) {
      /* component slitMS4=Slit() [20] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_slitMS4->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _slitMS4->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_slitMS4->_position_relative, _slitMS4->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Slit_trace(_slitMS4, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component slitMS4 [20] */
    if (!ABSORBED && _particle->_index == 21) {
      /* component slitMS5=Slit() [21] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_slitMS5->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _slitMS5->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_slitMS5->_position_relative, _slitMS5->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Slit_trace(_slitMS5, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component slitMS5 [21] */
    if (!ABSORBED && _particle->_index == 22) {
      /* component mon=Monitor() [22] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_mon->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _mon->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_mon->_position_relative, _mon->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Monitor_trace(_mon, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component mon [22] */
    if (!ABSORBED && _particle->_index == 23) {
      /* component emon1=E_monitor() [23] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_emon1->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _emon1->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_emon1->_position_relative, _emon1->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_E_monitor_trace(_emon1, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component emon1 [23] */
    if (!ABSORBED && _particle->_index == 24) {
      /* component sample=Slit() [24] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_sample->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _sample->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_sample->_position_relative, _sample->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Slit_trace(_sample, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component sample [24] */
    if (!ABSORBED && _particle->_index == 25) {
      /* component slit1mm=Slit() [25] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_slit1mm->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _slit1mm->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_slit1mm->_position_relative, _slit1mm->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Slit_trace(_slit1mm, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component slit1mm [25] */
    if (!ABSORBED && _particle->_index == 26) {
      /* component a3=Arm() [26] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_a3->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _a3->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_a3->_position_relative, _a3->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      _particle->_index++;
    } /* end component a3 [26] */
    if (!ABSORBED && _particle->_index == 27) {
      /* component c2=Collimator_linear() [27] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_c2->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _c2->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_c2->_position_relative, _c2->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Collimator_linear_trace(_c2, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component c2 [27] */
    if (!ABSORBED && _particle->_index == 28) {
      /* component ana=Arm() [28] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_ana->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _ana->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_ana->_position_relative, _ana->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      _particle->_index++;
    } /* end component ana [28] */
    if (!ABSORBED && _particle->_index == 29) {
      /* component a4=Arm() [29] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_a4->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _a4->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_a4->_position_relative, _a4->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      _particle->_index++;
    } /* end component a4 [29] */
    if (!ABSORBED && _particle->_index == 30) {
      /* component c3=Collimator_linear() [30] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_c3->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _c3->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_c3->_position_relative, _c3->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Collimator_linear_trace(_c3, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component c3 [30] */
    if (!ABSORBED && _particle->_index == 31) {
      /* component sng=Monitor() [31] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_sng->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _sng->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_sng->_position_relative, _sng->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Monitor_trace(_sng, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component sng [31] */
    if (!ABSORBED && _particle->_index == 32) {
      /* component Emon2=E_monitor() [32] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Emon2->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Emon2->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Emon2->_position_relative, _Emon2->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_E_monitor_trace(_Emon2, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component Emon2 [32] */
    if (_particle->_index > 32)
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
* instrument 'TAS1_Diff_Slit' and components SAVE
***************************************************************************** */

_class_Monitor *class_Monitor_save(_class_Monitor *_comp
) {
  #define xmin (_comp->_parameters.xmin)
  #define xmax (_comp->_parameters.xmax)
  #define ymin (_comp->_parameters.ymin)
  #define ymax (_comp->_parameters.ymax)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define restore_neutron (_comp->_parameters.restore_neutron)
  #define Nsum (_comp->_parameters.Nsum)
  #define psum (_comp->_parameters.psum)
  #define p2sum (_comp->_parameters.p2sum)
    char title[1024];
    sprintf(title, "Single monitor %s", NAME_CURRENT_COMP);
    DETECTOR_OUT_0D(title, Nsum, psum, p2sum);
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef xwidth
  #undef yheight
  #undef restore_neutron
  #undef Nsum
  #undef psum
  #undef p2sum
  return(_comp);
} /* class_Monitor_save */

_class_E_monitor *class_E_monitor_save(_class_E_monitor *_comp
) {
  #define nE (_comp->_parameters.nE)
  #define filename (_comp->_parameters.filename)
  #define xmin (_comp->_parameters.xmin)
  #define xmax (_comp->_parameters.xmax)
  #define ymin (_comp->_parameters.ymin)
  #define ymax (_comp->_parameters.ymax)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define Emin (_comp->_parameters.Emin)
  #define Emax (_comp->_parameters.Emax)
  #define restore_neutron (_comp->_parameters.restore_neutron)
  #define E_N (_comp->_parameters.E_N)
  #define E_p (_comp->_parameters.E_p)
  #define E_p2 (_comp->_parameters.E_p2)
  #define S_p (_comp->_parameters.S_p)
  #define S_pE (_comp->_parameters.S_pE)
  #define S_pE2 (_comp->_parameters.S_pE2)
  DETECTOR_OUT_1D(
      "Energy monitor",
      "Energy [meV]",
      "Intensity",
      "E", Emin, Emax, nE,
      &E_N[0],&E_p[0],&E_p2[0],
      filename);
  if (S_p) printf("<E> : %g meV , E-width : %g meV \n",
   S_pE/S_p,sqrt(S_pE2/S_p - S_pE*S_pE/(S_p*S_p)) );
  #undef nE
  #undef filename
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef xwidth
  #undef yheight
  #undef Emin
  #undef Emax
  #undef restore_neutron
  #undef E_N
  #undef E_p
  #undef E_p2
  #undef S_p
  #undef S_pE
  #undef S_pE2
  return(_comp);
} /* class_E_monitor_save */



int save(FILE *handle) { /* called by mccode_main for TAS1_Diff_Slit:SAVE */
  if (!handle) siminfo_init(NULL);

  /* call iteratively all components SAVE */





















  class_Monitor_save(_mon);

  class_E_monitor_save(_emon1);








  class_Monitor_save(_sng);

  class_E_monitor_save(_Emon2);

  if (!handle) siminfo_close(); 

  return(0);
} /* save */

/* *****************************************************************************
* instrument 'TAS1_Diff_Slit' and components FINALLY
***************************************************************************** */

_class_E_monitor *class_E_monitor_finally(_class_E_monitor *_comp
) {
  #define nE (_comp->_parameters.nE)
  #define filename (_comp->_parameters.filename)
  #define xmin (_comp->_parameters.xmin)
  #define xmax (_comp->_parameters.xmax)
  #define ymin (_comp->_parameters.ymin)
  #define ymax (_comp->_parameters.ymax)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define Emin (_comp->_parameters.Emin)
  #define Emax (_comp->_parameters.Emax)
  #define restore_neutron (_comp->_parameters.restore_neutron)
  #define E_N (_comp->_parameters.E_N)
  #define E_p (_comp->_parameters.E_p)
  #define E_p2 (_comp->_parameters.E_p2)
  #define S_p (_comp->_parameters.S_p)
  #define S_pE (_comp->_parameters.S_pE)
  #define S_pE2 (_comp->_parameters.S_pE2)
  destroy_darr1d(E_N);
  destroy_darr1d(E_p);
  destroy_darr1d(E_p2);
  #undef nE
  #undef filename
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef xwidth
  #undef yheight
  #undef Emin
  #undef Emax
  #undef restore_neutron
  #undef E_N
  #undef E_p
  #undef E_p2
  #undef S_p
  #undef S_pE
  #undef S_pE2
  return(_comp);
} /* class_E_monitor_finally */



int finally(void) { /* called by mccode_main for TAS1_Diff_Slit:FINALLY */
  siminfo_init(NULL);
  save(siminfo_file); /* save data when simulation ends */

  /* call iteratively all components FINALLY */






















  class_E_monitor_finally(_emon1);









  class_E_monitor_finally(_Emon2);

  siminfo_close(); 

  return(0);
} /* finally */

/* *****************************************************************************
* instrument 'TAS1_Diff_Slit' and components DISPLAY
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

_class_Collimator_linear *class_Collimator_linear_display(_class_Collimator_linear *_comp
) {
  #define xmin (_comp->_parameters.xmin)
  #define xmax (_comp->_parameters.xmax)
  #define ymin (_comp->_parameters.ymin)
  #define ymax (_comp->_parameters.ymax)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define length (_comp->_parameters.length)
  #define divergence (_comp->_parameters.divergence)
  #define transmission (_comp->_parameters.transmission)
  #define divergenceV (_comp->_parameters.divergenceV)
  #define slope (_comp->_parameters.slope)
  #define slopeV (_comp->_parameters.slopeV)
  double x;
  int i;

  
  for(x = xmin, i = 0; i <= 3; i++, x += (xmax - xmin)/3.0)
    multiline(5, x, (double)ymin, 0.0, x, (double)ymax, 0.0,
              x, (double)ymax, (double)length, x, (double)ymin, (double)length,
              x, (double)ymin, 0.0);
  line(xmin, ymin, 0,   xmax, ymin, 0);
  line(xmin, ymax, 0,   xmax, ymax, 0);
  line(xmin, ymin, length, xmax, ymin, length);
  line(xmin, ymax, length, xmax, ymax, length);
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef xwidth
  #undef yheight
  #undef length
  #undef divergence
  #undef transmission
  #undef divergenceV
  #undef slope
  #undef slopeV
  return(_comp);
} /* class_Collimator_linear_display */

_class_Monitor *class_Monitor_display(_class_Monitor *_comp
) {
  #define xmin (_comp->_parameters.xmin)
  #define xmax (_comp->_parameters.xmax)
  #define ymin (_comp->_parameters.ymin)
  #define ymax (_comp->_parameters.ymax)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define restore_neutron (_comp->_parameters.restore_neutron)
  #define Nsum (_comp->_parameters.Nsum)
  #define psum (_comp->_parameters.psum)
  #define p2sum (_comp->_parameters.p2sum)
  
  multiline(5, (double)xmin, (double)ymin, 0.0,
               (double)xmax, (double)ymin, 0.0,
               (double)xmax, (double)ymax, 0.0,
               (double)xmin, (double)ymax, 0.0,
               (double)xmin, (double)ymin, 0.0);
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef xwidth
  #undef yheight
  #undef restore_neutron
  #undef Nsum
  #undef psum
  #undef p2sum
  return(_comp);
} /* class_Monitor_display */

_class_E_monitor *class_E_monitor_display(_class_E_monitor *_comp
) {
  #define nE (_comp->_parameters.nE)
  #define filename (_comp->_parameters.filename)
  #define xmin (_comp->_parameters.xmin)
  #define xmax (_comp->_parameters.xmax)
  #define ymin (_comp->_parameters.ymin)
  #define ymax (_comp->_parameters.ymax)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define Emin (_comp->_parameters.Emin)
  #define Emax (_comp->_parameters.Emax)
  #define restore_neutron (_comp->_parameters.restore_neutron)
  #define E_N (_comp->_parameters.E_N)
  #define E_p (_comp->_parameters.E_p)
  #define E_p2 (_comp->_parameters.E_p2)
  #define S_p (_comp->_parameters.S_p)
  #define S_pE (_comp->_parameters.S_pE)
  #define S_pE2 (_comp->_parameters.S_pE2)

  multiline(5, (double)xmin, (double)ymin, 0.0,
               (double)xmax, (double)ymin, 0.0,
               (double)xmax, (double)ymax, 0.0,
               (double)xmin, (double)ymax, 0.0,
               (double)xmin, (double)ymin, 0.0);
  #undef nE
  #undef filename
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef xwidth
  #undef yheight
  #undef Emin
  #undef Emax
  #undef restore_neutron
  #undef E_N
  #undef E_p
  #undef E_p2
  #undef S_p
  #undef S_pE
  #undef S_pE2
  return(_comp);
} /* class_E_monitor_display */


  #undef magnify
  #undef line
  #undef dashed_line
  #undef multiline
  #undef rectangle
  #undef box
  #undef circle
  #undef cylinder
  #undef sphere

int display(void) { /* called by mccode_main for TAS1_Diff_Slit:DISPLAY */
  printf("MCDISPLAY: start\n");

  /* call iteratively all components DISPLAY */
  class_Arm_display(_a1);

  class_Source_simple_display(_source);

  class_Slit_display(_slit1);

  class_Slit_display(_slit2);

  class_Slit_display(_slit3);

  class_Arm_display(_focus_mono);

  class_Monochromator_flat_display(_m0);

  class_Monochromator_flat_display(_m1);

  class_Monochromator_flat_display(_m2);

  class_Monochromator_flat_display(_m3);

  class_Monochromator_flat_display(_m4);

  class_Monochromator_flat_display(_m5);

  class_Monochromator_flat_display(_m6);

  class_Monochromator_flat_display(_m7);

  class_Arm_display(_a2);

  class_Slit_display(_slitMS1);

  class_Slit_display(_slitMS2);

  class_Collimator_linear_display(_c1);

  class_Slit_display(_slitMS3);

  class_Slit_display(_slitMS4);

  class_Slit_display(_slitMS5);

  class_Monitor_display(_mon);

  class_E_monitor_display(_emon1);

  class_Slit_display(_sample);

  class_Slit_display(_slit1mm);

  class_Arm_display(_a3);

  class_Collimator_linear_display(_c2);

  class_Arm_display(_ana);

  class_Arm_display(_a4);

  class_Collimator_linear_display(_c3);

  class_Monitor_display(_sng);

  class_E_monitor_display(_Emon2);

  printf("MCDISPLAY: end\n");

  return(0);
} /* display */

void* _getvar_parameters(char* compname)
/* enables settings parameters based use of the GETPAR macro */
{
  if (!strcmp(compname, "a1")) return (void *) &(_a1->_parameters);
  if (!strcmp(compname, "source")) return (void *) &(_source->_parameters);
  if (!strcmp(compname, "slit1")) return (void *) &(_slit1->_parameters);
  if (!strcmp(compname, "slit2")) return (void *) &(_slit2->_parameters);
  if (!strcmp(compname, "slit3")) return (void *) &(_slit3->_parameters);
  if (!strcmp(compname, "focus_mono")) return (void *) &(_focus_mono->_parameters);
  if (!strcmp(compname, "m0")) return (void *) &(_m0->_parameters);
  if (!strcmp(compname, "m1")) return (void *) &(_m1->_parameters);
  if (!strcmp(compname, "m2")) return (void *) &(_m2->_parameters);
  if (!strcmp(compname, "m3")) return (void *) &(_m3->_parameters);
  if (!strcmp(compname, "m4")) return (void *) &(_m4->_parameters);
  if (!strcmp(compname, "m5")) return (void *) &(_m5->_parameters);
  if (!strcmp(compname, "m6")) return (void *) &(_m6->_parameters);
  if (!strcmp(compname, "m7")) return (void *) &(_m7->_parameters);
  if (!strcmp(compname, "a2")) return (void *) &(_a2->_parameters);
  if (!strcmp(compname, "slitMS1")) return (void *) &(_slitMS1->_parameters);
  if (!strcmp(compname, "slitMS2")) return (void *) &(_slitMS2->_parameters);
  if (!strcmp(compname, "c1")) return (void *) &(_c1->_parameters);
  if (!strcmp(compname, "slitMS3")) return (void *) &(_slitMS3->_parameters);
  if (!strcmp(compname, "slitMS4")) return (void *) &(_slitMS4->_parameters);
  if (!strcmp(compname, "slitMS5")) return (void *) &(_slitMS5->_parameters);
  if (!strcmp(compname, "mon")) return (void *) &(_mon->_parameters);
  if (!strcmp(compname, "emon1")) return (void *) &(_emon1->_parameters);
  if (!strcmp(compname, "sample")) return (void *) &(_sample->_parameters);
  if (!strcmp(compname, "slit1mm")) return (void *) &(_slit1mm->_parameters);
  if (!strcmp(compname, "a3")) return (void *) &(_a3->_parameters);
  if (!strcmp(compname, "c2")) return (void *) &(_c2->_parameters);
  if (!strcmp(compname, "ana")) return (void *) &(_ana->_parameters);
  if (!strcmp(compname, "a4")) return (void *) &(_a4->_parameters);
  if (!strcmp(compname, "c3")) return (void *) &(_c3->_parameters);
  if (!strcmp(compname, "sng")) return (void *) &(_sng->_parameters);
  if (!strcmp(compname, "Emon2")) return (void *) &(_Emon2->_parameters);
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

/* end of generated C code linup-3.c */
