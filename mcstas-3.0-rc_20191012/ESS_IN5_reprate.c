/* Automatically generated file. Do not edit. 
 * Format:     ANSI C source code
 * Creator:    McStas <http://www.mcstas.org>
 * Instrument: ESS_IN5_reprate.instr (ESS_IN5_reprate)
 * Date:       Sat Oct 12 09:05:01 2019
 * File:       ESS_IN5_reprate.c
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
* Start of instrument 'ESS_IN5_reprate' generated code
***************************************************************************** */

#ifdef MC_TRACE_ENABLED
int traceenabled = 1;
#else
int traceenabled = 0;
#endif
#define MCSTAS "/usr/share/mcstas/3.0-dev/"
int   defaultmain         = 1;
char  instrument_name[]   = "ESS_IN5_reprate";
char  instrument_source[] = "ESS_IN5_reprate.instr";
char *instrument_exe      = NULL; /* will be set to argv[0] in main */
char  instrument_code[]   = "Instrument ESS_IN5_reprate source code ESS_IN5_reprate.instr is not embedded in this executable.\n  Use --source option when running McStas.\n";

int main(int argc, char *argv[]){return mccode_main(argc, argv);}

/* *****************************************************************************
* instrument 'ESS_IN5_reprate' and components DECLARE
***************************************************************************** */

/* Instrument parameters: structure and a table for the initialisation
   (Used in e.g. inputparse and I/O function (e.g. detector_out) */

struct _struct_instrument_parameters {
  MCNUM _Lmin;
  MCNUM _Lmax;
  MCNUM _lambda0;
  MCNUM _Pulse_width;
  MCNUM _Num_pulses;
  MCNUM _GUI_start;
  MCNUM _FO1_DIST;
  MCNUM _L_ballistic_begin;
  MCNUM _L_ballistic_end;
  MCNUM _Length;
  MCNUM _SAMPLE_DIST;
  MCNUM _DETECTOR_DIST;
  MCNUM _GUI_h;
  MCNUM _GUI_w;
  MCNUM _GUI_GAP;
  MCNUM _H1;
  MCNUM _W1;
  MCNUM _H2;
  MCNUM _W2;
  MCNUM _H3;
  MCNUM _W3;
  MCNUM _H4;
  MCNUM _W4;
  MCNUM _H_chop;
  MCNUM _W_chop;
  MCNUM _H_end;
  MCNUM _W_end;
  MCNUM _ALPHA;
  MCNUM _M;
  MCNUM _F_slow1;
  MCNUM _F_slow2;
  MCNUM _F_fast1;
  MCNUM _F_fast2;
  MCNUM _N_fast;
  MCNUM _SLOW1_THETA;
  MCNUM _FO3;
  MCNUM _THETA_fast1;
  MCNUM _FAST_THETA;
  MCNUM _Gamma;
  MCNUM _Etun;
  MCNUM _V_HOLE;
  MCNUM _FRAC_QUASIEL;
  MCNUM _FRAC_TUNNEL;
  MCNUM _TT;
  MCNUM _RES_DE;
  MCNUM _port;
  MCNUM _cold;
};
typedef struct _struct_instrument_parameters _class_instrument_parameters;

struct _instrument_struct {
  char   _name[256]; /* the name of this instrument e.g. 'ESS_IN5_reprate' */
/* Counters per component instance */
  double counter_AbsorbProp[52]; /* absorbed events in PROP routines */
  double counter_N[52], counter_P[52], counter_P2[52]; /* event counters after each component instance */
  _class_particle _trajectory[52]; /* current trajectory for STORE/RESTORE */
/* Components position table (absolute and relative coords) */
  Coords _position_relative[52]; /* positions of all components */
  Coords _position_absolute[52];
  _class_instrument_parameters _parameters; /* instrument parameters */
} _instrument_var;
struct _instrument_struct *instrument = & _instrument_var;
#pragma acc declare create ( _instrument_var )
#pragma acc declare create ( instrument )

int numipar = 47;
struct mcinputtable_struct mcinputtable[] = {
  "Lmin", &(_instrument_var._parameters._Lmin), instr_type_double, "4.9", 
  "Lmax", &(_instrument_var._parameters._Lmax), instr_type_double, "5.1", 
  "lambda0", &(_instrument_var._parameters._lambda0), instr_type_double, "5", 
  "Pulse_width", &(_instrument_var._parameters._Pulse_width), instr_type_double, "0.002", 
  "Num_pulses", &(_instrument_var._parameters._Num_pulses), instr_type_double, "1", 
  "GUI_start", &(_instrument_var._parameters._GUI_start), instr_type_double, "2.0", 
  "FO1_DIST", &(_instrument_var._parameters._FO1_DIST), instr_type_double, "6", 
  "L_ballistic_begin", &(_instrument_var._parameters._L_ballistic_begin), instr_type_double, "19.5", 
  "L_ballistic_end", &(_instrument_var._parameters._L_ballistic_end), instr_type_double, "17", 
  "Length", &(_instrument_var._parameters._Length), instr_type_double, "100", 
  "SAMPLE_DIST", &(_instrument_var._parameters._SAMPLE_DIST), instr_type_double, "1.2", 
  "DETECTOR_DIST", &(_instrument_var._parameters._DETECTOR_DIST), instr_type_double, "4", 
  "GUI_h", &(_instrument_var._parameters._GUI_h), instr_type_double, "0.105", 
  "GUI_w", &(_instrument_var._parameters._GUI_w), instr_type_double, "0.1", 
  "GUI_GAP", &(_instrument_var._parameters._GUI_GAP), instr_type_double, "0.05", 
  "H1", &(_instrument_var._parameters._H1), instr_type_double, "0.167", 
  "W1", &(_instrument_var._parameters._W1), instr_type_double, "0.116", 
  "H2", &(_instrument_var._parameters._H2), instr_type_double, "0.185", 
  "W2", &(_instrument_var._parameters._W2), instr_type_double, "0.15", 
  "H3", &(_instrument_var._parameters._H3), instr_type_double, "0.19", 
  "W3", &(_instrument_var._parameters._W3), instr_type_double, "0.15", 
  "H4", &(_instrument_var._parameters._H4), instr_type_double, "0.213", 
  "W4", &(_instrument_var._parameters._W4), instr_type_double, "0.14", 
  "H_chop", &(_instrument_var._parameters._H_chop), instr_type_double, "0.075", 
  "W_chop", &(_instrument_var._parameters._W_chop), instr_type_double, "0.03", 
  "H_end", &(_instrument_var._parameters._H_end), instr_type_double, "0.042", 
  "W_end", &(_instrument_var._parameters._W_end), instr_type_double, "0.0215", 
  "ALPHA", &(_instrument_var._parameters._ALPHA), instr_type_double, "3.4", 
  "M", &(_instrument_var._parameters._M), instr_type_double, "3.5", 
  "F_slow1", &(_instrument_var._parameters._F_slow1), instr_type_double, "16.6667", 
  "F_slow2", &(_instrument_var._parameters._F_slow2), instr_type_double, "0", 
  "F_fast1", &(_instrument_var._parameters._F_fast1), instr_type_double, "0", 
  "F_fast2", &(_instrument_var._parameters._F_fast2), instr_type_double, "200", 
  "N_fast", &(_instrument_var._parameters._N_fast), instr_type_double, "1", 
  "SLOW1_THETA", &(_instrument_var._parameters._SLOW1_THETA), instr_type_double, "120", 
  "FO3", &(_instrument_var._parameters._FO3), instr_type_double, "1", 
  "THETA_fast1", &(_instrument_var._parameters._THETA_fast1), instr_type_double, "180", 
  "FAST_THETA", &(_instrument_var._parameters._FAST_THETA), instr_type_double, "5", 
  "Gamma", &(_instrument_var._parameters._Gamma), instr_type_double, "0", 
  "Etun", &(_instrument_var._parameters._Etun), instr_type_double, "1", 
  "V_HOLE", &(_instrument_var._parameters._V_HOLE), instr_type_double, "0", 
  "FRAC_QUASIEL", &(_instrument_var._parameters._FRAC_QUASIEL), instr_type_double, "0", 
  "FRAC_TUNNEL", &(_instrument_var._parameters._FRAC_TUNNEL), instr_type_double, "0", 
  "TT", &(_instrument_var._parameters._TT), instr_type_double, "50", 
  "RES_DE", &(_instrument_var._parameters._RES_DE), instr_type_double, "0.5", 
  "port", &(_instrument_var._parameters._port), instr_type_double, "30", 
  "cold", &(_instrument_var._parameters._cold), instr_type_double, "0.95", 
  NULL, NULL, instr_type_double, ""
};


/* ************************************************************************** */
/*             SHARE user declarations for all components                     */
/* ************************************************************************** */

/* Shared user declarations for all components types 'ESS_moderator_long'. */
double Mezei_M_fct(double l, double temp)
  {
    double a=949.0/temp;
    return 2*a*a*exp(-a/(l*l))/(l*l*l*l*l);
  }

  double Mezei_F_fct(double t, double tau, int n)
  {
    return (exp(-t/tau)-exp(-n*t/tau))*n/(n-1)/tau;
  }

/* Shared user declarations for all components types 'Guide'. */
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


/* Shared user declarations for all components types 'Tunneling_sample'. */
struct StructVarsV
{
double sigma_a; /* Absorption cross section per atom (barns) */
    double sigma_i; /* Incoherent scattering cross section per atom (barns) */
    double rho;     /* Density of atoms (AA-3) */
    double my_s;
    double my_a_v;
    char   isrect;      /* true when sample is a box */
    double distance;    /* when non zero, gives rect target distance */
    double aw,ah;       /* rectangular angular dimensions */
    double xw,yh;       /* rectangular metrical dimensions */
    double tx,ty,tz;    /* target coords */
  };


/* ************************************************************************** */
/*             End of SHARE user declarations for all components              */
/* ************************************************************************** */


/* ********************** component definition declarations. **************** */

/* component source=ESS_moderator_long() [1] DECLARE */
/* Parameter definition for component type 'ESS_moderator_long' */
struct _struct_ESS_moderator_long_parameters {
  /* Component type 'ESS_moderator_long' setting parameters */
  MCNUM width_c;
  MCNUM yheight;
  MCNUM Lmin;
  MCNUM Lmax;
  MCNUM dist;
  MCNUM focus_xw;
  MCNUM focus_yh;
  MCNUM nu;
  MCNUM T;
  MCNUM tau;
  MCNUM tau1;
  MCNUM tau2;
  MCNUM d;
  MCNUM n;
  MCNUM cold_frac;
  MCNUM n2;
  MCNUM chi2;
  MCNUM I0;
  MCNUM I2;
  long target_index;
  MCNUM cyl_radius;
  MCNUM branch1;
  MCNUM branch2;
  MCNUM branch_tail;
  long n_pulses;
  MCNUM width_t;
  MCNUM T_t;
  MCNUM tau_t;
  MCNUM tau1_t;
  MCNUM tau2_t;
  MCNUM chi2_t;
  MCNUM I0_t;
  MCNUM I2_t;
  MCNUM branch1_t;
  MCNUM branch2_t;
  long src_2012;
  MCNUM tfocus_dist;
  MCNUM tfocus_time;
  MCNUM tfocus_width;
  MCNUM beamport_angle;
  /* Component type 'ESS_moderator_long' private parameters */
  /* Component type 'ESS_moderator_long' DECLARE code stored as structure members */
  double l_range, w_mult, w_geom, w_geom_c, w_geom_t;
  double tx,ty,tz;
  double t1x,t1y,t1z,t2x,t2y,t2z;
                                                     
  double T_n, tau_n, tau1_n, tau2_n, chi2_n, I0_n, I2_n, branch1_n, branch2_n;

                                  
  double r_empty;                                                   
  double r_optics;
}; /* _struct_ESS_moderator_long_parameters */
typedef struct _struct_ESS_moderator_long_parameters _class_ESS_moderator_long_parameters;

/* Parameters for component type 'ESS_moderator_long' */
struct _struct_ESS_moderator_long {
  char     _name[256]; /* e.g. source */
  char     _type[256]; /* ESS_moderator_long */
  long     _index; /* e.g. 2 index in TRACE list */
  Coords   _position_absolute;
  Coords   _position_relative; /* wrt PREVIOUS */
  Rotation _rotation_absolute;
  Rotation _rotation_relative; /* wrt PREVIOUS */
  int      _rotation_is_identity;
  _class_ESS_moderator_long_parameters _parameters;
};
typedef struct _struct_ESS_moderator_long _class_ESS_moderator_long;
_class_ESS_moderator_long _source_var;
_class_ESS_moderator_long *_source = &_source_var;
#pragma acc declare create ( _source_var )
#pragma acc declare create ( _source )

/* component Origin=Progress_bar() [2] DECLARE */
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

/* component TOFmoderator_zoom=TOF_monitor() [3] DECLARE */
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
  char     _name[256]; /* e.g. TOFmoderator_zoom */
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
_class_TOF_monitor _TOFmoderator_zoom_var;
_class_TOF_monitor *_TOFmoderator_zoom = &_TOFmoderator_zoom_var;
#pragma acc declare create ( _TOFmoderator_zoom_var )
#pragma acc declare create ( _TOFmoderator_zoom )

_class_TOF_monitor _TOFmoderator_var;
_class_TOF_monitor *_TOFmoderator = &_TOFmoderator_var;
#pragma acc declare create ( _TOFmoderator_var )
#pragma acc declare create ( _TOFmoderator )

/* component Lmon_guistart=L_monitor() [5] DECLARE */
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
  char     _name[256]; /* e.g. Lmon_guistart */
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
_class_L_monitor _Lmon_guistart_var;
_class_L_monitor *_Lmon_guistart = &_Lmon_guistart_var;
#pragma acc declare create ( _Lmon_guistart_var )
#pragma acc declare create ( _Lmon_guistart )

_class_L_monitor _Lmon_normalize_var;
_class_L_monitor *_Lmon_normalize = &_Lmon_normalize_var;
#pragma acc declare create ( _Lmon_normalize_var )
#pragma acc declare create ( _Lmon_normalize )

/* component Guide1=Guide() [7] DECLARE */
/* Parameter definition for component type 'Guide' */
struct _struct_Guide_parameters {
  /* Component type 'Guide' setting parameters */
  char reflect[16384];
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
  /* Component type 'Guide' private parameters */
  /* Component type 'Guide' DECLARE code stored as structure members */
t_Table pTable;
}; /* _struct_Guide_parameters */
typedef struct _struct_Guide_parameters _class_Guide_parameters;

/* Parameters for component type 'Guide' */
struct _struct_Guide {
  char     _name[256]; /* e.g. Guide1 */
  char     _type[256]; /* Guide */
  long     _index; /* e.g. 2 index in TRACE list */
  Coords   _position_absolute;
  Coords   _position_relative; /* wrt PREVIOUS */
  Rotation _rotation_absolute;
  Rotation _rotation_relative; /* wrt PREVIOUS */
  int      _rotation_is_identity;
  _class_Guide_parameters _parameters;
};
typedef struct _struct_Guide _class_Guide;
_class_Guide _Guide1_var;
_class_Guide *_Guide1 = &_Guide1_var;
#pragma acc declare create ( _Guide1_var )
#pragma acc declare create ( _Guide1 )

_class_L_monitor _Lmonslow1_var;
_class_L_monitor *_Lmonslow1 = &_Lmonslow1_var;
#pragma acc declare create ( _Lmonslow1_var )
#pragma acc declare create ( _Lmonslow1 )

/* component PSDslow1=PSD_monitor() [9] DECLARE */
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
  char     _name[256]; /* e.g. PSDslow1 */
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
_class_PSD_monitor _PSDslow1_var;
_class_PSD_monitor *_PSDslow1 = &_PSDslow1_var;
#pragma acc declare create ( _PSDslow1_var )
#pragma acc declare create ( _PSDslow1 )

/* component FOchop1=DiskChopper() [10] DECLARE */
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
  char     _name[256]; /* e.g. FOchop1 */
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
_class_DiskChopper _FOchop1_var;
_class_DiskChopper *_FOchop1 = &_FOchop1_var;
#pragma acc declare create ( _FOchop1_var )
#pragma acc declare create ( _FOchop1 )

/* component TOFLmon1=TOFLambda_monitor() [11] DECLARE */
/* Parameter definition for component type 'TOFLambda_monitor' */
struct _struct_TOFLambda_monitor_parameters {
  /* Component type 'TOFLambda_monitor' setting parameters */
  MCNUM nL;
  MCNUM nt;
  MCNUM tmin;
  MCNUM tmax;
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
  /* Component type 'TOFLambda_monitor' private parameters */
  /* Component type 'TOFLambda_monitor' DECLARE code stored as structure members */
  DArray2d TOFL_N;
  DArray2d TOFL_p;
  DArray2d TOFL_p2;
  double tt_0, tt_1;
}; /* _struct_TOFLambda_monitor_parameters */
typedef struct _struct_TOFLambda_monitor_parameters _class_TOFLambda_monitor_parameters;

/* Parameters for component type 'TOFLambda_monitor' */
struct _struct_TOFLambda_monitor {
  char     _name[256]; /* e.g. TOFLmon1 */
  char     _type[256]; /* TOFLambda_monitor */
  long     _index; /* e.g. 2 index in TRACE list */
  Coords   _position_absolute;
  Coords   _position_relative; /* wrt PREVIOUS */
  Rotation _rotation_absolute;
  Rotation _rotation_relative; /* wrt PREVIOUS */
  int      _rotation_is_identity;
  _class_TOFLambda_monitor_parameters _parameters;
};
typedef struct _struct_TOFLambda_monitor _class_TOFLambda_monitor;
_class_TOFLambda_monitor _TOFLmon1_var;
_class_TOFLambda_monitor *_TOFLmon1 = &_TOFLmon1_var;
#pragma acc declare create ( _TOFLmon1_var )
#pragma acc declare create ( _TOFLmon1 )

_class_L_monitor _Lmon_afterslow1_var;
_class_L_monitor *_Lmon_afterslow1 = &_Lmon_afterslow1_var;
#pragma acc declare create ( _Lmon_afterslow1_var )
#pragma acc declare create ( _Lmon_afterslow1 )

_class_PSD_monitor _PSD_afterslow1_var;
_class_PSD_monitor *_PSD_afterslow1 = &_PSD_afterslow1_var;
#pragma acc declare create ( _PSD_afterslow1_var )
#pragma acc declare create ( _PSD_afterslow1 )

_class_Guide _Guidelong1_var;
_class_Guide *_Guidelong1 = &_Guidelong1_var;
#pragma acc declare create ( _Guidelong1_var )
#pragma acc declare create ( _Guidelong1 )

_class_Guide _Guidelong1b_var;
_class_Guide *_Guidelong1b = &_Guidelong1b_var;
#pragma acc declare create ( _Guidelong1b_var )
#pragma acc declare create ( _Guidelong1b )

_class_L_monitor _Lmon_slow2_var;
_class_L_monitor *_Lmon_slow2 = &_Lmon_slow2_var;
#pragma acc declare create ( _Lmon_slow2_var )
#pragma acc declare create ( _Lmon_slow2 )

_class_DiskChopper _FOchop2_var;
_class_DiskChopper *_FOchop2 = &_FOchop2_var;
#pragma acc declare create ( _FOchop2_var )
#pragma acc declare create ( _FOchop2 )

_class_DiskChopper _Fastchop1_var;
_class_DiskChopper *_Fastchop1 = &_Fastchop1_var;
#pragma acc declare create ( _Fastchop1_var )
#pragma acc declare create ( _Fastchop1 )

_class_PSD_monitor _PSD_afterslow2_var;
_class_PSD_monitor *_PSD_afterslow2 = &_PSD_afterslow2_var;
#pragma acc declare create ( _PSD_afterslow2_var )
#pragma acc declare create ( _PSD_afterslow2 )

_class_L_monitor _Lmon_afterslow2_var;
_class_L_monitor *_Lmon_afterslow2 = &_Lmon_afterslow2_var;
#pragma acc declare create ( _Lmon_afterslow2_var )
#pragma acc declare create ( _Lmon_afterslow2 )

_class_TOFLambda_monitor _TOFL_afterslow2_var;
_class_TOFLambda_monitor *_TOFL_afterslow2 = &_TOFL_afterslow2_var;
#pragma acc declare create ( _TOFL_afterslow2_var )
#pragma acc declare create ( _TOFL_afterslow2 )

_class_Guide _Guidelong2_var;
_class_Guide *_Guidelong2 = &_Guidelong2_var;
#pragma acc declare create ( _Guidelong2_var )
#pragma acc declare create ( _Guidelong2 )

_class_L_monitor _Lmon_beforeballistic_var;
_class_L_monitor *_Lmon_beforeballistic = &_Lmon_beforeballistic_var;
#pragma acc declare create ( _Lmon_beforeballistic_var )
#pragma acc declare create ( _Lmon_beforeballistic )

_class_PSD_monitor _PSD_beforeballistic_var;
_class_PSD_monitor *_PSD_beforeballistic = &_PSD_beforeballistic_var;
#pragma acc declare create ( _PSD_beforeballistic_var )
#pragma acc declare create ( _PSD_beforeballistic )

_class_Guide _Guidelong2a_var;
_class_Guide *_Guidelong2a = &_Guidelong2a_var;
#pragma acc declare create ( _Guidelong2a_var )
#pragma acc declare create ( _Guidelong2a )

_class_L_monitor _Lmonfast2_var;
_class_L_monitor *_Lmonfast2 = &_Lmonfast2_var;
#pragma acc declare create ( _Lmonfast2_var )
#pragma acc declare create ( _Lmonfast2 )

_class_L_monitor _Lmonfast2_zoom_var;
_class_L_monitor *_Lmonfast2_zoom = &_Lmonfast2_zoom_var;
#pragma acc declare create ( _Lmonfast2_zoom_var )
#pragma acc declare create ( _Lmonfast2_zoom )

_class_TOFLambda_monitor _TOFLfast2_var;
_class_TOFLambda_monitor *_TOFLfast2 = &_TOFLfast2_var;
#pragma acc declare create ( _TOFLfast2_var )
#pragma acc declare create ( _TOFLfast2 )

_class_TOFLambda_monitor _TOFLfast2zoom_var;
_class_TOFLambda_monitor *_TOFLfast2zoom = &_TOFLfast2zoom_var;
#pragma acc declare create ( _TOFLfast2zoom_var )
#pragma acc declare create ( _TOFLfast2zoom )

_class_PSD_monitor _PSDfast2_var;
_class_PSD_monitor *_PSDfast2 = &_PSDfast2_var;
#pragma acc declare create ( _PSDfast2_var )
#pragma acc declare create ( _PSDfast2 )

_class_DiskChopper _Fastchop2_var;
_class_DiskChopper *_Fastchop2 = &_Fastchop2_var;
#pragma acc declare create ( _Fastchop2_var )
#pragma acc declare create ( _Fastchop2 )

_class_DiskChopper _Fastchop2counter_var;
_class_DiskChopper *_Fastchop2counter = &_Fastchop2counter_var;
#pragma acc declare create ( _Fastchop2counter_var )
#pragma acc declare create ( _Fastchop2counter )

_class_DiskChopper _FOchop3_var;
_class_DiskChopper *_FOchop3 = &_FOchop3_var;
#pragma acc declare create ( _FOchop3_var )
#pragma acc declare create ( _FOchop3 )

_class_TOF_monitor _TOFfast2_zoom_var;
_class_TOF_monitor *_TOFfast2_zoom = &_TOFfast2_zoom_var;
#pragma acc declare create ( _TOFfast2_zoom_var )
#pragma acc declare create ( _TOFfast2_zoom )

_class_L_monitor _Lmon_afterfast2_var;
_class_L_monitor *_Lmon_afterfast2 = &_Lmon_afterfast2_var;
#pragma acc declare create ( _Lmon_afterfast2_var )
#pragma acc declare create ( _Lmon_afterfast2 )

_class_TOFLambda_monitor _TOFL_afterfast2_var;
_class_TOFLambda_monitor *_TOFL_afterfast2 = &_TOFL_afterfast2_var;
#pragma acc declare create ( _TOFL_afterfast2_var )
#pragma acc declare create ( _TOFL_afterfast2 )

_class_TOFLambda_monitor _TOFL_afterfast2_zoom_var;
_class_TOFLambda_monitor *_TOFL_afterfast2_zoom = &_TOFL_afterfast2_zoom_var;
#pragma acc declare create ( _TOFL_afterfast2_zoom_var )
#pragma acc declare create ( _TOFL_afterfast2_zoom )

_class_PSD_monitor _PSD_afterfast2_var;
_class_PSD_monitor *_PSD_afterfast2 = &_PSD_afterfast2_var;
#pragma acc declare create ( _PSD_afterfast2_var )
#pragma acc declare create ( _PSD_afterfast2 )

_class_Guide _Guidesample_var;
_class_Guide *_Guidesample = &_Guidesample_var;
#pragma acc declare create ( _Guidesample_var )
#pragma acc declare create ( _Guidesample )

_class_L_monitor _Lmon_guideend_var;
_class_L_monitor *_Lmon_guideend = &_Lmon_guideend_var;
#pragma acc declare create ( _Lmon_guideend_var )
#pragma acc declare create ( _Lmon_guideend )

_class_PSD_monitor _PSDsample_var;
_class_PSD_monitor *_PSDsample = &_PSDsample_var;
#pragma acc declare create ( _PSDsample_var )
#pragma acc declare create ( _PSDsample )

_class_TOF_monitor _TOFsample_zoom_var;
_class_TOF_monitor *_TOFsample_zoom = &_TOFsample_zoom_var;
#pragma acc declare create ( _TOFsample_zoom_var )
#pragma acc declare create ( _TOFsample_zoom )

/* component Esample=E_monitor() [43] DECLARE */
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
  char     _name[256]; /* e.g. Esample */
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
_class_E_monitor _Esample_var;
_class_E_monitor *_Esample = &_Esample_var;
#pragma acc declare create ( _Esample_var )
#pragma acc declare create ( _Esample )

_class_L_monitor _Lmon_sample_zoom_var;
_class_L_monitor *_Lmon_sample_zoom = &_Lmon_sample_zoom_var;
#pragma acc declare create ( _Lmon_sample_zoom_var )
#pragma acc declare create ( _Lmon_sample_zoom )

/* component sample=Tunneling_sample() [45] DECLARE */
/* Parameter definition for component type 'Tunneling_sample' */
struct _struct_Tunneling_sample_parameters {
  /* Component type 'Tunneling_sample' setting parameters */
  MCNUM thickness;
  MCNUM radius;
  MCNUM focus_r;
  MCNUM p_interact;
  MCNUM f_QE;
  MCNUM f_tun;
  MCNUM gamma;
  MCNUM E_tun;
  MCNUM target_x;
  MCNUM target_y;
  MCNUM target_z;
  MCNUM focus_xw;
  MCNUM focus_yh;
  MCNUM focus_aw;
  MCNUM focus_ah;
  MCNUM xwidth;
  MCNUM yheight;
  MCNUM zdepth;
  MCNUM sigma_abs;
  MCNUM sigma_inc;
  MCNUM Vc;
  long target_index;
  /* Component type 'Tunneling_sample' private parameters */
  /* Component type 'Tunneling_sample' DECLARE code stored as structure members */
  struct StructVarsV VarsV;
  double ftun, fQE;
}; /* _struct_Tunneling_sample_parameters */
typedef struct _struct_Tunneling_sample_parameters _class_Tunneling_sample_parameters;

/* Parameters for component type 'Tunneling_sample' */
struct _struct_Tunneling_sample {
  char     _name[256]; /* e.g. sample */
  char     _type[256]; /* Tunneling_sample */
  long     _index; /* e.g. 2 index in TRACE list */
  Coords   _position_absolute;
  Coords   _position_relative; /* wrt PREVIOUS */
  Rotation _rotation_absolute;
  Rotation _rotation_relative; /* wrt PREVIOUS */
  int      _rotation_is_identity;
  _class_Tunneling_sample_parameters _parameters;
};
typedef struct _struct_Tunneling_sample _class_Tunneling_sample;
_class_Tunneling_sample _sample_var;
_class_Tunneling_sample *_sample = &_sample_var;
#pragma acc declare create ( _sample_var )
#pragma acc declare create ( _sample )

/* component detectorarm=Arm() [46] DECLARE */
/* Parameter definition for component type 'Arm' */
struct _struct_Arm_parameters {
  char Arm_has_no_parameters;
}; /* _struct_Arm_parameters */
typedef struct _struct_Arm_parameters _class_Arm_parameters;

/* Parameters for component type 'Arm' */
struct _struct_Arm {
  char     _name[256]; /* e.g. detectorarm */
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
_class_Arm _detectorarm_var;
_class_Arm *_detectorarm = &_detectorarm_var;
#pragma acc declare create ( _detectorarm_var )
#pragma acc declare create ( _detectorarm )

_class_TOF_monitor _TOFdetector_var;
_class_TOF_monitor *_TOFdetector = &_TOFdetector_var;
#pragma acc declare create ( _TOFdetector_var )
#pragma acc declare create ( _TOFdetector )

_class_TOF_monitor _TOFdetector_zoom_var;
_class_TOF_monitor *_TOFdetector_zoom = &_TOFdetector_zoom_var;
#pragma acc declare create ( _TOFdetector_zoom_var )
#pragma acc declare create ( _TOFdetector_zoom )

_class_E_monitor _Edetector_var;
_class_E_monitor *_Edetector = &_Edetector_var;
#pragma acc declare create ( _Edetector_var )
#pragma acc declare create ( _Edetector )

/* component TOF2Edetector=TOF2E_monitor() [50] DECLARE */
/* Parameter definition for component type 'TOF2E_monitor' */
struct _struct_TOF2E_monitor_parameters {
  /* Component type 'TOF2E_monitor' setting parameters */
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
  MCNUM T_zero;
  MCNUM L_flight;
  MCNUM restore_neutron;
  /* Component type 'TOF2E_monitor' private parameters */
  /* Component type 'TOF2E_monitor' DECLARE code stored as structure members */
  DArray1d E_N;
  DArray1d E_p;
  DArray1d E_p2;
  double S_p, S_pE, S_pE2;
}; /* _struct_TOF2E_monitor_parameters */
typedef struct _struct_TOF2E_monitor_parameters _class_TOF2E_monitor_parameters;

/* Parameters for component type 'TOF2E_monitor' */
struct _struct_TOF2E_monitor {
  char     _name[256]; /* e.g. TOF2Edetector */
  char     _type[256]; /* TOF2E_monitor */
  long     _index; /* e.g. 2 index in TRACE list */
  Coords   _position_absolute;
  Coords   _position_relative; /* wrt PREVIOUS */
  Rotation _rotation_absolute;
  Rotation _rotation_relative; /* wrt PREVIOUS */
  int      _rotation_is_identity;
  _class_TOF2E_monitor_parameters _parameters;
};
typedef struct _struct_TOF2E_monitor _class_TOF2E_monitor;
_class_TOF2E_monitor _TOF2Edetector_var;
_class_TOF2E_monitor *_TOF2Edetector = &_TOF2Edetector_var;
#pragma acc declare create ( _TOF2Edetector_var )
#pragma acc declare create ( _TOF2Edetector )

int mcNUMCOMP = 50;

/* User declarations from instrument definition. Can define functions. */
        double FREQ;
        double t_FO1, t_FO2, t_fast1, t_fast2, t_fast2a, t_fast3, t_sample, t_detector;
        double tmin_zoom, tmax_zoom, t_offset;
        double E_target;

#undef compcurname
#undef compcurtype
#undef compcurindex
/* end of instrument 'ESS_IN5_reprate' and components DECLARE */

/* *****************************************************************************
* instrument 'ESS_IN5_reprate' and components INITIALISE
***************************************************************************** */

/* component source=ESS_moderator_long() SETTING, POSITION/ROTATION */
int _source_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_source_setpos] component source=ESS_moderator_long() SETTING [/usr/share/mcstas/3.0-dev/obsolete/ESS_moderator_long.comp:166]");
  stracpy(_source->_name, "source", 16384);
  stracpy(_source->_type, "ESS_moderator_long", 16384);
  _source->_index=1;
  _source->_parameters.width_c = 0;
  #define width_c (_source->_parameters.width_c)
  _source->_parameters.yheight = 0.12;
  #define yheight (_source->_parameters.yheight)
  _source->_parameters.Lmin = instrument->_parameters._Lmin;
  #define Lmin (_source->_parameters.Lmin)
  _source->_parameters.Lmax = instrument->_parameters._Lmax;
  #define Lmax (_source->_parameters.Lmax)
  _source->_parameters.dist = instrument->_parameters._GUI_start;
  #define dist (_source->_parameters.dist)
  _source->_parameters.focus_xw = instrument->_parameters._GUI_w;
  #define focus_xw (_source->_parameters.focus_xw)
  _source->_parameters.focus_yh = instrument->_parameters._GUI_h;
  #define focus_yh (_source->_parameters.focus_yh)
  _source->_parameters.nu = FREQ;
  #define nu (_source->_parameters.nu)
  _source->_parameters.T = 50;
  #define T (_source->_parameters.T)
  _source->_parameters.tau = 287e-6;
  #define tau (_source->_parameters.tau)
  _source->_parameters.tau1 = 0;
  #define tau1 (_source->_parameters.tau1)
  _source->_parameters.tau2 = 20e-6;
  #define tau2 (_source->_parameters.tau2)
  _source->_parameters.d = instrument->_parameters._Pulse_width;
  #define d (_source->_parameters.d)
  _source->_parameters.n = 20;
  #define n (_source->_parameters.n)
  _source->_parameters.cold_frac = instrument->_parameters._cold;
  #define cold_frac (_source->_parameters.cold_frac)
  _source->_parameters.n2 = 5;
  #define n2 (_source->_parameters.n2)
  _source->_parameters.chi2 = 0.9;
  #define chi2 (_source->_parameters.chi2)
  _source->_parameters.I0 = 6.9e11;
  #define I0 (_source->_parameters.I0)
  _source->_parameters.I2 = 27.6e10;
  #define I2 (_source->_parameters.I2)
  _source->_parameters.target_index = 0;
  #define target_index (_source->_parameters.target_index)
  _source->_parameters.cyl_radius = 0.085;
  #define cyl_radius (_source->_parameters.cyl_radius)
  _source->_parameters.branch1 = 1;
  #define branch1 (_source->_parameters.branch1)
  _source->_parameters.branch2 = 0.5;
  #define branch2 (_source->_parameters.branch2)
  _source->_parameters.branch_tail = 0.1;
  #define branch_tail (_source->_parameters.branch_tail)
  _source->_parameters.n_pulses = instrument->_parameters._Num_pulses;
  #define n_pulses (_source->_parameters.n_pulses)
  _source->_parameters.width_t = 0.12;
  #define width_t (_source->_parameters.width_t)
  _source->_parameters.T_t = 325;
  #define T_t (_source->_parameters.T_t)
  _source->_parameters.tau_t = 80e-6;
  #define tau_t (_source->_parameters.tau_t)
  _source->_parameters.tau1_t = 400e-6;
  #define tau1_t (_source->_parameters.tau1_t)
  _source->_parameters.tau2_t = 12e-6;
  #define tau2_t (_source->_parameters.tau2_t)
  _source->_parameters.chi2_t = 2.5;
  #define chi2_t (_source->_parameters.chi2_t)
  _source->_parameters.I0_t = 13.5e11;
  #define I0_t (_source->_parameters.I0_t)
  _source->_parameters.I2_t = 27.6e10;
  #define I2_t (_source->_parameters.I2_t)
  _source->_parameters.branch1_t = 0.5;
  #define branch1_t (_source->_parameters.branch1_t)
  _source->_parameters.branch2_t = 0.5;
  #define branch2_t (_source->_parameters.branch2_t)
  _source->_parameters.src_2012 = 1;
  #define src_2012 (_source->_parameters.src_2012)
  _source->_parameters.tfocus_dist = 0.1;
  #define tfocus_dist (_source->_parameters.tfocus_dist)
  _source->_parameters.tfocus_time = 0.0;
  #define tfocus_time (_source->_parameters.tfocus_time)
  _source->_parameters.tfocus_width = 0.0;
  #define tfocus_width (_source->_parameters.tfocus_width)
  _source->_parameters.beamport_angle = instrument->_parameters._port;
  #define beamport_angle (_source->_parameters.beamport_angle)

  #define l_range (_source->_parameters.l_range)
  #define w_mult (_source->_parameters.w_mult)
  #define w_geom (_source->_parameters.w_geom)
  #define w_geom_c (_source->_parameters.w_geom_c)
  #define w_geom_t (_source->_parameters.w_geom_t)
  #define tx (_source->_parameters.tx)
  #define ty (_source->_parameters.ty)
  #define tz (_source->_parameters.tz)
  #define t1x (_source->_parameters.t1x)
  #define t1y (_source->_parameters.t1y)
  #define t1z (_source->_parameters.t1z)
  #define t2x (_source->_parameters.t2x)
  #define t2y (_source->_parameters.t2y)
  #define t2z (_source->_parameters.t2z)
  #define T_n (_source->_parameters.T_n)
  #define tau_n (_source->_parameters.tau_n)
  #define tau1_n (_source->_parameters.tau1_n)
  #define tau2_n (_source->_parameters.tau2_n)
  #define chi2_n (_source->_parameters.chi2_n)
  #define I0_n (_source->_parameters.I0_n)
  #define I2_n (_source->_parameters.I2_n)
  #define branch1_n (_source->_parameters.branch1_n)
  #define branch2_n (_source->_parameters.branch2_n)
  #define r_empty (_source->_parameters.r_empty)
  #define r_optics (_source->_parameters.r_optics)

  #undef width_c
  #undef yheight
  #undef Lmin
  #undef Lmax
  #undef dist
  #undef focus_xw
  #undef focus_yh
  #undef nu
  #undef T
  #undef tau
  #undef tau1
  #undef tau2
  #undef d
  #undef n
  #undef cold_frac
  #undef n2
  #undef chi2
  #undef I0
  #undef I2
  #undef target_index
  #undef cyl_radius
  #undef branch1
  #undef branch2
  #undef branch_tail
  #undef n_pulses
  #undef width_t
  #undef T_t
  #undef tau_t
  #undef tau1_t
  #undef tau2_t
  #undef chi2_t
  #undef I0_t
  #undef I2_t
  #undef branch1_t
  #undef branch2_t
  #undef src_2012
  #undef tfocus_dist
  #undef tfocus_time
  #undef tfocus_width
  #undef beamport_angle
  #undef l_range
  #undef w_mult
  #undef w_geom
  #undef w_geom_c
  #undef w_geom_t
  #undef tx
  #undef ty
  #undef tz
  #undef t1x
  #undef t1y
  #undef t1z
  #undef t2x
  #undef t2y
  #undef t2z
  #undef T_n
  #undef tau_n
  #undef tau1_n
  #undef tau2_n
  #undef chi2_n
  #undef I0_n
  #undef I2_n
  #undef branch1_n
  #undef branch2_n
  #undef r_empty
  #undef r_optics
  /* component source=ESS_moderator_long() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(_source->_rotation_absolute,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_copy(_source->_rotation_relative, _source->_rotation_absolute);
    _source->_rotation_is_identity =  rot_test_identity(_source->_rotation_relative);
    _source->_position_absolute = coords_set(
      0, 0, -0.1);
    tc1 = coords_neg(_source->_position_absolute);
    _source->_position_relative = rot_apply(_source->_rotation_absolute, tc1);
  } /* source=ESS_moderator_long() AT ROTATED */
  DEBUG_COMPONENT("source", _source->_position_absolute, _source->_rotation_absolute);
  instrument->_position_absolute[1] = _source->_position_absolute;
  instrument->_position_relative[1] = _source->_position_relative;
  instrument->counter_N[1]  = instrument->counter_P[1] = instrument->counter_P2[1] = 0;
  instrument->counter_AbsorbProp[1]= 0;
  return(0);
} /* _source_setpos */

/* component Origin=Progress_bar() SETTING, POSITION/ROTATION */
int _Origin_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Origin_setpos] component Origin=Progress_bar() SETTING [/usr/share/mcstas/3.0-dev/misc/Progress_bar.comp:57]");
  stracpy(_Origin->_name, "Origin", 16384);
  stracpy(_Origin->_type, "Progress_bar", 16384);
  _Origin->_index=2;
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
    rot_transpose(_source->_rotation_absolute, tr1);
    rot_mul(_Origin->_rotation_absolute, tr1, _Origin->_rotation_relative);
    _Origin->_rotation_is_identity =  rot_test_identity(_Origin->_rotation_relative);
    _Origin->_position_absolute = coords_set(
      0, 0, 0);
    tc1 = coords_sub(_source->_position_absolute, _Origin->_position_absolute);
    _Origin->_position_relative = rot_apply(_Origin->_rotation_absolute, tc1);
  } /* Origin=Progress_bar() AT ROTATED */
  DEBUG_COMPONENT("Origin", _Origin->_position_absolute, _Origin->_rotation_absolute);
  instrument->_position_absolute[2] = _Origin->_position_absolute;
  instrument->_position_relative[2] = _Origin->_position_relative;
  instrument->counter_N[2]  = instrument->counter_P[2] = instrument->counter_P2[2] = 0;
  instrument->counter_AbsorbProp[2]= 0;
  return(0);
} /* _Origin_setpos */

/* component TOFmoderator_zoom=TOF_monitor() SETTING, POSITION/ROTATION */
int _TOFmoderator_zoom_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_TOFmoderator_zoom_setpos] component TOFmoderator_zoom=TOF_monitor() SETTING [/usr/share/mcstas/3.0-dev/monitors/TOF_monitor.comp:63]");
  stracpy(_TOFmoderator_zoom->_name, "TOFmoderator_zoom", 16384);
  stracpy(_TOFmoderator_zoom->_type, "TOF_monitor", 16384);
  _TOFmoderator_zoom->_index=3;
  _TOFmoderator_zoom->_parameters.nt = 1000;
  #define nt (_TOFmoderator_zoom->_parameters.nt)
  if("TOFmoderator_zoom.dat" && strlen("TOFmoderator_zoom.dat"))
    stracpy(_TOFmoderator_zoom->_parameters.filename, "TOFmoderator_zoom.dat" ? "TOFmoderator_zoom.dat" : "", 16384);
  else 
  _TOFmoderator_zoom->_parameters.filename[0]='\0';
  #define filename (_TOFmoderator_zoom->_parameters.filename)
  _TOFmoderator_zoom->_parameters.xmin = -0.05;
  #define xmin (_TOFmoderator_zoom->_parameters.xmin)
  _TOFmoderator_zoom->_parameters.xmax = 0.05;
  #define xmax (_TOFmoderator_zoom->_parameters.xmax)
  _TOFmoderator_zoom->_parameters.ymin = -0.05;
  #define ymin (_TOFmoderator_zoom->_parameters.ymin)
  _TOFmoderator_zoom->_parameters.ymax = 0.05;
  #define ymax (_TOFmoderator_zoom->_parameters.ymax)
  _TOFmoderator_zoom->_parameters.xwidth = 0.12;
  #define xwidth (_TOFmoderator_zoom->_parameters.xwidth)
  _TOFmoderator_zoom->_parameters.yheight = 0.12;
  #define yheight (_TOFmoderator_zoom->_parameters.yheight)
  _TOFmoderator_zoom->_parameters.tmin = 0;
  #define tmin (_TOFmoderator_zoom->_parameters.tmin)
  _TOFmoderator_zoom->_parameters.tmax = 5000;
  #define tmax (_TOFmoderator_zoom->_parameters.tmax)
  _TOFmoderator_zoom->_parameters.dt = 1.0;
  #define dt (_TOFmoderator_zoom->_parameters.dt)
  _TOFmoderator_zoom->_parameters.restore_neutron = 1;
  #define restore_neutron (_TOFmoderator_zoom->_parameters.restore_neutron)

  #define TOF_N (_TOFmoderator_zoom->_parameters.TOF_N)
  #define TOF_p (_TOFmoderator_zoom->_parameters.TOF_p)
  #define TOF_p2 (_TOFmoderator_zoom->_parameters.TOF_p2)
  #define t_min (_TOFmoderator_zoom->_parameters.t_min)
  #define t_max (_TOFmoderator_zoom->_parameters.t_max)
  #define delta_t (_TOFmoderator_zoom->_parameters.delta_t)

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
  /* component TOFmoderator_zoom=TOF_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Origin->_rotation_absolute, _TOFmoderator_zoom->_rotation_absolute);
    rot_transpose(_Origin->_rotation_absolute, tr1);
    rot_mul(_TOFmoderator_zoom->_rotation_absolute, tr1, _TOFmoderator_zoom->_rotation_relative);
    _TOFmoderator_zoom->_rotation_is_identity =  rot_test_identity(_TOFmoderator_zoom->_rotation_relative);
    tc1 = coords_set(
      0, 0, 1E-6);
    rot_transpose(_Origin->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _TOFmoderator_zoom->_position_absolute = coords_add(_Origin->_position_absolute, tc2);
    tc1 = coords_sub(_Origin->_position_absolute, _TOFmoderator_zoom->_position_absolute);
    _TOFmoderator_zoom->_position_relative = rot_apply(_TOFmoderator_zoom->_rotation_absolute, tc1);
  } /* TOFmoderator_zoom=TOF_monitor() AT ROTATED */
  DEBUG_COMPONENT("TOFmoderator_zoom", _TOFmoderator_zoom->_position_absolute, _TOFmoderator_zoom->_rotation_absolute);
  instrument->_position_absolute[3] = _TOFmoderator_zoom->_position_absolute;
  instrument->_position_relative[3] = _TOFmoderator_zoom->_position_relative;
  instrument->counter_N[3]  = instrument->counter_P[3] = instrument->counter_P2[3] = 0;
  instrument->counter_AbsorbProp[3]= 0;
  return(0);
} /* _TOFmoderator_zoom_setpos */

/* component TOFmoderator=TOF_monitor() SETTING, POSITION/ROTATION */
int _TOFmoderator_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_TOFmoderator_setpos] component TOFmoderator=TOF_monitor() SETTING [/usr/share/mcstas/3.0-dev/monitors/TOF_monitor.comp:63]");
  stracpy(_TOFmoderator->_name, "TOFmoderator", 16384);
  stracpy(_TOFmoderator->_type, "TOF_monitor", 16384);
  _TOFmoderator->_index=4;
  _TOFmoderator->_parameters.nt = 1000;
  #define nt (_TOFmoderator->_parameters.nt)
  if("TOFmoderator.dat" && strlen("TOFmoderator.dat"))
    stracpy(_TOFmoderator->_parameters.filename, "TOFmoderator.dat" ? "TOFmoderator.dat" : "", 16384);
  else 
  _TOFmoderator->_parameters.filename[0]='\0';
  #define filename (_TOFmoderator->_parameters.filename)
  _TOFmoderator->_parameters.xmin = -0.05;
  #define xmin (_TOFmoderator->_parameters.xmin)
  _TOFmoderator->_parameters.xmax = 0.05;
  #define xmax (_TOFmoderator->_parameters.xmax)
  _TOFmoderator->_parameters.ymin = -0.05;
  #define ymin (_TOFmoderator->_parameters.ymin)
  _TOFmoderator->_parameters.ymax = 0.05;
  #define ymax (_TOFmoderator->_parameters.ymax)
  _TOFmoderator->_parameters.xwidth = 0.12;
  #define xwidth (_TOFmoderator->_parameters.xwidth)
  _TOFmoderator->_parameters.yheight = 0.12;
  #define yheight (_TOFmoderator->_parameters.yheight)
  _TOFmoderator->_parameters.tmin = 0;
  #define tmin (_TOFmoderator->_parameters.tmin)
  _TOFmoderator->_parameters.tmax = 3.0e5;
  #define tmax (_TOFmoderator->_parameters.tmax)
  _TOFmoderator->_parameters.dt = 1.0;
  #define dt (_TOFmoderator->_parameters.dt)
  _TOFmoderator->_parameters.restore_neutron = 1;
  #define restore_neutron (_TOFmoderator->_parameters.restore_neutron)

  #define TOF_N (_TOFmoderator->_parameters.TOF_N)
  #define TOF_p (_TOFmoderator->_parameters.TOF_p)
  #define TOF_p2 (_TOFmoderator->_parameters.TOF_p2)
  #define t_min (_TOFmoderator->_parameters.t_min)
  #define t_max (_TOFmoderator->_parameters.t_max)
  #define delta_t (_TOFmoderator->_parameters.delta_t)

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
  /* component TOFmoderator=TOF_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Origin->_rotation_absolute, _TOFmoderator->_rotation_absolute);
    rot_transpose(_TOFmoderator_zoom->_rotation_absolute, tr1);
    rot_mul(_TOFmoderator->_rotation_absolute, tr1, _TOFmoderator->_rotation_relative);
    _TOFmoderator->_rotation_is_identity =  rot_test_identity(_TOFmoderator->_rotation_relative);
    tc1 = coords_set(
      0, 0, 1E-6);
    rot_transpose(_Origin->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _TOFmoderator->_position_absolute = coords_add(_Origin->_position_absolute, tc2);
    tc1 = coords_sub(_TOFmoderator_zoom->_position_absolute, _TOFmoderator->_position_absolute);
    _TOFmoderator->_position_relative = rot_apply(_TOFmoderator->_rotation_absolute, tc1);
  } /* TOFmoderator=TOF_monitor() AT ROTATED */
  DEBUG_COMPONENT("TOFmoderator", _TOFmoderator->_position_absolute, _TOFmoderator->_rotation_absolute);
  instrument->_position_absolute[4] = _TOFmoderator->_position_absolute;
  instrument->_position_relative[4] = _TOFmoderator->_position_relative;
  instrument->counter_N[4]  = instrument->counter_P[4] = instrument->counter_P2[4] = 0;
  instrument->counter_AbsorbProp[4]= 0;
  return(0);
} /* _TOFmoderator_setpos */

/* component Lmon_guistart=L_monitor() SETTING, POSITION/ROTATION */
int _Lmon_guistart_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Lmon_guistart_setpos] component Lmon_guistart=L_monitor() SETTING [/usr/share/mcstas/3.0-dev/monitors/L_monitor.comp:65]");
  stracpy(_Lmon_guistart->_name, "Lmon_guistart", 16384);
  stracpy(_Lmon_guistart->_type, "L_monitor", 16384);
  _Lmon_guistart->_index=5;
  _Lmon_guistart->_parameters.nL = 1000;
  #define nL (_Lmon_guistart->_parameters.nL)
  if("Lmon_guistart.dat" && strlen("Lmon_guistart.dat"))
    stracpy(_Lmon_guistart->_parameters.filename, "Lmon_guistart.dat" ? "Lmon_guistart.dat" : "", 16384);
  else 
  _Lmon_guistart->_parameters.filename[0]='\0';
  #define filename (_Lmon_guistart->_parameters.filename)
  _Lmon_guistart->_parameters.xmin = -0.05;
  #define xmin (_Lmon_guistart->_parameters.xmin)
  _Lmon_guistart->_parameters.xmax = 0.05;
  #define xmax (_Lmon_guistart->_parameters.xmax)
  _Lmon_guistart->_parameters.ymin = -0.05;
  #define ymin (_Lmon_guistart->_parameters.ymin)
  _Lmon_guistart->_parameters.ymax = 0.05;
  #define ymax (_Lmon_guistart->_parameters.ymax)
  _Lmon_guistart->_parameters.xwidth = instrument->_parameters._GUI_w + 0.01;
  #define xwidth (_Lmon_guistart->_parameters.xwidth)
  _Lmon_guistart->_parameters.yheight = instrument->_parameters._GUI_h + 0.01;
  #define yheight (_Lmon_guistart->_parameters.yheight)
  _Lmon_guistart->_parameters.Lmin = 0;
  #define Lmin (_Lmon_guistart->_parameters.Lmin)
  _Lmon_guistart->_parameters.Lmax = instrument->_parameters._Lmax + 1;
  #define Lmax (_Lmon_guistart->_parameters.Lmax)
  _Lmon_guistart->_parameters.restore_neutron = 0;
  #define restore_neutron (_Lmon_guistart->_parameters.restore_neutron)

  #define L_N (_Lmon_guistart->_parameters.L_N)
  #define L_p (_Lmon_guistart->_parameters.L_p)
  #define L_p2 (_Lmon_guistart->_parameters.L_p2)

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
  /* component Lmon_guistart=L_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Origin->_rotation_absolute, _Lmon_guistart->_rotation_absolute);
    rot_transpose(_TOFmoderator->_rotation_absolute, tr1);
    rot_mul(_Lmon_guistart->_rotation_absolute, tr1, _Lmon_guistart->_rotation_relative);
    _Lmon_guistart->_rotation_is_identity =  rot_test_identity(_Lmon_guistart->_rotation_relative);
    tc1 = coords_set(
      0, 0, instrument->_parameters._GUI_start -2e-6);
    rot_transpose(_Origin->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Lmon_guistart->_position_absolute = coords_add(_Origin->_position_absolute, tc2);
    tc1 = coords_sub(_TOFmoderator->_position_absolute, _Lmon_guistart->_position_absolute);
    _Lmon_guistart->_position_relative = rot_apply(_Lmon_guistart->_rotation_absolute, tc1);
  } /* Lmon_guistart=L_monitor() AT ROTATED */
  DEBUG_COMPONENT("Lmon_guistart", _Lmon_guistart->_position_absolute, _Lmon_guistart->_rotation_absolute);
  instrument->_position_absolute[5] = _Lmon_guistart->_position_absolute;
  instrument->_position_relative[5] = _Lmon_guistart->_position_relative;
  instrument->counter_N[5]  = instrument->counter_P[5] = instrument->counter_P2[5] = 0;
  instrument->counter_AbsorbProp[5]= 0;
  return(0);
} /* _Lmon_guistart_setpos */

/* component Lmon_normalize=L_monitor() SETTING, POSITION/ROTATION */
int _Lmon_normalize_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Lmon_normalize_setpos] component Lmon_normalize=L_monitor() SETTING [/usr/share/mcstas/3.0-dev/monitors/L_monitor.comp:65]");
  stracpy(_Lmon_normalize->_name, "Lmon_normalize", 16384);
  stracpy(_Lmon_normalize->_type, "L_monitor", 16384);
  _Lmon_normalize->_index=6;
  _Lmon_normalize->_parameters.nL = 2880;
  #define nL (_Lmon_normalize->_parameters.nL)
  if("Lmon_guistart_normalize.dat" && strlen("Lmon_guistart_normalize.dat"))
    stracpy(_Lmon_normalize->_parameters.filename, "Lmon_guistart_normalize.dat" ? "Lmon_guistart_normalize.dat" : "", 16384);
  else 
  _Lmon_normalize->_parameters.filename[0]='\0';
  #define filename (_Lmon_normalize->_parameters.filename)
  _Lmon_normalize->_parameters.xmin = -0.05;
  #define xmin (_Lmon_normalize->_parameters.xmin)
  _Lmon_normalize->_parameters.xmax = 0.05;
  #define xmax (_Lmon_normalize->_parameters.xmax)
  _Lmon_normalize->_parameters.ymin = -0.05;
  #define ymin (_Lmon_normalize->_parameters.ymin)
  _Lmon_normalize->_parameters.ymax = 0.05;
  #define ymax (_Lmon_normalize->_parameters.ymax)
  _Lmon_normalize->_parameters.xwidth = 0.10;
  #define xwidth (_Lmon_normalize->_parameters.xwidth)
  _Lmon_normalize->_parameters.yheight = 0.10;
  #define yheight (_Lmon_normalize->_parameters.yheight)
  _Lmon_normalize->_parameters.Lmin = 0;
  #define Lmin (_Lmon_normalize->_parameters.Lmin)
  _Lmon_normalize->_parameters.Lmax = 20;
  #define Lmax (_Lmon_normalize->_parameters.Lmax)
  _Lmon_normalize->_parameters.restore_neutron = 0;
  #define restore_neutron (_Lmon_normalize->_parameters.restore_neutron)

  #define L_N (_Lmon_normalize->_parameters.L_N)
  #define L_p (_Lmon_normalize->_parameters.L_p)
  #define L_p2 (_Lmon_normalize->_parameters.L_p2)

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
  /* component Lmon_normalize=L_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Origin->_rotation_absolute, _Lmon_normalize->_rotation_absolute);
    rot_transpose(_Lmon_guistart->_rotation_absolute, tr1);
    rot_mul(_Lmon_normalize->_rotation_absolute, tr1, _Lmon_normalize->_rotation_relative);
    _Lmon_normalize->_rotation_is_identity =  rot_test_identity(_Lmon_normalize->_rotation_relative);
    tc1 = coords_set(
      0, 0, instrument->_parameters._GUI_start -1e-6);
    rot_transpose(_Origin->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Lmon_normalize->_position_absolute = coords_add(_Origin->_position_absolute, tc2);
    tc1 = coords_sub(_Lmon_guistart->_position_absolute, _Lmon_normalize->_position_absolute);
    _Lmon_normalize->_position_relative = rot_apply(_Lmon_normalize->_rotation_absolute, tc1);
  } /* Lmon_normalize=L_monitor() AT ROTATED */
  DEBUG_COMPONENT("Lmon_normalize", _Lmon_normalize->_position_absolute, _Lmon_normalize->_rotation_absolute);
  instrument->_position_absolute[6] = _Lmon_normalize->_position_absolute;
  instrument->_position_relative[6] = _Lmon_normalize->_position_relative;
  instrument->counter_N[6]  = instrument->counter_P[6] = instrument->counter_P2[6] = 0;
  instrument->counter_AbsorbProp[6]= 0;
  return(0);
} /* _Lmon_normalize_setpos */

/* component Guide1=Guide() SETTING, POSITION/ROTATION */
int _Guide1_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Guide1_setpos] component Guide1=Guide() SETTING [/usr/share/mcstas/3.0-dev/optics/Guide.comp:73]");
  stracpy(_Guide1->_name, "Guide1", 16384);
  stracpy(_Guide1->_type, "Guide", 16384);
  _Guide1->_index=7;
  _Guide1->_parameters.reflect[0]='\0';
  #define reflect (_Guide1->_parameters.reflect)
  _Guide1->_parameters.w1 = instrument->_parameters._GUI_w;
  #define w1 (_Guide1->_parameters.w1)
  _Guide1->_parameters.h1 = instrument->_parameters._GUI_h;
  #define h1 (_Guide1->_parameters.h1)
  _Guide1->_parameters.w2 = instrument->_parameters._W1;
  #define w2 (_Guide1->_parameters.w2)
  _Guide1->_parameters.h2 = instrument->_parameters._H1;
  #define h2 (_Guide1->_parameters.h2)
  _Guide1->_parameters.l = instrument->_parameters._FO1_DIST - instrument->_parameters._GUI_start - instrument->_parameters._GUI_GAP / 2;
  #define l (_Guide1->_parameters.l)
  _Guide1->_parameters.R0 = 1;
  #define R0 (_Guide1->_parameters.R0)
  _Guide1->_parameters.Qc = 0.0219;
  #define Qc (_Guide1->_parameters.Qc)
  _Guide1->_parameters.alpha = instrument->_parameters._ALPHA;
  #define alpha (_Guide1->_parameters.alpha)
  _Guide1->_parameters.m = instrument->_parameters._M;
  #define m (_Guide1->_parameters.m)
  _Guide1->_parameters.W = 0.003;
  #define W (_Guide1->_parameters.W)

  #define pTable (_Guide1->_parameters.pTable)

  #undef reflect
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
  #undef pTable
  /* component Guide1=Guide() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Origin->_rotation_absolute, _Guide1->_rotation_absolute);
    rot_transpose(_Lmon_normalize->_rotation_absolute, tr1);
    rot_mul(_Guide1->_rotation_absolute, tr1, _Guide1->_rotation_relative);
    _Guide1->_rotation_is_identity =  rot_test_identity(_Guide1->_rotation_relative);
    tc1 = coords_set(
      0, 0, instrument->_parameters._GUI_start);
    rot_transpose(_Origin->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Guide1->_position_absolute = coords_add(_Origin->_position_absolute, tc2);
    tc1 = coords_sub(_Lmon_normalize->_position_absolute, _Guide1->_position_absolute);
    _Guide1->_position_relative = rot_apply(_Guide1->_rotation_absolute, tc1);
  } /* Guide1=Guide() AT ROTATED */
  DEBUG_COMPONENT("Guide1", _Guide1->_position_absolute, _Guide1->_rotation_absolute);
  instrument->_position_absolute[7] = _Guide1->_position_absolute;
  instrument->_position_relative[7] = _Guide1->_position_relative;
  instrument->counter_N[7]  = instrument->counter_P[7] = instrument->counter_P2[7] = 0;
  instrument->counter_AbsorbProp[7]= 0;
  return(0);
} /* _Guide1_setpos */

/* component Lmonslow1=L_monitor() SETTING, POSITION/ROTATION */
int _Lmonslow1_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Lmonslow1_setpos] component Lmonslow1=L_monitor() SETTING [/usr/share/mcstas/3.0-dev/monitors/L_monitor.comp:65]");
  stracpy(_Lmonslow1->_name, "Lmonslow1", 16384);
  stracpy(_Lmonslow1->_type, "L_monitor", 16384);
  _Lmonslow1->_index=8;
  _Lmonslow1->_parameters.nL = 200;
  #define nL (_Lmonslow1->_parameters.nL)
  if("Lmonslow1.dat" && strlen("Lmonslow1.dat"))
    stracpy(_Lmonslow1->_parameters.filename, "Lmonslow1.dat" ? "Lmonslow1.dat" : "", 16384);
  else 
  _Lmonslow1->_parameters.filename[0]='\0';
  #define filename (_Lmonslow1->_parameters.filename)
  _Lmonslow1->_parameters.xmin = -0.05;
  #define xmin (_Lmonslow1->_parameters.xmin)
  _Lmonslow1->_parameters.xmax = 0.05;
  #define xmax (_Lmonslow1->_parameters.xmax)
  _Lmonslow1->_parameters.ymin = -0.05;
  #define ymin (_Lmonslow1->_parameters.ymin)
  _Lmonslow1->_parameters.ymax = 0.05;
  #define ymax (_Lmonslow1->_parameters.ymax)
  _Lmonslow1->_parameters.xwidth = 0.06;
  #define xwidth (_Lmonslow1->_parameters.xwidth)
  _Lmonslow1->_parameters.yheight = 0.21;
  #define yheight (_Lmonslow1->_parameters.yheight)
  _Lmonslow1->_parameters.Lmin = 0;
  #define Lmin (_Lmonslow1->_parameters.Lmin)
  _Lmonslow1->_parameters.Lmax = instrument->_parameters._Lmax + 1;
  #define Lmax (_Lmonslow1->_parameters.Lmax)
  _Lmonslow1->_parameters.restore_neutron = 0;
  #define restore_neutron (_Lmonslow1->_parameters.restore_neutron)

  #define L_N (_Lmonslow1->_parameters.L_N)
  #define L_p (_Lmonslow1->_parameters.L_p)
  #define L_p2 (_Lmonslow1->_parameters.L_p2)

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
  /* component Lmonslow1=L_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Origin->_rotation_absolute, _Lmonslow1->_rotation_absolute);
    rot_transpose(_Guide1->_rotation_absolute, tr1);
    rot_mul(_Lmonslow1->_rotation_absolute, tr1, _Lmonslow1->_rotation_relative);
    _Lmonslow1->_rotation_is_identity =  rot_test_identity(_Lmonslow1->_rotation_relative);
    tc1 = coords_set(
      0, 0, instrument->_parameters._FO1_DIST -2e-6);
    rot_transpose(_Origin->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Lmonslow1->_position_absolute = coords_add(_Origin->_position_absolute, tc2);
    tc1 = coords_sub(_Guide1->_position_absolute, _Lmonslow1->_position_absolute);
    _Lmonslow1->_position_relative = rot_apply(_Lmonslow1->_rotation_absolute, tc1);
  } /* Lmonslow1=L_monitor() AT ROTATED */
  DEBUG_COMPONENT("Lmonslow1", _Lmonslow1->_position_absolute, _Lmonslow1->_rotation_absolute);
  instrument->_position_absolute[8] = _Lmonslow1->_position_absolute;
  instrument->_position_relative[8] = _Lmonslow1->_position_relative;
  instrument->counter_N[8]  = instrument->counter_P[8] = instrument->counter_P2[8] = 0;
  instrument->counter_AbsorbProp[8]= 0;
  return(0);
} /* _Lmonslow1_setpos */

/* component PSDslow1=PSD_monitor() SETTING, POSITION/ROTATION */
int _PSDslow1_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_PSDslow1_setpos] component PSDslow1=PSD_monitor() SETTING [/usr/share/mcstas/3.0-dev/monitors/PSD_monitor.comp:64]");
  stracpy(_PSDslow1->_name, "PSDslow1", 16384);
  stracpy(_PSDslow1->_type, "PSD_monitor", 16384);
  _PSDslow1->_index=9;
  _PSDslow1->_parameters.nx = 90;
  #define nx (_PSDslow1->_parameters.nx)
  _PSDslow1->_parameters.ny = 90;
  #define ny (_PSDslow1->_parameters.ny)
  if("PSDslow1.dat" && strlen("PSDslow1.dat"))
    stracpy(_PSDslow1->_parameters.filename, "PSDslow1.dat" ? "PSDslow1.dat" : "", 16384);
  else 
  _PSDslow1->_parameters.filename[0]='\0';
  #define filename (_PSDslow1->_parameters.filename)
  _PSDslow1->_parameters.xmin = -0.05;
  #define xmin (_PSDslow1->_parameters.xmin)
  _PSDslow1->_parameters.xmax = 0.05;
  #define xmax (_PSDslow1->_parameters.xmax)
  _PSDslow1->_parameters.ymin = -0.05;
  #define ymin (_PSDslow1->_parameters.ymin)
  _PSDslow1->_parameters.ymax = 0.05;
  #define ymax (_PSDslow1->_parameters.ymax)
  _PSDslow1->_parameters.xwidth = 0.1;
  #define xwidth (_PSDslow1->_parameters.xwidth)
  _PSDslow1->_parameters.yheight = 0.25;
  #define yheight (_PSDslow1->_parameters.yheight)
  _PSDslow1->_parameters.restore_neutron = 0;
  #define restore_neutron (_PSDslow1->_parameters.restore_neutron)

  #define PSD_N (_PSDslow1->_parameters.PSD_N)
  #define PSD_p (_PSDslow1->_parameters.PSD_p)
  #define PSD_p2 (_PSDslow1->_parameters.PSD_p2)

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
  /* component PSDslow1=PSD_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Origin->_rotation_absolute, _PSDslow1->_rotation_absolute);
    rot_transpose(_Lmonslow1->_rotation_absolute, tr1);
    rot_mul(_PSDslow1->_rotation_absolute, tr1, _PSDslow1->_rotation_relative);
    _PSDslow1->_rotation_is_identity =  rot_test_identity(_PSDslow1->_rotation_relative);
    tc1 = coords_set(
      0, 0, instrument->_parameters._FO1_DIST -1e-6);
    rot_transpose(_Origin->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _PSDslow1->_position_absolute = coords_add(_Origin->_position_absolute, tc2);
    tc1 = coords_sub(_Lmonslow1->_position_absolute, _PSDslow1->_position_absolute);
    _PSDslow1->_position_relative = rot_apply(_PSDslow1->_rotation_absolute, tc1);
  } /* PSDslow1=PSD_monitor() AT ROTATED */
  DEBUG_COMPONENT("PSDslow1", _PSDslow1->_position_absolute, _PSDslow1->_rotation_absolute);
  instrument->_position_absolute[9] = _PSDslow1->_position_absolute;
  instrument->_position_relative[9] = _PSDslow1->_position_relative;
  instrument->counter_N[9]  = instrument->counter_P[9] = instrument->counter_P2[9] = 0;
  instrument->counter_AbsorbProp[9]= 0;
  return(0);
} /* _PSDslow1_setpos */

/* component FOchop1=DiskChopper() SETTING, POSITION/ROTATION */
int _FOchop1_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_FOchop1_setpos] component FOchop1=DiskChopper() SETTING [/usr/share/mcstas/3.0-dev/optics/DiskChopper.comp:67]");
  stracpy(_FOchop1->_name, "FOchop1", 16384);
  stracpy(_FOchop1->_type, "DiskChopper", 16384);
  _FOchop1->_index=10;
  _FOchop1->_parameters.theta_0 = instrument->_parameters._SLOW1_THETA;
  #define theta_0 (_FOchop1->_parameters.theta_0)
  _FOchop1->_parameters.radius = 0.6;
  #define radius (_FOchop1->_parameters.radius)
  _FOchop1->_parameters.yheight = 0.20;
  #define yheight (_FOchop1->_parameters.yheight)
  _FOchop1->_parameters.nu = instrument->_parameters._F_slow1;
  #define nu (_FOchop1->_parameters.nu)
  _FOchop1->_parameters.nslit = 1;
  #define nslit (_FOchop1->_parameters.nslit)
  _FOchop1->_parameters.jitter = 0;
  #define jitter (_FOchop1->_parameters.jitter)
  _FOchop1->_parameters.delay = t_FO1;
  #define delay (_FOchop1->_parameters.delay)
  _FOchop1->_parameters.isfirst = 0;
  #define isfirst (_FOchop1->_parameters.isfirst)
  _FOchop1->_parameters.n_pulse = 1;
  #define n_pulse (_FOchop1->_parameters.n_pulse)
  _FOchop1->_parameters.abs_out = 1;
  #define abs_out (_FOchop1->_parameters.abs_out)
  _FOchop1->_parameters.phase = 0;
  #define phase (_FOchop1->_parameters.phase)
  _FOchop1->_parameters.xwidth = 0;
  #define xwidth (_FOchop1->_parameters.xwidth)
  _FOchop1->_parameters.verbose = 0;
  #define verbose (_FOchop1->_parameters.verbose)

  #define Tg (_FOchop1->_parameters.Tg)
  #define To (_FOchop1->_parameters.To)
  #define delta_y (_FOchop1->_parameters.delta_y)
  #define height (_FOchop1->_parameters.height)
  #define omega (_FOchop1->_parameters.omega)

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
  /* component FOchop1=DiskChopper() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Origin->_rotation_absolute, _FOchop1->_rotation_absolute);
    rot_transpose(_PSDslow1->_rotation_absolute, tr1);
    rot_mul(_FOchop1->_rotation_absolute, tr1, _FOchop1->_rotation_relative);
    _FOchop1->_rotation_is_identity =  rot_test_identity(_FOchop1->_rotation_relative);
    tc1 = coords_set(
      0, 0, instrument->_parameters._FO1_DIST);
    rot_transpose(_Origin->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _FOchop1->_position_absolute = coords_add(_Origin->_position_absolute, tc2);
    tc1 = coords_sub(_PSDslow1->_position_absolute, _FOchop1->_position_absolute);
    _FOchop1->_position_relative = rot_apply(_FOchop1->_rotation_absolute, tc1);
  } /* FOchop1=DiskChopper() AT ROTATED */
  DEBUG_COMPONENT("FOchop1", _FOchop1->_position_absolute, _FOchop1->_rotation_absolute);
  instrument->_position_absolute[10] = _FOchop1->_position_absolute;
  instrument->_position_relative[10] = _FOchop1->_position_relative;
  instrument->counter_N[10]  = instrument->counter_P[10] = instrument->counter_P2[10] = 0;
  instrument->counter_AbsorbProp[10]= 0;
  return(0);
} /* _FOchop1_setpos */

/* component TOFLmon1=TOFLambda_monitor() SETTING, POSITION/ROTATION */
int _TOFLmon1_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_TOFLmon1_setpos] component TOFLmon1=TOFLambda_monitor() SETTING [/usr/share/mcstas/3.0-dev/monitors/TOFLambda_monitor.comp:68]");
  stracpy(_TOFLmon1->_name, "TOFLmon1", 16384);
  stracpy(_TOFLmon1->_type, "TOFLambda_monitor", 16384);
  _TOFLmon1->_index=11;
  _TOFLmon1->_parameters.nL = 200;
  #define nL (_TOFLmon1->_parameters.nL)
  _TOFLmon1->_parameters.nt = 200;
  #define nt (_TOFLmon1->_parameters.nt)
  _TOFLmon1->_parameters.tmin = 0;
  #define tmin (_TOFLmon1->_parameters.tmin)
  _TOFLmon1->_parameters.tmax = 3e5;
  #define tmax (_TOFLmon1->_parameters.tmax)
  if("TOFLmon1.dat" && strlen("TOFLmon1.dat"))
    stracpy(_TOFLmon1->_parameters.filename, "TOFLmon1.dat" ? "TOFLmon1.dat" : "", 16384);
  else 
  _TOFLmon1->_parameters.filename[0]='\0';
  #define filename (_TOFLmon1->_parameters.filename)
  _TOFLmon1->_parameters.xmin = -0.05;
  #define xmin (_TOFLmon1->_parameters.xmin)
  _TOFLmon1->_parameters.xmax = 0.05;
  #define xmax (_TOFLmon1->_parameters.xmax)
  _TOFLmon1->_parameters.ymin = -0.05;
  #define ymin (_TOFLmon1->_parameters.ymin)
  _TOFLmon1->_parameters.ymax = 0.05;
  #define ymax (_TOFLmon1->_parameters.ymax)
  _TOFLmon1->_parameters.xwidth = 0.05;
  #define xwidth (_TOFLmon1->_parameters.xwidth)
  _TOFLmon1->_parameters.yheight = 0.21;
  #define yheight (_TOFLmon1->_parameters.yheight)
  _TOFLmon1->_parameters.Lmin = 0;
  #define Lmin (_TOFLmon1->_parameters.Lmin)
  _TOFLmon1->_parameters.Lmax = instrument->_parameters._Lmax + 1;
  #define Lmax (_TOFLmon1->_parameters.Lmax)
  _TOFLmon1->_parameters.restore_neutron = 0;
  #define restore_neutron (_TOFLmon1->_parameters.restore_neutron)

  #define tt_0 (_TOFLmon1->_parameters.tt_0)
  #define tt_1 (_TOFLmon1->_parameters.tt_1)
  #define TOFL_N (_TOFLmon1->_parameters.TOFL_N)
  #define TOFL_p (_TOFLmon1->_parameters.TOFL_p)
  #define TOFL_p2 (_TOFLmon1->_parameters.TOFL_p2)

  #undef nL
  #undef nt
  #undef tmin
  #undef tmax
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
  #undef tt_0
  #undef tt_1
  #undef TOFL_N
  #undef TOFL_p
  #undef TOFL_p2
  /* component TOFLmon1=TOFLambda_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _FOchop1->_rotation_absolute, _TOFLmon1->_rotation_absolute);
    rot_transpose(_FOchop1->_rotation_absolute, tr1);
    rot_mul(_TOFLmon1->_rotation_absolute, tr1, _TOFLmon1->_rotation_relative);
    _TOFLmon1->_rotation_is_identity =  rot_test_identity(_TOFLmon1->_rotation_relative);
    tc1 = coords_set(
      0, 0, 1e-6);
    rot_transpose(_FOchop1->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _TOFLmon1->_position_absolute = coords_add(_FOchop1->_position_absolute, tc2);
    tc1 = coords_sub(_FOchop1->_position_absolute, _TOFLmon1->_position_absolute);
    _TOFLmon1->_position_relative = rot_apply(_TOFLmon1->_rotation_absolute, tc1);
  } /* TOFLmon1=TOFLambda_monitor() AT ROTATED */
  DEBUG_COMPONENT("TOFLmon1", _TOFLmon1->_position_absolute, _TOFLmon1->_rotation_absolute);
  instrument->_position_absolute[11] = _TOFLmon1->_position_absolute;
  instrument->_position_relative[11] = _TOFLmon1->_position_relative;
  instrument->counter_N[11]  = instrument->counter_P[11] = instrument->counter_P2[11] = 0;
  instrument->counter_AbsorbProp[11]= 0;
  return(0);
} /* _TOFLmon1_setpos */

/* component Lmon_afterslow1=L_monitor() SETTING, POSITION/ROTATION */
int _Lmon_afterslow1_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Lmon_afterslow1_setpos] component Lmon_afterslow1=L_monitor() SETTING [/usr/share/mcstas/3.0-dev/monitors/L_monitor.comp:65]");
  stracpy(_Lmon_afterslow1->_name, "Lmon_afterslow1", 16384);
  stracpy(_Lmon_afterslow1->_type, "L_monitor", 16384);
  _Lmon_afterslow1->_index=12;
  _Lmon_afterslow1->_parameters.nL = 200;
  #define nL (_Lmon_afterslow1->_parameters.nL)
  if("Lmon_afterslow1.dat" && strlen("Lmon_afterslow1.dat"))
    stracpy(_Lmon_afterslow1->_parameters.filename, "Lmon_afterslow1.dat" ? "Lmon_afterslow1.dat" : "", 16384);
  else 
  _Lmon_afterslow1->_parameters.filename[0]='\0';
  #define filename (_Lmon_afterslow1->_parameters.filename)
  _Lmon_afterslow1->_parameters.xmin = -0.05;
  #define xmin (_Lmon_afterslow1->_parameters.xmin)
  _Lmon_afterslow1->_parameters.xmax = 0.05;
  #define xmax (_Lmon_afterslow1->_parameters.xmax)
  _Lmon_afterslow1->_parameters.ymin = -0.05;
  #define ymin (_Lmon_afterslow1->_parameters.ymin)
  _Lmon_afterslow1->_parameters.ymax = 0.05;
  #define ymax (_Lmon_afterslow1->_parameters.ymax)
  _Lmon_afterslow1->_parameters.xwidth = 0.06;
  #define xwidth (_Lmon_afterslow1->_parameters.xwidth)
  _Lmon_afterslow1->_parameters.yheight = 0.21;
  #define yheight (_Lmon_afterslow1->_parameters.yheight)
  _Lmon_afterslow1->_parameters.Lmin = 0;
  #define Lmin (_Lmon_afterslow1->_parameters.Lmin)
  _Lmon_afterslow1->_parameters.Lmax = instrument->_parameters._Lmax + 1;
  #define Lmax (_Lmon_afterslow1->_parameters.Lmax)
  _Lmon_afterslow1->_parameters.restore_neutron = 0;
  #define restore_neutron (_Lmon_afterslow1->_parameters.restore_neutron)

  #define L_N (_Lmon_afterslow1->_parameters.L_N)
  #define L_p (_Lmon_afterslow1->_parameters.L_p)
  #define L_p2 (_Lmon_afterslow1->_parameters.L_p2)

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
  /* component Lmon_afterslow1=L_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _FOchop1->_rotation_absolute, _Lmon_afterslow1->_rotation_absolute);
    rot_transpose(_TOFLmon1->_rotation_absolute, tr1);
    rot_mul(_Lmon_afterslow1->_rotation_absolute, tr1, _Lmon_afterslow1->_rotation_relative);
    _Lmon_afterslow1->_rotation_is_identity =  rot_test_identity(_Lmon_afterslow1->_rotation_relative);
    tc1 = coords_set(
      0, 0, 2e-6);
    rot_transpose(_FOchop1->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Lmon_afterslow1->_position_absolute = coords_add(_FOchop1->_position_absolute, tc2);
    tc1 = coords_sub(_TOFLmon1->_position_absolute, _Lmon_afterslow1->_position_absolute);
    _Lmon_afterslow1->_position_relative = rot_apply(_Lmon_afterslow1->_rotation_absolute, tc1);
  } /* Lmon_afterslow1=L_monitor() AT ROTATED */
  DEBUG_COMPONENT("Lmon_afterslow1", _Lmon_afterslow1->_position_absolute, _Lmon_afterslow1->_rotation_absolute);
  instrument->_position_absolute[12] = _Lmon_afterslow1->_position_absolute;
  instrument->_position_relative[12] = _Lmon_afterslow1->_position_relative;
  instrument->counter_N[12]  = instrument->counter_P[12] = instrument->counter_P2[12] = 0;
  instrument->counter_AbsorbProp[12]= 0;
  return(0);
} /* _Lmon_afterslow1_setpos */

/* component PSD_afterslow1=PSD_monitor() SETTING, POSITION/ROTATION */
int _PSD_afterslow1_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_PSD_afterslow1_setpos] component PSD_afterslow1=PSD_monitor() SETTING [/usr/share/mcstas/3.0-dev/monitors/PSD_monitor.comp:64]");
  stracpy(_PSD_afterslow1->_name, "PSD_afterslow1", 16384);
  stracpy(_PSD_afterslow1->_type, "PSD_monitor", 16384);
  _PSD_afterslow1->_index=13;
  _PSD_afterslow1->_parameters.nx = 90;
  #define nx (_PSD_afterslow1->_parameters.nx)
  _PSD_afterslow1->_parameters.ny = 90;
  #define ny (_PSD_afterslow1->_parameters.ny)
  if("PSD_afterslow1.dat" && strlen("PSD_afterslow1.dat"))
    stracpy(_PSD_afterslow1->_parameters.filename, "PSD_afterslow1.dat" ? "PSD_afterslow1.dat" : "", 16384);
  else 
  _PSD_afterslow1->_parameters.filename[0]='\0';
  #define filename (_PSD_afterslow1->_parameters.filename)
  _PSD_afterslow1->_parameters.xmin = -0.05;
  #define xmin (_PSD_afterslow1->_parameters.xmin)
  _PSD_afterslow1->_parameters.xmax = 0.05;
  #define xmax (_PSD_afterslow1->_parameters.xmax)
  _PSD_afterslow1->_parameters.ymin = -0.05;
  #define ymin (_PSD_afterslow1->_parameters.ymin)
  _PSD_afterslow1->_parameters.ymax = 0.05;
  #define ymax (_PSD_afterslow1->_parameters.ymax)
  _PSD_afterslow1->_parameters.xwidth = 0.1;
  #define xwidth (_PSD_afterslow1->_parameters.xwidth)
  _PSD_afterslow1->_parameters.yheight = 0.25;
  #define yheight (_PSD_afterslow1->_parameters.yheight)
  _PSD_afterslow1->_parameters.restore_neutron = 0;
  #define restore_neutron (_PSD_afterslow1->_parameters.restore_neutron)

  #define PSD_N (_PSD_afterslow1->_parameters.PSD_N)
  #define PSD_p (_PSD_afterslow1->_parameters.PSD_p)
  #define PSD_p2 (_PSD_afterslow1->_parameters.PSD_p2)

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
  /* component PSD_afterslow1=PSD_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _FOchop1->_rotation_absolute, _PSD_afterslow1->_rotation_absolute);
    rot_transpose(_Lmon_afterslow1->_rotation_absolute, tr1);
    rot_mul(_PSD_afterslow1->_rotation_absolute, tr1, _PSD_afterslow1->_rotation_relative);
    _PSD_afterslow1->_rotation_is_identity =  rot_test_identity(_PSD_afterslow1->_rotation_relative);
    tc1 = coords_set(
      0, 0, 3e-6);
    rot_transpose(_FOchop1->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _PSD_afterslow1->_position_absolute = coords_add(_FOchop1->_position_absolute, tc2);
    tc1 = coords_sub(_Lmon_afterslow1->_position_absolute, _PSD_afterslow1->_position_absolute);
    _PSD_afterslow1->_position_relative = rot_apply(_PSD_afterslow1->_rotation_absolute, tc1);
  } /* PSD_afterslow1=PSD_monitor() AT ROTATED */
  DEBUG_COMPONENT("PSD_afterslow1", _PSD_afterslow1->_position_absolute, _PSD_afterslow1->_rotation_absolute);
  instrument->_position_absolute[13] = _PSD_afterslow1->_position_absolute;
  instrument->_position_relative[13] = _PSD_afterslow1->_position_relative;
  instrument->counter_N[13]  = instrument->counter_P[13] = instrument->counter_P2[13] = 0;
  instrument->counter_AbsorbProp[13]= 0;
  return(0);
} /* _PSD_afterslow1_setpos */

/* component Guidelong1=Guide() SETTING, POSITION/ROTATION */
int _Guidelong1_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Guidelong1_setpos] component Guidelong1=Guide() SETTING [/usr/share/mcstas/3.0-dev/optics/Guide.comp:73]");
  stracpy(_Guidelong1->_name, "Guidelong1", 16384);
  stracpy(_Guidelong1->_type, "Guide", 16384);
  _Guidelong1->_index=14;
  _Guidelong1->_parameters.reflect[0]='\0';
  #define reflect (_Guidelong1->_parameters.reflect)
  _Guidelong1->_parameters.w1 = instrument->_parameters._W1;
  #define w1 (_Guidelong1->_parameters.w1)
  _Guidelong1->_parameters.h1 = instrument->_parameters._H1;
  #define h1 (_Guidelong1->_parameters.h1)
  _Guidelong1->_parameters.w2 = instrument->_parameters._W2;
  #define w2 (_Guidelong1->_parameters.w2)
  _Guidelong1->_parameters.h2 = instrument->_parameters._H2;
  #define h2 (_Guidelong1->_parameters.h2)
  _Guidelong1->_parameters.l = instrument->_parameters._L_ballistic_begin + instrument->_parameters._GUI_start - instrument->_parameters._FO1_DIST - instrument->_parameters._GUI_GAP / 2;
  #define l (_Guidelong1->_parameters.l)
  _Guidelong1->_parameters.R0 = 1;
  #define R0 (_Guidelong1->_parameters.R0)
  _Guidelong1->_parameters.Qc = 0.0219;
  #define Qc (_Guidelong1->_parameters.Qc)
  _Guidelong1->_parameters.alpha = instrument->_parameters._ALPHA;
  #define alpha (_Guidelong1->_parameters.alpha)
  _Guidelong1->_parameters.m = instrument->_parameters._M;
  #define m (_Guidelong1->_parameters.m)
  _Guidelong1->_parameters.W = 0.003;
  #define W (_Guidelong1->_parameters.W)

  #define pTable (_Guidelong1->_parameters.pTable)

  #undef reflect
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
  #undef pTable
  /* component Guidelong1=Guide() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _FOchop1->_rotation_absolute, _Guidelong1->_rotation_absolute);
    rot_transpose(_PSD_afterslow1->_rotation_absolute, tr1);
    rot_mul(_Guidelong1->_rotation_absolute, tr1, _Guidelong1->_rotation_relative);
    _Guidelong1->_rotation_is_identity =  rot_test_identity(_Guidelong1->_rotation_relative);
    tc1 = coords_set(
      0, 0, instrument->_parameters._GUI_GAP / 2);
    rot_transpose(_FOchop1->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Guidelong1->_position_absolute = coords_add(_FOchop1->_position_absolute, tc2);
    tc1 = coords_sub(_PSD_afterslow1->_position_absolute, _Guidelong1->_position_absolute);
    _Guidelong1->_position_relative = rot_apply(_Guidelong1->_rotation_absolute, tc1);
  } /* Guidelong1=Guide() AT ROTATED */
  DEBUG_COMPONENT("Guidelong1", _Guidelong1->_position_absolute, _Guidelong1->_rotation_absolute);
  instrument->_position_absolute[14] = _Guidelong1->_position_absolute;
  instrument->_position_relative[14] = _Guidelong1->_position_relative;
  instrument->counter_N[14]  = instrument->counter_P[14] = instrument->counter_P2[14] = 0;
  instrument->counter_AbsorbProp[14]= 0;
  return(0);
} /* _Guidelong1_setpos */

/* component Guidelong1b=Guide() SETTING, POSITION/ROTATION */
int _Guidelong1b_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Guidelong1b_setpos] component Guidelong1b=Guide() SETTING [/usr/share/mcstas/3.0-dev/optics/Guide.comp:73]");
  stracpy(_Guidelong1b->_name, "Guidelong1b", 16384);
  stracpy(_Guidelong1b->_type, "Guide", 16384);
  _Guidelong1b->_index=15;
  _Guidelong1b->_parameters.reflect[0]='\0';
  #define reflect (_Guidelong1b->_parameters.reflect)
  _Guidelong1b->_parameters.w1 = instrument->_parameters._W2;
  #define w1 (_Guidelong1b->_parameters.w1)
  _Guidelong1b->_parameters.h1 = instrument->_parameters._H2;
  #define h1 (_Guidelong1b->_parameters.h1)
  _Guidelong1b->_parameters.w2 = instrument->_parameters._W3;
  #define w2 (_Guidelong1b->_parameters.w2)
  _Guidelong1b->_parameters.h2 = instrument->_parameters._H3;
  #define h2 (_Guidelong1b->_parameters.h2)
  _Guidelong1b->_parameters.l = instrument->_parameters._Length / 2 - instrument->_parameters._L_ballistic_begin - instrument->_parameters._GUI_start - instrument->_parameters._GUI_GAP / 2;
  #define l (_Guidelong1b->_parameters.l)
  _Guidelong1b->_parameters.R0 = 1;
  #define R0 (_Guidelong1b->_parameters.R0)
  _Guidelong1b->_parameters.Qc = 0.0219;
  #define Qc (_Guidelong1b->_parameters.Qc)
  _Guidelong1b->_parameters.alpha = instrument->_parameters._ALPHA;
  #define alpha (_Guidelong1b->_parameters.alpha)
  _Guidelong1b->_parameters.m = instrument->_parameters._M;
  #define m (_Guidelong1b->_parameters.m)
  _Guidelong1b->_parameters.W = 0.003;
  #define W (_Guidelong1b->_parameters.W)

  #define pTable (_Guidelong1b->_parameters.pTable)

  #undef reflect
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
  #undef pTable
  /* component Guidelong1b=Guide() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Origin->_rotation_absolute, _Guidelong1b->_rotation_absolute);
    rot_transpose(_Guidelong1->_rotation_absolute, tr1);
    rot_mul(_Guidelong1b->_rotation_absolute, tr1, _Guidelong1b->_rotation_relative);
    _Guidelong1b->_rotation_is_identity =  rot_test_identity(_Guidelong1b->_rotation_relative);
    tc1 = coords_set(
      0, 0, instrument->_parameters._L_ballistic_begin + instrument->_parameters._GUI_start);
    rot_transpose(_Origin->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Guidelong1b->_position_absolute = coords_add(_Origin->_position_absolute, tc2);
    tc1 = coords_sub(_Guidelong1->_position_absolute, _Guidelong1b->_position_absolute);
    _Guidelong1b->_position_relative = rot_apply(_Guidelong1b->_rotation_absolute, tc1);
  } /* Guidelong1b=Guide() AT ROTATED */
  DEBUG_COMPONENT("Guidelong1b", _Guidelong1b->_position_absolute, _Guidelong1b->_rotation_absolute);
  instrument->_position_absolute[15] = _Guidelong1b->_position_absolute;
  instrument->_position_relative[15] = _Guidelong1b->_position_relative;
  instrument->counter_N[15]  = instrument->counter_P[15] = instrument->counter_P2[15] = 0;
  instrument->counter_AbsorbProp[15]= 0;
  return(0);
} /* _Guidelong1b_setpos */

/* component Lmon_slow2=L_monitor() SETTING, POSITION/ROTATION */
int _Lmon_slow2_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Lmon_slow2_setpos] component Lmon_slow2=L_monitor() SETTING [/usr/share/mcstas/3.0-dev/monitors/L_monitor.comp:65]");
  stracpy(_Lmon_slow2->_name, "Lmon_slow2", 16384);
  stracpy(_Lmon_slow2->_type, "L_monitor", 16384);
  _Lmon_slow2->_index=16;
  _Lmon_slow2->_parameters.nL = 200;
  #define nL (_Lmon_slow2->_parameters.nL)
  if("Lmon_slow2.dat" && strlen("Lmon_slow2.dat"))
    stracpy(_Lmon_slow2->_parameters.filename, "Lmon_slow2.dat" ? "Lmon_slow2.dat" : "", 16384);
  else 
  _Lmon_slow2->_parameters.filename[0]='\0';
  #define filename (_Lmon_slow2->_parameters.filename)
  _Lmon_slow2->_parameters.xmin = -0.05;
  #define xmin (_Lmon_slow2->_parameters.xmin)
  _Lmon_slow2->_parameters.xmax = 0.05;
  #define xmax (_Lmon_slow2->_parameters.xmax)
  _Lmon_slow2->_parameters.ymin = -0.05;
  #define ymin (_Lmon_slow2->_parameters.ymin)
  _Lmon_slow2->_parameters.ymax = 0.05;
  #define ymax (_Lmon_slow2->_parameters.ymax)
  _Lmon_slow2->_parameters.xwidth = 0.06;
  #define xwidth (_Lmon_slow2->_parameters.xwidth)
  _Lmon_slow2->_parameters.yheight = 0.21;
  #define yheight (_Lmon_slow2->_parameters.yheight)
  _Lmon_slow2->_parameters.Lmin = 0;
  #define Lmin (_Lmon_slow2->_parameters.Lmin)
  _Lmon_slow2->_parameters.Lmax = instrument->_parameters._Lmax + 1;
  #define Lmax (_Lmon_slow2->_parameters.Lmax)
  _Lmon_slow2->_parameters.restore_neutron = 0;
  #define restore_neutron (_Lmon_slow2->_parameters.restore_neutron)

  #define L_N (_Lmon_slow2->_parameters.L_N)
  #define L_p (_Lmon_slow2->_parameters.L_p)
  #define L_p2 (_Lmon_slow2->_parameters.L_p2)

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
  /* component Lmon_slow2=L_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Origin->_rotation_absolute, _Lmon_slow2->_rotation_absolute);
    rot_transpose(_Guidelong1b->_rotation_absolute, tr1);
    rot_mul(_Lmon_slow2->_rotation_absolute, tr1, _Lmon_slow2->_rotation_relative);
    _Lmon_slow2->_rotation_is_identity =  rot_test_identity(_Lmon_slow2->_rotation_relative);
    tc1 = coords_set(
      0, 0, instrument->_parameters._Length / 2 -1e-6);
    rot_transpose(_Origin->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Lmon_slow2->_position_absolute = coords_add(_Origin->_position_absolute, tc2);
    tc1 = coords_sub(_Guidelong1b->_position_absolute, _Lmon_slow2->_position_absolute);
    _Lmon_slow2->_position_relative = rot_apply(_Lmon_slow2->_rotation_absolute, tc1);
  } /* Lmon_slow2=L_monitor() AT ROTATED */
  DEBUG_COMPONENT("Lmon_slow2", _Lmon_slow2->_position_absolute, _Lmon_slow2->_rotation_absolute);
  instrument->_position_absolute[16] = _Lmon_slow2->_position_absolute;
  instrument->_position_relative[16] = _Lmon_slow2->_position_relative;
  instrument->counter_N[16]  = instrument->counter_P[16] = instrument->counter_P2[16] = 0;
  instrument->counter_AbsorbProp[16]= 0;
  return(0);
} /* _Lmon_slow2_setpos */

/* component FOchop2=DiskChopper() SETTING, POSITION/ROTATION */
int _FOchop2_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_FOchop2_setpos] component FOchop2=DiskChopper() SETTING [/usr/share/mcstas/3.0-dev/optics/DiskChopper.comp:67]");
  stracpy(_FOchop2->_name, "FOchop2", 16384);
  stracpy(_FOchop2->_type, "DiskChopper", 16384);
  _FOchop2->_index=17;
  _FOchop2->_parameters.theta_0 = 155;
  #define theta_0 (_FOchop2->_parameters.theta_0)
  _FOchop2->_parameters.radius = 0.6;
  #define radius (_FOchop2->_parameters.radius)
  _FOchop2->_parameters.yheight = 0.6;
  #define yheight (_FOchop2->_parameters.yheight)
  _FOchop2->_parameters.nu = instrument->_parameters._F_slow2;
  #define nu (_FOchop2->_parameters.nu)
  _FOchop2->_parameters.nslit = 1;
  #define nslit (_FOchop2->_parameters.nslit)
  _FOchop2->_parameters.jitter = 0;
  #define jitter (_FOchop2->_parameters.jitter)
  _FOchop2->_parameters.delay = t_FO2;
  #define delay (_FOchop2->_parameters.delay)
  _FOchop2->_parameters.isfirst = 0;
  #define isfirst (_FOchop2->_parameters.isfirst)
  _FOchop2->_parameters.n_pulse = 1;
  #define n_pulse (_FOchop2->_parameters.n_pulse)
  _FOchop2->_parameters.abs_out = 1;
  #define abs_out (_FOchop2->_parameters.abs_out)
  _FOchop2->_parameters.phase = 0;
  #define phase (_FOchop2->_parameters.phase)
  _FOchop2->_parameters.xwidth = 0;
  #define xwidth (_FOchop2->_parameters.xwidth)
  _FOchop2->_parameters.verbose = 0;
  #define verbose (_FOchop2->_parameters.verbose)

  #define Tg (_FOchop2->_parameters.Tg)
  #define To (_FOchop2->_parameters.To)
  #define delta_y (_FOchop2->_parameters.delta_y)
  #define height (_FOchop2->_parameters.height)
  #define omega (_FOchop2->_parameters.omega)

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
  /* component FOchop2=DiskChopper() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Origin->_rotation_absolute, _FOchop2->_rotation_absolute);
    rot_transpose(_Lmon_slow2->_rotation_absolute, tr1);
    rot_mul(_FOchop2->_rotation_absolute, tr1, _FOchop2->_rotation_relative);
    _FOchop2->_rotation_is_identity =  rot_test_identity(_FOchop2->_rotation_relative);
    tc1 = coords_set(
      0, 0.1, instrument->_parameters._Length / 2);
    rot_transpose(_Origin->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _FOchop2->_position_absolute = coords_add(_Origin->_position_absolute, tc2);
    tc1 = coords_sub(_Lmon_slow2->_position_absolute, _FOchop2->_position_absolute);
    _FOchop2->_position_relative = rot_apply(_FOchop2->_rotation_absolute, tc1);
  } /* FOchop2=DiskChopper() AT ROTATED */
  DEBUG_COMPONENT("FOchop2", _FOchop2->_position_absolute, _FOchop2->_rotation_absolute);
  instrument->_position_absolute[17] = _FOchop2->_position_absolute;
  instrument->_position_relative[17] = _FOchop2->_position_relative;
  instrument->counter_N[17]  = instrument->counter_P[17] = instrument->counter_P2[17] = 0;
  instrument->counter_AbsorbProp[17]= 0;
  return(0);
} /* _FOchop2_setpos */

/* component Fastchop1=DiskChopper() SETTING, POSITION/ROTATION */
int _Fastchop1_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Fastchop1_setpos] component Fastchop1=DiskChopper() SETTING [/usr/share/mcstas/3.0-dev/optics/DiskChopper.comp:67]");
  stracpy(_Fastchop1->_name, "Fastchop1", 16384);
  stracpy(_Fastchop1->_type, "DiskChopper", 16384);
  _Fastchop1->_index=18;
  _Fastchop1->_parameters.theta_0 = instrument->_parameters._THETA_fast1;
  #define theta_0 (_Fastchop1->_parameters.theta_0)
  _Fastchop1->_parameters.radius = 0.5;
  #define radius (_Fastchop1->_parameters.radius)
  _Fastchop1->_parameters.yheight = 0.5;
  #define yheight (_Fastchop1->_parameters.yheight)
  _Fastchop1->_parameters.nu = instrument->_parameters._F_fast1;
  #define nu (_Fastchop1->_parameters.nu)
  _Fastchop1->_parameters.nslit = instrument->_parameters._N_fast;
  #define nslit (_Fastchop1->_parameters.nslit)
  _Fastchop1->_parameters.jitter = 0;
  #define jitter (_Fastchop1->_parameters.jitter)
  _Fastchop1->_parameters.delay = t_fast1;
  #define delay (_Fastchop1->_parameters.delay)
  _Fastchop1->_parameters.isfirst = 0;
  #define isfirst (_Fastchop1->_parameters.isfirst)
  _Fastchop1->_parameters.n_pulse = 1;
  #define n_pulse (_Fastchop1->_parameters.n_pulse)
  _Fastchop1->_parameters.abs_out = 1;
  #define abs_out (_Fastchop1->_parameters.abs_out)
  _Fastchop1->_parameters.phase = 0;
  #define phase (_Fastchop1->_parameters.phase)
  _Fastchop1->_parameters.xwidth = 0;
  #define xwidth (_Fastchop1->_parameters.xwidth)
  _Fastchop1->_parameters.verbose = 0;
  #define verbose (_Fastchop1->_parameters.verbose)

  #define Tg (_Fastchop1->_parameters.Tg)
  #define To (_Fastchop1->_parameters.To)
  #define delta_y (_Fastchop1->_parameters.delta_y)
  #define height (_Fastchop1->_parameters.height)
  #define omega (_Fastchop1->_parameters.omega)

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
  /* component Fastchop1=DiskChopper() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Origin->_rotation_absolute, _Fastchop1->_rotation_absolute);
    rot_transpose(_FOchop2->_rotation_absolute, tr1);
    rot_mul(_Fastchop1->_rotation_absolute, tr1, _Fastchop1->_rotation_relative);
    _Fastchop1->_rotation_is_identity =  rot_test_identity(_Fastchop1->_rotation_relative);
    tc1 = coords_set(
      0, 0.1, instrument->_parameters._Length / 2 + instrument->_parameters._GUI_GAP);
    rot_transpose(_Origin->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Fastchop1->_position_absolute = coords_add(_Origin->_position_absolute, tc2);
    tc1 = coords_sub(_FOchop2->_position_absolute, _Fastchop1->_position_absolute);
    _Fastchop1->_position_relative = rot_apply(_Fastchop1->_rotation_absolute, tc1);
  } /* Fastchop1=DiskChopper() AT ROTATED */
  DEBUG_COMPONENT("Fastchop1", _Fastchop1->_position_absolute, _Fastchop1->_rotation_absolute);
  instrument->_position_absolute[18] = _Fastchop1->_position_absolute;
  instrument->_position_relative[18] = _Fastchop1->_position_relative;
  instrument->counter_N[18]  = instrument->counter_P[18] = instrument->counter_P2[18] = 0;
  instrument->counter_AbsorbProp[18]= 0;
  return(0);
} /* _Fastchop1_setpos */

/* component PSD_afterslow2=PSD_monitor() SETTING, POSITION/ROTATION */
int _PSD_afterslow2_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_PSD_afterslow2_setpos] component PSD_afterslow2=PSD_monitor() SETTING [/usr/share/mcstas/3.0-dev/monitors/PSD_monitor.comp:64]");
  stracpy(_PSD_afterslow2->_name, "PSD_afterslow2", 16384);
  stracpy(_PSD_afterslow2->_type, "PSD_monitor", 16384);
  _PSD_afterslow2->_index=19;
  _PSD_afterslow2->_parameters.nx = 90;
  #define nx (_PSD_afterslow2->_parameters.nx)
  _PSD_afterslow2->_parameters.ny = 90;
  #define ny (_PSD_afterslow2->_parameters.ny)
  if("PSD_afterslow2.dat" && strlen("PSD_afterslow2.dat"))
    stracpy(_PSD_afterslow2->_parameters.filename, "PSD_afterslow2.dat" ? "PSD_afterslow2.dat" : "", 16384);
  else 
  _PSD_afterslow2->_parameters.filename[0]='\0';
  #define filename (_PSD_afterslow2->_parameters.filename)
  _PSD_afterslow2->_parameters.xmin = -0.05;
  #define xmin (_PSD_afterslow2->_parameters.xmin)
  _PSD_afterslow2->_parameters.xmax = 0.05;
  #define xmax (_PSD_afterslow2->_parameters.xmax)
  _PSD_afterslow2->_parameters.ymin = -0.05;
  #define ymin (_PSD_afterslow2->_parameters.ymin)
  _PSD_afterslow2->_parameters.ymax = 0.05;
  #define ymax (_PSD_afterslow2->_parameters.ymax)
  _PSD_afterslow2->_parameters.xwidth = 0.1;
  #define xwidth (_PSD_afterslow2->_parameters.xwidth)
  _PSD_afterslow2->_parameters.yheight = 0.25;
  #define yheight (_PSD_afterslow2->_parameters.yheight)
  _PSD_afterslow2->_parameters.restore_neutron = 0;
  #define restore_neutron (_PSD_afterslow2->_parameters.restore_neutron)

  #define PSD_N (_PSD_afterslow2->_parameters.PSD_N)
  #define PSD_p (_PSD_afterslow2->_parameters.PSD_p)
  #define PSD_p2 (_PSD_afterslow2->_parameters.PSD_p2)

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
  /* component PSD_afterslow2=PSD_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Fastchop1->_rotation_absolute, _PSD_afterslow2->_rotation_absolute);
    rot_transpose(_Fastchop1->_rotation_absolute, tr1);
    rot_mul(_PSD_afterslow2->_rotation_absolute, tr1, _PSD_afterslow2->_rotation_relative);
    _PSD_afterslow2->_rotation_is_identity =  rot_test_identity(_PSD_afterslow2->_rotation_relative);
    tc1 = coords_set(
      0, -0.1, 1e-6);
    rot_transpose(_Fastchop1->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _PSD_afterslow2->_position_absolute = coords_add(_Fastchop1->_position_absolute, tc2);
    tc1 = coords_sub(_Fastchop1->_position_absolute, _PSD_afterslow2->_position_absolute);
    _PSD_afterslow2->_position_relative = rot_apply(_PSD_afterslow2->_rotation_absolute, tc1);
  } /* PSD_afterslow2=PSD_monitor() AT ROTATED */
  DEBUG_COMPONENT("PSD_afterslow2", _PSD_afterslow2->_position_absolute, _PSD_afterslow2->_rotation_absolute);
  instrument->_position_absolute[19] = _PSD_afterslow2->_position_absolute;
  instrument->_position_relative[19] = _PSD_afterslow2->_position_relative;
  instrument->counter_N[19]  = instrument->counter_P[19] = instrument->counter_P2[19] = 0;
  instrument->counter_AbsorbProp[19]= 0;
  return(0);
} /* _PSD_afterslow2_setpos */

/* component Lmon_afterslow2=L_monitor() SETTING, POSITION/ROTATION */
int _Lmon_afterslow2_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Lmon_afterslow2_setpos] component Lmon_afterslow2=L_monitor() SETTING [/usr/share/mcstas/3.0-dev/monitors/L_monitor.comp:65]");
  stracpy(_Lmon_afterslow2->_name, "Lmon_afterslow2", 16384);
  stracpy(_Lmon_afterslow2->_type, "L_monitor", 16384);
  _Lmon_afterslow2->_index=20;
  _Lmon_afterslow2->_parameters.nL = 200;
  #define nL (_Lmon_afterslow2->_parameters.nL)
  if("Lmon_afterslow2.dat" && strlen("Lmon_afterslow2.dat"))
    stracpy(_Lmon_afterslow2->_parameters.filename, "Lmon_afterslow2.dat" ? "Lmon_afterslow2.dat" : "", 16384);
  else 
  _Lmon_afterslow2->_parameters.filename[0]='\0';
  #define filename (_Lmon_afterslow2->_parameters.filename)
  _Lmon_afterslow2->_parameters.xmin = -0.05;
  #define xmin (_Lmon_afterslow2->_parameters.xmin)
  _Lmon_afterslow2->_parameters.xmax = 0.05;
  #define xmax (_Lmon_afterslow2->_parameters.xmax)
  _Lmon_afterslow2->_parameters.ymin = -0.05;
  #define ymin (_Lmon_afterslow2->_parameters.ymin)
  _Lmon_afterslow2->_parameters.ymax = 0.05;
  #define ymax (_Lmon_afterslow2->_parameters.ymax)
  _Lmon_afterslow2->_parameters.xwidth = 0.06;
  #define xwidth (_Lmon_afterslow2->_parameters.xwidth)
  _Lmon_afterslow2->_parameters.yheight = 0.21;
  #define yheight (_Lmon_afterslow2->_parameters.yheight)
  _Lmon_afterslow2->_parameters.Lmin = 0;
  #define Lmin (_Lmon_afterslow2->_parameters.Lmin)
  _Lmon_afterslow2->_parameters.Lmax = instrument->_parameters._Lmax + 1;
  #define Lmax (_Lmon_afterslow2->_parameters.Lmax)
  _Lmon_afterslow2->_parameters.restore_neutron = 0;
  #define restore_neutron (_Lmon_afterslow2->_parameters.restore_neutron)

  #define L_N (_Lmon_afterslow2->_parameters.L_N)
  #define L_p (_Lmon_afterslow2->_parameters.L_p)
  #define L_p2 (_Lmon_afterslow2->_parameters.L_p2)

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
  /* component Lmon_afterslow2=L_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _PSD_afterslow2->_rotation_absolute, _Lmon_afterslow2->_rotation_absolute);
    rot_transpose(_PSD_afterslow2->_rotation_absolute, tr1);
    rot_mul(_Lmon_afterslow2->_rotation_absolute, tr1, _Lmon_afterslow2->_rotation_relative);
    _Lmon_afterslow2->_rotation_is_identity =  rot_test_identity(_Lmon_afterslow2->_rotation_relative);
    tc1 = coords_set(
      0, 0, 1e-6);
    rot_transpose(_PSD_afterslow2->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Lmon_afterslow2->_position_absolute = coords_add(_PSD_afterslow2->_position_absolute, tc2);
    tc1 = coords_sub(_PSD_afterslow2->_position_absolute, _Lmon_afterslow2->_position_absolute);
    _Lmon_afterslow2->_position_relative = rot_apply(_Lmon_afterslow2->_rotation_absolute, tc1);
  } /* Lmon_afterslow2=L_monitor() AT ROTATED */
  DEBUG_COMPONENT("Lmon_afterslow2", _Lmon_afterslow2->_position_absolute, _Lmon_afterslow2->_rotation_absolute);
  instrument->_position_absolute[20] = _Lmon_afterslow2->_position_absolute;
  instrument->_position_relative[20] = _Lmon_afterslow2->_position_relative;
  instrument->counter_N[20]  = instrument->counter_P[20] = instrument->counter_P2[20] = 0;
  instrument->counter_AbsorbProp[20]= 0;
  return(0);
} /* _Lmon_afterslow2_setpos */

/* component TOFL_afterslow2=TOFLambda_monitor() SETTING, POSITION/ROTATION */
int _TOFL_afterslow2_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_TOFL_afterslow2_setpos] component TOFL_afterslow2=TOFLambda_monitor() SETTING [/usr/share/mcstas/3.0-dev/monitors/TOFLambda_monitor.comp:68]");
  stracpy(_TOFL_afterslow2->_name, "TOFL_afterslow2", 16384);
  stracpy(_TOFL_afterslow2->_type, "TOFLambda_monitor", 16384);
  _TOFL_afterslow2->_index=21;
  _TOFL_afterslow2->_parameters.nL = 200;
  #define nL (_TOFL_afterslow2->_parameters.nL)
  _TOFL_afterslow2->_parameters.nt = 200;
  #define nt (_TOFL_afterslow2->_parameters.nt)
  _TOFL_afterslow2->_parameters.tmin = 0;
  #define tmin (_TOFL_afterslow2->_parameters.tmin)
  _TOFL_afterslow2->_parameters.tmax = 3.0e5;
  #define tmax (_TOFL_afterslow2->_parameters.tmax)
  if("TOFL_afterslow2.dat" && strlen("TOFL_afterslow2.dat"))
    stracpy(_TOFL_afterslow2->_parameters.filename, "TOFL_afterslow2.dat" ? "TOFL_afterslow2.dat" : "", 16384);
  else 
  _TOFL_afterslow2->_parameters.filename[0]='\0';
  #define filename (_TOFL_afterslow2->_parameters.filename)
  _TOFL_afterslow2->_parameters.xmin = -0.05;
  #define xmin (_TOFL_afterslow2->_parameters.xmin)
  _TOFL_afterslow2->_parameters.xmax = 0.05;
  #define xmax (_TOFL_afterslow2->_parameters.xmax)
  _TOFL_afterslow2->_parameters.ymin = -0.05;
  #define ymin (_TOFL_afterslow2->_parameters.ymin)
  _TOFL_afterslow2->_parameters.ymax = 0.05;
  #define ymax (_TOFL_afterslow2->_parameters.ymax)
  _TOFL_afterslow2->_parameters.xwidth = 0.05;
  #define xwidth (_TOFL_afterslow2->_parameters.xwidth)
  _TOFL_afterslow2->_parameters.yheight = 0.21;
  #define yheight (_TOFL_afterslow2->_parameters.yheight)
  _TOFL_afterslow2->_parameters.Lmin = 0;
  #define Lmin (_TOFL_afterslow2->_parameters.Lmin)
  _TOFL_afterslow2->_parameters.Lmax = instrument->_parameters._Lmax + 1;
  #define Lmax (_TOFL_afterslow2->_parameters.Lmax)
  _TOFL_afterslow2->_parameters.restore_neutron = 0;
  #define restore_neutron (_TOFL_afterslow2->_parameters.restore_neutron)

  #define tt_0 (_TOFL_afterslow2->_parameters.tt_0)
  #define tt_1 (_TOFL_afterslow2->_parameters.tt_1)
  #define TOFL_N (_TOFL_afterslow2->_parameters.TOFL_N)
  #define TOFL_p (_TOFL_afterslow2->_parameters.TOFL_p)
  #define TOFL_p2 (_TOFL_afterslow2->_parameters.TOFL_p2)

  #undef nL
  #undef nt
  #undef tmin
  #undef tmax
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
  #undef tt_0
  #undef tt_1
  #undef TOFL_N
  #undef TOFL_p
  #undef TOFL_p2
  /* component TOFL_afterslow2=TOFLambda_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Lmon_afterslow2->_rotation_absolute, _TOFL_afterslow2->_rotation_absolute);
    rot_transpose(_Lmon_afterslow2->_rotation_absolute, tr1);
    rot_mul(_TOFL_afterslow2->_rotation_absolute, tr1, _TOFL_afterslow2->_rotation_relative);
    _TOFL_afterslow2->_rotation_is_identity =  rot_test_identity(_TOFL_afterslow2->_rotation_relative);
    tc1 = coords_set(
      0, 0, 1e-6);
    rot_transpose(_Lmon_afterslow2->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _TOFL_afterslow2->_position_absolute = coords_add(_Lmon_afterslow2->_position_absolute, tc2);
    tc1 = coords_sub(_Lmon_afterslow2->_position_absolute, _TOFL_afterslow2->_position_absolute);
    _TOFL_afterslow2->_position_relative = rot_apply(_TOFL_afterslow2->_rotation_absolute, tc1);
  } /* TOFL_afterslow2=TOFLambda_monitor() AT ROTATED */
  DEBUG_COMPONENT("TOFL_afterslow2", _TOFL_afterslow2->_position_absolute, _TOFL_afterslow2->_rotation_absolute);
  instrument->_position_absolute[21] = _TOFL_afterslow2->_position_absolute;
  instrument->_position_relative[21] = _TOFL_afterslow2->_position_relative;
  instrument->counter_N[21]  = instrument->counter_P[21] = instrument->counter_P2[21] = 0;
  instrument->counter_AbsorbProp[21]= 0;
  return(0);
} /* _TOFL_afterslow2_setpos */

/* component Guidelong2=Guide() SETTING, POSITION/ROTATION */
int _Guidelong2_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Guidelong2_setpos] component Guidelong2=Guide() SETTING [/usr/share/mcstas/3.0-dev/optics/Guide.comp:73]");
  stracpy(_Guidelong2->_name, "Guidelong2", 16384);
  stracpy(_Guidelong2->_type, "Guide", 16384);
  _Guidelong2->_index=22;
  _Guidelong2->_parameters.reflect[0]='\0';
  #define reflect (_Guidelong2->_parameters.reflect)
  _Guidelong2->_parameters.w1 = instrument->_parameters._W3;
  #define w1 (_Guidelong2->_parameters.w1)
  _Guidelong2->_parameters.h1 = instrument->_parameters._H3;
  #define h1 (_Guidelong2->_parameters.h1)
  _Guidelong2->_parameters.w2 = instrument->_parameters._W4;
  #define w2 (_Guidelong2->_parameters.w2)
  _Guidelong2->_parameters.h2 = instrument->_parameters._H4;
  #define h2 (_Guidelong2->_parameters.h2)
  _Guidelong2->_parameters.l = instrument->_parameters._Length / 2 -2 * instrument->_parameters._GUI_GAP - instrument->_parameters._L_ballistic_end;
  #define l (_Guidelong2->_parameters.l)
  _Guidelong2->_parameters.R0 = 1;
  #define R0 (_Guidelong2->_parameters.R0)
  _Guidelong2->_parameters.Qc = 0.0219;
  #define Qc (_Guidelong2->_parameters.Qc)
  _Guidelong2->_parameters.alpha = instrument->_parameters._ALPHA;
  #define alpha (_Guidelong2->_parameters.alpha)
  _Guidelong2->_parameters.m = instrument->_parameters._M;
  #define m (_Guidelong2->_parameters.m)
  _Guidelong2->_parameters.W = 0.003;
  #define W (_Guidelong2->_parameters.W)

  #define pTable (_Guidelong2->_parameters.pTable)

  #undef reflect
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
  #undef pTable
  /* component Guidelong2=Guide() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Fastchop1->_rotation_absolute, _Guidelong2->_rotation_absolute);
    rot_transpose(_TOFL_afterslow2->_rotation_absolute, tr1);
    rot_mul(_Guidelong2->_rotation_absolute, tr1, _Guidelong2->_rotation_relative);
    _Guidelong2->_rotation_is_identity =  rot_test_identity(_Guidelong2->_rotation_relative);
    tc1 = coords_set(
      0, -0.1, instrument->_parameters._GUI_GAP / 2);
    rot_transpose(_Fastchop1->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Guidelong2->_position_absolute = coords_add(_Fastchop1->_position_absolute, tc2);
    tc1 = coords_sub(_TOFL_afterslow2->_position_absolute, _Guidelong2->_position_absolute);
    _Guidelong2->_position_relative = rot_apply(_Guidelong2->_rotation_absolute, tc1);
  } /* Guidelong2=Guide() AT ROTATED */
  DEBUG_COMPONENT("Guidelong2", _Guidelong2->_position_absolute, _Guidelong2->_rotation_absolute);
  instrument->_position_absolute[22] = _Guidelong2->_position_absolute;
  instrument->_position_relative[22] = _Guidelong2->_position_relative;
  instrument->counter_N[22]  = instrument->counter_P[22] = instrument->counter_P2[22] = 0;
  instrument->counter_AbsorbProp[22]= 0;
  return(0);
} /* _Guidelong2_setpos */

/* component Lmon_beforeballistic=L_monitor() SETTING, POSITION/ROTATION */
int _Lmon_beforeballistic_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Lmon_beforeballistic_setpos] component Lmon_beforeballistic=L_monitor() SETTING [/usr/share/mcstas/3.0-dev/monitors/L_monitor.comp:65]");
  stracpy(_Lmon_beforeballistic->_name, "Lmon_beforeballistic", 16384);
  stracpy(_Lmon_beforeballistic->_type, "L_monitor", 16384);
  _Lmon_beforeballistic->_index=23;
  _Lmon_beforeballistic->_parameters.nL = 200;
  #define nL (_Lmon_beforeballistic->_parameters.nL)
  if("Lmon_before_ballistic.dat" && strlen("Lmon_before_ballistic.dat"))
    stracpy(_Lmon_beforeballistic->_parameters.filename, "Lmon_before_ballistic.dat" ? "Lmon_before_ballistic.dat" : "", 16384);
  else 
  _Lmon_beforeballistic->_parameters.filename[0]='\0';
  #define filename (_Lmon_beforeballistic->_parameters.filename)
  _Lmon_beforeballistic->_parameters.xmin = -0.05;
  #define xmin (_Lmon_beforeballistic->_parameters.xmin)
  _Lmon_beforeballistic->_parameters.xmax = 0.05;
  #define xmax (_Lmon_beforeballistic->_parameters.xmax)
  _Lmon_beforeballistic->_parameters.ymin = -0.05;
  #define ymin (_Lmon_beforeballistic->_parameters.ymin)
  _Lmon_beforeballistic->_parameters.ymax = 0.05;
  #define ymax (_Lmon_beforeballistic->_parameters.ymax)
  _Lmon_beforeballistic->_parameters.xwidth = 0.06;
  #define xwidth (_Lmon_beforeballistic->_parameters.xwidth)
  _Lmon_beforeballistic->_parameters.yheight = 0.18;
  #define yheight (_Lmon_beforeballistic->_parameters.yheight)
  _Lmon_beforeballistic->_parameters.Lmin = 0;
  #define Lmin (_Lmon_beforeballistic->_parameters.Lmin)
  _Lmon_beforeballistic->_parameters.Lmax = instrument->_parameters._Lmax + 1;
  #define Lmax (_Lmon_beforeballistic->_parameters.Lmax)
  _Lmon_beforeballistic->_parameters.restore_neutron = 0;
  #define restore_neutron (_Lmon_beforeballistic->_parameters.restore_neutron)

  #define L_N (_Lmon_beforeballistic->_parameters.L_N)
  #define L_p (_Lmon_beforeballistic->_parameters.L_p)
  #define L_p2 (_Lmon_beforeballistic->_parameters.L_p2)

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
  /* component Lmon_beforeballistic=L_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Guidelong2->_rotation_absolute, _Lmon_beforeballistic->_rotation_absolute);
    rot_transpose(_Guidelong2->_rotation_absolute, tr1);
    rot_mul(_Lmon_beforeballistic->_rotation_absolute, tr1, _Lmon_beforeballistic->_rotation_relative);
    _Lmon_beforeballistic->_rotation_is_identity =  rot_test_identity(_Lmon_beforeballistic->_rotation_relative);
    tc1 = coords_set(
      0, 0, instrument->_parameters._Length / 2 -2 * instrument->_parameters._GUI_GAP - instrument->_parameters._L_ballistic_end + 1e-6);
    rot_transpose(_Guidelong2->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Lmon_beforeballistic->_position_absolute = coords_add(_Guidelong2->_position_absolute, tc2);
    tc1 = coords_sub(_Guidelong2->_position_absolute, _Lmon_beforeballistic->_position_absolute);
    _Lmon_beforeballistic->_position_relative = rot_apply(_Lmon_beforeballistic->_rotation_absolute, tc1);
  } /* Lmon_beforeballistic=L_monitor() AT ROTATED */
  DEBUG_COMPONENT("Lmon_beforeballistic", _Lmon_beforeballistic->_position_absolute, _Lmon_beforeballistic->_rotation_absolute);
  instrument->_position_absolute[23] = _Lmon_beforeballistic->_position_absolute;
  instrument->_position_relative[23] = _Lmon_beforeballistic->_position_relative;
  instrument->counter_N[23]  = instrument->counter_P[23] = instrument->counter_P2[23] = 0;
  instrument->counter_AbsorbProp[23]= 0;
  return(0);
} /* _Lmon_beforeballistic_setpos */

/* component PSD_beforeballistic=PSD_monitor() SETTING, POSITION/ROTATION */
int _PSD_beforeballistic_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_PSD_beforeballistic_setpos] component PSD_beforeballistic=PSD_monitor() SETTING [/usr/share/mcstas/3.0-dev/monitors/PSD_monitor.comp:64]");
  stracpy(_PSD_beforeballistic->_name, "PSD_beforeballistic", 16384);
  stracpy(_PSD_beforeballistic->_type, "PSD_monitor", 16384);
  _PSD_beforeballistic->_index=24;
  _PSD_beforeballistic->_parameters.nx = 90;
  #define nx (_PSD_beforeballistic->_parameters.nx)
  _PSD_beforeballistic->_parameters.ny = 90;
  #define ny (_PSD_beforeballistic->_parameters.ny)
  if("PSD_beforeballistic.dat" && strlen("PSD_beforeballistic.dat"))
    stracpy(_PSD_beforeballistic->_parameters.filename, "PSD_beforeballistic.dat" ? "PSD_beforeballistic.dat" : "", 16384);
  else 
  _PSD_beforeballistic->_parameters.filename[0]='\0';
  #define filename (_PSD_beforeballistic->_parameters.filename)
  _PSD_beforeballistic->_parameters.xmin = -0.05;
  #define xmin (_PSD_beforeballistic->_parameters.xmin)
  _PSD_beforeballistic->_parameters.xmax = 0.05;
  #define xmax (_PSD_beforeballistic->_parameters.xmax)
  _PSD_beforeballistic->_parameters.ymin = -0.05;
  #define ymin (_PSD_beforeballistic->_parameters.ymin)
  _PSD_beforeballistic->_parameters.ymax = 0.05;
  #define ymax (_PSD_beforeballistic->_parameters.ymax)
  _PSD_beforeballistic->_parameters.xwidth = 0.1;
  #define xwidth (_PSD_beforeballistic->_parameters.xwidth)
  _PSD_beforeballistic->_parameters.yheight = 0.25;
  #define yheight (_PSD_beforeballistic->_parameters.yheight)
  _PSD_beforeballistic->_parameters.restore_neutron = 0;
  #define restore_neutron (_PSD_beforeballistic->_parameters.restore_neutron)

  #define PSD_N (_PSD_beforeballistic->_parameters.PSD_N)
  #define PSD_p (_PSD_beforeballistic->_parameters.PSD_p)
  #define PSD_p2 (_PSD_beforeballistic->_parameters.PSD_p2)

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
  /* component PSD_beforeballistic=PSD_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Lmon_beforeballistic->_rotation_absolute, _PSD_beforeballistic->_rotation_absolute);
    rot_transpose(_Lmon_beforeballistic->_rotation_absolute, tr1);
    rot_mul(_PSD_beforeballistic->_rotation_absolute, tr1, _PSD_beforeballistic->_rotation_relative);
    _PSD_beforeballistic->_rotation_is_identity =  rot_test_identity(_PSD_beforeballistic->_rotation_relative);
    tc1 = coords_set(
      0, 0, 1e-6);
    rot_transpose(_Lmon_beforeballistic->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _PSD_beforeballistic->_position_absolute = coords_add(_Lmon_beforeballistic->_position_absolute, tc2);
    tc1 = coords_sub(_Lmon_beforeballistic->_position_absolute, _PSD_beforeballistic->_position_absolute);
    _PSD_beforeballistic->_position_relative = rot_apply(_PSD_beforeballistic->_rotation_absolute, tc1);
  } /* PSD_beforeballistic=PSD_monitor() AT ROTATED */
  DEBUG_COMPONENT("PSD_beforeballistic", _PSD_beforeballistic->_position_absolute, _PSD_beforeballistic->_rotation_absolute);
  instrument->_position_absolute[24] = _PSD_beforeballistic->_position_absolute;
  instrument->_position_relative[24] = _PSD_beforeballistic->_position_relative;
  instrument->counter_N[24]  = instrument->counter_P[24] = instrument->counter_P2[24] = 0;
  instrument->counter_AbsorbProp[24]= 0;
  return(0);
} /* _PSD_beforeballistic_setpos */

/* component Guidelong2a=Guide() SETTING, POSITION/ROTATION */
int _Guidelong2a_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Guidelong2a_setpos] component Guidelong2a=Guide() SETTING [/usr/share/mcstas/3.0-dev/optics/Guide.comp:73]");
  stracpy(_Guidelong2a->_name, "Guidelong2a", 16384);
  stracpy(_Guidelong2a->_type, "Guide", 16384);
  _Guidelong2a->_index=25;
  _Guidelong2a->_parameters.reflect[0]='\0';
  #define reflect (_Guidelong2a->_parameters.reflect)
  _Guidelong2a->_parameters.w1 = instrument->_parameters._W4;
  #define w1 (_Guidelong2a->_parameters.w1)
  _Guidelong2a->_parameters.h1 = instrument->_parameters._H4;
  #define h1 (_Guidelong2a->_parameters.h1)
  _Guidelong2a->_parameters.w2 = instrument->_parameters._W_chop;
  #define w2 (_Guidelong2a->_parameters.w2)
  _Guidelong2a->_parameters.h2 = instrument->_parameters._H_chop;
  #define h2 (_Guidelong2a->_parameters.h2)
  _Guidelong2a->_parameters.l = instrument->_parameters._L_ballistic_end;
  #define l (_Guidelong2a->_parameters.l)
  _Guidelong2a->_parameters.R0 = 1;
  #define R0 (_Guidelong2a->_parameters.R0)
  _Guidelong2a->_parameters.Qc = 0.0219;
  #define Qc (_Guidelong2a->_parameters.Qc)
  _Guidelong2a->_parameters.alpha = instrument->_parameters._ALPHA;
  #define alpha (_Guidelong2a->_parameters.alpha)
  _Guidelong2a->_parameters.m = instrument->_parameters._M;
  #define m (_Guidelong2a->_parameters.m)
  _Guidelong2a->_parameters.W = 0.003;
  #define W (_Guidelong2a->_parameters.W)

  #define pTable (_Guidelong2a->_parameters.pTable)

  #undef reflect
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
  #undef pTable
  /* component Guidelong2a=Guide() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _PSD_beforeballistic->_rotation_absolute, _Guidelong2a->_rotation_absolute);
    rot_transpose(_PSD_beforeballistic->_rotation_absolute, tr1);
    rot_mul(_Guidelong2a->_rotation_absolute, tr1, _Guidelong2a->_rotation_relative);
    _Guidelong2a->_rotation_is_identity =  rot_test_identity(_Guidelong2a->_rotation_relative);
    tc1 = coords_set(
      0, 0, 1e-6);
    rot_transpose(_PSD_beforeballistic->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Guidelong2a->_position_absolute = coords_add(_PSD_beforeballistic->_position_absolute, tc2);
    tc1 = coords_sub(_PSD_beforeballistic->_position_absolute, _Guidelong2a->_position_absolute);
    _Guidelong2a->_position_relative = rot_apply(_Guidelong2a->_rotation_absolute, tc1);
  } /* Guidelong2a=Guide() AT ROTATED */
  DEBUG_COMPONENT("Guidelong2a", _Guidelong2a->_position_absolute, _Guidelong2a->_rotation_absolute);
  instrument->_position_absolute[25] = _Guidelong2a->_position_absolute;
  instrument->_position_relative[25] = _Guidelong2a->_position_relative;
  instrument->counter_N[25]  = instrument->counter_P[25] = instrument->counter_P2[25] = 0;
  instrument->counter_AbsorbProp[25]= 0;
  return(0);
} /* _Guidelong2a_setpos */

/* component Lmonfast2=L_monitor() SETTING, POSITION/ROTATION */
int _Lmonfast2_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Lmonfast2_setpos] component Lmonfast2=L_monitor() SETTING [/usr/share/mcstas/3.0-dev/monitors/L_monitor.comp:65]");
  stracpy(_Lmonfast2->_name, "Lmonfast2", 16384);
  stracpy(_Lmonfast2->_type, "L_monitor", 16384);
  _Lmonfast2->_index=26;
  _Lmonfast2->_parameters.nL = 200;
  #define nL (_Lmonfast2->_parameters.nL)
  if("Lmonfast2.dat" && strlen("Lmonfast2.dat"))
    stracpy(_Lmonfast2->_parameters.filename, "Lmonfast2.dat" ? "Lmonfast2.dat" : "", 16384);
  else 
  _Lmonfast2->_parameters.filename[0]='\0';
  #define filename (_Lmonfast2->_parameters.filename)
  _Lmonfast2->_parameters.xmin = -0.05;
  #define xmin (_Lmonfast2->_parameters.xmin)
  _Lmonfast2->_parameters.xmax = 0.05;
  #define xmax (_Lmonfast2->_parameters.xmax)
  _Lmonfast2->_parameters.ymin = -0.05;
  #define ymin (_Lmonfast2->_parameters.ymin)
  _Lmonfast2->_parameters.ymax = 0.05;
  #define ymax (_Lmonfast2->_parameters.ymax)
  _Lmonfast2->_parameters.xwidth = 0.06;
  #define xwidth (_Lmonfast2->_parameters.xwidth)
  _Lmonfast2->_parameters.yheight = 0.18;
  #define yheight (_Lmonfast2->_parameters.yheight)
  _Lmonfast2->_parameters.Lmin = 0;
  #define Lmin (_Lmonfast2->_parameters.Lmin)
  _Lmonfast2->_parameters.Lmax = instrument->_parameters._Lmax + 1;
  #define Lmax (_Lmonfast2->_parameters.Lmax)
  _Lmonfast2->_parameters.restore_neutron = 0;
  #define restore_neutron (_Lmonfast2->_parameters.restore_neutron)

  #define L_N (_Lmonfast2->_parameters.L_N)
  #define L_p (_Lmonfast2->_parameters.L_p)
  #define L_p2 (_Lmonfast2->_parameters.L_p2)

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
  /* component Lmonfast2=L_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Guidelong2a->_rotation_absolute, _Lmonfast2->_rotation_absolute);
    rot_transpose(_Guidelong2a->_rotation_absolute, tr1);
    rot_mul(_Lmonfast2->_rotation_absolute, tr1, _Lmonfast2->_rotation_relative);
    _Lmonfast2->_rotation_is_identity =  rot_test_identity(_Lmonfast2->_rotation_relative);
    tc1 = coords_set(
      0, 0, instrument->_parameters._L_ballistic_end + 1e-6);
    rot_transpose(_Guidelong2a->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Lmonfast2->_position_absolute = coords_add(_Guidelong2a->_position_absolute, tc2);
    tc1 = coords_sub(_Guidelong2a->_position_absolute, _Lmonfast2->_position_absolute);
    _Lmonfast2->_position_relative = rot_apply(_Lmonfast2->_rotation_absolute, tc1);
  } /* Lmonfast2=L_monitor() AT ROTATED */
  DEBUG_COMPONENT("Lmonfast2", _Lmonfast2->_position_absolute, _Lmonfast2->_rotation_absolute);
  instrument->_position_absolute[26] = _Lmonfast2->_position_absolute;
  instrument->_position_relative[26] = _Lmonfast2->_position_relative;
  instrument->counter_N[26]  = instrument->counter_P[26] = instrument->counter_P2[26] = 0;
  instrument->counter_AbsorbProp[26]= 0;
  return(0);
} /* _Lmonfast2_setpos */

/* component Lmonfast2_zoom=L_monitor() SETTING, POSITION/ROTATION */
int _Lmonfast2_zoom_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Lmonfast2_zoom_setpos] component Lmonfast2_zoom=L_monitor() SETTING [/usr/share/mcstas/3.0-dev/monitors/L_monitor.comp:65]");
  stracpy(_Lmonfast2_zoom->_name, "Lmonfast2_zoom", 16384);
  stracpy(_Lmonfast2_zoom->_type, "L_monitor", 16384);
  _Lmonfast2_zoom->_index=27;
  _Lmonfast2_zoom->_parameters.nL = 200;
  #define nL (_Lmonfast2_zoom->_parameters.nL)
  if("Lmonfast2_zoom.dat" && strlen("Lmonfast2_zoom.dat"))
    stracpy(_Lmonfast2_zoom->_parameters.filename, "Lmonfast2_zoom.dat" ? "Lmonfast2_zoom.dat" : "", 16384);
  else 
  _Lmonfast2_zoom->_parameters.filename[0]='\0';
  #define filename (_Lmonfast2_zoom->_parameters.filename)
  _Lmonfast2_zoom->_parameters.xmin = -0.05;
  #define xmin (_Lmonfast2_zoom->_parameters.xmin)
  _Lmonfast2_zoom->_parameters.xmax = 0.05;
  #define xmax (_Lmonfast2_zoom->_parameters.xmax)
  _Lmonfast2_zoom->_parameters.ymin = -0.05;
  #define ymin (_Lmonfast2_zoom->_parameters.ymin)
  _Lmonfast2_zoom->_parameters.ymax = 0.05;
  #define ymax (_Lmonfast2_zoom->_parameters.ymax)
  _Lmonfast2_zoom->_parameters.xwidth = 0.06;
  #define xwidth (_Lmonfast2_zoom->_parameters.xwidth)
  _Lmonfast2_zoom->_parameters.yheight = 0.18;
  #define yheight (_Lmonfast2_zoom->_parameters.yheight)
  _Lmonfast2_zoom->_parameters.Lmin = instrument->_parameters._lambda0 -0.2;
  #define Lmin (_Lmonfast2_zoom->_parameters.Lmin)
  _Lmonfast2_zoom->_parameters.Lmax = instrument->_parameters._lambda0 + 0.2;
  #define Lmax (_Lmonfast2_zoom->_parameters.Lmax)
  _Lmonfast2_zoom->_parameters.restore_neutron = 0;
  #define restore_neutron (_Lmonfast2_zoom->_parameters.restore_neutron)

  #define L_N (_Lmonfast2_zoom->_parameters.L_N)
  #define L_p (_Lmonfast2_zoom->_parameters.L_p)
  #define L_p2 (_Lmonfast2_zoom->_parameters.L_p2)

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
  /* component Lmonfast2_zoom=L_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Lmonfast2->_rotation_absolute, _Lmonfast2_zoom->_rotation_absolute);
    rot_transpose(_Lmonfast2->_rotation_absolute, tr1);
    rot_mul(_Lmonfast2_zoom->_rotation_absolute, tr1, _Lmonfast2_zoom->_rotation_relative);
    _Lmonfast2_zoom->_rotation_is_identity =  rot_test_identity(_Lmonfast2_zoom->_rotation_relative);
    tc1 = coords_set(
      0, 0, 1e-6);
    rot_transpose(_Lmonfast2->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Lmonfast2_zoom->_position_absolute = coords_add(_Lmonfast2->_position_absolute, tc2);
    tc1 = coords_sub(_Lmonfast2->_position_absolute, _Lmonfast2_zoom->_position_absolute);
    _Lmonfast2_zoom->_position_relative = rot_apply(_Lmonfast2_zoom->_rotation_absolute, tc1);
  } /* Lmonfast2_zoom=L_monitor() AT ROTATED */
  DEBUG_COMPONENT("Lmonfast2_zoom", _Lmonfast2_zoom->_position_absolute, _Lmonfast2_zoom->_rotation_absolute);
  instrument->_position_absolute[27] = _Lmonfast2_zoom->_position_absolute;
  instrument->_position_relative[27] = _Lmonfast2_zoom->_position_relative;
  instrument->counter_N[27]  = instrument->counter_P[27] = instrument->counter_P2[27] = 0;
  instrument->counter_AbsorbProp[27]= 0;
  return(0);
} /* _Lmonfast2_zoom_setpos */

/* component TOFLfast2=TOFLambda_monitor() SETTING, POSITION/ROTATION */
int _TOFLfast2_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_TOFLfast2_setpos] component TOFLfast2=TOFLambda_monitor() SETTING [/usr/share/mcstas/3.0-dev/monitors/TOFLambda_monitor.comp:68]");
  stracpy(_TOFLfast2->_name, "TOFLfast2", 16384);
  stracpy(_TOFLfast2->_type, "TOFLambda_monitor", 16384);
  _TOFLfast2->_index=28;
  _TOFLfast2->_parameters.nL = 200;
  #define nL (_TOFLfast2->_parameters.nL)
  _TOFLfast2->_parameters.nt = 200;
  #define nt (_TOFLfast2->_parameters.nt)
  _TOFLfast2->_parameters.tmin = 0;
  #define tmin (_TOFLfast2->_parameters.tmin)
  _TOFLfast2->_parameters.tmax = 3.0e5;
  #define tmax (_TOFLfast2->_parameters.tmax)
  if("TOFLfast2.dat" && strlen("TOFLfast2.dat"))
    stracpy(_TOFLfast2->_parameters.filename, "TOFLfast2.dat" ? "TOFLfast2.dat" : "", 16384);
  else 
  _TOFLfast2->_parameters.filename[0]='\0';
  #define filename (_TOFLfast2->_parameters.filename)
  _TOFLfast2->_parameters.xmin = -0.05;
  #define xmin (_TOFLfast2->_parameters.xmin)
  _TOFLfast2->_parameters.xmax = 0.05;
  #define xmax (_TOFLfast2->_parameters.xmax)
  _TOFLfast2->_parameters.ymin = -0.05;
  #define ymin (_TOFLfast2->_parameters.ymin)
  _TOFLfast2->_parameters.ymax = 0.05;
  #define ymax (_TOFLfast2->_parameters.ymax)
  _TOFLfast2->_parameters.xwidth = 0.05;
  #define xwidth (_TOFLfast2->_parameters.xwidth)
  _TOFLfast2->_parameters.yheight = 0.12;
  #define yheight (_TOFLfast2->_parameters.yheight)
  _TOFLfast2->_parameters.Lmin = 0;
  #define Lmin (_TOFLfast2->_parameters.Lmin)
  _TOFLfast2->_parameters.Lmax = instrument->_parameters._Lmax + 1;
  #define Lmax (_TOFLfast2->_parameters.Lmax)
  _TOFLfast2->_parameters.restore_neutron = 0;
  #define restore_neutron (_TOFLfast2->_parameters.restore_neutron)

  #define tt_0 (_TOFLfast2->_parameters.tt_0)
  #define tt_1 (_TOFLfast2->_parameters.tt_1)
  #define TOFL_N (_TOFLfast2->_parameters.TOFL_N)
  #define TOFL_p (_TOFLfast2->_parameters.TOFL_p)
  #define TOFL_p2 (_TOFLfast2->_parameters.TOFL_p2)

  #undef nL
  #undef nt
  #undef tmin
  #undef tmax
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
  #undef tt_0
  #undef tt_1
  #undef TOFL_N
  #undef TOFL_p
  #undef TOFL_p2
  /* component TOFLfast2=TOFLambda_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Lmonfast2_zoom->_rotation_absolute, _TOFLfast2->_rotation_absolute);
    rot_transpose(_Lmonfast2_zoom->_rotation_absolute, tr1);
    rot_mul(_TOFLfast2->_rotation_absolute, tr1, _TOFLfast2->_rotation_relative);
    _TOFLfast2->_rotation_is_identity =  rot_test_identity(_TOFLfast2->_rotation_relative);
    tc1 = coords_set(
      0, 0, 1e-6);
    rot_transpose(_Lmonfast2_zoom->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _TOFLfast2->_position_absolute = coords_add(_Lmonfast2_zoom->_position_absolute, tc2);
    tc1 = coords_sub(_Lmonfast2_zoom->_position_absolute, _TOFLfast2->_position_absolute);
    _TOFLfast2->_position_relative = rot_apply(_TOFLfast2->_rotation_absolute, tc1);
  } /* TOFLfast2=TOFLambda_monitor() AT ROTATED */
  DEBUG_COMPONENT("TOFLfast2", _TOFLfast2->_position_absolute, _TOFLfast2->_rotation_absolute);
  instrument->_position_absolute[28] = _TOFLfast2->_position_absolute;
  instrument->_position_relative[28] = _TOFLfast2->_position_relative;
  instrument->counter_N[28]  = instrument->counter_P[28] = instrument->counter_P2[28] = 0;
  instrument->counter_AbsorbProp[28]= 0;
  return(0);
} /* _TOFLfast2_setpos */

/* component TOFLfast2zoom=TOFLambda_monitor() SETTING, POSITION/ROTATION */
int _TOFLfast2zoom_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_TOFLfast2zoom_setpos] component TOFLfast2zoom=TOFLambda_monitor() SETTING [/usr/share/mcstas/3.0-dev/monitors/TOFLambda_monitor.comp:68]");
  stracpy(_TOFLfast2zoom->_name, "TOFLfast2zoom", 16384);
  stracpy(_TOFLfast2zoom->_type, "TOFLambda_monitor", 16384);
  _TOFLfast2zoom->_index=29;
  _TOFLfast2zoom->_parameters.nL = 200;
  #define nL (_TOFLfast2zoom->_parameters.nL)
  _TOFLfast2zoom->_parameters.nt = 200;
  #define nt (_TOFLfast2zoom->_parameters.nt)
  _TOFLfast2zoom->_parameters.tmin = tmin_zoom;
  #define tmin (_TOFLfast2zoom->_parameters.tmin)
  _TOFLfast2zoom->_parameters.tmax = tmax_zoom;
  #define tmax (_TOFLfast2zoom->_parameters.tmax)
  if("TOFLfast2_zoom.dat" && strlen("TOFLfast2_zoom.dat"))
    stracpy(_TOFLfast2zoom->_parameters.filename, "TOFLfast2_zoom.dat" ? "TOFLfast2_zoom.dat" : "", 16384);
  else 
  _TOFLfast2zoom->_parameters.filename[0]='\0';
  #define filename (_TOFLfast2zoom->_parameters.filename)
  _TOFLfast2zoom->_parameters.xmin = -0.05;
  #define xmin (_TOFLfast2zoom->_parameters.xmin)
  _TOFLfast2zoom->_parameters.xmax = 0.05;
  #define xmax (_TOFLfast2zoom->_parameters.xmax)
  _TOFLfast2zoom->_parameters.ymin = -0.05;
  #define ymin (_TOFLfast2zoom->_parameters.ymin)
  _TOFLfast2zoom->_parameters.ymax = 0.05;
  #define ymax (_TOFLfast2zoom->_parameters.ymax)
  _TOFLfast2zoom->_parameters.xwidth = 0.05;
  #define xwidth (_TOFLfast2zoom->_parameters.xwidth)
  _TOFLfast2zoom->_parameters.yheight = 0.12;
  #define yheight (_TOFLfast2zoom->_parameters.yheight)
  _TOFLfast2zoom->_parameters.Lmin = instrument->_parameters._lambda0 -0.2;
  #define Lmin (_TOFLfast2zoom->_parameters.Lmin)
  _TOFLfast2zoom->_parameters.Lmax = instrument->_parameters._lambda0 + 0.2;
  #define Lmax (_TOFLfast2zoom->_parameters.Lmax)
  _TOFLfast2zoom->_parameters.restore_neutron = 0;
  #define restore_neutron (_TOFLfast2zoom->_parameters.restore_neutron)

  #define tt_0 (_TOFLfast2zoom->_parameters.tt_0)
  #define tt_1 (_TOFLfast2zoom->_parameters.tt_1)
  #define TOFL_N (_TOFLfast2zoom->_parameters.TOFL_N)
  #define TOFL_p (_TOFLfast2zoom->_parameters.TOFL_p)
  #define TOFL_p2 (_TOFLfast2zoom->_parameters.TOFL_p2)

  #undef nL
  #undef nt
  #undef tmin
  #undef tmax
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
  #undef tt_0
  #undef tt_1
  #undef TOFL_N
  #undef TOFL_p
  #undef TOFL_p2
  /* component TOFLfast2zoom=TOFLambda_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _TOFLfast2->_rotation_absolute, _TOFLfast2zoom->_rotation_absolute);
    rot_transpose(_TOFLfast2->_rotation_absolute, tr1);
    rot_mul(_TOFLfast2zoom->_rotation_absolute, tr1, _TOFLfast2zoom->_rotation_relative);
    _TOFLfast2zoom->_rotation_is_identity =  rot_test_identity(_TOFLfast2zoom->_rotation_relative);
    tc1 = coords_set(
      0, 0, 1e-6);
    rot_transpose(_TOFLfast2->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _TOFLfast2zoom->_position_absolute = coords_add(_TOFLfast2->_position_absolute, tc2);
    tc1 = coords_sub(_TOFLfast2->_position_absolute, _TOFLfast2zoom->_position_absolute);
    _TOFLfast2zoom->_position_relative = rot_apply(_TOFLfast2zoom->_rotation_absolute, tc1);
  } /* TOFLfast2zoom=TOFLambda_monitor() AT ROTATED */
  DEBUG_COMPONENT("TOFLfast2zoom", _TOFLfast2zoom->_position_absolute, _TOFLfast2zoom->_rotation_absolute);
  instrument->_position_absolute[29] = _TOFLfast2zoom->_position_absolute;
  instrument->_position_relative[29] = _TOFLfast2zoom->_position_relative;
  instrument->counter_N[29]  = instrument->counter_P[29] = instrument->counter_P2[29] = 0;
  instrument->counter_AbsorbProp[29]= 0;
  return(0);
} /* _TOFLfast2zoom_setpos */

/* component PSDfast2=PSD_monitor() SETTING, POSITION/ROTATION */
int _PSDfast2_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_PSDfast2_setpos] component PSDfast2=PSD_monitor() SETTING [/usr/share/mcstas/3.0-dev/monitors/PSD_monitor.comp:64]");
  stracpy(_PSDfast2->_name, "PSDfast2", 16384);
  stracpy(_PSDfast2->_type, "PSD_monitor", 16384);
  _PSDfast2->_index=30;
  _PSDfast2->_parameters.nx = 90;
  #define nx (_PSDfast2->_parameters.nx)
  _PSDfast2->_parameters.ny = 90;
  #define ny (_PSDfast2->_parameters.ny)
  if("PSDfast2.dat" && strlen("PSDfast2.dat"))
    stracpy(_PSDfast2->_parameters.filename, "PSDfast2.dat" ? "PSDfast2.dat" : "", 16384);
  else 
  _PSDfast2->_parameters.filename[0]='\0';
  #define filename (_PSDfast2->_parameters.filename)
  _PSDfast2->_parameters.xmin = -0.05;
  #define xmin (_PSDfast2->_parameters.xmin)
  _PSDfast2->_parameters.xmax = 0.05;
  #define xmax (_PSDfast2->_parameters.xmax)
  _PSDfast2->_parameters.ymin = -0.05;
  #define ymin (_PSDfast2->_parameters.ymin)
  _PSDfast2->_parameters.ymax = 0.05;
  #define ymax (_PSDfast2->_parameters.ymax)
  _PSDfast2->_parameters.xwidth = 0.1;
  #define xwidth (_PSDfast2->_parameters.xwidth)
  _PSDfast2->_parameters.yheight = 0.25;
  #define yheight (_PSDfast2->_parameters.yheight)
  _PSDfast2->_parameters.restore_neutron = 0;
  #define restore_neutron (_PSDfast2->_parameters.restore_neutron)

  #define PSD_N (_PSDfast2->_parameters.PSD_N)
  #define PSD_p (_PSDfast2->_parameters.PSD_p)
  #define PSD_p2 (_PSDfast2->_parameters.PSD_p2)

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
  /* component PSDfast2=PSD_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _TOFLfast2zoom->_rotation_absolute, _PSDfast2->_rotation_absolute);
    rot_transpose(_TOFLfast2zoom->_rotation_absolute, tr1);
    rot_mul(_PSDfast2->_rotation_absolute, tr1, _PSDfast2->_rotation_relative);
    _PSDfast2->_rotation_is_identity =  rot_test_identity(_PSDfast2->_rotation_relative);
    tc1 = coords_set(
      0, 0, 1e-6);
    rot_transpose(_TOFLfast2zoom->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _PSDfast2->_position_absolute = coords_add(_TOFLfast2zoom->_position_absolute, tc2);
    tc1 = coords_sub(_TOFLfast2zoom->_position_absolute, _PSDfast2->_position_absolute);
    _PSDfast2->_position_relative = rot_apply(_PSDfast2->_rotation_absolute, tc1);
  } /* PSDfast2=PSD_monitor() AT ROTATED */
  DEBUG_COMPONENT("PSDfast2", _PSDfast2->_position_absolute, _PSDfast2->_rotation_absolute);
  instrument->_position_absolute[30] = _PSDfast2->_position_absolute;
  instrument->_position_relative[30] = _PSDfast2->_position_relative;
  instrument->counter_N[30]  = instrument->counter_P[30] = instrument->counter_P2[30] = 0;
  instrument->counter_AbsorbProp[30]= 0;
  return(0);
} /* _PSDfast2_setpos */

/* component Fastchop2=DiskChopper() SETTING, POSITION/ROTATION */
int _Fastchop2_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Fastchop2_setpos] component Fastchop2=DiskChopper() SETTING [/usr/share/mcstas/3.0-dev/optics/DiskChopper.comp:67]");
  stracpy(_Fastchop2->_name, "Fastchop2", 16384);
  stracpy(_Fastchop2->_type, "DiskChopper", 16384);
  _Fastchop2->_index=31;
  _Fastchop2->_parameters.theta_0 = instrument->_parameters._FAST_THETA;
  #define theta_0 (_Fastchop2->_parameters.theta_0)
  _Fastchop2->_parameters.radius = 0.35;
  #define radius (_Fastchop2->_parameters.radius)
  _Fastchop2->_parameters.yheight = 0.35;
  #define yheight (_Fastchop2->_parameters.yheight)
  _Fastchop2->_parameters.nu = instrument->_parameters._F_fast2;
  #define nu (_Fastchop2->_parameters.nu)
  _Fastchop2->_parameters.nslit = instrument->_parameters._N_fast;
  #define nslit (_Fastchop2->_parameters.nslit)
  _Fastchop2->_parameters.jitter = 0;
  #define jitter (_Fastchop2->_parameters.jitter)
  _Fastchop2->_parameters.delay = t_fast2;
  #define delay (_Fastchop2->_parameters.delay)
  _Fastchop2->_parameters.isfirst = 0;
  #define isfirst (_Fastchop2->_parameters.isfirst)
  _Fastchop2->_parameters.n_pulse = 1;
  #define n_pulse (_Fastchop2->_parameters.n_pulse)
  _Fastchop2->_parameters.abs_out = 1;
  #define abs_out (_Fastchop2->_parameters.abs_out)
  _Fastchop2->_parameters.phase = 0;
  #define phase (_Fastchop2->_parameters.phase)
  _Fastchop2->_parameters.xwidth = 0;
  #define xwidth (_Fastchop2->_parameters.xwidth)
  _Fastchop2->_parameters.verbose = 0;
  #define verbose (_Fastchop2->_parameters.verbose)

  #define Tg (_Fastchop2->_parameters.Tg)
  #define To (_Fastchop2->_parameters.To)
  #define delta_y (_Fastchop2->_parameters.delta_y)
  #define height (_Fastchop2->_parameters.height)
  #define omega (_Fastchop2->_parameters.omega)

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
  /* component Fastchop2=DiskChopper() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Origin->_rotation_absolute, _Fastchop2->_rotation_absolute);
    rot_transpose(_PSDfast2->_rotation_absolute, tr1);
    rot_mul(_Fastchop2->_rotation_absolute, tr1, _Fastchop2->_rotation_relative);
    _Fastchop2->_rotation_is_identity =  rot_test_identity(_Fastchop2->_rotation_relative);
    tc1 = coords_set(
      0, 0.04, instrument->_parameters._Length);
    rot_transpose(_Origin->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Fastchop2->_position_absolute = coords_add(_Origin->_position_absolute, tc2);
    tc1 = coords_sub(_PSDfast2->_position_absolute, _Fastchop2->_position_absolute);
    _Fastchop2->_position_relative = rot_apply(_Fastchop2->_rotation_absolute, tc1);
  } /* Fastchop2=DiskChopper() AT ROTATED */
  DEBUG_COMPONENT("Fastchop2", _Fastchop2->_position_absolute, _Fastchop2->_rotation_absolute);
  instrument->_position_absolute[31] = _Fastchop2->_position_absolute;
  instrument->_position_relative[31] = _Fastchop2->_position_relative;
  instrument->counter_N[31]  = instrument->counter_P[31] = instrument->counter_P2[31] = 0;
  instrument->counter_AbsorbProp[31]= 0;
  return(0);
} /* _Fastchop2_setpos */

/* component Fastchop2counter=DiskChopper() SETTING, POSITION/ROTATION */
int _Fastchop2counter_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Fastchop2counter_setpos] component Fastchop2counter=DiskChopper() SETTING [/usr/share/mcstas/3.0-dev/optics/DiskChopper.comp:67]");
  stracpy(_Fastchop2counter->_name, "Fastchop2counter", 16384);
  stracpy(_Fastchop2counter->_type, "DiskChopper", 16384);
  _Fastchop2counter->_index=32;
  _Fastchop2counter->_parameters.theta_0 = instrument->_parameters._FAST_THETA;
  #define theta_0 (_Fastchop2counter->_parameters.theta_0)
  _Fastchop2counter->_parameters.radius = 0.35;
  #define radius (_Fastchop2counter->_parameters.radius)
  _Fastchop2counter->_parameters.yheight = 0.35;
  #define yheight (_Fastchop2counter->_parameters.yheight)
  _Fastchop2counter->_parameters.nu = - instrument->_parameters._F_fast2;
  #define nu (_Fastchop2counter->_parameters.nu)
  _Fastchop2counter->_parameters.nslit = instrument->_parameters._N_fast;
  #define nslit (_Fastchop2counter->_parameters.nslit)
  _Fastchop2counter->_parameters.jitter = 0;
  #define jitter (_Fastchop2counter->_parameters.jitter)
  _Fastchop2counter->_parameters.delay = - t_fast2a;
  #define delay (_Fastchop2counter->_parameters.delay)
  _Fastchop2counter->_parameters.isfirst = 0;
  #define isfirst (_Fastchop2counter->_parameters.isfirst)
  _Fastchop2counter->_parameters.n_pulse = 1;
  #define n_pulse (_Fastchop2counter->_parameters.n_pulse)
  _Fastchop2counter->_parameters.abs_out = 1;
  #define abs_out (_Fastchop2counter->_parameters.abs_out)
  _Fastchop2counter->_parameters.phase = 0;
  #define phase (_Fastchop2counter->_parameters.phase)
  _Fastchop2counter->_parameters.xwidth = 0;
  #define xwidth (_Fastchop2counter->_parameters.xwidth)
  _Fastchop2counter->_parameters.verbose = 0;
  #define verbose (_Fastchop2counter->_parameters.verbose)

  #define Tg (_Fastchop2counter->_parameters.Tg)
  #define To (_Fastchop2counter->_parameters.To)
  #define delta_y (_Fastchop2counter->_parameters.delta_y)
  #define height (_Fastchop2counter->_parameters.height)
  #define omega (_Fastchop2counter->_parameters.omega)

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
  /* component Fastchop2counter=DiskChopper() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Origin->_rotation_absolute, _Fastchop2counter->_rotation_absolute);
    rot_transpose(_Fastchop2->_rotation_absolute, tr1);
    rot_mul(_Fastchop2counter->_rotation_absolute, tr1, _Fastchop2counter->_rotation_relative);
    _Fastchop2counter->_rotation_is_identity =  rot_test_identity(_Fastchop2counter->_rotation_relative);
    tc1 = coords_set(
      0, 0.04, instrument->_parameters._Length + instrument->_parameters._GUI_GAP);
    rot_transpose(_Origin->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Fastchop2counter->_position_absolute = coords_add(_Origin->_position_absolute, tc2);
    tc1 = coords_sub(_Fastchop2->_position_absolute, _Fastchop2counter->_position_absolute);
    _Fastchop2counter->_position_relative = rot_apply(_Fastchop2counter->_rotation_absolute, tc1);
  } /* Fastchop2counter=DiskChopper() AT ROTATED */
  DEBUG_COMPONENT("Fastchop2counter", _Fastchop2counter->_position_absolute, _Fastchop2counter->_rotation_absolute);
  instrument->_position_absolute[32] = _Fastchop2counter->_position_absolute;
  instrument->_position_relative[32] = _Fastchop2counter->_position_relative;
  instrument->counter_N[32]  = instrument->counter_P[32] = instrument->counter_P2[32] = 0;
  instrument->counter_AbsorbProp[32]= 0;
  return(0);
} /* _Fastchop2counter_setpos */

/* component FOchop3=DiskChopper() SETTING, POSITION/ROTATION */
int _FOchop3_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_FOchop3_setpos] component FOchop3=DiskChopper() SETTING [/usr/share/mcstas/3.0-dev/optics/DiskChopper.comp:67]");
  stracpy(_FOchop3->_name, "FOchop3", 16384);
  stracpy(_FOchop3->_type, "DiskChopper", 16384);
  _FOchop3->_index=33;
  _FOchop3->_parameters.theta_0 = 2 * instrument->_parameters._FAST_THETA;
  #define theta_0 (_FOchop3->_parameters.theta_0)
  _FOchop3->_parameters.radius = 0.35;
  #define radius (_FOchop3->_parameters.radius)
  _FOchop3->_parameters.yheight = 0.35;
  #define yheight (_FOchop3->_parameters.yheight)
  _FOchop3->_parameters.nu = instrument->_parameters._F_fast2 / instrument->_parameters._FO3;
  #define nu (_FOchop3->_parameters.nu)
  _FOchop3->_parameters.nslit = instrument->_parameters._N_fast;
  #define nslit (_FOchop3->_parameters.nslit)
  _FOchop3->_parameters.jitter = 0;
  #define jitter (_FOchop3->_parameters.jitter)
  _FOchop3->_parameters.delay = t_fast3;
  #define delay (_FOchop3->_parameters.delay)
  _FOchop3->_parameters.isfirst = 0;
  #define isfirst (_FOchop3->_parameters.isfirst)
  _FOchop3->_parameters.n_pulse = 1;
  #define n_pulse (_FOchop3->_parameters.n_pulse)
  _FOchop3->_parameters.abs_out = 1;
  #define abs_out (_FOchop3->_parameters.abs_out)
  _FOchop3->_parameters.phase = 0;
  #define phase (_FOchop3->_parameters.phase)
  _FOchop3->_parameters.xwidth = 0;
  #define xwidth (_FOchop3->_parameters.xwidth)
  _FOchop3->_parameters.verbose = 0;
  #define verbose (_FOchop3->_parameters.verbose)

  #define Tg (_FOchop3->_parameters.Tg)
  #define To (_FOchop3->_parameters.To)
  #define delta_y (_FOchop3->_parameters.delta_y)
  #define height (_FOchop3->_parameters.height)
  #define omega (_FOchop3->_parameters.omega)

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
  /* component FOchop3=DiskChopper() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Origin->_rotation_absolute, _FOchop3->_rotation_absolute);
    rot_transpose(_Fastchop2counter->_rotation_absolute, tr1);
    rot_mul(_FOchop3->_rotation_absolute, tr1, _FOchop3->_rotation_relative);
    _FOchop3->_rotation_is_identity =  rot_test_identity(_FOchop3->_rotation_relative);
    tc1 = coords_set(
      0, 0.04, instrument->_parameters._Length + 2 * instrument->_parameters._GUI_GAP);
    rot_transpose(_Origin->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _FOchop3->_position_absolute = coords_add(_Origin->_position_absolute, tc2);
    tc1 = coords_sub(_Fastchop2counter->_position_absolute, _FOchop3->_position_absolute);
    _FOchop3->_position_relative = rot_apply(_FOchop3->_rotation_absolute, tc1);
  } /* FOchop3=DiskChopper() AT ROTATED */
  DEBUG_COMPONENT("FOchop3", _FOchop3->_position_absolute, _FOchop3->_rotation_absolute);
  instrument->_position_absolute[33] = _FOchop3->_position_absolute;
  instrument->_position_relative[33] = _FOchop3->_position_relative;
  instrument->counter_N[33]  = instrument->counter_P[33] = instrument->counter_P2[33] = 0;
  instrument->counter_AbsorbProp[33]= 0;
  return(0);
} /* _FOchop3_setpos */

/* component TOFfast2_zoom=TOF_monitor() SETTING, POSITION/ROTATION */
int _TOFfast2_zoom_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_TOFfast2_zoom_setpos] component TOFfast2_zoom=TOF_monitor() SETTING [/usr/share/mcstas/3.0-dev/monitors/TOF_monitor.comp:63]");
  stracpy(_TOFfast2_zoom->_name, "TOFfast2_zoom", 16384);
  stracpy(_TOFfast2_zoom->_type, "TOF_monitor", 16384);
  _TOFfast2_zoom->_index=34;
  _TOFfast2_zoom->_parameters.nt = 100;
  #define nt (_TOFfast2_zoom->_parameters.nt)
  if("TOF_fast2.dat" && strlen("TOF_fast2.dat"))
    stracpy(_TOFfast2_zoom->_parameters.filename, "TOF_fast2.dat" ? "TOF_fast2.dat" : "", 16384);
  else 
  _TOFfast2_zoom->_parameters.filename[0]='\0';
  #define filename (_TOFfast2_zoom->_parameters.filename)
  _TOFfast2_zoom->_parameters.xmin = -0.05;
  #define xmin (_TOFfast2_zoom->_parameters.xmin)
  _TOFfast2_zoom->_parameters.xmax = 0.05;
  #define xmax (_TOFfast2_zoom->_parameters.xmax)
  _TOFfast2_zoom->_parameters.ymin = -0.05;
  #define ymin (_TOFfast2_zoom->_parameters.ymin)
  _TOFfast2_zoom->_parameters.ymax = 0.05;
  #define ymax (_TOFfast2_zoom->_parameters.ymax)
  _TOFfast2_zoom->_parameters.xwidth = 1.1;
  #define xwidth (_TOFfast2_zoom->_parameters.xwidth)
  _TOFfast2_zoom->_parameters.yheight = 1.2;
  #define yheight (_TOFfast2_zoom->_parameters.yheight)
  _TOFfast2_zoom->_parameters.tmin = 1e6 * t_fast3 -2e2;
  #define tmin (_TOFfast2_zoom->_parameters.tmin)
  _TOFfast2_zoom->_parameters.tmax = 1e6 * t_fast3 + 2e2;
  #define tmax (_TOFfast2_zoom->_parameters.tmax)
  _TOFfast2_zoom->_parameters.dt = 1.0;
  #define dt (_TOFfast2_zoom->_parameters.dt)
  _TOFfast2_zoom->_parameters.restore_neutron = 0;
  #define restore_neutron (_TOFfast2_zoom->_parameters.restore_neutron)

  #define TOF_N (_TOFfast2_zoom->_parameters.TOF_N)
  #define TOF_p (_TOFfast2_zoom->_parameters.TOF_p)
  #define TOF_p2 (_TOFfast2_zoom->_parameters.TOF_p2)
  #define t_min (_TOFfast2_zoom->_parameters.t_min)
  #define t_max (_TOFfast2_zoom->_parameters.t_max)
  #define delta_t (_TOFfast2_zoom->_parameters.delta_t)

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
  /* component TOFfast2_zoom=TOF_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _FOchop3->_rotation_absolute, _TOFfast2_zoom->_rotation_absolute);
    rot_transpose(_FOchop3->_rotation_absolute, tr1);
    rot_mul(_TOFfast2_zoom->_rotation_absolute, tr1, _TOFfast2_zoom->_rotation_relative);
    _TOFfast2_zoom->_rotation_is_identity =  rot_test_identity(_TOFfast2_zoom->_rotation_relative);
    tc1 = coords_set(
      0, -0.04, 1e-6);
    rot_transpose(_FOchop3->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _TOFfast2_zoom->_position_absolute = coords_add(_FOchop3->_position_absolute, tc2);
    tc1 = coords_sub(_FOchop3->_position_absolute, _TOFfast2_zoom->_position_absolute);
    _TOFfast2_zoom->_position_relative = rot_apply(_TOFfast2_zoom->_rotation_absolute, tc1);
  } /* TOFfast2_zoom=TOF_monitor() AT ROTATED */
  DEBUG_COMPONENT("TOFfast2_zoom", _TOFfast2_zoom->_position_absolute, _TOFfast2_zoom->_rotation_absolute);
  instrument->_position_absolute[34] = _TOFfast2_zoom->_position_absolute;
  instrument->_position_relative[34] = _TOFfast2_zoom->_position_relative;
  instrument->counter_N[34]  = instrument->counter_P[34] = instrument->counter_P2[34] = 0;
  instrument->counter_AbsorbProp[34]= 0;
  return(0);
} /* _TOFfast2_zoom_setpos */

/* component Lmon_afterfast2=L_monitor() SETTING, POSITION/ROTATION */
int _Lmon_afterfast2_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Lmon_afterfast2_setpos] component Lmon_afterfast2=L_monitor() SETTING [/usr/share/mcstas/3.0-dev/monitors/L_monitor.comp:65]");
  stracpy(_Lmon_afterfast2->_name, "Lmon_afterfast2", 16384);
  stracpy(_Lmon_afterfast2->_type, "L_monitor", 16384);
  _Lmon_afterfast2->_index=35;
  _Lmon_afterfast2->_parameters.nL = 500;
  #define nL (_Lmon_afterfast2->_parameters.nL)
  if("Lmon_afterfast2.dat" && strlen("Lmon_afterfast2.dat"))
    stracpy(_Lmon_afterfast2->_parameters.filename, "Lmon_afterfast2.dat" ? "Lmon_afterfast2.dat" : "", 16384);
  else 
  _Lmon_afterfast2->_parameters.filename[0]='\0';
  #define filename (_Lmon_afterfast2->_parameters.filename)
  _Lmon_afterfast2->_parameters.xmin = -0.05;
  #define xmin (_Lmon_afterfast2->_parameters.xmin)
  _Lmon_afterfast2->_parameters.xmax = 0.05;
  #define xmax (_Lmon_afterfast2->_parameters.xmax)
  _Lmon_afterfast2->_parameters.ymin = -0.05;
  #define ymin (_Lmon_afterfast2->_parameters.ymin)
  _Lmon_afterfast2->_parameters.ymax = 0.05;
  #define ymax (_Lmon_afterfast2->_parameters.ymax)
  _Lmon_afterfast2->_parameters.xwidth = 0.06;
  #define xwidth (_Lmon_afterfast2->_parameters.xwidth)
  _Lmon_afterfast2->_parameters.yheight = 0.18;
  #define yheight (_Lmon_afterfast2->_parameters.yheight)
  _Lmon_afterfast2->_parameters.Lmin = 0;
  #define Lmin (_Lmon_afterfast2->_parameters.Lmin)
  _Lmon_afterfast2->_parameters.Lmax = instrument->_parameters._Lmax + 1;
  #define Lmax (_Lmon_afterfast2->_parameters.Lmax)
  _Lmon_afterfast2->_parameters.restore_neutron = 0;
  #define restore_neutron (_Lmon_afterfast2->_parameters.restore_neutron)

  #define L_N (_Lmon_afterfast2->_parameters.L_N)
  #define L_p (_Lmon_afterfast2->_parameters.L_p)
  #define L_p2 (_Lmon_afterfast2->_parameters.L_p2)

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
  /* component Lmon_afterfast2=L_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _TOFfast2_zoom->_rotation_absolute, _Lmon_afterfast2->_rotation_absolute);
    rot_transpose(_TOFfast2_zoom->_rotation_absolute, tr1);
    rot_mul(_Lmon_afterfast2->_rotation_absolute, tr1, _Lmon_afterfast2->_rotation_relative);
    _Lmon_afterfast2->_rotation_is_identity =  rot_test_identity(_Lmon_afterfast2->_rotation_relative);
    tc1 = coords_set(
      0, 0, 1e-6);
    rot_transpose(_TOFfast2_zoom->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Lmon_afterfast2->_position_absolute = coords_add(_TOFfast2_zoom->_position_absolute, tc2);
    tc1 = coords_sub(_TOFfast2_zoom->_position_absolute, _Lmon_afterfast2->_position_absolute);
    _Lmon_afterfast2->_position_relative = rot_apply(_Lmon_afterfast2->_rotation_absolute, tc1);
  } /* Lmon_afterfast2=L_monitor() AT ROTATED */
  DEBUG_COMPONENT("Lmon_afterfast2", _Lmon_afterfast2->_position_absolute, _Lmon_afterfast2->_rotation_absolute);
  instrument->_position_absolute[35] = _Lmon_afterfast2->_position_absolute;
  instrument->_position_relative[35] = _Lmon_afterfast2->_position_relative;
  instrument->counter_N[35]  = instrument->counter_P[35] = instrument->counter_P2[35] = 0;
  instrument->counter_AbsorbProp[35]= 0;
  return(0);
} /* _Lmon_afterfast2_setpos */

/* component TOFL_afterfast2=TOFLambda_monitor() SETTING, POSITION/ROTATION */
int _TOFL_afterfast2_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_TOFL_afterfast2_setpos] component TOFL_afterfast2=TOFLambda_monitor() SETTING [/usr/share/mcstas/3.0-dev/monitors/TOFLambda_monitor.comp:68]");
  stracpy(_TOFL_afterfast2->_name, "TOFL_afterfast2", 16384);
  stracpy(_TOFL_afterfast2->_type, "TOFLambda_monitor", 16384);
  _TOFL_afterfast2->_index=36;
  _TOFL_afterfast2->_parameters.nL = 200;
  #define nL (_TOFL_afterfast2->_parameters.nL)
  _TOFL_afterfast2->_parameters.nt = 200;
  #define nt (_TOFL_afterfast2->_parameters.nt)
  _TOFL_afterfast2->_parameters.tmin = 0;
  #define tmin (_TOFL_afterfast2->_parameters.tmin)
  _TOFL_afterfast2->_parameters.tmax = 3.0e5;
  #define tmax (_TOFL_afterfast2->_parameters.tmax)
  if("TOF_afterfast2.dat" && strlen("TOF_afterfast2.dat"))
    stracpy(_TOFL_afterfast2->_parameters.filename, "TOF_afterfast2.dat" ? "TOF_afterfast2.dat" : "", 16384);
  else 
  _TOFL_afterfast2->_parameters.filename[0]='\0';
  #define filename (_TOFL_afterfast2->_parameters.filename)
  _TOFL_afterfast2->_parameters.xmin = -0.05;
  #define xmin (_TOFL_afterfast2->_parameters.xmin)
  _TOFL_afterfast2->_parameters.xmax = 0.05;
  #define xmax (_TOFL_afterfast2->_parameters.xmax)
  _TOFL_afterfast2->_parameters.ymin = -0.05;
  #define ymin (_TOFL_afterfast2->_parameters.ymin)
  _TOFL_afterfast2->_parameters.ymax = 0.05;
  #define ymax (_TOFL_afterfast2->_parameters.ymax)
  _TOFL_afterfast2->_parameters.xwidth = 0.05;
  #define xwidth (_TOFL_afterfast2->_parameters.xwidth)
  _TOFL_afterfast2->_parameters.yheight = 0.12;
  #define yheight (_TOFL_afterfast2->_parameters.yheight)
  _TOFL_afterfast2->_parameters.Lmin = 0;
  #define Lmin (_TOFL_afterfast2->_parameters.Lmin)
  _TOFL_afterfast2->_parameters.Lmax = instrument->_parameters._Lmax + 1;
  #define Lmax (_TOFL_afterfast2->_parameters.Lmax)
  _TOFL_afterfast2->_parameters.restore_neutron = 0;
  #define restore_neutron (_TOFL_afterfast2->_parameters.restore_neutron)

  #define tt_0 (_TOFL_afterfast2->_parameters.tt_0)
  #define tt_1 (_TOFL_afterfast2->_parameters.tt_1)
  #define TOFL_N (_TOFL_afterfast2->_parameters.TOFL_N)
  #define TOFL_p (_TOFL_afterfast2->_parameters.TOFL_p)
  #define TOFL_p2 (_TOFL_afterfast2->_parameters.TOFL_p2)

  #undef nL
  #undef nt
  #undef tmin
  #undef tmax
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
  #undef tt_0
  #undef tt_1
  #undef TOFL_N
  #undef TOFL_p
  #undef TOFL_p2
  /* component TOFL_afterfast2=TOFLambda_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Lmon_afterfast2->_rotation_absolute, _TOFL_afterfast2->_rotation_absolute);
    rot_transpose(_Lmon_afterfast2->_rotation_absolute, tr1);
    rot_mul(_TOFL_afterfast2->_rotation_absolute, tr1, _TOFL_afterfast2->_rotation_relative);
    _TOFL_afterfast2->_rotation_is_identity =  rot_test_identity(_TOFL_afterfast2->_rotation_relative);
    tc1 = coords_set(
      0, 0, 1e-6);
    rot_transpose(_Lmon_afterfast2->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _TOFL_afterfast2->_position_absolute = coords_add(_Lmon_afterfast2->_position_absolute, tc2);
    tc1 = coords_sub(_Lmon_afterfast2->_position_absolute, _TOFL_afterfast2->_position_absolute);
    _TOFL_afterfast2->_position_relative = rot_apply(_TOFL_afterfast2->_rotation_absolute, tc1);
  } /* TOFL_afterfast2=TOFLambda_monitor() AT ROTATED */
  DEBUG_COMPONENT("TOFL_afterfast2", _TOFL_afterfast2->_position_absolute, _TOFL_afterfast2->_rotation_absolute);
  instrument->_position_absolute[36] = _TOFL_afterfast2->_position_absolute;
  instrument->_position_relative[36] = _TOFL_afterfast2->_position_relative;
  instrument->counter_N[36]  = instrument->counter_P[36] = instrument->counter_P2[36] = 0;
  instrument->counter_AbsorbProp[36]= 0;
  return(0);
} /* _TOFL_afterfast2_setpos */

/* component TOFL_afterfast2_zoom=TOFLambda_monitor() SETTING, POSITION/ROTATION */
int _TOFL_afterfast2_zoom_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_TOFL_afterfast2_zoom_setpos] component TOFL_afterfast2_zoom=TOFLambda_monitor() SETTING [/usr/share/mcstas/3.0-dev/monitors/TOFLambda_monitor.comp:68]");
  stracpy(_TOFL_afterfast2_zoom->_name, "TOFL_afterfast2_zoom", 16384);
  stracpy(_TOFL_afterfast2_zoom->_type, "TOFLambda_monitor", 16384);
  _TOFL_afterfast2_zoom->_index=37;
  _TOFL_afterfast2_zoom->_parameters.nL = 200;
  #define nL (_TOFL_afterfast2_zoom->_parameters.nL)
  _TOFL_afterfast2_zoom->_parameters.nt = 200;
  #define nt (_TOFL_afterfast2_zoom->_parameters.nt)
  _TOFL_afterfast2_zoom->_parameters.tmin = tmin_zoom;
  #define tmin (_TOFL_afterfast2_zoom->_parameters.tmin)
  _TOFL_afterfast2_zoom->_parameters.tmax = tmax_zoom;
  #define tmax (_TOFL_afterfast2_zoom->_parameters.tmax)
  if("TOFL_afterfast2_zoom.dat" && strlen("TOFL_afterfast2_zoom.dat"))
    stracpy(_TOFL_afterfast2_zoom->_parameters.filename, "TOFL_afterfast2_zoom.dat" ? "TOFL_afterfast2_zoom.dat" : "", 16384);
  else 
  _TOFL_afterfast2_zoom->_parameters.filename[0]='\0';
  #define filename (_TOFL_afterfast2_zoom->_parameters.filename)
  _TOFL_afterfast2_zoom->_parameters.xmin = -0.05;
  #define xmin (_TOFL_afterfast2_zoom->_parameters.xmin)
  _TOFL_afterfast2_zoom->_parameters.xmax = 0.05;
  #define xmax (_TOFL_afterfast2_zoom->_parameters.xmax)
  _TOFL_afterfast2_zoom->_parameters.ymin = -0.05;
  #define ymin (_TOFL_afterfast2_zoom->_parameters.ymin)
  _TOFL_afterfast2_zoom->_parameters.ymax = 0.05;
  #define ymax (_TOFL_afterfast2_zoom->_parameters.ymax)
  _TOFL_afterfast2_zoom->_parameters.xwidth = 0.05;
  #define xwidth (_TOFL_afterfast2_zoom->_parameters.xwidth)
  _TOFL_afterfast2_zoom->_parameters.yheight = 0.12;
  #define yheight (_TOFL_afterfast2_zoom->_parameters.yheight)
  _TOFL_afterfast2_zoom->_parameters.Lmin = instrument->_parameters._lambda0 -0.2;
  #define Lmin (_TOFL_afterfast2_zoom->_parameters.Lmin)
  _TOFL_afterfast2_zoom->_parameters.Lmax = instrument->_parameters._lambda0 + 0.2;
  #define Lmax (_TOFL_afterfast2_zoom->_parameters.Lmax)
  _TOFL_afterfast2_zoom->_parameters.restore_neutron = 0;
  #define restore_neutron (_TOFL_afterfast2_zoom->_parameters.restore_neutron)

  #define tt_0 (_TOFL_afterfast2_zoom->_parameters.tt_0)
  #define tt_1 (_TOFL_afterfast2_zoom->_parameters.tt_1)
  #define TOFL_N (_TOFL_afterfast2_zoom->_parameters.TOFL_N)
  #define TOFL_p (_TOFL_afterfast2_zoom->_parameters.TOFL_p)
  #define TOFL_p2 (_TOFL_afterfast2_zoom->_parameters.TOFL_p2)

  #undef nL
  #undef nt
  #undef tmin
  #undef tmax
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
  #undef tt_0
  #undef tt_1
  #undef TOFL_N
  #undef TOFL_p
  #undef TOFL_p2
  /* component TOFL_afterfast2_zoom=TOFLambda_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _TOFL_afterfast2->_rotation_absolute, _TOFL_afterfast2_zoom->_rotation_absolute);
    rot_transpose(_TOFL_afterfast2->_rotation_absolute, tr1);
    rot_mul(_TOFL_afterfast2_zoom->_rotation_absolute, tr1, _TOFL_afterfast2_zoom->_rotation_relative);
    _TOFL_afterfast2_zoom->_rotation_is_identity =  rot_test_identity(_TOFL_afterfast2_zoom->_rotation_relative);
    tc1 = coords_set(
      0, 0, 1e-6);
    rot_transpose(_TOFL_afterfast2->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _TOFL_afterfast2_zoom->_position_absolute = coords_add(_TOFL_afterfast2->_position_absolute, tc2);
    tc1 = coords_sub(_TOFL_afterfast2->_position_absolute, _TOFL_afterfast2_zoom->_position_absolute);
    _TOFL_afterfast2_zoom->_position_relative = rot_apply(_TOFL_afterfast2_zoom->_rotation_absolute, tc1);
  } /* TOFL_afterfast2_zoom=TOFLambda_monitor() AT ROTATED */
  DEBUG_COMPONENT("TOFL_afterfast2_zoom", _TOFL_afterfast2_zoom->_position_absolute, _TOFL_afterfast2_zoom->_rotation_absolute);
  instrument->_position_absolute[37] = _TOFL_afterfast2_zoom->_position_absolute;
  instrument->_position_relative[37] = _TOFL_afterfast2_zoom->_position_relative;
  instrument->counter_N[37]  = instrument->counter_P[37] = instrument->counter_P2[37] = 0;
  instrument->counter_AbsorbProp[37]= 0;
  return(0);
} /* _TOFL_afterfast2_zoom_setpos */

/* component PSD_afterfast2=PSD_monitor() SETTING, POSITION/ROTATION */
int _PSD_afterfast2_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_PSD_afterfast2_setpos] component PSD_afterfast2=PSD_monitor() SETTING [/usr/share/mcstas/3.0-dev/monitors/PSD_monitor.comp:64]");
  stracpy(_PSD_afterfast2->_name, "PSD_afterfast2", 16384);
  stracpy(_PSD_afterfast2->_type, "PSD_monitor", 16384);
  _PSD_afterfast2->_index=38;
  _PSD_afterfast2->_parameters.nx = 90;
  #define nx (_PSD_afterfast2->_parameters.nx)
  _PSD_afterfast2->_parameters.ny = 90;
  #define ny (_PSD_afterfast2->_parameters.ny)
  if("PSD_afterfast2.dat" && strlen("PSD_afterfast2.dat"))
    stracpy(_PSD_afterfast2->_parameters.filename, "PSD_afterfast2.dat" ? "PSD_afterfast2.dat" : "", 16384);
  else 
  _PSD_afterfast2->_parameters.filename[0]='\0';
  #define filename (_PSD_afterfast2->_parameters.filename)
  _PSD_afterfast2->_parameters.xmin = -0.05;
  #define xmin (_PSD_afterfast2->_parameters.xmin)
  _PSD_afterfast2->_parameters.xmax = 0.05;
  #define xmax (_PSD_afterfast2->_parameters.xmax)
  _PSD_afterfast2->_parameters.ymin = -0.05;
  #define ymin (_PSD_afterfast2->_parameters.ymin)
  _PSD_afterfast2->_parameters.ymax = 0.05;
  #define ymax (_PSD_afterfast2->_parameters.ymax)
  _PSD_afterfast2->_parameters.xwidth = 0.1;
  #define xwidth (_PSD_afterfast2->_parameters.xwidth)
  _PSD_afterfast2->_parameters.yheight = 0.25;
  #define yheight (_PSD_afterfast2->_parameters.yheight)
  _PSD_afterfast2->_parameters.restore_neutron = 0;
  #define restore_neutron (_PSD_afterfast2->_parameters.restore_neutron)

  #define PSD_N (_PSD_afterfast2->_parameters.PSD_N)
  #define PSD_p (_PSD_afterfast2->_parameters.PSD_p)
  #define PSD_p2 (_PSD_afterfast2->_parameters.PSD_p2)

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
  /* component PSD_afterfast2=PSD_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _TOFL_afterfast2_zoom->_rotation_absolute, _PSD_afterfast2->_rotation_absolute);
    rot_transpose(_TOFL_afterfast2_zoom->_rotation_absolute, tr1);
    rot_mul(_PSD_afterfast2->_rotation_absolute, tr1, _PSD_afterfast2->_rotation_relative);
    _PSD_afterfast2->_rotation_is_identity =  rot_test_identity(_PSD_afterfast2->_rotation_relative);
    tc1 = coords_set(
      0, 0, 1e-6);
    rot_transpose(_TOFL_afterfast2_zoom->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _PSD_afterfast2->_position_absolute = coords_add(_TOFL_afterfast2_zoom->_position_absolute, tc2);
    tc1 = coords_sub(_TOFL_afterfast2_zoom->_position_absolute, _PSD_afterfast2->_position_absolute);
    _PSD_afterfast2->_position_relative = rot_apply(_PSD_afterfast2->_rotation_absolute, tc1);
  } /* PSD_afterfast2=PSD_monitor() AT ROTATED */
  DEBUG_COMPONENT("PSD_afterfast2", _PSD_afterfast2->_position_absolute, _PSD_afterfast2->_rotation_absolute);
  instrument->_position_absolute[38] = _PSD_afterfast2->_position_absolute;
  instrument->_position_relative[38] = _PSD_afterfast2->_position_relative;
  instrument->counter_N[38]  = instrument->counter_P[38] = instrument->counter_P2[38] = 0;
  instrument->counter_AbsorbProp[38]= 0;
  return(0);
} /* _PSD_afterfast2_setpos */

/* component Guidesample=Guide() SETTING, POSITION/ROTATION */
int _Guidesample_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Guidesample_setpos] component Guidesample=Guide() SETTING [/usr/share/mcstas/3.0-dev/optics/Guide.comp:73]");
  stracpy(_Guidesample->_name, "Guidesample", 16384);
  stracpy(_Guidesample->_type, "Guide", 16384);
  _Guidesample->_index=39;
  _Guidesample->_parameters.reflect[0]='\0';
  #define reflect (_Guidesample->_parameters.reflect)
  _Guidesample->_parameters.w1 = instrument->_parameters._W_chop;
  #define w1 (_Guidesample->_parameters.w1)
  _Guidesample->_parameters.h1 = instrument->_parameters._H_chop;
  #define h1 (_Guidesample->_parameters.h1)
  _Guidesample->_parameters.w2 = instrument->_parameters._W_end;
  #define w2 (_Guidesample->_parameters.w2)
  _Guidesample->_parameters.h2 = instrument->_parameters._H_end;
  #define h2 (_Guidesample->_parameters.h2)
  _Guidesample->_parameters.l = instrument->_parameters._SAMPLE_DIST -4 * instrument->_parameters._GUI_GAP;
  #define l (_Guidesample->_parameters.l)
  _Guidesample->_parameters.R0 = 1;
  #define R0 (_Guidesample->_parameters.R0)
  _Guidesample->_parameters.Qc = 0.0219;
  #define Qc (_Guidesample->_parameters.Qc)
  _Guidesample->_parameters.alpha = instrument->_parameters._ALPHA;
  #define alpha (_Guidesample->_parameters.alpha)
  _Guidesample->_parameters.m = instrument->_parameters._M;
  #define m (_Guidesample->_parameters.m)
  _Guidesample->_parameters.W = 0.003;
  #define W (_Guidesample->_parameters.W)

  #define pTable (_Guidesample->_parameters.pTable)

  #undef reflect
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
  #undef pTable
  /* component Guidesample=Guide() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _PSD_afterfast2->_rotation_absolute, _Guidesample->_rotation_absolute);
    rot_transpose(_PSD_afterfast2->_rotation_absolute, tr1);
    rot_mul(_Guidesample->_rotation_absolute, tr1, _Guidesample->_rotation_relative);
    _Guidesample->_rotation_is_identity =  rot_test_identity(_Guidesample->_rotation_relative);
    tc1 = coords_set(
      0, 0, instrument->_parameters._GUI_GAP / 2);
    rot_transpose(_PSD_afterfast2->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Guidesample->_position_absolute = coords_add(_PSD_afterfast2->_position_absolute, tc2);
    tc1 = coords_sub(_PSD_afterfast2->_position_absolute, _Guidesample->_position_absolute);
    _Guidesample->_position_relative = rot_apply(_Guidesample->_rotation_absolute, tc1);
  } /* Guidesample=Guide() AT ROTATED */
  DEBUG_COMPONENT("Guidesample", _Guidesample->_position_absolute, _Guidesample->_rotation_absolute);
  instrument->_position_absolute[39] = _Guidesample->_position_absolute;
  instrument->_position_relative[39] = _Guidesample->_position_relative;
  instrument->counter_N[39]  = instrument->counter_P[39] = instrument->counter_P2[39] = 0;
  instrument->counter_AbsorbProp[39]= 0;
  return(0);
} /* _Guidesample_setpos */

/* component Lmon_guideend=L_monitor() SETTING, POSITION/ROTATION */
int _Lmon_guideend_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Lmon_guideend_setpos] component Lmon_guideend=L_monitor() SETTING [/usr/share/mcstas/3.0-dev/monitors/L_monitor.comp:65]");
  stracpy(_Lmon_guideend->_name, "Lmon_guideend", 16384);
  stracpy(_Lmon_guideend->_type, "L_monitor", 16384);
  _Lmon_guideend->_index=40;
  _Lmon_guideend->_parameters.nL = 1000;
  #define nL (_Lmon_guideend->_parameters.nL)
  if("Lmon_guideend.dat" && strlen("Lmon_guideend.dat"))
    stracpy(_Lmon_guideend->_parameters.filename, "Lmon_guideend.dat" ? "Lmon_guideend.dat" : "", 16384);
  else 
  _Lmon_guideend->_parameters.filename[0]='\0';
  #define filename (_Lmon_guideend->_parameters.filename)
  _Lmon_guideend->_parameters.xmin = -0.05;
  #define xmin (_Lmon_guideend->_parameters.xmin)
  _Lmon_guideend->_parameters.xmax = 0.05;
  #define xmax (_Lmon_guideend->_parameters.xmax)
  _Lmon_guideend->_parameters.ymin = -0.05;
  #define ymin (_Lmon_guideend->_parameters.ymin)
  _Lmon_guideend->_parameters.ymax = 0.05;
  #define ymax (_Lmon_guideend->_parameters.ymax)
  _Lmon_guideend->_parameters.xwidth = instrument->_parameters._W_end + 0.01;
  #define xwidth (_Lmon_guideend->_parameters.xwidth)
  _Lmon_guideend->_parameters.yheight = instrument->_parameters._H_end + 0.01;
  #define yheight (_Lmon_guideend->_parameters.yheight)
  _Lmon_guideend->_parameters.Lmin = 0;
  #define Lmin (_Lmon_guideend->_parameters.Lmin)
  _Lmon_guideend->_parameters.Lmax = instrument->_parameters._Lmax + 1;
  #define Lmax (_Lmon_guideend->_parameters.Lmax)
  _Lmon_guideend->_parameters.restore_neutron = 0;
  #define restore_neutron (_Lmon_guideend->_parameters.restore_neutron)

  #define L_N (_Lmon_guideend->_parameters.L_N)
  #define L_p (_Lmon_guideend->_parameters.L_p)
  #define L_p2 (_Lmon_guideend->_parameters.L_p2)

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
  /* component Lmon_guideend=L_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Guidesample->_rotation_absolute, _Lmon_guideend->_rotation_absolute);
    rot_transpose(_Guidesample->_rotation_absolute, tr1);
    rot_mul(_Lmon_guideend->_rotation_absolute, tr1, _Lmon_guideend->_rotation_relative);
    _Lmon_guideend->_rotation_is_identity =  rot_test_identity(_Lmon_guideend->_rotation_relative);
    tc1 = coords_set(
      0, 0, instrument->_parameters._SAMPLE_DIST -4 * instrument->_parameters._GUI_GAP + 1e-6);
    rot_transpose(_Guidesample->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Lmon_guideend->_position_absolute = coords_add(_Guidesample->_position_absolute, tc2);
    tc1 = coords_sub(_Guidesample->_position_absolute, _Lmon_guideend->_position_absolute);
    _Lmon_guideend->_position_relative = rot_apply(_Lmon_guideend->_rotation_absolute, tc1);
  } /* Lmon_guideend=L_monitor() AT ROTATED */
  DEBUG_COMPONENT("Lmon_guideend", _Lmon_guideend->_position_absolute, _Lmon_guideend->_rotation_absolute);
  instrument->_position_absolute[40] = _Lmon_guideend->_position_absolute;
  instrument->_position_relative[40] = _Lmon_guideend->_position_relative;
  instrument->counter_N[40]  = instrument->counter_P[40] = instrument->counter_P2[40] = 0;
  instrument->counter_AbsorbProp[40]= 0;
  return(0);
} /* _Lmon_guideend_setpos */

/* component PSDsample=PSD_monitor() SETTING, POSITION/ROTATION */
int _PSDsample_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_PSDsample_setpos] component PSDsample=PSD_monitor() SETTING [/usr/share/mcstas/3.0-dev/monitors/PSD_monitor.comp:64]");
  stracpy(_PSDsample->_name, "PSDsample", 16384);
  stracpy(_PSDsample->_type, "PSD_monitor", 16384);
  _PSDsample->_index=41;
  _PSDsample->_parameters.nx = 90;
  #define nx (_PSDsample->_parameters.nx)
  _PSDsample->_parameters.ny = 90;
  #define ny (_PSDsample->_parameters.ny)
  if("PSDsample.dat" && strlen("PSDsample.dat"))
    stracpy(_PSDsample->_parameters.filename, "PSDsample.dat" ? "PSDsample.dat" : "", 16384);
  else 
  _PSDsample->_parameters.filename[0]='\0';
  #define filename (_PSDsample->_parameters.filename)
  _PSDsample->_parameters.xmin = -0.05;
  #define xmin (_PSDsample->_parameters.xmin)
  _PSDsample->_parameters.xmax = 0.05;
  #define xmax (_PSDsample->_parameters.xmax)
  _PSDsample->_parameters.ymin = -0.05;
  #define ymin (_PSDsample->_parameters.ymin)
  _PSDsample->_parameters.ymax = 0.05;
  #define ymax (_PSDsample->_parameters.ymax)
  _PSDsample->_parameters.xwidth = 0.1;
  #define xwidth (_PSDsample->_parameters.xwidth)
  _PSDsample->_parameters.yheight = 0.25;
  #define yheight (_PSDsample->_parameters.yheight)
  _PSDsample->_parameters.restore_neutron = 0;
  #define restore_neutron (_PSDsample->_parameters.restore_neutron)

  #define PSD_N (_PSDsample->_parameters.PSD_N)
  #define PSD_p (_PSDsample->_parameters.PSD_p)
  #define PSD_p2 (_PSDsample->_parameters.PSD_p2)

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
  /* component PSDsample=PSD_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Lmon_guideend->_rotation_absolute, _PSDsample->_rotation_absolute);
    rot_transpose(_Lmon_guideend->_rotation_absolute, tr1);
    rot_mul(_PSDsample->_rotation_absolute, tr1, _PSDsample->_rotation_relative);
    _PSDsample->_rotation_is_identity =  rot_test_identity(_PSDsample->_rotation_relative);
    tc1 = coords_set(
      0, 0, instrument->_parameters._SAMPLE_DIST -3 * instrument->_parameters._GUI_GAP + 1e-6);
    rot_transpose(_Lmon_guideend->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _PSDsample->_position_absolute = coords_add(_Lmon_guideend->_position_absolute, tc2);
    tc1 = coords_sub(_Lmon_guideend->_position_absolute, _PSDsample->_position_absolute);
    _PSDsample->_position_relative = rot_apply(_PSDsample->_rotation_absolute, tc1);
  } /* PSDsample=PSD_monitor() AT ROTATED */
  DEBUG_COMPONENT("PSDsample", _PSDsample->_position_absolute, _PSDsample->_rotation_absolute);
  instrument->_position_absolute[41] = _PSDsample->_position_absolute;
  instrument->_position_relative[41] = _PSDsample->_position_relative;
  instrument->counter_N[41]  = instrument->counter_P[41] = instrument->counter_P2[41] = 0;
  instrument->counter_AbsorbProp[41]= 0;
  return(0);
} /* _PSDsample_setpos */

/* component TOFsample_zoom=TOF_monitor() SETTING, POSITION/ROTATION */
int _TOFsample_zoom_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_TOFsample_zoom_setpos] component TOFsample_zoom=TOF_monitor() SETTING [/usr/share/mcstas/3.0-dev/monitors/TOF_monitor.comp:63]");
  stracpy(_TOFsample_zoom->_name, "TOFsample_zoom", 16384);
  stracpy(_TOFsample_zoom->_type, "TOF_monitor", 16384);
  _TOFsample_zoom->_index=42;
  _TOFsample_zoom->_parameters.nt = 500;
  #define nt (_TOFsample_zoom->_parameters.nt)
  if("TOF_sample.dat" && strlen("TOF_sample.dat"))
    stracpy(_TOFsample_zoom->_parameters.filename, "TOF_sample.dat" ? "TOF_sample.dat" : "", 16384);
  else 
  _TOFsample_zoom->_parameters.filename[0]='\0';
  #define filename (_TOFsample_zoom->_parameters.filename)
  _TOFsample_zoom->_parameters.xmin = -0.05;
  #define xmin (_TOFsample_zoom->_parameters.xmin)
  _TOFsample_zoom->_parameters.xmax = 0.05;
  #define xmax (_TOFsample_zoom->_parameters.xmax)
  _TOFsample_zoom->_parameters.ymin = -0.05;
  #define ymin (_TOFsample_zoom->_parameters.ymin)
  _TOFsample_zoom->_parameters.ymax = 0.05;
  #define ymax (_TOFsample_zoom->_parameters.ymax)
  _TOFsample_zoom->_parameters.xwidth = 0.02;
  #define xwidth (_TOFsample_zoom->_parameters.xwidth)
  _TOFsample_zoom->_parameters.yheight = 0.04;
  #define yheight (_TOFsample_zoom->_parameters.yheight)
  _TOFsample_zoom->_parameters.tmin = 1e6 * t_sample -5e4;
  #define tmin (_TOFsample_zoom->_parameters.tmin)
  _TOFsample_zoom->_parameters.tmax = 1e6 * t_sample + 5e4;
  #define tmax (_TOFsample_zoom->_parameters.tmax)
  _TOFsample_zoom->_parameters.dt = 1.0;
  #define dt (_TOFsample_zoom->_parameters.dt)
  _TOFsample_zoom->_parameters.restore_neutron = 0;
  #define restore_neutron (_TOFsample_zoom->_parameters.restore_neutron)

  #define TOF_N (_TOFsample_zoom->_parameters.TOF_N)
  #define TOF_p (_TOFsample_zoom->_parameters.TOF_p)
  #define TOF_p2 (_TOFsample_zoom->_parameters.TOF_p2)
  #define t_min (_TOFsample_zoom->_parameters.t_min)
  #define t_max (_TOFsample_zoom->_parameters.t_max)
  #define delta_t (_TOFsample_zoom->_parameters.delta_t)

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
  /* component TOFsample_zoom=TOF_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _PSDsample->_rotation_absolute, _TOFsample_zoom->_rotation_absolute);
    rot_transpose(_PSDsample->_rotation_absolute, tr1);
    rot_mul(_TOFsample_zoom->_rotation_absolute, tr1, _TOFsample_zoom->_rotation_relative);
    _TOFsample_zoom->_rotation_is_identity =  rot_test_identity(_TOFsample_zoom->_rotation_relative);
    tc1 = coords_set(
      0, 0, 1e-6);
    rot_transpose(_PSDsample->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _TOFsample_zoom->_position_absolute = coords_add(_PSDsample->_position_absolute, tc2);
    tc1 = coords_sub(_PSDsample->_position_absolute, _TOFsample_zoom->_position_absolute);
    _TOFsample_zoom->_position_relative = rot_apply(_TOFsample_zoom->_rotation_absolute, tc1);
  } /* TOFsample_zoom=TOF_monitor() AT ROTATED */
  DEBUG_COMPONENT("TOFsample_zoom", _TOFsample_zoom->_position_absolute, _TOFsample_zoom->_rotation_absolute);
  instrument->_position_absolute[42] = _TOFsample_zoom->_position_absolute;
  instrument->_position_relative[42] = _TOFsample_zoom->_position_relative;
  instrument->counter_N[42]  = instrument->counter_P[42] = instrument->counter_P2[42] = 0;
  instrument->counter_AbsorbProp[42]= 0;
  return(0);
} /* _TOFsample_zoom_setpos */

/* component Esample=E_monitor() SETTING, POSITION/ROTATION */
int _Esample_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Esample_setpos] component Esample=E_monitor() SETTING [/usr/share/mcstas/3.0-dev/monitors/E_monitor.comp:66]");
  stracpy(_Esample->_name, "Esample", 16384);
  stracpy(_Esample->_type, "E_monitor", 16384);
  _Esample->_index=43;
  _Esample->_parameters.nE = 400;
  #define nE (_Esample->_parameters.nE)
  if("Esample" && strlen("Esample"))
    stracpy(_Esample->_parameters.filename, "Esample" ? "Esample" : "", 16384);
  else 
  _Esample->_parameters.filename[0]='\0';
  #define filename (_Esample->_parameters.filename)
  _Esample->_parameters.xmin = -0.05;
  #define xmin (_Esample->_parameters.xmin)
  _Esample->_parameters.xmax = 0.05;
  #define xmax (_Esample->_parameters.xmax)
  _Esample->_parameters.ymin = -0.05;
  #define ymin (_Esample->_parameters.ymin)
  _Esample->_parameters.ymax = 0.05;
  #define ymax (_Esample->_parameters.ymax)
  _Esample->_parameters.xwidth = 0.02;
  #define xwidth (_Esample->_parameters.xwidth)
  _Esample->_parameters.yheight = 0.04;
  #define yheight (_Esample->_parameters.yheight)
  _Esample->_parameters.Emin = E_target - instrument->_parameters._RES_DE;
  #define Emin (_Esample->_parameters.Emin)
  _Esample->_parameters.Emax = E_target + instrument->_parameters._RES_DE;
  #define Emax (_Esample->_parameters.Emax)
  _Esample->_parameters.restore_neutron = 0;
  #define restore_neutron (_Esample->_parameters.restore_neutron)

  #define E_N (_Esample->_parameters.E_N)
  #define E_p (_Esample->_parameters.E_p)
  #define E_p2 (_Esample->_parameters.E_p2)
  #define S_p (_Esample->_parameters.S_p)
  #define S_pE (_Esample->_parameters.S_pE)
  #define S_pE2 (_Esample->_parameters.S_pE2)

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
  /* component Esample=E_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _TOFsample_zoom->_rotation_absolute, _Esample->_rotation_absolute);
    rot_transpose(_TOFsample_zoom->_rotation_absolute, tr1);
    rot_mul(_Esample->_rotation_absolute, tr1, _Esample->_rotation_relative);
    _Esample->_rotation_is_identity =  rot_test_identity(_Esample->_rotation_relative);
    tc1 = coords_set(
      0, 0, 1e-6);
    rot_transpose(_TOFsample_zoom->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Esample->_position_absolute = coords_add(_TOFsample_zoom->_position_absolute, tc2);
    tc1 = coords_sub(_TOFsample_zoom->_position_absolute, _Esample->_position_absolute);
    _Esample->_position_relative = rot_apply(_Esample->_rotation_absolute, tc1);
  } /* Esample=E_monitor() AT ROTATED */
  DEBUG_COMPONENT("Esample", _Esample->_position_absolute, _Esample->_rotation_absolute);
  instrument->_position_absolute[43] = _Esample->_position_absolute;
  instrument->_position_relative[43] = _Esample->_position_relative;
  instrument->counter_N[43]  = instrument->counter_P[43] = instrument->counter_P2[43] = 0;
  instrument->counter_AbsorbProp[43]= 0;
  return(0);
} /* _Esample_setpos */

/* component Lmon_sample_zoom=L_monitor() SETTING, POSITION/ROTATION */
int _Lmon_sample_zoom_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Lmon_sample_zoom_setpos] component Lmon_sample_zoom=L_monitor() SETTING [/usr/share/mcstas/3.0-dev/monitors/L_monitor.comp:65]");
  stracpy(_Lmon_sample_zoom->_name, "Lmon_sample_zoom", 16384);
  stracpy(_Lmon_sample_zoom->_type, "L_monitor", 16384);
  _Lmon_sample_zoom->_index=44;
  _Lmon_sample_zoom->_parameters.nL = 100;
  #define nL (_Lmon_sample_zoom->_parameters.nL)
  if("LMON_sample_zoom.dat" && strlen("LMON_sample_zoom.dat"))
    stracpy(_Lmon_sample_zoom->_parameters.filename, "LMON_sample_zoom.dat" ? "LMON_sample_zoom.dat" : "", 16384);
  else 
  _Lmon_sample_zoom->_parameters.filename[0]='\0';
  #define filename (_Lmon_sample_zoom->_parameters.filename)
  _Lmon_sample_zoom->_parameters.xmin = -0.05;
  #define xmin (_Lmon_sample_zoom->_parameters.xmin)
  _Lmon_sample_zoom->_parameters.xmax = 0.05;
  #define xmax (_Lmon_sample_zoom->_parameters.xmax)
  _Lmon_sample_zoom->_parameters.ymin = -0.05;
  #define ymin (_Lmon_sample_zoom->_parameters.ymin)
  _Lmon_sample_zoom->_parameters.ymax = 0.05;
  #define ymax (_Lmon_sample_zoom->_parameters.ymax)
  _Lmon_sample_zoom->_parameters.xwidth = 0.02;
  #define xwidth (_Lmon_sample_zoom->_parameters.xwidth)
  _Lmon_sample_zoom->_parameters.yheight = 0.04;
  #define yheight (_Lmon_sample_zoom->_parameters.yheight)
  _Lmon_sample_zoom->_parameters.Lmin = instrument->_parameters._lambda0 -0.1;
  #define Lmin (_Lmon_sample_zoom->_parameters.Lmin)
  _Lmon_sample_zoom->_parameters.Lmax = instrument->_parameters._lambda0 + 0.1;
  #define Lmax (_Lmon_sample_zoom->_parameters.Lmax)
  _Lmon_sample_zoom->_parameters.restore_neutron = 0;
  #define restore_neutron (_Lmon_sample_zoom->_parameters.restore_neutron)

  #define L_N (_Lmon_sample_zoom->_parameters.L_N)
  #define L_p (_Lmon_sample_zoom->_parameters.L_p)
  #define L_p2 (_Lmon_sample_zoom->_parameters.L_p2)

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
  /* component Lmon_sample_zoom=L_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Esample->_rotation_absolute, _Lmon_sample_zoom->_rotation_absolute);
    rot_transpose(_Esample->_rotation_absolute, tr1);
    rot_mul(_Lmon_sample_zoom->_rotation_absolute, tr1, _Lmon_sample_zoom->_rotation_relative);
    _Lmon_sample_zoom->_rotation_is_identity =  rot_test_identity(_Lmon_sample_zoom->_rotation_relative);
    tc1 = coords_set(
      0, 0, 1e-6);
    rot_transpose(_Esample->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Lmon_sample_zoom->_position_absolute = coords_add(_Esample->_position_absolute, tc2);
    tc1 = coords_sub(_Esample->_position_absolute, _Lmon_sample_zoom->_position_absolute);
    _Lmon_sample_zoom->_position_relative = rot_apply(_Lmon_sample_zoom->_rotation_absolute, tc1);
  } /* Lmon_sample_zoom=L_monitor() AT ROTATED */
  DEBUG_COMPONENT("Lmon_sample_zoom", _Lmon_sample_zoom->_position_absolute, _Lmon_sample_zoom->_rotation_absolute);
  instrument->_position_absolute[44] = _Lmon_sample_zoom->_position_absolute;
  instrument->_position_relative[44] = _Lmon_sample_zoom->_position_relative;
  instrument->counter_N[44]  = instrument->counter_P[44] = instrument->counter_P2[44] = 0;
  instrument->counter_AbsorbProp[44]= 0;
  return(0);
} /* _Lmon_sample_zoom_setpos */

/* component sample=Tunneling_sample() SETTING, POSITION/ROTATION */
int _sample_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_sample_setpos] component sample=Tunneling_sample() SETTING [/usr/share/mcstas/3.0-dev/samples/Tunneling_sample.comp:114]");
  stracpy(_sample->_name, "sample", 16384);
  stracpy(_sample->_type, "Tunneling_sample", 16384);
  _sample->_index=45;
  _sample->_parameters.thickness = 0.01 - instrument->_parameters._V_HOLE;
  #define thickness (_sample->_parameters.thickness)
  _sample->_parameters.radius = 0.01;
  #define radius (_sample->_parameters.radius)
  _sample->_parameters.focus_r = 0;
  #define focus_r (_sample->_parameters.focus_r)
  _sample->_parameters.p_interact = 1;
  #define p_interact (_sample->_parameters.p_interact)
  _sample->_parameters.f_QE = instrument->_parameters._FRAC_QUASIEL;
  #define f_QE (_sample->_parameters.f_QE)
  _sample->_parameters.f_tun = instrument->_parameters._FRAC_TUNNEL;
  #define f_tun (_sample->_parameters.f_tun)
  _sample->_parameters.gamma = instrument->_parameters._Gamma;
  #define gamma (_sample->_parameters.gamma)
  _sample->_parameters.E_tun = instrument->_parameters._Etun;
  #define E_tun (_sample->_parameters.E_tun)
  _sample->_parameters.target_x = 0;
  #define target_x (_sample->_parameters.target_x)
  _sample->_parameters.target_y = 0;
  #define target_y (_sample->_parameters.target_y)
  _sample->_parameters.target_z = 0.235;
  #define target_z (_sample->_parameters.target_z)
  _sample->_parameters.focus_xw = 0.015;
  #define focus_xw (_sample->_parameters.focus_xw)
  _sample->_parameters.focus_yh = 0.015;
  #define focus_yh (_sample->_parameters.focus_yh)
  _sample->_parameters.focus_aw = 0;
  #define focus_aw (_sample->_parameters.focus_aw)
  _sample->_parameters.focus_ah = 0;
  #define focus_ah (_sample->_parameters.focus_ah)
  _sample->_parameters.xwidth = 0;
  #define xwidth (_sample->_parameters.xwidth)
  _sample->_parameters.yheight = 0.04;
  #define yheight (_sample->_parameters.yheight)
  _sample->_parameters.zdepth = 0;
  #define zdepth (_sample->_parameters.zdepth)
  _sample->_parameters.sigma_abs = 5.08;
  #define sigma_abs (_sample->_parameters.sigma_abs)
  _sample->_parameters.sigma_inc = 4.935;
  #define sigma_inc (_sample->_parameters.sigma_inc)
  _sample->_parameters.Vc = 13.827;
  #define Vc (_sample->_parameters.Vc)
  _sample->_parameters.target_index = 2;
  #define target_index (_sample->_parameters.target_index)

  #define ftun (_sample->_parameters.ftun)
  #define fQE (_sample->_parameters.fQE)
  #define VarsV (_sample->_parameters.VarsV)

  #undef thickness
  #undef radius
  #undef focus_r
  #undef p_interact
  #undef f_QE
  #undef f_tun
  #undef gamma
  #undef E_tun
  #undef target_x
  #undef target_y
  #undef target_z
  #undef focus_xw
  #undef focus_yh
  #undef focus_aw
  #undef focus_ah
  #undef xwidth
  #undef yheight
  #undef zdepth
  #undef sigma_abs
  #undef sigma_inc
  #undef Vc
  #undef target_index
  #undef ftun
  #undef fQE
  #undef VarsV
  /* component sample=Tunneling_sample() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Fastchop2->_rotation_absolute, _sample->_rotation_absolute);
    rot_transpose(_Lmon_sample_zoom->_rotation_absolute, tr1);
    rot_mul(_sample->_rotation_absolute, tr1, _sample->_rotation_relative);
    _sample->_rotation_is_identity =  rot_test_identity(_sample->_rotation_relative);
    tc1 = coords_set(
      0, -0.04, instrument->_parameters._SAMPLE_DIST);
    rot_transpose(_Fastchop2->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _sample->_position_absolute = coords_add(_Fastchop2->_position_absolute, tc2);
    tc1 = coords_sub(_Lmon_sample_zoom->_position_absolute, _sample->_position_absolute);
    _sample->_position_relative = rot_apply(_sample->_rotation_absolute, tc1);
  } /* sample=Tunneling_sample() AT ROTATED */
  DEBUG_COMPONENT("sample", _sample->_position_absolute, _sample->_rotation_absolute);
  instrument->_position_absolute[45] = _sample->_position_absolute;
  instrument->_position_relative[45] = _sample->_position_relative;
  instrument->counter_N[45]  = instrument->counter_P[45] = instrument->counter_P2[45] = 0;
  instrument->counter_AbsorbProp[45]= 0;
  return(0);
} /* _sample_setpos */

/* component detectorarm=Arm() SETTING, POSITION/ROTATION */
int _detectorarm_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_detectorarm_setpos] component detectorarm=Arm() SETTING [Arm:0]");
  stracpy(_detectorarm->_name, "detectorarm", 16384);
  stracpy(_detectorarm->_type, "Arm", 16384);
  _detectorarm->_index=46;
  /* component detectorarm=Arm() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (instrument->_parameters._TT)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _sample->_rotation_absolute, _detectorarm->_rotation_absolute);
    rot_transpose(_sample->_rotation_absolute, tr1);
    rot_mul(_detectorarm->_rotation_absolute, tr1, _detectorarm->_rotation_relative);
    _detectorarm->_rotation_is_identity =  rot_test_identity(_detectorarm->_rotation_relative);
    tc1 = coords_set(
      0, 0, 0);
    rot_transpose(_sample->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _detectorarm->_position_absolute = coords_add(_sample->_position_absolute, tc2);
    tc1 = coords_sub(_sample->_position_absolute, _detectorarm->_position_absolute);
    _detectorarm->_position_relative = rot_apply(_detectorarm->_rotation_absolute, tc1);
  } /* detectorarm=Arm() AT ROTATED */
  DEBUG_COMPONENT("detectorarm", _detectorarm->_position_absolute, _detectorarm->_rotation_absolute);
  instrument->_position_absolute[46] = _detectorarm->_position_absolute;
  instrument->_position_relative[46] = _detectorarm->_position_relative;
  instrument->counter_N[46]  = instrument->counter_P[46] = instrument->counter_P2[46] = 0;
  instrument->counter_AbsorbProp[46]= 0;
  return(0);
} /* _detectorarm_setpos */

/* component TOFdetector=TOF_monitor() SETTING, POSITION/ROTATION */
int _TOFdetector_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_TOFdetector_setpos] component TOFdetector=TOF_monitor() SETTING [/usr/share/mcstas/3.0-dev/monitors/TOF_monitor.comp:63]");
  stracpy(_TOFdetector->_name, "TOFdetector", 16384);
  stracpy(_TOFdetector->_type, "TOF_monitor", 16384);
  _TOFdetector->_index=47;
  _TOFdetector->_parameters.nt = 512;
  #define nt (_TOFdetector->_parameters.nt)
  if("TOF.dat" && strlen("TOF.dat"))
    stracpy(_TOFdetector->_parameters.filename, "TOF.dat" ? "TOF.dat" : "", 16384);
  else 
  _TOFdetector->_parameters.filename[0]='\0';
  #define filename (_TOFdetector->_parameters.filename)
  _TOFdetector->_parameters.xmin = -0.05;
  #define xmin (_TOFdetector->_parameters.xmin)
  _TOFdetector->_parameters.xmax = 0.05;
  #define xmax (_TOFdetector->_parameters.xmax)
  _TOFdetector->_parameters.ymin = -0.05;
  #define ymin (_TOFdetector->_parameters.ymin)
  _TOFdetector->_parameters.ymax = 0.05;
  #define ymax (_TOFdetector->_parameters.ymax)
  _TOFdetector->_parameters.xwidth = 0.015;
  #define xwidth (_TOFdetector->_parameters.xwidth)
  _TOFdetector->_parameters.yheight = 0.015;
  #define yheight (_TOFdetector->_parameters.yheight)
  _TOFdetector->_parameters.tmin = 0;
  #define tmin (_TOFdetector->_parameters.tmin)
  _TOFdetector->_parameters.tmax = 2e5;
  #define tmax (_TOFdetector->_parameters.tmax)
  _TOFdetector->_parameters.dt = 1.0;
  #define dt (_TOFdetector->_parameters.dt)
  _TOFdetector->_parameters.restore_neutron = 0;
  #define restore_neutron (_TOFdetector->_parameters.restore_neutron)

  #define TOF_N (_TOFdetector->_parameters.TOF_N)
  #define TOF_p (_TOFdetector->_parameters.TOF_p)
  #define TOF_p2 (_TOFdetector->_parameters.TOF_p2)
  #define t_min (_TOFdetector->_parameters.t_min)
  #define t_max (_TOFdetector->_parameters.t_max)
  #define delta_t (_TOFdetector->_parameters.delta_t)

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
  /* component TOFdetector=TOF_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _detectorarm->_rotation_absolute, _TOFdetector->_rotation_absolute);
    rot_transpose(_detectorarm->_rotation_absolute, tr1);
    rot_mul(_TOFdetector->_rotation_absolute, tr1, _TOFdetector->_rotation_relative);
    _TOFdetector->_rotation_is_identity =  rot_test_identity(_TOFdetector->_rotation_relative);
    tc1 = coords_set(
      0, 0, instrument->_parameters._DETECTOR_DIST);
    rot_transpose(_detectorarm->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _TOFdetector->_position_absolute = coords_add(_detectorarm->_position_absolute, tc2);
    tc1 = coords_sub(_detectorarm->_position_absolute, _TOFdetector->_position_absolute);
    _TOFdetector->_position_relative = rot_apply(_TOFdetector->_rotation_absolute, tc1);
  } /* TOFdetector=TOF_monitor() AT ROTATED */
  DEBUG_COMPONENT("TOFdetector", _TOFdetector->_position_absolute, _TOFdetector->_rotation_absolute);
  instrument->_position_absolute[47] = _TOFdetector->_position_absolute;
  instrument->_position_relative[47] = _TOFdetector->_position_relative;
  instrument->counter_N[47]  = instrument->counter_P[47] = instrument->counter_P2[47] = 0;
  instrument->counter_AbsorbProp[47]= 0;
  return(0);
} /* _TOFdetector_setpos */

/* component TOFdetector_zoom=TOF_monitor() SETTING, POSITION/ROTATION */
int _TOFdetector_zoom_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_TOFdetector_zoom_setpos] component TOFdetector_zoom=TOF_monitor() SETTING [/usr/share/mcstas/3.0-dev/monitors/TOF_monitor.comp:63]");
  stracpy(_TOFdetector_zoom->_name, "TOFdetector_zoom", 16384);
  stracpy(_TOFdetector_zoom->_type, "TOF_monitor", 16384);
  _TOFdetector_zoom->_index=48;
  _TOFdetector_zoom->_parameters.nt = 100;
  #define nt (_TOFdetector_zoom->_parameters.nt)
  if("TOF_zoom.dat" && strlen("TOF_zoom.dat"))
    stracpy(_TOFdetector_zoom->_parameters.filename, "TOF_zoom.dat" ? "TOF_zoom.dat" : "", 16384);
  else 
  _TOFdetector_zoom->_parameters.filename[0]='\0';
  #define filename (_TOFdetector_zoom->_parameters.filename)
  _TOFdetector_zoom->_parameters.xmin = -0.05;
  #define xmin (_TOFdetector_zoom->_parameters.xmin)
  _TOFdetector_zoom->_parameters.xmax = 0.05;
  #define xmax (_TOFdetector_zoom->_parameters.xmax)
  _TOFdetector_zoom->_parameters.ymin = -0.05;
  #define ymin (_TOFdetector_zoom->_parameters.ymin)
  _TOFdetector_zoom->_parameters.ymax = 0.05;
  #define ymax (_TOFdetector_zoom->_parameters.ymax)
  _TOFdetector_zoom->_parameters.xwidth = 0.015;
  #define xwidth (_TOFdetector_zoom->_parameters.xwidth)
  _TOFdetector_zoom->_parameters.yheight = 0.015;
  #define yheight (_TOFdetector_zoom->_parameters.yheight)
  _TOFdetector_zoom->_parameters.tmin = 1e6 * t_detector -10e2;
  #define tmin (_TOFdetector_zoom->_parameters.tmin)
  _TOFdetector_zoom->_parameters.tmax = 1e6 * t_detector + 10e2;
  #define tmax (_TOFdetector_zoom->_parameters.tmax)
  _TOFdetector_zoom->_parameters.dt = 1.0;
  #define dt (_TOFdetector_zoom->_parameters.dt)
  _TOFdetector_zoom->_parameters.restore_neutron = 0;
  #define restore_neutron (_TOFdetector_zoom->_parameters.restore_neutron)

  #define TOF_N (_TOFdetector_zoom->_parameters.TOF_N)
  #define TOF_p (_TOFdetector_zoom->_parameters.TOF_p)
  #define TOF_p2 (_TOFdetector_zoom->_parameters.TOF_p2)
  #define t_min (_TOFdetector_zoom->_parameters.t_min)
  #define t_max (_TOFdetector_zoom->_parameters.t_max)
  #define delta_t (_TOFdetector_zoom->_parameters.delta_t)

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
  /* component TOFdetector_zoom=TOF_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _TOFdetector->_rotation_absolute, _TOFdetector_zoom->_rotation_absolute);
    rot_transpose(_TOFdetector->_rotation_absolute, tr1);
    rot_mul(_TOFdetector_zoom->_rotation_absolute, tr1, _TOFdetector_zoom->_rotation_relative);
    _TOFdetector_zoom->_rotation_is_identity =  rot_test_identity(_TOFdetector_zoom->_rotation_relative);
    tc1 = coords_set(
      0, 0, 1e-6);
    rot_transpose(_TOFdetector->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _TOFdetector_zoom->_position_absolute = coords_add(_TOFdetector->_position_absolute, tc2);
    tc1 = coords_sub(_TOFdetector->_position_absolute, _TOFdetector_zoom->_position_absolute);
    _TOFdetector_zoom->_position_relative = rot_apply(_TOFdetector_zoom->_rotation_absolute, tc1);
  } /* TOFdetector_zoom=TOF_monitor() AT ROTATED */
  DEBUG_COMPONENT("TOFdetector_zoom", _TOFdetector_zoom->_position_absolute, _TOFdetector_zoom->_rotation_absolute);
  instrument->_position_absolute[48] = _TOFdetector_zoom->_position_absolute;
  instrument->_position_relative[48] = _TOFdetector_zoom->_position_relative;
  instrument->counter_N[48]  = instrument->counter_P[48] = instrument->counter_P2[48] = 0;
  instrument->counter_AbsorbProp[48]= 0;
  return(0);
} /* _TOFdetector_zoom_setpos */

/* component Edetector=E_monitor() SETTING, POSITION/ROTATION */
int _Edetector_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Edetector_setpos] component Edetector=E_monitor() SETTING [/usr/share/mcstas/3.0-dev/monitors/E_monitor.comp:66]");
  stracpy(_Edetector->_name, "Edetector", 16384);
  stracpy(_Edetector->_type, "E_monitor", 16384);
  _Edetector->_index=49;
  _Edetector->_parameters.nE = 400;
  #define nE (_Edetector->_parameters.nE)
  if("Edet" && strlen("Edet"))
    stracpy(_Edetector->_parameters.filename, "Edet" ? "Edet" : "", 16384);
  else 
  _Edetector->_parameters.filename[0]='\0';
  #define filename (_Edetector->_parameters.filename)
  _Edetector->_parameters.xmin = -0.05;
  #define xmin (_Edetector->_parameters.xmin)
  _Edetector->_parameters.xmax = 0.05;
  #define xmax (_Edetector->_parameters.xmax)
  _Edetector->_parameters.ymin = -0.05;
  #define ymin (_Edetector->_parameters.ymin)
  _Edetector->_parameters.ymax = 0.05;
  #define ymax (_Edetector->_parameters.ymax)
  _Edetector->_parameters.xwidth = 0.015;
  #define xwidth (_Edetector->_parameters.xwidth)
  _Edetector->_parameters.yheight = 0.015;
  #define yheight (_Edetector->_parameters.yheight)
  _Edetector->_parameters.Emin = E_target - instrument->_parameters._RES_DE;
  #define Emin (_Edetector->_parameters.Emin)
  _Edetector->_parameters.Emax = E_target + instrument->_parameters._RES_DE;
  #define Emax (_Edetector->_parameters.Emax)
  _Edetector->_parameters.restore_neutron = 0;
  #define restore_neutron (_Edetector->_parameters.restore_neutron)

  #define E_N (_Edetector->_parameters.E_N)
  #define E_p (_Edetector->_parameters.E_p)
  #define E_p2 (_Edetector->_parameters.E_p2)
  #define S_p (_Edetector->_parameters.S_p)
  #define S_pE (_Edetector->_parameters.S_pE)
  #define S_pE2 (_Edetector->_parameters.S_pE2)

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
  /* component Edetector=E_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _TOFdetector_zoom->_rotation_absolute, _Edetector->_rotation_absolute);
    rot_transpose(_TOFdetector_zoom->_rotation_absolute, tr1);
    rot_mul(_Edetector->_rotation_absolute, tr1, _Edetector->_rotation_relative);
    _Edetector->_rotation_is_identity =  rot_test_identity(_Edetector->_rotation_relative);
    tc1 = coords_set(
      0, 0, 1e-6);
    rot_transpose(_TOFdetector_zoom->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Edetector->_position_absolute = coords_add(_TOFdetector_zoom->_position_absolute, tc2);
    tc1 = coords_sub(_TOFdetector_zoom->_position_absolute, _Edetector->_position_absolute);
    _Edetector->_position_relative = rot_apply(_Edetector->_rotation_absolute, tc1);
  } /* Edetector=E_monitor() AT ROTATED */
  DEBUG_COMPONENT("Edetector", _Edetector->_position_absolute, _Edetector->_rotation_absolute);
  instrument->_position_absolute[49] = _Edetector->_position_absolute;
  instrument->_position_relative[49] = _Edetector->_position_relative;
  instrument->counter_N[49]  = instrument->counter_P[49] = instrument->counter_P2[49] = 0;
  instrument->counter_AbsorbProp[49]= 0;
  return(0);
} /* _Edetector_setpos */

/* component TOF2Edetector=TOF2E_monitor() SETTING, POSITION/ROTATION */
int _TOF2Edetector_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_TOF2Edetector_setpos] component TOF2Edetector=TOF2E_monitor() SETTING [/usr/share/mcstas/3.0-dev/monitors/TOF2E_monitor.comp:67]");
  stracpy(_TOF2Edetector->_name, "TOF2Edetector", 16384);
  stracpy(_TOF2Edetector->_type, "TOF2E_monitor", 16384);
  _TOF2Edetector->_index=50;
  _TOF2Edetector->_parameters.nE = 200;
  #define nE (_TOF2Edetector->_parameters.nE)
  if("TOF2E.dat" && strlen("TOF2E.dat"))
    stracpy(_TOF2Edetector->_parameters.filename, "TOF2E.dat" ? "TOF2E.dat" : "", 16384);
  else 
  _TOF2Edetector->_parameters.filename[0]='\0';
  #define filename (_TOF2Edetector->_parameters.filename)
  _TOF2Edetector->_parameters.xmin = -0.05;
  #define xmin (_TOF2Edetector->_parameters.xmin)
  _TOF2Edetector->_parameters.xmax = 0.05;
  #define xmax (_TOF2Edetector->_parameters.xmax)
  _TOF2Edetector->_parameters.ymin = -0.05;
  #define ymin (_TOF2Edetector->_parameters.ymin)
  _TOF2Edetector->_parameters.ymax = 0.05;
  #define ymax (_TOF2Edetector->_parameters.ymax)
  _TOF2Edetector->_parameters.xwidth = 0.015;
  #define xwidth (_TOF2Edetector->_parameters.xwidth)
  _TOF2Edetector->_parameters.yheight = 0.015;
  #define yheight (_TOF2Edetector->_parameters.yheight)
  _TOF2Edetector->_parameters.Emin = E_target - instrument->_parameters._RES_DE;
  #define Emin (_TOF2Edetector->_parameters.Emin)
  _TOF2Edetector->_parameters.Emax = E_target + instrument->_parameters._RES_DE;
  #define Emax (_TOF2Edetector->_parameters.Emax)
  _TOF2Edetector->_parameters.T_zero = t_sample;
  #define T_zero (_TOF2Edetector->_parameters.T_zero)
  _TOF2Edetector->_parameters.L_flight = instrument->_parameters._DETECTOR_DIST;
  #define L_flight (_TOF2Edetector->_parameters.L_flight)
  _TOF2Edetector->_parameters.restore_neutron = 0;
  #define restore_neutron (_TOF2Edetector->_parameters.restore_neutron)

  #define E_N (_TOF2Edetector->_parameters.E_N)
  #define E_p (_TOF2Edetector->_parameters.E_p)
  #define E_p2 (_TOF2Edetector->_parameters.E_p2)
  #define S_p (_TOF2Edetector->_parameters.S_p)
  #define S_pE (_TOF2Edetector->_parameters.S_pE)
  #define S_pE2 (_TOF2Edetector->_parameters.S_pE2)

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
  #undef T_zero
  #undef L_flight
  #undef restore_neutron
  #undef E_N
  #undef E_p
  #undef E_p2
  #undef S_p
  #undef S_pE
  #undef S_pE2
  /* component TOF2Edetector=TOF2E_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Edetector->_rotation_absolute, _TOF2Edetector->_rotation_absolute);
    rot_transpose(_Edetector->_rotation_absolute, tr1);
    rot_mul(_TOF2Edetector->_rotation_absolute, tr1, _TOF2Edetector->_rotation_relative);
    _TOF2Edetector->_rotation_is_identity =  rot_test_identity(_TOF2Edetector->_rotation_relative);
    tc1 = coords_set(
      0, 0, 1e-6);
    rot_transpose(_Edetector->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _TOF2Edetector->_position_absolute = coords_add(_Edetector->_position_absolute, tc2);
    tc1 = coords_sub(_Edetector->_position_absolute, _TOF2Edetector->_position_absolute);
    _TOF2Edetector->_position_relative = rot_apply(_TOF2Edetector->_rotation_absolute, tc1);
  } /* TOF2Edetector=TOF2E_monitor() AT ROTATED */
  DEBUG_COMPONENT("TOF2Edetector", _TOF2Edetector->_position_absolute, _TOF2Edetector->_rotation_absolute);
  instrument->_position_absolute[50] = _TOF2Edetector->_position_absolute;
  instrument->_position_relative[50] = _TOF2Edetector->_position_relative;
  instrument->counter_N[50]  = instrument->counter_P[50] = instrument->counter_P2[50] = 0;
  instrument->counter_AbsorbProp[50]= 0;
  return(0);
} /* _TOF2Edetector_setpos */

_class_ESS_moderator_long *class_ESS_moderator_long_init(_class_ESS_moderator_long *_comp
) {
  #define width_c (_comp->_parameters.width_c)
  #define yheight (_comp->_parameters.yheight)
  #define Lmin (_comp->_parameters.Lmin)
  #define Lmax (_comp->_parameters.Lmax)
  #define dist (_comp->_parameters.dist)
  #define focus_xw (_comp->_parameters.focus_xw)
  #define focus_yh (_comp->_parameters.focus_yh)
  #define nu (_comp->_parameters.nu)
  #define T (_comp->_parameters.T)
  #define tau (_comp->_parameters.tau)
  #define tau1 (_comp->_parameters.tau1)
  #define tau2 (_comp->_parameters.tau2)
  #define d (_comp->_parameters.d)
  #define n (_comp->_parameters.n)
  #define cold_frac (_comp->_parameters.cold_frac)
  #define n2 (_comp->_parameters.n2)
  #define chi2 (_comp->_parameters.chi2)
  #define I0 (_comp->_parameters.I0)
  #define I2 (_comp->_parameters.I2)
  #define target_index (_comp->_parameters.target_index)
  #define cyl_radius (_comp->_parameters.cyl_radius)
  #define branch1 (_comp->_parameters.branch1)
  #define branch2 (_comp->_parameters.branch2)
  #define branch_tail (_comp->_parameters.branch_tail)
  #define n_pulses (_comp->_parameters.n_pulses)
  #define width_t (_comp->_parameters.width_t)
  #define T_t (_comp->_parameters.T_t)
  #define tau_t (_comp->_parameters.tau_t)
  #define tau1_t (_comp->_parameters.tau1_t)
  #define tau2_t (_comp->_parameters.tau2_t)
  #define chi2_t (_comp->_parameters.chi2_t)
  #define I0_t (_comp->_parameters.I0_t)
  #define I2_t (_comp->_parameters.I2_t)
  #define branch1_t (_comp->_parameters.branch1_t)
  #define branch2_t (_comp->_parameters.branch2_t)
  #define src_2012 (_comp->_parameters.src_2012)
  #define tfocus_dist (_comp->_parameters.tfocus_dist)
  #define tfocus_time (_comp->_parameters.tfocus_time)
  #define tfocus_width (_comp->_parameters.tfocus_width)
  #define beamport_angle (_comp->_parameters.beamport_angle)
  #define l_range (_comp->_parameters.l_range)
  #define w_mult (_comp->_parameters.w_mult)
  #define w_geom (_comp->_parameters.w_geom)
  #define w_geom_c (_comp->_parameters.w_geom_c)
  #define w_geom_t (_comp->_parameters.w_geom_t)
  #define tx (_comp->_parameters.tx)
  #define ty (_comp->_parameters.ty)
  #define tz (_comp->_parameters.tz)
  #define t1x (_comp->_parameters.t1x)
  #define t1y (_comp->_parameters.t1y)
  #define t1z (_comp->_parameters.t1z)
  #define t2x (_comp->_parameters.t2x)
  #define t2y (_comp->_parameters.t2y)
  #define t2z (_comp->_parameters.t2z)
  #define T_n (_comp->_parameters.T_n)
  #define tau_n (_comp->_parameters.tau_n)
  #define tau1_n (_comp->_parameters.tau1_n)
  #define tau2_n (_comp->_parameters.tau2_n)
  #define chi2_n (_comp->_parameters.chi2_n)
  #define I0_n (_comp->_parameters.I0_n)
  #define I2_n (_comp->_parameters.I2_n)
  #define branch1_n (_comp->_parameters.branch1_n)
  #define branch2_n (_comp->_parameters.branch2_n)
  #define r_empty (_comp->_parameters.r_empty)
  #define r_optics (_comp->_parameters.r_optics)
  if (Lmin>=Lmax) {
    printf("ESS_moderator_long: %s: Unmeaningful definition of wavelength range!\n ERROR - Exiting\n",
           NAME_CURRENT_COMP);
      exit(0);
  }

  n_pulses=(double)floor(n_pulses);
  r_empty = 2.0;
  if (n_pulses == 0) n_pulses=1;

  if (target_index && !dist)
  {
    Coords ToTarget;
    ToTarget = coords_sub(POS_A_COMP_INDEX(INDEX_CURRENT_COMP+target_index),POS_A_CURRENT_COMP);
    ToTarget = rot_apply(ROT_A_CURRENT_COMP, ToTarget);
    coords_get(ToTarget, &tx, &ty, &tz);
    dist=sqrt(tx*tx+ty*ty+tz*tz);
  } else if (target_index && !dist) {
    printf("ESS_moderator_long: %s: Please choose to set either the dist parameter or specify a target_index.\nExit\n", NAME_CURRENT_COMP);
    exit(-1);
  } else {
    tx=0, ty=0, tz=dist;
  }

  if (focus_xw < 0 || focus_yh < 0)
  {
    printf("ESS_moderator_long: %s: Please specify both focus_xw and focus_yh as positive numbers.\nExit\n", NAME_CURRENT_COMP);
    exit(-1);
  }

  if (dist < r_empty && dist > 0)
  {
    printf("ESS_moderator_long: %s WARNING: Provided dist parameter is %g and hence inside the vacated zone of the beam extraction system!\nYou might be placing optics in a restricted area!!!\n", NAME_CURRENT_COMP, dist);
  }

  if (beamport_angle < 0 || beamport_angle > 60)
  {
    printf("ESS_moderator_long: %s: Please select a beamport_angle between 0 and 60 degrees!\nExit\n", NAME_CURRENT_COMP);
    exit(-1);
  }

  if (width_c && cyl_radius) {
    printf("ESS_moderator_long: %s: Please specify EITHER cold-moderator radius (cyl_radius) or length of visible arch (width_c)!\nExit\n", NAME_CURRENT_COMP);
    exit(-1);
  } else if (cyl_radius) {
    width_c = 2*PI*cyl_radius*60/360;
  } else {
    cyl_radius = 360*width_c/(2*PI*60);
  }
  r_optics = 6.0 - r_empty - cyl_radius;

  if (n == 1 || n2 == 1 || Lmin<=0 || Lmax <=0 || dist == 0
    || branch2 == 0 || branch_tail == 0 || tau == 0)
  {
    printf("ESS_moderator_long: %s: Check parameters (lead to Math Error).\n Avoid 0 value for {Lmin Lmax dist d tau branch1/2/tail} and 1 value for {n n2 branch1/2/tail}\n", NAME_CURRENT_COMP);
    exit(-1);
  }

  if (tau1==0 && !(branch1==1)) {
    branch1=1;
    printf("ESS_moderator_long: %s: WARNING: Setting tau1 to zero implies branch 1=1.\n", NAME_CURRENT_COMP);
  }

  l_range = Lmax-Lmin;
  w_geom_c  = width_c*yheight*1.0e4;     /* source area correction */
  w_geom_t  = width_t*yheight*1.0e4;
  w_mult  = l_range;            /* wavelength range correction */
  w_mult *= 1.0/mcget_ncount();   /* Correct for number of rays */
  w_mult *= nu;               /* Correct for frequency */

  /* Calculate location of thermal wings wrt beamport_angle (z) direction */
  /* Wing 1 (left) is at -beamport_angle */
  t1z = cyl_radius*cos(-DEG2RAD*beamport_angle);
  t1x = cyl_radius*sin(-DEG2RAD*beamport_angle);
  t1y = 0;
  /* Wing 2 (right) is at 60-beamport_angle */
  t2z = cyl_radius*cos(DEG2RAD*(60-beamport_angle));
  t2x = cyl_radius*sin(DEG2RAD*(60-beamport_angle));
  t2y = 0;
  /* We want unit vectors... */
  NORM(t1x,t1y,t1z);
  NORM(t2x,t2y,t2z);

  #undef width_c
  #undef yheight
  #undef Lmin
  #undef Lmax
  #undef dist
  #undef focus_xw
  #undef focus_yh
  #undef nu
  #undef T
  #undef tau
  #undef tau1
  #undef tau2
  #undef d
  #undef n
  #undef cold_frac
  #undef n2
  #undef chi2
  #undef I0
  #undef I2
  #undef target_index
  #undef cyl_radius
  #undef branch1
  #undef branch2
  #undef branch_tail
  #undef n_pulses
  #undef width_t
  #undef T_t
  #undef tau_t
  #undef tau1_t
  #undef tau2_t
  #undef chi2_t
  #undef I0_t
  #undef I2_t
  #undef branch1_t
  #undef branch2_t
  #undef src_2012
  #undef tfocus_dist
  #undef tfocus_time
  #undef tfocus_width
  #undef beamport_angle
  #undef l_range
  #undef w_mult
  #undef w_geom
  #undef w_geom_c
  #undef w_geom_t
  #undef tx
  #undef ty
  #undef tz
  #undef t1x
  #undef t1y
  #undef t1z
  #undef t2x
  #undef t2y
  #undef t2z
  #undef T_n
  #undef tau_n
  #undef tau1_n
  #undef tau2_n
  #undef chi2_n
  #undef I0_n
  #undef I2_n
  #undef branch1_n
  #undef branch2_n
  #undef r_empty
  #undef r_optics
  return(_comp);
} /* class_ESS_moderator_long_init */

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

_class_Guide *class_Guide_init(_class_Guide *_comp
) {
  #define reflect (_comp->_parameters.reflect)
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
  #define pTable (_comp->_parameters.pTable)
if (mcgravitation) fprintf(stderr,"WARNING: Guide: %s: "
    "This component produces wrong results with gravitation !\n"
    "Use Guide_gravity.\n",
    NAME_CURRENT_COMP);

  if (!w2) w2=w1;
  if (!h2) h2=h1;

  if (reflect && strlen(reflect) && strcmp(reflect,"NULL") && strcmp(reflect,"0")) {
    if (Table_Read(&pTable, reflect, 1) <= 0) /* read 1st block data from file into pTable */
      exit(fprintf(stderr,"Guide: %s: can not read file %s\n", NAME_CURRENT_COMP, reflect));
  } else {
    if (W < 0 || R0 < 0 || Qc < 0 || m < 0)
    { fprintf(stderr,"Guide: %s: W R0 Qc must be >0.\n", NAME_CURRENT_COMP);
      exit(-1); }
  }
  #undef reflect
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
  #undef pTable
  return(_comp);
} /* class_Guide_init */

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

_class_TOFLambda_monitor *class_TOFLambda_monitor_init(_class_TOFLambda_monitor *_comp
) {
  #define nL (_comp->_parameters.nL)
  #define nt (_comp->_parameters.nt)
  #define tmin (_comp->_parameters.tmin)
  #define tmax (_comp->_parameters.tmax)
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
  #define tt_0 (_comp->_parameters.tt_0)
  #define tt_1 (_comp->_parameters.tt_1)
  #define TOFL_N (_comp->_parameters.TOFL_N)
  #define TOFL_p (_comp->_parameters.TOFL_p)
  #define TOFL_p2 (_comp->_parameters.TOFL_p2)
  int i,j;

  if (xwidth  > 0) { xmax = xwidth/2;  xmin = -xmax; }
  if (yheight > 0) { ymax = yheight/2; ymin = -ymax; }

  if ((xmin >= xmax) || (ymin >= ymax)) {
    printf("TOFlambda_monitor: %s: Null detection area !\n"
           "ERROR              (xwidth,yheight,xmin,xmax,ymin,ymax). Exiting",
           NAME_CURRENT_COMP);
    exit(0);
  }

  TOFL_N = create_darr2d(nt, nL);
  TOFL_p = create_darr2d(nt, nL);
  TOFL_p2 = create_darr2d(nt, nL);
  tt_0 = tmin*1e-6;
  tt_1 = tmax*1e-6;
  for (i=0; i<nL; i++){
    for (j=0; j<nt; j++){
      TOFL_N[j][i] = 0;
      TOFL_p[j][i] = 0;
      TOFL_p2[j][i] = 0;
    }
  }
  #undef nL
  #undef nt
  #undef tmin
  #undef tmax
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
  #undef tt_0
  #undef tt_1
  #undef TOFL_N
  #undef TOFL_p
  #undef TOFL_p2
  return(_comp);
} /* class_TOFLambda_monitor_init */

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

_class_Tunneling_sample *class_Tunneling_sample_init(_class_Tunneling_sample *_comp
) {
  #define thickness (_comp->_parameters.thickness)
  #define radius (_comp->_parameters.radius)
  #define focus_r (_comp->_parameters.focus_r)
  #define p_interact (_comp->_parameters.p_interact)
  #define f_QE (_comp->_parameters.f_QE)
  #define f_tun (_comp->_parameters.f_tun)
  #define gamma (_comp->_parameters.gamma)
  #define E_tun (_comp->_parameters.E_tun)
  #define target_x (_comp->_parameters.target_x)
  #define target_y (_comp->_parameters.target_y)
  #define target_z (_comp->_parameters.target_z)
  #define focus_xw (_comp->_parameters.focus_xw)
  #define focus_yh (_comp->_parameters.focus_yh)
  #define focus_aw (_comp->_parameters.focus_aw)
  #define focus_ah (_comp->_parameters.focus_ah)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define zdepth (_comp->_parameters.zdepth)
  #define sigma_abs (_comp->_parameters.sigma_abs)
  #define sigma_inc (_comp->_parameters.sigma_inc)
  #define Vc (_comp->_parameters.Vc)
  #define target_index (_comp->_parameters.target_index)
  #define ftun (_comp->_parameters.ftun)
  #define fQE (_comp->_parameters.fQE)
  #define VarsV (_comp->_parameters.VarsV)
  if (!xwidth || !yheight || !zdepth) /* Cannot define a rectangle */
    if (!radius || !yheight)              /* Cannot define a cylinder either */
      exit(fprintf(stderr,"V_sample: %s: sample has no volume (zero dimensions)\n", NAME_CURRENT_COMP));
    else                              /* It is a cylinder */
      VarsV.isrect=0;
  else                                /* It is a rectangle */
    VarsV.isrect=1;

/*  if(VarsV.isrect) printf("isrect"); else printf("not isrect"); */

  VarsV.sigma_a=sigma_abs;
  VarsV.sigma_i=sigma_inc;
  VarsV.rho = (1/Vc);
  VarsV.my_s=(VarsV.rho * 100 * VarsV.sigma_i);
  VarsV.my_a_v=(VarsV.rho * 100 * VarsV.sigma_a);

  /* now compute target coords if a component index is supplied */
  VarsV.tx= VarsV.ty=VarsV.tz=0;
  if (target_index)
  {
    Coords ToTarget;
    ToTarget = coords_sub(POS_A_COMP_INDEX(INDEX_CURRENT_COMP+target_index),POS_A_CURRENT_COMP);
    ToTarget = rot_apply(ROT_A_CURRENT_COMP, ToTarget);
    coords_get(ToTarget, &VarsV.tx, &VarsV.ty, &VarsV.tz);
  }
  else
  { VarsV.tx = target_x; VarsV.ty = target_y; VarsV.tz = target_z; }

  if (!(VarsV.tx || VarsV.ty || VarsV.tz))
    printf("Tunneling_sample: %s: The target is not defined. Using direct beam (Z-axis).\n",
      NAME_CURRENT_COMP);

  VarsV.distance=sqrt(VarsV.tx*VarsV.tx+VarsV.ty*VarsV.ty+VarsV.tz*VarsV.tz);

  /* different ways of setting rectangular area */
  VarsV.aw  = VarsV.ah = 0;
  if (focus_xw) {
  VarsV.xw = focus_xw;
  }
  if (focus_yh) {
    VarsV.yh = focus_yh;
  }
  if (focus_aw) {
    VarsV.aw = DEG2RAD*focus_aw;
  }
  if (focus_ah) {
    VarsV.ah = DEG2RAD*focus_ah;
  }

  /* Check that probabilities are positive and do not exceed unity */
  if (f_tun<0)
    ftun=0;
  else
    ftun=f_tun;
  if(f_QE<0)
    fQE=0;
  else
    fQE=f_QE;
  //  printf("Tunneling_sample: ftun %g %g fQE %g %g \n",f_tun,ftun,f_QE,fQE);
  if ((ftun+fQE)>1) {
    ftun=0;
    printf("Tunneling_sample: Sum of inelastic probabilities > 1. Setting f_tun=0");
    if (fQE>1) {
      fQE=0;
      printf("Tunneling_sample: Probability fQE > 1. Setting fQE=0.");
    }
  }
  #undef thickness
  #undef radius
  #undef focus_r
  #undef p_interact
  #undef f_QE
  #undef f_tun
  #undef gamma
  #undef E_tun
  #undef target_x
  #undef target_y
  #undef target_z
  #undef focus_xw
  #undef focus_yh
  #undef focus_aw
  #undef focus_ah
  #undef xwidth
  #undef yheight
  #undef zdepth
  #undef sigma_abs
  #undef sigma_inc
  #undef Vc
  #undef target_index
  #undef ftun
  #undef fQE
  #undef VarsV
  return(_comp);
} /* class_Tunneling_sample_init */

_class_TOF2E_monitor *class_TOF2E_monitor_init(_class_TOF2E_monitor *_comp
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
  #define T_zero (_comp->_parameters.T_zero)
  #define L_flight (_comp->_parameters.L_flight)
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
  #undef T_zero
  #undef L_flight
  #undef restore_neutron
  #undef E_N
  #undef E_p
  #undef E_p2
  #undef S_p
  #undef S_pE
  #undef S_pE2
  return(_comp);
} /* class_TOF2E_monitor_init */



int init(void) { /* called by mccode_main for ESS_IN5_reprate:INITIALISE */
  DEBUG_INSTR();

  /* code_main/parseoptions/readparams sets instrument parameters value */
  stracpy(instrument->_name, "ESS_IN5_reprate", 256);

  /* Instrument 'ESS_IN5_reprate' INITIALISE */
  SIG_MESSAGE("[ESS_IN5_reprate] INITIALISE [ESS_IN5_reprate.instr:95]");
  #define Lmin (instrument->_parameters._Lmin)
  #define Lmax (instrument->_parameters._Lmax)
  #define lambda0 (instrument->_parameters._lambda0)
  #define Pulse_width (instrument->_parameters._Pulse_width)
  #define Num_pulses (instrument->_parameters._Num_pulses)
  #define GUI_start (instrument->_parameters._GUI_start)
  #define FO1_DIST (instrument->_parameters._FO1_DIST)
  #define L_ballistic_begin (instrument->_parameters._L_ballistic_begin)
  #define L_ballistic_end (instrument->_parameters._L_ballistic_end)
  #define Length (instrument->_parameters._Length)
  #define SAMPLE_DIST (instrument->_parameters._SAMPLE_DIST)
  #define DETECTOR_DIST (instrument->_parameters._DETECTOR_DIST)
  #define GUI_h (instrument->_parameters._GUI_h)
  #define GUI_w (instrument->_parameters._GUI_w)
  #define GUI_GAP (instrument->_parameters._GUI_GAP)
  #define H1 (instrument->_parameters._H1)
  #define W1 (instrument->_parameters._W1)
  #define H2 (instrument->_parameters._H2)
  #define W2 (instrument->_parameters._W2)
  #define H3 (instrument->_parameters._H3)
  #define W3 (instrument->_parameters._W3)
  #define H4 (instrument->_parameters._H4)
  #define W4 (instrument->_parameters._W4)
  #define H_chop (instrument->_parameters._H_chop)
  #define W_chop (instrument->_parameters._W_chop)
  #define H_end (instrument->_parameters._H_end)
  #define W_end (instrument->_parameters._W_end)
  #define ALPHA (instrument->_parameters._ALPHA)
  #define M (instrument->_parameters._M)
  #define F_slow1 (instrument->_parameters._F_slow1)
  #define F_slow2 (instrument->_parameters._F_slow2)
  #define F_fast1 (instrument->_parameters._F_fast1)
  #define F_fast2 (instrument->_parameters._F_fast2)
  #define N_fast (instrument->_parameters._N_fast)
  #define SLOW1_THETA (instrument->_parameters._SLOW1_THETA)
  #define FO3 (instrument->_parameters._FO3)
  #define THETA_fast1 (instrument->_parameters._THETA_fast1)
  #define FAST_THETA (instrument->_parameters._FAST_THETA)
  #define Gamma (instrument->_parameters._Gamma)
  #define Etun (instrument->_parameters._Etun)
  #define V_HOLE (instrument->_parameters._V_HOLE)
  #define FRAC_QUASIEL (instrument->_parameters._FRAC_QUASIEL)
  #define FRAC_TUNNEL (instrument->_parameters._FRAC_TUNNEL)
  #define TT (instrument->_parameters._TT)
  #define RES_DE (instrument->_parameters._RES_DE)
  #define port (instrument->_parameters._port)
  #define cold (instrument->_parameters._cold)
{
        FREQ = 1000.0/60.0;
        t_offset=Pulse_width/2.0+170e-6;
        t_FO1 = FO1_DIST*lambda0/(2*PI*K2V)+t_offset;
        t_FO2 = Length/2*lambda0/(2*PI*K2V)+t_offset;
        t_fast1 = (Length/2+0.1)*lambda0/(2*PI*K2V)+t_offset;
        t_fast2 = Length*lambda0/(2*PI*K2V)+t_offset;
        t_fast2a = (Length+GUI_GAP)*lambda0/(2*PI*K2V)+t_offset;
        t_fast3 = (Length+2*GUI_GAP)*lambda0/(2*PI*K2V)+t_offset;
        t_sample = (Length+SAMPLE_DIST)*lambda0/(2*PI*K2V)+t_offset;
        t_detector = (Length+SAMPLE_DIST+DETECTOR_DIST)*lambda0/(2*PI*K2V)+t_offset;
        tmin_zoom = t_fast2*1e6-1e3;
        tmax_zoom = t_fast2*1e6+1e3;
        E_target = VS2E*(K2V*2*PI/lambda0)*(K2V*2*PI/lambda0);
}
  #undef Lmin
  #undef Lmax
  #undef lambda0
  #undef Pulse_width
  #undef Num_pulses
  #undef GUI_start
  #undef FO1_DIST
  #undef L_ballistic_begin
  #undef L_ballistic_end
  #undef Length
  #undef SAMPLE_DIST
  #undef DETECTOR_DIST
  #undef GUI_h
  #undef GUI_w
  #undef GUI_GAP
  #undef H1
  #undef W1
  #undef H2
  #undef W2
  #undef H3
  #undef W3
  #undef H4
  #undef W4
  #undef H_chop
  #undef W_chop
  #undef H_end
  #undef W_end
  #undef ALPHA
  #undef M
  #undef F_slow1
  #undef F_slow2
  #undef F_fast1
  #undef F_fast2
  #undef N_fast
  #undef SLOW1_THETA
  #undef FO3
  #undef THETA_fast1
  #undef FAST_THETA
  #undef Gamma
  #undef Etun
  #undef V_HOLE
  #undef FRAC_QUASIEL
  #undef FRAC_TUNNEL
  #undef TT
  #undef RES_DE
  #undef port
  #undef cold
  _source_setpos(); /* type ESS_moderator_long */
  _Origin_setpos(); /* type Progress_bar */
  _TOFmoderator_zoom_setpos(); /* type TOF_monitor */
  _TOFmoderator_setpos(); /* type TOF_monitor */
  _Lmon_guistart_setpos(); /* type L_monitor */
  _Lmon_normalize_setpos(); /* type L_monitor */
  _Guide1_setpos(); /* type Guide */
  _Lmonslow1_setpos(); /* type L_monitor */
  _PSDslow1_setpos(); /* type PSD_monitor */
  _FOchop1_setpos(); /* type DiskChopper */
  _TOFLmon1_setpos(); /* type TOFLambda_monitor */
  _Lmon_afterslow1_setpos(); /* type L_monitor */
  _PSD_afterslow1_setpos(); /* type PSD_monitor */
  _Guidelong1_setpos(); /* type Guide */
  _Guidelong1b_setpos(); /* type Guide */
  _Lmon_slow2_setpos(); /* type L_monitor */
  _FOchop2_setpos(); /* type DiskChopper */
  _Fastchop1_setpos(); /* type DiskChopper */
  _PSD_afterslow2_setpos(); /* type PSD_monitor */
  _Lmon_afterslow2_setpos(); /* type L_monitor */
  _TOFL_afterslow2_setpos(); /* type TOFLambda_monitor */
  _Guidelong2_setpos(); /* type Guide */
  _Lmon_beforeballistic_setpos(); /* type L_monitor */
  _PSD_beforeballistic_setpos(); /* type PSD_monitor */
  _Guidelong2a_setpos(); /* type Guide */
  _Lmonfast2_setpos(); /* type L_monitor */
  _Lmonfast2_zoom_setpos(); /* type L_monitor */
  _TOFLfast2_setpos(); /* type TOFLambda_monitor */
  _TOFLfast2zoom_setpos(); /* type TOFLambda_monitor */
  _PSDfast2_setpos(); /* type PSD_monitor */
  _Fastchop2_setpos(); /* type DiskChopper */
  _Fastchop2counter_setpos(); /* type DiskChopper */
  _FOchop3_setpos(); /* type DiskChopper */
  _TOFfast2_zoom_setpos(); /* type TOF_monitor */
  _Lmon_afterfast2_setpos(); /* type L_monitor */
  _TOFL_afterfast2_setpos(); /* type TOFLambda_monitor */
  _TOFL_afterfast2_zoom_setpos(); /* type TOFLambda_monitor */
  _PSD_afterfast2_setpos(); /* type PSD_monitor */
  _Guidesample_setpos(); /* type Guide */
  _Lmon_guideend_setpos(); /* type L_monitor */
  _PSDsample_setpos(); /* type PSD_monitor */
  _TOFsample_zoom_setpos(); /* type TOF_monitor */
  _Esample_setpos(); /* type E_monitor */
  _Lmon_sample_zoom_setpos(); /* type L_monitor */
  _sample_setpos(); /* type Tunneling_sample */
  _detectorarm_setpos(); /* type Arm */
  _TOFdetector_setpos(); /* type TOF_monitor */
  _TOFdetector_zoom_setpos(); /* type TOF_monitor */
  _Edetector_setpos(); /* type E_monitor */
  _TOF2Edetector_setpos(); /* type TOF2E_monitor */

  /* call iteratively all components INITIALISE */
  class_ESS_moderator_long_init(_source);

  class_Progress_bar_init(_Origin);

  class_TOF_monitor_init(_TOFmoderator_zoom);

  class_TOF_monitor_init(_TOFmoderator);

  class_L_monitor_init(_Lmon_guistart);

  class_L_monitor_init(_Lmon_normalize);

  class_Guide_init(_Guide1);

  class_L_monitor_init(_Lmonslow1);

  class_PSD_monitor_init(_PSDslow1);

  class_DiskChopper_init(_FOchop1);

  class_TOFLambda_monitor_init(_TOFLmon1);

  class_L_monitor_init(_Lmon_afterslow1);

  class_PSD_monitor_init(_PSD_afterslow1);

  class_Guide_init(_Guidelong1);

  class_Guide_init(_Guidelong1b);

  class_L_monitor_init(_Lmon_slow2);

  class_DiskChopper_init(_FOchop2);

  class_DiskChopper_init(_Fastchop1);

  class_PSD_monitor_init(_PSD_afterslow2);

  class_L_monitor_init(_Lmon_afterslow2);

  class_TOFLambda_monitor_init(_TOFL_afterslow2);

  class_Guide_init(_Guidelong2);

  class_L_monitor_init(_Lmon_beforeballistic);

  class_PSD_monitor_init(_PSD_beforeballistic);

  class_Guide_init(_Guidelong2a);

  class_L_monitor_init(_Lmonfast2);

  class_L_monitor_init(_Lmonfast2_zoom);

  class_TOFLambda_monitor_init(_TOFLfast2);

  class_TOFLambda_monitor_init(_TOFLfast2zoom);

  class_PSD_monitor_init(_PSDfast2);

  class_DiskChopper_init(_Fastchop2);

  class_DiskChopper_init(_Fastchop2counter);

  class_DiskChopper_init(_FOchop3);

  class_TOF_monitor_init(_TOFfast2_zoom);

  class_L_monitor_init(_Lmon_afterfast2);

  class_TOFLambda_monitor_init(_TOFL_afterfast2);

  class_TOFLambda_monitor_init(_TOFL_afterfast2_zoom);

  class_PSD_monitor_init(_PSD_afterfast2);

  class_Guide_init(_Guidesample);

  class_L_monitor_init(_Lmon_guideend);

  class_PSD_monitor_init(_PSDsample);

  class_TOF_monitor_init(_TOFsample_zoom);

  class_E_monitor_init(_Esample);

  class_L_monitor_init(_Lmon_sample_zoom);

  class_Tunneling_sample_init(_sample);


  class_TOF_monitor_init(_TOFdetector);

  class_TOF_monitor_init(_TOFdetector_zoom);

  class_E_monitor_init(_Edetector);

  class_TOF2E_monitor_init(_TOF2Edetector);

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
_class_ESS_moderator_long *class_ESS_moderator_long_trace(_class_ESS_moderator_long *_comp
  , _class_particle *_particle) {
  ABSORBED=SCATTERED=RESTORE=0;

  #define width_c (_comp->_parameters.width_c)
  #define yheight (_comp->_parameters.yheight)
  #define Lmin (_comp->_parameters.Lmin)
  #define Lmax (_comp->_parameters.Lmax)
  #define dist (_comp->_parameters.dist)
  #define focus_xw (_comp->_parameters.focus_xw)
  #define focus_yh (_comp->_parameters.focus_yh)
  #define nu (_comp->_parameters.nu)
  #define T (_comp->_parameters.T)
  #define tau (_comp->_parameters.tau)
  #define tau1 (_comp->_parameters.tau1)
  #define tau2 (_comp->_parameters.tau2)
  #define d (_comp->_parameters.d)
  #define n (_comp->_parameters.n)
  #define cold_frac (_comp->_parameters.cold_frac)
  #define n2 (_comp->_parameters.n2)
  #define chi2 (_comp->_parameters.chi2)
  #define I0 (_comp->_parameters.I0)
  #define I2 (_comp->_parameters.I2)
  #define target_index (_comp->_parameters.target_index)
  #define cyl_radius (_comp->_parameters.cyl_radius)
  #define branch1 (_comp->_parameters.branch1)
  #define branch2 (_comp->_parameters.branch2)
  #define branch_tail (_comp->_parameters.branch_tail)
  #define n_pulses (_comp->_parameters.n_pulses)
  #define width_t (_comp->_parameters.width_t)
  #define T_t (_comp->_parameters.T_t)
  #define tau_t (_comp->_parameters.tau_t)
  #define tau1_t (_comp->_parameters.tau1_t)
  #define tau2_t (_comp->_parameters.tau2_t)
  #define chi2_t (_comp->_parameters.chi2_t)
  #define I0_t (_comp->_parameters.I0_t)
  #define I2_t (_comp->_parameters.I2_t)
  #define branch1_t (_comp->_parameters.branch1_t)
  #define branch2_t (_comp->_parameters.branch2_t)
  #define src_2012 (_comp->_parameters.src_2012)
  #define tfocus_dist (_comp->_parameters.tfocus_dist)
  #define tfocus_time (_comp->_parameters.tfocus_time)
  #define tfocus_width (_comp->_parameters.tfocus_width)
  #define beamport_angle (_comp->_parameters.beamport_angle)
  #define l_range (_comp->_parameters.l_range)
  #define w_mult (_comp->_parameters.w_mult)
  #define w_geom (_comp->_parameters.w_geom)
  #define w_geom_c (_comp->_parameters.w_geom_c)
  #define w_geom_t (_comp->_parameters.w_geom_t)
  #define tx (_comp->_parameters.tx)
  #define ty (_comp->_parameters.ty)
  #define tz (_comp->_parameters.tz)
  #define t1x (_comp->_parameters.t1x)
  #define t1y (_comp->_parameters.t1y)
  #define t1z (_comp->_parameters.t1z)
  #define t2x (_comp->_parameters.t2x)
  #define t2y (_comp->_parameters.t2y)
  #define t2z (_comp->_parameters.t2z)
  #define T_n (_comp->_parameters.T_n)
  #define tau_n (_comp->_parameters.tau_n)
  #define tau1_n (_comp->_parameters.tau1_n)
  #define tau2_n (_comp->_parameters.tau2_n)
  #define chi2_n (_comp->_parameters.chi2_n)
  #define I0_n (_comp->_parameters.I0_n)
  #define I2_n (_comp->_parameters.I2_n)
  #define branch1_n (_comp->_parameters.branch1_n)
  #define branch2_n (_comp->_parameters.branch2_n)
  #define r_empty (_comp->_parameters.r_empty)
  #define r_optics (_comp->_parameters.r_optics)
  double v,tau_l,E,lambda,k,r,xf,yf,dx,dy,w_focus,tail_flag,cor,dt,xprime,yprime,zprime;

  /* Bispectral source - choice of spectrum and initial position */
  int cold = ( rand01() < cold_frac ) ? 1 : 0;

  /* Geometry adapted from ESS MCNPX model, mid 2012 */
  if (cold) {          //case: cold moderator
    double theta_tmp;

    //choose random point on cylinder surface
    theta_tmp = randpm1()*PI/6 + PI/2 + (- 30 + beamport_angle)*DEG2RAD;
        x     = cyl_radius * cos(theta_tmp);
        y     = 0.5*randpm1()*yheight;
        z     = cyl_radius * sin(theta_tmp);

    //spectrum related constants - ESS 2001 Cold moderator
    //T=50, tau=287e-6, tau1=0, tau2=20e-6, chi2=0.9, I0=6.9e11, I2=27.6e10, branch1=0, branch2=0.5;
    T_n=T; tau_n=tau; tau1_n=tau1; tau2_n=tau2; chi2_n=chi2; I0_n=I0; I2_n=I2; branch1_n=branch1; branch2_n=branch2;
    w_geom = w_geom_c;
  }
  else //case: thermal moderator
  {
    /* choose "left" or "right" thermal wing */
    int isleft = ( rand01() < 0.5 ) ? 1 : 0;
    double poshorz, posvert;

    poshorz = cyl_radius+rand01()*width_t;
    posvert = 0.5*randpm1()*yheight;

    if (isleft) {
      x = t1x * poshorz;
      z = t1z * poshorz;
    } else {
      x = t2x * poshorz;
      z = t2z * poshorz;
    }
    y = posvert;

    /* x = cyl_radius + width_t*rand01(); */
    /* y = 0.5*randpm1()*width_t; */
    /* z = cyl_radius; */

    //spectrum related constants - ESS 2001 Thermal moderator
    //T_t=325, tau_t=80e-6, tau1_t=400e-6, tau2_t=12e-6, chi2_t=2.5, I0_t=13.5e11, I2_t=27.6e10, branch1_t=0.5, branch2_t=0.5;
    T_n=T_t; tau_n=tau_t; tau1_n=tau1_t; tau2_n=tau2_t; chi2_n=chi2_t; I0_n=I0_t; I2_n=I2_t; branch1_n=branch1_t; branch2_n=branch2_t;
    w_geom = w_geom_t;
  }

  randvec_target_rect_real(&xf, &yf, &r, &w_focus, tx, ty, tz, focus_xw, focus_yh, ROT_A_CURRENT_COMP, x, y, z, 2);

  dx = xf-x;
  dy = yf-y;
  r = sqrt(dx*dx+dy*dy+dist*dist);

  lambda = Lmin+l_range*rand01();    /* Choose from uniform distribution */
  k = 2*PI/lambda;
  v = K2V*k;

  vz = v*dist/r;
  vy = v*dy/r;
  vx = v*dx/r;

  /* Determine delta-t needed to reach first chopper */
  if (tfocus_width>0) {
    dt = tfocus_dist/vz;			/* Flight time to time window (chopper) */
  }
  tail_flag = (rand01()<branch_tail);   /* Choose tail/bulk */
  if (tail_flag)
  {
    if (rand01() < branch2_n)
    {
      if (tau1_n>0)
      {
        if (rand01() < branch1_n)     /* Quick and dirty non-general solution */
        {  /* FIRST CASE a */
          tau_l = tau_n;
          p = 1/(branch1_n*branch2_n*branch_tail); /* Correct for switching prob. */
        }
        else
        {  /* FIRST CASE b */
          tau_l = tau1_n;
          p = 1/((1-branch1_n)*branch2_n*branch_tail); /* Correct for switching prob. */
        }
      }
      else
      {
        tau_l = tau_n;
        p = 1/(branch2_n*branch_tail); /* Correct for switching prob. */
      }
      t = -tau_l*log(1e-12+rand01());       /* Sample from long-time tail a */
      /* Correct for true pulse shape */
      p *= w_focus;                         /* Correct for target focusing */
      p *= tau_l/d;                         /* Correct for tail part */
      p *= I0_n*w_mult*w_geom*Mezei_M_fct(lambda,T_n);           /* Calculate true intensity */
    }
    else
    {
      /* SECOND CASE */
      tau_l = tau2_n*lambda;
      t = -tau_l*log(1e-12+rand01());       /* Sample from long-time tail */
      p = n2/(n2-1)*((1-exp(-d/tau_l))-(1-exp(-n2*d/tau_l))*exp(-(n2-1)*t/tau_l)/n);
                                            /* Correct for true pulse shape */
      p /= (1-branch2_n)*branch_tail;          /* Correct for switching prob. */
      p *= tau_l/d;                         /* Correct for tail part */
      p *= w_focus;                         /* Correct for target focusing */
      p *= I2_n*w_mult*w_geom/(1+exp(chi2_n*lambda-2.2))/lambda; /* Calculate true intensity */
    }
    t += d;                                 /* Add pulse length */
  }
  else /* Tail-flag */
  {
    if (tfocus_width>0)
    {
      t = tfocus_time-dt;                    /* Set time to hit time window center */
      t += randpm1()*tfocus_width/2.0;       /* Add random time within window width */
    }
    else
    {
      t = d*rand01();                        /* Sample from bulk pulse */
    }
    if (t<0) ABSORB;                       /* Kill neutron if outside pulse duration */
    if (t>d) ABSORB;
    if (rand01() < branch2_n)
    {
      if (rand01() < branch1_n)     /* Quick and dirty non-general solution */
      {  /* FIRST CASE a */
        tau_l = tau_n;
        p = 1/(branch1_n*branch2_n*(1-branch_tail)); /* Correct for switching prob. */
      }
      else
      {  /* FIRST CASE b */
        tau_l = tau1_n;
        p = 1/((1-branch1_n)*branch2_n*(1-branch_tail)); /* Correct for switching prob. */
      }
      p *= 1-n/(n-1)*(exp(-t/tau_l)-exp(-n*t/tau_l)/n); /* Correct for true pulse shape */
      p *= w_focus;                         /* Correct for target focusing */
      if (tfocus_width>0) {
        p *= tfocus_width/d;    	  	  /* Correct for time focusing */
      }
      p *= I0_n*w_mult*w_geom*Mezei_M_fct(lambda,T_n);       /* Calculate true intensity */
    }
    else
    {
      /* SECOND CASE */
      tau_l = tau2_n*lambda;
      p = 1-n2/(n2-1)*(exp(-t/tau_l)-exp(-n2*t/tau_l)/n2); /* Correct for true pulse shape */
      p /= (1-branch2_n)*(1-branch_tail);   /* Correct for switching prob. */
      p *= w_focus;                         /* Correct for target focusing */
      if (tfocus_width) {
        p *= tfocus_width/d;    		  /* Correct for time focusing */
      }
      p *= I2_n*w_mult*w_geom/(1+exp(chi2_n*lambda-2.2))/lambda;    /* Calculate true intensity */
    }
  }

  if (cold && src_2012)
  {
    /* Correction factors to converts 'predicted' spectrum from cold moderator to the one observed in MCNPX */
    if (lambda<=2.5)
      cor=log(1.402+0.898*lambda)*(2.0776-4.1093*lambda+4.8836*pow(lambda,2)-2.4715*pow(lambda,3)+0.4521*pow(lambda,4));
    else if (lambda <= 3.5)
      cor = log(1.402 + 0.898*lambda)*(4.3369 - 1.8367*lambda + 0.2524*pow(lambda,2) );
    else if (lambda  > 3.5)
      cor = log(1.402 + 0.898*lambda);
  }
  else
  {
    /* Thermal (pre-)moderator, i.e. no correction */
    cor = 1.0;
  }
  p *= cor;

  /* Correct weight for sampling of cold vs. thermal events. */
  if (cold) {
    p /=cold_frac;
  } else {
    p/=(1-cold_frac);
  }

  t+=(double)floor((n_pulses)*rand01())/nu;   /* Select a random pulse */
  #undef width_c
  #undef yheight
  #undef Lmin
  #undef Lmax
  #undef dist
  #undef focus_xw
  #undef focus_yh
  #undef nu
  #undef T
  #undef tau
  #undef tau1
  #undef tau2
  #undef d
  #undef n
  #undef cold_frac
  #undef n2
  #undef chi2
  #undef I0
  #undef I2
  #undef target_index
  #undef cyl_radius
  #undef branch1
  #undef branch2
  #undef branch_tail
  #undef n_pulses
  #undef width_t
  #undef T_t
  #undef tau_t
  #undef tau1_t
  #undef tau2_t
  #undef chi2_t
  #undef I0_t
  #undef I2_t
  #undef branch1_t
  #undef branch2_t
  #undef src_2012
  #undef tfocus_dist
  #undef tfocus_time
  #undef tfocus_width
  #undef beamport_angle
  #undef l_range
  #undef w_mult
  #undef w_geom
  #undef w_geom_c
  #undef w_geom_t
  #undef tx
  #undef ty
  #undef tz
  #undef t1x
  #undef t1y
  #undef t1z
  #undef t2x
  #undef t2y
  #undef t2z
  #undef T_n
  #undef tau_n
  #undef tau1_n
  #undef tau2_n
  #undef chi2_n
  #undef I0_n
  #undef I2_n
  #undef branch1_n
  #undef branch2_n
  #undef r_empty
  #undef r_optics
  return(_comp);
} /* class_ESS_moderator_long_trace */

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
_class_Guide *class_Guide_trace(_class_Guide *_comp
  , _class_particle *_particle) {
  ABSORBED=SCATTERED=RESTORE=0;

  #define reflect (_comp->_parameters.reflect)
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
  #define pTable (_comp->_parameters.pTable)
  double t1,t2;                                 /* Intersection times. */
  double av,ah,bv,bh,cv1,cv2,ch1,ch2,d;         /* Intermediate values */
  double weight;                                /* Internal probability weight */
  double vdotn_v1,vdotn_v2,vdotn_h1,vdotn_h2;   /* Dot products. */
  int i;                                        /* Which mirror hit? */
  double q;                                     /* Q [1/AA] of reflection */
  double nlen2;                                 /* Vector lengths squared */

  /* ToDo: These could be precalculated. */
  double ww = .5*(w2 - w1), hh = .5*(h2 - h1);
  double whalf = .5*w1, hhalf = .5*h1;

  /* Propagate neutron to guide entrance. */
  PROP_Z0;
  /* Scatter here to ensure that fully transmitted neutrons will not be
     absorbed in a GROUP construction, e.g. all neutrons - even the
     later absorbed ones are scattered at the guide entry. */
  SCATTER;
  if(x <= -whalf || x >= whalf || y <= -hhalf || y >= hhalf)
    ABSORB;
  for(;;)
  {
    /* Compute the dot products of v and n for the four mirrors. */
    av = l*vx; bv = ww*vz;
    ah = l*vy; bh = hh*vz;
    vdotn_v1 = bv + av;         /* Left vertical */
    vdotn_v2 = bv - av;         /* Right vertical */
    vdotn_h1 = bh + ah;         /* Lower horizontal */
    vdotn_h2 = bh - ah;         /* Upper horizontal */
    /* Compute the dot products of (O - r) and n as c1+c2 and c1-c2 */
    cv1 = -whalf*l - z*ww; cv2 = x*l;
    ch1 = -hhalf*l - z*hh; ch2 = y*l;
    /* Compute intersection times. */
    t1 = (l - z)/vz;
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
      break;                    /* Neutron left guide. */
    PROP_DT(t1);
    switch(i)
    {
      case 1:                   /* Left vertical mirror */
        nlen2 = l*l + ww*ww;
        q = V2Q*(-2)*vdotn_v1/sqrt(nlen2);
        d = 2*vdotn_v1/nlen2;
        vx = vx - d*l;
        vz = vz - d*ww;
        break;
      case 2:                   /* Right vertical mirror */
        nlen2 = l*l + ww*ww;
        q = V2Q*(-2)*vdotn_v2/sqrt(nlen2);
        d = 2*vdotn_v2/nlen2;
        vx = vx + d*l;
        vz = vz - d*ww;
        break;
      case 3:                   /* Lower horizontal mirror */
        nlen2 = l*l + hh*hh;
        q = V2Q*(-2)*vdotn_h1/sqrt(nlen2);
        d = 2*vdotn_h1/nlen2;
        vy = vy - d*l;
        vz = vz - d*hh;
        break;
      case 4:                   /* Upper horizontal mirror */
        nlen2 = l*l + hh*hh;
        q = V2Q*(-2)*vdotn_h2/sqrt(nlen2);
        d = 2*vdotn_h2/nlen2;
        vy = vy + d*l;
        vz = vz - d*hh;
        break;
    }
    /* Now compute reflectivity. */
    weight = 1.0; /* Initial internal weight factor */
    if(m == 0)
      ABSORB;
    if (reflect && strlen(reflect) && strcmp(reflect,"NULL") && strcmp(reflect,"0"))
       TableReflecFunc(q, &pTable, &weight);
    else {
      double par[] = {R0, Qc, alpha, m, W};
      StdReflecFunc(q, par, &weight);
    }
    if (weight > 0)
      p *= weight;
    else ABSORB;
    SCATTER;
  }
  #undef reflect
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
  #undef pTable
  return(_comp);
} /* class_Guide_trace */

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
_class_TOFLambda_monitor *class_TOFLambda_monitor_trace(_class_TOFLambda_monitor *_comp
  , _class_particle *_particle) {
  ABSORBED=SCATTERED=RESTORE=0;

  #define nL (_comp->_parameters.nL)
  #define nt (_comp->_parameters.nt)
  #define tmin (_comp->_parameters.tmin)
  #define tmax (_comp->_parameters.tmax)
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
  #define tt_0 (_comp->_parameters.tt_0)
  #define tt_1 (_comp->_parameters.tt_1)
  #define TOFL_N (_comp->_parameters.TOFL_N)
  #define TOFL_p (_comp->_parameters.TOFL_p)
  #define TOFL_p2 (_comp->_parameters.TOFL_p2)
  int i,j;
  double div;
  double lambda;

  PROP_Z0;
  lambda = (2*PI/V2K)/sqrt(vx*vx + vy*vy + vz*vz);
  if (x>xmin && x<xmax && y>ymin && y<ymax && lambda > Lmin && lambda < Lmax){
    if (t < tt_1 && t > tt_0)
    {
      i = floor((lambda - Lmin)*nL/(Lmax - Lmin));
      j = floor((t-tt_0)*nt/(tt_1-tt_0));
      /*  printf("tt_0, tt_1, nt %g %g %i t j %g %i \n",tt_0,tt_1,nt,t,j);
      */        TOFL_N[j][i]++;
      TOFL_p[j][i] += p;
      TOFL_p2[j][i] += p*p;
    }
  }
  if (restore_neutron) {
    RESTORE_NEUTRON(INDEX_CURRENT_COMP, x, y, z, vx, vy, vz, t, sx, sy, sz, p);
  }
  #undef nL
  #undef nt
  #undef tmin
  #undef tmax
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
  #undef tt_0
  #undef tt_1
  #undef TOFL_N
  #undef TOFL_p
  #undef TOFL_p2
  return(_comp);
} /* class_TOFLambda_monitor_trace */

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

#pragma acc routine seq
_class_Tunneling_sample *class_Tunneling_sample_trace(_class_Tunneling_sample *_comp
  , _class_particle *_particle) {
  ABSORBED=SCATTERED=RESTORE=0;

  #define thickness (_comp->_parameters.thickness)
  #define radius (_comp->_parameters.radius)
  #define focus_r (_comp->_parameters.focus_r)
  #define p_interact (_comp->_parameters.p_interact)
  #define f_QE (_comp->_parameters.f_QE)
  #define f_tun (_comp->_parameters.f_tun)
  #define gamma (_comp->_parameters.gamma)
  #define E_tun (_comp->_parameters.E_tun)
  #define target_x (_comp->_parameters.target_x)
  #define target_y (_comp->_parameters.target_y)
  #define target_z (_comp->_parameters.target_z)
  #define focus_xw (_comp->_parameters.focus_xw)
  #define focus_yh (_comp->_parameters.focus_yh)
  #define focus_aw (_comp->_parameters.focus_aw)
  #define focus_ah (_comp->_parameters.focus_ah)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define zdepth (_comp->_parameters.zdepth)
  #define sigma_abs (_comp->_parameters.sigma_abs)
  #define sigma_inc (_comp->_parameters.sigma_inc)
  #define Vc (_comp->_parameters.Vc)
  #define target_index (_comp->_parameters.target_index)
  #define ftun (_comp->_parameters.ftun)
  #define fQE (_comp->_parameters.fQE)
  #define VarsV (_comp->_parameters.VarsV)
  double t0, t3;                /* Entry/exit time for outer cylinder */
  double t1, t2;                /* Entry/exit time for inner cylinder */
  double v;                     /* Neutron velocity */
  double dt0, dt1, dt2, dt;     /* Flight times through sample */
  double l_full;                /* Flight path length for non-scattered neutron */
  double l_i, l_o=0;            /* Flight path lenght in/out for scattered neutron */
  double my_a=0;                /* Velocity-dependent attenuation factor */
  double solid_angle=0;         /* Solid angle of target as seen from scattering point */
  double aim_x=0, aim_y=0, aim_z=1;   /* Position of target relative to scattering point */
  double v_i, v_f, E_i, E_f;    /* initial and final energies and velocities */
  double dE;                    /* Energy transfer */
  double scatt_choice;          /* Representing random choice of scattering type */
  int    intersect=0;

  if (VarsV.isrect)
    intersect = box_intersect(&t0, &t3, x, y, z, vx, vy, vz, xwidth, yheight, zdepth);
  else
    intersect = cylinder_intersect(&t0, &t3, x, y, z, vx, vy, vz, radius, yheight);
  if(intersect)
  {
    if(t0 < 0) ABSORB; /* we already passed the sample; this is illegal */
    /* Neutron enters at t=t0. */
    if(VarsV.isrect)
      t1 = t2 = t3;
    else
      if(!thickness || !cylinder_intersect(&t1, &t2, x, y, z, vx, vy, vz, radius-thickness, yheight))
        t1 = t2 = t3;

    dt0 = t1-t0;                /* Time in sample, ingoing */
    dt1 = t2-t1;                /* Time in hole */
    dt2 = t3-t2;                /* Time in sample, outgoing */
    v = sqrt(vx*vx + vy*vy + vz*vz);
    l_full = v * (dt0 + dt2);   /* Length of full path through sample */
    if (v) my_a = VarsV.my_a_v*(2200/v);

    if (p_interact >= 1 || rand01()<p_interact)          /* Scattering */
    {
      dt = rand01()*(dt0+dt2);    /* Time of scattering (relative to t0) */
      l_i = v*dt;                 /* Penetration in sample: scattering+abs */
      if (dt > dt0)
        dt += dt1;                /* jump to 2nd side of cylinder */

      PROP_DT(dt+t0);             /* Point of scattering */

      if ((VarsV.tx || VarsV.ty || VarsV.tz)) {
        aim_x = VarsV.tx-x;       /* Vector pointing at target (anal./det.) */
        aim_y = VarsV.ty-y;
        aim_z = VarsV.tz-z;
      }
      if(VarsV.aw && VarsV.ah) {
        randvec_target_rect_angular(&vx, &vy, &vz, &solid_angle,
          aim_x, aim_y, aim_z, VarsV.aw, VarsV.ah, ROT_A_CURRENT_COMP);
      } else if(VarsV.xw && VarsV.yh) {
        randvec_target_rect(&vx, &vy, &vz, &solid_angle,
          aim_x, aim_y, aim_z, VarsV.xw, VarsV.yh, ROT_A_CURRENT_COMP);
      } else {
        randvec_target_circle(&vx, &vy, &vz, &solid_angle, aim_x, aim_y, aim_z, focus_r);
      }
      NORM(vx, vy, vz);
      
      scatt_choice = rand01();  /* chooses type of scattering */
      v_i = v;                  /* Store initial velocity in case of inel. */
      E_i = VS2E*v_i*v_i;
      //      printf("Tunneling_sample: scatt_choice %g fQE %g ftun %g \n",scatt_choice, fQE, ftun);
      if (scatt_choice<(fQE+ftun))    /* Inelastic choices */
	{
          if (scatt_choice<fQE) /* Quasielastic */
          { dE = gamma*tan(PI/2*randpm1());
	  //            printf("Tunneling_sample: Inside Quasielastic code \n");
	  }
	  else
          { if (randpm1()>0)
	      dE = E_tun;
		else
              dE = -E_tun;
	  //	  printf("Tunneling_sample: Inside tunnel code \n");
	  }
          E_f = E_i + dE;
          if (E_f <= 0) 
            ABSORB;
	  v_f = SE2V*sqrt(E_f);
          v = v_f;
	  /*          printf("vi: %g Ei: %g dE: %g Ef %g vf: %g v: %g \n",
		      v_i,E_i,dE,E_f,v_f,v); */
	}

      vx *= v;
      vy *= v;
      vz *= v;

      if(!VarsV.isrect) {
        if(!cylinder_intersect(&t0, &t3, x, y, z, vx, vy, vz, radius, yheight))
        {
          /* ??? did not hit cylinder */
          printf("FATAL ERROR: Did not hit cylinder from inside.\n");
          exit(1);
        }
        dt = t3; /* outgoing point */
        if(thickness && cylinder_intersect(&t1, &t2, x, y, z, vx, vy, vz, radius-thickness, yheight) &&
           t2 > 0)
          dt -= (t2-t1);            /* Subtract hollow part */
      }
      else
      {
      if(!box_intersect(&t0, &t3, x, y, z, vx, vy, vz, xwidth, yheight, zdepth))
        {
          /* ??? did not hit box */
          printf("FATAL ERROR: Did not hit box from inside.\n");
          exit(1);
        }
        dt = t3;
      }
      l_o = v*dt; /* trajectory after scattering point: absorption only */

      p *= v/v_i*l_full*VarsV.my_s*exp(-my_a*(l_i+v_i/v*l_o)-VarsV.my_s*l_i);
      /* We do not consider scattering from 2nd part (outgoing) */
      p /= 4*PI/solid_angle;
      p /= p_interact;

      /* Polarisation part (1/3 NSF, 2/3 SF) */
      sx *= -1.0/3.0;
      sy *= -1.0/3.0;
      sz *= -1.0/3.0;
	
      SCATTER;
    }
  else /* Transmitting; always elastic */
    {
      p *= exp(-(my_a+VarsV.my_s)*l_full);
      p /= (1-p_interact);
    }
  }
  #undef thickness
  #undef radius
  #undef focus_r
  #undef p_interact
  #undef f_QE
  #undef f_tun
  #undef gamma
  #undef E_tun
  #undef target_x
  #undef target_y
  #undef target_z
  #undef focus_xw
  #undef focus_yh
  #undef focus_aw
  #undef focus_ah
  #undef xwidth
  #undef yheight
  #undef zdepth
  #undef sigma_abs
  #undef sigma_inc
  #undef Vc
  #undef target_index
  #undef ftun
  #undef fQE
  #undef VarsV
  return(_comp);
} /* class_Tunneling_sample_trace */

#pragma acc routine seq
_class_TOF2E_monitor *class_TOF2E_monitor_trace(_class_TOF2E_monitor *_comp
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
  #define T_zero (_comp->_parameters.T_zero)
  #define L_flight (_comp->_parameters.L_flight)
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
    E = VS2E*(L_flight/(t-T_zero))*(L_flight/(t-T_zero));

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
  #undef T_zero
  #undef L_flight
  #undef restore_neutron
  #undef E_N
  #undef E_p
  #undef E_p2
  #undef S_p
  #undef S_pE
  #undef S_pE2
  return(_comp);
} /* class_TOF2E_monitor_trace */

/* *****************************************************************************
* instrument 'ESS_IN5_reprate' TRACE
***************************************************************************** */

#pragma acc routine seq
int raytrace(_class_particle* _particle) { /* called by mccode_main for ESS_IN5_reprate:TRACE */

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
      /* component source=ESS_moderator_long() [1] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_source->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _source->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_source->_position_relative, _source->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_ESS_moderator_long_trace(_source, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component source [1] */
    if (!ABSORBED && _particle->_index == 2) {
      /* component Origin=Progress_bar() [2] */
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
    } /* end component Origin [2] */
    if (!ABSORBED && _particle->_index == 3) {
      /* component TOFmoderator_zoom=TOF_monitor() [3] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_TOFmoderator_zoom->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _TOFmoderator_zoom->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_TOFmoderator_zoom->_position_relative, _TOFmoderator_zoom->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_TOF_monitor_trace(_TOFmoderator_zoom, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component TOFmoderator_zoom [3] */
    if (!ABSORBED && _particle->_index == 4) {
      /* component TOFmoderator=TOF_monitor() [4] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_TOFmoderator->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _TOFmoderator->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_TOFmoderator->_position_relative, _TOFmoderator->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_TOF_monitor_trace(_TOFmoderator, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component TOFmoderator [4] */
    if (!ABSORBED && _particle->_index == 5) {
      /* component Lmon_guistart=L_monitor() [5] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Lmon_guistart->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Lmon_guistart->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Lmon_guistart->_position_relative, _Lmon_guistart->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_L_monitor_trace(_Lmon_guistart, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component Lmon_guistart [5] */
    if (!ABSORBED && _particle->_index == 6) {
      /* component Lmon_normalize=L_monitor() [6] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Lmon_normalize->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Lmon_normalize->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Lmon_normalize->_position_relative, _Lmon_normalize->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_L_monitor_trace(_Lmon_normalize, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component Lmon_normalize [6] */
    if (!ABSORBED && _particle->_index == 7) {
      /* component Guide1=Guide() [7] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Guide1->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Guide1->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Guide1->_position_relative, _Guide1->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Guide_trace(_Guide1, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component Guide1 [7] */
    if (!ABSORBED && _particle->_index == 8) {
      /* component Lmonslow1=L_monitor() [8] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Lmonslow1->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Lmonslow1->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Lmonslow1->_position_relative, _Lmonslow1->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_L_monitor_trace(_Lmonslow1, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component Lmonslow1 [8] */
    if (!ABSORBED && _particle->_index == 9) {
      /* component PSDslow1=PSD_monitor() [9] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_PSDslow1->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _PSDslow1->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_PSDslow1->_position_relative, _PSDslow1->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_PSD_monitor_trace(_PSDslow1, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component PSDslow1 [9] */
    if (!ABSORBED && _particle->_index == 10) {
      /* component FOchop1=DiskChopper() [10] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_FOchop1->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _FOchop1->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_FOchop1->_position_relative, _FOchop1->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_DiskChopper_trace(_FOchop1, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component FOchop1 [10] */
    if (!ABSORBED && _particle->_index == 11) {
      /* component TOFLmon1=TOFLambda_monitor() [11] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_TOFLmon1->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _TOFLmon1->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_TOFLmon1->_position_relative, _TOFLmon1->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_TOFLambda_monitor_trace(_TOFLmon1, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component TOFLmon1 [11] */
    if (!ABSORBED && _particle->_index == 12) {
      /* component Lmon_afterslow1=L_monitor() [12] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Lmon_afterslow1->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Lmon_afterslow1->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Lmon_afterslow1->_position_relative, _Lmon_afterslow1->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_L_monitor_trace(_Lmon_afterslow1, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component Lmon_afterslow1 [12] */
    if (!ABSORBED && _particle->_index == 13) {
      /* component PSD_afterslow1=PSD_monitor() [13] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_PSD_afterslow1->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _PSD_afterslow1->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_PSD_afterslow1->_position_relative, _PSD_afterslow1->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_PSD_monitor_trace(_PSD_afterslow1, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component PSD_afterslow1 [13] */
    if (!ABSORBED && _particle->_index == 14) {
      /* component Guidelong1=Guide() [14] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Guidelong1->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Guidelong1->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Guidelong1->_position_relative, _Guidelong1->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Guide_trace(_Guidelong1, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component Guidelong1 [14] */
    if (!ABSORBED && _particle->_index == 15) {
      /* component Guidelong1b=Guide() [15] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Guidelong1b->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Guidelong1b->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Guidelong1b->_position_relative, _Guidelong1b->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Guide_trace(_Guidelong1b, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component Guidelong1b [15] */
    if (!ABSORBED && _particle->_index == 16) {
      /* component Lmon_slow2=L_monitor() [16] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Lmon_slow2->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Lmon_slow2->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Lmon_slow2->_position_relative, _Lmon_slow2->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_L_monitor_trace(_Lmon_slow2, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component Lmon_slow2 [16] */
    if (!ABSORBED && _particle->_index == 17) {
      /* component FOchop2=DiskChopper() [17] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_FOchop2->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _FOchop2->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_FOchop2->_position_relative, _FOchop2->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_DiskChopper_trace(_FOchop2, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component FOchop2 [17] */
    if (!ABSORBED && _particle->_index == 18) {
      /* component Fastchop1=DiskChopper() [18] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Fastchop1->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Fastchop1->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Fastchop1->_position_relative, _Fastchop1->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_DiskChopper_trace(_Fastchop1, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component Fastchop1 [18] */
    if (!ABSORBED && _particle->_index == 19) {
      /* component PSD_afterslow2=PSD_monitor() [19] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_PSD_afterslow2->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _PSD_afterslow2->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_PSD_afterslow2->_position_relative, _PSD_afterslow2->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_PSD_monitor_trace(_PSD_afterslow2, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component PSD_afterslow2 [19] */
    if (!ABSORBED && _particle->_index == 20) {
      /* component Lmon_afterslow2=L_monitor() [20] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Lmon_afterslow2->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Lmon_afterslow2->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Lmon_afterslow2->_position_relative, _Lmon_afterslow2->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_L_monitor_trace(_Lmon_afterslow2, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component Lmon_afterslow2 [20] */
    if (!ABSORBED && _particle->_index == 21) {
      /* component TOFL_afterslow2=TOFLambda_monitor() [21] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_TOFL_afterslow2->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _TOFL_afterslow2->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_TOFL_afterslow2->_position_relative, _TOFL_afterslow2->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_TOFLambda_monitor_trace(_TOFL_afterslow2, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component TOFL_afterslow2 [21] */
    if (!ABSORBED && _particle->_index == 22) {
      /* component Guidelong2=Guide() [22] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Guidelong2->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Guidelong2->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Guidelong2->_position_relative, _Guidelong2->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Guide_trace(_Guidelong2, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component Guidelong2 [22] */
    if (!ABSORBED && _particle->_index == 23) {
      /* component Lmon_beforeballistic=L_monitor() [23] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Lmon_beforeballistic->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Lmon_beforeballistic->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Lmon_beforeballistic->_position_relative, _Lmon_beforeballistic->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_L_monitor_trace(_Lmon_beforeballistic, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component Lmon_beforeballistic [23] */
    if (!ABSORBED && _particle->_index == 24) {
      /* component PSD_beforeballistic=PSD_monitor() [24] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_PSD_beforeballistic->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _PSD_beforeballistic->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_PSD_beforeballistic->_position_relative, _PSD_beforeballistic->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_PSD_monitor_trace(_PSD_beforeballistic, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component PSD_beforeballistic [24] */
    if (!ABSORBED && _particle->_index == 25) {
      /* component Guidelong2a=Guide() [25] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Guidelong2a->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Guidelong2a->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Guidelong2a->_position_relative, _Guidelong2a->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Guide_trace(_Guidelong2a, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component Guidelong2a [25] */
    if (!ABSORBED && _particle->_index == 26) {
      /* component Lmonfast2=L_monitor() [26] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Lmonfast2->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Lmonfast2->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Lmonfast2->_position_relative, _Lmonfast2->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_L_monitor_trace(_Lmonfast2, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component Lmonfast2 [26] */
    if (!ABSORBED && _particle->_index == 27) {
      /* component Lmonfast2_zoom=L_monitor() [27] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Lmonfast2_zoom->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Lmonfast2_zoom->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Lmonfast2_zoom->_position_relative, _Lmonfast2_zoom->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_L_monitor_trace(_Lmonfast2_zoom, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component Lmonfast2_zoom [27] */
    if (!ABSORBED && _particle->_index == 28) {
      /* component TOFLfast2=TOFLambda_monitor() [28] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_TOFLfast2->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _TOFLfast2->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_TOFLfast2->_position_relative, _TOFLfast2->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_TOFLambda_monitor_trace(_TOFLfast2, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component TOFLfast2 [28] */
    if (!ABSORBED && _particle->_index == 29) {
      /* component TOFLfast2zoom=TOFLambda_monitor() [29] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_TOFLfast2zoom->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _TOFLfast2zoom->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_TOFLfast2zoom->_position_relative, _TOFLfast2zoom->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_TOFLambda_monitor_trace(_TOFLfast2zoom, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component TOFLfast2zoom [29] */
    if (!ABSORBED && _particle->_index == 30) {
      /* component PSDfast2=PSD_monitor() [30] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_PSDfast2->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _PSDfast2->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_PSDfast2->_position_relative, _PSDfast2->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_PSD_monitor_trace(_PSDfast2, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component PSDfast2 [30] */
    if (!ABSORBED && _particle->_index == 31) {
      /* component Fastchop2=DiskChopper() [31] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Fastchop2->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Fastchop2->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Fastchop2->_position_relative, _Fastchop2->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_DiskChopper_trace(_Fastchop2, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component Fastchop2 [31] */
    if (!ABSORBED && _particle->_index == 32) {
      /* component Fastchop2counter=DiskChopper() [32] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Fastchop2counter->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Fastchop2counter->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Fastchop2counter->_position_relative, _Fastchop2counter->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_DiskChopper_trace(_Fastchop2counter, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component Fastchop2counter [32] */
    if (!ABSORBED && _particle->_index == 33) {
      /* component FOchop3=DiskChopper() [33] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_FOchop3->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _FOchop3->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_FOchop3->_position_relative, _FOchop3->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_DiskChopper_trace(_FOchop3, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component FOchop3 [33] */
    if (!ABSORBED && _particle->_index == 34) {
      /* component TOFfast2_zoom=TOF_monitor() [34] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_TOFfast2_zoom->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _TOFfast2_zoom->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_TOFfast2_zoom->_position_relative, _TOFfast2_zoom->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_TOF_monitor_trace(_TOFfast2_zoom, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component TOFfast2_zoom [34] */
    if (!ABSORBED && _particle->_index == 35) {
      /* component Lmon_afterfast2=L_monitor() [35] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Lmon_afterfast2->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Lmon_afterfast2->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Lmon_afterfast2->_position_relative, _Lmon_afterfast2->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_L_monitor_trace(_Lmon_afterfast2, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component Lmon_afterfast2 [35] */
    if (!ABSORBED && _particle->_index == 36) {
      /* component TOFL_afterfast2=TOFLambda_monitor() [36] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_TOFL_afterfast2->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _TOFL_afterfast2->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_TOFL_afterfast2->_position_relative, _TOFL_afterfast2->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_TOFLambda_monitor_trace(_TOFL_afterfast2, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component TOFL_afterfast2 [36] */
    if (!ABSORBED && _particle->_index == 37) {
      /* component TOFL_afterfast2_zoom=TOFLambda_monitor() [37] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_TOFL_afterfast2_zoom->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _TOFL_afterfast2_zoom->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_TOFL_afterfast2_zoom->_position_relative, _TOFL_afterfast2_zoom->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_TOFLambda_monitor_trace(_TOFL_afterfast2_zoom, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component TOFL_afterfast2_zoom [37] */
    if (!ABSORBED && _particle->_index == 38) {
      /* component PSD_afterfast2=PSD_monitor() [38] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_PSD_afterfast2->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _PSD_afterfast2->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_PSD_afterfast2->_position_relative, _PSD_afterfast2->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_PSD_monitor_trace(_PSD_afterfast2, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component PSD_afterfast2 [38] */
    if (!ABSORBED && _particle->_index == 39) {
      /* component Guidesample=Guide() [39] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Guidesample->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Guidesample->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Guidesample->_position_relative, _Guidesample->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Guide_trace(_Guidesample, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component Guidesample [39] */
    if (!ABSORBED && _particle->_index == 40) {
      /* component Lmon_guideend=L_monitor() [40] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Lmon_guideend->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Lmon_guideend->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Lmon_guideend->_position_relative, _Lmon_guideend->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_L_monitor_trace(_Lmon_guideend, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component Lmon_guideend [40] */
    if (!ABSORBED && _particle->_index == 41) {
      /* component PSDsample=PSD_monitor() [41] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_PSDsample->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _PSDsample->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_PSDsample->_position_relative, _PSDsample->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_PSD_monitor_trace(_PSDsample, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component PSDsample [41] */
    if (!ABSORBED && _particle->_index == 42) {
      /* component TOFsample_zoom=TOF_monitor() [42] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_TOFsample_zoom->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _TOFsample_zoom->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_TOFsample_zoom->_position_relative, _TOFsample_zoom->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_TOF_monitor_trace(_TOFsample_zoom, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component TOFsample_zoom [42] */
    if (!ABSORBED && _particle->_index == 43) {
      /* component Esample=E_monitor() [43] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Esample->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Esample->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Esample->_position_relative, _Esample->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_E_monitor_trace(_Esample, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component Esample [43] */
    if (!ABSORBED && _particle->_index == 44) {
      /* component Lmon_sample_zoom=L_monitor() [44] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Lmon_sample_zoom->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Lmon_sample_zoom->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Lmon_sample_zoom->_position_relative, _Lmon_sample_zoom->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_L_monitor_trace(_Lmon_sample_zoom, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component Lmon_sample_zoom [44] */
    if (!ABSORBED && _particle->_index == 45) {
      /* component sample=Tunneling_sample() [45] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_sample->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _sample->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_sample->_position_relative, _sample->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Tunneling_sample_trace(_sample, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component sample [45] */
    if (!ABSORBED && _particle->_index == 46) {
      /* component detectorarm=Arm() [46] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_detectorarm->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _detectorarm->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_detectorarm->_position_relative, _detectorarm->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      _particle->_index++;
    } /* end component detectorarm [46] */
    if (!ABSORBED && _particle->_index == 47) {
      /* component TOFdetector=TOF_monitor() [47] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_TOFdetector->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _TOFdetector->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_TOFdetector->_position_relative, _TOFdetector->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_TOF_monitor_trace(_TOFdetector, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component TOFdetector [47] */
    if (!ABSORBED && _particle->_index == 48) {
      /* component TOFdetector_zoom=TOF_monitor() [48] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_TOFdetector_zoom->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _TOFdetector_zoom->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_TOFdetector_zoom->_position_relative, _TOFdetector_zoom->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_TOF_monitor_trace(_TOFdetector_zoom, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component TOFdetector_zoom [48] */
    if (!ABSORBED && _particle->_index == 49) {
      /* component Edetector=E_monitor() [49] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Edetector->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Edetector->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Edetector->_position_relative, _Edetector->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_E_monitor_trace(_Edetector, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component Edetector [49] */
    if (!ABSORBED && _particle->_index == 50) {
      /* component TOF2Edetector=TOF2E_monitor() [50] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_TOF2Edetector->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _TOF2Edetector->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_TOF2Edetector->_position_relative, _TOF2Edetector->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_TOF2E_monitor_trace(_TOF2Edetector, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component TOF2Edetector [50] */
    if (_particle->_index > 50)
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
* instrument 'ESS_IN5_reprate' and components SAVE
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

_class_TOFLambda_monitor *class_TOFLambda_monitor_save(_class_TOFLambda_monitor *_comp
) {
  #define nL (_comp->_parameters.nL)
  #define nt (_comp->_parameters.nt)
  #define tmin (_comp->_parameters.tmin)
  #define tmax (_comp->_parameters.tmax)
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
  #define tt_0 (_comp->_parameters.tt_0)
  #define tt_1 (_comp->_parameters.tt_1)
  #define TOFL_N (_comp->_parameters.TOFL_N)
  #define TOFL_p (_comp->_parameters.TOFL_p)
  #define TOFL_p2 (_comp->_parameters.TOFL_p2)
  DETECTOR_OUT_2D(
    "TOF-wavelength monitor",
    "Time-of-flight [\\gms]", "Wavelength [AA]",
    tmin, tmax, Lmin, Lmax,
    nt, nL,
    &TOFL_N[0][0],&TOFL_p[0][0],&TOFL_p2[0][0],
    filename);
  #undef nL
  #undef nt
  #undef tmin
  #undef tmax
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
  #undef tt_0
  #undef tt_1
  #undef TOFL_N
  #undef TOFL_p
  #undef TOFL_p2
  return(_comp);
} /* class_TOFLambda_monitor_save */

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

_class_TOF2E_monitor *class_TOF2E_monitor_save(_class_TOF2E_monitor *_comp
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
  #define T_zero (_comp->_parameters.T_zero)
  #define L_flight (_comp->_parameters.L_flight)
  #define restore_neutron (_comp->_parameters.restore_neutron)
  #define E_N (_comp->_parameters.E_N)
  #define E_p (_comp->_parameters.E_p)
  #define E_p2 (_comp->_parameters.E_p2)
  #define S_p (_comp->_parameters.S_p)
  #define S_pE (_comp->_parameters.S_pE)
  #define S_pE2 (_comp->_parameters.S_pE2)
  DETECTOR_OUT_1D(
      "TOF-to-Energy monitor",
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
  #undef T_zero
  #undef L_flight
  #undef restore_neutron
  #undef E_N
  #undef E_p
  #undef E_p2
  #undef S_p
  #undef S_pE
  #undef S_pE2
  return(_comp);
} /* class_TOF2E_monitor_save */



int save(FILE *handle) { /* called by mccode_main for ESS_IN5_reprate:SAVE */
  if (!handle) siminfo_init(NULL);

  /* call iteratively all components SAVE */

  class_Progress_bar_save(_Origin);

  class_TOF_monitor_save(_TOFmoderator_zoom);

  class_TOF_monitor_save(_TOFmoderator);

  class_L_monitor_save(_Lmon_guistart);

  class_L_monitor_save(_Lmon_normalize);


  class_L_monitor_save(_Lmonslow1);

  class_PSD_monitor_save(_PSDslow1);


  class_TOFLambda_monitor_save(_TOFLmon1);

  class_L_monitor_save(_Lmon_afterslow1);

  class_PSD_monitor_save(_PSD_afterslow1);



  class_L_monitor_save(_Lmon_slow2);



  class_PSD_monitor_save(_PSD_afterslow2);

  class_L_monitor_save(_Lmon_afterslow2);

  class_TOFLambda_monitor_save(_TOFL_afterslow2);


  class_L_monitor_save(_Lmon_beforeballistic);

  class_PSD_monitor_save(_PSD_beforeballistic);


  class_L_monitor_save(_Lmonfast2);

  class_L_monitor_save(_Lmonfast2_zoom);

  class_TOFLambda_monitor_save(_TOFLfast2);

  class_TOFLambda_monitor_save(_TOFLfast2zoom);

  class_PSD_monitor_save(_PSDfast2);




  class_TOF_monitor_save(_TOFfast2_zoom);

  class_L_monitor_save(_Lmon_afterfast2);

  class_TOFLambda_monitor_save(_TOFL_afterfast2);

  class_TOFLambda_monitor_save(_TOFL_afterfast2_zoom);

  class_PSD_monitor_save(_PSD_afterfast2);


  class_L_monitor_save(_Lmon_guideend);

  class_PSD_monitor_save(_PSDsample);

  class_TOF_monitor_save(_TOFsample_zoom);

  class_E_monitor_save(_Esample);

  class_L_monitor_save(_Lmon_sample_zoom);



  class_TOF_monitor_save(_TOFdetector);

  class_TOF_monitor_save(_TOFdetector_zoom);

  class_E_monitor_save(_Edetector);

  class_TOF2E_monitor_save(_TOF2Edetector);

  if (!handle) siminfo_close(); 

  return(0);
} /* save */

/* *****************************************************************************
* instrument 'ESS_IN5_reprate' and components FINALLY
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

_class_TOFLambda_monitor *class_TOFLambda_monitor_finally(_class_TOFLambda_monitor *_comp
) {
  #define nL (_comp->_parameters.nL)
  #define nt (_comp->_parameters.nt)
  #define tmin (_comp->_parameters.tmin)
  #define tmax (_comp->_parameters.tmax)
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
  #define tt_0 (_comp->_parameters.tt_0)
  #define tt_1 (_comp->_parameters.tt_1)
  #define TOFL_N (_comp->_parameters.TOFL_N)
  #define TOFL_p (_comp->_parameters.TOFL_p)
  #define TOFL_p2 (_comp->_parameters.TOFL_p2)
  destroy_darr2d(TOFL_N);
  destroy_darr2d(TOFL_p);
  destroy_darr2d(TOFL_p2);
  #undef nL
  #undef nt
  #undef tmin
  #undef tmax
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
  #undef tt_0
  #undef tt_1
  #undef TOFL_N
  #undef TOFL_p
  #undef TOFL_p2
  return(_comp);
} /* class_TOFLambda_monitor_finally */

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

_class_TOF2E_monitor *class_TOF2E_monitor_finally(_class_TOF2E_monitor *_comp
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
  #define T_zero (_comp->_parameters.T_zero)
  #define L_flight (_comp->_parameters.L_flight)
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
  #undef T_zero
  #undef L_flight
  #undef restore_neutron
  #undef E_N
  #undef E_p
  #undef E_p2
  #undef S_p
  #undef S_pE
  #undef S_pE2
  return(_comp);
} /* class_TOF2E_monitor_finally */



int finally(void) { /* called by mccode_main for ESS_IN5_reprate:FINALLY */
  siminfo_init(NULL);
  save(siminfo_file); /* save data when simulation ends */

  /* call iteratively all components FINALLY */

  class_Progress_bar_finally(_Origin);

  class_TOF_monitor_finally(_TOFmoderator_zoom);

  class_TOF_monitor_finally(_TOFmoderator);

  class_L_monitor_finally(_Lmon_guistart);

  class_L_monitor_finally(_Lmon_normalize);


  class_L_monitor_finally(_Lmonslow1);

  class_PSD_monitor_finally(_PSDslow1);


  class_TOFLambda_monitor_finally(_TOFLmon1);

  class_L_monitor_finally(_Lmon_afterslow1);

  class_PSD_monitor_finally(_PSD_afterslow1);



  class_L_monitor_finally(_Lmon_slow2);



  class_PSD_monitor_finally(_PSD_afterslow2);

  class_L_monitor_finally(_Lmon_afterslow2);

  class_TOFLambda_monitor_finally(_TOFL_afterslow2);


  class_L_monitor_finally(_Lmon_beforeballistic);

  class_PSD_monitor_finally(_PSD_beforeballistic);


  class_L_monitor_finally(_Lmonfast2);

  class_L_monitor_finally(_Lmonfast2_zoom);

  class_TOFLambda_monitor_finally(_TOFLfast2);

  class_TOFLambda_monitor_finally(_TOFLfast2zoom);

  class_PSD_monitor_finally(_PSDfast2);




  class_TOF_monitor_finally(_TOFfast2_zoom);

  class_L_monitor_finally(_Lmon_afterfast2);

  class_TOFLambda_monitor_finally(_TOFL_afterfast2);

  class_TOFLambda_monitor_finally(_TOFL_afterfast2_zoom);

  class_PSD_monitor_finally(_PSD_afterfast2);


  class_L_monitor_finally(_Lmon_guideend);

  class_PSD_monitor_finally(_PSDsample);

  class_TOF_monitor_finally(_TOFsample_zoom);

  class_E_monitor_finally(_Esample);

  class_L_monitor_finally(_Lmon_sample_zoom);



  class_TOF_monitor_finally(_TOFdetector);

  class_TOF_monitor_finally(_TOFdetector_zoom);

  class_E_monitor_finally(_Edetector);

  class_TOF2E_monitor_finally(_TOF2Edetector);

  siminfo_close(); 

  return(0);
} /* finally */

/* *****************************************************************************
* instrument 'ESS_IN5_reprate' and components DISPLAY
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
_class_ESS_moderator_long *class_ESS_moderator_long_display(_class_ESS_moderator_long *_comp
) {
  #define width_c (_comp->_parameters.width_c)
  #define yheight (_comp->_parameters.yheight)
  #define Lmin (_comp->_parameters.Lmin)
  #define Lmax (_comp->_parameters.Lmax)
  #define dist (_comp->_parameters.dist)
  #define focus_xw (_comp->_parameters.focus_xw)
  #define focus_yh (_comp->_parameters.focus_yh)
  #define nu (_comp->_parameters.nu)
  #define T (_comp->_parameters.T)
  #define tau (_comp->_parameters.tau)
  #define tau1 (_comp->_parameters.tau1)
  #define tau2 (_comp->_parameters.tau2)
  #define d (_comp->_parameters.d)
  #define n (_comp->_parameters.n)
  #define cold_frac (_comp->_parameters.cold_frac)
  #define n2 (_comp->_parameters.n2)
  #define chi2 (_comp->_parameters.chi2)
  #define I0 (_comp->_parameters.I0)
  #define I2 (_comp->_parameters.I2)
  #define target_index (_comp->_parameters.target_index)
  #define cyl_radius (_comp->_parameters.cyl_radius)
  #define branch1 (_comp->_parameters.branch1)
  #define branch2 (_comp->_parameters.branch2)
  #define branch_tail (_comp->_parameters.branch_tail)
  #define n_pulses (_comp->_parameters.n_pulses)
  #define width_t (_comp->_parameters.width_t)
  #define T_t (_comp->_parameters.T_t)
  #define tau_t (_comp->_parameters.tau_t)
  #define tau1_t (_comp->_parameters.tau1_t)
  #define tau2_t (_comp->_parameters.tau2_t)
  #define chi2_t (_comp->_parameters.chi2_t)
  #define I0_t (_comp->_parameters.I0_t)
  #define I2_t (_comp->_parameters.I2_t)
  #define branch1_t (_comp->_parameters.branch1_t)
  #define branch2_t (_comp->_parameters.branch2_t)
  #define src_2012 (_comp->_parameters.src_2012)
  #define tfocus_dist (_comp->_parameters.tfocus_dist)
  #define tfocus_time (_comp->_parameters.tfocus_time)
  #define tfocus_width (_comp->_parameters.tfocus_width)
  #define beamport_angle (_comp->_parameters.beamport_angle)
  #define l_range (_comp->_parameters.l_range)
  #define w_mult (_comp->_parameters.w_mult)
  #define w_geom (_comp->_parameters.w_geom)
  #define w_geom_c (_comp->_parameters.w_geom_c)
  #define w_geom_t (_comp->_parameters.w_geom_t)
  #define tx (_comp->_parameters.tx)
  #define ty (_comp->_parameters.ty)
  #define tz (_comp->_parameters.tz)
  #define t1x (_comp->_parameters.t1x)
  #define t1y (_comp->_parameters.t1y)
  #define t1z (_comp->_parameters.t1z)
  #define t2x (_comp->_parameters.t2x)
  #define t2y (_comp->_parameters.t2y)
  #define t2z (_comp->_parameters.t2z)
  #define T_n (_comp->_parameters.T_n)
  #define tau_n (_comp->_parameters.tau_n)
  #define tau1_n (_comp->_parameters.tau1_n)
  #define tau2_n (_comp->_parameters.tau2_n)
  #define chi2_n (_comp->_parameters.chi2_n)
  #define I0_n (_comp->_parameters.I0_n)
  #define I2_n (_comp->_parameters.I2_n)
  #define branch1_n (_comp->_parameters.branch1_n)
  #define branch2_n (_comp->_parameters.branch2_n)
  #define r_empty (_comp->_parameters.r_empty)
  #define r_optics (_comp->_parameters.r_optics)
  /* Draw cold moderator as cylinder */
  
  circle("xz", 0,  yheight/2.0, 0, cyl_radius);
  circle("xz", 0,  -yheight/2.0, 0, cyl_radius);
  line(0, -yheight/2.0, cyl_radius, 0, yheight/2.0, cyl_radius);
  line(0, -yheight/2.0, -cyl_radius, 0, yheight/2.0, -cyl_radius);
  line(cyl_radius, yheight/2.0, 0, cyl_radius, yheight/2.0, 0);
  line(-cyl_radius, -yheight/2.0, 0, -cyl_radius, yheight/2.0, 0);
  /* Draw thermal moderators as a couple of squares + some lines */
  // Left
  multiline(4, t1x*cyl_radius, -yheight/2.0, t1z*cyl_radius,
	    t1x*(cyl_radius + width_t), -yheight/2.0, t1z*(cyl_radius + width_t),
      	    t1x*(cyl_radius + width_t), yheight/2.0,  t1z*(cyl_radius + width_t),
	    t1x*cyl_radius, yheight/2.0, t1z*cyl_radius);
	    // Right
  multiline(4, t2x*cyl_radius, -yheight/2.0, t2z*cyl_radius,
	    t2x*(cyl_radius + width_t), -yheight/2.0, t2z*(cyl_radius + width_t),
      	    t2x*(cyl_radius + width_t), yheight/2.0,  t2z*(cyl_radius + width_t),
	    t2x*cyl_radius, yheight/2.0, t2z*cyl_radius);

  /* Dashed lines for indicating "beam extraction" area... */
  dashed_line(t1x*cyl_radius, -yheight/2.0, t1z*cyl_radius, t1x*r_empty, -yheight/2.0, t1z*r_empty,10);
  dashed_line(t1x*cyl_radius, yheight/2.0, t1z*cyl_radius, t1x*r_empty, yheight/2.0, t1z*r_empty,10);
  dashed_line(t2x*cyl_radius, -yheight/2.0, t2z*cyl_radius, t2x*r_empty, -yheight/2.0, t2z*r_empty,5);
  dashed_line(t2x*cyl_radius, yheight/2.0, t2z*cyl_radius, t2x*r_empty, yheight/2.0, t2z*r_empty,5);

  /* Circles indicating extent of the "empty" zone where optics is not allowed */
  circle("xz", 0,  yheight/2.0, 0, r_empty);
  circle("xz", 0,  -yheight/2.0, 0, r_empty);

  /* Circles indicating the builk shielding of the target monolith at 6 m */
  circle("xz", 0,  focus_yh/2.0 , 0, 6);
  circle("xz", 0, -focus_yh/2.0 , 0, 6);
  circle("xz", 0,  2, 0, 6);
  circle("xz", 0, -2, 0, 6);

  /* Rectangle indicating the chosen focus rectangle - where the optics starts... */
  rectangle("xy",tx,ty,tz,focus_xw,focus_yh);
  #undef width_c
  #undef yheight
  #undef Lmin
  #undef Lmax
  #undef dist
  #undef focus_xw
  #undef focus_yh
  #undef nu
  #undef T
  #undef tau
  #undef tau1
  #undef tau2
  #undef d
  #undef n
  #undef cold_frac
  #undef n2
  #undef chi2
  #undef I0
  #undef I2
  #undef target_index
  #undef cyl_radius
  #undef branch1
  #undef branch2
  #undef branch_tail
  #undef n_pulses
  #undef width_t
  #undef T_t
  #undef tau_t
  #undef tau1_t
  #undef tau2_t
  #undef chi2_t
  #undef I0_t
  #undef I2_t
  #undef branch1_t
  #undef branch2_t
  #undef src_2012
  #undef tfocus_dist
  #undef tfocus_time
  #undef tfocus_width
  #undef beamport_angle
  #undef l_range
  #undef w_mult
  #undef w_geom
  #undef w_geom_c
  #undef w_geom_t
  #undef tx
  #undef ty
  #undef tz
  #undef t1x
  #undef t1y
  #undef t1z
  #undef t2x
  #undef t2y
  #undef t2z
  #undef T_n
  #undef tau_n
  #undef tau1_n
  #undef tau2_n
  #undef chi2_n
  #undef I0_n
  #undef I2_n
  #undef branch1_n
  #undef branch2_n
  #undef r_empty
  #undef r_optics
  return(_comp);
} /* class_ESS_moderator_long_display */

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

_class_Guide *class_Guide_display(_class_Guide *_comp
) {
  #define reflect (_comp->_parameters.reflect)
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
  #define pTable (_comp->_parameters.pTable)
  
  multiline(5,
            -w1/2.0, -h1/2.0, 0.0,
             w1/2.0, -h1/2.0, 0.0,
             w1/2.0,  h1/2.0, 0.0,
            -w1/2.0,  h1/2.0, 0.0,
            -w1/2.0, -h1/2.0, 0.0);
  multiline(5,
            -w2/2.0, -h2/2.0, (double)l,
             w2/2.0, -h2/2.0, (double)l,
             w2/2.0,  h2/2.0, (double)l,
            -w2/2.0,  h2/2.0, (double)l,
            -w2/2.0, -h2/2.0, (double)l);
  line(-w1/2.0, -h1/2.0, 0, -w2/2.0, -h2/2.0, (double)l);
  line( w1/2.0, -h1/2.0, 0,  w2/2.0, -h2/2.0, (double)l);
  line( w1/2.0,  h1/2.0, 0,  w2/2.0,  h2/2.0, (double)l);
  line(-w1/2.0,  h1/2.0, 0, -w2/2.0,  h2/2.0, (double)l);
  #undef reflect
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
  #undef pTable
  return(_comp);
} /* class_Guide_display */

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

_class_TOFLambda_monitor *class_TOFLambda_monitor_display(_class_TOFLambda_monitor *_comp
) {
  #define nL (_comp->_parameters.nL)
  #define nt (_comp->_parameters.nt)
  #define tmin (_comp->_parameters.tmin)
  #define tmax (_comp->_parameters.tmax)
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
  #define tt_0 (_comp->_parameters.tt_0)
  #define tt_1 (_comp->_parameters.tt_1)
  #define TOFL_N (_comp->_parameters.TOFL_N)
  #define TOFL_p (_comp->_parameters.TOFL_p)
  #define TOFL_p2 (_comp->_parameters.TOFL_p2)
  multiline(5, (double)xmin, (double)ymin, 0.0,
               (double)xmax, (double)ymin, 0.0,
               (double)xmax, (double)ymax, 0.0,
               (double)xmin, (double)ymax, 0.0,
               (double)xmin, (double)ymin, 0.0);
  #undef nL
  #undef nt
  #undef tmin
  #undef tmax
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
  #undef tt_0
  #undef tt_1
  #undef TOFL_N
  #undef TOFL_p
  #undef TOFL_p2
  return(_comp);
} /* class_TOFLambda_monitor_display */

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

_class_Tunneling_sample *class_Tunneling_sample_display(_class_Tunneling_sample *_comp
) {
  #define thickness (_comp->_parameters.thickness)
  #define radius (_comp->_parameters.radius)
  #define focus_r (_comp->_parameters.focus_r)
  #define p_interact (_comp->_parameters.p_interact)
  #define f_QE (_comp->_parameters.f_QE)
  #define f_tun (_comp->_parameters.f_tun)
  #define gamma (_comp->_parameters.gamma)
  #define E_tun (_comp->_parameters.E_tun)
  #define target_x (_comp->_parameters.target_x)
  #define target_y (_comp->_parameters.target_y)
  #define target_z (_comp->_parameters.target_z)
  #define focus_xw (_comp->_parameters.focus_xw)
  #define focus_yh (_comp->_parameters.focus_yh)
  #define focus_aw (_comp->_parameters.focus_aw)
  #define focus_ah (_comp->_parameters.focus_ah)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define zdepth (_comp->_parameters.zdepth)
  #define sigma_abs (_comp->_parameters.sigma_abs)
  #define sigma_inc (_comp->_parameters.sigma_inc)
  #define Vc (_comp->_parameters.Vc)
  #define target_index (_comp->_parameters.target_index)
  #define ftun (_comp->_parameters.ftun)
  #define fQE (_comp->_parameters.fQE)
  #define VarsV (_comp->_parameters.VarsV)
  
  if (!VarsV.isrect) 
  {
    circle("xz", 0,  yheight/2.0, 0, radius);
    circle("xz", 0, -yheight/2.0, 0, radius);
    line(-radius, -yheight/2.0, 0, -radius, +yheight/2.0, 0);
    line(+radius, -yheight/2.0, 0, +radius, +yheight/2.0, 0);
    line(0, -yheight/2.0, -radius, 0, +yheight/2.0, -radius);
    line(0, -yheight/2.0, +radius, 0, +yheight/2.0, +radius);
    if (thickness) 
    {
      double radius_i=radius-thickness;
      circle("xz", 0,  yheight/2.0, 0, radius_i);
      circle("xz", 0, -yheight/2.0, 0, radius_i);
      line(-radius_i, -yheight/2.0, 0, -radius_i, +yheight/2.0, 0);
      line(+radius_i, -yheight/2.0, 0, +radius_i, +yheight/2.0, 0);
      line(0, -yheight/2.0, -radius_i, 0, +yheight/2.0, -radius_i);
      line(0, -yheight/2.0, +radius_i, 0, +yheight/2.0, +radius_i);
    }
  }
  else
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
  }
  #undef thickness
  #undef radius
  #undef focus_r
  #undef p_interact
  #undef f_QE
  #undef f_tun
  #undef gamma
  #undef E_tun
  #undef target_x
  #undef target_y
  #undef target_z
  #undef focus_xw
  #undef focus_yh
  #undef focus_aw
  #undef focus_ah
  #undef xwidth
  #undef yheight
  #undef zdepth
  #undef sigma_abs
  #undef sigma_inc
  #undef Vc
  #undef target_index
  #undef ftun
  #undef fQE
  #undef VarsV
  return(_comp);
} /* class_Tunneling_sample_display */

_class_Arm *class_Arm_display(_class_Arm *_comp
) {
  /* A bit ugly; hard-coded dimensions. */
  
  line(0,0,0,0.2,0,0);
  line(0,0,0,0,0.2,0);
  line(0,0,0,0,0,0.2);
  return(_comp);
} /* class_Arm_display */

_class_TOF2E_monitor *class_TOF2E_monitor_display(_class_TOF2E_monitor *_comp
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
  #define T_zero (_comp->_parameters.T_zero)
  #define L_flight (_comp->_parameters.L_flight)
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
  #undef T_zero
  #undef L_flight
  #undef restore_neutron
  #undef E_N
  #undef E_p
  #undef E_p2
  #undef S_p
  #undef S_pE
  #undef S_pE2
  return(_comp);
} /* class_TOF2E_monitor_display */


  #undef magnify
  #undef line
  #undef dashed_line
  #undef multiline
  #undef rectangle
  #undef box
  #undef circle
  #undef cylinder
  #undef sphere

int display(void) { /* called by mccode_main for ESS_IN5_reprate:DISPLAY */
  printf("MCDISPLAY: start\n");

  /* call iteratively all components DISPLAY */
  class_ESS_moderator_long_display(_source);

  class_Progress_bar_display(_Origin);

  class_TOF_monitor_display(_TOFmoderator_zoom);

  class_TOF_monitor_display(_TOFmoderator);

  class_L_monitor_display(_Lmon_guistart);

  class_L_monitor_display(_Lmon_normalize);

  class_Guide_display(_Guide1);

  class_L_monitor_display(_Lmonslow1);

  class_PSD_monitor_display(_PSDslow1);

  class_DiskChopper_display(_FOchop1);

  class_TOFLambda_monitor_display(_TOFLmon1);

  class_L_monitor_display(_Lmon_afterslow1);

  class_PSD_monitor_display(_PSD_afterslow1);

  class_Guide_display(_Guidelong1);

  class_Guide_display(_Guidelong1b);

  class_L_monitor_display(_Lmon_slow2);

  class_DiskChopper_display(_FOchop2);

  class_DiskChopper_display(_Fastchop1);

  class_PSD_monitor_display(_PSD_afterslow2);

  class_L_monitor_display(_Lmon_afterslow2);

  class_TOFLambda_monitor_display(_TOFL_afterslow2);

  class_Guide_display(_Guidelong2);

  class_L_monitor_display(_Lmon_beforeballistic);

  class_PSD_monitor_display(_PSD_beforeballistic);

  class_Guide_display(_Guidelong2a);

  class_L_monitor_display(_Lmonfast2);

  class_L_monitor_display(_Lmonfast2_zoom);

  class_TOFLambda_monitor_display(_TOFLfast2);

  class_TOFLambda_monitor_display(_TOFLfast2zoom);

  class_PSD_monitor_display(_PSDfast2);

  class_DiskChopper_display(_Fastchop2);

  class_DiskChopper_display(_Fastchop2counter);

  class_DiskChopper_display(_FOchop3);

  class_TOF_monitor_display(_TOFfast2_zoom);

  class_L_monitor_display(_Lmon_afterfast2);

  class_TOFLambda_monitor_display(_TOFL_afterfast2);

  class_TOFLambda_monitor_display(_TOFL_afterfast2_zoom);

  class_PSD_monitor_display(_PSD_afterfast2);

  class_Guide_display(_Guidesample);

  class_L_monitor_display(_Lmon_guideend);

  class_PSD_monitor_display(_PSDsample);

  class_TOF_monitor_display(_TOFsample_zoom);

  class_E_monitor_display(_Esample);

  class_L_monitor_display(_Lmon_sample_zoom);

  class_Tunneling_sample_display(_sample);

  class_Arm_display(_detectorarm);

  class_TOF_monitor_display(_TOFdetector);

  class_TOF_monitor_display(_TOFdetector_zoom);

  class_E_monitor_display(_Edetector);

  class_TOF2E_monitor_display(_TOF2Edetector);

  printf("MCDISPLAY: end\n");

  return(0);
} /* display */

void* _getvar_parameters(char* compname)
/* enables settings parameters based use of the GETPAR macro */
{
  if (!strcmp(compname, "source")) return (void *) &(_source->_parameters);
  if (!strcmp(compname, "Origin")) return (void *) &(_Origin->_parameters);
  if (!strcmp(compname, "TOFmoderator_zoom")) return (void *) &(_TOFmoderator_zoom->_parameters);
  if (!strcmp(compname, "TOFmoderator")) return (void *) &(_TOFmoderator->_parameters);
  if (!strcmp(compname, "Lmon_guistart")) return (void *) &(_Lmon_guistart->_parameters);
  if (!strcmp(compname, "Lmon_normalize")) return (void *) &(_Lmon_normalize->_parameters);
  if (!strcmp(compname, "Guide1")) return (void *) &(_Guide1->_parameters);
  if (!strcmp(compname, "Lmonslow1")) return (void *) &(_Lmonslow1->_parameters);
  if (!strcmp(compname, "PSDslow1")) return (void *) &(_PSDslow1->_parameters);
  if (!strcmp(compname, "FOchop1")) return (void *) &(_FOchop1->_parameters);
  if (!strcmp(compname, "TOFLmon1")) return (void *) &(_TOFLmon1->_parameters);
  if (!strcmp(compname, "Lmon_afterslow1")) return (void *) &(_Lmon_afterslow1->_parameters);
  if (!strcmp(compname, "PSD_afterslow1")) return (void *) &(_PSD_afterslow1->_parameters);
  if (!strcmp(compname, "Guidelong1")) return (void *) &(_Guidelong1->_parameters);
  if (!strcmp(compname, "Guidelong1b")) return (void *) &(_Guidelong1b->_parameters);
  if (!strcmp(compname, "Lmon_slow2")) return (void *) &(_Lmon_slow2->_parameters);
  if (!strcmp(compname, "FOchop2")) return (void *) &(_FOchop2->_parameters);
  if (!strcmp(compname, "Fastchop1")) return (void *) &(_Fastchop1->_parameters);
  if (!strcmp(compname, "PSD_afterslow2")) return (void *) &(_PSD_afterslow2->_parameters);
  if (!strcmp(compname, "Lmon_afterslow2")) return (void *) &(_Lmon_afterslow2->_parameters);
  if (!strcmp(compname, "TOFL_afterslow2")) return (void *) &(_TOFL_afterslow2->_parameters);
  if (!strcmp(compname, "Guidelong2")) return (void *) &(_Guidelong2->_parameters);
  if (!strcmp(compname, "Lmon_beforeballistic")) return (void *) &(_Lmon_beforeballistic->_parameters);
  if (!strcmp(compname, "PSD_beforeballistic")) return (void *) &(_PSD_beforeballistic->_parameters);
  if (!strcmp(compname, "Guidelong2a")) return (void *) &(_Guidelong2a->_parameters);
  if (!strcmp(compname, "Lmonfast2")) return (void *) &(_Lmonfast2->_parameters);
  if (!strcmp(compname, "Lmonfast2_zoom")) return (void *) &(_Lmonfast2_zoom->_parameters);
  if (!strcmp(compname, "TOFLfast2")) return (void *) &(_TOFLfast2->_parameters);
  if (!strcmp(compname, "TOFLfast2zoom")) return (void *) &(_TOFLfast2zoom->_parameters);
  if (!strcmp(compname, "PSDfast2")) return (void *) &(_PSDfast2->_parameters);
  if (!strcmp(compname, "Fastchop2")) return (void *) &(_Fastchop2->_parameters);
  if (!strcmp(compname, "Fastchop2counter")) return (void *) &(_Fastchop2counter->_parameters);
  if (!strcmp(compname, "FOchop3")) return (void *) &(_FOchop3->_parameters);
  if (!strcmp(compname, "TOFfast2_zoom")) return (void *) &(_TOFfast2_zoom->_parameters);
  if (!strcmp(compname, "Lmon_afterfast2")) return (void *) &(_Lmon_afterfast2->_parameters);
  if (!strcmp(compname, "TOFL_afterfast2")) return (void *) &(_TOFL_afterfast2->_parameters);
  if (!strcmp(compname, "TOFL_afterfast2_zoom")) return (void *) &(_TOFL_afterfast2_zoom->_parameters);
  if (!strcmp(compname, "PSD_afterfast2")) return (void *) &(_PSD_afterfast2->_parameters);
  if (!strcmp(compname, "Guidesample")) return (void *) &(_Guidesample->_parameters);
  if (!strcmp(compname, "Lmon_guideend")) return (void *) &(_Lmon_guideend->_parameters);
  if (!strcmp(compname, "PSDsample")) return (void *) &(_PSDsample->_parameters);
  if (!strcmp(compname, "TOFsample_zoom")) return (void *) &(_TOFsample_zoom->_parameters);
  if (!strcmp(compname, "Esample")) return (void *) &(_Esample->_parameters);
  if (!strcmp(compname, "Lmon_sample_zoom")) return (void *) &(_Lmon_sample_zoom->_parameters);
  if (!strcmp(compname, "sample")) return (void *) &(_sample->_parameters);
  if (!strcmp(compname, "detectorarm")) return (void *) &(_detectorarm->_parameters);
  if (!strcmp(compname, "TOFdetector")) return (void *) &(_TOFdetector->_parameters);
  if (!strcmp(compname, "TOFdetector_zoom")) return (void *) &(_TOFdetector_zoom->_parameters);
  if (!strcmp(compname, "Edetector")) return (void *) &(_Edetector->_parameters);
  if (!strcmp(compname, "TOF2Edetector")) return (void *) &(_TOF2Edetector->_parameters);
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

/* end of generated C code ESS_IN5_reprate.c */
