/* Automatically generated file. Do not edit. 
 * Format:     ANSI C source code
 * Creator:    McStas <http://www.mcstas.org>
 * Instrument: ISIS_OSIRIS.instr (ISIS_OSIRIS)
 * Date:       Sat Oct 12 09:07:05 2019
 * File:       ISIS_OSIRIS.c
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
* Start of instrument 'ISIS_OSIRIS' generated code
***************************************************************************** */

#ifdef MC_TRACE_ENABLED
int traceenabled = 1;
#else
int traceenabled = 0;
#endif
#define MCSTAS "/usr/share/mcstas/3.0-dev/"
int   defaultmain         = 1;
char  instrument_name[]   = "ISIS_OSIRIS";
char  instrument_source[] = "ISIS_OSIRIS.instr";
char *instrument_exe      = NULL; /* will be set to argv[0] in main */
char  instrument_code[]   = "Instrument ISIS_OSIRIS source code ISIS_OSIRIS.instr is not embedded in this executable.\n  Use --source option when running McStas.\n";

int main(int argc, char *argv[]){return mccode_main(argc, argv);}

/* *****************************************************************************
* instrument 'ISIS_OSIRIS' and components DECLARE
***************************************************************************** */

/* Instrument parameters: structure and a table for the initialisation
   (Used in e.g. inputparse and I/O function (e.g. detector_out) */

struct _struct_instrument_parameters {
  MCNUM _LAMBDA;
  MCNUM _DLAMBDA;
  MCNUM _GUIDEREFLECTIVITY;
};
typedef struct _struct_instrument_parameters _class_instrument_parameters;

/* instrument SPLIT, GROUP and JUMP control logic */
struct instrument_logic_struct {
  long Group_graphiteAnalyzer; /* equals index of scattering comp when in group */
};

struct _instrument_struct {
  char   _name[256]; /* the name of this instrument e.g. 'ISIS_OSIRIS' */
/* Counters per component instance */
  double counter_AbsorbProp[74]; /* absorbed events in PROP routines */
  double counter_N[74], counter_P[74], counter_P2[74]; /* event counters after each component instance */
  _class_particle _trajectory[74]; /* current trajectory for STORE/RESTORE */
/* Components position table (absolute and relative coords) */
  Coords _position_relative[74]; /* positions of all components */
  Coords _position_absolute[74];
  _class_instrument_parameters _parameters; /* instrument parameters */
  struct instrument_logic_struct logic; /* instrument logic */
} _instrument_var;
struct _instrument_struct *instrument = & _instrument_var;
#pragma acc declare create ( _instrument_var )
#pragma acc declare create ( instrument )

int numipar = 3;
struct mcinputtable_struct mcinputtable[] = {
  "LAMBDA", &(_instrument_var._parameters._LAMBDA), instr_type_double, "6.66", 
  "DLAMBDA", &(_instrument_var._parameters._DLAMBDA), instr_type_double, "0.1", 
  "GUIDEREFLECTIVITY", &(_instrument_var._parameters._GUIDEREFLECTIVITY), instr_type_double, "1.0", 
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


/* Shared user declarations for all components types 'Guide_curved'. */


/* Shared user declarations for all components types 'Incoherent'. */

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


struct StructVarsInc
{
double  sigma_a; /* Absorption cross section per atom (barns) */
    double  sigma_i; /* Incoherent scattering cross section per atom (barns) */
    double  rho;     /* Density of atoms (AA-3) */
    double  my_s;
    double  my_a_v;
    int     shape;       /* 0 cylinder, 1 box, 2 sphere, 3 OFF file */
    double  aw,ah;       /* rectangular angular dimensions */
    double  xw,yh;       /* rectangular metrical dimensions */
    double  tx,ty,tz;    /* target coords */
    };

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



/* ************************************************************************** */
/*             End of SHARE user declarations for all components              */
/* ************************************************************************** */


/* ********************** component definition declarations. **************** */

/* component armStart=Arm() [1] DECLARE */
/* Parameter definition for component type 'Arm' */
struct _struct_Arm_parameters {
  char Arm_has_no_parameters;
}; /* _struct_Arm_parameters */
typedef struct _struct_Arm_parameters _class_Arm_parameters;

/* Parameters for component type 'Arm' */
struct _struct_Arm {
  char     _name[256]; /* e.g. armStart */
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
_class_Arm _armStart_var;
_class_Arm *_armStart = &_armStart_var;
#pragma acc declare create ( _armStart_var )
#pragma acc declare create ( _armStart )

/* component isis_mod=ISIS_moderator() [2] DECLARE */
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
  char     _name[256]; /* e.g. isis_mod */
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
_class_ISIS_moderator _isis_mod_var;
_class_ISIS_moderator *_isis_mod = &_isis_mod_var;
#pragma acc declare create ( _isis_mod_var )
#pragma acc declare create ( _isis_mod )

/* component tofSource=TOF_monitor() [3] DECLARE */
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
  char     _name[256]; /* e.g. tofSource */
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
_class_TOF_monitor _tofSource_var;
_class_TOF_monitor *_tofSource = &_tofSource_var;
#pragma acc declare create ( _tofSource_var )
#pragma acc declare create ( _tofSource )

_class_Arm _armStraight1_var;
_class_Arm *_armStraight1 = &_armStraight1_var;
#pragma acc declare create ( _armStraight1_var )
#pragma acc declare create ( _armStraight1 )

/* component lamStart=L_monitor() [5] DECLARE */
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

/* component psdStart=PSD_monitor() [6] DECLARE */
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
  char     _name[256]; /* e.g. psdStart */
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
_class_PSD_monitor _psdStart_var;
_class_PSD_monitor *_psdStart = &_psdStart_var;
#pragma acc declare create ( _psdStart_var )
#pragma acc declare create ( _psdStart )

/* component guide_straight_1=Guide() [7] DECLARE */
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
  char     _name[256]; /* e.g. guide_straight_1 */
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
_class_Guide _guide_straight_1_var;
_class_Guide *_guide_straight_1 = &_guide_straight_1_var;
#pragma acc declare create ( _guide_straight_1_var )
#pragma acc declare create ( _guide_straight_1 )

_class_Arm _armChopper1_var;
_class_Arm *_armChopper1 = &_armChopper1_var;
#pragma acc declare create ( _armChopper1_var )
#pragma acc declare create ( _armChopper1 )

/* component chopper1=DiskChopper() [9] DECLARE */
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
  char     _name[256]; /* e.g. chopper1 */
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
_class_DiskChopper _chopper1_var;
_class_DiskChopper *_chopper1 = &_chopper1_var;
#pragma acc declare create ( _chopper1_var )
#pragma acc declare create ( _chopper1 )

_class_Arm _armCurved1_var;
_class_Arm *_armCurved1 = &_armCurved1_var;
#pragma acc declare create ( _armCurved1_var )
#pragma acc declare create ( _armCurved1 )

/* component guide_curved1=Guide_curved() [11] DECLARE */
/* Parameter definition for component type 'Guide_curved' */
struct _struct_Guide_curved_parameters {
  /* Component type 'Guide_curved' setting parameters */
  MCNUM w1;
  MCNUM h1;
  MCNUM l;
  MCNUM R0;
  MCNUM Qc;
  MCNUM alpha;
  MCNUM m;
  MCNUM W;
  MCNUM curvature;
}; /* _struct_Guide_curved_parameters */
typedef struct _struct_Guide_curved_parameters _class_Guide_curved_parameters;

/* Parameters for component type 'Guide_curved' */
struct _struct_Guide_curved {
  char     _name[256]; /* e.g. guide_curved1 */
  char     _type[256]; /* Guide_curved */
  long     _index; /* e.g. 2 index in TRACE list */
  Coords   _position_absolute;
  Coords   _position_relative; /* wrt PREVIOUS */
  Rotation _rotation_absolute;
  Rotation _rotation_relative; /* wrt PREVIOUS */
  int      _rotation_is_identity;
  _class_Guide_curved_parameters _parameters;
};
typedef struct _struct_Guide_curved _class_Guide_curved;
_class_Guide_curved _guide_curved1_var;
_class_Guide_curved *_guide_curved1 = &_guide_curved1_var;
#pragma acc declare create ( _guide_curved1_var )
#pragma acc declare create ( _guide_curved1 )

_class_Arm _armChopper2_var;
_class_Arm *_armChopper2 = &_armChopper2_var;
#pragma acc declare create ( _armChopper2_var )
#pragma acc declare create ( _armChopper2 )

_class_DiskChopper _chopper2_var;
_class_DiskChopper *_chopper2 = &_chopper2_var;
#pragma acc declare create ( _chopper2_var )
#pragma acc declare create ( _chopper2 )

_class_Arm _armCurved2_var;
_class_Arm *_armCurved2 = &_armCurved2_var;
#pragma acc declare create ( _armCurved2_var )
#pragma acc declare create ( _armCurved2 )

_class_Guide_curved _guide_curved2_var;
_class_Guide_curved *_guide_curved2 = &_guide_curved2_var;
#pragma acc declare create ( _guide_curved2_var )
#pragma acc declare create ( _guide_curved2 )

_class_Arm _armStraight2_var;
_class_Arm *_armStraight2 = &_armStraight2_var;
#pragma acc declare create ( _armStraight2_var )
#pragma acc declare create ( _armStraight2 )

_class_Guide _guide_straight_2_var;
_class_Guide *_guide_straight_2 = &_guide_straight_2_var;
#pragma acc declare create ( _guide_straight_2_var )
#pragma acc declare create ( _guide_straight_2 )

_class_Arm _armSuper_var;
_class_Arm *_armSuper = &_armSuper_var;
#pragma acc declare create ( _armSuper_var )
#pragma acc declare create ( _armSuper )

_class_Guide _guide_super_var;
_class_Guide *_guide_super = &_guide_super_var;
#pragma acc declare create ( _guide_super_var )
#pragma acc declare create ( _guide_super )

_class_Arm _armTarget_var;
_class_Arm *_armTarget = &_armTarget_var;
#pragma acc declare create ( _armTarget_var )
#pragma acc declare create ( _armTarget )

_class_L_monitor _lamTarget_var;
_class_L_monitor *_lamTarget = &_lamTarget_var;
#pragma acc declare create ( _lamTarget_var )
#pragma acc declare create ( _lamTarget )

_class_PSD_monitor _psdTarget_var;
_class_PSD_monitor *_psdTarget = &_psdTarget_var;
#pragma acc declare create ( _psdTarget_var )
#pragma acc declare create ( _psdTarget )

/* component vsample1=Incoherent() [23] DECLARE */
/* Parameter definition for component type 'Incoherent' */
struct _struct_Incoherent_parameters {
  /* Component type 'Incoherent' setting parameters */
  char geometry[16384];
  MCNUM radius;
  MCNUM xwidth;
  MCNUM yheight;
  MCNUM zdepth;
  MCNUM thickness;
  MCNUM target_x;
  MCNUM target_y;
  MCNUM target_z;
  MCNUM focus_r;
  MCNUM focus_xw;
  MCNUM focus_yh;
  MCNUM focus_aw;
  MCNUM focus_ah;
  long target_index;
  MCNUM pack;
  MCNUM p_interact;
  MCNUM f_QE;
  MCNUM gamma;
  MCNUM sigma_abs;
  MCNUM sigma_inc;
  MCNUM Vc;
  MCNUM concentric;
  MCNUM order;
  /* Component type 'Incoherent' private parameters */
  /* Component type 'Incoherent' DECLARE code stored as structure members */

  struct StructVarsInc VarsInc;
  off_struct offdata;

}; /* _struct_Incoherent_parameters */
typedef struct _struct_Incoherent_parameters _class_Incoherent_parameters;

/* Parameters for component type 'Incoherent' */
struct _struct_Incoherent {
  char     _name[256]; /* e.g. vsample1 */
  char     _type[256]; /* Incoherent */
  long     _index; /* e.g. 2 index in TRACE list */
  Coords   _position_absolute;
  Coords   _position_relative; /* wrt PREVIOUS */
  Rotation _rotation_absolute;
  Rotation _rotation_relative; /* wrt PREVIOUS */
  int      _rotation_is_identity;
  _class_Incoherent_parameters _parameters;
};
typedef struct _struct_Incoherent _class_Incoherent;
_class_Incoherent _vsample1_var;
_class_Incoherent *_vsample1 = &_vsample1_var;
#pragma acc declare create ( _vsample1_var )
#pragma acc declare create ( _vsample1 )

_class_Incoherent _vsample2_var;
_class_Incoherent *_vsample2 = &_vsample2_var;
#pragma acc declare create ( _vsample2_var )
#pragma acc declare create ( _vsample2 )

_class_Incoherent _vsample3_var;
_class_Incoherent *_vsample3 = &_vsample3_var;
#pragma acc declare create ( _vsample3_var )
#pragma acc declare create ( _vsample3 )

/* component beamStop=Beamstop() [26] DECLARE */
/* Parameter definition for component type 'Beamstop' */
struct _struct_Beamstop_parameters {
  /* Component type 'Beamstop' setting parameters */
  MCNUM xmin;
  MCNUM xmax;
  MCNUM ymin;
  MCNUM ymax;
  MCNUM xwidth;
  MCNUM yheight;
  MCNUM radius;
}; /* _struct_Beamstop_parameters */
typedef struct _struct_Beamstop_parameters _class_Beamstop_parameters;

/* Parameters for component type 'Beamstop' */
struct _struct_Beamstop {
  char     _name[256]; /* e.g. beamStop */
  char     _type[256]; /* Beamstop */
  long     _index; /* e.g. 2 index in TRACE list */
  Coords   _position_absolute;
  Coords   _position_relative; /* wrt PREVIOUS */
  Rotation _rotation_absolute;
  Rotation _rotation_relative; /* wrt PREVIOUS */
  int      _rotation_is_identity;
  _class_Beamstop_parameters _parameters;
};
typedef struct _struct_Beamstop _class_Beamstop;
_class_Beamstop _beamStop_var;
_class_Beamstop *_beamStop = &_beamStop_var;
#pragma acc declare create ( _beamStop_var )
#pragma acc declare create ( _beamStop )

_class_Arm _armAnalyzer_var;
_class_Arm *_armAnalyzer = &_armAnalyzer_var;
#pragma acc declare create ( _armAnalyzer_var )
#pragma acc declare create ( _armAnalyzer )

/* component graphite_analyser0=Monochromator_pol() [28] DECLARE */
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
  char     _name[256]; /* e.g. graphite_analyser0 */
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
_class_Monochromator_pol _graphite_analyser0_var;
_class_Monochromator_pol *_graphite_analyser0 = &_graphite_analyser0_var;
#pragma acc declare create ( _graphite_analyser0_var )
#pragma acc declare create ( _graphite_analyser0 )

_class_Monochromator_pol _graphite_analyser1_var;
_class_Monochromator_pol *_graphite_analyser1 = &_graphite_analyser1_var;
#pragma acc declare create ( _graphite_analyser1_var )
#pragma acc declare create ( _graphite_analyser1 )

_class_Monochromator_pol _graphite_analyser2_var;
_class_Monochromator_pol *_graphite_analyser2 = &_graphite_analyser2_var;
#pragma acc declare create ( _graphite_analyser2_var )
#pragma acc declare create ( _graphite_analyser2 )

_class_Monochromator_pol _graphite_analyser3_var;
_class_Monochromator_pol *_graphite_analyser3 = &_graphite_analyser3_var;
#pragma acc declare create ( _graphite_analyser3_var )
#pragma acc declare create ( _graphite_analyser3 )

_class_Monochromator_pol _graphite_analyser4_var;
_class_Monochromator_pol *_graphite_analyser4 = &_graphite_analyser4_var;
#pragma acc declare create ( _graphite_analyser4_var )
#pragma acc declare create ( _graphite_analyser4 )

_class_Monochromator_pol _graphite_analyser5_var;
_class_Monochromator_pol *_graphite_analyser5 = &_graphite_analyser5_var;
#pragma acc declare create ( _graphite_analyser5_var )
#pragma acc declare create ( _graphite_analyser5 )

_class_Monochromator_pol _graphite_analyser6_var;
_class_Monochromator_pol *_graphite_analyser6 = &_graphite_analyser6_var;
#pragma acc declare create ( _graphite_analyser6_var )
#pragma acc declare create ( _graphite_analyser6 )

_class_Monochromator_pol _graphite_analyser7_var;
_class_Monochromator_pol *_graphite_analyser7 = &_graphite_analyser7_var;
#pragma acc declare create ( _graphite_analyser7_var )
#pragma acc declare create ( _graphite_analyser7 )

_class_Monochromator_pol _graphite_analyser8_var;
_class_Monochromator_pol *_graphite_analyser8 = &_graphite_analyser8_var;
#pragma acc declare create ( _graphite_analyser8_var )
#pragma acc declare create ( _graphite_analyser8 )

_class_Monochromator_pol _graphite_analyser9_var;
_class_Monochromator_pol *_graphite_analyser9 = &_graphite_analyser9_var;
#pragma acc declare create ( _graphite_analyser9_var )
#pragma acc declare create ( _graphite_analyser9 )

_class_Monochromator_pol _graphite_analyser10_var;
_class_Monochromator_pol *_graphite_analyser10 = &_graphite_analyser10_var;
#pragma acc declare create ( _graphite_analyser10_var )
#pragma acc declare create ( _graphite_analyser10 )

_class_Monochromator_pol _graphite_analyser11_var;
_class_Monochromator_pol *_graphite_analyser11 = &_graphite_analyser11_var;
#pragma acc declare create ( _graphite_analyser11_var )
#pragma acc declare create ( _graphite_analyser11 )

_class_Monochromator_pol _graphite_analyser12_var;
_class_Monochromator_pol *_graphite_analyser12 = &_graphite_analyser12_var;
#pragma acc declare create ( _graphite_analyser12_var )
#pragma acc declare create ( _graphite_analyser12 )

_class_Monochromator_pol _graphite_analyser13_var;
_class_Monochromator_pol *_graphite_analyser13 = &_graphite_analyser13_var;
#pragma acc declare create ( _graphite_analyser13_var )
#pragma acc declare create ( _graphite_analyser13 )

_class_Monochromator_pol _graphite_analyser14_var;
_class_Monochromator_pol *_graphite_analyser14 = &_graphite_analyser14_var;
#pragma acc declare create ( _graphite_analyser14_var )
#pragma acc declare create ( _graphite_analyser14 )

_class_Monochromator_pol _graphite_analyser15_var;
_class_Monochromator_pol *_graphite_analyser15 = &_graphite_analyser15_var;
#pragma acc declare create ( _graphite_analyser15_var )
#pragma acc declare create ( _graphite_analyser15 )

_class_Monochromator_pol _graphite_analyser16_var;
_class_Monochromator_pol *_graphite_analyser16 = &_graphite_analyser16_var;
#pragma acc declare create ( _graphite_analyser16_var )
#pragma acc declare create ( _graphite_analyser16 )

_class_Monochromator_pol _graphite_analyser17_var;
_class_Monochromator_pol *_graphite_analyser17 = &_graphite_analyser17_var;
#pragma acc declare create ( _graphite_analyser17_var )
#pragma acc declare create ( _graphite_analyser17 )

_class_Monochromator_pol _graphite_analyser18_var;
_class_Monochromator_pol *_graphite_analyser18 = &_graphite_analyser18_var;
#pragma acc declare create ( _graphite_analyser18_var )
#pragma acc declare create ( _graphite_analyser18 )

_class_Monochromator_pol _graphite_analyser19_var;
_class_Monochromator_pol *_graphite_analyser19 = &_graphite_analyser19_var;
#pragma acc declare create ( _graphite_analyser19_var )
#pragma acc declare create ( _graphite_analyser19 )

_class_Monochromator_pol _graphite_analyser20_var;
_class_Monochromator_pol *_graphite_analyser20 = &_graphite_analyser20_var;
#pragma acc declare create ( _graphite_analyser20_var )
#pragma acc declare create ( _graphite_analyser20 )

_class_Monochromator_pol _graphite_analyser21_var;
_class_Monochromator_pol *_graphite_analyser21 = &_graphite_analyser21_var;
#pragma acc declare create ( _graphite_analyser21_var )
#pragma acc declare create ( _graphite_analyser21 )

_class_Monochromator_pol _graphite_analyser22_var;
_class_Monochromator_pol *_graphite_analyser22 = &_graphite_analyser22_var;
#pragma acc declare create ( _graphite_analyser22_var )
#pragma acc declare create ( _graphite_analyser22 )

_class_Monochromator_pol _graphite_analyser23_var;
_class_Monochromator_pol *_graphite_analyser23 = &_graphite_analyser23_var;
#pragma acc declare create ( _graphite_analyser23_var )
#pragma acc declare create ( _graphite_analyser23 )

_class_Monochromator_pol _graphite_analyser24_var;
_class_Monochromator_pol *_graphite_analyser24 = &_graphite_analyser24_var;
#pragma acc declare create ( _graphite_analyser24_var )
#pragma acc declare create ( _graphite_analyser24 )

_class_Monochromator_pol _graphite_analyser25_var;
_class_Monochromator_pol *_graphite_analyser25 = &_graphite_analyser25_var;
#pragma acc declare create ( _graphite_analyser25_var )
#pragma acc declare create ( _graphite_analyser25 )

_class_Monochromator_pol _graphite_analyser26_var;
_class_Monochromator_pol *_graphite_analyser26 = &_graphite_analyser26_var;
#pragma acc declare create ( _graphite_analyser26_var )
#pragma acc declare create ( _graphite_analyser26 )

_class_Monochromator_pol _graphite_analyser27_var;
_class_Monochromator_pol *_graphite_analyser27 = &_graphite_analyser27_var;
#pragma acc declare create ( _graphite_analyser27_var )
#pragma acc declare create ( _graphite_analyser27 )

_class_Monochromator_pol _graphite_analyser28_var;
_class_Monochromator_pol *_graphite_analyser28 = &_graphite_analyser28_var;
#pragma acc declare create ( _graphite_analyser28_var )
#pragma acc declare create ( _graphite_analyser28 )

_class_Monochromator_pol _graphite_analyser29_var;
_class_Monochromator_pol *_graphite_analyser29 = &_graphite_analyser29_var;
#pragma acc declare create ( _graphite_analyser29_var )
#pragma acc declare create ( _graphite_analyser29 )

_class_Monochromator_pol _graphite_analyser30_var;
_class_Monochromator_pol *_graphite_analyser30 = &_graphite_analyser30_var;
#pragma acc declare create ( _graphite_analyser30_var )
#pragma acc declare create ( _graphite_analyser30 )

_class_Monochromator_pol _graphite_analyser31_var;
_class_Monochromator_pol *_graphite_analyser31 = &_graphite_analyser31_var;
#pragma acc declare create ( _graphite_analyser31_var )
#pragma acc declare create ( _graphite_analyser31 )

_class_Monochromator_pol _graphite_analyser32_var;
_class_Monochromator_pol *_graphite_analyser32 = &_graphite_analyser32_var;
#pragma acc declare create ( _graphite_analyser32_var )
#pragma acc declare create ( _graphite_analyser32 )

_class_Monochromator_pol _graphite_analyser33_var;
_class_Monochromator_pol *_graphite_analyser33 = &_graphite_analyser33_var;
#pragma acc declare create ( _graphite_analyser33_var )
#pragma acc declare create ( _graphite_analyser33 )

_class_Monochromator_pol _graphite_analyser34_var;
_class_Monochromator_pol *_graphite_analyser34 = &_graphite_analyser34_var;
#pragma acc declare create ( _graphite_analyser34_var )
#pragma acc declare create ( _graphite_analyser34 )

_class_Monochromator_pol _graphite_analyser35_var;
_class_Monochromator_pol *_graphite_analyser35 = &_graphite_analyser35_var;
#pragma acc declare create ( _graphite_analyser35_var )
#pragma acc declare create ( _graphite_analyser35 )

_class_Monochromator_pol _graphite_analyser36_var;
_class_Monochromator_pol *_graphite_analyser36 = &_graphite_analyser36_var;
#pragma acc declare create ( _graphite_analyser36_var )
#pragma acc declare create ( _graphite_analyser36 )

_class_Monochromator_pol _graphite_analyser37_var;
_class_Monochromator_pol *_graphite_analyser37 = &_graphite_analyser37_var;
#pragma acc declare create ( _graphite_analyser37_var )
#pragma acc declare create ( _graphite_analyser37 )

_class_Monochromator_pol _graphite_analyser38_var;
_class_Monochromator_pol *_graphite_analyser38 = &_graphite_analyser38_var;
#pragma acc declare create ( _graphite_analyser38_var )
#pragma acc declare create ( _graphite_analyser38 )

_class_Monochromator_pol _graphite_analyser39_var;
_class_Monochromator_pol *_graphite_analyser39 = &_graphite_analyser39_var;
#pragma acc declare create ( _graphite_analyser39_var )
#pragma acc declare create ( _graphite_analyser39 )

_class_Arm _armHelp_var;
_class_Arm *_armHelp = &_armHelp_var;
#pragma acc declare create ( _armHelp_var )
#pragma acc declare create ( _armHelp )

_class_Arm _armDetector_var;
_class_Arm *_armDetector = &_armDetector_var;
#pragma acc declare create ( _armDetector_var )
#pragma acc declare create ( _armDetector )

_class_L_monitor _lamDetector_var;
_class_L_monitor *_lamDetector = &_lamDetector_var;
#pragma acc declare create ( _lamDetector_var )
#pragma acc declare create ( _lamDetector )

_class_PSD_monitor _psdDetector_var;
_class_PSD_monitor *_psdDetector = &_psdDetector_var;
#pragma acc declare create ( _psdDetector_var )
#pragma acc declare create ( _psdDetector )

_class_TOF_monitor _tofDetector_var;
_class_TOF_monitor *_tofDetector = &_tofDetector_var;
#pragma acc declare create ( _tofDetector_var )
#pragma acc declare create ( _tofDetector )

int mcNUMCOMP = 72;

/* User declarations from instrument definition. Can define functions. */
  double calcAlpha(double length, double radius) {
    // calculate angle of arm after curved guide
    return RAD2DEG * length/radius;
  }

  double calcX(double length, double radius) {
    // calculate position and angle of arm after curved guide
    double alpha = DEG2RAD * calcAlpha(length, radius);
    return radius*(1.0-cos(alpha));
  }

  double calcZ(double length, double radius) {
    // calculate position and angle of arm after curved guide
    double alpha = DEG2RAD * calcAlpha(length, radius);
    return radius*sin(alpha);
  }

  double calcT0(double lambda, double distance) {
    // calculate phaseoffset in s of chopper
    double k = 0;

    if(lambda<=0)
      return 0;

    k = 2*PI/lambda;
    return distance/(k*K2V);
  }

  // BEAM PARAMETERS
  const double MODSIZEX = 0.; // m (use default)
  const double MODSIZEY = 0.; // m (use default)
  double LAMBDAMIN  =   1.0; // - means use wavelength (AA)
  double LAMBDAMAX  =  22.0; // - means use wavelength (AA)

  // GUIDE SECTION PARAMETERS
  const double GUIDEWIDTH      =     0.043; // m
  const double GUIDEHEIGHT     =     0.065; // m
  const double GUIDEENDWIDTH   =     0.022; // m
  const double GUIDEENDHEIGHT  =     0.044; // m

  const double DMODGUIDE       =    1.6980; // m
  const double DSTRAIGHTGUIDE1 =    4.5440; // m (include 6.6 cm space window)
  const double DCHOPPER1       =    0.0404; // m
  const double DCURVEDGUIDE1   =    3.6800; // m
  const double DCHOPPER2       =    0.0404; // m
  const double DCURVEDGUIDE2   =   21.2977; // m
  const double CURVERADIUS     = 2050.0000; // m
  const double DSTRAIGHTGUIDE2 =    0.8656; // m (include 8.76cm+17.8cm space windows)
  const double DSUPERGUIDE     =    1.5063; // m
  const double DGUIDETARGET    =    0.2500; // m
  double DMODTARGET;

  // chopper parameters
  const double CHOPPERRADIUS  =     0.308; // m
  const double CHOPPERHEIGHT  =     0.091; // m

  double CHOPPER1OMEGA  =  50*2*PI; // rad/s
  double CHOPPER1WINDOW =       66; // deg

  double CHOPPER2OMEGA  = -50*2*PI; // rad/s
  double CHOPPER2WINDOW =       98; // deg

  // TARGET SECTION PARAMETERS
  const double VANADIUMHEIGHT      =  0.07; // m
  /* This variable helps us select which target to interact with*/
  double probTarget;

  // ANALYZER SECTION PARAMETERS
  const double ANARADIUS     = 0.9068; // m
  const double ANAHEIGHT     = 0.01; // m
  const double ANAWIDTH      = 0.01; // m
  const double ANAMOSAIC     = 0.8*60; // arc minutes
  //  const double anad          = 3.354; //
  const double ANAQ          = 2*PI/3.354; // AA-1

  // He DETECTOR PARAMETERS
  const double DETECTORDISTANCE = 0.67; //m
  // effective detectorwidth and height: sqrt(3.141*((1.27/2.0)**2))=1.125cm
  const double DETECTORWIDTH  = 0.01125; //m (half inch)
  const double DETECTORHEIGHT = 0.01125; //m (half inch)

  // Variable to help classify hits from different analyzer crystals
  int groupNumber = 0;

#undef compcurname
#undef compcurtype
#undef compcurindex
/* end of instrument 'ISIS_OSIRIS' and components DECLARE */

/* *****************************************************************************
* instrument 'ISIS_OSIRIS' and components INITIALISE
***************************************************************************** */

/* component armStart=Arm() SETTING, POSITION/ROTATION */
int _armStart_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_armStart_setpos] component armStart=Arm() SETTING [Arm:0]");
  stracpy(_armStart->_name, "armStart", 16384);
  stracpy(_armStart->_type, "Arm", 16384);
  _armStart->_index=1;
  /* component armStart=Arm() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(_armStart->_rotation_absolute,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_copy(_armStart->_rotation_relative, _armStart->_rotation_absolute);
    _armStart->_rotation_is_identity =  rot_test_identity(_armStart->_rotation_relative);
    _armStart->_position_absolute = coords_set(
      0, 0, 0);
    tc1 = coords_neg(_armStart->_position_absolute);
    _armStart->_position_relative = rot_apply(_armStart->_rotation_absolute, tc1);
  } /* armStart=Arm() AT ROTATED */
  DEBUG_COMPONENT("armStart", _armStart->_position_absolute, _armStart->_rotation_absolute);
  instrument->_position_absolute[1] = _armStart->_position_absolute;
  instrument->_position_relative[1] = _armStart->_position_relative;
  instrument->counter_N[1]  = instrument->counter_P[1] = instrument->counter_P2[1] = 0;
  instrument->counter_AbsorbProp[1]= 0;
  return(0);
} /* _armStart_setpos */

/* component isis_mod=ISIS_moderator() SETTING, POSITION/ROTATION */
int _isis_mod_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_isis_mod_setpos] component isis_mod=ISIS_moderator() SETTING [/usr/share/mcstas/3.0-dev/contrib/ISIS_moderator.comp:836]");
  stracpy(_isis_mod->_name, "isis_mod", 16384);
  stracpy(_isis_mod->_type, "ISIS_moderator", 16384);
  _isis_mod->_index=2;
  if("iris" && strlen("iris"))
    stracpy(_isis_mod->_parameters.Face, "iris" ? "iris" : "", 16384);
  else 
  _isis_mod->_parameters.Face[0]='\0';
  #define Face (_isis_mod->_parameters.Face)
  _isis_mod->_parameters.Emin = - LAMBDAMIN;
  #define Emin (_isis_mod->_parameters.Emin)
  _isis_mod->_parameters.Emax = - LAMBDAMAX;
  #define Emax (_isis_mod->_parameters.Emax)
  _isis_mod->_parameters.dist = DMODGUIDE;
  #define dist (_isis_mod->_parameters.dist)
  _isis_mod->_parameters.focus_xw = GUIDEWIDTH;
  #define focus_xw (_isis_mod->_parameters.focus_xw)
  _isis_mod->_parameters.focus_yh = GUIDEHEIGHT;
  #define focus_yh (_isis_mod->_parameters.focus_yh)
  _isis_mod->_parameters.xwidth = MODSIZEX;
  #define xwidth (_isis_mod->_parameters.xwidth)
  _isis_mod->_parameters.yheight = MODSIZEY;
  #define yheight (_isis_mod->_parameters.yheight)
  _isis_mod->_parameters.CAngle = 0.0;
  #define CAngle (_isis_mod->_parameters.CAngle)
  _isis_mod->_parameters.SAC = 1;
  #define SAC (_isis_mod->_parameters.SAC)
  _isis_mod->_parameters.Lmin = 0;
  #define Lmin (_isis_mod->_parameters.Lmin)
  _isis_mod->_parameters.Lmax = 0;
  #define Lmax (_isis_mod->_parameters.Lmax)
  _isis_mod->_parameters.target_index = + 1;
  #define target_index (_isis_mod->_parameters.target_index)
  _isis_mod->_parameters.verbose = 0;
  #define verbose (_isis_mod->_parameters.verbose)

  #define p_in (_isis_mod->_parameters.p_in)
  #define Tnpts (_isis_mod->_parameters.Tnpts)
  #define scaleSize (_isis_mod->_parameters.scaleSize)
  #define angleArea (_isis_mod->_parameters.angleArea)
  #define Nsim (_isis_mod->_parameters.Nsim)
  #define Ncount (_isis_mod->_parameters.Ncount)
  #define TS (_isis_mod->_parameters.TS)
  #define rtE0 (_isis_mod->_parameters.rtE0)
  #define rtE1 (_isis_mod->_parameters.rtE1)
  #define rtmodX (_isis_mod->_parameters.rtmodX)
  #define rtmodY (_isis_mod->_parameters.rtmodY)
  #define TargetStation (_isis_mod->_parameters.TargetStation)
  #define CurrentWeight (_isis_mod->_parameters.CurrentWeight)

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
  /* component isis_mod=ISIS_moderator() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _armStart->_rotation_absolute, _isis_mod->_rotation_absolute);
    rot_transpose(_armStart->_rotation_absolute, tr1);
    rot_mul(_isis_mod->_rotation_absolute, tr1, _isis_mod->_rotation_relative);
    _isis_mod->_rotation_is_identity =  rot_test_identity(_isis_mod->_rotation_relative);
    tc1 = coords_set(
      0, 0, 0);
    rot_transpose(_armStart->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _isis_mod->_position_absolute = coords_add(_armStart->_position_absolute, tc2);
    tc1 = coords_sub(_armStart->_position_absolute, _isis_mod->_position_absolute);
    _isis_mod->_position_relative = rot_apply(_isis_mod->_rotation_absolute, tc1);
  } /* isis_mod=ISIS_moderator() AT ROTATED */
  DEBUG_COMPONENT("isis_mod", _isis_mod->_position_absolute, _isis_mod->_rotation_absolute);
  instrument->_position_absolute[2] = _isis_mod->_position_absolute;
  instrument->_position_relative[2] = _isis_mod->_position_relative;
  instrument->counter_N[2]  = instrument->counter_P[2] = instrument->counter_P2[2] = 0;
  instrument->counter_AbsorbProp[2]= 0;
  return(0);
} /* _isis_mod_setpos */

/* component tofSource=TOF_monitor() SETTING, POSITION/ROTATION */
int _tofSource_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_tofSource_setpos] component tofSource=TOF_monitor() SETTING [/usr/share/mcstas/3.0-dev/monitors/TOF_monitor.comp:63]");
  stracpy(_tofSource->_name, "tofSource", 16384);
  stracpy(_tofSource->_type, "TOF_monitor", 16384);
  _tofSource->_index=3;
  _tofSource->_parameters.nt = 200;
  #define nt (_tofSource->_parameters.nt)
  if("tofSource.dat" && strlen("tofSource.dat"))
    stracpy(_tofSource->_parameters.filename, "tofSource.dat" ? "tofSource.dat" : "", 16384);
  else 
  _tofSource->_parameters.filename[0]='\0';
  #define filename (_tofSource->_parameters.filename)
  _tofSource->_parameters.xmin = -0.05;
  #define xmin (_tofSource->_parameters.xmin)
  _tofSource->_parameters.xmax = 0.05;
  #define xmax (_tofSource->_parameters.xmax)
  _tofSource->_parameters.ymin = -0.05;
  #define ymin (_tofSource->_parameters.ymin)
  _tofSource->_parameters.ymax = 0.05;
  #define ymax (_tofSource->_parameters.ymax)
  _tofSource->_parameters.xwidth = 0.20;
  #define xwidth (_tofSource->_parameters.xwidth)
  _tofSource->_parameters.yheight = 0.20;
  #define yheight (_tofSource->_parameters.yheight)
  _tofSource->_parameters.tmin = 0;
  #define tmin (_tofSource->_parameters.tmin)
  _tofSource->_parameters.tmax = 2000;
  #define tmax (_tofSource->_parameters.tmax)
  _tofSource->_parameters.dt = 1.0;
  #define dt (_tofSource->_parameters.dt)
  _tofSource->_parameters.restore_neutron = 0;
  #define restore_neutron (_tofSource->_parameters.restore_neutron)

  #define TOF_N (_tofSource->_parameters.TOF_N)
  #define TOF_p (_tofSource->_parameters.TOF_p)
  #define TOF_p2 (_tofSource->_parameters.TOF_p2)
  #define t_min (_tofSource->_parameters.t_min)
  #define t_max (_tofSource->_parameters.t_max)
  #define delta_t (_tofSource->_parameters.delta_t)

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
  /* component tofSource=TOF_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _isis_mod->_rotation_absolute, _tofSource->_rotation_absolute);
    rot_transpose(_isis_mod->_rotation_absolute, tr1);
    rot_mul(_tofSource->_rotation_absolute, tr1, _tofSource->_rotation_relative);
    _tofSource->_rotation_is_identity =  rot_test_identity(_tofSource->_rotation_relative);
    tc1 = coords_set(
      0, 0, 0);
    rot_transpose(_isis_mod->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _tofSource->_position_absolute = coords_add(_isis_mod->_position_absolute, tc2);
    tc1 = coords_sub(_isis_mod->_position_absolute, _tofSource->_position_absolute);
    _tofSource->_position_relative = rot_apply(_tofSource->_rotation_absolute, tc1);
  } /* tofSource=TOF_monitor() AT ROTATED */
  DEBUG_COMPONENT("tofSource", _tofSource->_position_absolute, _tofSource->_rotation_absolute);
  instrument->_position_absolute[3] = _tofSource->_position_absolute;
  instrument->_position_relative[3] = _tofSource->_position_relative;
  instrument->counter_N[3]  = instrument->counter_P[3] = instrument->counter_P2[3] = 0;
  instrument->counter_AbsorbProp[3]= 0;
  return(0);
} /* _tofSource_setpos */

/* component armStraight1=Arm() SETTING, POSITION/ROTATION */
int _armStraight1_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_armStraight1_setpos] component armStraight1=Arm() SETTING [Arm:0]");
  stracpy(_armStraight1->_name, "armStraight1", 16384);
  stracpy(_armStraight1->_type, "Arm", 16384);
  _armStraight1->_index=4;
  /* component armStraight1=Arm() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _armStart->_rotation_absolute, _armStraight1->_rotation_absolute);
    rot_transpose(_tofSource->_rotation_absolute, tr1);
    rot_mul(_armStraight1->_rotation_absolute, tr1, _armStraight1->_rotation_relative);
    _armStraight1->_rotation_is_identity =  rot_test_identity(_armStraight1->_rotation_relative);
    tc1 = coords_set(
      0, 0, DMODGUIDE);
    rot_transpose(_armStart->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _armStraight1->_position_absolute = coords_add(_armStart->_position_absolute, tc2);
    tc1 = coords_sub(_tofSource->_position_absolute, _armStraight1->_position_absolute);
    _armStraight1->_position_relative = rot_apply(_armStraight1->_rotation_absolute, tc1);
  } /* armStraight1=Arm() AT ROTATED */
  DEBUG_COMPONENT("armStraight1", _armStraight1->_position_absolute, _armStraight1->_rotation_absolute);
  instrument->_position_absolute[4] = _armStraight1->_position_absolute;
  instrument->_position_relative[4] = _armStraight1->_position_relative;
  instrument->counter_N[4]  = instrument->counter_P[4] = instrument->counter_P2[4] = 0;
  instrument->counter_AbsorbProp[4]= 0;
  return(0);
} /* _armStraight1_setpos */

/* component lamStart=L_monitor() SETTING, POSITION/ROTATION */
int _lamStart_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_lamStart_setpos] component lamStart=L_monitor() SETTING [/usr/share/mcstas/3.0-dev/monitors/L_monitor.comp:65]");
  stracpy(_lamStart->_name, "lamStart", 16384);
  stracpy(_lamStart->_type, "L_monitor", 16384);
  _lamStart->_index=5;
  _lamStart->_parameters.nL = 120;
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
  _lamStart->_parameters.Lmin = 0.98 * LAMBDAMIN;
  #define Lmin (_lamStart->_parameters.Lmin)
  _lamStart->_parameters.Lmax = 1.02 * LAMBDAMAX;
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
    rot_mul(tr1, _armStraight1->_rotation_absolute, _lamStart->_rotation_absolute);
    rot_transpose(_armStraight1->_rotation_absolute, tr1);
    rot_mul(_lamStart->_rotation_absolute, tr1, _lamStart->_rotation_relative);
    _lamStart->_rotation_is_identity =  rot_test_identity(_lamStart->_rotation_relative);
    tc1 = coords_set(
      0, 0, 0);
    rot_transpose(_armStraight1->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _lamStart->_position_absolute = coords_add(_armStraight1->_position_absolute, tc2);
    tc1 = coords_sub(_armStraight1->_position_absolute, _lamStart->_position_absolute);
    _lamStart->_position_relative = rot_apply(_lamStart->_rotation_absolute, tc1);
  } /* lamStart=L_monitor() AT ROTATED */
  DEBUG_COMPONENT("lamStart", _lamStart->_position_absolute, _lamStart->_rotation_absolute);
  instrument->_position_absolute[5] = _lamStart->_position_absolute;
  instrument->_position_relative[5] = _lamStart->_position_relative;
  instrument->counter_N[5]  = instrument->counter_P[5] = instrument->counter_P2[5] = 0;
  instrument->counter_AbsorbProp[5]= 0;
  return(0);
} /* _lamStart_setpos */

/* component psdStart=PSD_monitor() SETTING, POSITION/ROTATION */
int _psdStart_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_psdStart_setpos] component psdStart=PSD_monitor() SETTING [/usr/share/mcstas/3.0-dev/monitors/PSD_monitor.comp:64]");
  stracpy(_psdStart->_name, "psdStart", 16384);
  stracpy(_psdStart->_type, "PSD_monitor", 16384);
  _psdStart->_index=6;
  _psdStart->_parameters.nx = 60;
  #define nx (_psdStart->_parameters.nx)
  _psdStart->_parameters.ny = 60;
  #define ny (_psdStart->_parameters.ny)
  if("psdStart.dat" && strlen("psdStart.dat"))
    stracpy(_psdStart->_parameters.filename, "psdStart.dat" ? "psdStart.dat" : "", 16384);
  else 
  _psdStart->_parameters.filename[0]='\0';
  #define filename (_psdStart->_parameters.filename)
  _psdStart->_parameters.xmin = -0.05;
  #define xmin (_psdStart->_parameters.xmin)
  _psdStart->_parameters.xmax = 0.05;
  #define xmax (_psdStart->_parameters.xmax)
  _psdStart->_parameters.ymin = -0.05;
  #define ymin (_psdStart->_parameters.ymin)
  _psdStart->_parameters.ymax = 0.05;
  #define ymax (_psdStart->_parameters.ymax)
  _psdStart->_parameters.xwidth = 0.10;
  #define xwidth (_psdStart->_parameters.xwidth)
  _psdStart->_parameters.yheight = 0.10;
  #define yheight (_psdStart->_parameters.yheight)
  _psdStart->_parameters.restore_neutron = 0;
  #define restore_neutron (_psdStart->_parameters.restore_neutron)

  #define PSD_N (_psdStart->_parameters.PSD_N)
  #define PSD_p (_psdStart->_parameters.PSD_p)
  #define PSD_p2 (_psdStart->_parameters.PSD_p2)

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
  /* component psdStart=PSD_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _armStraight1->_rotation_absolute, _psdStart->_rotation_absolute);
    rot_transpose(_lamStart->_rotation_absolute, tr1);
    rot_mul(_psdStart->_rotation_absolute, tr1, _psdStart->_rotation_relative);
    _psdStart->_rotation_is_identity =  rot_test_identity(_psdStart->_rotation_relative);
    tc1 = coords_set(
      0, 0, 0);
    rot_transpose(_armStraight1->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _psdStart->_position_absolute = coords_add(_armStraight1->_position_absolute, tc2);
    tc1 = coords_sub(_lamStart->_position_absolute, _psdStart->_position_absolute);
    _psdStart->_position_relative = rot_apply(_psdStart->_rotation_absolute, tc1);
  } /* psdStart=PSD_monitor() AT ROTATED */
  DEBUG_COMPONENT("psdStart", _psdStart->_position_absolute, _psdStart->_rotation_absolute);
  instrument->_position_absolute[6] = _psdStart->_position_absolute;
  instrument->_position_relative[6] = _psdStart->_position_relative;
  instrument->counter_N[6]  = instrument->counter_P[6] = instrument->counter_P2[6] = 0;
  instrument->counter_AbsorbProp[6]= 0;
  return(0);
} /* _psdStart_setpos */

/* component guide_straight_1=Guide() SETTING, POSITION/ROTATION */
int _guide_straight_1_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_guide_straight_1_setpos] component guide_straight_1=Guide() SETTING [/usr/share/mcstas/3.0-dev/optics/Guide.comp:73]");
  stracpy(_guide_straight_1->_name, "guide_straight_1", 16384);
  stracpy(_guide_straight_1->_type, "Guide", 16384);
  _guide_straight_1->_index=7;
  _guide_straight_1->_parameters.reflect[0]='\0';
  #define reflect (_guide_straight_1->_parameters.reflect)
  _guide_straight_1->_parameters.w1 = GUIDEWIDTH;
  #define w1 (_guide_straight_1->_parameters.w1)
  _guide_straight_1->_parameters.h1 = GUIDEHEIGHT;
  #define h1 (_guide_straight_1->_parameters.h1)
  _guide_straight_1->_parameters.w2 = GUIDEWIDTH;
  #define w2 (_guide_straight_1->_parameters.w2)
  _guide_straight_1->_parameters.h2 = GUIDEHEIGHT;
  #define h2 (_guide_straight_1->_parameters.h2)
  _guide_straight_1->_parameters.l = DSTRAIGHTGUIDE1;
  #define l (_guide_straight_1->_parameters.l)
  _guide_straight_1->_parameters.R0 = instrument->_parameters._GUIDEREFLECTIVITY;
  #define R0 (_guide_straight_1->_parameters.R0)
  _guide_straight_1->_parameters.Qc = 0.0219;
  #define Qc (_guide_straight_1->_parameters.Qc)
  _guide_straight_1->_parameters.alpha = 6.07;
  #define alpha (_guide_straight_1->_parameters.alpha)
  _guide_straight_1->_parameters.m = 2.0;
  #define m (_guide_straight_1->_parameters.m)
  _guide_straight_1->_parameters.W = 0.003;
  #define W (_guide_straight_1->_parameters.W)

  #define pTable (_guide_straight_1->_parameters.pTable)

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
  /* component guide_straight_1=Guide() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _armStraight1->_rotation_absolute, _guide_straight_1->_rotation_absolute);
    rot_transpose(_psdStart->_rotation_absolute, tr1);
    rot_mul(_guide_straight_1->_rotation_absolute, tr1, _guide_straight_1->_rotation_relative);
    _guide_straight_1->_rotation_is_identity =  rot_test_identity(_guide_straight_1->_rotation_relative);
    tc1 = coords_set(
      0, 0, 0);
    rot_transpose(_armStraight1->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _guide_straight_1->_position_absolute = coords_add(_armStraight1->_position_absolute, tc2);
    tc1 = coords_sub(_psdStart->_position_absolute, _guide_straight_1->_position_absolute);
    _guide_straight_1->_position_relative = rot_apply(_guide_straight_1->_rotation_absolute, tc1);
  } /* guide_straight_1=Guide() AT ROTATED */
  DEBUG_COMPONENT("guide_straight_1", _guide_straight_1->_position_absolute, _guide_straight_1->_rotation_absolute);
  instrument->_position_absolute[7] = _guide_straight_1->_position_absolute;
  instrument->_position_relative[7] = _guide_straight_1->_position_relative;
  instrument->counter_N[7]  = instrument->counter_P[7] = instrument->counter_P2[7] = 0;
  instrument->counter_AbsorbProp[7]= 0;
  return(0);
} /* _guide_straight_1_setpos */

/* component armChopper1=Arm() SETTING, POSITION/ROTATION */
int _armChopper1_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_armChopper1_setpos] component armChopper1=Arm() SETTING [Arm:0]");
  stracpy(_armChopper1->_name, "armChopper1", 16384);
  stracpy(_armChopper1->_type, "Arm", 16384);
  _armChopper1->_index=8;
  /* component armChopper1=Arm() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _armStraight1->_rotation_absolute, _armChopper1->_rotation_absolute);
    rot_transpose(_guide_straight_1->_rotation_absolute, tr1);
    rot_mul(_armChopper1->_rotation_absolute, tr1, _armChopper1->_rotation_relative);
    _armChopper1->_rotation_is_identity =  rot_test_identity(_armChopper1->_rotation_relative);
    tc1 = coords_set(
      0, 0, DSTRAIGHTGUIDE1);
    rot_transpose(_armStraight1->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _armChopper1->_position_absolute = coords_add(_armStraight1->_position_absolute, tc2);
    tc1 = coords_sub(_guide_straight_1->_position_absolute, _armChopper1->_position_absolute);
    _armChopper1->_position_relative = rot_apply(_armChopper1->_rotation_absolute, tc1);
  } /* armChopper1=Arm() AT ROTATED */
  DEBUG_COMPONENT("armChopper1", _armChopper1->_position_absolute, _armChopper1->_rotation_absolute);
  instrument->_position_absolute[8] = _armChopper1->_position_absolute;
  instrument->_position_relative[8] = _armChopper1->_position_relative;
  instrument->counter_N[8]  = instrument->counter_P[8] = instrument->counter_P2[8] = 0;
  instrument->counter_AbsorbProp[8]= 0;
  return(0);
} /* _armChopper1_setpos */

/* component chopper1=DiskChopper() SETTING, POSITION/ROTATION */
int _chopper1_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_chopper1_setpos] component chopper1=DiskChopper() SETTING [/usr/share/mcstas/3.0-dev/optics/DiskChopper.comp:67]");
  stracpy(_chopper1->_name, "chopper1", 16384);
  stracpy(_chopper1->_type, "DiskChopper", 16384);
  _chopper1->_index=9;
  _chopper1->_parameters.theta_0 = CHOPPER1WINDOW;
  #define theta_0 (_chopper1->_parameters.theta_0)
  _chopper1->_parameters.radius = CHOPPERRADIUS;
  #define radius (_chopper1->_parameters.radius)
  _chopper1->_parameters.yheight = CHOPPERHEIGHT;
  #define yheight (_chopper1->_parameters.yheight)
  _chopper1->_parameters.nu = CHOPPER1OMEGA / ( 2 * PI );
  #define nu (_chopper1->_parameters.nu)
  _chopper1->_parameters.nslit = 1;
  #define nslit (_chopper1->_parameters.nslit)
  _chopper1->_parameters.jitter = 0;
  #define jitter (_chopper1->_parameters.jitter)
  _chopper1->_parameters.delay = calcT0 ( instrument->_parameters._LAMBDA , DMODGUIDE + DSTRAIGHTGUIDE1 );
  #define delay (_chopper1->_parameters.delay)
  _chopper1->_parameters.isfirst = 0;
  #define isfirst (_chopper1->_parameters.isfirst)
  _chopper1->_parameters.n_pulse = 1;
  #define n_pulse (_chopper1->_parameters.n_pulse)
  _chopper1->_parameters.abs_out = 1;
  #define abs_out (_chopper1->_parameters.abs_out)
  _chopper1->_parameters.phase = 0;
  #define phase (_chopper1->_parameters.phase)
  _chopper1->_parameters.xwidth = 0;
  #define xwidth (_chopper1->_parameters.xwidth)
  _chopper1->_parameters.verbose = 0;
  #define verbose (_chopper1->_parameters.verbose)

  #define Tg (_chopper1->_parameters.Tg)
  #define To (_chopper1->_parameters.To)
  #define delta_y (_chopper1->_parameters.delta_y)
  #define height (_chopper1->_parameters.height)
  #define omega (_chopper1->_parameters.omega)

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
  /* component chopper1=DiskChopper() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _armChopper1->_rotation_absolute, _chopper1->_rotation_absolute);
    rot_transpose(_armChopper1->_rotation_absolute, tr1);
    rot_mul(_chopper1->_rotation_absolute, tr1, _chopper1->_rotation_relative);
    _chopper1->_rotation_is_identity =  rot_test_identity(_chopper1->_rotation_relative);
    tc1 = coords_set(
      0, 0, 0);
    rot_transpose(_armChopper1->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _chopper1->_position_absolute = coords_add(_armChopper1->_position_absolute, tc2);
    tc1 = coords_sub(_armChopper1->_position_absolute, _chopper1->_position_absolute);
    _chopper1->_position_relative = rot_apply(_chopper1->_rotation_absolute, tc1);
  } /* chopper1=DiskChopper() AT ROTATED */
  DEBUG_COMPONENT("chopper1", _chopper1->_position_absolute, _chopper1->_rotation_absolute);
  instrument->_position_absolute[9] = _chopper1->_position_absolute;
  instrument->_position_relative[9] = _chopper1->_position_relative;
  instrument->counter_N[9]  = instrument->counter_P[9] = instrument->counter_P2[9] = 0;
  instrument->counter_AbsorbProp[9]= 0;
  return(0);
} /* _chopper1_setpos */

/* component armCurved1=Arm() SETTING, POSITION/ROTATION */
int _armCurved1_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_armCurved1_setpos] component armCurved1=Arm() SETTING [Arm:0]");
  stracpy(_armCurved1->_name, "armCurved1", 16384);
  stracpy(_armCurved1->_type, "Arm", 16384);
  _armCurved1->_index=10;
  /* component armCurved1=Arm() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _armChopper1->_rotation_absolute, _armCurved1->_rotation_absolute);
    rot_transpose(_chopper1->_rotation_absolute, tr1);
    rot_mul(_armCurved1->_rotation_absolute, tr1, _armCurved1->_rotation_relative);
    _armCurved1->_rotation_is_identity =  rot_test_identity(_armCurved1->_rotation_relative);
    tc1 = coords_set(
      0, 0, DCHOPPER1);
    rot_transpose(_armChopper1->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _armCurved1->_position_absolute = coords_add(_armChopper1->_position_absolute, tc2);
    tc1 = coords_sub(_chopper1->_position_absolute, _armCurved1->_position_absolute);
    _armCurved1->_position_relative = rot_apply(_armCurved1->_rotation_absolute, tc1);
  } /* armCurved1=Arm() AT ROTATED */
  DEBUG_COMPONENT("armCurved1", _armCurved1->_position_absolute, _armCurved1->_rotation_absolute);
  instrument->_position_absolute[10] = _armCurved1->_position_absolute;
  instrument->_position_relative[10] = _armCurved1->_position_relative;
  instrument->counter_N[10]  = instrument->counter_P[10] = instrument->counter_P2[10] = 0;
  instrument->counter_AbsorbProp[10]= 0;
  return(0);
} /* _armCurved1_setpos */

/* component guide_curved1=Guide_curved() SETTING, POSITION/ROTATION */
int _guide_curved1_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_guide_curved1_setpos] component guide_curved1=Guide_curved() SETTING [/usr/share/mcstas/3.0-dev/contrib/Guide_curved.comp:60]");
  stracpy(_guide_curved1->_name, "guide_curved1", 16384);
  stracpy(_guide_curved1->_type, "Guide_curved", 16384);
  _guide_curved1->_index=11;
  _guide_curved1->_parameters.w1 = GUIDEWIDTH;
  #define w1 (_guide_curved1->_parameters.w1)
  _guide_curved1->_parameters.h1 = GUIDEHEIGHT;
  #define h1 (_guide_curved1->_parameters.h1)
  _guide_curved1->_parameters.l = DCURVEDGUIDE1;
  #define l (_guide_curved1->_parameters.l)
  _guide_curved1->_parameters.R0 = instrument->_parameters._GUIDEREFLECTIVITY;
  #define R0 (_guide_curved1->_parameters.R0)
  _guide_curved1->_parameters.Qc = 0.0218;
  #define Qc (_guide_curved1->_parameters.Qc)
  _guide_curved1->_parameters.alpha = 4.38;
  #define alpha (_guide_curved1->_parameters.alpha)
  _guide_curved1->_parameters.m = 2.0;
  #define m (_guide_curved1->_parameters.m)
  _guide_curved1->_parameters.W = 0.003;
  #define W (_guide_curved1->_parameters.W)
  _guide_curved1->_parameters.curvature = CURVERADIUS;
  #define curvature (_guide_curved1->_parameters.curvature)

  #undef w1
  #undef h1
  #undef l
  #undef R0
  #undef Qc
  #undef alpha
  #undef m
  #undef W
  #undef curvature
  /* component guide_curved1=Guide_curved() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _armCurved1->_rotation_absolute, _guide_curved1->_rotation_absolute);
    rot_transpose(_armCurved1->_rotation_absolute, tr1);
    rot_mul(_guide_curved1->_rotation_absolute, tr1, _guide_curved1->_rotation_relative);
    _guide_curved1->_rotation_is_identity =  rot_test_identity(_guide_curved1->_rotation_relative);
    tc1 = coords_set(
      0, 0, 0);
    rot_transpose(_armCurved1->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _guide_curved1->_position_absolute = coords_add(_armCurved1->_position_absolute, tc2);
    tc1 = coords_sub(_armCurved1->_position_absolute, _guide_curved1->_position_absolute);
    _guide_curved1->_position_relative = rot_apply(_guide_curved1->_rotation_absolute, tc1);
  } /* guide_curved1=Guide_curved() AT ROTATED */
  DEBUG_COMPONENT("guide_curved1", _guide_curved1->_position_absolute, _guide_curved1->_rotation_absolute);
  instrument->_position_absolute[11] = _guide_curved1->_position_absolute;
  instrument->_position_relative[11] = _guide_curved1->_position_relative;
  instrument->counter_N[11]  = instrument->counter_P[11] = instrument->counter_P2[11] = 0;
  instrument->counter_AbsorbProp[11]= 0;
  return(0);
} /* _guide_curved1_setpos */

/* component armChopper2=Arm() SETTING, POSITION/ROTATION */
int _armChopper2_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_armChopper2_setpos] component armChopper2=Arm() SETTING [Arm:0]");
  stracpy(_armChopper2->_name, "armChopper2", 16384);
  stracpy(_armChopper2->_type, "Arm", 16384);
  _armChopper2->_index=12;
  /* component armChopper2=Arm() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (calcAlpha ( DCURVEDGUIDE1 , CURVERADIUS ))*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _armCurved1->_rotation_absolute, _armChopper2->_rotation_absolute);
    rot_transpose(_guide_curved1->_rotation_absolute, tr1);
    rot_mul(_armChopper2->_rotation_absolute, tr1, _armChopper2->_rotation_relative);
    _armChopper2->_rotation_is_identity =  rot_test_identity(_armChopper2->_rotation_relative);
    tc1 = coords_set(
      calcX ( DCURVEDGUIDE1 , CURVERADIUS ), 0, calcZ ( DCURVEDGUIDE1 , CURVERADIUS ));
    rot_transpose(_armCurved1->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _armChopper2->_position_absolute = coords_add(_armCurved1->_position_absolute, tc2);
    tc1 = coords_sub(_guide_curved1->_position_absolute, _armChopper2->_position_absolute);
    _armChopper2->_position_relative = rot_apply(_armChopper2->_rotation_absolute, tc1);
  } /* armChopper2=Arm() AT ROTATED */
  DEBUG_COMPONENT("armChopper2", _armChopper2->_position_absolute, _armChopper2->_rotation_absolute);
  instrument->_position_absolute[12] = _armChopper2->_position_absolute;
  instrument->_position_relative[12] = _armChopper2->_position_relative;
  instrument->counter_N[12]  = instrument->counter_P[12] = instrument->counter_P2[12] = 0;
  instrument->counter_AbsorbProp[12]= 0;
  return(0);
} /* _armChopper2_setpos */

/* component chopper2=DiskChopper() SETTING, POSITION/ROTATION */
int _chopper2_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_chopper2_setpos] component chopper2=DiskChopper() SETTING [/usr/share/mcstas/3.0-dev/optics/DiskChopper.comp:67]");
  stracpy(_chopper2->_name, "chopper2", 16384);
  stracpy(_chopper2->_type, "DiskChopper", 16384);
  _chopper2->_index=13;
  _chopper2->_parameters.theta_0 = CHOPPER2WINDOW;
  #define theta_0 (_chopper2->_parameters.theta_0)
  _chopper2->_parameters.radius = CHOPPERRADIUS;
  #define radius (_chopper2->_parameters.radius)
  _chopper2->_parameters.yheight = CHOPPERHEIGHT;
  #define yheight (_chopper2->_parameters.yheight)
  _chopper2->_parameters.nu = CHOPPER2OMEGA / ( 2 * PI );
  #define nu (_chopper2->_parameters.nu)
  _chopper2->_parameters.nslit = 1;
  #define nslit (_chopper2->_parameters.nslit)
  _chopper2->_parameters.jitter = 0;
  #define jitter (_chopper2->_parameters.jitter)
  _chopper2->_parameters.delay = calcT0 ( instrument->_parameters._LAMBDA , DMODGUIDE + DSTRAIGHTGUIDE1 + DCHOPPER1 + DCURVEDGUIDE1 );
  #define delay (_chopper2->_parameters.delay)
  _chopper2->_parameters.isfirst = 0;
  #define isfirst (_chopper2->_parameters.isfirst)
  _chopper2->_parameters.n_pulse = 1;
  #define n_pulse (_chopper2->_parameters.n_pulse)
  _chopper2->_parameters.abs_out = 1;
  #define abs_out (_chopper2->_parameters.abs_out)
  _chopper2->_parameters.phase = 0;
  #define phase (_chopper2->_parameters.phase)
  _chopper2->_parameters.xwidth = 0;
  #define xwidth (_chopper2->_parameters.xwidth)
  _chopper2->_parameters.verbose = 0;
  #define verbose (_chopper2->_parameters.verbose)

  #define Tg (_chopper2->_parameters.Tg)
  #define To (_chopper2->_parameters.To)
  #define delta_y (_chopper2->_parameters.delta_y)
  #define height (_chopper2->_parameters.height)
  #define omega (_chopper2->_parameters.omega)

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
  /* component chopper2=DiskChopper() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _armChopper2->_rotation_absolute, _chopper2->_rotation_absolute);
    rot_transpose(_armChopper2->_rotation_absolute, tr1);
    rot_mul(_chopper2->_rotation_absolute, tr1, _chopper2->_rotation_relative);
    _chopper2->_rotation_is_identity =  rot_test_identity(_chopper2->_rotation_relative);
    tc1 = coords_set(
      0, 0, 0);
    rot_transpose(_armChopper2->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _chopper2->_position_absolute = coords_add(_armChopper2->_position_absolute, tc2);
    tc1 = coords_sub(_armChopper2->_position_absolute, _chopper2->_position_absolute);
    _chopper2->_position_relative = rot_apply(_chopper2->_rotation_absolute, tc1);
  } /* chopper2=DiskChopper() AT ROTATED */
  DEBUG_COMPONENT("chopper2", _chopper2->_position_absolute, _chopper2->_rotation_absolute);
  instrument->_position_absolute[13] = _chopper2->_position_absolute;
  instrument->_position_relative[13] = _chopper2->_position_relative;
  instrument->counter_N[13]  = instrument->counter_P[13] = instrument->counter_P2[13] = 0;
  instrument->counter_AbsorbProp[13]= 0;
  return(0);
} /* _chopper2_setpos */

/* component armCurved2=Arm() SETTING, POSITION/ROTATION */
int _armCurved2_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_armCurved2_setpos] component armCurved2=Arm() SETTING [Arm:0]");
  stracpy(_armCurved2->_name, "armCurved2", 16384);
  stracpy(_armCurved2->_type, "Arm", 16384);
  _armCurved2->_index=14;
  /* component armCurved2=Arm() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _armChopper2->_rotation_absolute, _armCurved2->_rotation_absolute);
    rot_transpose(_chopper2->_rotation_absolute, tr1);
    rot_mul(_armCurved2->_rotation_absolute, tr1, _armCurved2->_rotation_relative);
    _armCurved2->_rotation_is_identity =  rot_test_identity(_armCurved2->_rotation_relative);
    tc1 = coords_set(
      0, 0, DCHOPPER2);
    rot_transpose(_armChopper2->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _armCurved2->_position_absolute = coords_add(_armChopper2->_position_absolute, tc2);
    tc1 = coords_sub(_chopper2->_position_absolute, _armCurved2->_position_absolute);
    _armCurved2->_position_relative = rot_apply(_armCurved2->_rotation_absolute, tc1);
  } /* armCurved2=Arm() AT ROTATED */
  DEBUG_COMPONENT("armCurved2", _armCurved2->_position_absolute, _armCurved2->_rotation_absolute);
  instrument->_position_absolute[14] = _armCurved2->_position_absolute;
  instrument->_position_relative[14] = _armCurved2->_position_relative;
  instrument->counter_N[14]  = instrument->counter_P[14] = instrument->counter_P2[14] = 0;
  instrument->counter_AbsorbProp[14]= 0;
  return(0);
} /* _armCurved2_setpos */

/* component guide_curved2=Guide_curved() SETTING, POSITION/ROTATION */
int _guide_curved2_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_guide_curved2_setpos] component guide_curved2=Guide_curved() SETTING [/usr/share/mcstas/3.0-dev/contrib/Guide_curved.comp:60]");
  stracpy(_guide_curved2->_name, "guide_curved2", 16384);
  stracpy(_guide_curved2->_type, "Guide_curved", 16384);
  _guide_curved2->_index=15;
  _guide_curved2->_parameters.w1 = GUIDEWIDTH;
  #define w1 (_guide_curved2->_parameters.w1)
  _guide_curved2->_parameters.h1 = GUIDEHEIGHT;
  #define h1 (_guide_curved2->_parameters.h1)
  _guide_curved2->_parameters.l = DCURVEDGUIDE2;
  #define l (_guide_curved2->_parameters.l)
  _guide_curved2->_parameters.R0 = instrument->_parameters._GUIDEREFLECTIVITY;
  #define R0 (_guide_curved2->_parameters.R0)
  _guide_curved2->_parameters.Qc = 0.0218;
  #define Qc (_guide_curved2->_parameters.Qc)
  _guide_curved2->_parameters.alpha = 4.38;
  #define alpha (_guide_curved2->_parameters.alpha)
  _guide_curved2->_parameters.m = 2.0;
  #define m (_guide_curved2->_parameters.m)
  _guide_curved2->_parameters.W = 0.003;
  #define W (_guide_curved2->_parameters.W)
  _guide_curved2->_parameters.curvature = CURVERADIUS;
  #define curvature (_guide_curved2->_parameters.curvature)

  #undef w1
  #undef h1
  #undef l
  #undef R0
  #undef Qc
  #undef alpha
  #undef m
  #undef W
  #undef curvature
  /* component guide_curved2=Guide_curved() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _armCurved2->_rotation_absolute, _guide_curved2->_rotation_absolute);
    rot_transpose(_armCurved2->_rotation_absolute, tr1);
    rot_mul(_guide_curved2->_rotation_absolute, tr1, _guide_curved2->_rotation_relative);
    _guide_curved2->_rotation_is_identity =  rot_test_identity(_guide_curved2->_rotation_relative);
    tc1 = coords_set(
      0, 0, 0);
    rot_transpose(_armCurved2->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _guide_curved2->_position_absolute = coords_add(_armCurved2->_position_absolute, tc2);
    tc1 = coords_sub(_armCurved2->_position_absolute, _guide_curved2->_position_absolute);
    _guide_curved2->_position_relative = rot_apply(_guide_curved2->_rotation_absolute, tc1);
  } /* guide_curved2=Guide_curved() AT ROTATED */
  DEBUG_COMPONENT("guide_curved2", _guide_curved2->_position_absolute, _guide_curved2->_rotation_absolute);
  instrument->_position_absolute[15] = _guide_curved2->_position_absolute;
  instrument->_position_relative[15] = _guide_curved2->_position_relative;
  instrument->counter_N[15]  = instrument->counter_P[15] = instrument->counter_P2[15] = 0;
  instrument->counter_AbsorbProp[15]= 0;
  return(0);
} /* _guide_curved2_setpos */

/* component armStraight2=Arm() SETTING, POSITION/ROTATION */
int _armStraight2_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_armStraight2_setpos] component armStraight2=Arm() SETTING [Arm:0]");
  stracpy(_armStraight2->_name, "armStraight2", 16384);
  stracpy(_armStraight2->_type, "Arm", 16384);
  _armStraight2->_index=16;
  /* component armStraight2=Arm() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (calcAlpha ( DCURVEDGUIDE2 , CURVERADIUS ))*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _armCurved2->_rotation_absolute, _armStraight2->_rotation_absolute);
    rot_transpose(_guide_curved2->_rotation_absolute, tr1);
    rot_mul(_armStraight2->_rotation_absolute, tr1, _armStraight2->_rotation_relative);
    _armStraight2->_rotation_is_identity =  rot_test_identity(_armStraight2->_rotation_relative);
    tc1 = coords_set(
      calcX ( DCURVEDGUIDE2 , CURVERADIUS ), 0, calcZ ( DCURVEDGUIDE2 , CURVERADIUS ));
    rot_transpose(_armCurved2->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _armStraight2->_position_absolute = coords_add(_armCurved2->_position_absolute, tc2);
    tc1 = coords_sub(_guide_curved2->_position_absolute, _armStraight2->_position_absolute);
    _armStraight2->_position_relative = rot_apply(_armStraight2->_rotation_absolute, tc1);
  } /* armStraight2=Arm() AT ROTATED */
  DEBUG_COMPONENT("armStraight2", _armStraight2->_position_absolute, _armStraight2->_rotation_absolute);
  instrument->_position_absolute[16] = _armStraight2->_position_absolute;
  instrument->_position_relative[16] = _armStraight2->_position_relative;
  instrument->counter_N[16]  = instrument->counter_P[16] = instrument->counter_P2[16] = 0;
  instrument->counter_AbsorbProp[16]= 0;
  return(0);
} /* _armStraight2_setpos */

/* component guide_straight_2=Guide() SETTING, POSITION/ROTATION */
int _guide_straight_2_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_guide_straight_2_setpos] component guide_straight_2=Guide() SETTING [/usr/share/mcstas/3.0-dev/optics/Guide.comp:73]");
  stracpy(_guide_straight_2->_name, "guide_straight_2", 16384);
  stracpy(_guide_straight_2->_type, "Guide", 16384);
  _guide_straight_2->_index=17;
  _guide_straight_2->_parameters.reflect[0]='\0';
  #define reflect (_guide_straight_2->_parameters.reflect)
  _guide_straight_2->_parameters.w1 = GUIDEWIDTH;
  #define w1 (_guide_straight_2->_parameters.w1)
  _guide_straight_2->_parameters.h1 = GUIDEHEIGHT;
  #define h1 (_guide_straight_2->_parameters.h1)
  _guide_straight_2->_parameters.w2 = GUIDEWIDTH;
  #define w2 (_guide_straight_2->_parameters.w2)
  _guide_straight_2->_parameters.h2 = GUIDEHEIGHT;
  #define h2 (_guide_straight_2->_parameters.h2)
  _guide_straight_2->_parameters.l = DSTRAIGHTGUIDE2;
  #define l (_guide_straight_2->_parameters.l)
  _guide_straight_2->_parameters.R0 = instrument->_parameters._GUIDEREFLECTIVITY;
  #define R0 (_guide_straight_2->_parameters.R0)
  _guide_straight_2->_parameters.Qc = 0.0219;
  #define Qc (_guide_straight_2->_parameters.Qc)
  _guide_straight_2->_parameters.alpha = 6.07;
  #define alpha (_guide_straight_2->_parameters.alpha)
  _guide_straight_2->_parameters.m = 2.0;
  #define m (_guide_straight_2->_parameters.m)
  _guide_straight_2->_parameters.W = 0.003;
  #define W (_guide_straight_2->_parameters.W)

  #define pTable (_guide_straight_2->_parameters.pTable)

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
  /* component guide_straight_2=Guide() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _armStraight2->_rotation_absolute, _guide_straight_2->_rotation_absolute);
    rot_transpose(_armStraight2->_rotation_absolute, tr1);
    rot_mul(_guide_straight_2->_rotation_absolute, tr1, _guide_straight_2->_rotation_relative);
    _guide_straight_2->_rotation_is_identity =  rot_test_identity(_guide_straight_2->_rotation_relative);
    tc1 = coords_set(
      0, 0, 0);
    rot_transpose(_armStraight2->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _guide_straight_2->_position_absolute = coords_add(_armStraight2->_position_absolute, tc2);
    tc1 = coords_sub(_armStraight2->_position_absolute, _guide_straight_2->_position_absolute);
    _guide_straight_2->_position_relative = rot_apply(_guide_straight_2->_rotation_absolute, tc1);
  } /* guide_straight_2=Guide() AT ROTATED */
  DEBUG_COMPONENT("guide_straight_2", _guide_straight_2->_position_absolute, _guide_straight_2->_rotation_absolute);
  instrument->_position_absolute[17] = _guide_straight_2->_position_absolute;
  instrument->_position_relative[17] = _guide_straight_2->_position_relative;
  instrument->counter_N[17]  = instrument->counter_P[17] = instrument->counter_P2[17] = 0;
  instrument->counter_AbsorbProp[17]= 0;
  return(0);
} /* _guide_straight_2_setpos */

/* component armSuper=Arm() SETTING, POSITION/ROTATION */
int _armSuper_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_armSuper_setpos] component armSuper=Arm() SETTING [Arm:0]");
  stracpy(_armSuper->_name, "armSuper", 16384);
  stracpy(_armSuper->_type, "Arm", 16384);
  _armSuper->_index=18;
  /* component armSuper=Arm() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _armStraight2->_rotation_absolute, _armSuper->_rotation_absolute);
    rot_transpose(_guide_straight_2->_rotation_absolute, tr1);
    rot_mul(_armSuper->_rotation_absolute, tr1, _armSuper->_rotation_relative);
    _armSuper->_rotation_is_identity =  rot_test_identity(_armSuper->_rotation_relative);
    tc1 = coords_set(
      0, 0, DSTRAIGHTGUIDE2);
    rot_transpose(_armStraight2->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _armSuper->_position_absolute = coords_add(_armStraight2->_position_absolute, tc2);
    tc1 = coords_sub(_guide_straight_2->_position_absolute, _armSuper->_position_absolute);
    _armSuper->_position_relative = rot_apply(_armSuper->_rotation_absolute, tc1);
  } /* armSuper=Arm() AT ROTATED */
  DEBUG_COMPONENT("armSuper", _armSuper->_position_absolute, _armSuper->_rotation_absolute);
  instrument->_position_absolute[18] = _armSuper->_position_absolute;
  instrument->_position_relative[18] = _armSuper->_position_relative;
  instrument->counter_N[18]  = instrument->counter_P[18] = instrument->counter_P2[18] = 0;
  instrument->counter_AbsorbProp[18]= 0;
  return(0);
} /* _armSuper_setpos */

/* component guide_super=Guide() SETTING, POSITION/ROTATION */
int _guide_super_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_guide_super_setpos] component guide_super=Guide() SETTING [/usr/share/mcstas/3.0-dev/optics/Guide.comp:73]");
  stracpy(_guide_super->_name, "guide_super", 16384);
  stracpy(_guide_super->_type, "Guide", 16384);
  _guide_super->_index=19;
  _guide_super->_parameters.reflect[0]='\0';
  #define reflect (_guide_super->_parameters.reflect)
  _guide_super->_parameters.w1 = GUIDEWIDTH;
  #define w1 (_guide_super->_parameters.w1)
  _guide_super->_parameters.h1 = GUIDEHEIGHT;
  #define h1 (_guide_super->_parameters.h1)
  _guide_super->_parameters.w2 = GUIDEENDWIDTH;
  #define w2 (_guide_super->_parameters.w2)
  _guide_super->_parameters.h2 = GUIDEENDHEIGHT;
  #define h2 (_guide_super->_parameters.h2)
  _guide_super->_parameters.l = DSUPERGUIDE;
  #define l (_guide_super->_parameters.l)
  _guide_super->_parameters.R0 = instrument->_parameters._GUIDEREFLECTIVITY;
  #define R0 (_guide_super->_parameters.R0)
  _guide_super->_parameters.Qc = 0.0219;
  #define Qc (_guide_super->_parameters.Qc)
  _guide_super->_parameters.alpha = 6.07;
  #define alpha (_guide_super->_parameters.alpha)
  _guide_super->_parameters.m = 3.6;
  #define m (_guide_super->_parameters.m)
  _guide_super->_parameters.W = 0.003;
  #define W (_guide_super->_parameters.W)

  #define pTable (_guide_super->_parameters.pTable)

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
  /* component guide_super=Guide() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _armSuper->_rotation_absolute, _guide_super->_rotation_absolute);
    rot_transpose(_armSuper->_rotation_absolute, tr1);
    rot_mul(_guide_super->_rotation_absolute, tr1, _guide_super->_rotation_relative);
    _guide_super->_rotation_is_identity =  rot_test_identity(_guide_super->_rotation_relative);
    tc1 = coords_set(
      0, 0, 0);
    rot_transpose(_armSuper->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _guide_super->_position_absolute = coords_add(_armSuper->_position_absolute, tc2);
    tc1 = coords_sub(_armSuper->_position_absolute, _guide_super->_position_absolute);
    _guide_super->_position_relative = rot_apply(_guide_super->_rotation_absolute, tc1);
  } /* guide_super=Guide() AT ROTATED */
  DEBUG_COMPONENT("guide_super", _guide_super->_position_absolute, _guide_super->_rotation_absolute);
  instrument->_position_absolute[19] = _guide_super->_position_absolute;
  instrument->_position_relative[19] = _guide_super->_position_relative;
  instrument->counter_N[19]  = instrument->counter_P[19] = instrument->counter_P2[19] = 0;
  instrument->counter_AbsorbProp[19]= 0;
  return(0);
} /* _guide_super_setpos */

/* component armTarget=Arm() SETTING, POSITION/ROTATION */
int _armTarget_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_armTarget_setpos] component armTarget=Arm() SETTING [Arm:0]");
  stracpy(_armTarget->_name, "armTarget", 16384);
  stracpy(_armTarget->_type, "Arm", 16384);
  _armTarget->_index=20;
  /* component armTarget=Arm() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _armSuper->_rotation_absolute, _armTarget->_rotation_absolute);
    rot_transpose(_guide_super->_rotation_absolute, tr1);
    rot_mul(_armTarget->_rotation_absolute, tr1, _armTarget->_rotation_relative);
    _armTarget->_rotation_is_identity =  rot_test_identity(_armTarget->_rotation_relative);
    tc1 = coords_set(
      0, 0, DSUPERGUIDE + DGUIDETARGET);
    rot_transpose(_armSuper->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _armTarget->_position_absolute = coords_add(_armSuper->_position_absolute, tc2);
    tc1 = coords_sub(_guide_super->_position_absolute, _armTarget->_position_absolute);
    _armTarget->_position_relative = rot_apply(_armTarget->_rotation_absolute, tc1);
  } /* armTarget=Arm() AT ROTATED */
  DEBUG_COMPONENT("armTarget", _armTarget->_position_absolute, _armTarget->_rotation_absolute);
  instrument->_position_absolute[20] = _armTarget->_position_absolute;
  instrument->_position_relative[20] = _armTarget->_position_relative;
  instrument->counter_N[20]  = instrument->counter_P[20] = instrument->counter_P2[20] = 0;
  instrument->counter_AbsorbProp[20]= 0;
  return(0);
} /* _armTarget_setpos */

/* component lamTarget=L_monitor() SETTING, POSITION/ROTATION */
int _lamTarget_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_lamTarget_setpos] component lamTarget=L_monitor() SETTING [/usr/share/mcstas/3.0-dev/monitors/L_monitor.comp:65]");
  stracpy(_lamTarget->_name, "lamTarget", 16384);
  stracpy(_lamTarget->_type, "L_monitor", 16384);
  _lamTarget->_index=21;
  _lamTarget->_parameters.nL = 120;
  #define nL (_lamTarget->_parameters.nL)
  if("lambdaTarget.dat" && strlen("lambdaTarget.dat"))
    stracpy(_lamTarget->_parameters.filename, "lambdaTarget.dat" ? "lambdaTarget.dat" : "", 16384);
  else 
  _lamTarget->_parameters.filename[0]='\0';
  #define filename (_lamTarget->_parameters.filename)
  _lamTarget->_parameters.xmin = -0.05;
  #define xmin (_lamTarget->_parameters.xmin)
  _lamTarget->_parameters.xmax = 0.05;
  #define xmax (_lamTarget->_parameters.xmax)
  _lamTarget->_parameters.ymin = -0.05;
  #define ymin (_lamTarget->_parameters.ymin)
  _lamTarget->_parameters.ymax = 0.05;
  #define ymax (_lamTarget->_parameters.ymax)
  _lamTarget->_parameters.xwidth = 0.10;
  #define xwidth (_lamTarget->_parameters.xwidth)
  _lamTarget->_parameters.yheight = 0.10;
  #define yheight (_lamTarget->_parameters.yheight)
  _lamTarget->_parameters.Lmin = 0.98 * LAMBDAMIN;
  #define Lmin (_lamTarget->_parameters.Lmin)
  _lamTarget->_parameters.Lmax = 1.02 * LAMBDAMAX;
  #define Lmax (_lamTarget->_parameters.Lmax)
  _lamTarget->_parameters.restore_neutron = 0;
  #define restore_neutron (_lamTarget->_parameters.restore_neutron)

  #define L_N (_lamTarget->_parameters.L_N)
  #define L_p (_lamTarget->_parameters.L_p)
  #define L_p2 (_lamTarget->_parameters.L_p2)

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
  /* component lamTarget=L_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _armTarget->_rotation_absolute, _lamTarget->_rotation_absolute);
    rot_transpose(_armTarget->_rotation_absolute, tr1);
    rot_mul(_lamTarget->_rotation_absolute, tr1, _lamTarget->_rotation_relative);
    _lamTarget->_rotation_is_identity =  rot_test_identity(_lamTarget->_rotation_relative);
    tc1 = coords_set(
      0, 0, 0);
    rot_transpose(_armTarget->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _lamTarget->_position_absolute = coords_add(_armTarget->_position_absolute, tc2);
    tc1 = coords_sub(_armTarget->_position_absolute, _lamTarget->_position_absolute);
    _lamTarget->_position_relative = rot_apply(_lamTarget->_rotation_absolute, tc1);
  } /* lamTarget=L_monitor() AT ROTATED */
  DEBUG_COMPONENT("lamTarget", _lamTarget->_position_absolute, _lamTarget->_rotation_absolute);
  instrument->_position_absolute[21] = _lamTarget->_position_absolute;
  instrument->_position_relative[21] = _lamTarget->_position_relative;
  instrument->counter_N[21]  = instrument->counter_P[21] = instrument->counter_P2[21] = 0;
  instrument->counter_AbsorbProp[21]= 0;
  return(0);
} /* _lamTarget_setpos */

/* component psdTarget=PSD_monitor() SETTING, POSITION/ROTATION */
int _psdTarget_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_psdTarget_setpos] component psdTarget=PSD_monitor() SETTING [/usr/share/mcstas/3.0-dev/monitors/PSD_monitor.comp:64]");
  stracpy(_psdTarget->_name, "psdTarget", 16384);
  stracpy(_psdTarget->_type, "PSD_monitor", 16384);
  _psdTarget->_index=22;
  _psdTarget->_parameters.nx = 100;
  #define nx (_psdTarget->_parameters.nx)
  _psdTarget->_parameters.ny = 100;
  #define ny (_psdTarget->_parameters.ny)
  if("psdTarget.dat" && strlen("psdTarget.dat"))
    stracpy(_psdTarget->_parameters.filename, "psdTarget.dat" ? "psdTarget.dat" : "", 16384);
  else 
  _psdTarget->_parameters.filename[0]='\0';
  #define filename (_psdTarget->_parameters.filename)
  _psdTarget->_parameters.xmin = -0.05;
  #define xmin (_psdTarget->_parameters.xmin)
  _psdTarget->_parameters.xmax = 0.05;
  #define xmax (_psdTarget->_parameters.xmax)
  _psdTarget->_parameters.ymin = -0.05;
  #define ymin (_psdTarget->_parameters.ymin)
  _psdTarget->_parameters.ymax = 0.05;
  #define ymax (_psdTarget->_parameters.ymax)
  _psdTarget->_parameters.xwidth = 0.10;
  #define xwidth (_psdTarget->_parameters.xwidth)
  _psdTarget->_parameters.yheight = 0.10;
  #define yheight (_psdTarget->_parameters.yheight)
  _psdTarget->_parameters.restore_neutron = 0;
  #define restore_neutron (_psdTarget->_parameters.restore_neutron)

  #define PSD_N (_psdTarget->_parameters.PSD_N)
  #define PSD_p (_psdTarget->_parameters.PSD_p)
  #define PSD_p2 (_psdTarget->_parameters.PSD_p2)

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
  /* component psdTarget=PSD_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _armTarget->_rotation_absolute, _psdTarget->_rotation_absolute);
    rot_transpose(_lamTarget->_rotation_absolute, tr1);
    rot_mul(_psdTarget->_rotation_absolute, tr1, _psdTarget->_rotation_relative);
    _psdTarget->_rotation_is_identity =  rot_test_identity(_psdTarget->_rotation_relative);
    tc1 = coords_set(
      0, 0, 0);
    rot_transpose(_armTarget->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _psdTarget->_position_absolute = coords_add(_armTarget->_position_absolute, tc2);
    tc1 = coords_sub(_lamTarget->_position_absolute, _psdTarget->_position_absolute);
    _psdTarget->_position_relative = rot_apply(_psdTarget->_rotation_absolute, tc1);
  } /* psdTarget=PSD_monitor() AT ROTATED */
  DEBUG_COMPONENT("psdTarget", _psdTarget->_position_absolute, _psdTarget->_rotation_absolute);
  instrument->_position_absolute[22] = _psdTarget->_position_absolute;
  instrument->_position_relative[22] = _psdTarget->_position_relative;
  instrument->counter_N[22]  = instrument->counter_P[22] = instrument->counter_P2[22] = 0;
  instrument->counter_AbsorbProp[22]= 0;
  return(0);
} /* _psdTarget_setpos */

/* component vsample1=Incoherent() SETTING, POSITION/ROTATION */
int _vsample1_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_vsample1_setpos] component vsample1=Incoherent() SETTING [/usr/share/mcstas/3.0-dev/samples/Incoherent.comp:153]");
  stracpy(_vsample1->_name, "vsample1", 16384);
  stracpy(_vsample1->_type, "Incoherent", 16384);
  _vsample1->_index=23;
  _vsample1->_parameters.geometry[0]='\0';
  #define geometry (_vsample1->_parameters.geometry)
  _vsample1->_parameters.radius = 0.02;
  #define radius (_vsample1->_parameters.radius)
  _vsample1->_parameters.xwidth = 0;
  #define xwidth (_vsample1->_parameters.xwidth)
  _vsample1->_parameters.yheight = VANADIUMHEIGHT;
  #define yheight (_vsample1->_parameters.yheight)
  _vsample1->_parameters.zdepth = 0;
  #define zdepth (_vsample1->_parameters.zdepth)
  _vsample1->_parameters.thickness = 0.00004;
  #define thickness (_vsample1->_parameters.thickness)
  _vsample1->_parameters.target_x = - ANARADIUS;
  #define target_x (_vsample1->_parameters.target_x)
  _vsample1->_parameters.target_y = 0;
  #define target_y (_vsample1->_parameters.target_y)
  _vsample1->_parameters.target_z = 0;
  #define target_z (_vsample1->_parameters.target_z)
  _vsample1->_parameters.focus_r = 0;
  #define focus_r (_vsample1->_parameters.focus_r)
  _vsample1->_parameters.focus_xw = 0;
  #define focus_xw (_vsample1->_parameters.focus_xw)
  _vsample1->_parameters.focus_yh = 0;
  #define focus_yh (_vsample1->_parameters.focus_yh)
  _vsample1->_parameters.focus_aw = 0.7;
  #define focus_aw (_vsample1->_parameters.focus_aw)
  _vsample1->_parameters.focus_ah = 30.0;
  #define focus_ah (_vsample1->_parameters.focus_ah)
  _vsample1->_parameters.target_index = 0;
  #define target_index (_vsample1->_parameters.target_index)
  _vsample1->_parameters.pack = 1;
  #define pack (_vsample1->_parameters.pack)
  _vsample1->_parameters.p_interact = 1;
  #define p_interact (_vsample1->_parameters.p_interact)
  _vsample1->_parameters.f_QE = 0;
  #define f_QE (_vsample1->_parameters.f_QE)
  _vsample1->_parameters.gamma = 0;
  #define gamma (_vsample1->_parameters.gamma)
  _vsample1->_parameters.sigma_abs = 5.08;
  #define sigma_abs (_vsample1->_parameters.sigma_abs)
  _vsample1->_parameters.sigma_inc = 5.08;
  #define sigma_inc (_vsample1->_parameters.sigma_inc)
  _vsample1->_parameters.Vc = 13.827;
  #define Vc (_vsample1->_parameters.Vc)
  _vsample1->_parameters.concentric = 0;
  #define concentric (_vsample1->_parameters.concentric)
  _vsample1->_parameters.order = 0;
  #define order (_vsample1->_parameters.order)

  #define VarsInc (_vsample1->_parameters.VarsInc)
  #define offdata (_vsample1->_parameters.offdata)

  #undef geometry
  #undef radius
  #undef xwidth
  #undef yheight
  #undef zdepth
  #undef thickness
  #undef target_x
  #undef target_y
  #undef target_z
  #undef focus_r
  #undef focus_xw
  #undef focus_yh
  #undef focus_aw
  #undef focus_ah
  #undef target_index
  #undef pack
  #undef p_interact
  #undef f_QE
  #undef gamma
  #undef sigma_abs
  #undef sigma_inc
  #undef Vc
  #undef concentric
  #undef order
  #undef VarsInc
  #undef offdata
  /* component vsample1=Incoherent() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _armTarget->_rotation_absolute, _vsample1->_rotation_absolute);
    rot_transpose(_psdTarget->_rotation_absolute, tr1);
    rot_mul(_vsample1->_rotation_absolute, tr1, _vsample1->_rotation_relative);
    _vsample1->_rotation_is_identity =  rot_test_identity(_vsample1->_rotation_relative);
    tc1 = coords_set(
      0, 0, 0);
    rot_transpose(_armTarget->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _vsample1->_position_absolute = coords_add(_armTarget->_position_absolute, tc2);
    tc1 = coords_sub(_psdTarget->_position_absolute, _vsample1->_position_absolute);
    _vsample1->_position_relative = rot_apply(_vsample1->_rotation_absolute, tc1);
  } /* vsample1=Incoherent() AT ROTATED */
  DEBUG_COMPONENT("vsample1", _vsample1->_position_absolute, _vsample1->_rotation_absolute);
  instrument->_position_absolute[23] = _vsample1->_position_absolute;
  instrument->_position_relative[23] = _vsample1->_position_relative;
  instrument->counter_N[23]  = instrument->counter_P[23] = instrument->counter_P2[23] = 0;
  instrument->counter_AbsorbProp[23]= 0;
  return(0);
} /* _vsample1_setpos */

/* component vsample2=Incoherent() SETTING, POSITION/ROTATION */
int _vsample2_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_vsample2_setpos] component vsample2=Incoherent() SETTING [/usr/share/mcstas/3.0-dev/samples/Incoherent.comp:153]");
  stracpy(_vsample2->_name, "vsample2", 16384);
  stracpy(_vsample2->_type, "Incoherent", 16384);
  _vsample2->_index=24;
  _vsample2->_parameters.geometry[0]='\0';
  #define geometry (_vsample2->_parameters.geometry)
  _vsample2->_parameters.radius = 0.016;
  #define radius (_vsample2->_parameters.radius)
  _vsample2->_parameters.xwidth = 0;
  #define xwidth (_vsample2->_parameters.xwidth)
  _vsample2->_parameters.yheight = VANADIUMHEIGHT;
  #define yheight (_vsample2->_parameters.yheight)
  _vsample2->_parameters.zdepth = 0;
  #define zdepth (_vsample2->_parameters.zdepth)
  _vsample2->_parameters.thickness = 0.00004;
  #define thickness (_vsample2->_parameters.thickness)
  _vsample2->_parameters.target_x = - ANARADIUS;
  #define target_x (_vsample2->_parameters.target_x)
  _vsample2->_parameters.target_y = 0;
  #define target_y (_vsample2->_parameters.target_y)
  _vsample2->_parameters.target_z = 0;
  #define target_z (_vsample2->_parameters.target_z)
  _vsample2->_parameters.focus_r = 0;
  #define focus_r (_vsample2->_parameters.focus_r)
  _vsample2->_parameters.focus_xw = 0;
  #define focus_xw (_vsample2->_parameters.focus_xw)
  _vsample2->_parameters.focus_yh = 0;
  #define focus_yh (_vsample2->_parameters.focus_yh)
  _vsample2->_parameters.focus_aw = 0.7;
  #define focus_aw (_vsample2->_parameters.focus_aw)
  _vsample2->_parameters.focus_ah = 30.0;
  #define focus_ah (_vsample2->_parameters.focus_ah)
  _vsample2->_parameters.target_index = 0;
  #define target_index (_vsample2->_parameters.target_index)
  _vsample2->_parameters.pack = 1;
  #define pack (_vsample2->_parameters.pack)
  _vsample2->_parameters.p_interact = 1;
  #define p_interact (_vsample2->_parameters.p_interact)
  _vsample2->_parameters.f_QE = 0;
  #define f_QE (_vsample2->_parameters.f_QE)
  _vsample2->_parameters.gamma = 0;
  #define gamma (_vsample2->_parameters.gamma)
  _vsample2->_parameters.sigma_abs = 5.08;
  #define sigma_abs (_vsample2->_parameters.sigma_abs)
  _vsample2->_parameters.sigma_inc = 5.08;
  #define sigma_inc (_vsample2->_parameters.sigma_inc)
  _vsample2->_parameters.Vc = 13.827;
  #define Vc (_vsample2->_parameters.Vc)
  _vsample2->_parameters.concentric = 0;
  #define concentric (_vsample2->_parameters.concentric)
  _vsample2->_parameters.order = 0;
  #define order (_vsample2->_parameters.order)

  #define VarsInc (_vsample2->_parameters.VarsInc)
  #define offdata (_vsample2->_parameters.offdata)

  #undef geometry
  #undef radius
  #undef xwidth
  #undef yheight
  #undef zdepth
  #undef thickness
  #undef target_x
  #undef target_y
  #undef target_z
  #undef focus_r
  #undef focus_xw
  #undef focus_yh
  #undef focus_aw
  #undef focus_ah
  #undef target_index
  #undef pack
  #undef p_interact
  #undef f_QE
  #undef gamma
  #undef sigma_abs
  #undef sigma_inc
  #undef Vc
  #undef concentric
  #undef order
  #undef VarsInc
  #undef offdata
  /* component vsample2=Incoherent() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _armTarget->_rotation_absolute, _vsample2->_rotation_absolute);
    rot_transpose(_vsample1->_rotation_absolute, tr1);
    rot_mul(_vsample2->_rotation_absolute, tr1, _vsample2->_rotation_relative);
    _vsample2->_rotation_is_identity =  rot_test_identity(_vsample2->_rotation_relative);
    tc1 = coords_set(
      0, 0, 0);
    rot_transpose(_armTarget->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _vsample2->_position_absolute = coords_add(_armTarget->_position_absolute, tc2);
    tc1 = coords_sub(_vsample1->_position_absolute, _vsample2->_position_absolute);
    _vsample2->_position_relative = rot_apply(_vsample2->_rotation_absolute, tc1);
  } /* vsample2=Incoherent() AT ROTATED */
  DEBUG_COMPONENT("vsample2", _vsample2->_position_absolute, _vsample2->_rotation_absolute);
  instrument->_position_absolute[24] = _vsample2->_position_absolute;
  instrument->_position_relative[24] = _vsample2->_position_relative;
  instrument->counter_N[24]  = instrument->counter_P[24] = instrument->counter_P2[24] = 0;
  instrument->counter_AbsorbProp[24]= 0;
  return(0);
} /* _vsample2_setpos */

/* component vsample3=Incoherent() SETTING, POSITION/ROTATION */
int _vsample3_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_vsample3_setpos] component vsample3=Incoherent() SETTING [/usr/share/mcstas/3.0-dev/samples/Incoherent.comp:153]");
  stracpy(_vsample3->_name, "vsample3", 16384);
  stracpy(_vsample3->_type, "Incoherent", 16384);
  _vsample3->_index=25;
  _vsample3->_parameters.geometry[0]='\0';
  #define geometry (_vsample3->_parameters.geometry)
  _vsample3->_parameters.radius = 0.01;
  #define radius (_vsample3->_parameters.radius)
  _vsample3->_parameters.xwidth = 0;
  #define xwidth (_vsample3->_parameters.xwidth)
  _vsample3->_parameters.yheight = VANADIUMHEIGHT;
  #define yheight (_vsample3->_parameters.yheight)
  _vsample3->_parameters.zdepth = 0;
  #define zdepth (_vsample3->_parameters.zdepth)
  _vsample3->_parameters.thickness = 0.00004;
  #define thickness (_vsample3->_parameters.thickness)
  _vsample3->_parameters.target_x = - ANARADIUS;
  #define target_x (_vsample3->_parameters.target_x)
  _vsample3->_parameters.target_y = 0;
  #define target_y (_vsample3->_parameters.target_y)
  _vsample3->_parameters.target_z = 0;
  #define target_z (_vsample3->_parameters.target_z)
  _vsample3->_parameters.focus_r = 0;
  #define focus_r (_vsample3->_parameters.focus_r)
  _vsample3->_parameters.focus_xw = 0;
  #define focus_xw (_vsample3->_parameters.focus_xw)
  _vsample3->_parameters.focus_yh = 0;
  #define focus_yh (_vsample3->_parameters.focus_yh)
  _vsample3->_parameters.focus_aw = 0.7;
  #define focus_aw (_vsample3->_parameters.focus_aw)
  _vsample3->_parameters.focus_ah = 30.0;
  #define focus_ah (_vsample3->_parameters.focus_ah)
  _vsample3->_parameters.target_index = 0;
  #define target_index (_vsample3->_parameters.target_index)
  _vsample3->_parameters.pack = 1;
  #define pack (_vsample3->_parameters.pack)
  _vsample3->_parameters.p_interact = 1;
  #define p_interact (_vsample3->_parameters.p_interact)
  _vsample3->_parameters.f_QE = 0;
  #define f_QE (_vsample3->_parameters.f_QE)
  _vsample3->_parameters.gamma = 0;
  #define gamma (_vsample3->_parameters.gamma)
  _vsample3->_parameters.sigma_abs = 5.08;
  #define sigma_abs (_vsample3->_parameters.sigma_abs)
  _vsample3->_parameters.sigma_inc = 5.08;
  #define sigma_inc (_vsample3->_parameters.sigma_inc)
  _vsample3->_parameters.Vc = 13.827;
  #define Vc (_vsample3->_parameters.Vc)
  _vsample3->_parameters.concentric = 0;
  #define concentric (_vsample3->_parameters.concentric)
  _vsample3->_parameters.order = 0;
  #define order (_vsample3->_parameters.order)

  #define VarsInc (_vsample3->_parameters.VarsInc)
  #define offdata (_vsample3->_parameters.offdata)

  #undef geometry
  #undef radius
  #undef xwidth
  #undef yheight
  #undef zdepth
  #undef thickness
  #undef target_x
  #undef target_y
  #undef target_z
  #undef focus_r
  #undef focus_xw
  #undef focus_yh
  #undef focus_aw
  #undef focus_ah
  #undef target_index
  #undef pack
  #undef p_interact
  #undef f_QE
  #undef gamma
  #undef sigma_abs
  #undef sigma_inc
  #undef Vc
  #undef concentric
  #undef order
  #undef VarsInc
  #undef offdata
  /* component vsample3=Incoherent() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _armTarget->_rotation_absolute, _vsample3->_rotation_absolute);
    rot_transpose(_vsample2->_rotation_absolute, tr1);
    rot_mul(_vsample3->_rotation_absolute, tr1, _vsample3->_rotation_relative);
    _vsample3->_rotation_is_identity =  rot_test_identity(_vsample3->_rotation_relative);
    tc1 = coords_set(
      0, 0, 0);
    rot_transpose(_armTarget->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _vsample3->_position_absolute = coords_add(_armTarget->_position_absolute, tc2);
    tc1 = coords_sub(_vsample2->_position_absolute, _vsample3->_position_absolute);
    _vsample3->_position_relative = rot_apply(_vsample3->_rotation_absolute, tc1);
  } /* vsample3=Incoherent() AT ROTATED */
  DEBUG_COMPONENT("vsample3", _vsample3->_position_absolute, _vsample3->_rotation_absolute);
  instrument->_position_absolute[25] = _vsample3->_position_absolute;
  instrument->_position_relative[25] = _vsample3->_position_relative;
  instrument->counter_N[25]  = instrument->counter_P[25] = instrument->counter_P2[25] = 0;
  instrument->counter_AbsorbProp[25]= 0;
  return(0);
} /* _vsample3_setpos */

/* component beamStop=Beamstop() SETTING, POSITION/ROTATION */
int _beamStop_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_beamStop_setpos] component beamStop=Beamstop() SETTING [/usr/share/mcstas/3.0-dev/optics/Beamstop.comp:50]");
  stracpy(_beamStop->_name, "beamStop", 16384);
  stracpy(_beamStop->_type, "Beamstop", 16384);
  _beamStop->_index=26;
  _beamStop->_parameters.xmin = -0.05;
  #define xmin (_beamStop->_parameters.xmin)
  _beamStop->_parameters.xmax = 0.05;
  #define xmax (_beamStop->_parameters.xmax)
  _beamStop->_parameters.ymin = -0.05;
  #define ymin (_beamStop->_parameters.ymin)
  _beamStop->_parameters.ymax = 0.05;
  #define ymax (_beamStop->_parameters.ymax)
  _beamStop->_parameters.xwidth = 0;
  #define xwidth (_beamStop->_parameters.xwidth)
  _beamStop->_parameters.yheight = 0;
  #define yheight (_beamStop->_parameters.yheight)
  _beamStop->_parameters.radius = 0.1;
  #define radius (_beamStop->_parameters.radius)

  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef xwidth
  #undef yheight
  #undef radius
  /* component beamStop=Beamstop() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _armTarget->_rotation_absolute, _beamStop->_rotation_absolute);
    rot_transpose(_vsample3->_rotation_absolute, tr1);
    rot_mul(_beamStop->_rotation_absolute, tr1, _beamStop->_rotation_relative);
    _beamStop->_rotation_is_identity =  rot_test_identity(_beamStop->_rotation_relative);
    tc1 = coords_set(
      0, 0, 0.1);
    rot_transpose(_armTarget->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _beamStop->_position_absolute = coords_add(_armTarget->_position_absolute, tc2);
    tc1 = coords_sub(_vsample3->_position_absolute, _beamStop->_position_absolute);
    _beamStop->_position_relative = rot_apply(_beamStop->_rotation_absolute, tc1);
  } /* beamStop=Beamstop() AT ROTATED */
  DEBUG_COMPONENT("beamStop", _beamStop->_position_absolute, _beamStop->_rotation_absolute);
  instrument->_position_absolute[26] = _beamStop->_position_absolute;
  instrument->_position_relative[26] = _beamStop->_position_relative;
  instrument->counter_N[26]  = instrument->counter_P[26] = instrument->counter_P2[26] = 0;
  instrument->counter_AbsorbProp[26]= 0;
  return(0);
} /* _beamStop_setpos */

/* component armAnalyzer=Arm() SETTING, POSITION/ROTATION */
int _armAnalyzer_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_armAnalyzer_setpos] component armAnalyzer=Arm() SETTING [Arm:0]");
  stracpy(_armAnalyzer->_name, "armAnalyzer", 16384);
  stracpy(_armAnalyzer->_type, "Arm", 16384);
  _armAnalyzer->_index=27;
  /* component armAnalyzer=Arm() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (-90)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _armTarget->_rotation_absolute, _armAnalyzer->_rotation_absolute);
    rot_transpose(_beamStop->_rotation_absolute, tr1);
    rot_mul(_armAnalyzer->_rotation_absolute, tr1, _armAnalyzer->_rotation_relative);
    _armAnalyzer->_rotation_is_identity =  rot_test_identity(_armAnalyzer->_rotation_relative);
    tc1 = coords_set(
      0, 0, 0);
    rot_transpose(_armTarget->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _armAnalyzer->_position_absolute = coords_add(_armTarget->_position_absolute, tc2);
    tc1 = coords_sub(_beamStop->_position_absolute, _armAnalyzer->_position_absolute);
    _armAnalyzer->_position_relative = rot_apply(_armAnalyzer->_rotation_absolute, tc1);
  } /* armAnalyzer=Arm() AT ROTATED */
  DEBUG_COMPONENT("armAnalyzer", _armAnalyzer->_position_absolute, _armAnalyzer->_rotation_absolute);
  instrument->_position_absolute[27] = _armAnalyzer->_position_absolute;
  instrument->_position_relative[27] = _armAnalyzer->_position_relative;
  instrument->counter_N[27]  = instrument->counter_P[27] = instrument->counter_P2[27] = 0;
  instrument->counter_AbsorbProp[27]= 0;
  return(0);
} /* _armAnalyzer_setpos */

/* component graphite_analyser0=Monochromator_pol() SETTING, POSITION/ROTATION */
int _graphite_analyser0_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_graphite_analyser0_setpos] component graphite_analyser0=Monochromator_pol() SETTING [/usr/share/mcstas/3.0-dev/optics/Monochromator_pol.comp:109]");
  stracpy(_graphite_analyser0->_name, "graphite_analyser0", 16384);
  stracpy(_graphite_analyser0->_type, "Monochromator_pol", 16384);
  _graphite_analyser0->_index=28;
  _graphite_analyser0->_parameters.zwidth = ANAWIDTH;
  #define zwidth (_graphite_analyser0->_parameters.zwidth)
  _graphite_analyser0->_parameters.yheight = ANAHEIGHT;
  #define yheight (_graphite_analyser0->_parameters.yheight)
  _graphite_analyser0->_parameters.mosaic = ANAMOSAIC;
  #define mosaic (_graphite_analyser0->_parameters.mosaic)
  _graphite_analyser0->_parameters.dspread = 0.0025;
  #define dspread (_graphite_analyser0->_parameters.dspread)
  _graphite_analyser0->_parameters.Q = ANAQ;
  #define Q (_graphite_analyser0->_parameters.Q)
  _graphite_analyser0->_parameters.DM = 0;
  #define DM (_graphite_analyser0->_parameters.DM)
  _graphite_analyser0->_parameters.pThreshold = 0.001;
  #define pThreshold (_graphite_analyser0->_parameters.pThreshold)
  _graphite_analyser0->_parameters.Rup = 1.0;
  #define Rup (_graphite_analyser0->_parameters.Rup)
  _graphite_analyser0->_parameters.Rdown = 1.0;
  #define Rdown (_graphite_analyser0->_parameters.Rdown)
  _graphite_analyser0->_parameters.debug = 0;
  #define debug (_graphite_analyser0->_parameters.debug)

  #define mos_rms (_graphite_analyser0->_parameters.mos_rms)
  #define d_rms (_graphite_analyser0->_parameters.d_rms)
  #define mono_Q (_graphite_analyser0->_parameters.mono_Q)
  #define FN (_graphite_analyser0->_parameters.FN)
  #define FM (_graphite_analyser0->_parameters.FM)

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
  /* component graphite_analyser0=Monochromator_pol() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (90)*DEG2RAD, (8.928)*DEG2RAD);
    rot_mul(tr1, _armAnalyzer->_rotation_absolute, _graphite_analyser0->_rotation_absolute);
    rot_transpose(_armAnalyzer->_rotation_absolute, tr1);
    rot_mul(_graphite_analyser0->_rotation_absolute, tr1, _graphite_analyser0->_rotation_relative);
    _graphite_analyser0->_rotation_is_identity =  rot_test_identity(_graphite_analyser0->_rotation_relative);
    tc1 = coords_set(
      0, -0.195, 0.899629);
    rot_transpose(_armAnalyzer->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _graphite_analyser0->_position_absolute = coords_add(_armAnalyzer->_position_absolute, tc2);
    tc1 = coords_sub(_armAnalyzer->_position_absolute, _graphite_analyser0->_position_absolute);
    _graphite_analyser0->_position_relative = rot_apply(_graphite_analyser0->_rotation_absolute, tc1);
  } /* graphite_analyser0=Monochromator_pol() AT ROTATED */
  DEBUG_COMPONENT("graphite_analyser0", _graphite_analyser0->_position_absolute, _graphite_analyser0->_rotation_absolute);
  instrument->_position_absolute[28] = _graphite_analyser0->_position_absolute;
  instrument->_position_relative[28] = _graphite_analyser0->_position_relative;
  instrument->counter_N[28]  = instrument->counter_P[28] = instrument->counter_P2[28] = 0;
  instrument->counter_AbsorbProp[28]= 0;
  return(0);
} /* _graphite_analyser0_setpos */

/* component graphite_analyser1=Monochromator_pol() SETTING, POSITION/ROTATION */
int _graphite_analyser1_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_graphite_analyser1_setpos] component graphite_analyser1=Monochromator_pol() SETTING [/usr/share/mcstas/3.0-dev/optics/Monochromator_pol.comp:109]");
  stracpy(_graphite_analyser1->_name, "graphite_analyser1", 16384);
  stracpy(_graphite_analyser1->_type, "Monochromator_pol", 16384);
  _graphite_analyser1->_index=29;
  _graphite_analyser1->_parameters.zwidth = ANAWIDTH;
  #define zwidth (_graphite_analyser1->_parameters.zwidth)
  _graphite_analyser1->_parameters.yheight = ANAHEIGHT;
  #define yheight (_graphite_analyser1->_parameters.yheight)
  _graphite_analyser1->_parameters.mosaic = ANAMOSAIC;
  #define mosaic (_graphite_analyser1->_parameters.mosaic)
  _graphite_analyser1->_parameters.dspread = 0.0025;
  #define dspread (_graphite_analyser1->_parameters.dspread)
  _graphite_analyser1->_parameters.Q = ANAQ;
  #define Q (_graphite_analyser1->_parameters.Q)
  _graphite_analyser1->_parameters.DM = 0;
  #define DM (_graphite_analyser1->_parameters.DM)
  _graphite_analyser1->_parameters.pThreshold = 0.001;
  #define pThreshold (_graphite_analyser1->_parameters.pThreshold)
  _graphite_analyser1->_parameters.Rup = 1.0;
  #define Rup (_graphite_analyser1->_parameters.Rup)
  _graphite_analyser1->_parameters.Rdown = 1.0;
  #define Rdown (_graphite_analyser1->_parameters.Rdown)
  _graphite_analyser1->_parameters.debug = 0;
  #define debug (_graphite_analyser1->_parameters.debug)

  #define mos_rms (_graphite_analyser1->_parameters.mos_rms)
  #define d_rms (_graphite_analyser1->_parameters.d_rms)
  #define mono_Q (_graphite_analyser1->_parameters.mono_Q)
  #define FN (_graphite_analyser1->_parameters.FN)
  #define FM (_graphite_analyser1->_parameters.FM)

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
  /* component graphite_analyser1=Monochromator_pol() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (90)*DEG2RAD, (8.928)*DEG2RAD);
    rot_mul(tr1, _armAnalyzer->_rotation_absolute, _graphite_analyser1->_rotation_absolute);
    rot_transpose(_graphite_analyser0->_rotation_absolute, tr1);
    rot_mul(_graphite_analyser1->_rotation_absolute, tr1, _graphite_analyser1->_rotation_relative);
    _graphite_analyser1->_rotation_is_identity =  rot_test_identity(_graphite_analyser1->_rotation_relative);
    tc1 = coords_set(
      0, -0.185, 0.9012);
    rot_transpose(_armAnalyzer->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _graphite_analyser1->_position_absolute = coords_add(_armAnalyzer->_position_absolute, tc2);
    tc1 = coords_sub(_graphite_analyser0->_position_absolute, _graphite_analyser1->_position_absolute);
    _graphite_analyser1->_position_relative = rot_apply(_graphite_analyser1->_rotation_absolute, tc1);
  } /* graphite_analyser1=Monochromator_pol() AT ROTATED */
  DEBUG_COMPONENT("graphite_analyser1", _graphite_analyser1->_position_absolute, _graphite_analyser1->_rotation_absolute);
  instrument->_position_absolute[29] = _graphite_analyser1->_position_absolute;
  instrument->_position_relative[29] = _graphite_analyser1->_position_relative;
  instrument->counter_N[29]  = instrument->counter_P[29] = instrument->counter_P2[29] = 0;
  instrument->counter_AbsorbProp[29]= 0;
  return(0);
} /* _graphite_analyser1_setpos */

/* component graphite_analyser2=Monochromator_pol() SETTING, POSITION/ROTATION */
int _graphite_analyser2_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_graphite_analyser2_setpos] component graphite_analyser2=Monochromator_pol() SETTING [/usr/share/mcstas/3.0-dev/optics/Monochromator_pol.comp:109]");
  stracpy(_graphite_analyser2->_name, "graphite_analyser2", 16384);
  stracpy(_graphite_analyser2->_type, "Monochromator_pol", 16384);
  _graphite_analyser2->_index=30;
  _graphite_analyser2->_parameters.zwidth = ANAWIDTH;
  #define zwidth (_graphite_analyser2->_parameters.zwidth)
  _graphite_analyser2->_parameters.yheight = ANAHEIGHT;
  #define yheight (_graphite_analyser2->_parameters.yheight)
  _graphite_analyser2->_parameters.mosaic = ANAMOSAIC;
  #define mosaic (_graphite_analyser2->_parameters.mosaic)
  _graphite_analyser2->_parameters.dspread = 0.0025;
  #define dspread (_graphite_analyser2->_parameters.dspread)
  _graphite_analyser2->_parameters.Q = ANAQ;
  #define Q (_graphite_analyser2->_parameters.Q)
  _graphite_analyser2->_parameters.DM = 0;
  #define DM (_graphite_analyser2->_parameters.DM)
  _graphite_analyser2->_parameters.pThreshold = 0.001;
  #define pThreshold (_graphite_analyser2->_parameters.pThreshold)
  _graphite_analyser2->_parameters.Rup = 1.0;
  #define Rup (_graphite_analyser2->_parameters.Rup)
  _graphite_analyser2->_parameters.Rdown = 1.0;
  #define Rdown (_graphite_analyser2->_parameters.Rdown)
  _graphite_analyser2->_parameters.debug = 0;
  #define debug (_graphite_analyser2->_parameters.debug)

  #define mos_rms (_graphite_analyser2->_parameters.mos_rms)
  #define d_rms (_graphite_analyser2->_parameters.d_rms)
  #define mono_Q (_graphite_analyser2->_parameters.mono_Q)
  #define FN (_graphite_analyser2->_parameters.FN)
  #define FM (_graphite_analyser2->_parameters.FM)

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
  /* component graphite_analyser2=Monochromator_pol() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (90)*DEG2RAD, (8.1305)*DEG2RAD);
    rot_mul(tr1, _armAnalyzer->_rotation_absolute, _graphite_analyser2->_rotation_absolute);
    rot_transpose(_graphite_analyser1->_rotation_absolute, tr1);
    rot_mul(_graphite_analyser2->_rotation_absolute, tr1, _graphite_analyser2->_rotation_relative);
    _graphite_analyser2->_rotation_is_identity =  rot_test_identity(_graphite_analyser2->_rotation_relative);
    tc1 = coords_set(
      0, -0.175, 0.902637);
    rot_transpose(_armAnalyzer->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _graphite_analyser2->_position_absolute = coords_add(_armAnalyzer->_position_absolute, tc2);
    tc1 = coords_sub(_graphite_analyser1->_position_absolute, _graphite_analyser2->_position_absolute);
    _graphite_analyser2->_position_relative = rot_apply(_graphite_analyser2->_rotation_absolute, tc1);
  } /* graphite_analyser2=Monochromator_pol() AT ROTATED */
  DEBUG_COMPONENT("graphite_analyser2", _graphite_analyser2->_position_absolute, _graphite_analyser2->_rotation_absolute);
  instrument->_position_absolute[30] = _graphite_analyser2->_position_absolute;
  instrument->_position_relative[30] = _graphite_analyser2->_position_relative;
  instrument->counter_N[30]  = instrument->counter_P[30] = instrument->counter_P2[30] = 0;
  instrument->counter_AbsorbProp[30]= 0;
  return(0);
} /* _graphite_analyser2_setpos */

/* component graphite_analyser3=Monochromator_pol() SETTING, POSITION/ROTATION */
int _graphite_analyser3_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_graphite_analyser3_setpos] component graphite_analyser3=Monochromator_pol() SETTING [/usr/share/mcstas/3.0-dev/optics/Monochromator_pol.comp:109]");
  stracpy(_graphite_analyser3->_name, "graphite_analyser3", 16384);
  stracpy(_graphite_analyser3->_type, "Monochromator_pol", 16384);
  _graphite_analyser3->_index=31;
  _graphite_analyser3->_parameters.zwidth = ANAWIDTH;
  #define zwidth (_graphite_analyser3->_parameters.zwidth)
  _graphite_analyser3->_parameters.yheight = ANAHEIGHT;
  #define yheight (_graphite_analyser3->_parameters.yheight)
  _graphite_analyser3->_parameters.mosaic = ANAMOSAIC;
  #define mosaic (_graphite_analyser3->_parameters.mosaic)
  _graphite_analyser3->_parameters.dspread = 0.0025;
  #define dspread (_graphite_analyser3->_parameters.dspread)
  _graphite_analyser3->_parameters.Q = ANAQ;
  #define Q (_graphite_analyser3->_parameters.Q)
  _graphite_analyser3->_parameters.DM = 0;
  #define DM (_graphite_analyser3->_parameters.DM)
  _graphite_analyser3->_parameters.pThreshold = 0.001;
  #define pThreshold (_graphite_analyser3->_parameters.pThreshold)
  _graphite_analyser3->_parameters.Rup = 1.0;
  #define Rup (_graphite_analyser3->_parameters.Rup)
  _graphite_analyser3->_parameters.Rdown = 1.0;
  #define Rdown (_graphite_analyser3->_parameters.Rdown)
  _graphite_analyser3->_parameters.debug = 0;
  #define debug (_graphite_analyser3->_parameters.debug)

  #define mos_rms (_graphite_analyser3->_parameters.mos_rms)
  #define d_rms (_graphite_analyser3->_parameters.d_rms)
  #define mono_Q (_graphite_analyser3->_parameters.mono_Q)
  #define FN (_graphite_analyser3->_parameters.FN)
  #define FM (_graphite_analyser3->_parameters.FM)

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
  /* component graphite_analyser3=Monochromator_pol() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (90)*DEG2RAD, (7.062)*DEG2RAD);
    rot_mul(tr1, _armAnalyzer->_rotation_absolute, _graphite_analyser3->_rotation_absolute);
    rot_transpose(_graphite_analyser2->_rotation_absolute, tr1);
    rot_mul(_graphite_analyser3->_rotation_absolute, tr1, _graphite_analyser3->_rotation_relative);
    _graphite_analyser3->_rotation_is_identity =  rot_test_identity(_graphite_analyser3->_rotation_relative);
    tc1 = coords_set(
      0, -0.165, 0.903887);
    rot_transpose(_armAnalyzer->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _graphite_analyser3->_position_absolute = coords_add(_armAnalyzer->_position_absolute, tc2);
    tc1 = coords_sub(_graphite_analyser2->_position_absolute, _graphite_analyser3->_position_absolute);
    _graphite_analyser3->_position_relative = rot_apply(_graphite_analyser3->_rotation_absolute, tc1);
  } /* graphite_analyser3=Monochromator_pol() AT ROTATED */
  DEBUG_COMPONENT("graphite_analyser3", _graphite_analyser3->_position_absolute, _graphite_analyser3->_rotation_absolute);
  instrument->_position_absolute[31] = _graphite_analyser3->_position_absolute;
  instrument->_position_relative[31] = _graphite_analyser3->_position_relative;
  instrument->counter_N[31]  = instrument->counter_P[31] = instrument->counter_P2[31] = 0;
  instrument->counter_AbsorbProp[31]= 0;
  return(0);
} /* _graphite_analyser3_setpos */

/* component graphite_analyser4=Monochromator_pol() SETTING, POSITION/ROTATION */
int _graphite_analyser4_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_graphite_analyser4_setpos] component graphite_analyser4=Monochromator_pol() SETTING [/usr/share/mcstas/3.0-dev/optics/Monochromator_pol.comp:109]");
  stracpy(_graphite_analyser4->_name, "graphite_analyser4", 16384);
  stracpy(_graphite_analyser4->_type, "Monochromator_pol", 16384);
  _graphite_analyser4->_index=32;
  _graphite_analyser4->_parameters.zwidth = ANAWIDTH;
  #define zwidth (_graphite_analyser4->_parameters.zwidth)
  _graphite_analyser4->_parameters.yheight = ANAHEIGHT;
  #define yheight (_graphite_analyser4->_parameters.yheight)
  _graphite_analyser4->_parameters.mosaic = ANAMOSAIC;
  #define mosaic (_graphite_analyser4->_parameters.mosaic)
  _graphite_analyser4->_parameters.dspread = 0.0025;
  #define dspread (_graphite_analyser4->_parameters.dspread)
  _graphite_analyser4->_parameters.Q = ANAQ;
  #define Q (_graphite_analyser4->_parameters.Q)
  _graphite_analyser4->_parameters.DM = 0;
  #define DM (_graphite_analyser4->_parameters.DM)
  _graphite_analyser4->_parameters.pThreshold = 0.001;
  #define pThreshold (_graphite_analyser4->_parameters.pThreshold)
  _graphite_analyser4->_parameters.Rup = 1.0;
  #define Rup (_graphite_analyser4->_parameters.Rup)
  _graphite_analyser4->_parameters.Rdown = 1.0;
  #define Rdown (_graphite_analyser4->_parameters.Rdown)
  _graphite_analyser4->_parameters.debug = 0;
  #define debug (_graphite_analyser4->_parameters.debug)

  #define mos_rms (_graphite_analyser4->_parameters.mos_rms)
  #define d_rms (_graphite_analyser4->_parameters.d_rms)
  #define mono_Q (_graphite_analyser4->_parameters.mono_Q)
  #define FN (_graphite_analyser4->_parameters.FN)
  #define FM (_graphite_analyser4->_parameters.FM)

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
  /* component graphite_analyser4=Monochromator_pol() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (90)*DEG2RAD, (6.791)*DEG2RAD);
    rot_mul(tr1, _armAnalyzer->_rotation_absolute, _graphite_analyser4->_rotation_absolute);
    rot_transpose(_graphite_analyser3->_rotation_absolute, tr1);
    rot_mul(_graphite_analyser4->_rotation_absolute, tr1, _graphite_analyser4->_rotation_relative);
    _graphite_analyser4->_rotation_is_identity =  rot_test_identity(_graphite_analyser4->_rotation_relative);
    tc1 = coords_set(
      0, -0.155, 0.905081);
    rot_transpose(_armAnalyzer->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _graphite_analyser4->_position_absolute = coords_add(_armAnalyzer->_position_absolute, tc2);
    tc1 = coords_sub(_graphite_analyser3->_position_absolute, _graphite_analyser4->_position_absolute);
    _graphite_analyser4->_position_relative = rot_apply(_graphite_analyser4->_rotation_absolute, tc1);
  } /* graphite_analyser4=Monochromator_pol() AT ROTATED */
  DEBUG_COMPONENT("graphite_analyser4", _graphite_analyser4->_position_absolute, _graphite_analyser4->_rotation_absolute);
  instrument->_position_absolute[32] = _graphite_analyser4->_position_absolute;
  instrument->_position_relative[32] = _graphite_analyser4->_position_relative;
  instrument->counter_N[32]  = instrument->counter_P[32] = instrument->counter_P2[32] = 0;
  instrument->counter_AbsorbProp[32]= 0;
  return(0);
} /* _graphite_analyser4_setpos */

/* component graphite_analyser5=Monochromator_pol() SETTING, POSITION/ROTATION */
int _graphite_analyser5_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_graphite_analyser5_setpos] component graphite_analyser5=Monochromator_pol() SETTING [/usr/share/mcstas/3.0-dev/optics/Monochromator_pol.comp:109]");
  stracpy(_graphite_analyser5->_name, "graphite_analyser5", 16384);
  stracpy(_graphite_analyser5->_type, "Monochromator_pol", 16384);
  _graphite_analyser5->_index=33;
  _graphite_analyser5->_parameters.zwidth = ANAWIDTH;
  #define zwidth (_graphite_analyser5->_parameters.zwidth)
  _graphite_analyser5->_parameters.yheight = ANAHEIGHT;
  #define yheight (_graphite_analyser5->_parameters.yheight)
  _graphite_analyser5->_parameters.mosaic = ANAMOSAIC;
  #define mosaic (_graphite_analyser5->_parameters.mosaic)
  _graphite_analyser5->_parameters.dspread = 0.0025;
  #define dspread (_graphite_analyser5->_parameters.dspread)
  _graphite_analyser5->_parameters.Q = ANAQ;
  #define Q (_graphite_analyser5->_parameters.Q)
  _graphite_analyser5->_parameters.DM = 0;
  #define DM (_graphite_analyser5->_parameters.DM)
  _graphite_analyser5->_parameters.pThreshold = 0.001;
  #define pThreshold (_graphite_analyser5->_parameters.pThreshold)
  _graphite_analyser5->_parameters.Rup = 1.0;
  #define Rup (_graphite_analyser5->_parameters.Rup)
  _graphite_analyser5->_parameters.Rdown = 1.0;
  #define Rdown (_graphite_analyser5->_parameters.Rdown)
  _graphite_analyser5->_parameters.debug = 0;
  #define debug (_graphite_analyser5->_parameters.debug)

  #define mos_rms (_graphite_analyser5->_parameters.mos_rms)
  #define d_rms (_graphite_analyser5->_parameters.d_rms)
  #define mono_Q (_graphite_analyser5->_parameters.mono_Q)
  #define FN (_graphite_analyser5->_parameters.FN)
  #define FM (_graphite_analyser5->_parameters.FM)

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
  /* component graphite_analyser5=Monochromator_pol() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (90)*DEG2RAD, (5.898)*DEG2RAD);
    rot_mul(tr1, _armAnalyzer->_rotation_absolute, _graphite_analyser5->_rotation_absolute);
    rot_transpose(_graphite_analyser4->_rotation_absolute, tr1);
    rot_mul(_graphite_analyser5->_rotation_absolute, tr1, _graphite_analyser5->_rotation_relative);
    _graphite_analyser5->_rotation_is_identity =  rot_test_identity(_graphite_analyser5->_rotation_relative);
    tc1 = coords_set(
      0, -0.145, 0.906123);
    rot_transpose(_armAnalyzer->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _graphite_analyser5->_position_absolute = coords_add(_armAnalyzer->_position_absolute, tc2);
    tc1 = coords_sub(_graphite_analyser4->_position_absolute, _graphite_analyser5->_position_absolute);
    _graphite_analyser5->_position_relative = rot_apply(_graphite_analyser5->_rotation_absolute, tc1);
  } /* graphite_analyser5=Monochromator_pol() AT ROTATED */
  DEBUG_COMPONENT("graphite_analyser5", _graphite_analyser5->_position_absolute, _graphite_analyser5->_rotation_absolute);
  instrument->_position_absolute[33] = _graphite_analyser5->_position_absolute;
  instrument->_position_relative[33] = _graphite_analyser5->_position_relative;
  instrument->counter_N[33]  = instrument->counter_P[33] = instrument->counter_P2[33] = 0;
  instrument->counter_AbsorbProp[33]= 0;
  return(0);
} /* _graphite_analyser5_setpos */

/* component graphite_analyser6=Monochromator_pol() SETTING, POSITION/ROTATION */
int _graphite_analyser6_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_graphite_analyser6_setpos] component graphite_analyser6=Monochromator_pol() SETTING [/usr/share/mcstas/3.0-dev/optics/Monochromator_pol.comp:109]");
  stracpy(_graphite_analyser6->_name, "graphite_analyser6", 16384);
  stracpy(_graphite_analyser6->_type, "Monochromator_pol", 16384);
  _graphite_analyser6->_index=34;
  _graphite_analyser6->_parameters.zwidth = ANAWIDTH;
  #define zwidth (_graphite_analyser6->_parameters.zwidth)
  _graphite_analyser6->_parameters.yheight = ANAHEIGHT;
  #define yheight (_graphite_analyser6->_parameters.yheight)
  _graphite_analyser6->_parameters.mosaic = ANAMOSAIC;
  #define mosaic (_graphite_analyser6->_parameters.mosaic)
  _graphite_analyser6->_parameters.dspread = 0.0025;
  #define dspread (_graphite_analyser6->_parameters.dspread)
  _graphite_analyser6->_parameters.Q = ANAQ;
  #define Q (_graphite_analyser6->_parameters.Q)
  _graphite_analyser6->_parameters.DM = 0;
  #define DM (_graphite_analyser6->_parameters.DM)
  _graphite_analyser6->_parameters.pThreshold = 0.001;
  #define pThreshold (_graphite_analyser6->_parameters.pThreshold)
  _graphite_analyser6->_parameters.Rup = 1.0;
  #define Rup (_graphite_analyser6->_parameters.Rup)
  _graphite_analyser6->_parameters.Rdown = 1.0;
  #define Rdown (_graphite_analyser6->_parameters.Rdown)
  _graphite_analyser6->_parameters.debug = 0;
  #define debug (_graphite_analyser6->_parameters.debug)

  #define mos_rms (_graphite_analyser6->_parameters.mos_rms)
  #define d_rms (_graphite_analyser6->_parameters.d_rms)
  #define mono_Q (_graphite_analyser6->_parameters.mono_Q)
  #define FN (_graphite_analyser6->_parameters.FN)
  #define FM (_graphite_analyser6->_parameters.FM)

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
  /* component graphite_analyser6=Monochromator_pol() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (90)*DEG2RAD, (4.8325)*DEG2RAD);
    rot_mul(tr1, _armAnalyzer->_rotation_absolute, _graphite_analyser6->_rotation_absolute);
    rot_transpose(_graphite_analyser5->_rotation_absolute, tr1);
    rot_mul(_graphite_analyser6->_rotation_absolute, tr1, _graphite_analyser6->_rotation_relative);
    _graphite_analyser6->_rotation_is_identity =  rot_test_identity(_graphite_analyser6->_rotation_relative);
    tc1 = coords_set(
      0, -0.135, 0.90698);
    rot_transpose(_armAnalyzer->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _graphite_analyser6->_position_absolute = coords_add(_armAnalyzer->_position_absolute, tc2);
    tc1 = coords_sub(_graphite_analyser5->_position_absolute, _graphite_analyser6->_position_absolute);
    _graphite_analyser6->_position_relative = rot_apply(_graphite_analyser6->_rotation_absolute, tc1);
  } /* graphite_analyser6=Monochromator_pol() AT ROTATED */
  DEBUG_COMPONENT("graphite_analyser6", _graphite_analyser6->_position_absolute, _graphite_analyser6->_rotation_absolute);
  instrument->_position_absolute[34] = _graphite_analyser6->_position_absolute;
  instrument->_position_relative[34] = _graphite_analyser6->_position_relative;
  instrument->counter_N[34]  = instrument->counter_P[34] = instrument->counter_P2[34] = 0;
  instrument->counter_AbsorbProp[34]= 0;
  return(0);
} /* _graphite_analyser6_setpos */

/* component graphite_analyser7=Monochromator_pol() SETTING, POSITION/ROTATION */
int _graphite_analyser7_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_graphite_analyser7_setpos] component graphite_analyser7=Monochromator_pol() SETTING [/usr/share/mcstas/3.0-dev/optics/Monochromator_pol.comp:109]");
  stracpy(_graphite_analyser7->_name, "graphite_analyser7", 16384);
  stracpy(_graphite_analyser7->_type, "Monochromator_pol", 16384);
  _graphite_analyser7->_index=35;
  _graphite_analyser7->_parameters.zwidth = ANAWIDTH;
  #define zwidth (_graphite_analyser7->_parameters.zwidth)
  _graphite_analyser7->_parameters.yheight = ANAHEIGHT;
  #define yheight (_graphite_analyser7->_parameters.yheight)
  _graphite_analyser7->_parameters.mosaic = ANAMOSAIC;
  #define mosaic (_graphite_analyser7->_parameters.mosaic)
  _graphite_analyser7->_parameters.dspread = 0.0025;
  #define dspread (_graphite_analyser7->_parameters.dspread)
  _graphite_analyser7->_parameters.Q = ANAQ;
  #define Q (_graphite_analyser7->_parameters.Q)
  _graphite_analyser7->_parameters.DM = 0;
  #define DM (_graphite_analyser7->_parameters.DM)
  _graphite_analyser7->_parameters.pThreshold = 0.001;
  #define pThreshold (_graphite_analyser7->_parameters.pThreshold)
  _graphite_analyser7->_parameters.Rup = 1.0;
  #define Rup (_graphite_analyser7->_parameters.Rup)
  _graphite_analyser7->_parameters.Rdown = 1.0;
  #define Rdown (_graphite_analyser7->_parameters.Rdown)
  _graphite_analyser7->_parameters.debug = 0;
  #define debug (_graphite_analyser7->_parameters.debug)

  #define mos_rms (_graphite_analyser7->_parameters.mos_rms)
  #define d_rms (_graphite_analyser7->_parameters.d_rms)
  #define mono_Q (_graphite_analyser7->_parameters.mono_Q)
  #define FN (_graphite_analyser7->_parameters.FN)
  #define FM (_graphite_analyser7->_parameters.FM)

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
  /* component graphite_analyser7=Monochromator_pol() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (90)*DEG2RAD, (4.66)*DEG2RAD);
    rot_mul(tr1, _armAnalyzer->_rotation_absolute, _graphite_analyser7->_rotation_absolute);
    rot_transpose(_graphite_analyser6->_rotation_absolute, tr1);
    rot_mul(_graphite_analyser7->_rotation_absolute, tr1, _graphite_analyser7->_rotation_relative);
    _graphite_analyser7->_rotation_is_identity =  rot_test_identity(_graphite_analyser7->_rotation_relative);
    tc1 = coords_set(
      0, -0.125, 0.907797);
    rot_transpose(_armAnalyzer->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _graphite_analyser7->_position_absolute = coords_add(_armAnalyzer->_position_absolute, tc2);
    tc1 = coords_sub(_graphite_analyser6->_position_absolute, _graphite_analyser7->_position_absolute);
    _graphite_analyser7->_position_relative = rot_apply(_graphite_analyser7->_rotation_absolute, tc1);
  } /* graphite_analyser7=Monochromator_pol() AT ROTATED */
  DEBUG_COMPONENT("graphite_analyser7", _graphite_analyser7->_position_absolute, _graphite_analyser7->_rotation_absolute);
  instrument->_position_absolute[35] = _graphite_analyser7->_position_absolute;
  instrument->_position_relative[35] = _graphite_analyser7->_position_relative;
  instrument->counter_N[35]  = instrument->counter_P[35] = instrument->counter_P2[35] = 0;
  instrument->counter_AbsorbProp[35]= 0;
  return(0);
} /* _graphite_analyser7_setpos */

/* component graphite_analyser8=Monochromator_pol() SETTING, POSITION/ROTATION */
int _graphite_analyser8_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_graphite_analyser8_setpos] component graphite_analyser8=Monochromator_pol() SETTING [/usr/share/mcstas/3.0-dev/optics/Monochromator_pol.comp:109]");
  stracpy(_graphite_analyser8->_name, "graphite_analyser8", 16384);
  stracpy(_graphite_analyser8->_type, "Monochromator_pol", 16384);
  _graphite_analyser8->_index=36;
  _graphite_analyser8->_parameters.zwidth = ANAWIDTH;
  #define zwidth (_graphite_analyser8->_parameters.zwidth)
  _graphite_analyser8->_parameters.yheight = ANAHEIGHT;
  #define yheight (_graphite_analyser8->_parameters.yheight)
  _graphite_analyser8->_parameters.mosaic = ANAMOSAIC;
  #define mosaic (_graphite_analyser8->_parameters.mosaic)
  _graphite_analyser8->_parameters.dspread = 0.0025;
  #define dspread (_graphite_analyser8->_parameters.dspread)
  _graphite_analyser8->_parameters.Q = ANAQ;
  #define Q (_graphite_analyser8->_parameters.Q)
  _graphite_analyser8->_parameters.DM = 0;
  #define DM (_graphite_analyser8->_parameters.DM)
  _graphite_analyser8->_parameters.pThreshold = 0.001;
  #define pThreshold (_graphite_analyser8->_parameters.pThreshold)
  _graphite_analyser8->_parameters.Rup = 1.0;
  #define Rup (_graphite_analyser8->_parameters.Rup)
  _graphite_analyser8->_parameters.Rdown = 1.0;
  #define Rdown (_graphite_analyser8->_parameters.Rdown)
  _graphite_analyser8->_parameters.debug = 0;
  #define debug (_graphite_analyser8->_parameters.debug)

  #define mos_rms (_graphite_analyser8->_parameters.mos_rms)
  #define d_rms (_graphite_analyser8->_parameters.d_rms)
  #define mono_Q (_graphite_analyser8->_parameters.mono_Q)
  #define FN (_graphite_analyser8->_parameters.FN)
  #define FM (_graphite_analyser8->_parameters.FM)

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
  /* component graphite_analyser8=Monochromator_pol() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (90)*DEG2RAD, (3.6815)*DEG2RAD);
    rot_mul(tr1, _armAnalyzer->_rotation_absolute, _graphite_analyser8->_rotation_absolute);
    rot_transpose(_graphite_analyser7->_rotation_absolute, tr1);
    rot_mul(_graphite_analyser8->_rotation_absolute, tr1, _graphite_analyser8->_rotation_relative);
    _graphite_analyser8->_rotation_is_identity =  rot_test_identity(_graphite_analyser8->_rotation_relative);
    tc1 = coords_set(
      0, -0.115, 0.90845);
    rot_transpose(_armAnalyzer->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _graphite_analyser8->_position_absolute = coords_add(_armAnalyzer->_position_absolute, tc2);
    tc1 = coords_sub(_graphite_analyser7->_position_absolute, _graphite_analyser8->_position_absolute);
    _graphite_analyser8->_position_relative = rot_apply(_graphite_analyser8->_rotation_absolute, tc1);
  } /* graphite_analyser8=Monochromator_pol() AT ROTATED */
  DEBUG_COMPONENT("graphite_analyser8", _graphite_analyser8->_position_absolute, _graphite_analyser8->_rotation_absolute);
  instrument->_position_absolute[36] = _graphite_analyser8->_position_absolute;
  instrument->_position_relative[36] = _graphite_analyser8->_position_relative;
  instrument->counter_N[36]  = instrument->counter_P[36] = instrument->counter_P2[36] = 0;
  instrument->counter_AbsorbProp[36]= 0;
  return(0);
} /* _graphite_analyser8_setpos */

/* component graphite_analyser9=Monochromator_pol() SETTING, POSITION/ROTATION */
int _graphite_analyser9_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_graphite_analyser9_setpos] component graphite_analyser9=Monochromator_pol() SETTING [/usr/share/mcstas/3.0-dev/optics/Monochromator_pol.comp:109]");
  stracpy(_graphite_analyser9->_name, "graphite_analyser9", 16384);
  stracpy(_graphite_analyser9->_type, "Monochromator_pol", 16384);
  _graphite_analyser9->_index=37;
  _graphite_analyser9->_parameters.zwidth = ANAWIDTH;
  #define zwidth (_graphite_analyser9->_parameters.zwidth)
  _graphite_analyser9->_parameters.yheight = ANAHEIGHT;
  #define yheight (_graphite_analyser9->_parameters.yheight)
  _graphite_analyser9->_parameters.mosaic = ANAMOSAIC;
  #define mosaic (_graphite_analyser9->_parameters.mosaic)
  _graphite_analyser9->_parameters.dspread = 0.0025;
  #define dspread (_graphite_analyser9->_parameters.dspread)
  _graphite_analyser9->_parameters.Q = ANAQ;
  #define Q (_graphite_analyser9->_parameters.Q)
  _graphite_analyser9->_parameters.DM = 0;
  #define DM (_graphite_analyser9->_parameters.DM)
  _graphite_analyser9->_parameters.pThreshold = 0.001;
  #define pThreshold (_graphite_analyser9->_parameters.pThreshold)
  _graphite_analyser9->_parameters.Rup = 1.0;
  #define Rup (_graphite_analyser9->_parameters.Rup)
  _graphite_analyser9->_parameters.Rdown = 1.0;
  #define Rdown (_graphite_analyser9->_parameters.Rdown)
  _graphite_analyser9->_parameters.debug = 0;
  #define debug (_graphite_analyser9->_parameters.debug)

  #define mos_rms (_graphite_analyser9->_parameters.mos_rms)
  #define d_rms (_graphite_analyser9->_parameters.d_rms)
  #define mono_Q (_graphite_analyser9->_parameters.mono_Q)
  #define FN (_graphite_analyser9->_parameters.FN)
  #define FM (_graphite_analyser9->_parameters.FM)

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
  /* component graphite_analyser9=Monochromator_pol() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (90)*DEG2RAD, (2.6165)*DEG2RAD);
    rot_mul(tr1, _armAnalyzer->_rotation_absolute, _graphite_analyser9->_rotation_absolute);
    rot_transpose(_graphite_analyser8->_rotation_absolute, tr1);
    rot_mul(_graphite_analyser9->_rotation_absolute, tr1, _graphite_analyser9->_rotation_relative);
    _graphite_analyser9->_rotation_is_identity =  rot_test_identity(_graphite_analyser9->_rotation_relative);
    tc1 = coords_set(
      0, -0.105, 0.908918);
    rot_transpose(_armAnalyzer->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _graphite_analyser9->_position_absolute = coords_add(_armAnalyzer->_position_absolute, tc2);
    tc1 = coords_sub(_graphite_analyser8->_position_absolute, _graphite_analyser9->_position_absolute);
    _graphite_analyser9->_position_relative = rot_apply(_graphite_analyser9->_rotation_absolute, tc1);
  } /* graphite_analyser9=Monochromator_pol() AT ROTATED */
  DEBUG_COMPONENT("graphite_analyser9", _graphite_analyser9->_position_absolute, _graphite_analyser9->_rotation_absolute);
  instrument->_position_absolute[37] = _graphite_analyser9->_position_absolute;
  instrument->_position_relative[37] = _graphite_analyser9->_position_relative;
  instrument->counter_N[37]  = instrument->counter_P[37] = instrument->counter_P2[37] = 0;
  instrument->counter_AbsorbProp[37]= 0;
  return(0);
} /* _graphite_analyser9_setpos */

/* component graphite_analyser10=Monochromator_pol() SETTING, POSITION/ROTATION */
int _graphite_analyser10_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_graphite_analyser10_setpos] component graphite_analyser10=Monochromator_pol() SETTING [/usr/share/mcstas/3.0-dev/optics/Monochromator_pol.comp:109]");
  stracpy(_graphite_analyser10->_name, "graphite_analyser10", 16384);
  stracpy(_graphite_analyser10->_type, "Monochromator_pol", 16384);
  _graphite_analyser10->_index=38;
  _graphite_analyser10->_parameters.zwidth = ANAWIDTH;
  #define zwidth (_graphite_analyser10->_parameters.zwidth)
  _graphite_analyser10->_parameters.yheight = ANAHEIGHT;
  #define yheight (_graphite_analyser10->_parameters.yheight)
  _graphite_analyser10->_parameters.mosaic = ANAMOSAIC;
  #define mosaic (_graphite_analyser10->_parameters.mosaic)
  _graphite_analyser10->_parameters.dspread = 0.0025;
  #define dspread (_graphite_analyser10->_parameters.dspread)
  _graphite_analyser10->_parameters.Q = ANAQ;
  #define Q (_graphite_analyser10->_parameters.Q)
  _graphite_analyser10->_parameters.DM = 0;
  #define DM (_graphite_analyser10->_parameters.DM)
  _graphite_analyser10->_parameters.pThreshold = 0.001;
  #define pThreshold (_graphite_analyser10->_parameters.pThreshold)
  _graphite_analyser10->_parameters.Rup = 1.0;
  #define Rup (_graphite_analyser10->_parameters.Rup)
  _graphite_analyser10->_parameters.Rdown = 1.0;
  #define Rdown (_graphite_analyser10->_parameters.Rdown)
  _graphite_analyser10->_parameters.debug = 0;
  #define debug (_graphite_analyser10->_parameters.debug)

  #define mos_rms (_graphite_analyser10->_parameters.mos_rms)
  #define d_rms (_graphite_analyser10->_parameters.d_rms)
  #define mono_Q (_graphite_analyser10->_parameters.mono_Q)
  #define FN (_graphite_analyser10->_parameters.FN)
  #define FM (_graphite_analyser10->_parameters.FM)

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
  /* component graphite_analyser10=Monochromator_pol() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (90)*DEG2RAD, (2.53)*DEG2RAD);
    rot_mul(tr1, _armAnalyzer->_rotation_absolute, _graphite_analyser10->_rotation_absolute);
    rot_transpose(_graphite_analyser9->_rotation_absolute, tr1);
    rot_mul(_graphite_analyser10->_rotation_absolute, tr1, _graphite_analyser10->_rotation_relative);
    _graphite_analyser10->_rotation_is_identity =  rot_test_identity(_graphite_analyser10->_rotation_relative);
    tc1 = coords_set(
      0, -0.095, 0.90936);
    rot_transpose(_armAnalyzer->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _graphite_analyser10->_position_absolute = coords_add(_armAnalyzer->_position_absolute, tc2);
    tc1 = coords_sub(_graphite_analyser9->_position_absolute, _graphite_analyser10->_position_absolute);
    _graphite_analyser10->_position_relative = rot_apply(_graphite_analyser10->_rotation_absolute, tc1);
  } /* graphite_analyser10=Monochromator_pol() AT ROTATED */
  DEBUG_COMPONENT("graphite_analyser10", _graphite_analyser10->_position_absolute, _graphite_analyser10->_rotation_absolute);
  instrument->_position_absolute[38] = _graphite_analyser10->_position_absolute;
  instrument->_position_relative[38] = _graphite_analyser10->_position_relative;
  instrument->counter_N[38]  = instrument->counter_P[38] = instrument->counter_P2[38] = 0;
  instrument->counter_AbsorbProp[38]= 0;
  return(0);
} /* _graphite_analyser10_setpos */

/* component graphite_analyser11=Monochromator_pol() SETTING, POSITION/ROTATION */
int _graphite_analyser11_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_graphite_analyser11_setpos] component graphite_analyser11=Monochromator_pol() SETTING [/usr/share/mcstas/3.0-dev/optics/Monochromator_pol.comp:109]");
  stracpy(_graphite_analyser11->_name, "graphite_analyser11", 16384);
  stracpy(_graphite_analyser11->_type, "Monochromator_pol", 16384);
  _graphite_analyser11->_index=39;
  _graphite_analyser11->_parameters.zwidth = ANAWIDTH;
  #define zwidth (_graphite_analyser11->_parameters.zwidth)
  _graphite_analyser11->_parameters.yheight = ANAHEIGHT;
  #define yheight (_graphite_analyser11->_parameters.yheight)
  _graphite_analyser11->_parameters.mosaic = ANAMOSAIC;
  #define mosaic (_graphite_analyser11->_parameters.mosaic)
  _graphite_analyser11->_parameters.dspread = 0.0025;
  #define dspread (_graphite_analyser11->_parameters.dspread)
  _graphite_analyser11->_parameters.Q = ANAQ;
  #define Q (_graphite_analyser11->_parameters.Q)
  _graphite_analyser11->_parameters.DM = 0;
  #define DM (_graphite_analyser11->_parameters.DM)
  _graphite_analyser11->_parameters.pThreshold = 0.001;
  #define pThreshold (_graphite_analyser11->_parameters.pThreshold)
  _graphite_analyser11->_parameters.Rup = 1.0;
  #define Rup (_graphite_analyser11->_parameters.Rup)
  _graphite_analyser11->_parameters.Rdown = 1.0;
  #define Rdown (_graphite_analyser11->_parameters.Rdown)
  _graphite_analyser11->_parameters.debug = 0;
  #define debug (_graphite_analyser11->_parameters.debug)

  #define mos_rms (_graphite_analyser11->_parameters.mos_rms)
  #define d_rms (_graphite_analyser11->_parameters.d_rms)
  #define mono_Q (_graphite_analyser11->_parameters.mono_Q)
  #define FN (_graphite_analyser11->_parameters.FN)
  #define FM (_graphite_analyser11->_parameters.FM)

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
  /* component graphite_analyser11=Monochromator_pol() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (90)*DEG2RAD, (1.465)*DEG2RAD);
    rot_mul(tr1, _armAnalyzer->_rotation_absolute, _graphite_analyser11->_rotation_absolute);
    rot_transpose(_graphite_analyser10->_rotation_absolute, tr1);
    rot_mul(_graphite_analyser11->_rotation_absolute, tr1, _graphite_analyser11->_rotation_relative);
    _graphite_analyser11->_rotation_is_identity =  rot_test_identity(_graphite_analyser11->_rotation_relative);
    tc1 = coords_set(
      0, -0.085, 0.909628);
    rot_transpose(_armAnalyzer->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _graphite_analyser11->_position_absolute = coords_add(_armAnalyzer->_position_absolute, tc2);
    tc1 = coords_sub(_graphite_analyser10->_position_absolute, _graphite_analyser11->_position_absolute);
    _graphite_analyser11->_position_relative = rot_apply(_graphite_analyser11->_rotation_absolute, tc1);
  } /* graphite_analyser11=Monochromator_pol() AT ROTATED */
  DEBUG_COMPONENT("graphite_analyser11", _graphite_analyser11->_position_absolute, _graphite_analyser11->_rotation_absolute);
  instrument->_position_absolute[39] = _graphite_analyser11->_position_absolute;
  instrument->_position_relative[39] = _graphite_analyser11->_position_relative;
  instrument->counter_N[39]  = instrument->counter_P[39] = instrument->counter_P2[39] = 0;
  instrument->counter_AbsorbProp[39]= 0;
  return(0);
} /* _graphite_analyser11_setpos */

/* component graphite_analyser12=Monochromator_pol() SETTING, POSITION/ROTATION */
int _graphite_analyser12_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_graphite_analyser12_setpos] component graphite_analyser12=Monochromator_pol() SETTING [/usr/share/mcstas/3.0-dev/optics/Monochromator_pol.comp:109]");
  stracpy(_graphite_analyser12->_name, "graphite_analyser12", 16384);
  stracpy(_graphite_analyser12->_type, "Monochromator_pol", 16384);
  _graphite_analyser12->_index=40;
  _graphite_analyser12->_parameters.zwidth = ANAWIDTH;
  #define zwidth (_graphite_analyser12->_parameters.zwidth)
  _graphite_analyser12->_parameters.yheight = ANAHEIGHT;
  #define yheight (_graphite_analyser12->_parameters.yheight)
  _graphite_analyser12->_parameters.mosaic = ANAMOSAIC;
  #define mosaic (_graphite_analyser12->_parameters.mosaic)
  _graphite_analyser12->_parameters.dspread = 0.0025;
  #define dspread (_graphite_analyser12->_parameters.dspread)
  _graphite_analyser12->_parameters.Q = ANAQ;
  #define Q (_graphite_analyser12->_parameters.Q)
  _graphite_analyser12->_parameters.DM = 0;
  #define DM (_graphite_analyser12->_parameters.DM)
  _graphite_analyser12->_parameters.pThreshold = 0.001;
  #define pThreshold (_graphite_analyser12->_parameters.pThreshold)
  _graphite_analyser12->_parameters.Rup = 1.0;
  #define Rup (_graphite_analyser12->_parameters.Rup)
  _graphite_analyser12->_parameters.Rdown = 1.0;
  #define Rdown (_graphite_analyser12->_parameters.Rdown)
  _graphite_analyser12->_parameters.debug = 0;
  #define debug (_graphite_analyser12->_parameters.debug)

  #define mos_rms (_graphite_analyser12->_parameters.mos_rms)
  #define d_rms (_graphite_analyser12->_parameters.d_rms)
  #define mono_Q (_graphite_analyser12->_parameters.mono_Q)
  #define FN (_graphite_analyser12->_parameters.FN)
  #define FM (_graphite_analyser12->_parameters.FM)

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
  /* component graphite_analyser12=Monochromator_pol() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (90)*DEG2RAD, (0.4)*DEG2RAD);
    rot_mul(tr1, _armAnalyzer->_rotation_absolute, _graphite_analyser12->_rotation_absolute);
    rot_transpose(_graphite_analyser11->_rotation_absolute, tr1);
    rot_mul(_graphite_analyser12->_rotation_absolute, tr1, _graphite_analyser12->_rotation_relative);
    _graphite_analyser12->_rotation_is_identity =  rot_test_identity(_graphite_analyser12->_rotation_relative);
    tc1 = coords_set(
      0, -0.075, 0.90971);
    rot_transpose(_armAnalyzer->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _graphite_analyser12->_position_absolute = coords_add(_armAnalyzer->_position_absolute, tc2);
    tc1 = coords_sub(_graphite_analyser11->_position_absolute, _graphite_analyser12->_position_absolute);
    _graphite_analyser12->_position_relative = rot_apply(_graphite_analyser12->_rotation_absolute, tc1);
  } /* graphite_analyser12=Monochromator_pol() AT ROTATED */
  DEBUG_COMPONENT("graphite_analyser12", _graphite_analyser12->_position_absolute, _graphite_analyser12->_rotation_absolute);
  instrument->_position_absolute[40] = _graphite_analyser12->_position_absolute;
  instrument->_position_relative[40] = _graphite_analyser12->_position_relative;
  instrument->counter_N[40]  = instrument->counter_P[40] = instrument->counter_P2[40] = 0;
  instrument->counter_AbsorbProp[40]= 0;
  return(0);
} /* _graphite_analyser12_setpos */

/* component graphite_analyser13=Monochromator_pol() SETTING, POSITION/ROTATION */
int _graphite_analyser13_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_graphite_analyser13_setpos] component graphite_analyser13=Monochromator_pol() SETTING [/usr/share/mcstas/3.0-dev/optics/Monochromator_pol.comp:109]");
  stracpy(_graphite_analyser13->_name, "graphite_analyser13", 16384);
  stracpy(_graphite_analyser13->_type, "Monochromator_pol", 16384);
  _graphite_analyser13->_index=41;
  _graphite_analyser13->_parameters.zwidth = ANAWIDTH;
  #define zwidth (_graphite_analyser13->_parameters.zwidth)
  _graphite_analyser13->_parameters.yheight = ANAHEIGHT;
  #define yheight (_graphite_analyser13->_parameters.yheight)
  _graphite_analyser13->_parameters.mosaic = ANAMOSAIC;
  #define mosaic (_graphite_analyser13->_parameters.mosaic)
  _graphite_analyser13->_parameters.dspread = 0.0025;
  #define dspread (_graphite_analyser13->_parameters.dspread)
  _graphite_analyser13->_parameters.Q = ANAQ;
  #define Q (_graphite_analyser13->_parameters.Q)
  _graphite_analyser13->_parameters.DM = 0;
  #define DM (_graphite_analyser13->_parameters.DM)
  _graphite_analyser13->_parameters.pThreshold = 0.001;
  #define pThreshold (_graphite_analyser13->_parameters.pThreshold)
  _graphite_analyser13->_parameters.Rup = 1.0;
  #define Rup (_graphite_analyser13->_parameters.Rup)
  _graphite_analyser13->_parameters.Rdown = 1.0;
  #define Rdown (_graphite_analyser13->_parameters.Rdown)
  _graphite_analyser13->_parameters.debug = 0;
  #define debug (_graphite_analyser13->_parameters.debug)

  #define mos_rms (_graphite_analyser13->_parameters.mos_rms)
  #define d_rms (_graphite_analyser13->_parameters.d_rms)
  #define mono_Q (_graphite_analyser13->_parameters.mono_Q)
  #define FN (_graphite_analyser13->_parameters.FN)
  #define FM (_graphite_analyser13->_parameters.FM)

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
  /* component graphite_analyser13=Monochromator_pol() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (90)*DEG2RAD, (0.339)*DEG2RAD);
    rot_mul(tr1, _armAnalyzer->_rotation_absolute, _graphite_analyser13->_rotation_absolute);
    rot_transpose(_graphite_analyser12->_rotation_absolute, tr1);
    rot_mul(_graphite_analyser13->_rotation_absolute, tr1, _graphite_analyser13->_rotation_relative);
    _graphite_analyser13->_rotation_is_identity =  rot_test_identity(_graphite_analyser13->_rotation_relative);
    tc1 = coords_set(
      0, -0.065, 0.90977);
    rot_transpose(_armAnalyzer->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _graphite_analyser13->_position_absolute = coords_add(_armAnalyzer->_position_absolute, tc2);
    tc1 = coords_sub(_graphite_analyser12->_position_absolute, _graphite_analyser13->_position_absolute);
    _graphite_analyser13->_position_relative = rot_apply(_graphite_analyser13->_rotation_absolute, tc1);
  } /* graphite_analyser13=Monochromator_pol() AT ROTATED */
  DEBUG_COMPONENT("graphite_analyser13", _graphite_analyser13->_position_absolute, _graphite_analyser13->_rotation_absolute);
  instrument->_position_absolute[41] = _graphite_analyser13->_position_absolute;
  instrument->_position_relative[41] = _graphite_analyser13->_position_relative;
  instrument->counter_N[41]  = instrument->counter_P[41] = instrument->counter_P2[41] = 0;
  instrument->counter_AbsorbProp[41]= 0;
  return(0);
} /* _graphite_analyser13_setpos */

/* component graphite_analyser14=Monochromator_pol() SETTING, POSITION/ROTATION */
int _graphite_analyser14_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_graphite_analyser14_setpos] component graphite_analyser14=Monochromator_pol() SETTING [/usr/share/mcstas/3.0-dev/optics/Monochromator_pol.comp:109]");
  stracpy(_graphite_analyser14->_name, "graphite_analyser14", 16384);
  stracpy(_graphite_analyser14->_type, "Monochromator_pol", 16384);
  _graphite_analyser14->_index=42;
  _graphite_analyser14->_parameters.zwidth = ANAWIDTH;
  #define zwidth (_graphite_analyser14->_parameters.zwidth)
  _graphite_analyser14->_parameters.yheight = ANAHEIGHT;
  #define yheight (_graphite_analyser14->_parameters.yheight)
  _graphite_analyser14->_parameters.mosaic = ANAMOSAIC;
  #define mosaic (_graphite_analyser14->_parameters.mosaic)
  _graphite_analyser14->_parameters.dspread = 0.0025;
  #define dspread (_graphite_analyser14->_parameters.dspread)
  _graphite_analyser14->_parameters.Q = ANAQ;
  #define Q (_graphite_analyser14->_parameters.Q)
  _graphite_analyser14->_parameters.DM = 0;
  #define DM (_graphite_analyser14->_parameters.DM)
  _graphite_analyser14->_parameters.pThreshold = 0.001;
  #define pThreshold (_graphite_analyser14->_parameters.pThreshold)
  _graphite_analyser14->_parameters.Rup = 1.0;
  #define Rup (_graphite_analyser14->_parameters.Rup)
  _graphite_analyser14->_parameters.Rdown = 1.0;
  #define Rdown (_graphite_analyser14->_parameters.Rdown)
  _graphite_analyser14->_parameters.debug = 0;
  #define debug (_graphite_analyser14->_parameters.debug)

  #define mos_rms (_graphite_analyser14->_parameters.mos_rms)
  #define d_rms (_graphite_analyser14->_parameters.d_rms)
  #define mono_Q (_graphite_analyser14->_parameters.mono_Q)
  #define FN (_graphite_analyser14->_parameters.FN)
  #define FM (_graphite_analyser14->_parameters.FM)

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
  /* component graphite_analyser14=Monochromator_pol() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (90)*DEG2RAD, (-0.724)*DEG2RAD);
    rot_mul(tr1, _armAnalyzer->_rotation_absolute, _graphite_analyser14->_rotation_absolute);
    rot_transpose(_graphite_analyser13->_rotation_absolute, tr1);
    rot_mul(_graphite_analyser14->_rotation_absolute, tr1, _graphite_analyser14->_rotation_relative);
    _graphite_analyser14->_rotation_is_identity =  rot_test_identity(_graphite_analyser14->_rotation_relative);
    tc1 = coords_set(
      0, -0.055, 0.909655);
    rot_transpose(_armAnalyzer->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _graphite_analyser14->_position_absolute = coords_add(_armAnalyzer->_position_absolute, tc2);
    tc1 = coords_sub(_graphite_analyser13->_position_absolute, _graphite_analyser14->_position_absolute);
    _graphite_analyser14->_position_relative = rot_apply(_graphite_analyser14->_rotation_absolute, tc1);
  } /* graphite_analyser14=Monochromator_pol() AT ROTATED */
  DEBUG_COMPONENT("graphite_analyser14", _graphite_analyser14->_position_absolute, _graphite_analyser14->_rotation_absolute);
  instrument->_position_absolute[42] = _graphite_analyser14->_position_absolute;
  instrument->_position_relative[42] = _graphite_analyser14->_position_relative;
  instrument->counter_N[42]  = instrument->counter_P[42] = instrument->counter_P2[42] = 0;
  instrument->counter_AbsorbProp[42]= 0;
  return(0);
} /* _graphite_analyser14_setpos */

/* component graphite_analyser15=Monochromator_pol() SETTING, POSITION/ROTATION */
int _graphite_analyser15_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_graphite_analyser15_setpos] component graphite_analyser15=Monochromator_pol() SETTING [/usr/share/mcstas/3.0-dev/optics/Monochromator_pol.comp:109]");
  stracpy(_graphite_analyser15->_name, "graphite_analyser15", 16384);
  stracpy(_graphite_analyser15->_type, "Monochromator_pol", 16384);
  _graphite_analyser15->_index=43;
  _graphite_analyser15->_parameters.zwidth = ANAWIDTH;
  #define zwidth (_graphite_analyser15->_parameters.zwidth)
  _graphite_analyser15->_parameters.yheight = ANAHEIGHT;
  #define yheight (_graphite_analyser15->_parameters.yheight)
  _graphite_analyser15->_parameters.mosaic = ANAMOSAIC;
  #define mosaic (_graphite_analyser15->_parameters.mosaic)
  _graphite_analyser15->_parameters.dspread = 0.0025;
  #define dspread (_graphite_analyser15->_parameters.dspread)
  _graphite_analyser15->_parameters.Q = ANAQ;
  #define Q (_graphite_analyser15->_parameters.Q)
  _graphite_analyser15->_parameters.DM = 0;
  #define DM (_graphite_analyser15->_parameters.DM)
  _graphite_analyser15->_parameters.pThreshold = 0.001;
  #define pThreshold (_graphite_analyser15->_parameters.pThreshold)
  _graphite_analyser15->_parameters.Rup = 1.0;
  #define Rup (_graphite_analyser15->_parameters.Rup)
  _graphite_analyser15->_parameters.Rdown = 1.0;
  #define Rdown (_graphite_analyser15->_parameters.Rdown)
  _graphite_analyser15->_parameters.debug = 0;
  #define debug (_graphite_analyser15->_parameters.debug)

  #define mos_rms (_graphite_analyser15->_parameters.mos_rms)
  #define d_rms (_graphite_analyser15->_parameters.d_rms)
  #define mono_Q (_graphite_analyser15->_parameters.mono_Q)
  #define FN (_graphite_analyser15->_parameters.FN)
  #define FM (_graphite_analyser15->_parameters.FM)

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
  /* component graphite_analyser15=Monochromator_pol() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (90)*DEG2RAD, (-1.726)*DEG2RAD);
    rot_mul(tr1, _armAnalyzer->_rotation_absolute, _graphite_analyser15->_rotation_absolute);
    rot_transpose(_graphite_analyser14->_rotation_absolute, tr1);
    rot_mul(_graphite_analyser15->_rotation_absolute, tr1, _graphite_analyser15->_rotation_relative);
    _graphite_analyser15->_rotation_is_identity =  rot_test_identity(_graphite_analyser15->_rotation_relative);
    tc1 = coords_set(
      0, -0.045, 0.909363);
    rot_transpose(_armAnalyzer->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _graphite_analyser15->_position_absolute = coords_add(_armAnalyzer->_position_absolute, tc2);
    tc1 = coords_sub(_graphite_analyser14->_position_absolute, _graphite_analyser15->_position_absolute);
    _graphite_analyser15->_position_relative = rot_apply(_graphite_analyser15->_rotation_absolute, tc1);
  } /* graphite_analyser15=Monochromator_pol() AT ROTATED */
  DEBUG_COMPONENT("graphite_analyser15", _graphite_analyser15->_position_absolute, _graphite_analyser15->_rotation_absolute);
  instrument->_position_absolute[43] = _graphite_analyser15->_position_absolute;
  instrument->_position_relative[43] = _graphite_analyser15->_position_relative;
  instrument->counter_N[43]  = instrument->counter_P[43] = instrument->counter_P2[43] = 0;
  instrument->counter_AbsorbProp[43]= 0;
  return(0);
} /* _graphite_analyser15_setpos */

/* component graphite_analyser16=Monochromator_pol() SETTING, POSITION/ROTATION */
int _graphite_analyser16_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_graphite_analyser16_setpos] component graphite_analyser16=Monochromator_pol() SETTING [/usr/share/mcstas/3.0-dev/optics/Monochromator_pol.comp:109]");
  stracpy(_graphite_analyser16->_name, "graphite_analyser16", 16384);
  stracpy(_graphite_analyser16->_type, "Monochromator_pol", 16384);
  _graphite_analyser16->_index=44;
  _graphite_analyser16->_parameters.zwidth = ANAWIDTH;
  #define zwidth (_graphite_analyser16->_parameters.zwidth)
  _graphite_analyser16->_parameters.yheight = ANAHEIGHT;
  #define yheight (_graphite_analyser16->_parameters.yheight)
  _graphite_analyser16->_parameters.mosaic = ANAMOSAIC;
  #define mosaic (_graphite_analyser16->_parameters.mosaic)
  _graphite_analyser16->_parameters.dspread = 0.0025;
  #define dspread (_graphite_analyser16->_parameters.dspread)
  _graphite_analyser16->_parameters.Q = ANAQ;
  #define Q (_graphite_analyser16->_parameters.Q)
  _graphite_analyser16->_parameters.DM = 0;
  #define DM (_graphite_analyser16->_parameters.DM)
  _graphite_analyser16->_parameters.pThreshold = 0.001;
  #define pThreshold (_graphite_analyser16->_parameters.pThreshold)
  _graphite_analyser16->_parameters.Rup = 1.0;
  #define Rup (_graphite_analyser16->_parameters.Rup)
  _graphite_analyser16->_parameters.Rdown = 1.0;
  #define Rdown (_graphite_analyser16->_parameters.Rdown)
  _graphite_analyser16->_parameters.debug = 0;
  #define debug (_graphite_analyser16->_parameters.debug)

  #define mos_rms (_graphite_analyser16->_parameters.mos_rms)
  #define d_rms (_graphite_analyser16->_parameters.d_rms)
  #define mono_Q (_graphite_analyser16->_parameters.mono_Q)
  #define FN (_graphite_analyser16->_parameters.FN)
  #define FM (_graphite_analyser16->_parameters.FM)

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
  /* component graphite_analyser16=Monochromator_pol() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (90)*DEG2RAD, (-1.854)*DEG2RAD);
    rot_mul(tr1, _armAnalyzer->_rotation_absolute, _graphite_analyser16->_rotation_absolute);
    rot_transpose(_graphite_analyser15->_rotation_absolute, tr1);
    rot_mul(_graphite_analyser16->_rotation_absolute, tr1, _graphite_analyser16->_rotation_relative);
    _graphite_analyser16->_rotation_is_identity =  rot_test_identity(_graphite_analyser16->_rotation_relative);
    tc1 = coords_set(
      0, -0.035, 0.909041);
    rot_transpose(_armAnalyzer->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _graphite_analyser16->_position_absolute = coords_add(_armAnalyzer->_position_absolute, tc2);
    tc1 = coords_sub(_graphite_analyser15->_position_absolute, _graphite_analyser16->_position_absolute);
    _graphite_analyser16->_position_relative = rot_apply(_graphite_analyser16->_rotation_absolute, tc1);
  } /* graphite_analyser16=Monochromator_pol() AT ROTATED */
  DEBUG_COMPONENT("graphite_analyser16", _graphite_analyser16->_position_absolute, _graphite_analyser16->_rotation_absolute);
  instrument->_position_absolute[44] = _graphite_analyser16->_position_absolute;
  instrument->_position_relative[44] = _graphite_analyser16->_position_relative;
  instrument->counter_N[44]  = instrument->counter_P[44] = instrument->counter_P2[44] = 0;
  instrument->counter_AbsorbProp[44]= 0;
  return(0);
} /* _graphite_analyser16_setpos */

/* component graphite_analyser17=Monochromator_pol() SETTING, POSITION/ROTATION */
int _graphite_analyser17_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_graphite_analyser17_setpos] component graphite_analyser17=Monochromator_pol() SETTING [/usr/share/mcstas/3.0-dev/optics/Monochromator_pol.comp:109]");
  stracpy(_graphite_analyser17->_name, "graphite_analyser17", 16384);
  stracpy(_graphite_analyser17->_type, "Monochromator_pol", 16384);
  _graphite_analyser17->_index=45;
  _graphite_analyser17->_parameters.zwidth = ANAWIDTH;
  #define zwidth (_graphite_analyser17->_parameters.zwidth)
  _graphite_analyser17->_parameters.yheight = ANAHEIGHT;
  #define yheight (_graphite_analyser17->_parameters.yheight)
  _graphite_analyser17->_parameters.mosaic = ANAMOSAIC;
  #define mosaic (_graphite_analyser17->_parameters.mosaic)
  _graphite_analyser17->_parameters.dspread = 0.0025;
  #define dspread (_graphite_analyser17->_parameters.dspread)
  _graphite_analyser17->_parameters.Q = ANAQ;
  #define Q (_graphite_analyser17->_parameters.Q)
  _graphite_analyser17->_parameters.DM = 0;
  #define DM (_graphite_analyser17->_parameters.DM)
  _graphite_analyser17->_parameters.pThreshold = 0.001;
  #define pThreshold (_graphite_analyser17->_parameters.pThreshold)
  _graphite_analyser17->_parameters.Rup = 1.0;
  #define Rup (_graphite_analyser17->_parameters.Rup)
  _graphite_analyser17->_parameters.Rdown = 1.0;
  #define Rdown (_graphite_analyser17->_parameters.Rdown)
  _graphite_analyser17->_parameters.debug = 0;
  #define debug (_graphite_analyser17->_parameters.debug)

  #define mos_rms (_graphite_analyser17->_parameters.mos_rms)
  #define d_rms (_graphite_analyser17->_parameters.d_rms)
  #define mono_Q (_graphite_analyser17->_parameters.mono_Q)
  #define FN (_graphite_analyser17->_parameters.FN)
  #define FM (_graphite_analyser17->_parameters.FM)

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
  /* component graphite_analyser17=Monochromator_pol() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (90)*DEG2RAD, (-2.9185)*DEG2RAD);
    rot_mul(tr1, _armAnalyzer->_rotation_absolute, _graphite_analyser17->_rotation_absolute);
    rot_transpose(_graphite_analyser16->_rotation_absolute, tr1);
    rot_mul(_graphite_analyser17->_rotation_absolute, tr1, _graphite_analyser17->_rotation_relative);
    _graphite_analyser17->_rotation_is_identity =  rot_test_identity(_graphite_analyser17->_rotation_relative);
    tc1 = coords_set(
      0, -0.025, 0.908542);
    rot_transpose(_armAnalyzer->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _graphite_analyser17->_position_absolute = coords_add(_armAnalyzer->_position_absolute, tc2);
    tc1 = coords_sub(_graphite_analyser16->_position_absolute, _graphite_analyser17->_position_absolute);
    _graphite_analyser17->_position_relative = rot_apply(_graphite_analyser17->_rotation_absolute, tc1);
  } /* graphite_analyser17=Monochromator_pol() AT ROTATED */
  DEBUG_COMPONENT("graphite_analyser17", _graphite_analyser17->_position_absolute, _graphite_analyser17->_rotation_absolute);
  instrument->_position_absolute[45] = _graphite_analyser17->_position_absolute;
  instrument->_position_relative[45] = _graphite_analyser17->_position_relative;
  instrument->counter_N[45]  = instrument->counter_P[45] = instrument->counter_P2[45] = 0;
  instrument->counter_AbsorbProp[45]= 0;
  return(0);
} /* _graphite_analyser17_setpos */

/* component graphite_analyser18=Monochromator_pol() SETTING, POSITION/ROTATION */
int _graphite_analyser18_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_graphite_analyser18_setpos] component graphite_analyser18=Monochromator_pol() SETTING [/usr/share/mcstas/3.0-dev/optics/Monochromator_pol.comp:109]");
  stracpy(_graphite_analyser18->_name, "graphite_analyser18", 16384);
  stracpy(_graphite_analyser18->_type, "Monochromator_pol", 16384);
  _graphite_analyser18->_index=46;
  _graphite_analyser18->_parameters.zwidth = ANAWIDTH;
  #define zwidth (_graphite_analyser18->_parameters.zwidth)
  _graphite_analyser18->_parameters.yheight = ANAHEIGHT;
  #define yheight (_graphite_analyser18->_parameters.yheight)
  _graphite_analyser18->_parameters.mosaic = ANAMOSAIC;
  #define mosaic (_graphite_analyser18->_parameters.mosaic)
  _graphite_analyser18->_parameters.dspread = 0.0025;
  #define dspread (_graphite_analyser18->_parameters.dspread)
  _graphite_analyser18->_parameters.Q = ANAQ;
  #define Q (_graphite_analyser18->_parameters.Q)
  _graphite_analyser18->_parameters.DM = 0;
  #define DM (_graphite_analyser18->_parameters.DM)
  _graphite_analyser18->_parameters.pThreshold = 0.001;
  #define pThreshold (_graphite_analyser18->_parameters.pThreshold)
  _graphite_analyser18->_parameters.Rup = 1.0;
  #define Rup (_graphite_analyser18->_parameters.Rup)
  _graphite_analyser18->_parameters.Rdown = 1.0;
  #define Rdown (_graphite_analyser18->_parameters.Rdown)
  _graphite_analyser18->_parameters.debug = 0;
  #define debug (_graphite_analyser18->_parameters.debug)

  #define mos_rms (_graphite_analyser18->_parameters.mos_rms)
  #define d_rms (_graphite_analyser18->_parameters.d_rms)
  #define mono_Q (_graphite_analyser18->_parameters.mono_Q)
  #define FN (_graphite_analyser18->_parameters.FN)
  #define FM (_graphite_analyser18->_parameters.FM)

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
  /* component graphite_analyser18=Monochromator_pol() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (90)*DEG2RAD, (-3.855)*DEG2RAD);
    rot_mul(tr1, _armAnalyzer->_rotation_absolute, _graphite_analyser18->_rotation_absolute);
    rot_transpose(_graphite_analyser17->_rotation_absolute, tr1);
    rot_mul(_graphite_analyser18->_rotation_absolute, tr1, _graphite_analyser18->_rotation_relative);
    _graphite_analyser18->_rotation_is_identity =  rot_test_identity(_graphite_analyser18->_rotation_relative);
    tc1 = coords_set(
      0, -0.015, 0.907878);
    rot_transpose(_armAnalyzer->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _graphite_analyser18->_position_absolute = coords_add(_armAnalyzer->_position_absolute, tc2);
    tc1 = coords_sub(_graphite_analyser17->_position_absolute, _graphite_analyser18->_position_absolute);
    _graphite_analyser18->_position_relative = rot_apply(_graphite_analyser18->_rotation_absolute, tc1);
  } /* graphite_analyser18=Monochromator_pol() AT ROTATED */
  DEBUG_COMPONENT("graphite_analyser18", _graphite_analyser18->_position_absolute, _graphite_analyser18->_rotation_absolute);
  instrument->_position_absolute[46] = _graphite_analyser18->_position_absolute;
  instrument->_position_relative[46] = _graphite_analyser18->_position_relative;
  instrument->counter_N[46]  = instrument->counter_P[46] = instrument->counter_P2[46] = 0;
  instrument->counter_AbsorbProp[46]= 0;
  return(0);
} /* _graphite_analyser18_setpos */

/* component graphite_analyser19=Monochromator_pol() SETTING, POSITION/ROTATION */
int _graphite_analyser19_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_graphite_analyser19_setpos] component graphite_analyser19=Monochromator_pol() SETTING [/usr/share/mcstas/3.0-dev/optics/Monochromator_pol.comp:109]");
  stracpy(_graphite_analyser19->_name, "graphite_analyser19", 16384);
  stracpy(_graphite_analyser19->_type, "Monochromator_pol", 16384);
  _graphite_analyser19->_index=47;
  _graphite_analyser19->_parameters.zwidth = ANAWIDTH;
  #define zwidth (_graphite_analyser19->_parameters.zwidth)
  _graphite_analyser19->_parameters.yheight = ANAHEIGHT;
  #define yheight (_graphite_analyser19->_parameters.yheight)
  _graphite_analyser19->_parameters.mosaic = ANAMOSAIC;
  #define mosaic (_graphite_analyser19->_parameters.mosaic)
  _graphite_analyser19->_parameters.dspread = 0.0025;
  #define dspread (_graphite_analyser19->_parameters.dspread)
  _graphite_analyser19->_parameters.Q = ANAQ;
  #define Q (_graphite_analyser19->_parameters.Q)
  _graphite_analyser19->_parameters.DM = 0;
  #define DM (_graphite_analyser19->_parameters.DM)
  _graphite_analyser19->_parameters.pThreshold = 0.001;
  #define pThreshold (_graphite_analyser19->_parameters.pThreshold)
  _graphite_analyser19->_parameters.Rup = 1.0;
  #define Rup (_graphite_analyser19->_parameters.Rup)
  _graphite_analyser19->_parameters.Rdown = 1.0;
  #define Rdown (_graphite_analyser19->_parameters.Rdown)
  _graphite_analyser19->_parameters.debug = 0;
  #define debug (_graphite_analyser19->_parameters.debug)

  #define mos_rms (_graphite_analyser19->_parameters.mos_rms)
  #define d_rms (_graphite_analyser19->_parameters.d_rms)
  #define mono_Q (_graphite_analyser19->_parameters.mono_Q)
  #define FN (_graphite_analyser19->_parameters.FN)
  #define FM (_graphite_analyser19->_parameters.FM)

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
  /* component graphite_analyser19=Monochromator_pol() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (90)*DEG2RAD, (-4.0525)*DEG2RAD);
    rot_mul(tr1, _armAnalyzer->_rotation_absolute, _graphite_analyser19->_rotation_absolute);
    rot_transpose(_graphite_analyser18->_rotation_absolute, tr1);
    rot_mul(_graphite_analyser19->_rotation_absolute, tr1, _graphite_analyser19->_rotation_relative);
    _graphite_analyser19->_rotation_is_identity =  rot_test_identity(_graphite_analyser19->_rotation_relative);
    tc1 = coords_set(
      0, -0.005, 0.907171);
    rot_transpose(_armAnalyzer->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _graphite_analyser19->_position_absolute = coords_add(_armAnalyzer->_position_absolute, tc2);
    tc1 = coords_sub(_graphite_analyser18->_position_absolute, _graphite_analyser19->_position_absolute);
    _graphite_analyser19->_position_relative = rot_apply(_graphite_analyser19->_rotation_absolute, tc1);
  } /* graphite_analyser19=Monochromator_pol() AT ROTATED */
  DEBUG_COMPONENT("graphite_analyser19", _graphite_analyser19->_position_absolute, _graphite_analyser19->_rotation_absolute);
  instrument->_position_absolute[47] = _graphite_analyser19->_position_absolute;
  instrument->_position_relative[47] = _graphite_analyser19->_position_relative;
  instrument->counter_N[47]  = instrument->counter_P[47] = instrument->counter_P2[47] = 0;
  instrument->counter_AbsorbProp[47]= 0;
  return(0);
} /* _graphite_analyser19_setpos */

/* component graphite_analyser20=Monochromator_pol() SETTING, POSITION/ROTATION */
int _graphite_analyser20_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_graphite_analyser20_setpos] component graphite_analyser20=Monochromator_pol() SETTING [/usr/share/mcstas/3.0-dev/optics/Monochromator_pol.comp:109]");
  stracpy(_graphite_analyser20->_name, "graphite_analyser20", 16384);
  stracpy(_graphite_analyser20->_type, "Monochromator_pol", 16384);
  _graphite_analyser20->_index=48;
  _graphite_analyser20->_parameters.zwidth = ANAWIDTH;
  #define zwidth (_graphite_analyser20->_parameters.zwidth)
  _graphite_analyser20->_parameters.yheight = ANAHEIGHT;
  #define yheight (_graphite_analyser20->_parameters.yheight)
  _graphite_analyser20->_parameters.mosaic = ANAMOSAIC;
  #define mosaic (_graphite_analyser20->_parameters.mosaic)
  _graphite_analyser20->_parameters.dspread = 0.0025;
  #define dspread (_graphite_analyser20->_parameters.dspread)
  _graphite_analyser20->_parameters.Q = ANAQ;
  #define Q (_graphite_analyser20->_parameters.Q)
  _graphite_analyser20->_parameters.DM = 0;
  #define DM (_graphite_analyser20->_parameters.DM)
  _graphite_analyser20->_parameters.pThreshold = 0.001;
  #define pThreshold (_graphite_analyser20->_parameters.pThreshold)
  _graphite_analyser20->_parameters.Rup = 1.0;
  #define Rup (_graphite_analyser20->_parameters.Rup)
  _graphite_analyser20->_parameters.Rdown = 1.0;
  #define Rdown (_graphite_analyser20->_parameters.Rdown)
  _graphite_analyser20->_parameters.debug = 0;
  #define debug (_graphite_analyser20->_parameters.debug)

  #define mos_rms (_graphite_analyser20->_parameters.mos_rms)
  #define d_rms (_graphite_analyser20->_parameters.d_rms)
  #define mono_Q (_graphite_analyser20->_parameters.mono_Q)
  #define FN (_graphite_analyser20->_parameters.FN)
  #define FM (_graphite_analyser20->_parameters.FM)

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
  /* component graphite_analyser20=Monochromator_pol() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (90)*DEG2RAD, (-5.1175)*DEG2RAD);
    rot_mul(tr1, _armAnalyzer->_rotation_absolute, _graphite_analyser20->_rotation_absolute);
    rot_transpose(_graphite_analyser19->_rotation_absolute, tr1);
    rot_mul(_graphite_analyser20->_rotation_absolute, tr1, _graphite_analyser20->_rotation_relative);
    _graphite_analyser20->_rotation_is_identity =  rot_test_identity(_graphite_analyser20->_rotation_relative);
    tc1 = coords_set(
      0, 0.005, 0.906285);
    rot_transpose(_armAnalyzer->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _graphite_analyser20->_position_absolute = coords_add(_armAnalyzer->_position_absolute, tc2);
    tc1 = coords_sub(_graphite_analyser19->_position_absolute, _graphite_analyser20->_position_absolute);
    _graphite_analyser20->_position_relative = rot_apply(_graphite_analyser20->_rotation_absolute, tc1);
  } /* graphite_analyser20=Monochromator_pol() AT ROTATED */
  DEBUG_COMPONENT("graphite_analyser20", _graphite_analyser20->_position_absolute, _graphite_analyser20->_rotation_absolute);
  instrument->_position_absolute[48] = _graphite_analyser20->_position_absolute;
  instrument->_position_relative[48] = _graphite_analyser20->_position_relative;
  instrument->counter_N[48]  = instrument->counter_P[48] = instrument->counter_P2[48] = 0;
  instrument->counter_AbsorbProp[48]= 0;
  return(0);
} /* _graphite_analyser20_setpos */

/* component graphite_analyser21=Monochromator_pol() SETTING, POSITION/ROTATION */
int _graphite_analyser21_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_graphite_analyser21_setpos] component graphite_analyser21=Monochromator_pol() SETTING [/usr/share/mcstas/3.0-dev/optics/Monochromator_pol.comp:109]");
  stracpy(_graphite_analyser21->_name, "graphite_analyser21", 16384);
  stracpy(_graphite_analyser21->_type, "Monochromator_pol", 16384);
  _graphite_analyser21->_index=49;
  _graphite_analyser21->_parameters.zwidth = ANAWIDTH;
  #define zwidth (_graphite_analyser21->_parameters.zwidth)
  _graphite_analyser21->_parameters.yheight = ANAHEIGHT;
  #define yheight (_graphite_analyser21->_parameters.yheight)
  _graphite_analyser21->_parameters.mosaic = ANAMOSAIC;
  #define mosaic (_graphite_analyser21->_parameters.mosaic)
  _graphite_analyser21->_parameters.dspread = 0.0025;
  #define dspread (_graphite_analyser21->_parameters.dspread)
  _graphite_analyser21->_parameters.Q = ANAQ;
  #define Q (_graphite_analyser21->_parameters.Q)
  _graphite_analyser21->_parameters.DM = 0;
  #define DM (_graphite_analyser21->_parameters.DM)
  _graphite_analyser21->_parameters.pThreshold = 0.001;
  #define pThreshold (_graphite_analyser21->_parameters.pThreshold)
  _graphite_analyser21->_parameters.Rup = 1.0;
  #define Rup (_graphite_analyser21->_parameters.Rup)
  _graphite_analyser21->_parameters.Rdown = 1.0;
  #define Rdown (_graphite_analyser21->_parameters.Rdown)
  _graphite_analyser21->_parameters.debug = 0;
  #define debug (_graphite_analyser21->_parameters.debug)

  #define mos_rms (_graphite_analyser21->_parameters.mos_rms)
  #define d_rms (_graphite_analyser21->_parameters.d_rms)
  #define mono_Q (_graphite_analyser21->_parameters.mono_Q)
  #define FN (_graphite_analyser21->_parameters.FN)
  #define FM (_graphite_analyser21->_parameters.FM)

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
  /* component graphite_analyser21=Monochromator_pol() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (90)*DEG2RAD, (-5.985)*DEG2RAD);
    rot_mul(tr1, _armAnalyzer->_rotation_absolute, _graphite_analyser21->_rotation_absolute);
    rot_transpose(_graphite_analyser20->_rotation_absolute, tr1);
    rot_mul(_graphite_analyser21->_rotation_absolute, tr1, _graphite_analyser21->_rotation_relative);
    _graphite_analyser21->_rotation_is_identity =  rot_test_identity(_graphite_analyser21->_rotation_relative);
    tc1 = coords_set(
      0, 0.015, 0.905246);
    rot_transpose(_armAnalyzer->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _graphite_analyser21->_position_absolute = coords_add(_armAnalyzer->_position_absolute, tc2);
    tc1 = coords_sub(_graphite_analyser20->_position_absolute, _graphite_analyser21->_position_absolute);
    _graphite_analyser21->_position_relative = rot_apply(_graphite_analyser21->_rotation_absolute, tc1);
  } /* graphite_analyser21=Monochromator_pol() AT ROTATED */
  DEBUG_COMPONENT("graphite_analyser21", _graphite_analyser21->_position_absolute, _graphite_analyser21->_rotation_absolute);
  instrument->_position_absolute[49] = _graphite_analyser21->_position_absolute;
  instrument->_position_relative[49] = _graphite_analyser21->_position_relative;
  instrument->counter_N[49]  = instrument->counter_P[49] = instrument->counter_P2[49] = 0;
  instrument->counter_AbsorbProp[49]= 0;
  return(0);
} /* _graphite_analyser21_setpos */

/* component graphite_analyser22=Monochromator_pol() SETTING, POSITION/ROTATION */
int _graphite_analyser22_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_graphite_analyser22_setpos] component graphite_analyser22=Monochromator_pol() SETTING [/usr/share/mcstas/3.0-dev/optics/Monochromator_pol.comp:109]");
  stracpy(_graphite_analyser22->_name, "graphite_analyser22", 16384);
  stracpy(_graphite_analyser22->_type, "Monochromator_pol", 16384);
  _graphite_analyser22->_index=50;
  _graphite_analyser22->_parameters.zwidth = ANAWIDTH;
  #define zwidth (_graphite_analyser22->_parameters.zwidth)
  _graphite_analyser22->_parameters.yheight = ANAHEIGHT;
  #define yheight (_graphite_analyser22->_parameters.yheight)
  _graphite_analyser22->_parameters.mosaic = ANAMOSAIC;
  #define mosaic (_graphite_analyser22->_parameters.mosaic)
  _graphite_analyser22->_parameters.dspread = 0.0025;
  #define dspread (_graphite_analyser22->_parameters.dspread)
  _graphite_analyser22->_parameters.Q = ANAQ;
  #define Q (_graphite_analyser22->_parameters.Q)
  _graphite_analyser22->_parameters.DM = 0;
  #define DM (_graphite_analyser22->_parameters.DM)
  _graphite_analyser22->_parameters.pThreshold = 0.001;
  #define pThreshold (_graphite_analyser22->_parameters.pThreshold)
  _graphite_analyser22->_parameters.Rup = 1.0;
  #define Rup (_graphite_analyser22->_parameters.Rup)
  _graphite_analyser22->_parameters.Rdown = 1.0;
  #define Rdown (_graphite_analyser22->_parameters.Rdown)
  _graphite_analyser22->_parameters.debug = 0;
  #define debug (_graphite_analyser22->_parameters.debug)

  #define mos_rms (_graphite_analyser22->_parameters.mos_rms)
  #define d_rms (_graphite_analyser22->_parameters.d_rms)
  #define mono_Q (_graphite_analyser22->_parameters.mono_Q)
  #define FN (_graphite_analyser22->_parameters.FN)
  #define FM (_graphite_analyser22->_parameters.FM)

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
  /* component graphite_analyser22=Monochromator_pol() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (90)*DEG2RAD, (-6.2575)*DEG2RAD);
    rot_mul(tr1, _armAnalyzer->_rotation_absolute, _graphite_analyser22->_rotation_absolute);
    rot_transpose(_graphite_analyser21->_rotation_absolute, tr1);
    rot_mul(_graphite_analyser22->_rotation_absolute, tr1, _graphite_analyser22->_rotation_relative);
    _graphite_analyser22->_rotation_is_identity =  rot_test_identity(_graphite_analyser22->_rotation_relative);
    tc1 = coords_set(
      0, 0.025, 0.904153);
    rot_transpose(_armAnalyzer->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _graphite_analyser22->_position_absolute = coords_add(_armAnalyzer->_position_absolute, tc2);
    tc1 = coords_sub(_graphite_analyser21->_position_absolute, _graphite_analyser22->_position_absolute);
    _graphite_analyser22->_position_relative = rot_apply(_graphite_analyser22->_rotation_absolute, tc1);
  } /* graphite_analyser22=Monochromator_pol() AT ROTATED */
  DEBUG_COMPONENT("graphite_analyser22", _graphite_analyser22->_position_absolute, _graphite_analyser22->_rotation_absolute);
  instrument->_position_absolute[50] = _graphite_analyser22->_position_absolute;
  instrument->_position_relative[50] = _graphite_analyser22->_position_relative;
  instrument->counter_N[50]  = instrument->counter_P[50] = instrument->counter_P2[50] = 0;
  instrument->counter_AbsorbProp[50]= 0;
  return(0);
} /* _graphite_analyser22_setpos */

/* component graphite_analyser23=Monochromator_pol() SETTING, POSITION/ROTATION */
int _graphite_analyser23_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_graphite_analyser23_setpos] component graphite_analyser23=Monochromator_pol() SETTING [/usr/share/mcstas/3.0-dev/optics/Monochromator_pol.comp:109]");
  stracpy(_graphite_analyser23->_name, "graphite_analyser23", 16384);
  stracpy(_graphite_analyser23->_type, "Monochromator_pol", 16384);
  _graphite_analyser23->_index=51;
  _graphite_analyser23->_parameters.zwidth = ANAWIDTH;
  #define zwidth (_graphite_analyser23->_parameters.zwidth)
  _graphite_analyser23->_parameters.yheight = ANAHEIGHT;
  #define yheight (_graphite_analyser23->_parameters.yheight)
  _graphite_analyser23->_parameters.mosaic = ANAMOSAIC;
  #define mosaic (_graphite_analyser23->_parameters.mosaic)
  _graphite_analyser23->_parameters.dspread = 0.0025;
  #define dspread (_graphite_analyser23->_parameters.dspread)
  _graphite_analyser23->_parameters.Q = ANAQ;
  #define Q (_graphite_analyser23->_parameters.Q)
  _graphite_analyser23->_parameters.DM = 0;
  #define DM (_graphite_analyser23->_parameters.DM)
  _graphite_analyser23->_parameters.pThreshold = 0.001;
  #define pThreshold (_graphite_analyser23->_parameters.pThreshold)
  _graphite_analyser23->_parameters.Rup = 1.0;
  #define Rup (_graphite_analyser23->_parameters.Rup)
  _graphite_analyser23->_parameters.Rdown = 1.0;
  #define Rdown (_graphite_analyser23->_parameters.Rdown)
  _graphite_analyser23->_parameters.debug = 0;
  #define debug (_graphite_analyser23->_parameters.debug)

  #define mos_rms (_graphite_analyser23->_parameters.mos_rms)
  #define d_rms (_graphite_analyser23->_parameters.d_rms)
  #define mono_Q (_graphite_analyser23->_parameters.mono_Q)
  #define FN (_graphite_analyser23->_parameters.FN)
  #define FM (_graphite_analyser23->_parameters.FM)

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
  /* component graphite_analyser23=Monochromator_pol() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (90)*DEG2RAD, (-7.325)*DEG2RAD);
    rot_mul(tr1, _armAnalyzer->_rotation_absolute, _graphite_analyser23->_rotation_absolute);
    rot_transpose(_graphite_analyser22->_rotation_absolute, tr1);
    rot_mul(_graphite_analyser23->_rotation_absolute, tr1, _graphite_analyser23->_rotation_relative);
    _graphite_analyser23->_rotation_is_identity =  rot_test_identity(_graphite_analyser23->_rotation_relative);
    tc1 = coords_set(
      0, 0.035, 0.902877);
    rot_transpose(_armAnalyzer->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _graphite_analyser23->_position_absolute = coords_add(_armAnalyzer->_position_absolute, tc2);
    tc1 = coords_sub(_graphite_analyser22->_position_absolute, _graphite_analyser23->_position_absolute);
    _graphite_analyser23->_position_relative = rot_apply(_graphite_analyser23->_rotation_absolute, tc1);
  } /* graphite_analyser23=Monochromator_pol() AT ROTATED */
  DEBUG_COMPONENT("graphite_analyser23", _graphite_analyser23->_position_absolute, _graphite_analyser23->_rotation_absolute);
  instrument->_position_absolute[51] = _graphite_analyser23->_position_absolute;
  instrument->_position_relative[51] = _graphite_analyser23->_position_relative;
  instrument->counter_N[51]  = instrument->counter_P[51] = instrument->counter_P2[51] = 0;
  instrument->counter_AbsorbProp[51]= 0;
  return(0);
} /* _graphite_analyser23_setpos */

/* component graphite_analyser24=Monochromator_pol() SETTING, POSITION/ROTATION */
int _graphite_analyser24_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_graphite_analyser24_setpos] component graphite_analyser24=Monochromator_pol() SETTING [/usr/share/mcstas/3.0-dev/optics/Monochromator_pol.comp:109]");
  stracpy(_graphite_analyser24->_name, "graphite_analyser24", 16384);
  stracpy(_graphite_analyser24->_type, "Monochromator_pol", 16384);
  _graphite_analyser24->_index=52;
  _graphite_analyser24->_parameters.zwidth = ANAWIDTH;
  #define zwidth (_graphite_analyser24->_parameters.zwidth)
  _graphite_analyser24->_parameters.yheight = ANAHEIGHT;
  #define yheight (_graphite_analyser24->_parameters.yheight)
  _graphite_analyser24->_parameters.mosaic = ANAMOSAIC;
  #define mosaic (_graphite_analyser24->_parameters.mosaic)
  _graphite_analyser24->_parameters.dspread = 0.0025;
  #define dspread (_graphite_analyser24->_parameters.dspread)
  _graphite_analyser24->_parameters.Q = ANAQ;
  #define Q (_graphite_analyser24->_parameters.Q)
  _graphite_analyser24->_parameters.DM = 0;
  #define DM (_graphite_analyser24->_parameters.DM)
  _graphite_analyser24->_parameters.pThreshold = 0.001;
  #define pThreshold (_graphite_analyser24->_parameters.pThreshold)
  _graphite_analyser24->_parameters.Rup = 1.0;
  #define Rup (_graphite_analyser24->_parameters.Rup)
  _graphite_analyser24->_parameters.Rdown = 1.0;
  #define Rdown (_graphite_analyser24->_parameters.Rdown)
  _graphite_analyser24->_parameters.debug = 0;
  #define debug (_graphite_analyser24->_parameters.debug)

  #define mos_rms (_graphite_analyser24->_parameters.mos_rms)
  #define d_rms (_graphite_analyser24->_parameters.d_rms)
  #define mono_Q (_graphite_analyser24->_parameters.mono_Q)
  #define FN (_graphite_analyser24->_parameters.FN)
  #define FM (_graphite_analyser24->_parameters.FM)

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
  /* component graphite_analyser24=Monochromator_pol() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (90)*DEG2RAD, (-8.12)*DEG2RAD);
    rot_mul(tr1, _armAnalyzer->_rotation_absolute, _graphite_analyser24->_rotation_absolute);
    rot_transpose(_graphite_analyser23->_rotation_absolute, tr1);
    rot_mul(_graphite_analyser24->_rotation_absolute, tr1, _graphite_analyser24->_rotation_relative);
    _graphite_analyser24->_rotation_is_identity =  rot_test_identity(_graphite_analyser24->_rotation_relative);
    tc1 = coords_set(
      0, 0.045, 0.901459);
    rot_transpose(_armAnalyzer->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _graphite_analyser24->_position_absolute = coords_add(_armAnalyzer->_position_absolute, tc2);
    tc1 = coords_sub(_graphite_analyser23->_position_absolute, _graphite_analyser24->_position_absolute);
    _graphite_analyser24->_position_relative = rot_apply(_graphite_analyser24->_rotation_absolute, tc1);
  } /* graphite_analyser24=Monochromator_pol() AT ROTATED */
  DEBUG_COMPONENT("graphite_analyser24", _graphite_analyser24->_position_absolute, _graphite_analyser24->_rotation_absolute);
  instrument->_position_absolute[52] = _graphite_analyser24->_position_absolute;
  instrument->_position_relative[52] = _graphite_analyser24->_position_relative;
  instrument->counter_N[52]  = instrument->counter_P[52] = instrument->counter_P2[52] = 0;
  instrument->counter_AbsorbProp[52]= 0;
  return(0);
} /* _graphite_analyser24_setpos */

/* component graphite_analyser25=Monochromator_pol() SETTING, POSITION/ROTATION */
int _graphite_analyser25_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_graphite_analyser25_setpos] component graphite_analyser25=Monochromator_pol() SETTING [/usr/share/mcstas/3.0-dev/optics/Monochromator_pol.comp:109]");
  stracpy(_graphite_analyser25->_name, "graphite_analyser25", 16384);
  stracpy(_graphite_analyser25->_type, "Monochromator_pol", 16384);
  _graphite_analyser25->_index=53;
  _graphite_analyser25->_parameters.zwidth = ANAWIDTH;
  #define zwidth (_graphite_analyser25->_parameters.zwidth)
  _graphite_analyser25->_parameters.yheight = ANAHEIGHT;
  #define yheight (_graphite_analyser25->_parameters.yheight)
  _graphite_analyser25->_parameters.mosaic = ANAMOSAIC;
  #define mosaic (_graphite_analyser25->_parameters.mosaic)
  _graphite_analyser25->_parameters.dspread = 0.0025;
  #define dspread (_graphite_analyser25->_parameters.dspread)
  _graphite_analyser25->_parameters.Q = ANAQ;
  #define Q (_graphite_analyser25->_parameters.Q)
  _graphite_analyser25->_parameters.DM = 0;
  #define DM (_graphite_analyser25->_parameters.DM)
  _graphite_analyser25->_parameters.pThreshold = 0.001;
  #define pThreshold (_graphite_analyser25->_parameters.pThreshold)
  _graphite_analyser25->_parameters.Rup = 1.0;
  #define Rup (_graphite_analyser25->_parameters.Rup)
  _graphite_analyser25->_parameters.Rdown = 1.0;
  #define Rdown (_graphite_analyser25->_parameters.Rdown)
  _graphite_analyser25->_parameters.debug = 0;
  #define debug (_graphite_analyser25->_parameters.debug)

  #define mos_rms (_graphite_analyser25->_parameters.mos_rms)
  #define d_rms (_graphite_analyser25->_parameters.d_rms)
  #define mono_Q (_graphite_analyser25->_parameters.mono_Q)
  #define FN (_graphite_analyser25->_parameters.FN)
  #define FM (_graphite_analyser25->_parameters.FM)

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
  /* component graphite_analyser25=Monochromator_pol() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (90)*DEG2RAD, (-8.4745)*DEG2RAD);
    rot_mul(tr1, _armAnalyzer->_rotation_absolute, _graphite_analyser25->_rotation_absolute);
    rot_transpose(_graphite_analyser24->_rotation_absolute, tr1);
    rot_mul(_graphite_analyser25->_rotation_absolute, tr1, _graphite_analyser25->_rotation_relative);
    _graphite_analyser25->_rotation_is_identity =  rot_test_identity(_graphite_analyser25->_rotation_relative);
    tc1 = coords_set(
      0, 0.055, 0.899973);
    rot_transpose(_armAnalyzer->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _graphite_analyser25->_position_absolute = coords_add(_armAnalyzer->_position_absolute, tc2);
    tc1 = coords_sub(_graphite_analyser24->_position_absolute, _graphite_analyser25->_position_absolute);
    _graphite_analyser25->_position_relative = rot_apply(_graphite_analyser25->_rotation_absolute, tc1);
  } /* graphite_analyser25=Monochromator_pol() AT ROTATED */
  DEBUG_COMPONENT("graphite_analyser25", _graphite_analyser25->_position_absolute, _graphite_analyser25->_rotation_absolute);
  instrument->_position_absolute[53] = _graphite_analyser25->_position_absolute;
  instrument->_position_relative[53] = _graphite_analyser25->_position_relative;
  instrument->counter_N[53]  = instrument->counter_P[53] = instrument->counter_P2[53] = 0;
  instrument->counter_AbsorbProp[53]= 0;
  return(0);
} /* _graphite_analyser25_setpos */

/* component graphite_analyser26=Monochromator_pol() SETTING, POSITION/ROTATION */
int _graphite_analyser26_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_graphite_analyser26_setpos] component graphite_analyser26=Monochromator_pol() SETTING [/usr/share/mcstas/3.0-dev/optics/Monochromator_pol.comp:109]");
  stracpy(_graphite_analyser26->_name, "graphite_analyser26", 16384);
  stracpy(_graphite_analyser26->_type, "Monochromator_pol", 16384);
  _graphite_analyser26->_index=54;
  _graphite_analyser26->_parameters.zwidth = ANAWIDTH;
  #define zwidth (_graphite_analyser26->_parameters.zwidth)
  _graphite_analyser26->_parameters.yheight = ANAHEIGHT;
  #define yheight (_graphite_analyser26->_parameters.yheight)
  _graphite_analyser26->_parameters.mosaic = ANAMOSAIC;
  #define mosaic (_graphite_analyser26->_parameters.mosaic)
  _graphite_analyser26->_parameters.dspread = 0.0025;
  #define dspread (_graphite_analyser26->_parameters.dspread)
  _graphite_analyser26->_parameters.Q = ANAQ;
  #define Q (_graphite_analyser26->_parameters.Q)
  _graphite_analyser26->_parameters.DM = 0;
  #define DM (_graphite_analyser26->_parameters.DM)
  _graphite_analyser26->_parameters.pThreshold = 0.001;
  #define pThreshold (_graphite_analyser26->_parameters.pThreshold)
  _graphite_analyser26->_parameters.Rup = 1.0;
  #define Rup (_graphite_analyser26->_parameters.Rup)
  _graphite_analyser26->_parameters.Rdown = 1.0;
  #define Rdown (_graphite_analyser26->_parameters.Rdown)
  _graphite_analyser26->_parameters.debug = 0;
  #define debug (_graphite_analyser26->_parameters.debug)

  #define mos_rms (_graphite_analyser26->_parameters.mos_rms)
  #define d_rms (_graphite_analyser26->_parameters.d_rms)
  #define mono_Q (_graphite_analyser26->_parameters.mono_Q)
  #define FN (_graphite_analyser26->_parameters.FN)
  #define FM (_graphite_analyser26->_parameters.FM)

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
  /* component graphite_analyser26=Monochromator_pol() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (90)*DEG2RAD, (-9.5425)*DEG2RAD);
    rot_mul(tr1, _armAnalyzer->_rotation_absolute, _graphite_analyser26->_rotation_absolute);
    rot_transpose(_graphite_analyser25->_rotation_absolute, tr1);
    rot_mul(_graphite_analyser26->_rotation_absolute, tr1, _graphite_analyser26->_rotation_relative);
    _graphite_analyser26->_rotation_is_identity =  rot_test_identity(_graphite_analyser26->_rotation_relative);
    tc1 = coords_set(
      0, 0.065, 0.898303);
    rot_transpose(_armAnalyzer->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _graphite_analyser26->_position_absolute = coords_add(_armAnalyzer->_position_absolute, tc2);
    tc1 = coords_sub(_graphite_analyser25->_position_absolute, _graphite_analyser26->_position_absolute);
    _graphite_analyser26->_position_relative = rot_apply(_graphite_analyser26->_rotation_absolute, tc1);
  } /* graphite_analyser26=Monochromator_pol() AT ROTATED */
  DEBUG_COMPONENT("graphite_analyser26", _graphite_analyser26->_position_absolute, _graphite_analyser26->_rotation_absolute);
  instrument->_position_absolute[54] = _graphite_analyser26->_position_absolute;
  instrument->_position_relative[54] = _graphite_analyser26->_position_relative;
  instrument->counter_N[54]  = instrument->counter_P[54] = instrument->counter_P2[54] = 0;
  instrument->counter_AbsorbProp[54]= 0;
  return(0);
} /* _graphite_analyser26_setpos */

/* component graphite_analyser27=Monochromator_pol() SETTING, POSITION/ROTATION */
int _graphite_analyser27_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_graphite_analyser27_setpos] component graphite_analyser27=Monochromator_pol() SETTING [/usr/share/mcstas/3.0-dev/optics/Monochromator_pol.comp:109]");
  stracpy(_graphite_analyser27->_name, "graphite_analyser27", 16384);
  stracpy(_graphite_analyser27->_type, "Monochromator_pol", 16384);
  _graphite_analyser27->_index=55;
  _graphite_analyser27->_parameters.zwidth = ANAWIDTH;
  #define zwidth (_graphite_analyser27->_parameters.zwidth)
  _graphite_analyser27->_parameters.yheight = ANAHEIGHT;
  #define yheight (_graphite_analyser27->_parameters.yheight)
  _graphite_analyser27->_parameters.mosaic = ANAMOSAIC;
  #define mosaic (_graphite_analyser27->_parameters.mosaic)
  _graphite_analyser27->_parameters.dspread = 0.0025;
  #define dspread (_graphite_analyser27->_parameters.dspread)
  _graphite_analyser27->_parameters.Q = ANAQ;
  #define Q (_graphite_analyser27->_parameters.Q)
  _graphite_analyser27->_parameters.DM = 0;
  #define DM (_graphite_analyser27->_parameters.DM)
  _graphite_analyser27->_parameters.pThreshold = 0.001;
  #define pThreshold (_graphite_analyser27->_parameters.pThreshold)
  _graphite_analyser27->_parameters.Rup = 1.0;
  #define Rup (_graphite_analyser27->_parameters.Rup)
  _graphite_analyser27->_parameters.Rdown = 1.0;
  #define Rdown (_graphite_analyser27->_parameters.Rdown)
  _graphite_analyser27->_parameters.debug = 0;
  #define debug (_graphite_analyser27->_parameters.debug)

  #define mos_rms (_graphite_analyser27->_parameters.mos_rms)
  #define d_rms (_graphite_analyser27->_parameters.d_rms)
  #define mono_Q (_graphite_analyser27->_parameters.mono_Q)
  #define FN (_graphite_analyser27->_parameters.FN)
  #define FM (_graphite_analyser27->_parameters.FM)

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
  /* component graphite_analyser27=Monochromator_pol() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (90)*DEG2RAD, (-10.256)*DEG2RAD);
    rot_mul(tr1, _armAnalyzer->_rotation_absolute, _graphite_analyser27->_rotation_absolute);
    rot_transpose(_graphite_analyser26->_rotation_absolute, tr1);
    rot_mul(_graphite_analyser27->_rotation_absolute, tr1, _graphite_analyser27->_rotation_relative);
    _graphite_analyser27->_rotation_is_identity =  rot_test_identity(_graphite_analyser27->_rotation_relative);
    tc1 = coords_set(
      0, 0.075, 0.8965);
    rot_transpose(_armAnalyzer->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _graphite_analyser27->_position_absolute = coords_add(_armAnalyzer->_position_absolute, tc2);
    tc1 = coords_sub(_graphite_analyser26->_position_absolute, _graphite_analyser27->_position_absolute);
    _graphite_analyser27->_position_relative = rot_apply(_graphite_analyser27->_rotation_absolute, tc1);
  } /* graphite_analyser27=Monochromator_pol() AT ROTATED */
  DEBUG_COMPONENT("graphite_analyser27", _graphite_analyser27->_position_absolute, _graphite_analyser27->_rotation_absolute);
  instrument->_position_absolute[55] = _graphite_analyser27->_position_absolute;
  instrument->_position_relative[55] = _graphite_analyser27->_position_relative;
  instrument->counter_N[55]  = instrument->counter_P[55] = instrument->counter_P2[55] = 0;
  instrument->counter_AbsorbProp[55]= 0;
  return(0);
} /* _graphite_analyser27_setpos */

/* component graphite_analyser28=Monochromator_pol() SETTING, POSITION/ROTATION */
int _graphite_analyser28_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_graphite_analyser28_setpos] component graphite_analyser28=Monochromator_pol() SETTING [/usr/share/mcstas/3.0-dev/optics/Monochromator_pol.comp:109]");
  stracpy(_graphite_analyser28->_name, "graphite_analyser28", 16384);
  stracpy(_graphite_analyser28->_type, "Monochromator_pol", 16384);
  _graphite_analyser28->_index=56;
  _graphite_analyser28->_parameters.zwidth = ANAWIDTH;
  #define zwidth (_graphite_analyser28->_parameters.zwidth)
  _graphite_analyser28->_parameters.yheight = ANAHEIGHT;
  #define yheight (_graphite_analyser28->_parameters.yheight)
  _graphite_analyser28->_parameters.mosaic = ANAMOSAIC;
  #define mosaic (_graphite_analyser28->_parameters.mosaic)
  _graphite_analyser28->_parameters.dspread = 0.0025;
  #define dspread (_graphite_analyser28->_parameters.dspread)
  _graphite_analyser28->_parameters.Q = ANAQ;
  #define Q (_graphite_analyser28->_parameters.Q)
  _graphite_analyser28->_parameters.DM = 0;
  #define DM (_graphite_analyser28->_parameters.DM)
  _graphite_analyser28->_parameters.pThreshold = 0.001;
  #define pThreshold (_graphite_analyser28->_parameters.pThreshold)
  _graphite_analyser28->_parameters.Rup = 1.0;
  #define Rup (_graphite_analyser28->_parameters.Rup)
  _graphite_analyser28->_parameters.Rdown = 1.0;
  #define Rdown (_graphite_analyser28->_parameters.Rdown)
  _graphite_analyser28->_parameters.debug = 0;
  #define debug (_graphite_analyser28->_parameters.debug)

  #define mos_rms (_graphite_analyser28->_parameters.mos_rms)
  #define d_rms (_graphite_analyser28->_parameters.d_rms)
  #define mono_Q (_graphite_analyser28->_parameters.mono_Q)
  #define FN (_graphite_analyser28->_parameters.FN)
  #define FM (_graphite_analyser28->_parameters.FM)

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
  /* component graphite_analyser28=Monochromator_pol() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (90)*DEG2RAD, (-10.704)*DEG2RAD);
    rot_mul(tr1, _armAnalyzer->_rotation_absolute, _graphite_analyser28->_rotation_absolute);
    rot_transpose(_graphite_analyser27->_rotation_absolute, tr1);
    rot_mul(_graphite_analyser28->_rotation_absolute, tr1, _graphite_analyser28->_rotation_relative);
    _graphite_analyser28->_rotation_is_identity =  rot_test_identity(_graphite_analyser28->_rotation_relative);
    tc1 = coords_set(
      0, 0.085, 0.894614);
    rot_transpose(_armAnalyzer->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _graphite_analyser28->_position_absolute = coords_add(_armAnalyzer->_position_absolute, tc2);
    tc1 = coords_sub(_graphite_analyser27->_position_absolute, _graphite_analyser28->_position_absolute);
    _graphite_analyser28->_position_relative = rot_apply(_graphite_analyser28->_rotation_absolute, tc1);
  } /* graphite_analyser28=Monochromator_pol() AT ROTATED */
  DEBUG_COMPONENT("graphite_analyser28", _graphite_analyser28->_position_absolute, _graphite_analyser28->_rotation_absolute);
  instrument->_position_absolute[56] = _graphite_analyser28->_position_absolute;
  instrument->_position_relative[56] = _graphite_analyser28->_position_relative;
  instrument->counter_N[56]  = instrument->counter_P[56] = instrument->counter_P2[56] = 0;
  instrument->counter_AbsorbProp[56]= 0;
  return(0);
} /* _graphite_analyser28_setpos */

/* component graphite_analyser29=Monochromator_pol() SETTING, POSITION/ROTATION */
int _graphite_analyser29_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_graphite_analyser29_setpos] component graphite_analyser29=Monochromator_pol() SETTING [/usr/share/mcstas/3.0-dev/optics/Monochromator_pol.comp:109]");
  stracpy(_graphite_analyser29->_name, "graphite_analyser29", 16384);
  stracpy(_graphite_analyser29->_type, "Monochromator_pol", 16384);
  _graphite_analyser29->_index=57;
  _graphite_analyser29->_parameters.zwidth = ANAWIDTH;
  #define zwidth (_graphite_analyser29->_parameters.zwidth)
  _graphite_analyser29->_parameters.yheight = ANAHEIGHT;
  #define yheight (_graphite_analyser29->_parameters.yheight)
  _graphite_analyser29->_parameters.mosaic = ANAMOSAIC;
  #define mosaic (_graphite_analyser29->_parameters.mosaic)
  _graphite_analyser29->_parameters.dspread = 0.0025;
  #define dspread (_graphite_analyser29->_parameters.dspread)
  _graphite_analyser29->_parameters.Q = ANAQ;
  #define Q (_graphite_analyser29->_parameters.Q)
  _graphite_analyser29->_parameters.DM = 0;
  #define DM (_graphite_analyser29->_parameters.DM)
  _graphite_analyser29->_parameters.pThreshold = 0.001;
  #define pThreshold (_graphite_analyser29->_parameters.pThreshold)
  _graphite_analyser29->_parameters.Rup = 1.0;
  #define Rup (_graphite_analyser29->_parameters.Rup)
  _graphite_analyser29->_parameters.Rdown = 1.0;
  #define Rdown (_graphite_analyser29->_parameters.Rdown)
  _graphite_analyser29->_parameters.debug = 0;
  #define debug (_graphite_analyser29->_parameters.debug)

  #define mos_rms (_graphite_analyser29->_parameters.mos_rms)
  #define d_rms (_graphite_analyser29->_parameters.d_rms)
  #define mono_Q (_graphite_analyser29->_parameters.mono_Q)
  #define FN (_graphite_analyser29->_parameters.FN)
  #define FM (_graphite_analyser29->_parameters.FM)

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
  /* component graphite_analyser29=Monochromator_pol() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (90)*DEG2RAD, (-11.768)*DEG2RAD);
    rot_mul(tr1, _armAnalyzer->_rotation_absolute, _graphite_analyser29->_rotation_absolute);
    rot_transpose(_graphite_analyser28->_rotation_absolute, tr1);
    rot_mul(_graphite_analyser29->_rotation_absolute, tr1, _graphite_analyser29->_rotation_relative);
    _graphite_analyser29->_rotation_is_identity =  rot_test_identity(_graphite_analyser29->_rotation_relative);
    tc1 = coords_set(
      0, 0.095, 0.892542);
    rot_transpose(_armAnalyzer->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _graphite_analyser29->_position_absolute = coords_add(_armAnalyzer->_position_absolute, tc2);
    tc1 = coords_sub(_graphite_analyser28->_position_absolute, _graphite_analyser29->_position_absolute);
    _graphite_analyser29->_position_relative = rot_apply(_graphite_analyser29->_rotation_absolute, tc1);
  } /* graphite_analyser29=Monochromator_pol() AT ROTATED */
  DEBUG_COMPONENT("graphite_analyser29", _graphite_analyser29->_position_absolute, _graphite_analyser29->_rotation_absolute);
  instrument->_position_absolute[57] = _graphite_analyser29->_position_absolute;
  instrument->_position_relative[57] = _graphite_analyser29->_position_relative;
  instrument->counter_N[57]  = instrument->counter_P[57] = instrument->counter_P2[57] = 0;
  instrument->counter_AbsorbProp[57]= 0;
  return(0);
} /* _graphite_analyser29_setpos */

/* component graphite_analyser30=Monochromator_pol() SETTING, POSITION/ROTATION */
int _graphite_analyser30_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_graphite_analyser30_setpos] component graphite_analyser30=Monochromator_pol() SETTING [/usr/share/mcstas/3.0-dev/optics/Monochromator_pol.comp:109]");
  stracpy(_graphite_analyser30->_name, "graphite_analyser30", 16384);
  stracpy(_graphite_analyser30->_type, "Monochromator_pol", 16384);
  _graphite_analyser30->_index=58;
  _graphite_analyser30->_parameters.zwidth = ANAWIDTH;
  #define zwidth (_graphite_analyser30->_parameters.zwidth)
  _graphite_analyser30->_parameters.yheight = ANAHEIGHT;
  #define yheight (_graphite_analyser30->_parameters.yheight)
  _graphite_analyser30->_parameters.mosaic = ANAMOSAIC;
  #define mosaic (_graphite_analyser30->_parameters.mosaic)
  _graphite_analyser30->_parameters.dspread = 0.0025;
  #define dspread (_graphite_analyser30->_parameters.dspread)
  _graphite_analyser30->_parameters.Q = ANAQ;
  #define Q (_graphite_analyser30->_parameters.Q)
  _graphite_analyser30->_parameters.DM = 0;
  #define DM (_graphite_analyser30->_parameters.DM)
  _graphite_analyser30->_parameters.pThreshold = 0.001;
  #define pThreshold (_graphite_analyser30->_parameters.pThreshold)
  _graphite_analyser30->_parameters.Rup = 1.0;
  #define Rup (_graphite_analyser30->_parameters.Rup)
  _graphite_analyser30->_parameters.Rdown = 1.0;
  #define Rdown (_graphite_analyser30->_parameters.Rdown)
  _graphite_analyser30->_parameters.debug = 0;
  #define debug (_graphite_analyser30->_parameters.debug)

  #define mos_rms (_graphite_analyser30->_parameters.mos_rms)
  #define d_rms (_graphite_analyser30->_parameters.d_rms)
  #define mono_Q (_graphite_analyser30->_parameters.mono_Q)
  #define FN (_graphite_analyser30->_parameters.FN)
  #define FM (_graphite_analyser30->_parameters.FM)

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
  /* component graphite_analyser30=Monochromator_pol() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (90)*DEG2RAD, (-12.384)*DEG2RAD);
    rot_mul(tr1, _armAnalyzer->_rotation_absolute, _graphite_analyser30->_rotation_absolute);
    rot_transpose(_graphite_analyser29->_rotation_absolute, tr1);
    rot_mul(_graphite_analyser30->_rotation_absolute, tr1, _graphite_analyser30->_rotation_relative);
    _graphite_analyser30->_rotation_is_identity =  rot_test_identity(_graphite_analyser30->_rotation_relative);
    tc1 = coords_set(
      0, 0.105, 0.890353);
    rot_transpose(_armAnalyzer->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _graphite_analyser30->_position_absolute = coords_add(_armAnalyzer->_position_absolute, tc2);
    tc1 = coords_sub(_graphite_analyser29->_position_absolute, _graphite_analyser30->_position_absolute);
    _graphite_analyser30->_position_relative = rot_apply(_graphite_analyser30->_rotation_absolute, tc1);
  } /* graphite_analyser30=Monochromator_pol() AT ROTATED */
  DEBUG_COMPONENT("graphite_analyser30", _graphite_analyser30->_position_absolute, _graphite_analyser30->_rotation_absolute);
  instrument->_position_absolute[58] = _graphite_analyser30->_position_absolute;
  instrument->_position_relative[58] = _graphite_analyser30->_position_relative;
  instrument->counter_N[58]  = instrument->counter_P[58] = instrument->counter_P2[58] = 0;
  instrument->counter_AbsorbProp[58]= 0;
  return(0);
} /* _graphite_analyser30_setpos */

/* component graphite_analyser31=Monochromator_pol() SETTING, POSITION/ROTATION */
int _graphite_analyser31_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_graphite_analyser31_setpos] component graphite_analyser31=Monochromator_pol() SETTING [/usr/share/mcstas/3.0-dev/optics/Monochromator_pol.comp:109]");
  stracpy(_graphite_analyser31->_name, "graphite_analyser31", 16384);
  stracpy(_graphite_analyser31->_type, "Monochromator_pol", 16384);
  _graphite_analyser31->_index=59;
  _graphite_analyser31->_parameters.zwidth = ANAWIDTH;
  #define zwidth (_graphite_analyser31->_parameters.zwidth)
  _graphite_analyser31->_parameters.yheight = ANAHEIGHT;
  #define yheight (_graphite_analyser31->_parameters.yheight)
  _graphite_analyser31->_parameters.mosaic = ANAMOSAIC;
  #define mosaic (_graphite_analyser31->_parameters.mosaic)
  _graphite_analyser31->_parameters.dspread = 0.0025;
  #define dspread (_graphite_analyser31->_parameters.dspread)
  _graphite_analyser31->_parameters.Q = ANAQ;
  #define Q (_graphite_analyser31->_parameters.Q)
  _graphite_analyser31->_parameters.DM = 0;
  #define DM (_graphite_analyser31->_parameters.DM)
  _graphite_analyser31->_parameters.pThreshold = 0.001;
  #define pThreshold (_graphite_analyser31->_parameters.pThreshold)
  _graphite_analyser31->_parameters.Rup = 1.0;
  #define Rup (_graphite_analyser31->_parameters.Rup)
  _graphite_analyser31->_parameters.Rdown = 1.0;
  #define Rdown (_graphite_analyser31->_parameters.Rdown)
  _graphite_analyser31->_parameters.debug = 0;
  #define debug (_graphite_analyser31->_parameters.debug)

  #define mos_rms (_graphite_analyser31->_parameters.mos_rms)
  #define d_rms (_graphite_analyser31->_parameters.d_rms)
  #define mono_Q (_graphite_analyser31->_parameters.mono_Q)
  #define FN (_graphite_analyser31->_parameters.FN)
  #define FM (_graphite_analyser31->_parameters.FM)

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
  /* component graphite_analyser31=Monochromator_pol() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (90)*DEG2RAD, (-12.948)*DEG2RAD);
    rot_mul(tr1, _armAnalyzer->_rotation_absolute, _graphite_analyser31->_rotation_absolute);
    rot_transpose(_graphite_analyser30->_rotation_absolute, tr1);
    rot_mul(_graphite_analyser31->_rotation_absolute, tr1, _graphite_analyser31->_rotation_relative);
    _graphite_analyser31->_rotation_is_identity =  rot_test_identity(_graphite_analyser31->_rotation_relative);
    tc1 = coords_set(
      0, 0.115, 0.88806);
    rot_transpose(_armAnalyzer->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _graphite_analyser31->_position_absolute = coords_add(_armAnalyzer->_position_absolute, tc2);
    tc1 = coords_sub(_graphite_analyser30->_position_absolute, _graphite_analyser31->_position_absolute);
    _graphite_analyser31->_position_relative = rot_apply(_graphite_analyser31->_rotation_absolute, tc1);
  } /* graphite_analyser31=Monochromator_pol() AT ROTATED */
  DEBUG_COMPONENT("graphite_analyser31", _graphite_analyser31->_position_absolute, _graphite_analyser31->_rotation_absolute);
  instrument->_position_absolute[59] = _graphite_analyser31->_position_absolute;
  instrument->_position_relative[59] = _graphite_analyser31->_position_relative;
  instrument->counter_N[59]  = instrument->counter_P[59] = instrument->counter_P2[59] = 0;
  instrument->counter_AbsorbProp[59]= 0;
  return(0);
} /* _graphite_analyser31_setpos */

/* component graphite_analyser32=Monochromator_pol() SETTING, POSITION/ROTATION */
int _graphite_analyser32_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_graphite_analyser32_setpos] component graphite_analyser32=Monochromator_pol() SETTING [/usr/share/mcstas/3.0-dev/optics/Monochromator_pol.comp:109]");
  stracpy(_graphite_analyser32->_name, "graphite_analyser32", 16384);
  stracpy(_graphite_analyser32->_type, "Monochromator_pol", 16384);
  _graphite_analyser32->_index=60;
  _graphite_analyser32->_parameters.zwidth = ANAWIDTH;
  #define zwidth (_graphite_analyser32->_parameters.zwidth)
  _graphite_analyser32->_parameters.yheight = ANAHEIGHT;
  #define yheight (_graphite_analyser32->_parameters.yheight)
  _graphite_analyser32->_parameters.mosaic = ANAMOSAIC;
  #define mosaic (_graphite_analyser32->_parameters.mosaic)
  _graphite_analyser32->_parameters.dspread = 0.0025;
  #define dspread (_graphite_analyser32->_parameters.dspread)
  _graphite_analyser32->_parameters.Q = ANAQ;
  #define Q (_graphite_analyser32->_parameters.Q)
  _graphite_analyser32->_parameters.DM = 0;
  #define DM (_graphite_analyser32->_parameters.DM)
  _graphite_analyser32->_parameters.pThreshold = 0.001;
  #define pThreshold (_graphite_analyser32->_parameters.pThreshold)
  _graphite_analyser32->_parameters.Rup = 1.0;
  #define Rup (_graphite_analyser32->_parameters.Rup)
  _graphite_analyser32->_parameters.Rdown = 1.0;
  #define Rdown (_graphite_analyser32->_parameters.Rdown)
  _graphite_analyser32->_parameters.debug = 0;
  #define debug (_graphite_analyser32->_parameters.debug)

  #define mos_rms (_graphite_analyser32->_parameters.mos_rms)
  #define d_rms (_graphite_analyser32->_parameters.d_rms)
  #define mono_Q (_graphite_analyser32->_parameters.mono_Q)
  #define FN (_graphite_analyser32->_parameters.FN)
  #define FM (_graphite_analyser32->_parameters.FM)

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
  /* component graphite_analyser32=Monochromator_pol() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (90)*DEG2RAD, (-14.015)*DEG2RAD);
    rot_mul(tr1, _armAnalyzer->_rotation_absolute, _graphite_analyser32->_rotation_absolute);
    rot_transpose(_graphite_analyser31->_rotation_absolute, tr1);
    rot_mul(_graphite_analyser32->_rotation_absolute, tr1, _graphite_analyser32->_rotation_relative);
    _graphite_analyser32->_rotation_is_identity =  rot_test_identity(_graphite_analyser32->_rotation_relative);
    tc1 = coords_set(
      0, 0.125, 0.885575);
    rot_transpose(_armAnalyzer->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _graphite_analyser32->_position_absolute = coords_add(_armAnalyzer->_position_absolute, tc2);
    tc1 = coords_sub(_graphite_analyser31->_position_absolute, _graphite_analyser32->_position_absolute);
    _graphite_analyser32->_position_relative = rot_apply(_graphite_analyser32->_rotation_absolute, tc1);
  } /* graphite_analyser32=Monochromator_pol() AT ROTATED */
  DEBUG_COMPONENT("graphite_analyser32", _graphite_analyser32->_position_absolute, _graphite_analyser32->_rotation_absolute);
  instrument->_position_absolute[60] = _graphite_analyser32->_position_absolute;
  instrument->_position_relative[60] = _graphite_analyser32->_position_relative;
  instrument->counter_N[60]  = instrument->counter_P[60] = instrument->counter_P2[60] = 0;
  instrument->counter_AbsorbProp[60]= 0;
  return(0);
} /* _graphite_analyser32_setpos */

/* component graphite_analyser33=Monochromator_pol() SETTING, POSITION/ROTATION */
int _graphite_analyser33_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_graphite_analyser33_setpos] component graphite_analyser33=Monochromator_pol() SETTING [/usr/share/mcstas/3.0-dev/optics/Monochromator_pol.comp:109]");
  stracpy(_graphite_analyser33->_name, "graphite_analyser33", 16384);
  stracpy(_graphite_analyser33->_type, "Monochromator_pol", 16384);
  _graphite_analyser33->_index=61;
  _graphite_analyser33->_parameters.zwidth = ANAWIDTH;
  #define zwidth (_graphite_analyser33->_parameters.zwidth)
  _graphite_analyser33->_parameters.yheight = ANAHEIGHT;
  #define yheight (_graphite_analyser33->_parameters.yheight)
  _graphite_analyser33->_parameters.mosaic = ANAMOSAIC;
  #define mosaic (_graphite_analyser33->_parameters.mosaic)
  _graphite_analyser33->_parameters.dspread = 0.0025;
  #define dspread (_graphite_analyser33->_parameters.dspread)
  _graphite_analyser33->_parameters.Q = ANAQ;
  #define Q (_graphite_analyser33->_parameters.Q)
  _graphite_analyser33->_parameters.DM = 0;
  #define DM (_graphite_analyser33->_parameters.DM)
  _graphite_analyser33->_parameters.pThreshold = 0.001;
  #define pThreshold (_graphite_analyser33->_parameters.pThreshold)
  _graphite_analyser33->_parameters.Rup = 1.0;
  #define Rup (_graphite_analyser33->_parameters.Rup)
  _graphite_analyser33->_parameters.Rdown = 1.0;
  #define Rdown (_graphite_analyser33->_parameters.Rdown)
  _graphite_analyser33->_parameters.debug = 0;
  #define debug (_graphite_analyser33->_parameters.debug)

  #define mos_rms (_graphite_analyser33->_parameters.mos_rms)
  #define d_rms (_graphite_analyser33->_parameters.d_rms)
  #define mono_Q (_graphite_analyser33->_parameters.mono_Q)
  #define FN (_graphite_analyser33->_parameters.FN)
  #define FM (_graphite_analyser33->_parameters.FM)

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
  /* component graphite_analyser33=Monochromator_pol() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (90)*DEG2RAD, (-14.518)*DEG2RAD);
    rot_mul(tr1, _armAnalyzer->_rotation_absolute, _graphite_analyser33->_rotation_absolute);
    rot_transpose(_graphite_analyser32->_rotation_absolute, tr1);
    rot_mul(_graphite_analyser33->_rotation_absolute, tr1, _graphite_analyser33->_rotation_relative);
    _graphite_analyser33->_rotation_is_identity =  rot_test_identity(_graphite_analyser33->_rotation_relative);
    tc1 = coords_set(
      0, 0.135, 0.882991);
    rot_transpose(_armAnalyzer->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _graphite_analyser33->_position_absolute = coords_add(_armAnalyzer->_position_absolute, tc2);
    tc1 = coords_sub(_graphite_analyser32->_position_absolute, _graphite_analyser33->_position_absolute);
    _graphite_analyser33->_position_relative = rot_apply(_graphite_analyser33->_rotation_absolute, tc1);
  } /* graphite_analyser33=Monochromator_pol() AT ROTATED */
  DEBUG_COMPONENT("graphite_analyser33", _graphite_analyser33->_position_absolute, _graphite_analyser33->_rotation_absolute);
  instrument->_position_absolute[61] = _graphite_analyser33->_position_absolute;
  instrument->_position_relative[61] = _graphite_analyser33->_position_relative;
  instrument->counter_N[61]  = instrument->counter_P[61] = instrument->counter_P2[61] = 0;
  instrument->counter_AbsorbProp[61]= 0;
  return(0);
} /* _graphite_analyser33_setpos */

/* component graphite_analyser34=Monochromator_pol() SETTING, POSITION/ROTATION */
int _graphite_analyser34_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_graphite_analyser34_setpos] component graphite_analyser34=Monochromator_pol() SETTING [/usr/share/mcstas/3.0-dev/optics/Monochromator_pol.comp:109]");
  stracpy(_graphite_analyser34->_name, "graphite_analyser34", 16384);
  stracpy(_graphite_analyser34->_type, "Monochromator_pol", 16384);
  _graphite_analyser34->_index=62;
  _graphite_analyser34->_parameters.zwidth = ANAWIDTH;
  #define zwidth (_graphite_analyser34->_parameters.zwidth)
  _graphite_analyser34->_parameters.yheight = ANAHEIGHT;
  #define yheight (_graphite_analyser34->_parameters.yheight)
  _graphite_analyser34->_parameters.mosaic = ANAMOSAIC;
  #define mosaic (_graphite_analyser34->_parameters.mosaic)
  _graphite_analyser34->_parameters.dspread = 0.0025;
  #define dspread (_graphite_analyser34->_parameters.dspread)
  _graphite_analyser34->_parameters.Q = ANAQ;
  #define Q (_graphite_analyser34->_parameters.Q)
  _graphite_analyser34->_parameters.DM = 0;
  #define DM (_graphite_analyser34->_parameters.DM)
  _graphite_analyser34->_parameters.pThreshold = 0.001;
  #define pThreshold (_graphite_analyser34->_parameters.pThreshold)
  _graphite_analyser34->_parameters.Rup = 1.0;
  #define Rup (_graphite_analyser34->_parameters.Rup)
  _graphite_analyser34->_parameters.Rdown = 1.0;
  #define Rdown (_graphite_analyser34->_parameters.Rdown)
  _graphite_analyser34->_parameters.debug = 0;
  #define debug (_graphite_analyser34->_parameters.debug)

  #define mos_rms (_graphite_analyser34->_parameters.mos_rms)
  #define d_rms (_graphite_analyser34->_parameters.d_rms)
  #define mono_Q (_graphite_analyser34->_parameters.mono_Q)
  #define FN (_graphite_analyser34->_parameters.FN)
  #define FM (_graphite_analyser34->_parameters.FM)

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
  /* component graphite_analyser34=Monochromator_pol() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (90)*DEG2RAD, (-15.2215)*DEG2RAD);
    rot_mul(tr1, _armAnalyzer->_rotation_absolute, _graphite_analyser34->_rotation_absolute);
    rot_transpose(_graphite_analyser33->_rotation_absolute, tr1);
    rot_mul(_graphite_analyser34->_rotation_absolute, tr1, _graphite_analyser34->_rotation_relative);
    _graphite_analyser34->_rotation_is_identity =  rot_test_identity(_graphite_analyser34->_rotation_relative);
    tc1 = coords_set(
      0, 0.145, 0.880194);
    rot_transpose(_armAnalyzer->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _graphite_analyser34->_position_absolute = coords_add(_armAnalyzer->_position_absolute, tc2);
    tc1 = coords_sub(_graphite_analyser33->_position_absolute, _graphite_analyser34->_position_absolute);
    _graphite_analyser34->_position_relative = rot_apply(_graphite_analyser34->_rotation_absolute, tc1);
  } /* graphite_analyser34=Monochromator_pol() AT ROTATED */
  DEBUG_COMPONENT("graphite_analyser34", _graphite_analyser34->_position_absolute, _graphite_analyser34->_rotation_absolute);
  instrument->_position_absolute[62] = _graphite_analyser34->_position_absolute;
  instrument->_position_relative[62] = _graphite_analyser34->_position_relative;
  instrument->counter_N[62]  = instrument->counter_P[62] = instrument->counter_P2[62] = 0;
  instrument->counter_AbsorbProp[62]= 0;
  return(0);
} /* _graphite_analyser34_setpos */

/* component graphite_analyser35=Monochromator_pol() SETTING, POSITION/ROTATION */
int _graphite_analyser35_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_graphite_analyser35_setpos] component graphite_analyser35=Monochromator_pol() SETTING [/usr/share/mcstas/3.0-dev/optics/Monochromator_pol.comp:109]");
  stracpy(_graphite_analyser35->_name, "graphite_analyser35", 16384);
  stracpy(_graphite_analyser35->_type, "Monochromator_pol", 16384);
  _graphite_analyser35->_index=63;
  _graphite_analyser35->_parameters.zwidth = ANAWIDTH;
  #define zwidth (_graphite_analyser35->_parameters.zwidth)
  _graphite_analyser35->_parameters.yheight = ANAHEIGHT;
  #define yheight (_graphite_analyser35->_parameters.yheight)
  _graphite_analyser35->_parameters.mosaic = ANAMOSAIC;
  #define mosaic (_graphite_analyser35->_parameters.mosaic)
  _graphite_analyser35->_parameters.dspread = 0.0025;
  #define dspread (_graphite_analyser35->_parameters.dspread)
  _graphite_analyser35->_parameters.Q = ANAQ;
  #define Q (_graphite_analyser35->_parameters.Q)
  _graphite_analyser35->_parameters.DM = 0;
  #define DM (_graphite_analyser35->_parameters.DM)
  _graphite_analyser35->_parameters.pThreshold = 0.001;
  #define pThreshold (_graphite_analyser35->_parameters.pThreshold)
  _graphite_analyser35->_parameters.Rup = 1.0;
  #define Rup (_graphite_analyser35->_parameters.Rup)
  _graphite_analyser35->_parameters.Rdown = 1.0;
  #define Rdown (_graphite_analyser35->_parameters.Rdown)
  _graphite_analyser35->_parameters.debug = 0;
  #define debug (_graphite_analyser35->_parameters.debug)

  #define mos_rms (_graphite_analyser35->_parameters.mos_rms)
  #define d_rms (_graphite_analyser35->_parameters.d_rms)
  #define mono_Q (_graphite_analyser35->_parameters.mono_Q)
  #define FN (_graphite_analyser35->_parameters.FN)
  #define FM (_graphite_analyser35->_parameters.FM)

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
  /* component graphite_analyser35=Monochromator_pol() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (90)*DEG2RAD, (-16.2955)*DEG2RAD);
    rot_mul(tr1, _armAnalyzer->_rotation_absolute, _graphite_analyser35->_rotation_absolute);
    rot_transpose(_graphite_analyser34->_rotation_absolute, tr1);
    rot_mul(_graphite_analyser35->_rotation_absolute, tr1, _graphite_analyser35->_rotation_relative);
    _graphite_analyser35->_rotation_is_identity =  rot_test_identity(_graphite_analyser35->_rotation_relative);
    tc1 = coords_set(
      0, 0.155, 0.877283);
    rot_transpose(_armAnalyzer->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _graphite_analyser35->_position_absolute = coords_add(_armAnalyzer->_position_absolute, tc2);
    tc1 = coords_sub(_graphite_analyser34->_position_absolute, _graphite_analyser35->_position_absolute);
    _graphite_analyser35->_position_relative = rot_apply(_graphite_analyser35->_rotation_absolute, tc1);
  } /* graphite_analyser35=Monochromator_pol() AT ROTATED */
  DEBUG_COMPONENT("graphite_analyser35", _graphite_analyser35->_position_absolute, _graphite_analyser35->_rotation_absolute);
  instrument->_position_absolute[63] = _graphite_analyser35->_position_absolute;
  instrument->_position_relative[63] = _graphite_analyser35->_position_relative;
  instrument->counter_N[63]  = instrument->counter_P[63] = instrument->counter_P2[63] = 0;
  instrument->counter_AbsorbProp[63]= 0;
  return(0);
} /* _graphite_analyser35_setpos */

/* component graphite_analyser36=Monochromator_pol() SETTING, POSITION/ROTATION */
int _graphite_analyser36_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_graphite_analyser36_setpos] component graphite_analyser36=Monochromator_pol() SETTING [/usr/share/mcstas/3.0-dev/optics/Monochromator_pol.comp:109]");
  stracpy(_graphite_analyser36->_name, "graphite_analyser36", 16384);
  stracpy(_graphite_analyser36->_type, "Monochromator_pol", 16384);
  _graphite_analyser36->_index=64;
  _graphite_analyser36->_parameters.zwidth = ANAWIDTH;
  #define zwidth (_graphite_analyser36->_parameters.zwidth)
  _graphite_analyser36->_parameters.yheight = ANAHEIGHT;
  #define yheight (_graphite_analyser36->_parameters.yheight)
  _graphite_analyser36->_parameters.mosaic = ANAMOSAIC;
  #define mosaic (_graphite_analyser36->_parameters.mosaic)
  _graphite_analyser36->_parameters.dspread = 0.0025;
  #define dspread (_graphite_analyser36->_parameters.dspread)
  _graphite_analyser36->_parameters.Q = ANAQ;
  #define Q (_graphite_analyser36->_parameters.Q)
  _graphite_analyser36->_parameters.DM = 0;
  #define DM (_graphite_analyser36->_parameters.DM)
  _graphite_analyser36->_parameters.pThreshold = 0.001;
  #define pThreshold (_graphite_analyser36->_parameters.pThreshold)
  _graphite_analyser36->_parameters.Rup = 1.0;
  #define Rup (_graphite_analyser36->_parameters.Rup)
  _graphite_analyser36->_parameters.Rdown = 1.0;
  #define Rdown (_graphite_analyser36->_parameters.Rdown)
  _graphite_analyser36->_parameters.debug = 0;
  #define debug (_graphite_analyser36->_parameters.debug)

  #define mos_rms (_graphite_analyser36->_parameters.mos_rms)
  #define d_rms (_graphite_analyser36->_parameters.d_rms)
  #define mono_Q (_graphite_analyser36->_parameters.mono_Q)
  #define FN (_graphite_analyser36->_parameters.FN)
  #define FM (_graphite_analyser36->_parameters.FM)

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
  /* component graphite_analyser36=Monochromator_pol() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (90)*DEG2RAD, (-16.666)*DEG2RAD);
    rot_mul(tr1, _armAnalyzer->_rotation_absolute, _graphite_analyser36->_rotation_absolute);
    rot_transpose(_graphite_analyser35->_rotation_absolute, tr1);
    rot_mul(_graphite_analyser36->_rotation_absolute, tr1, _graphite_analyser36->_rotation_relative);
    _graphite_analyser36->_rotation_is_identity =  rot_test_identity(_graphite_analyser36->_rotation_relative);
    tc1 = coords_set(
      0, 0.165, 0.874375);
    rot_transpose(_armAnalyzer->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _graphite_analyser36->_position_absolute = coords_add(_armAnalyzer->_position_absolute, tc2);
    tc1 = coords_sub(_graphite_analyser35->_position_absolute, _graphite_analyser36->_position_absolute);
    _graphite_analyser36->_position_relative = rot_apply(_graphite_analyser36->_rotation_absolute, tc1);
  } /* graphite_analyser36=Monochromator_pol() AT ROTATED */
  DEBUG_COMPONENT("graphite_analyser36", _graphite_analyser36->_position_absolute, _graphite_analyser36->_rotation_absolute);
  instrument->_position_absolute[64] = _graphite_analyser36->_position_absolute;
  instrument->_position_relative[64] = _graphite_analyser36->_position_relative;
  instrument->counter_N[64]  = instrument->counter_P[64] = instrument->counter_P2[64] = 0;
  instrument->counter_AbsorbProp[64]= 0;
  return(0);
} /* _graphite_analyser36_setpos */

/* component graphite_analyser37=Monochromator_pol() SETTING, POSITION/ROTATION */
int _graphite_analyser37_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_graphite_analyser37_setpos] component graphite_analyser37=Monochromator_pol() SETTING [/usr/share/mcstas/3.0-dev/optics/Monochromator_pol.comp:109]");
  stracpy(_graphite_analyser37->_name, "graphite_analyser37", 16384);
  stracpy(_graphite_analyser37->_type, "Monochromator_pol", 16384);
  _graphite_analyser37->_index=65;
  _graphite_analyser37->_parameters.zwidth = ANAWIDTH;
  #define zwidth (_graphite_analyser37->_parameters.zwidth)
  _graphite_analyser37->_parameters.yheight = ANAHEIGHT;
  #define yheight (_graphite_analyser37->_parameters.yheight)
  _graphite_analyser37->_parameters.mosaic = ANAMOSAIC;
  #define mosaic (_graphite_analyser37->_parameters.mosaic)
  _graphite_analyser37->_parameters.dspread = 0.0025;
  #define dspread (_graphite_analyser37->_parameters.dspread)
  _graphite_analyser37->_parameters.Q = ANAQ;
  #define Q (_graphite_analyser37->_parameters.Q)
  _graphite_analyser37->_parameters.DM = 0;
  #define DM (_graphite_analyser37->_parameters.DM)
  _graphite_analyser37->_parameters.pThreshold = 0.001;
  #define pThreshold (_graphite_analyser37->_parameters.pThreshold)
  _graphite_analyser37->_parameters.Rup = 1.0;
  #define Rup (_graphite_analyser37->_parameters.Rup)
  _graphite_analyser37->_parameters.Rdown = 1.0;
  #define Rdown (_graphite_analyser37->_parameters.Rdown)
  _graphite_analyser37->_parameters.debug = 0;
  #define debug (_graphite_analyser37->_parameters.debug)

  #define mos_rms (_graphite_analyser37->_parameters.mos_rms)
  #define d_rms (_graphite_analyser37->_parameters.d_rms)
  #define mono_Q (_graphite_analyser37->_parameters.mono_Q)
  #define FN (_graphite_analyser37->_parameters.FN)
  #define FM (_graphite_analyser37->_parameters.FM)

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
  /* component graphite_analyser37=Monochromator_pol() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (90)*DEG2RAD, (-17.5315)*DEG2RAD);
    rot_mul(tr1, _armAnalyzer->_rotation_absolute, _graphite_analyser37->_rotation_absolute);
    rot_transpose(_graphite_analyser36->_rotation_absolute, tr1);
    rot_mul(_graphite_analyser37->_rotation_absolute, tr1, _graphite_analyser37->_rotation_relative);
    _graphite_analyser37->_rotation_is_identity =  rot_test_identity(_graphite_analyser37->_rotation_relative);
    tc1 = coords_set(
      0, 0.175, 0.871225);
    rot_transpose(_armAnalyzer->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _graphite_analyser37->_position_absolute = coords_add(_armAnalyzer->_position_absolute, tc2);
    tc1 = coords_sub(_graphite_analyser36->_position_absolute, _graphite_analyser37->_position_absolute);
    _graphite_analyser37->_position_relative = rot_apply(_graphite_analyser37->_rotation_absolute, tc1);
  } /* graphite_analyser37=Monochromator_pol() AT ROTATED */
  DEBUG_COMPONENT("graphite_analyser37", _graphite_analyser37->_position_absolute, _graphite_analyser37->_rotation_absolute);
  instrument->_position_absolute[65] = _graphite_analyser37->_position_absolute;
  instrument->_position_relative[65] = _graphite_analyser37->_position_relative;
  instrument->counter_N[65]  = instrument->counter_P[65] = instrument->counter_P2[65] = 0;
  instrument->counter_AbsorbProp[65]= 0;
  return(0);
} /* _graphite_analyser37_setpos */

/* component graphite_analyser38=Monochromator_pol() SETTING, POSITION/ROTATION */
int _graphite_analyser38_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_graphite_analyser38_setpos] component graphite_analyser38=Monochromator_pol() SETTING [/usr/share/mcstas/3.0-dev/optics/Monochromator_pol.comp:109]");
  stracpy(_graphite_analyser38->_name, "graphite_analyser38", 16384);
  stracpy(_graphite_analyser38->_type, "Monochromator_pol", 16384);
  _graphite_analyser38->_index=66;
  _graphite_analyser38->_parameters.zwidth = ANAWIDTH;
  #define zwidth (_graphite_analyser38->_parameters.zwidth)
  _graphite_analyser38->_parameters.yheight = ANAHEIGHT;
  #define yheight (_graphite_analyser38->_parameters.yheight)
  _graphite_analyser38->_parameters.mosaic = ANAMOSAIC;
  #define mosaic (_graphite_analyser38->_parameters.mosaic)
  _graphite_analyser38->_parameters.dspread = 0.0025;
  #define dspread (_graphite_analyser38->_parameters.dspread)
  _graphite_analyser38->_parameters.Q = ANAQ;
  #define Q (_graphite_analyser38->_parameters.Q)
  _graphite_analyser38->_parameters.DM = 0;
  #define DM (_graphite_analyser38->_parameters.DM)
  _graphite_analyser38->_parameters.pThreshold = 0.001;
  #define pThreshold (_graphite_analyser38->_parameters.pThreshold)
  _graphite_analyser38->_parameters.Rup = 1.0;
  #define Rup (_graphite_analyser38->_parameters.Rup)
  _graphite_analyser38->_parameters.Rdown = 1.0;
  #define Rdown (_graphite_analyser38->_parameters.Rdown)
  _graphite_analyser38->_parameters.debug = 0;
  #define debug (_graphite_analyser38->_parameters.debug)

  #define mos_rms (_graphite_analyser38->_parameters.mos_rms)
  #define d_rms (_graphite_analyser38->_parameters.d_rms)
  #define mono_Q (_graphite_analyser38->_parameters.mono_Q)
  #define FN (_graphite_analyser38->_parameters.FN)
  #define FM (_graphite_analyser38->_parameters.FM)

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
  /* component graphite_analyser38=Monochromator_pol() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (90)*DEG2RAD, (-18.604)*DEG2RAD);
    rot_mul(tr1, _armAnalyzer->_rotation_absolute, _graphite_analyser38->_rotation_absolute);
    rot_transpose(_graphite_analyser37->_rotation_absolute, tr1);
    rot_mul(_graphite_analyser38->_rotation_absolute, tr1, _graphite_analyser38->_rotation_relative);
    _graphite_analyser38->_rotation_is_identity =  rot_test_identity(_graphite_analyser38->_rotation_relative);
    tc1 = coords_set(
      0, 0.185, 0.867897);
    rot_transpose(_armAnalyzer->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _graphite_analyser38->_position_absolute = coords_add(_armAnalyzer->_position_absolute, tc2);
    tc1 = coords_sub(_graphite_analyser37->_position_absolute, _graphite_analyser38->_position_absolute);
    _graphite_analyser38->_position_relative = rot_apply(_graphite_analyser38->_rotation_absolute, tc1);
  } /* graphite_analyser38=Monochromator_pol() AT ROTATED */
  DEBUG_COMPONENT("graphite_analyser38", _graphite_analyser38->_position_absolute, _graphite_analyser38->_rotation_absolute);
  instrument->_position_absolute[66] = _graphite_analyser38->_position_absolute;
  instrument->_position_relative[66] = _graphite_analyser38->_position_relative;
  instrument->counter_N[66]  = instrument->counter_P[66] = instrument->counter_P2[66] = 0;
  instrument->counter_AbsorbProp[66]= 0;
  return(0);
} /* _graphite_analyser38_setpos */

/* component graphite_analyser39=Monochromator_pol() SETTING, POSITION/ROTATION */
int _graphite_analyser39_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_graphite_analyser39_setpos] component graphite_analyser39=Monochromator_pol() SETTING [/usr/share/mcstas/3.0-dev/optics/Monochromator_pol.comp:109]");
  stracpy(_graphite_analyser39->_name, "graphite_analyser39", 16384);
  stracpy(_graphite_analyser39->_type, "Monochromator_pol", 16384);
  _graphite_analyser39->_index=67;
  _graphite_analyser39->_parameters.zwidth = ANAWIDTH;
  #define zwidth (_graphite_analyser39->_parameters.zwidth)
  _graphite_analyser39->_parameters.yheight = ANAHEIGHT;
  #define yheight (_graphite_analyser39->_parameters.yheight)
  _graphite_analyser39->_parameters.mosaic = ANAMOSAIC;
  #define mosaic (_graphite_analyser39->_parameters.mosaic)
  _graphite_analyser39->_parameters.dspread = 0.0025;
  #define dspread (_graphite_analyser39->_parameters.dspread)
  _graphite_analyser39->_parameters.Q = ANAQ;
  #define Q (_graphite_analyser39->_parameters.Q)
  _graphite_analyser39->_parameters.DM = 0;
  #define DM (_graphite_analyser39->_parameters.DM)
  _graphite_analyser39->_parameters.pThreshold = 0.001;
  #define pThreshold (_graphite_analyser39->_parameters.pThreshold)
  _graphite_analyser39->_parameters.Rup = 1.0;
  #define Rup (_graphite_analyser39->_parameters.Rup)
  _graphite_analyser39->_parameters.Rdown = 1.0;
  #define Rdown (_graphite_analyser39->_parameters.Rdown)
  _graphite_analyser39->_parameters.debug = 0;
  #define debug (_graphite_analyser39->_parameters.debug)

  #define mos_rms (_graphite_analyser39->_parameters.mos_rms)
  #define d_rms (_graphite_analyser39->_parameters.d_rms)
  #define mono_Q (_graphite_analyser39->_parameters.mono_Q)
  #define FN (_graphite_analyser39->_parameters.FN)
  #define FM (_graphite_analyser39->_parameters.FM)

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
  /* component graphite_analyser39=Monochromator_pol() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (90)*DEG2RAD, (-18.811)*DEG2RAD);
    rot_mul(tr1, _armAnalyzer->_rotation_absolute, _graphite_analyser39->_rotation_absolute);
    rot_transpose(_graphite_analyser38->_rotation_absolute, tr1);
    rot_mul(_graphite_analyser39->_rotation_absolute, tr1, _graphite_analyser39->_rotation_relative);
    _graphite_analyser39->_rotation_is_identity =  rot_test_identity(_graphite_analyser39->_rotation_relative);
    tc1 = coords_set(
      0, 0.195, 0.864474);
    rot_transpose(_armAnalyzer->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _graphite_analyser39->_position_absolute = coords_add(_armAnalyzer->_position_absolute, tc2);
    tc1 = coords_sub(_graphite_analyser38->_position_absolute, _graphite_analyser39->_position_absolute);
    _graphite_analyser39->_position_relative = rot_apply(_graphite_analyser39->_rotation_absolute, tc1);
  } /* graphite_analyser39=Monochromator_pol() AT ROTATED */
  DEBUG_COMPONENT("graphite_analyser39", _graphite_analyser39->_position_absolute, _graphite_analyser39->_rotation_absolute);
  instrument->_position_absolute[67] = _graphite_analyser39->_position_absolute;
  instrument->_position_relative[67] = _graphite_analyser39->_position_relative;
  instrument->counter_N[67]  = instrument->counter_P[67] = instrument->counter_P2[67] = 0;
  instrument->counter_AbsorbProp[67]= 0;
  return(0);
} /* _graphite_analyser39_setpos */

/* component armHelp=Arm() SETTING, POSITION/ROTATION */
int _armHelp_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_armHelp_setpos] component armHelp=Arm() SETTING [Arm:0]");
  stracpy(_armHelp->_name, "armHelp", 16384);
  stracpy(_armHelp->_type, "Arm", 16384);
  _armHelp->_index=68;
  /* component armHelp=Arm() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (180)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _armAnalyzer->_rotation_absolute, _armHelp->_rotation_absolute);
    rot_transpose(_graphite_analyser39->_rotation_absolute, tr1);
    rot_mul(_armHelp->_rotation_absolute, tr1, _armHelp->_rotation_relative);
    _armHelp->_rotation_is_identity =  rot_test_identity(_armHelp->_rotation_relative);
    tc1 = coords_set(
      0, 0, ANARADIUS);
    rot_transpose(_armAnalyzer->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _armHelp->_position_absolute = coords_add(_armAnalyzer->_position_absolute, tc2);
    tc1 = coords_sub(_graphite_analyser39->_position_absolute, _armHelp->_position_absolute);
    _armHelp->_position_relative = rot_apply(_armHelp->_rotation_absolute, tc1);
  } /* armHelp=Arm() AT ROTATED */
  DEBUG_COMPONENT("armHelp", _armHelp->_position_absolute, _armHelp->_rotation_absolute);
  instrument->_position_absolute[68] = _armHelp->_position_absolute;
  instrument->_position_relative[68] = _armHelp->_position_relative;
  instrument->counter_N[68]  = instrument->counter_P[68] = instrument->counter_P2[68] = 0;
  instrument->counter_AbsorbProp[68]= 0;
  return(0);
} /* _armHelp_setpos */

/* component armDetector=Arm() SETTING, POSITION/ROTATION */
int _armDetector_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_armDetector_setpos] component armDetector=Arm() SETTING [Arm:0]");
  stracpy(_armDetector->_name, "armDetector", 16384);
  stracpy(_armDetector->_type, "Arm", 16384);
  _armDetector->_index=69;
  /* component armDetector=Arm() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (2.0 * 4.6526)*DEG2RAD, (0)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _armHelp->_rotation_absolute, _armDetector->_rotation_absolute);
    rot_transpose(_armHelp->_rotation_absolute, tr1);
    rot_mul(_armDetector->_rotation_absolute, tr1, _armDetector->_rotation_relative);
    _armDetector->_rotation_is_identity =  rot_test_identity(_armDetector->_rotation_relative);
    tc1 = coords_set(
      0, 0, 0);
    rot_transpose(_armHelp->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _armDetector->_position_absolute = coords_add(_armHelp->_position_absolute, tc2);
    tc1 = coords_sub(_armHelp->_position_absolute, _armDetector->_position_absolute);
    _armDetector->_position_relative = rot_apply(_armDetector->_rotation_absolute, tc1);
  } /* armDetector=Arm() AT ROTATED */
  DEBUG_COMPONENT("armDetector", _armDetector->_position_absolute, _armDetector->_rotation_absolute);
  instrument->_position_absolute[69] = _armDetector->_position_absolute;
  instrument->_position_relative[69] = _armDetector->_position_relative;
  instrument->counter_N[69]  = instrument->counter_P[69] = instrument->counter_P2[69] = 0;
  instrument->counter_AbsorbProp[69]= 0;
  return(0);
} /* _armDetector_setpos */

/* component lamDetector=L_monitor() SETTING, POSITION/ROTATION */
int _lamDetector_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_lamDetector_setpos] component lamDetector=L_monitor() SETTING [/usr/share/mcstas/3.0-dev/monitors/L_monitor.comp:65]");
  stracpy(_lamDetector->_name, "lamDetector", 16384);
  stracpy(_lamDetector->_type, "L_monitor", 16384);
  _lamDetector->_index=70;
  _lamDetector->_parameters.nL = 120;
  #define nL (_lamDetector->_parameters.nL)
  if("lambdaDetector.dat" && strlen("lambdaDetector.dat"))
    stracpy(_lamDetector->_parameters.filename, "lambdaDetector.dat" ? "lambdaDetector.dat" : "", 16384);
  else 
  _lamDetector->_parameters.filename[0]='\0';
  #define filename (_lamDetector->_parameters.filename)
  _lamDetector->_parameters.xmin = -0.05;
  #define xmin (_lamDetector->_parameters.xmin)
  _lamDetector->_parameters.xmax = 0.05;
  #define xmax (_lamDetector->_parameters.xmax)
  _lamDetector->_parameters.ymin = -0.05;
  #define ymin (_lamDetector->_parameters.ymin)
  _lamDetector->_parameters.ymax = 0.05;
  #define ymax (_lamDetector->_parameters.ymax)
  _lamDetector->_parameters.xwidth = 0.10;
  #define xwidth (_lamDetector->_parameters.xwidth)
  _lamDetector->_parameters.yheight = 0.10;
  #define yheight (_lamDetector->_parameters.yheight)
  _lamDetector->_parameters.Lmin = 0.99 * instrument->_parameters._LAMBDA;
  #define Lmin (_lamDetector->_parameters.Lmin)
  _lamDetector->_parameters.Lmax = 1.01 * instrument->_parameters._LAMBDA;
  #define Lmax (_lamDetector->_parameters.Lmax)
  _lamDetector->_parameters.restore_neutron = 0;
  #define restore_neutron (_lamDetector->_parameters.restore_neutron)

  #define L_N (_lamDetector->_parameters.L_N)
  #define L_p (_lamDetector->_parameters.L_p)
  #define L_p2 (_lamDetector->_parameters.L_p2)

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
  /* component lamDetector=L_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _armDetector->_rotation_absolute, _lamDetector->_rotation_absolute);
    rot_transpose(_armDetector->_rotation_absolute, tr1);
    rot_mul(_lamDetector->_rotation_absolute, tr1, _lamDetector->_rotation_relative);
    _lamDetector->_rotation_is_identity =  rot_test_identity(_lamDetector->_rotation_relative);
    tc1 = coords_set(
      0, 0, DETECTORDISTANCE);
    rot_transpose(_armDetector->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _lamDetector->_position_absolute = coords_add(_armDetector->_position_absolute, tc2);
    tc1 = coords_sub(_armDetector->_position_absolute, _lamDetector->_position_absolute);
    _lamDetector->_position_relative = rot_apply(_lamDetector->_rotation_absolute, tc1);
  } /* lamDetector=L_monitor() AT ROTATED */
  DEBUG_COMPONENT("lamDetector", _lamDetector->_position_absolute, _lamDetector->_rotation_absolute);
  instrument->_position_absolute[70] = _lamDetector->_position_absolute;
  instrument->_position_relative[70] = _lamDetector->_position_relative;
  instrument->counter_N[70]  = instrument->counter_P[70] = instrument->counter_P2[70] = 0;
  instrument->counter_AbsorbProp[70]= 0;
  return(0);
} /* _lamDetector_setpos */

/* component psdDetector=PSD_monitor() SETTING, POSITION/ROTATION */
int _psdDetector_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_psdDetector_setpos] component psdDetector=PSD_monitor() SETTING [/usr/share/mcstas/3.0-dev/monitors/PSD_monitor.comp:64]");
  stracpy(_psdDetector->_name, "psdDetector", 16384);
  stracpy(_psdDetector->_type, "PSD_monitor", 16384);
  _psdDetector->_index=71;
  _psdDetector->_parameters.nx = 50;
  #define nx (_psdDetector->_parameters.nx)
  _psdDetector->_parameters.ny = 50;
  #define ny (_psdDetector->_parameters.ny)
  if("psdDetector.dat" && strlen("psdDetector.dat"))
    stracpy(_psdDetector->_parameters.filename, "psdDetector.dat" ? "psdDetector.dat" : "", 16384);
  else 
  _psdDetector->_parameters.filename[0]='\0';
  #define filename (_psdDetector->_parameters.filename)
  _psdDetector->_parameters.xmin = -0.05;
  #define xmin (_psdDetector->_parameters.xmin)
  _psdDetector->_parameters.xmax = 0.05;
  #define xmax (_psdDetector->_parameters.xmax)
  _psdDetector->_parameters.ymin = -0.05;
  #define ymin (_psdDetector->_parameters.ymin)
  _psdDetector->_parameters.ymax = 0.05;
  #define ymax (_psdDetector->_parameters.ymax)
  _psdDetector->_parameters.xwidth = 0.10;
  #define xwidth (_psdDetector->_parameters.xwidth)
  _psdDetector->_parameters.yheight = 0.10;
  #define yheight (_psdDetector->_parameters.yheight)
  _psdDetector->_parameters.restore_neutron = 0;
  #define restore_neutron (_psdDetector->_parameters.restore_neutron)

  #define PSD_N (_psdDetector->_parameters.PSD_N)
  #define PSD_p (_psdDetector->_parameters.PSD_p)
  #define PSD_p2 (_psdDetector->_parameters.PSD_p2)

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
  /* component psdDetector=PSD_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _armDetector->_rotation_absolute, _psdDetector->_rotation_absolute);
    rot_transpose(_lamDetector->_rotation_absolute, tr1);
    rot_mul(_psdDetector->_rotation_absolute, tr1, _psdDetector->_rotation_relative);
    _psdDetector->_rotation_is_identity =  rot_test_identity(_psdDetector->_rotation_relative);
    tc1 = coords_set(
      0, 0, DETECTORDISTANCE);
    rot_transpose(_armDetector->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _psdDetector->_position_absolute = coords_add(_armDetector->_position_absolute, tc2);
    tc1 = coords_sub(_lamDetector->_position_absolute, _psdDetector->_position_absolute);
    _psdDetector->_position_relative = rot_apply(_psdDetector->_rotation_absolute, tc1);
  } /* psdDetector=PSD_monitor() AT ROTATED */
  DEBUG_COMPONENT("psdDetector", _psdDetector->_position_absolute, _psdDetector->_rotation_absolute);
  instrument->_position_absolute[71] = _psdDetector->_position_absolute;
  instrument->_position_relative[71] = _psdDetector->_position_relative;
  instrument->counter_N[71]  = instrument->counter_P[71] = instrument->counter_P2[71] = 0;
  instrument->counter_AbsorbProp[71]= 0;
  return(0);
} /* _psdDetector_setpos */

/* component tofDetector=TOF_monitor() SETTING, POSITION/ROTATION */
int _tofDetector_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_tofDetector_setpos] component tofDetector=TOF_monitor() SETTING [/usr/share/mcstas/3.0-dev/monitors/TOF_monitor.comp:63]");
  stracpy(_tofDetector->_name, "tofDetector", 16384);
  stracpy(_tofDetector->_type, "TOF_monitor", 16384);
  _tofDetector->_index=72;
  _tofDetector->_parameters.nt = 120;
  #define nt (_tofDetector->_parameters.nt)
  if("tofDetector.dat" && strlen("tofDetector.dat"))
    stracpy(_tofDetector->_parameters.filename, "tofDetector.dat" ? "tofDetector.dat" : "", 16384);
  else 
  _tofDetector->_parameters.filename[0]='\0';
  #define filename (_tofDetector->_parameters.filename)
  _tofDetector->_parameters.xmin = -0.05;
  #define xmin (_tofDetector->_parameters.xmin)
  _tofDetector->_parameters.xmax = 0.05;
  #define xmax (_tofDetector->_parameters.xmax)
  _tofDetector->_parameters.ymin = -0.05;
  #define ymin (_tofDetector->_parameters.ymin)
  _tofDetector->_parameters.ymax = 0.05;
  #define ymax (_tofDetector->_parameters.ymax)
  _tofDetector->_parameters.xwidth = 0.10;
  #define xwidth (_tofDetector->_parameters.xwidth)
  _tofDetector->_parameters.yheight = 0.10;
  #define yheight (_tofDetector->_parameters.yheight)
  _tofDetector->_parameters.tmin = calcT0 ( instrument->_parameters._LAMBDA , DMODTARGET + ANARADIUS + DETECTORDISTANCE -1.5 ) * 1e6;
  #define tmin (_tofDetector->_parameters.tmin)
  _tofDetector->_parameters.tmax = calcT0 ( instrument->_parameters._LAMBDA , DMODTARGET + ANARADIUS + DETECTORDISTANCE + 1.5 ) * 1e6;
  #define tmax (_tofDetector->_parameters.tmax)
  _tofDetector->_parameters.dt = 1.0;
  #define dt (_tofDetector->_parameters.dt)
  _tofDetector->_parameters.restore_neutron = 0;
  #define restore_neutron (_tofDetector->_parameters.restore_neutron)

  #define TOF_N (_tofDetector->_parameters.TOF_N)
  #define TOF_p (_tofDetector->_parameters.TOF_p)
  #define TOF_p2 (_tofDetector->_parameters.TOF_p2)
  #define t_min (_tofDetector->_parameters.t_min)
  #define t_max (_tofDetector->_parameters.t_max)
  #define delta_t (_tofDetector->_parameters.delta_t)

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
  /* component tofDetector=TOF_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _armDetector->_rotation_absolute, _tofDetector->_rotation_absolute);
    rot_transpose(_psdDetector->_rotation_absolute, tr1);
    rot_mul(_tofDetector->_rotation_absolute, tr1, _tofDetector->_rotation_relative);
    _tofDetector->_rotation_is_identity =  rot_test_identity(_tofDetector->_rotation_relative);
    tc1 = coords_set(
      0, 0, DETECTORDISTANCE);
    rot_transpose(_armDetector->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _tofDetector->_position_absolute = coords_add(_armDetector->_position_absolute, tc2);
    tc1 = coords_sub(_psdDetector->_position_absolute, _tofDetector->_position_absolute);
    _tofDetector->_position_relative = rot_apply(_tofDetector->_rotation_absolute, tc1);
  } /* tofDetector=TOF_monitor() AT ROTATED */
  DEBUG_COMPONENT("tofDetector", _tofDetector->_position_absolute, _tofDetector->_rotation_absolute);
  instrument->_position_absolute[72] = _tofDetector->_position_absolute;
  instrument->_position_relative[72] = _tofDetector->_position_relative;
  instrument->counter_N[72]  = instrument->counter_P[72] = instrument->counter_P2[72] = 0;
  instrument->counter_AbsorbProp[72]= 0;
  return(0);
} /* _tofDetector_setpos */

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

_class_Guide_curved *class_Guide_curved_init(_class_Guide_curved *_comp
) {
  #define w1 (_comp->_parameters.w1)
  #define h1 (_comp->_parameters.h1)
  #define l (_comp->_parameters.l)
  #define R0 (_comp->_parameters.R0)
  #define Qc (_comp->_parameters.Qc)
  #define alpha (_comp->_parameters.alpha)
  #define m (_comp->_parameters.m)
  #define W (_comp->_parameters.W)
  #define curvature (_comp->_parameters.curvature)
if (mcgravitation) fprintf(stderr,"WARNING: Guide_curved: %s: "
    "This component produces wrong results with gravitation !\n"
    "Use Guide_gravity.\n",
    NAME_CURRENT_COMP);
  #undef w1
  #undef h1
  #undef l
  #undef R0
  #undef Qc
  #undef alpha
  #undef m
  #undef W
  #undef curvature
  return(_comp);
} /* class_Guide_curved_init */

_class_Incoherent *class_Incoherent_init(_class_Incoherent *_comp
) {
  #define geometry (_comp->_parameters.geometry)
  #define radius (_comp->_parameters.radius)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define zdepth (_comp->_parameters.zdepth)
  #define thickness (_comp->_parameters.thickness)
  #define target_x (_comp->_parameters.target_x)
  #define target_y (_comp->_parameters.target_y)
  #define target_z (_comp->_parameters.target_z)
  #define focus_r (_comp->_parameters.focus_r)
  #define focus_xw (_comp->_parameters.focus_xw)
  #define focus_yh (_comp->_parameters.focus_yh)
  #define focus_aw (_comp->_parameters.focus_aw)
  #define focus_ah (_comp->_parameters.focus_ah)
  #define target_index (_comp->_parameters.target_index)
  #define pack (_comp->_parameters.pack)
  #define p_interact (_comp->_parameters.p_interact)
  #define f_QE (_comp->_parameters.f_QE)
  #define gamma (_comp->_parameters.gamma)
  #define sigma_abs (_comp->_parameters.sigma_abs)
  #define sigma_inc (_comp->_parameters.sigma_inc)
  #define Vc (_comp->_parameters.Vc)
  #define concentric (_comp->_parameters.concentric)
  #define order (_comp->_parameters.order)
  #define VarsInc (_comp->_parameters.VarsInc)
  #define offdata (_comp->_parameters.offdata)
  VarsInc.shape=-1; /* -1:no shape, 0:cyl, 1:box, 2:sphere, 3:any-shape  */
  if (geometry && strlen(geometry) && strcmp(geometry, "NULL") && strcmp(geometry, "0")) {
	  if (off_init(geometry, xwidth, yheight, zdepth, 0, &offdata)) {
      VarsInc.shape=3; thickness=0; concentric=0;
    }
  }
  else if (xwidth && yheight && zdepth)  VarsInc.shape=1; /* box */
  else if (radius > 0 &&  yheight)       VarsInc.shape=0; /* cylinder */
  else if (radius > 0 && !yheight)       VarsInc.shape=2; /* sphere */

  if (VarsInc.shape < 0)
    exit(fprintf(stderr,"Incoherent: %s: sample has invalid dimensions.\n"
                        "ERROR       Please check parameter values (xwidth, yheight, zdepth, radius).\n", NAME_CURRENT_COMP));
  if (thickness) {
    if (radius && (radius < thickness || ( yheight && (yheight < 2*thickness)))) {
      fprintf(stderr,"Incoherent: %s: hollow sample thickness is larger than its volume (sphere/cylinder).\n"
                     "WARNING     Please check parameter values. Using bulk sample (thickness=0).\n", NAME_CURRENT_COMP);
      thickness=0;
    }
    else if (!radius && (xwidth < 2*thickness || yheight < 2*thickness || zdepth < 2*thickness)) {
      fprintf(stderr,"Incoherent: %s: hollow sample thickness is larger than its volume (box).\n"
                     "WARNING     Please check parameter values. Using bulk sample (thickness=0).\n", NAME_CURRENT_COMP);
      thickness=0;
    }
  }

  if (concentric && thickness<=0) {
    printf("Incoherent: %s:Can not use concentric mode\n"
           "WARNING     on non hollow shape. Ignoring.\n",
           NAME_CURRENT_COMP);
    concentric=0;
  }

  VarsInc.sigma_a= sigma_abs;
  VarsInc.sigma_i= sigma_inc;
  VarsInc.rho    = (pack/Vc);
  VarsInc.my_s   = (VarsInc.rho * 100 * VarsInc.sigma_i);
  VarsInc.my_a_v = (VarsInc.rho * 100 * VarsInc.sigma_a);

  /* now compute target coords if a component index is supplied */
  VarsInc.tx= VarsInc.ty=VarsInc.tz=0;
  if (!target_index && !target_x && !target_y && !target_z) target_index=1;
  if (target_index)
  {
    Coords ToTarget;
    ToTarget = coords_sub(POS_A_COMP_INDEX(INDEX_CURRENT_COMP+target_index),POS_A_CURRENT_COMP);
    ToTarget = rot_apply(ROT_A_CURRENT_COMP, ToTarget);
    coords_get(ToTarget, &VarsInc.tx, &VarsInc.ty, &VarsInc.tz);
  }
  else
  { VarsInc.tx = target_x; VarsInc.ty = target_y; VarsInc.tz = target_z; }

  if (!(VarsInc.tx || VarsInc.ty || VarsInc.tz)) {
    MPI_MASTER(
    printf("Incoherent: %s: The target is not defined. Using direct beam (Z-axis).\n",
      NAME_CURRENT_COMP);
    );
    VarsInc.tz=1;
  }

  /* different ways of setting rectangular area */
  VarsInc.aw  = VarsInc.ah = 0;
  if (focus_xw) { VarsInc.xw = focus_xw; }
  if (focus_yh) { VarsInc.yh = focus_yh; }
  if (focus_aw) { VarsInc.aw = DEG2RAD*focus_aw; }
  if (focus_ah) { VarsInc.ah = DEG2RAD*focus_ah; }

  MPI_MASTER(
  printf("Incoherent: %s: Vc=%g [Angs] sigma_abs=%g [barn] sigma_inc=%g [barn]\n",
      NAME_CURRENT_COMP, Vc, VarsInc.sigma_a, VarsInc.sigma_i);
  );

  #undef geometry
  #undef radius
  #undef xwidth
  #undef yheight
  #undef zdepth
  #undef thickness
  #undef target_x
  #undef target_y
  #undef target_z
  #undef focus_r
  #undef focus_xw
  #undef focus_yh
  #undef focus_aw
  #undef focus_ah
  #undef target_index
  #undef pack
  #undef p_interact
  #undef f_QE
  #undef gamma
  #undef sigma_abs
  #undef sigma_inc
  #undef Vc
  #undef concentric
  #undef order
  #undef VarsInc
  #undef offdata
  return(_comp);
} /* class_Incoherent_init */

_class_Beamstop *class_Beamstop_init(_class_Beamstop *_comp
) {
  #define xmin (_comp->_parameters.xmin)
  #define xmax (_comp->_parameters.xmax)
  #define ymin (_comp->_parameters.ymin)
  #define ymax (_comp->_parameters.ymax)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define radius (_comp->_parameters.radius)
if (xwidth  > 0) { xmax = xwidth/2;  xmin = -xmax; }
  if (yheight > 0) { ymax = yheight/2; ymin = -ymax; }

  if (xmin == 0 && xmax == 0 && ymin == 0 & ymax == 0 && radius == 0)
  { fprintf(stderr,"Beamstop: %s: Error: give geometry\n", NAME_CURRENT_COMP); exit(-1); }
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef xwidth
  #undef yheight
  #undef radius
  return(_comp);
} /* class_Beamstop_init */

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



int init(void) { /* called by mccode_main for ISIS_OSIRIS:INITIALISE */
  DEBUG_INSTR();

  /* code_main/parseoptions/readparams sets instrument parameters value */
  stracpy(instrument->_name, "ISIS_OSIRIS", 256);

  /* Instrument 'ISIS_OSIRIS' INITIALISE */
  SIG_MESSAGE("[ISIS_OSIRIS] INITIALISE [ISIS_OSIRIS.instr:128]");
  #define LAMBDA (instrument->_parameters._LAMBDA)
  #define DLAMBDA (instrument->_parameters._DLAMBDA)
  #define GUIDEREFLECTIVITY (instrument->_parameters._GUIDEREFLECTIVITY)
{

  DMODTARGET = DMODGUIDE+DSTRAIGHTGUIDE1+DCHOPPER1+DCURVEDGUIDE1+
    DCHOPPER2+DCURVEDGUIDE2+DSTRAIGHTGUIDE2+DSUPERGUIDE+DGUIDETARGET; // m

  if(LAMBDA > 0) {

    // make a narrow interval around the selected LAMBDA

    LAMBDAMIN =  LAMBDA - DLAMBDA;
    LAMBDAMAX =  LAMBDA + DLAMBDA;
    if (LAMBDAMIN < 1.0)
      LAMBDAMIN = 1.0;

  }

  printf("\nParameters:\nLambda: %f AA", LAMBDA);
  printf("\nLambdamin, max: %f AA, %f AA", LAMBDAMIN, LAMBDAMAX);
  printf("\nDistance from moderator to target: %f", DMODTARGET);
  printf("\nChopper1 - omega: %f rad/s, window: %f deg",
	 CHOPPER1OMEGA, CHOPPER1WINDOW);
  printf("\nChopper2 - omega: %f rad/s, window: %f deg\n\n",
	 CHOPPER2OMEGA, CHOPPER2WINDOW);

}
  #undef LAMBDA
  #undef DLAMBDA
  #undef GUIDEREFLECTIVITY
  _armStart_setpos(); /* type Arm */
  _isis_mod_setpos(); /* type ISIS_moderator */
  _tofSource_setpos(); /* type TOF_monitor */
  _armStraight1_setpos(); /* type Arm */
  _lamStart_setpos(); /* type L_monitor */
  _psdStart_setpos(); /* type PSD_monitor */
  _guide_straight_1_setpos(); /* type Guide */
  _armChopper1_setpos(); /* type Arm */
  _chopper1_setpos(); /* type DiskChopper */
  _armCurved1_setpos(); /* type Arm */
  _guide_curved1_setpos(); /* type Guide_curved */
  _armChopper2_setpos(); /* type Arm */
  _chopper2_setpos(); /* type DiskChopper */
  _armCurved2_setpos(); /* type Arm */
  _guide_curved2_setpos(); /* type Guide_curved */
  _armStraight2_setpos(); /* type Arm */
  _guide_straight_2_setpos(); /* type Guide */
  _armSuper_setpos(); /* type Arm */
  _guide_super_setpos(); /* type Guide */
  _armTarget_setpos(); /* type Arm */
  _lamTarget_setpos(); /* type L_monitor */
  _psdTarget_setpos(); /* type PSD_monitor */
  _vsample1_setpos(); /* type Incoherent */
  _vsample2_setpos(); /* type Incoherent */
  _vsample3_setpos(); /* type Incoherent */
  _beamStop_setpos(); /* type Beamstop */
  _armAnalyzer_setpos(); /* type Arm */
  _graphite_analyser0_setpos(); /* type Monochromator_pol */
  _graphite_analyser1_setpos(); /* type Monochromator_pol */
  _graphite_analyser2_setpos(); /* type Monochromator_pol */
  _graphite_analyser3_setpos(); /* type Monochromator_pol */
  _graphite_analyser4_setpos(); /* type Monochromator_pol */
  _graphite_analyser5_setpos(); /* type Monochromator_pol */
  _graphite_analyser6_setpos(); /* type Monochromator_pol */
  _graphite_analyser7_setpos(); /* type Monochromator_pol */
  _graphite_analyser8_setpos(); /* type Monochromator_pol */
  _graphite_analyser9_setpos(); /* type Monochromator_pol */
  _graphite_analyser10_setpos(); /* type Monochromator_pol */
  _graphite_analyser11_setpos(); /* type Monochromator_pol */
  _graphite_analyser12_setpos(); /* type Monochromator_pol */
  _graphite_analyser13_setpos(); /* type Monochromator_pol */
  _graphite_analyser14_setpos(); /* type Monochromator_pol */
  _graphite_analyser15_setpos(); /* type Monochromator_pol */
  _graphite_analyser16_setpos(); /* type Monochromator_pol */
  _graphite_analyser17_setpos(); /* type Monochromator_pol */
  _graphite_analyser18_setpos(); /* type Monochromator_pol */
  _graphite_analyser19_setpos(); /* type Monochromator_pol */
  _graphite_analyser20_setpos(); /* type Monochromator_pol */
  _graphite_analyser21_setpos(); /* type Monochromator_pol */
  _graphite_analyser22_setpos(); /* type Monochromator_pol */
  _graphite_analyser23_setpos(); /* type Monochromator_pol */
  _graphite_analyser24_setpos(); /* type Monochromator_pol */
  _graphite_analyser25_setpos(); /* type Monochromator_pol */
  _graphite_analyser26_setpos(); /* type Monochromator_pol */
  _graphite_analyser27_setpos(); /* type Monochromator_pol */
  _graphite_analyser28_setpos(); /* type Monochromator_pol */
  _graphite_analyser29_setpos(); /* type Monochromator_pol */
  _graphite_analyser30_setpos(); /* type Monochromator_pol */
  _graphite_analyser31_setpos(); /* type Monochromator_pol */
  _graphite_analyser32_setpos(); /* type Monochromator_pol */
  _graphite_analyser33_setpos(); /* type Monochromator_pol */
  _graphite_analyser34_setpos(); /* type Monochromator_pol */
  _graphite_analyser35_setpos(); /* type Monochromator_pol */
  _graphite_analyser36_setpos(); /* type Monochromator_pol */
  _graphite_analyser37_setpos(); /* type Monochromator_pol */
  _graphite_analyser38_setpos(); /* type Monochromator_pol */
  _graphite_analyser39_setpos(); /* type Monochromator_pol */
  _armHelp_setpos(); /* type Arm */
  _armDetector_setpos(); /* type Arm */
  _lamDetector_setpos(); /* type L_monitor */
  _psdDetector_setpos(); /* type PSD_monitor */
  _tofDetector_setpos(); /* type TOF_monitor */

  /* call iteratively all components INITIALISE */

  class_ISIS_moderator_init(_isis_mod);

  class_TOF_monitor_init(_tofSource);


  class_L_monitor_init(_lamStart);

  class_PSD_monitor_init(_psdStart);

  class_Guide_init(_guide_straight_1);


  class_DiskChopper_init(_chopper1);


  class_Guide_curved_init(_guide_curved1);


  class_DiskChopper_init(_chopper2);


  class_Guide_curved_init(_guide_curved2);


  class_Guide_init(_guide_straight_2);


  class_Guide_init(_guide_super);


  class_L_monitor_init(_lamTarget);

  class_PSD_monitor_init(_psdTarget);

  class_Incoherent_init(_vsample1);

  class_Incoherent_init(_vsample2);

  class_Incoherent_init(_vsample3);

  class_Beamstop_init(_beamStop);


  class_Monochromator_pol_init(_graphite_analyser0);

  class_Monochromator_pol_init(_graphite_analyser1);

  class_Monochromator_pol_init(_graphite_analyser2);

  class_Monochromator_pol_init(_graphite_analyser3);

  class_Monochromator_pol_init(_graphite_analyser4);

  class_Monochromator_pol_init(_graphite_analyser5);

  class_Monochromator_pol_init(_graphite_analyser6);

  class_Monochromator_pol_init(_graphite_analyser7);

  class_Monochromator_pol_init(_graphite_analyser8);

  class_Monochromator_pol_init(_graphite_analyser9);

  class_Monochromator_pol_init(_graphite_analyser10);

  class_Monochromator_pol_init(_graphite_analyser11);

  class_Monochromator_pol_init(_graphite_analyser12);

  class_Monochromator_pol_init(_graphite_analyser13);

  class_Monochromator_pol_init(_graphite_analyser14);

  class_Monochromator_pol_init(_graphite_analyser15);

  class_Monochromator_pol_init(_graphite_analyser16);

  class_Monochromator_pol_init(_graphite_analyser17);

  class_Monochromator_pol_init(_graphite_analyser18);

  class_Monochromator_pol_init(_graphite_analyser19);

  class_Monochromator_pol_init(_graphite_analyser20);

  class_Monochromator_pol_init(_graphite_analyser21);

  class_Monochromator_pol_init(_graphite_analyser22);

  class_Monochromator_pol_init(_graphite_analyser23);

  class_Monochromator_pol_init(_graphite_analyser24);

  class_Monochromator_pol_init(_graphite_analyser25);

  class_Monochromator_pol_init(_graphite_analyser26);

  class_Monochromator_pol_init(_graphite_analyser27);

  class_Monochromator_pol_init(_graphite_analyser28);

  class_Monochromator_pol_init(_graphite_analyser29);

  class_Monochromator_pol_init(_graphite_analyser30);

  class_Monochromator_pol_init(_graphite_analyser31);

  class_Monochromator_pol_init(_graphite_analyser32);

  class_Monochromator_pol_init(_graphite_analyser33);

  class_Monochromator_pol_init(_graphite_analyser34);

  class_Monochromator_pol_init(_graphite_analyser35);

  class_Monochromator_pol_init(_graphite_analyser36);

  class_Monochromator_pol_init(_graphite_analyser37);

  class_Monochromator_pol_init(_graphite_analyser38);

  class_Monochromator_pol_init(_graphite_analyser39);



  class_L_monitor_init(_lamDetector);

  class_PSD_monitor_init(_psdDetector);

  class_TOF_monitor_init(_tofDetector);

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
_class_Arm *class_Arm_trace(_class_Arm *_comp
  , _class_particle *_particle) {
  ABSORBED=SCATTERED=RESTORE=0;


  // EXTEND code here
  if (!strcmp(_comp->_name, "armTarget")) {
  probTarget = rand01();
  }
  if (!strcmp(_comp->_name, "armAnalyzer")) {
  // Multiply probability by 3 because we have 3 targets
  p*=3.0;
  }

  return(_comp);
} /* class_Arm_trace */

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
_class_Guide_curved *class_Guide_curved_trace(_class_Guide_curved *_comp
  , _class_particle *_particle) {
  ABSORBED=SCATTERED=RESTORE=0;

  #define w1 (_comp->_parameters.w1)
  #define h1 (_comp->_parameters.h1)
  #define l (_comp->_parameters.l)
  #define R0 (_comp->_parameters.R0)
  #define Qc (_comp->_parameters.Qc)
  #define alpha (_comp->_parameters.alpha)
  #define m (_comp->_parameters.m)
  #define W (_comp->_parameters.W)
  #define curvature (_comp->_parameters.curvature)
  double t11, t12, t21, t22, theta, alphaAng, endtime, phi;
  double time, time1, time2, q, R;
  int ii, i_bounce;

  double whalf  = 0.5*w1, hhalf = 0.5*h1;   /* half width and height of guide */
  double z_off  = curvature*sin(l/curvature);       /* z-component of total guide length */
  double R1     = curvature - whalf;        /* radius of curvature of inside mirror */
  double R2     = curvature + whalf;        /* radius of curvature of outside mirror */
  double vel    = sqrt(vx*vx + vy*vy + vz*vz);  /* neutron velocity */
  double vel_xz = sqrt(vx*vx + vz*vz);      /* in plane velocity */
  double K      = V2K*vel;        /* neutron wavevector */
  double lambda = 2.0*PI/K;       /* neutron wavelength */

/* Propagate neutron to guide entrance. */

  PROP_Z0;
  if(x <= -whalf || x >= whalf || y <= -hhalf || y >= hhalf)
    ABSORB;
  SCATTER;
  for(;;)
  {
    double par[]={R0, Qc, alpha, m, W};
    /* Find itersection points of neutron with inside and outside guide walls */
    ii = cylinder_intersect(&t11, &t12 ,x - curvature, y, z, vx, vy, vz, R1, h1);
    ii = cylinder_intersect(&t21, &t22 ,x - curvature, y, z, vx, vy, vz, R2, h1);

    /* Choose appropriate reflection time */
    time1 = (t11 < 1e-7) ? t12 : t11;
    time2 = (t21 < 1e-7) ? t22 : t21;
    time  = (time1 < 1e-7 || time2 < time1) ? time2 : time1;

    /* Has neutron left the guide? */
    endtime = (z_off - z)/vz;
    if (time > endtime || time <= 1e-7) break;

    PROP_DT(time);

    /* Find reflection surface */
    R = (time == time1) ? R1 : R2;
    i_bounce = (fabs(y - hhalf) < 1e-7 || fabs(y + hhalf) < 1e-7) ? 2 : 1;
    switch(i_bounce) {
    case 1:           /* Inside or Outside wall */
      phi   = atan(vx/vz);        /* angle of neutron trajectory */
      alphaAng = asin(z/R);      /* angle of guide wall */
      theta = fabs(phi - alphaAng);    /* angle of reflection */
              vz    = vel_xz*cos(2.0*alphaAng - phi);
      vx    = vel_xz*sin(2.0*alphaAng - phi);
      break;
    case 2:       /* Top or Bottom wall */
      theta = fabs(atan(vy/vz));
      vy    = -vy;
      break;
    }
    /* Now compute reflectivity. */
    if (m == 0 || !R0) ABSORB;

    q = 4.0*PI*sin(theta)/lambda;
    StdReflecFunc(q, par, &R);
    if (R >= 0) p *= R; else ABSORB;
    SCATTER;
  }
  #undef w1
  #undef h1
  #undef l
  #undef R0
  #undef Qc
  #undef alpha
  #undef m
  #undef W
  #undef curvature
  return(_comp);
} /* class_Guide_curved_trace */

#pragma acc routine seq
_class_Incoherent *class_Incoherent_trace(_class_Incoherent *_comp
  , _class_particle *_particle) {
  ABSORBED=SCATTERED=RESTORE=0;

  #define geometry (_comp->_parameters.geometry)
  #define radius (_comp->_parameters.radius)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define zdepth (_comp->_parameters.zdepth)
  #define thickness (_comp->_parameters.thickness)
  #define target_x (_comp->_parameters.target_x)
  #define target_y (_comp->_parameters.target_y)
  #define target_z (_comp->_parameters.target_z)
  #define focus_r (_comp->_parameters.focus_r)
  #define focus_xw (_comp->_parameters.focus_xw)
  #define focus_yh (_comp->_parameters.focus_yh)
  #define focus_aw (_comp->_parameters.focus_aw)
  #define focus_ah (_comp->_parameters.focus_ah)
  #define target_index (_comp->_parameters.target_index)
  #define pack (_comp->_parameters.pack)
  #define p_interact (_comp->_parameters.p_interact)
  #define f_QE (_comp->_parameters.f_QE)
  #define gamma (_comp->_parameters.gamma)
  #define sigma_abs (_comp->_parameters.sigma_abs)
  #define sigma_inc (_comp->_parameters.sigma_inc)
  #define Vc (_comp->_parameters.Vc)
  #define concentric (_comp->_parameters.concentric)
  #define order (_comp->_parameters.order)
  #define VarsInc (_comp->_parameters.VarsInc)
  #define offdata (_comp->_parameters.offdata)
  double t0, t3;                /* Entry/exit time for outer surface */
  double t1, t2;                /* Entry/exit time for inner surface */
  double dt0, dt1, dt2, dt;     /* Flight times through sample */
  double v=0;                   /* Neutron velocity */
  double d_path;                /* Flight path length for non-scattered neutron */
  double l_i, l_o=0;            /* Flight path lenght in/out for scattered neutron */
  double my_a=0,my_t=0;         /* Velocity-dependent attenuation factor and total Xsec */
  double solid_angle=0;         /* Solid angle of target as seen from scattering point */
  double aim_x=0, aim_y=0, aim_z=1;   /* Position of target relative to scattering point */
  double v_i, v_f, E_i, E_f; /* initial and final energies and velocities */
  double dE;                 /* Energy transfer */
  int    intersect=0;
  int    flag_concentric=0;
  int    flag=0;
  double mc_trans, p_trans, mc_scatt, p_scatt, ws;
  double p_mult=1;

  do { /* Main interaction loop. Ends with intersect=0 */

    /* Intersection neutron trajectory / sample (sample surface) */
    if (VarsInc.shape == 0)
      intersect = cylinder_intersect(&t0, &t3, x, y, z, vx, vy, vz, radius, yheight);
    else if (VarsInc.shape == 1)
      intersect = box_intersect(&t0, &t3, x, y, z, vx, vy, vz, xwidth, yheight, zdepth);
    else if (VarsInc.shape == 2)
      intersect = sphere_intersect(&t0, &t3, x, y, z, vx, vy, vz, radius);
    else if (VarsInc.shape == 3)
      intersect = off_intersect(&t0, &t3, NULL, NULL, x, y, z, vx, vy, vz, offdata );
    if (intersect) {
      int flag_ishollow = 0;
      if (thickness>0) {
        if (VarsInc.shape==0 && cylinder_intersect(&t1,&t2, x,y,z,vx,vy,vz, radius-thickness,yheight-2*thickness))
          flag_ishollow=1;
        else if (VarsInc.shape==2 && sphere_intersect   (&t1,&t2, x,y,z,vx,vy,vz, radius-thickness))
          flag_ishollow=1;
        else if (VarsInc.shape==1 && box_intersect(&t1,&t2, x,y,z,vx,vy,vz, xwidth-2*thickness, yheight-2*thickness, zdepth-2*thickness))
          flag_ishollow = 1;
      }
      if (!flag_ishollow) t1 = t2 = t3; /* no empty space inside */

      dt0 = t1-t0;                /* Time in sample, ingoing */
      dt1 = t2-t1;                /* Time in hole */
      dt2 = t3-t2;                /* Time in sample, outgoing */

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
       && VarsInc.shape==0 && thickness>0) {
        flag_concentric=1;
      }

      if (flag_concentric == 1) {
        dt1=dt2=0; /* force exit when reaching hole/2nd part */
      }

      if (!dt0 && !dt2) {
        intersect = 0; /* the sample was passed entirely */
        break;
      }

      p_mult = 1;
      if (!v) v = sqrt(vx*vx + vy*vy + vz*vz);
      if (v) my_a = VarsInc.my_a_v*(2200/v);
      else {
        printf("Incoherent: %s: ERROR: Null velocity\n",NAME_CURRENT_COMP);
        ABSORB; /* should never occur */
      }

      my_t = my_a + VarsInc.my_s;  /* total scattering Xsect (tmp var) */
      if (my_t <= 0) {
        printf("Incoherent: %s: ERROR: Null total cross section %g. Removing event.\n",
          NAME_CURRENT_COMP, my_t);
        ABSORB; /* should never occur */
      }
      d_path = v * (dt0 + dt2);   /* Length of full path through sample */
      /* Proba of scattering vs absorption (integrating along the whole trajectory) */
      ws = VarsInc.my_s/my_t;  /* (inc+coh)/(inc+coh+abs) */
      /* Proba of transmission along length d_path */
      p_trans = exp(-my_t*d_path);
      p_scatt = 1 - p_trans; /* portion of beam which scatters */
      flag = 0; /* flag used for propagation to exit point before ending */
      /* are we next to the exit ? probably no scattering (avoid rounding errors) */
      if (VarsInc.my_s*d_path <= 4e-7) {
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
      if (!flag && mc_scatt > 0 && (mc_scatt >= 1 || (rand01()) < mc_scatt)) { /* Interaction neutron/sample */
        p_mult *= ws; /* Update weight ; account for absorption and retain scattered fraction */
        if (!mc_scatt) ABSORB;
        /* we have chosen portion mc_scatt of beam instead of p_scatt, so we compensate */
        p_mult *= fabs(p_scatt/mc_scatt); /* lower than 1 */
      } else {
        flag = 1; /* Transmission : no interaction neutron/sample */
        if (!mc_trans) ABSORB;
        p_mult *= fabs(p_trans/mc_trans);  /* attenuate beam by portion which is scattered (and left along) */
      }

      if (flag) { /* propagate to exit of sample and finish */
        intersect = 0;
        p *= p_mult; /* apply absorption correction */
        PROP_DT(dt0+dt2);
        break; /* exit main multi scatt while loop */
      }
      if (my_t*d_path < 1e-6)
      /* For very weak scattering, use simple uniform sampling of scattering
         point to avoid rounding errors. */
        dt = rand0max(d_path); /* length */
      else
        dt = -log(1 - rand0max((1 - exp(-my_t*d_path)))) / my_t; /* length */
      l_i = dt;/* Penetration in sample: scattering+abs */
      dt /= v; /* Time from present position to scattering point */

      /* If t0 is in hole, propagate to next part of the hollow cylinder */
      if (dt1 > 0 && dt0 > 0 && dt > dt0) dt += dt1;

      PROP_DT(dt);                /* Point of scattering */

      if ((VarsInc.tx || VarsInc.ty || VarsInc.tz)) {
        aim_x = VarsInc.tx-x;       /* Vector pointing at target (anal./det.) */
        aim_y = VarsInc.ty-y;
        aim_z = VarsInc.tz-z;
      }
      if(VarsInc.aw && VarsInc.ah) {
        randvec_target_rect_angular(&vx, &vy, &vz, &solid_angle,
          aim_x, aim_y, aim_z, VarsInc.aw, VarsInc.ah, ROT_A_CURRENT_COMP);
      } else if(VarsInc.xw && VarsInc.yh) {
        randvec_target_rect(&vx, &vy, &vz, &solid_angle,
          aim_x, aim_y, aim_z, VarsInc.xw, VarsInc.yh, ROT_A_CURRENT_COMP);
      } else {
        randvec_target_circle(&vx, &vy, &vz, &solid_angle, aim_x, aim_y, aim_z, focus_r);
      }
      NORM(vx, vy, vz);

      v_i = v;          /* Store initial velocity in case of quasielastic */
      if (rand01()<f_QE)	/* Quasielastic contribution */
      {
          E_i = VS2E*v_i*v_i;
          dE = gamma*tan(PI/2*randpm1());
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

      /* We do not consider scattering from 2nd part (outgoing) */
      p_mult *= solid_angle/4/PI;
      p *= p_mult;

      /* Polarisation part (1/3 NSF, 2/3 SF) */
      sx *= -1.0/3.0;
      sy *= -1.0/3.0;
      sz *= -1.0/3.0;

      SCATTER;

      /* test for a given multiple order */
      if (order && SCATTERED >= order) {
        intersect=0; /* reached required number of SCATTERing */
        break;       /* finish multiple scattering loop */
      }
    } /* end if intersect */
  } while (intersect); /* end do (intersect) (multiple scattering loop) */
  #undef geometry
  #undef radius
  #undef xwidth
  #undef yheight
  #undef zdepth
  #undef thickness
  #undef target_x
  #undef target_y
  #undef target_z
  #undef focus_r
  #undef focus_xw
  #undef focus_yh
  #undef focus_aw
  #undef focus_ah
  #undef target_index
  #undef pack
  #undef p_interact
  #undef f_QE
  #undef gamma
  #undef sigma_abs
  #undef sigma_inc
  #undef Vc
  #undef concentric
  #undef order
  #undef VarsInc
  #undef offdata
  return(_comp);
} /* class_Incoherent_trace */

#pragma acc routine seq
_class_Beamstop *class_Beamstop_trace(_class_Beamstop *_comp
  , _class_particle *_particle) {
  ABSORBED=SCATTERED=RESTORE=0;

  #define xmin (_comp->_parameters.xmin)
  #define xmax (_comp->_parameters.xmax)
  #define ymin (_comp->_parameters.ymin)
  #define ymax (_comp->_parameters.ymax)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define radius (_comp->_parameters.radius)
    double Time = t;
    ALLOW_BACKPROP;
    PROP_Z0;
    Time = t - Time;
    if ((Time>=0) && ((radius!=0) && (x*x + y*y <= radius*radius))
    || ((Time>=0) && (radius==0) && (x>xmin && x<xmax && y>ymin && y<ymax)))
      ABSORB;
    else
      RESTORE_NEUTRON(INDEX_CURRENT_COMP, x, y, z, vx, vy, vz, t, sx, sy, sz, p);
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef xwidth
  #undef yheight
  #undef radius
  return(_comp);
} /* class_Beamstop_trace */

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


  // EXTEND code here
  if (!strcmp(_comp->_name, "graphite_analyser0")) {
   groupNumber = 0;
  }
  if (!strcmp(_comp->_name, "graphite_analyser1")) {
   groupNumber = 0;
  }
  if (!strcmp(_comp->_name, "graphite_analyser2")) {
   groupNumber = 0;
  }
  if (!strcmp(_comp->_name, "graphite_analyser3")) {
   groupNumber = 0;
  }
  if (!strcmp(_comp->_name, "graphite_analyser4")) {
   groupNumber = 0;
  }
  if (!strcmp(_comp->_name, "graphite_analyser5")) {
   groupNumber = 0;
  }
  if (!strcmp(_comp->_name, "graphite_analyser6")) {
   groupNumber = 0;
  }
  if (!strcmp(_comp->_name, "graphite_analyser7")) {
   groupNumber = 0;
  }
  if (!strcmp(_comp->_name, "graphite_analyser8")) {
   groupNumber = 0;
  }
  if (!strcmp(_comp->_name, "graphite_analyser9")) {
   groupNumber = 0;
  }
  if (!strcmp(_comp->_name, "graphite_analyser10")) {
   groupNumber = 1;
  }
  if (!strcmp(_comp->_name, "graphite_analyser11")) {
   groupNumber = 1;
  }
  if (!strcmp(_comp->_name, "graphite_analyser12")) {
   groupNumber = 1;
  }
  if (!strcmp(_comp->_name, "graphite_analyser13")) {
   groupNumber = 1;
  }
  if (!strcmp(_comp->_name, "graphite_analyser14")) {
   groupNumber = 1;
  }
  if (!strcmp(_comp->_name, "graphite_analyser15")) {
   groupNumber = 1;
  }
  if (!strcmp(_comp->_name, "graphite_analyser16")) {
   groupNumber = 1;
  }
  if (!strcmp(_comp->_name, "graphite_analyser17")) {
   groupNumber = 1;
  }
  if (!strcmp(_comp->_name, "graphite_analyser18")) {
   groupNumber = 1;
  }
  if (!strcmp(_comp->_name, "graphite_analyser19")) {
   groupNumber = 1;
  }
  if (!strcmp(_comp->_name, "graphite_analyser20")) {
   groupNumber = 2;
  }
  if (!strcmp(_comp->_name, "graphite_analyser21")) {
   groupNumber = 2;
  }
  if (!strcmp(_comp->_name, "graphite_analyser22")) {
   groupNumber = 2;
  }
  if (!strcmp(_comp->_name, "graphite_analyser23")) {
   groupNumber = 2;
  }
  if (!strcmp(_comp->_name, "graphite_analyser24")) {
   groupNumber = 2;
  }
  if (!strcmp(_comp->_name, "graphite_analyser25")) {
   groupNumber = 2;
  }
  if (!strcmp(_comp->_name, "graphite_analyser26")) {
   groupNumber = 2;
  }
  if (!strcmp(_comp->_name, "graphite_analyser27")) {
   groupNumber = 2;
  }
  if (!strcmp(_comp->_name, "graphite_analyser28")) {
   groupNumber = 2;
  }
  if (!strcmp(_comp->_name, "graphite_analyser29")) {
   groupNumber = 2;
  }
  if (!strcmp(_comp->_name, "graphite_analyser30")) {
   groupNumber = 3;
  }
  if (!strcmp(_comp->_name, "graphite_analyser31")) {
   groupNumber = 3;
  }
  if (!strcmp(_comp->_name, "graphite_analyser32")) {
   groupNumber = 3;
  }
  if (!strcmp(_comp->_name, "graphite_analyser33")) {
   groupNumber = 3;
  }
  if (!strcmp(_comp->_name, "graphite_analyser34")) {
   groupNumber = 3;
  }
  if (!strcmp(_comp->_name, "graphite_analyser35")) {
   groupNumber = 3;
  }
  if (!strcmp(_comp->_name, "graphite_analyser36")) {
   groupNumber = 3;
  }
  if (!strcmp(_comp->_name, "graphite_analyser37")) {
   groupNumber = 3;
  }
  if (!strcmp(_comp->_name, "graphite_analyser38")) {
   groupNumber = 3;
  }
  if (!strcmp(_comp->_name, "graphite_analyser39")) {
   groupNumber = 3;
  }

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

/* *****************************************************************************
* instrument 'ISIS_OSIRIS' TRACE
***************************************************************************** */

#pragma acc routine seq
int raytrace(_class_particle* _particle) { /* called by mccode_main for ISIS_OSIRIS:TRACE */

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
      /* component armStart=Arm() [1] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_armStart->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _armStart->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_armStart->_position_relative, _armStart->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      _particle->_index++;
    } /* end component armStart [1] */
    if (!ABSORBED && _particle->_index == 2) {
      /* component isis_mod=ISIS_moderator() [2] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_isis_mod->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _isis_mod->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_isis_mod->_position_relative, _isis_mod->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_ISIS_moderator_trace(_isis_mod, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component isis_mod [2] */
    if (!ABSORBED && _particle->_index == 3) {
      /* component tofSource=TOF_monitor() [3] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_tofSource->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _tofSource->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_tofSource->_position_relative, _tofSource->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_TOF_monitor_trace(_tofSource, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component tofSource [3] */
    if (!ABSORBED && _particle->_index == 4) {
      /* component armStraight1=Arm() [4] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_armStraight1->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _armStraight1->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_armStraight1->_position_relative, _armStraight1->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      _particle->_index++;
    } /* end component armStraight1 [4] */
    if (!ABSORBED && _particle->_index == 5) {
      /* component lamStart=L_monitor() [5] */
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
    } /* end component lamStart [5] */
    if (!ABSORBED && _particle->_index == 6) {
      /* component psdStart=PSD_monitor() [6] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_psdStart->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _psdStart->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_psdStart->_position_relative, _psdStart->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_PSD_monitor_trace(_psdStart, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component psdStart [6] */
    if (!ABSORBED && _particle->_index == 7) {
      /* component guide_straight_1=Guide() [7] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_guide_straight_1->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _guide_straight_1->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_guide_straight_1->_position_relative, _guide_straight_1->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Guide_trace(_guide_straight_1, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component guide_straight_1 [7] */
    if (!ABSORBED && _particle->_index == 8) {
      /* component armChopper1=Arm() [8] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_armChopper1->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _armChopper1->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_armChopper1->_position_relative, _armChopper1->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      _particle->_index++;
    } /* end component armChopper1 [8] */
    if (!ABSORBED && _particle->_index == 9) {
      /* component chopper1=DiskChopper() [9] */
  #define theta_0 (_chopper1->_parameters.theta_0)
  #define radius (_chopper1->_parameters.radius)
  #define yheight (_chopper1->_parameters.yheight)
  #define nu (_chopper1->_parameters.nu)
  #define nslit (_chopper1->_parameters.nslit)
  #define jitter (_chopper1->_parameters.jitter)
  #define delay (_chopper1->_parameters.delay)
  #define isfirst (_chopper1->_parameters.isfirst)
  #define n_pulse (_chopper1->_parameters.n_pulse)
  #define abs_out (_chopper1->_parameters.abs_out)
  #define phase (_chopper1->_parameters.phase)
  #define xwidth (_chopper1->_parameters.xwidth)
  #define verbose (_chopper1->_parameters.verbose)
  #define Tg (_chopper1->_parameters.Tg)
  #define To (_chopper1->_parameters.To)
  #define delta_y (_chopper1->_parameters.delta_y)
  #define height (_chopper1->_parameters.height)
  #define omega (_chopper1->_parameters.omega)
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_chopper1->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _chopper1->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_chopper1->_position_relative, _chopper1->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      if (( instrument->_parameters._LAMBDA > 0 )) // conditional WHEN execution
      class_DiskChopper_trace(_chopper1, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
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
      _particle->_index++;
    } /* end component chopper1 [9] */
    if (!ABSORBED && _particle->_index == 10) {
      /* component armCurved1=Arm() [10] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_armCurved1->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _armCurved1->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_armCurved1->_position_relative, _armCurved1->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      _particle->_index++;
    } /* end component armCurved1 [10] */
    if (!ABSORBED && _particle->_index == 11) {
      /* component guide_curved1=Guide_curved() [11] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_guide_curved1->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _guide_curved1->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_guide_curved1->_position_relative, _guide_curved1->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Guide_curved_trace(_guide_curved1, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component guide_curved1 [11] */
    if (!ABSORBED && _particle->_index == 12) {
      /* component armChopper2=Arm() [12] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_armChopper2->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _armChopper2->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_armChopper2->_position_relative, _armChopper2->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      _particle->_index++;
    } /* end component armChopper2 [12] */
    if (!ABSORBED && _particle->_index == 13) {
      /* component chopper2=DiskChopper() [13] */
  #define theta_0 (_chopper2->_parameters.theta_0)
  #define radius (_chopper2->_parameters.radius)
  #define yheight (_chopper2->_parameters.yheight)
  #define nu (_chopper2->_parameters.nu)
  #define nslit (_chopper2->_parameters.nslit)
  #define jitter (_chopper2->_parameters.jitter)
  #define delay (_chopper2->_parameters.delay)
  #define isfirst (_chopper2->_parameters.isfirst)
  #define n_pulse (_chopper2->_parameters.n_pulse)
  #define abs_out (_chopper2->_parameters.abs_out)
  #define phase (_chopper2->_parameters.phase)
  #define xwidth (_chopper2->_parameters.xwidth)
  #define verbose (_chopper2->_parameters.verbose)
  #define Tg (_chopper2->_parameters.Tg)
  #define To (_chopper2->_parameters.To)
  #define delta_y (_chopper2->_parameters.delta_y)
  #define height (_chopper2->_parameters.height)
  #define omega (_chopper2->_parameters.omega)
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_chopper2->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _chopper2->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_chopper2->_position_relative, _chopper2->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      if (( instrument->_parameters._LAMBDA > 0 )) // conditional WHEN execution
      class_DiskChopper_trace(_chopper2, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
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
      _particle->_index++;
    } /* end component chopper2 [13] */
    if (!ABSORBED && _particle->_index == 14) {
      /* component armCurved2=Arm() [14] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_armCurved2->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _armCurved2->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_armCurved2->_position_relative, _armCurved2->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      _particle->_index++;
    } /* end component armCurved2 [14] */
    if (!ABSORBED && _particle->_index == 15) {
      /* component guide_curved2=Guide_curved() [15] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_guide_curved2->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _guide_curved2->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_guide_curved2->_position_relative, _guide_curved2->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Guide_curved_trace(_guide_curved2, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component guide_curved2 [15] */
    if (!ABSORBED && _particle->_index == 16) {
      /* component armStraight2=Arm() [16] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_armStraight2->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _armStraight2->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_armStraight2->_position_relative, _armStraight2->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      _particle->_index++;
    } /* end component armStraight2 [16] */
    if (!ABSORBED && _particle->_index == 17) {
      /* component guide_straight_2=Guide() [17] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_guide_straight_2->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _guide_straight_2->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_guide_straight_2->_position_relative, _guide_straight_2->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Guide_trace(_guide_straight_2, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component guide_straight_2 [17] */
    if (!ABSORBED && _particle->_index == 18) {
      /* component armSuper=Arm() [18] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_armSuper->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _armSuper->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_armSuper->_position_relative, _armSuper->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      _particle->_index++;
    } /* end component armSuper [18] */
    if (!ABSORBED && _particle->_index == 19) {
      /* component guide_super=Guide() [19] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_guide_super->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _guide_super->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_guide_super->_position_relative, _guide_super->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Guide_trace(_guide_super, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component guide_super [19] */
    if (!ABSORBED && _particle->_index == 20) {
      /* component armTarget=Arm() [20] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_armTarget->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _armTarget->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_armTarget->_position_relative, _armTarget->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Arm_trace(_armTarget, _particle); /* contains EXTEND code */
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component armTarget [20] */
    if (!ABSORBED && _particle->_index == 21) {
      /* component lamTarget=L_monitor() [21] */
  #define nL (_lamTarget->_parameters.nL)
  #define filename (_lamTarget->_parameters.filename)
  #define xmin (_lamTarget->_parameters.xmin)
  #define xmax (_lamTarget->_parameters.xmax)
  #define ymin (_lamTarget->_parameters.ymin)
  #define ymax (_lamTarget->_parameters.ymax)
  #define xwidth (_lamTarget->_parameters.xwidth)
  #define yheight (_lamTarget->_parameters.yheight)
  #define Lmin (_lamTarget->_parameters.Lmin)
  #define Lmax (_lamTarget->_parameters.Lmax)
  #define restore_neutron (_lamTarget->_parameters.restore_neutron)
  #define L_N (_lamTarget->_parameters.L_N)
  #define L_p (_lamTarget->_parameters.L_p)
  #define L_p2 (_lamTarget->_parameters.L_p2)
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_lamTarget->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _lamTarget->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_lamTarget->_position_relative, _lamTarget->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      if (( instrument->_parameters._LAMBDA == 0 )) // conditional WHEN execution
      class_L_monitor_trace(_lamTarget, _particle);
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
    } /* end component lamTarget [21] */
    if (!ABSORBED && _particle->_index == 22) {
      /* component psdTarget=PSD_monitor() [22] */
  #define nx (_psdTarget->_parameters.nx)
  #define ny (_psdTarget->_parameters.ny)
  #define filename (_psdTarget->_parameters.filename)
  #define xmin (_psdTarget->_parameters.xmin)
  #define xmax (_psdTarget->_parameters.xmax)
  #define ymin (_psdTarget->_parameters.ymin)
  #define ymax (_psdTarget->_parameters.ymax)
  #define xwidth (_psdTarget->_parameters.xwidth)
  #define yheight (_psdTarget->_parameters.yheight)
  #define restore_neutron (_psdTarget->_parameters.restore_neutron)
  #define PSD_N (_psdTarget->_parameters.PSD_N)
  #define PSD_p (_psdTarget->_parameters.PSD_p)
  #define PSD_p2 (_psdTarget->_parameters.PSD_p2)
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_psdTarget->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _psdTarget->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_psdTarget->_position_relative, _psdTarget->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      if (( instrument->_parameters._LAMBDA == 0 )) // conditional WHEN execution
      class_PSD_monitor_trace(_psdTarget, _particle);
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
    } /* end component psdTarget [22] */
    if (!ABSORBED && _particle->_index == 23) {
      /* component vsample1=Incoherent() [23] */
  #define geometry (_vsample1->_parameters.geometry)
  #define radius (_vsample1->_parameters.radius)
  #define xwidth (_vsample1->_parameters.xwidth)
  #define yheight (_vsample1->_parameters.yheight)
  #define zdepth (_vsample1->_parameters.zdepth)
  #define thickness (_vsample1->_parameters.thickness)
  #define target_x (_vsample1->_parameters.target_x)
  #define target_y (_vsample1->_parameters.target_y)
  #define target_z (_vsample1->_parameters.target_z)
  #define focus_r (_vsample1->_parameters.focus_r)
  #define focus_xw (_vsample1->_parameters.focus_xw)
  #define focus_yh (_vsample1->_parameters.focus_yh)
  #define focus_aw (_vsample1->_parameters.focus_aw)
  #define focus_ah (_vsample1->_parameters.focus_ah)
  #define target_index (_vsample1->_parameters.target_index)
  #define pack (_vsample1->_parameters.pack)
  #define p_interact (_vsample1->_parameters.p_interact)
  #define f_QE (_vsample1->_parameters.f_QE)
  #define gamma (_vsample1->_parameters.gamma)
  #define sigma_abs (_vsample1->_parameters.sigma_abs)
  #define sigma_inc (_vsample1->_parameters.sigma_inc)
  #define Vc (_vsample1->_parameters.Vc)
  #define concentric (_vsample1->_parameters.concentric)
  #define order (_vsample1->_parameters.order)
  #define VarsInc (_vsample1->_parameters.VarsInc)
  #define offdata (_vsample1->_parameters.offdata)
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_vsample1->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _vsample1->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_vsample1->_position_relative, _vsample1->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      if (( instrument->_parameters._LAMBDA > 0 && probTarget < 1.0 / 3.0 )) // conditional WHEN execution
      class_Incoherent_trace(_vsample1, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
  #undef geometry
  #undef radius
  #undef xwidth
  #undef yheight
  #undef zdepth
  #undef thickness
  #undef target_x
  #undef target_y
  #undef target_z
  #undef focus_r
  #undef focus_xw
  #undef focus_yh
  #undef focus_aw
  #undef focus_ah
  #undef target_index
  #undef pack
  #undef p_interact
  #undef f_QE
  #undef gamma
  #undef sigma_abs
  #undef sigma_inc
  #undef Vc
  #undef concentric
  #undef order
  #undef VarsInc
  #undef offdata
      _particle->_index++;
    } /* end component vsample1 [23] */
    if (!ABSORBED && _particle->_index == 24) {
      /* component vsample2=Incoherent() [24] */
  #define geometry (_vsample2->_parameters.geometry)
  #define radius (_vsample2->_parameters.radius)
  #define xwidth (_vsample2->_parameters.xwidth)
  #define yheight (_vsample2->_parameters.yheight)
  #define zdepth (_vsample2->_parameters.zdepth)
  #define thickness (_vsample2->_parameters.thickness)
  #define target_x (_vsample2->_parameters.target_x)
  #define target_y (_vsample2->_parameters.target_y)
  #define target_z (_vsample2->_parameters.target_z)
  #define focus_r (_vsample2->_parameters.focus_r)
  #define focus_xw (_vsample2->_parameters.focus_xw)
  #define focus_yh (_vsample2->_parameters.focus_yh)
  #define focus_aw (_vsample2->_parameters.focus_aw)
  #define focus_ah (_vsample2->_parameters.focus_ah)
  #define target_index (_vsample2->_parameters.target_index)
  #define pack (_vsample2->_parameters.pack)
  #define p_interact (_vsample2->_parameters.p_interact)
  #define f_QE (_vsample2->_parameters.f_QE)
  #define gamma (_vsample2->_parameters.gamma)
  #define sigma_abs (_vsample2->_parameters.sigma_abs)
  #define sigma_inc (_vsample2->_parameters.sigma_inc)
  #define Vc (_vsample2->_parameters.Vc)
  #define concentric (_vsample2->_parameters.concentric)
  #define order (_vsample2->_parameters.order)
  #define VarsInc (_vsample2->_parameters.VarsInc)
  #define offdata (_vsample2->_parameters.offdata)
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_vsample2->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _vsample2->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_vsample2->_position_relative, _vsample2->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      if (( instrument->_parameters._LAMBDA > 0 && probTarget >= 1.0 / 3.0 && probTarget < 2.0 / 3.0 )) // conditional WHEN execution
      class_Incoherent_trace(_vsample2, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
  #undef geometry
  #undef radius
  #undef xwidth
  #undef yheight
  #undef zdepth
  #undef thickness
  #undef target_x
  #undef target_y
  #undef target_z
  #undef focus_r
  #undef focus_xw
  #undef focus_yh
  #undef focus_aw
  #undef focus_ah
  #undef target_index
  #undef pack
  #undef p_interact
  #undef f_QE
  #undef gamma
  #undef sigma_abs
  #undef sigma_inc
  #undef Vc
  #undef concentric
  #undef order
  #undef VarsInc
  #undef offdata
      _particle->_index++;
    } /* end component vsample2 [24] */
    if (!ABSORBED && _particle->_index == 25) {
      /* component vsample3=Incoherent() [25] */
  #define geometry (_vsample3->_parameters.geometry)
  #define radius (_vsample3->_parameters.radius)
  #define xwidth (_vsample3->_parameters.xwidth)
  #define yheight (_vsample3->_parameters.yheight)
  #define zdepth (_vsample3->_parameters.zdepth)
  #define thickness (_vsample3->_parameters.thickness)
  #define target_x (_vsample3->_parameters.target_x)
  #define target_y (_vsample3->_parameters.target_y)
  #define target_z (_vsample3->_parameters.target_z)
  #define focus_r (_vsample3->_parameters.focus_r)
  #define focus_xw (_vsample3->_parameters.focus_xw)
  #define focus_yh (_vsample3->_parameters.focus_yh)
  #define focus_aw (_vsample3->_parameters.focus_aw)
  #define focus_ah (_vsample3->_parameters.focus_ah)
  #define target_index (_vsample3->_parameters.target_index)
  #define pack (_vsample3->_parameters.pack)
  #define p_interact (_vsample3->_parameters.p_interact)
  #define f_QE (_vsample3->_parameters.f_QE)
  #define gamma (_vsample3->_parameters.gamma)
  #define sigma_abs (_vsample3->_parameters.sigma_abs)
  #define sigma_inc (_vsample3->_parameters.sigma_inc)
  #define Vc (_vsample3->_parameters.Vc)
  #define concentric (_vsample3->_parameters.concentric)
  #define order (_vsample3->_parameters.order)
  #define VarsInc (_vsample3->_parameters.VarsInc)
  #define offdata (_vsample3->_parameters.offdata)
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_vsample3->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _vsample3->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_vsample3->_position_relative, _vsample3->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      if (( instrument->_parameters._LAMBDA > 0 && probTarget >= 2.0 / 3.0 )) // conditional WHEN execution
      class_Incoherent_trace(_vsample3, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
  #undef geometry
  #undef radius
  #undef xwidth
  #undef yheight
  #undef zdepth
  #undef thickness
  #undef target_x
  #undef target_y
  #undef target_z
  #undef focus_r
  #undef focus_xw
  #undef focus_yh
  #undef focus_aw
  #undef focus_ah
  #undef target_index
  #undef pack
  #undef p_interact
  #undef f_QE
  #undef gamma
  #undef sigma_abs
  #undef sigma_inc
  #undef Vc
  #undef concentric
  #undef order
  #undef VarsInc
  #undef offdata
      _particle->_index++;
    } /* end component vsample3 [25] */
    if (!ABSORBED && _particle->_index == 26) {
      /* component beamStop=Beamstop() [26] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_beamStop->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _beamStop->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_beamStop->_position_relative, _beamStop->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Beamstop_trace(_beamStop, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component beamStop [26] */
    if (!ABSORBED && _particle->_index == 27) {
      /* component armAnalyzer=Arm() [27] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_armAnalyzer->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _armAnalyzer->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_armAnalyzer->_position_relative, _armAnalyzer->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Arm_trace(_armAnalyzer, _particle); /* contains EXTEND code */
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component armAnalyzer [27] */
    if (!ABSORBED && _particle->_index == 28) {
      /* component graphite_analyser0=Monochromator_pol() [28] */
  #define zwidth (_graphite_analyser0->_parameters.zwidth)
  #define yheight (_graphite_analyser0->_parameters.yheight)
  #define mosaic (_graphite_analyser0->_parameters.mosaic)
  #define dspread (_graphite_analyser0->_parameters.dspread)
  #define Q (_graphite_analyser0->_parameters.Q)
  #define DM (_graphite_analyser0->_parameters.DM)
  #define pThreshold (_graphite_analyser0->_parameters.pThreshold)
  #define Rup (_graphite_analyser0->_parameters.Rup)
  #define Rdown (_graphite_analyser0->_parameters.Rdown)
  #define debug (_graphite_analyser0->_parameters.debug)
  #define mos_rms (_graphite_analyser0->_parameters.mos_rms)
  #define d_rms (_graphite_analyser0->_parameters.d_rms)
  #define mono_Q (_graphite_analyser0->_parameters.mono_Q)
  #define FN (_graphite_analyser0->_parameters.FN)
  #define FM (_graphite_analyser0->_parameters.FM)
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_graphite_analyser0->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _graphite_analyser0->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_graphite_analyser0->_position_relative, _graphite_analyser0->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      if (( instrument->_parameters._LAMBDA > 0 )) // conditional WHEN execution
      class_Monochromator_pol_trace(_graphite_analyser0, _particle); /* contains EXTEND code */
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
      // GROUP graphiteAnalyzer: from graphite_analyser0 [28] to graphite_analyser39 [67]
      if (SCATTERED) _particle->_index = 67; // when SCATTERED in GROUP: reach exit of GROUP after graphite_analyser39
      else ABSORBED=0; // not SCATTERED within GROUP: always tries next
      _particle->_index++;
    } /* end component graphite_analyser0 [28] */
    if (!ABSORBED && _particle->_index == 29) {
      /* component graphite_analyser1=Monochromator_pol() [29] */
  #define zwidth (_graphite_analyser1->_parameters.zwidth)
  #define yheight (_graphite_analyser1->_parameters.yheight)
  #define mosaic (_graphite_analyser1->_parameters.mosaic)
  #define dspread (_graphite_analyser1->_parameters.dspread)
  #define Q (_graphite_analyser1->_parameters.Q)
  #define DM (_graphite_analyser1->_parameters.DM)
  #define pThreshold (_graphite_analyser1->_parameters.pThreshold)
  #define Rup (_graphite_analyser1->_parameters.Rup)
  #define Rdown (_graphite_analyser1->_parameters.Rdown)
  #define debug (_graphite_analyser1->_parameters.debug)
  #define mos_rms (_graphite_analyser1->_parameters.mos_rms)
  #define d_rms (_graphite_analyser1->_parameters.d_rms)
  #define mono_Q (_graphite_analyser1->_parameters.mono_Q)
  #define FN (_graphite_analyser1->_parameters.FN)
  #define FM (_graphite_analyser1->_parameters.FM)
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_graphite_analyser1->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _graphite_analyser1->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_graphite_analyser1->_position_relative, _graphite_analyser1->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      if (( instrument->_parameters._LAMBDA > 0 )) // conditional WHEN execution
      class_Monochromator_pol_trace(_graphite_analyser1, _particle); /* contains EXTEND code */
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
      // GROUP graphiteAnalyzer: from graphite_analyser0 [28] to graphite_analyser39 [67]
      if (SCATTERED) _particle->_index = 67; // when SCATTERED in GROUP: reach exit of GROUP after graphite_analyser39
      else ABSORBED=0; // not SCATTERED within GROUP: always tries next
      _particle->_index++;
    } /* end component graphite_analyser1 [29] */
    if (!ABSORBED && _particle->_index == 30) {
      /* component graphite_analyser2=Monochromator_pol() [30] */
  #define zwidth (_graphite_analyser2->_parameters.zwidth)
  #define yheight (_graphite_analyser2->_parameters.yheight)
  #define mosaic (_graphite_analyser2->_parameters.mosaic)
  #define dspread (_graphite_analyser2->_parameters.dspread)
  #define Q (_graphite_analyser2->_parameters.Q)
  #define DM (_graphite_analyser2->_parameters.DM)
  #define pThreshold (_graphite_analyser2->_parameters.pThreshold)
  #define Rup (_graphite_analyser2->_parameters.Rup)
  #define Rdown (_graphite_analyser2->_parameters.Rdown)
  #define debug (_graphite_analyser2->_parameters.debug)
  #define mos_rms (_graphite_analyser2->_parameters.mos_rms)
  #define d_rms (_graphite_analyser2->_parameters.d_rms)
  #define mono_Q (_graphite_analyser2->_parameters.mono_Q)
  #define FN (_graphite_analyser2->_parameters.FN)
  #define FM (_graphite_analyser2->_parameters.FM)
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_graphite_analyser2->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _graphite_analyser2->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_graphite_analyser2->_position_relative, _graphite_analyser2->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      if (( instrument->_parameters._LAMBDA > 0 )) // conditional WHEN execution
      class_Monochromator_pol_trace(_graphite_analyser2, _particle); /* contains EXTEND code */
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
      // GROUP graphiteAnalyzer: from graphite_analyser0 [28] to graphite_analyser39 [67]
      if (SCATTERED) _particle->_index = 67; // when SCATTERED in GROUP: reach exit of GROUP after graphite_analyser39
      else ABSORBED=0; // not SCATTERED within GROUP: always tries next
      _particle->_index++;
    } /* end component graphite_analyser2 [30] */
    if (!ABSORBED && _particle->_index == 31) {
      /* component graphite_analyser3=Monochromator_pol() [31] */
  #define zwidth (_graphite_analyser3->_parameters.zwidth)
  #define yheight (_graphite_analyser3->_parameters.yheight)
  #define mosaic (_graphite_analyser3->_parameters.mosaic)
  #define dspread (_graphite_analyser3->_parameters.dspread)
  #define Q (_graphite_analyser3->_parameters.Q)
  #define DM (_graphite_analyser3->_parameters.DM)
  #define pThreshold (_graphite_analyser3->_parameters.pThreshold)
  #define Rup (_graphite_analyser3->_parameters.Rup)
  #define Rdown (_graphite_analyser3->_parameters.Rdown)
  #define debug (_graphite_analyser3->_parameters.debug)
  #define mos_rms (_graphite_analyser3->_parameters.mos_rms)
  #define d_rms (_graphite_analyser3->_parameters.d_rms)
  #define mono_Q (_graphite_analyser3->_parameters.mono_Q)
  #define FN (_graphite_analyser3->_parameters.FN)
  #define FM (_graphite_analyser3->_parameters.FM)
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_graphite_analyser3->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _graphite_analyser3->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_graphite_analyser3->_position_relative, _graphite_analyser3->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      if (( instrument->_parameters._LAMBDA > 0 )) // conditional WHEN execution
      class_Monochromator_pol_trace(_graphite_analyser3, _particle); /* contains EXTEND code */
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
      // GROUP graphiteAnalyzer: from graphite_analyser0 [28] to graphite_analyser39 [67]
      if (SCATTERED) _particle->_index = 67; // when SCATTERED in GROUP: reach exit of GROUP after graphite_analyser39
      else ABSORBED=0; // not SCATTERED within GROUP: always tries next
      _particle->_index++;
    } /* end component graphite_analyser3 [31] */
    if (!ABSORBED && _particle->_index == 32) {
      /* component graphite_analyser4=Monochromator_pol() [32] */
  #define zwidth (_graphite_analyser4->_parameters.zwidth)
  #define yheight (_graphite_analyser4->_parameters.yheight)
  #define mosaic (_graphite_analyser4->_parameters.mosaic)
  #define dspread (_graphite_analyser4->_parameters.dspread)
  #define Q (_graphite_analyser4->_parameters.Q)
  #define DM (_graphite_analyser4->_parameters.DM)
  #define pThreshold (_graphite_analyser4->_parameters.pThreshold)
  #define Rup (_graphite_analyser4->_parameters.Rup)
  #define Rdown (_graphite_analyser4->_parameters.Rdown)
  #define debug (_graphite_analyser4->_parameters.debug)
  #define mos_rms (_graphite_analyser4->_parameters.mos_rms)
  #define d_rms (_graphite_analyser4->_parameters.d_rms)
  #define mono_Q (_graphite_analyser4->_parameters.mono_Q)
  #define FN (_graphite_analyser4->_parameters.FN)
  #define FM (_graphite_analyser4->_parameters.FM)
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_graphite_analyser4->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _graphite_analyser4->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_graphite_analyser4->_position_relative, _graphite_analyser4->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      if (( instrument->_parameters._LAMBDA > 0 )) // conditional WHEN execution
      class_Monochromator_pol_trace(_graphite_analyser4, _particle); /* contains EXTEND code */
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
      // GROUP graphiteAnalyzer: from graphite_analyser0 [28] to graphite_analyser39 [67]
      if (SCATTERED) _particle->_index = 67; // when SCATTERED in GROUP: reach exit of GROUP after graphite_analyser39
      else ABSORBED=0; // not SCATTERED within GROUP: always tries next
      _particle->_index++;
    } /* end component graphite_analyser4 [32] */
    if (!ABSORBED && _particle->_index == 33) {
      /* component graphite_analyser5=Monochromator_pol() [33] */
  #define zwidth (_graphite_analyser5->_parameters.zwidth)
  #define yheight (_graphite_analyser5->_parameters.yheight)
  #define mosaic (_graphite_analyser5->_parameters.mosaic)
  #define dspread (_graphite_analyser5->_parameters.dspread)
  #define Q (_graphite_analyser5->_parameters.Q)
  #define DM (_graphite_analyser5->_parameters.DM)
  #define pThreshold (_graphite_analyser5->_parameters.pThreshold)
  #define Rup (_graphite_analyser5->_parameters.Rup)
  #define Rdown (_graphite_analyser5->_parameters.Rdown)
  #define debug (_graphite_analyser5->_parameters.debug)
  #define mos_rms (_graphite_analyser5->_parameters.mos_rms)
  #define d_rms (_graphite_analyser5->_parameters.d_rms)
  #define mono_Q (_graphite_analyser5->_parameters.mono_Q)
  #define FN (_graphite_analyser5->_parameters.FN)
  #define FM (_graphite_analyser5->_parameters.FM)
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_graphite_analyser5->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _graphite_analyser5->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_graphite_analyser5->_position_relative, _graphite_analyser5->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      if (( instrument->_parameters._LAMBDA > 0 )) // conditional WHEN execution
      class_Monochromator_pol_trace(_graphite_analyser5, _particle); /* contains EXTEND code */
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
      // GROUP graphiteAnalyzer: from graphite_analyser0 [28] to graphite_analyser39 [67]
      if (SCATTERED) _particle->_index = 67; // when SCATTERED in GROUP: reach exit of GROUP after graphite_analyser39
      else ABSORBED=0; // not SCATTERED within GROUP: always tries next
      _particle->_index++;
    } /* end component graphite_analyser5 [33] */
    if (!ABSORBED && _particle->_index == 34) {
      /* component graphite_analyser6=Monochromator_pol() [34] */
  #define zwidth (_graphite_analyser6->_parameters.zwidth)
  #define yheight (_graphite_analyser6->_parameters.yheight)
  #define mosaic (_graphite_analyser6->_parameters.mosaic)
  #define dspread (_graphite_analyser6->_parameters.dspread)
  #define Q (_graphite_analyser6->_parameters.Q)
  #define DM (_graphite_analyser6->_parameters.DM)
  #define pThreshold (_graphite_analyser6->_parameters.pThreshold)
  #define Rup (_graphite_analyser6->_parameters.Rup)
  #define Rdown (_graphite_analyser6->_parameters.Rdown)
  #define debug (_graphite_analyser6->_parameters.debug)
  #define mos_rms (_graphite_analyser6->_parameters.mos_rms)
  #define d_rms (_graphite_analyser6->_parameters.d_rms)
  #define mono_Q (_graphite_analyser6->_parameters.mono_Q)
  #define FN (_graphite_analyser6->_parameters.FN)
  #define FM (_graphite_analyser6->_parameters.FM)
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_graphite_analyser6->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _graphite_analyser6->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_graphite_analyser6->_position_relative, _graphite_analyser6->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      if (( instrument->_parameters._LAMBDA > 0 )) // conditional WHEN execution
      class_Monochromator_pol_trace(_graphite_analyser6, _particle); /* contains EXTEND code */
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
      // GROUP graphiteAnalyzer: from graphite_analyser0 [28] to graphite_analyser39 [67]
      if (SCATTERED) _particle->_index = 67; // when SCATTERED in GROUP: reach exit of GROUP after graphite_analyser39
      else ABSORBED=0; // not SCATTERED within GROUP: always tries next
      _particle->_index++;
    } /* end component graphite_analyser6 [34] */
    if (!ABSORBED && _particle->_index == 35) {
      /* component graphite_analyser7=Monochromator_pol() [35] */
  #define zwidth (_graphite_analyser7->_parameters.zwidth)
  #define yheight (_graphite_analyser7->_parameters.yheight)
  #define mosaic (_graphite_analyser7->_parameters.mosaic)
  #define dspread (_graphite_analyser7->_parameters.dspread)
  #define Q (_graphite_analyser7->_parameters.Q)
  #define DM (_graphite_analyser7->_parameters.DM)
  #define pThreshold (_graphite_analyser7->_parameters.pThreshold)
  #define Rup (_graphite_analyser7->_parameters.Rup)
  #define Rdown (_graphite_analyser7->_parameters.Rdown)
  #define debug (_graphite_analyser7->_parameters.debug)
  #define mos_rms (_graphite_analyser7->_parameters.mos_rms)
  #define d_rms (_graphite_analyser7->_parameters.d_rms)
  #define mono_Q (_graphite_analyser7->_parameters.mono_Q)
  #define FN (_graphite_analyser7->_parameters.FN)
  #define FM (_graphite_analyser7->_parameters.FM)
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_graphite_analyser7->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _graphite_analyser7->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_graphite_analyser7->_position_relative, _graphite_analyser7->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      if (( instrument->_parameters._LAMBDA > 0 )) // conditional WHEN execution
      class_Monochromator_pol_trace(_graphite_analyser7, _particle); /* contains EXTEND code */
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
      // GROUP graphiteAnalyzer: from graphite_analyser0 [28] to graphite_analyser39 [67]
      if (SCATTERED) _particle->_index = 67; // when SCATTERED in GROUP: reach exit of GROUP after graphite_analyser39
      else ABSORBED=0; // not SCATTERED within GROUP: always tries next
      _particle->_index++;
    } /* end component graphite_analyser7 [35] */
    if (!ABSORBED && _particle->_index == 36) {
      /* component graphite_analyser8=Monochromator_pol() [36] */
  #define zwidth (_graphite_analyser8->_parameters.zwidth)
  #define yheight (_graphite_analyser8->_parameters.yheight)
  #define mosaic (_graphite_analyser8->_parameters.mosaic)
  #define dspread (_graphite_analyser8->_parameters.dspread)
  #define Q (_graphite_analyser8->_parameters.Q)
  #define DM (_graphite_analyser8->_parameters.DM)
  #define pThreshold (_graphite_analyser8->_parameters.pThreshold)
  #define Rup (_graphite_analyser8->_parameters.Rup)
  #define Rdown (_graphite_analyser8->_parameters.Rdown)
  #define debug (_graphite_analyser8->_parameters.debug)
  #define mos_rms (_graphite_analyser8->_parameters.mos_rms)
  #define d_rms (_graphite_analyser8->_parameters.d_rms)
  #define mono_Q (_graphite_analyser8->_parameters.mono_Q)
  #define FN (_graphite_analyser8->_parameters.FN)
  #define FM (_graphite_analyser8->_parameters.FM)
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_graphite_analyser8->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _graphite_analyser8->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_graphite_analyser8->_position_relative, _graphite_analyser8->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      if (( instrument->_parameters._LAMBDA > 0 )) // conditional WHEN execution
      class_Monochromator_pol_trace(_graphite_analyser8, _particle); /* contains EXTEND code */
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
      // GROUP graphiteAnalyzer: from graphite_analyser0 [28] to graphite_analyser39 [67]
      if (SCATTERED) _particle->_index = 67; // when SCATTERED in GROUP: reach exit of GROUP after graphite_analyser39
      else ABSORBED=0; // not SCATTERED within GROUP: always tries next
      _particle->_index++;
    } /* end component graphite_analyser8 [36] */
    if (!ABSORBED && _particle->_index == 37) {
      /* component graphite_analyser9=Monochromator_pol() [37] */
  #define zwidth (_graphite_analyser9->_parameters.zwidth)
  #define yheight (_graphite_analyser9->_parameters.yheight)
  #define mosaic (_graphite_analyser9->_parameters.mosaic)
  #define dspread (_graphite_analyser9->_parameters.dspread)
  #define Q (_graphite_analyser9->_parameters.Q)
  #define DM (_graphite_analyser9->_parameters.DM)
  #define pThreshold (_graphite_analyser9->_parameters.pThreshold)
  #define Rup (_graphite_analyser9->_parameters.Rup)
  #define Rdown (_graphite_analyser9->_parameters.Rdown)
  #define debug (_graphite_analyser9->_parameters.debug)
  #define mos_rms (_graphite_analyser9->_parameters.mos_rms)
  #define d_rms (_graphite_analyser9->_parameters.d_rms)
  #define mono_Q (_graphite_analyser9->_parameters.mono_Q)
  #define FN (_graphite_analyser9->_parameters.FN)
  #define FM (_graphite_analyser9->_parameters.FM)
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_graphite_analyser9->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _graphite_analyser9->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_graphite_analyser9->_position_relative, _graphite_analyser9->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      if (( instrument->_parameters._LAMBDA > 0 )) // conditional WHEN execution
      class_Monochromator_pol_trace(_graphite_analyser9, _particle); /* contains EXTEND code */
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
      // GROUP graphiteAnalyzer: from graphite_analyser0 [28] to graphite_analyser39 [67]
      if (SCATTERED) _particle->_index = 67; // when SCATTERED in GROUP: reach exit of GROUP after graphite_analyser39
      else ABSORBED=0; // not SCATTERED within GROUP: always tries next
      _particle->_index++;
    } /* end component graphite_analyser9 [37] */
    if (!ABSORBED && _particle->_index == 38) {
      /* component graphite_analyser10=Monochromator_pol() [38] */
  #define zwidth (_graphite_analyser10->_parameters.zwidth)
  #define yheight (_graphite_analyser10->_parameters.yheight)
  #define mosaic (_graphite_analyser10->_parameters.mosaic)
  #define dspread (_graphite_analyser10->_parameters.dspread)
  #define Q (_graphite_analyser10->_parameters.Q)
  #define DM (_graphite_analyser10->_parameters.DM)
  #define pThreshold (_graphite_analyser10->_parameters.pThreshold)
  #define Rup (_graphite_analyser10->_parameters.Rup)
  #define Rdown (_graphite_analyser10->_parameters.Rdown)
  #define debug (_graphite_analyser10->_parameters.debug)
  #define mos_rms (_graphite_analyser10->_parameters.mos_rms)
  #define d_rms (_graphite_analyser10->_parameters.d_rms)
  #define mono_Q (_graphite_analyser10->_parameters.mono_Q)
  #define FN (_graphite_analyser10->_parameters.FN)
  #define FM (_graphite_analyser10->_parameters.FM)
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_graphite_analyser10->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _graphite_analyser10->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_graphite_analyser10->_position_relative, _graphite_analyser10->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      if (( instrument->_parameters._LAMBDA > 0 )) // conditional WHEN execution
      class_Monochromator_pol_trace(_graphite_analyser10, _particle); /* contains EXTEND code */
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
      // GROUP graphiteAnalyzer: from graphite_analyser0 [28] to graphite_analyser39 [67]
      if (SCATTERED) _particle->_index = 67; // when SCATTERED in GROUP: reach exit of GROUP after graphite_analyser39
      else ABSORBED=0; // not SCATTERED within GROUP: always tries next
      _particle->_index++;
    } /* end component graphite_analyser10 [38] */
    if (!ABSORBED && _particle->_index == 39) {
      /* component graphite_analyser11=Monochromator_pol() [39] */
  #define zwidth (_graphite_analyser11->_parameters.zwidth)
  #define yheight (_graphite_analyser11->_parameters.yheight)
  #define mosaic (_graphite_analyser11->_parameters.mosaic)
  #define dspread (_graphite_analyser11->_parameters.dspread)
  #define Q (_graphite_analyser11->_parameters.Q)
  #define DM (_graphite_analyser11->_parameters.DM)
  #define pThreshold (_graphite_analyser11->_parameters.pThreshold)
  #define Rup (_graphite_analyser11->_parameters.Rup)
  #define Rdown (_graphite_analyser11->_parameters.Rdown)
  #define debug (_graphite_analyser11->_parameters.debug)
  #define mos_rms (_graphite_analyser11->_parameters.mos_rms)
  #define d_rms (_graphite_analyser11->_parameters.d_rms)
  #define mono_Q (_graphite_analyser11->_parameters.mono_Q)
  #define FN (_graphite_analyser11->_parameters.FN)
  #define FM (_graphite_analyser11->_parameters.FM)
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_graphite_analyser11->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _graphite_analyser11->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_graphite_analyser11->_position_relative, _graphite_analyser11->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      if (( instrument->_parameters._LAMBDA > 0 )) // conditional WHEN execution
      class_Monochromator_pol_trace(_graphite_analyser11, _particle); /* contains EXTEND code */
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
      // GROUP graphiteAnalyzer: from graphite_analyser0 [28] to graphite_analyser39 [67]
      if (SCATTERED) _particle->_index = 67; // when SCATTERED in GROUP: reach exit of GROUP after graphite_analyser39
      else ABSORBED=0; // not SCATTERED within GROUP: always tries next
      _particle->_index++;
    } /* end component graphite_analyser11 [39] */
    if (!ABSORBED && _particle->_index == 40) {
      /* component graphite_analyser12=Monochromator_pol() [40] */
  #define zwidth (_graphite_analyser12->_parameters.zwidth)
  #define yheight (_graphite_analyser12->_parameters.yheight)
  #define mosaic (_graphite_analyser12->_parameters.mosaic)
  #define dspread (_graphite_analyser12->_parameters.dspread)
  #define Q (_graphite_analyser12->_parameters.Q)
  #define DM (_graphite_analyser12->_parameters.DM)
  #define pThreshold (_graphite_analyser12->_parameters.pThreshold)
  #define Rup (_graphite_analyser12->_parameters.Rup)
  #define Rdown (_graphite_analyser12->_parameters.Rdown)
  #define debug (_graphite_analyser12->_parameters.debug)
  #define mos_rms (_graphite_analyser12->_parameters.mos_rms)
  #define d_rms (_graphite_analyser12->_parameters.d_rms)
  #define mono_Q (_graphite_analyser12->_parameters.mono_Q)
  #define FN (_graphite_analyser12->_parameters.FN)
  #define FM (_graphite_analyser12->_parameters.FM)
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_graphite_analyser12->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _graphite_analyser12->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_graphite_analyser12->_position_relative, _graphite_analyser12->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      if (( instrument->_parameters._LAMBDA > 0 )) // conditional WHEN execution
      class_Monochromator_pol_trace(_graphite_analyser12, _particle); /* contains EXTEND code */
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
      // GROUP graphiteAnalyzer: from graphite_analyser0 [28] to graphite_analyser39 [67]
      if (SCATTERED) _particle->_index = 67; // when SCATTERED in GROUP: reach exit of GROUP after graphite_analyser39
      else ABSORBED=0; // not SCATTERED within GROUP: always tries next
      _particle->_index++;
    } /* end component graphite_analyser12 [40] */
    if (!ABSORBED && _particle->_index == 41) {
      /* component graphite_analyser13=Monochromator_pol() [41] */
  #define zwidth (_graphite_analyser13->_parameters.zwidth)
  #define yheight (_graphite_analyser13->_parameters.yheight)
  #define mosaic (_graphite_analyser13->_parameters.mosaic)
  #define dspread (_graphite_analyser13->_parameters.dspread)
  #define Q (_graphite_analyser13->_parameters.Q)
  #define DM (_graphite_analyser13->_parameters.DM)
  #define pThreshold (_graphite_analyser13->_parameters.pThreshold)
  #define Rup (_graphite_analyser13->_parameters.Rup)
  #define Rdown (_graphite_analyser13->_parameters.Rdown)
  #define debug (_graphite_analyser13->_parameters.debug)
  #define mos_rms (_graphite_analyser13->_parameters.mos_rms)
  #define d_rms (_graphite_analyser13->_parameters.d_rms)
  #define mono_Q (_graphite_analyser13->_parameters.mono_Q)
  #define FN (_graphite_analyser13->_parameters.FN)
  #define FM (_graphite_analyser13->_parameters.FM)
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_graphite_analyser13->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _graphite_analyser13->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_graphite_analyser13->_position_relative, _graphite_analyser13->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      if (( instrument->_parameters._LAMBDA > 0 )) // conditional WHEN execution
      class_Monochromator_pol_trace(_graphite_analyser13, _particle); /* contains EXTEND code */
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
      // GROUP graphiteAnalyzer: from graphite_analyser0 [28] to graphite_analyser39 [67]
      if (SCATTERED) _particle->_index = 67; // when SCATTERED in GROUP: reach exit of GROUP after graphite_analyser39
      else ABSORBED=0; // not SCATTERED within GROUP: always tries next
      _particle->_index++;
    } /* end component graphite_analyser13 [41] */
    if (!ABSORBED && _particle->_index == 42) {
      /* component graphite_analyser14=Monochromator_pol() [42] */
  #define zwidth (_graphite_analyser14->_parameters.zwidth)
  #define yheight (_graphite_analyser14->_parameters.yheight)
  #define mosaic (_graphite_analyser14->_parameters.mosaic)
  #define dspread (_graphite_analyser14->_parameters.dspread)
  #define Q (_graphite_analyser14->_parameters.Q)
  #define DM (_graphite_analyser14->_parameters.DM)
  #define pThreshold (_graphite_analyser14->_parameters.pThreshold)
  #define Rup (_graphite_analyser14->_parameters.Rup)
  #define Rdown (_graphite_analyser14->_parameters.Rdown)
  #define debug (_graphite_analyser14->_parameters.debug)
  #define mos_rms (_graphite_analyser14->_parameters.mos_rms)
  #define d_rms (_graphite_analyser14->_parameters.d_rms)
  #define mono_Q (_graphite_analyser14->_parameters.mono_Q)
  #define FN (_graphite_analyser14->_parameters.FN)
  #define FM (_graphite_analyser14->_parameters.FM)
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_graphite_analyser14->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _graphite_analyser14->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_graphite_analyser14->_position_relative, _graphite_analyser14->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      if (( instrument->_parameters._LAMBDA > 0 )) // conditional WHEN execution
      class_Monochromator_pol_trace(_graphite_analyser14, _particle); /* contains EXTEND code */
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
      // GROUP graphiteAnalyzer: from graphite_analyser0 [28] to graphite_analyser39 [67]
      if (SCATTERED) _particle->_index = 67; // when SCATTERED in GROUP: reach exit of GROUP after graphite_analyser39
      else ABSORBED=0; // not SCATTERED within GROUP: always tries next
      _particle->_index++;
    } /* end component graphite_analyser14 [42] */
    if (!ABSORBED && _particle->_index == 43) {
      /* component graphite_analyser15=Monochromator_pol() [43] */
  #define zwidth (_graphite_analyser15->_parameters.zwidth)
  #define yheight (_graphite_analyser15->_parameters.yheight)
  #define mosaic (_graphite_analyser15->_parameters.mosaic)
  #define dspread (_graphite_analyser15->_parameters.dspread)
  #define Q (_graphite_analyser15->_parameters.Q)
  #define DM (_graphite_analyser15->_parameters.DM)
  #define pThreshold (_graphite_analyser15->_parameters.pThreshold)
  #define Rup (_graphite_analyser15->_parameters.Rup)
  #define Rdown (_graphite_analyser15->_parameters.Rdown)
  #define debug (_graphite_analyser15->_parameters.debug)
  #define mos_rms (_graphite_analyser15->_parameters.mos_rms)
  #define d_rms (_graphite_analyser15->_parameters.d_rms)
  #define mono_Q (_graphite_analyser15->_parameters.mono_Q)
  #define FN (_graphite_analyser15->_parameters.FN)
  #define FM (_graphite_analyser15->_parameters.FM)
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_graphite_analyser15->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _graphite_analyser15->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_graphite_analyser15->_position_relative, _graphite_analyser15->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      if (( instrument->_parameters._LAMBDA > 0 )) // conditional WHEN execution
      class_Monochromator_pol_trace(_graphite_analyser15, _particle); /* contains EXTEND code */
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
      // GROUP graphiteAnalyzer: from graphite_analyser0 [28] to graphite_analyser39 [67]
      if (SCATTERED) _particle->_index = 67; // when SCATTERED in GROUP: reach exit of GROUP after graphite_analyser39
      else ABSORBED=0; // not SCATTERED within GROUP: always tries next
      _particle->_index++;
    } /* end component graphite_analyser15 [43] */
    if (!ABSORBED && _particle->_index == 44) {
      /* component graphite_analyser16=Monochromator_pol() [44] */
  #define zwidth (_graphite_analyser16->_parameters.zwidth)
  #define yheight (_graphite_analyser16->_parameters.yheight)
  #define mosaic (_graphite_analyser16->_parameters.mosaic)
  #define dspread (_graphite_analyser16->_parameters.dspread)
  #define Q (_graphite_analyser16->_parameters.Q)
  #define DM (_graphite_analyser16->_parameters.DM)
  #define pThreshold (_graphite_analyser16->_parameters.pThreshold)
  #define Rup (_graphite_analyser16->_parameters.Rup)
  #define Rdown (_graphite_analyser16->_parameters.Rdown)
  #define debug (_graphite_analyser16->_parameters.debug)
  #define mos_rms (_graphite_analyser16->_parameters.mos_rms)
  #define d_rms (_graphite_analyser16->_parameters.d_rms)
  #define mono_Q (_graphite_analyser16->_parameters.mono_Q)
  #define FN (_graphite_analyser16->_parameters.FN)
  #define FM (_graphite_analyser16->_parameters.FM)
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_graphite_analyser16->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _graphite_analyser16->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_graphite_analyser16->_position_relative, _graphite_analyser16->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      if (( instrument->_parameters._LAMBDA > 0 )) // conditional WHEN execution
      class_Monochromator_pol_trace(_graphite_analyser16, _particle); /* contains EXTEND code */
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
      // GROUP graphiteAnalyzer: from graphite_analyser0 [28] to graphite_analyser39 [67]
      if (SCATTERED) _particle->_index = 67; // when SCATTERED in GROUP: reach exit of GROUP after graphite_analyser39
      else ABSORBED=0; // not SCATTERED within GROUP: always tries next
      _particle->_index++;
    } /* end component graphite_analyser16 [44] */
    if (!ABSORBED && _particle->_index == 45) {
      /* component graphite_analyser17=Monochromator_pol() [45] */
  #define zwidth (_graphite_analyser17->_parameters.zwidth)
  #define yheight (_graphite_analyser17->_parameters.yheight)
  #define mosaic (_graphite_analyser17->_parameters.mosaic)
  #define dspread (_graphite_analyser17->_parameters.dspread)
  #define Q (_graphite_analyser17->_parameters.Q)
  #define DM (_graphite_analyser17->_parameters.DM)
  #define pThreshold (_graphite_analyser17->_parameters.pThreshold)
  #define Rup (_graphite_analyser17->_parameters.Rup)
  #define Rdown (_graphite_analyser17->_parameters.Rdown)
  #define debug (_graphite_analyser17->_parameters.debug)
  #define mos_rms (_graphite_analyser17->_parameters.mos_rms)
  #define d_rms (_graphite_analyser17->_parameters.d_rms)
  #define mono_Q (_graphite_analyser17->_parameters.mono_Q)
  #define FN (_graphite_analyser17->_parameters.FN)
  #define FM (_graphite_analyser17->_parameters.FM)
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_graphite_analyser17->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _graphite_analyser17->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_graphite_analyser17->_position_relative, _graphite_analyser17->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      if (( instrument->_parameters._LAMBDA > 0 )) // conditional WHEN execution
      class_Monochromator_pol_trace(_graphite_analyser17, _particle); /* contains EXTEND code */
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
      // GROUP graphiteAnalyzer: from graphite_analyser0 [28] to graphite_analyser39 [67]
      if (SCATTERED) _particle->_index = 67; // when SCATTERED in GROUP: reach exit of GROUP after graphite_analyser39
      else ABSORBED=0; // not SCATTERED within GROUP: always tries next
      _particle->_index++;
    } /* end component graphite_analyser17 [45] */
    if (!ABSORBED && _particle->_index == 46) {
      /* component graphite_analyser18=Monochromator_pol() [46] */
  #define zwidth (_graphite_analyser18->_parameters.zwidth)
  #define yheight (_graphite_analyser18->_parameters.yheight)
  #define mosaic (_graphite_analyser18->_parameters.mosaic)
  #define dspread (_graphite_analyser18->_parameters.dspread)
  #define Q (_graphite_analyser18->_parameters.Q)
  #define DM (_graphite_analyser18->_parameters.DM)
  #define pThreshold (_graphite_analyser18->_parameters.pThreshold)
  #define Rup (_graphite_analyser18->_parameters.Rup)
  #define Rdown (_graphite_analyser18->_parameters.Rdown)
  #define debug (_graphite_analyser18->_parameters.debug)
  #define mos_rms (_graphite_analyser18->_parameters.mos_rms)
  #define d_rms (_graphite_analyser18->_parameters.d_rms)
  #define mono_Q (_graphite_analyser18->_parameters.mono_Q)
  #define FN (_graphite_analyser18->_parameters.FN)
  #define FM (_graphite_analyser18->_parameters.FM)
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_graphite_analyser18->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _graphite_analyser18->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_graphite_analyser18->_position_relative, _graphite_analyser18->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      if (( instrument->_parameters._LAMBDA > 0 )) // conditional WHEN execution
      class_Monochromator_pol_trace(_graphite_analyser18, _particle); /* contains EXTEND code */
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
      // GROUP graphiteAnalyzer: from graphite_analyser0 [28] to graphite_analyser39 [67]
      if (SCATTERED) _particle->_index = 67; // when SCATTERED in GROUP: reach exit of GROUP after graphite_analyser39
      else ABSORBED=0; // not SCATTERED within GROUP: always tries next
      _particle->_index++;
    } /* end component graphite_analyser18 [46] */
    if (!ABSORBED && _particle->_index == 47) {
      /* component graphite_analyser19=Monochromator_pol() [47] */
  #define zwidth (_graphite_analyser19->_parameters.zwidth)
  #define yheight (_graphite_analyser19->_parameters.yheight)
  #define mosaic (_graphite_analyser19->_parameters.mosaic)
  #define dspread (_graphite_analyser19->_parameters.dspread)
  #define Q (_graphite_analyser19->_parameters.Q)
  #define DM (_graphite_analyser19->_parameters.DM)
  #define pThreshold (_graphite_analyser19->_parameters.pThreshold)
  #define Rup (_graphite_analyser19->_parameters.Rup)
  #define Rdown (_graphite_analyser19->_parameters.Rdown)
  #define debug (_graphite_analyser19->_parameters.debug)
  #define mos_rms (_graphite_analyser19->_parameters.mos_rms)
  #define d_rms (_graphite_analyser19->_parameters.d_rms)
  #define mono_Q (_graphite_analyser19->_parameters.mono_Q)
  #define FN (_graphite_analyser19->_parameters.FN)
  #define FM (_graphite_analyser19->_parameters.FM)
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_graphite_analyser19->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _graphite_analyser19->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_graphite_analyser19->_position_relative, _graphite_analyser19->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      if (( instrument->_parameters._LAMBDA > 0 )) // conditional WHEN execution
      class_Monochromator_pol_trace(_graphite_analyser19, _particle); /* contains EXTEND code */
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
      // GROUP graphiteAnalyzer: from graphite_analyser0 [28] to graphite_analyser39 [67]
      if (SCATTERED) _particle->_index = 67; // when SCATTERED in GROUP: reach exit of GROUP after graphite_analyser39
      else ABSORBED=0; // not SCATTERED within GROUP: always tries next
      _particle->_index++;
    } /* end component graphite_analyser19 [47] */
    if (!ABSORBED && _particle->_index == 48) {
      /* component graphite_analyser20=Monochromator_pol() [48] */
  #define zwidth (_graphite_analyser20->_parameters.zwidth)
  #define yheight (_graphite_analyser20->_parameters.yheight)
  #define mosaic (_graphite_analyser20->_parameters.mosaic)
  #define dspread (_graphite_analyser20->_parameters.dspread)
  #define Q (_graphite_analyser20->_parameters.Q)
  #define DM (_graphite_analyser20->_parameters.DM)
  #define pThreshold (_graphite_analyser20->_parameters.pThreshold)
  #define Rup (_graphite_analyser20->_parameters.Rup)
  #define Rdown (_graphite_analyser20->_parameters.Rdown)
  #define debug (_graphite_analyser20->_parameters.debug)
  #define mos_rms (_graphite_analyser20->_parameters.mos_rms)
  #define d_rms (_graphite_analyser20->_parameters.d_rms)
  #define mono_Q (_graphite_analyser20->_parameters.mono_Q)
  #define FN (_graphite_analyser20->_parameters.FN)
  #define FM (_graphite_analyser20->_parameters.FM)
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_graphite_analyser20->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _graphite_analyser20->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_graphite_analyser20->_position_relative, _graphite_analyser20->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      if (( instrument->_parameters._LAMBDA > 0 )) // conditional WHEN execution
      class_Monochromator_pol_trace(_graphite_analyser20, _particle); /* contains EXTEND code */
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
      // GROUP graphiteAnalyzer: from graphite_analyser0 [28] to graphite_analyser39 [67]
      if (SCATTERED) _particle->_index = 67; // when SCATTERED in GROUP: reach exit of GROUP after graphite_analyser39
      else ABSORBED=0; // not SCATTERED within GROUP: always tries next
      _particle->_index++;
    } /* end component graphite_analyser20 [48] */
    if (!ABSORBED && _particle->_index == 49) {
      /* component graphite_analyser21=Monochromator_pol() [49] */
  #define zwidth (_graphite_analyser21->_parameters.zwidth)
  #define yheight (_graphite_analyser21->_parameters.yheight)
  #define mosaic (_graphite_analyser21->_parameters.mosaic)
  #define dspread (_graphite_analyser21->_parameters.dspread)
  #define Q (_graphite_analyser21->_parameters.Q)
  #define DM (_graphite_analyser21->_parameters.DM)
  #define pThreshold (_graphite_analyser21->_parameters.pThreshold)
  #define Rup (_graphite_analyser21->_parameters.Rup)
  #define Rdown (_graphite_analyser21->_parameters.Rdown)
  #define debug (_graphite_analyser21->_parameters.debug)
  #define mos_rms (_graphite_analyser21->_parameters.mos_rms)
  #define d_rms (_graphite_analyser21->_parameters.d_rms)
  #define mono_Q (_graphite_analyser21->_parameters.mono_Q)
  #define FN (_graphite_analyser21->_parameters.FN)
  #define FM (_graphite_analyser21->_parameters.FM)
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_graphite_analyser21->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _graphite_analyser21->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_graphite_analyser21->_position_relative, _graphite_analyser21->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      if (( instrument->_parameters._LAMBDA > 0 )) // conditional WHEN execution
      class_Monochromator_pol_trace(_graphite_analyser21, _particle); /* contains EXTEND code */
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
      // GROUP graphiteAnalyzer: from graphite_analyser0 [28] to graphite_analyser39 [67]
      if (SCATTERED) _particle->_index = 67; // when SCATTERED in GROUP: reach exit of GROUP after graphite_analyser39
      else ABSORBED=0; // not SCATTERED within GROUP: always tries next
      _particle->_index++;
    } /* end component graphite_analyser21 [49] */
    if (!ABSORBED && _particle->_index == 50) {
      /* component graphite_analyser22=Monochromator_pol() [50] */
  #define zwidth (_graphite_analyser22->_parameters.zwidth)
  #define yheight (_graphite_analyser22->_parameters.yheight)
  #define mosaic (_graphite_analyser22->_parameters.mosaic)
  #define dspread (_graphite_analyser22->_parameters.dspread)
  #define Q (_graphite_analyser22->_parameters.Q)
  #define DM (_graphite_analyser22->_parameters.DM)
  #define pThreshold (_graphite_analyser22->_parameters.pThreshold)
  #define Rup (_graphite_analyser22->_parameters.Rup)
  #define Rdown (_graphite_analyser22->_parameters.Rdown)
  #define debug (_graphite_analyser22->_parameters.debug)
  #define mos_rms (_graphite_analyser22->_parameters.mos_rms)
  #define d_rms (_graphite_analyser22->_parameters.d_rms)
  #define mono_Q (_graphite_analyser22->_parameters.mono_Q)
  #define FN (_graphite_analyser22->_parameters.FN)
  #define FM (_graphite_analyser22->_parameters.FM)
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_graphite_analyser22->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _graphite_analyser22->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_graphite_analyser22->_position_relative, _graphite_analyser22->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      if (( instrument->_parameters._LAMBDA > 0 )) // conditional WHEN execution
      class_Monochromator_pol_trace(_graphite_analyser22, _particle); /* contains EXTEND code */
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
      // GROUP graphiteAnalyzer: from graphite_analyser0 [28] to graphite_analyser39 [67]
      if (SCATTERED) _particle->_index = 67; // when SCATTERED in GROUP: reach exit of GROUP after graphite_analyser39
      else ABSORBED=0; // not SCATTERED within GROUP: always tries next
      _particle->_index++;
    } /* end component graphite_analyser22 [50] */
    if (!ABSORBED && _particle->_index == 51) {
      /* component graphite_analyser23=Monochromator_pol() [51] */
  #define zwidth (_graphite_analyser23->_parameters.zwidth)
  #define yheight (_graphite_analyser23->_parameters.yheight)
  #define mosaic (_graphite_analyser23->_parameters.mosaic)
  #define dspread (_graphite_analyser23->_parameters.dspread)
  #define Q (_graphite_analyser23->_parameters.Q)
  #define DM (_graphite_analyser23->_parameters.DM)
  #define pThreshold (_graphite_analyser23->_parameters.pThreshold)
  #define Rup (_graphite_analyser23->_parameters.Rup)
  #define Rdown (_graphite_analyser23->_parameters.Rdown)
  #define debug (_graphite_analyser23->_parameters.debug)
  #define mos_rms (_graphite_analyser23->_parameters.mos_rms)
  #define d_rms (_graphite_analyser23->_parameters.d_rms)
  #define mono_Q (_graphite_analyser23->_parameters.mono_Q)
  #define FN (_graphite_analyser23->_parameters.FN)
  #define FM (_graphite_analyser23->_parameters.FM)
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_graphite_analyser23->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _graphite_analyser23->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_graphite_analyser23->_position_relative, _graphite_analyser23->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      if (( instrument->_parameters._LAMBDA > 0 )) // conditional WHEN execution
      class_Monochromator_pol_trace(_graphite_analyser23, _particle); /* contains EXTEND code */
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
      // GROUP graphiteAnalyzer: from graphite_analyser0 [28] to graphite_analyser39 [67]
      if (SCATTERED) _particle->_index = 67; // when SCATTERED in GROUP: reach exit of GROUP after graphite_analyser39
      else ABSORBED=0; // not SCATTERED within GROUP: always tries next
      _particle->_index++;
    } /* end component graphite_analyser23 [51] */
    if (!ABSORBED && _particle->_index == 52) {
      /* component graphite_analyser24=Monochromator_pol() [52] */
  #define zwidth (_graphite_analyser24->_parameters.zwidth)
  #define yheight (_graphite_analyser24->_parameters.yheight)
  #define mosaic (_graphite_analyser24->_parameters.mosaic)
  #define dspread (_graphite_analyser24->_parameters.dspread)
  #define Q (_graphite_analyser24->_parameters.Q)
  #define DM (_graphite_analyser24->_parameters.DM)
  #define pThreshold (_graphite_analyser24->_parameters.pThreshold)
  #define Rup (_graphite_analyser24->_parameters.Rup)
  #define Rdown (_graphite_analyser24->_parameters.Rdown)
  #define debug (_graphite_analyser24->_parameters.debug)
  #define mos_rms (_graphite_analyser24->_parameters.mos_rms)
  #define d_rms (_graphite_analyser24->_parameters.d_rms)
  #define mono_Q (_graphite_analyser24->_parameters.mono_Q)
  #define FN (_graphite_analyser24->_parameters.FN)
  #define FM (_graphite_analyser24->_parameters.FM)
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_graphite_analyser24->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _graphite_analyser24->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_graphite_analyser24->_position_relative, _graphite_analyser24->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      if (( instrument->_parameters._LAMBDA > 0 )) // conditional WHEN execution
      class_Monochromator_pol_trace(_graphite_analyser24, _particle); /* contains EXTEND code */
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
      // GROUP graphiteAnalyzer: from graphite_analyser0 [28] to graphite_analyser39 [67]
      if (SCATTERED) _particle->_index = 67; // when SCATTERED in GROUP: reach exit of GROUP after graphite_analyser39
      else ABSORBED=0; // not SCATTERED within GROUP: always tries next
      _particle->_index++;
    } /* end component graphite_analyser24 [52] */
    if (!ABSORBED && _particle->_index == 53) {
      /* component graphite_analyser25=Monochromator_pol() [53] */
  #define zwidth (_graphite_analyser25->_parameters.zwidth)
  #define yheight (_graphite_analyser25->_parameters.yheight)
  #define mosaic (_graphite_analyser25->_parameters.mosaic)
  #define dspread (_graphite_analyser25->_parameters.dspread)
  #define Q (_graphite_analyser25->_parameters.Q)
  #define DM (_graphite_analyser25->_parameters.DM)
  #define pThreshold (_graphite_analyser25->_parameters.pThreshold)
  #define Rup (_graphite_analyser25->_parameters.Rup)
  #define Rdown (_graphite_analyser25->_parameters.Rdown)
  #define debug (_graphite_analyser25->_parameters.debug)
  #define mos_rms (_graphite_analyser25->_parameters.mos_rms)
  #define d_rms (_graphite_analyser25->_parameters.d_rms)
  #define mono_Q (_graphite_analyser25->_parameters.mono_Q)
  #define FN (_graphite_analyser25->_parameters.FN)
  #define FM (_graphite_analyser25->_parameters.FM)
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_graphite_analyser25->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _graphite_analyser25->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_graphite_analyser25->_position_relative, _graphite_analyser25->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      if (( instrument->_parameters._LAMBDA > 0 )) // conditional WHEN execution
      class_Monochromator_pol_trace(_graphite_analyser25, _particle); /* contains EXTEND code */
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
      // GROUP graphiteAnalyzer: from graphite_analyser0 [28] to graphite_analyser39 [67]
      if (SCATTERED) _particle->_index = 67; // when SCATTERED in GROUP: reach exit of GROUP after graphite_analyser39
      else ABSORBED=0; // not SCATTERED within GROUP: always tries next
      _particle->_index++;
    } /* end component graphite_analyser25 [53] */
    if (!ABSORBED && _particle->_index == 54) {
      /* component graphite_analyser26=Monochromator_pol() [54] */
  #define zwidth (_graphite_analyser26->_parameters.zwidth)
  #define yheight (_graphite_analyser26->_parameters.yheight)
  #define mosaic (_graphite_analyser26->_parameters.mosaic)
  #define dspread (_graphite_analyser26->_parameters.dspread)
  #define Q (_graphite_analyser26->_parameters.Q)
  #define DM (_graphite_analyser26->_parameters.DM)
  #define pThreshold (_graphite_analyser26->_parameters.pThreshold)
  #define Rup (_graphite_analyser26->_parameters.Rup)
  #define Rdown (_graphite_analyser26->_parameters.Rdown)
  #define debug (_graphite_analyser26->_parameters.debug)
  #define mos_rms (_graphite_analyser26->_parameters.mos_rms)
  #define d_rms (_graphite_analyser26->_parameters.d_rms)
  #define mono_Q (_graphite_analyser26->_parameters.mono_Q)
  #define FN (_graphite_analyser26->_parameters.FN)
  #define FM (_graphite_analyser26->_parameters.FM)
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_graphite_analyser26->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _graphite_analyser26->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_graphite_analyser26->_position_relative, _graphite_analyser26->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      if (( instrument->_parameters._LAMBDA > 0 )) // conditional WHEN execution
      class_Monochromator_pol_trace(_graphite_analyser26, _particle); /* contains EXTEND code */
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
      // GROUP graphiteAnalyzer: from graphite_analyser0 [28] to graphite_analyser39 [67]
      if (SCATTERED) _particle->_index = 67; // when SCATTERED in GROUP: reach exit of GROUP after graphite_analyser39
      else ABSORBED=0; // not SCATTERED within GROUP: always tries next
      _particle->_index++;
    } /* end component graphite_analyser26 [54] */
    if (!ABSORBED && _particle->_index == 55) {
      /* component graphite_analyser27=Monochromator_pol() [55] */
  #define zwidth (_graphite_analyser27->_parameters.zwidth)
  #define yheight (_graphite_analyser27->_parameters.yheight)
  #define mosaic (_graphite_analyser27->_parameters.mosaic)
  #define dspread (_graphite_analyser27->_parameters.dspread)
  #define Q (_graphite_analyser27->_parameters.Q)
  #define DM (_graphite_analyser27->_parameters.DM)
  #define pThreshold (_graphite_analyser27->_parameters.pThreshold)
  #define Rup (_graphite_analyser27->_parameters.Rup)
  #define Rdown (_graphite_analyser27->_parameters.Rdown)
  #define debug (_graphite_analyser27->_parameters.debug)
  #define mos_rms (_graphite_analyser27->_parameters.mos_rms)
  #define d_rms (_graphite_analyser27->_parameters.d_rms)
  #define mono_Q (_graphite_analyser27->_parameters.mono_Q)
  #define FN (_graphite_analyser27->_parameters.FN)
  #define FM (_graphite_analyser27->_parameters.FM)
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_graphite_analyser27->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _graphite_analyser27->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_graphite_analyser27->_position_relative, _graphite_analyser27->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      if (( instrument->_parameters._LAMBDA > 0 )) // conditional WHEN execution
      class_Monochromator_pol_trace(_graphite_analyser27, _particle); /* contains EXTEND code */
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
      // GROUP graphiteAnalyzer: from graphite_analyser0 [28] to graphite_analyser39 [67]
      if (SCATTERED) _particle->_index = 67; // when SCATTERED in GROUP: reach exit of GROUP after graphite_analyser39
      else ABSORBED=0; // not SCATTERED within GROUP: always tries next
      _particle->_index++;
    } /* end component graphite_analyser27 [55] */
    if (!ABSORBED && _particle->_index == 56) {
      /* component graphite_analyser28=Monochromator_pol() [56] */
  #define zwidth (_graphite_analyser28->_parameters.zwidth)
  #define yheight (_graphite_analyser28->_parameters.yheight)
  #define mosaic (_graphite_analyser28->_parameters.mosaic)
  #define dspread (_graphite_analyser28->_parameters.dspread)
  #define Q (_graphite_analyser28->_parameters.Q)
  #define DM (_graphite_analyser28->_parameters.DM)
  #define pThreshold (_graphite_analyser28->_parameters.pThreshold)
  #define Rup (_graphite_analyser28->_parameters.Rup)
  #define Rdown (_graphite_analyser28->_parameters.Rdown)
  #define debug (_graphite_analyser28->_parameters.debug)
  #define mos_rms (_graphite_analyser28->_parameters.mos_rms)
  #define d_rms (_graphite_analyser28->_parameters.d_rms)
  #define mono_Q (_graphite_analyser28->_parameters.mono_Q)
  #define FN (_graphite_analyser28->_parameters.FN)
  #define FM (_graphite_analyser28->_parameters.FM)
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_graphite_analyser28->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _graphite_analyser28->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_graphite_analyser28->_position_relative, _graphite_analyser28->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      if (( instrument->_parameters._LAMBDA > 0 )) // conditional WHEN execution
      class_Monochromator_pol_trace(_graphite_analyser28, _particle); /* contains EXTEND code */
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
      // GROUP graphiteAnalyzer: from graphite_analyser0 [28] to graphite_analyser39 [67]
      if (SCATTERED) _particle->_index = 67; // when SCATTERED in GROUP: reach exit of GROUP after graphite_analyser39
      else ABSORBED=0; // not SCATTERED within GROUP: always tries next
      _particle->_index++;
    } /* end component graphite_analyser28 [56] */
    if (!ABSORBED && _particle->_index == 57) {
      /* component graphite_analyser29=Monochromator_pol() [57] */
  #define zwidth (_graphite_analyser29->_parameters.zwidth)
  #define yheight (_graphite_analyser29->_parameters.yheight)
  #define mosaic (_graphite_analyser29->_parameters.mosaic)
  #define dspread (_graphite_analyser29->_parameters.dspread)
  #define Q (_graphite_analyser29->_parameters.Q)
  #define DM (_graphite_analyser29->_parameters.DM)
  #define pThreshold (_graphite_analyser29->_parameters.pThreshold)
  #define Rup (_graphite_analyser29->_parameters.Rup)
  #define Rdown (_graphite_analyser29->_parameters.Rdown)
  #define debug (_graphite_analyser29->_parameters.debug)
  #define mos_rms (_graphite_analyser29->_parameters.mos_rms)
  #define d_rms (_graphite_analyser29->_parameters.d_rms)
  #define mono_Q (_graphite_analyser29->_parameters.mono_Q)
  #define FN (_graphite_analyser29->_parameters.FN)
  #define FM (_graphite_analyser29->_parameters.FM)
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_graphite_analyser29->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _graphite_analyser29->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_graphite_analyser29->_position_relative, _graphite_analyser29->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      if (( instrument->_parameters._LAMBDA > 0 )) // conditional WHEN execution
      class_Monochromator_pol_trace(_graphite_analyser29, _particle); /* contains EXTEND code */
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
      // GROUP graphiteAnalyzer: from graphite_analyser0 [28] to graphite_analyser39 [67]
      if (SCATTERED) _particle->_index = 67; // when SCATTERED in GROUP: reach exit of GROUP after graphite_analyser39
      else ABSORBED=0; // not SCATTERED within GROUP: always tries next
      _particle->_index++;
    } /* end component graphite_analyser29 [57] */
    if (!ABSORBED && _particle->_index == 58) {
      /* component graphite_analyser30=Monochromator_pol() [58] */
  #define zwidth (_graphite_analyser30->_parameters.zwidth)
  #define yheight (_graphite_analyser30->_parameters.yheight)
  #define mosaic (_graphite_analyser30->_parameters.mosaic)
  #define dspread (_graphite_analyser30->_parameters.dspread)
  #define Q (_graphite_analyser30->_parameters.Q)
  #define DM (_graphite_analyser30->_parameters.DM)
  #define pThreshold (_graphite_analyser30->_parameters.pThreshold)
  #define Rup (_graphite_analyser30->_parameters.Rup)
  #define Rdown (_graphite_analyser30->_parameters.Rdown)
  #define debug (_graphite_analyser30->_parameters.debug)
  #define mos_rms (_graphite_analyser30->_parameters.mos_rms)
  #define d_rms (_graphite_analyser30->_parameters.d_rms)
  #define mono_Q (_graphite_analyser30->_parameters.mono_Q)
  #define FN (_graphite_analyser30->_parameters.FN)
  #define FM (_graphite_analyser30->_parameters.FM)
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_graphite_analyser30->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _graphite_analyser30->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_graphite_analyser30->_position_relative, _graphite_analyser30->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      if (( instrument->_parameters._LAMBDA > 0 )) // conditional WHEN execution
      class_Monochromator_pol_trace(_graphite_analyser30, _particle); /* contains EXTEND code */
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
      // GROUP graphiteAnalyzer: from graphite_analyser0 [28] to graphite_analyser39 [67]
      if (SCATTERED) _particle->_index = 67; // when SCATTERED in GROUP: reach exit of GROUP after graphite_analyser39
      else ABSORBED=0; // not SCATTERED within GROUP: always tries next
      _particle->_index++;
    } /* end component graphite_analyser30 [58] */
    if (!ABSORBED && _particle->_index == 59) {
      /* component graphite_analyser31=Monochromator_pol() [59] */
  #define zwidth (_graphite_analyser31->_parameters.zwidth)
  #define yheight (_graphite_analyser31->_parameters.yheight)
  #define mosaic (_graphite_analyser31->_parameters.mosaic)
  #define dspread (_graphite_analyser31->_parameters.dspread)
  #define Q (_graphite_analyser31->_parameters.Q)
  #define DM (_graphite_analyser31->_parameters.DM)
  #define pThreshold (_graphite_analyser31->_parameters.pThreshold)
  #define Rup (_graphite_analyser31->_parameters.Rup)
  #define Rdown (_graphite_analyser31->_parameters.Rdown)
  #define debug (_graphite_analyser31->_parameters.debug)
  #define mos_rms (_graphite_analyser31->_parameters.mos_rms)
  #define d_rms (_graphite_analyser31->_parameters.d_rms)
  #define mono_Q (_graphite_analyser31->_parameters.mono_Q)
  #define FN (_graphite_analyser31->_parameters.FN)
  #define FM (_graphite_analyser31->_parameters.FM)
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_graphite_analyser31->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _graphite_analyser31->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_graphite_analyser31->_position_relative, _graphite_analyser31->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      if (( instrument->_parameters._LAMBDA > 0 )) // conditional WHEN execution
      class_Monochromator_pol_trace(_graphite_analyser31, _particle); /* contains EXTEND code */
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
      // GROUP graphiteAnalyzer: from graphite_analyser0 [28] to graphite_analyser39 [67]
      if (SCATTERED) _particle->_index = 67; // when SCATTERED in GROUP: reach exit of GROUP after graphite_analyser39
      else ABSORBED=0; // not SCATTERED within GROUP: always tries next
      _particle->_index++;
    } /* end component graphite_analyser31 [59] */
    if (!ABSORBED && _particle->_index == 60) {
      /* component graphite_analyser32=Monochromator_pol() [60] */
  #define zwidth (_graphite_analyser32->_parameters.zwidth)
  #define yheight (_graphite_analyser32->_parameters.yheight)
  #define mosaic (_graphite_analyser32->_parameters.mosaic)
  #define dspread (_graphite_analyser32->_parameters.dspread)
  #define Q (_graphite_analyser32->_parameters.Q)
  #define DM (_graphite_analyser32->_parameters.DM)
  #define pThreshold (_graphite_analyser32->_parameters.pThreshold)
  #define Rup (_graphite_analyser32->_parameters.Rup)
  #define Rdown (_graphite_analyser32->_parameters.Rdown)
  #define debug (_graphite_analyser32->_parameters.debug)
  #define mos_rms (_graphite_analyser32->_parameters.mos_rms)
  #define d_rms (_graphite_analyser32->_parameters.d_rms)
  #define mono_Q (_graphite_analyser32->_parameters.mono_Q)
  #define FN (_graphite_analyser32->_parameters.FN)
  #define FM (_graphite_analyser32->_parameters.FM)
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_graphite_analyser32->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _graphite_analyser32->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_graphite_analyser32->_position_relative, _graphite_analyser32->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      if (( instrument->_parameters._LAMBDA > 0 )) // conditional WHEN execution
      class_Monochromator_pol_trace(_graphite_analyser32, _particle); /* contains EXTEND code */
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
      // GROUP graphiteAnalyzer: from graphite_analyser0 [28] to graphite_analyser39 [67]
      if (SCATTERED) _particle->_index = 67; // when SCATTERED in GROUP: reach exit of GROUP after graphite_analyser39
      else ABSORBED=0; // not SCATTERED within GROUP: always tries next
      _particle->_index++;
    } /* end component graphite_analyser32 [60] */
    if (!ABSORBED && _particle->_index == 61) {
      /* component graphite_analyser33=Monochromator_pol() [61] */
  #define zwidth (_graphite_analyser33->_parameters.zwidth)
  #define yheight (_graphite_analyser33->_parameters.yheight)
  #define mosaic (_graphite_analyser33->_parameters.mosaic)
  #define dspread (_graphite_analyser33->_parameters.dspread)
  #define Q (_graphite_analyser33->_parameters.Q)
  #define DM (_graphite_analyser33->_parameters.DM)
  #define pThreshold (_graphite_analyser33->_parameters.pThreshold)
  #define Rup (_graphite_analyser33->_parameters.Rup)
  #define Rdown (_graphite_analyser33->_parameters.Rdown)
  #define debug (_graphite_analyser33->_parameters.debug)
  #define mos_rms (_graphite_analyser33->_parameters.mos_rms)
  #define d_rms (_graphite_analyser33->_parameters.d_rms)
  #define mono_Q (_graphite_analyser33->_parameters.mono_Q)
  #define FN (_graphite_analyser33->_parameters.FN)
  #define FM (_graphite_analyser33->_parameters.FM)
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_graphite_analyser33->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _graphite_analyser33->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_graphite_analyser33->_position_relative, _graphite_analyser33->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      if (( instrument->_parameters._LAMBDA > 0 )) // conditional WHEN execution
      class_Monochromator_pol_trace(_graphite_analyser33, _particle); /* contains EXTEND code */
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
      // GROUP graphiteAnalyzer: from graphite_analyser0 [28] to graphite_analyser39 [67]
      if (SCATTERED) _particle->_index = 67; // when SCATTERED in GROUP: reach exit of GROUP after graphite_analyser39
      else ABSORBED=0; // not SCATTERED within GROUP: always tries next
      _particle->_index++;
    } /* end component graphite_analyser33 [61] */
    if (!ABSORBED && _particle->_index == 62) {
      /* component graphite_analyser34=Monochromator_pol() [62] */
  #define zwidth (_graphite_analyser34->_parameters.zwidth)
  #define yheight (_graphite_analyser34->_parameters.yheight)
  #define mosaic (_graphite_analyser34->_parameters.mosaic)
  #define dspread (_graphite_analyser34->_parameters.dspread)
  #define Q (_graphite_analyser34->_parameters.Q)
  #define DM (_graphite_analyser34->_parameters.DM)
  #define pThreshold (_graphite_analyser34->_parameters.pThreshold)
  #define Rup (_graphite_analyser34->_parameters.Rup)
  #define Rdown (_graphite_analyser34->_parameters.Rdown)
  #define debug (_graphite_analyser34->_parameters.debug)
  #define mos_rms (_graphite_analyser34->_parameters.mos_rms)
  #define d_rms (_graphite_analyser34->_parameters.d_rms)
  #define mono_Q (_graphite_analyser34->_parameters.mono_Q)
  #define FN (_graphite_analyser34->_parameters.FN)
  #define FM (_graphite_analyser34->_parameters.FM)
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_graphite_analyser34->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _graphite_analyser34->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_graphite_analyser34->_position_relative, _graphite_analyser34->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      if (( instrument->_parameters._LAMBDA > 0 )) // conditional WHEN execution
      class_Monochromator_pol_trace(_graphite_analyser34, _particle); /* contains EXTEND code */
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
      // GROUP graphiteAnalyzer: from graphite_analyser0 [28] to graphite_analyser39 [67]
      if (SCATTERED) _particle->_index = 67; // when SCATTERED in GROUP: reach exit of GROUP after graphite_analyser39
      else ABSORBED=0; // not SCATTERED within GROUP: always tries next
      _particle->_index++;
    } /* end component graphite_analyser34 [62] */
    if (!ABSORBED && _particle->_index == 63) {
      /* component graphite_analyser35=Monochromator_pol() [63] */
  #define zwidth (_graphite_analyser35->_parameters.zwidth)
  #define yheight (_graphite_analyser35->_parameters.yheight)
  #define mosaic (_graphite_analyser35->_parameters.mosaic)
  #define dspread (_graphite_analyser35->_parameters.dspread)
  #define Q (_graphite_analyser35->_parameters.Q)
  #define DM (_graphite_analyser35->_parameters.DM)
  #define pThreshold (_graphite_analyser35->_parameters.pThreshold)
  #define Rup (_graphite_analyser35->_parameters.Rup)
  #define Rdown (_graphite_analyser35->_parameters.Rdown)
  #define debug (_graphite_analyser35->_parameters.debug)
  #define mos_rms (_graphite_analyser35->_parameters.mos_rms)
  #define d_rms (_graphite_analyser35->_parameters.d_rms)
  #define mono_Q (_graphite_analyser35->_parameters.mono_Q)
  #define FN (_graphite_analyser35->_parameters.FN)
  #define FM (_graphite_analyser35->_parameters.FM)
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_graphite_analyser35->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _graphite_analyser35->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_graphite_analyser35->_position_relative, _graphite_analyser35->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      if (( instrument->_parameters._LAMBDA > 0 )) // conditional WHEN execution
      class_Monochromator_pol_trace(_graphite_analyser35, _particle); /* contains EXTEND code */
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
      // GROUP graphiteAnalyzer: from graphite_analyser0 [28] to graphite_analyser39 [67]
      if (SCATTERED) _particle->_index = 67; // when SCATTERED in GROUP: reach exit of GROUP after graphite_analyser39
      else ABSORBED=0; // not SCATTERED within GROUP: always tries next
      _particle->_index++;
    } /* end component graphite_analyser35 [63] */
    if (!ABSORBED && _particle->_index == 64) {
      /* component graphite_analyser36=Monochromator_pol() [64] */
  #define zwidth (_graphite_analyser36->_parameters.zwidth)
  #define yheight (_graphite_analyser36->_parameters.yheight)
  #define mosaic (_graphite_analyser36->_parameters.mosaic)
  #define dspread (_graphite_analyser36->_parameters.dspread)
  #define Q (_graphite_analyser36->_parameters.Q)
  #define DM (_graphite_analyser36->_parameters.DM)
  #define pThreshold (_graphite_analyser36->_parameters.pThreshold)
  #define Rup (_graphite_analyser36->_parameters.Rup)
  #define Rdown (_graphite_analyser36->_parameters.Rdown)
  #define debug (_graphite_analyser36->_parameters.debug)
  #define mos_rms (_graphite_analyser36->_parameters.mos_rms)
  #define d_rms (_graphite_analyser36->_parameters.d_rms)
  #define mono_Q (_graphite_analyser36->_parameters.mono_Q)
  #define FN (_graphite_analyser36->_parameters.FN)
  #define FM (_graphite_analyser36->_parameters.FM)
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_graphite_analyser36->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _graphite_analyser36->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_graphite_analyser36->_position_relative, _graphite_analyser36->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      if (( instrument->_parameters._LAMBDA > 0 )) // conditional WHEN execution
      class_Monochromator_pol_trace(_graphite_analyser36, _particle); /* contains EXTEND code */
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
      // GROUP graphiteAnalyzer: from graphite_analyser0 [28] to graphite_analyser39 [67]
      if (SCATTERED) _particle->_index = 67; // when SCATTERED in GROUP: reach exit of GROUP after graphite_analyser39
      else ABSORBED=0; // not SCATTERED within GROUP: always tries next
      _particle->_index++;
    } /* end component graphite_analyser36 [64] */
    if (!ABSORBED && _particle->_index == 65) {
      /* component graphite_analyser37=Monochromator_pol() [65] */
  #define zwidth (_graphite_analyser37->_parameters.zwidth)
  #define yheight (_graphite_analyser37->_parameters.yheight)
  #define mosaic (_graphite_analyser37->_parameters.mosaic)
  #define dspread (_graphite_analyser37->_parameters.dspread)
  #define Q (_graphite_analyser37->_parameters.Q)
  #define DM (_graphite_analyser37->_parameters.DM)
  #define pThreshold (_graphite_analyser37->_parameters.pThreshold)
  #define Rup (_graphite_analyser37->_parameters.Rup)
  #define Rdown (_graphite_analyser37->_parameters.Rdown)
  #define debug (_graphite_analyser37->_parameters.debug)
  #define mos_rms (_graphite_analyser37->_parameters.mos_rms)
  #define d_rms (_graphite_analyser37->_parameters.d_rms)
  #define mono_Q (_graphite_analyser37->_parameters.mono_Q)
  #define FN (_graphite_analyser37->_parameters.FN)
  #define FM (_graphite_analyser37->_parameters.FM)
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_graphite_analyser37->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _graphite_analyser37->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_graphite_analyser37->_position_relative, _graphite_analyser37->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      if (( instrument->_parameters._LAMBDA > 0 )) // conditional WHEN execution
      class_Monochromator_pol_trace(_graphite_analyser37, _particle); /* contains EXTEND code */
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
      // GROUP graphiteAnalyzer: from graphite_analyser0 [28] to graphite_analyser39 [67]
      if (SCATTERED) _particle->_index = 67; // when SCATTERED in GROUP: reach exit of GROUP after graphite_analyser39
      else ABSORBED=0; // not SCATTERED within GROUP: always tries next
      _particle->_index++;
    } /* end component graphite_analyser37 [65] */
    if (!ABSORBED && _particle->_index == 66) {
      /* component graphite_analyser38=Monochromator_pol() [66] */
  #define zwidth (_graphite_analyser38->_parameters.zwidth)
  #define yheight (_graphite_analyser38->_parameters.yheight)
  #define mosaic (_graphite_analyser38->_parameters.mosaic)
  #define dspread (_graphite_analyser38->_parameters.dspread)
  #define Q (_graphite_analyser38->_parameters.Q)
  #define DM (_graphite_analyser38->_parameters.DM)
  #define pThreshold (_graphite_analyser38->_parameters.pThreshold)
  #define Rup (_graphite_analyser38->_parameters.Rup)
  #define Rdown (_graphite_analyser38->_parameters.Rdown)
  #define debug (_graphite_analyser38->_parameters.debug)
  #define mos_rms (_graphite_analyser38->_parameters.mos_rms)
  #define d_rms (_graphite_analyser38->_parameters.d_rms)
  #define mono_Q (_graphite_analyser38->_parameters.mono_Q)
  #define FN (_graphite_analyser38->_parameters.FN)
  #define FM (_graphite_analyser38->_parameters.FM)
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_graphite_analyser38->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _graphite_analyser38->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_graphite_analyser38->_position_relative, _graphite_analyser38->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      if (( instrument->_parameters._LAMBDA > 0 )) // conditional WHEN execution
      class_Monochromator_pol_trace(_graphite_analyser38, _particle); /* contains EXTEND code */
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
      // GROUP graphiteAnalyzer: from graphite_analyser0 [28] to graphite_analyser39 [67]
      if (SCATTERED) _particle->_index = 67; // when SCATTERED in GROUP: reach exit of GROUP after graphite_analyser39
      else ABSORBED=0; // not SCATTERED within GROUP: always tries next
      _particle->_index++;
    } /* end component graphite_analyser38 [66] */
    if (!ABSORBED && _particle->_index == 67) {
      /* component graphite_analyser39=Monochromator_pol() [67] */
  #define zwidth (_graphite_analyser39->_parameters.zwidth)
  #define yheight (_graphite_analyser39->_parameters.yheight)
  #define mosaic (_graphite_analyser39->_parameters.mosaic)
  #define dspread (_graphite_analyser39->_parameters.dspread)
  #define Q (_graphite_analyser39->_parameters.Q)
  #define DM (_graphite_analyser39->_parameters.DM)
  #define pThreshold (_graphite_analyser39->_parameters.pThreshold)
  #define Rup (_graphite_analyser39->_parameters.Rup)
  #define Rdown (_graphite_analyser39->_parameters.Rdown)
  #define debug (_graphite_analyser39->_parameters.debug)
  #define mos_rms (_graphite_analyser39->_parameters.mos_rms)
  #define d_rms (_graphite_analyser39->_parameters.d_rms)
  #define mono_Q (_graphite_analyser39->_parameters.mono_Q)
  #define FN (_graphite_analyser39->_parameters.FN)
  #define FM (_graphite_analyser39->_parameters.FM)
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_graphite_analyser39->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _graphite_analyser39->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_graphite_analyser39->_position_relative, _graphite_analyser39->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      if (( instrument->_parameters._LAMBDA > 0 )) // conditional WHEN execution
      class_Monochromator_pol_trace(_graphite_analyser39, _particle); /* contains EXTEND code */
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
      // GROUP graphiteAnalyzer: from graphite_analyser0 [28] to graphite_analyser39 [67]
      if (SCATTERED) _particle->_index = 67; // when SCATTERED in GROUP: reach exit of GROUP after graphite_analyser39
      else ABSORB;     // not SCATTERED at end of GROUP: removes left events
      _particle->_index++;
    } /* end component graphite_analyser39 [67] */
    if (!ABSORBED && _particle->_index == 68) {
      /* component armHelp=Arm() [68] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_armHelp->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _armHelp->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_armHelp->_position_relative, _armHelp->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      _particle->_index++;
    } /* end component armHelp [68] */
    if (!ABSORBED && _particle->_index == 69) {
      /* component armDetector=Arm() [69] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_armDetector->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _armDetector->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_armDetector->_position_relative, _armDetector->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      _particle->_index++;
    } /* end component armDetector [69] */
    if (!ABSORBED && _particle->_index == 70) {
      /* component lamDetector=L_monitor() [70] */
  #define nL (_lamDetector->_parameters.nL)
  #define filename (_lamDetector->_parameters.filename)
  #define xmin (_lamDetector->_parameters.xmin)
  #define xmax (_lamDetector->_parameters.xmax)
  #define ymin (_lamDetector->_parameters.ymin)
  #define ymax (_lamDetector->_parameters.ymax)
  #define xwidth (_lamDetector->_parameters.xwidth)
  #define yheight (_lamDetector->_parameters.yheight)
  #define Lmin (_lamDetector->_parameters.Lmin)
  #define Lmax (_lamDetector->_parameters.Lmax)
  #define restore_neutron (_lamDetector->_parameters.restore_neutron)
  #define L_N (_lamDetector->_parameters.L_N)
  #define L_p (_lamDetector->_parameters.L_p)
  #define L_p2 (_lamDetector->_parameters.L_p2)
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_lamDetector->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _lamDetector->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_lamDetector->_position_relative, _lamDetector->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      if (( instrument->_parameters._LAMBDA > 0 )) // conditional WHEN execution
      class_L_monitor_trace(_lamDetector, _particle);
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
    } /* end component lamDetector [70] */
    if (!ABSORBED && _particle->_index == 71) {
      /* component psdDetector=PSD_monitor() [71] */
  #define nx (_psdDetector->_parameters.nx)
  #define ny (_psdDetector->_parameters.ny)
  #define filename (_psdDetector->_parameters.filename)
  #define xmin (_psdDetector->_parameters.xmin)
  #define xmax (_psdDetector->_parameters.xmax)
  #define ymin (_psdDetector->_parameters.ymin)
  #define ymax (_psdDetector->_parameters.ymax)
  #define xwidth (_psdDetector->_parameters.xwidth)
  #define yheight (_psdDetector->_parameters.yheight)
  #define restore_neutron (_psdDetector->_parameters.restore_neutron)
  #define PSD_N (_psdDetector->_parameters.PSD_N)
  #define PSD_p (_psdDetector->_parameters.PSD_p)
  #define PSD_p2 (_psdDetector->_parameters.PSD_p2)
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_psdDetector->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _psdDetector->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_psdDetector->_position_relative, _psdDetector->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      if (( instrument->_parameters._LAMBDA > 0 )) // conditional WHEN execution
      class_PSD_monitor_trace(_psdDetector, _particle);
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
    } /* end component psdDetector [71] */
    if (!ABSORBED && _particle->_index == 72) {
      /* component tofDetector=TOF_monitor() [72] */
  #define nt (_tofDetector->_parameters.nt)
  #define filename (_tofDetector->_parameters.filename)
  #define xmin (_tofDetector->_parameters.xmin)
  #define xmax (_tofDetector->_parameters.xmax)
  #define ymin (_tofDetector->_parameters.ymin)
  #define ymax (_tofDetector->_parameters.ymax)
  #define xwidth (_tofDetector->_parameters.xwidth)
  #define yheight (_tofDetector->_parameters.yheight)
  #define tmin (_tofDetector->_parameters.tmin)
  #define tmax (_tofDetector->_parameters.tmax)
  #define dt (_tofDetector->_parameters.dt)
  #define restore_neutron (_tofDetector->_parameters.restore_neutron)
  #define TOF_N (_tofDetector->_parameters.TOF_N)
  #define TOF_p (_tofDetector->_parameters.TOF_p)
  #define TOF_p2 (_tofDetector->_parameters.TOF_p2)
  #define t_min (_tofDetector->_parameters.t_min)
  #define t_max (_tofDetector->_parameters.t_max)
  #define delta_t (_tofDetector->_parameters.delta_t)
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_tofDetector->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _tofDetector->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_tofDetector->_position_relative, _tofDetector->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      if (( instrument->_parameters._LAMBDA > 0 )) // conditional WHEN execution
      class_TOF_monitor_trace(_tofDetector, _particle);
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
    } /* end component tofDetector [72] */
    if (_particle->_index > 72)
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
* instrument 'ISIS_OSIRIS' and components SAVE
***************************************************************************** */

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



int save(FILE *handle) { /* called by mccode_main for ISIS_OSIRIS:SAVE */
  if (!handle) siminfo_init(NULL);

  /* call iteratively all components SAVE */


  class_TOF_monitor_save(_tofSource);


  class_L_monitor_save(_lamStart);

  class_PSD_monitor_save(_psdStart);















  class_L_monitor_save(_lamTarget);

  class_PSD_monitor_save(_psdTarget);
















































  class_L_monitor_save(_lamDetector);

  class_PSD_monitor_save(_psdDetector);

  class_TOF_monitor_save(_tofDetector);

  if (!handle) siminfo_close(); 

  return(0);
} /* save */

/* *****************************************************************************
* instrument 'ISIS_OSIRIS' and components FINALLY
***************************************************************************** */

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



int finally(void) { /* called by mccode_main for ISIS_OSIRIS:FINALLY */
  siminfo_init(NULL);
  save(siminfo_file); /* save data when simulation ends */

  /* call iteratively all components FINALLY */


  class_TOF_monitor_finally(_tofSource);


  class_L_monitor_finally(_lamStart);

  class_PSD_monitor_finally(_psdStart);















  class_L_monitor_finally(_lamTarget);

  class_PSD_monitor_finally(_psdTarget);
















































  class_L_monitor_finally(_lamDetector);

  class_PSD_monitor_finally(_psdDetector);

  class_TOF_monitor_finally(_tofDetector);

  siminfo_close(); 

  return(0);
} /* finally */

/* *****************************************************************************
* instrument 'ISIS_OSIRIS' and components DISPLAY
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

_class_Guide_curved *class_Guide_curved_display(_class_Guide_curved *_comp
) {
  #define w1 (_comp->_parameters.w1)
  #define h1 (_comp->_parameters.h1)
  #define l (_comp->_parameters.l)
  #define R0 (_comp->_parameters.R0)
  #define Qc (_comp->_parameters.Qc)
  #define alpha (_comp->_parameters.alpha)
  #define m (_comp->_parameters.m)
  #define W (_comp->_parameters.W)
  #define curvature (_comp->_parameters.curvature)
  double x1, x2, z1, z2;
  double xplot1[100], xplot2[100], zplot1[100], zplot2[100];
  int n = 100;
  int j = 1;
  double R1 = (curvature - 0.5*w1);    /* radius of inside arc */
  double R2 = (curvature + 0.5*w1);    /* radius of outside arc */

  

  for(j=0;j<n;j++) {
    z1 = ((double)j)*(R1*l/curvature)/(double)(n - 1);
    z2 = ((double)j)*(R2*l/curvature)/(double)(n - 1);
    x1 = curvature - sqrt(R1*R1 - z1*z1);
    x2 = curvature - sqrt(R2*R2 - z2*z2);
    xplot1[j] = x1;
    xplot2[j] = x2;
    zplot1[j] = z1;
    zplot2[j] = z2;
  }
  line(xplot1[0], 0.5*h1,zplot1[0],xplot2[0], 0.5*h1,zplot2[0]);
  line(xplot1[0], 0.5*h1,zplot1[0],xplot1[0],-0.5*h1,zplot1[0]);
  line(xplot2[0],-0.5*h1,zplot2[0],xplot2[0], 0.5*h1,zplot2[0]);
  line(xplot1[0],-0.5*h1,zplot1[0],xplot2[0],-0.5*h1,zplot2[0]);
  for(j=0;j<n-1;j++) {
    line(xplot1[j],  0.5*h1, zplot1[j], xplot1[j+1],  0.5*h1, zplot1[j+1]);
    line(xplot2[j],  0.5*h1, zplot2[j], xplot2[j+1],  0.5*h1, zplot2[j+1]);
    line(xplot1[j], -0.5*h1, zplot1[j], xplot1[j+1], -0.5*h1, zplot1[j+1]);
    line(xplot2[j], -0.5*h1, zplot2[j], xplot2[j+1], -0.5*h1, zplot2[j+1]);
  }
  line(xplot1[n-1], 0.5*h1,zplot1[n-1],xplot2[n-1], 0.5*h1,zplot2[n-1]);
  line(xplot1[n-1], 0.5*h1,zplot1[n-1],xplot1[n-1],-0.5*h1,zplot1[n-1]);
  line(xplot2[n-1],-0.5*h1,zplot2[n-1],xplot2[n-1], 0.5*h1,zplot2[n-1]);
  line(xplot1[n-1],-0.5*h1,zplot1[n-1],xplot2[n-1],-0.5*h1,zplot2[n-1]);
  #undef w1
  #undef h1
  #undef l
  #undef R0
  #undef Qc
  #undef alpha
  #undef m
  #undef W
  #undef curvature
  return(_comp);
} /* class_Guide_curved_display */

_class_Incoherent *class_Incoherent_display(_class_Incoherent *_comp
) {
  #define geometry (_comp->_parameters.geometry)
  #define radius (_comp->_parameters.radius)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define zdepth (_comp->_parameters.zdepth)
  #define thickness (_comp->_parameters.thickness)
  #define target_x (_comp->_parameters.target_x)
  #define target_y (_comp->_parameters.target_y)
  #define target_z (_comp->_parameters.target_z)
  #define focus_r (_comp->_parameters.focus_r)
  #define focus_xw (_comp->_parameters.focus_xw)
  #define focus_yh (_comp->_parameters.focus_yh)
  #define focus_aw (_comp->_parameters.focus_aw)
  #define focus_ah (_comp->_parameters.focus_ah)
  #define target_index (_comp->_parameters.target_index)
  #define pack (_comp->_parameters.pack)
  #define p_interact (_comp->_parameters.p_interact)
  #define f_QE (_comp->_parameters.f_QE)
  #define gamma (_comp->_parameters.gamma)
  #define sigma_abs (_comp->_parameters.sigma_abs)
  #define sigma_inc (_comp->_parameters.sigma_inc)
  #define Vc (_comp->_parameters.Vc)
  #define concentric (_comp->_parameters.concentric)
  #define order (_comp->_parameters.order)
  #define VarsInc (_comp->_parameters.VarsInc)
  #define offdata (_comp->_parameters.offdata)
  
  if (geometry && strlen(geometry) && strcmp(geometry, "NULL") && strcmp(geometry, "0")) {	/* OFF file */
    off_display(offdata);
  }
  else if (radius > 0 &&  yheight) {	/* cylinder */
      cylinder(0,0,0,radius,yheight,0, 0,1,0);
      cylinder(0,0,0,radius-thickness,yheight,0, 0,1,0);
  }
  else if (xwidth && yheight) { 	/* box/rectangle XY */
      box(0,0,0,xwidth, yheight, zdepth);
  }
  else if (radius > 0 && !yheight) {	/* sphere */
      sphere(0,0,0,radius,0);
  }
  #undef geometry
  #undef radius
  #undef xwidth
  #undef yheight
  #undef zdepth
  #undef thickness
  #undef target_x
  #undef target_y
  #undef target_z
  #undef focus_r
  #undef focus_xw
  #undef focus_yh
  #undef focus_aw
  #undef focus_ah
  #undef target_index
  #undef pack
  #undef p_interact
  #undef f_QE
  #undef gamma
  #undef sigma_abs
  #undef sigma_inc
  #undef Vc
  #undef concentric
  #undef order
  #undef VarsInc
  #undef offdata
  return(_comp);
} /* class_Incoherent_display */

_class_Beamstop *class_Beamstop_display(_class_Beamstop *_comp
) {
  #define xmin (_comp->_parameters.xmin)
  #define xmax (_comp->_parameters.xmax)
  #define ymin (_comp->_parameters.ymin)
  #define ymax (_comp->_parameters.ymax)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define radius (_comp->_parameters.radius)
  
  if (radius != 0)
    circle("xy", 0, 0, 0, radius);
  else
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
  #undef radius
  return(_comp);
} /* class_Beamstop_display */

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


  #undef magnify
  #undef line
  #undef dashed_line
  #undef multiline
  #undef rectangle
  #undef box
  #undef circle
  #undef cylinder
  #undef sphere

int display(void) { /* called by mccode_main for ISIS_OSIRIS:DISPLAY */
  printf("MCDISPLAY: start\n");

  /* call iteratively all components DISPLAY */
  class_Arm_display(_armStart);

  class_ISIS_moderator_display(_isis_mod);

  class_TOF_monitor_display(_tofSource);

  class_Arm_display(_armStraight1);

  class_L_monitor_display(_lamStart);

  class_PSD_monitor_display(_psdStart);

  class_Guide_display(_guide_straight_1);

  class_Arm_display(_armChopper1);

  class_DiskChopper_display(_chopper1);

  class_Arm_display(_armCurved1);

  class_Guide_curved_display(_guide_curved1);

  class_Arm_display(_armChopper2);

  class_DiskChopper_display(_chopper2);

  class_Arm_display(_armCurved2);

  class_Guide_curved_display(_guide_curved2);

  class_Arm_display(_armStraight2);

  class_Guide_display(_guide_straight_2);

  class_Arm_display(_armSuper);

  class_Guide_display(_guide_super);

  class_Arm_display(_armTarget);

  class_L_monitor_display(_lamTarget);

  class_PSD_monitor_display(_psdTarget);

  class_Incoherent_display(_vsample1);

  class_Incoherent_display(_vsample2);

  class_Incoherent_display(_vsample3);

  class_Beamstop_display(_beamStop);

  class_Arm_display(_armAnalyzer);

  class_Monochromator_pol_display(_graphite_analyser0);

  class_Monochromator_pol_display(_graphite_analyser1);

  class_Monochromator_pol_display(_graphite_analyser2);

  class_Monochromator_pol_display(_graphite_analyser3);

  class_Monochromator_pol_display(_graphite_analyser4);

  class_Monochromator_pol_display(_graphite_analyser5);

  class_Monochromator_pol_display(_graphite_analyser6);

  class_Monochromator_pol_display(_graphite_analyser7);

  class_Monochromator_pol_display(_graphite_analyser8);

  class_Monochromator_pol_display(_graphite_analyser9);

  class_Monochromator_pol_display(_graphite_analyser10);

  class_Monochromator_pol_display(_graphite_analyser11);

  class_Monochromator_pol_display(_graphite_analyser12);

  class_Monochromator_pol_display(_graphite_analyser13);

  class_Monochromator_pol_display(_graphite_analyser14);

  class_Monochromator_pol_display(_graphite_analyser15);

  class_Monochromator_pol_display(_graphite_analyser16);

  class_Monochromator_pol_display(_graphite_analyser17);

  class_Monochromator_pol_display(_graphite_analyser18);

  class_Monochromator_pol_display(_graphite_analyser19);

  class_Monochromator_pol_display(_graphite_analyser20);

  class_Monochromator_pol_display(_graphite_analyser21);

  class_Monochromator_pol_display(_graphite_analyser22);

  class_Monochromator_pol_display(_graphite_analyser23);

  class_Monochromator_pol_display(_graphite_analyser24);

  class_Monochromator_pol_display(_graphite_analyser25);

  class_Monochromator_pol_display(_graphite_analyser26);

  class_Monochromator_pol_display(_graphite_analyser27);

  class_Monochromator_pol_display(_graphite_analyser28);

  class_Monochromator_pol_display(_graphite_analyser29);

  class_Monochromator_pol_display(_graphite_analyser30);

  class_Monochromator_pol_display(_graphite_analyser31);

  class_Monochromator_pol_display(_graphite_analyser32);

  class_Monochromator_pol_display(_graphite_analyser33);

  class_Monochromator_pol_display(_graphite_analyser34);

  class_Monochromator_pol_display(_graphite_analyser35);

  class_Monochromator_pol_display(_graphite_analyser36);

  class_Monochromator_pol_display(_graphite_analyser37);

  class_Monochromator_pol_display(_graphite_analyser38);

  class_Monochromator_pol_display(_graphite_analyser39);

  class_Arm_display(_armHelp);

  class_Arm_display(_armDetector);

  class_L_monitor_display(_lamDetector);

  class_PSD_monitor_display(_psdDetector);

  class_TOF_monitor_display(_tofDetector);

  printf("MCDISPLAY: end\n");

  return(0);
} /* display */

void* _getvar_parameters(char* compname)
/* enables settings parameters based use of the GETPAR macro */
{
  if (!strcmp(compname, "armStart")) return (void *) &(_armStart->_parameters);
  if (!strcmp(compname, "isis_mod")) return (void *) &(_isis_mod->_parameters);
  if (!strcmp(compname, "tofSource")) return (void *) &(_tofSource->_parameters);
  if (!strcmp(compname, "armStraight1")) return (void *) &(_armStraight1->_parameters);
  if (!strcmp(compname, "lamStart")) return (void *) &(_lamStart->_parameters);
  if (!strcmp(compname, "psdStart")) return (void *) &(_psdStart->_parameters);
  if (!strcmp(compname, "guide_straight_1")) return (void *) &(_guide_straight_1->_parameters);
  if (!strcmp(compname, "armChopper1")) return (void *) &(_armChopper1->_parameters);
  if (!strcmp(compname, "chopper1")) return (void *) &(_chopper1->_parameters);
  if (!strcmp(compname, "armCurved1")) return (void *) &(_armCurved1->_parameters);
  if (!strcmp(compname, "guide_curved1")) return (void *) &(_guide_curved1->_parameters);
  if (!strcmp(compname, "armChopper2")) return (void *) &(_armChopper2->_parameters);
  if (!strcmp(compname, "chopper2")) return (void *) &(_chopper2->_parameters);
  if (!strcmp(compname, "armCurved2")) return (void *) &(_armCurved2->_parameters);
  if (!strcmp(compname, "guide_curved2")) return (void *) &(_guide_curved2->_parameters);
  if (!strcmp(compname, "armStraight2")) return (void *) &(_armStraight2->_parameters);
  if (!strcmp(compname, "guide_straight_2")) return (void *) &(_guide_straight_2->_parameters);
  if (!strcmp(compname, "armSuper")) return (void *) &(_armSuper->_parameters);
  if (!strcmp(compname, "guide_super")) return (void *) &(_guide_super->_parameters);
  if (!strcmp(compname, "armTarget")) return (void *) &(_armTarget->_parameters);
  if (!strcmp(compname, "lamTarget")) return (void *) &(_lamTarget->_parameters);
  if (!strcmp(compname, "psdTarget")) return (void *) &(_psdTarget->_parameters);
  if (!strcmp(compname, "vsample1")) return (void *) &(_vsample1->_parameters);
  if (!strcmp(compname, "vsample2")) return (void *) &(_vsample2->_parameters);
  if (!strcmp(compname, "vsample3")) return (void *) &(_vsample3->_parameters);
  if (!strcmp(compname, "beamStop")) return (void *) &(_beamStop->_parameters);
  if (!strcmp(compname, "armAnalyzer")) return (void *) &(_armAnalyzer->_parameters);
  if (!strcmp(compname, "graphite_analyser0")) return (void *) &(_graphite_analyser0->_parameters);
  if (!strcmp(compname, "graphite_analyser1")) return (void *) &(_graphite_analyser1->_parameters);
  if (!strcmp(compname, "graphite_analyser2")) return (void *) &(_graphite_analyser2->_parameters);
  if (!strcmp(compname, "graphite_analyser3")) return (void *) &(_graphite_analyser3->_parameters);
  if (!strcmp(compname, "graphite_analyser4")) return (void *) &(_graphite_analyser4->_parameters);
  if (!strcmp(compname, "graphite_analyser5")) return (void *) &(_graphite_analyser5->_parameters);
  if (!strcmp(compname, "graphite_analyser6")) return (void *) &(_graphite_analyser6->_parameters);
  if (!strcmp(compname, "graphite_analyser7")) return (void *) &(_graphite_analyser7->_parameters);
  if (!strcmp(compname, "graphite_analyser8")) return (void *) &(_graphite_analyser8->_parameters);
  if (!strcmp(compname, "graphite_analyser9")) return (void *) &(_graphite_analyser9->_parameters);
  if (!strcmp(compname, "graphite_analyser10")) return (void *) &(_graphite_analyser10->_parameters);
  if (!strcmp(compname, "graphite_analyser11")) return (void *) &(_graphite_analyser11->_parameters);
  if (!strcmp(compname, "graphite_analyser12")) return (void *) &(_graphite_analyser12->_parameters);
  if (!strcmp(compname, "graphite_analyser13")) return (void *) &(_graphite_analyser13->_parameters);
  if (!strcmp(compname, "graphite_analyser14")) return (void *) &(_graphite_analyser14->_parameters);
  if (!strcmp(compname, "graphite_analyser15")) return (void *) &(_graphite_analyser15->_parameters);
  if (!strcmp(compname, "graphite_analyser16")) return (void *) &(_graphite_analyser16->_parameters);
  if (!strcmp(compname, "graphite_analyser17")) return (void *) &(_graphite_analyser17->_parameters);
  if (!strcmp(compname, "graphite_analyser18")) return (void *) &(_graphite_analyser18->_parameters);
  if (!strcmp(compname, "graphite_analyser19")) return (void *) &(_graphite_analyser19->_parameters);
  if (!strcmp(compname, "graphite_analyser20")) return (void *) &(_graphite_analyser20->_parameters);
  if (!strcmp(compname, "graphite_analyser21")) return (void *) &(_graphite_analyser21->_parameters);
  if (!strcmp(compname, "graphite_analyser22")) return (void *) &(_graphite_analyser22->_parameters);
  if (!strcmp(compname, "graphite_analyser23")) return (void *) &(_graphite_analyser23->_parameters);
  if (!strcmp(compname, "graphite_analyser24")) return (void *) &(_graphite_analyser24->_parameters);
  if (!strcmp(compname, "graphite_analyser25")) return (void *) &(_graphite_analyser25->_parameters);
  if (!strcmp(compname, "graphite_analyser26")) return (void *) &(_graphite_analyser26->_parameters);
  if (!strcmp(compname, "graphite_analyser27")) return (void *) &(_graphite_analyser27->_parameters);
  if (!strcmp(compname, "graphite_analyser28")) return (void *) &(_graphite_analyser28->_parameters);
  if (!strcmp(compname, "graphite_analyser29")) return (void *) &(_graphite_analyser29->_parameters);
  if (!strcmp(compname, "graphite_analyser30")) return (void *) &(_graphite_analyser30->_parameters);
  if (!strcmp(compname, "graphite_analyser31")) return (void *) &(_graphite_analyser31->_parameters);
  if (!strcmp(compname, "graphite_analyser32")) return (void *) &(_graphite_analyser32->_parameters);
  if (!strcmp(compname, "graphite_analyser33")) return (void *) &(_graphite_analyser33->_parameters);
  if (!strcmp(compname, "graphite_analyser34")) return (void *) &(_graphite_analyser34->_parameters);
  if (!strcmp(compname, "graphite_analyser35")) return (void *) &(_graphite_analyser35->_parameters);
  if (!strcmp(compname, "graphite_analyser36")) return (void *) &(_graphite_analyser36->_parameters);
  if (!strcmp(compname, "graphite_analyser37")) return (void *) &(_graphite_analyser37->_parameters);
  if (!strcmp(compname, "graphite_analyser38")) return (void *) &(_graphite_analyser38->_parameters);
  if (!strcmp(compname, "graphite_analyser39")) return (void *) &(_graphite_analyser39->_parameters);
  if (!strcmp(compname, "armHelp")) return (void *) &(_armHelp->_parameters);
  if (!strcmp(compname, "armDetector")) return (void *) &(_armDetector->_parameters);
  if (!strcmp(compname, "lamDetector")) return (void *) &(_lamDetector->_parameters);
  if (!strcmp(compname, "psdDetector")) return (void *) &(_psdDetector->_parameters);
  if (!strcmp(compname, "tofDetector")) return (void *) &(_tofDetector->_parameters);
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

/* end of generated C code ISIS_OSIRIS.c */
