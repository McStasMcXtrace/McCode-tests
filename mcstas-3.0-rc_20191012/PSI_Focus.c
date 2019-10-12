/* Automatically generated file. Do not edit. 
 * Format:     ANSI C source code
 * Creator:    McStas <http://www.mcstas.org>
 * Instrument: PSI_Focus.instr (PSI_Focus)
 * Date:       Sat Oct 12 09:04:27 2019
 * File:       PSI_Focus.c
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
* Start of instrument 'PSI_Focus' generated code
***************************************************************************** */

#ifdef MC_TRACE_ENABLED
int traceenabled = 1;
#else
int traceenabled = 0;
#endif
#define MCSTAS "/usr/share/mcstas/3.0-dev/"
int   defaultmain         = 1;
char  instrument_name[]   = "PSI_Focus";
char  instrument_source[] = "PSI_Focus.instr";
char *instrument_exe      = NULL; /* will be set to argv[0] in main */
char  instrument_code[]   = "Instrument PSI_Focus source code PSI_Focus.instr is not embedded in this executable.\n  Use --source option when running McStas.\n";

int main(int argc, char *argv[]){return mccode_main(argc, argv);}

/* *****************************************************************************
* instrument 'PSI_Focus' and components DECLARE
***************************************************************************** */

/* Instrument parameters: structure and a table for the initialisation
   (Used in e.g. inputparse and I/O function (e.g. detector_out) */

struct _struct_instrument_parameters {
  MCNUM _lambda;
  MCNUM _chopp_ratio;
  MCNUM _DET;
};
typedef struct _struct_instrument_parameters _class_instrument_parameters;

/* instrument SPLIT, GROUP and JUMP control logic */
struct instrument_logic_struct {
  long Split_mono; /* this is the SPLIT counter decremented down to 0 */
  _class_particle Split_mono_particle; /* this is the particle to duplicate */
  long Split_Sample; /* this is the SPLIT counter decremented down to 0 */
  _class_particle Split_Sample_particle; /* this is the particle to duplicate */
};

struct _instrument_struct {
  char   _name[256]; /* the name of this instrument e.g. 'PSI_Focus' */
/* Counters per component instance */
  double counter_AbsorbProp[46]; /* absorbed events in PROP routines */
  double counter_N[46], counter_P[46], counter_P2[46]; /* event counters after each component instance */
  _class_particle _trajectory[46]; /* current trajectory for STORE/RESTORE */
/* Components position table (absolute and relative coords) */
  Coords _position_relative[46]; /* positions of all components */
  Coords _position_absolute[46];
  _class_instrument_parameters _parameters; /* instrument parameters */
  struct instrument_logic_struct logic; /* instrument logic */
} _instrument_var;
struct _instrument_struct *instrument = & _instrument_var;
#pragma acc declare create ( _instrument_var )
#pragma acc declare create ( instrument )

int numipar = 3;
struct mcinputtable_struct mcinputtable[] = {
  "lambda", &(_instrument_var._parameters._lambda), instr_type_double, "3.4", 
  "chopp_ratio", &(_instrument_var._parameters._chopp_ratio), instr_type_double, "1", 
  "DET", &(_instrument_var._parameters._DET), instr_type_double, "-69.9", 
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


/* Shared user declarations for all components types 'Guide'. */

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


/* Shared user declarations for all components types 'Bender'. */


/* Shared user declarations for all components types 'Guide_channeled'. */


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


/* Shared user declarations for all components types 'Monochromator_2foc'. */

#ifndef DIV_CUTOFF
#define DIV_CUTOFF 2            /* ~ 10^-5 cutoff. */
#endif

/* Shared user declarations for all components types 'FermiChopper'. */

#ifndef FermiChopper_TimeAccuracy
#define FermiChopper_TimeAccuracy 1e-9
#define FermiChopper_MAXITER      100

/* Definition of internal variable structure: all counters */
struct FermiChopper_struct
{
double omega, ph0, t0;  /* chopper rotation */
double C_slit;          /* slit curvature radius in [m] */
double L_slit;          /* slit package length [m] */
double sum_t;
double sum_v;
double sum_N;
double sum_N_pass;
/* events */
long absorb_alreadyinside;
long absorb_topbottom;
long absorb_cylentrance;
long absorb_sideentrance;
long absorb_notreachentrance;
long absorb_packentrance;
long absorb_slitcoating;
long warn_notreachslitwall;
long absorb_exitslitpack;
long absorb_maxiterations;
long absorb_wrongdirection;
long absorb_nocontrol;
long absorb_cylexit;
long warn_notreachslitoutput;
char compcurname[256];
double X[5],Z[5]; /* position of interaction locations in the Fermi chopper rotating frame
[0]=cylinder input
[1]=slit input
[2]=slit wall reflection
[3]=slit exit
[4]=cylinder exit
*/
};

/*****************************************************************************
* FC_zrot: returns Z' in rotating frame, from X,Z and t,omega,ph0
****************************************************************************/
double FC_zrot(double X, double Z, double T, struct FermiChopper_struct FCs)
{
  double omega =FCs.omega;
  double ph0   =FCs.ph0;

  return( Z*cos(omega*T+ph0)-X*sin(omega*T+ph0) );
}

/*****************************************************************************
 * FC_xrot: returns X' in rotating frame, from X,Z and omega,t,ph0
 *          additional coordinate shift in case of curved slits
 ****************************************************************************/
double FC_xrot(double X, double Z, double T, struct FermiChopper_struct FCs)
{
  double omega =FCs.omega;
  double ph0   =FCs.ph0;
  double C_slit=FCs.C_slit;
  double ret, tmp;

  ret = X*cos(omega*T+ph0)+Z*sin(omega*T+ph0);

  if (C_slit) {
    tmp  = fabs(FC_zrot(X, Z, T, FCs));
    if (tmp < FCs.L_slit/2) {
      tmp  = (FCs.L_slit/2 - tmp)*C_slit;
      ret += (1-sqrt(1-tmp*tmp))/C_slit;
    }
  }
  return( ret );
}

/*****************************************************************************
 * FC_xzrot_dt(x,z,vx,vz, t,dt, type='x' or 'z', FCs)
 *   returns X' or Z' in rotating frame, from X,Z and t,omega,ph0
 *              taking into account propagation with velocity during time dt
 ****************************************************************************/
double FC_xzrot_dt(double x, double z, double vx, double vz,
                   double t, double dt, char type, struct FermiChopper_struct FCs)
{
  if (dt) /* with propagation */
    return( (type == 'x' ? FC_xrot(x+vx*dt, z+vz*dt, t+dt, FCs)
                         : FC_zrot(x+vx*dt, z+vz*dt, t+dt, FCs)) );
  else    /* without propagation */
    return( (type == 'x' ? FC_xrot(x,z,t,FCs)
                         : FC_zrot(x,z,t,FCs)) );
}

/*****************************************************************************
 * FC_xzbrent(x,z,vx,vz, t,dt, type='x' or 'z', d, FCs)
 *   solves X'=d and Z'=d with Brent algorithm in time interval [0, dt].
 *           Returns time within [0,dt], from NumRecip in C, chap 9, p360 (zbrent)
 *           ERRORS: return -1 not used
 *                          -2 if exceed MAX iteration
 *                          -3 no sign change in range
 ****************************************************************************/
double FC_xzbrent(double x, double z, double vx, double vz,
                  double t, double dt,
                  char type, double d, struct FermiChopper_struct FCs)
{
  int iter;
  double a=0,b=dt;
  double c,dd,e,min1,min2;
  double tol=FermiChopper_TimeAccuracy;
  double EPS=FermiChopper_TimeAccuracy;
  double fa=FC_xzrot_dt(x,z,vx,vz, t,a, type, FCs) - d;
  double fb=FC_xzrot_dt(x,z,vx,vz, t,b, type, FCs) - d;
  double fc,p,q,r,s,tol1,xm;

  if (fb*fa > 0.0) return -3;
  fc=fb;
  for (iter=1;iter<=FermiChopper_MAXITER;iter++) {
    if (fb*fc > 0.0) {
      c=a;
      fc=fa;
      e=dd=b-a;
    }
    if (fabs(fc) < fabs(fb)) {
      a=b;
      b=c;
      c=a;
      fa=fb;
      fb=fc;
      fc=fa;
    }
    tol1=2.0*EPS*fabs(b)+0.5*tol;
    xm=0.5*(c-b);
    if (fabs(xm) <= tol1 || fb == 0.0) return b;
    if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
      s=fb/fa;
      if (a == c) {
        p=2.0*xm*s;
        q=1.0-s;
      } else {
        q=fa/fc;
        r=fb/fc;
        p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
        q=(q-1.0)*(r-1.0)*(s-1.0);
      }
      if (p > 0.0)  q = -q;
      p=fabs(p);
      min1=3.0*xm*q-fabs(tol1*q);
      min2=fabs(e*q);
      if (2.0*p < (min1 < min2 ? min1 : min2)) {
        e=dd;
        dd=p/q;
      } else {
        dd=xm;
        e=dd;
      }
    } else {
      dd=xm;
      e=dd;
    }
    a=b;
    fa=fb;
    if (fabs(dd) > tol1)
      b += dd;
    else
      b += (xm > 0.0 ? fabs(tol1) : -fabs(tol1));
    fb=FC_xzrot_dt(x,z,vx,vz, t,b, type, FCs) - d;
  }
  return -2;
} /* FC_xzbrent */

/*****************************************************************************
 * Wrappers to intersection algorithms
 ****************************************************************************/
double FC_xintersect(double x, double z, double vx, double vz,
                   double t, double dt,
                   double d, struct FermiChopper_struct FCs)
{
  return(FC_xzbrent(x, z, vx, vz, t, dt, 'x', d, FCs));
}
double FC_zintersect(double x, double z, double vx, double vz,
                   double t, double dt,
                   double d, struct FermiChopper_struct FCs)
{
  return(FC_xzbrent(x, z, vx, vz, t, dt, 'z', d, FCs));
}

#endif

/* Shared user declarations for all components types 'V_sample'. */
struct StructVarsV
{
double  sigma_a; /* Absorption cross section per atom (barns) */
    double  sigma_i; /* Incoherent scattering cross section per atom (barns) */
    double  rho;     /* Density of atoms (AA-3) */
    double  my_s;
    double  my_a_v;
    int     shapetyp;    /* 0 double well cylynder, 1 box,  3 sphere */
    double  distance;    /* when non zero, gives rect target distance */
    double  aw,ah;       /* rectangular angular dimensions */
    double  xw,yh;       /* rectangular metrical dimensions */
    double  tx,ty,tz;    /* target coords */
  };


/* ************************************************************************** */
/*             End of SHARE user declarations for all components              */
/* ************************************************************************** */


/* ********************** component definition declarations. **************** */

/* component a1=Progress_bar() [1] DECLARE */
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
  char     _name[256]; /* e.g. a1 */
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
_class_Progress_bar _a1_var;
_class_Progress_bar *_a1 = &_a1_var;
#pragma acc declare create ( _a1_var )
#pragma acc declare create ( _a1 )

/* component csource=Source_gen() [2] DECLARE */
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
  char     _name[256]; /* e.g. csource */
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
_class_Source_gen _csource_var;
_class_Source_gen *_csource = &_csource_var;
#pragma acc declare create ( _csource_var )
#pragma acc declare create ( _csource )

/* component guide1=Guide() [3] DECLARE */
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
  char     _name[256]; /* e.g. guide1 */
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
_class_Guide _guide1_var;
_class_Guide *_guide1 = &_guide1_var;
#pragma acc declare create ( _guide1_var )
#pragma acc declare create ( _guide1 )

/* component guide2=Bender() [4] DECLARE */
/* Parameter definition for component type 'Bender' */
struct _struct_Bender_parameters {
  /* Component type 'Bender' setting parameters */
  MCNUM w;
  MCNUM h;
  MCNUM r;
  MCNUM Win;
  MCNUM k;
  MCNUM d;
  MCNUM l;
  MCNUM R0a;
  MCNUM Qca;
  MCNUM alphaa;
  MCNUM ma;
  MCNUM Wa;
  MCNUM R0i;
  MCNUM Qci;
  MCNUM alphai;
  MCNUM mi;
  MCNUM Wi;
  MCNUM R0s;
  MCNUM Qcs;
  MCNUM alphas;
  MCNUM ms;
  MCNUM Ws;
  /* Component type 'Bender' private parameters */
  /* Component type 'Bender' DECLARE code stored as structure members */
double bk, mWin;
}; /* _struct_Bender_parameters */
typedef struct _struct_Bender_parameters _class_Bender_parameters;

/* Parameters for component type 'Bender' */
struct _struct_Bender {
  char     _name[256]; /* e.g. guide2 */
  char     _type[256]; /* Bender */
  long     _index; /* e.g. 2 index in TRACE list */
  Coords   _position_absolute;
  Coords   _position_relative; /* wrt PREVIOUS */
  Rotation _rotation_absolute;
  Rotation _rotation_relative; /* wrt PREVIOUS */
  int      _rotation_is_identity;
  _class_Bender_parameters _parameters;
};
typedef struct _struct_Bender _class_Bender;
_class_Bender _guide2_var;
_class_Bender *_guide2 = &_guide2_var;
#pragma acc declare create ( _guide2_var )
#pragma acc declare create ( _guide2 )

_class_Guide _bunker_var;
_class_Guide *_bunker = &_bunker_var;
#pragma acc declare create ( _bunker_var )
#pragma acc declare create ( _bunker )

_class_Guide _guide3_var;
_class_Guide *_guide3 = &_guide3_var;
#pragma acc declare create ( _guide3_var )
#pragma acc declare create ( _guide3 )

/* component lambdaGuideExit=L_monitor() [7] DECLARE */
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
  char     _name[256]; /* e.g. lambdaGuideExit */
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
_class_L_monitor _lambdaGuideExit_var;
_class_L_monitor *_lambdaGuideExit = &_lambdaGuideExit_var;
#pragma acc declare create ( _lambdaGuideExit_var )
#pragma acc declare create ( _lambdaGuideExit )

/* component DivMonGuideExit=Divergence_monitor() [8] DECLARE */
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
  char     _name[256]; /* e.g. DivMonGuideExit */
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
_class_Divergence_monitor _DivMonGuideExit_var;
_class_Divergence_monitor *_DivMonGuideExit = &_DivMonGuideExit_var;
#pragma acc declare create ( _DivMonGuideExit_var )
#pragma acc declare create ( _DivMonGuideExit )

/* component PSDGuideExit=PSD_monitor() [9] DECLARE */
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
  char     _name[256]; /* e.g. PSDGuideExit */
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
_class_PSD_monitor _PSDGuideExit_var;
_class_PSD_monitor *_PSDGuideExit = &_PSDGuideExit_var;
#pragma acc declare create ( _PSDGuideExit_var )
#pragma acc declare create ( _PSDGuideExit )

/* component FOCUSguide=Guide_channeled() [10] DECLARE */
/* Parameter definition for component type 'Guide_channeled' */
struct _struct_Guide_channeled_parameters {
  /* Component type 'Guide_channeled' setting parameters */
  MCNUM w1;
  MCNUM h1;
  MCNUM w2;
  MCNUM h2;
  MCNUM l;
  MCNUM R0;
  MCNUM Qc;
  MCNUM alpha;
  MCNUM m;
  MCNUM nslit;
  MCNUM d;
  MCNUM Qcx;
  MCNUM Qcy;
  MCNUM alphax;
  MCNUM alphay;
  MCNUM W;
  MCNUM mx;
  MCNUM my;
  MCNUM nu;
  MCNUM phase;
  /* Component type 'Guide_channeled' private parameters */
  /* Component type 'Guide_channeled' DECLARE code stored as structure members */
double w1c;
double w2c;
double ww, hh;
double whalf, hhalf;
double lwhalf, lhhalf;
}; /* _struct_Guide_channeled_parameters */
typedef struct _struct_Guide_channeled_parameters _class_Guide_channeled_parameters;

/* Parameters for component type 'Guide_channeled' */
struct _struct_Guide_channeled {
  char     _name[256]; /* e.g. FOCUSguide */
  char     _type[256]; /* Guide_channeled */
  long     _index; /* e.g. 2 index in TRACE list */
  Coords   _position_absolute;
  Coords   _position_relative; /* wrt PREVIOUS */
  Rotation _rotation_absolute;
  Rotation _rotation_relative; /* wrt PREVIOUS */
  int      _rotation_is_identity;
  _class_Guide_channeled_parameters _parameters;
};
typedef struct _struct_Guide_channeled _class_Guide_channeled;
_class_Guide_channeled _FOCUSguide_var;
_class_Guide_channeled *_FOCUSguide = &_FOCUSguide_var;
#pragma acc declare create ( _FOCUSguide_var )
#pragma acc declare create ( _FOCUSguide )

/* component FirstChopper=DiskChopper() [11] DECLARE */
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
  char     _name[256]; /* e.g. FirstChopper */
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
_class_DiskChopper _FirstChopper_var;
_class_DiskChopper *_FirstChopper = &_FirstChopper_var;
#pragma acc declare create ( _FirstChopper_var )
#pragma acc declare create ( _FirstChopper )

/* component DISCTOF=Monitor_nD() [12] DECLARE */
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
  char     _name[256]; /* e.g. DISCTOF */
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
_class_Monitor_nD _DISCTOF_var;
_class_Monitor_nD *_DISCTOF = &_DISCTOF_var;
#pragma acc declare create ( _DISCTOF_var )
#pragma acc declare create ( _DISCTOF )

_class_PSD_monitor _PSDmon1Chopper_var;
_class_PSD_monitor *_PSDmon1Chopper = &_PSDmon1Chopper_var;
#pragma acc declare create ( _PSDmon1Chopper_var )
#pragma acc declare create ( _PSDmon1Chopper )

/* component VacuumTube1entry=Slit() [14] DECLARE */
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
  char     _name[256]; /* e.g. VacuumTube1entry */
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
_class_Slit _VacuumTube1entry_var;
_class_Slit *_VacuumTube1entry = &_VacuumTube1entry_var;
#pragma acc declare create ( _VacuumTube1entry_var )
#pragma acc declare create ( _VacuumTube1entry )

_class_Slit _VacuumTube1exit_var;
_class_Slit *_VacuumTube1exit = &_VacuumTube1exit_var;
#pragma acc declare create ( _VacuumTube1exit_var )
#pragma acc declare create ( _VacuumTube1exit )

_class_Slit _VacuumTube2entry_var;
_class_Slit *_VacuumTube2entry = &_VacuumTube2entry_var;
#pragma acc declare create ( _VacuumTube2entry_var )
#pragma acc declare create ( _VacuumTube2entry )

_class_Slit _VacuumTube2exit_var;
_class_Slit *_VacuumTube2exit = &_VacuumTube2exit_var;
#pragma acc declare create ( _VacuumTube2exit_var )
#pragma acc declare create ( _VacuumTube2exit )

_class_Slit _VacuumTube3entry_var;
_class_Slit *_VacuumTube3entry = &_VacuumTube3entry_var;
#pragma acc declare create ( _VacuumTube3entry_var )
#pragma acc declare create ( _VacuumTube3entry )

_class_Slit _VacuumTube3exit_var;
_class_Slit *_VacuumTube3exit = &_VacuumTube3exit_var;
#pragma acc declare create ( _VacuumTube3exit_var )
#pragma acc declare create ( _VacuumTube3exit )

_class_PSD_monitor _PSDmonMono_var;
_class_PSD_monitor *_PSDmonMono = &_PSDmonMono_var;
#pragma acc declare create ( _PSDmonMono_var )
#pragma acc declare create ( _PSDmonMono )

_class_Monitor_nD _MONOTOF_var;
_class_Monitor_nD *_MONOTOF = &_MONOTOF_var;
#pragma acc declare create ( _MONOTOF_var )
#pragma acc declare create ( _MONOTOF )

_class_Divergence_monitor _DivMonMono_var;
_class_Divergence_monitor *_DivMonMono = &_DivMonMono_var;
#pragma acc declare create ( _DivMonMono_var )
#pragma acc declare create ( _DivMonMono )

/* component focus_mono=Arm() [23] DECLARE */
/* Parameter definition for component type 'Arm' */
struct _struct_Arm_parameters {
  char Arm_has_no_parameters;
}; /* _struct_Arm_parameters */
typedef struct _struct_Arm_parameters _class_Arm_parameters;

/* Parameters for component type 'Arm' */
struct _struct_Arm {
  char     _name[256]; /* e.g. focus_mono */
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
_class_Arm _focus_mono_var;
_class_Arm *_focus_mono = &_focus_mono_var;
#pragma acc declare create ( _focus_mono_var )
#pragma acc declare create ( _focus_mono )

/* component mono=Monochromator_2foc() [24] DECLARE */
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
  char     _name[256]; /* e.g. mono */
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
_class_Monochromator_2foc _mono_var;
_class_Monochromator_2foc *_mono = &_mono_var;
#pragma acc declare create ( _mono_var )
#pragma acc declare create ( _mono )

_class_Arm _a2_var;
_class_Arm *_a2 = &_a2_var;
#pragma acc declare create ( _a2_var )
#pragma acc declare create ( _a2 )

_class_Monitor_nD _FERMITOF_before_var;
_class_Monitor_nD *_FERMITOF_before = &_FERMITOF_before_var;
#pragma acc declare create ( _FERMITOF_before_var )
#pragma acc declare create ( _FERMITOF_before )

_class_L_monitor _lambdaFermi_var;
_class_L_monitor *_lambdaFermi = &_lambdaFermi_var;
#pragma acc declare create ( _lambdaFermi_var )
#pragma acc declare create ( _lambdaFermi )

/* component EMON_Fermi=E_monitor() [28] DECLARE */
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
  char     _name[256]; /* e.g. EMON_Fermi */
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
_class_E_monitor _EMON_Fermi_var;
_class_E_monitor *_EMON_Fermi = &_EMON_Fermi_var;
#pragma acc declare create ( _EMON_Fermi_var )
#pragma acc declare create ( _EMON_Fermi )

_class_Divergence_monitor _DivMonfermi1_var;
_class_Divergence_monitor *_DivMonfermi1 = &_DivMonfermi1_var;
#pragma acc declare create ( _DivMonfermi1_var )
#pragma acc declare create ( _DivMonfermi1 )

_class_PSD_monitor _PSD_Fermi1_var;
_class_PSD_monitor *_PSD_Fermi1 = &_PSD_Fermi1_var;
#pragma acc declare create ( _PSD_Fermi1_var )
#pragma acc declare create ( _PSD_Fermi1 )

/* component FoChopper=FermiChopper() [31] DECLARE */
/* Parameter definition for component type 'FermiChopper' */
struct _struct_FermiChopper_parameters {
  /* Component type 'FermiChopper' setting parameters */
  MCNUM phase;
  MCNUM radius;
  MCNUM nu;
  MCNUM w;
  MCNUM nslit;
  MCNUM R0;
  MCNUM Qc;
  MCNUM alpha;
  MCNUM m;
  MCNUM W;
  MCNUM length;
  MCNUM eff;
  MCNUM zero_time;
  MCNUM xwidth;
  MCNUM verbose;
  MCNUM yheight;
  MCNUM curvature;
  MCNUM delay;
  /* Component type 'FermiChopper' private parameters */
  /* Component type 'FermiChopper' DECLARE code stored as structure members */
  struct FermiChopper_struct FCVars;
}; /* _struct_FermiChopper_parameters */
typedef struct _struct_FermiChopper_parameters _class_FermiChopper_parameters;

/* Parameters for component type 'FermiChopper' */
struct _struct_FermiChopper {
  char     _name[256]; /* e.g. FoChopper */
  char     _type[256]; /* FermiChopper */
  long     _index; /* e.g. 2 index in TRACE list */
  Coords   _position_absolute;
  Coords   _position_relative; /* wrt PREVIOUS */
  Rotation _rotation_absolute;
  Rotation _rotation_relative; /* wrt PREVIOUS */
  int      _rotation_is_identity;
  _class_FermiChopper_parameters _parameters;
};
typedef struct _struct_FermiChopper _class_FermiChopper;
_class_FermiChopper _FoChopper_var;
_class_FermiChopper *_FoChopper = &_FoChopper_var;
#pragma acc declare create ( _FoChopper_var )
#pragma acc declare create ( _FoChopper )

_class_PSD_monitor _PSD_Fermi2_var;
_class_PSD_monitor *_PSD_Fermi2 = &_PSD_Fermi2_var;
#pragma acc declare create ( _PSD_Fermi2_var )
#pragma acc declare create ( _PSD_Fermi2 )

_class_Divergence_monitor _DivMonfermi2_var;
_class_Divergence_monitor *_DivMonfermi2 = &_DivMonfermi2_var;
#pragma acc declare create ( _DivMonfermi2_var )
#pragma acc declare create ( _DivMonfermi2 )

_class_Monitor_nD _FERMITOF1_var;
_class_Monitor_nD *_FERMITOF1 = &_FERMITOF1_var;
#pragma acc declare create ( _FERMITOF1_var )
#pragma acc declare create ( _FERMITOF1 )

_class_Slit _SAMPLE_SLIT_var;
_class_Slit *_SAMPLE_SLIT = &_SAMPLE_SLIT_var;
#pragma acc declare create ( _SAMPLE_SLIT_var )
#pragma acc declare create ( _SAMPLE_SLIT )

_class_Monitor_nD _FERMITOF2_var;
_class_Monitor_nD *_FERMITOF2 = &_FERMITOF2_var;
#pragma acc declare create ( _FERMITOF2_var )
#pragma acc declare create ( _FERMITOF2 )

_class_PSD_monitor _PSD_SAMPLE_var;
_class_PSD_monitor *_PSD_SAMPLE = &_PSD_SAMPLE_var;
#pragma acc declare create ( _PSD_SAMPLE_var )
#pragma acc declare create ( _PSD_SAMPLE )

_class_Divergence_monitor _DivMon_Sample_var;
_class_Divergence_monitor *_DivMon_Sample = &_DivMon_Sample_var;
#pragma acc declare create ( _DivMon_Sample_var )
#pragma acc declare create ( _DivMon_Sample )

_class_E_monitor _EMON_SAMPLE_var;
_class_E_monitor *_EMON_SAMPLE = &_EMON_SAMPLE_var;
#pragma acc declare create ( _EMON_SAMPLE_var )
#pragma acc declare create ( _EMON_SAMPLE )

_class_Arm _a3_var;
_class_Arm *_a3 = &_a3_var;
#pragma acc declare create ( _a3_var )
#pragma acc declare create ( _a3 )

/* component Sample=V_sample() [41] DECLARE */
/* Parameter definition for component type 'V_sample' */
struct _struct_V_sample_parameters {
  /* Component type 'V_sample' setting parameters */
  MCNUM radius;
  MCNUM thickness;
  MCNUM zdepth;
  MCNUM Vc;
  MCNUM sigma_abs;
  MCNUM sigma_inc;
  MCNUM radius_i;
  MCNUM radius_o;
  MCNUM h;
  MCNUM focus_r;
  MCNUM pack;
  MCNUM frac;
  MCNUM f_QE;
  MCNUM gamma;
  MCNUM target_x;
  MCNUM target_y;
  MCNUM target_z;
  MCNUM focus_xw;
  MCNUM focus_yh;
  MCNUM focus_aw;
  MCNUM focus_ah;
  MCNUM xwidth;
  MCNUM yheight;
  MCNUM zthick;
  MCNUM rad_sphere;
  MCNUM sig_a;
  MCNUM sig_i;
  MCNUM V0;
  long target_index;
  MCNUM multiples;
  /* Component type 'V_sample' private parameters */
  /* Component type 'V_sample' DECLARE code stored as structure members */
  struct StructVarsV VarsV;
}; /* _struct_V_sample_parameters */
typedef struct _struct_V_sample_parameters _class_V_sample_parameters;

/* Parameters for component type 'V_sample' */
struct _struct_V_sample {
  char     _name[256]; /* e.g. Sample */
  char     _type[256]; /* V_sample */
  long     _index; /* e.g. 2 index in TRACE list */
  Coords   _position_absolute;
  Coords   _position_relative; /* wrt PREVIOUS */
  Rotation _rotation_absolute;
  Rotation _rotation_relative; /* wrt PREVIOUS */
  int      _rotation_is_identity;
  _class_V_sample_parameters _parameters;
};
typedef struct _struct_V_sample _class_V_sample;
_class_V_sample _Sample_var;
_class_V_sample *_Sample = &_Sample_var;
#pragma acc declare create ( _Sample_var )
#pragma acc declare create ( _Sample )

_class_Monitor_nD _TOF_Det_var;
_class_Monitor_nD *_TOF_Det = &_TOF_Det_var;
#pragma acc declare create ( _TOF_Det_var )
#pragma acc declare create ( _TOF_Det )

_class_Monitor_nD _FoDet_var;
_class_Monitor_nD *_FoDet = &_FoDet_var;
#pragma acc declare create ( _FoDet_var )
#pragma acc declare create ( _FoDet )

_class_E_monitor _EMON_DET_var;
_class_E_monitor *_EMON_DET = &_EMON_DET_var;
#pragma acc declare create ( _EMON_DET_var )
#pragma acc declare create ( _EMON_DET )

int mcNUMCOMP = 44;

/* User declarations from instrument definition. Can define functions. */
 double lambda;
 double chopp_ratio;
 double DET;
 double PHM, TTM;
 double RV_2, RH_2;
 double FO_PHA;
 double DISC_SPEED, FERMI_SPEED;
 /* mono d-spacing */
 double monodd = 3.355;
 /* time constant [mus/A/m] */
 double c = 252.78;
 /* distance sample - detector */
 double dsd = 2.5;
 /* distance guide - monochromator */
 double dgm = 2.996;
 /* distance fermi-chopper - sample */
 double dfcs = 0.4997;
 /* distance disk chopper - monochromator */
 double ddcm = 2.944;
 /* distance  monochromator - fermi-chopper */
 double mfc = 1.002;
 /* distance  monochromator - sample */
 double dms = 1.4997;
 /* setting range of sources */
 double dL = 0.7;
 double LMIN, LMAX;
 double EMIN, EMAX;


#undef compcurname
#undef compcurtype
#undef compcurindex
/* end of instrument 'PSI_Focus' and components DECLARE */

/* *****************************************************************************
* instrument 'PSI_Focus' and components INITIALISE
***************************************************************************** */

/* component a1=Progress_bar() SETTING, POSITION/ROTATION */
int _a1_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_a1_setpos] component a1=Progress_bar() SETTING [/usr/share/mcstas/3.0-dev/misc/Progress_bar.comp:57]");
  stracpy(_a1->_name, "a1", 16384);
  stracpy(_a1->_type, "Progress_bar", 16384);
  _a1->_index=1;
  if("NULL" && strlen("NULL"))
    stracpy(_a1->_parameters.profile, "NULL" ? "NULL" : "", 16384);
  else 
  _a1->_parameters.profile[0]='\0';
  #define profile (_a1->_parameters.profile)
  _a1->_parameters.percent = 10;
  #define percent (_a1->_parameters.percent)
  _a1->_parameters.flag_save = 0;
  #define flag_save (_a1->_parameters.flag_save)
  _a1->_parameters.minutes = 0;
  #define minutes (_a1->_parameters.minutes)

  #define IntermediateCnts (_a1->_parameters.IntermediateCnts)
  #define StartTime (_a1->_parameters.StartTime)
  #define EndTime (_a1->_parameters.EndTime)
  #define CurrentTime (_a1->_parameters.CurrentTime)

  #undef profile
  #undef percent
  #undef flag_save
  #undef minutes
  #undef IntermediateCnts
  #undef StartTime
  #undef EndTime
  #undef CurrentTime
  /* component a1=Progress_bar() AT ROTATED */
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
  } /* a1=Progress_bar() AT ROTATED */
  DEBUG_COMPONENT("a1", _a1->_position_absolute, _a1->_rotation_absolute);
  instrument->_position_absolute[1] = _a1->_position_absolute;
  instrument->_position_relative[1] = _a1->_position_relative;
  instrument->counter_N[1]  = instrument->counter_P[1] = instrument->counter_P2[1] = 0;
  instrument->counter_AbsorbProp[1]= 0;
  return(0);
} /* _a1_setpos */

/* component csource=Source_gen() SETTING, POSITION/ROTATION */
int _csource_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_csource_setpos] component csource=Source_gen() SETTING [/usr/share/mcstas/3.0-dev/sources/Source_gen.comp:206]");
  stracpy(_csource->_name, "csource", 16384);
  stracpy(_csource->_type, "Source_gen", 16384);
  _csource->_index=2;
  if("NULL" && strlen("NULL"))
    stracpy(_csource->_parameters.flux_file, "NULL" ? "NULL" : "", 16384);
  else 
  _csource->_parameters.flux_file[0]='\0';
  #define flux_file (_csource->_parameters.flux_file)
  if("NULL" && strlen("NULL"))
    stracpy(_csource->_parameters.xdiv_file, "NULL" ? "NULL" : "", 16384);
  else 
  _csource->_parameters.xdiv_file[0]='\0';
  #define xdiv_file (_csource->_parameters.xdiv_file)
  if("NULL" && strlen("NULL"))
    stracpy(_csource->_parameters.ydiv_file, "NULL" ? "NULL" : "", 16384);
  else 
  _csource->_parameters.ydiv_file[0]='\0';
  #define ydiv_file (_csource->_parameters.ydiv_file)
  _csource->_parameters.radius = 0.0;
  #define radius (_csource->_parameters.radius)
  _csource->_parameters.dist = 1.77;
  #define dist (_csource->_parameters.dist)
  _csource->_parameters.focus_xw = 0.05;
  #define focus_xw (_csource->_parameters.focus_xw)
  _csource->_parameters.focus_yh = 0.12;
  #define focus_yh (_csource->_parameters.focus_yh)
  _csource->_parameters.focus_aw = 0;
  #define focus_aw (_csource->_parameters.focus_aw)
  _csource->_parameters.focus_ah = 0;
  #define focus_ah (_csource->_parameters.focus_ah)
  _csource->_parameters.E0 = 0;
  #define E0 (_csource->_parameters.E0)
  _csource->_parameters.dE = 0;
  #define dE (_csource->_parameters.dE)
  _csource->_parameters.lambda0 = 3.4;
  #define lambda0 (_csource->_parameters.lambda0)
  _csource->_parameters.dlambda = dL;
  #define dlambda (_csource->_parameters.dlambda)
  _csource->_parameters.I1 = 8.5E11;
  #define I1 (_csource->_parameters.I1)
  _csource->_parameters.yheight = 0.135;
  #define yheight (_csource->_parameters.yheight)
  _csource->_parameters.xwidth = 0.08;
  #define xwidth (_csource->_parameters.xwidth)
  _csource->_parameters.verbose = 0;
  #define verbose (_csource->_parameters.verbose)
  _csource->_parameters.T1 = 296.16;
  #define T1 (_csource->_parameters.T1)
  _csource->_parameters.flux_file_perAA = 0;
  #define flux_file_perAA (_csource->_parameters.flux_file_perAA)
  _csource->_parameters.flux_file_log = 0;
  #define flux_file_log (_csource->_parameters.flux_file_log)
  _csource->_parameters.Lmin = 0;
  #define Lmin (_csource->_parameters.Lmin)
  _csource->_parameters.Lmax = 0;
  #define Lmax (_csource->_parameters.Lmax)
  _csource->_parameters.Emin = 0;
  #define Emin (_csource->_parameters.Emin)
  _csource->_parameters.Emax = 0;
  #define Emax (_csource->_parameters.Emax)
  _csource->_parameters.T2 = 40.68;
  #define T2 (_csource->_parameters.T2)
  _csource->_parameters.I2 = 5.2E11;
  #define I2 (_csource->_parameters.I2)
  _csource->_parameters.T3 = 0;
  #define T3 (_csource->_parameters.T3)
  _csource->_parameters.I3 = 0;
  #define I3 (_csource->_parameters.I3)
  _csource->_parameters.zdepth = 0;
  #define zdepth (_csource->_parameters.zdepth)
  _csource->_parameters.target_index = + 1;
  #define target_index (_csource->_parameters.target_index)

  #define p_in (_csource->_parameters.p_in)
  #define lambda1 (_csource->_parameters.lambda1)
  #define lambda2 (_csource->_parameters.lambda2)
  #define lambda3 (_csource->_parameters.lambda3)
  #define pTable (_csource->_parameters.pTable)
  #define pTable_x (_csource->_parameters.pTable_x)
  #define pTable_y (_csource->_parameters.pTable_y)
  #define pTable_xmin (_csource->_parameters.pTable_xmin)
  #define pTable_xmax (_csource->_parameters.pTable_xmax)
  #define pTable_xsum (_csource->_parameters.pTable_xsum)
  #define pTable_ymin (_csource->_parameters.pTable_ymin)
  #define pTable_ymax (_csource->_parameters.pTable_ymax)
  #define pTable_ysum (_csource->_parameters.pTable_ysum)
  #define pTable_dxmin (_csource->_parameters.pTable_dxmin)
  #define pTable_dxmax (_csource->_parameters.pTable_dxmax)
  #define pTable_dymin (_csource->_parameters.pTable_dymin)
  #define pTable_dymax (_csource->_parameters.pTable_dymax)

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
  /* component csource=Source_gen() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (0)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _a1->_rotation_absolute, _csource->_rotation_absolute);
    rot_transpose(_a1->_rotation_absolute, tr1);
    rot_mul(_csource->_rotation_absolute, tr1, _csource->_rotation_relative);
    _csource->_rotation_is_identity =  rot_test_identity(_csource->_rotation_relative);
    tc1 = coords_set(
      0, 0, 0);
    rot_transpose(_a1->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _csource->_position_absolute = coords_add(_a1->_position_absolute, tc2);
    tc1 = coords_sub(_a1->_position_absolute, _csource->_position_absolute);
    _csource->_position_relative = rot_apply(_csource->_rotation_absolute, tc1);
  } /* csource=Source_gen() AT ROTATED */
  DEBUG_COMPONENT("csource", _csource->_position_absolute, _csource->_rotation_absolute);
  instrument->_position_absolute[2] = _csource->_position_absolute;
  instrument->_position_relative[2] = _csource->_position_relative;
  instrument->counter_N[2]  = instrument->counter_P[2] = instrument->counter_P2[2] = 0;
  instrument->counter_AbsorbProp[2]= 0;
  return(0);
} /* _csource_setpos */

/* component guide1=Guide() SETTING, POSITION/ROTATION */
int _guide1_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_guide1_setpos] component guide1=Guide() SETTING [/usr/share/mcstas/3.0-dev/optics/Guide.comp:74]");
  stracpy(_guide1->_name, "guide1", 16384);
  stracpy(_guide1->_type, "Guide", 16384);
  _guide1->_index=3;
  _guide1->_parameters.reflect[0]='\0';
  #define reflect (_guide1->_parameters.reflect)
  _guide1->_parameters.w1 = 0.05;
  #define w1 (_guide1->_parameters.w1)
  _guide1->_parameters.h1 = 0.12;
  #define h1 (_guide1->_parameters.h1)
  _guide1->_parameters.w2 = 0.05;
  #define w2 (_guide1->_parameters.w2)
  _guide1->_parameters.h2 = 0.12;
  #define h2 (_guide1->_parameters.h2)
  _guide1->_parameters.l = 4.66;
  #define l (_guide1->_parameters.l)
  _guide1->_parameters.R0 = 0.995;
  #define R0 (_guide1->_parameters.R0)
  _guide1->_parameters.Qc = 0.0217;
  #define Qc (_guide1->_parameters.Qc)
  _guide1->_parameters.alpha = 5.76;
  #define alpha (_guide1->_parameters.alpha)
  _guide1->_parameters.m = 2.0;
  #define m (_guide1->_parameters.m)
  _guide1->_parameters.W = 0.0033;
  #define W (_guide1->_parameters.W)

  #define pTable (_guide1->_parameters.pTable)

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
  /* component guide1=Guide() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (0)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _a1->_rotation_absolute, _guide1->_rotation_absolute);
    rot_transpose(_csource->_rotation_absolute, tr1);
    rot_mul(_guide1->_rotation_absolute, tr1, _guide1->_rotation_relative);
    _guide1->_rotation_is_identity =  rot_test_identity(_guide1->_rotation_relative);
    tc1 = coords_set(
      0, 0, 1.5);
    rot_transpose(_a1->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _guide1->_position_absolute = coords_add(_a1->_position_absolute, tc2);
    tc1 = coords_sub(_csource->_position_absolute, _guide1->_position_absolute);
    _guide1->_position_relative = rot_apply(_guide1->_rotation_absolute, tc1);
  } /* guide1=Guide() AT ROTATED */
  DEBUG_COMPONENT("guide1", _guide1->_position_absolute, _guide1->_rotation_absolute);
  instrument->_position_absolute[3] = _guide1->_position_absolute;
  instrument->_position_relative[3] = _guide1->_position_relative;
  instrument->counter_N[3]  = instrument->counter_P[3] = instrument->counter_P2[3] = 0;
  instrument->counter_AbsorbProp[3]= 0;
  return(0);
} /* _guide1_setpos */

/* component guide2=Bender() SETTING, POSITION/ROTATION */
int _guide2_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_guide2_setpos] component guide2=Bender() SETTING [/usr/share/mcstas/3.0-dev/optics/Bender.comp:112]");
  stracpy(_guide2->_name, "guide2", 16384);
  stracpy(_guide2->_type, "Bender", 16384);
  _guide2->_index=4;
  _guide2->_parameters.w = 0.05;
  #define w (_guide2->_parameters.w)
  _guide2->_parameters.h = 0.12;
  #define h (_guide2->_parameters.h)
  _guide2->_parameters.r = 1445;
  #define r (_guide2->_parameters.r)
  _guide2->_parameters.Win = 0.04;
  #define Win (_guide2->_parameters.Win)
  _guide2->_parameters.k = 1;
  #define k (_guide2->_parameters.k)
  _guide2->_parameters.d = 0.001;
  #define d (_guide2->_parameters.d)
  _guide2->_parameters.l = 24.5;
  #define l (_guide2->_parameters.l)
  _guide2->_parameters.R0a = 0.995;
  #define R0a (_guide2->_parameters.R0a)
  _guide2->_parameters.Qca = 0.0217;
  #define Qca (_guide2->_parameters.Qca)
  _guide2->_parameters.alphaa = 5.76;
  #define alphaa (_guide2->_parameters.alphaa)
  _guide2->_parameters.ma = 2;
  #define ma (_guide2->_parameters.ma)
  _guide2->_parameters.Wa = 0.0033;
  #define Wa (_guide2->_parameters.Wa)
  _guide2->_parameters.R0i = 0.995;
  #define R0i (_guide2->_parameters.R0i)
  _guide2->_parameters.Qci = 0.0217;
  #define Qci (_guide2->_parameters.Qci)
  _guide2->_parameters.alphai = 5.76;
  #define alphai (_guide2->_parameters.alphai)
  _guide2->_parameters.mi = 2;
  #define mi (_guide2->_parameters.mi)
  _guide2->_parameters.Wi = 0.0033;
  #define Wi (_guide2->_parameters.Wi)
  _guide2->_parameters.R0s = 0.995;
  #define R0s (_guide2->_parameters.R0s)
  _guide2->_parameters.Qcs = 0.0217;
  #define Qcs (_guide2->_parameters.Qcs)
  _guide2->_parameters.alphas = 5.76;
  #define alphas (_guide2->_parameters.alphas)
  _guide2->_parameters.ms = 2;
  #define ms (_guide2->_parameters.ms)
  _guide2->_parameters.Ws = 0.0033;
  #define Ws (_guide2->_parameters.Ws)

  #define bk (_guide2->_parameters.bk)
  #define mWin (_guide2->_parameters.mWin)

  #undef w
  #undef h
  #undef r
  #undef Win
  #undef k
  #undef d
  #undef l
  #undef R0a
  #undef Qca
  #undef alphaa
  #undef ma
  #undef Wa
  #undef R0i
  #undef Qci
  #undef alphai
  #undef mi
  #undef Wi
  #undef R0s
  #undef Qcs
  #undef alphas
  #undef ms
  #undef Ws
  #undef bk
  #undef mWin
  /* component guide2=Bender() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (0)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _a1->_rotation_absolute, _guide2->_rotation_absolute);
    rot_transpose(_guide1->_rotation_absolute, tr1);
    rot_mul(_guide2->_rotation_absolute, tr1, _guide2->_rotation_relative);
    _guide2->_rotation_is_identity =  rot_test_identity(_guide2->_rotation_relative);
    tc1 = coords_set(
      0, 0, 6.16);
    rot_transpose(_a1->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _guide2->_position_absolute = coords_add(_a1->_position_absolute, tc2);
    tc1 = coords_sub(_guide1->_position_absolute, _guide2->_position_absolute);
    _guide2->_position_relative = rot_apply(_guide2->_rotation_absolute, tc1);
  } /* guide2=Bender() AT ROTATED */
  DEBUG_COMPONENT("guide2", _guide2->_position_absolute, _guide2->_rotation_absolute);
  instrument->_position_absolute[4] = _guide2->_position_absolute;
  instrument->_position_relative[4] = _guide2->_position_relative;
  instrument->counter_N[4]  = instrument->counter_P[4] = instrument->counter_P2[4] = 0;
  instrument->counter_AbsorbProp[4]= 0;
  return(0);
} /* _guide2_setpos */

/* component bunker=Guide() SETTING, POSITION/ROTATION */
int _bunker_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_bunker_setpos] component bunker=Guide() SETTING [/usr/share/mcstas/3.0-dev/optics/Guide.comp:74]");
  stracpy(_bunker->_name, "bunker", 16384);
  stracpy(_bunker->_type, "Guide", 16384);
  _bunker->_index=5;
  _bunker->_parameters.reflect[0]='\0';
  #define reflect (_bunker->_parameters.reflect)
  _bunker->_parameters.w1 = 0.05;
  #define w1 (_bunker->_parameters.w1)
  _bunker->_parameters.h1 = .12;
  #define h1 (_bunker->_parameters.h1)
  _bunker->_parameters.w2 = 0.05;
  #define w2 (_bunker->_parameters.w2)
  _bunker->_parameters.h2 = .12;
  #define h2 (_bunker->_parameters.h2)
  _bunker->_parameters.l = 3.0;
  #define l (_bunker->_parameters.l)
  _bunker->_parameters.R0 = 0.995;
  #define R0 (_bunker->_parameters.R0)
  _bunker->_parameters.Qc = 0.0217;
  #define Qc (_bunker->_parameters.Qc)
  _bunker->_parameters.alpha = 5.76;
  #define alpha (_bunker->_parameters.alpha)
  _bunker->_parameters.m = 2.0;
  #define m (_bunker->_parameters.m)
  _bunker->_parameters.W = 0.0033;
  #define W (_bunker->_parameters.W)

  #define pTable (_bunker->_parameters.pTable)

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
  /* component bunker=Guide() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (0)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _a1->_rotation_absolute, _bunker->_rotation_absolute);
    rot_transpose(_guide2->_rotation_absolute, tr1);
    rot_mul(_bunker->_rotation_absolute, tr1, _bunker->_rotation_relative);
    _bunker->_rotation_is_identity =  rot_test_identity(_bunker->_rotation_relative);
    tc1 = coords_set(
      0, 0, 30.66);
    rot_transpose(_a1->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _bunker->_position_absolute = coords_add(_a1->_position_absolute, tc2);
    tc1 = coords_sub(_guide2->_position_absolute, _bunker->_position_absolute);
    _bunker->_position_relative = rot_apply(_bunker->_rotation_absolute, tc1);
  } /* bunker=Guide() AT ROTATED */
  DEBUG_COMPONENT("bunker", _bunker->_position_absolute, _bunker->_rotation_absolute);
  instrument->_position_absolute[5] = _bunker->_position_absolute;
  instrument->_position_relative[5] = _bunker->_position_relative;
  instrument->counter_N[5]  = instrument->counter_P[5] = instrument->counter_P2[5] = 0;
  instrument->counter_AbsorbProp[5]= 0;
  return(0);
} /* _bunker_setpos */

/* component guide3=Guide() SETTING, POSITION/ROTATION */
int _guide3_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_guide3_setpos] component guide3=Guide() SETTING [/usr/share/mcstas/3.0-dev/optics/Guide.comp:74]");
  stracpy(_guide3->_name, "guide3", 16384);
  stracpy(_guide3->_type, "Guide", 16384);
  _guide3->_index=6;
  _guide3->_parameters.reflect[0]='\0';
  #define reflect (_guide3->_parameters.reflect)
  _guide3->_parameters.w1 = 0.05;
  #define w1 (_guide3->_parameters.w1)
  _guide3->_parameters.h1 = .12;
  #define h1 (_guide3->_parameters.h1)
  _guide3->_parameters.w2 = 0.05;
  #define w2 (_guide3->_parameters.w2)
  _guide3->_parameters.h2 = .12;
  #define h2 (_guide3->_parameters.h2)
  _guide3->_parameters.l = 32.95;
  #define l (_guide3->_parameters.l)
  _guide3->_parameters.R0 = 0.995;
  #define R0 (_guide3->_parameters.R0)
  _guide3->_parameters.Qc = 0.0217;
  #define Qc (_guide3->_parameters.Qc)
  _guide3->_parameters.alpha = 5.76;
  #define alpha (_guide3->_parameters.alpha)
  _guide3->_parameters.m = 2.0;
  #define m (_guide3->_parameters.m)
  _guide3->_parameters.W = 0.0033;
  #define W (_guide3->_parameters.W)

  #define pTable (_guide3->_parameters.pTable)

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
  /* component guide3=Guide() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (0)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _a1->_rotation_absolute, _guide3->_rotation_absolute);
    rot_transpose(_bunker->_rotation_absolute, tr1);
    rot_mul(_guide3->_rotation_absolute, tr1, _guide3->_rotation_relative);
    _guide3->_rotation_is_identity =  rot_test_identity(_guide3->_rotation_relative);
    tc1 = coords_set(
      0, 0, 33.66);
    rot_transpose(_a1->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _guide3->_position_absolute = coords_add(_a1->_position_absolute, tc2);
    tc1 = coords_sub(_bunker->_position_absolute, _guide3->_position_absolute);
    _guide3->_position_relative = rot_apply(_guide3->_rotation_absolute, tc1);
  } /* guide3=Guide() AT ROTATED */
  DEBUG_COMPONENT("guide3", _guide3->_position_absolute, _guide3->_rotation_absolute);
  instrument->_position_absolute[6] = _guide3->_position_absolute;
  instrument->_position_relative[6] = _guide3->_position_relative;
  instrument->counter_N[6]  = instrument->counter_P[6] = instrument->counter_P2[6] = 0;
  instrument->counter_AbsorbProp[6]= 0;
  return(0);
} /* _guide3_setpos */

/* component lambdaGuideExit=L_monitor() SETTING, POSITION/ROTATION */
int _lambdaGuideExit_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_lambdaGuideExit_setpos] component lambdaGuideExit=L_monitor() SETTING [/usr/share/mcstas/3.0-dev/monitors/L_monitor.comp:65]");
  stracpy(_lambdaGuideExit->_name, "lambdaGuideExit", 16384);
  stracpy(_lambdaGuideExit->_type, "L_monitor", 16384);
  _lambdaGuideExit->_index=7;
  _lambdaGuideExit->_parameters.nL = 60;
  #define nL (_lambdaGuideExit->_parameters.nL)
  if("lambdaguide.dat" && strlen("lambdaguide.dat"))
    stracpy(_lambdaGuideExit->_parameters.filename, "lambdaguide.dat" ? "lambdaguide.dat" : "", 16384);
  else 
  _lambdaGuideExit->_parameters.filename[0]='\0';
  #define filename (_lambdaGuideExit->_parameters.filename)
  _lambdaGuideExit->_parameters.xmin = -0.05;
  #define xmin (_lambdaGuideExit->_parameters.xmin)
  _lambdaGuideExit->_parameters.xmax = 0.05;
  #define xmax (_lambdaGuideExit->_parameters.xmax)
  _lambdaGuideExit->_parameters.ymin = -0.05;
  #define ymin (_lambdaGuideExit->_parameters.ymin)
  _lambdaGuideExit->_parameters.ymax = 0.05;
  #define ymax (_lambdaGuideExit->_parameters.ymax)
  _lambdaGuideExit->_parameters.xwidth = 0.05;
  #define xwidth (_lambdaGuideExit->_parameters.xwidth)
  _lambdaGuideExit->_parameters.yheight = 0.12;
  #define yheight (_lambdaGuideExit->_parameters.yheight)
  _lambdaGuideExit->_parameters.Lmin = LMIN;
  #define Lmin (_lambdaGuideExit->_parameters.Lmin)
  _lambdaGuideExit->_parameters.Lmax = LMAX;
  #define Lmax (_lambdaGuideExit->_parameters.Lmax)
  _lambdaGuideExit->_parameters.restore_neutron = 0;
  #define restore_neutron (_lambdaGuideExit->_parameters.restore_neutron)

  #define L_N (_lambdaGuideExit->_parameters.L_N)
  #define L_p (_lambdaGuideExit->_parameters.L_p)
  #define L_p2 (_lambdaGuideExit->_parameters.L_p2)

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
  /* component lambdaGuideExit=L_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (0)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _a1->_rotation_absolute, _lambdaGuideExit->_rotation_absolute);
    rot_transpose(_guide3->_rotation_absolute, tr1);
    rot_mul(_lambdaGuideExit->_rotation_absolute, tr1, _lambdaGuideExit->_rotation_relative);
    _lambdaGuideExit->_rotation_is_identity =  rot_test_identity(_lambdaGuideExit->_rotation_relative);
    tc1 = coords_set(
      0, 0, 66.611);
    rot_transpose(_a1->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _lambdaGuideExit->_position_absolute = coords_add(_a1->_position_absolute, tc2);
    tc1 = coords_sub(_guide3->_position_absolute, _lambdaGuideExit->_position_absolute);
    _lambdaGuideExit->_position_relative = rot_apply(_lambdaGuideExit->_rotation_absolute, tc1);
  } /* lambdaGuideExit=L_monitor() AT ROTATED */
  DEBUG_COMPONENT("lambdaGuideExit", _lambdaGuideExit->_position_absolute, _lambdaGuideExit->_rotation_absolute);
  instrument->_position_absolute[7] = _lambdaGuideExit->_position_absolute;
  instrument->_position_relative[7] = _lambdaGuideExit->_position_relative;
  instrument->counter_N[7]  = instrument->counter_P[7] = instrument->counter_P2[7] = 0;
  instrument->counter_AbsorbProp[7]= 0;
  return(0);
} /* _lambdaGuideExit_setpos */

/* component DivMonGuideExit=Divergence_monitor() SETTING, POSITION/ROTATION */
int _DivMonGuideExit_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_DivMonGuideExit_setpos] component DivMonGuideExit=Divergence_monitor() SETTING [/usr/share/mcstas/3.0-dev/monitors/Divergence_monitor.comp:71]");
  stracpy(_DivMonGuideExit->_name, "DivMonGuideExit", 16384);
  stracpy(_DivMonGuideExit->_type, "Divergence_monitor", 16384);
  _DivMonGuideExit->_index=8;
  _DivMonGuideExit->_parameters.nh = 60;
  #define nh (_DivMonGuideExit->_parameters.nh)
  _DivMonGuideExit->_parameters.nv = 60;
  #define nv (_DivMonGuideExit->_parameters.nv)
  if("divguide.dat" && strlen("divguide.dat"))
    stracpy(_DivMonGuideExit->_parameters.filename, "divguide.dat" ? "divguide.dat" : "", 16384);
  else 
  _DivMonGuideExit->_parameters.filename[0]='\0';
  #define filename (_DivMonGuideExit->_parameters.filename)
  _DivMonGuideExit->_parameters.xmin = -0.05;
  #define xmin (_DivMonGuideExit->_parameters.xmin)
  _DivMonGuideExit->_parameters.xmax = 0.05;
  #define xmax (_DivMonGuideExit->_parameters.xmax)
  _DivMonGuideExit->_parameters.ymin = -0.05;
  #define ymin (_DivMonGuideExit->_parameters.ymin)
  _DivMonGuideExit->_parameters.ymax = 0.05;
  #define ymax (_DivMonGuideExit->_parameters.ymax)
  _DivMonGuideExit->_parameters.xwidth = 0.05;
  #define xwidth (_DivMonGuideExit->_parameters.xwidth)
  _DivMonGuideExit->_parameters.yheight = 0.12;
  #define yheight (_DivMonGuideExit->_parameters.yheight)
  _DivMonGuideExit->_parameters.maxdiv_h = 2;
  #define maxdiv_h (_DivMonGuideExit->_parameters.maxdiv_h)
  _DivMonGuideExit->_parameters.maxdiv_v = 2;
  #define maxdiv_v (_DivMonGuideExit->_parameters.maxdiv_v)
  _DivMonGuideExit->_parameters.restore_neutron = 0;
  #define restore_neutron (_DivMonGuideExit->_parameters.restore_neutron)
  _DivMonGuideExit->_parameters.nx = 0;
  #define nx (_DivMonGuideExit->_parameters.nx)
  _DivMonGuideExit->_parameters.ny = 0;
  #define ny (_DivMonGuideExit->_parameters.ny)
  _DivMonGuideExit->_parameters.nz = 1;
  #define nz (_DivMonGuideExit->_parameters.nz)

  #define Div_N (_DivMonGuideExit->_parameters.Div_N)
  #define Div_p (_DivMonGuideExit->_parameters.Div_p)
  #define Div_p2 (_DivMonGuideExit->_parameters.Div_p2)

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
  /* component DivMonGuideExit=Divergence_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (0)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _a1->_rotation_absolute, _DivMonGuideExit->_rotation_absolute);
    rot_transpose(_lambdaGuideExit->_rotation_absolute, tr1);
    rot_mul(_DivMonGuideExit->_rotation_absolute, tr1, _DivMonGuideExit->_rotation_relative);
    _DivMonGuideExit->_rotation_is_identity =  rot_test_identity(_DivMonGuideExit->_rotation_relative);
    tc1 = coords_set(
      0, 0, 66.612);
    rot_transpose(_a1->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _DivMonGuideExit->_position_absolute = coords_add(_a1->_position_absolute, tc2);
    tc1 = coords_sub(_lambdaGuideExit->_position_absolute, _DivMonGuideExit->_position_absolute);
    _DivMonGuideExit->_position_relative = rot_apply(_DivMonGuideExit->_rotation_absolute, tc1);
  } /* DivMonGuideExit=Divergence_monitor() AT ROTATED */
  DEBUG_COMPONENT("DivMonGuideExit", _DivMonGuideExit->_position_absolute, _DivMonGuideExit->_rotation_absolute);
  instrument->_position_absolute[8] = _DivMonGuideExit->_position_absolute;
  instrument->_position_relative[8] = _DivMonGuideExit->_position_relative;
  instrument->counter_N[8]  = instrument->counter_P[8] = instrument->counter_P2[8] = 0;
  instrument->counter_AbsorbProp[8]= 0;
  return(0);
} /* _DivMonGuideExit_setpos */

/* component PSDGuideExit=PSD_monitor() SETTING, POSITION/ROTATION */
int _PSDGuideExit_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_PSDGuideExit_setpos] component PSDGuideExit=PSD_monitor() SETTING [/usr/share/mcstas/3.0-dev/monitors/PSD_monitor.comp:64]");
  stracpy(_PSDGuideExit->_name, "PSDGuideExit", 16384);
  stracpy(_PSDGuideExit->_type, "PSD_monitor", 16384);
  _PSDGuideExit->_index=9;
  _PSDGuideExit->_parameters.nx = 60;
  #define nx (_PSDGuideExit->_parameters.nx)
  _PSDGuideExit->_parameters.ny = 60;
  #define ny (_PSDGuideExit->_parameters.ny)
  if("psdguide.dat" && strlen("psdguide.dat"))
    stracpy(_PSDGuideExit->_parameters.filename, "psdguide.dat" ? "psdguide.dat" : "", 16384);
  else 
  _PSDGuideExit->_parameters.filename[0]='\0';
  #define filename (_PSDGuideExit->_parameters.filename)
  _PSDGuideExit->_parameters.xmin = -0.05;
  #define xmin (_PSDGuideExit->_parameters.xmin)
  _PSDGuideExit->_parameters.xmax = 0.05;
  #define xmax (_PSDGuideExit->_parameters.xmax)
  _PSDGuideExit->_parameters.ymin = -0.05;
  #define ymin (_PSDGuideExit->_parameters.ymin)
  _PSDGuideExit->_parameters.ymax = 0.05;
  #define ymax (_PSDGuideExit->_parameters.ymax)
  _PSDGuideExit->_parameters.xwidth = 0.05;
  #define xwidth (_PSDGuideExit->_parameters.xwidth)
  _PSDGuideExit->_parameters.yheight = 0.12;
  #define yheight (_PSDGuideExit->_parameters.yheight)
  _PSDGuideExit->_parameters.restore_neutron = 0;
  #define restore_neutron (_PSDGuideExit->_parameters.restore_neutron)

  #define PSD_N (_PSDGuideExit->_parameters.PSD_N)
  #define PSD_p (_PSDGuideExit->_parameters.PSD_p)
  #define PSD_p2 (_PSDGuideExit->_parameters.PSD_p2)

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
  /* component PSDGuideExit=PSD_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (0)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _a1->_rotation_absolute, _PSDGuideExit->_rotation_absolute);
    rot_transpose(_DivMonGuideExit->_rotation_absolute, tr1);
    rot_mul(_PSDGuideExit->_rotation_absolute, tr1, _PSDGuideExit->_rotation_relative);
    _PSDGuideExit->_rotation_is_identity =  rot_test_identity(_PSDGuideExit->_rotation_relative);
    tc1 = coords_set(
      0, 0, 66.613);
    rot_transpose(_a1->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _PSDGuideExit->_position_absolute = coords_add(_a1->_position_absolute, tc2);
    tc1 = coords_sub(_DivMonGuideExit->_position_absolute, _PSDGuideExit->_position_absolute);
    _PSDGuideExit->_position_relative = rot_apply(_PSDGuideExit->_rotation_absolute, tc1);
  } /* PSDGuideExit=PSD_monitor() AT ROTATED */
  DEBUG_COMPONENT("PSDGuideExit", _PSDGuideExit->_position_absolute, _PSDGuideExit->_rotation_absolute);
  instrument->_position_absolute[9] = _PSDGuideExit->_position_absolute;
  instrument->_position_relative[9] = _PSDGuideExit->_position_relative;
  instrument->counter_N[9]  = instrument->counter_P[9] = instrument->counter_P2[9] = 0;
  instrument->counter_AbsorbProp[9]= 0;
  return(0);
} /* _PSDGuideExit_setpos */

/* component FOCUSguide=Guide_channeled() SETTING, POSITION/ROTATION */
int _FOCUSguide_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_FOCUSguide_setpos] component FOCUSguide=Guide_channeled() SETTING [/usr/share/mcstas/3.0-dev/optics/Guide_channeled.comp:89]");
  stracpy(_FOCUSguide->_name, "FOCUSguide", 16384);
  stracpy(_FOCUSguide->_type, "Guide_channeled", 16384);
  _FOCUSguide->_index=10;
  _FOCUSguide->_parameters.w1 = 0.05;
  #define w1 (_FOCUSguide->_parameters.w1)
  _FOCUSguide->_parameters.h1 = 0.12;
  #define h1 (_FOCUSguide->_parameters.h1)
  _FOCUSguide->_parameters.w2 = 0.05;
  #define w2 (_FOCUSguide->_parameters.w2)
  _FOCUSguide->_parameters.h2 = 0.095;
  #define h2 (_FOCUSguide->_parameters.h2)
  _FOCUSguide->_parameters.l = 2.95;
  #define l (_FOCUSguide->_parameters.l)
  _FOCUSguide->_parameters.R0 = 0.995;
  #define R0 (_FOCUSguide->_parameters.R0)
  _FOCUSguide->_parameters.Qc = 0;
  #define Qc (_FOCUSguide->_parameters.Qc)
  _FOCUSguide->_parameters.alpha = 0;
  #define alpha (_FOCUSguide->_parameters.alpha)
  _FOCUSguide->_parameters.m = 0;
  #define m (_FOCUSguide->_parameters.m)
  _FOCUSguide->_parameters.nslit = 1;
  #define nslit (_FOCUSguide->_parameters.nslit)
  _FOCUSguide->_parameters.d = 0.002;
  #define d (_FOCUSguide->_parameters.d)
  _FOCUSguide->_parameters.Qcx = 0.0217;
  #define Qcx (_FOCUSguide->_parameters.Qcx)
  _FOCUSguide->_parameters.Qcy = 0.0217;
  #define Qcy (_FOCUSguide->_parameters.Qcy)
  _FOCUSguide->_parameters.alphax = 5.76;
  #define alphax (_FOCUSguide->_parameters.alphax)
  _FOCUSguide->_parameters.alphay = 5.64;
  #define alphay (_FOCUSguide->_parameters.alphay)
  _FOCUSguide->_parameters.W = 0.0033;
  #define W (_FOCUSguide->_parameters.W)
  _FOCUSguide->_parameters.mx = 2.4;
  #define mx (_FOCUSguide->_parameters.mx)
  _FOCUSguide->_parameters.my = 3;
  #define my (_FOCUSguide->_parameters.my)
  _FOCUSguide->_parameters.nu = 0;
  #define nu (_FOCUSguide->_parameters.nu)
  _FOCUSguide->_parameters.phase = 0;
  #define phase (_FOCUSguide->_parameters.phase)

  #define w1c (_FOCUSguide->_parameters.w1c)
  #define w2c (_FOCUSguide->_parameters.w2c)
  #define ww (_FOCUSguide->_parameters.ww)
  #define hh (_FOCUSguide->_parameters.hh)
  #define whalf (_FOCUSguide->_parameters.whalf)
  #define hhalf (_FOCUSguide->_parameters.hhalf)
  #define lwhalf (_FOCUSguide->_parameters.lwhalf)
  #define lhhalf (_FOCUSguide->_parameters.lhhalf)

  #undef w1
  #undef h1
  #undef w2
  #undef h2
  #undef l
  #undef R0
  #undef Qc
  #undef alpha
  #undef m
  #undef nslit
  #undef d
  #undef Qcx
  #undef Qcy
  #undef alphax
  #undef alphay
  #undef W
  #undef mx
  #undef my
  #undef nu
  #undef phase
  #undef w1c
  #undef w2c
  #undef ww
  #undef hh
  #undef whalf
  #undef hhalf
  #undef lwhalf
  #undef lhhalf
  /* component FOCUSguide=Guide_channeled() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (0)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _a1->_rotation_absolute, _FOCUSguide->_rotation_absolute);
    rot_transpose(_PSDGuideExit->_rotation_absolute, tr1);
    rot_mul(_FOCUSguide->_rotation_absolute, tr1, _FOCUSguide->_rotation_relative);
    _FOCUSguide->_rotation_is_identity =  rot_test_identity(_FOCUSguide->_rotation_relative);
    tc1 = coords_set(
      0, 0, 66.67);
    rot_transpose(_a1->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _FOCUSguide->_position_absolute = coords_add(_a1->_position_absolute, tc2);
    tc1 = coords_sub(_PSDGuideExit->_position_absolute, _FOCUSguide->_position_absolute);
    _FOCUSguide->_position_relative = rot_apply(_FOCUSguide->_rotation_absolute, tc1);
  } /* FOCUSguide=Guide_channeled() AT ROTATED */
  DEBUG_COMPONENT("FOCUSguide", _FOCUSguide->_position_absolute, _FOCUSguide->_rotation_absolute);
  instrument->_position_absolute[10] = _FOCUSguide->_position_absolute;
  instrument->_position_relative[10] = _FOCUSguide->_position_relative;
  instrument->counter_N[10]  = instrument->counter_P[10] = instrument->counter_P2[10] = 0;
  instrument->counter_AbsorbProp[10]= 0;
  return(0);
} /* _FOCUSguide_setpos */

/* component FirstChopper=DiskChopper() SETTING, POSITION/ROTATION */
int _FirstChopper_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_FirstChopper_setpos] component FirstChopper=DiskChopper() SETTING [/usr/share/mcstas/3.0-dev/optics/DiskChopper.comp:67]");
  stracpy(_FirstChopper->_name, "FirstChopper", 16384);
  stracpy(_FirstChopper->_type, "DiskChopper", 16384);
  _FirstChopper->_index=11;
  _FirstChopper->_parameters.theta_0 = 0;
  #define theta_0 (_FirstChopper->_parameters.theta_0)
  _FirstChopper->_parameters.radius = 0.70;
  #define radius (_FirstChopper->_parameters.radius)
  _FirstChopper->_parameters.yheight = 0.70;
  #define yheight (_FirstChopper->_parameters.yheight)
  _FirstChopper->_parameters.nu = DISC_SPEED / 2 / PI;
  #define nu (_FirstChopper->_parameters.nu)
  _FirstChopper->_parameters.nslit = 2;
  #define nslit (_FirstChopper->_parameters.nslit)
  _FirstChopper->_parameters.jitter = 0;
  #define jitter (_FirstChopper->_parameters.jitter)
  _FirstChopper->_parameters.delay = FO_PHA;
  #define delay (_FirstChopper->_parameters.delay)
  _FirstChopper->_parameters.isfirst = 1;
  #define isfirst (_FirstChopper->_parameters.isfirst)
  _FirstChopper->_parameters.n_pulse = 1;
  #define n_pulse (_FirstChopper->_parameters.n_pulse)
  _FirstChopper->_parameters.abs_out = 0;
  #define abs_out (_FirstChopper->_parameters.abs_out)
  _FirstChopper->_parameters.phase = 0;
  #define phase (_FirstChopper->_parameters.phase)
  _FirstChopper->_parameters.xwidth = 0.07;
  #define xwidth (_FirstChopper->_parameters.xwidth)
  _FirstChopper->_parameters.verbose = 1;
  #define verbose (_FirstChopper->_parameters.verbose)

  #define Tg (_FirstChopper->_parameters.Tg)
  #define To (_FirstChopper->_parameters.To)
  #define delta_y (_FirstChopper->_parameters.delta_y)
  #define height (_FirstChopper->_parameters.height)
  #define omega (_FirstChopper->_parameters.omega)

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
  /* component FirstChopper=DiskChopper() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _a1->_rotation_absolute, _FirstChopper->_rotation_absolute);
    rot_transpose(_FOCUSguide->_rotation_absolute, tr1);
    rot_mul(_FirstChopper->_rotation_absolute, tr1, _FirstChopper->_rotation_relative);
    _FirstChopper->_rotation_is_identity =  rot_test_identity(_FirstChopper->_rotation_relative);
    tc1 = coords_set(
      0, 0, 69.674);
    rot_transpose(_a1->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _FirstChopper->_position_absolute = coords_add(_a1->_position_absolute, tc2);
    tc1 = coords_sub(_FOCUSguide->_position_absolute, _FirstChopper->_position_absolute);
    _FirstChopper->_position_relative = rot_apply(_FirstChopper->_rotation_absolute, tc1);
  } /* FirstChopper=DiskChopper() AT ROTATED */
  DEBUG_COMPONENT("FirstChopper", _FirstChopper->_position_absolute, _FirstChopper->_rotation_absolute);
  instrument->_position_absolute[11] = _FirstChopper->_position_absolute;
  instrument->_position_relative[11] = _FirstChopper->_position_relative;
  instrument->counter_N[11]  = instrument->counter_P[11] = instrument->counter_P2[11] = 0;
  instrument->counter_AbsorbProp[11]= 0;
  return(0);
} /* _FirstChopper_setpos */

/* component DISCTOF=Monitor_nD() SETTING, POSITION/ROTATION */
int _DISCTOF_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_DISCTOF_setpos] component DISCTOF=Monitor_nD() SETTING [/usr/share/mcstas/3.0-dev/monitors/Monitor_nD.comp:228]");
  stracpy(_DISCTOF->_name, "DISCTOF", 16384);
  stracpy(_DISCTOF->_type, "Monitor_nD", 16384);
  _DISCTOF->_index=12;
  _DISCTOF->_parameters.user1 = '\0';
  #define user1 (_DISCTOF->_parameters.user1)
  _DISCTOF->_parameters.user2 = '\0';
  #define user2 (_DISCTOF->_parameters.user2)
  _DISCTOF->_parameters.user3 = '\0';
  #define user3 (_DISCTOF->_parameters.user3)
  _DISCTOF->_parameters.xwidth = 0;
  #define xwidth (_DISCTOF->_parameters.xwidth)
  _DISCTOF->_parameters.yheight = 0;
  #define yheight (_DISCTOF->_parameters.yheight)
  _DISCTOF->_parameters.zdepth = 0;
  #define zdepth (_DISCTOF->_parameters.zdepth)
  _DISCTOF->_parameters.xmin = -0.05;
  #define xmin (_DISCTOF->_parameters.xmin)
  _DISCTOF->_parameters.xmax = 0.05;
  #define xmax (_DISCTOF->_parameters.xmax)
  _DISCTOF->_parameters.ymin = -0.1;
  #define ymin (_DISCTOF->_parameters.ymin)
  _DISCTOF->_parameters.ymax = 0.1;
  #define ymax (_DISCTOF->_parameters.ymax)
  _DISCTOF->_parameters.zmin = 0;
  #define zmin (_DISCTOF->_parameters.zmin)
  _DISCTOF->_parameters.zmax = 0;
  #define zmax (_DISCTOF->_parameters.zmax)
  _DISCTOF->_parameters.bins = 30;
  #define bins (_DISCTOF->_parameters.bins)
  _DISCTOF->_parameters.min = -1e40;
  #define min (_DISCTOF->_parameters.min)
  _DISCTOF->_parameters.max = 1e40;
  #define max (_DISCTOF->_parameters.max)
  _DISCTOF->_parameters.restore_neutron = 0;
  #define restore_neutron (_DISCTOF->_parameters.restore_neutron)
  _DISCTOF->_parameters.radius = 0;
  #define radius (_DISCTOF->_parameters.radius)
  if("auto, time" && strlen("auto, time"))
    stracpy(_DISCTOF->_parameters.options, "auto, time" ? "auto, time" : "", 16384);
  else 
  _DISCTOF->_parameters.options[0]='\0';
  #define options (_DISCTOF->_parameters.options)
  if("DISC_TOF.dat" && strlen("DISC_TOF.dat"))
    stracpy(_DISCTOF->_parameters.filename, "DISC_TOF.dat" ? "DISC_TOF.dat" : "", 16384);
  else 
  _DISCTOF->_parameters.filename[0]='\0';
  #define filename (_DISCTOF->_parameters.filename)
  if("NULL" && strlen("NULL"))
    stracpy(_DISCTOF->_parameters.geometry, "NULL" ? "NULL" : "", 16384);
  else 
  _DISCTOF->_parameters.geometry[0]='\0';
  #define geometry (_DISCTOF->_parameters.geometry)
  if("NULL" && strlen("NULL"))
    stracpy(_DISCTOF->_parameters.username1, "NULL" ? "NULL" : "", 16384);
  else 
  _DISCTOF->_parameters.username1[0]='\0';
  #define username1 (_DISCTOF->_parameters.username1)
  if("NULL" && strlen("NULL"))
    stracpy(_DISCTOF->_parameters.username2, "NULL" ? "NULL" : "", 16384);
  else 
  _DISCTOF->_parameters.username2[0]='\0';
  #define username2 (_DISCTOF->_parameters.username2)
  if("NULL" && strlen("NULL"))
    stracpy(_DISCTOF->_parameters.username3, "NULL" ? "NULL" : "", 16384);
  else 
  _DISCTOF->_parameters.username3[0]='\0';
  #define username3 (_DISCTOF->_parameters.username3)

  #define DEFS (_DISCTOF->_parameters.DEFS)
  #define Vars (_DISCTOF->_parameters.Vars)
  #define detector (_DISCTOF->_parameters.detector)
  #define offdata (_DISCTOF->_parameters.offdata)

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
  /* component DISCTOF=Monitor_nD() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _a1->_rotation_absolute, _DISCTOF->_rotation_absolute);
    rot_transpose(_FirstChopper->_rotation_absolute, tr1);
    rot_mul(_DISCTOF->_rotation_absolute, tr1, _DISCTOF->_rotation_relative);
    _DISCTOF->_rotation_is_identity =  rot_test_identity(_DISCTOF->_rotation_relative);
    tc1 = coords_set(
      0, 0, 69.718);
    rot_transpose(_a1->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _DISCTOF->_position_absolute = coords_add(_a1->_position_absolute, tc2);
    tc1 = coords_sub(_FirstChopper->_position_absolute, _DISCTOF->_position_absolute);
    _DISCTOF->_position_relative = rot_apply(_DISCTOF->_rotation_absolute, tc1);
  } /* DISCTOF=Monitor_nD() AT ROTATED */
  DEBUG_COMPONENT("DISCTOF", _DISCTOF->_position_absolute, _DISCTOF->_rotation_absolute);
  instrument->_position_absolute[12] = _DISCTOF->_position_absolute;
  instrument->_position_relative[12] = _DISCTOF->_position_relative;
  instrument->counter_N[12]  = instrument->counter_P[12] = instrument->counter_P2[12] = 0;
  instrument->counter_AbsorbProp[12]= 0;
  return(0);
} /* _DISCTOF_setpos */

/* component PSDmon1Chopper=PSD_monitor() SETTING, POSITION/ROTATION */
int _PSDmon1Chopper_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_PSDmon1Chopper_setpos] component PSDmon1Chopper=PSD_monitor() SETTING [/usr/share/mcstas/3.0-dev/monitors/PSD_monitor.comp:64]");
  stracpy(_PSDmon1Chopper->_name, "PSDmon1Chopper", 16384);
  stracpy(_PSDmon1Chopper->_type, "PSD_monitor", 16384);
  _PSDmon1Chopper->_index=13;
  _PSDmon1Chopper->_parameters.nx = 60;
  #define nx (_PSDmon1Chopper->_parameters.nx)
  _PSDmon1Chopper->_parameters.ny = 60;
  #define ny (_PSDmon1Chopper->_parameters.ny)
  if("psdchopper.dat" && strlen("psdchopper.dat"))
    stracpy(_PSDmon1Chopper->_parameters.filename, "psdchopper.dat" ? "psdchopper.dat" : "", 16384);
  else 
  _PSDmon1Chopper->_parameters.filename[0]='\0';
  #define filename (_PSDmon1Chopper->_parameters.filename)
  _PSDmon1Chopper->_parameters.xmin = -0.05;
  #define xmin (_PSDmon1Chopper->_parameters.xmin)
  _PSDmon1Chopper->_parameters.xmax = 0.05;
  #define xmax (_PSDmon1Chopper->_parameters.xmax)
  _PSDmon1Chopper->_parameters.ymin = -0.05;
  #define ymin (_PSDmon1Chopper->_parameters.ymin)
  _PSDmon1Chopper->_parameters.ymax = 0.05;
  #define ymax (_PSDmon1Chopper->_parameters.ymax)
  _PSDmon1Chopper->_parameters.xwidth = 0.1;
  #define xwidth (_PSDmon1Chopper->_parameters.xwidth)
  _PSDmon1Chopper->_parameters.yheight = 0.2;
  #define yheight (_PSDmon1Chopper->_parameters.yheight)
  _PSDmon1Chopper->_parameters.restore_neutron = 0;
  #define restore_neutron (_PSDmon1Chopper->_parameters.restore_neutron)

  #define PSD_N (_PSDmon1Chopper->_parameters.PSD_N)
  #define PSD_p (_PSDmon1Chopper->_parameters.PSD_p)
  #define PSD_p2 (_PSDmon1Chopper->_parameters.PSD_p2)

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
  /* component PSDmon1Chopper=PSD_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (0)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _a1->_rotation_absolute, _PSDmon1Chopper->_rotation_absolute);
    rot_transpose(_DISCTOF->_rotation_absolute, tr1);
    rot_mul(_PSDmon1Chopper->_rotation_absolute, tr1, _PSDmon1Chopper->_rotation_relative);
    _PSDmon1Chopper->_rotation_is_identity =  rot_test_identity(_PSDmon1Chopper->_rotation_relative);
    tc1 = coords_set(
      0, 0, 69.719);
    rot_transpose(_a1->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _PSDmon1Chopper->_position_absolute = coords_add(_a1->_position_absolute, tc2);
    tc1 = coords_sub(_DISCTOF->_position_absolute, _PSDmon1Chopper->_position_absolute);
    _PSDmon1Chopper->_position_relative = rot_apply(_PSDmon1Chopper->_rotation_absolute, tc1);
  } /* PSDmon1Chopper=PSD_monitor() AT ROTATED */
  DEBUG_COMPONENT("PSDmon1Chopper", _PSDmon1Chopper->_position_absolute, _PSDmon1Chopper->_rotation_absolute);
  instrument->_position_absolute[13] = _PSDmon1Chopper->_position_absolute;
  instrument->_position_relative[13] = _PSDmon1Chopper->_position_relative;
  instrument->counter_N[13]  = instrument->counter_P[13] = instrument->counter_P2[13] = 0;
  instrument->counter_AbsorbProp[13]= 0;
  return(0);
} /* _PSDmon1Chopper_setpos */

/* component VacuumTube1entry=Slit() SETTING, POSITION/ROTATION */
int _VacuumTube1entry_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_VacuumTube1entry_setpos] component VacuumTube1entry=Slit() SETTING [/usr/share/mcstas/3.0-dev/optics/Slit.comp:47]");
  stracpy(_VacuumTube1entry->_name, "VacuumTube1entry", 16384);
  stracpy(_VacuumTube1entry->_type, "Slit", 16384);
  _VacuumTube1entry->_index=14;
  _VacuumTube1entry->_parameters.xmin = -0.01;
  #define xmin (_VacuumTube1entry->_parameters.xmin)
  _VacuumTube1entry->_parameters.xmax = 0.01;
  #define xmax (_VacuumTube1entry->_parameters.xmax)
  _VacuumTube1entry->_parameters.ymin = -0.01;
  #define ymin (_VacuumTube1entry->_parameters.ymin)
  _VacuumTube1entry->_parameters.ymax = 0.01;
  #define ymax (_VacuumTube1entry->_parameters.ymax)
  _VacuumTube1entry->_parameters.radius = 0.085;
  #define radius (_VacuumTube1entry->_parameters.radius)
  _VacuumTube1entry->_parameters.xwidth = 0;
  #define xwidth (_VacuumTube1entry->_parameters.xwidth)
  _VacuumTube1entry->_parameters.yheight = 0;
  #define yheight (_VacuumTube1entry->_parameters.yheight)

  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef radius
  #undef xwidth
  #undef yheight
  /* component VacuumTube1entry=Slit() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _a1->_rotation_absolute, _VacuumTube1entry->_rotation_absolute);
    rot_transpose(_PSDmon1Chopper->_rotation_absolute, tr1);
    rot_mul(_VacuumTube1entry->_rotation_absolute, tr1, _VacuumTube1entry->_rotation_relative);
    _VacuumTube1entry->_rotation_is_identity =  rot_test_identity(_VacuumTube1entry->_rotation_relative);
    tc1 = coords_set(
      0, 0, 70.072);
    rot_transpose(_a1->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _VacuumTube1entry->_position_absolute = coords_add(_a1->_position_absolute, tc2);
    tc1 = coords_sub(_PSDmon1Chopper->_position_absolute, _VacuumTube1entry->_position_absolute);
    _VacuumTube1entry->_position_relative = rot_apply(_VacuumTube1entry->_rotation_absolute, tc1);
  } /* VacuumTube1entry=Slit() AT ROTATED */
  DEBUG_COMPONENT("VacuumTube1entry", _VacuumTube1entry->_position_absolute, _VacuumTube1entry->_rotation_absolute);
  instrument->_position_absolute[14] = _VacuumTube1entry->_position_absolute;
  instrument->_position_relative[14] = _VacuumTube1entry->_position_relative;
  instrument->counter_N[14]  = instrument->counter_P[14] = instrument->counter_P2[14] = 0;
  instrument->counter_AbsorbProp[14]= 0;
  return(0);
} /* _VacuumTube1entry_setpos */

/* component VacuumTube1exit=Slit() SETTING, POSITION/ROTATION */
int _VacuumTube1exit_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_VacuumTube1exit_setpos] component VacuumTube1exit=Slit() SETTING [/usr/share/mcstas/3.0-dev/optics/Slit.comp:47]");
  stracpy(_VacuumTube1exit->_name, "VacuumTube1exit", 16384);
  stracpy(_VacuumTube1exit->_type, "Slit", 16384);
  _VacuumTube1exit->_index=15;
  _VacuumTube1exit->_parameters.xmin = -0.01;
  #define xmin (_VacuumTube1exit->_parameters.xmin)
  _VacuumTube1exit->_parameters.xmax = 0.01;
  #define xmax (_VacuumTube1exit->_parameters.xmax)
  _VacuumTube1exit->_parameters.ymin = -0.01;
  #define ymin (_VacuumTube1exit->_parameters.ymin)
  _VacuumTube1exit->_parameters.ymax = 0.01;
  #define ymax (_VacuumTube1exit->_parameters.ymax)
  _VacuumTube1exit->_parameters.radius = 0.085;
  #define radius (_VacuumTube1exit->_parameters.radius)
  _VacuumTube1exit->_parameters.xwidth = 0;
  #define xwidth (_VacuumTube1exit->_parameters.xwidth)
  _VacuumTube1exit->_parameters.yheight = 0;
  #define yheight (_VacuumTube1exit->_parameters.yheight)

  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef radius
  #undef xwidth
  #undef yheight
  /* component VacuumTube1exit=Slit() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _a1->_rotation_absolute, _VacuumTube1exit->_rotation_absolute);
    rot_transpose(_VacuumTube1entry->_rotation_absolute, tr1);
    rot_mul(_VacuumTube1exit->_rotation_absolute, tr1, _VacuumTube1exit->_rotation_relative);
    _VacuumTube1exit->_rotation_is_identity =  rot_test_identity(_VacuumTube1exit->_rotation_relative);
    tc1 = coords_set(
      0, 0, 70.607);
    rot_transpose(_a1->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _VacuumTube1exit->_position_absolute = coords_add(_a1->_position_absolute, tc2);
    tc1 = coords_sub(_VacuumTube1entry->_position_absolute, _VacuumTube1exit->_position_absolute);
    _VacuumTube1exit->_position_relative = rot_apply(_VacuumTube1exit->_rotation_absolute, tc1);
  } /* VacuumTube1exit=Slit() AT ROTATED */
  DEBUG_COMPONENT("VacuumTube1exit", _VacuumTube1exit->_position_absolute, _VacuumTube1exit->_rotation_absolute);
  instrument->_position_absolute[15] = _VacuumTube1exit->_position_absolute;
  instrument->_position_relative[15] = _VacuumTube1exit->_position_relative;
  instrument->counter_N[15]  = instrument->counter_P[15] = instrument->counter_P2[15] = 0;
  instrument->counter_AbsorbProp[15]= 0;
  return(0);
} /* _VacuumTube1exit_setpos */

/* component VacuumTube2entry=Slit() SETTING, POSITION/ROTATION */
int _VacuumTube2entry_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_VacuumTube2entry_setpos] component VacuumTube2entry=Slit() SETTING [/usr/share/mcstas/3.0-dev/optics/Slit.comp:47]");
  stracpy(_VacuumTube2entry->_name, "VacuumTube2entry", 16384);
  stracpy(_VacuumTube2entry->_type, "Slit", 16384);
  _VacuumTube2entry->_index=16;
  _VacuumTube2entry->_parameters.xmin = -0.01;
  #define xmin (_VacuumTube2entry->_parameters.xmin)
  _VacuumTube2entry->_parameters.xmax = 0.01;
  #define xmax (_VacuumTube2entry->_parameters.xmax)
  _VacuumTube2entry->_parameters.ymin = -0.01;
  #define ymin (_VacuumTube2entry->_parameters.ymin)
  _VacuumTube2entry->_parameters.ymax = 0.01;
  #define ymax (_VacuumTube2entry->_parameters.ymax)
  _VacuumTube2entry->_parameters.radius = 0.1;
  #define radius (_VacuumTube2entry->_parameters.radius)
  _VacuumTube2entry->_parameters.xwidth = 0;
  #define xwidth (_VacuumTube2entry->_parameters.xwidth)
  _VacuumTube2entry->_parameters.yheight = 0;
  #define yheight (_VacuumTube2entry->_parameters.yheight)

  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef radius
  #undef xwidth
  #undef yheight
  /* component VacuumTube2entry=Slit() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _a1->_rotation_absolute, _VacuumTube2entry->_rotation_absolute);
    rot_transpose(_VacuumTube1exit->_rotation_absolute, tr1);
    rot_mul(_VacuumTube2entry->_rotation_absolute, tr1, _VacuumTube2entry->_rotation_relative);
    _VacuumTube2entry->_rotation_is_identity =  rot_test_identity(_VacuumTube2entry->_rotation_relative);
    tc1 = coords_set(
      0, 0, 70.607);
    rot_transpose(_a1->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _VacuumTube2entry->_position_absolute = coords_add(_a1->_position_absolute, tc2);
    tc1 = coords_sub(_VacuumTube1exit->_position_absolute, _VacuumTube2entry->_position_absolute);
    _VacuumTube2entry->_position_relative = rot_apply(_VacuumTube2entry->_rotation_absolute, tc1);
  } /* VacuumTube2entry=Slit() AT ROTATED */
  DEBUG_COMPONENT("VacuumTube2entry", _VacuumTube2entry->_position_absolute, _VacuumTube2entry->_rotation_absolute);
  instrument->_position_absolute[16] = _VacuumTube2entry->_position_absolute;
  instrument->_position_relative[16] = _VacuumTube2entry->_position_relative;
  instrument->counter_N[16]  = instrument->counter_P[16] = instrument->counter_P2[16] = 0;
  instrument->counter_AbsorbProp[16]= 0;
  return(0);
} /* _VacuumTube2entry_setpos */

/* component VacuumTube2exit=Slit() SETTING, POSITION/ROTATION */
int _VacuumTube2exit_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_VacuumTube2exit_setpos] component VacuumTube2exit=Slit() SETTING [/usr/share/mcstas/3.0-dev/optics/Slit.comp:47]");
  stracpy(_VacuumTube2exit->_name, "VacuumTube2exit", 16384);
  stracpy(_VacuumTube2exit->_type, "Slit", 16384);
  _VacuumTube2exit->_index=17;
  _VacuumTube2exit->_parameters.xmin = -0.01;
  #define xmin (_VacuumTube2exit->_parameters.xmin)
  _VacuumTube2exit->_parameters.xmax = 0.01;
  #define xmax (_VacuumTube2exit->_parameters.xmax)
  _VacuumTube2exit->_parameters.ymin = -0.01;
  #define ymin (_VacuumTube2exit->_parameters.ymin)
  _VacuumTube2exit->_parameters.ymax = 0.01;
  #define ymax (_VacuumTube2exit->_parameters.ymax)
  _VacuumTube2exit->_parameters.radius = 0.1;
  #define radius (_VacuumTube2exit->_parameters.radius)
  _VacuumTube2exit->_parameters.xwidth = 0;
  #define xwidth (_VacuumTube2exit->_parameters.xwidth)
  _VacuumTube2exit->_parameters.yheight = 0;
  #define yheight (_VacuumTube2exit->_parameters.yheight)

  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef radius
  #undef xwidth
  #undef yheight
  /* component VacuumTube2exit=Slit() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _a1->_rotation_absolute, _VacuumTube2exit->_rotation_absolute);
    rot_transpose(_VacuumTube2entry->_rotation_absolute, tr1);
    rot_mul(_VacuumTube2exit->_rotation_absolute, tr1, _VacuumTube2exit->_rotation_relative);
    _VacuumTube2exit->_rotation_is_identity =  rot_test_identity(_VacuumTube2exit->_rotation_relative);
    tc1 = coords_set(
      0, 0, 71.402);
    rot_transpose(_a1->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _VacuumTube2exit->_position_absolute = coords_add(_a1->_position_absolute, tc2);
    tc1 = coords_sub(_VacuumTube2entry->_position_absolute, _VacuumTube2exit->_position_absolute);
    _VacuumTube2exit->_position_relative = rot_apply(_VacuumTube2exit->_rotation_absolute, tc1);
  } /* VacuumTube2exit=Slit() AT ROTATED */
  DEBUG_COMPONENT("VacuumTube2exit", _VacuumTube2exit->_position_absolute, _VacuumTube2exit->_rotation_absolute);
  instrument->_position_absolute[17] = _VacuumTube2exit->_position_absolute;
  instrument->_position_relative[17] = _VacuumTube2exit->_position_relative;
  instrument->counter_N[17]  = instrument->counter_P[17] = instrument->counter_P2[17] = 0;
  instrument->counter_AbsorbProp[17]= 0;
  return(0);
} /* _VacuumTube2exit_setpos */

/* component VacuumTube3entry=Slit() SETTING, POSITION/ROTATION */
int _VacuumTube3entry_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_VacuumTube3entry_setpos] component VacuumTube3entry=Slit() SETTING [/usr/share/mcstas/3.0-dev/optics/Slit.comp:47]");
  stracpy(_VacuumTube3entry->_name, "VacuumTube3entry", 16384);
  stracpy(_VacuumTube3entry->_type, "Slit", 16384);
  _VacuumTube3entry->_index=18;
  _VacuumTube3entry->_parameters.xmin = -0.01;
  #define xmin (_VacuumTube3entry->_parameters.xmin)
  _VacuumTube3entry->_parameters.xmax = 0.01;
  #define xmax (_VacuumTube3entry->_parameters.xmax)
  _VacuumTube3entry->_parameters.ymin = -0.01;
  #define ymin (_VacuumTube3entry->_parameters.ymin)
  _VacuumTube3entry->_parameters.ymax = 0.01;
  #define ymax (_VacuumTube3entry->_parameters.ymax)
  _VacuumTube3entry->_parameters.radius = 0.1335;
  #define radius (_VacuumTube3entry->_parameters.radius)
  _VacuumTube3entry->_parameters.xwidth = 0;
  #define xwidth (_VacuumTube3entry->_parameters.xwidth)
  _VacuumTube3entry->_parameters.yheight = 0;
  #define yheight (_VacuumTube3entry->_parameters.yheight)

  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef radius
  #undef xwidth
  #undef yheight
  /* component VacuumTube3entry=Slit() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _a1->_rotation_absolute, _VacuumTube3entry->_rotation_absolute);
    rot_transpose(_VacuumTube2exit->_rotation_absolute, tr1);
    rot_mul(_VacuumTube3entry->_rotation_absolute, tr1, _VacuumTube3entry->_rotation_relative);
    _VacuumTube3entry->_rotation_is_identity =  rot_test_identity(_VacuumTube3entry->_rotation_relative);
    tc1 = coords_set(
      0, 0, 71.402);
    rot_transpose(_a1->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _VacuumTube3entry->_position_absolute = coords_add(_a1->_position_absolute, tc2);
    tc1 = coords_sub(_VacuumTube2exit->_position_absolute, _VacuumTube3entry->_position_absolute);
    _VacuumTube3entry->_position_relative = rot_apply(_VacuumTube3entry->_rotation_absolute, tc1);
  } /* VacuumTube3entry=Slit() AT ROTATED */
  DEBUG_COMPONENT("VacuumTube3entry", _VacuumTube3entry->_position_absolute, _VacuumTube3entry->_rotation_absolute);
  instrument->_position_absolute[18] = _VacuumTube3entry->_position_absolute;
  instrument->_position_relative[18] = _VacuumTube3entry->_position_relative;
  instrument->counter_N[18]  = instrument->counter_P[18] = instrument->counter_P2[18] = 0;
  instrument->counter_AbsorbProp[18]= 0;
  return(0);
} /* _VacuumTube3entry_setpos */

/* component VacuumTube3exit=Slit() SETTING, POSITION/ROTATION */
int _VacuumTube3exit_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_VacuumTube3exit_setpos] component VacuumTube3exit=Slit() SETTING [/usr/share/mcstas/3.0-dev/optics/Slit.comp:47]");
  stracpy(_VacuumTube3exit->_name, "VacuumTube3exit", 16384);
  stracpy(_VacuumTube3exit->_type, "Slit", 16384);
  _VacuumTube3exit->_index=19;
  _VacuumTube3exit->_parameters.xmin = -0.01;
  #define xmin (_VacuumTube3exit->_parameters.xmin)
  _VacuumTube3exit->_parameters.xmax = 0.01;
  #define xmax (_VacuumTube3exit->_parameters.xmax)
  _VacuumTube3exit->_parameters.ymin = -0.01;
  #define ymin (_VacuumTube3exit->_parameters.ymin)
  _VacuumTube3exit->_parameters.ymax = 0.01;
  #define ymax (_VacuumTube3exit->_parameters.ymax)
  _VacuumTube3exit->_parameters.radius = 0.1335;
  #define radius (_VacuumTube3exit->_parameters.radius)
  _VacuumTube3exit->_parameters.xwidth = 0;
  #define xwidth (_VacuumTube3exit->_parameters.xwidth)
  _VacuumTube3exit->_parameters.yheight = 0;
  #define yheight (_VacuumTube3exit->_parameters.yheight)

  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef radius
  #undef xwidth
  #undef yheight
  /* component VacuumTube3exit=Slit() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _a1->_rotation_absolute, _VacuumTube3exit->_rotation_absolute);
    rot_transpose(_VacuumTube3entry->_rotation_absolute, tr1);
    rot_mul(_VacuumTube3exit->_rotation_absolute, tr1, _VacuumTube3exit->_rotation_relative);
    _VacuumTube3exit->_rotation_is_identity =  rot_test_identity(_VacuumTube3exit->_rotation_relative);
    tc1 = coords_set(
      0, 0, 72.007);
    rot_transpose(_a1->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _VacuumTube3exit->_position_absolute = coords_add(_a1->_position_absolute, tc2);
    tc1 = coords_sub(_VacuumTube3entry->_position_absolute, _VacuumTube3exit->_position_absolute);
    _VacuumTube3exit->_position_relative = rot_apply(_VacuumTube3exit->_rotation_absolute, tc1);
  } /* VacuumTube3exit=Slit() AT ROTATED */
  DEBUG_COMPONENT("VacuumTube3exit", _VacuumTube3exit->_position_absolute, _VacuumTube3exit->_rotation_absolute);
  instrument->_position_absolute[19] = _VacuumTube3exit->_position_absolute;
  instrument->_position_relative[19] = _VacuumTube3exit->_position_relative;
  instrument->counter_N[19]  = instrument->counter_P[19] = instrument->counter_P2[19] = 0;
  instrument->counter_AbsorbProp[19]= 0;
  return(0);
} /* _VacuumTube3exit_setpos */

/* component PSDmonMono=PSD_monitor() SETTING, POSITION/ROTATION */
int _PSDmonMono_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_PSDmonMono_setpos] component PSDmonMono=PSD_monitor() SETTING [/usr/share/mcstas/3.0-dev/monitors/PSD_monitor.comp:64]");
  stracpy(_PSDmonMono->_name, "PSDmonMono", 16384);
  stracpy(_PSDmonMono->_type, "PSD_monitor", 16384);
  _PSDmonMono->_index=20;
  _PSDmonMono->_parameters.nx = 60;
  #define nx (_PSDmonMono->_parameters.nx)
  _PSDmonMono->_parameters.ny = 60;
  #define ny (_PSDmonMono->_parameters.ny)
  if("psdmono.dat" && strlen("psdmono.dat"))
    stracpy(_PSDmonMono->_parameters.filename, "psdmono.dat" ? "psdmono.dat" : "", 16384);
  else 
  _PSDmonMono->_parameters.filename[0]='\0';
  #define filename (_PSDmonMono->_parameters.filename)
  _PSDmonMono->_parameters.xmin = -0.05;
  #define xmin (_PSDmonMono->_parameters.xmin)
  _PSDmonMono->_parameters.xmax = 0.05;
  #define xmax (_PSDmonMono->_parameters.xmax)
  _PSDmonMono->_parameters.ymin = -0.05;
  #define ymin (_PSDmonMono->_parameters.ymin)
  _PSDmonMono->_parameters.ymax = 0.05;
  #define ymax (_PSDmonMono->_parameters.ymax)
  _PSDmonMono->_parameters.xwidth = 0.1;
  #define xwidth (_PSDmonMono->_parameters.xwidth)
  _PSDmonMono->_parameters.yheight = 0.2;
  #define yheight (_PSDmonMono->_parameters.yheight)
  _PSDmonMono->_parameters.restore_neutron = 0;
  #define restore_neutron (_PSDmonMono->_parameters.restore_neutron)

  #define PSD_N (_PSDmonMono->_parameters.PSD_N)
  #define PSD_p (_PSDmonMono->_parameters.PSD_p)
  #define PSD_p2 (_PSDmonMono->_parameters.PSD_p2)

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
  /* component PSDmonMono=PSD_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (0)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _a1->_rotation_absolute, _PSDmonMono->_rotation_absolute);
    rot_transpose(_VacuumTube3exit->_rotation_absolute, tr1);
    rot_mul(_PSDmonMono->_rotation_absolute, tr1, _PSDmonMono->_rotation_relative);
    _PSDmonMono->_rotation_is_identity =  rot_test_identity(_PSDmonMono->_rotation_relative);
    tc1 = coords_set(
      0, 0, 72.4);
    rot_transpose(_a1->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _PSDmonMono->_position_absolute = coords_add(_a1->_position_absolute, tc2);
    tc1 = coords_sub(_VacuumTube3exit->_position_absolute, _PSDmonMono->_position_absolute);
    _PSDmonMono->_position_relative = rot_apply(_PSDmonMono->_rotation_absolute, tc1);
  } /* PSDmonMono=PSD_monitor() AT ROTATED */
  DEBUG_COMPONENT("PSDmonMono", _PSDmonMono->_position_absolute, _PSDmonMono->_rotation_absolute);
  instrument->_position_absolute[20] = _PSDmonMono->_position_absolute;
  instrument->_position_relative[20] = _PSDmonMono->_position_relative;
  instrument->counter_N[20]  = instrument->counter_P[20] = instrument->counter_P2[20] = 0;
  instrument->counter_AbsorbProp[20]= 0;
  return(0);
} /* _PSDmonMono_setpos */

/* component MONOTOF=Monitor_nD() SETTING, POSITION/ROTATION */
int _MONOTOF_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_MONOTOF_setpos] component MONOTOF=Monitor_nD() SETTING [/usr/share/mcstas/3.0-dev/monitors/Monitor_nD.comp:228]");
  stracpy(_MONOTOF->_name, "MONOTOF", 16384);
  stracpy(_MONOTOF->_type, "Monitor_nD", 16384);
  _MONOTOF->_index=21;
  _MONOTOF->_parameters.user1 = '\0';
  #define user1 (_MONOTOF->_parameters.user1)
  _MONOTOF->_parameters.user2 = '\0';
  #define user2 (_MONOTOF->_parameters.user2)
  _MONOTOF->_parameters.user3 = '\0';
  #define user3 (_MONOTOF->_parameters.user3)
  _MONOTOF->_parameters.xwidth = 0;
  #define xwidth (_MONOTOF->_parameters.xwidth)
  _MONOTOF->_parameters.yheight = 0;
  #define yheight (_MONOTOF->_parameters.yheight)
  _MONOTOF->_parameters.zdepth = 0;
  #define zdepth (_MONOTOF->_parameters.zdepth)
  _MONOTOF->_parameters.xmin = -0.05;
  #define xmin (_MONOTOF->_parameters.xmin)
  _MONOTOF->_parameters.xmax = 0.05;
  #define xmax (_MONOTOF->_parameters.xmax)
  _MONOTOF->_parameters.ymin = -0.1;
  #define ymin (_MONOTOF->_parameters.ymin)
  _MONOTOF->_parameters.ymax = 0.1;
  #define ymax (_MONOTOF->_parameters.ymax)
  _MONOTOF->_parameters.zmin = 0;
  #define zmin (_MONOTOF->_parameters.zmin)
  _MONOTOF->_parameters.zmax = 0;
  #define zmax (_MONOTOF->_parameters.zmax)
  _MONOTOF->_parameters.bins = 30;
  #define bins (_MONOTOF->_parameters.bins)
  _MONOTOF->_parameters.min = -1e40;
  #define min (_MONOTOF->_parameters.min)
  _MONOTOF->_parameters.max = 1e40;
  #define max (_MONOTOF->_parameters.max)
  _MONOTOF->_parameters.restore_neutron = 0;
  #define restore_neutron (_MONOTOF->_parameters.restore_neutron)
  _MONOTOF->_parameters.radius = 0;
  #define radius (_MONOTOF->_parameters.radius)
  if("auto, time" && strlen("auto, time"))
    stracpy(_MONOTOF->_parameters.options, "auto, time" ? "auto, time" : "", 16384);
  else 
  _MONOTOF->_parameters.options[0]='\0';
  #define options (_MONOTOF->_parameters.options)
  if("MONO_TOF.dat" && strlen("MONO_TOF.dat"))
    stracpy(_MONOTOF->_parameters.filename, "MONO_TOF.dat" ? "MONO_TOF.dat" : "", 16384);
  else 
  _MONOTOF->_parameters.filename[0]='\0';
  #define filename (_MONOTOF->_parameters.filename)
  if("NULL" && strlen("NULL"))
    stracpy(_MONOTOF->_parameters.geometry, "NULL" ? "NULL" : "", 16384);
  else 
  _MONOTOF->_parameters.geometry[0]='\0';
  #define geometry (_MONOTOF->_parameters.geometry)
  if("NULL" && strlen("NULL"))
    stracpy(_MONOTOF->_parameters.username1, "NULL" ? "NULL" : "", 16384);
  else 
  _MONOTOF->_parameters.username1[0]='\0';
  #define username1 (_MONOTOF->_parameters.username1)
  if("NULL" && strlen("NULL"))
    stracpy(_MONOTOF->_parameters.username2, "NULL" ? "NULL" : "", 16384);
  else 
  _MONOTOF->_parameters.username2[0]='\0';
  #define username2 (_MONOTOF->_parameters.username2)
  if("NULL" && strlen("NULL"))
    stracpy(_MONOTOF->_parameters.username3, "NULL" ? "NULL" : "", 16384);
  else 
  _MONOTOF->_parameters.username3[0]='\0';
  #define username3 (_MONOTOF->_parameters.username3)

  #define DEFS (_MONOTOF->_parameters.DEFS)
  #define Vars (_MONOTOF->_parameters.Vars)
  #define detector (_MONOTOF->_parameters.detector)
  #define offdata (_MONOTOF->_parameters.offdata)

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
  /* component MONOTOF=Monitor_nD() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _a1->_rotation_absolute, _MONOTOF->_rotation_absolute);
    rot_transpose(_PSDmonMono->_rotation_absolute, tr1);
    rot_mul(_MONOTOF->_rotation_absolute, tr1, _MONOTOF->_rotation_relative);
    _MONOTOF->_rotation_is_identity =  rot_test_identity(_MONOTOF->_rotation_relative);
    tc1 = coords_set(
      0, 0, 72.41);
    rot_transpose(_a1->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _MONOTOF->_position_absolute = coords_add(_a1->_position_absolute, tc2);
    tc1 = coords_sub(_PSDmonMono->_position_absolute, _MONOTOF->_position_absolute);
    _MONOTOF->_position_relative = rot_apply(_MONOTOF->_rotation_absolute, tc1);
  } /* MONOTOF=Monitor_nD() AT ROTATED */
  DEBUG_COMPONENT("MONOTOF", _MONOTOF->_position_absolute, _MONOTOF->_rotation_absolute);
  instrument->_position_absolute[21] = _MONOTOF->_position_absolute;
  instrument->_position_relative[21] = _MONOTOF->_position_relative;
  instrument->counter_N[21]  = instrument->counter_P[21] = instrument->counter_P2[21] = 0;
  instrument->counter_AbsorbProp[21]= 0;
  return(0);
} /* _MONOTOF_setpos */

/* component DivMonMono=Divergence_monitor() SETTING, POSITION/ROTATION */
int _DivMonMono_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_DivMonMono_setpos] component DivMonMono=Divergence_monitor() SETTING [/usr/share/mcstas/3.0-dev/monitors/Divergence_monitor.comp:71]");
  stracpy(_DivMonMono->_name, "DivMonMono", 16384);
  stracpy(_DivMonMono->_type, "Divergence_monitor", 16384);
  _DivMonMono->_index=22;
  _DivMonMono->_parameters.nh = 60;
  #define nh (_DivMonMono->_parameters.nh)
  _DivMonMono->_parameters.nv = 60;
  #define nv (_DivMonMono->_parameters.nv)
  if("divmono.dat" && strlen("divmono.dat"))
    stracpy(_DivMonMono->_parameters.filename, "divmono.dat" ? "divmono.dat" : "", 16384);
  else 
  _DivMonMono->_parameters.filename[0]='\0';
  #define filename (_DivMonMono->_parameters.filename)
  _DivMonMono->_parameters.xmin = -0.05;
  #define xmin (_DivMonMono->_parameters.xmin)
  _DivMonMono->_parameters.xmax = 0.05;
  #define xmax (_DivMonMono->_parameters.xmax)
  _DivMonMono->_parameters.ymin = -0.05;
  #define ymin (_DivMonMono->_parameters.ymin)
  _DivMonMono->_parameters.ymax = 0.05;
  #define ymax (_DivMonMono->_parameters.ymax)
  _DivMonMono->_parameters.xwidth = 0.1;
  #define xwidth (_DivMonMono->_parameters.xwidth)
  _DivMonMono->_parameters.yheight = 0.2;
  #define yheight (_DivMonMono->_parameters.yheight)
  _DivMonMono->_parameters.maxdiv_h = 3;
  #define maxdiv_h (_DivMonMono->_parameters.maxdiv_h)
  _DivMonMono->_parameters.maxdiv_v = 3;
  #define maxdiv_v (_DivMonMono->_parameters.maxdiv_v)
  _DivMonMono->_parameters.restore_neutron = 0;
  #define restore_neutron (_DivMonMono->_parameters.restore_neutron)
  _DivMonMono->_parameters.nx = 0;
  #define nx (_DivMonMono->_parameters.nx)
  _DivMonMono->_parameters.ny = 0;
  #define ny (_DivMonMono->_parameters.ny)
  _DivMonMono->_parameters.nz = 1;
  #define nz (_DivMonMono->_parameters.nz)

  #define Div_N (_DivMonMono->_parameters.Div_N)
  #define Div_p (_DivMonMono->_parameters.Div_p)
  #define Div_p2 (_DivMonMono->_parameters.Div_p2)

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
  /* component DivMonMono=Divergence_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (0)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _a1->_rotation_absolute, _DivMonMono->_rotation_absolute);
    rot_transpose(_MONOTOF->_rotation_absolute, tr1);
    rot_mul(_DivMonMono->_rotation_absolute, tr1, _DivMonMono->_rotation_relative);
    _DivMonMono->_rotation_is_identity =  rot_test_identity(_DivMonMono->_rotation_relative);
    tc1 = coords_set(
      0, 0, 72.42);
    rot_transpose(_a1->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _DivMonMono->_position_absolute = coords_add(_a1->_position_absolute, tc2);
    tc1 = coords_sub(_MONOTOF->_position_absolute, _DivMonMono->_position_absolute);
    _DivMonMono->_position_relative = rot_apply(_DivMonMono->_rotation_absolute, tc1);
  } /* DivMonMono=Divergence_monitor() AT ROTATED */
  DEBUG_COMPONENT("DivMonMono", _DivMonMono->_position_absolute, _DivMonMono->_rotation_absolute);
  instrument->_position_absolute[22] = _DivMonMono->_position_absolute;
  instrument->_position_relative[22] = _DivMonMono->_position_relative;
  instrument->counter_N[22]  = instrument->counter_P[22] = instrument->counter_P2[22] = 0;
  instrument->counter_AbsorbProp[22]= 0;
  return(0);
} /* _DivMonMono_setpos */

/* component focus_mono=Arm() SETTING, POSITION/ROTATION */
int _focus_mono_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_focus_mono_setpos] component focus_mono=Arm() SETTING [Arm:0]");
  stracpy(_focus_mono->_name, "focus_mono", 16384);
  stracpy(_focus_mono->_type, "Arm", 16384);
  _focus_mono->_index=23;
  /* component focus_mono=Arm() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (PHM)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _a1->_rotation_absolute, _focus_mono->_rotation_absolute);
    rot_transpose(_DivMonMono->_rotation_absolute, tr1);
    rot_mul(_focus_mono->_rotation_absolute, tr1, _focus_mono->_rotation_relative);
    _focus_mono->_rotation_is_identity =  rot_test_identity(_focus_mono->_rotation_relative);
    tc1 = coords_set(
      0, 0, 72.617);
    rot_transpose(_a1->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _focus_mono->_position_absolute = coords_add(_a1->_position_absolute, tc2);
    tc1 = coords_sub(_DivMonMono->_position_absolute, _focus_mono->_position_absolute);
    _focus_mono->_position_relative = rot_apply(_focus_mono->_rotation_absolute, tc1);
  } /* focus_mono=Arm() AT ROTATED */
  DEBUG_COMPONENT("focus_mono", _focus_mono->_position_absolute, _focus_mono->_rotation_absolute);
  instrument->_position_absolute[23] = _focus_mono->_position_absolute;
  instrument->_position_relative[23] = _focus_mono->_position_relative;
  instrument->counter_N[23]  = instrument->counter_P[23] = instrument->counter_P2[23] = 0;
  instrument->counter_AbsorbProp[23]= 0;
  return(0);
} /* _focus_mono_setpos */

/* component mono=Monochromator_2foc() SETTING, POSITION/ROTATION */
int _mono_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_mono_setpos] component mono=Monochromator_2foc() SETTING [/usr/share/mcstas/3.0-dev/contrib/Monochromator_2foc.comp:100]");
  stracpy(_mono->_name, "mono", 16384);
  stracpy(_mono->_type, "Monochromator_2foc", 16384);
  _mono->_index=24;
  _mono->_parameters.reflect[0]='\0';
  #define reflect (_mono->_parameters.reflect)
  _mono->_parameters.zwidth = 0.019;
  #define zwidth (_mono->_parameters.zwidth)
  _mono->_parameters.yheight = 0.025;
  #define yheight (_mono->_parameters.yheight)
  _mono->_parameters.gap = 0.001;
  #define gap (_mono->_parameters.gap)
  _mono->_parameters.NH = 9;
  #define NH (_mono->_parameters.NH)
  _mono->_parameters.NV = 7;
  #define NV (_mono->_parameters.NV)
  _mono->_parameters.mosaich = 48;
  #define mosaich (_mono->_parameters.mosaich)
  _mono->_parameters.mosaicv = 48;
  #define mosaicv (_mono->_parameters.mosaicv)
  _mono->_parameters.r0 = 0.99;
  #define r0 (_mono->_parameters.r0)
  _mono->_parameters.Q = 1.873;
  #define Q (_mono->_parameters.Q)
  _mono->_parameters.RV = RV_2;
  #define RV (_mono->_parameters.RV)
  _mono->_parameters.RH = RH_2;
  #define RH (_mono->_parameters.RH)
  _mono->_parameters.DM = 0;
  #define DM (_mono->_parameters.DM)
  _mono->_parameters.mosaic = 0;
  #define mosaic (_mono->_parameters.mosaic)
  _mono->_parameters.width = 0;
  #define width (_mono->_parameters.width)
  _mono->_parameters.height = 0;
  #define height (_mono->_parameters.height)
  _mono->_parameters.verbose = 0;
  #define verbose (_mono->_parameters.verbose)

  #define mos_y (_mono->_parameters.mos_y)
  #define mos_z (_mono->_parameters.mos_z)
  #define mono_Q (_mono->_parameters.mono_Q)
  #define SlabWidth (_mono->_parameters.SlabWidth)
  #define SlabHeight (_mono->_parameters.SlabHeight)
  #define rTable (_mono->_parameters.rTable)

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
  /* component mono=Monochromator_2foc() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _focus_mono->_rotation_absolute, _mono->_rotation_absolute);
    rot_transpose(_focus_mono->_rotation_absolute, tr1);
    rot_mul(_mono->_rotation_absolute, tr1, _mono->_rotation_relative);
    _mono->_rotation_is_identity =  rot_test_identity(_mono->_rotation_relative);
    tc1 = coords_set(
      0, 0, 0);
    rot_transpose(_focus_mono->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _mono->_position_absolute = coords_add(_focus_mono->_position_absolute, tc2);
    tc1 = coords_sub(_focus_mono->_position_absolute, _mono->_position_absolute);
    _mono->_position_relative = rot_apply(_mono->_rotation_absolute, tc1);
  } /* mono=Monochromator_2foc() AT ROTATED */
  DEBUG_COMPONENT("mono", _mono->_position_absolute, _mono->_rotation_absolute);
  instrument->_position_absolute[24] = _mono->_position_absolute;
  instrument->_position_relative[24] = _mono->_position_relative;
  instrument->counter_N[24]  = instrument->counter_P[24] = instrument->counter_P2[24] = 0;
  instrument->counter_AbsorbProp[24]= 0;
  return(0);
} /* _mono_setpos */

/* component a2=Arm() SETTING, POSITION/ROTATION */
int _a2_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_a2_setpos] component a2=Arm() SETTING [Arm:0]");
  stracpy(_a2->_name, "a2", 16384);
  stracpy(_a2->_type, "Arm", 16384);
  _a2->_index=25;
  /* component a2=Arm() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (TTM)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _a1->_rotation_absolute, _a2->_rotation_absolute);
    rot_transpose(_mono->_rotation_absolute, tr1);
    rot_mul(_a2->_rotation_absolute, tr1, _a2->_rotation_relative);
    _a2->_rotation_is_identity =  rot_test_identity(_a2->_rotation_relative);
    tc1 = coords_set(
      0, 0, 0);
    rot_transpose(_focus_mono->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _a2->_position_absolute = coords_add(_focus_mono->_position_absolute, tc2);
    tc1 = coords_sub(_mono->_position_absolute, _a2->_position_absolute);
    _a2->_position_relative = rot_apply(_a2->_rotation_absolute, tc1);
  } /* a2=Arm() AT ROTATED */
  DEBUG_COMPONENT("a2", _a2->_position_absolute, _a2->_rotation_absolute);
  instrument->_position_absolute[25] = _a2->_position_absolute;
  instrument->_position_relative[25] = _a2->_position_relative;
  instrument->counter_N[25]  = instrument->counter_P[25] = instrument->counter_P2[25] = 0;
  instrument->counter_AbsorbProp[25]= 0;
  return(0);
} /* _a2_setpos */

/* component FERMITOF_before=Monitor_nD() SETTING, POSITION/ROTATION */
int _FERMITOF_before_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_FERMITOF_before_setpos] component FERMITOF_before=Monitor_nD() SETTING [/usr/share/mcstas/3.0-dev/monitors/Monitor_nD.comp:228]");
  stracpy(_FERMITOF_before->_name, "FERMITOF_before", 16384);
  stracpy(_FERMITOF_before->_type, "Monitor_nD", 16384);
  _FERMITOF_before->_index=26;
  _FERMITOF_before->_parameters.user1 = '\0';
  #define user1 (_FERMITOF_before->_parameters.user1)
  _FERMITOF_before->_parameters.user2 = '\0';
  #define user2 (_FERMITOF_before->_parameters.user2)
  _FERMITOF_before->_parameters.user3 = '\0';
  #define user3 (_FERMITOF_before->_parameters.user3)
  _FERMITOF_before->_parameters.xwidth = 0;
  #define xwidth (_FERMITOF_before->_parameters.xwidth)
  _FERMITOF_before->_parameters.yheight = 0;
  #define yheight (_FERMITOF_before->_parameters.yheight)
  _FERMITOF_before->_parameters.zdepth = 0;
  #define zdepth (_FERMITOF_before->_parameters.zdepth)
  _FERMITOF_before->_parameters.xmin = -0.05;
  #define xmin (_FERMITOF_before->_parameters.xmin)
  _FERMITOF_before->_parameters.xmax = 0.05;
  #define xmax (_FERMITOF_before->_parameters.xmax)
  _FERMITOF_before->_parameters.ymin = -0.1;
  #define ymin (_FERMITOF_before->_parameters.ymin)
  _FERMITOF_before->_parameters.ymax = 0.1;
  #define ymax (_FERMITOF_before->_parameters.ymax)
  _FERMITOF_before->_parameters.zmin = 0;
  #define zmin (_FERMITOF_before->_parameters.zmin)
  _FERMITOF_before->_parameters.zmax = 0;
  #define zmax (_FERMITOF_before->_parameters.zmax)
  _FERMITOF_before->_parameters.bins = 30;
  #define bins (_FERMITOF_before->_parameters.bins)
  _FERMITOF_before->_parameters.min = -1e40;
  #define min (_FERMITOF_before->_parameters.min)
  _FERMITOF_before->_parameters.max = 1e40;
  #define max (_FERMITOF_before->_parameters.max)
  _FERMITOF_before->_parameters.restore_neutron = 0;
  #define restore_neutron (_FERMITOF_before->_parameters.restore_neutron)
  _FERMITOF_before->_parameters.radius = 0;
  #define radius (_FERMITOF_before->_parameters.radius)
  if("auto, time" && strlen("auto, time"))
    stracpy(_FERMITOF_before->_parameters.options, "auto, time" ? "auto, time" : "", 16384);
  else 
  _FERMITOF_before->_parameters.options[0]='\0';
  #define options (_FERMITOF_before->_parameters.options)
  if("FERMI_TOF_before.dat" && strlen("FERMI_TOF_before.dat"))
    stracpy(_FERMITOF_before->_parameters.filename, "FERMI_TOF_before.dat" ? "FERMI_TOF_before.dat" : "", 16384);
  else 
  _FERMITOF_before->_parameters.filename[0]='\0';
  #define filename (_FERMITOF_before->_parameters.filename)
  if("NULL" && strlen("NULL"))
    stracpy(_FERMITOF_before->_parameters.geometry, "NULL" ? "NULL" : "", 16384);
  else 
  _FERMITOF_before->_parameters.geometry[0]='\0';
  #define geometry (_FERMITOF_before->_parameters.geometry)
  if("NULL" && strlen("NULL"))
    stracpy(_FERMITOF_before->_parameters.username1, "NULL" ? "NULL" : "", 16384);
  else 
  _FERMITOF_before->_parameters.username1[0]='\0';
  #define username1 (_FERMITOF_before->_parameters.username1)
  if("NULL" && strlen("NULL"))
    stracpy(_FERMITOF_before->_parameters.username2, "NULL" ? "NULL" : "", 16384);
  else 
  _FERMITOF_before->_parameters.username2[0]='\0';
  #define username2 (_FERMITOF_before->_parameters.username2)
  if("NULL" && strlen("NULL"))
    stracpy(_FERMITOF_before->_parameters.username3, "NULL" ? "NULL" : "", 16384);
  else 
  _FERMITOF_before->_parameters.username3[0]='\0';
  #define username3 (_FERMITOF_before->_parameters.username3)

  #define DEFS (_FERMITOF_before->_parameters.DEFS)
  #define Vars (_FERMITOF_before->_parameters.Vars)
  #define detector (_FERMITOF_before->_parameters.detector)
  #define offdata (_FERMITOF_before->_parameters.offdata)

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
  /* component FERMITOF_before=Monitor_nD() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _a2->_rotation_absolute, _FERMITOF_before->_rotation_absolute);
    rot_transpose(_a2->_rotation_absolute, tr1);
    rot_mul(_FERMITOF_before->_rotation_absolute, tr1, _FERMITOF_before->_rotation_relative);
    _FERMITOF_before->_rotation_is_identity =  rot_test_identity(_FERMITOF_before->_rotation_relative);
    tc1 = coords_set(
      0, 0, 0.9);
    rot_transpose(_a2->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _FERMITOF_before->_position_absolute = coords_add(_a2->_position_absolute, tc2);
    tc1 = coords_sub(_a2->_position_absolute, _FERMITOF_before->_position_absolute);
    _FERMITOF_before->_position_relative = rot_apply(_FERMITOF_before->_rotation_absolute, tc1);
  } /* FERMITOF_before=Monitor_nD() AT ROTATED */
  DEBUG_COMPONENT("FERMITOF_before", _FERMITOF_before->_position_absolute, _FERMITOF_before->_rotation_absolute);
  instrument->_position_absolute[26] = _FERMITOF_before->_position_absolute;
  instrument->_position_relative[26] = _FERMITOF_before->_position_relative;
  instrument->counter_N[26]  = instrument->counter_P[26] = instrument->counter_P2[26] = 0;
  instrument->counter_AbsorbProp[26]= 0;
  return(0);
} /* _FERMITOF_before_setpos */

/* component lambdaFermi=L_monitor() SETTING, POSITION/ROTATION */
int _lambdaFermi_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_lambdaFermi_setpos] component lambdaFermi=L_monitor() SETTING [/usr/share/mcstas/3.0-dev/monitors/L_monitor.comp:65]");
  stracpy(_lambdaFermi->_name, "lambdaFermi", 16384);
  stracpy(_lambdaFermi->_type, "L_monitor", 16384);
  _lambdaFermi->_index=27;
  _lambdaFermi->_parameters.nL = 60;
  #define nL (_lambdaFermi->_parameters.nL)
  if("lambdafermi.dat" && strlen("lambdafermi.dat"))
    stracpy(_lambdaFermi->_parameters.filename, "lambdafermi.dat" ? "lambdafermi.dat" : "", 16384);
  else 
  _lambdaFermi->_parameters.filename[0]='\0';
  #define filename (_lambdaFermi->_parameters.filename)
  _lambdaFermi->_parameters.xmin = -0.05;
  #define xmin (_lambdaFermi->_parameters.xmin)
  _lambdaFermi->_parameters.xmax = 0.05;
  #define xmax (_lambdaFermi->_parameters.xmax)
  _lambdaFermi->_parameters.ymin = -0.05;
  #define ymin (_lambdaFermi->_parameters.ymin)
  _lambdaFermi->_parameters.ymax = 0.05;
  #define ymax (_lambdaFermi->_parameters.ymax)
  _lambdaFermi->_parameters.xwidth = 0.1;
  #define xwidth (_lambdaFermi->_parameters.xwidth)
  _lambdaFermi->_parameters.yheight = 0.2;
  #define yheight (_lambdaFermi->_parameters.yheight)
  _lambdaFermi->_parameters.Lmin = LMIN;
  #define Lmin (_lambdaFermi->_parameters.Lmin)
  _lambdaFermi->_parameters.Lmax = LMAX;
  #define Lmax (_lambdaFermi->_parameters.Lmax)
  _lambdaFermi->_parameters.restore_neutron = 0;
  #define restore_neutron (_lambdaFermi->_parameters.restore_neutron)

  #define L_N (_lambdaFermi->_parameters.L_N)
  #define L_p (_lambdaFermi->_parameters.L_p)
  #define L_p2 (_lambdaFermi->_parameters.L_p2)

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
  /* component lambdaFermi=L_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _a2->_rotation_absolute, _lambdaFermi->_rotation_absolute);
    rot_transpose(_FERMITOF_before->_rotation_absolute, tr1);
    rot_mul(_lambdaFermi->_rotation_absolute, tr1, _lambdaFermi->_rotation_relative);
    _lambdaFermi->_rotation_is_identity =  rot_test_identity(_lambdaFermi->_rotation_relative);
    tc1 = coords_set(
      0, 0, 0.901);
    rot_transpose(_a2->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _lambdaFermi->_position_absolute = coords_add(_a2->_position_absolute, tc2);
    tc1 = coords_sub(_FERMITOF_before->_position_absolute, _lambdaFermi->_position_absolute);
    _lambdaFermi->_position_relative = rot_apply(_lambdaFermi->_rotation_absolute, tc1);
  } /* lambdaFermi=L_monitor() AT ROTATED */
  DEBUG_COMPONENT("lambdaFermi", _lambdaFermi->_position_absolute, _lambdaFermi->_rotation_absolute);
  instrument->_position_absolute[27] = _lambdaFermi->_position_absolute;
  instrument->_position_relative[27] = _lambdaFermi->_position_relative;
  instrument->counter_N[27]  = instrument->counter_P[27] = instrument->counter_P2[27] = 0;
  instrument->counter_AbsorbProp[27]= 0;
  return(0);
} /* _lambdaFermi_setpos */

/* component EMON_Fermi=E_monitor() SETTING, POSITION/ROTATION */
int _EMON_Fermi_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_EMON_Fermi_setpos] component EMON_Fermi=E_monitor() SETTING [/usr/share/mcstas/3.0-dev/monitors/E_monitor.comp:66]");
  stracpy(_EMON_Fermi->_name, "EMON_Fermi", 16384);
  stracpy(_EMON_Fermi->_type, "E_monitor", 16384);
  _EMON_Fermi->_index=28;
  _EMON_Fermi->_parameters.nE = 60;
  #define nE (_EMON_Fermi->_parameters.nE)
  if("emon_fermi.dat" && strlen("emon_fermi.dat"))
    stracpy(_EMON_Fermi->_parameters.filename, "emon_fermi.dat" ? "emon_fermi.dat" : "", 16384);
  else 
  _EMON_Fermi->_parameters.filename[0]='\0';
  #define filename (_EMON_Fermi->_parameters.filename)
  _EMON_Fermi->_parameters.xmin = -0.05;
  #define xmin (_EMON_Fermi->_parameters.xmin)
  _EMON_Fermi->_parameters.xmax = 0.05;
  #define xmax (_EMON_Fermi->_parameters.xmax)
  _EMON_Fermi->_parameters.ymin = -0.05;
  #define ymin (_EMON_Fermi->_parameters.ymin)
  _EMON_Fermi->_parameters.ymax = 0.05;
  #define ymax (_EMON_Fermi->_parameters.ymax)
  _EMON_Fermi->_parameters.xwidth = 0.06;
  #define xwidth (_EMON_Fermi->_parameters.xwidth)
  _EMON_Fermi->_parameters.yheight = 0.1;
  #define yheight (_EMON_Fermi->_parameters.yheight)
  _EMON_Fermi->_parameters.Emin = EMIN;
  #define Emin (_EMON_Fermi->_parameters.Emin)
  _EMON_Fermi->_parameters.Emax = EMAX;
  #define Emax (_EMON_Fermi->_parameters.Emax)
  _EMON_Fermi->_parameters.restore_neutron = 0;
  #define restore_neutron (_EMON_Fermi->_parameters.restore_neutron)

  #define E_N (_EMON_Fermi->_parameters.E_N)
  #define E_p (_EMON_Fermi->_parameters.E_p)
  #define E_p2 (_EMON_Fermi->_parameters.E_p2)
  #define S_p (_EMON_Fermi->_parameters.S_p)
  #define S_pE (_EMON_Fermi->_parameters.S_pE)
  #define S_pE2 (_EMON_Fermi->_parameters.S_pE2)

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
  /* component EMON_Fermi=E_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (0)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _a2->_rotation_absolute, _EMON_Fermi->_rotation_absolute);
    rot_transpose(_lambdaFermi->_rotation_absolute, tr1);
    rot_mul(_EMON_Fermi->_rotation_absolute, tr1, _EMON_Fermi->_rotation_relative);
    _EMON_Fermi->_rotation_is_identity =  rot_test_identity(_EMON_Fermi->_rotation_relative);
    tc1 = coords_set(
      0, 0, 0.9397);
    rot_transpose(_a2->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _EMON_Fermi->_position_absolute = coords_add(_a2->_position_absolute, tc2);
    tc1 = coords_sub(_lambdaFermi->_position_absolute, _EMON_Fermi->_position_absolute);
    _EMON_Fermi->_position_relative = rot_apply(_EMON_Fermi->_rotation_absolute, tc1);
  } /* EMON_Fermi=E_monitor() AT ROTATED */
  DEBUG_COMPONENT("EMON_Fermi", _EMON_Fermi->_position_absolute, _EMON_Fermi->_rotation_absolute);
  instrument->_position_absolute[28] = _EMON_Fermi->_position_absolute;
  instrument->_position_relative[28] = _EMON_Fermi->_position_relative;
  instrument->counter_N[28]  = instrument->counter_P[28] = instrument->counter_P2[28] = 0;
  instrument->counter_AbsorbProp[28]= 0;
  return(0);
} /* _EMON_Fermi_setpos */

/* component DivMonfermi1=Divergence_monitor() SETTING, POSITION/ROTATION */
int _DivMonfermi1_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_DivMonfermi1_setpos] component DivMonfermi1=Divergence_monitor() SETTING [/usr/share/mcstas/3.0-dev/monitors/Divergence_monitor.comp:71]");
  stracpy(_DivMonfermi1->_name, "DivMonfermi1", 16384);
  stracpy(_DivMonfermi1->_type, "Divergence_monitor", 16384);
  _DivMonfermi1->_index=29;
  _DivMonfermi1->_parameters.nh = 30;
  #define nh (_DivMonfermi1->_parameters.nh)
  _DivMonfermi1->_parameters.nv = 30;
  #define nv (_DivMonfermi1->_parameters.nv)
  if("divfermi1.dat" && strlen("divfermi1.dat"))
    stracpy(_DivMonfermi1->_parameters.filename, "divfermi1.dat" ? "divfermi1.dat" : "", 16384);
  else 
  _DivMonfermi1->_parameters.filename[0]='\0';
  #define filename (_DivMonfermi1->_parameters.filename)
  _DivMonfermi1->_parameters.xmin = -0.05;
  #define xmin (_DivMonfermi1->_parameters.xmin)
  _DivMonfermi1->_parameters.xmax = 0.05;
  #define xmax (_DivMonfermi1->_parameters.xmax)
  _DivMonfermi1->_parameters.ymin = -0.05;
  #define ymin (_DivMonfermi1->_parameters.ymin)
  _DivMonfermi1->_parameters.ymax = 0.05;
  #define ymax (_DivMonfermi1->_parameters.ymax)
  _DivMonfermi1->_parameters.xwidth = 0.06;
  #define xwidth (_DivMonfermi1->_parameters.xwidth)
  _DivMonfermi1->_parameters.yheight = 0.1;
  #define yheight (_DivMonfermi1->_parameters.yheight)
  _DivMonfermi1->_parameters.maxdiv_h = 2;
  #define maxdiv_h (_DivMonfermi1->_parameters.maxdiv_h)
  _DivMonfermi1->_parameters.maxdiv_v = 2;
  #define maxdiv_v (_DivMonfermi1->_parameters.maxdiv_v)
  _DivMonfermi1->_parameters.restore_neutron = 0;
  #define restore_neutron (_DivMonfermi1->_parameters.restore_neutron)
  _DivMonfermi1->_parameters.nx = 0;
  #define nx (_DivMonfermi1->_parameters.nx)
  _DivMonfermi1->_parameters.ny = 0;
  #define ny (_DivMonfermi1->_parameters.ny)
  _DivMonfermi1->_parameters.nz = 1;
  #define nz (_DivMonfermi1->_parameters.nz)

  #define Div_N (_DivMonfermi1->_parameters.Div_N)
  #define Div_p (_DivMonfermi1->_parameters.Div_p)
  #define Div_p2 (_DivMonfermi1->_parameters.Div_p2)

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
  /* component DivMonfermi1=Divergence_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _a2->_rotation_absolute, _DivMonfermi1->_rotation_absolute);
    rot_transpose(_EMON_Fermi->_rotation_absolute, tr1);
    rot_mul(_DivMonfermi1->_rotation_absolute, tr1, _DivMonfermi1->_rotation_relative);
    _DivMonfermi1->_rotation_is_identity =  rot_test_identity(_DivMonfermi1->_rotation_relative);
    tc1 = coords_set(
      0, 0, 0.9398);
    rot_transpose(_a2->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _DivMonfermi1->_position_absolute = coords_add(_a2->_position_absolute, tc2);
    tc1 = coords_sub(_EMON_Fermi->_position_absolute, _DivMonfermi1->_position_absolute);
    _DivMonfermi1->_position_relative = rot_apply(_DivMonfermi1->_rotation_absolute, tc1);
  } /* DivMonfermi1=Divergence_monitor() AT ROTATED */
  DEBUG_COMPONENT("DivMonfermi1", _DivMonfermi1->_position_absolute, _DivMonfermi1->_rotation_absolute);
  instrument->_position_absolute[29] = _DivMonfermi1->_position_absolute;
  instrument->_position_relative[29] = _DivMonfermi1->_position_relative;
  instrument->counter_N[29]  = instrument->counter_P[29] = instrument->counter_P2[29] = 0;
  instrument->counter_AbsorbProp[29]= 0;
  return(0);
} /* _DivMonfermi1_setpos */

/* component PSD_Fermi1=PSD_monitor() SETTING, POSITION/ROTATION */
int _PSD_Fermi1_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_PSD_Fermi1_setpos] component PSD_Fermi1=PSD_monitor() SETTING [/usr/share/mcstas/3.0-dev/monitors/PSD_monitor.comp:64]");
  stracpy(_PSD_Fermi1->_name, "PSD_Fermi1", 16384);
  stracpy(_PSD_Fermi1->_type, "PSD_monitor", 16384);
  _PSD_Fermi1->_index=30;
  _PSD_Fermi1->_parameters.nx = 30;
  #define nx (_PSD_Fermi1->_parameters.nx)
  _PSD_Fermi1->_parameters.ny = 30;
  #define ny (_PSD_Fermi1->_parameters.ny)
  if("psdfermi1.dat" && strlen("psdfermi1.dat"))
    stracpy(_PSD_Fermi1->_parameters.filename, "psdfermi1.dat" ? "psdfermi1.dat" : "", 16384);
  else 
  _PSD_Fermi1->_parameters.filename[0]='\0';
  #define filename (_PSD_Fermi1->_parameters.filename)
  _PSD_Fermi1->_parameters.xmin = -0.05;
  #define xmin (_PSD_Fermi1->_parameters.xmin)
  _PSD_Fermi1->_parameters.xmax = 0.05;
  #define xmax (_PSD_Fermi1->_parameters.xmax)
  _PSD_Fermi1->_parameters.ymin = -0.05;
  #define ymin (_PSD_Fermi1->_parameters.ymin)
  _PSD_Fermi1->_parameters.ymax = 0.05;
  #define ymax (_PSD_Fermi1->_parameters.ymax)
  _PSD_Fermi1->_parameters.xwidth = 0.06;
  #define xwidth (_PSD_Fermi1->_parameters.xwidth)
  _PSD_Fermi1->_parameters.yheight = 0.11;
  #define yheight (_PSD_Fermi1->_parameters.yheight)
  _PSD_Fermi1->_parameters.restore_neutron = 0;
  #define restore_neutron (_PSD_Fermi1->_parameters.restore_neutron)

  #define PSD_N (_PSD_Fermi1->_parameters.PSD_N)
  #define PSD_p (_PSD_Fermi1->_parameters.PSD_p)
  #define PSD_p2 (_PSD_Fermi1->_parameters.PSD_p2)

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
  /* component PSD_Fermi1=PSD_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (0)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _a2->_rotation_absolute, _PSD_Fermi1->_rotation_absolute);
    rot_transpose(_DivMonfermi1->_rotation_absolute, tr1);
    rot_mul(_PSD_Fermi1->_rotation_absolute, tr1, _PSD_Fermi1->_rotation_relative);
    _PSD_Fermi1->_rotation_is_identity =  rot_test_identity(_PSD_Fermi1->_rotation_relative);
    tc1 = coords_set(
      0, 0, 0.9399);
    rot_transpose(_a2->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _PSD_Fermi1->_position_absolute = coords_add(_a2->_position_absolute, tc2);
    tc1 = coords_sub(_DivMonfermi1->_position_absolute, _PSD_Fermi1->_position_absolute);
    _PSD_Fermi1->_position_relative = rot_apply(_PSD_Fermi1->_rotation_absolute, tc1);
  } /* PSD_Fermi1=PSD_monitor() AT ROTATED */
  DEBUG_COMPONENT("PSD_Fermi1", _PSD_Fermi1->_position_absolute, _PSD_Fermi1->_rotation_absolute);
  instrument->_position_absolute[30] = _PSD_Fermi1->_position_absolute;
  instrument->_position_relative[30] = _PSD_Fermi1->_position_relative;
  instrument->counter_N[30]  = instrument->counter_P[30] = instrument->counter_P2[30] = 0;
  instrument->counter_AbsorbProp[30]= 0;
  return(0);
} /* _PSD_Fermi1_setpos */

/* component FoChopper=FermiChopper() SETTING, POSITION/ROTATION */
int _FoChopper_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_FoChopper_setpos] component FoChopper=FermiChopper() SETTING [/usr/share/mcstas/3.0-dev/optics/FermiChopper.comp:298]");
  stracpy(_FoChopper->_name, "FoChopper", 16384);
  stracpy(_FoChopper->_type, "FermiChopper", 16384);
  _FoChopper->_index=31;
  _FoChopper->_parameters.phase = 0;
  #define phase (_FoChopper->_parameters.phase)
  _FoChopper->_parameters.radius = 0.06;
  #define radius (_FoChopper->_parameters.radius)
  _FoChopper->_parameters.nu = - FERMI_SPEED;
  #define nu (_FoChopper->_parameters.nu)
  _FoChopper->_parameters.w = 0.0005;
  #define w (_FoChopper->_parameters.w)
  _FoChopper->_parameters.nslit = 120;
  #define nslit (_FoChopper->_parameters.nslit)
  _FoChopper->_parameters.R0 = 0;
  #define R0 (_FoChopper->_parameters.R0)
  _FoChopper->_parameters.Qc = 0;
  #define Qc (_FoChopper->_parameters.Qc)
  _FoChopper->_parameters.alpha = 0;
  #define alpha (_FoChopper->_parameters.alpha)
  _FoChopper->_parameters.m = 0;
  #define m (_FoChopper->_parameters.m)
  _FoChopper->_parameters.W = 0.0001;
  #define W (_FoChopper->_parameters.W)
  _FoChopper->_parameters.length = 0.012;
  #define length (_FoChopper->_parameters.length)
  _FoChopper->_parameters.eff = 0.95;
  #define eff (_FoChopper->_parameters.eff)
  _FoChopper->_parameters.zero_time = 0;
  #define zero_time (_FoChopper->_parameters.zero_time)
  _FoChopper->_parameters.xwidth = 0;
  #define xwidth (_FoChopper->_parameters.xwidth)
  _FoChopper->_parameters.verbose = 0;
  #define verbose (_FoChopper->_parameters.verbose)
  _FoChopper->_parameters.yheight = 0.11;
  #define yheight (_FoChopper->_parameters.yheight)
  _FoChopper->_parameters.curvature = 0;
  #define curvature (_FoChopper->_parameters.curvature)
  _FoChopper->_parameters.delay = 0;
  #define delay (_FoChopper->_parameters.delay)

  #define FCVars (_FoChopper->_parameters.FCVars)

  #undef phase
  #undef radius
  #undef nu
  #undef w
  #undef nslit
  #undef R0
  #undef Qc
  #undef alpha
  #undef m
  #undef W
  #undef length
  #undef eff
  #undef zero_time
  #undef xwidth
  #undef verbose
  #undef yheight
  #undef curvature
  #undef delay
  #undef FCVars
  /* component FoChopper=FermiChopper() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _a2->_rotation_absolute, _FoChopper->_rotation_absolute);
    rot_transpose(_PSD_Fermi1->_rotation_absolute, tr1);
    rot_mul(_FoChopper->_rotation_absolute, tr1, _FoChopper->_rotation_relative);
    _FoChopper->_rotation_is_identity =  rot_test_identity(_FoChopper->_rotation_relative);
    tc1 = coords_set(
      0, 0, 1.002);
    rot_transpose(_a2->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _FoChopper->_position_absolute = coords_add(_a2->_position_absolute, tc2);
    tc1 = coords_sub(_PSD_Fermi1->_position_absolute, _FoChopper->_position_absolute);
    _FoChopper->_position_relative = rot_apply(_FoChopper->_rotation_absolute, tc1);
  } /* FoChopper=FermiChopper() AT ROTATED */
  DEBUG_COMPONENT("FoChopper", _FoChopper->_position_absolute, _FoChopper->_rotation_absolute);
  instrument->_position_absolute[31] = _FoChopper->_position_absolute;
  instrument->_position_relative[31] = _FoChopper->_position_relative;
  instrument->counter_N[31]  = instrument->counter_P[31] = instrument->counter_P2[31] = 0;
  instrument->counter_AbsorbProp[31]= 0;
  return(0);
} /* _FoChopper_setpos */

/* component PSD_Fermi2=PSD_monitor() SETTING, POSITION/ROTATION */
int _PSD_Fermi2_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_PSD_Fermi2_setpos] component PSD_Fermi2=PSD_monitor() SETTING [/usr/share/mcstas/3.0-dev/monitors/PSD_monitor.comp:64]");
  stracpy(_PSD_Fermi2->_name, "PSD_Fermi2", 16384);
  stracpy(_PSD_Fermi2->_type, "PSD_monitor", 16384);
  _PSD_Fermi2->_index=32;
  _PSD_Fermi2->_parameters.nx = 30;
  #define nx (_PSD_Fermi2->_parameters.nx)
  _PSD_Fermi2->_parameters.ny = 30;
  #define ny (_PSD_Fermi2->_parameters.ny)
  if("psdfermi2.dat" && strlen("psdfermi2.dat"))
    stracpy(_PSD_Fermi2->_parameters.filename, "psdfermi2.dat" ? "psdfermi2.dat" : "", 16384);
  else 
  _PSD_Fermi2->_parameters.filename[0]='\0';
  #define filename (_PSD_Fermi2->_parameters.filename)
  _PSD_Fermi2->_parameters.xmin = -0.05;
  #define xmin (_PSD_Fermi2->_parameters.xmin)
  _PSD_Fermi2->_parameters.xmax = 0.05;
  #define xmax (_PSD_Fermi2->_parameters.xmax)
  _PSD_Fermi2->_parameters.ymin = -0.05;
  #define ymin (_PSD_Fermi2->_parameters.ymin)
  _PSD_Fermi2->_parameters.ymax = 0.05;
  #define ymax (_PSD_Fermi2->_parameters.ymax)
  _PSD_Fermi2->_parameters.xwidth = 0.06;
  #define xwidth (_PSD_Fermi2->_parameters.xwidth)
  _PSD_Fermi2->_parameters.yheight = 0.11;
  #define yheight (_PSD_Fermi2->_parameters.yheight)
  _PSD_Fermi2->_parameters.restore_neutron = 0;
  #define restore_neutron (_PSD_Fermi2->_parameters.restore_neutron)

  #define PSD_N (_PSD_Fermi2->_parameters.PSD_N)
  #define PSD_p (_PSD_Fermi2->_parameters.PSD_p)
  #define PSD_p2 (_PSD_Fermi2->_parameters.PSD_p2)

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
  /* component PSD_Fermi2=PSD_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (0)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _a2->_rotation_absolute, _PSD_Fermi2->_rotation_absolute);
    rot_transpose(_FoChopper->_rotation_absolute, tr1);
    rot_mul(_PSD_Fermi2->_rotation_absolute, tr1, _PSD_Fermi2->_rotation_relative);
    _PSD_Fermi2->_rotation_is_identity =  rot_test_identity(_PSD_Fermi2->_rotation_relative);
    tc1 = coords_set(
      0, 0, 1.063);
    rot_transpose(_a2->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _PSD_Fermi2->_position_absolute = coords_add(_a2->_position_absolute, tc2);
    tc1 = coords_sub(_FoChopper->_position_absolute, _PSD_Fermi2->_position_absolute);
    _PSD_Fermi2->_position_relative = rot_apply(_PSD_Fermi2->_rotation_absolute, tc1);
  } /* PSD_Fermi2=PSD_monitor() AT ROTATED */
  DEBUG_COMPONENT("PSD_Fermi2", _PSD_Fermi2->_position_absolute, _PSD_Fermi2->_rotation_absolute);
  instrument->_position_absolute[32] = _PSD_Fermi2->_position_absolute;
  instrument->_position_relative[32] = _PSD_Fermi2->_position_relative;
  instrument->counter_N[32]  = instrument->counter_P[32] = instrument->counter_P2[32] = 0;
  instrument->counter_AbsorbProp[32]= 0;
  return(0);
} /* _PSD_Fermi2_setpos */

/* component DivMonfermi2=Divergence_monitor() SETTING, POSITION/ROTATION */
int _DivMonfermi2_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_DivMonfermi2_setpos] component DivMonfermi2=Divergence_monitor() SETTING [/usr/share/mcstas/3.0-dev/monitors/Divergence_monitor.comp:71]");
  stracpy(_DivMonfermi2->_name, "DivMonfermi2", 16384);
  stracpy(_DivMonfermi2->_type, "Divergence_monitor", 16384);
  _DivMonfermi2->_index=33;
  _DivMonfermi2->_parameters.nh = 30;
  #define nh (_DivMonfermi2->_parameters.nh)
  _DivMonfermi2->_parameters.nv = 30;
  #define nv (_DivMonfermi2->_parameters.nv)
  if("divfermi2.dat" && strlen("divfermi2.dat"))
    stracpy(_DivMonfermi2->_parameters.filename, "divfermi2.dat" ? "divfermi2.dat" : "", 16384);
  else 
  _DivMonfermi2->_parameters.filename[0]='\0';
  #define filename (_DivMonfermi2->_parameters.filename)
  _DivMonfermi2->_parameters.xmin = -0.05;
  #define xmin (_DivMonfermi2->_parameters.xmin)
  _DivMonfermi2->_parameters.xmax = 0.05;
  #define xmax (_DivMonfermi2->_parameters.xmax)
  _DivMonfermi2->_parameters.ymin = -0.05;
  #define ymin (_DivMonfermi2->_parameters.ymin)
  _DivMonfermi2->_parameters.ymax = 0.05;
  #define ymax (_DivMonfermi2->_parameters.ymax)
  _DivMonfermi2->_parameters.xwidth = 0.06;
  #define xwidth (_DivMonfermi2->_parameters.xwidth)
  _DivMonfermi2->_parameters.yheight = 0.1;
  #define yheight (_DivMonfermi2->_parameters.yheight)
  _DivMonfermi2->_parameters.maxdiv_h = 2;
  #define maxdiv_h (_DivMonfermi2->_parameters.maxdiv_h)
  _DivMonfermi2->_parameters.maxdiv_v = 2;
  #define maxdiv_v (_DivMonfermi2->_parameters.maxdiv_v)
  _DivMonfermi2->_parameters.restore_neutron = 0;
  #define restore_neutron (_DivMonfermi2->_parameters.restore_neutron)
  _DivMonfermi2->_parameters.nx = 0;
  #define nx (_DivMonfermi2->_parameters.nx)
  _DivMonfermi2->_parameters.ny = 0;
  #define ny (_DivMonfermi2->_parameters.ny)
  _DivMonfermi2->_parameters.nz = 1;
  #define nz (_DivMonfermi2->_parameters.nz)

  #define Div_N (_DivMonfermi2->_parameters.Div_N)
  #define Div_p (_DivMonfermi2->_parameters.Div_p)
  #define Div_p2 (_DivMonfermi2->_parameters.Div_p2)

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
  /* component DivMonfermi2=Divergence_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _a2->_rotation_absolute, _DivMonfermi2->_rotation_absolute);
    rot_transpose(_PSD_Fermi2->_rotation_absolute, tr1);
    rot_mul(_DivMonfermi2->_rotation_absolute, tr1, _DivMonfermi2->_rotation_relative);
    _DivMonfermi2->_rotation_is_identity =  rot_test_identity(_DivMonfermi2->_rotation_relative);
    tc1 = coords_set(
      0, 0, 1.064);
    rot_transpose(_a2->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _DivMonfermi2->_position_absolute = coords_add(_a2->_position_absolute, tc2);
    tc1 = coords_sub(_PSD_Fermi2->_position_absolute, _DivMonfermi2->_position_absolute);
    _DivMonfermi2->_position_relative = rot_apply(_DivMonfermi2->_rotation_absolute, tc1);
  } /* DivMonfermi2=Divergence_monitor() AT ROTATED */
  DEBUG_COMPONENT("DivMonfermi2", _DivMonfermi2->_position_absolute, _DivMonfermi2->_rotation_absolute);
  instrument->_position_absolute[33] = _DivMonfermi2->_position_absolute;
  instrument->_position_relative[33] = _DivMonfermi2->_position_relative;
  instrument->counter_N[33]  = instrument->counter_P[33] = instrument->counter_P2[33] = 0;
  instrument->counter_AbsorbProp[33]= 0;
  return(0);
} /* _DivMonfermi2_setpos */

/* component FERMITOF1=Monitor_nD() SETTING, POSITION/ROTATION */
int _FERMITOF1_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_FERMITOF1_setpos] component FERMITOF1=Monitor_nD() SETTING [/usr/share/mcstas/3.0-dev/monitors/Monitor_nD.comp:228]");
  stracpy(_FERMITOF1->_name, "FERMITOF1", 16384);
  stracpy(_FERMITOF1->_type, "Monitor_nD", 16384);
  _FERMITOF1->_index=34;
  _FERMITOF1->_parameters.user1 = '\0';
  #define user1 (_FERMITOF1->_parameters.user1)
  _FERMITOF1->_parameters.user2 = '\0';
  #define user2 (_FERMITOF1->_parameters.user2)
  _FERMITOF1->_parameters.user3 = '\0';
  #define user3 (_FERMITOF1->_parameters.user3)
  _FERMITOF1->_parameters.xwidth = 0;
  #define xwidth (_FERMITOF1->_parameters.xwidth)
  _FERMITOF1->_parameters.yheight = 0;
  #define yheight (_FERMITOF1->_parameters.yheight)
  _FERMITOF1->_parameters.zdepth = 0;
  #define zdepth (_FERMITOF1->_parameters.zdepth)
  _FERMITOF1->_parameters.xmin = -0.05;
  #define xmin (_FERMITOF1->_parameters.xmin)
  _FERMITOF1->_parameters.xmax = 0.05;
  #define xmax (_FERMITOF1->_parameters.xmax)
  _FERMITOF1->_parameters.ymin = -0.1;
  #define ymin (_FERMITOF1->_parameters.ymin)
  _FERMITOF1->_parameters.ymax = 0.1;
  #define ymax (_FERMITOF1->_parameters.ymax)
  _FERMITOF1->_parameters.zmin = 0;
  #define zmin (_FERMITOF1->_parameters.zmin)
  _FERMITOF1->_parameters.zmax = 0;
  #define zmax (_FERMITOF1->_parameters.zmax)
  _FERMITOF1->_parameters.bins = 30;
  #define bins (_FERMITOF1->_parameters.bins)
  _FERMITOF1->_parameters.min = -1e40;
  #define min (_FERMITOF1->_parameters.min)
  _FERMITOF1->_parameters.max = 1e40;
  #define max (_FERMITOF1->_parameters.max)
  _FERMITOF1->_parameters.restore_neutron = 0;
  #define restore_neutron (_FERMITOF1->_parameters.restore_neutron)
  _FERMITOF1->_parameters.radius = 0;
  #define radius (_FERMITOF1->_parameters.radius)
  if("auto, time" && strlen("auto, time"))
    stracpy(_FERMITOF1->_parameters.options, "auto, time" ? "auto, time" : "", 16384);
  else 
  _FERMITOF1->_parameters.options[0]='\0';
  #define options (_FERMITOF1->_parameters.options)
  if("FERMI_TOF1.dat.dat" && strlen("FERMI_TOF1.dat.dat"))
    stracpy(_FERMITOF1->_parameters.filename, "FERMI_TOF1.dat.dat" ? "FERMI_TOF1.dat.dat" : "", 16384);
  else 
  _FERMITOF1->_parameters.filename[0]='\0';
  #define filename (_FERMITOF1->_parameters.filename)
  if("NULL" && strlen("NULL"))
    stracpy(_FERMITOF1->_parameters.geometry, "NULL" ? "NULL" : "", 16384);
  else 
  _FERMITOF1->_parameters.geometry[0]='\0';
  #define geometry (_FERMITOF1->_parameters.geometry)
  if("NULL" && strlen("NULL"))
    stracpy(_FERMITOF1->_parameters.username1, "NULL" ? "NULL" : "", 16384);
  else 
  _FERMITOF1->_parameters.username1[0]='\0';
  #define username1 (_FERMITOF1->_parameters.username1)
  if("NULL" && strlen("NULL"))
    stracpy(_FERMITOF1->_parameters.username2, "NULL" ? "NULL" : "", 16384);
  else 
  _FERMITOF1->_parameters.username2[0]='\0';
  #define username2 (_FERMITOF1->_parameters.username2)
  if("NULL" && strlen("NULL"))
    stracpy(_FERMITOF1->_parameters.username3, "NULL" ? "NULL" : "", 16384);
  else 
  _FERMITOF1->_parameters.username3[0]='\0';
  #define username3 (_FERMITOF1->_parameters.username3)

  #define DEFS (_FERMITOF1->_parameters.DEFS)
  #define Vars (_FERMITOF1->_parameters.Vars)
  #define detector (_FERMITOF1->_parameters.detector)
  #define offdata (_FERMITOF1->_parameters.offdata)

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
  /* component FERMITOF1=Monitor_nD() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _a2->_rotation_absolute, _FERMITOF1->_rotation_absolute);
    rot_transpose(_DivMonfermi2->_rotation_absolute, tr1);
    rot_mul(_FERMITOF1->_rotation_absolute, tr1, _FERMITOF1->_rotation_relative);
    _FERMITOF1->_rotation_is_identity =  rot_test_identity(_FERMITOF1->_rotation_relative);
    tc1 = coords_set(
      0, 0, 1.154);
    rot_transpose(_a2->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _FERMITOF1->_position_absolute = coords_add(_a2->_position_absolute, tc2);
    tc1 = coords_sub(_DivMonfermi2->_position_absolute, _FERMITOF1->_position_absolute);
    _FERMITOF1->_position_relative = rot_apply(_FERMITOF1->_rotation_absolute, tc1);
  } /* FERMITOF1=Monitor_nD() AT ROTATED */
  DEBUG_COMPONENT("FERMITOF1", _FERMITOF1->_position_absolute, _FERMITOF1->_rotation_absolute);
  instrument->_position_absolute[34] = _FERMITOF1->_position_absolute;
  instrument->_position_relative[34] = _FERMITOF1->_position_relative;
  instrument->counter_N[34]  = instrument->counter_P[34] = instrument->counter_P2[34] = 0;
  instrument->counter_AbsorbProp[34]= 0;
  return(0);
} /* _FERMITOF1_setpos */

/* component SAMPLE_SLIT=Slit() SETTING, POSITION/ROTATION */
int _SAMPLE_SLIT_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_SAMPLE_SLIT_setpos] component SAMPLE_SLIT=Slit() SETTING [/usr/share/mcstas/3.0-dev/optics/Slit.comp:47]");
  stracpy(_SAMPLE_SLIT->_name, "SAMPLE_SLIT", 16384);
  stracpy(_SAMPLE_SLIT->_type, "Slit", 16384);
  _SAMPLE_SLIT->_index=35;
  _SAMPLE_SLIT->_parameters.xmin = -0.01;
  #define xmin (_SAMPLE_SLIT->_parameters.xmin)
  _SAMPLE_SLIT->_parameters.xmax = 0.01;
  #define xmax (_SAMPLE_SLIT->_parameters.xmax)
  _SAMPLE_SLIT->_parameters.ymin = -0.01;
  #define ymin (_SAMPLE_SLIT->_parameters.ymin)
  _SAMPLE_SLIT->_parameters.ymax = 0.01;
  #define ymax (_SAMPLE_SLIT->_parameters.ymax)
  _SAMPLE_SLIT->_parameters.radius = 0;
  #define radius (_SAMPLE_SLIT->_parameters.radius)
  _SAMPLE_SLIT->_parameters.xwidth = 0.02;
  #define xwidth (_SAMPLE_SLIT->_parameters.xwidth)
  _SAMPLE_SLIT->_parameters.yheight = 0.06;
  #define yheight (_SAMPLE_SLIT->_parameters.yheight)

  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef radius
  #undef xwidth
  #undef yheight
  /* component SAMPLE_SLIT=Slit() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _a2->_rotation_absolute, _SAMPLE_SLIT->_rotation_absolute);
    rot_transpose(_FERMITOF1->_rotation_absolute, tr1);
    rot_mul(_SAMPLE_SLIT->_rotation_absolute, tr1, _SAMPLE_SLIT->_rotation_relative);
    _SAMPLE_SLIT->_rotation_is_identity =  rot_test_identity(_SAMPLE_SLIT->_rotation_relative);
    tc1 = coords_set(
      0, 0, 1.155);
    rot_transpose(_a2->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _SAMPLE_SLIT->_position_absolute = coords_add(_a2->_position_absolute, tc2);
    tc1 = coords_sub(_FERMITOF1->_position_absolute, _SAMPLE_SLIT->_position_absolute);
    _SAMPLE_SLIT->_position_relative = rot_apply(_SAMPLE_SLIT->_rotation_absolute, tc1);
  } /* SAMPLE_SLIT=Slit() AT ROTATED */
  DEBUG_COMPONENT("SAMPLE_SLIT", _SAMPLE_SLIT->_position_absolute, _SAMPLE_SLIT->_rotation_absolute);
  instrument->_position_absolute[35] = _SAMPLE_SLIT->_position_absolute;
  instrument->_position_relative[35] = _SAMPLE_SLIT->_position_relative;
  instrument->counter_N[35]  = instrument->counter_P[35] = instrument->counter_P2[35] = 0;
  instrument->counter_AbsorbProp[35]= 0;
  return(0);
} /* _SAMPLE_SLIT_setpos */

/* component FERMITOF2=Monitor_nD() SETTING, POSITION/ROTATION */
int _FERMITOF2_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_FERMITOF2_setpos] component FERMITOF2=Monitor_nD() SETTING [/usr/share/mcstas/3.0-dev/monitors/Monitor_nD.comp:228]");
  stracpy(_FERMITOF2->_name, "FERMITOF2", 16384);
  stracpy(_FERMITOF2->_type, "Monitor_nD", 16384);
  _FERMITOF2->_index=36;
  _FERMITOF2->_parameters.user1 = '\0';
  #define user1 (_FERMITOF2->_parameters.user1)
  _FERMITOF2->_parameters.user2 = '\0';
  #define user2 (_FERMITOF2->_parameters.user2)
  _FERMITOF2->_parameters.user3 = '\0';
  #define user3 (_FERMITOF2->_parameters.user3)
  _FERMITOF2->_parameters.xwidth = 0;
  #define xwidth (_FERMITOF2->_parameters.xwidth)
  _FERMITOF2->_parameters.yheight = 0;
  #define yheight (_FERMITOF2->_parameters.yheight)
  _FERMITOF2->_parameters.zdepth = 0;
  #define zdepth (_FERMITOF2->_parameters.zdepth)
  _FERMITOF2->_parameters.xmin = -0.05;
  #define xmin (_FERMITOF2->_parameters.xmin)
  _FERMITOF2->_parameters.xmax = 0.05;
  #define xmax (_FERMITOF2->_parameters.xmax)
  _FERMITOF2->_parameters.ymin = -0.1;
  #define ymin (_FERMITOF2->_parameters.ymin)
  _FERMITOF2->_parameters.ymax = 0.1;
  #define ymax (_FERMITOF2->_parameters.ymax)
  _FERMITOF2->_parameters.zmin = 0;
  #define zmin (_FERMITOF2->_parameters.zmin)
  _FERMITOF2->_parameters.zmax = 0;
  #define zmax (_FERMITOF2->_parameters.zmax)
  _FERMITOF2->_parameters.bins = 30;
  #define bins (_FERMITOF2->_parameters.bins)
  _FERMITOF2->_parameters.min = -1e40;
  #define min (_FERMITOF2->_parameters.min)
  _FERMITOF2->_parameters.max = 1e40;
  #define max (_FERMITOF2->_parameters.max)
  _FERMITOF2->_parameters.restore_neutron = 0;
  #define restore_neutron (_FERMITOF2->_parameters.restore_neutron)
  _FERMITOF2->_parameters.radius = 0;
  #define radius (_FERMITOF2->_parameters.radius)
  if("auto, time" && strlen("auto, time"))
    stracpy(_FERMITOF2->_parameters.options, "auto, time" ? "auto, time" : "", 16384);
  else 
  _FERMITOF2->_parameters.options[0]='\0';
  #define options (_FERMITOF2->_parameters.options)
  if("FERMI_TOF2.dat.dat" && strlen("FERMI_TOF2.dat.dat"))
    stracpy(_FERMITOF2->_parameters.filename, "FERMI_TOF2.dat.dat" ? "FERMI_TOF2.dat.dat" : "", 16384);
  else 
  _FERMITOF2->_parameters.filename[0]='\0';
  #define filename (_FERMITOF2->_parameters.filename)
  if("NULL" && strlen("NULL"))
    stracpy(_FERMITOF2->_parameters.geometry, "NULL" ? "NULL" : "", 16384);
  else 
  _FERMITOF2->_parameters.geometry[0]='\0';
  #define geometry (_FERMITOF2->_parameters.geometry)
  if("NULL" && strlen("NULL"))
    stracpy(_FERMITOF2->_parameters.username1, "NULL" ? "NULL" : "", 16384);
  else 
  _FERMITOF2->_parameters.username1[0]='\0';
  #define username1 (_FERMITOF2->_parameters.username1)
  if("NULL" && strlen("NULL"))
    stracpy(_FERMITOF2->_parameters.username2, "NULL" ? "NULL" : "", 16384);
  else 
  _FERMITOF2->_parameters.username2[0]='\0';
  #define username2 (_FERMITOF2->_parameters.username2)
  if("NULL" && strlen("NULL"))
    stracpy(_FERMITOF2->_parameters.username3, "NULL" ? "NULL" : "", 16384);
  else 
  _FERMITOF2->_parameters.username3[0]='\0';
  #define username3 (_FERMITOF2->_parameters.username3)

  #define DEFS (_FERMITOF2->_parameters.DEFS)
  #define Vars (_FERMITOF2->_parameters.Vars)
  #define detector (_FERMITOF2->_parameters.detector)
  #define offdata (_FERMITOF2->_parameters.offdata)

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
  /* component FERMITOF2=Monitor_nD() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _a2->_rotation_absolute, _FERMITOF2->_rotation_absolute);
    rot_transpose(_SAMPLE_SLIT->_rotation_absolute, tr1);
    rot_mul(_FERMITOF2->_rotation_absolute, tr1, _FERMITOF2->_rotation_relative);
    _FERMITOF2->_rotation_is_identity =  rot_test_identity(_FERMITOF2->_rotation_relative);
    tc1 = coords_set(
      0, 0, 1.2);
    rot_transpose(_a2->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _FERMITOF2->_position_absolute = coords_add(_a2->_position_absolute, tc2);
    tc1 = coords_sub(_SAMPLE_SLIT->_position_absolute, _FERMITOF2->_position_absolute);
    _FERMITOF2->_position_relative = rot_apply(_FERMITOF2->_rotation_absolute, tc1);
  } /* FERMITOF2=Monitor_nD() AT ROTATED */
  DEBUG_COMPONENT("FERMITOF2", _FERMITOF2->_position_absolute, _FERMITOF2->_rotation_absolute);
  instrument->_position_absolute[36] = _FERMITOF2->_position_absolute;
  instrument->_position_relative[36] = _FERMITOF2->_position_relative;
  instrument->counter_N[36]  = instrument->counter_P[36] = instrument->counter_P2[36] = 0;
  instrument->counter_AbsorbProp[36]= 0;
  return(0);
} /* _FERMITOF2_setpos */

/* component PSD_SAMPLE=PSD_monitor() SETTING, POSITION/ROTATION */
int _PSD_SAMPLE_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_PSD_SAMPLE_setpos] component PSD_SAMPLE=PSD_monitor() SETTING [/usr/share/mcstas/3.0-dev/monitors/PSD_monitor.comp:64]");
  stracpy(_PSD_SAMPLE->_name, "PSD_SAMPLE", 16384);
  stracpy(_PSD_SAMPLE->_type, "PSD_monitor", 16384);
  _PSD_SAMPLE->_index=37;
  _PSD_SAMPLE->_parameters.nx = 30;
  #define nx (_PSD_SAMPLE->_parameters.nx)
  _PSD_SAMPLE->_parameters.ny = 30;
  #define ny (_PSD_SAMPLE->_parameters.ny)
  if("psdsample.dat" && strlen("psdsample.dat"))
    stracpy(_PSD_SAMPLE->_parameters.filename, "psdsample.dat" ? "psdsample.dat" : "", 16384);
  else 
  _PSD_SAMPLE->_parameters.filename[0]='\0';
  #define filename (_PSD_SAMPLE->_parameters.filename)
  _PSD_SAMPLE->_parameters.xmin = -0.05;
  #define xmin (_PSD_SAMPLE->_parameters.xmin)
  _PSD_SAMPLE->_parameters.xmax = 0.05;
  #define xmax (_PSD_SAMPLE->_parameters.xmax)
  _PSD_SAMPLE->_parameters.ymin = -0.05;
  #define ymin (_PSD_SAMPLE->_parameters.ymin)
  _PSD_SAMPLE->_parameters.ymax = 0.05;
  #define ymax (_PSD_SAMPLE->_parameters.ymax)
  _PSD_SAMPLE->_parameters.xwidth = 0.02;
  #define xwidth (_PSD_SAMPLE->_parameters.xwidth)
  _PSD_SAMPLE->_parameters.yheight = 0.1;
  #define yheight (_PSD_SAMPLE->_parameters.yheight)
  _PSD_SAMPLE->_parameters.restore_neutron = 0;
  #define restore_neutron (_PSD_SAMPLE->_parameters.restore_neutron)

  #define PSD_N (_PSD_SAMPLE->_parameters.PSD_N)
  #define PSD_p (_PSD_SAMPLE->_parameters.PSD_p)
  #define PSD_p2 (_PSD_SAMPLE->_parameters.PSD_p2)

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
  /* component PSD_SAMPLE=PSD_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (0)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _a2->_rotation_absolute, _PSD_SAMPLE->_rotation_absolute);
    rot_transpose(_FERMITOF2->_rotation_absolute, tr1);
    rot_mul(_PSD_SAMPLE->_rotation_absolute, tr1, _PSD_SAMPLE->_rotation_relative);
    _PSD_SAMPLE->_rotation_is_identity =  rot_test_identity(_PSD_SAMPLE->_rotation_relative);
    tc1 = coords_set(
      0, 0, 1.45);
    rot_transpose(_a2->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _PSD_SAMPLE->_position_absolute = coords_add(_a2->_position_absolute, tc2);
    tc1 = coords_sub(_FERMITOF2->_position_absolute, _PSD_SAMPLE->_position_absolute);
    _PSD_SAMPLE->_position_relative = rot_apply(_PSD_SAMPLE->_rotation_absolute, tc1);
  } /* PSD_SAMPLE=PSD_monitor() AT ROTATED */
  DEBUG_COMPONENT("PSD_SAMPLE", _PSD_SAMPLE->_position_absolute, _PSD_SAMPLE->_rotation_absolute);
  instrument->_position_absolute[37] = _PSD_SAMPLE->_position_absolute;
  instrument->_position_relative[37] = _PSD_SAMPLE->_position_relative;
  instrument->counter_N[37]  = instrument->counter_P[37] = instrument->counter_P2[37] = 0;
  instrument->counter_AbsorbProp[37]= 0;
  return(0);
} /* _PSD_SAMPLE_setpos */

/* component DivMon_Sample=Divergence_monitor() SETTING, POSITION/ROTATION */
int _DivMon_Sample_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_DivMon_Sample_setpos] component DivMon_Sample=Divergence_monitor() SETTING [/usr/share/mcstas/3.0-dev/monitors/Divergence_monitor.comp:71]");
  stracpy(_DivMon_Sample->_name, "DivMon_Sample", 16384);
  stracpy(_DivMon_Sample->_type, "Divergence_monitor", 16384);
  _DivMon_Sample->_index=38;
  _DivMon_Sample->_parameters.nh = 30;
  #define nh (_DivMon_Sample->_parameters.nh)
  _DivMon_Sample->_parameters.nv = 30;
  #define nv (_DivMon_Sample->_parameters.nv)
  if("div2.dat" && strlen("div2.dat"))
    stracpy(_DivMon_Sample->_parameters.filename, "div2.dat" ? "div2.dat" : "", 16384);
  else 
  _DivMon_Sample->_parameters.filename[0]='\0';
  #define filename (_DivMon_Sample->_parameters.filename)
  _DivMon_Sample->_parameters.xmin = -0.05;
  #define xmin (_DivMon_Sample->_parameters.xmin)
  _DivMon_Sample->_parameters.xmax = 0.05;
  #define xmax (_DivMon_Sample->_parameters.xmax)
  _DivMon_Sample->_parameters.ymin = -0.05;
  #define ymin (_DivMon_Sample->_parameters.ymin)
  _DivMon_Sample->_parameters.ymax = 0.05;
  #define ymax (_DivMon_Sample->_parameters.ymax)
  _DivMon_Sample->_parameters.xwidth = 0.01;
  #define xwidth (_DivMon_Sample->_parameters.xwidth)
  _DivMon_Sample->_parameters.yheight = 0.06;
  #define yheight (_DivMon_Sample->_parameters.yheight)
  _DivMon_Sample->_parameters.maxdiv_h = 3;
  #define maxdiv_h (_DivMon_Sample->_parameters.maxdiv_h)
  _DivMon_Sample->_parameters.maxdiv_v = 3;
  #define maxdiv_v (_DivMon_Sample->_parameters.maxdiv_v)
  _DivMon_Sample->_parameters.restore_neutron = 0;
  #define restore_neutron (_DivMon_Sample->_parameters.restore_neutron)
  _DivMon_Sample->_parameters.nx = 0;
  #define nx (_DivMon_Sample->_parameters.nx)
  _DivMon_Sample->_parameters.ny = 0;
  #define ny (_DivMon_Sample->_parameters.ny)
  _DivMon_Sample->_parameters.nz = 1;
  #define nz (_DivMon_Sample->_parameters.nz)

  #define Div_N (_DivMon_Sample->_parameters.Div_N)
  #define Div_p (_DivMon_Sample->_parameters.Div_p)
  #define Div_p2 (_DivMon_Sample->_parameters.Div_p2)

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
  /* component DivMon_Sample=Divergence_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (0)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _a2->_rotation_absolute, _DivMon_Sample->_rotation_absolute);
    rot_transpose(_PSD_SAMPLE->_rotation_absolute, tr1);
    rot_mul(_DivMon_Sample->_rotation_absolute, tr1, _DivMon_Sample->_rotation_relative);
    _DivMon_Sample->_rotation_is_identity =  rot_test_identity(_DivMon_Sample->_rotation_relative);
    tc1 = coords_set(
      0, 0, 1.451);
    rot_transpose(_a2->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _DivMon_Sample->_position_absolute = coords_add(_a2->_position_absolute, tc2);
    tc1 = coords_sub(_PSD_SAMPLE->_position_absolute, _DivMon_Sample->_position_absolute);
    _DivMon_Sample->_position_relative = rot_apply(_DivMon_Sample->_rotation_absolute, tc1);
  } /* DivMon_Sample=Divergence_monitor() AT ROTATED */
  DEBUG_COMPONENT("DivMon_Sample", _DivMon_Sample->_position_absolute, _DivMon_Sample->_rotation_absolute);
  instrument->_position_absolute[38] = _DivMon_Sample->_position_absolute;
  instrument->_position_relative[38] = _DivMon_Sample->_position_relative;
  instrument->counter_N[38]  = instrument->counter_P[38] = instrument->counter_P2[38] = 0;
  instrument->counter_AbsorbProp[38]= 0;
  return(0);
} /* _DivMon_Sample_setpos */

/* component EMON_SAMPLE=E_monitor() SETTING, POSITION/ROTATION */
int _EMON_SAMPLE_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_EMON_SAMPLE_setpos] component EMON_SAMPLE=E_monitor() SETTING [/usr/share/mcstas/3.0-dev/monitors/E_monitor.comp:66]");
  stracpy(_EMON_SAMPLE->_name, "EMON_SAMPLE", 16384);
  stracpy(_EMON_SAMPLE->_type, "E_monitor", 16384);
  _EMON_SAMPLE->_index=39;
  _EMON_SAMPLE->_parameters.nE = 60;
  #define nE (_EMON_SAMPLE->_parameters.nE)
  if("emon_sample.dat" && strlen("emon_sample.dat"))
    stracpy(_EMON_SAMPLE->_parameters.filename, "emon_sample.dat" ? "emon_sample.dat" : "", 16384);
  else 
  _EMON_SAMPLE->_parameters.filename[0]='\0';
  #define filename (_EMON_SAMPLE->_parameters.filename)
  _EMON_SAMPLE->_parameters.xmin = -0.05;
  #define xmin (_EMON_SAMPLE->_parameters.xmin)
  _EMON_SAMPLE->_parameters.xmax = 0.05;
  #define xmax (_EMON_SAMPLE->_parameters.xmax)
  _EMON_SAMPLE->_parameters.ymin = -0.05;
  #define ymin (_EMON_SAMPLE->_parameters.ymin)
  _EMON_SAMPLE->_parameters.ymax = 0.05;
  #define ymax (_EMON_SAMPLE->_parameters.ymax)
  _EMON_SAMPLE->_parameters.xwidth = 0.01;
  #define xwidth (_EMON_SAMPLE->_parameters.xwidth)
  _EMON_SAMPLE->_parameters.yheight = 0.06;
  #define yheight (_EMON_SAMPLE->_parameters.yheight)
  _EMON_SAMPLE->_parameters.Emin = EMIN;
  #define Emin (_EMON_SAMPLE->_parameters.Emin)
  _EMON_SAMPLE->_parameters.Emax = EMAX;
  #define Emax (_EMON_SAMPLE->_parameters.Emax)
  _EMON_SAMPLE->_parameters.restore_neutron = 0;
  #define restore_neutron (_EMON_SAMPLE->_parameters.restore_neutron)

  #define E_N (_EMON_SAMPLE->_parameters.E_N)
  #define E_p (_EMON_SAMPLE->_parameters.E_p)
  #define E_p2 (_EMON_SAMPLE->_parameters.E_p2)
  #define S_p (_EMON_SAMPLE->_parameters.S_p)
  #define S_pE (_EMON_SAMPLE->_parameters.S_pE)
  #define S_pE2 (_EMON_SAMPLE->_parameters.S_pE2)

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
  /* component EMON_SAMPLE=E_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (0)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _a2->_rotation_absolute, _EMON_SAMPLE->_rotation_absolute);
    rot_transpose(_DivMon_Sample->_rotation_absolute, tr1);
    rot_mul(_EMON_SAMPLE->_rotation_absolute, tr1, _EMON_SAMPLE->_rotation_relative);
    _EMON_SAMPLE->_rotation_is_identity =  rot_test_identity(_EMON_SAMPLE->_rotation_relative);
    tc1 = coords_set(
      0, 0, 1.452);
    rot_transpose(_a2->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _EMON_SAMPLE->_position_absolute = coords_add(_a2->_position_absolute, tc2);
    tc1 = coords_sub(_DivMon_Sample->_position_absolute, _EMON_SAMPLE->_position_absolute);
    _EMON_SAMPLE->_position_relative = rot_apply(_EMON_SAMPLE->_rotation_absolute, tc1);
  } /* EMON_SAMPLE=E_monitor() AT ROTATED */
  DEBUG_COMPONENT("EMON_SAMPLE", _EMON_SAMPLE->_position_absolute, _EMON_SAMPLE->_rotation_absolute);
  instrument->_position_absolute[39] = _EMON_SAMPLE->_position_absolute;
  instrument->_position_relative[39] = _EMON_SAMPLE->_position_relative;
  instrument->counter_N[39]  = instrument->counter_P[39] = instrument->counter_P2[39] = 0;
  instrument->counter_AbsorbProp[39]= 0;
  return(0);
} /* _EMON_SAMPLE_setpos */

/* component a3=Arm() SETTING, POSITION/ROTATION */
int _a3_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_a3_setpos] component a3=Arm() SETTING [Arm:0]");
  stracpy(_a3->_name, "a3", 16384);
  stracpy(_a3->_type, "Arm", 16384);
  _a3->_index=40;
  /* component a3=Arm() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _a2->_rotation_absolute, _a3->_rotation_absolute);
    rot_transpose(_EMON_SAMPLE->_rotation_absolute, tr1);
    rot_mul(_a3->_rotation_absolute, tr1, _a3->_rotation_relative);
    _a3->_rotation_is_identity =  rot_test_identity(_a3->_rotation_relative);
    tc1 = coords_set(
      0, 0, 1.5);
    rot_transpose(_a2->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _a3->_position_absolute = coords_add(_a2->_position_absolute, tc2);
    tc1 = coords_sub(_EMON_SAMPLE->_position_absolute, _a3->_position_absolute);
    _a3->_position_relative = rot_apply(_a3->_rotation_absolute, tc1);
  } /* a3=Arm() AT ROTATED */
  DEBUG_COMPONENT("a3", _a3->_position_absolute, _a3->_rotation_absolute);
  instrument->_position_absolute[40] = _a3->_position_absolute;
  instrument->_position_relative[40] = _a3->_position_relative;
  instrument->counter_N[40]  = instrument->counter_P[40] = instrument->counter_P2[40] = 0;
  instrument->counter_AbsorbProp[40]= 0;
  return(0);
} /* _a3_setpos */

/* component Sample=V_sample() SETTING, POSITION/ROTATION */
int _Sample_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Sample_setpos] component Sample=V_sample() SETTING [/usr/share/mcstas/3.0-dev/obsolete/V_sample.comp:121]");
  stracpy(_Sample->_name, "Sample", 16384);
  stracpy(_Sample->_type, "V_sample", 16384);
  _Sample->_index=41;
  _Sample->_parameters.radius = 0.008;
  #define radius (_Sample->_parameters.radius)
  _Sample->_parameters.thickness = 0;
  #define thickness (_Sample->_parameters.thickness)
  _Sample->_parameters.zdepth = 0;
  #define zdepth (_Sample->_parameters.zdepth)
  _Sample->_parameters.Vc = 13.827;
  #define Vc (_Sample->_parameters.Vc)
  _Sample->_parameters.sigma_abs = 5.08;
  #define sigma_abs (_Sample->_parameters.sigma_abs)
  _Sample->_parameters.sigma_inc = 5.08;
  #define sigma_inc (_Sample->_parameters.sigma_inc)
  _Sample->_parameters.radius_i = 0;
  #define radius_i (_Sample->_parameters.radius_i)
  _Sample->_parameters.radius_o = 0;
  #define radius_o (_Sample->_parameters.radius_o)
  _Sample->_parameters.h = 0;
  #define h (_Sample->_parameters.h)
  _Sample->_parameters.focus_r = 0;
  #define focus_r (_Sample->_parameters.focus_r)
  _Sample->_parameters.pack = 1;
  #define pack (_Sample->_parameters.pack)
  _Sample->_parameters.frac = 1;
  #define frac (_Sample->_parameters.frac)
  _Sample->_parameters.f_QE = 0;
  #define f_QE (_Sample->_parameters.f_QE)
  _Sample->_parameters.gamma = 0;
  #define gamma (_Sample->_parameters.gamma)
  _Sample->_parameters.target_x = 0;
  #define target_x (_Sample->_parameters.target_x)
  _Sample->_parameters.target_y = 0;
  #define target_y (_Sample->_parameters.target_y)
  _Sample->_parameters.target_z = 0;
  #define target_z (_Sample->_parameters.target_z)
  _Sample->_parameters.focus_xw = 0.336;
  #define focus_xw (_Sample->_parameters.focus_xw)
  _Sample->_parameters.focus_yh = 0.4;
  #define focus_yh (_Sample->_parameters.focus_yh)
  _Sample->_parameters.focus_aw = 0;
  #define focus_aw (_Sample->_parameters.focus_aw)
  _Sample->_parameters.focus_ah = 0;
  #define focus_ah (_Sample->_parameters.focus_ah)
  _Sample->_parameters.xwidth = 0;
  #define xwidth (_Sample->_parameters.xwidth)
  _Sample->_parameters.yheight = 0.055;
  #define yheight (_Sample->_parameters.yheight)
  _Sample->_parameters.zthick = 0;
  #define zthick (_Sample->_parameters.zthick)
  _Sample->_parameters.rad_sphere = 0;
  #define rad_sphere (_Sample->_parameters.rad_sphere)
  _Sample->_parameters.sig_a = 0;
  #define sig_a (_Sample->_parameters.sig_a)
  _Sample->_parameters.sig_i = 0;
  #define sig_i (_Sample->_parameters.sig_i)
  _Sample->_parameters.V0 = 0;
  #define V0 (_Sample->_parameters.V0)
  _Sample->_parameters.target_index = + 3;
  #define target_index (_Sample->_parameters.target_index)
  _Sample->_parameters.multiples = 1;
  #define multiples (_Sample->_parameters.multiples)

  #define VarsV (_Sample->_parameters.VarsV)

  #undef radius
  #undef thickness
  #undef zdepth
  #undef Vc
  #undef sigma_abs
  #undef sigma_inc
  #undef radius_i
  #undef radius_o
  #undef h
  #undef focus_r
  #undef pack
  #undef frac
  #undef f_QE
  #undef gamma
  #undef target_x
  #undef target_y
  #undef target_z
  #undef focus_xw
  #undef focus_yh
  #undef focus_aw
  #undef focus_ah
  #undef xwidth
  #undef yheight
  #undef zthick
  #undef rad_sphere
  #undef sig_a
  #undef sig_i
  #undef V0
  #undef target_index
  #undef multiples
  #undef VarsV
  /* component Sample=V_sample() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _a3->_rotation_absolute, _Sample->_rotation_absolute);
    rot_transpose(_a3->_rotation_absolute, tr1);
    rot_mul(_Sample->_rotation_absolute, tr1, _Sample->_rotation_relative);
    _Sample->_rotation_is_identity =  rot_test_identity(_Sample->_rotation_relative);
    tc1 = coords_set(
      0, 0, 0);
    rot_transpose(_a3->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Sample->_position_absolute = coords_add(_a3->_position_absolute, tc2);
    tc1 = coords_sub(_a3->_position_absolute, _Sample->_position_absolute);
    _Sample->_position_relative = rot_apply(_Sample->_rotation_absolute, tc1);
  } /* Sample=V_sample() AT ROTATED */
  DEBUG_COMPONENT("Sample", _Sample->_position_absolute, _Sample->_rotation_absolute);
  instrument->_position_absolute[41] = _Sample->_position_absolute;
  instrument->_position_relative[41] = _Sample->_position_relative;
  instrument->counter_N[41]  = instrument->counter_P[41] = instrument->counter_P2[41] = 0;
  instrument->counter_AbsorbProp[41]= 0;
  return(0);
} /* _Sample_setpos */

/* component TOF_Det=Monitor_nD() SETTING, POSITION/ROTATION */
int _TOF_Det_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_TOF_Det_setpos] component TOF_Det=Monitor_nD() SETTING [/usr/share/mcstas/3.0-dev/monitors/Monitor_nD.comp:228]");
  stracpy(_TOF_Det->_name, "TOF_Det", 16384);
  stracpy(_TOF_Det->_type, "Monitor_nD", 16384);
  _TOF_Det->_index=42;
  _TOF_Det->_parameters.user1 = '\0';
  #define user1 (_TOF_Det->_parameters.user1)
  _TOF_Det->_parameters.user2 = '\0';
  #define user2 (_TOF_Det->_parameters.user2)
  _TOF_Det->_parameters.user3 = '\0';
  #define user3 (_TOF_Det->_parameters.user3)
  _TOF_Det->_parameters.xwidth = 1.5;
  #define xwidth (_TOF_Det->_parameters.xwidth)
  _TOF_Det->_parameters.yheight = 0.2;
  #define yheight (_TOF_Det->_parameters.yheight)
  _TOF_Det->_parameters.zdepth = 0;
  #define zdepth (_TOF_Det->_parameters.zdepth)
  _TOF_Det->_parameters.xmin = 0;
  #define xmin (_TOF_Det->_parameters.xmin)
  _TOF_Det->_parameters.xmax = 0;
  #define xmax (_TOF_Det->_parameters.xmax)
  _TOF_Det->_parameters.ymin = 0;
  #define ymin (_TOF_Det->_parameters.ymin)
  _TOF_Det->_parameters.ymax = 0;
  #define ymax (_TOF_Det->_parameters.ymax)
  _TOF_Det->_parameters.zmin = 0;
  #define zmin (_TOF_Det->_parameters.zmin)
  _TOF_Det->_parameters.zmax = 0;
  #define zmax (_TOF_Det->_parameters.zmax)
  _TOF_Det->_parameters.bins = 100;
  #define bins (_TOF_Det->_parameters.bins)
  _TOF_Det->_parameters.min = -1e40;
  #define min (_TOF_Det->_parameters.min)
  _TOF_Det->_parameters.max = 1e40;
  #define max (_TOF_Det->_parameters.max)
  _TOF_Det->_parameters.restore_neutron = 0;
  #define restore_neutron (_TOF_Det->_parameters.restore_neutron)
  _TOF_Det->_parameters.radius = 0;
  #define radius (_TOF_Det->_parameters.radius)
  if("auto time, abs angle limits=[10 180], banana, parallel" && strlen("auto time, abs angle limits=[10 180], banana, parallel"))
    stracpy(_TOF_Det->_parameters.options, "auto time, abs angle limits=[10 180], banana, parallel" ? "auto time, abs angle limits=[10 180], banana, parallel" : "", 16384);
  else 
  _TOF_Det->_parameters.options[0]='\0';
  #define options (_TOF_Det->_parameters.options)
  if("NULL" && strlen("NULL"))
    stracpy(_TOF_Det->_parameters.filename, "NULL" ? "NULL" : "", 16384);
  else 
  _TOF_Det->_parameters.filename[0]='\0';
  #define filename (_TOF_Det->_parameters.filename)
  if("NULL" && strlen("NULL"))
    stracpy(_TOF_Det->_parameters.geometry, "NULL" ? "NULL" : "", 16384);
  else 
  _TOF_Det->_parameters.geometry[0]='\0';
  #define geometry (_TOF_Det->_parameters.geometry)
  if("NULL" && strlen("NULL"))
    stracpy(_TOF_Det->_parameters.username1, "NULL" ? "NULL" : "", 16384);
  else 
  _TOF_Det->_parameters.username1[0]='\0';
  #define username1 (_TOF_Det->_parameters.username1)
  if("NULL" && strlen("NULL"))
    stracpy(_TOF_Det->_parameters.username2, "NULL" ? "NULL" : "", 16384);
  else 
  _TOF_Det->_parameters.username2[0]='\0';
  #define username2 (_TOF_Det->_parameters.username2)
  if("NULL" && strlen("NULL"))
    stracpy(_TOF_Det->_parameters.username3, "NULL" ? "NULL" : "", 16384);
  else 
  _TOF_Det->_parameters.username3[0]='\0';
  #define username3 (_TOF_Det->_parameters.username3)

  #define DEFS (_TOF_Det->_parameters.DEFS)
  #define Vars (_TOF_Det->_parameters.Vars)
  #define detector (_TOF_Det->_parameters.detector)
  #define offdata (_TOF_Det->_parameters.offdata)

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
  /* component TOF_Det=Monitor_nD() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Sample->_rotation_absolute, _TOF_Det->_rotation_absolute);
    rot_transpose(_Sample->_rotation_absolute, tr1);
    rot_mul(_TOF_Det->_rotation_absolute, tr1, _TOF_Det->_rotation_relative);
    _TOF_Det->_rotation_is_identity =  rot_test_identity(_TOF_Det->_rotation_relative);
    tc1 = coords_set(
      0, 0, 0);
    rot_transpose(_Sample->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _TOF_Det->_position_absolute = coords_add(_Sample->_position_absolute, tc2);
    tc1 = coords_sub(_Sample->_position_absolute, _TOF_Det->_position_absolute);
    _TOF_Det->_position_relative = rot_apply(_TOF_Det->_rotation_absolute, tc1);
  } /* TOF_Det=Monitor_nD() AT ROTATED */
  DEBUG_COMPONENT("TOF_Det", _TOF_Det->_position_absolute, _TOF_Det->_rotation_absolute);
  instrument->_position_absolute[42] = _TOF_Det->_position_absolute;
  instrument->_position_relative[42] = _TOF_Det->_position_relative;
  instrument->counter_N[42]  = instrument->counter_P[42] = instrument->counter_P2[42] = 0;
  instrument->counter_AbsorbProp[42]= 0;
  return(0);
} /* _TOF_Det_setpos */

/* component FoDet=Monitor_nD() SETTING, POSITION/ROTATION */
int _FoDet_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_FoDet_setpos] component FoDet=Monitor_nD() SETTING [/usr/share/mcstas/3.0-dev/monitors/Monitor_nD.comp:228]");
  stracpy(_FoDet->_name, "FoDet", 16384);
  stracpy(_FoDet->_type, "Monitor_nD", 16384);
  _FoDet->_index=43;
  _FoDet->_parameters.user1 = '\0';
  #define user1 (_FoDet->_parameters.user1)
  _FoDet->_parameters.user2 = '\0';
  #define user2 (_FoDet->_parameters.user2)
  _FoDet->_parameters.user3 = '\0';
  #define user3 (_FoDet->_parameters.user3)
  _FoDet->_parameters.xwidth = 0.336;
  #define xwidth (_FoDet->_parameters.xwidth)
  _FoDet->_parameters.yheight = 0.4;
  #define yheight (_FoDet->_parameters.yheight)
  _FoDet->_parameters.zdepth = 0;
  #define zdepth (_FoDet->_parameters.zdepth)
  _FoDet->_parameters.xmin = 0;
  #define xmin (_FoDet->_parameters.xmin)
  _FoDet->_parameters.xmax = 0;
  #define xmax (_FoDet->_parameters.xmax)
  _FoDet->_parameters.ymin = 0;
  #define ymin (_FoDet->_parameters.ymin)
  _FoDet->_parameters.ymax = 0;
  #define ymax (_FoDet->_parameters.ymax)
  _FoDet->_parameters.zmin = 0;
  #define zmin (_FoDet->_parameters.zmin)
  _FoDet->_parameters.zmax = 0;
  #define zmax (_FoDet->_parameters.zmax)
  _FoDet->_parameters.bins = 0;
  #define bins (_FoDet->_parameters.bins)
  _FoDet->_parameters.min = -1e40;
  #define min (_FoDet->_parameters.min)
  _FoDet->_parameters.max = 1e40;
  #define max (_FoDet->_parameters.max)
  _FoDet->_parameters.restore_neutron = 0;
  #define restore_neutron (_FoDet->_parameters.restore_neutron)
  _FoDet->_parameters.radius = 0;
  #define radius (_FoDet->_parameters.radius)
  if("t auto file=TofDet.dat" && strlen("t auto file=TofDet.dat"))
    stracpy(_FoDet->_parameters.options, "t auto file=TofDet.dat" ? "t auto file=TofDet.dat" : "", 16384);
  else 
  _FoDet->_parameters.options[0]='\0';
  #define options (_FoDet->_parameters.options)
  if("NULL" && strlen("NULL"))
    stracpy(_FoDet->_parameters.filename, "NULL" ? "NULL" : "", 16384);
  else 
  _FoDet->_parameters.filename[0]='\0';
  #define filename (_FoDet->_parameters.filename)
  if("NULL" && strlen("NULL"))
    stracpy(_FoDet->_parameters.geometry, "NULL" ? "NULL" : "", 16384);
  else 
  _FoDet->_parameters.geometry[0]='\0';
  #define geometry (_FoDet->_parameters.geometry)
  if("NULL" && strlen("NULL"))
    stracpy(_FoDet->_parameters.username1, "NULL" ? "NULL" : "", 16384);
  else 
  _FoDet->_parameters.username1[0]='\0';
  #define username1 (_FoDet->_parameters.username1)
  if("NULL" && strlen("NULL"))
    stracpy(_FoDet->_parameters.username2, "NULL" ? "NULL" : "", 16384);
  else 
  _FoDet->_parameters.username2[0]='\0';
  #define username2 (_FoDet->_parameters.username2)
  if("NULL" && strlen("NULL"))
    stracpy(_FoDet->_parameters.username3, "NULL" ? "NULL" : "", 16384);
  else 
  _FoDet->_parameters.username3[0]='\0';
  #define username3 (_FoDet->_parameters.username3)

  #define DEFS (_FoDet->_parameters.DEFS)
  #define Vars (_FoDet->_parameters.Vars)
  #define detector (_FoDet->_parameters.detector)
  #define offdata (_FoDet->_parameters.offdata)

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
  /* component FoDet=Monitor_nD() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (instrument->_parameters._DET)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _a3->_rotation_absolute, _FoDet->_rotation_absolute);
    rot_transpose(_TOF_Det->_rotation_absolute, tr1);
    rot_mul(_FoDet->_rotation_absolute, tr1, _FoDet->_rotation_relative);
    _FoDet->_rotation_is_identity =  rot_test_identity(_FoDet->_rotation_relative);
    tc1 = coords_set(
      0, 0, 2.5);
    rot_transpose(_a3->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _FoDet->_position_absolute = coords_add(_a3->_position_absolute, tc2);
    tc1 = coords_sub(_TOF_Det->_position_absolute, _FoDet->_position_absolute);
    _FoDet->_position_relative = rot_apply(_FoDet->_rotation_absolute, tc1);
  } /* FoDet=Monitor_nD() AT ROTATED */
  DEBUG_COMPONENT("FoDet", _FoDet->_position_absolute, _FoDet->_rotation_absolute);
  instrument->_position_absolute[43] = _FoDet->_position_absolute;
  instrument->_position_relative[43] = _FoDet->_position_relative;
  instrument->counter_N[43]  = instrument->counter_P[43] = instrument->counter_P2[43] = 0;
  instrument->counter_AbsorbProp[43]= 0;
  return(0);
} /* _FoDet_setpos */

/* component EMON_DET=E_monitor() SETTING, POSITION/ROTATION */
int _EMON_DET_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_EMON_DET_setpos] component EMON_DET=E_monitor() SETTING [/usr/share/mcstas/3.0-dev/monitors/E_monitor.comp:66]");
  stracpy(_EMON_DET->_name, "EMON_DET", 16384);
  stracpy(_EMON_DET->_type, "E_monitor", 16384);
  _EMON_DET->_index=44;
  _EMON_DET->_parameters.nE = 80;
  #define nE (_EMON_DET->_parameters.nE)
  if("emon_det.dat" && strlen("emon_det.dat"))
    stracpy(_EMON_DET->_parameters.filename, "emon_det.dat" ? "emon_det.dat" : "", 16384);
  else 
  _EMON_DET->_parameters.filename[0]='\0';
  #define filename (_EMON_DET->_parameters.filename)
  _EMON_DET->_parameters.xmin = -0.05;
  #define xmin (_EMON_DET->_parameters.xmin)
  _EMON_DET->_parameters.xmax = 0.05;
  #define xmax (_EMON_DET->_parameters.xmax)
  _EMON_DET->_parameters.ymin = -0.05;
  #define ymin (_EMON_DET->_parameters.ymin)
  _EMON_DET->_parameters.ymax = 0.05;
  #define ymax (_EMON_DET->_parameters.ymax)
  _EMON_DET->_parameters.xwidth = 0.4;
  #define xwidth (_EMON_DET->_parameters.xwidth)
  _EMON_DET->_parameters.yheight = 0.4;
  #define yheight (_EMON_DET->_parameters.yheight)
  _EMON_DET->_parameters.Emin = EMIN;
  #define Emin (_EMON_DET->_parameters.Emin)
  _EMON_DET->_parameters.Emax = EMAX;
  #define Emax (_EMON_DET->_parameters.Emax)
  _EMON_DET->_parameters.restore_neutron = 0;
  #define restore_neutron (_EMON_DET->_parameters.restore_neutron)

  #define E_N (_EMON_DET->_parameters.E_N)
  #define E_p (_EMON_DET->_parameters.E_p)
  #define E_p2 (_EMON_DET->_parameters.E_p2)
  #define S_p (_EMON_DET->_parameters.S_p)
  #define S_pE (_EMON_DET->_parameters.S_pE)
  #define S_pE2 (_EMON_DET->_parameters.S_pE2)

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
  /* component EMON_DET=E_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    Rotation tr1;
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (instrument->_parameters._DET)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _a3->_rotation_absolute, _EMON_DET->_rotation_absolute);
    rot_transpose(_FoDet->_rotation_absolute, tr1);
    rot_mul(_EMON_DET->_rotation_absolute, tr1, _EMON_DET->_rotation_relative);
    _EMON_DET->_rotation_is_identity =  rot_test_identity(_EMON_DET->_rotation_relative);
    tc1 = coords_set(
      0, 0, 2.501);
    rot_transpose(_a3->_rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _EMON_DET->_position_absolute = coords_add(_a3->_position_absolute, tc2);
    tc1 = coords_sub(_FoDet->_position_absolute, _EMON_DET->_position_absolute);
    _EMON_DET->_position_relative = rot_apply(_EMON_DET->_rotation_absolute, tc1);
  } /* EMON_DET=E_monitor() AT ROTATED */
  DEBUG_COMPONENT("EMON_DET", _EMON_DET->_position_absolute, _EMON_DET->_rotation_absolute);
  instrument->_position_absolute[44] = _EMON_DET->_position_absolute;
  instrument->_position_relative[44] = _EMON_DET->_position_relative;
  instrument->counter_N[44]  = instrument->counter_P[44] = instrument->counter_P2[44] = 0;
  instrument->counter_AbsorbProp[44]= 0;
  return(0);
} /* _EMON_DET_setpos */

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

_class_Bender *class_Bender_init(_class_Bender *_comp
) {
  #define w (_comp->_parameters.w)
  #define h (_comp->_parameters.h)
  #define r (_comp->_parameters.r)
  #define Win (_comp->_parameters.Win)
  #define k (_comp->_parameters.k)
  #define d (_comp->_parameters.d)
  #define l (_comp->_parameters.l)
  #define R0a (_comp->_parameters.R0a)
  #define Qca (_comp->_parameters.Qca)
  #define alphaa (_comp->_parameters.alphaa)
  #define ma (_comp->_parameters.ma)
  #define Wa (_comp->_parameters.Wa)
  #define R0i (_comp->_parameters.R0i)
  #define Qci (_comp->_parameters.Qci)
  #define alphai (_comp->_parameters.alphai)
  #define mi (_comp->_parameters.mi)
  #define Wi (_comp->_parameters.Wi)
  #define R0s (_comp->_parameters.R0s)
  #define Qcs (_comp->_parameters.Qcs)
  #define alphas (_comp->_parameters.alphas)
  #define ms (_comp->_parameters.ms)
  #define Ws (_comp->_parameters.Ws)
  #define bk (_comp->_parameters.bk)
  #define mWin (_comp->_parameters.mWin)
if (r <0)
      { fprintf(stderr,"Bender: error: %s: to bend in the other direction\n", NAME_CURRENT_COMP);
        fprintf(stderr,"        rotate comp on z-axis by 180 deg.\n"); exit(-1); }

      if (k*d > w)
      { fprintf(stderr,"Bender: error: %s has (k*d > w).\n", NAME_CURRENT_COMP);
        exit(-1); }
      if (w*h*r*Win*k == 0)
      { fprintf(stderr,"Bender: error: %s has one of w,h,r,Win,k null.\n", NAME_CURRENT_COMP);
        exit(-1); }
      /* width of one channel + thickness d of partition */
      mWin = Win;
      if (l!= 0 && r != 0) mWin = (double)l/(double)r;
      bk=(w+d)/k;
      if (mcgravitation) fprintf(stderr,"WARNING: Bender: %s: "
        "This component produces wrong results with gravitation !\n",
        NAME_CURRENT_COMP);
  #undef w
  #undef h
  #undef r
  #undef Win
  #undef k
  #undef d
  #undef l
  #undef R0a
  #undef Qca
  #undef alphaa
  #undef ma
  #undef Wa
  #undef R0i
  #undef Qci
  #undef alphai
  #undef mi
  #undef Wi
  #undef R0s
  #undef Qcs
  #undef alphas
  #undef ms
  #undef Ws
  #undef bk
  #undef mWin
  return(_comp);
} /* class_Bender_init */

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

_class_Guide_channeled *class_Guide_channeled_init(_class_Guide_channeled *_comp
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
  #define nslit (_comp->_parameters.nslit)
  #define d (_comp->_parameters.d)
  #define Qcx (_comp->_parameters.Qcx)
  #define Qcy (_comp->_parameters.Qcy)
  #define alphax (_comp->_parameters.alphax)
  #define alphay (_comp->_parameters.alphay)
  #define W (_comp->_parameters.W)
  #define mx (_comp->_parameters.mx)
  #define my (_comp->_parameters.my)
  #define nu (_comp->_parameters.nu)
  #define phase (_comp->_parameters.phase)
  #define w1c (_comp->_parameters.w1c)
  #define w2c (_comp->_parameters.w2c)
  #define ww (_comp->_parameters.ww)
  #define hh (_comp->_parameters.hh)
  #define whalf (_comp->_parameters.whalf)
  #define hhalf (_comp->_parameters.hhalf)
  #define lwhalf (_comp->_parameters.lwhalf)
  #define lhhalf (_comp->_parameters.lhhalf)
if (!w2) w2=w1;
  if (!h2) h2=h1;
  if (nslit <= 0 || W <=0)
  { fprintf(stderr,"Guide_channeled: %s: nslit and W must be positive\n", NAME_CURRENT_COMP);
    exit(-1); }
  w1c = (w1 + d)/(double)nslit;
  w2c = (w2 + d)/(double)nslit;
  ww = .5*(w2c - w1c);
  hh = .5*(h2 - h1);
  whalf = .5*(w1c - d);
  hhalf = .5*h1;
  lwhalf = l*whalf;
  lhhalf = l*hhalf;

  if (m)     { mx=my=m; }
  if (Qc)    { Qcx=Qcy=Qc; }
  if (alpha) { alphax=alphay=alpha; }

  if ((nslit > 1) && (w1 != w2))
  {
    fprintf(stderr,"WARNING: Guide_channeled: %s:"
    "This component does not work with multichannel focusing guide\n"
    "Use Guide_gravity for that.\n", NAME_CURRENT_COMP);
    exit(-1);
  }

  if (d*nslit > w1) exit(fprintf(stderr, "Guide_channeled: %s: absorbing walls fill input window. No space left for transmission (d*nslit > w1).\n", NAME_CURRENT_COMP));

  if (mcgravitation) fprintf(stderr,"WARNING: Guide_channeled: %s: "
    "This component produces wrong results with gravitation !\n"
    "Use Guide_gravity.\n",
    NAME_CURRENT_COMP);
  if (nu != 0 || phase != 0) {
      if (w1 != w2 || h1 != h2)
      exit(fprintf(stderr,"Guide_channeled: %s: rotating slit pack must be straight (w1=w2 and h1=h2).\n", NAME_CURRENT_COMP));
      printf("Guide_channeled: %s: Fermi Chopper mode: frequency=%g [Hz] phase=%g [deg]\n",
        NAME_CURRENT_COMP, nu, phase);
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
  #undef nslit
  #undef d
  #undef Qcx
  #undef Qcy
  #undef alphax
  #undef alphay
  #undef W
  #undef mx
  #undef my
  #undef nu
  #undef phase
  #undef w1c
  #undef w2c
  #undef ww
  #undef hh
  #undef whalf
  #undef hhalf
  #undef lwhalf
  #undef lhhalf
  return(_comp);
} /* class_Guide_channeled_init */

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

_class_FermiChopper *class_FermiChopper_init(_class_FermiChopper *_comp
) {
  #define phase (_comp->_parameters.phase)
  #define radius (_comp->_parameters.radius)
  #define nu (_comp->_parameters.nu)
  #define w (_comp->_parameters.w)
  #define nslit (_comp->_parameters.nslit)
  #define R0 (_comp->_parameters.R0)
  #define Qc (_comp->_parameters.Qc)
  #define alpha (_comp->_parameters.alpha)
  #define m (_comp->_parameters.m)
  #define W (_comp->_parameters.W)
  #define length (_comp->_parameters.length)
  #define eff (_comp->_parameters.eff)
  #define zero_time (_comp->_parameters.zero_time)
  #define xwidth (_comp->_parameters.xwidth)
  #define verbose (_comp->_parameters.verbose)
  #define yheight (_comp->_parameters.yheight)
  #define curvature (_comp->_parameters.curvature)
  #define delay (_comp->_parameters.delay)
  #define FCVars (_comp->_parameters.FCVars)

/************************ CALCULATION CONSTANTS *****************************/
  strcpy(FCVars.compcurname, NAME_CURRENT_COMP);

  FCVars.omega    = 2*PI*nu;
  if (!phase && delay) {
     FCVars.ph0= fmod(-delay*nu*360,360)*DEG2RAD;
  } else FCVars.ph0      = phase*DEG2RAD;
  FCVars.sum_t=FCVars.sum_v=FCVars.sum_N=FCVars.sum_N_pass=0;

  /* check of input parameters */
  if (nslit < 1) nslit=1;
  if (yheight <= 0) exit(printf("FermiChopper: %s: FATAL: unrealistic cylinder yheight =%g [m]\n", NAME_CURRENT_COMP, yheight));

  if (m <= 0) { m=0; R0=0; }
  if (radius <= 0) {
    printf("FermiChopper: %s: FATAL: Unrealistic cylinder radius radius=%g [m]\n", NAME_CURRENT_COMP, radius);
    exit(-1);
  }
  if (xwidth > 0 && xwidth < radius*2 && nslit > 0) {
    w = xwidth/nslit;
  }
  if (w <= 0) {
    printf("FermiChopper: %s: FATAL: Slits in the package have unrealistic width w=%g [m]\n", NAME_CURRENT_COMP, w);
    exit(-1);
  }
  if (nslit*w > radius*2) {
    nslit = floor(radius/w);
    printf("FermiChopper: %s: Too many slits to fit in the cylinder\n"
           "    Adjusting nslit=%f\n", NAME_CURRENT_COMP, nslit);
  }
  if (length > radius*2) {
    length = 2*sqrt(radius*radius - nslit*w*nslit*w/4);
    printf("FermiChopper: %s: Slit package is longer than the whole\n"
           "    chopper cylinder. Adjusting length=%g [m]\n", NAME_CURRENT_COMP, length);
  }

  if (eff <= 0 || eff > 1) {
    eff = 0.95;
    printf("FermiChopper: %s: Efficiency is unrealistic\n"
           "    Adjusting eff=%f\n", NAME_CURRENT_COMP, eff);
  }
  if (Qc <= 0) { Qc = 0.02176; m = 0; R0 = 0; }
  if (W <= 0) W=1e-6;

  if (curvature) {
    FCVars.C_slit = curvature;
    if (1 < fabs(radius*curvature))
      exit(printf("FermiChopper: %s: Slit curvature is unrealistic\n",
           NAME_CURRENT_COMP));
  }
  FCVars.L_slit = length;
  if (verbose && nu)
    printf("FermiChopper: %s: Frequency nu=%g [Hz] %g [rpm], time frame=%g [s] phase=%g [deg]\n"
      , NAME_CURRENT_COMP, nu, nu*60, 2/nu, FCVars.ph0*RAD2DEG);

  FCVars.absorb_alreadyinside    = 0;
  FCVars.absorb_topbottom        = 0;
  FCVars.absorb_cylentrance      = 0;
  FCVars.absorb_sideentrance     = 0;
  FCVars.absorb_notreachentrance = 0;
  FCVars.absorb_packentrance     = 0;
  FCVars.absorb_slitcoating      = 0;
  FCVars.warn_notreachslitwall   = 0;
  FCVars.absorb_exitslitpack     = 0;
  FCVars.absorb_maxiterations    = 0;
  FCVars.absorb_wrongdirection   = 0;
  FCVars.absorb_nocontrol        = 0;
  FCVars.absorb_cylexit          = 0;
  FCVars.warn_notreachslitoutput = 0;

  /* fix for the wrong coordinate frame orientation to come back to McStas XYZ system */
  FCVars.omega *= -1;
  FCVars.ph0   *= -1;
  FCVars.t0     = -FCVars.ph0/FCVars.omega;
  #undef phase
  #undef radius
  #undef nu
  #undef w
  #undef nslit
  #undef R0
  #undef Qc
  #undef alpha
  #undef m
  #undef W
  #undef length
  #undef eff
  #undef zero_time
  #undef xwidth
  #undef verbose
  #undef yheight
  #undef curvature
  #undef delay
  #undef FCVars
  return(_comp);
} /* class_FermiChopper_init */

_class_V_sample *class_V_sample_init(_class_V_sample *_comp
) {
  #define radius (_comp->_parameters.radius)
  #define thickness (_comp->_parameters.thickness)
  #define zdepth (_comp->_parameters.zdepth)
  #define Vc (_comp->_parameters.Vc)
  #define sigma_abs (_comp->_parameters.sigma_abs)
  #define sigma_inc (_comp->_parameters.sigma_inc)
  #define radius_i (_comp->_parameters.radius_i)
  #define radius_o (_comp->_parameters.radius_o)
  #define h (_comp->_parameters.h)
  #define focus_r (_comp->_parameters.focus_r)
  #define pack (_comp->_parameters.pack)
  #define frac (_comp->_parameters.frac)
  #define f_QE (_comp->_parameters.f_QE)
  #define gamma (_comp->_parameters.gamma)
  #define target_x (_comp->_parameters.target_x)
  #define target_y (_comp->_parameters.target_y)
  #define target_z (_comp->_parameters.target_z)
  #define focus_xw (_comp->_parameters.focus_xw)
  #define focus_yh (_comp->_parameters.focus_yh)
  #define focus_aw (_comp->_parameters.focus_aw)
  #define focus_ah (_comp->_parameters.focus_ah)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define zthick (_comp->_parameters.zthick)
  #define rad_sphere (_comp->_parameters.rad_sphere)
  #define sig_a (_comp->_parameters.sig_a)
  #define sig_i (_comp->_parameters.sig_i)
  #define V0 (_comp->_parameters.V0)
  #define target_index (_comp->_parameters.target_index)
  #define multiples (_comp->_parameters.multiples)
  #define VarsV (_comp->_parameters.VarsV)
  /* Backward compatibility */
  if (radius) radius_o = radius;
  if (thickness) radius_i = radius_o - thickness;
  if (zdepth) zthick = zdepth;
  if (yheight) h = yheight;
  if (Vc) V0 = Vc;
  if (sigma_abs) sig_a = sigma_abs;
  if (sigma_inc) sig_i = sigma_inc;

  VarsV.shapetyp = -1;
  if (xwidth && yheight && zdepth)  VarsV.shapetyp=1; /* box */
  else if (radius > 0 && yheight)        VarsV.shapetyp=0; /* cylinder */
  else if (radius && !yheight)           VarsV.shapetyp=2; /* sphere */
  
  if (VarsV.shapetyp < 0)
    exit(fprintf(stderr,"V_sample: %s: sample has invalid dimensions. Please check parameter values.\n", NAME_CURRENT_COMP));

  VarsV.sigma_a=sig_a;
  VarsV.sigma_i=sig_i;
  VarsV.rho = (pack/V0);
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
    printf("V_sample: %s: The target is not defined. Using direct beam (Z-axis).\n",
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
  #undef radius
  #undef thickness
  #undef zdepth
  #undef Vc
  #undef sigma_abs
  #undef sigma_inc
  #undef radius_i
  #undef radius_o
  #undef h
  #undef focus_r
  #undef pack
  #undef frac
  #undef f_QE
  #undef gamma
  #undef target_x
  #undef target_y
  #undef target_z
  #undef focus_xw
  #undef focus_yh
  #undef focus_aw
  #undef focus_ah
  #undef xwidth
  #undef yheight
  #undef zthick
  #undef rad_sphere
  #undef sig_a
  #undef sig_i
  #undef V0
  #undef target_index
  #undef multiples
  #undef VarsV
  return(_comp);
} /* class_V_sample_init */



int init(void) { /* called by mccode_main for PSI_Focus:INITIALISE */
  DEBUG_INSTR();

  /* code_main/parseoptions/readparams sets instrument parameters value */
  stracpy(instrument->_name, "PSI_Focus", 256);

  /* Instrument 'PSI_Focus' INITIALISE */
  SIG_MESSAGE("[PSI_Focus] INITIALISE [PSI_Focus.instr:62]");
  #define lambda (instrument->_parameters._lambda)
  #define chopp_ratio (instrument->_parameters._chopp_ratio)
  #define DET (instrument->_parameters._DET)
{

 /* setting theta and 2Theta */
 PHM = asin(lambda/2/monodd)*180/PI;
 TTM = 2*PHM;

 /* setting speed of Fermi chopper, Disk chopper */
 FERMI_SPEED=1.0e6/(PI*lambda*c*(1./tan(PHM*PI/180))*(dsd*pow((1),1.5)+dfcs)*(1.0-(dms/dgm)));
 DISC_SPEED = FERMI_SPEED/chopp_ratio;

 /* calculate the phase between Fermi and Disk chopper */
 FO_PHA = (c*lambda*(ddcm+mfc))/1e6;

 /* setting monochromator curvature */
 RV_2 = 2*dms*(sin(DEG2RAD*PHM));
 RH_2 = 2*dms/(sin(DEG2RAD*PHM));

 /* setting E and L monitors */
 LMIN = lambda - dL;
 LMAX = lambda + dL;
 EMIN = 81.81/(LMAX*LMAX);
 EMAX = 81.81/(LMIN*LMIN);

}
  #undef lambda
  #undef chopp_ratio
  #undef DET
  _a1_setpos(); /* type Progress_bar */
  _csource_setpos(); /* type Source_gen */
  _guide1_setpos(); /* type Guide */
  _guide2_setpos(); /* type Bender */
  _bunker_setpos(); /* type Guide */
  _guide3_setpos(); /* type Guide */
  _lambdaGuideExit_setpos(); /* type L_monitor */
  _DivMonGuideExit_setpos(); /* type Divergence_monitor */
  _PSDGuideExit_setpos(); /* type PSD_monitor */
  _FOCUSguide_setpos(); /* type Guide_channeled */
  _FirstChopper_setpos(); /* type DiskChopper */
  _DISCTOF_setpos(); /* type Monitor_nD */
  _PSDmon1Chopper_setpos(); /* type PSD_monitor */
  _VacuumTube1entry_setpos(); /* type Slit */
  _VacuumTube1exit_setpos(); /* type Slit */
  _VacuumTube2entry_setpos(); /* type Slit */
  _VacuumTube2exit_setpos(); /* type Slit */
  _VacuumTube3entry_setpos(); /* type Slit */
  _VacuumTube3exit_setpos(); /* type Slit */
  _PSDmonMono_setpos(); /* type PSD_monitor */
  _MONOTOF_setpos(); /* type Monitor_nD */
  _DivMonMono_setpos(); /* type Divergence_monitor */
  _focus_mono_setpos(); /* type Arm */
  _mono_setpos(); /* type Monochromator_2foc */
  _a2_setpos(); /* type Arm */
  _FERMITOF_before_setpos(); /* type Monitor_nD */
  _lambdaFermi_setpos(); /* type L_monitor */
  _EMON_Fermi_setpos(); /* type E_monitor */
  _DivMonfermi1_setpos(); /* type Divergence_monitor */
  _PSD_Fermi1_setpos(); /* type PSD_monitor */
  _FoChopper_setpos(); /* type FermiChopper */
  _PSD_Fermi2_setpos(); /* type PSD_monitor */
  _DivMonfermi2_setpos(); /* type Divergence_monitor */
  _FERMITOF1_setpos(); /* type Monitor_nD */
  _SAMPLE_SLIT_setpos(); /* type Slit */
  _FERMITOF2_setpos(); /* type Monitor_nD */
  _PSD_SAMPLE_setpos(); /* type PSD_monitor */
  _DivMon_Sample_setpos(); /* type Divergence_monitor */
  _EMON_SAMPLE_setpos(); /* type E_monitor */
  _a3_setpos(); /* type Arm */
  _Sample_setpos(); /* type V_sample */
  _TOF_Det_setpos(); /* type Monitor_nD */
  _FoDet_setpos(); /* type Monitor_nD */
  _EMON_DET_setpos(); /* type E_monitor */

  /* call iteratively all components INITIALISE */
  class_Progress_bar_init(_a1);

  class_Source_gen_init(_csource);

  class_Guide_init(_guide1);

  class_Bender_init(_guide2);

  class_Guide_init(_bunker);

  class_Guide_init(_guide3);

  class_L_monitor_init(_lambdaGuideExit);

  class_Divergence_monitor_init(_DivMonGuideExit);

  class_PSD_monitor_init(_PSDGuideExit);

  class_Guide_channeled_init(_FOCUSguide);

  class_DiskChopper_init(_FirstChopper);

  class_Monitor_nD_init(_DISCTOF);

  class_PSD_monitor_init(_PSDmon1Chopper);

  class_Slit_init(_VacuumTube1entry);

  class_Slit_init(_VacuumTube1exit);

  class_Slit_init(_VacuumTube2entry);

  class_Slit_init(_VacuumTube2exit);

  class_Slit_init(_VacuumTube3entry);

  class_Slit_init(_VacuumTube3exit);

  class_PSD_monitor_init(_PSDmonMono);

  class_Monitor_nD_init(_MONOTOF);

  class_Divergence_monitor_init(_DivMonMono);


  class_Monochromator_2foc_init(_mono);


  class_Monitor_nD_init(_FERMITOF_before);

  class_L_monitor_init(_lambdaFermi);

  class_E_monitor_init(_EMON_Fermi);

  class_Divergence_monitor_init(_DivMonfermi1);

  class_PSD_monitor_init(_PSD_Fermi1);

  class_FermiChopper_init(_FoChopper);

  class_PSD_monitor_init(_PSD_Fermi2);

  class_Divergence_monitor_init(_DivMonfermi2);

  class_Monitor_nD_init(_FERMITOF1);

  class_Slit_init(_SAMPLE_SLIT);

  class_Monitor_nD_init(_FERMITOF2);

  class_PSD_monitor_init(_PSD_SAMPLE);

  class_Divergence_monitor_init(_DivMon_Sample);

  class_E_monitor_init(_EMON_SAMPLE);


  class_V_sample_init(_Sample);

  class_Monitor_nD_init(_TOF_Det);

  class_Monitor_nD_init(_FoDet);

  class_E_monitor_init(_EMON_DET);

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
_class_Bender *class_Bender_trace(_class_Bender *_comp
  , _class_particle *_particle) {
  ABSORBED=SCATTERED=RESTORE=0;

  #define w (_comp->_parameters.w)
  #define h (_comp->_parameters.h)
  #define r (_comp->_parameters.r)
  #define Win (_comp->_parameters.Win)
  #define k (_comp->_parameters.k)
  #define d (_comp->_parameters.d)
  #define l (_comp->_parameters.l)
  #define R0a (_comp->_parameters.R0a)
  #define Qca (_comp->_parameters.Qca)
  #define alphaa (_comp->_parameters.alphaa)
  #define ma (_comp->_parameters.ma)
  #define Wa (_comp->_parameters.Wa)
  #define R0i (_comp->_parameters.R0i)
  #define Qci (_comp->_parameters.Qci)
  #define alphai (_comp->_parameters.alphai)
  #define mi (_comp->_parameters.mi)
  #define Wi (_comp->_parameters.Wi)
  #define R0s (_comp->_parameters.R0s)
  #define Qcs (_comp->_parameters.Qcs)
  #define alphas (_comp->_parameters.alphas)
  #define ms (_comp->_parameters.ms)
  #define Ws (_comp->_parameters.Ws)
  #define bk (_comp->_parameters.bk)
  #define mWin (_comp->_parameters.mWin)
    int i,num,numa,numi;
    double dru,ab,dab,R,Q,Ta,vpl;
    double einmWin,ausmWin,zykmWin,aeumWin,innmWin,ref,innref,aeuref;
    double einzei,auszei,zykzei;

    /* does the neutron hit the bender at the entrance? */
    PROP_Z0;
    if ((fabs(x)<w/2) && (fabs(y)<h/2))
    {
      /*** reflections in the XZ-plane ***/

      /* distance between neutron and concave side of the channel at the entrance */
      dru=floor((w/2-x)/bk)*bk;
      ab=w/2.0-x-dru;

      /* radius of the channel */
      R=r-dru;

      /* does the neutron hit the partition at the entrance? */
      if (ab<bk-d)
      {
        double aeu[] = {R0a, Qca, alphaa, ma, Wa};
        /* velocity in the XZ-plane */
        vpl=sqrt(vx*vx+vz*vz);

        /* divergence of the neutron at the entrance */
        einmWin=atan(vx/vz);

        /* maximal distance between neutron and concave side of the channel */
        dab=R-cos(einmWin)*(R-ab);

        /* reflection angle at the concave side */
        aeumWin=acos((R-dab)/R);

        /* reflection coefficient at the concave side */
        Q=2.0*V2K*vpl*sin(aeumWin);
        StdReflecFunc(Q, aeu, &aeuref);

        /* does the neutron hit the convex side of the channel? */
        innmWin=0.0;
        innref=1.0;
        if (dab>bk-d)
        {
           double inn[] = {R0i, Qci, alphai, mi, Wi};
           /* reflection coefficient at the convex side */
           innmWin=acos((R-dab)/(R-bk+d));
           Q=2.0*V2K*vpl*sin(innmWin);
           StdReflecFunc(Q, inn, &innref);
        }

        /* divergence of the neutron at the exit */
        zykmWin=2.0*(aeumWin-innmWin);
        ausmWin=fmod(mWin+einmWin+aeumWin-innmWin
          *(1.0+SIGN(einmWin)),zykmWin)-zykmWin/2.0;
        ausmWin+=innmWin*SIGN(ausmWin);

        /* number of reflections at the concave side */
        numa=(mWin+einmWin+aeumWin-innmWin*(1.0+SIGN(einmWin)))/zykmWin;

        /* number of reflections at the convex side */
        numi=numa;
        if (ausmWin*einmWin<0)
        {
           if (ausmWin-einmWin>0)
              numi++;
           else
              numi--;
        }

        /* is the reflection coefficient too small? */
        if (((numa>0) && (aeuref<=0)) || ((numi>0) && (innref<=0)))
           ABSORB;

        /* calculation of the neutron probability weight p */
        for (i=1;i<=numa;i++)
            p*=aeuref;
        for (i=1;i<=numi;i++)
            p*=innref;

        /* time to cross the bender */
        Ta=(2*numa*(tan(aeumWin)-tan(innmWin))
          +tan(ausmWin)-tan(einmWin)
          -tan(innmWin)*(SIGN(ausmWin)-SIGN(einmWin)))
          *(R-dab)/vpl;
        t+=Ta;

        /* distance between neutron and concave side of channel at the exit */
        ab=R-(R-dab)/cos(ausmWin);

        /* calculation of the exit coordinates in the XZ-plane */
        x=w/2.0-ab-dru;
        z=r*mWin;
        vx=sin(ausmWin)*vpl;
        vz=cos(ausmWin)*vpl;

        /*** reflections at top and bottom side (Y axis) ***/

        if (vy!=0.0)
        {
          double s[] = {R0s, Qcs, alphas, ms, Ws};
          /* reflection coefficent at the top and bottom side */
          Q=2.0*V2K*fabs(vy);
          StdReflecFunc(Q, s, &ref);

          /* number of reflections at top and bottom */
          einzei=h/2.0/fabs(vy)+y/vy;
          zykzei=h/fabs(vy);
          num=(Ta+einzei)/zykzei;

          /* time between the last reflection and the exit */
          auszei=fmod(Ta+einzei,zykzei);

          /* is the reflection coefficient too small? */
          if ((num>0) && (ref<=0))
             ABSORB;

          /* calculation of the probability weight p */
          for (i=1;i<=num;i++) {
               p*=ref;
               vy*=-1.0; }

          /* calculation of the exit coordinate */
          y=auszei*vy-vy*h/fabs(vy)/2.0;
        } /* if (vy!=0.0) */
        SCATTER;
      } /* if (dab>bk-d)  */
      else
        ABSORB; /* hit separating walls */
    }
    else /* if ((fabs(x)<w/2) && (fabs(y)<h/2))   */
      ABSORB; /* miss entry window */

  #undef w
  #undef h
  #undef r
  #undef Win
  #undef k
  #undef d
  #undef l
  #undef R0a
  #undef Qca
  #undef alphaa
  #undef ma
  #undef Wa
  #undef R0i
  #undef Qci
  #undef alphai
  #undef mi
  #undef Wi
  #undef R0s
  #undef Qcs
  #undef alphas
  #undef ms
  #undef Ws
  #undef bk
  #undef mWin
  return(_comp);
} /* class_Bender_trace */

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
_class_Guide_channeled *class_Guide_channeled_trace(_class_Guide_channeled *_comp
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
  #define nslit (_comp->_parameters.nslit)
  #define d (_comp->_parameters.d)
  #define Qcx (_comp->_parameters.Qcx)
  #define Qcy (_comp->_parameters.Qcy)
  #define alphax (_comp->_parameters.alphax)
  #define alphay (_comp->_parameters.alphay)
  #define W (_comp->_parameters.W)
  #define mx (_comp->_parameters.mx)
  #define my (_comp->_parameters.my)
  #define nu (_comp->_parameters.nu)
  #define phase (_comp->_parameters.phase)
  #define w1c (_comp->_parameters.w1c)
  #define w2c (_comp->_parameters.w2c)
  #define ww (_comp->_parameters.ww)
  #define hh (_comp->_parameters.hh)
  #define whalf (_comp->_parameters.whalf)
  #define hhalf (_comp->_parameters.hhalf)
  #define lwhalf (_comp->_parameters.lwhalf)
  #define lhhalf (_comp->_parameters.lhhalf)
  double t1,t2;                                 /* Intersection times. */
  double av,ah,bv,bh,cv1,cv2,ch1,ch2,dd;        /* Intermediate values */
  double vdotn_v1,vdotn_v2,vdotn_h1,vdotn_h2;   /* Dot products. */
  int i;                                        /* Which mirror hit? */
  double q;                                     /* Q [1/AA] of reflection */
  double nlen2;                                 /* Vector lengths squared */
  double edge;
  double hadj;                                  /* Channel displacement */
  double angle=0;

  if (nu != 0 || phase != 0) { /* rotate neutron w/r to guide element */
    /* approximation of rotating straight Fermi Chopper */
    Coords   X = coords_set(x,y,z-l/2);  /* current coordinates of neutron in centered static frame */
    Rotation R;
    double dt=(-z+l/2)/vz; /* time shift to each center of slit package */
    angle=fmod(360*nu*(t+dt)+phase, 360); /* in deg */
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

  /* Propagate neutron to guide entrance. */
  PROP_Z0;
  /* Scatter here to ensure that fully transmitted neutrons will not be
     absorbed in a GROUP construction, e.g. all neutrons - even the
     later absorbed ones are scattered at the guide entry. */
  SCATTER;
  if(x <= w1/-2.0 || x >= w1/2.0 || y <= -hhalf || y >= hhalf)
    ABSORB;
  /* Shift origin to center of channel hit (absorb if hit dividing walls) */
  x += w1/2.0;
  edge = floor(x/w1c)*w1c;
  if(x - edge > w1c - d)
  {
    x -= w1/2.0; /* Re-adjust origin */
    ABSORB;
  }
  x -= (edge + (w1c - d)/2.0);
  hadj = edge + (w1c - d)/2.0 - w1/2.0;
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
        dd = 2*vdotn_v1/nlen2;
        vx = vx - dd*l;
        vz = vz - dd*ww;
        break;
      case 2:                   /* Right vertical mirror */
        nlen2 = l*l + ww*ww;
        q = V2Q*(-2)*vdotn_v2/sqrt(nlen2);
        dd = 2*vdotn_v2/nlen2;
        vx = vx + dd*l;
        vz = vz - dd*ww;
        break;
      case 3:                   /* Lower horizontal mirror */
        nlen2 = l*l + hh*hh;
        q = V2Q*(-2)*vdotn_h1/sqrt(nlen2);
        dd = 2*vdotn_h1/nlen2;
        vy = vy - dd*l;
        vz = vz - dd*hh;
        break;
      case 4:                   /* Upper horizontal mirror */
        nlen2 = l*l + hh*hh;
        q = V2Q*(-2)*vdotn_h2/sqrt(nlen2);
        dd = 2*vdotn_h2/nlen2;
        vy = vy + dd*l;
        vz = vz - dd*hh;
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
        double par[] = {R0, Qcx, alphax, mx, W};
        StdReflecFunc(q, par, &ref);
        if (ref > 0)
          p *= ref;
        else {
          x += hadj; /* Re-adjust origin */
          ABSORB;                               /* Cutoff ~ 1E-10 */
        }
      } else {
        double par[] = {R0, Qcy, alphay, my, W};
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
  } /* end for */
  x += hadj; /* Re-adjust origin */
  if (nu != 0 || phase != 0) { /* rotate back neutron w/r to guide element */
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
  #undef w1
  #undef h1
  #undef w2
  #undef h2
  #undef l
  #undef R0
  #undef Qc
  #undef alpha
  #undef m
  #undef nslit
  #undef d
  #undef Qcx
  #undef Qcy
  #undef alphax
  #undef alphay
  #undef W
  #undef mx
  #undef my
  #undef nu
  #undef phase
  #undef w1c
  #undef w2c
  #undef ww
  #undef hh
  #undef whalf
  #undef hhalf
  #undef lwhalf
  #undef lhhalf
  return(_comp);
} /* class_Guide_channeled_trace */

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
  if(!strcmp("DISCTOF", _comp->_name)) {
    user1 = FLT_MAX;
    user2 = FLT_MAX;
    user3 = FLT_MAX;
  }
  if(!strcmp("MONOTOF", _comp->_name)) {
    user1 = FLT_MAX;
    user2 = FLT_MAX;
    user3 = FLT_MAX;
  }
  if(!strcmp("FERMITOF_before", _comp->_name)) {
    user1 = FLT_MAX;
    user2 = FLT_MAX;
    user3 = FLT_MAX;
  }
  if(!strcmp("FERMITOF1", _comp->_name)) {
    user1 = FLT_MAX;
    user2 = FLT_MAX;
    user3 = FLT_MAX;
  }
  if(!strcmp("FERMITOF2", _comp->_name)) {
    user1 = FLT_MAX;
    user2 = FLT_MAX;
    user3 = FLT_MAX;
  }
  if(!strcmp("TOF_Det", _comp->_name)) {
    user1 = FLT_MAX;
    user2 = FLT_MAX;
    user3 = FLT_MAX;
  }
  if(!strcmp("FoDet", _comp->_name)) {
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
_class_FermiChopper *class_FermiChopper_trace(_class_FermiChopper *_comp
  , _class_particle *_particle) {
  ABSORBED=SCATTERED=RESTORE=0;

  #define phase (_comp->_parameters.phase)
  #define radius (_comp->_parameters.radius)
  #define nu (_comp->_parameters.nu)
  #define w (_comp->_parameters.w)
  #define nslit (_comp->_parameters.nslit)
  #define R0 (_comp->_parameters.R0)
  #define Qc (_comp->_parameters.Qc)
  #define alpha (_comp->_parameters.alpha)
  #define m (_comp->_parameters.m)
  #define W (_comp->_parameters.W)
  #define length (_comp->_parameters.length)
  #define eff (_comp->_parameters.eff)
  #define zero_time (_comp->_parameters.zero_time)
  #define xwidth (_comp->_parameters.xwidth)
  #define verbose (_comp->_parameters.verbose)
  #define yheight (_comp->_parameters.yheight)
  #define curvature (_comp->_parameters.curvature)
  #define delay (_comp->_parameters.delay)
  #define FCVars (_comp->_parameters.FCVars)

  /** local CALCULATION VARIABLES**************************************/

  /** Interaction with slit package ***************************/
  double slit_input; /* length of the slits */

  /** Variables for calculating interaction with blades ***************/
  double xp1, zp1, vxp1;

  /**  Reflections ***********************************************/
  double n1;

  /**  Multiple Reflections  ******************************/
  int loopcounter=0;   /* How many reflections happen? */

  /** Time variables *********************************/
  double t3=0;      /* interaction time at n3 position (slit wall) */
  double t1=0,t2=0; /* cylinder intersection time (at entry and exit of slit pack - n1 n2 - or cylinder) */
  double dt;        /* interaction intervals (e.g. time for exiting slit pack) */

  /************************ TIME OF FLIGHT RESET ************************/
  if (zero_time == 1 && nu)
    t -= floor( (t+1/(4*nu))*(2*nu) )/(2*nu);

  /* zero arrays used to store positions */
  for (loopcounter=0; loopcounter < 5; loopcounter++)
    FCVars.X[loopcounter] = FCVars.Z[loopcounter]=0;

  /************** test, if the neutron interacts with the cylinder ***/
  if (cylinder_intersect (&t1, &t2, x, y, z, vx, vy, vz, radius, yheight))
  {
    if (t1 <= 0) { /* Neutron started inside the cylinder */
      if (verbose > 0 && FCVars.absorb_alreadyinside<FermiChopper_MAXITER) printf("FermiChopper: %s: ABSORB Neutron started inside the cylinder, t1=%8.3g (enter).\n",
        NAME_CURRENT_COMP, t1);
      FCVars.absorb_alreadyinside++;
      ABSORB;
    }
    if (verbose > 2)
      printf("FermiChopper: %s:         t1=%8.3g t2=%8.3g xyz=[%8.3g %8.3g %8.3g] v=[%8.3g %8.3g %8.3g] t=%8.3g (init).\n",
           NAME_CURRENT_COMP, t1, t2, x,y,z,vx,vy,vz,t);

    dt=t2-t1;     /* total time of flight inside the cylinder  */
    PROP_DT(t1);  /* Propagates neutron to entrance of the cylinder */
    SCATTER;

    if (verbose > 2)
      printf("FermiChopper: %s: PROP_DT t1=%8.3g t2=%8.3g xyz=[%8.3g %8.3g %8.3g] v=[%8.3g %8.3g %8.3g] t=%8.3g (IN cyl).\n",
           NAME_CURRENT_COMP, t1, t2, x,y,z,vx,vy,vz,t);

    /* neutron must not enter or leave from top or bottom of cylinder. */
    if (fabs(y) >= yheight/2.0 || fabs(y+vy*dt) >= yheight/2.0) {
      if (verbose > 2) printf("FermiChopper: %s: ABSORB Neutron hits top/bottom of cylinder, y=%8.3g (enter).\n",
        NAME_CURRENT_COMP, y);
      FCVars.absorb_topbottom++; ABSORB;
    }

    vxp1 = sqrt(vx*vx+vy*vy+vz*vz);
    FCVars.sum_v += p*vxp1;
    FCVars.sum_t += p*t;
    FCVars.sum_N += p;

    if (zero_time > 1 && FCVars.sum_N) { /* automatic phase adjustment */
      double mean_t, mean_phase;
      mean_t     = FCVars.sum_t/FCVars.sum_N;
      mean_phase = fmod(mean_t*nu*2*PI, 2*PI);
      /* now we shift the phase so that the neutron pulse is centered on the slit pack */
      mean_phase+= radius/vxp1*2*PI*nu;
      FCVars.ph0 = mean_phase;
    }

    /* neutron must enter the cylinder opening: |X'| < full package width*/
    xp1 = FC_xrot(x,z, t,FCVars); /* X'(t) */
    if (fabs(xp1) >= nslit*w/2) {
      if (verbose > 2)
        printf("FermiChopper: %s: ABSORB Neutron X is outside cylinder aperture, x'=%8.3g > %g (enter).\n",
        NAME_CURRENT_COMP, xp1, nslit*w/2);
      FCVars.absorb_cylentrance++; ABSORB;
    }

/*********************** PROPAGATE TO SLIT PACKAGE **************************/

    /* zp1 = Z' at entrance of cylinder Z'(t) */
    zp1  = FC_zrot(x,z, t, FCVars);

    FCVars.X[0] = xp1; FCVars.Z[0] = zp1;

    /* here we should have sqrt(x*x+z*z) == sqrt(xp1*xp1+zp1*zp1) */
    t3 = sqrt(x*x+z*z) / sqrt(xp1*xp1+zp1*zp1);
    if (t3 < 0.99 || 1.01 < t3) {
      if (verbose > 1 && FCVars.absorb_cylentrance < FermiChopper_MAXITER)
        printf("FermiChopper: %s: ABSORB Neutron radius on cylinder in static frame r=%g does not match that of rotating frame r'=%g.\n",
        NAME_CURRENT_COMP, sqrt(x*x+z*z), sqrt(xp1*xp1+zp1*zp1));
      FCVars.absorb_cylentrance++; ABSORB;
    }

    /* Checking on which side of the Chopper the Neutron enters: sign(Z') */
    slit_input = (zp1 > 0 ? length/2 : -length/2);

    /* time shift to reach slit package in [0,time to exit cylinder]: Z'=slit_input */
    /* t3 is used here as a tmp variable, will be redefined in for loop  */
    t3 = FC_zintersect(x,z,vx,vz, t,dt, slit_input, FCVars);

    if( (t3 < 0)||(t3 > dt) ) {
      if (verbose > 2 && FCVars.absorb_notreachentrance < FermiChopper_MAXITER) {
        printf("FermiChopper: %s: Can not reach entrance of slits. dt=%8.3g t3=%8.3g (intersection:1).\n",
        NAME_CURRENT_COMP, dt, t3);
        if (t3 == -3)
            printf("          No sign change to determine intersection\n");
        else if (t3 == -2)
            printf("          Max iterations reached\n");
        else if (t3 < 0)
            printf("          Error when solving intersection\n");
      }
      FCVars.absorb_notreachentrance++;
      ABSORB; /* neutron can not reach slit entrance */
    }

    /* Propagating to the slit package entrance */
    PROP_DT(t3); /* dt = t2-t1: time in cylinder */
    dt -= t3; /* remaining time from slit pack entry to exit of cylinder */
    xp1 = FC_xrot(x,z, t, FCVars); /* should be slit_input */
    zp1 = FC_zrot(x,z, t, FCVars);
    FCVars.X[1] = xp1; FCVars.Z[1] = zp1;

    if (mcdotrace) {
      /* indicate position of neutron in mcdisplay */
      double xp2 = x; double zp2 = z; x = xp1; z=zp1; SCATTER; x=xp2; z=zp2;
    } else SCATTER;

    if (verbose > 2)
      printf("FermiChopper: %s: PROP_DT t=%8.3g dt=%8.3g x'=%8.3g z'=%8.3g length=%g (slit enter).\n",
           NAME_CURRENT_COMP, t, dt, xp1, zp1, slit_input);

    /* must have X'< slit package width at package Z'=slit_input */
    if (fabs(xp1) >= nslit*w/2) {
      if (verbose > 2) printf("FermiChopper: %s: ABSORB Neutron X is outside slit package, x'=%8.3g > %g (enter).\n",
        NAME_CURRENT_COMP, xp1, nslit*w/2);
      FCVars.absorb_packentrance++; ABSORB;
    }

    /* solve Z'=-slit_input for time of exit of slit package */
    /* t3 is used here as a tmp variable, will be redefined in for loop  */
    t3 = FC_zintersect(x,z,vx,vz, t,dt*1.1, -slit_input, FCVars);

    if((t3 < FermiChopper_TimeAccuracy)||(t3 > dt)) {
      if (verbose > 1 && FCVars.warn_notreachslitoutput < FermiChopper_MAXITER) {
        printf("FermiChopper: %s: Can not reach exit of slits. dt=%8.3g t3=%8.3g (intersection:2).\n",
          NAME_CURRENT_COMP, dt, t3);
        if (t3 == -3)
            printf("          No sign change to determine intersection\n");
        else if (t3 == -2)
            printf("          Max iterations reached\n");
        else if (t3 < 0)
            printf("          Error when solving intersection\n");
      }
      FCVars.warn_notreachslitoutput++;
      /* we estimate analytically the time to the slit exit */
      t3 = length/vxp1;
    }
    dt = t3; /* reduce time interval to [0, time of slit exit] */

    /* here we should have dt*v = length (slit entrance -> exit) */
    t3 = fabs(FC_xzrot_dt(x,z,vx,vz,    t,   0 , 'z', FCVars) - FC_xzrot_dt(x,z,vx,vz,    t,   dt, 'z', FCVars))/length;
    if (fabs(t3-1) > 0.02) {
      if (verbose > 0 && FCVars.warn_notreachslitoutput < FermiChopper_MAXITER)
        printf("FermiChopper: %s: ABSORB Neutron propagation time v*dt/length=%g in slit does not match its length=%g (slit exit expected).\n",
          NAME_CURRENT_COMP, t3, length);
      FCVars.warn_notreachslitoutput++;
      ABSORB;
    }

    /* dt= time shift to go to exit of slit package (or exit of cylinder in case of error) */
    /*
      |---------|
      | /       |
      |o        * (dt)
      |         |
      |---------|
     */

/*********************PROPAGATION INSIDE THE SLIT PACKAGE *******************/

    /* index of slit hit at entrance n1=-N/2 (xp1=-) ... N/2 (xp1=+) */
    n1 = floor(xp1/w);

/******************* BEGIN LOOP FOR MULTIPLE REFLECTIONS ********************/

    for (loopcounter=0; loopcounter<=FermiChopper_MAXITER; loopcounter++) {
      double dt_to_tangent=0; /* time shift to go to tangent intersection */
      double n2,n3;           /* slit indices */
      double xp2, zp2,xp3=0,vxp2=0,vzp1=0; /* position, velocity */
      double q;                   /* used for calculating new velocities after reflection */
      int    i;

      /* compute trajectory tangents: m1=Vz'+w.X'(t), m2=Vz'+w.X'(t+dt) */
      xp1 = FC_xrot    (x,z,          t,   FCVars);          /* X'(t)    current position */
      xp2 = FC_xzrot_dt(x,z,vx,vz,    t,   dt, 'x', FCVars); /* X'(t+dt) slit exit */
      zp2 = FC_xzrot_dt(x,z,vx,vz,    t,   dt, 'z', FCVars); /* Z'(t+dt) slit exit */

      /* slit index at the end of the slit: */
      n2 = floor(xp2/w);

      /* quick exit for absorbing walls when neutron changes slit index */
      if (n2 != n1 && (m <= 0 || R0 <= 0 || Qc <= 0)) {
        if (verbose > 2) printf("FermiChopper: %s: ABSORB Neutron hits absorbing coating (change slit).\n",
            NAME_CURRENT_COMP);
        FCVars.absorb_slitcoating++; ABSORB;
      }

      /* compute transversal velocity to determine their intersection */
      vxp1= FC_xrot    (vx+z*FCVars.omega,vz-x*FCVars.omega,
                                      t,   FCVars);          /* dX'(t)/dt slope at current position*/

      vxp2= FC_xrot    (vx+(z+vz*dt)*FCVars.omega,vz-(x+vx*dt)*FCVars.omega,
                                      t+dt,FCVars);          /* dX'(t+dt)/dt slope at slit exit */

      /* absolute time at tangent intersection, changed to time shift below */
      dt_to_tangent = (vxp1 - vxp2 ? (xp2 - xp1 - dt*vxp2)/(vxp1 - vxp2) : -1);

      /* If algorithm fails, take the middle of the interval*/
      if (dt_to_tangent < 0 || dt_to_tangent > dt) {
        if (verbose > 1)
          printf("FermiChopper: %s: WARNING could not determine tangent intersection dt=%g. Using middle.\n",
          NAME_CURRENT_COMP, dt_to_tangent);
        dt_to_tangent=dt*0.5;
      }

      /*
           *(dt_to_tangent, xp3)
      |---------|
      | /     \ |
  xp1 |o       \|(dt) xp2
      |         |
      |---------|
     */

      /* point coordinates at tangent intersection/middle point (max deviation from optical axis) */
      xp3 = FC_xzrot_dt(x,z,vx,vz, t, dt_to_tangent, 'x', FCVars); /* X'(t+dt_to_tangent) */

      /* slit index at the tangent intersection/middle point */
      n3 = floor(xp3/w);

      if (verbose > 2)
        printf("FermiChopper: %s: t3=%8.3g slit_indices=[%g %g %g] (time at tangent intersection).\n",
           NAME_CURRENT_COMP, dt_to_tangent, n1, n2, n3);

      /* change slit means there is a reflection/intersection inside */
      if ( n2!=n1 || n3!=n1 ) {

        double t3a, t3b, distance_Wa, distance_Wb;
        if (m <= 0 || R0 <= 0 || Qc <= 0) {
          if (verbose > 2) printf("FermiChopper: %s: ABSORB Neutron hits absorbing coating (change slit).\n",
            NAME_CURRENT_COMP);
          FCVars.absorb_slitcoating++; ABSORB;
        }

        /* Choosing the first time it isn't in the slit anymore */
        if(n2 == n1 && n3 != n1){
          n2 = n3;
        }

        /* get position of slit wall towards which neutron is propagating */
        if (n2 > n1) {     /* X' positive side of slit is in principle the first intersection to test*/
          distance_Wa = n1*w+w;
          distance_Wb = n1*w;
        } else {            /* X' negative side of slit */
          distance_Wb = n1*w+w;
          distance_Wa = n1*w;
        }

        /* time shift to reach slit wall point in [0,dt_to_tangent]: X'=distance_W slit_wall */
        for (i=0; i< 2; i++) {
          /* first  attempt: [0,dt_to_tangent] (to max deviation location)
             second attempt: [0, dt]           (to slit exit) */

          double dt_search = (i == 0 ? dt_to_tangent : dt);
          t3a = FC_xintersect(x,z,vx,vz, t,dt_to_tangent, distance_Wa, FCVars);
          t3b = FC_xintersect(x,z,vx,vz, t,dt_to_tangent, distance_Wb, FCVars);
          if      (t3b < 0)             t3 = t3a;
          else if (t3a < 0 && t3b >= 0) t3 = t3b;
          else                          t3 = (t3a < t3b ? t3a : t3b);
          if (FermiChopper_TimeAccuracy < t3 && t3 < dt_search)
            break; /* we found the intersection location */
        }
        /* handle case where intersection search fails */
        if (t3 < FermiChopper_TimeAccuracy || t3 >= dt_to_tangent) {
          if (verbose > 1 && FCVars.warn_notreachslitwall < FermiChopper_MAXITER) {
            printf("FermiChopper: %s: Can not reach slit wall (iteration %i). dt=%8.3g t3=%8.3g (intersection:3).\n",
              NAME_CURRENT_COMP, loopcounter, dt_to_tangent, t3);
            if (t3 == -3)
              printf("        No sign change to determine intersection\n");
            else if (t3 == -2)
              printf("        Max iterations reached\n");
            else if (t3 < 0 || t3 >= dt_to_tangent)
              printf("        Error when solving intersection\n");
          }
          /* neutron can not reach slit wall. */
          FCVars.warn_notreachslitwall++;
          ABSORB;
        }

        /* Propagate to slit wall point (t+t3) on slit n3 wall */
        PROP_DT(t3); dt -= t3; /* dt: time remaining to slit exit after propagation */
        xp1 = FC_xrot(x,z, t, FCVars); /* X'(t+t3) : on slit wall */
        zp1 = FC_zrot(x,z, t, FCVars); /* Z'(t+t3) : on slit wall */
        FCVars.X[2] = xp1; FCVars.Z[2] = zp1;

        if (verbose > 2)
          printf("FermiChopper: %s: PROP_DT t3=%8.3g dt=%8.3g xyz=[%8.3g %8.3g %8.3g] (on wall). z'=%g\n",
           NAME_CURRENT_COMP, t3, dt, x,y,z, zp1);

        /* check if intersection point is still in the slit package, else exit loop */
        if (fabs(zp1) > length/2 || dt <= 0) {
          if (verbose > 2)
            printf("FermiChopper: %s: Neutron is outside slit pack (on slit wall).\n",
                NAME_CURRENT_COMP);
          break;
        }

    /*    here
      |---o-----|
      | /   \   |
      |/     \  |
      |       \ |(dt)
      |---------|
     */

        /* get velocity in rotating frame, on slit wall */
        vxp1 = FC_xrot(vx,vz, t, FCVars);
        vzp1 = FC_zrot(vx,vz, t, FCVars);

        q    = 2*V2Q*(fabs(vxp1));

        {
          double ref = 0;
          double par[] = {R0, Qc, alpha, m, W};
          StdReflecFunc(q, par, &ref);
          if (ref>0) p *= ref;
          else {
            if (verbose > 2) printf("FermiChopper: %s: ABSORB Neutron hits absorbing coating (on slit wall).\n",
              NAME_CURRENT_COMP);
            FCVars.absorb_slitcoating++; ABSORB;
          } /* Cutoff ~ 1E-10 */
        }
        if (mcdotrace) {
          double xp2 = x; double zp2 = z;
          /* indicate position of neutron in mcdisplay */
          x = FC_xrot(x,z, t,FCVars); z= FC_zrot(x,z, t,FCVars); SCATTER; x=xp2; z=zp2;
        } else SCATTER;

        /* reflect perpendicular velocity and compute new velocity in static frame */
        vxp1 *= -1;
        /* apply transposed Transformation matrix */
        vx = FC_xrot( vxp1,-vzp1, t,FCVars);
        vz = FC_zrot(-vxp1, vzp1, t,FCVars);

        /* recompute time to slit exit */
        /* solve Z'=-slit_input for time of exit of slit package */
        t3 = FC_zintersect(x,z,vx,vz, t,dt, -slit_input, FCVars);

        if(t3 < 0 || t3 > dt) {
          if (verbose > 1 && FCVars.warn_notreachslitoutput < FermiChopper_MAXITER) {
            printf("FermiChopper: %s: Can not reach exit of slits. dt=%8.3g t3=%8.3g (intersection:4).\n",
              NAME_CURRENT_COMP, dt, t3);
            if (t3 == -3)
              printf("              No sign change to determine intersection\n");
            else if (t3 == -2)
              printf("              Max iterations reached\n");
            else if (t3 < 0)
              printf("              Error when solving intersection\n");
          }
          FCVars.warn_notreachslitoutput++;
          ABSORB;
          /* neutron can not reach slit output. */
        } else dt = t3; /* reduce time interval to [0, time of slit exit] */

      } /* end if changed slit index */
      else { /* neutron remains in the same slit: direct flight */
        if (dt > 0) PROP_DT(dt); /* go to slit exit */
        break;
      }
    } /* end for */

    xp1 = FC_xrot(x,z, t,FCVars);
    zp1 = FC_zrot(x,z, t,FCVars);
    FCVars.X[3] = xp1; FCVars.Z[3] = zp1; /* slit exit */

    if (fabs(xp1) >= nslit*w/2)
    {
      if (verbose > 2) printf("FermiChopper: %s: ABSORB Neutron X is outside slit package, x'=%8.3g (slit exit).\n",
        NAME_CURRENT_COMP, xp1);
      FCVars.absorb_slitcoating++; ABSORB;
    }

    if (loopcounter >= FermiChopper_MAXITER)
    {
      if (verbose > 1 && FCVars.absorb_maxiterations < FermiChopper_MAXITER)
      printf("FermiChopper: %s: Max iterations %i reached inside slit. Absorb.\n",
        NAME_CURRENT_COMP, FermiChopper_MAXITER);
      FCVars.absorb_maxiterations++;
      ABSORB;
    }


    /********************* EXIT SLIT PACKAGE ********************************/

    /* propagation times to cylinder. Should use t2 to exit */
    if (!cylinder_intersect (&t1, &t2, x, y, z, vx, vy, vz, radius, yheight)) {
      /* must be inside the cylinder */
      if (verbose > 1) printf("FermiChopper: %s: ABSORB Neutron has unexpectidely exited cylinder ! (exiting)\n",
        NAME_CURRENT_COMP);
      FCVars.absorb_exitslitpack++;
      ABSORB; }

    if (t1 > 0)
    {
      if (verbose > 1 && FCVars.absorb_wrongdirection < FermiChopper_MAXITER)
      printf("FermiChopper: %s: Neutrons are leaving chopper\n"
             "              in the wrong direction. Absorb.\n", NAME_CURRENT_COMP);
      FCVars.absorb_wrongdirection++;
      ABSORB;
    }

    if (t2 <= 0 && FCVars.absorb_nocontrol < FermiChopper_MAXITER)
    {
      if (verbose > 1)
      printf("FermiChopper: %s: Neutrons are leaving chopper\n"
             "              without any control. Absorb.\n", NAME_CURRENT_COMP);
      FCVars.absorb_nocontrol++;
      ABSORB;
    }

    /* propagate to cylinder surface (exit) */
    PROP_DT(t2);
    SCATTER;

    xp1 = FC_xrot(x,z, t,FCVars);
    zp1 = FC_zrot(x,z, t,FCVars);
    FCVars.X[4] = xp1; FCVars.Z[4] = zp1;

    if (verbose > 2)
      printf("FermiChopper: %s: t1=%8.3g PROP_DT t2=%8.3g xyz=[%8.3g %8.3g %8.3g] (OUT cyl). z'=%g\n",
           NAME_CURRENT_COMP, t1, t2, x,y,z, zp1);

    /* Check if the neutron left the cylinder by its top or bottom */
    if (fabs(y) >= yheight/2)
    {
      if (verbose > 2) printf("FermiChopper: %s: ABSORB Neutron hits top/bottom of cylinder, y=%8.3g (exiting)\n",
        NAME_CURRENT_COMP, y);
      FCVars.absorb_topbottom++; ABSORB;
    }

    /* must have X'< slit package width at package Z'=cylinder output */
    if (fabs(xp1) >= nslit*w/2)
    {
      if (verbose > 2) printf("FermiChopper: %s: ABSORB Neutron X is outside slit package cylinder, xp1=%8.3g (exiting).\n",
        NAME_CURRENT_COMP, xp1);
      FCVars.absorb_cylexit++; ABSORB;
    }

    /* Transmission coefficent */
    p = p*eff;          //finite cross section + transmission

    FCVars.sum_N_pass += p;

  } /* end if cylinder_intersect */

  #undef phase
  #undef radius
  #undef nu
  #undef w
  #undef nslit
  #undef R0
  #undef Qc
  #undef alpha
  #undef m
  #undef W
  #undef length
  #undef eff
  #undef zero_time
  #undef xwidth
  #undef verbose
  #undef yheight
  #undef curvature
  #undef delay
  #undef FCVars
  return(_comp);
} /* class_FermiChopper_trace */

#pragma acc routine seq
_class_V_sample *class_V_sample_trace(_class_V_sample *_comp
  , _class_particle *_particle) {
  ABSORBED=SCATTERED=RESTORE=0;

  #define radius (_comp->_parameters.radius)
  #define thickness (_comp->_parameters.thickness)
  #define zdepth (_comp->_parameters.zdepth)
  #define Vc (_comp->_parameters.Vc)
  #define sigma_abs (_comp->_parameters.sigma_abs)
  #define sigma_inc (_comp->_parameters.sigma_inc)
  #define radius_i (_comp->_parameters.radius_i)
  #define radius_o (_comp->_parameters.radius_o)
  #define h (_comp->_parameters.h)
  #define focus_r (_comp->_parameters.focus_r)
  #define pack (_comp->_parameters.pack)
  #define frac (_comp->_parameters.frac)
  #define f_QE (_comp->_parameters.f_QE)
  #define gamma (_comp->_parameters.gamma)
  #define target_x (_comp->_parameters.target_x)
  #define target_y (_comp->_parameters.target_y)
  #define target_z (_comp->_parameters.target_z)
  #define focus_xw (_comp->_parameters.focus_xw)
  #define focus_yh (_comp->_parameters.focus_yh)
  #define focus_aw (_comp->_parameters.focus_aw)
  #define focus_ah (_comp->_parameters.focus_ah)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define zthick (_comp->_parameters.zthick)
  #define rad_sphere (_comp->_parameters.rad_sphere)
  #define sig_a (_comp->_parameters.sig_a)
  #define sig_i (_comp->_parameters.sig_i)
  #define V0 (_comp->_parameters.V0)
  #define target_index (_comp->_parameters.target_index)
  #define multiples (_comp->_parameters.multiples)
  #define VarsV (_comp->_parameters.VarsV)
  double t0, t3;                /* Entry/exit time for outer cylinder */
  double t1, t2;                /* Entry/exit time for inner cylinder */
  double v;                     /* Neutron velocity */
  double dt0, dt1, dt2, dt;     /* Flight times through sample */
  double l_full;                /* Flight path length for non-scattered neutron */
  double l_i, l_o=0;            /* Flight path lenght in/out for scattered neutron */
  double my_a=0;                  /* Velocity-dependent attenuation factor */
  double solid_angle=0;         /* Solid angle of target as seen from scattering point */
  double aim_x=0, aim_y=0, aim_z=1;   /* Position of target relative to scattering point */
  double v_i, v_f, E_i, E_f; /* initial and final energies and velocities */
  double dE;                 /* Energy transfer */
  int    intersect=0;

  if (VarsV.shapetyp == 2)
    intersect = sphere_intersect(&t0, &t3, x, y, z, vx, vy, vz, rad_sphere);
  else
    if (VarsV.shapetyp == 1)
      intersect = box_intersect(&t0, &t3, x, y, z, vx, vy, vz, xwidth, yheight, zthick);
  else
    intersect = cylinder_intersect(&t0, &t3, x, y, z, vx, vy, vz, radius_o, h);
  if(intersect)
  {
    if(t0 < 0) ABSORB; /* we already passed the sample; this is illegal */
    /* Neutron enters at t=t0. */
    if(VarsV.shapetyp == 1 || VarsV.shapetyp == 2)
      t1 = t2 = t3;
    else
      if(!radius_i || !cylinder_intersect(&t1, &t2, x, y, z, vx, vy, vz, radius_i, h))
        t1 = t2 = t3;

    dt0 = t1-t0;                /* Time in sample, ingoing */
    dt1 = t2-t1;                /* Time in hole */
    dt2 = t3-t2;                /* Time in sample, outgoing */
    v = sqrt(vx*vx + vy*vy + vz*vz);
    l_full = v * (dt0 + dt2);   /* Length of full path through sample */
    if (v) my_a = VarsV.my_a_v*(2200/v);

    if (frac >= 1 || rand01()<frac)          /* Scattering */
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

      if(VarsV.shapetyp == 0) {
        if(!cylinder_intersect(&t0, &t3, x, y, z, vx, vy, vz, radius_o, h)) {
          /* ??? did not hit cylinder */
          printf("FATAL ERROR: Did not hit cylinder from inside.\n");
          exit(1);
        }
        dt = t3; /* outgoing point */
        if(cylinder_intersect(&t1, &t2, x, y, z, vx, vy, vz, radius_i, h) &&
           t2 > 0)
          dt -= (t2-t1);            /* Subtract hollow part */
      }
      else {
        if(VarsV.shapetyp == 1) {
	      if(!box_intersect(&t0, &t3, x, y, z, vx, vy, vz, xwidth, yheight, zthick)) {
            /* ??? did not hit box */
            printf("FATAL ERROR: Did not hit box from inside.\n");
            exit(1);
          }
          dt = t3;
        }
        else {
	      if(!sphere_intersect(&t0, &t3, x, y, z, vx, vy, vz, rad_sphere)) {
            /* ??? did not hit sphere */
            printf("FATAL ERROR: Did not hit sphere from inside.\n");
            exit(1);
          }
          dt = t3;  
        }
      }
      l_o = v*dt; /* trajectory after scattering point: absorption only */

      p *= v/v_i*l_full*VarsV.my_s*exp(-my_a*(l_i+v_i/v*l_o)-VarsV.my_s*l_i);
      if (!multiples) {
	/* If no "multiples", correct by applying scattering cross-sec and
	   implicitly "absorb" further scattering (as in PowderN) 
	   We are currently (august 2007) having a debate on which solution 
	   is the most reasonable */
	p *= exp(-VarsV.my_s*l_o);
      }
      /* We do not consider scattering from 2nd part (outgoing) */
      p /= 4*PI/solid_angle;
      p /= frac;

      /* Polarisation part (1/3 NSF, 2/3 SF) */
      sx *= -1.0/3.0;
      sy *= -1.0/3.0;
      sz *= -1.0/3.0;

      SCATTER;
    }
    else /* Transmitting; always elastic */
    {
      p *= exp(-(my_a+VarsV.my_s)*l_full);
      p /= (1-frac);
    }
  }
  #undef radius
  #undef thickness
  #undef zdepth
  #undef Vc
  #undef sigma_abs
  #undef sigma_inc
  #undef radius_i
  #undef radius_o
  #undef h
  #undef focus_r
  #undef pack
  #undef frac
  #undef f_QE
  #undef gamma
  #undef target_x
  #undef target_y
  #undef target_z
  #undef focus_xw
  #undef focus_yh
  #undef focus_aw
  #undef focus_ah
  #undef xwidth
  #undef yheight
  #undef zthick
  #undef rad_sphere
  #undef sig_a
  #undef sig_i
  #undef V0
  #undef target_index
  #undef multiples
  #undef VarsV
  return(_comp);
} /* class_V_sample_trace */

/* *****************************************************************************
* instrument 'PSI_Focus' TRACE
***************************************************************************** */

#pragma acc routine seq
int raytrace(_class_particle* _particle) { /* called by mccode_main for PSI_Focus:TRACE */

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
      /* component a1=Progress_bar() [1] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_a1->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _a1->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_a1->_position_relative, _a1->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Progress_bar_trace(_a1, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component a1 [1] */
    if (!ABSORBED && _particle->_index == 2) {
      /* component csource=Source_gen() [2] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_csource->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _csource->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_csource->_position_relative, _csource->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Source_gen_trace(_csource, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component csource [2] */
    if (!ABSORBED && _particle->_index == 3) {
      /* component guide1=Guide() [3] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_guide1->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _guide1->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_guide1->_position_relative, _guide1->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Guide_trace(_guide1, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component guide1 [3] */
    if (!ABSORBED && _particle->_index == 4) {
      /* component guide2=Bender() [4] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_guide2->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _guide2->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_guide2->_position_relative, _guide2->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Bender_trace(_guide2, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component guide2 [4] */
    if (!ABSORBED && _particle->_index == 5) {
      /* component bunker=Guide() [5] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_bunker->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _bunker->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_bunker->_position_relative, _bunker->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Guide_trace(_bunker, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component bunker [5] */
    if (!ABSORBED && _particle->_index == 6) {
      /* component guide3=Guide() [6] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_guide3->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _guide3->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_guide3->_position_relative, _guide3->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Guide_trace(_guide3, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component guide3 [6] */
    if (!ABSORBED && _particle->_index == 7) {
      /* component lambdaGuideExit=L_monitor() [7] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_lambdaGuideExit->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _lambdaGuideExit->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_lambdaGuideExit->_position_relative, _lambdaGuideExit->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_L_monitor_trace(_lambdaGuideExit, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component lambdaGuideExit [7] */
    if (!ABSORBED && _particle->_index == 8) {
      /* component DivMonGuideExit=Divergence_monitor() [8] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_DivMonGuideExit->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _DivMonGuideExit->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_DivMonGuideExit->_position_relative, _DivMonGuideExit->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Divergence_monitor_trace(_DivMonGuideExit, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component DivMonGuideExit [8] */
    if (!ABSORBED && _particle->_index == 9) {
      /* component PSDGuideExit=PSD_monitor() [9] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_PSDGuideExit->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _PSDGuideExit->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_PSDGuideExit->_position_relative, _PSDGuideExit->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_PSD_monitor_trace(_PSDGuideExit, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component PSDGuideExit [9] */
    if (!ABSORBED && _particle->_index == 10) {
      /* component FOCUSguide=Guide_channeled() [10] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_FOCUSguide->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _FOCUSguide->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_FOCUSguide->_position_relative, _FOCUSguide->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Guide_channeled_trace(_FOCUSguide, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component FOCUSguide [10] */
    if (!ABSORBED && _particle->_index == 11) {
      /* component FirstChopper=DiskChopper() [11] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_FirstChopper->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _FirstChopper->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_FirstChopper->_position_relative, _FirstChopper->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_DiskChopper_trace(_FirstChopper, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component FirstChopper [11] */
    if (!ABSORBED && _particle->_index == 12) {
      /* component DISCTOF=Monitor_nD() [12] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_DISCTOF->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _DISCTOF->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_DISCTOF->_position_relative, _DISCTOF->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Monitor_nD_trace(_DISCTOF, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component DISCTOF [12] */
    if (!ABSORBED && _particle->_index == 13) {
      /* component PSDmon1Chopper=PSD_monitor() [13] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_PSDmon1Chopper->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _PSDmon1Chopper->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_PSDmon1Chopper->_position_relative, _PSDmon1Chopper->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_PSD_monitor_trace(_PSDmon1Chopper, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component PSDmon1Chopper [13] */
    if (!ABSORBED && _particle->_index == 14) {
      /* component VacuumTube1entry=Slit() [14] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_VacuumTube1entry->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _VacuumTube1entry->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_VacuumTube1entry->_position_relative, _VacuumTube1entry->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Slit_trace(_VacuumTube1entry, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component VacuumTube1entry [14] */
    if (!ABSORBED && _particle->_index == 15) {
      /* component VacuumTube1exit=Slit() [15] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_VacuumTube1exit->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _VacuumTube1exit->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_VacuumTube1exit->_position_relative, _VacuumTube1exit->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Slit_trace(_VacuumTube1exit, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component VacuumTube1exit [15] */
    if (!ABSORBED && _particle->_index == 16) {
      /* component VacuumTube2entry=Slit() [16] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_VacuumTube2entry->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _VacuumTube2entry->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_VacuumTube2entry->_position_relative, _VacuumTube2entry->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Slit_trace(_VacuumTube2entry, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component VacuumTube2entry [16] */
    if (!ABSORBED && _particle->_index == 17) {
      /* component VacuumTube2exit=Slit() [17] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_VacuumTube2exit->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _VacuumTube2exit->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_VacuumTube2exit->_position_relative, _VacuumTube2exit->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Slit_trace(_VacuumTube2exit, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component VacuumTube2exit [17] */
    if (!ABSORBED && _particle->_index == 18) {
      /* component VacuumTube3entry=Slit() [18] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_VacuumTube3entry->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _VacuumTube3entry->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_VacuumTube3entry->_position_relative, _VacuumTube3entry->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Slit_trace(_VacuumTube3entry, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component VacuumTube3entry [18] */
    if (!ABSORBED && _particle->_index == 19) {
      /* component VacuumTube3exit=Slit() [19] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_VacuumTube3exit->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _VacuumTube3exit->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_VacuumTube3exit->_position_relative, _VacuumTube3exit->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Slit_trace(_VacuumTube3exit, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component VacuumTube3exit [19] */
    if (!ABSORBED && _particle->_index == 20) {
      /* component PSDmonMono=PSD_monitor() [20] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_PSDmonMono->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _PSDmonMono->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_PSDmonMono->_position_relative, _PSDmonMono->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_PSD_monitor_trace(_PSDmonMono, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component PSDmonMono [20] */
    if (!ABSORBED && _particle->_index == 21) {
      /* component MONOTOF=Monitor_nD() [21] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_MONOTOF->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _MONOTOF->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_MONOTOF->_position_relative, _MONOTOF->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Monitor_nD_trace(_MONOTOF, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component MONOTOF [21] */
    if (!ABSORBED && _particle->_index == 22) {
      /* component DivMonMono=Divergence_monitor() [22] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_DivMonMono->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _DivMonMono->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_DivMonMono->_position_relative, _DivMonMono->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Divergence_monitor_trace(_DivMonMono, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component DivMonMono [22] */
    if (!ABSORBED && _particle->_index == 23) {
      /* component focus_mono=Arm() [23] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_focus_mono->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _focus_mono->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_focus_mono->_position_relative, _focus_mono->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      _particle->_index++;
    } /* end component focus_mono [23] */
    /* start SPLIT at mono */
    instrument->logic.Split_mono_particle=_particle;
    if (!ABSORBED && _particle->_index == 24) 
    for (instrument->logic.Split_mono = 0; instrument->logic.Split_mono< 10; instrument->logic.Split_mono++) {
    {
      /* component mono=Monochromator_2foc() [24] */
      _particle=instrument->logic.Split_mono_particle;
      p /= 10 > 0 ? 10 : 1;
  #define reflect (_mono->_parameters.reflect)
  #define zwidth (_mono->_parameters.zwidth)
  #define yheight (_mono->_parameters.yheight)
  #define gap (_mono->_parameters.gap)
  #define NH (_mono->_parameters.NH)
  #define NV (_mono->_parameters.NV)
  #define mosaich (_mono->_parameters.mosaich)
  #define mosaicv (_mono->_parameters.mosaicv)
  #define r0 (_mono->_parameters.r0)
  #define Q (_mono->_parameters.Q)
  #define RV (_mono->_parameters.RV)
  #define RH (_mono->_parameters.RH)
  #define DM (_mono->_parameters.DM)
  #define mosaic (_mono->_parameters.mosaic)
  #define width (_mono->_parameters.width)
  #define height (_mono->_parameters.height)
  #define verbose (_mono->_parameters.verbose)
  #define mos_y (_mono->_parameters.mos_y)
  #define mos_z (_mono->_parameters.mos_z)
  #define mono_Q (_mono->_parameters.mono_Q)
  #define SlabWidth (_mono->_parameters.SlabWidth)
  #define SlabHeight (_mono->_parameters.SlabHeight)
  #define rTable (_mono->_parameters.rTable)
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_mono->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _mono->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_mono->_position_relative, _mono->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Monochromator_2foc_trace(_mono, _particle);
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
    } /* end component mono [24] */
    if (!ABSORBED && _particle->_index == 25) {
      /* component a2=Arm() [25] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_a2->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _a2->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_a2->_position_relative, _a2->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      _particle->_index++;
    } /* end component a2 [25] */
    if (!ABSORBED && _particle->_index == 26) {
      /* component FERMITOF_before=Monitor_nD() [26] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_FERMITOF_before->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _FERMITOF_before->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_FERMITOF_before->_position_relative, _FERMITOF_before->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Monitor_nD_trace(_FERMITOF_before, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component FERMITOF_before [26] */
    if (!ABSORBED && _particle->_index == 27) {
      /* component lambdaFermi=L_monitor() [27] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_lambdaFermi->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _lambdaFermi->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_lambdaFermi->_position_relative, _lambdaFermi->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_L_monitor_trace(_lambdaFermi, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component lambdaFermi [27] */
    if (!ABSORBED && _particle->_index == 28) {
      /* component EMON_Fermi=E_monitor() [28] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_EMON_Fermi->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _EMON_Fermi->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_EMON_Fermi->_position_relative, _EMON_Fermi->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_E_monitor_trace(_EMON_Fermi, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component EMON_Fermi [28] */
    if (!ABSORBED && _particle->_index == 29) {
      /* component DivMonfermi1=Divergence_monitor() [29] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_DivMonfermi1->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _DivMonfermi1->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_DivMonfermi1->_position_relative, _DivMonfermi1->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Divergence_monitor_trace(_DivMonfermi1, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component DivMonfermi1 [29] */
    if (!ABSORBED && _particle->_index == 30) {
      /* component PSD_Fermi1=PSD_monitor() [30] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_PSD_Fermi1->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _PSD_Fermi1->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_PSD_Fermi1->_position_relative, _PSD_Fermi1->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_PSD_monitor_trace(_PSD_Fermi1, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component PSD_Fermi1 [30] */
    if (!ABSORBED && _particle->_index == 31) {
      /* component FoChopper=FermiChopper() [31] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_FoChopper->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _FoChopper->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_FoChopper->_position_relative, _FoChopper->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_FermiChopper_trace(_FoChopper, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component FoChopper [31] */
    if (!ABSORBED && _particle->_index == 32) {
      /* component PSD_Fermi2=PSD_monitor() [32] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_PSD_Fermi2->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _PSD_Fermi2->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_PSD_Fermi2->_position_relative, _PSD_Fermi2->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_PSD_monitor_trace(_PSD_Fermi2, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component PSD_Fermi2 [32] */
    if (!ABSORBED && _particle->_index == 33) {
      /* component DivMonfermi2=Divergence_monitor() [33] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_DivMonfermi2->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _DivMonfermi2->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_DivMonfermi2->_position_relative, _DivMonfermi2->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Divergence_monitor_trace(_DivMonfermi2, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component DivMonfermi2 [33] */
    if (!ABSORBED && _particle->_index == 34) {
      /* component FERMITOF1=Monitor_nD() [34] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_FERMITOF1->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _FERMITOF1->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_FERMITOF1->_position_relative, _FERMITOF1->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Monitor_nD_trace(_FERMITOF1, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component FERMITOF1 [34] */
    if (!ABSORBED && _particle->_index == 35) {
      /* component SAMPLE_SLIT=Slit() [35] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_SAMPLE_SLIT->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _SAMPLE_SLIT->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_SAMPLE_SLIT->_position_relative, _SAMPLE_SLIT->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Slit_trace(_SAMPLE_SLIT, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component SAMPLE_SLIT [35] */
    if (!ABSORBED && _particle->_index == 36) {
      /* component FERMITOF2=Monitor_nD() [36] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_FERMITOF2->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _FERMITOF2->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_FERMITOF2->_position_relative, _FERMITOF2->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Monitor_nD_trace(_FERMITOF2, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component FERMITOF2 [36] */
    if (!ABSORBED && _particle->_index == 37) {
      /* component PSD_SAMPLE=PSD_monitor() [37] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_PSD_SAMPLE->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _PSD_SAMPLE->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_PSD_SAMPLE->_position_relative, _PSD_SAMPLE->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_PSD_monitor_trace(_PSD_SAMPLE, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component PSD_SAMPLE [37] */
    if (!ABSORBED && _particle->_index == 38) {
      /* component DivMon_Sample=Divergence_monitor() [38] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_DivMon_Sample->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _DivMon_Sample->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_DivMon_Sample->_position_relative, _DivMon_Sample->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Divergence_monitor_trace(_DivMon_Sample, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component DivMon_Sample [38] */
    if (!ABSORBED && _particle->_index == 39) {
      /* component EMON_SAMPLE=E_monitor() [39] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_EMON_SAMPLE->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _EMON_SAMPLE->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_EMON_SAMPLE->_position_relative, _EMON_SAMPLE->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_E_monitor_trace(_EMON_SAMPLE, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component EMON_SAMPLE [39] */
    if (!ABSORBED && _particle->_index == 40) {
      /* component a3=Arm() [40] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_a3->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _a3->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_a3->_position_relative, _a3->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      _particle->_index++;
    } /* end component a3 [40] */
    /* start SPLIT at Sample */
    instrument->logic.Split_Sample_particle=_particle;
    if (!ABSORBED && _particle->_index == 41) 
    for (instrument->logic.Split_Sample = 0; instrument->logic.Split_Sample< 10; instrument->logic.Split_Sample++) {
    {
      /* component Sample=V_sample() [41] */
      _particle=instrument->logic.Split_Sample_particle;
      p /= 10 > 0 ? 10 : 1;
  #define radius (_Sample->_parameters.radius)
  #define thickness (_Sample->_parameters.thickness)
  #define zdepth (_Sample->_parameters.zdepth)
  #define Vc (_Sample->_parameters.Vc)
  #define sigma_abs (_Sample->_parameters.sigma_abs)
  #define sigma_inc (_Sample->_parameters.sigma_inc)
  #define radius_i (_Sample->_parameters.radius_i)
  #define radius_o (_Sample->_parameters.radius_o)
  #define h (_Sample->_parameters.h)
  #define focus_r (_Sample->_parameters.focus_r)
  #define pack (_Sample->_parameters.pack)
  #define frac (_Sample->_parameters.frac)
  #define f_QE (_Sample->_parameters.f_QE)
  #define gamma (_Sample->_parameters.gamma)
  #define target_x (_Sample->_parameters.target_x)
  #define target_y (_Sample->_parameters.target_y)
  #define target_z (_Sample->_parameters.target_z)
  #define focus_xw (_Sample->_parameters.focus_xw)
  #define focus_yh (_Sample->_parameters.focus_yh)
  #define focus_aw (_Sample->_parameters.focus_aw)
  #define focus_ah (_Sample->_parameters.focus_ah)
  #define xwidth (_Sample->_parameters.xwidth)
  #define yheight (_Sample->_parameters.yheight)
  #define zthick (_Sample->_parameters.zthick)
  #define rad_sphere (_Sample->_parameters.rad_sphere)
  #define sig_a (_Sample->_parameters.sig_a)
  #define sig_i (_Sample->_parameters.sig_i)
  #define V0 (_Sample->_parameters.V0)
  #define target_index (_Sample->_parameters.target_index)
  #define multiples (_Sample->_parameters.multiples)
  #define VarsV (_Sample->_parameters.VarsV)
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_Sample->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _Sample->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_Sample->_position_relative, _Sample->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_V_sample_trace(_Sample, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
  #undef radius
  #undef thickness
  #undef zdepth
  #undef Vc
  #undef sigma_abs
  #undef sigma_inc
  #undef radius_i
  #undef radius_o
  #undef h
  #undef focus_r
  #undef pack
  #undef frac
  #undef f_QE
  #undef gamma
  #undef target_x
  #undef target_y
  #undef target_z
  #undef focus_xw
  #undef focus_yh
  #undef focus_aw
  #undef focus_ah
  #undef xwidth
  #undef yheight
  #undef zthick
  #undef rad_sphere
  #undef sig_a
  #undef sig_i
  #undef V0
  #undef target_index
  #undef multiples
  #undef VarsV
      _particle->_index++;
    } /* end component Sample [41] */
    if (!ABSORBED && _particle->_index == 42) {
      /* component TOF_Det=Monitor_nD() [42] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_TOF_Det->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _TOF_Det->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_TOF_Det->_position_relative, _TOF_Det->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Monitor_nD_trace(_TOF_Det, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component TOF_Det [42] */
    if (!ABSORBED && _particle->_index == 43) {
      /* component FoDet=Monitor_nD() [43] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_FoDet->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _FoDet->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_FoDet->_position_relative, _FoDet->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_Monitor_nD_trace(_FoDet, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component FoDet [43] */
    if (!ABSORBED && _particle->_index == 44) {
      /* component EMON_DET=E_monitor() [44] */
      if (!flag_nocoordschange) { // flag activated by JUMP to pass coords change
        if (_EMON_DET->_rotation_is_identity) {
          coords_get(coords_add(coords_set(x,y,z), _EMON_DET->_position_relative),&x, &y, &z);
        } else
          mccoordschange(_EMON_DET->_position_relative, _EMON_DET->_rotation_relative, _particle);
      } else flag_nocoordschange=0;
      _particle_save = *_particle;
      class_E_monitor_trace(_EMON_DET, _particle);
      if (_particle->_restore)
        *_particle = _particle_save;
      _particle->_index++;
    } /* end component EMON_DET [44] */
    } /* end SPLIT at Sample */
    } /* end SPLIT at mono */
    if (_particle->_index > 44)
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
* instrument 'PSI_Focus' and components SAVE
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

_class_FermiChopper *class_FermiChopper_save(_class_FermiChopper *_comp
) {
  #define phase (_comp->_parameters.phase)
  #define radius (_comp->_parameters.radius)
  #define nu (_comp->_parameters.nu)
  #define w (_comp->_parameters.w)
  #define nslit (_comp->_parameters.nslit)
  #define R0 (_comp->_parameters.R0)
  #define Qc (_comp->_parameters.Qc)
  #define alpha (_comp->_parameters.alpha)
  #define m (_comp->_parameters.m)
  #define W (_comp->_parameters.W)
  #define length (_comp->_parameters.length)
  #define eff (_comp->_parameters.eff)
  #define zero_time (_comp->_parameters.zero_time)
  #define xwidth (_comp->_parameters.xwidth)
  #define verbose (_comp->_parameters.verbose)
  #define yheight (_comp->_parameters.yheight)
  #define curvature (_comp->_parameters.curvature)
  #define delay (_comp->_parameters.delay)
  #define FCVars (_comp->_parameters.FCVars)
  double mean_k, mean_v, mean_t, mean_w=0, mean_L=0;

  if (FCVars.sum_N) {
    mean_v = FCVars.sum_v/FCVars.sum_N;
    mean_t = FCVars.sum_t/FCVars.sum_N;
    mean_k = V2K*mean_v;
    if (mean_k) mean_L = 2*PI/mean_k;
    mean_w = VS2E*mean_v*mean_v;
    /* base opening time */
    double div, mean_phase;
    div        = atan2(w,length)/PI*180;
    mean_phase = fmod(mean_t*nu*360, 360);
    mean_phase+=radius/mean_v*360*nu;
    if (mean_phase > 180) mean_phase -= 360;

    if (!FCVars.sum_N_pass)
      printf("FermiChopper: %s: No neutron can pass the chopper.\n", NAME_CURRENT_COMP);
    if (!FCVars.sum_N_pass || verbose)
      printf("FermiChopper: %s\n"
           "              Mean velocity     v     = %g [m/s]\n"
           "              Mean wavelength   lambda= %g [Angs]\n"
           "              Mean energy       omega = %g [meV]\n"
           "              Mean arrival time t     = %g [s]\n"
           "              Mean phase              = %g [deg] (%s)\n"
           "              Slit pack divergence    = %g [deg] (full width)\n"
           "              Opening time      dt    = %g [s]\n"
           "              Intensity reaching FC   = %g [n/s]\n"
           "              Intensity passing  FC   = %g [n/s]\n"
           , NAME_CURRENT_COMP,
           mean_v, mean_L, mean_w, mean_t, mean_phase,
           (zero_time > 1 ? "set automatically" : "use phase=-(this value) to optimize"),
           2*div,
           (nu ? fabs(div/PI/nu) : 1),
           FCVars.sum_N,
           FCVars.sum_N_pass);
    if (!FCVars.sum_N_pass || verbose) {
      printf("FermiChopper: %s: Lost events anaylsis\n"
             "              Already inside:            %li\n"
             "              By Top/Bottom of cylinder: %li\n"
             "              At cylinder entrance:      %li\n"
             "              Hit cyl. entrance sides:   %li\n"
             "              Can't prop. to slit pack:  %li\n"
             "              At slit pack entrance:     %li\n"
             "              On absorbing slit coating: %li\n"
             "              Exiting slit pack:         %li\n"
             "              Too many iterations:       %li\n"
             "              Prop. in wrong direction : %li\n"
             "              Mad neutron (no control):  %li\n"
             "              At cylinder exit:          %li\n"
      , NAME_CURRENT_COMP,
      FCVars.absorb_alreadyinside,
      FCVars.absorb_topbottom,
      FCVars.absorb_cylentrance,
      FCVars.absorb_sideentrance,
      FCVars.absorb_notreachentrance,
      FCVars.absorb_packentrance,
      FCVars.absorb_slitcoating,

      FCVars.absorb_exitslitpack,
      FCVars.absorb_maxiterations,
      FCVars.absorb_wrongdirection,
      FCVars.absorb_nocontrol,
      FCVars.absorb_cylexit);

      if (FCVars.warn_notreachslitwall || FCVars.warn_notreachslitoutput)
        printf("Warning:      Can not reach slit wall:   %li\n"
               "Warning:      Can not reach slit output: %li\n",
        FCVars.warn_notreachslitwall,
        FCVars.warn_notreachslitoutput);
    }

  } else {
    printf("FermiChopper: %s: No neutron can reach the chopper.\n", NAME_CURRENT_COMP);
  }
  #undef phase
  #undef radius
  #undef nu
  #undef w
  #undef nslit
  #undef R0
  #undef Qc
  #undef alpha
  #undef m
  #undef W
  #undef length
  #undef eff
  #undef zero_time
  #undef xwidth
  #undef verbose
  #undef yheight
  #undef curvature
  #undef delay
  #undef FCVars
  return(_comp);
} /* class_FermiChopper_save */



int save(FILE *handle) { /* called by mccode_main for PSI_Focus:SAVE */
  if (!handle) siminfo_init(NULL);

  /* call iteratively all components SAVE */
  class_Progress_bar_save(_a1);






  class_L_monitor_save(_lambdaGuideExit);

  class_Divergence_monitor_save(_DivMonGuideExit);

  class_PSD_monitor_save(_PSDGuideExit);



  class_Monitor_nD_save(_DISCTOF);

  class_PSD_monitor_save(_PSDmon1Chopper);







  class_PSD_monitor_save(_PSDmonMono);

  class_Monitor_nD_save(_MONOTOF);

  class_Divergence_monitor_save(_DivMonMono);




  class_Monitor_nD_save(_FERMITOF_before);

  class_L_monitor_save(_lambdaFermi);

  class_E_monitor_save(_EMON_Fermi);

  class_Divergence_monitor_save(_DivMonfermi1);

  class_PSD_monitor_save(_PSD_Fermi1);

  class_FermiChopper_save(_FoChopper);

  class_PSD_monitor_save(_PSD_Fermi2);

  class_Divergence_monitor_save(_DivMonfermi2);

  class_Monitor_nD_save(_FERMITOF1);


  class_Monitor_nD_save(_FERMITOF2);

  class_PSD_monitor_save(_PSD_SAMPLE);

  class_Divergence_monitor_save(_DivMon_Sample);

  class_E_monitor_save(_EMON_SAMPLE);



  class_Monitor_nD_save(_TOF_Det);

  class_Monitor_nD_save(_FoDet);

  class_E_monitor_save(_EMON_DET);

  if (!handle) siminfo_close(); 

  return(0);
} /* save */

/* *****************************************************************************
* instrument 'PSI_Focus' and components FINALLY
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



int finally(void) { /* called by mccode_main for PSI_Focus:FINALLY */
  siminfo_init(NULL);
  save(siminfo_file); /* save data when simulation ends */

  /* call iteratively all components FINALLY */
  class_Progress_bar_finally(_a1);

  class_Source_gen_finally(_csource);





  class_L_monitor_finally(_lambdaGuideExit);

  class_Divergence_monitor_finally(_DivMonGuideExit);

  class_PSD_monitor_finally(_PSDGuideExit);



  class_Monitor_nD_finally(_DISCTOF);

  class_PSD_monitor_finally(_PSDmon1Chopper);







  class_PSD_monitor_finally(_PSDmonMono);

  class_Monitor_nD_finally(_MONOTOF);

  class_Divergence_monitor_finally(_DivMonMono);




  class_Monitor_nD_finally(_FERMITOF_before);

  class_L_monitor_finally(_lambdaFermi);

  class_E_monitor_finally(_EMON_Fermi);

  class_Divergence_monitor_finally(_DivMonfermi1);

  class_PSD_monitor_finally(_PSD_Fermi1);


  class_PSD_monitor_finally(_PSD_Fermi2);

  class_Divergence_monitor_finally(_DivMonfermi2);

  class_Monitor_nD_finally(_FERMITOF1);


  class_Monitor_nD_finally(_FERMITOF2);

  class_PSD_monitor_finally(_PSD_SAMPLE);

  class_Divergence_monitor_finally(_DivMon_Sample);

  class_E_monitor_finally(_EMON_SAMPLE);



  class_Monitor_nD_finally(_TOF_Det);

  class_Monitor_nD_finally(_FoDet);

  class_E_monitor_finally(_EMON_DET);

  siminfo_close(); 

  return(0);
} /* finally */

/* *****************************************************************************
* instrument 'PSI_Focus' and components DISPLAY
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

_class_Bender *class_Bender_display(_class_Bender *_comp
) {
  #define w (_comp->_parameters.w)
  #define h (_comp->_parameters.h)
  #define r (_comp->_parameters.r)
  #define Win (_comp->_parameters.Win)
  #define k (_comp->_parameters.k)
  #define d (_comp->_parameters.d)
  #define l (_comp->_parameters.l)
  #define R0a (_comp->_parameters.R0a)
  #define Qca (_comp->_parameters.Qca)
  #define alphaa (_comp->_parameters.alphaa)
  #define ma (_comp->_parameters.ma)
  #define Wa (_comp->_parameters.Wa)
  #define R0i (_comp->_parameters.R0i)
  #define Qci (_comp->_parameters.Qci)
  #define alphai (_comp->_parameters.alphai)
  #define mi (_comp->_parameters.mi)
  #define Wi (_comp->_parameters.Wi)
  #define R0s (_comp->_parameters.R0s)
  #define Qcs (_comp->_parameters.Qcs)
  #define alphas (_comp->_parameters.alphas)
  #define ms (_comp->_parameters.ms)
  #define Ws (_comp->_parameters.Ws)
  #define bk (_comp->_parameters.bk)
  #define mWin (_comp->_parameters.mWin)
  int i;
  double w1c, w2c, h1, h2, L, w1, w2;

  w1c = (w + d)/(double)k;
  w2c = w1c; h1 = h; h2 = h;
  L = r*mWin; w1 = w; w2 = w;

  
  for(i = 0; i < k; i++)
  {
    multiline(5,
              i*w1c - w1/2.0, -h1/2.0, 0.0,
              i*w2c - w2/2.0, -h2/2.0, (double)L,
              i*w2c - w2/2.0,  h2/2.0, (double)L,
              i*w1c - w1/2.0,  h1/2.0, 0.0,
              i*w1c - w1/2.0, -h1/2.0, 0.0);
    multiline(5,
              (i+1)*w1c - d - w1/2.0, -h1/2.0, 0.0,
              (i+1)*w2c - d - w2/2.0, -h2/2.0, (double)L,
              (i+1)*w2c - d - w2/2.0,  h2/2.0, (double)L,
              (i+1)*w1c - d - w1/2.0,  h1/2.0, 0.0,
              (i+1)*w1c - d - w1/2.0, -h1/2.0, 0.0);
  }
  line(-w1/2.0, -h1/2.0, 0.0, w1/2.0, -h1/2.0, 0.0);
  line(-w2/2.0, -h2/2.0, (double)L, w2/2.0, -h2/2.0, (double)L);
  #undef w
  #undef h
  #undef r
  #undef Win
  #undef k
  #undef d
  #undef l
  #undef R0a
  #undef Qca
  #undef alphaa
  #undef ma
  #undef Wa
  #undef R0i
  #undef Qci
  #undef alphai
  #undef mi
  #undef Wi
  #undef R0s
  #undef Qcs
  #undef alphas
  #undef ms
  #undef Ws
  #undef bk
  #undef mWin
  return(_comp);
} /* class_Bender_display */

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

_class_Guide_channeled *class_Guide_channeled_display(_class_Guide_channeled *_comp
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
  #define nslit (_comp->_parameters.nslit)
  #define d (_comp->_parameters.d)
  #define Qcx (_comp->_parameters.Qcx)
  #define Qcy (_comp->_parameters.Qcy)
  #define alphax (_comp->_parameters.alphax)
  #define alphay (_comp->_parameters.alphay)
  #define W (_comp->_parameters.W)
  #define mx (_comp->_parameters.mx)
  #define my (_comp->_parameters.my)
  #define nu (_comp->_parameters.nu)
  #define phase (_comp->_parameters.phase)
  #define w1c (_comp->_parameters.w1c)
  #define w2c (_comp->_parameters.w2c)
  #define ww (_comp->_parameters.ww)
  #define hh (_comp->_parameters.hh)
  #define whalf (_comp->_parameters.whalf)
  #define hhalf (_comp->_parameters.hhalf)
  #define lwhalf (_comp->_parameters.lwhalf)
  #define lhhalf (_comp->_parameters.lhhalf)
  int i;

  
  for(i = 0; i < nslit; i++)
  {
    multiline(5,
              i*w1c - w1/2.0, -h1/2.0, 0.0,
              i*w2c - w2/2.0, -h2/2.0, (double)l,
              i*w2c - w2/2.0,  h2/2.0, (double)l,
              i*w1c - w1/2.0,  h1/2.0, 0.0,
              i*w1c - w1/2.0, -h1/2.0, 0.0);
    multiline(5,
              (i+1)*w1c - d - w1/2.0, -h1/2.0, 0.0,
              (i+1)*w2c - d - w2/2.0, -h2/2.0, (double)l,
              (i+1)*w2c - d - w2/2.0,  h2/2.0, (double)l,
              (i+1)*w1c - d - w1/2.0,  h1/2.0, 0.0,
              (i+1)*w1c - d - w1/2.0, -h1/2.0, 0.0);
  }
  line(-w1/2.0, -h1/2.0, 0.0, w1/2.0, -h1/2.0, 0.0);
  line(-w2/2.0, -h2/2.0, (double)l, w2/2.0, -h2/2.0, (double)l);

  if (nu || phase) {
    double radius = sqrt(w1*w1+l*l);
    /* cylinder top/center/bottom  */
    circle("xz", 0,-h1/2,l/2,radius);
    circle("xz", 0,0    ,l/2,radius);
    circle("xz", 0, h1/2,l/2,radius);
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
  #undef nslit
  #undef d
  #undef Qcx
  #undef Qcy
  #undef alphax
  #undef alphay
  #undef W
  #undef mx
  #undef my
  #undef nu
  #undef phase
  #undef w1c
  #undef w2c
  #undef ww
  #undef hh
  #undef whalf
  #undef hhalf
  #undef lwhalf
  #undef lhhalf
  return(_comp);
} /* class_Guide_channeled_display */

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

_class_Arm *class_Arm_display(_class_Arm *_comp
) {
  /* A bit ugly; hard-coded dimensions. */
  
  line(0,0,0,0.2,0,0);
  line(0,0,0,0,0.2,0);
  line(0,0,0,0,0,0.2);
  return(_comp);
} /* class_Arm_display */

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

_class_FermiChopper *class_FermiChopper_display(_class_FermiChopper *_comp
) {
  #define phase (_comp->_parameters.phase)
  #define radius (_comp->_parameters.radius)
  #define nu (_comp->_parameters.nu)
  #define w (_comp->_parameters.w)
  #define nslit (_comp->_parameters.nslit)
  #define R0 (_comp->_parameters.R0)
  #define Qc (_comp->_parameters.Qc)
  #define alpha (_comp->_parameters.alpha)
  #define m (_comp->_parameters.m)
  #define W (_comp->_parameters.W)
  #define length (_comp->_parameters.length)
  #define eff (_comp->_parameters.eff)
  #define zero_time (_comp->_parameters.zero_time)
  #define xwidth (_comp->_parameters.xwidth)
  #define verbose (_comp->_parameters.verbose)
  #define yheight (_comp->_parameters.yheight)
  #define curvature (_comp->_parameters.curvature)
  #define delay (_comp->_parameters.delay)
  #define FCVars (_comp->_parameters.FCVars)
  double index_x=0;
  double index_z=0;
  double xpos, zpos;
  double Nz,Nx;
  double ymin=-yheight/2.0;
  double ymax= yheight/2.0;

  double omega=FCVars.omega;
  FCVars.omega=0;
  // FCVars.ph0  =0;
  Nz = (FCVars.C_slit ?  4 : 1);
  Nx = (nslit > 11    ? 11 : nslit);
  FCVars.C_slit *= -1;
  
  /* cylinder top/center/bottom  */
  circle("xz", 0,ymax,0,radius);
  circle("xz", 0,0   ,0,radius);
  circle("xz", 0,ymin,0,radius);
  /* vertical lines to make a kind of volume */
  line( 0  ,ymin,-radius, 0  ,ymax,-radius);
  line( 0  ,ymin, radius, 0  ,ymax, radius);
  line(-radius,ymin, 0  ,-radius,ymax, 0  );
  line( radius,ymin, 0  , radius,ymax, 0  );
  /* slit package */
  for (index_x = -Nx/2; index_x < Nx/2; index_x++) {
    for (index_z = -Nz/2; index_z < Nz/2; index_z++) {
      double xs1, zs1, zs2;
      double xp1, xp2, zp1, zp2;
      zs1 = length*index_z/Nz;
      zs2 = length*(index_z+1)/Nz;
      xs1 = w*nslit*index_x/Nx;
      xp1 = FC_xrot(xs1, zs1, 0, FCVars);
      xp2 = FC_xrot(xs1, zs2, 0, FCVars);
      zp1 = FC_zrot(xs1, zs1, 0, FCVars);
      zp2 = FC_zrot(xs1, zs2, 0, FCVars);
      multiline(5, xp1, ymin, zp1,
                   xp1, ymax, zp1,
                   xp2, ymax, zp2,
                   xp2, ymin, zp2,
                   xp1, ymin, zp1);
    }
  }
  /* cylinder inner sides containing slit package */
  double xp1, xp2, zp1, zp2;
  xpos = nslit*w/2;
  zpos = sqrt(radius*radius - xpos*xpos);
  xp1 = FC_xrot(xpos, -zpos, 0, FCVars);
  xp2 = FC_xrot(xpos, +zpos, 0, FCVars);
  zp1 = FC_zrot(xpos, -zpos, 0, FCVars);
  zp2 = FC_zrot(xpos, +zpos, 0, FCVars);
  multiline(5,  xp1, ymin, zp1,
                xp1, ymax, zp1,
                xp2, ymax, zp2,
                xp2, ymin, zp2,
                xp1, ymin, zp1);
  xpos *= -1;
  xp1 = FC_xrot(xpos, -zpos, 0, FCVars);
  xp2 = FC_xrot(xpos, +zpos, 0, FCVars);
  zp1 = FC_zrot(xpos, -zpos, 0, FCVars);
  zp2 = FC_zrot(xpos, +zpos, 0, FCVars);
  multiline(5,  xp1, ymin, zp1,
                xp1, ymax, zp1,
                xp2, ymax, zp2,
                xp2, ymin, zp2,
                xp1, ymin, zp1);
  FCVars.omega=omega;
  #undef phase
  #undef radius
  #undef nu
  #undef w
  #undef nslit
  #undef R0
  #undef Qc
  #undef alpha
  #undef m
  #undef W
  #undef length
  #undef eff
  #undef zero_time
  #undef xwidth
  #undef verbose
  #undef yheight
  #undef curvature
  #undef delay
  #undef FCVars
  return(_comp);
} /* class_FermiChopper_display */

_class_V_sample *class_V_sample_display(_class_V_sample *_comp
) {
  #define radius (_comp->_parameters.radius)
  #define thickness (_comp->_parameters.thickness)
  #define zdepth (_comp->_parameters.zdepth)
  #define Vc (_comp->_parameters.Vc)
  #define sigma_abs (_comp->_parameters.sigma_abs)
  #define sigma_inc (_comp->_parameters.sigma_inc)
  #define radius_i (_comp->_parameters.radius_i)
  #define radius_o (_comp->_parameters.radius_o)
  #define h (_comp->_parameters.h)
  #define focus_r (_comp->_parameters.focus_r)
  #define pack (_comp->_parameters.pack)
  #define frac (_comp->_parameters.frac)
  #define f_QE (_comp->_parameters.f_QE)
  #define gamma (_comp->_parameters.gamma)
  #define target_x (_comp->_parameters.target_x)
  #define target_y (_comp->_parameters.target_y)
  #define target_z (_comp->_parameters.target_z)
  #define focus_xw (_comp->_parameters.focus_xw)
  #define focus_yh (_comp->_parameters.focus_yh)
  #define focus_aw (_comp->_parameters.focus_aw)
  #define focus_ah (_comp->_parameters.focus_ah)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define zthick (_comp->_parameters.zthick)
  #define rad_sphere (_comp->_parameters.rad_sphere)
  #define sig_a (_comp->_parameters.sig_a)
  #define sig_i (_comp->_parameters.sig_i)
  #define V0 (_comp->_parameters.V0)
  #define target_index (_comp->_parameters.target_index)
  #define multiples (_comp->_parameters.multiples)
  #define VarsV (_comp->_parameters.VarsV)
  
  if (VarsV.shapetyp == 0) {
    circle("xz", 0,  h/2.0, 0, radius_i);
    circle("xz", 0,  h/2.0, 0, radius_o);
    circle("xz", 0, -h/2.0, 0, radius_i);
    circle("xz", 0, -h/2.0, 0, radius_o);
    line(-radius_i, -h/2.0, 0, -radius_i, +h/2.0, 0);
    line(+radius_i, -h/2.0, 0, +radius_i, +h/2.0, 0);
    line(0, -h/2.0, -radius_i, 0, +h/2.0, -radius_i);
    line(0, -h/2.0, +radius_i, 0, +h/2.0, +radius_i);
    line(-radius_o, -h/2.0, 0, -radius_o, +h/2.0, 0);
    line(+radius_o, -h/2.0, 0, +radius_o, +h/2.0, 0);
    line(0, -h/2.0, -radius_o, 0, +h/2.0, -radius_o);
    line(0, -h/2.0, +radius_o, 0, +h/2.0, +radius_o);
  }
  else { 
	if (VarsV.shapetyp == 1) {
      double xmin = -0.5*xwidth;
      double xmax =  0.5*xwidth;
      double ymin = -0.5*yheight;
      double ymax =  0.5*yheight;
      double zmin = -0.5*zthick;
      double zmax =  0.5*zthick;
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
    else {
      circle("xy", 0,  0.0, 0, rad_sphere);
      circle("xz", 0,  0.0, 0, rad_sphere);
      circle("yz", 0,  0.0, 0, rad_sphere);        
    }
  }
  #undef radius
  #undef thickness
  #undef zdepth
  #undef Vc
  #undef sigma_abs
  #undef sigma_inc
  #undef radius_i
  #undef radius_o
  #undef h
  #undef focus_r
  #undef pack
  #undef frac
  #undef f_QE
  #undef gamma
  #undef target_x
  #undef target_y
  #undef target_z
  #undef focus_xw
  #undef focus_yh
  #undef focus_aw
  #undef focus_ah
  #undef xwidth
  #undef yheight
  #undef zthick
  #undef rad_sphere
  #undef sig_a
  #undef sig_i
  #undef V0
  #undef target_index
  #undef multiples
  #undef VarsV
  return(_comp);
} /* class_V_sample_display */


  #undef magnify
  #undef line
  #undef dashed_line
  #undef multiline
  #undef rectangle
  #undef box
  #undef circle
  #undef cylinder
  #undef sphere

int display(void) { /* called by mccode_main for PSI_Focus:DISPLAY */
  printf("MCDISPLAY: start\n");

  /* call iteratively all components DISPLAY */
  class_Progress_bar_display(_a1);

  class_Source_gen_display(_csource);

  class_Guide_display(_guide1);

  class_Bender_display(_guide2);

  class_Guide_display(_bunker);

  class_Guide_display(_guide3);

  class_L_monitor_display(_lambdaGuideExit);

  class_Divergence_monitor_display(_DivMonGuideExit);

  class_PSD_monitor_display(_PSDGuideExit);

  class_Guide_channeled_display(_FOCUSguide);

  class_DiskChopper_display(_FirstChopper);

  class_Monitor_nD_display(_DISCTOF);

  class_PSD_monitor_display(_PSDmon1Chopper);

  class_Slit_display(_VacuumTube1entry);

  class_Slit_display(_VacuumTube1exit);

  class_Slit_display(_VacuumTube2entry);

  class_Slit_display(_VacuumTube2exit);

  class_Slit_display(_VacuumTube3entry);

  class_Slit_display(_VacuumTube3exit);

  class_PSD_monitor_display(_PSDmonMono);

  class_Monitor_nD_display(_MONOTOF);

  class_Divergence_monitor_display(_DivMonMono);

  class_Arm_display(_focus_mono);

  class_Monochromator_2foc_display(_mono);

  class_Arm_display(_a2);

  class_Monitor_nD_display(_FERMITOF_before);

  class_L_monitor_display(_lambdaFermi);

  class_E_monitor_display(_EMON_Fermi);

  class_Divergence_monitor_display(_DivMonfermi1);

  class_PSD_monitor_display(_PSD_Fermi1);

  class_FermiChopper_display(_FoChopper);

  class_PSD_monitor_display(_PSD_Fermi2);

  class_Divergence_monitor_display(_DivMonfermi2);

  class_Monitor_nD_display(_FERMITOF1);

  class_Slit_display(_SAMPLE_SLIT);

  class_Monitor_nD_display(_FERMITOF2);

  class_PSD_monitor_display(_PSD_SAMPLE);

  class_Divergence_monitor_display(_DivMon_Sample);

  class_E_monitor_display(_EMON_SAMPLE);

  class_Arm_display(_a3);

  class_V_sample_display(_Sample);

  class_Monitor_nD_display(_TOF_Det);

  class_Monitor_nD_display(_FoDet);

  class_E_monitor_display(_EMON_DET);

  printf("MCDISPLAY: end\n");

  return(0);
} /* display */

void* _getvar_parameters(char* compname)
/* enables settings parameters based use of the GETPAR macro */
{
  if (!strcmp(compname, "a1")) return (void *) &(_a1->_parameters);
  if (!strcmp(compname, "csource")) return (void *) &(_csource->_parameters);
  if (!strcmp(compname, "guide1")) return (void *) &(_guide1->_parameters);
  if (!strcmp(compname, "guide2")) return (void *) &(_guide2->_parameters);
  if (!strcmp(compname, "bunker")) return (void *) &(_bunker->_parameters);
  if (!strcmp(compname, "guide3")) return (void *) &(_guide3->_parameters);
  if (!strcmp(compname, "lambdaGuideExit")) return (void *) &(_lambdaGuideExit->_parameters);
  if (!strcmp(compname, "DivMonGuideExit")) return (void *) &(_DivMonGuideExit->_parameters);
  if (!strcmp(compname, "PSDGuideExit")) return (void *) &(_PSDGuideExit->_parameters);
  if (!strcmp(compname, "FOCUSguide")) return (void *) &(_FOCUSguide->_parameters);
  if (!strcmp(compname, "FirstChopper")) return (void *) &(_FirstChopper->_parameters);
  if (!strcmp(compname, "DISCTOF")) return (void *) &(_DISCTOF->_parameters);
  if (!strcmp(compname, "PSDmon1Chopper")) return (void *) &(_PSDmon1Chopper->_parameters);
  if (!strcmp(compname, "VacuumTube1entry")) return (void *) &(_VacuumTube1entry->_parameters);
  if (!strcmp(compname, "VacuumTube1exit")) return (void *) &(_VacuumTube1exit->_parameters);
  if (!strcmp(compname, "VacuumTube2entry")) return (void *) &(_VacuumTube2entry->_parameters);
  if (!strcmp(compname, "VacuumTube2exit")) return (void *) &(_VacuumTube2exit->_parameters);
  if (!strcmp(compname, "VacuumTube3entry")) return (void *) &(_VacuumTube3entry->_parameters);
  if (!strcmp(compname, "VacuumTube3exit")) return (void *) &(_VacuumTube3exit->_parameters);
  if (!strcmp(compname, "PSDmonMono")) return (void *) &(_PSDmonMono->_parameters);
  if (!strcmp(compname, "MONOTOF")) return (void *) &(_MONOTOF->_parameters);
  if (!strcmp(compname, "DivMonMono")) return (void *) &(_DivMonMono->_parameters);
  if (!strcmp(compname, "focus_mono")) return (void *) &(_focus_mono->_parameters);
  if (!strcmp(compname, "mono")) return (void *) &(_mono->_parameters);
  if (!strcmp(compname, "a2")) return (void *) &(_a2->_parameters);
  if (!strcmp(compname, "FERMITOF_before")) return (void *) &(_FERMITOF_before->_parameters);
  if (!strcmp(compname, "lambdaFermi")) return (void *) &(_lambdaFermi->_parameters);
  if (!strcmp(compname, "EMON_Fermi")) return (void *) &(_EMON_Fermi->_parameters);
  if (!strcmp(compname, "DivMonfermi1")) return (void *) &(_DivMonfermi1->_parameters);
  if (!strcmp(compname, "PSD_Fermi1")) return (void *) &(_PSD_Fermi1->_parameters);
  if (!strcmp(compname, "FoChopper")) return (void *) &(_FoChopper->_parameters);
  if (!strcmp(compname, "PSD_Fermi2")) return (void *) &(_PSD_Fermi2->_parameters);
  if (!strcmp(compname, "DivMonfermi2")) return (void *) &(_DivMonfermi2->_parameters);
  if (!strcmp(compname, "FERMITOF1")) return (void *) &(_FERMITOF1->_parameters);
  if (!strcmp(compname, "SAMPLE_SLIT")) return (void *) &(_SAMPLE_SLIT->_parameters);
  if (!strcmp(compname, "FERMITOF2")) return (void *) &(_FERMITOF2->_parameters);
  if (!strcmp(compname, "PSD_SAMPLE")) return (void *) &(_PSD_SAMPLE->_parameters);
  if (!strcmp(compname, "DivMon_Sample")) return (void *) &(_DivMon_Sample->_parameters);
  if (!strcmp(compname, "EMON_SAMPLE")) return (void *) &(_EMON_SAMPLE->_parameters);
  if (!strcmp(compname, "a3")) return (void *) &(_a3->_parameters);
  if (!strcmp(compname, "Sample")) return (void *) &(_Sample->_parameters);
  if (!strcmp(compname, "TOF_Det")) return (void *) &(_TOF_Det->_parameters);
  if (!strcmp(compname, "FoDet")) return (void *) &(_FoDet->_parameters);
  if (!strcmp(compname, "EMON_DET")) return (void *) &(_EMON_DET->_parameters);
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

/* end of generated C code PSI_Focus.c */
