/* Automatically generated file. Do not edit. 
 * Format:     ANSI C source code
 * Creator:    McStas <http://www.mcstas.org>
 * Instrument: /zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr (Test_Monochromators)
 * Date:       Wed Feb 26 19:20:59 2020
 * File:       ./Test_Monochromators.c
 * Compile:    cc -o Test_Monochromators.out ./Test_Monochromators.c 
 * CFLAGS=
 */


#define MCCODE_STRING "McStas 2.5 - Feb. 26, 2020"
#define FLAVOR "mcstas"
#define FLAVOR_UPPER "MCSTAS"
#define MC_USE_DEFAULT_MAIN
#define MC_TRACE_ENABLED
#define MC_EMBEDDED_RUNTIME

#line 1 "mccode-r.h"
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
* Release: McStas 2.5
* Version: $Revision$
*
* Runtime system header for McStas/McXtrace.
*
* In order to use this library as an external library, the following variables
* and macros must be declared (see details in the code)
*
*   struct mcinputtable_struct mcinputtable[];
*   int mcnumipar;
*   char mcinstrument_name[], mcinstrument_source[];
*   int mctraceenabled, mcdefaultmain;
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

#include <math.h>
#include <string.h>
#include <strings.h>
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
#define mcstatic static
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

#ifndef WIN32
#ifndef _POSIX_SOURCE
#define _POSIX_SOURCE 1
#endif
#endif

/* the version string is replaced when building distribution with mkdist */
#ifndef MCCODE_STRING
#define MCCODE_STRING "McStas 2.5 - Feb. 26, 2020"
#endif

#ifndef MCCODE_DATE
#define MCCODE_DATE "Feb. 26, 2020"
#endif

#ifndef MCCODE_VERSION
#define MCCODE_VERSION "2.5"
#endif

#ifndef MCCODE_NAME
#define MCCODE_NAME "McStas"
#endif

#ifndef MCCODE_PARTICLE
#define MCCODE_PARTICLE "neutron"
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
    instr_type_double, instr_type_int, instr_type_string
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
extern struct mcinputtable_struct mcinputtable[]; /* list of instrument parameters */
extern int    mcnumipar;                          /* number of instrument parameters */
extern char   mcinstrument_name[], mcinstrument_source[]; /* instrument name and filename */
extern char  *mcinstrument_exe;                           /* executable path = argv[0] or NULL */
extern MCNUM  mccomp_storein[]; /* 11 coords * number of components in instrument */
extern MCNUM  mcAbsorbProp[];
extern MCNUM  mcScattered;      /* number of SCATTER calls in current component */
extern MCNUM  mcRestore;        /* Flag to indicate if neutron needs to be restored */
#ifndef MC_ANCIENT_COMPATIBILITY
extern int mctraceenabled, mcdefaultmain;
#endif
#endif


/* Useful macros ============================================================ */

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
#define SIG_MESSAGE(msg) strcpy(mcsig_message, msg);
#else
#define SIG_MESSAGE(msg)
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
#  define M_2_SQRTPI 2/sqrt(M_PI)
#  define M_SQRT2 sqrt(2)
#  define M_SQRT1_2 sqrt(1/2)
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
#define POS_A_COMP_INDEX(index) \
    (mccomp_posa[index])
#define POS_R_COMP_INDEX(index) \
    (mccomp_posr[index])
/* number of SCATTER calls in current comp: mcScattered defined in generated C code */
#define SCATTERED mcScattered
/* Flag to indicate if neutron needs to be restored: mcRestore defined in generated C code */
#define RESTORE mcRestore


/* Retrieve component information from the kernel */
/* Name, position and orientation (both absolute and relative)  */
/* Any component: For "redundancy", see comment by KN */
#define tmp_name_comp(comp) #comp
#define NAME_COMP(comp) tmp_name_comp(comp)
#define tmp_pos_a_comp(comp) (mcposa ## comp)
#define POS_A_COMP(comp) tmp_pos_a_comp(comp)
#define tmp_pos_r_comp(comp) (mcposr ## comp)
#define POS_R_COMP(comp) tmp_pos_r_comp(comp)
#define tmp_rot_a_comp(comp) (mcrota ## comp)
#define ROT_A_COMP(comp) tmp_rot_a_comp(comp)
#define tmp_rot_r_comp(comp) (mcrotr ## comp)
#define ROT_R_COMP(comp) tmp_rot_r_comp(comp)

/* Current component name, index, position and orientation */
#define NAME_CURRENT_COMP  NAME_COMP(mccompcurname)
#define INDEX_CURRENT_COMP mccompcurindex
#define POS_A_CURRENT_COMP POS_A_COMP(mccompcurname)
#define POS_R_CURRENT_COMP POS_R_COMP(mccompcurname)
#define ROT_A_CURRENT_COMP ROT_A_COMP(mccompcurname)
#define ROT_R_CURRENT_COMP ROT_R_COMP(mccompcurname)

/* Note: The two-stage approach to MC_GETPAR is NOT redundant; without it,
* after #define C sample, MC_GETPAR(C,x) would refer to component C, not to
* component sample. Such are the joys of ANSI C.

* Anyway the usage of MCGETPAR requires that we use sometimes bare names...
*/
#define MC_GETPAR2(comp, par) (mcc ## comp ## _ ## par)
#define MC_GETPAR(comp, par) MC_GETPAR2(comp,par)

/* MCDISPLAY/trace and debugging message sent to stdout */
#ifdef MC_TRACE_ENABLED
#define DEBUG
#endif

#ifdef DEBUG
#define mcDEBUG_INSTR() if(!mcdotrace); else { printf("\nINSTRUMENT:\n"); printf("Instrument '%s' (%s)\n", mcinstrument_name, mcinstrument_source); }
#define mcDEBUG_COMPONENT(name,c,t) if(!mcdotrace); else {\
  printf("COMPONENT: \"%s\"\n" \
         "POS: %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g\n", \
         name, c.x, c.y, c.z, t[0][0], t[0][1], t[0][2], \
         t[1][0], t[1][1], t[1][2], t[2][0], t[2][1], t[2][2]); \
  mcAccumulatedILength += coords_len(coords_sub(mcLastComp,c)); \
  printf("Component %30s AT (%g,%g,%g)    %g m from origin\n", name, c.x, c.y, c.z, mcAccumulatedILength); \
  mcLastComp=c;\
  }
#define mcDEBUG_INSTR_END() if(!mcdotrace); else printf("INSTRUMENT END:\n");
#define mcDEBUG_ENTER() if(!mcdotrace); else printf("ENTER:\n");
#define mcDEBUG_COMP(c) if(!mcdotrace); else printf("COMP: \"%s\"\n", c);
#define mcDEBUG_LEAVE() if(!mcdotrace); else printf("LEAVE:\n");
#define mcDEBUG_ABSORB() if(!mcdotrace); else printf("ABSORB:\n");
#else
#define mcDEBUG_INSTR()
#define mcDEBUG_COMPONENT(name,c,t)
#define mcDEBUG_INSTR_END()
#define mcDEBUG_ENTER()
#define mcDEBUG_COMP(c)
#define mcDEBUG_LEAVE()
#define mcDEBUG_ABSORB()
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
void mcdis_dashed_linemcdis_dashed_line(double x1, double y1, double z1,
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

/* selection of random number generator. default is MT */
#ifndef MC_RAND_ALG
#define MC_RAND_ALG 1
#endif

#if MC_RAND_ALG == 0
   /* Use system random() (not recommended). */
#  define MC_RAND_MAX RAND_MAX
#elif MC_RAND_ALG == 1
   /* "Mersenne Twister", by Makoto Matsumoto and Takuji Nishimura. */
#  define MC_RAND_MAX ((unsigned long)0xffffffff)
#  define random mt_random
#  define srandom mt_srandom
#elif MC_RAND_ALG == 2
   /* Algorithm used in McStas CVS-080208 and earlier (not recommended). */
#  define MC_RAND_MAX 0x7fffffff
#  define random mc_random
#  define srandom mc_srandom
#else
#  error "Bad value for random number generator choice."
#endif

typedef int mc_int32_t;
mc_int32_t mc_random(void);
void mc_srandom (unsigned int x);
unsigned long mt_random(void);
void mt_srandom (unsigned long x);

double rand01();
double randpm1();
double rand0max(double max);
double randminmax(double min, double max);

double randnorm(void);
double randtriangle(void);

#ifndef DANSE
void mcinit(void);
void mcraytrace(void);
void mcsave(FILE *);
void mcfinally(void);
void mcdisplay(void);
#endif

/* simple vector algebra ==================================================== */
#define vec_prod(x, y, z, x1, y1, z1, x2, y2, z2) \
	vec_prod_func(&x, &y, &z, x1, y1, z1, x2, y2, z2)
mcstatic void vec_prod_func(double *x, double *y, double *z,
		double x1, double y1, double z1, double x2, double y2, double z2);

mcstatic double scalar_prod(
		double x1, double y1, double z1, double x2, double y2, double z2);

#define NORM(x,y,z) \
	norm_func(&x, &y, &z)
mcstatic void norm_func(double *x, double *y, double *z) {
	double temp = (*x * *x) + (*y * *y) + (*z * *z);
	if (temp != 0) {
		temp = sqrt(temp);
		*x /= temp;
		*y /= temp;
		*z /= temp;
	}
}
#define normal_vec(nx, ny, nz, x, y, z) \
    normal_vec_func(&(nx), &(ny), &(nz), x, y, z)
mcstatic void normal_vec_func(double *nx, double *ny, double *nz,
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

void mccoordschange(Coords a, Rotation t, double *x, double *y, double *z,
    double *vx, double *vy, double *vz, double *sx, double *sy, double *sz);
void
mccoordschange_polarisation(Rotation t, double *sx, double *sy, double *sz);

double mcestimate_error(double N, double p1, double p2);
void mcreadparams(void);

/* this is now in mcstas-r.h and mcxtrace-r.h as the number of state parameters is no longer equal*/
/* void mcsetstate(double x, double y, double z, double vx, double vy, double vz,
                double t, double sx, double sy, double sz, double p);
*/
void mcgenstate(void);

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
void randvec_target_circle(double *xo, double *yo, double *zo,
    double *solid_angle, double xi, double yi, double zi, double radius);
#define randvec_target_sphere randvec_target_circle
void randvec_target_rect_angular(double *xo, double *yo, double *zo,
    double *solid_angle,
               double xi, double yi, double zi, double height, double width, Rotation A);
#define randvec_target_rect(p0,p1,p2,p3,p4,p5,p6,p7,p8,p9)  randvec_target_rect_real(p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,0,0,0,1)
void randvec_target_rect_real(double *xo, double *yo, double *zo,
    double *solid_angle,
	       double xi, double yi, double zi, double height, double width, Rotation A,
			 double lx, double ly, double lz, int order);

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

static   char *mcdirname             = NULL;      /* name of output directory */
static   char *mcsiminfo_name        = "mccode";  /* default output sim file name */
char    *mcformat                    = NULL;      /* NULL (default) or a specific format */

/* file I/O definitions and function prototypes */

#ifndef MC_EMBEDDED_RUNTIME /* the mcstatic variables (from mccode-r.c) */
extern FILE * mcsiminfo_file;     /* handle to the output siminfo file */
extern int    mcgravitation;      /* flag to enable gravitation */
extern int    mcdotrace;          /* flag to print MCDISPLAY messages */
#else
mcstatic FILE *mcsiminfo_file        = NULL;
#endif

/* I/O function prototypes ================================================== */

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

#line 712 "./Test_Monochromators.c"

#line 1 "mcstas-r.h"
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
*   char mcinstrument_name[], mcinstrument_source[];
*   int mctraceenabled, mcdefaultmain;
*   extern MCNUM  mccomp_storein[];
*   extern MCNUM  mcAbsorbProp[];
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

#define SCATTER do {mcDEBUG_SCATTER(mcnlx, mcnly, mcnlz, mcnlvx, mcnlvy, mcnlvz, \
        mcnlt,mcnlsx,mcnlsy,mcnlsz, mcnlp); mcScattered++;} while(0)
#define ABSORB do {mcDEBUG_STATE(mcnlx, mcnly, mcnlz, mcnlvx, mcnlvy, mcnlvz, \
        mcnlt,mcnlsx,mcnlsy,mcnlsz, mcnlp); mcDEBUG_ABSORB(); MAGNET_OFF; goto mcabsorb;} while(0)

#define STORE_NEUTRON(index, x, y, z, vx, vy, vz, t, sx, sy, sz, p) \
  mcstore_neutron(mccomp_storein,index, x, y, z, vx, vy, vz, t, sx, sy, sz, p);
#define RESTORE_NEUTRON(index, x, y, z, vx, vy, vz, t, sx, sy, sz, p) \
  mcrestore_neutron(mccomp_storein,index, &x, &y, &z, &vx, &vy, &vz, &t, &sx, &sy, &sz, &p);

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
  }while (0)
    /* change coordinates from local system to magnet system */
/*    Rotation rotLM, rotTemp; \
      Coords   posLM = coords_sub(POS_A_CURRENT_COMP, mcMagnetPos); \
      rot_transpose(ROT_A_CURRENT_COMP, rotTemp); \
      rot_mul(rotTemp, mcMagnetRot, rotLM); \
      mcMagnetPrecession(mcnlx, mcnly, mcnlz, mcnlt, mcnlvx, mcnlvy, mcnlvz, \
               &mcnlsx, &mcnlsy, &mcnlsz, dt, posLM, rotLM); \
      } while(0)
*/

#define mcPROP_DT(dt) \
  do { \
    if (mcMagnet && dt > 0) PROP_MAGNET(dt);\
    mcnlx += mcnlvx*(dt); \
    mcnly += mcnlvy*(dt); \
    mcnlz += mcnlvz*(dt); \
    mcnlt += (dt); \
    if (isnan(p) || isinf(p)) { mcAbsorbProp[INDEX_CURRENT_COMP]++; ABSORB; }\
  } while(0)

/* ADD: E. Farhi, Aug 6th, 2001 PROP_GRAV_DT propagation with acceleration */
#define PROP_GRAV_DT(dt, Ax, Ay, Az) \
  do { \
    if(dt < 0 && mcallowbackprop == 0) { mcAbsorbProp[INDEX_CURRENT_COMP]++; ABSORB; }\
    if (mcMagnet) printf("Spin precession gravity\n"); \
    mcnlx  += mcnlvx*(dt) + (Ax)*(dt)*(dt)/2; \
    mcnly  += mcnlvy*(dt) + (Ay)*(dt)*(dt)/2; \
    mcnlz  += mcnlvz*(dt) + (Az)*(dt)*(dt)/2; \
    mcnlvx += (Ax)*(dt); \
    mcnlvy += (Ay)*(dt); \
    mcnlvz += (Az)*(dt); \
    mcnlt  += (dt); \
    DISALLOW_BACKPROP;\
  } while(0)


#define PROP_DT(dt) \
  do { \
    if(dt < 0) { RESTORE=1; goto mcabsorbComp; }; \
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
    mc_ret = solve_2nd_order(&mc_dt, NULL, -mc_gz/2, -mcnlvz, -mcnlz); \
    if (mc_ret && mc_dt>=0) {PROP_GRAV_DT(mc_dt, mc_gx, mc_gy, mc_gz); mcnlz=0;}\
    else { if (mcallowbackprop ==0) {mcAbsorbProp[INDEX_CURRENT_COMP]++; ABSORB; }}; }\
    else mcPROP_Z0; \
    DISALLOW_BACKPROP;\
  } while(0)

#define mcPROP_Z0 \
  do { \
    double mc_dt; \
    if(mcnlvz == 0) { mcAbsorbProp[INDEX_CURRENT_COMP]++; ABSORB; }; \
    mc_dt = -mcnlz/mcnlvz; \
    if(mc_dt < 0 && mcallowbackprop == 0) { mcAbsorbProp[INDEX_CURRENT_COMP]++; ABSORB; }; \
    mcPROP_DT(mc_dt); \
    mcnlz = 0; \
    DISALLOW_BACKPROP;\
  } while(0)

#define PROP_X0 \
  do { \
    if (mcgravitation) { Coords mcLocG; int mc_ret; \
    double mc_dt, mc_gx, mc_gy, mc_gz; \
    mcLocG = rot_apply(ROT_A_CURRENT_COMP, coords_set(0,-GRAVITY,0)); \
    coords_get(mcLocG, &mc_gx, &mc_gy, &mc_gz); \
    mc_ret = solve_2nd_order(&mc_dt, NULL, -mc_gx/2, -mcnlvx, -mcnlx); \
    if (mc_ret && mc_dt>=0) PROP_GRAV_DT(mc_dt, mc_gx, mc_gy, mc_gz); \
    else { if (mcallowbackprop ==0) {mcAbsorbProp[INDEX_CURRENT_COMP]++; ABSORB; }}; }\
    else mcPROP_X0; \
    DISALLOW_BACKPROP;\
  } while(0)

#define mcPROP_X0 \
  do { \
    double mc_dt; \
    if(mcnlvx == 0) { mcAbsorbProp[INDEX_CURRENT_COMP]++; ABSORB; }; \
    mc_dt = -mcnlx/mcnlvx; \
    if(mc_dt < 0 && mcallowbackprop == 0) { mcAbsorbProp[INDEX_CURRENT_COMP]++; ABSORB; }; \
    mcPROP_DT(mc_dt); \
    mcnlx = 0; \
    DISALLOW_BACKPROP;\
  } while(0)

#define PROP_Y0 \
  do { \
    if (mcgravitation) { Coords mcLocG; int mc_ret; \
    double mc_dt, mc_gx, mc_gy, mc_gz; \
    mcLocG = rot_apply(ROT_A_CURRENT_COMP, coords_set(0,-GRAVITY,0)); \
    coords_get(mcLocG, &mc_gx, &mc_gy, &mc_gz); \
    mc_ret = solve_2nd_order(&mc_dt, NULL, -mc_gy/2, -mcnlvy, -mcnly); \
    if (mc_ret && mc_dt>=0) PROP_GRAV_DT(mc_dt, mc_gx, mc_gy, mc_gz); \
    else { if (mcallowbackprop ==0) {mcAbsorbProp[INDEX_CURRENT_COMP]++; ABSORB; }}; }\
    else mcPROP_Y0; \
    DISALLOW_BACKPROP;\
  } while(0)


#define mcPROP_Y0 \
  do { \
    double mc_dt; \
    if(mcnlvy == 0) { mcAbsorbProp[INDEX_CURRENT_COMP]++; ABSORB; }; \
    mc_dt = -mcnly/mcnlvy; \
    if(mc_dt < 0 && mcallowbackprop == 0) { mcAbsorbProp[INDEX_CURRENT_COMP]++; ABSORB; }; \
    mcPROP_DT(mc_dt); \
    mcnly = 0; \
    DISALLOW_BACKPROP; \
  } while(0)

/*moved from mccode-r.h*/
void mcsetstate(double x, double y, double z, double vx, double vy, double vz,
                double t, double sx, double sy, double sz, double p);

#ifdef DEBUG

#define mcDEBUG_STATE(x,y,z,vx,vy,vz,t,sx,sy,sz,p) if(!mcdotrace); else \
  printf("STATE: %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g\n", \
         x,y,z,vx,vy,vz,t,sx,sy,sz,p);
#define mcDEBUG_SCATTER(x,y,z,vx,vy,vz,t,sx,sy,sz,p) if(!mcdotrace); else \
  printf("SCATTER: %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g\n", \
         x,y,z,vx,vy,vz,t,sx,sy,sz,p);

#else

#define mcDEBUG_STATE(x,y,z,vx,vy,vz,t,sx,sy,sz,p)
#define mcDEBUG_SCATTER(x,y,z,vx,vy,vz,t,sx,sy,sz,p)

#endif

#endif /* !MCCODE_H */

#endif /* MCSTAS_R_H */
/* End of file "mcstas-r.h". */

#line 945 "./Test_Monochromators.c"

#line 1 "mccode-r.c"
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
#endif

#include <sys/stat.h>

#ifdef _WIN32 
#include <direct.h>
# define  mkdir( D, M )   _mkdir( D ) 
#endif 

#ifndef DANSE
#ifdef MC_ANCIENT_COMPATIBILITY
int mctraceenabled = 0;
int mcdefaultmain  = 0;
#endif
/* else defined directly in the McCode generated C code */

static   long mcseed                 = 0; /* seed for random generator */
static   long mcstartdate            = 0; /* start simulation time */
static   int  mcdisable_output_files = 0; /* --no-output-files */
mcstatic int  mcgravitation          = 0; /* use gravitation flag, for PROP macros */
int      mcMagnet                    = 0; /* magnet stack flag */
mcstatic int  mcdotrace              = 0; /* flag for --trace and messages for DISPLAY */
int      mcallowbackprop             = 0;         /* flag to enable negative/backprop */

/* Number of particle histories to simulate. */
#ifdef NEUTRONICS
mcstatic unsigned long long int mcncount             = 1;
mcstatic unsigned long long int mcrun_num            = 0;
#else
mcstatic unsigned long long int mcncount             = 1000000;
mcstatic unsigned long long int mcrun_num            = 0;
#endif /* NEUTRONICS */

#else
#include "mcstas-globals.h"
#endif /* !DANSE */

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
      if (MPI_Reduce((double*)(sbuf+offset), (double*)(rbuf+offset),
              length, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD) != MPI_SUCCESS)
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
    mcparm_double, mcparminfo_double, mcparmerror_double,
    mcparmprinter_double
  }, {
    mcparm_int, mcparminfo_int, mcparmerror_int,
    mcparmprinter_int
  }, {
    mcparm_string, mcparminfo_string, mcparmerror_string,
    mcparmprinter_string
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

/*******************************************************************************
* mcfull_file: allocates a full file name=mcdirname+file. Catenate extension if missing.
*******************************************************************************/
char *mcfull_file(char *name, char *ext)
{
  int   dirlen=0;
  char *mem   =NULL;

  dirlen = mcdirname ? strlen(mcdirname) : 0;
  mem = (char*)malloc(dirlen + strlen(name) + CHAR_BUF_LENGTH);
  if(!mem) {
    exit(-fprintf(stderr, "Error: Out of memory %li (mcfull_file)\n", (long)(dirlen + strlen(name) + 256)));
  }
  strcpy(mem, "");

  /* prepend directory name to path if name does not contain a path */
  if (dirlen > 0 && !strchr(name, MC_PATHSEP_C)) {
    strcat(mem, mcdirname);
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
* mcnew_file: opens a new file within mcdirname if non NULL
*             the file is opened in "a" (append, create if does not exist)
*             the extension 'ext' is added if the file name does not include one.
*             the last argument is set to 0 if file did not exist, else to 1.
*******************************************************************************/
FILE *mcnew_file(char *name, char *ext, int *exists)
{
  char *mem;
  FILE *file=NULL;

  if (!name || strlen(name) == 0 || mcdisable_output_files) return(NULL);
  
  mem  = mcfull_file(name, ext); /* create mcdirname/name.ext */
  
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
* Used by: mcdetector_import
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
      exit(-fprintf(stderr, "Error: Out of memory creating %li 1D " MCCODE_STRING " data set for file '%s' (mcdetector_import)\n",
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
* mcdetector_import: build detector structure, merge non-lists from MPI
*                    compute basic stat, write "Detector:" line
* RETURN:            detector structure. Invalid data if detector.p1 == NULL
*                    Invalid detector sets m=0 and filename=""
*                    Simulation data  sets m=0 and filename=mcsiminfo_name
* This function is equivalent to the old 'mcdetector_out', returning a structure
*******************************************************************************/
MCDETECTOR mcdetector_import(
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

  snprintf(detector.instrument, CHAR_BUF_LENGTH, "%s (%s)", mcinstrument_name, mcinstrument_source);
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

  m=labs(m); n=labs(n); p=labs(p); /* make sure dimensions are positive */
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
    if (!strcmp(detector.component, mcinstrument_name)) {
      if (strlen(detector.filename))  /* we name it from its filename, or from its title */
        strncpy(c, detector.filename, CHAR_BUF_LENGTH);
      else
        snprintf(c, CHAR_BUF_LENGTH, "%s", mcinstrument_name);
    } else
      strncpy(c, detector.component, CHAR_BUF_LENGTH);  /* usual detectors written by components */

    printf("Detector: %s_I=%g %s_ERR=%g %s_N=%g",
           c, detector.intensity,
           c, detector.error,
           c, detector.events);
    printf(" \"%s\"\n", strlen(detector.filename) ? detector.filename : detector.component);
  }
  

  return(detector);
} /* mcdetector_import */

/* end MCDETECTOR import section ============================================ */

















/* ========================================================================== */

/*                               ASCII output                                 */
/*     The SIM file is YAML based, the data files have '#' headers            */

/* ========================================================================== */


/*******************************************************************************
* mcinfo_out: output instrument tags/info (only in SIM)
* Used in: mcsiminfo_init (ascii), mcinfo(stdout)
*******************************************************************************/
static void mcinfo_out(char *pre, FILE *f)
{
  char Parameters[CHAR_BUF_LENGTH] = "";
  int  i;

  if (!f || mcdisable_output_files) return;

  /* create parameter string ================================================ */
  for(i = 0; i < mcnumipar; i++)
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
    fprintf(f, "%sFile: %s%c%s\n",    pre, mcdirname, MC_PATHSEP_C, mcsiminfo_name);
  else
    fprintf(f, "%sCreator: %s\n",     pre, MCCODE_STRING);

  fprintf(f, "%sSource: %s\n",   pre, mcinstrument_source);
  fprintf(f, "%sParameters: %s\n",    pre, Parameters);
  
  fprintf(f, "%sTrace_enabled: %s\n", pre, mctraceenabled ? "yes" : "no");
  fprintf(f, "%sDefault_main: %s\n",  pre, mcdefaultmain ?  "yes" : "no");
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
* mcruninfo_out_backend: output simulation tags/info (both in SIM and data files)
* Used in: mcsiminfo_init (ascii case), mcdetector_out_xD_ascii, mcinfo(stdout)
*******************************************************************************/
static void mcruninfo_out_backend(char *pre, FILE *f, int info)
{
  int i;
  char Parameters[CHAR_BUF_LENGTH];

  if (!f || mcdisable_output_files) return;

  fprintf(f, "%sFormat: %s%s\n",      pre, 
    mcformat && strlen(mcformat) ? mcformat : MCCODE_NAME,
    mcformat && strcasestr(mcformat,"McCode") ? " with text headers" : "");
  fprintf(f, "%sURL: %s\n",         pre, "http://www.mccode.org");
  fprintf(f, "%sCreator: %s\n",     pre, MCCODE_STRING);
  fprintf(f, "%sInstrument: %s\n", pre, mcinstrument_source);
  fprintf(f, "%sNcount: %llu\n",        pre, mcget_ncount());
  fprintf(f, "%sTrace: %s\n",       pre, mcdotrace ? "yes" : "no");
  fprintf(f, "%sGravitation: %s\n", pre, mcgravitation ? "yes" : "no");
  snprintf(Parameters, CHAR_BUF_LENGTH, "%ld", mcseed);
  fprintf(f, "%sSeed: %s\n",        pre, Parameters);
  fprintf(f, "%sDirectory: %s\n",        pre, mcdirname ? mcdirname : ".");
#ifdef USE_MPI
  if (mpi_node_count > 1)
    fprintf(f, "%sNodes: %i\n",        pre, mpi_node_count);
#endif

  /* output parameter string ================================================ */
  for(i = 0; i < mcnumipar; i++) {
      if (!info){
          (*mcinputtypes[mcinputtable[i].type].printer)(Parameters, mcinputtable[i].par);
          fprintf(f, "%sParam: %s=%s\n", pre, mcinputtable[i].name, Parameters);
      }else{
        /*if an info run, some variables might not have values. Flag these by "NULL"*/
	if(mcinputtable[i].val && strlen(mcinputtable[i].val)){
            /* ... those with defautl values*/
            (*mcinputtypes[mcinputtable[i].type].printer)(Parameters, mcinputtable[i].par);
            fprintf(f, "%sParam: %s=%s\n", pre, mcinputtable[i].name, Parameters);
        }else{
            /* ... and those without */
            fprintf(f, "%sParam: %s=NULL\n", pre, mcinputtable[i].name);
	}
      }
  }
} /* mcruninfo_out_backend */

/************************
* wrapper function to mcruninfo_out_backend
*  Regular runs use this whereas the single call from mcinfo is directly to the backend
*************************/
static void mcruninfo_out(char *pre, FILE *f){
    mcruninfo_out_backend(pre,f,0);
}

/*******************************************************************************
* mcsiminfo_out:    wrapper to fprintf(mcsiminfo_file)
*******************************************************************************/
void mcsiminfo_out(char *format, ...)
{
  va_list ap;

  if(mcsiminfo_file && !mcdisable_output_files)
  {
    va_start(ap, format);
    vfprintf(mcsiminfo_file, format, ap);
    va_end(ap);
  }
} /* mcsiminfo_out */


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
    mcsiminfo_out("\nbegin data\n"); // detector.component
    mcdatainfo_out("  ", mcsiminfo_file, detector);
    mcsiminfo_out("end data\n");
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
    mcsiminfo_out("\nbegin data\n"); // detector.filename
    mcdatainfo_out("  ", mcsiminfo_file, detector);
    mcsiminfo_out("end data\n");
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
        mcsiminfo_out("\nbegin data\n"); // detector.filename
        mcdatainfo_out("  ", mcsiminfo_file, detector);
        mcsiminfo_out("end data\n");
      
        mcruninfo_out( "# ", outfile);
        mcdatainfo_out("# ", outfile,   detector);
      }
      fprintf(outfile, "# Data [%s/%s] %s:\n", detector.component, detector.filename, detector.zvar);
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
* Used in: mcsiminfo_init (nexus)
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
  nxprintattr(f, "creator",   "%s generated with " MCCODE_STRING, mcinstrument_name);
  
  /* count the number of existing NXentry and create the next one */
  NXgetgroupinfo(f, &count, name, class);
  sprintf(entry0, "entry%i", count+1);

  /* create the main NXentry (mandatory in NeXus) */
  if (NXmakegroup(f, entry0, "NXentry") == NX_OK) 
  if (NXopengroup(f, entry0, "NXentry") == NX_OK) {
    
    nxprintf(nxhandle, "program_name", MCCODE_STRING);
    nxprintf(f, "start_time", ctime(&t));
    nxprintf(f, "title", "%s%s%s simulation generated by instrument %s", 
      mcdirname && strlen(mcdirname) ? mcdirname : ".", MC_PATHSEP_S, mcsiminfo_name,
      mcinstrument_name);
    nxprintattr(f, "program_name", MCCODE_STRING);
    nxprintattr(f, "instrument",   mcinstrument_name);
    nxprintattr(f, "simulation",   "%s%s%s",
        mcdirname && strlen(mcdirname) ? mcdirname : ".", MC_PATHSEP_S, mcsiminfo_name);

    /* write NeXus instrument group */
    if (NXmakegroup(f, "instrument", "NXinstrument") == NX_OK)
    if (NXopengroup(f, "instrument", "NXinstrument") == NX_OK) {
      int   i;
      char *string=NULL;

      /* write NeXus parameters(types) data =================================== */
      string = (char*)malloc(CHAR_BUF_LENGTH);
      if (string) {
        strcpy(string, "");
        for(i = 0; i < mcnumipar; i++)
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
        
      nxprintattr(f, "name",          mcinstrument_name);
      nxprintf   (f, "name",          mcinstrument_name);
      nxprintattr(f, "Source",        mcinstrument_source);
      
      nxprintattr(f, "Trace_enabled", mctraceenabled ? "yes" : "no");
      nxprintattr(f, "Default_main",  mcdefaultmain ?  "yes" : "no");
      nxprintattr(f, "Embedded_runtime",  
  #ifdef MC_EMBEDDED_RUNTIME
           "yes"
  #else
           "no"
  #endif
           );
           
      /* add instrument source code when available */
      buffer = mcinfo_readfile(mcinstrument_source);
      if (buffer && strlen(buffer)) {
        long length=strlen(buffer);
        nxprintf (f, "description", buffer);
        NXopendata(f,"description");
        nxprintattr(f, "file_name", mcinstrument_source);
        nxprintattr(f, "file_size", "%li", length);
        nxprintattr(f, "MCCODE_STRING", MCCODE_STRING);
        NXclosedata(f);
        nxprintf (f,"instrument_source", "%s " MCCODE_NAME " " MCCODE_PARTICLE " Monte Carlo simulation", mcinstrument_name);
        free(buffer);
      } else
        nxprintf (f, "description", "File %s not found (instrument description %s is missing)", 
          mcinstrument_source, mcinstrument_name);
      
      /* add Mantid/IDF.xml when available */
      char *IDFfile=NULL;
      IDFfile = (char*)malloc(CHAR_BUF_LENGTH);
      sprintf(IDFfile,"%s%s",mcinstrument_source,".xml");
      buffer = mcinfo_readfile(IDFfile);
      if (buffer && strlen(buffer)) {
        NXmakegroup (nxhandle, "instrument_xml", "NXnote");
        NXopengroup (nxhandle, "instrument_xml", "NXnote");
        nxprintf(f, "data", buffer);
        nxprintf(f, "description", "IDF.xml file found with instrument %s", mcinstrument_source);
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
        mcdirname && strlen(mcdirname) ? mcdirname : ".", MC_PATHSEP_S, mcsiminfo_name);
      
      nxprintf   (f, "name",      "%s",     mcsiminfo_name);
      nxprintattr(f, "Format",    mcformat && strlen(mcformat) ? mcformat : MCCODE_NAME);
      nxprintattr(f, "URL",       "http://www.mccode.org");
      nxprintattr(f, "program",   MCCODE_STRING);
      nxprintattr(f, "Instrument",mcinstrument_source);
      nxprintattr(f, "Trace",     mcdotrace ?     "yes" : "no");
      nxprintattr(f, "Gravitation",mcgravitation ? "yes" : "no");
      nxprintattr(f, "Seed",      "%li", mcseed);
      nxprintattr(f, "Directory", mcdirname);
    #ifdef USE_MPI
      if (mpi_node_count > 1)
        nxprintf(f, "Nodes", "%i",        mpi_node_count);
    #endif
    
      /* output parameter string ================================================ */
      if (NXmakegroup(f, "Param", "NXparameters") == NX_OK)
      if (NXopengroup(f, "Param", "NXparameters") == NX_OK) {
        int i;
        char string[CHAR_BUF_LENGTH];
        for(i = 0; i < mcnumipar; i++) {
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

  /* the NXdetector group has been created in mcinfo_out_nexus (mcsiminfo_init) */
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

  /* the NXdetector group has been created in mcinfo_out_nexus (mcsiminfo_init) */
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

MCDETECTOR mcdetector_out_1D_nexus(MCDETECTOR detector_inc)
{
  MCDETECTOR detector = detector_inc;
  MPI_MASTER(
  mcdatainfo_out_nexus(nxhandle, detector);
  mcdetector_out_data_nexus(nxhandle, detector);
  );
  return(detector);
} /* mcdetector_out_1D_ascii */

MCDETECTOR mcdetector_out_2D_nexus(MCDETECTOR detector_inc)
{
  MCDETECTOR detector = detector_inc;
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
* mcsiminfo_init:   open SIM and write header
*******************************************************************************/
FILE *mcsiminfo_init(FILE *f)
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
  if (mcsiminfo_file || mcdisable_output_files) 
    return (mcsiminfo_file);
    
#ifdef USE_NEXUS
  /* only master writes NeXus header: calls NXopen(nxhandle) */
  if (mcformat && strcasestr(mcformat, "NeXus")) {
	  MPI_MASTER(
	  mcsiminfo_file = mcnew_file(mcsiminfo_name, "h5", &exists);
    if(!mcsiminfo_file)
      fprintf(stderr,
	      "Warning: could not open simulation description file '%s'\n",
	      mcsiminfo_name);
	  else
	    mcinfo_out_nexus(nxhandle);
	  );
    return(mcsiminfo_file); /* points to nxhandle */
  }
#endif
  
  /* write main description file (only MASTER) */
  MPI_MASTER(

  mcsiminfo_file = mcnew_file(mcsiminfo_name, "sim", &exists);
  if(!mcsiminfo_file)
    fprintf(stderr,
	    "Warning: could not open simulation description file '%s'\n",
	    mcsiminfo_name);
  else
  {
    /* write SIM header */
    time_t t=time(NULL);
    mcsiminfo_out("%s simulation description file for %s.\n", 
      MCCODE_NAME, mcinstrument_name);
    mcsiminfo_out("Date:    %s", ctime(&t)); /* includes \n */
    mcsiminfo_out("Program: %s\n\n", MCCODE_STRING);
    
    mcsiminfo_out("begin instrument: %s\n", mcinstrument_name);
    mcinfo_out(   "  ", mcsiminfo_file);
    mcsiminfo_out("end instrument\n");

    mcsiminfo_out("\nbegin simulation: %s\n", mcdirname);
    mcruninfo_out("  ", mcsiminfo_file);
    mcsiminfo_out("end simulation\n");

  }
  return (mcsiminfo_file);
  
  ); /* MPI_MASTER */
  
} /* mcsiminfo_init */

/*******************************************************************************
*   mcsiminfo_close:  close SIM
*******************************************************************************/
void mcsiminfo_close()
{
  MPI_MASTER(
  if(mcsiminfo_file && !mcdisable_output_files) {
#ifdef USE_NEXUS
    if (mcformat && strcasestr(mcformat, "NeXus")) {
      time_t t=time(NULL);
      nxprintf(nxhandle, "end_time", ctime(&t));
      nxprintf(nxhandle, "duration", "%li", (long)t-mcstartdate);
      NXclosegroup(nxhandle); /* NXentry */
      NXclose(&nxhandle);
    } else
#endif
      fclose(mcsiminfo_file);
    );
    mcsiminfo_file = NULL;
  }
} /* mcsiminfo_close */

/*******************************************************************************
* mcdetector_out_0D: wrapper for 0D (single value).
*   Output single detector/monitor data (p0, p1, p2).
*   Title is t, component name is c.
*******************************************************************************/
MCDETECTOR mcdetector_out_0D(char *t, double p0, double p1, double p2,
                         char *c, Coords posa)
{
  /* import and perform basic detector analysis (and handle MPI reduce) */
  MCDETECTOR detector = mcdetector_import(mcformat,
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
  MCDETECTOR detector = mcdetector_import(mcformat,
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
  if (xl && strlen(xl)) { strncpy(xvar, xl, CHAR_BUF_LENGTH); xvar[strcspn(xvar,"\n\r ")]='\0'; }
  else strcpy(xvar, "x");
  if (yl && strlen(yl)) { strncpy(yvar, yl, CHAR_BUF_LENGTH); yvar[strcspn(yvar,"\n\r ")]='\0'; }
  else strcpy(yvar, "y");

  MCDETECTOR detector;

  /* import and perform basic detector analysis (and handle MPI_Reduce) */
  if (labs(m) == 1) {/* n>1 on Y, m==1 on X: 1D, no X axis*/
    detector = mcdetector_import(mcformat,
      c, (t ? t : MCCODE_STRING " 1D data"),
      n, 1, 1,
      yl, "", "Signal per bin",
      yvar, "(I,Ierr)", "I",
      y1, y2, x1, x2, 0, 0, f,
      p0, p1, p2, posa); /* write Detector: line */
  } else if (labs(n)==1) {/* m>1 on X, n==1 on Y: 1D, no Y axis*/
    detector = mcdetector_import(mcformat,
      c, (t ? t : MCCODE_STRING " 1D data"),
      m, 1, 1,
      xl, "", "Signal per bin",
      xvar, "(I,Ierr)", "I",
      x1, x2, y1, y2, 0, 0, f,
      p0, p1, p2, posa); /* write Detector: line */
  }else {
    detector = mcdetector_import(mcformat,
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
                  1,labs(m),1,labs(n),
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
    mcdirname = dir;
  else
    mcdirname = dir+strlen("file://");
  
  
  
  MPI_MASTER(
    if(mkdir(mcdirname, 0777)) {
#ifndef DANSE
      fprintf(stderr, "Error: unable to create directory '%s' (mcuse_dir)\n", dir);
      fprintf(stderr, "(Maybe the directory already exists?)\n");
#endif
#ifdef USE_MPI
    MPI_Abort(MPI_COMM_WORLD, -1);
#endif
    exit(-1);
    }
  ); /* MPI_MASTER */
  
  /* remove trailing PATHSEP (if any) */
  while (strlen(mcdirname) && mcdirname[strlen(mcdirname) - 1] == MC_PATHSEP_C)
    mcdirname[strlen(mcdirname) - 1]='\0';
#endif /* !MC_PORTABLE */
} /* mcuse_dir */

/*******************************************************************************
* mcinfo: display instrument simulation info to stdout and exit
*******************************************************************************/
static void
mcinfo(void)
{
  fprintf(stdout, "begin instrument: %s\n", mcinstrument_name);
  mcinfo_out("  ", stdout);
  fprintf(stdout, "end instrument\n");
  fprintf(stdout, "begin simulation: %s\n", mcdirname ? mcdirname : ".");
  mcruninfo_out_backend("  ", stdout,1);
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
unsigned long long int mcget_ncount(void)
{
  return mcncount;
}

/* mcget_run_num: get curent number of rays in TRACE */
unsigned long long int mcget_run_num(void)
{
  return mcrun_num;
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
            mcdis_line(x+r*sin(i*2*M_PI/24),y,z+r*cos(i*2*M_PI/24),
                    x+r*sin((i+1)*2*M_PI/24),y,z+r*cos((i+1)*2*M_PI/24));
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
            rotate(ux,uy,uz, mx,my,mz, i*2*M_PI/24, nx,ny,nz);
            rotate(wx,wy,wz, mx,my,mz, (i+1)*2*M_PI/24, nx,ny,nz);
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
        rotate(ux,uy,uz, mx,my,mz, i*2*M_PI/24, nx,ny,nz);
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
        rotate(nx,ny,nz, nx,ny,nz, M_PI/N, 0,1,0);
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
Coords
coords_set(MCNUM x, MCNUM y, MCNUM z)
{
  Coords a;

  a.x = x;
  a.y = y;
  a.z = z;
  return a;
}

/* coords_get: get coordinates. Required when 'x','y','z' are #defined as ray pars */
Coords
coords_get(Coords a, MCNUM *x, MCNUM *y, MCNUM *z)
{
  *x = a.x;
  *y = a.y;
  *z = a.z;
  return a;
}

/* coords_add: Add two coordinates. */
Coords
coords_add(Coords a, Coords b)
{
  Coords c;

  c.x = a.x + b.x;
  c.y = a.y + b.y;
  c.z = a.z + b.z;
  if (fabs(c.z) < 1e-14) c.z=0.0;
  return c;
}

/* coords_sub: Subtract two coordinates. */
Coords
coords_sub(Coords a, Coords b)
{
  Coords c;

  c.x = a.x - b.x;
  c.y = a.y - b.y;
  c.z = a.z - b.z;
  if (fabs(c.z) < 1e-14) c.z=0.0;
  return c;
}

/* coords_neg: Negate coordinates. */
Coords
coords_neg(Coords a)
{
  Coords b;

  b.x = -a.x;
  b.y = -a.y;
  b.z = -a.z;
  return b;
}

/* coords_scale: Scale a vector. */
Coords coords_scale(Coords b, double scale) {
  Coords a;

  a.x = b.x*scale;
  a.y = b.y*scale;
  a.z = b.z*scale;
  return a;
}

/* coords_sp: Scalar product: a . b */
double coords_sp(Coords a, Coords b) {
  double value;

  value = a.x*b.x + a.y*b.y + a.z*b.z;
  return value;
}

/* coords_xp: Cross product: a = b x c. */
Coords coords_xp(Coords b, Coords c) {
  Coords a;

  a.x = b.y*c.z - c.y*b.z;
  a.y = b.z*c.x - c.z*b.x;
  a.z = b.x*c.y - c.x*b.y;
  return a;
}

/* coords_len: Gives length of coords set. */
double coords_len(Coords a) {
  return sqrt(a.x*a.x + a.y*a.y + a.z*a.z);
}

/* coords_mirror: Mirror a in plane (through the origin) defined by normal n*/
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
void coords_print(Coords a) {

  fprintf(stdout, "(%f, %f, %f)\n", a.x, a.y, a.z);
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
void
rot_set_rotation(Rotation t, double phx, double phy, double phz)
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
int
rot_test_identity(Rotation t)
{
  return (t[0][0] + t[1][1] + t[2][2] == 3);
}

/*******************************************************************************
* rot_mul: Matrix multiplication of transformations (this corresponds to
* combining transformations). After rot_mul(T1, T2, T3), doing T3 is
* equal to doing first T2, then T1.
* Note that T3 must not alias (use the same array as) T1 or T2.
*******************************************************************************/
void
rot_mul(Rotation t1, Rotation t2, Rotation t3)
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
void
rot_copy(Rotation dest, Rotation src)
{
  int i,j;
  for(i = 0; i < 3; i++)
    for(j = 0; j < 3; j++)
      dest[i][j] = src[i][j];
}

/*******************************************************************************
* rot_transpose: Matrix transposition, which is inversion for Rotation matrices
*******************************************************************************/
void
rot_transpose(Rotation src, Rotation dst)
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
Coords
rot_apply(Rotation t, Coords a)
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
mcstatic void vec_prod_func(double *x, double *y, double *z,
		double x1, double y1, double z1,
		double x2, double y2, double z2) {
    *x = (y1)*(z2) - (y2)*(z1);
    *y = (z1)*(x2) - (z2)*(x1);
    *z = (x1)*(y2) - (x2)*(y1);
}

/**
 * Scalar product: use coords_sp for Coords.
 */
mcstatic double scalar_prod(
		double x1, double y1, double z1,
		double x2, double y2, double z2) {
	return ((x1 * x2) + (y1 * y2) + (z1 * z2));
}

/*******************************************************************************
* mccoordschange: applies rotation to (x y z) and (vx vy vz) and Spin (sx,sy,sz)
*******************************************************************************/
void
mccoordschange(Coords a, Rotation t, double *x, double *y, double *z,
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

  if ( (vz && vy  && vx) && (*vz != 0.0 || *vx != 0.0 || *vy != 0.0) ) mccoordschange_polarisation(t, vx, vy, vz);

  if ( (sz && sy  && sx) && (*sz != 0.0 || *sx != 0.0 || *sy != 0.0) ) mccoordschange_polarisation(t, sx, sy, sz);

}

/*******************************************************************************
* mccoordschange_polarisation: applies rotation to vector (sx sy sz)
*******************************************************************************/
void
mccoordschange_polarisation(Rotation t, double *sx, double *sy, double *sz)
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
mcstatic void normal_vec_func(double *nx, double *ny, double *nz,
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
void
randvec_target_circle(double *xo, double *yo, double *zo, double *solid_angle,
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
} /* randvec_target_circle */

/*******************************************************************************
 * randvec_target_rect_angular: Choose random direction towards target at
 * (xi,yi,zi) with given ANGULAR dimension height x width. height=phi_x=[0,PI],
 * width=phi_y=[0,2*PI] (radians)
 * If height or width is zero, choose random direction in full 4PI, no target.
 *******************************************************************************/
void
randvec_target_rect_angular(double *xo, double *yo, double *zo, double *solid_angle,
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

} /* randvec_target_rect_angular */

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

void
randvec_target_rect_real(double *xo, double *yo, double *zo, double *solid_angle,
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
} /* randvec_target_rect_real */

/* SECTION: random numbers ================================================== */

/*
 * Copyright (c) 1983 Regents of the University of California.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms are permitted
 * provided that the above copyright notice and this paragraph are
 * duplicated in all such forms and that any documentation,
 * advertising materials, and other materials related to such
 * distribution and use acknowledge that the software was developed
 * by the University of California, Berkeley.  The name of the
 * University may not be used to endorse or promote products derived
 * from this software without specific prior written permission.
 * THIS SOFTWARE IS PROVIDED ``AS IS'' AND WITHOUT ANY EXPRESS OR
 * IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED
 * WARRANTIES OF MERCHANTIBILITY AND FITNESS FOR A PARTICULAR PURPOSE.
 */

/*
 * This is derived from the Berkeley source:
 *        @(#)random.c        5.5 (Berkeley) 7/6/88
 * It was reworked for the GNU C Library by Roland McGrath.
 * Rewritten to use reentrant functions by Ulrich Drepper, 1995.
 */

/*******************************************************************************
* Modified for McStas from glibc 2.0.7pre1 stdlib/random.c and
* stdlib/random_r.c.
*
* This way random() is more than four times faster compared to calling
* standard glibc random() on ix86 Linux, probably due to multithread support,
* ELF shared library overhead, etc. It also makes McStas generated
* simulations more portable (more likely to behave identically across
* platforms, important for parrallel computations).
*******************************************************************************/


#define        TYPE_3                3
#define        BREAK_3                128
#define        DEG_3                31
#define        SEP_3                3

static mc_int32_t randtbl[DEG_3 + 1] =
  {
    TYPE_3,

    -1726662223, 379960547, 1735697613, 1040273694, 1313901226,
    1627687941, -179304937, -2073333483, 1780058412, -1989503057,
    -615974602, 344556628, 939512070, -1249116260, 1507946756,
    -812545463, 154635395, 1388815473, -1926676823, 525320961,
    -1009028674, 968117788, -123449607, 1284210865, 435012392,
    -2017506339, -911064859, -370259173, 1132637927, 1398500161,
    -205601318,
  };

static mc_int32_t *fptr = &randtbl[SEP_3 + 1];
static mc_int32_t *rptr = &randtbl[1];
static mc_int32_t *state = &randtbl[1];
#define rand_deg DEG_3
#define rand_sep SEP_3
static mc_int32_t *end_ptr = &randtbl[sizeof (randtbl) / sizeof (randtbl[0])];

mc_int32_t
mc_random (void)
{
  mc_int32_t result;

  *fptr += *rptr;
  /* Chucking least random bit.  */
  result = (*fptr >> 1) & 0x7fffffff;
  ++fptr;
  if (fptr >= end_ptr)
  {
    fptr = state;
    ++rptr;
  }
  else
  {
    ++rptr;
    if (rptr >= end_ptr)
      rptr = state;
  }
  return result;
}

void
mc_srandom (unsigned int x)
{
  /* We must make sure the seed is not 0.  Take arbitrarily 1 in this case.  */
  state[0] = x ? x : 1;
  {
    long int i;
    for (i = 1; i < rand_deg; ++i)
    {
      /* This does:
         state[i] = (16807 * state[i - 1]) % 2147483647;
         but avoids overflowing 31 bits.  */
      long int hi = state[i - 1] / 127773;
      long int lo = state[i - 1] % 127773;
      long int test = 16807 * lo - 2836 * hi;
      state[i] = test + (test < 0 ? 2147483647 : 0);
    }
    fptr = &state[rand_sep];
    rptr = &state[0];
    for (i = 0; i < 10 * rand_deg; ++i)
      random ();
  }
}

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

static unsigned long mt[N]; /* the array for the state vector  */
static int mti=N+1; /* mti==N+1 means mt[N] is not initialized */

/* initializes mt[N] with a seed */
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

/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
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
    static unsigned long mag01[2]={0x0UL, MATRIX_A};
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

/* End of McCode random number routine. */

/* randnorm: generate a random number from normal law */
double
randnorm(void)
{
  static double v1, v2, s;
  static int phase = 0;
  double X, u1, u2;

  if(phase == 0)
  {
    do
    {
      u1 = rand01();
      u2 = rand01();
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

/**
 * Generate a random number from -1 to 1 with triangle distribution
 */
double randtriangle(void) {
	double randnum = rand01();
	if (randnum>0.5) return(1-sqrt(2*(randnum-0.5)));
	else return(sqrt(2*randnum)-1);
}

/**
 * Random number between 0.0 and 1.0 (including?)
 */
double rand01() {
	double randnum;
	randnum = (double) random();
	randnum /= (double) MC_RAND_MAX + 1;
	return randnum;
}

/**
 * Return a random number between 1 and -1
 */
double randpm1() {
	double randnum;
	randnum = (double) random();
	randnum /= ((double) MC_RAND_MAX + 1) / 2;
	randnum -= 1;
	return randnum;
}

/**
 * Return a random number between 0 and max.
 */
double rand0max(double max) {
	double randnum;
	randnum = (double) random();
	randnum /= ((double) MC_RAND_MAX + 1) / max;
	return randnum;
}

/**
 * Return a random number between min and max.
 */
double randminmax(double min, double max) {
	return rand0max(max - min) + max;
}

/* SECTION: main and signal handlers ======================================== */

/*******************************************************************************
* mchelp: displays instrument executable help with possible options
*******************************************************************************/
static void
mchelp(char *pgmname)
{
  int i;

  fprintf(stderr, "%s (%s) instrument simulation, generated with " MCCODE_STRING " (" MCCODE_DATE ")\n", mcinstrument_name, mcinstrument_source);
  fprintf(stderr, "Usage: %s [options] [parm=value ...]\n", pgmname);
  fprintf(stderr,
"Options are:\n"
"  -s SEED   --seed=SEED      Set random seed (must be != 0)\n"
"  -n COUNT  --ncount=COUNT   Set number of " MCCODE_PARTICLE "s to simulate.\n"
"  -d DIR    --dir=DIR        Put all data files in directory DIR.\n"
"  -t        --trace          Enable trace of " MCCODE_PARTICLE "s through instrument.\n"
"  -g        --gravitation    Enable gravitation for all trajectories.\n"
"  --no-output-files          Do not write any data files.\n"
"  -h        --help           Show this help message.\n"
"  -i        --info           Detailed instrument information.\n"
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
  if(mcnumipar > 0)
  {
    fprintf(stderr, "Instrument parameters are:\n");
    for(i = 0; i < mcnumipar; i++)
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
 if(mctraceenabled)
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
                    mcinstrument_name, mcinstrument_source));

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

  /* Add one to mcnumipar to avoid allocating zero size memory block. */
  paramsetarray = (int*)malloc((mcnumipar + 1)*sizeof(*paramsetarray));
  if(paramsetarray == NULL)
  {
    fprintf(stderr, "Error: insufficient memory (mcparseoptions)\n");
    exit(1);
  }
  for(j = 0; j < mcnumipar; j++)
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
    else if(!strcmp("--help", argv[i]))
      mcshowhelp(argv[0]);
    else if(!strcmp("-i", argv[i])) {
      mcformat=FLAVOR_UPPER;
      mcinfo();
    }
    else if(!strcmp("--info", argv[i]))
      mcinfo();
    else if(!strcmp("-t", argv[i]))
      mcenabletrace();
    else if(!strcmp("--trace", argv[i]))
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
    else if(argv[i][0] != '-' && (p = strchr(argv[i], '=')) != NULL)
    {
      *p++ = '\0';

      for(j = 0; j < mcnumipar; j++)
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
      if(j == mcnumipar)
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
    for(j = 0; j < mcnumipar; j++)
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
mcstatic char  mcsig_message[256];


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
  printf("# Simulation: %s (%s) \n", mcinstrument_name, mcinstrument_source);
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
    mcsave(NULL);
    fflush(stdout);
    return;
  }
  else
  if (sig == SIG_TERM)
  {
    printf("# " MCCODE_STRING ": Finishing simulation (save results and exit)\n");
    mcfinally();
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

/*******************************************************************************
* mccode_main: McCode main() function.
*******************************************************************************/
int mccode_main(int argc, char *argv[])
{
/*  double run_num = 0; */
  time_t  t;
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
  MPI_Comm_set_name(MPI_COMM_WORLD, mcinstrument_name);
  MPI_Get_processor_name(mpi_node_name, &mpi_node_name_len);
#endif /* USE_MPI */

t = time(NULL);
mcseed = (long)t+(long)getpid();

#ifdef USE_MPI
/* *** print number of nodes *********************************************** */
  if (mpi_node_count > 1) {
    MPI_MASTER(
    printf("Simulation '%s' (%s): running on %i nodes (master is '%s', MPI version %i.%i).\n",
      mcinstrument_name, mcinstrument_source, mpi_node_count, mpi_node_name, MPI_VERSION, MPI_SUBVERSION);
    );
  }
#endif /* USE_MPI */
  
  mcstartdate = (long)t;  /* set start date before parsing options and creating sim file */

/* *** parse options ******************************************************* */
  SIG_MESSAGE("main (Start)");
  mcformat=getenv(FLAVOR_UPPER "_FORMAT") ?
           getenv(FLAVOR_UPPER "_FORMAT") : FLAVOR_UPPER;
  mcinstrument_exe = argv[0]; /* store the executable path */
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
  mcsiminfo_init(NULL); /* open SIM */
  SIG_MESSAGE("main (Init)");
  mcinit();
#ifndef NOSIGNALS
#ifdef SIGINT
  if (signal( SIGINT ,sighandler) == SIG_IGN)
    signal( SIGINT,SIG_IGN);    /* interrupt (rubout) only after INIT */
#endif
#endif /* !NOSIGNALS */

/* ================ main particle generation/propagation loop ================ */
#if defined (USE_MPI)
  /* sliced Ncount on each MPI node */
  mcncount = mpi_node_count > 1 ?
    floor(mcncount / mpi_node_count) :
    mcncount; /* number of rays per node */
#endif

/* main particle event loop */
while(mcrun_num < mcncount || mcrun_num < mcget_ncount())
  {
#ifndef NEUTRONICS
    mcgenstate();
#endif
    /* old init: mcsetstate(0, 0, 0, 0, 0, 1, 0, sx=0, sy=1, sz=0, 1); */
    mcraytrace();
    mcrun_num++;
  }

#ifdef USE_MPI
 /* merge run_num from MPI nodes */
  if (mpi_node_count > 1) {
  double mcrun_num_double = (double)mcrun_num;
  mc_MPI_Sum(&mcrun_num_double, 1);
  mcrun_num = (unsigned long long)mcrun_num_double;
  }
#endif

/* save/finally executed by master node/thread */
  mcfinally();

#ifdef USE_MPI
  MPI_Finalize();
#endif /* USE_MPI */

  return 0;
} /* mccode_main */

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
  mcinit();

  /* *** parse options *** */
  SIG_MESSAGE("main (Start)");
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

#line 4977 "./Test_Monochromators.c"

#line 1 "mcstas-r.c"
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
* mcstore_neutron: stores neutron coodinates into global array (per component)
*******************************************************************************/
void
mcstore_neutron(MCNUM *s, int index, double x, double y, double z,
               double vx, double vy, double vz, double t,
               double sx, double sy, double sz, double p)
{
    double *dptr = &s[11*index];
    *dptr++  = x;
    *dptr++  = y ;
    *dptr++  = z ;
    *dptr++  = vx;
    *dptr++  = vy;
    *dptr++  = vz;
    *dptr++  = t ;
    *dptr++  = sx;
    *dptr++  = sy;
    *dptr++  = sz;
    *dptr    = p ;
} /* mcstore_neutron */

/*******************************************************************************
* mcrestore_neutron: restores neutron coodinates from global array
*******************************************************************************/
void
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
} /* mcrestore_neutron */

/*******************************************************************************
* mcsetstate: transfer parameters into global McStas variables 
*******************************************************************************/
void
mcsetstate(double x, double y, double z, double vx, double vy, double vz,
           double t, double sx, double sy, double sz, double p)
{
  extern double mcnx, mcny, mcnz, mcnvx, mcnvy, mcnvz;
  extern double mcnt, mcnsx, mcnsy, mcnsz, mcnp;

  mcnx = x;
  mcny = y;
  mcnz = z;
  mcnvx = vx;
  mcnvy = vy;
  mcnvz = vz;
  mcnt = t;
  mcnsx = sx;
  mcnsy = sy;
  mcnsz = sz;
  mcnp = p;
} /* mcsetstate */

/*******************************************************************************
* mcgenstate: set default neutron parameters 
*******************************************************************************/
void
mcgenstate(void)
{
  mcsetstate(0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1);
  /* old initialisation: mcsetstate(0, 0, 0, 0, 0, 1, 0, sx=0, sy=1, sz=0, 1); */
}

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
int
cylinder_intersect(double *t0, double *t1, double x, double y, double z,
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
int
sphere_intersect(double *t0, double *t1, double x, double y, double z,
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
int
plane_intersect(double *t, double x, double y, double z,
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

#line 5337 "./Test_Monochromators.c"
#ifdef MC_TRACE_ENABLED
int mctraceenabled = 1;
#else
int mctraceenabled = 0;
#endif
#define MCSTAS "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../"
int mcdefaultmain = 1;
char mcinstrument_name[] = "Test_Monochromators";
char mcinstrument_source[] = "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr";
char *mcinstrument_exe=NULL; /* will be set to argv[0] in main */
int main(int argc, char *argv[]){return mccode_main(argc, argv);}
void mcinit(void);
void mcraytrace(void);
void mcsave(FILE *);
void mcfinally(void);
void mcdisplay(void);

/* Shared user declarations for all components 'Monochromator_flat'. */
#line 71 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Monochromator_flat.comp"
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
#line 5377 "./Test_Monochromators.c"

/* Shared user declarations for all components 'Monochromator_pol'. */
#line 95 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Monochromator_pol.comp"
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
#define mc_pol_mu0=4*M_PI*1e-7

/*example field functions should have a variable set of arguments*/
#include <stdarg.h>
#include <stddef.h>
/*macros for some stuff*/
#ifndef MCSTAS_R_H
#include <mcstas-r.h>
#endif

/*Redefine MAGNET_OFF to also clear the magnetic field stack.*/
#undef MAGNET_OFF
#define MAGNET_OFF \
  do { \
    mcMagnet = 0; \
    mcmagnet_pop_all(); \
  } while(0)

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
    mcMagnet=0; \
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
void *mcmagnet_pop_all(void);

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

/*maximum number of rows to rebin a table = 1M*/
enum { mcread_table_rebin_maxsize = 1000000 };

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
* Version: $Revision: 5052 $
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

    /* logic here is Read_Table should include a call to FIND. If found the return value should just be used as
     * if the table had been read from disk. If not found then read the table and STORE.
     * Table_Free should include a call to GC. If this returns non-NULL then we should proceed with freeing the memory
     * associated with the table item - otherwise only decrement the reference counter since there are more references
     * that may need it.*/

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
            tr->table_ref = ((t_Table *) item);
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
                        /*the item is found and no garbage collection needed*/
                        tr->ref_count--;
                        return NULL;
                    }else{
                        /* The item is found and the reference counter is 1.
                         * This means we should garbage collect. Move remaining list items up one slot,
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
            /* item not found, and so should be garbage collected. This could be the case if freeing a
             * Table that has been constructed from code - not read from file. Return 0x1 to flag it for
             * collection.*/
            return (void *) 0x1 ;
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
    return Table_File_List_Handler(STORE,tab,0);
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

      if (!hfile && mcinstrument_source[0] != '\0' && strlen(mcinstrument_source)) /* search in instrument source location */
      {
        char *path_pos   = NULL;
        /* extract path: searches for last file separator */
        path_pos    = strrchr(mcinstrument_source, MC_PATHSEP_C);  /* last PATHSEP */
        if (path_pos) {
          long path_length = path_pos +1 - mcinstrument_source;  /* from start to path+sep */
          if (path_length) {
            strncpy(dir, mcinstrument_source, path_length);
            dir[path_length] = '\0';
            snprintf(path, 1024, "%s%c%s", dir, MC_PATHSEP_C, File);
            hfile = fopen(path, Mode);
          }
        }
      }
      if (!hfile && mcinstrument_exe[0] != '\0' && strlen(mcinstrument_exe)) /* search in PWD instrument executable location */
      {
        char *path_pos   = NULL;
        /* extract path: searches for last file separator */
        path_pos    = strrchr(mcinstrument_exe, MC_PATHSEP_C);  /* last PATHSEP */
        if (path_pos) {
          long path_length = path_pos +1 - mcinstrument_exe;  /* from start to path+sep */
          if (path_length) {
            strncpy(dir, mcinstrument_exe, path_length);
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
                    malloc_size = count_in_array*1.5;
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
        (Table->filename[0] != '\0' ? Table->filename : ""),
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
      /*return early if the rebinned table will become too large*/
      if (Length_Table > mcread_table_rebin_maxsize){
        fprintf(stderr,"WARNING: (Table_Rebin): Rebinning table from %s would exceed 1M rows. Skipping.\n", Table->filename); 
        return(Table->rows*Table->columns);
      }
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
*   ACTION: free a single Table. First Call Table_File_list_gc. If this returns
*   non-zero it means there are more refernces to the table, and so the table
*   should not bee freed.
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
      Table.filename[0] != '\0' ? Table.filename : "", buffer);
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

    /* first allocate an initial empty t_Table array */
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
      /*if the block is empty - don't store it*/
      if (nelements>0){
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
      }
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
    long index;
    if (!Table) return;
    for (index=0;index < Table[0].array_length; index++){
            Table_Free(&Table[index]);
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
    for (i=0; i<interpolator->field_dimensionality; i++){
        field[i]=w->data[i];
    }
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
    mcMagnetPrecession(mcnlx, mcnly, mcnlz, mcnlt, mcnlvx, mcnlvy, mcnlvz, \
	   	       &mcnlsx, &mcnlsy, &mcnlsz, dt, posLM, rotLM); \
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
    mcMagnet=1;
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
    mcMagnet=1;
  }else{
    mcMagnet=0;
  }
  return (void*) t;
}

void *mcmagnet_pop_all(void){
  void *t=mcmagnet_pop();
  while (t!=NULL){
    t=mcmagnet_pop();
  }
  return NULL;
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
  *bx =  Bmagnitude * sin(PI/2*(z+magnetLength/2.0)/magnetLength);
  *by =  Bmagnitude * cos(PI/2*(z+magnetLength/2.0)/magnetLength);
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
  mccoordschange(mc_pol_posLM, mc_pol_rotLM,
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

#line 8215 "./Test_Monochromators.c"

/* Shared user declarations for all components 'Monochromator_curved'. */
#line 109 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Monochromator_curved.comp"
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


#line 8241 "./Test_Monochromators.c"

/* Shared user declarations for all components 'Monochromator_2foc'. */
#line 83 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../contrib/Monochromator_2foc.comp"

#line 8246 "./Test_Monochromators.c"

/* Shared user declarations for all components 'Single_crystal'. */
#line 296 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../samples/Single_crystal.comp"
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
              "       The first line is required to be exactly 'OFF', '3' or 'ply'.\n",filename,strlen(line));
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
      fgets(line, CHAR_BUF_LENGTH, f);
      continue; 
    }
    if (ret != 3) {
      fprintf(stderr, "Error: can not read [xyz] coordinates for vertex %ld in file %s (interoff/off_init). Read %ld values.\n", 
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
      fgets(line, CHAR_BUF_LENGTH, f);
      continue; 
    }
    if (ret != 1) {
      fprintf(stderr, "Error: can not read polygon %ld length in file %s (interoff/off_init)\n", 
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
      fscanf(f, "%lg", &vtx);
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
    sprintf(pixelinfo, "%lu,%lu,%lu,%i,%g,%g,%g,%g,%g,%g", data.mantidoffset+pixel, data.mantidoffset, data.mantidoffset+data.polySize-1, nbVertex, cmx, cmy, cmz, x1-cmx, y1-cmy, z1-cmz);
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
      normal_vec(b1[0], b1[1], b1[2],
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
          normal_vec(b1x, b1y, b1z, nx, ny, nz);
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
    void randrotate(double *nx, double *ny, double *nz, double a, double b, double c) {
      double x1, y1, z1, x2, y2, z2;
      rotate(x1, y1, z1, *nx,*ny,*nz, a, 1, 0, 0); /* <1> = rot(<n>,a) */
      rotate(x2, y2, z2,  x1, y1, z1, b, 0, 1, 0); /* <2> = rot(<1>,b) */
      rotate(*nx,*ny,*nz, x2, y2, z2, c, 0, 0, 1); /* <n> = rot(<2>,c) */
    }
    /* Powder, back */
    void randderotate(double *nx, double *ny, double *nz, double a, double b, double c) {
      double x1, y1, z1, x2, y2, z2;
      rotate(x1, y1, z1, *nx,*ny,*nz, -c, 0,0,1);
      rotate(x2, y2, z2,  x1, y1, z1, -b, 0,1,0);
      rotate(*nx,*ny,*nz, x2, y2, z2, -a, 1,0,0);
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



    /* rotate vector counterclockwise */
    void vec_rotate_2d(double* x, double* y, double angle) {
        double c, s;
        double newx, newy;

        c = cos(angle);
        s = sin(angle);

        newx = *x*c - *y*s;
        newy = *x*s + *y*c;

        *x = newx;
        *y = newy;
    }


#endif /* !SINGLE_CRYSTAL_DECL */
#line 9909 "./Test_Monochromators.c"

/* Shared user declarations for all components 'PSD_monitor'. */
#line 57 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"

#ifndef ARRAYS_H
#define ARRAYS_H
typedef double* DArray1d;
DArray1d create_darr1d(int n);
void destroy_darr1d(DArray1d a);

typedef double** DArray2d;
DArray2d create_darr2d(int nx, int ny);
void destroy_darr2d(DArray2d a);

typedef double*** DArray3d;
DArray3d create_darr3d(int nx, int ny, int nz);
void destroy_darr3d(DArray3d a);
#endif
#ifndef ARRAYS_C
#define ARRAYS_C
#include <stdlib.h>

typedef double* DArray1d;
typedef double** DArray2d;
typedef double*** DArray3d;

DArray1d create_darr1d(int n){
  DArray1d arr2d;
  arr2d = calloc(n, sizeof(double));
  return arr2d;
}
void destroy_darr1d(DArray1d a){
  free(a);
}

DArray2d create_darr2d(int nx, int ny){
  DArray2d arr2d;
  arr2d = calloc(nx, sizeof(double *));

  double *p1;
  p1 = calloc(nx*ny, sizeof(double));

  int i;
  for (i=0; i<nx; i++){
    arr2d[i] = &(p1[i*ny]);
  }
  return arr2d;
}
void destroy_darr2d(DArray2d a){
  free(a[0]);
  free(a);
}

DArray3d create_darr3d(int nx, int ny, int nz){
  DArray3d arr3d;
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
  }
  return arr3d;
}
void destroy_darr3d(DArray3d a){
  free(a[0][0]);
  free(a[0]);
  free(a);
}
#endif

#line 9995 "./Test_Monochromators.c"

/* Instrument parameters. */
int mcipMono;
MCNUM mciplambda;
MCNUM mcipMoz;
MCNUM mcipPG;
MCNUM mcippowder;
MCNUM mciporder;

#define mcNUMIPAR 6
int mcnumipar = 6;
struct mcinputtable_struct mcinputtable[mcNUMIPAR+1] = {
  "Mono", &mcipMono, instr_type_int, "1", 
  "lambda", &mciplambda, instr_type_double, "2.0", 
  "Moz", &mcipMoz, instr_type_double, "40", 
  "PG", &mcipPG, instr_type_double, "0", 
  "powder", &mcippowder, instr_type_double, "0", 
  "order", &mciporder, instr_type_double, "1", 
  NULL, NULL, instr_type_double, ""
};

/* User declarations from instrument definition. */
#define mccompcurname  Test_Monochromators
#define mccompcurtype  INSTRUMENT
#define mccompcurindex 0
#define mcposaTest_Monochromators coords_set(0,0,0)
#define Mono mcipMono
#define lambda mciplambda
#define Moz mcipMoz
#define PG mcipPG
#define powder mcippowder
#define order mciporder
#line 46 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  double DM  = 3.3539; /* Monochromator d-spacing in Angs */
                       /* PG002 Orders : 1st 3.355 2e 1.6775, 3e 1.1183 */

  /* to compute */
  double A1,A2;
  double mono_q;

  /* This variable helps us switch on and off the different setups*/
  double filterProb;
#line 10038 "./Test_Monochromators.c"
#undef order
#undef powder
#undef PG
#undef Moz
#undef lambda
#undef Mono
#undef mcposaTest_Monochromators
#undef mccompcurindex
#undef mccompcurtype
#undef mccompcurname

/* neutron state table at each component input (local coords) */
/* [x, y, z, vx, vy, vz, t, sx, sy, sz, p] */
MCNUM mccomp_storein[11*16];
/* Components position table (absolute and relative coords) */
Coords mccomp_posa[16];
Coords mccomp_posr[16];
/* Counter for each comp to check for inactive ones */
MCNUM  mcNCounter[16];
MCNUM  mcPCounter[16];
MCNUM  mcP2Counter[16];
#define mcNUMCOMP 15 /* number of components */
/* Counter for PROP ABSORB */
MCNUM  mcAbsorbProp[16];
/* Flag true when previous component acted on the neutron (SCATTER) */
MCNUM mcScattered=0;
/* Flag true when neutron should be restored (RESTORE) */
MCNUM mcRestore=0;
/* Declarations of component definition and setting parameters. */

/* Setting parameters for component 'Origin' [1]. */
char mccOrigin_profile[16384];
MCNUM mccOrigin_percent;
MCNUM mccOrigin_flag_save;
MCNUM mccOrigin_minutes;

/* Setting parameters for component 'Source' [2]. */
MCNUM mccSource_radius;
MCNUM mccSource_yheight;
MCNUM mccSource_xwidth;
MCNUM mccSource_dist;
MCNUM mccSource_focus_xw;
MCNUM mccSource_focus_yh;
MCNUM mccSource_E0;
MCNUM mccSource_dE;
MCNUM mccSource_lambda0;
MCNUM mccSource_dlambda;
MCNUM mccSource_flux;
MCNUM mccSource_gauss;
int mccSource_target_index;

/* Definition parameters for component 'lamStart' [3]. */
#define mcclamStart_nL 200
/* Setting parameters for component 'lamStart' [3]. */
char mcclamStart_filename[16384];
MCNUM mcclamStart_xmin;
MCNUM mcclamStart_xmax;
MCNUM mcclamStart_ymin;
MCNUM mcclamStart_ymax;
MCNUM mcclamStart_xwidth;
MCNUM mcclamStart_yheight;
MCNUM mcclamStart_Lmin;
MCNUM mcclamStart_Lmax;
MCNUM mcclamStart_restore_neutron;
int mcclamStart_nowritefile;

/* Setting parameters for component 'Mono1' [6]. */
MCNUM mccMono1_zmin;
MCNUM mccMono1_zmax;
MCNUM mccMono1_ymin;
MCNUM mccMono1_ymax;
MCNUM mccMono1_zwidth;
MCNUM mccMono1_yheight;
MCNUM mccMono1_mosaich;
MCNUM mccMono1_mosaicv;
MCNUM mccMono1_r0;
MCNUM mccMono1_Q;
MCNUM mccMono1_DM;

/* Setting parameters for component 'Mono2' [7]. */
MCNUM mccMono2_zwidth;
MCNUM mccMono2_yheight;
MCNUM mccMono2_mosaic;
MCNUM mccMono2_dspread;
MCNUM mccMono2_Q;
MCNUM mccMono2_DM;
MCNUM mccMono2_pThreshold;
MCNUM mccMono2_Rup;
MCNUM mccMono2_Rdown;
int mccMono2_debug;

/* Setting parameters for component 'Mono3' [8]. */
MCNUM mccMono3_zwidth;
MCNUM mccMono3_yheight;
MCNUM mccMono3_mosaic;
MCNUM mccMono3_dspread;
MCNUM mccMono3_Q;
MCNUM mccMono3_DM;
MCNUM mccMono3_pThreshold;
MCNUM mccMono3_Rup;
MCNUM mccMono3_Rdown;
int mccMono3_debug;

/* Setting parameters for component 'Mono4' [9]. */
char mccMono4_reflect[16384];
char mccMono4_transmit[16384];
MCNUM mccMono4_zwidth;
MCNUM mccMono4_yheight;
MCNUM mccMono4_gap;
MCNUM mccMono4_NH;
MCNUM mccMono4_NV;
MCNUM mccMono4_mosaich;
MCNUM mccMono4_mosaicv;
MCNUM mccMono4_r0;
MCNUM mccMono4_t0;
MCNUM mccMono4_Q;
MCNUM mccMono4_RV;
MCNUM mccMono4_RH;
MCNUM mccMono4_DM;
MCNUM mccMono4_mosaic;
MCNUM mccMono4_width;
MCNUM mccMono4_height;
MCNUM mccMono4_verbose;
MCNUM mccMono4_order;

/* Setting parameters for component 'Mono5' [10]. */
char mccMono5_reflect[16384];
MCNUM mccMono5_zwidth;
MCNUM mccMono5_yheight;
MCNUM mccMono5_gap;
MCNUM mccMono5_NH;
MCNUM mccMono5_NV;
MCNUM mccMono5_mosaich;
MCNUM mccMono5_mosaicv;
MCNUM mccMono5_r0;
MCNUM mccMono5_Q;
MCNUM mccMono5_RV;
MCNUM mccMono5_RH;
MCNUM mccMono5_DM;
MCNUM mccMono5_mosaic;
MCNUM mccMono5_width;
MCNUM mccMono5_height;
MCNUM mccMono5_verbose;

/* Definition parameters for component 'Mono6' [11]. */
#define mccMono6_mosaic_AB Mosaic_AB_Undefined
/* Setting parameters for component 'Mono6' [11]. */
char mccMono6_reflections[16384];
char mccMono6_geometry[16384];
MCNUM mccMono6_xwidth;
MCNUM mccMono6_yheight;
MCNUM mccMono6_zdepth;
MCNUM mccMono6_radius;
MCNUM mccMono6_delta_d_d;
MCNUM mccMono6_mosaic;
MCNUM mccMono6_mosaic_a;
MCNUM mccMono6_mosaic_b;
MCNUM mccMono6_mosaic_c;
MCNUM mccMono6_recip_cell;
MCNUM mccMono6_barns;
MCNUM mccMono6_ax;
MCNUM mccMono6_ay;
MCNUM mccMono6_az;
MCNUM mccMono6_bx;
MCNUM mccMono6_by;
MCNUM mccMono6_bz;
MCNUM mccMono6_cx;
MCNUM mccMono6_cy;
MCNUM mccMono6_cz;
MCNUM mccMono6_p_transmit;
MCNUM mccMono6_sigma_abs;
MCNUM mccMono6_sigma_inc;
MCNUM mccMono6_aa;
MCNUM mccMono6_bb;
MCNUM mccMono6_cc;
MCNUM mccMono6_order;
MCNUM mccMono6_RX;
MCNUM mccMono6_RY;
MCNUM mccMono6_powder;
MCNUM mccMono6_PG;
MCNUM mccMono6_deltak;

/* Definition parameters for component 'Sphere1' [12]. */
#define mccSphere1_nx 90
#define mccSphere1_ny 90
/* Setting parameters for component 'Sphere1' [12]. */
char mccSphere1_filename[16384];
MCNUM mccSphere1_radius;
MCNUM mccSphere1_restore_neutron;
int mccSphere1_nowritefile;

/* Definition parameters for component 'lam1' [13]. */
#define mcclam1_nL 200
/* Setting parameters for component 'lam1' [13]. */
char mcclam1_filename[16384];
MCNUM mcclam1_xmin;
MCNUM mcclam1_xmax;
MCNUM mcclam1_ymin;
MCNUM mcclam1_ymax;
MCNUM mcclam1_xwidth;
MCNUM mcclam1_yheight;
MCNUM mcclam1_Lmin;
MCNUM mcclam1_Lmax;
MCNUM mcclam1_restore_neutron;
int mcclam1_nowritefile;

/* Setting parameters for component 'psd1' [14]. */
int mccpsd1_nx;
int mccpsd1_ny;
char mccpsd1_filename[16384];
MCNUM mccpsd1_xmin;
MCNUM mccpsd1_xmax;
MCNUM mccpsd1_ymin;
MCNUM mccpsd1_ymax;
MCNUM mccpsd1_xwidth;
MCNUM mccpsd1_yheight;
MCNUM mccpsd1_restore_neutron;

/* User component declarations. */

/* User declarations for component 'Origin' [1]. */
#define mccompcurname  Origin
#define mccompcurtype  Progress_bar
#define mccompcurindex 1
#define IntermediateCnts mccOrigin_IntermediateCnts
#define StartTime mccOrigin_StartTime
#define EndTime mccOrigin_EndTime
#define CurrentTime mccOrigin_CurrentTime
#define profile mccOrigin_profile
#define percent mccOrigin_percent
#define flag_save mccOrigin_flag_save
#define minutes mccOrigin_minutes
#line 44 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../misc/Progress_bar.comp"
#ifndef PROGRESS_BAR
#define PROGRESS_BAR
#else
#error Only one Progress_bar component may be used in an instrument definition.
#endif

double IntermediateCnts;
time_t StartTime;
time_t EndTime;
time_t CurrentTime;
#line 10282 "./Test_Monochromators.c"
#undef minutes
#undef flag_save
#undef percent
#undef profile
#undef CurrentTime
#undef EndTime
#undef StartTime
#undef IntermediateCnts
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'Source' [2]. */
#define mccompcurname  Source
#define mccompcurtype  Source_simple
#define mccompcurindex 2
#define pmul mccSource_pmul
#define square mccSource_square
#define srcArea mccSource_srcArea
#define radius mccSource_radius
#define yheight mccSource_yheight
#define xwidth mccSource_xwidth
#define dist mccSource_dist
#define focus_xw mccSource_focus_xw
#define focus_yh mccSource_focus_yh
#define E0 mccSource_E0
#define dE mccSource_dE
#define lambda0 mccSource_lambda0
#define dlambda mccSource_dlambda
#define flux mccSource_flux
#define gauss mccSource_gauss
#define target_index mccSource_target_index
#line 60 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../sources/Source_simple.comp"
double pmul, srcArea;
int square;
double tx,ty,tz;
#line 10319 "./Test_Monochromators.c"
#undef target_index
#undef gauss
#undef flux
#undef dlambda
#undef lambda0
#undef dE
#undef E0
#undef focus_yh
#undef focus_xw
#undef dist
#undef xwidth
#undef yheight
#undef radius
#undef srcArea
#undef square
#undef pmul
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'lamStart' [3]. */
#define mccompcurname  lamStart
#define mccompcurtype  L_monitor
#define mccompcurindex 3
#define nL mcclamStart_nL
#define L_N mcclamStart_L_N
#define L_p mcclamStart_L_p
#define L_p2 mcclamStart_L_p2
#define filename mcclamStart_filename
#define xmin mcclamStart_xmin
#define xmax mcclamStart_xmax
#define ymin mcclamStart_ymin
#define ymax mcclamStart_ymax
#define xwidth mcclamStart_xwidth
#define yheight mcclamStart_yheight
#define Lmin mcclamStart_Lmin
#define Lmax mcclamStart_Lmax
#define restore_neutron mcclamStart_restore_neutron
#define nowritefile mcclamStart_nowritefile
#line 57 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/L_monitor.comp"
double L_N[nL];
double L_p[nL], L_p2[nL];
#line 10362 "./Test_Monochromators.c"
#undef nowritefile
#undef restore_neutron
#undef Lmax
#undef Lmin
#undef yheight
#undef xwidth
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef filename
#undef L_p2
#undef L_p
#undef L_N
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'Mono_Arm' [4]. */
#define mccompcurname  Mono_Arm
#define mccompcurtype  Arm
#define mccompcurindex 4
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'Mono_Out' [5]. */
#define mccompcurname  Mono_Out
#define mccompcurtype  Arm
#define mccompcurindex 5
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'Mono1' [6]. */
#define mccompcurname  Mono1
#define mccompcurtype  Monochromator_flat
#define mccompcurindex 6
#define mos_rms_y mccMono1_mos_rms_y
#define mos_rms_z mccMono1_mos_rms_z
#define mos_rms_max mccMono1_mos_rms_max
#define mono_Q mccMono1_mono_Q
#define zmin mccMono1_zmin
#define zmax mccMono1_zmax
#define ymin mccMono1_ymin
#define ymax mccMono1_ymax
#define zwidth mccMono1_zwidth
#define yheight mccMono1_yheight
#define mosaich mccMono1_mosaich
#define mosaicv mccMono1_mosaicv
#define r0 mccMono1_r0
#define Q mccMono1_Q
#define DM mccMono1_DM
#line 95 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Monochromator_flat.comp"
  double mos_rms_y; /* root-mean-square of mosaic, in radians */
  double mos_rms_z;
  double mos_rms_max;
  double mono_Q;
#line 10422 "./Test_Monochromators.c"
#undef DM
#undef Q
#undef r0
#undef mosaicv
#undef mosaich
#undef yheight
#undef zwidth
#undef ymax
#undef ymin
#undef zmax
#undef zmin
#undef mono_Q
#undef mos_rms_max
#undef mos_rms_z
#undef mos_rms_y
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'Mono2' [7]. */
#define mccompcurname  Mono2
#define mccompcurtype  Monochromator_pol
#define mccompcurindex 7
#define mos_rms mccMono2_mos_rms
#define d_rms mccMono2_d_rms
#define mono_Q mccMono2_mono_Q
#define FN mccMono2_FN
#define FM mccMono2_FM
#define zwidth mccMono2_zwidth
#define yheight mccMono2_yheight
#define mosaic mccMono2_mosaic
#define dspread mccMono2_dspread
#define Q mccMono2_Q
#define DM mccMono2_DM
#define pThreshold mccMono2_pThreshold
#define Rup mccMono2_Rup
#define Rdown mccMono2_Rdown
#define debug mccMono2_debug
#line 100 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Monochromator_pol.comp"
double mos_rms; /* root-mean-square of mosaic, in radians */
double d_rms;   /* root-mean-square of d-spread, in AA */
double mono_Q;

double FN; /* Unit cell nuclear structure factor */
double FM; /* Unit cell magnetic structure factor */
#line 10468 "./Test_Monochromators.c"
#undef debug
#undef Rdown
#undef Rup
#undef pThreshold
#undef DM
#undef Q
#undef dspread
#undef mosaic
#undef yheight
#undef zwidth
#undef FM
#undef FN
#undef mono_Q
#undef d_rms
#undef mos_rms
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'Mono3' [8]. */
#define mccompcurname  Mono3
#define mccompcurtype  Monochromator_pol
#define mccompcurindex 8
#define mos_rms mccMono3_mos_rms
#define d_rms mccMono3_d_rms
#define mono_Q mccMono3_mono_Q
#define FN mccMono3_FN
#define FM mccMono3_FM
#define zwidth mccMono3_zwidth
#define yheight mccMono3_yheight
#define mosaic mccMono3_mosaic
#define dspread mccMono3_dspread
#define Q mccMono3_Q
#define DM mccMono3_DM
#define pThreshold mccMono3_pThreshold
#define Rup mccMono3_Rup
#define Rdown mccMono3_Rdown
#define debug mccMono3_debug
#line 100 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Monochromator_pol.comp"
double mos_rms; /* root-mean-square of mosaic, in radians */
double d_rms;   /* root-mean-square of d-spread, in AA */
double mono_Q;

double FN; /* Unit cell nuclear structure factor */
double FM; /* Unit cell magnetic structure factor */
#line 10514 "./Test_Monochromators.c"
#undef debug
#undef Rdown
#undef Rup
#undef pThreshold
#undef DM
#undef Q
#undef dspread
#undef mosaic
#undef yheight
#undef zwidth
#undef FM
#undef FN
#undef mono_Q
#undef d_rms
#undef mos_rms
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'Mono4' [9]. */
#define mccompcurname  Mono4
#define mccompcurtype  Monochromator_curved
#define mccompcurindex 9
#define mos_rms_y mccMono4_mos_rms_y
#define mos_rms_z mccMono4_mos_rms_z
#define mos_rms_max mccMono4_mos_rms_max
#define mono_Q mccMono4_mono_Q
#define SlabWidth mccMono4_SlabWidth
#define SlabHeight mccMono4_SlabHeight
#define rTable mccMono4_rTable
#define tTable mccMono4_tTable
#define row mccMono4_row
#define col mccMono4_col
#define tiltH mccMono4_tiltH
#define tiltV mccMono4_tiltV
#define reflect mccMono4_reflect
#define transmit mccMono4_transmit
#define zwidth mccMono4_zwidth
#define yheight mccMono4_yheight
#define gap mccMono4_gap
#define NH mccMono4_NH
#define NV mccMono4_NV
#define mosaich mccMono4_mosaich
#define mosaicv mccMono4_mosaicv
#define r0 mccMono4_r0
#define t0 mccMono4_t0
#define Q mccMono4_Q
#define RV mccMono4_RV
#define RH mccMono4_RH
#define DM mccMono4_DM
#define mosaic mccMono4_mosaic
#define width mccMono4_width
#define height mccMono4_height
#define verbose mccMono4_verbose
#define order mccMono4_order
#line 136 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Monochromator_curved.comp"
  double mos_rms_y; /* root-mean-square of mosaic, in radians */
  double mos_rms_z;
  double mos_rms_max;
  double mono_Q;
  double SlabWidth, SlabHeight;
  t_Table rTable, tTable;
  double row,col;
  double* tiltH;
  double* tiltV;
#line 10580 "./Test_Monochromators.c"
#undef order
#undef verbose
#undef height
#undef width
#undef mosaic
#undef DM
#undef RH
#undef RV
#undef Q
#undef t0
#undef r0
#undef mosaicv
#undef mosaich
#undef NV
#undef NH
#undef gap
#undef yheight
#undef zwidth
#undef transmit
#undef reflect
#undef tiltV
#undef tiltH
#undef col
#undef row
#undef tTable
#undef rTable
#undef SlabHeight
#undef SlabWidth
#undef mono_Q
#undef mos_rms_max
#undef mos_rms_z
#undef mos_rms_y
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'Mono5' [10]. */
#define mccompcurname  Mono5
#define mccompcurtype  Monochromator_2foc
#define mccompcurindex 10
#define mos_y mccMono5_mos_y
#define mos_z mccMono5_mos_z
#define mono_Q mccMono5_mono_Q
#define SlabWidth mccMono5_SlabWidth
#define SlabHeight mccMono5_SlabHeight
#define rTable mccMono5_rTable
#define reflect mccMono5_reflect
#define zwidth mccMono5_zwidth
#define yheight mccMono5_yheight
#define gap mccMono5_gap
#define NH mccMono5_NH
#define NV mccMono5_NV
#define mosaich mccMono5_mosaich
#define mosaicv mccMono5_mosaicv
#define r0 mccMono5_r0
#define Q mccMono5_Q
#define RV mccMono5_RV
#define RH mccMono5_RH
#define DM mccMono5_DM
#define mosaic mccMono5_mosaic
#define width mccMono5_width
#define height mccMono5_height
#define verbose mccMono5_verbose
#line 88 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../contrib/Monochromator_2foc.comp"
#ifndef DIV_CUTOFF
#define DIV_CUTOFF 2            /* ~ 10^-5 cutoff. */
#endif
double mos_y; /* mosaic, in radians */
double mos_z;
double mono_Q;
double SlabWidth, SlabHeight;
t_Table rTable;
#line 10653 "./Test_Monochromators.c"
#undef verbose
#undef height
#undef width
#undef mosaic
#undef DM
#undef RH
#undef RV
#undef Q
#undef r0
#undef mosaicv
#undef mosaich
#undef NV
#undef NH
#undef gap
#undef yheight
#undef zwidth
#undef reflect
#undef rTable
#undef SlabHeight
#undef SlabWidth
#undef mono_Q
#undef mos_z
#undef mos_y
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'Mono6' [11]. */
#define mccompcurname  Mono6
#define mccompcurtype  Single_crystal
#define mccompcurindex 11
#define mosaic_AB mccMono6_mosaic_AB
#define hkl_info mccMono6_hkl_info
#define offdata mccMono6_offdata
#define reflections mccMono6_reflections
#define geometry mccMono6_geometry
#define xwidth mccMono6_xwidth
#define yheight mccMono6_yheight
#define zdepth mccMono6_zdepth
#define radius mccMono6_radius
#define delta_d_d mccMono6_delta_d_d
#define mosaic mccMono6_mosaic
#define mosaic_a mccMono6_mosaic_a
#define mosaic_b mccMono6_mosaic_b
#define mosaic_c mccMono6_mosaic_c
#define recip_cell mccMono6_recip_cell
#define barns mccMono6_barns
#define ax mccMono6_ax
#define ay mccMono6_ay
#define az mccMono6_az
#define bx mccMono6_bx
#define by mccMono6_by
#define bz mccMono6_bz
#define cx mccMono6_cx
#define cy mccMono6_cy
#define cz mccMono6_cz
#define p_transmit mccMono6_p_transmit
#define sigma_abs mccMono6_sigma_abs
#define sigma_inc mccMono6_sigma_inc
#define aa mccMono6_aa
#define bb mccMono6_bb
#define cc mccMono6_cc
#define order mccMono6_order
#define RX mccMono6_RX
#define RY mccMono6_RY
#define powder mccMono6_powder
#define PG mccMono6_PG
#define deltak mccMono6_deltak
#line 970 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../samples/Single_crystal.comp"
  struct hkl_info_struct hkl_info;
  off_struct             offdata;
#line 10725 "./Test_Monochromators.c"
#undef deltak
#undef PG
#undef powder
#undef RY
#undef RX
#undef order
#undef cc
#undef bb
#undef aa
#undef sigma_inc
#undef sigma_abs
#undef p_transmit
#undef cz
#undef cy
#undef cx
#undef bz
#undef by
#undef bx
#undef az
#undef ay
#undef ax
#undef barns
#undef recip_cell
#undef mosaic_c
#undef mosaic_b
#undef mosaic_a
#undef mosaic
#undef delta_d_d
#undef radius
#undef zdepth
#undef yheight
#undef xwidth
#undef geometry
#undef reflections
#undef offdata
#undef hkl_info
#undef mosaic_AB
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'Sphere1' [12]. */
#define mccompcurname  Sphere1
#define mccompcurtype  PSD_monitor_4PI
#define mccompcurindex 12
#define nx mccSphere1_nx
#define ny mccSphere1_ny
#define PSD_N mccSphere1_PSD_N
#define PSD_p mccSphere1_PSD_p
#define PSD_p2 mccSphere1_PSD_p2
#define filename mccSphere1_filename
#define radius mccSphere1_radius
#define restore_neutron mccSphere1_restore_neutron
#define nowritefile mccSphere1_nowritefile
#line 55 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor_4PI.comp"
double PSD_N[nx][ny];
double PSD_p[nx][ny];
double PSD_p2[nx][ny];
#line 10784 "./Test_Monochromators.c"
#undef nowritefile
#undef restore_neutron
#undef radius
#undef filename
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef ny
#undef nx
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'lam1' [13]. */
#define mccompcurname  lam1
#define mccompcurtype  L_monitor
#define mccompcurindex 13
#define nL mcclam1_nL
#define L_N mcclam1_L_N
#define L_p mcclam1_L_p
#define L_p2 mcclam1_L_p2
#define filename mcclam1_filename
#define xmin mcclam1_xmin
#define xmax mcclam1_xmax
#define ymin mcclam1_ymin
#define ymax mcclam1_ymax
#define xwidth mcclam1_xwidth
#define yheight mcclam1_yheight
#define Lmin mcclam1_Lmin
#define Lmax mcclam1_Lmax
#define restore_neutron mcclam1_restore_neutron
#define nowritefile mcclam1_nowritefile
#line 57 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/L_monitor.comp"
double L_N[nL];
double L_p[nL], L_p2[nL];
#line 10820 "./Test_Monochromators.c"
#undef nowritefile
#undef restore_neutron
#undef Lmax
#undef Lmin
#undef yheight
#undef xwidth
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef filename
#undef L_p2
#undef L_p
#undef L_N
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'psd1' [14]. */
#define mccompcurname  psd1
#define mccompcurtype  PSD_monitor
#define mccompcurindex 14
#define PSD_N mccpsd1_PSD_N
#define PSD_p mccpsd1_PSD_p
#define PSD_p2 mccpsd1_PSD_p2
#define nx mccpsd1_nx
#define ny mccpsd1_ny
#define filename mccpsd1_filename
#define xmin mccpsd1_xmin
#define xmax mccpsd1_xmax
#define ymin mccpsd1_ymin
#define ymax mccpsd1_ymax
#define xwidth mccpsd1_xwidth
#define yheight mccpsd1_yheight
#define restore_neutron mccpsd1_restore_neutron
#line 62 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
  DArray2d PSD_N;
  DArray2d PSD_p;
  DArray2d PSD_p2;
#line 10861 "./Test_Monochromators.c"
#undef restore_neutron
#undef yheight
#undef xwidth
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef filename
#undef ny
#undef nx
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

Coords mcposaOrigin, mcposrOrigin;
Rotation mcrotaOrigin, mcrotrOrigin;
Coords mcposaSource, mcposrSource;
Rotation mcrotaSource, mcrotrSource;
Coords mcposalamStart, mcposrlamStart;
Rotation mcrotalamStart, mcrotrlamStart;
Coords mcposaMono_Arm, mcposrMono_Arm;
Rotation mcrotaMono_Arm, mcrotrMono_Arm;
Coords mcposaMono_Out, mcposrMono_Out;
Rotation mcrotaMono_Out, mcrotrMono_Out;
Coords mcposaMono1, mcposrMono1;
Rotation mcrotaMono1, mcrotrMono1;
Coords mcposaMono2, mcposrMono2;
Rotation mcrotaMono2, mcrotrMono2;
Coords mcposaMono3, mcposrMono3;
Rotation mcrotaMono3, mcrotrMono3;
Coords mcposaMono4, mcposrMono4;
Rotation mcrotaMono4, mcrotrMono4;
Coords mcposaMono5, mcposrMono5;
Rotation mcrotaMono5, mcrotrMono5;
Coords mcposaMono6, mcposrMono6;
Rotation mcrotaMono6, mcrotrMono6;
Coords mcposaSphere1, mcposrSphere1;
Rotation mcrotaSphere1, mcrotrSphere1;
Coords mcposalam1, mcposrlam1;
Rotation mcrotalam1, mcrotrlam1;
Coords mcposapsd1, mcposrpsd1;
Rotation mcrotapsd1, mcrotrpsd1;

MCNUM mcnx, mcny, mcnz, mcnvx, mcnvy, mcnvz, mcnt, mcnsx, mcnsy, mcnsz, mcnp;

/* end declare */

void mcinit(void) {
#define mccompcurname  Test_Monochromators
#define mccompcurtype  INSTRUMENT
#define mccompcurindex 0
#define mcposaTest_Monochromators coords_set(0,0,0)
#define Mono mcipMono
#define lambda mciplambda
#define Moz mcipMoz
#define PG mcipPG
#define powder mcippowder
#define order mciporder
#line 59 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
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

  printf("\n%s: ", NAME_CURRENT_COMP);
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
#line 10957 "./Test_Monochromators.c"
#undef order
#undef powder
#undef PG
#undef Moz
#undef lambda
#undef Mono
#undef mcposaTest_Monochromators
#undef mccompcurindex
#undef mccompcurtype
#undef mccompcurname
  /* Computation of coordinate transformations. */
  {
    Coords mctc1, mctc2, mcLastComp;
    Rotation mctr1;
    double mcAccumulatedILength = 0;
    /* Initialize "last" component origin as (0,0,0) */
    mcLastComp = coords_set(0,0,0);

    mcDEBUG_INSTR()
  /* Component initializations. */
    /* Component Origin. */
  /* Setting parameters for component Origin. */
  SIG_MESSAGE("Origin (Init:SetPar)");
#line 39 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  if("NULL") strncpy(mccOrigin_profile, "NULL" ? "NULL" : "", 16384); else mccOrigin_profile[0]='\0';
#line 39 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccOrigin_percent = 10;
#line 39 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccOrigin_flag_save = 0;
#line 39 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccOrigin_minutes = 0;
#line 10989 "./Test_Monochromators.c"

  SIG_MESSAGE("Origin (Init:Place/Rotate)");
  rot_set_rotation(mcrotaOrigin,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 10996 "./Test_Monochromators.c"
  rot_copy(mcrotrOrigin, mcrotaOrigin);
  mcposaOrigin = coords_set(
#line 98 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
    0,
#line 98 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
    0,
#line 98 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
    0);
#line 11005 "./Test_Monochromators.c"
  mctc1 = coords_neg(mcposaOrigin);
  mcposrOrigin = rot_apply(mcrotaOrigin, mctc1);
  mcDEBUG_COMPONENT("Origin", mcposaOrigin, mcrotaOrigin)
  mccomp_posa[1] = mcposaOrigin;
  mccomp_posr[1] = mcposrOrigin;
  mcNCounter[1]  = mcPCounter[1] = mcP2Counter[1] = 0;
  mcAbsorbProp[1]= 0;
    /* Component Source. */
  /* Setting parameters for component Source. */
  SIG_MESSAGE("Source (Init:SetPar)");
#line 102 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccSource_radius = 0.01;
#line 52 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccSource_yheight = 0;
#line 52 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccSource_xwidth = 0;
#line 102 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccSource_dist = 1.0;
#line 103 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccSource_focus_xw = 0.02;
#line 103 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccSource_focus_yh = 0.02;
#line 54 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccSource_E0 = 0;
#line 54 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccSource_dE = 0;
#line 104 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccSource_lambda0 = mciplambda;
#line 104 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccSource_dlambda = 0.2 * mciplambda;
#line 55 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccSource_flux = 1;
#line 55 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccSource_gauss = 0;
#line 55 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccSource_target_index = + 1;
#line 11042 "./Test_Monochromators.c"

  SIG_MESSAGE("Source (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 11049 "./Test_Monochromators.c"
  rot_mul(mctr1, mcrotaOrigin, mcrotaSource);
  rot_transpose(mcrotaOrigin, mctr1);
  rot_mul(mcrotaSource, mctr1, mcrotrSource);
  mctc1 = coords_set(
#line 105 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
    0,
#line 105 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
    0,
#line 105 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
    0);
#line 11060 "./Test_Monochromators.c"
  rot_transpose(mcrotaOrigin, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaSource = coords_add(mcposaOrigin, mctc2);
  mctc1 = coords_sub(mcposaOrigin, mcposaSource);
  mcposrSource = rot_apply(mcrotaSource, mctc1);
  mcDEBUG_COMPONENT("Source", mcposaSource, mcrotaSource)
  mccomp_posa[2] = mcposaSource;
  mccomp_posr[2] = mcposrSource;
  mcNCounter[2]  = mcPCounter[2] = mcP2Counter[2] = 0;
  mcAbsorbProp[2]= 0;
    /* Component lamStart. */
  /* Setting parameters for component lamStart. */
  SIG_MESSAGE("lamStart (Init:SetPar)");
#line 108 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  if("lambdaStart.dat") strncpy(mcclamStart_filename, "lambdaStart.dat" ? "lambdaStart.dat" : "", 16384); else mcclamStart_filename[0]='\0';
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mcclamStart_xmin = -0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mcclamStart_xmax = 0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mcclamStart_ymin = -0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mcclamStart_ymax = 0.05;
#line 108 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mcclamStart_xwidth = 0.10;
#line 109 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mcclamStart_yheight = 0.10;
#line 109 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mcclamStart_Lmin = 0.7 * mciplambda;
#line 109 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mcclamStart_Lmax = 1.3 * mciplambda;
#line 51 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mcclamStart_restore_neutron = 0;
#line 51 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mcclamStart_nowritefile = 0;
#line 11096 "./Test_Monochromators.c"

  SIG_MESSAGE("lamStart (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 11103 "./Test_Monochromators.c"
  rot_mul(mctr1, mcrotaOrigin, mcrotalamStart);
  rot_transpose(mcrotaSource, mctr1);
  rot_mul(mcrotalamStart, mctr1, mcrotrlamStart);
  mctc1 = coords_set(
#line 110 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
    0,
#line 110 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
    0,
#line 110 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
    0.5);
#line 11114 "./Test_Monochromators.c"
  rot_transpose(mcrotaOrigin, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposalamStart = coords_add(mcposaOrigin, mctc2);
  mctc1 = coords_sub(mcposaSource, mcposalamStart);
  mcposrlamStart = rot_apply(mcrotalamStart, mctc1);
  mcDEBUG_COMPONENT("lamStart", mcposalamStart, mcrotalamStart)
  mccomp_posa[3] = mcposalamStart;
  mccomp_posr[3] = mcposrlamStart;
  mcNCounter[3]  = mcPCounter[3] = mcP2Counter[3] = 0;
  mcAbsorbProp[3]= 0;
    /* Component Mono_Arm. */
  /* Setting parameters for component Mono_Arm. */
  SIG_MESSAGE("Mono_Arm (Init:SetPar)");

  SIG_MESSAGE("Mono_Arm (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 114 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
    (0)*DEG2RAD,
#line 114 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
    (A1)*DEG2RAD,
#line 114 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
    (0)*DEG2RAD);
#line 11137 "./Test_Monochromators.c"
  rot_mul(mctr1, mcrotaSource, mcrotaMono_Arm);
  rot_transpose(mcrotalamStart, mctr1);
  rot_mul(mcrotaMono_Arm, mctr1, mcrotrMono_Arm);
  mctc1 = coords_set(
#line 114 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
    0,
#line 114 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
    0,
#line 114 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
    1.0);
#line 11148 "./Test_Monochromators.c"
  rot_transpose(mcrotaSource, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaMono_Arm = coords_add(mcposaSource, mctc2);
  mctc1 = coords_sub(mcposalamStart, mcposaMono_Arm);
  mcposrMono_Arm = rot_apply(mcrotaMono_Arm, mctc1);
  mcDEBUG_COMPONENT("Mono_Arm", mcposaMono_Arm, mcrotaMono_Arm)
  mccomp_posa[4] = mcposaMono_Arm;
  mccomp_posr[4] = mcposrMono_Arm;
  mcNCounter[4]  = mcPCounter[4] = mcP2Counter[4] = 0;
  mcAbsorbProp[4]= 0;
    /* Component Mono_Out. */
  /* Setting parameters for component Mono_Out. */
  SIG_MESSAGE("Mono_Out (Init:SetPar)");

  SIG_MESSAGE("Mono_Out (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 117 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
    (0)*DEG2RAD,
#line 117 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
    (A2)*DEG2RAD,
#line 117 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
    (0)*DEG2RAD);
#line 11171 "./Test_Monochromators.c"
  rot_mul(mctr1, mcrotaSource, mcrotaMono_Out);
  rot_transpose(mcrotaMono_Arm, mctr1);
  rot_mul(mcrotaMono_Out, mctr1, mcrotrMono_Out);
  mctc1 = coords_set(
#line 117 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
    0,
#line 117 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
    0,
#line 117 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
    0);
#line 11182 "./Test_Monochromators.c"
  rot_transpose(mcrotaMono_Arm, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaMono_Out = coords_add(mcposaMono_Arm, mctc2);
  mctc1 = coords_sub(mcposaMono_Arm, mcposaMono_Out);
  mcposrMono_Out = rot_apply(mcrotaMono_Out, mctc1);
  mcDEBUG_COMPONENT("Mono_Out", mcposaMono_Out, mcrotaMono_Out)
  mccomp_posa[5] = mcposaMono_Out;
  mccomp_posr[5] = mcposrMono_Out;
  mcNCounter[5]  = mcPCounter[5] = mcP2Counter[5] = 0;
  mcAbsorbProp[5]= 0;
    /* Component Mono1. */
  /* Setting parameters for component Mono1. */
  SIG_MESSAGE("Mono1 (Init:SetPar)");
#line 64 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccMono1_zmin = -0.05;
#line 64 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccMono1_zmax = 0.05;
#line 64 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccMono1_ymin = -0.05;
#line 64 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccMono1_ymax = 0.05;
#line 120 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccMono1_zwidth = 0.10;
#line 120 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccMono1_yheight = 0.10;
#line 121 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccMono1_mosaich = mcipMoz;
#line 121 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccMono1_mosaicv = mcipMoz;
#line 122 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccMono1_r0 = 0.9;
#line 122 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccMono1_Q = mono_q;
#line 66 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccMono1_DM = 0;
#line 11218 "./Test_Monochromators.c"

  SIG_MESSAGE("Mono1 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 11225 "./Test_Monochromators.c"
  rot_mul(mctr1, mcrotaMono_Arm, mcrotaMono1);
  rot_transpose(mcrotaMono_Out, mctr1);
  rot_mul(mcrotaMono1, mctr1, mcrotrMono1);
  mctc1 = coords_set(
#line 123 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
    0,
#line 123 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
    0,
#line 123 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
    0);
#line 11236 "./Test_Monochromators.c"
  rot_transpose(mcrotaMono_Arm, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaMono1 = coords_add(mcposaMono_Arm, mctc2);
  mctc1 = coords_sub(mcposaMono_Out, mcposaMono1);
  mcposrMono1 = rot_apply(mcrotaMono1, mctc1);
  mcDEBUG_COMPONENT("Mono1", mcposaMono1, mcrotaMono1)
  mccomp_posa[6] = mcposaMono1;
  mccomp_posr[6] = mcposrMono1;
  mcNCounter[6]  = mcPCounter[6] = mcP2Counter[6] = 0;
  mcAbsorbProp[6]= 0;
    /* Component Mono2. */
  /* Setting parameters for component Mono2. */
  SIG_MESSAGE("Mono2 (Init:SetPar)");
#line 126 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccMono2_zwidth = 0.10;
#line 126 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccMono2_yheight = 0.10;
#line 127 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccMono2_mosaic = mcipMoz;
#line 127 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccMono2_dspread = 0.0;
#line 128 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccMono2_Q = mono_q;
#line 89 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccMono2_DM = 0;
#line 89 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccMono2_pThreshold = 0;
#line 128 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccMono2_Rup = 0.9;
#line 128 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccMono2_Rdown = 0.9;
#line 89 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccMono2_debug = 0;
#line 11270 "./Test_Monochromators.c"

  SIG_MESSAGE("Mono2 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 11277 "./Test_Monochromators.c"
  rot_mul(mctr1, mcrotaMono_Arm, mcrotaMono2);
  rot_transpose(mcrotaMono1, mctr1);
  rot_mul(mcrotaMono2, mctr1, mcrotrMono2);
  mctc1 = coords_set(
#line 129 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
    0,
#line 129 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
    0,
#line 129 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
    0);
#line 11288 "./Test_Monochromators.c"
  rot_transpose(mcrotaMono_Arm, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaMono2 = coords_add(mcposaMono_Arm, mctc2);
  mctc1 = coords_sub(mcposaMono1, mcposaMono2);
  mcposrMono2 = rot_apply(mcrotaMono2, mctc1);
  mcDEBUG_COMPONENT("Mono2", mcposaMono2, mcrotaMono2)
  mccomp_posa[7] = mcposaMono2;
  mccomp_posr[7] = mcposrMono2;
  mcNCounter[7]  = mcPCounter[7] = mcP2Counter[7] = 0;
  mcAbsorbProp[7]= 0;
    /* Component Mono3. */
  /* Setting parameters for component Mono3. */
  SIG_MESSAGE("Mono3 (Init:SetPar)");
#line 132 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccMono3_zwidth = 0.10;
#line 132 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccMono3_yheight = 0.10;
#line 133 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccMono3_mosaic = mcipMoz;
#line 133 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccMono3_dspread = 0.0;
#line 134 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccMono3_Q = mono_q;
#line 89 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccMono3_DM = 0;
#line 134 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccMono3_pThreshold = 0.01;
#line 134 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccMono3_Rup = 0.9;
#line 134 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccMono3_Rdown = 0.9;
#line 89 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccMono3_debug = 0;
#line 11322 "./Test_Monochromators.c"

  SIG_MESSAGE("Mono3 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 11329 "./Test_Monochromators.c"
  rot_mul(mctr1, mcrotaMono_Arm, mcrotaMono3);
  rot_transpose(mcrotaMono2, mctr1);
  rot_mul(mcrotaMono3, mctr1, mcrotrMono3);
  mctc1 = coords_set(
#line 135 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
    0,
#line 135 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
    0,
#line 135 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
    0);
#line 11340 "./Test_Monochromators.c"
  rot_transpose(mcrotaMono_Arm, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaMono3 = coords_add(mcposaMono_Arm, mctc2);
  mctc1 = coords_sub(mcposaMono2, mcposaMono3);
  mcposrMono3 = rot_apply(mcrotaMono3, mctc1);
  mcDEBUG_COMPONENT("Mono3", mcposaMono3, mcrotaMono3)
  mccomp_posa[8] = mcposaMono3;
  mccomp_posr[8] = mcposrMono3;
  mcNCounter[8]  = mcPCounter[8] = mcP2Counter[8] = 0;
  mcAbsorbProp[8]= 0;
    /* Component Mono4. */
  /* Setting parameters for component Mono4. */
  SIG_MESSAGE("Mono4 (Init:SetPar)");
#line 99 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  if("NULL") strncpy(mccMono4_reflect, "NULL" ? "NULL" : "", 16384); else mccMono4_reflect[0]='\0';
#line 99 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  if("NULL") strncpy(mccMono4_transmit, "NULL" ? "NULL" : "", 16384); else mccMono4_transmit[0]='\0';
#line 100 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccMono4_zwidth = 0.01;
#line 100 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccMono4_yheight = 0.01;
#line 140 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccMono4_gap = 0;
#line 101 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccMono4_NH = 11;
#line 101 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccMono4_NV = 11;
#line 139 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccMono4_mosaich = mcipMoz;
#line 139 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccMono4_mosaicv = mcipMoz;
#line 140 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccMono4_r0 = 0.9;
#line 101 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccMono4_t0 = 1.0;
#line 140 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccMono4_Q = mono_q;
#line 102 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccMono4_RV = 0;
#line 102 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccMono4_RH = 0;
#line 102 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccMono4_DM = 0;
#line 102 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccMono4_mosaic = 0;
#line 138 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccMono4_width = 0.10;
#line 138 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccMono4_height = 0.10;
#line 102 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccMono4_verbose = 0;
#line 102 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccMono4_order = 0;
#line 11394 "./Test_Monochromators.c"

  SIG_MESSAGE("Mono4 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 11401 "./Test_Monochromators.c"
  rot_mul(mctr1, mcrotaMono_Arm, mcrotaMono4);
  rot_transpose(mcrotaMono3, mctr1);
  rot_mul(mcrotaMono4, mctr1, mcrotrMono4);
  mctc1 = coords_set(
#line 141 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
    0,
#line 141 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
    0,
#line 141 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
    0);
#line 11412 "./Test_Monochromators.c"
  rot_transpose(mcrotaMono_Arm, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaMono4 = coords_add(mcposaMono_Arm, mctc2);
  mctc1 = coords_sub(mcposaMono3, mcposaMono4);
  mcposrMono4 = rot_apply(mcrotaMono4, mctc1);
  mcDEBUG_COMPONENT("Mono4", mcposaMono4, mcrotaMono4)
  mccomp_posa[9] = mcposaMono4;
  mccomp_posr[9] = mcposrMono4;
  mcNCounter[9]  = mcPCounter[9] = mcP2Counter[9] = 0;
  mcAbsorbProp[9]= 0;
    /* Component Mono5. */
  /* Setting parameters for component Mono5. */
  SIG_MESSAGE("Mono5 (Init:SetPar)");
#line 78 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  if(0) strncpy(mccMono5_reflect, 0 ? 0 : "", 16384); else mccMono5_reflect[0]='\0';
#line 78 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccMono5_zwidth = 0.01;
#line 78 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccMono5_yheight = 0.01;
#line 146 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccMono5_gap = 0;
#line 78 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccMono5_NH = 11;
#line 78 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccMono5_NV = 11;
#line 145 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccMono5_mosaich = mcipMoz;
#line 145 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccMono5_mosaicv = mcipMoz;
#line 146 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccMono5_r0 = 0.9;
#line 146 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccMono5_Q = mono_q;
#line 78 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccMono5_RV = 0;
#line 78 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccMono5_RH = 0;
#line 78 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccMono5_DM = 0;
#line 78 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccMono5_mosaic = 0;
#line 144 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccMono5_width = 0.10;
#line 144 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccMono5_height = 0.10;
#line 78 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccMono5_verbose = 0;
#line 11460 "./Test_Monochromators.c"

  SIG_MESSAGE("Mono5 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 11467 "./Test_Monochromators.c"
  rot_mul(mctr1, mcrotaMono_Arm, mcrotaMono5);
  rot_transpose(mcrotaMono4, mctr1);
  rot_mul(mcrotaMono5, mctr1, mcrotrMono5);
  mctc1 = coords_set(
#line 147 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
    0,
#line 147 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
    0,
#line 147 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
    0);
#line 11478 "./Test_Monochromators.c"
  rot_transpose(mcrotaMono_Arm, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaMono5 = coords_add(mcposaMono_Arm, mctc2);
  mctc1 = coords_sub(mcposaMono4, mcposaMono5);
  mcposrMono5 = rot_apply(mcrotaMono5, mctc1);
  mcDEBUG_COMPONENT("Mono5", mcposaMono5, mcrotaMono5)
  mccomp_posa[10] = mcposaMono5;
  mccomp_posr[10] = mcposrMono5;
  mcNCounter[10]  = mcPCounter[10] = mcP2Counter[10] = 0;
  mcAbsorbProp[10]= 0;
    /* Component Mono6. */
  /* Setting parameters for component Mono6. */
  SIG_MESSAGE("Mono6 (Init:SetPar)");
#line 151 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  if("C_graphite.lau") strncpy(mccMono6_reflections, "C_graphite.lau" ? "C_graphite.lau" : "", 16384); else mccMono6_reflections[0]='\0';
#line 282 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  if(0) strncpy(mccMono6_geometry, 0 ? 0 : "", 16384); else mccMono6_geometry[0]='\0';
#line 150 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccMono6_xwidth = 0.002;
#line 150 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccMono6_yheight = 0.1;
#line 150 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccMono6_zdepth = 0.1;
#line 283 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccMono6_radius = 0;
#line 151 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccMono6_delta_d_d = 1e-2;
#line 151 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccMono6_mosaic = mcipMoz;
#line 284 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccMono6_mosaic_a = -1;
#line 284 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccMono6_mosaic_b = -1;
#line 284 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccMono6_mosaic_c = -1;
#line 285 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccMono6_recip_cell = 0;
#line 151 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccMono6_barns = 0;
#line 152 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccMono6_ax = 0;
#line 152 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccMono6_ay = 2.14;
#line 152 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccMono6_az = -1.24;
#line 153 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccMono6_bx = 0;
#line 153 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccMono6_by = 0;
#line 153 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccMono6_bz = 2.47;
#line 154 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccMono6_cx = 6.71;
#line 154 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccMono6_cy = 0;
#line 154 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccMono6_cz = 0;
#line 289 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccMono6_p_transmit = 0.001;
#line 155 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccMono6_sigma_abs = 0.014;
#line 155 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccMono6_sigma_inc = 0.004;
#line 290 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccMono6_aa = 0;
#line 290 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccMono6_bb = 0;
#line 290 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccMono6_cc = 0;
#line 155 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccMono6_order = mciporder;
#line 290 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccMono6_RX = 0;
#line 290 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccMono6_RY = 0;
#line 155 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccMono6_powder = mcippowder;
#line 155 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccMono6_PG = mcipPG;
#line 291 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccMono6_deltak = 1e-6;
#line 11560 "./Test_Monochromators.c"

  SIG_MESSAGE("Mono6 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 11567 "./Test_Monochromators.c"
  rot_mul(mctr1, mcrotaMono_Arm, mcrotaMono6);
  rot_transpose(mcrotaMono5, mctr1);
  rot_mul(mcrotaMono6, mctr1, mcrotrMono6);
  mctc1 = coords_set(
#line 156 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
    0,
#line 156 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
    0,
#line 156 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
    0);
#line 11578 "./Test_Monochromators.c"
  rot_transpose(mcrotaMono_Arm, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaMono6 = coords_add(mcposaMono_Arm, mctc2);
  mctc1 = coords_sub(mcposaMono5, mcposaMono6);
  mcposrMono6 = rot_apply(mcrotaMono6, mctc1);
  mcDEBUG_COMPONENT("Mono6", mcposaMono6, mcrotaMono6)
  mccomp_posa[11] = mcposaMono6;
  mccomp_posr[11] = mcposrMono6;
  mcNCounter[11]  = mcPCounter[11] = mcP2Counter[11] = 0;
  mcAbsorbProp[11]= 0;
    /* Component Sphere1. */
  /* Setting parameters for component Sphere1. */
  SIG_MESSAGE("Sphere1 (Init:SetPar)");
#line 161 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  if("sphere") strncpy(mccSphere1_filename, "sphere" ? "sphere" : "", 16384); else mccSphere1_filename[0]='\0';
#line 161 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccSphere1_radius = 1;
#line 161 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccSphere1_restore_neutron = 1;
#line 49 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccSphere1_nowritefile = 0;
#line 11600 "./Test_Monochromators.c"

  SIG_MESSAGE("Sphere1 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 11607 "./Test_Monochromators.c"
  rot_mul(mctr1, mcrotaMono_Out, mcrotaSphere1);
  rot_transpose(mcrotaMono6, mctr1);
  rot_mul(mcrotaSphere1, mctr1, mcrotrSphere1);
  mctc1 = coords_set(
#line 162 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
    0,
#line 162 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
    0,
#line 162 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
    0);
#line 11618 "./Test_Monochromators.c"
  rot_transpose(mcrotaMono_Out, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaSphere1 = coords_add(mcposaMono_Out, mctc2);
  mctc1 = coords_sub(mcposaMono6, mcposaSphere1);
  mcposrSphere1 = rot_apply(mcrotaSphere1, mctc1);
  mcDEBUG_COMPONENT("Sphere1", mcposaSphere1, mcrotaSphere1)
  mccomp_posa[12] = mcposaSphere1;
  mccomp_posr[12] = mcposrSphere1;
  mcNCounter[12]  = mcPCounter[12] = mcP2Counter[12] = 0;
  mcAbsorbProp[12]= 0;
    /* Component lam1. */
  /* Setting parameters for component lam1. */
  SIG_MESSAGE("lam1 (Init:SetPar)");
#line 165 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  if("lambda1.dat") strncpy(mcclam1_filename, "lambda1.dat" ? "lambda1.dat" : "", 16384); else mcclam1_filename[0]='\0';
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mcclam1_xmin = -0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mcclam1_xmax = 0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mcclam1_ymin = -0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mcclam1_ymax = 0.05;
#line 165 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mcclam1_xwidth = 0.10;
#line 166 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mcclam1_yheight = 0.10;
#line 166 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mcclam1_Lmin = 0.7 * mciplambda;
#line 166 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mcclam1_Lmax = 1.3 * mciplambda;
#line 51 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mcclam1_restore_neutron = 0;
#line 51 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mcclam1_nowritefile = 0;
#line 11654 "./Test_Monochromators.c"

  SIG_MESSAGE("lam1 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 11661 "./Test_Monochromators.c"
  rot_mul(mctr1, mcrotaMono_Out, mcrotalam1);
  rot_transpose(mcrotaSphere1, mctr1);
  rot_mul(mcrotalam1, mctr1, mcrotrlam1);
  mctc1 = coords_set(
#line 167 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
    0,
#line 167 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
    0,
#line 167 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
    0.25);
#line 11672 "./Test_Monochromators.c"
  rot_transpose(mcrotaMono_Out, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposalam1 = coords_add(mcposaMono_Out, mctc2);
  mctc1 = coords_sub(mcposaSphere1, mcposalam1);
  mcposrlam1 = rot_apply(mcrotalam1, mctc1);
  mcDEBUG_COMPONENT("lam1", mcposalam1, mcrotalam1)
  mccomp_posa[13] = mcposalam1;
  mccomp_posr[13] = mcposrlam1;
  mcNCounter[13]  = mcPCounter[13] = mcP2Counter[13] = 0;
  mcAbsorbProp[13]= 0;
    /* Component psd1. */
  /* Setting parameters for component psd1. */
  SIG_MESSAGE("psd1 (Init:SetPar)");
#line 172 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccpsd1_nx = 20;
#line 172 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccpsd1_ny = 20;
#line 172 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  if("psd1.dat") strncpy(mccpsd1_filename, "psd1.dat" ? "psd1.dat" : "", 16384); else mccpsd1_filename[0]='\0';
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccpsd1_xmin = -0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccpsd1_xmax = 0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccpsd1_ymin = -0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccpsd1_ymax = 0.05;
#line 171 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccpsd1_xwidth = 0.10;
#line 171 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccpsd1_yheight = 0.10;
#line 51 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  mccpsd1_restore_neutron = 0;
#line 11706 "./Test_Monochromators.c"

  SIG_MESSAGE("psd1 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 11713 "./Test_Monochromators.c"
  rot_mul(mctr1, mcrotaMono_Out, mcrotapsd1);
  rot_transpose(mcrotalam1, mctr1);
  rot_mul(mcrotapsd1, mctr1, mcrotrpsd1);
  mctc1 = coords_set(
#line 173 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
    0,
#line 173 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
    0,
#line 173 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
    0.5);
#line 11724 "./Test_Monochromators.c"
  rot_transpose(mcrotaMono_Out, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposapsd1 = coords_add(mcposaMono_Out, mctc2);
  mctc1 = coords_sub(mcposalam1, mcposapsd1);
  mcposrpsd1 = rot_apply(mcrotapsd1, mctc1);
  mcDEBUG_COMPONENT("psd1", mcposapsd1, mcrotapsd1)
  mccomp_posa[14] = mcposapsd1;
  mccomp_posr[14] = mcposrpsd1;
  mcNCounter[14]  = mcPCounter[14] = mcP2Counter[14] = 0;
  mcAbsorbProp[14]= 0;
  /* Component initializations. */
  /* Initializations for component Origin. */
  SIG_MESSAGE("Origin (Init)");
#define mccompcurname  Origin
#define mccompcurtype  Progress_bar
#define mccompcurindex 1
#define IntermediateCnts mccOrigin_IntermediateCnts
#define StartTime mccOrigin_StartTime
#define EndTime mccOrigin_EndTime
#define CurrentTime mccOrigin_CurrentTime
#define profile mccOrigin_profile
#define percent mccOrigin_percent
#define flag_save mccOrigin_flag_save
#define minutes mccOrigin_minutes
#line 57 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../misc/Progress_bar.comp"
{
IntermediateCnts=0;
StartTime=0;
EndTime=0;
CurrentTime=0;

fprintf(stdout, "[%s] Initialize\n", mcinstrument_name);
  if (percent*mcget_ncount()/100 < 1e5) {
    percent=1e5*100.0/mcget_ncount();
  }
}
#line 11761 "./Test_Monochromators.c"
#undef minutes
#undef flag_save
#undef percent
#undef profile
#undef CurrentTime
#undef EndTime
#undef StartTime
#undef IntermediateCnts
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component Source. */
  SIG_MESSAGE("Source (Init)");
#define mccompcurname  Source
#define mccompcurtype  Source_simple
#define mccompcurindex 2
#define pmul mccSource_pmul
#define square mccSource_square
#define srcArea mccSource_srcArea
#define radius mccSource_radius
#define yheight mccSource_yheight
#define xwidth mccSource_xwidth
#define dist mccSource_dist
#define focus_xw mccSource_focus_xw
#define focus_yh mccSource_focus_yh
#define E0 mccSource_E0
#define dE mccSource_dE
#define lambda0 mccSource_lambda0
#define dlambda mccSource_dlambda
#define flux mccSource_flux
#define gauss mccSource_gauss
#define target_index mccSource_target_index
#line 65 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../sources/Source_simple.comp"
{
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
    ToTarget = coords_sub(POS_A_COMP_INDEX(INDEX_CURRENT_COMP+target_index),POS_A_CURRENT_COMP);
    ToTarget = rot_apply(ROT_A_CURRENT_COMP, ToTarget);
    coords_get(ToTarget, &tx, &ty, &tz);
    dist=sqrt(tx*tx+ty*ty+tz*tz);
  } else if (dist) {
    tx = 0;
    ty = 0;
    tz = dist;
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
}
#line 11855 "./Test_Monochromators.c"
#undef target_index
#undef gauss
#undef flux
#undef dlambda
#undef lambda0
#undef dE
#undef E0
#undef focus_yh
#undef focus_xw
#undef dist
#undef xwidth
#undef yheight
#undef radius
#undef srcArea
#undef square
#undef pmul
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component lamStart. */
  SIG_MESSAGE("lamStart (Init)");
#define mccompcurname  lamStart
#define mccompcurtype  L_monitor
#define mccompcurindex 3
#define nL mcclamStart_nL
#define L_N mcclamStart_L_N
#define L_p mcclamStart_L_p
#define L_p2 mcclamStart_L_p2
#define filename mcclamStart_filename
#define xmin mcclamStart_xmin
#define xmax mcclamStart_xmax
#define ymin mcclamStart_ymin
#define ymax mcclamStart_ymax
#define xwidth mcclamStart_xwidth
#define yheight mcclamStart_yheight
#define Lmin mcclamStart_Lmin
#define Lmax mcclamStart_Lmax
#define restore_neutron mcclamStart_restore_neutron
#define nowritefile mcclamStart_nowritefile
#line 62 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/L_monitor.comp"
{
int i;

if (xwidth  > 0) { xmax = xwidth/2;  xmin = -xmax; }
    if (yheight > 0) { ymax = yheight/2; ymin = -ymax; }

    if ((xmin >= xmax) || (ymin >= ymax)) {
            printf("L_monitor: %s: Null detection area !\n"
                   "ERROR      (xwidth,yheight,xmin,xmax,ymin,ymax). Exiting",
           NAME_CURRENT_COMP);
      exit(0);
    }

    for (i=0; i<nL; i++)
    {
      L_N[i] = 0;
      L_p[i] = 0;
      L_p2[i] = 0;
    }
}
#line 11917 "./Test_Monochromators.c"
#undef nowritefile
#undef restore_neutron
#undef Lmax
#undef Lmin
#undef yheight
#undef xwidth
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef filename
#undef L_p2
#undef L_p
#undef L_N
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component Mono_Arm. */
  SIG_MESSAGE("Mono_Arm (Init)");

  /* Initializations for component Mono_Out. */
  SIG_MESSAGE("Mono_Out (Init)");

  /* Initializations for component Mono1. */
  SIG_MESSAGE("Mono1 (Init)");
#define mccompcurname  Mono1
#define mccompcurtype  Monochromator_flat
#define mccompcurindex 6
#define mos_rms_y mccMono1_mos_rms_y
#define mos_rms_z mccMono1_mos_rms_z
#define mos_rms_max mccMono1_mos_rms_max
#define mono_Q mccMono1_mono_Q
#define zmin mccMono1_zmin
#define zmax mccMono1_zmax
#define ymin mccMono1_ymin
#define ymax mccMono1_ymax
#define zwidth mccMono1_zwidth
#define yheight mccMono1_yheight
#define mosaich mccMono1_mosaich
#define mosaicv mccMono1_mosaicv
#define r0 mccMono1_r0
#define Q mccMono1_Q
#define DM mccMono1_DM
#line 102 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Monochromator_flat.comp"
{
  mos_rms_y = MIN2RAD*mosaicv/sqrt(8*log(2));
  mos_rms_z = MIN2RAD*mosaich/sqrt(8*log(2));
  mos_rms_max = mos_rms_y > mos_rms_z ? mos_rms_y : mos_rms_z;

  mono_Q = Q;
  if (DM != 0) mono_Q = 2*PI/DM;

  if (zwidth>0)  { zmax = zwidth/2;  zmin=-zmax; }
  if (yheight>0) { ymax = yheight/2; ymin=-ymax; }

  if (zmin==zmax || ymin==ymax)
    exit(fprintf(stderr, "Monochromator_flat: %s : Surface is null (zmin,zmax,ymin,ymax)\n", NAME_CURRENT_COMP));
}
#line 11978 "./Test_Monochromators.c"
#undef DM
#undef Q
#undef r0
#undef mosaicv
#undef mosaich
#undef yheight
#undef zwidth
#undef ymax
#undef ymin
#undef zmax
#undef zmin
#undef mono_Q
#undef mos_rms_max
#undef mos_rms_z
#undef mos_rms_y
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component Mono2. */
  SIG_MESSAGE("Mono2 (Init)");
#define mccompcurname  Mono2
#define mccompcurtype  Monochromator_pol
#define mccompcurindex 7
#define mos_rms mccMono2_mos_rms
#define d_rms mccMono2_d_rms
#define mono_Q mccMono2_mono_Q
#define FN mccMono2_FN
#define FM mccMono2_FM
#define zwidth mccMono2_zwidth
#define yheight mccMono2_yheight
#define mosaic mccMono2_mosaic
#define dspread mccMono2_dspread
#define Q mccMono2_Q
#define DM mccMono2_DM
#define pThreshold mccMono2_pThreshold
#define Rup mccMono2_Rup
#define Rdown mccMono2_Rdown
#define debug mccMono2_debug
#line 109 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Monochromator_pol.comp"
{
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
}
#line 12038 "./Test_Monochromators.c"
#undef debug
#undef Rdown
#undef Rup
#undef pThreshold
#undef DM
#undef Q
#undef dspread
#undef mosaic
#undef yheight
#undef zwidth
#undef FM
#undef FN
#undef mono_Q
#undef d_rms
#undef mos_rms
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component Mono3. */
  SIG_MESSAGE("Mono3 (Init)");
#define mccompcurname  Mono3
#define mccompcurtype  Monochromator_pol
#define mccompcurindex 8
#define mos_rms mccMono3_mos_rms
#define d_rms mccMono3_d_rms
#define mono_Q mccMono3_mono_Q
#define FN mccMono3_FN
#define FM mccMono3_FM
#define zwidth mccMono3_zwidth
#define yheight mccMono3_yheight
#define mosaic mccMono3_mosaic
#define dspread mccMono3_dspread
#define Q mccMono3_Q
#define DM mccMono3_DM
#define pThreshold mccMono3_pThreshold
#define Rup mccMono3_Rup
#define Rdown mccMono3_Rdown
#define debug mccMono3_debug
#line 109 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Monochromator_pol.comp"
{
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
}
#line 12098 "./Test_Monochromators.c"
#undef debug
#undef Rdown
#undef Rup
#undef pThreshold
#undef DM
#undef Q
#undef dspread
#undef mosaic
#undef yheight
#undef zwidth
#undef FM
#undef FN
#undef mono_Q
#undef d_rms
#undef mos_rms
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component Mono4. */
  SIG_MESSAGE("Mono4 (Init)");
#define mccompcurname  Mono4
#define mccompcurtype  Monochromator_curved
#define mccompcurindex 9
#define mos_rms_y mccMono4_mos_rms_y
#define mos_rms_z mccMono4_mos_rms_z
#define mos_rms_max mccMono4_mos_rms_max
#define mono_Q mccMono4_mono_Q
#define SlabWidth mccMono4_SlabWidth
#define SlabHeight mccMono4_SlabHeight
#define rTable mccMono4_rTable
#define tTable mccMono4_tTable
#define row mccMono4_row
#define col mccMono4_col
#define tiltH mccMono4_tiltH
#define tiltV mccMono4_tiltV
#define reflect mccMono4_reflect
#define transmit mccMono4_transmit
#define zwidth mccMono4_zwidth
#define yheight mccMono4_yheight
#define gap mccMono4_gap
#define NH mccMono4_NH
#define NV mccMono4_NV
#define mosaich mccMono4_mosaich
#define mosaicv mccMono4_mosaicv
#define r0 mccMono4_r0
#define t0 mccMono4_t0
#define Q mccMono4_Q
#define RV mccMono4_RV
#define RH mccMono4_RH
#define DM mccMono4_DM
#define mosaic mccMono4_mosaic
#define width mccMono4_width
#define height mccMono4_height
#define verbose mccMono4_verbose
#define order mccMono4_order
#line 148 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Monochromator_curved.comp"
{
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

}
#line 12236 "./Test_Monochromators.c"
#undef order
#undef verbose
#undef height
#undef width
#undef mosaic
#undef DM
#undef RH
#undef RV
#undef Q
#undef t0
#undef r0
#undef mosaicv
#undef mosaich
#undef NV
#undef NH
#undef gap
#undef yheight
#undef zwidth
#undef transmit
#undef reflect
#undef tiltV
#undef tiltH
#undef col
#undef row
#undef tTable
#undef rTable
#undef SlabHeight
#undef SlabWidth
#undef mono_Q
#undef mos_rms_max
#undef mos_rms_z
#undef mos_rms_y
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component Mono5. */
  SIG_MESSAGE("Mono5 (Init)");
#define mccompcurname  Mono5
#define mccompcurtype  Monochromator_2foc
#define mccompcurindex 10
#define mos_y mccMono5_mos_y
#define mos_z mccMono5_mos_z
#define mono_Q mccMono5_mono_Q
#define SlabWidth mccMono5_SlabWidth
#define SlabHeight mccMono5_SlabHeight
#define rTable mccMono5_rTable
#define reflect mccMono5_reflect
#define zwidth mccMono5_zwidth
#define yheight mccMono5_yheight
#define gap mccMono5_gap
#define NH mccMono5_NH
#define NV mccMono5_NV
#define mosaich mccMono5_mosaich
#define mosaicv mccMono5_mosaicv
#define r0 mccMono5_r0
#define Q mccMono5_Q
#define RV mccMono5_RV
#define RH mccMono5_RH
#define DM mccMono5_DM
#define mosaic mccMono5_mosaic
#define width mccMono5_width
#define height mccMono5_height
#define verbose mccMono5_verbose
#line 99 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../contrib/Monochromator_2foc.comp"
{
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
}
#line 12349 "./Test_Monochromators.c"
#undef verbose
#undef height
#undef width
#undef mosaic
#undef DM
#undef RH
#undef RV
#undef Q
#undef r0
#undef mosaicv
#undef mosaich
#undef NV
#undef NH
#undef gap
#undef yheight
#undef zwidth
#undef reflect
#undef rTable
#undef SlabHeight
#undef SlabWidth
#undef mono_Q
#undef mos_z
#undef mos_y
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component Mono6. */
  SIG_MESSAGE("Mono6 (Init)");
#define mccompcurname  Mono6
#define mccompcurtype  Single_crystal
#define mccompcurindex 11
#define mosaic_AB mccMono6_mosaic_AB
#define hkl_info mccMono6_hkl_info
#define offdata mccMono6_offdata
#define reflections mccMono6_reflections
#define geometry mccMono6_geometry
#define xwidth mccMono6_xwidth
#define yheight mccMono6_yheight
#define zdepth mccMono6_zdepth
#define radius mccMono6_radius
#define delta_d_d mccMono6_delta_d_d
#define mosaic mccMono6_mosaic
#define mosaic_a mccMono6_mosaic_a
#define mosaic_b mccMono6_mosaic_b
#define mosaic_c mccMono6_mosaic_c
#define recip_cell mccMono6_recip_cell
#define barns mccMono6_barns
#define ax mccMono6_ax
#define ay mccMono6_ay
#define az mccMono6_az
#define bx mccMono6_bx
#define by mccMono6_by
#define bz mccMono6_bz
#define cx mccMono6_cx
#define cy mccMono6_cy
#define cz mccMono6_cz
#define p_transmit mccMono6_p_transmit
#define sigma_abs mccMono6_sigma_abs
#define sigma_inc mccMono6_sigma_inc
#define aa mccMono6_aa
#define bb mccMono6_bb
#define cc mccMono6_cc
#define order mccMono6_order
#define RX mccMono6_RX
#define RY mccMono6_RY
#define powder mccMono6_powder
#define PG mccMono6_PG
#define deltak mccMono6_deltak
#line 975 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../samples/Single_crystal.comp"
{
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

  /*this is necessary to allow a numerical array to be passed through as a DEFINITION parameter*/
  double mosaic_ABin[]=mosaic_AB;
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

  if (powder && PG)
    exit(fprintf(stderr,"Single_crystal: %s: powder and PG modes can not be used together!\n"
             "ERROR           Please use EITHER powder or PG mode.\n", NAME_CURRENT_COMP));

  if (powder && !(order==1)) {
    fprintf(stderr,"Single_crystal: %s: powder mode means implicit choice of no multiple scattering!\n"
            "WARNING setting order=1\n", NAME_CURRENT_COMP);
    order=1;
  }

  if (PG && !(order==1)) {
    fprintf(stderr,"Single_crystal: %s: PG mode means implicit choice of no multiple scattering!\n"
            "WARNING setting order=1\n", NAME_CURRENT_COMP);
    order=1;
  }


}
#line 12508 "./Test_Monochromators.c"
#undef deltak
#undef PG
#undef powder
#undef RY
#undef RX
#undef order
#undef cc
#undef bb
#undef aa
#undef sigma_inc
#undef sigma_abs
#undef p_transmit
#undef cz
#undef cy
#undef cx
#undef bz
#undef by
#undef bx
#undef az
#undef ay
#undef ax
#undef barns
#undef recip_cell
#undef mosaic_c
#undef mosaic_b
#undef mosaic_a
#undef mosaic
#undef delta_d_d
#undef radius
#undef zdepth
#undef yheight
#undef xwidth
#undef geometry
#undef reflections
#undef offdata
#undef hkl_info
#undef mosaic_AB
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component Sphere1. */
  SIG_MESSAGE("Sphere1 (Init)");
#define mccompcurname  Sphere1
#define mccompcurtype  PSD_monitor_4PI
#define mccompcurindex 12
#define nx mccSphere1_nx
#define ny mccSphere1_ny
#define PSD_N mccSphere1_PSD_N
#define PSD_p mccSphere1_PSD_p
#define PSD_p2 mccSphere1_PSD_p2
#define filename mccSphere1_filename
#define radius mccSphere1_radius
#define restore_neutron mccSphere1_restore_neutron
#define nowritefile mccSphere1_nowritefile
#line 61 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor_4PI.comp"
{
int i,j;

for (i=0; i<nx; i++)
    for (j=0; j<ny; j++)
    {
      PSD_N[i][j] = 0;
      PSD_p[i][j] = 0;
      PSD_p2[i][j] = 0;
    }
}
#line 12576 "./Test_Monochromators.c"
#undef nowritefile
#undef restore_neutron
#undef radius
#undef filename
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef ny
#undef nx
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component lam1. */
  SIG_MESSAGE("lam1 (Init)");
#define mccompcurname  lam1
#define mccompcurtype  L_monitor
#define mccompcurindex 13
#define nL mcclam1_nL
#define L_N mcclam1_L_N
#define L_p mcclam1_L_p
#define L_p2 mcclam1_L_p2
#define filename mcclam1_filename
#define xmin mcclam1_xmin
#define xmax mcclam1_xmax
#define ymin mcclam1_ymin
#define ymax mcclam1_ymax
#define xwidth mcclam1_xwidth
#define yheight mcclam1_yheight
#define Lmin mcclam1_Lmin
#define Lmax mcclam1_Lmax
#define restore_neutron mcclam1_restore_neutron
#define nowritefile mcclam1_nowritefile
#line 62 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/L_monitor.comp"
{
int i;

if (xwidth  > 0) { xmax = xwidth/2;  xmin = -xmax; }
    if (yheight > 0) { ymax = yheight/2; ymin = -ymax; }

    if ((xmin >= xmax) || (ymin >= ymax)) {
            printf("L_monitor: %s: Null detection area !\n"
                   "ERROR      (xwidth,yheight,xmin,xmax,ymin,ymax). Exiting",
           NAME_CURRENT_COMP);
      exit(0);
    }

    for (i=0; i<nL; i++)
    {
      L_N[i] = 0;
      L_p[i] = 0;
      L_p2[i] = 0;
    }
}
#line 12631 "./Test_Monochromators.c"
#undef nowritefile
#undef restore_neutron
#undef Lmax
#undef Lmin
#undef yheight
#undef xwidth
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef filename
#undef L_p2
#undef L_p
#undef L_N
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component psd1. */
  SIG_MESSAGE("psd1 (Init)");
#define mccompcurname  psd1
#define mccompcurtype  PSD_monitor
#define mccompcurindex 14
#define PSD_N mccpsd1_PSD_N
#define PSD_p mccpsd1_PSD_p
#define PSD_p2 mccpsd1_PSD_p2
#define nx mccpsd1_nx
#define ny mccpsd1_ny
#define filename mccpsd1_filename
#define xmin mccpsd1_xmin
#define xmax mccpsd1_xmax
#define ymin mccpsd1_ymin
#define ymax mccpsd1_ymax
#define xwidth mccpsd1_xwidth
#define yheight mccpsd1_yheight
#define restore_neutron mccpsd1_restore_neutron
#line 68 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
{
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
}
#line 12694 "./Test_Monochromators.c"
#undef restore_neutron
#undef yheight
#undef xwidth
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef filename
#undef ny
#undef nx
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if(mcdotrace) mcdisplay();
    mcDEBUG_INSTR_END()
  }

} /* end init */

void mcraytrace(void) {
  /* Neutronics-specific defines */
#ifdef NEUTRONICS
extern double mcnx, mcny, mcnz, mcnvx, mcnvy, mcnvz;
extern double mcnt, mcnsx, mcnsy, mcnsz, mcnp;
#endif
  /* End of Neutronics-specific defines */
  /* Copy neutron state to local variables. */
  MCNUM mcnlx = mcnx;
  MCNUM mcnly = mcny;
  MCNUM mcnlz = mcnz;
  MCNUM mcnlvx = mcnvx;
  MCNUM mcnlvy = mcnvy;
  MCNUM mcnlvz = mcnvz;
  MCNUM mcnlt = mcnt;
  MCNUM mcnlsx = mcnsx;
  MCNUM mcnlsy = mcnsy;
  MCNUM mcnlsz = mcnsz;
  MCNUM mcnlp = mcnp;

  mcDEBUG_ENTER()
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define mcabsorb mcabsorbAll
  /* TRACE Component Origin [1] */
  mccoordschange(mcposrOrigin, mcrotrOrigin,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component Origin (without coords transformations) */
  mcJumpTrace_Origin:
  SIG_MESSAGE("Origin (Trace)");
  mcDEBUG_COMP("Origin")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompOrigin
  STORE_NEUTRON(1,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[1]++;
  mcPCounter[1] += p;
  mcP2Counter[1] += p*p;
#define mccompcurname  Origin
#define mccompcurtype  Progress_bar
#define mccompcurindex 1
#define IntermediateCnts mccOrigin_IntermediateCnts
#define StartTime mccOrigin_StartTime
#define EndTime mccOrigin_EndTime
#define CurrentTime mccOrigin_CurrentTime
{   /* Declarations of Origin=Progress_bar() SETTING parameters. */
char* profile = mccOrigin_profile;
MCNUM percent = mccOrigin_percent;
MCNUM flag_save = mccOrigin_flag_save;
MCNUM minutes = mccOrigin_minutes;
#line 70 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../misc/Progress_bar.comp"
{
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
    if (flag_save) mcsave(NULL);
  }
}
#line 12865 "./Test_Monochromators.c"
}   /* End of Origin=Progress_bar() SETTING parameter declarations. */
#undef CurrentTime
#undef EndTime
#undef StartTime
#undef IntermediateCnts
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompOrigin:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(1,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlvx,
      mcnlvy,
      mcnlvz,
      mcnlt,
      mcnlsx,
      mcnlsy,
      mcnlsz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef sz
#undef sy
#undef sx
#undef t
#undef vz
#undef vy
#undef vx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlvx,
mcnlvy,
mcnlvz,
mcnlt,
mcnlsx,
mcnlsy,
mcnlsz,
mcnlp)

  /* TRACE Component Source [2] */
  mccoordschange(mcposrSource, mcrotrSource,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component Source (without coords transformations) */
  mcJumpTrace_Source:
  SIG_MESSAGE("Source (Trace)");
  mcDEBUG_COMP("Source")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompSource
  STORE_NEUTRON(2,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[2]++;
  mcPCounter[2] += p;
  mcP2Counter[2] += p*p;
#define mccompcurname  Source
#define mccompcurtype  Source_simple
#define mccompcurindex 2
#define pmul mccSource_pmul
#define square mccSource_square
#define srcArea mccSource_srcArea
{   /* Declarations of Source=Source_simple() SETTING parameters. */
MCNUM radius = mccSource_radius;
MCNUM yheight = mccSource_yheight;
MCNUM xwidth = mccSource_xwidth;
MCNUM dist = mccSource_dist;
MCNUM focus_xw = mccSource_focus_xw;
MCNUM focus_yh = mccSource_focus_yh;
MCNUM E0 = mccSource_E0;
MCNUM dE = mccSource_dE;
MCNUM lambda0 = mccSource_lambda0;
MCNUM dlambda = mccSource_dlambda;
MCNUM flux = mccSource_flux;
MCNUM gauss = mccSource_gauss;
int target_index = mccSource_target_index;
#line 125 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../sources/Source_simple.comp"
{
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
			  tx, ty, tz, focus_xw, focus_yh, ROT_A_CURRENT_COMP, x, y, z, 2);

 dx = xf-x;
 dy = yf-y;
 rf = sqrt(dx*dx+dy*dy+rf*rf);

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
}
#line 13036 "./Test_Monochromators.c"
}   /* End of Source=Source_simple() SETTING parameter declarations. */
#undef srcArea
#undef square
#undef pmul
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompSource:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(2,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlvx,
      mcnlvy,
      mcnlvz,
      mcnlt,
      mcnlsx,
      mcnlsy,
      mcnlsz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef sz
#undef sy
#undef sx
#undef t
#undef vz
#undef vy
#undef vx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlvx,
mcnlvy,
mcnlvz,
mcnlt,
mcnlsx,
mcnlsy,
mcnlsz,
mcnlp)

  /* TRACE Component lamStart [3] */
  mccoordschange(mcposrlamStart, mcrotrlamStart,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component lamStart (without coords transformations) */
  mcJumpTrace_lamStart:
  SIG_MESSAGE("lamStart (Trace)");
  mcDEBUG_COMP("lamStart")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbComplamStart
  STORE_NEUTRON(3,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[3]++;
  mcPCounter[3] += p;
  mcP2Counter[3] += p*p;
#define mccompcurname  lamStart
#define mccompcurtype  L_monitor
#define mccompcurindex 3
#define nL mcclamStart_nL
#define L_N mcclamStart_L_N
#define L_p mcclamStart_L_p
#define L_p2 mcclamStart_L_p2
{   /* Declarations of lamStart=L_monitor() SETTING parameters. */
char* filename = mcclamStart_filename;
MCNUM xmin = mcclamStart_xmin;
MCNUM xmax = mcclamStart_xmax;
MCNUM ymin = mcclamStart_ymin;
MCNUM ymax = mcclamStart_ymax;
MCNUM xwidth = mcclamStart_xwidth;
MCNUM yheight = mcclamStart_yheight;
MCNUM Lmin = mcclamStart_Lmin;
MCNUM Lmax = mcclamStart_Lmax;
MCNUM restore_neutron = mcclamStart_restore_neutron;
int nowritefile = mcclamStart_nowritefile;
#line 84 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/L_monitor.comp"
{
    int i;
    double L;

    PROP_Z0;
    if (x>xmin && x<xmax && y>ymin && y<ymax)
    {
      L = (2*PI/V2K)/sqrt(vx*vx + vy*vy + vz*vz);
      i = floor((L-Lmin)*nL/(Lmax-Lmin));
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
}
#line 13182 "./Test_Monochromators.c"
}   /* End of lamStart=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbComplamStart:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(3,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlvx,
      mcnlvy,
      mcnlvz,
      mcnlt,
      mcnlsx,
      mcnlsy,
      mcnlsz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef sz
#undef sy
#undef sx
#undef t
#undef vz
#undef vy
#undef vx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlvx,
mcnlvy,
mcnlvz,
mcnlt,
mcnlsx,
mcnlsy,
mcnlsz,
mcnlp)

  /* TRACE Component Mono_Arm [4] */
  mccoordschange(mcposrMono_Arm, mcrotrMono_Arm,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component Mono_Arm (without coords transformations) */
  mcJumpTrace_Mono_Arm:
  SIG_MESSAGE("Mono_Arm (Trace)");
  mcDEBUG_COMP("Mono_Arm")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompMono_Arm
  STORE_NEUTRON(4,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[4]++;
  mcPCounter[4] += p;
  mcP2Counter[4] += p*p;
#define mccompcurname  Mono_Arm
#define mccompcurtype  Arm
#define mccompcurindex 4
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompMono_Arm:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(4,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlvx,
      mcnlvy,
      mcnlvz,
      mcnlt,
      mcnlsx,
      mcnlsy,
      mcnlsz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef sz
#undef sy
#undef sx
#undef t
#undef vz
#undef vy
#undef vx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlvx,
mcnlvy,
mcnlvz,
mcnlt,
mcnlsx,
mcnlsy,
mcnlsz,
mcnlp)

  /* TRACE Component Mono_Out [5] */
  mccoordschange(mcposrMono_Out, mcrotrMono_Out,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component Mono_Out (without coords transformations) */
  mcJumpTrace_Mono_Out:
  SIG_MESSAGE("Mono_Out (Trace)");
  mcDEBUG_COMP("Mono_Out")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompMono_Out
  STORE_NEUTRON(5,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[5]++;
  mcPCounter[5] += p;
  mcP2Counter[5] += p*p;
#define mccompcurname  Mono_Out
#define mccompcurtype  Arm
#define mccompcurindex 5
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompMono_Out:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(5,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlvx,
      mcnlvy,
      mcnlvz,
      mcnlt,
      mcnlsx,
      mcnlsy,
      mcnlsz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef sz
#undef sy
#undef sx
#undef t
#undef vz
#undef vy
#undef vx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlvx,
mcnlvy,
mcnlvz,
mcnlt,
mcnlsx,
mcnlsy,
mcnlsz,
mcnlp)

  /* TRACE Component Mono1 [6] */
  mccoordschange(mcposrMono1, mcrotrMono1,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component Mono1 (without coords transformations) */
  mcJumpTrace_Mono1:
  SIG_MESSAGE("Mono1 (Trace)");
  mcDEBUG_COMP("Mono1")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompMono1
  STORE_NEUTRON(6,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[6]++;
  mcPCounter[6] += p;
  mcP2Counter[6] += p*p;
#define mccompcurname  Mono1
#define mccompcurtype  Monochromator_flat
#define mccompcurindex 6
#define mos_rms_y mccMono1_mos_rms_y
#define mos_rms_z mccMono1_mos_rms_z
#define mos_rms_max mccMono1_mos_rms_max
#define mono_Q mccMono1_mono_Q
{   /* Declarations of Mono1=Monochromator_flat() SETTING parameters. */
MCNUM zmin = mccMono1_zmin;
MCNUM zmax = mccMono1_zmax;
MCNUM ymin = mccMono1_ymin;
MCNUM ymax = mccMono1_ymax;
MCNUM zwidth = mccMono1_zwidth;
MCNUM yheight = mccMono1_yheight;
MCNUM mosaich = mccMono1_mosaich;
MCNUM mosaicv = mccMono1_mosaicv;
MCNUM r0 = mccMono1_r0;
MCNUM Q = mccMono1_Q;
MCNUM DM = mccMono1_DM;
/* 'Mono1=Monochromator_flat()' component instance has conditional execution */
if (( mcipMono == 1 ))

#line 118 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Monochromator_flat.comp"
{
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
}
#line 13650 "./Test_Monochromators.c"
}   /* End of Mono1=Monochromator_flat() SETTING parameter declarations. */
#undef mono_Q
#undef mos_rms_max
#undef mos_rms_z
#undef mos_rms_y
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompMono1:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(6,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlvx,
      mcnlvy,
      mcnlvz,
      mcnlt,
      mcnlsx,
      mcnlsy,
      mcnlsz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef sz
#undef sy
#undef sx
#undef t
#undef vz
#undef vy
#undef vx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlvx,
mcnlvy,
mcnlvz,
mcnlt,
mcnlsx,
mcnlsy,
mcnlsz,
mcnlp)

  /* TRACE Component Mono2 [7] */
  mccoordschange(mcposrMono2, mcrotrMono2,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component Mono2 (without coords transformations) */
  mcJumpTrace_Mono2:
  SIG_MESSAGE("Mono2 (Trace)");
  mcDEBUG_COMP("Mono2")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompMono2
  STORE_NEUTRON(7,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[7]++;
  mcPCounter[7] += p;
  mcP2Counter[7] += p*p;
#define mccompcurname  Mono2
#define mccompcurtype  Monochromator_pol
#define mccompcurindex 7
#define mos_rms mccMono2_mos_rms
#define d_rms mccMono2_d_rms
#define mono_Q mccMono2_mono_Q
#define FN mccMono2_FN
#define FM mccMono2_FM
{   /* Declarations of Mono2=Monochromator_pol() SETTING parameters. */
MCNUM zwidth = mccMono2_zwidth;
MCNUM yheight = mccMono2_yheight;
MCNUM mosaic = mccMono2_mosaic;
MCNUM dspread = mccMono2_dspread;
MCNUM Q = mccMono2_Q;
MCNUM DM = mccMono2_DM;
MCNUM pThreshold = mccMono2_pThreshold;
MCNUM Rup = mccMono2_Rup;
MCNUM Rdown = mccMono2_Rdown;
int debug = mccMono2_debug;
/* 'Mono2=Monochromator_pol()' component instance has conditional execution */
if (( mcipMono == 2 ))

#line 130 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Monochromator_pol.comp"
{
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

}
#line 13845 "./Test_Monochromators.c"
}   /* End of Mono2=Monochromator_pol() SETTING parameter declarations. */
#undef FM
#undef FN
#undef mono_Q
#undef d_rms
#undef mos_rms
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompMono2:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(7,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlvx,
      mcnlvy,
      mcnlvz,
      mcnlt,
      mcnlsx,
      mcnlsy,
      mcnlsz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef sz
#undef sy
#undef sx
#undef t
#undef vz
#undef vy
#undef vx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlvx,
mcnlvy,
mcnlvz,
mcnlt,
mcnlsx,
mcnlsy,
mcnlsz,
mcnlp)

  /* TRACE Component Mono3 [8] */
  mccoordschange(mcposrMono3, mcrotrMono3,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component Mono3 (without coords transformations) */
  mcJumpTrace_Mono3:
  SIG_MESSAGE("Mono3 (Trace)");
  mcDEBUG_COMP("Mono3")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompMono3
  STORE_NEUTRON(8,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[8]++;
  mcPCounter[8] += p;
  mcP2Counter[8] += p*p;
#define mccompcurname  Mono3
#define mccompcurtype  Monochromator_pol
#define mccompcurindex 8
#define mos_rms mccMono3_mos_rms
#define d_rms mccMono3_d_rms
#define mono_Q mccMono3_mono_Q
#define FN mccMono3_FN
#define FM mccMono3_FM
{   /* Declarations of Mono3=Monochromator_pol() SETTING parameters. */
MCNUM zwidth = mccMono3_zwidth;
MCNUM yheight = mccMono3_yheight;
MCNUM mosaic = mccMono3_mosaic;
MCNUM dspread = mccMono3_dspread;
MCNUM Q = mccMono3_Q;
MCNUM DM = mccMono3_DM;
MCNUM pThreshold = mccMono3_pThreshold;
MCNUM Rup = mccMono3_Rup;
MCNUM Rdown = mccMono3_Rdown;
int debug = mccMono3_debug;
/* 'Mono3=Monochromator_pol()' component instance has conditional execution */
if (( mcipMono == 3 ))

#line 130 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Monochromator_pol.comp"
{
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

}
#line 14041 "./Test_Monochromators.c"
}   /* End of Mono3=Monochromator_pol() SETTING parameter declarations. */
#undef FM
#undef FN
#undef mono_Q
#undef d_rms
#undef mos_rms
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompMono3:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(8,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlvx,
      mcnlvy,
      mcnlvz,
      mcnlt,
      mcnlsx,
      mcnlsy,
      mcnlsz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef sz
#undef sy
#undef sx
#undef t
#undef vz
#undef vy
#undef vx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlvx,
mcnlvy,
mcnlvz,
mcnlt,
mcnlsx,
mcnlsy,
mcnlsz,
mcnlp)

  /* TRACE Component Mono4 [9] */
  mccoordschange(mcposrMono4, mcrotrMono4,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component Mono4 (without coords transformations) */
  mcJumpTrace_Mono4:
  SIG_MESSAGE("Mono4 (Trace)");
  mcDEBUG_COMP("Mono4")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompMono4
  STORE_NEUTRON(9,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[9]++;
  mcPCounter[9] += p;
  mcP2Counter[9] += p*p;
#define mccompcurname  Mono4
#define mccompcurtype  Monochromator_curved
#define mccompcurindex 9
#define mos_rms_y mccMono4_mos_rms_y
#define mos_rms_z mccMono4_mos_rms_z
#define mos_rms_max mccMono4_mos_rms_max
#define mono_Q mccMono4_mono_Q
#define SlabWidth mccMono4_SlabWidth
#define SlabHeight mccMono4_SlabHeight
#define rTable mccMono4_rTable
#define tTable mccMono4_tTable
#define row mccMono4_row
#define col mccMono4_col
#define tiltH mccMono4_tiltH
#define tiltV mccMono4_tiltV
{   /* Declarations of Mono4=Monochromator_curved() SETTING parameters. */
char* reflect = mccMono4_reflect;
char* transmit = mccMono4_transmit;
MCNUM zwidth = mccMono4_zwidth;
MCNUM yheight = mccMono4_yheight;
MCNUM gap = mccMono4_gap;
MCNUM NH = mccMono4_NH;
MCNUM NV = mccMono4_NV;
MCNUM mosaich = mccMono4_mosaich;
MCNUM mosaicv = mccMono4_mosaicv;
MCNUM r0 = mccMono4_r0;
MCNUM t0 = mccMono4_t0;
MCNUM Q = mccMono4_Q;
MCNUM RV = mccMono4_RV;
MCNUM RH = mccMono4_RH;
MCNUM DM = mccMono4_DM;
MCNUM mosaic = mccMono4_mosaic;
MCNUM width = mccMono4_width;
MCNUM height = mccMono4_height;
MCNUM verbose = mccMono4_verbose;
MCNUM order = mccMono4_order;
/* 'Mono4=Monochromator_curved()' component instance has conditional execution */
if (( mcipMono == 4 ))

#line 230 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Monochromator_curved.comp"
{
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
      /* Visualise scattering point in proper, component frame 
	 - but only if the neutron is reflected, that is none of:
	 * transmitted
	 * falling outside the slab material */	
      if(!do_transmit) SCATTER;

      /* mccoordschange_polarisation(tt, &sx, &sy, &sz); */
    } /* End intersect the crystal (if z) */
    else {
      /* restore neutron state when no interaction */
      RESTORE_NEUTRON(INDEX_CURRENT_COMP, x, y, z, vx, vy, vz, t, sx, sy, sz, p);
    }
  } /* End neutron moving towards crystal (if vx)*/
}
#line 14415 "./Test_Monochromators.c"
}   /* End of Mono4=Monochromator_curved() SETTING parameter declarations. */
#undef tiltV
#undef tiltH
#undef col
#undef row
#undef tTable
#undef rTable
#undef SlabHeight
#undef SlabWidth
#undef mono_Q
#undef mos_rms_max
#undef mos_rms_z
#undef mos_rms_y
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompMono4:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(9,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlvx,
      mcnlvy,
      mcnlvz,
      mcnlt,
      mcnlsx,
      mcnlsy,
      mcnlsz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef sz
#undef sy
#undef sx
#undef t
#undef vz
#undef vy
#undef vx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlvx,
mcnlvy,
mcnlvz,
mcnlt,
mcnlsx,
mcnlsy,
mcnlsz,
mcnlp)

  /* TRACE Component Mono5 [10] */
  mccoordschange(mcposrMono5, mcrotrMono5,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component Mono5 (without coords transformations) */
  mcJumpTrace_Mono5:
  SIG_MESSAGE("Mono5 (Trace)");
  mcDEBUG_COMP("Mono5")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompMono5
  STORE_NEUTRON(10,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[10]++;
  mcPCounter[10] += p;
  mcP2Counter[10] += p*p;
#define mccompcurname  Mono5
#define mccompcurtype  Monochromator_2foc
#define mccompcurindex 10
#define mos_y mccMono5_mos_y
#define mos_z mccMono5_mos_z
#define mono_Q mccMono5_mono_Q
#define SlabWidth mccMono5_SlabWidth
#define SlabHeight mccMono5_SlabHeight
#define rTable mccMono5_rTable
{   /* Declarations of Mono5=Monochromator_2foc() SETTING parameters. */
char* reflect = mccMono5_reflect;
MCNUM zwidth = mccMono5_zwidth;
MCNUM yheight = mccMono5_yheight;
MCNUM gap = mccMono5_gap;
MCNUM NH = mccMono5_NH;
MCNUM NV = mccMono5_NV;
MCNUM mosaich = mccMono5_mosaich;
MCNUM mosaicv = mccMono5_mosaicv;
MCNUM r0 = mccMono5_r0;
MCNUM Q = mccMono5_Q;
MCNUM RV = mccMono5_RV;
MCNUM RH = mccMono5_RH;
MCNUM DM = mccMono5_DM;
MCNUM mosaic = mccMono5_mosaic;
MCNUM width = mccMono5_width;
MCNUM height = mccMono5_height;
MCNUM verbose = mccMono5_verbose;
/* 'Mono5=Monochromator_2foc()' component instance has conditional execution */
if (( mcipMono == 5 ))

#line 148 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../contrib/Monochromator_2foc.comp"
{
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
}
#line 14679 "./Test_Monochromators.c"
}   /* End of Mono5=Monochromator_2foc() SETTING parameter declarations. */
#undef rTable
#undef SlabHeight
#undef SlabWidth
#undef mono_Q
#undef mos_z
#undef mos_y
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompMono5:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(10,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlvx,
      mcnlvy,
      mcnlvz,
      mcnlt,
      mcnlsx,
      mcnlsy,
      mcnlsz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef sz
#undef sy
#undef sx
#undef t
#undef vz
#undef vy
#undef vx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlvx,
mcnlvy,
mcnlvz,
mcnlt,
mcnlsx,
mcnlsy,
mcnlsz,
mcnlp)

  /* TRACE Component Mono6 [11] */
  mccoordschange(mcposrMono6, mcrotrMono6,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component Mono6 (without coords transformations) */
  mcJumpTrace_Mono6:
  SIG_MESSAGE("Mono6 (Trace)");
  mcDEBUG_COMP("Mono6")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompMono6
  STORE_NEUTRON(11,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[11]++;
  mcPCounter[11] += p;
  mcP2Counter[11] += p*p;
#define mccompcurname  Mono6
#define mccompcurtype  Single_crystal
#define mccompcurindex 11
#define mosaic_AB mccMono6_mosaic_AB
#define hkl_info mccMono6_hkl_info
#define offdata mccMono6_offdata
{   /* Declarations of Mono6=Single_crystal() SETTING parameters. */
char* reflections = mccMono6_reflections;
char* geometry = mccMono6_geometry;
MCNUM xwidth = mccMono6_xwidth;
MCNUM yheight = mccMono6_yheight;
MCNUM zdepth = mccMono6_zdepth;
MCNUM radius = mccMono6_radius;
MCNUM delta_d_d = mccMono6_delta_d_d;
MCNUM mosaic = mccMono6_mosaic;
MCNUM mosaic_a = mccMono6_mosaic_a;
MCNUM mosaic_b = mccMono6_mosaic_b;
MCNUM mosaic_c = mccMono6_mosaic_c;
MCNUM recip_cell = mccMono6_recip_cell;
MCNUM barns = mccMono6_barns;
MCNUM ax = mccMono6_ax;
MCNUM ay = mccMono6_ay;
MCNUM az = mccMono6_az;
MCNUM bx = mccMono6_bx;
MCNUM by = mccMono6_by;
MCNUM bz = mccMono6_bz;
MCNUM cx = mccMono6_cx;
MCNUM cy = mccMono6_cy;
MCNUM cz = mccMono6_cz;
MCNUM p_transmit = mccMono6_p_transmit;
MCNUM sigma_abs = mccMono6_sigma_abs;
MCNUM sigma_inc = mccMono6_sigma_inc;
MCNUM aa = mccMono6_aa;
MCNUM bb = mccMono6_bb;
MCNUM cc = mccMono6_cc;
MCNUM order = mccMono6_order;
MCNUM RX = mccMono6_RX;
MCNUM RY = mccMono6_RY;
MCNUM powder = mccMono6_powder;
MCNUM PG = mccMono6_PG;
MCNUM deltak = mccMono6_deltak;
/* 'Mono6=Single_crystal()' component instance has conditional execution */
if (( mcipMono == 6 ))

#line 1065 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../samples/Single_crystal.comp"
{
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
  double adjust, r, sum;        /* Temporaries */

  double p_trans;               /* Transmission probability */
  double mc_trans, mc_interact; /* Transmission, interaction MC choices */
  int    intersect=0;
  double theta, phi;            /* rotation angles for curved lattice option */

  double curv_xangle;
  double curv_yangle;

  double _vx;
  double _vy;
  double _vz;


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
      double Alpha, Beta, Gamma;

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
        Beta  = randpm1()*PI/2;
        Gamma = randpm1()*PI;
        randrotate(&vx, &vy, &vz, Alpha, Beta, Gamma);
      }
      if (PG) { /* orientation of crystallite is random along <c> axis */
        Alpha = rand01()*2*PI*PG;
        PGrotate(&vx, &vy, &vz, Alpha, hkl_info.csx, hkl_info.csy, hkl_info.csz);
      }



      /* ------------------------------------------------------------------------- */
      /* lattice curvature option: rotate neutron velocity */
      /* WARNING: cannot be used together with the PG c-rotation! */
      curv_xangle = 0;
      curv_yangle = 0;

      _vx = vx;
      _vy = vy;
      _vz = vz;

      if(RY) { /* rotate v around x axis based on y pos, for vertical focus */
          curv_yangle = atan2(y, RY);
          vec_rotate_2d(&vy,&vz, curv_yangle);
          vec_rotate_2d(&sy,&sz, curv_yangle);

          /*changing y,z actually curves the crystal, not only the planes*/
          /*comment out if only curvature of the lattice planes is needed*/
          vec_rotate_2d(&y,&z, curv_yangle);
      }
      if(RX) { /* rotate v around y axis based on x pos, for horizontal focus */
          curv_xangle = atan2(x, RX);
          vec_rotate_2d(&vx,&vz, curv_xangle);
          vec_rotate_2d(&sx,&sz, curv_xangle);

          /*changing x,z actually curves the crystal, not only the planes*/
          /*comment out if only curvature of the lattice planes is needed*/
          vec_rotate_2d(&x,&z, curv_xangle);
      }

      kix = V2K*vx;
      kiy = V2K*vy;
      kiz = V2K*vz;
      vx = _vx;
      vy = _vy;
      vz = _vz;
      /* ------------------------------------------------------------------------- */



      /* (2). Intersection of Ewald sphere with reciprocal lattice points */

      /* in case we use 'SPLIT' then consecutive neutrons can be identical when entering here
         and we may skip the hkl_search call */
      if ( fabs(kix - hkl_info.kix) < deltak
        && fabs(kiy - hkl_info.kiy) < deltak
        && fabs(kiz - hkl_info.kiz) < deltak) {
        hkl_info.nb_reuses++;
      } else {
        /* Max possible tau for this ki with 5*sigma delta-d/d cutoff. */
        double tau_max   = 2*ki/(1 - 5*hkl_info.m_delta_d_d);
        double coh_xsect = 0, coh_refl = 0;

        /* call hkl_search */
        hkl_info.tau_count = hkl_search(L, T, hkl_info.count, hkl_info.V0,
              kix, kiy, kiz, tau_max,
              &coh_refl, &coh_xsect); /* CPU consuming */

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
          randderotate(&vx, &vy, &vz, Alpha, Beta, Gamma);
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
            randderotate(&vx, &vy, &vz, Alpha, Beta, Gamma);
          }
          if (PG) { /* orientation of crystallite is no longer random, rotation around <c> */
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
              "(r=%g, sum=%g, j=%i, tau_count=%i).\n", r, sum, j , hkl_info.tau_count);
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



      /* ------------------------------------------------------------------------- */
      /* lattice curvature option: rotate back neutron velocity */
      if(RX) {
          vec_rotate_2d(&vx,&vz, -curv_xangle);
          vec_rotate_2d(&sx,&sz, -curv_xangle);

          /*changing x,z actually curves the crystal, not only the planes*/
          /*comment out if only curvature of the lattice planes is needed*/
          vec_rotate_2d(&x,&z, -curv_xangle);
      }
      if(RY) {
          vec_rotate_2d(&vy,&vz, -curv_yangle);
          vec_rotate_2d(&sy,&sz, -curv_yangle);

          /*changing y,z actually curves the crystal, not only the planes*/
          /*comment out if only curvature of the lattice planes is needed*/
          vec_rotate_2d(&y,&z, -curv_yangle);
      }
      /* ------------------------------------------------------------------------- */



      SCATTER;
      if (powder) { /* orientation of crystallite is no longer random */
        randderotate(&vx, &vy, &vz, Alpha, Beta, Gamma);
      }
      if (PG) { /* orientation of crystallite is longer random */
        PGderotate(&vx, &vy, &vz, Alpha, hkl_info.csx, hkl_info.csy, hkl_info.csz);
      }
      /* exit if multiple scattering order has been reached */
      if (order && event_counter >= order) { intersect=0; break; }
      /* Repeat loop for next scattering event. */
    } while (intersect); /* end do (intersect) (multiple scattering loop) */
  } /* if intersect */
}
#line 15173 "./Test_Monochromators.c"
/* 'Mono6=Single_crystal()' component instance extend code */
    SIG_MESSAGE("Mono6 (Trace:Extend)");
if (( mcipMono == 6 )) {

#line 158 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Test_Monochromators/Test_Monochromators.instr"
  if(!SCATTERED) ABSORB;
#line 15179 "./Test_Monochromators.c"
}

}   /* End of Mono6=Single_crystal() SETTING parameter declarations. */
#undef offdata
#undef hkl_info
#undef mosaic_AB
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompMono6:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(11,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlvx,
      mcnlvy,
      mcnlvz,
      mcnlt,
      mcnlsx,
      mcnlsy,
      mcnlsz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef sz
#undef sy
#undef sx
#undef t
#undef vz
#undef vy
#undef vx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlvx,
mcnlvy,
mcnlvz,
mcnlt,
mcnlsx,
mcnlsy,
mcnlsz,
mcnlp)

  /* TRACE Component Sphere1 [12] */
  mccoordschange(mcposrSphere1, mcrotrSphere1,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component Sphere1 (without coords transformations) */
  mcJumpTrace_Sphere1:
  SIG_MESSAGE("Sphere1 (Trace)");
  mcDEBUG_COMP("Sphere1")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompSphere1
  STORE_NEUTRON(12,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[12]++;
  mcPCounter[12] += p;
  mcP2Counter[12] += p*p;
#define mccompcurname  Sphere1
#define mccompcurtype  PSD_monitor_4PI
#define mccompcurindex 12
#define nx mccSphere1_nx
#define ny mccSphere1_ny
#define PSD_N mccSphere1_PSD_N
#define PSD_p mccSphere1_PSD_p
#define PSD_p2 mccSphere1_PSD_p2
{   /* Declarations of Sphere1=PSD_monitor_4PI() SETTING parameters. */
char* filename = mccSphere1_filename;
MCNUM radius = mccSphere1_radius;
MCNUM restore_neutron = mccSphere1_restore_neutron;
int nowritefile = mccSphere1_nowritefile;
#line 74 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor_4PI.comp"
{
  double t0, t1, phi, theta;
  int i,j;

  if(sphere_intersect(&t0, &t1, x, y, z, vx, vy, vz, radius) && t1 > 0)
  {
    if(t0 < 0)
      t0 = t1;
    /* t0 is now time of intersection with the sphere. */
    PROP_DT(t0);
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
}
#line 15330 "./Test_Monochromators.c"
}   /* End of Sphere1=PSD_monitor_4PI() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef ny
#undef nx
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompSphere1:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(12,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlvx,
      mcnlvy,
      mcnlvz,
      mcnlt,
      mcnlsx,
      mcnlsy,
      mcnlsz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef sz
#undef sy
#undef sx
#undef t
#undef vz
#undef vy
#undef vx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlvx,
mcnlvy,
mcnlvz,
mcnlt,
mcnlsx,
mcnlsy,
mcnlsz,
mcnlp)

  /* TRACE Component lam1 [13] */
  mccoordschange(mcposrlam1, mcrotrlam1,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component lam1 (without coords transformations) */
  mcJumpTrace_lam1:
  SIG_MESSAGE("lam1 (Trace)");
  mcDEBUG_COMP("lam1")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbComplam1
  STORE_NEUTRON(13,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[13]++;
  mcPCounter[13] += p;
  mcP2Counter[13] += p*p;
#define mccompcurname  lam1
#define mccompcurtype  L_monitor
#define mccompcurindex 13
#define nL mcclam1_nL
#define L_N mcclam1_L_N
#define L_p mcclam1_L_p
#define L_p2 mcclam1_L_p2
{   /* Declarations of lam1=L_monitor() SETTING parameters. */
char* filename = mcclam1_filename;
MCNUM xmin = mcclam1_xmin;
MCNUM xmax = mcclam1_xmax;
MCNUM ymin = mcclam1_ymin;
MCNUM ymax = mcclam1_ymax;
MCNUM xwidth = mcclam1_xwidth;
MCNUM yheight = mcclam1_yheight;
MCNUM Lmin = mcclam1_Lmin;
MCNUM Lmax = mcclam1_Lmax;
MCNUM restore_neutron = mcclam1_restore_neutron;
int nowritefile = mcclam1_nowritefile;
#line 84 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/L_monitor.comp"
{
    int i;
    double L;

    PROP_Z0;
    if (x>xmin && x<xmax && y>ymin && y<ymax)
    {
      L = (2*PI/V2K)/sqrt(vx*vx + vy*vy + vz*vz);
      i = floor((L-Lmin)*nL/(Lmax-Lmin));
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
}
#line 15478 "./Test_Monochromators.c"
}   /* End of lam1=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbComplam1:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(13,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlvx,
      mcnlvy,
      mcnlvz,
      mcnlt,
      mcnlsx,
      mcnlsy,
      mcnlsz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef sz
#undef sy
#undef sx
#undef t
#undef vz
#undef vy
#undef vx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlvx,
mcnlvy,
mcnlvz,
mcnlt,
mcnlsx,
mcnlsy,
mcnlsz,
mcnlp)

  /* TRACE Component psd1 [14] */
  mccoordschange(mcposrpsd1, mcrotrpsd1,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component psd1 (without coords transformations) */
  mcJumpTrace_psd1:
  SIG_MESSAGE("psd1 (Trace)");
  mcDEBUG_COMP("psd1")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbComppsd1
  STORE_NEUTRON(14,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[14]++;
  mcPCounter[14] += p;
  mcP2Counter[14] += p*p;
#define mccompcurname  psd1
#define mccompcurtype  PSD_monitor
#define mccompcurindex 14
#define PSD_N mccpsd1_PSD_N
#define PSD_p mccpsd1_PSD_p
#define PSD_p2 mccpsd1_PSD_p2
{   /* Declarations of psd1=PSD_monitor() SETTING parameters. */
int nx = mccpsd1_nx;
int ny = mccpsd1_ny;
char* filename = mccpsd1_filename;
MCNUM xmin = mccpsd1_xmin;
MCNUM xmax = mccpsd1_xmax;
MCNUM ymin = mccpsd1_ymin;
MCNUM ymax = mccpsd1_ymax;
MCNUM xwidth = mccpsd1_xwidth;
MCNUM yheight = mccpsd1_yheight;
MCNUM restore_neutron = mccpsd1_restore_neutron;
#line 94 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
{
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
}
#line 15616 "./Test_Monochromators.c"
}   /* End of psd1=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbComppsd1:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(14,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlvx,
      mcnlvy,
      mcnlvz,
      mcnlt,
      mcnlsx,
      mcnlsy,
      mcnlsz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef sz
#undef sy
#undef sx
#undef t
#undef vz
#undef vy
#undef vx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlvx,
mcnlvy,
mcnlvz,
mcnlt,
mcnlsx,
mcnlsy,
mcnlsz,
mcnlp)

  mcabsorbAll:
  mcDEBUG_LEAVE()
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlvx,
mcnlvy,
mcnlvz,
mcnlt,
mcnlsx,
mcnlsy,
mcnlsz,
mcnlp)
  /* Copy neutron state to global variables. */
  mcnx = mcnlx;
  mcny = mcnly;
  mcnz = mcnlz;
  mcnvx = mcnlvx;
  mcnvy = mcnlvy;
  mcnvz = mcnlvz;
  mcnt = mcnlt;
  mcnsx = mcnlsx;
  mcnsy = mcnlsy;
  mcnsz = mcnlsz;
  mcnp = mcnlp;

} /* end trace */

void mcsave(FILE *handle) {
  if (!handle) mcsiminfo_init(NULL);
  /* User component SAVE code. */

  /* User SAVE code for component 'Origin'. */
  SIG_MESSAGE("Origin (Save)");
#define mccompcurname  Origin
#define mccompcurtype  Progress_bar
#define mccompcurindex 1
#define IntermediateCnts mccOrigin_IntermediateCnts
#define StartTime mccOrigin_StartTime
#define EndTime mccOrigin_EndTime
#define CurrentTime mccOrigin_CurrentTime
{   /* Declarations of Origin=Progress_bar() SETTING parameters. */
char* profile = mccOrigin_profile;
MCNUM percent = mccOrigin_percent;
MCNUM flag_save = mccOrigin_flag_save;
MCNUM minutes = mccOrigin_minutes;
#line 115 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../misc/Progress_bar.comp"
{
  MPI_MASTER(fprintf(stdout, "\nSave [%s]\n", mcinstrument_name););
  if (profile && strlen(profile) && strcmp(profile,"NULL") && strcmp(profile,"0")) {
    char filename[256];
    if (!strlen(profile) || !strcmp(profile,"NULL") || !strcmp(profile,"0")) strcpy(filename, mcinstrument_name);
    else strcpy(filename, profile);
    DETECTOR_OUT_1D(
        "Intensity profiler",
        "Component index [1]",
        "Intensity",
        "prof", 1, mcNUMCOMP, mcNUMCOMP-1,
        &mcNCounter[1],&mcPCounter[1],&mcP2Counter[1],
        filename);

  }
}
#line 15728 "./Test_Monochromators.c"
}   /* End of Origin=Progress_bar() SETTING parameter declarations. */
#undef CurrentTime
#undef EndTime
#undef StartTime
#undef IntermediateCnts
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'lamStart'. */
  SIG_MESSAGE("lamStart (Save)");
#define mccompcurname  lamStart
#define mccompcurtype  L_monitor
#define mccompcurindex 3
#define nL mcclamStart_nL
#define L_N mcclamStart_L_N
#define L_p mcclamStart_L_p
#define L_p2 mcclamStart_L_p2
{   /* Declarations of lamStart=L_monitor() SETTING parameters. */
char* filename = mcclamStart_filename;
MCNUM xmin = mcclamStart_xmin;
MCNUM xmax = mcclamStart_xmax;
MCNUM ymin = mcclamStart_ymin;
MCNUM ymax = mcclamStart_ymax;
MCNUM xwidth = mcclamStart_xwidth;
MCNUM yheight = mcclamStart_yheight;
MCNUM Lmin = mcclamStart_Lmin;
MCNUM Lmax = mcclamStart_Lmax;
MCNUM restore_neutron = mcclamStart_restore_neutron;
int nowritefile = mcclamStart_nowritefile;
#line 107 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/L_monitor.comp"
{
    if (!nowritefile) {
    DETECTOR_OUT_1D(
        "Wavelength monitor",
        "Wavelength [AA]",
        "Intensity",
        "L", Lmin, Lmax, nL,
        &L_N[0],&L_p[0],&L_p2[0],
        filename);
    }
}
#line 15771 "./Test_Monochromators.c"
}   /* End of lamStart=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'Sphere1'. */
  SIG_MESSAGE("Sphere1 (Save)");
#define mccompcurname  Sphere1
#define mccompcurtype  PSD_monitor_4PI
#define mccompcurindex 12
#define nx mccSphere1_nx
#define ny mccSphere1_ny
#define PSD_N mccSphere1_PSD_N
#define PSD_p mccSphere1_PSD_p
#define PSD_p2 mccSphere1_PSD_p2
{   /* Declarations of Sphere1=PSD_monitor_4PI() SETTING parameters. */
char* filename = mccSphere1_filename;
MCNUM radius = mccSphere1_radius;
MCNUM restore_neutron = mccSphere1_restore_neutron;
int nowritefile = mccSphere1_nowritefile;
#line 107 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor_4PI.comp"
{
  if (!nowritefile) {
    DETECTOR_OUT_2D(
    "4PI PSD monitor",
    "Longitude [deg]",
    "Lattitude [deg]",
    -180, 180, -90, 90,
    nx, ny,
    &PSD_N[0][0],&PSD_p[0][0],&PSD_p2[0][0],
    filename);
  }
}
#line 15809 "./Test_Monochromators.c"
}   /* End of Sphere1=PSD_monitor_4PI() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef ny
#undef nx
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'lam1'. */
  SIG_MESSAGE("lam1 (Save)");
#define mccompcurname  lam1
#define mccompcurtype  L_monitor
#define mccompcurindex 13
#define nL mcclam1_nL
#define L_N mcclam1_L_N
#define L_p mcclam1_L_p
#define L_p2 mcclam1_L_p2
{   /* Declarations of lam1=L_monitor() SETTING parameters. */
char* filename = mcclam1_filename;
MCNUM xmin = mcclam1_xmin;
MCNUM xmax = mcclam1_xmax;
MCNUM ymin = mcclam1_ymin;
MCNUM ymax = mcclam1_ymax;
MCNUM xwidth = mcclam1_xwidth;
MCNUM yheight = mcclam1_yheight;
MCNUM Lmin = mcclam1_Lmin;
MCNUM Lmax = mcclam1_Lmax;
MCNUM restore_neutron = mcclam1_restore_neutron;
int nowritefile = mcclam1_nowritefile;
#line 107 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/L_monitor.comp"
{
    if (!nowritefile) {
    DETECTOR_OUT_1D(
        "Wavelength monitor",
        "Wavelength [AA]",
        "Intensity",
        "L", Lmin, Lmax, nL,
        &L_N[0],&L_p[0],&L_p2[0],
        filename);
    }
}
#line 15853 "./Test_Monochromators.c"
}   /* End of lam1=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'psd1'. */
  SIG_MESSAGE("psd1 (Save)");
#define mccompcurname  psd1
#define mccompcurtype  PSD_monitor
#define mccompcurindex 14
#define PSD_N mccpsd1_PSD_N
#define PSD_p mccpsd1_PSD_p
#define PSD_p2 mccpsd1_PSD_p2
{   /* Declarations of psd1=PSD_monitor() SETTING parameters. */
int nx = mccpsd1_nx;
int ny = mccpsd1_ny;
char* filename = mccpsd1_filename;
MCNUM xmin = mccpsd1_xmin;
MCNUM xmax = mccpsd1_xmax;
MCNUM ymin = mccpsd1_ymin;
MCNUM ymax = mccpsd1_ymax;
MCNUM xwidth = mccpsd1_xwidth;
MCNUM yheight = mccpsd1_yheight;
MCNUM restore_neutron = mccpsd1_restore_neutron;
#line 110 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
{
  DETECTOR_OUT_2D(
    "PSD monitor",
    "X position [cm]",
    "Y position [cm]",
    xmin*100.0, xmax*100.0, ymin*100.0, ymax*100.0,
    nx, ny,
    &PSD_N[0][0],&PSD_p[0][0],&PSD_p2[0][0],
    filename);
}
#line 15893 "./Test_Monochromators.c"
}   /* End of psd1=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  if (!handle) mcsiminfo_close(); 
} /* end save */
void mcfinally(void) {
  /* User component FINALLY code. */
  mcsiminfo_init(NULL);
  mcsave(mcsiminfo_file); /* save data when simulation ends */

  /* User FINALLY code for component 'Origin'. */
  SIG_MESSAGE("Origin (Finally)");
#define mccompcurname  Origin
#define mccompcurtype  Progress_bar
#define mccompcurindex 1
#define IntermediateCnts mccOrigin_IntermediateCnts
#define StartTime mccOrigin_StartTime
#define EndTime mccOrigin_EndTime
#define CurrentTime mccOrigin_CurrentTime
{   /* Declarations of Origin=Progress_bar() SETTING parameters. */
char* profile = mccOrigin_profile;
MCNUM percent = mccOrigin_percent;
MCNUM flag_save = mccOrigin_flag_save;
MCNUM minutes = mccOrigin_minutes;
#line 133 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../misc/Progress_bar.comp"
{
  time_t NowTime;
  time(&NowTime);
  fprintf(stdout, "\nFinally [%s: %s]. Time: ", mcinstrument_name, mcdirname ? mcdirname : ".");
  if (difftime(NowTime,StartTime) < 60.0)
    fprintf(stdout, "%g [s] ", difftime(NowTime,StartTime));
  else if (difftime(NowTime,StartTime) > 3600.0)
    fprintf(stdout, "%g [h] ", difftime(NowTime,StartTime)/3660.0);
  else
    fprintf(stdout, "%g [min] ", difftime(NowTime,StartTime)/60.0);
  fprintf(stdout, "\n");
}
#line 15936 "./Test_Monochromators.c"
}   /* End of Origin=Progress_bar() SETTING parameter declarations. */
#undef CurrentTime
#undef EndTime
#undef StartTime
#undef IntermediateCnts
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[1]) fprintf(stderr, "Warning: No neutron could reach Component[1] Origin\n");
    if (mcAbsorbProp[1]) fprintf(stderr, "Warning: %g events were removed in Component[1] Origin=Progress_bar()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[1]);
    if (!mcNCounter[2]) fprintf(stderr, "Warning: No neutron could reach Component[2] Source\n");
    if (mcAbsorbProp[2]) fprintf(stderr, "Warning: %g events were removed in Component[2] Source=Source_simple()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[2]);
    if (!mcNCounter[3]) fprintf(stderr, "Warning: No neutron could reach Component[3] lamStart\n");
    if (mcAbsorbProp[3]) fprintf(stderr, "Warning: %g events were removed in Component[3] lamStart=L_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[3]);
    if (!mcNCounter[4]) fprintf(stderr, "Warning: No neutron could reach Component[4] Mono_Arm\n");
    if (mcAbsorbProp[4]) fprintf(stderr, "Warning: %g events were removed in Component[4] Mono_Arm=Arm()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[4]);
    if (!mcNCounter[5]) fprintf(stderr, "Warning: No neutron could reach Component[5] Mono_Out\n");
    if (mcAbsorbProp[5]) fprintf(stderr, "Warning: %g events were removed in Component[5] Mono_Out=Arm()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[5]);
    if (!mcNCounter[6]) fprintf(stderr, "Warning: No neutron could reach Component[6] Mono1\n");
    if (mcAbsorbProp[6]) fprintf(stderr, "Warning: %g events were removed in Component[6] Mono1=Monochromator_flat()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[6]);
    if (!mcNCounter[7]) fprintf(stderr, "Warning: No neutron could reach Component[7] Mono2\n");
    if (mcAbsorbProp[7]) fprintf(stderr, "Warning: %g events were removed in Component[7] Mono2=Monochromator_pol()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[7]);
    if (!mcNCounter[8]) fprintf(stderr, "Warning: No neutron could reach Component[8] Mono3\n");
    if (mcAbsorbProp[8]) fprintf(stderr, "Warning: %g events were removed in Component[8] Mono3=Monochromator_pol()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[8]);
  /* User FINALLY code for component 'Mono4'. */
  SIG_MESSAGE("Mono4 (Finally)");
#define mccompcurname  Mono4
#define mccompcurtype  Monochromator_curved
#define mccompcurindex 9
#define mos_rms_y mccMono4_mos_rms_y
#define mos_rms_z mccMono4_mos_rms_z
#define mos_rms_max mccMono4_mos_rms_max
#define mono_Q mccMono4_mono_Q
#define SlabWidth mccMono4_SlabWidth
#define SlabHeight mccMono4_SlabHeight
#define rTable mccMono4_rTable
#define tTable mccMono4_tTable
#define row mccMono4_row
#define col mccMono4_col
#define tiltH mccMono4_tiltH
#define tiltV mccMono4_tiltV
{   /* Declarations of Mono4=Monochromator_curved() SETTING parameters. */
char* reflect = mccMono4_reflect;
char* transmit = mccMono4_transmit;
MCNUM zwidth = mccMono4_zwidth;
MCNUM yheight = mccMono4_yheight;
MCNUM gap = mccMono4_gap;
MCNUM NH = mccMono4_NH;
MCNUM NV = mccMono4_NV;
MCNUM mosaich = mccMono4_mosaich;
MCNUM mosaicv = mccMono4_mosaicv;
MCNUM r0 = mccMono4_r0;
MCNUM t0 = mccMono4_t0;
MCNUM Q = mccMono4_Q;
MCNUM RV = mccMono4_RV;
MCNUM RH = mccMono4_RH;
MCNUM DM = mccMono4_DM;
MCNUM mosaic = mccMono4_mosaic;
MCNUM width = mccMono4_width;
MCNUM height = mccMono4_height;
MCNUM verbose = mccMono4_verbose;
MCNUM order = mccMono4_order;
#line 460 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Monochromator_curved.comp"
{
  Table_Free(&rTable);
  Table_Free(&tTable);
  if (tiltH) free(tiltH);
  if (tiltV) free(tiltV);
}
#line 16007 "./Test_Monochromators.c"
}   /* End of Mono4=Monochromator_curved() SETTING parameter declarations. */
#undef tiltV
#undef tiltH
#undef col
#undef row
#undef tTable
#undef rTable
#undef SlabHeight
#undef SlabWidth
#undef mono_Q
#undef mos_rms_max
#undef mos_rms_z
#undef mos_rms_y
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[9]) fprintf(stderr, "Warning: No neutron could reach Component[9] Mono4\n");
    if (mcAbsorbProp[9]) fprintf(stderr, "Warning: %g events were removed in Component[9] Mono4=Monochromator_curved()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[9]);
    if (!mcNCounter[10]) fprintf(stderr, "Warning: No neutron could reach Component[10] Mono5\n");
    if (mcAbsorbProp[10]) fprintf(stderr, "Warning: %g events were removed in Component[10] Mono5=Monochromator_2foc()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[10]);
  /* User FINALLY code for component 'Mono6'. */
  SIG_MESSAGE("Mono6 (Finally)");
#define mccompcurname  Mono6
#define mccompcurtype  Single_crystal
#define mccompcurindex 11
#define mosaic_AB mccMono6_mosaic_AB
#define hkl_info mccMono6_hkl_info
#define offdata mccMono6_offdata
{   /* Declarations of Mono6=Single_crystal() SETTING parameters. */
char* reflections = mccMono6_reflections;
char* geometry = mccMono6_geometry;
MCNUM xwidth = mccMono6_xwidth;
MCNUM yheight = mccMono6_yheight;
MCNUM zdepth = mccMono6_zdepth;
MCNUM radius = mccMono6_radius;
MCNUM delta_d_d = mccMono6_delta_d_d;
MCNUM mosaic = mccMono6_mosaic;
MCNUM mosaic_a = mccMono6_mosaic_a;
MCNUM mosaic_b = mccMono6_mosaic_b;
MCNUM mosaic_c = mccMono6_mosaic_c;
MCNUM recip_cell = mccMono6_recip_cell;
MCNUM barns = mccMono6_barns;
MCNUM ax = mccMono6_ax;
MCNUM ay = mccMono6_ay;
MCNUM az = mccMono6_az;
MCNUM bx = mccMono6_bx;
MCNUM by = mccMono6_by;
MCNUM bz = mccMono6_bz;
MCNUM cx = mccMono6_cx;
MCNUM cy = mccMono6_cy;
MCNUM cz = mccMono6_cz;
MCNUM p_transmit = mccMono6_p_transmit;
MCNUM sigma_abs = mccMono6_sigma_abs;
MCNUM sigma_inc = mccMono6_sigma_inc;
MCNUM aa = mccMono6_aa;
MCNUM bb = mccMono6_bb;
MCNUM cc = mccMono6_cc;
MCNUM order = mccMono6_order;
MCNUM RX = mccMono6_RX;
MCNUM RY = mccMono6_RY;
MCNUM powder = mccMono6_powder;
MCNUM PG = mccMono6_PG;
MCNUM deltak = mccMono6_deltak;
#line 1409 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../samples/Single_crystal.comp"
{
  MPI_MASTER(
  if (hkl_info.flag_warning)
    fprintf(stderr, "Single_crystal: %s: Error message was repeated %i times with absorbed neutrons.\n",
      NAME_CURRENT_COMP, hkl_info.flag_warning);

  /* in case this instance is used in a SPLIT, we can recommend the
     optimal iteration value */
  if (hkl_info.nb_refl_count) {
    double split_iterations = (double)hkl_info.nb_reuses/hkl_info.nb_refl_count + 1;
    double split_optimal    = (double)hkl_info.nb_refl/hkl_info.nb_refl_count;
    if (split_optimal > split_iterations + 2)
      printf("Single_crystal: %s: Info: you may highly improve the computation efficiency by using\n"
        "    SPLIT %i COMPONENT %s=Single_crystal(order=1, ...)\n"
        "  in the instrument description %s.\n",
        NAME_CURRENT_COMP, (int)split_optimal, NAME_CURRENT_COMP, mcinstrument_source);
  }
  );
}
#line 16092 "./Test_Monochromators.c"
}   /* End of Mono6=Single_crystal() SETTING parameter declarations. */
#undef offdata
#undef hkl_info
#undef mosaic_AB
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[11]) fprintf(stderr, "Warning: No neutron could reach Component[11] Mono6\n");
    if (mcAbsorbProp[11]) fprintf(stderr, "Warning: %g events were removed in Component[11] Mono6=Single_crystal()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[11]);
    if (!mcNCounter[12]) fprintf(stderr, "Warning: No neutron could reach Component[12] Sphere1\n");
    if (mcAbsorbProp[12]) fprintf(stderr, "Warning: %g events were removed in Component[12] Sphere1=PSD_monitor_4PI()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[12]);
    if (!mcNCounter[13]) fprintf(stderr, "Warning: No neutron could reach Component[13] lam1\n");
    if (mcAbsorbProp[13]) fprintf(stderr, "Warning: %g events were removed in Component[13] lam1=L_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[13]);
  /* User FINALLY code for component 'psd1'. */
  SIG_MESSAGE("psd1 (Finally)");
#define mccompcurname  psd1
#define mccompcurtype  PSD_monitor
#define mccompcurindex 14
#define PSD_N mccpsd1_PSD_N
#define PSD_p mccpsd1_PSD_p
#define PSD_p2 mccpsd1_PSD_p2
{   /* Declarations of psd1=PSD_monitor() SETTING parameters. */
int nx = mccpsd1_nx;
int ny = mccpsd1_ny;
char* filename = mccpsd1_filename;
MCNUM xmin = mccpsd1_xmin;
MCNUM xmax = mccpsd1_xmax;
MCNUM ymin = mccpsd1_ymin;
MCNUM ymax = mccpsd1_ymax;
MCNUM xwidth = mccpsd1_xwidth;
MCNUM yheight = mccpsd1_yheight;
MCNUM restore_neutron = mccpsd1_restore_neutron;
#line 122 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
{
  destroy_darr2d(PSD_N);
  destroy_darr2d(PSD_p);
  destroy_darr2d(PSD_p2);
}
#line 16132 "./Test_Monochromators.c"
}   /* End of psd1=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[14]) fprintf(stderr, "Warning: No neutron could reach Component[14] psd1\n");
    if (mcAbsorbProp[14]) fprintf(stderr, "Warning: %g events were removed in Component[14] psd1=PSD_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[14]);
  mcsiminfo_close(); 
} /* end finally */
#define magnify mcdis_magnify
#define line mcdis_line
#define dashed_line mcdis_dashed_line
#define multiline mcdis_multiline
#define rectangle mcdis_rectangle
#define box mcdis_box
#define circle mcdis_circle
#define cylinder mcdis_cylinder
#define sphere mcdis_sphere
void mcdisplay(void) {
  printf("MCDISPLAY: start\n");
  /* Components MCDISPLAY code. */

  /* MCDISPLAY code for component 'Origin'. */
  SIG_MESSAGE("Origin (McDisplay)");
  printf("MCDISPLAY: component %s\n", "Origin");
#define mccompcurname  Origin
#define mccompcurtype  Progress_bar
#define mccompcurindex 1
#define IntermediateCnts mccOrigin_IntermediateCnts
#define StartTime mccOrigin_StartTime
#define EndTime mccOrigin_EndTime
#define CurrentTime mccOrigin_CurrentTime
{   /* Declarations of Origin=Progress_bar() SETTING parameters. */
char* profile = mccOrigin_profile;
MCNUM percent = mccOrigin_percent;
MCNUM flag_save = mccOrigin_flag_save;
MCNUM minutes = mccOrigin_minutes;
#line 147 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../misc/Progress_bar.comp"
{
  
}
#line 16177 "./Test_Monochromators.c"
}   /* End of Origin=Progress_bar() SETTING parameter declarations. */
#undef CurrentTime
#undef EndTime
#undef StartTime
#undef IntermediateCnts
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'Source'. */
  SIG_MESSAGE("Source (McDisplay)");
  printf("MCDISPLAY: component %s\n", "Source");
#define mccompcurname  Source
#define mccompcurtype  Source_simple
#define mccompcurindex 2
#define pmul mccSource_pmul
#define square mccSource_square
#define srcArea mccSource_srcArea
{   /* Declarations of Source=Source_simple() SETTING parameters. */
MCNUM radius = mccSource_radius;
MCNUM yheight = mccSource_yheight;
MCNUM xwidth = mccSource_xwidth;
MCNUM dist = mccSource_dist;
MCNUM focus_xw = mccSource_focus_xw;
MCNUM focus_yh = mccSource_focus_yh;
MCNUM E0 = mccSource_E0;
MCNUM dE = mccSource_dE;
MCNUM lambda0 = mccSource_lambda0;
MCNUM dlambda = mccSource_dlambda;
MCNUM flux = mccSource_flux;
MCNUM gauss = mccSource_gauss;
int target_index = mccSource_target_index;
#line 171 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../sources/Source_simple.comp"
{
  if (square == 1) {
    
    rectangle("xy",0,0,0,xwidth,yheight);
  } else {
    
    circle("xy",0,0,0,radius);
  }
  if (dist) {
    dashed_line(0,0,0, -focus_xw/2+tx,-focus_yh/2+ty,tz, 4);
    dashed_line(0,0,0,  focus_xw/2+tx,-focus_yh/2+ty,tz, 4);
    dashed_line(0,0,0,  focus_xw/2+tx, focus_yh/2+ty,tz, 4);
    dashed_line(0,0,0, -focus_xw/2+tx, focus_yh/2+ty,tz, 4);
  }
}
#line 16226 "./Test_Monochromators.c"
}   /* End of Source=Source_simple() SETTING parameter declarations. */
#undef srcArea
#undef square
#undef pmul
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'lamStart'. */
  SIG_MESSAGE("lamStart (McDisplay)");
  printf("MCDISPLAY: component %s\n", "lamStart");
#define mccompcurname  lamStart
#define mccompcurtype  L_monitor
#define mccompcurindex 3
#define nL mcclamStart_nL
#define L_N mcclamStart_L_N
#define L_p mcclamStart_L_p
#define L_p2 mcclamStart_L_p2
{   /* Declarations of lamStart=L_monitor() SETTING parameters. */
char* filename = mcclamStart_filename;
MCNUM xmin = mcclamStart_xmin;
MCNUM xmax = mcclamStart_xmax;
MCNUM ymin = mcclamStart_ymin;
MCNUM ymax = mcclamStart_ymax;
MCNUM xwidth = mcclamStart_xwidth;
MCNUM yheight = mcclamStart_yheight;
MCNUM Lmin = mcclamStart_Lmin;
MCNUM Lmax = mcclamStart_Lmax;
MCNUM restore_neutron = mcclamStart_restore_neutron;
int nowritefile = mcclamStart_nowritefile;
#line 120 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/L_monitor.comp"
{
  
  multiline(5, (double)xmin, (double)ymin, 0.0,
               (double)xmax, (double)ymin, 0.0,
               (double)xmax, (double)ymax, 0.0,
               (double)xmin, (double)ymax, 0.0,
               (double)xmin, (double)ymin, 0.0);
}
#line 16266 "./Test_Monochromators.c"
}   /* End of lamStart=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'Mono_Arm'. */
  SIG_MESSAGE("Mono_Arm (McDisplay)");
  printf("MCDISPLAY: component %s\n", "Mono_Arm");
#define mccompcurname  Mono_Arm
#define mccompcurtype  Arm
#define mccompcurindex 4
#line 40 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Arm.comp"
{
  /* A bit ugly; hard-coded dimensions. */
  
  line(0,0,0,0.2,0,0);
  line(0,0,0,0,0.2,0);
  line(0,0,0,0,0,0.2);
}
#line 16290 "./Test_Monochromators.c"
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'Mono_Out'. */
  SIG_MESSAGE("Mono_Out (McDisplay)");
  printf("MCDISPLAY: component %s\n", "Mono_Out");
#define mccompcurname  Mono_Out
#define mccompcurtype  Arm
#define mccompcurindex 5
#line 40 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Arm.comp"
{
  /* A bit ugly; hard-coded dimensions. */
  
  line(0,0,0,0.2,0,0);
  line(0,0,0,0,0.2,0);
  line(0,0,0,0,0,0.2);
}
#line 16309 "./Test_Monochromators.c"
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'Mono1'. */
  SIG_MESSAGE("Mono1 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "Mono1");
#define mccompcurname  Mono1
#define mccompcurtype  Monochromator_flat
#define mccompcurindex 6
#define mos_rms_y mccMono1_mos_rms_y
#define mos_rms_z mccMono1_mos_rms_z
#define mos_rms_max mccMono1_mos_rms_max
#define mono_Q mccMono1_mono_Q
{   /* Declarations of Mono1=Monochromator_flat() SETTING parameters. */
MCNUM zmin = mccMono1_zmin;
MCNUM zmax = mccMono1_zmax;
MCNUM ymin = mccMono1_ymin;
MCNUM ymax = mccMono1_ymax;
MCNUM zwidth = mccMono1_zwidth;
MCNUM yheight = mccMono1_yheight;
MCNUM mosaich = mccMono1_mosaich;
MCNUM mosaicv = mccMono1_mosaicv;
MCNUM r0 = mccMono1_r0;
MCNUM Q = mccMono1_Q;
MCNUM DM = mccMono1_DM;
#line 254 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Monochromator_flat.comp"
{
  
  multiline(5, 0.0, (double)ymin, (double)zmin,
               0.0, (double)ymax, (double)zmin,
               0.0, (double)ymax, (double)zmax,
               0.0, (double)ymin, (double)zmax,
               0.0, (double)ymin, (double)zmin);
}
#line 16345 "./Test_Monochromators.c"
}   /* End of Mono1=Monochromator_flat() SETTING parameter declarations. */
#undef mono_Q
#undef mos_rms_max
#undef mos_rms_z
#undef mos_rms_y
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'Mono2'. */
  SIG_MESSAGE("Mono2 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "Mono2");
#define mccompcurname  Mono2
#define mccompcurtype  Monochromator_pol
#define mccompcurindex 7
#define mos_rms mccMono2_mos_rms
#define d_rms mccMono2_d_rms
#define mono_Q mccMono2_mono_Q
#define FN mccMono2_FN
#define FM mccMono2_FM
{   /* Declarations of Mono2=Monochromator_pol() SETTING parameters. */
MCNUM zwidth = mccMono2_zwidth;
MCNUM yheight = mccMono2_yheight;
MCNUM mosaic = mccMono2_mosaic;
MCNUM dspread = mccMono2_dspread;
MCNUM Q = mccMono2_Q;
MCNUM DM = mccMono2_DM;
MCNUM pThreshold = mccMono2_pThreshold;
MCNUM Rup = mccMono2_Rup;
MCNUM Rdown = mccMono2_Rdown;
int debug = mccMono2_debug;
#line 199 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Monochromator_pol.comp"
{
  
  rectangle("yz", 0, 0, 0, zwidth, yheight);
}
#line 16382 "./Test_Monochromators.c"
}   /* End of Mono2=Monochromator_pol() SETTING parameter declarations. */
#undef FM
#undef FN
#undef mono_Q
#undef d_rms
#undef mos_rms
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'Mono3'. */
  SIG_MESSAGE("Mono3 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "Mono3");
#define mccompcurname  Mono3
#define mccompcurtype  Monochromator_pol
#define mccompcurindex 8
#define mos_rms mccMono3_mos_rms
#define d_rms mccMono3_d_rms
#define mono_Q mccMono3_mono_Q
#define FN mccMono3_FN
#define FM mccMono3_FM
{   /* Declarations of Mono3=Monochromator_pol() SETTING parameters. */
MCNUM zwidth = mccMono3_zwidth;
MCNUM yheight = mccMono3_yheight;
MCNUM mosaic = mccMono3_mosaic;
MCNUM dspread = mccMono3_dspread;
MCNUM Q = mccMono3_Q;
MCNUM DM = mccMono3_DM;
MCNUM pThreshold = mccMono3_pThreshold;
MCNUM Rup = mccMono3_Rup;
MCNUM Rdown = mccMono3_Rdown;
int debug = mccMono3_debug;
#line 199 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Monochromator_pol.comp"
{
  
  rectangle("yz", 0, 0, 0, zwidth, yheight);
}
#line 16420 "./Test_Monochromators.c"
}   /* End of Mono3=Monochromator_pol() SETTING parameter declarations. */
#undef FM
#undef FN
#undef mono_Q
#undef d_rms
#undef mos_rms
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'Mono4'. */
  SIG_MESSAGE("Mono4 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "Mono4");
#define mccompcurname  Mono4
#define mccompcurtype  Monochromator_curved
#define mccompcurindex 9
#define mos_rms_y mccMono4_mos_rms_y
#define mos_rms_z mccMono4_mos_rms_z
#define mos_rms_max mccMono4_mos_rms_max
#define mono_Q mccMono4_mono_Q
#define SlabWidth mccMono4_SlabWidth
#define SlabHeight mccMono4_SlabHeight
#define rTable mccMono4_rTable
#define tTable mccMono4_tTable
#define row mccMono4_row
#define col mccMono4_col
#define tiltH mccMono4_tiltH
#define tiltV mccMono4_tiltV
{   /* Declarations of Mono4=Monochromator_curved() SETTING parameters. */
char* reflect = mccMono4_reflect;
char* transmit = mccMono4_transmit;
MCNUM zwidth = mccMono4_zwidth;
MCNUM yheight = mccMono4_yheight;
MCNUM gap = mccMono4_gap;
MCNUM NH = mccMono4_NH;
MCNUM NV = mccMono4_NV;
MCNUM mosaich = mccMono4_mosaich;
MCNUM mosaicv = mccMono4_mosaicv;
MCNUM r0 = mccMono4_r0;
MCNUM t0 = mccMono4_t0;
MCNUM Q = mccMono4_Q;
MCNUM RV = mccMono4_RV;
MCNUM RH = mccMono4_RH;
MCNUM DM = mccMono4_DM;
MCNUM mosaic = mccMono4_mosaic;
MCNUM width = mccMono4_width;
MCNUM height = mccMono4_height;
MCNUM verbose = mccMono4_verbose;
MCNUM order = mccMono4_order;
#line 468 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Monochromator_curved.comp"
{
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
}
#line 16500 "./Test_Monochromators.c"
}   /* End of Mono4=Monochromator_curved() SETTING parameter declarations. */
#undef tiltV
#undef tiltH
#undef col
#undef row
#undef tTable
#undef rTable
#undef SlabHeight
#undef SlabWidth
#undef mono_Q
#undef mos_rms_max
#undef mos_rms_z
#undef mos_rms_y
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'Mono5'. */
  SIG_MESSAGE("Mono5 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "Mono5");
#define mccompcurname  Mono5
#define mccompcurtype  Monochromator_2foc
#define mccompcurindex 10
#define mos_y mccMono5_mos_y
#define mos_z mccMono5_mos_z
#define mono_Q mccMono5_mono_Q
#define SlabWidth mccMono5_SlabWidth
#define SlabHeight mccMono5_SlabHeight
#define rTable mccMono5_rTable
{   /* Declarations of Mono5=Monochromator_2foc() SETTING parameters. */
char* reflect = mccMono5_reflect;
MCNUM zwidth = mccMono5_zwidth;
MCNUM yheight = mccMono5_yheight;
MCNUM gap = mccMono5_gap;
MCNUM NH = mccMono5_NH;
MCNUM NV = mccMono5_NV;
MCNUM mosaich = mccMono5_mosaich;
MCNUM mosaicv = mccMono5_mosaicv;
MCNUM r0 = mccMono5_r0;
MCNUM Q = mccMono5_Q;
MCNUM RV = mccMono5_RV;
MCNUM RH = mccMono5_RH;
MCNUM DM = mccMono5_DM;
MCNUM mosaic = mccMono5_mosaic;
MCNUM width = mccMono5_width;
MCNUM height = mccMono5_height;
MCNUM verbose = mccMono5_verbose;
#line 270 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../contrib/Monochromator_2foc.comp"
{
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
}
#line 16583 "./Test_Monochromators.c"
}   /* End of Mono5=Monochromator_2foc() SETTING parameter declarations. */
#undef rTable
#undef SlabHeight
#undef SlabWidth
#undef mono_Q
#undef mos_z
#undef mos_y
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'Mono6'. */
  SIG_MESSAGE("Mono6 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "Mono6");
#define mccompcurname  Mono6
#define mccompcurtype  Single_crystal
#define mccompcurindex 11
#define mosaic_AB mccMono6_mosaic_AB
#define hkl_info mccMono6_hkl_info
#define offdata mccMono6_offdata
{   /* Declarations of Mono6=Single_crystal() SETTING parameters. */
char* reflections = mccMono6_reflections;
char* geometry = mccMono6_geometry;
MCNUM xwidth = mccMono6_xwidth;
MCNUM yheight = mccMono6_yheight;
MCNUM zdepth = mccMono6_zdepth;
MCNUM radius = mccMono6_radius;
MCNUM delta_d_d = mccMono6_delta_d_d;
MCNUM mosaic = mccMono6_mosaic;
MCNUM mosaic_a = mccMono6_mosaic_a;
MCNUM mosaic_b = mccMono6_mosaic_b;
MCNUM mosaic_c = mccMono6_mosaic_c;
MCNUM recip_cell = mccMono6_recip_cell;
MCNUM barns = mccMono6_barns;
MCNUM ax = mccMono6_ax;
MCNUM ay = mccMono6_ay;
MCNUM az = mccMono6_az;
MCNUM bx = mccMono6_bx;
MCNUM by = mccMono6_by;
MCNUM bz = mccMono6_bz;
MCNUM cx = mccMono6_cx;
MCNUM cy = mccMono6_cy;
MCNUM cz = mccMono6_cz;
MCNUM p_transmit = mccMono6_p_transmit;
MCNUM sigma_abs = mccMono6_sigma_abs;
MCNUM sigma_inc = mccMono6_sigma_inc;
MCNUM aa = mccMono6_aa;
MCNUM bb = mccMono6_bb;
MCNUM cc = mccMono6_cc;
MCNUM order = mccMono6_order;
MCNUM RX = mccMono6_RX;
MCNUM RY = mccMono6_RY;
MCNUM powder = mccMono6_powder;
MCNUM PG = mccMono6_PG;
MCNUM deltak = mccMono6_deltak;
#line 1430 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../samples/Single_crystal.comp"
{
  
  if (hkl_info.shape == 0) {        /* cylinder */
    circle("xz", 0,  yheight/2.0, 0, radius);
    circle("xz", 0, -yheight/2.0, 0, radius);
    line(-radius, -yheight/2.0, 0, -radius, +yheight/2.0, 0);
    line(+radius, -yheight/2.0, 0, +radius, +yheight/2.0, 0);
    line(0, -yheight/2.0, -radius, 0, +yheight/2.0, -radius);
    line(0, -yheight/2.0, +radius, 0, +yheight/2.0, +radius);
  }
  else if (hkl_info.shape == 1) {         /* box */
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
  else if (hkl_info.shape == 2) {        /* sphere */
    circle("xy", 0,  0.0, 0, radius);
    circle("xz", 0,  0.0, 0, radius);
    circle("yz", 0,  0.0, 0, radius);
  }
  else if (hkl_info.shape == 3) {        /* OFF file */
    off_display(offdata);
  }
}
#line 16681 "./Test_Monochromators.c"
}   /* End of Mono6=Single_crystal() SETTING parameter declarations. */
#undef offdata
#undef hkl_info
#undef mosaic_AB
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'Sphere1'. */
  SIG_MESSAGE("Sphere1 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "Sphere1");
#define mccompcurname  Sphere1
#define mccompcurtype  PSD_monitor_4PI
#define mccompcurindex 12
#define nx mccSphere1_nx
#define ny mccSphere1_ny
#define PSD_N mccSphere1_PSD_N
#define PSD_p mccSphere1_PSD_p
#define PSD_p2 mccSphere1_PSD_p2
{   /* Declarations of Sphere1=PSD_monitor_4PI() SETTING parameters. */
char* filename = mccSphere1_filename;
MCNUM radius = mccSphere1_radius;
MCNUM restore_neutron = mccSphere1_restore_neutron;
int nowritefile = mccSphere1_nowritefile;
#line 121 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor_4PI.comp"
{
  
  circle("xy",0,0,0,radius);
  circle("xz",0,0,0,radius);
  circle("yz",0,0,0,radius);
}
#line 16713 "./Test_Monochromators.c"
}   /* End of Sphere1=PSD_monitor_4PI() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef ny
#undef nx
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'lam1'. */
  SIG_MESSAGE("lam1 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "lam1");
#define mccompcurname  lam1
#define mccompcurtype  L_monitor
#define mccompcurindex 13
#define nL mcclam1_nL
#define L_N mcclam1_L_N
#define L_p mcclam1_L_p
#define L_p2 mcclam1_L_p2
{   /* Declarations of lam1=L_monitor() SETTING parameters. */
char* filename = mcclam1_filename;
MCNUM xmin = mcclam1_xmin;
MCNUM xmax = mcclam1_xmax;
MCNUM ymin = mcclam1_ymin;
MCNUM ymax = mcclam1_ymax;
MCNUM xwidth = mcclam1_xwidth;
MCNUM yheight = mcclam1_yheight;
MCNUM Lmin = mcclam1_Lmin;
MCNUM Lmax = mcclam1_Lmax;
MCNUM restore_neutron = mcclam1_restore_neutron;
int nowritefile = mcclam1_nowritefile;
#line 120 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/L_monitor.comp"
{
  
  multiline(5, (double)xmin, (double)ymin, 0.0,
               (double)xmax, (double)ymin, 0.0,
               (double)xmax, (double)ymax, 0.0,
               (double)xmin, (double)ymax, 0.0,
               (double)xmin, (double)ymin, 0.0);
}
#line 16755 "./Test_Monochromators.c"
}   /* End of lam1=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'psd1'. */
  SIG_MESSAGE("psd1 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "psd1");
#define mccompcurname  psd1
#define mccompcurtype  PSD_monitor
#define mccompcurindex 14
#define PSD_N mccpsd1_PSD_N
#define PSD_p mccpsd1_PSD_p
#define PSD_p2 mccpsd1_PSD_p2
{   /* Declarations of psd1=PSD_monitor() SETTING parameters. */
int nx = mccpsd1_nx;
int ny = mccpsd1_ny;
char* filename = mccpsd1_filename;
MCNUM xmin = mccpsd1_xmin;
MCNUM xmax = mccpsd1_xmax;
MCNUM ymin = mccpsd1_ymin;
MCNUM ymax = mccpsd1_ymax;
MCNUM xwidth = mccpsd1_xwidth;
MCNUM yheight = mccpsd1_yheight;
MCNUM restore_neutron = mccpsd1_restore_neutron;
#line 129 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
{
  multiline(5,
    (double)xmin, (double)ymin, 0.0,
    (double)xmax, (double)ymin, 0.0,
    (double)xmax, (double)ymax, 0.0,
    (double)xmin, (double)ymax, 0.0,
    (double)xmin, (double)ymin, 0.0);
}
#line 16794 "./Test_Monochromators.c"
}   /* End of psd1=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  printf("MCDISPLAY: end\n");
} /* end display */
#undef magnify
#undef line
#undef dashed_line
#undef multiline
#undef rectangle
#undef box
#undef circle
#undef cylinder
#undef sphere
/* end of generated C code ./Test_Monochromators.c */
