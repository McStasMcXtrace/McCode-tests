/* Automatically generated file. Do not edit. 
 * Format:     ANSI C source code
 * Creator:    McStas <http://www.mcstas.org>
 * Instrument: /zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr (ESS_2001_bispectral)
 * Date:       Tue Feb 25 20:22:45 2020
 * File:       ./ESS_2001_bispectral.c
 * Compile:    cc -o ESS_2001_bispectral.out ./ESS_2001_bispectral.c 
 * CFLAGS=
 */


#define MCCODE_STRING "McStas 2.5 - Feb. 24, 2020"
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
#define MCCODE_STRING "McStas 2.5 - Feb. 24, 2020"
#endif

#ifndef MCCODE_DATE
#define MCCODE_DATE "Feb. 24, 2020"
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

#line 712 "./ESS_2001_bispectral.c"

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

#line 945 "./ESS_2001_bispectral.c"

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

#line 4977 "./ESS_2001_bispectral.c"

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

#line 5337 "./ESS_2001_bispectral.c"
#ifdef MC_TRACE_ENABLED
int mctraceenabled = 1;
#else
int mctraceenabled = 0;
#endif
#define MCSTAS "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../"
int mcdefaultmain = 1;
char mcinstrument_name[] = "ESS_2001_bispectral";
char mcinstrument_source[] = "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr";
char *mcinstrument_exe=NULL; /* will be set to argv[0] in main */
int main(int argc, char *argv[]){return mccode_main(argc, argv);}
void mcinit(void);
void mcraytrace(void);
void mcsave(FILE *);
void mcfinally(void);
void mcdisplay(void);

/* Shared user declarations for all components 'ESS_moderator_long_2001'. */
#line 96 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../obsolete/ESS_moderator_long_2001.comp"
double Mezei_M_fct(double l, double temp)
  {
    double a=949.0/temp;
    return 2*a*a*exp(-a/(l*l))/(l*l*l*l*l);
  }

  double Mezei_F_fct(double t, double tau, int n)
  {
    return (exp(-t/tau)-exp(-n*t/tau))*n/(n-1)/tau;
  }
#line 5367 "./ESS_2001_bispectral.c"

/* Shared user declarations for all components 'Mirror_Curved_Bispectral'. */
#line 55 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../contrib/Mirror_Curved_Bispectral.comp"
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

#line 6814 "./ESS_2001_bispectral.c"

/* Shared user declarations for all components 'Mirror_Elliptic_Bispectral'. */
#line 57 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../contrib/Mirror_Elliptic_Bispectral.comp"

#line 6819 "./ESS_2001_bispectral.c"

/* Instrument parameters. */
int mcipman_lam;
int mcipthermal;
MCNUM mcipcold_Lam_min;
MCNUM mcipcold_Lam_max;
MCNUM mcipthermal_Lam_min;
MCNUM mcipthermal_Lam_max;
MCNUM mcipthermal_hotspot_fac;
MCNUM mcipcold_hotspot_fac;
int mcipuse_full_guide_flag;
MCNUM mcipguide_start;
MCNUM mcipLength;
MCNUM mcipfocus_start_w;
MCNUM mcipfocus_end_w;
MCNUM mcipsmallaxis_w;
MCNUM mcipfocus_start_h;
MCNUM mcipfocus_end_h;
MCNUM mcipsmallaxis_h;
MCNUM mcipmaxdiv;
MCNUM mcipcoldperthermal;
int mcipmirror_type;
int mcipmirror_coating_type;
MCNUM mcipmirror_offset;
MCNUM mciptheta1;
MCNUM mciptheta2;
MCNUM mciptheta3;
MCNUM mcipm_mirror;
MCNUM mciph_mirror;
MCNUM mcipPulse_width;
MCNUM mcipfrequency;
MCNUM mcipgravity;
MCNUM mcipsubstrate_thickness;
MCNUM mcipcoating_thickness;
MCNUM mcipm1;
MCNUM mcipm2;

#define mcNUMIPAR 34
int mcnumipar = 34;
struct mcinputtable_struct mcinputtable[mcNUMIPAR+1] = {
  "man_lam", &mcipman_lam, instr_type_int, "1", 
  "thermal", &mcipthermal, instr_type_int, "2", 
  "cold_Lam_min", &mcipcold_Lam_min, instr_type_double, "0.1", 
  "cold_Lam_max", &mcipcold_Lam_max, instr_type_double, "10", 
  "thermal_Lam_min", &mcipthermal_Lam_min, instr_type_double, "0.1", 
  "thermal_Lam_max", &mcipthermal_Lam_max, instr_type_double, "10", 
  "thermal_hotspot_fac", &mcipthermal_hotspot_fac, instr_type_double, "1.0", 
  "cold_hotspot_fac", &mcipcold_hotspot_fac, instr_type_double, "1.0", 
  "use_full_guide_flag", &mcipuse_full_guide_flag, instr_type_int, "0", 
  "guide_start", &mcipguide_start, instr_type_double, "2", 
  "Length", &mcipLength, instr_type_double, "150", 
  "focus_start_w", &mcipfocus_start_w, instr_type_double, "-1.2620", 
  "focus_end_w", &mcipfocus_end_w, instr_type_double, "150.3187", 
  "smallaxis_w", &mcipsmallaxis_w, instr_type_double, "0.2550", 
  "focus_start_h", &mcipfocus_start_h, instr_type_double, "-2.1415", 
  "focus_end_h", &mcipfocus_end_h, instr_type_double, "150.0778", 
  "smallaxis_h", &mcipsmallaxis_h, instr_type_double, "0.3847", 
  "maxdiv", &mcipmaxdiv, instr_type_double, "2.0", 
  "coldperthermal", &mcipcoldperthermal, instr_type_double, "30", 
  "mirror_type", &mcipmirror_type, instr_type_int, "2", 
  "mirror_coating_type", &mcipmirror_coating_type, instr_type_int, "1", 
  "mirror_offset", &mcipmirror_offset, instr_type_double, "0", 
  "theta1", &mciptheta1, instr_type_double, "1.25", 
  "theta2", &mciptheta2, instr_type_double, "1.25", 
  "theta3", &mciptheta3, instr_type_double, "1.25", 
  "m_mirror", &mcipm_mirror, instr_type_double, "5", 
  "h_mirror", &mciph_mirror, instr_type_double, "0.15", 
  "Pulse_width", &mcipPulse_width, instr_type_double, "0.00286", 
  "frequency", &mcipfrequency, instr_type_double, "14.0", 
  "gravity", &mcipgravity, instr_type_double, "-9.81", 
  "substrate_thickness", &mcipsubstrate_thickness, instr_type_double, "0.0005", 
  "coating_thickness", &mcipcoating_thickness, instr_type_double, "10e-6", 
  "m1", &mcipm1, instr_type_double, "5", 
  "m2", &mcipm2, instr_type_double, "4", 
  NULL, NULL, instr_type_double, ""
};

/* User declarations from instrument definition. */
#define mccompcurname  ESS_2001_bispectral
#define mccompcurtype  INSTRUMENT
#define mccompcurindex 0
#define mcposaESS_2001_bispectral coords_set(0,0,0)
#define man_lam mcipman_lam
#define thermal mcipthermal
#define cold_Lam_min mcipcold_Lam_min
#define cold_Lam_max mcipcold_Lam_max
#define thermal_Lam_min mcipthermal_Lam_min
#define thermal_Lam_max mcipthermal_Lam_max
#define thermal_hotspot_fac mcipthermal_hotspot_fac
#define cold_hotspot_fac mcipcold_hotspot_fac
#define use_full_guide_flag mcipuse_full_guide_flag
#define guide_start mcipguide_start
#define Length mcipLength
#define focus_start_w mcipfocus_start_w
#define focus_end_w mcipfocus_end_w
#define smallaxis_w mcipsmallaxis_w
#define focus_start_h mcipfocus_start_h
#define focus_end_h mcipfocus_end_h
#define smallaxis_h mcipsmallaxis_h
#define maxdiv mcipmaxdiv
#define coldperthermal mcipcoldperthermal
#define mirror_type mcipmirror_type
#define mirror_coating_type mcipmirror_coating_type
#define mirror_offset mcipmirror_offset
#define theta1 mciptheta1
#define theta2 mciptheta2
#define theta3 mciptheta3
#define m_mirror mcipm_mirror
#define h_mirror mciph_mirror
#define Pulse_width mcipPulse_width
#define frequency mcipfrequency
#define gravity mcipgravity
#define substrate_thickness mcipsubstrate_thickness
#define coating_thickness mcipcoating_thickness
#define m1 mcipm1
#define m2 mcipm2
#line 67 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
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


#line 7108 "./ESS_2001_bispectral.c"
#undef m2
#undef m1
#undef coating_thickness
#undef substrate_thickness
#undef gravity
#undef frequency
#undef Pulse_width
#undef h_mirror
#undef m_mirror
#undef theta3
#undef theta2
#undef theta1
#undef mirror_offset
#undef mirror_coating_type
#undef mirror_type
#undef coldperthermal
#undef maxdiv
#undef smallaxis_h
#undef focus_end_h
#undef focus_start_h
#undef smallaxis_w
#undef focus_end_w
#undef focus_start_w
#undef Length
#undef guide_start
#undef use_full_guide_flag
#undef cold_hotspot_fac
#undef thermal_hotspot_fac
#undef thermal_Lam_max
#undef thermal_Lam_min
#undef cold_Lam_max
#undef cold_Lam_min
#undef thermal
#undef man_lam
#undef mcposaESS_2001_bispectral
#undef mccompcurindex
#undef mccompcurtype
#undef mccompcurname

/* neutron state table at each component input (local coords) */
/* [x, y, z, vx, vy, vz, t, sx, sy, sz, p] */
MCNUM mccomp_storein[11*28];
/* Components position table (absolute and relative coords) */
Coords mccomp_posa[28];
Coords mccomp_posr[28];
/* Counter for each comp to check for inactive ones */
MCNUM  mcNCounter[28];
MCNUM  mcPCounter[28];
MCNUM  mcP2Counter[28];
#define mcNUMCOMP 27 /* number of components */
/* Counter for PROP ABSORB */
MCNUM  mcAbsorbProp[28];
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

/* Setting parameters for component 'cold_source' [6]. */
MCNUM mcccold_source_size;
MCNUM mcccold_source_l_low;
MCNUM mcccold_source_l_high;
MCNUM mcccold_source_dist;
MCNUM mcccold_source_xw;
MCNUM mcccold_source_yh;
MCNUM mcccold_source_freq;
MCNUM mcccold_source_T;
MCNUM mcccold_source_tau;
MCNUM mcccold_source_tau1;
MCNUM mcccold_source_tau2;
MCNUM mcccold_source_d;
MCNUM mcccold_source_n;
MCNUM mcccold_source_n2;
MCNUM mcccold_source_chi2;
MCNUM mcccold_source_I0;
MCNUM mcccold_source_I2;
MCNUM mcccold_source_branch1;
MCNUM mcccold_source_branch2;
MCNUM mcccold_source_branch_tail;
MCNUM mcccold_source_twopulses;
int mcccold_source_target_index;

/* Setting parameters for component 'thermal_source' [7]. */
MCNUM mccthermal_source_size;
MCNUM mccthermal_source_l_low;
MCNUM mccthermal_source_l_high;
MCNUM mccthermal_source_dist;
MCNUM mccthermal_source_xw;
MCNUM mccthermal_source_yh;
MCNUM mccthermal_source_freq;
MCNUM mccthermal_source_T;
MCNUM mccthermal_source_tau;
MCNUM mccthermal_source_tau1;
MCNUM mccthermal_source_tau2;
MCNUM mccthermal_source_d;
MCNUM mccthermal_source_n;
MCNUM mccthermal_source_n2;
MCNUM mccthermal_source_chi2;
MCNUM mccthermal_source_I0;
MCNUM mccthermal_source_I2;
MCNUM mccthermal_source_branch1;
MCNUM mccthermal_source_branch2;
MCNUM mccthermal_source_branch_tail;
MCNUM mccthermal_source_twopulses;
int mccthermal_source_target_index;

/* Definition parameters for component 'mirror_full_center' [10]. */
#define mccmirror_full_center_reflect 0 /* declared as a string. May produce warnings at compile */
/* Setting parameters for component 'mirror_full_center' [10]. */
MCNUM mccmirror_full_center_focus_s;
MCNUM mccmirror_full_center_focus_e;
MCNUM mccmirror_full_center_mirror_start;
MCNUM mccmirror_full_center_guide_start;
MCNUM mccmirror_full_center_yheight;
MCNUM mccmirror_full_center_smallaxis;
MCNUM mccmirror_full_center_length;
MCNUM mccmirror_full_center_m;
MCNUM mccmirror_full_center_transmit;
MCNUM mccmirror_full_center_substrate_thickness;
MCNUM mccmirror_full_center_coating_thickness;
MCNUM mccmirror_full_center_theta_1;
MCNUM mccmirror_full_center_theta_2;
MCNUM mccmirror_full_center_theta_3;

/* Definition parameters for component 'guide_right' [12]. */
#define mccguide_right_reflect 0 /* declared as a string. May produce warnings at compile */
/* Setting parameters for component 'guide_right' [12]. */
MCNUM mccguide_right_focus_start_w;
MCNUM mccguide_right_focus_end_w;
MCNUM mccguide_right_focus_start_h;
MCNUM mccguide_right_focus_end_h;
MCNUM mccguide_right_mirror_start;
MCNUM mccguide_right_m;
MCNUM mccguide_right_smallaxis_w;
MCNUM mccguide_right_smallaxis_h;
MCNUM mccguide_right_length;
MCNUM mccguide_right_transmit;
MCNUM mccguide_right_substrate_thickness;
MCNUM mccguide_right_coating_thickness;

/* Definition parameters for component 'guide_bottom' [14]. */
#define mccguide_bottom_reflect 0 /* declared as a string. May produce warnings at compile */
/* Setting parameters for component 'guide_bottom' [14]. */
MCNUM mccguide_bottom_focus_start_w;
MCNUM mccguide_bottom_focus_end_w;
MCNUM mccguide_bottom_focus_start_h;
MCNUM mccguide_bottom_focus_end_h;
MCNUM mccguide_bottom_mirror_start;
MCNUM mccguide_bottom_m;
MCNUM mccguide_bottom_smallaxis_w;
MCNUM mccguide_bottom_smallaxis_h;
MCNUM mccguide_bottom_length;
MCNUM mccguide_bottom_transmit;
MCNUM mccguide_bottom_substrate_thickness;
MCNUM mccguide_bottom_coating_thickness;

/* Definition parameters for component 'cold_lambda_guidestart' [16]. */
#define mcccold_lambda_guidestart_nL 100
/* Setting parameters for component 'cold_lambda_guidestart' [16]. */
char mcccold_lambda_guidestart_filename[16384];
MCNUM mcccold_lambda_guidestart_xmin;
MCNUM mcccold_lambda_guidestart_xmax;
MCNUM mcccold_lambda_guidestart_ymin;
MCNUM mcccold_lambda_guidestart_ymax;
MCNUM mcccold_lambda_guidestart_xwidth;
MCNUM mcccold_lambda_guidestart_yheight;
MCNUM mcccold_lambda_guidestart_Lmin;
MCNUM mcccold_lambda_guidestart_Lmax;
MCNUM mcccold_lambda_guidestart_restore_neutron;
int mcccold_lambda_guidestart_nowritefile;

/* Definition parameters for component 'thermal_lambda_guidestart' [17]. */
#define mccthermal_lambda_guidestart_nL 100
/* Setting parameters for component 'thermal_lambda_guidestart' [17]. */
char mccthermal_lambda_guidestart_filename[16384];
MCNUM mccthermal_lambda_guidestart_xmin;
MCNUM mccthermal_lambda_guidestart_xmax;
MCNUM mccthermal_lambda_guidestart_ymin;
MCNUM mccthermal_lambda_guidestart_ymax;
MCNUM mccthermal_lambda_guidestart_xwidth;
MCNUM mccthermal_lambda_guidestart_yheight;
MCNUM mccthermal_lambda_guidestart_Lmin;
MCNUM mccthermal_lambda_guidestart_Lmax;
MCNUM mccthermal_lambda_guidestart_restore_neutron;
int mccthermal_lambda_guidestart_nowritefile;

/* Definition parameters for component 'lambda_guidestart' [18]. */
#define mcclambda_guidestart_nL 100
/* Setting parameters for component 'lambda_guidestart' [18]. */
char mcclambda_guidestart_filename[16384];
MCNUM mcclambda_guidestart_xmin;
MCNUM mcclambda_guidestart_xmax;
MCNUM mcclambda_guidestart_ymin;
MCNUM mcclambda_guidestart_ymax;
MCNUM mcclambda_guidestart_xwidth;
MCNUM mcclambda_guidestart_yheight;
MCNUM mcclambda_guidestart_Lmin;
MCNUM mcclambda_guidestart_Lmax;
MCNUM mcclambda_guidestart_restore_neutron;
int mcclambda_guidestart_nowritefile;

/* Definition parameters for component 'guide_top' [19]. */
#define mccguide_top_reflect 0 /* declared as a string. May produce warnings at compile */
/* Setting parameters for component 'guide_top' [19]. */
MCNUM mccguide_top_focus_start_w;
MCNUM mccguide_top_focus_end_w;
MCNUM mccguide_top_focus_start_h;
MCNUM mccguide_top_focus_end_h;
MCNUM mccguide_top_mirror_start;
MCNUM mccguide_top_m;
MCNUM mccguide_top_smallaxis_w;
MCNUM mccguide_top_smallaxis_h;
MCNUM mccguide_top_length;
MCNUM mccguide_top_transmit;
MCNUM mccguide_top_substrate_thickness;
MCNUM mccguide_top_coating_thickness;

/* Definition parameters for component 'guide_Left' [21]. */
#define mccguide_Left_reflect 0 /* declared as a string. May produce warnings at compile */
/* Setting parameters for component 'guide_Left' [21]. */
MCNUM mccguide_Left_focus_start_w;
MCNUM mccguide_Left_focus_end_w;
MCNUM mccguide_Left_focus_start_h;
MCNUM mccguide_Left_focus_end_h;
MCNUM mccguide_Left_mirror_start;
MCNUM mccguide_Left_m;
MCNUM mccguide_Left_smallaxis_w;
MCNUM mccguide_Left_smallaxis_h;
MCNUM mccguide_Left_length;
MCNUM mccguide_Left_transmit;
MCNUM mccguide_Left_substrate_thickness;
MCNUM mccguide_Left_coating_thickness;

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
#line 7373 "./ESS_2001_bispectral.c"
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

/* User declarations for component 'ArmForGuideRight' [2]. */
#define mccompcurname  ArmForGuideRight
#define mccompcurtype  Arm
#define mccompcurindex 2
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'ArmForGuideBottom' [3]. */
#define mccompcurname  ArmForGuideBottom
#define mccompcurtype  Arm
#define mccompcurindex 3
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'ArmForGuideTop' [4]. */
#define mccompcurname  ArmForGuideTop
#define mccompcurtype  Arm
#define mccompcurindex 4
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'ArmForGuideLeft' [5]. */
#define mccompcurname  ArmForGuideLeft
#define mccompcurtype  Arm
#define mccompcurindex 5
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'cold_source' [6]. */
#define mccompcurname  cold_source
#define mccompcurtype  ESS_moderator_long_2001
#define mccompcurindex 6
#define l_range mcccold_source_l_range
#define w_mult mcccold_source_w_mult
#define size mcccold_source_size
#define l_low mcccold_source_l_low
#define l_high mcccold_source_l_high
#define dist mcccold_source_dist
#define xw mcccold_source_xw
#define yh mcccold_source_yh
#define freq mcccold_source_freq
#define T mcccold_source_T
#define tau mcccold_source_tau
#define tau1 mcccold_source_tau1
#define tau2 mcccold_source_tau2
#define d mcccold_source_d
#define n mcccold_source_n
#define n2 mcccold_source_n2
#define chi2 mcccold_source_chi2
#define I0 mcccold_source_I0
#define I2 mcccold_source_I2
#define branch1 mcccold_source_branch1
#define branch2 mcccold_source_branch2
#define branch_tail mcccold_source_branch_tail
#define twopulses mcccold_source_twopulses
#define target_index mcccold_source_target_index
#line 110 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../obsolete/ESS_moderator_long_2001.comp"
  double l_range, w_mult, branchframe, tx, ty, tz;
#line 7448 "./ESS_2001_bispectral.c"
#undef target_index
#undef twopulses
#undef branch_tail
#undef branch2
#undef branch1
#undef I2
#undef I0
#undef chi2
#undef n2
#undef n
#undef d
#undef tau2
#undef tau1
#undef tau
#undef T
#undef freq
#undef yh
#undef xw
#undef dist
#undef l_high
#undef l_low
#undef size
#undef w_mult
#undef l_range
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'thermal_source' [7]. */
#define mccompcurname  thermal_source
#define mccompcurtype  ESS_moderator_long_2001
#define mccompcurindex 7
#define l_range mccthermal_source_l_range
#define w_mult mccthermal_source_w_mult
#define size mccthermal_source_size
#define l_low mccthermal_source_l_low
#define l_high mccthermal_source_l_high
#define dist mccthermal_source_dist
#define xw mccthermal_source_xw
#define yh mccthermal_source_yh
#define freq mccthermal_source_freq
#define T mccthermal_source_T
#define tau mccthermal_source_tau
#define tau1 mccthermal_source_tau1
#define tau2 mccthermal_source_tau2
#define d mccthermal_source_d
#define n mccthermal_source_n
#define n2 mccthermal_source_n2
#define chi2 mccthermal_source_chi2
#define I0 mccthermal_source_I0
#define I2 mccthermal_source_I2
#define branch1 mccthermal_source_branch1
#define branch2 mccthermal_source_branch2
#define branch_tail mccthermal_source_branch_tail
#define twopulses mccthermal_source_twopulses
#define target_index mccthermal_source_target_index
#line 110 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../obsolete/ESS_moderator_long_2001.comp"
  double l_range, w_mult, branchframe, tx, ty, tz;
#line 7507 "./ESS_2001_bispectral.c"
#undef target_index
#undef twopulses
#undef branch_tail
#undef branch2
#undef branch1
#undef I2
#undef I0
#undef chi2
#undef n2
#undef n
#undef d
#undef tau2
#undef tau1
#undef tau
#undef T
#undef freq
#undef yh
#undef xw
#undef dist
#undef l_high
#undef l_low
#undef size
#undef w_mult
#undef l_range
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'ColdFocus' [8]. */
#define mccompcurname  ColdFocus
#define mccompcurtype  Arm
#define mccompcurindex 8
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'ArmMidOne' [9]. */
#define mccompcurname  ArmMidOne
#define mccompcurtype  Arm
#define mccompcurindex 9
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'mirror_full_center' [10]. */
#define mccompcurname  mirror_full_center
#define mccompcurtype  Mirror_Curved_Bispectral
#define mccompcurindex 10
#define reflect mccmirror_full_center_reflect
#define pTable mccmirror_full_center_pTable
#define focus_s mccmirror_full_center_focus_s
#define focus_e mccmirror_full_center_focus_e
#define mirror_start mccmirror_full_center_mirror_start
#define guide_start mccmirror_full_center_guide_start
#define yheight mccmirror_full_center_yheight
#define smallaxis mccmirror_full_center_smallaxis
#define length mccmirror_full_center_length
#define m mccmirror_full_center_m
#define transmit mccmirror_full_center_transmit
#define substrate_thickness mccmirror_full_center_substrate_thickness
#define coating_thickness mccmirror_full_center_coating_thickness
#define theta_1 mccmirror_full_center_theta_1
#define theta_2 mccmirror_full_center_theta_2
#define theta_3 mccmirror_full_center_theta_3
#line 60 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../contrib/Mirror_Curved_Bispectral.comp"
t_Table pTable;

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



#line 7639 "./ESS_2001_bispectral.c"
#undef theta_3
#undef theta_2
#undef theta_1
#undef coating_thickness
#undef substrate_thickness
#undef transmit
#undef m
#undef length
#undef smallaxis
#undef yheight
#undef guide_start
#undef mirror_start
#undef focus_e
#undef focus_s
#undef pTable
#undef reflect
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'ArmForNeutronPropState_2' [11]. */
#define mccompcurname  ArmForNeutronPropState_2
#define mccompcurtype  Arm
#define mccompcurindex 11
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'guide_right' [12]. */
#define mccompcurname  guide_right
#define mccompcurtype  Mirror_Elliptic_Bispectral
#define mccompcurindex 12
#define reflect mccguide_right_reflect
#define pTable mccguide_right_pTable
#define focus_start_w mccguide_right_focus_start_w
#define focus_end_w mccguide_right_focus_end_w
#define focus_start_h mccguide_right_focus_start_h
#define focus_end_h mccguide_right_focus_end_h
#define mirror_start mccguide_right_mirror_start
#define m mccguide_right_m
#define smallaxis_w mccguide_right_smallaxis_w
#define smallaxis_h mccguide_right_smallaxis_h
#define length mccguide_right_length
#define transmit mccguide_right_transmit
#define substrate_thickness mccguide_right_substrate_thickness
#define coating_thickness mccguide_right_coating_thickness
#line 63 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../contrib/Mirror_Elliptic_Bispectral.comp"
t_Table pTable;

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

double zprime_w;


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

double b_w;
double f_w;
double asquared_w;
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


#line 7772 "./ESS_2001_bispectral.c"
#undef coating_thickness
#undef substrate_thickness
#undef transmit
#undef length
#undef smallaxis_h
#undef smallaxis_w
#undef m
#undef mirror_start
#undef focus_end_h
#undef focus_start_h
#undef focus_end_w
#undef focus_start_w
#undef pTable
#undef reflect
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'ArmForNeutronPropState_4' [13]. */
#define mccompcurname  ArmForNeutronPropState_4
#define mccompcurtype  Arm
#define mccompcurindex 13
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'guide_bottom' [14]. */
#define mccompcurname  guide_bottom
#define mccompcurtype  Mirror_Elliptic_Bispectral
#define mccompcurindex 14
#define reflect mccguide_bottom_reflect
#define pTable mccguide_bottom_pTable
#define focus_start_w mccguide_bottom_focus_start_w
#define focus_end_w mccguide_bottom_focus_end_w
#define focus_start_h mccguide_bottom_focus_start_h
#define focus_end_h mccguide_bottom_focus_end_h
#define mirror_start mccguide_bottom_mirror_start
#define m mccguide_bottom_m
#define smallaxis_w mccguide_bottom_smallaxis_w
#define smallaxis_h mccguide_bottom_smallaxis_h
#define length mccguide_bottom_length
#define transmit mccguide_bottom_transmit
#define substrate_thickness mccguide_bottom_substrate_thickness
#define coating_thickness mccguide_bottom_coating_thickness
#line 63 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../contrib/Mirror_Elliptic_Bispectral.comp"
t_Table pTable;

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

double zprime_w;


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

double b_w;
double f_w;
double asquared_w;
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


#line 7903 "./ESS_2001_bispectral.c"
#undef coating_thickness
#undef substrate_thickness
#undef transmit
#undef length
#undef smallaxis_h
#undef smallaxis_w
#undef m
#undef mirror_start
#undef focus_end_h
#undef focus_start_h
#undef focus_end_w
#undef focus_start_w
#undef pTable
#undef reflect
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'ArmForNeutronPropState_5' [15]. */
#define mccompcurname  ArmForNeutronPropState_5
#define mccompcurtype  Arm
#define mccompcurindex 15
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'cold_lambda_guidestart' [16]. */
#define mccompcurname  cold_lambda_guidestart
#define mccompcurtype  L_monitor
#define mccompcurindex 16
#define nL mcccold_lambda_guidestart_nL
#define L_N mcccold_lambda_guidestart_L_N
#define L_p mcccold_lambda_guidestart_L_p
#define L_p2 mcccold_lambda_guidestart_L_p2
#define filename mcccold_lambda_guidestart_filename
#define xmin mcccold_lambda_guidestart_xmin
#define xmax mcccold_lambda_guidestart_xmax
#define ymin mcccold_lambda_guidestart_ymin
#define ymax mcccold_lambda_guidestart_ymax
#define xwidth mcccold_lambda_guidestart_xwidth
#define yheight mcccold_lambda_guidestart_yheight
#define Lmin mcccold_lambda_guidestart_Lmin
#define Lmax mcccold_lambda_guidestart_Lmax
#define restore_neutron mcccold_lambda_guidestart_restore_neutron
#define nowritefile mcccold_lambda_guidestart_nowritefile
#line 57 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/L_monitor.comp"
double L_N[nL];
double L_p[nL], L_p2[nL];
#line 7952 "./ESS_2001_bispectral.c"
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

/* User declarations for component 'thermal_lambda_guidestart' [17]. */
#define mccompcurname  thermal_lambda_guidestart
#define mccompcurtype  L_monitor
#define mccompcurindex 17
#define nL mccthermal_lambda_guidestart_nL
#define L_N mccthermal_lambda_guidestart_L_N
#define L_p mccthermal_lambda_guidestart_L_p
#define L_p2 mccthermal_lambda_guidestart_L_p2
#define filename mccthermal_lambda_guidestart_filename
#define xmin mccthermal_lambda_guidestart_xmin
#define xmax mccthermal_lambda_guidestart_xmax
#define ymin mccthermal_lambda_guidestart_ymin
#define ymax mccthermal_lambda_guidestart_ymax
#define xwidth mccthermal_lambda_guidestart_xwidth
#define yheight mccthermal_lambda_guidestart_yheight
#define Lmin mccthermal_lambda_guidestart_Lmin
#define Lmax mccthermal_lambda_guidestart_Lmax
#define restore_neutron mccthermal_lambda_guidestart_restore_neutron
#define nowritefile mccthermal_lambda_guidestart_nowritefile
#line 57 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/L_monitor.comp"
double L_N[nL];
double L_p[nL], L_p2[nL];
#line 7994 "./ESS_2001_bispectral.c"
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

/* User declarations for component 'lambda_guidestart' [18]. */
#define mccompcurname  lambda_guidestart
#define mccompcurtype  L_monitor
#define mccompcurindex 18
#define nL mcclambda_guidestart_nL
#define L_N mcclambda_guidestart_L_N
#define L_p mcclambda_guidestart_L_p
#define L_p2 mcclambda_guidestart_L_p2
#define filename mcclambda_guidestart_filename
#define xmin mcclambda_guidestart_xmin
#define xmax mcclambda_guidestart_xmax
#define ymin mcclambda_guidestart_ymin
#define ymax mcclambda_guidestart_ymax
#define xwidth mcclambda_guidestart_xwidth
#define yheight mcclambda_guidestart_yheight
#define Lmin mcclambda_guidestart_Lmin
#define Lmax mcclambda_guidestart_Lmax
#define restore_neutron mcclambda_guidestart_restore_neutron
#define nowritefile mcclambda_guidestart_nowritefile
#line 57 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/L_monitor.comp"
double L_N[nL];
double L_p[nL], L_p2[nL];
#line 8036 "./ESS_2001_bispectral.c"
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

/* User declarations for component 'guide_top' [19]. */
#define mccompcurname  guide_top
#define mccompcurtype  Mirror_Elliptic_Bispectral
#define mccompcurindex 19
#define reflect mccguide_top_reflect
#define pTable mccguide_top_pTable
#define focus_start_w mccguide_top_focus_start_w
#define focus_end_w mccguide_top_focus_end_w
#define focus_start_h mccguide_top_focus_start_h
#define focus_end_h mccguide_top_focus_end_h
#define mirror_start mccguide_top_mirror_start
#define m mccguide_top_m
#define smallaxis_w mccguide_top_smallaxis_w
#define smallaxis_h mccguide_top_smallaxis_h
#define length mccguide_top_length
#define transmit mccguide_top_transmit
#define substrate_thickness mccguide_top_substrate_thickness
#define coating_thickness mccguide_top_coating_thickness
#line 63 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../contrib/Mirror_Elliptic_Bispectral.comp"
t_Table pTable;

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

double zprime_w;


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

double b_w;
double f_w;
double asquared_w;
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


#line 8160 "./ESS_2001_bispectral.c"
#undef coating_thickness
#undef substrate_thickness
#undef transmit
#undef length
#undef smallaxis_h
#undef smallaxis_w
#undef m
#undef mirror_start
#undef focus_end_h
#undef focus_start_h
#undef focus_end_w
#undef focus_start_w
#undef pTable
#undef reflect
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'ArmForNeutronPropState_6' [20]. */
#define mccompcurname  ArmForNeutronPropState_6
#define mccompcurtype  Arm
#define mccompcurindex 20
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'guide_Left' [21]. */
#define mccompcurname  guide_Left
#define mccompcurtype  Mirror_Elliptic_Bispectral
#define mccompcurindex 21
#define reflect mccguide_Left_reflect
#define pTable mccguide_Left_pTable
#define focus_start_w mccguide_Left_focus_start_w
#define focus_end_w mccguide_Left_focus_end_w
#define focus_start_h mccguide_Left_focus_start_h
#define focus_end_h mccguide_Left_focus_end_h
#define mirror_start mccguide_Left_mirror_start
#define m mccguide_Left_m
#define smallaxis_w mccguide_Left_smallaxis_w
#define smallaxis_h mccguide_Left_smallaxis_h
#define length mccguide_Left_length
#define transmit mccguide_Left_transmit
#define substrate_thickness mccguide_Left_substrate_thickness
#define coating_thickness mccguide_Left_coating_thickness
#line 63 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../contrib/Mirror_Elliptic_Bispectral.comp"
t_Table pTable;

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

double zprime_w;


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

double b_w;
double f_w;
double asquared_w;
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


#line 8291 "./ESS_2001_bispectral.c"
#undef coating_thickness
#undef substrate_thickness
#undef transmit
#undef length
#undef smallaxis_h
#undef smallaxis_w
#undef m
#undef mirror_start
#undef focus_end_h
#undef focus_start_h
#undef focus_end_w
#undef focus_start_w
#undef pTable
#undef reflect
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'ArmForNeutronPropState_7' [22]. */
#define mccompcurname  ArmForNeutronPropState_7
#define mccompcurtype  Arm
#define mccompcurindex 22
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'ArmMidTwo' [23]. */
#define mccompcurname  ArmMidTwo
#define mccompcurtype  Arm
#define mccompcurindex 23
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'ArmForNeutronPropState_8' [24]. */
#define mccompcurname  ArmForNeutronPropState_8
#define mccompcurtype  Arm
#define mccompcurindex 24
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'ArmMidThree' [25]. */
#define mccompcurname  ArmMidThree
#define mccompcurtype  Arm
#define mccompcurindex 25
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'ArmExit' [26]. */
#define mccompcurname  ArmExit
#define mccompcurtype  Arm
#define mccompcurindex 26
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

Coords mcposaOrigin, mcposrOrigin;
Rotation mcrotaOrigin, mcrotrOrigin;
Coords mcposaArmForGuideRight, mcposrArmForGuideRight;
Rotation mcrotaArmForGuideRight, mcrotrArmForGuideRight;
Coords mcposaArmForGuideBottom, mcposrArmForGuideBottom;
Rotation mcrotaArmForGuideBottom, mcrotrArmForGuideBottom;
Coords mcposaArmForGuideTop, mcposrArmForGuideTop;
Rotation mcrotaArmForGuideTop, mcrotrArmForGuideTop;
Coords mcposaArmForGuideLeft, mcposrArmForGuideLeft;
Rotation mcrotaArmForGuideLeft, mcrotrArmForGuideLeft;
Coords mcposacold_source, mcposrcold_source;
Rotation mcrotacold_source, mcrotrcold_source;
Coords mcposathermal_source, mcposrthermal_source;
Rotation mcrotathermal_source, mcrotrthermal_source;
Coords mcposaColdFocus, mcposrColdFocus;
Rotation mcrotaColdFocus, mcrotrColdFocus;
Coords mcposaArmMidOne, mcposrArmMidOne;
Rotation mcrotaArmMidOne, mcrotrArmMidOne;
Coords mcposamirror_full_center, mcposrmirror_full_center;
Rotation mcrotamirror_full_center, mcrotrmirror_full_center;
Coords mcposaArmForNeutronPropState_2, mcposrArmForNeutronPropState_2;
Rotation mcrotaArmForNeutronPropState_2, mcrotrArmForNeutronPropState_2;
Coords mcposaguide_right, mcposrguide_right;
Rotation mcrotaguide_right, mcrotrguide_right;
Coords mcposaArmForNeutronPropState_4, mcposrArmForNeutronPropState_4;
Rotation mcrotaArmForNeutronPropState_4, mcrotrArmForNeutronPropState_4;
Coords mcposaguide_bottom, mcposrguide_bottom;
Rotation mcrotaguide_bottom, mcrotrguide_bottom;
Coords mcposaArmForNeutronPropState_5, mcposrArmForNeutronPropState_5;
Rotation mcrotaArmForNeutronPropState_5, mcrotrArmForNeutronPropState_5;
Coords mcposacold_lambda_guidestart, mcposrcold_lambda_guidestart;
Rotation mcrotacold_lambda_guidestart, mcrotrcold_lambda_guidestart;
Coords mcposathermal_lambda_guidestart, mcposrthermal_lambda_guidestart;
Rotation mcrotathermal_lambda_guidestart, mcrotrthermal_lambda_guidestart;
Coords mcposalambda_guidestart, mcposrlambda_guidestart;
Rotation mcrotalambda_guidestart, mcrotrlambda_guidestart;
Coords mcposaguide_top, mcposrguide_top;
Rotation mcrotaguide_top, mcrotrguide_top;
Coords mcposaArmForNeutronPropState_6, mcposrArmForNeutronPropState_6;
Rotation mcrotaArmForNeutronPropState_6, mcrotrArmForNeutronPropState_6;
Coords mcposaguide_Left, mcposrguide_Left;
Rotation mcrotaguide_Left, mcrotrguide_Left;
Coords mcposaArmForNeutronPropState_7, mcposrArmForNeutronPropState_7;
Rotation mcrotaArmForNeutronPropState_7, mcrotrArmForNeutronPropState_7;
Coords mcposaArmMidTwo, mcposrArmMidTwo;
Rotation mcrotaArmMidTwo, mcrotrArmMidTwo;
Coords mcposaArmForNeutronPropState_8, mcposrArmForNeutronPropState_8;
Rotation mcrotaArmForNeutronPropState_8, mcrotrArmForNeutronPropState_8;
Coords mcposaArmMidThree, mcposrArmMidThree;
Rotation mcrotaArmMidThree, mcrotrArmMidThree;
Coords mcposaArmExit, mcposrArmExit;
Rotation mcrotaArmExit, mcrotrArmExit;

MCNUM mcnx, mcny, mcnz, mcnvx, mcnvy, mcnvz, mcnt, mcnsx, mcnsy, mcnsz, mcnp;

/* end declare */

void mcinit(void) {
#define mccompcurname  ESS_2001_bispectral
#define mccompcurtype  INSTRUMENT
#define mccompcurindex 0
#define mcposaESS_2001_bispectral coords_set(0,0,0)
#define man_lam mcipman_lam
#define thermal mcipthermal
#define cold_Lam_min mcipcold_Lam_min
#define cold_Lam_max mcipcold_Lam_max
#define thermal_Lam_min mcipthermal_Lam_min
#define thermal_Lam_max mcipthermal_Lam_max
#define thermal_hotspot_fac mcipthermal_hotspot_fac
#define cold_hotspot_fac mcipcold_hotspot_fac
#define use_full_guide_flag mcipuse_full_guide_flag
#define guide_start mcipguide_start
#define Length mcipLength
#define focus_start_w mcipfocus_start_w
#define focus_end_w mcipfocus_end_w
#define smallaxis_w mcipsmallaxis_w
#define focus_start_h mcipfocus_start_h
#define focus_end_h mcipfocus_end_h
#define smallaxis_h mcipsmallaxis_h
#define maxdiv mcipmaxdiv
#define coldperthermal mcipcoldperthermal
#define mirror_type mcipmirror_type
#define mirror_coating_type mcipmirror_coating_type
#define mirror_offset mcipmirror_offset
#define theta1 mciptheta1
#define theta2 mciptheta2
#define theta3 mciptheta3
#define m_mirror mcipm_mirror
#define h_mirror mciph_mirror
#define Pulse_width mcipPulse_width
#define frequency mcipfrequency
#define gravity mcipgravity
#define substrate_thickness mcipsubstrate_thickness
#define coating_thickness mcipcoating_thickness
#define m1 mcipm1
#define m2 mcipm2
#line 241 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
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
#line 8738 "./ESS_2001_bispectral.c"
#undef m2
#undef m1
#undef coating_thickness
#undef substrate_thickness
#undef gravity
#undef frequency
#undef Pulse_width
#undef h_mirror
#undef m_mirror
#undef theta3
#undef theta2
#undef theta1
#undef mirror_offset
#undef mirror_coating_type
#undef mirror_type
#undef coldperthermal
#undef maxdiv
#undef smallaxis_h
#undef focus_end_h
#undef focus_start_h
#undef smallaxis_w
#undef focus_end_w
#undef focus_start_w
#undef Length
#undef guide_start
#undef use_full_guide_flag
#undef cold_hotspot_fac
#undef thermal_hotspot_fac
#undef thermal_Lam_max
#undef thermal_Lam_min
#undef cold_Lam_max
#undef cold_Lam_min
#undef thermal
#undef man_lam
#undef mcposaESS_2001_bispectral
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
#line 39 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  if("NULL") strncpy(mccOrigin_profile, "NULL" ? "NULL" : "", 16384); else mccOrigin_profile[0]='\0';
#line 39 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mccOrigin_percent = 10;
#line 39 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mccOrigin_flag_save = 0;
#line 39 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mccOrigin_minutes = 0;
#line 8798 "./ESS_2001_bispectral.c"

  SIG_MESSAGE("Origin (Init:Place/Rotate)");
  rot_set_rotation(mcrotaOrigin,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 8805 "./ESS_2001_bispectral.c"
  rot_copy(mcrotrOrigin, mcrotaOrigin);
  mcposaOrigin = coords_set(
#line 538 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    0,
#line 538 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    0,
#line 538 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    0);
#line 8814 "./ESS_2001_bispectral.c"
  mctc1 = coords_neg(mcposaOrigin);
  mcposrOrigin = rot_apply(mcrotaOrigin, mctc1);
  mcDEBUG_COMPONENT("Origin", mcposaOrigin, mcrotaOrigin)
  mccomp_posa[1] = mcposaOrigin;
  mccomp_posr[1] = mcposrOrigin;
  mcNCounter[1]  = mcPCounter[1] = mcP2Counter[1] = 0;
  mcAbsorbProp[1]= 0;
    /* Component ArmForGuideRight. */
  /* Setting parameters for component ArmForGuideRight. */
  SIG_MESSAGE("ArmForGuideRight (Init:SetPar)");

  SIG_MESSAGE("ArmForGuideRight (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 553 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    (0)*DEG2RAD,
#line 553 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    (0)*DEG2RAD,
#line 553 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    (90)*DEG2RAD);
#line 8834 "./ESS_2001_bispectral.c"
  rot_mul(mctr1, mcrotaOrigin, mcrotaArmForGuideRight);
  rot_transpose(mcrotaOrigin, mctr1);
  rot_mul(mcrotaArmForGuideRight, mctr1, mcrotrArmForGuideRight);
  mctc1 = coords_set(
#line 552 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    0,
#line 552 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    0,
#line 552 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    0);
#line 8845 "./ESS_2001_bispectral.c"
  rot_transpose(mcrotaOrigin, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaArmForGuideRight = coords_add(mcposaOrigin, mctc2);
  mctc1 = coords_sub(mcposaOrigin, mcposaArmForGuideRight);
  mcposrArmForGuideRight = rot_apply(mcrotaArmForGuideRight, mctc1);
  mcDEBUG_COMPONENT("ArmForGuideRight", mcposaArmForGuideRight, mcrotaArmForGuideRight)
  mccomp_posa[2] = mcposaArmForGuideRight;
  mccomp_posr[2] = mcposrArmForGuideRight;
  mcNCounter[2]  = mcPCounter[2] = mcP2Counter[2] = 0;
  mcAbsorbProp[2]= 0;
    /* Component ArmForGuideBottom. */
  /* Setting parameters for component ArmForGuideBottom. */
  SIG_MESSAGE("ArmForGuideBottom (Init:SetPar)");

  SIG_MESSAGE("ArmForGuideBottom (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 555 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    (0)*DEG2RAD,
#line 555 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    (0)*DEG2RAD,
#line 555 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    (90)*DEG2RAD);
#line 8868 "./ESS_2001_bispectral.c"
  rot_mul(mctr1, mcrotaOrigin, mcrotaArmForGuideBottom);
  rot_transpose(mcrotaArmForGuideRight, mctr1);
  rot_mul(mcrotaArmForGuideBottom, mctr1, mcrotrArmForGuideBottom);
  mctc1 = coords_set(
#line 555 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    0,
#line 555 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    0,
#line 555 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    0);
#line 8879 "./ESS_2001_bispectral.c"
  rot_transpose(mcrotaOrigin, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaArmForGuideBottom = coords_add(mcposaOrigin, mctc2);
  mctc1 = coords_sub(mcposaArmForGuideRight, mcposaArmForGuideBottom);
  mcposrArmForGuideBottom = rot_apply(mcrotaArmForGuideBottom, mctc1);
  mcDEBUG_COMPONENT("ArmForGuideBottom", mcposaArmForGuideBottom, mcrotaArmForGuideBottom)
  mccomp_posa[3] = mcposaArmForGuideBottom;
  mccomp_posr[3] = mcposrArmForGuideBottom;
  mcNCounter[3]  = mcPCounter[3] = mcP2Counter[3] = 0;
  mcAbsorbProp[3]= 0;
    /* Component ArmForGuideTop. */
  /* Setting parameters for component ArmForGuideTop. */
  SIG_MESSAGE("ArmForGuideTop (Init:SetPar)");

  SIG_MESSAGE("ArmForGuideTop (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 557 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    (0)*DEG2RAD,
#line 557 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    (0)*DEG2RAD,
#line 557 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    (-90)*DEG2RAD);
#line 8902 "./ESS_2001_bispectral.c"
  rot_mul(mctr1, mcrotaOrigin, mcrotaArmForGuideTop);
  rot_transpose(mcrotaArmForGuideBottom, mctr1);
  rot_mul(mcrotaArmForGuideTop, mctr1, mcrotrArmForGuideTop);
  mctc1 = coords_set(
#line 557 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    0,
#line 557 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    0,
#line 557 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    0);
#line 8913 "./ESS_2001_bispectral.c"
  rot_transpose(mcrotaOrigin, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaArmForGuideTop = coords_add(mcposaOrigin, mctc2);
  mctc1 = coords_sub(mcposaArmForGuideBottom, mcposaArmForGuideTop);
  mcposrArmForGuideTop = rot_apply(mcrotaArmForGuideTop, mctc1);
  mcDEBUG_COMPONENT("ArmForGuideTop", mcposaArmForGuideTop, mcrotaArmForGuideTop)
  mccomp_posa[4] = mcposaArmForGuideTop;
  mccomp_posr[4] = mcposrArmForGuideTop;
  mcNCounter[4]  = mcPCounter[4] = mcP2Counter[4] = 0;
  mcAbsorbProp[4]= 0;
    /* Component ArmForGuideLeft. */
  /* Setting parameters for component ArmForGuideLeft. */
  SIG_MESSAGE("ArmForGuideLeft (Init:SetPar)");

  SIG_MESSAGE("ArmForGuideLeft (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 559 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    (0)*DEG2RAD,
#line 559 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    (0)*DEG2RAD,
#line 559 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    (180)*DEG2RAD);
#line 8936 "./ESS_2001_bispectral.c"
  rot_mul(mctr1, mcrotaOrigin, mcrotaArmForGuideLeft);
  rot_transpose(mcrotaArmForGuideTop, mctr1);
  rot_mul(mcrotaArmForGuideLeft, mctr1, mcrotrArmForGuideLeft);
  mctc1 = coords_set(
#line 559 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    0,
#line 559 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    0,
#line 559 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    0);
#line 8947 "./ESS_2001_bispectral.c"
  rot_transpose(mcrotaOrigin, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaArmForGuideLeft = coords_add(mcposaOrigin, mctc2);
  mctc1 = coords_sub(mcposaArmForGuideTop, mcposaArmForGuideLeft);
  mcposrArmForGuideLeft = rot_apply(mcrotaArmForGuideLeft, mctc1);
  mcDEBUG_COMPONENT("ArmForGuideLeft", mcposaArmForGuideLeft, mcrotaArmForGuideLeft)
  mccomp_posa[5] = mcposaArmForGuideLeft;
  mccomp_posr[5] = mcposrArmForGuideLeft;
  mcNCounter[5]  = mcPCounter[5] = mcP2Counter[5] = 0;
  mcAbsorbProp[5]= 0;
    /* Component cold_source. */
  /* Setting parameters for component cold_source. */
  SIG_MESSAGE("cold_source (Init:SetPar)");
#line 563 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mcccold_source_size = 0.12;
#line 563 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mcccold_source_l_low = mcipcold_Lam_min;
#line 563 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mcccold_source_l_high = mcipcold_Lam_max;
#line 88 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mcccold_source_dist = 0;
#line 564 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mcccold_source_xw = x_focus + 0.15;
#line 564 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mcccold_source_yh = y_focus + 0.15;
#line 564 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mcccold_source_freq = mcipfrequency;
#line 564 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mcccold_source_T = 50;
#line 564 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mcccold_source_tau = 287e-6;
#line 564 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mcccold_source_tau1 = 0;
#line 566 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mcccold_source_tau2 = 20e-6;
#line 566 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mcccold_source_d = mcipPulse_width;
#line 566 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mcccold_source_n = 20;
#line 566 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mcccold_source_n2 = 5;
#line 566 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mcccold_source_chi2 = 0.9;
#line 567 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mcccold_source_I0 = 6.9e11;
#line 567 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mcccold_source_I2 = 27.6e10;
#line 567 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mcccold_source_branch1 = 1.0;
#line 567 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mcccold_source_branch2 = 0.5;
#line 568 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mcccold_source_branch_tail = 0.1;
#line 91 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mcccold_source_twopulses = 0;
#line 568 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mcccold_source_target_index = 2;
#line 9005 "./ESS_2001_bispectral.c"

  SIG_MESSAGE("cold_source (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 9012 "./ESS_2001_bispectral.c"
  rot_mul(mctr1, mcrotaOrigin, mcrotacold_source);
  rot_transpose(mcrotaArmForGuideLeft, mctr1);
  rot_mul(mcrotacold_source, mctr1, mcrotrcold_source);
  mctc1 = coords_set(
#line 569 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    cold_moderator_x,
#line 569 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    0,
#line 569 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    0);
#line 9023 "./ESS_2001_bispectral.c"
  rot_transpose(mcrotaOrigin, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposacold_source = coords_add(mcposaOrigin, mctc2);
  mctc1 = coords_sub(mcposaArmForGuideLeft, mcposacold_source);
  mcposrcold_source = rot_apply(mcrotacold_source, mctc1);
  mcDEBUG_COMPONENT("cold_source", mcposacold_source, mcrotacold_source)
  mccomp_posa[6] = mcposacold_source;
  mccomp_posr[6] = mcposrcold_source;
  mcNCounter[6]  = mcPCounter[6] = mcP2Counter[6] = 0;
  mcAbsorbProp[6]= 0;
    /* Component thermal_source. */
  /* Setting parameters for component thermal_source. */
  SIG_MESSAGE("thermal_source (Init:SetPar)");
#line 590 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mccthermal_source_size = 0.12;
#line 590 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mccthermal_source_l_low = mcipthermal_Lam_min;
#line 590 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mccthermal_source_l_high = mcipthermal_Lam_max;
#line 590 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mccthermal_source_dist = extraction_start_pos;
#line 591 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mccthermal_source_xw = 2 * w_guide_leftstart + 0.1;
#line 591 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mccthermal_source_yh = 2 * h_guide_leftstart + 0.1;
#line 591 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mccthermal_source_freq = mcipfrequency;
#line 591 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mccthermal_source_T = 325;
#line 591 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mccthermal_source_tau = 80e-6;
#line 591 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mccthermal_source_tau1 = 400e-6;
#line 592 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mccthermal_source_tau2 = 12e-6;
#line 592 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mccthermal_source_d = mcipPulse_width;
#line 592 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mccthermal_source_n = 20;
#line 592 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mccthermal_source_n2 = 5;
#line 592 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mccthermal_source_chi2 = 2.5;
#line 593 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mccthermal_source_I0 = 13.5e11;
#line 593 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mccthermal_source_I2 = 27.6e10;
#line 593 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mccthermal_source_branch1 = 0.5;
#line 593 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mccthermal_source_branch2 = 0.5;
#line 594 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mccthermal_source_branch_tail = 0.1;
#line 91 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mccthermal_source_twopulses = 0;
#line 91 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mccthermal_source_target_index = 0;
#line 9081 "./ESS_2001_bispectral.c"

  SIG_MESSAGE("thermal_source (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 9088 "./ESS_2001_bispectral.c"
  rot_mul(mctr1, mcrotaOrigin, mcrotathermal_source);
  rot_transpose(mcrotacold_source, mctr1);
  rot_mul(mcrotathermal_source, mctr1, mcrotrthermal_source);
  mctc1 = coords_set(
#line 595 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    0,
#line 595 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    0,
#line 595 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    0);
#line 9099 "./ESS_2001_bispectral.c"
  rot_transpose(mcrotaOrigin, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposathermal_source = coords_add(mcposaOrigin, mctc2);
  mctc1 = coords_sub(mcposacold_source, mcposathermal_source);
  mcposrthermal_source = rot_apply(mcrotathermal_source, mctc1);
  mcDEBUG_COMPONENT("thermal_source", mcposathermal_source, mcrotathermal_source)
  mccomp_posa[7] = mcposathermal_source;
  mccomp_posr[7] = mcposrthermal_source;
  mcNCounter[7]  = mcPCounter[7] = mcP2Counter[7] = 0;
  mcAbsorbProp[7]= 0;
    /* Component ColdFocus. */
  /* Setting parameters for component ColdFocus. */
  SIG_MESSAGE("ColdFocus (Init:SetPar)");

  SIG_MESSAGE("ColdFocus (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 9119 "./ESS_2001_bispectral.c"
  rot_mul(mctr1, mcrotaOrigin, mcrotaColdFocus);
  rot_transpose(mcrotathermal_source, mctr1);
  rot_mul(mcrotaColdFocus, mctr1, mcrotrColdFocus);
  mctc1 = coords_set(
#line 614 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    x_mid,
#line 614 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    0,
#line 614 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    extraction_start_pos);
#line 9130 "./ESS_2001_bispectral.c"
  rot_transpose(mcrotaOrigin, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaColdFocus = coords_add(mcposaOrigin, mctc2);
  mctc1 = coords_sub(mcposathermal_source, mcposaColdFocus);
  mcposrColdFocus = rot_apply(mcrotaColdFocus, mctc1);
  mcDEBUG_COMPONENT("ColdFocus", mcposaColdFocus, mcrotaColdFocus)
  mccomp_posa[8] = mcposaColdFocus;
  mccomp_posr[8] = mcposrColdFocus;
  mcNCounter[8]  = mcPCounter[8] = mcP2Counter[8] = 0;
  mcAbsorbProp[8]= 0;
    /* Component ArmMidOne. */
  /* Setting parameters for component ArmMidOne. */
  SIG_MESSAGE("ArmMidOne (Init:SetPar)");

  SIG_MESSAGE("ArmMidOne (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 9150 "./ESS_2001_bispectral.c"
  rot_mul(mctr1, mcrotaOrigin, mcrotaArmMidOne);
  rot_transpose(mcrotaColdFocus, mctr1);
  rot_mul(mcrotaArmMidOne, mctr1, mcrotrArmMidOne);
  mctc1 = coords_set(
#line 625 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    0,
#line 625 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    0,
#line 625 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    ArmMidOnePos);
#line 9161 "./ESS_2001_bispectral.c"
  rot_transpose(mcrotaOrigin, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaArmMidOne = coords_add(mcposaOrigin, mctc2);
  mctc1 = coords_sub(mcposaColdFocus, mcposaArmMidOne);
  mcposrArmMidOne = rot_apply(mcrotaArmMidOne, mctc1);
  mcDEBUG_COMPONENT("ArmMidOne", mcposaArmMidOne, mcrotaArmMidOne)
  mccomp_posa[9] = mcposaArmMidOne;
  mccomp_posr[9] = mcposrArmMidOne;
  mcNCounter[9]  = mcPCounter[9] = mcP2Counter[9] = 0;
  mcAbsorbProp[9]= 0;
    /* Component mirror_full_center. */
  /* Setting parameters for component mirror_full_center. */
  SIG_MESSAGE("mirror_full_center (Init:SetPar)");
#line 657 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mccmirror_full_center_focus_s = mcipfocus_start_h;
#line 657 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mccmirror_full_center_focus_e = mcipfocus_end_h;
#line 657 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mccmirror_full_center_mirror_start = extraction_start_pos;
#line 657 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mccmirror_full_center_guide_start = mcipguide_start;
#line 657 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mccmirror_full_center_yheight = mciph_mirror;
#line 657 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mccmirror_full_center_smallaxis = mcipsmallaxis_h;
#line 657 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mccmirror_full_center_length = 2.5;
#line 658 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mccmirror_full_center_m = mcipm_mirror;
#line 658 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mccmirror_full_center_transmit = 1;
#line 656 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mccmirror_full_center_substrate_thickness = mcipsubstrate_thickness;
#line 656 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mccmirror_full_center_coating_thickness = mcipcoating_thickness;
#line 658 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mccmirror_full_center_theta_1 = mciptheta1;
#line 658 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mccmirror_full_center_theta_2 = mciptheta1;
#line 658 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mccmirror_full_center_theta_3 = mciptheta1;
#line 9203 "./ESS_2001_bispectral.c"

  SIG_MESSAGE("mirror_full_center (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 661 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    (0)*DEG2RAD,
#line 661 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    (-90)*DEG2RAD,
#line 661 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    (0)*DEG2RAD);
#line 9213 "./ESS_2001_bispectral.c"
  rot_mul(mctr1, mcrotaOrigin, mcrotamirror_full_center);
  rot_transpose(mcrotaArmMidOne, mctr1);
  rot_mul(mcrotamirror_full_center, mctr1, mcrotrmirror_full_center);
  mctc1 = coords_set(
#line 660 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    mcipmirror_offset,
#line 660 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    0,
#line 660 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    extraction_start_pos + mirror_full_length * 0.5);
#line 9224 "./ESS_2001_bispectral.c"
  rot_transpose(mcrotaOrigin, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposamirror_full_center = coords_add(mcposaOrigin, mctc2);
  mctc1 = coords_sub(mcposaArmMidOne, mcposamirror_full_center);
  mcposrmirror_full_center = rot_apply(mcrotamirror_full_center, mctc1);
  mcDEBUG_COMPONENT("mirror_full_center", mcposamirror_full_center, mcrotamirror_full_center)
  mccomp_posa[10] = mcposamirror_full_center;
  mccomp_posr[10] = mcposrmirror_full_center;
  mcNCounter[10]  = mcPCounter[10] = mcP2Counter[10] = 0;
  mcAbsorbProp[10]= 0;
    /* Component ArmForNeutronPropState_2. */
  /* Setting parameters for component ArmForNeutronPropState_2. */
  SIG_MESSAGE("ArmForNeutronPropState_2 (Init:SetPar)");

  SIG_MESSAGE("ArmForNeutronPropState_2 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 9244 "./ESS_2001_bispectral.c"
  rot_mul(mctr1, mcrotaOrigin, mcrotaArmForNeutronPropState_2);
  rot_transpose(mcrotamirror_full_center, mctr1);
  rot_mul(mcrotaArmForNeutronPropState_2, mctr1, mcrotrArmForNeutronPropState_2);
  mctc1 = coords_set(
#line 668 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    0,
#line 668 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    0,
#line 668 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    ArmMidOnePos);
#line 9255 "./ESS_2001_bispectral.c"
  rot_transpose(mcrotaOrigin, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaArmForNeutronPropState_2 = coords_add(mcposaOrigin, mctc2);
  mctc1 = coords_sub(mcposamirror_full_center, mcposaArmForNeutronPropState_2);
  mcposrArmForNeutronPropState_2 = rot_apply(mcrotaArmForNeutronPropState_2, mctc1);
  mcDEBUG_COMPONENT("ArmForNeutronPropState_2", mcposaArmForNeutronPropState_2, mcrotaArmForNeutronPropState_2)
  mccomp_posa[11] = mcposaArmForNeutronPropState_2;
  mccomp_posr[11] = mcposrArmForNeutronPropState_2;
  mcNCounter[11]  = mcPCounter[11] = mcP2Counter[11] = 0;
  mcAbsorbProp[11]= 0;
    /* Component guide_right. */
  /* Setting parameters for component guide_right. */
  SIG_MESSAGE("guide_right (Init:SetPar)");
#line 703 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mccguide_right_focus_start_w = mcipfocus_start_w;
#line 703 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mccguide_right_focus_end_w = mcipfocus_end_w;
#line 702 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mccguide_right_focus_start_h = mcipfocus_start_h;
#line 702 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mccguide_right_focus_end_h = mcipfocus_end_h;
#line 702 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mccguide_right_mirror_start = mcipguide_start;
#line 704 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mccguide_right_m = m [ 0 ];
#line 703 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mccguide_right_smallaxis_w = mcipsmallaxis_w;
#line 702 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mccguide_right_smallaxis_h = mcipsmallaxis_h;
#line 702 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mccguide_right_length = extraction_start_pos + mirror_full_length - mcipguide_start;
#line 52 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mccguide_right_transmit = 0;
#line 701 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mccguide_right_substrate_thickness = 0;
#line 701 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mccguide_right_coating_thickness = 0;
#line 9293 "./ESS_2001_bispectral.c"

  SIG_MESSAGE("guide_right (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 707 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    (0)*DEG2RAD,
#line 707 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    (-90)*DEG2RAD,
#line 707 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    (0)*DEG2RAD);
#line 9303 "./ESS_2001_bispectral.c"
  rot_mul(mctr1, mcrotaOrigin, mcrotaguide_right);
  rot_transpose(mcrotaArmForNeutronPropState_2, mctr1);
  rot_mul(mcrotaguide_right, mctr1, mcrotrguide_right);
  mctc1 = coords_set(
#line 706 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    0,
#line 706 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    0,
#line 706 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    mcipguide_start);
#line 9314 "./ESS_2001_bispectral.c"
  rot_transpose(mcrotaOrigin, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaguide_right = coords_add(mcposaOrigin, mctc2);
  mctc1 = coords_sub(mcposaArmForNeutronPropState_2, mcposaguide_right);
  mcposrguide_right = rot_apply(mcrotaguide_right, mctc1);
  mcDEBUG_COMPONENT("guide_right", mcposaguide_right, mcrotaguide_right)
  mccomp_posa[12] = mcposaguide_right;
  mccomp_posr[12] = mcposrguide_right;
  mcNCounter[12]  = mcPCounter[12] = mcP2Counter[12] = 0;
  mcAbsorbProp[12]= 0;
    /* Component ArmForNeutronPropState_4. */
  /* Setting parameters for component ArmForNeutronPropState_4. */
  SIG_MESSAGE("ArmForNeutronPropState_4 (Init:SetPar)");

  SIG_MESSAGE("ArmForNeutronPropState_4 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 9334 "./ESS_2001_bispectral.c"
  rot_mul(mctr1, mcrotaOrigin, mcrotaArmForNeutronPropState_4);
  rot_transpose(mcrotaguide_right, mctr1);
  rot_mul(mcrotaArmForNeutronPropState_4, mctr1, mcrotrArmForNeutronPropState_4);
  mctc1 = coords_set(
#line 714 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    0,
#line 714 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    0,
#line 714 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    ArmMidOnePos);
#line 9345 "./ESS_2001_bispectral.c"
  rot_transpose(mcrotaOrigin, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaArmForNeutronPropState_4 = coords_add(mcposaOrigin, mctc2);
  mctc1 = coords_sub(mcposaguide_right, mcposaArmForNeutronPropState_4);
  mcposrArmForNeutronPropState_4 = rot_apply(mcrotaArmForNeutronPropState_4, mctc1);
  mcDEBUG_COMPONENT("ArmForNeutronPropState_4", mcposaArmForNeutronPropState_4, mcrotaArmForNeutronPropState_4)
  mccomp_posa[13] = mcposaArmForNeutronPropState_4;
  mccomp_posr[13] = mcposrArmForNeutronPropState_4;
  mcNCounter[13]  = mcPCounter[13] = mcP2Counter[13] = 0;
  mcAbsorbProp[13]= 0;
    /* Component guide_bottom. */
  /* Setting parameters for component guide_bottom. */
  SIG_MESSAGE("guide_bottom (Init:SetPar)");
#line 753 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mccguide_bottom_focus_start_w = mcipfocus_start_h;
#line 753 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mccguide_bottom_focus_end_w = mcipfocus_end_h;
#line 752 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mccguide_bottom_focus_start_h = mcipfocus_start_w;
#line 752 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mccguide_bottom_focus_end_h = mcipfocus_end_w;
#line 752 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mccguide_bottom_mirror_start = mcipguide_start;
#line 754 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mccguide_bottom_m = m [ 0 ];
#line 753 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mccguide_bottom_smallaxis_w = mcipsmallaxis_h;
#line 752 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mccguide_bottom_smallaxis_h = mcipsmallaxis_w;
#line 752 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mccguide_bottom_length = extraction_start_pos + mirror_full_length - mcipguide_start;
#line 52 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mccguide_bottom_transmit = 0;
#line 751 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mccguide_bottom_substrate_thickness = 0;
#line 751 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mccguide_bottom_coating_thickness = 0;
#line 9383 "./ESS_2001_bispectral.c"

  SIG_MESSAGE("guide_bottom (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 757 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    (0)*DEG2RAD,
#line 757 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    (-90)*DEG2RAD,
#line 757 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    (0)*DEG2RAD);
#line 9393 "./ESS_2001_bispectral.c"
  rot_mul(mctr1, mcrotaArmForGuideBottom, mcrotaguide_bottom);
  rot_transpose(mcrotaArmForNeutronPropState_4, mctr1);
  rot_mul(mcrotaguide_bottom, mctr1, mcrotrguide_bottom);
  mctc1 = coords_set(
#line 756 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    0,
#line 756 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    0,
#line 756 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    mcipguide_start);
#line 9404 "./ESS_2001_bispectral.c"
  rot_transpose(mcrotaArmForGuideBottom, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaguide_bottom = coords_add(mcposaArmForGuideBottom, mctc2);
  mctc1 = coords_sub(mcposaArmForNeutronPropState_4, mcposaguide_bottom);
  mcposrguide_bottom = rot_apply(mcrotaguide_bottom, mctc1);
  mcDEBUG_COMPONENT("guide_bottom", mcposaguide_bottom, mcrotaguide_bottom)
  mccomp_posa[14] = mcposaguide_bottom;
  mccomp_posr[14] = mcposrguide_bottom;
  mcNCounter[14]  = mcPCounter[14] = mcP2Counter[14] = 0;
  mcAbsorbProp[14]= 0;
    /* Component ArmForNeutronPropState_5. */
  /* Setting parameters for component ArmForNeutronPropState_5. */
  SIG_MESSAGE("ArmForNeutronPropState_5 (Init:SetPar)");

  SIG_MESSAGE("ArmForNeutronPropState_5 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 9424 "./ESS_2001_bispectral.c"
  rot_mul(mctr1, mcrotaOrigin, mcrotaArmForNeutronPropState_5);
  rot_transpose(mcrotaguide_bottom, mctr1);
  rot_mul(mcrotaArmForNeutronPropState_5, mctr1, mcrotrArmForNeutronPropState_5);
  mctc1 = coords_set(
#line 764 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    0,
#line 764 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    0,
#line 764 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    ArmMidOnePos);
#line 9435 "./ESS_2001_bispectral.c"
  rot_transpose(mcrotaOrigin, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaArmForNeutronPropState_5 = coords_add(mcposaOrigin, mctc2);
  mctc1 = coords_sub(mcposaguide_bottom, mcposaArmForNeutronPropState_5);
  mcposrArmForNeutronPropState_5 = rot_apply(mcrotaArmForNeutronPropState_5, mctc1);
  mcDEBUG_COMPONENT("ArmForNeutronPropState_5", mcposaArmForNeutronPropState_5, mcrotaArmForNeutronPropState_5)
  mccomp_posa[15] = mcposaArmForNeutronPropState_5;
  mccomp_posr[15] = mcposrArmForNeutronPropState_5;
  mcNCounter[15]  = mcPCounter[15] = mcP2Counter[15] = 0;
  mcAbsorbProp[15]= 0;
    /* Component cold_lambda_guidestart. */
  /* Setting parameters for component cold_lambda_guidestart. */
  SIG_MESSAGE("cold_lambda_guidestart (Init:SetPar)");
#line 800 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  if("cold_lambda_guidestart") strncpy(mcccold_lambda_guidestart_filename, "cold_lambda_guidestart" ? "cold_lambda_guidestart" : "", 16384); else mcccold_lambda_guidestart_filename[0]='\0';
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mcccold_lambda_guidestart_xmin = -0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mcccold_lambda_guidestart_xmax = 0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mcccold_lambda_guidestart_ymin = -0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mcccold_lambda_guidestart_ymax = 0.05;
#line 801 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mcccold_lambda_guidestart_xwidth = w [ 0 ];
#line 801 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mcccold_lambda_guidestart_yheight = h [ 0 ];
#line 801 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mcccold_lambda_guidestart_Lmin = 0.01;
#line 801 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mcccold_lambda_guidestart_Lmax = 20;
#line 800 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mcccold_lambda_guidestart_restore_neutron = 1;
#line 51 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mcccold_lambda_guidestart_nowritefile = 0;
#line 9471 "./ESS_2001_bispectral.c"

  SIG_MESSAGE("cold_lambda_guidestart (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 9478 "./ESS_2001_bispectral.c"
  rot_mul(mctr1, mcrotaOrigin, mcrotacold_lambda_guidestart);
  rot_transpose(mcrotaArmForNeutronPropState_5, mctr1);
  rot_mul(mcrotacold_lambda_guidestart, mctr1, mcrotrcold_lambda_guidestart);
  mctc1 = coords_set(
#line 802 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    0,
#line 802 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    0,
#line 802 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    guide_z_pos [ 0 ]);
#line 9489 "./ESS_2001_bispectral.c"
  rot_transpose(mcrotaOrigin, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposacold_lambda_guidestart = coords_add(mcposaOrigin, mctc2);
  mctc1 = coords_sub(mcposaArmForNeutronPropState_5, mcposacold_lambda_guidestart);
  mcposrcold_lambda_guidestart = rot_apply(mcrotacold_lambda_guidestart, mctc1);
  mcDEBUG_COMPONENT("cold_lambda_guidestart", mcposacold_lambda_guidestart, mcrotacold_lambda_guidestart)
  mccomp_posa[16] = mcposacold_lambda_guidestart;
  mccomp_posr[16] = mcposrcold_lambda_guidestart;
  mcNCounter[16]  = mcPCounter[16] = mcP2Counter[16] = 0;
  mcAbsorbProp[16]= 0;
    /* Component thermal_lambda_guidestart. */
  /* Setting parameters for component thermal_lambda_guidestart. */
  SIG_MESSAGE("thermal_lambda_guidestart (Init:SetPar)");
#line 805 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  if("thermal_lambda_guidestart") strncpy(mccthermal_lambda_guidestart_filename, "thermal_lambda_guidestart" ? "thermal_lambda_guidestart" : "", 16384); else mccthermal_lambda_guidestart_filename[0]='\0';
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mccthermal_lambda_guidestart_xmin = -0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mccthermal_lambda_guidestart_xmax = 0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mccthermal_lambda_guidestart_ymin = -0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mccthermal_lambda_guidestart_ymax = 0.05;
#line 806 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mccthermal_lambda_guidestart_xwidth = w [ 0 ];
#line 806 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mccthermal_lambda_guidestart_yheight = h [ 0 ];
#line 806 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mccthermal_lambda_guidestart_Lmin = 0.01;
#line 806 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mccthermal_lambda_guidestart_Lmax = 20;
#line 805 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mccthermal_lambda_guidestart_restore_neutron = 1;
#line 51 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mccthermal_lambda_guidestart_nowritefile = 0;
#line 9525 "./ESS_2001_bispectral.c"

  SIG_MESSAGE("thermal_lambda_guidestart (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 9532 "./ESS_2001_bispectral.c"
  rot_mul(mctr1, mcrotaOrigin, mcrotathermal_lambda_guidestart);
  rot_transpose(mcrotacold_lambda_guidestart, mctr1);
  rot_mul(mcrotathermal_lambda_guidestart, mctr1, mcrotrthermal_lambda_guidestart);
  mctc1 = coords_set(
#line 807 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    0,
#line 807 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    0,
#line 807 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    guide_z_pos [ 0 ]);
#line 9543 "./ESS_2001_bispectral.c"
  rot_transpose(mcrotaOrigin, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposathermal_lambda_guidestart = coords_add(mcposaOrigin, mctc2);
  mctc1 = coords_sub(mcposacold_lambda_guidestart, mcposathermal_lambda_guidestart);
  mcposrthermal_lambda_guidestart = rot_apply(mcrotathermal_lambda_guidestart, mctc1);
  mcDEBUG_COMPONENT("thermal_lambda_guidestart", mcposathermal_lambda_guidestart, mcrotathermal_lambda_guidestart)
  mccomp_posa[17] = mcposathermal_lambda_guidestart;
  mccomp_posr[17] = mcposrthermal_lambda_guidestart;
  mcNCounter[17]  = mcPCounter[17] = mcP2Counter[17] = 0;
  mcAbsorbProp[17]= 0;
    /* Component lambda_guidestart. */
  /* Setting parameters for component lambda_guidestart. */
  SIG_MESSAGE("lambda_guidestart (Init:SetPar)");
#line 810 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  if("lambda_guidestart") strncpy(mcclambda_guidestart_filename, "lambda_guidestart" ? "lambda_guidestart" : "", 16384); else mcclambda_guidestart_filename[0]='\0';
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mcclambda_guidestart_xmin = -0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mcclambda_guidestart_xmax = 0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mcclambda_guidestart_ymin = -0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mcclambda_guidestart_ymax = 0.05;
#line 811 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mcclambda_guidestart_xwidth = w [ 0 ];
#line 811 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mcclambda_guidestart_yheight = h [ 0 ];
#line 811 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mcclambda_guidestart_Lmin = 0.01;
#line 811 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mcclambda_guidestart_Lmax = 20;
#line 810 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mcclambda_guidestart_restore_neutron = 1;
#line 51 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mcclambda_guidestart_nowritefile = 0;
#line 9579 "./ESS_2001_bispectral.c"

  SIG_MESSAGE("lambda_guidestart (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 9586 "./ESS_2001_bispectral.c"
  rot_mul(mctr1, mcrotaOrigin, mcrotalambda_guidestart);
  rot_transpose(mcrotathermal_lambda_guidestart, mctr1);
  rot_mul(mcrotalambda_guidestart, mctr1, mcrotrlambda_guidestart);
  mctc1 = coords_set(
#line 812 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    0,
#line 812 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    0,
#line 812 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    guide_z_pos [ 0 ]);
#line 9597 "./ESS_2001_bispectral.c"
  rot_transpose(mcrotaOrigin, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposalambda_guidestart = coords_add(mcposaOrigin, mctc2);
  mctc1 = coords_sub(mcposathermal_lambda_guidestart, mcposalambda_guidestart);
  mcposrlambda_guidestart = rot_apply(mcrotalambda_guidestart, mctc1);
  mcDEBUG_COMPONENT("lambda_guidestart", mcposalambda_guidestart, mcrotalambda_guidestart)
  mccomp_posa[18] = mcposalambda_guidestart;
  mccomp_posr[18] = mcposrlambda_guidestart;
  mcNCounter[18]  = mcPCounter[18] = mcP2Counter[18] = 0;
  mcAbsorbProp[18]= 0;
    /* Component guide_top. */
  /* Setting parameters for component guide_top. */
  SIG_MESSAGE("guide_top (Init:SetPar)");
#line 818 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mccguide_top_focus_start_w = mcipfocus_start_h;
#line 818 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mccguide_top_focus_end_w = mcipfocus_end_h;
#line 817 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mccguide_top_focus_start_h = mcipfocus_start_w;
#line 817 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mccguide_top_focus_end_h = mcipfocus_end_w;
#line 817 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mccguide_top_mirror_start = mcipguide_start;
#line 819 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mccguide_top_m = m [ 0 ];
#line 818 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mccguide_top_smallaxis_w = mcipsmallaxis_h;
#line 817 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mccguide_top_smallaxis_h = mcipsmallaxis_w;
#line 817 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mccguide_top_length = extraction_start_pos + mirror_full_length - mcipguide_start;
#line 52 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mccguide_top_transmit = 0;
#line 816 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mccguide_top_substrate_thickness = 0;
#line 816 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mccguide_top_coating_thickness = 0;
#line 9635 "./ESS_2001_bispectral.c"

  SIG_MESSAGE("guide_top (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 822 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    (0)*DEG2RAD,
#line 822 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    (-90)*DEG2RAD,
#line 822 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    (0)*DEG2RAD);
#line 9645 "./ESS_2001_bispectral.c"
  rot_mul(mctr1, mcrotaArmForGuideTop, mcrotaguide_top);
  rot_transpose(mcrotalambda_guidestart, mctr1);
  rot_mul(mcrotaguide_top, mctr1, mcrotrguide_top);
  mctc1 = coords_set(
#line 821 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    0,
#line 821 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    0,
#line 821 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    mcipguide_start);
#line 9656 "./ESS_2001_bispectral.c"
  rot_transpose(mcrotaArmForGuideTop, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaguide_top = coords_add(mcposaArmForGuideTop, mctc2);
  mctc1 = coords_sub(mcposalambda_guidestart, mcposaguide_top);
  mcposrguide_top = rot_apply(mcrotaguide_top, mctc1);
  mcDEBUG_COMPONENT("guide_top", mcposaguide_top, mcrotaguide_top)
  mccomp_posa[19] = mcposaguide_top;
  mccomp_posr[19] = mcposrguide_top;
  mcNCounter[19]  = mcPCounter[19] = mcP2Counter[19] = 0;
  mcAbsorbProp[19]= 0;
    /* Component ArmForNeutronPropState_6. */
  /* Setting parameters for component ArmForNeutronPropState_6. */
  SIG_MESSAGE("ArmForNeutronPropState_6 (Init:SetPar)");

  SIG_MESSAGE("ArmForNeutronPropState_6 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 9676 "./ESS_2001_bispectral.c"
  rot_mul(mctr1, mcrotaOrigin, mcrotaArmForNeutronPropState_6);
  rot_transpose(mcrotaguide_top, mctr1);
  rot_mul(mcrotaArmForNeutronPropState_6, mctr1, mcrotrArmForNeutronPropState_6);
  mctc1 = coords_set(
#line 829 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    0,
#line 829 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    0,
#line 829 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    ArmMidOnePos);
#line 9687 "./ESS_2001_bispectral.c"
  rot_transpose(mcrotaOrigin, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaArmForNeutronPropState_6 = coords_add(mcposaOrigin, mctc2);
  mctc1 = coords_sub(mcposaguide_top, mcposaArmForNeutronPropState_6);
  mcposrArmForNeutronPropState_6 = rot_apply(mcrotaArmForNeutronPropState_6, mctc1);
  mcDEBUG_COMPONENT("ArmForNeutronPropState_6", mcposaArmForNeutronPropState_6, mcrotaArmForNeutronPropState_6)
  mccomp_posa[20] = mcposaArmForNeutronPropState_6;
  mccomp_posr[20] = mcposrArmForNeutronPropState_6;
  mcNCounter[20]  = mcPCounter[20] = mcP2Counter[20] = 0;
  mcAbsorbProp[20]= 0;
    /* Component guide_Left. */
  /* Setting parameters for component guide_Left. */
  SIG_MESSAGE("guide_Left (Init:SetPar)");
#line 868 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mccguide_Left_focus_start_w = mcipfocus_start_w;
#line 868 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mccguide_Left_focus_end_w = mcipfocus_end_w;
#line 867 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mccguide_Left_focus_start_h = mcipfocus_start_h;
#line 867 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mccguide_Left_focus_end_h = mcipfocus_end_h;
#line 867 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mccguide_Left_mirror_start = guide_left_start;
#line 869 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mccguide_Left_m = m [ 0 ];
#line 868 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mccguide_Left_smallaxis_w = mcipsmallaxis_w;
#line 867 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mccguide_Left_smallaxis_h = mcipsmallaxis_h;
#line 867 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mccguide_Left_length = extraction_start_pos + mirror_full_length - guide_left_start;
#line 52 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mccguide_Left_transmit = 0;
#line 866 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mccguide_Left_substrate_thickness = 0;
#line 866 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
  mccguide_Left_coating_thickness = 0;
#line 9725 "./ESS_2001_bispectral.c"

  SIG_MESSAGE("guide_Left (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 872 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    (0)*DEG2RAD,
#line 872 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    (-90)*DEG2RAD,
#line 872 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    (0)*DEG2RAD);
#line 9735 "./ESS_2001_bispectral.c"
  rot_mul(mctr1, mcrotaArmForGuideLeft, mcrotaguide_Left);
  rot_transpose(mcrotaArmForNeutronPropState_6, mctr1);
  rot_mul(mcrotaguide_Left, mctr1, mcrotrguide_Left);
  mctc1 = coords_set(
#line 871 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    0,
#line 871 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    0,
#line 871 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    guide_left_start);
#line 9746 "./ESS_2001_bispectral.c"
  rot_transpose(mcrotaArmForGuideLeft, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaguide_Left = coords_add(mcposaArmForGuideLeft, mctc2);
  mctc1 = coords_sub(mcposaArmForNeutronPropState_6, mcposaguide_Left);
  mcposrguide_Left = rot_apply(mcrotaguide_Left, mctc1);
  mcDEBUG_COMPONENT("guide_Left", mcposaguide_Left, mcrotaguide_Left)
  mccomp_posa[21] = mcposaguide_Left;
  mccomp_posr[21] = mcposrguide_Left;
  mcNCounter[21]  = mcPCounter[21] = mcP2Counter[21] = 0;
  mcAbsorbProp[21]= 0;
    /* Component ArmForNeutronPropState_7. */
  /* Setting parameters for component ArmForNeutronPropState_7. */
  SIG_MESSAGE("ArmForNeutronPropState_7 (Init:SetPar)");

  SIG_MESSAGE("ArmForNeutronPropState_7 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 9766 "./ESS_2001_bispectral.c"
  rot_mul(mctr1, mcrotaOrigin, mcrotaArmForNeutronPropState_7);
  rot_transpose(mcrotaguide_Left, mctr1);
  rot_mul(mcrotaArmForNeutronPropState_7, mctr1, mcrotrArmForNeutronPropState_7);
  mctc1 = coords_set(
#line 879 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    0,
#line 879 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    0,
#line 879 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    ArmMidOnePos);
#line 9777 "./ESS_2001_bispectral.c"
  rot_transpose(mcrotaOrigin, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaArmForNeutronPropState_7 = coords_add(mcposaOrigin, mctc2);
  mctc1 = coords_sub(mcposaguide_Left, mcposaArmForNeutronPropState_7);
  mcposrArmForNeutronPropState_7 = rot_apply(mcrotaArmForNeutronPropState_7, mctc1);
  mcDEBUG_COMPONENT("ArmForNeutronPropState_7", mcposaArmForNeutronPropState_7, mcrotaArmForNeutronPropState_7)
  mccomp_posa[22] = mcposaArmForNeutronPropState_7;
  mccomp_posr[22] = mcposrArmForNeutronPropState_7;
  mcNCounter[22]  = mcPCounter[22] = mcP2Counter[22] = 0;
  mcAbsorbProp[22]= 0;
    /* Component ArmMidTwo. */
  /* Setting parameters for component ArmMidTwo. */
  SIG_MESSAGE("ArmMidTwo (Init:SetPar)");

  SIG_MESSAGE("ArmMidTwo (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 9797 "./ESS_2001_bispectral.c"
  rot_mul(mctr1, mcrotaArmMidOne, mcrotaArmMidTwo);
  rot_transpose(mcrotaArmForNeutronPropState_7, mctr1);
  rot_mul(mcrotaArmMidTwo, mctr1, mcrotrArmMidTwo);
  mctc1 = coords_set(
#line 918 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    0,
#line 918 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    0,
#line 918 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    0);
#line 9808 "./ESS_2001_bispectral.c"
  rot_transpose(mcrotaArmMidOne, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaArmMidTwo = coords_add(mcposaArmMidOne, mctc2);
  mctc1 = coords_sub(mcposaArmForNeutronPropState_7, mcposaArmMidTwo);
  mcposrArmMidTwo = rot_apply(mcrotaArmMidTwo, mctc1);
  mcDEBUG_COMPONENT("ArmMidTwo", mcposaArmMidTwo, mcrotaArmMidTwo)
  mccomp_posa[23] = mcposaArmMidTwo;
  mccomp_posr[23] = mcposrArmMidTwo;
  mcNCounter[23]  = mcPCounter[23] = mcP2Counter[23] = 0;
  mcAbsorbProp[23]= 0;
    /* Component ArmForNeutronPropState_8. */
  /* Setting parameters for component ArmForNeutronPropState_8. */
  SIG_MESSAGE("ArmForNeutronPropState_8 (Init:SetPar)");

  SIG_MESSAGE("ArmForNeutronPropState_8 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 9828 "./ESS_2001_bispectral.c"
  rot_mul(mctr1, mcrotaOrigin, mcrotaArmForNeutronPropState_8);
  rot_transpose(mcrotaArmMidTwo, mctr1);
  rot_mul(mcrotaArmForNeutronPropState_8, mctr1, mcrotrArmForNeutronPropState_8);
  mctc1 = coords_set(
#line 924 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    0,
#line 924 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    0,
#line 924 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    ArmMidOnePos);
#line 9839 "./ESS_2001_bispectral.c"
  rot_transpose(mcrotaOrigin, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaArmForNeutronPropState_8 = coords_add(mcposaOrigin, mctc2);
  mctc1 = coords_sub(mcposaArmMidTwo, mcposaArmForNeutronPropState_8);
  mcposrArmForNeutronPropState_8 = rot_apply(mcrotaArmForNeutronPropState_8, mctc1);
  mcDEBUG_COMPONENT("ArmForNeutronPropState_8", mcposaArmForNeutronPropState_8, mcrotaArmForNeutronPropState_8)
  mccomp_posa[24] = mcposaArmForNeutronPropState_8;
  mccomp_posr[24] = mcposrArmForNeutronPropState_8;
  mcNCounter[24]  = mcPCounter[24] = mcP2Counter[24] = 0;
  mcAbsorbProp[24]= 0;
    /* Component ArmMidThree. */
  /* Setting parameters for component ArmMidThree. */
  SIG_MESSAGE("ArmMidThree (Init:SetPar)");

  SIG_MESSAGE("ArmMidThree (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 9859 "./ESS_2001_bispectral.c"
  rot_mul(mctr1, mcrotaArmMidOne, mcrotaArmMidThree);
  rot_transpose(mcrotaArmForNeutronPropState_8, mctr1);
  rot_mul(mcrotaArmMidThree, mctr1, mcrotrArmMidThree);
  mctc1 = coords_set(
#line 952 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    0,
#line 952 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    0,
#line 952 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    0);
#line 9870 "./ESS_2001_bispectral.c"
  rot_transpose(mcrotaArmMidOne, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaArmMidThree = coords_add(mcposaArmMidOne, mctc2);
  mctc1 = coords_sub(mcposaArmForNeutronPropState_8, mcposaArmMidThree);
  mcposrArmMidThree = rot_apply(mcrotaArmMidThree, mctc1);
  mcDEBUG_COMPONENT("ArmMidThree", mcposaArmMidThree, mcrotaArmMidThree)
  mccomp_posa[25] = mcposaArmMidThree;
  mccomp_posr[25] = mcposrArmMidThree;
  mcNCounter[25]  = mcPCounter[25] = mcP2Counter[25] = 0;
  mcAbsorbProp[25]= 0;
    /* Component ArmExit. */
  /* Setting parameters for component ArmExit. */
  SIG_MESSAGE("ArmExit (Init:SetPar)");

  SIG_MESSAGE("ArmExit (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 9890 "./ESS_2001_bispectral.c"
  rot_mul(mctr1, mcrotaOrigin, mcrotaArmExit);
  rot_transpose(mcrotaArmMidThree, mctr1);
  rot_mul(mcrotaArmExit, mctr1, mcrotrArmExit);
  mctc1 = coords_set(
#line 955 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    0,
#line 955 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    0,
#line 955 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
    ArmExitPos);
#line 9901 "./ESS_2001_bispectral.c"
  rot_transpose(mcrotaOrigin, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaArmExit = coords_add(mcposaOrigin, mctc2);
  mctc1 = coords_sub(mcposaArmMidThree, mcposaArmExit);
  mcposrArmExit = rot_apply(mcrotaArmExit, mctc1);
  mcDEBUG_COMPONENT("ArmExit", mcposaArmExit, mcrotaArmExit)
  mccomp_posa[26] = mcposaArmExit;
  mccomp_posr[26] = mcposrArmExit;
  mcNCounter[26]  = mcPCounter[26] = mcP2Counter[26] = 0;
  mcAbsorbProp[26]= 0;
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
#line 9938 "./ESS_2001_bispectral.c"
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

  /* Initializations for component ArmForGuideRight. */
  SIG_MESSAGE("ArmForGuideRight (Init)");

  /* Initializations for component ArmForGuideBottom. */
  SIG_MESSAGE("ArmForGuideBottom (Init)");

  /* Initializations for component ArmForGuideTop. */
  SIG_MESSAGE("ArmForGuideTop (Init)");

  /* Initializations for component ArmForGuideLeft. */
  SIG_MESSAGE("ArmForGuideLeft (Init)");

  /* Initializations for component cold_source. */
  SIG_MESSAGE("cold_source (Init)");
#define mccompcurname  cold_source
#define mccompcurtype  ESS_moderator_long_2001
#define mccompcurindex 6
#define l_range mcccold_source_l_range
#define w_mult mcccold_source_w_mult
#define size mcccold_source_size
#define l_low mcccold_source_l_low
#define l_high mcccold_source_l_high
#define dist mcccold_source_dist
#define xw mcccold_source_xw
#define yh mcccold_source_yh
#define freq mcccold_source_freq
#define T mcccold_source_T
#define tau mcccold_source_tau
#define tau1 mcccold_source_tau1
#define tau2 mcccold_source_tau2
#define d mcccold_source_d
#define n mcccold_source_n
#define n2 mcccold_source_n2
#define chi2 mcccold_source_chi2
#define I0 mcccold_source_I0
#define I2 mcccold_source_I2
#define branch1 mcccold_source_branch1
#define branch2 mcccold_source_branch2
#define branch_tail mcccold_source_branch_tail
#define twopulses mcccold_source_twopulses
#define target_index mcccold_source_target_index
#line 114 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../obsolete/ESS_moderator_long_2001.comp"
{
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
}
#line 10038 "./ESS_2001_bispectral.c"
#undef target_index
#undef twopulses
#undef branch_tail
#undef branch2
#undef branch1
#undef I2
#undef I0
#undef chi2
#undef n2
#undef n
#undef d
#undef tau2
#undef tau1
#undef tau
#undef T
#undef freq
#undef yh
#undef xw
#undef dist
#undef l_high
#undef l_low
#undef size
#undef w_mult
#undef l_range
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component thermal_source. */
  SIG_MESSAGE("thermal_source (Init)");
#define mccompcurname  thermal_source
#define mccompcurtype  ESS_moderator_long_2001
#define mccompcurindex 7
#define l_range mccthermal_source_l_range
#define w_mult mccthermal_source_w_mult
#define size mccthermal_source_size
#define l_low mccthermal_source_l_low
#define l_high mccthermal_source_l_high
#define dist mccthermal_source_dist
#define xw mccthermal_source_xw
#define yh mccthermal_source_yh
#define freq mccthermal_source_freq
#define T mccthermal_source_T
#define tau mccthermal_source_tau
#define tau1 mccthermal_source_tau1
#define tau2 mccthermal_source_tau2
#define d mccthermal_source_d
#define n mccthermal_source_n
#define n2 mccthermal_source_n2
#define chi2 mccthermal_source_chi2
#define I0 mccthermal_source_I0
#define I2 mccthermal_source_I2
#define branch1 mccthermal_source_branch1
#define branch2 mccthermal_source_branch2
#define branch_tail mccthermal_source_branch_tail
#define twopulses mccthermal_source_twopulses
#define target_index mccthermal_source_target_index
#line 114 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../obsolete/ESS_moderator_long_2001.comp"
{
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
}
#line 10142 "./ESS_2001_bispectral.c"
#undef target_index
#undef twopulses
#undef branch_tail
#undef branch2
#undef branch1
#undef I2
#undef I0
#undef chi2
#undef n2
#undef n
#undef d
#undef tau2
#undef tau1
#undef tau
#undef T
#undef freq
#undef yh
#undef xw
#undef dist
#undef l_high
#undef l_low
#undef size
#undef w_mult
#undef l_range
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component ColdFocus. */
  SIG_MESSAGE("ColdFocus (Init)");

  /* Initializations for component ArmMidOne. */
  SIG_MESSAGE("ArmMidOne (Init)");

  /* Initializations for component mirror_full_center. */
  SIG_MESSAGE("mirror_full_center (Init)");
#define mccompcurname  mirror_full_center
#define mccompcurtype  Mirror_Curved_Bispectral
#define mccompcurindex 10
#define reflect mccmirror_full_center_reflect
#define pTable mccmirror_full_center_pTable
#define focus_s mccmirror_full_center_focus_s
#define focus_e mccmirror_full_center_focus_e
#define mirror_start mccmirror_full_center_mirror_start
#define guide_start mccmirror_full_center_guide_start
#define yheight mccmirror_full_center_yheight
#define smallaxis mccmirror_full_center_smallaxis
#define length mccmirror_full_center_length
#define m mccmirror_full_center_m
#define transmit mccmirror_full_center_transmit
#define substrate_thickness mccmirror_full_center_substrate_thickness
#define coating_thickness mccmirror_full_center_coating_thickness
#define theta_1 mccmirror_full_center_theta_1
#define theta_2 mccmirror_full_center_theta_2
#define theta_3 mccmirror_full_center_theta_3
#line 129 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../contrib/Mirror_Curved_Bispectral.comp"
{
  if (reflect && strlen(reflect) && strcmp(reflect,"NULL") && strcmp(reflect,"0")) {
    if (Table_Read(&pTable, reflect, 1) <= 0) /* read 1st block data from file into pTable */
      exit(fprintf(stderr,"Mirror: %s: can not read file %s\n", NAME_CURRENT_COMP, reflect));
  }
}
#line 10205 "./ESS_2001_bispectral.c"
#undef theta_3
#undef theta_2
#undef theta_1
#undef coating_thickness
#undef substrate_thickness
#undef transmit
#undef m
#undef length
#undef smallaxis
#undef yheight
#undef guide_start
#undef mirror_start
#undef focus_e
#undef focus_s
#undef pTable
#undef reflect
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component ArmForNeutronPropState_2. */
  SIG_MESSAGE("ArmForNeutronPropState_2 (Init)");

  /* Initializations for component guide_right. */
  SIG_MESSAGE("guide_right (Init)");
#define mccompcurname  guide_right
#define mccompcurtype  Mirror_Elliptic_Bispectral
#define mccompcurindex 12
#define reflect mccguide_right_reflect
#define pTable mccguide_right_pTable
#define focus_start_w mccguide_right_focus_start_w
#define focus_end_w mccguide_right_focus_end_w
#define focus_start_h mccguide_right_focus_start_h
#define focus_end_h mccguide_right_focus_end_h
#define mirror_start mccguide_right_mirror_start
#define m mccguide_right_m
#define smallaxis_w mccguide_right_smallaxis_w
#define smallaxis_h mccguide_right_smallaxis_h
#define length mccguide_right_length
#define transmit mccguide_right_transmit
#define substrate_thickness mccguide_right_substrate_thickness
#define coating_thickness mccguide_right_coating_thickness
#line 151 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../contrib/Mirror_Elliptic_Bispectral.comp"
{
  if (reflect && strlen(reflect) && strcmp(reflect,"NULL") && strcmp(reflect,"0")) {
    if (Table_Read(&pTable, reflect, 1) <= 0) /* read 1st block data from file into pTable */
      exit(fprintf(stderr,"Mirror: %s: can not read file %s\n", NAME_CURRENT_COMP, reflect));
  }
}
#line 10255 "./ESS_2001_bispectral.c"
#undef coating_thickness
#undef substrate_thickness
#undef transmit
#undef length
#undef smallaxis_h
#undef smallaxis_w
#undef m
#undef mirror_start
#undef focus_end_h
#undef focus_start_h
#undef focus_end_w
#undef focus_start_w
#undef pTable
#undef reflect
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component ArmForNeutronPropState_4. */
  SIG_MESSAGE("ArmForNeutronPropState_4 (Init)");

  /* Initializations for component guide_bottom. */
  SIG_MESSAGE("guide_bottom (Init)");
#define mccompcurname  guide_bottom
#define mccompcurtype  Mirror_Elliptic_Bispectral
#define mccompcurindex 14
#define reflect mccguide_bottom_reflect
#define pTable mccguide_bottom_pTable
#define focus_start_w mccguide_bottom_focus_start_w
#define focus_end_w mccguide_bottom_focus_end_w
#define focus_start_h mccguide_bottom_focus_start_h
#define focus_end_h mccguide_bottom_focus_end_h
#define mirror_start mccguide_bottom_mirror_start
#define m mccguide_bottom_m
#define smallaxis_w mccguide_bottom_smallaxis_w
#define smallaxis_h mccguide_bottom_smallaxis_h
#define length mccguide_bottom_length
#define transmit mccguide_bottom_transmit
#define substrate_thickness mccguide_bottom_substrate_thickness
#define coating_thickness mccguide_bottom_coating_thickness
#line 151 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../contrib/Mirror_Elliptic_Bispectral.comp"
{
  if (reflect && strlen(reflect) && strcmp(reflect,"NULL") && strcmp(reflect,"0")) {
    if (Table_Read(&pTable, reflect, 1) <= 0) /* read 1st block data from file into pTable */
      exit(fprintf(stderr,"Mirror: %s: can not read file %s\n", NAME_CURRENT_COMP, reflect));
  }
}
#line 10303 "./ESS_2001_bispectral.c"
#undef coating_thickness
#undef substrate_thickness
#undef transmit
#undef length
#undef smallaxis_h
#undef smallaxis_w
#undef m
#undef mirror_start
#undef focus_end_h
#undef focus_start_h
#undef focus_end_w
#undef focus_start_w
#undef pTable
#undef reflect
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component ArmForNeutronPropState_5. */
  SIG_MESSAGE("ArmForNeutronPropState_5 (Init)");

  /* Initializations for component cold_lambda_guidestart. */
  SIG_MESSAGE("cold_lambda_guidestart (Init)");
#define mccompcurname  cold_lambda_guidestart
#define mccompcurtype  L_monitor
#define mccompcurindex 16
#define nL mcccold_lambda_guidestart_nL
#define L_N mcccold_lambda_guidestart_L_N
#define L_p mcccold_lambda_guidestart_L_p
#define L_p2 mcccold_lambda_guidestart_L_p2
#define filename mcccold_lambda_guidestart_filename
#define xmin mcccold_lambda_guidestart_xmin
#define xmax mcccold_lambda_guidestart_xmax
#define ymin mcccold_lambda_guidestart_ymin
#define ymax mcccold_lambda_guidestart_ymax
#define xwidth mcccold_lambda_guidestart_xwidth
#define yheight mcccold_lambda_guidestart_yheight
#define Lmin mcccold_lambda_guidestart_Lmin
#define Lmax mcccold_lambda_guidestart_Lmax
#define restore_neutron mcccold_lambda_guidestart_restore_neutron
#define nowritefile mcccold_lambda_guidestart_nowritefile
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
#line 10366 "./ESS_2001_bispectral.c"
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

  /* Initializations for component thermal_lambda_guidestart. */
  SIG_MESSAGE("thermal_lambda_guidestart (Init)");
#define mccompcurname  thermal_lambda_guidestart
#define mccompcurtype  L_monitor
#define mccompcurindex 17
#define nL mccthermal_lambda_guidestart_nL
#define L_N mccthermal_lambda_guidestart_L_N
#define L_p mccthermal_lambda_guidestart_L_p
#define L_p2 mccthermal_lambda_guidestart_L_p2
#define filename mccthermal_lambda_guidestart_filename
#define xmin mccthermal_lambda_guidestart_xmin
#define xmax mccthermal_lambda_guidestart_xmax
#define ymin mccthermal_lambda_guidestart_ymin
#define ymax mccthermal_lambda_guidestart_ymax
#define xwidth mccthermal_lambda_guidestart_xwidth
#define yheight mccthermal_lambda_guidestart_yheight
#define Lmin mccthermal_lambda_guidestart_Lmin
#define Lmax mccthermal_lambda_guidestart_Lmax
#define restore_neutron mccthermal_lambda_guidestart_restore_neutron
#define nowritefile mccthermal_lambda_guidestart_nowritefile
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
#line 10427 "./ESS_2001_bispectral.c"
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

  /* Initializations for component lambda_guidestart. */
  SIG_MESSAGE("lambda_guidestart (Init)");
#define mccompcurname  lambda_guidestart
#define mccompcurtype  L_monitor
#define mccompcurindex 18
#define nL mcclambda_guidestart_nL
#define L_N mcclambda_guidestart_L_N
#define L_p mcclambda_guidestart_L_p
#define L_p2 mcclambda_guidestart_L_p2
#define filename mcclambda_guidestart_filename
#define xmin mcclambda_guidestart_xmin
#define xmax mcclambda_guidestart_xmax
#define ymin mcclambda_guidestart_ymin
#define ymax mcclambda_guidestart_ymax
#define xwidth mcclambda_guidestart_xwidth
#define yheight mcclambda_guidestart_yheight
#define Lmin mcclambda_guidestart_Lmin
#define Lmax mcclambda_guidestart_Lmax
#define restore_neutron mcclambda_guidestart_restore_neutron
#define nowritefile mcclambda_guidestart_nowritefile
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
#line 10488 "./ESS_2001_bispectral.c"
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

  /* Initializations for component guide_top. */
  SIG_MESSAGE("guide_top (Init)");
#define mccompcurname  guide_top
#define mccompcurtype  Mirror_Elliptic_Bispectral
#define mccompcurindex 19
#define reflect mccguide_top_reflect
#define pTable mccguide_top_pTable
#define focus_start_w mccguide_top_focus_start_w
#define focus_end_w mccguide_top_focus_end_w
#define focus_start_h mccguide_top_focus_start_h
#define focus_end_h mccguide_top_focus_end_h
#define mirror_start mccguide_top_mirror_start
#define m mccguide_top_m
#define smallaxis_w mccguide_top_smallaxis_w
#define smallaxis_h mccguide_top_smallaxis_h
#define length mccguide_top_length
#define transmit mccguide_top_transmit
#define substrate_thickness mccguide_top_substrate_thickness
#define coating_thickness mccguide_top_coating_thickness
#line 151 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../contrib/Mirror_Elliptic_Bispectral.comp"
{
  if (reflect && strlen(reflect) && strcmp(reflect,"NULL") && strcmp(reflect,"0")) {
    if (Table_Read(&pTable, reflect, 1) <= 0) /* read 1st block data from file into pTable */
      exit(fprintf(stderr,"Mirror: %s: can not read file %s\n", NAME_CURRENT_COMP, reflect));
  }
}
#line 10534 "./ESS_2001_bispectral.c"
#undef coating_thickness
#undef substrate_thickness
#undef transmit
#undef length
#undef smallaxis_h
#undef smallaxis_w
#undef m
#undef mirror_start
#undef focus_end_h
#undef focus_start_h
#undef focus_end_w
#undef focus_start_w
#undef pTable
#undef reflect
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component ArmForNeutronPropState_6. */
  SIG_MESSAGE("ArmForNeutronPropState_6 (Init)");

  /* Initializations for component guide_Left. */
  SIG_MESSAGE("guide_Left (Init)");
#define mccompcurname  guide_Left
#define mccompcurtype  Mirror_Elliptic_Bispectral
#define mccompcurindex 21
#define reflect mccguide_Left_reflect
#define pTable mccguide_Left_pTable
#define focus_start_w mccguide_Left_focus_start_w
#define focus_end_w mccguide_Left_focus_end_w
#define focus_start_h mccguide_Left_focus_start_h
#define focus_end_h mccguide_Left_focus_end_h
#define mirror_start mccguide_Left_mirror_start
#define m mccguide_Left_m
#define smallaxis_w mccguide_Left_smallaxis_w
#define smallaxis_h mccguide_Left_smallaxis_h
#define length mccguide_Left_length
#define transmit mccguide_Left_transmit
#define substrate_thickness mccguide_Left_substrate_thickness
#define coating_thickness mccguide_Left_coating_thickness
#line 151 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../contrib/Mirror_Elliptic_Bispectral.comp"
{
  if (reflect && strlen(reflect) && strcmp(reflect,"NULL") && strcmp(reflect,"0")) {
    if (Table_Read(&pTable, reflect, 1) <= 0) /* read 1st block data from file into pTable */
      exit(fprintf(stderr,"Mirror: %s: can not read file %s\n", NAME_CURRENT_COMP, reflect));
  }
}
#line 10582 "./ESS_2001_bispectral.c"
#undef coating_thickness
#undef substrate_thickness
#undef transmit
#undef length
#undef smallaxis_h
#undef smallaxis_w
#undef m
#undef mirror_start
#undef focus_end_h
#undef focus_start_h
#undef focus_end_w
#undef focus_start_w
#undef pTable
#undef reflect
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component ArmForNeutronPropState_7. */
  SIG_MESSAGE("ArmForNeutronPropState_7 (Init)");

  /* Initializations for component ArmMidTwo. */
  SIG_MESSAGE("ArmMidTwo (Init)");

  /* Initializations for component ArmForNeutronPropState_8. */
  SIG_MESSAGE("ArmForNeutronPropState_8 (Init)");

  /* Initializations for component ArmMidThree. */
  SIG_MESSAGE("ArmMidThree (Init)");

  /* Initializations for component ArmExit. */
  SIG_MESSAGE("ArmExit (Init)");

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
#line 10769 "./ESS_2001_bispectral.c"
/* 'Origin=Progress_bar()' component instance extend code */
    SIG_MESSAGE("Origin (Trace:Extend)");
#line 541 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
	if(two_sources!=0){                //switch between sources
		dummy=thermalmultiplier*rand01();
		if (dummy>1){
		flag=1;
		}
	else {flag=0;}
	}
#line 10780 "./ESS_2001_bispectral.c"
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

  /* TRACE Component ArmForGuideRight [2] */
  mccoordschange(mcposrArmForGuideRight, mcrotrArmForGuideRight,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component ArmForGuideRight (without coords transformations) */
  mcJumpTrace_ArmForGuideRight:
  SIG_MESSAGE("ArmForGuideRight (Trace)");
  mcDEBUG_COMP("ArmForGuideRight")
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

#define mcabsorbComp mcabsorbCompArmForGuideRight
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
#define mccompcurname  ArmForGuideRight
#define mccompcurtype  Arm
#define mccompcurindex 2
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompArmForGuideRight:
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

  /* TRACE Component ArmForGuideBottom [3] */
  mccoordschange(mcposrArmForGuideBottom, mcrotrArmForGuideBottom,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component ArmForGuideBottom (without coords transformations) */
  mcJumpTrace_ArmForGuideBottom:
  SIG_MESSAGE("ArmForGuideBottom (Trace)");
  mcDEBUG_COMP("ArmForGuideBottom")
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

#define mcabsorbComp mcabsorbCompArmForGuideBottom
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
#define mccompcurname  ArmForGuideBottom
#define mccompcurtype  Arm
#define mccompcurindex 3
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompArmForGuideBottom:
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

  /* TRACE Component ArmForGuideTop [4] */
  mccoordschange(mcposrArmForGuideTop, mcrotrArmForGuideTop,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component ArmForGuideTop (without coords transformations) */
  mcJumpTrace_ArmForGuideTop:
  SIG_MESSAGE("ArmForGuideTop (Trace)");
  mcDEBUG_COMP("ArmForGuideTop")
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

#define mcabsorbComp mcabsorbCompArmForGuideTop
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
#define mccompcurname  ArmForGuideTop
#define mccompcurtype  Arm
#define mccompcurindex 4
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompArmForGuideTop:
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

  /* TRACE Component ArmForGuideLeft [5] */
  mccoordschange(mcposrArmForGuideLeft, mcrotrArmForGuideLeft,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component ArmForGuideLeft (without coords transformations) */
  mcJumpTrace_ArmForGuideLeft:
  SIG_MESSAGE("ArmForGuideLeft (Trace)");
  mcDEBUG_COMP("ArmForGuideLeft")
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

#define mcabsorbComp mcabsorbCompArmForGuideLeft
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
#define mccompcurname  ArmForGuideLeft
#define mccompcurtype  Arm
#define mccompcurindex 5
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompArmForGuideLeft:
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

  /* TRACE Component cold_source [6] */
  mccoordschange(mcposrcold_source, mcrotrcold_source,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component cold_source (without coords transformations) */
  mcJumpTrace_cold_source:
  SIG_MESSAGE("cold_source (Trace)");
  mcDEBUG_COMP("cold_source")
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

#define mcabsorbComp mcabsorbCompcold_source
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
#define mccompcurname  cold_source
#define mccompcurtype  ESS_moderator_long_2001
#define mccompcurindex 6
#define l_range mcccold_source_l_range
#define w_mult mcccold_source_w_mult
{   /* Declarations of cold_source=ESS_moderator_long_2001() SETTING parameters. */
MCNUM size = mcccold_source_size;
MCNUM l_low = mcccold_source_l_low;
MCNUM l_high = mcccold_source_l_high;
MCNUM dist = mcccold_source_dist;
MCNUM xw = mcccold_source_xw;
MCNUM yh = mcccold_source_yh;
MCNUM freq = mcccold_source_freq;
MCNUM T = mcccold_source_T;
MCNUM tau = mcccold_source_tau;
MCNUM tau1 = mcccold_source_tau1;
MCNUM tau2 = mcccold_source_tau2;
MCNUM d = mcccold_source_d;
MCNUM n = mcccold_source_n;
MCNUM n2 = mcccold_source_n2;
MCNUM chi2 = mcccold_source_chi2;
MCNUM I0 = mcccold_source_I0;
MCNUM I2 = mcccold_source_I2;
MCNUM branch1 = mcccold_source_branch1;
MCNUM branch2 = mcccold_source_branch2;
MCNUM branch_tail = mcccold_source_branch_tail;
MCNUM twopulses = mcccold_source_twopulses;
int target_index = mcccold_source_target_index;
/* 'cold_source=ESS_moderator_long_2001()' component instance has conditional execution */
if (( flag == 1 ))

#line 160 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../obsolete/ESS_moderator_long_2001.comp"
{
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
}
#line 11433 "./ESS_2001_bispectral.c"
/* 'cold_source=ESS_moderator_long_2001()' component instance extend code */
    SIG_MESSAGE("cold_source (Trace:Extend)");
if (( flag == 1 )) {

#line 572 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
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
#line 11452 "./ESS_2001_bispectral.c"
}

}   /* End of cold_source=ESS_moderator_long_2001() SETTING parameter declarations. */
#undef w_mult
#undef l_range
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompcold_source:
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

  /* TRACE Component thermal_source [7] */
  mccoordschange(mcposrthermal_source, mcrotrthermal_source,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component thermal_source (without coords transformations) */
  mcJumpTrace_thermal_source:
  SIG_MESSAGE("thermal_source (Trace)");
  mcDEBUG_COMP("thermal_source")
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

#define mcabsorbComp mcabsorbCompthermal_source
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
#define mccompcurname  thermal_source
#define mccompcurtype  ESS_moderator_long_2001
#define mccompcurindex 7
#define l_range mccthermal_source_l_range
#define w_mult mccthermal_source_w_mult
{   /* Declarations of thermal_source=ESS_moderator_long_2001() SETTING parameters. */
MCNUM size = mccthermal_source_size;
MCNUM l_low = mccthermal_source_l_low;
MCNUM l_high = mccthermal_source_l_high;
MCNUM dist = mccthermal_source_dist;
MCNUM xw = mccthermal_source_xw;
MCNUM yh = mccthermal_source_yh;
MCNUM freq = mccthermal_source_freq;
MCNUM T = mccthermal_source_T;
MCNUM tau = mccthermal_source_tau;
MCNUM tau1 = mccthermal_source_tau1;
MCNUM tau2 = mccthermal_source_tau2;
MCNUM d = mccthermal_source_d;
MCNUM n = mccthermal_source_n;
MCNUM n2 = mccthermal_source_n2;
MCNUM chi2 = mccthermal_source_chi2;
MCNUM I0 = mccthermal_source_I0;
MCNUM I2 = mccthermal_source_I2;
MCNUM branch1 = mccthermal_source_branch1;
MCNUM branch2 = mccthermal_source_branch2;
MCNUM branch_tail = mccthermal_source_branch_tail;
MCNUM twopulses = mccthermal_source_twopulses;
int target_index = mccthermal_source_target_index;
/* 'thermal_source=ESS_moderator_long_2001()' component instance has conditional execution */
if (( flag == 0 ))

#line 160 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../obsolete/ESS_moderator_long_2001.comp"
{
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
}
#line 11692 "./ESS_2001_bispectral.c"
/* 'thermal_source=ESS_moderator_long_2001()' component instance extend code */
    SIG_MESSAGE("thermal_source (Trace:Extend)");
if (( flag == 0 )) {

#line 598 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
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
#line 11709 "./ESS_2001_bispectral.c"
}

}   /* End of thermal_source=ESS_moderator_long_2001() SETTING parameter declarations. */
#undef w_mult
#undef l_range
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompthermal_source:
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

  /* TRACE Component ColdFocus [8] */
  mccoordschange(mcposrColdFocus, mcrotrColdFocus,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component ColdFocus (without coords transformations) */
  mcJumpTrace_ColdFocus:
  SIG_MESSAGE("ColdFocus (Trace)");
  mcDEBUG_COMP("ColdFocus")
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

#define mcabsorbComp mcabsorbCompColdFocus
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
#define mccompcurname  ColdFocus
#define mccompcurtype  Arm
#define mccompcurindex 8
/* 'ColdFocus=Arm()' component instance extend code */
    SIG_MESSAGE("ColdFocus (Trace:Extend)");
#line 616 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
PROP_Z0;
#line 11821 "./ESS_2001_bispectral.c"
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompColdFocus:
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

  /* TRACE Component ArmMidOne [9] */
  mccoordschange(mcposrArmMidOne, mcrotrArmMidOne,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component ArmMidOne (without coords transformations) */
  mcJumpTrace_ArmMidOne:
  SIG_MESSAGE("ArmMidOne (Trace)");
  mcDEBUG_COMP("ArmMidOne")
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

#define mcabsorbComp mcabsorbCompArmMidOne
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
#define mccompcurname  ArmMidOne
#define mccompcurtype  Arm
#define mccompcurindex 9
/* 'ArmMidOne=Arm()' component instance extend code */
    SIG_MESSAGE("ArmMidOne (Trace:Extend)");
#line 627 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
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
#line 11953 "./ESS_2001_bispectral.c"
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompArmMidOne:
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

  /* TRACE Component mirror_full_center [10] */
  mccoordschange(mcposrmirror_full_center, mcrotrmirror_full_center,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component mirror_full_center (without coords transformations) */
  mcJumpTrace_mirror_full_center:
  SIG_MESSAGE("mirror_full_center (Trace)");
  mcDEBUG_COMP("mirror_full_center")
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

#define mcabsorbComp mcabsorbCompmirror_full_center
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
#define mccompcurname  mirror_full_center
#define mccompcurtype  Mirror_Curved_Bispectral
#define mccompcurindex 10
#define reflect mccmirror_full_center_reflect
#define pTable mccmirror_full_center_pTable
{   /* Declarations of mirror_full_center=Mirror_Curved_Bispectral() SETTING parameters. */
MCNUM focus_s = mccmirror_full_center_focus_s;
MCNUM focus_e = mccmirror_full_center_focus_e;
MCNUM mirror_start = mccmirror_full_center_mirror_start;
MCNUM guide_start = mccmirror_full_center_guide_start;
MCNUM yheight = mccmirror_full_center_yheight;
MCNUM smallaxis = mccmirror_full_center_smallaxis;
MCNUM length = mccmirror_full_center_length;
MCNUM m = mccmirror_full_center_m;
MCNUM transmit = mccmirror_full_center_transmit;
MCNUM substrate_thickness = mccmirror_full_center_substrate_thickness;
MCNUM coating_thickness = mccmirror_full_center_coating_thickness;
MCNUM theta_1 = mccmirror_full_center_theta_1;
MCNUM theta_2 = mccmirror_full_center_theta_2;
MCNUM theta_3 = mccmirror_full_center_theta_3;
#line 137 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../contrib/Mirror_Curved_Bispectral.comp"
{
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


}
#line 12574 "./ESS_2001_bispectral.c"
/* 'mirror_full_center=Mirror_Curved_Bispectral()' component instance extend code */
    SIG_MESSAGE("mirror_full_center (Trace:Extend)");
#line 663 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
if (SCATTERED) {
//{printf("I scatter\n");
guide_scatt=5; PROP_DT(1e-9); SCATTER; }
#line 12581 "./ESS_2001_bispectral.c"
}   /* End of mirror_full_center=Mirror_Curved_Bispectral() SETTING parameter declarations. */
#undef pTable
#undef reflect
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompmirror_full_center:
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

  /* TRACE Component ArmForNeutronPropState_2 [11] */
  mccoordschange(mcposrArmForNeutronPropState_2, mcrotrArmForNeutronPropState_2,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component ArmForNeutronPropState_2 (without coords transformations) */
  mcJumpTrace_ArmForNeutronPropState_2:
  SIG_MESSAGE("ArmForNeutronPropState_2 (Trace)");
  mcDEBUG_COMP("ArmForNeutronPropState_2")
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

#define mcabsorbComp mcabsorbCompArmForNeutronPropState_2
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
#define mccompcurname  ArmForNeutronPropState_2
#define mccompcurtype  Arm
#define mccompcurindex 11
/* 'ArmForNeutronPropState_2=Arm()' component instance extend code */
    SIG_MESSAGE("ArmForNeutronPropState_2 (Trace:Extend)");
#line 670 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
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

#line 12719 "./ESS_2001_bispectral.c"
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompArmForNeutronPropState_2:
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

  /* TRACE Component guide_right [12] */
  mccoordschange(mcposrguide_right, mcrotrguide_right,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component guide_right (without coords transformations) */
  mcJumpTrace_guide_right:
  SIG_MESSAGE("guide_right (Trace)");
  mcDEBUG_COMP("guide_right")
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

#define mcabsorbComp mcabsorbCompguide_right
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
#define mccompcurname  guide_right
#define mccompcurtype  Mirror_Elliptic_Bispectral
#define mccompcurindex 12
#define reflect mccguide_right_reflect
#define pTable mccguide_right_pTable
{   /* Declarations of guide_right=Mirror_Elliptic_Bispectral() SETTING parameters. */
MCNUM focus_start_w = mccguide_right_focus_start_w;
MCNUM focus_end_w = mccguide_right_focus_end_w;
MCNUM focus_start_h = mccguide_right_focus_start_h;
MCNUM focus_end_h = mccguide_right_focus_end_h;
MCNUM mirror_start = mccguide_right_mirror_start;
MCNUM m = mccguide_right_m;
MCNUM smallaxis_w = mccguide_right_smallaxis_w;
MCNUM smallaxis_h = mccguide_right_smallaxis_h;
MCNUM length = mccguide_right_length;
MCNUM transmit = mccguide_right_transmit;
MCNUM substrate_thickness = mccguide_right_substrate_thickness;
MCNUM coating_thickness = mccguide_right_coating_thickness;
#line 159 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../contrib/Mirror_Elliptic_Bispectral.comp"
{
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


}
#line 13113 "./ESS_2001_bispectral.c"
/* 'guide_right=Mirror_Elliptic_Bispectral()' component instance extend code */
    SIG_MESSAGE("guide_right (Trace:Extend)");
#line 709 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
if (SCATTERED){
//{printf("I scatter\n");
guide_scatt=3; PROP_DT(1e-9); SCATTER; }
#line 13120 "./ESS_2001_bispectral.c"
}   /* End of guide_right=Mirror_Elliptic_Bispectral() SETTING parameter declarations. */
#undef pTable
#undef reflect
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompguide_right:
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

  /* TRACE Component ArmForNeutronPropState_4 [13] */
  mccoordschange(mcposrArmForNeutronPropState_4, mcrotrArmForNeutronPropState_4,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component ArmForNeutronPropState_4 (without coords transformations) */
  mcJumpTrace_ArmForNeutronPropState_4:
  SIG_MESSAGE("ArmForNeutronPropState_4 (Trace)");
  mcDEBUG_COMP("ArmForNeutronPropState_4")
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

#define mcabsorbComp mcabsorbCompArmForNeutronPropState_4
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
#define mccompcurname  ArmForNeutronPropState_4
#define mccompcurtype  Arm
#define mccompcurindex 13
/* 'ArmForNeutronPropState_4=Arm()' component instance extend code */
    SIG_MESSAGE("ArmForNeutronPropState_4 (Trace:Extend)");
#line 716 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"

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

#line 13263 "./ESS_2001_bispectral.c"
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompArmForNeutronPropState_4:
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

  /* TRACE Component guide_bottom [14] */
  mccoordschange(mcposrguide_bottom, mcrotrguide_bottom,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component guide_bottom (without coords transformations) */
  mcJumpTrace_guide_bottom:
  SIG_MESSAGE("guide_bottom (Trace)");
  mcDEBUG_COMP("guide_bottom")
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

#define mcabsorbComp mcabsorbCompguide_bottom
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
#define mccompcurname  guide_bottom
#define mccompcurtype  Mirror_Elliptic_Bispectral
#define mccompcurindex 14
#define reflect mccguide_bottom_reflect
#define pTable mccguide_bottom_pTable
{   /* Declarations of guide_bottom=Mirror_Elliptic_Bispectral() SETTING parameters. */
MCNUM focus_start_w = mccguide_bottom_focus_start_w;
MCNUM focus_end_w = mccguide_bottom_focus_end_w;
MCNUM focus_start_h = mccguide_bottom_focus_start_h;
MCNUM focus_end_h = mccguide_bottom_focus_end_h;
MCNUM mirror_start = mccguide_bottom_mirror_start;
MCNUM m = mccguide_bottom_m;
MCNUM smallaxis_w = mccguide_bottom_smallaxis_w;
MCNUM smallaxis_h = mccguide_bottom_smallaxis_h;
MCNUM length = mccguide_bottom_length;
MCNUM transmit = mccguide_bottom_transmit;
MCNUM substrate_thickness = mccguide_bottom_substrate_thickness;
MCNUM coating_thickness = mccguide_bottom_coating_thickness;
#line 159 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../contrib/Mirror_Elliptic_Bispectral.comp"
{
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


}
#line 13657 "./ESS_2001_bispectral.c"
/* 'guide_bottom=Mirror_Elliptic_Bispectral()' component instance extend code */
    SIG_MESSAGE("guide_bottom (Trace:Extend)");
#line 759 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
if (SCATTERED){
//printf("I scatter\n");
guide_scatt=2; PROP_DT(1e-9); SCATTER; }
#line 13664 "./ESS_2001_bispectral.c"
}   /* End of guide_bottom=Mirror_Elliptic_Bispectral() SETTING parameter declarations. */
#undef pTable
#undef reflect
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompguide_bottom:
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

  /* TRACE Component ArmForNeutronPropState_5 [15] */
  mccoordschange(mcposrArmForNeutronPropState_5, mcrotrArmForNeutronPropState_5,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component ArmForNeutronPropState_5 (without coords transformations) */
  mcJumpTrace_ArmForNeutronPropState_5:
  SIG_MESSAGE("ArmForNeutronPropState_5 (Trace)");
  mcDEBUG_COMP("ArmForNeutronPropState_5")
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

#define mcabsorbComp mcabsorbCompArmForNeutronPropState_5
  STORE_NEUTRON(15,
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
  mcNCounter[15]++;
  mcPCounter[15] += p;
  mcP2Counter[15] += p*p;
#define mccompcurname  ArmForNeutronPropState_5
#define mccompcurtype  Arm
#define mccompcurindex 15
/* 'ArmForNeutronPropState_5=Arm()' component instance extend code */
    SIG_MESSAGE("ArmForNeutronPropState_5 (Trace:Extend)");
#line 766 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"

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

#line 13805 "./ESS_2001_bispectral.c"
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompArmForNeutronPropState_5:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(15,
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

  /* TRACE Component cold_lambda_guidestart [16] */
  mccoordschange(mcposrcold_lambda_guidestart, mcrotrcold_lambda_guidestart,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component cold_lambda_guidestart (without coords transformations) */
  mcJumpTrace_cold_lambda_guidestart:
  SIG_MESSAGE("cold_lambda_guidestart (Trace)");
  mcDEBUG_COMP("cold_lambda_guidestart")
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

#define mcabsorbComp mcabsorbCompcold_lambda_guidestart
  STORE_NEUTRON(16,
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
  mcNCounter[16]++;
  mcPCounter[16] += p;
  mcP2Counter[16] += p*p;
#define mccompcurname  cold_lambda_guidestart
#define mccompcurtype  L_monitor
#define mccompcurindex 16
#define nL mcccold_lambda_guidestart_nL
#define L_N mcccold_lambda_guidestart_L_N
#define L_p mcccold_lambda_guidestart_L_p
#define L_p2 mcccold_lambda_guidestart_L_p2
{   /* Declarations of cold_lambda_guidestart=L_monitor() SETTING parameters. */
char* filename = mcccold_lambda_guidestart_filename;
MCNUM xmin = mcccold_lambda_guidestart_xmin;
MCNUM xmax = mcccold_lambda_guidestart_xmax;
MCNUM ymin = mcccold_lambda_guidestart_ymin;
MCNUM ymax = mcccold_lambda_guidestart_ymax;
MCNUM xwidth = mcccold_lambda_guidestart_xwidth;
MCNUM yheight = mcccold_lambda_guidestart_yheight;
MCNUM Lmin = mcccold_lambda_guidestart_Lmin;
MCNUM Lmax = mcccold_lambda_guidestart_Lmax;
MCNUM restore_neutron = mcccold_lambda_guidestart_restore_neutron;
int nowritefile = mcccold_lambda_guidestart_nowritefile;
/* 'cold_lambda_guidestart=L_monitor()' component instance has conditional execution */
if (( flag == 1 ))

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
#line 13949 "./ESS_2001_bispectral.c"
}   /* End of cold_lambda_guidestart=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompcold_lambda_guidestart:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(16,
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

  /* TRACE Component thermal_lambda_guidestart [17] */
  mccoordschange(mcposrthermal_lambda_guidestart, mcrotrthermal_lambda_guidestart,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component thermal_lambda_guidestart (without coords transformations) */
  mcJumpTrace_thermal_lambda_guidestart:
  SIG_MESSAGE("thermal_lambda_guidestart (Trace)");
  mcDEBUG_COMP("thermal_lambda_guidestart")
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

#define mcabsorbComp mcabsorbCompthermal_lambda_guidestart
  STORE_NEUTRON(17,
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
  mcNCounter[17]++;
  mcPCounter[17] += p;
  mcP2Counter[17] += p*p;
#define mccompcurname  thermal_lambda_guidestart
#define mccompcurtype  L_monitor
#define mccompcurindex 17
#define nL mccthermal_lambda_guidestart_nL
#define L_N mccthermal_lambda_guidestart_L_N
#define L_p mccthermal_lambda_guidestart_L_p
#define L_p2 mccthermal_lambda_guidestart_L_p2
{   /* Declarations of thermal_lambda_guidestart=L_monitor() SETTING parameters. */
char* filename = mccthermal_lambda_guidestart_filename;
MCNUM xmin = mccthermal_lambda_guidestart_xmin;
MCNUM xmax = mccthermal_lambda_guidestart_xmax;
MCNUM ymin = mccthermal_lambda_guidestart_ymin;
MCNUM ymax = mccthermal_lambda_guidestart_ymax;
MCNUM xwidth = mccthermal_lambda_guidestart_xwidth;
MCNUM yheight = mccthermal_lambda_guidestart_yheight;
MCNUM Lmin = mccthermal_lambda_guidestart_Lmin;
MCNUM Lmax = mccthermal_lambda_guidestart_Lmax;
MCNUM restore_neutron = mccthermal_lambda_guidestart_restore_neutron;
int nowritefile = mccthermal_lambda_guidestart_nowritefile;
/* 'thermal_lambda_guidestart=L_monitor()' component instance has conditional execution */
if (( flag == 0 ))

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
#line 14098 "./ESS_2001_bispectral.c"
}   /* End of thermal_lambda_guidestart=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompthermal_lambda_guidestart:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(17,
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

  /* TRACE Component lambda_guidestart [18] */
  mccoordschange(mcposrlambda_guidestart, mcrotrlambda_guidestart,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component lambda_guidestart (without coords transformations) */
  mcJumpTrace_lambda_guidestart:
  SIG_MESSAGE("lambda_guidestart (Trace)");
  mcDEBUG_COMP("lambda_guidestart")
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

#define mcabsorbComp mcabsorbComplambda_guidestart
  STORE_NEUTRON(18,
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
  mcNCounter[18]++;
  mcPCounter[18] += p;
  mcP2Counter[18] += p*p;
#define mccompcurname  lambda_guidestart
#define mccompcurtype  L_monitor
#define mccompcurindex 18
#define nL mcclambda_guidestart_nL
#define L_N mcclambda_guidestart_L_N
#define L_p mcclambda_guidestart_L_p
#define L_p2 mcclambda_guidestart_L_p2
{   /* Declarations of lambda_guidestart=L_monitor() SETTING parameters. */
char* filename = mcclambda_guidestart_filename;
MCNUM xmin = mcclambda_guidestart_xmin;
MCNUM xmax = mcclambda_guidestart_xmax;
MCNUM ymin = mcclambda_guidestart_ymin;
MCNUM ymax = mcclambda_guidestart_ymax;
MCNUM xwidth = mcclambda_guidestart_xwidth;
MCNUM yheight = mcclambda_guidestart_yheight;
MCNUM Lmin = mcclambda_guidestart_Lmin;
MCNUM Lmax = mcclambda_guidestart_Lmax;
MCNUM restore_neutron = mcclambda_guidestart_restore_neutron;
int nowritefile = mcclambda_guidestart_nowritefile;
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
#line 14245 "./ESS_2001_bispectral.c"
}   /* End of lambda_guidestart=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbComplambda_guidestart:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(18,
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

  /* TRACE Component guide_top [19] */
  mccoordschange(mcposrguide_top, mcrotrguide_top,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component guide_top (without coords transformations) */
  mcJumpTrace_guide_top:
  SIG_MESSAGE("guide_top (Trace)");
  mcDEBUG_COMP("guide_top")
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

#define mcabsorbComp mcabsorbCompguide_top
  STORE_NEUTRON(19,
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
  mcNCounter[19]++;
  mcPCounter[19] += p;
  mcP2Counter[19] += p*p;
#define mccompcurname  guide_top
#define mccompcurtype  Mirror_Elliptic_Bispectral
#define mccompcurindex 19
#define reflect mccguide_top_reflect
#define pTable mccguide_top_pTable
{   /* Declarations of guide_top=Mirror_Elliptic_Bispectral() SETTING parameters. */
MCNUM focus_start_w = mccguide_top_focus_start_w;
MCNUM focus_end_w = mccguide_top_focus_end_w;
MCNUM focus_start_h = mccguide_top_focus_start_h;
MCNUM focus_end_h = mccguide_top_focus_end_h;
MCNUM mirror_start = mccguide_top_mirror_start;
MCNUM m = mccguide_top_m;
MCNUM smallaxis_w = mccguide_top_smallaxis_w;
MCNUM smallaxis_h = mccguide_top_smallaxis_h;
MCNUM length = mccguide_top_length;
MCNUM transmit = mccguide_top_transmit;
MCNUM substrate_thickness = mccguide_top_substrate_thickness;
MCNUM coating_thickness = mccguide_top_coating_thickness;
#line 159 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../contrib/Mirror_Elliptic_Bispectral.comp"
{
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


}
#line 14644 "./ESS_2001_bispectral.c"
/* 'guide_top=Mirror_Elliptic_Bispectral()' component instance extend code */
    SIG_MESSAGE("guide_top (Trace:Extend)");
#line 824 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
if (SCATTERED){
//{printf("I scatter\n");
guide_scatt=1; PROP_DT(1e-9); SCATTER; }
#line 14651 "./ESS_2001_bispectral.c"
}   /* End of guide_top=Mirror_Elliptic_Bispectral() SETTING parameter declarations. */
#undef pTable
#undef reflect
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompguide_top:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(19,
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

  /* TRACE Component ArmForNeutronPropState_6 [20] */
  mccoordschange(mcposrArmForNeutronPropState_6, mcrotrArmForNeutronPropState_6,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component ArmForNeutronPropState_6 (without coords transformations) */
  mcJumpTrace_ArmForNeutronPropState_6:
  SIG_MESSAGE("ArmForNeutronPropState_6 (Trace)");
  mcDEBUG_COMP("ArmForNeutronPropState_6")
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

#define mcabsorbComp mcabsorbCompArmForNeutronPropState_6
  STORE_NEUTRON(20,
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
  mcNCounter[20]++;
  mcPCounter[20] += p;
  mcP2Counter[20] += p*p;
#define mccompcurname  ArmForNeutronPropState_6
#define mccompcurtype  Arm
#define mccompcurindex 20
/* 'ArmForNeutronPropState_6=Arm()' component instance extend code */
    SIG_MESSAGE("ArmForNeutronPropState_6 (Trace:Extend)");
#line 831 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"

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


#line 14793 "./ESS_2001_bispectral.c"
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompArmForNeutronPropState_6:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(20,
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

  /* TRACE Component guide_Left [21] */
  mccoordschange(mcposrguide_Left, mcrotrguide_Left,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component guide_Left (without coords transformations) */
  mcJumpTrace_guide_Left:
  SIG_MESSAGE("guide_Left (Trace)");
  mcDEBUG_COMP("guide_Left")
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

#define mcabsorbComp mcabsorbCompguide_Left
  STORE_NEUTRON(21,
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
  mcNCounter[21]++;
  mcPCounter[21] += p;
  mcP2Counter[21] += p*p;
#define mccompcurname  guide_Left
#define mccompcurtype  Mirror_Elliptic_Bispectral
#define mccompcurindex 21
#define reflect mccguide_Left_reflect
#define pTable mccguide_Left_pTable
{   /* Declarations of guide_Left=Mirror_Elliptic_Bispectral() SETTING parameters. */
MCNUM focus_start_w = mccguide_Left_focus_start_w;
MCNUM focus_end_w = mccguide_Left_focus_end_w;
MCNUM focus_start_h = mccguide_Left_focus_start_h;
MCNUM focus_end_h = mccguide_Left_focus_end_h;
MCNUM mirror_start = mccguide_Left_mirror_start;
MCNUM m = mccguide_Left_m;
MCNUM smallaxis_w = mccguide_Left_smallaxis_w;
MCNUM smallaxis_h = mccguide_Left_smallaxis_h;
MCNUM length = mccguide_Left_length;
MCNUM transmit = mccguide_Left_transmit;
MCNUM substrate_thickness = mccguide_Left_substrate_thickness;
MCNUM coating_thickness = mccguide_Left_coating_thickness;
#line 159 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../contrib/Mirror_Elliptic_Bispectral.comp"
{
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


}
#line 15187 "./ESS_2001_bispectral.c"
/* 'guide_Left=Mirror_Elliptic_Bispectral()' component instance extend code */
    SIG_MESSAGE("guide_Left (Trace:Extend)");
#line 874 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
if (SCATTERED){
//{printf("I scatter\n");
guide_scatt=4; PROP_DT(1e-9); SCATTER; }
#line 15194 "./ESS_2001_bispectral.c"
}   /* End of guide_Left=Mirror_Elliptic_Bispectral() SETTING parameter declarations. */
#undef pTable
#undef reflect
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompguide_Left:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(21,
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

  /* TRACE Component ArmForNeutronPropState_7 [22] */
  mccoordschange(mcposrArmForNeutronPropState_7, mcrotrArmForNeutronPropState_7,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component ArmForNeutronPropState_7 (without coords transformations) */
  mcJumpTrace_ArmForNeutronPropState_7:
  SIG_MESSAGE("ArmForNeutronPropState_7 (Trace)");
  mcDEBUG_COMP("ArmForNeutronPropState_7")
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

#define mcabsorbComp mcabsorbCompArmForNeutronPropState_7
  STORE_NEUTRON(22,
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
  mcNCounter[22]++;
  mcPCounter[22] += p;
  mcP2Counter[22] += p*p;
#define mccompcurname  ArmForNeutronPropState_7
#define mccompcurtype  Arm
#define mccompcurindex 22
/* 'ArmForNeutronPropState_7=Arm()' component instance extend code */
    SIG_MESSAGE("ArmForNeutronPropState_7 (Trace:Extend)");
#line 881 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"

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


#line 15336 "./ESS_2001_bispectral.c"
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompArmForNeutronPropState_7:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(22,
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

  /* TRACE Component ArmMidTwo [23] */
  mccoordschange(mcposrArmMidTwo, mcrotrArmMidTwo,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component ArmMidTwo (without coords transformations) */
  mcJumpTrace_ArmMidTwo:
  SIG_MESSAGE("ArmMidTwo (Trace)");
  mcDEBUG_COMP("ArmMidTwo")
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

#define mcabsorbComp mcabsorbCompArmMidTwo
  STORE_NEUTRON(23,
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
  mcNCounter[23]++;
  mcPCounter[23] += p;
  mcP2Counter[23] += p*p;
#define mccompcurname  ArmMidTwo
#define mccompcurtype  Arm
#define mccompcurindex 23
/* 'ArmMidTwo=Arm()' component instance extend code */
    SIG_MESSAGE("ArmMidTwo (Trace:Extend)");
#line 920 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
 if (guide_scatt==0)
 {SCATTER; }
#line 15445 "./ESS_2001_bispectral.c"
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompArmMidTwo:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(23,
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

  /* TRACE Component ArmForNeutronPropState_8 [24] */
  mccoordschange(mcposrArmForNeutronPropState_8, mcrotrArmForNeutronPropState_8,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component ArmForNeutronPropState_8 (without coords transformations) */
  mcJumpTrace_ArmForNeutronPropState_8:
  SIG_MESSAGE("ArmForNeutronPropState_8 (Trace)");
  mcDEBUG_COMP("ArmForNeutronPropState_8")
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

#define mcabsorbComp mcabsorbCompArmForNeutronPropState_8
  STORE_NEUTRON(24,
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
  mcNCounter[24]++;
  mcPCounter[24] += p;
  mcP2Counter[24] += p*p;
#define mccompcurname  ArmForNeutronPropState_8
#define mccompcurtype  Arm
#define mccompcurindex 24
/* 'ArmForNeutronPropState_8=Arm()' component instance extend code */
    SIG_MESSAGE("ArmForNeutronPropState_8 (Trace:Extend)");
#line 926 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_2001_bispectral/ESS_2001_bispectral.instr"
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

#line 15576 "./ESS_2001_bispectral.c"
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompArmForNeutronPropState_8:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(24,
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

  /* TRACE Component ArmMidThree [25] */
  mccoordschange(mcposrArmMidThree, mcrotrArmMidThree,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component ArmMidThree (without coords transformations) */
  mcJumpTrace_ArmMidThree:
  SIG_MESSAGE("ArmMidThree (Trace)");
  mcDEBUG_COMP("ArmMidThree")
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

#define mcabsorbComp mcabsorbCompArmMidThree
  STORE_NEUTRON(25,
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
  mcNCounter[25]++;
  mcPCounter[25] += p;
  mcP2Counter[25] += p*p;
#define mccompcurname  ArmMidThree
#define mccompcurtype  Arm
#define mccompcurindex 25
if (( guide_scatt > 0 )) goto mcJumpTrace_ArmMidOne;
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompArmMidThree:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(25,
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

  /* TRACE Component ArmExit [26] */
  mccoordschange(mcposrArmExit, mcrotrArmExit,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component ArmExit (without coords transformations) */
  mcJumpTrace_ArmExit:
  SIG_MESSAGE("ArmExit (Trace)");
  mcDEBUG_COMP("ArmExit")
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

#define mcabsorbComp mcabsorbCompArmExit
  STORE_NEUTRON(26,
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
  mcNCounter[26]++;
  mcPCounter[26] += p;
  mcP2Counter[26] += p*p;
#define mccompcurname  ArmExit
#define mccompcurtype  Arm
#define mccompcurindex 26
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompArmExit:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(26,
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
#line 15891 "./ESS_2001_bispectral.c"
}   /* End of Origin=Progress_bar() SETTING parameter declarations. */
#undef CurrentTime
#undef EndTime
#undef StartTime
#undef IntermediateCnts
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'cold_lambda_guidestart'. */
  SIG_MESSAGE("cold_lambda_guidestart (Save)");
#define mccompcurname  cold_lambda_guidestart
#define mccompcurtype  L_monitor
#define mccompcurindex 16
#define nL mcccold_lambda_guidestart_nL
#define L_N mcccold_lambda_guidestart_L_N
#define L_p mcccold_lambda_guidestart_L_p
#define L_p2 mcccold_lambda_guidestart_L_p2
{   /* Declarations of cold_lambda_guidestart=L_monitor() SETTING parameters. */
char* filename = mcccold_lambda_guidestart_filename;
MCNUM xmin = mcccold_lambda_guidestart_xmin;
MCNUM xmax = mcccold_lambda_guidestart_xmax;
MCNUM ymin = mcccold_lambda_guidestart_ymin;
MCNUM ymax = mcccold_lambda_guidestart_ymax;
MCNUM xwidth = mcccold_lambda_guidestart_xwidth;
MCNUM yheight = mcccold_lambda_guidestart_yheight;
MCNUM Lmin = mcccold_lambda_guidestart_Lmin;
MCNUM Lmax = mcccold_lambda_guidestart_Lmax;
MCNUM restore_neutron = mcccold_lambda_guidestart_restore_neutron;
int nowritefile = mcccold_lambda_guidestart_nowritefile;
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
#line 15934 "./ESS_2001_bispectral.c"
}   /* End of cold_lambda_guidestart=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'thermal_lambda_guidestart'. */
  SIG_MESSAGE("thermal_lambda_guidestart (Save)");
#define mccompcurname  thermal_lambda_guidestart
#define mccompcurtype  L_monitor
#define mccompcurindex 17
#define nL mccthermal_lambda_guidestart_nL
#define L_N mccthermal_lambda_guidestart_L_N
#define L_p mccthermal_lambda_guidestart_L_p
#define L_p2 mccthermal_lambda_guidestart_L_p2
{   /* Declarations of thermal_lambda_guidestart=L_monitor() SETTING parameters. */
char* filename = mccthermal_lambda_guidestart_filename;
MCNUM xmin = mccthermal_lambda_guidestart_xmin;
MCNUM xmax = mccthermal_lambda_guidestart_xmax;
MCNUM ymin = mccthermal_lambda_guidestart_ymin;
MCNUM ymax = mccthermal_lambda_guidestart_ymax;
MCNUM xwidth = mccthermal_lambda_guidestart_xwidth;
MCNUM yheight = mccthermal_lambda_guidestart_yheight;
MCNUM Lmin = mccthermal_lambda_guidestart_Lmin;
MCNUM Lmax = mccthermal_lambda_guidestart_Lmax;
MCNUM restore_neutron = mccthermal_lambda_guidestart_restore_neutron;
int nowritefile = mccthermal_lambda_guidestart_nowritefile;
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
#line 15977 "./ESS_2001_bispectral.c"
}   /* End of thermal_lambda_guidestart=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'lambda_guidestart'. */
  SIG_MESSAGE("lambda_guidestart (Save)");
#define mccompcurname  lambda_guidestart
#define mccompcurtype  L_monitor
#define mccompcurindex 18
#define nL mcclambda_guidestart_nL
#define L_N mcclambda_guidestart_L_N
#define L_p mcclambda_guidestart_L_p
#define L_p2 mcclambda_guidestart_L_p2
{   /* Declarations of lambda_guidestart=L_monitor() SETTING parameters. */
char* filename = mcclambda_guidestart_filename;
MCNUM xmin = mcclambda_guidestart_xmin;
MCNUM xmax = mcclambda_guidestart_xmax;
MCNUM ymin = mcclambda_guidestart_ymin;
MCNUM ymax = mcclambda_guidestart_ymax;
MCNUM xwidth = mcclambda_guidestart_xwidth;
MCNUM yheight = mcclambda_guidestart_yheight;
MCNUM Lmin = mcclambda_guidestart_Lmin;
MCNUM Lmax = mcclambda_guidestart_Lmax;
MCNUM restore_neutron = mcclambda_guidestart_restore_neutron;
int nowritefile = mcclambda_guidestart_nowritefile;
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
#line 16020 "./ESS_2001_bispectral.c"
}   /* End of lambda_guidestart=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef nL
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
#line 16064 "./ESS_2001_bispectral.c"
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
    if (!mcNCounter[2]) fprintf(stderr, "Warning: No neutron could reach Component[2] ArmForGuideRight\n");
    if (mcAbsorbProp[2]) fprintf(stderr, "Warning: %g events were removed in Component[2] ArmForGuideRight=Arm()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[2]);
    if (!mcNCounter[3]) fprintf(stderr, "Warning: No neutron could reach Component[3] ArmForGuideBottom\n");
    if (mcAbsorbProp[3]) fprintf(stderr, "Warning: %g events were removed in Component[3] ArmForGuideBottom=Arm()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[3]);
    if (!mcNCounter[4]) fprintf(stderr, "Warning: No neutron could reach Component[4] ArmForGuideTop\n");
    if (mcAbsorbProp[4]) fprintf(stderr, "Warning: %g events were removed in Component[4] ArmForGuideTop=Arm()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[4]);
    if (!mcNCounter[5]) fprintf(stderr, "Warning: No neutron could reach Component[5] ArmForGuideLeft\n");
    if (mcAbsorbProp[5]) fprintf(stderr, "Warning: %g events were removed in Component[5] ArmForGuideLeft=Arm()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[5]);
    if (!mcNCounter[6]) fprintf(stderr, "Warning: No neutron could reach Component[6] cold_source\n");
    if (mcAbsorbProp[6]) fprintf(stderr, "Warning: %g events were removed in Component[6] cold_source=ESS_moderator_long_2001()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[6]);
    if (!mcNCounter[7]) fprintf(stderr, "Warning: No neutron could reach Component[7] thermal_source\n");
    if (mcAbsorbProp[7]) fprintf(stderr, "Warning: %g events were removed in Component[7] thermal_source=ESS_moderator_long_2001()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[7]);
    if (!mcNCounter[8]) fprintf(stderr, "Warning: No neutron could reach Component[8] ColdFocus\n");
    if (mcAbsorbProp[8]) fprintf(stderr, "Warning: %g events were removed in Component[8] ColdFocus=Arm()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[8]);
    if (!mcNCounter[9]) fprintf(stderr, "Warning: No neutron could reach Component[9] ArmMidOne\n");
    if (mcAbsorbProp[9]) fprintf(stderr, "Warning: %g events were removed in Component[9] ArmMidOne=Arm()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[9]);
    if (!mcNCounter[10]) fprintf(stderr, "Warning: No neutron could reach Component[10] mirror_full_center\n");
    if (mcAbsorbProp[10]) fprintf(stderr, "Warning: %g events were removed in Component[10] mirror_full_center=Mirror_Curved_Bispectral()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[10]);
    if (!mcNCounter[11]) fprintf(stderr, "Warning: No neutron could reach Component[11] ArmForNeutronPropState_2\n");
    if (mcAbsorbProp[11]) fprintf(stderr, "Warning: %g events were removed in Component[11] ArmForNeutronPropState_2=Arm()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[11]);
    if (!mcNCounter[12]) fprintf(stderr, "Warning: No neutron could reach Component[12] guide_right\n");
    if (mcAbsorbProp[12]) fprintf(stderr, "Warning: %g events were removed in Component[12] guide_right=Mirror_Elliptic_Bispectral()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[12]);
    if (!mcNCounter[13]) fprintf(stderr, "Warning: No neutron could reach Component[13] ArmForNeutronPropState_4\n");
    if (mcAbsorbProp[13]) fprintf(stderr, "Warning: %g events were removed in Component[13] ArmForNeutronPropState_4=Arm()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[13]);
    if (!mcNCounter[14]) fprintf(stderr, "Warning: No neutron could reach Component[14] guide_bottom\n");
    if (mcAbsorbProp[14]) fprintf(stderr, "Warning: %g events were removed in Component[14] guide_bottom=Mirror_Elliptic_Bispectral()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[14]);
    if (!mcNCounter[15]) fprintf(stderr, "Warning: No neutron could reach Component[15] ArmForNeutronPropState_5\n");
    if (mcAbsorbProp[15]) fprintf(stderr, "Warning: %g events were removed in Component[15] ArmForNeutronPropState_5=Arm()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[15]);
    if (!mcNCounter[16]) fprintf(stderr, "Warning: No neutron could reach Component[16] cold_lambda_guidestart\n");
    if (mcAbsorbProp[16]) fprintf(stderr, "Warning: %g events were removed in Component[16] cold_lambda_guidestart=L_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[16]);
    if (!mcNCounter[17]) fprintf(stderr, "Warning: No neutron could reach Component[17] thermal_lambda_guidestart\n");
    if (mcAbsorbProp[17]) fprintf(stderr, "Warning: %g events were removed in Component[17] thermal_lambda_guidestart=L_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[17]);
    if (!mcNCounter[18]) fprintf(stderr, "Warning: No neutron could reach Component[18] lambda_guidestart\n");
    if (mcAbsorbProp[18]) fprintf(stderr, "Warning: %g events were removed in Component[18] lambda_guidestart=L_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[18]);
    if (!mcNCounter[19]) fprintf(stderr, "Warning: No neutron could reach Component[19] guide_top\n");
    if (mcAbsorbProp[19]) fprintf(stderr, "Warning: %g events were removed in Component[19] guide_top=Mirror_Elliptic_Bispectral()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[19]);
    if (!mcNCounter[20]) fprintf(stderr, "Warning: No neutron could reach Component[20] ArmForNeutronPropState_6\n");
    if (mcAbsorbProp[20]) fprintf(stderr, "Warning: %g events were removed in Component[20] ArmForNeutronPropState_6=Arm()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[20]);
    if (!mcNCounter[21]) fprintf(stderr, "Warning: No neutron could reach Component[21] guide_Left\n");
    if (mcAbsorbProp[21]) fprintf(stderr, "Warning: %g events were removed in Component[21] guide_Left=Mirror_Elliptic_Bispectral()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[21]);
    if (!mcNCounter[22]) fprintf(stderr, "Warning: No neutron could reach Component[22] ArmForNeutronPropState_7\n");
    if (mcAbsorbProp[22]) fprintf(stderr, "Warning: %g events were removed in Component[22] ArmForNeutronPropState_7=Arm()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[22]);
    if (!mcNCounter[23]) fprintf(stderr, "Warning: No neutron could reach Component[23] ArmMidTwo\n");
    if (mcAbsorbProp[23]) fprintf(stderr, "Warning: %g events were removed in Component[23] ArmMidTwo=Arm()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[23]);
    if (!mcNCounter[24]) fprintf(stderr, "Warning: No neutron could reach Component[24] ArmForNeutronPropState_8\n");
    if (mcAbsorbProp[24]) fprintf(stderr, "Warning: %g events were removed in Component[24] ArmForNeutronPropState_8=Arm()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[24]);
    if (!mcNCounter[25]) fprintf(stderr, "Warning: No neutron could reach Component[25] ArmMidThree\n");
    if (mcAbsorbProp[25]) fprintf(stderr, "Warning: %g events were removed in Component[25] ArmMidThree=Arm()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[25]);
    if (!mcNCounter[26]) fprintf(stderr, "Warning: No neutron could reach Component[26] ArmExit\n");
    if (mcAbsorbProp[26]) fprintf(stderr, "Warning: %g events were removed in Component[26] ArmExit=Arm()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[26]);
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
#line 16160 "./ESS_2001_bispectral.c"
}   /* End of Origin=Progress_bar() SETTING parameter declarations. */
#undef CurrentTime
#undef EndTime
#undef StartTime
#undef IntermediateCnts
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'ArmForGuideRight'. */
  SIG_MESSAGE("ArmForGuideRight (McDisplay)");
  printf("MCDISPLAY: component %s\n", "ArmForGuideRight");
#define mccompcurname  ArmForGuideRight
#define mccompcurtype  Arm
#define mccompcurindex 2
#line 40 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Arm.comp"
{
  /* A bit ugly; hard-coded dimensions. */
  
  line(0,0,0,0.2,0,0);
  line(0,0,0,0,0.2,0);
  line(0,0,0,0,0,0.2);
}
#line 16184 "./ESS_2001_bispectral.c"
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'ArmForGuideBottom'. */
  SIG_MESSAGE("ArmForGuideBottom (McDisplay)");
  printf("MCDISPLAY: component %s\n", "ArmForGuideBottom");
#define mccompcurname  ArmForGuideBottom
#define mccompcurtype  Arm
#define mccompcurindex 3
#line 40 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Arm.comp"
{
  /* A bit ugly; hard-coded dimensions. */
  
  line(0,0,0,0.2,0,0);
  line(0,0,0,0,0.2,0);
  line(0,0,0,0,0,0.2);
}
#line 16203 "./ESS_2001_bispectral.c"
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'ArmForGuideTop'. */
  SIG_MESSAGE("ArmForGuideTop (McDisplay)");
  printf("MCDISPLAY: component %s\n", "ArmForGuideTop");
#define mccompcurname  ArmForGuideTop
#define mccompcurtype  Arm
#define mccompcurindex 4
#line 40 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Arm.comp"
{
  /* A bit ugly; hard-coded dimensions. */
  
  line(0,0,0,0.2,0,0);
  line(0,0,0,0,0.2,0);
  line(0,0,0,0,0,0.2);
}
#line 16222 "./ESS_2001_bispectral.c"
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'ArmForGuideLeft'. */
  SIG_MESSAGE("ArmForGuideLeft (McDisplay)");
  printf("MCDISPLAY: component %s\n", "ArmForGuideLeft");
#define mccompcurname  ArmForGuideLeft
#define mccompcurtype  Arm
#define mccompcurindex 5
#line 40 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Arm.comp"
{
  /* A bit ugly; hard-coded dimensions. */
  
  line(0,0,0,0.2,0,0);
  line(0,0,0,0,0.2,0);
  line(0,0,0,0,0,0.2);
}
#line 16241 "./ESS_2001_bispectral.c"
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'cold_source'. */
  SIG_MESSAGE("cold_source (McDisplay)");
  printf("MCDISPLAY: component %s\n", "cold_source");
#define mccompcurname  cold_source
#define mccompcurtype  ESS_moderator_long_2001
#define mccompcurindex 6
#define l_range mcccold_source_l_range
#define w_mult mcccold_source_w_mult
{   /* Declarations of cold_source=ESS_moderator_long_2001() SETTING parameters. */
MCNUM size = mcccold_source_size;
MCNUM l_low = mcccold_source_l_low;
MCNUM l_high = mcccold_source_l_high;
MCNUM dist = mcccold_source_dist;
MCNUM xw = mcccold_source_xw;
MCNUM yh = mcccold_source_yh;
MCNUM freq = mcccold_source_freq;
MCNUM T = mcccold_source_T;
MCNUM tau = mcccold_source_tau;
MCNUM tau1 = mcccold_source_tau1;
MCNUM tau2 = mcccold_source_tau2;
MCNUM d = mcccold_source_d;
MCNUM n = mcccold_source_n;
MCNUM n2 = mcccold_source_n2;
MCNUM chi2 = mcccold_source_chi2;
MCNUM I0 = mcccold_source_I0;
MCNUM I2 = mcccold_source_I2;
MCNUM branch1 = mcccold_source_branch1;
MCNUM branch2 = mcccold_source_branch2;
MCNUM branch_tail = mcccold_source_branch_tail;
MCNUM twopulses = mcccold_source_twopulses;
int target_index = mcccold_source_target_index;
#line 266 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../obsolete/ESS_moderator_long_2001.comp"
{
  
  rectangle("xy", 0, 0, 0, size, size);
}
#line 16282 "./ESS_2001_bispectral.c"
}   /* End of cold_source=ESS_moderator_long_2001() SETTING parameter declarations. */
#undef w_mult
#undef l_range
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'thermal_source'. */
  SIG_MESSAGE("thermal_source (McDisplay)");
  printf("MCDISPLAY: component %s\n", "thermal_source");
#define mccompcurname  thermal_source
#define mccompcurtype  ESS_moderator_long_2001
#define mccompcurindex 7
#define l_range mccthermal_source_l_range
#define w_mult mccthermal_source_w_mult
{   /* Declarations of thermal_source=ESS_moderator_long_2001() SETTING parameters. */
MCNUM size = mccthermal_source_size;
MCNUM l_low = mccthermal_source_l_low;
MCNUM l_high = mccthermal_source_l_high;
MCNUM dist = mccthermal_source_dist;
MCNUM xw = mccthermal_source_xw;
MCNUM yh = mccthermal_source_yh;
MCNUM freq = mccthermal_source_freq;
MCNUM T = mccthermal_source_T;
MCNUM tau = mccthermal_source_tau;
MCNUM tau1 = mccthermal_source_tau1;
MCNUM tau2 = mccthermal_source_tau2;
MCNUM d = mccthermal_source_d;
MCNUM n = mccthermal_source_n;
MCNUM n2 = mccthermal_source_n2;
MCNUM chi2 = mccthermal_source_chi2;
MCNUM I0 = mccthermal_source_I0;
MCNUM I2 = mccthermal_source_I2;
MCNUM branch1 = mccthermal_source_branch1;
MCNUM branch2 = mccthermal_source_branch2;
MCNUM branch_tail = mccthermal_source_branch_tail;
MCNUM twopulses = mccthermal_source_twopulses;
int target_index = mccthermal_source_target_index;
#line 266 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../obsolete/ESS_moderator_long_2001.comp"
{
  
  rectangle("xy", 0, 0, 0, size, size);
}
#line 16326 "./ESS_2001_bispectral.c"
}   /* End of thermal_source=ESS_moderator_long_2001() SETTING parameter declarations. */
#undef w_mult
#undef l_range
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'ColdFocus'. */
  SIG_MESSAGE("ColdFocus (McDisplay)");
  printf("MCDISPLAY: component %s\n", "ColdFocus");
#define mccompcurname  ColdFocus
#define mccompcurtype  Arm
#define mccompcurindex 8
#line 40 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Arm.comp"
{
  /* A bit ugly; hard-coded dimensions. */
  
  line(0,0,0,0.2,0,0);
  line(0,0,0,0,0.2,0);
  line(0,0,0,0,0,0.2);
}
#line 16348 "./ESS_2001_bispectral.c"
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'ArmMidOne'. */
  SIG_MESSAGE("ArmMidOne (McDisplay)");
  printf("MCDISPLAY: component %s\n", "ArmMidOne");
#define mccompcurname  ArmMidOne
#define mccompcurtype  Arm
#define mccompcurindex 9
#line 40 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Arm.comp"
{
  /* A bit ugly; hard-coded dimensions. */
  
  line(0,0,0,0.2,0,0);
  line(0,0,0,0,0.2,0);
  line(0,0,0,0,0,0.2);
}
#line 16367 "./ESS_2001_bispectral.c"
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'mirror_full_center'. */
  SIG_MESSAGE("mirror_full_center (McDisplay)");
  printf("MCDISPLAY: component %s\n", "mirror_full_center");
#define mccompcurname  mirror_full_center
#define mccompcurtype  Mirror_Curved_Bispectral
#define mccompcurindex 10
#define reflect mccmirror_full_center_reflect
#define pTable mccmirror_full_center_pTable
{   /* Declarations of mirror_full_center=Mirror_Curved_Bispectral() SETTING parameters. */
MCNUM focus_s = mccmirror_full_center_focus_s;
MCNUM focus_e = mccmirror_full_center_focus_e;
MCNUM mirror_start = mccmirror_full_center_mirror_start;
MCNUM guide_start = mccmirror_full_center_guide_start;
MCNUM yheight = mccmirror_full_center_yheight;
MCNUM smallaxis = mccmirror_full_center_smallaxis;
MCNUM length = mccmirror_full_center_length;
MCNUM m = mccmirror_full_center_m;
MCNUM transmit = mccmirror_full_center_transmit;
MCNUM substrate_thickness = mccmirror_full_center_substrate_thickness;
MCNUM coating_thickness = mccmirror_full_center_coating_thickness;
MCNUM theta_1 = mccmirror_full_center_theta_1;
MCNUM theta_2 = mccmirror_full_center_theta_2;
MCNUM theta_3 = mccmirror_full_center_theta_3;
#line 638 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../contrib/Mirror_Curved_Bispectral.comp"
{



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
}
#line 16555 "./ESS_2001_bispectral.c"
}   /* End of mirror_full_center=Mirror_Curved_Bispectral() SETTING parameter declarations. */
#undef pTable
#undef reflect
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'ArmForNeutronPropState_2'. */
  SIG_MESSAGE("ArmForNeutronPropState_2 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "ArmForNeutronPropState_2");
#define mccompcurname  ArmForNeutronPropState_2
#define mccompcurtype  Arm
#define mccompcurindex 11
#line 40 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Arm.comp"
{
  /* A bit ugly; hard-coded dimensions. */
  
  line(0,0,0,0.2,0,0);
  line(0,0,0,0,0.2,0);
  line(0,0,0,0,0,0.2);
}
#line 16577 "./ESS_2001_bispectral.c"
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'guide_right'. */
  SIG_MESSAGE("guide_right (McDisplay)");
  printf("MCDISPLAY: component %s\n", "guide_right");
#define mccompcurname  guide_right
#define mccompcurtype  Mirror_Elliptic_Bispectral
#define mccompcurindex 12
#define reflect mccguide_right_reflect
#define pTable mccguide_right_pTable
{   /* Declarations of guide_right=Mirror_Elliptic_Bispectral() SETTING parameters. */
MCNUM focus_start_w = mccguide_right_focus_start_w;
MCNUM focus_end_w = mccguide_right_focus_end_w;
MCNUM focus_start_h = mccguide_right_focus_start_h;
MCNUM focus_end_h = mccguide_right_focus_end_h;
MCNUM mirror_start = mccguide_right_mirror_start;
MCNUM m = mccguide_right_m;
MCNUM smallaxis_w = mccguide_right_smallaxis_w;
MCNUM smallaxis_h = mccguide_right_smallaxis_h;
MCNUM length = mccguide_right_length;
MCNUM transmit = mccguide_right_transmit;
MCNUM substrate_thickness = mccguide_right_substrate_thickness;
MCNUM coating_thickness = mccguide_right_coating_thickness;
#line 435 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../contrib/Mirror_Elliptic_Bispectral.comp"
{




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
}
#line 16729 "./ESS_2001_bispectral.c"
}   /* End of guide_right=Mirror_Elliptic_Bispectral() SETTING parameter declarations. */
#undef pTable
#undef reflect
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'ArmForNeutronPropState_4'. */
  SIG_MESSAGE("ArmForNeutronPropState_4 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "ArmForNeutronPropState_4");
#define mccompcurname  ArmForNeutronPropState_4
#define mccompcurtype  Arm
#define mccompcurindex 13
#line 40 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Arm.comp"
{
  /* A bit ugly; hard-coded dimensions. */
  
  line(0,0,0,0.2,0,0);
  line(0,0,0,0,0.2,0);
  line(0,0,0,0,0,0.2);
}
#line 16751 "./ESS_2001_bispectral.c"
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'guide_bottom'. */
  SIG_MESSAGE("guide_bottom (McDisplay)");
  printf("MCDISPLAY: component %s\n", "guide_bottom");
#define mccompcurname  guide_bottom
#define mccompcurtype  Mirror_Elliptic_Bispectral
#define mccompcurindex 14
#define reflect mccguide_bottom_reflect
#define pTable mccguide_bottom_pTable
{   /* Declarations of guide_bottom=Mirror_Elliptic_Bispectral() SETTING parameters. */
MCNUM focus_start_w = mccguide_bottom_focus_start_w;
MCNUM focus_end_w = mccguide_bottom_focus_end_w;
MCNUM focus_start_h = mccguide_bottom_focus_start_h;
MCNUM focus_end_h = mccguide_bottom_focus_end_h;
MCNUM mirror_start = mccguide_bottom_mirror_start;
MCNUM m = mccguide_bottom_m;
MCNUM smallaxis_w = mccguide_bottom_smallaxis_w;
MCNUM smallaxis_h = mccguide_bottom_smallaxis_h;
MCNUM length = mccguide_bottom_length;
MCNUM transmit = mccguide_bottom_transmit;
MCNUM substrate_thickness = mccguide_bottom_substrate_thickness;
MCNUM coating_thickness = mccguide_bottom_coating_thickness;
#line 435 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../contrib/Mirror_Elliptic_Bispectral.comp"
{




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
}
#line 16903 "./ESS_2001_bispectral.c"
}   /* End of guide_bottom=Mirror_Elliptic_Bispectral() SETTING parameter declarations. */
#undef pTable
#undef reflect
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'ArmForNeutronPropState_5'. */
  SIG_MESSAGE("ArmForNeutronPropState_5 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "ArmForNeutronPropState_5");
#define mccompcurname  ArmForNeutronPropState_5
#define mccompcurtype  Arm
#define mccompcurindex 15
#line 40 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Arm.comp"
{
  /* A bit ugly; hard-coded dimensions. */
  
  line(0,0,0,0.2,0,0);
  line(0,0,0,0,0.2,0);
  line(0,0,0,0,0,0.2);
}
#line 16925 "./ESS_2001_bispectral.c"
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'cold_lambda_guidestart'. */
  SIG_MESSAGE("cold_lambda_guidestart (McDisplay)");
  printf("MCDISPLAY: component %s\n", "cold_lambda_guidestart");
#define mccompcurname  cold_lambda_guidestart
#define mccompcurtype  L_monitor
#define mccompcurindex 16
#define nL mcccold_lambda_guidestart_nL
#define L_N mcccold_lambda_guidestart_L_N
#define L_p mcccold_lambda_guidestart_L_p
#define L_p2 mcccold_lambda_guidestart_L_p2
{   /* Declarations of cold_lambda_guidestart=L_monitor() SETTING parameters. */
char* filename = mcccold_lambda_guidestart_filename;
MCNUM xmin = mcccold_lambda_guidestart_xmin;
MCNUM xmax = mcccold_lambda_guidestart_xmax;
MCNUM ymin = mcccold_lambda_guidestart_ymin;
MCNUM ymax = mcccold_lambda_guidestart_ymax;
MCNUM xwidth = mcccold_lambda_guidestart_xwidth;
MCNUM yheight = mcccold_lambda_guidestart_yheight;
MCNUM Lmin = mcccold_lambda_guidestart_Lmin;
MCNUM Lmax = mcccold_lambda_guidestart_Lmax;
MCNUM restore_neutron = mcccold_lambda_guidestart_restore_neutron;
int nowritefile = mcccold_lambda_guidestart_nowritefile;
#line 120 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/L_monitor.comp"
{
  
  multiline(5, (double)xmin, (double)ymin, 0.0,
               (double)xmax, (double)ymin, 0.0,
               (double)xmax, (double)ymax, 0.0,
               (double)xmin, (double)ymax, 0.0,
               (double)xmin, (double)ymin, 0.0);
}
#line 16961 "./ESS_2001_bispectral.c"
}   /* End of cold_lambda_guidestart=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'thermal_lambda_guidestart'. */
  SIG_MESSAGE("thermal_lambda_guidestart (McDisplay)");
  printf("MCDISPLAY: component %s\n", "thermal_lambda_guidestart");
#define mccompcurname  thermal_lambda_guidestart
#define mccompcurtype  L_monitor
#define mccompcurindex 17
#define nL mccthermal_lambda_guidestart_nL
#define L_N mccthermal_lambda_guidestart_L_N
#define L_p mccthermal_lambda_guidestart_L_p
#define L_p2 mccthermal_lambda_guidestart_L_p2
{   /* Declarations of thermal_lambda_guidestart=L_monitor() SETTING parameters. */
char* filename = mccthermal_lambda_guidestart_filename;
MCNUM xmin = mccthermal_lambda_guidestart_xmin;
MCNUM xmax = mccthermal_lambda_guidestart_xmax;
MCNUM ymin = mccthermal_lambda_guidestart_ymin;
MCNUM ymax = mccthermal_lambda_guidestart_ymax;
MCNUM xwidth = mccthermal_lambda_guidestart_xwidth;
MCNUM yheight = mccthermal_lambda_guidestart_yheight;
MCNUM Lmin = mccthermal_lambda_guidestart_Lmin;
MCNUM Lmax = mccthermal_lambda_guidestart_Lmax;
MCNUM restore_neutron = mccthermal_lambda_guidestart_restore_neutron;
int nowritefile = mccthermal_lambda_guidestart_nowritefile;
#line 120 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/L_monitor.comp"
{
  
  multiline(5, (double)xmin, (double)ymin, 0.0,
               (double)xmax, (double)ymin, 0.0,
               (double)xmax, (double)ymax, 0.0,
               (double)xmin, (double)ymax, 0.0,
               (double)xmin, (double)ymin, 0.0);
}
#line 17002 "./ESS_2001_bispectral.c"
}   /* End of thermal_lambda_guidestart=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'lambda_guidestart'. */
  SIG_MESSAGE("lambda_guidestart (McDisplay)");
  printf("MCDISPLAY: component %s\n", "lambda_guidestart");
#define mccompcurname  lambda_guidestart
#define mccompcurtype  L_monitor
#define mccompcurindex 18
#define nL mcclambda_guidestart_nL
#define L_N mcclambda_guidestart_L_N
#define L_p mcclambda_guidestart_L_p
#define L_p2 mcclambda_guidestart_L_p2
{   /* Declarations of lambda_guidestart=L_monitor() SETTING parameters. */
char* filename = mcclambda_guidestart_filename;
MCNUM xmin = mcclambda_guidestart_xmin;
MCNUM xmax = mcclambda_guidestart_xmax;
MCNUM ymin = mcclambda_guidestart_ymin;
MCNUM ymax = mcclambda_guidestart_ymax;
MCNUM xwidth = mcclambda_guidestart_xwidth;
MCNUM yheight = mcclambda_guidestart_yheight;
MCNUM Lmin = mcclambda_guidestart_Lmin;
MCNUM Lmax = mcclambda_guidestart_Lmax;
MCNUM restore_neutron = mcclambda_guidestart_restore_neutron;
int nowritefile = mcclambda_guidestart_nowritefile;
#line 120 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/L_monitor.comp"
{
  
  multiline(5, (double)xmin, (double)ymin, 0.0,
               (double)xmax, (double)ymin, 0.0,
               (double)xmax, (double)ymax, 0.0,
               (double)xmin, (double)ymax, 0.0,
               (double)xmin, (double)ymin, 0.0);
}
#line 17043 "./ESS_2001_bispectral.c"
}   /* End of lambda_guidestart=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'guide_top'. */
  SIG_MESSAGE("guide_top (McDisplay)");
  printf("MCDISPLAY: component %s\n", "guide_top");
#define mccompcurname  guide_top
#define mccompcurtype  Mirror_Elliptic_Bispectral
#define mccompcurindex 19
#define reflect mccguide_top_reflect
#define pTable mccguide_top_pTable
{   /* Declarations of guide_top=Mirror_Elliptic_Bispectral() SETTING parameters. */
MCNUM focus_start_w = mccguide_top_focus_start_w;
MCNUM focus_end_w = mccguide_top_focus_end_w;
MCNUM focus_start_h = mccguide_top_focus_start_h;
MCNUM focus_end_h = mccguide_top_focus_end_h;
MCNUM mirror_start = mccguide_top_mirror_start;
MCNUM m = mccguide_top_m;
MCNUM smallaxis_w = mccguide_top_smallaxis_w;
MCNUM smallaxis_h = mccguide_top_smallaxis_h;
MCNUM length = mccguide_top_length;
MCNUM transmit = mccguide_top_transmit;
MCNUM substrate_thickness = mccguide_top_substrate_thickness;
MCNUM coating_thickness = mccguide_top_coating_thickness;
#line 435 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../contrib/Mirror_Elliptic_Bispectral.comp"
{




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
}
#line 17200 "./ESS_2001_bispectral.c"
}   /* End of guide_top=Mirror_Elliptic_Bispectral() SETTING parameter declarations. */
#undef pTable
#undef reflect
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'ArmForNeutronPropState_6'. */
  SIG_MESSAGE("ArmForNeutronPropState_6 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "ArmForNeutronPropState_6");
#define mccompcurname  ArmForNeutronPropState_6
#define mccompcurtype  Arm
#define mccompcurindex 20
#line 40 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Arm.comp"
{
  /* A bit ugly; hard-coded dimensions. */
  
  line(0,0,0,0.2,0,0);
  line(0,0,0,0,0.2,0);
  line(0,0,0,0,0,0.2);
}
#line 17222 "./ESS_2001_bispectral.c"
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'guide_Left'. */
  SIG_MESSAGE("guide_Left (McDisplay)");
  printf("MCDISPLAY: component %s\n", "guide_Left");
#define mccompcurname  guide_Left
#define mccompcurtype  Mirror_Elliptic_Bispectral
#define mccompcurindex 21
#define reflect mccguide_Left_reflect
#define pTable mccguide_Left_pTable
{   /* Declarations of guide_Left=Mirror_Elliptic_Bispectral() SETTING parameters. */
MCNUM focus_start_w = mccguide_Left_focus_start_w;
MCNUM focus_end_w = mccguide_Left_focus_end_w;
MCNUM focus_start_h = mccguide_Left_focus_start_h;
MCNUM focus_end_h = mccguide_Left_focus_end_h;
MCNUM mirror_start = mccguide_Left_mirror_start;
MCNUM m = mccguide_Left_m;
MCNUM smallaxis_w = mccguide_Left_smallaxis_w;
MCNUM smallaxis_h = mccguide_Left_smallaxis_h;
MCNUM length = mccguide_Left_length;
MCNUM transmit = mccguide_Left_transmit;
MCNUM substrate_thickness = mccguide_Left_substrate_thickness;
MCNUM coating_thickness = mccguide_Left_coating_thickness;
#line 435 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../contrib/Mirror_Elliptic_Bispectral.comp"
{




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
}
#line 17374 "./ESS_2001_bispectral.c"
}   /* End of guide_Left=Mirror_Elliptic_Bispectral() SETTING parameter declarations. */
#undef pTable
#undef reflect
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'ArmForNeutronPropState_7'. */
  SIG_MESSAGE("ArmForNeutronPropState_7 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "ArmForNeutronPropState_7");
#define mccompcurname  ArmForNeutronPropState_7
#define mccompcurtype  Arm
#define mccompcurindex 22
#line 40 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Arm.comp"
{
  /* A bit ugly; hard-coded dimensions. */
  
  line(0,0,0,0.2,0,0);
  line(0,0,0,0,0.2,0);
  line(0,0,0,0,0,0.2);
}
#line 17396 "./ESS_2001_bispectral.c"
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'ArmMidTwo'. */
  SIG_MESSAGE("ArmMidTwo (McDisplay)");
  printf("MCDISPLAY: component %s\n", "ArmMidTwo");
#define mccompcurname  ArmMidTwo
#define mccompcurtype  Arm
#define mccompcurindex 23
#line 40 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Arm.comp"
{
  /* A bit ugly; hard-coded dimensions. */
  
  line(0,0,0,0.2,0,0);
  line(0,0,0,0,0.2,0);
  line(0,0,0,0,0,0.2);
}
#line 17415 "./ESS_2001_bispectral.c"
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'ArmForNeutronPropState_8'. */
  SIG_MESSAGE("ArmForNeutronPropState_8 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "ArmForNeutronPropState_8");
#define mccompcurname  ArmForNeutronPropState_8
#define mccompcurtype  Arm
#define mccompcurindex 24
#line 40 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Arm.comp"
{
  /* A bit ugly; hard-coded dimensions. */
  
  line(0,0,0,0.2,0,0);
  line(0,0,0,0,0.2,0);
  line(0,0,0,0,0,0.2);
}
#line 17434 "./ESS_2001_bispectral.c"
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'ArmMidThree'. */
  SIG_MESSAGE("ArmMidThree (McDisplay)");
  printf("MCDISPLAY: component %s\n", "ArmMidThree");
#define mccompcurname  ArmMidThree
#define mccompcurtype  Arm
#define mccompcurindex 25
#line 40 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Arm.comp"
{
  /* A bit ugly; hard-coded dimensions. */
  
  line(0,0,0,0.2,0,0);
  line(0,0,0,0,0.2,0);
  line(0,0,0,0,0,0.2);
}
#line 17453 "./ESS_2001_bispectral.c"
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'ArmExit'. */
  SIG_MESSAGE("ArmExit (McDisplay)");
  printf("MCDISPLAY: component %s\n", "ArmExit");
#define mccompcurname  ArmExit
#define mccompcurtype  Arm
#define mccompcurindex 26
#line 40 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Arm.comp"
{
  /* A bit ugly; hard-coded dimensions. */
  
  line(0,0,0,0.2,0,0);
  line(0,0,0,0,0.2,0);
  line(0,0,0,0,0,0.2);
}
#line 17472 "./ESS_2001_bispectral.c"
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
/* end of generated C code ./ESS_2001_bispectral.c */
