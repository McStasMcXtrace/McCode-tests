/* Automatically generated file. Do not edit. 
 * Format:     ANSI C source code
 * Creator:    McStas <http://www.mcstas.org>
 * Instrument: /zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr (Union_test_mask)
 * Date:       Tue Feb 25 21:11:07 2020
 * File:       ./Union_test_mask.c
 * Compile:    cc -o Union_test_mask.out ./Union_test_mask.c  -I@MCCODE_LIB@/share/
 * CFLAGS= -I@MCCODE_LIB@/share/
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

#line 712 "./Union_test_mask.c"

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

#line 945 "./Union_test_mask.c"

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

#line 4977 "./Union_test_mask.c"

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

#line 5337 "./Union_test_mask.c"
#ifdef MC_TRACE_ENABLED
int mctraceenabled = 1;
#else
int mctraceenabled = 0;
#endif
#define MCSTAS "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../"
int mcdefaultmain = 1;
char mcinstrument_name[] = "Union_test_mask";
char mcinstrument_source[] = "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr";
char *mcinstrument_exe=NULL; /* will be set to argv[0] in main */
int main(int argc, char *argv[]){return mccode_main(argc, argv);}
void mcinit(void);
void mcraytrace(void);
void mcsave(FILE *);
void mcfinally(void);
void mcdisplay(void);

/* Shared user declarations for all components 'Incoherent_process'. */
#line 62 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../contrib/union/Incoherent_process.comp"
#ifndef Union
#define Union $Revision: 0.8 $

#include "Union_functions.c"
#include "Union_initialization.c"

#endif


struct Incoherent_physics_storage_struct{
    // Variables that needs to be transfered between any of the following places:
    // The initialize in this component
    // The function for calculating my
    // The function for calculating scattering
    
    double my_scattering;
    double QE_sampling_frequency;
    double lorentzian_width;
    
};

// Function for calculating my in Incoherent case
int Incoherent_physics_my(double *my,double *k_initial, union data_transfer_union data_transfer, struct focus_data_struct *focus_data) {
    *my = data_transfer.pointer_to_a_Incoherent_physics_storage_struct->my_scattering;
    return 1;
};

// Function for basic incoherent scattering event
int Incoherent_physics_scattering(double *k_final, double *k_initial, double *weight, union data_transfer_union data_transfer, struct focus_data_struct *focus_data) {

    //New version of incoherent scattering
    double k_length = sqrt(k_initial[0]*k_initial[0]+k_initial[1]*k_initial[1]+k_initial[2]*k_initial[2]);

    Coords k_out;
    // Here is the focusing system in action, get a vector
    double solid_angle;
    focus_data->focusing_function(&k_out,&solid_angle,focus_data);
    NORM(k_out.x,k_out.y,k_out.z);
    *weight *= solid_angle*0.25/PI;
    
    double v_i,v_f,E_i,dE,E_f;
    
    if (rand01() < data_transfer.pointer_to_a_Incoherent_physics_storage_struct->QE_sampling_frequency) {
      v_i = k_length * K2V;
      E_i = VS2E*v_i*v_i;
      dE = data_transfer.pointer_to_a_Incoherent_physics_storage_struct->lorentzian_width*tan(PI/2*randpm1());
      E_f = E_i + dE;
      if (E_f <= 0)
        return 0;
      v_f = SE2V*sqrt(E_f);
      k_length = v_f*V2K;
    }
    
    k_final[0] = k_out.x*k_length; k_final[1] = k_out.y*k_length; k_final[2] = k_out.z*k_length;
    return 1;
};

#line 5414 "./Union_test_mask.c"

/* Shared user declarations for all components 'Union_make_material'. */
#line 57 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../contrib/union/Union_make_material.comp"
#ifndef Union
#define Union $Revision: 0.8 $

#include "Union_functions.c"
#include "Union_initialization.c"

#endif

// This function checks if global_process_element should be included in this material when using automatic linking, returns 1 if yes, 0 if no.
int automatic_linking_materials_function(struct global_process_element_struct global_process_element, struct pointer_to_global_material_list global_material_list,int current_index) {
    // Remember this function is used before the current material is added to global_material_list
    // debug info
    //MPI_MASTER(
    //printf("Checking if process with index %d should be automatically linked to material with index %d\n",global_process_element.component_index,current_index);
    //)

    // Check if this is the first make_material, which makes the problem simpler.
    if (global_material_list.num_elements == 0) {
       if (global_process_element.component_index < current_index) return 1;
       else return 0;
    }
    // In case there are more than 1 make_material, global_material_list.elements[global_material_list.num_elements-1].component_index makes sense.
    if (global_process_element.component_index < current_index && global_process_element.component_index > global_material_list.elements[global_material_list.num_elements-1].component_index) return 1;
    else return 0;
}

void manual_linking_function_material(char *input_string, struct pointer_to_global_process_list *global_process_list, struct pointer_to_1d_int_list *accepted_processes, char *component_name) {
    // Need to check a input_string of text for an occurance of name. If it is in the inputstring, yes return 1, otherwise 0.
   char *token;
   int loop_index;
   char local_string[256];
   
   strcpy(local_string,input_string);
   // get the first token
   token = strtok(local_string,",");
   
   // walk through other tokens
   while( token != NULL ) 
   {
      //printf( " %s\n", token );
      for (loop_index=0;loop_index<global_process_list->num_elements;loop_index++) {
        if (strcmp(token,global_process_list->elements[loop_index].name) == 0) {
          add_element_to_int_list(accepted_processes,loop_index);
          break;
        }
        
        if (loop_index == global_process_list->num_elements - 1) {
          // All possible process names have been looked through, and the break was not executed.
          // Alert the user to this problem by showing the process name that was not found and the currently available processes
            printf("\n");
            printf("ERROR: The process string \"%s\" in Union material \"%s\" had an entry that did not match a specified process. \n",input_string,component_name);
            printf("       The unrecoignized process name was: \"%s\" \n",token);
            printf("       The processes available at this point (need to be defined before the material): \n");
            for (loop_index=0;loop_index<global_process_list->num_elements;loop_index++)
              printf("         %s\n",global_process_list->elements[loop_index].name);
            exit(1);
        }
      }
      
      // Updates the token
      token = strtok(NULL,",");
   }
}

// This function is needed in initialize of all geometry components
// Possible to insert these functions in make material, as they are only compiled once instead of many times
int manual_linking_function(char *name, char *input_string) {
    // Need to check a input_string of text for an occurance of name. If it is in the inputstring, yes return 1, otherwise 0.
   char *token;
   int return_integer=0;
   char local_string[124];
   
   strcpy(local_string,input_string);
   /* get the first token */
   token = strtok(local_string,",");
   
   /* walk through other tokens */
   while( token != NULL ) 
   {
      //printf( " %s\n", token );
      if (strcmp(token,name) == 0) return_integer=1;
      
      token = strtok(NULL,",");
   }
   
   return return_integer;
}



/*
int count_commas(char *string) {
  int return_value = 0;
  
  int index;
  for (index=0;index<strlen(string);index++) {
    printf("%c \n",string[index]);
    if (string[index]==',') return_value++;
  }
    
  //printf("number_of_commas = %d \n",return_value);
  return return_value;
}
*/


#line 5524 "./Union_test_mask.c"

/* Shared user declarations for all components 'Union_box'. */
#line 77 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../contrib/union/Union_box.comp"
#ifndef Union
#define Union $Revision: 0.8 $

#include "Union_functions.c"
#include "Union_initialization.c"

#endif


void mcdisplay_box_function(struct lines_to_draw *lines_to_draw_output,int index, struct geometry_struct **Geometries,int number_of_volumes) {
    // Function to call in mcdisplay section of the sample component for this volume
    // One can assume that Volumes[index] refers to a volume with the geometry described in this file
    
    double depth = Geometries[index]->geometry_parameters.p_box_storage->z_depth;
    double width1 = Geometries[index]->geometry_parameters.p_box_storage->x_width1;
    double width2 = Geometries[index]->geometry_parameters.p_box_storage->x_width2;
    double height1 = Geometries[index]->geometry_parameters.p_box_storage->y_height1;
    double height2 = Geometries[index]->geometry_parameters.p_box_storage->y_height2;
    
    Coords x_vector = Geometries[index]->geometry_parameters.p_box_storage->x_vector;
    Coords y_vector = Geometries[index]->geometry_parameters.p_box_storage->y_vector;
    Coords z_vector = Geometries[index]->geometry_parameters.p_box_storage->z_vector;
    
    Coords center = Geometries[index]->center;
    
    Coords square1[4],square2[4];
    
    square1[0] = coords_add(coords_add(coords_add(center,coords_scalar_mult(z_vector,-0.5*depth)),coords_scalar_mult(x_vector,-0.5*width1)),coords_scalar_mult(y_vector,-0.5*height1));
    
    square1[1] = coords_add(square1[0],coords_scalar_mult(x_vector,width1));
    square1[2] = coords_add(square1[1],coords_scalar_mult(y_vector,height1));
    square1[3] = coords_add(square1[0],coords_scalar_mult(y_vector,height1));
    
    square2[0] = coords_add(coords_add(coords_add(center,coords_scalar_mult(z_vector,0.5*depth)),coords_scalar_mult(x_vector,-0.5*width2)),coords_scalar_mult(y_vector,-0.5*height2));
    
    square2[1] = coords_add(square2[0],coords_scalar_mult(x_vector,width2));
    square2[2] = coords_add(square2[1],coords_scalar_mult(y_vector,height2));
    square2[3] = coords_add(square2[0],coords_scalar_mult(y_vector,height2));
    
    struct lines_to_draw lines_to_draw_temp;
    lines_to_draw_temp.number_of_lines = 0;
    
    int iterate;
    for (iterate=0;iterate<3;iterate++) {
        lines_to_draw_temp = draw_line_with_highest_priority(square1[iterate],square1[iterate+1],index,Geometries,number_of_volumes,2);
        merge_lines_to_draw(lines_to_draw_output,&lines_to_draw_temp);
    }
    lines_to_draw_temp = draw_line_with_highest_priority(square1[3],square1[0],index,Geometries,number_of_volumes,2);
    merge_lines_to_draw(lines_to_draw_output,&lines_to_draw_temp);

    for (iterate=0;iterate<3;iterate++) {
        lines_to_draw_temp = draw_line_with_highest_priority(square2[iterate],square2[iterate+1],index,Geometries,number_of_volumes,2);
        merge_lines_to_draw(lines_to_draw_output,&lines_to_draw_temp);
    }
    lines_to_draw_temp = draw_line_with_highest_priority(square2[3],square2[0],index,Geometries,number_of_volumes,2);
    merge_lines_to_draw(lines_to_draw_output,&lines_to_draw_temp);

    for (iterate=0;iterate<4;iterate++) {
        lines_to_draw_temp = draw_line_with_highest_priority(square1[iterate],square2[iterate],index,Geometries,number_of_volumes,2);
        merge_lines_to_draw(lines_to_draw_output,&lines_to_draw_temp);
    }
};

void initialize_box_geometry_from_main_component(struct geometry_struct *box) {
    // Function to be called in initialize of the main component
    // This is done as the rotation matrix needs to be relative to the main component instead of global
    // Everything done in initialize in this component file has the rotation matrix relative to global
    Coords simple_vector = coords_set(1,0,0);
    Coords rotated_vector;
    
    rotated_vector = rot_apply(box->rotation_matrix,simple_vector);
    NORM(rotated_vector.x,rotated_vector.y,rotated_vector.z);
    box->geometry_parameters.p_box_storage->x_vector = rotated_vector;
    
    simple_vector = coords_set(0,1,0);
    rotated_vector = rot_apply(box->rotation_matrix,simple_vector);
    NORM(rotated_vector.x,rotated_vector.y,rotated_vector.z);
    box->geometry_parameters.p_box_storage->y_vector = rotated_vector;
    
    simple_vector = coords_set(0,0,1);
    rotated_vector = rot_apply(box->rotation_matrix,simple_vector);
    NORM(rotated_vector.x,rotated_vector.y,rotated_vector.z);
    box->geometry_parameters.p_box_storage->z_vector = rotated_vector;
};

struct pointer_to_1d_coords_list box_shell_points(struct geometry_struct *geometry,int max_number_of_points) {
  // This function returns an array of corner positions for the box in the main coordinate system.
  // Normally one would limit it to a maximum number of points, but as there are only 8 for the box,
  //  it is hardcoded to 8. Other geometries can be approximated with a variable number of points.
  
  struct pointer_to_1d_coords_list corner_points;
  corner_points.elements = malloc(8*sizeof(Coords));
  corner_points.num_elements = 8;
  
  double depth = geometry->geometry_parameters.p_box_storage->z_depth;
  double width1 = geometry->geometry_parameters.p_box_storage->x_width1;
  double width2 = geometry->geometry_parameters.p_box_storage->x_width2;
  double height1 = geometry->geometry_parameters.p_box_storage->y_height1;
  double height2 = geometry->geometry_parameters.p_box_storage->y_height2;
    
  Coords x_vector = geometry->geometry_parameters.p_box_storage->x_vector;
  Coords y_vector = geometry->geometry_parameters.p_box_storage->y_vector;
  Coords z_vector = geometry->geometry_parameters.p_box_storage->z_vector;
    
  Coords center = geometry->center;
    
  corner_points.elements[0] = coords_add(coords_add(coords_add(center,coords_scalar_mult(z_vector,-0.5*depth)),coords_scalar_mult(x_vector,-0.5*width1)),coords_scalar_mult(y_vector,-0.5*height1));
    
  corner_points.elements[1] = coords_add(corner_points.elements[0],coords_scalar_mult(x_vector,width1));
  corner_points.elements[2] = coords_add(corner_points.elements[1],coords_scalar_mult(y_vector,height1));
  corner_points.elements[3] = coords_add(corner_points.elements[0],coords_scalar_mult(y_vector,height1));
    
  corner_points.elements[4] = coords_add(coords_add(coords_add(center,coords_scalar_mult(z_vector,0.5*depth)),coords_scalar_mult(x_vector,-0.5*width2)),coords_scalar_mult(y_vector,-0.5*height2));
    
  corner_points.elements[5] = coords_add(corner_points.elements[4],coords_scalar_mult(x_vector,width2));
  corner_points.elements[6] = coords_add(corner_points.elements[5],coords_scalar_mult(y_vector,height2));
  corner_points.elements[7] = coords_add(corner_points.elements[4],coords_scalar_mult(y_vector,height2));
  
  return corner_points;

}

#line 5650 "./Union_test_mask.c"

/* Shared user declarations for all components 'Union_master'. */
#line 54 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../contrib/union/Union_master.comp"
#ifndef Union
#define Union $Revision: 0.9 $

#include "Union_functions.c"
#include "Union_initialization.c"

#endif
// TEST
struct logger_with_data_struct loggers_with_data_array;
#line 5663 "./Union_test_mask.c"

/* Instrument parameters. */

#define mcNUMIPAR 0
int mcnumipar = 0;
struct mcinputtable_struct mcinputtable[mcNUMIPAR+1] = {
  NULL, NULL, instr_type_double, ""
};

/* User declarations from instrument definition. */
#define mccompcurname  Union_test_mask
#define mccompcurtype  INSTRUMENT
#define mccompcurindex 0
#define mcposaUnion_test_mask coords_set(0,0,0)
#line 30 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
int scattered_1;
int scattered_2;
int sample_1_index = 1;
int sample_2_index = 2;
#line 5683 "./Union_test_mask.c"
#undef mcposaUnion_test_mask
#undef mccompcurindex
#undef mccompcurtype
#undef mccompcurname

/* neutron state table at each component input (local coords) */
/* [x, y, z, vx, vy, vz, t, sx, sy, sz, p] */
MCNUM mccomp_storein[11*19];
/* Components position table (absolute and relative coords) */
Coords mccomp_posa[19];
Coords mccomp_posr[19];
/* Counter for each comp to check for inactive ones */
MCNUM  mcNCounter[19];
MCNUM  mcPCounter[19];
MCNUM  mcP2Counter[19];
#define mcNUMCOMP 18 /* number of components */
/* Counter for PROP ABSORB */
MCNUM  mcAbsorbProp[19];
/* Flag true when previous component acted on the neutron (SCATTER) */
MCNUM mcScattered=0;
/* Flag true when neutron should be restored (RESTORE) */
MCNUM mcRestore=0;
/* Declarations of component definition and setting parameters. */

/* Setting parameters for component 'test_incoherent' [1]. */
MCNUM mcctest_incoherent_sigma;
MCNUM mcctest_incoherent_f_QE;
MCNUM mcctest_incoherent_gamma;
MCNUM mcctest_incoherent_packing_factor;
MCNUM mcctest_incoherent_unit_cell_volume;
MCNUM mcctest_incoherent_interact_fraction;

/* Setting parameters for component 'test_material' [2]. */
char mcctest_material_process_string[16384];
MCNUM mcctest_material_my_absorption;
MCNUM mcctest_material_absorber;

/* Setting parameters for component 'a1' [3]. */
char mcca1_profile[16384];
MCNUM mcca1_percent;
MCNUM mcca1_flag_save;
MCNUM mcca1_minutes;

/* Setting parameters for component 'source' [4]. */
MCNUM mccsource_xwidth;
MCNUM mccsource_yheight;
MCNUM mccsource_focus_aw;
MCNUM mccsource_focus_ah;
MCNUM mccsource_E0;
MCNUM mccsource_dE;
MCNUM mccsource_lambda0;
MCNUM mccsource_dlambda;
MCNUM mccsource_gauss;
MCNUM mccsource_flux;

/* Setting parameters for component 'volume_1' [6]. */
char mccvolume_1_material_string[16384];
MCNUM mccvolume_1_priority;
MCNUM mccvolume_1_xwidth;
MCNUM mccvolume_1_yheight;
MCNUM mccvolume_1_zdepth;
MCNUM mccvolume_1_xwidth2;
MCNUM mccvolume_1_yheight2;
MCNUM mccvolume_1_visualize;
int mccvolume_1_target_index;
MCNUM mccvolume_1_target_x;
MCNUM mccvolume_1_target_y;
MCNUM mccvolume_1_target_z;
MCNUM mccvolume_1_focus_aw;
MCNUM mccvolume_1_focus_ah;
MCNUM mccvolume_1_focus_xw;
MCNUM mccvolume_1_focus_xh;
MCNUM mccvolume_1_focus_r;
MCNUM mccvolume_1_p_interact;
char mccvolume_1_mask_string[16384];
char mccvolume_1_mask_setting[16384];
MCNUM mccvolume_1_number_of_activations;

/* Setting parameters for component 'volume_2' [7]. */
char mccvolume_2_material_string[16384];
MCNUM mccvolume_2_priority;
MCNUM mccvolume_2_xwidth;
MCNUM mccvolume_2_yheight;
MCNUM mccvolume_2_zdepth;
MCNUM mccvolume_2_xwidth2;
MCNUM mccvolume_2_yheight2;
MCNUM mccvolume_2_visualize;
int mccvolume_2_target_index;
MCNUM mccvolume_2_target_x;
MCNUM mccvolume_2_target_y;
MCNUM mccvolume_2_target_z;
MCNUM mccvolume_2_focus_aw;
MCNUM mccvolume_2_focus_ah;
MCNUM mccvolume_2_focus_xw;
MCNUM mccvolume_2_focus_xh;
MCNUM mccvolume_2_focus_r;
MCNUM mccvolume_2_p_interact;
char mccvolume_2_mask_string[16384];
char mccvolume_2_mask_setting[16384];
MCNUM mccvolume_2_number_of_activations;

/* Setting parameters for component 'mask_1_1' [8]. */
char mccmask_1_1_material_string[16384];
MCNUM mccmask_1_1_priority;
MCNUM mccmask_1_1_xwidth;
MCNUM mccmask_1_1_yheight;
MCNUM mccmask_1_1_zdepth;
MCNUM mccmask_1_1_xwidth2;
MCNUM mccmask_1_1_yheight2;
MCNUM mccmask_1_1_visualize;
int mccmask_1_1_target_index;
MCNUM mccmask_1_1_target_x;
MCNUM mccmask_1_1_target_y;
MCNUM mccmask_1_1_target_z;
MCNUM mccmask_1_1_focus_aw;
MCNUM mccmask_1_1_focus_ah;
MCNUM mccmask_1_1_focus_xw;
MCNUM mccmask_1_1_focus_xh;
MCNUM mccmask_1_1_focus_r;
MCNUM mccmask_1_1_p_interact;
char mccmask_1_1_mask_string[16384];
char mccmask_1_1_mask_setting[16384];
MCNUM mccmask_1_1_number_of_activations;

/* Setting parameters for component 'mask_1_2' [9]. */
char mccmask_1_2_material_string[16384];
MCNUM mccmask_1_2_priority;
MCNUM mccmask_1_2_xwidth;
MCNUM mccmask_1_2_yheight;
MCNUM mccmask_1_2_zdepth;
MCNUM mccmask_1_2_xwidth2;
MCNUM mccmask_1_2_yheight2;
MCNUM mccmask_1_2_visualize;
int mccmask_1_2_target_index;
MCNUM mccmask_1_2_target_x;
MCNUM mccmask_1_2_target_y;
MCNUM mccmask_1_2_target_z;
MCNUM mccmask_1_2_focus_aw;
MCNUM mccmask_1_2_focus_ah;
MCNUM mccmask_1_2_focus_xw;
MCNUM mccmask_1_2_focus_xh;
MCNUM mccmask_1_2_focus_r;
MCNUM mccmask_1_2_p_interact;
char mccmask_1_2_mask_string[16384];
char mccmask_1_2_mask_setting[16384];
MCNUM mccmask_1_2_number_of_activations;

/* Setting parameters for component 'mask_2_1' [10]. */
char mccmask_2_1_material_string[16384];
MCNUM mccmask_2_1_priority;
MCNUM mccmask_2_1_xwidth;
MCNUM mccmask_2_1_yheight;
MCNUM mccmask_2_1_zdepth;
MCNUM mccmask_2_1_xwidth2;
MCNUM mccmask_2_1_yheight2;
MCNUM mccmask_2_1_visualize;
int mccmask_2_1_target_index;
MCNUM mccmask_2_1_target_x;
MCNUM mccmask_2_1_target_y;
MCNUM mccmask_2_1_target_z;
MCNUM mccmask_2_1_focus_aw;
MCNUM mccmask_2_1_focus_ah;
MCNUM mccmask_2_1_focus_xw;
MCNUM mccmask_2_1_focus_xh;
MCNUM mccmask_2_1_focus_r;
MCNUM mccmask_2_1_p_interact;
char mccmask_2_1_mask_string[16384];
char mccmask_2_1_mask_setting[16384];
MCNUM mccmask_2_1_number_of_activations;

/* Setting parameters for component 'mask_2_2' [11]. */
char mccmask_2_2_material_string[16384];
MCNUM mccmask_2_2_priority;
MCNUM mccmask_2_2_xwidth;
MCNUM mccmask_2_2_yheight;
MCNUM mccmask_2_2_zdepth;
MCNUM mccmask_2_2_xwidth2;
MCNUM mccmask_2_2_yheight2;
MCNUM mccmask_2_2_visualize;
int mccmask_2_2_target_index;
MCNUM mccmask_2_2_target_x;
MCNUM mccmask_2_2_target_y;
MCNUM mccmask_2_2_target_z;
MCNUM mccmask_2_2_focus_aw;
MCNUM mccmask_2_2_focus_ah;
MCNUM mccmask_2_2_focus_xw;
MCNUM mccmask_2_2_focus_xh;
MCNUM mccmask_2_2_focus_r;
MCNUM mccmask_2_2_p_interact;
char mccmask_2_2_mask_string[16384];
char mccmask_2_2_mask_setting[16384];
MCNUM mccmask_2_2_number_of_activations;

/* Setting parameters for component 'test_sample' [12]. */
MCNUM mcctest_sample_allow_inside_start;
MCNUM mcctest_sample_history_limit;
MCNUM mcctest_sample_enable_conditionals;
MCNUM mcctest_sample_inherit_number_of_scattering_events;

/* Definition parameters for component 'm4pi' [13]. */
#define mccm4pi_nx 180
#define mccm4pi_ny 180
/* Setting parameters for component 'm4pi' [13]. */
char mccm4pi_filename[16384];
MCNUM mccm4pi_radius;
MCNUM mccm4pi_restore_neutron;
int mccm4pi_nowritefile;

/* Definition parameters for component 'm4pi_1' [14]. */
#define mccm4pi_1_nx 180
#define mccm4pi_1_ny 180
/* Setting parameters for component 'm4pi_1' [14]. */
char mccm4pi_1_filename[16384];
MCNUM mccm4pi_1_radius;
MCNUM mccm4pi_1_restore_neutron;
int mccm4pi_1_nowritefile;

/* Definition parameters for component 'm4pi_2' [15]. */
#define mccm4pi_2_nx 180
#define mccm4pi_2_ny 180
/* Setting parameters for component 'm4pi_2' [15]. */
char mccm4pi_2_filename[16384];
MCNUM mccm4pi_2_radius;
MCNUM mccm4pi_2_restore_neutron;
int mccm4pi_2_nowritefile;

/* Definition parameters for component 'm4pi_both' [17]. */
#define mccm4pi_both_nx 180
#define mccm4pi_both_ny 180
/* Setting parameters for component 'm4pi_both' [17]. */
char mccm4pi_both_filename[16384];
MCNUM mccm4pi_both_radius;
MCNUM mccm4pi_both_restore_neutron;
int mccm4pi_both_nowritefile;

/* User component declarations. */

/* User declarations for component 'test_incoherent' [1]. */
#define mccompcurname  test_incoherent
#define mccompcurtype  Incoherent_process
#define mccompcurindex 1
#define This_process mcctest_incoherent_This_process
#define Incoherent_storage mcctest_incoherent_Incoherent_storage
#define effective_my_scattering mcctest_incoherent_effective_my_scattering
#define sigma mcctest_incoherent_sigma
#define f_QE mcctest_incoherent_f_QE
#define gamma mcctest_incoherent_gamma
#define packing_factor mcctest_incoherent_packing_factor
#define unit_cell_volume mcctest_incoherent_unit_cell_volume
#define interact_fraction mcctest_incoherent_interact_fraction
#line 123 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../contrib/union/Incoherent_process.comp"
// Needed for transport to the main component
struct global_process_element_struct global_process_element;
struct scattering_process_struct This_process;

#ifndef PROCESS_DETECTOR
	//struct pointer_to_global_process_list global_process_list = {0,NULL};
	#define PROCESS_DETECTOR dummy
#endif

// Declare for this component, to do calculations on the input / store in the transported data
struct Incoherent_physics_storage_struct Incoherent_storage;
double effective_my_scattering;

#line 5948 "./Union_test_mask.c"
#undef interact_fraction
#undef unit_cell_volume
#undef packing_factor
#undef gamma
#undef f_QE
#undef sigma
#undef effective_my_scattering
#undef Incoherent_storage
#undef This_process
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'test_material' [2]. */
#define mccompcurname  test_material
#define mccompcurtype  Union_make_material
#define mccompcurindex 2
#define loop_index mcctest_material_loop_index
#define this_material mcctest_material_this_material
#define accepted_processes mcctest_material_accepted_processes
#define global_material_element mcctest_material_global_material_element
#define process_string mcctest_material_process_string
#define my_absorption mcctest_material_my_absorption
#define absorber mcctest_material_absorber
#line 167 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../contrib/union/Union_make_material.comp"
// Needed for transport to the main component
struct global_material_element_struct global_material_element;
struct physics_struct this_material;

#ifndef MATERIAL_DETECTOR
	//struct pointer_to_global_material_list global_material_list = {0,NULL};
	#define MATERIAL_DETECTOR dummy
#endif

int loop_index;
int found_process;
int specified_processes;
char local_string[256];

struct pointer_to_1d_int_list accepted_processes = {0,NULL};

// Add setup for loggers since make_material is called before any volume / master
#ifndef UNION_LOGGER_DECLARE
	//struct pointer_to_global_logger_list global_all_volume_logger_list = {0,NULL};
    //struct pointer_to_global_logger_list global_specific_volumes_logger_list = {0,NULL};
	#define UNION_LOGGER_DECLARE dummy
#endif

#line 5997 "./Union_test_mask.c"
#undef absorber
#undef my_absorption
#undef process_string
#undef global_material_element
#undef accepted_processes
#undef this_material
#undef loop_index
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'a1' [3]. */
#define mccompcurname  a1
#define mccompcurtype  Progress_bar
#define mccompcurindex 3
#define IntermediateCnts mcca1_IntermediateCnts
#define StartTime mcca1_StartTime
#define EndTime mcca1_EndTime
#define CurrentTime mcca1_CurrentTime
#define profile mcca1_profile
#define percent mcca1_percent
#define flag_save mcca1_flag_save
#define minutes mcca1_minutes
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
#line 6032 "./Union_test_mask.c"
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

/* User declarations for component 'source' [4]. */
#define mccompcurname  source
#define mccompcurtype  Source_div
#define mccompcurindex 4
#define thetah mccsource_thetah
#define thetav mccsource_thetav
#define sigmah mccsource_sigmah
#define sigmav mccsource_sigmav
#define tan_h mccsource_tan_h
#define tan_v mccsource_tan_v
#define p_init mccsource_p_init
#define dist mccsource_dist
#define focus_xw mccsource_focus_xw
#define focus_yh mccsource_focus_yh
#define xwidth mccsource_xwidth
#define yheight mccsource_yheight
#define focus_aw mccsource_focus_aw
#define focus_ah mccsource_focus_ah
#define E0 mccsource_E0
#define dE mccsource_dE
#define lambda0 mccsource_lambda0
#define dlambda mccsource_dlambda
#define gauss mccsource_gauss
#define flux mccsource_flux
#line 69 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../sources/Source_div.comp"
double thetah, thetav, sigmah, sigmav, tan_h, tan_v, p_init, dist, focus_xw, focus_yh;
#line 6071 "./Union_test_mask.c"
#undef flux
#undef gauss
#undef dlambda
#undef lambda0
#undef dE
#undef E0
#undef focus_ah
#undef focus_aw
#undef yheight
#undef xwidth
#undef focus_yh
#undef focus_xw
#undef dist
#undef p_init
#undef tan_v
#undef tan_h
#undef sigmav
#undef sigmah
#undef thetav
#undef thetah
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'beam_center' [5]. */
#define mccompcurname  beam_center
#define mccompcurtype  Arm
#define mccompcurindex 5
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'volume_1' [6]. */
#define mccompcurname  volume_1
#define mccompcurtype  Union_box
#define mccompcurindex 6
#define loop_index mccvolume_1_loop_index
#define this_box_volume mccvolume_1_this_box_volume
#define global_geometry_element mccvolume_1_global_geometry_element
#define this_box_storage mccvolume_1_this_box_storage
#define material_string mccvolume_1_material_string
#define priority mccvolume_1_priority
#define xwidth mccvolume_1_xwidth
#define yheight mccvolume_1_yheight
#define zdepth mccvolume_1_zdepth
#define xwidth2 mccvolume_1_xwidth2
#define yheight2 mccvolume_1_yheight2
#define visualize mccvolume_1_visualize
#define target_index mccvolume_1_target_index
#define target_x mccvolume_1_target_x
#define target_y mccvolume_1_target_y
#define target_z mccvolume_1_target_z
#define focus_aw mccvolume_1_focus_aw
#define focus_ah mccvolume_1_focus_ah
#define focus_xw mccvolume_1_focus_xw
#define focus_xh mccvolume_1_focus_xh
#define focus_r mccvolume_1_focus_r
#define p_interact mccvolume_1_p_interact
#define mask_string mccvolume_1_mask_string
#define mask_setting mccvolume_1_mask_setting
#define number_of_activations mccvolume_1_number_of_activations
#line 203 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../contrib/union/Union_box.comp"
// Needed for transport to the main component
struct global_geometry_element_struct global_geometry_element;

#ifndef ANY_GEOMETRY_DETECTOR_DECLARE
    #define ANY_GEOMETRY_DETECTOR_DECLARE dummy
	//struct pointer_to_global_geometry_list global_geometry_list = {0,NULL};
#endif

int loop_index,found_geometries;

double x_component;
double y_component;
double z_component;

struct Volume_struct this_box_volume;
struct box_storage this_box_storage;

#line 6151 "./Union_test_mask.c"
#undef number_of_activations
#undef mask_setting
#undef mask_string
#undef p_interact
#undef focus_r
#undef focus_xh
#undef focus_xw
#undef focus_ah
#undef focus_aw
#undef target_z
#undef target_y
#undef target_x
#undef target_index
#undef visualize
#undef yheight2
#undef xwidth2
#undef zdepth
#undef yheight
#undef xwidth
#undef priority
#undef material_string
#undef this_box_storage
#undef global_geometry_element
#undef this_box_volume
#undef loop_index
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'volume_2' [7]. */
#define mccompcurname  volume_2
#define mccompcurtype  Union_box
#define mccompcurindex 7
#define loop_index mccvolume_2_loop_index
#define this_box_volume mccvolume_2_this_box_volume
#define global_geometry_element mccvolume_2_global_geometry_element
#define this_box_storage mccvolume_2_this_box_storage
#define material_string mccvolume_2_material_string
#define priority mccvolume_2_priority
#define xwidth mccvolume_2_xwidth
#define yheight mccvolume_2_yheight
#define zdepth mccvolume_2_zdepth
#define xwidth2 mccvolume_2_xwidth2
#define yheight2 mccvolume_2_yheight2
#define visualize mccvolume_2_visualize
#define target_index mccvolume_2_target_index
#define target_x mccvolume_2_target_x
#define target_y mccvolume_2_target_y
#define target_z mccvolume_2_target_z
#define focus_aw mccvolume_2_focus_aw
#define focus_ah mccvolume_2_focus_ah
#define focus_xw mccvolume_2_focus_xw
#define focus_xh mccvolume_2_focus_xh
#define focus_r mccvolume_2_focus_r
#define p_interact mccvolume_2_p_interact
#define mask_string mccvolume_2_mask_string
#define mask_setting mccvolume_2_mask_setting
#define number_of_activations mccvolume_2_number_of_activations
#line 203 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../contrib/union/Union_box.comp"
// Needed for transport to the main component
struct global_geometry_element_struct global_geometry_element;

#ifndef ANY_GEOMETRY_DETECTOR_DECLARE
    #define ANY_GEOMETRY_DETECTOR_DECLARE dummy
	//struct pointer_to_global_geometry_list global_geometry_list = {0,NULL};
#endif

int loop_index,found_geometries;

double x_component;
double y_component;
double z_component;

struct Volume_struct this_box_volume;
struct box_storage this_box_storage;

#line 6228 "./Union_test_mask.c"
#undef number_of_activations
#undef mask_setting
#undef mask_string
#undef p_interact
#undef focus_r
#undef focus_xh
#undef focus_xw
#undef focus_ah
#undef focus_aw
#undef target_z
#undef target_y
#undef target_x
#undef target_index
#undef visualize
#undef yheight2
#undef xwidth2
#undef zdepth
#undef yheight
#undef xwidth
#undef priority
#undef material_string
#undef this_box_storage
#undef global_geometry_element
#undef this_box_volume
#undef loop_index
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'mask_1_1' [8]. */
#define mccompcurname  mask_1_1
#define mccompcurtype  Union_box
#define mccompcurindex 8
#define loop_index mccmask_1_1_loop_index
#define this_box_volume mccmask_1_1_this_box_volume
#define global_geometry_element mccmask_1_1_global_geometry_element
#define this_box_storage mccmask_1_1_this_box_storage
#define material_string mccmask_1_1_material_string
#define priority mccmask_1_1_priority
#define xwidth mccmask_1_1_xwidth
#define yheight mccmask_1_1_yheight
#define zdepth mccmask_1_1_zdepth
#define xwidth2 mccmask_1_1_xwidth2
#define yheight2 mccmask_1_1_yheight2
#define visualize mccmask_1_1_visualize
#define target_index mccmask_1_1_target_index
#define target_x mccmask_1_1_target_x
#define target_y mccmask_1_1_target_y
#define target_z mccmask_1_1_target_z
#define focus_aw mccmask_1_1_focus_aw
#define focus_ah mccmask_1_1_focus_ah
#define focus_xw mccmask_1_1_focus_xw
#define focus_xh mccmask_1_1_focus_xh
#define focus_r mccmask_1_1_focus_r
#define p_interact mccmask_1_1_p_interact
#define mask_string mccmask_1_1_mask_string
#define mask_setting mccmask_1_1_mask_setting
#define number_of_activations mccmask_1_1_number_of_activations
#line 203 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../contrib/union/Union_box.comp"
// Needed for transport to the main component
struct global_geometry_element_struct global_geometry_element;

#ifndef ANY_GEOMETRY_DETECTOR_DECLARE
    #define ANY_GEOMETRY_DETECTOR_DECLARE dummy
	//struct pointer_to_global_geometry_list global_geometry_list = {0,NULL};
#endif

int loop_index,found_geometries;

double x_component;
double y_component;
double z_component;

struct Volume_struct this_box_volume;
struct box_storage this_box_storage;

#line 6305 "./Union_test_mask.c"
#undef number_of_activations
#undef mask_setting
#undef mask_string
#undef p_interact
#undef focus_r
#undef focus_xh
#undef focus_xw
#undef focus_ah
#undef focus_aw
#undef target_z
#undef target_y
#undef target_x
#undef target_index
#undef visualize
#undef yheight2
#undef xwidth2
#undef zdepth
#undef yheight
#undef xwidth
#undef priority
#undef material_string
#undef this_box_storage
#undef global_geometry_element
#undef this_box_volume
#undef loop_index
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'mask_1_2' [9]. */
#define mccompcurname  mask_1_2
#define mccompcurtype  Union_box
#define mccompcurindex 9
#define loop_index mccmask_1_2_loop_index
#define this_box_volume mccmask_1_2_this_box_volume
#define global_geometry_element mccmask_1_2_global_geometry_element
#define this_box_storage mccmask_1_2_this_box_storage
#define material_string mccmask_1_2_material_string
#define priority mccmask_1_2_priority
#define xwidth mccmask_1_2_xwidth
#define yheight mccmask_1_2_yheight
#define zdepth mccmask_1_2_zdepth
#define xwidth2 mccmask_1_2_xwidth2
#define yheight2 mccmask_1_2_yheight2
#define visualize mccmask_1_2_visualize
#define target_index mccmask_1_2_target_index
#define target_x mccmask_1_2_target_x
#define target_y mccmask_1_2_target_y
#define target_z mccmask_1_2_target_z
#define focus_aw mccmask_1_2_focus_aw
#define focus_ah mccmask_1_2_focus_ah
#define focus_xw mccmask_1_2_focus_xw
#define focus_xh mccmask_1_2_focus_xh
#define focus_r mccmask_1_2_focus_r
#define p_interact mccmask_1_2_p_interact
#define mask_string mccmask_1_2_mask_string
#define mask_setting mccmask_1_2_mask_setting
#define number_of_activations mccmask_1_2_number_of_activations
#line 203 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../contrib/union/Union_box.comp"
// Needed for transport to the main component
struct global_geometry_element_struct global_geometry_element;

#ifndef ANY_GEOMETRY_DETECTOR_DECLARE
    #define ANY_GEOMETRY_DETECTOR_DECLARE dummy
	//struct pointer_to_global_geometry_list global_geometry_list = {0,NULL};
#endif

int loop_index,found_geometries;

double x_component;
double y_component;
double z_component;

struct Volume_struct this_box_volume;
struct box_storage this_box_storage;

#line 6382 "./Union_test_mask.c"
#undef number_of_activations
#undef mask_setting
#undef mask_string
#undef p_interact
#undef focus_r
#undef focus_xh
#undef focus_xw
#undef focus_ah
#undef focus_aw
#undef target_z
#undef target_y
#undef target_x
#undef target_index
#undef visualize
#undef yheight2
#undef xwidth2
#undef zdepth
#undef yheight
#undef xwidth
#undef priority
#undef material_string
#undef this_box_storage
#undef global_geometry_element
#undef this_box_volume
#undef loop_index
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'mask_2_1' [10]. */
#define mccompcurname  mask_2_1
#define mccompcurtype  Union_box
#define mccompcurindex 10
#define loop_index mccmask_2_1_loop_index
#define this_box_volume mccmask_2_1_this_box_volume
#define global_geometry_element mccmask_2_1_global_geometry_element
#define this_box_storage mccmask_2_1_this_box_storage
#define material_string mccmask_2_1_material_string
#define priority mccmask_2_1_priority
#define xwidth mccmask_2_1_xwidth
#define yheight mccmask_2_1_yheight
#define zdepth mccmask_2_1_zdepth
#define xwidth2 mccmask_2_1_xwidth2
#define yheight2 mccmask_2_1_yheight2
#define visualize mccmask_2_1_visualize
#define target_index mccmask_2_1_target_index
#define target_x mccmask_2_1_target_x
#define target_y mccmask_2_1_target_y
#define target_z mccmask_2_1_target_z
#define focus_aw mccmask_2_1_focus_aw
#define focus_ah mccmask_2_1_focus_ah
#define focus_xw mccmask_2_1_focus_xw
#define focus_xh mccmask_2_1_focus_xh
#define focus_r mccmask_2_1_focus_r
#define p_interact mccmask_2_1_p_interact
#define mask_string mccmask_2_1_mask_string
#define mask_setting mccmask_2_1_mask_setting
#define number_of_activations mccmask_2_1_number_of_activations
#line 203 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../contrib/union/Union_box.comp"
// Needed for transport to the main component
struct global_geometry_element_struct global_geometry_element;

#ifndef ANY_GEOMETRY_DETECTOR_DECLARE
    #define ANY_GEOMETRY_DETECTOR_DECLARE dummy
	//struct pointer_to_global_geometry_list global_geometry_list = {0,NULL};
#endif

int loop_index,found_geometries;

double x_component;
double y_component;
double z_component;

struct Volume_struct this_box_volume;
struct box_storage this_box_storage;

#line 6459 "./Union_test_mask.c"
#undef number_of_activations
#undef mask_setting
#undef mask_string
#undef p_interact
#undef focus_r
#undef focus_xh
#undef focus_xw
#undef focus_ah
#undef focus_aw
#undef target_z
#undef target_y
#undef target_x
#undef target_index
#undef visualize
#undef yheight2
#undef xwidth2
#undef zdepth
#undef yheight
#undef xwidth
#undef priority
#undef material_string
#undef this_box_storage
#undef global_geometry_element
#undef this_box_volume
#undef loop_index
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'mask_2_2' [11]. */
#define mccompcurname  mask_2_2
#define mccompcurtype  Union_box
#define mccompcurindex 11
#define loop_index mccmask_2_2_loop_index
#define this_box_volume mccmask_2_2_this_box_volume
#define global_geometry_element mccmask_2_2_global_geometry_element
#define this_box_storage mccmask_2_2_this_box_storage
#define material_string mccmask_2_2_material_string
#define priority mccmask_2_2_priority
#define xwidth mccmask_2_2_xwidth
#define yheight mccmask_2_2_yheight
#define zdepth mccmask_2_2_zdepth
#define xwidth2 mccmask_2_2_xwidth2
#define yheight2 mccmask_2_2_yheight2
#define visualize mccmask_2_2_visualize
#define target_index mccmask_2_2_target_index
#define target_x mccmask_2_2_target_x
#define target_y mccmask_2_2_target_y
#define target_z mccmask_2_2_target_z
#define focus_aw mccmask_2_2_focus_aw
#define focus_ah mccmask_2_2_focus_ah
#define focus_xw mccmask_2_2_focus_xw
#define focus_xh mccmask_2_2_focus_xh
#define focus_r mccmask_2_2_focus_r
#define p_interact mccmask_2_2_p_interact
#define mask_string mccmask_2_2_mask_string
#define mask_setting mccmask_2_2_mask_setting
#define number_of_activations mccmask_2_2_number_of_activations
#line 203 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../contrib/union/Union_box.comp"
// Needed for transport to the main component
struct global_geometry_element_struct global_geometry_element;

#ifndef ANY_GEOMETRY_DETECTOR_DECLARE
    #define ANY_GEOMETRY_DETECTOR_DECLARE dummy
	//struct pointer_to_global_geometry_list global_geometry_list = {0,NULL};
#endif

int loop_index,found_geometries;

double x_component;
double y_component;
double z_component;

struct Volume_struct this_box_volume;
struct box_storage this_box_storage;

#line 6536 "./Union_test_mask.c"
#undef number_of_activations
#undef mask_setting
#undef mask_string
#undef p_interact
#undef focus_r
#undef focus_xh
#undef focus_xw
#undef focus_ah
#undef focus_aw
#undef target_z
#undef target_y
#undef target_x
#undef target_index
#undef visualize
#undef yheight2
#undef xwidth2
#undef zdepth
#undef yheight
#undef xwidth
#undef priority
#undef material_string
#undef this_box_storage
#undef global_geometry_element
#undef this_box_volume
#undef loop_index
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'test_sample' [12]. */
#define mccompcurname  test_sample
#define mccompcurtype  Union_master
#define mccompcurindex 12
#define verbal mcctest_sample_verbal
#define list_verbal mcctest_sample_list_verbal
#define trace_verbal mcctest_sample_trace_verbal
#define finally_verbal mcctest_sample_finally_verbal
#define starting_volume_warning mcctest_sample_starting_volume_warning
#define global_master_element mcctest_sample_global_master_element
#define this_global_master_index mcctest_sample_this_global_master_index
#define previous_master_index mcctest_sample_previous_master_index
#define geometry_list_index mcctest_sample_geometry_list_index
#define intersection_time_table mcctest_sample_intersection_time_table
#define Volumes mcctest_sample_Volumes
#define Geometries mcctest_sample_Geometries
#define starting_lists mcctest_sample_starting_lists
#define r mcctest_sample_r
#define r_start mcctest_sample_r_start
#define v mcctest_sample_v
#define error_msg mcctest_sample_error_msg
#define component_error_msg mcctest_sample_component_error_msg
#define string_output mcctest_sample_string_output
#define number_of_volumes mcctest_sample_number_of_volumes
#define volume_index mcctest_sample_volume_index
#define process_index mcctest_sample_process_index
#define solutions mcctest_sample_solutions
#define max_number_of_processes mcctest_sample_max_number_of_processes
#define limit mcctest_sample_limit
#define solution mcctest_sample_solution
#define min_solution mcctest_sample_min_solution
#define min_volume mcctest_sample_min_volume
#define time_found mcctest_sample_time_found
#define intersection_time mcctest_sample_intersection_time
#define min_intersection_time mcctest_sample_min_intersection_time
#define process mcctest_sample_process
#define process_start mcctest_sample_process_start
#define my_trace mcctest_sample_my_trace
#define p_my_trace mcctest_sample_p_my_trace
#define my_trace_fraction_control mcctest_sample_my_trace_fraction_control
#define k mcctest_sample_k
#define k_new mcctest_sample_k_new
#define k_old mcctest_sample_k_old
#define v_length mcctest_sample_v_length
#define my_sum mcctest_sample_my_sum
#define my_sum_plus_abs mcctest_sample_my_sum_plus_abs
#define culmative_probability mcctest_sample_culmative_probability
#define mc_prop mcctest_sample_mc_prop
#define time_to_scattering mcctest_sample_time_to_scattering
#define length_to_scattering mcctest_sample_length_to_scattering
#define length_to_boundery mcctest_sample_length_to_boundery
#define time_to_boundery mcctest_sample_time_to_boundery
#define selected_process mcctest_sample_selected_process
#define scattering_event mcctest_sample_scattering_event
#define time_propagated_without_scattering mcctest_sample_time_propagated_without_scattering
#define a_next_volume_found mcctest_sample_a_next_volume_found
#define next_volume mcctest_sample_next_volume
#define next_volume_priority mcctest_sample_next_volume_priority
#define done mcctest_sample_done
#define current_volume mcctest_sample_current_volume
#define number_of_solutions mcctest_sample_number_of_solutions
#define number_of_solutions_static mcctest_sample_number_of_solutions_static
#define check mcctest_sample_check
#define start mcctest_sample_start
#define intersection_with_children mcctest_sample_intersection_with_children
#define geometry_output mcctest_sample_geometry_output
#define tree_next_volume mcctest_sample_tree_next_volume
#define pre_allocated1 mcctest_sample_pre_allocated1
#define pre_allocated2 mcctest_sample_pre_allocated2
#define pre_allocated3 mcctest_sample_pre_allocated3
#define ray_position mcctest_sample_ray_position
#define ray_velocity mcctest_sample_ray_velocity
#define ray_velocity_final mcctest_sample_ray_velocity_final
#define volume_0_found mcctest_sample_volume_0_found
#define scattered_flag mcctest_sample_scattered_flag
#define scattered_flag_VP mcctest_sample_scattered_flag_VP
#define master_transposed_rotation_matrix mcctest_sample_master_transposed_rotation_matrix
#define temp_rotation_matrix mcctest_sample_temp_rotation_matrix
#define non_rotated_position mcctest_sample_non_rotated_position
#define rotated_position mcctest_sample_rotated_position
#define enable_tagging mcctest_sample_enable_tagging
#define stop_tagging_ray mcctest_sample_stop_tagging_ray
#define stop_creating_nodes mcctest_sample_stop_creating_nodes
#define enable_tagging_check mcctest_sample_enable_tagging_check
#define master_tagging_node_list mcctest_sample_master_tagging_node_list
#define current_tagging_node mcctest_sample_current_tagging_node
#define tagging_leaf_counter mcctest_sample_tagging_leaf_counter
#define number_of_scattering_events mcctest_sample_number_of_scattering_events
#define real_transmission_probability mcctest_sample_real_transmission_probability
#define mc_transmission_probability mcctest_sample_mc_transmission_probability
#define number_of_masks mcctest_sample_number_of_masks
#define number_of_masked_volumes mcctest_sample_number_of_masked_volumes
#define need_to_run_within_which_volume mcctest_sample_need_to_run_within_which_volume
#define mask_index_main mcctest_sample_mask_index_main
#define mask_iterate mcctest_sample_mask_iterate
#define mask_status_list mcctest_sample_mask_status_list
#define current_mask_intersect_list_status mcctest_sample_current_mask_intersect_list_status
#define mask_volume_index_list mcctest_sample_mask_volume_index_list
#define geometry_component_index_list mcctest_sample_geometry_component_index_list
#define Volume_copies_allocated mcctest_sample_Volume_copies_allocated
#define p_old mcctest_sample_p_old
#define this_logger mcctest_sample_this_logger
#define conditional_status mcctest_sample_conditional_status
#define tagging_conditional_list mcctest_sample_tagging_conditional_list
#define free_tagging_conditioanl_list mcctest_sample_free_tagging_conditioanl_list
#define logger_conditional_extend_array mcctest_sample_logger_conditional_extend_array
#define tagging_conditional_extend mcctest_sample_tagging_conditional_extend
#define max_conditional_extend_index mcctest_sample_max_conditional_extend_index
#define safty_distance mcctest_sample_safty_distance
#define safty_distance2 mcctest_sample_safty_distance2
#define number_of_processes_array mcctest_sample_number_of_processes_array
#define temporary_focus_data mcctest_sample_temporary_focus_data
#define focus_data_index mcctest_sample_focus_data_index
#define allow_inside_start mcctest_sample_allow_inside_start
#define history_limit mcctest_sample_history_limit
#define enable_conditionals mcctest_sample_enable_conditionals
#define inherit_number_of_scattering_events mcctest_sample_inherit_number_of_scattering_events
#line 68 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../contrib/union/Union_master.comp"
  // Settings that can be moved to input.
  int verbal = 1;
  int list_verbal = 0;
  int finally_verbal = 0;
  
  // Performance intensive, obsolete, now done in precompiler
  int trace_verbal = 0;
  int enable_tagging = 0;

  // It is possible to surpress warnings on starting volume by setting this to 1
  int starting_volume_warning = 0;
  
  // New precompiler settings for verbal / tagging
  //#define Union_trace_verbal_setting
  //#define Union_enable_tagging_setting

  // Declare the global variables (not to be in output parameters)
  struct global_master_element_struct global_master_element;
  int this_global_master_index;
  
  #ifndef MASTER_DETECTOR
	//struct pointer_to_global_master_list global_master_list = {0,NULL};
	#define MASTER_DETECTOR dummy
  #endif
  
  // variables used for assigning global information to local variables
  int previous_master_index,geometry_list_index;

  // The main structures used in this component
  struct intersection_time_table_struct intersection_time_table;
  struct Volume_struct **Volumes;
  struct geometry_struct **Geometries;
  struct Volume_struct **Volume_copies;
  struct starting_lists_struct starting_lists;
  
  // garbage collection for volume_copies
  struct pointer_to_1d_int_list Volume_copies_allocated;

  // Vectors in old format (still used by intersect function, will go to Coords in future)
  double r[3],r_start[3],v[3];

  // Error handling
  int error_msg, component_error_msg=0;
  
  // For verbal output
  char string_output[128];
  
  // Variables for ray-tracing algorithm
  int number_of_volumes, volume_index, process_index;
  int iterate,solutions,max_number_of_processes,limit;
  int solution,min_solution,min_volume,time_found;
  double intersection_time,min_intersection_time;
  
  struct scattering_process_struct *process,*process_start;
  double *my_trace,*p_my_trace,*my_trace_fraction_control;
  double k[3],k_new[3],k_old[3],k_rotated[3];
  double v_length,my_sum,my_sum_plus_abs,culmative_probability,mc_prop,time_to_scattering;
  double length_to_scattering,length_to_boundery,length_to_boundery_fp,time_to_boundery;
  int selected_process,scattering_event;
  double time_propagated_without_scattering;
  
  int a_next_volume_found,next_volume;
  double next_volume_priority;
  
  int done,current_volume,ray_sucseeded;
  int *number_of_solutions;
  int number_of_solutions_static;
  int *check,*start;
  int intersection_with_children,geometry_output;
  
  // For within_which_volume
  int tree_next_volume;
  int *pre_allocated1,*pre_allocated2,*pre_allocated3;
  Coords ray_position,ray_velocity,ray_velocity_rotated,ray_velocity_final,wavevector,wavevector_rotated;
  int volume_0_found=0;
  
  int *scattered_flag;
  int **scattered_flag_VP;
  
  // For coordinate transformations
  Rotation master_transposed_rotation_matrix;
  Rotation temp_rotation_matrix;
  Rotation temp_transpose_rotation_matrix;
  Coords non_rotated_position;
  Coords rotated_position;
  int non_isotropic_found;
  
  // For tagging
  struct list_of_tagging_tree_node_pointers master_tagging_node_list;
  struct tagging_tree_node_struct *current_tagging_node;
  
  int tagging_leaf_counter=0,stop_tagging_ray,stop_creating_nodes;
  int number_of_scattering_events;
  
  // For geometry p interact
  double real_transmission_probability,mc_transmission_probability;
  
  // Process p interact
  int number_of_process_interacts_set,index_of_lacking_process;
  double total_process_interact;
  
  // Volume nr -> component index
  struct pointer_to_1d_int_list geometry_component_index_list;
  
  // Masks
  struct pointer_to_1d_int_list mask_volume_index_list;
  int number_of_masks=0;
  int number_of_masked_volumes=0;
  struct pointer_to_1d_int_list mask_status_list;
  struct pointer_to_1d_int_list current_mask_intersect_list_status;
  int mask_index_main,mask_iterate;
  int *mask_start,*mask_check;
  int need_to_run_within_which_volume;
  
  // Loggers
  //struct logger_with_data_struct loggers_with_data_array;
  int *number_of_processes_array;
  double p_old;
  int log_index,conditional_status;
  struct logger_struct *this_logger;
  //  union detector_pointer_union detector_pointer;
  
  // Conditionals
  struct conditional_list_struct *tagging_conditional_list;
  int *logger_conditional_extend_array;
  int max_conditional_extend_index;
  int tagging_conditional_extend;
  int free_tagging_conditioanl_list;
  
  // Reliability control
  // Safty distance is needed to avoid having ray positions closer to a wall than the precision of intersection functions
  double safty_distance,safty_distance2;
  
  // Focusing
  struct focus_data_struct temporary_focus_data;
  int focus_data_index;
  
#line 6821 "./Union_test_mask.c"
#undef inherit_number_of_scattering_events
#undef enable_conditionals
#undef history_limit
#undef allow_inside_start
#undef focus_data_index
#undef temporary_focus_data
#undef number_of_processes_array
#undef safty_distance2
#undef safty_distance
#undef max_conditional_extend_index
#undef tagging_conditional_extend
#undef logger_conditional_extend_array
#undef free_tagging_conditioanl_list
#undef tagging_conditional_list
#undef conditional_status
#undef this_logger
#undef p_old
#undef Volume_copies_allocated
#undef geometry_component_index_list
#undef mask_volume_index_list
#undef current_mask_intersect_list_status
#undef mask_status_list
#undef mask_iterate
#undef mask_index_main
#undef need_to_run_within_which_volume
#undef number_of_masked_volumes
#undef number_of_masks
#undef mc_transmission_probability
#undef real_transmission_probability
#undef number_of_scattering_events
#undef tagging_leaf_counter
#undef current_tagging_node
#undef master_tagging_node_list
#undef enable_tagging_check
#undef stop_creating_nodes
#undef stop_tagging_ray
#undef enable_tagging
#undef rotated_position
#undef non_rotated_position
#undef temp_rotation_matrix
#undef master_transposed_rotation_matrix
#undef scattered_flag_VP
#undef scattered_flag
#undef volume_0_found
#undef ray_velocity_final
#undef ray_velocity
#undef ray_position
#undef pre_allocated3
#undef pre_allocated2
#undef pre_allocated1
#undef tree_next_volume
#undef geometry_output
#undef intersection_with_children
#undef start
#undef check
#undef number_of_solutions_static
#undef number_of_solutions
#undef current_volume
#undef done
#undef next_volume_priority
#undef next_volume
#undef a_next_volume_found
#undef time_propagated_without_scattering
#undef scattering_event
#undef selected_process
#undef time_to_boundery
#undef length_to_boundery
#undef length_to_scattering
#undef time_to_scattering
#undef mc_prop
#undef culmative_probability
#undef my_sum_plus_abs
#undef my_sum
#undef v_length
#undef k_old
#undef k_new
#undef k
#undef my_trace_fraction_control
#undef p_my_trace
#undef my_trace
#undef process_start
#undef process
#undef min_intersection_time
#undef intersection_time
#undef time_found
#undef min_volume
#undef min_solution
#undef solution
#undef limit
#undef max_number_of_processes
#undef solutions
#undef process_index
#undef volume_index
#undef number_of_volumes
#undef string_output
#undef component_error_msg
#undef error_msg
#undef v
#undef r_start
#undef r
#undef starting_lists
#undef Geometries
#undef Volumes
#undef intersection_time_table
#undef geometry_list_index
#undef previous_master_index
#undef this_global_master_index
#undef global_master_element
#undef starting_volume_warning
#undef finally_verbal
#undef trace_verbal
#undef list_verbal
#undef verbal
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'm4pi' [13]. */
#define mccompcurname  m4pi
#define mccompcurtype  PSD_monitor_4PI
#define mccompcurindex 13
#define nx mccm4pi_nx
#define ny mccm4pi_ny
#define PSD_N mccm4pi_PSD_N
#define PSD_p mccm4pi_PSD_p
#define PSD_p2 mccm4pi_PSD_p2
#define filename mccm4pi_filename
#define radius mccm4pi_radius
#define restore_neutron mccm4pi_restore_neutron
#define nowritefile mccm4pi_nowritefile
#line 55 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor_4PI.comp"
double PSD_N[nx][ny];
double PSD_p[nx][ny];
double PSD_p2[nx][ny];
#line 6956 "./Union_test_mask.c"
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

/* User declarations for component 'm4pi_1' [14]. */
#define mccompcurname  m4pi_1
#define mccompcurtype  PSD_monitor_4PI
#define mccompcurindex 14
#define nx mccm4pi_1_nx
#define ny mccm4pi_1_ny
#define PSD_N mccm4pi_1_PSD_N
#define PSD_p mccm4pi_1_PSD_p
#define PSD_p2 mccm4pi_1_PSD_p2
#define filename mccm4pi_1_filename
#define radius mccm4pi_1_radius
#define restore_neutron mccm4pi_1_restore_neutron
#define nowritefile mccm4pi_1_nowritefile
#line 55 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor_4PI.comp"
double PSD_N[nx][ny];
double PSD_p[nx][ny];
double PSD_p2[nx][ny];
#line 6987 "./Union_test_mask.c"
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

/* User declarations for component 'm4pi_2' [15]. */
#define mccompcurname  m4pi_2
#define mccompcurtype  PSD_monitor_4PI
#define mccompcurindex 15
#define nx mccm4pi_2_nx
#define ny mccm4pi_2_ny
#define PSD_N mccm4pi_2_PSD_N
#define PSD_p mccm4pi_2_PSD_p
#define PSD_p2 mccm4pi_2_PSD_p2
#define filename mccm4pi_2_filename
#define radius mccm4pi_2_radius
#define restore_neutron mccm4pi_2_restore_neutron
#define nowritefile mccm4pi_2_nowritefile
#line 55 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor_4PI.comp"
double PSD_N[nx][ny];
double PSD_p[nx][ny];
double PSD_p2[nx][ny];
#line 7018 "./Union_test_mask.c"
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

/* User declarations for component 'arm_1' [16]. */
#define mccompcurname  arm_1
#define mccompcurtype  Arm
#define mccompcurindex 16
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'm4pi_both' [17]. */
#define mccompcurname  m4pi_both
#define mccompcurtype  PSD_monitor_4PI
#define mccompcurindex 17
#define nx mccm4pi_both_nx
#define ny mccm4pi_both_ny
#define PSD_N mccm4pi_both_PSD_N
#define PSD_p mccm4pi_both_PSD_p
#define PSD_p2 mccm4pi_both_PSD_p2
#define filename mccm4pi_both_filename
#define radius mccm4pi_both_radius
#define restore_neutron mccm4pi_both_restore_neutron
#define nowritefile mccm4pi_both_nowritefile
#line 55 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor_4PI.comp"
double PSD_N[nx][ny];
double PSD_p[nx][ny];
double PSD_p2[nx][ny];
#line 7057 "./Union_test_mask.c"
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

Coords mcposatest_incoherent, mcposrtest_incoherent;
Rotation mcrotatest_incoherent, mcrotrtest_incoherent;
Coords mcposatest_material, mcposrtest_material;
Rotation mcrotatest_material, mcrotrtest_material;
Coords mcposaa1, mcposra1;
Rotation mcrotaa1, mcrotra1;
Coords mcposasource, mcposrsource;
Rotation mcrotasource, mcrotrsource;
Coords mcposabeam_center, mcposrbeam_center;
Rotation mcrotabeam_center, mcrotrbeam_center;
Coords mcposavolume_1, mcposrvolume_1;
Rotation mcrotavolume_1, mcrotrvolume_1;
Coords mcposavolume_2, mcposrvolume_2;
Rotation mcrotavolume_2, mcrotrvolume_2;
Coords mcposamask_1_1, mcposrmask_1_1;
Rotation mcrotamask_1_1, mcrotrmask_1_1;
Coords mcposamask_1_2, mcposrmask_1_2;
Rotation mcrotamask_1_2, mcrotrmask_1_2;
Coords mcposamask_2_1, mcposrmask_2_1;
Rotation mcrotamask_2_1, mcrotrmask_2_1;
Coords mcposamask_2_2, mcposrmask_2_2;
Rotation mcrotamask_2_2, mcrotrmask_2_2;
Coords mcposatest_sample, mcposrtest_sample;
Rotation mcrotatest_sample, mcrotrtest_sample;
Coords mcposam4pi, mcposrm4pi;
Rotation mcrotam4pi, mcrotrm4pi;
Coords mcposam4pi_1, mcposrm4pi_1;
Rotation mcrotam4pi_1, mcrotrm4pi_1;
Coords mcposam4pi_2, mcposrm4pi_2;
Rotation mcrotam4pi_2, mcrotrm4pi_2;
Coords mcposaarm_1, mcposrarm_1;
Rotation mcrotaarm_1, mcrotrarm_1;
Coords mcposam4pi_both, mcposrm4pi_both;
Rotation mcrotam4pi_both, mcrotrm4pi_both;

MCNUM mcnx, mcny, mcnz, mcnvx, mcnvy, mcnvz, mcnt, mcnsx, mcnsy, mcnsz, mcnp;

/* end declare */

void mcinit(void) {
#define mccompcurname  Union_test_mask
#define mccompcurtype  INSTRUMENT
#define mccompcurindex 0
#define mcposaUnion_test_mask coords_set(0,0,0)
#undef mcposaUnion_test_mask
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
    /* Component test_incoherent. */
  /* Setting parameters for component test_incoherent. */
  SIG_MESSAGE("test_incoherent (Init:SetPar)");
#line 38 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mcctest_incoherent_sigma = 20;
#line 55 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mcctest_incoherent_f_QE = 0;
#line 55 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mcctest_incoherent_gamma = 0;
#line 38 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mcctest_incoherent_packing_factor = 1;
#line 38 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mcctest_incoherent_unit_cell_volume = 66.4;
#line 55 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mcctest_incoherent_interact_fraction = -1;
#line 7144 "./Union_test_mask.c"

  SIG_MESSAGE("test_incoherent (Init:Place/Rotate)");
  rot_set_rotation(mcrotatest_incoherent,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 7151 "./Union_test_mask.c"
  rot_copy(mcrotrtest_incoherent, mcrotatest_incoherent);
  mcposatest_incoherent = coords_set(
#line 39 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
    0,
#line 39 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
    0,
#line 39 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
    0);
#line 7160 "./Union_test_mask.c"
  mctc1 = coords_neg(mcposatest_incoherent);
  mcposrtest_incoherent = rot_apply(mcrotatest_incoherent, mctc1);
  mcDEBUG_COMPONENT("test_incoherent", mcposatest_incoherent, mcrotatest_incoherent)
  mccomp_posa[1] = mcposatest_incoherent;
  mccomp_posr[1] = mcposrtest_incoherent;
  mcNCounter[1]  = mcPCounter[1] = mcP2Counter[1] = 0;
  mcAbsorbProp[1]= 0;
    /* Component test_material. */
  /* Setting parameters for component test_material. */
  SIG_MESSAGE("test_material (Init:SetPar)");
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  if("NULL") strncpy(mcctest_material_process_string, "NULL" ? "NULL" : "", 16384); else mcctest_material_process_string[0]='\0';
#line 41 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mcctest_material_my_absorption = 4;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mcctest_material_absorber = 0;
#line 7177 "./Union_test_mask.c"

  SIG_MESSAGE("test_material (Init:Place/Rotate)");
  rot_set_rotation(mcrotatest_material,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 7184 "./Union_test_mask.c"
  rot_transpose(mcrotatest_incoherent, mctr1);
  rot_mul(mcrotatest_material, mctr1, mcrotrtest_material);
  mcposatest_material = coords_set(
#line 42 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
    0,
#line 42 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
    0,
#line 42 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
    0);
#line 7194 "./Union_test_mask.c"
  mctc1 = coords_sub(mcposatest_incoherent, mcposatest_material);
  mcposrtest_material = rot_apply(mcrotatest_material, mctc1);
  mcDEBUG_COMPONENT("test_material", mcposatest_material, mcrotatest_material)
  mccomp_posa[2] = mcposatest_material;
  mccomp_posr[2] = mcposrtest_material;
  mcNCounter[2]  = mcPCounter[2] = mcP2Counter[2] = 0;
  mcAbsorbProp[2]= 0;
    /* Component a1. */
  /* Setting parameters for component a1. */
  SIG_MESSAGE("a1 (Init:SetPar)");
#line 39 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  if("NULL") strncpy(mcca1_profile, "NULL" ? "NULL" : "", 16384); else mcca1_profile[0]='\0';
#line 39 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mcca1_percent = 10;
#line 39 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mcca1_flag_save = 0;
#line 39 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mcca1_minutes = 0;
#line 7213 "./Union_test_mask.c"

  SIG_MESSAGE("a1 (Init:Place/Rotate)");
  rot_set_rotation(mcrotaa1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 7220 "./Union_test_mask.c"
  rot_transpose(mcrotatest_material, mctr1);
  rot_mul(mcrotaa1, mctr1, mcrotra1);
  mcposaa1 = coords_set(
#line 45 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
    0,
#line 45 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
    0,
#line 45 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
    0);
#line 7230 "./Union_test_mask.c"
  mctc1 = coords_sub(mcposatest_material, mcposaa1);
  mcposra1 = rot_apply(mcrotaa1, mctc1);
  mcDEBUG_COMPONENT("a1", mcposaa1, mcrotaa1)
  mccomp_posa[3] = mcposaa1;
  mccomp_posr[3] = mcposra1;
  mcNCounter[3]  = mcPCounter[3] = mcP2Counter[3] = 0;
  mcAbsorbProp[3]= 0;
    /* Component source. */
  /* Setting parameters for component source. */
  SIG_MESSAGE("source (Init:SetPar)");
#line 49 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccsource_xwidth = 0.0005;
#line 49 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccsource_yheight = 0.0005;
#line 49 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccsource_focus_aw = 0.3;
#line 49 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccsource_focus_ah = 1;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccsource_E0 = 10;
#line 51 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccsource_dE = 0;
#line 64 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccsource_lambda0 = 0.0;
#line 64 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccsource_dlambda = 0.0;
#line 64 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccsource_gauss = 0;
#line 64 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccsource_flux = 1;
#line 7261 "./Union_test_mask.c"

  SIG_MESSAGE("source (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 52 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
    (0)*DEG2RAD,
#line 52 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
    (0)*DEG2RAD,
#line 52 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
    (0)*DEG2RAD);
#line 7271 "./Union_test_mask.c"
  rot_mul(mctr1, mcrotaa1, mcrotasource);
  rot_transpose(mcrotaa1, mctr1);
  rot_mul(mcrotasource, mctr1, mcrotrsource);
  mctc1 = coords_set(
#line 52 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
    0,
#line 52 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
    0,
#line 52 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
    0);
#line 7282 "./Union_test_mask.c"
  rot_transpose(mcrotaa1, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposasource = coords_add(mcposaa1, mctc2);
  mctc1 = coords_sub(mcposaa1, mcposasource);
  mcposrsource = rot_apply(mcrotasource, mctc1);
  mcDEBUG_COMPONENT("source", mcposasource, mcrotasource)
  mccomp_posa[4] = mcposasource;
  mccomp_posr[4] = mcposrsource;
  mcNCounter[4]  = mcPCounter[4] = mcP2Counter[4] = 0;
  mcAbsorbProp[4]= 0;
    /* Component beam_center. */
  /* Setting parameters for component beam_center. */
  SIG_MESSAGE("beam_center (Init:SetPar)");

  SIG_MESSAGE("beam_center (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 63 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
    (0)*DEG2RAD,
#line 63 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
    (0)*DEG2RAD,
#line 63 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
    (0)*DEG2RAD);
#line 7305 "./Union_test_mask.c"
  rot_mul(mctr1, mcrotaa1, mcrotabeam_center);
  rot_transpose(mcrotasource, mctr1);
  rot_mul(mcrotabeam_center, mctr1, mcrotrbeam_center);
  mctc1 = coords_set(
#line 62 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
    0,
#line 62 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
    0,
#line 62 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
    10.4);
#line 7316 "./Union_test_mask.c"
  rot_transpose(mcrotaa1, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposabeam_center = coords_add(mcposaa1, mctc2);
  mctc1 = coords_sub(mcposasource, mcposabeam_center);
  mcposrbeam_center = rot_apply(mcrotabeam_center, mctc1);
  mcDEBUG_COMPONENT("beam_center", mcposabeam_center, mcrotabeam_center)
  mccomp_posa[5] = mcposabeam_center;
  mccomp_posr[5] = mcposrbeam_center;
  mcNCounter[5]  = mcPCounter[5] = mcP2Counter[5] = 0;
  mcAbsorbProp[5]= 0;
    /* Component volume_1. */
  /* Setting parameters for component volume_1. */
  SIG_MESSAGE("volume_1 (Init:SetPar)");
#line 67 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  if("test_material") strncpy(mccvolume_1_material_string, "test_material" ? "test_material" : "", 16384); else mccvolume_1_material_string[0]='\0';
#line 67 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccvolume_1_priority = 10;
#line 67 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccvolume_1_xwidth = 0.1;
#line 67 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccvolume_1_yheight = 0.1;
#line 67 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccvolume_1_zdepth = 0.1;
#line 70 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccvolume_1_xwidth2 = -1;
#line 70 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccvolume_1_yheight2 = -1;
#line 70 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccvolume_1_visualize = 1;
#line 70 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccvolume_1_target_index = 0;
#line 70 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccvolume_1_target_x = 0;
#line 70 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccvolume_1_target_y = 0;
#line 70 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccvolume_1_target_z = 0;
#line 70 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccvolume_1_focus_aw = 0;
#line 70 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccvolume_1_focus_ah = 0;
#line 70 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccvolume_1_focus_xw = 0;
#line 70 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccvolume_1_focus_xh = 0;
#line 70 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccvolume_1_focus_r = 0;
#line 70 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccvolume_1_p_interact = 0;
#line 70 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  if(0) strncpy(mccvolume_1_mask_string, 0 ? 0 : "", 16384); else mccvolume_1_mask_string[0]='\0';
#line 70 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  if(0) strncpy(mccvolume_1_mask_setting, 0 ? 0 : "", 16384); else mccvolume_1_mask_setting[0]='\0';
#line 70 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccvolume_1_number_of_activations = 1;
#line 7372 "./Union_test_mask.c"

  SIG_MESSAGE("volume_1 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 69 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
    (0)*DEG2RAD,
#line 69 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
    (0)*DEG2RAD,
#line 69 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
    (0)*DEG2RAD);
#line 7382 "./Union_test_mask.c"
  rot_mul(mctr1, mcrotabeam_center, mcrotavolume_1);
  rot_transpose(mcrotabeam_center, mctr1);
  rot_mul(mcrotavolume_1, mctr1, mcrotrvolume_1);
  mctc1 = coords_set(
#line 68 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
    0,
#line 68 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
    0,
#line 68 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
    0);
#line 7393 "./Union_test_mask.c"
  rot_transpose(mcrotabeam_center, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposavolume_1 = coords_add(mcposabeam_center, mctc2);
  mctc1 = coords_sub(mcposabeam_center, mcposavolume_1);
  mcposrvolume_1 = rot_apply(mcrotavolume_1, mctc1);
  mcDEBUG_COMPONENT("volume_1", mcposavolume_1, mcrotavolume_1)
  mccomp_posa[6] = mcposavolume_1;
  mccomp_posr[6] = mcposrvolume_1;
  mcNCounter[6]  = mcPCounter[6] = mcP2Counter[6] = 0;
  mcAbsorbProp[6]= 0;
    /* Component volume_2. */
  /* Setting parameters for component volume_2. */
  SIG_MESSAGE("volume_2 (Init:SetPar)");
#line 71 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  if("test_material") strncpy(mccvolume_2_material_string, "test_material" ? "test_material" : "", 16384); else mccvolume_2_material_string[0]='\0';
#line 71 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccvolume_2_priority = 11;
#line 71 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccvolume_2_xwidth = 0.09;
#line 71 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccvolume_2_yheight = 0.1;
#line 71 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccvolume_2_zdepth = 0.1;
#line 70 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccvolume_2_xwidth2 = -1;
#line 70 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccvolume_2_yheight2 = -1;
#line 70 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccvolume_2_visualize = 1;
#line 70 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccvolume_2_target_index = 0;
#line 70 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccvolume_2_target_x = 0;
#line 70 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccvolume_2_target_y = 0;
#line 70 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccvolume_2_target_z = 0;
#line 70 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccvolume_2_focus_aw = 0;
#line 70 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccvolume_2_focus_ah = 0;
#line 70 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccvolume_2_focus_xw = 0;
#line 70 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccvolume_2_focus_xh = 0;
#line 70 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccvolume_2_focus_r = 0;
#line 70 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccvolume_2_p_interact = 0;
#line 70 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  if(0) strncpy(mccvolume_2_mask_string, 0 ? 0 : "", 16384); else mccvolume_2_mask_string[0]='\0';
#line 70 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  if(0) strncpy(mccvolume_2_mask_setting, 0 ? 0 : "", 16384); else mccvolume_2_mask_setting[0]='\0';
#line 70 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccvolume_2_number_of_activations = 1;
#line 7449 "./Union_test_mask.c"

  SIG_MESSAGE("volume_2 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 73 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
    (0)*DEG2RAD,
#line 73 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
    (0)*DEG2RAD,
#line 73 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
    (0)*DEG2RAD);
#line 7459 "./Union_test_mask.c"
  rot_mul(mctr1, mcrotabeam_center, mcrotavolume_2);
  rot_transpose(mcrotavolume_1, mctr1);
  rot_mul(mcrotavolume_2, mctr1, mcrotrvolume_2);
  mctc1 = coords_set(
#line 72 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
    0,
#line 72 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
    0.01,
#line 72 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
    0.01);
#line 7470 "./Union_test_mask.c"
  rot_transpose(mcrotabeam_center, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposavolume_2 = coords_add(mcposabeam_center, mctc2);
  mctc1 = coords_sub(mcposavolume_1, mcposavolume_2);
  mcposrvolume_2 = rot_apply(mcrotavolume_2, mctc1);
  mcDEBUG_COMPONENT("volume_2", mcposavolume_2, mcrotavolume_2)
  mccomp_posa[7] = mcposavolume_2;
  mccomp_posr[7] = mcposrvolume_2;
  mcNCounter[7]  = mcPCounter[7] = mcP2Counter[7] = 0;
  mcAbsorbProp[7]= 0;
    /* Component mask_1_1. */
  /* Setting parameters for component mask_1_1. */
  SIG_MESSAGE("mask_1_1 (Init:SetPar)");
#line 70 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  if(0) strncpy(mccmask_1_1_material_string, 0 ? 0 : "", 16384); else mccmask_1_1_material_string[0]='\0';
#line 75 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccmask_1_1_priority = 12;
#line 75 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccmask_1_1_xwidth = 0.03;
#line 75 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccmask_1_1_yheight = 0.03;
#line 75 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccmask_1_1_zdepth = 0.03;
#line 70 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccmask_1_1_xwidth2 = -1;
#line 70 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccmask_1_1_yheight2 = -1;
#line 70 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccmask_1_1_visualize = 1;
#line 70 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccmask_1_1_target_index = 0;
#line 70 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccmask_1_1_target_x = 0;
#line 70 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccmask_1_1_target_y = 0;
#line 70 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccmask_1_1_target_z = 0;
#line 70 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccmask_1_1_focus_aw = 0;
#line 70 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccmask_1_1_focus_ah = 0;
#line 70 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccmask_1_1_focus_xw = 0;
#line 70 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccmask_1_1_focus_xh = 0;
#line 70 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccmask_1_1_focus_r = 0;
#line 70 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccmask_1_1_p_interact = 0;
#line 75 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  if("volume_1") strncpy(mccmask_1_1_mask_string, "volume_1" ? "volume_1" : "", 16384); else mccmask_1_1_mask_string[0]='\0';
#line 75 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  if("ALL") strncpy(mccmask_1_1_mask_setting, "ALL" ? "ALL" : "", 16384); else mccmask_1_1_mask_setting[0]='\0';
#line 70 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccmask_1_1_number_of_activations = 1;
#line 7526 "./Union_test_mask.c"

  SIG_MESSAGE("mask_1_1 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 77 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
    (0)*DEG2RAD,
#line 77 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
    (0)*DEG2RAD,
#line 77 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
    (0)*DEG2RAD);
#line 7536 "./Union_test_mask.c"
  rot_mul(mctr1, mcrotabeam_center, mcrotamask_1_1);
  rot_transpose(mcrotavolume_2, mctr1);
  rot_mul(mcrotamask_1_1, mctr1, mcrotrmask_1_1);
  mctc1 = coords_set(
#line 76 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
    0,
#line 76 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
    0.0,
#line 76 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
    0.0);
#line 7547 "./Union_test_mask.c"
  rot_transpose(mcrotabeam_center, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposamask_1_1 = coords_add(mcposabeam_center, mctc2);
  mctc1 = coords_sub(mcposavolume_2, mcposamask_1_1);
  mcposrmask_1_1 = rot_apply(mcrotamask_1_1, mctc1);
  mcDEBUG_COMPONENT("mask_1_1", mcposamask_1_1, mcrotamask_1_1)
  mccomp_posa[8] = mcposamask_1_1;
  mccomp_posr[8] = mcposrmask_1_1;
  mcNCounter[8]  = mcPCounter[8] = mcP2Counter[8] = 0;
  mcAbsorbProp[8]= 0;
    /* Component mask_1_2. */
  /* Setting parameters for component mask_1_2. */
  SIG_MESSAGE("mask_1_2 (Init:SetPar)");
#line 70 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  if(0) strncpy(mccmask_1_2_material_string, 0 ? 0 : "", 16384); else mccmask_1_2_material_string[0]='\0';
#line 79 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccmask_1_2_priority = 13;
#line 79 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccmask_1_2_xwidth = 0.01;
#line 79 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccmask_1_2_yheight = 0.01;
#line 79 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccmask_1_2_zdepth = 0.05;
#line 70 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccmask_1_2_xwidth2 = -1;
#line 70 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccmask_1_2_yheight2 = -1;
#line 70 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccmask_1_2_visualize = 1;
#line 70 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccmask_1_2_target_index = 0;
#line 70 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccmask_1_2_target_x = 0;
#line 70 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccmask_1_2_target_y = 0;
#line 70 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccmask_1_2_target_z = 0;
#line 70 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccmask_1_2_focus_aw = 0;
#line 70 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccmask_1_2_focus_ah = 0;
#line 70 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccmask_1_2_focus_xw = 0;
#line 70 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccmask_1_2_focus_xh = 0;
#line 70 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccmask_1_2_focus_r = 0;
#line 70 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccmask_1_2_p_interact = 0;
#line 79 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  if("volume_1") strncpy(mccmask_1_2_mask_string, "volume_1" ? "volume_1" : "", 16384); else mccmask_1_2_mask_string[0]='\0';
#line 79 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  if("ALL") strncpy(mccmask_1_2_mask_setting, "ALL" ? "ALL" : "", 16384); else mccmask_1_2_mask_setting[0]='\0';
#line 70 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccmask_1_2_number_of_activations = 1;
#line 7603 "./Union_test_mask.c"

  SIG_MESSAGE("mask_1_2 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 81 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
    (0)*DEG2RAD,
#line 81 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
    (0)*DEG2RAD,
#line 81 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
    (0)*DEG2RAD);
#line 7613 "./Union_test_mask.c"
  rot_mul(mctr1, mcrotabeam_center, mcrotamask_1_2);
  rot_transpose(mcrotamask_1_1, mctr1);
  rot_mul(mcrotamask_1_2, mctr1, mcrotrmask_1_2);
  mctc1 = coords_set(
#line 80 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
    0,
#line 80 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
    0.0,
#line 80 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
    0);
#line 7624 "./Union_test_mask.c"
  rot_transpose(mcrotabeam_center, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposamask_1_2 = coords_add(mcposabeam_center, mctc2);
  mctc1 = coords_sub(mcposamask_1_1, mcposamask_1_2);
  mcposrmask_1_2 = rot_apply(mcrotamask_1_2, mctc1);
  mcDEBUG_COMPONENT("mask_1_2", mcposamask_1_2, mcrotamask_1_2)
  mccomp_posa[9] = mcposamask_1_2;
  mccomp_posr[9] = mcposrmask_1_2;
  mcNCounter[9]  = mcPCounter[9] = mcP2Counter[9] = 0;
  mcAbsorbProp[9]= 0;
    /* Component mask_2_1. */
  /* Setting parameters for component mask_2_1. */
  SIG_MESSAGE("mask_2_1 (Init:SetPar)");
#line 70 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  if(0) strncpy(mccmask_2_1_material_string, 0 ? 0 : "", 16384); else mccmask_2_1_material_string[0]='\0';
#line 83 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccmask_2_1_priority = 14;
#line 83 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccmask_2_1_xwidth = 0.04;
#line 83 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccmask_2_1_yheight = 0.04;
#line 83 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccmask_2_1_zdepth = 0.04;
#line 70 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccmask_2_1_xwidth2 = -1;
#line 70 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccmask_2_1_yheight2 = -1;
#line 70 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccmask_2_1_visualize = 1;
#line 70 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccmask_2_1_target_index = 0;
#line 70 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccmask_2_1_target_x = 0;
#line 70 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccmask_2_1_target_y = 0;
#line 70 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccmask_2_1_target_z = 0;
#line 70 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccmask_2_1_focus_aw = 0;
#line 70 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccmask_2_1_focus_ah = 0;
#line 70 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccmask_2_1_focus_xw = 0;
#line 70 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccmask_2_1_focus_xh = 0;
#line 70 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccmask_2_1_focus_r = 0;
#line 70 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccmask_2_1_p_interact = 0;
#line 83 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  if("volume_2") strncpy(mccmask_2_1_mask_string, "volume_2" ? "volume_2" : "", 16384); else mccmask_2_1_mask_string[0]='\0';
#line 83 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  if("ALL") strncpy(mccmask_2_1_mask_setting, "ALL" ? "ALL" : "", 16384); else mccmask_2_1_mask_setting[0]='\0';
#line 70 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccmask_2_1_number_of_activations = 1;
#line 7680 "./Union_test_mask.c"

  SIG_MESSAGE("mask_2_1 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 85 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
    (0)*DEG2RAD,
#line 85 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
    (0)*DEG2RAD,
#line 85 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
    (0)*DEG2RAD);
#line 7690 "./Union_test_mask.c"
  rot_mul(mctr1, mcrotabeam_center, mcrotamask_2_1);
  rot_transpose(mcrotamask_1_2, mctr1);
  rot_mul(mcrotamask_2_1, mctr1, mcrotrmask_2_1);
  mctc1 = coords_set(
#line 84 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
    0,
#line 84 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
    0.0,
#line 84 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
    0.0);
#line 7701 "./Union_test_mask.c"
  rot_transpose(mcrotabeam_center, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposamask_2_1 = coords_add(mcposabeam_center, mctc2);
  mctc1 = coords_sub(mcposamask_1_2, mcposamask_2_1);
  mcposrmask_2_1 = rot_apply(mcrotamask_2_1, mctc1);
  mcDEBUG_COMPONENT("mask_2_1", mcposamask_2_1, mcrotamask_2_1)
  mccomp_posa[10] = mcposamask_2_1;
  mccomp_posr[10] = mcposrmask_2_1;
  mcNCounter[10]  = mcPCounter[10] = mcP2Counter[10] = 0;
  mcAbsorbProp[10]= 0;
    /* Component mask_2_2. */
  /* Setting parameters for component mask_2_2. */
  SIG_MESSAGE("mask_2_2 (Init:SetPar)");
#line 70 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  if(0) strncpy(mccmask_2_2_material_string, 0 ? 0 : "", 16384); else mccmask_2_2_material_string[0]='\0';
#line 87 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccmask_2_2_priority = 15;
#line 87 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccmask_2_2_xwidth = 0.015;
#line 87 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccmask_2_2_yheight = 0.015;
#line 87 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccmask_2_2_zdepth = 0.055;
#line 70 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccmask_2_2_xwidth2 = -1;
#line 70 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccmask_2_2_yheight2 = -1;
#line 70 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccmask_2_2_visualize = 1;
#line 70 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccmask_2_2_target_index = 0;
#line 70 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccmask_2_2_target_x = 0;
#line 70 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccmask_2_2_target_y = 0;
#line 70 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccmask_2_2_target_z = 0;
#line 70 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccmask_2_2_focus_aw = 0;
#line 70 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccmask_2_2_focus_ah = 0;
#line 70 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccmask_2_2_focus_xw = 0;
#line 70 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccmask_2_2_focus_xh = 0;
#line 70 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccmask_2_2_focus_r = 0;
#line 70 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccmask_2_2_p_interact = 0;
#line 87 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  if("volume_2") strncpy(mccmask_2_2_mask_string, "volume_2" ? "volume_2" : "", 16384); else mccmask_2_2_mask_string[0]='\0';
#line 87 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  if("ALL") strncpy(mccmask_2_2_mask_setting, "ALL" ? "ALL" : "", 16384); else mccmask_2_2_mask_setting[0]='\0';
#line 70 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccmask_2_2_number_of_activations = 1;
#line 7757 "./Union_test_mask.c"

  SIG_MESSAGE("mask_2_2 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 89 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
    (0)*DEG2RAD,
#line 89 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
    (0)*DEG2RAD,
#line 89 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
    (0)*DEG2RAD);
#line 7767 "./Union_test_mask.c"
  rot_mul(mctr1, mcrotabeam_center, mcrotamask_2_2);
  rot_transpose(mcrotamask_2_1, mctr1);
  rot_mul(mcrotamask_2_2, mctr1, mcrotrmask_2_2);
  mctc1 = coords_set(
#line 88 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
    0,
#line 88 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
    0.0,
#line 88 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
    0.0);
#line 7778 "./Union_test_mask.c"
  rot_transpose(mcrotabeam_center, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposamask_2_2 = coords_add(mcposabeam_center, mctc2);
  mctc1 = coords_sub(mcposamask_2_1, mcposamask_2_2);
  mcposrmask_2_2 = rot_apply(mcrotamask_2_2, mctc1);
  mcDEBUG_COMPONENT("mask_2_2", mcposamask_2_2, mcrotamask_2_2)
  mccomp_posa[11] = mcposamask_2_2;
  mccomp_posr[11] = mcposrmask_2_2;
  mcNCounter[11]  = mcPCounter[11] = mcP2Counter[11] = 0;
  mcAbsorbProp[11]= 0;
    /* Component test_sample. */
  /* Setting parameters for component test_sample. */
  SIG_MESSAGE("test_sample (Init:SetPar)");
#line 45 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mcctest_sample_allow_inside_start = 0;
#line 45 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mcctest_sample_history_limit = 300000;
#line 45 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mcctest_sample_enable_conditionals = 1;
#line 45 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mcctest_sample_inherit_number_of_scattering_events = 0;
#line 7800 "./Union_test_mask.c"

  SIG_MESSAGE("test_sample (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 95 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
    (0)*DEG2RAD,
#line 95 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
    (0)*DEG2RAD,
#line 95 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
    (0)*DEG2RAD);
#line 7810 "./Union_test_mask.c"
  rot_mul(mctr1, mcrotabeam_center, mcrotatest_sample);
  rot_transpose(mcrotamask_2_2, mctr1);
  rot_mul(mcrotatest_sample, mctr1, mcrotrtest_sample);
  mctc1 = coords_set(
#line 94 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
    0,
#line 94 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
    0,
#line 94 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
    0);
#line 7821 "./Union_test_mask.c"
  rot_transpose(mcrotabeam_center, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposatest_sample = coords_add(mcposabeam_center, mctc2);
  mctc1 = coords_sub(mcposamask_2_2, mcposatest_sample);
  mcposrtest_sample = rot_apply(mcrotatest_sample, mctc1);
  mcDEBUG_COMPONENT("test_sample", mcposatest_sample, mcrotatest_sample)
  mccomp_posa[12] = mcposatest_sample;
  mccomp_posr[12] = mcposrtest_sample;
  mcNCounter[12]  = mcPCounter[12] = mcP2Counter[12] = 0;
  mcAbsorbProp[12]= 0;
    /* Component m4pi. */
  /* Setting parameters for component m4pi. */
  SIG_MESSAGE("m4pi (Init:SetPar)");
#line 103 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  if("Events") strncpy(mccm4pi_filename, "Events" ? "Events" : "", 16384); else mccm4pi_filename[0]='\0';
#line 103 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccm4pi_radius = 1;
#line 103 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccm4pi_restore_neutron = 1;
#line 49 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccm4pi_nowritefile = 0;
#line 7843 "./Union_test_mask.c"

  SIG_MESSAGE("m4pi (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 7850 "./Union_test_mask.c"
  rot_mul(mctr1, mcrotabeam_center, mcrotam4pi);
  rot_transpose(mcrotatest_sample, mctr1);
  rot_mul(mcrotam4pi, mctr1, mcrotrm4pi);
  mctc1 = coords_set(
#line 104 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
    0,
#line 104 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
    0,
#line 104 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
    0);
#line 7861 "./Union_test_mask.c"
  rot_transpose(mcrotabeam_center, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposam4pi = coords_add(mcposabeam_center, mctc2);
  mctc1 = coords_sub(mcposatest_sample, mcposam4pi);
  mcposrm4pi = rot_apply(mcrotam4pi, mctc1);
  mcDEBUG_COMPONENT("m4pi", mcposam4pi, mcrotam4pi)
  mccomp_posa[13] = mcposam4pi;
  mccomp_posr[13] = mcposrm4pi;
  mcNCounter[13]  = mcPCounter[13] = mcP2Counter[13] = 0;
  mcAbsorbProp[13]= 0;
    /* Component m4pi_1. */
  /* Setting parameters for component m4pi_1. */
  SIG_MESSAGE("m4pi_1 (Init:SetPar)");
#line 106 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  if("Events_1") strncpy(mccm4pi_1_filename, "Events_1" ? "Events_1" : "", 16384); else mccm4pi_1_filename[0]='\0';
#line 106 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccm4pi_1_radius = 1;
#line 106 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccm4pi_1_restore_neutron = 1;
#line 49 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccm4pi_1_nowritefile = 0;
#line 7883 "./Union_test_mask.c"

  SIG_MESSAGE("m4pi_1 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 7890 "./Union_test_mask.c"
  rot_mul(mctr1, mcrotabeam_center, mcrotam4pi_1);
  rot_transpose(mcrotam4pi, mctr1);
  rot_mul(mcrotam4pi_1, mctr1, mcrotrm4pi_1);
  mctc1 = coords_set(
#line 107 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
    0,
#line 107 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
    0,
#line 107 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
    0);
#line 7901 "./Union_test_mask.c"
  rot_transpose(mcrotabeam_center, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposam4pi_1 = coords_add(mcposabeam_center, mctc2);
  mctc1 = coords_sub(mcposam4pi, mcposam4pi_1);
  mcposrm4pi_1 = rot_apply(mcrotam4pi_1, mctc1);
  mcDEBUG_COMPONENT("m4pi_1", mcposam4pi_1, mcrotam4pi_1)
  mccomp_posa[14] = mcposam4pi_1;
  mccomp_posr[14] = mcposrm4pi_1;
  mcNCounter[14]  = mcPCounter[14] = mcP2Counter[14] = 0;
  mcAbsorbProp[14]= 0;
    /* Component m4pi_2. */
  /* Setting parameters for component m4pi_2. */
  SIG_MESSAGE("m4pi_2 (Init:SetPar)");
#line 109 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  if("Events_2") strncpy(mccm4pi_2_filename, "Events_2" ? "Events_2" : "", 16384); else mccm4pi_2_filename[0]='\0';
#line 109 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccm4pi_2_radius = 1;
#line 109 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccm4pi_2_restore_neutron = 1;
#line 49 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccm4pi_2_nowritefile = 0;
#line 7923 "./Union_test_mask.c"

  SIG_MESSAGE("m4pi_2 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 7930 "./Union_test_mask.c"
  rot_mul(mctr1, mcrotabeam_center, mcrotam4pi_2);
  rot_transpose(mcrotam4pi_1, mctr1);
  rot_mul(mcrotam4pi_2, mctr1, mcrotrm4pi_2);
  mctc1 = coords_set(
#line 110 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
    0,
#line 110 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
    0,
#line 110 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
    0);
#line 7941 "./Union_test_mask.c"
  rot_transpose(mcrotabeam_center, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposam4pi_2 = coords_add(mcposabeam_center, mctc2);
  mctc1 = coords_sub(mcposam4pi_1, mcposam4pi_2);
  mcposrm4pi_2 = rot_apply(mcrotam4pi_2, mctc1);
  mcDEBUG_COMPONENT("m4pi_2", mcposam4pi_2, mcrotam4pi_2)
  mccomp_posa[15] = mcposam4pi_2;
  mccomp_posr[15] = mcposrm4pi_2;
  mcNCounter[15]  = mcPCounter[15] = mcP2Counter[15] = 0;
  mcAbsorbProp[15]= 0;
    /* Component arm_1. */
  /* Setting parameters for component arm_1. */
  SIG_MESSAGE("arm_1 (Init:SetPar)");

  SIG_MESSAGE("arm_1 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 7961 "./Union_test_mask.c"
  rot_mul(mctr1, mcrotabeam_center, mcrotaarm_1);
  rot_transpose(mcrotam4pi_2, mctr1);
  rot_mul(mcrotaarm_1, mctr1, mcrotrarm_1);
  mctc1 = coords_set(
#line 113 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
    0,
#line 113 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
    0,
#line 113 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
    0);
#line 7972 "./Union_test_mask.c"
  rot_transpose(mcrotabeam_center, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaarm_1 = coords_add(mcposabeam_center, mctc2);
  mctc1 = coords_sub(mcposam4pi_2, mcposaarm_1);
  mcposrarm_1 = rot_apply(mcrotaarm_1, mctc1);
  mcDEBUG_COMPONENT("arm_1", mcposaarm_1, mcrotaarm_1)
  mccomp_posa[16] = mcposaarm_1;
  mccomp_posr[16] = mcposrarm_1;
  mcNCounter[16]  = mcPCounter[16] = mcP2Counter[16] = 0;
  mcAbsorbProp[16]= 0;
    /* Component m4pi_both. */
  /* Setting parameters for component m4pi_both. */
  SIG_MESSAGE("m4pi_both (Init:SetPar)");
#line 119 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  if("Events_12") strncpy(mccm4pi_both_filename, "Events_12" ? "Events_12" : "", 16384); else mccm4pi_both_filename[0]='\0';
#line 119 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccm4pi_both_radius = 1;
#line 49 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccm4pi_both_restore_neutron = 0;
#line 49 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
  mccm4pi_both_nowritefile = 0;
#line 7994 "./Union_test_mask.c"

  SIG_MESSAGE("m4pi_both (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 8001 "./Union_test_mask.c"
  rot_mul(mctr1, mcrotabeam_center, mcrotam4pi_both);
  rot_transpose(mcrotaarm_1, mctr1);
  rot_mul(mcrotam4pi_both, mctr1, mcrotrm4pi_both);
  mctc1 = coords_set(
#line 120 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
    0,
#line 120 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
    0,
#line 120 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
    0);
#line 8012 "./Union_test_mask.c"
  rot_transpose(mcrotabeam_center, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposam4pi_both = coords_add(mcposabeam_center, mctc2);
  mctc1 = coords_sub(mcposaarm_1, mcposam4pi_both);
  mcposrm4pi_both = rot_apply(mcrotam4pi_both, mctc1);
  mcDEBUG_COMPONENT("m4pi_both", mcposam4pi_both, mcrotam4pi_both)
  mccomp_posa[17] = mcposam4pi_both;
  mccomp_posr[17] = mcposrm4pi_both;
  mcNCounter[17]  = mcPCounter[17] = mcP2Counter[17] = 0;
  mcAbsorbProp[17]= 0;
  /* Component initializations. */
  /* Initializations for component test_incoherent. */
  SIG_MESSAGE("test_incoherent (Init)");
#define mccompcurname  test_incoherent
#define mccompcurtype  Incoherent_process
#define mccompcurindex 1
#define This_process mcctest_incoherent_This_process
#define Incoherent_storage mcctest_incoherent_Incoherent_storage
#define effective_my_scattering mcctest_incoherent_effective_my_scattering
#define sigma mcctest_incoherent_sigma
#define f_QE mcctest_incoherent_f_QE
#define gamma mcctest_incoherent_gamma
#define packing_factor mcctest_incoherent_packing_factor
#define unit_cell_volume mcctest_incoherent_unit_cell_volume
#define interact_fraction mcctest_incoherent_interact_fraction
#line 139 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../contrib/union/Incoherent_process.comp"
{
  // Initialize done in the component
  effective_my_scattering = ((packing_factor/unit_cell_volume) * 100 * sigma);
  Incoherent_storage.my_scattering = effective_my_scattering;
  
  Incoherent_storage.QE_sampling_frequency = f_QE;
  Incoherent_storage.lorentzian_width = gamma;

  // Need to specify if this process is isotropic
  This_process.non_isotropic_rot_index = -1; // Yes (powder)
  //This_process.non_isotropic_rot_index =  1;  // No (single crystal)

  // Packing the data into a structure that is transported to the main component
  sprintf(This_process.name,NAME_CURRENT_COMP);
  This_process.process_p_interact = interact_fraction;
  This_process.data_transfer.pointer_to_a_Incoherent_physics_storage_struct = &Incoherent_storage;
  //This_process.data_transfer.pointer_to_a_Incoherent_physics_storage_struct->my_scattering = effective_my_scattering;
  This_process.probability_for_scattering_function = &Incoherent_physics_my;
  This_process.scattering_function = &Incoherent_physics_scattering;

  // This will be the same for all process's, and can thus be moved to an include.
  sprintf(global_process_element.name,NAME_CURRENT_COMP);
  global_process_element.component_index = INDEX_CURRENT_COMP;
  global_process_element.p_scattering_process = &This_process;
  add_element_to_process_list(&global_process_list,global_process_element);
}
#line 8065 "./Union_test_mask.c"
#undef interact_fraction
#undef unit_cell_volume
#undef packing_factor
#undef gamma
#undef f_QE
#undef sigma
#undef effective_my_scattering
#undef Incoherent_storage
#undef This_process
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component test_material. */
  SIG_MESSAGE("test_material (Init)");
#define mccompcurname  test_material
#define mccompcurtype  Union_make_material
#define mccompcurindex 2
#define loop_index mcctest_material_loop_index
#define this_material mcctest_material_this_material
#define accepted_processes mcctest_material_accepted_processes
#define global_material_element mcctest_material_global_material_element
#define process_string mcctest_material_process_string
#define my_absorption mcctest_material_my_absorption
#define absorber mcctest_material_absorber
#line 193 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../contrib/union/Union_make_material.comp"
{

  /*
  // Comma test
  printf("Starting comma test on string: %s \n",process_string);
  printf("Number of commas in string: %d \n",count_commas(process_string));
  exit(1);
  */
  

  if (0 == strcmp(NAME_CURRENT_COMP,"vacuum") || 0 == strcmp(NAME_CURRENT_COMP,"Vacuum")) {
    printf("ERROR, a Union material may not be called Vacuum. A vacuum volume may be created by material=\"Vacuum\" in a geometry component.\n");
    exit(1);
  }
  if (0 == strcmp(NAME_CURRENT_COMP,"exit") || 0 == strcmp(NAME_CURRENT_COMP,"Exit")) {
    printf("ERROR, a Union material may not be called Exit. A exit volume may be created by material=\"Exit\" in a geometry component.\n");
    exit(1);
  }
  if (my_absorption < 0) {
    printf("ERROR, Union make material named %s have a negative absorption cross section!.\n",NAME_CURRENT_COMP);
    exit(1);
  }
  
  if (absorber == 0) {
    if (process_string && strlen(process_string) && strcmp(process_string,"NULL") && strcmp(process_string, "0")) {
        manual_linking_function_material(process_string, &global_process_list, &accepted_processes, NAME_CURRENT_COMP);
    } else {
      for (loop_index=0;loop_index<global_process_list.num_elements;loop_index++) {
        // printf("Automatic linking chosen [loop index = %d] with process_string = %s  \n",loop_index,process_string);
        // automatic linking
        // accept a process if index is between current and former index of make_material components
          if (1 == automatic_linking_materials_function(global_process_list.elements[loop_index],global_material_list,INDEX_CURRENT_COMP))
              add_element_to_int_list(&accepted_processes,loop_index);
      }
    }
  }

  this_material.number_of_processes = accepted_processes.num_elements; // Add number of processes
  this_material.is_vacuum = 0; // This material is not vacuum
  
  if (this_material.number_of_processes == 0 && my_absorption == 0) {
    printf("ERROR, the material named %s has no processes assigned and no absorption cross section, making it eqvialent to vacuum. Vacuums are assigned by setting material=\"Vacuum\" in a geometry component.\n",NAME_CURRENT_COMP);
    exit(1);
  }
  
  // add process element to this_material, building an array of processes called p_scattering_array
  if (this_material.number_of_processes > 0) this_material.p_scattering_array = malloc(this_material.number_of_processes * sizeof(struct scattering_process_struct));
  for (loop_index=0;loop_index<accepted_processes.num_elements;loop_index++) {
        this_material.p_scattering_array[loop_index]=*global_process_list.elements[accepted_processes.elements[loop_index]].p_scattering_process;
  }
  
  this_material.my_a = my_absorption;  // add the absorption to this material
  sprintf(this_material.name,NAME_CURRENT_COMP);
  
  // packing the information into the global_material_element, which is then included in the global_material_list.
  sprintf(global_material_element.name,NAME_CURRENT_COMP);
  global_material_element.component_index = INDEX_CURRENT_COMP;
  global_material_element.physics = &this_material; // Would be nicer if this material was a pointer, now we have the (small) data two places
  add_element_to_material_list(&global_material_list,global_material_element);
}
#line 8152 "./Union_test_mask.c"
#undef absorber
#undef my_absorption
#undef process_string
#undef global_material_element
#undef accepted_processes
#undef this_material
#undef loop_index
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component a1. */
  SIG_MESSAGE("a1 (Init)");
#define mccompcurname  a1
#define mccompcurtype  Progress_bar
#define mccompcurindex 3
#define IntermediateCnts mcca1_IntermediateCnts
#define StartTime mcca1_StartTime
#define EndTime mcca1_EndTime
#define CurrentTime mcca1_CurrentTime
#define profile mcca1_profile
#define percent mcca1_percent
#define flag_save mcca1_flag_save
#define minutes mcca1_minutes
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
#line 8189 "./Union_test_mask.c"
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

  /* Initializations for component source. */
  SIG_MESSAGE("source (Init)");
#define mccompcurname  source
#define mccompcurtype  Source_div
#define mccompcurindex 4
#define thetah mccsource_thetah
#define thetav mccsource_thetav
#define sigmah mccsource_sigmah
#define sigmav mccsource_sigmav
#define tan_h mccsource_tan_h
#define tan_v mccsource_tan_v
#define p_init mccsource_p_init
#define dist mccsource_dist
#define focus_xw mccsource_focus_xw
#define focus_yh mccsource_focus_yh
#define xwidth mccsource_xwidth
#define yheight mccsource_yheight
#define focus_aw mccsource_focus_aw
#define focus_ah mccsource_focus_ah
#define E0 mccsource_E0
#define dE mccsource_dE
#define lambda0 mccsource_lambda0
#define dlambda mccsource_dlambda
#define gauss mccsource_gauss
#define flux mccsource_flux
#line 72 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../sources/Source_div.comp"
{
sigmah = DEG2RAD*focus_aw/(sqrt(8.0*log(2.0)));
  sigmav = DEG2RAD*focus_ah/(sqrt(8.0*log(2.0)));

  if (xwidth < 0 || yheight < 0 || focus_aw < 0 || focus_ah < 0) {
      printf("Source_div: %s: Error in input parameter values!\n"
             "ERROR       Exiting\n",
           NAME_CURRENT_COMP);
      exit(0);
  }
  if ((!lambda0 && !E0 && !dE && !dlambda)) {
    printf("Source_div: %s: You must specify either a wavelength or energy range!\n ERROR - Exiting\n",
           NAME_CURRENT_COMP);
    exit(0);
  }
  if ((!lambda0 && !dlambda && (E0 <= 0 || dE < 0 || E0-dE <= 0))
    || (!E0 && !dE && (lambda0 <= 0 || dlambda < 0 || lambda0-dlambda <= 0))) {
    printf("Source_div: %s: Unmeaningful definition of wavelength or energy range!\n ERROR - Exiting\n",
           NAME_CURRENT_COMP);
      exit(0);
  }
  /* compute distance to next component */
  Coords ToTarget;
  double tx,ty,tz;
  ToTarget = coords_sub(POS_A_COMP_INDEX(INDEX_CURRENT_COMP+1),POS_A_CURRENT_COMP);
  ToTarget = rot_apply(ROT_A_CURRENT_COMP, ToTarget);
  coords_get(ToTarget, &tx, &ty, &tz);
  dist=sqrt(tx*tx+ty*ty+tz*tz);
  /* compute target area */
  if (dist) {
    focus_xw=dist*tan(focus_aw*DEG2RAD);
    focus_yh=dist*tan(focus_ah*DEG2RAD);
  }

  p_init  = flux*1e4*xwidth*yheight/mcget_ncount();
  if (!focus_aw || !focus_ah)
    exit(printf("Source_div: %s: Zero divergence defined. \n"
                "ERROR       Use non zero values for focus_aw and focus_ah.\n",
           NAME_CURRENT_COMP));
  p_init *= 2*fabs(DEG2RAD*focus_aw*sin(DEG2RAD*focus_ah/2));  /* solid angle */
  if (dlambda)
    p_init *= 2*dlambda;
  else if (dE)
    p_init *= 2*dE;
}
#line 8273 "./Union_test_mask.c"
#undef flux
#undef gauss
#undef dlambda
#undef lambda0
#undef dE
#undef E0
#undef focus_ah
#undef focus_aw
#undef yheight
#undef xwidth
#undef focus_yh
#undef focus_xw
#undef dist
#undef p_init
#undef tan_v
#undef tan_h
#undef sigmav
#undef sigmah
#undef thetav
#undef thetah
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component beam_center. */
  SIG_MESSAGE("beam_center (Init)");

  /* Initializations for component volume_1. */
  SIG_MESSAGE("volume_1 (Init)");
#define mccompcurname  volume_1
#define mccompcurtype  Union_box
#define mccompcurindex 6
#define loop_index mccvolume_1_loop_index
#define this_box_volume mccvolume_1_this_box_volume
#define global_geometry_element mccvolume_1_global_geometry_element
#define this_box_storage mccvolume_1_this_box_storage
#define material_string mccvolume_1_material_string
#define priority mccvolume_1_priority
#define xwidth mccvolume_1_xwidth
#define yheight mccvolume_1_yheight
#define zdepth mccvolume_1_zdepth
#define xwidth2 mccvolume_1_xwidth2
#define yheight2 mccvolume_1_yheight2
#define visualize mccvolume_1_visualize
#define target_index mccvolume_1_target_index
#define target_x mccvolume_1_target_x
#define target_y mccvolume_1_target_y
#define target_z mccvolume_1_target_z
#define focus_aw mccvolume_1_focus_aw
#define focus_ah mccvolume_1_focus_ah
#define focus_xw mccvolume_1_focus_xw
#define focus_xh mccvolume_1_focus_xh
#define focus_r mccvolume_1_focus_r
#define p_interact mccvolume_1_p_interact
#define mask_string mccvolume_1_mask_string
#define mask_setting mccvolume_1_mask_setting
#define number_of_activations mccvolume_1_number_of_activations
#line 223 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../contrib/union/Union_box.comp"
{
// Initializes the focusing system for this volume including input sanitation.
focus_initialize(&this_box_volume.geometry, POS_A_COMP_INDEX(INDEX_CURRENT_COMP+target_index), POS_A_CURRENT_COMP, ROT_A_CURRENT_COMP, target_index, target_x, target_y, target_z, focus_aw, focus_ah, focus_xw, focus_xh, focus_r, NAME_CURRENT_COMP);

// Input sanitation for this geometry
if (xwidth <= 0) {
  printf("\nERROR in Union_box named %s, the xwidth is <= 0. \n",NAME_CURRENT_COMP);
  exit(1);
}
if (yheight <= 0) {
  printf("\nERROR in Union_box named %s, yheight is <= 0. \n",NAME_CURRENT_COMP);
  exit(1);
}
if (zdepth <= 0) {
  printf("\nERROR in Union_box named %s, zdepth is <= 0. \n",NAME_CURRENT_COMP);
  exit(1);
}
if (xwidth2 <= 0 && xwidth2 != -1) {
  printf("\nERROR in Union_box named %s, the xwidth2 is <= 0. \n",NAME_CURRENT_COMP);
  exit(1);
}
if (yheight2 <= 0 && yheight2 != -1) {
  printf("\nERROR in Union_box named %s, yheight2 is <= 0. \n",NAME_CURRENT_COMP);
  exit(1);
}

// Use sanitation
if (global_material_list.num_elements == 0) {
  printf("\nERROR: Need to define a material using Union_make_material before using a Union geometry component. \n");
  printf("       %s was defined before first use of Union_make_material.\n",NAME_CURRENT_COMP);
  exit(1);
}

this_box_volume.geometry.is_masked_volume = 0;
this_box_volume.geometry.is_exit_volume = 0;
this_box_volume.geometry.is_mask_volume = 0;
// check if the volume is a mask, if it is the material string is irelevant.
if (mask_string && strlen(mask_string) && strcmp(mask_string, "NULL") && strcmp(mask_string, "0")) {
    // A mask volume is used to limit the extend of other volumes, called the masked volumes. These are specified in the mask_string.
    // In order for a ray to enter a masked volume, it needs to be both in the region covered by that volume AND the mask volume.
    // When more than
    this_box_volume.geometry.mask_mode = 1; // Default is mask mode is ALL
    if (mask_setting && strlen(mask_setting) && strcmp(mask_setting, "NULL") && strcmp(mask_setting, "0")) {
        if (strcmp(mask_setting,"ALL") == 0 || strcmp(mask_setting,"All") == 0) this_box_volume.geometry.mask_mode = 1;
        else if (strcmp(mask_setting,"ANY") == 0 || strcmp(mask_setting,"Any") == 0) this_box_volume.geometry.mask_mode = 2;
        else {
            printf("The mask_mode of component %s is set to %s, but must be either ALL or ANY.\n",NAME_CURRENT_COMP,mask_setting);
            exit(1);
        }
    }
    
    for (loop_index=0;loop_index<global_geometry_list.num_elements;loop_index++) {
        // Add mask list
        if (1 == manual_linking_function(global_geometry_list.elements[loop_index].name,mask_string)) {
            add_element_to_int_list(&this_box_volume.geometry.mask_list,global_geometry_list.elements[loop_index].component_index);
            add_element_to_int_list(&global_geometry_list.elements[loop_index].Volume->geometry.masked_by_list,INDEX_CURRENT_COMP);
            global_geometry_list.elements[loop_index].Volume->geometry.is_masked_volume = 1;
            if (this_box_volume.geometry.mask_mode == 2)
                global_geometry_list.elements[loop_index].Volume->geometry.mask_mode = 2;
            if (this_box_volume.geometry.mask_mode == 1) {
                if (global_geometry_list.elements[loop_index].Volume->geometry.is_masked_volume == 1 && global_geometry_list.elements[loop_index].Volume->geometry.mask_mode != 2)
                    // If more than one mask is added to one volume, the ANY mode overwrites the (default) ALL mode.
                    global_geometry_list.elements[loop_index].Volume->geometry.mask_mode = 1;
            }
            
            found_geometries = 1;
        }
    }
    if (found_geometries == 0) {
        printf("The mask_string in geometry: %s did not find any of the specified volumes in the mask_string %s \n",NAME_CURRENT_COMP,mask_string);
        exit(1);
    }
    this_box_volume.p_physics = malloc(sizeof(struct physics_struct));
    this_box_volume.p_physics->is_vacuum = 0; // Makes this volume a vacuum
    this_box_volume.p_physics->number_of_processes = (int) 0; // Should not be used.
    this_box_volume.p_physics->my_a = 0; // Should not be used.
    sprintf(this_box_volume.p_physics->name,"Mask");
    this_box_volume.geometry.is_mask_volume = 1;
    
    
// Read the material input, or if it lacks, use automatic linking.
} else if (material_string && strlen(material_string) && strcmp(material_string, "NULL") && strcmp(material_string, "0")) {
    // A geometry string was given, use it to determine which material
    if (0 == strcmp(material_string,"vacuum") || 0 == strcmp(material_string,"Vacuum")) {
        // One could have a global physics struct for vacuum instead of creating one for each
        this_box_volume.p_physics = malloc(sizeof(struct physics_struct));
        this_box_volume.p_physics->is_vacuum = 1; // Makes this volume a vacuum
        this_box_volume.p_physics->number_of_processes = (int) 0; // Should not be used.
        this_box_volume.p_physics->my_a = 0; // Should not be used.
        sprintf(this_box_volume.p_physics->name,"Vacuum");
    } else if (0 == strcmp(material_string,"exit") || 0 == strcmp(material_string,"Exit")) {
        // One could have a global physics struct for vacuum instead of creating one for each
        this_box_volume.p_physics = malloc(sizeof(struct physics_struct));
        this_box_volume.p_physics->is_vacuum = 1; // Makes this volume a vacuum
        this_box_volume.p_physics->number_of_processes = (int) 0; // Should not be used.
        this_box_volume.p_physics->my_a = 0; // Should not be used.
        this_box_volume.geometry.is_exit_volume = 1;
        sprintf(this_box_volume.p_physics->name,"Exit");
    } else {
        #ifndef MATERIAL_DETECTOR
            printf("Need to define a material before refering to it in a geometry %s.\n",NAME_CURRENT_COMP);
            exit(1);
        #endif
        for (loop_index=0;loop_index<global_material_list.num_elements;loop_index++) {
            if (0 == strcmp(material_string,global_material_list.elements[loop_index].name)) {
                this_box_volume.p_physics = global_material_list.elements[loop_index].physics;
                break;
            }
            if (loop_index == global_material_list.num_elements-1) {
            printf("\n");
            printf("ERROR: The material string \"%s\" in Union geometry \"%s\" did not match a specified material. \n",material_string,NAME_CURRENT_COMP);
            printf("       The materials available at this point (need to be defined before the geometry): \n");
            for (loop_index=0;loop_index<global_material_list.num_elements;loop_index++)
              printf("         %s\n",global_material_list.elements[loop_index].name);
            printf("\n");
            printf("       It is also possible to use one of the defualt materials avaiable: \n");
            printf("           Vacuum (for a Volume without scattering or absorption)\n");
            printf("           Exit (for a Volume where the ray exits the component if it enters)\n");
            printf("           Mask (for a Volume that masks existing volumes specified in the mask_string\n");
            exit(1);
          }
        }
    }
} else {
    // Automatic linking, simply using the last defined material.
    #ifndef MATERIAL_DETECTOR
        printf("Need to define a material before the geometry to use automatic linking %s.\n",NAME_CURRENT_COMP);
        exit(1);
    #endif
    this_box_volume.p_physics = global_material_list.elements[global_material_list.num_elements-1].physics;
}

sprintf(this_box_volume.name,NAME_CURRENT_COMP);
sprintf(this_box_volume.geometry.shape,"box");
this_box_volume.geometry.priority_value = priority;
this_box_volume.geometry.geometry_p_interact = p_interact;
// Currently the coordinates will be in absolute space.
this_box_volume.geometry.center = POS_A_CURRENT_COMP;

this_box_storage.z_depth = zdepth;
this_box_storage.x_width1 = xwidth;
this_box_storage.y_height1 = yheight;

this_box_storage.is_rectangle = 0;
if (xwidth2 < 0 && yheight2 < 0) this_box_storage.is_rectangle = 1;
if (xwidth == xwidth2 && yheight == yheight2) this_box_storage.is_rectangle = 1;

if (xwidth2 < 0) {
    this_box_storage.x_width2 = xwidth;
    xwidth2 = xwidth;
} else this_box_storage.x_width2 = xwidth2;

if (yheight2 < 0) {
    this_box_storage.y_height2 = yheight;
    yheight2 = yheight;
} else this_box_storage.y_height2 = yheight2;


this_box_storage.normal_vectors[0] = coords_set(0,0,1);
this_box_storage.normal_vectors[1] = coords_set(0,0,1);

// for sides with y component = 0
x_component = 2*zdepth/sqrt((xwidth-xwidth2)*(xwidth-xwidth2)+4*zdepth*zdepth);
z_component = (xwidth-xwidth2)/sqrt(4*zdepth*zdepth+(xwidth-xwidth2)*(xwidth-xwidth2));

this_box_storage.normal_vectors[2] = coords_set(x_component,0,z_component);
this_box_storage.normal_vectors[3] = coords_set(-x_component,0,z_component);

// for sides with x component = 0
y_component = 2*zdepth/sqrt((yheight-yheight2)*(yheight-yheight2)+4*zdepth*zdepth);
z_component = (yheight-yheight2)/sqrt(4*zdepth*zdepth+(yheight-yheight2)*(yheight-yheight2));

this_box_storage.normal_vectors[4] = coords_set(0,y_component,z_component);
this_box_storage.normal_vectors[5] = coords_set(0,-y_component,z_component);

this_box_volume.geometry.visualization_on = visualize;

this_box_volume.geometry.geometry_parameters.p_box_storage = &this_box_storage;

// Assign pointers to functions for intersection with the shape, checking if a point is inside the shape
if (this_box_storage.is_rectangle == 1) {
this_box_volume.geometry.intersect_function = &sample_box_intersect_simple;
this_box_volume.geometry.within_function = &r_within_box_simple;
}
else {
this_box_volume.geometry.intersect_function = &sample_box_intersect_advanced;
this_box_volume.geometry.within_function = &r_within_box_advanced;
}

this_box_volume.geometry.shell_points = &box_shell_points;
this_box_volume.geometry.mcdisplay_function = &mcdisplay_box_function;
this_box_volume.geometry.initialize_from_main_function = &initialize_box_geometry_from_main_component;
this_box_volume.geometry.process_rot_allocated = 0;

this_box_volume.geometry.copy_geometry_parameters = &allocate_box_storage_copy;

rot_copy(this_box_volume.geometry.rotation_matrix,ROT_A_CURRENT_COMP); //check how ROT_R_CURRENT_COMP would work
rot_transpose(ROT_A_CURRENT_COMP,this_box_volume.geometry.transpose_rotation_matrix);

// Initialize loggers
this_box_volume.loggers.num_elements = 0;

// packing the information into the global_geometry_element, which is then included in the global_geometry_list.
sprintf(global_geometry_element.name,NAME_CURRENT_COMP);
global_geometry_element.activation_counter = number_of_activations;
global_geometry_element.component_index = INDEX_CURRENT_COMP;
global_geometry_element.Volume = &this_box_volume; // Would be nicer if this m was a pointer, now we have the (small) data two places
add_element_to_geometry_list(&global_geometry_list,global_geometry_element);
}
#line 8541 "./Union_test_mask.c"
#undef number_of_activations
#undef mask_setting
#undef mask_string
#undef p_interact
#undef focus_r
#undef focus_xh
#undef focus_xw
#undef focus_ah
#undef focus_aw
#undef target_z
#undef target_y
#undef target_x
#undef target_index
#undef visualize
#undef yheight2
#undef xwidth2
#undef zdepth
#undef yheight
#undef xwidth
#undef priority
#undef material_string
#undef this_box_storage
#undef global_geometry_element
#undef this_box_volume
#undef loop_index
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component volume_2. */
  SIG_MESSAGE("volume_2 (Init)");
#define mccompcurname  volume_2
#define mccompcurtype  Union_box
#define mccompcurindex 7
#define loop_index mccvolume_2_loop_index
#define this_box_volume mccvolume_2_this_box_volume
#define global_geometry_element mccvolume_2_global_geometry_element
#define this_box_storage mccvolume_2_this_box_storage
#define material_string mccvolume_2_material_string
#define priority mccvolume_2_priority
#define xwidth mccvolume_2_xwidth
#define yheight mccvolume_2_yheight
#define zdepth mccvolume_2_zdepth
#define xwidth2 mccvolume_2_xwidth2
#define yheight2 mccvolume_2_yheight2
#define visualize mccvolume_2_visualize
#define target_index mccvolume_2_target_index
#define target_x mccvolume_2_target_x
#define target_y mccvolume_2_target_y
#define target_z mccvolume_2_target_z
#define focus_aw mccvolume_2_focus_aw
#define focus_ah mccvolume_2_focus_ah
#define focus_xw mccvolume_2_focus_xw
#define focus_xh mccvolume_2_focus_xh
#define focus_r mccvolume_2_focus_r
#define p_interact mccvolume_2_p_interact
#define mask_string mccvolume_2_mask_string
#define mask_setting mccvolume_2_mask_setting
#define number_of_activations mccvolume_2_number_of_activations
#line 223 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../contrib/union/Union_box.comp"
{
// Initializes the focusing system for this volume including input sanitation.
focus_initialize(&this_box_volume.geometry, POS_A_COMP_INDEX(INDEX_CURRENT_COMP+target_index), POS_A_CURRENT_COMP, ROT_A_CURRENT_COMP, target_index, target_x, target_y, target_z, focus_aw, focus_ah, focus_xw, focus_xh, focus_r, NAME_CURRENT_COMP);

// Input sanitation for this geometry
if (xwidth <= 0) {
  printf("\nERROR in Union_box named %s, the xwidth is <= 0. \n",NAME_CURRENT_COMP);
  exit(1);
}
if (yheight <= 0) {
  printf("\nERROR in Union_box named %s, yheight is <= 0. \n",NAME_CURRENT_COMP);
  exit(1);
}
if (zdepth <= 0) {
  printf("\nERROR in Union_box named %s, zdepth is <= 0. \n",NAME_CURRENT_COMP);
  exit(1);
}
if (xwidth2 <= 0 && xwidth2 != -1) {
  printf("\nERROR in Union_box named %s, the xwidth2 is <= 0. \n",NAME_CURRENT_COMP);
  exit(1);
}
if (yheight2 <= 0 && yheight2 != -1) {
  printf("\nERROR in Union_box named %s, yheight2 is <= 0. \n",NAME_CURRENT_COMP);
  exit(1);
}

// Use sanitation
if (global_material_list.num_elements == 0) {
  printf("\nERROR: Need to define a material using Union_make_material before using a Union geometry component. \n");
  printf("       %s was defined before first use of Union_make_material.\n",NAME_CURRENT_COMP);
  exit(1);
}

this_box_volume.geometry.is_masked_volume = 0;
this_box_volume.geometry.is_exit_volume = 0;
this_box_volume.geometry.is_mask_volume = 0;
// check if the volume is a mask, if it is the material string is irelevant.
if (mask_string && strlen(mask_string) && strcmp(mask_string, "NULL") && strcmp(mask_string, "0")) {
    // A mask volume is used to limit the extend of other volumes, called the masked volumes. These are specified in the mask_string.
    // In order for a ray to enter a masked volume, it needs to be both in the region covered by that volume AND the mask volume.
    // When more than
    this_box_volume.geometry.mask_mode = 1; // Default is mask mode is ALL
    if (mask_setting && strlen(mask_setting) && strcmp(mask_setting, "NULL") && strcmp(mask_setting, "0")) {
        if (strcmp(mask_setting,"ALL") == 0 || strcmp(mask_setting,"All") == 0) this_box_volume.geometry.mask_mode = 1;
        else if (strcmp(mask_setting,"ANY") == 0 || strcmp(mask_setting,"Any") == 0) this_box_volume.geometry.mask_mode = 2;
        else {
            printf("The mask_mode of component %s is set to %s, but must be either ALL or ANY.\n",NAME_CURRENT_COMP,mask_setting);
            exit(1);
        }
    }
    
    for (loop_index=0;loop_index<global_geometry_list.num_elements;loop_index++) {
        // Add mask list
        if (1 == manual_linking_function(global_geometry_list.elements[loop_index].name,mask_string)) {
            add_element_to_int_list(&this_box_volume.geometry.mask_list,global_geometry_list.elements[loop_index].component_index);
            add_element_to_int_list(&global_geometry_list.elements[loop_index].Volume->geometry.masked_by_list,INDEX_CURRENT_COMP);
            global_geometry_list.elements[loop_index].Volume->geometry.is_masked_volume = 1;
            if (this_box_volume.geometry.mask_mode == 2)
                global_geometry_list.elements[loop_index].Volume->geometry.mask_mode = 2;
            if (this_box_volume.geometry.mask_mode == 1) {
                if (global_geometry_list.elements[loop_index].Volume->geometry.is_masked_volume == 1 && global_geometry_list.elements[loop_index].Volume->geometry.mask_mode != 2)
                    // If more than one mask is added to one volume, the ANY mode overwrites the (default) ALL mode.
                    global_geometry_list.elements[loop_index].Volume->geometry.mask_mode = 1;
            }
            
            found_geometries = 1;
        }
    }
    if (found_geometries == 0) {
        printf("The mask_string in geometry: %s did not find any of the specified volumes in the mask_string %s \n",NAME_CURRENT_COMP,mask_string);
        exit(1);
    }
    this_box_volume.p_physics = malloc(sizeof(struct physics_struct));
    this_box_volume.p_physics->is_vacuum = 0; // Makes this volume a vacuum
    this_box_volume.p_physics->number_of_processes = (int) 0; // Should not be used.
    this_box_volume.p_physics->my_a = 0; // Should not be used.
    sprintf(this_box_volume.p_physics->name,"Mask");
    this_box_volume.geometry.is_mask_volume = 1;
    
    
// Read the material input, or if it lacks, use automatic linking.
} else if (material_string && strlen(material_string) && strcmp(material_string, "NULL") && strcmp(material_string, "0")) {
    // A geometry string was given, use it to determine which material
    if (0 == strcmp(material_string,"vacuum") || 0 == strcmp(material_string,"Vacuum")) {
        // One could have a global physics struct for vacuum instead of creating one for each
        this_box_volume.p_physics = malloc(sizeof(struct physics_struct));
        this_box_volume.p_physics->is_vacuum = 1; // Makes this volume a vacuum
        this_box_volume.p_physics->number_of_processes = (int) 0; // Should not be used.
        this_box_volume.p_physics->my_a = 0; // Should not be used.
        sprintf(this_box_volume.p_physics->name,"Vacuum");
    } else if (0 == strcmp(material_string,"exit") || 0 == strcmp(material_string,"Exit")) {
        // One could have a global physics struct for vacuum instead of creating one for each
        this_box_volume.p_physics = malloc(sizeof(struct physics_struct));
        this_box_volume.p_physics->is_vacuum = 1; // Makes this volume a vacuum
        this_box_volume.p_physics->number_of_processes = (int) 0; // Should not be used.
        this_box_volume.p_physics->my_a = 0; // Should not be used.
        this_box_volume.geometry.is_exit_volume = 1;
        sprintf(this_box_volume.p_physics->name,"Exit");
    } else {
        #ifndef MATERIAL_DETECTOR
            printf("Need to define a material before refering to it in a geometry %s.\n",NAME_CURRENT_COMP);
            exit(1);
        #endif
        for (loop_index=0;loop_index<global_material_list.num_elements;loop_index++) {
            if (0 == strcmp(material_string,global_material_list.elements[loop_index].name)) {
                this_box_volume.p_physics = global_material_list.elements[loop_index].physics;
                break;
            }
            if (loop_index == global_material_list.num_elements-1) {
            printf("\n");
            printf("ERROR: The material string \"%s\" in Union geometry \"%s\" did not match a specified material. \n",material_string,NAME_CURRENT_COMP);
            printf("       The materials available at this point (need to be defined before the geometry): \n");
            for (loop_index=0;loop_index<global_material_list.num_elements;loop_index++)
              printf("         %s\n",global_material_list.elements[loop_index].name);
            printf("\n");
            printf("       It is also possible to use one of the defualt materials avaiable: \n");
            printf("           Vacuum (for a Volume without scattering or absorption)\n");
            printf("           Exit (for a Volume where the ray exits the component if it enters)\n");
            printf("           Mask (for a Volume that masks existing volumes specified in the mask_string\n");
            exit(1);
          }
        }
    }
} else {
    // Automatic linking, simply using the last defined material.
    #ifndef MATERIAL_DETECTOR
        printf("Need to define a material before the geometry to use automatic linking %s.\n",NAME_CURRENT_COMP);
        exit(1);
    #endif
    this_box_volume.p_physics = global_material_list.elements[global_material_list.num_elements-1].physics;
}

sprintf(this_box_volume.name,NAME_CURRENT_COMP);
sprintf(this_box_volume.geometry.shape,"box");
this_box_volume.geometry.priority_value = priority;
this_box_volume.geometry.geometry_p_interact = p_interact;
// Currently the coordinates will be in absolute space.
this_box_volume.geometry.center = POS_A_CURRENT_COMP;

this_box_storage.z_depth = zdepth;
this_box_storage.x_width1 = xwidth;
this_box_storage.y_height1 = yheight;

this_box_storage.is_rectangle = 0;
if (xwidth2 < 0 && yheight2 < 0) this_box_storage.is_rectangle = 1;
if (xwidth == xwidth2 && yheight == yheight2) this_box_storage.is_rectangle = 1;

if (xwidth2 < 0) {
    this_box_storage.x_width2 = xwidth;
    xwidth2 = xwidth;
} else this_box_storage.x_width2 = xwidth2;

if (yheight2 < 0) {
    this_box_storage.y_height2 = yheight;
    yheight2 = yheight;
} else this_box_storage.y_height2 = yheight2;


this_box_storage.normal_vectors[0] = coords_set(0,0,1);
this_box_storage.normal_vectors[1] = coords_set(0,0,1);

// for sides with y component = 0
x_component = 2*zdepth/sqrt((xwidth-xwidth2)*(xwidth-xwidth2)+4*zdepth*zdepth);
z_component = (xwidth-xwidth2)/sqrt(4*zdepth*zdepth+(xwidth-xwidth2)*(xwidth-xwidth2));

this_box_storage.normal_vectors[2] = coords_set(x_component,0,z_component);
this_box_storage.normal_vectors[3] = coords_set(-x_component,0,z_component);

// for sides with x component = 0
y_component = 2*zdepth/sqrt((yheight-yheight2)*(yheight-yheight2)+4*zdepth*zdepth);
z_component = (yheight-yheight2)/sqrt(4*zdepth*zdepth+(yheight-yheight2)*(yheight-yheight2));

this_box_storage.normal_vectors[4] = coords_set(0,y_component,z_component);
this_box_storage.normal_vectors[5] = coords_set(0,-y_component,z_component);

this_box_volume.geometry.visualization_on = visualize;

this_box_volume.geometry.geometry_parameters.p_box_storage = &this_box_storage;

// Assign pointers to functions for intersection with the shape, checking if a point is inside the shape
if (this_box_storage.is_rectangle == 1) {
this_box_volume.geometry.intersect_function = &sample_box_intersect_simple;
this_box_volume.geometry.within_function = &r_within_box_simple;
}
else {
this_box_volume.geometry.intersect_function = &sample_box_intersect_advanced;
this_box_volume.geometry.within_function = &r_within_box_advanced;
}

this_box_volume.geometry.shell_points = &box_shell_points;
this_box_volume.geometry.mcdisplay_function = &mcdisplay_box_function;
this_box_volume.geometry.initialize_from_main_function = &initialize_box_geometry_from_main_component;
this_box_volume.geometry.process_rot_allocated = 0;

this_box_volume.geometry.copy_geometry_parameters = &allocate_box_storage_copy;

rot_copy(this_box_volume.geometry.rotation_matrix,ROT_A_CURRENT_COMP); //check how ROT_R_CURRENT_COMP would work
rot_transpose(ROT_A_CURRENT_COMP,this_box_volume.geometry.transpose_rotation_matrix);

// Initialize loggers
this_box_volume.loggers.num_elements = 0;

// packing the information into the global_geometry_element, which is then included in the global_geometry_list.
sprintf(global_geometry_element.name,NAME_CURRENT_COMP);
global_geometry_element.activation_counter = number_of_activations;
global_geometry_element.component_index = INDEX_CURRENT_COMP;
global_geometry_element.Volume = &this_box_volume; // Would be nicer if this m was a pointer, now we have the (small) data two places
add_element_to_geometry_list(&global_geometry_list,global_geometry_element);
}
#line 8811 "./Union_test_mask.c"
#undef number_of_activations
#undef mask_setting
#undef mask_string
#undef p_interact
#undef focus_r
#undef focus_xh
#undef focus_xw
#undef focus_ah
#undef focus_aw
#undef target_z
#undef target_y
#undef target_x
#undef target_index
#undef visualize
#undef yheight2
#undef xwidth2
#undef zdepth
#undef yheight
#undef xwidth
#undef priority
#undef material_string
#undef this_box_storage
#undef global_geometry_element
#undef this_box_volume
#undef loop_index
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component mask_1_1. */
  SIG_MESSAGE("mask_1_1 (Init)");
#define mccompcurname  mask_1_1
#define mccompcurtype  Union_box
#define mccompcurindex 8
#define loop_index mccmask_1_1_loop_index
#define this_box_volume mccmask_1_1_this_box_volume
#define global_geometry_element mccmask_1_1_global_geometry_element
#define this_box_storage mccmask_1_1_this_box_storage
#define material_string mccmask_1_1_material_string
#define priority mccmask_1_1_priority
#define xwidth mccmask_1_1_xwidth
#define yheight mccmask_1_1_yheight
#define zdepth mccmask_1_1_zdepth
#define xwidth2 mccmask_1_1_xwidth2
#define yheight2 mccmask_1_1_yheight2
#define visualize mccmask_1_1_visualize
#define target_index mccmask_1_1_target_index
#define target_x mccmask_1_1_target_x
#define target_y mccmask_1_1_target_y
#define target_z mccmask_1_1_target_z
#define focus_aw mccmask_1_1_focus_aw
#define focus_ah mccmask_1_1_focus_ah
#define focus_xw mccmask_1_1_focus_xw
#define focus_xh mccmask_1_1_focus_xh
#define focus_r mccmask_1_1_focus_r
#define p_interact mccmask_1_1_p_interact
#define mask_string mccmask_1_1_mask_string
#define mask_setting mccmask_1_1_mask_setting
#define number_of_activations mccmask_1_1_number_of_activations
#line 223 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../contrib/union/Union_box.comp"
{
// Initializes the focusing system for this volume including input sanitation.
focus_initialize(&this_box_volume.geometry, POS_A_COMP_INDEX(INDEX_CURRENT_COMP+target_index), POS_A_CURRENT_COMP, ROT_A_CURRENT_COMP, target_index, target_x, target_y, target_z, focus_aw, focus_ah, focus_xw, focus_xh, focus_r, NAME_CURRENT_COMP);

// Input sanitation for this geometry
if (xwidth <= 0) {
  printf("\nERROR in Union_box named %s, the xwidth is <= 0. \n",NAME_CURRENT_COMP);
  exit(1);
}
if (yheight <= 0) {
  printf("\nERROR in Union_box named %s, yheight is <= 0. \n",NAME_CURRENT_COMP);
  exit(1);
}
if (zdepth <= 0) {
  printf("\nERROR in Union_box named %s, zdepth is <= 0. \n",NAME_CURRENT_COMP);
  exit(1);
}
if (xwidth2 <= 0 && xwidth2 != -1) {
  printf("\nERROR in Union_box named %s, the xwidth2 is <= 0. \n",NAME_CURRENT_COMP);
  exit(1);
}
if (yheight2 <= 0 && yheight2 != -1) {
  printf("\nERROR in Union_box named %s, yheight2 is <= 0. \n",NAME_CURRENT_COMP);
  exit(1);
}

// Use sanitation
if (global_material_list.num_elements == 0) {
  printf("\nERROR: Need to define a material using Union_make_material before using a Union geometry component. \n");
  printf("       %s was defined before first use of Union_make_material.\n",NAME_CURRENT_COMP);
  exit(1);
}

this_box_volume.geometry.is_masked_volume = 0;
this_box_volume.geometry.is_exit_volume = 0;
this_box_volume.geometry.is_mask_volume = 0;
// check if the volume is a mask, if it is the material string is irelevant.
if (mask_string && strlen(mask_string) && strcmp(mask_string, "NULL") && strcmp(mask_string, "0")) {
    // A mask volume is used to limit the extend of other volumes, called the masked volumes. These are specified in the mask_string.
    // In order for a ray to enter a masked volume, it needs to be both in the region covered by that volume AND the mask volume.
    // When more than
    this_box_volume.geometry.mask_mode = 1; // Default is mask mode is ALL
    if (mask_setting && strlen(mask_setting) && strcmp(mask_setting, "NULL") && strcmp(mask_setting, "0")) {
        if (strcmp(mask_setting,"ALL") == 0 || strcmp(mask_setting,"All") == 0) this_box_volume.geometry.mask_mode = 1;
        else if (strcmp(mask_setting,"ANY") == 0 || strcmp(mask_setting,"Any") == 0) this_box_volume.geometry.mask_mode = 2;
        else {
            printf("The mask_mode of component %s is set to %s, but must be either ALL or ANY.\n",NAME_CURRENT_COMP,mask_setting);
            exit(1);
        }
    }
    
    for (loop_index=0;loop_index<global_geometry_list.num_elements;loop_index++) {
        // Add mask list
        if (1 == manual_linking_function(global_geometry_list.elements[loop_index].name,mask_string)) {
            add_element_to_int_list(&this_box_volume.geometry.mask_list,global_geometry_list.elements[loop_index].component_index);
            add_element_to_int_list(&global_geometry_list.elements[loop_index].Volume->geometry.masked_by_list,INDEX_CURRENT_COMP);
            global_geometry_list.elements[loop_index].Volume->geometry.is_masked_volume = 1;
            if (this_box_volume.geometry.mask_mode == 2)
                global_geometry_list.elements[loop_index].Volume->geometry.mask_mode = 2;
            if (this_box_volume.geometry.mask_mode == 1) {
                if (global_geometry_list.elements[loop_index].Volume->geometry.is_masked_volume == 1 && global_geometry_list.elements[loop_index].Volume->geometry.mask_mode != 2)
                    // If more than one mask is added to one volume, the ANY mode overwrites the (default) ALL mode.
                    global_geometry_list.elements[loop_index].Volume->geometry.mask_mode = 1;
            }
            
            found_geometries = 1;
        }
    }
    if (found_geometries == 0) {
        printf("The mask_string in geometry: %s did not find any of the specified volumes in the mask_string %s \n",NAME_CURRENT_COMP,mask_string);
        exit(1);
    }
    this_box_volume.p_physics = malloc(sizeof(struct physics_struct));
    this_box_volume.p_physics->is_vacuum = 0; // Makes this volume a vacuum
    this_box_volume.p_physics->number_of_processes = (int) 0; // Should not be used.
    this_box_volume.p_physics->my_a = 0; // Should not be used.
    sprintf(this_box_volume.p_physics->name,"Mask");
    this_box_volume.geometry.is_mask_volume = 1;
    
    
// Read the material input, or if it lacks, use automatic linking.
} else if (material_string && strlen(material_string) && strcmp(material_string, "NULL") && strcmp(material_string, "0")) {
    // A geometry string was given, use it to determine which material
    if (0 == strcmp(material_string,"vacuum") || 0 == strcmp(material_string,"Vacuum")) {
        // One could have a global physics struct for vacuum instead of creating one for each
        this_box_volume.p_physics = malloc(sizeof(struct physics_struct));
        this_box_volume.p_physics->is_vacuum = 1; // Makes this volume a vacuum
        this_box_volume.p_physics->number_of_processes = (int) 0; // Should not be used.
        this_box_volume.p_physics->my_a = 0; // Should not be used.
        sprintf(this_box_volume.p_physics->name,"Vacuum");
    } else if (0 == strcmp(material_string,"exit") || 0 == strcmp(material_string,"Exit")) {
        // One could have a global physics struct for vacuum instead of creating one for each
        this_box_volume.p_physics = malloc(sizeof(struct physics_struct));
        this_box_volume.p_physics->is_vacuum = 1; // Makes this volume a vacuum
        this_box_volume.p_physics->number_of_processes = (int) 0; // Should not be used.
        this_box_volume.p_physics->my_a = 0; // Should not be used.
        this_box_volume.geometry.is_exit_volume = 1;
        sprintf(this_box_volume.p_physics->name,"Exit");
    } else {
        #ifndef MATERIAL_DETECTOR
            printf("Need to define a material before refering to it in a geometry %s.\n",NAME_CURRENT_COMP);
            exit(1);
        #endif
        for (loop_index=0;loop_index<global_material_list.num_elements;loop_index++) {
            if (0 == strcmp(material_string,global_material_list.elements[loop_index].name)) {
                this_box_volume.p_physics = global_material_list.elements[loop_index].physics;
                break;
            }
            if (loop_index == global_material_list.num_elements-1) {
            printf("\n");
            printf("ERROR: The material string \"%s\" in Union geometry \"%s\" did not match a specified material. \n",material_string,NAME_CURRENT_COMP);
            printf("       The materials available at this point (need to be defined before the geometry): \n");
            for (loop_index=0;loop_index<global_material_list.num_elements;loop_index++)
              printf("         %s\n",global_material_list.elements[loop_index].name);
            printf("\n");
            printf("       It is also possible to use one of the defualt materials avaiable: \n");
            printf("           Vacuum (for a Volume without scattering or absorption)\n");
            printf("           Exit (for a Volume where the ray exits the component if it enters)\n");
            printf("           Mask (for a Volume that masks existing volumes specified in the mask_string\n");
            exit(1);
          }
        }
    }
} else {
    // Automatic linking, simply using the last defined material.
    #ifndef MATERIAL_DETECTOR
        printf("Need to define a material before the geometry to use automatic linking %s.\n",NAME_CURRENT_COMP);
        exit(1);
    #endif
    this_box_volume.p_physics = global_material_list.elements[global_material_list.num_elements-1].physics;
}

sprintf(this_box_volume.name,NAME_CURRENT_COMP);
sprintf(this_box_volume.geometry.shape,"box");
this_box_volume.geometry.priority_value = priority;
this_box_volume.geometry.geometry_p_interact = p_interact;
// Currently the coordinates will be in absolute space.
this_box_volume.geometry.center = POS_A_CURRENT_COMP;

this_box_storage.z_depth = zdepth;
this_box_storage.x_width1 = xwidth;
this_box_storage.y_height1 = yheight;

this_box_storage.is_rectangle = 0;
if (xwidth2 < 0 && yheight2 < 0) this_box_storage.is_rectangle = 1;
if (xwidth == xwidth2 && yheight == yheight2) this_box_storage.is_rectangle = 1;

if (xwidth2 < 0) {
    this_box_storage.x_width2 = xwidth;
    xwidth2 = xwidth;
} else this_box_storage.x_width2 = xwidth2;

if (yheight2 < 0) {
    this_box_storage.y_height2 = yheight;
    yheight2 = yheight;
} else this_box_storage.y_height2 = yheight2;


this_box_storage.normal_vectors[0] = coords_set(0,0,1);
this_box_storage.normal_vectors[1] = coords_set(0,0,1);

// for sides with y component = 0
x_component = 2*zdepth/sqrt((xwidth-xwidth2)*(xwidth-xwidth2)+4*zdepth*zdepth);
z_component = (xwidth-xwidth2)/sqrt(4*zdepth*zdepth+(xwidth-xwidth2)*(xwidth-xwidth2));

this_box_storage.normal_vectors[2] = coords_set(x_component,0,z_component);
this_box_storage.normal_vectors[3] = coords_set(-x_component,0,z_component);

// for sides with x component = 0
y_component = 2*zdepth/sqrt((yheight-yheight2)*(yheight-yheight2)+4*zdepth*zdepth);
z_component = (yheight-yheight2)/sqrt(4*zdepth*zdepth+(yheight-yheight2)*(yheight-yheight2));

this_box_storage.normal_vectors[4] = coords_set(0,y_component,z_component);
this_box_storage.normal_vectors[5] = coords_set(0,-y_component,z_component);

this_box_volume.geometry.visualization_on = visualize;

this_box_volume.geometry.geometry_parameters.p_box_storage = &this_box_storage;

// Assign pointers to functions for intersection with the shape, checking if a point is inside the shape
if (this_box_storage.is_rectangle == 1) {
this_box_volume.geometry.intersect_function = &sample_box_intersect_simple;
this_box_volume.geometry.within_function = &r_within_box_simple;
}
else {
this_box_volume.geometry.intersect_function = &sample_box_intersect_advanced;
this_box_volume.geometry.within_function = &r_within_box_advanced;
}

this_box_volume.geometry.shell_points = &box_shell_points;
this_box_volume.geometry.mcdisplay_function = &mcdisplay_box_function;
this_box_volume.geometry.initialize_from_main_function = &initialize_box_geometry_from_main_component;
this_box_volume.geometry.process_rot_allocated = 0;

this_box_volume.geometry.copy_geometry_parameters = &allocate_box_storage_copy;

rot_copy(this_box_volume.geometry.rotation_matrix,ROT_A_CURRENT_COMP); //check how ROT_R_CURRENT_COMP would work
rot_transpose(ROT_A_CURRENT_COMP,this_box_volume.geometry.transpose_rotation_matrix);

// Initialize loggers
this_box_volume.loggers.num_elements = 0;

// packing the information into the global_geometry_element, which is then included in the global_geometry_list.
sprintf(global_geometry_element.name,NAME_CURRENT_COMP);
global_geometry_element.activation_counter = number_of_activations;
global_geometry_element.component_index = INDEX_CURRENT_COMP;
global_geometry_element.Volume = &this_box_volume; // Would be nicer if this m was a pointer, now we have the (small) data two places
add_element_to_geometry_list(&global_geometry_list,global_geometry_element);
}
#line 9081 "./Union_test_mask.c"
#undef number_of_activations
#undef mask_setting
#undef mask_string
#undef p_interact
#undef focus_r
#undef focus_xh
#undef focus_xw
#undef focus_ah
#undef focus_aw
#undef target_z
#undef target_y
#undef target_x
#undef target_index
#undef visualize
#undef yheight2
#undef xwidth2
#undef zdepth
#undef yheight
#undef xwidth
#undef priority
#undef material_string
#undef this_box_storage
#undef global_geometry_element
#undef this_box_volume
#undef loop_index
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component mask_1_2. */
  SIG_MESSAGE("mask_1_2 (Init)");
#define mccompcurname  mask_1_2
#define mccompcurtype  Union_box
#define mccompcurindex 9
#define loop_index mccmask_1_2_loop_index
#define this_box_volume mccmask_1_2_this_box_volume
#define global_geometry_element mccmask_1_2_global_geometry_element
#define this_box_storage mccmask_1_2_this_box_storage
#define material_string mccmask_1_2_material_string
#define priority mccmask_1_2_priority
#define xwidth mccmask_1_2_xwidth
#define yheight mccmask_1_2_yheight
#define zdepth mccmask_1_2_zdepth
#define xwidth2 mccmask_1_2_xwidth2
#define yheight2 mccmask_1_2_yheight2
#define visualize mccmask_1_2_visualize
#define target_index mccmask_1_2_target_index
#define target_x mccmask_1_2_target_x
#define target_y mccmask_1_2_target_y
#define target_z mccmask_1_2_target_z
#define focus_aw mccmask_1_2_focus_aw
#define focus_ah mccmask_1_2_focus_ah
#define focus_xw mccmask_1_2_focus_xw
#define focus_xh mccmask_1_2_focus_xh
#define focus_r mccmask_1_2_focus_r
#define p_interact mccmask_1_2_p_interact
#define mask_string mccmask_1_2_mask_string
#define mask_setting mccmask_1_2_mask_setting
#define number_of_activations mccmask_1_2_number_of_activations
#line 223 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../contrib/union/Union_box.comp"
{
// Initializes the focusing system for this volume including input sanitation.
focus_initialize(&this_box_volume.geometry, POS_A_COMP_INDEX(INDEX_CURRENT_COMP+target_index), POS_A_CURRENT_COMP, ROT_A_CURRENT_COMP, target_index, target_x, target_y, target_z, focus_aw, focus_ah, focus_xw, focus_xh, focus_r, NAME_CURRENT_COMP);

// Input sanitation for this geometry
if (xwidth <= 0) {
  printf("\nERROR in Union_box named %s, the xwidth is <= 0. \n",NAME_CURRENT_COMP);
  exit(1);
}
if (yheight <= 0) {
  printf("\nERROR in Union_box named %s, yheight is <= 0. \n",NAME_CURRENT_COMP);
  exit(1);
}
if (zdepth <= 0) {
  printf("\nERROR in Union_box named %s, zdepth is <= 0. \n",NAME_CURRENT_COMP);
  exit(1);
}
if (xwidth2 <= 0 && xwidth2 != -1) {
  printf("\nERROR in Union_box named %s, the xwidth2 is <= 0. \n",NAME_CURRENT_COMP);
  exit(1);
}
if (yheight2 <= 0 && yheight2 != -1) {
  printf("\nERROR in Union_box named %s, yheight2 is <= 0. \n",NAME_CURRENT_COMP);
  exit(1);
}

// Use sanitation
if (global_material_list.num_elements == 0) {
  printf("\nERROR: Need to define a material using Union_make_material before using a Union geometry component. \n");
  printf("       %s was defined before first use of Union_make_material.\n",NAME_CURRENT_COMP);
  exit(1);
}

this_box_volume.geometry.is_masked_volume = 0;
this_box_volume.geometry.is_exit_volume = 0;
this_box_volume.geometry.is_mask_volume = 0;
// check if the volume is a mask, if it is the material string is irelevant.
if (mask_string && strlen(mask_string) && strcmp(mask_string, "NULL") && strcmp(mask_string, "0")) {
    // A mask volume is used to limit the extend of other volumes, called the masked volumes. These are specified in the mask_string.
    // In order for a ray to enter a masked volume, it needs to be both in the region covered by that volume AND the mask volume.
    // When more than
    this_box_volume.geometry.mask_mode = 1; // Default is mask mode is ALL
    if (mask_setting && strlen(mask_setting) && strcmp(mask_setting, "NULL") && strcmp(mask_setting, "0")) {
        if (strcmp(mask_setting,"ALL") == 0 || strcmp(mask_setting,"All") == 0) this_box_volume.geometry.mask_mode = 1;
        else if (strcmp(mask_setting,"ANY") == 0 || strcmp(mask_setting,"Any") == 0) this_box_volume.geometry.mask_mode = 2;
        else {
            printf("The mask_mode of component %s is set to %s, but must be either ALL or ANY.\n",NAME_CURRENT_COMP,mask_setting);
            exit(1);
        }
    }
    
    for (loop_index=0;loop_index<global_geometry_list.num_elements;loop_index++) {
        // Add mask list
        if (1 == manual_linking_function(global_geometry_list.elements[loop_index].name,mask_string)) {
            add_element_to_int_list(&this_box_volume.geometry.mask_list,global_geometry_list.elements[loop_index].component_index);
            add_element_to_int_list(&global_geometry_list.elements[loop_index].Volume->geometry.masked_by_list,INDEX_CURRENT_COMP);
            global_geometry_list.elements[loop_index].Volume->geometry.is_masked_volume = 1;
            if (this_box_volume.geometry.mask_mode == 2)
                global_geometry_list.elements[loop_index].Volume->geometry.mask_mode = 2;
            if (this_box_volume.geometry.mask_mode == 1) {
                if (global_geometry_list.elements[loop_index].Volume->geometry.is_masked_volume == 1 && global_geometry_list.elements[loop_index].Volume->geometry.mask_mode != 2)
                    // If more than one mask is added to one volume, the ANY mode overwrites the (default) ALL mode.
                    global_geometry_list.elements[loop_index].Volume->geometry.mask_mode = 1;
            }
            
            found_geometries = 1;
        }
    }
    if (found_geometries == 0) {
        printf("The mask_string in geometry: %s did not find any of the specified volumes in the mask_string %s \n",NAME_CURRENT_COMP,mask_string);
        exit(1);
    }
    this_box_volume.p_physics = malloc(sizeof(struct physics_struct));
    this_box_volume.p_physics->is_vacuum = 0; // Makes this volume a vacuum
    this_box_volume.p_physics->number_of_processes = (int) 0; // Should not be used.
    this_box_volume.p_physics->my_a = 0; // Should not be used.
    sprintf(this_box_volume.p_physics->name,"Mask");
    this_box_volume.geometry.is_mask_volume = 1;
    
    
// Read the material input, or if it lacks, use automatic linking.
} else if (material_string && strlen(material_string) && strcmp(material_string, "NULL") && strcmp(material_string, "0")) {
    // A geometry string was given, use it to determine which material
    if (0 == strcmp(material_string,"vacuum") || 0 == strcmp(material_string,"Vacuum")) {
        // One could have a global physics struct for vacuum instead of creating one for each
        this_box_volume.p_physics = malloc(sizeof(struct physics_struct));
        this_box_volume.p_physics->is_vacuum = 1; // Makes this volume a vacuum
        this_box_volume.p_physics->number_of_processes = (int) 0; // Should not be used.
        this_box_volume.p_physics->my_a = 0; // Should not be used.
        sprintf(this_box_volume.p_physics->name,"Vacuum");
    } else if (0 == strcmp(material_string,"exit") || 0 == strcmp(material_string,"Exit")) {
        // One could have a global physics struct for vacuum instead of creating one for each
        this_box_volume.p_physics = malloc(sizeof(struct physics_struct));
        this_box_volume.p_physics->is_vacuum = 1; // Makes this volume a vacuum
        this_box_volume.p_physics->number_of_processes = (int) 0; // Should not be used.
        this_box_volume.p_physics->my_a = 0; // Should not be used.
        this_box_volume.geometry.is_exit_volume = 1;
        sprintf(this_box_volume.p_physics->name,"Exit");
    } else {
        #ifndef MATERIAL_DETECTOR
            printf("Need to define a material before refering to it in a geometry %s.\n",NAME_CURRENT_COMP);
            exit(1);
        #endif
        for (loop_index=0;loop_index<global_material_list.num_elements;loop_index++) {
            if (0 == strcmp(material_string,global_material_list.elements[loop_index].name)) {
                this_box_volume.p_physics = global_material_list.elements[loop_index].physics;
                break;
            }
            if (loop_index == global_material_list.num_elements-1) {
            printf("\n");
            printf("ERROR: The material string \"%s\" in Union geometry \"%s\" did not match a specified material. \n",material_string,NAME_CURRENT_COMP);
            printf("       The materials available at this point (need to be defined before the geometry): \n");
            for (loop_index=0;loop_index<global_material_list.num_elements;loop_index++)
              printf("         %s\n",global_material_list.elements[loop_index].name);
            printf("\n");
            printf("       It is also possible to use one of the defualt materials avaiable: \n");
            printf("           Vacuum (for a Volume without scattering or absorption)\n");
            printf("           Exit (for a Volume where the ray exits the component if it enters)\n");
            printf("           Mask (for a Volume that masks existing volumes specified in the mask_string\n");
            exit(1);
          }
        }
    }
} else {
    // Automatic linking, simply using the last defined material.
    #ifndef MATERIAL_DETECTOR
        printf("Need to define a material before the geometry to use automatic linking %s.\n",NAME_CURRENT_COMP);
        exit(1);
    #endif
    this_box_volume.p_physics = global_material_list.elements[global_material_list.num_elements-1].physics;
}

sprintf(this_box_volume.name,NAME_CURRENT_COMP);
sprintf(this_box_volume.geometry.shape,"box");
this_box_volume.geometry.priority_value = priority;
this_box_volume.geometry.geometry_p_interact = p_interact;
// Currently the coordinates will be in absolute space.
this_box_volume.geometry.center = POS_A_CURRENT_COMP;

this_box_storage.z_depth = zdepth;
this_box_storage.x_width1 = xwidth;
this_box_storage.y_height1 = yheight;

this_box_storage.is_rectangle = 0;
if (xwidth2 < 0 && yheight2 < 0) this_box_storage.is_rectangle = 1;
if (xwidth == xwidth2 && yheight == yheight2) this_box_storage.is_rectangle = 1;

if (xwidth2 < 0) {
    this_box_storage.x_width2 = xwidth;
    xwidth2 = xwidth;
} else this_box_storage.x_width2 = xwidth2;

if (yheight2 < 0) {
    this_box_storage.y_height2 = yheight;
    yheight2 = yheight;
} else this_box_storage.y_height2 = yheight2;


this_box_storage.normal_vectors[0] = coords_set(0,0,1);
this_box_storage.normal_vectors[1] = coords_set(0,0,1);

// for sides with y component = 0
x_component = 2*zdepth/sqrt((xwidth-xwidth2)*(xwidth-xwidth2)+4*zdepth*zdepth);
z_component = (xwidth-xwidth2)/sqrt(4*zdepth*zdepth+(xwidth-xwidth2)*(xwidth-xwidth2));

this_box_storage.normal_vectors[2] = coords_set(x_component,0,z_component);
this_box_storage.normal_vectors[3] = coords_set(-x_component,0,z_component);

// for sides with x component = 0
y_component = 2*zdepth/sqrt((yheight-yheight2)*(yheight-yheight2)+4*zdepth*zdepth);
z_component = (yheight-yheight2)/sqrt(4*zdepth*zdepth+(yheight-yheight2)*(yheight-yheight2));

this_box_storage.normal_vectors[4] = coords_set(0,y_component,z_component);
this_box_storage.normal_vectors[5] = coords_set(0,-y_component,z_component);

this_box_volume.geometry.visualization_on = visualize;

this_box_volume.geometry.geometry_parameters.p_box_storage = &this_box_storage;

// Assign pointers to functions for intersection with the shape, checking if a point is inside the shape
if (this_box_storage.is_rectangle == 1) {
this_box_volume.geometry.intersect_function = &sample_box_intersect_simple;
this_box_volume.geometry.within_function = &r_within_box_simple;
}
else {
this_box_volume.geometry.intersect_function = &sample_box_intersect_advanced;
this_box_volume.geometry.within_function = &r_within_box_advanced;
}

this_box_volume.geometry.shell_points = &box_shell_points;
this_box_volume.geometry.mcdisplay_function = &mcdisplay_box_function;
this_box_volume.geometry.initialize_from_main_function = &initialize_box_geometry_from_main_component;
this_box_volume.geometry.process_rot_allocated = 0;

this_box_volume.geometry.copy_geometry_parameters = &allocate_box_storage_copy;

rot_copy(this_box_volume.geometry.rotation_matrix,ROT_A_CURRENT_COMP); //check how ROT_R_CURRENT_COMP would work
rot_transpose(ROT_A_CURRENT_COMP,this_box_volume.geometry.transpose_rotation_matrix);

// Initialize loggers
this_box_volume.loggers.num_elements = 0;

// packing the information into the global_geometry_element, which is then included in the global_geometry_list.
sprintf(global_geometry_element.name,NAME_CURRENT_COMP);
global_geometry_element.activation_counter = number_of_activations;
global_geometry_element.component_index = INDEX_CURRENT_COMP;
global_geometry_element.Volume = &this_box_volume; // Would be nicer if this m was a pointer, now we have the (small) data two places
add_element_to_geometry_list(&global_geometry_list,global_geometry_element);
}
#line 9351 "./Union_test_mask.c"
#undef number_of_activations
#undef mask_setting
#undef mask_string
#undef p_interact
#undef focus_r
#undef focus_xh
#undef focus_xw
#undef focus_ah
#undef focus_aw
#undef target_z
#undef target_y
#undef target_x
#undef target_index
#undef visualize
#undef yheight2
#undef xwidth2
#undef zdepth
#undef yheight
#undef xwidth
#undef priority
#undef material_string
#undef this_box_storage
#undef global_geometry_element
#undef this_box_volume
#undef loop_index
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component mask_2_1. */
  SIG_MESSAGE("mask_2_1 (Init)");
#define mccompcurname  mask_2_1
#define mccompcurtype  Union_box
#define mccompcurindex 10
#define loop_index mccmask_2_1_loop_index
#define this_box_volume mccmask_2_1_this_box_volume
#define global_geometry_element mccmask_2_1_global_geometry_element
#define this_box_storage mccmask_2_1_this_box_storage
#define material_string mccmask_2_1_material_string
#define priority mccmask_2_1_priority
#define xwidth mccmask_2_1_xwidth
#define yheight mccmask_2_1_yheight
#define zdepth mccmask_2_1_zdepth
#define xwidth2 mccmask_2_1_xwidth2
#define yheight2 mccmask_2_1_yheight2
#define visualize mccmask_2_1_visualize
#define target_index mccmask_2_1_target_index
#define target_x mccmask_2_1_target_x
#define target_y mccmask_2_1_target_y
#define target_z mccmask_2_1_target_z
#define focus_aw mccmask_2_1_focus_aw
#define focus_ah mccmask_2_1_focus_ah
#define focus_xw mccmask_2_1_focus_xw
#define focus_xh mccmask_2_1_focus_xh
#define focus_r mccmask_2_1_focus_r
#define p_interact mccmask_2_1_p_interact
#define mask_string mccmask_2_1_mask_string
#define mask_setting mccmask_2_1_mask_setting
#define number_of_activations mccmask_2_1_number_of_activations
#line 223 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../contrib/union/Union_box.comp"
{
// Initializes the focusing system for this volume including input sanitation.
focus_initialize(&this_box_volume.geometry, POS_A_COMP_INDEX(INDEX_CURRENT_COMP+target_index), POS_A_CURRENT_COMP, ROT_A_CURRENT_COMP, target_index, target_x, target_y, target_z, focus_aw, focus_ah, focus_xw, focus_xh, focus_r, NAME_CURRENT_COMP);

// Input sanitation for this geometry
if (xwidth <= 0) {
  printf("\nERROR in Union_box named %s, the xwidth is <= 0. \n",NAME_CURRENT_COMP);
  exit(1);
}
if (yheight <= 0) {
  printf("\nERROR in Union_box named %s, yheight is <= 0. \n",NAME_CURRENT_COMP);
  exit(1);
}
if (zdepth <= 0) {
  printf("\nERROR in Union_box named %s, zdepth is <= 0. \n",NAME_CURRENT_COMP);
  exit(1);
}
if (xwidth2 <= 0 && xwidth2 != -1) {
  printf("\nERROR in Union_box named %s, the xwidth2 is <= 0. \n",NAME_CURRENT_COMP);
  exit(1);
}
if (yheight2 <= 0 && yheight2 != -1) {
  printf("\nERROR in Union_box named %s, yheight2 is <= 0. \n",NAME_CURRENT_COMP);
  exit(1);
}

// Use sanitation
if (global_material_list.num_elements == 0) {
  printf("\nERROR: Need to define a material using Union_make_material before using a Union geometry component. \n");
  printf("       %s was defined before first use of Union_make_material.\n",NAME_CURRENT_COMP);
  exit(1);
}

this_box_volume.geometry.is_masked_volume = 0;
this_box_volume.geometry.is_exit_volume = 0;
this_box_volume.geometry.is_mask_volume = 0;
// check if the volume is a mask, if it is the material string is irelevant.
if (mask_string && strlen(mask_string) && strcmp(mask_string, "NULL") && strcmp(mask_string, "0")) {
    // A mask volume is used to limit the extend of other volumes, called the masked volumes. These are specified in the mask_string.
    // In order for a ray to enter a masked volume, it needs to be both in the region covered by that volume AND the mask volume.
    // When more than
    this_box_volume.geometry.mask_mode = 1; // Default is mask mode is ALL
    if (mask_setting && strlen(mask_setting) && strcmp(mask_setting, "NULL") && strcmp(mask_setting, "0")) {
        if (strcmp(mask_setting,"ALL") == 0 || strcmp(mask_setting,"All") == 0) this_box_volume.geometry.mask_mode = 1;
        else if (strcmp(mask_setting,"ANY") == 0 || strcmp(mask_setting,"Any") == 0) this_box_volume.geometry.mask_mode = 2;
        else {
            printf("The mask_mode of component %s is set to %s, but must be either ALL or ANY.\n",NAME_CURRENT_COMP,mask_setting);
            exit(1);
        }
    }
    
    for (loop_index=0;loop_index<global_geometry_list.num_elements;loop_index++) {
        // Add mask list
        if (1 == manual_linking_function(global_geometry_list.elements[loop_index].name,mask_string)) {
            add_element_to_int_list(&this_box_volume.geometry.mask_list,global_geometry_list.elements[loop_index].component_index);
            add_element_to_int_list(&global_geometry_list.elements[loop_index].Volume->geometry.masked_by_list,INDEX_CURRENT_COMP);
            global_geometry_list.elements[loop_index].Volume->geometry.is_masked_volume = 1;
            if (this_box_volume.geometry.mask_mode == 2)
                global_geometry_list.elements[loop_index].Volume->geometry.mask_mode = 2;
            if (this_box_volume.geometry.mask_mode == 1) {
                if (global_geometry_list.elements[loop_index].Volume->geometry.is_masked_volume == 1 && global_geometry_list.elements[loop_index].Volume->geometry.mask_mode != 2)
                    // If more than one mask is added to one volume, the ANY mode overwrites the (default) ALL mode.
                    global_geometry_list.elements[loop_index].Volume->geometry.mask_mode = 1;
            }
            
            found_geometries = 1;
        }
    }
    if (found_geometries == 0) {
        printf("The mask_string in geometry: %s did not find any of the specified volumes in the mask_string %s \n",NAME_CURRENT_COMP,mask_string);
        exit(1);
    }
    this_box_volume.p_physics = malloc(sizeof(struct physics_struct));
    this_box_volume.p_physics->is_vacuum = 0; // Makes this volume a vacuum
    this_box_volume.p_physics->number_of_processes = (int) 0; // Should not be used.
    this_box_volume.p_physics->my_a = 0; // Should not be used.
    sprintf(this_box_volume.p_physics->name,"Mask");
    this_box_volume.geometry.is_mask_volume = 1;
    
    
// Read the material input, or if it lacks, use automatic linking.
} else if (material_string && strlen(material_string) && strcmp(material_string, "NULL") && strcmp(material_string, "0")) {
    // A geometry string was given, use it to determine which material
    if (0 == strcmp(material_string,"vacuum") || 0 == strcmp(material_string,"Vacuum")) {
        // One could have a global physics struct for vacuum instead of creating one for each
        this_box_volume.p_physics = malloc(sizeof(struct physics_struct));
        this_box_volume.p_physics->is_vacuum = 1; // Makes this volume a vacuum
        this_box_volume.p_physics->number_of_processes = (int) 0; // Should not be used.
        this_box_volume.p_physics->my_a = 0; // Should not be used.
        sprintf(this_box_volume.p_physics->name,"Vacuum");
    } else if (0 == strcmp(material_string,"exit") || 0 == strcmp(material_string,"Exit")) {
        // One could have a global physics struct for vacuum instead of creating one for each
        this_box_volume.p_physics = malloc(sizeof(struct physics_struct));
        this_box_volume.p_physics->is_vacuum = 1; // Makes this volume a vacuum
        this_box_volume.p_physics->number_of_processes = (int) 0; // Should not be used.
        this_box_volume.p_physics->my_a = 0; // Should not be used.
        this_box_volume.geometry.is_exit_volume = 1;
        sprintf(this_box_volume.p_physics->name,"Exit");
    } else {
        #ifndef MATERIAL_DETECTOR
            printf("Need to define a material before refering to it in a geometry %s.\n",NAME_CURRENT_COMP);
            exit(1);
        #endif
        for (loop_index=0;loop_index<global_material_list.num_elements;loop_index++) {
            if (0 == strcmp(material_string,global_material_list.elements[loop_index].name)) {
                this_box_volume.p_physics = global_material_list.elements[loop_index].physics;
                break;
            }
            if (loop_index == global_material_list.num_elements-1) {
            printf("\n");
            printf("ERROR: The material string \"%s\" in Union geometry \"%s\" did not match a specified material. \n",material_string,NAME_CURRENT_COMP);
            printf("       The materials available at this point (need to be defined before the geometry): \n");
            for (loop_index=0;loop_index<global_material_list.num_elements;loop_index++)
              printf("         %s\n",global_material_list.elements[loop_index].name);
            printf("\n");
            printf("       It is also possible to use one of the defualt materials avaiable: \n");
            printf("           Vacuum (for a Volume without scattering or absorption)\n");
            printf("           Exit (for a Volume where the ray exits the component if it enters)\n");
            printf("           Mask (for a Volume that masks existing volumes specified in the mask_string\n");
            exit(1);
          }
        }
    }
} else {
    // Automatic linking, simply using the last defined material.
    #ifndef MATERIAL_DETECTOR
        printf("Need to define a material before the geometry to use automatic linking %s.\n",NAME_CURRENT_COMP);
        exit(1);
    #endif
    this_box_volume.p_physics = global_material_list.elements[global_material_list.num_elements-1].physics;
}

sprintf(this_box_volume.name,NAME_CURRENT_COMP);
sprintf(this_box_volume.geometry.shape,"box");
this_box_volume.geometry.priority_value = priority;
this_box_volume.geometry.geometry_p_interact = p_interact;
// Currently the coordinates will be in absolute space.
this_box_volume.geometry.center = POS_A_CURRENT_COMP;

this_box_storage.z_depth = zdepth;
this_box_storage.x_width1 = xwidth;
this_box_storage.y_height1 = yheight;

this_box_storage.is_rectangle = 0;
if (xwidth2 < 0 && yheight2 < 0) this_box_storage.is_rectangle = 1;
if (xwidth == xwidth2 && yheight == yheight2) this_box_storage.is_rectangle = 1;

if (xwidth2 < 0) {
    this_box_storage.x_width2 = xwidth;
    xwidth2 = xwidth;
} else this_box_storage.x_width2 = xwidth2;

if (yheight2 < 0) {
    this_box_storage.y_height2 = yheight;
    yheight2 = yheight;
} else this_box_storage.y_height2 = yheight2;


this_box_storage.normal_vectors[0] = coords_set(0,0,1);
this_box_storage.normal_vectors[1] = coords_set(0,0,1);

// for sides with y component = 0
x_component = 2*zdepth/sqrt((xwidth-xwidth2)*(xwidth-xwidth2)+4*zdepth*zdepth);
z_component = (xwidth-xwidth2)/sqrt(4*zdepth*zdepth+(xwidth-xwidth2)*(xwidth-xwidth2));

this_box_storage.normal_vectors[2] = coords_set(x_component,0,z_component);
this_box_storage.normal_vectors[3] = coords_set(-x_component,0,z_component);

// for sides with x component = 0
y_component = 2*zdepth/sqrt((yheight-yheight2)*(yheight-yheight2)+4*zdepth*zdepth);
z_component = (yheight-yheight2)/sqrt(4*zdepth*zdepth+(yheight-yheight2)*(yheight-yheight2));

this_box_storage.normal_vectors[4] = coords_set(0,y_component,z_component);
this_box_storage.normal_vectors[5] = coords_set(0,-y_component,z_component);

this_box_volume.geometry.visualization_on = visualize;

this_box_volume.geometry.geometry_parameters.p_box_storage = &this_box_storage;

// Assign pointers to functions for intersection with the shape, checking if a point is inside the shape
if (this_box_storage.is_rectangle == 1) {
this_box_volume.geometry.intersect_function = &sample_box_intersect_simple;
this_box_volume.geometry.within_function = &r_within_box_simple;
}
else {
this_box_volume.geometry.intersect_function = &sample_box_intersect_advanced;
this_box_volume.geometry.within_function = &r_within_box_advanced;
}

this_box_volume.geometry.shell_points = &box_shell_points;
this_box_volume.geometry.mcdisplay_function = &mcdisplay_box_function;
this_box_volume.geometry.initialize_from_main_function = &initialize_box_geometry_from_main_component;
this_box_volume.geometry.process_rot_allocated = 0;

this_box_volume.geometry.copy_geometry_parameters = &allocate_box_storage_copy;

rot_copy(this_box_volume.geometry.rotation_matrix,ROT_A_CURRENT_COMP); //check how ROT_R_CURRENT_COMP would work
rot_transpose(ROT_A_CURRENT_COMP,this_box_volume.geometry.transpose_rotation_matrix);

// Initialize loggers
this_box_volume.loggers.num_elements = 0;

// packing the information into the global_geometry_element, which is then included in the global_geometry_list.
sprintf(global_geometry_element.name,NAME_CURRENT_COMP);
global_geometry_element.activation_counter = number_of_activations;
global_geometry_element.component_index = INDEX_CURRENT_COMP;
global_geometry_element.Volume = &this_box_volume; // Would be nicer if this m was a pointer, now we have the (small) data two places
add_element_to_geometry_list(&global_geometry_list,global_geometry_element);
}
#line 9621 "./Union_test_mask.c"
#undef number_of_activations
#undef mask_setting
#undef mask_string
#undef p_interact
#undef focus_r
#undef focus_xh
#undef focus_xw
#undef focus_ah
#undef focus_aw
#undef target_z
#undef target_y
#undef target_x
#undef target_index
#undef visualize
#undef yheight2
#undef xwidth2
#undef zdepth
#undef yheight
#undef xwidth
#undef priority
#undef material_string
#undef this_box_storage
#undef global_geometry_element
#undef this_box_volume
#undef loop_index
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component mask_2_2. */
  SIG_MESSAGE("mask_2_2 (Init)");
#define mccompcurname  mask_2_2
#define mccompcurtype  Union_box
#define mccompcurindex 11
#define loop_index mccmask_2_2_loop_index
#define this_box_volume mccmask_2_2_this_box_volume
#define global_geometry_element mccmask_2_2_global_geometry_element
#define this_box_storage mccmask_2_2_this_box_storage
#define material_string mccmask_2_2_material_string
#define priority mccmask_2_2_priority
#define xwidth mccmask_2_2_xwidth
#define yheight mccmask_2_2_yheight
#define zdepth mccmask_2_2_zdepth
#define xwidth2 mccmask_2_2_xwidth2
#define yheight2 mccmask_2_2_yheight2
#define visualize mccmask_2_2_visualize
#define target_index mccmask_2_2_target_index
#define target_x mccmask_2_2_target_x
#define target_y mccmask_2_2_target_y
#define target_z mccmask_2_2_target_z
#define focus_aw mccmask_2_2_focus_aw
#define focus_ah mccmask_2_2_focus_ah
#define focus_xw mccmask_2_2_focus_xw
#define focus_xh mccmask_2_2_focus_xh
#define focus_r mccmask_2_2_focus_r
#define p_interact mccmask_2_2_p_interact
#define mask_string mccmask_2_2_mask_string
#define mask_setting mccmask_2_2_mask_setting
#define number_of_activations mccmask_2_2_number_of_activations
#line 223 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../contrib/union/Union_box.comp"
{
// Initializes the focusing system for this volume including input sanitation.
focus_initialize(&this_box_volume.geometry, POS_A_COMP_INDEX(INDEX_CURRENT_COMP+target_index), POS_A_CURRENT_COMP, ROT_A_CURRENT_COMP, target_index, target_x, target_y, target_z, focus_aw, focus_ah, focus_xw, focus_xh, focus_r, NAME_CURRENT_COMP);

// Input sanitation for this geometry
if (xwidth <= 0) {
  printf("\nERROR in Union_box named %s, the xwidth is <= 0. \n",NAME_CURRENT_COMP);
  exit(1);
}
if (yheight <= 0) {
  printf("\nERROR in Union_box named %s, yheight is <= 0. \n",NAME_CURRENT_COMP);
  exit(1);
}
if (zdepth <= 0) {
  printf("\nERROR in Union_box named %s, zdepth is <= 0. \n",NAME_CURRENT_COMP);
  exit(1);
}
if (xwidth2 <= 0 && xwidth2 != -1) {
  printf("\nERROR in Union_box named %s, the xwidth2 is <= 0. \n",NAME_CURRENT_COMP);
  exit(1);
}
if (yheight2 <= 0 && yheight2 != -1) {
  printf("\nERROR in Union_box named %s, yheight2 is <= 0. \n",NAME_CURRENT_COMP);
  exit(1);
}

// Use sanitation
if (global_material_list.num_elements == 0) {
  printf("\nERROR: Need to define a material using Union_make_material before using a Union geometry component. \n");
  printf("       %s was defined before first use of Union_make_material.\n",NAME_CURRENT_COMP);
  exit(1);
}

this_box_volume.geometry.is_masked_volume = 0;
this_box_volume.geometry.is_exit_volume = 0;
this_box_volume.geometry.is_mask_volume = 0;
// check if the volume is a mask, if it is the material string is irelevant.
if (mask_string && strlen(mask_string) && strcmp(mask_string, "NULL") && strcmp(mask_string, "0")) {
    // A mask volume is used to limit the extend of other volumes, called the masked volumes. These are specified in the mask_string.
    // In order for a ray to enter a masked volume, it needs to be both in the region covered by that volume AND the mask volume.
    // When more than
    this_box_volume.geometry.mask_mode = 1; // Default is mask mode is ALL
    if (mask_setting && strlen(mask_setting) && strcmp(mask_setting, "NULL") && strcmp(mask_setting, "0")) {
        if (strcmp(mask_setting,"ALL") == 0 || strcmp(mask_setting,"All") == 0) this_box_volume.geometry.mask_mode = 1;
        else if (strcmp(mask_setting,"ANY") == 0 || strcmp(mask_setting,"Any") == 0) this_box_volume.geometry.mask_mode = 2;
        else {
            printf("The mask_mode of component %s is set to %s, but must be either ALL or ANY.\n",NAME_CURRENT_COMP,mask_setting);
            exit(1);
        }
    }
    
    for (loop_index=0;loop_index<global_geometry_list.num_elements;loop_index++) {
        // Add mask list
        if (1 == manual_linking_function(global_geometry_list.elements[loop_index].name,mask_string)) {
            add_element_to_int_list(&this_box_volume.geometry.mask_list,global_geometry_list.elements[loop_index].component_index);
            add_element_to_int_list(&global_geometry_list.elements[loop_index].Volume->geometry.masked_by_list,INDEX_CURRENT_COMP);
            global_geometry_list.elements[loop_index].Volume->geometry.is_masked_volume = 1;
            if (this_box_volume.geometry.mask_mode == 2)
                global_geometry_list.elements[loop_index].Volume->geometry.mask_mode = 2;
            if (this_box_volume.geometry.mask_mode == 1) {
                if (global_geometry_list.elements[loop_index].Volume->geometry.is_masked_volume == 1 && global_geometry_list.elements[loop_index].Volume->geometry.mask_mode != 2)
                    // If more than one mask is added to one volume, the ANY mode overwrites the (default) ALL mode.
                    global_geometry_list.elements[loop_index].Volume->geometry.mask_mode = 1;
            }
            
            found_geometries = 1;
        }
    }
    if (found_geometries == 0) {
        printf("The mask_string in geometry: %s did not find any of the specified volumes in the mask_string %s \n",NAME_CURRENT_COMP,mask_string);
        exit(1);
    }
    this_box_volume.p_physics = malloc(sizeof(struct physics_struct));
    this_box_volume.p_physics->is_vacuum = 0; // Makes this volume a vacuum
    this_box_volume.p_physics->number_of_processes = (int) 0; // Should not be used.
    this_box_volume.p_physics->my_a = 0; // Should not be used.
    sprintf(this_box_volume.p_physics->name,"Mask");
    this_box_volume.geometry.is_mask_volume = 1;
    
    
// Read the material input, or if it lacks, use automatic linking.
} else if (material_string && strlen(material_string) && strcmp(material_string, "NULL") && strcmp(material_string, "0")) {
    // A geometry string was given, use it to determine which material
    if (0 == strcmp(material_string,"vacuum") || 0 == strcmp(material_string,"Vacuum")) {
        // One could have a global physics struct for vacuum instead of creating one for each
        this_box_volume.p_physics = malloc(sizeof(struct physics_struct));
        this_box_volume.p_physics->is_vacuum = 1; // Makes this volume a vacuum
        this_box_volume.p_physics->number_of_processes = (int) 0; // Should not be used.
        this_box_volume.p_physics->my_a = 0; // Should not be used.
        sprintf(this_box_volume.p_physics->name,"Vacuum");
    } else if (0 == strcmp(material_string,"exit") || 0 == strcmp(material_string,"Exit")) {
        // One could have a global physics struct for vacuum instead of creating one for each
        this_box_volume.p_physics = malloc(sizeof(struct physics_struct));
        this_box_volume.p_physics->is_vacuum = 1; // Makes this volume a vacuum
        this_box_volume.p_physics->number_of_processes = (int) 0; // Should not be used.
        this_box_volume.p_physics->my_a = 0; // Should not be used.
        this_box_volume.geometry.is_exit_volume = 1;
        sprintf(this_box_volume.p_physics->name,"Exit");
    } else {
        #ifndef MATERIAL_DETECTOR
            printf("Need to define a material before refering to it in a geometry %s.\n",NAME_CURRENT_COMP);
            exit(1);
        #endif
        for (loop_index=0;loop_index<global_material_list.num_elements;loop_index++) {
            if (0 == strcmp(material_string,global_material_list.elements[loop_index].name)) {
                this_box_volume.p_physics = global_material_list.elements[loop_index].physics;
                break;
            }
            if (loop_index == global_material_list.num_elements-1) {
            printf("\n");
            printf("ERROR: The material string \"%s\" in Union geometry \"%s\" did not match a specified material. \n",material_string,NAME_CURRENT_COMP);
            printf("       The materials available at this point (need to be defined before the geometry): \n");
            for (loop_index=0;loop_index<global_material_list.num_elements;loop_index++)
              printf("         %s\n",global_material_list.elements[loop_index].name);
            printf("\n");
            printf("       It is also possible to use one of the defualt materials avaiable: \n");
            printf("           Vacuum (for a Volume without scattering or absorption)\n");
            printf("           Exit (for a Volume where the ray exits the component if it enters)\n");
            printf("           Mask (for a Volume that masks existing volumes specified in the mask_string\n");
            exit(1);
          }
        }
    }
} else {
    // Automatic linking, simply using the last defined material.
    #ifndef MATERIAL_DETECTOR
        printf("Need to define a material before the geometry to use automatic linking %s.\n",NAME_CURRENT_COMP);
        exit(1);
    #endif
    this_box_volume.p_physics = global_material_list.elements[global_material_list.num_elements-1].physics;
}

sprintf(this_box_volume.name,NAME_CURRENT_COMP);
sprintf(this_box_volume.geometry.shape,"box");
this_box_volume.geometry.priority_value = priority;
this_box_volume.geometry.geometry_p_interact = p_interact;
// Currently the coordinates will be in absolute space.
this_box_volume.geometry.center = POS_A_CURRENT_COMP;

this_box_storage.z_depth = zdepth;
this_box_storage.x_width1 = xwidth;
this_box_storage.y_height1 = yheight;

this_box_storage.is_rectangle = 0;
if (xwidth2 < 0 && yheight2 < 0) this_box_storage.is_rectangle = 1;
if (xwidth == xwidth2 && yheight == yheight2) this_box_storage.is_rectangle = 1;

if (xwidth2 < 0) {
    this_box_storage.x_width2 = xwidth;
    xwidth2 = xwidth;
} else this_box_storage.x_width2 = xwidth2;

if (yheight2 < 0) {
    this_box_storage.y_height2 = yheight;
    yheight2 = yheight;
} else this_box_storage.y_height2 = yheight2;


this_box_storage.normal_vectors[0] = coords_set(0,0,1);
this_box_storage.normal_vectors[1] = coords_set(0,0,1);

// for sides with y component = 0
x_component = 2*zdepth/sqrt((xwidth-xwidth2)*(xwidth-xwidth2)+4*zdepth*zdepth);
z_component = (xwidth-xwidth2)/sqrt(4*zdepth*zdepth+(xwidth-xwidth2)*(xwidth-xwidth2));

this_box_storage.normal_vectors[2] = coords_set(x_component,0,z_component);
this_box_storage.normal_vectors[3] = coords_set(-x_component,0,z_component);

// for sides with x component = 0
y_component = 2*zdepth/sqrt((yheight-yheight2)*(yheight-yheight2)+4*zdepth*zdepth);
z_component = (yheight-yheight2)/sqrt(4*zdepth*zdepth+(yheight-yheight2)*(yheight-yheight2));

this_box_storage.normal_vectors[4] = coords_set(0,y_component,z_component);
this_box_storage.normal_vectors[5] = coords_set(0,-y_component,z_component);

this_box_volume.geometry.visualization_on = visualize;

this_box_volume.geometry.geometry_parameters.p_box_storage = &this_box_storage;

// Assign pointers to functions for intersection with the shape, checking if a point is inside the shape
if (this_box_storage.is_rectangle == 1) {
this_box_volume.geometry.intersect_function = &sample_box_intersect_simple;
this_box_volume.geometry.within_function = &r_within_box_simple;
}
else {
this_box_volume.geometry.intersect_function = &sample_box_intersect_advanced;
this_box_volume.geometry.within_function = &r_within_box_advanced;
}

this_box_volume.geometry.shell_points = &box_shell_points;
this_box_volume.geometry.mcdisplay_function = &mcdisplay_box_function;
this_box_volume.geometry.initialize_from_main_function = &initialize_box_geometry_from_main_component;
this_box_volume.geometry.process_rot_allocated = 0;

this_box_volume.geometry.copy_geometry_parameters = &allocate_box_storage_copy;

rot_copy(this_box_volume.geometry.rotation_matrix,ROT_A_CURRENT_COMP); //check how ROT_R_CURRENT_COMP would work
rot_transpose(ROT_A_CURRENT_COMP,this_box_volume.geometry.transpose_rotation_matrix);

// Initialize loggers
this_box_volume.loggers.num_elements = 0;

// packing the information into the global_geometry_element, which is then included in the global_geometry_list.
sprintf(global_geometry_element.name,NAME_CURRENT_COMP);
global_geometry_element.activation_counter = number_of_activations;
global_geometry_element.component_index = INDEX_CURRENT_COMP;
global_geometry_element.Volume = &this_box_volume; // Would be nicer if this m was a pointer, now we have the (small) data two places
add_element_to_geometry_list(&global_geometry_list,global_geometry_element);
}
#line 9891 "./Union_test_mask.c"
#undef number_of_activations
#undef mask_setting
#undef mask_string
#undef p_interact
#undef focus_r
#undef focus_xh
#undef focus_xw
#undef focus_ah
#undef focus_aw
#undef target_z
#undef target_y
#undef target_x
#undef target_index
#undef visualize
#undef yheight2
#undef xwidth2
#undef zdepth
#undef yheight
#undef xwidth
#undef priority
#undef material_string
#undef this_box_storage
#undef global_geometry_element
#undef this_box_volume
#undef loop_index
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component test_sample. */
  SIG_MESSAGE("test_sample (Init)");
#define mccompcurname  test_sample
#define mccompcurtype  Union_master
#define mccompcurindex 12
#define verbal mcctest_sample_verbal
#define list_verbal mcctest_sample_list_verbal
#define trace_verbal mcctest_sample_trace_verbal
#define finally_verbal mcctest_sample_finally_verbal
#define starting_volume_warning mcctest_sample_starting_volume_warning
#define global_master_element mcctest_sample_global_master_element
#define this_global_master_index mcctest_sample_this_global_master_index
#define previous_master_index mcctest_sample_previous_master_index
#define geometry_list_index mcctest_sample_geometry_list_index
#define intersection_time_table mcctest_sample_intersection_time_table
#define Volumes mcctest_sample_Volumes
#define Geometries mcctest_sample_Geometries
#define starting_lists mcctest_sample_starting_lists
#define r mcctest_sample_r
#define r_start mcctest_sample_r_start
#define v mcctest_sample_v
#define error_msg mcctest_sample_error_msg
#define component_error_msg mcctest_sample_component_error_msg
#define string_output mcctest_sample_string_output
#define number_of_volumes mcctest_sample_number_of_volumes
#define volume_index mcctest_sample_volume_index
#define process_index mcctest_sample_process_index
#define solutions mcctest_sample_solutions
#define max_number_of_processes mcctest_sample_max_number_of_processes
#define limit mcctest_sample_limit
#define solution mcctest_sample_solution
#define min_solution mcctest_sample_min_solution
#define min_volume mcctest_sample_min_volume
#define time_found mcctest_sample_time_found
#define intersection_time mcctest_sample_intersection_time
#define min_intersection_time mcctest_sample_min_intersection_time
#define process mcctest_sample_process
#define process_start mcctest_sample_process_start
#define my_trace mcctest_sample_my_trace
#define p_my_trace mcctest_sample_p_my_trace
#define my_trace_fraction_control mcctest_sample_my_trace_fraction_control
#define k mcctest_sample_k
#define k_new mcctest_sample_k_new
#define k_old mcctest_sample_k_old
#define v_length mcctest_sample_v_length
#define my_sum mcctest_sample_my_sum
#define my_sum_plus_abs mcctest_sample_my_sum_plus_abs
#define culmative_probability mcctest_sample_culmative_probability
#define mc_prop mcctest_sample_mc_prop
#define time_to_scattering mcctest_sample_time_to_scattering
#define length_to_scattering mcctest_sample_length_to_scattering
#define length_to_boundery mcctest_sample_length_to_boundery
#define time_to_boundery mcctest_sample_time_to_boundery
#define selected_process mcctest_sample_selected_process
#define scattering_event mcctest_sample_scattering_event
#define time_propagated_without_scattering mcctest_sample_time_propagated_without_scattering
#define a_next_volume_found mcctest_sample_a_next_volume_found
#define next_volume mcctest_sample_next_volume
#define next_volume_priority mcctest_sample_next_volume_priority
#define done mcctest_sample_done
#define current_volume mcctest_sample_current_volume
#define number_of_solutions mcctest_sample_number_of_solutions
#define number_of_solutions_static mcctest_sample_number_of_solutions_static
#define check mcctest_sample_check
#define start mcctest_sample_start
#define intersection_with_children mcctest_sample_intersection_with_children
#define geometry_output mcctest_sample_geometry_output
#define tree_next_volume mcctest_sample_tree_next_volume
#define pre_allocated1 mcctest_sample_pre_allocated1
#define pre_allocated2 mcctest_sample_pre_allocated2
#define pre_allocated3 mcctest_sample_pre_allocated3
#define ray_position mcctest_sample_ray_position
#define ray_velocity mcctest_sample_ray_velocity
#define ray_velocity_final mcctest_sample_ray_velocity_final
#define volume_0_found mcctest_sample_volume_0_found
#define scattered_flag mcctest_sample_scattered_flag
#define scattered_flag_VP mcctest_sample_scattered_flag_VP
#define master_transposed_rotation_matrix mcctest_sample_master_transposed_rotation_matrix
#define temp_rotation_matrix mcctest_sample_temp_rotation_matrix
#define non_rotated_position mcctest_sample_non_rotated_position
#define rotated_position mcctest_sample_rotated_position
#define enable_tagging mcctest_sample_enable_tagging
#define stop_tagging_ray mcctest_sample_stop_tagging_ray
#define stop_creating_nodes mcctest_sample_stop_creating_nodes
#define enable_tagging_check mcctest_sample_enable_tagging_check
#define master_tagging_node_list mcctest_sample_master_tagging_node_list
#define current_tagging_node mcctest_sample_current_tagging_node
#define tagging_leaf_counter mcctest_sample_tagging_leaf_counter
#define number_of_scattering_events mcctest_sample_number_of_scattering_events
#define real_transmission_probability mcctest_sample_real_transmission_probability
#define mc_transmission_probability mcctest_sample_mc_transmission_probability
#define number_of_masks mcctest_sample_number_of_masks
#define number_of_masked_volumes mcctest_sample_number_of_masked_volumes
#define need_to_run_within_which_volume mcctest_sample_need_to_run_within_which_volume
#define mask_index_main mcctest_sample_mask_index_main
#define mask_iterate mcctest_sample_mask_iterate
#define mask_status_list mcctest_sample_mask_status_list
#define current_mask_intersect_list_status mcctest_sample_current_mask_intersect_list_status
#define mask_volume_index_list mcctest_sample_mask_volume_index_list
#define geometry_component_index_list mcctest_sample_geometry_component_index_list
#define Volume_copies_allocated mcctest_sample_Volume_copies_allocated
#define p_old mcctest_sample_p_old
#define this_logger mcctest_sample_this_logger
#define conditional_status mcctest_sample_conditional_status
#define tagging_conditional_list mcctest_sample_tagging_conditional_list
#define free_tagging_conditioanl_list mcctest_sample_free_tagging_conditioanl_list
#define logger_conditional_extend_array mcctest_sample_logger_conditional_extend_array
#define tagging_conditional_extend mcctest_sample_tagging_conditional_extend
#define max_conditional_extend_index mcctest_sample_max_conditional_extend_index
#define safty_distance mcctest_sample_safty_distance
#define safty_distance2 mcctest_sample_safty_distance2
#define number_of_processes_array mcctest_sample_number_of_processes_array
#define temporary_focus_data mcctest_sample_temporary_focus_data
#define focus_data_index mcctest_sample_focus_data_index
#define allow_inside_start mcctest_sample_allow_inside_start
#define history_limit mcctest_sample_history_limit
#define enable_conditionals mcctest_sample_enable_conditionals
#define inherit_number_of_scattering_events mcctest_sample_inherit_number_of_scattering_events
#line 208 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../contrib/union/Union_master.comp"
{
  // Use sanitation
  #ifndef ANY_GEOMETRY_DETECTOR_DECLARE
    printf("\nERROR: Need to define at least one Volume using Union_cylinder or Union_box before using the Union_master component. \n");
    exit(1);
  #endif
  #ifdef ANY_GEOMETRY_DETECTOR_DECLARE
    if (global_geometry_list.num_elements == 0) {
      printf("\nERROR: Need to define at least one Volume using Union_cylinder or Union_box before using the Union_master component. \n");
      printf("       Union_master component named \"%s\" is before any Volumes in the instrument file. At least one Volume need to be defined before\n",NAME_CURRENT_COMP);
    
      exit(1);
    }
  #endif
  
  // Parameters describing the safety distances close to surfaces, as scattering should not occur closer to a surface than the
  //  accuracy of the intersection calculation.
  safty_distance = 1E-11;
  safty_distance2 = safty_distance*2;
  
  // Write information to the global_master_list about the current Union_master
  sprintf(global_master_element.name,NAME_CURRENT_COMP);
  global_master_element.component_index = INDEX_CURRENT_COMP;
  add_element_to_master_list(&global_master_list,global_master_element);
  if (inherit_number_of_scattering_events == 1 && global_master_list.num_elements == 1) {
     printf("ERROR in Union_master with name %s. Inherit_number_of_scattering_events set to 1 for first Union_master component, but there is no preceeding Union_master component. Aborting.\n",NAME_CURRENT_COMP);
     exit(1);
  }
  this_global_master_index = global_master_list.num_elements - 1; // Save the index for this master in global master list
  
  // Set the component index of the previous Union_master component if one exists
  if (global_master_list.num_elements == 1) previous_master_index = 0; // no previous index
  else previous_master_index = global_master_list.elements[global_master_list.num_elements-2].component_index; // -2 because of zero indexing and needing the previous index.
  //printf("Assigned previous_master_index = %d \n",previous_master_index);
  
  // All volumes in the global_geometry_list is being check for activity using the number_of_activations input made for each geometry (default is 1)
  // In addition it is counted how many volumes, mask volumes and masked volumes are active in this Union_master.
  number_of_volumes = 1; // Starting with 1 as the surrounding vacuum is considered a volume
  number_of_masks = 0;   // Starting with 0 mask volumes
  number_of_masked_volumes = 0; // Starting with 0 masked volumes
  for (iterate=0;iterate<global_geometry_list.num_elements;iterate++) {
    if (global_geometry_list.elements[iterate].component_index < INDEX_CURRENT_COMP && global_geometry_list.elements[iterate].activation_counter > 0) {
        global_geometry_list.elements[iterate].active = 1;
        global_geometry_list.elements[iterate].activation_counter--;
        number_of_volumes++;
        if (global_geometry_list.elements[iterate].Volume->geometry.is_mask_volume == 1) number_of_masks++;
        if (global_geometry_list.elements[iterate].Volume->geometry.is_masked_volume == 1) number_of_masked_volumes++;
    } else global_geometry_list.elements[iterate].active = 0;
  }
  //printf("Found number of volumes to be %d \n",number_of_volumes);
  
  // Allocation of global lists
  geometry_component_index_list.num_elements = number_of_volumes;
  geometry_component_index_list.elements = malloc( geometry_component_index_list.num_elements * sizeof(int));
  mask_volume_index_list.num_elements = number_of_masks;
  if (number_of_masks >0) mask_volume_index_list.elements = malloc( number_of_masks * sizeof(int));
  mask_status_list.num_elements = number_of_masks;
  if (number_of_masks >0) mask_status_list.elements = malloc( number_of_masks * sizeof(int));
  current_mask_intersect_list_status.num_elements = number_of_masked_volumes;
  if (number_of_masked_volumes >0) current_mask_intersect_list_status.elements = malloc( number_of_masked_volumes * sizeof(int));
  
  // Make a list of component index from each volume index
  volume_index = 0;
  for (iterate=0;iterate<global_geometry_list.num_elements;iterate++) {
    if (global_geometry_list.elements[iterate].active == 1)
        geometry_component_index_list.elements[++volume_index] = global_geometry_list.elements[iterate].component_index;
      
  }
  geometry_component_index_list.elements[0] = 0; // Volume 0 is never set in the above code, but should never be used.
  
  // The input for this component is done through a series of input components
  // All information needed is stored in global lists, some of which is printed here for an overview to the user.
  MPI_MASTER( // MPI_MASTER ensures just one thread output this information to the user
  if (verbal == 1) {
      printf("---------------------------------------------------------------------\n");
      printf("global_process_list.num_elements: %d\n",global_process_list.num_elements);
      for (iterate=0;iterate<global_process_list.num_elements;iterate++) {
      printf("name of process [%d]: %s \n",iterate,global_process_list.elements[iterate].name);
      printf("component index [%d]: %d \n",iterate,global_process_list.elements[iterate].component_index);
      }
      
      printf("---------------------------------------------------------------------\n");
      printf("global_material_list.num_elements: %d\n",global_material_list.num_elements);
      for (iterate=0;iterate<global_material_list.num_elements;iterate++) {
      printf("name of material    [%d]: %s \n",iterate,global_material_list.elements[iterate].name);
      printf("component index     [%d]: %d \n",iterate,global_material_list.elements[iterate].component_index);
      printf("my_absoprtion       [%d]: %f \n",iterate,global_material_list.elements[iterate].physics->my_a);
      printf("number of processes [%d]: %d \n",iterate,global_material_list.elements[iterate].physics->number_of_processes);
      }

      printf("---------------------------------------------------------------------\n");
      printf("global_geometry_list.num_elements: %d\n",global_material_list.num_elements);
      for (iterate=0;iterate<global_geometry_list.num_elements;iterate++) {
        if (global_geometry_list.elements[iterate].active == 1) {
          printf("\n");
          printf("name of geometry    [%d]: %s \n",iterate,global_geometry_list.elements[iterate].name);
          printf("component index     [%d]: %d \n",iterate,global_geometry_list.elements[iterate].component_index);
          printf("Volume.name         [%d]: %s \n",iterate,global_geometry_list.elements[iterate].Volume->name);
          if (global_geometry_list.elements[iterate].Volume->geometry.is_mask_volume == 0) {
          printf("Volume.p_physics.is_vacuum           [%d]: %d \n",iterate,global_geometry_list.elements[iterate].Volume->p_physics->is_vacuum);
          printf("Volume.p_physics.my_absorption       [%d]: %f \n",iterate,global_geometry_list.elements[iterate].Volume->p_physics->my_a);
          printf("Volume.p_physics.number of processes [%d]: %d \n",iterate,global_geometry_list.elements[iterate].Volume->p_physics->number_of_processes);
          }
          printf("Volume.geometry.shape                [%d]: %s \n",iterate,global_geometry_list.elements[iterate].Volume->geometry.shape);
          printf("Volume.geometry.center.x             [%d]: %f \n",iterate,global_geometry_list.elements[iterate].Volume->geometry.center.x);
          printf("Volume.geometry.center.y             [%d]: %f \n",iterate,global_geometry_list.elements[iterate].Volume->geometry.center.y);
          printf("Volume.geometry.center.z             [%d]: %f \n",iterate,global_geometry_list.elements[iterate].Volume->geometry.center.z);
          printf("Volume.geometry.rotation_matrix[0]           [%d]: [%f %f %f] \n",iterate,global_geometry_list.elements[iterate].Volume->geometry.rotation_matrix[0][0],global_geometry_list.elements[iterate].Volume->geometry.rotation_matrix[0][1],global_geometry_list.elements[iterate].Volume->geometry.rotation_matrix[0][2]);
          printf("Volume.geometry.rotation_matrix[1]           [%d]: [%f %f %f] \n",iterate,global_geometry_list.elements[iterate].Volume->geometry.rotation_matrix[1][0],global_geometry_list.elements[iterate].Volume->geometry.rotation_matrix[1][1],global_geometry_list.elements[iterate].Volume->geometry.rotation_matrix[1][2]);
          printf("Volume.geometry.rotation_matrix[2]           [%d]: [%f %f %f] \n",iterate,global_geometry_list.elements[iterate].Volume->geometry.rotation_matrix[2][0],global_geometry_list.elements[iterate].Volume->geometry.rotation_matrix[2][1],global_geometry_list.elements[iterate].Volume->geometry.rotation_matrix[2][2]);
          if (strcmp(global_geometry_list.elements[iterate].Volume->geometry.shape,"cylinder") == 0) {
          printf("Volume.geometry.geometry_parameters.cyl_radius [%d]: %f \n",iterate,global_geometry_list.elements[iterate].Volume->geometry.geometry_parameters.p_cylinder_storage->cyl_radius);
          printf("Volume.geometry.geometry_parameters.height [%d]: %f \n",iterate,global_geometry_list.elements[iterate].Volume->geometry.geometry_parameters.p_cylinder_storage->height);
          }
          printf("Volume.geometry.focus_data_array.elements[0].Aim             [%d]: [%f %f %f] \n",iterate,global_geometry_list.elements[iterate].Volume->geometry.focus_data_array.elements[0].Aim.x,global_geometry_list.elements[iterate].Volume->geometry.focus_data_array.elements[0].Aim.y,global_geometry_list.elements[iterate].Volume->geometry.focus_data_array.elements[0].Aim.z);
        }
      }
      printf("---------------------------------------------------------------------\n");
      printf("number_of_volumes = %d\n",number_of_volumes);
      printf("number_of_masks = %d\n",number_of_masks);
      printf("number_of_masked_volumes = %d\n",number_of_masked_volumes);
  }
  ); // End MPI_MASTER
  
  
  // --- Initialization tasks independent of volume stucture -----------------------

  // Store a pointer to the conditional list and update the current index in that structure
  // If no tagging_conditionals were defined between this and the previous master, a dummy is allocated instead
  if (global_tagging_conditional_list.num_elements == global_tagging_conditional_list.current_index + 1) {
    tagging_conditional_list = &global_tagging_conditional_list.elements[global_tagging_conditional_list.current_index++].conditional_list;
    free_tagging_conditioanl_list = 0;
  } else {
    tagging_conditional_list = malloc(sizeof(struct conditional_list_struct));
    tagging_conditional_list->num_elements = 0;
    free_tagging_conditioanl_list = 1;
  }
  
  // Find the maximum logger extend index so that the correct memory allocation can be performed later
  // Here the loggers applied to all volumes are searched, later this result is compared to volume specific loggers and updated
  max_conditional_extend_index = -1;
  for (iterate=0;iterate<global_all_volume_logger_list.num_elements;iterate++) {
    if (global_all_volume_logger_list.elements[iterate].logger->logger_extend_index > max_conditional_extend_index) {
      max_conditional_extend_index = global_all_volume_logger_list.elements[iterate].logger->logger_extend_index;
    }
  }
  
  // The absolute rotation of this component is saved for use in initialization
  rot_transpose(ROT_A_CURRENT_COMP,master_transposed_rotation_matrix);
  
  // Preceeding componnets can add coordinates and rotations to global_positions_to_transform and global_rotations_to_transform
  //  in order to have these transformed into the coordinate system of the next master compoent in the instrument file.
  // Here these transformations are performed, and the lists are cleared so no transformed information is further altered by
  //  next master components.
  
  // Position transformation
  for (iterate=0;iterate<global_positions_to_transform_list.num_elements;iterate++) {
    //print_position(*global_positions_to_transform_list.positions[iterate],"Position to transform before");
    non_rotated_position = coords_sub(*(global_positions_to_transform_list.positions[iterate]),POS_A_CURRENT_COMP);
    *(global_positions_to_transform_list.positions[iterate]) = rot_apply(ROT_A_CURRENT_COMP,non_rotated_position);
    //print_position(*global_positions_to_transform_list.positions[iterate],"Position to transform after");
  }
  if (global_positions_to_transform_list.num_elements > 0) {
    global_positions_to_transform_list.num_elements = 0;
    free(global_positions_to_transform_list.positions);
  }
  // Rotation transformation
  for (iterate=0;iterate<global_rotations_to_transform_list.num_elements;iterate++) {
    //print_rotation(*(global_rotations_to_transform_list.rotations[iterate]),"rotation matrix to be updated");
    rot_mul(master_transposed_rotation_matrix,*(global_rotations_to_transform_list.rotations[iterate]),temp_rotation_matrix);
    rot_copy(*(global_rotations_to_transform_list.rotations[iterate]),temp_rotation_matrix);
  }
  if (global_rotations_to_transform_list.num_elements > 0) {
    global_rotations_to_transform_list.num_elements = 0;
    free(global_rotations_to_transform_list.rotations);
  }

  
  // --- Definition of volumes and loading of appropriate data  -----------------------
  
  // The information stored in global lists is to be stored in one array of structures that is allocated here
  Volumes = malloc(number_of_volumes * sizeof(struct Volume_struct*));
  scattered_flag = malloc(number_of_volumes*sizeof(int));
  scattered_flag_VP = (int**) malloc(number_of_volumes * sizeof(int*));
  number_of_processes_array = malloc(number_of_volumes*sizeof(int));
  
  // The mcdisplay functions need access to the other geomtries, but can not use the Volumes struct because of order of definition.
  // A separate list of pointers to the geometry structures is thus allocated
  Geometries = malloc(number_of_volumes * sizeof(struct geometry_struct *));
  
  // When activation counter is used to have several copies of one volume, it can become necessary to have soft copies of volumes
  // Not all of these will necessarily be allocated or used.
  Volume_copies = malloc(number_of_volumes * sizeof(struct Volume_struct *));
  Volume_copies_allocated.num_elements = 0;
  
  // The central structure is called a "Volume", it describes a region in space with certain scattering processes and absorption cross section

  // ---  Volume 0 ------------------------------------------------------------------------------------------------
  // Volume 0 is the vacuum surrounding the experiment (infinite, everywhere) and its properties are hardcoded here
  Volumes[0] = malloc(sizeof(struct Volume_struct));
  strcpy(Volumes[0]->name,"Surrounding vacuum");
  // Assign geometry
  
  // This information is meaningless for volume 0, and is never be acsessed in the logic.
  Volumes[0]->geometry.priority_value = 0.0;
  Volumes[0]->geometry.center.x = 0;
  Volumes[0]->geometry.center.y = 0;
  Volumes[0]->geometry.center.z = 0;
  strcpy(Volumes[0]->geometry.shape,"vacuum");
  Volumes[0]->geometry.within_function = &r_within_surroundings; // Always returns 1
  // No physics struct allocated
  Volumes[0]->p_physics = NULL;
  number_of_processes_array[volume_index] = 0;
  
  // These are never used for volume 0, but by setting the length to 0 it is automatically skipped in many forloops without the need for an if statement
  Volumes[0]->geometry.children.num_elements=0;
  Volumes[0]->geometry.direct_children.num_elements=0;
  Volumes[0]->geometry.destinations_list.num_elements=0;
  Volumes[0]->geometry.reduced_destinations_list.num_elements=0;
  
  Volumes[0]->geometry.masked_by_list.num_elements = 0;
  Volumes[0]->geometry.mask_list.num_elements = 0;
  Volumes[0]->geometry.masked_by_mask_index_list.num_elements = 0;
  Volumes[0]->geometry.mask_mode=0;
  Volumes[0]->geometry.is_mask_volume=0;
  Volumes[0]->geometry.is_masked_volume=0;
  
  // A pointer to the geometry structure
  Geometries[0] = &Volumes[0]->geometry;
  
  // Logging initialization
  Volumes[0]->loggers.num_elements = 0;
  
  
  // --- Loop over user defined volumes ------------------------------------------------------------------------
  // Here the user defined volumes are loaded into the volume structure that is used in the ray-tracing
  //  algorithm. Not all user defined volumes are used, some could be used by a previous master, some
  //  could be used by the previous master, this one, and perhaps more. This is controlled by the
  //  activation counter input for geometries, and is here condensed to the active variable.
  // Volumes that were used before
  
  
  max_number_of_processes = 0; // The maximum number of processes in a volume is assumed 0 and updated during the following loop
  
  volume_index = 0;
  mask_index_main = 0;
  for (geometry_list_index=0;geometry_list_index<global_geometry_list.num_elements;geometry_list_index++) {
    if (global_geometry_list.elements[geometry_list_index].active == 1) { // Only include the volume if it is active
      volume_index++;
      // Connect a volume for each of the geometry.comp instances in the McStas instrument files
      if (global_geometry_list.elements[geometry_list_index].activation_counter == 0) {
        // This is the last time this volume is used, use the hard copy from the geometry component
        Volumes[volume_index] = global_geometry_list.elements[geometry_list_index].Volume;
        //printf("used hard copy of volume %d \n",volume_index);
      } else {
        // Since this volume is still needed more than this once, we need to make a shallow copy and use instead
        
        //printf("\n");
        //printf("making local copy of volume %d \n",volume_index);
        Volume_copies[volume_index] = malloc(sizeof(struct Volume_struct));
        *(Volume_copies[volume_index]) = *global_geometry_list.elements[geometry_list_index].Volume; // Makes shallow copy
        Volumes[volume_index] = Volume_copies[volume_index];
        add_element_to_int_list(&Volume_copies_allocated,volume_index); // Keep track of dynamically allocated volumes in order to free them in FINALLY.
        
        // The geometry storage needs a shallow copy as well (hard copy not necessary for any current geometries), may need changes in future
        // A simple copy_geometry_parameters function is added to the geometry in each geometry component
        Volumes[volume_index]->geometry.geometry_parameters = Volumes[volume_index]->geometry.copy_geometry_parameters(&global_geometry_list.elements[geometry_list_index].Volume->geometry.geometry_parameters);
          
      }
      
      // This section identifies the different non isotropic processes in the current volume and give them appropriate transformation matrices
      // Identify the number of non isotropic processes in a material (this code can be safely executed for the same material many times)
      // A setting of -1 means no transformation necessary, other settings are assigned a unique identifier instead
      non_isotropic_found = 0;
      for (iterate=0;iterate<Volumes[volume_index]->p_physics->number_of_processes;iterate++) {
        if (Volumes[volume_index]->p_physics->p_scattering_array[iterate].non_isotropic_rot_index != -1) {
            Volumes[volume_index]->p_physics->p_scattering_array[iterate].non_isotropic_rot_index = non_isotropic_found;
            non_isotropic_found++;
        }
      }
      
      Volumes[volume_index]->geometry.focus_array_indices.num_elements=0;
      // For the non_isotropic volumes found, rotation matrices need to be allocated and calculated
      if (non_isotropic_found > 0) {
        // Allocation of rotation and transpose rotation matrices
        if (Volumes[volume_index]->geometry.process_rot_allocated == 0) {
          Volumes[volume_index]->geometry.process_rot_matrix_array = malloc(non_isotropic_found * sizeof(Rotation));
          Volumes[volume_index]->geometry.transpose_process_rot_matrix_array = malloc(non_isotropic_found * sizeof(Rotation));
          Volumes[volume_index]->geometry.process_rot_allocated = 1;
        }
      
        // Calculation of the appropriate rotation matrices for transformation between Union_master and the process in a given volume.
        non_isotropic_found = 0;
        for (iterate=0;iterate<Volumes[volume_index]->p_physics->number_of_processes;iterate++) {
          if (Volumes[volume_index]->p_physics->p_scattering_array[iterate].non_isotropic_rot_index != -1) {
            // Transformation for each process / geometry combination
            
            // The focus vector is given in relation to the geometry and needs to be transformed to the process
            // Work on temporary_focus_data_element which is added to the focus_data_array_at the end
            temporary_focus_data = Volumes[volume_index]->geometry.focus_data_array.elements[0];
            
            // Correct for process rotation
            temporary_focus_data.Aim = rot_apply(Volumes[volume_index]->p_physics->p_scattering_array[iterate].rotation_matrix,temporary_focus_data.Aim);
            
            // Add element to focus_array_indices
            // focus_array_indices refers to the correct element in focus_data_array for this volume/process combination
            // focus_data_array[0] is the isotropic version in all cases, so the first non_isotropic goes to focus_data_array[1]
            //  and so forth. When a process is isotropic, this array is appended with a zero.
            // The focus_array_indices maps process numbers to the correct focus_data_array index.
            add_element_to_int_list(&Volumes[volume_index]->geometry.focus_array_indices,non_isotropic_found+1);
            
            // Add the new focus_data element to this volumes focus_data_array.
            add_element_to_focus_data_array(&Volumes[volume_index]->geometry.focus_data_array,temporary_focus_data);
            
            // Quick error check to see the length is correct which indirectly confirms the indices are correct
            if (Volumes[volume_index]->geometry.focus_data_array.num_elements != non_isotropic_found + 2) {
              printf("ERROR, focus_data_array length for volume %s inconsistent with number of non isotropic processes found!\n",Volumes[volume_index]->name);
              exit(1);
            }
            
            // Create rotation matrix for this specific volume / process combination to transform from master coordinate system to the non-isotropics process coordinate system
            // This is done by multipling the transpose master component roration matrix, the volume rotation, and then the process rotation matrix onto the velocity / wavevector
            rot_mul(Volumes[volume_index]->geometry.rotation_matrix,master_transposed_rotation_matrix,temp_rotation_matrix);
            rot_mul(Volumes[volume_index]->p_physics->p_scattering_array[iterate].rotation_matrix,temp_rotation_matrix,Volumes[volume_index]->geometry.process_rot_matrix_array[non_isotropic_found]);
            
            // Need to transpose as well to transform back to the master coordinate system
            rot_transpose(Volumes[volume_index]->geometry.process_rot_matrix_array[non_isotropic_found],Volumes[volume_index]->geometry.transpose_process_rot_matrix_array[non_isotropic_found]);

            // Debug print
            //print_rotation(Volumes[volume_index]->geometry.process_rot_matrix_array[non_isotropic_found],"Process rotation matrix");
            //print_rotation(Volumes[volume_index]->geometry.transpose_process_rot_matrix_array[non_isotropic_found],"Transpose process rotation matrix");
            
            non_isotropic_found++;
          } else {
            // This process can use the standard isotropic focus_data_array which is indexed zero.
            add_element_to_int_list(&Volumes[volume_index]->geometry.focus_array_indices,0);
          }
        }
      } else {
      // No non isotropic volumes found, focus_array_indices should just be a list of 0's of same length as the number of processes.
      // In this way all processes use the isotropic focus_data structure
        Volumes[volume_index]->geometry.focus_array_indices.elements = malloc(Volumes[volume_index]->p_physics->number_of_processes * sizeof(int));
        for (iterate=0;iterate<Volumes[volume_index]->p_physics->number_of_processes;iterate++)
          Volumes[volume_index]->geometry.focus_array_indices.elements[iterate] = 0;
          
      }
      //print_1d_int_list(Volumes[volume_index]->geometry.focus_array_indices,"focus_array_indices");
      
      
      // This component works in its local coordinate system, and thus all information from the input components should be transformed to its coordinate system.
      // All the input components saved their absolute rotation/position into their Volume structure, and the absolute rotation of the current component is known.
      // The next section finds the relative rotation and translation of all the volumes and the master component.
      
      // Transform the rotation matrices for each volume
      rot_mul(ROT_A_CURRENT_COMP,Volumes[volume_index]->geometry.transpose_rotation_matrix,temp_rotation_matrix);
      // Copy the result back to the volumes structure
      rot_copy(Volumes[volume_index]->geometry.rotation_matrix,temp_rotation_matrix);
      // Now update the transpose as well
      rot_transpose(Volumes[volume_index]->geometry.rotation_matrix,temp_rotation_matrix);
      rot_copy(Volumes[volume_index]->geometry.transpose_rotation_matrix,temp_rotation_matrix);
        
      // Transform the position for each volume
      non_rotated_position.x = Volumes[volume_index]->geometry.center.x - POS_A_CURRENT_COMP.x;
      non_rotated_position.y = Volumes[volume_index]->geometry.center.y - POS_A_CURRENT_COMP.y;
      non_rotated_position.z = Volumes[volume_index]->geometry.center.z - POS_A_CURRENT_COMP.z;

      rot_transpose(ROT_A_CURRENT_COMP,temp_rotation_matrix); // REVIEW LINE
      rotated_position = rot_apply(ROT_A_CURRENT_COMP,non_rotated_position);

      Volumes[volume_index]->geometry.center.x = rotated_position.x;
      Volumes[volume_index]->geometry.center.y = rotated_position.y;
      Volumes[volume_index]->geometry.center.z = rotated_position.z;
      
      // The focus_data information need to be updated as well
      rot_mul(ROT_A_CURRENT_COMP,Volumes[volume_index]->geometry.focus_data_array.elements[0].absolute_rotation,temp_rotation_matrix);
      // Copy the result back to the volumes structure
      rot_copy(Volumes[volume_index]->geometry.focus_data_array.elements[0].absolute_rotation,temp_rotation_matrix);
      
      // Use same rotation on the aim vector of the isotropic focus_data element
      Volumes[volume_index]->geometry.focus_data_array.elements[0].Aim = rot_apply(Volumes[volume_index]->geometry.rotation_matrix,Volumes[volume_index]->geometry.focus_data_array.elements[0].Aim);
      
      // To allocate enough memory to hold information on all processes, the maximum of these is updated if this volume has more
      if (Volumes[volume_index]->p_physics->number_of_processes > max_number_of_processes)
          max_number_of_processes = Volumes[volume_index]->p_physics->number_of_processes;
        
      // Allocate memory to scattered_flag_VP (holds statistics for scatterings in each process of the volume)
      scattered_flag_VP[volume_index] = malloc(Volumes[volume_index]->p_physics->number_of_processes * sizeof(int));
      number_of_processes_array[volume_index] = Volumes[volume_index]->p_physics->number_of_processes;
      
      // Normalizing and error checking process interact fraction
      number_of_process_interacts_set = 0; total_process_interact=0;
      for (process_index=0;process_index<Volumes[volume_index]->p_physics->number_of_processes;process_index++) {
        if (Volumes[volume_index]->p_physics->p_scattering_array[process_index].process_p_interact != -1) {
          number_of_process_interacts_set++;
          total_process_interact += Volumes[volume_index]->p_physics->p_scattering_array[process_index].process_p_interact;
        } else {
          index_of_lacking_process = process_index;
        }
      }
      
      if (number_of_process_interacts_set == 0) Volumes[volume_index]->p_physics->interact_control = 0;
      else Volumes[volume_index]->p_physics->interact_control = 1;
      
      // If all are set, check if they need renormalization so that the sum is one.
      if (number_of_process_interacts_set == Volumes[volume_index]->p_physics->number_of_processes) {
        if (total_process_interact > 1.001 || total_process_interact < 0.999) {
          for (process_index=0;process_index<Volumes[volume_index]->p_physics->number_of_processes;process_index++) {
            Volumes[volume_index]->p_physics->p_scattering_array[process_index].process_p_interact = Volumes[volume_index]->p_physics->p_scattering_array[process_index].process_p_interact/total_process_interact;
          }
        }
      } else if ( number_of_process_interacts_set != 0) {
        if (number_of_process_interacts_set == Volumes[volume_index]->p_physics->number_of_processes - 1) {// If all but one is set, it is an easy fix
          Volumes[volume_index]->p_physics->p_scattering_array[index_of_lacking_process].process_p_interact = 1 - total_process_interact;
          if (total_process_interact >= 1) {
            printf("ERROR, material %s has a total interact_fraction above 1 and a process without an interact_fraction. Either set all so they can be renormalized, or have a sum below 1, so that the last can have 1 - sum.\n",Volumes[volume_index]->p_physics->name);
            exit(1);
          }
        } else {
            printf("ERROR, material %s needs to have all, all minus one or none of its processes with an interact_fraction \n",Volumes[volume_index]->p_physics->name);
            exit(1);
        }
      }
        
      // Some initialization can only happen after the rotation matrix relative to the master is known
      // Such initialization is placed in the geometry component, and executed here through a function pointer
      Volumes[volume_index]->geometry.initialize_from_main_function(&Volumes[volume_index]->geometry);
      
      // Add pointer to geometry to Geometries
      Geometries[volume_index] = &Volumes[volume_index]->geometry;
      
      // Initialize mask intersect list
      Volumes[volume_index]->geometry.mask_intersect_list.num_elements = 0;
      
      // Here the mask_list and masked_by_list for the volume is updated from component index values to volume indexes
      for (iterate=0;iterate<Volumes[volume_index]->geometry.mask_list.num_elements;iterate++)
        Volumes[volume_index]->geometry.mask_list.elements[iterate] = find_on_int_list(geometry_component_index_list,Volumes[volume_index]->geometry.mask_list.elements[iterate]);
      
      for (iterate=0;iterate<Volumes[volume_index]->geometry.masked_by_list.num_elements;iterate++)
        Volumes[volume_index]->geometry.masked_by_list.elements[iterate] = find_on_int_list(geometry_component_index_list,Volumes[volume_index]->geometry.masked_by_list.elements[iterate]);
      
      // If the volume is a mask, its volume number is added to the mask_volume_index list so volume index can be converted to mask_index.
      if (Volumes[volume_index]->geometry.is_mask_volume == 1) Volumes[volume_index]->geometry.mask_index = mask_index_main;
      if (Volumes[volume_index]->geometry.is_mask_volume == 1) mask_volume_index_list.elements[mask_index_main++] = volume_index;
     
      // Check all loggers assosiated with this volume and update the max_conditional_extend_index if necessary
      //printf("reached max_test for volume %d \n",volume_index);
      for (iterate=0;iterate<Volumes[volume_index]->loggers.num_elements;iterate++) {
        //printf("iterate = %d \n",iterate);
        for (process_index=0;process_index<Volumes[volume_index]->loggers.p_logger_volume[iterate].num_elements;process_index++) {
        //printf("process_index = %d \n",process_index);
          if (Volumes[volume_index]->loggers.p_logger_volume[iterate].p_logger_process[process_index] != NULL) {
            if (Volumes[volume_index]->loggers.p_logger_volume[iterate].p_logger_process[process_index]->logger_extend_index > max_conditional_extend_index)
              max_conditional_extend_index = Volumes[volume_index]->loggers.p_logger_volume[iterate].p_logger_process[process_index]->logger_extend_index;
          }
        }
      }
      //printf("did max_test for volume %d\n",volume_index);
      
    
    }
  } // Initialization for each volume done
  
  // ------- Initialization of ray-tracing algorithm ------------------------------------
   
  my_trace = malloc(max_number_of_processes*sizeof(double));
  my_trace_fraction_control = malloc(max_number_of_processes*sizeof(double));
  
  // All geometries can have 2 intersections currently, when this changes the maximum number of solutions need to be reported to the Union_master.
  number_of_solutions = &number_of_solutions_static;
  component_error_msg = 0;
  
  // Pre allocated memory for destination list search
  pre_allocated1 = malloc(number_of_volumes * sizeof(int));
  pre_allocated2 = malloc(number_of_volumes * sizeof(int));
  pre_allocated3 = malloc(number_of_volumes * sizeof(int));
  
  // Allocate memory for logger_conditional_extend_array used in the extend section of the master component, if it is needed.
  if (max_conditional_extend_index > -1) {
    logger_conditional_extend_array = malloc((max_conditional_extend_index + 1)*sizeof(int));
  }
  
  // In this function different lists of volume indecies are generated. They are the key to the speed of the component and central for the logic.
  // They use simple set algebra to generate these lists for each volume:
  // Children list for volume n: Indicies of volumes that are entirely within the set described by volume n
  // Overlap list for volume n: Indicies of volume that contains some of the set described by volume n (excluding volume n)
  // Intersect check list for volume n: Indicies of volumes to check for intersection if a ray originates from volume n (is generated from the children and overlap lists)
  // Parents list for volume n: Indicies of volumes that contain the entire set of volume n
  // Grandparents lists for volume n: Indicies of volumes that contain the entire set of at least one parent of volume n
  // Destination list for volume n: Indicies of volumes that could be the destination volume when a ray leaves volume n
  // The overlap, parents and grandparents lists are local variables in the function, and not in the main scope.
  
  generate_lists(Volumes, &starting_lists, number_of_volumes, list_verbal);
  
  // Generate "safe starting list", which contains all volumes that the ray may enter from other components
  // These are all volumes without scattering or absorption
  
  // Updating mask lists from volume index to global_mask_indices
  // Filling out the masked_by list that uses mask indices
  for (volume_index=0;volume_index<number_of_volumes;volume_index++) {
    Volumes[volume_index]->geometry.masked_by_mask_index_list.num_elements = Volumes[volume_index]->geometry.masked_by_list.num_elements;
    Volumes[volume_index]->geometry.masked_by_mask_index_list.elements = malloc(Volumes[volume_index]->geometry.masked_by_mask_index_list.num_elements * sizeof(int));
    for (iterate=0;iterate<Volumes[volume_index]->geometry.masked_by_list.num_elements;iterate++)
        Volumes[volume_index]->geometry.masked_by_mask_index_list.elements[iterate] = find_on_int_list(mask_volume_index_list,Volumes[volume_index]->geometry.masked_by_list.elements[iterate]);
  }
  
  // Optimizing speed of the within_which_volume search algorithm
  int volume_index_main; // REVIEW_LINE
  /*
  for (volume_index_main=0;volume_index_main<number_of_volumes;volume_index_main++) {
    Volumes[volume_index_main]->geometry.destinations_logic_list.elements[0] = 0;
  }
  */
  
  // Checking for equal priorities in order to alert the user to a potential input error
  for (volume_index_main=0;volume_index_main<number_of_volumes;volume_index_main++) {
    for (volume_index=0;volume_index<number_of_volumes;volume_index++)
        if (Volumes[volume_index_main]->geometry.priority_value == Volumes[volume_index]->geometry.priority_value && volume_index_main != volume_index) {
            if (Volumes[volume_index_main]->geometry.is_mask_volume == 0 && Volumes[volume_index]->geometry.is_mask_volume == 0) {
              // Priority of masks do not matter
              printf("ERROR in Union_master with name %s. The volumes named %s and %s have the same priority. Change the priorities so the one present in case of overlap has highest priority.\n",NAME_CURRENT_COMP,Volumes[volume_index_main]->name,Volumes[volume_index]->name);
              exit(1);
            }
        }
  }
  
  
  // Printing the generated lists for all volumes.
  MPI_MASTER(
  if (verbal) printf("\n ---- Overview of the lists generated for each volume ---- \n");
  
  
  printf("List overview for surrounding vacuum\n");
  for (volume_index_main=0;volume_index_main<number_of_volumes;volume_index_main++) {
      
      if (verbal) {
        if (volume_index_main != 0) {
          if (volume_index_main,Volumes[volume_index_main]->geometry.is_mask_volume == 0 ||
              volume_index_main,Volumes[volume_index_main]->geometry.is_masked_volume == 0 ||
              volume_index_main,Volumes[volume_index_main]->geometry.is_exit_volume == 0) {
            printf("List overview for %s with %s shape made of %s\n",
                   Volumes[volume_index_main]->name,
                   Volumes[volume_index_main]->geometry.shape,
                   Volumes[volume_index_main]->p_physics->name);
          } else {
            printf("List overview for %s with shape %s\n",
                   Volumes[volume_index_main]->name,
                   Volumes[volume_index_main]->geometry.shape);
          }
        }
      }
  
      if (verbal) sprintf(string_output,"Children for Volume                  %d",volume_index_main);
      if (verbal) print_1d_int_list(Volumes[volume_index_main]->geometry.children,string_output);
      
      if (verbal) sprintf(string_output,"Direct_children for Volume           %d",volume_index_main);
      if (verbal) print_1d_int_list(Volumes[volume_index_main]->geometry.direct_children,string_output);
      
      if (verbal) sprintf(string_output,"Intersect_check_list for Volume      %d",volume_index_main);
      if (verbal) print_1d_int_list(Volumes[volume_index_main]->geometry.intersect_check_list,string_output);
      
      if (verbal) sprintf(string_output,"Mask_intersect_list for Volume       %d",volume_index_main);
      if (verbal) print_1d_int_list(Volumes[volume_index_main]->geometry.mask_intersect_list,string_output);
      
      if (verbal) sprintf(string_output,"Destinations_list for Volume         %d",volume_index_main);
      if (verbal) print_1d_int_list(Volumes[volume_index_main]->geometry.destinations_list,string_output);
      
      //if (verbal) sprintf(string_output,"Destinations_logic_list for Volume   %d",volume_index_main);
      //if (verbal) print_1d_int_list(Volumes[volume_index_main]->geometry.destinations_logic_list,string_output);
      
      if (verbal) sprintf(string_output,"Reduced_destinations_list for Volume %d",volume_index_main);
      if (verbal) print_1d_int_list(Volumes[volume_index_main]->geometry.reduced_destinations_list,string_output);
      
      if (verbal) sprintf(string_output,"Next_volume_list for Volume          %d",volume_index_main);
      if (verbal) print_1d_int_list(Volumes[volume_index_main]->geometry.next_volume_list,string_output);
      
      if (verbal) {
        if (volume_index_main != 0)
                           printf("      Is_vacuum for Volume                 %d = %d\n",volume_index_main,Volumes[volume_index_main]->p_physics->is_vacuum);
      }
      if (verbal) {
        if (volume_index_main != 0)
                           printf("      is_mask_volume for Volume            %d = %d\n",volume_index_main,Volumes[volume_index_main]->geometry.is_mask_volume);
      }
      if (verbal) {
        if (volume_index_main != 0)
                           printf("      is_masked_volume for Volume          %d = %d\n",volume_index_main,Volumes[volume_index_main]->geometry.is_masked_volume);
      }
      if (verbal) {
        if (volume_index_main != 0)
                           printf("      is_exit_volume for Volume            %d = %d\n",volume_index_main,Volumes[volume_index_main]->geometry.is_exit_volume);
      }
      
      if (verbal) sprintf(string_output,"mask_list for Volume                 %d",volume_index_main);
      if (verbal) print_1d_int_list(Volumes[volume_index_main]->geometry.mask_list,string_output);
      
      if (verbal) sprintf(string_output,"masked_by_list for Volume            %d",volume_index_main);
      if (verbal) print_1d_int_list(Volumes[volume_index_main]->geometry.masked_by_list,string_output);
      
      if (verbal) sprintf(string_output,"masked_by_mask_index_list for Volume %d",volume_index_main);
      if (verbal) print_1d_int_list(Volumes[volume_index_main]->geometry.masked_by_mask_index_list,string_output);
      
      if (verbal)          printf("      mask_mode for Volume                 %d = %d\n",volume_index_main,Volumes[volume_index_main]->geometry.mask_mode);
      if (verbal) printf("\n");
  }
  ) // End of MPI_MASTER
  

  // Initializing intersection_time_table
  // The intersection time table contains all information on intersection times for the current position/direction, and is cleared everytime a ray changes direction.
  // Not all entries needs to be calculated, so there is a variable that keeps track of which intersection times have been calculated in order to avoid redoing that.
  // When the intersections times are calculated for a volume, all future intersections are kept in the time table.
  // Thus the memory allocation have to take into account how many intersections there can be with each volume, but it is currently set to 2, but can easily be changed. This may need to be reported by individual geometry components in the future.
  
  intersection_time_table.num_volumes = number_of_volumes;
  
  intersection_time_table.n_elements = (int*) malloc(intersection_time_table.num_volumes * sizeof(int));
  intersection_time_table.calculated = (int*) malloc(intersection_time_table.num_volumes * sizeof(int));
  intersection_time_table.intersection_times = (double**) malloc(intersection_time_table.num_volumes * sizeof(double));
  for (iterate = 0;iterate < intersection_time_table.num_volumes;iterate++){
      if (strcmp(Volumes[iterate]->geometry.shape, "mesh") == 0) {
        intersection_time_table.n_elements[iterate] = (int) 100; // Meshes can have any number of intersections, here we allocate room for 100
      } else {
        intersection_time_table.n_elements[iterate] = (int) 2; // number of intersection for all other geometries
      }
      if (iterate == 0) intersection_time_table.n_elements[iterate] = (int) 0; // number of intersection solutions
      intersection_time_table.calculated[iterate] = (int) 0; // Initializing calculated logic
      
      if (iterate == 0) {
           intersection_time_table.intersection_times[0] = NULL;
      }
      else {
        intersection_time_table.intersection_times[iterate] = (double*) malloc(intersection_time_table.n_elements[iterate]*sizeof(double));
        for (solutions = 0;solutions < intersection_time_table.n_elements[iterate];solutions++) {
              intersection_time_table.intersection_times[iterate][solutions] = -1.0;
        }
      }
  }
  
  // If enabled, the tagging system tracks all different histories sampled by the program.

  // Initialize the tagging tree
  // Allocate a list of host nodes with the same length as the number of volumes
  
  stop_creating_nodes = 0; stop_tagging_ray = 0; tagging_leaf_counter = 0;
  if (enable_tagging) {
    master_tagging_node_list.num_elements = number_of_volumes;
    master_tagging_node_list.elements = malloc(master_tagging_node_list.num_elements * sizeof(struct tagging_tree_node_struct*));
  
    // Initialize
    for (volume_index=0;volume_index<number_of_volumes;volume_index++) {
      //if (verbal) printf("Allocating master tagging node for volume number %d \n",volume_index);
      master_tagging_node_list.elements[volume_index] = initialize_tagging_tree_node(master_tagging_node_list.elements[volume_index], NULL, Volumes[volume_index]);
      //if (verbal) printf("Allocated master tagging node for volume number %d \n",volume_index);
    }
  }
  
  // Initialize loggers
  loggers_with_data_array.allocated_elements = 0;
  loggers_with_data_array.used_elements = 0;
  
  
  // Signal initialization complete
  MPI_MASTER(
    printf("Union_master component %s initialized sucessfully\n",NAME_CURRENT_COMP);
  )

}
#line 10707 "./Union_test_mask.c"
#undef inherit_number_of_scattering_events
#undef enable_conditionals
#undef history_limit
#undef allow_inside_start
#undef focus_data_index
#undef temporary_focus_data
#undef number_of_processes_array
#undef safty_distance2
#undef safty_distance
#undef max_conditional_extend_index
#undef tagging_conditional_extend
#undef logger_conditional_extend_array
#undef free_tagging_conditioanl_list
#undef tagging_conditional_list
#undef conditional_status
#undef this_logger
#undef p_old
#undef Volume_copies_allocated
#undef geometry_component_index_list
#undef mask_volume_index_list
#undef current_mask_intersect_list_status
#undef mask_status_list
#undef mask_iterate
#undef mask_index_main
#undef need_to_run_within_which_volume
#undef number_of_masked_volumes
#undef number_of_masks
#undef mc_transmission_probability
#undef real_transmission_probability
#undef number_of_scattering_events
#undef tagging_leaf_counter
#undef current_tagging_node
#undef master_tagging_node_list
#undef enable_tagging_check
#undef stop_creating_nodes
#undef stop_tagging_ray
#undef enable_tagging
#undef rotated_position
#undef non_rotated_position
#undef temp_rotation_matrix
#undef master_transposed_rotation_matrix
#undef scattered_flag_VP
#undef scattered_flag
#undef volume_0_found
#undef ray_velocity_final
#undef ray_velocity
#undef ray_position
#undef pre_allocated3
#undef pre_allocated2
#undef pre_allocated1
#undef tree_next_volume
#undef geometry_output
#undef intersection_with_children
#undef start
#undef check
#undef number_of_solutions_static
#undef number_of_solutions
#undef current_volume
#undef done
#undef next_volume_priority
#undef next_volume
#undef a_next_volume_found
#undef time_propagated_without_scattering
#undef scattering_event
#undef selected_process
#undef time_to_boundery
#undef length_to_boundery
#undef length_to_scattering
#undef time_to_scattering
#undef mc_prop
#undef culmative_probability
#undef my_sum_plus_abs
#undef my_sum
#undef v_length
#undef k_old
#undef k_new
#undef k
#undef my_trace_fraction_control
#undef p_my_trace
#undef my_trace
#undef process_start
#undef process
#undef min_intersection_time
#undef intersection_time
#undef time_found
#undef min_volume
#undef min_solution
#undef solution
#undef limit
#undef max_number_of_processes
#undef solutions
#undef process_index
#undef volume_index
#undef number_of_volumes
#undef string_output
#undef component_error_msg
#undef error_msg
#undef v
#undef r_start
#undef r
#undef starting_lists
#undef Geometries
#undef Volumes
#undef intersection_time_table
#undef geometry_list_index
#undef previous_master_index
#undef this_global_master_index
#undef global_master_element
#undef starting_volume_warning
#undef finally_verbal
#undef trace_verbal
#undef list_verbal
#undef verbal
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component m4pi. */
  SIG_MESSAGE("m4pi (Init)");
#define mccompcurname  m4pi
#define mccompcurtype  PSD_monitor_4PI
#define mccompcurindex 13
#define nx mccm4pi_nx
#define ny mccm4pi_ny
#define PSD_N mccm4pi_PSD_N
#define PSD_p mccm4pi_PSD_p
#define PSD_p2 mccm4pi_PSD_p2
#define filename mccm4pi_filename
#define radius mccm4pi_radius
#define restore_neutron mccm4pi_restore_neutron
#define nowritefile mccm4pi_nowritefile
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
#line 10851 "./Union_test_mask.c"
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

  /* Initializations for component m4pi_1. */
  SIG_MESSAGE("m4pi_1 (Init)");
#define mccompcurname  m4pi_1
#define mccompcurtype  PSD_monitor_4PI
#define mccompcurindex 14
#define nx mccm4pi_1_nx
#define ny mccm4pi_1_ny
#define PSD_N mccm4pi_1_PSD_N
#define PSD_p mccm4pi_1_PSD_p
#define PSD_p2 mccm4pi_1_PSD_p2
#define filename mccm4pi_1_filename
#define radius mccm4pi_1_radius
#define restore_neutron mccm4pi_1_restore_neutron
#define nowritefile mccm4pi_1_nowritefile
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
#line 10891 "./Union_test_mask.c"
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

  /* Initializations for component m4pi_2. */
  SIG_MESSAGE("m4pi_2 (Init)");
#define mccompcurname  m4pi_2
#define mccompcurtype  PSD_monitor_4PI
#define mccompcurindex 15
#define nx mccm4pi_2_nx
#define ny mccm4pi_2_ny
#define PSD_N mccm4pi_2_PSD_N
#define PSD_p mccm4pi_2_PSD_p
#define PSD_p2 mccm4pi_2_PSD_p2
#define filename mccm4pi_2_filename
#define radius mccm4pi_2_radius
#define restore_neutron mccm4pi_2_restore_neutron
#define nowritefile mccm4pi_2_nowritefile
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
#line 10931 "./Union_test_mask.c"
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

  /* Initializations for component arm_1. */
  SIG_MESSAGE("arm_1 (Init)");

  /* Initializations for component m4pi_both. */
  SIG_MESSAGE("m4pi_both (Init)");
#define mccompcurname  m4pi_both
#define mccompcurtype  PSD_monitor_4PI
#define mccompcurindex 17
#define nx mccm4pi_both_nx
#define ny mccm4pi_both_ny
#define PSD_N mccm4pi_both_PSD_N
#define PSD_p mccm4pi_both_PSD_p
#define PSD_p2 mccm4pi_both_PSD_p2
#define filename mccm4pi_both_filename
#define radius mccm4pi_both_radius
#define restore_neutron mccm4pi_both_restore_neutron
#define nowritefile mccm4pi_both_nowritefile
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
#line 10974 "./Union_test_mask.c"
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
  /* TRACE Component test_incoherent [1] */
  mccoordschange(mcposrtest_incoherent, mcrotrtest_incoherent,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component test_incoherent (without coords transformations) */
  mcJumpTrace_test_incoherent:
  SIG_MESSAGE("test_incoherent (Trace)");
  mcDEBUG_COMP("test_incoherent")
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

#define mcabsorbComp mcabsorbComptest_incoherent
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
#define mccompcurname  test_incoherent
#define mccompcurtype  Incoherent_process
#define mccompcurindex 1
#define This_process mcctest_incoherent_This_process
#define Incoherent_storage mcctest_incoherent_Incoherent_storage
#define effective_my_scattering mcctest_incoherent_effective_my_scattering
{   /* Declarations of test_incoherent=Incoherent_process() SETTING parameters. */
MCNUM sigma = mcctest_incoherent_sigma;
MCNUM f_QE = mcctest_incoherent_f_QE;
MCNUM gamma = mcctest_incoherent_gamma;
MCNUM packing_factor = mcctest_incoherent_packing_factor;
MCNUM unit_cell_volume = mcctest_incoherent_unit_cell_volume;
MCNUM interact_fraction = mcctest_incoherent_interact_fraction;
}   /* End of test_incoherent=Incoherent_process() SETTING parameter declarations. */
#undef effective_my_scattering
#undef Incoherent_storage
#undef This_process
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbComptest_incoherent:
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

  /* TRACE Component test_material [2] */
  mccoordschange(mcposrtest_material, mcrotrtest_material,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component test_material (without coords transformations) */
  mcJumpTrace_test_material:
  SIG_MESSAGE("test_material (Trace)");
  mcDEBUG_COMP("test_material")
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

#define mcabsorbComp mcabsorbComptest_material
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
#define mccompcurname  test_material
#define mccompcurtype  Union_make_material
#define mccompcurindex 2
#define loop_index mcctest_material_loop_index
#define this_material mcctest_material_this_material
#define accepted_processes mcctest_material_accepted_processes
#define global_material_element mcctest_material_global_material_element
{   /* Declarations of test_material=Union_make_material() SETTING parameters. */
char* process_string = mcctest_material_process_string;
MCNUM my_absorption = mcctest_material_my_absorption;
MCNUM absorber = mcctest_material_absorber;
}   /* End of test_material=Union_make_material() SETTING parameter declarations. */
#undef global_material_element
#undef accepted_processes
#undef this_material
#undef loop_index
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbComptest_material:
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

  /* TRACE Component a1 [3] */
  mccoordschange(mcposra1, mcrotra1,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component a1 (without coords transformations) */
  mcJumpTrace_a1:
  SIG_MESSAGE("a1 (Trace)");
  mcDEBUG_COMP("a1")
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

#define mcabsorbComp mcabsorbCompa1
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
#define mccompcurname  a1
#define mccompcurtype  Progress_bar
#define mccompcurindex 3
#define IntermediateCnts mcca1_IntermediateCnts
#define StartTime mcca1_StartTime
#define EndTime mcca1_EndTime
#define CurrentTime mcca1_CurrentTime
{   /* Declarations of a1=Progress_bar() SETTING parameters. */
char* profile = mcca1_profile;
MCNUM percent = mcca1_percent;
MCNUM flag_save = mcca1_flag_save;
MCNUM minutes = mcca1_minutes;
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
#line 11374 "./Union_test_mask.c"
}   /* End of a1=Progress_bar() SETTING parameter declarations. */
#undef CurrentTime
#undef EndTime
#undef StartTime
#undef IntermediateCnts
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompa1:
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

  /* TRACE Component source [4] */
  mccoordschange(mcposrsource, mcrotrsource,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component source (without coords transformations) */
  mcJumpTrace_source:
  SIG_MESSAGE("source (Trace)");
  mcDEBUG_COMP("source")
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

#define mcabsorbComp mcabsorbCompsource
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
#define mccompcurname  source
#define mccompcurtype  Source_div
#define mccompcurindex 4
#define thetah mccsource_thetah
#define thetav mccsource_thetav
#define sigmah mccsource_sigmah
#define sigmav mccsource_sigmav
#define tan_h mccsource_tan_h
#define tan_v mccsource_tan_v
#define p_init mccsource_p_init
#define dist mccsource_dist
#define focus_xw mccsource_focus_xw
#define focus_yh mccsource_focus_yh
{   /* Declarations of source=Source_div() SETTING parameters. */
MCNUM xwidth = mccsource_xwidth;
MCNUM yheight = mccsource_yheight;
MCNUM focus_aw = mccsource_focus_aw;
MCNUM focus_ah = mccsource_focus_ah;
MCNUM E0 = mccsource_E0;
MCNUM dE = mccsource_dE;
MCNUM lambda0 = mccsource_lambda0;
MCNUM dlambda = mccsource_dlambda;
MCNUM gauss = mccsource_gauss;
MCNUM flux = mccsource_flux;
#line 118 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../sources/Source_div.comp"
{
  double E,lambda,v;

  p=p_init;
  z=0;
  t=0;

  x=randpm1()*xwidth/2.0;
  y=randpm1()*yheight/2.0;
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

  if (gauss==1) {
      thetah = randnorm()*sigmah;
      thetav = randnorm()*sigmav;
  } else {
      /*find limits of uniform sampling scheme for vertical divergence.
        thetav should be acos(1-2*U) for U\in[0,1]. for theta measured from vertical axis
        we only use a sub-interval for U and measure from horizontal plane.*/
      double sample_lim1,u2;
      sample_lim1=(1-cos(M_PI_2 - focus_ah/2.0*DEG2RAD))*0.5;
      u2=randpm1()*(sample_lim1-0.5) + 0.5;
      thetav = acos(1-2*u2) - M_PI_2;
      thetah = randpm1()*focus_aw*DEG2RAD/2;
  }

  tan_h = tan(thetah);
  tan_v = tan(thetav);

  /* Perform the correct treatment - no small angle approx. here! */
  vz = v / sqrt(1 + tan_v*tan_v + tan_h*tan_h);
  vy = tan_v * vz;
  vx = tan_h * vz;
}
#line 11552 "./Union_test_mask.c"
}   /* End of source=Source_div() SETTING parameter declarations. */
#undef focus_yh
#undef focus_xw
#undef dist
#undef p_init
#undef tan_v
#undef tan_h
#undef sigmav
#undef sigmah
#undef thetav
#undef thetah
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompsource:
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

  /* TRACE Component beam_center [5] */
  mccoordschange(mcposrbeam_center, mcrotrbeam_center,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component beam_center (without coords transformations) */
  mcJumpTrace_beam_center:
  SIG_MESSAGE("beam_center (Trace)");
  mcDEBUG_COMP("beam_center")
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

#define mcabsorbComp mcabsorbCompbeam_center
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
#define mccompcurname  beam_center
#define mccompcurtype  Arm
#define mccompcurindex 5
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompbeam_center:
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

  /* TRACE Component volume_1 [6] */
  mccoordschange(mcposrvolume_1, mcrotrvolume_1,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component volume_1 (without coords transformations) */
  mcJumpTrace_volume_1:
  SIG_MESSAGE("volume_1 (Trace)");
  mcDEBUG_COMP("volume_1")
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

#define mcabsorbComp mcabsorbCompvolume_1
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
#define mccompcurname  volume_1
#define mccompcurtype  Union_box
#define mccompcurindex 6
#define loop_index mccvolume_1_loop_index
#define this_box_volume mccvolume_1_this_box_volume
#define global_geometry_element mccvolume_1_global_geometry_element
#define this_box_storage mccvolume_1_this_box_storage
{   /* Declarations of volume_1=Union_box() SETTING parameters. */
char* material_string = mccvolume_1_material_string;
MCNUM priority = mccvolume_1_priority;
MCNUM xwidth = mccvolume_1_xwidth;
MCNUM yheight = mccvolume_1_yheight;
MCNUM zdepth = mccvolume_1_zdepth;
MCNUM xwidth2 = mccvolume_1_xwidth2;
MCNUM yheight2 = mccvolume_1_yheight2;
MCNUM visualize = mccvolume_1_visualize;
int target_index = mccvolume_1_target_index;
MCNUM target_x = mccvolume_1_target_x;
MCNUM target_y = mccvolume_1_target_y;
MCNUM target_z = mccvolume_1_target_z;
MCNUM focus_aw = mccvolume_1_focus_aw;
MCNUM focus_ah = mccvolume_1_focus_ah;
MCNUM focus_xw = mccvolume_1_focus_xw;
MCNUM focus_xh = mccvolume_1_focus_xh;
MCNUM focus_r = mccvolume_1_focus_r;
MCNUM p_interact = mccvolume_1_p_interact;
char* mask_string = mccvolume_1_mask_string;
char* mask_setting = mccvolume_1_mask_setting;
MCNUM number_of_activations = mccvolume_1_number_of_activations;
}   /* End of volume_1=Union_box() SETTING parameter declarations. */
#undef this_box_storage
#undef global_geometry_element
#undef this_box_volume
#undef loop_index
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompvolume_1:
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

  /* TRACE Component volume_2 [7] */
  mccoordschange(mcposrvolume_2, mcrotrvolume_2,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component volume_2 (without coords transformations) */
  mcJumpTrace_volume_2:
  SIG_MESSAGE("volume_2 (Trace)");
  mcDEBUG_COMP("volume_2")
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

#define mcabsorbComp mcabsorbCompvolume_2
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
#define mccompcurname  volume_2
#define mccompcurtype  Union_box
#define mccompcurindex 7
#define loop_index mccvolume_2_loop_index
#define this_box_volume mccvolume_2_this_box_volume
#define global_geometry_element mccvolume_2_global_geometry_element
#define this_box_storage mccvolume_2_this_box_storage
{   /* Declarations of volume_2=Union_box() SETTING parameters. */
char* material_string = mccvolume_2_material_string;
MCNUM priority = mccvolume_2_priority;
MCNUM xwidth = mccvolume_2_xwidth;
MCNUM yheight = mccvolume_2_yheight;
MCNUM zdepth = mccvolume_2_zdepth;
MCNUM xwidth2 = mccvolume_2_xwidth2;
MCNUM yheight2 = mccvolume_2_yheight2;
MCNUM visualize = mccvolume_2_visualize;
int target_index = mccvolume_2_target_index;
MCNUM target_x = mccvolume_2_target_x;
MCNUM target_y = mccvolume_2_target_y;
MCNUM target_z = mccvolume_2_target_z;
MCNUM focus_aw = mccvolume_2_focus_aw;
MCNUM focus_ah = mccvolume_2_focus_ah;
MCNUM focus_xw = mccvolume_2_focus_xw;
MCNUM focus_xh = mccvolume_2_focus_xh;
MCNUM focus_r = mccvolume_2_focus_r;
MCNUM p_interact = mccvolume_2_p_interact;
char* mask_string = mccvolume_2_mask_string;
char* mask_setting = mccvolume_2_mask_setting;
MCNUM number_of_activations = mccvolume_2_number_of_activations;
}   /* End of volume_2=Union_box() SETTING parameter declarations. */
#undef this_box_storage
#undef global_geometry_element
#undef this_box_volume
#undef loop_index
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompvolume_2:
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

  /* TRACE Component mask_1_1 [8] */
  mccoordschange(mcposrmask_1_1, mcrotrmask_1_1,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component mask_1_1 (without coords transformations) */
  mcJumpTrace_mask_1_1:
  SIG_MESSAGE("mask_1_1 (Trace)");
  mcDEBUG_COMP("mask_1_1")
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

#define mcabsorbComp mcabsorbCompmask_1_1
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
#define mccompcurname  mask_1_1
#define mccompcurtype  Union_box
#define mccompcurindex 8
#define loop_index mccmask_1_1_loop_index
#define this_box_volume mccmask_1_1_this_box_volume
#define global_geometry_element mccmask_1_1_global_geometry_element
#define this_box_storage mccmask_1_1_this_box_storage
{   /* Declarations of mask_1_1=Union_box() SETTING parameters. */
char* material_string = mccmask_1_1_material_string;
MCNUM priority = mccmask_1_1_priority;
MCNUM xwidth = mccmask_1_1_xwidth;
MCNUM yheight = mccmask_1_1_yheight;
MCNUM zdepth = mccmask_1_1_zdepth;
MCNUM xwidth2 = mccmask_1_1_xwidth2;
MCNUM yheight2 = mccmask_1_1_yheight2;
MCNUM visualize = mccmask_1_1_visualize;
int target_index = mccmask_1_1_target_index;
MCNUM target_x = mccmask_1_1_target_x;
MCNUM target_y = mccmask_1_1_target_y;
MCNUM target_z = mccmask_1_1_target_z;
MCNUM focus_aw = mccmask_1_1_focus_aw;
MCNUM focus_ah = mccmask_1_1_focus_ah;
MCNUM focus_xw = mccmask_1_1_focus_xw;
MCNUM focus_xh = mccmask_1_1_focus_xh;
MCNUM focus_r = mccmask_1_1_focus_r;
MCNUM p_interact = mccmask_1_1_p_interact;
char* mask_string = mccmask_1_1_mask_string;
char* mask_setting = mccmask_1_1_mask_setting;
MCNUM number_of_activations = mccmask_1_1_number_of_activations;
}   /* End of mask_1_1=Union_box() SETTING parameter declarations. */
#undef this_box_storage
#undef global_geometry_element
#undef this_box_volume
#undef loop_index
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompmask_1_1:
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

  /* TRACE Component mask_1_2 [9] */
  mccoordschange(mcposrmask_1_2, mcrotrmask_1_2,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component mask_1_2 (without coords transformations) */
  mcJumpTrace_mask_1_2:
  SIG_MESSAGE("mask_1_2 (Trace)");
  mcDEBUG_COMP("mask_1_2")
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

#define mcabsorbComp mcabsorbCompmask_1_2
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
#define mccompcurname  mask_1_2
#define mccompcurtype  Union_box
#define mccompcurindex 9
#define loop_index mccmask_1_2_loop_index
#define this_box_volume mccmask_1_2_this_box_volume
#define global_geometry_element mccmask_1_2_global_geometry_element
#define this_box_storage mccmask_1_2_this_box_storage
{   /* Declarations of mask_1_2=Union_box() SETTING parameters. */
char* material_string = mccmask_1_2_material_string;
MCNUM priority = mccmask_1_2_priority;
MCNUM xwidth = mccmask_1_2_xwidth;
MCNUM yheight = mccmask_1_2_yheight;
MCNUM zdepth = mccmask_1_2_zdepth;
MCNUM xwidth2 = mccmask_1_2_xwidth2;
MCNUM yheight2 = mccmask_1_2_yheight2;
MCNUM visualize = mccmask_1_2_visualize;
int target_index = mccmask_1_2_target_index;
MCNUM target_x = mccmask_1_2_target_x;
MCNUM target_y = mccmask_1_2_target_y;
MCNUM target_z = mccmask_1_2_target_z;
MCNUM focus_aw = mccmask_1_2_focus_aw;
MCNUM focus_ah = mccmask_1_2_focus_ah;
MCNUM focus_xw = mccmask_1_2_focus_xw;
MCNUM focus_xh = mccmask_1_2_focus_xh;
MCNUM focus_r = mccmask_1_2_focus_r;
MCNUM p_interact = mccmask_1_2_p_interact;
char* mask_string = mccmask_1_2_mask_string;
char* mask_setting = mccmask_1_2_mask_setting;
MCNUM number_of_activations = mccmask_1_2_number_of_activations;
}   /* End of mask_1_2=Union_box() SETTING parameter declarations. */
#undef this_box_storage
#undef global_geometry_element
#undef this_box_volume
#undef loop_index
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompmask_1_2:
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

  /* TRACE Component mask_2_1 [10] */
  mccoordschange(mcposrmask_2_1, mcrotrmask_2_1,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component mask_2_1 (without coords transformations) */
  mcJumpTrace_mask_2_1:
  SIG_MESSAGE("mask_2_1 (Trace)");
  mcDEBUG_COMP("mask_2_1")
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

#define mcabsorbComp mcabsorbCompmask_2_1
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
#define mccompcurname  mask_2_1
#define mccompcurtype  Union_box
#define mccompcurindex 10
#define loop_index mccmask_2_1_loop_index
#define this_box_volume mccmask_2_1_this_box_volume
#define global_geometry_element mccmask_2_1_global_geometry_element
#define this_box_storage mccmask_2_1_this_box_storage
{   /* Declarations of mask_2_1=Union_box() SETTING parameters. */
char* material_string = mccmask_2_1_material_string;
MCNUM priority = mccmask_2_1_priority;
MCNUM xwidth = mccmask_2_1_xwidth;
MCNUM yheight = mccmask_2_1_yheight;
MCNUM zdepth = mccmask_2_1_zdepth;
MCNUM xwidth2 = mccmask_2_1_xwidth2;
MCNUM yheight2 = mccmask_2_1_yheight2;
MCNUM visualize = mccmask_2_1_visualize;
int target_index = mccmask_2_1_target_index;
MCNUM target_x = mccmask_2_1_target_x;
MCNUM target_y = mccmask_2_1_target_y;
MCNUM target_z = mccmask_2_1_target_z;
MCNUM focus_aw = mccmask_2_1_focus_aw;
MCNUM focus_ah = mccmask_2_1_focus_ah;
MCNUM focus_xw = mccmask_2_1_focus_xw;
MCNUM focus_xh = mccmask_2_1_focus_xh;
MCNUM focus_r = mccmask_2_1_focus_r;
MCNUM p_interact = mccmask_2_1_p_interact;
char* mask_string = mccmask_2_1_mask_string;
char* mask_setting = mccmask_2_1_mask_setting;
MCNUM number_of_activations = mccmask_2_1_number_of_activations;
}   /* End of mask_2_1=Union_box() SETTING parameter declarations. */
#undef this_box_storage
#undef global_geometry_element
#undef this_box_volume
#undef loop_index
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompmask_2_1:
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

  /* TRACE Component mask_2_2 [11] */
  mccoordschange(mcposrmask_2_2, mcrotrmask_2_2,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component mask_2_2 (without coords transformations) */
  mcJumpTrace_mask_2_2:
  SIG_MESSAGE("mask_2_2 (Trace)");
  mcDEBUG_COMP("mask_2_2")
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

#define mcabsorbComp mcabsorbCompmask_2_2
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
#define mccompcurname  mask_2_2
#define mccompcurtype  Union_box
#define mccompcurindex 11
#define loop_index mccmask_2_2_loop_index
#define this_box_volume mccmask_2_2_this_box_volume
#define global_geometry_element mccmask_2_2_global_geometry_element
#define this_box_storage mccmask_2_2_this_box_storage
{   /* Declarations of mask_2_2=Union_box() SETTING parameters. */
char* material_string = mccmask_2_2_material_string;
MCNUM priority = mccmask_2_2_priority;
MCNUM xwidth = mccmask_2_2_xwidth;
MCNUM yheight = mccmask_2_2_yheight;
MCNUM zdepth = mccmask_2_2_zdepth;
MCNUM xwidth2 = mccmask_2_2_xwidth2;
MCNUM yheight2 = mccmask_2_2_yheight2;
MCNUM visualize = mccmask_2_2_visualize;
int target_index = mccmask_2_2_target_index;
MCNUM target_x = mccmask_2_2_target_x;
MCNUM target_y = mccmask_2_2_target_y;
MCNUM target_z = mccmask_2_2_target_z;
MCNUM focus_aw = mccmask_2_2_focus_aw;
MCNUM focus_ah = mccmask_2_2_focus_ah;
MCNUM focus_xw = mccmask_2_2_focus_xw;
MCNUM focus_xh = mccmask_2_2_focus_xh;
MCNUM focus_r = mccmask_2_2_focus_r;
MCNUM p_interact = mccmask_2_2_p_interact;
char* mask_string = mccmask_2_2_mask_string;
char* mask_setting = mccmask_2_2_mask_setting;
MCNUM number_of_activations = mccmask_2_2_number_of_activations;
}   /* End of mask_2_2=Union_box() SETTING parameter declarations. */
#undef this_box_storage
#undef global_geometry_element
#undef this_box_volume
#undef loop_index
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompmask_2_2:
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

  /* TRACE Component test_sample [12] */
  mccoordschange(mcposrtest_sample, mcrotrtest_sample,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component test_sample (without coords transformations) */
  mcJumpTrace_test_sample:
  SIG_MESSAGE("test_sample (Trace)");
  mcDEBUG_COMP("test_sample")
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

#define mcabsorbComp mcabsorbComptest_sample
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
#define mccompcurname  test_sample
#define mccompcurtype  Union_master
#define mccompcurindex 12
#define verbal mcctest_sample_verbal
#define list_verbal mcctest_sample_list_verbal
#define trace_verbal mcctest_sample_trace_verbal
#define finally_verbal mcctest_sample_finally_verbal
#define starting_volume_warning mcctest_sample_starting_volume_warning
#define global_master_element mcctest_sample_global_master_element
#define this_global_master_index mcctest_sample_this_global_master_index
#define previous_master_index mcctest_sample_previous_master_index
#define geometry_list_index mcctest_sample_geometry_list_index
#define intersection_time_table mcctest_sample_intersection_time_table
#define Volumes mcctest_sample_Volumes
#define Geometries mcctest_sample_Geometries
#define starting_lists mcctest_sample_starting_lists
#define r mcctest_sample_r
#define r_start mcctest_sample_r_start
#define v mcctest_sample_v
#define error_msg mcctest_sample_error_msg
#define component_error_msg mcctest_sample_component_error_msg
#define string_output mcctest_sample_string_output
#define number_of_volumes mcctest_sample_number_of_volumes
#define volume_index mcctest_sample_volume_index
#define process_index mcctest_sample_process_index
#define solutions mcctest_sample_solutions
#define max_number_of_processes mcctest_sample_max_number_of_processes
#define limit mcctest_sample_limit
#define solution mcctest_sample_solution
#define min_solution mcctest_sample_min_solution
#define min_volume mcctest_sample_min_volume
#define time_found mcctest_sample_time_found
#define intersection_time mcctest_sample_intersection_time
#define min_intersection_time mcctest_sample_min_intersection_time
#define process mcctest_sample_process
#define process_start mcctest_sample_process_start
#define my_trace mcctest_sample_my_trace
#define p_my_trace mcctest_sample_p_my_trace
#define my_trace_fraction_control mcctest_sample_my_trace_fraction_control
#define k mcctest_sample_k
#define k_new mcctest_sample_k_new
#define k_old mcctest_sample_k_old
#define v_length mcctest_sample_v_length
#define my_sum mcctest_sample_my_sum
#define my_sum_plus_abs mcctest_sample_my_sum_plus_abs
#define culmative_probability mcctest_sample_culmative_probability
#define mc_prop mcctest_sample_mc_prop
#define time_to_scattering mcctest_sample_time_to_scattering
#define length_to_scattering mcctest_sample_length_to_scattering
#define length_to_boundery mcctest_sample_length_to_boundery
#define time_to_boundery mcctest_sample_time_to_boundery
#define selected_process mcctest_sample_selected_process
#define scattering_event mcctest_sample_scattering_event
#define time_propagated_without_scattering mcctest_sample_time_propagated_without_scattering
#define a_next_volume_found mcctest_sample_a_next_volume_found
#define next_volume mcctest_sample_next_volume
#define next_volume_priority mcctest_sample_next_volume_priority
#define done mcctest_sample_done
#define current_volume mcctest_sample_current_volume
#define number_of_solutions mcctest_sample_number_of_solutions
#define number_of_solutions_static mcctest_sample_number_of_solutions_static
#define check mcctest_sample_check
#define start mcctest_sample_start
#define intersection_with_children mcctest_sample_intersection_with_children
#define geometry_output mcctest_sample_geometry_output
#define tree_next_volume mcctest_sample_tree_next_volume
#define pre_allocated1 mcctest_sample_pre_allocated1
#define pre_allocated2 mcctest_sample_pre_allocated2
#define pre_allocated3 mcctest_sample_pre_allocated3
#define ray_position mcctest_sample_ray_position
#define ray_velocity mcctest_sample_ray_velocity
#define ray_velocity_final mcctest_sample_ray_velocity_final
#define volume_0_found mcctest_sample_volume_0_found
#define scattered_flag mcctest_sample_scattered_flag
#define scattered_flag_VP mcctest_sample_scattered_flag_VP
#define master_transposed_rotation_matrix mcctest_sample_master_transposed_rotation_matrix
#define temp_rotation_matrix mcctest_sample_temp_rotation_matrix
#define non_rotated_position mcctest_sample_non_rotated_position
#define rotated_position mcctest_sample_rotated_position
#define enable_tagging mcctest_sample_enable_tagging
#define stop_tagging_ray mcctest_sample_stop_tagging_ray
#define stop_creating_nodes mcctest_sample_stop_creating_nodes
#define enable_tagging_check mcctest_sample_enable_tagging_check
#define master_tagging_node_list mcctest_sample_master_tagging_node_list
#define current_tagging_node mcctest_sample_current_tagging_node
#define tagging_leaf_counter mcctest_sample_tagging_leaf_counter
#define number_of_scattering_events mcctest_sample_number_of_scattering_events
#define real_transmission_probability mcctest_sample_real_transmission_probability
#define mc_transmission_probability mcctest_sample_mc_transmission_probability
#define number_of_masks mcctest_sample_number_of_masks
#define number_of_masked_volumes mcctest_sample_number_of_masked_volumes
#define need_to_run_within_which_volume mcctest_sample_need_to_run_within_which_volume
#define mask_index_main mcctest_sample_mask_index_main
#define mask_iterate mcctest_sample_mask_iterate
#define mask_status_list mcctest_sample_mask_status_list
#define current_mask_intersect_list_status mcctest_sample_current_mask_intersect_list_status
#define mask_volume_index_list mcctest_sample_mask_volume_index_list
#define geometry_component_index_list mcctest_sample_geometry_component_index_list
#define Volume_copies_allocated mcctest_sample_Volume_copies_allocated
#define p_old mcctest_sample_p_old
#define this_logger mcctest_sample_this_logger
#define conditional_status mcctest_sample_conditional_status
#define tagging_conditional_list mcctest_sample_tagging_conditional_list
#define free_tagging_conditioanl_list mcctest_sample_free_tagging_conditioanl_list
#define logger_conditional_extend_array mcctest_sample_logger_conditional_extend_array
#define tagging_conditional_extend mcctest_sample_tagging_conditional_extend
#define max_conditional_extend_index mcctest_sample_max_conditional_extend_index
#define safty_distance mcctest_sample_safty_distance
#define safty_distance2 mcctest_sample_safty_distance2
#define number_of_processes_array mcctest_sample_number_of_processes_array
#define temporary_focus_data mcctest_sample_temporary_focus_data
#define focus_data_index mcctest_sample_focus_data_index
{   /* Declarations of test_sample=Union_master() SETTING parameters. */
MCNUM allow_inside_start = mcctest_sample_allow_inside_start;
MCNUM history_limit = mcctest_sample_history_limit;
MCNUM enable_conditionals = mcctest_sample_enable_conditionals;
MCNUM inherit_number_of_scattering_events = mcctest_sample_inherit_number_of_scattering_events;
#line 877 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../contrib/union/Union_master.comp"
{
  
  #ifdef Union_trace_verbal_setting
    printf("\n\n\n\n\n----------- NEW RAY -------------------------------------------------\n");
    printf("Union_master component name: %s \n \n",NAME_CURRENT_COMP);
  #endif
  
  // Initialize logic
  done = 0;
  error_msg = 0;
  clear_intersection_table(&intersection_time_table);
  
  time_propagated_without_scattering = 0;
  v_length = sqrt(vx*vx+vy*vy+vz*vz);
  
  // Initialize logger system / Statistics
  number_of_scattering_events = 0;
  
  if (inherit_number_of_scattering_events==1) // Continue number of scattering from previous Union_master
    number_of_scattering_events = global_master_list.elements[this_global_master_index-1].stored_number_of_scattering_events;
  
  // Zero scattered_flag_VP data
  for (volume_index = 1;volume_index<number_of_volumes;volume_index++) { // No reason to update volume 0, as scattering doesn't happen there
    scattered_flag[volume_index] = 0;
    for (process_index=0;process_index<number_of_processes_array[volume_index];process_index++)
      scattered_flag_VP[volume_index][process_index] = 0;
  }
  
  // If first Union_master in instrument, reset loggers_with_data_array and clean unused data.
  // Unused data happens when logging data is passed to the next Union_master, but the ray is absorbed on the way.
  // Could be improved by using the precompiler instead as ncount times the number of Union_masters could be avoided.
  if (global_master_list.elements[0].component_index == INDEX_CURRENT_COMP) {
    // If this is the first Union master, clean up logger data for rays that did not make it through Union components
    for (log_index=loggers_with_data_array.used_elements-1;log_index>-1;log_index--) {
      loggers_with_data_array.logger_pointers[log_index]->function_pointers.clear_temp(&loggers_with_data_array.logger_pointers[log_index]->data_union);
    }
    loggers_with_data_array.used_elements = 0;
  }
    
  tagging_conditional_extend = 0;
  for (iterate=0;iterate<max_conditional_extend_index+1;iterate++) {
    logger_conditional_extend_array[iterate] = 0;
  }

  // Need to clean up the double notation for position and velocity. // REVIEW_LINE
  r_start[0] = x; r_start[1] = y; r_start[2] = z;
  r[0]=x;r[1]=y;r[2]=z;v[0]=vx;v[1]=vy;v[2]=vz; // REVIEW_LINE r and v are bad names
  
  ray_position = coords_set(x,y,z);
  ray_velocity = coords_set(vx,vy,vz);
  
  // Mask update: need to check the mask status for the initial position
  // mask status for a mask is 1 if the ray position is inside, 0 if it is outside
  for (iterate=0;iterate<number_of_masks;iterate++) {
    if(Volumes[mask_volume_index_list.elements[iterate]]->geometry.within_function(ray_position,&Volumes[mask_volume_index_list.elements[iterate]]->geometry) == 1) {
      mask_status_list.elements[iterate] = 1;
    } else {
      mask_status_list.elements[iterate] = 0;
    }
  }
  
  #ifdef Union_trace_verbal_setting
    print_1d_int_list(mask_status_list,"Initial mask status list");
  #endif
  
  // Now the initial current_volume can be found, which requires the up to date mask_status_list
  current_volume = within_which_volume(ray_position,starting_lists.reduced_start_list,starting_lists.starting_destinations_list,Volumes,&mask_status_list,number_of_volumes,pre_allocated1,pre_allocated2,pre_allocated3);
  
  // Using the mask_status_list and the current volume, the current_mask_intersect_list_status can be made
  //  it contains the effective mask status of all volumes on the current volumes mask intersect list, which needs to be calculated,
  //  but only when the current volume or mask status changes, not under for example scattering inside the current volume
  update_current_mask_intersect_status(&current_mask_intersect_list_status, &mask_status_list, Volumes, &current_volume);
  
  #ifdef Union_trace_verbal_setting
    printf("Starting current_volume = %d\n",current_volume);
  #endif
  
  // Check if the ray appeared in an allowed starting volume, unless this check is disabled by the user for advanced cases
  if (allow_inside_start == 0 && starting_lists.allowed_starting_volume_logic_list.elements[current_volume] == 0) {
    printf("ERROR, ray ''teleported'' into Union component %s, if intentional, set allow_inside_start=1\n",NAME_CURRENT_COMP);
    exit(1);
  }
  // Warn the user that rays have appeared inside a volume instead of outside as expected
  if (starting_volume_warning == 0 && current_volume != 0) {
        printf("WARNING: Ray started in volume ''%s'' rather than the surrounding vacuum in component %s. This warning is only shown once.\n",Volumes[current_volume]->name,NAME_CURRENT_COMP);
        starting_volume_warning = 1;
  }
  
  // Placing the new ray at the start of the tagging tree corresponding to current volume
  // A history limit can be imposed so that no new nodes are created after this limit (may be necessary to fit in memory)
  // Rays can still follow the nodes created before even when no additional nodes are created, but if a situation that
  //  requires a new node is encountered, stop_tagging_ray is set to 1, stopping further tagging and preventing the data
  //  for that ray to be used further.
  if (enable_tagging) {
    current_tagging_node = master_tagging_node_list.elements[current_volume];
    stop_tagging_ray = 0; // Allow this ray to be tracked
    if (tagging_leaf_counter > history_limit) stop_creating_nodes = 1;
  }
  
  #ifdef Union_trace_verbal_setting
    if (enable_tagging) printf("current_tagging_node->intensity = %f\n",current_tagging_node->intensity);
    if (enable_tagging) printf("current_tagging_node->number_of_rays = %d \n",current_tagging_node->number_of_rays);
  #endif
  
  // Propagation loop including scattering
  // This while loop continues until the ray leaves the ensamble of user defined volumes either through volume 0
  //  or a dedicated exit volume. The loop is cancelled after a large number of iterations as a failsafe for errors.
  // A single run of the loop will either be a propagation to the next volume along the path of the ray, or a
  //  scattering event at some point along the path of the ray in the current volume.
  limit = 100000;
  while (done == 0) {
    limit--;
    
    #ifdef Union_trace_verbal_setting
      printf("----------- START OF WHILE LOOP --------------------------------------\n");
      print_intersection_table(intersection_time_table);
      printf("current_volume = %d \n",current_volume);
    #endif
    
    // Calculating intersections with the necessary volumes. The relevant set of volumes depend on the current volume and the mask status array.
    // First the volumes on the current volumes intersect list is checked, then its mask interset list. Before checking the volume itself, it is
    //  checked if any children of the current volume is intersected, in which case the intersection calculation with the current volume can be
    //  skipped.
    
    // Checking intersections for all volumes in the intersect list.
    for (start=check=Volumes[current_volume]->geometry.intersect_check_list.elements;check-start<Volumes[current_volume]->geometry.intersect_check_list.num_elements;check++) {
    // This will leave check as a pointer to the intergers in the intersect_check_list and iccrement nicely
        #ifdef Union_trace_verbal_setting
          printf("Intersect_list = %d being checked \n",*check);
        #endif
    
        if (intersection_time_table.calculated[*check] == 0) {
             #ifdef Union_trace_verbal_setting
               printf("running intersection for intersect_list with *check = %d \n",*check);
               // if (trace_verbal) printf("r = (%f,%f,%f) v = (%f,%f,%f) \n",r[0],r[1],r[2],v[0],v[1],v[2]);
             #endif
            // Calculate intersections using intersect function imbedded in the relevant volume structure using parameters
            //  that are also imbedded in the structure.
            geometry_output = Volumes[*check]->geometry.intersect_function(intersection_time_table.intersection_times[*check],number_of_solutions,r_start,v,&Volumes[*check]->geometry);
            intersection_time_table.calculated[*check] = 1;
            // if (trace_verbal) printf("succesfully calculated intersection times for volume *check = %d \n",*check);
        }
    }
    
    // Mask update: add additional loop for checking intersections with masked volumes depending on mask statuses
    for (mask_iterate=0;mask_iterate<Volumes[current_volume]->geometry.mask_intersect_list.num_elements;mask_iterate++) {
      if (current_mask_intersect_list_status.elements[mask_iterate] == 1) { // Only check if the mask is active
        #ifdef Union_trace_verbal_setting
          printf("Mask Intersect_list = %d being checked \n",Volumes[current_volume]->geometry.mask_intersect_list.elements[mask_iterate]);
        #endif
        if (intersection_time_table.calculated[Volumes[current_volume]->geometry.mask_intersect_list.elements[mask_iterate]] == 0) {
          #ifdef Union_trace_verbal_setting
            printf("running intersection for mask_intersect_list element = %d \n",Volumes[current_volume]->geometry.mask_intersect_list.elements[mask_iterate]);
          //  printf("r = (%f,%f,%f) v = (%f,%f,%f) \n",r[0],r[1],r[2],v[0],v[1],v[2]);
          #endif
          // Calculate intersections using intersect function imbedded in the relevant volume structure using parameters
            //  that are also imbedded in the structure.
          geometry_output = Volumes[Volumes[current_volume]->geometry.mask_intersect_list.elements[mask_iterate]]->geometry.intersect_function(intersection_time_table.intersection_times[Volumes[current_volume]->geometry.mask_intersect_list.elements[mask_iterate]],number_of_solutions,r_start,v,&Volumes[Volumes[current_volume]->geometry.mask_intersect_list.elements[mask_iterate]]->geometry);
          intersection_time_table.calculated[Volumes[current_volume]->geometry.mask_intersect_list.elements[mask_iterate]] = 1;
          // if (trace_verbal) printf("succesfully calculated intersection times for volume *check = %d \n",*check);
        }
      }
    }
    
    // Checking if there are intersections with children of current volume, which means there is an intersection before current_volume, and thus can be skipped. But only if they have not been overwritten. In case current_volume is 0, there is no need to do this regardless.
    if (current_volume != 0 && intersection_time_table.calculated[current_volume] == 0) {
        #ifdef Union_trace_verbal_setting
          printf("Checking if children of current_volume = %d have intersections. \n",current_volume);
        #endif
        intersection_with_children = 0;
        //for (start = check = Volumes[current_volume]->geometry.direct_children.elements;check - start < Volumes[current_volume]->geometry.children.num_elements;check++) { // REVIEW LINE. Caused bug with masks.
        for (start = check = Volumes[current_volume]->geometry.children.elements;check - start < Volumes[current_volume]->geometry.children.num_elements;check++) {
            #ifdef Union_trace_verbal_setting
              printf("Checking if child %d of current_volume = %d have intersections. \n",*check,current_volume);
            #endif
            // Only check the first of the two results in the intersection table, as they are ordered, and the second is of no interest
            if (intersection_time_table.calculated[*check] == 1 && intersection_time_table.intersection_times[*check][0] > time_propagated_without_scattering) {
                // If this child is masked, its mask status need to be 1 in order to be taken into account
                if (Volumes[*check]->geometry.is_masked_volume == 0) {
                  #ifdef Union_trace_verbal_setting
                    printf("Found an child of current_volume with an intersection. Skips calculating for current_volume \n");
                  #endif
                  intersection_with_children = 1;
                  break; // No need to check more, if there is just one it is not necessary to calculate intersection with current_volume yet
                } else {
                  #ifdef Union_trace_verbal_setting
                    printf("Found an child of current_volume with an intersection, but it is masked. Check to see if it can skip calculating for current_volume \n");
                  #endif
                  
                  if (Volumes[*check]->geometry.mask_mode == 2) { // ANY mask mode
                    tree_next_volume = 0;
                    for (mask_start=mask_check=Volumes[*check]->geometry.masked_by_mask_index_list.elements;mask_check-mask_start<Volumes[*check]->geometry.masked_by_mask_index_list.num_elements;mask_check++) {
                       if (mask_status_list.elements[*mask_check] == 1) {
                         intersection_with_children = 1;
                         break;
                       }
                    }
                  } else { // ALL mask mode
                    intersection_with_children = 1;
                    for (mask_start=mask_check=Volumes[*check]->geometry.masked_by_mask_index_list.elements;mask_check-mask_start<Volumes[*check]->geometry.masked_by_mask_index_list.num_elements;mask_check++) {
                      if (mask_status_list.elements[*mask_check] == 0) {
                        intersection_with_children = 0;
                        break;
                      }
                    }
                  }
                  #ifdef Union_trace_verbal_setting
                    printf("The mask status was 1, can actually skip intersection calculation for current volume \n");
                  #endif
                  if (intersection_with_children == 1) break;
                }
            }
        }
        #ifdef Union_trace_verbal_setting
          printf("intersection_with_children = %d \n",intersection_with_children);
        #endif
        if (intersection_with_children == 0) {
            geometry_output = Volumes[current_volume]->geometry.intersect_function(intersection_time_table.intersection_times[current_volume],number_of_solutions,r_start,v,&Volumes[current_volume]->geometry);
            intersection_time_table.calculated[current_volume] = 1;
        }
    }

    // At this point, intersection_time_table is updated with intersection times of all possible intersections.
    #ifdef Union_trace_verbal_setting
      print_intersection_table(intersection_time_table);
    #endif
    
    // Next task is to find the next intersection time. The next intersection must be greater than the time_propagated_without_scattering (0 at start of loop)
    // Loops are eqvialent to the 3 intersection calculation loops already completed
    
    // First loop for checking intersect_check_list
    time_found = 0;
    for (start=check=Volumes[current_volume]->geometry.intersect_check_list.elements;check-start<Volumes[current_volume]->geometry.intersect_check_list.num_elements;check++) {
        for (solution = 0;solution<intersection_time_table.n_elements[*check];solution++) {
            if (time_found) {
                if ((intersection_time = intersection_time_table.intersection_times[*check][solution]) > time_propagated_without_scattering &&  intersection_time < min_intersection_time) {
                    min_intersection_time = intersection_time;min_solution = solution;min_volume = *check;
                }
            } else {
                if ((intersection_time = intersection_time_table.intersection_times[*check][solution]) > time_propagated_without_scattering) {
                    min_intersection_time = intersection_time;min_solution = solution;min_volume = *check;
                    time_found = 1;
                }
            }
        }
    }
    
    // Now check the masked_intersect_list, but only the ones that are currently active
    for (mask_iterate=0;mask_iterate<Volumes[current_volume]->geometry.mask_intersect_list.num_elements;mask_iterate++) {
      if (current_mask_intersect_list_status.elements[mask_iterate] == 1) {
        for (solution = 0;solution<intersection_time_table.n_elements[Volumes[current_volume]->geometry.mask_intersect_list.elements[mask_iterate]];solution++) {
            if (time_found) {
                if ((intersection_time = intersection_time_table.intersection_times[Volumes[current_volume]->geometry.mask_intersect_list.elements[mask_iterate]][solution]) > time_propagated_without_scattering &&  intersection_time < min_intersection_time) {
                    min_intersection_time = intersection_time;min_solution = solution;min_volume = Volumes[current_volume]->geometry.mask_intersect_list.elements[mask_iterate];
                }
            } else {
                if ((intersection_time = intersection_time_table.intersection_times[Volumes[current_volume]->geometry.mask_intersect_list.elements[mask_iterate]][solution]) > time_propagated_without_scattering) {
                    min_intersection_time = intersection_time;min_solution = solution;min_volume = Volumes[current_volume]->geometry.mask_intersect_list.elements[mask_iterate];
                    time_found = 1;
                }
            }
        }
      }
    }
    
    // And check the current_volume
    for (solution = 0;solution<intersection_time_table.n_elements[current_volume];solution++) {
        if (time_found) {
            if ((intersection_time = intersection_time_table.intersection_times[current_volume][solution]) > time_propagated_without_scattering && intersection_time < min_intersection_time) {
                min_intersection_time = intersection_time;min_solution = solution;min_volume = current_volume;
            }
        } else {
            if ((intersection_time = intersection_time_table.intersection_times[current_volume][solution]) > time_propagated_without_scattering) {
                min_intersection_time = intersection_time;min_solution = solution;min_volume = current_volume;
                time_found = 1;
            }
        }
    }
    
    #ifdef Union_trace_verbal_setting
      printf("min_intersection_time = %f \n",min_intersection_time);
      printf("min_solution          = %d \n",min_solution);
      printf("min_volume            = %d \n",min_volume);
      printf("time_found            = %d \n",time_found);
    #endif

    // If a time is found, propagation continues, and it will be checked if a scattering occurs before the next intersection.
    // If a time is not found, the ray must be leaving the ensamble of volumes and the loop will be concluded
    if (time_found) {
        time_to_boundery = min_intersection_time - time_propagated_without_scattering; // calculate the time remaining before the next intersection
        scattering_event = 0; // Assume a scattering event will not occur
        
        // Check if a scattering event should occur
        if (current_volume != 0) { // Volume 0 is always vacuum, and if this is the current volume, an event will not occur
          if (Volumes[current_volume]->p_physics->number_of_processes == 0) { // If there are no processes, the volume could be vacuum or an absorber
            if (Volumes[current_volume]->p_physics->is_vacuum == 0)
              // This volume does not have physical processes but does have an absorption cross section, so the ray weight is reduced accordingly
              p *= exp(-Volumes[current_volume]->p_physics->my_a*2200*time_to_boundery);
    
    
              //#ifdef Union_trace_verbal_setting
              //printf("name of material: %s \n",Volumes[current_volume]->name);
              //printf("length to boundery  = %f\n",length_to_boundery);
              //printf("absorption cross section  = %f\n",Volumes[current_volume]->p_physics->my_a);
              //printf("chance to get through this length of absorber: %f \%\n",100*exp(-Volumes[current_volume]->p_physics->my_a*length_to_boundery));
              //#endif
                    
          } else {
            // Since there is a non-zero number of processes in this material, all the scattering cross section for these are calculated
            my_sum = 0; k[0] = V2K*vx; k[1] = V2K*vy; k[2] = V2K*vz; p_my_trace = my_trace; wavevector = coords_set(k[0],k[1],k[2]);
            for (process_start = process = Volumes[current_volume]->p_physics->p_scattering_array;process - process_start < Volumes[current_volume]->p_physics->number_of_processes;process++) {
            
              if (Volumes[current_volume]->p_physics->p_scattering_array[process - process_start].non_isotropic_rot_index != -1) {
                // If the process is not isotropic, the wavevector is transformed into the local coordinate system of the process
                wavevector_rotated = rot_apply(Volumes[current_volume]->geometry.process_rot_matrix_array[Volumes[current_volume]->p_physics->p_scattering_array[process - process_start].non_isotropic_rot_index],wavevector);
                
                coords_get(wavevector_rotated,&k_rotated[0],&k_rotated[1],&k_rotated[2]);

              } else {
                k_rotated[0] = k[0]; k_rotated[1] = k[1]; k_rotated[2] = k[2];
              }
              
              // Find correct focus_data_array index for this volume/process
              focus_data_index = Volumes[current_volume]->geometry.focus_array_indices.elements[process - process_start];
              
              // Call the probability for scattering function assighed to this specific procress (the process pointer is updated in the for loop)
              process->probability_for_scattering_function(p_my_trace,k_rotated,process->data_transfer,&Volumes[current_volume]->geometry.focus_data_array.elements[focus_data_index]);
              
              my_sum += *p_my_trace;
              #ifdef Union_trace_verbal_setting
                printf("my_trace = %f, my_sum = %f\n",*p_my_trace,my_sum);
              #endif
              p_my_trace++; // increment the pointer so that it point to the next element (max number of process in any material is allocated)
            }
            
            #ifdef Union_trace_verbal_setting
              printf("time_propagated_without_scattering = %f.\n",time_propagated_without_scattering);
              printf("v_length                           = %f.\n",v_length);
            #endif
            
            length_to_boundery = time_to_boundery * v_length;
            
            #ifdef Union_trace_verbal_setting
              printf("exp(- length_to_boundery*my_sum) = %f. length_to_boundery = %f. my_sum = %f.\n",exp(-length_to_boundery*my_sum),length_to_boundery,my_sum);
            #endif
            
            // Selecting if a scattering takes place, and what scattering process.
            // This section have too many if statements, and unessecary calculations
            // Could make seperate functions for p_interact on/off and interact_fraction on/off,
            //   and set function pointers to these in initialize, thus avoiding many unessecary if statements and calculations of x/x.
            
            my_sum_plus_abs = my_sum + Volumes[current_volume]->p_physics->my_a*(2200/v_length);
            
            if (my_sum < 1E-18) {
                // The scattering cross section is basicly zero, no scattering should occur.
                scattering_event = 0;
                p *= exp(-length_to_boundery*my_sum_plus_abs); // Correct for absorption and the almost zero scattering
            } else if (length_to_boundery < safty_distance2) {
                // Too close to boundery to safly make another scattering, attenuate
                p *= exp(-length_to_boundery*my_sum_plus_abs); // Attentuate the beam for the small distance
                scattering_event = 0;
            } else {
                // The scattering cross section is above zero and the distance to the boundery is sufficient for a scattering
                if (Volumes[current_volume]->geometry.geometry_p_interact != 0) {
                    // a fraction of the beam (geometry_p_interact) is forced to scatter
                    real_transmission_probability = exp(-length_to_boundery*my_sum_plus_abs);
                    mc_transmission_probability = (1.0 - Volumes[current_volume]->geometry.geometry_p_interact);
                    if ((scattering_event = (rand01() > mc_transmission_probability))) {
                        // Scattering event happens, this is the correction for the weight
                        p *= (1.0-real_transmission_probability)/(1.0-mc_transmission_probability); // Absorption simulated in weight
                        // Find length to next scattering knowing the ray will scatter.
                        length_to_scattering = safty_distance -log(1.0 - rand0max((1.0 - exp(-my_sum_plus_abs*(length_to_boundery-safty_distance2))))) / my_sum_plus_abs;
                    } else {
                        // Scattering event does not happen, this is the appropriate correction
                        p *= real_transmission_probability/mc_transmission_probability; // Absorption simulated in weight
                    }
                } else {
                    // probability to scatter is the natural value
                    if(my_sum*length_to_boundery < 1e-6) { // Scattering probability very small, linear method is used as exponential is unreliable
                      if (length_to_boundery > safty_distance2) {
                        if (rand01() < exp(-length_to_boundery*my_sum_plus_abs)) {
                          // Scattering happens, use linear description to select scattering position
                          length_to_scattering = safty_distance + rand0max(length_to_boundery - safty_distance2);
                          // Weight factor necessary to correct for using the linear scattering position distribution
                          p *= length_to_boundery*my_sum*exp(-length_to_scattering*my_sum_plus_abs); // Absorption simulated in weight
                          scattering_event = 1;
                        } else scattering_event = 0;
                      } else {
                        // The distance is too short to reliably make a scattering event (in comparison to accuraccy of intersect functions)
                        p *= exp(-length_to_boundery*my_sum_plus_abs); // Attentuate the beam for the small distance
                        scattering_event = 0;
                      }
                    } else {
                        // Strong scattering, use exponential description to select scattering position between safetydistance and infinity
                        length_to_scattering = safty_distance -log(1 - rand01() ) / my_sum_plus_abs;
                        // Scattering happens if the scattering position is before the boundery (and safty distance)
                        if (length_to_scattering < length_to_boundery - safty_distance) scattering_event = 1;
                        else scattering_event = 0;
                    }
                }
                
                if (scattering_event == 1) {
                  // Adjust weight for absorption
                  p *= my_sum/my_sum_plus_abs;
                  // Safety feature, alert in case of nonsense my results / negative absorption
                  if (my_sum/my_sum_plus_abs > 1.0) printf("WARNING: Absorption weight factor above 1! Should not happen! \n");
                  // Select process
                  if (Volumes[current_volume]->p_physics->number_of_processes == 1) { // trivial case
                    // Select the only available process, which will always have index 0
                    selected_process = 0;
                  } else {
                    if (Volumes[current_volume]->p_physics->interact_control == 1) {
                      // Interact_fraction is used to influence the choice of process in this material
                      mc_prop = rand01();culmative_probability=0;total_process_interact=1.0;
                  
                      // If any of the processes have probability 0, they are excluded from the selection
                      for (iterate = 0;iterate < Volumes[current_volume]->p_physics->number_of_processes;iterate++) {
                        if (my_trace[iterate] < 1E-18) {
                          // When this happens, the total force probability is corrected and the probability for this particular instance is set to 0
                          total_process_interact -= Volumes[current_volume]->p_physics->p_scattering_array[iterate].process_p_interact;
                          my_trace_fraction_control[iterate] = 0;
                          // In cases where my_trace is not zero, the forced fraction is still used.
                        } else my_trace_fraction_control[iterate] = Volumes[current_volume]->p_physics->p_scattering_array[iterate].process_p_interact;
                      }
                      // Randomly select a process using the weights stored in my_trace_fraction_control divided by total_process_interact
                      for (iterate = 0;iterate < Volumes[current_volume]->p_physics->number_of_processes;iterate++) {
                        culmative_probability += my_trace_fraction_control[iterate]/total_process_interact;
                        if (culmative_probability > mc_prop) {
                          selected_process = iterate;
                          p *= (my_trace[iterate]/my_sum)*(total_process_interact/my_trace_fraction_control[iterate]);
                          break;
                        }
                      }
                    
                    } else {
                      // Select a process based on their relative attenuations factors
                      mc_prop = rand01();culmative_probability=0;
                      for (iterate = 0;iterate < Volumes[current_volume]->p_physics->number_of_processes;iterate++) {
                        culmative_probability += my_trace[iterate]/my_sum;
                        if (culmative_probability > mc_prop) {
                          selected_process = iterate;
                          break;
                        }
                      }
                    }
                  }
                } // end of select process
            }
            
          }
        } // Done checking for scttering event and in case of scattering selecting a process
     
        if (scattering_event == 1) {
            #ifdef Union_trace_verbal_setting
              printf("SCATTERING EVENT \n");
              printf("current_volume            = %d \n",current_volume);
              printf("r = (%f,%f,%f) v = (%f,%f,%f) \n",r[0],r[1],r[2],v[0],v[1],v[2]);
            // printf("did scatter: my_trace = %E = %f \n",my_trace[selected_process],my_trace[selected_process]);
            #endif
            
            // Calculate the time to scattering
            time_to_scattering = length_to_scattering/v_length;
            
            #ifdef Union_trace_verbal_setting
              printf("time to scattering        = %2.20f \n",time_to_scattering);
            #endif
            
            //#ifdef Union_trace_verbal_setting
            //printf("length to boundery = %f, length to scattering = %f \n",length_to_boundery,length_to_scattering);
            //#endif
            
            //PROP_DT(time_to_scattering); // May be replace by version without gravity
            
            // Reduce the double book keeping done here // REVIEW LINE
            x += time_to_scattering*vx; y += time_to_scattering*vy; z += time_to_scattering*vz; t += time_to_scattering;
            r_start[0] = x; r_start[1] = y; r_start[2] = z;
            r[0] = x; r[1] = y; r[2] = z;
            ray_position = coords_set(x,y,z);
            ray_velocity = coords_set(vx,vy,vz);
            
            // Safe check that should be unecessary. Used to fine tune how close to the edge of a volume a scattering event is allowed to take place (1E-14 m away currently).
            if (Volumes[current_volume]->geometry.within_function(ray_position,&Volumes[current_volume]->geometry) == 0) {
              printf("\nERROR, propagated out of current volume instead of to a point within!\n");
              printf("length_to_scattering_specified = %2.20f\n             length propagated = %2.20f\n            length_to_boundery = %2.20f \n   current_position = (%lf,%lf,%lf) \n",length_to_scattering,sqrt(time_to_scattering*time_to_scattering*vx*vx+time_to_scattering*time_to_scattering*vy*vy+time_to_scattering*time_to_scattering*vz*vz),length_to_boundery,x,y,z);
              
              volume_index = within_which_volume(ray_position,starting_lists.reduced_start_list,starting_lists.starting_destinations_list,Volumes,&mask_status_list,number_of_volumes,pre_allocated1,pre_allocated2,pre_allocated3);
              
              printf("Debug info: Volumes[current_volume]->name = %s, but now inside volume number %d named %s.\n",Volumes[current_volume]->name,volume_index,Volumes[volume_index]->name);
              printf("Ray absorbed \n");
              ABSORB;
            }
            
            // Save information before scattering event needed in logging section
            p_old = p;
            k_old[0] = k[0];k_old[1] = k[1];k_old[2] = k[2];
            
            // Find correct focus_data_array index for this volume/process
            focus_data_index = Volumes[current_volume]->geometry.focus_array_indices.elements[selected_process];
        
            // Rotation to local process coordinate system (for non isotropic processes)
            if (Volumes[current_volume]->p_physics->p_scattering_array[selected_process].non_isotropic_rot_index != -1) {
                ray_velocity_rotated = rot_apply(Volumes[current_volume]->geometry.process_rot_matrix_array[Volumes[current_volume]->p_physics->p_scattering_array[selected_process].non_isotropic_rot_index],ray_velocity);
            } else {
                ray_velocity_rotated = ray_velocity;
            }
                
            // test_physics_scattering(double *k_final, double *k_initial, union data_transfer_union data_transfer) {
            //k[0] = V2K*ray_velocity.x; k[1] = V2K*ray_velocity.y; k[2] = V2K*ray_velocity.z;
            coords_get(coords_scalar_mult(ray_velocity_rotated,V2K),&k[0],&k[1],&k[2]);
            
            // I may replace a intial and final k with one instance that serves as both input and output
            if (0 == Volumes[current_volume]->p_physics->p_scattering_array[selected_process].scattering_function(k_new,k,&p,Volumes[current_volume]->p_physics->p_scattering_array[selected_process].data_transfer,&Volumes[current_volume]->geometry.focus_data_array.elements[0])) {
              /*
              // PowderN and Single_crystal requires the option of absorbing the neutron, which is weird. If there is a scattering probability, there should be a new direction.
              // It can arise from need to simplify sampling process and end up in cases where weight factor is 0, and the ray should be absorbed in these cases
              printf("ERROR: Union_master: %s.Absorbed ray because scattering function returned 0 (error/absorb)\n",NAME_CURRENT_COMP);
              component_error_msg++;
              if (component_error_msg > 100) {
                printf("To many errors encountered, exiting. \n");
                exit(1);
              }
              */
              ABSORB;
            }
            
            // Update velocity using k
            ray_velocity_rotated = coords_set(K2V*k_new[0],K2V*k_new[1],K2V*k_new[2]);
            
            // Transformation back to main coordinate system (maybe one should only do this when multiple scattering in that volume was over, especially if there is only one non isotropic frame)
            if (Volumes[current_volume]->p_physics->p_scattering_array[selected_process].non_isotropic_rot_index != -1) {
                ray_velocity_final = rot_apply(Volumes[current_volume]->geometry.transpose_process_rot_matrix_array[Volumes[current_volume]->p_physics->p_scattering_array[selected_process].non_isotropic_rot_index],ray_velocity_rotated);
            } else {
               ray_velocity_final = ray_velocity_rotated;
            }
            
            // Write velocity to global variable (temp, only really necessary at final)
            //vx = ray_velocity.x; vy = ray_velocity.y; vz = ray_velocity.z;
            coords_get(ray_velocity_final,&vx,&vy,&vz);
            
            // Write velocity in array format as it is still used by intersect functions (temp, they need to be updated to ray_position / ray_velocity)
            v[0] = vx; v[1] = vy; v[2] = vz;
            v_length = sqrt(vx*vx+vy*vy+vz*vz);
            k_new[0] = V2K*vx; k_new[1] = V2K*vy; k_new[2] = V2K*vz;
            if (verbal) if (v_length < 1) printf("velocity set to less than 1\n");
            
            #ifdef Union_trace_verbal_setting
              printf("Running logger system for specific volumes \n");
            #endif
            // Logging for detector components assosiated with this volume
            for (log_index=0;log_index<Volumes[current_volume]->loggers.num_elements;log_index++) {
              //printf("logging time! Volume specific version. Volume name = %s \n",Volumes[current_volume]->name);
              //printf("  log_index = %d \n",log_index);
              if (Volumes[current_volume]->loggers.p_logger_volume[log_index].p_logger_process[selected_process] != NULL) {
                // Technically the scattering function could edit k, the wavevector before the scattering, even though there would be little point to doing that.
                // Could save a secure copy and pass that instead to be certain that no scattering process accidently tampers with the logging.
                // printf("  the logging function pointer was not NULL \n");
                // PV (number of time scattered in this volume/process combination is not recorded. Need to expand scattering_flag to contain 2D, volume and process
                
                // This function calls a logger function which in turn stores some data among the passed, and possibly performs some basic data analysis
                Volumes[current_volume]->loggers.p_logger_volume[log_index].p_logger_process[selected_process]->function_pointers.active_record_function(&ray_position,k_new,k_old,p,p_old,t,scattered_flag[current_volume],scattered_flag_VP[current_volume][selected_process],number_of_scattering_events,Volumes[current_volume]->loggers.p_logger_volume[log_index].p_logger_process[selected_process],&loggers_with_data_array);
                // If the logging component have a conditional attatched, the collected data will be written to a temporary place
                // At the end of the rays life, it will be checked if the condition is met
                //  if it is met, the temporary data is transfered to permanent, and temp is cleared.
                //  if it is not met, the temporary data is cleared.
              }
            }
            
            #ifdef Union_trace_verbal_setting
              printf("Running logger system for all volumes \n");
            #endif
            for (log_index=0;log_index<global_all_volume_logger_list.num_elements;log_index++) {
              //printf("logging time! Global version. log_index = %d \n",log_index);
              // As above, but on a global scale, meaning scattering in all volumes are logged
              
              // Problems with VN, PV, as there is no assosiated volume or process. The functions however need to have the same input to make the logger components general.
              // Could be interesting to have a monitor that just globally measurres the second scattering event in any volume (must be two in the same). Weird but not meaningless.
              global_all_volume_logger_list.elements[log_index].logger->function_pointers.active_record_function(&ray_position,k_new,k_old,p,p_old,t,scattered_flag[current_volume],scattered_flag_VP[current_volume][selected_process],number_of_scattering_events,global_all_volume_logger_list.elements[log_index].logger,&loggers_with_data_array);
            }
            
            
            SCATTER;
            ++number_of_scattering_events;
            ++scattered_flag[current_volume];
            ++scattered_flag_VP[current_volume][selected_process];
            

            // Empty intersection time lists
            clear_intersection_table(&intersection_time_table);
            time_propagated_without_scattering = 0.0;
            #ifdef Union_trace_verbal_setting
              printf("SCATTERED SUCSSESFULLY \n");
              printf("r = (%f,%f,%f) v = (%f,%f,%f) \n",x,y,z,vx,vy,vz);
            
              if (enable_tagging && stop_tagging_ray == 0) printf("Before new process node: current_tagging_node->intensity = %f\n",current_tagging_node->intensity);
              if (enable_tagging && stop_tagging_ray == 0) printf("Before new process node: current_tagging_node->number_of_rays = %d\n",current_tagging_node->number_of_rays);
            #endif
            
            if (enable_tagging && stop_tagging_ray == 0)
                current_tagging_node = goto_process_node(current_tagging_node,selected_process,Volumes[current_volume],&stop_tagging_ray,stop_creating_nodes);
            
            #ifdef Union_trace_verbal_setting
              if (enable_tagging && stop_tagging_ray == 0) printf("After new process node: current_tagging_node->intensity = %f\n",current_tagging_node->intensity);
              if (enable_tagging && stop_tagging_ray == 0) printf("After new process node: current_tagging_node->number_of_rays = %d\n",current_tagging_node->number_of_rays);
            #endif
            
        } else {
            #ifdef Union_trace_verbal_setting
              printf("Propagate out of volume %d\n", current_volume);
              printf("r = (%f,%f,%f) v = (%f,%f,%f) \n",x,y,z,vx,vy,vz);
            #endif
            // Propagate neutron to found minimum time
            // PROP_DT(min_intersection_time - time_propagated_without_scattering);
            //time_to_boundery = min_intersection_time - time_propagated_without_scattering;
            x += time_to_boundery*vx;
            y += time_to_boundery*vy;
            z += time_to_boundery*vz;
            t += time_to_boundery;
            r[0] = x; r[1] = y; r[2] = z;
            ray_position = coords_set(x,y,z);
            ray_velocity = coords_set(vx,vy,vz);
            
            /*
            // Absorption moved to before testing if scattering occurs
            if (current_volume != 0) {
                if (Volumes[current_volume]->p_physics->is_vacuum == 0) {
                    // Absorption is done explicitly when propagating out of a volume, but between all scattering events is done implicitly
                   
                    // Old version
                    //length_to_boundery = (min_intersection_time - time_propagated_without_scattering) * v_length;
                    //p *= exp(-Volumes[current_volume]->p_physics->my_a*(2200/v_length)*length_to_boundery);
                    
                    if (Volumes[current_volume]->p_physics->number_of_processes == 0) {
                      // Optimized version
                      //p *= exp(-Volumes[current_volume]->p_physics->my_a*2200*time_to_boundery);
                    
                    //printf("name of material: %s \n",Volumes[current_volume]->name);
                    //printf("length to boundery  = %f\n",length_to_boundery);
                    //printf("absorption cross section  = %f\n",Volumes[current_volume]->p_physics->my_a);
                    //printf("chance to get through this length of absorber: %f \%\n",100*exp(-Volumes[current_volume]->p_physics->my_a*length_to_boundery));
                    }
                }
            }
            */
            
            time_propagated_without_scattering = min_intersection_time;
            SCATTER; // For debugging purposes
            #ifdef Union_trace_verbal_setting
              printf("r = (%f,%f,%f) v = (%f,%f,%f) \n",x,y,z,vx,vy,vz);
            #endif
            // Remove this entry from the intersection_time_table
            intersection_time_table.intersection_times[min_volume][min_solution] = -1;
            
            // Use destination list for corresponding intersection entry n,i) to find next volume
            #ifdef Union_trace_verbal_setting
              printf("PROPAGATION FROM VOLUME %d \n",current_volume);
            #endif
            if (min_volume == current_volume) {
                #ifdef Union_trace_verbal_setting
                  printf("min_volume == current_volume \n");
                #endif
                // List approach to finding the next volume.
                // When the ray intersects the current volume, the next volume must be on the destination list of the current volume
                // However, the reduced_destination_list can be investigated first, and depending on the results, the
                //  direct children of the volumes on the reduced destination list are investigated.
                // In the worst case, all direct children are investigated, which is eqvivalent to the entire destination list.
                // There is however a certain overhead in the logic needed to set up this tree, avoid duplicates of direct children, and so on.
                // This method is only faster than just checking the destination list when there are direct children (nested structures),
                //  but in general the tree method scales better with complexity, and is only slightly slower in simple cases.
                
                if (Volumes[current_volume]->geometry.destinations_list.num_elements == 1)
                    tree_next_volume = Volumes[current_volume]->geometry.destinations_list.elements[0];
                else {
                    ray_position = coords_set(x,y,z);
                    ray_velocity = coords_set(vx,vy,vz);
                    tree_next_volume = within_which_volume(ray_position,Volumes[current_volume]->geometry.reduced_destinations_list,Volumes[current_volume]->geometry.destinations_list,Volumes,&mask_status_list,number_of_volumes,pre_allocated1,pre_allocated2,pre_allocated3);
                }

                #ifdef Union_trace_verbal_setting
                  if (trace_verbal) printf("tree method moves from %d to %d\n",current_volume,tree_next_volume);
                
                  if (enable_tagging && stop_tagging_ray == 0) printf("Before new tree volume node: current_tagging_node->intensity = %f\n",current_tagging_node->intensity);
                  if (enable_tagging && stop_tagging_ray == 0) printf("Before new tree volume node: current_tagging_node->number_of_rays = %d\n",current_tagging_node->number_of_rays);
                #endif
                
                if (enable_tagging && stop_tagging_ray == 0)
                    current_tagging_node = goto_volume_node(current_tagging_node, current_volume, tree_next_volume, Volumes,&stop_tagging_ray,stop_creating_nodes);
                
                #ifdef Union_trace_verbal_setting
                  if (enable_tagging && stop_tagging_ray == 0) printf("After new tree volume node: current_tagging_node->intensity = %f\n",current_tagging_node->intensity);
                  if (enable_tagging && stop_tagging_ray == 0) printf("After new tree volume node: current_tagging_node->number_of_rays = %d\n",current_tagging_node->number_of_rays);
                #endif
                
                // Set next volume to the solution found in the tree method
                current_volume = tree_next_volume;
                update_current_mask_intersect_status(&current_mask_intersect_list_status, &mask_status_list, Volumes, &current_volume);
                #ifdef Union_trace_verbal_setting
                  print_1d_int_list(current_mask_intersect_list_status,"Updated current_mask_intersect_list_status");
                #endif
                
                
                // Debugging phase
                /*
                if (tree_next_volume == 0) {
                    volume_0_found=0;
                    for (start = check = Volumes[current_volume]->geometry.destinations_list.elements;check - start < Volumes[current_volume]->geometry.destinations_list.num_elements;check++) {
                        if (*check == 0) {
                            volume_0_found = 1;
                        }
                    }
                    if (volume_0_found==0) printf("ERROR The within_which_volume function returned 0 for a volume where volume 0 is not on the destination list!\n");
                }
                */
                
            } else {
                #ifdef Union_trace_verbal_setting
                  if (enable_tagging && stop_tagging_ray == 0) printf("Before new intersection volume node: current_tagging_node->intensity = %f\n",current_tagging_node->intensity);
                  if (enable_tagging && stop_tagging_ray == 0) printf("Before new intersection volume node: current_tagging_node->number_of_rays = %d\n",current_tagging_node->number_of_rays);
                #endif
            
                //if (enable_tagging) current_tagging_node = goto_volume_node(current_tagging_node, current_volume, min_volume, Volumes);
                
                
                
                // Mask update: If the min_volume is not a mask, things are simple, current_volume = min_volume.
                //               however, if it is a mask, the mask status will switch.
                //               if the mask status becomes one, the masked volumes inside may be the next volume (unless they are children of the mask)
                //               if the mask status becomes zero (and the current volume is masked by min_volume), the destinations list of the mask is searched
                //               if the mask status becomes zero (and the current volume is NOT masked by min volume), the current volume doesn't change
                
                
                if (Volumes[min_volume]->geometry.is_mask_volume == 0) {
                  #ifdef Union_trace_verbal_setting
                    printf("Min volume is not a mask, next volume = min volume\n");
                  #endif
                  if (enable_tagging && stop_tagging_ray == 0) current_tagging_node = goto_volume_node(current_tagging_node, current_volume, min_volume, Volumes,&stop_tagging_ray,stop_creating_nodes);
                  current_volume = min_volume;
                  //update_current_mask_intersect_status(&current_mask_intersect_list_status, &mask_status_list, Volumes, &current_volume);
                } else {
                  #ifdef Union_trace_verbal_setting
                    printf("Current volume is not a mask, complex decision tree\n");
                  #endif
                  if (mask_status_list.elements[Volumes[min_volume]->geometry.mask_index] == 1) {
                    // We are leaving the mask, change the status
                    #ifdef Union_trace_verbal_setting
                      printf("mask status changed from 1 to 0 as a mask is left\n");
                    #endif
                    mask_status_list.elements[Volumes[min_volume]->geometry.mask_index] = 0;
                    // If the current volume is masked by this mask, run within_which_volume using the masks destination list, otherwise keep the current volume
                    //if (on_int_list(Volumes[min_volume]->geometry.mask_list,current_volume))
                    if (on_int_list(Volumes[current_volume]->geometry.masked_by_list,min_volume) == 1) {
                      #ifdef Union_trace_verbal_setting
                        printf("The current volume was masked by this mask, and my need updating\n");
                      #endif
                      // In case of ANY mode, need to see if another mask on the masked_by list of the current volume is active, and if so, nothing happens
                      need_to_run_within_which_volume = 1;
                      if (Volumes[current_volume]->geometry.mask_mode == 2) {
                        for (mask_start=mask_check=Volumes[current_volume]->geometry.masked_by_mask_index_list.elements;mask_check-mask_start<Volumes[current_volume]->geometry.masked_by_mask_index_list.num_elements;mask_check++) {
                          if (mask_status_list.elements[*mask_check] == 1) {
                            // Nothing needs to be done, the effective mask status of the current volume is still 1
                            need_to_run_within_which_volume = 0;
                            break;
                          }
                        }
                      }
                      if (need_to_run_within_which_volume == 1) {
                        #ifdef Union_trace_verbal_setting
                          printf("The current volume was masked by this mask, and does need updating\n");
                        #endif
                        if (Volumes[min_volume]->geometry.destinations_list.num_elements == 1) {
                          #ifdef Union_trace_verbal_setting
                            printf("Only one element in the destination tree of the mask\n");
                          #endif
                          // If there is only one element on the destinations list (quite common) there is no reason to run within_which_volume
                          // Instead the mask status is calculated here
                          if (Volumes[Volumes[min_volume]->geometry.destinations_list.elements[0]]->geometry.is_masked_volume == 1) {
                            #ifdef Union_trace_verbal_setting
                              printf("The one element is however masked, so the mask status need to be calculated\n");
                            #endif
                            // figure out the effective mask status of this volume
                            if (Volumes[Volumes[min_volume]->geometry.destinations_list.elements[0]]->geometry.mask_mode == 2) { // ANY mask mode
                              tree_next_volume = 0;
                              for (mask_start=mask_check=Volumes[Volumes[min_volume]->geometry.destinations_list.elements[0]]->geometry.masked_by_mask_index_list.elements;mask_check-mask_start<Volumes[Volumes[min_volume]->geometry.destinations_list.elements[0]]->geometry.masked_by_mask_index_list.num_elements;mask_check++) {
                                if (mask_status_list.elements[*mask_check] == 1) {
                                  tree_next_volume = Volumes[min_volume]->geometry.destinations_list.elements[0];
                                  break;
                                }
                              }
                            } else { // ALL mask mode
                              tree_next_volume = Volumes[min_volume]->geometry.destinations_list.elements[0];
                              for (mask_start=mask_check=Volumes[Volumes[min_volume]->geometry.destinations_list.elements[0]]->geometry.masked_by_mask_index_list.elements;mask_check-mask_start<Volumes[Volumes[min_volume]->geometry.destinations_list.elements[0]]->geometry.masked_by_mask_index_list.num_elements;mask_check++) {
                                if (mask_status_list.elements[*mask_check] == 0) {
                                  tree_next_volume = 0;
                                  break;
                                }
                              }
                            }
                          } else tree_next_volume = Volumes[min_volume]->geometry.destinations_list.elements[0];
                          #ifdef Union_trace_verbal_setting
                            printf("The method found the next tree volume to be %d\n",tree_next_volume);
                          #endif
                          if (enable_tagging && stop_tagging_ray == 0) current_tagging_node = goto_volume_node(current_tagging_node, current_volume, tree_next_volume, Volumes,&stop_tagging_ray,stop_creating_nodes);
                          current_volume = tree_next_volume;
                          //update_current_mask_intersect_status(&current_mask_intersect_list_status, &mask_status_list, Volumes, &current_volume);
                        } else {
                          #ifdef Union_trace_verbal_setting
                            printf("Many elements in destinations list, use within_which_volume\n");
                          #endif
                          ray_position = coords_set(x,y,z);
                          ray_velocity = coords_set(vx,vy,vz);
                          tree_next_volume = within_which_volume(ray_position,Volumes[min_volume]->geometry.reduced_destinations_list,Volumes[min_volume]->geometry.destinations_list,Volumes,&mask_status_list,number_of_volumes,pre_allocated1,pre_allocated2,pre_allocated3);
                        // } Bug fixed on 27/11/2016
                          if (enable_tagging && stop_tagging_ray == 0) current_tagging_node = goto_volume_node(current_tagging_node, current_volume, tree_next_volume, Volumes,&stop_tagging_ray,stop_creating_nodes);
                          current_volume = tree_next_volume;
                          #ifdef Union_trace_verbal_setting
                            printf("Set new new volume to %d\n",tree_next_volume);
                          #endif
                        } // Moved here on 27/11/2016, problem was two calls to current_tagging_node when only one element in destinations list
                        //update_current_mask_intersect_status(&current_mask_intersect_list_status, &mask_status_list, Volumes, &current_volume);
                      } else {
                        #ifdef Union_trace_verbal_setting
                          printf("Did not need updating, as another mask was covering the volume\n");
                        #endif
                      }
                    }

                  } else {
                    // Here beccause the mask status of the mask that is intersected was 0, and it is thus switched to 1
                    mask_status_list.elements[Volumes[min_volume]->geometry.mask_index] = 1;
                    // When entering a mask, the new highest priority volume may be one of the masked volumes, if not we keep the current volume
                    ray_position = coords_set(x,y,z);
                    ray_velocity = coords_set(vx,vy,vz);
                    //tree_next_volume = within_which_volume(ray_position,Volumes[min_volume]->geometry.mask_list,Volumes[min_volume]->geometry.destinations_list,Volumes,&mask_status_list,number_of_volumes,pre_allocated1,pre_allocated2,pre_allocated3);
                    // Bug found on the 2/9 2016, the destinations_list of a mask does not contain the volumes inside it. Could make an additional list for this.
                    // The temporary fix will be to use the mask list for both reduced destinations list and destinations list.
                    tree_next_volume = within_which_volume(ray_position,Volumes[min_volume]->geometry.mask_list,Volumes[min_volume]->geometry.mask_list,Volumes,&mask_status_list,number_of_volumes,pre_allocated1,pre_allocated2,pre_allocated3);
                    // if within_which_volume returns 0, no result was found (volume 0 can not be masked, so it could not be on the mask list)
                    if (tree_next_volume != 0) {
                      if (Volumes[tree_next_volume]->geometry.priority_value > Volumes[current_volume]->geometry.priority_value) {
                        // In case the current volume has a higher priority, nothing happens, otherwise change current volume
                        if (enable_tagging && stop_tagging_ray == 0) current_tagging_node = goto_volume_node(current_tagging_node, current_volume, tree_next_volume, Volumes,&stop_tagging_ray,stop_creating_nodes);
                        current_volume = tree_next_volume;
                      }
                    }
                    //update_current_mask_intersect_status(&current_mask_intersect_list_status, &mask_status_list, Volumes, &current_volume);
                  }
                }
                
                // Regardless of the outcome of the above code, either the mask status or current volume have changed, and thus a effective mask update is needed.
                update_current_mask_intersect_status(&current_mask_intersect_list_status, &mask_status_list, Volumes, &current_volume);
                #ifdef Union_trace_verbal_setting
                  print_1d_int_list(mask_status_list,"Updated mask status list");
                  print_1d_int_list(current_mask_intersect_list_status,"Updated current_mask_intersect_list_status");
                
                  if (enable_tagging) printf("After new intersection volume node: current_tagging_node->intensity = %f\n",current_tagging_node->intensity);
                  if (enable_tagging) printf("After new intersection volume node: current_tagging_node->number_of_rays = %d\n",current_tagging_node->number_of_rays);
                #endif
                
            }
            if (Volumes[current_volume]->geometry.is_exit_volume==1) {
                    done = 1; // Exit volumes allow the ray to escape the component
                    ray_sucseeded = 1; // Allows the ray to
                    /*
                    // Moved to after while loop to collect code
                    if (enable_tagging && stop_tagging_ray == 0)
                        add_statistics_to_node(current_tagging_node,&ray_position, &ray_velocity, &p, &tagging_leaf_counter);
                
                    x += vx*1E-9; y += vy*1E-9; z += vz*1E-9; t += 1E-9;
                    */
            }
            #ifdef Union_trace_verbal_setting
              printf(" TO VOLUME %d \n",current_volume);
            #endif
        }
    } else { // Here because a shortest time is not found
        if (current_volume == 0) {
            done = 1;
            ray_sucseeded = 1;
            
        } else { // Check for errors (debugging phase)
            if (error_msg == 0) {
              component_error_msg++;
              ray_sucseeded = 0;
              done = 1; // stop the loop
              printf("\n----------------------------------------------------------------------------------------------------\n");
              printf("Union_master %s: Somehow reached a situation with no intersection time found, but still inside volume %d instead of 0\n",NAME_CURRENT_COMP,current_volume);
              for (volume_index = 1; volume_index < number_of_volumes; volume_index++) {
                if (Volumes[volume_index]->geometry.within_function(ray_position,&Volumes[volume_index]->geometry) == 1)
                    printf("The ray is in volume       %d\n",volume_index);
              }
              
              print_1d_int_list(mask_status_list,"mask status list");
              for (iterate=0;iterate<number_of_volumes;iterate++)
                 printf("%d:%d - ",iterate,scattered_flag[iterate]);
              printf("\n");
              printf("r = (%f,%f,%f) v = (%f,%f,%f) \n",x,y,z,vx,vy,vz);
              
              printf("Trace error number (%d/100) \n",component_error_msg);
              //print_intersection_table(intersection_time_table);
              //printf("Cluster is loosing the thread, debug for this behavior\n");
              //printf("global_all_volume_logger_list.num_elements = %d\n",global_all_volume_logger_list.num_elements);
              //printf("global_specific_volumes_logger_list.num_elements = %d\n",global_specific_volumes_logger_list.num_elements);
              
              //exit(1); // Temporary
            }
            error_msg++;
            
            //exit(1); // temp for debug
            
            if (component_error_msg > 100) {
                printf("To many errors encountered, exiting. \n");
                exit(1);
            }
        }
    }
    /*
    */
    if (limit == 0) {done = 1; ray_sucseeded = 0; printf("Reached limit on number of interactions, and discarded the neutron, was in volume %d\n",current_volume);ABSORB;}
    #ifdef Union_trace_verbal_setting
      printf("----------- END OF WHILE LOOP --------------------------------------\n");
    #endif
    //printf("This ray did %d iterations in the while loop\n",1000-limit);
    
  }
  // Could move all add_statistics and similar to this point, but need to filter for failed rays
  if (ray_sucseeded == 1) {
    
    /*
    // Instead of keeping global and specific loggers apart, lets do them in one loop using the loggers_with_data_array
    // Only needed for conditionals
    // Loop over global loggers, may or may not have a conditional, only necessary to do anything here when they do.
    for (log_index=0;log_index<global_all_volume_logger_list.num_elements;log_index++) {
      if (global_all_volume_logger_list.elements[log_index].logger->has_conditional == 1) {
        if (1 == conditional_return = global_all_volume_logger_list.elements[log_index].logger->conditional(scattered_flag,scattered_flag_VP,global_loggers.elements[log_index].conditional_data_union,k_final,current_volume,x,y,z)) {
          global_all_volume_logger_list.elements[log_index].logger->function_pointers->temp_to_perm(global_loggers.elements[log_index].logger_data_union);
          
        }
        if (global_all_volume_logger_list.elements[log_index].logger->conditional_extend_index != -1)
          // The user can set a condtional_extend_index, so that the evaluation of this specific conditional can be taken easily from extend
          conditional_extend_array.elements[global_all_volume_logger_list.elements[log_index].logger->conditional_extend_index] = conditional_return;
          // Do not need to reset these to 0, as they will be manually set to 0 if not fulfilled
      }
    }
    */
    #ifdef Union_trace_verbal_setting
      printf("----------- logger loop --------------------------------------\n");
    #endif
    // Loggers attatched to specific volumes need to be handled with care to avoid looping over all loggers for every ray
    //for (log_index=0;log_index<loggers_with_data_array.used_elements;log_index++) {
    // TEST
    if (enable_conditionals == 1) {
      for (log_index=loggers_with_data_array.used_elements-1;log_index>-1;log_index--) {
        // Check all conditionals attatched to the current logger
        this_logger = loggers_with_data_array.logger_pointers[log_index];
        conditional_status = 1;
        for (iterate=0;iterate<loggers_with_data_array.logger_pointers[log_index]->conditional_list.num_elements;iterate++) {
          // Call this particular conditional. If it fails, report the status and break
          //printf("checking conditional! \n");
          #ifdef Union_trace_verbal_setting
            printf("Checking conditional number %d for logger named %s \n",iterate,loggers_with_data_array.logger_pointers[log_index]->name);
          #endif
          if (0 == this_logger->conditional_list.conditional_functions[iterate](
                         this_logger->conditional_list.p_data_unions[iterate],
                         &ray_position,&ray_velocity,&p,&t,&current_volume,
                         &number_of_scattering_events,scattered_flag,scattered_flag_VP)) {
            conditional_status = 0;
            break;
          }
        }
        if (conditional_status == 1) {
          // If a logger does not have a conditional, it will write directly to perm, and not even add it to the loggers_with_data_array, thus we know the temp_to_perm function needs to be called
          // The input for the temp_to_perm function is a pointer to the logger_data_union for the appropriate logger
          
          if (loggers_with_data_array.logger_pointers[log_index]->function_pointers.select_t_to_p == 1) {
            loggers_with_data_array.logger_pointers[log_index]->function_pointers.temp_to_perm(&loggers_with_data_array.logger_pointers[log_index]->data_union);
          }
          else if (loggers_with_data_array.logger_pointers[log_index]->function_pointers.select_t_to_p == 2) {
            loggers_with_data_array.logger_pointers[log_index]->function_pointers.temp_to_perm_final_p(&loggers_with_data_array.logger_pointers[log_index]->data_union,p);
          }
        
          // Luxury feature to be added later
          if (loggers_with_data_array.logger_pointers[log_index]->logger_extend_index != -1) {
            #ifdef Union_trace_verbal_setting
              printf("Updating logger_conditional_extend_array[%d] to 1 (max length = %d)\n",loggers_with_data_array.logger_pointers[log_index]->logger_extend_index,max_conditional_extend_index);
            #endif
            // The user can set a condtional_extend_index, so that the evaluation of this specific conditional can be taken easily from extend
            logger_conditional_extend_array[loggers_with_data_array.logger_pointers[log_index]->logger_extend_index] = 1;
            // Are all reset to 0 for each new ray
            #ifdef Union_trace_verbal_setting
              printf("Updated extend index sucessfully\n");
            #endif
          
            //printf("extend_array[%d] = 1 \n",loggers_with_data_array.logger_pointers[log_index]->logger_extend_index);
          }
        
          // Need to remove the current element from logger_with_data as it has been cleared and written to disk
          // The remaining elements is passed on to the next Union_master as it may fulfill the conditional after that master
          if (global_master_list.elements[global_master_list.num_elements-1].component_index != INDEX_CURRENT_COMP) {
            // Move current logger pointer in logger_with_data to end position
            loggers_with_data_array.logger_pointers[log_index] = loggers_with_data_array.logger_pointers[loggers_with_data_array.used_elements-1];
            // Decrease logger_with_data.used_elements with 1
            loggers_with_data_array.used_elements--;
          
          }
        
        
        } else {
           // Conditional status was 0, clear temp data for this logger, but only if this is the last Union_master,
           //  as the logger data can be written if one of the ray fulfills the conditional afer one of the
           //  subsequent Union masters.
           // The job of cleaning was moved to the start of trace on 15/5/2017
           //if (global_master_list.elements[global_master_list.num_elements-1].component_index == INDEX_CURRENT_COMP)
           //  loggers_with_data_array.logger_pointers[log_index]->function_pointers.clear_temp(&loggers_with_data_array.logger_pointers[log_index]->data_union);
        }
      }
    }
    
    if (enable_tagging && stop_tagging_ray == 0) {
      conditional_status = 1;
      for (iterate=0;iterate<tagging_conditional_list->num_elements;iterate++) {
        // Call this particular conditional. If it fails, report the status and break
        // Since a conditional can work for a logger and master_tagging at the same time, it may be evaluated twice
        //printf("checking conditional! \n");
        #ifdef Union_trace_verbal_setting
          printf("Checking tagging conditional number %d\n",iterate);
        #endif
        if (0 == tagging_conditional_list->conditional_functions[iterate](
                         tagging_conditional_list->p_data_unions[iterate],
                         &ray_position,&ray_velocity,&p,&t,&current_volume,
                         &number_of_scattering_events,scattered_flag,scattered_flag_VP)) {
          conditional_status = 0;
          break;
        }
      }
      if (conditional_status == 1) {
        tagging_conditional_extend = 1;
        #ifdef Union_trace_verbal_setting
          printf("Before adding statistics to node: current_tagging_nodbe->intensity = %f\n",current_tagging_node->intensity);
          printf("Before adding statistics to node: current_tagging_node->number_of_rays = %d\n",current_tagging_node->number_of_rays);
        #endif
        
          add_statistics_to_node(current_tagging_node,&ray_position, &ray_velocity, &p, &tagging_leaf_counter);
        
        #ifdef Union_trace_verbal_setting
          printf("After adding statistics to node: current_tagging_node->intensity = %f\n",current_tagging_node->intensity);
          printf("After adding statistics to node: current_tagging_node->number_of_rays = %d\n",current_tagging_node->number_of_rays);
        #endif
      }
    }
    
    // Move the rays a nano meter away from the surface it left, in case activation counter > 1, this will prevent the ray from starting on a volume boundery
    x += vx*1E-9; y += vy*1E-9; z += vz*1E-9; t += 1E-9;
        
  } else {
    ABSORB; // Absorb rays that didn't exit correctly for whatever reason
    // Could error log here
  }
  
  // TEST
  // Stores nubmer of scattering events in global master list so that another master with inherit_number_of_scattering_events can continue
  global_master_list.elements[this_global_master_index].stored_number_of_scattering_events = number_of_scattering_events;
  
  
}
#line 13752 "./Union_test_mask.c"
/* 'test_sample=Union_master()' component instance extend code */
    SIG_MESSAGE("test_sample (Trace:Extend)");
#line 98 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
	if (scattered_flag[sample_1_index] > 0) scattered_1 = 1; else scattered_1 = 0;
	if (scattered_flag[sample_2_index] > 0) scattered_2 = 1; else scattered_2 = 0;
#line 13758 "./Union_test_mask.c"
}   /* End of test_sample=Union_master() SETTING parameter declarations. */
#undef focus_data_index
#undef temporary_focus_data
#undef number_of_processes_array
#undef safty_distance2
#undef safty_distance
#undef max_conditional_extend_index
#undef tagging_conditional_extend
#undef logger_conditional_extend_array
#undef free_tagging_conditioanl_list
#undef tagging_conditional_list
#undef conditional_status
#undef this_logger
#undef p_old
#undef Volume_copies_allocated
#undef geometry_component_index_list
#undef mask_volume_index_list
#undef current_mask_intersect_list_status
#undef mask_status_list
#undef mask_iterate
#undef mask_index_main
#undef need_to_run_within_which_volume
#undef number_of_masked_volumes
#undef number_of_masks
#undef mc_transmission_probability
#undef real_transmission_probability
#undef number_of_scattering_events
#undef tagging_leaf_counter
#undef current_tagging_node
#undef master_tagging_node_list
#undef enable_tagging_check
#undef stop_creating_nodes
#undef stop_tagging_ray
#undef enable_tagging
#undef rotated_position
#undef non_rotated_position
#undef temp_rotation_matrix
#undef master_transposed_rotation_matrix
#undef scattered_flag_VP
#undef scattered_flag
#undef volume_0_found
#undef ray_velocity_final
#undef ray_velocity
#undef ray_position
#undef pre_allocated3
#undef pre_allocated2
#undef pre_allocated1
#undef tree_next_volume
#undef geometry_output
#undef intersection_with_children
#undef start
#undef check
#undef number_of_solutions_static
#undef number_of_solutions
#undef current_volume
#undef done
#undef next_volume_priority
#undef next_volume
#undef a_next_volume_found
#undef time_propagated_without_scattering
#undef scattering_event
#undef selected_process
#undef time_to_boundery
#undef length_to_boundery
#undef length_to_scattering
#undef time_to_scattering
#undef mc_prop
#undef culmative_probability
#undef my_sum_plus_abs
#undef my_sum
#undef v_length
#undef k_old
#undef k_new
#undef k
#undef my_trace_fraction_control
#undef p_my_trace
#undef my_trace
#undef process_start
#undef process
#undef min_intersection_time
#undef intersection_time
#undef time_found
#undef min_volume
#undef min_solution
#undef solution
#undef limit
#undef max_number_of_processes
#undef solutions
#undef process_index
#undef volume_index
#undef number_of_volumes
#undef string_output
#undef component_error_msg
#undef error_msg
#undef v
#undef r_start
#undef r
#undef starting_lists
#undef Geometries
#undef Volumes
#undef intersection_time_table
#undef geometry_list_index
#undef previous_master_index
#undef this_global_master_index
#undef global_master_element
#undef starting_volume_warning
#undef finally_verbal
#undef trace_verbal
#undef list_verbal
#undef verbal
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbComptest_sample:
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

  /* TRACE Component m4pi [13] */
  mccoordschange(mcposrm4pi, mcrotrm4pi,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component m4pi (without coords transformations) */
  mcJumpTrace_m4pi:
  SIG_MESSAGE("m4pi (Trace)");
  mcDEBUG_COMP("m4pi")
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

#define mcabsorbComp mcabsorbCompm4pi
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
#define mccompcurname  m4pi
#define mccompcurtype  PSD_monitor_4PI
#define mccompcurindex 13
#define nx mccm4pi_nx
#define ny mccm4pi_ny
#define PSD_N mccm4pi_PSD_N
#define PSD_p mccm4pi_PSD_p
#define PSD_p2 mccm4pi_PSD_p2
{   /* Declarations of m4pi=PSD_monitor_4PI() SETTING parameters. */
char* filename = mccm4pi_filename;
MCNUM radius = mccm4pi_radius;
MCNUM restore_neutron = mccm4pi_restore_neutron;
int nowritefile = mccm4pi_nowritefile;
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
#line 14014 "./Union_test_mask.c"
}   /* End of m4pi=PSD_monitor_4PI() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef ny
#undef nx
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompm4pi:
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

  /* TRACE Component m4pi_1 [14] */
  mccoordschange(mcposrm4pi_1, mcrotrm4pi_1,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component m4pi_1 (without coords transformations) */
  mcJumpTrace_m4pi_1:
  SIG_MESSAGE("m4pi_1 (Trace)");
  mcDEBUG_COMP("m4pi_1")
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

#define mcabsorbComp mcabsorbCompm4pi_1
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
#define mccompcurname  m4pi_1
#define mccompcurtype  PSD_monitor_4PI
#define mccompcurindex 14
#define nx mccm4pi_1_nx
#define ny mccm4pi_1_ny
#define PSD_N mccm4pi_1_PSD_N
#define PSD_p mccm4pi_1_PSD_p
#define PSD_p2 mccm4pi_1_PSD_p2
{   /* Declarations of m4pi_1=PSD_monitor_4PI() SETTING parameters. */
char* filename = mccm4pi_1_filename;
MCNUM radius = mccm4pi_1_radius;
MCNUM restore_neutron = mccm4pi_1_restore_neutron;
int nowritefile = mccm4pi_1_nowritefile;
/* 'm4pi_1=PSD_monitor_4PI()' component instance has conditional execution */
if (( scattered_1 ))

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
#line 14168 "./Union_test_mask.c"
}   /* End of m4pi_1=PSD_monitor_4PI() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef ny
#undef nx
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompm4pi_1:
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

  /* TRACE Component m4pi_2 [15] */
  mccoordschange(mcposrm4pi_2, mcrotrm4pi_2,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component m4pi_2 (without coords transformations) */
  mcJumpTrace_m4pi_2:
  SIG_MESSAGE("m4pi_2 (Trace)");
  mcDEBUG_COMP("m4pi_2")
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

#define mcabsorbComp mcabsorbCompm4pi_2
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
#define mccompcurname  m4pi_2
#define mccompcurtype  PSD_monitor_4PI
#define mccompcurindex 15
#define nx mccm4pi_2_nx
#define ny mccm4pi_2_ny
#define PSD_N mccm4pi_2_PSD_N
#define PSD_p mccm4pi_2_PSD_p
#define PSD_p2 mccm4pi_2_PSD_p2
{   /* Declarations of m4pi_2=PSD_monitor_4PI() SETTING parameters. */
char* filename = mccm4pi_2_filename;
MCNUM radius = mccm4pi_2_radius;
MCNUM restore_neutron = mccm4pi_2_restore_neutron;
int nowritefile = mccm4pi_2_nowritefile;
/* 'm4pi_2=PSD_monitor_4PI()' component instance has conditional execution */
if (( scattered_2 ))

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
#line 14322 "./Union_test_mask.c"
}   /* End of m4pi_2=PSD_monitor_4PI() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef ny
#undef nx
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompm4pi_2:
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

  /* TRACE Component arm_1 [16] */
  mccoordschange(mcposrarm_1, mcrotrarm_1,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component arm_1 (without coords transformations) */
  mcJumpTrace_arm_1:
  SIG_MESSAGE("arm_1 (Trace)");
  mcDEBUG_COMP("arm_1")
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

#define mcabsorbComp mcabsorbComparm_1
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
#define mccompcurname  arm_1
#define mccompcurtype  Arm
#define mccompcurindex 16
/* 'arm_1=Arm()' component instance extend code */
    SIG_MESSAGE("arm_1 (Trace:Extend)");
#line 116 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/Union_test_mask/Union_test_mask.instr"
	if (scattered_1 + scattered_2 != 2) ABSORB;
#line 14436 "./Union_test_mask.c"
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbComparm_1:
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

  /* TRACE Component m4pi_both [17] */
  mccoordschange(mcposrm4pi_both, mcrotrm4pi_both,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component m4pi_both (without coords transformations) */
  mcJumpTrace_m4pi_both:
  SIG_MESSAGE("m4pi_both (Trace)");
  mcDEBUG_COMP("m4pi_both")
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

#define mcabsorbComp mcabsorbCompm4pi_both
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
#define mccompcurname  m4pi_both
#define mccompcurtype  PSD_monitor_4PI
#define mccompcurindex 17
#define nx mccm4pi_both_nx
#define ny mccm4pi_both_ny
#define PSD_N mccm4pi_both_PSD_N
#define PSD_p mccm4pi_both_PSD_p
#define PSD_p2 mccm4pi_both_PSD_p2
{   /* Declarations of m4pi_both=PSD_monitor_4PI() SETTING parameters. */
char* filename = mccm4pi_both_filename;
MCNUM radius = mccm4pi_both_radius;
MCNUM restore_neutron = mccm4pi_both_restore_neutron;
int nowritefile = mccm4pi_both_nowritefile;
/* 'm4pi_both=PSD_monitor_4PI()' component instance has conditional execution */
if (( scattered_1 && scattered_2 ))

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
#line 14584 "./Union_test_mask.c"
}   /* End of m4pi_both=PSD_monitor_4PI() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef ny
#undef nx
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompm4pi_both:
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

  /* User SAVE code for component 'a1'. */
  SIG_MESSAGE("a1 (Save)");
#define mccompcurname  a1
#define mccompcurtype  Progress_bar
#define mccompcurindex 3
#define IntermediateCnts mcca1_IntermediateCnts
#define StartTime mcca1_StartTime
#define EndTime mcca1_EndTime
#define CurrentTime mcca1_CurrentTime
{   /* Declarations of a1=Progress_bar() SETTING parameters. */
char* profile = mcca1_profile;
MCNUM percent = mcca1_percent;
MCNUM flag_save = mcca1_flag_save;
MCNUM minutes = mcca1_minutes;
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
#line 14698 "./Union_test_mask.c"
}   /* End of a1=Progress_bar() SETTING parameter declarations. */
#undef CurrentTime
#undef EndTime
#undef StartTime
#undef IntermediateCnts
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'm4pi'. */
  SIG_MESSAGE("m4pi (Save)");
#define mccompcurname  m4pi
#define mccompcurtype  PSD_monitor_4PI
#define mccompcurindex 13
#define nx mccm4pi_nx
#define ny mccm4pi_ny
#define PSD_N mccm4pi_PSD_N
#define PSD_p mccm4pi_PSD_p
#define PSD_p2 mccm4pi_PSD_p2
{   /* Declarations of m4pi=PSD_monitor_4PI() SETTING parameters. */
char* filename = mccm4pi_filename;
MCNUM radius = mccm4pi_radius;
MCNUM restore_neutron = mccm4pi_restore_neutron;
int nowritefile = mccm4pi_nowritefile;
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
#line 14736 "./Union_test_mask.c"
}   /* End of m4pi=PSD_monitor_4PI() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef ny
#undef nx
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'm4pi_1'. */
  SIG_MESSAGE("m4pi_1 (Save)");
#define mccompcurname  m4pi_1
#define mccompcurtype  PSD_monitor_4PI
#define mccompcurindex 14
#define nx mccm4pi_1_nx
#define ny mccm4pi_1_ny
#define PSD_N mccm4pi_1_PSD_N
#define PSD_p mccm4pi_1_PSD_p
#define PSD_p2 mccm4pi_1_PSD_p2
{   /* Declarations of m4pi_1=PSD_monitor_4PI() SETTING parameters. */
char* filename = mccm4pi_1_filename;
MCNUM radius = mccm4pi_1_radius;
MCNUM restore_neutron = mccm4pi_1_restore_neutron;
int nowritefile = mccm4pi_1_nowritefile;
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
#line 14775 "./Union_test_mask.c"
}   /* End of m4pi_1=PSD_monitor_4PI() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef ny
#undef nx
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'm4pi_2'. */
  SIG_MESSAGE("m4pi_2 (Save)");
#define mccompcurname  m4pi_2
#define mccompcurtype  PSD_monitor_4PI
#define mccompcurindex 15
#define nx mccm4pi_2_nx
#define ny mccm4pi_2_ny
#define PSD_N mccm4pi_2_PSD_N
#define PSD_p mccm4pi_2_PSD_p
#define PSD_p2 mccm4pi_2_PSD_p2
{   /* Declarations of m4pi_2=PSD_monitor_4PI() SETTING parameters. */
char* filename = mccm4pi_2_filename;
MCNUM radius = mccm4pi_2_radius;
MCNUM restore_neutron = mccm4pi_2_restore_neutron;
int nowritefile = mccm4pi_2_nowritefile;
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
#line 14814 "./Union_test_mask.c"
}   /* End of m4pi_2=PSD_monitor_4PI() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef ny
#undef nx
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'm4pi_both'. */
  SIG_MESSAGE("m4pi_both (Save)");
#define mccompcurname  m4pi_both
#define mccompcurtype  PSD_monitor_4PI
#define mccompcurindex 17
#define nx mccm4pi_both_nx
#define ny mccm4pi_both_ny
#define PSD_N mccm4pi_both_PSD_N
#define PSD_p mccm4pi_both_PSD_p
#define PSD_p2 mccm4pi_both_PSD_p2
{   /* Declarations of m4pi_both=PSD_monitor_4PI() SETTING parameters. */
char* filename = mccm4pi_both_filename;
MCNUM radius = mccm4pi_both_radius;
MCNUM restore_neutron = mccm4pi_both_restore_neutron;
int nowritefile = mccm4pi_both_nowritefile;
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
#line 14853 "./Union_test_mask.c"
}   /* End of m4pi_both=PSD_monitor_4PI() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef ny
#undef nx
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  if (!handle) mcsiminfo_close(); 
} /* end save */
void mcfinally(void) {
  /* User component FINALLY code. */
  mcsiminfo_init(NULL);
  mcsave(mcsiminfo_file); /* save data when simulation ends */

  /* User FINALLY code for component 'test_incoherent'. */
  SIG_MESSAGE("test_incoherent (Finally)");
#define mccompcurname  test_incoherent
#define mccompcurtype  Incoherent_process
#define mccompcurindex 1
#define This_process mcctest_incoherent_This_process
#define Incoherent_storage mcctest_incoherent_Incoherent_storage
#define effective_my_scattering mcctest_incoherent_effective_my_scattering
{   /* Declarations of test_incoherent=Incoherent_process() SETTING parameters. */
MCNUM sigma = mcctest_incoherent_sigma;
MCNUM f_QE = mcctest_incoherent_f_QE;
MCNUM gamma = mcctest_incoherent_gamma;
MCNUM packing_factor = mcctest_incoherent_packing_factor;
MCNUM unit_cell_volume = mcctest_incoherent_unit_cell_volume;
MCNUM interact_fraction = mcctest_incoherent_interact_fraction;
#line 171 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../contrib/union/Incoherent_process.comp"
{
// Since the process and it's storage is a static allocation, there is nothing to deallocate

}
#line 14891 "./Union_test_mask.c"
}   /* End of test_incoherent=Incoherent_process() SETTING parameter declarations. */
#undef effective_my_scattering
#undef Incoherent_storage
#undef This_process
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[1]) fprintf(stderr, "Warning: No neutron could reach Component[1] test_incoherent\n");
    if (mcAbsorbProp[1]) fprintf(stderr, "Warning: %g events were removed in Component[1] test_incoherent=Incoherent_process()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[1]);
  /* User FINALLY code for component 'test_material'. */
  SIG_MESSAGE("test_material (Finally)");
#define mccompcurname  test_material
#define mccompcurtype  Union_make_material
#define mccompcurindex 2
#define loop_index mcctest_material_loop_index
#define this_material mcctest_material_this_material
#define accepted_processes mcctest_material_accepted_processes
#define global_material_element mcctest_material_global_material_element
{   /* Declarations of test_material=Union_make_material() SETTING parameters. */
char* process_string = mcctest_material_process_string;
MCNUM my_absorption = mcctest_material_my_absorption;
MCNUM absorber = mcctest_material_absorber;
#line 259 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../contrib/union/Union_make_material.comp"
{
// The elements of the scattering array used static allocation and is thus deallocated automatically
if (this_material.number_of_processes > 0) free(this_material.p_scattering_array);
if (accepted_processes.num_elements > 0) free(accepted_processes.elements);

// Checking if any Union volumes are defined after the master component
#ifdef MASTER_DETECTOR
  #ifdef ANY_GEOMETRY_DETECTOR_DECLARE
    #ifndef MASTER_DETECTOR_WARNING
      for (loop_index=0;loop_index<global_geometry_list.num_elements;loop_index++) {
        if (global_geometry_list.elements[loop_index].component_index > global_master_list.elements[global_master_list.num_elements-1].component_index) {
          printf("WARNING: No Union_master component defined after Union volume named %s, this components did not affect the simulation in any way.\n",global_geometry_list.elements[loop_index].name);
        }
      }
      // Decided to have this as a warning without exiting the simulation
      // In order to only show this warning once, the MASTER_DETECTOR_WARNING is defined
      #define MASTER_DETECTOR_WARNING dummy
    #endif
  #endif
#endif

// Checking if the user remembered to put in a Union_master
#ifndef MASTER_DETECTOR
  #ifdef ANY_GEOMETRY_DETECTOR_DECLARE
    #ifndef MASTER_DETECTOR_WARNING
      printf("\nWARNING: No Union_master component used, these components did not affect the simulation in any way:\n");
      for (loop_index=0;loop_index<global_geometry_list.num_elements;loop_index++)
        printf("  %s\n",global_geometry_list.elements[loop_index].name);
      printf("\n");
      // Decided to have this as a warning without exiting the simulation
      // In order to only show this warning once, the MASTER_DETECTOR_WARNING is defined
      #define MASTER_DETECTOR_WARNING dummy
    #endif
  #endif
#endif


}
#line 14954 "./Union_test_mask.c"
}   /* End of test_material=Union_make_material() SETTING parameter declarations. */
#undef global_material_element
#undef accepted_processes
#undef this_material
#undef loop_index
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[2]) fprintf(stderr, "Warning: No neutron could reach Component[2] test_material\n");
    if (mcAbsorbProp[2]) fprintf(stderr, "Warning: %g events were removed in Component[2] test_material=Union_make_material()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[2]);
  /* User FINALLY code for component 'a1'. */
  SIG_MESSAGE("a1 (Finally)");
#define mccompcurname  a1
#define mccompcurtype  Progress_bar
#define mccompcurindex 3
#define IntermediateCnts mcca1_IntermediateCnts
#define StartTime mcca1_StartTime
#define EndTime mcca1_EndTime
#define CurrentTime mcca1_CurrentTime
{   /* Declarations of a1=Progress_bar() SETTING parameters. */
char* profile = mcca1_profile;
MCNUM percent = mcca1_percent;
MCNUM flag_save = mcca1_flag_save;
MCNUM minutes = mcca1_minutes;
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
#line 14993 "./Union_test_mask.c"
}   /* End of a1=Progress_bar() SETTING parameter declarations. */
#undef CurrentTime
#undef EndTime
#undef StartTime
#undef IntermediateCnts
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[3]) fprintf(stderr, "Warning: No neutron could reach Component[3] a1\n");
    if (mcAbsorbProp[3]) fprintf(stderr, "Warning: %g events were removed in Component[3] a1=Progress_bar()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[3]);
    if (!mcNCounter[4]) fprintf(stderr, "Warning: No neutron could reach Component[4] source\n");
    if (mcAbsorbProp[4]) fprintf(stderr, "Warning: %g events were removed in Component[4] source=Source_div()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[4]);
    if (!mcNCounter[5]) fprintf(stderr, "Warning: No neutron could reach Component[5] beam_center\n");
    if (mcAbsorbProp[5]) fprintf(stderr, "Warning: %g events were removed in Component[5] beam_center=Arm()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[5]);
    if (!mcNCounter[6]) fprintf(stderr, "Warning: No neutron could reach Component[6] volume_1\n");
    if (mcAbsorbProp[6]) fprintf(stderr, "Warning: %g events were removed in Component[6] volume_1=Union_box()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[6]);
    if (!mcNCounter[7]) fprintf(stderr, "Warning: No neutron could reach Component[7] volume_2\n");
    if (mcAbsorbProp[7]) fprintf(stderr, "Warning: %g events were removed in Component[7] volume_2=Union_box()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[7]);
    if (!mcNCounter[8]) fprintf(stderr, "Warning: No neutron could reach Component[8] mask_1_1\n");
    if (mcAbsorbProp[8]) fprintf(stderr, "Warning: %g events were removed in Component[8] mask_1_1=Union_box()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[8]);
    if (!mcNCounter[9]) fprintf(stderr, "Warning: No neutron could reach Component[9] mask_1_2\n");
    if (mcAbsorbProp[9]) fprintf(stderr, "Warning: %g events were removed in Component[9] mask_1_2=Union_box()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[9]);
    if (!mcNCounter[10]) fprintf(stderr, "Warning: No neutron could reach Component[10] mask_2_1\n");
    if (mcAbsorbProp[10]) fprintf(stderr, "Warning: %g events were removed in Component[10] mask_2_1=Union_box()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[10]);
    if (!mcNCounter[11]) fprintf(stderr, "Warning: No neutron could reach Component[11] mask_2_2\n");
    if (mcAbsorbProp[11]) fprintf(stderr, "Warning: %g events were removed in Component[11] mask_2_2=Union_box()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[11]);
  /* User FINALLY code for component 'test_sample'. */
  SIG_MESSAGE("test_sample (Finally)");
#define mccompcurname  test_sample
#define mccompcurtype  Union_master
#define mccompcurindex 12
#define verbal mcctest_sample_verbal
#define list_verbal mcctest_sample_list_verbal
#define trace_verbal mcctest_sample_trace_verbal
#define finally_verbal mcctest_sample_finally_verbal
#define starting_volume_warning mcctest_sample_starting_volume_warning
#define global_master_element mcctest_sample_global_master_element
#define this_global_master_index mcctest_sample_this_global_master_index
#define previous_master_index mcctest_sample_previous_master_index
#define geometry_list_index mcctest_sample_geometry_list_index
#define intersection_time_table mcctest_sample_intersection_time_table
#define Volumes mcctest_sample_Volumes
#define Geometries mcctest_sample_Geometries
#define starting_lists mcctest_sample_starting_lists
#define r mcctest_sample_r
#define r_start mcctest_sample_r_start
#define v mcctest_sample_v
#define error_msg mcctest_sample_error_msg
#define component_error_msg mcctest_sample_component_error_msg
#define string_output mcctest_sample_string_output
#define number_of_volumes mcctest_sample_number_of_volumes
#define volume_index mcctest_sample_volume_index
#define process_index mcctest_sample_process_index
#define solutions mcctest_sample_solutions
#define max_number_of_processes mcctest_sample_max_number_of_processes
#define limit mcctest_sample_limit
#define solution mcctest_sample_solution
#define min_solution mcctest_sample_min_solution
#define min_volume mcctest_sample_min_volume
#define time_found mcctest_sample_time_found
#define intersection_time mcctest_sample_intersection_time
#define min_intersection_time mcctest_sample_min_intersection_time
#define process mcctest_sample_process
#define process_start mcctest_sample_process_start
#define my_trace mcctest_sample_my_trace
#define p_my_trace mcctest_sample_p_my_trace
#define my_trace_fraction_control mcctest_sample_my_trace_fraction_control
#define k mcctest_sample_k
#define k_new mcctest_sample_k_new
#define k_old mcctest_sample_k_old
#define v_length mcctest_sample_v_length
#define my_sum mcctest_sample_my_sum
#define my_sum_plus_abs mcctest_sample_my_sum_plus_abs
#define culmative_probability mcctest_sample_culmative_probability
#define mc_prop mcctest_sample_mc_prop
#define time_to_scattering mcctest_sample_time_to_scattering
#define length_to_scattering mcctest_sample_length_to_scattering
#define length_to_boundery mcctest_sample_length_to_boundery
#define time_to_boundery mcctest_sample_time_to_boundery
#define selected_process mcctest_sample_selected_process
#define scattering_event mcctest_sample_scattering_event
#define time_propagated_without_scattering mcctest_sample_time_propagated_without_scattering
#define a_next_volume_found mcctest_sample_a_next_volume_found
#define next_volume mcctest_sample_next_volume
#define next_volume_priority mcctest_sample_next_volume_priority
#define done mcctest_sample_done
#define current_volume mcctest_sample_current_volume
#define number_of_solutions mcctest_sample_number_of_solutions
#define number_of_solutions_static mcctest_sample_number_of_solutions_static
#define check mcctest_sample_check
#define start mcctest_sample_start
#define intersection_with_children mcctest_sample_intersection_with_children
#define geometry_output mcctest_sample_geometry_output
#define tree_next_volume mcctest_sample_tree_next_volume
#define pre_allocated1 mcctest_sample_pre_allocated1
#define pre_allocated2 mcctest_sample_pre_allocated2
#define pre_allocated3 mcctest_sample_pre_allocated3
#define ray_position mcctest_sample_ray_position
#define ray_velocity mcctest_sample_ray_velocity
#define ray_velocity_final mcctest_sample_ray_velocity_final
#define volume_0_found mcctest_sample_volume_0_found
#define scattered_flag mcctest_sample_scattered_flag
#define scattered_flag_VP mcctest_sample_scattered_flag_VP
#define master_transposed_rotation_matrix mcctest_sample_master_transposed_rotation_matrix
#define temp_rotation_matrix mcctest_sample_temp_rotation_matrix
#define non_rotated_position mcctest_sample_non_rotated_position
#define rotated_position mcctest_sample_rotated_position
#define enable_tagging mcctest_sample_enable_tagging
#define stop_tagging_ray mcctest_sample_stop_tagging_ray
#define stop_creating_nodes mcctest_sample_stop_creating_nodes
#define enable_tagging_check mcctest_sample_enable_tagging_check
#define master_tagging_node_list mcctest_sample_master_tagging_node_list
#define current_tagging_node mcctest_sample_current_tagging_node
#define tagging_leaf_counter mcctest_sample_tagging_leaf_counter
#define number_of_scattering_events mcctest_sample_number_of_scattering_events
#define real_transmission_probability mcctest_sample_real_transmission_probability
#define mc_transmission_probability mcctest_sample_mc_transmission_probability
#define number_of_masks mcctest_sample_number_of_masks
#define number_of_masked_volumes mcctest_sample_number_of_masked_volumes
#define need_to_run_within_which_volume mcctest_sample_need_to_run_within_which_volume
#define mask_index_main mcctest_sample_mask_index_main
#define mask_iterate mcctest_sample_mask_iterate
#define mask_status_list mcctest_sample_mask_status_list
#define current_mask_intersect_list_status mcctest_sample_current_mask_intersect_list_status
#define mask_volume_index_list mcctest_sample_mask_volume_index_list
#define geometry_component_index_list mcctest_sample_geometry_component_index_list
#define Volume_copies_allocated mcctest_sample_Volume_copies_allocated
#define p_old mcctest_sample_p_old
#define this_logger mcctest_sample_this_logger
#define conditional_status mcctest_sample_conditional_status
#define tagging_conditional_list mcctest_sample_tagging_conditional_list
#define free_tagging_conditioanl_list mcctest_sample_free_tagging_conditioanl_list
#define logger_conditional_extend_array mcctest_sample_logger_conditional_extend_array
#define tagging_conditional_extend mcctest_sample_tagging_conditional_extend
#define max_conditional_extend_index mcctest_sample_max_conditional_extend_index
#define safty_distance mcctest_sample_safty_distance
#define safty_distance2 mcctest_sample_safty_distance2
#define number_of_processes_array mcctest_sample_number_of_processes_array
#define temporary_focus_data mcctest_sample_temporary_focus_data
#define focus_data_index mcctest_sample_focus_data_index
{   /* Declarations of test_sample=Union_master() SETTING parameters. */
MCNUM allow_inside_start = mcctest_sample_allow_inside_start;
MCNUM history_limit = mcctest_sample_history_limit;
MCNUM enable_conditionals = mcctest_sample_enable_conditionals;
MCNUM inherit_number_of_scattering_events = mcctest_sample_inherit_number_of_scattering_events;
#line 1948 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../contrib/union/Union_master.comp"
{
// write out histories from tagging system if enabled
if (enable_tagging) {
    if (finally_verbal) printf("Writing tagging tree to disk \n");
    if (finally_verbal) printf("Number of leafs = %d \n",tagging_leaf_counter);
    // While writing the tagging tree to disk, all the leafs are deallocated
    write_tagging_tree(&master_tagging_node_list, Volumes, tagging_leaf_counter, number_of_volumes);
}
if (master_tagging_node_list.num_elements > 0) free(master_tagging_node_list.elements);




if (finally_verbal) printf("Freeing variables which are always allocated \n");
// free allocated arrays specific to this master union component
free(scattered_flag);
free(my_trace);
free(pre_allocated1);
free(pre_allocated2);
free(pre_allocated3);
free(number_of_processes_array);

if (finally_verbal) printf("Freeing intersection_time_table \n");
for (iterate = 1;iterate < intersection_time_table.num_volumes;iterate++){
    free(intersection_time_table.intersection_times[iterate]);
}

free(intersection_time_table.n_elements);
free(intersection_time_table.calculated);
free(intersection_time_table.intersection_times);

if (free_tagging_conditioanl_list == 1) free(tagging_conditional_list);

/*
if (tagging_conditional_list->num_elements > 0) {
  free(tagging_conditional_list.conditional_functions);
  free(tagging_conditional_list.p_data_unions);
}
*/

if (finally_verbal) printf("Freeing lists for individual volumes \n");
for (volume_index=0;volume_index<number_of_volumes;volume_index++) {
  
  if (Volumes[volume_index]->geometry.intersect_check_list.num_elements > 0) free(Volumes[volume_index]->geometry.intersect_check_list.elements);
  if (Volumes[volume_index]->geometry.destinations_list.num_elements > 0) free(Volumes[volume_index]->geometry.destinations_list.elements);
  if (Volumes[volume_index]->geometry.reduced_destinations_list.num_elements > 0) free(Volumes[volume_index]->geometry.reduced_destinations_list.elements);
  if (Volumes[volume_index]->geometry.children.num_elements > 0) free(Volumes[volume_index]->geometry.children.elements);
  if (Volumes[volume_index]->geometry.direct_children.num_elements > 0) free(Volumes[volume_index]->geometry.direct_children.elements);
  if (Volumes[volume_index]->geometry.masked_by_list.num_elements > 0) free(Volumes[volume_index]->geometry.masked_by_list.elements);
  if (Volumes[volume_index]->geometry.masked_by_mask_index_list.num_elements > 0) free(Volumes[volume_index]->geometry.masked_by_mask_index_list.elements);
  if (Volumes[volume_index]->geometry.mask_list.num_elements > 0) free(Volumes[volume_index]->geometry.mask_list.elements);
  if (Volumes[volume_index]->geometry.mask_intersect_list.num_elements > 0) free(Volumes[volume_index]->geometry.mask_intersect_list.elements);
  if (enable_tagging)
    if (Volumes[volume_index]->geometry.next_volume_list.num_elements > 0) free(Volumes[volume_index]->geometry.next_volume_list.elements);
  // Add dealocation of logging
  

  if (volume_index > 0) { // Volume 0 does not have physical properties allocated
    free(scattered_flag_VP[volume_index]);
    if (Volumes[volume_index]->geometry.process_rot_allocated == 1) {
          free(Volumes[volume_index]->geometry.process_rot_matrix_array);
          free(Volumes[volume_index]->geometry.transpose_process_rot_matrix_array);
          Volumes[volume_index]->geometry.process_rot_allocated = 0;
    }
    if (on_int_list(Volume_copies_allocated,volume_index))
      // This is a local copy of a volume, deallocate that local copy (all the allocated memory attachted to it was just deallocated, so this should not leave any leaks)
      free(Volumes[volume_index]);
    else
      // Only free p_physics for vacuum volumes for the original at the end (there is a p_physics allocated for each vacuum volume)
      if (Volumes[volume_index]->p_physics->is_vacuum == 1 ) free(Volumes[volume_index]->p_physics);
  }
  
  if (Volumes[volume_index]->loggers.num_elements >0) {
    for (iterate=0;iterate<Volumes[volume_index]->loggers.num_elements;iterate++) {
      free(Volumes[volume_index]->loggers.p_logger_volume[iterate].p_logger_process);
    }
    free(Volumes[volume_index]->loggers.p_logger_volume);
  }
  
}

free(scattered_flag_VP);

if (finally_verbal) printf("Freeing starting lists \n");
if (starting_lists.allowed_starting_volume_logic_list.num_elements > 0) free(starting_lists.allowed_starting_volume_logic_list.elements);
if (starting_lists.reduced_start_list.num_elements > 0) free(starting_lists.reduced_start_list.elements);
if (starting_lists.start_logic_list.num_elements > 0) free(starting_lists.start_logic_list.elements);

if (finally_verbal) printf("Freeing mask lists \n");
if (mask_status_list.num_elements>0) free(mask_status_list.elements);
if (current_mask_intersect_list_status.num_elements>0) free(current_mask_intersect_list_status.elements);
if (mask_volume_index_list.num_elements>0) free(mask_volume_index_list.elements);

if (finally_verbal) printf("Freeing component index list \n");
if (geometry_component_index_list.num_elements>0) free(geometry_component_index_list.elements);


if (finally_verbal) printf("Freeing Volumes \n");
free(Volumes);

// Free global allocated arrays if this is the last master union component in the instrument file

if (global_master_list.elements[global_master_list.num_elements-1].component_index == INDEX_CURRENT_COMP) {
    if (finally_verbal) printf("Freeing global arrays because this is the last Union master component\n");
    
    // Freeing lists allocated in Union_initialization
    //#ifdef PROCESS_DETECTOR
        if (finally_verbal) printf("Freeing global process list \n");
        if (global_process_list.num_elements > 0) free(global_process_list.elements);
    //#endif

    //#ifdef MATERIAL_DETECTOR
        if (finally_verbal) printf("Freeing global material list \n");
        if (global_material_list.num_elements > 0) free(global_material_list.elements);
    //#endif

    //#ifdef ANY_GEOMETRY_DETECTOR_DECLARE
        if (finally_verbal) printf("Freeing global geometry list \n");
        if (global_geometry_list.num_elements > 0) free(global_geometry_list.elements);
    //#endif
    
    //#ifdef MASTER_DETECTOR
        if (finally_verbal) printf("Freeing global master list \n");
        if (global_master_list.num_elements > 0) free(global_master_list.elements);
    //#endif
    
    
    //#ifdef UNION_LOGGER_DECLARE
    if (finally_verbal) printf("Freeing global logger lists \n");
    for (iterate=0;iterate<global_all_volume_logger_list.num_elements;iterate++) {
      if (global_all_volume_logger_list.elements[iterate].logger->conditional_list.num_elements > 0) {
        free(global_all_volume_logger_list.elements[iterate].logger->conditional_list.conditional_functions);
        free(global_all_volume_logger_list.elements[iterate].logger->conditional_list.p_data_unions);
      }
    }
    if (global_all_volume_logger_list.num_elements > 0) free(global_all_volume_logger_list.elements);
    
    
    for (iterate=0;iterate<global_specific_volumes_logger_list.num_elements;iterate++) {
      if (global_specific_volumes_logger_list.elements[iterate].logger->conditional_list.num_elements > 0) {
        free(global_specific_volumes_logger_list.elements[iterate].logger->conditional_list.conditional_functions);
        free(global_specific_volumes_logger_list.elements[iterate].logger->conditional_list.p_data_unions);
      }
    }
    if (global_specific_volumes_logger_list.num_elements > 0) free(global_specific_volumes_logger_list.elements);
    //#endif
    
    for (iterate=0;iterate<global_tagging_conditional_list.num_elements;iterate++) {
      if (global_tagging_conditional_list.elements[iterate].conditional_list.num_elements > 0) {
        free(global_tagging_conditional_list.elements[iterate].conditional_list.conditional_functions);
        free(global_tagging_conditional_list.elements[iterate].conditional_list.p_data_unions);
      }
    }
    if (global_tagging_conditional_list.num_elements>0) free(global_tagging_conditional_list.elements);
    
    /*
    if (finally_verbal) printf("Freeing global tagging conditional list \n");
    if (global_tagging_conditional_list.num_elements > 0) free(global_tagging_conditional_list.elements);
    */
}

}
#line 15303 "./Union_test_mask.c"
}   /* End of test_sample=Union_master() SETTING parameter declarations. */
#undef focus_data_index
#undef temporary_focus_data
#undef number_of_processes_array
#undef safty_distance2
#undef safty_distance
#undef max_conditional_extend_index
#undef tagging_conditional_extend
#undef logger_conditional_extend_array
#undef free_tagging_conditioanl_list
#undef tagging_conditional_list
#undef conditional_status
#undef this_logger
#undef p_old
#undef Volume_copies_allocated
#undef geometry_component_index_list
#undef mask_volume_index_list
#undef current_mask_intersect_list_status
#undef mask_status_list
#undef mask_iterate
#undef mask_index_main
#undef need_to_run_within_which_volume
#undef number_of_masked_volumes
#undef number_of_masks
#undef mc_transmission_probability
#undef real_transmission_probability
#undef number_of_scattering_events
#undef tagging_leaf_counter
#undef current_tagging_node
#undef master_tagging_node_list
#undef enable_tagging_check
#undef stop_creating_nodes
#undef stop_tagging_ray
#undef enable_tagging
#undef rotated_position
#undef non_rotated_position
#undef temp_rotation_matrix
#undef master_transposed_rotation_matrix
#undef scattered_flag_VP
#undef scattered_flag
#undef volume_0_found
#undef ray_velocity_final
#undef ray_velocity
#undef ray_position
#undef pre_allocated3
#undef pre_allocated2
#undef pre_allocated1
#undef tree_next_volume
#undef geometry_output
#undef intersection_with_children
#undef start
#undef check
#undef number_of_solutions_static
#undef number_of_solutions
#undef current_volume
#undef done
#undef next_volume_priority
#undef next_volume
#undef a_next_volume_found
#undef time_propagated_without_scattering
#undef scattering_event
#undef selected_process
#undef time_to_boundery
#undef length_to_boundery
#undef length_to_scattering
#undef time_to_scattering
#undef mc_prop
#undef culmative_probability
#undef my_sum_plus_abs
#undef my_sum
#undef v_length
#undef k_old
#undef k_new
#undef k
#undef my_trace_fraction_control
#undef p_my_trace
#undef my_trace
#undef process_start
#undef process
#undef min_intersection_time
#undef intersection_time
#undef time_found
#undef min_volume
#undef min_solution
#undef solution
#undef limit
#undef max_number_of_processes
#undef solutions
#undef process_index
#undef volume_index
#undef number_of_volumes
#undef string_output
#undef component_error_msg
#undef error_msg
#undef v
#undef r_start
#undef r
#undef starting_lists
#undef Geometries
#undef Volumes
#undef intersection_time_table
#undef geometry_list_index
#undef previous_master_index
#undef this_global_master_index
#undef global_master_element
#undef starting_volume_warning
#undef finally_verbal
#undef trace_verbal
#undef list_verbal
#undef verbal
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[12]) fprintf(stderr, "Warning: No neutron could reach Component[12] test_sample\n");
    if (mcAbsorbProp[12]) fprintf(stderr, "Warning: %g events were removed in Component[12] test_sample=Union_master()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[12]);
    if (!mcNCounter[13]) fprintf(stderr, "Warning: No neutron could reach Component[13] m4pi\n");
    if (mcAbsorbProp[13]) fprintf(stderr, "Warning: %g events were removed in Component[13] m4pi=PSD_monitor_4PI()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[13]);
    if (!mcNCounter[14]) fprintf(stderr, "Warning: No neutron could reach Component[14] m4pi_1\n");
    if (mcAbsorbProp[14]) fprintf(stderr, "Warning: %g events were removed in Component[14] m4pi_1=PSD_monitor_4PI()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[14]);
    if (!mcNCounter[15]) fprintf(stderr, "Warning: No neutron could reach Component[15] m4pi_2\n");
    if (mcAbsorbProp[15]) fprintf(stderr, "Warning: %g events were removed in Component[15] m4pi_2=PSD_monitor_4PI()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[15]);
    if (!mcNCounter[16]) fprintf(stderr, "Warning: No neutron could reach Component[16] arm_1\n");
    if (mcAbsorbProp[16]) fprintf(stderr, "Warning: %g events were removed in Component[16] arm_1=Arm()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[16]);
    if (!mcNCounter[17]) fprintf(stderr, "Warning: No neutron could reach Component[17] m4pi_both\n");
    if (mcAbsorbProp[17]) fprintf(stderr, "Warning: %g events were removed in Component[17] m4pi_both=PSD_monitor_4PI()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[17]);
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

  /* MCDISPLAY code for component 'a1'. */
  SIG_MESSAGE("a1 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "a1");
#define mccompcurname  a1
#define mccompcurtype  Progress_bar
#define mccompcurindex 3
#define IntermediateCnts mcca1_IntermediateCnts
#define StartTime mcca1_StartTime
#define EndTime mcca1_EndTime
#define CurrentTime mcca1_CurrentTime
{   /* Declarations of a1=Progress_bar() SETTING parameters. */
char* profile = mcca1_profile;
MCNUM percent = mcca1_percent;
MCNUM flag_save = mcca1_flag_save;
MCNUM minutes = mcca1_minutes;
#line 147 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../misc/Progress_bar.comp"
{
  
}
#line 15464 "./Union_test_mask.c"
}   /* End of a1=Progress_bar() SETTING parameter declarations. */
#undef CurrentTime
#undef EndTime
#undef StartTime
#undef IntermediateCnts
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'source'. */
  SIG_MESSAGE("source (McDisplay)");
  printf("MCDISPLAY: component %s\n", "source");
#define mccompcurname  source
#define mccompcurtype  Source_div
#define mccompcurindex 4
#define thetah mccsource_thetah
#define thetav mccsource_thetav
#define sigmah mccsource_sigmah
#define sigmav mccsource_sigmav
#define tan_h mccsource_tan_h
#define tan_v mccsource_tan_v
#define p_init mccsource_p_init
#define dist mccsource_dist
#define focus_xw mccsource_focus_xw
#define focus_yh mccsource_focus_yh
{   /* Declarations of source=Source_div() SETTING parameters. */
MCNUM xwidth = mccsource_xwidth;
MCNUM yheight = mccsource_yheight;
MCNUM focus_aw = mccsource_focus_aw;
MCNUM focus_ah = mccsource_focus_ah;
MCNUM E0 = mccsource_E0;
MCNUM dE = mccsource_dE;
MCNUM lambda0 = mccsource_lambda0;
MCNUM dlambda = mccsource_dlambda;
MCNUM gauss = mccsource_gauss;
MCNUM flux = mccsource_flux;
#line 167 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../sources/Source_div.comp"
{
  
  multiline(5, -xwidth/2.0, -yheight/2.0, 0.0,
                xwidth/2.0, -yheight/2.0, 0.0,
                xwidth/2.0,  yheight/2.0, 0.0,
               -xwidth/2.0,  yheight/2.0, 0.0,
               -xwidth/2.0, -yheight/2.0, 0.0);
  if (dist) {
    dashed_line(0,0,0, -focus_xw/2,-focus_yh/2,dist, 4);
    dashed_line(0,0,0,  focus_xw/2,-focus_yh/2,dist, 4);
    dashed_line(0,0,0,  focus_xw/2, focus_yh/2,dist, 4);
    dashed_line(0,0,0, -focus_xw/2, focus_yh/2,dist, 4);
  }
}
#line 15516 "./Union_test_mask.c"
}   /* End of source=Source_div() SETTING parameter declarations. */
#undef focus_yh
#undef focus_xw
#undef dist
#undef p_init
#undef tan_v
#undef tan_h
#undef sigmav
#undef sigmah
#undef thetav
#undef thetah
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'beam_center'. */
  SIG_MESSAGE("beam_center (McDisplay)");
  printf("MCDISPLAY: component %s\n", "beam_center");
#define mccompcurname  beam_center
#define mccompcurtype  Arm
#define mccompcurindex 5
#line 40 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Arm.comp"
{
  /* A bit ugly; hard-coded dimensions. */
  
  line(0,0,0,0.2,0,0);
  line(0,0,0,0,0.2,0);
  line(0,0,0,0,0,0.2);
}
#line 15546 "./Union_test_mask.c"
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'test_sample'. */
  SIG_MESSAGE("test_sample (McDisplay)");
  printf("MCDISPLAY: component %s\n", "test_sample");
#define mccompcurname  test_sample
#define mccompcurtype  Union_master
#define mccompcurindex 12
#define verbal mcctest_sample_verbal
#define list_verbal mcctest_sample_list_verbal
#define trace_verbal mcctest_sample_trace_verbal
#define finally_verbal mcctest_sample_finally_verbal
#define starting_volume_warning mcctest_sample_starting_volume_warning
#define global_master_element mcctest_sample_global_master_element
#define this_global_master_index mcctest_sample_this_global_master_index
#define previous_master_index mcctest_sample_previous_master_index
#define geometry_list_index mcctest_sample_geometry_list_index
#define intersection_time_table mcctest_sample_intersection_time_table
#define Volumes mcctest_sample_Volumes
#define Geometries mcctest_sample_Geometries
#define starting_lists mcctest_sample_starting_lists
#define r mcctest_sample_r
#define r_start mcctest_sample_r_start
#define v mcctest_sample_v
#define error_msg mcctest_sample_error_msg
#define component_error_msg mcctest_sample_component_error_msg
#define string_output mcctest_sample_string_output
#define number_of_volumes mcctest_sample_number_of_volumes
#define volume_index mcctest_sample_volume_index
#define process_index mcctest_sample_process_index
#define solutions mcctest_sample_solutions
#define max_number_of_processes mcctest_sample_max_number_of_processes
#define limit mcctest_sample_limit
#define solution mcctest_sample_solution
#define min_solution mcctest_sample_min_solution
#define min_volume mcctest_sample_min_volume
#define time_found mcctest_sample_time_found
#define intersection_time mcctest_sample_intersection_time
#define min_intersection_time mcctest_sample_min_intersection_time
#define process mcctest_sample_process
#define process_start mcctest_sample_process_start
#define my_trace mcctest_sample_my_trace
#define p_my_trace mcctest_sample_p_my_trace
#define my_trace_fraction_control mcctest_sample_my_trace_fraction_control
#define k mcctest_sample_k
#define k_new mcctest_sample_k_new
#define k_old mcctest_sample_k_old
#define v_length mcctest_sample_v_length
#define my_sum mcctest_sample_my_sum
#define my_sum_plus_abs mcctest_sample_my_sum_plus_abs
#define culmative_probability mcctest_sample_culmative_probability
#define mc_prop mcctest_sample_mc_prop
#define time_to_scattering mcctest_sample_time_to_scattering
#define length_to_scattering mcctest_sample_length_to_scattering
#define length_to_boundery mcctest_sample_length_to_boundery
#define time_to_boundery mcctest_sample_time_to_boundery
#define selected_process mcctest_sample_selected_process
#define scattering_event mcctest_sample_scattering_event
#define time_propagated_without_scattering mcctest_sample_time_propagated_without_scattering
#define a_next_volume_found mcctest_sample_a_next_volume_found
#define next_volume mcctest_sample_next_volume
#define next_volume_priority mcctest_sample_next_volume_priority
#define done mcctest_sample_done
#define current_volume mcctest_sample_current_volume
#define number_of_solutions mcctest_sample_number_of_solutions
#define number_of_solutions_static mcctest_sample_number_of_solutions_static
#define check mcctest_sample_check
#define start mcctest_sample_start
#define intersection_with_children mcctest_sample_intersection_with_children
#define geometry_output mcctest_sample_geometry_output
#define tree_next_volume mcctest_sample_tree_next_volume
#define pre_allocated1 mcctest_sample_pre_allocated1
#define pre_allocated2 mcctest_sample_pre_allocated2
#define pre_allocated3 mcctest_sample_pre_allocated3
#define ray_position mcctest_sample_ray_position
#define ray_velocity mcctest_sample_ray_velocity
#define ray_velocity_final mcctest_sample_ray_velocity_final
#define volume_0_found mcctest_sample_volume_0_found
#define scattered_flag mcctest_sample_scattered_flag
#define scattered_flag_VP mcctest_sample_scattered_flag_VP
#define master_transposed_rotation_matrix mcctest_sample_master_transposed_rotation_matrix
#define temp_rotation_matrix mcctest_sample_temp_rotation_matrix
#define non_rotated_position mcctest_sample_non_rotated_position
#define rotated_position mcctest_sample_rotated_position
#define enable_tagging mcctest_sample_enable_tagging
#define stop_tagging_ray mcctest_sample_stop_tagging_ray
#define stop_creating_nodes mcctest_sample_stop_creating_nodes
#define enable_tagging_check mcctest_sample_enable_tagging_check
#define master_tagging_node_list mcctest_sample_master_tagging_node_list
#define current_tagging_node mcctest_sample_current_tagging_node
#define tagging_leaf_counter mcctest_sample_tagging_leaf_counter
#define number_of_scattering_events mcctest_sample_number_of_scattering_events
#define real_transmission_probability mcctest_sample_real_transmission_probability
#define mc_transmission_probability mcctest_sample_mc_transmission_probability
#define number_of_masks mcctest_sample_number_of_masks
#define number_of_masked_volumes mcctest_sample_number_of_masked_volumes
#define need_to_run_within_which_volume mcctest_sample_need_to_run_within_which_volume
#define mask_index_main mcctest_sample_mask_index_main
#define mask_iterate mcctest_sample_mask_iterate
#define mask_status_list mcctest_sample_mask_status_list
#define current_mask_intersect_list_status mcctest_sample_current_mask_intersect_list_status
#define mask_volume_index_list mcctest_sample_mask_volume_index_list
#define geometry_component_index_list mcctest_sample_geometry_component_index_list
#define Volume_copies_allocated mcctest_sample_Volume_copies_allocated
#define p_old mcctest_sample_p_old
#define this_logger mcctest_sample_this_logger
#define conditional_status mcctest_sample_conditional_status
#define tagging_conditional_list mcctest_sample_tagging_conditional_list
#define free_tagging_conditioanl_list mcctest_sample_free_tagging_conditioanl_list
#define logger_conditional_extend_array mcctest_sample_logger_conditional_extend_array
#define tagging_conditional_extend mcctest_sample_tagging_conditional_extend
#define max_conditional_extend_index mcctest_sample_max_conditional_extend_index
#define safty_distance mcctest_sample_safty_distance
#define safty_distance2 mcctest_sample_safty_distance2
#define number_of_processes_array mcctest_sample_number_of_processes_array
#define temporary_focus_data mcctest_sample_temporary_focus_data
#define focus_data_index mcctest_sample_focus_data_index
{   /* Declarations of test_sample=Union_master() SETTING parameters. */
MCNUM allow_inside_start = mcctest_sample_allow_inside_start;
MCNUM history_limit = mcctest_sample_history_limit;
MCNUM enable_conditionals = mcctest_sample_enable_conditionals;
MCNUM inherit_number_of_scattering_events = mcctest_sample_inherit_number_of_scattering_events;
#line 2113 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../contrib/union/Union_master.comp"
{
  // mcdisplay is handled in the component files for each geometry and called here. The line function is only available in this section, and not through functions,
  //   so all the lines to be drawn for each volume are collected in a structure that is then drawn here.
  magnify("xyz");
  struct lines_to_draw lines_to_draw_master;
  for (volume_index = 1; volume_index < number_of_volumes; volume_index++) {
        if (Volumes[volume_index]->geometry.visualization_on == 1) {
            lines_to_draw_master.number_of_lines = 0;
            
            Volumes[volume_index]->geometry.mcdisplay_function(&lines_to_draw_master,volume_index,Geometries,number_of_volumes);
            
            for (iterate = 0;iterate<lines_to_draw_master.number_of_lines;iterate++) {
               if (lines_to_draw_master.lines[iterate].number_of_dashes == 1) {
                 line(lines_to_draw_master.lines[iterate].point1.x,lines_to_draw_master.lines[iterate].point1.y,lines_to_draw_master.lines[iterate].point1.z,lines_to_draw_master.lines[iterate].point2.x,lines_to_draw_master.lines[iterate].point2.y,lines_to_draw_master.lines[iterate].point2.z);
               }
               else {
                 dashed_line(lines_to_draw_master.lines[iterate].point1.x,lines_to_draw_master.lines[iterate].point1.y,lines_to_draw_master.lines[iterate].point1.z,lines_to_draw_master.lines[iterate].point2.x,lines_to_draw_master.lines[iterate].point2.y,lines_to_draw_master.lines[iterate].point2.z,lines_to_draw_master.lines[iterate].number_of_dashes);
               }
            }
            
            if (lines_to_draw_master.number_of_lines>0) free(lines_to_draw_master.lines);
        }
   }

}
#line 15697 "./Union_test_mask.c"
}   /* End of test_sample=Union_master() SETTING parameter declarations. */
#undef focus_data_index
#undef temporary_focus_data
#undef number_of_processes_array
#undef safty_distance2
#undef safty_distance
#undef max_conditional_extend_index
#undef tagging_conditional_extend
#undef logger_conditional_extend_array
#undef free_tagging_conditioanl_list
#undef tagging_conditional_list
#undef conditional_status
#undef this_logger
#undef p_old
#undef Volume_copies_allocated
#undef geometry_component_index_list
#undef mask_volume_index_list
#undef current_mask_intersect_list_status
#undef mask_status_list
#undef mask_iterate
#undef mask_index_main
#undef need_to_run_within_which_volume
#undef number_of_masked_volumes
#undef number_of_masks
#undef mc_transmission_probability
#undef real_transmission_probability
#undef number_of_scattering_events
#undef tagging_leaf_counter
#undef current_tagging_node
#undef master_tagging_node_list
#undef enable_tagging_check
#undef stop_creating_nodes
#undef stop_tagging_ray
#undef enable_tagging
#undef rotated_position
#undef non_rotated_position
#undef temp_rotation_matrix
#undef master_transposed_rotation_matrix
#undef scattered_flag_VP
#undef scattered_flag
#undef volume_0_found
#undef ray_velocity_final
#undef ray_velocity
#undef ray_position
#undef pre_allocated3
#undef pre_allocated2
#undef pre_allocated1
#undef tree_next_volume
#undef geometry_output
#undef intersection_with_children
#undef start
#undef check
#undef number_of_solutions_static
#undef number_of_solutions
#undef current_volume
#undef done
#undef next_volume_priority
#undef next_volume
#undef a_next_volume_found
#undef time_propagated_without_scattering
#undef scattering_event
#undef selected_process
#undef time_to_boundery
#undef length_to_boundery
#undef length_to_scattering
#undef time_to_scattering
#undef mc_prop
#undef culmative_probability
#undef my_sum_plus_abs
#undef my_sum
#undef v_length
#undef k_old
#undef k_new
#undef k
#undef my_trace_fraction_control
#undef p_my_trace
#undef my_trace
#undef process_start
#undef process
#undef min_intersection_time
#undef intersection_time
#undef time_found
#undef min_volume
#undef min_solution
#undef solution
#undef limit
#undef max_number_of_processes
#undef solutions
#undef process_index
#undef volume_index
#undef number_of_volumes
#undef string_output
#undef component_error_msg
#undef error_msg
#undef v
#undef r_start
#undef r
#undef starting_lists
#undef Geometries
#undef Volumes
#undef intersection_time_table
#undef geometry_list_index
#undef previous_master_index
#undef this_global_master_index
#undef global_master_element
#undef starting_volume_warning
#undef finally_verbal
#undef trace_verbal
#undef list_verbal
#undef verbal
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'm4pi'. */
  SIG_MESSAGE("m4pi (McDisplay)");
  printf("MCDISPLAY: component %s\n", "m4pi");
#define mccompcurname  m4pi
#define mccompcurtype  PSD_monitor_4PI
#define mccompcurindex 13
#define nx mccm4pi_nx
#define ny mccm4pi_ny
#define PSD_N mccm4pi_PSD_N
#define PSD_p mccm4pi_PSD_p
#define PSD_p2 mccm4pi_PSD_p2
{   /* Declarations of m4pi=PSD_monitor_4PI() SETTING parameters. */
char* filename = mccm4pi_filename;
MCNUM radius = mccm4pi_radius;
MCNUM restore_neutron = mccm4pi_restore_neutron;
int nowritefile = mccm4pi_nowritefile;
#line 121 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor_4PI.comp"
{
  
  circle("xy",0,0,0,radius);
  circle("xz",0,0,0,radius);
  circle("yz",0,0,0,radius);
}
#line 15835 "./Union_test_mask.c"
}   /* End of m4pi=PSD_monitor_4PI() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef ny
#undef nx
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'm4pi_1'. */
  SIG_MESSAGE("m4pi_1 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "m4pi_1");
#define mccompcurname  m4pi_1
#define mccompcurtype  PSD_monitor_4PI
#define mccompcurindex 14
#define nx mccm4pi_1_nx
#define ny mccm4pi_1_ny
#define PSD_N mccm4pi_1_PSD_N
#define PSD_p mccm4pi_1_PSD_p
#define PSD_p2 mccm4pi_1_PSD_p2
{   /* Declarations of m4pi_1=PSD_monitor_4PI() SETTING parameters. */
char* filename = mccm4pi_1_filename;
MCNUM radius = mccm4pi_1_radius;
MCNUM restore_neutron = mccm4pi_1_restore_neutron;
int nowritefile = mccm4pi_1_nowritefile;
#line 121 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor_4PI.comp"
{
  
  circle("xy",0,0,0,radius);
  circle("xz",0,0,0,radius);
  circle("yz",0,0,0,radius);
}
#line 15869 "./Union_test_mask.c"
}   /* End of m4pi_1=PSD_monitor_4PI() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef ny
#undef nx
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'm4pi_2'. */
  SIG_MESSAGE("m4pi_2 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "m4pi_2");
#define mccompcurname  m4pi_2
#define mccompcurtype  PSD_monitor_4PI
#define mccompcurindex 15
#define nx mccm4pi_2_nx
#define ny mccm4pi_2_ny
#define PSD_N mccm4pi_2_PSD_N
#define PSD_p mccm4pi_2_PSD_p
#define PSD_p2 mccm4pi_2_PSD_p2
{   /* Declarations of m4pi_2=PSD_monitor_4PI() SETTING parameters. */
char* filename = mccm4pi_2_filename;
MCNUM radius = mccm4pi_2_radius;
MCNUM restore_neutron = mccm4pi_2_restore_neutron;
int nowritefile = mccm4pi_2_nowritefile;
#line 121 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor_4PI.comp"
{
  
  circle("xy",0,0,0,radius);
  circle("xz",0,0,0,radius);
  circle("yz",0,0,0,radius);
}
#line 15903 "./Union_test_mask.c"
}   /* End of m4pi_2=PSD_monitor_4PI() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef ny
#undef nx
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'arm_1'. */
  SIG_MESSAGE("arm_1 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "arm_1");
#define mccompcurname  arm_1
#define mccompcurtype  Arm
#define mccompcurindex 16
#line 40 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Arm.comp"
{
  /* A bit ugly; hard-coded dimensions. */
  
  line(0,0,0,0.2,0,0);
  line(0,0,0,0,0.2,0);
  line(0,0,0,0,0,0.2);
}
#line 15928 "./Union_test_mask.c"
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'm4pi_both'. */
  SIG_MESSAGE("m4pi_both (McDisplay)");
  printf("MCDISPLAY: component %s\n", "m4pi_both");
#define mccompcurname  m4pi_both
#define mccompcurtype  PSD_monitor_4PI
#define mccompcurindex 17
#define nx mccm4pi_both_nx
#define ny mccm4pi_both_ny
#define PSD_N mccm4pi_both_PSD_N
#define PSD_p mccm4pi_both_PSD_p
#define PSD_p2 mccm4pi_both_PSD_p2
{   /* Declarations of m4pi_both=PSD_monitor_4PI() SETTING parameters. */
char* filename = mccm4pi_both_filename;
MCNUM radius = mccm4pi_both_radius;
MCNUM restore_neutron = mccm4pi_both_restore_neutron;
int nowritefile = mccm4pi_both_nowritefile;
#line 121 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor_4PI.comp"
{
  
  circle("xy",0,0,0,radius);
  circle("xz",0,0,0,radius);
  circle("yz",0,0,0,radius);
}
#line 15956 "./Union_test_mask.c"
}   /* End of m4pi_both=PSD_monitor_4PI() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef ny
#undef nx
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
/* end of generated C code ./Union_test_mask.c */
