/* Automatically generated file. Do not edit. 
 * Format:     ANSI C source code
 * Creator:    McStas <http://www.mcstas.org>
 * Instrument: ISIS_CRISP.instr (ISIS_CRISP)
 * Date:       Sat Oct 12 09:47:53 2019
 * File:       ISIS_CRISP.c
 * Compile:    cc -o ISIS_CRISP.out ISIS_CRISP.c  -lgsl -lgslcblas -L@MCCODE_LIB@/miniconda3/lib -I@MCCODE_LIB@/miniconda3/include
 * CFLAGS= -lgsl -lgslcblas -L@MCCODE_LIB@/miniconda3/lib -I@MCCODE_LIB@/miniconda3/include
 */


#define MCCODE_STRING "McStas 2.6rc1 - Oct. 12, 2019"
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
* Release: McStas 2.6rc1
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
#define MCCODE_STRING "McStas 2.6rc1 - Oct. 12, 2019"
#endif

#ifndef MCCODE_DATE
#define MCCODE_DATE "Oct. 12, 2019"
#endif

#ifndef MCCODE_VERSION
#define MCCODE_VERSION "2.6rc1"
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

#line 712 "ISIS_CRISP.c"

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

#line 945 "ISIS_CRISP.c"

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

  if (!detector.p1 || !detector.m || detector.filename[0] == '\0')
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

#line 4977 "ISIS_CRISP.c"

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

#line 5337 "ISIS_CRISP.c"
#ifdef MC_TRACE_ENABLED
int mctraceenabled = 1;
#else
int mctraceenabled = 0;
#endif
#define MCSTAS "/usr/share/mcstas/2.6rc1/"
int mcdefaultmain = 1;
char mcinstrument_name[] = "ISIS_CRISP";
char mcinstrument_source[] = "ISIS_CRISP.instr";
char *mcinstrument_exe=NULL; /* will be set to argv[0] in main */
int main(int argc, char *argv[]){return mccode_main(argc, argv);}
void mcinit(void);
void mcraytrace(void);
void mcsave(FILE *);
void mcfinally(void);
void mcdisplay(void);

/* Shared user declarations for all components 'ISIS_moderator'. */
#line 64 "/usr/share/mcstas/2.6rc1/contrib/ISIS_moderator.comp"
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
#line 6104 "ISIS_CRISP.c"

/* Shared user declarations for all components 'PSD_monitor'. */
#line 57 "/usr/share/mcstas/2.6rc1/monitors/PSD_monitor.comp"

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

#line 6190 "ISIS_CRISP.c"

/* Shared user declarations for all components 'Guide_tapering'. */
#line 91 "/usr/share/mcstas/2.6rc1/optics/Guide_tapering.comp"
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

#line 7784 "ISIS_CRISP.c"

/* Shared user declarations for all components 'Multilayer_Sample'. */
#line 60 "/usr/share/mcstas/2.6rc1/contrib/Multilayer_Sample.comp"
#ifndef GSL_VERSION
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_cblas.h>
#endif
#line 7794 "ISIS_CRISP.c"

/* Instrument parameters. */
MCNUM mcipglen;
MCNUM mcipflen;
MCNUM mcipw1;
MCNUM mcipvm;
MCNUM mcipFRAC;

#define mcNUMIPAR 5
int mcnumipar = 5;
struct mcinputtable_struct mcinputtable[mcNUMIPAR+1] = {
  "glen", &mcipglen, instr_type_double, "1.4", 
  "flen", &mcipflen, instr_type_double, "0.4", 
  "w1", &mcipw1, instr_type_double, "0.05", 
  "vm", &mcipvm, instr_type_double, "3.0", 
  "FRAC", &mcipFRAC, instr_type_double, "0", 
  NULL, NULL, instr_type_double, ""
};

/* User declarations from instrument definition. */
#define mccompcurname  ISIS_CRISP
#define mccompcurtype  INSTRUMENT
#define mccompcurindex 0
#define mcposaISIS_CRISP coords_set(0,0,0)
#define glen mcipglen
#define flen mcipflen
#define w1 mcipw1
#define vm mcipvm
#define FRAC mcipFRAC
#line 45 "ISIS_CRISP.instr"
  double glen,flen,spos,tend,t15,fp1;
  double linw,loutw,l,u1,u2,div1,b_ell_q,w1,w12,a_ell_q,lbw,w2;
#line 7827 "ISIS_CRISP.c"
#undef FRAC
#undef vm
#undef w1
#undef flen
#undef glen
#undef mcposaISIS_CRISP
#undef mccompcurindex
#undef mccompcurtype
#undef mccompcurname

/* neutron state table at each component input (local coords) */
/* [x, y, z, vx, vy, vz, t, sx, sy, sz, p] */
MCNUM mccomp_storein[11*20];
/* Components position table (absolute and relative coords) */
Coords mccomp_posa[20];
Coords mccomp_posr[20];
/* Counter for each comp to check for inactive ones */
MCNUM  mcNCounter[20];
MCNUM  mcPCounter[20];
MCNUM  mcP2Counter[20];
#define mcNUMCOMP 19 /* number of components */
/* Counter for PROP ABSORB */
MCNUM  mcAbsorbProp[20];
/* Flag true when previous component acted on the neutron (SCATTER) */
MCNUM mcScattered=0;
/* Flag true when neutron should be restored (RESTORE) */
MCNUM mcRestore=0;
/* Declarations of component definition and setting parameters. */

/* Definition parameters for component 'isis_source' [2]. */
#define mccisis_source_Face "crisp" /* declared as a string. May produce warnings at compile */
/* Setting parameters for component 'isis_source' [2]. */
MCNUM mccisis_source_Emin;
MCNUM mccisis_source_Emax;
MCNUM mccisis_source_dist;
MCNUM mccisis_source_focus_xw;
MCNUM mccisis_source_focus_yh;
MCNUM mccisis_source_xwidth;
MCNUM mccisis_source_yheight;
MCNUM mccisis_source_CAngle;
MCNUM mccisis_source_SAC;
MCNUM mccisis_source_Lmin;
MCNUM mccisis_source_Lmax;
int mccisis_source_target_index;
MCNUM mccisis_source_verbose;

/* Setting parameters for component 'defslit1' [3]. */
MCNUM mccdefslit1_xmin;
MCNUM mccdefslit1_xmax;
MCNUM mccdefslit1_ymin;
MCNUM mccdefslit1_ymax;
MCNUM mccdefslit1_radius;
MCNUM mccdefslit1_xwidth;
MCNUM mccdefslit1_yheight;

/* Setting parameters for component 'defslit2' [4]. */
MCNUM mccdefslit2_xmin;
MCNUM mccdefslit2_xmax;
MCNUM mccdefslit2_ymin;
MCNUM mccdefslit2_ymax;
MCNUM mccdefslit2_radius;
MCNUM mccdefslit2_xwidth;
MCNUM mccdefslit2_yheight;

/* Setting parameters for component 'coarseslit1' [5]. */
MCNUM mcccoarseslit1_xmin;
MCNUM mcccoarseslit1_xmax;
MCNUM mcccoarseslit1_ymin;
MCNUM mcccoarseslit1_ymax;
MCNUM mcccoarseslit1_radius;
MCNUM mcccoarseslit1_xwidth;
MCNUM mcccoarseslit1_yheight;

/* Setting parameters for component 'coarseslit2' [6]. */
MCNUM mcccoarseslit2_xmin;
MCNUM mcccoarseslit2_xmax;
MCNUM mcccoarseslit2_ymin;
MCNUM mcccoarseslit2_ymax;
MCNUM mcccoarseslit2_radius;
MCNUM mcccoarseslit2_xwidth;
MCNUM mcccoarseslit2_yheight;

/* Definition parameters for component 'lmon1' [7]. */
#define mcclmon1_nL 500
/* Setting parameters for component 'lmon1' [7]. */
char mcclmon1_filename[16384];
MCNUM mcclmon1_xmin;
MCNUM mcclmon1_xmax;
MCNUM mcclmon1_ymin;
MCNUM mcclmon1_ymax;
MCNUM mcclmon1_xwidth;
MCNUM mcclmon1_yheight;
MCNUM mcclmon1_Lmin;
MCNUM mcclmon1_Lmax;
MCNUM mcclmon1_restore_neutron;
int mcclmon1_nowritefile;

/* Setting parameters for component 'PSDmon1' [8]. */
int mccPSDmon1_nx;
int mccPSDmon1_ny;
char mccPSDmon1_filename[16384];
MCNUM mccPSDmon1_xmin;
MCNUM mccPSDmon1_xmax;
MCNUM mccPSDmon1_ymin;
MCNUM mccPSDmon1_ymax;
MCNUM mccPSDmon1_xwidth;
MCNUM mccPSDmon1_yheight;
MCNUM mccPSDmon1_restore_neutron;

/* Setting parameters for component 'slit1' [9]. */
MCNUM mccslit1_xmin;
MCNUM mccslit1_xmax;
MCNUM mccslit1_ymin;
MCNUM mccslit1_ymax;
MCNUM mccslit1_radius;
MCNUM mccslit1_xwidth;
MCNUM mccslit1_yheight;

/* Setting parameters for component 'PSDmons1' [10]. */
int mccPSDmons1_nx;
int mccPSDmons1_ny;
char mccPSDmons1_filename[16384];
MCNUM mccPSDmons1_xmin;
MCNUM mccPSDmons1_xmax;
MCNUM mccPSDmons1_ymin;
MCNUM mccPSDmons1_ymax;
MCNUM mccPSDmons1_xwidth;
MCNUM mccPSDmons1_yheight;
MCNUM mccPSDmons1_restore_neutron;

/* Setting parameters for component 'eguide1' [11]. */
char mcceguide1_option[16384];
MCNUM mcceguide1_w1;
MCNUM mcceguide1_h1;
MCNUM mcceguide1_l;
MCNUM mcceguide1_linw;
MCNUM mcceguide1_loutw;
MCNUM mcceguide1_linh;
MCNUM mcceguide1_louth;
MCNUM mcceguide1_R0;
MCNUM mcceguide1_Qcx;
MCNUM mcceguide1_Qcy;
MCNUM mcceguide1_alphax;
MCNUM mcceguide1_alphay;
MCNUM mcceguide1_W;
MCNUM mcceguide1_mx;
MCNUM mcceguide1_my;
MCNUM mcceguide1_segno;
MCNUM mcceguide1_curvature;
MCNUM mcceguide1_curvature_v;

/* Setting parameters for component 'slit2' [12]. */
MCNUM mccslit2_xmin;
MCNUM mccslit2_xmax;
MCNUM mccslit2_ymin;
MCNUM mccslit2_ymax;
MCNUM mccslit2_radius;
MCNUM mccslit2_xwidth;
MCNUM mccslit2_yheight;

/* Definition parameters for component 'lmon2' [13]. */
#define mcclmon2_nL 500
/* Setting parameters for component 'lmon2' [13]. */
char mcclmon2_filename[16384];
MCNUM mcclmon2_xmin;
MCNUM mcclmon2_xmax;
MCNUM mcclmon2_ymin;
MCNUM mcclmon2_ymax;
MCNUM mcclmon2_xwidth;
MCNUM mcclmon2_yheight;
MCNUM mcclmon2_Lmin;
MCNUM mcclmon2_Lmax;
MCNUM mcclmon2_restore_neutron;
int mcclmon2_nowritefile;

/* Definition parameters for component 'samp1' [14]. */
#define mccsamp1_xwidth 0.05
#define mccsamp1_zlength 0.15
#define mccsamp1_nlayer 0
#define mccsamp1_sldPar { 0.0 , 6.35e-6 }
#define mccsamp1_dPar { 0.0 }
#define mccsamp1_sigmaPar { 5.0 }
#define mccsamp1_frac_inc mcipFRAC
#define mccsamp1_ythick 0.01
#define mccsamp1_mu_inc 0.138
#define mccsamp1_target_index 1
#define mccsamp1_focus_xw 2 * tend
#define mccsamp1_focus_yh 0.01

/* Setting parameters for component 'slit3' [15]. */
MCNUM mccslit3_xmin;
MCNUM mccslit3_xmax;
MCNUM mccslit3_ymin;
MCNUM mccslit3_ymax;
MCNUM mccslit3_radius;
MCNUM mccslit3_xwidth;
MCNUM mccslit3_yheight;

/* Setting parameters for component 'eguide2' [16]. */
char mcceguide2_option[16384];
MCNUM mcceguide2_w1;
MCNUM mcceguide2_h1;
MCNUM mcceguide2_l;
MCNUM mcceguide2_linw;
MCNUM mcceguide2_loutw;
MCNUM mcceguide2_linh;
MCNUM mcceguide2_louth;
MCNUM mcceguide2_R0;
MCNUM mcceguide2_Qcx;
MCNUM mcceguide2_Qcy;
MCNUM mcceguide2_alphax;
MCNUM mcceguide2_alphay;
MCNUM mcceguide2_W;
MCNUM mcceguide2_mx;
MCNUM mcceguide2_my;
MCNUM mcceguide2_segno;
MCNUM mcceguide2_curvature;
MCNUM mcceguide2_curvature_v;

/* Definition parameters for component 'lmon3' [17]. */
#define mcclmon3_nL 500
/* Setting parameters for component 'lmon3' [17]. */
char mcclmon3_filename[16384];
MCNUM mcclmon3_xmin;
MCNUM mcclmon3_xmax;
MCNUM mcclmon3_ymin;
MCNUM mcclmon3_ymax;
MCNUM mcclmon3_xwidth;
MCNUM mcclmon3_yheight;
MCNUM mcclmon3_Lmin;
MCNUM mcclmon3_Lmax;
MCNUM mcclmon3_restore_neutron;
int mcclmon3_nowritefile;

/* Setting parameters for component 'PSDdet1' [18]. */
int mccPSDdet1_nx;
int mccPSDdet1_ny;
char mccPSDdet1_filename[16384];
MCNUM mccPSDdet1_xmin;
MCNUM mccPSDdet1_xmax;
MCNUM mccPSDdet1_ymin;
MCNUM mccPSDdet1_ymax;
MCNUM mccPSDdet1_xwidth;
MCNUM mccPSDdet1_yheight;
MCNUM mccPSDdet1_restore_neutron;

/* User component declarations. */

/* User declarations for component 'a1' [1]. */
#define mccompcurname  a1
#define mccompcurtype  Arm
#define mccompcurindex 1
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'isis_source' [2]. */
#define mccompcurname  isis_source
#define mccompcurtype  ISIS_moderator
#define mccompcurindex 2
#define Face mccisis_source_Face
#define p_in mccisis_source_p_in
#define Tnpts mccisis_source_Tnpts
#define scaleSize mccisis_source_scaleSize
#define angleArea mccisis_source_angleArea
#define Nsim mccisis_source_Nsim
#define Ncount mccisis_source_Ncount
#define TS mccisis_source_TS
#define rtE0 mccisis_source_rtE0
#define rtE1 mccisis_source_rtE1
#define rtmodX mccisis_source_rtmodX
#define rtmodY mccisis_source_rtmodY
#define TargetStation mccisis_source_TargetStation
#define CurrentWeight mccisis_source_CurrentWeight
#define Emin mccisis_source_Emin
#define Emax mccisis_source_Emax
#define dist mccisis_source_dist
#define focus_xw mccisis_source_focus_xw
#define focus_yh mccisis_source_focus_yh
#define xwidth mccisis_source_xwidth
#define yheight mccisis_source_yheight
#define CAngle mccisis_source_CAngle
#define SAC mccisis_source_SAC
#define Lmin mccisis_source_Lmin
#define Lmax mccisis_source_Lmax
#define target_index mccisis_source_target_index
#define verbose mccisis_source_verbose
#line 815 "/usr/share/mcstas/2.6rc1/contrib/ISIS_moderator.comp"
  #include <ctype.h>
  /* global variables */

  double p_in;         /* Polorization term (from McSTAS) */
  int Tnpts;           /* Number of points in parameteriation */
  double scaleSize;    /* correction for the actual area of the moderator viewed */
  double angleArea;    /* Area seen by the window  */
  double Nsim;	       /* Total number of neutrons to be simulated */
  int Ncount;          /* Number of neutron simulate so far*/
  Source TS;

  /* runtime variables*/

  double rtE0,rtE1;       /* runtime Energy minima and maxima so we can use angstroms as negative input */
  double rtmodX,rtmodY;    /* runtime moderator sizes, so that a negative argument may give a default size */

  int TargetStation;
  double CurrentWeight;
#line 8134 "ISIS_CRISP.c"
#undef verbose
#undef target_index
#undef Lmax
#undef Lmin
#undef SAC
#undef CAngle
#undef yheight
#undef xwidth
#undef focus_yh
#undef focus_xw
#undef dist
#undef Emax
#undef Emin
#undef CurrentWeight
#undef TargetStation
#undef rtmodY
#undef rtmodX
#undef rtE1
#undef rtE0
#undef TS
#undef Ncount
#undef Nsim
#undef angleArea
#undef scaleSize
#undef Tnpts
#undef p_in
#undef Face
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'defslit1' [3]. */
#define mccompcurname  defslit1
#define mccompcurtype  Slit
#define mccompcurindex 3
#define xmin mccdefslit1_xmin
#define xmax mccdefslit1_xmax
#define ymin mccdefslit1_ymin
#define ymax mccdefslit1_ymax
#define radius mccdefslit1_radius
#define xwidth mccdefslit1_xwidth
#define yheight mccdefslit1_yheight
#undef yheight
#undef xwidth
#undef radius
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'defslit2' [4]. */
#define mccompcurname  defslit2
#define mccompcurtype  Slit
#define mccompcurindex 4
#define xmin mccdefslit2_xmin
#define xmax mccdefslit2_xmax
#define ymin mccdefslit2_ymin
#define ymax mccdefslit2_ymax
#define radius mccdefslit2_radius
#define xwidth mccdefslit2_xwidth
#define yheight mccdefslit2_yheight
#undef yheight
#undef xwidth
#undef radius
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'coarseslit1' [5]. */
#define mccompcurname  coarseslit1
#define mccompcurtype  Slit
#define mccompcurindex 5
#define xmin mcccoarseslit1_xmin
#define xmax mcccoarseslit1_xmax
#define ymin mcccoarseslit1_ymin
#define ymax mcccoarseslit1_ymax
#define radius mcccoarseslit1_radius
#define xwidth mcccoarseslit1_xwidth
#define yheight mcccoarseslit1_yheight
#undef yheight
#undef xwidth
#undef radius
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'coarseslit2' [6]. */
#define mccompcurname  coarseslit2
#define mccompcurtype  Slit
#define mccompcurindex 6
#define xmin mcccoarseslit2_xmin
#define xmax mcccoarseslit2_xmax
#define ymin mcccoarseslit2_ymin
#define ymax mcccoarseslit2_ymax
#define radius mcccoarseslit2_radius
#define xwidth mcccoarseslit2_xwidth
#define yheight mcccoarseslit2_yheight
#undef yheight
#undef xwidth
#undef radius
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'lmon1' [7]. */
#define mccompcurname  lmon1
#define mccompcurtype  L_monitor
#define mccompcurindex 7
#define nL mcclmon1_nL
#define L_N mcclmon1_L_N
#define L_p mcclmon1_L_p
#define L_p2 mcclmon1_L_p2
#define filename mcclmon1_filename
#define xmin mcclmon1_xmin
#define xmax mcclmon1_xmax
#define ymin mcclmon1_ymin
#define ymax mcclmon1_ymax
#define xwidth mcclmon1_xwidth
#define yheight mcclmon1_yheight
#define Lmin mcclmon1_Lmin
#define Lmax mcclmon1_Lmax
#define restore_neutron mcclmon1_restore_neutron
#define nowritefile mcclmon1_nowritefile
#line 57 "/usr/share/mcstas/2.6rc1/monitors/L_monitor.comp"
double L_N[nL];
double L_p[nL], L_p2[nL];
#line 8276 "ISIS_CRISP.c"
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

/* User declarations for component 'PSDmon1' [8]. */
#define mccompcurname  PSDmon1
#define mccompcurtype  PSD_monitor
#define mccompcurindex 8
#define PSD_N mccPSDmon1_PSD_N
#define PSD_p mccPSDmon1_PSD_p
#define PSD_p2 mccPSDmon1_PSD_p2
#define nx mccPSDmon1_nx
#define ny mccPSDmon1_ny
#define filename mccPSDmon1_filename
#define xmin mccPSDmon1_xmin
#define xmax mccPSDmon1_xmax
#define ymin mccPSDmon1_ymin
#define ymax mccPSDmon1_ymax
#define xwidth mccPSDmon1_xwidth
#define yheight mccPSDmon1_yheight
#define restore_neutron mccPSDmon1_restore_neutron
#line 62 "/usr/share/mcstas/2.6rc1/monitors/PSD_monitor.comp"
  DArray2d PSD_N;
  DArray2d PSD_p;
  DArray2d PSD_p2;
#line 8317 "ISIS_CRISP.c"
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

/* User declarations for component 'slit1' [9]. */
#define mccompcurname  slit1
#define mccompcurtype  Slit
#define mccompcurindex 9
#define xmin mccslit1_xmin
#define xmax mccslit1_xmax
#define ymin mccslit1_ymin
#define ymax mccslit1_ymax
#define radius mccslit1_radius
#define xwidth mccslit1_xwidth
#define yheight mccslit1_yheight
#undef yheight
#undef xwidth
#undef radius
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'PSDmons1' [10]. */
#define mccompcurname  PSDmons1
#define mccompcurtype  PSD_monitor
#define mccompcurindex 10
#define PSD_N mccPSDmons1_PSD_N
#define PSD_p mccPSDmons1_PSD_p
#define PSD_p2 mccPSDmons1_PSD_p2
#define nx mccPSDmons1_nx
#define ny mccPSDmons1_ny
#define filename mccPSDmons1_filename
#define xmin mccPSDmons1_xmin
#define xmax mccPSDmons1_xmax
#define ymin mccPSDmons1_ymin
#define ymax mccPSDmons1_ymax
#define xwidth mccPSDmons1_xwidth
#define yheight mccPSDmons1_yheight
#define restore_neutron mccPSDmons1_restore_neutron
#line 62 "/usr/share/mcstas/2.6rc1/monitors/PSD_monitor.comp"
  DArray2d PSD_N;
  DArray2d PSD_p;
  DArray2d PSD_p2;
#line 8378 "ISIS_CRISP.c"
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

/* User declarations for component 'eguide1' [11]. */
#define mccompcurname  eguide1
#define mccompcurtype  Guide_tapering
#define mccompcurindex 11
#define w1c mcceguide1_w1c
#define w2c mcceguide1_w2c
#define ww mcceguide1_ww
#define hh mcceguide1_hh
#define whalf mcceguide1_whalf
#define hhalf mcceguide1_hhalf
#define lwhalf mcceguide1_lwhalf
#define lhhalf mcceguide1_lhhalf
#define h1_in mcceguide1_h1_in
#define h2_out mcceguide1_h2_out
#define w1_in mcceguide1_w1_in
#define w2_out mcceguide1_w2_out
#define l_seg mcceguide1_l_seg
#define seg mcceguide1_seg
#define h12 mcceguide1_h12
#define h2 mcceguide1_h2
#define w12 mcceguide1_w12
#define w2 mcceguide1_w2
#define a_ell_q mcceguide1_a_ell_q
#define b_ell_q mcceguide1_b_ell_q
#define lbw mcceguide1_lbw
#define lbh mcceguide1_lbh
#define mxi mcceguide1_mxi
#define u1 mcceguide1_u1
#define u2 mcceguide1_u2
#define div1 mcceguide1_div1
#define p2_para mcceguide1_p2_para
#define test mcceguide1_test
#define Div1 mcceguide1_Div1
#define i mcceguide1_i
#define ii mcceguide1_ii
#define seg mcceguide1_seg
#define fu mcceguide1_fu
#define pos mcceguide1_pos
#define file_name mcceguide1_file_name
#define ep mcceguide1_ep
#define num mcceguide1_num
#define rotation_h mcceguide1_rotation_h
#define rotation_v mcceguide1_rotation_v
#define option mcceguide1_option
#define w1 mcceguide1_w1
#define h1 mcceguide1_h1
#define l mcceguide1_l
#define linw mcceguide1_linw
#define loutw mcceguide1_loutw
#define linh mcceguide1_linh
#define louth mcceguide1_louth
#define R0 mcceguide1_R0
#define Qcx mcceguide1_Qcx
#define Qcy mcceguide1_Qcy
#define alphax mcceguide1_alphax
#define alphay mcceguide1_alphay
#define W mcceguide1_W
#define mx mcceguide1_mx
#define my mcceguide1_my
#define segno mcceguide1_segno
#define curvature mcceguide1_curvature
#define curvature_v mcceguide1_curvature_v
#line 96 "/usr/share/mcstas/2.6rc1/optics/Guide_tapering.comp"
double *w1c;
double *w2c;
double *ww, *hh;
double *whalf, *hhalf;
double *lwhalf, *lhhalf;
double *h1_in, *h2_out, *w1_in, *w2_out;
double l_seg, h12, h2, w12, w2, a_ell_q, b_ell_q, lbw, lbh;
double mxi ,u1 ,u2 ,div1, p2_para;
double test,Div1;
int i,ii,seg;
char *fu;
char *pos;
char file_name[1024];
char *ep;
FILE *num;
double rotation_h, rotation_v;
#line 8475 "ISIS_CRISP.c"
#undef curvature_v
#undef curvature
#undef segno
#undef my
#undef mx
#undef W
#undef alphay
#undef alphax
#undef Qcy
#undef Qcx
#undef R0
#undef louth
#undef linh
#undef loutw
#undef linw
#undef l
#undef h1
#undef w1
#undef option
#undef rotation_v
#undef rotation_h
#undef num
#undef ep
#undef file_name
#undef pos
#undef fu
#undef seg
#undef ii
#undef i
#undef Div1
#undef test
#undef p2_para
#undef div1
#undef u2
#undef u1
#undef mxi
#undef lbh
#undef lbw
#undef b_ell_q
#undef a_ell_q
#undef w2
#undef w12
#undef h2
#undef h12
#undef seg
#undef l_seg
#undef w2_out
#undef w1_in
#undef h2_out
#undef h1_in
#undef lhhalf
#undef lwhalf
#undef hhalf
#undef whalf
#undef hh
#undef ww
#undef w2c
#undef w1c
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'slit2' [12]. */
#define mccompcurname  slit2
#define mccompcurtype  Slit
#define mccompcurindex 12
#define xmin mccslit2_xmin
#define xmax mccslit2_xmax
#define ymin mccslit2_ymin
#define ymax mccslit2_ymax
#define radius mccslit2_radius
#define xwidth mccslit2_xwidth
#define yheight mccslit2_yheight
#undef yheight
#undef xwidth
#undef radius
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'lmon2' [13]. */
#define mccompcurname  lmon2
#define mccompcurtype  L_monitor
#define mccompcurindex 13
#define nL mcclmon2_nL
#define L_N mcclmon2_L_N
#define L_p mcclmon2_L_p
#define L_p2 mcclmon2_L_p2
#define filename mcclmon2_filename
#define xmin mcclmon2_xmin
#define xmax mcclmon2_xmax
#define ymin mcclmon2_ymin
#define ymax mcclmon2_ymax
#define xwidth mcclmon2_xwidth
#define yheight mcclmon2_yheight
#define Lmin mcclmon2_Lmin
#define Lmax mcclmon2_Lmax
#define restore_neutron mcclmon2_restore_neutron
#define nowritefile mcclmon2_nowritefile
#line 57 "/usr/share/mcstas/2.6rc1/monitors/L_monitor.comp"
double L_N[nL];
double L_p[nL], L_p2[nL];
#line 8582 "ISIS_CRISP.c"
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

/* User declarations for component 'samp1' [14]. */
#define mccompcurname  samp1
#define mccompcurtype  Multilayer_Sample
#define mccompcurindex 14
#define xwidth mccsamp1_xwidth
#define zlength mccsamp1_zlength
#define nlayer mccsamp1_nlayer
#define sldPar mccsamp1_sldPar
#define dPar mccsamp1_dPar
#define sigmaPar mccsamp1_sigmaPar
#define frac_inc mccsamp1_frac_inc
#define ythick mccsamp1_ythick
#define mu_inc mccsamp1_mu_inc
#define target_index mccsamp1_target_index
#define focus_xw mccsamp1_focus_xw
#define focus_yh mccsamp1_focus_yh
#define sldParPtr mccsamp1_sldParPtr
#define dParPtr mccsamp1_dParPtr
#define sigmaParPtr mccsamp1_sigmaParPtr
#define tx mccsamp1_tx
#define ty mccsamp1_ty
#define tz mccsamp1_tz
#line 70 "/usr/share/mcstas/2.6rc1/contrib/Multilayer_Sample.comp"
double sldParPtr[]   = sldPar;
double dParPtr[] = dPar;
double sigmaParPtr[] = sigmaPar;
double xmin, xmax, zmin, zmax;
double tx, ty, tz;
#line 8630 "ISIS_CRISP.c"
#undef tz
#undef ty
#undef tx
#undef sigmaParPtr
#undef dParPtr
#undef sldParPtr
#undef focus_yh
#undef focus_xw
#undef target_index
#undef mu_inc
#undef ythick
#undef frac_inc
#undef sigmaPar
#undef dPar
#undef sldPar
#undef nlayer
#undef zlength
#undef xwidth
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'slit3' [15]. */
#define mccompcurname  slit3
#define mccompcurtype  Slit
#define mccompcurindex 15
#define xmin mccslit3_xmin
#define xmax mccslit3_xmax
#define ymin mccslit3_ymin
#define ymax mccslit3_ymax
#define radius mccslit3_radius
#define xwidth mccslit3_xwidth
#define yheight mccslit3_yheight
#undef yheight
#undef xwidth
#undef radius
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'eguide2' [16]. */
#define mccompcurname  eguide2
#define mccompcurtype  Guide_tapering
#define mccompcurindex 16
#define w1c mcceguide2_w1c
#define w2c mcceguide2_w2c
#define ww mcceguide2_ww
#define hh mcceguide2_hh
#define whalf mcceguide2_whalf
#define hhalf mcceguide2_hhalf
#define lwhalf mcceguide2_lwhalf
#define lhhalf mcceguide2_lhhalf
#define h1_in mcceguide2_h1_in
#define h2_out mcceguide2_h2_out
#define w1_in mcceguide2_w1_in
#define w2_out mcceguide2_w2_out
#define l_seg mcceguide2_l_seg
#define seg mcceguide2_seg
#define h12 mcceguide2_h12
#define h2 mcceguide2_h2
#define w12 mcceguide2_w12
#define w2 mcceguide2_w2
#define a_ell_q mcceguide2_a_ell_q
#define b_ell_q mcceguide2_b_ell_q
#define lbw mcceguide2_lbw
#define lbh mcceguide2_lbh
#define mxi mcceguide2_mxi
#define u1 mcceguide2_u1
#define u2 mcceguide2_u2
#define div1 mcceguide2_div1
#define p2_para mcceguide2_p2_para
#define test mcceguide2_test
#define Div1 mcceguide2_Div1
#define i mcceguide2_i
#define ii mcceguide2_ii
#define seg mcceguide2_seg
#define fu mcceguide2_fu
#define pos mcceguide2_pos
#define file_name mcceguide2_file_name
#define ep mcceguide2_ep
#define num mcceguide2_num
#define rotation_h mcceguide2_rotation_h
#define rotation_v mcceguide2_rotation_v
#define option mcceguide2_option
#define w1 mcceguide2_w1
#define h1 mcceguide2_h1
#define l mcceguide2_l
#define linw mcceguide2_linw
#define loutw mcceguide2_loutw
#define linh mcceguide2_linh
#define louth mcceguide2_louth
#define R0 mcceguide2_R0
#define Qcx mcceguide2_Qcx
#define Qcy mcceguide2_Qcy
#define alphax mcceguide2_alphax
#define alphay mcceguide2_alphay
#define W mcceguide2_W
#define mx mcceguide2_mx
#define my mcceguide2_my
#define segno mcceguide2_segno
#define curvature mcceguide2_curvature
#define curvature_v mcceguide2_curvature_v
#line 96 "/usr/share/mcstas/2.6rc1/optics/Guide_tapering.comp"
double *w1c;
double *w2c;
double *ww, *hh;
double *whalf, *hhalf;
double *lwhalf, *lhhalf;
double *h1_in, *h2_out, *w1_in, *w2_out;
double l_seg, h12, h2, w12, w2, a_ell_q, b_ell_q, lbw, lbh;
double mxi ,u1 ,u2 ,div1, p2_para;
double test,Div1;
int i,ii,seg;
char *fu;
char *pos;
char file_name[1024];
char *ep;
FILE *num;
double rotation_h, rotation_v;
#line 8754 "ISIS_CRISP.c"
#undef curvature_v
#undef curvature
#undef segno
#undef my
#undef mx
#undef W
#undef alphay
#undef alphax
#undef Qcy
#undef Qcx
#undef R0
#undef louth
#undef linh
#undef loutw
#undef linw
#undef l
#undef h1
#undef w1
#undef option
#undef rotation_v
#undef rotation_h
#undef num
#undef ep
#undef file_name
#undef pos
#undef fu
#undef seg
#undef ii
#undef i
#undef Div1
#undef test
#undef p2_para
#undef div1
#undef u2
#undef u1
#undef mxi
#undef lbh
#undef lbw
#undef b_ell_q
#undef a_ell_q
#undef w2
#undef w12
#undef h2
#undef h12
#undef seg
#undef l_seg
#undef w2_out
#undef w1_in
#undef h2_out
#undef h1_in
#undef lhhalf
#undef lwhalf
#undef hhalf
#undef whalf
#undef hh
#undef ww
#undef w2c
#undef w1c
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'lmon3' [17]. */
#define mccompcurname  lmon3
#define mccompcurtype  L_monitor
#define mccompcurindex 17
#define nL mcclmon3_nL
#define L_N mcclmon3_L_N
#define L_p mcclmon3_L_p
#define L_p2 mcclmon3_L_p2
#define filename mcclmon3_filename
#define xmin mcclmon3_xmin
#define xmax mcclmon3_xmax
#define ymin mcclmon3_ymin
#define ymax mcclmon3_ymax
#define xwidth mcclmon3_xwidth
#define yheight mcclmon3_yheight
#define Lmin mcclmon3_Lmin
#define Lmax mcclmon3_Lmax
#define restore_neutron mcclmon3_restore_neutron
#define nowritefile mcclmon3_nowritefile
#line 57 "/usr/share/mcstas/2.6rc1/monitors/L_monitor.comp"
double L_N[nL];
double L_p[nL], L_p2[nL];
#line 8839 "ISIS_CRISP.c"
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

/* User declarations for component 'PSDdet1' [18]. */
#define mccompcurname  PSDdet1
#define mccompcurtype  PSD_monitor
#define mccompcurindex 18
#define PSD_N mccPSDdet1_PSD_N
#define PSD_p mccPSDdet1_PSD_p
#define PSD_p2 mccPSDdet1_PSD_p2
#define nx mccPSDdet1_nx
#define ny mccPSDdet1_ny
#define filename mccPSDdet1_filename
#define xmin mccPSDdet1_xmin
#define xmax mccPSDdet1_xmax
#define ymin mccPSDdet1_ymin
#define ymax mccPSDdet1_ymax
#define xwidth mccPSDdet1_xwidth
#define yheight mccPSDdet1_yheight
#define restore_neutron mccPSDdet1_restore_neutron
#line 62 "/usr/share/mcstas/2.6rc1/monitors/PSD_monitor.comp"
  DArray2d PSD_N;
  DArray2d PSD_p;
  DArray2d PSD_p2;
#line 8880 "ISIS_CRISP.c"
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

Coords mcposaa1, mcposra1;
Rotation mcrotaa1, mcrotra1;
Coords mcposaisis_source, mcposrisis_source;
Rotation mcrotaisis_source, mcrotrisis_source;
Coords mcposadefslit1, mcposrdefslit1;
Rotation mcrotadefslit1, mcrotrdefslit1;
Coords mcposadefslit2, mcposrdefslit2;
Rotation mcrotadefslit2, mcrotrdefslit2;
Coords mcposacoarseslit1, mcposrcoarseslit1;
Rotation mcrotacoarseslit1, mcrotrcoarseslit1;
Coords mcposacoarseslit2, mcposrcoarseslit2;
Rotation mcrotacoarseslit2, mcrotrcoarseslit2;
Coords mcposalmon1, mcposrlmon1;
Rotation mcrotalmon1, mcrotrlmon1;
Coords mcposaPSDmon1, mcposrPSDmon1;
Rotation mcrotaPSDmon1, mcrotrPSDmon1;
Coords mcposaslit1, mcposrslit1;
Rotation mcrotaslit1, mcrotrslit1;
Coords mcposaPSDmons1, mcposrPSDmons1;
Rotation mcrotaPSDmons1, mcrotrPSDmons1;
Coords mcposaeguide1, mcposreguide1;
Rotation mcrotaeguide1, mcrotreguide1;
Coords mcposaslit2, mcposrslit2;
Rotation mcrotaslit2, mcrotrslit2;
Coords mcposalmon2, mcposrlmon2;
Rotation mcrotalmon2, mcrotrlmon2;
Coords mcposasamp1, mcposrsamp1;
Rotation mcrotasamp1, mcrotrsamp1;
Coords mcposaslit3, mcposrslit3;
Rotation mcrotaslit3, mcrotrslit3;
Coords mcposaeguide2, mcposreguide2;
Rotation mcrotaeguide2, mcrotreguide2;
Coords mcposalmon3, mcposrlmon3;
Rotation mcrotalmon3, mcrotrlmon3;
Coords mcposaPSDdet1, mcposrPSDdet1;
Rotation mcrotaPSDdet1, mcrotrPSDdet1;

MCNUM mcnx, mcny, mcnz, mcnvx, mcnvy, mcnvz, mcnt, mcnsx, mcnsy, mcnsz, mcnp;

/* end declare */

void mcinit(void) {
#define mccompcurname  ISIS_CRISP
#define mccompcurtype  INSTRUMENT
#define mccompcurindex 0
#define mcposaISIS_CRISP coords_set(0,0,0)
#define glen mcipglen
#define flen mcipflen
#define w1 mcipw1
#define vm mcipvm
#define FRAC mcipFRAC
#line 50 "ISIS_CRISP.instr"
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
#line 8972 "ISIS_CRISP.c"
#undef FRAC
#undef vm
#undef w1
#undef flen
#undef glen
#undef mcposaISIS_CRISP
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
    /* Component a1. */
  /* Setting parameters for component a1. */
  SIG_MESSAGE("a1 (Init:SetPar)");

  SIG_MESSAGE("a1 (Init:Place/Rotate)");
  rot_set_rotation(mcrotaa1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 9001 "ISIS_CRISP.c"
  rot_copy(mcrotra1, mcrotaa1);
  mcposaa1 = coords_set(
#line 76 "ISIS_CRISP.instr"
    0,
#line 76 "ISIS_CRISP.instr"
    0,
#line 76 "ISIS_CRISP.instr"
    0);
#line 9010 "ISIS_CRISP.c"
  mctc1 = coords_neg(mcposaa1);
  mcposra1 = rot_apply(mcrotaa1, mctc1);
  mcDEBUG_COMPONENT("a1", mcposaa1, mcrotaa1)
  mccomp_posa[1] = mcposaa1;
  mccomp_posr[1] = mcposra1;
  mcNCounter[1]  = mcPCounter[1] = mcP2Counter[1] = 0;
  mcAbsorbProp[1]= 0;
    /* Component isis_source. */
  /* Setting parameters for component isis_source. */
  SIG_MESSAGE("isis_source (Init:SetPar)");
#line 79 "ISIS_CRISP.instr"
  mccisis_source_Emin = -6.5;
#line 79 "ISIS_CRISP.instr"
  mccisis_source_Emax = -0.55;
#line 79 "ISIS_CRISP.instr"
  mccisis_source_dist = 7.2695;
#line 80 "ISIS_CRISP.instr"
  mccisis_source_focus_xw = 0.045;
#line 80 "ISIS_CRISP.instr"
  mccisis_source_focus_yh = 0.0045;
#line 80 "ISIS_CRISP.instr"
  mccisis_source_xwidth = -1;
#line 80 "ISIS_CRISP.instr"
  mccisis_source_yheight = -1;
#line 81 "ISIS_CRISP.instr"
  mccisis_source_CAngle = 0.0;
#line 81 "ISIS_CRISP.instr"
  mccisis_source_SAC = 1;
#line 58 "ISIS_CRISP.instr"
  mccisis_source_Lmin = 0;
#line 58 "ISIS_CRISP.instr"
  mccisis_source_Lmax = 0;
#line 58 "ISIS_CRISP.instr"
  mccisis_source_target_index = + 1;
#line 58 "ISIS_CRISP.instr"
  mccisis_source_verbose = 0;
#line 9047 "ISIS_CRISP.c"

  SIG_MESSAGE("isis_source (Init:Place/Rotate)");
  rot_set_rotation(mcrotaisis_source,
#line 83 "ISIS_CRISP.instr"
    (1.5)*DEG2RAD,
#line 83 "ISIS_CRISP.instr"
    (0.0)*DEG2RAD,
#line 83 "ISIS_CRISP.instr"
    (0.0)*DEG2RAD);
#line 9057 "ISIS_CRISP.c"
  rot_transpose(mcrotaa1, mctr1);
  rot_mul(mcrotaisis_source, mctr1, mcrotrisis_source);
  mctc1 = coords_set(
#line 82 "ISIS_CRISP.instr"
    0.0,
#line 82 "ISIS_CRISP.instr"
    0.0,
#line 82 "ISIS_CRISP.instr"
    0.00001);
#line 9067 "ISIS_CRISP.c"
  rot_transpose(mcrotaa1, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaisis_source = coords_add(mcposaa1, mctc2);
  mctc1 = coords_sub(mcposaa1, mcposaisis_source);
  mcposrisis_source = rot_apply(mcrotaisis_source, mctc1);
  mcDEBUG_COMPONENT("isis_source", mcposaisis_source, mcrotaisis_source)
  mccomp_posa[2] = mcposaisis_source;
  mccomp_posr[2] = mcposrisis_source;
  mcNCounter[2]  = mcPCounter[2] = mcP2Counter[2] = 0;
  mcAbsorbProp[2]= 0;
    /* Component defslit1. */
  /* Setting parameters for component defslit1. */
  SIG_MESSAGE("defslit1 (Init:SetPar)");
#line 86 "ISIS_CRISP.instr"
  mccdefslit1_xmin = -0.0345;
#line 86 "ISIS_CRISP.instr"
  mccdefslit1_xmax = 0.0345;
#line 86 "ISIS_CRISP.instr"
  mccdefslit1_ymin = -0.01731;
#line 87 "ISIS_CRISP.instr"
  mccdefslit1_ymax = 0.01731;
#line 46 "ISIS_CRISP.instr"
  mccdefslit1_radius = 0;
#line 46 "ISIS_CRISP.instr"
  mccdefslit1_xwidth = 0;
#line 46 "ISIS_CRISP.instr"
  mccdefslit1_yheight = 0;
#line 9095 "ISIS_CRISP.c"

  SIG_MESSAGE("defslit1 (Init:Place/Rotate)");
  rot_set_rotation(mcrotadefslit1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 9102 "ISIS_CRISP.c"
  rot_transpose(mcrotaisis_source, mctr1);
  rot_mul(mcrotadefslit1, mctr1, mcrotrdefslit1);
  mcposadefslit1 = coords_set(
#line 88 "ISIS_CRISP.instr"
    0.0,
#line 88 "ISIS_CRISP.instr"
    -0.09829,
#line 88 "ISIS_CRISP.instr"
    3.7537);
#line 9112 "ISIS_CRISP.c"
  mctc1 = coords_sub(mcposaisis_source, mcposadefslit1);
  mcposrdefslit1 = rot_apply(mcrotadefslit1, mctc1);
  mcDEBUG_COMPONENT("defslit1", mcposadefslit1, mcrotadefslit1)
  mccomp_posa[3] = mcposadefslit1;
  mccomp_posr[3] = mcposrdefslit1;
  mcNCounter[3]  = mcPCounter[3] = mcP2Counter[3] = 0;
  mcAbsorbProp[3]= 0;
    /* Component defslit2. */
  /* Setting parameters for component defslit2. */
  SIG_MESSAGE("defslit2 (Init:SetPar)");
#line 91 "ISIS_CRISP.instr"
  mccdefslit2_xmin = -0.0289;
#line 91 "ISIS_CRISP.instr"
  mccdefslit2_xmax = 0.0289;
#line 91 "ISIS_CRISP.instr"
  mccdefslit2_ymin = -0.0123;
#line 92 "ISIS_CRISP.instr"
  mccdefslit2_ymax = 0.0123;
#line 46 "ISIS_CRISP.instr"
  mccdefslit2_radius = 0;
#line 46 "ISIS_CRISP.instr"
  mccdefslit2_xwidth = 0;
#line 46 "ISIS_CRISP.instr"
  mccdefslit2_yheight = 0;
#line 9137 "ISIS_CRISP.c"

  SIG_MESSAGE("defslit2 (Init:Place/Rotate)");
  rot_set_rotation(mcrotadefslit2,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 9144 "ISIS_CRISP.c"
  rot_transpose(mcrotadefslit1, mctr1);
  rot_mul(mcrotadefslit2, mctr1, mcrotrdefslit2);
  mcposadefslit2 = coords_set(
#line 93 "ISIS_CRISP.instr"
    0.0,
#line 93 "ISIS_CRISP.instr"
    -0.15976,
#line 93 "ISIS_CRISP.instr"
    6.101);
#line 9154 "ISIS_CRISP.c"
  mctc1 = coords_sub(mcposadefslit1, mcposadefslit2);
  mcposrdefslit2 = rot_apply(mcrotadefslit2, mctc1);
  mcDEBUG_COMPONENT("defslit2", mcposadefslit2, mcrotadefslit2)
  mccomp_posa[4] = mcposadefslit2;
  mccomp_posr[4] = mcposrdefslit2;
  mcNCounter[4]  = mcPCounter[4] = mcP2Counter[4] = 0;
  mcAbsorbProp[4]= 0;
    /* Component coarseslit1. */
  /* Setting parameters for component coarseslit1. */
  SIG_MESSAGE("coarseslit1 (Init:SetPar)");
#line 96 "ISIS_CRISP.instr"
  mcccoarseslit1_xmin = -0.03;
#line 96 "ISIS_CRISP.instr"
  mcccoarseslit1_xmax = 0.03;
#line 96 "ISIS_CRISP.instr"
  mcccoarseslit1_ymin = -0.1;
#line 97 "ISIS_CRISP.instr"
  mcccoarseslit1_ymax = 0.1;
#line 46 "ISIS_CRISP.instr"
  mcccoarseslit1_radius = 0;
#line 46 "ISIS_CRISP.instr"
  mcccoarseslit1_xwidth = 0;
#line 46 "ISIS_CRISP.instr"
  mcccoarseslit1_yheight = 0;
#line 9179 "ISIS_CRISP.c"

  SIG_MESSAGE("coarseslit1 (Init:Place/Rotate)");
  rot_set_rotation(mcrotacoarseslit1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 9186 "ISIS_CRISP.c"
  rot_transpose(mcrotadefslit2, mctr1);
  rot_mul(mcrotacoarseslit1, mctr1, mcrotrcoarseslit1);
  mcposacoarseslit1 = coords_set(
#line 98 "ISIS_CRISP.instr"
    0.0,
#line 98 "ISIS_CRISP.instr"
    -0.1757,
#line 98 "ISIS_CRISP.instr"
    6.7096);
#line 9196 "ISIS_CRISP.c"
  mctc1 = coords_sub(mcposadefslit2, mcposacoarseslit1);
  mcposrcoarseslit1 = rot_apply(mcrotacoarseslit1, mctc1);
  mcDEBUG_COMPONENT("coarseslit1", mcposacoarseslit1, mcrotacoarseslit1)
  mccomp_posa[5] = mcposacoarseslit1;
  mccomp_posr[5] = mcposrcoarseslit1;
  mcNCounter[5]  = mcPCounter[5] = mcP2Counter[5] = 0;
  mcAbsorbProp[5]= 0;
    /* Component coarseslit2. */
  /* Setting parameters for component coarseslit2. */
  SIG_MESSAGE("coarseslit2 (Init:SetPar)");
#line 101 "ISIS_CRISP.instr"
  mcccoarseslit2_xmin = -0.1;
#line 101 "ISIS_CRISP.instr"
  mcccoarseslit2_xmax = 0.1;
#line 101 "ISIS_CRISP.instr"
  mcccoarseslit2_ymin = -0.01;
#line 102 "ISIS_CRISP.instr"
  mcccoarseslit2_ymax = 0.01;
#line 46 "ISIS_CRISP.instr"
  mcccoarseslit2_radius = 0;
#line 46 "ISIS_CRISP.instr"
  mcccoarseslit2_xwidth = 0;
#line 46 "ISIS_CRISP.instr"
  mcccoarseslit2_yheight = 0;
#line 9221 "ISIS_CRISP.c"

  SIG_MESSAGE("coarseslit2 (Init:Place/Rotate)");
  rot_set_rotation(mcrotacoarseslit2,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 9228 "ISIS_CRISP.c"
  rot_transpose(mcrotacoarseslit1, mctr1);
  rot_mul(mcrotacoarseslit2, mctr1, mcrotrcoarseslit2);
  mcposacoarseslit2 = coords_set(
#line 103 "ISIS_CRISP.instr"
    0.0,
#line 103 "ISIS_CRISP.instr"
    -0.17832,
#line 103 "ISIS_CRISP.instr"
    6.8096);
#line 9238 "ISIS_CRISP.c"
  mctc1 = coords_sub(mcposacoarseslit1, mcposacoarseslit2);
  mcposrcoarseslit2 = rot_apply(mcrotacoarseslit2, mctc1);
  mcDEBUG_COMPONENT("coarseslit2", mcposacoarseslit2, mcrotacoarseslit2)
  mccomp_posa[6] = mcposacoarseslit2;
  mccomp_posr[6] = mcposrcoarseslit2;
  mcNCounter[6]  = mcPCounter[6] = mcP2Counter[6] = 0;
  mcAbsorbProp[6]= 0;
    /* Component lmon1. */
  /* Setting parameters for component lmon1. */
  SIG_MESSAGE("lmon1 (Init:SetPar)");
#line 106 "ISIS_CRISP.instr"
  if("lmon1.dat") strncpy(mcclmon1_filename, "lmon1.dat" ? "lmon1.dat" : "", 16384); else mcclmon1_filename[0]='\0';
#line 106 "ISIS_CRISP.instr"
  mcclmon1_xmin = -0.06;
#line 107 "ISIS_CRISP.instr"
  mcclmon1_xmax = 0.06;
#line 107 "ISIS_CRISP.instr"
  mcclmon1_ymin = -0.04;
#line 107 "ISIS_CRISP.instr"
  mcclmon1_ymax = 0.04;
#line 51 "ISIS_CRISP.instr"
  mcclmon1_xwidth = 0;
#line 51 "ISIS_CRISP.instr"
  mcclmon1_yheight = 0;
#line 107 "ISIS_CRISP.instr"
  mcclmon1_Lmin = 0.0;
#line 108 "ISIS_CRISP.instr"
  mcclmon1_Lmax = 10.0;
#line 51 "ISIS_CRISP.instr"
  mcclmon1_restore_neutron = 0;
#line 51 "ISIS_CRISP.instr"
  mcclmon1_nowritefile = 0;
#line 9271 "ISIS_CRISP.c"

  SIG_MESSAGE("lmon1 (Init:Place/Rotate)");
  rot_set_rotation(mcrotalmon1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 9278 "ISIS_CRISP.c"
  rot_transpose(mcrotacoarseslit2, mctr1);
  rot_mul(mcrotalmon1, mctr1, mcrotrlmon1);
  mcposalmon1 = coords_set(
#line 109 "ISIS_CRISP.instr"
    0.0,
#line 109 "ISIS_CRISP.instr"
    -0.1833,
#line 109 "ISIS_CRISP.instr"
    7.0);
#line 9288 "ISIS_CRISP.c"
  mctc1 = coords_sub(mcposacoarseslit2, mcposalmon1);
  mcposrlmon1 = rot_apply(mcrotalmon1, mctc1);
  mcDEBUG_COMPONENT("lmon1", mcposalmon1, mcrotalmon1)
  mccomp_posa[7] = mcposalmon1;
  mccomp_posr[7] = mcposrlmon1;
  mcNCounter[7]  = mcPCounter[7] = mcP2Counter[7] = 0;
  mcAbsorbProp[7]= 0;
    /* Component PSDmon1. */
  /* Setting parameters for component PSDmon1. */
  SIG_MESSAGE("PSDmon1 (Init:SetPar)");
#line 112 "ISIS_CRISP.instr"
  mccPSDmon1_nx = 100;
#line 112 "ISIS_CRISP.instr"
  mccPSDmon1_ny = 100;
#line 112 "ISIS_CRISP.instr"
  if("PSD1.dat") strncpy(mccPSDmon1_filename, "PSD1.dat" ? "PSD1.dat" : "", 16384); else mccPSDmon1_filename[0]='\0';
#line 113 "ISIS_CRISP.instr"
  mccPSDmon1_xmin = -0.05;
#line 113 "ISIS_CRISP.instr"
  mccPSDmon1_xmax = 0.05;
#line 113 "ISIS_CRISP.instr"
  mccPSDmon1_ymin = -0.01;
#line 113 "ISIS_CRISP.instr"
  mccPSDmon1_ymax = 0.01;
#line 50 "ISIS_CRISP.instr"
  mccPSDmon1_xwidth = 0;
#line 50 "ISIS_CRISP.instr"
  mccPSDmon1_yheight = 0;
#line 51 "ISIS_CRISP.instr"
  mccPSDmon1_restore_neutron = 0;
#line 9319 "ISIS_CRISP.c"

  SIG_MESSAGE("PSDmon1 (Init:Place/Rotate)");
  rot_set_rotation(mcrotaPSDmon1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 9326 "ISIS_CRISP.c"
  rot_transpose(mcrotalmon1, mctr1);
  rot_mul(mcrotaPSDmon1, mctr1, mcrotrPSDmon1);
  mcposaPSDmon1 = coords_set(
#line 114 "ISIS_CRISP.instr"
    0.0,
#line 114 "ISIS_CRISP.instr"
    -0.1833,
#line 114 "ISIS_CRISP.instr"
    7.01);
#line 9336 "ISIS_CRISP.c"
  mctc1 = coords_sub(mcposalmon1, mcposaPSDmon1);
  mcposrPSDmon1 = rot_apply(mcrotaPSDmon1, mctc1);
  mcDEBUG_COMPONENT("PSDmon1", mcposaPSDmon1, mcrotaPSDmon1)
  mccomp_posa[8] = mcposaPSDmon1;
  mccomp_posr[8] = mcposrPSDmon1;
  mcNCounter[8]  = mcPCounter[8] = mcP2Counter[8] = 0;
  mcAbsorbProp[8]= 0;
    /* Component slit1. */
  /* Setting parameters for component slit1. */
  SIG_MESSAGE("slit1 (Init:SetPar)");
#line 46 "ISIS_CRISP.instr"
  mccslit1_xmin = 0;
#line 46 "ISIS_CRISP.instr"
  mccslit1_xmax = 0;
#line 46 "ISIS_CRISP.instr"
  mccslit1_ymin = 0;
#line 46 "ISIS_CRISP.instr"
  mccslit1_ymax = 0;
#line 46 "ISIS_CRISP.instr"
  mccslit1_radius = 0;
#line 117 "ISIS_CRISP.instr"
  mccslit1_xwidth = 0.04;
#line 117 "ISIS_CRISP.instr"
  mccslit1_yheight = 0.004;
#line 9361 "ISIS_CRISP.c"

  SIG_MESSAGE("slit1 (Init:Place/Rotate)");
  rot_set_rotation(mcrotaslit1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 9368 "ISIS_CRISP.c"
  rot_transpose(mcrotaPSDmon1, mctr1);
  rot_mul(mcrotaslit1, mctr1, mcrotrslit1);
  mcposaslit1 = coords_set(
#line 118 "ISIS_CRISP.instr"
    0.0,
#line 118 "ISIS_CRISP.instr"
    -0.1904,
#line 118 "ISIS_CRISP.instr"
    7.2695);
#line 9378 "ISIS_CRISP.c"
  mctc1 = coords_sub(mcposaPSDmon1, mcposaslit1);
  mcposrslit1 = rot_apply(mcrotaslit1, mctc1);
  mcDEBUG_COMPONENT("slit1", mcposaslit1, mcrotaslit1)
  mccomp_posa[9] = mcposaslit1;
  mccomp_posr[9] = mcposrslit1;
  mcNCounter[9]  = mcPCounter[9] = mcP2Counter[9] = 0;
  mcAbsorbProp[9]= 0;
    /* Component PSDmons1. */
  /* Setting parameters for component PSDmons1. */
  SIG_MESSAGE("PSDmons1 (Init:SetPar)");
#line 121 "ISIS_CRISP.instr"
  mccPSDmons1_nx = 100;
#line 121 "ISIS_CRISP.instr"
  mccPSDmons1_ny = 100;
#line 121 "ISIS_CRISP.instr"
  if("PSDs1.dat") strncpy(mccPSDmons1_filename, "PSDs1.dat" ? "PSDs1.dat" : "", 16384); else mccPSDmons1_filename[0]='\0';
#line 122 "ISIS_CRISP.instr"
  mccPSDmons1_xmin = -0.05;
#line 122 "ISIS_CRISP.instr"
  mccPSDmons1_xmax = 0.05;
#line 122 "ISIS_CRISP.instr"
  mccPSDmons1_ymin = -0.01;
#line 122 "ISIS_CRISP.instr"
  mccPSDmons1_ymax = 0.01;
#line 50 "ISIS_CRISP.instr"
  mccPSDmons1_xwidth = 0;
#line 50 "ISIS_CRISP.instr"
  mccPSDmons1_yheight = 0;
#line 51 "ISIS_CRISP.instr"
  mccPSDmons1_restore_neutron = 0;
#line 9409 "ISIS_CRISP.c"

  SIG_MESSAGE("PSDmons1 (Init:Place/Rotate)");
  rot_set_rotation(mcrotaPSDmons1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 9416 "ISIS_CRISP.c"
  rot_transpose(mcrotaslit1, mctr1);
  rot_mul(mcrotaPSDmons1, mctr1, mcrotrPSDmons1);
  mcposaPSDmons1 = coords_set(
#line 123 "ISIS_CRISP.instr"
    0.0,
#line 123 "ISIS_CRISP.instr"
    7.27 * t15,
#line 123 "ISIS_CRISP.instr"
    7.27);
#line 9426 "ISIS_CRISP.c"
  mctc1 = coords_sub(mcposaslit1, mcposaPSDmons1);
  mcposrPSDmons1 = rot_apply(mcrotaPSDmons1, mctc1);
  mcDEBUG_COMPONENT("PSDmons1", mcposaPSDmons1, mcrotaPSDmons1)
  mccomp_posa[10] = mcposaPSDmons1;
  mccomp_posr[10] = mcposrPSDmons1;
  mcNCounter[10]  = mcPCounter[10] = mcP2Counter[10] = 0;
  mcAbsorbProp[10]= 0;
    /* Component eguide1. */
  /* Setting parameters for component eguide1. */
  SIG_MESSAGE("eguide1 (Init:SetPar)");
#line 126 "ISIS_CRISP.instr"
  if("elliptical") strncpy(mcceguide1_option, "elliptical" ? "elliptical" : "", 16384); else mcceguide1_option[0]='\0';
#line 126 "ISIS_CRISP.instr"
  mcceguide1_w1 = mcipw1;
#line 126 "ISIS_CRISP.instr"
  mcceguide1_h1 = 0.1;
#line 126 "ISIS_CRISP.instr"
  mcceguide1_l = mcipglen;
#line 127 "ISIS_CRISP.instr"
  mcceguide1_linw = 10.0;
#line 127 "ISIS_CRISP.instr"
  mcceguide1_loutw = mcipflen;
#line 81 "ISIS_CRISP.instr"
  mcceguide1_linh = 0;
#line 81 "ISIS_CRISP.instr"
  mcceguide1_louth = 0;
#line 81 "ISIS_CRISP.instr"
  mcceguide1_R0 = 0.99;
#line 82 "ISIS_CRISP.instr"
  mcceguide1_Qcx = 0.021;
#line 82 "ISIS_CRISP.instr"
  mcceguide1_Qcy = 0.021;
#line 82 "ISIS_CRISP.instr"
  mcceguide1_alphax = 6.07;
#line 82 "ISIS_CRISP.instr"
  mcceguide1_alphay = 6.07;
#line 82 "ISIS_CRISP.instr"
  mcceguide1_W = 0.003;
#line 127 "ISIS_CRISP.instr"
  mcceguide1_mx = mcipvm;
#line 127 "ISIS_CRISP.instr"
  mcceguide1_my = 0;
#line 126 "ISIS_CRISP.instr"
  mcceguide1_segno = 10;
#line 83 "ISIS_CRISP.instr"
  mcceguide1_curvature = 0;
#line 83 "ISIS_CRISP.instr"
  mcceguide1_curvature_v = 0;
#line 9475 "ISIS_CRISP.c"

  SIG_MESSAGE("eguide1 (Init:Place/Rotate)");
  rot_set_rotation(mcrotaeguide1,
#line 129 "ISIS_CRISP.instr"
    (1.5)*DEG2RAD,
#line 129 "ISIS_CRISP.instr"
    (0.0)*DEG2RAD,
#line 129 "ISIS_CRISP.instr"
    (0.0)*DEG2RAD);
#line 9485 "ISIS_CRISP.c"
  rot_transpose(mcrotaPSDmons1, mctr1);
  rot_mul(mcrotaeguide1, mctr1, mcrotreguide1);
  mcposaeguide1 = coords_set(
#line 128 "ISIS_CRISP.instr"
    0.0,
#line 128 "ISIS_CRISP.instr"
    ( spos - mcipflen - mcipglen ) * t15,
#line 128 "ISIS_CRISP.instr"
    spos - mcipflen - mcipglen);
#line 9495 "ISIS_CRISP.c"
  mctc1 = coords_sub(mcposaPSDmons1, mcposaeguide1);
  mcposreguide1 = rot_apply(mcrotaeguide1, mctc1);
  mcDEBUG_COMPONENT("eguide1", mcposaeguide1, mcrotaeguide1)
  mccomp_posa[11] = mcposaeguide1;
  mccomp_posr[11] = mcposreguide1;
  mcNCounter[11]  = mcPCounter[11] = mcP2Counter[11] = 0;
  mcAbsorbProp[11]= 0;
    /* Component slit2. */
  /* Setting parameters for component slit2. */
  SIG_MESSAGE("slit2 (Init:SetPar)");
#line 46 "ISIS_CRISP.instr"
  mccslit2_xmin = 0;
#line 46 "ISIS_CRISP.instr"
  mccslit2_xmax = 0;
#line 46 "ISIS_CRISP.instr"
  mccslit2_ymin = 0;
#line 46 "ISIS_CRISP.instr"
  mccslit2_ymax = 0;
#line 46 "ISIS_CRISP.instr"
  mccslit2_radius = 0;
#line 132 "ISIS_CRISP.instr"
  mccslit2_xwidth = tend;
#line 132 "ISIS_CRISP.instr"
  mccslit2_yheight = 0.0025;
#line 9520 "ISIS_CRISP.c"

  SIG_MESSAGE("slit2 (Init:Place/Rotate)");
  rot_set_rotation(mcrotaslit2,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 9527 "ISIS_CRISP.c"
  rot_transpose(mcrotaeguide1, mctr1);
  rot_mul(mcrotaslit2, mctr1, mcrotrslit2);
  mcposaslit2 = coords_set(
#line 133 "ISIS_CRISP.instr"
    0.0,
#line 133 "ISIS_CRISP.instr"
    -0.2581,
#line 133 "ISIS_CRISP.instr"
    9.8555);
#line 9537 "ISIS_CRISP.c"
  mctc1 = coords_sub(mcposaeguide1, mcposaslit2);
  mcposrslit2 = rot_apply(mcrotaslit2, mctc1);
  mcDEBUG_COMPONENT("slit2", mcposaslit2, mcrotaslit2)
  mccomp_posa[12] = mcposaslit2;
  mccomp_posr[12] = mcposrslit2;
  mcNCounter[12]  = mcPCounter[12] = mcP2Counter[12] = 0;
  mcAbsorbProp[12]= 0;
    /* Component lmon2. */
  /* Setting parameters for component lmon2. */
  SIG_MESSAGE("lmon2 (Init:SetPar)");
#line 136 "ISIS_CRISP.instr"
  if("lmon2.dat") strncpy(mcclmon2_filename, "lmon2.dat" ? "lmon2.dat" : "", 16384); else mcclmon2_filename[0]='\0';
#line 136 "ISIS_CRISP.instr"
  mcclmon2_xmin = -0.06;
#line 137 "ISIS_CRISP.instr"
  mcclmon2_xmax = 0.06;
#line 137 "ISIS_CRISP.instr"
  mcclmon2_ymin = -0.04;
#line 137 "ISIS_CRISP.instr"
  mcclmon2_ymax = 0.04;
#line 51 "ISIS_CRISP.instr"
  mcclmon2_xwidth = 0;
#line 51 "ISIS_CRISP.instr"
  mcclmon2_yheight = 0;
#line 137 "ISIS_CRISP.instr"
  mcclmon2_Lmin = 0.0;
#line 138 "ISIS_CRISP.instr"
  mcclmon2_Lmax = 10.0;
#line 51 "ISIS_CRISP.instr"
  mcclmon2_restore_neutron = 0;
#line 51 "ISIS_CRISP.instr"
  mcclmon2_nowritefile = 0;
#line 9570 "ISIS_CRISP.c"

  SIG_MESSAGE("lmon2 (Init:Place/Rotate)");
  rot_set_rotation(mcrotalmon2,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 9577 "ISIS_CRISP.c"
  rot_transpose(mcrotaslit2, mctr1);
  rot_mul(mcrotalmon2, mctr1, mcrotrlmon2);
  mcposalmon2 = coords_set(
#line 139 "ISIS_CRISP.instr"
    0.0,
#line 139 "ISIS_CRISP.instr"
    -0.2608,
#line 139 "ISIS_CRISP.instr"
    9.96);
#line 9587 "ISIS_CRISP.c"
  mctc1 = coords_sub(mcposaslit2, mcposalmon2);
  mcposrlmon2 = rot_apply(mcrotalmon2, mctc1);
  mcDEBUG_COMPONENT("lmon2", mcposalmon2, mcrotalmon2)
  mccomp_posa[13] = mcposalmon2;
  mccomp_posr[13] = mcposrlmon2;
  mcNCounter[13]  = mcPCounter[13] = mcP2Counter[13] = 0;
  mcAbsorbProp[13]= 0;
    /* Component samp1. */
  /* Setting parameters for component samp1. */
  SIG_MESSAGE("samp1 (Init:SetPar)");

  SIG_MESSAGE("samp1 (Init:Place/Rotate)");
  rot_set_rotation(mcrotasamp1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 9604 "ISIS_CRISP.c"
  rot_transpose(mcrotalmon2, mctr1);
  rot_mul(mcrotasamp1, mctr1, mcrotrsamp1);
  mcposasamp1 = coords_set(
#line 146 "ISIS_CRISP.instr"
    0.0,
#line 146 "ISIS_CRISP.instr"
    spos * t15,
#line 146 "ISIS_CRISP.instr"
    spos);
#line 9614 "ISIS_CRISP.c"
  mctc1 = coords_sub(mcposalmon2, mcposasamp1);
  mcposrsamp1 = rot_apply(mcrotasamp1, mctc1);
  mcDEBUG_COMPONENT("samp1", mcposasamp1, mcrotasamp1)
  mccomp_posa[14] = mcposasamp1;
  mccomp_posr[14] = mcposrsamp1;
  mcNCounter[14]  = mcPCounter[14] = mcP2Counter[14] = 0;
  mcAbsorbProp[14]= 0;
    /* Component slit3. */
  /* Setting parameters for component slit3. */
  SIG_MESSAGE("slit3 (Init:SetPar)");
#line 46 "ISIS_CRISP.instr"
  mccslit3_xmin = 0;
#line 46 "ISIS_CRISP.instr"
  mccslit3_xmax = 0;
#line 46 "ISIS_CRISP.instr"
  mccslit3_ymin = 0;
#line 46 "ISIS_CRISP.instr"
  mccslit3_ymax = 0;
#line 46 "ISIS_CRISP.instr"
  mccslit3_radius = 0;
#line 149 "ISIS_CRISP.instr"
  mccslit3_xwidth = tend;
#line 149 "ISIS_CRISP.instr"
  mccslit3_yheight = 0.003;
#line 9639 "ISIS_CRISP.c"

  SIG_MESSAGE("slit3 (Init:Place/Rotate)");
  rot_set_rotation(mcrotaslit3,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 9646 "ISIS_CRISP.c"
  rot_transpose(mcrotasamp1, mctr1);
  rot_mul(mcrotaslit3, mctr1, mcrotrslit3);
  mcposaslit3 = coords_set(
#line 150 "ISIS_CRISP.instr"
    0.0,
#line 150 "ISIS_CRISP.instr"
    ( spos + mcipflen -0.01 ) * t15 -2.0 * ( mcipflen -0.01 ) * t15,
#line 150 "ISIS_CRISP.instr"
    spos + mcipflen -0.01);
#line 9656 "ISIS_CRISP.c"
  mctc1 = coords_sub(mcposasamp1, mcposaslit3);
  mcposrslit3 = rot_apply(mcrotaslit3, mctc1);
  mcDEBUG_COMPONENT("slit3", mcposaslit3, mcrotaslit3)
  mccomp_posa[15] = mcposaslit3;
  mccomp_posr[15] = mcposrslit3;
  mcNCounter[15]  = mcPCounter[15] = mcP2Counter[15] = 0;
  mcAbsorbProp[15]= 0;
    /* Component eguide2. */
  /* Setting parameters for component eguide2. */
  SIG_MESSAGE("eguide2 (Init:SetPar)");
#line 153 "ISIS_CRISP.instr"
  if("elliptical") strncpy(mcceguide2_option, "elliptical" ? "elliptical" : "", 16384); else mcceguide2_option[0]='\0';
#line 153 "ISIS_CRISP.instr"
  mcceguide2_w1 = tend;
#line 153 "ISIS_CRISP.instr"
  mcceguide2_h1 = 0.1;
#line 153 "ISIS_CRISP.instr"
  mcceguide2_l = mcipglen;
#line 154 "ISIS_CRISP.instr"
  mcceguide2_linw = mcipflen;
#line 154 "ISIS_CRISP.instr"
  mcceguide2_loutw = 10.0;
#line 81 "ISIS_CRISP.instr"
  mcceguide2_linh = 0;
#line 81 "ISIS_CRISP.instr"
  mcceguide2_louth = 0;
#line 81 "ISIS_CRISP.instr"
  mcceguide2_R0 = 0.99;
#line 82 "ISIS_CRISP.instr"
  mcceguide2_Qcx = 0.021;
#line 82 "ISIS_CRISP.instr"
  mcceguide2_Qcy = 0.021;
#line 82 "ISIS_CRISP.instr"
  mcceguide2_alphax = 6.07;
#line 82 "ISIS_CRISP.instr"
  mcceguide2_alphay = 6.07;
#line 82 "ISIS_CRISP.instr"
  mcceguide2_W = 0.003;
#line 154 "ISIS_CRISP.instr"
  mcceguide2_mx = mcipvm;
#line 154 "ISIS_CRISP.instr"
  mcceguide2_my = 0;
#line 153 "ISIS_CRISP.instr"
  mcceguide2_segno = 10;
#line 83 "ISIS_CRISP.instr"
  mcceguide2_curvature = 0;
#line 83 "ISIS_CRISP.instr"
  mcceguide2_curvature_v = 0;
#line 9705 "ISIS_CRISP.c"

  SIG_MESSAGE("eguide2 (Init:Place/Rotate)");
  rot_set_rotation(mcrotaeguide2,
#line 156 "ISIS_CRISP.instr"
    (-1.5)*DEG2RAD,
#line 156 "ISIS_CRISP.instr"
    (0.0)*DEG2RAD,
#line 156 "ISIS_CRISP.instr"
    (0.0)*DEG2RAD);
#line 9715 "ISIS_CRISP.c"
  rot_transpose(mcrotaslit3, mctr1);
  rot_mul(mcrotaeguide2, mctr1, mcrotreguide2);
  mcposaeguide2 = coords_set(
#line 155 "ISIS_CRISP.instr"
    0.0,
#line 155 "ISIS_CRISP.instr"
    ( spos + mcipflen ) * t15 -2.0 * mcipflen * t15,
#line 155 "ISIS_CRISP.instr"
    spos + mcipflen);
#line 9725 "ISIS_CRISP.c"
  mctc1 = coords_sub(mcposaslit3, mcposaeguide2);
  mcposreguide2 = rot_apply(mcrotaeguide2, mctc1);
  mcDEBUG_COMPONENT("eguide2", mcposaeguide2, mcrotaeguide2)
  mccomp_posa[16] = mcposaeguide2;
  mccomp_posr[16] = mcposreguide2;
  mcNCounter[16]  = mcPCounter[16] = mcP2Counter[16] = 0;
  mcAbsorbProp[16]= 0;
    /* Component lmon3. */
  /* Setting parameters for component lmon3. */
  SIG_MESSAGE("lmon3 (Init:SetPar)");
#line 159 "ISIS_CRISP.instr"
  if("lmon3.dat") strncpy(mcclmon3_filename, "lmon3.dat" ? "lmon3.dat" : "", 16384); else mcclmon3_filename[0]='\0';
#line 159 "ISIS_CRISP.instr"
  mcclmon3_xmin = -0.05;
#line 160 "ISIS_CRISP.instr"
  mcclmon3_xmax = 0.05;
#line 160 "ISIS_CRISP.instr"
  mcclmon3_ymin = -0.05;
#line 160 "ISIS_CRISP.instr"
  mcclmon3_ymax = 0.05;
#line 51 "ISIS_CRISP.instr"
  mcclmon3_xwidth = 0;
#line 51 "ISIS_CRISP.instr"
  mcclmon3_yheight = 0;
#line 160 "ISIS_CRISP.instr"
  mcclmon3_Lmin = 0.0;
#line 161 "ISIS_CRISP.instr"
  mcclmon3_Lmax = 10.0;
#line 51 "ISIS_CRISP.instr"
  mcclmon3_restore_neutron = 0;
#line 51 "ISIS_CRISP.instr"
  mcclmon3_nowritefile = 0;
#line 9758 "ISIS_CRISP.c"

  SIG_MESSAGE("lmon3 (Init:Place/Rotate)");
  rot_set_rotation(mcrotalmon3,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 9765 "ISIS_CRISP.c"
  rot_transpose(mcrotaeguide2, mctr1);
  rot_mul(mcrotalmon3, mctr1, mcrotrlmon3);
  mcposalmon3 = coords_set(
#line 162 "ISIS_CRISP.instr"
    0.0,
#line 162 "ISIS_CRISP.instr"
    12.11 * t15 -2.0 * ( 12.11 - spos ) * t15,
#line 162 "ISIS_CRISP.instr"
    12.11);
#line 9775 "ISIS_CRISP.c"
  mctc1 = coords_sub(mcposaeguide2, mcposalmon3);
  mcposrlmon3 = rot_apply(mcrotalmon3, mctc1);
  mcDEBUG_COMPONENT("lmon3", mcposalmon3, mcrotalmon3)
  mccomp_posa[17] = mcposalmon3;
  mccomp_posr[17] = mcposrlmon3;
  mcNCounter[17]  = mcPCounter[17] = mcP2Counter[17] = 0;
  mcAbsorbProp[17]= 0;
    /* Component PSDdet1. */
  /* Setting parameters for component PSDdet1. */
  SIG_MESSAGE("PSDdet1 (Init:SetPar)");
#line 165 "ISIS_CRISP.instr"
  mccPSDdet1_nx = 100;
#line 165 "ISIS_CRISP.instr"
  mccPSDdet1_ny = 100;
#line 165 "ISIS_CRISP.instr"
  if("PSD3.dat") strncpy(mccPSDdet1_filename, "PSD3.dat" ? "PSD3.dat" : "", 16384); else mccPSDdet1_filename[0]='\0';
#line 166 "ISIS_CRISP.instr"
  mccPSDdet1_xmin = -0.05;
#line 166 "ISIS_CRISP.instr"
  mccPSDdet1_xmax = 0.05;
#line 166 "ISIS_CRISP.instr"
  mccPSDdet1_ymin = -0.05;
#line 166 "ISIS_CRISP.instr"
  mccPSDdet1_ymax = 0.05;
#line 50 "ISIS_CRISP.instr"
  mccPSDdet1_xwidth = 0;
#line 50 "ISIS_CRISP.instr"
  mccPSDdet1_yheight = 0;
#line 51 "ISIS_CRISP.instr"
  mccPSDdet1_restore_neutron = 0;
#line 9806 "ISIS_CRISP.c"

  SIG_MESSAGE("PSDdet1 (Init:Place/Rotate)");
  rot_set_rotation(mcrotaPSDdet1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 9813 "ISIS_CRISP.c"
  rot_transpose(mcrotalmon3, mctr1);
  rot_mul(mcrotaPSDdet1, mctr1, mcrotrPSDdet1);
  mcposaPSDdet1 = coords_set(
#line 167 "ISIS_CRISP.instr"
    0.0,
#line 167 "ISIS_CRISP.instr"
    12.111 * t15 -2.0 * ( 12.111 - spos ) * t15,
#line 167 "ISIS_CRISP.instr"
    12.111);
#line 9823 "ISIS_CRISP.c"
  mctc1 = coords_sub(mcposalmon3, mcposaPSDdet1);
  mcposrPSDdet1 = rot_apply(mcrotaPSDdet1, mctc1);
  mcDEBUG_COMPONENT("PSDdet1", mcposaPSDdet1, mcrotaPSDdet1)
  mccomp_posa[18] = mcposaPSDdet1;
  mccomp_posr[18] = mcposrPSDdet1;
  mcNCounter[18]  = mcPCounter[18] = mcP2Counter[18] = 0;
  mcAbsorbProp[18]= 0;
  /* Component initializations. */
  /* Initializations for component a1. */
  SIG_MESSAGE("a1 (Init)");

  /* Initializations for component isis_source. */
  SIG_MESSAGE("isis_source (Init)");
#define mccompcurname  isis_source
#define mccompcurtype  ISIS_moderator
#define mccompcurindex 2
#define Face mccisis_source_Face
#define p_in mccisis_source_p_in
#define Tnpts mccisis_source_Tnpts
#define scaleSize mccisis_source_scaleSize
#define angleArea mccisis_source_angleArea
#define Nsim mccisis_source_Nsim
#define Ncount mccisis_source_Ncount
#define TS mccisis_source_TS
#define rtE0 mccisis_source_rtE0
#define rtE1 mccisis_source_rtE1
#define rtmodX mccisis_source_rtmodX
#define rtmodY mccisis_source_rtmodY
#define TargetStation mccisis_source_TargetStation
#define CurrentWeight mccisis_source_CurrentWeight
#define Emin mccisis_source_Emin
#define Emax mccisis_source_Emax
#define dist mccisis_source_dist
#define focus_xw mccisis_source_focus_xw
#define focus_yh mccisis_source_focus_yh
#define xwidth mccisis_source_xwidth
#define yheight mccisis_source_yheight
#define CAngle mccisis_source_CAngle
#define SAC mccisis_source_SAC
#define Lmin mccisis_source_Lmin
#define Lmax mccisis_source_Lmax
#define target_index mccisis_source_target_index
#define verbose mccisis_source_verbose
#line 836 "/usr/share/mcstas/2.6rc1/contrib/ISIS_moderator.comp"
{
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
  
}
#line 10090 "ISIS_CRISP.c"
#undef verbose
#undef target_index
#undef Lmax
#undef Lmin
#undef SAC
#undef CAngle
#undef yheight
#undef xwidth
#undef focus_yh
#undef focus_xw
#undef dist
#undef Emax
#undef Emin
#undef CurrentWeight
#undef TargetStation
#undef rtmodY
#undef rtmodX
#undef rtE1
#undef rtE0
#undef TS
#undef Ncount
#undef Nsim
#undef angleArea
#undef scaleSize
#undef Tnpts
#undef p_in
#undef Face
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component defslit1. */
  SIG_MESSAGE("defslit1 (Init)");
#define mccompcurname  defslit1
#define mccompcurtype  Slit
#define mccompcurindex 3
#define xmin mccdefslit1_xmin
#define xmax mccdefslit1_xmax
#define ymin mccdefslit1_ymin
#define ymax mccdefslit1_ymax
#define radius mccdefslit1_radius
#define xwidth mccdefslit1_xwidth
#define yheight mccdefslit1_yheight
#line 50 "/usr/share/mcstas/2.6rc1/optics/Slit.comp"
{
if (xwidth > 0)  { 
  if (!xmin && !xmax) {
    xmax=xwidth/2;  xmin=-xmax;
  } else {
    fprintf(stderr,"Slit: %s: Error: please specify EITHER xmin & xmax or xwidth\n", NAME_CURRENT_COMP); exit(-1);
  }
 }
 if (yheight > 0) { 
   if (!ymin && !ymax) {
     ymax=yheight/2; ymin=-ymax; 
   } else {
     fprintf(stderr,"Slit: %s: Error: please specify EITHER ymin & ymax or ywidth\n", NAME_CURRENT_COMP); exit(-1);
   }
 }
 if (xmin == 0 && xmax == 0 && ymin == 0 && ymax == 0 && radius == 0)
    { fprintf(stderr,"Slit: %s: Warning: Running with CLOSED slit - is this intentional?? \n", NAME_CURRENT_COMP); }

}
#line 10154 "ISIS_CRISP.c"
#undef yheight
#undef xwidth
#undef radius
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component defslit2. */
  SIG_MESSAGE("defslit2 (Init)");
#define mccompcurname  defslit2
#define mccompcurtype  Slit
#define mccompcurindex 4
#define xmin mccdefslit2_xmin
#define xmax mccdefslit2_xmax
#define ymin mccdefslit2_ymin
#define ymax mccdefslit2_ymax
#define radius mccdefslit2_radius
#define xwidth mccdefslit2_xwidth
#define yheight mccdefslit2_yheight
#line 50 "/usr/share/mcstas/2.6rc1/optics/Slit.comp"
{
if (xwidth > 0)  { 
  if (!xmin && !xmax) {
    xmax=xwidth/2;  xmin=-xmax;
  } else {
    fprintf(stderr,"Slit: %s: Error: please specify EITHER xmin & xmax or xwidth\n", NAME_CURRENT_COMP); exit(-1);
  }
 }
 if (yheight > 0) { 
   if (!ymin && !ymax) {
     ymax=yheight/2; ymin=-ymax; 
   } else {
     fprintf(stderr,"Slit: %s: Error: please specify EITHER ymin & ymax or ywidth\n", NAME_CURRENT_COMP); exit(-1);
   }
 }
 if (xmin == 0 && xmax == 0 && ymin == 0 && ymax == 0 && radius == 0)
    { fprintf(stderr,"Slit: %s: Warning: Running with CLOSED slit - is this intentional?? \n", NAME_CURRENT_COMP); }

}
#line 10198 "ISIS_CRISP.c"
#undef yheight
#undef xwidth
#undef radius
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component coarseslit1. */
  SIG_MESSAGE("coarseslit1 (Init)");
#define mccompcurname  coarseslit1
#define mccompcurtype  Slit
#define mccompcurindex 5
#define xmin mcccoarseslit1_xmin
#define xmax mcccoarseslit1_xmax
#define ymin mcccoarseslit1_ymin
#define ymax mcccoarseslit1_ymax
#define radius mcccoarseslit1_radius
#define xwidth mcccoarseslit1_xwidth
#define yheight mcccoarseslit1_yheight
#line 50 "/usr/share/mcstas/2.6rc1/optics/Slit.comp"
{
if (xwidth > 0)  { 
  if (!xmin && !xmax) {
    xmax=xwidth/2;  xmin=-xmax;
  } else {
    fprintf(stderr,"Slit: %s: Error: please specify EITHER xmin & xmax or xwidth\n", NAME_CURRENT_COMP); exit(-1);
  }
 }
 if (yheight > 0) { 
   if (!ymin && !ymax) {
     ymax=yheight/2; ymin=-ymax; 
   } else {
     fprintf(stderr,"Slit: %s: Error: please specify EITHER ymin & ymax or ywidth\n", NAME_CURRENT_COMP); exit(-1);
   }
 }
 if (xmin == 0 && xmax == 0 && ymin == 0 && ymax == 0 && radius == 0)
    { fprintf(stderr,"Slit: %s: Warning: Running with CLOSED slit - is this intentional?? \n", NAME_CURRENT_COMP); }

}
#line 10242 "ISIS_CRISP.c"
#undef yheight
#undef xwidth
#undef radius
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component coarseslit2. */
  SIG_MESSAGE("coarseslit2 (Init)");
#define mccompcurname  coarseslit2
#define mccompcurtype  Slit
#define mccompcurindex 6
#define xmin mcccoarseslit2_xmin
#define xmax mcccoarseslit2_xmax
#define ymin mcccoarseslit2_ymin
#define ymax mcccoarseslit2_ymax
#define radius mcccoarseslit2_radius
#define xwidth mcccoarseslit2_xwidth
#define yheight mcccoarseslit2_yheight
#line 50 "/usr/share/mcstas/2.6rc1/optics/Slit.comp"
{
if (xwidth > 0)  { 
  if (!xmin && !xmax) {
    xmax=xwidth/2;  xmin=-xmax;
  } else {
    fprintf(stderr,"Slit: %s: Error: please specify EITHER xmin & xmax or xwidth\n", NAME_CURRENT_COMP); exit(-1);
  }
 }
 if (yheight > 0) { 
   if (!ymin && !ymax) {
     ymax=yheight/2; ymin=-ymax; 
   } else {
     fprintf(stderr,"Slit: %s: Error: please specify EITHER ymin & ymax or ywidth\n", NAME_CURRENT_COMP); exit(-1);
   }
 }
 if (xmin == 0 && xmax == 0 && ymin == 0 && ymax == 0 && radius == 0)
    { fprintf(stderr,"Slit: %s: Warning: Running with CLOSED slit - is this intentional?? \n", NAME_CURRENT_COMP); }

}
#line 10286 "ISIS_CRISP.c"
#undef yheight
#undef xwidth
#undef radius
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component lmon1. */
  SIG_MESSAGE("lmon1 (Init)");
#define mccompcurname  lmon1
#define mccompcurtype  L_monitor
#define mccompcurindex 7
#define nL mcclmon1_nL
#define L_N mcclmon1_L_N
#define L_p mcclmon1_L_p
#define L_p2 mcclmon1_L_p2
#define filename mcclmon1_filename
#define xmin mcclmon1_xmin
#define xmax mcclmon1_xmax
#define ymin mcclmon1_ymin
#define ymax mcclmon1_ymax
#define xwidth mcclmon1_xwidth
#define yheight mcclmon1_yheight
#define Lmin mcclmon1_Lmin
#define Lmax mcclmon1_Lmax
#define restore_neutron mcclmon1_restore_neutron
#define nowritefile mcclmon1_nowritefile
#line 62 "/usr/share/mcstas/2.6rc1/monitors/L_monitor.comp"
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
#line 10339 "ISIS_CRISP.c"
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

  /* Initializations for component PSDmon1. */
  SIG_MESSAGE("PSDmon1 (Init)");
#define mccompcurname  PSDmon1
#define mccompcurtype  PSD_monitor
#define mccompcurindex 8
#define PSD_N mccPSDmon1_PSD_N
#define PSD_p mccPSDmon1_PSD_p
#define PSD_p2 mccPSDmon1_PSD_p2
#define nx mccPSDmon1_nx
#define ny mccPSDmon1_ny
#define filename mccPSDmon1_filename
#define xmin mccPSDmon1_xmin
#define xmax mccPSDmon1_xmax
#define ymin mccPSDmon1_ymin
#define ymax mccPSDmon1_ymax
#define xwidth mccPSDmon1_xwidth
#define yheight mccPSDmon1_yheight
#define restore_neutron mccPSDmon1_restore_neutron
#line 68 "/usr/share/mcstas/2.6rc1/monitors/PSD_monitor.comp"
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
#line 10402 "ISIS_CRISP.c"
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

  /* Initializations for component slit1. */
  SIG_MESSAGE("slit1 (Init)");
#define mccompcurname  slit1
#define mccompcurtype  Slit
#define mccompcurindex 9
#define xmin mccslit1_xmin
#define xmax mccslit1_xmax
#define ymin mccslit1_ymin
#define ymax mccslit1_ymax
#define radius mccslit1_radius
#define xwidth mccslit1_xwidth
#define yheight mccslit1_yheight
#line 50 "/usr/share/mcstas/2.6rc1/optics/Slit.comp"
{
if (xwidth > 0)  { 
  if (!xmin && !xmax) {
    xmax=xwidth/2;  xmin=-xmax;
  } else {
    fprintf(stderr,"Slit: %s: Error: please specify EITHER xmin & xmax or xwidth\n", NAME_CURRENT_COMP); exit(-1);
  }
 }
 if (yheight > 0) { 
   if (!ymin && !ymax) {
     ymax=yheight/2; ymin=-ymax; 
   } else {
     fprintf(stderr,"Slit: %s: Error: please specify EITHER ymin & ymax or ywidth\n", NAME_CURRENT_COMP); exit(-1);
   }
 }
 if (xmin == 0 && xmax == 0 && ymin == 0 && ymax == 0 && radius == 0)
    { fprintf(stderr,"Slit: %s: Warning: Running with CLOSED slit - is this intentional?? \n", NAME_CURRENT_COMP); }

}
#line 10452 "ISIS_CRISP.c"
#undef yheight
#undef xwidth
#undef radius
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component PSDmons1. */
  SIG_MESSAGE("PSDmons1 (Init)");
#define mccompcurname  PSDmons1
#define mccompcurtype  PSD_monitor
#define mccompcurindex 10
#define PSD_N mccPSDmons1_PSD_N
#define PSD_p mccPSDmons1_PSD_p
#define PSD_p2 mccPSDmons1_PSD_p2
#define nx mccPSDmons1_nx
#define ny mccPSDmons1_ny
#define filename mccPSDmons1_filename
#define xmin mccPSDmons1_xmin
#define xmax mccPSDmons1_xmax
#define ymin mccPSDmons1_ymin
#define ymax mccPSDmons1_ymax
#define xwidth mccPSDmons1_xwidth
#define yheight mccPSDmons1_yheight
#define restore_neutron mccPSDmons1_restore_neutron
#line 68 "/usr/share/mcstas/2.6rc1/monitors/PSD_monitor.comp"
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
#line 10507 "ISIS_CRISP.c"
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

  /* Initializations for component eguide1. */
  SIG_MESSAGE("eguide1 (Init)");
#define mccompcurname  eguide1
#define mccompcurtype  Guide_tapering
#define mccompcurindex 11
#define w1c mcceguide1_w1c
#define w2c mcceguide1_w2c
#define ww mcceguide1_ww
#define hh mcceguide1_hh
#define whalf mcceguide1_whalf
#define hhalf mcceguide1_hhalf
#define lwhalf mcceguide1_lwhalf
#define lhhalf mcceguide1_lhhalf
#define h1_in mcceguide1_h1_in
#define h2_out mcceguide1_h2_out
#define w1_in mcceguide1_w1_in
#define w2_out mcceguide1_w2_out
#define l_seg mcceguide1_l_seg
#define seg mcceguide1_seg
#define h12 mcceguide1_h12
#define h2 mcceguide1_h2
#define w12 mcceguide1_w12
#define w2 mcceguide1_w2
#define a_ell_q mcceguide1_a_ell_q
#define b_ell_q mcceguide1_b_ell_q
#define lbw mcceguide1_lbw
#define lbh mcceguide1_lbh
#define mxi mcceguide1_mxi
#define u1 mcceguide1_u1
#define u2 mcceguide1_u2
#define div1 mcceguide1_div1
#define p2_para mcceguide1_p2_para
#define test mcceguide1_test
#define Div1 mcceguide1_Div1
#define i mcceguide1_i
#define ii mcceguide1_ii
#define seg mcceguide1_seg
#define fu mcceguide1_fu
#define pos mcceguide1_pos
#define file_name mcceguide1_file_name
#define ep mcceguide1_ep
#define num mcceguide1_num
#define rotation_h mcceguide1_rotation_h
#define rotation_v mcceguide1_rotation_v
#define option mcceguide1_option
#define w1 mcceguide1_w1
#define h1 mcceguide1_h1
#define l mcceguide1_l
#define linw mcceguide1_linw
#define loutw mcceguide1_loutw
#define linh mcceguide1_linh
#define louth mcceguide1_louth
#define R0 mcceguide1_R0
#define Qcx mcceguide1_Qcx
#define Qcy mcceguide1_Qcy
#define alphax mcceguide1_alphax
#define alphay mcceguide1_alphay
#define W mcceguide1_W
#define mx mcceguide1_mx
#define my mcceguide1_my
#define segno mcceguide1_segno
#define curvature mcceguide1_curvature
#define curvature_v mcceguide1_curvature_v
#line 115 "/usr/share/mcstas/2.6rc1/optics/Guide_tapering.comp"
{
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
          fgets(segment[ii].st,128,num);
          if (ii >  799) {
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
}
#line 10924 "ISIS_CRISP.c"
#undef curvature_v
#undef curvature
#undef segno
#undef my
#undef mx
#undef W
#undef alphay
#undef alphax
#undef Qcy
#undef Qcx
#undef R0
#undef louth
#undef linh
#undef loutw
#undef linw
#undef l
#undef h1
#undef w1
#undef option
#undef rotation_v
#undef rotation_h
#undef num
#undef ep
#undef file_name
#undef pos
#undef fu
#undef seg
#undef ii
#undef i
#undef Div1
#undef test
#undef p2_para
#undef div1
#undef u2
#undef u1
#undef mxi
#undef lbh
#undef lbw
#undef b_ell_q
#undef a_ell_q
#undef w2
#undef w12
#undef h2
#undef h12
#undef seg
#undef l_seg
#undef w2_out
#undef w1_in
#undef h2_out
#undef h1_in
#undef lhhalf
#undef lwhalf
#undef hhalf
#undef whalf
#undef hh
#undef ww
#undef w2c
#undef w1c
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component slit2. */
  SIG_MESSAGE("slit2 (Init)");
#define mccompcurname  slit2
#define mccompcurtype  Slit
#define mccompcurindex 12
#define xmin mccslit2_xmin
#define xmax mccslit2_xmax
#define ymin mccslit2_ymin
#define ymax mccslit2_ymax
#define radius mccslit2_radius
#define xwidth mccslit2_xwidth
#define yheight mccslit2_yheight
#line 50 "/usr/share/mcstas/2.6rc1/optics/Slit.comp"
{
if (xwidth > 0)  { 
  if (!xmin && !xmax) {
    xmax=xwidth/2;  xmin=-xmax;
  } else {
    fprintf(stderr,"Slit: %s: Error: please specify EITHER xmin & xmax or xwidth\n", NAME_CURRENT_COMP); exit(-1);
  }
 }
 if (yheight > 0) { 
   if (!ymin && !ymax) {
     ymax=yheight/2; ymin=-ymax; 
   } else {
     fprintf(stderr,"Slit: %s: Error: please specify EITHER ymin & ymax or ywidth\n", NAME_CURRENT_COMP); exit(-1);
   }
 }
 if (xmin == 0 && xmax == 0 && ymin == 0 && ymax == 0 && radius == 0)
    { fprintf(stderr,"Slit: %s: Warning: Running with CLOSED slit - is this intentional?? \n", NAME_CURRENT_COMP); }

}
#line 11019 "ISIS_CRISP.c"
#undef yheight
#undef xwidth
#undef radius
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component lmon2. */
  SIG_MESSAGE("lmon2 (Init)");
#define mccompcurname  lmon2
#define mccompcurtype  L_monitor
#define mccompcurindex 13
#define nL mcclmon2_nL
#define L_N mcclmon2_L_N
#define L_p mcclmon2_L_p
#define L_p2 mcclmon2_L_p2
#define filename mcclmon2_filename
#define xmin mcclmon2_xmin
#define xmax mcclmon2_xmax
#define ymin mcclmon2_ymin
#define ymax mcclmon2_ymax
#define xwidth mcclmon2_xwidth
#define yheight mcclmon2_yheight
#define Lmin mcclmon2_Lmin
#define Lmax mcclmon2_Lmax
#define restore_neutron mcclmon2_restore_neutron
#define nowritefile mcclmon2_nowritefile
#line 62 "/usr/share/mcstas/2.6rc1/monitors/L_monitor.comp"
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
#line 11072 "ISIS_CRISP.c"
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

  /* Initializations for component samp1. */
  SIG_MESSAGE("samp1 (Init)");
#define mccompcurname  samp1
#define mccompcurtype  Multilayer_Sample
#define mccompcurindex 14
#define xwidth mccsamp1_xwidth
#define zlength mccsamp1_zlength
#define nlayer mccsamp1_nlayer
#define sldPar mccsamp1_sldPar
#define dPar mccsamp1_dPar
#define sigmaPar mccsamp1_sigmaPar
#define frac_inc mccsamp1_frac_inc
#define ythick mccsamp1_ythick
#define mu_inc mccsamp1_mu_inc
#define target_index mccsamp1_target_index
#define focus_xw mccsamp1_focus_xw
#define focus_yh mccsamp1_focus_yh
#define sldParPtr mccsamp1_sldParPtr
#define dParPtr mccsamp1_dParPtr
#define sigmaParPtr mccsamp1_sigmaParPtr
#define tx mccsamp1_tx
#define ty mccsamp1_ty
#define tz mccsamp1_tz
#line 78 "/usr/share/mcstas/2.6rc1/contrib/Multilayer_Sample.comp"
{
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

}
#line 11137 "ISIS_CRISP.c"
#undef tz
#undef ty
#undef tx
#undef sigmaParPtr
#undef dParPtr
#undef sldParPtr
#undef focus_yh
#undef focus_xw
#undef target_index
#undef mu_inc
#undef ythick
#undef frac_inc
#undef sigmaPar
#undef dPar
#undef sldPar
#undef nlayer
#undef zlength
#undef xwidth
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component slit3. */
  SIG_MESSAGE("slit3 (Init)");
#define mccompcurname  slit3
#define mccompcurtype  Slit
#define mccompcurindex 15
#define xmin mccslit3_xmin
#define xmax mccslit3_xmax
#define ymin mccslit3_ymin
#define ymax mccslit3_ymax
#define radius mccslit3_radius
#define xwidth mccslit3_xwidth
#define yheight mccslit3_yheight
#line 50 "/usr/share/mcstas/2.6rc1/optics/Slit.comp"
{
if (xwidth > 0)  { 
  if (!xmin && !xmax) {
    xmax=xwidth/2;  xmin=-xmax;
  } else {
    fprintf(stderr,"Slit: %s: Error: please specify EITHER xmin & xmax or xwidth\n", NAME_CURRENT_COMP); exit(-1);
  }
 }
 if (yheight > 0) { 
   if (!ymin && !ymax) {
     ymax=yheight/2; ymin=-ymax; 
   } else {
     fprintf(stderr,"Slit: %s: Error: please specify EITHER ymin & ymax or ywidth\n", NAME_CURRENT_COMP); exit(-1);
   }
 }
 if (xmin == 0 && xmax == 0 && ymin == 0 && ymax == 0 && radius == 0)
    { fprintf(stderr,"Slit: %s: Warning: Running with CLOSED slit - is this intentional?? \n", NAME_CURRENT_COMP); }

}
#line 11192 "ISIS_CRISP.c"
#undef yheight
#undef xwidth
#undef radius
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component eguide2. */
  SIG_MESSAGE("eguide2 (Init)");
#define mccompcurname  eguide2
#define mccompcurtype  Guide_tapering
#define mccompcurindex 16
#define w1c mcceguide2_w1c
#define w2c mcceguide2_w2c
#define ww mcceguide2_ww
#define hh mcceguide2_hh
#define whalf mcceguide2_whalf
#define hhalf mcceguide2_hhalf
#define lwhalf mcceguide2_lwhalf
#define lhhalf mcceguide2_lhhalf
#define h1_in mcceguide2_h1_in
#define h2_out mcceguide2_h2_out
#define w1_in mcceguide2_w1_in
#define w2_out mcceguide2_w2_out
#define l_seg mcceguide2_l_seg
#define seg mcceguide2_seg
#define h12 mcceguide2_h12
#define h2 mcceguide2_h2
#define w12 mcceguide2_w12
#define w2 mcceguide2_w2
#define a_ell_q mcceguide2_a_ell_q
#define b_ell_q mcceguide2_b_ell_q
#define lbw mcceguide2_lbw
#define lbh mcceguide2_lbh
#define mxi mcceguide2_mxi
#define u1 mcceguide2_u1
#define u2 mcceguide2_u2
#define div1 mcceguide2_div1
#define p2_para mcceguide2_p2_para
#define test mcceguide2_test
#define Div1 mcceguide2_Div1
#define i mcceguide2_i
#define ii mcceguide2_ii
#define seg mcceguide2_seg
#define fu mcceguide2_fu
#define pos mcceguide2_pos
#define file_name mcceguide2_file_name
#define ep mcceguide2_ep
#define num mcceguide2_num
#define rotation_h mcceguide2_rotation_h
#define rotation_v mcceguide2_rotation_v
#define option mcceguide2_option
#define w1 mcceguide2_w1
#define h1 mcceguide2_h1
#define l mcceguide2_l
#define linw mcceguide2_linw
#define loutw mcceguide2_loutw
#define linh mcceguide2_linh
#define louth mcceguide2_louth
#define R0 mcceguide2_R0
#define Qcx mcceguide2_Qcx
#define Qcy mcceguide2_Qcy
#define alphax mcceguide2_alphax
#define alphay mcceguide2_alphay
#define W mcceguide2_W
#define mx mcceguide2_mx
#define my mcceguide2_my
#define segno mcceguide2_segno
#define curvature mcceguide2_curvature
#define curvature_v mcceguide2_curvature_v
#line 115 "/usr/share/mcstas/2.6rc1/optics/Guide_tapering.comp"
{
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
          fgets(segment[ii].st,128,num);
          if (ii >  799) {
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
}
#line 11603 "ISIS_CRISP.c"
#undef curvature_v
#undef curvature
#undef segno
#undef my
#undef mx
#undef W
#undef alphay
#undef alphax
#undef Qcy
#undef Qcx
#undef R0
#undef louth
#undef linh
#undef loutw
#undef linw
#undef l
#undef h1
#undef w1
#undef option
#undef rotation_v
#undef rotation_h
#undef num
#undef ep
#undef file_name
#undef pos
#undef fu
#undef seg
#undef ii
#undef i
#undef Div1
#undef test
#undef p2_para
#undef div1
#undef u2
#undef u1
#undef mxi
#undef lbh
#undef lbw
#undef b_ell_q
#undef a_ell_q
#undef w2
#undef w12
#undef h2
#undef h12
#undef seg
#undef l_seg
#undef w2_out
#undef w1_in
#undef h2_out
#undef h1_in
#undef lhhalf
#undef lwhalf
#undef hhalf
#undef whalf
#undef hh
#undef ww
#undef w2c
#undef w1c
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component lmon3. */
  SIG_MESSAGE("lmon3 (Init)");
#define mccompcurname  lmon3
#define mccompcurtype  L_monitor
#define mccompcurindex 17
#define nL mcclmon3_nL
#define L_N mcclmon3_L_N
#define L_p mcclmon3_L_p
#define L_p2 mcclmon3_L_p2
#define filename mcclmon3_filename
#define xmin mcclmon3_xmin
#define xmax mcclmon3_xmax
#define ymin mcclmon3_ymin
#define ymax mcclmon3_ymax
#define xwidth mcclmon3_xwidth
#define yheight mcclmon3_yheight
#define Lmin mcclmon3_Lmin
#define Lmax mcclmon3_Lmax
#define restore_neutron mcclmon3_restore_neutron
#define nowritefile mcclmon3_nowritefile
#line 62 "/usr/share/mcstas/2.6rc1/monitors/L_monitor.comp"
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
#line 11707 "ISIS_CRISP.c"
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

  /* Initializations for component PSDdet1. */
  SIG_MESSAGE("PSDdet1 (Init)");
#define mccompcurname  PSDdet1
#define mccompcurtype  PSD_monitor
#define mccompcurindex 18
#define PSD_N mccPSDdet1_PSD_N
#define PSD_p mccPSDdet1_PSD_p
#define PSD_p2 mccPSDdet1_PSD_p2
#define nx mccPSDdet1_nx
#define ny mccPSDdet1_ny
#define filename mccPSDdet1_filename
#define xmin mccPSDdet1_xmin
#define xmax mccPSDdet1_xmax
#define ymin mccPSDdet1_ymin
#define ymax mccPSDdet1_ymax
#define xwidth mccPSDdet1_xwidth
#define yheight mccPSDdet1_yheight
#define restore_neutron mccPSDdet1_restore_neutron
#line 68 "/usr/share/mcstas/2.6rc1/monitors/PSD_monitor.comp"
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
#line 11770 "ISIS_CRISP.c"
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
  /* TRACE Component a1 [1] */
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
#define mccompcurname  a1
#define mccompcurtype  Arm
#define mccompcurindex 1
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompa1:
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

  /* TRACE Component isis_source [2] */
  mccoordschange(mcposrisis_source, mcrotrisis_source,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component isis_source (without coords transformations) */
  mcJumpTrace_isis_source:
  SIG_MESSAGE("isis_source (Trace)");
  mcDEBUG_COMP("isis_source")
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

#define mcabsorbComp mcabsorbCompisis_source
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
#define mccompcurname  isis_source
#define mccompcurtype  ISIS_moderator
#define mccompcurindex 2
#define Face mccisis_source_Face
#define p_in mccisis_source_p_in
#define Tnpts mccisis_source_Tnpts
#define scaleSize mccisis_source_scaleSize
#define angleArea mccisis_source_angleArea
#define Nsim mccisis_source_Nsim
#define Ncount mccisis_source_Ncount
#define TS mccisis_source_TS
#define rtE0 mccisis_source_rtE0
#define rtE1 mccisis_source_rtE1
#define rtmodX mccisis_source_rtmodX
#define rtmodY mccisis_source_rtmodY
#define TargetStation mccisis_source_TargetStation
#define CurrentWeight mccisis_source_CurrentWeight
{   /* Declarations of isis_source=ISIS_moderator() SETTING parameters. */
MCNUM Emin = mccisis_source_Emin;
MCNUM Emax = mccisis_source_Emax;
MCNUM dist = mccisis_source_dist;
MCNUM focus_xw = mccisis_source_focus_xw;
MCNUM focus_yh = mccisis_source_focus_yh;
MCNUM xwidth = mccisis_source_xwidth;
MCNUM yheight = mccisis_source_yheight;
MCNUM CAngle = mccisis_source_CAngle;
MCNUM SAC = mccisis_source_SAC;
MCNUM Lmin = mccisis_source_Lmin;
MCNUM Lmax = mccisis_source_Lmax;
int target_index = mccisis_source_target_index;
MCNUM verbose = mccisis_source_verbose;
#line 1060 "/usr/share/mcstas/2.6rc1/contrib/ISIS_moderator.comp"
{
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
}
#line 12076 "ISIS_CRISP.c"
}   /* End of isis_source=ISIS_moderator() SETTING parameter declarations. */
#undef CurrentWeight
#undef TargetStation
#undef rtmodY
#undef rtmodX
#undef rtE1
#undef rtE0
#undef TS
#undef Ncount
#undef Nsim
#undef angleArea
#undef scaleSize
#undef Tnpts
#undef p_in
#undef Face
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompisis_source:
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

  /* TRACE Component defslit1 [3] */
  mccoordschange(mcposrdefslit1, mcrotrdefslit1,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component defslit1 (without coords transformations) */
  mcJumpTrace_defslit1:
  SIG_MESSAGE("defslit1 (Trace)");
  mcDEBUG_COMP("defslit1")
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

#define mcabsorbComp mcabsorbCompdefslit1
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
#define mccompcurname  defslit1
#define mccompcurtype  Slit
#define mccompcurindex 3
{   /* Declarations of defslit1=Slit() SETTING parameters. */
MCNUM xmin = mccdefslit1_xmin;
MCNUM xmax = mccdefslit1_xmax;
MCNUM ymin = mccdefslit1_ymin;
MCNUM ymax = mccdefslit1_ymax;
MCNUM radius = mccdefslit1_radius;
MCNUM xwidth = mccdefslit1_xwidth;
MCNUM yheight = mccdefslit1_yheight;
#line 71 "/usr/share/mcstas/2.6rc1/optics/Slit.comp"
{
    PROP_Z0;
    if (((radius == 0) && (x<xmin || x>xmax || y<ymin || y>ymax))
	|| ((radius != 0) && (x*x + y*y > radius*radius))) {
      RESTORE_NEUTRON(INDEX_CURRENT_COMP, x, y, z, vx, vy, vz, t, sx, sy, sz, p);
      ABSORB;
    }
    else
        SCATTER;
}
#line 12214 "ISIS_CRISP.c"
}   /* End of defslit1=Slit() SETTING parameter declarations. */
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompdefslit1:
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

  /* TRACE Component defslit2 [4] */
  mccoordschange(mcposrdefslit2, mcrotrdefslit2,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component defslit2 (without coords transformations) */
  mcJumpTrace_defslit2:
  SIG_MESSAGE("defslit2 (Trace)");
  mcDEBUG_COMP("defslit2")
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

#define mcabsorbComp mcabsorbCompdefslit2
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
#define mccompcurname  defslit2
#define mccompcurtype  Slit
#define mccompcurindex 4
{   /* Declarations of defslit2=Slit() SETTING parameters. */
MCNUM xmin = mccdefslit2_xmin;
MCNUM xmax = mccdefslit2_xmax;
MCNUM ymin = mccdefslit2_ymin;
MCNUM ymax = mccdefslit2_ymax;
MCNUM radius = mccdefslit2_radius;
MCNUM xwidth = mccdefslit2_xwidth;
MCNUM yheight = mccdefslit2_yheight;
#line 71 "/usr/share/mcstas/2.6rc1/optics/Slit.comp"
{
    PROP_Z0;
    if (((radius == 0) && (x<xmin || x>xmax || y<ymin || y>ymax))
	|| ((radius != 0) && (x*x + y*y > radius*radius))) {
      RESTORE_NEUTRON(INDEX_CURRENT_COMP, x, y, z, vx, vy, vz, t, sx, sy, sz, p);
      ABSORB;
    }
    else
        SCATTER;
}
#line 12338 "ISIS_CRISP.c"
}   /* End of defslit2=Slit() SETTING parameter declarations. */
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompdefslit2:
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

  /* TRACE Component coarseslit1 [5] */
  mccoordschange(mcposrcoarseslit1, mcrotrcoarseslit1,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component coarseslit1 (without coords transformations) */
  mcJumpTrace_coarseslit1:
  SIG_MESSAGE("coarseslit1 (Trace)");
  mcDEBUG_COMP("coarseslit1")
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

#define mcabsorbComp mcabsorbCompcoarseslit1
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
#define mccompcurname  coarseslit1
#define mccompcurtype  Slit
#define mccompcurindex 5
{   /* Declarations of coarseslit1=Slit() SETTING parameters. */
MCNUM xmin = mcccoarseslit1_xmin;
MCNUM xmax = mcccoarseslit1_xmax;
MCNUM ymin = mcccoarseslit1_ymin;
MCNUM ymax = mcccoarseslit1_ymax;
MCNUM radius = mcccoarseslit1_radius;
MCNUM xwidth = mcccoarseslit1_xwidth;
MCNUM yheight = mcccoarseslit1_yheight;
#line 71 "/usr/share/mcstas/2.6rc1/optics/Slit.comp"
{
    PROP_Z0;
    if (((radius == 0) && (x<xmin || x>xmax || y<ymin || y>ymax))
	|| ((radius != 0) && (x*x + y*y > radius*radius))) {
      RESTORE_NEUTRON(INDEX_CURRENT_COMP, x, y, z, vx, vy, vz, t, sx, sy, sz, p);
      ABSORB;
    }
    else
        SCATTER;
}
#line 12462 "ISIS_CRISP.c"
}   /* End of coarseslit1=Slit() SETTING parameter declarations. */
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompcoarseslit1:
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

  /* TRACE Component coarseslit2 [6] */
  mccoordschange(mcposrcoarseslit2, mcrotrcoarseslit2,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component coarseslit2 (without coords transformations) */
  mcJumpTrace_coarseslit2:
  SIG_MESSAGE("coarseslit2 (Trace)");
  mcDEBUG_COMP("coarseslit2")
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

#define mcabsorbComp mcabsorbCompcoarseslit2
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
#define mccompcurname  coarseslit2
#define mccompcurtype  Slit
#define mccompcurindex 6
{   /* Declarations of coarseslit2=Slit() SETTING parameters. */
MCNUM xmin = mcccoarseslit2_xmin;
MCNUM xmax = mcccoarseslit2_xmax;
MCNUM ymin = mcccoarseslit2_ymin;
MCNUM ymax = mcccoarseslit2_ymax;
MCNUM radius = mcccoarseslit2_radius;
MCNUM xwidth = mcccoarseslit2_xwidth;
MCNUM yheight = mcccoarseslit2_yheight;
#line 71 "/usr/share/mcstas/2.6rc1/optics/Slit.comp"
{
    PROP_Z0;
    if (((radius == 0) && (x<xmin || x>xmax || y<ymin || y>ymax))
	|| ((radius != 0) && (x*x + y*y > radius*radius))) {
      RESTORE_NEUTRON(INDEX_CURRENT_COMP, x, y, z, vx, vy, vz, t, sx, sy, sz, p);
      ABSORB;
    }
    else
        SCATTER;
}
#line 12586 "ISIS_CRISP.c"
}   /* End of coarseslit2=Slit() SETTING parameter declarations. */
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompcoarseslit2:
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

  /* TRACE Component lmon1 [7] */
  mccoordschange(mcposrlmon1, mcrotrlmon1,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component lmon1 (without coords transformations) */
  mcJumpTrace_lmon1:
  SIG_MESSAGE("lmon1 (Trace)");
  mcDEBUG_COMP("lmon1")
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

#define mcabsorbComp mcabsorbComplmon1
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
#define mccompcurname  lmon1
#define mccompcurtype  L_monitor
#define mccompcurindex 7
#define nL mcclmon1_nL
#define L_N mcclmon1_L_N
#define L_p mcclmon1_L_p
#define L_p2 mcclmon1_L_p2
{   /* Declarations of lmon1=L_monitor() SETTING parameters. */
char* filename = mcclmon1_filename;
MCNUM xmin = mcclmon1_xmin;
MCNUM xmax = mcclmon1_xmax;
MCNUM ymin = mcclmon1_ymin;
MCNUM ymax = mcclmon1_ymax;
MCNUM xwidth = mcclmon1_xwidth;
MCNUM yheight = mcclmon1_yheight;
MCNUM Lmin = mcclmon1_Lmin;
MCNUM Lmax = mcclmon1_Lmax;
MCNUM restore_neutron = mcclmon1_restore_neutron;
int nowritefile = mcclmon1_nowritefile;
#line 84 "/usr/share/mcstas/2.6rc1/monitors/L_monitor.comp"
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
#line 12729 "ISIS_CRISP.c"
}   /* End of lmon1=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbComplmon1:
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

  /* TRACE Component PSDmon1 [8] */
  mccoordschange(mcposrPSDmon1, mcrotrPSDmon1,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component PSDmon1 (without coords transformations) */
  mcJumpTrace_PSDmon1:
  SIG_MESSAGE("PSDmon1 (Trace)");
  mcDEBUG_COMP("PSDmon1")
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

#define mcabsorbComp mcabsorbCompPSDmon1
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
#define mccompcurname  PSDmon1
#define mccompcurtype  PSD_monitor
#define mccompcurindex 8
#define PSD_N mccPSDmon1_PSD_N
#define PSD_p mccPSDmon1_PSD_p
#define PSD_p2 mccPSDmon1_PSD_p2
{   /* Declarations of PSDmon1=PSD_monitor() SETTING parameters. */
int nx = mccPSDmon1_nx;
int ny = mccPSDmon1_ny;
char* filename = mccPSDmon1_filename;
MCNUM xmin = mccPSDmon1_xmin;
MCNUM xmax = mccPSDmon1_xmax;
MCNUM ymin = mccPSDmon1_ymin;
MCNUM ymax = mccPSDmon1_ymax;
MCNUM xwidth = mccPSDmon1_xwidth;
MCNUM yheight = mccPSDmon1_yheight;
MCNUM restore_neutron = mccPSDmon1_restore_neutron;
#line 94 "/usr/share/mcstas/2.6rc1/monitors/PSD_monitor.comp"
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
#line 12867 "ISIS_CRISP.c"
}   /* End of PSDmon1=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompPSDmon1:
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

  /* TRACE Component slit1 [9] */
  mccoordschange(mcposrslit1, mcrotrslit1,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component slit1 (without coords transformations) */
  mcJumpTrace_slit1:
  SIG_MESSAGE("slit1 (Trace)");
  mcDEBUG_COMP("slit1")
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

#define mcabsorbComp mcabsorbCompslit1
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
#define mccompcurname  slit1
#define mccompcurtype  Slit
#define mccompcurindex 9
{   /* Declarations of slit1=Slit() SETTING parameters. */
MCNUM xmin = mccslit1_xmin;
MCNUM xmax = mccslit1_xmax;
MCNUM ymin = mccslit1_ymin;
MCNUM ymax = mccslit1_ymax;
MCNUM radius = mccslit1_radius;
MCNUM xwidth = mccslit1_xwidth;
MCNUM yheight = mccslit1_yheight;
#line 71 "/usr/share/mcstas/2.6rc1/optics/Slit.comp"
{
    PROP_Z0;
    if (((radius == 0) && (x<xmin || x>xmax || y<ymin || y>ymax))
	|| ((radius != 0) && (x*x + y*y > radius*radius))) {
      RESTORE_NEUTRON(INDEX_CURRENT_COMP, x, y, z, vx, vy, vz, t, sx, sy, sz, p);
      ABSORB;
    }
    else
        SCATTER;
}
#line 12994 "ISIS_CRISP.c"
}   /* End of slit1=Slit() SETTING parameter declarations. */
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompslit1:
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

  /* TRACE Component PSDmons1 [10] */
  mccoordschange(mcposrPSDmons1, mcrotrPSDmons1,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component PSDmons1 (without coords transformations) */
  mcJumpTrace_PSDmons1:
  SIG_MESSAGE("PSDmons1 (Trace)");
  mcDEBUG_COMP("PSDmons1")
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

#define mcabsorbComp mcabsorbCompPSDmons1
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
#define mccompcurname  PSDmons1
#define mccompcurtype  PSD_monitor
#define mccompcurindex 10
#define PSD_N mccPSDmons1_PSD_N
#define PSD_p mccPSDmons1_PSD_p
#define PSD_p2 mccPSDmons1_PSD_p2
{   /* Declarations of PSDmons1=PSD_monitor() SETTING parameters. */
int nx = mccPSDmons1_nx;
int ny = mccPSDmons1_ny;
char* filename = mccPSDmons1_filename;
MCNUM xmin = mccPSDmons1_xmin;
MCNUM xmax = mccPSDmons1_xmax;
MCNUM ymin = mccPSDmons1_ymin;
MCNUM ymax = mccPSDmons1_ymax;
MCNUM xwidth = mccPSDmons1_xwidth;
MCNUM yheight = mccPSDmons1_yheight;
MCNUM restore_neutron = mccPSDmons1_restore_neutron;
#line 94 "/usr/share/mcstas/2.6rc1/monitors/PSD_monitor.comp"
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
#line 13128 "ISIS_CRISP.c"
}   /* End of PSDmons1=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompPSDmons1:
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

  /* TRACE Component eguide1 [11] */
  mccoordschange(mcposreguide1, mcrotreguide1,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component eguide1 (without coords transformations) */
  mcJumpTrace_eguide1:
  SIG_MESSAGE("eguide1 (Trace)");
  mcDEBUG_COMP("eguide1")
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

#define mcabsorbComp mcabsorbCompeguide1
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
#define mccompcurname  eguide1
#define mccompcurtype  Guide_tapering
#define mccompcurindex 11
#define w1c mcceguide1_w1c
#define w2c mcceguide1_w2c
#define ww mcceguide1_ww
#define hh mcceguide1_hh
#define whalf mcceguide1_whalf
#define hhalf mcceguide1_hhalf
#define lwhalf mcceguide1_lwhalf
#define lhhalf mcceguide1_lhhalf
#define h1_in mcceguide1_h1_in
#define h2_out mcceguide1_h2_out
#define w1_in mcceguide1_w1_in
#define w2_out mcceguide1_w2_out
#define l_seg mcceguide1_l_seg
#define seg mcceguide1_seg
#define h12 mcceguide1_h12
#define h2 mcceguide1_h2
#define w12 mcceguide1_w12
#define w2 mcceguide1_w2
#define a_ell_q mcceguide1_a_ell_q
#define b_ell_q mcceguide1_b_ell_q
#define lbw mcceguide1_lbw
#define lbh mcceguide1_lbh
#define mxi mcceguide1_mxi
#define u1 mcceguide1_u1
#define u2 mcceguide1_u2
#define div1 mcceguide1_div1
#define p2_para mcceguide1_p2_para
#define test mcceguide1_test
#define Div1 mcceguide1_Div1
#define i mcceguide1_i
#define ii mcceguide1_ii
#define seg mcceguide1_seg
#define fu mcceguide1_fu
#define pos mcceguide1_pos
#define file_name mcceguide1_file_name
#define ep mcceguide1_ep
#define num mcceguide1_num
#define rotation_h mcceguide1_rotation_h
#define rotation_v mcceguide1_rotation_v
{   /* Declarations of eguide1=Guide_tapering() SETTING parameters. */
char* option = mcceguide1_option;
MCNUM w1 = mcceguide1_w1;
MCNUM h1 = mcceguide1_h1;
MCNUM l = mcceguide1_l;
MCNUM linw = mcceguide1_linw;
MCNUM loutw = mcceguide1_loutw;
MCNUM linh = mcceguide1_linh;
MCNUM louth = mcceguide1_louth;
MCNUM R0 = mcceguide1_R0;
MCNUM Qcx = mcceguide1_Qcx;
MCNUM Qcy = mcceguide1_Qcy;
MCNUM alphax = mcceguide1_alphax;
MCNUM alphay = mcceguide1_alphay;
MCNUM W = mcceguide1_W;
MCNUM mx = mcceguide1_mx;
MCNUM my = mcceguide1_my;
MCNUM segno = mcceguide1_segno;
MCNUM curvature = mcceguide1_curvature;
MCNUM curvature_v = mcceguide1_curvature_v;
#line 452 "/usr/share/mcstas/2.6rc1/optics/Guide_tapering.comp"
{
  double t1,t2,ts,zr;                           /* Intersection times. */
  double av,ah,bv,bh,cv1,cv2,ch1,ch2,dd;        /* Intermediate values */
  double vdotn_v1,vdotn_v2,vdotn_h1,vdotn_h2;   /* Dot products. */
  int i;                                        /* Which mirror hit? */
  double q;                                     /* Q [1/AA] of reflection */
  double vlen2,nlen2;                           /* Vector lengths squared */
  double edge;
  double hadj;                                  /* Channel displacement */
  double sz=0;
  int ii;
  Coords mctc1,mctc2;
  Rotation mctr1;

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

}
#line 13450 "ISIS_CRISP.c"
}   /* End of eguide1=Guide_tapering() SETTING parameter declarations. */
#undef rotation_v
#undef rotation_h
#undef num
#undef ep
#undef file_name
#undef pos
#undef fu
#undef seg
#undef ii
#undef i
#undef Div1
#undef test
#undef p2_para
#undef div1
#undef u2
#undef u1
#undef mxi
#undef lbh
#undef lbw
#undef b_ell_q
#undef a_ell_q
#undef w2
#undef w12
#undef h2
#undef h12
#undef seg
#undef l_seg
#undef w2_out
#undef w1_in
#undef h2_out
#undef h1_in
#undef lhhalf
#undef lwhalf
#undef hhalf
#undef whalf
#undef hh
#undef ww
#undef w2c
#undef w1c
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompeguide1:
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

  /* TRACE Component slit2 [12] */
  mccoordschange(mcposrslit2, mcrotrslit2,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component slit2 (without coords transformations) */
  mcJumpTrace_slit2:
  SIG_MESSAGE("slit2 (Trace)");
  mcDEBUG_COMP("slit2")
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

#define mcabsorbComp mcabsorbCompslit2
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
#define mccompcurname  slit2
#define mccompcurtype  Slit
#define mccompcurindex 12
{   /* Declarations of slit2=Slit() SETTING parameters. */
MCNUM xmin = mccslit2_xmin;
MCNUM xmax = mccslit2_xmax;
MCNUM ymin = mccslit2_ymin;
MCNUM ymax = mccslit2_ymax;
MCNUM radius = mccslit2_radius;
MCNUM xwidth = mccslit2_xwidth;
MCNUM yheight = mccslit2_yheight;
#line 71 "/usr/share/mcstas/2.6rc1/optics/Slit.comp"
{
    PROP_Z0;
    if (((radius == 0) && (x<xmin || x>xmax || y<ymin || y>ymax))
	|| ((radius != 0) && (x*x + y*y > radius*radius))) {
      RESTORE_NEUTRON(INDEX_CURRENT_COMP, x, y, z, vx, vy, vz, t, sx, sy, sz, p);
      ABSORB;
    }
    else
        SCATTER;
}
#line 13613 "ISIS_CRISP.c"
}   /* End of slit2=Slit() SETTING parameter declarations. */
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompslit2:
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

  /* TRACE Component lmon2 [13] */
  mccoordschange(mcposrlmon2, mcrotrlmon2,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component lmon2 (without coords transformations) */
  mcJumpTrace_lmon2:
  SIG_MESSAGE("lmon2 (Trace)");
  mcDEBUG_COMP("lmon2")
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

#define mcabsorbComp mcabsorbComplmon2
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
#define mccompcurname  lmon2
#define mccompcurtype  L_monitor
#define mccompcurindex 13
#define nL mcclmon2_nL
#define L_N mcclmon2_L_N
#define L_p mcclmon2_L_p
#define L_p2 mcclmon2_L_p2
{   /* Declarations of lmon2=L_monitor() SETTING parameters. */
char* filename = mcclmon2_filename;
MCNUM xmin = mcclmon2_xmin;
MCNUM xmax = mcclmon2_xmax;
MCNUM ymin = mcclmon2_ymin;
MCNUM ymax = mcclmon2_ymax;
MCNUM xwidth = mcclmon2_xwidth;
MCNUM yheight = mcclmon2_yheight;
MCNUM Lmin = mcclmon2_Lmin;
MCNUM Lmax = mcclmon2_Lmax;
MCNUM restore_neutron = mcclmon2_restore_neutron;
int nowritefile = mcclmon2_nowritefile;
#line 84 "/usr/share/mcstas/2.6rc1/monitors/L_monitor.comp"
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
#line 13756 "ISIS_CRISP.c"
}   /* End of lmon2=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbComplmon2:
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

  /* TRACE Component samp1 [14] */
  mccoordschange(mcposrsamp1, mcrotrsamp1,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component samp1 (without coords transformations) */
  mcJumpTrace_samp1:
  SIG_MESSAGE("samp1 (Trace)");
  mcDEBUG_COMP("samp1")
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

#define mcabsorbComp mcabsorbCompsamp1
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
#define mccompcurname  samp1
#define mccompcurtype  Multilayer_Sample
#define mccompcurindex 14
#define xwidth mccsamp1_xwidth
#define zlength mccsamp1_zlength
#define nlayer mccsamp1_nlayer
#define sldPar mccsamp1_sldPar
#define dPar mccsamp1_dPar
#define sigmaPar mccsamp1_sigmaPar
#define frac_inc mccsamp1_frac_inc
#define ythick mccsamp1_ythick
#define mu_inc mccsamp1_mu_inc
#define target_index mccsamp1_target_index
#define focus_xw mccsamp1_focus_xw
#define focus_yh mccsamp1_focus_yh
#define sldParPtr mccsamp1_sldParPtr
#define dParPtr mccsamp1_dParPtr
#define sigmaParPtr mccsamp1_sigmaParPtr
#define tx mccsamp1_tx
#define ty mccsamp1_ty
#define tz mccsamp1_tz
#line 101 "/usr/share/mcstas/2.6rc1/contrib/Multilayer_Sample.comp"
{
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
        tcvar1=gsl_complex_rect(1.0 - (lsq * sldParPtr[i] / tpi),dbl0);
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
		gsl_matrix_complex_set(betan,i+1,0,gsl_complex_mul_real(tcvar1,tpi*dParPtr[i]/lambda));
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
      tcvar4=gsl_complex_mul_real(tcvar3,-1.0*tlc*sigmaParPtr[0]*sigmaParPtr[0]);
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
      		tcvar4=gsl_complex_mul_real(tcvar3,-1.0*tlc*sigmaParPtr[i]*sigmaParPtr[i]);
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
}
#line 14113 "ISIS_CRISP.c"
#undef tz
#undef ty
#undef tx
#undef sigmaParPtr
#undef dParPtr
#undef sldParPtr
#undef focus_yh
#undef focus_xw
#undef target_index
#undef mu_inc
#undef ythick
#undef frac_inc
#undef sigmaPar
#undef dPar
#undef sldPar
#undef nlayer
#undef zlength
#undef xwidth
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompsamp1:
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

  /* TRACE Component slit3 [15] */
  mccoordschange(mcposrslit3, mcrotrslit3,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component slit3 (without coords transformations) */
  mcJumpTrace_slit3:
  SIG_MESSAGE("slit3 (Trace)");
  mcDEBUG_COMP("slit3")
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

#define mcabsorbComp mcabsorbCompslit3
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
#define mccompcurname  slit3
#define mccompcurtype  Slit
#define mccompcurindex 15
{   /* Declarations of slit3=Slit() SETTING parameters. */
MCNUM xmin = mccslit3_xmin;
MCNUM xmax = mccslit3_xmax;
MCNUM ymin = mccslit3_ymin;
MCNUM ymax = mccslit3_ymax;
MCNUM radius = mccslit3_radius;
MCNUM xwidth = mccslit3_xwidth;
MCNUM yheight = mccslit3_yheight;
#line 71 "/usr/share/mcstas/2.6rc1/optics/Slit.comp"
{
    PROP_Z0;
    if (((radius == 0) && (x<xmin || x>xmax || y<ymin || y>ymax))
	|| ((radius != 0) && (x*x + y*y > radius*radius))) {
      RESTORE_NEUTRON(INDEX_CURRENT_COMP, x, y, z, vx, vy, vz, t, sx, sy, sz, p);
      ABSORB;
    }
    else
        SCATTER;
}
#line 14254 "ISIS_CRISP.c"
}   /* End of slit3=Slit() SETTING parameter declarations. */
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompslit3:
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

  /* TRACE Component eguide2 [16] */
  mccoordschange(mcposreguide2, mcrotreguide2,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component eguide2 (without coords transformations) */
  mcJumpTrace_eguide2:
  SIG_MESSAGE("eguide2 (Trace)");
  mcDEBUG_COMP("eguide2")
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

#define mcabsorbComp mcabsorbCompeguide2
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
#define mccompcurname  eguide2
#define mccompcurtype  Guide_tapering
#define mccompcurindex 16
#define w1c mcceguide2_w1c
#define w2c mcceguide2_w2c
#define ww mcceguide2_ww
#define hh mcceguide2_hh
#define whalf mcceguide2_whalf
#define hhalf mcceguide2_hhalf
#define lwhalf mcceguide2_lwhalf
#define lhhalf mcceguide2_lhhalf
#define h1_in mcceguide2_h1_in
#define h2_out mcceguide2_h2_out
#define w1_in mcceguide2_w1_in
#define w2_out mcceguide2_w2_out
#define l_seg mcceguide2_l_seg
#define seg mcceguide2_seg
#define h12 mcceguide2_h12
#define h2 mcceguide2_h2
#define w12 mcceguide2_w12
#define w2 mcceguide2_w2
#define a_ell_q mcceguide2_a_ell_q
#define b_ell_q mcceguide2_b_ell_q
#define lbw mcceguide2_lbw
#define lbh mcceguide2_lbh
#define mxi mcceguide2_mxi
#define u1 mcceguide2_u1
#define u2 mcceguide2_u2
#define div1 mcceguide2_div1
#define p2_para mcceguide2_p2_para
#define test mcceguide2_test
#define Div1 mcceguide2_Div1
#define i mcceguide2_i
#define ii mcceguide2_ii
#define seg mcceguide2_seg
#define fu mcceguide2_fu
#define pos mcceguide2_pos
#define file_name mcceguide2_file_name
#define ep mcceguide2_ep
#define num mcceguide2_num
#define rotation_h mcceguide2_rotation_h
#define rotation_v mcceguide2_rotation_v
{   /* Declarations of eguide2=Guide_tapering() SETTING parameters. */
char* option = mcceguide2_option;
MCNUM w1 = mcceguide2_w1;
MCNUM h1 = mcceguide2_h1;
MCNUM l = mcceguide2_l;
MCNUM linw = mcceguide2_linw;
MCNUM loutw = mcceguide2_loutw;
MCNUM linh = mcceguide2_linh;
MCNUM louth = mcceguide2_louth;
MCNUM R0 = mcceguide2_R0;
MCNUM Qcx = mcceguide2_Qcx;
MCNUM Qcy = mcceguide2_Qcy;
MCNUM alphax = mcceguide2_alphax;
MCNUM alphay = mcceguide2_alphay;
MCNUM W = mcceguide2_W;
MCNUM mx = mcceguide2_mx;
MCNUM my = mcceguide2_my;
MCNUM segno = mcceguide2_segno;
MCNUM curvature = mcceguide2_curvature;
MCNUM curvature_v = mcceguide2_curvature_v;
#line 452 "/usr/share/mcstas/2.6rc1/optics/Guide_tapering.comp"
{
  double t1,t2,ts,zr;                           /* Intersection times. */
  double av,ah,bv,bh,cv1,cv2,ch1,ch2,dd;        /* Intermediate values */
  double vdotn_v1,vdotn_v2,vdotn_h1,vdotn_h2;   /* Dot products. */
  int i;                                        /* Which mirror hit? */
  double q;                                     /* Q [1/AA] of reflection */
  double vlen2,nlen2;                           /* Vector lengths squared */
  double edge;
  double hadj;                                  /* Channel displacement */
  double sz=0;
  int ii;
  Coords mctc1,mctc2;
  Rotation mctr1;

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

}
#line 14573 "ISIS_CRISP.c"
}   /* End of eguide2=Guide_tapering() SETTING parameter declarations. */
#undef rotation_v
#undef rotation_h
#undef num
#undef ep
#undef file_name
#undef pos
#undef fu
#undef seg
#undef ii
#undef i
#undef Div1
#undef test
#undef p2_para
#undef div1
#undef u2
#undef u1
#undef mxi
#undef lbh
#undef lbw
#undef b_ell_q
#undef a_ell_q
#undef w2
#undef w12
#undef h2
#undef h12
#undef seg
#undef l_seg
#undef w2_out
#undef w1_in
#undef h2_out
#undef h1_in
#undef lhhalf
#undef lwhalf
#undef hhalf
#undef whalf
#undef hh
#undef ww
#undef w2c
#undef w1c
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompeguide2:
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

  /* TRACE Component lmon3 [17] */
  mccoordschange(mcposrlmon3, mcrotrlmon3,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component lmon3 (without coords transformations) */
  mcJumpTrace_lmon3:
  SIG_MESSAGE("lmon3 (Trace)");
  mcDEBUG_COMP("lmon3")
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

#define mcabsorbComp mcabsorbComplmon3
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
#define mccompcurname  lmon3
#define mccompcurtype  L_monitor
#define mccompcurindex 17
#define nL mcclmon3_nL
#define L_N mcclmon3_L_N
#define L_p mcclmon3_L_p
#define L_p2 mcclmon3_L_p2
{   /* Declarations of lmon3=L_monitor() SETTING parameters. */
char* filename = mcclmon3_filename;
MCNUM xmin = mcclmon3_xmin;
MCNUM xmax = mcclmon3_xmax;
MCNUM ymin = mcclmon3_ymin;
MCNUM ymax = mcclmon3_ymax;
MCNUM xwidth = mcclmon3_xwidth;
MCNUM yheight = mcclmon3_yheight;
MCNUM Lmin = mcclmon3_Lmin;
MCNUM Lmax = mcclmon3_Lmax;
MCNUM restore_neutron = mcclmon3_restore_neutron;
int nowritefile = mcclmon3_nowritefile;
#line 84 "/usr/share/mcstas/2.6rc1/monitors/L_monitor.comp"
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
#line 14755 "ISIS_CRISP.c"
}   /* End of lmon3=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbComplmon3:
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

  /* TRACE Component PSDdet1 [18] */
  mccoordschange(mcposrPSDdet1, mcrotrPSDdet1,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component PSDdet1 (without coords transformations) */
  mcJumpTrace_PSDdet1:
  SIG_MESSAGE("PSDdet1 (Trace)");
  mcDEBUG_COMP("PSDdet1")
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

#define mcabsorbComp mcabsorbCompPSDdet1
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
#define mccompcurname  PSDdet1
#define mccompcurtype  PSD_monitor
#define mccompcurindex 18
#define PSD_N mccPSDdet1_PSD_N
#define PSD_p mccPSDdet1_PSD_p
#define PSD_p2 mccPSDdet1_PSD_p2
{   /* Declarations of PSDdet1=PSD_monitor() SETTING parameters. */
int nx = mccPSDdet1_nx;
int ny = mccPSDdet1_ny;
char* filename = mccPSDdet1_filename;
MCNUM xmin = mccPSDdet1_xmin;
MCNUM xmax = mccPSDdet1_xmax;
MCNUM ymin = mccPSDdet1_ymin;
MCNUM ymax = mccPSDdet1_ymax;
MCNUM xwidth = mccPSDdet1_xwidth;
MCNUM yheight = mccPSDdet1_yheight;
MCNUM restore_neutron = mccPSDdet1_restore_neutron;
#line 94 "/usr/share/mcstas/2.6rc1/monitors/PSD_monitor.comp"
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
#line 14893 "ISIS_CRISP.c"
}   /* End of PSDdet1=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompPSDdet1:
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

  /* User SAVE code for component 'lmon1'. */
  SIG_MESSAGE("lmon1 (Save)");
#define mccompcurname  lmon1
#define mccompcurtype  L_monitor
#define mccompcurindex 7
#define nL mcclmon1_nL
#define L_N mcclmon1_L_N
#define L_p mcclmon1_L_p
#define L_p2 mcclmon1_L_p2
{   /* Declarations of lmon1=L_monitor() SETTING parameters. */
char* filename = mcclmon1_filename;
MCNUM xmin = mcclmon1_xmin;
MCNUM xmax = mcclmon1_xmax;
MCNUM ymin = mcclmon1_ymin;
MCNUM ymax = mcclmon1_ymax;
MCNUM xwidth = mcclmon1_xwidth;
MCNUM yheight = mcclmon1_yheight;
MCNUM Lmin = mcclmon1_Lmin;
MCNUM Lmax = mcclmon1_Lmax;
MCNUM restore_neutron = mcclmon1_restore_neutron;
int nowritefile = mcclmon1_nowritefile;
#line 107 "/usr/share/mcstas/2.6rc1/monitors/L_monitor.comp"
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
#line 15007 "ISIS_CRISP.c"
}   /* End of lmon1=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'PSDmon1'. */
  SIG_MESSAGE("PSDmon1 (Save)");
#define mccompcurname  PSDmon1
#define mccompcurtype  PSD_monitor
#define mccompcurindex 8
#define PSD_N mccPSDmon1_PSD_N
#define PSD_p mccPSDmon1_PSD_p
#define PSD_p2 mccPSDmon1_PSD_p2
{   /* Declarations of PSDmon1=PSD_monitor() SETTING parameters. */
int nx = mccPSDmon1_nx;
int ny = mccPSDmon1_ny;
char* filename = mccPSDmon1_filename;
MCNUM xmin = mccPSDmon1_xmin;
MCNUM xmax = mccPSDmon1_xmax;
MCNUM ymin = mccPSDmon1_ymin;
MCNUM ymax = mccPSDmon1_ymax;
MCNUM xwidth = mccPSDmon1_xwidth;
MCNUM yheight = mccPSDmon1_yheight;
MCNUM restore_neutron = mccPSDmon1_restore_neutron;
#line 110 "/usr/share/mcstas/2.6rc1/monitors/PSD_monitor.comp"
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
#line 15047 "ISIS_CRISP.c"
}   /* End of PSDmon1=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'PSDmons1'. */
  SIG_MESSAGE("PSDmons1 (Save)");
#define mccompcurname  PSDmons1
#define mccompcurtype  PSD_monitor
#define mccompcurindex 10
#define PSD_N mccPSDmons1_PSD_N
#define PSD_p mccPSDmons1_PSD_p
#define PSD_p2 mccPSDmons1_PSD_p2
{   /* Declarations of PSDmons1=PSD_monitor() SETTING parameters. */
int nx = mccPSDmons1_nx;
int ny = mccPSDmons1_ny;
char* filename = mccPSDmons1_filename;
MCNUM xmin = mccPSDmons1_xmin;
MCNUM xmax = mccPSDmons1_xmax;
MCNUM ymin = mccPSDmons1_ymin;
MCNUM ymax = mccPSDmons1_ymax;
MCNUM xwidth = mccPSDmons1_xwidth;
MCNUM yheight = mccPSDmons1_yheight;
MCNUM restore_neutron = mccPSDmons1_restore_neutron;
#line 110 "/usr/share/mcstas/2.6rc1/monitors/PSD_monitor.comp"
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
#line 15086 "ISIS_CRISP.c"
}   /* End of PSDmons1=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'lmon2'. */
  SIG_MESSAGE("lmon2 (Save)");
#define mccompcurname  lmon2
#define mccompcurtype  L_monitor
#define mccompcurindex 13
#define nL mcclmon2_nL
#define L_N mcclmon2_L_N
#define L_p mcclmon2_L_p
#define L_p2 mcclmon2_L_p2
{   /* Declarations of lmon2=L_monitor() SETTING parameters. */
char* filename = mcclmon2_filename;
MCNUM xmin = mcclmon2_xmin;
MCNUM xmax = mcclmon2_xmax;
MCNUM ymin = mcclmon2_ymin;
MCNUM ymax = mcclmon2_ymax;
MCNUM xwidth = mcclmon2_xwidth;
MCNUM yheight = mcclmon2_yheight;
MCNUM Lmin = mcclmon2_Lmin;
MCNUM Lmax = mcclmon2_Lmax;
MCNUM restore_neutron = mcclmon2_restore_neutron;
int nowritefile = mcclmon2_nowritefile;
#line 107 "/usr/share/mcstas/2.6rc1/monitors/L_monitor.comp"
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
#line 15128 "ISIS_CRISP.c"
}   /* End of lmon2=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'lmon3'. */
  SIG_MESSAGE("lmon3 (Save)");
#define mccompcurname  lmon3
#define mccompcurtype  L_monitor
#define mccompcurindex 17
#define nL mcclmon3_nL
#define L_N mcclmon3_L_N
#define L_p mcclmon3_L_p
#define L_p2 mcclmon3_L_p2
{   /* Declarations of lmon3=L_monitor() SETTING parameters. */
char* filename = mcclmon3_filename;
MCNUM xmin = mcclmon3_xmin;
MCNUM xmax = mcclmon3_xmax;
MCNUM ymin = mcclmon3_ymin;
MCNUM ymax = mcclmon3_ymax;
MCNUM xwidth = mcclmon3_xwidth;
MCNUM yheight = mcclmon3_yheight;
MCNUM Lmin = mcclmon3_Lmin;
MCNUM Lmax = mcclmon3_Lmax;
MCNUM restore_neutron = mcclmon3_restore_neutron;
int nowritefile = mcclmon3_nowritefile;
#line 107 "/usr/share/mcstas/2.6rc1/monitors/L_monitor.comp"
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
#line 15171 "ISIS_CRISP.c"
}   /* End of lmon3=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'PSDdet1'. */
  SIG_MESSAGE("PSDdet1 (Save)");
#define mccompcurname  PSDdet1
#define mccompcurtype  PSD_monitor
#define mccompcurindex 18
#define PSD_N mccPSDdet1_PSD_N
#define PSD_p mccPSDdet1_PSD_p
#define PSD_p2 mccPSDdet1_PSD_p2
{   /* Declarations of PSDdet1=PSD_monitor() SETTING parameters. */
int nx = mccPSDdet1_nx;
int ny = mccPSDdet1_ny;
char* filename = mccPSDdet1_filename;
MCNUM xmin = mccPSDdet1_xmin;
MCNUM xmax = mccPSDdet1_xmax;
MCNUM ymin = mccPSDdet1_ymin;
MCNUM ymax = mccPSDdet1_ymax;
MCNUM xwidth = mccPSDdet1_xwidth;
MCNUM yheight = mccPSDdet1_yheight;
MCNUM restore_neutron = mccPSDdet1_restore_neutron;
#line 110 "/usr/share/mcstas/2.6rc1/monitors/PSD_monitor.comp"
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
#line 15211 "ISIS_CRISP.c"
}   /* End of PSDdet1=PSD_monitor() SETTING parameter declarations. */
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

    if (!mcNCounter[1]) fprintf(stderr, "Warning: No neutron could reach Component[1] a1\n");
    if (mcAbsorbProp[1]) fprintf(stderr, "Warning: %g events were removed in Component[1] a1=Arm()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[1]);
    if (!mcNCounter[2]) fprintf(stderr, "Warning: No neutron could reach Component[2] isis_source\n");
    if (mcAbsorbProp[2]) fprintf(stderr, "Warning: %g events were removed in Component[2] isis_source=ISIS_moderator()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[2]);
    if (!mcNCounter[3]) fprintf(stderr, "Warning: No neutron could reach Component[3] defslit1\n");
    if (mcAbsorbProp[3]) fprintf(stderr, "Warning: %g events were removed in Component[3] defslit1=Slit()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[3]);
    if (!mcNCounter[4]) fprintf(stderr, "Warning: No neutron could reach Component[4] defslit2\n");
    if (mcAbsorbProp[4]) fprintf(stderr, "Warning: %g events were removed in Component[4] defslit2=Slit()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[4]);
    if (!mcNCounter[5]) fprintf(stderr, "Warning: No neutron could reach Component[5] coarseslit1\n");
    if (mcAbsorbProp[5]) fprintf(stderr, "Warning: %g events were removed in Component[5] coarseslit1=Slit()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[5]);
    if (!mcNCounter[6]) fprintf(stderr, "Warning: No neutron could reach Component[6] coarseslit2\n");
    if (mcAbsorbProp[6]) fprintf(stderr, "Warning: %g events were removed in Component[6] coarseslit2=Slit()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[6]);
    if (!mcNCounter[7]) fprintf(stderr, "Warning: No neutron could reach Component[7] lmon1\n");
    if (mcAbsorbProp[7]) fprintf(stderr, "Warning: %g events were removed in Component[7] lmon1=L_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[7]);
  /* User FINALLY code for component 'PSDmon1'. */
  SIG_MESSAGE("PSDmon1 (Finally)");
#define mccompcurname  PSDmon1
#define mccompcurtype  PSD_monitor
#define mccompcurindex 8
#define PSD_N mccPSDmon1_PSD_N
#define PSD_p mccPSDmon1_PSD_p
#define PSD_p2 mccPSDmon1_PSD_p2
{   /* Declarations of PSDmon1=PSD_monitor() SETTING parameters. */
int nx = mccPSDmon1_nx;
int ny = mccPSDmon1_ny;
char* filename = mccPSDmon1_filename;
MCNUM xmin = mccPSDmon1_xmin;
MCNUM xmax = mccPSDmon1_xmax;
MCNUM ymin = mccPSDmon1_ymin;
MCNUM ymax = mccPSDmon1_ymax;
MCNUM xwidth = mccPSDmon1_xwidth;
MCNUM yheight = mccPSDmon1_yheight;
MCNUM restore_neutron = mccPSDmon1_restore_neutron;
#line 122 "/usr/share/mcstas/2.6rc1/monitors/PSD_monitor.comp"
{
  destroy_darr2d(PSD_N);
  destroy_darr2d(PSD_p);
  destroy_darr2d(PSD_p2);
}
#line 15266 "ISIS_CRISP.c"
}   /* End of PSDmon1=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[8]) fprintf(stderr, "Warning: No neutron could reach Component[8] PSDmon1\n");
    if (mcAbsorbProp[8]) fprintf(stderr, "Warning: %g events were removed in Component[8] PSDmon1=PSD_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[8]);
    if (!mcNCounter[9]) fprintf(stderr, "Warning: No neutron could reach Component[9] slit1\n");
    if (mcAbsorbProp[9]) fprintf(stderr, "Warning: %g events were removed in Component[9] slit1=Slit()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[9]);
  /* User FINALLY code for component 'PSDmons1'. */
  SIG_MESSAGE("PSDmons1 (Finally)");
#define mccompcurname  PSDmons1
#define mccompcurtype  PSD_monitor
#define mccompcurindex 10
#define PSD_N mccPSDmons1_PSD_N
#define PSD_p mccPSDmons1_PSD_p
#define PSD_p2 mccPSDmons1_PSD_p2
{   /* Declarations of PSDmons1=PSD_monitor() SETTING parameters. */
int nx = mccPSDmons1_nx;
int ny = mccPSDmons1_ny;
char* filename = mccPSDmons1_filename;
MCNUM xmin = mccPSDmons1_xmin;
MCNUM xmax = mccPSDmons1_xmax;
MCNUM ymin = mccPSDmons1_ymin;
MCNUM ymax = mccPSDmons1_ymax;
MCNUM xwidth = mccPSDmons1_xwidth;
MCNUM yheight = mccPSDmons1_yheight;
MCNUM restore_neutron = mccPSDmons1_restore_neutron;
#line 122 "/usr/share/mcstas/2.6rc1/monitors/PSD_monitor.comp"
{
  destroy_darr2d(PSD_N);
  destroy_darr2d(PSD_p);
  destroy_darr2d(PSD_p2);
}
#line 15304 "ISIS_CRISP.c"
}   /* End of PSDmons1=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[10]) fprintf(stderr, "Warning: No neutron could reach Component[10] PSDmons1\n");
    if (mcAbsorbProp[10]) fprintf(stderr, "Warning: %g events were removed in Component[10] PSDmons1=PSD_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[10]);
  /* User FINALLY code for component 'eguide1'. */
  SIG_MESSAGE("eguide1 (Finally)");
#define mccompcurname  eguide1
#define mccompcurtype  Guide_tapering
#define mccompcurindex 11
#define w1c mcceguide1_w1c
#define w2c mcceguide1_w2c
#define ww mcceguide1_ww
#define hh mcceguide1_hh
#define whalf mcceguide1_whalf
#define hhalf mcceguide1_hhalf
#define lwhalf mcceguide1_lwhalf
#define lhhalf mcceguide1_lhhalf
#define h1_in mcceguide1_h1_in
#define h2_out mcceguide1_h2_out
#define w1_in mcceguide1_w1_in
#define w2_out mcceguide1_w2_out
#define l_seg mcceguide1_l_seg
#define seg mcceguide1_seg
#define h12 mcceguide1_h12
#define h2 mcceguide1_h2
#define w12 mcceguide1_w12
#define w2 mcceguide1_w2
#define a_ell_q mcceguide1_a_ell_q
#define b_ell_q mcceguide1_b_ell_q
#define lbw mcceguide1_lbw
#define lbh mcceguide1_lbh
#define mxi mcceguide1_mxi
#define u1 mcceguide1_u1
#define u2 mcceguide1_u2
#define div1 mcceguide1_div1
#define p2_para mcceguide1_p2_para
#define test mcceguide1_test
#define Div1 mcceguide1_Div1
#define i mcceguide1_i
#define ii mcceguide1_ii
#define seg mcceguide1_seg
#define fu mcceguide1_fu
#define pos mcceguide1_pos
#define file_name mcceguide1_file_name
#define ep mcceguide1_ep
#define num mcceguide1_num
#define rotation_h mcceguide1_rotation_h
#define rotation_v mcceguide1_rotation_v
{   /* Declarations of eguide1=Guide_tapering() SETTING parameters. */
char* option = mcceguide1_option;
MCNUM w1 = mcceguide1_w1;
MCNUM h1 = mcceguide1_h1;
MCNUM l = mcceguide1_l;
MCNUM linw = mcceguide1_linw;
MCNUM loutw = mcceguide1_loutw;
MCNUM linh = mcceguide1_linh;
MCNUM louth = mcceguide1_louth;
MCNUM R0 = mcceguide1_R0;
MCNUM Qcx = mcceguide1_Qcx;
MCNUM Qcy = mcceguide1_Qcy;
MCNUM alphax = mcceguide1_alphax;
MCNUM alphay = mcceguide1_alphay;
MCNUM W = mcceguide1_W;
MCNUM mx = mcceguide1_mx;
MCNUM my = mcceguide1_my;
MCNUM segno = mcceguide1_segno;
MCNUM curvature = mcceguide1_curvature;
MCNUM curvature_v = mcceguide1_curvature_v;
#line 608 "/usr/share/mcstas/2.6rc1/optics/Guide_tapering.comp"
{
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
}
#line 15394 "ISIS_CRISP.c"
}   /* End of eguide1=Guide_tapering() SETTING parameter declarations. */
#undef rotation_v
#undef rotation_h
#undef num
#undef ep
#undef file_name
#undef pos
#undef fu
#undef seg
#undef ii
#undef i
#undef Div1
#undef test
#undef p2_para
#undef div1
#undef u2
#undef u1
#undef mxi
#undef lbh
#undef lbw
#undef b_ell_q
#undef a_ell_q
#undef w2
#undef w12
#undef h2
#undef h12
#undef seg
#undef l_seg
#undef w2_out
#undef w1_in
#undef h2_out
#undef h1_in
#undef lhhalf
#undef lwhalf
#undef hhalf
#undef whalf
#undef hh
#undef ww
#undef w2c
#undef w1c
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[11]) fprintf(stderr, "Warning: No neutron could reach Component[11] eguide1\n");
    if (mcAbsorbProp[11]) fprintf(stderr, "Warning: %g events were removed in Component[11] eguide1=Guide_tapering()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[11]);
    if (!mcNCounter[12]) fprintf(stderr, "Warning: No neutron could reach Component[12] slit2\n");
    if (mcAbsorbProp[12]) fprintf(stderr, "Warning: %g events were removed in Component[12] slit2=Slit()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[12]);
    if (!mcNCounter[13]) fprintf(stderr, "Warning: No neutron could reach Component[13] lmon2\n");
    if (mcAbsorbProp[13]) fprintf(stderr, "Warning: %g events were removed in Component[13] lmon2=L_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[13]);
    if (!mcNCounter[14]) fprintf(stderr, "Warning: No neutron could reach Component[14] samp1\n");
    if (mcAbsorbProp[14]) fprintf(stderr, "Warning: %g events were removed in Component[14] samp1=Multilayer_Sample()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[14]);
    if (!mcNCounter[15]) fprintf(stderr, "Warning: No neutron could reach Component[15] slit3\n");
    if (mcAbsorbProp[15]) fprintf(stderr, "Warning: %g events were removed in Component[15] slit3=Slit()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[15]);
  /* User FINALLY code for component 'eguide2'. */
  SIG_MESSAGE("eguide2 (Finally)");
#define mccompcurname  eguide2
#define mccompcurtype  Guide_tapering
#define mccompcurindex 16
#define w1c mcceguide2_w1c
#define w2c mcceguide2_w2c
#define ww mcceguide2_ww
#define hh mcceguide2_hh
#define whalf mcceguide2_whalf
#define hhalf mcceguide2_hhalf
#define lwhalf mcceguide2_lwhalf
#define lhhalf mcceguide2_lhhalf
#define h1_in mcceguide2_h1_in
#define h2_out mcceguide2_h2_out
#define w1_in mcceguide2_w1_in
#define w2_out mcceguide2_w2_out
#define l_seg mcceguide2_l_seg
#define seg mcceguide2_seg
#define h12 mcceguide2_h12
#define h2 mcceguide2_h2
#define w12 mcceguide2_w12
#define w2 mcceguide2_w2
#define a_ell_q mcceguide2_a_ell_q
#define b_ell_q mcceguide2_b_ell_q
#define lbw mcceguide2_lbw
#define lbh mcceguide2_lbh
#define mxi mcceguide2_mxi
#define u1 mcceguide2_u1
#define u2 mcceguide2_u2
#define div1 mcceguide2_div1
#define p2_para mcceguide2_p2_para
#define test mcceguide2_test
#define Div1 mcceguide2_Div1
#define i mcceguide2_i
#define ii mcceguide2_ii
#define seg mcceguide2_seg
#define fu mcceguide2_fu
#define pos mcceguide2_pos
#define file_name mcceguide2_file_name
#define ep mcceguide2_ep
#define num mcceguide2_num
#define rotation_h mcceguide2_rotation_h
#define rotation_v mcceguide2_rotation_v
{   /* Declarations of eguide2=Guide_tapering() SETTING parameters. */
char* option = mcceguide2_option;
MCNUM w1 = mcceguide2_w1;
MCNUM h1 = mcceguide2_h1;
MCNUM l = mcceguide2_l;
MCNUM linw = mcceguide2_linw;
MCNUM loutw = mcceguide2_loutw;
MCNUM linh = mcceguide2_linh;
MCNUM louth = mcceguide2_louth;
MCNUM R0 = mcceguide2_R0;
MCNUM Qcx = mcceguide2_Qcx;
MCNUM Qcy = mcceguide2_Qcy;
MCNUM alphax = mcceguide2_alphax;
MCNUM alphay = mcceguide2_alphay;
MCNUM W = mcceguide2_W;
MCNUM mx = mcceguide2_mx;
MCNUM my = mcceguide2_my;
MCNUM segno = mcceguide2_segno;
MCNUM curvature = mcceguide2_curvature;
MCNUM curvature_v = mcceguide2_curvature_v;
#line 608 "/usr/share/mcstas/2.6rc1/optics/Guide_tapering.comp"
{
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
}
#line 15528 "ISIS_CRISP.c"
}   /* End of eguide2=Guide_tapering() SETTING parameter declarations. */
#undef rotation_v
#undef rotation_h
#undef num
#undef ep
#undef file_name
#undef pos
#undef fu
#undef seg
#undef ii
#undef i
#undef Div1
#undef test
#undef p2_para
#undef div1
#undef u2
#undef u1
#undef mxi
#undef lbh
#undef lbw
#undef b_ell_q
#undef a_ell_q
#undef w2
#undef w12
#undef h2
#undef h12
#undef seg
#undef l_seg
#undef w2_out
#undef w1_in
#undef h2_out
#undef h1_in
#undef lhhalf
#undef lwhalf
#undef hhalf
#undef whalf
#undef hh
#undef ww
#undef w2c
#undef w1c
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[16]) fprintf(stderr, "Warning: No neutron could reach Component[16] eguide2\n");
    if (mcAbsorbProp[16]) fprintf(stderr, "Warning: %g events were removed in Component[16] eguide2=Guide_tapering()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[16]);
    if (!mcNCounter[17]) fprintf(stderr, "Warning: No neutron could reach Component[17] lmon3\n");
    if (mcAbsorbProp[17]) fprintf(stderr, "Warning: %g events were removed in Component[17] lmon3=L_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[17]);
  /* User FINALLY code for component 'PSDdet1'. */
  SIG_MESSAGE("PSDdet1 (Finally)");
#define mccompcurname  PSDdet1
#define mccompcurtype  PSD_monitor
#define mccompcurindex 18
#define PSD_N mccPSDdet1_PSD_N
#define PSD_p mccPSDdet1_PSD_p
#define PSD_p2 mccPSDdet1_PSD_p2
{   /* Declarations of PSDdet1=PSD_monitor() SETTING parameters. */
int nx = mccPSDdet1_nx;
int ny = mccPSDdet1_ny;
char* filename = mccPSDdet1_filename;
MCNUM xmin = mccPSDdet1_xmin;
MCNUM xmax = mccPSDdet1_xmax;
MCNUM ymin = mccPSDdet1_ymin;
MCNUM ymax = mccPSDdet1_ymax;
MCNUM xwidth = mccPSDdet1_xwidth;
MCNUM yheight = mccPSDdet1_yheight;
MCNUM restore_neutron = mccPSDdet1_restore_neutron;
#line 122 "/usr/share/mcstas/2.6rc1/monitors/PSD_monitor.comp"
{
  destroy_darr2d(PSD_N);
  destroy_darr2d(PSD_p);
  destroy_darr2d(PSD_p2);
}
#line 15602 "ISIS_CRISP.c"
}   /* End of PSDdet1=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[18]) fprintf(stderr, "Warning: No neutron could reach Component[18] PSDdet1\n");
    if (mcAbsorbProp[18]) fprintf(stderr, "Warning: %g events were removed in Component[18] PSDdet1=PSD_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[18]);
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
#define mccompcurtype  Arm
#define mccompcurindex 1
#line 40 "/usr/share/mcstas/2.6rc1/optics/Arm.comp"
{
  /* A bit ugly; hard-coded dimensions. */
  
  line(0,0,0,0.2,0,0);
  line(0,0,0,0,0.2,0);
  line(0,0,0,0,0,0.2);
}
#line 15642 "ISIS_CRISP.c"
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'isis_source'. */
  SIG_MESSAGE("isis_source (McDisplay)");
  printf("MCDISPLAY: component %s\n", "isis_source");
#define mccompcurname  isis_source
#define mccompcurtype  ISIS_moderator
#define mccompcurindex 2
#define Face mccisis_source_Face
#define p_in mccisis_source_p_in
#define Tnpts mccisis_source_Tnpts
#define scaleSize mccisis_source_scaleSize
#define angleArea mccisis_source_angleArea
#define Nsim mccisis_source_Nsim
#define Ncount mccisis_source_Ncount
#define TS mccisis_source_TS
#define rtE0 mccisis_source_rtE0
#define rtE1 mccisis_source_rtE1
#define rtmodX mccisis_source_rtmodX
#define rtmodY mccisis_source_rtmodY
#define TargetStation mccisis_source_TargetStation
#define CurrentWeight mccisis_source_CurrentWeight
{   /* Declarations of isis_source=ISIS_moderator() SETTING parameters. */
MCNUM Emin = mccisis_source_Emin;
MCNUM Emax = mccisis_source_Emax;
MCNUM dist = mccisis_source_dist;
MCNUM focus_xw = mccisis_source_focus_xw;
MCNUM focus_yh = mccisis_source_focus_yh;
MCNUM xwidth = mccisis_source_xwidth;
MCNUM yheight = mccisis_source_yheight;
MCNUM CAngle = mccisis_source_CAngle;
MCNUM SAC = mccisis_source_SAC;
MCNUM Lmin = mccisis_source_Lmin;
MCNUM Lmax = mccisis_source_Lmax;
int target_index = mccisis_source_target_index;
MCNUM verbose = mccisis_source_verbose;
#line 1119 "/usr/share/mcstas/2.6rc1/contrib/ISIS_moderator.comp"
{
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

}
#line 15708 "ISIS_CRISP.c"
}   /* End of isis_source=ISIS_moderator() SETTING parameter declarations. */
#undef CurrentWeight
#undef TargetStation
#undef rtmodY
#undef rtmodX
#undef rtE1
#undef rtE0
#undef TS
#undef Ncount
#undef Nsim
#undef angleArea
#undef scaleSize
#undef Tnpts
#undef p_in
#undef Face
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'defslit1'. */
  SIG_MESSAGE("defslit1 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "defslit1");
#define mccompcurname  defslit1
#define mccompcurtype  Slit
#define mccompcurindex 3
{   /* Declarations of defslit1=Slit() SETTING parameters. */
MCNUM xmin = mccdefslit1_xmin;
MCNUM xmax = mccdefslit1_xmax;
MCNUM ymin = mccdefslit1_ymin;
MCNUM ymax = mccdefslit1_ymax;
MCNUM radius = mccdefslit1_radius;
MCNUM xwidth = mccdefslit1_xwidth;
MCNUM yheight = mccdefslit1_yheight;
#line 83 "/usr/share/mcstas/2.6rc1/optics/Slit.comp"
{
  
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
}
#line 15765 "ISIS_CRISP.c"
}   /* End of defslit1=Slit() SETTING parameter declarations. */
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'defslit2'. */
  SIG_MESSAGE("defslit2 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "defslit2");
#define mccompcurname  defslit2
#define mccompcurtype  Slit
#define mccompcurindex 4
{   /* Declarations of defslit2=Slit() SETTING parameters. */
MCNUM xmin = mccdefslit2_xmin;
MCNUM xmax = mccdefslit2_xmax;
MCNUM ymin = mccdefslit2_ymin;
MCNUM ymax = mccdefslit2_ymax;
MCNUM radius = mccdefslit2_radius;
MCNUM xwidth = mccdefslit2_xwidth;
MCNUM yheight = mccdefslit2_yheight;
#line 83 "/usr/share/mcstas/2.6rc1/optics/Slit.comp"
{
  
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
}
#line 15808 "ISIS_CRISP.c"
}   /* End of defslit2=Slit() SETTING parameter declarations. */
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'coarseslit1'. */
  SIG_MESSAGE("coarseslit1 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "coarseslit1");
#define mccompcurname  coarseslit1
#define mccompcurtype  Slit
#define mccompcurindex 5
{   /* Declarations of coarseslit1=Slit() SETTING parameters. */
MCNUM xmin = mcccoarseslit1_xmin;
MCNUM xmax = mcccoarseslit1_xmax;
MCNUM ymin = mcccoarseslit1_ymin;
MCNUM ymax = mcccoarseslit1_ymax;
MCNUM radius = mcccoarseslit1_radius;
MCNUM xwidth = mcccoarseslit1_xwidth;
MCNUM yheight = mcccoarseslit1_yheight;
#line 83 "/usr/share/mcstas/2.6rc1/optics/Slit.comp"
{
  
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
}
#line 15851 "ISIS_CRISP.c"
}   /* End of coarseslit1=Slit() SETTING parameter declarations. */
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'coarseslit2'. */
  SIG_MESSAGE("coarseslit2 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "coarseslit2");
#define mccompcurname  coarseslit2
#define mccompcurtype  Slit
#define mccompcurindex 6
{   /* Declarations of coarseslit2=Slit() SETTING parameters. */
MCNUM xmin = mcccoarseslit2_xmin;
MCNUM xmax = mcccoarseslit2_xmax;
MCNUM ymin = mcccoarseslit2_ymin;
MCNUM ymax = mcccoarseslit2_ymax;
MCNUM radius = mcccoarseslit2_radius;
MCNUM xwidth = mcccoarseslit2_xwidth;
MCNUM yheight = mcccoarseslit2_yheight;
#line 83 "/usr/share/mcstas/2.6rc1/optics/Slit.comp"
{
  
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
}
#line 15894 "ISIS_CRISP.c"
}   /* End of coarseslit2=Slit() SETTING parameter declarations. */
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'lmon1'. */
  SIG_MESSAGE("lmon1 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "lmon1");
#define mccompcurname  lmon1
#define mccompcurtype  L_monitor
#define mccompcurindex 7
#define nL mcclmon1_nL
#define L_N mcclmon1_L_N
#define L_p mcclmon1_L_p
#define L_p2 mcclmon1_L_p2
{   /* Declarations of lmon1=L_monitor() SETTING parameters. */
char* filename = mcclmon1_filename;
MCNUM xmin = mcclmon1_xmin;
MCNUM xmax = mcclmon1_xmax;
MCNUM ymin = mcclmon1_ymin;
MCNUM ymax = mcclmon1_ymax;
MCNUM xwidth = mcclmon1_xwidth;
MCNUM yheight = mcclmon1_yheight;
MCNUM Lmin = mcclmon1_Lmin;
MCNUM Lmax = mcclmon1_Lmax;
MCNUM restore_neutron = mcclmon1_restore_neutron;
int nowritefile = mcclmon1_nowritefile;
#line 120 "/usr/share/mcstas/2.6rc1/monitors/L_monitor.comp"
{
  
  multiline(5, (double)xmin, (double)ymin, 0.0,
               (double)xmax, (double)ymin, 0.0,
               (double)xmax, (double)ymax, 0.0,
               (double)xmin, (double)ymax, 0.0,
               (double)xmin, (double)ymin, 0.0);
}
#line 15931 "ISIS_CRISP.c"
}   /* End of lmon1=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'PSDmon1'. */
  SIG_MESSAGE("PSDmon1 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "PSDmon1");
#define mccompcurname  PSDmon1
#define mccompcurtype  PSD_monitor
#define mccompcurindex 8
#define PSD_N mccPSDmon1_PSD_N
#define PSD_p mccPSDmon1_PSD_p
#define PSD_p2 mccPSDmon1_PSD_p2
{   /* Declarations of PSDmon1=PSD_monitor() SETTING parameters. */
int nx = mccPSDmon1_nx;
int ny = mccPSDmon1_ny;
char* filename = mccPSDmon1_filename;
MCNUM xmin = mccPSDmon1_xmin;
MCNUM xmax = mccPSDmon1_xmax;
MCNUM ymin = mccPSDmon1_ymin;
MCNUM ymax = mccPSDmon1_ymax;
MCNUM xwidth = mccPSDmon1_xwidth;
MCNUM yheight = mccPSDmon1_yheight;
MCNUM restore_neutron = mccPSDmon1_restore_neutron;
#line 129 "/usr/share/mcstas/2.6rc1/monitors/PSD_monitor.comp"
{
  multiline(5,
    (double)xmin, (double)ymin, 0.0,
    (double)xmax, (double)ymin, 0.0,
    (double)xmax, (double)ymax, 0.0,
    (double)xmin, (double)ymax, 0.0,
    (double)xmin, (double)ymin, 0.0);
}
#line 15970 "ISIS_CRISP.c"
}   /* End of PSDmon1=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'slit1'. */
  SIG_MESSAGE("slit1 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "slit1");
#define mccompcurname  slit1
#define mccompcurtype  Slit
#define mccompcurindex 9
{   /* Declarations of slit1=Slit() SETTING parameters. */
MCNUM xmin = mccslit1_xmin;
MCNUM xmax = mccslit1_xmax;
MCNUM ymin = mccslit1_ymin;
MCNUM ymax = mccslit1_ymax;
MCNUM radius = mccslit1_radius;
MCNUM xwidth = mccslit1_xwidth;
MCNUM yheight = mccslit1_yheight;
#line 83 "/usr/share/mcstas/2.6rc1/optics/Slit.comp"
{
  
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
}
#line 16016 "ISIS_CRISP.c"
}   /* End of slit1=Slit() SETTING parameter declarations. */
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'PSDmons1'. */
  SIG_MESSAGE("PSDmons1 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "PSDmons1");
#define mccompcurname  PSDmons1
#define mccompcurtype  PSD_monitor
#define mccompcurindex 10
#define PSD_N mccPSDmons1_PSD_N
#define PSD_p mccPSDmons1_PSD_p
#define PSD_p2 mccPSDmons1_PSD_p2
{   /* Declarations of PSDmons1=PSD_monitor() SETTING parameters. */
int nx = mccPSDmons1_nx;
int ny = mccPSDmons1_ny;
char* filename = mccPSDmons1_filename;
MCNUM xmin = mccPSDmons1_xmin;
MCNUM xmax = mccPSDmons1_xmax;
MCNUM ymin = mccPSDmons1_ymin;
MCNUM ymax = mccPSDmons1_ymax;
MCNUM xwidth = mccPSDmons1_xwidth;
MCNUM yheight = mccPSDmons1_yheight;
MCNUM restore_neutron = mccPSDmons1_restore_neutron;
#line 129 "/usr/share/mcstas/2.6rc1/monitors/PSD_monitor.comp"
{
  multiline(5,
    (double)xmin, (double)ymin, 0.0,
    (double)xmax, (double)ymin, 0.0,
    (double)xmax, (double)ymax, 0.0,
    (double)xmin, (double)ymax, 0.0,
    (double)xmin, (double)ymin, 0.0);
}
#line 16051 "ISIS_CRISP.c"
}   /* End of PSDmons1=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'eguide1'. */
  SIG_MESSAGE("eguide1 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "eguide1");
#define mccompcurname  eguide1
#define mccompcurtype  Guide_tapering
#define mccompcurindex 11
#define w1c mcceguide1_w1c
#define w2c mcceguide1_w2c
#define ww mcceguide1_ww
#define hh mcceguide1_hh
#define whalf mcceguide1_whalf
#define hhalf mcceguide1_hhalf
#define lwhalf mcceguide1_lwhalf
#define lhhalf mcceguide1_lhhalf
#define h1_in mcceguide1_h1_in
#define h2_out mcceguide1_h2_out
#define w1_in mcceguide1_w1_in
#define w2_out mcceguide1_w2_out
#define l_seg mcceguide1_l_seg
#define seg mcceguide1_seg
#define h12 mcceguide1_h12
#define h2 mcceguide1_h2
#define w12 mcceguide1_w12
#define w2 mcceguide1_w2
#define a_ell_q mcceguide1_a_ell_q
#define b_ell_q mcceguide1_b_ell_q
#define lbw mcceguide1_lbw
#define lbh mcceguide1_lbh
#define mxi mcceguide1_mxi
#define u1 mcceguide1_u1
#define u2 mcceguide1_u2
#define div1 mcceguide1_div1
#define p2_para mcceguide1_p2_para
#define test mcceguide1_test
#define Div1 mcceguide1_Div1
#define i mcceguide1_i
#define ii mcceguide1_ii
#define seg mcceguide1_seg
#define fu mcceguide1_fu
#define pos mcceguide1_pos
#define file_name mcceguide1_file_name
#define ep mcceguide1_ep
#define num mcceguide1_num
#define rotation_h mcceguide1_rotation_h
#define rotation_v mcceguide1_rotation_v
{   /* Declarations of eguide1=Guide_tapering() SETTING parameters. */
char* option = mcceguide1_option;
MCNUM w1 = mcceguide1_w1;
MCNUM h1 = mcceguide1_h1;
MCNUM l = mcceguide1_l;
MCNUM linw = mcceguide1_linw;
MCNUM loutw = mcceguide1_loutw;
MCNUM linh = mcceguide1_linh;
MCNUM louth = mcceguide1_louth;
MCNUM R0 = mcceguide1_R0;
MCNUM Qcx = mcceguide1_Qcx;
MCNUM Qcy = mcceguide1_Qcy;
MCNUM alphax = mcceguide1_alphax;
MCNUM alphay = mcceguide1_alphay;
MCNUM W = mcceguide1_W;
MCNUM mx = mcceguide1_mx;
MCNUM my = mcceguide1_my;
MCNUM segno = mcceguide1_segno;
MCNUM curvature = mcceguide1_curvature;
MCNUM curvature_v = mcceguide1_curvature_v;
#line 624 "/usr/share/mcstas/2.6rc1/optics/Guide_tapering.comp"
{
  double x;
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

}
#line 16158 "ISIS_CRISP.c"
}   /* End of eguide1=Guide_tapering() SETTING parameter declarations. */
#undef rotation_v
#undef rotation_h
#undef num
#undef ep
#undef file_name
#undef pos
#undef fu
#undef seg
#undef ii
#undef i
#undef Div1
#undef test
#undef p2_para
#undef div1
#undef u2
#undef u1
#undef mxi
#undef lbh
#undef lbw
#undef b_ell_q
#undef a_ell_q
#undef w2
#undef w12
#undef h2
#undef h12
#undef seg
#undef l_seg
#undef w2_out
#undef w1_in
#undef h2_out
#undef h1_in
#undef lhhalf
#undef lwhalf
#undef hhalf
#undef whalf
#undef hh
#undef ww
#undef w2c
#undef w1c
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'slit2'. */
  SIG_MESSAGE("slit2 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "slit2");
#define mccompcurname  slit2
#define mccompcurtype  Slit
#define mccompcurindex 12
{   /* Declarations of slit2=Slit() SETTING parameters. */
MCNUM xmin = mccslit2_xmin;
MCNUM xmax = mccslit2_xmax;
MCNUM ymin = mccslit2_ymin;
MCNUM ymax = mccslit2_ymax;
MCNUM radius = mccslit2_radius;
MCNUM xwidth = mccslit2_xwidth;
MCNUM yheight = mccslit2_yheight;
#line 83 "/usr/share/mcstas/2.6rc1/optics/Slit.comp"
{
  
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
}
#line 16240 "ISIS_CRISP.c"
}   /* End of slit2=Slit() SETTING parameter declarations. */
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'lmon2'. */
  SIG_MESSAGE("lmon2 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "lmon2");
#define mccompcurname  lmon2
#define mccompcurtype  L_monitor
#define mccompcurindex 13
#define nL mcclmon2_nL
#define L_N mcclmon2_L_N
#define L_p mcclmon2_L_p
#define L_p2 mcclmon2_L_p2
{   /* Declarations of lmon2=L_monitor() SETTING parameters. */
char* filename = mcclmon2_filename;
MCNUM xmin = mcclmon2_xmin;
MCNUM xmax = mcclmon2_xmax;
MCNUM ymin = mcclmon2_ymin;
MCNUM ymax = mcclmon2_ymax;
MCNUM xwidth = mcclmon2_xwidth;
MCNUM yheight = mcclmon2_yheight;
MCNUM Lmin = mcclmon2_Lmin;
MCNUM Lmax = mcclmon2_Lmax;
MCNUM restore_neutron = mcclmon2_restore_neutron;
int nowritefile = mcclmon2_nowritefile;
#line 120 "/usr/share/mcstas/2.6rc1/monitors/L_monitor.comp"
{
  
  multiline(5, (double)xmin, (double)ymin, 0.0,
               (double)xmax, (double)ymin, 0.0,
               (double)xmax, (double)ymax, 0.0,
               (double)xmin, (double)ymax, 0.0,
               (double)xmin, (double)ymin, 0.0);
}
#line 16277 "ISIS_CRISP.c"
}   /* End of lmon2=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'samp1'. */
  SIG_MESSAGE("samp1 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "samp1");
#define mccompcurname  samp1
#define mccompcurtype  Multilayer_Sample
#define mccompcurindex 14
#define xwidth mccsamp1_xwidth
#define zlength mccsamp1_zlength
#define nlayer mccsamp1_nlayer
#define sldPar mccsamp1_sldPar
#define dPar mccsamp1_dPar
#define sigmaPar mccsamp1_sigmaPar
#define frac_inc mccsamp1_frac_inc
#define ythick mccsamp1_ythick
#define mu_inc mccsamp1_mu_inc
#define target_index mccsamp1_target_index
#define focus_xw mccsamp1_focus_xw
#define focus_yh mccsamp1_focus_yh
#define sldParPtr mccsamp1_sldParPtr
#define dParPtr mccsamp1_dParPtr
#define sigmaParPtr mccsamp1_sigmaParPtr
#define tx mccsamp1_tx
#define ty mccsamp1_ty
#define tz mccsamp1_tz
#line 332 "/usr/share/mcstas/2.6rc1/contrib/Multilayer_Sample.comp"
{
  box(0.0, (double)-ythick/2.0, 0.0, (double)xwidth, (double)ythick, (double)zlength);
}
#line 16315 "ISIS_CRISP.c"
#undef tz
#undef ty
#undef tx
#undef sigmaParPtr
#undef dParPtr
#undef sldParPtr
#undef focus_yh
#undef focus_xw
#undef target_index
#undef mu_inc
#undef ythick
#undef frac_inc
#undef sigmaPar
#undef dPar
#undef sldPar
#undef nlayer
#undef zlength
#undef xwidth
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'slit3'. */
  SIG_MESSAGE("slit3 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "slit3");
#define mccompcurname  slit3
#define mccompcurtype  Slit
#define mccompcurindex 15
{   /* Declarations of slit3=Slit() SETTING parameters. */
MCNUM xmin = mccslit3_xmin;
MCNUM xmax = mccslit3_xmax;
MCNUM ymin = mccslit3_ymin;
MCNUM ymax = mccslit3_ymax;
MCNUM radius = mccslit3_radius;
MCNUM xwidth = mccslit3_xwidth;
MCNUM yheight = mccslit3_yheight;
#line 83 "/usr/share/mcstas/2.6rc1/optics/Slit.comp"
{
  
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
}
#line 16375 "ISIS_CRISP.c"
}   /* End of slit3=Slit() SETTING parameter declarations. */
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'eguide2'. */
  SIG_MESSAGE("eguide2 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "eguide2");
#define mccompcurname  eguide2
#define mccompcurtype  Guide_tapering
#define mccompcurindex 16
#define w1c mcceguide2_w1c
#define w2c mcceguide2_w2c
#define ww mcceguide2_ww
#define hh mcceguide2_hh
#define whalf mcceguide2_whalf
#define hhalf mcceguide2_hhalf
#define lwhalf mcceguide2_lwhalf
#define lhhalf mcceguide2_lhhalf
#define h1_in mcceguide2_h1_in
#define h2_out mcceguide2_h2_out
#define w1_in mcceguide2_w1_in
#define w2_out mcceguide2_w2_out
#define l_seg mcceguide2_l_seg
#define seg mcceguide2_seg
#define h12 mcceguide2_h12
#define h2 mcceguide2_h2
#define w12 mcceguide2_w12
#define w2 mcceguide2_w2
#define a_ell_q mcceguide2_a_ell_q
#define b_ell_q mcceguide2_b_ell_q
#define lbw mcceguide2_lbw
#define lbh mcceguide2_lbh
#define mxi mcceguide2_mxi
#define u1 mcceguide2_u1
#define u2 mcceguide2_u2
#define div1 mcceguide2_div1
#define p2_para mcceguide2_p2_para
#define test mcceguide2_test
#define Div1 mcceguide2_Div1
#define i mcceguide2_i
#define ii mcceguide2_ii
#define seg mcceguide2_seg
#define fu mcceguide2_fu
#define pos mcceguide2_pos
#define file_name mcceguide2_file_name
#define ep mcceguide2_ep
#define num mcceguide2_num
#define rotation_h mcceguide2_rotation_h
#define rotation_v mcceguide2_rotation_v
{   /* Declarations of eguide2=Guide_tapering() SETTING parameters. */
char* option = mcceguide2_option;
MCNUM w1 = mcceguide2_w1;
MCNUM h1 = mcceguide2_h1;
MCNUM l = mcceguide2_l;
MCNUM linw = mcceguide2_linw;
MCNUM loutw = mcceguide2_loutw;
MCNUM linh = mcceguide2_linh;
MCNUM louth = mcceguide2_louth;
MCNUM R0 = mcceguide2_R0;
MCNUM Qcx = mcceguide2_Qcx;
MCNUM Qcy = mcceguide2_Qcy;
MCNUM alphax = mcceguide2_alphax;
MCNUM alphay = mcceguide2_alphay;
MCNUM W = mcceguide2_W;
MCNUM mx = mcceguide2_mx;
MCNUM my = mcceguide2_my;
MCNUM segno = mcceguide2_segno;
MCNUM curvature = mcceguide2_curvature;
MCNUM curvature_v = mcceguide2_curvature_v;
#line 624 "/usr/share/mcstas/2.6rc1/optics/Guide_tapering.comp"
{
  double x;
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

}
#line 16479 "ISIS_CRISP.c"
}   /* End of eguide2=Guide_tapering() SETTING parameter declarations. */
#undef rotation_v
#undef rotation_h
#undef num
#undef ep
#undef file_name
#undef pos
#undef fu
#undef seg
#undef ii
#undef i
#undef Div1
#undef test
#undef p2_para
#undef div1
#undef u2
#undef u1
#undef mxi
#undef lbh
#undef lbw
#undef b_ell_q
#undef a_ell_q
#undef w2
#undef w12
#undef h2
#undef h12
#undef seg
#undef l_seg
#undef w2_out
#undef w1_in
#undef h2_out
#undef h1_in
#undef lhhalf
#undef lwhalf
#undef hhalf
#undef whalf
#undef hh
#undef ww
#undef w2c
#undef w1c
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'lmon3'. */
  SIG_MESSAGE("lmon3 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "lmon3");
#define mccompcurname  lmon3
#define mccompcurtype  L_monitor
#define mccompcurindex 17
#define nL mcclmon3_nL
#define L_N mcclmon3_L_N
#define L_p mcclmon3_L_p
#define L_p2 mcclmon3_L_p2
{   /* Declarations of lmon3=L_monitor() SETTING parameters. */
char* filename = mcclmon3_filename;
MCNUM xmin = mcclmon3_xmin;
MCNUM xmax = mcclmon3_xmax;
MCNUM ymin = mcclmon3_ymin;
MCNUM ymax = mcclmon3_ymax;
MCNUM xwidth = mcclmon3_xwidth;
MCNUM yheight = mcclmon3_yheight;
MCNUM Lmin = mcclmon3_Lmin;
MCNUM Lmax = mcclmon3_Lmax;
MCNUM restore_neutron = mcclmon3_restore_neutron;
int nowritefile = mcclmon3_nowritefile;
#line 120 "/usr/share/mcstas/2.6rc1/monitors/L_monitor.comp"
{
  
  multiline(5, (double)xmin, (double)ymin, 0.0,
               (double)xmax, (double)ymin, 0.0,
               (double)xmax, (double)ymax, 0.0,
               (double)xmin, (double)ymax, 0.0,
               (double)xmin, (double)ymin, 0.0);
}
#line 16555 "ISIS_CRISP.c"
}   /* End of lmon3=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'PSDdet1'. */
  SIG_MESSAGE("PSDdet1 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "PSDdet1");
#define mccompcurname  PSDdet1
#define mccompcurtype  PSD_monitor
#define mccompcurindex 18
#define PSD_N mccPSDdet1_PSD_N
#define PSD_p mccPSDdet1_PSD_p
#define PSD_p2 mccPSDdet1_PSD_p2
{   /* Declarations of PSDdet1=PSD_monitor() SETTING parameters. */
int nx = mccPSDdet1_nx;
int ny = mccPSDdet1_ny;
char* filename = mccPSDdet1_filename;
MCNUM xmin = mccPSDdet1_xmin;
MCNUM xmax = mccPSDdet1_xmax;
MCNUM ymin = mccPSDdet1_ymin;
MCNUM ymax = mccPSDdet1_ymax;
MCNUM xwidth = mccPSDdet1_xwidth;
MCNUM yheight = mccPSDdet1_yheight;
MCNUM restore_neutron = mccPSDdet1_restore_neutron;
#line 129 "/usr/share/mcstas/2.6rc1/monitors/PSD_monitor.comp"
{
  multiline(5,
    (double)xmin, (double)ymin, 0.0,
    (double)xmax, (double)ymin, 0.0,
    (double)xmax, (double)ymax, 0.0,
    (double)xmin, (double)ymax, 0.0,
    (double)xmin, (double)ymin, 0.0);
}
#line 16594 "ISIS_CRISP.c"
}   /* End of PSDdet1=PSD_monitor() SETTING parameter declarations. */
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
/* end of generated C code ISIS_CRISP.c */
