/* Automatically generated file. Do not edit. 
 * Format:     ANSI C source code
 * Creator:    McStas <http://www.mcstas.org>
 * Instrument: /zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr (HZB_FLEX)
 * Date:       Wed Nov 20 00:13:52 2019
 * File:       ./HZB_FLEX.c
 * Compile:    cc -o HZB_FLEX.out ./HZB_FLEX.c 
 * CFLAGS=
 */


#define MCCODE_STRING "McStas 2.5 - Nov. 19, 2019"
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
#define MCCODE_STRING "McStas 2.5 - Nov. 19, 2019"
#endif

#ifndef MCCODE_DATE
#define MCCODE_DATE "Nov. 19, 2019"
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

#line 712 "./HZB_FLEX.c"

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

#line 945 "./HZB_FLEX.c"

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

#line 4977 "./HZB_FLEX.c"

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

#line 5337 "./HZB_FLEX.c"
#ifdef MC_TRACE_ENABLED
int mctraceenabled = 1;
#else
int mctraceenabled = 0;
#endif
#define MCSTAS "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../"
int mcdefaultmain = 1;
char mcinstrument_name[] = "HZB_FLEX";
char mcinstrument_source[] = "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr";
char *mcinstrument_exe=NULL; /* will be set to argv[0] in main */
int main(int argc, char *argv[]){return mccode_main(argc, argv);}
void mcinit(void);
void mcraytrace(void);
void mcsave(FILE *);
void mcfinally(void);
void mcdisplay(void);

/* Shared user declarations for all components 'Source_gen'. */
#line 140 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../sources/Source_gen.comp"
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

#line 6839 "./HZB_FLEX.c"

/* Shared user declarations for all components 'Guide'. */
#line 63 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide.comp"

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

#line 6992 "./HZB_FLEX.c"

/* Shared user declarations for all components 'Guide_curved'. */
#line 55 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../contrib/Guide_curved.comp"

#line 6997 "./HZB_FLEX.c"

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

#line 7083 "./HZB_FLEX.c"

/* Shared user declarations for all components 'Guide_tapering'. */
#line 91 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide_tapering.comp"

#line 7088 "./HZB_FLEX.c"

/* Shared user declarations for all components 'Guide_channeled'. */
#line 76 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide_channeled.comp"

#line 7093 "./HZB_FLEX.c"

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


#line 7119 "./HZB_FLEX.c"

/* Instrument parameters. */
MCNUM mcipkI;
MCNUM mcipkF;
MCNUM mcipwVS;
MCNUM mciptilt;
MCNUM mcipSA;
MCNUM mcipA3;
MCNUM mcipA4;
MCNUM mcipL3;
MCNUM mcipL4;
int mcipMono_flatswitch;

#define mcNUMIPAR 10
int mcnumipar = 10;
struct mcinputtable_struct mcinputtable[mcNUMIPAR+1] = {
  "kI", &mcipkI, instr_type_double, "1.55", 
  "kF", &mcipkF, instr_type_double, "1.55", 
  "wVS", &mcipwVS, instr_type_double, "0.03", 
  "tilt", &mciptilt, instr_type_double, "0", 
  "SA", &mcipSA, instr_type_double, "-1", 
  "A3", &mcipA3, instr_type_double, "0", 
  "A4", &mcipA4, instr_type_double, "70", 
  "L3", &mcipL3, instr_type_double, "1.00", 
  "L4", &mcipL4, instr_type_double, "1.00", 
  "Mono_flatswitch", &mcipMono_flatswitch, instr_type_int, "0", 
  NULL, NULL, instr_type_double, ""
};

/* User declarations from instrument definition. */
#define mccompcurname  HZB_FLEX
#define mccompcurtype  INSTRUMENT
#define mccompcurindex 0
#define mcposaHZB_FLEX coords_set(0,0,0)
#define kI mcipkI
#define kF mcipkF
#define wVS mcipwVS
#define tilt mciptilt
#define SA mcipSA
#define A3 mcipA3
#define A4 mcipA4
#define L3 mcipL3
#define L4 mcipL4
#define Mono_flatswitch mcipMono_flatswitch
#line 37 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
#include <math.h>

// unit conversions
double m_n=1.674928e-27;     // in kg
double hbar=1.054572e-34;    // in J sec
double meV_to_SI=1.6022e-22; // conversion
double SI_to_meV=6.2414e21;

// guide coating properties
double MGUIDE=3.2;
double W_para = 0.0025;
double R0_para = 0.99;
double Qc_para = 0.0217;
double alpha_para=3.90;

// curved guide properties
double LGUIDE_1=18.65;
double LGUIDE_2=17.8;
double beta_guide_1, s_guide_1;
double beta_guide_2, s_guide_2;

// monochromator properties
double RMH, RMV, rho_MH, rho_MV;
double lambdaI, EI;

// analyser properties
// double RAH, RAV, rho_AH, rho_AV;
// rho_AV=1.00 for kf=0.952
// rho_AV=1.23 for kf=1.55
// rho_AV=1.40 for kf=2.62
// rho_AV=1.97 for kf=3.696
double RAH, rho_AV=1.23, RAV;
double lambdaF, EF;

// angles
double A1, A2, A5, A6;

// distances
double L1=1.75;
double L2=2.2; // This is an actual variable in the real experiment
// double L3=1.10;
// double L4=1.00;


//velocity selector frequency
double SelFreq;

// monitor/source properties
double dlambdaI,dkI,dEI,l_min,l_max,e_min,e_max;

#line 7215 "./HZB_FLEX.c"
#undef Mono_flatswitch
#undef L4
#undef L3
#undef A4
#undef A3
#undef SA
#undef tilt
#undef wVS
#undef kF
#undef kI
#undef mcposaHZB_FLEX
#undef mccompcurindex
#undef mccompcurtype
#undef mccompcurname

/* neutron state table at each component input (local coords) */
/* [x, y, z, vx, vy, vz, t, sx, sy, sz, p] */
MCNUM mccomp_storein[11*36];
/* Components position table (absolute and relative coords) */
Coords mccomp_posa[36];
Coords mccomp_posr[36];
/* Counter for each comp to check for inactive ones */
MCNUM  mcNCounter[36];
MCNUM  mcPCounter[36];
MCNUM  mcP2Counter[36];
#define mcNUMCOMP 35 /* number of components */
/* Counter for PROP ABSORB */
MCNUM  mcAbsorbProp[36];
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

/* Setting parameters for component 'Gen_Source' [2]. */
char mccGen_Source_flux_file[16384];
char mccGen_Source_xdiv_file[16384];
char mccGen_Source_ydiv_file[16384];
MCNUM mccGen_Source_radius;
MCNUM mccGen_Source_dist;
MCNUM mccGen_Source_focus_xw;
MCNUM mccGen_Source_focus_yh;
MCNUM mccGen_Source_focus_aw;
MCNUM mccGen_Source_focus_ah;
MCNUM mccGen_Source_E0;
MCNUM mccGen_Source_dE;
MCNUM mccGen_Source_lambda0;
MCNUM mccGen_Source_dlambda;
MCNUM mccGen_Source_I1;
MCNUM mccGen_Source_yheight;
MCNUM mccGen_Source_xwidth;
MCNUM mccGen_Source_verbose;
MCNUM mccGen_Source_T1;
MCNUM mccGen_Source_flux_file_perAA;
MCNUM mccGen_Source_flux_file_log;
MCNUM mccGen_Source_Lmin;
MCNUM mccGen_Source_Lmax;
MCNUM mccGen_Source_Emin;
MCNUM mccGen_Source_Emax;
MCNUM mccGen_Source_T2;
MCNUM mccGen_Source_I2;
MCNUM mccGen_Source_T3;
MCNUM mccGen_Source_I3;
MCNUM mccGen_Source_zdepth;
int mccGen_Source_target_index;

/* Setting parameters for component 'NL1A_NL1B_NL2_NL3_InPile' [4]. */
char mccNL1A_NL1B_NL2_NL3_InPile_reflect[16384];
MCNUM mccNL1A_NL1B_NL2_NL3_InPile_w1;
MCNUM mccNL1A_NL1B_NL2_NL3_InPile_h1;
MCNUM mccNL1A_NL1B_NL2_NL3_InPile_w2;
MCNUM mccNL1A_NL1B_NL2_NL3_InPile_h2;
MCNUM mccNL1A_NL1B_NL2_NL3_InPile_l;
MCNUM mccNL1A_NL1B_NL2_NL3_InPile_R0;
MCNUM mccNL1A_NL1B_NL2_NL3_InPile_Qc;
MCNUM mccNL1A_NL1B_NL2_NL3_InPile_alpha;
MCNUM mccNL1A_NL1B_NL2_NL3_InPile_m;
MCNUM mccNL1A_NL1B_NL2_NL3_InPile_W;

/* Setting parameters for component 'NL1B_Straight1' [6]. */
char mccNL1B_Straight1_reflect[16384];
MCNUM mccNL1B_Straight1_w1;
MCNUM mccNL1B_Straight1_h1;
MCNUM mccNL1B_Straight1_w2;
MCNUM mccNL1B_Straight1_h2;
MCNUM mccNL1B_Straight1_l;
MCNUM mccNL1B_Straight1_R0;
MCNUM mccNL1B_Straight1_Qc;
MCNUM mccNL1B_Straight1_alpha;
MCNUM mccNL1B_Straight1_m;
MCNUM mccNL1B_Straight1_W;

/* Setting parameters for component 'NL1B_Curved_Guide_1' [8]. */
MCNUM mccNL1B_Curved_Guide_1_w1;
MCNUM mccNL1B_Curved_Guide_1_h1;
MCNUM mccNL1B_Curved_Guide_1_l;
MCNUM mccNL1B_Curved_Guide_1_R0;
MCNUM mccNL1B_Curved_Guide_1_Qc;
MCNUM mccNL1B_Curved_Guide_1_alpha;
MCNUM mccNL1B_Curved_Guide_1_m;
MCNUM mccNL1B_Curved_Guide_1_W;
MCNUM mccNL1B_Curved_Guide_1_curvature;

/* Setting parameters for component 'Before_selec' [10]. */
int mccBefore_selec_nx;
int mccBefore_selec_ny;
char mccBefore_selec_filename[16384];
MCNUM mccBefore_selec_xmin;
MCNUM mccBefore_selec_xmax;
MCNUM mccBefore_selec_ymin;
MCNUM mccBefore_selec_ymax;
MCNUM mccBefore_selec_xwidth;
MCNUM mccBefore_selec_yheight;
MCNUM mccBefore_selec_restore_neutron;

/* Setting parameters for component 'SELEC' [11]. */
MCNUM mccSELEC_xmin;
MCNUM mccSELEC_xmax;
MCNUM mccSELEC_ymin;
MCNUM mccSELEC_ymax;
MCNUM mccSELEC_length;
MCNUM mccSELEC_xwidth;
MCNUM mccSELEC_yheight;
MCNUM mccSELEC_nslit;
MCNUM mccSELEC_d;
MCNUM mccSELEC_radius;
MCNUM mccSELEC_alpha;
MCNUM mccSELEC_nu;

/* Setting parameters for component 'After_selec' [12]. */
int mccAfter_selec_nx;
int mccAfter_selec_ny;
char mccAfter_selec_filename[16384];
MCNUM mccAfter_selec_xmin;
MCNUM mccAfter_selec_xmax;
MCNUM mccAfter_selec_ymin;
MCNUM mccAfter_selec_ymax;
MCNUM mccAfter_selec_xwidth;
MCNUM mccAfter_selec_yheight;
MCNUM mccAfter_selec_restore_neutron;

/* Setting parameters for component 'NL1B_Curved_Guide_2' [14]. */
MCNUM mccNL1B_Curved_Guide_2_w1;
MCNUM mccNL1B_Curved_Guide_2_h1;
MCNUM mccNL1B_Curved_Guide_2_l;
MCNUM mccNL1B_Curved_Guide_2_R0;
MCNUM mccNL1B_Curved_Guide_2_Qc;
MCNUM mccNL1B_Curved_Guide_2_alpha;
MCNUM mccNL1B_Curved_Guide_2_m;
MCNUM mccNL1B_Curved_Guide_2_W;
MCNUM mccNL1B_Curved_Guide_2_curvature;

/* Setting parameters for component 'NL1B_Straight2' [16]. */
char mccNL1B_Straight2_reflect[16384];
MCNUM mccNL1B_Straight2_w1;
MCNUM mccNL1B_Straight2_h1;
MCNUM mccNL1B_Straight2_w2;
MCNUM mccNL1B_Straight2_h2;
MCNUM mccNL1B_Straight2_l;
MCNUM mccNL1B_Straight2_R0;
MCNUM mccNL1B_Straight2_Qc;
MCNUM mccNL1B_Straight2_alpha;
MCNUM mccNL1B_Straight2_m;
MCNUM mccNL1B_Straight2_W;

/* Setting parameters for component 'elliptical_piece' [18]. */
char mccelliptical_piece_option[16384];
MCNUM mccelliptical_piece_w1;
MCNUM mccelliptical_piece_h1;
MCNUM mccelliptical_piece_l;
MCNUM mccelliptical_piece_linw;
MCNUM mccelliptical_piece_loutw;
MCNUM mccelliptical_piece_linh;
MCNUM mccelliptical_piece_louth;
MCNUM mccelliptical_piece_R0;
MCNUM mccelliptical_piece_Qcx;
MCNUM mccelliptical_piece_Qcy;
MCNUM mccelliptical_piece_alphax;
MCNUM mccelliptical_piece_alphay;
MCNUM mccelliptical_piece_W;
MCNUM mccelliptical_piece_mx;
MCNUM mccelliptical_piece_my;
MCNUM mccelliptical_piece_segno;
MCNUM mccelliptical_piece_curvature;
MCNUM mccelliptical_piece_curvature_v;

/* Setting parameters for component 'Virtual_source' [19]. */
MCNUM mccVirtual_source_xmin;
MCNUM mccVirtual_source_xmax;
MCNUM mccVirtual_source_ymin;
MCNUM mccVirtual_source_ymax;
MCNUM mccVirtual_source_radius;
MCNUM mccVirtual_source_xwidth;
MCNUM mccVirtual_source_yheight;

/* Setting parameters for component 'NL1B_Vertical_Guide' [21]. */
MCNUM mccNL1B_Vertical_Guide_w1;
MCNUM mccNL1B_Vertical_Guide_h1;
MCNUM mccNL1B_Vertical_Guide_w2;
MCNUM mccNL1B_Vertical_Guide_h2;
MCNUM mccNL1B_Vertical_Guide_l;
MCNUM mccNL1B_Vertical_Guide_R0;
MCNUM mccNL1B_Vertical_Guide_Qc;
MCNUM mccNL1B_Vertical_Guide_alpha;
MCNUM mccNL1B_Vertical_Guide_m;
MCNUM mccNL1B_Vertical_Guide_nslit;
MCNUM mccNL1B_Vertical_Guide_d;
MCNUM mccNL1B_Vertical_Guide_Qcx;
MCNUM mccNL1B_Vertical_Guide_Qcy;
MCNUM mccNL1B_Vertical_Guide_alphax;
MCNUM mccNL1B_Vertical_Guide_alphay;
MCNUM mccNL1B_Vertical_Guide_W;
MCNUM mccNL1B_Vertical_Guide_mx;
MCNUM mccNL1B_Vertical_Guide_my;
MCNUM mccNL1B_Vertical_Guide_nu;
MCNUM mccNL1B_Vertical_Guide_phase;

/* Definition parameters for component 'energy_endguide' [23]. */
#define mccenergy_endguide_nE 200
/* Setting parameters for component 'energy_endguide' [23]. */
char mccenergy_endguide_filename[16384];
MCNUM mccenergy_endguide_xmin;
MCNUM mccenergy_endguide_xmax;
MCNUM mccenergy_endguide_ymin;
MCNUM mccenergy_endguide_ymax;
MCNUM mccenergy_endguide_xwidth;
MCNUM mccenergy_endguide_yheight;
MCNUM mccenergy_endguide_Emin;
MCNUM mccenergy_endguide_Emax;
MCNUM mccenergy_endguide_restore_neutron;
int mccenergy_endguide_nowritefile;

/* Setting parameters for component 'Monochromator' [25]. */
char mccMonochromator_reflect[16384];
char mccMonochromator_transmit[16384];
MCNUM mccMonochromator_zwidth;
MCNUM mccMonochromator_yheight;
MCNUM mccMonochromator_gap;
MCNUM mccMonochromator_NH;
MCNUM mccMonochromator_NV;
MCNUM mccMonochromator_mosaich;
MCNUM mccMonochromator_mosaicv;
MCNUM mccMonochromator_r0;
MCNUM mccMonochromator_t0;
MCNUM mccMonochromator_Q;
MCNUM mccMonochromator_RV;
MCNUM mccMonochromator_RH;
MCNUM mccMonochromator_DM;
MCNUM mccMonochromator_mosaic;
MCNUM mccMonochromator_width;
MCNUM mccMonochromator_height;
MCNUM mccMonochromator_verbose;
MCNUM mccMonochromator_order;

/* Definition parameters for component 'energy_mono' [27]. */
#define mccenergy_mono_nE 200
/* Setting parameters for component 'energy_mono' [27]. */
char mccenergy_mono_filename[16384];
MCNUM mccenergy_mono_xmin;
MCNUM mccenergy_mono_xmax;
MCNUM mccenergy_mono_ymin;
MCNUM mccenergy_mono_ymax;
MCNUM mccenergy_mono_xwidth;
MCNUM mccenergy_mono_yheight;
MCNUM mccenergy_mono_Emin;
MCNUM mccenergy_mono_Emax;
MCNUM mccenergy_mono_restore_neutron;
int mccenergy_mono_nowritefile;

/* Definition parameters for component 'energy_pre_sample' [28]. */
#define mccenergy_pre_sample_nE 50
/* Setting parameters for component 'energy_pre_sample' [28]. */
char mccenergy_pre_sample_filename[16384];
MCNUM mccenergy_pre_sample_xmin;
MCNUM mccenergy_pre_sample_xmax;
MCNUM mccenergy_pre_sample_ymin;
MCNUM mccenergy_pre_sample_ymax;
MCNUM mccenergy_pre_sample_xwidth;
MCNUM mccenergy_pre_sample_yheight;
MCNUM mccenergy_pre_sample_Emin;
MCNUM mccenergy_pre_sample_Emax;
MCNUM mccenergy_pre_sample_restore_neutron;
int mccenergy_pre_sample_nowritefile;

/* Definition parameters for component 'div_mono' [31]. */
#define mccdiv_mono_nh 20
#define mccdiv_mono_nv 20
/* Setting parameters for component 'div_mono' [31]. */
char mccdiv_mono_filename[16384];
MCNUM mccdiv_mono_xmin;
MCNUM mccdiv_mono_xmax;
MCNUM mccdiv_mono_ymin;
MCNUM mccdiv_mono_ymax;
MCNUM mccdiv_mono_xwidth;
MCNUM mccdiv_mono_yheight;
MCNUM mccdiv_mono_maxdiv_h;
MCNUM mccdiv_mono_maxdiv_v;
MCNUM mccdiv_mono_restore_neutron;
MCNUM mccdiv_mono_nx;
MCNUM mccdiv_mono_ny;
MCNUM mccdiv_mono_nz;
int mccdiv_mono_nowritefile;

/* Definition parameters for component 'div_mono_H' [32]. */
#define mccdiv_mono_H_nh 35
/* Setting parameters for component 'div_mono_H' [32]. */
char mccdiv_mono_H_filename[16384];
MCNUM mccdiv_mono_H_xmin;
MCNUM mccdiv_mono_H_xmax;
MCNUM mccdiv_mono_H_ymin;
MCNUM mccdiv_mono_H_ymax;
MCNUM mccdiv_mono_H_xwidth;
MCNUM mccdiv_mono_H_yheight;
MCNUM mccdiv_mono_H_h_maxdiv;
MCNUM mccdiv_mono_H_restore_neutron;
int mccdiv_mono_H_nowritefile;

/* Setting parameters for component 'psd_sam' [33]. */
int mccpsd_sam_nx;
int mccpsd_sam_ny;
char mccpsd_sam_filename[16384];
MCNUM mccpsd_sam_xmin;
MCNUM mccpsd_sam_xmax;
MCNUM mccpsd_sam_ymin;
MCNUM mccpsd_sam_ymax;
MCNUM mccpsd_sam_xwidth;
MCNUM mccpsd_sam_yheight;
MCNUM mccpsd_sam_restore_neutron;

/* Definition parameters for component '1dpsd' [34]. */
#define mcc1dpsd_nx 30
/* Setting parameters for component '1dpsd' [34]. */
char mcc1dpsd_filename[16384];
MCNUM mcc1dpsd_xmin;
MCNUM mcc1dpsd_xmax;
MCNUM mcc1dpsd_ymin;
MCNUM mcc1dpsd_ymax;
MCNUM mcc1dpsd_xwidth;
MCNUM mcc1dpsd_yheight;
MCNUM mcc1dpsd_restore_neutron;
int mcc1dpsd_nowritefile;

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
#line 7589 "./HZB_FLEX.c"
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

/* User declarations for component 'Gen_Source' [2]. */
#define mccompcurname  Gen_Source
#define mccompcurtype  Source_gen
#define mccompcurindex 2
#define p_in mccGen_Source_p_in
#define lambda1 mccGen_Source_lambda1
#define lambda2 mccGen_Source_lambda2
#define lambda3 mccGen_Source_lambda3
#define pTable mccGen_Source_pTable
#define pTable_x mccGen_Source_pTable_x
#define pTable_y mccGen_Source_pTable_y
#define pTable_xmin mccGen_Source_pTable_xmin
#define pTable_xmax mccGen_Source_pTable_xmax
#define pTable_xsum mccGen_Source_pTable_xsum
#define pTable_ymin mccGen_Source_pTable_ymin
#define pTable_ymax mccGen_Source_pTable_ymax
#define pTable_ysum mccGen_Source_pTable_ysum
#define pTable_dxmin mccGen_Source_pTable_dxmin
#define pTable_dxmax mccGen_Source_pTable_dxmax
#define pTable_dymin mccGen_Source_pTable_dymin
#define pTable_dymax mccGen_Source_pTable_dymax
#define flux_file mccGen_Source_flux_file
#define xdiv_file mccGen_Source_xdiv_file
#define ydiv_file mccGen_Source_ydiv_file
#define radius mccGen_Source_radius
#define dist mccGen_Source_dist
#define focus_xw mccGen_Source_focus_xw
#define focus_yh mccGen_Source_focus_yh
#define focus_aw mccGen_Source_focus_aw
#define focus_ah mccGen_Source_focus_ah
#define E0 mccGen_Source_E0
#define dE mccGen_Source_dE
#define lambda0 mccGen_Source_lambda0
#define dlambda mccGen_Source_dlambda
#define I1 mccGen_Source_I1
#define yheight mccGen_Source_yheight
#define xwidth mccGen_Source_xwidth
#define verbose mccGen_Source_verbose
#define T1 mccGen_Source_T1
#define flux_file_perAA mccGen_Source_flux_file_perAA
#define flux_file_log mccGen_Source_flux_file_log
#define Lmin mccGen_Source_Lmin
#define Lmax mccGen_Source_Lmax
#define Emin mccGen_Source_Emin
#define Emax mccGen_Source_Emax
#define T2 mccGen_Source_T2
#define I2 mccGen_Source_I2
#define T3 mccGen_Source_T3
#define I3 mccGen_Source_I3
#define zdepth mccGen_Source_zdepth
#define target_index mccGen_Source_target_index
#line 184 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../sources/Source_gen.comp"

  double p_in;
  double lambda1;  /* first Maxwellian source */
  double lambda2;  /* second Maxwellian source */
  double lambda3;  /* third Maxwellian source */
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

#line 7673 "./HZB_FLEX.c"
#undef target_index
#undef zdepth
#undef I3
#undef T3
#undef I2
#undef T2
#undef Emax
#undef Emin
#undef Lmax
#undef Lmin
#undef flux_file_log
#undef flux_file_perAA
#undef T1
#undef verbose
#undef xwidth
#undef yheight
#undef I1
#undef dlambda
#undef lambda0
#undef dE
#undef E0
#undef focus_ah
#undef focus_aw
#undef focus_yh
#undef focus_xw
#undef dist
#undef radius
#undef ydiv_file
#undef xdiv_file
#undef flux_file
#undef pTable_dymax
#undef pTable_dymin
#undef pTable_dxmax
#undef pTable_dxmin
#undef pTable_ysum
#undef pTable_ymax
#undef pTable_ymin
#undef pTable_xsum
#undef pTable_xmax
#undef pTable_xmin
#undef pTable_y
#undef pTable_x
#undef pTable
#undef lambda3
#undef lambda2
#undef lambda1
#undef p_in
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'NL1A_NL1B_NL2_NL3_InPile_Entrance_Window' [3]. */
#define mccompcurname  NL1A_NL1B_NL2_NL3_InPile_Entrance_Window
#define mccompcurtype  Arm
#define mccompcurindex 3
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'NL1A_NL1B_NL2_NL3_InPile' [4]. */
#define mccompcurname  NL1A_NL1B_NL2_NL3_InPile
#define mccompcurtype  Guide
#define mccompcurindex 4
#define pTable mccNL1A_NL1B_NL2_NL3_InPile_pTable
#define reflect mccNL1A_NL1B_NL2_NL3_InPile_reflect
#define w1 mccNL1A_NL1B_NL2_NL3_InPile_w1
#define h1 mccNL1A_NL1B_NL2_NL3_InPile_h1
#define w2 mccNL1A_NL1B_NL2_NL3_InPile_w2
#define h2 mccNL1A_NL1B_NL2_NL3_InPile_h2
#define l mccNL1A_NL1B_NL2_NL3_InPile_l
#define R0 mccNL1A_NL1B_NL2_NL3_InPile_R0
#define Qc mccNL1A_NL1B_NL2_NL3_InPile_Qc
#define alpha mccNL1A_NL1B_NL2_NL3_InPile_alpha
#define m mccNL1A_NL1B_NL2_NL3_InPile_m
#define W mccNL1A_NL1B_NL2_NL3_InPile_W
#line 70 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide.comp"
t_Table pTable;
#line 7751 "./HZB_FLEX.c"
#undef W
#undef m
#undef alpha
#undef Qc
#undef R0
#undef l
#undef h2
#undef w2
#undef h1
#undef w1
#undef reflect
#undef pTable
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'NL1B_Straight1_Entrance_Window' [5]. */
#define mccompcurname  NL1B_Straight1_Entrance_Window
#define mccompcurtype  Arm
#define mccompcurindex 5
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'NL1B_Straight1' [6]. */
#define mccompcurname  NL1B_Straight1
#define mccompcurtype  Guide
#define mccompcurindex 6
#define pTable mccNL1B_Straight1_pTable
#define reflect mccNL1B_Straight1_reflect
#define w1 mccNL1B_Straight1_w1
#define h1 mccNL1B_Straight1_h1
#define w2 mccNL1B_Straight1_w2
#define h2 mccNL1B_Straight1_h2
#define l mccNL1B_Straight1_l
#define R0 mccNL1B_Straight1_R0
#define Qc mccNL1B_Straight1_Qc
#define alpha mccNL1B_Straight1_alpha
#define m mccNL1B_Straight1_m
#define W mccNL1B_Straight1_W
#line 70 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide.comp"
t_Table pTable;
#line 7794 "./HZB_FLEX.c"
#undef W
#undef m
#undef alpha
#undef Qc
#undef R0
#undef l
#undef h2
#undef w2
#undef h1
#undef w1
#undef reflect
#undef pTable
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'NL1B_Curved_Entrance_Window' [7]. */
#define mccompcurname  NL1B_Curved_Entrance_Window
#define mccompcurtype  Arm
#define mccompcurindex 7
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'NL1B_Velocity_Selector_Gap_Entrance_Window' [9]. */
#define mccompcurname  NL1B_Velocity_Selector_Gap_Entrance_Window
#define mccompcurtype  Arm
#define mccompcurindex 9
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'Before_selec' [10]. */
#define mccompcurname  Before_selec
#define mccompcurtype  PSD_monitor
#define mccompcurindex 10
#define PSD_N mccBefore_selec_PSD_N
#define PSD_p mccBefore_selec_PSD_p
#define PSD_p2 mccBefore_selec_PSD_p2
#define nx mccBefore_selec_nx
#define ny mccBefore_selec_ny
#define filename mccBefore_selec_filename
#define xmin mccBefore_selec_xmin
#define xmax mccBefore_selec_xmax
#define ymin mccBefore_selec_ymin
#define ymax mccBefore_selec_ymax
#define xwidth mccBefore_selec_xwidth
#define yheight mccBefore_selec_yheight
#define restore_neutron mccBefore_selec_restore_neutron
#line 62 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
  DArray2d PSD_N;
  DArray2d PSD_p;
  DArray2d PSD_p2;
#line 7848 "./HZB_FLEX.c"
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

/* User declarations for component 'SELEC' [11]. */
#define mccompcurname  SELEC
#define mccompcurtype  Selector
#define mccompcurindex 11
#define xmin mccSELEC_xmin
#define xmax mccSELEC_xmax
#define ymin mccSELEC_ymin
#define ymax mccSELEC_ymax
#define length mccSELEC_length
#define xwidth mccSELEC_xwidth
#define yheight mccSELEC_yheight
#define nslit mccSELEC_nslit
#define d mccSELEC_d
#define radius mccSELEC_radius
#define alpha mccSELEC_alpha
#define nu mccSELEC_nu
#undef nu
#undef alpha
#undef radius
#undef d
#undef nslit
#undef yheight
#undef xwidth
#undef length
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'After_selec' [12]. */
#define mccompcurname  After_selec
#define mccompcurtype  PSD_monitor
#define mccompcurindex 12
#define PSD_N mccAfter_selec_PSD_N
#define PSD_p mccAfter_selec_PSD_p
#define PSD_p2 mccAfter_selec_PSD_p2
#define nx mccAfter_selec_nx
#define ny mccAfter_selec_ny
#define filename mccAfter_selec_filename
#define xmin mccAfter_selec_xmin
#define xmax mccAfter_selec_xmax
#define ymin mccAfter_selec_ymin
#define ymax mccAfter_selec_ymax
#define xwidth mccAfter_selec_xwidth
#define yheight mccAfter_selec_yheight
#define restore_neutron mccAfter_selec_restore_neutron
#line 62 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
  DArray2d PSD_N;
  DArray2d PSD_p;
  DArray2d PSD_p2;
#line 7919 "./HZB_FLEX.c"
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

/* User declarations for component 'NL1B_Curved_2_Entrance_Window' [13]. */
#define mccompcurname  NL1B_Curved_2_Entrance_Window
#define mccompcurtype  Arm
#define mccompcurindex 13
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'NL1B_Straight2_Entrance_Window' [15]. */
#define mccompcurname  NL1B_Straight2_Entrance_Window
#define mccompcurtype  Arm
#define mccompcurindex 15
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'NL1B_Straight2' [16]. */
#define mccompcurname  NL1B_Straight2
#define mccompcurtype  Guide
#define mccompcurindex 16
#define pTable mccNL1B_Straight2_pTable
#define reflect mccNL1B_Straight2_reflect
#define w1 mccNL1B_Straight2_w1
#define h1 mccNL1B_Straight2_h1
#define w2 mccNL1B_Straight2_w2
#define h2 mccNL1B_Straight2_h2
#define l mccNL1B_Straight2_l
#define R0 mccNL1B_Straight2_R0
#define Qc mccNL1B_Straight2_Qc
#define alpha mccNL1B_Straight2_alpha
#define m mccNL1B_Straight2_m
#define W mccNL1B_Straight2_W
#line 70 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide.comp"
t_Table pTable;
#line 7971 "./HZB_FLEX.c"
#undef W
#undef m
#undef alpha
#undef Qc
#undef R0
#undef l
#undef h2
#undef w2
#undef h1
#undef w1
#undef reflect
#undef pTable
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'NL1B_Elliptical_Entrance_Window' [17]. */
#define mccompcurname  NL1B_Elliptical_Entrance_Window
#define mccompcurtype  Arm
#define mccompcurindex 17
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'elliptical_piece' [18]. */
#define mccompcurname  elliptical_piece
#define mccompcurtype  Guide_tapering
#define mccompcurindex 18
#define w1c mccelliptical_piece_w1c
#define w2c mccelliptical_piece_w2c
#define ww mccelliptical_piece_ww
#define hh mccelliptical_piece_hh
#define whalf mccelliptical_piece_whalf
#define hhalf mccelliptical_piece_hhalf
#define lwhalf mccelliptical_piece_lwhalf
#define lhhalf mccelliptical_piece_lhhalf
#define h1_in mccelliptical_piece_h1_in
#define h2_out mccelliptical_piece_h2_out
#define w1_in mccelliptical_piece_w1_in
#define w2_out mccelliptical_piece_w2_out
#define l_seg mccelliptical_piece_l_seg
#define seg mccelliptical_piece_seg
#define h12 mccelliptical_piece_h12
#define h2 mccelliptical_piece_h2
#define w12 mccelliptical_piece_w12
#define w2 mccelliptical_piece_w2
#define a_ell_q mccelliptical_piece_a_ell_q
#define b_ell_q mccelliptical_piece_b_ell_q
#define lbw mccelliptical_piece_lbw
#define lbh mccelliptical_piece_lbh
#define mxi mccelliptical_piece_mxi
#define u1 mccelliptical_piece_u1
#define u2 mccelliptical_piece_u2
#define div1 mccelliptical_piece_div1
#define p2_para mccelliptical_piece_p2_para
#define test mccelliptical_piece_test
#define Div1 mccelliptical_piece_Div1
#define i mccelliptical_piece_i
#define ii mccelliptical_piece_ii
#define seg mccelliptical_piece_seg
#define fu mccelliptical_piece_fu
#define pos mccelliptical_piece_pos
#define file_name mccelliptical_piece_file_name
#define ep mccelliptical_piece_ep
#define num mccelliptical_piece_num
#define rotation_h mccelliptical_piece_rotation_h
#define rotation_v mccelliptical_piece_rotation_v
#define option mccelliptical_piece_option
#define w1 mccelliptical_piece_w1
#define h1 mccelliptical_piece_h1
#define l mccelliptical_piece_l
#define linw mccelliptical_piece_linw
#define loutw mccelliptical_piece_loutw
#define linh mccelliptical_piece_linh
#define louth mccelliptical_piece_louth
#define R0 mccelliptical_piece_R0
#define Qcx mccelliptical_piece_Qcx
#define Qcy mccelliptical_piece_Qcy
#define alphax mccelliptical_piece_alphax
#define alphay mccelliptical_piece_alphay
#define W mccelliptical_piece_W
#define mx mccelliptical_piece_mx
#define my mccelliptical_piece_my
#define segno mccelliptical_piece_segno
#define curvature mccelliptical_piece_curvature
#define curvature_v mccelliptical_piece_curvature_v
#line 97 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide_tapering.comp"
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
#line 8075 "./HZB_FLEX.c"
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

/* User declarations for component 'Virtual_source' [19]. */
#define mccompcurname  Virtual_source
#define mccompcurtype  Slit
#define mccompcurindex 19
#define xmin mccVirtual_source_xmin
#define xmax mccVirtual_source_xmax
#define ymin mccVirtual_source_ymin
#define ymax mccVirtual_source_ymax
#define radius mccVirtual_source_radius
#define xwidth mccVirtual_source_xwidth
#define yheight mccVirtual_source_yheight
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

/* User declarations for component 'NL1B_Vertical_Guide_Entrance_Window' [20]. */
#define mccompcurname  NL1B_Vertical_Guide_Entrance_Window
#define mccompcurtype  Arm
#define mccompcurindex 20
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'NL1B_Vertical_Guide' [21]. */
#define mccompcurname  NL1B_Vertical_Guide
#define mccompcurtype  Guide_channeled
#define mccompcurindex 21
#define w1c mccNL1B_Vertical_Guide_w1c
#define w2c mccNL1B_Vertical_Guide_w2c
#define ww mccNL1B_Vertical_Guide_ww
#define hh mccNL1B_Vertical_Guide_hh
#define whalf mccNL1B_Vertical_Guide_whalf
#define hhalf mccNL1B_Vertical_Guide_hhalf
#define lwhalf mccNL1B_Vertical_Guide_lwhalf
#define lhhalf mccNL1B_Vertical_Guide_lhhalf
#define w1 mccNL1B_Vertical_Guide_w1
#define h1 mccNL1B_Vertical_Guide_h1
#define w2 mccNL1B_Vertical_Guide_w2
#define h2 mccNL1B_Vertical_Guide_h2
#define l mccNL1B_Vertical_Guide_l
#define R0 mccNL1B_Vertical_Guide_R0
#define Qc mccNL1B_Vertical_Guide_Qc
#define alpha mccNL1B_Vertical_Guide_alpha
#define m mccNL1B_Vertical_Guide_m
#define nslit mccNL1B_Vertical_Guide_nslit
#define d mccNL1B_Vertical_Guide_d
#define Qcx mccNL1B_Vertical_Guide_Qcx
#define Qcy mccNL1B_Vertical_Guide_Qcy
#define alphax mccNL1B_Vertical_Guide_alphax
#define alphay mccNL1B_Vertical_Guide_alphay
#define W mccNL1B_Vertical_Guide_W
#define mx mccNL1B_Vertical_Guide_mx
#define my mccNL1B_Vertical_Guide_my
#define nu mccNL1B_Vertical_Guide_nu
#define phase mccNL1B_Vertical_Guide_phase
#line 81 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide_channeled.comp"
double w1c;
double w2c;
double ww, hh;
double whalf, hhalf;
double lwhalf, lhhalf;
#line 8206 "./HZB_FLEX.c"
#undef phase
#undef nu
#undef my
#undef mx
#undef W
#undef alphay
#undef alphax
#undef Qcy
#undef Qcx
#undef d
#undef nslit
#undef m
#undef alpha
#undef Qc
#undef R0
#undef l
#undef h2
#undef w2
#undef h1
#undef w1
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

/* User declarations for component 'NL1B_Guide_Exit' [22]. */
#define mccompcurname  NL1B_Guide_Exit
#define mccompcurtype  Arm
#define mccompcurindex 22
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'energy_endguide' [23]. */
#define mccompcurname  energy_endguide
#define mccompcurtype  E_monitor
#define mccompcurindex 23
#define nE mccenergy_endguide_nE
#define E_N mccenergy_endguide_E_N
#define E_p mccenergy_endguide_E_p
#define E_p2 mccenergy_endguide_E_p2
#define S_p mccenergy_endguide_S_p
#define S_pE mccenergy_endguide_S_pE
#define S_pE2 mccenergy_endguide_S_pE2
#define filename mccenergy_endguide_filename
#define xmin mccenergy_endguide_xmin
#define xmax mccenergy_endguide_xmax
#define ymin mccenergy_endguide_ymin
#define ymax mccenergy_endguide_ymax
#define xwidth mccenergy_endguide_xwidth
#define yheight mccenergy_endguide_yheight
#define Emin mccenergy_endguide_Emin
#define Emax mccenergy_endguide_Emax
#define restore_neutron mccenergy_endguide_restore_neutron
#define nowritefile mccenergy_endguide_nowritefile
#line 60 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/E_monitor.comp"
double E_N[nE];
double E_p[nE], E_p2[nE];
double S_p, S_pE, S_pE2;
#line 8273 "./HZB_FLEX.c"
#undef nowritefile
#undef restore_neutron
#undef Emax
#undef Emin
#undef yheight
#undef xwidth
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef filename
#undef S_pE2
#undef S_pE
#undef S_p
#undef E_p2
#undef E_p
#undef E_N
#undef nE
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'Mono_center' [24]. */
#define mccompcurname  Mono_center
#define mccompcurtype  Arm
#define mccompcurindex 24
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'Monochromator' [25]. */
#define mccompcurname  Monochromator
#define mccompcurtype  Monochromator_curved
#define mccompcurindex 25
#define mos_rms_y mccMonochromator_mos_rms_y
#define mos_rms_z mccMonochromator_mos_rms_z
#define mos_rms_max mccMonochromator_mos_rms_max
#define mono_Q mccMonochromator_mono_Q
#define SlabWidth mccMonochromator_SlabWidth
#define SlabHeight mccMonochromator_SlabHeight
#define rTable mccMonochromator_rTable
#define tTable mccMonochromator_tTable
#define row mccMonochromator_row
#define col mccMonochromator_col
#define tiltH mccMonochromator_tiltH
#define tiltV mccMonochromator_tiltV
#define reflect mccMonochromator_reflect
#define transmit mccMonochromator_transmit
#define zwidth mccMonochromator_zwidth
#define yheight mccMonochromator_yheight
#define gap mccMonochromator_gap
#define NH mccMonochromator_NH
#define NV mccMonochromator_NV
#define mosaich mccMonochromator_mosaich
#define mosaicv mccMonochromator_mosaicv
#define r0 mccMonochromator_r0
#define t0 mccMonochromator_t0
#define Q mccMonochromator_Q
#define RV mccMonochromator_RV
#define RH mccMonochromator_RH
#define DM mccMonochromator_DM
#define mosaic mccMonochromator_mosaic
#define width mccMonochromator_width
#define height mccMonochromator_height
#define verbose mccMonochromator_verbose
#define order mccMonochromator_order
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
#line 8350 "./HZB_FLEX.c"
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

/* User declarations for component 'Mono_sample_arm' [26]. */
#define mccompcurname  Mono_sample_arm
#define mccompcurtype  Arm
#define mccompcurindex 26
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'energy_mono' [27]. */
#define mccompcurname  energy_mono
#define mccompcurtype  E_monitor
#define mccompcurindex 27
#define nE mccenergy_mono_nE
#define E_N mccenergy_mono_E_N
#define E_p mccenergy_mono_E_p
#define E_p2 mccenergy_mono_E_p2
#define S_p mccenergy_mono_S_p
#define S_pE mccenergy_mono_S_pE
#define S_pE2 mccenergy_mono_S_pE2
#define filename mccenergy_mono_filename
#define xmin mccenergy_mono_xmin
#define xmax mccenergy_mono_xmax
#define ymin mccenergy_mono_ymin
#define ymax mccenergy_mono_ymax
#define xwidth mccenergy_mono_xwidth
#define yheight mccenergy_mono_yheight
#define Emin mccenergy_mono_Emin
#define Emax mccenergy_mono_Emax
#define restore_neutron mccenergy_mono_restore_neutron
#define nowritefile mccenergy_mono_nowritefile
#line 60 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/E_monitor.comp"
double E_N[nE];
double E_p[nE], E_p2[nE];
double S_p, S_pE, S_pE2;
#line 8421 "./HZB_FLEX.c"
#undef nowritefile
#undef restore_neutron
#undef Emax
#undef Emin
#undef yheight
#undef xwidth
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef filename
#undef S_pE2
#undef S_pE
#undef S_p
#undef E_p2
#undef E_p
#undef E_N
#undef nE
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'energy_pre_sample' [28]. */
#define mccompcurname  energy_pre_sample
#define mccompcurtype  E_monitor
#define mccompcurindex 28
#define nE mccenergy_pre_sample_nE
#define E_N mccenergy_pre_sample_E_N
#define E_p mccenergy_pre_sample_E_p
#define E_p2 mccenergy_pre_sample_E_p2
#define S_p mccenergy_pre_sample_S_p
#define S_pE mccenergy_pre_sample_S_pE
#define S_pE2 mccenergy_pre_sample_S_pE2
#define filename mccenergy_pre_sample_filename
#define xmin mccenergy_pre_sample_xmin
#define xmax mccenergy_pre_sample_xmax
#define ymin mccenergy_pre_sample_ymin
#define ymax mccenergy_pre_sample_ymax
#define xwidth mccenergy_pre_sample_xwidth
#define yheight mccenergy_pre_sample_yheight
#define Emin mccenergy_pre_sample_Emin
#define Emax mccenergy_pre_sample_Emax
#define restore_neutron mccenergy_pre_sample_restore_neutron
#define nowritefile mccenergy_pre_sample_nowritefile
#line 60 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/E_monitor.comp"
double E_N[nE];
double E_p[nE], E_p2[nE];
double S_p, S_pE, S_pE2;
#line 8470 "./HZB_FLEX.c"
#undef nowritefile
#undef restore_neutron
#undef Emax
#undef Emin
#undef yheight
#undef xwidth
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef filename
#undef S_pE2
#undef S_pE
#undef S_p
#undef E_p2
#undef E_p
#undef E_N
#undef nE
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'Sample_center' [29]. */
#define mccompcurname  Sample_center
#define mccompcurtype  Arm
#define mccompcurindex 29
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'Sample_analyser_arm' [30]. */
#define mccompcurname  Sample_analyser_arm
#define mccompcurtype  Arm
#define mccompcurindex 30
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'div_mono' [31]. */
#define mccompcurname  div_mono
#define mccompcurtype  Divergence_monitor
#define mccompcurindex 31
#define nh mccdiv_mono_nh
#define nv mccdiv_mono_nv
#define Div_N mccdiv_mono_Div_N
#define Div_p mccdiv_mono_Div_p
#define Div_p2 mccdiv_mono_Div_p2
#define filename mccdiv_mono_filename
#define xmin mccdiv_mono_xmin
#define xmax mccdiv_mono_xmax
#define ymin mccdiv_mono_ymin
#define ymax mccdiv_mono_ymax
#define xwidth mccdiv_mono_xwidth
#define yheight mccdiv_mono_yheight
#define maxdiv_h mccdiv_mono_maxdiv_h
#define maxdiv_v mccdiv_mono_maxdiv_v
#define restore_neutron mccdiv_mono_restore_neutron
#define nx mccdiv_mono_nx
#define ny mccdiv_mono_ny
#define nz mccdiv_mono_nz
#define nowritefile mccdiv_mono_nowritefile
#line 61 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/Divergence_monitor.comp"
double Div_N[nh][nv];
double Div_p[nh][nv];
double Div_p2[nh][nv];
#line 8536 "./HZB_FLEX.c"
#undef nowritefile
#undef nz
#undef ny
#undef nx
#undef restore_neutron
#undef maxdiv_v
#undef maxdiv_h
#undef yheight
#undef xwidth
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef filename
#undef Div_p2
#undef Div_p
#undef Div_N
#undef nv
#undef nh
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'div_mono_H' [32]. */
#define mccompcurname  div_mono_H
#define mccompcurtype  Hdiv_monitor
#define mccompcurindex 32
#define nh mccdiv_mono_H_nh
#define Div_N mccdiv_mono_H_Div_N
#define Div_p mccdiv_mono_H_Div_p
#define Div_p2 mccdiv_mono_H_Div_p2
#define filename mccdiv_mono_H_filename
#define xmin mccdiv_mono_H_xmin
#define xmax mccdiv_mono_H_xmax
#define ymin mccdiv_mono_H_ymin
#define ymax mccdiv_mono_H_ymax
#define xwidth mccdiv_mono_H_xwidth
#define yheight mccdiv_mono_H_yheight
#define h_maxdiv mccdiv_mono_H_h_maxdiv
#define restore_neutron mccdiv_mono_H_restore_neutron
#define nowritefile mccdiv_mono_H_nowritefile
#line 57 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/Hdiv_monitor.comp"
double Div_N[nh];
double Div_p[nh];
double Div_p2[nh];
#line 8582 "./HZB_FLEX.c"
#undef nowritefile
#undef restore_neutron
#undef h_maxdiv
#undef yheight
#undef xwidth
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef filename
#undef Div_p2
#undef Div_p
#undef Div_N
#undef nh
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'psd_sam' [33]. */
#define mccompcurname  psd_sam
#define mccompcurtype  PSD_monitor
#define mccompcurindex 33
#define PSD_N mccpsd_sam_PSD_N
#define PSD_p mccpsd_sam_PSD_p
#define PSD_p2 mccpsd_sam_PSD_p2
#define nx mccpsd_sam_nx
#define ny mccpsd_sam_ny
#define filename mccpsd_sam_filename
#define xmin mccpsd_sam_xmin
#define xmax mccpsd_sam_xmax
#define ymin mccpsd_sam_ymin
#define ymax mccpsd_sam_ymax
#define xwidth mccpsd_sam_xwidth
#define yheight mccpsd_sam_yheight
#define restore_neutron mccpsd_sam_restore_neutron
#line 62 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
  DArray2d PSD_N;
  DArray2d PSD_p;
  DArray2d PSD_p2;
#line 8622 "./HZB_FLEX.c"
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

/* User declarations for component '1dpsd' [34]. */
#define mccompcurname  1dpsd
#define mccompcurtype  PSDlin_monitor
#define mccompcurindex 34
#define nx mcc1dpsd_nx
#define PSDlin_N mcc1dpsd_PSDlin_N
#define PSDlin_p mcc1dpsd_PSDlin_p
#define PSDlin_p2 mcc1dpsd_PSDlin_p2
#define filename mcc1dpsd_filename
#define xmin mcc1dpsd_xmin
#define xmax mcc1dpsd_xmax
#define ymin mcc1dpsd_ymin
#define ymax mcc1dpsd_ymax
#define xwidth mcc1dpsd_xwidth
#define yheight mcc1dpsd_yheight
#define restore_neutron mcc1dpsd_restore_neutron
#define nowritefile mcc1dpsd_nowritefile
#line 53 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/PSDlin_monitor.comp"
    double PSDlin_N[nx];
    double PSDlin_p[nx];
    double PSDlin_p2[nx];
#line 8661 "./HZB_FLEX.c"
#undef nowritefile
#undef restore_neutron
#undef yheight
#undef xwidth
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef filename
#undef PSDlin_p2
#undef PSDlin_p
#undef PSDlin_N
#undef nx
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

Coords mcposaOrigin, mcposrOrigin;
Rotation mcrotaOrigin, mcrotrOrigin;
Coords mcposaGen_Source, mcposrGen_Source;
Rotation mcrotaGen_Source, mcrotrGen_Source;
Coords mcposaNL1A_NL1B_NL2_NL3_InPile_Entrance_Window, mcposrNL1A_NL1B_NL2_NL3_InPile_Entrance_Window;
Rotation mcrotaNL1A_NL1B_NL2_NL3_InPile_Entrance_Window, mcrotrNL1A_NL1B_NL2_NL3_InPile_Entrance_Window;
Coords mcposaNL1A_NL1B_NL2_NL3_InPile, mcposrNL1A_NL1B_NL2_NL3_InPile;
Rotation mcrotaNL1A_NL1B_NL2_NL3_InPile, mcrotrNL1A_NL1B_NL2_NL3_InPile;
Coords mcposaNL1B_Straight1_Entrance_Window, mcposrNL1B_Straight1_Entrance_Window;
Rotation mcrotaNL1B_Straight1_Entrance_Window, mcrotrNL1B_Straight1_Entrance_Window;
Coords mcposaNL1B_Straight1, mcposrNL1B_Straight1;
Rotation mcrotaNL1B_Straight1, mcrotrNL1B_Straight1;
Coords mcposaNL1B_Curved_Entrance_Window, mcposrNL1B_Curved_Entrance_Window;
Rotation mcrotaNL1B_Curved_Entrance_Window, mcrotrNL1B_Curved_Entrance_Window;
Coords mcposaNL1B_Curved_Guide_1, mcposrNL1B_Curved_Guide_1;
Rotation mcrotaNL1B_Curved_Guide_1, mcrotrNL1B_Curved_Guide_1;
Coords mcposaNL1B_Velocity_Selector_Gap_Entrance_Window, mcposrNL1B_Velocity_Selector_Gap_Entrance_Window;
Rotation mcrotaNL1B_Velocity_Selector_Gap_Entrance_Window, mcrotrNL1B_Velocity_Selector_Gap_Entrance_Window;
Coords mcposaBefore_selec, mcposrBefore_selec;
Rotation mcrotaBefore_selec, mcrotrBefore_selec;
Coords mcposaSELEC, mcposrSELEC;
Rotation mcrotaSELEC, mcrotrSELEC;
Coords mcposaAfter_selec, mcposrAfter_selec;
Rotation mcrotaAfter_selec, mcrotrAfter_selec;
Coords mcposaNL1B_Curved_2_Entrance_Window, mcposrNL1B_Curved_2_Entrance_Window;
Rotation mcrotaNL1B_Curved_2_Entrance_Window, mcrotrNL1B_Curved_2_Entrance_Window;
Coords mcposaNL1B_Curved_Guide_2, mcposrNL1B_Curved_Guide_2;
Rotation mcrotaNL1B_Curved_Guide_2, mcrotrNL1B_Curved_Guide_2;
Coords mcposaNL1B_Straight2_Entrance_Window, mcposrNL1B_Straight2_Entrance_Window;
Rotation mcrotaNL1B_Straight2_Entrance_Window, mcrotrNL1B_Straight2_Entrance_Window;
Coords mcposaNL1B_Straight2, mcposrNL1B_Straight2;
Rotation mcrotaNL1B_Straight2, mcrotrNL1B_Straight2;
Coords mcposaNL1B_Elliptical_Entrance_Window, mcposrNL1B_Elliptical_Entrance_Window;
Rotation mcrotaNL1B_Elliptical_Entrance_Window, mcrotrNL1B_Elliptical_Entrance_Window;
Coords mcposaelliptical_piece, mcposrelliptical_piece;
Rotation mcrotaelliptical_piece, mcrotrelliptical_piece;
Coords mcposaVirtual_source, mcposrVirtual_source;
Rotation mcrotaVirtual_source, mcrotrVirtual_source;
Coords mcposaNL1B_Vertical_Guide_Entrance_Window, mcposrNL1B_Vertical_Guide_Entrance_Window;
Rotation mcrotaNL1B_Vertical_Guide_Entrance_Window, mcrotrNL1B_Vertical_Guide_Entrance_Window;
Coords mcposaNL1B_Vertical_Guide, mcposrNL1B_Vertical_Guide;
Rotation mcrotaNL1B_Vertical_Guide, mcrotrNL1B_Vertical_Guide;
Coords mcposaNL1B_Guide_Exit, mcposrNL1B_Guide_Exit;
Rotation mcrotaNL1B_Guide_Exit, mcrotrNL1B_Guide_Exit;
Coords mcposaenergy_endguide, mcposrenergy_endguide;
Rotation mcrotaenergy_endguide, mcrotrenergy_endguide;
Coords mcposaMono_center, mcposrMono_center;
Rotation mcrotaMono_center, mcrotrMono_center;
Coords mcposaMonochromator, mcposrMonochromator;
Rotation mcrotaMonochromator, mcrotrMonochromator;
Coords mcposaMono_sample_arm, mcposrMono_sample_arm;
Rotation mcrotaMono_sample_arm, mcrotrMono_sample_arm;
Coords mcposaenergy_mono, mcposrenergy_mono;
Rotation mcrotaenergy_mono, mcrotrenergy_mono;
Coords mcposaenergy_pre_sample, mcposrenergy_pre_sample;
Rotation mcrotaenergy_pre_sample, mcrotrenergy_pre_sample;
Coords mcposaSample_center, mcposrSample_center;
Rotation mcrotaSample_center, mcrotrSample_center;
Coords mcposaSample_analyser_arm, mcposrSample_analyser_arm;
Rotation mcrotaSample_analyser_arm, mcrotrSample_analyser_arm;
Coords mcposadiv_mono, mcposrdiv_mono;
Rotation mcrotadiv_mono, mcrotrdiv_mono;
Coords mcposadiv_mono_H, mcposrdiv_mono_H;
Rotation mcrotadiv_mono_H, mcrotrdiv_mono_H;
Coords mcposapsd_sam, mcposrpsd_sam;
Rotation mcrotapsd_sam, mcrotrpsd_sam;
Coords mcposa1dpsd, mcposr1dpsd;
Rotation mcrota1dpsd, mcrotr1dpsd;

MCNUM mcnx, mcny, mcnz, mcnvx, mcnvy, mcnvz, mcnt, mcnsx, mcnsy, mcnsz, mcnp;

/* end declare */

void mcinit(void) {
#define mccompcurname  HZB_FLEX
#define mccompcurtype  INSTRUMENT
#define mccompcurindex 0
#define mcposaHZB_FLEX coords_set(0,0,0)
#define kI mcipkI
#define kF mcipkF
#define wVS mcipwVS
#define tilt mciptilt
#define SA mcipSA
#define A3 mcipA3
#define A4 mcipA4
#define L3 mcipL3
#define L4 mcipL4
#define Mono_flatswitch mcipMono_flatswitch
#line 92 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
{
#include <math.h>

int    SM=-1;      // monochromator scattering sense is always negative
int    SS=+1;
// int    SA=-1;
double DM = 3.355; // PG 002 monochromator d-spacing in inv AA
double DA = 3.355; // PG 002 analyser d-spacing in inv AA

/* -------calculate spectrometer angles and monochromator/analyser curvatures-------*/
// ----- monochromator
  A1 = SM*RAD2DEG*asin(PI/(DM*kI));
  A2 = 2*A1;

//  RMH = L2/sin(A1*DEG2RAD);
  RMH = (L1+L2)/(2*sin(A1*DEG2RAD)); //additional factor done by RTP
//  rho_MH=1/RMH;
if (Mono_flatswitch == 1) RMH = 0.0;




//  RMV = 2*L2*sin(A1*DEG2RAD);
  RMV = (L1+L2)*sin(A1*DEG2RAD); //Facotr intropduced by RTP
//  rho_MV=1/RMV;

  printf("monochromator angles are:\n A1=%7.3f, A2=%7.3f [deg]\n",A1,A2);
  printf("radius of curvature RMH=%7.3f [m]\nradius of curvature RMV=%7.3f [m] \n",RMH,RMV);

//  printf("curvature rhoMH=%7.3f [m^-1]\ncurvature rhoMV=%7.3f [m] \n",rho_MH,rho_MV);

// ----- sample
//  A3 = 0;
//  A4 = 70;

printf("sample angles are:\n A3=%7.3f, A4=%7.3f [deg]\n",A3,A4);

// ----- analyser
  A5 = SA*RAD2DEG*asin(PI/(DA*kF));
  A6 = 2*A5;
  RAH = (L3+L4)/(2*sin(A5*DEG2RAD)); //Additional factor multiplied by RTP
  RAV = SA*1/rho_AV;
//  RAV should be (L3+L4)*sin(A5*DEG2RAD) but above equation gives better focusing for these specific conditions of kf=1.55 and the divergence...
  printf("analyser angles are:\n A5=%7.3f, A6=%7.3f [deg]\n",A5,A6);
  printf("radius of curvature RAH=%7.3f [m]\nradius of curvature RAV=%7.3f [m] \n",RAH,RAV);

/* ------------------------------calculate EI, lambdaI and so on ------------------*/
  lambdaI=2*PI/kI;  // in AA
  EI=SI_to_meV*hbar*hbar*kI*kI*1e20/2/m_n;

  SelFreq = 3956.06/(lambdaI*(360/(19.7+(2.1*tilt)))*0.25);
//  empirical formulas for source emmitance to save time on simulations!
//  1st attempt: linear
//  dlambdaI=0.01*(7-lambdaI)*lambdaI;
//  2nd attempt, 1/x, better one:
//  dlambdaI=((0.082/(lambdaI+0.4))-0.002)*lambdaI;
//  3rd attempt, 6th order polynomial, even better (with enough parameters we can model an elephant):
 dlambdaI = (1.75004E-05*pow(lambdaI,6) - 5.13903E-04*pow(lambdaI,5) + 6.04404E-03*pow(lambdaI,4) - 3.66989E-02*pow(lambdaI,3) + 1.22399E-01*pow(lambdaI,2) - 2.19626E-01*lambdaI + 2.0E-01)*lambdaI;

  l_min=lambdaI-dlambdaI;
  l_max=lambdaI+dlambdaI;
  e_min=81.81/l_max/l_max;
  e_max=81.81/l_min/l_min;

  printf("velocity selector frequency is %7.3f [Hz]\n",SelFreq);
  printf("lambda=%7.3f [A]\nkI=%7.3f [A^-1] \nEI=%7.3f [meV]\n",lambdaI,kI,EI);
  printf("Maxwellian source forced to emit only between l_min=%7.3f and l_max=%7.3f [A] in order to save simulation time\n",l_min,l_max);
  printf("correspondingly, the energy monitors at sample position are detecting between %7.3f and %7.3f [meV]\n",e_min,e_max);


// guide properties
  beta_guide_1=LGUIDE_1/2800*180/PI;
  beta_guide_2=LGUIDE_2/2800*180/PI;
  s_guide_1 = 2800-sqrt(2800*2800-LGUIDE_1*LGUIDE_1);
  s_guide_2 = 2800-sqrt(2800*2800-LGUIDE_2*LGUIDE_2);
}
#line 8844 "./HZB_FLEX.c"
#undef Mono_flatswitch
#undef L4
#undef L3
#undef A4
#undef A3
#undef SA
#undef tilt
#undef wVS
#undef kF
#undef kI
#undef mcposaHZB_FLEX
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
#line 39 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  if("NULL") strncpy(mccOrigin_profile, "NULL" ? "NULL" : "", 16384); else mccOrigin_profile[0]='\0';
#line 39 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccOrigin_percent = 10;
#line 39 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccOrigin_flag_save = 0;
#line 39 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccOrigin_minutes = 0;
#line 8880 "./HZB_FLEX.c"

  SIG_MESSAGE("Origin (Init:Place/Rotate)");
  rot_set_rotation(mcrotaOrigin,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 8887 "./HZB_FLEX.c"
  rot_copy(mcrotrOrigin, mcrotaOrigin);
  mcposaOrigin = coords_set(
#line 182 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    0,
#line 182 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    0,
#line 182 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    0);
#line 8896 "./HZB_FLEX.c"
  mctc1 = coords_neg(mcposaOrigin);
  mcposrOrigin = rot_apply(mcrotaOrigin, mctc1);
  mcDEBUG_COMPONENT("Origin", mcposaOrigin, mcrotaOrigin)
  mccomp_posa[1] = mcposaOrigin;
  mccomp_posr[1] = mcposrOrigin;
  mcNCounter[1]  = mcPCounter[1] = mcP2Counter[1] = 0;
  mcAbsorbProp[1]= 0;
    /* Component Gen_Source. */
  /* Setting parameters for component Gen_Source. */
  SIG_MESSAGE("Gen_Source (Init:SetPar)");
#line 129 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  if("NULL") strncpy(mccGen_Source_flux_file, "NULL" ? "NULL" : "", 16384); else mccGen_Source_flux_file[0]='\0';
#line 129 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  if("NULL") strncpy(mccGen_Source_xdiv_file, "NULL" ? "NULL" : "", 16384); else mccGen_Source_xdiv_file[0]='\0';
#line 129 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  if("NULL") strncpy(mccGen_Source_ydiv_file, "NULL" ? "NULL" : "", 16384); else mccGen_Source_ydiv_file[0]='\0';
#line 187 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccGen_Source_radius = 0.0775;
#line 187 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccGen_Source_dist = 1.53;
#line 187 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccGen_Source_focus_xw = 0.2;
#line 187 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccGen_Source_focus_yh = 0.125;
#line 130 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccGen_Source_focus_aw = 0;
#line 130 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccGen_Source_focus_ah = 0;
#line 131 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccGen_Source_E0 = 0;
#line 131 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccGen_Source_dE = 0;
#line 131 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccGen_Source_lambda0 = 0;
#line 131 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccGen_Source_dlambda = 0;
#line 188 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccGen_Source_I1 = 1e10;
#line 132 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccGen_Source_yheight = 0.1;
#line 132 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccGen_Source_xwidth = 0.1;
#line 132 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccGen_Source_verbose = 0;
#line 189 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccGen_Source_T1 = 45;
#line 133 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccGen_Source_flux_file_perAA = 0;
#line 133 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccGen_Source_flux_file_log = 0;
#line 188 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccGen_Source_Lmin = l_min;
#line 188 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccGen_Source_Lmax = l_max;
#line 134 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccGen_Source_Emin = 0;
#line 134 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccGen_Source_Emax = 0;
#line 134 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccGen_Source_T2 = 0;
#line 134 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccGen_Source_I2 = 0;
#line 134 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccGen_Source_T3 = 0;
#line 134 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccGen_Source_I3 = 0;
#line 134 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccGen_Source_zdepth = 0;
#line 134 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccGen_Source_target_index = + 1;
#line 8967 "./HZB_FLEX.c"

  SIG_MESSAGE("Gen_Source (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 191 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    (0)*DEG2RAD,
#line 191 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    (0)*DEG2RAD,
#line 191 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    (0)*DEG2RAD);
#line 8977 "./HZB_FLEX.c"
  rot_mul(mctr1, mcrotaOrigin, mcrotaGen_Source);
  rot_transpose(mcrotaOrigin, mctr1);
  rot_mul(mcrotaGen_Source, mctr1, mcrotrGen_Source);
  mctc1 = coords_set(
#line 190 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    0,
#line 190 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    0,
#line 190 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    0);
#line 8988 "./HZB_FLEX.c"
  rot_transpose(mcrotaOrigin, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaGen_Source = coords_add(mcposaOrigin, mctc2);
  mctc1 = coords_sub(mcposaOrigin, mcposaGen_Source);
  mcposrGen_Source = rot_apply(mcrotaGen_Source, mctc1);
  mcDEBUG_COMPONENT("Gen_Source", mcposaGen_Source, mcrotaGen_Source)
  mccomp_posa[2] = mcposaGen_Source;
  mccomp_posr[2] = mcposrGen_Source;
  mcNCounter[2]  = mcPCounter[2] = mcP2Counter[2] = 0;
  mcAbsorbProp[2]= 0;
    /* Component NL1A_NL1B_NL2_NL3_InPile_Entrance_Window. */
  /* Setting parameters for component NL1A_NL1B_NL2_NL3_InPile_Entrance_Window. */
  SIG_MESSAGE("NL1A_NL1B_NL2_NL3_InPile_Entrance_Window (Init:SetPar)");

  SIG_MESSAGE("NL1A_NL1B_NL2_NL3_InPile_Entrance_Window (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 197 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    (0)*DEG2RAD,
#line 197 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    (0)*DEG2RAD,
#line 197 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    (0)*DEG2RAD);
#line 9011 "./HZB_FLEX.c"
  rot_mul(mctr1, mcrotaOrigin, mcrotaNL1A_NL1B_NL2_NL3_InPile_Entrance_Window);
  rot_transpose(mcrotaGen_Source, mctr1);
  rot_mul(mcrotaNL1A_NL1B_NL2_NL3_InPile_Entrance_Window, mctr1, mcrotrNL1A_NL1B_NL2_NL3_InPile_Entrance_Window);
  mctc1 = coords_set(
#line 197 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    0,
#line 197 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    0,
#line 197 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    1.530);
#line 9022 "./HZB_FLEX.c"
  rot_transpose(mcrotaOrigin, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaNL1A_NL1B_NL2_NL3_InPile_Entrance_Window = coords_add(mcposaOrigin, mctc2);
  mctc1 = coords_sub(mcposaGen_Source, mcposaNL1A_NL1B_NL2_NL3_InPile_Entrance_Window);
  mcposrNL1A_NL1B_NL2_NL3_InPile_Entrance_Window = rot_apply(mcrotaNL1A_NL1B_NL2_NL3_InPile_Entrance_Window, mctc1);
  mcDEBUG_COMPONENT("NL1A_NL1B_NL2_NL3_InPile_Entrance_Window", mcposaNL1A_NL1B_NL2_NL3_InPile_Entrance_Window, mcrotaNL1A_NL1B_NL2_NL3_InPile_Entrance_Window)
  mccomp_posa[3] = mcposaNL1A_NL1B_NL2_NL3_InPile_Entrance_Window;
  mccomp_posr[3] = mcposrNL1A_NL1B_NL2_NL3_InPile_Entrance_Window;
  mcNCounter[3]  = mcPCounter[3] = mcP2Counter[3] = 0;
  mcAbsorbProp[3]= 0;
    /* Component NL1A_NL1B_NL2_NL3_InPile. */
  /* Setting parameters for component NL1A_NL1B_NL2_NL3_InPile. */
  SIG_MESSAGE("NL1A_NL1B_NL2_NL3_InPile (Init:SetPar)");
#line 58 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  if(0) strncpy(mccNL1A_NL1B_NL2_NL3_InPile_reflect, 0 ? 0 : "", 16384); else mccNL1A_NL1B_NL2_NL3_InPile_reflect[0]='\0';
#line 201 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccNL1A_NL1B_NL2_NL3_InPile_w1 = 0.170696;
#line 201 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccNL1A_NL1B_NL2_NL3_InPile_h1 = 0.129;
#line 201 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccNL1A_NL1B_NL2_NL3_InPile_w2 = 0.340522;
#line 201 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccNL1A_NL1B_NL2_NL3_InPile_h2 = 0.129;
#line 201 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccNL1A_NL1B_NL2_NL3_InPile_l = 1.870;
#line 202 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccNL1A_NL1B_NL2_NL3_InPile_R0 = R0_para;
#line 202 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccNL1A_NL1B_NL2_NL3_InPile_Qc = Qc_para;
#line 202 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccNL1A_NL1B_NL2_NL3_InPile_alpha = alpha_para;
#line 202 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccNL1A_NL1B_NL2_NL3_InPile_m = MGUIDE;
#line 202 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccNL1A_NL1B_NL2_NL3_InPile_W = W_para;
#line 9058 "./HZB_FLEX.c"

  SIG_MESSAGE("NL1A_NL1B_NL2_NL3_InPile (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 203 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    (0)*DEG2RAD,
#line 203 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    (0)*DEG2RAD,
#line 203 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    (0)*DEG2RAD);
#line 9068 "./HZB_FLEX.c"
  rot_mul(mctr1, mcrotaNL1A_NL1B_NL2_NL3_InPile_Entrance_Window, mcrotaNL1A_NL1B_NL2_NL3_InPile);
  rot_transpose(mcrotaNL1A_NL1B_NL2_NL3_InPile_Entrance_Window, mctr1);
  rot_mul(mcrotaNL1A_NL1B_NL2_NL3_InPile, mctr1, mcrotrNL1A_NL1B_NL2_NL3_InPile);
  mctc1 = coords_set(
#line 203 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    0,
#line 203 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    0,
#line 203 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    0);
#line 9079 "./HZB_FLEX.c"
  rot_transpose(mcrotaNL1A_NL1B_NL2_NL3_InPile_Entrance_Window, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaNL1A_NL1B_NL2_NL3_InPile = coords_add(mcposaNL1A_NL1B_NL2_NL3_InPile_Entrance_Window, mctc2);
  mctc1 = coords_sub(mcposaNL1A_NL1B_NL2_NL3_InPile_Entrance_Window, mcposaNL1A_NL1B_NL2_NL3_InPile);
  mcposrNL1A_NL1B_NL2_NL3_InPile = rot_apply(mcrotaNL1A_NL1B_NL2_NL3_InPile, mctc1);
  mcDEBUG_COMPONENT("NL1A_NL1B_NL2_NL3_InPile", mcposaNL1A_NL1B_NL2_NL3_InPile, mcrotaNL1A_NL1B_NL2_NL3_InPile)
  mccomp_posa[4] = mcposaNL1A_NL1B_NL2_NL3_InPile;
  mccomp_posr[4] = mcposrNL1A_NL1B_NL2_NL3_InPile;
  mcNCounter[4]  = mcPCounter[4] = mcP2Counter[4] = 0;
  mcAbsorbProp[4]= 0;
    /* Component NL1B_Straight1_Entrance_Window. */
  /* Setting parameters for component NL1B_Straight1_Entrance_Window. */
  SIG_MESSAGE("NL1B_Straight1_Entrance_Window (Init:SetPar)");

  SIG_MESSAGE("NL1B_Straight1_Entrance_Window (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 207 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    (0)*DEG2RAD,
#line 207 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    (-2.15)*DEG2RAD,
#line 207 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    (0)*DEG2RAD);
#line 9102 "./HZB_FLEX.c"
  rot_mul(mctr1, mcrotaOrigin, mcrotaNL1B_Straight1_Entrance_Window);
  rot_transpose(mcrotaNL1A_NL1B_NL2_NL3_InPile, mctr1);
  rot_mul(mcrotaNL1B_Straight1_Entrance_Window, mctr1, mcrotrNL1B_Straight1_Entrance_Window);
  mctc1 = coords_set(
#line 207 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    -0.1277,
#line 207 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    0,
#line 207 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    3.400);
#line 9113 "./HZB_FLEX.c"
  rot_transpose(mcrotaOrigin, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaNL1B_Straight1_Entrance_Window = coords_add(mcposaOrigin, mctc2);
  mctc1 = coords_sub(mcposaNL1A_NL1B_NL2_NL3_InPile, mcposaNL1B_Straight1_Entrance_Window);
  mcposrNL1B_Straight1_Entrance_Window = rot_apply(mcrotaNL1B_Straight1_Entrance_Window, mctc1);
  mcDEBUG_COMPONENT("NL1B_Straight1_Entrance_Window", mcposaNL1B_Straight1_Entrance_Window, mcrotaNL1B_Straight1_Entrance_Window)
  mccomp_posa[5] = mcposaNL1B_Straight1_Entrance_Window;
  mccomp_posr[5] = mcposrNL1B_Straight1_Entrance_Window;
  mcNCounter[5]  = mcPCounter[5] = mcP2Counter[5] = 0;
  mcAbsorbProp[5]= 0;
    /* Component NL1B_Straight1. */
  /* Setting parameters for component NL1B_Straight1. */
  SIG_MESSAGE("NL1B_Straight1 (Init:SetPar)");
#line 58 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  if(0) strncpy(mccNL1B_Straight1_reflect, 0 ? 0 : "", 16384); else mccNL1B_Straight1_reflect[0]='\0';
#line 211 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccNL1B_Straight1_w1 = 0.06;
#line 211 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccNL1B_Straight1_h1 = 0.125;
#line 211 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccNL1B_Straight1_w2 = 0.06;
#line 211 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccNL1B_Straight1_h2 = 0.125;
#line 211 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccNL1B_Straight1_l = 1.549;
#line 212 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccNL1B_Straight1_R0 = R0_para;
#line 212 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccNL1B_Straight1_Qc = Qc_para;
#line 212 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccNL1B_Straight1_alpha = alpha_para;
#line 212 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccNL1B_Straight1_m = MGUIDE;
#line 212 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccNL1B_Straight1_W = W_para;
#line 9149 "./HZB_FLEX.c"

  SIG_MESSAGE("NL1B_Straight1 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 213 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    (0)*DEG2RAD,
#line 213 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    (0)*DEG2RAD,
#line 213 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    (0)*DEG2RAD);
#line 9159 "./HZB_FLEX.c"
  rot_mul(mctr1, mcrotaNL1B_Straight1_Entrance_Window, mcrotaNL1B_Straight1);
  rot_transpose(mcrotaNL1B_Straight1_Entrance_Window, mctr1);
  rot_mul(mcrotaNL1B_Straight1, mctr1, mcrotrNL1B_Straight1);
  mctc1 = coords_set(
#line 213 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    0,
#line 213 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    0,
#line 213 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    0);
#line 9170 "./HZB_FLEX.c"
  rot_transpose(mcrotaNL1B_Straight1_Entrance_Window, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaNL1B_Straight1 = coords_add(mcposaNL1B_Straight1_Entrance_Window, mctc2);
  mctc1 = coords_sub(mcposaNL1B_Straight1_Entrance_Window, mcposaNL1B_Straight1);
  mcposrNL1B_Straight1 = rot_apply(mcrotaNL1B_Straight1, mctc1);
  mcDEBUG_COMPONENT("NL1B_Straight1", mcposaNL1B_Straight1, mcrotaNL1B_Straight1)
  mccomp_posa[6] = mcposaNL1B_Straight1;
  mccomp_posr[6] = mcposrNL1B_Straight1;
  mcNCounter[6]  = mcPCounter[6] = mcP2Counter[6] = 0;
  mcAbsorbProp[6]= 0;
    /* Component NL1B_Curved_Entrance_Window. */
  /* Setting parameters for component NL1B_Curved_Entrance_Window. */
  SIG_MESSAGE("NL1B_Curved_Entrance_Window (Init:SetPar)");

  SIG_MESSAGE("NL1B_Curved_Entrance_Window (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 216 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    (0)*DEG2RAD,
#line 216 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    (0)*DEG2RAD,
#line 216 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    (0)*DEG2RAD);
#line 9193 "./HZB_FLEX.c"
  rot_mul(mctr1, mcrotaNL1B_Straight1_Entrance_Window, mcrotaNL1B_Curved_Entrance_Window);
  rot_transpose(mcrotaNL1B_Straight1, mctr1);
  rot_mul(mcrotaNL1B_Curved_Entrance_Window, mctr1, mcrotrNL1B_Curved_Entrance_Window);
  mctc1 = coords_set(
#line 216 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    0,
#line 216 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    0,
#line 216 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    1.549);
#line 9204 "./HZB_FLEX.c"
  rot_transpose(mcrotaNL1B_Straight1_Entrance_Window, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaNL1B_Curved_Entrance_Window = coords_add(mcposaNL1B_Straight1_Entrance_Window, mctc2);
  mctc1 = coords_sub(mcposaNL1B_Straight1, mcposaNL1B_Curved_Entrance_Window);
  mcposrNL1B_Curved_Entrance_Window = rot_apply(mcrotaNL1B_Curved_Entrance_Window, mctc1);
  mcDEBUG_COMPONENT("NL1B_Curved_Entrance_Window", mcposaNL1B_Curved_Entrance_Window, mcrotaNL1B_Curved_Entrance_Window)
  mccomp_posa[7] = mcposaNL1B_Curved_Entrance_Window;
  mccomp_posr[7] = mcposrNL1B_Curved_Entrance_Window;
  mcNCounter[7]  = mcPCounter[7] = mcP2Counter[7] = 0;
  mcAbsorbProp[7]= 0;
    /* Component NL1B_Curved_Guide_1. */
  /* Setting parameters for component NL1B_Curved_Guide_1. */
  SIG_MESSAGE("NL1B_Curved_Guide_1 (Init:SetPar)");
#line 222 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccNL1B_Curved_Guide_1_w1 = 0.06;
#line 222 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccNL1B_Curved_Guide_1_h1 = 0.125;
#line 222 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccNL1B_Curved_Guide_1_l = LGUIDE_1;
#line 222 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccNL1B_Curved_Guide_1_R0 = R0_para;
#line 222 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccNL1B_Curved_Guide_1_Qc = Qc_para;
#line 223 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccNL1B_Curved_Guide_1_alpha = alpha_para;
#line 223 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccNL1B_Curved_Guide_1_m = MGUIDE;
#line 223 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccNL1B_Curved_Guide_1_W = W_para;
#line 223 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccNL1B_Curved_Guide_1_curvature = 2800;
#line 9236 "./HZB_FLEX.c"

  SIG_MESSAGE("NL1B_Curved_Guide_1 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 224 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    (0)*DEG2RAD,
#line 224 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    (0)*DEG2RAD,
#line 224 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    (180)*DEG2RAD);
#line 9246 "./HZB_FLEX.c"
  rot_mul(mctr1, mcrotaNL1B_Curved_Entrance_Window, mcrotaNL1B_Curved_Guide_1);
  rot_transpose(mcrotaNL1B_Curved_Entrance_Window, mctr1);
  rot_mul(mcrotaNL1B_Curved_Guide_1, mctr1, mcrotrNL1B_Curved_Guide_1);
  mctc1 = coords_set(
#line 224 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    0,
#line 224 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    0,
#line 224 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    0);
#line 9257 "./HZB_FLEX.c"
  rot_transpose(mcrotaNL1B_Curved_Entrance_Window, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaNL1B_Curved_Guide_1 = coords_add(mcposaNL1B_Curved_Entrance_Window, mctc2);
  mctc1 = coords_sub(mcposaNL1B_Curved_Entrance_Window, mcposaNL1B_Curved_Guide_1);
  mcposrNL1B_Curved_Guide_1 = rot_apply(mcrotaNL1B_Curved_Guide_1, mctc1);
  mcDEBUG_COMPONENT("NL1B_Curved_Guide_1", mcposaNL1B_Curved_Guide_1, mcrotaNL1B_Curved_Guide_1)
  mccomp_posa[8] = mcposaNL1B_Curved_Guide_1;
  mccomp_posr[8] = mcposrNL1B_Curved_Guide_1;
  mcNCounter[8]  = mcPCounter[8] = mcP2Counter[8] = 0;
  mcAbsorbProp[8]= 0;
    /* Component NL1B_Velocity_Selector_Gap_Entrance_Window. */
  /* Setting parameters for component NL1B_Velocity_Selector_Gap_Entrance_Window. */
  SIG_MESSAGE("NL1B_Velocity_Selector_Gap_Entrance_Window (Init:SetPar)");

  SIG_MESSAGE("NL1B_Velocity_Selector_Gap_Entrance_Window (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 229 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    (0)*DEG2RAD,
#line 229 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    (- beta_guide_1)*DEG2RAD,
#line 229 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    (0)*DEG2RAD);
#line 9280 "./HZB_FLEX.c"
  rot_mul(mctr1, mcrotaNL1B_Curved_Entrance_Window, mcrotaNL1B_Velocity_Selector_Gap_Entrance_Window);
  rot_transpose(mcrotaNL1B_Curved_Guide_1, mctr1);
  rot_mul(mcrotaNL1B_Velocity_Selector_Gap_Entrance_Window, mctr1, mcrotrNL1B_Velocity_Selector_Gap_Entrance_Window);
  mctc1 = coords_set(
#line 229 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    - s_guide_1,
#line 229 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    0,
#line 229 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    LGUIDE_1);
#line 9291 "./HZB_FLEX.c"
  rot_transpose(mcrotaNL1B_Curved_Entrance_Window, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaNL1B_Velocity_Selector_Gap_Entrance_Window = coords_add(mcposaNL1B_Curved_Entrance_Window, mctc2);
  mctc1 = coords_sub(mcposaNL1B_Curved_Guide_1, mcposaNL1B_Velocity_Selector_Gap_Entrance_Window);
  mcposrNL1B_Velocity_Selector_Gap_Entrance_Window = rot_apply(mcrotaNL1B_Velocity_Selector_Gap_Entrance_Window, mctc1);
  mcDEBUG_COMPONENT("NL1B_Velocity_Selector_Gap_Entrance_Window", mcposaNL1B_Velocity_Selector_Gap_Entrance_Window, mcrotaNL1B_Velocity_Selector_Gap_Entrance_Window)
  mccomp_posa[9] = mcposaNL1B_Velocity_Selector_Gap_Entrance_Window;
  mccomp_posr[9] = mcposrNL1B_Velocity_Selector_Gap_Entrance_Window;
  mcNCounter[9]  = mcPCounter[9] = mcP2Counter[9] = 0;
  mcAbsorbProp[9]= 0;
    /* Component Before_selec. */
  /* Setting parameters for component Before_selec. */
  SIG_MESSAGE("Before_selec (Init:SetPar)");
#line 234 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccBefore_selec_nx = 50;
#line 234 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccBefore_selec_ny = 50;
#line 235 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  if("PSD_monitor_before.dat") strncpy(mccBefore_selec_filename, "PSD_monitor_before.dat" ? "PSD_monitor_before.dat" : "", 16384); else mccBefore_selec_filename[0]='\0';
#line 50 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccBefore_selec_xmin = -0.05;
#line 50 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccBefore_selec_xmax = 0.05;
#line 50 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccBefore_selec_ymin = -0.05;
#line 50 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccBefore_selec_ymax = 0.05;
#line 234 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccBefore_selec_xwidth = 0.20;
#line 235 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccBefore_selec_yheight = 0.20;
#line 234 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccBefore_selec_restore_neutron = 1;
#line 9325 "./HZB_FLEX.c"

  SIG_MESSAGE("Before_selec (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 9332 "./HZB_FLEX.c"
  rot_mul(mctr1, mcrotaNL1B_Velocity_Selector_Gap_Entrance_Window, mcrotaBefore_selec);
  rot_transpose(mcrotaNL1B_Velocity_Selector_Gap_Entrance_Window, mctr1);
  rot_mul(mcrotaBefore_selec, mctr1, mcrotrBefore_selec);
  mctc1 = coords_set(
#line 236 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    0,
#line 236 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    0,
#line 236 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    0.19);
#line 9343 "./HZB_FLEX.c"
  rot_transpose(mcrotaNL1B_Velocity_Selector_Gap_Entrance_Window, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaBefore_selec = coords_add(mcposaNL1B_Velocity_Selector_Gap_Entrance_Window, mctc2);
  mctc1 = coords_sub(mcposaNL1B_Velocity_Selector_Gap_Entrance_Window, mcposaBefore_selec);
  mcposrBefore_selec = rot_apply(mcrotaBefore_selec, mctc1);
  mcDEBUG_COMPONENT("Before_selec", mcposaBefore_selec, mcrotaBefore_selec)
  mccomp_posa[10] = mcposaBefore_selec;
  mccomp_posr[10] = mcposrBefore_selec;
  mcNCounter[10]  = mcPCounter[10] = mcP2Counter[10] = 0;
  mcAbsorbProp[10]= 0;
    /* Component SELEC. */
  /* Setting parameters for component SELEC. */
  SIG_MESSAGE("SELEC (Init:SetPar)");
#line 239 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccSELEC_xmin = -0.0625;
#line 240 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccSELEC_xmax = 0.0625;
#line 241 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccSELEC_ymin = -0.03;
#line 242 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccSELEC_ymax = 0.03;
#line 243 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccSELEC_length = 0.25;
#line 62 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccSELEC_xwidth = 0;
#line 62 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccSELEC_yheight = 0;
#line 244 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccSELEC_nslit = 72;
#line 245 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccSELEC_d = 0.0004;
#line 246 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccSELEC_radius = 0.123;
#line 247 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccSELEC_alpha = 19.7;
#line 248 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccSELEC_nu = SelFreq;
#line 9381 "./HZB_FLEX.c"

  SIG_MESSAGE("SELEC (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 252 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    (0)*DEG2RAD,
#line 252 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    (mciptilt)*DEG2RAD,
#line 252 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    (90)*DEG2RAD);
#line 9391 "./HZB_FLEX.c"
  rot_mul(mctr1, mcrotaNL1B_Velocity_Selector_Gap_Entrance_Window, mcrotaSELEC);
  rot_transpose(mcrotaBefore_selec, mctr1);
  rot_mul(mcrotaSELEC, mctr1, mcrotrSELEC);
  mctc1 = coords_set(
#line 251 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    0,
#line 251 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    0,
#line 251 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    0.2);
#line 9402 "./HZB_FLEX.c"
  rot_transpose(mcrotaNL1B_Velocity_Selector_Gap_Entrance_Window, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaSELEC = coords_add(mcposaNL1B_Velocity_Selector_Gap_Entrance_Window, mctc2);
  mctc1 = coords_sub(mcposaBefore_selec, mcposaSELEC);
  mcposrSELEC = rot_apply(mcrotaSELEC, mctc1);
  mcDEBUG_COMPONENT("SELEC", mcposaSELEC, mcrotaSELEC)
  mccomp_posa[11] = mcposaSELEC;
  mccomp_posr[11] = mcposrSELEC;
  mcNCounter[11]  = mcPCounter[11] = mcP2Counter[11] = 0;
  mcAbsorbProp[11]= 0;
    /* Component After_selec. */
  /* Setting parameters for component After_selec. */
  SIG_MESSAGE("After_selec (Init:SetPar)");
#line 255 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccAfter_selec_nx = 50;
#line 255 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccAfter_selec_ny = 50;
#line 256 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  if("PSD_monitor_after.dat") strncpy(mccAfter_selec_filename, "PSD_monitor_after.dat" ? "PSD_monitor_after.dat" : "", 16384); else mccAfter_selec_filename[0]='\0';
#line 50 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccAfter_selec_xmin = -0.05;
#line 50 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccAfter_selec_xmax = 0.05;
#line 50 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccAfter_selec_ymin = -0.05;
#line 50 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccAfter_selec_ymax = 0.05;
#line 255 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccAfter_selec_xwidth = 0.20;
#line 256 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccAfter_selec_yheight = 0.20;
#line 255 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccAfter_selec_restore_neutron = 1;
#line 9436 "./HZB_FLEX.c"

  SIG_MESSAGE("After_selec (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 9443 "./HZB_FLEX.c"
  rot_mul(mctr1, mcrotaNL1B_Velocity_Selector_Gap_Entrance_Window, mcrotaAfter_selec);
  rot_transpose(mcrotaSELEC, mctr1);
  rot_mul(mcrotaAfter_selec, mctr1, mcrotrAfter_selec);
  mctc1 = coords_set(
#line 257 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    0,
#line 257 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    0,
#line 257 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    0.46);
#line 9454 "./HZB_FLEX.c"
  rot_transpose(mcrotaNL1B_Velocity_Selector_Gap_Entrance_Window, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaAfter_selec = coords_add(mcposaNL1B_Velocity_Selector_Gap_Entrance_Window, mctc2);
  mctc1 = coords_sub(mcposaSELEC, mcposaAfter_selec);
  mcposrAfter_selec = rot_apply(mcrotaAfter_selec, mctc1);
  mcDEBUG_COMPONENT("After_selec", mcposaAfter_selec, mcrotaAfter_selec)
  mccomp_posa[12] = mcposaAfter_selec;
  mccomp_posr[12] = mcposrAfter_selec;
  mcNCounter[12]  = mcPCounter[12] = mcP2Counter[12] = 0;
  mcAbsorbProp[12]= 0;
    /* Component NL1B_Curved_2_Entrance_Window. */
  /* Setting parameters for component NL1B_Curved_2_Entrance_Window. */
  SIG_MESSAGE("NL1B_Curved_2_Entrance_Window (Init:SetPar)");

  SIG_MESSAGE("NL1B_Curved_2_Entrance_Window (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 260 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    (0)*DEG2RAD,
#line 260 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    (0)*DEG2RAD,
#line 260 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    (0)*DEG2RAD);
#line 9477 "./HZB_FLEX.c"
  rot_mul(mctr1, mcrotaNL1B_Velocity_Selector_Gap_Entrance_Window, mcrotaNL1B_Curved_2_Entrance_Window);
  rot_transpose(mcrotaAfter_selec, mctr1);
  rot_mul(mcrotaNL1B_Curved_2_Entrance_Window, mctr1, mcrotrNL1B_Curved_2_Entrance_Window);
  mctc1 = coords_set(
#line 260 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    0,
#line 260 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    0,
#line 260 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    0.47);
#line 9488 "./HZB_FLEX.c"
  rot_transpose(mcrotaNL1B_Velocity_Selector_Gap_Entrance_Window, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaNL1B_Curved_2_Entrance_Window = coords_add(mcposaNL1B_Velocity_Selector_Gap_Entrance_Window, mctc2);
  mctc1 = coords_sub(mcposaAfter_selec, mcposaNL1B_Curved_2_Entrance_Window);
  mcposrNL1B_Curved_2_Entrance_Window = rot_apply(mcrotaNL1B_Curved_2_Entrance_Window, mctc1);
  mcDEBUG_COMPONENT("NL1B_Curved_2_Entrance_Window", mcposaNL1B_Curved_2_Entrance_Window, mcrotaNL1B_Curved_2_Entrance_Window)
  mccomp_posa[13] = mcposaNL1B_Curved_2_Entrance_Window;
  mccomp_posr[13] = mcposrNL1B_Curved_2_Entrance_Window;
  mcNCounter[13]  = mcPCounter[13] = mcP2Counter[13] = 0;
  mcAbsorbProp[13]= 0;
    /* Component NL1B_Curved_Guide_2. */
  /* Setting parameters for component NL1B_Curved_Guide_2. */
  SIG_MESSAGE("NL1B_Curved_Guide_2 (Init:SetPar)");
#line 263 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccNL1B_Curved_Guide_2_w1 = 0.06;
#line 263 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccNL1B_Curved_Guide_2_h1 = 0.125;
#line 263 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccNL1B_Curved_Guide_2_l = LGUIDE_2;
#line 263 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccNL1B_Curved_Guide_2_R0 = R0_para;
#line 263 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccNL1B_Curved_Guide_2_Qc = Qc_para;
#line 264 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccNL1B_Curved_Guide_2_alpha = alpha_para;
#line 264 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccNL1B_Curved_Guide_2_m = MGUIDE;
#line 264 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccNL1B_Curved_Guide_2_W = W_para;
#line 264 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccNL1B_Curved_Guide_2_curvature = 2800;
#line 9520 "./HZB_FLEX.c"

  SIG_MESSAGE("NL1B_Curved_Guide_2 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 265 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    (0)*DEG2RAD,
#line 265 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    (0)*DEG2RAD,
#line 265 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    (180)*DEG2RAD);
#line 9530 "./HZB_FLEX.c"
  rot_mul(mctr1, mcrotaNL1B_Curved_2_Entrance_Window, mcrotaNL1B_Curved_Guide_2);
  rot_transpose(mcrotaNL1B_Curved_2_Entrance_Window, mctr1);
  rot_mul(mcrotaNL1B_Curved_Guide_2, mctr1, mcrotrNL1B_Curved_Guide_2);
  mctc1 = coords_set(
#line 265 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    0,
#line 265 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    0,
#line 265 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    0);
#line 9541 "./HZB_FLEX.c"
  rot_transpose(mcrotaNL1B_Curved_2_Entrance_Window, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaNL1B_Curved_Guide_2 = coords_add(mcposaNL1B_Curved_2_Entrance_Window, mctc2);
  mctc1 = coords_sub(mcposaNL1B_Curved_2_Entrance_Window, mcposaNL1B_Curved_Guide_2);
  mcposrNL1B_Curved_Guide_2 = rot_apply(mcrotaNL1B_Curved_Guide_2, mctc1);
  mcDEBUG_COMPONENT("NL1B_Curved_Guide_2", mcposaNL1B_Curved_Guide_2, mcrotaNL1B_Curved_Guide_2)
  mccomp_posa[14] = mcposaNL1B_Curved_Guide_2;
  mccomp_posr[14] = mcposrNL1B_Curved_Guide_2;
  mcNCounter[14]  = mcPCounter[14] = mcP2Counter[14] = 0;
  mcAbsorbProp[14]= 0;
    /* Component NL1B_Straight2_Entrance_Window. */
  /* Setting parameters for component NL1B_Straight2_Entrance_Window. */
  SIG_MESSAGE("NL1B_Straight2_Entrance_Window (Init:SetPar)");

  SIG_MESSAGE("NL1B_Straight2_Entrance_Window (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 268 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    (0)*DEG2RAD,
#line 268 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    (- beta_guide_2)*DEG2RAD,
#line 268 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    (0)*DEG2RAD);
#line 9564 "./HZB_FLEX.c"
  rot_mul(mctr1, mcrotaNL1B_Curved_2_Entrance_Window, mcrotaNL1B_Straight2_Entrance_Window);
  rot_transpose(mcrotaNL1B_Curved_Guide_2, mctr1);
  rot_mul(mcrotaNL1B_Straight2_Entrance_Window, mctr1, mcrotrNL1B_Straight2_Entrance_Window);
  mctc1 = coords_set(
#line 268 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    - s_guide_2,
#line 268 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    0,
#line 268 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    LGUIDE_2);
#line 9575 "./HZB_FLEX.c"
  rot_transpose(mcrotaNL1B_Curved_2_Entrance_Window, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaNL1B_Straight2_Entrance_Window = coords_add(mcposaNL1B_Curved_2_Entrance_Window, mctc2);
  mctc1 = coords_sub(mcposaNL1B_Curved_Guide_2, mcposaNL1B_Straight2_Entrance_Window);
  mcposrNL1B_Straight2_Entrance_Window = rot_apply(mcrotaNL1B_Straight2_Entrance_Window, mctc1);
  mcDEBUG_COMPONENT("NL1B_Straight2_Entrance_Window", mcposaNL1B_Straight2_Entrance_Window, mcrotaNL1B_Straight2_Entrance_Window)
  mccomp_posa[15] = mcposaNL1B_Straight2_Entrance_Window;
  mccomp_posr[15] = mcposrNL1B_Straight2_Entrance_Window;
  mcNCounter[15]  = mcPCounter[15] = mcP2Counter[15] = 0;
  mcAbsorbProp[15]= 0;
    /* Component NL1B_Straight2. */
  /* Setting parameters for component NL1B_Straight2. */
  SIG_MESSAGE("NL1B_Straight2 (Init:SetPar)");
#line 58 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  if(0) strncpy(mccNL1B_Straight2_reflect, 0 ? 0 : "", 16384); else mccNL1B_Straight2_reflect[0]='\0';
#line 271 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccNL1B_Straight2_w1 = 0.06;
#line 271 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccNL1B_Straight2_h1 = 0.125;
#line 271 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccNL1B_Straight2_w2 = 0.06;
#line 271 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccNL1B_Straight2_h2 = 0.125;
#line 271 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccNL1B_Straight2_l = 3.500;
#line 272 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccNL1B_Straight2_R0 = R0_para;
#line 272 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccNL1B_Straight2_Qc = Qc_para;
#line 272 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccNL1B_Straight2_alpha = alpha_para;
#line 272 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccNL1B_Straight2_m = MGUIDE;
#line 272 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccNL1B_Straight2_W = W_para;
#line 9611 "./HZB_FLEX.c"

  SIG_MESSAGE("NL1B_Straight2 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 273 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    (0)*DEG2RAD,
#line 273 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    (0)*DEG2RAD,
#line 273 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    (0)*DEG2RAD);
#line 9621 "./HZB_FLEX.c"
  rot_mul(mctr1, mcrotaNL1B_Straight2_Entrance_Window, mcrotaNL1B_Straight2);
  rot_transpose(mcrotaNL1B_Straight2_Entrance_Window, mctr1);
  rot_mul(mcrotaNL1B_Straight2, mctr1, mcrotrNL1B_Straight2);
  mctc1 = coords_set(
#line 273 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    0,
#line 273 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    0,
#line 273 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    0);
#line 9632 "./HZB_FLEX.c"
  rot_transpose(mcrotaNL1B_Straight2_Entrance_Window, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaNL1B_Straight2 = coords_add(mcposaNL1B_Straight2_Entrance_Window, mctc2);
  mctc1 = coords_sub(mcposaNL1B_Straight2_Entrance_Window, mcposaNL1B_Straight2);
  mcposrNL1B_Straight2 = rot_apply(mcrotaNL1B_Straight2, mctc1);
  mcDEBUG_COMPONENT("NL1B_Straight2", mcposaNL1B_Straight2, mcrotaNL1B_Straight2)
  mccomp_posa[16] = mcposaNL1B_Straight2;
  mccomp_posr[16] = mcposrNL1B_Straight2;
  mcNCounter[16]  = mcPCounter[16] = mcP2Counter[16] = 0;
  mcAbsorbProp[16]= 0;
    /* Component NL1B_Elliptical_Entrance_Window. */
  /* Setting parameters for component NL1B_Elliptical_Entrance_Window. */
  SIG_MESSAGE("NL1B_Elliptical_Entrance_Window (Init:SetPar)");

  SIG_MESSAGE("NL1B_Elliptical_Entrance_Window (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 279 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    (0)*DEG2RAD,
#line 279 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    (0.000)*DEG2RAD,
#line 279 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    (0.00)*DEG2RAD);
#line 9655 "./HZB_FLEX.c"
  rot_mul(mctr1, mcrotaNL1B_Straight2_Entrance_Window, mcrotaNL1B_Elliptical_Entrance_Window);
  rot_transpose(mcrotaNL1B_Straight2, mctr1);
  rot_mul(mcrotaNL1B_Elliptical_Entrance_Window, mctr1, mcrotrNL1B_Elliptical_Entrance_Window);
  mctc1 = coords_set(
#line 279 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    0.0000,
#line 279 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    0,
#line 279 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    3.500);
#line 9666 "./HZB_FLEX.c"
  rot_transpose(mcrotaNL1B_Straight2_Entrance_Window, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaNL1B_Elliptical_Entrance_Window = coords_add(mcposaNL1B_Straight2_Entrance_Window, mctc2);
  mctc1 = coords_sub(mcposaNL1B_Straight2, mcposaNL1B_Elliptical_Entrance_Window);
  mcposrNL1B_Elliptical_Entrance_Window = rot_apply(mcrotaNL1B_Elliptical_Entrance_Window, mctc1);
  mcDEBUG_COMPONENT("NL1B_Elliptical_Entrance_Window", mcposaNL1B_Elliptical_Entrance_Window, mcrotaNL1B_Elliptical_Entrance_Window)
  mccomp_posa[17] = mcposaNL1B_Elliptical_Entrance_Window;
  mccomp_posr[17] = mcposrNL1B_Elliptical_Entrance_Window;
  mcNCounter[17]  = mcPCounter[17] = mcP2Counter[17] = 0;
  mcAbsorbProp[17]= 0;
    /* Component elliptical_piece. */
  /* Setting parameters for component elliptical_piece. */
  SIG_MESSAGE("elliptical_piece (Init:SetPar)");
#line 283 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  if("elliptical") strncpy(mccelliptical_piece_option, "elliptical" ? "elliptical" : "", 16384); else mccelliptical_piece_option[0]='\0';
#line 283 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccelliptical_piece_w1 = 0.06;
#line 283 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccelliptical_piece_h1 = 0.125;
#line 284 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccelliptical_piece_l = 2.5;
#line 284 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccelliptical_piece_linw = 2.9;
#line 284 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccelliptical_piece_loutw = 0.4;
#line 81 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccelliptical_piece_linh = 0;
#line 81 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccelliptical_piece_louth = 0;
#line 284 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccelliptical_piece_R0 = 0.99;
#line 284 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccelliptical_piece_Qcx = 0.0217;
#line 285 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccelliptical_piece_Qcy = 0.0217;
#line 285 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccelliptical_piece_alphax = 4.00;
#line 285 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccelliptical_piece_alphay = alpha_para;
#line 285 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccelliptical_piece_W = 0.002;
#line 286 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccelliptical_piece_mx = 5.2;
#line 286 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccelliptical_piece_my = MGUIDE;
#line 283 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccelliptical_piece_segno = 800;
#line 83 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccelliptical_piece_curvature = 0;
#line 83 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccelliptical_piece_curvature_v = 0;
#line 9718 "./HZB_FLEX.c"

  SIG_MESSAGE("elliptical_piece (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 289 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    (0)*DEG2RAD,
#line 289 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    (0)*DEG2RAD,
#line 289 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    (0)*DEG2RAD);
#line 9728 "./HZB_FLEX.c"
  rot_mul(mctr1, mcrotaNL1B_Elliptical_Entrance_Window, mcrotaelliptical_piece);
  rot_transpose(mcrotaNL1B_Elliptical_Entrance_Window, mctr1);
  rot_mul(mcrotaelliptical_piece, mctr1, mcrotrelliptical_piece);
  mctc1 = coords_set(
#line 288 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    0,
#line 288 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    0,
#line 288 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    0);
#line 9739 "./HZB_FLEX.c"
  rot_transpose(mcrotaNL1B_Elliptical_Entrance_Window, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaelliptical_piece = coords_add(mcposaNL1B_Elliptical_Entrance_Window, mctc2);
  mctc1 = coords_sub(mcposaNL1B_Elliptical_Entrance_Window, mcposaelliptical_piece);
  mcposrelliptical_piece = rot_apply(mcrotaelliptical_piece, mctc1);
  mcDEBUG_COMPONENT("elliptical_piece", mcposaelliptical_piece, mcrotaelliptical_piece)
  mccomp_posa[18] = mcposaelliptical_piece;
  mccomp_posr[18] = mcposrelliptical_piece;
  mcNCounter[18]  = mcPCounter[18] = mcP2Counter[18] = 0;
  mcAbsorbProp[18]= 0;
    /* Component Virtual_source. */
  /* Setting parameters for component Virtual_source. */
  SIG_MESSAGE("Virtual_source (Init:SetPar)");
#line 46 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccVirtual_source_xmin = 0;
#line 46 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccVirtual_source_xmax = 0;
#line 46 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccVirtual_source_ymin = 0;
#line 46 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccVirtual_source_ymax = 0;
#line 46 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccVirtual_source_radius = 0;
#line 295 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccVirtual_source_xwidth = mcipwVS;
#line 295 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccVirtual_source_yheight = 0.2;
#line 9767 "./HZB_FLEX.c"

  SIG_MESSAGE("Virtual_source (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 9774 "./HZB_FLEX.c"
  rot_mul(mctr1, mcrotaNL1B_Elliptical_Entrance_Window, mcrotaVirtual_source);
  rot_transpose(mcrotaelliptical_piece, mctr1);
  rot_mul(mcrotaVirtual_source, mctr1, mcrotrVirtual_source);
  mctc1 = coords_set(
#line 296 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    0,
#line 296 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    0,
#line 296 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    2.900);
#line 9785 "./HZB_FLEX.c"
  rot_transpose(mcrotaNL1B_Elliptical_Entrance_Window, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaVirtual_source = coords_add(mcposaNL1B_Elliptical_Entrance_Window, mctc2);
  mctc1 = coords_sub(mcposaelliptical_piece, mcposaVirtual_source);
  mcposrVirtual_source = rot_apply(mcrotaVirtual_source, mctc1);
  mcDEBUG_COMPONENT("Virtual_source", mcposaVirtual_source, mcrotaVirtual_source)
  mccomp_posa[19] = mcposaVirtual_source;
  mccomp_posr[19] = mcposrVirtual_source;
  mcNCounter[19]  = mcPCounter[19] = mcP2Counter[19] = 0;
  mcAbsorbProp[19]= 0;
    /* Component NL1B_Vertical_Guide_Entrance_Window. */
  /* Setting parameters for component NL1B_Vertical_Guide_Entrance_Window. */
  SIG_MESSAGE("NL1B_Vertical_Guide_Entrance_Window (Init:SetPar)");

  SIG_MESSAGE("NL1B_Vertical_Guide_Entrance_Window (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 303 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    (0)*DEG2RAD,
#line 303 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    (0.000)*DEG2RAD,
#line 303 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    (0.00)*DEG2RAD);
#line 9808 "./HZB_FLEX.c"
  rot_mul(mctr1, mcrotaNL1B_Elliptical_Entrance_Window, mcrotaNL1B_Vertical_Guide_Entrance_Window);
  rot_transpose(mcrotaVirtual_source, mctr1);
  rot_mul(mcrotaNL1B_Vertical_Guide_Entrance_Window, mctr1, mcrotrNL1B_Vertical_Guide_Entrance_Window);
  mctc1 = coords_set(
#line 303 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    0.0000,
#line 303 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    0,
#line 303 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    2.925);
#line 9819 "./HZB_FLEX.c"
  rot_transpose(mcrotaNL1B_Elliptical_Entrance_Window, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaNL1B_Vertical_Guide_Entrance_Window = coords_add(mcposaNL1B_Elliptical_Entrance_Window, mctc2);
  mctc1 = coords_sub(mcposaVirtual_source, mcposaNL1B_Vertical_Guide_Entrance_Window);
  mcposrNL1B_Vertical_Guide_Entrance_Window = rot_apply(mcrotaNL1B_Vertical_Guide_Entrance_Window, mctc1);
  mcDEBUG_COMPONENT("NL1B_Vertical_Guide_Entrance_Window", mcposaNL1B_Vertical_Guide_Entrance_Window, mcrotaNL1B_Vertical_Guide_Entrance_Window)
  mccomp_posa[20] = mcposaNL1B_Vertical_Guide_Entrance_Window;
  mccomp_posr[20] = mcposrNL1B_Vertical_Guide_Entrance_Window;
  mcNCounter[20]  = mcPCounter[20] = mcP2Counter[20] = 0;
  mcAbsorbProp[20]= 0;
    /* Component NL1B_Vertical_Guide. */
  /* Setting parameters for component NL1B_Vertical_Guide. */
  SIG_MESSAGE("NL1B_Vertical_Guide (Init:SetPar)");
#line 307 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccNL1B_Vertical_Guide_w1 = 0.03;
#line 307 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccNL1B_Vertical_Guide_h1 = 0.125;
#line 307 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccNL1B_Vertical_Guide_w2 = 0.15;
#line 307 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccNL1B_Vertical_Guide_h2 = 0.125;
#line 307 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccNL1B_Vertical_Guide_l = 1.475;
#line 308 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccNL1B_Vertical_Guide_R0 = R0_para;
#line 71 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccNL1B_Vertical_Guide_Qc = 0;
#line 71 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccNL1B_Vertical_Guide_alpha = 0;
#line 71 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccNL1B_Vertical_Guide_m = 0;
#line 309 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccNL1B_Vertical_Guide_nslit = 1;
#line 71 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccNL1B_Vertical_Guide_d = 0.0005;
#line 308 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccNL1B_Vertical_Guide_Qcx = Qc_para;
#line 308 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccNL1B_Vertical_Guide_Qcy = Qc_para;
#line 309 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccNL1B_Vertical_Guide_alphax = alpha_para;
#line 309 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccNL1B_Vertical_Guide_alphay = alpha_para;
#line 309 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccNL1B_Vertical_Guide_W = W_para;
#line 310 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccNL1B_Vertical_Guide_mx = 0;
#line 310 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccNL1B_Vertical_Guide_my = MGUIDE;
#line 72 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccNL1B_Vertical_Guide_nu = 0;
#line 72 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccNL1B_Vertical_Guide_phase = 0;
#line 9873 "./HZB_FLEX.c"

  SIG_MESSAGE("NL1B_Vertical_Guide (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 313 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    (0)*DEG2RAD,
#line 313 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    (0)*DEG2RAD,
#line 313 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    (0)*DEG2RAD);
#line 9883 "./HZB_FLEX.c"
  rot_mul(mctr1, mcrotaNL1B_Vertical_Guide_Entrance_Window, mcrotaNL1B_Vertical_Guide);
  rot_transpose(mcrotaNL1B_Vertical_Guide_Entrance_Window, mctr1);
  rot_mul(mcrotaNL1B_Vertical_Guide, mctr1, mcrotrNL1B_Vertical_Guide);
  mctc1 = coords_set(
#line 312 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    0,
#line 312 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    0,
#line 312 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    0);
#line 9894 "./HZB_FLEX.c"
  rot_transpose(mcrotaNL1B_Vertical_Guide_Entrance_Window, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaNL1B_Vertical_Guide = coords_add(mcposaNL1B_Vertical_Guide_Entrance_Window, mctc2);
  mctc1 = coords_sub(mcposaNL1B_Vertical_Guide_Entrance_Window, mcposaNL1B_Vertical_Guide);
  mcposrNL1B_Vertical_Guide = rot_apply(mcrotaNL1B_Vertical_Guide, mctc1);
  mcDEBUG_COMPONENT("NL1B_Vertical_Guide", mcposaNL1B_Vertical_Guide, mcrotaNL1B_Vertical_Guide)
  mccomp_posa[21] = mcposaNL1B_Vertical_Guide;
  mccomp_posr[21] = mcposrNL1B_Vertical_Guide;
  mcNCounter[21]  = mcPCounter[21] = mcP2Counter[21] = 0;
  mcAbsorbProp[21]= 0;
    /* Component NL1B_Guide_Exit. */
  /* Setting parameters for component NL1B_Guide_Exit. */
  SIG_MESSAGE("NL1B_Guide_Exit (Init:SetPar)");

  SIG_MESSAGE("NL1B_Guide_Exit (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 316 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    (0)*DEG2RAD,
#line 316 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    (0)*DEG2RAD,
#line 316 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    (0)*DEG2RAD);
#line 9917 "./HZB_FLEX.c"
  rot_mul(mctr1, mcrotaNL1B_Vertical_Guide_Entrance_Window, mcrotaNL1B_Guide_Exit);
  rot_transpose(mcrotaNL1B_Vertical_Guide, mctr1);
  rot_mul(mcrotaNL1B_Guide_Exit, mctr1, mcrotrNL1B_Guide_Exit);
  mctc1 = coords_set(
#line 316 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    0,
#line 316 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    0,
#line 316 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    1.475);
#line 9928 "./HZB_FLEX.c"
  rot_transpose(mcrotaNL1B_Vertical_Guide_Entrance_Window, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaNL1B_Guide_Exit = coords_add(mcposaNL1B_Vertical_Guide_Entrance_Window, mctc2);
  mctc1 = coords_sub(mcposaNL1B_Vertical_Guide, mcposaNL1B_Guide_Exit);
  mcposrNL1B_Guide_Exit = rot_apply(mcrotaNL1B_Guide_Exit, mctc1);
  mcDEBUG_COMPONENT("NL1B_Guide_Exit", mcposaNL1B_Guide_Exit, mcrotaNL1B_Guide_Exit)
  mccomp_posa[22] = mcposaNL1B_Guide_Exit;
  mccomp_posr[22] = mcposrNL1B_Guide_Exit;
  mcNCounter[22]  = mcPCounter[22] = mcP2Counter[22] = 0;
  mcAbsorbProp[22]= 0;
    /* Component energy_endguide. */
  /* Setting parameters for component energy_endguide. */
  SIG_MESSAGE("energy_endguide (Init:SetPar)");
#line 324 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  if("energy-endguide.dat") strncpy(mccenergy_endguide_filename, "energy-endguide.dat" ? "energy-endguide.dat" : "", 16384); else mccenergy_endguide_filename[0]='\0';
#line 53 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccenergy_endguide_xmin = -0.05;
#line 53 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccenergy_endguide_xmax = 0.05;
#line 53 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccenergy_endguide_ymin = -0.05;
#line 53 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccenergy_endguide_ymax = 0.05;
#line 325 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccenergy_endguide_xwidth = 0.17;
#line 325 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccenergy_endguide_yheight = 0.135;
#line 325 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccenergy_endguide_Emin = e_min;
#line 325 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccenergy_endguide_Emax = e_max;
#line 324 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccenergy_endguide_restore_neutron = 1;
#line 54 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccenergy_endguide_nowritefile = 0;
#line 9964 "./HZB_FLEX.c"

  SIG_MESSAGE("energy_endguide (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 9971 "./HZB_FLEX.c"
  rot_mul(mctr1, mcrotaNL1B_Guide_Exit, mcrotaenergy_endguide);
  rot_transpose(mcrotaNL1B_Guide_Exit, mctr1);
  rot_mul(mcrotaenergy_endguide, mctr1, mcrotrenergy_endguide);
  mctc1 = coords_set(
#line 326 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    0,
#line 326 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    0,
#line 326 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    0);
#line 9982 "./HZB_FLEX.c"
  rot_transpose(mcrotaNL1B_Guide_Exit, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaenergy_endguide = coords_add(mcposaNL1B_Guide_Exit, mctc2);
  mctc1 = coords_sub(mcposaNL1B_Guide_Exit, mcposaenergy_endguide);
  mcposrenergy_endguide = rot_apply(mcrotaenergy_endguide, mctc1);
  mcDEBUG_COMPONENT("energy_endguide", mcposaenergy_endguide, mcrotaenergy_endguide)
  mccomp_posa[23] = mcposaenergy_endguide;
  mccomp_posr[23] = mcposrenergy_endguide;
  mcNCounter[23]  = mcPCounter[23] = mcP2Counter[23] = 0;
  mcAbsorbProp[23]= 0;
    /* Component Mono_center. */
  /* Setting parameters for component Mono_center. */
  SIG_MESSAGE("Mono_center (Init:SetPar)");

  SIG_MESSAGE("Mono_center (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 343 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    (0)*DEG2RAD,
#line 343 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    (A1)*DEG2RAD,
#line 343 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    (0)*DEG2RAD);
#line 10005 "./HZB_FLEX.c"
  rot_mul(mctr1, mcrotaNL1B_Guide_Exit, mcrotaMono_center);
  rot_transpose(mcrotaenergy_endguide, mctr1);
  rot_mul(mcrotaMono_center, mctr1, mcrotrMono_center);
  mctc1 = coords_set(
#line 343 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    0,
#line 343 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    0,
#line 343 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    0.25);
#line 10016 "./HZB_FLEX.c"
  rot_transpose(mcrotaNL1B_Guide_Exit, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaMono_center = coords_add(mcposaNL1B_Guide_Exit, mctc2);
  mctc1 = coords_sub(mcposaenergy_endguide, mcposaMono_center);
  mcposrMono_center = rot_apply(mcrotaMono_center, mctc1);
  mcDEBUG_COMPONENT("Mono_center", mcposaMono_center, mcrotaMono_center)
  mccomp_posa[24] = mcposaMono_center;
  mccomp_posr[24] = mcposrMono_center;
  mcNCounter[24]  = mcPCounter[24] = mcP2Counter[24] = 0;
  mcAbsorbProp[24]= 0;
    /* Component Monochromator. */
  /* Setting parameters for component Monochromator. */
  SIG_MESSAGE("Monochromator (Init:SetPar)");
#line 99 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  if("NULL") strncpy(mccMonochromator_reflect, "NULL" ? "NULL" : "", 16384); else mccMonochromator_reflect[0]='\0';
#line 99 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  if("NULL") strncpy(mccMonochromator_transmit, "NULL" ? "NULL" : "", 16384); else mccMonochromator_transmit[0]='\0';
#line 347 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccMonochromator_zwidth = 0.02;
#line 348 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccMonochromator_yheight = 0.02;
#line 349 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccMonochromator_gap = 0.0005;
#line 350 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccMonochromator_NH = 15;
#line 351 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccMonochromator_NV = 7;
#line 352 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccMonochromator_mosaich = 40;
#line 353 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccMonochromator_mosaicv = 40;
#line 354 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccMonochromator_r0 = 1.0;
#line 355 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccMonochromator_t0 = 0.00001;
#line 101 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccMonochromator_Q = 1.8734;
#line 356 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccMonochromator_RV = RMV;
#line 357 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccMonochromator_RH = RMH;
#line 102 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccMonochromator_DM = 0;
#line 102 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccMonochromator_mosaic = 0;
#line 102 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccMonochromator_width = 0;
#line 102 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccMonochromator_height = 0;
#line 102 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccMonochromator_verbose = 0;
#line 102 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccMonochromator_order = 0;
#line 10070 "./HZB_FLEX.c"

  SIG_MESSAGE("Monochromator (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 359 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    (0)*DEG2RAD,
#line 359 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    (0)*DEG2RAD,
#line 359 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    (0)*DEG2RAD);
#line 10080 "./HZB_FLEX.c"
  rot_mul(mctr1, mcrotaMono_center, mcrotaMonochromator);
  rot_transpose(mcrotaMono_center, mctr1);
  rot_mul(mcrotaMonochromator, mctr1, mcrotrMonochromator);
  mctc1 = coords_set(
#line 359 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    0,
#line 359 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    0,
#line 359 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    0);
#line 10091 "./HZB_FLEX.c"
  rot_transpose(mcrotaMono_center, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaMonochromator = coords_add(mcposaMono_center, mctc2);
  mctc1 = coords_sub(mcposaMono_center, mcposaMonochromator);
  mcposrMonochromator = rot_apply(mcrotaMonochromator, mctc1);
  mcDEBUG_COMPONENT("Monochromator", mcposaMonochromator, mcrotaMonochromator)
  mccomp_posa[25] = mcposaMonochromator;
  mccomp_posr[25] = mcposrMonochromator;
  mcNCounter[25]  = mcPCounter[25] = mcP2Counter[25] = 0;
  mcAbsorbProp[25]= 0;
    /* Component Mono_sample_arm. */
  /* Setting parameters for component Mono_sample_arm. */
  SIG_MESSAGE("Mono_sample_arm (Init:SetPar)");

  SIG_MESSAGE("Mono_sample_arm (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 362 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    (0)*DEG2RAD,
#line 362 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    (A2)*DEG2RAD,
#line 362 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    (0)*DEG2RAD);
#line 10114 "./HZB_FLEX.c"
  rot_mul(mctr1, mcrotaNL1B_Guide_Exit, mcrotaMono_sample_arm);
  rot_transpose(mcrotaMonochromator, mctr1);
  rot_mul(mcrotaMono_sample_arm, mctr1, mcrotrMono_sample_arm);
  mctc1 = coords_set(
#line 362 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    0,
#line 362 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    0,
#line 362 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    0);
#line 10125 "./HZB_FLEX.c"
  rot_transpose(mcrotaMono_center, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaMono_sample_arm = coords_add(mcposaMono_center, mctc2);
  mctc1 = coords_sub(mcposaMonochromator, mcposaMono_sample_arm);
  mcposrMono_sample_arm = rot_apply(mcrotaMono_sample_arm, mctc1);
  mcDEBUG_COMPONENT("Mono_sample_arm", mcposaMono_sample_arm, mcrotaMono_sample_arm)
  mccomp_posa[26] = mcposaMono_sample_arm;
  mccomp_posr[26] = mcposrMono_sample_arm;
  mcNCounter[26]  = mcPCounter[26] = mcP2Counter[26] = 0;
  mcAbsorbProp[26]= 0;
    /* Component energy_mono. */
  /* Setting parameters for component energy_mono. */
  SIG_MESSAGE("energy_mono (Init:SetPar)");
#line 370 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  if("energy-mono.dat") strncpy(mccenergy_mono_filename, "energy-mono.dat" ? "energy-mono.dat" : "", 16384); else mccenergy_mono_filename[0]='\0';
#line 53 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccenergy_mono_xmin = -0.05;
#line 53 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccenergy_mono_xmax = 0.05;
#line 53 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccenergy_mono_ymin = -0.05;
#line 53 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccenergy_mono_ymax = 0.05;
#line 371 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccenergy_mono_xwidth = 0.2;
#line 371 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccenergy_mono_yheight = 0.125;
#line 371 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccenergy_mono_Emin = e_min;
#line 371 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccenergy_mono_Emax = e_max;
#line 370 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccenergy_mono_restore_neutron = 1;
#line 54 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccenergy_mono_nowritefile = 0;
#line 10161 "./HZB_FLEX.c"

  SIG_MESSAGE("energy_mono (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 10168 "./HZB_FLEX.c"
  rot_mul(mctr1, mcrotaMono_sample_arm, mcrotaenergy_mono);
  rot_transpose(mcrotaMono_sample_arm, mctr1);
  rot_mul(mcrotaenergy_mono, mctr1, mcrotrenergy_mono);
  mctc1 = coords_set(
#line 372 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    0,
#line 372 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    0,
#line 372 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    0.15);
#line 10179 "./HZB_FLEX.c"
  rot_transpose(mcrotaMono_sample_arm, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaenergy_mono = coords_add(mcposaMono_sample_arm, mctc2);
  mctc1 = coords_sub(mcposaMono_sample_arm, mcposaenergy_mono);
  mcposrenergy_mono = rot_apply(mcrotaenergy_mono, mctc1);
  mcDEBUG_COMPONENT("energy_mono", mcposaenergy_mono, mcrotaenergy_mono)
  mccomp_posa[27] = mcposaenergy_mono;
  mccomp_posr[27] = mcposrenergy_mono;
  mcNCounter[27]  = mcPCounter[27] = mcP2Counter[27] = 0;
  mcAbsorbProp[27]= 0;
    /* Component energy_pre_sample. */
  /* Setting parameters for component energy_pre_sample. */
  SIG_MESSAGE("energy_pre_sample (Init:SetPar)");
#line 387 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  if("energy-pre-sample.dat") strncpy(mccenergy_pre_sample_filename, "energy-pre-sample.dat" ? "energy-pre-sample.dat" : "", 16384); else mccenergy_pre_sample_filename[0]='\0';
#line 53 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccenergy_pre_sample_xmin = -0.05;
#line 53 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccenergy_pre_sample_xmax = 0.05;
#line 53 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccenergy_pre_sample_ymin = -0.05;
#line 53 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccenergy_pre_sample_ymax = 0.05;
#line 388 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccenergy_pre_sample_xwidth = 0.08;
#line 388 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccenergy_pre_sample_yheight = 0.08;
#line 388 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccenergy_pre_sample_Emin = e_min;
#line 388 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccenergy_pre_sample_Emax = e_max;
#line 387 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccenergy_pre_sample_restore_neutron = 1;
#line 54 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccenergy_pre_sample_nowritefile = 0;
#line 10215 "./HZB_FLEX.c"

  SIG_MESSAGE("energy_pre_sample (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 10222 "./HZB_FLEX.c"
  rot_mul(mctr1, mcrotaMono_sample_arm, mcrotaenergy_pre_sample);
  rot_transpose(mcrotaenergy_mono, mctr1);
  rot_mul(mcrotaenergy_pre_sample, mctr1, mcrotrenergy_pre_sample);
  mctc1 = coords_set(
#line 389 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    0,
#line 389 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    0,
#line 389 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    L2 -0.1);
#line 10233 "./HZB_FLEX.c"
  rot_transpose(mcrotaMono_sample_arm, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaenergy_pre_sample = coords_add(mcposaMono_sample_arm, mctc2);
  mctc1 = coords_sub(mcposaenergy_mono, mcposaenergy_pre_sample);
  mcposrenergy_pre_sample = rot_apply(mcrotaenergy_pre_sample, mctc1);
  mcDEBUG_COMPONENT("energy_pre_sample", mcposaenergy_pre_sample, mcrotaenergy_pre_sample)
  mccomp_posa[28] = mcposaenergy_pre_sample;
  mccomp_posr[28] = mcposrenergy_pre_sample;
  mcNCounter[28]  = mcPCounter[28] = mcP2Counter[28] = 0;
  mcAbsorbProp[28]= 0;
    /* Component Sample_center. */
  /* Setting parameters for component Sample_center. */
  SIG_MESSAGE("Sample_center (Init:SetPar)");

  SIG_MESSAGE("Sample_center (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 400 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    (0)*DEG2RAD,
#line 400 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    (mcipA3)*DEG2RAD,
#line 400 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    (0)*DEG2RAD);
#line 10256 "./HZB_FLEX.c"
  rot_mul(mctr1, mcrotaMono_sample_arm, mcrotaSample_center);
  rot_transpose(mcrotaenergy_pre_sample, mctr1);
  rot_mul(mcrotaSample_center, mctr1, mcrotrSample_center);
  mctc1 = coords_set(
#line 400 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    0,
#line 400 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    0,
#line 400 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    L2);
#line 10267 "./HZB_FLEX.c"
  rot_transpose(mcrotaMono_sample_arm, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaSample_center = coords_add(mcposaMono_sample_arm, mctc2);
  mctc1 = coords_sub(mcposaenergy_pre_sample, mcposaSample_center);
  mcposrSample_center = rot_apply(mcrotaSample_center, mctc1);
  mcDEBUG_COMPONENT("Sample_center", mcposaSample_center, mcrotaSample_center)
  mccomp_posa[29] = mcposaSample_center;
  mccomp_posr[29] = mcposrSample_center;
  mcNCounter[29]  = mcPCounter[29] = mcP2Counter[29] = 0;
  mcAbsorbProp[29]= 0;
    /* Component Sample_analyser_arm. */
  /* Setting parameters for component Sample_analyser_arm. */
  SIG_MESSAGE("Sample_analyser_arm (Init:SetPar)");

  SIG_MESSAGE("Sample_analyser_arm (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 405 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    (0)*DEG2RAD,
#line 405 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    (mcipA4)*DEG2RAD,
#line 405 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    (0)*DEG2RAD);
#line 10290 "./HZB_FLEX.c"
  rot_mul(mctr1, mcrotaMono_sample_arm, mcrotaSample_analyser_arm);
  rot_transpose(mcrotaSample_center, mctr1);
  rot_mul(mcrotaSample_analyser_arm, mctr1, mcrotrSample_analyser_arm);
  mctc1 = coords_set(
#line 405 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    0,
#line 405 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    0,
#line 405 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    0);
#line 10301 "./HZB_FLEX.c"
  rot_transpose(mcrotaSample_center, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaSample_analyser_arm = coords_add(mcposaSample_center, mctc2);
  mctc1 = coords_sub(mcposaSample_center, mcposaSample_analyser_arm);
  mcposrSample_analyser_arm = rot_apply(mcrotaSample_analyser_arm, mctc1);
  mcDEBUG_COMPONENT("Sample_analyser_arm", mcposaSample_analyser_arm, mcrotaSample_analyser_arm)
  mccomp_posa[30] = mcposaSample_analyser_arm;
  mccomp_posr[30] = mcposrSample_analyser_arm;
  mcNCounter[30]  = mcPCounter[30] = mcP2Counter[30] = 0;
  mcAbsorbProp[30]= 0;
    /* Component div_mono. */
  /* Setting parameters for component div_mono. */
  SIG_MESSAGE("div_mono (Init:SetPar)");
#line 408 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  if("div-sample.dat") strncpy(mccdiv_mono_filename, "div-sample.dat" ? "div-sample.dat" : "", 16384); else mccdiv_mono_filename[0]='\0';
#line 55 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccdiv_mono_xmin = -0.05;
#line 55 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccdiv_mono_xmax = 0.05;
#line 55 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccdiv_mono_ymin = -0.05;
#line 55 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccdiv_mono_ymax = 0.05;
#line 409 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccdiv_mono_xwidth = 0.1;
#line 409 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccdiv_mono_yheight = 0.1;
#line 410 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccdiv_mono_maxdiv_h = 4;
#line 410 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccdiv_mono_maxdiv_v = 4;
#line 409 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccdiv_mono_restore_neutron = 1;
#line 56 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccdiv_mono_nx = 0;
#line 56 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccdiv_mono_ny = 0;
#line 56 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccdiv_mono_nz = 1;
#line 56 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccdiv_mono_nowritefile = 0;
#line 10343 "./HZB_FLEX.c"

  SIG_MESSAGE("div_mono (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 10350 "./HZB_FLEX.c"
  rot_mul(mctr1, mcrotaMono_sample_arm, mcrotadiv_mono);
  rot_transpose(mcrotaSample_analyser_arm, mctr1);
  rot_mul(mcrotadiv_mono, mctr1, mcrotrdiv_mono);
  mctc1 = coords_set(
#line 411 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    0,
#line 411 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    0,
#line 411 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    L2 -0.05);
#line 10361 "./HZB_FLEX.c"
  rot_transpose(mcrotaMono_sample_arm, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposadiv_mono = coords_add(mcposaMono_sample_arm, mctc2);
  mctc1 = coords_sub(mcposaSample_analyser_arm, mcposadiv_mono);
  mcposrdiv_mono = rot_apply(mcrotadiv_mono, mctc1);
  mcDEBUG_COMPONENT("div_mono", mcposadiv_mono, mcrotadiv_mono)
  mccomp_posa[31] = mcposadiv_mono;
  mccomp_posr[31] = mcposrdiv_mono;
  mcNCounter[31]  = mcPCounter[31] = mcP2Counter[31] = 0;
  mcAbsorbProp[31]= 0;
    /* Component div_mono_H. */
  /* Setting parameters for component div_mono_H. */
  SIG_MESSAGE("div_mono_H (Init:SetPar)");
#line 414 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  if("H-div-sample.dat") strncpy(mccdiv_mono_H_filename, "H-div-sample.dat" ? "H-div-sample.dat" : "", 16384); else mccdiv_mono_H_filename[0]='\0';
#line 51 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccdiv_mono_H_xmin = -0.05;
#line 51 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccdiv_mono_H_xmax = 0.05;
#line 51 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccdiv_mono_H_ymin = -0.05;
#line 51 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccdiv_mono_H_ymax = 0.05;
#line 415 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccdiv_mono_H_xwidth = 0.005;
#line 415 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccdiv_mono_H_yheight = 0.7;
#line 416 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccdiv_mono_H_h_maxdiv = 4;
#line 415 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccdiv_mono_H_restore_neutron = 1;
#line 52 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccdiv_mono_H_nowritefile = 0;
#line 10395 "./HZB_FLEX.c"

  SIG_MESSAGE("div_mono_H (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 10402 "./HZB_FLEX.c"
  rot_mul(mctr1, mcrotaMono_sample_arm, mcrotadiv_mono_H);
  rot_transpose(mcrotadiv_mono, mctr1);
  rot_mul(mcrotadiv_mono_H, mctr1, mcrotrdiv_mono_H);
  mctc1 = coords_set(
#line 417 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    0,
#line 417 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    0,
#line 417 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    L2 -0.04);
#line 10413 "./HZB_FLEX.c"
  rot_transpose(mcrotaMono_sample_arm, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposadiv_mono_H = coords_add(mcposaMono_sample_arm, mctc2);
  mctc1 = coords_sub(mcposadiv_mono, mcposadiv_mono_H);
  mcposrdiv_mono_H = rot_apply(mcrotadiv_mono_H, mctc1);
  mcDEBUG_COMPONENT("div_mono_H", mcposadiv_mono_H, mcrotadiv_mono_H)
  mccomp_posa[32] = mcposadiv_mono_H;
  mccomp_posr[32] = mcposrdiv_mono_H;
  mcNCounter[32]  = mcPCounter[32] = mcP2Counter[32] = 0;
  mcAbsorbProp[32]= 0;
    /* Component psd_sam. */
  /* Setting parameters for component psd_sam. */
  SIG_MESSAGE("psd_sam (Init:SetPar)");
#line 420 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccpsd_sam_nx = 60;
#line 420 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccpsd_sam_ny = 60;
#line 420 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  if("psd-sam.dat") strncpy(mccpsd_sam_filename, "psd-sam.dat" ? "psd-sam.dat" : "", 16384); else mccpsd_sam_filename[0]='\0';
#line 50 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccpsd_sam_xmin = -0.05;
#line 50 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccpsd_sam_xmax = 0.05;
#line 50 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccpsd_sam_ymin = -0.05;
#line 50 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccpsd_sam_ymax = 0.05;
#line 421 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccpsd_sam_xwidth = 0.20;
#line 421 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccpsd_sam_yheight = 0.20;
#line 421 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mccpsd_sam_restore_neutron = 1;
#line 10447 "./HZB_FLEX.c"

  SIG_MESSAGE("psd_sam (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 422 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    (0)*DEG2RAD,
#line 422 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    (0)*DEG2RAD,
#line 422 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    (0)*DEG2RAD);
#line 10457 "./HZB_FLEX.c"
  rot_mul(mctr1, mcrotaMono_sample_arm, mcrotapsd_sam);
  rot_transpose(mcrotadiv_mono_H, mctr1);
  rot_mul(mcrotapsd_sam, mctr1, mcrotrpsd_sam);
  mctc1 = coords_set(
#line 422 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    0,
#line 422 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    0,
#line 422 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    0);
#line 10468 "./HZB_FLEX.c"
  rot_transpose(mcrotaSample_center, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposapsd_sam = coords_add(mcposaSample_center, mctc2);
  mctc1 = coords_sub(mcposadiv_mono_H, mcposapsd_sam);
  mcposrpsd_sam = rot_apply(mcrotapsd_sam, mctc1);
  mcDEBUG_COMPONENT("psd_sam", mcposapsd_sam, mcrotapsd_sam)
  mccomp_posa[33] = mcposapsd_sam;
  mccomp_posr[33] = mcposrpsd_sam;
  mcNCounter[33]  = mcPCounter[33] = mcP2Counter[33] = 0;
  mcAbsorbProp[33]= 0;
    /* Component 1dpsd. */
  /* Setting parameters for component 1dpsd. */
  SIG_MESSAGE("1dpsd (Init:SetPar)");
#line 424 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  if("linpsdsam.dat") strncpy(mcc1dpsd_filename, "linpsdsam.dat" ? "linpsdsam.dat" : "", 16384); else mcc1dpsd_filename[0]='\0';
#line 425 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mcc1dpsd_xmin = -0.1;
#line 425 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mcc1dpsd_xmax = 0.1;
#line 425 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mcc1dpsd_ymin = -0.1;
#line 425 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mcc1dpsd_ymax = 0.1;
#line 47 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mcc1dpsd_xwidth = 0;
#line 47 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mcc1dpsd_yheight = 0;
#line 47 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mcc1dpsd_restore_neutron = 0;
#line 47 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
  mcc1dpsd_nowritefile = 0;
#line 10500 "./HZB_FLEX.c"

  SIG_MESSAGE("1dpsd (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 426 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    (0)*DEG2RAD,
#line 426 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    (0)*DEG2RAD,
#line 426 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    (0)*DEG2RAD);
#line 10510 "./HZB_FLEX.c"
  rot_mul(mctr1, mcrotaMono_sample_arm, mcrota1dpsd);
  rot_transpose(mcrotapsd_sam, mctr1);
  rot_mul(mcrota1dpsd, mctr1, mcrotr1dpsd);
  mctc1 = coords_set(
#line 426 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    0,
#line 426 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    0,
#line 426 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/HZB_FLEX/HZB_FLEX.instr"
    0);
#line 10521 "./HZB_FLEX.c"
  rot_transpose(mcrotaSample_center, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposa1dpsd = coords_add(mcposaSample_center, mctc2);
  mctc1 = coords_sub(mcposapsd_sam, mcposa1dpsd);
  mcposr1dpsd = rot_apply(mcrota1dpsd, mctc1);
  mcDEBUG_COMPONENT("1dpsd", mcposa1dpsd, mcrota1dpsd)
  mccomp_posa[34] = mcposa1dpsd;
  mccomp_posr[34] = mcposr1dpsd;
  mcNCounter[34]  = mcPCounter[34] = mcP2Counter[34] = 0;
  mcAbsorbProp[34]= 0;
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
#line 10558 "./HZB_FLEX.c"
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

  /* Initializations for component Gen_Source. */
  SIG_MESSAGE("Gen_Source (Init)");
#define mccompcurname  Gen_Source
#define mccompcurtype  Source_gen
#define mccompcurindex 2
#define p_in mccGen_Source_p_in
#define lambda1 mccGen_Source_lambda1
#define lambda2 mccGen_Source_lambda2
#define lambda3 mccGen_Source_lambda3
#define pTable mccGen_Source_pTable
#define pTable_x mccGen_Source_pTable_x
#define pTable_y mccGen_Source_pTable_y
#define pTable_xmin mccGen_Source_pTable_xmin
#define pTable_xmax mccGen_Source_pTable_xmax
#define pTable_xsum mccGen_Source_pTable_xsum
#define pTable_ymin mccGen_Source_pTable_ymin
#define pTable_ymax mccGen_Source_pTable_ymax
#define pTable_ysum mccGen_Source_pTable_ysum
#define pTable_dxmin mccGen_Source_pTable_dxmin
#define pTable_dxmax mccGen_Source_pTable_dxmax
#define pTable_dymin mccGen_Source_pTable_dymin
#define pTable_dymax mccGen_Source_pTable_dymax
#define flux_file mccGen_Source_flux_file
#define xdiv_file mccGen_Source_xdiv_file
#define ydiv_file mccGen_Source_ydiv_file
#define radius mccGen_Source_radius
#define dist mccGen_Source_dist
#define focus_xw mccGen_Source_focus_xw
#define focus_yh mccGen_Source_focus_yh
#define focus_aw mccGen_Source_focus_aw
#define focus_ah mccGen_Source_focus_ah
#define E0 mccGen_Source_E0
#define dE mccGen_Source_dE
#define lambda0 mccGen_Source_lambda0
#define dlambda mccGen_Source_dlambda
#define I1 mccGen_Source_I1
#define yheight mccGen_Source_yheight
#define xwidth mccGen_Source_xwidth
#define verbose mccGen_Source_verbose
#define T1 mccGen_Source_T1
#define flux_file_perAA mccGen_Source_flux_file_perAA
#define flux_file_log mccGen_Source_flux_file_log
#define Lmin mccGen_Source_Lmin
#define Lmax mccGen_Source_Lmax
#define Emin mccGen_Source_Emin
#define Emax mccGen_Source_Emax
#define T2 mccGen_Source_T2
#define I2 mccGen_Source_I2
#define T3 mccGen_Source_T3
#define I3 mccGen_Source_I3
#define zdepth mccGen_Source_zdepth
#define target_index mccGen_Source_target_index
#line 206 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../sources/Source_gen.comp"
{
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
}
#line 10895 "./HZB_FLEX.c"
#undef target_index
#undef zdepth
#undef I3
#undef T3
#undef I2
#undef T2
#undef Emax
#undef Emin
#undef Lmax
#undef Lmin
#undef flux_file_log
#undef flux_file_perAA
#undef T1
#undef verbose
#undef xwidth
#undef yheight
#undef I1
#undef dlambda
#undef lambda0
#undef dE
#undef E0
#undef focus_ah
#undef focus_aw
#undef focus_yh
#undef focus_xw
#undef dist
#undef radius
#undef ydiv_file
#undef xdiv_file
#undef flux_file
#undef pTable_dymax
#undef pTable_dymin
#undef pTable_dxmax
#undef pTable_dxmin
#undef pTable_ysum
#undef pTable_ymax
#undef pTable_ymin
#undef pTable_xsum
#undef pTable_xmax
#undef pTable_xmin
#undef pTable_y
#undef pTable_x
#undef pTable
#undef lambda3
#undef lambda2
#undef lambda1
#undef p_in
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component NL1A_NL1B_NL2_NL3_InPile_Entrance_Window. */
  SIG_MESSAGE("NL1A_NL1B_NL2_NL3_InPile_Entrance_Window (Init)");

  /* Initializations for component NL1A_NL1B_NL2_NL3_InPile. */
  SIG_MESSAGE("NL1A_NL1B_NL2_NL3_InPile (Init)");
#define mccompcurname  NL1A_NL1B_NL2_NL3_InPile
#define mccompcurtype  Guide
#define mccompcurindex 4
#define pTable mccNL1A_NL1B_NL2_NL3_InPile_pTable
#define reflect mccNL1A_NL1B_NL2_NL3_InPile_reflect
#define w1 mccNL1A_NL1B_NL2_NL3_InPile_w1
#define h1 mccNL1A_NL1B_NL2_NL3_InPile_h1
#define w2 mccNL1A_NL1B_NL2_NL3_InPile_w2
#define h2 mccNL1A_NL1B_NL2_NL3_InPile_h2
#define l mccNL1A_NL1B_NL2_NL3_InPile_l
#define R0 mccNL1A_NL1B_NL2_NL3_InPile_R0
#define Qc mccNL1A_NL1B_NL2_NL3_InPile_Qc
#define alpha mccNL1A_NL1B_NL2_NL3_InPile_alpha
#define m mccNL1A_NL1B_NL2_NL3_InPile_m
#define W mccNL1A_NL1B_NL2_NL3_InPile_W
#line 74 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide.comp"
{
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
}
#line 10986 "./HZB_FLEX.c"
#undef W
#undef m
#undef alpha
#undef Qc
#undef R0
#undef l
#undef h2
#undef w2
#undef h1
#undef w1
#undef reflect
#undef pTable
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component NL1B_Straight1_Entrance_Window. */
  SIG_MESSAGE("NL1B_Straight1_Entrance_Window (Init)");

  /* Initializations for component NL1B_Straight1. */
  SIG_MESSAGE("NL1B_Straight1 (Init)");
#define mccompcurname  NL1B_Straight1
#define mccompcurtype  Guide
#define mccompcurindex 6
#define pTable mccNL1B_Straight1_pTable
#define reflect mccNL1B_Straight1_reflect
#define w1 mccNL1B_Straight1_w1
#define h1 mccNL1B_Straight1_h1
#define w2 mccNL1B_Straight1_w2
#define h2 mccNL1B_Straight1_h2
#define l mccNL1B_Straight1_l
#define R0 mccNL1B_Straight1_R0
#define Qc mccNL1B_Straight1_Qc
#define alpha mccNL1B_Straight1_alpha
#define m mccNL1B_Straight1_m
#define W mccNL1B_Straight1_W
#line 74 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide.comp"
{
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
}
#line 11042 "./HZB_FLEX.c"
#undef W
#undef m
#undef alpha
#undef Qc
#undef R0
#undef l
#undef h2
#undef w2
#undef h1
#undef w1
#undef reflect
#undef pTable
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component NL1B_Curved_Entrance_Window. */
  SIG_MESSAGE("NL1B_Curved_Entrance_Window (Init)");

  /* Initializations for component NL1B_Curved_Guide_1. */
  SIG_MESSAGE("NL1B_Curved_Guide_1 (Init)");
#define mccompcurname  NL1B_Curved_Guide_1
#define mccompcurtype  Guide_curved
#define mccompcurindex 8
#define w1 mccNL1B_Curved_Guide_1_w1
#define h1 mccNL1B_Curved_Guide_1_h1
#define l mccNL1B_Curved_Guide_1_l
#define R0 mccNL1B_Curved_Guide_1_R0
#define Qc mccNL1B_Curved_Guide_1_Qc
#define alpha mccNL1B_Curved_Guide_1_alpha
#define m mccNL1B_Curved_Guide_1_m
#define W mccNL1B_Curved_Guide_1_W
#define curvature mccNL1B_Curved_Guide_1_curvature
#line 60 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../contrib/Guide_curved.comp"
{
if (mcgravitation) fprintf(stderr,"WARNING: Guide_curved: %s: "
    "This component produces wrong results with gravitation !\n"
    "Use Guide_gravity.\n",
    NAME_CURRENT_COMP);
}
#line 11083 "./HZB_FLEX.c"
#undef curvature
#undef W
#undef m
#undef alpha
#undef Qc
#undef R0
#undef l
#undef h1
#undef w1
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component NL1B_Velocity_Selector_Gap_Entrance_Window. */
  SIG_MESSAGE("NL1B_Velocity_Selector_Gap_Entrance_Window (Init)");

  /* Initializations for component Before_selec. */
  SIG_MESSAGE("Before_selec (Init)");
#define mccompcurname  Before_selec
#define mccompcurtype  PSD_monitor
#define mccompcurindex 10
#define PSD_N mccBefore_selec_PSD_N
#define PSD_p mccBefore_selec_PSD_p
#define PSD_p2 mccBefore_selec_PSD_p2
#define nx mccBefore_selec_nx
#define ny mccBefore_selec_ny
#define filename mccBefore_selec_filename
#define xmin mccBefore_selec_xmin
#define xmax mccBefore_selec_xmax
#define ymin mccBefore_selec_ymin
#define ymax mccBefore_selec_ymax
#define xwidth mccBefore_selec_xwidth
#define yheight mccBefore_selec_yheight
#define restore_neutron mccBefore_selec_restore_neutron
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
#line 11143 "./HZB_FLEX.c"
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

  /* Initializations for component SELEC. */
  SIG_MESSAGE("SELEC (Init)");
#define mccompcurname  SELEC
#define mccompcurtype  Selector
#define mccompcurindex 11
#define xmin mccSELEC_xmin
#define xmax mccSELEC_xmax
#define ymin mccSELEC_ymin
#define ymax mccSELEC_ymax
#define length mccSELEC_length
#define xwidth mccSELEC_xwidth
#define yheight mccSELEC_yheight
#define nslit mccSELEC_nslit
#define d mccSELEC_d
#define radius mccSELEC_radius
#define alpha mccSELEC_alpha
#define nu mccSELEC_nu
#line 68 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Selector.comp"
{
if (xwidth > 0)  { xmax=xwidth/2;  xmin=-xmax; }
  if (yheight > 0) { ymax=yheight/2; ymin=-ymax; }
}
#line 11183 "./HZB_FLEX.c"
#undef nu
#undef alpha
#undef radius
#undef d
#undef nslit
#undef yheight
#undef xwidth
#undef length
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component After_selec. */
  SIG_MESSAGE("After_selec (Init)");
#define mccompcurname  After_selec
#define mccompcurtype  PSD_monitor
#define mccompcurindex 12
#define PSD_N mccAfter_selec_PSD_N
#define PSD_p mccAfter_selec_PSD_p
#define PSD_p2 mccAfter_selec_PSD_p2
#define nx mccAfter_selec_nx
#define ny mccAfter_selec_ny
#define filename mccAfter_selec_filename
#define xmin mccAfter_selec_xmin
#define xmax mccAfter_selec_xmax
#define ymin mccAfter_selec_ymin
#define ymax mccAfter_selec_ymax
#define xwidth mccAfter_selec_xwidth
#define yheight mccAfter_selec_yheight
#define restore_neutron mccAfter_selec_restore_neutron
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
#line 11243 "./HZB_FLEX.c"
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

  /* Initializations for component NL1B_Curved_2_Entrance_Window. */
  SIG_MESSAGE("NL1B_Curved_2_Entrance_Window (Init)");

  /* Initializations for component NL1B_Curved_Guide_2. */
  SIG_MESSAGE("NL1B_Curved_Guide_2 (Init)");
#define mccompcurname  NL1B_Curved_Guide_2
#define mccompcurtype  Guide_curved
#define mccompcurindex 14
#define w1 mccNL1B_Curved_Guide_2_w1
#define h1 mccNL1B_Curved_Guide_2_h1
#define l mccNL1B_Curved_Guide_2_l
#define R0 mccNL1B_Curved_Guide_2_R0
#define Qc mccNL1B_Curved_Guide_2_Qc
#define alpha mccNL1B_Curved_Guide_2_alpha
#define m mccNL1B_Curved_Guide_2_m
#define W mccNL1B_Curved_Guide_2_W
#define curvature mccNL1B_Curved_Guide_2_curvature
#line 60 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../contrib/Guide_curved.comp"
{
if (mcgravitation) fprintf(stderr,"WARNING: Guide_curved: %s: "
    "This component produces wrong results with gravitation !\n"
    "Use Guide_gravity.\n",
    NAME_CURRENT_COMP);
}
#line 11285 "./HZB_FLEX.c"
#undef curvature
#undef W
#undef m
#undef alpha
#undef Qc
#undef R0
#undef l
#undef h1
#undef w1
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component NL1B_Straight2_Entrance_Window. */
  SIG_MESSAGE("NL1B_Straight2_Entrance_Window (Init)");

  /* Initializations for component NL1B_Straight2. */
  SIG_MESSAGE("NL1B_Straight2 (Init)");
#define mccompcurname  NL1B_Straight2
#define mccompcurtype  Guide
#define mccompcurindex 16
#define pTable mccNL1B_Straight2_pTable
#define reflect mccNL1B_Straight2_reflect
#define w1 mccNL1B_Straight2_w1
#define h1 mccNL1B_Straight2_h1
#define w2 mccNL1B_Straight2_w2
#define h2 mccNL1B_Straight2_h2
#define l mccNL1B_Straight2_l
#define R0 mccNL1B_Straight2_R0
#define Qc mccNL1B_Straight2_Qc
#define alpha mccNL1B_Straight2_alpha
#define m mccNL1B_Straight2_m
#define W mccNL1B_Straight2_W
#line 74 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide.comp"
{
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
}
#line 11338 "./HZB_FLEX.c"
#undef W
#undef m
#undef alpha
#undef Qc
#undef R0
#undef l
#undef h2
#undef w2
#undef h1
#undef w1
#undef reflect
#undef pTable
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component NL1B_Elliptical_Entrance_Window. */
  SIG_MESSAGE("NL1B_Elliptical_Entrance_Window (Init)");

  /* Initializations for component elliptical_piece. */
  SIG_MESSAGE("elliptical_piece (Init)");
#define mccompcurname  elliptical_piece
#define mccompcurtype  Guide_tapering
#define mccompcurindex 18
#define w1c mccelliptical_piece_w1c
#define w2c mccelliptical_piece_w2c
#define ww mccelliptical_piece_ww
#define hh mccelliptical_piece_hh
#define whalf mccelliptical_piece_whalf
#define hhalf mccelliptical_piece_hhalf
#define lwhalf mccelliptical_piece_lwhalf
#define lhhalf mccelliptical_piece_lhhalf
#define h1_in mccelliptical_piece_h1_in
#define h2_out mccelliptical_piece_h2_out
#define w1_in mccelliptical_piece_w1_in
#define w2_out mccelliptical_piece_w2_out
#define l_seg mccelliptical_piece_l_seg
#define seg mccelliptical_piece_seg
#define h12 mccelliptical_piece_h12
#define h2 mccelliptical_piece_h2
#define w12 mccelliptical_piece_w12
#define w2 mccelliptical_piece_w2
#define a_ell_q mccelliptical_piece_a_ell_q
#define b_ell_q mccelliptical_piece_b_ell_q
#define lbw mccelliptical_piece_lbw
#define lbh mccelliptical_piece_lbh
#define mxi mccelliptical_piece_mxi
#define u1 mccelliptical_piece_u1
#define u2 mccelliptical_piece_u2
#define div1 mccelliptical_piece_div1
#define p2_para mccelliptical_piece_p2_para
#define test mccelliptical_piece_test
#define Div1 mccelliptical_piece_Div1
#define i mccelliptical_piece_i
#define ii mccelliptical_piece_ii
#define seg mccelliptical_piece_seg
#define fu mccelliptical_piece_fu
#define pos mccelliptical_piece_pos
#define file_name mccelliptical_piece_file_name
#define ep mccelliptical_piece_ep
#define num mccelliptical_piece_num
#define rotation_h mccelliptical_piece_rotation_h
#define rotation_v mccelliptical_piece_rotation_v
#define option mccelliptical_piece_option
#define w1 mccelliptical_piece_w1
#define h1 mccelliptical_piece_h1
#define l mccelliptical_piece_l
#define linw mccelliptical_piece_linw
#define loutw mccelliptical_piece_loutw
#define linh mccelliptical_piece_linh
#define louth mccelliptical_piece_louth
#define R0 mccelliptical_piece_R0
#define Qcx mccelliptical_piece_Qcx
#define Qcy mccelliptical_piece_Qcy
#define alphax mccelliptical_piece_alphax
#define alphay mccelliptical_piece_alphay
#define W mccelliptical_piece_W
#define mx mccelliptical_piece_mx
#define my mccelliptical_piece_my
#define segno mccelliptical_piece_segno
#define curvature mccelliptical_piece_curvature
#define curvature_v mccelliptical_piece_curvature_v
#line 116 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide_tapering.comp"
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
#line 11757 "./HZB_FLEX.c"
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

  /* Initializations for component Virtual_source. */
  SIG_MESSAGE("Virtual_source (Init)");
#define mccompcurname  Virtual_source
#define mccompcurtype  Slit
#define mccompcurindex 19
#define xmin mccVirtual_source_xmin
#define xmax mccVirtual_source_xmax
#define ymin mccVirtual_source_ymin
#define ymax mccVirtual_source_ymax
#define radius mccVirtual_source_radius
#define xwidth mccVirtual_source_xwidth
#define yheight mccVirtual_source_yheight
#line 50 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Slit.comp"
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
#line 11852 "./HZB_FLEX.c"
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

  /* Initializations for component NL1B_Vertical_Guide_Entrance_Window. */
  SIG_MESSAGE("NL1B_Vertical_Guide_Entrance_Window (Init)");

  /* Initializations for component NL1B_Vertical_Guide. */
  SIG_MESSAGE("NL1B_Vertical_Guide (Init)");
#define mccompcurname  NL1B_Vertical_Guide
#define mccompcurtype  Guide_channeled
#define mccompcurindex 21
#define w1c mccNL1B_Vertical_Guide_w1c
#define w2c mccNL1B_Vertical_Guide_w2c
#define ww mccNL1B_Vertical_Guide_ww
#define hh mccNL1B_Vertical_Guide_hh
#define whalf mccNL1B_Vertical_Guide_whalf
#define hhalf mccNL1B_Vertical_Guide_hhalf
#define lwhalf mccNL1B_Vertical_Guide_lwhalf
#define lhhalf mccNL1B_Vertical_Guide_lhhalf
#define w1 mccNL1B_Vertical_Guide_w1
#define h1 mccNL1B_Vertical_Guide_h1
#define w2 mccNL1B_Vertical_Guide_w2
#define h2 mccNL1B_Vertical_Guide_h2
#define l mccNL1B_Vertical_Guide_l
#define R0 mccNL1B_Vertical_Guide_R0
#define Qc mccNL1B_Vertical_Guide_Qc
#define alpha mccNL1B_Vertical_Guide_alpha
#define m mccNL1B_Vertical_Guide_m
#define nslit mccNL1B_Vertical_Guide_nslit
#define d mccNL1B_Vertical_Guide_d
#define Qcx mccNL1B_Vertical_Guide_Qcx
#define Qcy mccNL1B_Vertical_Guide_Qcy
#define alphax mccNL1B_Vertical_Guide_alphax
#define alphay mccNL1B_Vertical_Guide_alphay
#define W mccNL1B_Vertical_Guide_W
#define mx mccNL1B_Vertical_Guide_mx
#define my mccNL1B_Vertical_Guide_my
#define nu mccNL1B_Vertical_Guide_nu
#define phase mccNL1B_Vertical_Guide_phase
#line 89 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide_channeled.comp"
{
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
}
#line 11941 "./HZB_FLEX.c"
#undef phase
#undef nu
#undef my
#undef mx
#undef W
#undef alphay
#undef alphax
#undef Qcy
#undef Qcx
#undef d
#undef nslit
#undef m
#undef alpha
#undef Qc
#undef R0
#undef l
#undef h2
#undef w2
#undef h1
#undef w1
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

  /* Initializations for component NL1B_Guide_Exit. */
  SIG_MESSAGE("NL1B_Guide_Exit (Init)");

  /* Initializations for component energy_endguide. */
  SIG_MESSAGE("energy_endguide (Init)");
#define mccompcurname  energy_endguide
#define mccompcurtype  E_monitor
#define mccompcurindex 23
#define nE mccenergy_endguide_nE
#define E_N mccenergy_endguide_E_N
#define E_p mccenergy_endguide_E_p
#define E_p2 mccenergy_endguide_E_p2
#define S_p mccenergy_endguide_S_p
#define S_pE mccenergy_endguide_S_pE
#define S_pE2 mccenergy_endguide_S_pE2
#define filename mccenergy_endguide_filename
#define xmin mccenergy_endguide_xmin
#define xmax mccenergy_endguide_xmax
#define ymin mccenergy_endguide_ymin
#define ymax mccenergy_endguide_ymax
#define xwidth mccenergy_endguide_xwidth
#define yheight mccenergy_endguide_yheight
#define Emin mccenergy_endguide_Emin
#define Emax mccenergy_endguide_Emax
#define restore_neutron mccenergy_endguide_restore_neutron
#define nowritefile mccenergy_endguide_nowritefile
#line 66 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/E_monitor.comp"
{
int i;

if (xwidth  > 0) { xmax = xwidth/2;  xmin = -xmax; }
    if (yheight > 0) { ymax = yheight/2; ymin = -ymax; }

    if ((xmin >= xmax) || (ymin >= ymax)) {
            printf("E_monitor: %s: Null detection area !\n"
                   "ERROR      (xwidth,yheight,xmin,xmax,ymin,ymax). Exiting",
           NAME_CURRENT_COMP);
      exit(0);
    }

    for (i=0; i<nE; i++)
    {
      E_N[i] = 0;
      E_p[i] = 0;
      E_p2[i] = 0;
    }
    S_p = S_pE = S_pE2 = 0;
}
#line 12022 "./HZB_FLEX.c"
#undef nowritefile
#undef restore_neutron
#undef Emax
#undef Emin
#undef yheight
#undef xwidth
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef filename
#undef S_pE2
#undef S_pE
#undef S_p
#undef E_p2
#undef E_p
#undef E_N
#undef nE
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component Mono_center. */
  SIG_MESSAGE("Mono_center (Init)");

  /* Initializations for component Monochromator. */
  SIG_MESSAGE("Monochromator (Init)");
#define mccompcurname  Monochromator
#define mccompcurtype  Monochromator_curved
#define mccompcurindex 25
#define mos_rms_y mccMonochromator_mos_rms_y
#define mos_rms_z mccMonochromator_mos_rms_z
#define mos_rms_max mccMonochromator_mos_rms_max
#define mono_Q mccMonochromator_mono_Q
#define SlabWidth mccMonochromator_SlabWidth
#define SlabHeight mccMonochromator_SlabHeight
#define rTable mccMonochromator_rTable
#define tTable mccMonochromator_tTable
#define row mccMonochromator_row
#define col mccMonochromator_col
#define tiltH mccMonochromator_tiltH
#define tiltV mccMonochromator_tiltV
#define reflect mccMonochromator_reflect
#define transmit mccMonochromator_transmit
#define zwidth mccMonochromator_zwidth
#define yheight mccMonochromator_yheight
#define gap mccMonochromator_gap
#define NH mccMonochromator_NH
#define NV mccMonochromator_NV
#define mosaich mccMonochromator_mosaich
#define mosaicv mccMonochromator_mosaicv
#define r0 mccMonochromator_r0
#define t0 mccMonochromator_t0
#define Q mccMonochromator_Q
#define RV mccMonochromator_RV
#define RH mccMonochromator_RH
#define DM mccMonochromator_DM
#define mosaic mccMonochromator_mosaic
#define width mccMonochromator_width
#define height mccMonochromator_height
#define verbose mccMonochromator_verbose
#define order mccMonochromator_order
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
#line 12166 "./HZB_FLEX.c"
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

  /* Initializations for component Mono_sample_arm. */
  SIG_MESSAGE("Mono_sample_arm (Init)");

  /* Initializations for component energy_mono. */
  SIG_MESSAGE("energy_mono (Init)");
#define mccompcurname  energy_mono
#define mccompcurtype  E_monitor
#define mccompcurindex 27
#define nE mccenergy_mono_nE
#define E_N mccenergy_mono_E_N
#define E_p mccenergy_mono_E_p
#define E_p2 mccenergy_mono_E_p2
#define S_p mccenergy_mono_S_p
#define S_pE mccenergy_mono_S_pE
#define S_pE2 mccenergy_mono_S_pE2
#define filename mccenergy_mono_filename
#define xmin mccenergy_mono_xmin
#define xmax mccenergy_mono_xmax
#define ymin mccenergy_mono_ymin
#define ymax mccenergy_mono_ymax
#define xwidth mccenergy_mono_xwidth
#define yheight mccenergy_mono_yheight
#define Emin mccenergy_mono_Emin
#define Emax mccenergy_mono_Emax
#define restore_neutron mccenergy_mono_restore_neutron
#define nowritefile mccenergy_mono_nowritefile
#line 66 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/E_monitor.comp"
{
int i;

if (xwidth  > 0) { xmax = xwidth/2;  xmin = -xmax; }
    if (yheight > 0) { ymax = yheight/2; ymin = -ymax; }

    if ((xmin >= xmax) || (ymin >= ymax)) {
            printf("E_monitor: %s: Null detection area !\n"
                   "ERROR      (xwidth,yheight,xmin,xmax,ymin,ymax). Exiting",
           NAME_CURRENT_COMP);
      exit(0);
    }

    for (i=0; i<nE; i++)
    {
      E_N[i] = 0;
      E_p[i] = 0;
      E_p2[i] = 0;
    }
    S_p = S_pE = S_pE2 = 0;
}
#line 12251 "./HZB_FLEX.c"
#undef nowritefile
#undef restore_neutron
#undef Emax
#undef Emin
#undef yheight
#undef xwidth
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef filename
#undef S_pE2
#undef S_pE
#undef S_p
#undef E_p2
#undef E_p
#undef E_N
#undef nE
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component energy_pre_sample. */
  SIG_MESSAGE("energy_pre_sample (Init)");
#define mccompcurname  energy_pre_sample
#define mccompcurtype  E_monitor
#define mccompcurindex 28
#define nE mccenergy_pre_sample_nE
#define E_N mccenergy_pre_sample_E_N
#define E_p mccenergy_pre_sample_E_p
#define E_p2 mccenergy_pre_sample_E_p2
#define S_p mccenergy_pre_sample_S_p
#define S_pE mccenergy_pre_sample_S_pE
#define S_pE2 mccenergy_pre_sample_S_pE2
#define filename mccenergy_pre_sample_filename
#define xmin mccenergy_pre_sample_xmin
#define xmax mccenergy_pre_sample_xmax
#define ymin mccenergy_pre_sample_ymin
#define ymax mccenergy_pre_sample_ymax
#define xwidth mccenergy_pre_sample_xwidth
#define yheight mccenergy_pre_sample_yheight
#define Emin mccenergy_pre_sample_Emin
#define Emax mccenergy_pre_sample_Emax
#define restore_neutron mccenergy_pre_sample_restore_neutron
#define nowritefile mccenergy_pre_sample_nowritefile
#line 66 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/E_monitor.comp"
{
int i;

if (xwidth  > 0) { xmax = xwidth/2;  xmin = -xmax; }
    if (yheight > 0) { ymax = yheight/2; ymin = -ymax; }

    if ((xmin >= xmax) || (ymin >= ymax)) {
            printf("E_monitor: %s: Null detection area !\n"
                   "ERROR      (xwidth,yheight,xmin,xmax,ymin,ymax). Exiting",
           NAME_CURRENT_COMP);
      exit(0);
    }

    for (i=0; i<nE; i++)
    {
      E_N[i] = 0;
      E_p[i] = 0;
      E_p2[i] = 0;
    }
    S_p = S_pE = S_pE2 = 0;
}
#line 12319 "./HZB_FLEX.c"
#undef nowritefile
#undef restore_neutron
#undef Emax
#undef Emin
#undef yheight
#undef xwidth
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef filename
#undef S_pE2
#undef S_pE
#undef S_p
#undef E_p2
#undef E_p
#undef E_N
#undef nE
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component Sample_center. */
  SIG_MESSAGE("Sample_center (Init)");

  /* Initializations for component Sample_analyser_arm. */
  SIG_MESSAGE("Sample_analyser_arm (Init)");

  /* Initializations for component div_mono. */
  SIG_MESSAGE("div_mono (Init)");
#define mccompcurname  div_mono
#define mccompcurtype  Divergence_monitor
#define mccompcurindex 31
#define nh mccdiv_mono_nh
#define nv mccdiv_mono_nv
#define Div_N mccdiv_mono_Div_N
#define Div_p mccdiv_mono_Div_p
#define Div_p2 mccdiv_mono_Div_p2
#define filename mccdiv_mono_filename
#define xmin mccdiv_mono_xmin
#define xmax mccdiv_mono_xmax
#define ymin mccdiv_mono_ymin
#define ymax mccdiv_mono_ymax
#define xwidth mccdiv_mono_xwidth
#define yheight mccdiv_mono_yheight
#define maxdiv_h mccdiv_mono_maxdiv_h
#define maxdiv_v mccdiv_mono_maxdiv_v
#define restore_neutron mccdiv_mono_restore_neutron
#define nx mccdiv_mono_nx
#define ny mccdiv_mono_ny
#define nz mccdiv_mono_nz
#define nowritefile mccdiv_mono_nowritefile
#line 66 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/Divergence_monitor.comp"
{
int i,j;

if (xwidth  > 0) { xmax = xwidth/2;  xmin = -xmax; }
    if (yheight > 0) { ymax = yheight/2; ymin = -ymax; }

    if ((xmin >= xmax) || (ymin >= ymax)) {
            printf("Divergence_monitor: %s: Null detection area !\n"
                   "ERROR               (xwidth,yheight,xmin,xmax,ymin,ymax). Exiting",
           NAME_CURRENT_COMP);
      exit(0);
    }

    for (i=0; i<nh; i++)
     for (j=0; j<nv; j++)
     {
      Div_N[i][j] = 0;
      Div_p[i][j] = 0;
      Div_p2[i][j] = 0;
     }
    NORM(nx,ny,nz);
}
#line 12395 "./HZB_FLEX.c"
#undef nowritefile
#undef nz
#undef ny
#undef nx
#undef restore_neutron
#undef maxdiv_v
#undef maxdiv_h
#undef yheight
#undef xwidth
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef filename
#undef Div_p2
#undef Div_p
#undef Div_N
#undef nv
#undef nh
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component div_mono_H. */
  SIG_MESSAGE("div_mono_H (Init)");
#define mccompcurname  div_mono_H
#define mccompcurtype  Hdiv_monitor
#define mccompcurindex 32
#define nh mccdiv_mono_H_nh
#define Div_N mccdiv_mono_H_Div_N
#define Div_p mccdiv_mono_H_Div_p
#define Div_p2 mccdiv_mono_H_Div_p2
#define filename mccdiv_mono_H_filename
#define xmin mccdiv_mono_H_xmin
#define xmax mccdiv_mono_H_xmax
#define ymin mccdiv_mono_H_ymin
#define ymax mccdiv_mono_H_ymax
#define xwidth mccdiv_mono_H_xwidth
#define yheight mccdiv_mono_H_yheight
#define h_maxdiv mccdiv_mono_H_h_maxdiv
#define restore_neutron mccdiv_mono_H_restore_neutron
#define nowritefile mccdiv_mono_H_nowritefile
#line 62 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/Hdiv_monitor.comp"
{
int i;

if (xwidth  > 0) { xmax = xwidth/2;  xmin = -xmax; }
    if (yheight > 0) { ymax = yheight/2; ymin = -ymax; }

    if ((xmin >= xmax) || (ymin >= ymax)) {
            printf("Hdiv_monitor: %s: Null detection area !\n"
                   "ERROR         (xwidth,yheight,xmin,xmax,ymin,ymax). Exiting",
           NAME_CURRENT_COMP);
      exit(0);
    }

    for (i=0; i<nh; i++)
     {
/*       printf("HDiv_monitor: %d\n",i); */
      Div_N[i] = 0;
      Div_p[i] = 0;
      Div_p2[i] = 0;
     }
/*     printf("%d %d %d\n",i,nh,h_maxdiv); */
}
#line 12461 "./HZB_FLEX.c"
#undef nowritefile
#undef restore_neutron
#undef h_maxdiv
#undef yheight
#undef xwidth
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef filename
#undef Div_p2
#undef Div_p
#undef Div_N
#undef nh
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component psd_sam. */
  SIG_MESSAGE("psd_sam (Init)");
#define mccompcurname  psd_sam
#define mccompcurtype  PSD_monitor
#define mccompcurindex 33
#define PSD_N mccpsd_sam_PSD_N
#define PSD_p mccpsd_sam_PSD_p
#define PSD_p2 mccpsd_sam_PSD_p2
#define nx mccpsd_sam_nx
#define ny mccpsd_sam_ny
#define filename mccpsd_sam_filename
#define xmin mccpsd_sam_xmin
#define xmax mccpsd_sam_xmax
#define ymin mccpsd_sam_ymin
#define ymax mccpsd_sam_ymax
#define xwidth mccpsd_sam_xwidth
#define yheight mccpsd_sam_yheight
#define restore_neutron mccpsd_sam_restore_neutron
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
#line 12523 "./HZB_FLEX.c"
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

  /* Initializations for component 1dpsd. */
  SIG_MESSAGE("1dpsd (Init)");
#define mccompcurname  1dpsd
#define mccompcurtype  PSDlin_monitor
#define mccompcurindex 34
#define nx mcc1dpsd_nx
#define PSDlin_N mcc1dpsd_PSDlin_N
#define PSDlin_p mcc1dpsd_PSDlin_p
#define PSDlin_p2 mcc1dpsd_PSDlin_p2
#define filename mcc1dpsd_filename
#define xmin mcc1dpsd_xmin
#define xmax mcc1dpsd_xmax
#define ymin mcc1dpsd_ymin
#define ymax mcc1dpsd_ymax
#define xwidth mcc1dpsd_xwidth
#define yheight mcc1dpsd_yheight
#define restore_neutron mcc1dpsd_restore_neutron
#define nowritefile mcc1dpsd_nowritefile
#line 59 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/PSDlin_monitor.comp"
{
    int i;

    if (xwidth  > 0) { xmax = xwidth/2;  xmin = -xmax; }
    if (yheight > 0) { ymax = yheight/2; ymin = -ymax; }

    if ((xmin >= xmax) || (ymin >= ymax)) {
            printf("PSDlin_monitor: %s: Null detection area !\n"
                   "ERROR           (xwidth,yheight,xmin,xmax,ymin,ymax). Exiting",
           NAME_CURRENT_COMP);
      exit(0);
    }

    for (i=0; i<nx; i++)
    {
      PSDlin_N[i] = 0;
      PSDlin_p[i] = 0;
      PSDlin_p2[i] = 0;
    }
}
#line 12580 "./HZB_FLEX.c"
#undef nowritefile
#undef restore_neutron
#undef yheight
#undef xwidth
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef filename
#undef PSDlin_p2
#undef PSDlin_p
#undef PSDlin_N
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
#line 12751 "./HZB_FLEX.c"
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

  /* TRACE Component Gen_Source [2] */
  mccoordschange(mcposrGen_Source, mcrotrGen_Source,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component Gen_Source (without coords transformations) */
  mcJumpTrace_Gen_Source:
  SIG_MESSAGE("Gen_Source (Trace)");
  mcDEBUG_COMP("Gen_Source")
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

#define mcabsorbComp mcabsorbCompGen_Source
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
#define mccompcurname  Gen_Source
#define mccompcurtype  Source_gen
#define mccompcurindex 2
#define p_in mccGen_Source_p_in
#define lambda1 mccGen_Source_lambda1
#define lambda2 mccGen_Source_lambda2
#define lambda3 mccGen_Source_lambda3
#define pTable mccGen_Source_pTable
#define pTable_x mccGen_Source_pTable_x
#define pTable_y mccGen_Source_pTable_y
#define pTable_xmin mccGen_Source_pTable_xmin
#define pTable_xmax mccGen_Source_pTable_xmax
#define pTable_xsum mccGen_Source_pTable_xsum
#define pTable_ymin mccGen_Source_pTable_ymin
#define pTable_ymax mccGen_Source_pTable_ymax
#define pTable_ysum mccGen_Source_pTable_ysum
#define pTable_dxmin mccGen_Source_pTable_dxmin
#define pTable_dxmax mccGen_Source_pTable_dxmax
#define pTable_dymin mccGen_Source_pTable_dymin
#define pTable_dymax mccGen_Source_pTable_dymax
{   /* Declarations of Gen_Source=Source_gen() SETTING parameters. */
char* flux_file = mccGen_Source_flux_file;
char* xdiv_file = mccGen_Source_xdiv_file;
char* ydiv_file = mccGen_Source_ydiv_file;
MCNUM radius = mccGen_Source_radius;
MCNUM dist = mccGen_Source_dist;
MCNUM focus_xw = mccGen_Source_focus_xw;
MCNUM focus_yh = mccGen_Source_focus_yh;
MCNUM focus_aw = mccGen_Source_focus_aw;
MCNUM focus_ah = mccGen_Source_focus_ah;
MCNUM E0 = mccGen_Source_E0;
MCNUM dE = mccGen_Source_dE;
MCNUM lambda0 = mccGen_Source_lambda0;
MCNUM dlambda = mccGen_Source_dlambda;
MCNUM I1 = mccGen_Source_I1;
MCNUM yheight = mccGen_Source_yheight;
MCNUM xwidth = mccGen_Source_xwidth;
MCNUM verbose = mccGen_Source_verbose;
MCNUM T1 = mccGen_Source_T1;
MCNUM flux_file_perAA = mccGen_Source_flux_file_perAA;
MCNUM flux_file_log = mccGen_Source_flux_file_log;
MCNUM Lmin = mccGen_Source_Lmin;
MCNUM Lmax = mccGen_Source_Lmax;
MCNUM Emin = mccGen_Source_Emin;
MCNUM Emax = mccGen_Source_Emax;
MCNUM T2 = mccGen_Source_T2;
MCNUM I2 = mccGen_Source_I2;
MCNUM T3 = mccGen_Source_T3;
MCNUM I3 = mccGen_Source_I3;
MCNUM zdepth = mccGen_Source_zdepth;
int target_index = mccGen_Source_target_index;
#line 479 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../sources/Source_gen.comp"
{
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
       double xwidth=Table_Value(pTable, lambda, 1);
       if (flux_file_log) xwidth=exp(xwidth);
       p *= xwidth;
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
}
#line 12999 "./HZB_FLEX.c"
}   /* End of Gen_Source=Source_gen() SETTING parameter declarations. */
#undef pTable_dymax
#undef pTable_dymin
#undef pTable_dxmax
#undef pTable_dxmin
#undef pTable_ysum
#undef pTable_ymax
#undef pTable_ymin
#undef pTable_xsum
#undef pTable_xmax
#undef pTable_xmin
#undef pTable_y
#undef pTable_x
#undef pTable
#undef lambda3
#undef lambda2
#undef lambda1
#undef p_in
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompGen_Source:
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

  /* TRACE Component NL1A_NL1B_NL2_NL3_InPile_Entrance_Window [3] */
  mccoordschange(mcposrNL1A_NL1B_NL2_NL3_InPile_Entrance_Window, mcrotrNL1A_NL1B_NL2_NL3_InPile_Entrance_Window,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component NL1A_NL1B_NL2_NL3_InPile_Entrance_Window (without coords transformations) */
  mcJumpTrace_NL1A_NL1B_NL2_NL3_InPile_Entrance_Window:
  SIG_MESSAGE("NL1A_NL1B_NL2_NL3_InPile_Entrance_Window (Trace)");
  mcDEBUG_COMP("NL1A_NL1B_NL2_NL3_InPile_Entrance_Window")
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

#define mcabsorbComp mcabsorbCompNL1A_NL1B_NL2_NL3_InPile_Entrance_Window
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
#define mccompcurname  NL1A_NL1B_NL2_NL3_InPile_Entrance_Window
#define mccompcurtype  Arm
#define mccompcurindex 3
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompNL1A_NL1B_NL2_NL3_InPile_Entrance_Window:
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

  /* TRACE Component NL1A_NL1B_NL2_NL3_InPile [4] */
  mccoordschange(mcposrNL1A_NL1B_NL2_NL3_InPile, mcrotrNL1A_NL1B_NL2_NL3_InPile,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component NL1A_NL1B_NL2_NL3_InPile (without coords transformations) */
  mcJumpTrace_NL1A_NL1B_NL2_NL3_InPile:
  SIG_MESSAGE("NL1A_NL1B_NL2_NL3_InPile (Trace)");
  mcDEBUG_COMP("NL1A_NL1B_NL2_NL3_InPile")
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

#define mcabsorbComp mcabsorbCompNL1A_NL1B_NL2_NL3_InPile
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
#define mccompcurname  NL1A_NL1B_NL2_NL3_InPile
#define mccompcurtype  Guide
#define mccompcurindex 4
#define pTable mccNL1A_NL1B_NL2_NL3_InPile_pTable
{   /* Declarations of NL1A_NL1B_NL2_NL3_InPile=Guide() SETTING parameters. */
char* reflect = mccNL1A_NL1B_NL2_NL3_InPile_reflect;
MCNUM w1 = mccNL1A_NL1B_NL2_NL3_InPile_w1;
MCNUM h1 = mccNL1A_NL1B_NL2_NL3_InPile_h1;
MCNUM w2 = mccNL1A_NL1B_NL2_NL3_InPile_w2;
MCNUM h2 = mccNL1A_NL1B_NL2_NL3_InPile_h2;
MCNUM l = mccNL1A_NL1B_NL2_NL3_InPile_l;
MCNUM R0 = mccNL1A_NL1B_NL2_NL3_InPile_R0;
MCNUM Qc = mccNL1A_NL1B_NL2_NL3_InPile_Qc;
MCNUM alpha = mccNL1A_NL1B_NL2_NL3_InPile_alpha;
MCNUM m = mccNL1A_NL1B_NL2_NL3_InPile_m;
MCNUM W = mccNL1A_NL1B_NL2_NL3_InPile_W;
#line 94 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide.comp"
{
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
}
#line 13344 "./HZB_FLEX.c"
}   /* End of NL1A_NL1B_NL2_NL3_InPile=Guide() SETTING parameter declarations. */
#undef pTable
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompNL1A_NL1B_NL2_NL3_InPile:
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

  /* TRACE Component NL1B_Straight1_Entrance_Window [5] */
  mccoordschange(mcposrNL1B_Straight1_Entrance_Window, mcrotrNL1B_Straight1_Entrance_Window,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component NL1B_Straight1_Entrance_Window (without coords transformations) */
  mcJumpTrace_NL1B_Straight1_Entrance_Window:
  SIG_MESSAGE("NL1B_Straight1_Entrance_Window (Trace)");
  mcDEBUG_COMP("NL1B_Straight1_Entrance_Window")
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

#define mcabsorbComp mcabsorbCompNL1B_Straight1_Entrance_Window
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
#define mccompcurname  NL1B_Straight1_Entrance_Window
#define mccompcurtype  Arm
#define mccompcurindex 5
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompNL1B_Straight1_Entrance_Window:
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

  /* TRACE Component NL1B_Straight1 [6] */
  mccoordschange(mcposrNL1B_Straight1, mcrotrNL1B_Straight1,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component NL1B_Straight1 (without coords transformations) */
  mcJumpTrace_NL1B_Straight1:
  SIG_MESSAGE("NL1B_Straight1 (Trace)");
  mcDEBUG_COMP("NL1B_Straight1")
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

#define mcabsorbComp mcabsorbCompNL1B_Straight1
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
#define mccompcurname  NL1B_Straight1
#define mccompcurtype  Guide
#define mccompcurindex 6
#define pTable mccNL1B_Straight1_pTable
{   /* Declarations of NL1B_Straight1=Guide() SETTING parameters. */
char* reflect = mccNL1B_Straight1_reflect;
MCNUM w1 = mccNL1B_Straight1_w1;
MCNUM h1 = mccNL1B_Straight1_h1;
MCNUM w2 = mccNL1B_Straight1_w2;
MCNUM h2 = mccNL1B_Straight1_h2;
MCNUM l = mccNL1B_Straight1_l;
MCNUM R0 = mccNL1B_Straight1_R0;
MCNUM Qc = mccNL1B_Straight1_Qc;
MCNUM alpha = mccNL1B_Straight1_alpha;
MCNUM m = mccNL1B_Straight1_m;
MCNUM W = mccNL1B_Straight1_W;
#line 94 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide.comp"
{
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
}
#line 13673 "./HZB_FLEX.c"
}   /* End of NL1B_Straight1=Guide() SETTING parameter declarations. */
#undef pTable
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompNL1B_Straight1:
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

  /* TRACE Component NL1B_Curved_Entrance_Window [7] */
  mccoordschange(mcposrNL1B_Curved_Entrance_Window, mcrotrNL1B_Curved_Entrance_Window,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component NL1B_Curved_Entrance_Window (without coords transformations) */
  mcJumpTrace_NL1B_Curved_Entrance_Window:
  SIG_MESSAGE("NL1B_Curved_Entrance_Window (Trace)");
  mcDEBUG_COMP("NL1B_Curved_Entrance_Window")
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

#define mcabsorbComp mcabsorbCompNL1B_Curved_Entrance_Window
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
#define mccompcurname  NL1B_Curved_Entrance_Window
#define mccompcurtype  Arm
#define mccompcurindex 7
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompNL1B_Curved_Entrance_Window:
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

  /* TRACE Component NL1B_Curved_Guide_1 [8] */
  mccoordschange(mcposrNL1B_Curved_Guide_1, mcrotrNL1B_Curved_Guide_1,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component NL1B_Curved_Guide_1 (without coords transformations) */
  mcJumpTrace_NL1B_Curved_Guide_1:
  SIG_MESSAGE("NL1B_Curved_Guide_1 (Trace)");
  mcDEBUG_COMP("NL1B_Curved_Guide_1")
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

#define mcabsorbComp mcabsorbCompNL1B_Curved_Guide_1
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
#define mccompcurname  NL1B_Curved_Guide_1
#define mccompcurtype  Guide_curved
#define mccompcurindex 8
{   /* Declarations of NL1B_Curved_Guide_1=Guide_curved() SETTING parameters. */
MCNUM w1 = mccNL1B_Curved_Guide_1_w1;
MCNUM h1 = mccNL1B_Curved_Guide_1_h1;
MCNUM l = mccNL1B_Curved_Guide_1_l;
MCNUM R0 = mccNL1B_Curved_Guide_1_R0;
MCNUM Qc = mccNL1B_Curved_Guide_1_Qc;
MCNUM alpha = mccNL1B_Curved_Guide_1_alpha;
MCNUM m = mccNL1B_Curved_Guide_1_m;
MCNUM W = mccNL1B_Curved_Guide_1_W;
MCNUM curvature = mccNL1B_Curved_Guide_1_curvature;
#line 68 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../contrib/Guide_curved.comp"
{
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
}
#line 13956 "./HZB_FLEX.c"
}   /* End of NL1B_Curved_Guide_1=Guide_curved() SETTING parameter declarations. */
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompNL1B_Curved_Guide_1:
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

  /* TRACE Component NL1B_Velocity_Selector_Gap_Entrance_Window [9] */
  mccoordschange(mcposrNL1B_Velocity_Selector_Gap_Entrance_Window, mcrotrNL1B_Velocity_Selector_Gap_Entrance_Window,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component NL1B_Velocity_Selector_Gap_Entrance_Window (without coords transformations) */
  mcJumpTrace_NL1B_Velocity_Selector_Gap_Entrance_Window:
  SIG_MESSAGE("NL1B_Velocity_Selector_Gap_Entrance_Window (Trace)");
  mcDEBUG_COMP("NL1B_Velocity_Selector_Gap_Entrance_Window")
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

#define mcabsorbComp mcabsorbCompNL1B_Velocity_Selector_Gap_Entrance_Window
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
#define mccompcurname  NL1B_Velocity_Selector_Gap_Entrance_Window
#define mccompcurtype  Arm
#define mccompcurindex 9
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompNL1B_Velocity_Selector_Gap_Entrance_Window:
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

  /* TRACE Component Before_selec [10] */
  mccoordschange(mcposrBefore_selec, mcrotrBefore_selec,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component Before_selec (without coords transformations) */
  mcJumpTrace_Before_selec:
  SIG_MESSAGE("Before_selec (Trace)");
  mcDEBUG_COMP("Before_selec")
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

#define mcabsorbComp mcabsorbCompBefore_selec
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
#define mccompcurname  Before_selec
#define mccompcurtype  PSD_monitor
#define mccompcurindex 10
#define PSD_N mccBefore_selec_PSD_N
#define PSD_p mccBefore_selec_PSD_p
#define PSD_p2 mccBefore_selec_PSD_p2
{   /* Declarations of Before_selec=PSD_monitor() SETTING parameters. */
int nx = mccBefore_selec_nx;
int ny = mccBefore_selec_ny;
char* filename = mccBefore_selec_filename;
MCNUM xmin = mccBefore_selec_xmin;
MCNUM xmax = mccBefore_selec_xmax;
MCNUM ymin = mccBefore_selec_ymin;
MCNUM ymax = mccBefore_selec_ymax;
MCNUM xwidth = mccBefore_selec_xwidth;
MCNUM yheight = mccBefore_selec_yheight;
MCNUM restore_neutron = mccBefore_selec_restore_neutron;
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
#line 14193 "./HZB_FLEX.c"
}   /* End of Before_selec=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompBefore_selec:
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

  /* TRACE Component SELEC [11] */
  mccoordschange(mcposrSELEC, mcrotrSELEC,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component SELEC (without coords transformations) */
  mcJumpTrace_SELEC:
  SIG_MESSAGE("SELEC (Trace)");
  mcDEBUG_COMP("SELEC")
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

#define mcabsorbComp mcabsorbCompSELEC
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
#define mccompcurname  SELEC
#define mccompcurtype  Selector
#define mccompcurindex 11
{   /* Declarations of SELEC=Selector() SETTING parameters. */
MCNUM xmin = mccSELEC_xmin;
MCNUM xmax = mccSELEC_xmax;
MCNUM ymin = mccSELEC_ymin;
MCNUM ymax = mccSELEC_ymax;
MCNUM length = mccSELEC_length;
MCNUM xwidth = mccSELEC_xwidth;
MCNUM yheight = mccSELEC_yheight;
MCNUM nslit = mccSELEC_nslit;
MCNUM d = mccSELEC_d;
MCNUM radius = mccSELEC_radius;
MCNUM alpha = mccSELEC_alpha;
MCNUM nu = mccSELEC_nu;
#line 74 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Selector.comp"
{
  double E;
  double dt;
  double open_angle, closed_angle, n_angle_in, n_angle_out;
  double sel_phase, act_radius; // distance between neutron and selector axle

  PROP_Z0;
  E=VS2E*(vx*vx+vy*vy+vz*vz);
  if (x<xmin || x>xmax || y<ymin || y>ymax)
    ABSORB;                                  /* because outside frame */
  dt = length/vz;
  /* get phase angle of selector rotor as MonteCarlo choice
   only the free space between two neighboring spokes is taken
   p is adjusted to transmission for parallel beam
  */
  n_angle_in = atan2( x,y+radius)*RAD2DEG;
  act_radius = sqrt(x*x + pow(radius+y,2));
  closed_angle = d/act_radius*RAD2DEG;
  open_angle = 360/nslit-closed_angle;
  sel_phase = open_angle*rand01();
  p *= (open_angle/(closed_angle+open_angle));

  PROP_DT(dt);

  if (x<xmin || x>xmax || y<ymin || y>ymax)
    ABSORB;                                  /* because outside frame */
  /* now let's look whether the neutron is still
 between the same two spokes or absorbed meanwhile */

  n_angle_out = atan2(x,y+radius)*RAD2DEG;  /* neutron beam might be divergent */

  sel_phase = sel_phase + nu*dt*360 - alpha;  /* rotor turned, but spokes are torsaded */

  if (n_angle_out<(n_angle_in-sel_phase) || n_angle_out>(n_angle_in-sel_phase+open_angle) )
    ABSORB;              /* because must have passed absorber */
  else
    SCATTER;
}
#line 14353 "./HZB_FLEX.c"
}   /* End of SELEC=Selector() SETTING parameter declarations. */
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompSELEC:
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

  /* TRACE Component After_selec [12] */
  mccoordschange(mcposrAfter_selec, mcrotrAfter_selec,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component After_selec (without coords transformations) */
  mcJumpTrace_After_selec:
  SIG_MESSAGE("After_selec (Trace)");
  mcDEBUG_COMP("After_selec")
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

#define mcabsorbComp mcabsorbCompAfter_selec
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
#define mccompcurname  After_selec
#define mccompcurtype  PSD_monitor
#define mccompcurindex 12
#define PSD_N mccAfter_selec_PSD_N
#define PSD_p mccAfter_selec_PSD_p
#define PSD_p2 mccAfter_selec_PSD_p2
{   /* Declarations of After_selec=PSD_monitor() SETTING parameters. */
int nx = mccAfter_selec_nx;
int ny = mccAfter_selec_ny;
char* filename = mccAfter_selec_filename;
MCNUM xmin = mccAfter_selec_xmin;
MCNUM xmax = mccAfter_selec_xmax;
MCNUM ymin = mccAfter_selec_ymin;
MCNUM ymax = mccAfter_selec_ymax;
MCNUM xwidth = mccAfter_selec_xwidth;
MCNUM yheight = mccAfter_selec_yheight;
MCNUM restore_neutron = mccAfter_selec_restore_neutron;
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
#line 14487 "./HZB_FLEX.c"
}   /* End of After_selec=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompAfter_selec:
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

  /* TRACE Component NL1B_Curved_2_Entrance_Window [13] */
  mccoordschange(mcposrNL1B_Curved_2_Entrance_Window, mcrotrNL1B_Curved_2_Entrance_Window,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component NL1B_Curved_2_Entrance_Window (without coords transformations) */
  mcJumpTrace_NL1B_Curved_2_Entrance_Window:
  SIG_MESSAGE("NL1B_Curved_2_Entrance_Window (Trace)");
  mcDEBUG_COMP("NL1B_Curved_2_Entrance_Window")
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

#define mcabsorbComp mcabsorbCompNL1B_Curved_2_Entrance_Window
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
#define mccompcurname  NL1B_Curved_2_Entrance_Window
#define mccompcurtype  Arm
#define mccompcurindex 13
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompNL1B_Curved_2_Entrance_Window:
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

  /* TRACE Component NL1B_Curved_Guide_2 [14] */
  mccoordschange(mcposrNL1B_Curved_Guide_2, mcrotrNL1B_Curved_Guide_2,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component NL1B_Curved_Guide_2 (without coords transformations) */
  mcJumpTrace_NL1B_Curved_Guide_2:
  SIG_MESSAGE("NL1B_Curved_Guide_2 (Trace)");
  mcDEBUG_COMP("NL1B_Curved_Guide_2")
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

#define mcabsorbComp mcabsorbCompNL1B_Curved_Guide_2
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
#define mccompcurname  NL1B_Curved_Guide_2
#define mccompcurtype  Guide_curved
#define mccompcurindex 14
{   /* Declarations of NL1B_Curved_Guide_2=Guide_curved() SETTING parameters. */
MCNUM w1 = mccNL1B_Curved_Guide_2_w1;
MCNUM h1 = mccNL1B_Curved_Guide_2_h1;
MCNUM l = mccNL1B_Curved_Guide_2_l;
MCNUM R0 = mccNL1B_Curved_Guide_2_R0;
MCNUM Qc = mccNL1B_Curved_Guide_2_Qc;
MCNUM alpha = mccNL1B_Curved_Guide_2_alpha;
MCNUM m = mccNL1B_Curved_Guide_2_m;
MCNUM W = mccNL1B_Curved_Guide_2_W;
MCNUM curvature = mccNL1B_Curved_Guide_2_curvature;
#line 68 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../contrib/Guide_curved.comp"
{
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
}
#line 14772 "./HZB_FLEX.c"
}   /* End of NL1B_Curved_Guide_2=Guide_curved() SETTING parameter declarations. */
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompNL1B_Curved_Guide_2:
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

  /* TRACE Component NL1B_Straight2_Entrance_Window [15] */
  mccoordschange(mcposrNL1B_Straight2_Entrance_Window, mcrotrNL1B_Straight2_Entrance_Window,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component NL1B_Straight2_Entrance_Window (without coords transformations) */
  mcJumpTrace_NL1B_Straight2_Entrance_Window:
  SIG_MESSAGE("NL1B_Straight2_Entrance_Window (Trace)");
  mcDEBUG_COMP("NL1B_Straight2_Entrance_Window")
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

#define mcabsorbComp mcabsorbCompNL1B_Straight2_Entrance_Window
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
#define mccompcurname  NL1B_Straight2_Entrance_Window
#define mccompcurtype  Arm
#define mccompcurindex 15
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompNL1B_Straight2_Entrance_Window:
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

  /* TRACE Component NL1B_Straight2 [16] */
  mccoordschange(mcposrNL1B_Straight2, mcrotrNL1B_Straight2,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component NL1B_Straight2 (without coords transformations) */
  mcJumpTrace_NL1B_Straight2:
  SIG_MESSAGE("NL1B_Straight2 (Trace)");
  mcDEBUG_COMP("NL1B_Straight2")
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

#define mcabsorbComp mcabsorbCompNL1B_Straight2
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
#define mccompcurname  NL1B_Straight2
#define mccompcurtype  Guide
#define mccompcurindex 16
#define pTable mccNL1B_Straight2_pTable
{   /* Declarations of NL1B_Straight2=Guide() SETTING parameters. */
char* reflect = mccNL1B_Straight2_reflect;
MCNUM w1 = mccNL1B_Straight2_w1;
MCNUM h1 = mccNL1B_Straight2_h1;
MCNUM w2 = mccNL1B_Straight2_w2;
MCNUM h2 = mccNL1B_Straight2_h2;
MCNUM l = mccNL1B_Straight2_l;
MCNUM R0 = mccNL1B_Straight2_R0;
MCNUM Qc = mccNL1B_Straight2_Qc;
MCNUM alpha = mccNL1B_Straight2_alpha;
MCNUM m = mccNL1B_Straight2_m;
MCNUM W = mccNL1B_Straight2_W;
#line 94 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide.comp"
{
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
}
#line 15100 "./HZB_FLEX.c"
}   /* End of NL1B_Straight2=Guide() SETTING parameter declarations. */
#undef pTable
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompNL1B_Straight2:
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

  /* TRACE Component NL1B_Elliptical_Entrance_Window [17] */
  mccoordschange(mcposrNL1B_Elliptical_Entrance_Window, mcrotrNL1B_Elliptical_Entrance_Window,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component NL1B_Elliptical_Entrance_Window (without coords transformations) */
  mcJumpTrace_NL1B_Elliptical_Entrance_Window:
  SIG_MESSAGE("NL1B_Elliptical_Entrance_Window (Trace)");
  mcDEBUG_COMP("NL1B_Elliptical_Entrance_Window")
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

#define mcabsorbComp mcabsorbCompNL1B_Elliptical_Entrance_Window
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
#define mccompcurname  NL1B_Elliptical_Entrance_Window
#define mccompcurtype  Arm
#define mccompcurindex 17
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompNL1B_Elliptical_Entrance_Window:
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

  /* TRACE Component elliptical_piece [18] */
  mccoordschange(mcposrelliptical_piece, mcrotrelliptical_piece,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component elliptical_piece (without coords transformations) */
  mcJumpTrace_elliptical_piece:
  SIG_MESSAGE("elliptical_piece (Trace)");
  mcDEBUG_COMP("elliptical_piece")
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

#define mcabsorbComp mcabsorbCompelliptical_piece
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
#define mccompcurname  elliptical_piece
#define mccompcurtype  Guide_tapering
#define mccompcurindex 18
#define w1c mccelliptical_piece_w1c
#define w2c mccelliptical_piece_w2c
#define ww mccelliptical_piece_ww
#define hh mccelliptical_piece_hh
#define whalf mccelliptical_piece_whalf
#define hhalf mccelliptical_piece_hhalf
#define lwhalf mccelliptical_piece_lwhalf
#define lhhalf mccelliptical_piece_lhhalf
#define h1_in mccelliptical_piece_h1_in
#define h2_out mccelliptical_piece_h2_out
#define w1_in mccelliptical_piece_w1_in
#define w2_out mccelliptical_piece_w2_out
#define l_seg mccelliptical_piece_l_seg
#define seg mccelliptical_piece_seg
#define h12 mccelliptical_piece_h12
#define h2 mccelliptical_piece_h2
#define w12 mccelliptical_piece_w12
#define w2 mccelliptical_piece_w2
#define a_ell_q mccelliptical_piece_a_ell_q
#define b_ell_q mccelliptical_piece_b_ell_q
#define lbw mccelliptical_piece_lbw
#define lbh mccelliptical_piece_lbh
#define mxi mccelliptical_piece_mxi
#define u1 mccelliptical_piece_u1
#define u2 mccelliptical_piece_u2
#define div1 mccelliptical_piece_div1
#define p2_para mccelliptical_piece_p2_para
#define test mccelliptical_piece_test
#define Div1 mccelliptical_piece_Div1
#define i mccelliptical_piece_i
#define ii mccelliptical_piece_ii
#define seg mccelliptical_piece_seg
#define fu mccelliptical_piece_fu
#define pos mccelliptical_piece_pos
#define file_name mccelliptical_piece_file_name
#define ep mccelliptical_piece_ep
#define num mccelliptical_piece_num
#define rotation_h mccelliptical_piece_rotation_h
#define rotation_v mccelliptical_piece_rotation_v
{   /* Declarations of elliptical_piece=Guide_tapering() SETTING parameters. */
char* option = mccelliptical_piece_option;
MCNUM w1 = mccelliptical_piece_w1;
MCNUM h1 = mccelliptical_piece_h1;
MCNUM l = mccelliptical_piece_l;
MCNUM linw = mccelliptical_piece_linw;
MCNUM loutw = mccelliptical_piece_loutw;
MCNUM linh = mccelliptical_piece_linh;
MCNUM louth = mccelliptical_piece_louth;
MCNUM R0 = mccelliptical_piece_R0;
MCNUM Qcx = mccelliptical_piece_Qcx;
MCNUM Qcy = mccelliptical_piece_Qcy;
MCNUM alphax = mccelliptical_piece_alphax;
MCNUM alphay = mccelliptical_piece_alphay;
MCNUM W = mccelliptical_piece_W;
MCNUM mx = mccelliptical_piece_mx;
MCNUM my = mccelliptical_piece_my;
MCNUM segno = mccelliptical_piece_segno;
MCNUM curvature = mccelliptical_piece_curvature;
MCNUM curvature_v = mccelliptical_piece_curvature_v;
#line 453 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide_tapering.comp"
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
#line 15523 "./HZB_FLEX.c"
}   /* End of elliptical_piece=Guide_tapering() SETTING parameter declarations. */
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
  mcabsorbCompelliptical_piece:
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

  /* TRACE Component Virtual_source [19] */
  mccoordschange(mcposrVirtual_source, mcrotrVirtual_source,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component Virtual_source (without coords transformations) */
  mcJumpTrace_Virtual_source:
  SIG_MESSAGE("Virtual_source (Trace)");
  mcDEBUG_COMP("Virtual_source")
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

#define mcabsorbComp mcabsorbCompVirtual_source
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
#define mccompcurname  Virtual_source
#define mccompcurtype  Slit
#define mccompcurindex 19
{   /* Declarations of Virtual_source=Slit() SETTING parameters. */
MCNUM xmin = mccVirtual_source_xmin;
MCNUM xmax = mccVirtual_source_xmax;
MCNUM ymin = mccVirtual_source_ymin;
MCNUM ymax = mccVirtual_source_ymax;
MCNUM radius = mccVirtual_source_radius;
MCNUM xwidth = mccVirtual_source_xwidth;
MCNUM yheight = mccVirtual_source_yheight;
#line 71 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Slit.comp"
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
#line 15686 "./HZB_FLEX.c"
}   /* End of Virtual_source=Slit() SETTING parameter declarations. */
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompVirtual_source:
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

  /* TRACE Component NL1B_Vertical_Guide_Entrance_Window [20] */
  mccoordschange(mcposrNL1B_Vertical_Guide_Entrance_Window, mcrotrNL1B_Vertical_Guide_Entrance_Window,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component NL1B_Vertical_Guide_Entrance_Window (without coords transformations) */
  mcJumpTrace_NL1B_Vertical_Guide_Entrance_Window:
  SIG_MESSAGE("NL1B_Vertical_Guide_Entrance_Window (Trace)");
  mcDEBUG_COMP("NL1B_Vertical_Guide_Entrance_Window")
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

#define mcabsorbComp mcabsorbCompNL1B_Vertical_Guide_Entrance_Window
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
#define mccompcurname  NL1B_Vertical_Guide_Entrance_Window
#define mccompcurtype  Arm
#define mccompcurindex 20
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompNL1B_Vertical_Guide_Entrance_Window:
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

  /* TRACE Component NL1B_Vertical_Guide [21] */
  mccoordschange(mcposrNL1B_Vertical_Guide, mcrotrNL1B_Vertical_Guide,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component NL1B_Vertical_Guide (without coords transformations) */
  mcJumpTrace_NL1B_Vertical_Guide:
  SIG_MESSAGE("NL1B_Vertical_Guide (Trace)");
  mcDEBUG_COMP("NL1B_Vertical_Guide")
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

#define mcabsorbComp mcabsorbCompNL1B_Vertical_Guide
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
#define mccompcurname  NL1B_Vertical_Guide
#define mccompcurtype  Guide_channeled
#define mccompcurindex 21
#define w1c mccNL1B_Vertical_Guide_w1c
#define w2c mccNL1B_Vertical_Guide_w2c
#define ww mccNL1B_Vertical_Guide_ww
#define hh mccNL1B_Vertical_Guide_hh
#define whalf mccNL1B_Vertical_Guide_whalf
#define hhalf mccNL1B_Vertical_Guide_hhalf
#define lwhalf mccNL1B_Vertical_Guide_lwhalf
#define lhhalf mccNL1B_Vertical_Guide_lhhalf
{   /* Declarations of NL1B_Vertical_Guide=Guide_channeled() SETTING parameters. */
MCNUM w1 = mccNL1B_Vertical_Guide_w1;
MCNUM h1 = mccNL1B_Vertical_Guide_h1;
MCNUM w2 = mccNL1B_Vertical_Guide_w2;
MCNUM h2 = mccNL1B_Vertical_Guide_h2;
MCNUM l = mccNL1B_Vertical_Guide_l;
MCNUM R0 = mccNL1B_Vertical_Guide_R0;
MCNUM Qc = mccNL1B_Vertical_Guide_Qc;
MCNUM alpha = mccNL1B_Vertical_Guide_alpha;
MCNUM m = mccNL1B_Vertical_Guide_m;
MCNUM nslit = mccNL1B_Vertical_Guide_nslit;
MCNUM d = mccNL1B_Vertical_Guide_d;
MCNUM Qcx = mccNL1B_Vertical_Guide_Qcx;
MCNUM Qcy = mccNL1B_Vertical_Guide_Qcy;
MCNUM alphax = mccNL1B_Vertical_Guide_alphax;
MCNUM alphay = mccNL1B_Vertical_Guide_alphay;
MCNUM W = mccNL1B_Vertical_Guide_W;
MCNUM mx = mccNL1B_Vertical_Guide_mx;
MCNUM my = mccNL1B_Vertical_Guide_my;
MCNUM nu = mccNL1B_Vertical_Guide_nu;
MCNUM phase = mccNL1B_Vertical_Guide_phase;
#line 131 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide_channeled.comp"
{
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
}
#line 16088 "./HZB_FLEX.c"
}   /* End of NL1B_Vertical_Guide=Guide_channeled() SETTING parameter declarations. */
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
  mcabsorbCompNL1B_Vertical_Guide:
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

  /* TRACE Component NL1B_Guide_Exit [22] */
  mccoordschange(mcposrNL1B_Guide_Exit, mcrotrNL1B_Guide_Exit,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component NL1B_Guide_Exit (without coords transformations) */
  mcJumpTrace_NL1B_Guide_Exit:
  SIG_MESSAGE("NL1B_Guide_Exit (Trace)");
  mcDEBUG_COMP("NL1B_Guide_Exit")
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

#define mcabsorbComp mcabsorbCompNL1B_Guide_Exit
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
#define mccompcurname  NL1B_Guide_Exit
#define mccompcurtype  Arm
#define mccompcurindex 22
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompNL1B_Guide_Exit:
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

  /* TRACE Component energy_endguide [23] */
  mccoordschange(mcposrenergy_endguide, mcrotrenergy_endguide,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component energy_endguide (without coords transformations) */
  mcJumpTrace_energy_endguide:
  SIG_MESSAGE("energy_endguide (Trace)");
  mcDEBUG_COMP("energy_endguide")
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

#define mcabsorbComp mcabsorbCompenergy_endguide
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
#define mccompcurname  energy_endguide
#define mccompcurtype  E_monitor
#define mccompcurindex 23
#define nE mccenergy_endguide_nE
#define E_N mccenergy_endguide_E_N
#define E_p mccenergy_endguide_E_p
#define E_p2 mccenergy_endguide_E_p2
#define S_p mccenergy_endguide_S_p
#define S_pE mccenergy_endguide_S_pE
#define S_pE2 mccenergy_endguide_S_pE2
{   /* Declarations of energy_endguide=E_monitor() SETTING parameters. */
char* filename = mccenergy_endguide_filename;
MCNUM xmin = mccenergy_endguide_xmin;
MCNUM xmax = mccenergy_endguide_xmax;
MCNUM ymin = mccenergy_endguide_ymin;
MCNUM ymax = mccenergy_endguide_ymax;
MCNUM xwidth = mccenergy_endguide_xwidth;
MCNUM yheight = mccenergy_endguide_yheight;
MCNUM Emin = mccenergy_endguide_Emin;
MCNUM Emax = mccenergy_endguide_Emax;
MCNUM restore_neutron = mccenergy_endguide_restore_neutron;
int nowritefile = mccenergy_endguide_nowritefile;
#line 89 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/E_monitor.comp"
{
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
}
#line 16350 "./HZB_FLEX.c"
}   /* End of energy_endguide=E_monitor() SETTING parameter declarations. */
#undef S_pE2
#undef S_pE
#undef S_p
#undef E_p2
#undef E_p
#undef E_N
#undef nE
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompenergy_endguide:
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

  /* TRACE Component Mono_center [24] */
  mccoordschange(mcposrMono_center, mcrotrMono_center,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component Mono_center (without coords transformations) */
  mcJumpTrace_Mono_center:
  SIG_MESSAGE("Mono_center (Trace)");
  mcDEBUG_COMP("Mono_center")
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

#define mcabsorbComp mcabsorbCompMono_center
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
#define mccompcurname  Mono_center
#define mccompcurtype  Arm
#define mccompcurindex 24
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompMono_center:
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

  /* TRACE Component Monochromator [25] */
  mccoordschange(mcposrMonochromator, mcrotrMonochromator,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component Monochromator (without coords transformations) */
  mcJumpTrace_Monochromator:
  SIG_MESSAGE("Monochromator (Trace)");
  mcDEBUG_COMP("Monochromator")
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

#define mcabsorbComp mcabsorbCompMonochromator
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
#define mccompcurname  Monochromator
#define mccompcurtype  Monochromator_curved
#define mccompcurindex 25
#define mos_rms_y mccMonochromator_mos_rms_y
#define mos_rms_z mccMonochromator_mos_rms_z
#define mos_rms_max mccMonochromator_mos_rms_max
#define mono_Q mccMonochromator_mono_Q
#define SlabWidth mccMonochromator_SlabWidth
#define SlabHeight mccMonochromator_SlabHeight
#define rTable mccMonochromator_rTable
#define tTable mccMonochromator_tTable
#define row mccMonochromator_row
#define col mccMonochromator_col
#define tiltH mccMonochromator_tiltH
#define tiltV mccMonochromator_tiltV
{   /* Declarations of Monochromator=Monochromator_curved() SETTING parameters. */
char* reflect = mccMonochromator_reflect;
char* transmit = mccMonochromator_transmit;
MCNUM zwidth = mccMonochromator_zwidth;
MCNUM yheight = mccMonochromator_yheight;
MCNUM gap = mccMonochromator_gap;
MCNUM NH = mccMonochromator_NH;
MCNUM NV = mccMonochromator_NV;
MCNUM mosaich = mccMonochromator_mosaich;
MCNUM mosaicv = mccMonochromator_mosaicv;
MCNUM r0 = mccMonochromator_r0;
MCNUM t0 = mccMonochromator_t0;
MCNUM Q = mccMonochromator_Q;
MCNUM RV = mccMonochromator_RV;
MCNUM RH = mccMonochromator_RH;
MCNUM DM = mccMonochromator_DM;
MCNUM mosaic = mccMonochromator_mosaic;
MCNUM width = mccMonochromator_width;
MCNUM height = mccMonochromator_height;
MCNUM verbose = mccMonochromator_verbose;
MCNUM order = mccMonochromator_order;
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
#line 16827 "./HZB_FLEX.c"
}   /* End of Monochromator=Monochromator_curved() SETTING parameter declarations. */
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
  mcabsorbCompMonochromator:
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

  /* TRACE Component Mono_sample_arm [26] */
  mccoordschange(mcposrMono_sample_arm, mcrotrMono_sample_arm,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component Mono_sample_arm (without coords transformations) */
  mcJumpTrace_Mono_sample_arm:
  SIG_MESSAGE("Mono_sample_arm (Trace)");
  mcDEBUG_COMP("Mono_sample_arm")
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

#define mcabsorbComp mcabsorbCompMono_sample_arm
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
#define mccompcurname  Mono_sample_arm
#define mccompcurtype  Arm
#define mccompcurindex 26
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompMono_sample_arm:
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

  /* TRACE Component energy_mono [27] */
  mccoordschange(mcposrenergy_mono, mcrotrenergy_mono,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component energy_mono (without coords transformations) */
  mcJumpTrace_energy_mono:
  SIG_MESSAGE("energy_mono (Trace)");
  mcDEBUG_COMP("energy_mono")
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

#define mcabsorbComp mcabsorbCompenergy_mono
  STORE_NEUTRON(27,
    mcnlx,
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
  mcNCounter[27]++;
  mcPCounter[27] += p;
  mcP2Counter[27] += p*p;
#define mccompcurname  energy_mono
#define mccompcurtype  E_monitor
#define mccompcurindex 27
#define nE mccenergy_mono_nE
#define E_N mccenergy_mono_E_N
#define E_p mccenergy_mono_E_p
#define E_p2 mccenergy_mono_E_p2
#define S_p mccenergy_mono_S_p
#define S_pE mccenergy_mono_S_pE
#define S_pE2 mccenergy_mono_S_pE2
{   /* Declarations of energy_mono=E_monitor() SETTING parameters. */
char* filename = mccenergy_mono_filename;
MCNUM xmin = mccenergy_mono_xmin;
MCNUM xmax = mccenergy_mono_xmax;
MCNUM ymin = mccenergy_mono_ymin;
MCNUM ymax = mccenergy_mono_ymax;
MCNUM xwidth = mccenergy_mono_xwidth;
MCNUM yheight = mccenergy_mono_yheight;
MCNUM Emin = mccenergy_mono_Emin;
MCNUM Emax = mccenergy_mono_Emax;
MCNUM restore_neutron = mccenergy_mono_restore_neutron;
int nowritefile = mccenergy_mono_nowritefile;
#line 89 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/E_monitor.comp"
{
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
}
#line 17093 "./HZB_FLEX.c"
}   /* End of energy_mono=E_monitor() SETTING parameter declarations. */
#undef S_pE2
#undef S_pE
#undef S_p
#undef E_p2
#undef E_p
#undef E_N
#undef nE
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompenergy_mono:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(27,
      mcnlx,
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

  /* TRACE Component energy_pre_sample [28] */
  mccoordschange(mcposrenergy_pre_sample, mcrotrenergy_pre_sample,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component energy_pre_sample (without coords transformations) */
  mcJumpTrace_energy_pre_sample:
  SIG_MESSAGE("energy_pre_sample (Trace)");
  mcDEBUG_COMP("energy_pre_sample")
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

#define mcabsorbComp mcabsorbCompenergy_pre_sample
  STORE_NEUTRON(28,
    mcnlx,
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
  mcNCounter[28]++;
  mcPCounter[28] += p;
  mcP2Counter[28] += p*p;
#define mccompcurname  energy_pre_sample
#define mccompcurtype  E_monitor
#define mccompcurindex 28
#define nE mccenergy_pre_sample_nE
#define E_N mccenergy_pre_sample_E_N
#define E_p mccenergy_pre_sample_E_p
#define E_p2 mccenergy_pre_sample_E_p2
#define S_p mccenergy_pre_sample_S_p
#define S_pE mccenergy_pre_sample_S_pE
#define S_pE2 mccenergy_pre_sample_S_pE2
{   /* Declarations of energy_pre_sample=E_monitor() SETTING parameters. */
char* filename = mccenergy_pre_sample_filename;
MCNUM xmin = mccenergy_pre_sample_xmin;
MCNUM xmax = mccenergy_pre_sample_xmax;
MCNUM ymin = mccenergy_pre_sample_ymin;
MCNUM ymax = mccenergy_pre_sample_ymax;
MCNUM xwidth = mccenergy_pre_sample_xwidth;
MCNUM yheight = mccenergy_pre_sample_yheight;
MCNUM Emin = mccenergy_pre_sample_Emin;
MCNUM Emax = mccenergy_pre_sample_Emax;
MCNUM restore_neutron = mccenergy_pre_sample_restore_neutron;
int nowritefile = mccenergy_pre_sample_nowritefile;
#line 89 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/E_monitor.comp"
{
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
}
#line 17251 "./HZB_FLEX.c"
}   /* End of energy_pre_sample=E_monitor() SETTING parameter declarations. */
#undef S_pE2
#undef S_pE
#undef S_p
#undef E_p2
#undef E_p
#undef E_N
#undef nE
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompenergy_pre_sample:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(28,
      mcnlx,
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

  /* TRACE Component Sample_center [29] */
  mccoordschange(mcposrSample_center, mcrotrSample_center,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component Sample_center (without coords transformations) */
  mcJumpTrace_Sample_center:
  SIG_MESSAGE("Sample_center (Trace)");
  mcDEBUG_COMP("Sample_center")
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

#define mcabsorbComp mcabsorbCompSample_center
  STORE_NEUTRON(29,
    mcnlx,
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
  mcNCounter[29]++;
  mcPCounter[29] += p;
  mcP2Counter[29] += p*p;
#define mccompcurname  Sample_center
#define mccompcurtype  Arm
#define mccompcurindex 29
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompSample_center:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(29,
      mcnlx,
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

  /* TRACE Component Sample_analyser_arm [30] */
  mccoordschange(mcposrSample_analyser_arm, mcrotrSample_analyser_arm,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component Sample_analyser_arm (without coords transformations) */
  mcJumpTrace_Sample_analyser_arm:
  SIG_MESSAGE("Sample_analyser_arm (Trace)");
  mcDEBUG_COMP("Sample_analyser_arm")
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

#define mcabsorbComp mcabsorbCompSample_analyser_arm
  STORE_NEUTRON(30,
    mcnlx,
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
  mcNCounter[30]++;
  mcPCounter[30] += p;
  mcP2Counter[30] += p*p;
#define mccompcurname  Sample_analyser_arm
#define mccompcurtype  Arm
#define mccompcurindex 30
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompSample_analyser_arm:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(30,
      mcnlx,
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

  /* TRACE Component div_mono [31] */
  mccoordschange(mcposrdiv_mono, mcrotrdiv_mono,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component div_mono (without coords transformations) */
  mcJumpTrace_div_mono:
  SIG_MESSAGE("div_mono (Trace)");
  mcDEBUG_COMP("div_mono")
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

#define mcabsorbComp mcabsorbCompdiv_mono
  STORE_NEUTRON(31,
    mcnlx,
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
  mcNCounter[31]++;
  mcPCounter[31] += p;
  mcP2Counter[31] += p*p;
#define mccompcurname  div_mono
#define mccompcurtype  Divergence_monitor
#define mccompcurindex 31
#define nh mccdiv_mono_nh
#define nv mccdiv_mono_nv
#define Div_N mccdiv_mono_Div_N
#define Div_p mccdiv_mono_Div_p
#define Div_p2 mccdiv_mono_Div_p2
{   /* Declarations of div_mono=Divergence_monitor() SETTING parameters. */
char* filename = mccdiv_mono_filename;
MCNUM xmin = mccdiv_mono_xmin;
MCNUM xmax = mccdiv_mono_xmax;
MCNUM ymin = mccdiv_mono_ymin;
MCNUM ymax = mccdiv_mono_ymax;
MCNUM xwidth = mccdiv_mono_xwidth;
MCNUM yheight = mccdiv_mono_yheight;
MCNUM maxdiv_h = mccdiv_mono_maxdiv_h;
MCNUM maxdiv_v = mccdiv_mono_maxdiv_v;
MCNUM restore_neutron = mccdiv_mono_restore_neutron;
MCNUM nx = mccdiv_mono_nx;
MCNUM ny = mccdiv_mono_ny;
MCNUM nz = mccdiv_mono_nz;
int nowritefile = mccdiv_mono_nowritefile;
#line 89 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/Divergence_monitor.comp"
{
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
}
#line 17617 "./HZB_FLEX.c"
}   /* End of div_mono=Divergence_monitor() SETTING parameter declarations. */
#undef Div_p2
#undef Div_p
#undef Div_N
#undef nv
#undef nh
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompdiv_mono:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(31,
      mcnlx,
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

  /* TRACE Component div_mono_H [32] */
  mccoordschange(mcposrdiv_mono_H, mcrotrdiv_mono_H,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component div_mono_H (without coords transformations) */
  mcJumpTrace_div_mono_H:
  SIG_MESSAGE("div_mono_H (Trace)");
  mcDEBUG_COMP("div_mono_H")
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

#define mcabsorbComp mcabsorbCompdiv_mono_H
  STORE_NEUTRON(32,
    mcnlx,
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
  mcNCounter[32]++;
  mcPCounter[32] += p;
  mcP2Counter[32] += p*p;
#define mccompcurname  div_mono_H
#define mccompcurtype  Hdiv_monitor
#define mccompcurindex 32
#define nh mccdiv_mono_H_nh
#define Div_N mccdiv_mono_H_Div_N
#define Div_p mccdiv_mono_H_Div_p
#define Div_p2 mccdiv_mono_H_Div_p2
{   /* Declarations of div_mono_H=Hdiv_monitor() SETTING parameters. */
char* filename = mccdiv_mono_H_filename;
MCNUM xmin = mccdiv_mono_H_xmin;
MCNUM xmax = mccdiv_mono_H_xmax;
MCNUM ymin = mccdiv_mono_H_ymin;
MCNUM ymax = mccdiv_mono_H_ymax;
MCNUM xwidth = mccdiv_mono_H_xwidth;
MCNUM yheight = mccdiv_mono_H_yheight;
MCNUM h_maxdiv = mccdiv_mono_H_h_maxdiv;
MCNUM restore_neutron = mccdiv_mono_H_restore_neutron;
int nowritefile = mccdiv_mono_H_nowritefile;
#line 85 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/Hdiv_monitor.comp"
{
    int i;
    double h_div;

    PROP_Z0;
    if (x>xmin && x<xmax && y>ymin && y<ymax)
    {
      h_div = RAD2DEG*atan2(vx,vz);
      if (h_div < (double)h_maxdiv && h_div > -(double)h_maxdiv)
      {
        i = floor((h_div + (double)h_maxdiv)*nh/(2.0*(double)h_maxdiv));
        Div_N[i]++;
        Div_p[i] += p;
        Div_p2[i] += p*p;
        SCATTER;
      }
    }
    if (restore_neutron) {
      RESTORE_NEUTRON(INDEX_CURRENT_COMP, x, y, z, vx, vy, vz, t, sx, sy, sz, p);
    }
}
#line 17764 "./HZB_FLEX.c"
}   /* End of div_mono_H=Hdiv_monitor() SETTING parameter declarations. */
#undef Div_p2
#undef Div_p
#undef Div_N
#undef nh
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompdiv_mono_H:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(32,
      mcnlx,
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

  /* TRACE Component psd_sam [33] */
  mccoordschange(mcposrpsd_sam, mcrotrpsd_sam,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component psd_sam (without coords transformations) */
  mcJumpTrace_psd_sam:
  SIG_MESSAGE("psd_sam (Trace)");
  mcDEBUG_COMP("psd_sam")
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

#define mcabsorbComp mcabsorbComppsd_sam
  STORE_NEUTRON(33,
    mcnlx,
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
  mcNCounter[33]++;
  mcPCounter[33] += p;
  mcP2Counter[33] += p*p;
#define mccompcurname  psd_sam
#define mccompcurtype  PSD_monitor
#define mccompcurindex 33
#define PSD_N mccpsd_sam_PSD_N
#define PSD_p mccpsd_sam_PSD_p
#define PSD_p2 mccpsd_sam_PSD_p2
{   /* Declarations of psd_sam=PSD_monitor() SETTING parameters. */
int nx = mccpsd_sam_nx;
int ny = mccpsd_sam_ny;
char* filename = mccpsd_sam_filename;
MCNUM xmin = mccpsd_sam_xmin;
MCNUM xmax = mccpsd_sam_xmax;
MCNUM ymin = mccpsd_sam_ymin;
MCNUM ymax = mccpsd_sam_ymax;
MCNUM xwidth = mccpsd_sam_xwidth;
MCNUM yheight = mccpsd_sam_yheight;
MCNUM restore_neutron = mccpsd_sam_restore_neutron;
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
#line 17902 "./HZB_FLEX.c"
}   /* End of psd_sam=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbComppsd_sam:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(33,
      mcnlx,
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

  /* TRACE Component 1dpsd [34] */
  mccoordschange(mcposr1dpsd, mcrotr1dpsd,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component 1dpsd (without coords transformations) */
  mcJumpTrace_1dpsd:
  SIG_MESSAGE("1dpsd (Trace)");
  mcDEBUG_COMP("1dpsd")
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

#define mcabsorbComp mcabsorbComp1dpsd
  STORE_NEUTRON(34,
    mcnlx,
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
  mcNCounter[34]++;
  mcPCounter[34] += p;
  mcP2Counter[34] += p*p;
#define mccompcurname  1dpsd
#define mccompcurtype  PSDlin_monitor
#define mccompcurindex 34
#define nx mcc1dpsd_nx
#define PSDlin_N mcc1dpsd_PSDlin_N
#define PSDlin_p mcc1dpsd_PSDlin_p
#define PSDlin_p2 mcc1dpsd_PSDlin_p2
{   /* Declarations of 1dpsd=PSDlin_monitor() SETTING parameters. */
char* filename = mcc1dpsd_filename;
MCNUM xmin = mcc1dpsd_xmin;
MCNUM xmax = mcc1dpsd_xmax;
MCNUM ymin = mcc1dpsd_ymin;
MCNUM ymax = mcc1dpsd_ymax;
MCNUM xwidth = mcc1dpsd_xwidth;
MCNUM yheight = mcc1dpsd_yheight;
MCNUM restore_neutron = mcc1dpsd_restore_neutron;
int nowritefile = mcc1dpsd_nowritefile;
#line 81 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/PSDlin_monitor.comp"
{
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
}
#line 18045 "./HZB_FLEX.c"
}   /* End of 1dpsd=PSDlin_monitor() SETTING parameter declarations. */
#undef PSDlin_p2
#undef PSDlin_p
#undef PSDlin_N
#undef nx
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbComp1dpsd:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(34,
      mcnlx,
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
#line 18158 "./HZB_FLEX.c"
}   /* End of Origin=Progress_bar() SETTING parameter declarations. */
#undef CurrentTime
#undef EndTime
#undef StartTime
#undef IntermediateCnts
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'Before_selec'. */
  SIG_MESSAGE("Before_selec (Save)");
#define mccompcurname  Before_selec
#define mccompcurtype  PSD_monitor
#define mccompcurindex 10
#define PSD_N mccBefore_selec_PSD_N
#define PSD_p mccBefore_selec_PSD_p
#define PSD_p2 mccBefore_selec_PSD_p2
{   /* Declarations of Before_selec=PSD_monitor() SETTING parameters. */
int nx = mccBefore_selec_nx;
int ny = mccBefore_selec_ny;
char* filename = mccBefore_selec_filename;
MCNUM xmin = mccBefore_selec_xmin;
MCNUM xmax = mccBefore_selec_xmax;
MCNUM ymin = mccBefore_selec_ymin;
MCNUM ymax = mccBefore_selec_ymax;
MCNUM xwidth = mccBefore_selec_xwidth;
MCNUM yheight = mccBefore_selec_yheight;
MCNUM restore_neutron = mccBefore_selec_restore_neutron;
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
#line 18198 "./HZB_FLEX.c"
}   /* End of Before_selec=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'After_selec'. */
  SIG_MESSAGE("After_selec (Save)");
#define mccompcurname  After_selec
#define mccompcurtype  PSD_monitor
#define mccompcurindex 12
#define PSD_N mccAfter_selec_PSD_N
#define PSD_p mccAfter_selec_PSD_p
#define PSD_p2 mccAfter_selec_PSD_p2
{   /* Declarations of After_selec=PSD_monitor() SETTING parameters. */
int nx = mccAfter_selec_nx;
int ny = mccAfter_selec_ny;
char* filename = mccAfter_selec_filename;
MCNUM xmin = mccAfter_selec_xmin;
MCNUM xmax = mccAfter_selec_xmax;
MCNUM ymin = mccAfter_selec_ymin;
MCNUM ymax = mccAfter_selec_ymax;
MCNUM xwidth = mccAfter_selec_xwidth;
MCNUM yheight = mccAfter_selec_yheight;
MCNUM restore_neutron = mccAfter_selec_restore_neutron;
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
#line 18237 "./HZB_FLEX.c"
}   /* End of After_selec=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'energy_endguide'. */
  SIG_MESSAGE("energy_endguide (Save)");
#define mccompcurname  energy_endguide
#define mccompcurtype  E_monitor
#define mccompcurindex 23
#define nE mccenergy_endguide_nE
#define E_N mccenergy_endguide_E_N
#define E_p mccenergy_endguide_E_p
#define E_p2 mccenergy_endguide_E_p2
#define S_p mccenergy_endguide_S_p
#define S_pE mccenergy_endguide_S_pE
#define S_pE2 mccenergy_endguide_S_pE2
{   /* Declarations of energy_endguide=E_monitor() SETTING parameters. */
char* filename = mccenergy_endguide_filename;
MCNUM xmin = mccenergy_endguide_xmin;
MCNUM xmax = mccenergy_endguide_xmax;
MCNUM ymin = mccenergy_endguide_ymin;
MCNUM ymax = mccenergy_endguide_ymax;
MCNUM xwidth = mccenergy_endguide_xwidth;
MCNUM yheight = mccenergy_endguide_yheight;
MCNUM Emin = mccenergy_endguide_Emin;
MCNUM Emax = mccenergy_endguide_Emax;
MCNUM restore_neutron = mccenergy_endguide_restore_neutron;
int nowritefile = mccenergy_endguide_nowritefile;
#line 117 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/E_monitor.comp"
{
    if (!nowritefile) {
    DETECTOR_OUT_1D(
        "Energy monitor",
        "Energy [meV]",
        "Intensity",
        "E", Emin, Emax, nE,
        &E_N[0],&E_p[0],&E_p2[0],
        filename);
    if (S_p) printf("<E> : %g meV , E-width : %g meV \n",
     S_pE/S_p,sqrt(S_pE2/S_p - S_pE*S_pE/(S_p*S_p)) );
    }
}
#line 18284 "./HZB_FLEX.c"
}   /* End of energy_endguide=E_monitor() SETTING parameter declarations. */
#undef S_pE2
#undef S_pE
#undef S_p
#undef E_p2
#undef E_p
#undef E_N
#undef nE
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'energy_mono'. */
  SIG_MESSAGE("energy_mono (Save)");
#define mccompcurname  energy_mono
#define mccompcurtype  E_monitor
#define mccompcurindex 27
#define nE mccenergy_mono_nE
#define E_N mccenergy_mono_E_N
#define E_p mccenergy_mono_E_p
#define E_p2 mccenergy_mono_E_p2
#define S_p mccenergy_mono_S_p
#define S_pE mccenergy_mono_S_pE
#define S_pE2 mccenergy_mono_S_pE2
{   /* Declarations of energy_mono=E_monitor() SETTING parameters. */
char* filename = mccenergy_mono_filename;
MCNUM xmin = mccenergy_mono_xmin;
MCNUM xmax = mccenergy_mono_xmax;
MCNUM ymin = mccenergy_mono_ymin;
MCNUM ymax = mccenergy_mono_ymax;
MCNUM xwidth = mccenergy_mono_xwidth;
MCNUM yheight = mccenergy_mono_yheight;
MCNUM Emin = mccenergy_mono_Emin;
MCNUM Emax = mccenergy_mono_Emax;
MCNUM restore_neutron = mccenergy_mono_restore_neutron;
int nowritefile = mccenergy_mono_nowritefile;
#line 117 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/E_monitor.comp"
{
    if (!nowritefile) {
    DETECTOR_OUT_1D(
        "Energy monitor",
        "Energy [meV]",
        "Intensity",
        "E", Emin, Emax, nE,
        &E_N[0],&E_p[0],&E_p2[0],
        filename);
    if (S_p) printf("<E> : %g meV , E-width : %g meV \n",
     S_pE/S_p,sqrt(S_pE2/S_p - S_pE*S_pE/(S_p*S_p)) );
    }
}
#line 18335 "./HZB_FLEX.c"
}   /* End of energy_mono=E_monitor() SETTING parameter declarations. */
#undef S_pE2
#undef S_pE
#undef S_p
#undef E_p2
#undef E_p
#undef E_N
#undef nE
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'energy_pre_sample'. */
  SIG_MESSAGE("energy_pre_sample (Save)");
#define mccompcurname  energy_pre_sample
#define mccompcurtype  E_monitor
#define mccompcurindex 28
#define nE mccenergy_pre_sample_nE
#define E_N mccenergy_pre_sample_E_N
#define E_p mccenergy_pre_sample_E_p
#define E_p2 mccenergy_pre_sample_E_p2
#define S_p mccenergy_pre_sample_S_p
#define S_pE mccenergy_pre_sample_S_pE
#define S_pE2 mccenergy_pre_sample_S_pE2
{   /* Declarations of energy_pre_sample=E_monitor() SETTING parameters. */
char* filename = mccenergy_pre_sample_filename;
MCNUM xmin = mccenergy_pre_sample_xmin;
MCNUM xmax = mccenergy_pre_sample_xmax;
MCNUM ymin = mccenergy_pre_sample_ymin;
MCNUM ymax = mccenergy_pre_sample_ymax;
MCNUM xwidth = mccenergy_pre_sample_xwidth;
MCNUM yheight = mccenergy_pre_sample_yheight;
MCNUM Emin = mccenergy_pre_sample_Emin;
MCNUM Emax = mccenergy_pre_sample_Emax;
MCNUM restore_neutron = mccenergy_pre_sample_restore_neutron;
int nowritefile = mccenergy_pre_sample_nowritefile;
#line 117 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/E_monitor.comp"
{
    if (!nowritefile) {
    DETECTOR_OUT_1D(
        "Energy monitor",
        "Energy [meV]",
        "Intensity",
        "E", Emin, Emax, nE,
        &E_N[0],&E_p[0],&E_p2[0],
        filename);
    if (S_p) printf("<E> : %g meV , E-width : %g meV \n",
     S_pE/S_p,sqrt(S_pE2/S_p - S_pE*S_pE/(S_p*S_p)) );
    }
}
#line 18386 "./HZB_FLEX.c"
}   /* End of energy_pre_sample=E_monitor() SETTING parameter declarations. */
#undef S_pE2
#undef S_pE
#undef S_p
#undef E_p2
#undef E_p
#undef E_N
#undef nE
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'div_mono'. */
  SIG_MESSAGE("div_mono (Save)");
#define mccompcurname  div_mono
#define mccompcurtype  Divergence_monitor
#define mccompcurindex 31
#define nh mccdiv_mono_nh
#define nv mccdiv_mono_nv
#define Div_N mccdiv_mono_Div_N
#define Div_p mccdiv_mono_Div_p
#define Div_p2 mccdiv_mono_Div_p2
{   /* Declarations of div_mono=Divergence_monitor() SETTING parameters. */
char* filename = mccdiv_mono_filename;
MCNUM xmin = mccdiv_mono_xmin;
MCNUM xmax = mccdiv_mono_xmax;
MCNUM ymin = mccdiv_mono_ymin;
MCNUM ymax = mccdiv_mono_ymax;
MCNUM xwidth = mccdiv_mono_xwidth;
MCNUM yheight = mccdiv_mono_yheight;
MCNUM maxdiv_h = mccdiv_mono_maxdiv_h;
MCNUM maxdiv_v = mccdiv_mono_maxdiv_v;
MCNUM restore_neutron = mccdiv_mono_restore_neutron;
MCNUM nx = mccdiv_mono_nx;
MCNUM ny = mccdiv_mono_ny;
MCNUM nz = mccdiv_mono_nz;
int nowritefile = mccdiv_mono_nowritefile;
#line 117 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/Divergence_monitor.comp"
{
    if (!nowritefile) {
    DETECTOR_OUT_2D(
        "Divergence monitor",
        "X divergence [deg]",
        "Y divergence [deg]",
        -maxdiv_h, maxdiv_h, -maxdiv_v, maxdiv_v,
        nh, nv,
        &Div_N[0][0],&Div_p[0][0],&Div_p2[0][0],
        filename);
    }
}
#line 18437 "./HZB_FLEX.c"
}   /* End of div_mono=Divergence_monitor() SETTING parameter declarations. */
#undef Div_p2
#undef Div_p
#undef Div_N
#undef nv
#undef nh
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'div_mono_H'. */
  SIG_MESSAGE("div_mono_H (Save)");
#define mccompcurname  div_mono_H
#define mccompcurtype  Hdiv_monitor
#define mccompcurindex 32
#define nh mccdiv_mono_H_nh
#define Div_N mccdiv_mono_H_Div_N
#define Div_p mccdiv_mono_H_Div_p
#define Div_p2 mccdiv_mono_H_Div_p2
{   /* Declarations of div_mono_H=Hdiv_monitor() SETTING parameters. */
char* filename = mccdiv_mono_H_filename;
MCNUM xmin = mccdiv_mono_H_xmin;
MCNUM xmax = mccdiv_mono_H_xmax;
MCNUM ymin = mccdiv_mono_H_ymin;
MCNUM ymax = mccdiv_mono_H_ymax;
MCNUM xwidth = mccdiv_mono_H_xwidth;
MCNUM yheight = mccdiv_mono_H_yheight;
MCNUM h_maxdiv = mccdiv_mono_H_h_maxdiv;
MCNUM restore_neutron = mccdiv_mono_H_restore_neutron;
int nowritefile = mccdiv_mono_H_nowritefile;
#line 107 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/Hdiv_monitor.comp"
{
    if (!nowritefile) {
    DETECTOR_OUT_1D(
        "horizontal divergence monitor",
        "horizontal divergence [deg]",
        "Intensity",
        "divergence", -h_maxdiv, h_maxdiv, nh,
        &Div_N[0],&Div_p[0],&Div_p2[0],
        filename);
    }
}
#line 18480 "./HZB_FLEX.c"
}   /* End of div_mono_H=Hdiv_monitor() SETTING parameter declarations. */
#undef Div_p2
#undef Div_p
#undef Div_N
#undef nh
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'psd_sam'. */
  SIG_MESSAGE("psd_sam (Save)");
#define mccompcurname  psd_sam
#define mccompcurtype  PSD_monitor
#define mccompcurindex 33
#define PSD_N mccpsd_sam_PSD_N
#define PSD_p mccpsd_sam_PSD_p
#define PSD_p2 mccpsd_sam_PSD_p2
{   /* Declarations of psd_sam=PSD_monitor() SETTING parameters. */
int nx = mccpsd_sam_nx;
int ny = mccpsd_sam_ny;
char* filename = mccpsd_sam_filename;
MCNUM xmin = mccpsd_sam_xmin;
MCNUM xmax = mccpsd_sam_xmax;
MCNUM ymin = mccpsd_sam_ymin;
MCNUM ymax = mccpsd_sam_ymax;
MCNUM xwidth = mccpsd_sam_xwidth;
MCNUM yheight = mccpsd_sam_yheight;
MCNUM restore_neutron = mccpsd_sam_restore_neutron;
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
#line 18520 "./HZB_FLEX.c"
}   /* End of psd_sam=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component '1dpsd'. */
  SIG_MESSAGE("1dpsd (Save)");
#define mccompcurname  1dpsd
#define mccompcurtype  PSDlin_monitor
#define mccompcurindex 34
#define nx mcc1dpsd_nx
#define PSDlin_N mcc1dpsd_PSDlin_N
#define PSDlin_p mcc1dpsd_PSDlin_p
#define PSDlin_p2 mcc1dpsd_PSDlin_p2
{   /* Declarations of 1dpsd=PSDlin_monitor() SETTING parameters. */
char* filename = mcc1dpsd_filename;
MCNUM xmin = mcc1dpsd_xmin;
MCNUM xmax = mcc1dpsd_xmax;
MCNUM ymin = mcc1dpsd_ymin;
MCNUM ymax = mcc1dpsd_ymax;
MCNUM xwidth = mcc1dpsd_xwidth;
MCNUM yheight = mcc1dpsd_yheight;
MCNUM restore_neutron = mcc1dpsd_restore_neutron;
int nowritefile = mcc1dpsd_nowritefile;
#line 103 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/PSDlin_monitor.comp"
{
    if (!nowritefile) {
    DETECTOR_OUT_1D(
        "Linear PSD monitor",
        "x-Position [m]",
        "Intensity",
        "x", xmin, xmax, nx,
        &PSDlin_N[0],&PSDlin_p[0],&PSDlin_p2[0],
        filename);
    }
}
#line 18560 "./HZB_FLEX.c"
}   /* End of 1dpsd=PSDlin_monitor() SETTING parameter declarations. */
#undef PSDlin_p2
#undef PSDlin_p
#undef PSDlin_N
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
#line 18604 "./HZB_FLEX.c"
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
  /* User FINALLY code for component 'Gen_Source'. */
  SIG_MESSAGE("Gen_Source (Finally)");
#define mccompcurname  Gen_Source
#define mccompcurtype  Source_gen
#define mccompcurindex 2
#define p_in mccGen_Source_p_in
#define lambda1 mccGen_Source_lambda1
#define lambda2 mccGen_Source_lambda2
#define lambda3 mccGen_Source_lambda3
#define pTable mccGen_Source_pTable
#define pTable_x mccGen_Source_pTable_x
#define pTable_y mccGen_Source_pTable_y
#define pTable_xmin mccGen_Source_pTable_xmin
#define pTable_xmax mccGen_Source_pTable_xmax
#define pTable_xsum mccGen_Source_pTable_xsum
#define pTable_ymin mccGen_Source_pTable_ymin
#define pTable_ymax mccGen_Source_pTable_ymax
#define pTable_ysum mccGen_Source_pTable_ysum
#define pTable_dxmin mccGen_Source_pTable_dxmin
#define pTable_dxmax mccGen_Source_pTable_dxmax
#define pTable_dymin mccGen_Source_pTable_dymin
#define pTable_dymax mccGen_Source_pTable_dymax
{   /* Declarations of Gen_Source=Source_gen() SETTING parameters. */
char* flux_file = mccGen_Source_flux_file;
char* xdiv_file = mccGen_Source_xdiv_file;
char* ydiv_file = mccGen_Source_ydiv_file;
MCNUM radius = mccGen_Source_radius;
MCNUM dist = mccGen_Source_dist;
MCNUM focus_xw = mccGen_Source_focus_xw;
MCNUM focus_yh = mccGen_Source_focus_yh;
MCNUM focus_aw = mccGen_Source_focus_aw;
MCNUM focus_ah = mccGen_Source_focus_ah;
MCNUM E0 = mccGen_Source_E0;
MCNUM dE = mccGen_Source_dE;
MCNUM lambda0 = mccGen_Source_lambda0;
MCNUM dlambda = mccGen_Source_dlambda;
MCNUM I1 = mccGen_Source_I1;
MCNUM yheight = mccGen_Source_yheight;
MCNUM xwidth = mccGen_Source_xwidth;
MCNUM verbose = mccGen_Source_verbose;
MCNUM T1 = mccGen_Source_T1;
MCNUM flux_file_perAA = mccGen_Source_flux_file_perAA;
MCNUM flux_file_log = mccGen_Source_flux_file_log;
MCNUM Lmin = mccGen_Source_Lmin;
MCNUM Lmax = mccGen_Source_Lmax;
MCNUM Emin = mccGen_Source_Emin;
MCNUM Emax = mccGen_Source_Emax;
MCNUM T2 = mccGen_Source_T2;
MCNUM I2 = mccGen_Source_I2;
MCNUM T3 = mccGen_Source_T3;
MCNUM I3 = mccGen_Source_I3;
MCNUM zdepth = mccGen_Source_zdepth;
int target_index = mccGen_Source_target_index;
#line 571 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../sources/Source_gen.comp"
{
  Table_Free(&pTable);
  Table_Free(&pTable_x);
  Table_Free(&pTable_y);
}
#line 18675 "./HZB_FLEX.c"
}   /* End of Gen_Source=Source_gen() SETTING parameter declarations. */
#undef pTable_dymax
#undef pTable_dymin
#undef pTable_dxmax
#undef pTable_dxmin
#undef pTable_ysum
#undef pTable_ymax
#undef pTable_ymin
#undef pTable_xsum
#undef pTable_xmax
#undef pTable_xmin
#undef pTable_y
#undef pTable_x
#undef pTable
#undef lambda3
#undef lambda2
#undef lambda1
#undef p_in
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[2]) fprintf(stderr, "Warning: No neutron could reach Component[2] Gen_Source\n");
    if (mcAbsorbProp[2]) fprintf(stderr, "Warning: %g events were removed in Component[2] Gen_Source=Source_gen()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[2]);
    if (!mcNCounter[3]) fprintf(stderr, "Warning: No neutron could reach Component[3] NL1A_NL1B_NL2_NL3_InPile_Entrance_Window\n");
    if (mcAbsorbProp[3]) fprintf(stderr, "Warning: %g events were removed in Component[3] NL1A_NL1B_NL2_NL3_InPile_Entrance_Window=Arm()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[3]);
    if (!mcNCounter[4]) fprintf(stderr, "Warning: No neutron could reach Component[4] NL1A_NL1B_NL2_NL3_InPile\n");
    if (mcAbsorbProp[4]) fprintf(stderr, "Warning: %g events were removed in Component[4] NL1A_NL1B_NL2_NL3_InPile=Guide()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[4]);
    if (!mcNCounter[5]) fprintf(stderr, "Warning: No neutron could reach Component[5] NL1B_Straight1_Entrance_Window\n");
    if (mcAbsorbProp[5]) fprintf(stderr, "Warning: %g events were removed in Component[5] NL1B_Straight1_Entrance_Window=Arm()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[5]);
    if (!mcNCounter[6]) fprintf(stderr, "Warning: No neutron could reach Component[6] NL1B_Straight1\n");
    if (mcAbsorbProp[6]) fprintf(stderr, "Warning: %g events were removed in Component[6] NL1B_Straight1=Guide()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[6]);
    if (!mcNCounter[7]) fprintf(stderr, "Warning: No neutron could reach Component[7] NL1B_Curved_Entrance_Window\n");
    if (mcAbsorbProp[7]) fprintf(stderr, "Warning: %g events were removed in Component[7] NL1B_Curved_Entrance_Window=Arm()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[7]);
    if (!mcNCounter[8]) fprintf(stderr, "Warning: No neutron could reach Component[8] NL1B_Curved_Guide_1\n");
    if (mcAbsorbProp[8]) fprintf(stderr, "Warning: %g events were removed in Component[8] NL1B_Curved_Guide_1=Guide_curved()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[8]);
    if (!mcNCounter[9]) fprintf(stderr, "Warning: No neutron could reach Component[9] NL1B_Velocity_Selector_Gap_Entrance_Window\n");
    if (mcAbsorbProp[9]) fprintf(stderr, "Warning: %g events were removed in Component[9] NL1B_Velocity_Selector_Gap_Entrance_Window=Arm()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[9]);
  /* User FINALLY code for component 'Before_selec'. */
  SIG_MESSAGE("Before_selec (Finally)");
#define mccompcurname  Before_selec
#define mccompcurtype  PSD_monitor
#define mccompcurindex 10
#define PSD_N mccBefore_selec_PSD_N
#define PSD_p mccBefore_selec_PSD_p
#define PSD_p2 mccBefore_selec_PSD_p2
{   /* Declarations of Before_selec=PSD_monitor() SETTING parameters. */
int nx = mccBefore_selec_nx;
int ny = mccBefore_selec_ny;
char* filename = mccBefore_selec_filename;
MCNUM xmin = mccBefore_selec_xmin;
MCNUM xmax = mccBefore_selec_xmax;
MCNUM ymin = mccBefore_selec_ymin;
MCNUM ymax = mccBefore_selec_ymax;
MCNUM xwidth = mccBefore_selec_xwidth;
MCNUM yheight = mccBefore_selec_yheight;
MCNUM restore_neutron = mccBefore_selec_restore_neutron;
#line 122 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
{
  destroy_darr2d(PSD_N);
  destroy_darr2d(PSD_p);
  destroy_darr2d(PSD_p2);
}
#line 18739 "./HZB_FLEX.c"
}   /* End of Before_selec=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[10]) fprintf(stderr, "Warning: No neutron could reach Component[10] Before_selec\n");
    if (mcAbsorbProp[10]) fprintf(stderr, "Warning: %g events were removed in Component[10] Before_selec=PSD_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[10]);
    if (!mcNCounter[11]) fprintf(stderr, "Warning: No neutron could reach Component[11] SELEC\n");
    if (mcAbsorbProp[11]) fprintf(stderr, "Warning: %g events were removed in Component[11] SELEC=Selector()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[11]);
  /* User FINALLY code for component 'After_selec'. */
  SIG_MESSAGE("After_selec (Finally)");
#define mccompcurname  After_selec
#define mccompcurtype  PSD_monitor
#define mccompcurindex 12
#define PSD_N mccAfter_selec_PSD_N
#define PSD_p mccAfter_selec_PSD_p
#define PSD_p2 mccAfter_selec_PSD_p2
{   /* Declarations of After_selec=PSD_monitor() SETTING parameters. */
int nx = mccAfter_selec_nx;
int ny = mccAfter_selec_ny;
char* filename = mccAfter_selec_filename;
MCNUM xmin = mccAfter_selec_xmin;
MCNUM xmax = mccAfter_selec_xmax;
MCNUM ymin = mccAfter_selec_ymin;
MCNUM ymax = mccAfter_selec_ymax;
MCNUM xwidth = mccAfter_selec_xwidth;
MCNUM yheight = mccAfter_selec_yheight;
MCNUM restore_neutron = mccAfter_selec_restore_neutron;
#line 122 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
{
  destroy_darr2d(PSD_N);
  destroy_darr2d(PSD_p);
  destroy_darr2d(PSD_p2);
}
#line 18777 "./HZB_FLEX.c"
}   /* End of After_selec=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[12]) fprintf(stderr, "Warning: No neutron could reach Component[12] After_selec\n");
    if (mcAbsorbProp[12]) fprintf(stderr, "Warning: %g events were removed in Component[12] After_selec=PSD_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[12]);
    if (!mcNCounter[13]) fprintf(stderr, "Warning: No neutron could reach Component[13] NL1B_Curved_2_Entrance_Window\n");
    if (mcAbsorbProp[13]) fprintf(stderr, "Warning: %g events were removed in Component[13] NL1B_Curved_2_Entrance_Window=Arm()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[13]);
    if (!mcNCounter[14]) fprintf(stderr, "Warning: No neutron could reach Component[14] NL1B_Curved_Guide_2\n");
    if (mcAbsorbProp[14]) fprintf(stderr, "Warning: %g events were removed in Component[14] NL1B_Curved_Guide_2=Guide_curved()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[14]);
    if (!mcNCounter[15]) fprintf(stderr, "Warning: No neutron could reach Component[15] NL1B_Straight2_Entrance_Window\n");
    if (mcAbsorbProp[15]) fprintf(stderr, "Warning: %g events were removed in Component[15] NL1B_Straight2_Entrance_Window=Arm()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[15]);
    if (!mcNCounter[16]) fprintf(stderr, "Warning: No neutron could reach Component[16] NL1B_Straight2\n");
    if (mcAbsorbProp[16]) fprintf(stderr, "Warning: %g events were removed in Component[16] NL1B_Straight2=Guide()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[16]);
    if (!mcNCounter[17]) fprintf(stderr, "Warning: No neutron could reach Component[17] NL1B_Elliptical_Entrance_Window\n");
    if (mcAbsorbProp[17]) fprintf(stderr, "Warning: %g events were removed in Component[17] NL1B_Elliptical_Entrance_Window=Arm()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[17]);
  /* User FINALLY code for component 'elliptical_piece'. */
  SIG_MESSAGE("elliptical_piece (Finally)");
#define mccompcurname  elliptical_piece
#define mccompcurtype  Guide_tapering
#define mccompcurindex 18
#define w1c mccelliptical_piece_w1c
#define w2c mccelliptical_piece_w2c
#define ww mccelliptical_piece_ww
#define hh mccelliptical_piece_hh
#define whalf mccelliptical_piece_whalf
#define hhalf mccelliptical_piece_hhalf
#define lwhalf mccelliptical_piece_lwhalf
#define lhhalf mccelliptical_piece_lhhalf
#define h1_in mccelliptical_piece_h1_in
#define h2_out mccelliptical_piece_h2_out
#define w1_in mccelliptical_piece_w1_in
#define w2_out mccelliptical_piece_w2_out
#define l_seg mccelliptical_piece_l_seg
#define seg mccelliptical_piece_seg
#define h12 mccelliptical_piece_h12
#define h2 mccelliptical_piece_h2
#define w12 mccelliptical_piece_w12
#define w2 mccelliptical_piece_w2
#define a_ell_q mccelliptical_piece_a_ell_q
#define b_ell_q mccelliptical_piece_b_ell_q
#define lbw mccelliptical_piece_lbw
#define lbh mccelliptical_piece_lbh
#define mxi mccelliptical_piece_mxi
#define u1 mccelliptical_piece_u1
#define u2 mccelliptical_piece_u2
#define div1 mccelliptical_piece_div1
#define p2_para mccelliptical_piece_p2_para
#define test mccelliptical_piece_test
#define Div1 mccelliptical_piece_Div1
#define i mccelliptical_piece_i
#define ii mccelliptical_piece_ii
#define seg mccelliptical_piece_seg
#define fu mccelliptical_piece_fu
#define pos mccelliptical_piece_pos
#define file_name mccelliptical_piece_file_name
#define ep mccelliptical_piece_ep
#define num mccelliptical_piece_num
#define rotation_h mccelliptical_piece_rotation_h
#define rotation_v mccelliptical_piece_rotation_v
{   /* Declarations of elliptical_piece=Guide_tapering() SETTING parameters. */
char* option = mccelliptical_piece_option;
MCNUM w1 = mccelliptical_piece_w1;
MCNUM h1 = mccelliptical_piece_h1;
MCNUM l = mccelliptical_piece_l;
MCNUM linw = mccelliptical_piece_linw;
MCNUM loutw = mccelliptical_piece_loutw;
MCNUM linh = mccelliptical_piece_linh;
MCNUM louth = mccelliptical_piece_louth;
MCNUM R0 = mccelliptical_piece_R0;
MCNUM Qcx = mccelliptical_piece_Qcx;
MCNUM Qcy = mccelliptical_piece_Qcy;
MCNUM alphax = mccelliptical_piece_alphax;
MCNUM alphay = mccelliptical_piece_alphay;
MCNUM W = mccelliptical_piece_W;
MCNUM mx = mccelliptical_piece_mx;
MCNUM my = mccelliptical_piece_my;
MCNUM segno = mccelliptical_piece_segno;
MCNUM curvature = mccelliptical_piece_curvature;
MCNUM curvature_v = mccelliptical_piece_curvature_v;
#line 609 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide_tapering.comp"
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
#line 18877 "./HZB_FLEX.c"
}   /* End of elliptical_piece=Guide_tapering() SETTING parameter declarations. */
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

    if (!mcNCounter[18]) fprintf(stderr, "Warning: No neutron could reach Component[18] elliptical_piece\n");
    if (mcAbsorbProp[18]) fprintf(stderr, "Warning: %g events were removed in Component[18] elliptical_piece=Guide_tapering()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[18]);
    if (!mcNCounter[19]) fprintf(stderr, "Warning: No neutron could reach Component[19] Virtual_source\n");
    if (mcAbsorbProp[19]) fprintf(stderr, "Warning: %g events were removed in Component[19] Virtual_source=Slit()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[19]);
    if (!mcNCounter[20]) fprintf(stderr, "Warning: No neutron could reach Component[20] NL1B_Vertical_Guide_Entrance_Window\n");
    if (mcAbsorbProp[20]) fprintf(stderr, "Warning: %g events were removed in Component[20] NL1B_Vertical_Guide_Entrance_Window=Arm()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[20]);
    if (!mcNCounter[21]) fprintf(stderr, "Warning: No neutron could reach Component[21] NL1B_Vertical_Guide\n");
    if (mcAbsorbProp[21]) fprintf(stderr, "Warning: %g events were removed in Component[21] NL1B_Vertical_Guide=Guide_channeled()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[21]);
    if (!mcNCounter[22]) fprintf(stderr, "Warning: No neutron could reach Component[22] NL1B_Guide_Exit\n");
    if (mcAbsorbProp[22]) fprintf(stderr, "Warning: %g events were removed in Component[22] NL1B_Guide_Exit=Arm()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[22]);
    if (!mcNCounter[23]) fprintf(stderr, "Warning: No neutron could reach Component[23] energy_endguide\n");
    if (mcAbsorbProp[23]) fprintf(stderr, "Warning: %g events were removed in Component[23] energy_endguide=E_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[23]);
    if (!mcNCounter[24]) fprintf(stderr, "Warning: No neutron could reach Component[24] Mono_center\n");
    if (mcAbsorbProp[24]) fprintf(stderr, "Warning: %g events were removed in Component[24] Mono_center=Arm()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[24]);
  /* User FINALLY code for component 'Monochromator'. */
  SIG_MESSAGE("Monochromator (Finally)");
#define mccompcurname  Monochromator
#define mccompcurtype  Monochromator_curved
#define mccompcurindex 25
#define mos_rms_y mccMonochromator_mos_rms_y
#define mos_rms_z mccMonochromator_mos_rms_z
#define mos_rms_max mccMonochromator_mos_rms_max
#define mono_Q mccMonochromator_mono_Q
#define SlabWidth mccMonochromator_SlabWidth
#define SlabHeight mccMonochromator_SlabHeight
#define rTable mccMonochromator_rTable
#define tTable mccMonochromator_tTable
#define row mccMonochromator_row
#define col mccMonochromator_col
#define tiltH mccMonochromator_tiltH
#define tiltV mccMonochromator_tiltV
{   /* Declarations of Monochromator=Monochromator_curved() SETTING parameters. */
char* reflect = mccMonochromator_reflect;
char* transmit = mccMonochromator_transmit;
MCNUM zwidth = mccMonochromator_zwidth;
MCNUM yheight = mccMonochromator_yheight;
MCNUM gap = mccMonochromator_gap;
MCNUM NH = mccMonochromator_NH;
MCNUM NV = mccMonochromator_NV;
MCNUM mosaich = mccMonochromator_mosaich;
MCNUM mosaicv = mccMonochromator_mosaicv;
MCNUM r0 = mccMonochromator_r0;
MCNUM t0 = mccMonochromator_t0;
MCNUM Q = mccMonochromator_Q;
MCNUM RV = mccMonochromator_RV;
MCNUM RH = mccMonochromator_RH;
MCNUM DM = mccMonochromator_DM;
MCNUM mosaic = mccMonochromator_mosaic;
MCNUM width = mccMonochromator_width;
MCNUM height = mccMonochromator_height;
MCNUM verbose = mccMonochromator_verbose;
MCNUM order = mccMonochromator_order;
#line 460 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Monochromator_curved.comp"
{
  Table_Free(&rTable);
  Table_Free(&tTable);
  if (tiltH) free(tiltH);
  if (tiltV) free(tiltV);
}
#line 18981 "./HZB_FLEX.c"
}   /* End of Monochromator=Monochromator_curved() SETTING parameter declarations. */
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

    if (!mcNCounter[25]) fprintf(stderr, "Warning: No neutron could reach Component[25] Monochromator\n");
    if (mcAbsorbProp[25]) fprintf(stderr, "Warning: %g events were removed in Component[25] Monochromator=Monochromator_curved()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[25]);
    if (!mcNCounter[26]) fprintf(stderr, "Warning: No neutron could reach Component[26] Mono_sample_arm\n");
    if (mcAbsorbProp[26]) fprintf(stderr, "Warning: %g events were removed in Component[26] Mono_sample_arm=Arm()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[26]);
    if (!mcNCounter[27]) fprintf(stderr, "Warning: No neutron could reach Component[27] energy_mono\n");
    if (mcAbsorbProp[27]) fprintf(stderr, "Warning: %g events were removed in Component[27] energy_mono=E_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[27]);
    if (!mcNCounter[28]) fprintf(stderr, "Warning: No neutron could reach Component[28] energy_pre_sample\n");
    if (mcAbsorbProp[28]) fprintf(stderr, "Warning: %g events were removed in Component[28] energy_pre_sample=E_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[28]);
    if (!mcNCounter[29]) fprintf(stderr, "Warning: No neutron could reach Component[29] Sample_center\n");
    if (mcAbsorbProp[29]) fprintf(stderr, "Warning: %g events were removed in Component[29] Sample_center=Arm()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[29]);
    if (!mcNCounter[30]) fprintf(stderr, "Warning: No neutron could reach Component[30] Sample_analyser_arm\n");
    if (mcAbsorbProp[30]) fprintf(stderr, "Warning: %g events were removed in Component[30] Sample_analyser_arm=Arm()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[30]);
    if (!mcNCounter[31]) fprintf(stderr, "Warning: No neutron could reach Component[31] div_mono\n");
    if (mcAbsorbProp[31]) fprintf(stderr, "Warning: %g events were removed in Component[31] div_mono=Divergence_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[31]);
    if (!mcNCounter[32]) fprintf(stderr, "Warning: No neutron could reach Component[32] div_mono_H\n");
    if (mcAbsorbProp[32]) fprintf(stderr, "Warning: %g events were removed in Component[32] div_mono_H=Hdiv_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[32]);
  /* User FINALLY code for component 'psd_sam'. */
  SIG_MESSAGE("psd_sam (Finally)");
#define mccompcurname  psd_sam
#define mccompcurtype  PSD_monitor
#define mccompcurindex 33
#define PSD_N mccpsd_sam_PSD_N
#define PSD_p mccpsd_sam_PSD_p
#define PSD_p2 mccpsd_sam_PSD_p2
{   /* Declarations of psd_sam=PSD_monitor() SETTING parameters. */
int nx = mccpsd_sam_nx;
int ny = mccpsd_sam_ny;
char* filename = mccpsd_sam_filename;
MCNUM xmin = mccpsd_sam_xmin;
MCNUM xmax = mccpsd_sam_xmax;
MCNUM ymin = mccpsd_sam_ymin;
MCNUM ymax = mccpsd_sam_ymax;
MCNUM xwidth = mccpsd_sam_xwidth;
MCNUM yheight = mccpsd_sam_yheight;
MCNUM restore_neutron = mccpsd_sam_restore_neutron;
#line 122 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
{
  destroy_darr2d(PSD_N);
  destroy_darr2d(PSD_p);
  destroy_darr2d(PSD_p2);
}
#line 19040 "./HZB_FLEX.c"
}   /* End of psd_sam=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[33]) fprintf(stderr, "Warning: No neutron could reach Component[33] psd_sam\n");
    if (mcAbsorbProp[33]) fprintf(stderr, "Warning: %g events were removed in Component[33] psd_sam=PSD_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[33]);
    if (!mcNCounter[34]) fprintf(stderr, "Warning: No neutron could reach Component[34] 1dpsd\n");
    if (mcAbsorbProp[34]) fprintf(stderr, "Warning: %g events were removed in Component[34] 1dpsd=PSDlin_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[34]);
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
#line 19087 "./HZB_FLEX.c"
}   /* End of Origin=Progress_bar() SETTING parameter declarations. */
#undef CurrentTime
#undef EndTime
#undef StartTime
#undef IntermediateCnts
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'Gen_Source'. */
  SIG_MESSAGE("Gen_Source (McDisplay)");
  printf("MCDISPLAY: component %s\n", "Gen_Source");
#define mccompcurname  Gen_Source
#define mccompcurtype  Source_gen
#define mccompcurindex 2
#define p_in mccGen_Source_p_in
#define lambda1 mccGen_Source_lambda1
#define lambda2 mccGen_Source_lambda2
#define lambda3 mccGen_Source_lambda3
#define pTable mccGen_Source_pTable
#define pTable_x mccGen_Source_pTable_x
#define pTable_y mccGen_Source_pTable_y
#define pTable_xmin mccGen_Source_pTable_xmin
#define pTable_xmax mccGen_Source_pTable_xmax
#define pTable_xsum mccGen_Source_pTable_xsum
#define pTable_ymin mccGen_Source_pTable_ymin
#define pTable_ymax mccGen_Source_pTable_ymax
#define pTable_ysum mccGen_Source_pTable_ysum
#define pTable_dxmin mccGen_Source_pTable_dxmin
#define pTable_dxmax mccGen_Source_pTable_dxmax
#define pTable_dymin mccGen_Source_pTable_dymin
#define pTable_dymax mccGen_Source_pTable_dymax
{   /* Declarations of Gen_Source=Source_gen() SETTING parameters. */
char* flux_file = mccGen_Source_flux_file;
char* xdiv_file = mccGen_Source_xdiv_file;
char* ydiv_file = mccGen_Source_ydiv_file;
MCNUM radius = mccGen_Source_radius;
MCNUM dist = mccGen_Source_dist;
MCNUM focus_xw = mccGen_Source_focus_xw;
MCNUM focus_yh = mccGen_Source_focus_yh;
MCNUM focus_aw = mccGen_Source_focus_aw;
MCNUM focus_ah = mccGen_Source_focus_ah;
MCNUM E0 = mccGen_Source_E0;
MCNUM dE = mccGen_Source_dE;
MCNUM lambda0 = mccGen_Source_lambda0;
MCNUM dlambda = mccGen_Source_dlambda;
MCNUM I1 = mccGen_Source_I1;
MCNUM yheight = mccGen_Source_yheight;
MCNUM xwidth = mccGen_Source_xwidth;
MCNUM verbose = mccGen_Source_verbose;
MCNUM T1 = mccGen_Source_T1;
MCNUM flux_file_perAA = mccGen_Source_flux_file_perAA;
MCNUM flux_file_log = mccGen_Source_flux_file_log;
MCNUM Lmin = mccGen_Source_Lmin;
MCNUM Lmax = mccGen_Source_Lmax;
MCNUM Emin = mccGen_Source_Emin;
MCNUM Emax = mccGen_Source_Emax;
MCNUM T2 = mccGen_Source_T2;
MCNUM I2 = mccGen_Source_I2;
MCNUM T3 = mccGen_Source_T3;
MCNUM I3 = mccGen_Source_I3;
MCNUM zdepth = mccGen_Source_zdepth;
int target_index = mccGen_Source_target_index;
#line 578 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../sources/Source_gen.comp"
{
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
}
#line 19200 "./HZB_FLEX.c"
}   /* End of Gen_Source=Source_gen() SETTING parameter declarations. */
#undef pTable_dymax
#undef pTable_dymin
#undef pTable_dxmax
#undef pTable_dxmin
#undef pTable_ysum
#undef pTable_ymax
#undef pTable_ymin
#undef pTable_xsum
#undef pTable_xmax
#undef pTable_xmin
#undef pTable_y
#undef pTable_x
#undef pTable
#undef lambda3
#undef lambda2
#undef lambda1
#undef p_in
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'NL1A_NL1B_NL2_NL3_InPile_Entrance_Window'. */
  SIG_MESSAGE("NL1A_NL1B_NL2_NL3_InPile_Entrance_Window (McDisplay)");
  printf("MCDISPLAY: component %s\n", "NL1A_NL1B_NL2_NL3_InPile_Entrance_Window");
#define mccompcurname  NL1A_NL1B_NL2_NL3_InPile_Entrance_Window
#define mccompcurtype  Arm
#define mccompcurindex 3
#line 40 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Arm.comp"
{
  /* A bit ugly; hard-coded dimensions. */
  
  line(0,0,0,0.2,0,0);
  line(0,0,0,0,0.2,0);
  line(0,0,0,0,0,0.2);
}
#line 19237 "./HZB_FLEX.c"
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'NL1A_NL1B_NL2_NL3_InPile'. */
  SIG_MESSAGE("NL1A_NL1B_NL2_NL3_InPile (McDisplay)");
  printf("MCDISPLAY: component %s\n", "NL1A_NL1B_NL2_NL3_InPile");
#define mccompcurname  NL1A_NL1B_NL2_NL3_InPile
#define mccompcurtype  Guide
#define mccompcurindex 4
#define pTable mccNL1A_NL1B_NL2_NL3_InPile_pTable
{   /* Declarations of NL1A_NL1B_NL2_NL3_InPile=Guide() SETTING parameters. */
char* reflect = mccNL1A_NL1B_NL2_NL3_InPile_reflect;
MCNUM w1 = mccNL1A_NL1B_NL2_NL3_InPile_w1;
MCNUM h1 = mccNL1A_NL1B_NL2_NL3_InPile_h1;
MCNUM w2 = mccNL1A_NL1B_NL2_NL3_InPile_w2;
MCNUM h2 = mccNL1A_NL1B_NL2_NL3_InPile_h2;
MCNUM l = mccNL1A_NL1B_NL2_NL3_InPile_l;
MCNUM R0 = mccNL1A_NL1B_NL2_NL3_InPile_R0;
MCNUM Qc = mccNL1A_NL1B_NL2_NL3_InPile_Qc;
MCNUM alpha = mccNL1A_NL1B_NL2_NL3_InPile_alpha;
MCNUM m = mccNL1A_NL1B_NL2_NL3_InPile_m;
MCNUM W = mccNL1A_NL1B_NL2_NL3_InPile_W;
#line 202 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide.comp"
{
  
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
}
#line 19281 "./HZB_FLEX.c"
}   /* End of NL1A_NL1B_NL2_NL3_InPile=Guide() SETTING parameter declarations. */
#undef pTable
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'NL1B_Straight1_Entrance_Window'. */
  SIG_MESSAGE("NL1B_Straight1_Entrance_Window (McDisplay)");
  printf("MCDISPLAY: component %s\n", "NL1B_Straight1_Entrance_Window");
#define mccompcurname  NL1B_Straight1_Entrance_Window
#define mccompcurtype  Arm
#define mccompcurindex 5
#line 40 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Arm.comp"
{
  /* A bit ugly; hard-coded dimensions. */
  
  line(0,0,0,0.2,0,0);
  line(0,0,0,0,0.2,0);
  line(0,0,0,0,0,0.2);
}
#line 19302 "./HZB_FLEX.c"
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'NL1B_Straight1'. */
  SIG_MESSAGE("NL1B_Straight1 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "NL1B_Straight1");
#define mccompcurname  NL1B_Straight1
#define mccompcurtype  Guide
#define mccompcurindex 6
#define pTable mccNL1B_Straight1_pTable
{   /* Declarations of NL1B_Straight1=Guide() SETTING parameters. */
char* reflect = mccNL1B_Straight1_reflect;
MCNUM w1 = mccNL1B_Straight1_w1;
MCNUM h1 = mccNL1B_Straight1_h1;
MCNUM w2 = mccNL1B_Straight1_w2;
MCNUM h2 = mccNL1B_Straight1_h2;
MCNUM l = mccNL1B_Straight1_l;
MCNUM R0 = mccNL1B_Straight1_R0;
MCNUM Qc = mccNL1B_Straight1_Qc;
MCNUM alpha = mccNL1B_Straight1_alpha;
MCNUM m = mccNL1B_Straight1_m;
MCNUM W = mccNL1B_Straight1_W;
#line 202 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide.comp"
{
  
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
}
#line 19346 "./HZB_FLEX.c"
}   /* End of NL1B_Straight1=Guide() SETTING parameter declarations. */
#undef pTable
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'NL1B_Curved_Entrance_Window'. */
  SIG_MESSAGE("NL1B_Curved_Entrance_Window (McDisplay)");
  printf("MCDISPLAY: component %s\n", "NL1B_Curved_Entrance_Window");
#define mccompcurname  NL1B_Curved_Entrance_Window
#define mccompcurtype  Arm
#define mccompcurindex 7
#line 40 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Arm.comp"
{
  /* A bit ugly; hard-coded dimensions. */
  
  line(0,0,0,0.2,0,0);
  line(0,0,0,0,0.2,0);
  line(0,0,0,0,0,0.2);
}
#line 19367 "./HZB_FLEX.c"
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'NL1B_Curved_Guide_1'. */
  SIG_MESSAGE("NL1B_Curved_Guide_1 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "NL1B_Curved_Guide_1");
#define mccompcurname  NL1B_Curved_Guide_1
#define mccompcurtype  Guide_curved
#define mccompcurindex 8
{   /* Declarations of NL1B_Curved_Guide_1=Guide_curved() SETTING parameters. */
MCNUM w1 = mccNL1B_Curved_Guide_1_w1;
MCNUM h1 = mccNL1B_Curved_Guide_1_h1;
MCNUM l = mccNL1B_Curved_Guide_1_l;
MCNUM R0 = mccNL1B_Curved_Guide_1_R0;
MCNUM Qc = mccNL1B_Curved_Guide_1_Qc;
MCNUM alpha = mccNL1B_Curved_Guide_1_alpha;
MCNUM m = mccNL1B_Curved_Guide_1_m;
MCNUM W = mccNL1B_Curved_Guide_1_W;
MCNUM curvature = mccNL1B_Curved_Guide_1_curvature;
#line 133 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../contrib/Guide_curved.comp"
{
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
}
#line 19424 "./HZB_FLEX.c"
}   /* End of NL1B_Curved_Guide_1=Guide_curved() SETTING parameter declarations. */
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'NL1B_Velocity_Selector_Gap_Entrance_Window'. */
  SIG_MESSAGE("NL1B_Velocity_Selector_Gap_Entrance_Window (McDisplay)");
  printf("MCDISPLAY: component %s\n", "NL1B_Velocity_Selector_Gap_Entrance_Window");
#define mccompcurname  NL1B_Velocity_Selector_Gap_Entrance_Window
#define mccompcurtype  Arm
#define mccompcurindex 9
#line 40 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Arm.comp"
{
  /* A bit ugly; hard-coded dimensions. */
  
  line(0,0,0,0.2,0,0);
  line(0,0,0,0,0.2,0);
  line(0,0,0,0,0,0.2);
}
#line 19444 "./HZB_FLEX.c"
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'Before_selec'. */
  SIG_MESSAGE("Before_selec (McDisplay)");
  printf("MCDISPLAY: component %s\n", "Before_selec");
#define mccompcurname  Before_selec
#define mccompcurtype  PSD_monitor
#define mccompcurindex 10
#define PSD_N mccBefore_selec_PSD_N
#define PSD_p mccBefore_selec_PSD_p
#define PSD_p2 mccBefore_selec_PSD_p2
{   /* Declarations of Before_selec=PSD_monitor() SETTING parameters. */
int nx = mccBefore_selec_nx;
int ny = mccBefore_selec_ny;
char* filename = mccBefore_selec_filename;
MCNUM xmin = mccBefore_selec_xmin;
MCNUM xmax = mccBefore_selec_xmax;
MCNUM ymin = mccBefore_selec_ymin;
MCNUM ymax = mccBefore_selec_ymax;
MCNUM xwidth = mccBefore_selec_xwidth;
MCNUM yheight = mccBefore_selec_yheight;
MCNUM restore_neutron = mccBefore_selec_restore_neutron;
#line 129 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
{
  multiline(5,
    (double)xmin, (double)ymin, 0.0,
    (double)xmax, (double)ymin, 0.0,
    (double)xmax, (double)ymax, 0.0,
    (double)xmin, (double)ymax, 0.0,
    (double)xmin, (double)ymin, 0.0);
}
#line 19478 "./HZB_FLEX.c"
}   /* End of Before_selec=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'SELEC'. */
  SIG_MESSAGE("SELEC (McDisplay)");
  printf("MCDISPLAY: component %s\n", "SELEC");
#define mccompcurname  SELEC
#define mccompcurtype  Selector
#define mccompcurindex 11
{   /* Declarations of SELEC=Selector() SETTING parameters. */
MCNUM xmin = mccSELEC_xmin;
MCNUM xmax = mccSELEC_xmax;
MCNUM ymin = mccSELEC_ymin;
MCNUM ymax = mccSELEC_ymax;
MCNUM length = mccSELEC_length;
MCNUM xwidth = mccSELEC_xwidth;
MCNUM yheight = mccSELEC_yheight;
MCNUM nslit = mccSELEC_nslit;
MCNUM d = mccSELEC_d;
MCNUM radius = mccSELEC_radius;
MCNUM alpha = mccSELEC_alpha;
MCNUM nu = mccSELEC_nu;
#line 114 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Selector.comp"
{
  double phi, r0, Width, height, l0, l1;
  double r;
  double x0;
  double x1;
  double y0;
  double y1;
  double z0;
  double z1;
  double z2;
  double z3;
  double a;
  double xw, yh;

  phi = alpha;
  Width  = (xmax-xmin)/2;
  height = (ymax-ymin)/2;
  x0 = xmin; x1 = xmax;
  y0 = ymin; y1 = ymax;
  l0 = length; l1 = l0;
  r0 = radius;

  r = r0 + height;
  x0 = -Width/2.0;
  x1 =  Width/2.0;
  y0 = -height/2.0;
  y1 =  height/2.0;
  z0 =  0;
  z1 =  0;
  z2 =  l1;
  z3 =  l0;

  
  xw = Width/2.0;
  yh = height/2.0;
  /* Draw apertures */
  for(a = z0;;)
  {
    multiline(3, x0-xw, (double)y1, a,
              (double)x0, (double)y1, a,
              (double)x0, y1+yh, a);
    multiline(3, x1+xw, (double)y1, a,
              (double)x1, (double)y1, a,
              (double)x1, y1+yh, a);
    multiline(3, x0-xw, (double)y0, a,
              (double)x0, (double)y0, a,
              (double)x0, y0-yh, a);
    multiline(3, x1+xw, (double)y0, a,
              (double)x1, (double)y0, a,
              (double)x1, y0-yh, a);
    if(a == z3)
      break;
    else
      a = z3;
  }

  /* Draw cylinder. */
  circle("xy", 0, -r0, z1, r);
  circle("xy", 0, -r0, z2, r);
  line(0, -r0, z1, 0, -r0, z2);
  for(a = 0; a < 2*PI; a += PI/8)
  {
    multiline(4,
              0.0, -r0, z1,
              r*cos(a), r*sin(a) - r0, z1,
              r*cos(a + DEG2RAD*phi), r*sin(a + DEG2RAD*phi) - r0, z2,
              0.0, -r0, z2);
  }
  /*
  
  multiline(5, (double)xmin, (double)ymin, 0.0,
               (double)xmax, (double)ymin, 0.0,
               (double)xmax, (double)ymax, 0.0,
               (double)xmin, (double)ymax, 0.0,
               (double)xmin, (double)ymin, 0.0); */
}
#line 19583 "./HZB_FLEX.c"
}   /* End of SELEC=Selector() SETTING parameter declarations. */
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'After_selec'. */
  SIG_MESSAGE("After_selec (McDisplay)");
  printf("MCDISPLAY: component %s\n", "After_selec");
#define mccompcurname  After_selec
#define mccompcurtype  PSD_monitor
#define mccompcurindex 12
#define PSD_N mccAfter_selec_PSD_N
#define PSD_p mccAfter_selec_PSD_p
#define PSD_p2 mccAfter_selec_PSD_p2
{   /* Declarations of After_selec=PSD_monitor() SETTING parameters. */
int nx = mccAfter_selec_nx;
int ny = mccAfter_selec_ny;
char* filename = mccAfter_selec_filename;
MCNUM xmin = mccAfter_selec_xmin;
MCNUM xmax = mccAfter_selec_xmax;
MCNUM ymin = mccAfter_selec_ymin;
MCNUM ymax = mccAfter_selec_ymax;
MCNUM xwidth = mccAfter_selec_xwidth;
MCNUM yheight = mccAfter_selec_yheight;
MCNUM restore_neutron = mccAfter_selec_restore_neutron;
#line 129 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
{
  multiline(5,
    (double)xmin, (double)ymin, 0.0,
    (double)xmax, (double)ymin, 0.0,
    (double)xmax, (double)ymax, 0.0,
    (double)xmin, (double)ymax, 0.0,
    (double)xmin, (double)ymin, 0.0);
}
#line 19618 "./HZB_FLEX.c"
}   /* End of After_selec=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'NL1B_Curved_2_Entrance_Window'. */
  SIG_MESSAGE("NL1B_Curved_2_Entrance_Window (McDisplay)");
  printf("MCDISPLAY: component %s\n", "NL1B_Curved_2_Entrance_Window");
#define mccompcurname  NL1B_Curved_2_Entrance_Window
#define mccompcurtype  Arm
#define mccompcurindex 13
#line 40 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Arm.comp"
{
  /* A bit ugly; hard-coded dimensions. */
  
  line(0,0,0,0.2,0,0);
  line(0,0,0,0,0.2,0);
  line(0,0,0,0,0,0.2);
}
#line 19641 "./HZB_FLEX.c"
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'NL1B_Curved_Guide_2'. */
  SIG_MESSAGE("NL1B_Curved_Guide_2 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "NL1B_Curved_Guide_2");
#define mccompcurname  NL1B_Curved_Guide_2
#define mccompcurtype  Guide_curved
#define mccompcurindex 14
{   /* Declarations of NL1B_Curved_Guide_2=Guide_curved() SETTING parameters. */
MCNUM w1 = mccNL1B_Curved_Guide_2_w1;
MCNUM h1 = mccNL1B_Curved_Guide_2_h1;
MCNUM l = mccNL1B_Curved_Guide_2_l;
MCNUM R0 = mccNL1B_Curved_Guide_2_R0;
MCNUM Qc = mccNL1B_Curved_Guide_2_Qc;
MCNUM alpha = mccNL1B_Curved_Guide_2_alpha;
MCNUM m = mccNL1B_Curved_Guide_2_m;
MCNUM W = mccNL1B_Curved_Guide_2_W;
MCNUM curvature = mccNL1B_Curved_Guide_2_curvature;
#line 133 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../contrib/Guide_curved.comp"
{
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
}
#line 19698 "./HZB_FLEX.c"
}   /* End of NL1B_Curved_Guide_2=Guide_curved() SETTING parameter declarations. */
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'NL1B_Straight2_Entrance_Window'. */
  SIG_MESSAGE("NL1B_Straight2_Entrance_Window (McDisplay)");
  printf("MCDISPLAY: component %s\n", "NL1B_Straight2_Entrance_Window");
#define mccompcurname  NL1B_Straight2_Entrance_Window
#define mccompcurtype  Arm
#define mccompcurindex 15
#line 40 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Arm.comp"
{
  /* A bit ugly; hard-coded dimensions. */
  
  line(0,0,0,0.2,0,0);
  line(0,0,0,0,0.2,0);
  line(0,0,0,0,0,0.2);
}
#line 19718 "./HZB_FLEX.c"
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'NL1B_Straight2'. */
  SIG_MESSAGE("NL1B_Straight2 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "NL1B_Straight2");
#define mccompcurname  NL1B_Straight2
#define mccompcurtype  Guide
#define mccompcurindex 16
#define pTable mccNL1B_Straight2_pTable
{   /* Declarations of NL1B_Straight2=Guide() SETTING parameters. */
char* reflect = mccNL1B_Straight2_reflect;
MCNUM w1 = mccNL1B_Straight2_w1;
MCNUM h1 = mccNL1B_Straight2_h1;
MCNUM w2 = mccNL1B_Straight2_w2;
MCNUM h2 = mccNL1B_Straight2_h2;
MCNUM l = mccNL1B_Straight2_l;
MCNUM R0 = mccNL1B_Straight2_R0;
MCNUM Qc = mccNL1B_Straight2_Qc;
MCNUM alpha = mccNL1B_Straight2_alpha;
MCNUM m = mccNL1B_Straight2_m;
MCNUM W = mccNL1B_Straight2_W;
#line 202 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide.comp"
{
  
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
}
#line 19762 "./HZB_FLEX.c"
}   /* End of NL1B_Straight2=Guide() SETTING parameter declarations. */
#undef pTable
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'NL1B_Elliptical_Entrance_Window'. */
  SIG_MESSAGE("NL1B_Elliptical_Entrance_Window (McDisplay)");
  printf("MCDISPLAY: component %s\n", "NL1B_Elliptical_Entrance_Window");
#define mccompcurname  NL1B_Elliptical_Entrance_Window
#define mccompcurtype  Arm
#define mccompcurindex 17
#line 40 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Arm.comp"
{
  /* A bit ugly; hard-coded dimensions. */
  
  line(0,0,0,0.2,0,0);
  line(0,0,0,0,0.2,0);
  line(0,0,0,0,0,0.2);
}
#line 19783 "./HZB_FLEX.c"
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'elliptical_piece'. */
  SIG_MESSAGE("elliptical_piece (McDisplay)");
  printf("MCDISPLAY: component %s\n", "elliptical_piece");
#define mccompcurname  elliptical_piece
#define mccompcurtype  Guide_tapering
#define mccompcurindex 18
#define w1c mccelliptical_piece_w1c
#define w2c mccelliptical_piece_w2c
#define ww mccelliptical_piece_ww
#define hh mccelliptical_piece_hh
#define whalf mccelliptical_piece_whalf
#define hhalf mccelliptical_piece_hhalf
#define lwhalf mccelliptical_piece_lwhalf
#define lhhalf mccelliptical_piece_lhhalf
#define h1_in mccelliptical_piece_h1_in
#define h2_out mccelliptical_piece_h2_out
#define w1_in mccelliptical_piece_w1_in
#define w2_out mccelliptical_piece_w2_out
#define l_seg mccelliptical_piece_l_seg
#define seg mccelliptical_piece_seg
#define h12 mccelliptical_piece_h12
#define h2 mccelliptical_piece_h2
#define w12 mccelliptical_piece_w12
#define w2 mccelliptical_piece_w2
#define a_ell_q mccelliptical_piece_a_ell_q
#define b_ell_q mccelliptical_piece_b_ell_q
#define lbw mccelliptical_piece_lbw
#define lbh mccelliptical_piece_lbh
#define mxi mccelliptical_piece_mxi
#define u1 mccelliptical_piece_u1
#define u2 mccelliptical_piece_u2
#define div1 mccelliptical_piece_div1
#define p2_para mccelliptical_piece_p2_para
#define test mccelliptical_piece_test
#define Div1 mccelliptical_piece_Div1
#define i mccelliptical_piece_i
#define ii mccelliptical_piece_ii
#define seg mccelliptical_piece_seg
#define fu mccelliptical_piece_fu
#define pos mccelliptical_piece_pos
#define file_name mccelliptical_piece_file_name
#define ep mccelliptical_piece_ep
#define num mccelliptical_piece_num
#define rotation_h mccelliptical_piece_rotation_h
#define rotation_v mccelliptical_piece_rotation_v
{   /* Declarations of elliptical_piece=Guide_tapering() SETTING parameters. */
char* option = mccelliptical_piece_option;
MCNUM w1 = mccelliptical_piece_w1;
MCNUM h1 = mccelliptical_piece_h1;
MCNUM l = mccelliptical_piece_l;
MCNUM linw = mccelliptical_piece_linw;
MCNUM loutw = mccelliptical_piece_loutw;
MCNUM linh = mccelliptical_piece_linh;
MCNUM louth = mccelliptical_piece_louth;
MCNUM R0 = mccelliptical_piece_R0;
MCNUM Qcx = mccelliptical_piece_Qcx;
MCNUM Qcy = mccelliptical_piece_Qcy;
MCNUM alphax = mccelliptical_piece_alphax;
MCNUM alphay = mccelliptical_piece_alphay;
MCNUM W = mccelliptical_piece_W;
MCNUM mx = mccelliptical_piece_mx;
MCNUM my = mccelliptical_piece_my;
MCNUM segno = mccelliptical_piece_segno;
MCNUM curvature = mccelliptical_piece_curvature;
MCNUM curvature_v = mccelliptical_piece_curvature_v;
#line 625 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide_tapering.comp"
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
#line 19886 "./HZB_FLEX.c"
}   /* End of elliptical_piece=Guide_tapering() SETTING parameter declarations. */
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

  /* MCDISPLAY code for component 'Virtual_source'. */
  SIG_MESSAGE("Virtual_source (McDisplay)");
  printf("MCDISPLAY: component %s\n", "Virtual_source");
#define mccompcurname  Virtual_source
#define mccompcurtype  Slit
#define mccompcurindex 19
{   /* Declarations of Virtual_source=Slit() SETTING parameters. */
MCNUM xmin = mccVirtual_source_xmin;
MCNUM xmax = mccVirtual_source_xmax;
MCNUM ymin = mccVirtual_source_ymin;
MCNUM ymax = mccVirtual_source_ymax;
MCNUM radius = mccVirtual_source_radius;
MCNUM xwidth = mccVirtual_source_xwidth;
MCNUM yheight = mccVirtual_source_yheight;
#line 83 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Slit.comp"
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
#line 19968 "./HZB_FLEX.c"
}   /* End of Virtual_source=Slit() SETTING parameter declarations. */
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'NL1B_Vertical_Guide_Entrance_Window'. */
  SIG_MESSAGE("NL1B_Vertical_Guide_Entrance_Window (McDisplay)");
  printf("MCDISPLAY: component %s\n", "NL1B_Vertical_Guide_Entrance_Window");
#define mccompcurname  NL1B_Vertical_Guide_Entrance_Window
#define mccompcurtype  Arm
#define mccompcurindex 20
#line 40 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Arm.comp"
{
  /* A bit ugly; hard-coded dimensions. */
  
  line(0,0,0,0.2,0,0);
  line(0,0,0,0,0.2,0);
  line(0,0,0,0,0,0.2);
}
#line 19988 "./HZB_FLEX.c"
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'NL1B_Vertical_Guide'. */
  SIG_MESSAGE("NL1B_Vertical_Guide (McDisplay)");
  printf("MCDISPLAY: component %s\n", "NL1B_Vertical_Guide");
#define mccompcurname  NL1B_Vertical_Guide
#define mccompcurtype  Guide_channeled
#define mccompcurindex 21
#define w1c mccNL1B_Vertical_Guide_w1c
#define w2c mccNL1B_Vertical_Guide_w2c
#define ww mccNL1B_Vertical_Guide_ww
#define hh mccNL1B_Vertical_Guide_hh
#define whalf mccNL1B_Vertical_Guide_whalf
#define hhalf mccNL1B_Vertical_Guide_hhalf
#define lwhalf mccNL1B_Vertical_Guide_lwhalf
#define lhhalf mccNL1B_Vertical_Guide_lhhalf
{   /* Declarations of NL1B_Vertical_Guide=Guide_channeled() SETTING parameters. */
MCNUM w1 = mccNL1B_Vertical_Guide_w1;
MCNUM h1 = mccNL1B_Vertical_Guide_h1;
MCNUM w2 = mccNL1B_Vertical_Guide_w2;
MCNUM h2 = mccNL1B_Vertical_Guide_h2;
MCNUM l = mccNL1B_Vertical_Guide_l;
MCNUM R0 = mccNL1B_Vertical_Guide_R0;
MCNUM Qc = mccNL1B_Vertical_Guide_Qc;
MCNUM alpha = mccNL1B_Vertical_Guide_alpha;
MCNUM m = mccNL1B_Vertical_Guide_m;
MCNUM nslit = mccNL1B_Vertical_Guide_nslit;
MCNUM d = mccNL1B_Vertical_Guide_d;
MCNUM Qcx = mccNL1B_Vertical_Guide_Qcx;
MCNUM Qcy = mccNL1B_Vertical_Guide_Qcy;
MCNUM alphax = mccNL1B_Vertical_Guide_alphax;
MCNUM alphay = mccNL1B_Vertical_Guide_alphay;
MCNUM W = mccNL1B_Vertical_Guide_W;
MCNUM mx = mccNL1B_Vertical_Guide_mx;
MCNUM my = mccNL1B_Vertical_Guide_my;
MCNUM nu = mccNL1B_Vertical_Guide_nu;
MCNUM phase = mccNL1B_Vertical_Guide_phase;
#line 297 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide_channeled.comp"
{
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
}
#line 20059 "./HZB_FLEX.c"
}   /* End of NL1B_Vertical_Guide=Guide_channeled() SETTING parameter declarations. */
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

  /* MCDISPLAY code for component 'NL1B_Guide_Exit'. */
  SIG_MESSAGE("NL1B_Guide_Exit (McDisplay)");
  printf("MCDISPLAY: component %s\n", "NL1B_Guide_Exit");
#define mccompcurname  NL1B_Guide_Exit
#define mccompcurtype  Arm
#define mccompcurindex 22
#line 40 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Arm.comp"
{
  /* A bit ugly; hard-coded dimensions. */
  
  line(0,0,0,0.2,0,0);
  line(0,0,0,0,0.2,0);
  line(0,0,0,0,0,0.2);
}
#line 20087 "./HZB_FLEX.c"
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'energy_endguide'. */
  SIG_MESSAGE("energy_endguide (McDisplay)");
  printf("MCDISPLAY: component %s\n", "energy_endguide");
#define mccompcurname  energy_endguide
#define mccompcurtype  E_monitor
#define mccompcurindex 23
#define nE mccenergy_endguide_nE
#define E_N mccenergy_endguide_E_N
#define E_p mccenergy_endguide_E_p
#define E_p2 mccenergy_endguide_E_p2
#define S_p mccenergy_endguide_S_p
#define S_pE mccenergy_endguide_S_pE
#define S_pE2 mccenergy_endguide_S_pE2
{   /* Declarations of energy_endguide=E_monitor() SETTING parameters. */
char* filename = mccenergy_endguide_filename;
MCNUM xmin = mccenergy_endguide_xmin;
MCNUM xmax = mccenergy_endguide_xmax;
MCNUM ymin = mccenergy_endguide_ymin;
MCNUM ymax = mccenergy_endguide_ymax;
MCNUM xwidth = mccenergy_endguide_xwidth;
MCNUM yheight = mccenergy_endguide_yheight;
MCNUM Emin = mccenergy_endguide_Emin;
MCNUM Emax = mccenergy_endguide_Emax;
MCNUM restore_neutron = mccenergy_endguide_restore_neutron;
int nowritefile = mccenergy_endguide_nowritefile;
#line 132 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/E_monitor.comp"
{
  
  multiline(5, (double)xmin, (double)ymin, 0.0,
               (double)xmax, (double)ymin, 0.0,
               (double)xmax, (double)ymax, 0.0,
               (double)xmin, (double)ymax, 0.0,
               (double)xmin, (double)ymin, 0.0);
}
#line 20126 "./HZB_FLEX.c"
}   /* End of energy_endguide=E_monitor() SETTING parameter declarations. */
#undef S_pE2
#undef S_pE
#undef S_p
#undef E_p2
#undef E_p
#undef E_N
#undef nE
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'Mono_center'. */
  SIG_MESSAGE("Mono_center (McDisplay)");
  printf("MCDISPLAY: component %s\n", "Mono_center");
#define mccompcurname  Mono_center
#define mccompcurtype  Arm
#define mccompcurindex 24
#line 40 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Arm.comp"
{
  /* A bit ugly; hard-coded dimensions. */
  
  line(0,0,0,0.2,0,0);
  line(0,0,0,0,0.2,0);
  line(0,0,0,0,0,0.2);
}
#line 20153 "./HZB_FLEX.c"
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'Monochromator'. */
  SIG_MESSAGE("Monochromator (McDisplay)");
  printf("MCDISPLAY: component %s\n", "Monochromator");
#define mccompcurname  Monochromator
#define mccompcurtype  Monochromator_curved
#define mccompcurindex 25
#define mos_rms_y mccMonochromator_mos_rms_y
#define mos_rms_z mccMonochromator_mos_rms_z
#define mos_rms_max mccMonochromator_mos_rms_max
#define mono_Q mccMonochromator_mono_Q
#define SlabWidth mccMonochromator_SlabWidth
#define SlabHeight mccMonochromator_SlabHeight
#define rTable mccMonochromator_rTable
#define tTable mccMonochromator_tTable
#define row mccMonochromator_row
#define col mccMonochromator_col
#define tiltH mccMonochromator_tiltH
#define tiltV mccMonochromator_tiltV
{   /* Declarations of Monochromator=Monochromator_curved() SETTING parameters. */
char* reflect = mccMonochromator_reflect;
char* transmit = mccMonochromator_transmit;
MCNUM zwidth = mccMonochromator_zwidth;
MCNUM yheight = mccMonochromator_yheight;
MCNUM gap = mccMonochromator_gap;
MCNUM NH = mccMonochromator_NH;
MCNUM NV = mccMonochromator_NV;
MCNUM mosaich = mccMonochromator_mosaich;
MCNUM mosaicv = mccMonochromator_mosaicv;
MCNUM r0 = mccMonochromator_r0;
MCNUM t0 = mccMonochromator_t0;
MCNUM Q = mccMonochromator_Q;
MCNUM RV = mccMonochromator_RV;
MCNUM RH = mccMonochromator_RH;
MCNUM DM = mccMonochromator_DM;
MCNUM mosaic = mccMonochromator_mosaic;
MCNUM width = mccMonochromator_width;
MCNUM height = mccMonochromator_height;
MCNUM verbose = mccMonochromator_verbose;
MCNUM order = mccMonochromator_order;
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
#line 20227 "./HZB_FLEX.c"
}   /* End of Monochromator=Monochromator_curved() SETTING parameter declarations. */
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

  /* MCDISPLAY code for component 'Mono_sample_arm'. */
  SIG_MESSAGE("Mono_sample_arm (McDisplay)");
  printf("MCDISPLAY: component %s\n", "Mono_sample_arm");
#define mccompcurname  Mono_sample_arm
#define mccompcurtype  Arm
#define mccompcurindex 26
#line 40 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Arm.comp"
{
  /* A bit ugly; hard-coded dimensions. */
  
  line(0,0,0,0.2,0,0);
  line(0,0,0,0,0.2,0);
  line(0,0,0,0,0,0.2);
}
#line 20259 "./HZB_FLEX.c"
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'energy_mono'. */
  SIG_MESSAGE("energy_mono (McDisplay)");
  printf("MCDISPLAY: component %s\n", "energy_mono");
#define mccompcurname  energy_mono
#define mccompcurtype  E_monitor
#define mccompcurindex 27
#define nE mccenergy_mono_nE
#define E_N mccenergy_mono_E_N
#define E_p mccenergy_mono_E_p
#define E_p2 mccenergy_mono_E_p2
#define S_p mccenergy_mono_S_p
#define S_pE mccenergy_mono_S_pE
#define S_pE2 mccenergy_mono_S_pE2
{   /* Declarations of energy_mono=E_monitor() SETTING parameters. */
char* filename = mccenergy_mono_filename;
MCNUM xmin = mccenergy_mono_xmin;
MCNUM xmax = mccenergy_mono_xmax;
MCNUM ymin = mccenergy_mono_ymin;
MCNUM ymax = mccenergy_mono_ymax;
MCNUM xwidth = mccenergy_mono_xwidth;
MCNUM yheight = mccenergy_mono_yheight;
MCNUM Emin = mccenergy_mono_Emin;
MCNUM Emax = mccenergy_mono_Emax;
MCNUM restore_neutron = mccenergy_mono_restore_neutron;
int nowritefile = mccenergy_mono_nowritefile;
#line 132 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/E_monitor.comp"
{
  
  multiline(5, (double)xmin, (double)ymin, 0.0,
               (double)xmax, (double)ymin, 0.0,
               (double)xmax, (double)ymax, 0.0,
               (double)xmin, (double)ymax, 0.0,
               (double)xmin, (double)ymin, 0.0);
}
#line 20298 "./HZB_FLEX.c"
}   /* End of energy_mono=E_monitor() SETTING parameter declarations. */
#undef S_pE2
#undef S_pE
#undef S_p
#undef E_p2
#undef E_p
#undef E_N
#undef nE
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'energy_pre_sample'. */
  SIG_MESSAGE("energy_pre_sample (McDisplay)");
  printf("MCDISPLAY: component %s\n", "energy_pre_sample");
#define mccompcurname  energy_pre_sample
#define mccompcurtype  E_monitor
#define mccompcurindex 28
#define nE mccenergy_pre_sample_nE
#define E_N mccenergy_pre_sample_E_N
#define E_p mccenergy_pre_sample_E_p
#define E_p2 mccenergy_pre_sample_E_p2
#define S_p mccenergy_pre_sample_S_p
#define S_pE mccenergy_pre_sample_S_pE
#define S_pE2 mccenergy_pre_sample_S_pE2
{   /* Declarations of energy_pre_sample=E_monitor() SETTING parameters. */
char* filename = mccenergy_pre_sample_filename;
MCNUM xmin = mccenergy_pre_sample_xmin;
MCNUM xmax = mccenergy_pre_sample_xmax;
MCNUM ymin = mccenergy_pre_sample_ymin;
MCNUM ymax = mccenergy_pre_sample_ymax;
MCNUM xwidth = mccenergy_pre_sample_xwidth;
MCNUM yheight = mccenergy_pre_sample_yheight;
MCNUM Emin = mccenergy_pre_sample_Emin;
MCNUM Emax = mccenergy_pre_sample_Emax;
MCNUM restore_neutron = mccenergy_pre_sample_restore_neutron;
int nowritefile = mccenergy_pre_sample_nowritefile;
#line 132 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/E_monitor.comp"
{
  
  multiline(5, (double)xmin, (double)ymin, 0.0,
               (double)xmax, (double)ymin, 0.0,
               (double)xmax, (double)ymax, 0.0,
               (double)xmin, (double)ymax, 0.0,
               (double)xmin, (double)ymin, 0.0);
}
#line 20345 "./HZB_FLEX.c"
}   /* End of energy_pre_sample=E_monitor() SETTING parameter declarations. */
#undef S_pE2
#undef S_pE
#undef S_p
#undef E_p2
#undef E_p
#undef E_N
#undef nE
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'Sample_center'. */
  SIG_MESSAGE("Sample_center (McDisplay)");
  printf("MCDISPLAY: component %s\n", "Sample_center");
#define mccompcurname  Sample_center
#define mccompcurtype  Arm
#define mccompcurindex 29
#line 40 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Arm.comp"
{
  /* A bit ugly; hard-coded dimensions. */
  
  line(0,0,0,0.2,0,0);
  line(0,0,0,0,0.2,0);
  line(0,0,0,0,0,0.2);
}
#line 20372 "./HZB_FLEX.c"
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'Sample_analyser_arm'. */
  SIG_MESSAGE("Sample_analyser_arm (McDisplay)");
  printf("MCDISPLAY: component %s\n", "Sample_analyser_arm");
#define mccompcurname  Sample_analyser_arm
#define mccompcurtype  Arm
#define mccompcurindex 30
#line 40 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Arm.comp"
{
  /* A bit ugly; hard-coded dimensions. */
  
  line(0,0,0,0.2,0,0);
  line(0,0,0,0,0.2,0);
  line(0,0,0,0,0,0.2);
}
#line 20391 "./HZB_FLEX.c"
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'div_mono'. */
  SIG_MESSAGE("div_mono (McDisplay)");
  printf("MCDISPLAY: component %s\n", "div_mono");
#define mccompcurname  div_mono
#define mccompcurtype  Divergence_monitor
#define mccompcurindex 31
#define nh mccdiv_mono_nh
#define nv mccdiv_mono_nv
#define Div_N mccdiv_mono_Div_N
#define Div_p mccdiv_mono_Div_p
#define Div_p2 mccdiv_mono_Div_p2
{   /* Declarations of div_mono=Divergence_monitor() SETTING parameters. */
char* filename = mccdiv_mono_filename;
MCNUM xmin = mccdiv_mono_xmin;
MCNUM xmax = mccdiv_mono_xmax;
MCNUM ymin = mccdiv_mono_ymin;
MCNUM ymax = mccdiv_mono_ymax;
MCNUM xwidth = mccdiv_mono_xwidth;
MCNUM yheight = mccdiv_mono_yheight;
MCNUM maxdiv_h = mccdiv_mono_maxdiv_h;
MCNUM maxdiv_v = mccdiv_mono_maxdiv_v;
MCNUM restore_neutron = mccdiv_mono_restore_neutron;
MCNUM nx = mccdiv_mono_nx;
MCNUM ny = mccdiv_mono_ny;
MCNUM nz = mccdiv_mono_nz;
int nowritefile = mccdiv_mono_nowritefile;
#line 131 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/Divergence_monitor.comp"
{
    
    multiline(5, (double)xmin, (double)ymin, 0.0,
                 (double)xmax, (double)ymin, 0.0,
                 (double)xmax, (double)ymax, 0.0,
                 (double)xmin, (double)ymax, 0.0,
                 (double)xmin, (double)ymin, 0.0);
}
#line 20431 "./HZB_FLEX.c"
}   /* End of div_mono=Divergence_monitor() SETTING parameter declarations. */
#undef Div_p2
#undef Div_p
#undef Div_N
#undef nv
#undef nh
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'div_mono_H'. */
  SIG_MESSAGE("div_mono_H (McDisplay)");
  printf("MCDISPLAY: component %s\n", "div_mono_H");
#define mccompcurname  div_mono_H
#define mccompcurtype  Hdiv_monitor
#define mccompcurindex 32
#define nh mccdiv_mono_H_nh
#define Div_N mccdiv_mono_H_Div_N
#define Div_p mccdiv_mono_H_Div_p
#define Div_p2 mccdiv_mono_H_Div_p2
{   /* Declarations of div_mono_H=Hdiv_monitor() SETTING parameters. */
char* filename = mccdiv_mono_H_filename;
MCNUM xmin = mccdiv_mono_H_xmin;
MCNUM xmax = mccdiv_mono_H_xmax;
MCNUM ymin = mccdiv_mono_H_ymin;
MCNUM ymax = mccdiv_mono_H_ymax;
MCNUM xwidth = mccdiv_mono_H_xwidth;
MCNUM yheight = mccdiv_mono_H_yheight;
MCNUM h_maxdiv = mccdiv_mono_H_h_maxdiv;
MCNUM restore_neutron = mccdiv_mono_H_restore_neutron;
int nowritefile = mccdiv_mono_H_nowritefile;
#line 120 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/Hdiv_monitor.comp"
{
    
    multiline(5, (double)xmin, (double)ymin, 0.0,
                 (double)xmax, (double)ymin, 0.0,
                 (double)xmax, (double)ymax, 0.0,
                 (double)xmin, (double)ymax, 0.0,
                 (double)xmin, (double)ymin, 0.0);
}
#line 20472 "./HZB_FLEX.c"
}   /* End of div_mono_H=Hdiv_monitor() SETTING parameter declarations. */
#undef Div_p2
#undef Div_p
#undef Div_N
#undef nh
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'psd_sam'. */
  SIG_MESSAGE("psd_sam (McDisplay)");
  printf("MCDISPLAY: component %s\n", "psd_sam");
#define mccompcurname  psd_sam
#define mccompcurtype  PSD_monitor
#define mccompcurindex 33
#define PSD_N mccpsd_sam_PSD_N
#define PSD_p mccpsd_sam_PSD_p
#define PSD_p2 mccpsd_sam_PSD_p2
{   /* Declarations of psd_sam=PSD_monitor() SETTING parameters. */
int nx = mccpsd_sam_nx;
int ny = mccpsd_sam_ny;
char* filename = mccpsd_sam_filename;
MCNUM xmin = mccpsd_sam_xmin;
MCNUM xmax = mccpsd_sam_xmax;
MCNUM ymin = mccpsd_sam_ymin;
MCNUM ymax = mccpsd_sam_ymax;
MCNUM xwidth = mccpsd_sam_xwidth;
MCNUM yheight = mccpsd_sam_yheight;
MCNUM restore_neutron = mccpsd_sam_restore_neutron;
#line 129 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
{
  multiline(5,
    (double)xmin, (double)ymin, 0.0,
    (double)xmax, (double)ymin, 0.0,
    (double)xmax, (double)ymax, 0.0,
    (double)xmin, (double)ymax, 0.0,
    (double)xmin, (double)ymin, 0.0);
}
#line 20511 "./HZB_FLEX.c"
}   /* End of psd_sam=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component '1dpsd'. */
  SIG_MESSAGE("1dpsd (McDisplay)");
  printf("MCDISPLAY: component %s\n", "1dpsd");
#define mccompcurname  1dpsd
#define mccompcurtype  PSDlin_monitor
#define mccompcurindex 34
#define nx mcc1dpsd_nx
#define PSDlin_N mcc1dpsd_PSDlin_N
#define PSDlin_p mcc1dpsd_PSDlin_p
#define PSDlin_p2 mcc1dpsd_PSDlin_p2
{   /* Declarations of 1dpsd=PSDlin_monitor() SETTING parameters. */
char* filename = mcc1dpsd_filename;
MCNUM xmin = mcc1dpsd_xmin;
MCNUM xmax = mcc1dpsd_xmax;
MCNUM ymin = mcc1dpsd_ymin;
MCNUM ymax = mcc1dpsd_ymax;
MCNUM xwidth = mcc1dpsd_xwidth;
MCNUM yheight = mcc1dpsd_yheight;
MCNUM restore_neutron = mcc1dpsd_restore_neutron;
int nowritefile = mcc1dpsd_nowritefile;
#line 116 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/PSDlin_monitor.comp"
{
  
  multiline(5, (double)xmin, (double)ymin, 0.0,
               (double)xmax, (double)ymin, 0.0,
               (double)xmax, (double)ymax, 0.0,
               (double)xmin, (double)ymax, 0.0,
               (double)xmin, (double)ymin, 0.0);
}
#line 20549 "./HZB_FLEX.c"
}   /* End of 1dpsd=PSDlin_monitor() SETTING parameter declarations. */
#undef PSDlin_p2
#undef PSDlin_p
#undef PSDlin_N
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
/* end of generated C code ./HZB_FLEX.c */