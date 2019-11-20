/* Automatically generated file. Do not edit. 
 * Format:     ANSI C source code
 * Creator:    McStas <http://www.mcstas.org>
 * Instrument: /zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr (ESS_2015_test)
 * Date:       Wed Nov 20 00:11:12 2019
 * File:       ./ESS_2015_test.c
 * Compile:    cc -o ESS_2015_test.out ./ESS_2015_test.c 
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

#line 712 "./ESS_2015_test.c"

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

#line 945 "./ESS_2015_test.c"

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

#line 4977 "./ESS_2015_test.c"

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

#line 5337 "./ESS_2015_test.c"
#ifdef MC_TRACE_ENABLED
int mctraceenabled = 1;
#else
int mctraceenabled = 0;
#endif
#define MCSTAS "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../"
int mcdefaultmain = 1;
char mcinstrument_name[] = "ESS_2015_test";
char mcinstrument_source[] = "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr";
char *mcinstrument_exe=NULL; /* will be set to argv[0] in main */
int main(int argc, char *argv[]){return mccode_main(argc, argv);}
void mcinit(void);
void mcraytrace(void);
void mcsave(FILE *);
void mcfinally(void);
void mcdisplay(void);

/* Shared user declarations for all components 'ESS_moderator'. */
#line 91 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../sources/ESS_moderator.comp"
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

#line 6365 "./ESS_2015_test.c"

/* Instrument parameters. */
MCNUM mcipfrac;
MCNUM mciplmin;
MCNUM mciplmax;
MCNUM mcipANGLE;
MCNUM mcipmon_shift;

#define mcNUMIPAR 5
int mcnumipar = 5;
struct mcinputtable_struct mcinputtable[mcNUMIPAR+1] = {
  "frac", &mcipfrac, instr_type_double, "0.5", 
  "lmin", &mciplmin, instr_type_double, "0.3", 
  "lmax", &mciplmax, instr_type_double, "8.3", 
  "ANGLE", &mcipANGLE, instr_type_double, "-55", 
  "mon_shift", &mcipmon_shift, instr_type_double, "0.0", 
  NULL, NULL, instr_type_double, ""
};

/* User declarations from instrument definition. */
#define mccompcurname  ESS_2015_test
#define mccompcurtype  INSTRUMENT
#define mccompcurindex 0
#define mcposaESS_2015_test coords_set(0,0,0)
#define frac mcipfrac
#define lmin mciplmin
#define lmax mciplmax
#define ANGLE mcipANGLE
#define mon_shift mcipmon_shift
#line 40 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
  //double lambdamin,lambdamax;
  int IsCold,np,sgn;
  char srcdef[128];
  double Yheight,power;
  double cosa,sina;
  double ctrX,ctrZ,BX0,BZ0;
  double mon_width;
#line 6403 "./ESS_2015_test.c"
#undef mon_shift
#undef ANGLE
#undef lmax
#undef lmin
#undef frac
#undef mcposaESS_2015_test
#undef mccompcurindex
#undef mccompcurtype
#undef mccompcurname

/* neutron state table at each component input (local coords) */
/* [x, y, z, vx, vy, vz, t, sx, sy, sz, p] */
MCNUM mccomp_storein[11*12];
/* Components position table (absolute and relative coords) */
Coords mccomp_posa[12];
Coords mccomp_posr[12];
/* Counter for each comp to check for inactive ones */
MCNUM  mcNCounter[12];
MCNUM  mcPCounter[12];
MCNUM  mcP2Counter[12];
#define mcNUMCOMP 11 /* number of components */
/* Counter for PROP ABSORB */
MCNUM  mcAbsorbProp[12];
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
MCNUM mccSource_isleft;
MCNUM mccSource_Lmin;
MCNUM mccSource_Lmax;
MCNUM mccSource_cold_frac;
MCNUM mccSource_dist;
MCNUM mccSource_focus_xw;
MCNUM mccSource_focus_yh;
int mccSource_target_index;
MCNUM mccSource_tmax_multiplier;
MCNUM mccSource_yheight_c;
MCNUM mccSource_yheight_t;
int mccSource_n_pulses;
MCNUM mccSource_acc_power;
MCNUM mccSource_beamport_angle;
char mccSource_sourcedef[16384];
MCNUM mccSource_xwidth_c;
MCNUM mccSource_xwidth_t;
MCNUM mccSource_extraction_opening;

/* Definition parameters for component 'XmonT' [4]. */
#define mccXmonT_nx 250
/* Setting parameters for component 'XmonT' [4]. */
char mccXmonT_filename[16384];
MCNUM mccXmonT_xmin;
MCNUM mccXmonT_xmax;
MCNUM mccXmonT_ymin;
MCNUM mccXmonT_ymax;
MCNUM mccXmonT_xwidth;
MCNUM mccXmonT_yheight;
MCNUM mccXmonT_restore_neutron;
int mccXmonT_nowritefile;

/* Definition parameters for component 'XmonC' [5]. */
#define mccXmonC_nx 250
/* Setting parameters for component 'XmonC' [5]. */
char mccXmonC_filename[16384];
MCNUM mccXmonC_xmin;
MCNUM mccXmonC_xmax;
MCNUM mccXmonC_ymin;
MCNUM mccXmonC_ymax;
MCNUM mccXmonC_xwidth;
MCNUM mccXmonC_yheight;
MCNUM mccXmonC_restore_neutron;
int mccXmonC_nowritefile;

/* Definition parameters for component 'tofT' [6]. */
#define mcctofT_nt 100
/* Setting parameters for component 'tofT' [6]. */
char mcctofT_filename[16384];
MCNUM mcctofT_xmin;
MCNUM mcctofT_xmax;
MCNUM mcctofT_ymin;
MCNUM mcctofT_ymax;
MCNUM mcctofT_xwidth;
MCNUM mcctofT_yheight;
MCNUM mcctofT_tmin;
MCNUM mcctofT_tmax;
MCNUM mcctofT_dt;
MCNUM mcctofT_restore_neutron;
int mcctofT_nowritefile;

/* Definition parameters for component 'tofC' [7]. */
#define mcctofC_nt 100
/* Setting parameters for component 'tofC' [7]. */
char mcctofC_filename[16384];
MCNUM mcctofC_xmin;
MCNUM mcctofC_xmax;
MCNUM mcctofC_ymin;
MCNUM mcctofC_ymax;
MCNUM mcctofC_xwidth;
MCNUM mcctofC_yheight;
MCNUM mcctofC_tmin;
MCNUM mcctofC_tmax;
MCNUM mcctofC_dt;
MCNUM mcctofC_restore_neutron;
int mcctofC_nowritefile;

/* Setting parameters for component 'Slit1' [8]. */
MCNUM mccSlit1_xmin;
MCNUM mccSlit1_xmax;
MCNUM mccSlit1_ymin;
MCNUM mccSlit1_ymax;
MCNUM mccSlit1_radius;
MCNUM mccSlit1_xwidth;
MCNUM mccSlit1_yheight;

/* Definition parameters for component 'LamT' [9]. */
#define mccLamT_nL 160
/* Setting parameters for component 'LamT' [9]. */
char mccLamT_filename[16384];
MCNUM mccLamT_xmin;
MCNUM mccLamT_xmax;
MCNUM mccLamT_ymin;
MCNUM mccLamT_ymax;
MCNUM mccLamT_xwidth;
MCNUM mccLamT_yheight;
MCNUM mccLamT_Lmin;
MCNUM mccLamT_Lmax;
MCNUM mccLamT_restore_neutron;
int mccLamT_nowritefile;

/* Definition parameters for component 'LamC' [10]. */
#define mccLamC_nL 160
/* Setting parameters for component 'LamC' [10]. */
char mccLamC_filename[16384];
MCNUM mccLamC_xmin;
MCNUM mccLamC_xmax;
MCNUM mccLamC_ymin;
MCNUM mccLamC_ymax;
MCNUM mccLamC_xwidth;
MCNUM mccLamC_yheight;
MCNUM mccLamC_Lmin;
MCNUM mccLamC_Lmax;
MCNUM mccLamC_restore_neutron;
int mccLamC_nowritefile;

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
#line 6581 "./ESS_2015_test.c"
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
#define mccompcurtype  ESS_moderator
#define mccompcurindex 2
#define l_range mccSource_l_range
#define w_mult mccSource_w_mult
#define w_geom mccSource_w_geom
#define w_geom_c mccSource_w_geom_c
#define w_geom_t mccSource_w_geom_t
#define w_stat mccSource_w_stat
#define cold mccSource_cold
#define wasleft mccSource_wasleft
#define cosine mccSource_cosine
#define tx mccSource_tx
#define ty mccSource_ty
#define tz mccSource_tz
#define cx mccSource_cx
#define cy mccSource_cy
#define cz mccSource_cz
#define cyl_radius mccSource_cyl_radius
#define t1x mccSource_t1x
#define t1y mccSource_t1y
#define t1z mccSource_t1z
#define t2x mccSource_t2x
#define t2y mccSource_t2y
#define t2z mccSource_t2z
#define modextras mccSource_modextras
#define cold_bril mccSource_cold_bril
#define thermal_bril mccSource_thermal_bril
#define tfocus_width mccSource_tfocus_width
#define tfocus_time mccSource_tfocus_time
#define dt mccSource_dt
#define r_empty mccSource_r_empty
#define lambda mccSource_lambda
#define xprime mccSource_xprime
#define yprime mccSource_yprime
#define zprime mccSource_zprime
#define sgnport mccSource_sgnport
#define abs_beamport_angle mccSource_abs_beamport_angle
#define cos_beamport_angle mccSource_cos_beamport_angle
#define edge_thermal mccSource_edge_thermal
#define cos_thermal mccSource_cos_thermal
#define cos_cold mccSource_cos_cold
#define sin_beamport_angle mccSource_sin_beamport_angle
#define isleft mccSource_isleft
#define Lmin mccSource_Lmin
#define Lmax mccSource_Lmax
#define cold_frac mccSource_cold_frac
#define dist mccSource_dist
#define focus_xw mccSource_focus_xw
#define focus_yh mccSource_focus_yh
#define target_index mccSource_target_index
#define tmax_multiplier mccSource_tmax_multiplier
#define yheight_c mccSource_yheight_c
#define yheight_t mccSource_yheight_t
#define n_pulses mccSource_n_pulses
#define acc_power mccSource_acc_power
#define beamport_angle mccSource_beamport_angle
#define sourcedef mccSource_sourcedef
#define xwidth_c mccSource_xwidth_c
#define xwidth_t mccSource_xwidth_t
#define extraction_opening mccSource_extraction_opening
#line 96 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../sources/ESS_moderator.comp"
// Weighting-oriented scalars
double l_range, w_mult, w_geom, w_geom_c, w_geom_t, w_stat;
// Flag to indicate if originating from cold or thermal moderator
int cold;
// Flag to indicate if originating from left or right thermal slab
int wasleft;
// Flag to indicate if we should apply the solid angle correction from randvec
int cosine;
// Target coordinates (focusing)
  double tx,ty,tz;
  // Displacement of cold moderator centre from comp origin, cold moderator radius 2014 only!
  double cx,cy,cz,cyl_radius;
  // Vectors giving direction of the thermal slabs
  double t1x,t1y,t1z,t2x,t2y,t2z;

  /* Struct for transfer to ess-source-lib funcs  */
  ess_moderator_struct modextras;
  /* Cold and thermal function pointers */
  functype cold_bril;
  functype thermal_bril;
  /* Time-focusing parameters - currently unused */
  double tfocus_width,  tfocus_time,  dt;
  /* radius of vacated zone around the moderator */
  double r_empty;
  double lambda,xprime,yprime,zprime;
  double abs_beamport_angle,cos_beamport_angle,sin_beamport_angle;
  /* variables needed to correct for the emission surface angle*/
  double cos_thermal,cos_cold, edge_thermal;
  int sgnport;
#line 6686 "./ESS_2015_test.c"
#undef extraction_opening
#undef xwidth_t
#undef xwidth_c
#undef sourcedef
#undef beamport_angle
#undef acc_power
#undef n_pulses
#undef yheight_t
#undef yheight_c
#undef tmax_multiplier
#undef target_index
#undef focus_yh
#undef focus_xw
#undef dist
#undef cold_frac
#undef Lmax
#undef Lmin
#undef isleft
#undef sin_beamport_angle
#undef cos_cold
#undef cos_thermal
#undef edge_thermal
#undef cos_beamport_angle
#undef abs_beamport_angle
#undef sgnport
#undef zprime
#undef yprime
#undef xprime
#undef lambda
#undef r_empty
#undef dt
#undef tfocus_time
#undef tfocus_width
#undef thermal_bril
#undef cold_bril
#undef modextras
#undef t2z
#undef t2y
#undef t2x
#undef t1z
#undef t1y
#undef t1x
#undef cyl_radius
#undef cz
#undef cy
#undef cx
#undef tz
#undef ty
#undef tx
#undef cosine
#undef wasleft
#undef cold
#undef w_stat
#undef w_geom_t
#undef w_geom_c
#undef w_geom
#undef w_mult
#undef l_range
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'PortOrig' [3]. */
#define mccompcurname  PortOrig
#define mccompcurtype  Arm
#define mccompcurindex 3
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'XmonT' [4]. */
#define mccompcurname  XmonT
#define mccompcurtype  PSDlin_monitor
#define mccompcurindex 4
#define nx mccXmonT_nx
#define PSDlin_N mccXmonT_PSDlin_N
#define PSDlin_p mccXmonT_PSDlin_p
#define PSDlin_p2 mccXmonT_PSDlin_p2
#define filename mccXmonT_filename
#define xmin mccXmonT_xmin
#define xmax mccXmonT_xmax
#define ymin mccXmonT_ymin
#define ymax mccXmonT_ymax
#define xwidth mccXmonT_xwidth
#define yheight mccXmonT_yheight
#define restore_neutron mccXmonT_restore_neutron
#define nowritefile mccXmonT_nowritefile
#line 53 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/PSDlin_monitor.comp"
    double PSDlin_N[nx];
    double PSDlin_p[nx];
    double PSDlin_p2[nx];
#line 6778 "./ESS_2015_test.c"
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

/* User declarations for component 'XmonC' [5]. */
#define mccompcurname  XmonC
#define mccompcurtype  PSDlin_monitor
#define mccompcurindex 5
#define nx mccXmonC_nx
#define PSDlin_N mccXmonC_PSDlin_N
#define PSDlin_p mccXmonC_PSDlin_p
#define PSDlin_p2 mccXmonC_PSDlin_p2
#define filename mccXmonC_filename
#define xmin mccXmonC_xmin
#define xmax mccXmonC_xmax
#define ymin mccXmonC_ymin
#define ymax mccXmonC_ymax
#define xwidth mccXmonC_xwidth
#define yheight mccXmonC_yheight
#define restore_neutron mccXmonC_restore_neutron
#define nowritefile mccXmonC_nowritefile
#line 53 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/PSDlin_monitor.comp"
    double PSDlin_N[nx];
    double PSDlin_p[nx];
    double PSDlin_p2[nx];
#line 6817 "./ESS_2015_test.c"
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

/* User declarations for component 'tofT' [6]. */
#define mccompcurname  tofT
#define mccompcurtype  TOF_monitor
#define mccompcurindex 6
#define nt mcctofT_nt
#define TOF_N mcctofT_TOF_N
#define TOF_p mcctofT_TOF_p
#define TOF_p2 mcctofT_TOF_p2
#define t_min mcctofT_t_min
#define t_max mcctofT_t_max
#define delta_t mcctofT_delta_t
#define filename mcctofT_filename
#define xmin mcctofT_xmin
#define xmax mcctofT_xmax
#define ymin mcctofT_ymin
#define ymax mcctofT_ymax
#define xwidth mcctofT_xwidth
#define yheight mcctofT_yheight
#define tmin mcctofT_tmin
#define tmax mcctofT_tmax
#define dt mcctofT_dt
#define restore_neutron mcctofT_restore_neutron
#define nowritefile mcctofT_nowritefile
#line 54 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/TOF_monitor.comp"
double TOF_N[nt];
double TOF_p[nt];
double TOF_p2[nt];
double t_min, t_max, delta_t;
#line 6863 "./ESS_2015_test.c"
#undef nowritefile
#undef restore_neutron
#undef dt
#undef tmax
#undef tmin
#undef yheight
#undef xwidth
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef filename
#undef delta_t
#undef t_max
#undef t_min
#undef TOF_p2
#undef TOF_p
#undef TOF_N
#undef nt
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'tofC' [7]. */
#define mccompcurname  tofC
#define mccompcurtype  TOF_monitor
#define mccompcurindex 7
#define nt mcctofC_nt
#define TOF_N mcctofC_TOF_N
#define TOF_p mcctofC_TOF_p
#define TOF_p2 mcctofC_TOF_p2
#define t_min mcctofC_t_min
#define t_max mcctofC_t_max
#define delta_t mcctofC_delta_t
#define filename mcctofC_filename
#define xmin mcctofC_xmin
#define xmax mcctofC_xmax
#define ymin mcctofC_ymin
#define ymax mcctofC_ymax
#define xwidth mcctofC_xwidth
#define yheight mcctofC_yheight
#define tmin mcctofC_tmin
#define tmax mcctofC_tmax
#define dt mcctofC_dt
#define restore_neutron mcctofC_restore_neutron
#define nowritefile mcctofC_nowritefile
#line 54 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/TOF_monitor.comp"
double TOF_N[nt];
double TOF_p[nt];
double TOF_p2[nt];
double t_min, t_max, delta_t;
#line 6915 "./ESS_2015_test.c"
#undef nowritefile
#undef restore_neutron
#undef dt
#undef tmax
#undef tmin
#undef yheight
#undef xwidth
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef filename
#undef delta_t
#undef t_max
#undef t_min
#undef TOF_p2
#undef TOF_p
#undef TOF_N
#undef nt
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'Slit1' [8]. */
#define mccompcurname  Slit1
#define mccompcurtype  Slit
#define mccompcurindex 8
#define xmin mccSlit1_xmin
#define xmax mccSlit1_xmax
#define ymin mccSlit1_ymin
#define ymax mccSlit1_ymax
#define radius mccSlit1_radius
#define xwidth mccSlit1_xwidth
#define yheight mccSlit1_yheight
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

/* User declarations for component 'LamT' [9]. */
#define mccompcurname  LamT
#define mccompcurtype  L_monitor
#define mccompcurindex 9
#define nL mccLamT_nL
#define L_N mccLamT_L_N
#define L_p mccLamT_L_p
#define L_p2 mccLamT_L_p2
#define filename mccLamT_filename
#define xmin mccLamT_xmin
#define xmax mccLamT_xmax
#define ymin mccLamT_ymin
#define ymax mccLamT_ymax
#define xwidth mccLamT_xwidth
#define yheight mccLamT_yheight
#define Lmin mccLamT_Lmin
#define Lmax mccLamT_Lmax
#define restore_neutron mccLamT_restore_neutron
#define nowritefile mccLamT_nowritefile
#line 57 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/L_monitor.comp"
double L_N[nL];
double L_p[nL], L_p2[nL];
#line 6983 "./ESS_2015_test.c"
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

/* User declarations for component 'LamC' [10]. */
#define mccompcurname  LamC
#define mccompcurtype  L_monitor
#define mccompcurindex 10
#define nL mccLamC_nL
#define L_N mccLamC_L_N
#define L_p mccLamC_L_p
#define L_p2 mccLamC_L_p2
#define filename mccLamC_filename
#define xmin mccLamC_xmin
#define xmax mccLamC_xmax
#define ymin mccLamC_ymin
#define ymax mccLamC_ymax
#define xwidth mccLamC_xwidth
#define yheight mccLamC_yheight
#define Lmin mccLamC_Lmin
#define Lmax mccLamC_Lmax
#define restore_neutron mccLamC_restore_neutron
#define nowritefile mccLamC_nowritefile
#line 57 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/L_monitor.comp"
double L_N[nL];
double L_p[nL], L_p2[nL];
#line 7025 "./ESS_2015_test.c"
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

Coords mcposaOrigin, mcposrOrigin;
Rotation mcrotaOrigin, mcrotrOrigin;
Coords mcposaSource, mcposrSource;
Rotation mcrotaSource, mcrotrSource;
Coords mcposaPortOrig, mcposrPortOrig;
Rotation mcrotaPortOrig, mcrotrPortOrig;
Coords mcposaXmonT, mcposrXmonT;
Rotation mcrotaXmonT, mcrotrXmonT;
Coords mcposaXmonC, mcposrXmonC;
Rotation mcrotaXmonC, mcrotrXmonC;
Coords mcposatofT, mcposrtofT;
Rotation mcrotatofT, mcrotrtofT;
Coords mcposatofC, mcposrtofC;
Rotation mcrotatofC, mcrotrtofC;
Coords mcposaSlit1, mcposrSlit1;
Rotation mcrotaSlit1, mcrotrSlit1;
Coords mcposaLamT, mcposrLamT;
Rotation mcrotaLamT, mcrotrLamT;
Coords mcposaLamC, mcposrLamC;
Rotation mcrotaLamC, mcrotrLamC;

MCNUM mcnx, mcny, mcnz, mcnvx, mcnvy, mcnvz, mcnt, mcnsx, mcnsy, mcnsz, mcnp;

/* end declare */

void mcinit(void) {
#define mccompcurname  ESS_2015_test
#define mccompcurtype  INSTRUMENT
#define mccompcurindex 0
#define mcposaESS_2015_test coords_set(0,0,0)
#define frac mcipfrac
#define lmin mciplmin
#define lmax mciplmax
#define ANGLE mcipANGLE
#define mon_shift mcipmon_shift
#line 52 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
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
#line 7100 "./ESS_2015_test.c"
#undef mon_shift
#undef ANGLE
#undef lmax
#undef lmin
#undef frac
#undef mcposaESS_2015_test
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
#line 39 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
  if("NULL") strncpy(mccOrigin_profile, "NULL" ? "NULL" : "", 16384); else mccOrigin_profile[0]='\0';
#line 39 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
  mccOrigin_percent = 10;
#line 39 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
  mccOrigin_flag_save = 0;
#line 39 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
  mccOrigin_minutes = 0;
#line 7131 "./ESS_2015_test.c"

  SIG_MESSAGE("Origin (Init:Place/Rotate)");
  rot_set_rotation(mcrotaOrigin,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 7138 "./ESS_2015_test.c"
  rot_copy(mcrotrOrigin, mcrotaOrigin);
  mcposaOrigin = coords_set(
#line 84 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
    0,
#line 84 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
    0,
#line 84 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
    0);
#line 7147 "./ESS_2015_test.c"
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
#line 81 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
  mccSource_isleft = 0.9;
#line 87 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
  mccSource_Lmin = mciplmin;
#line 87 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
  mccSource_Lmax = mciplmax;
#line 87 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
  mccSource_cold_frac = mcipfrac;
#line 82 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
  mccSource_dist = 0;
#line 86 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
  mccSource_focus_xw = 0.01;
#line 86 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
  mccSource_focus_yh = 0.01;
#line 86 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
  mccSource_target_index = + 6;
#line 82 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
  mccSource_tmax_multiplier = 3;
#line 86 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
  mccSource_yheight_c = Yheight;
#line 86 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
  mccSource_yheight_t = Yheight;
#line 87 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
  mccSource_n_pulses = np;
#line 87 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
  mccSource_acc_power = power;
#line 87 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
  mccSource_beamport_angle = mcipANGLE;
#line 87 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
  if(srcdef) strncpy(mccSource_sourcedef, srcdef ? srcdef : "", 16384); else mccSource_sourcedef[0]='\0';
#line 83 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
  mccSource_xwidth_c = 0.1;
#line 83 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
  mccSource_xwidth_t = 0.18;
#line 83 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
  mccSource_extraction_opening = 120;
#line 7194 "./ESS_2015_test.c"

  SIG_MESSAGE("Source (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 7201 "./ESS_2015_test.c"
  rot_mul(mctr1, mcrotaOrigin, mcrotaSource);
  rot_transpose(mcrotaOrigin, mctr1);
  rot_mul(mcrotaSource, mctr1, mcrotrSource);
  mctc1 = coords_set(
#line 88 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
    0,
#line 88 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
    0,
#line 88 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
    0);
#line 7212 "./ESS_2015_test.c"
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
    /* Component PortOrig. */
  /* Setting parameters for component PortOrig. */
  SIG_MESSAGE("PortOrig (Init:SetPar)");

  SIG_MESSAGE("PortOrig (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 95 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
    (0)*DEG2RAD,
#line 95 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
    (0)*DEG2RAD,
#line 95 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
    (0)*DEG2RAD);
#line 7235 "./ESS_2015_test.c"
  rot_mul(mctr1, mcrotaOrigin, mcrotaPortOrig);
  rot_transpose(mcrotaSource, mctr1);
  rot_mul(mcrotaPortOrig, mctr1, mcrotrPortOrig);
  mctc1 = coords_set(
#line 94 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
    BX0,
#line 94 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
    0,
#line 94 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
    BZ0);
#line 7246 "./ESS_2015_test.c"
  rot_transpose(mcrotaOrigin, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaPortOrig = coords_add(mcposaOrigin, mctc2);
  mctc1 = coords_sub(mcposaSource, mcposaPortOrig);
  mcposrPortOrig = rot_apply(mcrotaPortOrig, mctc1);
  mcDEBUG_COMPONENT("PortOrig", mcposaPortOrig, mcrotaPortOrig)
  mccomp_posa[3] = mcposaPortOrig;
  mccomp_posr[3] = mcposrPortOrig;
  mcNCounter[3]  = mcPCounter[3] = mcP2Counter[3] = 0;
  mcAbsorbProp[3]= 0;
    /* Component XmonT. */
  /* Setting parameters for component XmonT. */
  SIG_MESSAGE("XmonT (Init:SetPar)");
#line 98 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
  if("XmonT.dat") strncpy(mccXmonT_filename, "XmonT.dat" ? "XmonT.dat" : "", 16384); else mccXmonT_filename[0]='\0';
#line 98 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
  mccXmonT_xmin = - mon_width / 2 + mcipmon_shift;
#line 98 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
  mccXmonT_xmax = mon_width / 2 + mcipmon_shift;
#line 46 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
  mccXmonT_ymin = -0.05;
#line 46 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
  mccXmonT_ymax = 0.05;
#line 47 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
  mccXmonT_xwidth = 0;
#line 98 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
  mccXmonT_yheight = 0.04;
#line 99 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
  mccXmonT_restore_neutron = 1;
#line 47 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
  mccXmonT_nowritefile = 0;
#line 7278 "./ESS_2015_test.c"

  SIG_MESSAGE("XmonT (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 7285 "./ESS_2015_test.c"
  rot_mul(mctr1, mcrotaPortOrig, mcrotaXmonT);
  rot_transpose(mcrotaPortOrig, mctr1);
  rot_mul(mcrotaXmonT, mctr1, mcrotrXmonT);
  mctc1 = coords_set(
#line 100 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
    0,
#line 100 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
    0,
#line 100 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
    0.3);
#line 7296 "./ESS_2015_test.c"
  rot_transpose(mcrotaPortOrig, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaXmonT = coords_add(mcposaPortOrig, mctc2);
  mctc1 = coords_sub(mcposaPortOrig, mcposaXmonT);
  mcposrXmonT = rot_apply(mcrotaXmonT, mctc1);
  mcDEBUG_COMPONENT("XmonT", mcposaXmonT, mcrotaXmonT)
  mccomp_posa[4] = mcposaXmonT;
  mccomp_posr[4] = mcposrXmonT;
  mcNCounter[4]  = mcPCounter[4] = mcP2Counter[4] = 0;
  mcAbsorbProp[4]= 0;
    /* Component XmonC. */
  /* Setting parameters for component XmonC. */
  SIG_MESSAGE("XmonC (Init:SetPar)");
#line 103 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
  if("XmonC.dat") strncpy(mccXmonC_filename, "XmonC.dat" ? "XmonC.dat" : "", 16384); else mccXmonC_filename[0]='\0';
#line 103 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
  mccXmonC_xmin = - mon_width / 2 + mcipmon_shift;
#line 103 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
  mccXmonC_xmax = mon_width / 2 + mcipmon_shift;
#line 46 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
  mccXmonC_ymin = -0.05;
#line 46 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
  mccXmonC_ymax = 0.05;
#line 47 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
  mccXmonC_xwidth = 0;
#line 103 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
  mccXmonC_yheight = 0.04;
#line 104 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
  mccXmonC_restore_neutron = 1;
#line 47 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
  mccXmonC_nowritefile = 0;
#line 7328 "./ESS_2015_test.c"

  SIG_MESSAGE("XmonC (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 7335 "./ESS_2015_test.c"
  rot_mul(mctr1, mcrotaPortOrig, mcrotaXmonC);
  rot_transpose(mcrotaXmonT, mctr1);
  rot_mul(mcrotaXmonC, mctr1, mcrotrXmonC);
  mctc1 = coords_set(
#line 105 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
    0,
#line 105 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
    0,
#line 105 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
    0.3);
#line 7346 "./ESS_2015_test.c"
  rot_transpose(mcrotaPortOrig, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaXmonC = coords_add(mcposaPortOrig, mctc2);
  mctc1 = coords_sub(mcposaXmonT, mcposaXmonC);
  mcposrXmonC = rot_apply(mcrotaXmonC, mctc1);
  mcDEBUG_COMPONENT("XmonC", mcposaXmonC, mcrotaXmonC)
  mccomp_posa[5] = mcposaXmonC;
  mccomp_posr[5] = mcposrXmonC;
  mcNCounter[5]  = mcPCounter[5] = mcP2Counter[5] = 0;
  mcAbsorbProp[5]= 0;
    /* Component tofT. */
  /* Setting parameters for component tofT. */
  SIG_MESSAGE("tofT (Init:SetPar)");
#line 108 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
  if("tofT.dat") strncpy(mcctofT_filename, "tofT.dat" ? "tofT.dat" : "", 16384); else mcctofT_filename[0]='\0';
#line 108 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
  mcctofT_xmin = - mon_width / 2 + mcipmon_shift;
#line 108 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
  mcctofT_xmax = mon_width / 2 + mcipmon_shift;
#line 47 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
  mcctofT_ymin = -0.05;
#line 47 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
  mcctofT_ymax = 0.05;
#line 48 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
  mcctofT_xwidth = 0;
#line 108 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
  mcctofT_yheight = 0.04;
#line 109 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
  mcctofT_tmin = 0;
#line 109 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
  mcctofT_tmax = 10000;
#line 48 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
  mcctofT_dt = 1.0;
#line 109 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
  mcctofT_restore_neutron = 1;
#line 48 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
  mcctofT_nowritefile = 0;
#line 7384 "./ESS_2015_test.c"

  SIG_MESSAGE("tofT (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 7391 "./ESS_2015_test.c"
  rot_mul(mctr1, mcrotaPortOrig, mcrotatofT);
  rot_transpose(mcrotaXmonC, mctr1);
  rot_mul(mcrotatofT, mctr1, mcrotrtofT);
  mctc1 = coords_set(
#line 110 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
    0,
#line 110 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
    0,
#line 110 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
    0.3);
#line 7402 "./ESS_2015_test.c"
  rot_transpose(mcrotaPortOrig, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposatofT = coords_add(mcposaPortOrig, mctc2);
  mctc1 = coords_sub(mcposaXmonC, mcposatofT);
  mcposrtofT = rot_apply(mcrotatofT, mctc1);
  mcDEBUG_COMPONENT("tofT", mcposatofT, mcrotatofT)
  mccomp_posa[6] = mcposatofT;
  mccomp_posr[6] = mcposrtofT;
  mcNCounter[6]  = mcPCounter[6] = mcP2Counter[6] = 0;
  mcAbsorbProp[6]= 0;
    /* Component tofC. */
  /* Setting parameters for component tofC. */
  SIG_MESSAGE("tofC (Init:SetPar)");
#line 113 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
  if("tofC.dat") strncpy(mcctofC_filename, "tofC.dat" ? "tofC.dat" : "", 16384); else mcctofC_filename[0]='\0';
#line 113 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
  mcctofC_xmin = - mon_width / 2 + mcipmon_shift;
#line 113 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
  mcctofC_xmax = mon_width / 2 + mcipmon_shift;
#line 47 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
  mcctofC_ymin = -0.05;
#line 47 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
  mcctofC_ymax = 0.05;
#line 48 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
  mcctofC_xwidth = 0;
#line 113 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
  mcctofC_yheight = 0.04;
#line 114 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
  mcctofC_tmin = 0;
#line 114 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
  mcctofC_tmax = 10000;
#line 48 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
  mcctofC_dt = 1.0;
#line 114 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
  mcctofC_restore_neutron = 1;
#line 48 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
  mcctofC_nowritefile = 0;
#line 7440 "./ESS_2015_test.c"

  SIG_MESSAGE("tofC (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 7447 "./ESS_2015_test.c"
  rot_mul(mctr1, mcrotaPortOrig, mcrotatofC);
  rot_transpose(mcrotatofT, mctr1);
  rot_mul(mcrotatofC, mctr1, mcrotrtofC);
  mctc1 = coords_set(
#line 115 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
    0,
#line 115 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
    0,
#line 115 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
    0.3);
#line 7458 "./ESS_2015_test.c"
  rot_transpose(mcrotaPortOrig, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposatofC = coords_add(mcposaPortOrig, mctc2);
  mctc1 = coords_sub(mcposatofT, mcposatofC);
  mcposrtofC = rot_apply(mcrotatofC, mctc1);
  mcDEBUG_COMPONENT("tofC", mcposatofC, mcrotatofC)
  mccomp_posa[7] = mcposatofC;
  mccomp_posr[7] = mcposrtofC;
  mcNCounter[7]  = mcPCounter[7] = mcP2Counter[7] = 0;
  mcAbsorbProp[7]= 0;
    /* Component Slit1. */
  /* Setting parameters for component Slit1. */
  SIG_MESSAGE("Slit1 (Init:SetPar)");
#line 46 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
  mccSlit1_xmin = 0;
#line 46 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
  mccSlit1_xmax = 0;
#line 46 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
  mccSlit1_ymin = 0;
#line 46 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
  mccSlit1_ymax = 0;
#line 46 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
  mccSlit1_radius = 0;
#line 118 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
  mccSlit1_xwidth = 0.01;
#line 118 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
  mccSlit1_yheight = 0.01;
#line 7486 "./ESS_2015_test.c"

  SIG_MESSAGE("Slit1 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 7493 "./ESS_2015_test.c"
  rot_mul(mctr1, mcrotaPortOrig, mcrotaSlit1);
  rot_transpose(mcrotatofC, mctr1);
  rot_mul(mcrotaSlit1, mctr1, mcrotrSlit1);
  mctc1 = coords_set(
#line 119 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
    0,
#line 119 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
    0,
#line 119 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
    100);
#line 7504 "./ESS_2015_test.c"
  rot_transpose(mcrotaPortOrig, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaSlit1 = coords_add(mcposaPortOrig, mctc2);
  mctc1 = coords_sub(mcposatofC, mcposaSlit1);
  mcposrSlit1 = rot_apply(mcrotaSlit1, mctc1);
  mcDEBUG_COMPONENT("Slit1", mcposaSlit1, mcrotaSlit1)
  mccomp_posa[8] = mcposaSlit1;
  mccomp_posr[8] = mcposrSlit1;
  mcNCounter[8]  = mcPCounter[8] = mcP2Counter[8] = 0;
  mcAbsorbProp[8]= 0;
    /* Component LamT. */
  /* Setting parameters for component LamT. */
  SIG_MESSAGE("LamT (Init:SetPar)");
#line 122 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
  if("LamT.dat") strncpy(mccLamT_filename, "LamT.dat" ? "LamT.dat" : "", 16384); else mccLamT_filename[0]='\0';
#line 50 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
  mccLamT_xmin = -0.05;
#line 50 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
  mccLamT_xmax = 0.05;
#line 50 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
  mccLamT_ymin = -0.05;
#line 50 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
  mccLamT_ymax = 0.05;
#line 122 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
  mccLamT_xwidth = 0.02;
#line 122 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
  mccLamT_yheight = 0.02;
#line 123 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
  mccLamT_Lmin = mciplmin;
#line 123 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
  mccLamT_Lmax = mciplmax;
#line 123 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
  mccLamT_restore_neutron = 1;
#line 51 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
  mccLamT_nowritefile = 0;
#line 7540 "./ESS_2015_test.c"

  SIG_MESSAGE("LamT (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 7547 "./ESS_2015_test.c"
  rot_mul(mctr1, mcrotaSlit1, mcrotaLamT);
  rot_transpose(mcrotaSlit1, mctr1);
  rot_mul(mcrotaLamT, mctr1, mcrotrLamT);
  mctc1 = coords_set(
#line 124 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
    0,
#line 124 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
    0,
#line 124 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
    0.0001);
#line 7558 "./ESS_2015_test.c"
  rot_transpose(mcrotaSlit1, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaLamT = coords_add(mcposaSlit1, mctc2);
  mctc1 = coords_sub(mcposaSlit1, mcposaLamT);
  mcposrLamT = rot_apply(mcrotaLamT, mctc1);
  mcDEBUG_COMPONENT("LamT", mcposaLamT, mcrotaLamT)
  mccomp_posa[9] = mcposaLamT;
  mccomp_posr[9] = mcposrLamT;
  mcNCounter[9]  = mcPCounter[9] = mcP2Counter[9] = 0;
  mcAbsorbProp[9]= 0;
    /* Component LamC. */
  /* Setting parameters for component LamC. */
  SIG_MESSAGE("LamC (Init:SetPar)");
#line 128 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
  if("LamC.dat") strncpy(mccLamC_filename, "LamC.dat" ? "LamC.dat" : "", 16384); else mccLamC_filename[0]='\0';
#line 50 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
  mccLamC_xmin = -0.05;
#line 50 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
  mccLamC_xmax = 0.05;
#line 50 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
  mccLamC_ymin = -0.05;
#line 50 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
  mccLamC_ymax = 0.05;
#line 128 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
  mccLamC_xwidth = 0.02;
#line 128 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
  mccLamC_yheight = 0.02;
#line 129 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
  mccLamC_Lmin = mciplmin;
#line 129 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
  mccLamC_Lmax = mciplmax;
#line 129 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
  mccLamC_restore_neutron = 1;
#line 51 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
  mccLamC_nowritefile = 0;
#line 7594 "./ESS_2015_test.c"

  SIG_MESSAGE("LamC (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 7601 "./ESS_2015_test.c"
  rot_mul(mctr1, mcrotaSlit1, mcrotaLamC);
  rot_transpose(mcrotaLamT, mctr1);
  rot_mul(mcrotaLamC, mctr1, mcrotrLamC);
  mctc1 = coords_set(
#line 130 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
    0,
#line 130 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
    0,
#line 130 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
    0.0001);
#line 7612 "./ESS_2015_test.c"
  rot_transpose(mcrotaSlit1, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaLamC = coords_add(mcposaSlit1, mctc2);
  mctc1 = coords_sub(mcposaLamT, mcposaLamC);
  mcposrLamC = rot_apply(mcrotaLamC, mctc1);
  mcDEBUG_COMPONENT("LamC", mcposaLamC, mcrotaLamC)
  mccomp_posa[10] = mcposaLamC;
  mccomp_posr[10] = mcposrLamC;
  mcNCounter[10]  = mcPCounter[10] = mcP2Counter[10] = 0;
  mcAbsorbProp[10]= 0;
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
#line 7649 "./ESS_2015_test.c"
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
#define mccompcurtype  ESS_moderator
#define mccompcurindex 2
#define l_range mccSource_l_range
#define w_mult mccSource_w_mult
#define w_geom mccSource_w_geom
#define w_geom_c mccSource_w_geom_c
#define w_geom_t mccSource_w_geom_t
#define w_stat mccSource_w_stat
#define cold mccSource_cold
#define wasleft mccSource_wasleft
#define cosine mccSource_cosine
#define tx mccSource_tx
#define ty mccSource_ty
#define tz mccSource_tz
#define cx mccSource_cx
#define cy mccSource_cy
#define cz mccSource_cz
#define cyl_radius mccSource_cyl_radius
#define t1x mccSource_t1x
#define t1y mccSource_t1y
#define t1z mccSource_t1z
#define t2x mccSource_t2x
#define t2y mccSource_t2y
#define t2z mccSource_t2z
#define modextras mccSource_modextras
#define cold_bril mccSource_cold_bril
#define thermal_bril mccSource_thermal_bril
#define tfocus_width mccSource_tfocus_width
#define tfocus_time mccSource_tfocus_time
#define dt mccSource_dt
#define r_empty mccSource_r_empty
#define lambda mccSource_lambda
#define xprime mccSource_xprime
#define yprime mccSource_yprime
#define zprime mccSource_zprime
#define sgnport mccSource_sgnport
#define abs_beamport_angle mccSource_abs_beamport_angle
#define cos_beamport_angle mccSource_cos_beamport_angle
#define edge_thermal mccSource_edge_thermal
#define cos_thermal mccSource_cos_thermal
#define cos_cold mccSource_cos_cold
#define sin_beamport_angle mccSource_sin_beamport_angle
#define isleft mccSource_isleft
#define Lmin mccSource_Lmin
#define Lmax mccSource_Lmax
#define cold_frac mccSource_cold_frac
#define dist mccSource_dist
#define focus_xw mccSource_focus_xw
#define focus_yh mccSource_focus_yh
#define target_index mccSource_target_index
#define tmax_multiplier mccSource_tmax_multiplier
#define yheight_c mccSource_yheight_c
#define yheight_t mccSource_yheight_t
#define n_pulses mccSource_n_pulses
#define acc_power mccSource_acc_power
#define beamport_angle mccSource_beamport_angle
#define sourcedef mccSource_sourcedef
#define xwidth_c mccSource_xwidth_c
#define xwidth_t mccSource_xwidth_t
#define extraction_opening mccSource_extraction_opening
#line 128 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../sources/ESS_moderator.comp"
{

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

}
#line 7892 "./ESS_2015_test.c"
#undef extraction_opening
#undef xwidth_t
#undef xwidth_c
#undef sourcedef
#undef beamport_angle
#undef acc_power
#undef n_pulses
#undef yheight_t
#undef yheight_c
#undef tmax_multiplier
#undef target_index
#undef focus_yh
#undef focus_xw
#undef dist
#undef cold_frac
#undef Lmax
#undef Lmin
#undef isleft
#undef sin_beamport_angle
#undef cos_cold
#undef cos_thermal
#undef edge_thermal
#undef cos_beamport_angle
#undef abs_beamport_angle
#undef sgnport
#undef zprime
#undef yprime
#undef xprime
#undef lambda
#undef r_empty
#undef dt
#undef tfocus_time
#undef tfocus_width
#undef thermal_bril
#undef cold_bril
#undef modextras
#undef t2z
#undef t2y
#undef t2x
#undef t1z
#undef t1y
#undef t1x
#undef cyl_radius
#undef cz
#undef cy
#undef cx
#undef tz
#undef ty
#undef tx
#undef cosine
#undef wasleft
#undef cold
#undef w_stat
#undef w_geom_t
#undef w_geom_c
#undef w_geom
#undef w_mult
#undef l_range
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component PortOrig. */
  SIG_MESSAGE("PortOrig (Init)");

  /* Initializations for component XmonT. */
  SIG_MESSAGE("XmonT (Init)");
#define mccompcurname  XmonT
#define mccompcurtype  PSDlin_monitor
#define mccompcurindex 4
#define nx mccXmonT_nx
#define PSDlin_N mccXmonT_PSDlin_N
#define PSDlin_p mccXmonT_PSDlin_p
#define PSDlin_p2 mccXmonT_PSDlin_p2
#define filename mccXmonT_filename
#define xmin mccXmonT_xmin
#define xmax mccXmonT_xmax
#define ymin mccXmonT_ymin
#define ymax mccXmonT_ymax
#define xwidth mccXmonT_xwidth
#define yheight mccXmonT_yheight
#define restore_neutron mccXmonT_restore_neutron
#define nowritefile mccXmonT_nowritefile
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
#line 7997 "./ESS_2015_test.c"
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

  /* Initializations for component XmonC. */
  SIG_MESSAGE("XmonC (Init)");
#define mccompcurname  XmonC
#define mccompcurtype  PSDlin_monitor
#define mccompcurindex 5
#define nx mccXmonC_nx
#define PSDlin_N mccXmonC_PSDlin_N
#define PSDlin_p mccXmonC_PSDlin_p
#define PSDlin_p2 mccXmonC_PSDlin_p2
#define filename mccXmonC_filename
#define xmin mccXmonC_xmin
#define xmax mccXmonC_xmax
#define ymin mccXmonC_ymin
#define ymax mccXmonC_ymax
#define xwidth mccXmonC_xwidth
#define yheight mccXmonC_yheight
#define restore_neutron mccXmonC_restore_neutron
#define nowritefile mccXmonC_nowritefile
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
#line 8054 "./ESS_2015_test.c"
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

  /* Initializations for component tofT. */
  SIG_MESSAGE("tofT (Init)");
#define mccompcurname  tofT
#define mccompcurtype  TOF_monitor
#define mccompcurindex 6
#define nt mcctofT_nt
#define TOF_N mcctofT_TOF_N
#define TOF_p mcctofT_TOF_p
#define TOF_p2 mcctofT_TOF_p2
#define t_min mcctofT_t_min
#define t_max mcctofT_t_max
#define delta_t mcctofT_delta_t
#define filename mcctofT_filename
#define xmin mcctofT_xmin
#define xmax mcctofT_xmax
#define ymin mcctofT_ymin
#define ymax mcctofT_ymax
#define xwidth mcctofT_xwidth
#define yheight mcctofT_yheight
#define tmin mcctofT_tmin
#define tmax mcctofT_tmax
#define dt mcctofT_dt
#define restore_neutron mcctofT_restore_neutron
#define nowritefile mcctofT_nowritefile
#line 61 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/TOF_monitor.comp"
{
int i;

if (xwidth  > 0) { xmax = xwidth/2;  xmin = -xmax; }
    if (yheight > 0) { ymax = yheight/2; ymin = -ymax; }

    if ((xmin >= xmax) || (ymin >= ymax)) {
            printf("TOF_monitor: %s: Null detection area !\n"
                   "ERROR        (xwidth,yheight,xmin,xmax,ymin,ymax). Exiting",
           NAME_CURRENT_COMP);
      exit(0);
    }

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
}
#line 8129 "./ESS_2015_test.c"
#undef nowritefile
#undef restore_neutron
#undef dt
#undef tmax
#undef tmin
#undef yheight
#undef xwidth
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef filename
#undef delta_t
#undef t_max
#undef t_min
#undef TOF_p2
#undef TOF_p
#undef TOF_N
#undef nt
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component tofC. */
  SIG_MESSAGE("tofC (Init)");
#define mccompcurname  tofC
#define mccompcurtype  TOF_monitor
#define mccompcurindex 7
#define nt mcctofC_nt
#define TOF_N mcctofC_TOF_N
#define TOF_p mcctofC_TOF_p
#define TOF_p2 mcctofC_TOF_p2
#define t_min mcctofC_t_min
#define t_max mcctofC_t_max
#define delta_t mcctofC_delta_t
#define filename mcctofC_filename
#define xmin mcctofC_xmin
#define xmax mcctofC_xmax
#define ymin mcctofC_ymin
#define ymax mcctofC_ymax
#define xwidth mcctofC_xwidth
#define yheight mcctofC_yheight
#define tmin mcctofC_tmin
#define tmax mcctofC_tmax
#define dt mcctofC_dt
#define restore_neutron mcctofC_restore_neutron
#define nowritefile mcctofC_nowritefile
#line 61 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/TOF_monitor.comp"
{
int i;

if (xwidth  > 0) { xmax = xwidth/2;  xmin = -xmax; }
    if (yheight > 0) { ymax = yheight/2; ymin = -ymax; }

    if ((xmin >= xmax) || (ymin >= ymax)) {
            printf("TOF_monitor: %s: Null detection area !\n"
                   "ERROR        (xwidth,yheight,xmin,xmax,ymin,ymax). Exiting",
           NAME_CURRENT_COMP);
      exit(0);
    }

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
}
#line 8210 "./ESS_2015_test.c"
#undef nowritefile
#undef restore_neutron
#undef dt
#undef tmax
#undef tmin
#undef yheight
#undef xwidth
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef filename
#undef delta_t
#undef t_max
#undef t_min
#undef TOF_p2
#undef TOF_p
#undef TOF_N
#undef nt
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component Slit1. */
  SIG_MESSAGE("Slit1 (Init)");
#define mccompcurname  Slit1
#define mccompcurtype  Slit
#define mccompcurindex 8
#define xmin mccSlit1_xmin
#define xmax mccSlit1_xmax
#define ymin mccSlit1_ymin
#define ymax mccSlit1_ymax
#define radius mccSlit1_radius
#define xwidth mccSlit1_xwidth
#define yheight mccSlit1_yheight
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
#line 8266 "./ESS_2015_test.c"
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

  /* Initializations for component LamT. */
  SIG_MESSAGE("LamT (Init)");
#define mccompcurname  LamT
#define mccompcurtype  L_monitor
#define mccompcurindex 9
#define nL mccLamT_nL
#define L_N mccLamT_L_N
#define L_p mccLamT_L_p
#define L_p2 mccLamT_L_p2
#define filename mccLamT_filename
#define xmin mccLamT_xmin
#define xmax mccLamT_xmax
#define ymin mccLamT_ymin
#define ymax mccLamT_ymax
#define xwidth mccLamT_xwidth
#define yheight mccLamT_yheight
#define Lmin mccLamT_Lmin
#define Lmax mccLamT_Lmax
#define restore_neutron mccLamT_restore_neutron
#define nowritefile mccLamT_nowritefile
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
#line 8319 "./ESS_2015_test.c"
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

  /* Initializations for component LamC. */
  SIG_MESSAGE("LamC (Init)");
#define mccompcurname  LamC
#define mccompcurtype  L_monitor
#define mccompcurindex 10
#define nL mccLamC_nL
#define L_N mccLamC_L_N
#define L_p mccLamC_L_p
#define L_p2 mccLamC_L_p2
#define filename mccLamC_filename
#define xmin mccLamC_xmin
#define xmax mccLamC_xmax
#define ymin mccLamC_ymin
#define ymax mccLamC_ymax
#define xwidth mccLamC_xwidth
#define yheight mccLamC_yheight
#define Lmin mccLamC_Lmin
#define Lmax mccLamC_Lmax
#define restore_neutron mccLamC_restore_neutron
#define nowritefile mccLamC_nowritefile
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
#line 8380 "./ESS_2015_test.c"
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
#line 8553 "./ESS_2015_test.c"
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
#define mccompcurtype  ESS_moderator
#define mccompcurindex 2
#define l_range mccSource_l_range
#define w_mult mccSource_w_mult
#define w_geom mccSource_w_geom
#define w_geom_c mccSource_w_geom_c
#define w_geom_t mccSource_w_geom_t
#define w_stat mccSource_w_stat
#define cold mccSource_cold
#define wasleft mccSource_wasleft
#define cosine mccSource_cosine
#define tx mccSource_tx
#define ty mccSource_ty
#define tz mccSource_tz
#define cx mccSource_cx
#define cy mccSource_cy
#define cz mccSource_cz
#define cyl_radius mccSource_cyl_radius
#define t1x mccSource_t1x
#define t1y mccSource_t1y
#define t1z mccSource_t1z
#define t2x mccSource_t2x
#define t2y mccSource_t2y
#define t2z mccSource_t2z
#define modextras mccSource_modextras
#define cold_bril mccSource_cold_bril
#define thermal_bril mccSource_thermal_bril
#define tfocus_width mccSource_tfocus_width
#define tfocus_time mccSource_tfocus_time
#define dt mccSource_dt
#define r_empty mccSource_r_empty
#define lambda mccSource_lambda
#define xprime mccSource_xprime
#define yprime mccSource_yprime
#define zprime mccSource_zprime
#define sgnport mccSource_sgnport
#define abs_beamport_angle mccSource_abs_beamport_angle
#define cos_beamport_angle mccSource_cos_beamport_angle
#define edge_thermal mccSource_edge_thermal
#define cos_thermal mccSource_cos_thermal
#define cos_cold mccSource_cos_cold
#define sin_beamport_angle mccSource_sin_beamport_angle
{   /* Declarations of Source=ESS_moderator() SETTING parameters. */
MCNUM isleft = mccSource_isleft;
MCNUM Lmin = mccSource_Lmin;
MCNUM Lmax = mccSource_Lmax;
MCNUM cold_frac = mccSource_cold_frac;
MCNUM dist = mccSource_dist;
MCNUM focus_xw = mccSource_focus_xw;
MCNUM focus_yh = mccSource_focus_yh;
int target_index = mccSource_target_index;
MCNUM tmax_multiplier = mccSource_tmax_multiplier;
MCNUM yheight_c = mccSource_yheight_c;
MCNUM yheight_t = mccSource_yheight_t;
int n_pulses = mccSource_n_pulses;
MCNUM acc_power = mccSource_acc_power;
MCNUM beamport_angle = mccSource_beamport_angle;
char* sourcedef = mccSource_sourcedef;
MCNUM xwidth_c = mccSource_xwidth_c;
MCNUM xwidth_t = mccSource_xwidth_t;
MCNUM extraction_opening = mccSource_extraction_opening;
#line 295 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../sources/ESS_moderator.comp"
{
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

}
#line 8865 "./ESS_2015_test.c"
/* 'Source=ESS_moderator()' component instance extend code */
    SIG_MESSAGE("Source (Trace:Extend)");
#line 90 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/ESS_2015_test/ESS_2015_test.instr"
  IsCold=cold;
#line 8870 "./ESS_2015_test.c"
}   /* End of Source=ESS_moderator() SETTING parameter declarations. */
#undef sin_beamport_angle
#undef cos_cold
#undef cos_thermal
#undef edge_thermal
#undef cos_beamport_angle
#undef abs_beamport_angle
#undef sgnport
#undef zprime
#undef yprime
#undef xprime
#undef lambda
#undef r_empty
#undef dt
#undef tfocus_time
#undef tfocus_width
#undef thermal_bril
#undef cold_bril
#undef modextras
#undef t2z
#undef t2y
#undef t2x
#undef t1z
#undef t1y
#undef t1x
#undef cyl_radius
#undef cz
#undef cy
#undef cx
#undef tz
#undef ty
#undef tx
#undef cosine
#undef wasleft
#undef cold
#undef w_stat
#undef w_geom_t
#undef w_geom_c
#undef w_geom
#undef w_mult
#undef l_range
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

  /* TRACE Component PortOrig [3] */
  mccoordschange(mcposrPortOrig, mcrotrPortOrig,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component PortOrig (without coords transformations) */
  mcJumpTrace_PortOrig:
  SIG_MESSAGE("PortOrig (Trace)");
  mcDEBUG_COMP("PortOrig")
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

#define mcabsorbComp mcabsorbCompPortOrig
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
#define mccompcurname  PortOrig
#define mccompcurtype  Arm
#define mccompcurindex 3
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompPortOrig:
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

  /* TRACE Component XmonT [4] */
  mccoordschange(mcposrXmonT, mcrotrXmonT,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component XmonT (without coords transformations) */
  mcJumpTrace_XmonT:
  SIG_MESSAGE("XmonT (Trace)");
  mcDEBUG_COMP("XmonT")
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

#define mcabsorbComp mcabsorbCompXmonT
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
#define mccompcurname  XmonT
#define mccompcurtype  PSDlin_monitor
#define mccompcurindex 4
#define nx mccXmonT_nx
#define PSDlin_N mccXmonT_PSDlin_N
#define PSDlin_p mccXmonT_PSDlin_p
#define PSDlin_p2 mccXmonT_PSDlin_p2
{   /* Declarations of XmonT=PSDlin_monitor() SETTING parameters. */
char* filename = mccXmonT_filename;
MCNUM xmin = mccXmonT_xmin;
MCNUM xmax = mccXmonT_xmax;
MCNUM ymin = mccXmonT_ymin;
MCNUM ymax = mccXmonT_ymax;
MCNUM xwidth = mccXmonT_xwidth;
MCNUM yheight = mccXmonT_yheight;
MCNUM restore_neutron = mccXmonT_restore_neutron;
int nowritefile = mccXmonT_nowritefile;
/* 'XmonT=PSDlin_monitor()' component instance has conditional execution */
if (( ! IsCold ))

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
#line 9155 "./ESS_2015_test.c"
}   /* End of XmonT=PSDlin_monitor() SETTING parameter declarations. */
#undef PSDlin_p2
#undef PSDlin_p
#undef PSDlin_N
#undef nx
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompXmonT:
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

  /* TRACE Component XmonC [5] */
  mccoordschange(mcposrXmonC, mcrotrXmonC,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component XmonC (without coords transformations) */
  mcJumpTrace_XmonC:
  SIG_MESSAGE("XmonC (Trace)");
  mcDEBUG_COMP("XmonC")
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

#define mcabsorbComp mcabsorbCompXmonC
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
#define mccompcurname  XmonC
#define mccompcurtype  PSDlin_monitor
#define mccompcurindex 5
#define nx mccXmonC_nx
#define PSDlin_N mccXmonC_PSDlin_N
#define PSDlin_p mccXmonC_PSDlin_p
#define PSDlin_p2 mccXmonC_PSDlin_p2
{   /* Declarations of XmonC=PSDlin_monitor() SETTING parameters. */
char* filename = mccXmonC_filename;
MCNUM xmin = mccXmonC_xmin;
MCNUM xmax = mccXmonC_xmax;
MCNUM ymin = mccXmonC_ymin;
MCNUM ymax = mccXmonC_ymax;
MCNUM xwidth = mccXmonC_xwidth;
MCNUM yheight = mccXmonC_yheight;
MCNUM restore_neutron = mccXmonC_restore_neutron;
int nowritefile = mccXmonC_nowritefile;
/* 'XmonC=PSDlin_monitor()' component instance has conditional execution */
if (( IsCold ))

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
#line 9301 "./ESS_2015_test.c"
}   /* End of XmonC=PSDlin_monitor() SETTING parameter declarations. */
#undef PSDlin_p2
#undef PSDlin_p
#undef PSDlin_N
#undef nx
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompXmonC:
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

  /* TRACE Component tofT [6] */
  mccoordschange(mcposrtofT, mcrotrtofT,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component tofT (without coords transformations) */
  mcJumpTrace_tofT:
  SIG_MESSAGE("tofT (Trace)");
  mcDEBUG_COMP("tofT")
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

#define mcabsorbComp mcabsorbComptofT
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
#define mccompcurname  tofT
#define mccompcurtype  TOF_monitor
#define mccompcurindex 6
#define nt mcctofT_nt
#define TOF_N mcctofT_TOF_N
#define TOF_p mcctofT_TOF_p
#define TOF_p2 mcctofT_TOF_p2
#define t_min mcctofT_t_min
#define t_max mcctofT_t_max
#define delta_t mcctofT_delta_t
{   /* Declarations of tofT=TOF_monitor() SETTING parameters. */
char* filename = mcctofT_filename;
MCNUM xmin = mcctofT_xmin;
MCNUM xmax = mcctofT_xmax;
MCNUM ymin = mcctofT_ymin;
MCNUM ymax = mcctofT_ymax;
MCNUM xwidth = mcctofT_xwidth;
MCNUM yheight = mcctofT_yheight;
MCNUM tmin = mcctofT_tmin;
MCNUM tmax = mcctofT_tmax;
MCNUM dt = mcctofT_dt;
MCNUM restore_neutron = mcctofT_restore_neutron;
int nowritefile = mcctofT_nowritefile;
/* 'tofT=TOF_monitor()' component instance has conditional execution */
if (( ! IsCold ))

#line 95 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/TOF_monitor.comp"
{
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
}
#line 9451 "./ESS_2015_test.c"
}   /* End of tofT=TOF_monitor() SETTING parameter declarations. */
#undef delta_t
#undef t_max
#undef t_min
#undef TOF_p2
#undef TOF_p
#undef TOF_N
#undef nt
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbComptofT:
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

  /* TRACE Component tofC [7] */
  mccoordschange(mcposrtofC, mcrotrtofC,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component tofC (without coords transformations) */
  mcJumpTrace_tofC:
  SIG_MESSAGE("tofC (Trace)");
  mcDEBUG_COMP("tofC")
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

#define mcabsorbComp mcabsorbComptofC
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
#define mccompcurname  tofC
#define mccompcurtype  TOF_monitor
#define mccompcurindex 7
#define nt mcctofC_nt
#define TOF_N mcctofC_TOF_N
#define TOF_p mcctofC_TOF_p
#define TOF_p2 mcctofC_TOF_p2
#define t_min mcctofC_t_min
#define t_max mcctofC_t_max
#define delta_t mcctofC_delta_t
{   /* Declarations of tofC=TOF_monitor() SETTING parameters. */
char* filename = mcctofC_filename;
MCNUM xmin = mcctofC_xmin;
MCNUM xmax = mcctofC_xmax;
MCNUM ymin = mcctofC_ymin;
MCNUM ymax = mcctofC_ymax;
MCNUM xwidth = mcctofC_xwidth;
MCNUM yheight = mcctofC_yheight;
MCNUM tmin = mcctofC_tmin;
MCNUM tmax = mcctofC_tmax;
MCNUM dt = mcctofC_dt;
MCNUM restore_neutron = mcctofC_restore_neutron;
int nowritefile = mcctofC_nowritefile;
/* 'tofC=TOF_monitor()' component instance has conditional execution */
if (( IsCold ))

#line 95 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/TOF_monitor.comp"
{
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
}
#line 9604 "./ESS_2015_test.c"
}   /* End of tofC=TOF_monitor() SETTING parameter declarations. */
#undef delta_t
#undef t_max
#undef t_min
#undef TOF_p2
#undef TOF_p
#undef TOF_N
#undef nt
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbComptofC:
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

  /* TRACE Component Slit1 [8] */
  mccoordschange(mcposrSlit1, mcrotrSlit1,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component Slit1 (without coords transformations) */
  mcJumpTrace_Slit1:
  SIG_MESSAGE("Slit1 (Trace)");
  mcDEBUG_COMP("Slit1")
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

#define mcabsorbComp mcabsorbCompSlit1
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
#define mccompcurname  Slit1
#define mccompcurtype  Slit
#define mccompcurindex 8
{   /* Declarations of Slit1=Slit() SETTING parameters. */
MCNUM xmin = mccSlit1_xmin;
MCNUM xmax = mccSlit1_xmax;
MCNUM ymin = mccSlit1_ymin;
MCNUM ymax = mccSlit1_ymax;
MCNUM radius = mccSlit1_radius;
MCNUM xwidth = mccSlit1_xwidth;
MCNUM yheight = mccSlit1_yheight;
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
#line 9735 "./ESS_2015_test.c"
}   /* End of Slit1=Slit() SETTING parameter declarations. */
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompSlit1:
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

  /* TRACE Component LamT [9] */
  mccoordschange(mcposrLamT, mcrotrLamT,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component LamT (without coords transformations) */
  mcJumpTrace_LamT:
  SIG_MESSAGE("LamT (Trace)");
  mcDEBUG_COMP("LamT")
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

#define mcabsorbComp mcabsorbCompLamT
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
#define mccompcurname  LamT
#define mccompcurtype  L_monitor
#define mccompcurindex 9
#define nL mccLamT_nL
#define L_N mccLamT_L_N
#define L_p mccLamT_L_p
#define L_p2 mccLamT_L_p2
{   /* Declarations of LamT=L_monitor() SETTING parameters. */
char* filename = mccLamT_filename;
MCNUM xmin = mccLamT_xmin;
MCNUM xmax = mccLamT_xmax;
MCNUM ymin = mccLamT_ymin;
MCNUM ymax = mccLamT_ymax;
MCNUM xwidth = mccLamT_xwidth;
MCNUM yheight = mccLamT_yheight;
MCNUM Lmin = mccLamT_Lmin;
MCNUM Lmax = mccLamT_Lmax;
MCNUM restore_neutron = mccLamT_restore_neutron;
int nowritefile = mccLamT_nowritefile;
/* 'LamT=L_monitor()' component instance has conditional execution */
if (( ! IsCold ))

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
#line 9880 "./ESS_2015_test.c"
}   /* End of LamT=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompLamT:
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

  /* TRACE Component LamC [10] */
  mccoordschange(mcposrLamC, mcrotrLamC,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component LamC (without coords transformations) */
  mcJumpTrace_LamC:
  SIG_MESSAGE("LamC (Trace)");
  mcDEBUG_COMP("LamC")
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

#define mcabsorbComp mcabsorbCompLamC
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
#define mccompcurname  LamC
#define mccompcurtype  L_monitor
#define mccompcurindex 10
#define nL mccLamC_nL
#define L_N mccLamC_L_N
#define L_p mccLamC_L_p
#define L_p2 mccLamC_L_p2
{   /* Declarations of LamC=L_monitor() SETTING parameters. */
char* filename = mccLamC_filename;
MCNUM xmin = mccLamC_xmin;
MCNUM xmax = mccLamC_xmax;
MCNUM ymin = mccLamC_ymin;
MCNUM ymax = mccLamC_ymax;
MCNUM xwidth = mccLamC_xwidth;
MCNUM yheight = mccLamC_yheight;
MCNUM Lmin = mccLamC_Lmin;
MCNUM Lmax = mccLamC_Lmax;
MCNUM restore_neutron = mccLamC_restore_neutron;
int nowritefile = mccLamC_nowritefile;
/* 'LamC=L_monitor()' component instance has conditional execution */
if (( IsCold ))

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
#line 10029 "./ESS_2015_test.c"
}   /* End of LamC=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompLamC:
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
#line 10142 "./ESS_2015_test.c"
}   /* End of Origin=Progress_bar() SETTING parameter declarations. */
#undef CurrentTime
#undef EndTime
#undef StartTime
#undef IntermediateCnts
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'XmonT'. */
  SIG_MESSAGE("XmonT (Save)");
#define mccompcurname  XmonT
#define mccompcurtype  PSDlin_monitor
#define mccompcurindex 4
#define nx mccXmonT_nx
#define PSDlin_N mccXmonT_PSDlin_N
#define PSDlin_p mccXmonT_PSDlin_p
#define PSDlin_p2 mccXmonT_PSDlin_p2
{   /* Declarations of XmonT=PSDlin_monitor() SETTING parameters. */
char* filename = mccXmonT_filename;
MCNUM xmin = mccXmonT_xmin;
MCNUM xmax = mccXmonT_xmax;
MCNUM ymin = mccXmonT_ymin;
MCNUM ymax = mccXmonT_ymax;
MCNUM xwidth = mccXmonT_xwidth;
MCNUM yheight = mccXmonT_yheight;
MCNUM restore_neutron = mccXmonT_restore_neutron;
int nowritefile = mccXmonT_nowritefile;
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
#line 10183 "./ESS_2015_test.c"
}   /* End of XmonT=PSDlin_monitor() SETTING parameter declarations. */
#undef PSDlin_p2
#undef PSDlin_p
#undef PSDlin_N
#undef nx
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'XmonC'. */
  SIG_MESSAGE("XmonC (Save)");
#define mccompcurname  XmonC
#define mccompcurtype  PSDlin_monitor
#define mccompcurindex 5
#define nx mccXmonC_nx
#define PSDlin_N mccXmonC_PSDlin_N
#define PSDlin_p mccXmonC_PSDlin_p
#define PSDlin_p2 mccXmonC_PSDlin_p2
{   /* Declarations of XmonC=PSDlin_monitor() SETTING parameters. */
char* filename = mccXmonC_filename;
MCNUM xmin = mccXmonC_xmin;
MCNUM xmax = mccXmonC_xmax;
MCNUM ymin = mccXmonC_ymin;
MCNUM ymax = mccXmonC_ymax;
MCNUM xwidth = mccXmonC_xwidth;
MCNUM yheight = mccXmonC_yheight;
MCNUM restore_neutron = mccXmonC_restore_neutron;
int nowritefile = mccXmonC_nowritefile;
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
#line 10224 "./ESS_2015_test.c"
}   /* End of XmonC=PSDlin_monitor() SETTING parameter declarations. */
#undef PSDlin_p2
#undef PSDlin_p
#undef PSDlin_N
#undef nx
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'tofT'. */
  SIG_MESSAGE("tofT (Save)");
#define mccompcurname  tofT
#define mccompcurtype  TOF_monitor
#define mccompcurindex 6
#define nt mcctofT_nt
#define TOF_N mcctofT_TOF_N
#define TOF_p mcctofT_TOF_p
#define TOF_p2 mcctofT_TOF_p2
#define t_min mcctofT_t_min
#define t_max mcctofT_t_max
#define delta_t mcctofT_delta_t
{   /* Declarations of tofT=TOF_monitor() SETTING parameters. */
char* filename = mcctofT_filename;
MCNUM xmin = mcctofT_xmin;
MCNUM xmax = mcctofT_xmax;
MCNUM ymin = mcctofT_ymin;
MCNUM ymax = mcctofT_ymax;
MCNUM xwidth = mcctofT_xwidth;
MCNUM yheight = mcctofT_yheight;
MCNUM tmin = mcctofT_tmin;
MCNUM tmax = mcctofT_tmax;
MCNUM dt = mcctofT_dt;
MCNUM restore_neutron = mcctofT_restore_neutron;
int nowritefile = mcctofT_nowritefile;
#line 115 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/TOF_monitor.comp"
{
    if (!nowritefile) {
    DETECTOR_OUT_1D(
        "Time-of-flight monitor",
        "Time-of-flight [\\gms]",
        "Intensity",
        "t", t_min, t_max, nt,
        &TOF_N[0],&TOF_p[0],&TOF_p2[0],
        filename);
    }
}
#line 10271 "./ESS_2015_test.c"
}   /* End of tofT=TOF_monitor() SETTING parameter declarations. */
#undef delta_t
#undef t_max
#undef t_min
#undef TOF_p2
#undef TOF_p
#undef TOF_N
#undef nt
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'tofC'. */
  SIG_MESSAGE("tofC (Save)");
#define mccompcurname  tofC
#define mccompcurtype  TOF_monitor
#define mccompcurindex 7
#define nt mcctofC_nt
#define TOF_N mcctofC_TOF_N
#define TOF_p mcctofC_TOF_p
#define TOF_p2 mcctofC_TOF_p2
#define t_min mcctofC_t_min
#define t_max mcctofC_t_max
#define delta_t mcctofC_delta_t
{   /* Declarations of tofC=TOF_monitor() SETTING parameters. */
char* filename = mcctofC_filename;
MCNUM xmin = mcctofC_xmin;
MCNUM xmax = mcctofC_xmax;
MCNUM ymin = mcctofC_ymin;
MCNUM ymax = mcctofC_ymax;
MCNUM xwidth = mcctofC_xwidth;
MCNUM yheight = mcctofC_yheight;
MCNUM tmin = mcctofC_tmin;
MCNUM tmax = mcctofC_tmax;
MCNUM dt = mcctofC_dt;
MCNUM restore_neutron = mcctofC_restore_neutron;
int nowritefile = mcctofC_nowritefile;
#line 115 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/TOF_monitor.comp"
{
    if (!nowritefile) {
    DETECTOR_OUT_1D(
        "Time-of-flight monitor",
        "Time-of-flight [\\gms]",
        "Intensity",
        "t", t_min, t_max, nt,
        &TOF_N[0],&TOF_p[0],&TOF_p2[0],
        filename);
    }
}
#line 10321 "./ESS_2015_test.c"
}   /* End of tofC=TOF_monitor() SETTING parameter declarations. */
#undef delta_t
#undef t_max
#undef t_min
#undef TOF_p2
#undef TOF_p
#undef TOF_N
#undef nt
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'LamT'. */
  SIG_MESSAGE("LamT (Save)");
#define mccompcurname  LamT
#define mccompcurtype  L_monitor
#define mccompcurindex 9
#define nL mccLamT_nL
#define L_N mccLamT_L_N
#define L_p mccLamT_L_p
#define L_p2 mccLamT_L_p2
{   /* Declarations of LamT=L_monitor() SETTING parameters. */
char* filename = mccLamT_filename;
MCNUM xmin = mccLamT_xmin;
MCNUM xmax = mccLamT_xmax;
MCNUM ymin = mccLamT_ymin;
MCNUM ymax = mccLamT_ymax;
MCNUM xwidth = mccLamT_xwidth;
MCNUM yheight = mccLamT_yheight;
MCNUM Lmin = mccLamT_Lmin;
MCNUM Lmax = mccLamT_Lmax;
MCNUM restore_neutron = mccLamT_restore_neutron;
int nowritefile = mccLamT_nowritefile;
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
#line 10367 "./ESS_2015_test.c"
}   /* End of LamT=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'LamC'. */
  SIG_MESSAGE("LamC (Save)");
#define mccompcurname  LamC
#define mccompcurtype  L_monitor
#define mccompcurindex 10
#define nL mccLamC_nL
#define L_N mccLamC_L_N
#define L_p mccLamC_L_p
#define L_p2 mccLamC_L_p2
{   /* Declarations of LamC=L_monitor() SETTING parameters. */
char* filename = mccLamC_filename;
MCNUM xmin = mccLamC_xmin;
MCNUM xmax = mccLamC_xmax;
MCNUM ymin = mccLamC_ymin;
MCNUM ymax = mccLamC_ymax;
MCNUM xwidth = mccLamC_xwidth;
MCNUM yheight = mccLamC_yheight;
MCNUM Lmin = mccLamC_Lmin;
MCNUM Lmax = mccLamC_Lmax;
MCNUM restore_neutron = mccLamC_restore_neutron;
int nowritefile = mccLamC_nowritefile;
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
#line 10410 "./ESS_2015_test.c"
}   /* End of LamC=L_monitor() SETTING parameter declarations. */
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
#line 10454 "./ESS_2015_test.c"
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
    if (mcAbsorbProp[2]) fprintf(stderr, "Warning: %g events were removed in Component[2] Source=ESS_moderator()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[2]);
    if (!mcNCounter[3]) fprintf(stderr, "Warning: No neutron could reach Component[3] PortOrig\n");
    if (mcAbsorbProp[3]) fprintf(stderr, "Warning: %g events were removed in Component[3] PortOrig=Arm()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[3]);
    if (!mcNCounter[4]) fprintf(stderr, "Warning: No neutron could reach Component[4] XmonT\n");
    if (mcAbsorbProp[4]) fprintf(stderr, "Warning: %g events were removed in Component[4] XmonT=PSDlin_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[4]);
    if (!mcNCounter[5]) fprintf(stderr, "Warning: No neutron could reach Component[5] XmonC\n");
    if (mcAbsorbProp[5]) fprintf(stderr, "Warning: %g events were removed in Component[5] XmonC=PSDlin_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[5]);
    if (!mcNCounter[6]) fprintf(stderr, "Warning: No neutron could reach Component[6] tofT\n");
    if (mcAbsorbProp[6]) fprintf(stderr, "Warning: %g events were removed in Component[6] tofT=TOF_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[6]);
    if (!mcNCounter[7]) fprintf(stderr, "Warning: No neutron could reach Component[7] tofC\n");
    if (mcAbsorbProp[7]) fprintf(stderr, "Warning: %g events were removed in Component[7] tofC=TOF_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[7]);
    if (!mcNCounter[8]) fprintf(stderr, "Warning: No neutron could reach Component[8] Slit1\n");
    if (mcAbsorbProp[8]) fprintf(stderr, "Warning: %g events were removed in Component[8] Slit1=Slit()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[8]);
    if (!mcNCounter[9]) fprintf(stderr, "Warning: No neutron could reach Component[9] LamT\n");
    if (mcAbsorbProp[9]) fprintf(stderr, "Warning: %g events were removed in Component[9] LamT=L_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[9]);
    if (!mcNCounter[10]) fprintf(stderr, "Warning: No neutron could reach Component[10] LamC\n");
    if (mcAbsorbProp[10]) fprintf(stderr, "Warning: %g events were removed in Component[10] LamC=L_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[10]);
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
#line 10518 "./ESS_2015_test.c"
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
#define mccompcurtype  ESS_moderator
#define mccompcurindex 2
#define l_range mccSource_l_range
#define w_mult mccSource_w_mult
#define w_geom mccSource_w_geom
#define w_geom_c mccSource_w_geom_c
#define w_geom_t mccSource_w_geom_t
#define w_stat mccSource_w_stat
#define cold mccSource_cold
#define wasleft mccSource_wasleft
#define cosine mccSource_cosine
#define tx mccSource_tx
#define ty mccSource_ty
#define tz mccSource_tz
#define cx mccSource_cx
#define cy mccSource_cy
#define cz mccSource_cz
#define cyl_radius mccSource_cyl_radius
#define t1x mccSource_t1x
#define t1y mccSource_t1y
#define t1z mccSource_t1z
#define t2x mccSource_t2x
#define t2y mccSource_t2y
#define t2z mccSource_t2z
#define modextras mccSource_modextras
#define cold_bril mccSource_cold_bril
#define thermal_bril mccSource_thermal_bril
#define tfocus_width mccSource_tfocus_width
#define tfocus_time mccSource_tfocus_time
#define dt mccSource_dt
#define r_empty mccSource_r_empty
#define lambda mccSource_lambda
#define xprime mccSource_xprime
#define yprime mccSource_yprime
#define zprime mccSource_zprime
#define sgnport mccSource_sgnport
#define abs_beamport_angle mccSource_abs_beamport_angle
#define cos_beamport_angle mccSource_cos_beamport_angle
#define edge_thermal mccSource_edge_thermal
#define cos_thermal mccSource_cos_thermal
#define cos_cold mccSource_cos_cold
#define sin_beamport_angle mccSource_sin_beamport_angle
{   /* Declarations of Source=ESS_moderator() SETTING parameters. */
MCNUM isleft = mccSource_isleft;
MCNUM Lmin = mccSource_Lmin;
MCNUM Lmax = mccSource_Lmax;
MCNUM cold_frac = mccSource_cold_frac;
MCNUM dist = mccSource_dist;
MCNUM focus_xw = mccSource_focus_xw;
MCNUM focus_yh = mccSource_focus_yh;
int target_index = mccSource_target_index;
MCNUM tmax_multiplier = mccSource_tmax_multiplier;
MCNUM yheight_c = mccSource_yheight_c;
MCNUM yheight_t = mccSource_yheight_t;
int n_pulses = mccSource_n_pulses;
MCNUM acc_power = mccSource_acc_power;
MCNUM beamport_angle = mccSource_beamport_angle;
char* sourcedef = mccSource_sourcedef;
MCNUM xwidth_c = mccSource_xwidth_c;
MCNUM xwidth_t = mccSource_xwidth_t;
MCNUM extraction_opening = mccSource_extraction_opening;
#line 440 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../sources/ESS_moderator.comp"
{
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

}
#line 10751 "./ESS_2015_test.c"
}   /* End of Source=ESS_moderator() SETTING parameter declarations. */
#undef sin_beamport_angle
#undef cos_cold
#undef cos_thermal
#undef edge_thermal
#undef cos_beamport_angle
#undef abs_beamport_angle
#undef sgnport
#undef zprime
#undef yprime
#undef xprime
#undef lambda
#undef r_empty
#undef dt
#undef tfocus_time
#undef tfocus_width
#undef thermal_bril
#undef cold_bril
#undef modextras
#undef t2z
#undef t2y
#undef t2x
#undef t1z
#undef t1y
#undef t1x
#undef cyl_radius
#undef cz
#undef cy
#undef cx
#undef tz
#undef ty
#undef tx
#undef cosine
#undef wasleft
#undef cold
#undef w_stat
#undef w_geom_t
#undef w_geom_c
#undef w_geom
#undef w_mult
#undef l_range
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'PortOrig'. */
  SIG_MESSAGE("PortOrig (McDisplay)");
  printf("MCDISPLAY: component %s\n", "PortOrig");
#define mccompcurname  PortOrig
#define mccompcurtype  Arm
#define mccompcurindex 3
#line 40 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Arm.comp"
{
  /* A bit ugly; hard-coded dimensions. */
  
  line(0,0,0,0.2,0,0);
  line(0,0,0,0,0.2,0);
  line(0,0,0,0,0,0.2);
}
#line 10811 "./ESS_2015_test.c"
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'XmonT'. */
  SIG_MESSAGE("XmonT (McDisplay)");
  printf("MCDISPLAY: component %s\n", "XmonT");
#define mccompcurname  XmonT
#define mccompcurtype  PSDlin_monitor
#define mccompcurindex 4
#define nx mccXmonT_nx
#define PSDlin_N mccXmonT_PSDlin_N
#define PSDlin_p mccXmonT_PSDlin_p
#define PSDlin_p2 mccXmonT_PSDlin_p2
{   /* Declarations of XmonT=PSDlin_monitor() SETTING parameters. */
char* filename = mccXmonT_filename;
MCNUM xmin = mccXmonT_xmin;
MCNUM xmax = mccXmonT_xmax;
MCNUM ymin = mccXmonT_ymin;
MCNUM ymax = mccXmonT_ymax;
MCNUM xwidth = mccXmonT_xwidth;
MCNUM yheight = mccXmonT_yheight;
MCNUM restore_neutron = mccXmonT_restore_neutron;
int nowritefile = mccXmonT_nowritefile;
#line 116 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/PSDlin_monitor.comp"
{
  
  multiline(5, (double)xmin, (double)ymin, 0.0,
               (double)xmax, (double)ymin, 0.0,
               (double)xmax, (double)ymax, 0.0,
               (double)xmin, (double)ymax, 0.0,
               (double)xmin, (double)ymin, 0.0);
}
#line 10845 "./ESS_2015_test.c"
}   /* End of XmonT=PSDlin_monitor() SETTING parameter declarations. */
#undef PSDlin_p2
#undef PSDlin_p
#undef PSDlin_N
#undef nx
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'XmonC'. */
  SIG_MESSAGE("XmonC (McDisplay)");
  printf("MCDISPLAY: component %s\n", "XmonC");
#define mccompcurname  XmonC
#define mccompcurtype  PSDlin_monitor
#define mccompcurindex 5
#define nx mccXmonC_nx
#define PSDlin_N mccXmonC_PSDlin_N
#define PSDlin_p mccXmonC_PSDlin_p
#define PSDlin_p2 mccXmonC_PSDlin_p2
{   /* Declarations of XmonC=PSDlin_monitor() SETTING parameters. */
char* filename = mccXmonC_filename;
MCNUM xmin = mccXmonC_xmin;
MCNUM xmax = mccXmonC_xmax;
MCNUM ymin = mccXmonC_ymin;
MCNUM ymax = mccXmonC_ymax;
MCNUM xwidth = mccXmonC_xwidth;
MCNUM yheight = mccXmonC_yheight;
MCNUM restore_neutron = mccXmonC_restore_neutron;
int nowritefile = mccXmonC_nowritefile;
#line 116 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/PSDlin_monitor.comp"
{
  
  multiline(5, (double)xmin, (double)ymin, 0.0,
               (double)xmax, (double)ymin, 0.0,
               (double)xmax, (double)ymax, 0.0,
               (double)xmin, (double)ymax, 0.0,
               (double)xmin, (double)ymin, 0.0);
}
#line 10884 "./ESS_2015_test.c"
}   /* End of XmonC=PSDlin_monitor() SETTING parameter declarations. */
#undef PSDlin_p2
#undef PSDlin_p
#undef PSDlin_N
#undef nx
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'tofT'. */
  SIG_MESSAGE("tofT (McDisplay)");
  printf("MCDISPLAY: component %s\n", "tofT");
#define mccompcurname  tofT
#define mccompcurtype  TOF_monitor
#define mccompcurindex 6
#define nt mcctofT_nt
#define TOF_N mcctofT_TOF_N
#define TOF_p mcctofT_TOF_p
#define TOF_p2 mcctofT_TOF_p2
#define t_min mcctofT_t_min
#define t_max mcctofT_t_max
#define delta_t mcctofT_delta_t
{   /* Declarations of tofT=TOF_monitor() SETTING parameters. */
char* filename = mcctofT_filename;
MCNUM xmin = mcctofT_xmin;
MCNUM xmax = mcctofT_xmax;
MCNUM ymin = mcctofT_ymin;
MCNUM ymax = mcctofT_ymax;
MCNUM xwidth = mcctofT_xwidth;
MCNUM yheight = mcctofT_yheight;
MCNUM tmin = mcctofT_tmin;
MCNUM tmax = mcctofT_tmax;
MCNUM dt = mcctofT_dt;
MCNUM restore_neutron = mcctofT_restore_neutron;
int nowritefile = mcctofT_nowritefile;
#line 128 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/TOF_monitor.comp"
{
  
  multiline(5, (double)xmin, (double)ymin, 0.0,
               (double)xmax, (double)ymin, 0.0,
               (double)xmax, (double)ymax, 0.0,
               (double)xmin, (double)ymax, 0.0,
               (double)xmin, (double)ymin, 0.0);
}
#line 10929 "./ESS_2015_test.c"
}   /* End of tofT=TOF_monitor() SETTING parameter declarations. */
#undef delta_t
#undef t_max
#undef t_min
#undef TOF_p2
#undef TOF_p
#undef TOF_N
#undef nt
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'tofC'. */
  SIG_MESSAGE("tofC (McDisplay)");
  printf("MCDISPLAY: component %s\n", "tofC");
#define mccompcurname  tofC
#define mccompcurtype  TOF_monitor
#define mccompcurindex 7
#define nt mcctofC_nt
#define TOF_N mcctofC_TOF_N
#define TOF_p mcctofC_TOF_p
#define TOF_p2 mcctofC_TOF_p2
#define t_min mcctofC_t_min
#define t_max mcctofC_t_max
#define delta_t mcctofC_delta_t
{   /* Declarations of tofC=TOF_monitor() SETTING parameters. */
char* filename = mcctofC_filename;
MCNUM xmin = mcctofC_xmin;
MCNUM xmax = mcctofC_xmax;
MCNUM ymin = mcctofC_ymin;
MCNUM ymax = mcctofC_ymax;
MCNUM xwidth = mcctofC_xwidth;
MCNUM yheight = mcctofC_yheight;
MCNUM tmin = mcctofC_tmin;
MCNUM tmax = mcctofC_tmax;
MCNUM dt = mcctofC_dt;
MCNUM restore_neutron = mcctofC_restore_neutron;
int nowritefile = mcctofC_nowritefile;
#line 128 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/TOF_monitor.comp"
{
  
  multiline(5, (double)xmin, (double)ymin, 0.0,
               (double)xmax, (double)ymin, 0.0,
               (double)xmax, (double)ymax, 0.0,
               (double)xmin, (double)ymax, 0.0,
               (double)xmin, (double)ymin, 0.0);
}
#line 10977 "./ESS_2015_test.c"
}   /* End of tofC=TOF_monitor() SETTING parameter declarations. */
#undef delta_t
#undef t_max
#undef t_min
#undef TOF_p2
#undef TOF_p
#undef TOF_N
#undef nt
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'Slit1'. */
  SIG_MESSAGE("Slit1 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "Slit1");
#define mccompcurname  Slit1
#define mccompcurtype  Slit
#define mccompcurindex 8
{   /* Declarations of Slit1=Slit() SETTING parameters. */
MCNUM xmin = mccSlit1_xmin;
MCNUM xmax = mccSlit1_xmax;
MCNUM ymin = mccSlit1_ymin;
MCNUM ymax = mccSlit1_ymax;
MCNUM radius = mccSlit1_radius;
MCNUM xwidth = mccSlit1_xwidth;
MCNUM yheight = mccSlit1_yheight;
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
#line 11027 "./ESS_2015_test.c"
}   /* End of Slit1=Slit() SETTING parameter declarations. */
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'LamT'. */
  SIG_MESSAGE("LamT (McDisplay)");
  printf("MCDISPLAY: component %s\n", "LamT");
#define mccompcurname  LamT
#define mccompcurtype  L_monitor
#define mccompcurindex 9
#define nL mccLamT_nL
#define L_N mccLamT_L_N
#define L_p mccLamT_L_p
#define L_p2 mccLamT_L_p2
{   /* Declarations of LamT=L_monitor() SETTING parameters. */
char* filename = mccLamT_filename;
MCNUM xmin = mccLamT_xmin;
MCNUM xmax = mccLamT_xmax;
MCNUM ymin = mccLamT_ymin;
MCNUM ymax = mccLamT_ymax;
MCNUM xwidth = mccLamT_xwidth;
MCNUM yheight = mccLamT_yheight;
MCNUM Lmin = mccLamT_Lmin;
MCNUM Lmax = mccLamT_Lmax;
MCNUM restore_neutron = mccLamT_restore_neutron;
int nowritefile = mccLamT_nowritefile;
#line 120 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/L_monitor.comp"
{
  
  multiline(5, (double)xmin, (double)ymin, 0.0,
               (double)xmax, (double)ymin, 0.0,
               (double)xmax, (double)ymax, 0.0,
               (double)xmin, (double)ymax, 0.0,
               (double)xmin, (double)ymin, 0.0);
}
#line 11064 "./ESS_2015_test.c"
}   /* End of LamT=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'LamC'. */
  SIG_MESSAGE("LamC (McDisplay)");
  printf("MCDISPLAY: component %s\n", "LamC");
#define mccompcurname  LamC
#define mccompcurtype  L_monitor
#define mccompcurindex 10
#define nL mccLamC_nL
#define L_N mccLamC_L_N
#define L_p mccLamC_L_p
#define L_p2 mccLamC_L_p2
{   /* Declarations of LamC=L_monitor() SETTING parameters. */
char* filename = mccLamC_filename;
MCNUM xmin = mccLamC_xmin;
MCNUM xmax = mccLamC_xmax;
MCNUM ymin = mccLamC_ymin;
MCNUM ymax = mccLamC_ymax;
MCNUM xwidth = mccLamC_xwidth;
MCNUM yheight = mccLamC_yheight;
MCNUM Lmin = mccLamC_Lmin;
MCNUM Lmax = mccLamC_Lmax;
MCNUM restore_neutron = mccLamC_restore_neutron;
int nowritefile = mccLamC_nowritefile;
#line 120 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/L_monitor.comp"
{
  
  multiline(5, (double)xmin, (double)ymin, 0.0,
               (double)xmax, (double)ymin, 0.0,
               (double)xmax, (double)ymax, 0.0,
               (double)xmin, (double)ymax, 0.0,
               (double)xmin, (double)ymin, 0.0);
}
#line 11105 "./ESS_2015_test.c"
}   /* End of LamC=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef nL
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
/* end of generated C code ./ESS_2015_test.c */
