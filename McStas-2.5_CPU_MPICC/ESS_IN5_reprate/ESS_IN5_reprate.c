/* Automatically generated file. Do not edit. 
 * Format:     ANSI C source code
 * Creator:    McStas <http://www.mcstas.org>
 * Instrument: /zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr (ESS_IN5_reprate)
 * Date:       Tue Feb 25 20:23:35 2020
 * File:       ./ESS_IN5_reprate.c
 * Compile:    cc -o ESS_IN5_reprate.out ./ESS_IN5_reprate.c 
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

#line 712 "./ESS_IN5_reprate.c"

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

#line 945 "./ESS_IN5_reprate.c"

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

#line 4977 "./ESS_IN5_reprate.c"

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

#line 5337 "./ESS_IN5_reprate.c"
#ifdef MC_TRACE_ENABLED
int mctraceenabled = 1;
#else
int mctraceenabled = 0;
#endif
#define MCSTAS "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../"
int mcdefaultmain = 1;
char mcinstrument_name[] = "ESS_IN5_reprate";
char mcinstrument_source[] = "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr";
char *mcinstrument_exe=NULL; /* will be set to argv[0] in main */
int main(int argc, char *argv[]){return mccode_main(argc, argv);}
void mcinit(void);
void mcraytrace(void);
void mcsave(FILE *);
void mcfinally(void);
void mcdisplay(void);

/* Shared user declarations for all components 'ESS_moderator_long'. */
#line 140 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../obsolete/ESS_moderator_long.comp"
double Mezei_M_fct(double l, double temp)
  {
    double a=949.0/temp;
    return 2*a*a*exp(-a/(l*l))/(l*l*l*l*l);
  }

  double Mezei_F_fct(double t, double tau, int n)
  {
    return (exp(-t/tau)-exp(-n*t/tau))*n/(n-1)/tau;
  }
#line 5367 "./ESS_IN5_reprate.c"

/* Shared user declarations for all components 'Guide'. */
#line 63 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide.comp"
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

#line 6962 "./ESS_IN5_reprate.c"

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

#line 7048 "./ESS_IN5_reprate.c"

/* Shared user declarations for all components 'Tunneling_sample'. */
#line 92 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../samples/Tunneling_sample.comp"
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
#line 7065 "./ESS_IN5_reprate.c"

/* Instrument parameters. */
MCNUM mcipLmin;
MCNUM mcipLmax;
MCNUM mciplambda0;
MCNUM mcipPulse_width;
MCNUM mcipNum_pulses;
MCNUM mcipGUI_start;
MCNUM mcipFO1_DIST;
MCNUM mcipL_ballistic_begin;
MCNUM mcipL_ballistic_end;
MCNUM mcipLength;
MCNUM mcipSAMPLE_DIST;
MCNUM mcipDETECTOR_DIST;
MCNUM mcipGUI_h;
MCNUM mcipGUI_w;
MCNUM mcipGUI_GAP;
MCNUM mcipH1;
MCNUM mcipW1;
MCNUM mcipH2;
MCNUM mcipW2;
MCNUM mcipH3;
MCNUM mcipW3;
MCNUM mcipH4;
MCNUM mcipW4;
MCNUM mcipH_chop;
MCNUM mcipW_chop;
MCNUM mcipH_end;
MCNUM mcipW_end;
MCNUM mcipALPHA;
MCNUM mcipM;
MCNUM mcipF_slow1;
MCNUM mcipF_slow2;
MCNUM mcipF_fast1;
MCNUM mcipF_fast2;
MCNUM mcipN_fast;
MCNUM mcipSLOW1_THETA;
MCNUM mcipFO3;
MCNUM mcipTHETA_fast1;
MCNUM mcipFAST_THETA;
MCNUM mcipGamma;
MCNUM mcipEtun;
MCNUM mcipV_HOLE;
MCNUM mcipFRAC_QUASIEL;
MCNUM mcipFRAC_TUNNEL;
MCNUM mcipTT;
MCNUM mcipRES_DE;
MCNUM mcipport;
MCNUM mcipcold;

#define mcNUMIPAR 47
int mcnumipar = 47;
struct mcinputtable_struct mcinputtable[mcNUMIPAR+1] = {
  "Lmin", &mcipLmin, instr_type_double, "4.9", 
  "Lmax", &mcipLmax, instr_type_double, "5.1", 
  "lambda0", &mciplambda0, instr_type_double, "5", 
  "Pulse_width", &mcipPulse_width, instr_type_double, "0.002", 
  "Num_pulses", &mcipNum_pulses, instr_type_double, "1", 
  "GUI_start", &mcipGUI_start, instr_type_double, "2.0", 
  "FO1_DIST", &mcipFO1_DIST, instr_type_double, "6", 
  "L_ballistic_begin", &mcipL_ballistic_begin, instr_type_double, "19.5", 
  "L_ballistic_end", &mcipL_ballistic_end, instr_type_double, "17", 
  "Length", &mcipLength, instr_type_double, "100", 
  "SAMPLE_DIST", &mcipSAMPLE_DIST, instr_type_double, "1.2", 
  "DETECTOR_DIST", &mcipDETECTOR_DIST, instr_type_double, "4", 
  "GUI_h", &mcipGUI_h, instr_type_double, "0.105", 
  "GUI_w", &mcipGUI_w, instr_type_double, "0.1", 
  "GUI_GAP", &mcipGUI_GAP, instr_type_double, "0.05", 
  "H1", &mcipH1, instr_type_double, "0.167", 
  "W1", &mcipW1, instr_type_double, "0.116", 
  "H2", &mcipH2, instr_type_double, "0.185", 
  "W2", &mcipW2, instr_type_double, "0.15", 
  "H3", &mcipH3, instr_type_double, "0.19", 
  "W3", &mcipW3, instr_type_double, "0.15", 
  "H4", &mcipH4, instr_type_double, "0.213", 
  "W4", &mcipW4, instr_type_double, "0.14", 
  "H_chop", &mcipH_chop, instr_type_double, "0.075", 
  "W_chop", &mcipW_chop, instr_type_double, "0.03", 
  "H_end", &mcipH_end, instr_type_double, "0.042", 
  "W_end", &mcipW_end, instr_type_double, "0.0215", 
  "ALPHA", &mcipALPHA, instr_type_double, "3.4", 
  "M", &mcipM, instr_type_double, "3.5", 
  "F_slow1", &mcipF_slow1, instr_type_double, "16.6667", 
  "F_slow2", &mcipF_slow2, instr_type_double, "0", 
  "F_fast1", &mcipF_fast1, instr_type_double, "0", 
  "F_fast2", &mcipF_fast2, instr_type_double, "200", 
  "N_fast", &mcipN_fast, instr_type_double, "1", 
  "SLOW1_THETA", &mcipSLOW1_THETA, instr_type_double, "120", 
  "FO3", &mcipFO3, instr_type_double, "1", 
  "THETA_fast1", &mcipTHETA_fast1, instr_type_double, "180", 
  "FAST_THETA", &mcipFAST_THETA, instr_type_double, "5", 
  "Gamma", &mcipGamma, instr_type_double, "0", 
  "Etun", &mcipEtun, instr_type_double, "1", 
  "V_HOLE", &mcipV_HOLE, instr_type_double, "0", 
  "FRAC_QUASIEL", &mcipFRAC_QUASIEL, instr_type_double, "0", 
  "FRAC_TUNNEL", &mcipFRAC_TUNNEL, instr_type_double, "0", 
  "TT", &mcipTT, instr_type_double, "50", 
  "RES_DE", &mcipRES_DE, instr_type_double, "0.5", 
  "port", &mcipport, instr_type_double, "30", 
  "cold", &mcipcold, instr_type_double, "0.95", 
  NULL, NULL, instr_type_double, ""
};

/* User declarations from instrument definition. */
#define mccompcurname  ESS_IN5_reprate
#define mccompcurtype  INSTRUMENT
#define mccompcurindex 0
#define mcposaESS_IN5_reprate coords_set(0,0,0)
#define Lmin mcipLmin
#define Lmax mcipLmax
#define lambda0 mciplambda0
#define Pulse_width mcipPulse_width
#define Num_pulses mcipNum_pulses
#define GUI_start mcipGUI_start
#define FO1_DIST mcipFO1_DIST
#define L_ballistic_begin mcipL_ballistic_begin
#define L_ballistic_end mcipL_ballistic_end
#define Length mcipLength
#define SAMPLE_DIST mcipSAMPLE_DIST
#define DETECTOR_DIST mcipDETECTOR_DIST
#define GUI_h mcipGUI_h
#define GUI_w mcipGUI_w
#define GUI_GAP mcipGUI_GAP
#define H1 mcipH1
#define W1 mcipW1
#define H2 mcipH2
#define W2 mcipW2
#define H3 mcipH3
#define W3 mcipW3
#define H4 mcipH4
#define W4 mcipW4
#define H_chop mcipH_chop
#define W_chop mcipW_chop
#define H_end mcipH_end
#define W_end mcipW_end
#define ALPHA mcipALPHA
#define M mcipM
#define F_slow1 mcipF_slow1
#define F_slow2 mcipF_slow2
#define F_fast1 mcipF_fast1
#define F_fast2 mcipF_fast2
#define N_fast mcipN_fast
#define SLOW1_THETA mcipSLOW1_THETA
#define FO3 mcipFO3
#define THETA_fast1 mcipTHETA_fast1
#define FAST_THETA mcipFAST_THETA
#define Gamma mcipGamma
#define Etun mcipEtun
#define V_HOLE mcipV_HOLE
#define FRAC_QUASIEL mcipFRAC_QUASIEL
#define FRAC_TUNNEL mcipFRAC_TUNNEL
#define TT mcipTT
#define RES_DE mcipRES_DE
#define port mcipport
#define cold mcipcold
#line 88 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
        double FREQ;
        double t_FO1, t_FO2, t_fast1, t_fast2, t_fast2a, t_fast3, t_sample, t_detector;
        double tmin_zoom, tmax_zoom, t_offset;
        double E_target;
#line 7226 "./ESS_IN5_reprate.c"
#undef cold
#undef port
#undef RES_DE
#undef TT
#undef FRAC_TUNNEL
#undef FRAC_QUASIEL
#undef V_HOLE
#undef Etun
#undef Gamma
#undef FAST_THETA
#undef THETA_fast1
#undef FO3
#undef SLOW1_THETA
#undef N_fast
#undef F_fast2
#undef F_fast1
#undef F_slow2
#undef F_slow1
#undef M
#undef ALPHA
#undef W_end
#undef H_end
#undef W_chop
#undef H_chop
#undef W4
#undef H4
#undef W3
#undef H3
#undef W2
#undef H2
#undef W1
#undef H1
#undef GUI_GAP
#undef GUI_w
#undef GUI_h
#undef DETECTOR_DIST
#undef SAMPLE_DIST
#undef Length
#undef L_ballistic_end
#undef L_ballistic_begin
#undef FO1_DIST
#undef GUI_start
#undef Num_pulses
#undef Pulse_width
#undef lambda0
#undef Lmax
#undef Lmin
#undef mcposaESS_IN5_reprate
#undef mccompcurindex
#undef mccompcurtype
#undef mccompcurname

/* neutron state table at each component input (local coords) */
/* [x, y, z, vx, vy, vz, t, sx, sy, sz, p] */
MCNUM mccomp_storein[11*52];
/* Components position table (absolute and relative coords) */
Coords mccomp_posa[52];
Coords mccomp_posr[52];
/* Counter for each comp to check for inactive ones */
MCNUM  mcNCounter[52];
MCNUM  mcPCounter[52];
MCNUM  mcP2Counter[52];
#define mcNUMCOMP 51 /* number of components */
/* Counter for PROP ABSORB */
MCNUM  mcAbsorbProp[52];
/* Flag true when previous component acted on the neutron (SCATTER) */
MCNUM mcScattered=0;
/* Flag true when neutron should be restored (RESTORE) */
MCNUM mcRestore=0;
/* Declarations of component definition and setting parameters. */

/* Setting parameters for component 'source' [1]. */
MCNUM mccsource_width_c;
MCNUM mccsource_yheight;
MCNUM mccsource_Lmin;
MCNUM mccsource_Lmax;
MCNUM mccsource_dist;
MCNUM mccsource_focus_xw;
MCNUM mccsource_focus_yh;
MCNUM mccsource_nu;
MCNUM mccsource_T;
MCNUM mccsource_tau;
MCNUM mccsource_tau1;
MCNUM mccsource_tau2;
MCNUM mccsource_d;
MCNUM mccsource_n;
MCNUM mccsource_cold_frac;
MCNUM mccsource_n2;
MCNUM mccsource_chi2;
MCNUM mccsource_I0;
MCNUM mccsource_I2;
int mccsource_target_index;
MCNUM mccsource_cyl_radius;
MCNUM mccsource_branch1;
MCNUM mccsource_branch2;
MCNUM mccsource_branch_tail;
int mccsource_n_pulses;
MCNUM mccsource_width_t;
MCNUM mccsource_T_t;
MCNUM mccsource_tau_t;
MCNUM mccsource_tau1_t;
MCNUM mccsource_tau2_t;
MCNUM mccsource_chi2_t;
MCNUM mccsource_I0_t;
MCNUM mccsource_I2_t;
MCNUM mccsource_branch1_t;
MCNUM mccsource_branch2_t;
int mccsource_src_2012;
MCNUM mccsource_tfocus_dist;
MCNUM mccsource_tfocus_time;
MCNUM mccsource_tfocus_width;
MCNUM mccsource_beamport_angle;

/* Setting parameters for component 'Origin' [2]. */
char mccOrigin_profile[16384];
MCNUM mccOrigin_percent;
MCNUM mccOrigin_flag_save;
MCNUM mccOrigin_minutes;

/* Definition parameters for component 'TOFmoderator_zoom' [3]. */
#define mccTOFmoderator_zoom_nt 1000
/* Setting parameters for component 'TOFmoderator_zoom' [3]. */
char mccTOFmoderator_zoom_filename[16384];
MCNUM mccTOFmoderator_zoom_xmin;
MCNUM mccTOFmoderator_zoom_xmax;
MCNUM mccTOFmoderator_zoom_ymin;
MCNUM mccTOFmoderator_zoom_ymax;
MCNUM mccTOFmoderator_zoom_xwidth;
MCNUM mccTOFmoderator_zoom_yheight;
MCNUM mccTOFmoderator_zoom_tmin;
MCNUM mccTOFmoderator_zoom_tmax;
MCNUM mccTOFmoderator_zoom_dt;
MCNUM mccTOFmoderator_zoom_restore_neutron;
int mccTOFmoderator_zoom_nowritefile;

/* Definition parameters for component 'TOFmoderator' [4]. */
#define mccTOFmoderator_nt 1000
/* Setting parameters for component 'TOFmoderator' [4]. */
char mccTOFmoderator_filename[16384];
MCNUM mccTOFmoderator_xmin;
MCNUM mccTOFmoderator_xmax;
MCNUM mccTOFmoderator_ymin;
MCNUM mccTOFmoderator_ymax;
MCNUM mccTOFmoderator_xwidth;
MCNUM mccTOFmoderator_yheight;
MCNUM mccTOFmoderator_tmin;
MCNUM mccTOFmoderator_tmax;
MCNUM mccTOFmoderator_dt;
MCNUM mccTOFmoderator_restore_neutron;
int mccTOFmoderator_nowritefile;

/* Definition parameters for component 'Lmon_guistart' [5]. */
#define mccLmon_guistart_nL 1000
/* Setting parameters for component 'Lmon_guistart' [5]. */
char mccLmon_guistart_filename[16384];
MCNUM mccLmon_guistart_xmin;
MCNUM mccLmon_guistart_xmax;
MCNUM mccLmon_guistart_ymin;
MCNUM mccLmon_guistart_ymax;
MCNUM mccLmon_guistart_xwidth;
MCNUM mccLmon_guistart_yheight;
MCNUM mccLmon_guistart_Lmin;
MCNUM mccLmon_guistart_Lmax;
MCNUM mccLmon_guistart_restore_neutron;
int mccLmon_guistart_nowritefile;

/* Definition parameters for component 'Lmon_normalize' [6]. */
#define mccLmon_normalize_nL 2880
/* Setting parameters for component 'Lmon_normalize' [6]. */
char mccLmon_normalize_filename[16384];
MCNUM mccLmon_normalize_xmin;
MCNUM mccLmon_normalize_xmax;
MCNUM mccLmon_normalize_ymin;
MCNUM mccLmon_normalize_ymax;
MCNUM mccLmon_normalize_xwidth;
MCNUM mccLmon_normalize_yheight;
MCNUM mccLmon_normalize_Lmin;
MCNUM mccLmon_normalize_Lmax;
MCNUM mccLmon_normalize_restore_neutron;
int mccLmon_normalize_nowritefile;

/* Setting parameters for component 'Guide1' [7]. */
char mccGuide1_reflect[16384];
MCNUM mccGuide1_w1;
MCNUM mccGuide1_h1;
MCNUM mccGuide1_w2;
MCNUM mccGuide1_h2;
MCNUM mccGuide1_l;
MCNUM mccGuide1_R0;
MCNUM mccGuide1_Qc;
MCNUM mccGuide1_alpha;
MCNUM mccGuide1_m;
MCNUM mccGuide1_W;

/* Definition parameters for component 'Lmonslow1' [8]. */
#define mccLmonslow1_nL 200
/* Setting parameters for component 'Lmonslow1' [8]. */
char mccLmonslow1_filename[16384];
MCNUM mccLmonslow1_xmin;
MCNUM mccLmonslow1_xmax;
MCNUM mccLmonslow1_ymin;
MCNUM mccLmonslow1_ymax;
MCNUM mccLmonslow1_xwidth;
MCNUM mccLmonslow1_yheight;
MCNUM mccLmonslow1_Lmin;
MCNUM mccLmonslow1_Lmax;
MCNUM mccLmonslow1_restore_neutron;
int mccLmonslow1_nowritefile;

/* Setting parameters for component 'PSDslow1' [9]. */
int mccPSDslow1_nx;
int mccPSDslow1_ny;
char mccPSDslow1_filename[16384];
MCNUM mccPSDslow1_xmin;
MCNUM mccPSDslow1_xmax;
MCNUM mccPSDslow1_ymin;
MCNUM mccPSDslow1_ymax;
MCNUM mccPSDslow1_xwidth;
MCNUM mccPSDslow1_yheight;
MCNUM mccPSDslow1_restore_neutron;

/* Setting parameters for component 'FOchop1' [10]. */
MCNUM mccFOchop1_theta_0;
MCNUM mccFOchop1_radius;
MCNUM mccFOchop1_yheight;
MCNUM mccFOchop1_nu;
MCNUM mccFOchop1_nslit;
MCNUM mccFOchop1_jitter;
MCNUM mccFOchop1_delay;
MCNUM mccFOchop1_isfirst;
MCNUM mccFOchop1_n_pulse;
MCNUM mccFOchop1_abs_out;
MCNUM mccFOchop1_phase;
MCNUM mccFOchop1_xwidth;
MCNUM mccFOchop1_verbose;

/* Definition parameters for component 'TOFLmon1' [11]. */
#define mccTOFLmon1_nL 200
#define mccTOFLmon1_nt 200
#define mccTOFLmon1_tmin 0
#define mccTOFLmon1_tmax 3e5
/* Setting parameters for component 'TOFLmon1' [11]. */
char mccTOFLmon1_filename[16384];
MCNUM mccTOFLmon1_xmin;
MCNUM mccTOFLmon1_xmax;
MCNUM mccTOFLmon1_ymin;
MCNUM mccTOFLmon1_ymax;
MCNUM mccTOFLmon1_xwidth;
MCNUM mccTOFLmon1_yheight;
MCNUM mccTOFLmon1_Lmin;
MCNUM mccTOFLmon1_Lmax;
MCNUM mccTOFLmon1_restore_neutron;
int mccTOFLmon1_nowritefile;

/* Definition parameters for component 'Lmon_afterslow1' [12]. */
#define mccLmon_afterslow1_nL 200
/* Setting parameters for component 'Lmon_afterslow1' [12]. */
char mccLmon_afterslow1_filename[16384];
MCNUM mccLmon_afterslow1_xmin;
MCNUM mccLmon_afterslow1_xmax;
MCNUM mccLmon_afterslow1_ymin;
MCNUM mccLmon_afterslow1_ymax;
MCNUM mccLmon_afterslow1_xwidth;
MCNUM mccLmon_afterslow1_yheight;
MCNUM mccLmon_afterslow1_Lmin;
MCNUM mccLmon_afterslow1_Lmax;
MCNUM mccLmon_afterslow1_restore_neutron;
int mccLmon_afterslow1_nowritefile;

/* Setting parameters for component 'PSD_afterslow1' [13]. */
int mccPSD_afterslow1_nx;
int mccPSD_afterslow1_ny;
char mccPSD_afterslow1_filename[16384];
MCNUM mccPSD_afterslow1_xmin;
MCNUM mccPSD_afterslow1_xmax;
MCNUM mccPSD_afterslow1_ymin;
MCNUM mccPSD_afterslow1_ymax;
MCNUM mccPSD_afterslow1_xwidth;
MCNUM mccPSD_afterslow1_yheight;
MCNUM mccPSD_afterslow1_restore_neutron;

/* Setting parameters for component 'Guidelong1' [14]. */
char mccGuidelong1_reflect[16384];
MCNUM mccGuidelong1_w1;
MCNUM mccGuidelong1_h1;
MCNUM mccGuidelong1_w2;
MCNUM mccGuidelong1_h2;
MCNUM mccGuidelong1_l;
MCNUM mccGuidelong1_R0;
MCNUM mccGuidelong1_Qc;
MCNUM mccGuidelong1_alpha;
MCNUM mccGuidelong1_m;
MCNUM mccGuidelong1_W;

/* Setting parameters for component 'Guidelong1b' [15]. */
char mccGuidelong1b_reflect[16384];
MCNUM mccGuidelong1b_w1;
MCNUM mccGuidelong1b_h1;
MCNUM mccGuidelong1b_w2;
MCNUM mccGuidelong1b_h2;
MCNUM mccGuidelong1b_l;
MCNUM mccGuidelong1b_R0;
MCNUM mccGuidelong1b_Qc;
MCNUM mccGuidelong1b_alpha;
MCNUM mccGuidelong1b_m;
MCNUM mccGuidelong1b_W;

/* Definition parameters for component 'Lmon_slow2' [16]. */
#define mccLmon_slow2_nL 200
/* Setting parameters for component 'Lmon_slow2' [16]. */
char mccLmon_slow2_filename[16384];
MCNUM mccLmon_slow2_xmin;
MCNUM mccLmon_slow2_xmax;
MCNUM mccLmon_slow2_ymin;
MCNUM mccLmon_slow2_ymax;
MCNUM mccLmon_slow2_xwidth;
MCNUM mccLmon_slow2_yheight;
MCNUM mccLmon_slow2_Lmin;
MCNUM mccLmon_slow2_Lmax;
MCNUM mccLmon_slow2_restore_neutron;
int mccLmon_slow2_nowritefile;

/* Setting parameters for component 'FOchop2' [17]. */
MCNUM mccFOchop2_theta_0;
MCNUM mccFOchop2_radius;
MCNUM mccFOchop2_yheight;
MCNUM mccFOchop2_nu;
MCNUM mccFOchop2_nslit;
MCNUM mccFOchop2_jitter;
MCNUM mccFOchop2_delay;
MCNUM mccFOchop2_isfirst;
MCNUM mccFOchop2_n_pulse;
MCNUM mccFOchop2_abs_out;
MCNUM mccFOchop2_phase;
MCNUM mccFOchop2_xwidth;
MCNUM mccFOchop2_verbose;

/* Setting parameters for component 'Fastchop1' [18]. */
MCNUM mccFastchop1_theta_0;
MCNUM mccFastchop1_radius;
MCNUM mccFastchop1_yheight;
MCNUM mccFastchop1_nu;
MCNUM mccFastchop1_nslit;
MCNUM mccFastchop1_jitter;
MCNUM mccFastchop1_delay;
MCNUM mccFastchop1_isfirst;
MCNUM mccFastchop1_n_pulse;
MCNUM mccFastchop1_abs_out;
MCNUM mccFastchop1_phase;
MCNUM mccFastchop1_xwidth;
MCNUM mccFastchop1_verbose;

/* Setting parameters for component 'PSD_afterslow2' [19]. */
int mccPSD_afterslow2_nx;
int mccPSD_afterslow2_ny;
char mccPSD_afterslow2_filename[16384];
MCNUM mccPSD_afterslow2_xmin;
MCNUM mccPSD_afterslow2_xmax;
MCNUM mccPSD_afterslow2_ymin;
MCNUM mccPSD_afterslow2_ymax;
MCNUM mccPSD_afterslow2_xwidth;
MCNUM mccPSD_afterslow2_yheight;
MCNUM mccPSD_afterslow2_restore_neutron;

/* Definition parameters for component 'Lmon_afterslow2' [20]. */
#define mccLmon_afterslow2_nL 200
/* Setting parameters for component 'Lmon_afterslow2' [20]. */
char mccLmon_afterslow2_filename[16384];
MCNUM mccLmon_afterslow2_xmin;
MCNUM mccLmon_afterslow2_xmax;
MCNUM mccLmon_afterslow2_ymin;
MCNUM mccLmon_afterslow2_ymax;
MCNUM mccLmon_afterslow2_xwidth;
MCNUM mccLmon_afterslow2_yheight;
MCNUM mccLmon_afterslow2_Lmin;
MCNUM mccLmon_afterslow2_Lmax;
MCNUM mccLmon_afterslow2_restore_neutron;
int mccLmon_afterslow2_nowritefile;

/* Definition parameters for component 'TOFL_afterslow2' [21]. */
#define mccTOFL_afterslow2_nL 200
#define mccTOFL_afterslow2_nt 200
#define mccTOFL_afterslow2_tmin 0
#define mccTOFL_afterslow2_tmax 3.0e5
/* Setting parameters for component 'TOFL_afterslow2' [21]. */
char mccTOFL_afterslow2_filename[16384];
MCNUM mccTOFL_afterslow2_xmin;
MCNUM mccTOFL_afterslow2_xmax;
MCNUM mccTOFL_afterslow2_ymin;
MCNUM mccTOFL_afterslow2_ymax;
MCNUM mccTOFL_afterslow2_xwidth;
MCNUM mccTOFL_afterslow2_yheight;
MCNUM mccTOFL_afterslow2_Lmin;
MCNUM mccTOFL_afterslow2_Lmax;
MCNUM mccTOFL_afterslow2_restore_neutron;
int mccTOFL_afterslow2_nowritefile;

/* Setting parameters for component 'Guidelong2' [22]. */
char mccGuidelong2_reflect[16384];
MCNUM mccGuidelong2_w1;
MCNUM mccGuidelong2_h1;
MCNUM mccGuidelong2_w2;
MCNUM mccGuidelong2_h2;
MCNUM mccGuidelong2_l;
MCNUM mccGuidelong2_R0;
MCNUM mccGuidelong2_Qc;
MCNUM mccGuidelong2_alpha;
MCNUM mccGuidelong2_m;
MCNUM mccGuidelong2_W;

/* Definition parameters for component 'Lmon_beforeballistic' [23]. */
#define mccLmon_beforeballistic_nL 200
/* Setting parameters for component 'Lmon_beforeballistic' [23]. */
char mccLmon_beforeballistic_filename[16384];
MCNUM mccLmon_beforeballistic_xmin;
MCNUM mccLmon_beforeballistic_xmax;
MCNUM mccLmon_beforeballistic_ymin;
MCNUM mccLmon_beforeballistic_ymax;
MCNUM mccLmon_beforeballistic_xwidth;
MCNUM mccLmon_beforeballistic_yheight;
MCNUM mccLmon_beforeballistic_Lmin;
MCNUM mccLmon_beforeballistic_Lmax;
MCNUM mccLmon_beforeballistic_restore_neutron;
int mccLmon_beforeballistic_nowritefile;

/* Setting parameters for component 'PSD_beforeballistic' [24]. */
int mccPSD_beforeballistic_nx;
int mccPSD_beforeballistic_ny;
char mccPSD_beforeballistic_filename[16384];
MCNUM mccPSD_beforeballistic_xmin;
MCNUM mccPSD_beforeballistic_xmax;
MCNUM mccPSD_beforeballistic_ymin;
MCNUM mccPSD_beforeballistic_ymax;
MCNUM mccPSD_beforeballistic_xwidth;
MCNUM mccPSD_beforeballistic_yheight;
MCNUM mccPSD_beforeballistic_restore_neutron;

/* Setting parameters for component 'Guidelong2a' [25]. */
char mccGuidelong2a_reflect[16384];
MCNUM mccGuidelong2a_w1;
MCNUM mccGuidelong2a_h1;
MCNUM mccGuidelong2a_w2;
MCNUM mccGuidelong2a_h2;
MCNUM mccGuidelong2a_l;
MCNUM mccGuidelong2a_R0;
MCNUM mccGuidelong2a_Qc;
MCNUM mccGuidelong2a_alpha;
MCNUM mccGuidelong2a_m;
MCNUM mccGuidelong2a_W;

/* Definition parameters for component 'Lmonfast2' [26]. */
#define mccLmonfast2_nL 200
/* Setting parameters for component 'Lmonfast2' [26]. */
char mccLmonfast2_filename[16384];
MCNUM mccLmonfast2_xmin;
MCNUM mccLmonfast2_xmax;
MCNUM mccLmonfast2_ymin;
MCNUM mccLmonfast2_ymax;
MCNUM mccLmonfast2_xwidth;
MCNUM mccLmonfast2_yheight;
MCNUM mccLmonfast2_Lmin;
MCNUM mccLmonfast2_Lmax;
MCNUM mccLmonfast2_restore_neutron;
int mccLmonfast2_nowritefile;

/* Definition parameters for component 'Lmonfast2_zoom' [27]. */
#define mccLmonfast2_zoom_nL 200
/* Setting parameters for component 'Lmonfast2_zoom' [27]. */
char mccLmonfast2_zoom_filename[16384];
MCNUM mccLmonfast2_zoom_xmin;
MCNUM mccLmonfast2_zoom_xmax;
MCNUM mccLmonfast2_zoom_ymin;
MCNUM mccLmonfast2_zoom_ymax;
MCNUM mccLmonfast2_zoom_xwidth;
MCNUM mccLmonfast2_zoom_yheight;
MCNUM mccLmonfast2_zoom_Lmin;
MCNUM mccLmonfast2_zoom_Lmax;
MCNUM mccLmonfast2_zoom_restore_neutron;
int mccLmonfast2_zoom_nowritefile;

/* Definition parameters for component 'TOFLfast2' [28]. */
#define mccTOFLfast2_nL 200
#define mccTOFLfast2_nt 200
#define mccTOFLfast2_tmin 0
#define mccTOFLfast2_tmax 3.0e5
/* Setting parameters for component 'TOFLfast2' [28]. */
char mccTOFLfast2_filename[16384];
MCNUM mccTOFLfast2_xmin;
MCNUM mccTOFLfast2_xmax;
MCNUM mccTOFLfast2_ymin;
MCNUM mccTOFLfast2_ymax;
MCNUM mccTOFLfast2_xwidth;
MCNUM mccTOFLfast2_yheight;
MCNUM mccTOFLfast2_Lmin;
MCNUM mccTOFLfast2_Lmax;
MCNUM mccTOFLfast2_restore_neutron;
int mccTOFLfast2_nowritefile;

/* Definition parameters for component 'TOFLfast2zoom' [29]. */
#define mccTOFLfast2zoom_nL 200
#define mccTOFLfast2zoom_nt 200
#define mccTOFLfast2zoom_tmin tmin_zoom
#define mccTOFLfast2zoom_tmax tmax_zoom
/* Setting parameters for component 'TOFLfast2zoom' [29]. */
char mccTOFLfast2zoom_filename[16384];
MCNUM mccTOFLfast2zoom_xmin;
MCNUM mccTOFLfast2zoom_xmax;
MCNUM mccTOFLfast2zoom_ymin;
MCNUM mccTOFLfast2zoom_ymax;
MCNUM mccTOFLfast2zoom_xwidth;
MCNUM mccTOFLfast2zoom_yheight;
MCNUM mccTOFLfast2zoom_Lmin;
MCNUM mccTOFLfast2zoom_Lmax;
MCNUM mccTOFLfast2zoom_restore_neutron;
int mccTOFLfast2zoom_nowritefile;

/* Setting parameters for component 'PSDfast2' [30]. */
int mccPSDfast2_nx;
int mccPSDfast2_ny;
char mccPSDfast2_filename[16384];
MCNUM mccPSDfast2_xmin;
MCNUM mccPSDfast2_xmax;
MCNUM mccPSDfast2_ymin;
MCNUM mccPSDfast2_ymax;
MCNUM mccPSDfast2_xwidth;
MCNUM mccPSDfast2_yheight;
MCNUM mccPSDfast2_restore_neutron;

/* Setting parameters for component 'Fastchop2' [31]. */
MCNUM mccFastchop2_theta_0;
MCNUM mccFastchop2_radius;
MCNUM mccFastchop2_yheight;
MCNUM mccFastchop2_nu;
MCNUM mccFastchop2_nslit;
MCNUM mccFastchop2_jitter;
MCNUM mccFastchop2_delay;
MCNUM mccFastchop2_isfirst;
MCNUM mccFastchop2_n_pulse;
MCNUM mccFastchop2_abs_out;
MCNUM mccFastchop2_phase;
MCNUM mccFastchop2_xwidth;
MCNUM mccFastchop2_verbose;

/* Setting parameters for component 'Fastchop2counter' [32]. */
MCNUM mccFastchop2counter_theta_0;
MCNUM mccFastchop2counter_radius;
MCNUM mccFastchop2counter_yheight;
MCNUM mccFastchop2counter_nu;
MCNUM mccFastchop2counter_nslit;
MCNUM mccFastchop2counter_jitter;
MCNUM mccFastchop2counter_delay;
MCNUM mccFastchop2counter_isfirst;
MCNUM mccFastchop2counter_n_pulse;
MCNUM mccFastchop2counter_abs_out;
MCNUM mccFastchop2counter_phase;
MCNUM mccFastchop2counter_xwidth;
MCNUM mccFastchop2counter_verbose;

/* Setting parameters for component 'FOchop3' [33]. */
MCNUM mccFOchop3_theta_0;
MCNUM mccFOchop3_radius;
MCNUM mccFOchop3_yheight;
MCNUM mccFOchop3_nu;
MCNUM mccFOchop3_nslit;
MCNUM mccFOchop3_jitter;
MCNUM mccFOchop3_delay;
MCNUM mccFOchop3_isfirst;
MCNUM mccFOchop3_n_pulse;
MCNUM mccFOchop3_abs_out;
MCNUM mccFOchop3_phase;
MCNUM mccFOchop3_xwidth;
MCNUM mccFOchop3_verbose;

/* Definition parameters for component 'TOFfast2_zoom' [34]. */
#define mccTOFfast2_zoom_nt 100
/* Setting parameters for component 'TOFfast2_zoom' [34]. */
char mccTOFfast2_zoom_filename[16384];
MCNUM mccTOFfast2_zoom_xmin;
MCNUM mccTOFfast2_zoom_xmax;
MCNUM mccTOFfast2_zoom_ymin;
MCNUM mccTOFfast2_zoom_ymax;
MCNUM mccTOFfast2_zoom_xwidth;
MCNUM mccTOFfast2_zoom_yheight;
MCNUM mccTOFfast2_zoom_tmin;
MCNUM mccTOFfast2_zoom_tmax;
MCNUM mccTOFfast2_zoom_dt;
MCNUM mccTOFfast2_zoom_restore_neutron;
int mccTOFfast2_zoom_nowritefile;

/* Definition parameters for component 'Lmon_afterfast2' [35]. */
#define mccLmon_afterfast2_nL 500
/* Setting parameters for component 'Lmon_afterfast2' [35]. */
char mccLmon_afterfast2_filename[16384];
MCNUM mccLmon_afterfast2_xmin;
MCNUM mccLmon_afterfast2_xmax;
MCNUM mccLmon_afterfast2_ymin;
MCNUM mccLmon_afterfast2_ymax;
MCNUM mccLmon_afterfast2_xwidth;
MCNUM mccLmon_afterfast2_yheight;
MCNUM mccLmon_afterfast2_Lmin;
MCNUM mccLmon_afterfast2_Lmax;
MCNUM mccLmon_afterfast2_restore_neutron;
int mccLmon_afterfast2_nowritefile;

/* Definition parameters for component 'TOFL_afterfast2' [36]. */
#define mccTOFL_afterfast2_nL 200
#define mccTOFL_afterfast2_nt 200
#define mccTOFL_afterfast2_tmin 0
#define mccTOFL_afterfast2_tmax 3.0e5
/* Setting parameters for component 'TOFL_afterfast2' [36]. */
char mccTOFL_afterfast2_filename[16384];
MCNUM mccTOFL_afterfast2_xmin;
MCNUM mccTOFL_afterfast2_xmax;
MCNUM mccTOFL_afterfast2_ymin;
MCNUM mccTOFL_afterfast2_ymax;
MCNUM mccTOFL_afterfast2_xwidth;
MCNUM mccTOFL_afterfast2_yheight;
MCNUM mccTOFL_afterfast2_Lmin;
MCNUM mccTOFL_afterfast2_Lmax;
MCNUM mccTOFL_afterfast2_restore_neutron;
int mccTOFL_afterfast2_nowritefile;

/* Definition parameters for component 'TOFL_afterfast2_zoom' [37]. */
#define mccTOFL_afterfast2_zoom_nL 200
#define mccTOFL_afterfast2_zoom_nt 200
#define mccTOFL_afterfast2_zoom_tmin tmin_zoom
#define mccTOFL_afterfast2_zoom_tmax tmax_zoom
/* Setting parameters for component 'TOFL_afterfast2_zoom' [37]. */
char mccTOFL_afterfast2_zoom_filename[16384];
MCNUM mccTOFL_afterfast2_zoom_xmin;
MCNUM mccTOFL_afterfast2_zoom_xmax;
MCNUM mccTOFL_afterfast2_zoom_ymin;
MCNUM mccTOFL_afterfast2_zoom_ymax;
MCNUM mccTOFL_afterfast2_zoom_xwidth;
MCNUM mccTOFL_afterfast2_zoom_yheight;
MCNUM mccTOFL_afterfast2_zoom_Lmin;
MCNUM mccTOFL_afterfast2_zoom_Lmax;
MCNUM mccTOFL_afterfast2_zoom_restore_neutron;
int mccTOFL_afterfast2_zoom_nowritefile;

/* Setting parameters for component 'PSD_afterfast2' [38]. */
int mccPSD_afterfast2_nx;
int mccPSD_afterfast2_ny;
char mccPSD_afterfast2_filename[16384];
MCNUM mccPSD_afterfast2_xmin;
MCNUM mccPSD_afterfast2_xmax;
MCNUM mccPSD_afterfast2_ymin;
MCNUM mccPSD_afterfast2_ymax;
MCNUM mccPSD_afterfast2_xwidth;
MCNUM mccPSD_afterfast2_yheight;
MCNUM mccPSD_afterfast2_restore_neutron;

/* Setting parameters for component 'Guidesample' [39]. */
char mccGuidesample_reflect[16384];
MCNUM mccGuidesample_w1;
MCNUM mccGuidesample_h1;
MCNUM mccGuidesample_w2;
MCNUM mccGuidesample_h2;
MCNUM mccGuidesample_l;
MCNUM mccGuidesample_R0;
MCNUM mccGuidesample_Qc;
MCNUM mccGuidesample_alpha;
MCNUM mccGuidesample_m;
MCNUM mccGuidesample_W;

/* Definition parameters for component 'Lmon_guideend' [40]. */
#define mccLmon_guideend_nL 1000
/* Setting parameters for component 'Lmon_guideend' [40]. */
char mccLmon_guideend_filename[16384];
MCNUM mccLmon_guideend_xmin;
MCNUM mccLmon_guideend_xmax;
MCNUM mccLmon_guideend_ymin;
MCNUM mccLmon_guideend_ymax;
MCNUM mccLmon_guideend_xwidth;
MCNUM mccLmon_guideend_yheight;
MCNUM mccLmon_guideend_Lmin;
MCNUM mccLmon_guideend_Lmax;
MCNUM mccLmon_guideend_restore_neutron;
int mccLmon_guideend_nowritefile;

/* Setting parameters for component 'PSDsample' [41]. */
int mccPSDsample_nx;
int mccPSDsample_ny;
char mccPSDsample_filename[16384];
MCNUM mccPSDsample_xmin;
MCNUM mccPSDsample_xmax;
MCNUM mccPSDsample_ymin;
MCNUM mccPSDsample_ymax;
MCNUM mccPSDsample_xwidth;
MCNUM mccPSDsample_yheight;
MCNUM mccPSDsample_restore_neutron;

/* Definition parameters for component 'TOFsample_zoom' [42]. */
#define mccTOFsample_zoom_nt 500
/* Setting parameters for component 'TOFsample_zoom' [42]. */
char mccTOFsample_zoom_filename[16384];
MCNUM mccTOFsample_zoom_xmin;
MCNUM mccTOFsample_zoom_xmax;
MCNUM mccTOFsample_zoom_ymin;
MCNUM mccTOFsample_zoom_ymax;
MCNUM mccTOFsample_zoom_xwidth;
MCNUM mccTOFsample_zoom_yheight;
MCNUM mccTOFsample_zoom_tmin;
MCNUM mccTOFsample_zoom_tmax;
MCNUM mccTOFsample_zoom_dt;
MCNUM mccTOFsample_zoom_restore_neutron;
int mccTOFsample_zoom_nowritefile;

/* Definition parameters for component 'Esample' [43]. */
#define mccEsample_nE 400
/* Setting parameters for component 'Esample' [43]. */
char mccEsample_filename[16384];
MCNUM mccEsample_xmin;
MCNUM mccEsample_xmax;
MCNUM mccEsample_ymin;
MCNUM mccEsample_ymax;
MCNUM mccEsample_xwidth;
MCNUM mccEsample_yheight;
MCNUM mccEsample_Emin;
MCNUM mccEsample_Emax;
MCNUM mccEsample_restore_neutron;
int mccEsample_nowritefile;

/* Definition parameters for component 'Lmon_sample_zoom' [44]. */
#define mccLmon_sample_zoom_nL 100
/* Setting parameters for component 'Lmon_sample_zoom' [44]. */
char mccLmon_sample_zoom_filename[16384];
MCNUM mccLmon_sample_zoom_xmin;
MCNUM mccLmon_sample_zoom_xmax;
MCNUM mccLmon_sample_zoom_ymin;
MCNUM mccLmon_sample_zoom_ymax;
MCNUM mccLmon_sample_zoom_xwidth;
MCNUM mccLmon_sample_zoom_yheight;
MCNUM mccLmon_sample_zoom_Lmin;
MCNUM mccLmon_sample_zoom_Lmax;
MCNUM mccLmon_sample_zoom_restore_neutron;
int mccLmon_sample_zoom_nowritefile;

/* Setting parameters for component 'sample' [45]. */
MCNUM mccsample_thickness;
MCNUM mccsample_radius;
MCNUM mccsample_focus_r;
MCNUM mccsample_p_interact;
MCNUM mccsample_f_QE;
MCNUM mccsample_f_tun;
MCNUM mccsample_gamma;
MCNUM mccsample_E_tun;
MCNUM mccsample_target_x;
MCNUM mccsample_target_y;
MCNUM mccsample_target_z;
MCNUM mccsample_focus_xw;
MCNUM mccsample_focus_yh;
MCNUM mccsample_focus_aw;
MCNUM mccsample_focus_ah;
MCNUM mccsample_xwidth;
MCNUM mccsample_yheight;
MCNUM mccsample_zdepth;
MCNUM mccsample_sigma_abs;
MCNUM mccsample_sigma_inc;
MCNUM mccsample_Vc;
int mccsample_target_index;

/* Definition parameters for component 'TOFdetector' [47]. */
#define mccTOFdetector_nt 512
/* Setting parameters for component 'TOFdetector' [47]. */
char mccTOFdetector_filename[16384];
MCNUM mccTOFdetector_xmin;
MCNUM mccTOFdetector_xmax;
MCNUM mccTOFdetector_ymin;
MCNUM mccTOFdetector_ymax;
MCNUM mccTOFdetector_xwidth;
MCNUM mccTOFdetector_yheight;
MCNUM mccTOFdetector_tmin;
MCNUM mccTOFdetector_tmax;
MCNUM mccTOFdetector_dt;
MCNUM mccTOFdetector_restore_neutron;
int mccTOFdetector_nowritefile;

/* Definition parameters for component 'TOFdetector_zoom' [48]. */
#define mccTOFdetector_zoom_nt 100
/* Setting parameters for component 'TOFdetector_zoom' [48]. */
char mccTOFdetector_zoom_filename[16384];
MCNUM mccTOFdetector_zoom_xmin;
MCNUM mccTOFdetector_zoom_xmax;
MCNUM mccTOFdetector_zoom_ymin;
MCNUM mccTOFdetector_zoom_ymax;
MCNUM mccTOFdetector_zoom_xwidth;
MCNUM mccTOFdetector_zoom_yheight;
MCNUM mccTOFdetector_zoom_tmin;
MCNUM mccTOFdetector_zoom_tmax;
MCNUM mccTOFdetector_zoom_dt;
MCNUM mccTOFdetector_zoom_restore_neutron;
int mccTOFdetector_zoom_nowritefile;

/* Definition parameters for component 'Edetector' [49]. */
#define mccEdetector_nE 400
/* Setting parameters for component 'Edetector' [49]. */
char mccEdetector_filename[16384];
MCNUM mccEdetector_xmin;
MCNUM mccEdetector_xmax;
MCNUM mccEdetector_ymin;
MCNUM mccEdetector_ymax;
MCNUM mccEdetector_xwidth;
MCNUM mccEdetector_yheight;
MCNUM mccEdetector_Emin;
MCNUM mccEdetector_Emax;
MCNUM mccEdetector_restore_neutron;
int mccEdetector_nowritefile;

/* Definition parameters for component 'TOF2Edetector' [50]. */
#define mccTOF2Edetector_nE 200
/* Setting parameters for component 'TOF2Edetector' [50]. */
char mccTOF2Edetector_filename[16384];
MCNUM mccTOF2Edetector_xmin;
MCNUM mccTOF2Edetector_xmax;
MCNUM mccTOF2Edetector_ymin;
MCNUM mccTOF2Edetector_ymax;
MCNUM mccTOF2Edetector_xwidth;
MCNUM mccTOF2Edetector_yheight;
MCNUM mccTOF2Edetector_Emin;
MCNUM mccTOF2Edetector_Emax;
MCNUM mccTOF2Edetector_T_zero;
MCNUM mccTOF2Edetector_L_flight;
MCNUM mccTOF2Edetector_restore_neutron;
int mccTOF2Edetector_nowritefile;

/* User component declarations. */

/* User declarations for component 'source' [1]. */
#define mccompcurname  source
#define mccompcurtype  ESS_moderator_long
#define mccompcurindex 1
#define l_range mccsource_l_range
#define w_mult mccsource_w_mult
#define w_geom mccsource_w_geom
#define w_geom_c mccsource_w_geom_c
#define w_geom_t mccsource_w_geom_t
#define tx mccsource_tx
#define ty mccsource_ty
#define tz mccsource_tz
#define t1x mccsource_t1x
#define t1y mccsource_t1y
#define t1z mccsource_t1z
#define t2x mccsource_t2x
#define t2y mccsource_t2y
#define t2z mccsource_t2z
#define T_n mccsource_T_n
#define tau_n mccsource_tau_n
#define tau1_n mccsource_tau1_n
#define tau2_n mccsource_tau2_n
#define chi2_n mccsource_chi2_n
#define I0_n mccsource_I0_n
#define I2_n mccsource_I2_n
#define branch1_n mccsource_branch1_n
#define branch2_n mccsource_branch2_n
#define r_empty mccsource_r_empty
#define r_optics mccsource_r_optics
#define width_c mccsource_width_c
#define yheight mccsource_yheight
#define Lmin mccsource_Lmin
#define Lmax mccsource_Lmax
#define dist mccsource_dist
#define focus_xw mccsource_focus_xw
#define focus_yh mccsource_focus_yh
#define nu mccsource_nu
#define T mccsource_T
#define tau mccsource_tau
#define tau1 mccsource_tau1
#define tau2 mccsource_tau2
#define d mccsource_d
#define n mccsource_n
#define cold_frac mccsource_cold_frac
#define n2 mccsource_n2
#define chi2 mccsource_chi2
#define I0 mccsource_I0
#define I2 mccsource_I2
#define target_index mccsource_target_index
#define cyl_radius mccsource_cyl_radius
#define branch1 mccsource_branch1
#define branch2 mccsource_branch2
#define branch_tail mccsource_branch_tail
#define n_pulses mccsource_n_pulses
#define width_t mccsource_width_t
#define T_t mccsource_T_t
#define tau_t mccsource_tau_t
#define tau1_t mccsource_tau1_t
#define tau2_t mccsource_tau2_t
#define chi2_t mccsource_chi2_t
#define I0_t mccsource_I0_t
#define I2_t mccsource_I2_t
#define branch1_t mccsource_branch1_t
#define branch2_t mccsource_branch2_t
#define src_2012 mccsource_src_2012
#define tfocus_dist mccsource_tfocus_dist
#define tfocus_time mccsource_tfocus_time
#define tfocus_width mccsource_tfocus_width
#define beamport_angle mccsource_beamport_angle
#line 154 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../obsolete/ESS_moderator_long.comp"
  double l_range, w_mult, w_geom, w_geom_c, w_geom_t;
  double tx,ty,tz;
  double t1x,t1y,t1z,t2x,t2y,t2z;
  /* Neutron-specific distribution-shape variables */
  double T_n, tau_n, tau1_n, tau2_n, chi2_n, I0_n, I2_n, branch1_n, branch2_n;

  /* Target station geometry... */
  double r_empty; /* two meters from moderator surface and out... */
  double r_optics;
#line 8134 "./ESS_IN5_reprate.c"
#undef beamport_angle
#undef tfocus_width
#undef tfocus_time
#undef tfocus_dist
#undef src_2012
#undef branch2_t
#undef branch1_t
#undef I2_t
#undef I0_t
#undef chi2_t
#undef tau2_t
#undef tau1_t
#undef tau_t
#undef T_t
#undef width_t
#undef n_pulses
#undef branch_tail
#undef branch2
#undef branch1
#undef cyl_radius
#undef target_index
#undef I2
#undef I0
#undef chi2
#undef n2
#undef cold_frac
#undef n
#undef d
#undef tau2
#undef tau1
#undef tau
#undef T
#undef nu
#undef focus_yh
#undef focus_xw
#undef dist
#undef Lmax
#undef Lmin
#undef yheight
#undef width_c
#undef r_optics
#undef r_empty
#undef branch2_n
#undef branch1_n
#undef I2_n
#undef I0_n
#undef chi2_n
#undef tau2_n
#undef tau1_n
#undef tau_n
#undef T_n
#undef t2z
#undef t2y
#undef t2x
#undef t1z
#undef t1y
#undef t1x
#undef tz
#undef ty
#undef tx
#undef w_geom_t
#undef w_geom_c
#undef w_geom
#undef w_mult
#undef l_range
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'Origin' [2]. */
#define mccompcurname  Origin
#define mccompcurtype  Progress_bar
#define mccompcurindex 2
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
#line 8227 "./ESS_IN5_reprate.c"
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

/* User declarations for component 'TOFmoderator_zoom' [3]. */
#define mccompcurname  TOFmoderator_zoom
#define mccompcurtype  TOF_monitor
#define mccompcurindex 3
#define nt mccTOFmoderator_zoom_nt
#define TOF_N mccTOFmoderator_zoom_TOF_N
#define TOF_p mccTOFmoderator_zoom_TOF_p
#define TOF_p2 mccTOFmoderator_zoom_TOF_p2
#define t_min mccTOFmoderator_zoom_t_min
#define t_max mccTOFmoderator_zoom_t_max
#define delta_t mccTOFmoderator_zoom_delta_t
#define filename mccTOFmoderator_zoom_filename
#define xmin mccTOFmoderator_zoom_xmin
#define xmax mccTOFmoderator_zoom_xmax
#define ymin mccTOFmoderator_zoom_ymin
#define ymax mccTOFmoderator_zoom_ymax
#define xwidth mccTOFmoderator_zoom_xwidth
#define yheight mccTOFmoderator_zoom_yheight
#define tmin mccTOFmoderator_zoom_tmin
#define tmax mccTOFmoderator_zoom_tmax
#define dt mccTOFmoderator_zoom_dt
#define restore_neutron mccTOFmoderator_zoom_restore_neutron
#define nowritefile mccTOFmoderator_zoom_nowritefile
#line 54 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/TOF_monitor.comp"
double TOF_N[nt];
double TOF_p[nt];
double TOF_p2[nt];
double t_min, t_max, delta_t;
#line 8268 "./ESS_IN5_reprate.c"
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

/* User declarations for component 'TOFmoderator' [4]. */
#define mccompcurname  TOFmoderator
#define mccompcurtype  TOF_monitor
#define mccompcurindex 4
#define nt mccTOFmoderator_nt
#define TOF_N mccTOFmoderator_TOF_N
#define TOF_p mccTOFmoderator_TOF_p
#define TOF_p2 mccTOFmoderator_TOF_p2
#define t_min mccTOFmoderator_t_min
#define t_max mccTOFmoderator_t_max
#define delta_t mccTOFmoderator_delta_t
#define filename mccTOFmoderator_filename
#define xmin mccTOFmoderator_xmin
#define xmax mccTOFmoderator_xmax
#define ymin mccTOFmoderator_ymin
#define ymax mccTOFmoderator_ymax
#define xwidth mccTOFmoderator_xwidth
#define yheight mccTOFmoderator_yheight
#define tmin mccTOFmoderator_tmin
#define tmax mccTOFmoderator_tmax
#define dt mccTOFmoderator_dt
#define restore_neutron mccTOFmoderator_restore_neutron
#define nowritefile mccTOFmoderator_nowritefile
#line 54 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/TOF_monitor.comp"
double TOF_N[nt];
double TOF_p[nt];
double TOF_p2[nt];
double t_min, t_max, delta_t;
#line 8320 "./ESS_IN5_reprate.c"
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

/* User declarations for component 'Lmon_guistart' [5]. */
#define mccompcurname  Lmon_guistart
#define mccompcurtype  L_monitor
#define mccompcurindex 5
#define nL mccLmon_guistart_nL
#define L_N mccLmon_guistart_L_N
#define L_p mccLmon_guistart_L_p
#define L_p2 mccLmon_guistart_L_p2
#define filename mccLmon_guistart_filename
#define xmin mccLmon_guistart_xmin
#define xmax mccLmon_guistart_xmax
#define ymin mccLmon_guistart_ymin
#define ymax mccLmon_guistart_ymax
#define xwidth mccLmon_guistart_xwidth
#define yheight mccLmon_guistart_yheight
#define Lmin mccLmon_guistart_Lmin
#define Lmax mccLmon_guistart_Lmax
#define restore_neutron mccLmon_guistart_restore_neutron
#define nowritefile mccLmon_guistart_nowritefile
#line 57 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/L_monitor.comp"
double L_N[nL];
double L_p[nL], L_p2[nL];
#line 8366 "./ESS_IN5_reprate.c"
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

/* User declarations for component 'Lmon_normalize' [6]. */
#define mccompcurname  Lmon_normalize
#define mccompcurtype  L_monitor
#define mccompcurindex 6
#define nL mccLmon_normalize_nL
#define L_N mccLmon_normalize_L_N
#define L_p mccLmon_normalize_L_p
#define L_p2 mccLmon_normalize_L_p2
#define filename mccLmon_normalize_filename
#define xmin mccLmon_normalize_xmin
#define xmax mccLmon_normalize_xmax
#define ymin mccLmon_normalize_ymin
#define ymax mccLmon_normalize_ymax
#define xwidth mccLmon_normalize_xwidth
#define yheight mccLmon_normalize_yheight
#define Lmin mccLmon_normalize_Lmin
#define Lmax mccLmon_normalize_Lmax
#define restore_neutron mccLmon_normalize_restore_neutron
#define nowritefile mccLmon_normalize_nowritefile
#line 57 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/L_monitor.comp"
double L_N[nL];
double L_p[nL], L_p2[nL];
#line 8408 "./ESS_IN5_reprate.c"
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

/* User declarations for component 'Guide1' [7]. */
#define mccompcurname  Guide1
#define mccompcurtype  Guide
#define mccompcurindex 7
#define pTable mccGuide1_pTable
#define reflect mccGuide1_reflect
#define w1 mccGuide1_w1
#define h1 mccGuide1_h1
#define w2 mccGuide1_w2
#define h2 mccGuide1_h2
#define l mccGuide1_l
#define R0 mccGuide1_R0
#define Qc mccGuide1_Qc
#define alpha mccGuide1_alpha
#define m mccGuide1_m
#define W mccGuide1_W
#line 69 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide.comp"
t_Table pTable;
#line 8446 "./ESS_IN5_reprate.c"
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

/* User declarations for component 'Lmonslow1' [8]. */
#define mccompcurname  Lmonslow1
#define mccompcurtype  L_monitor
#define mccompcurindex 8
#define nL mccLmonslow1_nL
#define L_N mccLmonslow1_L_N
#define L_p mccLmonslow1_L_p
#define L_p2 mccLmonslow1_L_p2
#define filename mccLmonslow1_filename
#define xmin mccLmonslow1_xmin
#define xmax mccLmonslow1_xmax
#define ymin mccLmonslow1_ymin
#define ymax mccLmonslow1_ymax
#define xwidth mccLmonslow1_xwidth
#define yheight mccLmonslow1_yheight
#define Lmin mccLmonslow1_Lmin
#define Lmax mccLmonslow1_Lmax
#define restore_neutron mccLmonslow1_restore_neutron
#define nowritefile mccLmonslow1_nowritefile
#line 57 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/L_monitor.comp"
double L_N[nL];
double L_p[nL], L_p2[nL];
#line 8485 "./ESS_IN5_reprate.c"
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

/* User declarations for component 'PSDslow1' [9]. */
#define mccompcurname  PSDslow1
#define mccompcurtype  PSD_monitor
#define mccompcurindex 9
#define PSD_N mccPSDslow1_PSD_N
#define PSD_p mccPSDslow1_PSD_p
#define PSD_p2 mccPSDslow1_PSD_p2
#define nx mccPSDslow1_nx
#define ny mccPSDslow1_ny
#define filename mccPSDslow1_filename
#define xmin mccPSDslow1_xmin
#define xmax mccPSDslow1_xmax
#define ymin mccPSDslow1_ymin
#define ymax mccPSDslow1_ymax
#define xwidth mccPSDslow1_xwidth
#define yheight mccPSDslow1_yheight
#define restore_neutron mccPSDslow1_restore_neutron
#line 62 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
  DArray2d PSD_N;
  DArray2d PSD_p;
  DArray2d PSD_p2;
#line 8526 "./ESS_IN5_reprate.c"
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

/* User declarations for component 'FOchop1' [10]. */
#define mccompcurname  FOchop1
#define mccompcurtype  DiskChopper
#define mccompcurindex 10
#define Tg mccFOchop1_Tg
#define To mccFOchop1_To
#define delta_y mccFOchop1_delta_y
#define height mccFOchop1_height
#define omega mccFOchop1_omega
#define theta_0 mccFOchop1_theta_0
#define radius mccFOchop1_radius
#define yheight mccFOchop1_yheight
#define nu mccFOchop1_nu
#define nslit mccFOchop1_nslit
#define jitter mccFOchop1_jitter
#define delay mccFOchop1_delay
#define isfirst mccFOchop1_isfirst
#define n_pulse mccFOchop1_n_pulse
#define abs_out mccFOchop1_abs_out
#define phase mccFOchop1_phase
#define xwidth mccFOchop1_xwidth
#define verbose mccFOchop1_verbose
#line 66 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/DiskChopper.comp"
double Tg,To,delta_y,height,omega;
#line 8568 "./ESS_IN5_reprate.c"
#undef verbose
#undef xwidth
#undef phase
#undef abs_out
#undef n_pulse
#undef isfirst
#undef delay
#undef jitter
#undef nslit
#undef nu
#undef yheight
#undef radius
#undef theta_0
#undef omega
#undef height
#undef delta_y
#undef To
#undef Tg
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'TOFLmon1' [11]. */
#define mccompcurname  TOFLmon1
#define mccompcurtype  TOFLambda_monitor
#define mccompcurindex 11
#define nL mccTOFLmon1_nL
#define nt mccTOFLmon1_nt
#define tmin mccTOFLmon1_tmin
#define tmax mccTOFLmon1_tmax
#define tt_0 mccTOFLmon1_tt_0
#define tt_1 mccTOFLmon1_tt_1
#define TOFL_N mccTOFLmon1_TOFL_N
#define TOFL_p mccTOFLmon1_TOFL_p
#define TOFL_p2 mccTOFLmon1_TOFL_p2
#define filename mccTOFLmon1_filename
#define xmin mccTOFLmon1_xmin
#define xmax mccTOFLmon1_xmax
#define ymin mccTOFLmon1_ymin
#define ymax mccTOFLmon1_ymax
#define xwidth mccTOFLmon1_xwidth
#define yheight mccTOFLmon1_yheight
#define Lmin mccTOFLmon1_Lmin
#define Lmax mccTOFLmon1_Lmax
#define restore_neutron mccTOFLmon1_restore_neutron
#define nowritefile mccTOFLmon1_nowritefile
#line 55 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/TOFLambda_monitor.comp"
double TOFL_N[nt][nL];
double TOFL_p[nt][nL];
double TOFL_p2[nt][nL];
double tt_0, tt_1;
#line 8620 "./ESS_IN5_reprate.c"
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
#undef TOFL_p2
#undef TOFL_p
#undef TOFL_N
#undef tt_1
#undef tt_0
#undef tmax
#undef tmin
#undef nt
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'Lmon_afterslow1' [12]. */
#define mccompcurname  Lmon_afterslow1
#define mccompcurtype  L_monitor
#define mccompcurindex 12
#define nL mccLmon_afterslow1_nL
#define L_N mccLmon_afterslow1_L_N
#define L_p mccLmon_afterslow1_L_p
#define L_p2 mccLmon_afterslow1_L_p2
#define filename mccLmon_afterslow1_filename
#define xmin mccLmon_afterslow1_xmin
#define xmax mccLmon_afterslow1_xmax
#define ymin mccLmon_afterslow1_ymin
#define ymax mccLmon_afterslow1_ymax
#define xwidth mccLmon_afterslow1_xwidth
#define yheight mccLmon_afterslow1_yheight
#define Lmin mccLmon_afterslow1_Lmin
#define Lmax mccLmon_afterslow1_Lmax
#define restore_neutron mccLmon_afterslow1_restore_neutron
#define nowritefile mccLmon_afterslow1_nowritefile
#line 57 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/L_monitor.comp"
double L_N[nL];
double L_p[nL], L_p2[nL];
#line 8667 "./ESS_IN5_reprate.c"
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

/* User declarations for component 'PSD_afterslow1' [13]. */
#define mccompcurname  PSD_afterslow1
#define mccompcurtype  PSD_monitor
#define mccompcurindex 13
#define PSD_N mccPSD_afterslow1_PSD_N
#define PSD_p mccPSD_afterslow1_PSD_p
#define PSD_p2 mccPSD_afterslow1_PSD_p2
#define nx mccPSD_afterslow1_nx
#define ny mccPSD_afterslow1_ny
#define filename mccPSD_afterslow1_filename
#define xmin mccPSD_afterslow1_xmin
#define xmax mccPSD_afterslow1_xmax
#define ymin mccPSD_afterslow1_ymin
#define ymax mccPSD_afterslow1_ymax
#define xwidth mccPSD_afterslow1_xwidth
#define yheight mccPSD_afterslow1_yheight
#define restore_neutron mccPSD_afterslow1_restore_neutron
#line 62 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
  DArray2d PSD_N;
  DArray2d PSD_p;
  DArray2d PSD_p2;
#line 8708 "./ESS_IN5_reprate.c"
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

/* User declarations for component 'Guidelong1' [14]. */
#define mccompcurname  Guidelong1
#define mccompcurtype  Guide
#define mccompcurindex 14
#define pTable mccGuidelong1_pTable
#define reflect mccGuidelong1_reflect
#define w1 mccGuidelong1_w1
#define h1 mccGuidelong1_h1
#define w2 mccGuidelong1_w2
#define h2 mccGuidelong1_h2
#define l mccGuidelong1_l
#define R0 mccGuidelong1_R0
#define Qc mccGuidelong1_Qc
#define alpha mccGuidelong1_alpha
#define m mccGuidelong1_m
#define W mccGuidelong1_W
#line 69 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide.comp"
t_Table pTable;
#line 8744 "./ESS_IN5_reprate.c"
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

/* User declarations for component 'Guidelong1b' [15]. */
#define mccompcurname  Guidelong1b
#define mccompcurtype  Guide
#define mccompcurindex 15
#define pTable mccGuidelong1b_pTable
#define reflect mccGuidelong1b_reflect
#define w1 mccGuidelong1b_w1
#define h1 mccGuidelong1b_h1
#define w2 mccGuidelong1b_w2
#define h2 mccGuidelong1b_h2
#define l mccGuidelong1b_l
#define R0 mccGuidelong1b_R0
#define Qc mccGuidelong1b_Qc
#define alpha mccGuidelong1b_alpha
#define m mccGuidelong1b_m
#define W mccGuidelong1b_W
#line 69 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide.comp"
t_Table pTable;
#line 8779 "./ESS_IN5_reprate.c"
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

/* User declarations for component 'Lmon_slow2' [16]. */
#define mccompcurname  Lmon_slow2
#define mccompcurtype  L_monitor
#define mccompcurindex 16
#define nL mccLmon_slow2_nL
#define L_N mccLmon_slow2_L_N
#define L_p mccLmon_slow2_L_p
#define L_p2 mccLmon_slow2_L_p2
#define filename mccLmon_slow2_filename
#define xmin mccLmon_slow2_xmin
#define xmax mccLmon_slow2_xmax
#define ymin mccLmon_slow2_ymin
#define ymax mccLmon_slow2_ymax
#define xwidth mccLmon_slow2_xwidth
#define yheight mccLmon_slow2_yheight
#define Lmin mccLmon_slow2_Lmin
#define Lmax mccLmon_slow2_Lmax
#define restore_neutron mccLmon_slow2_restore_neutron
#define nowritefile mccLmon_slow2_nowritefile
#line 57 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/L_monitor.comp"
double L_N[nL];
double L_p[nL], L_p2[nL];
#line 8818 "./ESS_IN5_reprate.c"
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

/* User declarations for component 'FOchop2' [17]. */
#define mccompcurname  FOchop2
#define mccompcurtype  DiskChopper
#define mccompcurindex 17
#define Tg mccFOchop2_Tg
#define To mccFOchop2_To
#define delta_y mccFOchop2_delta_y
#define height mccFOchop2_height
#define omega mccFOchop2_omega
#define theta_0 mccFOchop2_theta_0
#define radius mccFOchop2_radius
#define yheight mccFOchop2_yheight
#define nu mccFOchop2_nu
#define nslit mccFOchop2_nslit
#define jitter mccFOchop2_jitter
#define delay mccFOchop2_delay
#define isfirst mccFOchop2_isfirst
#define n_pulse mccFOchop2_n_pulse
#define abs_out mccFOchop2_abs_out
#define phase mccFOchop2_phase
#define xwidth mccFOchop2_xwidth
#define verbose mccFOchop2_verbose
#line 66 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/DiskChopper.comp"
double Tg,To,delta_y,height,omega;
#line 8862 "./ESS_IN5_reprate.c"
#undef verbose
#undef xwidth
#undef phase
#undef abs_out
#undef n_pulse
#undef isfirst
#undef delay
#undef jitter
#undef nslit
#undef nu
#undef yheight
#undef radius
#undef theta_0
#undef omega
#undef height
#undef delta_y
#undef To
#undef Tg
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'Fastchop1' [18]. */
#define mccompcurname  Fastchop1
#define mccompcurtype  DiskChopper
#define mccompcurindex 18
#define Tg mccFastchop1_Tg
#define To mccFastchop1_To
#define delta_y mccFastchop1_delta_y
#define height mccFastchop1_height
#define omega mccFastchop1_omega
#define theta_0 mccFastchop1_theta_0
#define radius mccFastchop1_radius
#define yheight mccFastchop1_yheight
#define nu mccFastchop1_nu
#define nslit mccFastchop1_nslit
#define jitter mccFastchop1_jitter
#define delay mccFastchop1_delay
#define isfirst mccFastchop1_isfirst
#define n_pulse mccFastchop1_n_pulse
#define abs_out mccFastchop1_abs_out
#define phase mccFastchop1_phase
#define xwidth mccFastchop1_xwidth
#define verbose mccFastchop1_verbose
#line 66 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/DiskChopper.comp"
double Tg,To,delta_y,height,omega;
#line 8909 "./ESS_IN5_reprate.c"
#undef verbose
#undef xwidth
#undef phase
#undef abs_out
#undef n_pulse
#undef isfirst
#undef delay
#undef jitter
#undef nslit
#undef nu
#undef yheight
#undef radius
#undef theta_0
#undef omega
#undef height
#undef delta_y
#undef To
#undef Tg
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'PSD_afterslow2' [19]. */
#define mccompcurname  PSD_afterslow2
#define mccompcurtype  PSD_monitor
#define mccompcurindex 19
#define PSD_N mccPSD_afterslow2_PSD_N
#define PSD_p mccPSD_afterslow2_PSD_p
#define PSD_p2 mccPSD_afterslow2_PSD_p2
#define nx mccPSD_afterslow2_nx
#define ny mccPSD_afterslow2_ny
#define filename mccPSD_afterslow2_filename
#define xmin mccPSD_afterslow2_xmin
#define xmax mccPSD_afterslow2_xmax
#define ymin mccPSD_afterslow2_ymin
#define ymax mccPSD_afterslow2_ymax
#define xwidth mccPSD_afterslow2_xwidth
#define yheight mccPSD_afterslow2_yheight
#define restore_neutron mccPSD_afterslow2_restore_neutron
#line 62 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
  DArray2d PSD_N;
  DArray2d PSD_p;
  DArray2d PSD_p2;
#line 8953 "./ESS_IN5_reprate.c"
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

/* User declarations for component 'Lmon_afterslow2' [20]. */
#define mccompcurname  Lmon_afterslow2
#define mccompcurtype  L_monitor
#define mccompcurindex 20
#define nL mccLmon_afterslow2_nL
#define L_N mccLmon_afterslow2_L_N
#define L_p mccLmon_afterslow2_L_p
#define L_p2 mccLmon_afterslow2_L_p2
#define filename mccLmon_afterslow2_filename
#define xmin mccLmon_afterslow2_xmin
#define xmax mccLmon_afterslow2_xmax
#define ymin mccLmon_afterslow2_ymin
#define ymax mccLmon_afterslow2_ymax
#define xwidth mccLmon_afterslow2_xwidth
#define yheight mccLmon_afterslow2_yheight
#define Lmin mccLmon_afterslow2_Lmin
#define Lmax mccLmon_afterslow2_Lmax
#define restore_neutron mccLmon_afterslow2_restore_neutron
#define nowritefile mccLmon_afterslow2_nowritefile
#line 57 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/L_monitor.comp"
double L_N[nL];
double L_p[nL], L_p2[nL];
#line 8993 "./ESS_IN5_reprate.c"
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

/* User declarations for component 'TOFL_afterslow2' [21]. */
#define mccompcurname  TOFL_afterslow2
#define mccompcurtype  TOFLambda_monitor
#define mccompcurindex 21
#define nL mccTOFL_afterslow2_nL
#define nt mccTOFL_afterslow2_nt
#define tmin mccTOFL_afterslow2_tmin
#define tmax mccTOFL_afterslow2_tmax
#define tt_0 mccTOFL_afterslow2_tt_0
#define tt_1 mccTOFL_afterslow2_tt_1
#define TOFL_N mccTOFL_afterslow2_TOFL_N
#define TOFL_p mccTOFL_afterslow2_TOFL_p
#define TOFL_p2 mccTOFL_afterslow2_TOFL_p2
#define filename mccTOFL_afterslow2_filename
#define xmin mccTOFL_afterslow2_xmin
#define xmax mccTOFL_afterslow2_xmax
#define ymin mccTOFL_afterslow2_ymin
#define ymax mccTOFL_afterslow2_ymax
#define xwidth mccTOFL_afterslow2_xwidth
#define yheight mccTOFL_afterslow2_yheight
#define Lmin mccTOFL_afterslow2_Lmin
#define Lmax mccTOFL_afterslow2_Lmax
#define restore_neutron mccTOFL_afterslow2_restore_neutron
#define nowritefile mccTOFL_afterslow2_nowritefile
#line 55 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/TOFLambda_monitor.comp"
double TOFL_N[nt][nL];
double TOFL_p[nt][nL];
double TOFL_p2[nt][nL];
double tt_0, tt_1;
#line 9042 "./ESS_IN5_reprate.c"
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
#undef TOFL_p2
#undef TOFL_p
#undef TOFL_N
#undef tt_1
#undef tt_0
#undef tmax
#undef tmin
#undef nt
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'Guidelong2' [22]. */
#define mccompcurname  Guidelong2
#define mccompcurtype  Guide
#define mccompcurindex 22
#define pTable mccGuidelong2_pTable
#define reflect mccGuidelong2_reflect
#define w1 mccGuidelong2_w1
#define h1 mccGuidelong2_h1
#define w2 mccGuidelong2_w2
#define h2 mccGuidelong2_h2
#define l mccGuidelong2_l
#define R0 mccGuidelong2_R0
#define Qc mccGuidelong2_Qc
#define alpha mccGuidelong2_alpha
#define m mccGuidelong2_m
#define W mccGuidelong2_W
#line 69 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide.comp"
t_Table pTable;
#line 9085 "./ESS_IN5_reprate.c"
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

/* User declarations for component 'Lmon_beforeballistic' [23]. */
#define mccompcurname  Lmon_beforeballistic
#define mccompcurtype  L_monitor
#define mccompcurindex 23
#define nL mccLmon_beforeballistic_nL
#define L_N mccLmon_beforeballistic_L_N
#define L_p mccLmon_beforeballistic_L_p
#define L_p2 mccLmon_beforeballistic_L_p2
#define filename mccLmon_beforeballistic_filename
#define xmin mccLmon_beforeballistic_xmin
#define xmax mccLmon_beforeballistic_xmax
#define ymin mccLmon_beforeballistic_ymin
#define ymax mccLmon_beforeballistic_ymax
#define xwidth mccLmon_beforeballistic_xwidth
#define yheight mccLmon_beforeballistic_yheight
#define Lmin mccLmon_beforeballistic_Lmin
#define Lmax mccLmon_beforeballistic_Lmax
#define restore_neutron mccLmon_beforeballistic_restore_neutron
#define nowritefile mccLmon_beforeballistic_nowritefile
#line 57 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/L_monitor.comp"
double L_N[nL];
double L_p[nL], L_p2[nL];
#line 9124 "./ESS_IN5_reprate.c"
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

/* User declarations for component 'PSD_beforeballistic' [24]. */
#define mccompcurname  PSD_beforeballistic
#define mccompcurtype  PSD_monitor
#define mccompcurindex 24
#define PSD_N mccPSD_beforeballistic_PSD_N
#define PSD_p mccPSD_beforeballistic_PSD_p
#define PSD_p2 mccPSD_beforeballistic_PSD_p2
#define nx mccPSD_beforeballistic_nx
#define ny mccPSD_beforeballistic_ny
#define filename mccPSD_beforeballistic_filename
#define xmin mccPSD_beforeballistic_xmin
#define xmax mccPSD_beforeballistic_xmax
#define ymin mccPSD_beforeballistic_ymin
#define ymax mccPSD_beforeballistic_ymax
#define xwidth mccPSD_beforeballistic_xwidth
#define yheight mccPSD_beforeballistic_yheight
#define restore_neutron mccPSD_beforeballistic_restore_neutron
#line 62 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
  DArray2d PSD_N;
  DArray2d PSD_p;
  DArray2d PSD_p2;
#line 9165 "./ESS_IN5_reprate.c"
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

/* User declarations for component 'Guidelong2a' [25]. */
#define mccompcurname  Guidelong2a
#define mccompcurtype  Guide
#define mccompcurindex 25
#define pTable mccGuidelong2a_pTable
#define reflect mccGuidelong2a_reflect
#define w1 mccGuidelong2a_w1
#define h1 mccGuidelong2a_h1
#define w2 mccGuidelong2a_w2
#define h2 mccGuidelong2a_h2
#define l mccGuidelong2a_l
#define R0 mccGuidelong2a_R0
#define Qc mccGuidelong2a_Qc
#define alpha mccGuidelong2a_alpha
#define m mccGuidelong2a_m
#define W mccGuidelong2a_W
#line 69 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide.comp"
t_Table pTable;
#line 9201 "./ESS_IN5_reprate.c"
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

/* User declarations for component 'Lmonfast2' [26]. */
#define mccompcurname  Lmonfast2
#define mccompcurtype  L_monitor
#define mccompcurindex 26
#define nL mccLmonfast2_nL
#define L_N mccLmonfast2_L_N
#define L_p mccLmonfast2_L_p
#define L_p2 mccLmonfast2_L_p2
#define filename mccLmonfast2_filename
#define xmin mccLmonfast2_xmin
#define xmax mccLmonfast2_xmax
#define ymin mccLmonfast2_ymin
#define ymax mccLmonfast2_ymax
#define xwidth mccLmonfast2_xwidth
#define yheight mccLmonfast2_yheight
#define Lmin mccLmonfast2_Lmin
#define Lmax mccLmonfast2_Lmax
#define restore_neutron mccLmonfast2_restore_neutron
#define nowritefile mccLmonfast2_nowritefile
#line 57 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/L_monitor.comp"
double L_N[nL];
double L_p[nL], L_p2[nL];
#line 9240 "./ESS_IN5_reprate.c"
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

/* User declarations for component 'Lmonfast2_zoom' [27]. */
#define mccompcurname  Lmonfast2_zoom
#define mccompcurtype  L_monitor
#define mccompcurindex 27
#define nL mccLmonfast2_zoom_nL
#define L_N mccLmonfast2_zoom_L_N
#define L_p mccLmonfast2_zoom_L_p
#define L_p2 mccLmonfast2_zoom_L_p2
#define filename mccLmonfast2_zoom_filename
#define xmin mccLmonfast2_zoom_xmin
#define xmax mccLmonfast2_zoom_xmax
#define ymin mccLmonfast2_zoom_ymin
#define ymax mccLmonfast2_zoom_ymax
#define xwidth mccLmonfast2_zoom_xwidth
#define yheight mccLmonfast2_zoom_yheight
#define Lmin mccLmonfast2_zoom_Lmin
#define Lmax mccLmonfast2_zoom_Lmax
#define restore_neutron mccLmonfast2_zoom_restore_neutron
#define nowritefile mccLmonfast2_zoom_nowritefile
#line 57 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/L_monitor.comp"
double L_N[nL];
double L_p[nL], L_p2[nL];
#line 9282 "./ESS_IN5_reprate.c"
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

/* User declarations for component 'TOFLfast2' [28]. */
#define mccompcurname  TOFLfast2
#define mccompcurtype  TOFLambda_monitor
#define mccompcurindex 28
#define nL mccTOFLfast2_nL
#define nt mccTOFLfast2_nt
#define tmin mccTOFLfast2_tmin
#define tmax mccTOFLfast2_tmax
#define tt_0 mccTOFLfast2_tt_0
#define tt_1 mccTOFLfast2_tt_1
#define TOFL_N mccTOFLfast2_TOFL_N
#define TOFL_p mccTOFLfast2_TOFL_p
#define TOFL_p2 mccTOFLfast2_TOFL_p2
#define filename mccTOFLfast2_filename
#define xmin mccTOFLfast2_xmin
#define xmax mccTOFLfast2_xmax
#define ymin mccTOFLfast2_ymin
#define ymax mccTOFLfast2_ymax
#define xwidth mccTOFLfast2_xwidth
#define yheight mccTOFLfast2_yheight
#define Lmin mccTOFLfast2_Lmin
#define Lmax mccTOFLfast2_Lmax
#define restore_neutron mccTOFLfast2_restore_neutron
#define nowritefile mccTOFLfast2_nowritefile
#line 55 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/TOFLambda_monitor.comp"
double TOFL_N[nt][nL];
double TOFL_p[nt][nL];
double TOFL_p2[nt][nL];
double tt_0, tt_1;
#line 9331 "./ESS_IN5_reprate.c"
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
#undef TOFL_p2
#undef TOFL_p
#undef TOFL_N
#undef tt_1
#undef tt_0
#undef tmax
#undef tmin
#undef nt
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'TOFLfast2zoom' [29]. */
#define mccompcurname  TOFLfast2zoom
#define mccompcurtype  TOFLambda_monitor
#define mccompcurindex 29
#define nL mccTOFLfast2zoom_nL
#define nt mccTOFLfast2zoom_nt
#define tmin mccTOFLfast2zoom_tmin
#define tmax mccTOFLfast2zoom_tmax
#define tt_0 mccTOFLfast2zoom_tt_0
#define tt_1 mccTOFLfast2zoom_tt_1
#define TOFL_N mccTOFLfast2zoom_TOFL_N
#define TOFL_p mccTOFLfast2zoom_TOFL_p
#define TOFL_p2 mccTOFLfast2zoom_TOFL_p2
#define filename mccTOFLfast2zoom_filename
#define xmin mccTOFLfast2zoom_xmin
#define xmax mccTOFLfast2zoom_xmax
#define ymin mccTOFLfast2zoom_ymin
#define ymax mccTOFLfast2zoom_ymax
#define xwidth mccTOFLfast2zoom_xwidth
#define yheight mccTOFLfast2zoom_yheight
#define Lmin mccTOFLfast2zoom_Lmin
#define Lmax mccTOFLfast2zoom_Lmax
#define restore_neutron mccTOFLfast2zoom_restore_neutron
#define nowritefile mccTOFLfast2zoom_nowritefile
#line 55 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/TOFLambda_monitor.comp"
double TOFL_N[nt][nL];
double TOFL_p[nt][nL];
double TOFL_p2[nt][nL];
double tt_0, tt_1;
#line 9385 "./ESS_IN5_reprate.c"
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
#undef TOFL_p2
#undef TOFL_p
#undef TOFL_N
#undef tt_1
#undef tt_0
#undef tmax
#undef tmin
#undef nt
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'PSDfast2' [30]. */
#define mccompcurname  PSDfast2
#define mccompcurtype  PSD_monitor
#define mccompcurindex 30
#define PSD_N mccPSDfast2_PSD_N
#define PSD_p mccPSDfast2_PSD_p
#define PSD_p2 mccPSDfast2_PSD_p2
#define nx mccPSDfast2_nx
#define ny mccPSDfast2_ny
#define filename mccPSDfast2_filename
#define xmin mccPSDfast2_xmin
#define xmax mccPSDfast2_xmax
#define ymin mccPSDfast2_ymin
#define ymax mccPSDfast2_ymax
#define xwidth mccPSDfast2_xwidth
#define yheight mccPSDfast2_yheight
#define restore_neutron mccPSDfast2_restore_neutron
#line 62 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
  DArray2d PSD_N;
  DArray2d PSD_p;
  DArray2d PSD_p2;
#line 9431 "./ESS_IN5_reprate.c"
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

/* User declarations for component 'Fastchop2' [31]. */
#define mccompcurname  Fastchop2
#define mccompcurtype  DiskChopper
#define mccompcurindex 31
#define Tg mccFastchop2_Tg
#define To mccFastchop2_To
#define delta_y mccFastchop2_delta_y
#define height mccFastchop2_height
#define omega mccFastchop2_omega
#define theta_0 mccFastchop2_theta_0
#define radius mccFastchop2_radius
#define yheight mccFastchop2_yheight
#define nu mccFastchop2_nu
#define nslit mccFastchop2_nslit
#define jitter mccFastchop2_jitter
#define delay mccFastchop2_delay
#define isfirst mccFastchop2_isfirst
#define n_pulse mccFastchop2_n_pulse
#define abs_out mccFastchop2_abs_out
#define phase mccFastchop2_phase
#define xwidth mccFastchop2_xwidth
#define verbose mccFastchop2_verbose
#line 66 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/DiskChopper.comp"
double Tg,To,delta_y,height,omega;
#line 9473 "./ESS_IN5_reprate.c"
#undef verbose
#undef xwidth
#undef phase
#undef abs_out
#undef n_pulse
#undef isfirst
#undef delay
#undef jitter
#undef nslit
#undef nu
#undef yheight
#undef radius
#undef theta_0
#undef omega
#undef height
#undef delta_y
#undef To
#undef Tg
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'Fastchop2counter' [32]. */
#define mccompcurname  Fastchop2counter
#define mccompcurtype  DiskChopper
#define mccompcurindex 32
#define Tg mccFastchop2counter_Tg
#define To mccFastchop2counter_To
#define delta_y mccFastchop2counter_delta_y
#define height mccFastchop2counter_height
#define omega mccFastchop2counter_omega
#define theta_0 mccFastchop2counter_theta_0
#define radius mccFastchop2counter_radius
#define yheight mccFastchop2counter_yheight
#define nu mccFastchop2counter_nu
#define nslit mccFastchop2counter_nslit
#define jitter mccFastchop2counter_jitter
#define delay mccFastchop2counter_delay
#define isfirst mccFastchop2counter_isfirst
#define n_pulse mccFastchop2counter_n_pulse
#define abs_out mccFastchop2counter_abs_out
#define phase mccFastchop2counter_phase
#define xwidth mccFastchop2counter_xwidth
#define verbose mccFastchop2counter_verbose
#line 66 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/DiskChopper.comp"
double Tg,To,delta_y,height,omega;
#line 9520 "./ESS_IN5_reprate.c"
#undef verbose
#undef xwidth
#undef phase
#undef abs_out
#undef n_pulse
#undef isfirst
#undef delay
#undef jitter
#undef nslit
#undef nu
#undef yheight
#undef radius
#undef theta_0
#undef omega
#undef height
#undef delta_y
#undef To
#undef Tg
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'FOchop3' [33]. */
#define mccompcurname  FOchop3
#define mccompcurtype  DiskChopper
#define mccompcurindex 33
#define Tg mccFOchop3_Tg
#define To mccFOchop3_To
#define delta_y mccFOchop3_delta_y
#define height mccFOchop3_height
#define omega mccFOchop3_omega
#define theta_0 mccFOchop3_theta_0
#define radius mccFOchop3_radius
#define yheight mccFOchop3_yheight
#define nu mccFOchop3_nu
#define nslit mccFOchop3_nslit
#define jitter mccFOchop3_jitter
#define delay mccFOchop3_delay
#define isfirst mccFOchop3_isfirst
#define n_pulse mccFOchop3_n_pulse
#define abs_out mccFOchop3_abs_out
#define phase mccFOchop3_phase
#define xwidth mccFOchop3_xwidth
#define verbose mccFOchop3_verbose
#line 66 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/DiskChopper.comp"
double Tg,To,delta_y,height,omega;
#line 9567 "./ESS_IN5_reprate.c"
#undef verbose
#undef xwidth
#undef phase
#undef abs_out
#undef n_pulse
#undef isfirst
#undef delay
#undef jitter
#undef nslit
#undef nu
#undef yheight
#undef radius
#undef theta_0
#undef omega
#undef height
#undef delta_y
#undef To
#undef Tg
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'TOFfast2_zoom' [34]. */
#define mccompcurname  TOFfast2_zoom
#define mccompcurtype  TOF_monitor
#define mccompcurindex 34
#define nt mccTOFfast2_zoom_nt
#define TOF_N mccTOFfast2_zoom_TOF_N
#define TOF_p mccTOFfast2_zoom_TOF_p
#define TOF_p2 mccTOFfast2_zoom_TOF_p2
#define t_min mccTOFfast2_zoom_t_min
#define t_max mccTOFfast2_zoom_t_max
#define delta_t mccTOFfast2_zoom_delta_t
#define filename mccTOFfast2_zoom_filename
#define xmin mccTOFfast2_zoom_xmin
#define xmax mccTOFfast2_zoom_xmax
#define ymin mccTOFfast2_zoom_ymin
#define ymax mccTOFfast2_zoom_ymax
#define xwidth mccTOFfast2_zoom_xwidth
#define yheight mccTOFfast2_zoom_yheight
#define tmin mccTOFfast2_zoom_tmin
#define tmax mccTOFfast2_zoom_tmax
#define dt mccTOFfast2_zoom_dt
#define restore_neutron mccTOFfast2_zoom_restore_neutron
#define nowritefile mccTOFfast2_zoom_nowritefile
#line 54 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/TOF_monitor.comp"
double TOF_N[nt];
double TOF_p[nt];
double TOF_p2[nt];
double t_min, t_max, delta_t;
#line 9618 "./ESS_IN5_reprate.c"
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

/* User declarations for component 'Lmon_afterfast2' [35]. */
#define mccompcurname  Lmon_afterfast2
#define mccompcurtype  L_monitor
#define mccompcurindex 35
#define nL mccLmon_afterfast2_nL
#define L_N mccLmon_afterfast2_L_N
#define L_p mccLmon_afterfast2_L_p
#define L_p2 mccLmon_afterfast2_L_p2
#define filename mccLmon_afterfast2_filename
#define xmin mccLmon_afterfast2_xmin
#define xmax mccLmon_afterfast2_xmax
#define ymin mccLmon_afterfast2_ymin
#define ymax mccLmon_afterfast2_ymax
#define xwidth mccLmon_afterfast2_xwidth
#define yheight mccLmon_afterfast2_yheight
#define Lmin mccLmon_afterfast2_Lmin
#define Lmax mccLmon_afterfast2_Lmax
#define restore_neutron mccLmon_afterfast2_restore_neutron
#define nowritefile mccLmon_afterfast2_nowritefile
#line 57 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/L_monitor.comp"
double L_N[nL];
double L_p[nL], L_p2[nL];
#line 9664 "./ESS_IN5_reprate.c"
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

/* User declarations for component 'TOFL_afterfast2' [36]. */
#define mccompcurname  TOFL_afterfast2
#define mccompcurtype  TOFLambda_monitor
#define mccompcurindex 36
#define nL mccTOFL_afterfast2_nL
#define nt mccTOFL_afterfast2_nt
#define tmin mccTOFL_afterfast2_tmin
#define tmax mccTOFL_afterfast2_tmax
#define tt_0 mccTOFL_afterfast2_tt_0
#define tt_1 mccTOFL_afterfast2_tt_1
#define TOFL_N mccTOFL_afterfast2_TOFL_N
#define TOFL_p mccTOFL_afterfast2_TOFL_p
#define TOFL_p2 mccTOFL_afterfast2_TOFL_p2
#define filename mccTOFL_afterfast2_filename
#define xmin mccTOFL_afterfast2_xmin
#define xmax mccTOFL_afterfast2_xmax
#define ymin mccTOFL_afterfast2_ymin
#define ymax mccTOFL_afterfast2_ymax
#define xwidth mccTOFL_afterfast2_xwidth
#define yheight mccTOFL_afterfast2_yheight
#define Lmin mccTOFL_afterfast2_Lmin
#define Lmax mccTOFL_afterfast2_Lmax
#define restore_neutron mccTOFL_afterfast2_restore_neutron
#define nowritefile mccTOFL_afterfast2_nowritefile
#line 55 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/TOFLambda_monitor.comp"
double TOFL_N[nt][nL];
double TOFL_p[nt][nL];
double TOFL_p2[nt][nL];
double tt_0, tt_1;
#line 9713 "./ESS_IN5_reprate.c"
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
#undef TOFL_p2
#undef TOFL_p
#undef TOFL_N
#undef tt_1
#undef tt_0
#undef tmax
#undef tmin
#undef nt
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'TOFL_afterfast2_zoom' [37]. */
#define mccompcurname  TOFL_afterfast2_zoom
#define mccompcurtype  TOFLambda_monitor
#define mccompcurindex 37
#define nL mccTOFL_afterfast2_zoom_nL
#define nt mccTOFL_afterfast2_zoom_nt
#define tmin mccTOFL_afterfast2_zoom_tmin
#define tmax mccTOFL_afterfast2_zoom_tmax
#define tt_0 mccTOFL_afterfast2_zoom_tt_0
#define tt_1 mccTOFL_afterfast2_zoom_tt_1
#define TOFL_N mccTOFL_afterfast2_zoom_TOFL_N
#define TOFL_p mccTOFL_afterfast2_zoom_TOFL_p
#define TOFL_p2 mccTOFL_afterfast2_zoom_TOFL_p2
#define filename mccTOFL_afterfast2_zoom_filename
#define xmin mccTOFL_afterfast2_zoom_xmin
#define xmax mccTOFL_afterfast2_zoom_xmax
#define ymin mccTOFL_afterfast2_zoom_ymin
#define ymax mccTOFL_afterfast2_zoom_ymax
#define xwidth mccTOFL_afterfast2_zoom_xwidth
#define yheight mccTOFL_afterfast2_zoom_yheight
#define Lmin mccTOFL_afterfast2_zoom_Lmin
#define Lmax mccTOFL_afterfast2_zoom_Lmax
#define restore_neutron mccTOFL_afterfast2_zoom_restore_neutron
#define nowritefile mccTOFL_afterfast2_zoom_nowritefile
#line 55 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/TOFLambda_monitor.comp"
double TOFL_N[nt][nL];
double TOFL_p[nt][nL];
double TOFL_p2[nt][nL];
double tt_0, tt_1;
#line 9767 "./ESS_IN5_reprate.c"
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
#undef TOFL_p2
#undef TOFL_p
#undef TOFL_N
#undef tt_1
#undef tt_0
#undef tmax
#undef tmin
#undef nt
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'PSD_afterfast2' [38]. */
#define mccompcurname  PSD_afterfast2
#define mccompcurtype  PSD_monitor
#define mccompcurindex 38
#define PSD_N mccPSD_afterfast2_PSD_N
#define PSD_p mccPSD_afterfast2_PSD_p
#define PSD_p2 mccPSD_afterfast2_PSD_p2
#define nx mccPSD_afterfast2_nx
#define ny mccPSD_afterfast2_ny
#define filename mccPSD_afterfast2_filename
#define xmin mccPSD_afterfast2_xmin
#define xmax mccPSD_afterfast2_xmax
#define ymin mccPSD_afterfast2_ymin
#define ymax mccPSD_afterfast2_ymax
#define xwidth mccPSD_afterfast2_xwidth
#define yheight mccPSD_afterfast2_yheight
#define restore_neutron mccPSD_afterfast2_restore_neutron
#line 62 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
  DArray2d PSD_N;
  DArray2d PSD_p;
  DArray2d PSD_p2;
#line 9813 "./ESS_IN5_reprate.c"
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

/* User declarations for component 'Guidesample' [39]. */
#define mccompcurname  Guidesample
#define mccompcurtype  Guide
#define mccompcurindex 39
#define pTable mccGuidesample_pTable
#define reflect mccGuidesample_reflect
#define w1 mccGuidesample_w1
#define h1 mccGuidesample_h1
#define w2 mccGuidesample_w2
#define h2 mccGuidesample_h2
#define l mccGuidesample_l
#define R0 mccGuidesample_R0
#define Qc mccGuidesample_Qc
#define alpha mccGuidesample_alpha
#define m mccGuidesample_m
#define W mccGuidesample_W
#line 69 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide.comp"
t_Table pTable;
#line 9849 "./ESS_IN5_reprate.c"
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

/* User declarations for component 'Lmon_guideend' [40]. */
#define mccompcurname  Lmon_guideend
#define mccompcurtype  L_monitor
#define mccompcurindex 40
#define nL mccLmon_guideend_nL
#define L_N mccLmon_guideend_L_N
#define L_p mccLmon_guideend_L_p
#define L_p2 mccLmon_guideend_L_p2
#define filename mccLmon_guideend_filename
#define xmin mccLmon_guideend_xmin
#define xmax mccLmon_guideend_xmax
#define ymin mccLmon_guideend_ymin
#define ymax mccLmon_guideend_ymax
#define xwidth mccLmon_guideend_xwidth
#define yheight mccLmon_guideend_yheight
#define Lmin mccLmon_guideend_Lmin
#define Lmax mccLmon_guideend_Lmax
#define restore_neutron mccLmon_guideend_restore_neutron
#define nowritefile mccLmon_guideend_nowritefile
#line 57 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/L_monitor.comp"
double L_N[nL];
double L_p[nL], L_p2[nL];
#line 9888 "./ESS_IN5_reprate.c"
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

/* User declarations for component 'PSDsample' [41]. */
#define mccompcurname  PSDsample
#define mccompcurtype  PSD_monitor
#define mccompcurindex 41
#define PSD_N mccPSDsample_PSD_N
#define PSD_p mccPSDsample_PSD_p
#define PSD_p2 mccPSDsample_PSD_p2
#define nx mccPSDsample_nx
#define ny mccPSDsample_ny
#define filename mccPSDsample_filename
#define xmin mccPSDsample_xmin
#define xmax mccPSDsample_xmax
#define ymin mccPSDsample_ymin
#define ymax mccPSDsample_ymax
#define xwidth mccPSDsample_xwidth
#define yheight mccPSDsample_yheight
#define restore_neutron mccPSDsample_restore_neutron
#line 62 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
  DArray2d PSD_N;
  DArray2d PSD_p;
  DArray2d PSD_p2;
#line 9929 "./ESS_IN5_reprate.c"
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

/* User declarations for component 'TOFsample_zoom' [42]. */
#define mccompcurname  TOFsample_zoom
#define mccompcurtype  TOF_monitor
#define mccompcurindex 42
#define nt mccTOFsample_zoom_nt
#define TOF_N mccTOFsample_zoom_TOF_N
#define TOF_p mccTOFsample_zoom_TOF_p
#define TOF_p2 mccTOFsample_zoom_TOF_p2
#define t_min mccTOFsample_zoom_t_min
#define t_max mccTOFsample_zoom_t_max
#define delta_t mccTOFsample_zoom_delta_t
#define filename mccTOFsample_zoom_filename
#define xmin mccTOFsample_zoom_xmin
#define xmax mccTOFsample_zoom_xmax
#define ymin mccTOFsample_zoom_ymin
#define ymax mccTOFsample_zoom_ymax
#define xwidth mccTOFsample_zoom_xwidth
#define yheight mccTOFsample_zoom_yheight
#define tmin mccTOFsample_zoom_tmin
#define tmax mccTOFsample_zoom_tmax
#define dt mccTOFsample_zoom_dt
#define restore_neutron mccTOFsample_zoom_restore_neutron
#define nowritefile mccTOFsample_zoom_nowritefile
#line 54 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/TOF_monitor.comp"
double TOF_N[nt];
double TOF_p[nt];
double TOF_p2[nt];
double t_min, t_max, delta_t;
#line 9975 "./ESS_IN5_reprate.c"
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

/* User declarations for component 'Esample' [43]. */
#define mccompcurname  Esample
#define mccompcurtype  E_monitor
#define mccompcurindex 43
#define nE mccEsample_nE
#define E_N mccEsample_E_N
#define E_p mccEsample_E_p
#define E_p2 mccEsample_E_p2
#define S_p mccEsample_S_p
#define S_pE mccEsample_S_pE
#define S_pE2 mccEsample_S_pE2
#define filename mccEsample_filename
#define xmin mccEsample_xmin
#define xmax mccEsample_xmax
#define ymin mccEsample_ymin
#define ymax mccEsample_ymax
#define xwidth mccEsample_xwidth
#define yheight mccEsample_yheight
#define Emin mccEsample_Emin
#define Emax mccEsample_Emax
#define restore_neutron mccEsample_restore_neutron
#define nowritefile mccEsample_nowritefile
#line 60 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/E_monitor.comp"
double E_N[nE];
double E_p[nE], E_p2[nE];
double S_p, S_pE, S_pE2;
#line 10025 "./ESS_IN5_reprate.c"
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

/* User declarations for component 'Lmon_sample_zoom' [44]. */
#define mccompcurname  Lmon_sample_zoom
#define mccompcurtype  L_monitor
#define mccompcurindex 44
#define nL mccLmon_sample_zoom_nL
#define L_N mccLmon_sample_zoom_L_N
#define L_p mccLmon_sample_zoom_L_p
#define L_p2 mccLmon_sample_zoom_L_p2
#define filename mccLmon_sample_zoom_filename
#define xmin mccLmon_sample_zoom_xmin
#define xmax mccLmon_sample_zoom_xmax
#define ymin mccLmon_sample_zoom_ymin
#define ymax mccLmon_sample_zoom_ymax
#define xwidth mccLmon_sample_zoom_xwidth
#define yheight mccLmon_sample_zoom_yheight
#define Lmin mccLmon_sample_zoom_Lmin
#define Lmax mccLmon_sample_zoom_Lmax
#define restore_neutron mccLmon_sample_zoom_restore_neutron
#define nowritefile mccLmon_sample_zoom_nowritefile
#line 57 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/L_monitor.comp"
double L_N[nL];
double L_p[nL], L_p2[nL];
#line 10070 "./ESS_IN5_reprate.c"
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

/* User declarations for component 'sample' [45]. */
#define mccompcurname  sample
#define mccompcurtype  Tunneling_sample
#define mccompcurindex 45
#define ftun mccsample_ftun
#define fQE mccsample_fQE
#define VarsV mccsample_VarsV
#define thickness mccsample_thickness
#define radius mccsample_radius
#define focus_r mccsample_focus_r
#define p_interact mccsample_p_interact
#define f_QE mccsample_f_QE
#define f_tun mccsample_f_tun
#define gamma mccsample_gamma
#define E_tun mccsample_E_tun
#define target_x mccsample_target_x
#define target_y mccsample_target_y
#define target_z mccsample_target_z
#define focus_xw mccsample_focus_xw
#define focus_yh mccsample_focus_yh
#define focus_aw mccsample_focus_aw
#define focus_ah mccsample_focus_ah
#define xwidth mccsample_xwidth
#define yheight mccsample_yheight
#define zdepth mccsample_zdepth
#define sigma_abs mccsample_sigma_abs
#define sigma_inc mccsample_sigma_inc
#define Vc mccsample_Vc
#define target_index mccsample_target_index
#line 109 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../samples/Tunneling_sample.comp"
  struct StructVarsV VarsV;
  double ftun, fQE;
#line 10122 "./ESS_IN5_reprate.c"
#undef target_index
#undef Vc
#undef sigma_inc
#undef sigma_abs
#undef zdepth
#undef yheight
#undef xwidth
#undef focus_ah
#undef focus_aw
#undef focus_yh
#undef focus_xw
#undef target_z
#undef target_y
#undef target_x
#undef E_tun
#undef gamma
#undef f_tun
#undef f_QE
#undef p_interact
#undef focus_r
#undef radius
#undef thickness
#undef VarsV
#undef fQE
#undef ftun
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'detectorarm' [46]. */
#define mccompcurname  detectorarm
#define mccompcurtype  Arm
#define mccompcurindex 46
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'TOFdetector' [47]. */
#define mccompcurname  TOFdetector
#define mccompcurtype  TOF_monitor
#define mccompcurindex 47
#define nt mccTOFdetector_nt
#define TOF_N mccTOFdetector_TOF_N
#define TOF_p mccTOFdetector_TOF_p
#define TOF_p2 mccTOFdetector_TOF_p2
#define t_min mccTOFdetector_t_min
#define t_max mccTOFdetector_t_max
#define delta_t mccTOFdetector_delta_t
#define filename mccTOFdetector_filename
#define xmin mccTOFdetector_xmin
#define xmax mccTOFdetector_xmax
#define ymin mccTOFdetector_ymin
#define ymax mccTOFdetector_ymax
#define xwidth mccTOFdetector_xwidth
#define yheight mccTOFdetector_yheight
#define tmin mccTOFdetector_tmin
#define tmax mccTOFdetector_tmax
#define dt mccTOFdetector_dt
#define restore_neutron mccTOFdetector_restore_neutron
#define nowritefile mccTOFdetector_nowritefile
#line 54 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/TOF_monitor.comp"
double TOF_N[nt];
double TOF_p[nt];
double TOF_p2[nt];
double t_min, t_max, delta_t;
#line 10188 "./ESS_IN5_reprate.c"
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

/* User declarations for component 'TOFdetector_zoom' [48]. */
#define mccompcurname  TOFdetector_zoom
#define mccompcurtype  TOF_monitor
#define mccompcurindex 48
#define nt mccTOFdetector_zoom_nt
#define TOF_N mccTOFdetector_zoom_TOF_N
#define TOF_p mccTOFdetector_zoom_TOF_p
#define TOF_p2 mccTOFdetector_zoom_TOF_p2
#define t_min mccTOFdetector_zoom_t_min
#define t_max mccTOFdetector_zoom_t_max
#define delta_t mccTOFdetector_zoom_delta_t
#define filename mccTOFdetector_zoom_filename
#define xmin mccTOFdetector_zoom_xmin
#define xmax mccTOFdetector_zoom_xmax
#define ymin mccTOFdetector_zoom_ymin
#define ymax mccTOFdetector_zoom_ymax
#define xwidth mccTOFdetector_zoom_xwidth
#define yheight mccTOFdetector_zoom_yheight
#define tmin mccTOFdetector_zoom_tmin
#define tmax mccTOFdetector_zoom_tmax
#define dt mccTOFdetector_zoom_dt
#define restore_neutron mccTOFdetector_zoom_restore_neutron
#define nowritefile mccTOFdetector_zoom_nowritefile
#line 54 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/TOF_monitor.comp"
double TOF_N[nt];
double TOF_p[nt];
double TOF_p2[nt];
double t_min, t_max, delta_t;
#line 10240 "./ESS_IN5_reprate.c"
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

/* User declarations for component 'Edetector' [49]. */
#define mccompcurname  Edetector
#define mccompcurtype  E_monitor
#define mccompcurindex 49
#define nE mccEdetector_nE
#define E_N mccEdetector_E_N
#define E_p mccEdetector_E_p
#define E_p2 mccEdetector_E_p2
#define S_p mccEdetector_S_p
#define S_pE mccEdetector_S_pE
#define S_pE2 mccEdetector_S_pE2
#define filename mccEdetector_filename
#define xmin mccEdetector_xmin
#define xmax mccEdetector_xmax
#define ymin mccEdetector_ymin
#define ymax mccEdetector_ymax
#define xwidth mccEdetector_xwidth
#define yheight mccEdetector_yheight
#define Emin mccEdetector_Emin
#define Emax mccEdetector_Emax
#define restore_neutron mccEdetector_restore_neutron
#define nowritefile mccEdetector_nowritefile
#line 60 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/E_monitor.comp"
double E_N[nE];
double E_p[nE], E_p2[nE];
double S_p, S_pE, S_pE2;
#line 10290 "./ESS_IN5_reprate.c"
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

/* User declarations for component 'TOF2Edetector' [50]. */
#define mccompcurname  TOF2Edetector
#define mccompcurtype  TOF2E_monitor
#define mccompcurindex 50
#define nE mccTOF2Edetector_nE
#define E_N mccTOF2Edetector_E_N
#define E_p mccTOF2Edetector_E_p
#define E_p2 mccTOF2Edetector_E_p2
#define S_p mccTOF2Edetector_S_p
#define S_pE mccTOF2Edetector_S_pE
#define S_pE2 mccTOF2Edetector_S_pE2
#define filename mccTOF2Edetector_filename
#define xmin mccTOF2Edetector_xmin
#define xmax mccTOF2Edetector_xmax
#define ymin mccTOF2Edetector_ymin
#define ymax mccTOF2Edetector_ymax
#define xwidth mccTOF2Edetector_xwidth
#define yheight mccTOF2Edetector_yheight
#define Emin mccTOF2Edetector_Emin
#define Emax mccTOF2Edetector_Emax
#define T_zero mccTOF2Edetector_T_zero
#define L_flight mccTOF2Edetector_L_flight
#define restore_neutron mccTOF2Edetector_restore_neutron
#define nowritefile mccTOF2Edetector_nowritefile
#line 61 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/TOF2E_monitor.comp"
double E_N[nE];
double E_p[nE], E_p2[nE];
double S_p, S_pE, S_pE2;
#line 10341 "./ESS_IN5_reprate.c"
#undef nowritefile
#undef restore_neutron
#undef L_flight
#undef T_zero
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

Coords mcposasource, mcposrsource;
Rotation mcrotasource, mcrotrsource;
Coords mcposaOrigin, mcposrOrigin;
Rotation mcrotaOrigin, mcrotrOrigin;
Coords mcposaTOFmoderator_zoom, mcposrTOFmoderator_zoom;
Rotation mcrotaTOFmoderator_zoom, mcrotrTOFmoderator_zoom;
Coords mcposaTOFmoderator, mcposrTOFmoderator;
Rotation mcrotaTOFmoderator, mcrotrTOFmoderator;
Coords mcposaLmon_guistart, mcposrLmon_guistart;
Rotation mcrotaLmon_guistart, mcrotrLmon_guistart;
Coords mcposaLmon_normalize, mcposrLmon_normalize;
Rotation mcrotaLmon_normalize, mcrotrLmon_normalize;
Coords mcposaGuide1, mcposrGuide1;
Rotation mcrotaGuide1, mcrotrGuide1;
Coords mcposaLmonslow1, mcposrLmonslow1;
Rotation mcrotaLmonslow1, mcrotrLmonslow1;
Coords mcposaPSDslow1, mcposrPSDslow1;
Rotation mcrotaPSDslow1, mcrotrPSDslow1;
Coords mcposaFOchop1, mcposrFOchop1;
Rotation mcrotaFOchop1, mcrotrFOchop1;
Coords mcposaTOFLmon1, mcposrTOFLmon1;
Rotation mcrotaTOFLmon1, mcrotrTOFLmon1;
Coords mcposaLmon_afterslow1, mcposrLmon_afterslow1;
Rotation mcrotaLmon_afterslow1, mcrotrLmon_afterslow1;
Coords mcposaPSD_afterslow1, mcposrPSD_afterslow1;
Rotation mcrotaPSD_afterslow1, mcrotrPSD_afterslow1;
Coords mcposaGuidelong1, mcposrGuidelong1;
Rotation mcrotaGuidelong1, mcrotrGuidelong1;
Coords mcposaGuidelong1b, mcposrGuidelong1b;
Rotation mcrotaGuidelong1b, mcrotrGuidelong1b;
Coords mcposaLmon_slow2, mcposrLmon_slow2;
Rotation mcrotaLmon_slow2, mcrotrLmon_slow2;
Coords mcposaFOchop2, mcposrFOchop2;
Rotation mcrotaFOchop2, mcrotrFOchop2;
Coords mcposaFastchop1, mcposrFastchop1;
Rotation mcrotaFastchop1, mcrotrFastchop1;
Coords mcposaPSD_afterslow2, mcposrPSD_afterslow2;
Rotation mcrotaPSD_afterslow2, mcrotrPSD_afterslow2;
Coords mcposaLmon_afterslow2, mcposrLmon_afterslow2;
Rotation mcrotaLmon_afterslow2, mcrotrLmon_afterslow2;
Coords mcposaTOFL_afterslow2, mcposrTOFL_afterslow2;
Rotation mcrotaTOFL_afterslow2, mcrotrTOFL_afterslow2;
Coords mcposaGuidelong2, mcposrGuidelong2;
Rotation mcrotaGuidelong2, mcrotrGuidelong2;
Coords mcposaLmon_beforeballistic, mcposrLmon_beforeballistic;
Rotation mcrotaLmon_beforeballistic, mcrotrLmon_beforeballistic;
Coords mcposaPSD_beforeballistic, mcposrPSD_beforeballistic;
Rotation mcrotaPSD_beforeballistic, mcrotrPSD_beforeballistic;
Coords mcposaGuidelong2a, mcposrGuidelong2a;
Rotation mcrotaGuidelong2a, mcrotrGuidelong2a;
Coords mcposaLmonfast2, mcposrLmonfast2;
Rotation mcrotaLmonfast2, mcrotrLmonfast2;
Coords mcposaLmonfast2_zoom, mcposrLmonfast2_zoom;
Rotation mcrotaLmonfast2_zoom, mcrotrLmonfast2_zoom;
Coords mcposaTOFLfast2, mcposrTOFLfast2;
Rotation mcrotaTOFLfast2, mcrotrTOFLfast2;
Coords mcposaTOFLfast2zoom, mcposrTOFLfast2zoom;
Rotation mcrotaTOFLfast2zoom, mcrotrTOFLfast2zoom;
Coords mcposaPSDfast2, mcposrPSDfast2;
Rotation mcrotaPSDfast2, mcrotrPSDfast2;
Coords mcposaFastchop2, mcposrFastchop2;
Rotation mcrotaFastchop2, mcrotrFastchop2;
Coords mcposaFastchop2counter, mcposrFastchop2counter;
Rotation mcrotaFastchop2counter, mcrotrFastchop2counter;
Coords mcposaFOchop3, mcposrFOchop3;
Rotation mcrotaFOchop3, mcrotrFOchop3;
Coords mcposaTOFfast2_zoom, mcposrTOFfast2_zoom;
Rotation mcrotaTOFfast2_zoom, mcrotrTOFfast2_zoom;
Coords mcposaLmon_afterfast2, mcposrLmon_afterfast2;
Rotation mcrotaLmon_afterfast2, mcrotrLmon_afterfast2;
Coords mcposaTOFL_afterfast2, mcposrTOFL_afterfast2;
Rotation mcrotaTOFL_afterfast2, mcrotrTOFL_afterfast2;
Coords mcposaTOFL_afterfast2_zoom, mcposrTOFL_afterfast2_zoom;
Rotation mcrotaTOFL_afterfast2_zoom, mcrotrTOFL_afterfast2_zoom;
Coords mcposaPSD_afterfast2, mcposrPSD_afterfast2;
Rotation mcrotaPSD_afterfast2, mcrotrPSD_afterfast2;
Coords mcposaGuidesample, mcposrGuidesample;
Rotation mcrotaGuidesample, mcrotrGuidesample;
Coords mcposaLmon_guideend, mcposrLmon_guideend;
Rotation mcrotaLmon_guideend, mcrotrLmon_guideend;
Coords mcposaPSDsample, mcposrPSDsample;
Rotation mcrotaPSDsample, mcrotrPSDsample;
Coords mcposaTOFsample_zoom, mcposrTOFsample_zoom;
Rotation mcrotaTOFsample_zoom, mcrotrTOFsample_zoom;
Coords mcposaEsample, mcposrEsample;
Rotation mcrotaEsample, mcrotrEsample;
Coords mcposaLmon_sample_zoom, mcposrLmon_sample_zoom;
Rotation mcrotaLmon_sample_zoom, mcrotrLmon_sample_zoom;
Coords mcposasample, mcposrsample;
Rotation mcrotasample, mcrotrsample;
Coords mcposadetectorarm, mcposrdetectorarm;
Rotation mcrotadetectorarm, mcrotrdetectorarm;
Coords mcposaTOFdetector, mcposrTOFdetector;
Rotation mcrotaTOFdetector, mcrotrTOFdetector;
Coords mcposaTOFdetector_zoom, mcposrTOFdetector_zoom;
Rotation mcrotaTOFdetector_zoom, mcrotrTOFdetector_zoom;
Coords mcposaEdetector, mcposrEdetector;
Rotation mcrotaEdetector, mcrotrEdetector;
Coords mcposaTOF2Edetector, mcposrTOF2Edetector;
Rotation mcrotaTOF2Edetector, mcrotrTOF2Edetector;

MCNUM mcnx, mcny, mcnz, mcnvx, mcnvy, mcnvz, mcnt, mcnsx, mcnsy, mcnsz, mcnp;

/* end declare */

void mcinit(void) {
#define mccompcurname  ESS_IN5_reprate
#define mccompcurtype  INSTRUMENT
#define mccompcurindex 0
#define mcposaESS_IN5_reprate coords_set(0,0,0)
#define Lmin mcipLmin
#define Lmax mcipLmax
#define lambda0 mciplambda0
#define Pulse_width mcipPulse_width
#define Num_pulses mcipNum_pulses
#define GUI_start mcipGUI_start
#define FO1_DIST mcipFO1_DIST
#define L_ballistic_begin mcipL_ballistic_begin
#define L_ballistic_end mcipL_ballistic_end
#define Length mcipLength
#define SAMPLE_DIST mcipSAMPLE_DIST
#define DETECTOR_DIST mcipDETECTOR_DIST
#define GUI_h mcipGUI_h
#define GUI_w mcipGUI_w
#define GUI_GAP mcipGUI_GAP
#define H1 mcipH1
#define W1 mcipW1
#define H2 mcipH2
#define W2 mcipW2
#define H3 mcipH3
#define W3 mcipW3
#define H4 mcipH4
#define W4 mcipW4
#define H_chop mcipH_chop
#define W_chop mcipW_chop
#define H_end mcipH_end
#define W_end mcipW_end
#define ALPHA mcipALPHA
#define M mcipM
#define F_slow1 mcipF_slow1
#define F_slow2 mcipF_slow2
#define F_fast1 mcipF_fast1
#define F_fast2 mcipF_fast2
#define N_fast mcipN_fast
#define SLOW1_THETA mcipSLOW1_THETA
#define FO3 mcipFO3
#define THETA_fast1 mcipTHETA_fast1
#define FAST_THETA mcipFAST_THETA
#define Gamma mcipGamma
#define Etun mcipEtun
#define V_HOLE mcipV_HOLE
#define FRAC_QUASIEL mcipFRAC_QUASIEL
#define FRAC_TUNNEL mcipFRAC_TUNNEL
#define TT mcipTT
#define RES_DE mcipRES_DE
#define port mcipport
#define cold mcipcold
#line 95 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
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
#line 10539 "./ESS_IN5_reprate.c"
#undef cold
#undef port
#undef RES_DE
#undef TT
#undef FRAC_TUNNEL
#undef FRAC_QUASIEL
#undef V_HOLE
#undef Etun
#undef Gamma
#undef FAST_THETA
#undef THETA_fast1
#undef FO3
#undef SLOW1_THETA
#undef N_fast
#undef F_fast2
#undef F_fast1
#undef F_slow2
#undef F_slow1
#undef M
#undef ALPHA
#undef W_end
#undef H_end
#undef W_chop
#undef H_chop
#undef W4
#undef H4
#undef W3
#undef H3
#undef W2
#undef H2
#undef W1
#undef H1
#undef GUI_GAP
#undef GUI_w
#undef GUI_h
#undef DETECTOR_DIST
#undef SAMPLE_DIST
#undef Length
#undef L_ballistic_end
#undef L_ballistic_begin
#undef FO1_DIST
#undef GUI_start
#undef Num_pulses
#undef Pulse_width
#undef lambda0
#undef Lmax
#undef Lmin
#undef mcposaESS_IN5_reprate
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
    /* Component source. */
  /* Setting parameters for component source. */
  SIG_MESSAGE("source (Init:SetPar)");
#line 125 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccsource_width_c = 0;
#line 125 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccsource_yheight = 0.12;
#line 116 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccsource_Lmin = mcipLmin;
#line 116 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccsource_Lmax = mcipLmax;
#line 116 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccsource_dist = mcipGUI_start;
#line 116 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccsource_focus_xw = mcipGUI_w;
#line 117 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccsource_focus_yh = mcipGUI_h;
#line 117 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccsource_nu = FREQ;
#line 117 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccsource_T = 50;
#line 117 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccsource_tau = 287e-6;
#line 117 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccsource_tau1 = 0;
#line 118 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccsource_tau2 = 20e-6;
#line 118 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccsource_d = mcipPulse_width;
#line 118 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccsource_n = 20;
#line 120 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccsource_cold_frac = mcipcold;
#line 118 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccsource_n2 = 5;
#line 118 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccsource_chi2 = 0.9;
#line 119 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccsource_I0 = 6.9e11;
#line 119 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccsource_I2 = 27.6e10;
#line 127 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccsource_target_index = 0;
#line 116 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccsource_cyl_radius = 0.085;
#line 119 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccsource_branch1 = 1;
#line 119 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccsource_branch2 = 0.5;
#line 120 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccsource_branch_tail = 0.1;
#line 120 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccsource_n_pulses = mcipNum_pulses;
#line 129 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccsource_width_t = 0.12;
#line 129 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccsource_T_t = 325;
#line 129 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccsource_tau_t = 80e-6;
#line 129 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccsource_tau1_t = 400e-6;
#line 130 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccsource_tau2_t = 12e-6;
#line 130 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccsource_chi2_t = 2.5;
#line 130 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccsource_I0_t = 13.5e11;
#line 130 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccsource_I2_t = 27.6e10;
#line 130 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccsource_branch1_t = 0.5;
#line 131 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccsource_branch2_t = 0.5;
#line 131 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccsource_src_2012 = 1;
#line 131 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccsource_tfocus_dist = 0.1;
#line 131 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccsource_tfocus_time = 0.0;
#line 131 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccsource_tfocus_width = 0.0;
#line 120 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccsource_beamport_angle = mcipport;
#line 10684 "./ESS_IN5_reprate.c"

  SIG_MESSAGE("source (Init:Place/Rotate)");
  rot_set_rotation(mcrotasource,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 10691 "./ESS_IN5_reprate.c"
  rot_copy(mcrotrsource, mcrotasource);
  mcposasource = coords_set(
#line 121 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    0,
#line 121 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    0,
#line 121 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    -0.1);
#line 10700 "./ESS_IN5_reprate.c"
  mctc1 = coords_neg(mcposasource);
  mcposrsource = rot_apply(mcrotasource, mctc1);
  mcDEBUG_COMPONENT("source", mcposasource, mcrotasource)
  mccomp_posa[1] = mcposasource;
  mccomp_posr[1] = mcposrsource;
  mcNCounter[1]  = mcPCounter[1] = mcP2Counter[1] = 0;
  mcAbsorbProp[1]= 0;
    /* Component Origin. */
  /* Setting parameters for component Origin. */
  SIG_MESSAGE("Origin (Init:SetPar)");
#line 39 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  if("NULL") strncpy(mccOrigin_profile, "NULL" ? "NULL" : "", 16384); else mccOrigin_profile[0]='\0';
#line 39 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccOrigin_percent = 10;
#line 39 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccOrigin_flag_save = 0;
#line 39 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccOrigin_minutes = 0;
#line 10719 "./ESS_IN5_reprate.c"

  SIG_MESSAGE("Origin (Init:Place/Rotate)");
  rot_set_rotation(mcrotaOrigin,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 10726 "./ESS_IN5_reprate.c"
  rot_transpose(mcrotasource, mctr1);
  rot_mul(mcrotaOrigin, mctr1, mcrotrOrigin);
  mcposaOrigin = coords_set(
#line 125 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    0,
#line 125 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    0,
#line 125 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    0);
#line 10736 "./ESS_IN5_reprate.c"
  mctc1 = coords_sub(mcposasource, mcposaOrigin);
  mcposrOrigin = rot_apply(mcrotaOrigin, mctc1);
  mcDEBUG_COMPONENT("Origin", mcposaOrigin, mcrotaOrigin)
  mccomp_posa[2] = mcposaOrigin;
  mccomp_posr[2] = mcposrOrigin;
  mcNCounter[2]  = mcPCounter[2] = mcP2Counter[2] = 0;
  mcAbsorbProp[2]= 0;
    /* Component TOFmoderator_zoom. */
  /* Setting parameters for component TOFmoderator_zoom. */
  SIG_MESSAGE("TOFmoderator_zoom (Init:SetPar)");
#line 128 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  if("TOFmoderator_zoom.dat") strncpy(mccTOFmoderator_zoom_filename, "TOFmoderator_zoom.dat" ? "TOFmoderator_zoom.dat" : "", 16384); else mccTOFmoderator_zoom_filename[0]='\0';
#line 47 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFmoderator_zoom_xmin = -0.05;
#line 47 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFmoderator_zoom_xmax = 0.05;
#line 47 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFmoderator_zoom_ymin = -0.05;
#line 47 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFmoderator_zoom_ymax = 0.05;
#line 128 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFmoderator_zoom_xwidth = 0.12;
#line 129 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFmoderator_zoom_yheight = 0.12;
#line 129 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFmoderator_zoom_tmin = 0;
#line 129 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFmoderator_zoom_tmax = 5000;
#line 48 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFmoderator_zoom_dt = 1.0;
#line 129 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFmoderator_zoom_restore_neutron = 1;
#line 48 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFmoderator_zoom_nowritefile = 0;
#line 10771 "./ESS_IN5_reprate.c"

  SIG_MESSAGE("TOFmoderator_zoom (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 10778 "./ESS_IN5_reprate.c"
  rot_mul(mctr1, mcrotaOrigin, mcrotaTOFmoderator_zoom);
  rot_transpose(mcrotaOrigin, mctr1);
  rot_mul(mcrotaTOFmoderator_zoom, mctr1, mcrotrTOFmoderator_zoom);
  mctc1 = coords_set(
#line 130 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    0,
#line 130 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    0,
#line 130 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    1E-6);
#line 10789 "./ESS_IN5_reprate.c"
  rot_transpose(mcrotaOrigin, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaTOFmoderator_zoom = coords_add(mcposaOrigin, mctc2);
  mctc1 = coords_sub(mcposaOrigin, mcposaTOFmoderator_zoom);
  mcposrTOFmoderator_zoom = rot_apply(mcrotaTOFmoderator_zoom, mctc1);
  mcDEBUG_COMPONENT("TOFmoderator_zoom", mcposaTOFmoderator_zoom, mcrotaTOFmoderator_zoom)
  mccomp_posa[3] = mcposaTOFmoderator_zoom;
  mccomp_posr[3] = mcposrTOFmoderator_zoom;
  mcNCounter[3]  = mcPCounter[3] = mcP2Counter[3] = 0;
  mcAbsorbProp[3]= 0;
    /* Component TOFmoderator. */
  /* Setting parameters for component TOFmoderator. */
  SIG_MESSAGE("TOFmoderator (Init:SetPar)");
#line 133 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  if("TOFmoderator.dat") strncpy(mccTOFmoderator_filename, "TOFmoderator.dat" ? "TOFmoderator.dat" : "", 16384); else mccTOFmoderator_filename[0]='\0';
#line 47 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFmoderator_xmin = -0.05;
#line 47 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFmoderator_xmax = 0.05;
#line 47 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFmoderator_ymin = -0.05;
#line 47 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFmoderator_ymax = 0.05;
#line 133 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFmoderator_xwidth = 0.12;
#line 134 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFmoderator_yheight = 0.12;
#line 134 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFmoderator_tmin = 0;
#line 134 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFmoderator_tmax = 3.0e5;
#line 48 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFmoderator_dt = 1.0;
#line 134 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFmoderator_restore_neutron = 1;
#line 48 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFmoderator_nowritefile = 0;
#line 10827 "./ESS_IN5_reprate.c"

  SIG_MESSAGE("TOFmoderator (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 10834 "./ESS_IN5_reprate.c"
  rot_mul(mctr1, mcrotaOrigin, mcrotaTOFmoderator);
  rot_transpose(mcrotaTOFmoderator_zoom, mctr1);
  rot_mul(mcrotaTOFmoderator, mctr1, mcrotrTOFmoderator);
  mctc1 = coords_set(
#line 135 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    0,
#line 135 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    0,
#line 135 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    1E-6);
#line 10845 "./ESS_IN5_reprate.c"
  rot_transpose(mcrotaOrigin, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaTOFmoderator = coords_add(mcposaOrigin, mctc2);
  mctc1 = coords_sub(mcposaTOFmoderator_zoom, mcposaTOFmoderator);
  mcposrTOFmoderator = rot_apply(mcrotaTOFmoderator, mctc1);
  mcDEBUG_COMPONENT("TOFmoderator", mcposaTOFmoderator, mcrotaTOFmoderator)
  mccomp_posa[4] = mcposaTOFmoderator;
  mccomp_posr[4] = mcposrTOFmoderator;
  mcNCounter[4]  = mcPCounter[4] = mcP2Counter[4] = 0;
  mcAbsorbProp[4]= 0;
    /* Component Lmon_guistart. */
  /* Setting parameters for component Lmon_guistart. */
  SIG_MESSAGE("Lmon_guistart (Init:SetPar)");
#line 138 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  if("Lmon_guistart.dat") strncpy(mccLmon_guistart_filename, "Lmon_guistart.dat" ? "Lmon_guistart.dat" : "", 16384); else mccLmon_guistart_filename[0]='\0';
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmon_guistart_xmin = -0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmon_guistart_xmax = 0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmon_guistart_ymin = -0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmon_guistart_ymax = 0.05;
#line 138 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmon_guistart_xwidth = mcipGUI_w + 0.01;
#line 139 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmon_guistart_yheight = mcipGUI_h + 0.01;
#line 139 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmon_guistart_Lmin = 0;
#line 139 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmon_guistart_Lmax = mcipLmax + 1;
#line 51 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmon_guistart_restore_neutron = 0;
#line 51 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmon_guistart_nowritefile = 0;
#line 10881 "./ESS_IN5_reprate.c"

  SIG_MESSAGE("Lmon_guistart (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 10888 "./ESS_IN5_reprate.c"
  rot_mul(mctr1, mcrotaOrigin, mcrotaLmon_guistart);
  rot_transpose(mcrotaTOFmoderator, mctr1);
  rot_mul(mcrotaLmon_guistart, mctr1, mcrotrLmon_guistart);
  mctc1 = coords_set(
#line 140 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    0,
#line 140 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    0,
#line 140 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    mcipGUI_start -2e-6);
#line 10899 "./ESS_IN5_reprate.c"
  rot_transpose(mcrotaOrigin, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaLmon_guistart = coords_add(mcposaOrigin, mctc2);
  mctc1 = coords_sub(mcposaTOFmoderator, mcposaLmon_guistart);
  mcposrLmon_guistart = rot_apply(mcrotaLmon_guistart, mctc1);
  mcDEBUG_COMPONENT("Lmon_guistart", mcposaLmon_guistart, mcrotaLmon_guistart)
  mccomp_posa[5] = mcposaLmon_guistart;
  mccomp_posr[5] = mcposrLmon_guistart;
  mcNCounter[5]  = mcPCounter[5] = mcP2Counter[5] = 0;
  mcAbsorbProp[5]= 0;
    /* Component Lmon_normalize. */
  /* Setting parameters for component Lmon_normalize. */
  SIG_MESSAGE("Lmon_normalize (Init:SetPar)");
#line 143 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  if("Lmon_guistart_normalize.dat") strncpy(mccLmon_normalize_filename, "Lmon_guistart_normalize.dat" ? "Lmon_guistart_normalize.dat" : "", 16384); else mccLmon_normalize_filename[0]='\0';
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmon_normalize_xmin = -0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmon_normalize_xmax = 0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmon_normalize_ymin = -0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmon_normalize_ymax = 0.05;
#line 143 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmon_normalize_xwidth = 0.10;
#line 144 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmon_normalize_yheight = 0.10;
#line 144 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmon_normalize_Lmin = 0;
#line 144 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmon_normalize_Lmax = 20;
#line 51 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmon_normalize_restore_neutron = 0;
#line 51 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmon_normalize_nowritefile = 0;
#line 10935 "./ESS_IN5_reprate.c"

  SIG_MESSAGE("Lmon_normalize (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 10942 "./ESS_IN5_reprate.c"
  rot_mul(mctr1, mcrotaOrigin, mcrotaLmon_normalize);
  rot_transpose(mcrotaLmon_guistart, mctr1);
  rot_mul(mcrotaLmon_normalize, mctr1, mcrotrLmon_normalize);
  mctc1 = coords_set(
#line 145 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    0,
#line 145 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    0,
#line 145 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    mcipGUI_start -1e-6);
#line 10953 "./ESS_IN5_reprate.c"
  rot_transpose(mcrotaOrigin, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaLmon_normalize = coords_add(mcposaOrigin, mctc2);
  mctc1 = coords_sub(mcposaLmon_guistart, mcposaLmon_normalize);
  mcposrLmon_normalize = rot_apply(mcrotaLmon_normalize, mctc1);
  mcDEBUG_COMPONENT("Lmon_normalize", mcposaLmon_normalize, mcrotaLmon_normalize)
  mccomp_posa[6] = mcposaLmon_normalize;
  mccomp_posr[6] = mcposrLmon_normalize;
  mcNCounter[6]  = mcPCounter[6] = mcP2Counter[6] = 0;
  mcAbsorbProp[6]= 0;
    /* Component Guide1. */
  /* Setting parameters for component Guide1. */
  SIG_MESSAGE("Guide1 (Init:SetPar)");
#line 58 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  if(0) strncpy(mccGuide1_reflect, 0 ? 0 : "", 16384); else mccGuide1_reflect[0]='\0';
#line 150 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccGuide1_w1 = mcipGUI_w;
#line 150 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccGuide1_h1 = mcipGUI_h;
#line 150 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccGuide1_w2 = mcipW1;
#line 150 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccGuide1_h2 = mcipH1;
#line 150 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccGuide1_l = mcipFO1_DIST - mcipGUI_start - mcipGUI_GAP / 2;
#line 150 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccGuide1_R0 = 1;
#line 58 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccGuide1_Qc = 0.0219;
#line 151 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccGuide1_alpha = mcipALPHA;
#line 151 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccGuide1_m = mcipM;
#line 151 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccGuide1_W = 0.003;
#line 10989 "./ESS_IN5_reprate.c"

  SIG_MESSAGE("Guide1 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 10996 "./ESS_IN5_reprate.c"
  rot_mul(mctr1, mcrotaOrigin, mcrotaGuide1);
  rot_transpose(mcrotaLmon_normalize, mctr1);
  rot_mul(mcrotaGuide1, mctr1, mcrotrGuide1);
  mctc1 = coords_set(
#line 152 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    0,
#line 152 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    0,
#line 152 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    mcipGUI_start);
#line 11007 "./ESS_IN5_reprate.c"
  rot_transpose(mcrotaOrigin, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaGuide1 = coords_add(mcposaOrigin, mctc2);
  mctc1 = coords_sub(mcposaLmon_normalize, mcposaGuide1);
  mcposrGuide1 = rot_apply(mcrotaGuide1, mctc1);
  mcDEBUG_COMPONENT("Guide1", mcposaGuide1, mcrotaGuide1)
  mccomp_posa[7] = mcposaGuide1;
  mccomp_posr[7] = mcposrGuide1;
  mcNCounter[7]  = mcPCounter[7] = mcP2Counter[7] = 0;
  mcAbsorbProp[7]= 0;
    /* Component Lmonslow1. */
  /* Setting parameters for component Lmonslow1. */
  SIG_MESSAGE("Lmonslow1 (Init:SetPar)");
#line 155 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  if("Lmonslow1.dat") strncpy(mccLmonslow1_filename, "Lmonslow1.dat" ? "Lmonslow1.dat" : "", 16384); else mccLmonslow1_filename[0]='\0';
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmonslow1_xmin = -0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmonslow1_xmax = 0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmonslow1_ymin = -0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmonslow1_ymax = 0.05;
#line 155 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmonslow1_xwidth = 0.06;
#line 156 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmonslow1_yheight = 0.21;
#line 156 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmonslow1_Lmin = 0;
#line 156 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmonslow1_Lmax = mcipLmax + 1;
#line 51 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmonslow1_restore_neutron = 0;
#line 51 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmonslow1_nowritefile = 0;
#line 11043 "./ESS_IN5_reprate.c"

  SIG_MESSAGE("Lmonslow1 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 11050 "./ESS_IN5_reprate.c"
  rot_mul(mctr1, mcrotaOrigin, mcrotaLmonslow1);
  rot_transpose(mcrotaGuide1, mctr1);
  rot_mul(mcrotaLmonslow1, mctr1, mcrotrLmonslow1);
  mctc1 = coords_set(
#line 157 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    0,
#line 157 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    0,
#line 157 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    mcipFO1_DIST -2e-6);
#line 11061 "./ESS_IN5_reprate.c"
  rot_transpose(mcrotaOrigin, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaLmonslow1 = coords_add(mcposaOrigin, mctc2);
  mctc1 = coords_sub(mcposaGuide1, mcposaLmonslow1);
  mcposrLmonslow1 = rot_apply(mcrotaLmonslow1, mctc1);
  mcDEBUG_COMPONENT("Lmonslow1", mcposaLmonslow1, mcrotaLmonslow1)
  mccomp_posa[8] = mcposaLmonslow1;
  mccomp_posr[8] = mcposrLmonslow1;
  mcNCounter[8]  = mcPCounter[8] = mcP2Counter[8] = 0;
  mcAbsorbProp[8]= 0;
    /* Component PSDslow1. */
  /* Setting parameters for component PSDslow1. */
  SIG_MESSAGE("PSDslow1 (Init:SetPar)");
#line 49 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccPSDslow1_nx = 90;
#line 49 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccPSDslow1_ny = 90;
#line 160 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  if("PSDslow1.dat") strncpy(mccPSDslow1_filename, "PSDslow1.dat" ? "PSDslow1.dat" : "", 16384); else mccPSDslow1_filename[0]='\0';
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccPSDslow1_xmin = -0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccPSDslow1_xmax = 0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccPSDslow1_ymin = -0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccPSDslow1_ymax = 0.05;
#line 160 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccPSDslow1_xwidth = 0.1;
#line 160 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccPSDslow1_yheight = 0.25;
#line 51 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccPSDslow1_restore_neutron = 0;
#line 11095 "./ESS_IN5_reprate.c"

  SIG_MESSAGE("PSDslow1 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 11102 "./ESS_IN5_reprate.c"
  rot_mul(mctr1, mcrotaOrigin, mcrotaPSDslow1);
  rot_transpose(mcrotaLmonslow1, mctr1);
  rot_mul(mcrotaPSDslow1, mctr1, mcrotrPSDslow1);
  mctc1 = coords_set(
#line 161 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    0,
#line 161 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    0,
#line 161 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    mcipFO1_DIST -1e-6);
#line 11113 "./ESS_IN5_reprate.c"
  rot_transpose(mcrotaOrigin, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaPSDslow1 = coords_add(mcposaOrigin, mctc2);
  mctc1 = coords_sub(mcposaLmonslow1, mcposaPSDslow1);
  mcposrPSDslow1 = rot_apply(mcrotaPSDslow1, mctc1);
  mcDEBUG_COMPONENT("PSDslow1", mcposaPSDslow1, mcrotaPSDslow1)
  mccomp_posa[9] = mcposaPSDslow1;
  mccomp_posr[9] = mcposrPSDslow1;
  mcNCounter[9]  = mcPCounter[9] = mcP2Counter[9] = 0;
  mcAbsorbProp[9]= 0;
    /* Component FOchop1. */
  /* Setting parameters for component FOchop1. */
  SIG_MESSAGE("FOchop1 (Init:SetPar)");
#line 164 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccFOchop1_theta_0 = mcipSLOW1_THETA;
#line 164 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccFOchop1_radius = 0.6;
#line 164 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccFOchop1_yheight = 0.20;
#line 164 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccFOchop1_nu = mcipF_slow1;
#line 164 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccFOchop1_nslit = 1;
#line 59 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccFOchop1_jitter = 0;
#line 164 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccFOchop1_delay = t_FO1;
#line 59 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccFOchop1_isfirst = 0;
#line 59 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccFOchop1_n_pulse = 1;
#line 59 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccFOchop1_abs_out = 1;
#line 59 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccFOchop1_phase = 0;
#line 59 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccFOchop1_xwidth = 0;
#line 59 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccFOchop1_verbose = 0;
#line 11153 "./ESS_IN5_reprate.c"

  SIG_MESSAGE("FOchop1 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 11160 "./ESS_IN5_reprate.c"
  rot_mul(mctr1, mcrotaOrigin, mcrotaFOchop1);
  rot_transpose(mcrotaPSDslow1, mctr1);
  rot_mul(mcrotaFOchop1, mctr1, mcrotrFOchop1);
  mctc1 = coords_set(
#line 165 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    0,
#line 165 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    0,
#line 165 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    mcipFO1_DIST);
#line 11171 "./ESS_IN5_reprate.c"
  rot_transpose(mcrotaOrigin, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaFOchop1 = coords_add(mcposaOrigin, mctc2);
  mctc1 = coords_sub(mcposaPSDslow1, mcposaFOchop1);
  mcposrFOchop1 = rot_apply(mcrotaFOchop1, mctc1);
  mcDEBUG_COMPONENT("FOchop1", mcposaFOchop1, mcrotaFOchop1)
  mccomp_posa[10] = mcposaFOchop1;
  mccomp_posr[10] = mcposrFOchop1;
  mcNCounter[10]  = mcPCounter[10] = mcP2Counter[10] = 0;
  mcAbsorbProp[10]= 0;
    /* Component TOFLmon1. */
  /* Setting parameters for component TOFLmon1. */
  SIG_MESSAGE("TOFLmon1 (Init:SetPar)");
#line 168 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  if("TOFLmon1.dat") strncpy(mccTOFLmon1_filename, "TOFLmon1.dat" ? "TOFLmon1.dat" : "", 16384); else mccTOFLmon1_filename[0]='\0';
#line 49 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFLmon1_xmin = -0.05;
#line 49 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFLmon1_xmax = 0.05;
#line 49 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFLmon1_ymin = -0.05;
#line 49 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFLmon1_ymax = 0.05;
#line 169 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFLmon1_xwidth = 0.05;
#line 169 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFLmon1_yheight = 0.21;
#line 169 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFLmon1_Lmin = 0;
#line 170 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFLmon1_Lmax = mcipLmax + 1;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFLmon1_restore_neutron = 0;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFLmon1_nowritefile = 0;
#line 11207 "./ESS_IN5_reprate.c"

  SIG_MESSAGE("TOFLmon1 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 11214 "./ESS_IN5_reprate.c"
  rot_mul(mctr1, mcrotaFOchop1, mcrotaTOFLmon1);
  rot_transpose(mcrotaFOchop1, mctr1);
  rot_mul(mcrotaTOFLmon1, mctr1, mcrotrTOFLmon1);
  mctc1 = coords_set(
#line 171 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    0,
#line 171 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    0,
#line 171 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    1e-6);
#line 11225 "./ESS_IN5_reprate.c"
  rot_transpose(mcrotaFOchop1, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaTOFLmon1 = coords_add(mcposaFOchop1, mctc2);
  mctc1 = coords_sub(mcposaFOchop1, mcposaTOFLmon1);
  mcposrTOFLmon1 = rot_apply(mcrotaTOFLmon1, mctc1);
  mcDEBUG_COMPONENT("TOFLmon1", mcposaTOFLmon1, mcrotaTOFLmon1)
  mccomp_posa[11] = mcposaTOFLmon1;
  mccomp_posr[11] = mcposrTOFLmon1;
  mcNCounter[11]  = mcPCounter[11] = mcP2Counter[11] = 0;
  mcAbsorbProp[11]= 0;
    /* Component Lmon_afterslow1. */
  /* Setting parameters for component Lmon_afterslow1. */
  SIG_MESSAGE("Lmon_afterslow1 (Init:SetPar)");
#line 174 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  if("Lmon_afterslow1.dat") strncpy(mccLmon_afterslow1_filename, "Lmon_afterslow1.dat" ? "Lmon_afterslow1.dat" : "", 16384); else mccLmon_afterslow1_filename[0]='\0';
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmon_afterslow1_xmin = -0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmon_afterslow1_xmax = 0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmon_afterslow1_ymin = -0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmon_afterslow1_ymax = 0.05;
#line 174 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmon_afterslow1_xwidth = 0.06;
#line 175 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmon_afterslow1_yheight = 0.21;
#line 175 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmon_afterslow1_Lmin = 0;
#line 175 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmon_afterslow1_Lmax = mcipLmax + 1;
#line 51 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmon_afterslow1_restore_neutron = 0;
#line 51 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmon_afterslow1_nowritefile = 0;
#line 11261 "./ESS_IN5_reprate.c"

  SIG_MESSAGE("Lmon_afterslow1 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 11268 "./ESS_IN5_reprate.c"
  rot_mul(mctr1, mcrotaFOchop1, mcrotaLmon_afterslow1);
  rot_transpose(mcrotaTOFLmon1, mctr1);
  rot_mul(mcrotaLmon_afterslow1, mctr1, mcrotrLmon_afterslow1);
  mctc1 = coords_set(
#line 176 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    0,
#line 176 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    0,
#line 176 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    2e-6);
#line 11279 "./ESS_IN5_reprate.c"
  rot_transpose(mcrotaFOchop1, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaLmon_afterslow1 = coords_add(mcposaFOchop1, mctc2);
  mctc1 = coords_sub(mcposaTOFLmon1, mcposaLmon_afterslow1);
  mcposrLmon_afterslow1 = rot_apply(mcrotaLmon_afterslow1, mctc1);
  mcDEBUG_COMPONENT("Lmon_afterslow1", mcposaLmon_afterslow1, mcrotaLmon_afterslow1)
  mccomp_posa[12] = mcposaLmon_afterslow1;
  mccomp_posr[12] = mcposrLmon_afterslow1;
  mcNCounter[12]  = mcPCounter[12] = mcP2Counter[12] = 0;
  mcAbsorbProp[12]= 0;
    /* Component PSD_afterslow1. */
  /* Setting parameters for component PSD_afterslow1. */
  SIG_MESSAGE("PSD_afterslow1 (Init:SetPar)");
#line 49 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccPSD_afterslow1_nx = 90;
#line 49 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccPSD_afterslow1_ny = 90;
#line 179 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  if("PSD_afterslow1.dat") strncpy(mccPSD_afterslow1_filename, "PSD_afterslow1.dat" ? "PSD_afterslow1.dat" : "", 16384); else mccPSD_afterslow1_filename[0]='\0';
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccPSD_afterslow1_xmin = -0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccPSD_afterslow1_xmax = 0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccPSD_afterslow1_ymin = -0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccPSD_afterslow1_ymax = 0.05;
#line 179 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccPSD_afterslow1_xwidth = 0.1;
#line 179 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccPSD_afterslow1_yheight = 0.25;
#line 51 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccPSD_afterslow1_restore_neutron = 0;
#line 11313 "./ESS_IN5_reprate.c"

  SIG_MESSAGE("PSD_afterslow1 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 11320 "./ESS_IN5_reprate.c"
  rot_mul(mctr1, mcrotaFOchop1, mcrotaPSD_afterslow1);
  rot_transpose(mcrotaLmon_afterslow1, mctr1);
  rot_mul(mcrotaPSD_afterslow1, mctr1, mcrotrPSD_afterslow1);
  mctc1 = coords_set(
#line 180 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    0,
#line 180 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    0,
#line 180 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    3e-6);
#line 11331 "./ESS_IN5_reprate.c"
  rot_transpose(mcrotaFOchop1, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaPSD_afterslow1 = coords_add(mcposaFOchop1, mctc2);
  mctc1 = coords_sub(mcposaLmon_afterslow1, mcposaPSD_afterslow1);
  mcposrPSD_afterslow1 = rot_apply(mcrotaPSD_afterslow1, mctc1);
  mcDEBUG_COMPONENT("PSD_afterslow1", mcposaPSD_afterslow1, mcrotaPSD_afterslow1)
  mccomp_posa[13] = mcposaPSD_afterslow1;
  mccomp_posr[13] = mcposrPSD_afterslow1;
  mcNCounter[13]  = mcPCounter[13] = mcP2Counter[13] = 0;
  mcAbsorbProp[13]= 0;
    /* Component Guidelong1. */
  /* Setting parameters for component Guidelong1. */
  SIG_MESSAGE("Guidelong1 (Init:SetPar)");
#line 58 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  if(0) strncpy(mccGuidelong1_reflect, 0 ? 0 : "", 16384); else mccGuidelong1_reflect[0]='\0';
#line 183 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccGuidelong1_w1 = mcipW1;
#line 183 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccGuidelong1_h1 = mcipH1;
#line 183 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccGuidelong1_w2 = mcipW2;
#line 183 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccGuidelong1_h2 = mcipH2;
#line 184 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccGuidelong1_l = mcipL_ballistic_begin + mcipGUI_start - mcipFO1_DIST - mcipGUI_GAP / 2;
#line 185 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccGuidelong1_R0 = 1;
#line 58 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccGuidelong1_Qc = 0.0219;
#line 185 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccGuidelong1_alpha = mcipALPHA;
#line 185 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccGuidelong1_m = mcipM;
#line 185 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccGuidelong1_W = 0.003;
#line 11367 "./ESS_IN5_reprate.c"

  SIG_MESSAGE("Guidelong1 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 11374 "./ESS_IN5_reprate.c"
  rot_mul(mctr1, mcrotaFOchop1, mcrotaGuidelong1);
  rot_transpose(mcrotaPSD_afterslow1, mctr1);
  rot_mul(mcrotaGuidelong1, mctr1, mcrotrGuidelong1);
  mctc1 = coords_set(
#line 186 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    0,
#line 186 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    0,
#line 186 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    mcipGUI_GAP / 2);
#line 11385 "./ESS_IN5_reprate.c"
  rot_transpose(mcrotaFOchop1, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaGuidelong1 = coords_add(mcposaFOchop1, mctc2);
  mctc1 = coords_sub(mcposaPSD_afterslow1, mcposaGuidelong1);
  mcposrGuidelong1 = rot_apply(mcrotaGuidelong1, mctc1);
  mcDEBUG_COMPONENT("Guidelong1", mcposaGuidelong1, mcrotaGuidelong1)
  mccomp_posa[14] = mcposaGuidelong1;
  mccomp_posr[14] = mcposrGuidelong1;
  mcNCounter[14]  = mcPCounter[14] = mcP2Counter[14] = 0;
  mcAbsorbProp[14]= 0;
    /* Component Guidelong1b. */
  /* Setting parameters for component Guidelong1b. */
  SIG_MESSAGE("Guidelong1b (Init:SetPar)");
#line 58 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  if(0) strncpy(mccGuidelong1b_reflect, 0 ? 0 : "", 16384); else mccGuidelong1b_reflect[0]='\0';
#line 189 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccGuidelong1b_w1 = mcipW2;
#line 189 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccGuidelong1b_h1 = mcipH2;
#line 189 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccGuidelong1b_w2 = mcipW3;
#line 189 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccGuidelong1b_h2 = mcipH3;
#line 189 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccGuidelong1b_l = mcipLength / 2 - mcipL_ballistic_begin - mcipGUI_start - mcipGUI_GAP / 2;
#line 190 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccGuidelong1b_R0 = 1;
#line 58 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccGuidelong1b_Qc = 0.0219;
#line 190 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccGuidelong1b_alpha = mcipALPHA;
#line 190 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccGuidelong1b_m = mcipM;
#line 190 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccGuidelong1b_W = 0.003;
#line 11421 "./ESS_IN5_reprate.c"

  SIG_MESSAGE("Guidelong1b (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 11428 "./ESS_IN5_reprate.c"
  rot_mul(mctr1, mcrotaOrigin, mcrotaGuidelong1b);
  rot_transpose(mcrotaGuidelong1, mctr1);
  rot_mul(mcrotaGuidelong1b, mctr1, mcrotrGuidelong1b);
  mctc1 = coords_set(
#line 191 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    0,
#line 191 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    0,
#line 191 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    mcipL_ballistic_begin + mcipGUI_start);
#line 11439 "./ESS_IN5_reprate.c"
  rot_transpose(mcrotaOrigin, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaGuidelong1b = coords_add(mcposaOrigin, mctc2);
  mctc1 = coords_sub(mcposaGuidelong1, mcposaGuidelong1b);
  mcposrGuidelong1b = rot_apply(mcrotaGuidelong1b, mctc1);
  mcDEBUG_COMPONENT("Guidelong1b", mcposaGuidelong1b, mcrotaGuidelong1b)
  mccomp_posa[15] = mcposaGuidelong1b;
  mccomp_posr[15] = mcposrGuidelong1b;
  mcNCounter[15]  = mcPCounter[15] = mcP2Counter[15] = 0;
  mcAbsorbProp[15]= 0;
    /* Component Lmon_slow2. */
  /* Setting parameters for component Lmon_slow2. */
  SIG_MESSAGE("Lmon_slow2 (Init:SetPar)");
#line 194 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  if("Lmon_slow2.dat") strncpy(mccLmon_slow2_filename, "Lmon_slow2.dat" ? "Lmon_slow2.dat" : "", 16384); else mccLmon_slow2_filename[0]='\0';
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmon_slow2_xmin = -0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmon_slow2_xmax = 0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmon_slow2_ymin = -0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmon_slow2_ymax = 0.05;
#line 194 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmon_slow2_xwidth = 0.06;
#line 195 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmon_slow2_yheight = 0.21;
#line 195 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmon_slow2_Lmin = 0;
#line 195 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmon_slow2_Lmax = mcipLmax + 1;
#line 51 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmon_slow2_restore_neutron = 0;
#line 51 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmon_slow2_nowritefile = 0;
#line 11475 "./ESS_IN5_reprate.c"

  SIG_MESSAGE("Lmon_slow2 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 11482 "./ESS_IN5_reprate.c"
  rot_mul(mctr1, mcrotaOrigin, mcrotaLmon_slow2);
  rot_transpose(mcrotaGuidelong1b, mctr1);
  rot_mul(mcrotaLmon_slow2, mctr1, mcrotrLmon_slow2);
  mctc1 = coords_set(
#line 196 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    0,
#line 196 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    0,
#line 196 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    mcipLength / 2 -1e-6);
#line 11493 "./ESS_IN5_reprate.c"
  rot_transpose(mcrotaOrigin, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaLmon_slow2 = coords_add(mcposaOrigin, mctc2);
  mctc1 = coords_sub(mcposaGuidelong1b, mcposaLmon_slow2);
  mcposrLmon_slow2 = rot_apply(mcrotaLmon_slow2, mctc1);
  mcDEBUG_COMPONENT("Lmon_slow2", mcposaLmon_slow2, mcrotaLmon_slow2)
  mccomp_posa[16] = mcposaLmon_slow2;
  mccomp_posr[16] = mcposrLmon_slow2;
  mcNCounter[16]  = mcPCounter[16] = mcP2Counter[16] = 0;
  mcAbsorbProp[16]= 0;
    /* Component FOchop2. */
  /* Setting parameters for component FOchop2. */
  SIG_MESSAGE("FOchop2 (Init:SetPar)");
#line 199 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccFOchop2_theta_0 = 155;
#line 199 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccFOchop2_radius = 0.6;
#line 199 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccFOchop2_yheight = 0.6;
#line 199 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccFOchop2_nu = mcipF_slow2;
#line 199 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccFOchop2_nslit = 1;
#line 59 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccFOchop2_jitter = 0;
#line 199 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccFOchop2_delay = t_FO2;
#line 59 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccFOchop2_isfirst = 0;
#line 59 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccFOchop2_n_pulse = 1;
#line 59 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccFOchop2_abs_out = 1;
#line 59 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccFOchop2_phase = 0;
#line 59 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccFOchop2_xwidth = 0;
#line 59 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccFOchop2_verbose = 0;
#line 11533 "./ESS_IN5_reprate.c"

  SIG_MESSAGE("FOchop2 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 11540 "./ESS_IN5_reprate.c"
  rot_mul(mctr1, mcrotaOrigin, mcrotaFOchop2);
  rot_transpose(mcrotaLmon_slow2, mctr1);
  rot_mul(mcrotaFOchop2, mctr1, mcrotrFOchop2);
  mctc1 = coords_set(
#line 200 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    0,
#line 200 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    0.1,
#line 200 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    mcipLength / 2);
#line 11551 "./ESS_IN5_reprate.c"
  rot_transpose(mcrotaOrigin, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaFOchop2 = coords_add(mcposaOrigin, mctc2);
  mctc1 = coords_sub(mcposaLmon_slow2, mcposaFOchop2);
  mcposrFOchop2 = rot_apply(mcrotaFOchop2, mctc1);
  mcDEBUG_COMPONENT("FOchop2", mcposaFOchop2, mcrotaFOchop2)
  mccomp_posa[17] = mcposaFOchop2;
  mccomp_posr[17] = mcposrFOchop2;
  mcNCounter[17]  = mcPCounter[17] = mcP2Counter[17] = 0;
  mcAbsorbProp[17]= 0;
    /* Component Fastchop1. */
  /* Setting parameters for component Fastchop1. */
  SIG_MESSAGE("Fastchop1 (Init:SetPar)");
#line 203 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccFastchop1_theta_0 = mcipTHETA_fast1;
#line 203 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccFastchop1_radius = 0.5;
#line 203 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccFastchop1_yheight = 0.5;
#line 203 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccFastchop1_nu = mcipF_fast1;
#line 203 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccFastchop1_nslit = mcipN_fast;
#line 59 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccFastchop1_jitter = 0;
#line 203 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccFastchop1_delay = t_fast1;
#line 59 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccFastchop1_isfirst = 0;
#line 59 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccFastchop1_n_pulse = 1;
#line 59 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccFastchop1_abs_out = 1;
#line 59 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccFastchop1_phase = 0;
#line 59 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccFastchop1_xwidth = 0;
#line 59 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccFastchop1_verbose = 0;
#line 11591 "./ESS_IN5_reprate.c"

  SIG_MESSAGE("Fastchop1 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 11598 "./ESS_IN5_reprate.c"
  rot_mul(mctr1, mcrotaOrigin, mcrotaFastchop1);
  rot_transpose(mcrotaFOchop2, mctr1);
  rot_mul(mcrotaFastchop1, mctr1, mcrotrFastchop1);
  mctc1 = coords_set(
#line 204 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    0,
#line 204 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    0.1,
#line 204 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    mcipLength / 2 + mcipGUI_GAP);
#line 11609 "./ESS_IN5_reprate.c"
  rot_transpose(mcrotaOrigin, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaFastchop1 = coords_add(mcposaOrigin, mctc2);
  mctc1 = coords_sub(mcposaFOchop2, mcposaFastchop1);
  mcposrFastchop1 = rot_apply(mcrotaFastchop1, mctc1);
  mcDEBUG_COMPONENT("Fastchop1", mcposaFastchop1, mcrotaFastchop1)
  mccomp_posa[18] = mcposaFastchop1;
  mccomp_posr[18] = mcposrFastchop1;
  mcNCounter[18]  = mcPCounter[18] = mcP2Counter[18] = 0;
  mcAbsorbProp[18]= 0;
    /* Component PSD_afterslow2. */
  /* Setting parameters for component PSD_afterslow2. */
  SIG_MESSAGE("PSD_afterslow2 (Init:SetPar)");
#line 49 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccPSD_afterslow2_nx = 90;
#line 49 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccPSD_afterslow2_ny = 90;
#line 207 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  if("PSD_afterslow2.dat") strncpy(mccPSD_afterslow2_filename, "PSD_afterslow2.dat" ? "PSD_afterslow2.dat" : "", 16384); else mccPSD_afterslow2_filename[0]='\0';
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccPSD_afterslow2_xmin = -0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccPSD_afterslow2_xmax = 0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccPSD_afterslow2_ymin = -0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccPSD_afterslow2_ymax = 0.05;
#line 207 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccPSD_afterslow2_xwidth = 0.1;
#line 207 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccPSD_afterslow2_yheight = 0.25;
#line 51 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccPSD_afterslow2_restore_neutron = 0;
#line 11643 "./ESS_IN5_reprate.c"

  SIG_MESSAGE("PSD_afterslow2 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 11650 "./ESS_IN5_reprate.c"
  rot_mul(mctr1, mcrotaFastchop1, mcrotaPSD_afterslow2);
  rot_transpose(mcrotaFastchop1, mctr1);
  rot_mul(mcrotaPSD_afterslow2, mctr1, mcrotrPSD_afterslow2);
  mctc1 = coords_set(
#line 208 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    0,
#line 208 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    -0.1,
#line 208 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    1e-6);
#line 11661 "./ESS_IN5_reprate.c"
  rot_transpose(mcrotaFastchop1, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaPSD_afterslow2 = coords_add(mcposaFastchop1, mctc2);
  mctc1 = coords_sub(mcposaFastchop1, mcposaPSD_afterslow2);
  mcposrPSD_afterslow2 = rot_apply(mcrotaPSD_afterslow2, mctc1);
  mcDEBUG_COMPONENT("PSD_afterslow2", mcposaPSD_afterslow2, mcrotaPSD_afterslow2)
  mccomp_posa[19] = mcposaPSD_afterslow2;
  mccomp_posr[19] = mcposrPSD_afterslow2;
  mcNCounter[19]  = mcPCounter[19] = mcP2Counter[19] = 0;
  mcAbsorbProp[19]= 0;
    /* Component Lmon_afterslow2. */
  /* Setting parameters for component Lmon_afterslow2. */
  SIG_MESSAGE("Lmon_afterslow2 (Init:SetPar)");
#line 211 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  if("Lmon_afterslow2.dat") strncpy(mccLmon_afterslow2_filename, "Lmon_afterslow2.dat" ? "Lmon_afterslow2.dat" : "", 16384); else mccLmon_afterslow2_filename[0]='\0';
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmon_afterslow2_xmin = -0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmon_afterslow2_xmax = 0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmon_afterslow2_ymin = -0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmon_afterslow2_ymax = 0.05;
#line 211 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmon_afterslow2_xwidth = 0.06;
#line 212 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmon_afterslow2_yheight = 0.21;
#line 212 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmon_afterslow2_Lmin = 0;
#line 212 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmon_afterslow2_Lmax = mcipLmax + 1;
#line 51 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmon_afterslow2_restore_neutron = 0;
#line 51 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmon_afterslow2_nowritefile = 0;
#line 11697 "./ESS_IN5_reprate.c"

  SIG_MESSAGE("Lmon_afterslow2 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 11704 "./ESS_IN5_reprate.c"
  rot_mul(mctr1, mcrotaPSD_afterslow2, mcrotaLmon_afterslow2);
  rot_transpose(mcrotaPSD_afterslow2, mctr1);
  rot_mul(mcrotaLmon_afterslow2, mctr1, mcrotrLmon_afterslow2);
  mctc1 = coords_set(
#line 213 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    0,
#line 213 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    0,
#line 213 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    1e-6);
#line 11715 "./ESS_IN5_reprate.c"
  rot_transpose(mcrotaPSD_afterslow2, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaLmon_afterslow2 = coords_add(mcposaPSD_afterslow2, mctc2);
  mctc1 = coords_sub(mcposaPSD_afterslow2, mcposaLmon_afterslow2);
  mcposrLmon_afterslow2 = rot_apply(mcrotaLmon_afterslow2, mctc1);
  mcDEBUG_COMPONENT("Lmon_afterslow2", mcposaLmon_afterslow2, mcrotaLmon_afterslow2)
  mccomp_posa[20] = mcposaLmon_afterslow2;
  mccomp_posr[20] = mcposrLmon_afterslow2;
  mcNCounter[20]  = mcPCounter[20] = mcP2Counter[20] = 0;
  mcAbsorbProp[20]= 0;
    /* Component TOFL_afterslow2. */
  /* Setting parameters for component TOFL_afterslow2. */
  SIG_MESSAGE("TOFL_afterslow2 (Init:SetPar)");
#line 216 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  if("TOFL_afterslow2.dat") strncpy(mccTOFL_afterslow2_filename, "TOFL_afterslow2.dat" ? "TOFL_afterslow2.dat" : "", 16384); else mccTOFL_afterslow2_filename[0]='\0';
#line 49 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFL_afterslow2_xmin = -0.05;
#line 49 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFL_afterslow2_xmax = 0.05;
#line 49 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFL_afterslow2_ymin = -0.05;
#line 49 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFL_afterslow2_ymax = 0.05;
#line 217 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFL_afterslow2_xwidth = 0.05;
#line 217 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFL_afterslow2_yheight = 0.21;
#line 217 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFL_afterslow2_Lmin = 0;
#line 218 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFL_afterslow2_Lmax = mcipLmax + 1;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFL_afterslow2_restore_neutron = 0;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFL_afterslow2_nowritefile = 0;
#line 11751 "./ESS_IN5_reprate.c"

  SIG_MESSAGE("TOFL_afterslow2 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 11758 "./ESS_IN5_reprate.c"
  rot_mul(mctr1, mcrotaLmon_afterslow2, mcrotaTOFL_afterslow2);
  rot_transpose(mcrotaLmon_afterslow2, mctr1);
  rot_mul(mcrotaTOFL_afterslow2, mctr1, mcrotrTOFL_afterslow2);
  mctc1 = coords_set(
#line 219 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    0,
#line 219 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    0,
#line 219 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    1e-6);
#line 11769 "./ESS_IN5_reprate.c"
  rot_transpose(mcrotaLmon_afterslow2, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaTOFL_afterslow2 = coords_add(mcposaLmon_afterslow2, mctc2);
  mctc1 = coords_sub(mcposaLmon_afterslow2, mcposaTOFL_afterslow2);
  mcposrTOFL_afterslow2 = rot_apply(mcrotaTOFL_afterslow2, mctc1);
  mcDEBUG_COMPONENT("TOFL_afterslow2", mcposaTOFL_afterslow2, mcrotaTOFL_afterslow2)
  mccomp_posa[21] = mcposaTOFL_afterslow2;
  mccomp_posr[21] = mcposrTOFL_afterslow2;
  mcNCounter[21]  = mcPCounter[21] = mcP2Counter[21] = 0;
  mcAbsorbProp[21]= 0;
    /* Component Guidelong2. */
  /* Setting parameters for component Guidelong2. */
  SIG_MESSAGE("Guidelong2 (Init:SetPar)");
#line 58 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  if(0) strncpy(mccGuidelong2_reflect, 0 ? 0 : "", 16384); else mccGuidelong2_reflect[0]='\0';
#line 222 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccGuidelong2_w1 = mcipW3;
#line 222 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccGuidelong2_h1 = mcipH3;
#line 222 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccGuidelong2_w2 = mcipW4;
#line 222 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccGuidelong2_h2 = mcipH4;
#line 222 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccGuidelong2_l = mcipLength / 2 -2 * mcipGUI_GAP - mcipL_ballistic_end;
#line 223 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccGuidelong2_R0 = 1;
#line 58 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccGuidelong2_Qc = 0.0219;
#line 223 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccGuidelong2_alpha = mcipALPHA;
#line 223 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccGuidelong2_m = mcipM;
#line 223 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccGuidelong2_W = 0.003;
#line 11805 "./ESS_IN5_reprate.c"

  SIG_MESSAGE("Guidelong2 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 11812 "./ESS_IN5_reprate.c"
  rot_mul(mctr1, mcrotaFastchop1, mcrotaGuidelong2);
  rot_transpose(mcrotaTOFL_afterslow2, mctr1);
  rot_mul(mcrotaGuidelong2, mctr1, mcrotrGuidelong2);
  mctc1 = coords_set(
#line 224 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    0,
#line 224 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    -0.1,
#line 224 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    mcipGUI_GAP / 2);
#line 11823 "./ESS_IN5_reprate.c"
  rot_transpose(mcrotaFastchop1, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaGuidelong2 = coords_add(mcposaFastchop1, mctc2);
  mctc1 = coords_sub(mcposaTOFL_afterslow2, mcposaGuidelong2);
  mcposrGuidelong2 = rot_apply(mcrotaGuidelong2, mctc1);
  mcDEBUG_COMPONENT("Guidelong2", mcposaGuidelong2, mcrotaGuidelong2)
  mccomp_posa[22] = mcposaGuidelong2;
  mccomp_posr[22] = mcposrGuidelong2;
  mcNCounter[22]  = mcPCounter[22] = mcP2Counter[22] = 0;
  mcAbsorbProp[22]= 0;
    /* Component Lmon_beforeballistic. */
  /* Setting parameters for component Lmon_beforeballistic. */
  SIG_MESSAGE("Lmon_beforeballistic (Init:SetPar)");
#line 227 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  if("Lmon_before_ballistic.dat") strncpy(mccLmon_beforeballistic_filename, "Lmon_before_ballistic.dat" ? "Lmon_before_ballistic.dat" : "", 16384); else mccLmon_beforeballistic_filename[0]='\0';
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmon_beforeballistic_xmin = -0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmon_beforeballistic_xmax = 0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmon_beforeballistic_ymin = -0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmon_beforeballistic_ymax = 0.05;
#line 227 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmon_beforeballistic_xwidth = 0.06;
#line 228 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmon_beforeballistic_yheight = 0.18;
#line 228 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmon_beforeballistic_Lmin = 0;
#line 228 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmon_beforeballistic_Lmax = mcipLmax + 1;
#line 51 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmon_beforeballistic_restore_neutron = 0;
#line 51 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmon_beforeballistic_nowritefile = 0;
#line 11859 "./ESS_IN5_reprate.c"

  SIG_MESSAGE("Lmon_beforeballistic (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 11866 "./ESS_IN5_reprate.c"
  rot_mul(mctr1, mcrotaGuidelong2, mcrotaLmon_beforeballistic);
  rot_transpose(mcrotaGuidelong2, mctr1);
  rot_mul(mcrotaLmon_beforeballistic, mctr1, mcrotrLmon_beforeballistic);
  mctc1 = coords_set(
#line 229 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    0,
#line 229 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    0,
#line 229 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    mcipLength / 2 -2 * mcipGUI_GAP - mcipL_ballistic_end + 1e-6);
#line 11877 "./ESS_IN5_reprate.c"
  rot_transpose(mcrotaGuidelong2, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaLmon_beforeballistic = coords_add(mcposaGuidelong2, mctc2);
  mctc1 = coords_sub(mcposaGuidelong2, mcposaLmon_beforeballistic);
  mcposrLmon_beforeballistic = rot_apply(mcrotaLmon_beforeballistic, mctc1);
  mcDEBUG_COMPONENT("Lmon_beforeballistic", mcposaLmon_beforeballistic, mcrotaLmon_beforeballistic)
  mccomp_posa[23] = mcposaLmon_beforeballistic;
  mccomp_posr[23] = mcposrLmon_beforeballistic;
  mcNCounter[23]  = mcPCounter[23] = mcP2Counter[23] = 0;
  mcAbsorbProp[23]= 0;
    /* Component PSD_beforeballistic. */
  /* Setting parameters for component PSD_beforeballistic. */
  SIG_MESSAGE("PSD_beforeballistic (Init:SetPar)");
#line 49 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccPSD_beforeballistic_nx = 90;
#line 49 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccPSD_beforeballistic_ny = 90;
#line 232 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  if("PSD_beforeballistic.dat") strncpy(mccPSD_beforeballistic_filename, "PSD_beforeballistic.dat" ? "PSD_beforeballistic.dat" : "", 16384); else mccPSD_beforeballistic_filename[0]='\0';
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccPSD_beforeballistic_xmin = -0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccPSD_beforeballistic_xmax = 0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccPSD_beforeballistic_ymin = -0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccPSD_beforeballistic_ymax = 0.05;
#line 232 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccPSD_beforeballistic_xwidth = 0.1;
#line 232 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccPSD_beforeballistic_yheight = 0.25;
#line 51 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccPSD_beforeballistic_restore_neutron = 0;
#line 11911 "./ESS_IN5_reprate.c"

  SIG_MESSAGE("PSD_beforeballistic (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 11918 "./ESS_IN5_reprate.c"
  rot_mul(mctr1, mcrotaLmon_beforeballistic, mcrotaPSD_beforeballistic);
  rot_transpose(mcrotaLmon_beforeballistic, mctr1);
  rot_mul(mcrotaPSD_beforeballistic, mctr1, mcrotrPSD_beforeballistic);
  mctc1 = coords_set(
#line 233 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    0,
#line 233 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    0,
#line 233 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    1e-6);
#line 11929 "./ESS_IN5_reprate.c"
  rot_transpose(mcrotaLmon_beforeballistic, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaPSD_beforeballistic = coords_add(mcposaLmon_beforeballistic, mctc2);
  mctc1 = coords_sub(mcposaLmon_beforeballistic, mcposaPSD_beforeballistic);
  mcposrPSD_beforeballistic = rot_apply(mcrotaPSD_beforeballistic, mctc1);
  mcDEBUG_COMPONENT("PSD_beforeballistic", mcposaPSD_beforeballistic, mcrotaPSD_beforeballistic)
  mccomp_posa[24] = mcposaPSD_beforeballistic;
  mccomp_posr[24] = mcposrPSD_beforeballistic;
  mcNCounter[24]  = mcPCounter[24] = mcP2Counter[24] = 0;
  mcAbsorbProp[24]= 0;
    /* Component Guidelong2a. */
  /* Setting parameters for component Guidelong2a. */
  SIG_MESSAGE("Guidelong2a (Init:SetPar)");
#line 58 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  if(0) strncpy(mccGuidelong2a_reflect, 0 ? 0 : "", 16384); else mccGuidelong2a_reflect[0]='\0';
#line 236 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccGuidelong2a_w1 = mcipW4;
#line 236 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccGuidelong2a_h1 = mcipH4;
#line 236 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccGuidelong2a_w2 = mcipW_chop;
#line 236 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccGuidelong2a_h2 = mcipH_chop;
#line 236 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccGuidelong2a_l = mcipL_ballistic_end;
#line 237 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccGuidelong2a_R0 = 1;
#line 58 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccGuidelong2a_Qc = 0.0219;
#line 237 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccGuidelong2a_alpha = mcipALPHA;
#line 237 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccGuidelong2a_m = mcipM;
#line 237 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccGuidelong2a_W = 0.003;
#line 11965 "./ESS_IN5_reprate.c"

  SIG_MESSAGE("Guidelong2a (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 11972 "./ESS_IN5_reprate.c"
  rot_mul(mctr1, mcrotaPSD_beforeballistic, mcrotaGuidelong2a);
  rot_transpose(mcrotaPSD_beforeballistic, mctr1);
  rot_mul(mcrotaGuidelong2a, mctr1, mcrotrGuidelong2a);
  mctc1 = coords_set(
#line 238 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    0,
#line 238 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    0,
#line 238 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    1e-6);
#line 11983 "./ESS_IN5_reprate.c"
  rot_transpose(mcrotaPSD_beforeballistic, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaGuidelong2a = coords_add(mcposaPSD_beforeballistic, mctc2);
  mctc1 = coords_sub(mcposaPSD_beforeballistic, mcposaGuidelong2a);
  mcposrGuidelong2a = rot_apply(mcrotaGuidelong2a, mctc1);
  mcDEBUG_COMPONENT("Guidelong2a", mcposaGuidelong2a, mcrotaGuidelong2a)
  mccomp_posa[25] = mcposaGuidelong2a;
  mccomp_posr[25] = mcposrGuidelong2a;
  mcNCounter[25]  = mcPCounter[25] = mcP2Counter[25] = 0;
  mcAbsorbProp[25]= 0;
    /* Component Lmonfast2. */
  /* Setting parameters for component Lmonfast2. */
  SIG_MESSAGE("Lmonfast2 (Init:SetPar)");
#line 241 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  if("Lmonfast2.dat") strncpy(mccLmonfast2_filename, "Lmonfast2.dat" ? "Lmonfast2.dat" : "", 16384); else mccLmonfast2_filename[0]='\0';
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmonfast2_xmin = -0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmonfast2_xmax = 0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmonfast2_ymin = -0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmonfast2_ymax = 0.05;
#line 241 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmonfast2_xwidth = 0.06;
#line 242 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmonfast2_yheight = 0.18;
#line 242 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmonfast2_Lmin = 0;
#line 242 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmonfast2_Lmax = mcipLmax + 1;
#line 51 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmonfast2_restore_neutron = 0;
#line 51 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmonfast2_nowritefile = 0;
#line 12019 "./ESS_IN5_reprate.c"

  SIG_MESSAGE("Lmonfast2 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 12026 "./ESS_IN5_reprate.c"
  rot_mul(mctr1, mcrotaGuidelong2a, mcrotaLmonfast2);
  rot_transpose(mcrotaGuidelong2a, mctr1);
  rot_mul(mcrotaLmonfast2, mctr1, mcrotrLmonfast2);
  mctc1 = coords_set(
#line 243 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    0,
#line 243 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    0,
#line 243 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    mcipL_ballistic_end + 1e-6);
#line 12037 "./ESS_IN5_reprate.c"
  rot_transpose(mcrotaGuidelong2a, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaLmonfast2 = coords_add(mcposaGuidelong2a, mctc2);
  mctc1 = coords_sub(mcposaGuidelong2a, mcposaLmonfast2);
  mcposrLmonfast2 = rot_apply(mcrotaLmonfast2, mctc1);
  mcDEBUG_COMPONENT("Lmonfast2", mcposaLmonfast2, mcrotaLmonfast2)
  mccomp_posa[26] = mcposaLmonfast2;
  mccomp_posr[26] = mcposrLmonfast2;
  mcNCounter[26]  = mcPCounter[26] = mcP2Counter[26] = 0;
  mcAbsorbProp[26]= 0;
    /* Component Lmonfast2_zoom. */
  /* Setting parameters for component Lmonfast2_zoom. */
  SIG_MESSAGE("Lmonfast2_zoom (Init:SetPar)");
#line 246 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  if("Lmonfast2_zoom.dat") strncpy(mccLmonfast2_zoom_filename, "Lmonfast2_zoom.dat" ? "Lmonfast2_zoom.dat" : "", 16384); else mccLmonfast2_zoom_filename[0]='\0';
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmonfast2_zoom_xmin = -0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmonfast2_zoom_xmax = 0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmonfast2_zoom_ymin = -0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmonfast2_zoom_ymax = 0.05;
#line 246 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmonfast2_zoom_xwidth = 0.06;
#line 247 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmonfast2_zoom_yheight = 0.18;
#line 247 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmonfast2_zoom_Lmin = mciplambda0 -0.2;
#line 247 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmonfast2_zoom_Lmax = mciplambda0 + 0.2;
#line 51 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmonfast2_zoom_restore_neutron = 0;
#line 51 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmonfast2_zoom_nowritefile = 0;
#line 12073 "./ESS_IN5_reprate.c"

  SIG_MESSAGE("Lmonfast2_zoom (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 12080 "./ESS_IN5_reprate.c"
  rot_mul(mctr1, mcrotaLmonfast2, mcrotaLmonfast2_zoom);
  rot_transpose(mcrotaLmonfast2, mctr1);
  rot_mul(mcrotaLmonfast2_zoom, mctr1, mcrotrLmonfast2_zoom);
  mctc1 = coords_set(
#line 248 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    0,
#line 248 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    0,
#line 248 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    1e-6);
#line 12091 "./ESS_IN5_reprate.c"
  rot_transpose(mcrotaLmonfast2, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaLmonfast2_zoom = coords_add(mcposaLmonfast2, mctc2);
  mctc1 = coords_sub(mcposaLmonfast2, mcposaLmonfast2_zoom);
  mcposrLmonfast2_zoom = rot_apply(mcrotaLmonfast2_zoom, mctc1);
  mcDEBUG_COMPONENT("Lmonfast2_zoom", mcposaLmonfast2_zoom, mcrotaLmonfast2_zoom)
  mccomp_posa[27] = mcposaLmonfast2_zoom;
  mccomp_posr[27] = mcposrLmonfast2_zoom;
  mcNCounter[27]  = mcPCounter[27] = mcP2Counter[27] = 0;
  mcAbsorbProp[27]= 0;
    /* Component TOFLfast2. */
  /* Setting parameters for component TOFLfast2. */
  SIG_MESSAGE("TOFLfast2 (Init:SetPar)");
#line 251 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  if("TOFLfast2.dat") strncpy(mccTOFLfast2_filename, "TOFLfast2.dat" ? "TOFLfast2.dat" : "", 16384); else mccTOFLfast2_filename[0]='\0';
#line 49 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFLfast2_xmin = -0.05;
#line 49 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFLfast2_xmax = 0.05;
#line 49 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFLfast2_ymin = -0.05;
#line 49 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFLfast2_ymax = 0.05;
#line 252 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFLfast2_xwidth = 0.05;
#line 252 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFLfast2_yheight = 0.12;
#line 252 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFLfast2_Lmin = 0;
#line 253 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFLfast2_Lmax = mcipLmax + 1;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFLfast2_restore_neutron = 0;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFLfast2_nowritefile = 0;
#line 12127 "./ESS_IN5_reprate.c"

  SIG_MESSAGE("TOFLfast2 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 12134 "./ESS_IN5_reprate.c"
  rot_mul(mctr1, mcrotaLmonfast2_zoom, mcrotaTOFLfast2);
  rot_transpose(mcrotaLmonfast2_zoom, mctr1);
  rot_mul(mcrotaTOFLfast2, mctr1, mcrotrTOFLfast2);
  mctc1 = coords_set(
#line 254 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    0,
#line 254 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    0,
#line 254 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    1e-6);
#line 12145 "./ESS_IN5_reprate.c"
  rot_transpose(mcrotaLmonfast2_zoom, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaTOFLfast2 = coords_add(mcposaLmonfast2_zoom, mctc2);
  mctc1 = coords_sub(mcposaLmonfast2_zoom, mcposaTOFLfast2);
  mcposrTOFLfast2 = rot_apply(mcrotaTOFLfast2, mctc1);
  mcDEBUG_COMPONENT("TOFLfast2", mcposaTOFLfast2, mcrotaTOFLfast2)
  mccomp_posa[28] = mcposaTOFLfast2;
  mccomp_posr[28] = mcposrTOFLfast2;
  mcNCounter[28]  = mcPCounter[28] = mcP2Counter[28] = 0;
  mcAbsorbProp[28]= 0;
    /* Component TOFLfast2zoom. */
  /* Setting parameters for component TOFLfast2zoom. */
  SIG_MESSAGE("TOFLfast2zoom (Init:SetPar)");
#line 257 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  if("TOFLfast2_zoom.dat") strncpy(mccTOFLfast2zoom_filename, "TOFLfast2_zoom.dat" ? "TOFLfast2_zoom.dat" : "", 16384); else mccTOFLfast2zoom_filename[0]='\0';
#line 49 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFLfast2zoom_xmin = -0.05;
#line 49 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFLfast2zoom_xmax = 0.05;
#line 49 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFLfast2zoom_ymin = -0.05;
#line 49 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFLfast2zoom_ymax = 0.05;
#line 258 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFLfast2zoom_xwidth = 0.05;
#line 258 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFLfast2zoom_yheight = 0.12;
#line 258 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFLfast2zoom_Lmin = mciplambda0 -0.2;
#line 259 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFLfast2zoom_Lmax = mciplambda0 + 0.2;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFLfast2zoom_restore_neutron = 0;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFLfast2zoom_nowritefile = 0;
#line 12181 "./ESS_IN5_reprate.c"

  SIG_MESSAGE("TOFLfast2zoom (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 12188 "./ESS_IN5_reprate.c"
  rot_mul(mctr1, mcrotaTOFLfast2, mcrotaTOFLfast2zoom);
  rot_transpose(mcrotaTOFLfast2, mctr1);
  rot_mul(mcrotaTOFLfast2zoom, mctr1, mcrotrTOFLfast2zoom);
  mctc1 = coords_set(
#line 260 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    0,
#line 260 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    0,
#line 260 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    1e-6);
#line 12199 "./ESS_IN5_reprate.c"
  rot_transpose(mcrotaTOFLfast2, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaTOFLfast2zoom = coords_add(mcposaTOFLfast2, mctc2);
  mctc1 = coords_sub(mcposaTOFLfast2, mcposaTOFLfast2zoom);
  mcposrTOFLfast2zoom = rot_apply(mcrotaTOFLfast2zoom, mctc1);
  mcDEBUG_COMPONENT("TOFLfast2zoom", mcposaTOFLfast2zoom, mcrotaTOFLfast2zoom)
  mccomp_posa[29] = mcposaTOFLfast2zoom;
  mccomp_posr[29] = mcposrTOFLfast2zoom;
  mcNCounter[29]  = mcPCounter[29] = mcP2Counter[29] = 0;
  mcAbsorbProp[29]= 0;
    /* Component PSDfast2. */
  /* Setting parameters for component PSDfast2. */
  SIG_MESSAGE("PSDfast2 (Init:SetPar)");
#line 49 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccPSDfast2_nx = 90;
#line 49 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccPSDfast2_ny = 90;
#line 263 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  if("PSDfast2.dat") strncpy(mccPSDfast2_filename, "PSDfast2.dat" ? "PSDfast2.dat" : "", 16384); else mccPSDfast2_filename[0]='\0';
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccPSDfast2_xmin = -0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccPSDfast2_xmax = 0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccPSDfast2_ymin = -0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccPSDfast2_ymax = 0.05;
#line 263 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccPSDfast2_xwidth = 0.1;
#line 263 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccPSDfast2_yheight = 0.25;
#line 51 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccPSDfast2_restore_neutron = 0;
#line 12233 "./ESS_IN5_reprate.c"

  SIG_MESSAGE("PSDfast2 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 12240 "./ESS_IN5_reprate.c"
  rot_mul(mctr1, mcrotaTOFLfast2zoom, mcrotaPSDfast2);
  rot_transpose(mcrotaTOFLfast2zoom, mctr1);
  rot_mul(mcrotaPSDfast2, mctr1, mcrotrPSDfast2);
  mctc1 = coords_set(
#line 264 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    0,
#line 264 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    0,
#line 264 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    1e-6);
#line 12251 "./ESS_IN5_reprate.c"
  rot_transpose(mcrotaTOFLfast2zoom, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaPSDfast2 = coords_add(mcposaTOFLfast2zoom, mctc2);
  mctc1 = coords_sub(mcposaTOFLfast2zoom, mcposaPSDfast2);
  mcposrPSDfast2 = rot_apply(mcrotaPSDfast2, mctc1);
  mcDEBUG_COMPONENT("PSDfast2", mcposaPSDfast2, mcrotaPSDfast2)
  mccomp_posa[30] = mcposaPSDfast2;
  mccomp_posr[30] = mcposrPSDfast2;
  mcNCounter[30]  = mcPCounter[30] = mcP2Counter[30] = 0;
  mcAbsorbProp[30]= 0;
    /* Component Fastchop2. */
  /* Setting parameters for component Fastchop2. */
  SIG_MESSAGE("Fastchop2 (Init:SetPar)");
#line 267 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccFastchop2_theta_0 = mcipFAST_THETA;
#line 267 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccFastchop2_radius = 0.35;
#line 267 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccFastchop2_yheight = 0.35;
#line 267 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccFastchop2_nu = mcipF_fast2;
#line 267 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccFastchop2_nslit = mcipN_fast;
#line 59 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccFastchop2_jitter = 0;
#line 267 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccFastchop2_delay = t_fast2;
#line 59 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccFastchop2_isfirst = 0;
#line 59 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccFastchop2_n_pulse = 1;
#line 59 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccFastchop2_abs_out = 1;
#line 59 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccFastchop2_phase = 0;
#line 59 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccFastchop2_xwidth = 0;
#line 59 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccFastchop2_verbose = 0;
#line 12291 "./ESS_IN5_reprate.c"

  SIG_MESSAGE("Fastchop2 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 12298 "./ESS_IN5_reprate.c"
  rot_mul(mctr1, mcrotaOrigin, mcrotaFastchop2);
  rot_transpose(mcrotaPSDfast2, mctr1);
  rot_mul(mcrotaFastchop2, mctr1, mcrotrFastchop2);
  mctc1 = coords_set(
#line 268 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    0,
#line 268 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    0.04,
#line 268 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    mcipLength);
#line 12309 "./ESS_IN5_reprate.c"
  rot_transpose(mcrotaOrigin, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaFastchop2 = coords_add(mcposaOrigin, mctc2);
  mctc1 = coords_sub(mcposaPSDfast2, mcposaFastchop2);
  mcposrFastchop2 = rot_apply(mcrotaFastchop2, mctc1);
  mcDEBUG_COMPONENT("Fastchop2", mcposaFastchop2, mcrotaFastchop2)
  mccomp_posa[31] = mcposaFastchop2;
  mccomp_posr[31] = mcposrFastchop2;
  mcNCounter[31]  = mcPCounter[31] = mcP2Counter[31] = 0;
  mcAbsorbProp[31]= 0;
    /* Component Fastchop2counter. */
  /* Setting parameters for component Fastchop2counter. */
  SIG_MESSAGE("Fastchop2counter (Init:SetPar)");
#line 271 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccFastchop2counter_theta_0 = mcipFAST_THETA;
#line 271 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccFastchop2counter_radius = 0.35;
#line 271 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccFastchop2counter_yheight = 0.35;
#line 271 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccFastchop2counter_nu = - mcipF_fast2;
#line 271 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccFastchop2counter_nslit = mcipN_fast;
#line 59 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccFastchop2counter_jitter = 0;
#line 271 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccFastchop2counter_delay = - t_fast2a;
#line 59 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccFastchop2counter_isfirst = 0;
#line 59 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccFastchop2counter_n_pulse = 1;
#line 59 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccFastchop2counter_abs_out = 1;
#line 59 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccFastchop2counter_phase = 0;
#line 59 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccFastchop2counter_xwidth = 0;
#line 59 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccFastchop2counter_verbose = 0;
#line 12349 "./ESS_IN5_reprate.c"

  SIG_MESSAGE("Fastchop2counter (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 12356 "./ESS_IN5_reprate.c"
  rot_mul(mctr1, mcrotaOrigin, mcrotaFastchop2counter);
  rot_transpose(mcrotaFastchop2, mctr1);
  rot_mul(mcrotaFastchop2counter, mctr1, mcrotrFastchop2counter);
  mctc1 = coords_set(
#line 272 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    0,
#line 272 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    0.04,
#line 272 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    mcipLength + mcipGUI_GAP);
#line 12367 "./ESS_IN5_reprate.c"
  rot_transpose(mcrotaOrigin, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaFastchop2counter = coords_add(mcposaOrigin, mctc2);
  mctc1 = coords_sub(mcposaFastchop2, mcposaFastchop2counter);
  mcposrFastchop2counter = rot_apply(mcrotaFastchop2counter, mctc1);
  mcDEBUG_COMPONENT("Fastchop2counter", mcposaFastchop2counter, mcrotaFastchop2counter)
  mccomp_posa[32] = mcposaFastchop2counter;
  mccomp_posr[32] = mcposrFastchop2counter;
  mcNCounter[32]  = mcPCounter[32] = mcP2Counter[32] = 0;
  mcAbsorbProp[32]= 0;
    /* Component FOchop3. */
  /* Setting parameters for component FOchop3. */
  SIG_MESSAGE("FOchop3 (Init:SetPar)");
#line 275 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccFOchop3_theta_0 = 2 * mcipFAST_THETA;
#line 275 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccFOchop3_radius = 0.35;
#line 275 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccFOchop3_yheight = 0.35;
#line 275 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccFOchop3_nu = mcipF_fast2 / mcipFO3;
#line 275 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccFOchop3_nslit = mcipN_fast;
#line 59 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccFOchop3_jitter = 0;
#line 275 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccFOchop3_delay = t_fast3;
#line 59 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccFOchop3_isfirst = 0;
#line 59 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccFOchop3_n_pulse = 1;
#line 59 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccFOchop3_abs_out = 1;
#line 59 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccFOchop3_phase = 0;
#line 59 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccFOchop3_xwidth = 0;
#line 59 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccFOchop3_verbose = 0;
#line 12407 "./ESS_IN5_reprate.c"

  SIG_MESSAGE("FOchop3 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 12414 "./ESS_IN5_reprate.c"
  rot_mul(mctr1, mcrotaOrigin, mcrotaFOchop3);
  rot_transpose(mcrotaFastchop2counter, mctr1);
  rot_mul(mcrotaFOchop3, mctr1, mcrotrFOchop3);
  mctc1 = coords_set(
#line 276 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    0,
#line 276 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    0.04,
#line 276 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    mcipLength + 2 * mcipGUI_GAP);
#line 12425 "./ESS_IN5_reprate.c"
  rot_transpose(mcrotaOrigin, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaFOchop3 = coords_add(mcposaOrigin, mctc2);
  mctc1 = coords_sub(mcposaFastchop2counter, mcposaFOchop3);
  mcposrFOchop3 = rot_apply(mcrotaFOchop3, mctc1);
  mcDEBUG_COMPONENT("FOchop3", mcposaFOchop3, mcrotaFOchop3)
  mccomp_posa[33] = mcposaFOchop3;
  mccomp_posr[33] = mcposrFOchop3;
  mcNCounter[33]  = mcPCounter[33] = mcP2Counter[33] = 0;
  mcAbsorbProp[33]= 0;
    /* Component TOFfast2_zoom. */
  /* Setting parameters for component TOFfast2_zoom. */
  SIG_MESSAGE("TOFfast2_zoom (Init:SetPar)");
#line 279 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  if("TOF_fast2.dat") strncpy(mccTOFfast2_zoom_filename, "TOF_fast2.dat" ? "TOF_fast2.dat" : "", 16384); else mccTOFfast2_zoom_filename[0]='\0';
#line 47 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFfast2_zoom_xmin = -0.05;
#line 47 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFfast2_zoom_xmax = 0.05;
#line 47 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFfast2_zoom_ymin = -0.05;
#line 47 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFfast2_zoom_ymax = 0.05;
#line 279 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFfast2_zoom_xwidth = 1.1;
#line 280 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFfast2_zoom_yheight = 1.2;
#line 280 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFfast2_zoom_tmin = 1e6 * t_fast3 -2e2;
#line 280 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFfast2_zoom_tmax = 1e6 * t_fast3 + 2e2;
#line 48 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFfast2_zoom_dt = 1.0;
#line 48 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFfast2_zoom_restore_neutron = 0;
#line 48 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFfast2_zoom_nowritefile = 0;
#line 12463 "./ESS_IN5_reprate.c"

  SIG_MESSAGE("TOFfast2_zoom (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 12470 "./ESS_IN5_reprate.c"
  rot_mul(mctr1, mcrotaFOchop3, mcrotaTOFfast2_zoom);
  rot_transpose(mcrotaFOchop3, mctr1);
  rot_mul(mcrotaTOFfast2_zoom, mctr1, mcrotrTOFfast2_zoom);
  mctc1 = coords_set(
#line 281 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    0,
#line 281 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    -0.04,
#line 281 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    1e-6);
#line 12481 "./ESS_IN5_reprate.c"
  rot_transpose(mcrotaFOchop3, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaTOFfast2_zoom = coords_add(mcposaFOchop3, mctc2);
  mctc1 = coords_sub(mcposaFOchop3, mcposaTOFfast2_zoom);
  mcposrTOFfast2_zoom = rot_apply(mcrotaTOFfast2_zoom, mctc1);
  mcDEBUG_COMPONENT("TOFfast2_zoom", mcposaTOFfast2_zoom, mcrotaTOFfast2_zoom)
  mccomp_posa[34] = mcposaTOFfast2_zoom;
  mccomp_posr[34] = mcposrTOFfast2_zoom;
  mcNCounter[34]  = mcPCounter[34] = mcP2Counter[34] = 0;
  mcAbsorbProp[34]= 0;
    /* Component Lmon_afterfast2. */
  /* Setting parameters for component Lmon_afterfast2. */
  SIG_MESSAGE("Lmon_afterfast2 (Init:SetPar)");
#line 284 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  if("Lmon_afterfast2.dat") strncpy(mccLmon_afterfast2_filename, "Lmon_afterfast2.dat" ? "Lmon_afterfast2.dat" : "", 16384); else mccLmon_afterfast2_filename[0]='\0';
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmon_afterfast2_xmin = -0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmon_afterfast2_xmax = 0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmon_afterfast2_ymin = -0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmon_afterfast2_ymax = 0.05;
#line 284 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmon_afterfast2_xwidth = 0.06;
#line 285 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmon_afterfast2_yheight = 0.18;
#line 285 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmon_afterfast2_Lmin = 0;
#line 285 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmon_afterfast2_Lmax = mcipLmax + 1;
#line 51 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmon_afterfast2_restore_neutron = 0;
#line 51 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmon_afterfast2_nowritefile = 0;
#line 12517 "./ESS_IN5_reprate.c"

  SIG_MESSAGE("Lmon_afterfast2 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 12524 "./ESS_IN5_reprate.c"
  rot_mul(mctr1, mcrotaTOFfast2_zoom, mcrotaLmon_afterfast2);
  rot_transpose(mcrotaTOFfast2_zoom, mctr1);
  rot_mul(mcrotaLmon_afterfast2, mctr1, mcrotrLmon_afterfast2);
  mctc1 = coords_set(
#line 286 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    0,
#line 286 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    0,
#line 286 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    1e-6);
#line 12535 "./ESS_IN5_reprate.c"
  rot_transpose(mcrotaTOFfast2_zoom, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaLmon_afterfast2 = coords_add(mcposaTOFfast2_zoom, mctc2);
  mctc1 = coords_sub(mcposaTOFfast2_zoom, mcposaLmon_afterfast2);
  mcposrLmon_afterfast2 = rot_apply(mcrotaLmon_afterfast2, mctc1);
  mcDEBUG_COMPONENT("Lmon_afterfast2", mcposaLmon_afterfast2, mcrotaLmon_afterfast2)
  mccomp_posa[35] = mcposaLmon_afterfast2;
  mccomp_posr[35] = mcposrLmon_afterfast2;
  mcNCounter[35]  = mcPCounter[35] = mcP2Counter[35] = 0;
  mcAbsorbProp[35]= 0;
    /* Component TOFL_afterfast2. */
  /* Setting parameters for component TOFL_afterfast2. */
  SIG_MESSAGE("TOFL_afterfast2 (Init:SetPar)");
#line 289 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  if("TOF_afterfast2.dat") strncpy(mccTOFL_afterfast2_filename, "TOF_afterfast2.dat" ? "TOF_afterfast2.dat" : "", 16384); else mccTOFL_afterfast2_filename[0]='\0';
#line 49 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFL_afterfast2_xmin = -0.05;
#line 49 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFL_afterfast2_xmax = 0.05;
#line 49 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFL_afterfast2_ymin = -0.05;
#line 49 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFL_afterfast2_ymax = 0.05;
#line 290 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFL_afterfast2_xwidth = 0.05;
#line 290 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFL_afterfast2_yheight = 0.12;
#line 290 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFL_afterfast2_Lmin = 0;
#line 291 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFL_afterfast2_Lmax = mcipLmax + 1;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFL_afterfast2_restore_neutron = 0;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFL_afterfast2_nowritefile = 0;
#line 12571 "./ESS_IN5_reprate.c"

  SIG_MESSAGE("TOFL_afterfast2 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 12578 "./ESS_IN5_reprate.c"
  rot_mul(mctr1, mcrotaLmon_afterfast2, mcrotaTOFL_afterfast2);
  rot_transpose(mcrotaLmon_afterfast2, mctr1);
  rot_mul(mcrotaTOFL_afterfast2, mctr1, mcrotrTOFL_afterfast2);
  mctc1 = coords_set(
#line 292 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    0,
#line 292 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    0,
#line 292 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    1e-6);
#line 12589 "./ESS_IN5_reprate.c"
  rot_transpose(mcrotaLmon_afterfast2, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaTOFL_afterfast2 = coords_add(mcposaLmon_afterfast2, mctc2);
  mctc1 = coords_sub(mcposaLmon_afterfast2, mcposaTOFL_afterfast2);
  mcposrTOFL_afterfast2 = rot_apply(mcrotaTOFL_afterfast2, mctc1);
  mcDEBUG_COMPONENT("TOFL_afterfast2", mcposaTOFL_afterfast2, mcrotaTOFL_afterfast2)
  mccomp_posa[36] = mcposaTOFL_afterfast2;
  mccomp_posr[36] = mcposrTOFL_afterfast2;
  mcNCounter[36]  = mcPCounter[36] = mcP2Counter[36] = 0;
  mcAbsorbProp[36]= 0;
    /* Component TOFL_afterfast2_zoom. */
  /* Setting parameters for component TOFL_afterfast2_zoom. */
  SIG_MESSAGE("TOFL_afterfast2_zoom (Init:SetPar)");
#line 295 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  if("TOFL_afterfast2_zoom.dat") strncpy(mccTOFL_afterfast2_zoom_filename, "TOFL_afterfast2_zoom.dat" ? "TOFL_afterfast2_zoom.dat" : "", 16384); else mccTOFL_afterfast2_zoom_filename[0]='\0';
#line 49 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFL_afterfast2_zoom_xmin = -0.05;
#line 49 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFL_afterfast2_zoom_xmax = 0.05;
#line 49 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFL_afterfast2_zoom_ymin = -0.05;
#line 49 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFL_afterfast2_zoom_ymax = 0.05;
#line 296 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFL_afterfast2_zoom_xwidth = 0.05;
#line 296 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFL_afterfast2_zoom_yheight = 0.12;
#line 296 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFL_afterfast2_zoom_Lmin = mciplambda0 -0.2;
#line 297 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFL_afterfast2_zoom_Lmax = mciplambda0 + 0.2;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFL_afterfast2_zoom_restore_neutron = 0;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFL_afterfast2_zoom_nowritefile = 0;
#line 12625 "./ESS_IN5_reprate.c"

  SIG_MESSAGE("TOFL_afterfast2_zoom (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 12632 "./ESS_IN5_reprate.c"
  rot_mul(mctr1, mcrotaTOFL_afterfast2, mcrotaTOFL_afterfast2_zoom);
  rot_transpose(mcrotaTOFL_afterfast2, mctr1);
  rot_mul(mcrotaTOFL_afterfast2_zoom, mctr1, mcrotrTOFL_afterfast2_zoom);
  mctc1 = coords_set(
#line 298 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    0,
#line 298 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    0,
#line 298 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    1e-6);
#line 12643 "./ESS_IN5_reprate.c"
  rot_transpose(mcrotaTOFL_afterfast2, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaTOFL_afterfast2_zoom = coords_add(mcposaTOFL_afterfast2, mctc2);
  mctc1 = coords_sub(mcposaTOFL_afterfast2, mcposaTOFL_afterfast2_zoom);
  mcposrTOFL_afterfast2_zoom = rot_apply(mcrotaTOFL_afterfast2_zoom, mctc1);
  mcDEBUG_COMPONENT("TOFL_afterfast2_zoom", mcposaTOFL_afterfast2_zoom, mcrotaTOFL_afterfast2_zoom)
  mccomp_posa[37] = mcposaTOFL_afterfast2_zoom;
  mccomp_posr[37] = mcposrTOFL_afterfast2_zoom;
  mcNCounter[37]  = mcPCounter[37] = mcP2Counter[37] = 0;
  mcAbsorbProp[37]= 0;
    /* Component PSD_afterfast2. */
  /* Setting parameters for component PSD_afterfast2. */
  SIG_MESSAGE("PSD_afterfast2 (Init:SetPar)");
#line 49 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccPSD_afterfast2_nx = 90;
#line 49 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccPSD_afterfast2_ny = 90;
#line 301 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  if("PSD_afterfast2.dat") strncpy(mccPSD_afterfast2_filename, "PSD_afterfast2.dat" ? "PSD_afterfast2.dat" : "", 16384); else mccPSD_afterfast2_filename[0]='\0';
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccPSD_afterfast2_xmin = -0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccPSD_afterfast2_xmax = 0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccPSD_afterfast2_ymin = -0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccPSD_afterfast2_ymax = 0.05;
#line 301 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccPSD_afterfast2_xwidth = 0.1;
#line 301 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccPSD_afterfast2_yheight = 0.25;
#line 51 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccPSD_afterfast2_restore_neutron = 0;
#line 12677 "./ESS_IN5_reprate.c"

  SIG_MESSAGE("PSD_afterfast2 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 12684 "./ESS_IN5_reprate.c"
  rot_mul(mctr1, mcrotaTOFL_afterfast2_zoom, mcrotaPSD_afterfast2);
  rot_transpose(mcrotaTOFL_afterfast2_zoom, mctr1);
  rot_mul(mcrotaPSD_afterfast2, mctr1, mcrotrPSD_afterfast2);
  mctc1 = coords_set(
#line 302 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    0,
#line 302 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    0,
#line 302 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    1e-6);
#line 12695 "./ESS_IN5_reprate.c"
  rot_transpose(mcrotaTOFL_afterfast2_zoom, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaPSD_afterfast2 = coords_add(mcposaTOFL_afterfast2_zoom, mctc2);
  mctc1 = coords_sub(mcposaTOFL_afterfast2_zoom, mcposaPSD_afterfast2);
  mcposrPSD_afterfast2 = rot_apply(mcrotaPSD_afterfast2, mctc1);
  mcDEBUG_COMPONENT("PSD_afterfast2", mcposaPSD_afterfast2, mcrotaPSD_afterfast2)
  mccomp_posa[38] = mcposaPSD_afterfast2;
  mccomp_posr[38] = mcposrPSD_afterfast2;
  mcNCounter[38]  = mcPCounter[38] = mcP2Counter[38] = 0;
  mcAbsorbProp[38]= 0;
    /* Component Guidesample. */
  /* Setting parameters for component Guidesample. */
  SIG_MESSAGE("Guidesample (Init:SetPar)");
#line 58 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  if(0) strncpy(mccGuidesample_reflect, 0 ? 0 : "", 16384); else mccGuidesample_reflect[0]='\0';
#line 305 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccGuidesample_w1 = mcipW_chop;
#line 305 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccGuidesample_h1 = mcipH_chop;
#line 305 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccGuidesample_w2 = mcipW_end;
#line 305 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccGuidesample_h2 = mcipH_end;
#line 305 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccGuidesample_l = mcipSAMPLE_DIST -4 * mcipGUI_GAP;
#line 306 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccGuidesample_R0 = 1;
#line 58 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccGuidesample_Qc = 0.0219;
#line 306 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccGuidesample_alpha = mcipALPHA;
#line 306 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccGuidesample_m = mcipM;
#line 306 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccGuidesample_W = 0.003;
#line 12731 "./ESS_IN5_reprate.c"

  SIG_MESSAGE("Guidesample (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 12738 "./ESS_IN5_reprate.c"
  rot_mul(mctr1, mcrotaPSD_afterfast2, mcrotaGuidesample);
  rot_transpose(mcrotaPSD_afterfast2, mctr1);
  rot_mul(mcrotaGuidesample, mctr1, mcrotrGuidesample);
  mctc1 = coords_set(
#line 307 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    0,
#line 307 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    0,
#line 307 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    mcipGUI_GAP / 2);
#line 12749 "./ESS_IN5_reprate.c"
  rot_transpose(mcrotaPSD_afterfast2, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaGuidesample = coords_add(mcposaPSD_afterfast2, mctc2);
  mctc1 = coords_sub(mcposaPSD_afterfast2, mcposaGuidesample);
  mcposrGuidesample = rot_apply(mcrotaGuidesample, mctc1);
  mcDEBUG_COMPONENT("Guidesample", mcposaGuidesample, mcrotaGuidesample)
  mccomp_posa[39] = mcposaGuidesample;
  mccomp_posr[39] = mcposrGuidesample;
  mcNCounter[39]  = mcPCounter[39] = mcP2Counter[39] = 0;
  mcAbsorbProp[39]= 0;
    /* Component Lmon_guideend. */
  /* Setting parameters for component Lmon_guideend. */
  SIG_MESSAGE("Lmon_guideend (Init:SetPar)");
#line 310 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  if("Lmon_guideend.dat") strncpy(mccLmon_guideend_filename, "Lmon_guideend.dat" ? "Lmon_guideend.dat" : "", 16384); else mccLmon_guideend_filename[0]='\0';
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmon_guideend_xmin = -0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmon_guideend_xmax = 0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmon_guideend_ymin = -0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmon_guideend_ymax = 0.05;
#line 310 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmon_guideend_xwidth = mcipW_end + 0.01;
#line 311 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmon_guideend_yheight = mcipH_end + 0.01;
#line 311 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmon_guideend_Lmin = 0;
#line 311 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmon_guideend_Lmax = mcipLmax + 1;
#line 51 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmon_guideend_restore_neutron = 0;
#line 51 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmon_guideend_nowritefile = 0;
#line 12785 "./ESS_IN5_reprate.c"

  SIG_MESSAGE("Lmon_guideend (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 12792 "./ESS_IN5_reprate.c"
  rot_mul(mctr1, mcrotaGuidesample, mcrotaLmon_guideend);
  rot_transpose(mcrotaGuidesample, mctr1);
  rot_mul(mcrotaLmon_guideend, mctr1, mcrotrLmon_guideend);
  mctc1 = coords_set(
#line 312 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    0,
#line 312 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    0,
#line 312 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    mcipSAMPLE_DIST -4 * mcipGUI_GAP + 1e-6);
#line 12803 "./ESS_IN5_reprate.c"
  rot_transpose(mcrotaGuidesample, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaLmon_guideend = coords_add(mcposaGuidesample, mctc2);
  mctc1 = coords_sub(mcposaGuidesample, mcposaLmon_guideend);
  mcposrLmon_guideend = rot_apply(mcrotaLmon_guideend, mctc1);
  mcDEBUG_COMPONENT("Lmon_guideend", mcposaLmon_guideend, mcrotaLmon_guideend)
  mccomp_posa[40] = mcposaLmon_guideend;
  mccomp_posr[40] = mcposrLmon_guideend;
  mcNCounter[40]  = mcPCounter[40] = mcP2Counter[40] = 0;
  mcAbsorbProp[40]= 0;
    /* Component PSDsample. */
  /* Setting parameters for component PSDsample. */
  SIG_MESSAGE("PSDsample (Init:SetPar)");
#line 49 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccPSDsample_nx = 90;
#line 49 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccPSDsample_ny = 90;
#line 315 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  if("PSDsample.dat") strncpy(mccPSDsample_filename, "PSDsample.dat" ? "PSDsample.dat" : "", 16384); else mccPSDsample_filename[0]='\0';
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccPSDsample_xmin = -0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccPSDsample_xmax = 0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccPSDsample_ymin = -0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccPSDsample_ymax = 0.05;
#line 315 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccPSDsample_xwidth = 0.1;
#line 315 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccPSDsample_yheight = 0.25;
#line 51 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccPSDsample_restore_neutron = 0;
#line 12837 "./ESS_IN5_reprate.c"

  SIG_MESSAGE("PSDsample (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 12844 "./ESS_IN5_reprate.c"
  rot_mul(mctr1, mcrotaLmon_guideend, mcrotaPSDsample);
  rot_transpose(mcrotaLmon_guideend, mctr1);
  rot_mul(mcrotaPSDsample, mctr1, mcrotrPSDsample);
  mctc1 = coords_set(
#line 316 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    0,
#line 316 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    0,
#line 316 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    mcipSAMPLE_DIST -3 * mcipGUI_GAP + 1e-6);
#line 12855 "./ESS_IN5_reprate.c"
  rot_transpose(mcrotaLmon_guideend, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaPSDsample = coords_add(mcposaLmon_guideend, mctc2);
  mctc1 = coords_sub(mcposaLmon_guideend, mcposaPSDsample);
  mcposrPSDsample = rot_apply(mcrotaPSDsample, mctc1);
  mcDEBUG_COMPONENT("PSDsample", mcposaPSDsample, mcrotaPSDsample)
  mccomp_posa[41] = mcposaPSDsample;
  mccomp_posr[41] = mcposrPSDsample;
  mcNCounter[41]  = mcPCounter[41] = mcP2Counter[41] = 0;
  mcAbsorbProp[41]= 0;
    /* Component TOFsample_zoom. */
  /* Setting parameters for component TOFsample_zoom. */
  SIG_MESSAGE("TOFsample_zoom (Init:SetPar)");
#line 319 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  if("TOF_sample.dat") strncpy(mccTOFsample_zoom_filename, "TOF_sample.dat" ? "TOF_sample.dat" : "", 16384); else mccTOFsample_zoom_filename[0]='\0';
#line 47 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFsample_zoom_xmin = -0.05;
#line 47 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFsample_zoom_xmax = 0.05;
#line 47 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFsample_zoom_ymin = -0.05;
#line 47 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFsample_zoom_ymax = 0.05;
#line 319 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFsample_zoom_xwidth = 0.02;
#line 320 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFsample_zoom_yheight = 0.04;
#line 320 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFsample_zoom_tmin = 1e6 * t_sample -5e4;
#line 320 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFsample_zoom_tmax = 1e6 * t_sample + 5e4;
#line 48 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFsample_zoom_dt = 1.0;
#line 48 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFsample_zoom_restore_neutron = 0;
#line 48 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFsample_zoom_nowritefile = 0;
#line 12893 "./ESS_IN5_reprate.c"

  SIG_MESSAGE("TOFsample_zoom (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 12900 "./ESS_IN5_reprate.c"
  rot_mul(mctr1, mcrotaPSDsample, mcrotaTOFsample_zoom);
  rot_transpose(mcrotaPSDsample, mctr1);
  rot_mul(mcrotaTOFsample_zoom, mctr1, mcrotrTOFsample_zoom);
  mctc1 = coords_set(
#line 321 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    0,
#line 321 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    0,
#line 321 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    1e-6);
#line 12911 "./ESS_IN5_reprate.c"
  rot_transpose(mcrotaPSDsample, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaTOFsample_zoom = coords_add(mcposaPSDsample, mctc2);
  mctc1 = coords_sub(mcposaPSDsample, mcposaTOFsample_zoom);
  mcposrTOFsample_zoom = rot_apply(mcrotaTOFsample_zoom, mctc1);
  mcDEBUG_COMPONENT("TOFsample_zoom", mcposaTOFsample_zoom, mcrotaTOFsample_zoom)
  mccomp_posa[42] = mcposaTOFsample_zoom;
  mccomp_posr[42] = mcposrTOFsample_zoom;
  mcNCounter[42]  = mcPCounter[42] = mcP2Counter[42] = 0;
  mcAbsorbProp[42]= 0;
    /* Component Esample. */
  /* Setting parameters for component Esample. */
  SIG_MESSAGE("Esample (Init:SetPar)");
#line 324 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  if("Esample") strncpy(mccEsample_filename, "Esample" ? "Esample" : "", 16384); else mccEsample_filename[0]='\0';
#line 53 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccEsample_xmin = -0.05;
#line 53 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccEsample_xmax = 0.05;
#line 53 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccEsample_ymin = -0.05;
#line 53 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccEsample_ymax = 0.05;
#line 324 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccEsample_xwidth = 0.02;
#line 325 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccEsample_yheight = 0.04;
#line 325 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccEsample_Emin = E_target - mcipRES_DE;
#line 325 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccEsample_Emax = E_target + mcipRES_DE;
#line 54 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccEsample_restore_neutron = 0;
#line 54 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccEsample_nowritefile = 0;
#line 12947 "./ESS_IN5_reprate.c"

  SIG_MESSAGE("Esample (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 12954 "./ESS_IN5_reprate.c"
  rot_mul(mctr1, mcrotaTOFsample_zoom, mcrotaEsample);
  rot_transpose(mcrotaTOFsample_zoom, mctr1);
  rot_mul(mcrotaEsample, mctr1, mcrotrEsample);
  mctc1 = coords_set(
#line 326 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    0,
#line 326 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    0,
#line 326 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    1e-6);
#line 12965 "./ESS_IN5_reprate.c"
  rot_transpose(mcrotaTOFsample_zoom, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaEsample = coords_add(mcposaTOFsample_zoom, mctc2);
  mctc1 = coords_sub(mcposaTOFsample_zoom, mcposaEsample);
  mcposrEsample = rot_apply(mcrotaEsample, mctc1);
  mcDEBUG_COMPONENT("Esample", mcposaEsample, mcrotaEsample)
  mccomp_posa[43] = mcposaEsample;
  mccomp_posr[43] = mcposrEsample;
  mcNCounter[43]  = mcPCounter[43] = mcP2Counter[43] = 0;
  mcAbsorbProp[43]= 0;
    /* Component Lmon_sample_zoom. */
  /* Setting parameters for component Lmon_sample_zoom. */
  SIG_MESSAGE("Lmon_sample_zoom (Init:SetPar)");
#line 329 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  if("LMON_sample_zoom.dat") strncpy(mccLmon_sample_zoom_filename, "LMON_sample_zoom.dat" ? "LMON_sample_zoom.dat" : "", 16384); else mccLmon_sample_zoom_filename[0]='\0';
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmon_sample_zoom_xmin = -0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmon_sample_zoom_xmax = 0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmon_sample_zoom_ymin = -0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmon_sample_zoom_ymax = 0.05;
#line 329 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmon_sample_zoom_xwidth = 0.02;
#line 330 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmon_sample_zoom_yheight = 0.04;
#line 330 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmon_sample_zoom_Lmin = mciplambda0 -0.1;
#line 330 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmon_sample_zoom_Lmax = mciplambda0 + 0.1;
#line 51 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmon_sample_zoom_restore_neutron = 0;
#line 51 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccLmon_sample_zoom_nowritefile = 0;
#line 13001 "./ESS_IN5_reprate.c"

  SIG_MESSAGE("Lmon_sample_zoom (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 13008 "./ESS_IN5_reprate.c"
  rot_mul(mctr1, mcrotaEsample, mcrotaLmon_sample_zoom);
  rot_transpose(mcrotaEsample, mctr1);
  rot_mul(mcrotaLmon_sample_zoom, mctr1, mcrotrLmon_sample_zoom);
  mctc1 = coords_set(
#line 331 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    0,
#line 331 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    0,
#line 331 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    1e-6);
#line 13019 "./ESS_IN5_reprate.c"
  rot_transpose(mcrotaEsample, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaLmon_sample_zoom = coords_add(mcposaEsample, mctc2);
  mctc1 = coords_sub(mcposaEsample, mcposaLmon_sample_zoom);
  mcposrLmon_sample_zoom = rot_apply(mcrotaLmon_sample_zoom, mctc1);
  mcDEBUG_COMPONENT("Lmon_sample_zoom", mcposaLmon_sample_zoom, mcrotaLmon_sample_zoom)
  mccomp_posa[44] = mcposaLmon_sample_zoom;
  mccomp_posr[44] = mcposrLmon_sample_zoom;
  mcNCounter[44]  = mcPCounter[44] = mcP2Counter[44] = 0;
  mcAbsorbProp[44]= 0;
    /* Component sample. */
  /* Setting parameters for component sample. */
  SIG_MESSAGE("sample (Init:SetPar)");
#line 334 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccsample_thickness = 0.01 - mcipV_HOLE;
#line 334 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccsample_radius = 0.01;
#line 83 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccsample_focus_r = 0;
#line 334 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccsample_p_interact = 1;
#line 335 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccsample_f_QE = mcipFRAC_QUASIEL;
#line 335 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccsample_f_tun = mcipFRAC_TUNNEL;
#line 335 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccsample_gamma = mcipGamma;
#line 335 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccsample_E_tun = mcipEtun;
#line 85 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccsample_target_x = 0;
#line 85 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccsample_target_y = 0;
#line 85 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccsample_target_z = 0.235;
#line 336 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccsample_focus_xw = 0.015;
#line 336 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccsample_focus_yh = 0.015;
#line 86 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccsample_focus_aw = 0;
#line 86 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccsample_focus_ah = 0;
#line 86 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccsample_xwidth = 0;
#line 334 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccsample_yheight = 0.04;
#line 86 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccsample_zdepth = 0;
#line 86 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccsample_sigma_abs = 5.08;
#line 86 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccsample_sigma_inc = 4.935;
#line 86 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccsample_Vc = 13.827;
#line 336 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccsample_target_index = 2;
#line 13077 "./ESS_IN5_reprate.c"

  SIG_MESSAGE("sample (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 13084 "./ESS_IN5_reprate.c"
  rot_mul(mctr1, mcrotaFastchop2, mcrotasample);
  rot_transpose(mcrotaLmon_sample_zoom, mctr1);
  rot_mul(mcrotasample, mctr1, mcrotrsample);
  mctc1 = coords_set(
#line 337 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    0,
#line 337 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    -0.04,
#line 337 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    mcipSAMPLE_DIST);
#line 13095 "./ESS_IN5_reprate.c"
  rot_transpose(mcrotaFastchop2, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposasample = coords_add(mcposaFastchop2, mctc2);
  mctc1 = coords_sub(mcposaLmon_sample_zoom, mcposasample);
  mcposrsample = rot_apply(mcrotasample, mctc1);
  mcDEBUG_COMPONENT("sample", mcposasample, mcrotasample)
  mccomp_posa[45] = mcposasample;
  mccomp_posr[45] = mcposrsample;
  mcNCounter[45]  = mcPCounter[45] = mcP2Counter[45] = 0;
  mcAbsorbProp[45]= 0;
    /* Component detectorarm. */
  /* Setting parameters for component detectorarm. */
  SIG_MESSAGE("detectorarm (Init:SetPar)");

  SIG_MESSAGE("detectorarm (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 342 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    (0)*DEG2RAD,
#line 342 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    (mcipTT)*DEG2RAD,
#line 342 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    (0)*DEG2RAD);
#line 13118 "./ESS_IN5_reprate.c"
  rot_mul(mctr1, mcrotasample, mcrotadetectorarm);
  rot_transpose(mcrotasample, mctr1);
  rot_mul(mcrotadetectorarm, mctr1, mcrotrdetectorarm);
  mctc1 = coords_set(
#line 341 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    0,
#line 341 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    0,
#line 341 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    0);
#line 13129 "./ESS_IN5_reprate.c"
  rot_transpose(mcrotasample, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposadetectorarm = coords_add(mcposasample, mctc2);
  mctc1 = coords_sub(mcposasample, mcposadetectorarm);
  mcposrdetectorarm = rot_apply(mcrotadetectorarm, mctc1);
  mcDEBUG_COMPONENT("detectorarm", mcposadetectorarm, mcrotadetectorarm)
  mccomp_posa[46] = mcposadetectorarm;
  mccomp_posr[46] = mcposrdetectorarm;
  mcNCounter[46]  = mcPCounter[46] = mcP2Counter[46] = 0;
  mcAbsorbProp[46]= 0;
    /* Component TOFdetector. */
  /* Setting parameters for component TOFdetector. */
  SIG_MESSAGE("TOFdetector (Init:SetPar)");
#line 345 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  if("TOF.dat") strncpy(mccTOFdetector_filename, "TOF.dat" ? "TOF.dat" : "", 16384); else mccTOFdetector_filename[0]='\0';
#line 47 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFdetector_xmin = -0.05;
#line 47 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFdetector_xmax = 0.05;
#line 47 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFdetector_ymin = -0.05;
#line 47 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFdetector_ymax = 0.05;
#line 345 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFdetector_xwidth = 0.015;
#line 346 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFdetector_yheight = 0.015;
#line 346 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFdetector_tmin = 0;
#line 346 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFdetector_tmax = 2e5;
#line 48 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFdetector_dt = 1.0;
#line 48 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFdetector_restore_neutron = 0;
#line 48 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFdetector_nowritefile = 0;
#line 13167 "./ESS_IN5_reprate.c"

  SIG_MESSAGE("TOFdetector (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 13174 "./ESS_IN5_reprate.c"
  rot_mul(mctr1, mcrotadetectorarm, mcrotaTOFdetector);
  rot_transpose(mcrotadetectorarm, mctr1);
  rot_mul(mcrotaTOFdetector, mctr1, mcrotrTOFdetector);
  mctc1 = coords_set(
#line 347 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    0,
#line 347 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    0,
#line 347 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    mcipDETECTOR_DIST);
#line 13185 "./ESS_IN5_reprate.c"
  rot_transpose(mcrotadetectorarm, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaTOFdetector = coords_add(mcposadetectorarm, mctc2);
  mctc1 = coords_sub(mcposadetectorarm, mcposaTOFdetector);
  mcposrTOFdetector = rot_apply(mcrotaTOFdetector, mctc1);
  mcDEBUG_COMPONENT("TOFdetector", mcposaTOFdetector, mcrotaTOFdetector)
  mccomp_posa[47] = mcposaTOFdetector;
  mccomp_posr[47] = mcposrTOFdetector;
  mcNCounter[47]  = mcPCounter[47] = mcP2Counter[47] = 0;
  mcAbsorbProp[47]= 0;
    /* Component TOFdetector_zoom. */
  /* Setting parameters for component TOFdetector_zoom. */
  SIG_MESSAGE("TOFdetector_zoom (Init:SetPar)");
#line 350 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  if("TOF_zoom.dat") strncpy(mccTOFdetector_zoom_filename, "TOF_zoom.dat" ? "TOF_zoom.dat" : "", 16384); else mccTOFdetector_zoom_filename[0]='\0';
#line 47 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFdetector_zoom_xmin = -0.05;
#line 47 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFdetector_zoom_xmax = 0.05;
#line 47 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFdetector_zoom_ymin = -0.05;
#line 47 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFdetector_zoom_ymax = 0.05;
#line 350 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFdetector_zoom_xwidth = 0.015;
#line 351 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFdetector_zoom_yheight = 0.015;
#line 351 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFdetector_zoom_tmin = 1e6 * t_detector -10e2;
#line 351 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFdetector_zoom_tmax = 1e6 * t_detector + 10e2;
#line 48 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFdetector_zoom_dt = 1.0;
#line 48 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFdetector_zoom_restore_neutron = 0;
#line 48 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOFdetector_zoom_nowritefile = 0;
#line 13223 "./ESS_IN5_reprate.c"

  SIG_MESSAGE("TOFdetector_zoom (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 13230 "./ESS_IN5_reprate.c"
  rot_mul(mctr1, mcrotaTOFdetector, mcrotaTOFdetector_zoom);
  rot_transpose(mcrotaTOFdetector, mctr1);
  rot_mul(mcrotaTOFdetector_zoom, mctr1, mcrotrTOFdetector_zoom);
  mctc1 = coords_set(
#line 352 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    0,
#line 352 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    0,
#line 352 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    1e-6);
#line 13241 "./ESS_IN5_reprate.c"
  rot_transpose(mcrotaTOFdetector, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaTOFdetector_zoom = coords_add(mcposaTOFdetector, mctc2);
  mctc1 = coords_sub(mcposaTOFdetector, mcposaTOFdetector_zoom);
  mcposrTOFdetector_zoom = rot_apply(mcrotaTOFdetector_zoom, mctc1);
  mcDEBUG_COMPONENT("TOFdetector_zoom", mcposaTOFdetector_zoom, mcrotaTOFdetector_zoom)
  mccomp_posa[48] = mcposaTOFdetector_zoom;
  mccomp_posr[48] = mcposrTOFdetector_zoom;
  mcNCounter[48]  = mcPCounter[48] = mcP2Counter[48] = 0;
  mcAbsorbProp[48]= 0;
    /* Component Edetector. */
  /* Setting parameters for component Edetector. */
  SIG_MESSAGE("Edetector (Init:SetPar)");
#line 355 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  if("Edet") strncpy(mccEdetector_filename, "Edet" ? "Edet" : "", 16384); else mccEdetector_filename[0]='\0';
#line 53 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccEdetector_xmin = -0.05;
#line 53 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccEdetector_xmax = 0.05;
#line 53 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccEdetector_ymin = -0.05;
#line 53 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccEdetector_ymax = 0.05;
#line 355 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccEdetector_xwidth = 0.015;
#line 356 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccEdetector_yheight = 0.015;
#line 356 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccEdetector_Emin = E_target - mcipRES_DE;
#line 356 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccEdetector_Emax = E_target + mcipRES_DE;
#line 54 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccEdetector_restore_neutron = 0;
#line 54 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccEdetector_nowritefile = 0;
#line 13277 "./ESS_IN5_reprate.c"

  SIG_MESSAGE("Edetector (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 13284 "./ESS_IN5_reprate.c"
  rot_mul(mctr1, mcrotaTOFdetector_zoom, mcrotaEdetector);
  rot_transpose(mcrotaTOFdetector_zoom, mctr1);
  rot_mul(mcrotaEdetector, mctr1, mcrotrEdetector);
  mctc1 = coords_set(
#line 357 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    0,
#line 357 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    0,
#line 357 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    1e-6);
#line 13295 "./ESS_IN5_reprate.c"
  rot_transpose(mcrotaTOFdetector_zoom, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaEdetector = coords_add(mcposaTOFdetector_zoom, mctc2);
  mctc1 = coords_sub(mcposaTOFdetector_zoom, mcposaEdetector);
  mcposrEdetector = rot_apply(mcrotaEdetector, mctc1);
  mcDEBUG_COMPONENT("Edetector", mcposaEdetector, mcrotaEdetector)
  mccomp_posa[49] = mcposaEdetector;
  mccomp_posr[49] = mcposrEdetector;
  mcNCounter[49]  = mcPCounter[49] = mcP2Counter[49] = 0;
  mcAbsorbProp[49]= 0;
    /* Component TOF2Edetector. */
  /* Setting parameters for component TOF2Edetector. */
  SIG_MESSAGE("TOF2Edetector (Init:SetPar)");
#line 360 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  if("TOF2E.dat") strncpy(mccTOF2Edetector_filename, "TOF2E.dat" ? "TOF2E.dat" : "", 16384); else mccTOF2Edetector_filename[0]='\0';
#line 55 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOF2Edetector_xmin = -0.05;
#line 55 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOF2Edetector_xmax = 0.05;
#line 55 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOF2Edetector_ymin = -0.05;
#line 55 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOF2Edetector_ymax = 0.05;
#line 360 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOF2Edetector_xwidth = 0.015;
#line 361 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOF2Edetector_yheight = 0.015;
#line 361 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOF2Edetector_Emin = E_target - mcipRES_DE;
#line 361 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOF2Edetector_Emax = E_target + mcipRES_DE;
#line 361 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOF2Edetector_T_zero = t_sample;
#line 361 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOF2Edetector_L_flight = mcipDETECTOR_DIST;
#line 55 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOF2Edetector_restore_neutron = 0;
#line 55 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
  mccTOF2Edetector_nowritefile = 0;
#line 13335 "./ESS_IN5_reprate.c"

  SIG_MESSAGE("TOF2Edetector (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 13342 "./ESS_IN5_reprate.c"
  rot_mul(mctr1, mcrotaEdetector, mcrotaTOF2Edetector);
  rot_transpose(mcrotaEdetector, mctr1);
  rot_mul(mcrotaTOF2Edetector, mctr1, mcrotrTOF2Edetector);
  mctc1 = coords_set(
#line 362 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    0,
#line 362 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    0,
#line 362 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ESS_IN5_reprate/ESS_IN5_reprate.instr"
    1e-6);
#line 13353 "./ESS_IN5_reprate.c"
  rot_transpose(mcrotaEdetector, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaTOF2Edetector = coords_add(mcposaEdetector, mctc2);
  mctc1 = coords_sub(mcposaEdetector, mcposaTOF2Edetector);
  mcposrTOF2Edetector = rot_apply(mcrotaTOF2Edetector, mctc1);
  mcDEBUG_COMPONENT("TOF2Edetector", mcposaTOF2Edetector, mcrotaTOF2Edetector)
  mccomp_posa[50] = mcposaTOF2Edetector;
  mccomp_posr[50] = mcposrTOF2Edetector;
  mcNCounter[50]  = mcPCounter[50] = mcP2Counter[50] = 0;
  mcAbsorbProp[50]= 0;
  /* Component initializations. */
  /* Initializations for component source. */
  SIG_MESSAGE("source (Init)");
#define mccompcurname  source
#define mccompcurtype  ESS_moderator_long
#define mccompcurindex 1
#define l_range mccsource_l_range
#define w_mult mccsource_w_mult
#define w_geom mccsource_w_geom
#define w_geom_c mccsource_w_geom_c
#define w_geom_t mccsource_w_geom_t
#define tx mccsource_tx
#define ty mccsource_ty
#define tz mccsource_tz
#define t1x mccsource_t1x
#define t1y mccsource_t1y
#define t1z mccsource_t1z
#define t2x mccsource_t2x
#define t2y mccsource_t2y
#define t2z mccsource_t2z
#define T_n mccsource_T_n
#define tau_n mccsource_tau_n
#define tau1_n mccsource_tau1_n
#define tau2_n mccsource_tau2_n
#define chi2_n mccsource_chi2_n
#define I0_n mccsource_I0_n
#define I2_n mccsource_I2_n
#define branch1_n mccsource_branch1_n
#define branch2_n mccsource_branch2_n
#define r_empty mccsource_r_empty
#define r_optics mccsource_r_optics
#define width_c mccsource_width_c
#define yheight mccsource_yheight
#define Lmin mccsource_Lmin
#define Lmax mccsource_Lmax
#define dist mccsource_dist
#define focus_xw mccsource_focus_xw
#define focus_yh mccsource_focus_yh
#define nu mccsource_nu
#define T mccsource_T
#define tau mccsource_tau
#define tau1 mccsource_tau1
#define tau2 mccsource_tau2
#define d mccsource_d
#define n mccsource_n
#define cold_frac mccsource_cold_frac
#define n2 mccsource_n2
#define chi2 mccsource_chi2
#define I0 mccsource_I0
#define I2 mccsource_I2
#define target_index mccsource_target_index
#define cyl_radius mccsource_cyl_radius
#define branch1 mccsource_branch1
#define branch2 mccsource_branch2
#define branch_tail mccsource_branch_tail
#define n_pulses mccsource_n_pulses
#define width_t mccsource_width_t
#define T_t mccsource_T_t
#define tau_t mccsource_tau_t
#define tau1_t mccsource_tau1_t
#define tau2_t mccsource_tau2_t
#define chi2_t mccsource_chi2_t
#define I0_t mccsource_I0_t
#define I2_t mccsource_I2_t
#define branch1_t mccsource_branch1_t
#define branch2_t mccsource_branch2_t
#define src_2012 mccsource_src_2012
#define tfocus_dist mccsource_tfocus_dist
#define tfocus_time mccsource_tfocus_time
#define tfocus_width mccsource_tfocus_width
#define beamport_angle mccsource_beamport_angle
#line 166 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../obsolete/ESS_moderator_long.comp"
{
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

}
#line 13521 "./ESS_IN5_reprate.c"
#undef beamport_angle
#undef tfocus_width
#undef tfocus_time
#undef tfocus_dist
#undef src_2012
#undef branch2_t
#undef branch1_t
#undef I2_t
#undef I0_t
#undef chi2_t
#undef tau2_t
#undef tau1_t
#undef tau_t
#undef T_t
#undef width_t
#undef n_pulses
#undef branch_tail
#undef branch2
#undef branch1
#undef cyl_radius
#undef target_index
#undef I2
#undef I0
#undef chi2
#undef n2
#undef cold_frac
#undef n
#undef d
#undef tau2
#undef tau1
#undef tau
#undef T
#undef nu
#undef focus_yh
#undef focus_xw
#undef dist
#undef Lmax
#undef Lmin
#undef yheight
#undef width_c
#undef r_optics
#undef r_empty
#undef branch2_n
#undef branch1_n
#undef I2_n
#undef I0_n
#undef chi2_n
#undef tau2_n
#undef tau1_n
#undef tau_n
#undef T_n
#undef t2z
#undef t2y
#undef t2x
#undef t1z
#undef t1y
#undef t1x
#undef tz
#undef ty
#undef tx
#undef w_geom_t
#undef w_geom_c
#undef w_geom
#undef w_mult
#undef l_range
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component Origin. */
  SIG_MESSAGE("Origin (Init)");
#define mccompcurname  Origin
#define mccompcurtype  Progress_bar
#define mccompcurindex 2
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
#line 13616 "./ESS_IN5_reprate.c"
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

  /* Initializations for component TOFmoderator_zoom. */
  SIG_MESSAGE("TOFmoderator_zoom (Init)");
#define mccompcurname  TOFmoderator_zoom
#define mccompcurtype  TOF_monitor
#define mccompcurindex 3
#define nt mccTOFmoderator_zoom_nt
#define TOF_N mccTOFmoderator_zoom_TOF_N
#define TOF_p mccTOFmoderator_zoom_TOF_p
#define TOF_p2 mccTOFmoderator_zoom_TOF_p2
#define t_min mccTOFmoderator_zoom_t_min
#define t_max mccTOFmoderator_zoom_t_max
#define delta_t mccTOFmoderator_zoom_delta_t
#define filename mccTOFmoderator_zoom_filename
#define xmin mccTOFmoderator_zoom_xmin
#define xmax mccTOFmoderator_zoom_xmax
#define ymin mccTOFmoderator_zoom_ymin
#define ymax mccTOFmoderator_zoom_ymax
#define xwidth mccTOFmoderator_zoom_xwidth
#define yheight mccTOFmoderator_zoom_yheight
#define tmin mccTOFmoderator_zoom_tmin
#define tmax mccTOFmoderator_zoom_tmax
#define dt mccTOFmoderator_zoom_dt
#define restore_neutron mccTOFmoderator_zoom_restore_neutron
#define nowritefile mccTOFmoderator_zoom_nowritefile
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
#line 13686 "./ESS_IN5_reprate.c"
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

  /* Initializations for component TOFmoderator. */
  SIG_MESSAGE("TOFmoderator (Init)");
#define mccompcurname  TOFmoderator
#define mccompcurtype  TOF_monitor
#define mccompcurindex 4
#define nt mccTOFmoderator_nt
#define TOF_N mccTOFmoderator_TOF_N
#define TOF_p mccTOFmoderator_TOF_p
#define TOF_p2 mccTOFmoderator_TOF_p2
#define t_min mccTOFmoderator_t_min
#define t_max mccTOFmoderator_t_max
#define delta_t mccTOFmoderator_delta_t
#define filename mccTOFmoderator_filename
#define xmin mccTOFmoderator_xmin
#define xmax mccTOFmoderator_xmax
#define ymin mccTOFmoderator_ymin
#define ymax mccTOFmoderator_ymax
#define xwidth mccTOFmoderator_xwidth
#define yheight mccTOFmoderator_yheight
#define tmin mccTOFmoderator_tmin
#define tmax mccTOFmoderator_tmax
#define dt mccTOFmoderator_dt
#define restore_neutron mccTOFmoderator_restore_neutron
#define nowritefile mccTOFmoderator_nowritefile
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
#line 13767 "./ESS_IN5_reprate.c"
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

  /* Initializations for component Lmon_guistart. */
  SIG_MESSAGE("Lmon_guistart (Init)");
#define mccompcurname  Lmon_guistart
#define mccompcurtype  L_monitor
#define mccompcurindex 5
#define nL mccLmon_guistart_nL
#define L_N mccLmon_guistart_L_N
#define L_p mccLmon_guistart_L_p
#define L_p2 mccLmon_guistart_L_p2
#define filename mccLmon_guistart_filename
#define xmin mccLmon_guistart_xmin
#define xmax mccLmon_guistart_xmax
#define ymin mccLmon_guistart_ymin
#define ymax mccLmon_guistart_ymax
#define xwidth mccLmon_guistart_xwidth
#define yheight mccLmon_guistart_yheight
#define Lmin mccLmon_guistart_Lmin
#define Lmax mccLmon_guistart_Lmax
#define restore_neutron mccLmon_guistart_restore_neutron
#define nowritefile mccLmon_guistart_nowritefile
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
#line 13832 "./ESS_IN5_reprate.c"
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

  /* Initializations for component Lmon_normalize. */
  SIG_MESSAGE("Lmon_normalize (Init)");
#define mccompcurname  Lmon_normalize
#define mccompcurtype  L_monitor
#define mccompcurindex 6
#define nL mccLmon_normalize_nL
#define L_N mccLmon_normalize_L_N
#define L_p mccLmon_normalize_L_p
#define L_p2 mccLmon_normalize_L_p2
#define filename mccLmon_normalize_filename
#define xmin mccLmon_normalize_xmin
#define xmax mccLmon_normalize_xmax
#define ymin mccLmon_normalize_ymin
#define ymax mccLmon_normalize_ymax
#define xwidth mccLmon_normalize_xwidth
#define yheight mccLmon_normalize_yheight
#define Lmin mccLmon_normalize_Lmin
#define Lmax mccLmon_normalize_Lmax
#define restore_neutron mccLmon_normalize_restore_neutron
#define nowritefile mccLmon_normalize_nowritefile
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
#line 13893 "./ESS_IN5_reprate.c"
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

  /* Initializations for component Guide1. */
  SIG_MESSAGE("Guide1 (Init)");
#define mccompcurname  Guide1
#define mccompcurtype  Guide
#define mccompcurindex 7
#define pTable mccGuide1_pTable
#define reflect mccGuide1_reflect
#define w1 mccGuide1_w1
#define h1 mccGuide1_h1
#define w2 mccGuide1_w2
#define h2 mccGuide1_h2
#define l mccGuide1_l
#define R0 mccGuide1_R0
#define Qc mccGuide1_Qc
#define alpha mccGuide1_alpha
#define m mccGuide1_m
#define W mccGuide1_W
#line 73 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide.comp"
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
#line 13949 "./ESS_IN5_reprate.c"
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

  /* Initializations for component Lmonslow1. */
  SIG_MESSAGE("Lmonslow1 (Init)");
#define mccompcurname  Lmonslow1
#define mccompcurtype  L_monitor
#define mccompcurindex 8
#define nL mccLmonslow1_nL
#define L_N mccLmonslow1_L_N
#define L_p mccLmonslow1_L_p
#define L_p2 mccLmonslow1_L_p2
#define filename mccLmonslow1_filename
#define xmin mccLmonslow1_xmin
#define xmax mccLmonslow1_xmax
#define ymin mccLmonslow1_ymin
#define ymax mccLmonslow1_ymax
#define xwidth mccLmonslow1_xwidth
#define yheight mccLmonslow1_yheight
#define Lmin mccLmonslow1_Lmin
#define Lmax mccLmonslow1_Lmax
#define restore_neutron mccLmonslow1_restore_neutron
#define nowritefile mccLmonslow1_nowritefile
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
#line 14007 "./ESS_IN5_reprate.c"
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

  /* Initializations for component PSDslow1. */
  SIG_MESSAGE("PSDslow1 (Init)");
#define mccompcurname  PSDslow1
#define mccompcurtype  PSD_monitor
#define mccompcurindex 9
#define PSD_N mccPSDslow1_PSD_N
#define PSD_p mccPSDslow1_PSD_p
#define PSD_p2 mccPSDslow1_PSD_p2
#define nx mccPSDslow1_nx
#define ny mccPSDslow1_ny
#define filename mccPSDslow1_filename
#define xmin mccPSDslow1_xmin
#define xmax mccPSDslow1_xmax
#define ymin mccPSDslow1_ymin
#define ymax mccPSDslow1_ymax
#define xwidth mccPSDslow1_xwidth
#define yheight mccPSDslow1_yheight
#define restore_neutron mccPSDslow1_restore_neutron
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
#line 14070 "./ESS_IN5_reprate.c"
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

  /* Initializations for component FOchop1. */
  SIG_MESSAGE("FOchop1 (Init)");
#define mccompcurname  FOchop1
#define mccompcurtype  DiskChopper
#define mccompcurindex 10
#define Tg mccFOchop1_Tg
#define To mccFOchop1_To
#define delta_y mccFOchop1_delta_y
#define height mccFOchop1_height
#define omega mccFOchop1_omega
#define theta_0 mccFOchop1_theta_0
#define radius mccFOchop1_radius
#define yheight mccFOchop1_yheight
#define nu mccFOchop1_nu
#define nslit mccFOchop1_nslit
#define jitter mccFOchop1_jitter
#define delay mccFOchop1_delay
#define isfirst mccFOchop1_isfirst
#define n_pulse mccFOchop1_n_pulse
#define abs_out mccFOchop1_abs_out
#define phase mccFOchop1_phase
#define xwidth mccFOchop1_xwidth
#define verbose mccFOchop1_verbose
#line 70 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/DiskChopper.comp"
{
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
}
#line 14173 "./ESS_IN5_reprate.c"
#undef verbose
#undef xwidth
#undef phase
#undef abs_out
#undef n_pulse
#undef isfirst
#undef delay
#undef jitter
#undef nslit
#undef nu
#undef yheight
#undef radius
#undef theta_0
#undef omega
#undef height
#undef delta_y
#undef To
#undef Tg
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component TOFLmon1. */
  SIG_MESSAGE("TOFLmon1 (Init)");
#define mccompcurname  TOFLmon1
#define mccompcurtype  TOFLambda_monitor
#define mccompcurindex 11
#define nL mccTOFLmon1_nL
#define nt mccTOFLmon1_nt
#define tmin mccTOFLmon1_tmin
#define tmax mccTOFLmon1_tmax
#define tt_0 mccTOFLmon1_tt_0
#define tt_1 mccTOFLmon1_tt_1
#define TOFL_N mccTOFLmon1_TOFL_N
#define TOFL_p mccTOFLmon1_TOFL_p
#define TOFL_p2 mccTOFLmon1_TOFL_p2
#define filename mccTOFLmon1_filename
#define xmin mccTOFLmon1_xmin
#define xmax mccTOFLmon1_xmax
#define ymin mccTOFLmon1_ymin
#define ymax mccTOFLmon1_ymax
#define xwidth mccTOFLmon1_xwidth
#define yheight mccTOFLmon1_yheight
#define Lmin mccTOFLmon1_Lmin
#define Lmax mccTOFLmon1_Lmax
#define restore_neutron mccTOFLmon1_restore_neutron
#define nowritefile mccTOFLmon1_nowritefile
#line 61 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/TOFLambda_monitor.comp"
{
int i,j;

if (xwidth  > 0) { xmax = xwidth/2;  xmin = -xmax; }
    if (yheight > 0) { ymax = yheight/2; ymin = -ymax; }

    if ((xmin >= xmax) || (ymin >= ymax)) {
            printf("TOFlambda_monitor: %s: Null detection area !\n"
                   "ERROR              (xwidth,yheight,xmin,xmax,ymin,ymax). Exiting",
           NAME_CURRENT_COMP);
      exit(0);
    }

    tt_0 = tmin*1e-6;
    tt_1 = tmax*1e-6;
    for (i=0; i<nL; i++)
     for (j=0; j<nt; j++)
     {
      TOFL_N[j][i] = 0;
      TOFL_p[j][i] = 0;
      TOFL_p2[j][i] = 0;
     }
}
#line 14245 "./ESS_IN5_reprate.c"
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
#undef TOFL_p2
#undef TOFL_p
#undef TOFL_N
#undef tt_1
#undef tt_0
#undef tmax
#undef tmin
#undef nt
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component Lmon_afterslow1. */
  SIG_MESSAGE("Lmon_afterslow1 (Init)");
#define mccompcurname  Lmon_afterslow1
#define mccompcurtype  L_monitor
#define mccompcurindex 12
#define nL mccLmon_afterslow1_nL
#define L_N mccLmon_afterslow1_L_N
#define L_p mccLmon_afterslow1_L_p
#define L_p2 mccLmon_afterslow1_L_p2
#define filename mccLmon_afterslow1_filename
#define xmin mccLmon_afterslow1_xmin
#define xmax mccLmon_afterslow1_xmax
#define ymin mccLmon_afterslow1_ymin
#define ymax mccLmon_afterslow1_ymax
#define xwidth mccLmon_afterslow1_xwidth
#define yheight mccLmon_afterslow1_yheight
#define Lmin mccLmon_afterslow1_Lmin
#define Lmax mccLmon_afterslow1_Lmax
#define restore_neutron mccLmon_afterslow1_restore_neutron
#define nowritefile mccLmon_afterslow1_nowritefile
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
#line 14311 "./ESS_IN5_reprate.c"
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

  /* Initializations for component PSD_afterslow1. */
  SIG_MESSAGE("PSD_afterslow1 (Init)");
#define mccompcurname  PSD_afterslow1
#define mccompcurtype  PSD_monitor
#define mccompcurindex 13
#define PSD_N mccPSD_afterslow1_PSD_N
#define PSD_p mccPSD_afterslow1_PSD_p
#define PSD_p2 mccPSD_afterslow1_PSD_p2
#define nx mccPSD_afterslow1_nx
#define ny mccPSD_afterslow1_ny
#define filename mccPSD_afterslow1_filename
#define xmin mccPSD_afterslow1_xmin
#define xmax mccPSD_afterslow1_xmax
#define ymin mccPSD_afterslow1_ymin
#define ymax mccPSD_afterslow1_ymax
#define xwidth mccPSD_afterslow1_xwidth
#define yheight mccPSD_afterslow1_yheight
#define restore_neutron mccPSD_afterslow1_restore_neutron
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
#line 14374 "./ESS_IN5_reprate.c"
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

  /* Initializations for component Guidelong1. */
  SIG_MESSAGE("Guidelong1 (Init)");
#define mccompcurname  Guidelong1
#define mccompcurtype  Guide
#define mccompcurindex 14
#define pTable mccGuidelong1_pTable
#define reflect mccGuidelong1_reflect
#define w1 mccGuidelong1_w1
#define h1 mccGuidelong1_h1
#define w2 mccGuidelong1_w2
#define h2 mccGuidelong1_h2
#define l mccGuidelong1_l
#define R0 mccGuidelong1_R0
#define Qc mccGuidelong1_Qc
#define alpha mccGuidelong1_alpha
#define m mccGuidelong1_m
#define W mccGuidelong1_W
#line 73 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide.comp"
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
#line 14428 "./ESS_IN5_reprate.c"
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

  /* Initializations for component Guidelong1b. */
  SIG_MESSAGE("Guidelong1b (Init)");
#define mccompcurname  Guidelong1b
#define mccompcurtype  Guide
#define mccompcurindex 15
#define pTable mccGuidelong1b_pTable
#define reflect mccGuidelong1b_reflect
#define w1 mccGuidelong1b_w1
#define h1 mccGuidelong1b_h1
#define w2 mccGuidelong1b_w2
#define h2 mccGuidelong1b_h2
#define l mccGuidelong1b_l
#define R0 mccGuidelong1b_R0
#define Qc mccGuidelong1b_Qc
#define alpha mccGuidelong1b_alpha
#define m mccGuidelong1b_m
#define W mccGuidelong1b_W
#line 73 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide.comp"
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
#line 14481 "./ESS_IN5_reprate.c"
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

  /* Initializations for component Lmon_slow2. */
  SIG_MESSAGE("Lmon_slow2 (Init)");
#define mccompcurname  Lmon_slow2
#define mccompcurtype  L_monitor
#define mccompcurindex 16
#define nL mccLmon_slow2_nL
#define L_N mccLmon_slow2_L_N
#define L_p mccLmon_slow2_L_p
#define L_p2 mccLmon_slow2_L_p2
#define filename mccLmon_slow2_filename
#define xmin mccLmon_slow2_xmin
#define xmax mccLmon_slow2_xmax
#define ymin mccLmon_slow2_ymin
#define ymax mccLmon_slow2_ymax
#define xwidth mccLmon_slow2_xwidth
#define yheight mccLmon_slow2_yheight
#define Lmin mccLmon_slow2_Lmin
#define Lmax mccLmon_slow2_Lmax
#define restore_neutron mccLmon_slow2_restore_neutron
#define nowritefile mccLmon_slow2_nowritefile
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
#line 14539 "./ESS_IN5_reprate.c"
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

  /* Initializations for component FOchop2. */
  SIG_MESSAGE("FOchop2 (Init)");
#define mccompcurname  FOchop2
#define mccompcurtype  DiskChopper
#define mccompcurindex 17
#define Tg mccFOchop2_Tg
#define To mccFOchop2_To
#define delta_y mccFOchop2_delta_y
#define height mccFOchop2_height
#define omega mccFOchop2_omega
#define theta_0 mccFOchop2_theta_0
#define radius mccFOchop2_radius
#define yheight mccFOchop2_yheight
#define nu mccFOchop2_nu
#define nslit mccFOchop2_nslit
#define jitter mccFOchop2_jitter
#define delay mccFOchop2_delay
#define isfirst mccFOchop2_isfirst
#define n_pulse mccFOchop2_n_pulse
#define abs_out mccFOchop2_abs_out
#define phase mccFOchop2_phase
#define xwidth mccFOchop2_xwidth
#define verbose mccFOchop2_verbose
#line 70 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/DiskChopper.comp"
{
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
}
#line 14644 "./ESS_IN5_reprate.c"
#undef verbose
#undef xwidth
#undef phase
#undef abs_out
#undef n_pulse
#undef isfirst
#undef delay
#undef jitter
#undef nslit
#undef nu
#undef yheight
#undef radius
#undef theta_0
#undef omega
#undef height
#undef delta_y
#undef To
#undef Tg
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component Fastchop1. */
  SIG_MESSAGE("Fastchop1 (Init)");
#define mccompcurname  Fastchop1
#define mccompcurtype  DiskChopper
#define mccompcurindex 18
#define Tg mccFastchop1_Tg
#define To mccFastchop1_To
#define delta_y mccFastchop1_delta_y
#define height mccFastchop1_height
#define omega mccFastchop1_omega
#define theta_0 mccFastchop1_theta_0
#define radius mccFastchop1_radius
#define yheight mccFastchop1_yheight
#define nu mccFastchop1_nu
#define nslit mccFastchop1_nslit
#define jitter mccFastchop1_jitter
#define delay mccFastchop1_delay
#define isfirst mccFastchop1_isfirst
#define n_pulse mccFastchop1_n_pulse
#define abs_out mccFastchop1_abs_out
#define phase mccFastchop1_phase
#define xwidth mccFastchop1_xwidth
#define verbose mccFastchop1_verbose
#line 70 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/DiskChopper.comp"
{
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
}
#line 14752 "./ESS_IN5_reprate.c"
#undef verbose
#undef xwidth
#undef phase
#undef abs_out
#undef n_pulse
#undef isfirst
#undef delay
#undef jitter
#undef nslit
#undef nu
#undef yheight
#undef radius
#undef theta_0
#undef omega
#undef height
#undef delta_y
#undef To
#undef Tg
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component PSD_afterslow2. */
  SIG_MESSAGE("PSD_afterslow2 (Init)");
#define mccompcurname  PSD_afterslow2
#define mccompcurtype  PSD_monitor
#define mccompcurindex 19
#define PSD_N mccPSD_afterslow2_PSD_N
#define PSD_p mccPSD_afterslow2_PSD_p
#define PSD_p2 mccPSD_afterslow2_PSD_p2
#define nx mccPSD_afterslow2_nx
#define ny mccPSD_afterslow2_ny
#define filename mccPSD_afterslow2_filename
#define xmin mccPSD_afterslow2_xmin
#define xmax mccPSD_afterslow2_xmax
#define ymin mccPSD_afterslow2_ymin
#define ymax mccPSD_afterslow2_ymax
#define xwidth mccPSD_afterslow2_xwidth
#define yheight mccPSD_afterslow2_yheight
#define restore_neutron mccPSD_afterslow2_restore_neutron
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
#line 14818 "./ESS_IN5_reprate.c"
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

  /* Initializations for component Lmon_afterslow2. */
  SIG_MESSAGE("Lmon_afterslow2 (Init)");
#define mccompcurname  Lmon_afterslow2
#define mccompcurtype  L_monitor
#define mccompcurindex 20
#define nL mccLmon_afterslow2_nL
#define L_N mccLmon_afterslow2_L_N
#define L_p mccLmon_afterslow2_L_p
#define L_p2 mccLmon_afterslow2_L_p2
#define filename mccLmon_afterslow2_filename
#define xmin mccLmon_afterslow2_xmin
#define xmax mccLmon_afterslow2_xmax
#define ymin mccLmon_afterslow2_ymin
#define ymax mccLmon_afterslow2_ymax
#define xwidth mccLmon_afterslow2_xwidth
#define yheight mccLmon_afterslow2_yheight
#define Lmin mccLmon_afterslow2_Lmin
#define Lmax mccLmon_afterslow2_Lmax
#define restore_neutron mccLmon_afterslow2_restore_neutron
#define nowritefile mccLmon_afterslow2_nowritefile
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
#line 14877 "./ESS_IN5_reprate.c"
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

  /* Initializations for component TOFL_afterslow2. */
  SIG_MESSAGE("TOFL_afterslow2 (Init)");
#define mccompcurname  TOFL_afterslow2
#define mccompcurtype  TOFLambda_monitor
#define mccompcurindex 21
#define nL mccTOFL_afterslow2_nL
#define nt mccTOFL_afterslow2_nt
#define tmin mccTOFL_afterslow2_tmin
#define tmax mccTOFL_afterslow2_tmax
#define tt_0 mccTOFL_afterslow2_tt_0
#define tt_1 mccTOFL_afterslow2_tt_1
#define TOFL_N mccTOFL_afterslow2_TOFL_N
#define TOFL_p mccTOFL_afterslow2_TOFL_p
#define TOFL_p2 mccTOFL_afterslow2_TOFL_p2
#define filename mccTOFL_afterslow2_filename
#define xmin mccTOFL_afterslow2_xmin
#define xmax mccTOFL_afterslow2_xmax
#define ymin mccTOFL_afterslow2_ymin
#define ymax mccTOFL_afterslow2_ymax
#define xwidth mccTOFL_afterslow2_xwidth
#define yheight mccTOFL_afterslow2_yheight
#define Lmin mccTOFL_afterslow2_Lmin
#define Lmax mccTOFL_afterslow2_Lmax
#define restore_neutron mccTOFL_afterslow2_restore_neutron
#define nowritefile mccTOFL_afterslow2_nowritefile
#line 61 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/TOFLambda_monitor.comp"
{
int i,j;

if (xwidth  > 0) { xmax = xwidth/2;  xmin = -xmax; }
    if (yheight > 0) { ymax = yheight/2; ymin = -ymax; }

    if ((xmin >= xmax) || (ymin >= ymax)) {
            printf("TOFlambda_monitor: %s: Null detection area !\n"
                   "ERROR              (xwidth,yheight,xmin,xmax,ymin,ymax). Exiting",
           NAME_CURRENT_COMP);
      exit(0);
    }

    tt_0 = tmin*1e-6;
    tt_1 = tmax*1e-6;
    for (i=0; i<nL; i++)
     for (j=0; j<nt; j++)
     {
      TOFL_N[j][i] = 0;
      TOFL_p[j][i] = 0;
      TOFL_p2[j][i] = 0;
     }
}
#line 14946 "./ESS_IN5_reprate.c"
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
#undef TOFL_p2
#undef TOFL_p
#undef TOFL_N
#undef tt_1
#undef tt_0
#undef tmax
#undef tmin
#undef nt
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component Guidelong2. */
  SIG_MESSAGE("Guidelong2 (Init)");
#define mccompcurname  Guidelong2
#define mccompcurtype  Guide
#define mccompcurindex 22
#define pTable mccGuidelong2_pTable
#define reflect mccGuidelong2_reflect
#define w1 mccGuidelong2_w1
#define h1 mccGuidelong2_h1
#define w2 mccGuidelong2_w2
#define h2 mccGuidelong2_h2
#define l mccGuidelong2_l
#define R0 mccGuidelong2_R0
#define Qc mccGuidelong2_Qc
#define alpha mccGuidelong2_alpha
#define m mccGuidelong2_m
#define W mccGuidelong2_W
#line 73 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide.comp"
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
#line 15007 "./ESS_IN5_reprate.c"
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

  /* Initializations for component Lmon_beforeballistic. */
  SIG_MESSAGE("Lmon_beforeballistic (Init)");
#define mccompcurname  Lmon_beforeballistic
#define mccompcurtype  L_monitor
#define mccompcurindex 23
#define nL mccLmon_beforeballistic_nL
#define L_N mccLmon_beforeballistic_L_N
#define L_p mccLmon_beforeballistic_L_p
#define L_p2 mccLmon_beforeballistic_L_p2
#define filename mccLmon_beforeballistic_filename
#define xmin mccLmon_beforeballistic_xmin
#define xmax mccLmon_beforeballistic_xmax
#define ymin mccLmon_beforeballistic_ymin
#define ymax mccLmon_beforeballistic_ymax
#define xwidth mccLmon_beforeballistic_xwidth
#define yheight mccLmon_beforeballistic_yheight
#define Lmin mccLmon_beforeballistic_Lmin
#define Lmax mccLmon_beforeballistic_Lmax
#define restore_neutron mccLmon_beforeballistic_restore_neutron
#define nowritefile mccLmon_beforeballistic_nowritefile
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
#line 15065 "./ESS_IN5_reprate.c"
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

  /* Initializations for component PSD_beforeballistic. */
  SIG_MESSAGE("PSD_beforeballistic (Init)");
#define mccompcurname  PSD_beforeballistic
#define mccompcurtype  PSD_monitor
#define mccompcurindex 24
#define PSD_N mccPSD_beforeballistic_PSD_N
#define PSD_p mccPSD_beforeballistic_PSD_p
#define PSD_p2 mccPSD_beforeballistic_PSD_p2
#define nx mccPSD_beforeballistic_nx
#define ny mccPSD_beforeballistic_ny
#define filename mccPSD_beforeballistic_filename
#define xmin mccPSD_beforeballistic_xmin
#define xmax mccPSD_beforeballistic_xmax
#define ymin mccPSD_beforeballistic_ymin
#define ymax mccPSD_beforeballistic_ymax
#define xwidth mccPSD_beforeballistic_xwidth
#define yheight mccPSD_beforeballistic_yheight
#define restore_neutron mccPSD_beforeballistic_restore_neutron
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
#line 15128 "./ESS_IN5_reprate.c"
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

  /* Initializations for component Guidelong2a. */
  SIG_MESSAGE("Guidelong2a (Init)");
#define mccompcurname  Guidelong2a
#define mccompcurtype  Guide
#define mccompcurindex 25
#define pTable mccGuidelong2a_pTable
#define reflect mccGuidelong2a_reflect
#define w1 mccGuidelong2a_w1
#define h1 mccGuidelong2a_h1
#define w2 mccGuidelong2a_w2
#define h2 mccGuidelong2a_h2
#define l mccGuidelong2a_l
#define R0 mccGuidelong2a_R0
#define Qc mccGuidelong2a_Qc
#define alpha mccGuidelong2a_alpha
#define m mccGuidelong2a_m
#define W mccGuidelong2a_W
#line 73 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide.comp"
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
#line 15182 "./ESS_IN5_reprate.c"
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

  /* Initializations for component Lmonfast2. */
  SIG_MESSAGE("Lmonfast2 (Init)");
#define mccompcurname  Lmonfast2
#define mccompcurtype  L_monitor
#define mccompcurindex 26
#define nL mccLmonfast2_nL
#define L_N mccLmonfast2_L_N
#define L_p mccLmonfast2_L_p
#define L_p2 mccLmonfast2_L_p2
#define filename mccLmonfast2_filename
#define xmin mccLmonfast2_xmin
#define xmax mccLmonfast2_xmax
#define ymin mccLmonfast2_ymin
#define ymax mccLmonfast2_ymax
#define xwidth mccLmonfast2_xwidth
#define yheight mccLmonfast2_yheight
#define Lmin mccLmonfast2_Lmin
#define Lmax mccLmonfast2_Lmax
#define restore_neutron mccLmonfast2_restore_neutron
#define nowritefile mccLmonfast2_nowritefile
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
#line 15240 "./ESS_IN5_reprate.c"
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

  /* Initializations for component Lmonfast2_zoom. */
  SIG_MESSAGE("Lmonfast2_zoom (Init)");
#define mccompcurname  Lmonfast2_zoom
#define mccompcurtype  L_monitor
#define mccompcurindex 27
#define nL mccLmonfast2_zoom_nL
#define L_N mccLmonfast2_zoom_L_N
#define L_p mccLmonfast2_zoom_L_p
#define L_p2 mccLmonfast2_zoom_L_p2
#define filename mccLmonfast2_zoom_filename
#define xmin mccLmonfast2_zoom_xmin
#define xmax mccLmonfast2_zoom_xmax
#define ymin mccLmonfast2_zoom_ymin
#define ymax mccLmonfast2_zoom_ymax
#define xwidth mccLmonfast2_zoom_xwidth
#define yheight mccLmonfast2_zoom_yheight
#define Lmin mccLmonfast2_zoom_Lmin
#define Lmax mccLmonfast2_zoom_Lmax
#define restore_neutron mccLmonfast2_zoom_restore_neutron
#define nowritefile mccLmonfast2_zoom_nowritefile
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
#line 15301 "./ESS_IN5_reprate.c"
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

  /* Initializations for component TOFLfast2. */
  SIG_MESSAGE("TOFLfast2 (Init)");
#define mccompcurname  TOFLfast2
#define mccompcurtype  TOFLambda_monitor
#define mccompcurindex 28
#define nL mccTOFLfast2_nL
#define nt mccTOFLfast2_nt
#define tmin mccTOFLfast2_tmin
#define tmax mccTOFLfast2_tmax
#define tt_0 mccTOFLfast2_tt_0
#define tt_1 mccTOFLfast2_tt_1
#define TOFL_N mccTOFLfast2_TOFL_N
#define TOFL_p mccTOFLfast2_TOFL_p
#define TOFL_p2 mccTOFLfast2_TOFL_p2
#define filename mccTOFLfast2_filename
#define xmin mccTOFLfast2_xmin
#define xmax mccTOFLfast2_xmax
#define ymin mccTOFLfast2_ymin
#define ymax mccTOFLfast2_ymax
#define xwidth mccTOFLfast2_xwidth
#define yheight mccTOFLfast2_yheight
#define Lmin mccTOFLfast2_Lmin
#define Lmax mccTOFLfast2_Lmax
#define restore_neutron mccTOFLfast2_restore_neutron
#define nowritefile mccTOFLfast2_nowritefile
#line 61 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/TOFLambda_monitor.comp"
{
int i,j;

if (xwidth  > 0) { xmax = xwidth/2;  xmin = -xmax; }
    if (yheight > 0) { ymax = yheight/2; ymin = -ymax; }

    if ((xmin >= xmax) || (ymin >= ymax)) {
            printf("TOFlambda_monitor: %s: Null detection area !\n"
                   "ERROR              (xwidth,yheight,xmin,xmax,ymin,ymax). Exiting",
           NAME_CURRENT_COMP);
      exit(0);
    }

    tt_0 = tmin*1e-6;
    tt_1 = tmax*1e-6;
    for (i=0; i<nL; i++)
     for (j=0; j<nt; j++)
     {
      TOFL_N[j][i] = 0;
      TOFL_p[j][i] = 0;
      TOFL_p2[j][i] = 0;
     }
}
#line 15370 "./ESS_IN5_reprate.c"
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
#undef TOFL_p2
#undef TOFL_p
#undef TOFL_N
#undef tt_1
#undef tt_0
#undef tmax
#undef tmin
#undef nt
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component TOFLfast2zoom. */
  SIG_MESSAGE("TOFLfast2zoom (Init)");
#define mccompcurname  TOFLfast2zoom
#define mccompcurtype  TOFLambda_monitor
#define mccompcurindex 29
#define nL mccTOFLfast2zoom_nL
#define nt mccTOFLfast2zoom_nt
#define tmin mccTOFLfast2zoom_tmin
#define tmax mccTOFLfast2zoom_tmax
#define tt_0 mccTOFLfast2zoom_tt_0
#define tt_1 mccTOFLfast2zoom_tt_1
#define TOFL_N mccTOFLfast2zoom_TOFL_N
#define TOFL_p mccTOFLfast2zoom_TOFL_p
#define TOFL_p2 mccTOFLfast2zoom_TOFL_p2
#define filename mccTOFLfast2zoom_filename
#define xmin mccTOFLfast2zoom_xmin
#define xmax mccTOFLfast2zoom_xmax
#define ymin mccTOFLfast2zoom_ymin
#define ymax mccTOFLfast2zoom_ymax
#define xwidth mccTOFLfast2zoom_xwidth
#define yheight mccTOFLfast2zoom_yheight
#define Lmin mccTOFLfast2zoom_Lmin
#define Lmax mccTOFLfast2zoom_Lmax
#define restore_neutron mccTOFLfast2zoom_restore_neutron
#define nowritefile mccTOFLfast2zoom_nowritefile
#line 61 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/TOFLambda_monitor.comp"
{
int i,j;

if (xwidth  > 0) { xmax = xwidth/2;  xmin = -xmax; }
    if (yheight > 0) { ymax = yheight/2; ymin = -ymax; }

    if ((xmin >= xmax) || (ymin >= ymax)) {
            printf("TOFlambda_monitor: %s: Null detection area !\n"
                   "ERROR              (xwidth,yheight,xmin,xmax,ymin,ymax). Exiting",
           NAME_CURRENT_COMP);
      exit(0);
    }

    tt_0 = tmin*1e-6;
    tt_1 = tmax*1e-6;
    for (i=0; i<nL; i++)
     for (j=0; j<nt; j++)
     {
      TOFL_N[j][i] = 0;
      TOFL_p[j][i] = 0;
      TOFL_p2[j][i] = 0;
     }
}
#line 15444 "./ESS_IN5_reprate.c"
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
#undef TOFL_p2
#undef TOFL_p
#undef TOFL_N
#undef tt_1
#undef tt_0
#undef tmax
#undef tmin
#undef nt
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component PSDfast2. */
  SIG_MESSAGE("PSDfast2 (Init)");
#define mccompcurname  PSDfast2
#define mccompcurtype  PSD_monitor
#define mccompcurindex 30
#define PSD_N mccPSDfast2_PSD_N
#define PSD_p mccPSDfast2_PSD_p
#define PSD_p2 mccPSDfast2_PSD_p2
#define nx mccPSDfast2_nx
#define ny mccPSDfast2_ny
#define filename mccPSDfast2_filename
#define xmin mccPSDfast2_xmin
#define xmax mccPSDfast2_xmax
#define ymin mccPSDfast2_ymin
#define ymax mccPSDfast2_ymax
#define xwidth mccPSDfast2_xwidth
#define yheight mccPSDfast2_yheight
#define restore_neutron mccPSDfast2_restore_neutron
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
#line 15512 "./ESS_IN5_reprate.c"
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

  /* Initializations for component Fastchop2. */
  SIG_MESSAGE("Fastchop2 (Init)");
#define mccompcurname  Fastchop2
#define mccompcurtype  DiskChopper
#define mccompcurindex 31
#define Tg mccFastchop2_Tg
#define To mccFastchop2_To
#define delta_y mccFastchop2_delta_y
#define height mccFastchop2_height
#define omega mccFastchop2_omega
#define theta_0 mccFastchop2_theta_0
#define radius mccFastchop2_radius
#define yheight mccFastchop2_yheight
#define nu mccFastchop2_nu
#define nslit mccFastchop2_nslit
#define jitter mccFastchop2_jitter
#define delay mccFastchop2_delay
#define isfirst mccFastchop2_isfirst
#define n_pulse mccFastchop2_n_pulse
#define abs_out mccFastchop2_abs_out
#define phase mccFastchop2_phase
#define xwidth mccFastchop2_xwidth
#define verbose mccFastchop2_verbose
#line 70 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/DiskChopper.comp"
{
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
}
#line 15615 "./ESS_IN5_reprate.c"
#undef verbose
#undef xwidth
#undef phase
#undef abs_out
#undef n_pulse
#undef isfirst
#undef delay
#undef jitter
#undef nslit
#undef nu
#undef yheight
#undef radius
#undef theta_0
#undef omega
#undef height
#undef delta_y
#undef To
#undef Tg
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component Fastchop2counter. */
  SIG_MESSAGE("Fastchop2counter (Init)");
#define mccompcurname  Fastchop2counter
#define mccompcurtype  DiskChopper
#define mccompcurindex 32
#define Tg mccFastchop2counter_Tg
#define To mccFastchop2counter_To
#define delta_y mccFastchop2counter_delta_y
#define height mccFastchop2counter_height
#define omega mccFastchop2counter_omega
#define theta_0 mccFastchop2counter_theta_0
#define radius mccFastchop2counter_radius
#define yheight mccFastchop2counter_yheight
#define nu mccFastchop2counter_nu
#define nslit mccFastchop2counter_nslit
#define jitter mccFastchop2counter_jitter
#define delay mccFastchop2counter_delay
#define isfirst mccFastchop2counter_isfirst
#define n_pulse mccFastchop2counter_n_pulse
#define abs_out mccFastchop2counter_abs_out
#define phase mccFastchop2counter_phase
#define xwidth mccFastchop2counter_xwidth
#define verbose mccFastchop2counter_verbose
#line 70 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/DiskChopper.comp"
{
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
}
#line 15723 "./ESS_IN5_reprate.c"
#undef verbose
#undef xwidth
#undef phase
#undef abs_out
#undef n_pulse
#undef isfirst
#undef delay
#undef jitter
#undef nslit
#undef nu
#undef yheight
#undef radius
#undef theta_0
#undef omega
#undef height
#undef delta_y
#undef To
#undef Tg
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component FOchop3. */
  SIG_MESSAGE("FOchop3 (Init)");
#define mccompcurname  FOchop3
#define mccompcurtype  DiskChopper
#define mccompcurindex 33
#define Tg mccFOchop3_Tg
#define To mccFOchop3_To
#define delta_y mccFOchop3_delta_y
#define height mccFOchop3_height
#define omega mccFOchop3_omega
#define theta_0 mccFOchop3_theta_0
#define radius mccFOchop3_radius
#define yheight mccFOchop3_yheight
#define nu mccFOchop3_nu
#define nslit mccFOchop3_nslit
#define jitter mccFOchop3_jitter
#define delay mccFOchop3_delay
#define isfirst mccFOchop3_isfirst
#define n_pulse mccFOchop3_n_pulse
#define abs_out mccFOchop3_abs_out
#define phase mccFOchop3_phase
#define xwidth mccFOchop3_xwidth
#define verbose mccFOchop3_verbose
#line 70 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/DiskChopper.comp"
{
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
}
#line 15831 "./ESS_IN5_reprate.c"
#undef verbose
#undef xwidth
#undef phase
#undef abs_out
#undef n_pulse
#undef isfirst
#undef delay
#undef jitter
#undef nslit
#undef nu
#undef yheight
#undef radius
#undef theta_0
#undef omega
#undef height
#undef delta_y
#undef To
#undef Tg
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component TOFfast2_zoom. */
  SIG_MESSAGE("TOFfast2_zoom (Init)");
#define mccompcurname  TOFfast2_zoom
#define mccompcurtype  TOF_monitor
#define mccompcurindex 34
#define nt mccTOFfast2_zoom_nt
#define TOF_N mccTOFfast2_zoom_TOF_N
#define TOF_p mccTOFfast2_zoom_TOF_p
#define TOF_p2 mccTOFfast2_zoom_TOF_p2
#define t_min mccTOFfast2_zoom_t_min
#define t_max mccTOFfast2_zoom_t_max
#define delta_t mccTOFfast2_zoom_delta_t
#define filename mccTOFfast2_zoom_filename
#define xmin mccTOFfast2_zoom_xmin
#define xmax mccTOFfast2_zoom_xmax
#define ymin mccTOFfast2_zoom_ymin
#define ymax mccTOFfast2_zoom_ymax
#define xwidth mccTOFfast2_zoom_xwidth
#define yheight mccTOFfast2_zoom_yheight
#define tmin mccTOFfast2_zoom_tmin
#define tmax mccTOFfast2_zoom_tmax
#define dt mccTOFfast2_zoom_dt
#define restore_neutron mccTOFfast2_zoom_restore_neutron
#define nowritefile mccTOFfast2_zoom_nowritefile
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
#line 15911 "./ESS_IN5_reprate.c"
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

  /* Initializations for component Lmon_afterfast2. */
  SIG_MESSAGE("Lmon_afterfast2 (Init)");
#define mccompcurname  Lmon_afterfast2
#define mccompcurtype  L_monitor
#define mccompcurindex 35
#define nL mccLmon_afterfast2_nL
#define L_N mccLmon_afterfast2_L_N
#define L_p mccLmon_afterfast2_L_p
#define L_p2 mccLmon_afterfast2_L_p2
#define filename mccLmon_afterfast2_filename
#define xmin mccLmon_afterfast2_xmin
#define xmax mccLmon_afterfast2_xmax
#define ymin mccLmon_afterfast2_ymin
#define ymax mccLmon_afterfast2_ymax
#define xwidth mccLmon_afterfast2_xwidth
#define yheight mccLmon_afterfast2_yheight
#define Lmin mccLmon_afterfast2_Lmin
#define Lmax mccLmon_afterfast2_Lmax
#define restore_neutron mccLmon_afterfast2_restore_neutron
#define nowritefile mccLmon_afterfast2_nowritefile
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
#line 15976 "./ESS_IN5_reprate.c"
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

  /* Initializations for component TOFL_afterfast2. */
  SIG_MESSAGE("TOFL_afterfast2 (Init)");
#define mccompcurname  TOFL_afterfast2
#define mccompcurtype  TOFLambda_monitor
#define mccompcurindex 36
#define nL mccTOFL_afterfast2_nL
#define nt mccTOFL_afterfast2_nt
#define tmin mccTOFL_afterfast2_tmin
#define tmax mccTOFL_afterfast2_tmax
#define tt_0 mccTOFL_afterfast2_tt_0
#define tt_1 mccTOFL_afterfast2_tt_1
#define TOFL_N mccTOFL_afterfast2_TOFL_N
#define TOFL_p mccTOFL_afterfast2_TOFL_p
#define TOFL_p2 mccTOFL_afterfast2_TOFL_p2
#define filename mccTOFL_afterfast2_filename
#define xmin mccTOFL_afterfast2_xmin
#define xmax mccTOFL_afterfast2_xmax
#define ymin mccTOFL_afterfast2_ymin
#define ymax mccTOFL_afterfast2_ymax
#define xwidth mccTOFL_afterfast2_xwidth
#define yheight mccTOFL_afterfast2_yheight
#define Lmin mccTOFL_afterfast2_Lmin
#define Lmax mccTOFL_afterfast2_Lmax
#define restore_neutron mccTOFL_afterfast2_restore_neutron
#define nowritefile mccTOFL_afterfast2_nowritefile
#line 61 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/TOFLambda_monitor.comp"
{
int i,j;

if (xwidth  > 0) { xmax = xwidth/2;  xmin = -xmax; }
    if (yheight > 0) { ymax = yheight/2; ymin = -ymax; }

    if ((xmin >= xmax) || (ymin >= ymax)) {
            printf("TOFlambda_monitor: %s: Null detection area !\n"
                   "ERROR              (xwidth,yheight,xmin,xmax,ymin,ymax). Exiting",
           NAME_CURRENT_COMP);
      exit(0);
    }

    tt_0 = tmin*1e-6;
    tt_1 = tmax*1e-6;
    for (i=0; i<nL; i++)
     for (j=0; j<nt; j++)
     {
      TOFL_N[j][i] = 0;
      TOFL_p[j][i] = 0;
      TOFL_p2[j][i] = 0;
     }
}
#line 16045 "./ESS_IN5_reprate.c"
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
#undef TOFL_p2
#undef TOFL_p
#undef TOFL_N
#undef tt_1
#undef tt_0
#undef tmax
#undef tmin
#undef nt
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component TOFL_afterfast2_zoom. */
  SIG_MESSAGE("TOFL_afterfast2_zoom (Init)");
#define mccompcurname  TOFL_afterfast2_zoom
#define mccompcurtype  TOFLambda_monitor
#define mccompcurindex 37
#define nL mccTOFL_afterfast2_zoom_nL
#define nt mccTOFL_afterfast2_zoom_nt
#define tmin mccTOFL_afterfast2_zoom_tmin
#define tmax mccTOFL_afterfast2_zoom_tmax
#define tt_0 mccTOFL_afterfast2_zoom_tt_0
#define tt_1 mccTOFL_afterfast2_zoom_tt_1
#define TOFL_N mccTOFL_afterfast2_zoom_TOFL_N
#define TOFL_p mccTOFL_afterfast2_zoom_TOFL_p
#define TOFL_p2 mccTOFL_afterfast2_zoom_TOFL_p2
#define filename mccTOFL_afterfast2_zoom_filename
#define xmin mccTOFL_afterfast2_zoom_xmin
#define xmax mccTOFL_afterfast2_zoom_xmax
#define ymin mccTOFL_afterfast2_zoom_ymin
#define ymax mccTOFL_afterfast2_zoom_ymax
#define xwidth mccTOFL_afterfast2_zoom_xwidth
#define yheight mccTOFL_afterfast2_zoom_yheight
#define Lmin mccTOFL_afterfast2_zoom_Lmin
#define Lmax mccTOFL_afterfast2_zoom_Lmax
#define restore_neutron mccTOFL_afterfast2_zoom_restore_neutron
#define nowritefile mccTOFL_afterfast2_zoom_nowritefile
#line 61 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/TOFLambda_monitor.comp"
{
int i,j;

if (xwidth  > 0) { xmax = xwidth/2;  xmin = -xmax; }
    if (yheight > 0) { ymax = yheight/2; ymin = -ymax; }

    if ((xmin >= xmax) || (ymin >= ymax)) {
            printf("TOFlambda_monitor: %s: Null detection area !\n"
                   "ERROR              (xwidth,yheight,xmin,xmax,ymin,ymax). Exiting",
           NAME_CURRENT_COMP);
      exit(0);
    }

    tt_0 = tmin*1e-6;
    tt_1 = tmax*1e-6;
    for (i=0; i<nL; i++)
     for (j=0; j<nt; j++)
     {
      TOFL_N[j][i] = 0;
      TOFL_p[j][i] = 0;
      TOFL_p2[j][i] = 0;
     }
}
#line 16119 "./ESS_IN5_reprate.c"
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
#undef TOFL_p2
#undef TOFL_p
#undef TOFL_N
#undef tt_1
#undef tt_0
#undef tmax
#undef tmin
#undef nt
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component PSD_afterfast2. */
  SIG_MESSAGE("PSD_afterfast2 (Init)");
#define mccompcurname  PSD_afterfast2
#define mccompcurtype  PSD_monitor
#define mccompcurindex 38
#define PSD_N mccPSD_afterfast2_PSD_N
#define PSD_p mccPSD_afterfast2_PSD_p
#define PSD_p2 mccPSD_afterfast2_PSD_p2
#define nx mccPSD_afterfast2_nx
#define ny mccPSD_afterfast2_ny
#define filename mccPSD_afterfast2_filename
#define xmin mccPSD_afterfast2_xmin
#define xmax mccPSD_afterfast2_xmax
#define ymin mccPSD_afterfast2_ymin
#define ymax mccPSD_afterfast2_ymax
#define xwidth mccPSD_afterfast2_xwidth
#define yheight mccPSD_afterfast2_yheight
#define restore_neutron mccPSD_afterfast2_restore_neutron
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
#line 16187 "./ESS_IN5_reprate.c"
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

  /* Initializations for component Guidesample. */
  SIG_MESSAGE("Guidesample (Init)");
#define mccompcurname  Guidesample
#define mccompcurtype  Guide
#define mccompcurindex 39
#define pTable mccGuidesample_pTable
#define reflect mccGuidesample_reflect
#define w1 mccGuidesample_w1
#define h1 mccGuidesample_h1
#define w2 mccGuidesample_w2
#define h2 mccGuidesample_h2
#define l mccGuidesample_l
#define R0 mccGuidesample_R0
#define Qc mccGuidesample_Qc
#define alpha mccGuidesample_alpha
#define m mccGuidesample_m
#define W mccGuidesample_W
#line 73 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide.comp"
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
#line 16241 "./ESS_IN5_reprate.c"
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

  /* Initializations for component Lmon_guideend. */
  SIG_MESSAGE("Lmon_guideend (Init)");
#define mccompcurname  Lmon_guideend
#define mccompcurtype  L_monitor
#define mccompcurindex 40
#define nL mccLmon_guideend_nL
#define L_N mccLmon_guideend_L_N
#define L_p mccLmon_guideend_L_p
#define L_p2 mccLmon_guideend_L_p2
#define filename mccLmon_guideend_filename
#define xmin mccLmon_guideend_xmin
#define xmax mccLmon_guideend_xmax
#define ymin mccLmon_guideend_ymin
#define ymax mccLmon_guideend_ymax
#define xwidth mccLmon_guideend_xwidth
#define yheight mccLmon_guideend_yheight
#define Lmin mccLmon_guideend_Lmin
#define Lmax mccLmon_guideend_Lmax
#define restore_neutron mccLmon_guideend_restore_neutron
#define nowritefile mccLmon_guideend_nowritefile
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
#line 16299 "./ESS_IN5_reprate.c"
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

  /* Initializations for component PSDsample. */
  SIG_MESSAGE("PSDsample (Init)");
#define mccompcurname  PSDsample
#define mccompcurtype  PSD_monitor
#define mccompcurindex 41
#define PSD_N mccPSDsample_PSD_N
#define PSD_p mccPSDsample_PSD_p
#define PSD_p2 mccPSDsample_PSD_p2
#define nx mccPSDsample_nx
#define ny mccPSDsample_ny
#define filename mccPSDsample_filename
#define xmin mccPSDsample_xmin
#define xmax mccPSDsample_xmax
#define ymin mccPSDsample_ymin
#define ymax mccPSDsample_ymax
#define xwidth mccPSDsample_xwidth
#define yheight mccPSDsample_yheight
#define restore_neutron mccPSDsample_restore_neutron
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
#line 16362 "./ESS_IN5_reprate.c"
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

  /* Initializations for component TOFsample_zoom. */
  SIG_MESSAGE("TOFsample_zoom (Init)");
#define mccompcurname  TOFsample_zoom
#define mccompcurtype  TOF_monitor
#define mccompcurindex 42
#define nt mccTOFsample_zoom_nt
#define TOF_N mccTOFsample_zoom_TOF_N
#define TOF_p mccTOFsample_zoom_TOF_p
#define TOF_p2 mccTOFsample_zoom_TOF_p2
#define t_min mccTOFsample_zoom_t_min
#define t_max mccTOFsample_zoom_t_max
#define delta_t mccTOFsample_zoom_delta_t
#define filename mccTOFsample_zoom_filename
#define xmin mccTOFsample_zoom_xmin
#define xmax mccTOFsample_zoom_xmax
#define ymin mccTOFsample_zoom_ymin
#define ymax mccTOFsample_zoom_ymax
#define xwidth mccTOFsample_zoom_xwidth
#define yheight mccTOFsample_zoom_yheight
#define tmin mccTOFsample_zoom_tmin
#define tmax mccTOFsample_zoom_tmax
#define dt mccTOFsample_zoom_dt
#define restore_neutron mccTOFsample_zoom_restore_neutron
#define nowritefile mccTOFsample_zoom_nowritefile
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
#line 16437 "./ESS_IN5_reprate.c"
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

  /* Initializations for component Esample. */
  SIG_MESSAGE("Esample (Init)");
#define mccompcurname  Esample
#define mccompcurtype  E_monitor
#define mccompcurindex 43
#define nE mccEsample_nE
#define E_N mccEsample_E_N
#define E_p mccEsample_E_p
#define E_p2 mccEsample_E_p2
#define S_p mccEsample_S_p
#define S_pE mccEsample_S_pE
#define S_pE2 mccEsample_S_pE2
#define filename mccEsample_filename
#define xmin mccEsample_xmin
#define xmax mccEsample_xmax
#define ymin mccEsample_ymin
#define ymax mccEsample_ymax
#define xwidth mccEsample_xwidth
#define yheight mccEsample_yheight
#define Emin mccEsample_Emin
#define Emax mccEsample_Emax
#define restore_neutron mccEsample_restore_neutron
#define nowritefile mccEsample_nowritefile
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
#line 16506 "./ESS_IN5_reprate.c"
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

  /* Initializations for component Lmon_sample_zoom. */
  SIG_MESSAGE("Lmon_sample_zoom (Init)");
#define mccompcurname  Lmon_sample_zoom
#define mccompcurtype  L_monitor
#define mccompcurindex 44
#define nL mccLmon_sample_zoom_nL
#define L_N mccLmon_sample_zoom_L_N
#define L_p mccLmon_sample_zoom_L_p
#define L_p2 mccLmon_sample_zoom_L_p2
#define filename mccLmon_sample_zoom_filename
#define xmin mccLmon_sample_zoom_xmin
#define xmax mccLmon_sample_zoom_xmax
#define ymin mccLmon_sample_zoom_ymin
#define ymax mccLmon_sample_zoom_ymax
#define xwidth mccLmon_sample_zoom_xwidth
#define yheight mccLmon_sample_zoom_yheight
#define Lmin mccLmon_sample_zoom_Lmin
#define Lmax mccLmon_sample_zoom_Lmax
#define restore_neutron mccLmon_sample_zoom_restore_neutron
#define nowritefile mccLmon_sample_zoom_nowritefile
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
#line 16570 "./ESS_IN5_reprate.c"
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

  /* Initializations for component sample. */
  SIG_MESSAGE("sample (Init)");
#define mccompcurname  sample
#define mccompcurtype  Tunneling_sample
#define mccompcurindex 45
#define ftun mccsample_ftun
#define fQE mccsample_fQE
#define VarsV mccsample_VarsV
#define thickness mccsample_thickness
#define radius mccsample_radius
#define focus_r mccsample_focus_r
#define p_interact mccsample_p_interact
#define f_QE mccsample_f_QE
#define f_tun mccsample_f_tun
#define gamma mccsample_gamma
#define E_tun mccsample_E_tun
#define target_x mccsample_target_x
#define target_y mccsample_target_y
#define target_z mccsample_target_z
#define focus_xw mccsample_focus_xw
#define focus_yh mccsample_focus_yh
#define focus_aw mccsample_focus_aw
#define focus_ah mccsample_focus_ah
#define xwidth mccsample_xwidth
#define yheight mccsample_yheight
#define zdepth mccsample_zdepth
#define sigma_abs mccsample_sigma_abs
#define sigma_inc mccsample_sigma_inc
#define Vc mccsample_Vc
#define target_index mccsample_target_index
#line 114 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../samples/Tunneling_sample.comp"
{
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
}
#line 16690 "./ESS_IN5_reprate.c"
#undef target_index
#undef Vc
#undef sigma_inc
#undef sigma_abs
#undef zdepth
#undef yheight
#undef xwidth
#undef focus_ah
#undef focus_aw
#undef focus_yh
#undef focus_xw
#undef target_z
#undef target_y
#undef target_x
#undef E_tun
#undef gamma
#undef f_tun
#undef f_QE
#undef p_interact
#undef focus_r
#undef radius
#undef thickness
#undef VarsV
#undef fQE
#undef ftun
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component detectorarm. */
  SIG_MESSAGE("detectorarm (Init)");

  /* Initializations for component TOFdetector. */
  SIG_MESSAGE("TOFdetector (Init)");
#define mccompcurname  TOFdetector
#define mccompcurtype  TOF_monitor
#define mccompcurindex 47
#define nt mccTOFdetector_nt
#define TOF_N mccTOFdetector_TOF_N
#define TOF_p mccTOFdetector_TOF_p
#define TOF_p2 mccTOFdetector_TOF_p2
#define t_min mccTOFdetector_t_min
#define t_max mccTOFdetector_t_max
#define delta_t mccTOFdetector_delta_t
#define filename mccTOFdetector_filename
#define xmin mccTOFdetector_xmin
#define xmax mccTOFdetector_xmax
#define ymin mccTOFdetector_ymin
#define ymax mccTOFdetector_ymax
#define xwidth mccTOFdetector_xwidth
#define yheight mccTOFdetector_yheight
#define tmin mccTOFdetector_tmin
#define tmax mccTOFdetector_tmax
#define dt mccTOFdetector_dt
#define restore_neutron mccTOFdetector_restore_neutron
#define nowritefile mccTOFdetector_nowritefile
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
#line 16780 "./ESS_IN5_reprate.c"
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

  /* Initializations for component TOFdetector_zoom. */
  SIG_MESSAGE("TOFdetector_zoom (Init)");
#define mccompcurname  TOFdetector_zoom
#define mccompcurtype  TOF_monitor
#define mccompcurindex 48
#define nt mccTOFdetector_zoom_nt
#define TOF_N mccTOFdetector_zoom_TOF_N
#define TOF_p mccTOFdetector_zoom_TOF_p
#define TOF_p2 mccTOFdetector_zoom_TOF_p2
#define t_min mccTOFdetector_zoom_t_min
#define t_max mccTOFdetector_zoom_t_max
#define delta_t mccTOFdetector_zoom_delta_t
#define filename mccTOFdetector_zoom_filename
#define xmin mccTOFdetector_zoom_xmin
#define xmax mccTOFdetector_zoom_xmax
#define ymin mccTOFdetector_zoom_ymin
#define ymax mccTOFdetector_zoom_ymax
#define xwidth mccTOFdetector_zoom_xwidth
#define yheight mccTOFdetector_zoom_yheight
#define tmin mccTOFdetector_zoom_tmin
#define tmax mccTOFdetector_zoom_tmax
#define dt mccTOFdetector_zoom_dt
#define restore_neutron mccTOFdetector_zoom_restore_neutron
#define nowritefile mccTOFdetector_zoom_nowritefile
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
#line 16861 "./ESS_IN5_reprate.c"
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

  /* Initializations for component Edetector. */
  SIG_MESSAGE("Edetector (Init)");
#define mccompcurname  Edetector
#define mccompcurtype  E_monitor
#define mccompcurindex 49
#define nE mccEdetector_nE
#define E_N mccEdetector_E_N
#define E_p mccEdetector_E_p
#define E_p2 mccEdetector_E_p2
#define S_p mccEdetector_S_p
#define S_pE mccEdetector_S_pE
#define S_pE2 mccEdetector_S_pE2
#define filename mccEdetector_filename
#define xmin mccEdetector_xmin
#define xmax mccEdetector_xmax
#define ymin mccEdetector_ymin
#define ymax mccEdetector_ymax
#define xwidth mccEdetector_xwidth
#define yheight mccEdetector_yheight
#define Emin mccEdetector_Emin
#define Emax mccEdetector_Emax
#define restore_neutron mccEdetector_restore_neutron
#define nowritefile mccEdetector_nowritefile
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
#line 16930 "./ESS_IN5_reprate.c"
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

  /* Initializations for component TOF2Edetector. */
  SIG_MESSAGE("TOF2Edetector (Init)");
#define mccompcurname  TOF2Edetector
#define mccompcurtype  TOF2E_monitor
#define mccompcurindex 50
#define nE mccTOF2Edetector_nE
#define E_N mccTOF2Edetector_E_N
#define E_p mccTOF2Edetector_E_p
#define E_p2 mccTOF2Edetector_E_p2
#define S_p mccTOF2Edetector_S_p
#define S_pE mccTOF2Edetector_S_pE
#define S_pE2 mccTOF2Edetector_S_pE2
#define filename mccTOF2Edetector_filename
#define xmin mccTOF2Edetector_xmin
#define xmax mccTOF2Edetector_xmax
#define ymin mccTOF2Edetector_ymin
#define ymax mccTOF2Edetector_ymax
#define xwidth mccTOF2Edetector_xwidth
#define yheight mccTOF2Edetector_yheight
#define Emin mccTOF2Edetector_Emin
#define Emax mccTOF2Edetector_Emax
#define T_zero mccTOF2Edetector_T_zero
#define L_flight mccTOF2Edetector_L_flight
#define restore_neutron mccTOF2Edetector_restore_neutron
#define nowritefile mccTOF2Edetector_nowritefile
#line 67 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/TOF2E_monitor.comp"
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
#line 17000 "./ESS_IN5_reprate.c"
#undef nowritefile
#undef restore_neutron
#undef L_flight
#undef T_zero
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
  /* TRACE Component source [1] */
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
#define mccompcurname  source
#define mccompcurtype  ESS_moderator_long
#define mccompcurindex 1
#define l_range mccsource_l_range
#define w_mult mccsource_w_mult
#define w_geom mccsource_w_geom
#define w_geom_c mccsource_w_geom_c
#define w_geom_t mccsource_w_geom_t
#define tx mccsource_tx
#define ty mccsource_ty
#define tz mccsource_tz
#define t1x mccsource_t1x
#define t1y mccsource_t1y
#define t1z mccsource_t1z
#define t2x mccsource_t2x
#define t2y mccsource_t2y
#define t2z mccsource_t2z
#define T_n mccsource_T_n
#define tau_n mccsource_tau_n
#define tau1_n mccsource_tau1_n
#define tau2_n mccsource_tau2_n
#define chi2_n mccsource_chi2_n
#define I0_n mccsource_I0_n
#define I2_n mccsource_I2_n
#define branch1_n mccsource_branch1_n
#define branch2_n mccsource_branch2_n
#define r_empty mccsource_r_empty
#define r_optics mccsource_r_optics
{   /* Declarations of source=ESS_moderator_long() SETTING parameters. */
MCNUM width_c = mccsource_width_c;
MCNUM yheight = mccsource_yheight;
MCNUM Lmin = mccsource_Lmin;
MCNUM Lmax = mccsource_Lmax;
MCNUM dist = mccsource_dist;
MCNUM focus_xw = mccsource_focus_xw;
MCNUM focus_yh = mccsource_focus_yh;
MCNUM nu = mccsource_nu;
MCNUM T = mccsource_T;
MCNUM tau = mccsource_tau;
MCNUM tau1 = mccsource_tau1;
MCNUM tau2 = mccsource_tau2;
MCNUM d = mccsource_d;
MCNUM n = mccsource_n;
MCNUM cold_frac = mccsource_cold_frac;
MCNUM n2 = mccsource_n2;
MCNUM chi2 = mccsource_chi2;
MCNUM I0 = mccsource_I0;
MCNUM I2 = mccsource_I2;
int target_index = mccsource_target_index;
MCNUM cyl_radius = mccsource_cyl_radius;
MCNUM branch1 = mccsource_branch1;
MCNUM branch2 = mccsource_branch2;
MCNUM branch_tail = mccsource_branch_tail;
int n_pulses = mccsource_n_pulses;
MCNUM width_t = mccsource_width_t;
MCNUM T_t = mccsource_T_t;
MCNUM tau_t = mccsource_tau_t;
MCNUM tau1_t = mccsource_tau1_t;
MCNUM tau2_t = mccsource_tau2_t;
MCNUM chi2_t = mccsource_chi2_t;
MCNUM I0_t = mccsource_I0_t;
MCNUM I2_t = mccsource_I2_t;
MCNUM branch1_t = mccsource_branch1_t;
MCNUM branch2_t = mccsource_branch2_t;
int src_2012 = mccsource_src_2012;
MCNUM tfocus_dist = mccsource_tfocus_dist;
MCNUM tfocus_time = mccsource_tfocus_time;
MCNUM tfocus_width = mccsource_tfocus_width;
MCNUM beamport_angle = mccsource_beamport_angle;
#line 253 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../obsolete/ESS_moderator_long.comp"
{
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
}
#line 17374 "./ESS_IN5_reprate.c"
}   /* End of source=ESS_moderator_long() SETTING parameter declarations. */
#undef r_optics
#undef r_empty
#undef branch2_n
#undef branch1_n
#undef I2_n
#undef I0_n
#undef chi2_n
#undef tau2_n
#undef tau1_n
#undef tau_n
#undef T_n
#undef t2z
#undef t2y
#undef t2x
#undef t1z
#undef t1y
#undef t1x
#undef tz
#undef ty
#undef tx
#undef w_geom_t
#undef w_geom_c
#undef w_geom
#undef w_mult
#undef l_range
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompsource:
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

  /* TRACE Component Origin [2] */
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
#define mccompcurname  Origin
#define mccompcurtype  Progress_bar
#define mccompcurindex 2
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
#line 17557 "./ESS_IN5_reprate.c"
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

  /* TRACE Component TOFmoderator_zoom [3] */
  mccoordschange(mcposrTOFmoderator_zoom, mcrotrTOFmoderator_zoom,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component TOFmoderator_zoom (without coords transformations) */
  mcJumpTrace_TOFmoderator_zoom:
  SIG_MESSAGE("TOFmoderator_zoom (Trace)");
  mcDEBUG_COMP("TOFmoderator_zoom")
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

#define mcabsorbComp mcabsorbCompTOFmoderator_zoom
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
#define mccompcurname  TOFmoderator_zoom
#define mccompcurtype  TOF_monitor
#define mccompcurindex 3
#define nt mccTOFmoderator_zoom_nt
#define TOF_N mccTOFmoderator_zoom_TOF_N
#define TOF_p mccTOFmoderator_zoom_TOF_p
#define TOF_p2 mccTOFmoderator_zoom_TOF_p2
#define t_min mccTOFmoderator_zoom_t_min
#define t_max mccTOFmoderator_zoom_t_max
#define delta_t mccTOFmoderator_zoom_delta_t
{   /* Declarations of TOFmoderator_zoom=TOF_monitor() SETTING parameters. */
char* filename = mccTOFmoderator_zoom_filename;
MCNUM xmin = mccTOFmoderator_zoom_xmin;
MCNUM xmax = mccTOFmoderator_zoom_xmax;
MCNUM ymin = mccTOFmoderator_zoom_ymin;
MCNUM ymax = mccTOFmoderator_zoom_ymax;
MCNUM xwidth = mccTOFmoderator_zoom_xwidth;
MCNUM yheight = mccTOFmoderator_zoom_yheight;
MCNUM tmin = mccTOFmoderator_zoom_tmin;
MCNUM tmax = mccTOFmoderator_zoom_tmax;
MCNUM dt = mccTOFmoderator_zoom_dt;
MCNUM restore_neutron = mccTOFmoderator_zoom_restore_neutron;
int nowritefile = mccTOFmoderator_zoom_nowritefile;
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
#line 17705 "./ESS_IN5_reprate.c"
}   /* End of TOFmoderator_zoom=TOF_monitor() SETTING parameter declarations. */
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
  mcabsorbCompTOFmoderator_zoom:
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

  /* TRACE Component TOFmoderator [4] */
  mccoordschange(mcposrTOFmoderator, mcrotrTOFmoderator,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component TOFmoderator (without coords transformations) */
  mcJumpTrace_TOFmoderator:
  SIG_MESSAGE("TOFmoderator (Trace)");
  mcDEBUG_COMP("TOFmoderator")
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

#define mcabsorbComp mcabsorbCompTOFmoderator
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
#define mccompcurname  TOFmoderator
#define mccompcurtype  TOF_monitor
#define mccompcurindex 4
#define nt mccTOFmoderator_nt
#define TOF_N mccTOFmoderator_TOF_N
#define TOF_p mccTOFmoderator_TOF_p
#define TOF_p2 mccTOFmoderator_TOF_p2
#define t_min mccTOFmoderator_t_min
#define t_max mccTOFmoderator_t_max
#define delta_t mccTOFmoderator_delta_t
{   /* Declarations of TOFmoderator=TOF_monitor() SETTING parameters. */
char* filename = mccTOFmoderator_filename;
MCNUM xmin = mccTOFmoderator_xmin;
MCNUM xmax = mccTOFmoderator_xmax;
MCNUM ymin = mccTOFmoderator_ymin;
MCNUM ymax = mccTOFmoderator_ymax;
MCNUM xwidth = mccTOFmoderator_xwidth;
MCNUM yheight = mccTOFmoderator_yheight;
MCNUM tmin = mccTOFmoderator_tmin;
MCNUM tmax = mccTOFmoderator_tmax;
MCNUM dt = mccTOFmoderator_dt;
MCNUM restore_neutron = mccTOFmoderator_restore_neutron;
int nowritefile = mccTOFmoderator_nowritefile;
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
#line 17856 "./ESS_IN5_reprate.c"
}   /* End of TOFmoderator=TOF_monitor() SETTING parameter declarations. */
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
  mcabsorbCompTOFmoderator:
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

  /* TRACE Component Lmon_guistart [5] */
  mccoordschange(mcposrLmon_guistart, mcrotrLmon_guistart,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component Lmon_guistart (without coords transformations) */
  mcJumpTrace_Lmon_guistart:
  SIG_MESSAGE("Lmon_guistart (Trace)");
  mcDEBUG_COMP("Lmon_guistart")
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

#define mcabsorbComp mcabsorbCompLmon_guistart
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
#define mccompcurname  Lmon_guistart
#define mccompcurtype  L_monitor
#define mccompcurindex 5
#define nL mccLmon_guistart_nL
#define L_N mccLmon_guistart_L_N
#define L_p mccLmon_guistart_L_p
#define L_p2 mccLmon_guistart_L_p2
{   /* Declarations of Lmon_guistart=L_monitor() SETTING parameters. */
char* filename = mccLmon_guistart_filename;
MCNUM xmin = mccLmon_guistart_xmin;
MCNUM xmax = mccLmon_guistart_xmax;
MCNUM ymin = mccLmon_guistart_ymin;
MCNUM ymax = mccLmon_guistart_ymax;
MCNUM xwidth = mccLmon_guistart_xwidth;
MCNUM yheight = mccLmon_guistart_yheight;
MCNUM Lmin = mccLmon_guistart_Lmin;
MCNUM Lmax = mccLmon_guistart_Lmax;
MCNUM restore_neutron = mccLmon_guistart_restore_neutron;
int nowritefile = mccLmon_guistart_nowritefile;
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
#line 18006 "./ESS_IN5_reprate.c"
}   /* End of Lmon_guistart=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompLmon_guistart:
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

  /* TRACE Component Lmon_normalize [6] */
  mccoordschange(mcposrLmon_normalize, mcrotrLmon_normalize,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component Lmon_normalize (without coords transformations) */
  mcJumpTrace_Lmon_normalize:
  SIG_MESSAGE("Lmon_normalize (Trace)");
  mcDEBUG_COMP("Lmon_normalize")
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

#define mcabsorbComp mcabsorbCompLmon_normalize
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
#define mccompcurname  Lmon_normalize
#define mccompcurtype  L_monitor
#define mccompcurindex 6
#define nL mccLmon_normalize_nL
#define L_N mccLmon_normalize_L_N
#define L_p mccLmon_normalize_L_p
#define L_p2 mccLmon_normalize_L_p2
{   /* Declarations of Lmon_normalize=L_monitor() SETTING parameters. */
char* filename = mccLmon_normalize_filename;
MCNUM xmin = mccLmon_normalize_xmin;
MCNUM xmax = mccLmon_normalize_xmax;
MCNUM ymin = mccLmon_normalize_ymin;
MCNUM ymax = mccLmon_normalize_ymax;
MCNUM xwidth = mccLmon_normalize_xwidth;
MCNUM yheight = mccLmon_normalize_yheight;
MCNUM Lmin = mccLmon_normalize_Lmin;
MCNUM Lmax = mccLmon_normalize_Lmax;
MCNUM restore_neutron = mccLmon_normalize_restore_neutron;
int nowritefile = mccLmon_normalize_nowritefile;
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
#line 18153 "./ESS_IN5_reprate.c"
}   /* End of Lmon_normalize=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompLmon_normalize:
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

  /* TRACE Component Guide1 [7] */
  mccoordschange(mcposrGuide1, mcrotrGuide1,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component Guide1 (without coords transformations) */
  mcJumpTrace_Guide1:
  SIG_MESSAGE("Guide1 (Trace)");
  mcDEBUG_COMP("Guide1")
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

#define mcabsorbComp mcabsorbCompGuide1
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
#define mccompcurname  Guide1
#define mccompcurtype  Guide
#define mccompcurindex 7
#define pTable mccGuide1_pTable
{   /* Declarations of Guide1=Guide() SETTING parameters. */
char* reflect = mccGuide1_reflect;
MCNUM w1 = mccGuide1_w1;
MCNUM h1 = mccGuide1_h1;
MCNUM w2 = mccGuide1_w2;
MCNUM h2 = mccGuide1_h2;
MCNUM l = mccGuide1_l;
MCNUM R0 = mccGuide1_R0;
MCNUM Qc = mccGuide1_Qc;
MCNUM alpha = mccGuide1_alpha;
MCNUM m = mccGuide1_m;
MCNUM W = mccGuide1_W;
#line 93 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide.comp"
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
#line 18382 "./ESS_IN5_reprate.c"
}   /* End of Guide1=Guide() SETTING parameter declarations. */
#undef pTable
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompGuide1:
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

  /* TRACE Component Lmonslow1 [8] */
  mccoordschange(mcposrLmonslow1, mcrotrLmonslow1,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component Lmonslow1 (without coords transformations) */
  mcJumpTrace_Lmonslow1:
  SIG_MESSAGE("Lmonslow1 (Trace)");
  mcDEBUG_COMP("Lmonslow1")
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

#define mcabsorbComp mcabsorbCompLmonslow1
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
#define mccompcurname  Lmonslow1
#define mccompcurtype  L_monitor
#define mccompcurindex 8
#define nL mccLmonslow1_nL
#define L_N mccLmonslow1_L_N
#define L_p mccLmonslow1_L_p
#define L_p2 mccLmonslow1_L_p2
{   /* Declarations of Lmonslow1=L_monitor() SETTING parameters. */
char* filename = mccLmonslow1_filename;
MCNUM xmin = mccLmonslow1_xmin;
MCNUM xmax = mccLmonslow1_xmax;
MCNUM ymin = mccLmonslow1_ymin;
MCNUM ymax = mccLmonslow1_ymax;
MCNUM xwidth = mccLmonslow1_xwidth;
MCNUM yheight = mccLmonslow1_yheight;
MCNUM Lmin = mccLmonslow1_Lmin;
MCNUM Lmax = mccLmonslow1_Lmax;
MCNUM restore_neutron = mccLmonslow1_restore_neutron;
int nowritefile = mccLmonslow1_nowritefile;
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
#line 18526 "./ESS_IN5_reprate.c"
}   /* End of Lmonslow1=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompLmonslow1:
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

  /* TRACE Component PSDslow1 [9] */
  mccoordschange(mcposrPSDslow1, mcrotrPSDslow1,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component PSDslow1 (without coords transformations) */
  mcJumpTrace_PSDslow1:
  SIG_MESSAGE("PSDslow1 (Trace)");
  mcDEBUG_COMP("PSDslow1")
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

#define mcabsorbComp mcabsorbCompPSDslow1
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
#define mccompcurname  PSDslow1
#define mccompcurtype  PSD_monitor
#define mccompcurindex 9
#define PSD_N mccPSDslow1_PSD_N
#define PSD_p mccPSDslow1_PSD_p
#define PSD_p2 mccPSDslow1_PSD_p2
{   /* Declarations of PSDslow1=PSD_monitor() SETTING parameters. */
int nx = mccPSDslow1_nx;
int ny = mccPSDslow1_ny;
char* filename = mccPSDslow1_filename;
MCNUM xmin = mccPSDslow1_xmin;
MCNUM xmax = mccPSDslow1_xmax;
MCNUM ymin = mccPSDslow1_ymin;
MCNUM ymax = mccPSDslow1_ymax;
MCNUM xwidth = mccPSDslow1_xwidth;
MCNUM yheight = mccPSDslow1_yheight;
MCNUM restore_neutron = mccPSDslow1_restore_neutron;
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
#line 18664 "./ESS_IN5_reprate.c"
}   /* End of PSDslow1=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompPSDslow1:
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

  /* TRACE Component FOchop1 [10] */
  mccoordschange(mcposrFOchop1, mcrotrFOchop1,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component FOchop1 (without coords transformations) */
  mcJumpTrace_FOchop1:
  SIG_MESSAGE("FOchop1 (Trace)");
  mcDEBUG_COMP("FOchop1")
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

#define mcabsorbComp mcabsorbCompFOchop1
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
#define mccompcurname  FOchop1
#define mccompcurtype  DiskChopper
#define mccompcurindex 10
#define Tg mccFOchop1_Tg
#define To mccFOchop1_To
#define delta_y mccFOchop1_delta_y
#define height mccFOchop1_height
#define omega mccFOchop1_omega
{   /* Declarations of FOchop1=DiskChopper() SETTING parameters. */
MCNUM theta_0 = mccFOchop1_theta_0;
MCNUM radius = mccFOchop1_radius;
MCNUM yheight = mccFOchop1_yheight;
MCNUM nu = mccFOchop1_nu;
MCNUM nslit = mccFOchop1_nslit;
MCNUM jitter = mccFOchop1_jitter;
MCNUM delay = mccFOchop1_delay;
MCNUM isfirst = mccFOchop1_isfirst;
MCNUM n_pulse = mccFOchop1_n_pulse;
MCNUM abs_out = mccFOchop1_abs_out;
MCNUM phase = mccFOchop1_phase;
MCNUM xwidth = mccFOchop1_xwidth;
MCNUM verbose = mccFOchop1_verbose;
#line 133 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/DiskChopper.comp"
{
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

}
#line 18825 "./ESS_IN5_reprate.c"
}   /* End of FOchop1=DiskChopper() SETTING parameter declarations. */
#undef omega
#undef height
#undef delta_y
#undef To
#undef Tg
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompFOchop1:
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

  /* TRACE Component TOFLmon1 [11] */
  mccoordschange(mcposrTOFLmon1, mcrotrTOFLmon1,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component TOFLmon1 (without coords transformations) */
  mcJumpTrace_TOFLmon1:
  SIG_MESSAGE("TOFLmon1 (Trace)");
  mcDEBUG_COMP("TOFLmon1")
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

#define mcabsorbComp mcabsorbCompTOFLmon1
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
#define mccompcurname  TOFLmon1
#define mccompcurtype  TOFLambda_monitor
#define mccompcurindex 11
#define nL mccTOFLmon1_nL
#define nt mccTOFLmon1_nt
#define tmin mccTOFLmon1_tmin
#define tmax mccTOFLmon1_tmax
#define tt_0 mccTOFLmon1_tt_0
#define tt_1 mccTOFLmon1_tt_1
#define TOFL_N mccTOFLmon1_TOFL_N
#define TOFL_p mccTOFLmon1_TOFL_p
#define TOFL_p2 mccTOFLmon1_TOFL_p2
{   /* Declarations of TOFLmon1=TOFLambda_monitor() SETTING parameters. */
char* filename = mccTOFLmon1_filename;
MCNUM xmin = mccTOFLmon1_xmin;
MCNUM xmax = mccTOFLmon1_xmax;
MCNUM ymin = mccTOFLmon1_ymin;
MCNUM ymax = mccTOFLmon1_ymax;
MCNUM xwidth = mccTOFLmon1_xwidth;
MCNUM yheight = mccTOFLmon1_yheight;
MCNUM Lmin = mccTOFLmon1_Lmin;
MCNUM Lmax = mccTOFLmon1_Lmax;
MCNUM restore_neutron = mccTOFLmon1_restore_neutron;
int nowritefile = mccTOFLmon1_nowritefile;
#line 85 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/TOFLambda_monitor.comp"
{
    int i,j;
    double div;
    double lambda;

    PROP_Z0;
    lambda = (2*PI/V2K)/sqrt(vx*vx + vy*vy + vz*vz);
    if (x>xmin && x<xmax && y>ymin && y<ymax &&
        lambda > Lmin && lambda < Lmax)
    {
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
}
#line 18981 "./ESS_IN5_reprate.c"
}   /* End of TOFLmon1=TOFLambda_monitor() SETTING parameter declarations. */
#undef TOFL_p2
#undef TOFL_p
#undef TOFL_N
#undef tt_1
#undef tt_0
#undef tmax
#undef tmin
#undef nt
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompTOFLmon1:
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

  /* TRACE Component Lmon_afterslow1 [12] */
  mccoordschange(mcposrLmon_afterslow1, mcrotrLmon_afterslow1,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component Lmon_afterslow1 (without coords transformations) */
  mcJumpTrace_Lmon_afterslow1:
  SIG_MESSAGE("Lmon_afterslow1 (Trace)");
  mcDEBUG_COMP("Lmon_afterslow1")
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

#define mcabsorbComp mcabsorbCompLmon_afterslow1
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
#define mccompcurname  Lmon_afterslow1
#define mccompcurtype  L_monitor
#define mccompcurindex 12
#define nL mccLmon_afterslow1_nL
#define L_N mccLmon_afterslow1_L_N
#define L_p mccLmon_afterslow1_L_p
#define L_p2 mccLmon_afterslow1_L_p2
{   /* Declarations of Lmon_afterslow1=L_monitor() SETTING parameters. */
char* filename = mccLmon_afterslow1_filename;
MCNUM xmin = mccLmon_afterslow1_xmin;
MCNUM xmax = mccLmon_afterslow1_xmax;
MCNUM ymin = mccLmon_afterslow1_ymin;
MCNUM ymax = mccLmon_afterslow1_ymax;
MCNUM xwidth = mccLmon_afterslow1_xwidth;
MCNUM yheight = mccLmon_afterslow1_yheight;
MCNUM Lmin = mccLmon_afterslow1_Lmin;
MCNUM Lmax = mccLmon_afterslow1_Lmax;
MCNUM restore_neutron = mccLmon_afterslow1_restore_neutron;
int nowritefile = mccLmon_afterslow1_nowritefile;
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
#line 19133 "./ESS_IN5_reprate.c"
}   /* End of Lmon_afterslow1=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompLmon_afterslow1:
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

  /* TRACE Component PSD_afterslow1 [13] */
  mccoordschange(mcposrPSD_afterslow1, mcrotrPSD_afterslow1,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component PSD_afterslow1 (without coords transformations) */
  mcJumpTrace_PSD_afterslow1:
  SIG_MESSAGE("PSD_afterslow1 (Trace)");
  mcDEBUG_COMP("PSD_afterslow1")
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

#define mcabsorbComp mcabsorbCompPSD_afterslow1
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
#define mccompcurname  PSD_afterslow1
#define mccompcurtype  PSD_monitor
#define mccompcurindex 13
#define PSD_N mccPSD_afterslow1_PSD_N
#define PSD_p mccPSD_afterslow1_PSD_p
#define PSD_p2 mccPSD_afterslow1_PSD_p2
{   /* Declarations of PSD_afterslow1=PSD_monitor() SETTING parameters. */
int nx = mccPSD_afterslow1_nx;
int ny = mccPSD_afterslow1_ny;
char* filename = mccPSD_afterslow1_filename;
MCNUM xmin = mccPSD_afterslow1_xmin;
MCNUM xmax = mccPSD_afterslow1_xmax;
MCNUM ymin = mccPSD_afterslow1_ymin;
MCNUM ymax = mccPSD_afterslow1_ymax;
MCNUM xwidth = mccPSD_afterslow1_xwidth;
MCNUM yheight = mccPSD_afterslow1_yheight;
MCNUM restore_neutron = mccPSD_afterslow1_restore_neutron;
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
#line 19271 "./ESS_IN5_reprate.c"
}   /* End of PSD_afterslow1=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompPSD_afterslow1:
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

  /* TRACE Component Guidelong1 [14] */
  mccoordschange(mcposrGuidelong1, mcrotrGuidelong1,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component Guidelong1 (without coords transformations) */
  mcJumpTrace_Guidelong1:
  SIG_MESSAGE("Guidelong1 (Trace)");
  mcDEBUG_COMP("Guidelong1")
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

#define mcabsorbComp mcabsorbCompGuidelong1
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
#define mccompcurname  Guidelong1
#define mccompcurtype  Guide
#define mccompcurindex 14
#define pTable mccGuidelong1_pTable
{   /* Declarations of Guidelong1=Guide() SETTING parameters. */
char* reflect = mccGuidelong1_reflect;
MCNUM w1 = mccGuidelong1_w1;
MCNUM h1 = mccGuidelong1_h1;
MCNUM w2 = mccGuidelong1_w2;
MCNUM h2 = mccGuidelong1_h2;
MCNUM l = mccGuidelong1_l;
MCNUM R0 = mccGuidelong1_R0;
MCNUM Qc = mccGuidelong1_Qc;
MCNUM alpha = mccGuidelong1_alpha;
MCNUM m = mccGuidelong1_m;
MCNUM W = mccGuidelong1_W;
#line 93 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide.comp"
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
#line 19499 "./ESS_IN5_reprate.c"
}   /* End of Guidelong1=Guide() SETTING parameter declarations. */
#undef pTable
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompGuidelong1:
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

  /* TRACE Component Guidelong1b [15] */
  mccoordschange(mcposrGuidelong1b, mcrotrGuidelong1b,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component Guidelong1b (without coords transformations) */
  mcJumpTrace_Guidelong1b:
  SIG_MESSAGE("Guidelong1b (Trace)");
  mcDEBUG_COMP("Guidelong1b")
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

#define mcabsorbComp mcabsorbCompGuidelong1b
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
#define mccompcurname  Guidelong1b
#define mccompcurtype  Guide
#define mccompcurindex 15
#define pTable mccGuidelong1b_pTable
{   /* Declarations of Guidelong1b=Guide() SETTING parameters. */
char* reflect = mccGuidelong1b_reflect;
MCNUM w1 = mccGuidelong1b_w1;
MCNUM h1 = mccGuidelong1b_h1;
MCNUM w2 = mccGuidelong1b_w2;
MCNUM h2 = mccGuidelong1b_h2;
MCNUM l = mccGuidelong1b_l;
MCNUM R0 = mccGuidelong1b_R0;
MCNUM Qc = mccGuidelong1b_Qc;
MCNUM alpha = mccGuidelong1b_alpha;
MCNUM m = mccGuidelong1b_m;
MCNUM W = mccGuidelong1b_W;
#line 93 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide.comp"
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
#line 19725 "./ESS_IN5_reprate.c"
}   /* End of Guidelong1b=Guide() SETTING parameter declarations. */
#undef pTable
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompGuidelong1b:
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

  /* TRACE Component Lmon_slow2 [16] */
  mccoordschange(mcposrLmon_slow2, mcrotrLmon_slow2,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component Lmon_slow2 (without coords transformations) */
  mcJumpTrace_Lmon_slow2:
  SIG_MESSAGE("Lmon_slow2 (Trace)");
  mcDEBUG_COMP("Lmon_slow2")
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

#define mcabsorbComp mcabsorbCompLmon_slow2
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
#define mccompcurname  Lmon_slow2
#define mccompcurtype  L_monitor
#define mccompcurindex 16
#define nL mccLmon_slow2_nL
#define L_N mccLmon_slow2_L_N
#define L_p mccLmon_slow2_L_p
#define L_p2 mccLmon_slow2_L_p2
{   /* Declarations of Lmon_slow2=L_monitor() SETTING parameters. */
char* filename = mccLmon_slow2_filename;
MCNUM xmin = mccLmon_slow2_xmin;
MCNUM xmax = mccLmon_slow2_xmax;
MCNUM ymin = mccLmon_slow2_ymin;
MCNUM ymax = mccLmon_slow2_ymax;
MCNUM xwidth = mccLmon_slow2_xwidth;
MCNUM yheight = mccLmon_slow2_yheight;
MCNUM Lmin = mccLmon_slow2_Lmin;
MCNUM Lmax = mccLmon_slow2_Lmax;
MCNUM restore_neutron = mccLmon_slow2_restore_neutron;
int nowritefile = mccLmon_slow2_nowritefile;
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
#line 19869 "./ESS_IN5_reprate.c"
}   /* End of Lmon_slow2=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompLmon_slow2:
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

  /* TRACE Component FOchop2 [17] */
  mccoordschange(mcposrFOchop2, mcrotrFOchop2,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component FOchop2 (without coords transformations) */
  mcJumpTrace_FOchop2:
  SIG_MESSAGE("FOchop2 (Trace)");
  mcDEBUG_COMP("FOchop2")
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

#define mcabsorbComp mcabsorbCompFOchop2
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
#define mccompcurname  FOchop2
#define mccompcurtype  DiskChopper
#define mccompcurindex 17
#define Tg mccFOchop2_Tg
#define To mccFOchop2_To
#define delta_y mccFOchop2_delta_y
#define height mccFOchop2_height
#define omega mccFOchop2_omega
{   /* Declarations of FOchop2=DiskChopper() SETTING parameters. */
MCNUM theta_0 = mccFOchop2_theta_0;
MCNUM radius = mccFOchop2_radius;
MCNUM yheight = mccFOchop2_yheight;
MCNUM nu = mccFOchop2_nu;
MCNUM nslit = mccFOchop2_nslit;
MCNUM jitter = mccFOchop2_jitter;
MCNUM delay = mccFOchop2_delay;
MCNUM isfirst = mccFOchop2_isfirst;
MCNUM n_pulse = mccFOchop2_n_pulse;
MCNUM abs_out = mccFOchop2_abs_out;
MCNUM phase = mccFOchop2_phase;
MCNUM xwidth = mccFOchop2_xwidth;
MCNUM verbose = mccFOchop2_verbose;
#line 133 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/DiskChopper.comp"
{
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

}
#line 20031 "./ESS_IN5_reprate.c"
}   /* End of FOchop2=DiskChopper() SETTING parameter declarations. */
#undef omega
#undef height
#undef delta_y
#undef To
#undef Tg
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompFOchop2:
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

  /* TRACE Component Fastchop1 [18] */
  mccoordschange(mcposrFastchop1, mcrotrFastchop1,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component Fastchop1 (without coords transformations) */
  mcJumpTrace_Fastchop1:
  SIG_MESSAGE("Fastchop1 (Trace)");
  mcDEBUG_COMP("Fastchop1")
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

#define mcabsorbComp mcabsorbCompFastchop1
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
#define mccompcurname  Fastchop1
#define mccompcurtype  DiskChopper
#define mccompcurindex 18
#define Tg mccFastchop1_Tg
#define To mccFastchop1_To
#define delta_y mccFastchop1_delta_y
#define height mccFastchop1_height
#define omega mccFastchop1_omega
{   /* Declarations of Fastchop1=DiskChopper() SETTING parameters. */
MCNUM theta_0 = mccFastchop1_theta_0;
MCNUM radius = mccFastchop1_radius;
MCNUM yheight = mccFastchop1_yheight;
MCNUM nu = mccFastchop1_nu;
MCNUM nslit = mccFastchop1_nslit;
MCNUM jitter = mccFastchop1_jitter;
MCNUM delay = mccFastchop1_delay;
MCNUM isfirst = mccFastchop1_isfirst;
MCNUM n_pulse = mccFastchop1_n_pulse;
MCNUM abs_out = mccFastchop1_abs_out;
MCNUM phase = mccFastchop1_phase;
MCNUM xwidth = mccFastchop1_xwidth;
MCNUM verbose = mccFastchop1_verbose;
#line 133 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/DiskChopper.comp"
{
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

}
#line 20194 "./ESS_IN5_reprate.c"
}   /* End of Fastchop1=DiskChopper() SETTING parameter declarations. */
#undef omega
#undef height
#undef delta_y
#undef To
#undef Tg
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompFastchop1:
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

  /* TRACE Component PSD_afterslow2 [19] */
  mccoordschange(mcposrPSD_afterslow2, mcrotrPSD_afterslow2,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component PSD_afterslow2 (without coords transformations) */
  mcJumpTrace_PSD_afterslow2:
  SIG_MESSAGE("PSD_afterslow2 (Trace)");
  mcDEBUG_COMP("PSD_afterslow2")
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

#define mcabsorbComp mcabsorbCompPSD_afterslow2
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
#define mccompcurname  PSD_afterslow2
#define mccompcurtype  PSD_monitor
#define mccompcurindex 19
#define PSD_N mccPSD_afterslow2_PSD_N
#define PSD_p mccPSD_afterslow2_PSD_p
#define PSD_p2 mccPSD_afterslow2_PSD_p2
{   /* Declarations of PSD_afterslow2=PSD_monitor() SETTING parameters. */
int nx = mccPSD_afterslow2_nx;
int ny = mccPSD_afterslow2_ny;
char* filename = mccPSD_afterslow2_filename;
MCNUM xmin = mccPSD_afterslow2_xmin;
MCNUM xmax = mccPSD_afterslow2_xmax;
MCNUM ymin = mccPSD_afterslow2_ymin;
MCNUM ymax = mccPSD_afterslow2_ymax;
MCNUM xwidth = mccPSD_afterslow2_xwidth;
MCNUM yheight = mccPSD_afterslow2_yheight;
MCNUM restore_neutron = mccPSD_afterslow2_restore_neutron;
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
#line 20333 "./ESS_IN5_reprate.c"
}   /* End of PSD_afterslow2=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompPSD_afterslow2:
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

  /* TRACE Component Lmon_afterslow2 [20] */
  mccoordschange(mcposrLmon_afterslow2, mcrotrLmon_afterslow2,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component Lmon_afterslow2 (without coords transformations) */
  mcJumpTrace_Lmon_afterslow2:
  SIG_MESSAGE("Lmon_afterslow2 (Trace)");
  mcDEBUG_COMP("Lmon_afterslow2")
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

#define mcabsorbComp mcabsorbCompLmon_afterslow2
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
#define mccompcurname  Lmon_afterslow2
#define mccompcurtype  L_monitor
#define mccompcurindex 20
#define nL mccLmon_afterslow2_nL
#define L_N mccLmon_afterslow2_L_N
#define L_p mccLmon_afterslow2_L_p
#define L_p2 mccLmon_afterslow2_L_p2
{   /* Declarations of Lmon_afterslow2=L_monitor() SETTING parameters. */
char* filename = mccLmon_afterslow2_filename;
MCNUM xmin = mccLmon_afterslow2_xmin;
MCNUM xmax = mccLmon_afterslow2_xmax;
MCNUM ymin = mccLmon_afterslow2_ymin;
MCNUM ymax = mccLmon_afterslow2_ymax;
MCNUM xwidth = mccLmon_afterslow2_xwidth;
MCNUM yheight = mccLmon_afterslow2_yheight;
MCNUM Lmin = mccLmon_afterslow2_Lmin;
MCNUM Lmax = mccLmon_afterslow2_Lmax;
MCNUM restore_neutron = mccLmon_afterslow2_restore_neutron;
int nowritefile = mccLmon_afterslow2_nowritefile;
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
#line 20479 "./ESS_IN5_reprate.c"
}   /* End of Lmon_afterslow2=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompLmon_afterslow2:
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

  /* TRACE Component TOFL_afterslow2 [21] */
  mccoordschange(mcposrTOFL_afterslow2, mcrotrTOFL_afterslow2,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component TOFL_afterslow2 (without coords transformations) */
  mcJumpTrace_TOFL_afterslow2:
  SIG_MESSAGE("TOFL_afterslow2 (Trace)");
  mcDEBUG_COMP("TOFL_afterslow2")
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

#define mcabsorbComp mcabsorbCompTOFL_afterslow2
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
#define mccompcurname  TOFL_afterslow2
#define mccompcurtype  TOFLambda_monitor
#define mccompcurindex 21
#define nL mccTOFL_afterslow2_nL
#define nt mccTOFL_afterslow2_nt
#define tmin mccTOFL_afterslow2_tmin
#define tmax mccTOFL_afterslow2_tmax
#define tt_0 mccTOFL_afterslow2_tt_0
#define tt_1 mccTOFL_afterslow2_tt_1
#define TOFL_N mccTOFL_afterslow2_TOFL_N
#define TOFL_p mccTOFL_afterslow2_TOFL_p
#define TOFL_p2 mccTOFL_afterslow2_TOFL_p2
{   /* Declarations of TOFL_afterslow2=TOFLambda_monitor() SETTING parameters. */
char* filename = mccTOFL_afterslow2_filename;
MCNUM xmin = mccTOFL_afterslow2_xmin;
MCNUM xmax = mccTOFL_afterslow2_xmax;
MCNUM ymin = mccTOFL_afterslow2_ymin;
MCNUM ymax = mccTOFL_afterslow2_ymax;
MCNUM xwidth = mccTOFL_afterslow2_xwidth;
MCNUM yheight = mccTOFL_afterslow2_yheight;
MCNUM Lmin = mccTOFL_afterslow2_Lmin;
MCNUM Lmax = mccTOFL_afterslow2_Lmax;
MCNUM restore_neutron = mccTOFL_afterslow2_restore_neutron;
int nowritefile = mccTOFL_afterslow2_nowritefile;
#line 85 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/TOFLambda_monitor.comp"
{
    int i,j;
    double div;
    double lambda;

    PROP_Z0;
    lambda = (2*PI/V2K)/sqrt(vx*vx + vy*vy + vz*vz);
    if (x>xmin && x<xmax && y>ymin && y<ymax &&
        lambda > Lmin && lambda < Lmax)
    {
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
}
#line 20634 "./ESS_IN5_reprate.c"
}   /* End of TOFL_afterslow2=TOFLambda_monitor() SETTING parameter declarations. */
#undef TOFL_p2
#undef TOFL_p
#undef TOFL_N
#undef tt_1
#undef tt_0
#undef tmax
#undef tmin
#undef nt
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompTOFL_afterslow2:
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

  /* TRACE Component Guidelong2 [22] */
  mccoordschange(mcposrGuidelong2, mcrotrGuidelong2,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component Guidelong2 (without coords transformations) */
  mcJumpTrace_Guidelong2:
  SIG_MESSAGE("Guidelong2 (Trace)");
  mcDEBUG_COMP("Guidelong2")
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

#define mcabsorbComp mcabsorbCompGuidelong2
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
#define mccompcurname  Guidelong2
#define mccompcurtype  Guide
#define mccompcurindex 22
#define pTable mccGuidelong2_pTable
{   /* Declarations of Guidelong2=Guide() SETTING parameters. */
char* reflect = mccGuidelong2_reflect;
MCNUM w1 = mccGuidelong2_w1;
MCNUM h1 = mccGuidelong2_h1;
MCNUM w2 = mccGuidelong2_w2;
MCNUM h2 = mccGuidelong2_h2;
MCNUM l = mccGuidelong2_l;
MCNUM R0 = mccGuidelong2_R0;
MCNUM Qc = mccGuidelong2_Qc;
MCNUM alpha = mccGuidelong2_alpha;
MCNUM m = mccGuidelong2_m;
MCNUM W = mccGuidelong2_W;
#line 93 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide.comp"
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
#line 20868 "./ESS_IN5_reprate.c"
}   /* End of Guidelong2=Guide() SETTING parameter declarations. */
#undef pTable
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompGuidelong2:
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

  /* TRACE Component Lmon_beforeballistic [23] */
  mccoordschange(mcposrLmon_beforeballistic, mcrotrLmon_beforeballistic,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component Lmon_beforeballistic (without coords transformations) */
  mcJumpTrace_Lmon_beforeballistic:
  SIG_MESSAGE("Lmon_beforeballistic (Trace)");
  mcDEBUG_COMP("Lmon_beforeballistic")
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

#define mcabsorbComp mcabsorbCompLmon_beforeballistic
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
#define mccompcurname  Lmon_beforeballistic
#define mccompcurtype  L_monitor
#define mccompcurindex 23
#define nL mccLmon_beforeballistic_nL
#define L_N mccLmon_beforeballistic_L_N
#define L_p mccLmon_beforeballistic_L_p
#define L_p2 mccLmon_beforeballistic_L_p2
{   /* Declarations of Lmon_beforeballistic=L_monitor() SETTING parameters. */
char* filename = mccLmon_beforeballistic_filename;
MCNUM xmin = mccLmon_beforeballistic_xmin;
MCNUM xmax = mccLmon_beforeballistic_xmax;
MCNUM ymin = mccLmon_beforeballistic_ymin;
MCNUM ymax = mccLmon_beforeballistic_ymax;
MCNUM xwidth = mccLmon_beforeballistic_xwidth;
MCNUM yheight = mccLmon_beforeballistic_yheight;
MCNUM Lmin = mccLmon_beforeballistic_Lmin;
MCNUM Lmax = mccLmon_beforeballistic_Lmax;
MCNUM restore_neutron = mccLmon_beforeballistic_restore_neutron;
int nowritefile = mccLmon_beforeballistic_nowritefile;
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
#line 21012 "./ESS_IN5_reprate.c"
}   /* End of Lmon_beforeballistic=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompLmon_beforeballistic:
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

  /* TRACE Component PSD_beforeballistic [24] */
  mccoordschange(mcposrPSD_beforeballistic, mcrotrPSD_beforeballistic,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component PSD_beforeballistic (without coords transformations) */
  mcJumpTrace_PSD_beforeballistic:
  SIG_MESSAGE("PSD_beforeballistic (Trace)");
  mcDEBUG_COMP("PSD_beforeballistic")
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

#define mcabsorbComp mcabsorbCompPSD_beforeballistic
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
#define mccompcurname  PSD_beforeballistic
#define mccompcurtype  PSD_monitor
#define mccompcurindex 24
#define PSD_N mccPSD_beforeballistic_PSD_N
#define PSD_p mccPSD_beforeballistic_PSD_p
#define PSD_p2 mccPSD_beforeballistic_PSD_p2
{   /* Declarations of PSD_beforeballistic=PSD_monitor() SETTING parameters. */
int nx = mccPSD_beforeballistic_nx;
int ny = mccPSD_beforeballistic_ny;
char* filename = mccPSD_beforeballistic_filename;
MCNUM xmin = mccPSD_beforeballistic_xmin;
MCNUM xmax = mccPSD_beforeballistic_xmax;
MCNUM ymin = mccPSD_beforeballistic_ymin;
MCNUM ymax = mccPSD_beforeballistic_ymax;
MCNUM xwidth = mccPSD_beforeballistic_xwidth;
MCNUM yheight = mccPSD_beforeballistic_yheight;
MCNUM restore_neutron = mccPSD_beforeballistic_restore_neutron;
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
#line 21150 "./ESS_IN5_reprate.c"
}   /* End of PSD_beforeballistic=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompPSD_beforeballistic:
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

  /* TRACE Component Guidelong2a [25] */
  mccoordschange(mcposrGuidelong2a, mcrotrGuidelong2a,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component Guidelong2a (without coords transformations) */
  mcJumpTrace_Guidelong2a:
  SIG_MESSAGE("Guidelong2a (Trace)");
  mcDEBUG_COMP("Guidelong2a")
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

#define mcabsorbComp mcabsorbCompGuidelong2a
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
#define mccompcurname  Guidelong2a
#define mccompcurtype  Guide
#define mccompcurindex 25
#define pTable mccGuidelong2a_pTable
{   /* Declarations of Guidelong2a=Guide() SETTING parameters. */
char* reflect = mccGuidelong2a_reflect;
MCNUM w1 = mccGuidelong2a_w1;
MCNUM h1 = mccGuidelong2a_h1;
MCNUM w2 = mccGuidelong2a_w2;
MCNUM h2 = mccGuidelong2a_h2;
MCNUM l = mccGuidelong2a_l;
MCNUM R0 = mccGuidelong2a_R0;
MCNUM Qc = mccGuidelong2a_Qc;
MCNUM alpha = mccGuidelong2a_alpha;
MCNUM m = mccGuidelong2a_m;
MCNUM W = mccGuidelong2a_W;
#line 93 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide.comp"
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
#line 21378 "./ESS_IN5_reprate.c"
}   /* End of Guidelong2a=Guide() SETTING parameter declarations. */
#undef pTable
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompGuidelong2a:
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

  /* TRACE Component Lmonfast2 [26] */
  mccoordschange(mcposrLmonfast2, mcrotrLmonfast2,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component Lmonfast2 (without coords transformations) */
  mcJumpTrace_Lmonfast2:
  SIG_MESSAGE("Lmonfast2 (Trace)");
  mcDEBUG_COMP("Lmonfast2")
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

#define mcabsorbComp mcabsorbCompLmonfast2
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
#define mccompcurname  Lmonfast2
#define mccompcurtype  L_monitor
#define mccompcurindex 26
#define nL mccLmonfast2_nL
#define L_N mccLmonfast2_L_N
#define L_p mccLmonfast2_L_p
#define L_p2 mccLmonfast2_L_p2
{   /* Declarations of Lmonfast2=L_monitor() SETTING parameters. */
char* filename = mccLmonfast2_filename;
MCNUM xmin = mccLmonfast2_xmin;
MCNUM xmax = mccLmonfast2_xmax;
MCNUM ymin = mccLmonfast2_ymin;
MCNUM ymax = mccLmonfast2_ymax;
MCNUM xwidth = mccLmonfast2_xwidth;
MCNUM yheight = mccLmonfast2_yheight;
MCNUM Lmin = mccLmonfast2_Lmin;
MCNUM Lmax = mccLmonfast2_Lmax;
MCNUM restore_neutron = mccLmonfast2_restore_neutron;
int nowritefile = mccLmonfast2_nowritefile;
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
#line 21522 "./ESS_IN5_reprate.c"
}   /* End of Lmonfast2=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompLmonfast2:
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

  /* TRACE Component Lmonfast2_zoom [27] */
  mccoordschange(mcposrLmonfast2_zoom, mcrotrLmonfast2_zoom,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component Lmonfast2_zoom (without coords transformations) */
  mcJumpTrace_Lmonfast2_zoom:
  SIG_MESSAGE("Lmonfast2_zoom (Trace)");
  mcDEBUG_COMP("Lmonfast2_zoom")
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

#define mcabsorbComp mcabsorbCompLmonfast2_zoom
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
#define mccompcurname  Lmonfast2_zoom
#define mccompcurtype  L_monitor
#define mccompcurindex 27
#define nL mccLmonfast2_zoom_nL
#define L_N mccLmonfast2_zoom_L_N
#define L_p mccLmonfast2_zoom_L_p
#define L_p2 mccLmonfast2_zoom_L_p2
{   /* Declarations of Lmonfast2_zoom=L_monitor() SETTING parameters. */
char* filename = mccLmonfast2_zoom_filename;
MCNUM xmin = mccLmonfast2_zoom_xmin;
MCNUM xmax = mccLmonfast2_zoom_xmax;
MCNUM ymin = mccLmonfast2_zoom_ymin;
MCNUM ymax = mccLmonfast2_zoom_ymax;
MCNUM xwidth = mccLmonfast2_zoom_xwidth;
MCNUM yheight = mccLmonfast2_zoom_yheight;
MCNUM Lmin = mccLmonfast2_zoom_Lmin;
MCNUM Lmax = mccLmonfast2_zoom_Lmax;
MCNUM restore_neutron = mccLmonfast2_zoom_restore_neutron;
int nowritefile = mccLmonfast2_zoom_nowritefile;
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
#line 21669 "./ESS_IN5_reprate.c"
}   /* End of Lmonfast2_zoom=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompLmonfast2_zoom:
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

  /* TRACE Component TOFLfast2 [28] */
  mccoordschange(mcposrTOFLfast2, mcrotrTOFLfast2,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component TOFLfast2 (without coords transformations) */
  mcJumpTrace_TOFLfast2:
  SIG_MESSAGE("TOFLfast2 (Trace)");
  mcDEBUG_COMP("TOFLfast2")
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

#define mcabsorbComp mcabsorbCompTOFLfast2
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
#define mccompcurname  TOFLfast2
#define mccompcurtype  TOFLambda_monitor
#define mccompcurindex 28
#define nL mccTOFLfast2_nL
#define nt mccTOFLfast2_nt
#define tmin mccTOFLfast2_tmin
#define tmax mccTOFLfast2_tmax
#define tt_0 mccTOFLfast2_tt_0
#define tt_1 mccTOFLfast2_tt_1
#define TOFL_N mccTOFLfast2_TOFL_N
#define TOFL_p mccTOFLfast2_TOFL_p
#define TOFL_p2 mccTOFLfast2_TOFL_p2
{   /* Declarations of TOFLfast2=TOFLambda_monitor() SETTING parameters. */
char* filename = mccTOFLfast2_filename;
MCNUM xmin = mccTOFLfast2_xmin;
MCNUM xmax = mccTOFLfast2_xmax;
MCNUM ymin = mccTOFLfast2_ymin;
MCNUM ymax = mccTOFLfast2_ymax;
MCNUM xwidth = mccTOFLfast2_xwidth;
MCNUM yheight = mccTOFLfast2_yheight;
MCNUM Lmin = mccTOFLfast2_Lmin;
MCNUM Lmax = mccTOFLfast2_Lmax;
MCNUM restore_neutron = mccTOFLfast2_restore_neutron;
int nowritefile = mccTOFLfast2_nowritefile;
#line 85 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/TOFLambda_monitor.comp"
{
    int i,j;
    double div;
    double lambda;

    PROP_Z0;
    lambda = (2*PI/V2K)/sqrt(vx*vx + vy*vy + vz*vz);
    if (x>xmin && x<xmax && y>ymin && y<ymax &&
        lambda > Lmin && lambda < Lmax)
    {
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
}
#line 21824 "./ESS_IN5_reprate.c"
}   /* End of TOFLfast2=TOFLambda_monitor() SETTING parameter declarations. */
#undef TOFL_p2
#undef TOFL_p
#undef TOFL_N
#undef tt_1
#undef tt_0
#undef tmax
#undef tmin
#undef nt
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompTOFLfast2:
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

  /* TRACE Component TOFLfast2zoom [29] */
  mccoordschange(mcposrTOFLfast2zoom, mcrotrTOFLfast2zoom,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component TOFLfast2zoom (without coords transformations) */
  mcJumpTrace_TOFLfast2zoom:
  SIG_MESSAGE("TOFLfast2zoom (Trace)");
  mcDEBUG_COMP("TOFLfast2zoom")
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

#define mcabsorbComp mcabsorbCompTOFLfast2zoom
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
#define mccompcurname  TOFLfast2zoom
#define mccompcurtype  TOFLambda_monitor
#define mccompcurindex 29
#define nL mccTOFLfast2zoom_nL
#define nt mccTOFLfast2zoom_nt
#define tmin mccTOFLfast2zoom_tmin
#define tmax mccTOFLfast2zoom_tmax
#define tt_0 mccTOFLfast2zoom_tt_0
#define tt_1 mccTOFLfast2zoom_tt_1
#define TOFL_N mccTOFLfast2zoom_TOFL_N
#define TOFL_p mccTOFLfast2zoom_TOFL_p
#define TOFL_p2 mccTOFLfast2zoom_TOFL_p2
{   /* Declarations of TOFLfast2zoom=TOFLambda_monitor() SETTING parameters. */
char* filename = mccTOFLfast2zoom_filename;
MCNUM xmin = mccTOFLfast2zoom_xmin;
MCNUM xmax = mccTOFLfast2zoom_xmax;
MCNUM ymin = mccTOFLfast2zoom_ymin;
MCNUM ymax = mccTOFLfast2zoom_ymax;
MCNUM xwidth = mccTOFLfast2zoom_xwidth;
MCNUM yheight = mccTOFLfast2zoom_yheight;
MCNUM Lmin = mccTOFLfast2zoom_Lmin;
MCNUM Lmax = mccTOFLfast2zoom_Lmax;
MCNUM restore_neutron = mccTOFLfast2zoom_restore_neutron;
int nowritefile = mccTOFLfast2zoom_nowritefile;
#line 85 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/TOFLambda_monitor.comp"
{
    int i,j;
    double div;
    double lambda;

    PROP_Z0;
    lambda = (2*PI/V2K)/sqrt(vx*vx + vy*vy + vz*vz);
    if (x>xmin && x<xmax && y>ymin && y<ymax &&
        lambda > Lmin && lambda < Lmax)
    {
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
}
#line 21984 "./ESS_IN5_reprate.c"
}   /* End of TOFLfast2zoom=TOFLambda_monitor() SETTING parameter declarations. */
#undef TOFL_p2
#undef TOFL_p
#undef TOFL_N
#undef tt_1
#undef tt_0
#undef tmax
#undef tmin
#undef nt
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompTOFLfast2zoom:
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

  /* TRACE Component PSDfast2 [30] */
  mccoordschange(mcposrPSDfast2, mcrotrPSDfast2,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component PSDfast2 (without coords transformations) */
  mcJumpTrace_PSDfast2:
  SIG_MESSAGE("PSDfast2 (Trace)");
  mcDEBUG_COMP("PSDfast2")
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

#define mcabsorbComp mcabsorbCompPSDfast2
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
#define mccompcurname  PSDfast2
#define mccompcurtype  PSD_monitor
#define mccompcurindex 30
#define PSD_N mccPSDfast2_PSD_N
#define PSD_p mccPSDfast2_PSD_p
#define PSD_p2 mccPSDfast2_PSD_p2
{   /* Declarations of PSDfast2=PSD_monitor() SETTING parameters. */
int nx = mccPSDfast2_nx;
int ny = mccPSDfast2_ny;
char* filename = mccPSDfast2_filename;
MCNUM xmin = mccPSDfast2_xmin;
MCNUM xmax = mccPSDfast2_xmax;
MCNUM ymin = mccPSDfast2_ymin;
MCNUM ymax = mccPSDfast2_ymax;
MCNUM xwidth = mccPSDfast2_xwidth;
MCNUM yheight = mccPSDfast2_yheight;
MCNUM restore_neutron = mccPSDfast2_restore_neutron;
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
#line 22127 "./ESS_IN5_reprate.c"
}   /* End of PSDfast2=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompPSDfast2:
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

  /* TRACE Component Fastchop2 [31] */
  mccoordschange(mcposrFastchop2, mcrotrFastchop2,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component Fastchop2 (without coords transformations) */
  mcJumpTrace_Fastchop2:
  SIG_MESSAGE("Fastchop2 (Trace)");
  mcDEBUG_COMP("Fastchop2")
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

#define mcabsorbComp mcabsorbCompFastchop2
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
#define mccompcurname  Fastchop2
#define mccompcurtype  DiskChopper
#define mccompcurindex 31
#define Tg mccFastchop2_Tg
#define To mccFastchop2_To
#define delta_y mccFastchop2_delta_y
#define height mccFastchop2_height
#define omega mccFastchop2_omega
{   /* Declarations of Fastchop2=DiskChopper() SETTING parameters. */
MCNUM theta_0 = mccFastchop2_theta_0;
MCNUM radius = mccFastchop2_radius;
MCNUM yheight = mccFastchop2_yheight;
MCNUM nu = mccFastchop2_nu;
MCNUM nslit = mccFastchop2_nslit;
MCNUM jitter = mccFastchop2_jitter;
MCNUM delay = mccFastchop2_delay;
MCNUM isfirst = mccFastchop2_isfirst;
MCNUM n_pulse = mccFastchop2_n_pulse;
MCNUM abs_out = mccFastchop2_abs_out;
MCNUM phase = mccFastchop2_phase;
MCNUM xwidth = mccFastchop2_xwidth;
MCNUM verbose = mccFastchop2_verbose;
#line 133 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/DiskChopper.comp"
{
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

}
#line 22288 "./ESS_IN5_reprate.c"
}   /* End of Fastchop2=DiskChopper() SETTING parameter declarations. */
#undef omega
#undef height
#undef delta_y
#undef To
#undef Tg
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompFastchop2:
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

  /* TRACE Component Fastchop2counter [32] */
  mccoordschange(mcposrFastchop2counter, mcrotrFastchop2counter,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component Fastchop2counter (without coords transformations) */
  mcJumpTrace_Fastchop2counter:
  SIG_MESSAGE("Fastchop2counter (Trace)");
  mcDEBUG_COMP("Fastchop2counter")
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

#define mcabsorbComp mcabsorbCompFastchop2counter
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
#define mccompcurname  Fastchop2counter
#define mccompcurtype  DiskChopper
#define mccompcurindex 32
#define Tg mccFastchop2counter_Tg
#define To mccFastchop2counter_To
#define delta_y mccFastchop2counter_delta_y
#define height mccFastchop2counter_height
#define omega mccFastchop2counter_omega
{   /* Declarations of Fastchop2counter=DiskChopper() SETTING parameters. */
MCNUM theta_0 = mccFastchop2counter_theta_0;
MCNUM radius = mccFastchop2counter_radius;
MCNUM yheight = mccFastchop2counter_yheight;
MCNUM nu = mccFastchop2counter_nu;
MCNUM nslit = mccFastchop2counter_nslit;
MCNUM jitter = mccFastchop2counter_jitter;
MCNUM delay = mccFastchop2counter_delay;
MCNUM isfirst = mccFastchop2counter_isfirst;
MCNUM n_pulse = mccFastchop2counter_n_pulse;
MCNUM abs_out = mccFastchop2counter_abs_out;
MCNUM phase = mccFastchop2counter_phase;
MCNUM xwidth = mccFastchop2counter_xwidth;
MCNUM verbose = mccFastchop2counter_verbose;
#line 133 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/DiskChopper.comp"
{
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

}
#line 22451 "./ESS_IN5_reprate.c"
}   /* End of Fastchop2counter=DiskChopper() SETTING parameter declarations. */
#undef omega
#undef height
#undef delta_y
#undef To
#undef Tg
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompFastchop2counter:
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

  /* TRACE Component FOchop3 [33] */
  mccoordschange(mcposrFOchop3, mcrotrFOchop3,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component FOchop3 (without coords transformations) */
  mcJumpTrace_FOchop3:
  SIG_MESSAGE("FOchop3 (Trace)");
  mcDEBUG_COMP("FOchop3")
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

#define mcabsorbComp mcabsorbCompFOchop3
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
#define mccompcurname  FOchop3
#define mccompcurtype  DiskChopper
#define mccompcurindex 33
#define Tg mccFOchop3_Tg
#define To mccFOchop3_To
#define delta_y mccFOchop3_delta_y
#define height mccFOchop3_height
#define omega mccFOchop3_omega
{   /* Declarations of FOchop3=DiskChopper() SETTING parameters. */
MCNUM theta_0 = mccFOchop3_theta_0;
MCNUM radius = mccFOchop3_radius;
MCNUM yheight = mccFOchop3_yheight;
MCNUM nu = mccFOchop3_nu;
MCNUM nslit = mccFOchop3_nslit;
MCNUM jitter = mccFOchop3_jitter;
MCNUM delay = mccFOchop3_delay;
MCNUM isfirst = mccFOchop3_isfirst;
MCNUM n_pulse = mccFOchop3_n_pulse;
MCNUM abs_out = mccFOchop3_abs_out;
MCNUM phase = mccFOchop3_phase;
MCNUM xwidth = mccFOchop3_xwidth;
MCNUM verbose = mccFOchop3_verbose;
#line 133 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/DiskChopper.comp"
{
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

}
#line 22614 "./ESS_IN5_reprate.c"
}   /* End of FOchop3=DiskChopper() SETTING parameter declarations. */
#undef omega
#undef height
#undef delta_y
#undef To
#undef Tg
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompFOchop3:
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

  /* TRACE Component TOFfast2_zoom [34] */
  mccoordschange(mcposrTOFfast2_zoom, mcrotrTOFfast2_zoom,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component TOFfast2_zoom (without coords transformations) */
  mcJumpTrace_TOFfast2_zoom:
  SIG_MESSAGE("TOFfast2_zoom (Trace)");
  mcDEBUG_COMP("TOFfast2_zoom")
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

#define mcabsorbComp mcabsorbCompTOFfast2_zoom
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
#define mccompcurname  TOFfast2_zoom
#define mccompcurtype  TOF_monitor
#define mccompcurindex 34
#define nt mccTOFfast2_zoom_nt
#define TOF_N mccTOFfast2_zoom_TOF_N
#define TOF_p mccTOFfast2_zoom_TOF_p
#define TOF_p2 mccTOFfast2_zoom_TOF_p2
#define t_min mccTOFfast2_zoom_t_min
#define t_max mccTOFfast2_zoom_t_max
#define delta_t mccTOFfast2_zoom_delta_t
{   /* Declarations of TOFfast2_zoom=TOF_monitor() SETTING parameters. */
char* filename = mccTOFfast2_zoom_filename;
MCNUM xmin = mccTOFfast2_zoom_xmin;
MCNUM xmax = mccTOFfast2_zoom_xmax;
MCNUM ymin = mccTOFfast2_zoom_ymin;
MCNUM ymax = mccTOFfast2_zoom_ymax;
MCNUM xwidth = mccTOFfast2_zoom_xwidth;
MCNUM yheight = mccTOFfast2_zoom_yheight;
MCNUM tmin = mccTOFfast2_zoom_tmin;
MCNUM tmax = mccTOFfast2_zoom_tmax;
MCNUM dt = mccTOFfast2_zoom_dt;
MCNUM restore_neutron = mccTOFfast2_zoom_restore_neutron;
int nowritefile = mccTOFfast2_zoom_nowritefile;
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
#line 22763 "./ESS_IN5_reprate.c"
}   /* End of TOFfast2_zoom=TOF_monitor() SETTING parameter declarations. */
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
  mcabsorbCompTOFfast2_zoom:
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

  /* TRACE Component Lmon_afterfast2 [35] */
  mccoordschange(mcposrLmon_afterfast2, mcrotrLmon_afterfast2,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component Lmon_afterfast2 (without coords transformations) */
  mcJumpTrace_Lmon_afterfast2:
  SIG_MESSAGE("Lmon_afterfast2 (Trace)");
  mcDEBUG_COMP("Lmon_afterfast2")
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

#define mcabsorbComp mcabsorbCompLmon_afterfast2
  STORE_NEUTRON(35,
    mcnlx,
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
  mcNCounter[35]++;
  mcPCounter[35] += p;
  mcP2Counter[35] += p*p;
#define mccompcurname  Lmon_afterfast2
#define mccompcurtype  L_monitor
#define mccompcurindex 35
#define nL mccLmon_afterfast2_nL
#define L_N mccLmon_afterfast2_L_N
#define L_p mccLmon_afterfast2_L_p
#define L_p2 mccLmon_afterfast2_L_p2
{   /* Declarations of Lmon_afterfast2=L_monitor() SETTING parameters. */
char* filename = mccLmon_afterfast2_filename;
MCNUM xmin = mccLmon_afterfast2_xmin;
MCNUM xmax = mccLmon_afterfast2_xmax;
MCNUM ymin = mccLmon_afterfast2_ymin;
MCNUM ymax = mccLmon_afterfast2_ymax;
MCNUM xwidth = mccLmon_afterfast2_xwidth;
MCNUM yheight = mccLmon_afterfast2_yheight;
MCNUM Lmin = mccLmon_afterfast2_Lmin;
MCNUM Lmax = mccLmon_afterfast2_Lmax;
MCNUM restore_neutron = mccLmon_afterfast2_restore_neutron;
int nowritefile = mccLmon_afterfast2_nowritefile;
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
#line 22913 "./ESS_IN5_reprate.c"
}   /* End of Lmon_afterfast2=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompLmon_afterfast2:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(35,
      mcnlx,
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

  /* TRACE Component TOFL_afterfast2 [36] */
  mccoordschange(mcposrTOFL_afterfast2, mcrotrTOFL_afterfast2,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component TOFL_afterfast2 (without coords transformations) */
  mcJumpTrace_TOFL_afterfast2:
  SIG_MESSAGE("TOFL_afterfast2 (Trace)");
  mcDEBUG_COMP("TOFL_afterfast2")
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

#define mcabsorbComp mcabsorbCompTOFL_afterfast2
  STORE_NEUTRON(36,
    mcnlx,
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
  mcNCounter[36]++;
  mcPCounter[36] += p;
  mcP2Counter[36] += p*p;
#define mccompcurname  TOFL_afterfast2
#define mccompcurtype  TOFLambda_monitor
#define mccompcurindex 36
#define nL mccTOFL_afterfast2_nL
#define nt mccTOFL_afterfast2_nt
#define tmin mccTOFL_afterfast2_tmin
#define tmax mccTOFL_afterfast2_tmax
#define tt_0 mccTOFL_afterfast2_tt_0
#define tt_1 mccTOFL_afterfast2_tt_1
#define TOFL_N mccTOFL_afterfast2_TOFL_N
#define TOFL_p mccTOFL_afterfast2_TOFL_p
#define TOFL_p2 mccTOFL_afterfast2_TOFL_p2
{   /* Declarations of TOFL_afterfast2=TOFLambda_monitor() SETTING parameters. */
char* filename = mccTOFL_afterfast2_filename;
MCNUM xmin = mccTOFL_afterfast2_xmin;
MCNUM xmax = mccTOFL_afterfast2_xmax;
MCNUM ymin = mccTOFL_afterfast2_ymin;
MCNUM ymax = mccTOFL_afterfast2_ymax;
MCNUM xwidth = mccTOFL_afterfast2_xwidth;
MCNUM yheight = mccTOFL_afterfast2_yheight;
MCNUM Lmin = mccTOFL_afterfast2_Lmin;
MCNUM Lmax = mccTOFL_afterfast2_Lmax;
MCNUM restore_neutron = mccTOFL_afterfast2_restore_neutron;
int nowritefile = mccTOFL_afterfast2_nowritefile;
#line 85 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/TOFLambda_monitor.comp"
{
    int i,j;
    double div;
    double lambda;

    PROP_Z0;
    lambda = (2*PI/V2K)/sqrt(vx*vx + vy*vy + vz*vz);
    if (x>xmin && x<xmax && y>ymin && y<ymax &&
        lambda > Lmin && lambda < Lmax)
    {
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
}
#line 23068 "./ESS_IN5_reprate.c"
}   /* End of TOFL_afterfast2=TOFLambda_monitor() SETTING parameter declarations. */
#undef TOFL_p2
#undef TOFL_p
#undef TOFL_N
#undef tt_1
#undef tt_0
#undef tmax
#undef tmin
#undef nt
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompTOFL_afterfast2:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(36,
      mcnlx,
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

  /* TRACE Component TOFL_afterfast2_zoom [37] */
  mccoordschange(mcposrTOFL_afterfast2_zoom, mcrotrTOFL_afterfast2_zoom,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component TOFL_afterfast2_zoom (without coords transformations) */
  mcJumpTrace_TOFL_afterfast2_zoom:
  SIG_MESSAGE("TOFL_afterfast2_zoom (Trace)");
  mcDEBUG_COMP("TOFL_afterfast2_zoom")
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

#define mcabsorbComp mcabsorbCompTOFL_afterfast2_zoom
  STORE_NEUTRON(37,
    mcnlx,
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
  mcNCounter[37]++;
  mcPCounter[37] += p;
  mcP2Counter[37] += p*p;
#define mccompcurname  TOFL_afterfast2_zoom
#define mccompcurtype  TOFLambda_monitor
#define mccompcurindex 37
#define nL mccTOFL_afterfast2_zoom_nL
#define nt mccTOFL_afterfast2_zoom_nt
#define tmin mccTOFL_afterfast2_zoom_tmin
#define tmax mccTOFL_afterfast2_zoom_tmax
#define tt_0 mccTOFL_afterfast2_zoom_tt_0
#define tt_1 mccTOFL_afterfast2_zoom_tt_1
#define TOFL_N mccTOFL_afterfast2_zoom_TOFL_N
#define TOFL_p mccTOFL_afterfast2_zoom_TOFL_p
#define TOFL_p2 mccTOFL_afterfast2_zoom_TOFL_p2
{   /* Declarations of TOFL_afterfast2_zoom=TOFLambda_monitor() SETTING parameters. */
char* filename = mccTOFL_afterfast2_zoom_filename;
MCNUM xmin = mccTOFL_afterfast2_zoom_xmin;
MCNUM xmax = mccTOFL_afterfast2_zoom_xmax;
MCNUM ymin = mccTOFL_afterfast2_zoom_ymin;
MCNUM ymax = mccTOFL_afterfast2_zoom_ymax;
MCNUM xwidth = mccTOFL_afterfast2_zoom_xwidth;
MCNUM yheight = mccTOFL_afterfast2_zoom_yheight;
MCNUM Lmin = mccTOFL_afterfast2_zoom_Lmin;
MCNUM Lmax = mccTOFL_afterfast2_zoom_Lmax;
MCNUM restore_neutron = mccTOFL_afterfast2_zoom_restore_neutron;
int nowritefile = mccTOFL_afterfast2_zoom_nowritefile;
#line 85 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/TOFLambda_monitor.comp"
{
    int i,j;
    double div;
    double lambda;

    PROP_Z0;
    lambda = (2*PI/V2K)/sqrt(vx*vx + vy*vy + vz*vz);
    if (x>xmin && x<xmax && y>ymin && y<ymax &&
        lambda > Lmin && lambda < Lmax)
    {
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
}
#line 23228 "./ESS_IN5_reprate.c"
}   /* End of TOFL_afterfast2_zoom=TOFLambda_monitor() SETTING parameter declarations. */
#undef TOFL_p2
#undef TOFL_p
#undef TOFL_N
#undef tt_1
#undef tt_0
#undef tmax
#undef tmin
#undef nt
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompTOFL_afterfast2_zoom:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(37,
      mcnlx,
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

  /* TRACE Component PSD_afterfast2 [38] */
  mccoordschange(mcposrPSD_afterfast2, mcrotrPSD_afterfast2,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component PSD_afterfast2 (without coords transformations) */
  mcJumpTrace_PSD_afterfast2:
  SIG_MESSAGE("PSD_afterfast2 (Trace)");
  mcDEBUG_COMP("PSD_afterfast2")
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

#define mcabsorbComp mcabsorbCompPSD_afterfast2
  STORE_NEUTRON(38,
    mcnlx,
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
  mcNCounter[38]++;
  mcPCounter[38] += p;
  mcP2Counter[38] += p*p;
#define mccompcurname  PSD_afterfast2
#define mccompcurtype  PSD_monitor
#define mccompcurindex 38
#define PSD_N mccPSD_afterfast2_PSD_N
#define PSD_p mccPSD_afterfast2_PSD_p
#define PSD_p2 mccPSD_afterfast2_PSD_p2
{   /* Declarations of PSD_afterfast2=PSD_monitor() SETTING parameters. */
int nx = mccPSD_afterfast2_nx;
int ny = mccPSD_afterfast2_ny;
char* filename = mccPSD_afterfast2_filename;
MCNUM xmin = mccPSD_afterfast2_xmin;
MCNUM xmax = mccPSD_afterfast2_xmax;
MCNUM ymin = mccPSD_afterfast2_ymin;
MCNUM ymax = mccPSD_afterfast2_ymax;
MCNUM xwidth = mccPSD_afterfast2_xwidth;
MCNUM yheight = mccPSD_afterfast2_yheight;
MCNUM restore_neutron = mccPSD_afterfast2_restore_neutron;
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
#line 23371 "./ESS_IN5_reprate.c"
}   /* End of PSD_afterfast2=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompPSD_afterfast2:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(38,
      mcnlx,
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

  /* TRACE Component Guidesample [39] */
  mccoordschange(mcposrGuidesample, mcrotrGuidesample,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component Guidesample (without coords transformations) */
  mcJumpTrace_Guidesample:
  SIG_MESSAGE("Guidesample (Trace)");
  mcDEBUG_COMP("Guidesample")
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

#define mcabsorbComp mcabsorbCompGuidesample
  STORE_NEUTRON(39,
    mcnlx,
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
  mcNCounter[39]++;
  mcPCounter[39] += p;
  mcP2Counter[39] += p*p;
#define mccompcurname  Guidesample
#define mccompcurtype  Guide
#define mccompcurindex 39
#define pTable mccGuidesample_pTable
{   /* Declarations of Guidesample=Guide() SETTING parameters. */
char* reflect = mccGuidesample_reflect;
MCNUM w1 = mccGuidesample_w1;
MCNUM h1 = mccGuidesample_h1;
MCNUM w2 = mccGuidesample_w2;
MCNUM h2 = mccGuidesample_h2;
MCNUM l = mccGuidesample_l;
MCNUM R0 = mccGuidesample_R0;
MCNUM Qc = mccGuidesample_Qc;
MCNUM alpha = mccGuidesample_alpha;
MCNUM m = mccGuidesample_m;
MCNUM W = mccGuidesample_W;
#line 93 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide.comp"
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
#line 23599 "./ESS_IN5_reprate.c"
}   /* End of Guidesample=Guide() SETTING parameter declarations. */
#undef pTable
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompGuidesample:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(39,
      mcnlx,
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

  /* TRACE Component Lmon_guideend [40] */
  mccoordschange(mcposrLmon_guideend, mcrotrLmon_guideend,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component Lmon_guideend (without coords transformations) */
  mcJumpTrace_Lmon_guideend:
  SIG_MESSAGE("Lmon_guideend (Trace)");
  mcDEBUG_COMP("Lmon_guideend")
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

#define mcabsorbComp mcabsorbCompLmon_guideend
  STORE_NEUTRON(40,
    mcnlx,
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
  mcNCounter[40]++;
  mcPCounter[40] += p;
  mcP2Counter[40] += p*p;
#define mccompcurname  Lmon_guideend
#define mccompcurtype  L_monitor
#define mccompcurindex 40
#define nL mccLmon_guideend_nL
#define L_N mccLmon_guideend_L_N
#define L_p mccLmon_guideend_L_p
#define L_p2 mccLmon_guideend_L_p2
{   /* Declarations of Lmon_guideend=L_monitor() SETTING parameters. */
char* filename = mccLmon_guideend_filename;
MCNUM xmin = mccLmon_guideend_xmin;
MCNUM xmax = mccLmon_guideend_xmax;
MCNUM ymin = mccLmon_guideend_ymin;
MCNUM ymax = mccLmon_guideend_ymax;
MCNUM xwidth = mccLmon_guideend_xwidth;
MCNUM yheight = mccLmon_guideend_yheight;
MCNUM Lmin = mccLmon_guideend_Lmin;
MCNUM Lmax = mccLmon_guideend_Lmax;
MCNUM restore_neutron = mccLmon_guideend_restore_neutron;
int nowritefile = mccLmon_guideend_nowritefile;
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
#line 23743 "./ESS_IN5_reprate.c"
}   /* End of Lmon_guideend=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompLmon_guideend:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(40,
      mcnlx,
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

  /* TRACE Component PSDsample [41] */
  mccoordschange(mcposrPSDsample, mcrotrPSDsample,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component PSDsample (without coords transformations) */
  mcJumpTrace_PSDsample:
  SIG_MESSAGE("PSDsample (Trace)");
  mcDEBUG_COMP("PSDsample")
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

#define mcabsorbComp mcabsorbCompPSDsample
  STORE_NEUTRON(41,
    mcnlx,
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
  mcNCounter[41]++;
  mcPCounter[41] += p;
  mcP2Counter[41] += p*p;
#define mccompcurname  PSDsample
#define mccompcurtype  PSD_monitor
#define mccompcurindex 41
#define PSD_N mccPSDsample_PSD_N
#define PSD_p mccPSDsample_PSD_p
#define PSD_p2 mccPSDsample_PSD_p2
{   /* Declarations of PSDsample=PSD_monitor() SETTING parameters. */
int nx = mccPSDsample_nx;
int ny = mccPSDsample_ny;
char* filename = mccPSDsample_filename;
MCNUM xmin = mccPSDsample_xmin;
MCNUM xmax = mccPSDsample_xmax;
MCNUM ymin = mccPSDsample_ymin;
MCNUM ymax = mccPSDsample_ymax;
MCNUM xwidth = mccPSDsample_xwidth;
MCNUM yheight = mccPSDsample_yheight;
MCNUM restore_neutron = mccPSDsample_restore_neutron;
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
#line 23881 "./ESS_IN5_reprate.c"
}   /* End of PSDsample=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompPSDsample:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(41,
      mcnlx,
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

  /* TRACE Component TOFsample_zoom [42] */
  mccoordschange(mcposrTOFsample_zoom, mcrotrTOFsample_zoom,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component TOFsample_zoom (without coords transformations) */
  mcJumpTrace_TOFsample_zoom:
  SIG_MESSAGE("TOFsample_zoom (Trace)");
  mcDEBUG_COMP("TOFsample_zoom")
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

#define mcabsorbComp mcabsorbCompTOFsample_zoom
  STORE_NEUTRON(42,
    mcnlx,
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
  mcNCounter[42]++;
  mcPCounter[42] += p;
  mcP2Counter[42] += p*p;
#define mccompcurname  TOFsample_zoom
#define mccompcurtype  TOF_monitor
#define mccompcurindex 42
#define nt mccTOFsample_zoom_nt
#define TOF_N mccTOFsample_zoom_TOF_N
#define TOF_p mccTOFsample_zoom_TOF_p
#define TOF_p2 mccTOFsample_zoom_TOF_p2
#define t_min mccTOFsample_zoom_t_min
#define t_max mccTOFsample_zoom_t_max
#define delta_t mccTOFsample_zoom_delta_t
{   /* Declarations of TOFsample_zoom=TOF_monitor() SETTING parameters. */
char* filename = mccTOFsample_zoom_filename;
MCNUM xmin = mccTOFsample_zoom_xmin;
MCNUM xmax = mccTOFsample_zoom_xmax;
MCNUM ymin = mccTOFsample_zoom_ymin;
MCNUM ymax = mccTOFsample_zoom_ymax;
MCNUM xwidth = mccTOFsample_zoom_xwidth;
MCNUM yheight = mccTOFsample_zoom_yheight;
MCNUM tmin = mccTOFsample_zoom_tmin;
MCNUM tmax = mccTOFsample_zoom_tmax;
MCNUM dt = mccTOFsample_zoom_dt;
MCNUM restore_neutron = mccTOFsample_zoom_restore_neutron;
int nowritefile = mccTOFsample_zoom_nowritefile;
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
#line 24028 "./ESS_IN5_reprate.c"
}   /* End of TOFsample_zoom=TOF_monitor() SETTING parameter declarations. */
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
  mcabsorbCompTOFsample_zoom:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(42,
      mcnlx,
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

  /* TRACE Component Esample [43] */
  mccoordschange(mcposrEsample, mcrotrEsample,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component Esample (without coords transformations) */
  mcJumpTrace_Esample:
  SIG_MESSAGE("Esample (Trace)");
  mcDEBUG_COMP("Esample")
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

#define mcabsorbComp mcabsorbCompEsample
  STORE_NEUTRON(43,
    mcnlx,
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
  mcNCounter[43]++;
  mcPCounter[43] += p;
  mcP2Counter[43] += p*p;
#define mccompcurname  Esample
#define mccompcurtype  E_monitor
#define mccompcurindex 43
#define nE mccEsample_nE
#define E_N mccEsample_E_N
#define E_p mccEsample_E_p
#define E_p2 mccEsample_E_p2
#define S_p mccEsample_S_p
#define S_pE mccEsample_S_pE
#define S_pE2 mccEsample_S_pE2
{   /* Declarations of Esample=E_monitor() SETTING parameters. */
char* filename = mccEsample_filename;
MCNUM xmin = mccEsample_xmin;
MCNUM xmax = mccEsample_xmax;
MCNUM ymin = mccEsample_ymin;
MCNUM ymax = mccEsample_ymax;
MCNUM xwidth = mccEsample_xwidth;
MCNUM yheight = mccEsample_yheight;
MCNUM Emin = mccEsample_Emin;
MCNUM Emax = mccEsample_Emax;
MCNUM restore_neutron = mccEsample_restore_neutron;
int nowritefile = mccEsample_nowritefile;
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
#line 24186 "./ESS_IN5_reprate.c"
}   /* End of Esample=E_monitor() SETTING parameter declarations. */
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
  mcabsorbCompEsample:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(43,
      mcnlx,
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

  /* TRACE Component Lmon_sample_zoom [44] */
  mccoordschange(mcposrLmon_sample_zoom, mcrotrLmon_sample_zoom,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component Lmon_sample_zoom (without coords transformations) */
  mcJumpTrace_Lmon_sample_zoom:
  SIG_MESSAGE("Lmon_sample_zoom (Trace)");
  mcDEBUG_COMP("Lmon_sample_zoom")
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

#define mcabsorbComp mcabsorbCompLmon_sample_zoom
  STORE_NEUTRON(44,
    mcnlx,
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
  mcNCounter[44]++;
  mcPCounter[44] += p;
  mcP2Counter[44] += p*p;
#define mccompcurname  Lmon_sample_zoom
#define mccompcurtype  L_monitor
#define mccompcurindex 44
#define nL mccLmon_sample_zoom_nL
#define L_N mccLmon_sample_zoom_L_N
#define L_p mccLmon_sample_zoom_L_p
#define L_p2 mccLmon_sample_zoom_L_p2
{   /* Declarations of Lmon_sample_zoom=L_monitor() SETTING parameters. */
char* filename = mccLmon_sample_zoom_filename;
MCNUM xmin = mccLmon_sample_zoom_xmin;
MCNUM xmax = mccLmon_sample_zoom_xmax;
MCNUM ymin = mccLmon_sample_zoom_ymin;
MCNUM ymax = mccLmon_sample_zoom_ymax;
MCNUM xwidth = mccLmon_sample_zoom_xwidth;
MCNUM yheight = mccLmon_sample_zoom_yheight;
MCNUM Lmin = mccLmon_sample_zoom_Lmin;
MCNUM Lmax = mccLmon_sample_zoom_Lmax;
MCNUM restore_neutron = mccLmon_sample_zoom_restore_neutron;
int nowritefile = mccLmon_sample_zoom_nowritefile;
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
#line 24336 "./ESS_IN5_reprate.c"
}   /* End of Lmon_sample_zoom=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompLmon_sample_zoom:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(44,
      mcnlx,
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

  /* TRACE Component sample [45] */
  mccoordschange(mcposrsample, mcrotrsample,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component sample (without coords transformations) */
  mcJumpTrace_sample:
  SIG_MESSAGE("sample (Trace)");
  mcDEBUG_COMP("sample")
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

#define mcabsorbComp mcabsorbCompsample
  STORE_NEUTRON(45,
    mcnlx,
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
  mcNCounter[45]++;
  mcPCounter[45] += p;
  mcP2Counter[45] += p*p;
#define mccompcurname  sample
#define mccompcurtype  Tunneling_sample
#define mccompcurindex 45
#define ftun mccsample_ftun
#define fQE mccsample_fQE
#define VarsV mccsample_VarsV
{   /* Declarations of sample=Tunneling_sample() SETTING parameters. */
MCNUM thickness = mccsample_thickness;
MCNUM radius = mccsample_radius;
MCNUM focus_r = mccsample_focus_r;
MCNUM p_interact = mccsample_p_interact;
MCNUM f_QE = mccsample_f_QE;
MCNUM f_tun = mccsample_f_tun;
MCNUM gamma = mccsample_gamma;
MCNUM E_tun = mccsample_E_tun;
MCNUM target_x = mccsample_target_x;
MCNUM target_y = mccsample_target_y;
MCNUM target_z = mccsample_target_z;
MCNUM focus_xw = mccsample_focus_xw;
MCNUM focus_yh = mccsample_focus_yh;
MCNUM focus_aw = mccsample_focus_aw;
MCNUM focus_ah = mccsample_focus_ah;
MCNUM xwidth = mccsample_xwidth;
MCNUM yheight = mccsample_yheight;
MCNUM zdepth = mccsample_zdepth;
MCNUM sigma_abs = mccsample_sigma_abs;
MCNUM sigma_inc = mccsample_sigma_inc;
MCNUM Vc = mccsample_Vc;
int target_index = mccsample_target_index;
#line 185 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../samples/Tunneling_sample.comp"
{
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
}
#line 24606 "./ESS_IN5_reprate.c"
}   /* End of sample=Tunneling_sample() SETTING parameter declarations. */
#undef VarsV
#undef fQE
#undef ftun
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompsample:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(45,
      mcnlx,
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

  /* TRACE Component detectorarm [46] */
  mccoordschange(mcposrdetectorarm, mcrotrdetectorarm,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component detectorarm (without coords transformations) */
  mcJumpTrace_detectorarm:
  SIG_MESSAGE("detectorarm (Trace)");
  mcDEBUG_COMP("detectorarm")
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

#define mcabsorbComp mcabsorbCompdetectorarm
  STORE_NEUTRON(46,
    mcnlx,
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
  mcNCounter[46]++;
  mcPCounter[46] += p;
  mcP2Counter[46] += p*p;
#define mccompcurname  detectorarm
#define mccompcurtype  Arm
#define mccompcurindex 46
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompdetectorarm:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(46,
      mcnlx,
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

  /* TRACE Component TOFdetector [47] */
  mccoordschange(mcposrTOFdetector, mcrotrTOFdetector,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component TOFdetector (without coords transformations) */
  mcJumpTrace_TOFdetector:
  SIG_MESSAGE("TOFdetector (Trace)");
  mcDEBUG_COMP("TOFdetector")
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

#define mcabsorbComp mcabsorbCompTOFdetector
  STORE_NEUTRON(47,
    mcnlx,
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
  mcNCounter[47]++;
  mcPCounter[47] += p;
  mcP2Counter[47] += p*p;
#define mccompcurname  TOFdetector
#define mccompcurtype  TOF_monitor
#define mccompcurindex 47
#define nt mccTOFdetector_nt
#define TOF_N mccTOFdetector_TOF_N
#define TOF_p mccTOFdetector_TOF_p
#define TOF_p2 mccTOFdetector_TOF_p2
#define t_min mccTOFdetector_t_min
#define t_max mccTOFdetector_t_max
#define delta_t mccTOFdetector_delta_t
{   /* Declarations of TOFdetector=TOF_monitor() SETTING parameters. */
char* filename = mccTOFdetector_filename;
MCNUM xmin = mccTOFdetector_xmin;
MCNUM xmax = mccTOFdetector_xmax;
MCNUM ymin = mccTOFdetector_ymin;
MCNUM ymax = mccTOFdetector_ymax;
MCNUM xwidth = mccTOFdetector_xwidth;
MCNUM yheight = mccTOFdetector_yheight;
MCNUM tmin = mccTOFdetector_tmin;
MCNUM tmax = mccTOFdetector_tmax;
MCNUM dt = mccTOFdetector_dt;
MCNUM restore_neutron = mccTOFdetector_restore_neutron;
int nowritefile = mccTOFdetector_nowritefile;
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
#line 24856 "./ESS_IN5_reprate.c"
}   /* End of TOFdetector=TOF_monitor() SETTING parameter declarations. */
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
  mcabsorbCompTOFdetector:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(47,
      mcnlx,
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

  /* TRACE Component TOFdetector_zoom [48] */
  mccoordschange(mcposrTOFdetector_zoom, mcrotrTOFdetector_zoom,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component TOFdetector_zoom (without coords transformations) */
  mcJumpTrace_TOFdetector_zoom:
  SIG_MESSAGE("TOFdetector_zoom (Trace)");
  mcDEBUG_COMP("TOFdetector_zoom")
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

#define mcabsorbComp mcabsorbCompTOFdetector_zoom
  STORE_NEUTRON(48,
    mcnlx,
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
  mcNCounter[48]++;
  mcPCounter[48] += p;
  mcP2Counter[48] += p*p;
#define mccompcurname  TOFdetector_zoom
#define mccompcurtype  TOF_monitor
#define mccompcurindex 48
#define nt mccTOFdetector_zoom_nt
#define TOF_N mccTOFdetector_zoom_TOF_N
#define TOF_p mccTOFdetector_zoom_TOF_p
#define TOF_p2 mccTOFdetector_zoom_TOF_p2
#define t_min mccTOFdetector_zoom_t_min
#define t_max mccTOFdetector_zoom_t_max
#define delta_t mccTOFdetector_zoom_delta_t
{   /* Declarations of TOFdetector_zoom=TOF_monitor() SETTING parameters. */
char* filename = mccTOFdetector_zoom_filename;
MCNUM xmin = mccTOFdetector_zoom_xmin;
MCNUM xmax = mccTOFdetector_zoom_xmax;
MCNUM ymin = mccTOFdetector_zoom_ymin;
MCNUM ymax = mccTOFdetector_zoom_ymax;
MCNUM xwidth = mccTOFdetector_zoom_xwidth;
MCNUM yheight = mccTOFdetector_zoom_yheight;
MCNUM tmin = mccTOFdetector_zoom_tmin;
MCNUM tmax = mccTOFdetector_zoom_tmax;
MCNUM dt = mccTOFdetector_zoom_dt;
MCNUM restore_neutron = mccTOFdetector_zoom_restore_neutron;
int nowritefile = mccTOFdetector_zoom_nowritefile;
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
#line 25007 "./ESS_IN5_reprate.c"
}   /* End of TOFdetector_zoom=TOF_monitor() SETTING parameter declarations. */
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
  mcabsorbCompTOFdetector_zoom:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(48,
      mcnlx,
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

  /* TRACE Component Edetector [49] */
  mccoordschange(mcposrEdetector, mcrotrEdetector,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component Edetector (without coords transformations) */
  mcJumpTrace_Edetector:
  SIG_MESSAGE("Edetector (Trace)");
  mcDEBUG_COMP("Edetector")
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

#define mcabsorbComp mcabsorbCompEdetector
  STORE_NEUTRON(49,
    mcnlx,
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
  mcNCounter[49]++;
  mcPCounter[49] += p;
  mcP2Counter[49] += p*p;
#define mccompcurname  Edetector
#define mccompcurtype  E_monitor
#define mccompcurindex 49
#define nE mccEdetector_nE
#define E_N mccEdetector_E_N
#define E_p mccEdetector_E_p
#define E_p2 mccEdetector_E_p2
#define S_p mccEdetector_S_p
#define S_pE mccEdetector_S_pE
#define S_pE2 mccEdetector_S_pE2
{   /* Declarations of Edetector=E_monitor() SETTING parameters. */
char* filename = mccEdetector_filename;
MCNUM xmin = mccEdetector_xmin;
MCNUM xmax = mccEdetector_xmax;
MCNUM ymin = mccEdetector_ymin;
MCNUM ymax = mccEdetector_ymax;
MCNUM xwidth = mccEdetector_xwidth;
MCNUM yheight = mccEdetector_yheight;
MCNUM Emin = mccEdetector_Emin;
MCNUM Emax = mccEdetector_Emax;
MCNUM restore_neutron = mccEdetector_restore_neutron;
int nowritefile = mccEdetector_nowritefile;
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
#line 25165 "./ESS_IN5_reprate.c"
}   /* End of Edetector=E_monitor() SETTING parameter declarations. */
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
  mcabsorbCompEdetector:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(49,
      mcnlx,
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

  /* TRACE Component TOF2Edetector [50] */
  mccoordschange(mcposrTOF2Edetector, mcrotrTOF2Edetector,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component TOF2Edetector (without coords transformations) */
  mcJumpTrace_TOF2Edetector:
  SIG_MESSAGE("TOF2Edetector (Trace)");
  mcDEBUG_COMP("TOF2Edetector")
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

#define mcabsorbComp mcabsorbCompTOF2Edetector
  STORE_NEUTRON(50,
    mcnlx,
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
  mcNCounter[50]++;
  mcPCounter[50] += p;
  mcP2Counter[50] += p*p;
#define mccompcurname  TOF2Edetector
#define mccompcurtype  TOF2E_monitor
#define mccompcurindex 50
#define nE mccTOF2Edetector_nE
#define E_N mccTOF2Edetector_E_N
#define E_p mccTOF2Edetector_E_p
#define E_p2 mccTOF2Edetector_E_p2
#define S_p mccTOF2Edetector_S_p
#define S_pE mccTOF2Edetector_S_pE
#define S_pE2 mccTOF2Edetector_S_pE2
{   /* Declarations of TOF2Edetector=TOF2E_monitor() SETTING parameters. */
char* filename = mccTOF2Edetector_filename;
MCNUM xmin = mccTOF2Edetector_xmin;
MCNUM xmax = mccTOF2Edetector_xmax;
MCNUM ymin = mccTOF2Edetector_ymin;
MCNUM ymax = mccTOF2Edetector_ymax;
MCNUM xwidth = mccTOF2Edetector_xwidth;
MCNUM yheight = mccTOF2Edetector_yheight;
MCNUM Emin = mccTOF2Edetector_Emin;
MCNUM Emax = mccTOF2Edetector_Emax;
MCNUM T_zero = mccTOF2Edetector_T_zero;
MCNUM L_flight = mccTOF2Edetector_L_flight;
MCNUM restore_neutron = mccTOF2Edetector_restore_neutron;
int nowritefile = mccTOF2Edetector_nowritefile;
#line 90 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/TOF2E_monitor.comp"
{
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
}
#line 25325 "./ESS_IN5_reprate.c"
}   /* End of TOF2Edetector=TOF2E_monitor() SETTING parameter declarations. */
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
  mcabsorbCompTOF2Edetector:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(50,
      mcnlx,
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
#define mccompcurindex 2
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
#line 25441 "./ESS_IN5_reprate.c"
}   /* End of Origin=Progress_bar() SETTING parameter declarations. */
#undef CurrentTime
#undef EndTime
#undef StartTime
#undef IntermediateCnts
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'TOFmoderator_zoom'. */
  SIG_MESSAGE("TOFmoderator_zoom (Save)");
#define mccompcurname  TOFmoderator_zoom
#define mccompcurtype  TOF_monitor
#define mccompcurindex 3
#define nt mccTOFmoderator_zoom_nt
#define TOF_N mccTOFmoderator_zoom_TOF_N
#define TOF_p mccTOFmoderator_zoom_TOF_p
#define TOF_p2 mccTOFmoderator_zoom_TOF_p2
#define t_min mccTOFmoderator_zoom_t_min
#define t_max mccTOFmoderator_zoom_t_max
#define delta_t mccTOFmoderator_zoom_delta_t
{   /* Declarations of TOFmoderator_zoom=TOF_monitor() SETTING parameters. */
char* filename = mccTOFmoderator_zoom_filename;
MCNUM xmin = mccTOFmoderator_zoom_xmin;
MCNUM xmax = mccTOFmoderator_zoom_xmax;
MCNUM ymin = mccTOFmoderator_zoom_ymin;
MCNUM ymax = mccTOFmoderator_zoom_ymax;
MCNUM xwidth = mccTOFmoderator_zoom_xwidth;
MCNUM yheight = mccTOFmoderator_zoom_yheight;
MCNUM tmin = mccTOFmoderator_zoom_tmin;
MCNUM tmax = mccTOFmoderator_zoom_tmax;
MCNUM dt = mccTOFmoderator_zoom_dt;
MCNUM restore_neutron = mccTOFmoderator_zoom_restore_neutron;
int nowritefile = mccTOFmoderator_zoom_nowritefile;
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
#line 25488 "./ESS_IN5_reprate.c"
}   /* End of TOFmoderator_zoom=TOF_monitor() SETTING parameter declarations. */
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

  /* User SAVE code for component 'TOFmoderator'. */
  SIG_MESSAGE("TOFmoderator (Save)");
#define mccompcurname  TOFmoderator
#define mccompcurtype  TOF_monitor
#define mccompcurindex 4
#define nt mccTOFmoderator_nt
#define TOF_N mccTOFmoderator_TOF_N
#define TOF_p mccTOFmoderator_TOF_p
#define TOF_p2 mccTOFmoderator_TOF_p2
#define t_min mccTOFmoderator_t_min
#define t_max mccTOFmoderator_t_max
#define delta_t mccTOFmoderator_delta_t
{   /* Declarations of TOFmoderator=TOF_monitor() SETTING parameters. */
char* filename = mccTOFmoderator_filename;
MCNUM xmin = mccTOFmoderator_xmin;
MCNUM xmax = mccTOFmoderator_xmax;
MCNUM ymin = mccTOFmoderator_ymin;
MCNUM ymax = mccTOFmoderator_ymax;
MCNUM xwidth = mccTOFmoderator_xwidth;
MCNUM yheight = mccTOFmoderator_yheight;
MCNUM tmin = mccTOFmoderator_tmin;
MCNUM tmax = mccTOFmoderator_tmax;
MCNUM dt = mccTOFmoderator_dt;
MCNUM restore_neutron = mccTOFmoderator_restore_neutron;
int nowritefile = mccTOFmoderator_nowritefile;
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
#line 25538 "./ESS_IN5_reprate.c"
}   /* End of TOFmoderator=TOF_monitor() SETTING parameter declarations. */
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

  /* User SAVE code for component 'Lmon_guistart'. */
  SIG_MESSAGE("Lmon_guistart (Save)");
#define mccompcurname  Lmon_guistart
#define mccompcurtype  L_monitor
#define mccompcurindex 5
#define nL mccLmon_guistart_nL
#define L_N mccLmon_guistart_L_N
#define L_p mccLmon_guistart_L_p
#define L_p2 mccLmon_guistart_L_p2
{   /* Declarations of Lmon_guistart=L_monitor() SETTING parameters. */
char* filename = mccLmon_guistart_filename;
MCNUM xmin = mccLmon_guistart_xmin;
MCNUM xmax = mccLmon_guistart_xmax;
MCNUM ymin = mccLmon_guistart_ymin;
MCNUM ymax = mccLmon_guistart_ymax;
MCNUM xwidth = mccLmon_guistart_xwidth;
MCNUM yheight = mccLmon_guistart_yheight;
MCNUM Lmin = mccLmon_guistart_Lmin;
MCNUM Lmax = mccLmon_guistart_Lmax;
MCNUM restore_neutron = mccLmon_guistart_restore_neutron;
int nowritefile = mccLmon_guistart_nowritefile;
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
#line 25584 "./ESS_IN5_reprate.c"
}   /* End of Lmon_guistart=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'Lmon_normalize'. */
  SIG_MESSAGE("Lmon_normalize (Save)");
#define mccompcurname  Lmon_normalize
#define mccompcurtype  L_monitor
#define mccompcurindex 6
#define nL mccLmon_normalize_nL
#define L_N mccLmon_normalize_L_N
#define L_p mccLmon_normalize_L_p
#define L_p2 mccLmon_normalize_L_p2
{   /* Declarations of Lmon_normalize=L_monitor() SETTING parameters. */
char* filename = mccLmon_normalize_filename;
MCNUM xmin = mccLmon_normalize_xmin;
MCNUM xmax = mccLmon_normalize_xmax;
MCNUM ymin = mccLmon_normalize_ymin;
MCNUM ymax = mccLmon_normalize_ymax;
MCNUM xwidth = mccLmon_normalize_xwidth;
MCNUM yheight = mccLmon_normalize_yheight;
MCNUM Lmin = mccLmon_normalize_Lmin;
MCNUM Lmax = mccLmon_normalize_Lmax;
MCNUM restore_neutron = mccLmon_normalize_restore_neutron;
int nowritefile = mccLmon_normalize_nowritefile;
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
#line 25627 "./ESS_IN5_reprate.c"
}   /* End of Lmon_normalize=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'Lmonslow1'. */
  SIG_MESSAGE("Lmonslow1 (Save)");
#define mccompcurname  Lmonslow1
#define mccompcurtype  L_monitor
#define mccompcurindex 8
#define nL mccLmonslow1_nL
#define L_N mccLmonslow1_L_N
#define L_p mccLmonslow1_L_p
#define L_p2 mccLmonslow1_L_p2
{   /* Declarations of Lmonslow1=L_monitor() SETTING parameters. */
char* filename = mccLmonslow1_filename;
MCNUM xmin = mccLmonslow1_xmin;
MCNUM xmax = mccLmonslow1_xmax;
MCNUM ymin = mccLmonslow1_ymin;
MCNUM ymax = mccLmonslow1_ymax;
MCNUM xwidth = mccLmonslow1_xwidth;
MCNUM yheight = mccLmonslow1_yheight;
MCNUM Lmin = mccLmonslow1_Lmin;
MCNUM Lmax = mccLmonslow1_Lmax;
MCNUM restore_neutron = mccLmonslow1_restore_neutron;
int nowritefile = mccLmonslow1_nowritefile;
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
#line 25670 "./ESS_IN5_reprate.c"
}   /* End of Lmonslow1=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'PSDslow1'. */
  SIG_MESSAGE("PSDslow1 (Save)");
#define mccompcurname  PSDslow1
#define mccompcurtype  PSD_monitor
#define mccompcurindex 9
#define PSD_N mccPSDslow1_PSD_N
#define PSD_p mccPSDslow1_PSD_p
#define PSD_p2 mccPSDslow1_PSD_p2
{   /* Declarations of PSDslow1=PSD_monitor() SETTING parameters. */
int nx = mccPSDslow1_nx;
int ny = mccPSDslow1_ny;
char* filename = mccPSDslow1_filename;
MCNUM xmin = mccPSDslow1_xmin;
MCNUM xmax = mccPSDslow1_xmax;
MCNUM ymin = mccPSDslow1_ymin;
MCNUM ymax = mccPSDslow1_ymax;
MCNUM xwidth = mccPSDslow1_xwidth;
MCNUM yheight = mccPSDslow1_yheight;
MCNUM restore_neutron = mccPSDslow1_restore_neutron;
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
#line 25710 "./ESS_IN5_reprate.c"
}   /* End of PSDslow1=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'TOFLmon1'. */
  SIG_MESSAGE("TOFLmon1 (Save)");
#define mccompcurname  TOFLmon1
#define mccompcurtype  TOFLambda_monitor
#define mccompcurindex 11
#define nL mccTOFLmon1_nL
#define nt mccTOFLmon1_nt
#define tmin mccTOFLmon1_tmin
#define tmax mccTOFLmon1_tmax
#define tt_0 mccTOFLmon1_tt_0
#define tt_1 mccTOFLmon1_tt_1
#define TOFL_N mccTOFLmon1_TOFL_N
#define TOFL_p mccTOFLmon1_TOFL_p
#define TOFL_p2 mccTOFLmon1_TOFL_p2
{   /* Declarations of TOFLmon1=TOFLambda_monitor() SETTING parameters. */
char* filename = mccTOFLmon1_filename;
MCNUM xmin = mccTOFLmon1_xmin;
MCNUM xmax = mccTOFLmon1_xmax;
MCNUM ymin = mccTOFLmon1_ymin;
MCNUM ymax = mccTOFLmon1_ymax;
MCNUM xwidth = mccTOFLmon1_xwidth;
MCNUM yheight = mccTOFLmon1_yheight;
MCNUM Lmin = mccTOFLmon1_Lmin;
MCNUM Lmax = mccTOFLmon1_Lmax;
MCNUM restore_neutron = mccTOFLmon1_restore_neutron;
int nowritefile = mccTOFLmon1_nowritefile;
#line 110 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/TOFLambda_monitor.comp"
{
    if (!nowritefile) {
    DETECTOR_OUT_2D(
        "TOF-wavelength monitor",
        "Time-of-flight [\\gms]", "Wavelength [AA]",
        tmin, tmax, Lmin, Lmax,
        nt, nL,
        &TOFL_N[0][0],&TOFL_p[0][0],&TOFL_p2[0][0],
        filename);
    }
}
#line 25757 "./ESS_IN5_reprate.c"
}   /* End of TOFLmon1=TOFLambda_monitor() SETTING parameter declarations. */
#undef TOFL_p2
#undef TOFL_p
#undef TOFL_N
#undef tt_1
#undef tt_0
#undef tmax
#undef tmin
#undef nt
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'Lmon_afterslow1'. */
  SIG_MESSAGE("Lmon_afterslow1 (Save)");
#define mccompcurname  Lmon_afterslow1
#define mccompcurtype  L_monitor
#define mccompcurindex 12
#define nL mccLmon_afterslow1_nL
#define L_N mccLmon_afterslow1_L_N
#define L_p mccLmon_afterslow1_L_p
#define L_p2 mccLmon_afterslow1_L_p2
{   /* Declarations of Lmon_afterslow1=L_monitor() SETTING parameters. */
char* filename = mccLmon_afterslow1_filename;
MCNUM xmin = mccLmon_afterslow1_xmin;
MCNUM xmax = mccLmon_afterslow1_xmax;
MCNUM ymin = mccLmon_afterslow1_ymin;
MCNUM ymax = mccLmon_afterslow1_ymax;
MCNUM xwidth = mccLmon_afterslow1_xwidth;
MCNUM yheight = mccLmon_afterslow1_yheight;
MCNUM Lmin = mccLmon_afterslow1_Lmin;
MCNUM Lmax = mccLmon_afterslow1_Lmax;
MCNUM restore_neutron = mccLmon_afterslow1_restore_neutron;
int nowritefile = mccLmon_afterslow1_nowritefile;
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
#line 25805 "./ESS_IN5_reprate.c"
}   /* End of Lmon_afterslow1=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'PSD_afterslow1'. */
  SIG_MESSAGE("PSD_afterslow1 (Save)");
#define mccompcurname  PSD_afterslow1
#define mccompcurtype  PSD_monitor
#define mccompcurindex 13
#define PSD_N mccPSD_afterslow1_PSD_N
#define PSD_p mccPSD_afterslow1_PSD_p
#define PSD_p2 mccPSD_afterslow1_PSD_p2
{   /* Declarations of PSD_afterslow1=PSD_monitor() SETTING parameters. */
int nx = mccPSD_afterslow1_nx;
int ny = mccPSD_afterslow1_ny;
char* filename = mccPSD_afterslow1_filename;
MCNUM xmin = mccPSD_afterslow1_xmin;
MCNUM xmax = mccPSD_afterslow1_xmax;
MCNUM ymin = mccPSD_afterslow1_ymin;
MCNUM ymax = mccPSD_afterslow1_ymax;
MCNUM xwidth = mccPSD_afterslow1_xwidth;
MCNUM yheight = mccPSD_afterslow1_yheight;
MCNUM restore_neutron = mccPSD_afterslow1_restore_neutron;
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
#line 25845 "./ESS_IN5_reprate.c"
}   /* End of PSD_afterslow1=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'Lmon_slow2'. */
  SIG_MESSAGE("Lmon_slow2 (Save)");
#define mccompcurname  Lmon_slow2
#define mccompcurtype  L_monitor
#define mccompcurindex 16
#define nL mccLmon_slow2_nL
#define L_N mccLmon_slow2_L_N
#define L_p mccLmon_slow2_L_p
#define L_p2 mccLmon_slow2_L_p2
{   /* Declarations of Lmon_slow2=L_monitor() SETTING parameters. */
char* filename = mccLmon_slow2_filename;
MCNUM xmin = mccLmon_slow2_xmin;
MCNUM xmax = mccLmon_slow2_xmax;
MCNUM ymin = mccLmon_slow2_ymin;
MCNUM ymax = mccLmon_slow2_ymax;
MCNUM xwidth = mccLmon_slow2_xwidth;
MCNUM yheight = mccLmon_slow2_yheight;
MCNUM Lmin = mccLmon_slow2_Lmin;
MCNUM Lmax = mccLmon_slow2_Lmax;
MCNUM restore_neutron = mccLmon_slow2_restore_neutron;
int nowritefile = mccLmon_slow2_nowritefile;
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
#line 25887 "./ESS_IN5_reprate.c"
}   /* End of Lmon_slow2=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'PSD_afterslow2'. */
  SIG_MESSAGE("PSD_afterslow2 (Save)");
#define mccompcurname  PSD_afterslow2
#define mccompcurtype  PSD_monitor
#define mccompcurindex 19
#define PSD_N mccPSD_afterslow2_PSD_N
#define PSD_p mccPSD_afterslow2_PSD_p
#define PSD_p2 mccPSD_afterslow2_PSD_p2
{   /* Declarations of PSD_afterslow2=PSD_monitor() SETTING parameters. */
int nx = mccPSD_afterslow2_nx;
int ny = mccPSD_afterslow2_ny;
char* filename = mccPSD_afterslow2_filename;
MCNUM xmin = mccPSD_afterslow2_xmin;
MCNUM xmax = mccPSD_afterslow2_xmax;
MCNUM ymin = mccPSD_afterslow2_ymin;
MCNUM ymax = mccPSD_afterslow2_ymax;
MCNUM xwidth = mccPSD_afterslow2_xwidth;
MCNUM yheight = mccPSD_afterslow2_yheight;
MCNUM restore_neutron = mccPSD_afterslow2_restore_neutron;
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
#line 25927 "./ESS_IN5_reprate.c"
}   /* End of PSD_afterslow2=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'Lmon_afterslow2'. */
  SIG_MESSAGE("Lmon_afterslow2 (Save)");
#define mccompcurname  Lmon_afterslow2
#define mccompcurtype  L_monitor
#define mccompcurindex 20
#define nL mccLmon_afterslow2_nL
#define L_N mccLmon_afterslow2_L_N
#define L_p mccLmon_afterslow2_L_p
#define L_p2 mccLmon_afterslow2_L_p2
{   /* Declarations of Lmon_afterslow2=L_monitor() SETTING parameters. */
char* filename = mccLmon_afterslow2_filename;
MCNUM xmin = mccLmon_afterslow2_xmin;
MCNUM xmax = mccLmon_afterslow2_xmax;
MCNUM ymin = mccLmon_afterslow2_ymin;
MCNUM ymax = mccLmon_afterslow2_ymax;
MCNUM xwidth = mccLmon_afterslow2_xwidth;
MCNUM yheight = mccLmon_afterslow2_yheight;
MCNUM Lmin = mccLmon_afterslow2_Lmin;
MCNUM Lmax = mccLmon_afterslow2_Lmax;
MCNUM restore_neutron = mccLmon_afterslow2_restore_neutron;
int nowritefile = mccLmon_afterslow2_nowritefile;
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
#line 25969 "./ESS_IN5_reprate.c"
}   /* End of Lmon_afterslow2=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'TOFL_afterslow2'. */
  SIG_MESSAGE("TOFL_afterslow2 (Save)");
#define mccompcurname  TOFL_afterslow2
#define mccompcurtype  TOFLambda_monitor
#define mccompcurindex 21
#define nL mccTOFL_afterslow2_nL
#define nt mccTOFL_afterslow2_nt
#define tmin mccTOFL_afterslow2_tmin
#define tmax mccTOFL_afterslow2_tmax
#define tt_0 mccTOFL_afterslow2_tt_0
#define tt_1 mccTOFL_afterslow2_tt_1
#define TOFL_N mccTOFL_afterslow2_TOFL_N
#define TOFL_p mccTOFL_afterslow2_TOFL_p
#define TOFL_p2 mccTOFL_afterslow2_TOFL_p2
{   /* Declarations of TOFL_afterslow2=TOFLambda_monitor() SETTING parameters. */
char* filename = mccTOFL_afterslow2_filename;
MCNUM xmin = mccTOFL_afterslow2_xmin;
MCNUM xmax = mccTOFL_afterslow2_xmax;
MCNUM ymin = mccTOFL_afterslow2_ymin;
MCNUM ymax = mccTOFL_afterslow2_ymax;
MCNUM xwidth = mccTOFL_afterslow2_xwidth;
MCNUM yheight = mccTOFL_afterslow2_yheight;
MCNUM Lmin = mccTOFL_afterslow2_Lmin;
MCNUM Lmax = mccTOFL_afterslow2_Lmax;
MCNUM restore_neutron = mccTOFL_afterslow2_restore_neutron;
int nowritefile = mccTOFL_afterslow2_nowritefile;
#line 110 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/TOFLambda_monitor.comp"
{
    if (!nowritefile) {
    DETECTOR_OUT_2D(
        "TOF-wavelength monitor",
        "Time-of-flight [\\gms]", "Wavelength [AA]",
        tmin, tmax, Lmin, Lmax,
        nt, nL,
        &TOFL_N[0][0],&TOFL_p[0][0],&TOFL_p2[0][0],
        filename);
    }
}
#line 26017 "./ESS_IN5_reprate.c"
}   /* End of TOFL_afterslow2=TOFLambda_monitor() SETTING parameter declarations. */
#undef TOFL_p2
#undef TOFL_p
#undef TOFL_N
#undef tt_1
#undef tt_0
#undef tmax
#undef tmin
#undef nt
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'Lmon_beforeballistic'. */
  SIG_MESSAGE("Lmon_beforeballistic (Save)");
#define mccompcurname  Lmon_beforeballistic
#define mccompcurtype  L_monitor
#define mccompcurindex 23
#define nL mccLmon_beforeballistic_nL
#define L_N mccLmon_beforeballistic_L_N
#define L_p mccLmon_beforeballistic_L_p
#define L_p2 mccLmon_beforeballistic_L_p2
{   /* Declarations of Lmon_beforeballistic=L_monitor() SETTING parameters. */
char* filename = mccLmon_beforeballistic_filename;
MCNUM xmin = mccLmon_beforeballistic_xmin;
MCNUM xmax = mccLmon_beforeballistic_xmax;
MCNUM ymin = mccLmon_beforeballistic_ymin;
MCNUM ymax = mccLmon_beforeballistic_ymax;
MCNUM xwidth = mccLmon_beforeballistic_xwidth;
MCNUM yheight = mccLmon_beforeballistic_yheight;
MCNUM Lmin = mccLmon_beforeballistic_Lmin;
MCNUM Lmax = mccLmon_beforeballistic_Lmax;
MCNUM restore_neutron = mccLmon_beforeballistic_restore_neutron;
int nowritefile = mccLmon_beforeballistic_nowritefile;
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
#line 26065 "./ESS_IN5_reprate.c"
}   /* End of Lmon_beforeballistic=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'PSD_beforeballistic'. */
  SIG_MESSAGE("PSD_beforeballistic (Save)");
#define mccompcurname  PSD_beforeballistic
#define mccompcurtype  PSD_monitor
#define mccompcurindex 24
#define PSD_N mccPSD_beforeballistic_PSD_N
#define PSD_p mccPSD_beforeballistic_PSD_p
#define PSD_p2 mccPSD_beforeballistic_PSD_p2
{   /* Declarations of PSD_beforeballistic=PSD_monitor() SETTING parameters. */
int nx = mccPSD_beforeballistic_nx;
int ny = mccPSD_beforeballistic_ny;
char* filename = mccPSD_beforeballistic_filename;
MCNUM xmin = mccPSD_beforeballistic_xmin;
MCNUM xmax = mccPSD_beforeballistic_xmax;
MCNUM ymin = mccPSD_beforeballistic_ymin;
MCNUM ymax = mccPSD_beforeballistic_ymax;
MCNUM xwidth = mccPSD_beforeballistic_xwidth;
MCNUM yheight = mccPSD_beforeballistic_yheight;
MCNUM restore_neutron = mccPSD_beforeballistic_restore_neutron;
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
#line 26105 "./ESS_IN5_reprate.c"
}   /* End of PSD_beforeballistic=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'Lmonfast2'. */
  SIG_MESSAGE("Lmonfast2 (Save)");
#define mccompcurname  Lmonfast2
#define mccompcurtype  L_monitor
#define mccompcurindex 26
#define nL mccLmonfast2_nL
#define L_N mccLmonfast2_L_N
#define L_p mccLmonfast2_L_p
#define L_p2 mccLmonfast2_L_p2
{   /* Declarations of Lmonfast2=L_monitor() SETTING parameters. */
char* filename = mccLmonfast2_filename;
MCNUM xmin = mccLmonfast2_xmin;
MCNUM xmax = mccLmonfast2_xmax;
MCNUM ymin = mccLmonfast2_ymin;
MCNUM ymax = mccLmonfast2_ymax;
MCNUM xwidth = mccLmonfast2_xwidth;
MCNUM yheight = mccLmonfast2_yheight;
MCNUM Lmin = mccLmonfast2_Lmin;
MCNUM Lmax = mccLmonfast2_Lmax;
MCNUM restore_neutron = mccLmonfast2_restore_neutron;
int nowritefile = mccLmonfast2_nowritefile;
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
#line 26147 "./ESS_IN5_reprate.c"
}   /* End of Lmonfast2=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'Lmonfast2_zoom'. */
  SIG_MESSAGE("Lmonfast2_zoom (Save)");
#define mccompcurname  Lmonfast2_zoom
#define mccompcurtype  L_monitor
#define mccompcurindex 27
#define nL mccLmonfast2_zoom_nL
#define L_N mccLmonfast2_zoom_L_N
#define L_p mccLmonfast2_zoom_L_p
#define L_p2 mccLmonfast2_zoom_L_p2
{   /* Declarations of Lmonfast2_zoom=L_monitor() SETTING parameters. */
char* filename = mccLmonfast2_zoom_filename;
MCNUM xmin = mccLmonfast2_zoom_xmin;
MCNUM xmax = mccLmonfast2_zoom_xmax;
MCNUM ymin = mccLmonfast2_zoom_ymin;
MCNUM ymax = mccLmonfast2_zoom_ymax;
MCNUM xwidth = mccLmonfast2_zoom_xwidth;
MCNUM yheight = mccLmonfast2_zoom_yheight;
MCNUM Lmin = mccLmonfast2_zoom_Lmin;
MCNUM Lmax = mccLmonfast2_zoom_Lmax;
MCNUM restore_neutron = mccLmonfast2_zoom_restore_neutron;
int nowritefile = mccLmonfast2_zoom_nowritefile;
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
#line 26190 "./ESS_IN5_reprate.c"
}   /* End of Lmonfast2_zoom=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'TOFLfast2'. */
  SIG_MESSAGE("TOFLfast2 (Save)");
#define mccompcurname  TOFLfast2
#define mccompcurtype  TOFLambda_monitor
#define mccompcurindex 28
#define nL mccTOFLfast2_nL
#define nt mccTOFLfast2_nt
#define tmin mccTOFLfast2_tmin
#define tmax mccTOFLfast2_tmax
#define tt_0 mccTOFLfast2_tt_0
#define tt_1 mccTOFLfast2_tt_1
#define TOFL_N mccTOFLfast2_TOFL_N
#define TOFL_p mccTOFLfast2_TOFL_p
#define TOFL_p2 mccTOFLfast2_TOFL_p2
{   /* Declarations of TOFLfast2=TOFLambda_monitor() SETTING parameters. */
char* filename = mccTOFLfast2_filename;
MCNUM xmin = mccTOFLfast2_xmin;
MCNUM xmax = mccTOFLfast2_xmax;
MCNUM ymin = mccTOFLfast2_ymin;
MCNUM ymax = mccTOFLfast2_ymax;
MCNUM xwidth = mccTOFLfast2_xwidth;
MCNUM yheight = mccTOFLfast2_yheight;
MCNUM Lmin = mccTOFLfast2_Lmin;
MCNUM Lmax = mccTOFLfast2_Lmax;
MCNUM restore_neutron = mccTOFLfast2_restore_neutron;
int nowritefile = mccTOFLfast2_nowritefile;
#line 110 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/TOFLambda_monitor.comp"
{
    if (!nowritefile) {
    DETECTOR_OUT_2D(
        "TOF-wavelength monitor",
        "Time-of-flight [\\gms]", "Wavelength [AA]",
        tmin, tmax, Lmin, Lmax,
        nt, nL,
        &TOFL_N[0][0],&TOFL_p[0][0],&TOFL_p2[0][0],
        filename);
    }
}
#line 26238 "./ESS_IN5_reprate.c"
}   /* End of TOFLfast2=TOFLambda_monitor() SETTING parameter declarations. */
#undef TOFL_p2
#undef TOFL_p
#undef TOFL_N
#undef tt_1
#undef tt_0
#undef tmax
#undef tmin
#undef nt
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'TOFLfast2zoom'. */
  SIG_MESSAGE("TOFLfast2zoom (Save)");
#define mccompcurname  TOFLfast2zoom
#define mccompcurtype  TOFLambda_monitor
#define mccompcurindex 29
#define nL mccTOFLfast2zoom_nL
#define nt mccTOFLfast2zoom_nt
#define tmin mccTOFLfast2zoom_tmin
#define tmax mccTOFLfast2zoom_tmax
#define tt_0 mccTOFLfast2zoom_tt_0
#define tt_1 mccTOFLfast2zoom_tt_1
#define TOFL_N mccTOFLfast2zoom_TOFL_N
#define TOFL_p mccTOFLfast2zoom_TOFL_p
#define TOFL_p2 mccTOFLfast2zoom_TOFL_p2
{   /* Declarations of TOFLfast2zoom=TOFLambda_monitor() SETTING parameters. */
char* filename = mccTOFLfast2zoom_filename;
MCNUM xmin = mccTOFLfast2zoom_xmin;
MCNUM xmax = mccTOFLfast2zoom_xmax;
MCNUM ymin = mccTOFLfast2zoom_ymin;
MCNUM ymax = mccTOFLfast2zoom_ymax;
MCNUM xwidth = mccTOFLfast2zoom_xwidth;
MCNUM yheight = mccTOFLfast2zoom_yheight;
MCNUM Lmin = mccTOFLfast2zoom_Lmin;
MCNUM Lmax = mccTOFLfast2zoom_Lmax;
MCNUM restore_neutron = mccTOFLfast2zoom_restore_neutron;
int nowritefile = mccTOFLfast2zoom_nowritefile;
#line 110 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/TOFLambda_monitor.comp"
{
    if (!nowritefile) {
    DETECTOR_OUT_2D(
        "TOF-wavelength monitor",
        "Time-of-flight [\\gms]", "Wavelength [AA]",
        tmin, tmax, Lmin, Lmax,
        nt, nL,
        &TOFL_N[0][0],&TOFL_p[0][0],&TOFL_p2[0][0],
        filename);
    }
}
#line 26291 "./ESS_IN5_reprate.c"
}   /* End of TOFLfast2zoom=TOFLambda_monitor() SETTING parameter declarations. */
#undef TOFL_p2
#undef TOFL_p
#undef TOFL_N
#undef tt_1
#undef tt_0
#undef tmax
#undef tmin
#undef nt
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'PSDfast2'. */
  SIG_MESSAGE("PSDfast2 (Save)");
#define mccompcurname  PSDfast2
#define mccompcurtype  PSD_monitor
#define mccompcurindex 30
#define PSD_N mccPSDfast2_PSD_N
#define PSD_p mccPSDfast2_PSD_p
#define PSD_p2 mccPSDfast2_PSD_p2
{   /* Declarations of PSDfast2=PSD_monitor() SETTING parameters. */
int nx = mccPSDfast2_nx;
int ny = mccPSDfast2_ny;
char* filename = mccPSDfast2_filename;
MCNUM xmin = mccPSDfast2_xmin;
MCNUM xmax = mccPSDfast2_xmax;
MCNUM ymin = mccPSDfast2_ymin;
MCNUM ymax = mccPSDfast2_ymax;
MCNUM xwidth = mccPSDfast2_xwidth;
MCNUM yheight = mccPSDfast2_yheight;
MCNUM restore_neutron = mccPSDfast2_restore_neutron;
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
#line 26336 "./ESS_IN5_reprate.c"
}   /* End of PSDfast2=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'TOFfast2_zoom'. */
  SIG_MESSAGE("TOFfast2_zoom (Save)");
#define mccompcurname  TOFfast2_zoom
#define mccompcurtype  TOF_monitor
#define mccompcurindex 34
#define nt mccTOFfast2_zoom_nt
#define TOF_N mccTOFfast2_zoom_TOF_N
#define TOF_p mccTOFfast2_zoom_TOF_p
#define TOF_p2 mccTOFfast2_zoom_TOF_p2
#define t_min mccTOFfast2_zoom_t_min
#define t_max mccTOFfast2_zoom_t_max
#define delta_t mccTOFfast2_zoom_delta_t
{   /* Declarations of TOFfast2_zoom=TOF_monitor() SETTING parameters. */
char* filename = mccTOFfast2_zoom_filename;
MCNUM xmin = mccTOFfast2_zoom_xmin;
MCNUM xmax = mccTOFfast2_zoom_xmax;
MCNUM ymin = mccTOFfast2_zoom_ymin;
MCNUM ymax = mccTOFfast2_zoom_ymax;
MCNUM xwidth = mccTOFfast2_zoom_xwidth;
MCNUM yheight = mccTOFfast2_zoom_yheight;
MCNUM tmin = mccTOFfast2_zoom_tmin;
MCNUM tmax = mccTOFfast2_zoom_tmax;
MCNUM dt = mccTOFfast2_zoom_dt;
MCNUM restore_neutron = mccTOFfast2_zoom_restore_neutron;
int nowritefile = mccTOFfast2_zoom_nowritefile;
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
#line 26382 "./ESS_IN5_reprate.c"
}   /* End of TOFfast2_zoom=TOF_monitor() SETTING parameter declarations. */
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

  /* User SAVE code for component 'Lmon_afterfast2'. */
  SIG_MESSAGE("Lmon_afterfast2 (Save)");
#define mccompcurname  Lmon_afterfast2
#define mccompcurtype  L_monitor
#define mccompcurindex 35
#define nL mccLmon_afterfast2_nL
#define L_N mccLmon_afterfast2_L_N
#define L_p mccLmon_afterfast2_L_p
#define L_p2 mccLmon_afterfast2_L_p2
{   /* Declarations of Lmon_afterfast2=L_monitor() SETTING parameters. */
char* filename = mccLmon_afterfast2_filename;
MCNUM xmin = mccLmon_afterfast2_xmin;
MCNUM xmax = mccLmon_afterfast2_xmax;
MCNUM ymin = mccLmon_afterfast2_ymin;
MCNUM ymax = mccLmon_afterfast2_ymax;
MCNUM xwidth = mccLmon_afterfast2_xwidth;
MCNUM yheight = mccLmon_afterfast2_yheight;
MCNUM Lmin = mccLmon_afterfast2_Lmin;
MCNUM Lmax = mccLmon_afterfast2_Lmax;
MCNUM restore_neutron = mccLmon_afterfast2_restore_neutron;
int nowritefile = mccLmon_afterfast2_nowritefile;
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
#line 26428 "./ESS_IN5_reprate.c"
}   /* End of Lmon_afterfast2=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'TOFL_afterfast2'. */
  SIG_MESSAGE("TOFL_afterfast2 (Save)");
#define mccompcurname  TOFL_afterfast2
#define mccompcurtype  TOFLambda_monitor
#define mccompcurindex 36
#define nL mccTOFL_afterfast2_nL
#define nt mccTOFL_afterfast2_nt
#define tmin mccTOFL_afterfast2_tmin
#define tmax mccTOFL_afterfast2_tmax
#define tt_0 mccTOFL_afterfast2_tt_0
#define tt_1 mccTOFL_afterfast2_tt_1
#define TOFL_N mccTOFL_afterfast2_TOFL_N
#define TOFL_p mccTOFL_afterfast2_TOFL_p
#define TOFL_p2 mccTOFL_afterfast2_TOFL_p2
{   /* Declarations of TOFL_afterfast2=TOFLambda_monitor() SETTING parameters. */
char* filename = mccTOFL_afterfast2_filename;
MCNUM xmin = mccTOFL_afterfast2_xmin;
MCNUM xmax = mccTOFL_afterfast2_xmax;
MCNUM ymin = mccTOFL_afterfast2_ymin;
MCNUM ymax = mccTOFL_afterfast2_ymax;
MCNUM xwidth = mccTOFL_afterfast2_xwidth;
MCNUM yheight = mccTOFL_afterfast2_yheight;
MCNUM Lmin = mccTOFL_afterfast2_Lmin;
MCNUM Lmax = mccTOFL_afterfast2_Lmax;
MCNUM restore_neutron = mccTOFL_afterfast2_restore_neutron;
int nowritefile = mccTOFL_afterfast2_nowritefile;
#line 110 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/TOFLambda_monitor.comp"
{
    if (!nowritefile) {
    DETECTOR_OUT_2D(
        "TOF-wavelength monitor",
        "Time-of-flight [\\gms]", "Wavelength [AA]",
        tmin, tmax, Lmin, Lmax,
        nt, nL,
        &TOFL_N[0][0],&TOFL_p[0][0],&TOFL_p2[0][0],
        filename);
    }
}
#line 26476 "./ESS_IN5_reprate.c"
}   /* End of TOFL_afterfast2=TOFLambda_monitor() SETTING parameter declarations. */
#undef TOFL_p2
#undef TOFL_p
#undef TOFL_N
#undef tt_1
#undef tt_0
#undef tmax
#undef tmin
#undef nt
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'TOFL_afterfast2_zoom'. */
  SIG_MESSAGE("TOFL_afterfast2_zoom (Save)");
#define mccompcurname  TOFL_afterfast2_zoom
#define mccompcurtype  TOFLambda_monitor
#define mccompcurindex 37
#define nL mccTOFL_afterfast2_zoom_nL
#define nt mccTOFL_afterfast2_zoom_nt
#define tmin mccTOFL_afterfast2_zoom_tmin
#define tmax mccTOFL_afterfast2_zoom_tmax
#define tt_0 mccTOFL_afterfast2_zoom_tt_0
#define tt_1 mccTOFL_afterfast2_zoom_tt_1
#define TOFL_N mccTOFL_afterfast2_zoom_TOFL_N
#define TOFL_p mccTOFL_afterfast2_zoom_TOFL_p
#define TOFL_p2 mccTOFL_afterfast2_zoom_TOFL_p2
{   /* Declarations of TOFL_afterfast2_zoom=TOFLambda_monitor() SETTING parameters. */
char* filename = mccTOFL_afterfast2_zoom_filename;
MCNUM xmin = mccTOFL_afterfast2_zoom_xmin;
MCNUM xmax = mccTOFL_afterfast2_zoom_xmax;
MCNUM ymin = mccTOFL_afterfast2_zoom_ymin;
MCNUM ymax = mccTOFL_afterfast2_zoom_ymax;
MCNUM xwidth = mccTOFL_afterfast2_zoom_xwidth;
MCNUM yheight = mccTOFL_afterfast2_zoom_yheight;
MCNUM Lmin = mccTOFL_afterfast2_zoom_Lmin;
MCNUM Lmax = mccTOFL_afterfast2_zoom_Lmax;
MCNUM restore_neutron = mccTOFL_afterfast2_zoom_restore_neutron;
int nowritefile = mccTOFL_afterfast2_zoom_nowritefile;
#line 110 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/TOFLambda_monitor.comp"
{
    if (!nowritefile) {
    DETECTOR_OUT_2D(
        "TOF-wavelength monitor",
        "Time-of-flight [\\gms]", "Wavelength [AA]",
        tmin, tmax, Lmin, Lmax,
        nt, nL,
        &TOFL_N[0][0],&TOFL_p[0][0],&TOFL_p2[0][0],
        filename);
    }
}
#line 26529 "./ESS_IN5_reprate.c"
}   /* End of TOFL_afterfast2_zoom=TOFLambda_monitor() SETTING parameter declarations. */
#undef TOFL_p2
#undef TOFL_p
#undef TOFL_N
#undef tt_1
#undef tt_0
#undef tmax
#undef tmin
#undef nt
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'PSD_afterfast2'. */
  SIG_MESSAGE("PSD_afterfast2 (Save)");
#define mccompcurname  PSD_afterfast2
#define mccompcurtype  PSD_monitor
#define mccompcurindex 38
#define PSD_N mccPSD_afterfast2_PSD_N
#define PSD_p mccPSD_afterfast2_PSD_p
#define PSD_p2 mccPSD_afterfast2_PSD_p2
{   /* Declarations of PSD_afterfast2=PSD_monitor() SETTING parameters. */
int nx = mccPSD_afterfast2_nx;
int ny = mccPSD_afterfast2_ny;
char* filename = mccPSD_afterfast2_filename;
MCNUM xmin = mccPSD_afterfast2_xmin;
MCNUM xmax = mccPSD_afterfast2_xmax;
MCNUM ymin = mccPSD_afterfast2_ymin;
MCNUM ymax = mccPSD_afterfast2_ymax;
MCNUM xwidth = mccPSD_afterfast2_xwidth;
MCNUM yheight = mccPSD_afterfast2_yheight;
MCNUM restore_neutron = mccPSD_afterfast2_restore_neutron;
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
#line 26574 "./ESS_IN5_reprate.c"
}   /* End of PSD_afterfast2=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'Lmon_guideend'. */
  SIG_MESSAGE("Lmon_guideend (Save)");
#define mccompcurname  Lmon_guideend
#define mccompcurtype  L_monitor
#define mccompcurindex 40
#define nL mccLmon_guideend_nL
#define L_N mccLmon_guideend_L_N
#define L_p mccLmon_guideend_L_p
#define L_p2 mccLmon_guideend_L_p2
{   /* Declarations of Lmon_guideend=L_monitor() SETTING parameters. */
char* filename = mccLmon_guideend_filename;
MCNUM xmin = mccLmon_guideend_xmin;
MCNUM xmax = mccLmon_guideend_xmax;
MCNUM ymin = mccLmon_guideend_ymin;
MCNUM ymax = mccLmon_guideend_ymax;
MCNUM xwidth = mccLmon_guideend_xwidth;
MCNUM yheight = mccLmon_guideend_yheight;
MCNUM Lmin = mccLmon_guideend_Lmin;
MCNUM Lmax = mccLmon_guideend_Lmax;
MCNUM restore_neutron = mccLmon_guideend_restore_neutron;
int nowritefile = mccLmon_guideend_nowritefile;
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
#line 26616 "./ESS_IN5_reprate.c"
}   /* End of Lmon_guideend=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'PSDsample'. */
  SIG_MESSAGE("PSDsample (Save)");
#define mccompcurname  PSDsample
#define mccompcurtype  PSD_monitor
#define mccompcurindex 41
#define PSD_N mccPSDsample_PSD_N
#define PSD_p mccPSDsample_PSD_p
#define PSD_p2 mccPSDsample_PSD_p2
{   /* Declarations of PSDsample=PSD_monitor() SETTING parameters. */
int nx = mccPSDsample_nx;
int ny = mccPSDsample_ny;
char* filename = mccPSDsample_filename;
MCNUM xmin = mccPSDsample_xmin;
MCNUM xmax = mccPSDsample_xmax;
MCNUM ymin = mccPSDsample_ymin;
MCNUM ymax = mccPSDsample_ymax;
MCNUM xwidth = mccPSDsample_xwidth;
MCNUM yheight = mccPSDsample_yheight;
MCNUM restore_neutron = mccPSDsample_restore_neutron;
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
#line 26656 "./ESS_IN5_reprate.c"
}   /* End of PSDsample=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'TOFsample_zoom'. */
  SIG_MESSAGE("TOFsample_zoom (Save)");
#define mccompcurname  TOFsample_zoom
#define mccompcurtype  TOF_monitor
#define mccompcurindex 42
#define nt mccTOFsample_zoom_nt
#define TOF_N mccTOFsample_zoom_TOF_N
#define TOF_p mccTOFsample_zoom_TOF_p
#define TOF_p2 mccTOFsample_zoom_TOF_p2
#define t_min mccTOFsample_zoom_t_min
#define t_max mccTOFsample_zoom_t_max
#define delta_t mccTOFsample_zoom_delta_t
{   /* Declarations of TOFsample_zoom=TOF_monitor() SETTING parameters. */
char* filename = mccTOFsample_zoom_filename;
MCNUM xmin = mccTOFsample_zoom_xmin;
MCNUM xmax = mccTOFsample_zoom_xmax;
MCNUM ymin = mccTOFsample_zoom_ymin;
MCNUM ymax = mccTOFsample_zoom_ymax;
MCNUM xwidth = mccTOFsample_zoom_xwidth;
MCNUM yheight = mccTOFsample_zoom_yheight;
MCNUM tmin = mccTOFsample_zoom_tmin;
MCNUM tmax = mccTOFsample_zoom_tmax;
MCNUM dt = mccTOFsample_zoom_dt;
MCNUM restore_neutron = mccTOFsample_zoom_restore_neutron;
int nowritefile = mccTOFsample_zoom_nowritefile;
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
#line 26702 "./ESS_IN5_reprate.c"
}   /* End of TOFsample_zoom=TOF_monitor() SETTING parameter declarations. */
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

  /* User SAVE code for component 'Esample'. */
  SIG_MESSAGE("Esample (Save)");
#define mccompcurname  Esample
#define mccompcurtype  E_monitor
#define mccompcurindex 43
#define nE mccEsample_nE
#define E_N mccEsample_E_N
#define E_p mccEsample_E_p
#define E_p2 mccEsample_E_p2
#define S_p mccEsample_S_p
#define S_pE mccEsample_S_pE
#define S_pE2 mccEsample_S_pE2
{   /* Declarations of Esample=E_monitor() SETTING parameters. */
char* filename = mccEsample_filename;
MCNUM xmin = mccEsample_xmin;
MCNUM xmax = mccEsample_xmax;
MCNUM ymin = mccEsample_ymin;
MCNUM ymax = mccEsample_ymax;
MCNUM xwidth = mccEsample_xwidth;
MCNUM yheight = mccEsample_yheight;
MCNUM Emin = mccEsample_Emin;
MCNUM Emax = mccEsample_Emax;
MCNUM restore_neutron = mccEsample_restore_neutron;
int nowritefile = mccEsample_nowritefile;
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
#line 26753 "./ESS_IN5_reprate.c"
}   /* End of Esample=E_monitor() SETTING parameter declarations. */
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

  /* User SAVE code for component 'Lmon_sample_zoom'. */
  SIG_MESSAGE("Lmon_sample_zoom (Save)");
#define mccompcurname  Lmon_sample_zoom
#define mccompcurtype  L_monitor
#define mccompcurindex 44
#define nL mccLmon_sample_zoom_nL
#define L_N mccLmon_sample_zoom_L_N
#define L_p mccLmon_sample_zoom_L_p
#define L_p2 mccLmon_sample_zoom_L_p2
{   /* Declarations of Lmon_sample_zoom=L_monitor() SETTING parameters. */
char* filename = mccLmon_sample_zoom_filename;
MCNUM xmin = mccLmon_sample_zoom_xmin;
MCNUM xmax = mccLmon_sample_zoom_xmax;
MCNUM ymin = mccLmon_sample_zoom_ymin;
MCNUM ymax = mccLmon_sample_zoom_ymax;
MCNUM xwidth = mccLmon_sample_zoom_xwidth;
MCNUM yheight = mccLmon_sample_zoom_yheight;
MCNUM Lmin = mccLmon_sample_zoom_Lmin;
MCNUM Lmax = mccLmon_sample_zoom_Lmax;
MCNUM restore_neutron = mccLmon_sample_zoom_restore_neutron;
int nowritefile = mccLmon_sample_zoom_nowritefile;
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
#line 26799 "./ESS_IN5_reprate.c"
}   /* End of Lmon_sample_zoom=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'TOFdetector'. */
  SIG_MESSAGE("TOFdetector (Save)");
#define mccompcurname  TOFdetector
#define mccompcurtype  TOF_monitor
#define mccompcurindex 47
#define nt mccTOFdetector_nt
#define TOF_N mccTOFdetector_TOF_N
#define TOF_p mccTOFdetector_TOF_p
#define TOF_p2 mccTOFdetector_TOF_p2
#define t_min mccTOFdetector_t_min
#define t_max mccTOFdetector_t_max
#define delta_t mccTOFdetector_delta_t
{   /* Declarations of TOFdetector=TOF_monitor() SETTING parameters. */
char* filename = mccTOFdetector_filename;
MCNUM xmin = mccTOFdetector_xmin;
MCNUM xmax = mccTOFdetector_xmax;
MCNUM ymin = mccTOFdetector_ymin;
MCNUM ymax = mccTOFdetector_ymax;
MCNUM xwidth = mccTOFdetector_xwidth;
MCNUM yheight = mccTOFdetector_yheight;
MCNUM tmin = mccTOFdetector_tmin;
MCNUM tmax = mccTOFdetector_tmax;
MCNUM dt = mccTOFdetector_dt;
MCNUM restore_neutron = mccTOFdetector_restore_neutron;
int nowritefile = mccTOFdetector_nowritefile;
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
#line 26846 "./ESS_IN5_reprate.c"
}   /* End of TOFdetector=TOF_monitor() SETTING parameter declarations. */
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

  /* User SAVE code for component 'TOFdetector_zoom'. */
  SIG_MESSAGE("TOFdetector_zoom (Save)");
#define mccompcurname  TOFdetector_zoom
#define mccompcurtype  TOF_monitor
#define mccompcurindex 48
#define nt mccTOFdetector_zoom_nt
#define TOF_N mccTOFdetector_zoom_TOF_N
#define TOF_p mccTOFdetector_zoom_TOF_p
#define TOF_p2 mccTOFdetector_zoom_TOF_p2
#define t_min mccTOFdetector_zoom_t_min
#define t_max mccTOFdetector_zoom_t_max
#define delta_t mccTOFdetector_zoom_delta_t
{   /* Declarations of TOFdetector_zoom=TOF_monitor() SETTING parameters. */
char* filename = mccTOFdetector_zoom_filename;
MCNUM xmin = mccTOFdetector_zoom_xmin;
MCNUM xmax = mccTOFdetector_zoom_xmax;
MCNUM ymin = mccTOFdetector_zoom_ymin;
MCNUM ymax = mccTOFdetector_zoom_ymax;
MCNUM xwidth = mccTOFdetector_zoom_xwidth;
MCNUM yheight = mccTOFdetector_zoom_yheight;
MCNUM tmin = mccTOFdetector_zoom_tmin;
MCNUM tmax = mccTOFdetector_zoom_tmax;
MCNUM dt = mccTOFdetector_zoom_dt;
MCNUM restore_neutron = mccTOFdetector_zoom_restore_neutron;
int nowritefile = mccTOFdetector_zoom_nowritefile;
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
#line 26896 "./ESS_IN5_reprate.c"
}   /* End of TOFdetector_zoom=TOF_monitor() SETTING parameter declarations. */
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

  /* User SAVE code for component 'Edetector'. */
  SIG_MESSAGE("Edetector (Save)");
#define mccompcurname  Edetector
#define mccompcurtype  E_monitor
#define mccompcurindex 49
#define nE mccEdetector_nE
#define E_N mccEdetector_E_N
#define E_p mccEdetector_E_p
#define E_p2 mccEdetector_E_p2
#define S_p mccEdetector_S_p
#define S_pE mccEdetector_S_pE
#define S_pE2 mccEdetector_S_pE2
{   /* Declarations of Edetector=E_monitor() SETTING parameters. */
char* filename = mccEdetector_filename;
MCNUM xmin = mccEdetector_xmin;
MCNUM xmax = mccEdetector_xmax;
MCNUM ymin = mccEdetector_ymin;
MCNUM ymax = mccEdetector_ymax;
MCNUM xwidth = mccEdetector_xwidth;
MCNUM yheight = mccEdetector_yheight;
MCNUM Emin = mccEdetector_Emin;
MCNUM Emax = mccEdetector_Emax;
MCNUM restore_neutron = mccEdetector_restore_neutron;
int nowritefile = mccEdetector_nowritefile;
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
#line 26947 "./ESS_IN5_reprate.c"
}   /* End of Edetector=E_monitor() SETTING parameter declarations. */
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

  /* User SAVE code for component 'TOF2Edetector'. */
  SIG_MESSAGE("TOF2Edetector (Save)");
#define mccompcurname  TOF2Edetector
#define mccompcurtype  TOF2E_monitor
#define mccompcurindex 50
#define nE mccTOF2Edetector_nE
#define E_N mccTOF2Edetector_E_N
#define E_p mccTOF2Edetector_E_p
#define E_p2 mccTOF2Edetector_E_p2
#define S_p mccTOF2Edetector_S_p
#define S_pE mccTOF2Edetector_S_pE
#define S_pE2 mccTOF2Edetector_S_pE2
{   /* Declarations of TOF2Edetector=TOF2E_monitor() SETTING parameters. */
char* filename = mccTOF2Edetector_filename;
MCNUM xmin = mccTOF2Edetector_xmin;
MCNUM xmax = mccTOF2Edetector_xmax;
MCNUM ymin = mccTOF2Edetector_ymin;
MCNUM ymax = mccTOF2Edetector_ymax;
MCNUM xwidth = mccTOF2Edetector_xwidth;
MCNUM yheight = mccTOF2Edetector_yheight;
MCNUM Emin = mccTOF2Edetector_Emin;
MCNUM Emax = mccTOF2Edetector_Emax;
MCNUM T_zero = mccTOF2Edetector_T_zero;
MCNUM L_flight = mccTOF2Edetector_L_flight;
MCNUM restore_neutron = mccTOF2Edetector_restore_neutron;
int nowritefile = mccTOF2Edetector_nowritefile;
#line 118 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/TOF2E_monitor.comp"
{
    if (!nowritefile) {
    DETECTOR_OUT_1D(
        "TOF-to-Energy monitor",
        "Energy [meV]",
        "Intensity",
        "E", Emin, Emax, nE,
        &E_N[0],&E_p[0],&E_p2[0],
        filename);
    if (S_p) printf("<E> : %g meV , E-width : %g meV \n",
     S_pE/S_p,sqrt(S_pE2/S_p - S_pE*S_pE/(S_p*S_p)) );
    }
}
#line 27000 "./ESS_IN5_reprate.c"
}   /* End of TOF2Edetector=TOF2E_monitor() SETTING parameter declarations. */
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

  if (!handle) mcsiminfo_close(); 
} /* end save */
void mcfinally(void) {
  /* User component FINALLY code. */
  mcsiminfo_init(NULL);
  mcsave(mcsiminfo_file); /* save data when simulation ends */

    if (!mcNCounter[1]) fprintf(stderr, "Warning: No neutron could reach Component[1] source\n");
    if (mcAbsorbProp[1]) fprintf(stderr, "Warning: %g events were removed in Component[1] source=ESS_moderator_long()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[1]);
  /* User FINALLY code for component 'Origin'. */
  SIG_MESSAGE("Origin (Finally)");
#define mccompcurname  Origin
#define mccompcurtype  Progress_bar
#define mccompcurindex 2
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
#line 27049 "./ESS_IN5_reprate.c"
}   /* End of Origin=Progress_bar() SETTING parameter declarations. */
#undef CurrentTime
#undef EndTime
#undef StartTime
#undef IntermediateCnts
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[2]) fprintf(stderr, "Warning: No neutron could reach Component[2] Origin\n");
    if (mcAbsorbProp[2]) fprintf(stderr, "Warning: %g events were removed in Component[2] Origin=Progress_bar()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[2]);
    if (!mcNCounter[3]) fprintf(stderr, "Warning: No neutron could reach Component[3] TOFmoderator_zoom\n");
    if (mcAbsorbProp[3]) fprintf(stderr, "Warning: %g events were removed in Component[3] TOFmoderator_zoom=TOF_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[3]);
    if (!mcNCounter[4]) fprintf(stderr, "Warning: No neutron could reach Component[4] TOFmoderator\n");
    if (mcAbsorbProp[4]) fprintf(stderr, "Warning: %g events were removed in Component[4] TOFmoderator=TOF_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[4]);
    if (!mcNCounter[5]) fprintf(stderr, "Warning: No neutron could reach Component[5] Lmon_guistart\n");
    if (mcAbsorbProp[5]) fprintf(stderr, "Warning: %g events were removed in Component[5] Lmon_guistart=L_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[5]);
    if (!mcNCounter[6]) fprintf(stderr, "Warning: No neutron could reach Component[6] Lmon_normalize\n");
    if (mcAbsorbProp[6]) fprintf(stderr, "Warning: %g events were removed in Component[6] Lmon_normalize=L_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[6]);
    if (!mcNCounter[7]) fprintf(stderr, "Warning: No neutron could reach Component[7] Guide1\n");
    if (mcAbsorbProp[7]) fprintf(stderr, "Warning: %g events were removed in Component[7] Guide1=Guide()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[7]);
    if (!mcNCounter[8]) fprintf(stderr, "Warning: No neutron could reach Component[8] Lmonslow1\n");
    if (mcAbsorbProp[8]) fprintf(stderr, "Warning: %g events were removed in Component[8] Lmonslow1=L_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[8]);
  /* User FINALLY code for component 'PSDslow1'. */
  SIG_MESSAGE("PSDslow1 (Finally)");
#define mccompcurname  PSDslow1
#define mccompcurtype  PSD_monitor
#define mccompcurindex 9
#define PSD_N mccPSDslow1_PSD_N
#define PSD_p mccPSDslow1_PSD_p
#define PSD_p2 mccPSDslow1_PSD_p2
{   /* Declarations of PSDslow1=PSD_monitor() SETTING parameters. */
int nx = mccPSDslow1_nx;
int ny = mccPSDslow1_ny;
char* filename = mccPSDslow1_filename;
MCNUM xmin = mccPSDslow1_xmin;
MCNUM xmax = mccPSDslow1_xmax;
MCNUM ymin = mccPSDslow1_ymin;
MCNUM ymax = mccPSDslow1_ymax;
MCNUM xwidth = mccPSDslow1_xwidth;
MCNUM yheight = mccPSDslow1_yheight;
MCNUM restore_neutron = mccPSDslow1_restore_neutron;
#line 122 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
{
  destroy_darr2d(PSD_N);
  destroy_darr2d(PSD_p);
  destroy_darr2d(PSD_p2);
}
#line 27098 "./ESS_IN5_reprate.c"
}   /* End of PSDslow1=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[9]) fprintf(stderr, "Warning: No neutron could reach Component[9] PSDslow1\n");
    if (mcAbsorbProp[9]) fprintf(stderr, "Warning: %g events were removed in Component[9] PSDslow1=PSD_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[9]);
    if (!mcNCounter[10]) fprintf(stderr, "Warning: No neutron could reach Component[10] FOchop1\n");
    if (mcAbsorbProp[10]) fprintf(stderr, "Warning: %g events were removed in Component[10] FOchop1=DiskChopper()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[10]);
    if (!mcNCounter[11]) fprintf(stderr, "Warning: No neutron could reach Component[11] TOFLmon1\n");
    if (mcAbsorbProp[11]) fprintf(stderr, "Warning: %g events were removed in Component[11] TOFLmon1=TOFLambda_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[11]);
    if (!mcNCounter[12]) fprintf(stderr, "Warning: No neutron could reach Component[12] Lmon_afterslow1\n");
    if (mcAbsorbProp[12]) fprintf(stderr, "Warning: %g events were removed in Component[12] Lmon_afterslow1=L_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[12]);
  /* User FINALLY code for component 'PSD_afterslow1'. */
  SIG_MESSAGE("PSD_afterslow1 (Finally)");
#define mccompcurname  PSD_afterslow1
#define mccompcurtype  PSD_monitor
#define mccompcurindex 13
#define PSD_N mccPSD_afterslow1_PSD_N
#define PSD_p mccPSD_afterslow1_PSD_p
#define PSD_p2 mccPSD_afterslow1_PSD_p2
{   /* Declarations of PSD_afterslow1=PSD_monitor() SETTING parameters. */
int nx = mccPSD_afterslow1_nx;
int ny = mccPSD_afterslow1_ny;
char* filename = mccPSD_afterslow1_filename;
MCNUM xmin = mccPSD_afterslow1_xmin;
MCNUM xmax = mccPSD_afterslow1_xmax;
MCNUM ymin = mccPSD_afterslow1_ymin;
MCNUM ymax = mccPSD_afterslow1_ymax;
MCNUM xwidth = mccPSD_afterslow1_xwidth;
MCNUM yheight = mccPSD_afterslow1_yheight;
MCNUM restore_neutron = mccPSD_afterslow1_restore_neutron;
#line 122 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
{
  destroy_darr2d(PSD_N);
  destroy_darr2d(PSD_p);
  destroy_darr2d(PSD_p2);
}
#line 27140 "./ESS_IN5_reprate.c"
}   /* End of PSD_afterslow1=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[13]) fprintf(stderr, "Warning: No neutron could reach Component[13] PSD_afterslow1\n");
    if (mcAbsorbProp[13]) fprintf(stderr, "Warning: %g events were removed in Component[13] PSD_afterslow1=PSD_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[13]);
    if (!mcNCounter[14]) fprintf(stderr, "Warning: No neutron could reach Component[14] Guidelong1\n");
    if (mcAbsorbProp[14]) fprintf(stderr, "Warning: %g events were removed in Component[14] Guidelong1=Guide()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[14]);
    if (!mcNCounter[15]) fprintf(stderr, "Warning: No neutron could reach Component[15] Guidelong1b\n");
    if (mcAbsorbProp[15]) fprintf(stderr, "Warning: %g events were removed in Component[15] Guidelong1b=Guide()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[15]);
    if (!mcNCounter[16]) fprintf(stderr, "Warning: No neutron could reach Component[16] Lmon_slow2\n");
    if (mcAbsorbProp[16]) fprintf(stderr, "Warning: %g events were removed in Component[16] Lmon_slow2=L_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[16]);
    if (!mcNCounter[17]) fprintf(stderr, "Warning: No neutron could reach Component[17] FOchop2\n");
    if (mcAbsorbProp[17]) fprintf(stderr, "Warning: %g events were removed in Component[17] FOchop2=DiskChopper()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[17]);
    if (!mcNCounter[18]) fprintf(stderr, "Warning: No neutron could reach Component[18] Fastchop1\n");
    if (mcAbsorbProp[18]) fprintf(stderr, "Warning: %g events were removed in Component[18] Fastchop1=DiskChopper()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[18]);
  /* User FINALLY code for component 'PSD_afterslow2'. */
  SIG_MESSAGE("PSD_afterslow2 (Finally)");
#define mccompcurname  PSD_afterslow2
#define mccompcurtype  PSD_monitor
#define mccompcurindex 19
#define PSD_N mccPSD_afterslow2_PSD_N
#define PSD_p mccPSD_afterslow2_PSD_p
#define PSD_p2 mccPSD_afterslow2_PSD_p2
{   /* Declarations of PSD_afterslow2=PSD_monitor() SETTING parameters. */
int nx = mccPSD_afterslow2_nx;
int ny = mccPSD_afterslow2_ny;
char* filename = mccPSD_afterslow2_filename;
MCNUM xmin = mccPSD_afterslow2_xmin;
MCNUM xmax = mccPSD_afterslow2_xmax;
MCNUM ymin = mccPSD_afterslow2_ymin;
MCNUM ymax = mccPSD_afterslow2_ymax;
MCNUM xwidth = mccPSD_afterslow2_xwidth;
MCNUM yheight = mccPSD_afterslow2_yheight;
MCNUM restore_neutron = mccPSD_afterslow2_restore_neutron;
#line 122 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
{
  destroy_darr2d(PSD_N);
  destroy_darr2d(PSD_p);
  destroy_darr2d(PSD_p2);
}
#line 27186 "./ESS_IN5_reprate.c"
}   /* End of PSD_afterslow2=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[19]) fprintf(stderr, "Warning: No neutron could reach Component[19] PSD_afterslow2\n");
    if (mcAbsorbProp[19]) fprintf(stderr, "Warning: %g events were removed in Component[19] PSD_afterslow2=PSD_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[19]);
    if (!mcNCounter[20]) fprintf(stderr, "Warning: No neutron could reach Component[20] Lmon_afterslow2\n");
    if (mcAbsorbProp[20]) fprintf(stderr, "Warning: %g events were removed in Component[20] Lmon_afterslow2=L_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[20]);
    if (!mcNCounter[21]) fprintf(stderr, "Warning: No neutron could reach Component[21] TOFL_afterslow2\n");
    if (mcAbsorbProp[21]) fprintf(stderr, "Warning: %g events were removed in Component[21] TOFL_afterslow2=TOFLambda_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[21]);
    if (!mcNCounter[22]) fprintf(stderr, "Warning: No neutron could reach Component[22] Guidelong2\n");
    if (mcAbsorbProp[22]) fprintf(stderr, "Warning: %g events were removed in Component[22] Guidelong2=Guide()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[22]);
    if (!mcNCounter[23]) fprintf(stderr, "Warning: No neutron could reach Component[23] Lmon_beforeballistic\n");
    if (mcAbsorbProp[23]) fprintf(stderr, "Warning: %g events were removed in Component[23] Lmon_beforeballistic=L_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[23]);
  /* User FINALLY code for component 'PSD_beforeballistic'. */
  SIG_MESSAGE("PSD_beforeballistic (Finally)");
#define mccompcurname  PSD_beforeballistic
#define mccompcurtype  PSD_monitor
#define mccompcurindex 24
#define PSD_N mccPSD_beforeballistic_PSD_N
#define PSD_p mccPSD_beforeballistic_PSD_p
#define PSD_p2 mccPSD_beforeballistic_PSD_p2
{   /* Declarations of PSD_beforeballistic=PSD_monitor() SETTING parameters. */
int nx = mccPSD_beforeballistic_nx;
int ny = mccPSD_beforeballistic_ny;
char* filename = mccPSD_beforeballistic_filename;
MCNUM xmin = mccPSD_beforeballistic_xmin;
MCNUM xmax = mccPSD_beforeballistic_xmax;
MCNUM ymin = mccPSD_beforeballistic_ymin;
MCNUM ymax = mccPSD_beforeballistic_ymax;
MCNUM xwidth = mccPSD_beforeballistic_xwidth;
MCNUM yheight = mccPSD_beforeballistic_yheight;
MCNUM restore_neutron = mccPSD_beforeballistic_restore_neutron;
#line 122 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
{
  destroy_darr2d(PSD_N);
  destroy_darr2d(PSD_p);
  destroy_darr2d(PSD_p2);
}
#line 27230 "./ESS_IN5_reprate.c"
}   /* End of PSD_beforeballistic=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[24]) fprintf(stderr, "Warning: No neutron could reach Component[24] PSD_beforeballistic\n");
    if (mcAbsorbProp[24]) fprintf(stderr, "Warning: %g events were removed in Component[24] PSD_beforeballistic=PSD_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[24]);
    if (!mcNCounter[25]) fprintf(stderr, "Warning: No neutron could reach Component[25] Guidelong2a\n");
    if (mcAbsorbProp[25]) fprintf(stderr, "Warning: %g events were removed in Component[25] Guidelong2a=Guide()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[25]);
    if (!mcNCounter[26]) fprintf(stderr, "Warning: No neutron could reach Component[26] Lmonfast2\n");
    if (mcAbsorbProp[26]) fprintf(stderr, "Warning: %g events were removed in Component[26] Lmonfast2=L_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[26]);
    if (!mcNCounter[27]) fprintf(stderr, "Warning: No neutron could reach Component[27] Lmonfast2_zoom\n");
    if (mcAbsorbProp[27]) fprintf(stderr, "Warning: %g events were removed in Component[27] Lmonfast2_zoom=L_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[27]);
    if (!mcNCounter[28]) fprintf(stderr, "Warning: No neutron could reach Component[28] TOFLfast2\n");
    if (mcAbsorbProp[28]) fprintf(stderr, "Warning: %g events were removed in Component[28] TOFLfast2=TOFLambda_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[28]);
    if (!mcNCounter[29]) fprintf(stderr, "Warning: No neutron could reach Component[29] TOFLfast2zoom\n");
    if (mcAbsorbProp[29]) fprintf(stderr, "Warning: %g events were removed in Component[29] TOFLfast2zoom=TOFLambda_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[29]);
  /* User FINALLY code for component 'PSDfast2'. */
  SIG_MESSAGE("PSDfast2 (Finally)");
#define mccompcurname  PSDfast2
#define mccompcurtype  PSD_monitor
#define mccompcurindex 30
#define PSD_N mccPSDfast2_PSD_N
#define PSD_p mccPSDfast2_PSD_p
#define PSD_p2 mccPSDfast2_PSD_p2
{   /* Declarations of PSDfast2=PSD_monitor() SETTING parameters. */
int nx = mccPSDfast2_nx;
int ny = mccPSDfast2_ny;
char* filename = mccPSDfast2_filename;
MCNUM xmin = mccPSDfast2_xmin;
MCNUM xmax = mccPSDfast2_xmax;
MCNUM ymin = mccPSDfast2_ymin;
MCNUM ymax = mccPSDfast2_ymax;
MCNUM xwidth = mccPSDfast2_xwidth;
MCNUM yheight = mccPSDfast2_yheight;
MCNUM restore_neutron = mccPSDfast2_restore_neutron;
#line 122 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
{
  destroy_darr2d(PSD_N);
  destroy_darr2d(PSD_p);
  destroy_darr2d(PSD_p2);
}
#line 27276 "./ESS_IN5_reprate.c"
}   /* End of PSDfast2=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[30]) fprintf(stderr, "Warning: No neutron could reach Component[30] PSDfast2\n");
    if (mcAbsorbProp[30]) fprintf(stderr, "Warning: %g events were removed in Component[30] PSDfast2=PSD_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[30]);
    if (!mcNCounter[31]) fprintf(stderr, "Warning: No neutron could reach Component[31] Fastchop2\n");
    if (mcAbsorbProp[31]) fprintf(stderr, "Warning: %g events were removed in Component[31] Fastchop2=DiskChopper()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[31]);
    if (!mcNCounter[32]) fprintf(stderr, "Warning: No neutron could reach Component[32] Fastchop2counter\n");
    if (mcAbsorbProp[32]) fprintf(stderr, "Warning: %g events were removed in Component[32] Fastchop2counter=DiskChopper()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[32]);
    if (!mcNCounter[33]) fprintf(stderr, "Warning: No neutron could reach Component[33] FOchop3\n");
    if (mcAbsorbProp[33]) fprintf(stderr, "Warning: %g events were removed in Component[33] FOchop3=DiskChopper()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[33]);
    if (!mcNCounter[34]) fprintf(stderr, "Warning: No neutron could reach Component[34] TOFfast2_zoom\n");
    if (mcAbsorbProp[34]) fprintf(stderr, "Warning: %g events were removed in Component[34] TOFfast2_zoom=TOF_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[34]);
    if (!mcNCounter[35]) fprintf(stderr, "Warning: No neutron could reach Component[35] Lmon_afterfast2\n");
    if (mcAbsorbProp[35]) fprintf(stderr, "Warning: %g events were removed in Component[35] Lmon_afterfast2=L_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[35]);
    if (!mcNCounter[36]) fprintf(stderr, "Warning: No neutron could reach Component[36] TOFL_afterfast2\n");
    if (mcAbsorbProp[36]) fprintf(stderr, "Warning: %g events were removed in Component[36] TOFL_afterfast2=TOFLambda_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[36]);
    if (!mcNCounter[37]) fprintf(stderr, "Warning: No neutron could reach Component[37] TOFL_afterfast2_zoom\n");
    if (mcAbsorbProp[37]) fprintf(stderr, "Warning: %g events were removed in Component[37] TOFL_afterfast2_zoom=TOFLambda_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[37]);
  /* User FINALLY code for component 'PSD_afterfast2'. */
  SIG_MESSAGE("PSD_afterfast2 (Finally)");
#define mccompcurname  PSD_afterfast2
#define mccompcurtype  PSD_monitor
#define mccompcurindex 38
#define PSD_N mccPSD_afterfast2_PSD_N
#define PSD_p mccPSD_afterfast2_PSD_p
#define PSD_p2 mccPSD_afterfast2_PSD_p2
{   /* Declarations of PSD_afterfast2=PSD_monitor() SETTING parameters. */
int nx = mccPSD_afterfast2_nx;
int ny = mccPSD_afterfast2_ny;
char* filename = mccPSD_afterfast2_filename;
MCNUM xmin = mccPSD_afterfast2_xmin;
MCNUM xmax = mccPSD_afterfast2_xmax;
MCNUM ymin = mccPSD_afterfast2_ymin;
MCNUM ymax = mccPSD_afterfast2_ymax;
MCNUM xwidth = mccPSD_afterfast2_xwidth;
MCNUM yheight = mccPSD_afterfast2_yheight;
MCNUM restore_neutron = mccPSD_afterfast2_restore_neutron;
#line 122 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
{
  destroy_darr2d(PSD_N);
  destroy_darr2d(PSD_p);
  destroy_darr2d(PSD_p2);
}
#line 27326 "./ESS_IN5_reprate.c"
}   /* End of PSD_afterfast2=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[38]) fprintf(stderr, "Warning: No neutron could reach Component[38] PSD_afterfast2\n");
    if (mcAbsorbProp[38]) fprintf(stderr, "Warning: %g events were removed in Component[38] PSD_afterfast2=PSD_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[38]);
    if (!mcNCounter[39]) fprintf(stderr, "Warning: No neutron could reach Component[39] Guidesample\n");
    if (mcAbsorbProp[39]) fprintf(stderr, "Warning: %g events were removed in Component[39] Guidesample=Guide()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[39]);
    if (!mcNCounter[40]) fprintf(stderr, "Warning: No neutron could reach Component[40] Lmon_guideend\n");
    if (mcAbsorbProp[40]) fprintf(stderr, "Warning: %g events were removed in Component[40] Lmon_guideend=L_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[40]);
  /* User FINALLY code for component 'PSDsample'. */
  SIG_MESSAGE("PSDsample (Finally)");
#define mccompcurname  PSDsample
#define mccompcurtype  PSD_monitor
#define mccompcurindex 41
#define PSD_N mccPSDsample_PSD_N
#define PSD_p mccPSDsample_PSD_p
#define PSD_p2 mccPSDsample_PSD_p2
{   /* Declarations of PSDsample=PSD_monitor() SETTING parameters. */
int nx = mccPSDsample_nx;
int ny = mccPSDsample_ny;
char* filename = mccPSDsample_filename;
MCNUM xmin = mccPSDsample_xmin;
MCNUM xmax = mccPSDsample_xmax;
MCNUM ymin = mccPSDsample_ymin;
MCNUM ymax = mccPSDsample_ymax;
MCNUM xwidth = mccPSDsample_xwidth;
MCNUM yheight = mccPSDsample_yheight;
MCNUM restore_neutron = mccPSDsample_restore_neutron;
#line 122 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
{
  destroy_darr2d(PSD_N);
  destroy_darr2d(PSD_p);
  destroy_darr2d(PSD_p2);
}
#line 27366 "./ESS_IN5_reprate.c"
}   /* End of PSDsample=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[41]) fprintf(stderr, "Warning: No neutron could reach Component[41] PSDsample\n");
    if (mcAbsorbProp[41]) fprintf(stderr, "Warning: %g events were removed in Component[41] PSDsample=PSD_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[41]);
    if (!mcNCounter[42]) fprintf(stderr, "Warning: No neutron could reach Component[42] TOFsample_zoom\n");
    if (mcAbsorbProp[42]) fprintf(stderr, "Warning: %g events were removed in Component[42] TOFsample_zoom=TOF_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[42]);
    if (!mcNCounter[43]) fprintf(stderr, "Warning: No neutron could reach Component[43] Esample\n");
    if (mcAbsorbProp[43]) fprintf(stderr, "Warning: %g events were removed in Component[43] Esample=E_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[43]);
    if (!mcNCounter[44]) fprintf(stderr, "Warning: No neutron could reach Component[44] Lmon_sample_zoom\n");
    if (mcAbsorbProp[44]) fprintf(stderr, "Warning: %g events were removed in Component[44] Lmon_sample_zoom=L_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[44]);
    if (!mcNCounter[45]) fprintf(stderr, "Warning: No neutron could reach Component[45] sample\n");
    if (mcAbsorbProp[45]) fprintf(stderr, "Warning: %g events were removed in Component[45] sample=Tunneling_sample()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[45]);
    if (!mcNCounter[46]) fprintf(stderr, "Warning: No neutron could reach Component[46] detectorarm\n");
    if (mcAbsorbProp[46]) fprintf(stderr, "Warning: %g events were removed in Component[46] detectorarm=Arm()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[46]);
    if (!mcNCounter[47]) fprintf(stderr, "Warning: No neutron could reach Component[47] TOFdetector\n");
    if (mcAbsorbProp[47]) fprintf(stderr, "Warning: %g events were removed in Component[47] TOFdetector=TOF_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[47]);
    if (!mcNCounter[48]) fprintf(stderr, "Warning: No neutron could reach Component[48] TOFdetector_zoom\n");
    if (mcAbsorbProp[48]) fprintf(stderr, "Warning: %g events were removed in Component[48] TOFdetector_zoom=TOF_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[48]);
    if (!mcNCounter[49]) fprintf(stderr, "Warning: No neutron could reach Component[49] Edetector\n");
    if (mcAbsorbProp[49]) fprintf(stderr, "Warning: %g events were removed in Component[49] Edetector=E_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[49]);
    if (!mcNCounter[50]) fprintf(stderr, "Warning: No neutron could reach Component[50] TOF2Edetector\n");
    if (mcAbsorbProp[50]) fprintf(stderr, "Warning: %g events were removed in Component[50] TOF2Edetector=TOF2E_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[50]);
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

  /* MCDISPLAY code for component 'source'. */
  SIG_MESSAGE("source (McDisplay)");
  printf("MCDISPLAY: component %s\n", "source");
#define mccompcurname  source
#define mccompcurtype  ESS_moderator_long
#define mccompcurindex 1
#define l_range mccsource_l_range
#define w_mult mccsource_w_mult
#define w_geom mccsource_w_geom
#define w_geom_c mccsource_w_geom_c
#define w_geom_t mccsource_w_geom_t
#define tx mccsource_tx
#define ty mccsource_ty
#define tz mccsource_tz
#define t1x mccsource_t1x
#define t1y mccsource_t1y
#define t1z mccsource_t1z
#define t2x mccsource_t2x
#define t2y mccsource_t2y
#define t2z mccsource_t2z
#define T_n mccsource_T_n
#define tau_n mccsource_tau_n
#define tau1_n mccsource_tau1_n
#define tau2_n mccsource_tau2_n
#define chi2_n mccsource_chi2_n
#define I0_n mccsource_I0_n
#define I2_n mccsource_I2_n
#define branch1_n mccsource_branch1_n
#define branch2_n mccsource_branch2_n
#define r_empty mccsource_r_empty
#define r_optics mccsource_r_optics
{   /* Declarations of source=ESS_moderator_long() SETTING parameters. */
MCNUM width_c = mccsource_width_c;
MCNUM yheight = mccsource_yheight;
MCNUM Lmin = mccsource_Lmin;
MCNUM Lmax = mccsource_Lmax;
MCNUM dist = mccsource_dist;
MCNUM focus_xw = mccsource_focus_xw;
MCNUM focus_yh = mccsource_focus_yh;
MCNUM nu = mccsource_nu;
MCNUM T = mccsource_T;
MCNUM tau = mccsource_tau;
MCNUM tau1 = mccsource_tau1;
MCNUM tau2 = mccsource_tau2;
MCNUM d = mccsource_d;
MCNUM n = mccsource_n;
MCNUM cold_frac = mccsource_cold_frac;
MCNUM n2 = mccsource_n2;
MCNUM chi2 = mccsource_chi2;
MCNUM I0 = mccsource_I0;
MCNUM I2 = mccsource_I2;
int target_index = mccsource_target_index;
MCNUM cyl_radius = mccsource_cyl_radius;
MCNUM branch1 = mccsource_branch1;
MCNUM branch2 = mccsource_branch2;
MCNUM branch_tail = mccsource_branch_tail;
int n_pulses = mccsource_n_pulses;
MCNUM width_t = mccsource_width_t;
MCNUM T_t = mccsource_T_t;
MCNUM tau_t = mccsource_tau_t;
MCNUM tau1_t = mccsource_tau1_t;
MCNUM tau2_t = mccsource_tau2_t;
MCNUM chi2_t = mccsource_chi2_t;
MCNUM I0_t = mccsource_I0_t;
MCNUM I2_t = mccsource_I2_t;
MCNUM branch1_t = mccsource_branch1_t;
MCNUM branch2_t = mccsource_branch2_t;
int src_2012 = mccsource_src_2012;
MCNUM tfocus_dist = mccsource_tfocus_dist;
MCNUM tfocus_time = mccsource_tfocus_time;
MCNUM tfocus_width = mccsource_tfocus_width;
MCNUM beamport_angle = mccsource_beamport_angle;
#line 437 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../obsolete/ESS_moderator_long.comp"
{
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
}
#line 27523 "./ESS_IN5_reprate.c"
}   /* End of source=ESS_moderator_long() SETTING parameter declarations. */
#undef r_optics
#undef r_empty
#undef branch2_n
#undef branch1_n
#undef I2_n
#undef I0_n
#undef chi2_n
#undef tau2_n
#undef tau1_n
#undef tau_n
#undef T_n
#undef t2z
#undef t2y
#undef t2x
#undef t1z
#undef t1y
#undef t1x
#undef tz
#undef ty
#undef tx
#undef w_geom_t
#undef w_geom_c
#undef w_geom
#undef w_mult
#undef l_range
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'Origin'. */
  SIG_MESSAGE("Origin (McDisplay)");
  printf("MCDISPLAY: component %s\n", "Origin");
#define mccompcurname  Origin
#define mccompcurtype  Progress_bar
#define mccompcurindex 2
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
#line 27573 "./ESS_IN5_reprate.c"
}   /* End of Origin=Progress_bar() SETTING parameter declarations. */
#undef CurrentTime
#undef EndTime
#undef StartTime
#undef IntermediateCnts
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'TOFmoderator_zoom'. */
  SIG_MESSAGE("TOFmoderator_zoom (McDisplay)");
  printf("MCDISPLAY: component %s\n", "TOFmoderator_zoom");
#define mccompcurname  TOFmoderator_zoom
#define mccompcurtype  TOF_monitor
#define mccompcurindex 3
#define nt mccTOFmoderator_zoom_nt
#define TOF_N mccTOFmoderator_zoom_TOF_N
#define TOF_p mccTOFmoderator_zoom_TOF_p
#define TOF_p2 mccTOFmoderator_zoom_TOF_p2
#define t_min mccTOFmoderator_zoom_t_min
#define t_max mccTOFmoderator_zoom_t_max
#define delta_t mccTOFmoderator_zoom_delta_t
{   /* Declarations of TOFmoderator_zoom=TOF_monitor() SETTING parameters. */
char* filename = mccTOFmoderator_zoom_filename;
MCNUM xmin = mccTOFmoderator_zoom_xmin;
MCNUM xmax = mccTOFmoderator_zoom_xmax;
MCNUM ymin = mccTOFmoderator_zoom_ymin;
MCNUM ymax = mccTOFmoderator_zoom_ymax;
MCNUM xwidth = mccTOFmoderator_zoom_xwidth;
MCNUM yheight = mccTOFmoderator_zoom_yheight;
MCNUM tmin = mccTOFmoderator_zoom_tmin;
MCNUM tmax = mccTOFmoderator_zoom_tmax;
MCNUM dt = mccTOFmoderator_zoom_dt;
MCNUM restore_neutron = mccTOFmoderator_zoom_restore_neutron;
int nowritefile = mccTOFmoderator_zoom_nowritefile;
#line 128 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/TOF_monitor.comp"
{
  
  multiline(5, (double)xmin, (double)ymin, 0.0,
               (double)xmax, (double)ymin, 0.0,
               (double)xmax, (double)ymax, 0.0,
               (double)xmin, (double)ymax, 0.0,
               (double)xmin, (double)ymin, 0.0);
}
#line 27618 "./ESS_IN5_reprate.c"
}   /* End of TOFmoderator_zoom=TOF_monitor() SETTING parameter declarations. */
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

  /* MCDISPLAY code for component 'TOFmoderator'. */
  SIG_MESSAGE("TOFmoderator (McDisplay)");
  printf("MCDISPLAY: component %s\n", "TOFmoderator");
#define mccompcurname  TOFmoderator
#define mccompcurtype  TOF_monitor
#define mccompcurindex 4
#define nt mccTOFmoderator_nt
#define TOF_N mccTOFmoderator_TOF_N
#define TOF_p mccTOFmoderator_TOF_p
#define TOF_p2 mccTOFmoderator_TOF_p2
#define t_min mccTOFmoderator_t_min
#define t_max mccTOFmoderator_t_max
#define delta_t mccTOFmoderator_delta_t
{   /* Declarations of TOFmoderator=TOF_monitor() SETTING parameters. */
char* filename = mccTOFmoderator_filename;
MCNUM xmin = mccTOFmoderator_xmin;
MCNUM xmax = mccTOFmoderator_xmax;
MCNUM ymin = mccTOFmoderator_ymin;
MCNUM ymax = mccTOFmoderator_ymax;
MCNUM xwidth = mccTOFmoderator_xwidth;
MCNUM yheight = mccTOFmoderator_yheight;
MCNUM tmin = mccTOFmoderator_tmin;
MCNUM tmax = mccTOFmoderator_tmax;
MCNUM dt = mccTOFmoderator_dt;
MCNUM restore_neutron = mccTOFmoderator_restore_neutron;
int nowritefile = mccTOFmoderator_nowritefile;
#line 128 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/TOF_monitor.comp"
{
  
  multiline(5, (double)xmin, (double)ymin, 0.0,
               (double)xmax, (double)ymin, 0.0,
               (double)xmax, (double)ymax, 0.0,
               (double)xmin, (double)ymax, 0.0,
               (double)xmin, (double)ymin, 0.0);
}
#line 27666 "./ESS_IN5_reprate.c"
}   /* End of TOFmoderator=TOF_monitor() SETTING parameter declarations. */
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

  /* MCDISPLAY code for component 'Lmon_guistart'. */
  SIG_MESSAGE("Lmon_guistart (McDisplay)");
  printf("MCDISPLAY: component %s\n", "Lmon_guistart");
#define mccompcurname  Lmon_guistart
#define mccompcurtype  L_monitor
#define mccompcurindex 5
#define nL mccLmon_guistart_nL
#define L_N mccLmon_guistart_L_N
#define L_p mccLmon_guistart_L_p
#define L_p2 mccLmon_guistart_L_p2
{   /* Declarations of Lmon_guistart=L_monitor() SETTING parameters. */
char* filename = mccLmon_guistart_filename;
MCNUM xmin = mccLmon_guistart_xmin;
MCNUM xmax = mccLmon_guistart_xmax;
MCNUM ymin = mccLmon_guistart_ymin;
MCNUM ymax = mccLmon_guistart_ymax;
MCNUM xwidth = mccLmon_guistart_xwidth;
MCNUM yheight = mccLmon_guistart_yheight;
MCNUM Lmin = mccLmon_guistart_Lmin;
MCNUM Lmax = mccLmon_guistart_Lmax;
MCNUM restore_neutron = mccLmon_guistart_restore_neutron;
int nowritefile = mccLmon_guistart_nowritefile;
#line 120 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/L_monitor.comp"
{
  
  multiline(5, (double)xmin, (double)ymin, 0.0,
               (double)xmax, (double)ymin, 0.0,
               (double)xmax, (double)ymax, 0.0,
               (double)xmin, (double)ymax, 0.0,
               (double)xmin, (double)ymin, 0.0);
}
#line 27710 "./ESS_IN5_reprate.c"
}   /* End of Lmon_guistart=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'Lmon_normalize'. */
  SIG_MESSAGE("Lmon_normalize (McDisplay)");
  printf("MCDISPLAY: component %s\n", "Lmon_normalize");
#define mccompcurname  Lmon_normalize
#define mccompcurtype  L_monitor
#define mccompcurindex 6
#define nL mccLmon_normalize_nL
#define L_N mccLmon_normalize_L_N
#define L_p mccLmon_normalize_L_p
#define L_p2 mccLmon_normalize_L_p2
{   /* Declarations of Lmon_normalize=L_monitor() SETTING parameters. */
char* filename = mccLmon_normalize_filename;
MCNUM xmin = mccLmon_normalize_xmin;
MCNUM xmax = mccLmon_normalize_xmax;
MCNUM ymin = mccLmon_normalize_ymin;
MCNUM ymax = mccLmon_normalize_ymax;
MCNUM xwidth = mccLmon_normalize_xwidth;
MCNUM yheight = mccLmon_normalize_yheight;
MCNUM Lmin = mccLmon_normalize_Lmin;
MCNUM Lmax = mccLmon_normalize_Lmax;
MCNUM restore_neutron = mccLmon_normalize_restore_neutron;
int nowritefile = mccLmon_normalize_nowritefile;
#line 120 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/L_monitor.comp"
{
  
  multiline(5, (double)xmin, (double)ymin, 0.0,
               (double)xmax, (double)ymin, 0.0,
               (double)xmax, (double)ymax, 0.0,
               (double)xmin, (double)ymax, 0.0,
               (double)xmin, (double)ymin, 0.0);
}
#line 27751 "./ESS_IN5_reprate.c"
}   /* End of Lmon_normalize=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'Guide1'. */
  SIG_MESSAGE("Guide1 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "Guide1");
#define mccompcurname  Guide1
#define mccompcurtype  Guide
#define mccompcurindex 7
#define pTable mccGuide1_pTable
{   /* Declarations of Guide1=Guide() SETTING parameters. */
char* reflect = mccGuide1_reflect;
MCNUM w1 = mccGuide1_w1;
MCNUM h1 = mccGuide1_h1;
MCNUM w2 = mccGuide1_w2;
MCNUM h2 = mccGuide1_h2;
MCNUM l = mccGuide1_l;
MCNUM R0 = mccGuide1_R0;
MCNUM Qc = mccGuide1_Qc;
MCNUM alpha = mccGuide1_alpha;
MCNUM m = mccGuide1_m;
MCNUM W = mccGuide1_W;
#line 201 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide.comp"
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
#line 27800 "./ESS_IN5_reprate.c"
}   /* End of Guide1=Guide() SETTING parameter declarations. */
#undef pTable
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'Lmonslow1'. */
  SIG_MESSAGE("Lmonslow1 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "Lmonslow1");
#define mccompcurname  Lmonslow1
#define mccompcurtype  L_monitor
#define mccompcurindex 8
#define nL mccLmonslow1_nL
#define L_N mccLmonslow1_L_N
#define L_p mccLmonslow1_L_p
#define L_p2 mccLmonslow1_L_p2
{   /* Declarations of Lmonslow1=L_monitor() SETTING parameters. */
char* filename = mccLmonslow1_filename;
MCNUM xmin = mccLmonslow1_xmin;
MCNUM xmax = mccLmonslow1_xmax;
MCNUM ymin = mccLmonslow1_ymin;
MCNUM ymax = mccLmonslow1_ymax;
MCNUM xwidth = mccLmonslow1_xwidth;
MCNUM yheight = mccLmonslow1_yheight;
MCNUM Lmin = mccLmonslow1_Lmin;
MCNUM Lmax = mccLmonslow1_Lmax;
MCNUM restore_neutron = mccLmonslow1_restore_neutron;
int nowritefile = mccLmonslow1_nowritefile;
#line 120 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/L_monitor.comp"
{
  
  multiline(5, (double)xmin, (double)ymin, 0.0,
               (double)xmax, (double)ymin, 0.0,
               (double)xmax, (double)ymax, 0.0,
               (double)xmin, (double)ymax, 0.0,
               (double)xmin, (double)ymin, 0.0);
}
#line 27838 "./ESS_IN5_reprate.c"
}   /* End of Lmonslow1=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'PSDslow1'. */
  SIG_MESSAGE("PSDslow1 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "PSDslow1");
#define mccompcurname  PSDslow1
#define mccompcurtype  PSD_monitor
#define mccompcurindex 9
#define PSD_N mccPSDslow1_PSD_N
#define PSD_p mccPSDslow1_PSD_p
#define PSD_p2 mccPSDslow1_PSD_p2
{   /* Declarations of PSDslow1=PSD_monitor() SETTING parameters. */
int nx = mccPSDslow1_nx;
int ny = mccPSDslow1_ny;
char* filename = mccPSDslow1_filename;
MCNUM xmin = mccPSDslow1_xmin;
MCNUM xmax = mccPSDslow1_xmax;
MCNUM ymin = mccPSDslow1_ymin;
MCNUM ymax = mccPSDslow1_ymax;
MCNUM xwidth = mccPSDslow1_xwidth;
MCNUM yheight = mccPSDslow1_yheight;
MCNUM restore_neutron = mccPSDslow1_restore_neutron;
#line 129 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
{
  multiline(5,
    (double)xmin, (double)ymin, 0.0,
    (double)xmax, (double)ymin, 0.0,
    (double)xmax, (double)ymax, 0.0,
    (double)xmin, (double)ymax, 0.0,
    (double)xmin, (double)ymin, 0.0);
}
#line 27877 "./ESS_IN5_reprate.c"
}   /* End of PSDslow1=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'FOchop1'. */
  SIG_MESSAGE("FOchop1 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "FOchop1");
#define mccompcurname  FOchop1
#define mccompcurtype  DiskChopper
#define mccompcurindex 10
#define Tg mccFOchop1_Tg
#define To mccFOchop1_To
#define delta_y mccFOchop1_delta_y
#define height mccFOchop1_height
#define omega mccFOchop1_omega
{   /* Declarations of FOchop1=DiskChopper() SETTING parameters. */
MCNUM theta_0 = mccFOchop1_theta_0;
MCNUM radius = mccFOchop1_radius;
MCNUM yheight = mccFOchop1_yheight;
MCNUM nu = mccFOchop1_nu;
MCNUM nslit = mccFOchop1_nslit;
MCNUM jitter = mccFOchop1_jitter;
MCNUM delay = mccFOchop1_delay;
MCNUM isfirst = mccFOchop1_isfirst;
MCNUM n_pulse = mccFOchop1_n_pulse;
MCNUM abs_out = mccFOchop1_abs_out;
MCNUM phase = mccFOchop1_phase;
MCNUM xwidth = mccFOchop1_xwidth;
MCNUM verbose = mccFOchop1_verbose;
#line 168 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/DiskChopper.comp"
{

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
}
#line 27938 "./ESS_IN5_reprate.c"
}   /* End of FOchop1=DiskChopper() SETTING parameter declarations. */
#undef omega
#undef height
#undef delta_y
#undef To
#undef Tg
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'TOFLmon1'. */
  SIG_MESSAGE("TOFLmon1 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "TOFLmon1");
#define mccompcurname  TOFLmon1
#define mccompcurtype  TOFLambda_monitor
#define mccompcurindex 11
#define nL mccTOFLmon1_nL
#define nt mccTOFLmon1_nt
#define tmin mccTOFLmon1_tmin
#define tmax mccTOFLmon1_tmax
#define tt_0 mccTOFLmon1_tt_0
#define tt_1 mccTOFLmon1_tt_1
#define TOFL_N mccTOFLmon1_TOFL_N
#define TOFL_p mccTOFLmon1_TOFL_p
#define TOFL_p2 mccTOFLmon1_TOFL_p2
{   /* Declarations of TOFLmon1=TOFLambda_monitor() SETTING parameters. */
char* filename = mccTOFLmon1_filename;
MCNUM xmin = mccTOFLmon1_xmin;
MCNUM xmax = mccTOFLmon1_xmax;
MCNUM ymin = mccTOFLmon1_ymin;
MCNUM ymax = mccTOFLmon1_ymax;
MCNUM xwidth = mccTOFLmon1_xwidth;
MCNUM yheight = mccTOFLmon1_yheight;
MCNUM Lmin = mccTOFLmon1_Lmin;
MCNUM Lmax = mccTOFLmon1_Lmax;
MCNUM restore_neutron = mccTOFLmon1_restore_neutron;
int nowritefile = mccTOFLmon1_nowritefile;
#line 123 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/TOFLambda_monitor.comp"
{
    
    multiline(5, (double)xmin, (double)ymin, 0.0,
                 (double)xmax, (double)ymin, 0.0,
                 (double)xmax, (double)ymax, 0.0,
                 (double)xmin, (double)ymax, 0.0,
                 (double)xmin, (double)ymin, 0.0);
}
#line 27985 "./ESS_IN5_reprate.c"
}   /* End of TOFLmon1=TOFLambda_monitor() SETTING parameter declarations. */
#undef TOFL_p2
#undef TOFL_p
#undef TOFL_N
#undef tt_1
#undef tt_0
#undef tmax
#undef tmin
#undef nt
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'Lmon_afterslow1'. */
  SIG_MESSAGE("Lmon_afterslow1 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "Lmon_afterslow1");
#define mccompcurname  Lmon_afterslow1
#define mccompcurtype  L_monitor
#define mccompcurindex 12
#define nL mccLmon_afterslow1_nL
#define L_N mccLmon_afterslow1_L_N
#define L_p mccLmon_afterslow1_L_p
#define L_p2 mccLmon_afterslow1_L_p2
{   /* Declarations of Lmon_afterslow1=L_monitor() SETTING parameters. */
char* filename = mccLmon_afterslow1_filename;
MCNUM xmin = mccLmon_afterslow1_xmin;
MCNUM xmax = mccLmon_afterslow1_xmax;
MCNUM ymin = mccLmon_afterslow1_ymin;
MCNUM ymax = mccLmon_afterslow1_ymax;
MCNUM xwidth = mccLmon_afterslow1_xwidth;
MCNUM yheight = mccLmon_afterslow1_yheight;
MCNUM Lmin = mccLmon_afterslow1_Lmin;
MCNUM Lmax = mccLmon_afterslow1_Lmax;
MCNUM restore_neutron = mccLmon_afterslow1_restore_neutron;
int nowritefile = mccLmon_afterslow1_nowritefile;
#line 120 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/L_monitor.comp"
{
  
  multiline(5, (double)xmin, (double)ymin, 0.0,
               (double)xmax, (double)ymin, 0.0,
               (double)xmax, (double)ymax, 0.0,
               (double)xmin, (double)ymax, 0.0,
               (double)xmin, (double)ymin, 0.0);
}
#line 28031 "./ESS_IN5_reprate.c"
}   /* End of Lmon_afterslow1=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'PSD_afterslow1'. */
  SIG_MESSAGE("PSD_afterslow1 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "PSD_afterslow1");
#define mccompcurname  PSD_afterslow1
#define mccompcurtype  PSD_monitor
#define mccompcurindex 13
#define PSD_N mccPSD_afterslow1_PSD_N
#define PSD_p mccPSD_afterslow1_PSD_p
#define PSD_p2 mccPSD_afterslow1_PSD_p2
{   /* Declarations of PSD_afterslow1=PSD_monitor() SETTING parameters. */
int nx = mccPSD_afterslow1_nx;
int ny = mccPSD_afterslow1_ny;
char* filename = mccPSD_afterslow1_filename;
MCNUM xmin = mccPSD_afterslow1_xmin;
MCNUM xmax = mccPSD_afterslow1_xmax;
MCNUM ymin = mccPSD_afterslow1_ymin;
MCNUM ymax = mccPSD_afterslow1_ymax;
MCNUM xwidth = mccPSD_afterslow1_xwidth;
MCNUM yheight = mccPSD_afterslow1_yheight;
MCNUM restore_neutron = mccPSD_afterslow1_restore_neutron;
#line 129 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
{
  multiline(5,
    (double)xmin, (double)ymin, 0.0,
    (double)xmax, (double)ymin, 0.0,
    (double)xmax, (double)ymax, 0.0,
    (double)xmin, (double)ymax, 0.0,
    (double)xmin, (double)ymin, 0.0);
}
#line 28070 "./ESS_IN5_reprate.c"
}   /* End of PSD_afterslow1=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'Guidelong1'. */
  SIG_MESSAGE("Guidelong1 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "Guidelong1");
#define mccompcurname  Guidelong1
#define mccompcurtype  Guide
#define mccompcurindex 14
#define pTable mccGuidelong1_pTable
{   /* Declarations of Guidelong1=Guide() SETTING parameters. */
char* reflect = mccGuidelong1_reflect;
MCNUM w1 = mccGuidelong1_w1;
MCNUM h1 = mccGuidelong1_h1;
MCNUM w2 = mccGuidelong1_w2;
MCNUM h2 = mccGuidelong1_h2;
MCNUM l = mccGuidelong1_l;
MCNUM R0 = mccGuidelong1_R0;
MCNUM Qc = mccGuidelong1_Qc;
MCNUM alpha = mccGuidelong1_alpha;
MCNUM m = mccGuidelong1_m;
MCNUM W = mccGuidelong1_W;
#line 201 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide.comp"
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
#line 28118 "./ESS_IN5_reprate.c"
}   /* End of Guidelong1=Guide() SETTING parameter declarations. */
#undef pTable
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'Guidelong1b'. */
  SIG_MESSAGE("Guidelong1b (McDisplay)");
  printf("MCDISPLAY: component %s\n", "Guidelong1b");
#define mccompcurname  Guidelong1b
#define mccompcurtype  Guide
#define mccompcurindex 15
#define pTable mccGuidelong1b_pTable
{   /* Declarations of Guidelong1b=Guide() SETTING parameters. */
char* reflect = mccGuidelong1b_reflect;
MCNUM w1 = mccGuidelong1b_w1;
MCNUM h1 = mccGuidelong1b_h1;
MCNUM w2 = mccGuidelong1b_w2;
MCNUM h2 = mccGuidelong1b_h2;
MCNUM l = mccGuidelong1b_l;
MCNUM R0 = mccGuidelong1b_R0;
MCNUM Qc = mccGuidelong1b_Qc;
MCNUM alpha = mccGuidelong1b_alpha;
MCNUM m = mccGuidelong1b_m;
MCNUM W = mccGuidelong1b_W;
#line 201 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide.comp"
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
#line 28164 "./ESS_IN5_reprate.c"
}   /* End of Guidelong1b=Guide() SETTING parameter declarations. */
#undef pTable
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'Lmon_slow2'. */
  SIG_MESSAGE("Lmon_slow2 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "Lmon_slow2");
#define mccompcurname  Lmon_slow2
#define mccompcurtype  L_monitor
#define mccompcurindex 16
#define nL mccLmon_slow2_nL
#define L_N mccLmon_slow2_L_N
#define L_p mccLmon_slow2_L_p
#define L_p2 mccLmon_slow2_L_p2
{   /* Declarations of Lmon_slow2=L_monitor() SETTING parameters. */
char* filename = mccLmon_slow2_filename;
MCNUM xmin = mccLmon_slow2_xmin;
MCNUM xmax = mccLmon_slow2_xmax;
MCNUM ymin = mccLmon_slow2_ymin;
MCNUM ymax = mccLmon_slow2_ymax;
MCNUM xwidth = mccLmon_slow2_xwidth;
MCNUM yheight = mccLmon_slow2_yheight;
MCNUM Lmin = mccLmon_slow2_Lmin;
MCNUM Lmax = mccLmon_slow2_Lmax;
MCNUM restore_neutron = mccLmon_slow2_restore_neutron;
int nowritefile = mccLmon_slow2_nowritefile;
#line 120 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/L_monitor.comp"
{
  
  multiline(5, (double)xmin, (double)ymin, 0.0,
               (double)xmax, (double)ymin, 0.0,
               (double)xmax, (double)ymax, 0.0,
               (double)xmin, (double)ymax, 0.0,
               (double)xmin, (double)ymin, 0.0);
}
#line 28202 "./ESS_IN5_reprate.c"
}   /* End of Lmon_slow2=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'FOchop2'. */
  SIG_MESSAGE("FOchop2 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "FOchop2");
#define mccompcurname  FOchop2
#define mccompcurtype  DiskChopper
#define mccompcurindex 17
#define Tg mccFOchop2_Tg
#define To mccFOchop2_To
#define delta_y mccFOchop2_delta_y
#define height mccFOchop2_height
#define omega mccFOchop2_omega
{   /* Declarations of FOchop2=DiskChopper() SETTING parameters. */
MCNUM theta_0 = mccFOchop2_theta_0;
MCNUM radius = mccFOchop2_radius;
MCNUM yheight = mccFOchop2_yheight;
MCNUM nu = mccFOchop2_nu;
MCNUM nslit = mccFOchop2_nslit;
MCNUM jitter = mccFOchop2_jitter;
MCNUM delay = mccFOchop2_delay;
MCNUM isfirst = mccFOchop2_isfirst;
MCNUM n_pulse = mccFOchop2_n_pulse;
MCNUM abs_out = mccFOchop2_abs_out;
MCNUM phase = mccFOchop2_phase;
MCNUM xwidth = mccFOchop2_xwidth;
MCNUM verbose = mccFOchop2_verbose;
#line 168 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/DiskChopper.comp"
{

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
}
#line 28264 "./ESS_IN5_reprate.c"
}   /* End of FOchop2=DiskChopper() SETTING parameter declarations. */
#undef omega
#undef height
#undef delta_y
#undef To
#undef Tg
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'Fastchop1'. */
  SIG_MESSAGE("Fastchop1 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "Fastchop1");
#define mccompcurname  Fastchop1
#define mccompcurtype  DiskChopper
#define mccompcurindex 18
#define Tg mccFastchop1_Tg
#define To mccFastchop1_To
#define delta_y mccFastchop1_delta_y
#define height mccFastchop1_height
#define omega mccFastchop1_omega
{   /* Declarations of Fastchop1=DiskChopper() SETTING parameters. */
MCNUM theta_0 = mccFastchop1_theta_0;
MCNUM radius = mccFastchop1_radius;
MCNUM yheight = mccFastchop1_yheight;
MCNUM nu = mccFastchop1_nu;
MCNUM nslit = mccFastchop1_nslit;
MCNUM jitter = mccFastchop1_jitter;
MCNUM delay = mccFastchop1_delay;
MCNUM isfirst = mccFastchop1_isfirst;
MCNUM n_pulse = mccFastchop1_n_pulse;
MCNUM abs_out = mccFastchop1_abs_out;
MCNUM phase = mccFastchop1_phase;
MCNUM xwidth = mccFastchop1_xwidth;
MCNUM verbose = mccFastchop1_verbose;
#line 168 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/DiskChopper.comp"
{

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
}
#line 28327 "./ESS_IN5_reprate.c"
}   /* End of Fastchop1=DiskChopper() SETTING parameter declarations. */
#undef omega
#undef height
#undef delta_y
#undef To
#undef Tg
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'PSD_afterslow2'. */
  SIG_MESSAGE("PSD_afterslow2 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "PSD_afterslow2");
#define mccompcurname  PSD_afterslow2
#define mccompcurtype  PSD_monitor
#define mccompcurindex 19
#define PSD_N mccPSD_afterslow2_PSD_N
#define PSD_p mccPSD_afterslow2_PSD_p
#define PSD_p2 mccPSD_afterslow2_PSD_p2
{   /* Declarations of PSD_afterslow2=PSD_monitor() SETTING parameters. */
int nx = mccPSD_afterslow2_nx;
int ny = mccPSD_afterslow2_ny;
char* filename = mccPSD_afterslow2_filename;
MCNUM xmin = mccPSD_afterslow2_xmin;
MCNUM xmax = mccPSD_afterslow2_xmax;
MCNUM ymin = mccPSD_afterslow2_ymin;
MCNUM ymax = mccPSD_afterslow2_ymax;
MCNUM xwidth = mccPSD_afterslow2_xwidth;
MCNUM yheight = mccPSD_afterslow2_yheight;
MCNUM restore_neutron = mccPSD_afterslow2_restore_neutron;
#line 129 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
{
  multiline(5,
    (double)xmin, (double)ymin, 0.0,
    (double)xmax, (double)ymin, 0.0,
    (double)xmax, (double)ymax, 0.0,
    (double)xmin, (double)ymax, 0.0,
    (double)xmin, (double)ymin, 0.0);
}
#line 28367 "./ESS_IN5_reprate.c"
}   /* End of PSD_afterslow2=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'Lmon_afterslow2'. */
  SIG_MESSAGE("Lmon_afterslow2 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "Lmon_afterslow2");
#define mccompcurname  Lmon_afterslow2
#define mccompcurtype  L_monitor
#define mccompcurindex 20
#define nL mccLmon_afterslow2_nL
#define L_N mccLmon_afterslow2_L_N
#define L_p mccLmon_afterslow2_L_p
#define L_p2 mccLmon_afterslow2_L_p2
{   /* Declarations of Lmon_afterslow2=L_monitor() SETTING parameters. */
char* filename = mccLmon_afterslow2_filename;
MCNUM xmin = mccLmon_afterslow2_xmin;
MCNUM xmax = mccLmon_afterslow2_xmax;
MCNUM ymin = mccLmon_afterslow2_ymin;
MCNUM ymax = mccLmon_afterslow2_ymax;
MCNUM xwidth = mccLmon_afterslow2_xwidth;
MCNUM yheight = mccLmon_afterslow2_yheight;
MCNUM Lmin = mccLmon_afterslow2_Lmin;
MCNUM Lmax = mccLmon_afterslow2_Lmax;
MCNUM restore_neutron = mccLmon_afterslow2_restore_neutron;
int nowritefile = mccLmon_afterslow2_nowritefile;
#line 120 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/L_monitor.comp"
{
  
  multiline(5, (double)xmin, (double)ymin, 0.0,
               (double)xmax, (double)ymin, 0.0,
               (double)xmax, (double)ymax, 0.0,
               (double)xmin, (double)ymax, 0.0,
               (double)xmin, (double)ymin, 0.0);
}
#line 28407 "./ESS_IN5_reprate.c"
}   /* End of Lmon_afterslow2=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'TOFL_afterslow2'. */
  SIG_MESSAGE("TOFL_afterslow2 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "TOFL_afterslow2");
#define mccompcurname  TOFL_afterslow2
#define mccompcurtype  TOFLambda_monitor
#define mccompcurindex 21
#define nL mccTOFL_afterslow2_nL
#define nt mccTOFL_afterslow2_nt
#define tmin mccTOFL_afterslow2_tmin
#define tmax mccTOFL_afterslow2_tmax
#define tt_0 mccTOFL_afterslow2_tt_0
#define tt_1 mccTOFL_afterslow2_tt_1
#define TOFL_N mccTOFL_afterslow2_TOFL_N
#define TOFL_p mccTOFL_afterslow2_TOFL_p
#define TOFL_p2 mccTOFL_afterslow2_TOFL_p2
{   /* Declarations of TOFL_afterslow2=TOFLambda_monitor() SETTING parameters. */
char* filename = mccTOFL_afterslow2_filename;
MCNUM xmin = mccTOFL_afterslow2_xmin;
MCNUM xmax = mccTOFL_afterslow2_xmax;
MCNUM ymin = mccTOFL_afterslow2_ymin;
MCNUM ymax = mccTOFL_afterslow2_ymax;
MCNUM xwidth = mccTOFL_afterslow2_xwidth;
MCNUM yheight = mccTOFL_afterslow2_yheight;
MCNUM Lmin = mccTOFL_afterslow2_Lmin;
MCNUM Lmax = mccTOFL_afterslow2_Lmax;
MCNUM restore_neutron = mccTOFL_afterslow2_restore_neutron;
int nowritefile = mccTOFL_afterslow2_nowritefile;
#line 123 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/TOFLambda_monitor.comp"
{
    
    multiline(5, (double)xmin, (double)ymin, 0.0,
                 (double)xmax, (double)ymin, 0.0,
                 (double)xmax, (double)ymax, 0.0,
                 (double)xmin, (double)ymax, 0.0,
                 (double)xmin, (double)ymin, 0.0);
}
#line 28453 "./ESS_IN5_reprate.c"
}   /* End of TOFL_afterslow2=TOFLambda_monitor() SETTING parameter declarations. */
#undef TOFL_p2
#undef TOFL_p
#undef TOFL_N
#undef tt_1
#undef tt_0
#undef tmax
#undef tmin
#undef nt
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'Guidelong2'. */
  SIG_MESSAGE("Guidelong2 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "Guidelong2");
#define mccompcurname  Guidelong2
#define mccompcurtype  Guide
#define mccompcurindex 22
#define pTable mccGuidelong2_pTable
{   /* Declarations of Guidelong2=Guide() SETTING parameters. */
char* reflect = mccGuidelong2_reflect;
MCNUM w1 = mccGuidelong2_w1;
MCNUM h1 = mccGuidelong2_h1;
MCNUM w2 = mccGuidelong2_w2;
MCNUM h2 = mccGuidelong2_h2;
MCNUM l = mccGuidelong2_l;
MCNUM R0 = mccGuidelong2_R0;
MCNUM Qc = mccGuidelong2_Qc;
MCNUM alpha = mccGuidelong2_alpha;
MCNUM m = mccGuidelong2_m;
MCNUM W = mccGuidelong2_W;
#line 201 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide.comp"
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
#line 28507 "./ESS_IN5_reprate.c"
}   /* End of Guidelong2=Guide() SETTING parameter declarations. */
#undef pTable
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'Lmon_beforeballistic'. */
  SIG_MESSAGE("Lmon_beforeballistic (McDisplay)");
  printf("MCDISPLAY: component %s\n", "Lmon_beforeballistic");
#define mccompcurname  Lmon_beforeballistic
#define mccompcurtype  L_monitor
#define mccompcurindex 23
#define nL mccLmon_beforeballistic_nL
#define L_N mccLmon_beforeballistic_L_N
#define L_p mccLmon_beforeballistic_L_p
#define L_p2 mccLmon_beforeballistic_L_p2
{   /* Declarations of Lmon_beforeballistic=L_monitor() SETTING parameters. */
char* filename = mccLmon_beforeballistic_filename;
MCNUM xmin = mccLmon_beforeballistic_xmin;
MCNUM xmax = mccLmon_beforeballistic_xmax;
MCNUM ymin = mccLmon_beforeballistic_ymin;
MCNUM ymax = mccLmon_beforeballistic_ymax;
MCNUM xwidth = mccLmon_beforeballistic_xwidth;
MCNUM yheight = mccLmon_beforeballistic_yheight;
MCNUM Lmin = mccLmon_beforeballistic_Lmin;
MCNUM Lmax = mccLmon_beforeballistic_Lmax;
MCNUM restore_neutron = mccLmon_beforeballistic_restore_neutron;
int nowritefile = mccLmon_beforeballistic_nowritefile;
#line 120 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/L_monitor.comp"
{
  
  multiline(5, (double)xmin, (double)ymin, 0.0,
               (double)xmax, (double)ymin, 0.0,
               (double)xmax, (double)ymax, 0.0,
               (double)xmin, (double)ymax, 0.0,
               (double)xmin, (double)ymin, 0.0);
}
#line 28545 "./ESS_IN5_reprate.c"
}   /* End of Lmon_beforeballistic=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'PSD_beforeballistic'. */
  SIG_MESSAGE("PSD_beforeballistic (McDisplay)");
  printf("MCDISPLAY: component %s\n", "PSD_beforeballistic");
#define mccompcurname  PSD_beforeballistic
#define mccompcurtype  PSD_monitor
#define mccompcurindex 24
#define PSD_N mccPSD_beforeballistic_PSD_N
#define PSD_p mccPSD_beforeballistic_PSD_p
#define PSD_p2 mccPSD_beforeballistic_PSD_p2
{   /* Declarations of PSD_beforeballistic=PSD_monitor() SETTING parameters. */
int nx = mccPSD_beforeballistic_nx;
int ny = mccPSD_beforeballistic_ny;
char* filename = mccPSD_beforeballistic_filename;
MCNUM xmin = mccPSD_beforeballistic_xmin;
MCNUM xmax = mccPSD_beforeballistic_xmax;
MCNUM ymin = mccPSD_beforeballistic_ymin;
MCNUM ymax = mccPSD_beforeballistic_ymax;
MCNUM xwidth = mccPSD_beforeballistic_xwidth;
MCNUM yheight = mccPSD_beforeballistic_yheight;
MCNUM restore_neutron = mccPSD_beforeballistic_restore_neutron;
#line 129 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
{
  multiline(5,
    (double)xmin, (double)ymin, 0.0,
    (double)xmax, (double)ymin, 0.0,
    (double)xmax, (double)ymax, 0.0,
    (double)xmin, (double)ymax, 0.0,
    (double)xmin, (double)ymin, 0.0);
}
#line 28584 "./ESS_IN5_reprate.c"
}   /* End of PSD_beforeballistic=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'Guidelong2a'. */
  SIG_MESSAGE("Guidelong2a (McDisplay)");
  printf("MCDISPLAY: component %s\n", "Guidelong2a");
#define mccompcurname  Guidelong2a
#define mccompcurtype  Guide
#define mccompcurindex 25
#define pTable mccGuidelong2a_pTable
{   /* Declarations of Guidelong2a=Guide() SETTING parameters. */
char* reflect = mccGuidelong2a_reflect;
MCNUM w1 = mccGuidelong2a_w1;
MCNUM h1 = mccGuidelong2a_h1;
MCNUM w2 = mccGuidelong2a_w2;
MCNUM h2 = mccGuidelong2a_h2;
MCNUM l = mccGuidelong2a_l;
MCNUM R0 = mccGuidelong2a_R0;
MCNUM Qc = mccGuidelong2a_Qc;
MCNUM alpha = mccGuidelong2a_alpha;
MCNUM m = mccGuidelong2a_m;
MCNUM W = mccGuidelong2a_W;
#line 201 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide.comp"
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
#line 28632 "./ESS_IN5_reprate.c"
}   /* End of Guidelong2a=Guide() SETTING parameter declarations. */
#undef pTable
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'Lmonfast2'. */
  SIG_MESSAGE("Lmonfast2 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "Lmonfast2");
#define mccompcurname  Lmonfast2
#define mccompcurtype  L_monitor
#define mccompcurindex 26
#define nL mccLmonfast2_nL
#define L_N mccLmonfast2_L_N
#define L_p mccLmonfast2_L_p
#define L_p2 mccLmonfast2_L_p2
{   /* Declarations of Lmonfast2=L_monitor() SETTING parameters. */
char* filename = mccLmonfast2_filename;
MCNUM xmin = mccLmonfast2_xmin;
MCNUM xmax = mccLmonfast2_xmax;
MCNUM ymin = mccLmonfast2_ymin;
MCNUM ymax = mccLmonfast2_ymax;
MCNUM xwidth = mccLmonfast2_xwidth;
MCNUM yheight = mccLmonfast2_yheight;
MCNUM Lmin = mccLmonfast2_Lmin;
MCNUM Lmax = mccLmonfast2_Lmax;
MCNUM restore_neutron = mccLmonfast2_restore_neutron;
int nowritefile = mccLmonfast2_nowritefile;
#line 120 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/L_monitor.comp"
{
  
  multiline(5, (double)xmin, (double)ymin, 0.0,
               (double)xmax, (double)ymin, 0.0,
               (double)xmax, (double)ymax, 0.0,
               (double)xmin, (double)ymax, 0.0,
               (double)xmin, (double)ymin, 0.0);
}
#line 28670 "./ESS_IN5_reprate.c"
}   /* End of Lmonfast2=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'Lmonfast2_zoom'. */
  SIG_MESSAGE("Lmonfast2_zoom (McDisplay)");
  printf("MCDISPLAY: component %s\n", "Lmonfast2_zoom");
#define mccompcurname  Lmonfast2_zoom
#define mccompcurtype  L_monitor
#define mccompcurindex 27
#define nL mccLmonfast2_zoom_nL
#define L_N mccLmonfast2_zoom_L_N
#define L_p mccLmonfast2_zoom_L_p
#define L_p2 mccLmonfast2_zoom_L_p2
{   /* Declarations of Lmonfast2_zoom=L_monitor() SETTING parameters. */
char* filename = mccLmonfast2_zoom_filename;
MCNUM xmin = mccLmonfast2_zoom_xmin;
MCNUM xmax = mccLmonfast2_zoom_xmax;
MCNUM ymin = mccLmonfast2_zoom_ymin;
MCNUM ymax = mccLmonfast2_zoom_ymax;
MCNUM xwidth = mccLmonfast2_zoom_xwidth;
MCNUM yheight = mccLmonfast2_zoom_yheight;
MCNUM Lmin = mccLmonfast2_zoom_Lmin;
MCNUM Lmax = mccLmonfast2_zoom_Lmax;
MCNUM restore_neutron = mccLmonfast2_zoom_restore_neutron;
int nowritefile = mccLmonfast2_zoom_nowritefile;
#line 120 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/L_monitor.comp"
{
  
  multiline(5, (double)xmin, (double)ymin, 0.0,
               (double)xmax, (double)ymin, 0.0,
               (double)xmax, (double)ymax, 0.0,
               (double)xmin, (double)ymax, 0.0,
               (double)xmin, (double)ymin, 0.0);
}
#line 28711 "./ESS_IN5_reprate.c"
}   /* End of Lmonfast2_zoom=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'TOFLfast2'. */
  SIG_MESSAGE("TOFLfast2 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "TOFLfast2");
#define mccompcurname  TOFLfast2
#define mccompcurtype  TOFLambda_monitor
#define mccompcurindex 28
#define nL mccTOFLfast2_nL
#define nt mccTOFLfast2_nt
#define tmin mccTOFLfast2_tmin
#define tmax mccTOFLfast2_tmax
#define tt_0 mccTOFLfast2_tt_0
#define tt_1 mccTOFLfast2_tt_1
#define TOFL_N mccTOFLfast2_TOFL_N
#define TOFL_p mccTOFLfast2_TOFL_p
#define TOFL_p2 mccTOFLfast2_TOFL_p2
{   /* Declarations of TOFLfast2=TOFLambda_monitor() SETTING parameters. */
char* filename = mccTOFLfast2_filename;
MCNUM xmin = mccTOFLfast2_xmin;
MCNUM xmax = mccTOFLfast2_xmax;
MCNUM ymin = mccTOFLfast2_ymin;
MCNUM ymax = mccTOFLfast2_ymax;
MCNUM xwidth = mccTOFLfast2_xwidth;
MCNUM yheight = mccTOFLfast2_yheight;
MCNUM Lmin = mccTOFLfast2_Lmin;
MCNUM Lmax = mccTOFLfast2_Lmax;
MCNUM restore_neutron = mccTOFLfast2_restore_neutron;
int nowritefile = mccTOFLfast2_nowritefile;
#line 123 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/TOFLambda_monitor.comp"
{
    
    multiline(5, (double)xmin, (double)ymin, 0.0,
                 (double)xmax, (double)ymin, 0.0,
                 (double)xmax, (double)ymax, 0.0,
                 (double)xmin, (double)ymax, 0.0,
                 (double)xmin, (double)ymin, 0.0);
}
#line 28757 "./ESS_IN5_reprate.c"
}   /* End of TOFLfast2=TOFLambda_monitor() SETTING parameter declarations. */
#undef TOFL_p2
#undef TOFL_p
#undef TOFL_N
#undef tt_1
#undef tt_0
#undef tmax
#undef tmin
#undef nt
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'TOFLfast2zoom'. */
  SIG_MESSAGE("TOFLfast2zoom (McDisplay)");
  printf("MCDISPLAY: component %s\n", "TOFLfast2zoom");
#define mccompcurname  TOFLfast2zoom
#define mccompcurtype  TOFLambda_monitor
#define mccompcurindex 29
#define nL mccTOFLfast2zoom_nL
#define nt mccTOFLfast2zoom_nt
#define tmin mccTOFLfast2zoom_tmin
#define tmax mccTOFLfast2zoom_tmax
#define tt_0 mccTOFLfast2zoom_tt_0
#define tt_1 mccTOFLfast2zoom_tt_1
#define TOFL_N mccTOFLfast2zoom_TOFL_N
#define TOFL_p mccTOFLfast2zoom_TOFL_p
#define TOFL_p2 mccTOFLfast2zoom_TOFL_p2
{   /* Declarations of TOFLfast2zoom=TOFLambda_monitor() SETTING parameters. */
char* filename = mccTOFLfast2zoom_filename;
MCNUM xmin = mccTOFLfast2zoom_xmin;
MCNUM xmax = mccTOFLfast2zoom_xmax;
MCNUM ymin = mccTOFLfast2zoom_ymin;
MCNUM ymax = mccTOFLfast2zoom_ymax;
MCNUM xwidth = mccTOFLfast2zoom_xwidth;
MCNUM yheight = mccTOFLfast2zoom_yheight;
MCNUM Lmin = mccTOFLfast2zoom_Lmin;
MCNUM Lmax = mccTOFLfast2zoom_Lmax;
MCNUM restore_neutron = mccTOFLfast2zoom_restore_neutron;
int nowritefile = mccTOFLfast2zoom_nowritefile;
#line 123 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/TOFLambda_monitor.comp"
{
    
    multiline(5, (double)xmin, (double)ymin, 0.0,
                 (double)xmax, (double)ymin, 0.0,
                 (double)xmax, (double)ymax, 0.0,
                 (double)xmin, (double)ymax, 0.0,
                 (double)xmin, (double)ymin, 0.0);
}
#line 28808 "./ESS_IN5_reprate.c"
}   /* End of TOFLfast2zoom=TOFLambda_monitor() SETTING parameter declarations. */
#undef TOFL_p2
#undef TOFL_p
#undef TOFL_N
#undef tt_1
#undef tt_0
#undef tmax
#undef tmin
#undef nt
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'PSDfast2'. */
  SIG_MESSAGE("PSDfast2 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "PSDfast2");
#define mccompcurname  PSDfast2
#define mccompcurtype  PSD_monitor
#define mccompcurindex 30
#define PSD_N mccPSDfast2_PSD_N
#define PSD_p mccPSDfast2_PSD_p
#define PSD_p2 mccPSDfast2_PSD_p2
{   /* Declarations of PSDfast2=PSD_monitor() SETTING parameters. */
int nx = mccPSDfast2_nx;
int ny = mccPSDfast2_ny;
char* filename = mccPSDfast2_filename;
MCNUM xmin = mccPSDfast2_xmin;
MCNUM xmax = mccPSDfast2_xmax;
MCNUM ymin = mccPSDfast2_ymin;
MCNUM ymax = mccPSDfast2_ymax;
MCNUM xwidth = mccPSDfast2_xwidth;
MCNUM yheight = mccPSDfast2_yheight;
MCNUM restore_neutron = mccPSDfast2_restore_neutron;
#line 129 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
{
  multiline(5,
    (double)xmin, (double)ymin, 0.0,
    (double)xmax, (double)ymin, 0.0,
    (double)xmax, (double)ymax, 0.0,
    (double)xmin, (double)ymax, 0.0,
    (double)xmin, (double)ymin, 0.0);
}
#line 28852 "./ESS_IN5_reprate.c"
}   /* End of PSDfast2=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'Fastchop2'. */
  SIG_MESSAGE("Fastchop2 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "Fastchop2");
#define mccompcurname  Fastchop2
#define mccompcurtype  DiskChopper
#define mccompcurindex 31
#define Tg mccFastchop2_Tg
#define To mccFastchop2_To
#define delta_y mccFastchop2_delta_y
#define height mccFastchop2_height
#define omega mccFastchop2_omega
{   /* Declarations of Fastchop2=DiskChopper() SETTING parameters. */
MCNUM theta_0 = mccFastchop2_theta_0;
MCNUM radius = mccFastchop2_radius;
MCNUM yheight = mccFastchop2_yheight;
MCNUM nu = mccFastchop2_nu;
MCNUM nslit = mccFastchop2_nslit;
MCNUM jitter = mccFastchop2_jitter;
MCNUM delay = mccFastchop2_delay;
MCNUM isfirst = mccFastchop2_isfirst;
MCNUM n_pulse = mccFastchop2_n_pulse;
MCNUM abs_out = mccFastchop2_abs_out;
MCNUM phase = mccFastchop2_phase;
MCNUM xwidth = mccFastchop2_xwidth;
MCNUM verbose = mccFastchop2_verbose;
#line 168 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/DiskChopper.comp"
{

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
}
#line 28913 "./ESS_IN5_reprate.c"
}   /* End of Fastchop2=DiskChopper() SETTING parameter declarations. */
#undef omega
#undef height
#undef delta_y
#undef To
#undef Tg
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'Fastchop2counter'. */
  SIG_MESSAGE("Fastchop2counter (McDisplay)");
  printf("MCDISPLAY: component %s\n", "Fastchop2counter");
#define mccompcurname  Fastchop2counter
#define mccompcurtype  DiskChopper
#define mccompcurindex 32
#define Tg mccFastchop2counter_Tg
#define To mccFastchop2counter_To
#define delta_y mccFastchop2counter_delta_y
#define height mccFastchop2counter_height
#define omega mccFastchop2counter_omega
{   /* Declarations of Fastchop2counter=DiskChopper() SETTING parameters. */
MCNUM theta_0 = mccFastchop2counter_theta_0;
MCNUM radius = mccFastchop2counter_radius;
MCNUM yheight = mccFastchop2counter_yheight;
MCNUM nu = mccFastchop2counter_nu;
MCNUM nslit = mccFastchop2counter_nslit;
MCNUM jitter = mccFastchop2counter_jitter;
MCNUM delay = mccFastchop2counter_delay;
MCNUM isfirst = mccFastchop2counter_isfirst;
MCNUM n_pulse = mccFastchop2counter_n_pulse;
MCNUM abs_out = mccFastchop2counter_abs_out;
MCNUM phase = mccFastchop2counter_phase;
MCNUM xwidth = mccFastchop2counter_xwidth;
MCNUM verbose = mccFastchop2counter_verbose;
#line 168 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/DiskChopper.comp"
{

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
}
#line 28976 "./ESS_IN5_reprate.c"
}   /* End of Fastchop2counter=DiskChopper() SETTING parameter declarations. */
#undef omega
#undef height
#undef delta_y
#undef To
#undef Tg
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'FOchop3'. */
  SIG_MESSAGE("FOchop3 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "FOchop3");
#define mccompcurname  FOchop3
#define mccompcurtype  DiskChopper
#define mccompcurindex 33
#define Tg mccFOchop3_Tg
#define To mccFOchop3_To
#define delta_y mccFOchop3_delta_y
#define height mccFOchop3_height
#define omega mccFOchop3_omega
{   /* Declarations of FOchop3=DiskChopper() SETTING parameters. */
MCNUM theta_0 = mccFOchop3_theta_0;
MCNUM radius = mccFOchop3_radius;
MCNUM yheight = mccFOchop3_yheight;
MCNUM nu = mccFOchop3_nu;
MCNUM nslit = mccFOchop3_nslit;
MCNUM jitter = mccFOchop3_jitter;
MCNUM delay = mccFOchop3_delay;
MCNUM isfirst = mccFOchop3_isfirst;
MCNUM n_pulse = mccFOchop3_n_pulse;
MCNUM abs_out = mccFOchop3_abs_out;
MCNUM phase = mccFOchop3_phase;
MCNUM xwidth = mccFOchop3_xwidth;
MCNUM verbose = mccFOchop3_verbose;
#line 168 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/DiskChopper.comp"
{

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
}
#line 29039 "./ESS_IN5_reprate.c"
}   /* End of FOchop3=DiskChopper() SETTING parameter declarations. */
#undef omega
#undef height
#undef delta_y
#undef To
#undef Tg
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'TOFfast2_zoom'. */
  SIG_MESSAGE("TOFfast2_zoom (McDisplay)");
  printf("MCDISPLAY: component %s\n", "TOFfast2_zoom");
#define mccompcurname  TOFfast2_zoom
#define mccompcurtype  TOF_monitor
#define mccompcurindex 34
#define nt mccTOFfast2_zoom_nt
#define TOF_N mccTOFfast2_zoom_TOF_N
#define TOF_p mccTOFfast2_zoom_TOF_p
#define TOF_p2 mccTOFfast2_zoom_TOF_p2
#define t_min mccTOFfast2_zoom_t_min
#define t_max mccTOFfast2_zoom_t_max
#define delta_t mccTOFfast2_zoom_delta_t
{   /* Declarations of TOFfast2_zoom=TOF_monitor() SETTING parameters. */
char* filename = mccTOFfast2_zoom_filename;
MCNUM xmin = mccTOFfast2_zoom_xmin;
MCNUM xmax = mccTOFfast2_zoom_xmax;
MCNUM ymin = mccTOFfast2_zoom_ymin;
MCNUM ymax = mccTOFfast2_zoom_ymax;
MCNUM xwidth = mccTOFfast2_zoom_xwidth;
MCNUM yheight = mccTOFfast2_zoom_yheight;
MCNUM tmin = mccTOFfast2_zoom_tmin;
MCNUM tmax = mccTOFfast2_zoom_tmax;
MCNUM dt = mccTOFfast2_zoom_dt;
MCNUM restore_neutron = mccTOFfast2_zoom_restore_neutron;
int nowritefile = mccTOFfast2_zoom_nowritefile;
#line 128 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/TOF_monitor.comp"
{
  
  multiline(5, (double)xmin, (double)ymin, 0.0,
               (double)xmax, (double)ymin, 0.0,
               (double)xmax, (double)ymax, 0.0,
               (double)xmin, (double)ymax, 0.0,
               (double)xmin, (double)ymin, 0.0);
}
#line 29085 "./ESS_IN5_reprate.c"
}   /* End of TOFfast2_zoom=TOF_monitor() SETTING parameter declarations. */
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

  /* MCDISPLAY code for component 'Lmon_afterfast2'. */
  SIG_MESSAGE("Lmon_afterfast2 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "Lmon_afterfast2");
#define mccompcurname  Lmon_afterfast2
#define mccompcurtype  L_monitor
#define mccompcurindex 35
#define nL mccLmon_afterfast2_nL
#define L_N mccLmon_afterfast2_L_N
#define L_p mccLmon_afterfast2_L_p
#define L_p2 mccLmon_afterfast2_L_p2
{   /* Declarations of Lmon_afterfast2=L_monitor() SETTING parameters. */
char* filename = mccLmon_afterfast2_filename;
MCNUM xmin = mccLmon_afterfast2_xmin;
MCNUM xmax = mccLmon_afterfast2_xmax;
MCNUM ymin = mccLmon_afterfast2_ymin;
MCNUM ymax = mccLmon_afterfast2_ymax;
MCNUM xwidth = mccLmon_afterfast2_xwidth;
MCNUM yheight = mccLmon_afterfast2_yheight;
MCNUM Lmin = mccLmon_afterfast2_Lmin;
MCNUM Lmax = mccLmon_afterfast2_Lmax;
MCNUM restore_neutron = mccLmon_afterfast2_restore_neutron;
int nowritefile = mccLmon_afterfast2_nowritefile;
#line 120 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/L_monitor.comp"
{
  
  multiline(5, (double)xmin, (double)ymin, 0.0,
               (double)xmax, (double)ymin, 0.0,
               (double)xmax, (double)ymax, 0.0,
               (double)xmin, (double)ymax, 0.0,
               (double)xmin, (double)ymin, 0.0);
}
#line 29129 "./ESS_IN5_reprate.c"
}   /* End of Lmon_afterfast2=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'TOFL_afterfast2'. */
  SIG_MESSAGE("TOFL_afterfast2 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "TOFL_afterfast2");
#define mccompcurname  TOFL_afterfast2
#define mccompcurtype  TOFLambda_monitor
#define mccompcurindex 36
#define nL mccTOFL_afterfast2_nL
#define nt mccTOFL_afterfast2_nt
#define tmin mccTOFL_afterfast2_tmin
#define tmax mccTOFL_afterfast2_tmax
#define tt_0 mccTOFL_afterfast2_tt_0
#define tt_1 mccTOFL_afterfast2_tt_1
#define TOFL_N mccTOFL_afterfast2_TOFL_N
#define TOFL_p mccTOFL_afterfast2_TOFL_p
#define TOFL_p2 mccTOFL_afterfast2_TOFL_p2
{   /* Declarations of TOFL_afterfast2=TOFLambda_monitor() SETTING parameters. */
char* filename = mccTOFL_afterfast2_filename;
MCNUM xmin = mccTOFL_afterfast2_xmin;
MCNUM xmax = mccTOFL_afterfast2_xmax;
MCNUM ymin = mccTOFL_afterfast2_ymin;
MCNUM ymax = mccTOFL_afterfast2_ymax;
MCNUM xwidth = mccTOFL_afterfast2_xwidth;
MCNUM yheight = mccTOFL_afterfast2_yheight;
MCNUM Lmin = mccTOFL_afterfast2_Lmin;
MCNUM Lmax = mccTOFL_afterfast2_Lmax;
MCNUM restore_neutron = mccTOFL_afterfast2_restore_neutron;
int nowritefile = mccTOFL_afterfast2_nowritefile;
#line 123 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/TOFLambda_monitor.comp"
{
    
    multiline(5, (double)xmin, (double)ymin, 0.0,
                 (double)xmax, (double)ymin, 0.0,
                 (double)xmax, (double)ymax, 0.0,
                 (double)xmin, (double)ymax, 0.0,
                 (double)xmin, (double)ymin, 0.0);
}
#line 29175 "./ESS_IN5_reprate.c"
}   /* End of TOFL_afterfast2=TOFLambda_monitor() SETTING parameter declarations. */
#undef TOFL_p2
#undef TOFL_p
#undef TOFL_N
#undef tt_1
#undef tt_0
#undef tmax
#undef tmin
#undef nt
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'TOFL_afterfast2_zoom'. */
  SIG_MESSAGE("TOFL_afterfast2_zoom (McDisplay)");
  printf("MCDISPLAY: component %s\n", "TOFL_afterfast2_zoom");
#define mccompcurname  TOFL_afterfast2_zoom
#define mccompcurtype  TOFLambda_monitor
#define mccompcurindex 37
#define nL mccTOFL_afterfast2_zoom_nL
#define nt mccTOFL_afterfast2_zoom_nt
#define tmin mccTOFL_afterfast2_zoom_tmin
#define tmax mccTOFL_afterfast2_zoom_tmax
#define tt_0 mccTOFL_afterfast2_zoom_tt_0
#define tt_1 mccTOFL_afterfast2_zoom_tt_1
#define TOFL_N mccTOFL_afterfast2_zoom_TOFL_N
#define TOFL_p mccTOFL_afterfast2_zoom_TOFL_p
#define TOFL_p2 mccTOFL_afterfast2_zoom_TOFL_p2
{   /* Declarations of TOFL_afterfast2_zoom=TOFLambda_monitor() SETTING parameters. */
char* filename = mccTOFL_afterfast2_zoom_filename;
MCNUM xmin = mccTOFL_afterfast2_zoom_xmin;
MCNUM xmax = mccTOFL_afterfast2_zoom_xmax;
MCNUM ymin = mccTOFL_afterfast2_zoom_ymin;
MCNUM ymax = mccTOFL_afterfast2_zoom_ymax;
MCNUM xwidth = mccTOFL_afterfast2_zoom_xwidth;
MCNUM yheight = mccTOFL_afterfast2_zoom_yheight;
MCNUM Lmin = mccTOFL_afterfast2_zoom_Lmin;
MCNUM Lmax = mccTOFL_afterfast2_zoom_Lmax;
MCNUM restore_neutron = mccTOFL_afterfast2_zoom_restore_neutron;
int nowritefile = mccTOFL_afterfast2_zoom_nowritefile;
#line 123 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/TOFLambda_monitor.comp"
{
    
    multiline(5, (double)xmin, (double)ymin, 0.0,
                 (double)xmax, (double)ymin, 0.0,
                 (double)xmax, (double)ymax, 0.0,
                 (double)xmin, (double)ymax, 0.0,
                 (double)xmin, (double)ymin, 0.0);
}
#line 29226 "./ESS_IN5_reprate.c"
}   /* End of TOFL_afterfast2_zoom=TOFLambda_monitor() SETTING parameter declarations. */
#undef TOFL_p2
#undef TOFL_p
#undef TOFL_N
#undef tt_1
#undef tt_0
#undef tmax
#undef tmin
#undef nt
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'PSD_afterfast2'. */
  SIG_MESSAGE("PSD_afterfast2 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "PSD_afterfast2");
#define mccompcurname  PSD_afterfast2
#define mccompcurtype  PSD_monitor
#define mccompcurindex 38
#define PSD_N mccPSD_afterfast2_PSD_N
#define PSD_p mccPSD_afterfast2_PSD_p
#define PSD_p2 mccPSD_afterfast2_PSD_p2
{   /* Declarations of PSD_afterfast2=PSD_monitor() SETTING parameters. */
int nx = mccPSD_afterfast2_nx;
int ny = mccPSD_afterfast2_ny;
char* filename = mccPSD_afterfast2_filename;
MCNUM xmin = mccPSD_afterfast2_xmin;
MCNUM xmax = mccPSD_afterfast2_xmax;
MCNUM ymin = mccPSD_afterfast2_ymin;
MCNUM ymax = mccPSD_afterfast2_ymax;
MCNUM xwidth = mccPSD_afterfast2_xwidth;
MCNUM yheight = mccPSD_afterfast2_yheight;
MCNUM restore_neutron = mccPSD_afterfast2_restore_neutron;
#line 129 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
{
  multiline(5,
    (double)xmin, (double)ymin, 0.0,
    (double)xmax, (double)ymin, 0.0,
    (double)xmax, (double)ymax, 0.0,
    (double)xmin, (double)ymax, 0.0,
    (double)xmin, (double)ymin, 0.0);
}
#line 29270 "./ESS_IN5_reprate.c"
}   /* End of PSD_afterfast2=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'Guidesample'. */
  SIG_MESSAGE("Guidesample (McDisplay)");
  printf("MCDISPLAY: component %s\n", "Guidesample");
#define mccompcurname  Guidesample
#define mccompcurtype  Guide
#define mccompcurindex 39
#define pTable mccGuidesample_pTable
{   /* Declarations of Guidesample=Guide() SETTING parameters. */
char* reflect = mccGuidesample_reflect;
MCNUM w1 = mccGuidesample_w1;
MCNUM h1 = mccGuidesample_h1;
MCNUM w2 = mccGuidesample_w2;
MCNUM h2 = mccGuidesample_h2;
MCNUM l = mccGuidesample_l;
MCNUM R0 = mccGuidesample_R0;
MCNUM Qc = mccGuidesample_Qc;
MCNUM alpha = mccGuidesample_alpha;
MCNUM m = mccGuidesample_m;
MCNUM W = mccGuidesample_W;
#line 201 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide.comp"
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
#line 29318 "./ESS_IN5_reprate.c"
}   /* End of Guidesample=Guide() SETTING parameter declarations. */
#undef pTable
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'Lmon_guideend'. */
  SIG_MESSAGE("Lmon_guideend (McDisplay)");
  printf("MCDISPLAY: component %s\n", "Lmon_guideend");
#define mccompcurname  Lmon_guideend
#define mccompcurtype  L_monitor
#define mccompcurindex 40
#define nL mccLmon_guideend_nL
#define L_N mccLmon_guideend_L_N
#define L_p mccLmon_guideend_L_p
#define L_p2 mccLmon_guideend_L_p2
{   /* Declarations of Lmon_guideend=L_monitor() SETTING parameters. */
char* filename = mccLmon_guideend_filename;
MCNUM xmin = mccLmon_guideend_xmin;
MCNUM xmax = mccLmon_guideend_xmax;
MCNUM ymin = mccLmon_guideend_ymin;
MCNUM ymax = mccLmon_guideend_ymax;
MCNUM xwidth = mccLmon_guideend_xwidth;
MCNUM yheight = mccLmon_guideend_yheight;
MCNUM Lmin = mccLmon_guideend_Lmin;
MCNUM Lmax = mccLmon_guideend_Lmax;
MCNUM restore_neutron = mccLmon_guideend_restore_neutron;
int nowritefile = mccLmon_guideend_nowritefile;
#line 120 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/L_monitor.comp"
{
  
  multiline(5, (double)xmin, (double)ymin, 0.0,
               (double)xmax, (double)ymin, 0.0,
               (double)xmax, (double)ymax, 0.0,
               (double)xmin, (double)ymax, 0.0,
               (double)xmin, (double)ymin, 0.0);
}
#line 29356 "./ESS_IN5_reprate.c"
}   /* End of Lmon_guideend=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'PSDsample'. */
  SIG_MESSAGE("PSDsample (McDisplay)");
  printf("MCDISPLAY: component %s\n", "PSDsample");
#define mccompcurname  PSDsample
#define mccompcurtype  PSD_monitor
#define mccompcurindex 41
#define PSD_N mccPSDsample_PSD_N
#define PSD_p mccPSDsample_PSD_p
#define PSD_p2 mccPSDsample_PSD_p2
{   /* Declarations of PSDsample=PSD_monitor() SETTING parameters. */
int nx = mccPSDsample_nx;
int ny = mccPSDsample_ny;
char* filename = mccPSDsample_filename;
MCNUM xmin = mccPSDsample_xmin;
MCNUM xmax = mccPSDsample_xmax;
MCNUM ymin = mccPSDsample_ymin;
MCNUM ymax = mccPSDsample_ymax;
MCNUM xwidth = mccPSDsample_xwidth;
MCNUM yheight = mccPSDsample_yheight;
MCNUM restore_neutron = mccPSDsample_restore_neutron;
#line 129 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
{
  multiline(5,
    (double)xmin, (double)ymin, 0.0,
    (double)xmax, (double)ymin, 0.0,
    (double)xmax, (double)ymax, 0.0,
    (double)xmin, (double)ymax, 0.0,
    (double)xmin, (double)ymin, 0.0);
}
#line 29395 "./ESS_IN5_reprate.c"
}   /* End of PSDsample=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'TOFsample_zoom'. */
  SIG_MESSAGE("TOFsample_zoom (McDisplay)");
  printf("MCDISPLAY: component %s\n", "TOFsample_zoom");
#define mccompcurname  TOFsample_zoom
#define mccompcurtype  TOF_monitor
#define mccompcurindex 42
#define nt mccTOFsample_zoom_nt
#define TOF_N mccTOFsample_zoom_TOF_N
#define TOF_p mccTOFsample_zoom_TOF_p
#define TOF_p2 mccTOFsample_zoom_TOF_p2
#define t_min mccTOFsample_zoom_t_min
#define t_max mccTOFsample_zoom_t_max
#define delta_t mccTOFsample_zoom_delta_t
{   /* Declarations of TOFsample_zoom=TOF_monitor() SETTING parameters. */
char* filename = mccTOFsample_zoom_filename;
MCNUM xmin = mccTOFsample_zoom_xmin;
MCNUM xmax = mccTOFsample_zoom_xmax;
MCNUM ymin = mccTOFsample_zoom_ymin;
MCNUM ymax = mccTOFsample_zoom_ymax;
MCNUM xwidth = mccTOFsample_zoom_xwidth;
MCNUM yheight = mccTOFsample_zoom_yheight;
MCNUM tmin = mccTOFsample_zoom_tmin;
MCNUM tmax = mccTOFsample_zoom_tmax;
MCNUM dt = mccTOFsample_zoom_dt;
MCNUM restore_neutron = mccTOFsample_zoom_restore_neutron;
int nowritefile = mccTOFsample_zoom_nowritefile;
#line 128 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/TOF_monitor.comp"
{
  
  multiline(5, (double)xmin, (double)ymin, 0.0,
               (double)xmax, (double)ymin, 0.0,
               (double)xmax, (double)ymax, 0.0,
               (double)xmin, (double)ymax, 0.0,
               (double)xmin, (double)ymin, 0.0);
}
#line 29439 "./ESS_IN5_reprate.c"
}   /* End of TOFsample_zoom=TOF_monitor() SETTING parameter declarations. */
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

  /* MCDISPLAY code for component 'Esample'. */
  SIG_MESSAGE("Esample (McDisplay)");
  printf("MCDISPLAY: component %s\n", "Esample");
#define mccompcurname  Esample
#define mccompcurtype  E_monitor
#define mccompcurindex 43
#define nE mccEsample_nE
#define E_N mccEsample_E_N
#define E_p mccEsample_E_p
#define E_p2 mccEsample_E_p2
#define S_p mccEsample_S_p
#define S_pE mccEsample_S_pE
#define S_pE2 mccEsample_S_pE2
{   /* Declarations of Esample=E_monitor() SETTING parameters. */
char* filename = mccEsample_filename;
MCNUM xmin = mccEsample_xmin;
MCNUM xmax = mccEsample_xmax;
MCNUM ymin = mccEsample_ymin;
MCNUM ymax = mccEsample_ymax;
MCNUM xwidth = mccEsample_xwidth;
MCNUM yheight = mccEsample_yheight;
MCNUM Emin = mccEsample_Emin;
MCNUM Emax = mccEsample_Emax;
MCNUM restore_neutron = mccEsample_restore_neutron;
int nowritefile = mccEsample_nowritefile;
#line 132 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/E_monitor.comp"
{
  
  multiline(5, (double)xmin, (double)ymin, 0.0,
               (double)xmax, (double)ymin, 0.0,
               (double)xmax, (double)ymax, 0.0,
               (double)xmin, (double)ymax, 0.0,
               (double)xmin, (double)ymin, 0.0);
}
#line 29486 "./ESS_IN5_reprate.c"
}   /* End of Esample=E_monitor() SETTING parameter declarations. */
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

  /* MCDISPLAY code for component 'Lmon_sample_zoom'. */
  SIG_MESSAGE("Lmon_sample_zoom (McDisplay)");
  printf("MCDISPLAY: component %s\n", "Lmon_sample_zoom");
#define mccompcurname  Lmon_sample_zoom
#define mccompcurtype  L_monitor
#define mccompcurindex 44
#define nL mccLmon_sample_zoom_nL
#define L_N mccLmon_sample_zoom_L_N
#define L_p mccLmon_sample_zoom_L_p
#define L_p2 mccLmon_sample_zoom_L_p2
{   /* Declarations of Lmon_sample_zoom=L_monitor() SETTING parameters. */
char* filename = mccLmon_sample_zoom_filename;
MCNUM xmin = mccLmon_sample_zoom_xmin;
MCNUM xmax = mccLmon_sample_zoom_xmax;
MCNUM ymin = mccLmon_sample_zoom_ymin;
MCNUM ymax = mccLmon_sample_zoom_ymax;
MCNUM xwidth = mccLmon_sample_zoom_xwidth;
MCNUM yheight = mccLmon_sample_zoom_yheight;
MCNUM Lmin = mccLmon_sample_zoom_Lmin;
MCNUM Lmax = mccLmon_sample_zoom_Lmax;
MCNUM restore_neutron = mccLmon_sample_zoom_restore_neutron;
int nowritefile = mccLmon_sample_zoom_nowritefile;
#line 120 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/L_monitor.comp"
{
  
  multiline(5, (double)xmin, (double)ymin, 0.0,
               (double)xmax, (double)ymin, 0.0,
               (double)xmax, (double)ymax, 0.0,
               (double)xmin, (double)ymax, 0.0,
               (double)xmin, (double)ymin, 0.0);
}
#line 29530 "./ESS_IN5_reprate.c"
}   /* End of Lmon_sample_zoom=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'sample'. */
  SIG_MESSAGE("sample (McDisplay)");
  printf("MCDISPLAY: component %s\n", "sample");
#define mccompcurname  sample
#define mccompcurtype  Tunneling_sample
#define mccompcurindex 45
#define ftun mccsample_ftun
#define fQE mccsample_fQE
#define VarsV mccsample_VarsV
{   /* Declarations of sample=Tunneling_sample() SETTING parameters. */
MCNUM thickness = mccsample_thickness;
MCNUM radius = mccsample_radius;
MCNUM focus_r = mccsample_focus_r;
MCNUM p_interact = mccsample_p_interact;
MCNUM f_QE = mccsample_f_QE;
MCNUM f_tun = mccsample_f_tun;
MCNUM gamma = mccsample_gamma;
MCNUM E_tun = mccsample_E_tun;
MCNUM target_x = mccsample_target_x;
MCNUM target_y = mccsample_target_y;
MCNUM target_z = mccsample_target_z;
MCNUM focus_xw = mccsample_focus_xw;
MCNUM focus_yh = mccsample_focus_yh;
MCNUM focus_aw = mccsample_focus_aw;
MCNUM focus_ah = mccsample_focus_ah;
MCNUM xwidth = mccsample_xwidth;
MCNUM yheight = mccsample_yheight;
MCNUM zdepth = mccsample_zdepth;
MCNUM sigma_abs = mccsample_sigma_abs;
MCNUM sigma_inc = mccsample_sigma_inc;
MCNUM Vc = mccsample_Vc;
int target_index = mccsample_target_index;
#line 321 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../samples/Tunneling_sample.comp"
{
  
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
}
#line 29618 "./ESS_IN5_reprate.c"
}   /* End of sample=Tunneling_sample() SETTING parameter declarations. */
#undef VarsV
#undef fQE
#undef ftun
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'detectorarm'. */
  SIG_MESSAGE("detectorarm (McDisplay)");
  printf("MCDISPLAY: component %s\n", "detectorarm");
#define mccompcurname  detectorarm
#define mccompcurtype  Arm
#define mccompcurindex 46
#line 40 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Arm.comp"
{
  /* A bit ugly; hard-coded dimensions. */
  
  line(0,0,0,0.2,0,0);
  line(0,0,0,0,0.2,0);
  line(0,0,0,0,0,0.2);
}
#line 29641 "./ESS_IN5_reprate.c"
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'TOFdetector'. */
  SIG_MESSAGE("TOFdetector (McDisplay)");
  printf("MCDISPLAY: component %s\n", "TOFdetector");
#define mccompcurname  TOFdetector
#define mccompcurtype  TOF_monitor
#define mccompcurindex 47
#define nt mccTOFdetector_nt
#define TOF_N mccTOFdetector_TOF_N
#define TOF_p mccTOFdetector_TOF_p
#define TOF_p2 mccTOFdetector_TOF_p2
#define t_min mccTOFdetector_t_min
#define t_max mccTOFdetector_t_max
#define delta_t mccTOFdetector_delta_t
{   /* Declarations of TOFdetector=TOF_monitor() SETTING parameters. */
char* filename = mccTOFdetector_filename;
MCNUM xmin = mccTOFdetector_xmin;
MCNUM xmax = mccTOFdetector_xmax;
MCNUM ymin = mccTOFdetector_ymin;
MCNUM ymax = mccTOFdetector_ymax;
MCNUM xwidth = mccTOFdetector_xwidth;
MCNUM yheight = mccTOFdetector_yheight;
MCNUM tmin = mccTOFdetector_tmin;
MCNUM tmax = mccTOFdetector_tmax;
MCNUM dt = mccTOFdetector_dt;
MCNUM restore_neutron = mccTOFdetector_restore_neutron;
int nowritefile = mccTOFdetector_nowritefile;
#line 128 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/TOF_monitor.comp"
{
  
  multiline(5, (double)xmin, (double)ymin, 0.0,
               (double)xmax, (double)ymin, 0.0,
               (double)xmax, (double)ymax, 0.0,
               (double)xmin, (double)ymax, 0.0,
               (double)xmin, (double)ymin, 0.0);
}
#line 29681 "./ESS_IN5_reprate.c"
}   /* End of TOFdetector=TOF_monitor() SETTING parameter declarations. */
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

  /* MCDISPLAY code for component 'TOFdetector_zoom'. */
  SIG_MESSAGE("TOFdetector_zoom (McDisplay)");
  printf("MCDISPLAY: component %s\n", "TOFdetector_zoom");
#define mccompcurname  TOFdetector_zoom
#define mccompcurtype  TOF_monitor
#define mccompcurindex 48
#define nt mccTOFdetector_zoom_nt
#define TOF_N mccTOFdetector_zoom_TOF_N
#define TOF_p mccTOFdetector_zoom_TOF_p
#define TOF_p2 mccTOFdetector_zoom_TOF_p2
#define t_min mccTOFdetector_zoom_t_min
#define t_max mccTOFdetector_zoom_t_max
#define delta_t mccTOFdetector_zoom_delta_t
{   /* Declarations of TOFdetector_zoom=TOF_monitor() SETTING parameters. */
char* filename = mccTOFdetector_zoom_filename;
MCNUM xmin = mccTOFdetector_zoom_xmin;
MCNUM xmax = mccTOFdetector_zoom_xmax;
MCNUM ymin = mccTOFdetector_zoom_ymin;
MCNUM ymax = mccTOFdetector_zoom_ymax;
MCNUM xwidth = mccTOFdetector_zoom_xwidth;
MCNUM yheight = mccTOFdetector_zoom_yheight;
MCNUM tmin = mccTOFdetector_zoom_tmin;
MCNUM tmax = mccTOFdetector_zoom_tmax;
MCNUM dt = mccTOFdetector_zoom_dt;
MCNUM restore_neutron = mccTOFdetector_zoom_restore_neutron;
int nowritefile = mccTOFdetector_zoom_nowritefile;
#line 128 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/TOF_monitor.comp"
{
  
  multiline(5, (double)xmin, (double)ymin, 0.0,
               (double)xmax, (double)ymin, 0.0,
               (double)xmax, (double)ymax, 0.0,
               (double)xmin, (double)ymax, 0.0,
               (double)xmin, (double)ymin, 0.0);
}
#line 29729 "./ESS_IN5_reprate.c"
}   /* End of TOFdetector_zoom=TOF_monitor() SETTING parameter declarations. */
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

  /* MCDISPLAY code for component 'Edetector'. */
  SIG_MESSAGE("Edetector (McDisplay)");
  printf("MCDISPLAY: component %s\n", "Edetector");
#define mccompcurname  Edetector
#define mccompcurtype  E_monitor
#define mccompcurindex 49
#define nE mccEdetector_nE
#define E_N mccEdetector_E_N
#define E_p mccEdetector_E_p
#define E_p2 mccEdetector_E_p2
#define S_p mccEdetector_S_p
#define S_pE mccEdetector_S_pE
#define S_pE2 mccEdetector_S_pE2
{   /* Declarations of Edetector=E_monitor() SETTING parameters. */
char* filename = mccEdetector_filename;
MCNUM xmin = mccEdetector_xmin;
MCNUM xmax = mccEdetector_xmax;
MCNUM ymin = mccEdetector_ymin;
MCNUM ymax = mccEdetector_ymax;
MCNUM xwidth = mccEdetector_xwidth;
MCNUM yheight = mccEdetector_yheight;
MCNUM Emin = mccEdetector_Emin;
MCNUM Emax = mccEdetector_Emax;
MCNUM restore_neutron = mccEdetector_restore_neutron;
int nowritefile = mccEdetector_nowritefile;
#line 132 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/E_monitor.comp"
{
  
  multiline(5, (double)xmin, (double)ymin, 0.0,
               (double)xmax, (double)ymin, 0.0,
               (double)xmax, (double)ymax, 0.0,
               (double)xmin, (double)ymax, 0.0,
               (double)xmin, (double)ymin, 0.0);
}
#line 29776 "./ESS_IN5_reprate.c"
}   /* End of Edetector=E_monitor() SETTING parameter declarations. */
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

  /* MCDISPLAY code for component 'TOF2Edetector'. */
  SIG_MESSAGE("TOF2Edetector (McDisplay)");
  printf("MCDISPLAY: component %s\n", "TOF2Edetector");
#define mccompcurname  TOF2Edetector
#define mccompcurtype  TOF2E_monitor
#define mccompcurindex 50
#define nE mccTOF2Edetector_nE
#define E_N mccTOF2Edetector_E_N
#define E_p mccTOF2Edetector_E_p
#define E_p2 mccTOF2Edetector_E_p2
#define S_p mccTOF2Edetector_S_p
#define S_pE mccTOF2Edetector_S_pE
#define S_pE2 mccTOF2Edetector_S_pE2
{   /* Declarations of TOF2Edetector=TOF2E_monitor() SETTING parameters. */
char* filename = mccTOF2Edetector_filename;
MCNUM xmin = mccTOF2Edetector_xmin;
MCNUM xmax = mccTOF2Edetector_xmax;
MCNUM ymin = mccTOF2Edetector_ymin;
MCNUM ymax = mccTOF2Edetector_ymax;
MCNUM xwidth = mccTOF2Edetector_xwidth;
MCNUM yheight = mccTOF2Edetector_yheight;
MCNUM Emin = mccTOF2Edetector_Emin;
MCNUM Emax = mccTOF2Edetector_Emax;
MCNUM T_zero = mccTOF2Edetector_T_zero;
MCNUM L_flight = mccTOF2Edetector_L_flight;
MCNUM restore_neutron = mccTOF2Edetector_restore_neutron;
int nowritefile = mccTOF2Edetector_nowritefile;
#line 133 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/TOF2E_monitor.comp"
{
  
  multiline(5, (double)xmin, (double)ymin, 0.0,
               (double)xmax, (double)ymin, 0.0,
               (double)xmax, (double)ymax, 0.0,
               (double)xmin, (double)ymax, 0.0,
               (double)xmin, (double)ymin, 0.0);
}
#line 29825 "./ESS_IN5_reprate.c"
}   /* End of TOF2Edetector=TOF2E_monitor() SETTING parameter declarations. */
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
/* end of generated C code ./ESS_IN5_reprate.c */
