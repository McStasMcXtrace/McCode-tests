/* Automatically generated file. Do not edit. 
 * Format:     ANSI C source code
 * Creator:    McStas <http://www.mcstas.org>
 * Instrument: /zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr (ISIS_SANS2d)
 * Date:       Wed Feb 26 19:10:55 2020
 * File:       ./ISIS_SANS2d.c
 * Compile:    cc -o ISIS_SANS2d.out ./ISIS_SANS2d.c 
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

#line 712 "./ISIS_SANS2d.c"

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

#line 945 "./ISIS_SANS2d.c"

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

#line 4977 "./ISIS_SANS2d.c"

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

#line 5337 "./ISIS_SANS2d.c"
#ifdef MC_TRACE_ENABLED
int mctraceenabled = 1;
#else
int mctraceenabled = 0;
#endif
#define MCSTAS "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../"
int mcdefaultmain = 1;
char mcinstrument_name[] = "ISIS_SANS2d";
char mcinstrument_source[] = "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr";
char *mcinstrument_exe=NULL; /* will be set to argv[0] in main */
int main(int argc, char *argv[]){return mccode_main(argc, argv);}
void mcinit(void);
void mcraytrace(void);
void mcsave(FILE *);
void mcfinally(void);
void mcdisplay(void);

/* Shared user declarations for all components 'ISIS_moderator'. */
#line 64 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../contrib/ISIS_moderator.comp"
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
#line 6104 "./ISIS_SANS2d.c"

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

#line 6190 "./ISIS_SANS2d.c"

/* Shared user declarations for all components 'Guide_gravity'. */
#line 124 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide_gravity.comp"
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

#ifndef Gravity_guide_Version
#define Gravity_guide_Version "$Revision$"

#ifndef PROP_GRAV_DT
#error McStas : You need PROP_GRAV_DT (McStas >= 1.4.3) to run this component
#endif

/*
* G:       (m/s^2) Gravitation acceleration along y axis [-9.81]
* Gx:      (m/s^2) Gravitation acceleration along x axis [0]
* Gy:      (m/s^2) Gravitation acceleration along y axis [-9.81]
* Gz:      (m/s^2) Gravitation acceleration along z axis [0]
* mh:      (1)    m-value of material for left/right vert. mirrors
* mv:      (1)    m-value of material for top/bottom horz. mirrors
* mx:      (1)    m-value of material for left/right vert. mirrors
* my:      (1)    m-value of material for top/bottom horz. mirrors
*/

  typedef struct Gravity_guide_Vars
  {
    double gx;
    double gy;
    double gz;
    double nx[6], ny[6], nz[6];
    double wx[6], wy[6], wz[6];
    double A[6], norm_n2[6], norm_n[6];
    long   N_reflection[7];
    double w1c, h1c;
    double w2c, h2c;
    double M[5];
    double Alpha[5];
    double nzC[5], norm_n2xy[5], Axy[5];
    double wav_lr, wav_tb, wav_z;
    double chamfer_z, chamfer_lr, chamfer_tb;
    char   compcurname[256];
    double fc_freq, fc_phase;
    double warnings;
  } Gravity_guide_Vars_type;

  void Gravity_guide_Init(Gravity_guide_Vars_type *aVars,
    MCNUM a_w1, MCNUM a_h1, MCNUM a_w2, MCNUM a_h2, MCNUM a_l, MCNUM a_R0,
    MCNUM a_Qc, MCNUM a_alpha, MCNUM a_m, MCNUM a_W, MCNUM a_nslit, MCNUM a_d,
    MCNUM a_Gx, MCNUM a_Gy, MCNUM a_Gz,
    MCNUM a_mleft, MCNUM a_mright, MCNUM a_mtop, MCNUM a_mbottom, MCNUM a_nhslit,
    MCNUM a_wavy_lr, MCNUM a_wavy_tb, MCNUM a_wavy_z, MCNUM a_wavy,
    MCNUM a_chamfers_z, MCNUM a_chamfers_lr, MCNUM a_chamfers_tb, MCNUM a_chamfers,
    MCNUM a_nu, MCNUM a_phase, MCNUM a_aleft, MCNUM a_aright, MCNUM a_atop, MCNUM a_abottom)
  {
    int i;

    for (i=0; i<7; aVars->N_reflection[i++] = 0);
    for (i=0; i<5; aVars->M[i++] = 0);
    for (i=0; i<5; aVars->Alpha[i++] = 0);

    aVars->gx = a_Gx; /* The gravitation vector in the current component axis system */
    aVars->gy = a_Gy;
    aVars->gz = a_Gz;
    aVars->warnings=0;

    if (a_nslit <= 0 || a_nhslit <= 0) { fprintf(stderr,"%s: Fatal: no channel in this guide (nhslit or nslit=0).\n", aVars->compcurname); exit(-1); }
    if (a_d < 0) { fprintf(stderr,"%s: Fatal: subdividing walls have negative thickness in this guide (d<0).\n", aVars->compcurname); exit(-1); }
    aVars->w1c = (a_w1 - (a_nslit-1) *a_d)/(double)a_nslit;
    aVars->w2c = (a_w2 - (a_nslit-1) *a_d)/(double)a_nslit;
    aVars->h1c = (a_h1 - (a_nhslit-1)*a_d)/(double)a_nhslit;
    aVars->h2c = (a_h2 - (a_nhslit-1)*a_d)/(double)a_nhslit;

    for (i=0; i <= 4;   aVars->M[i++]=a_m);
    for (i=0; i <= 4;   aVars->Alpha[i++]=a_alpha);
    if (a_mleft   >= 0) aVars->M[1] =a_mleft  ;
    if (a_mright  >= 0) aVars->M[2] =a_mright ;
    if (a_mtop    >= 0) aVars->M[3] =a_mtop   ;
    if (a_mbottom >= 0) aVars->M[4] =a_mbottom;
    if (a_aleft   >= 0) aVars->Alpha[1] =a_aleft  ;
    if (a_aright  >= 0) aVars->Alpha[2] =a_aright ;
    if (a_atop    >= 0) aVars->Alpha[3] =a_atop   ;
    if (a_abottom >= 0) aVars->Alpha[4] =a_abottom;

    /* n: normal vectors to surfaces */
    aVars->nx[1] =  a_l; aVars->ny[1] =  0;   aVars->nz[1] =  0.5*(aVars->w2c-aVars->w1c);  /* 1:+X left       */
    aVars->nx[2] = -a_l; aVars->ny[2] =  0;   aVars->nz[2] = -aVars->nz[1];             /* 2:-X right      */
    aVars->nx[3] =  0;   aVars->ny[3] =  a_l; aVars->nz[3] =  0.5*(aVars->h2c-aVars->h1c);  /* 3:+Y top        */
    aVars->nx[4] =  0;   aVars->ny[4] = -a_l; aVars->nz[4] = -aVars->nz[3];             /* 4:-Y bottom     */
    aVars->nx[5] =  0;   aVars->ny[5] =  0;   aVars->nz[5] =  a_l;                      /* 5:+Z exit       */
    aVars->nx[0] =  0;   aVars->ny[0] =  0;   aVars->nz[0] = -a_l;                      /* 0:Z0 input      */
    /* w: a point on these surfaces */
    aVars->wx[1] = +(aVars->w1c)/2; aVars->wy[1] =  0;              aVars->wz[1] = 0;   /* 1:+X left       */
    aVars->wx[2] = -(aVars->w1c)/2; aVars->wy[2] =  0;              aVars->wz[2] = 0;   /* 2:-X right      */
    aVars->wx[3] =  0;              aVars->wy[3] = +(aVars->h1c)/2; aVars->wz[3] = 0;   /* 3:+Y top        */
    aVars->wx[4] =  0;              aVars->wy[4] = -(aVars->h1c)/2; aVars->wz[4] = 0;   /* 4:-Y bottom     */
    aVars->wx[5] =  0;              aVars->wy[5] =  0;              aVars->wz[5] = a_l; /* 5:+Z exit       */
    aVars->wx[0] =  0;              aVars->wy[0] =  0;              aVars->wz[0] = 0;   /* 0:Z0 input      */

    for (i=0; i <= 5; i++)
    {
      aVars->A[i] = scalar_prod(aVars->nx[i], aVars->ny[i], aVars->nz[i], aVars->gx, aVars->gy, aVars->gz)/2;
      aVars->norm_n2[i] = aVars->nx[i]*aVars->nx[i] + aVars->ny[i]*aVars->ny[i] + aVars->nz[i]*aVars->nz[i];
      if (aVars->norm_n2[i] <= 0)
        { fprintf(stderr,"%s: Fatal: normal vector norm %i is null/negative ! check guide dimensions.\n", aVars->compcurname, i); exit(-1); } /* should never occur */
      else
        aVars->norm_n[i] = sqrt(aVars->norm_n2[i]);
    }
    /* partial computations for l/r/t/b sides, to save computing time */
    for (i=1; i <= 4; i++)
    { /* stores nz that changes in case non box element (focus/defocus) */
      aVars->nzC[i]      =  aVars->nz[i]; /* partial xy terms */
      aVars->norm_n2xy[i]=  aVars->nx[i]*aVars->nx[i] + aVars->ny[i]*aVars->ny[i];
      aVars->Axy[i]      = (aVars->nx[i]*aVars->gx    + aVars->ny[i]*aVars->gy)/2;
    }
    /* handle waviness init */
    if (a_wavy && (!a_wavy_tb && !a_wavy_lr && !a_wavy_z))
    { aVars->wav_tb=aVars->wav_lr=aVars->wav_z=a_wavy; }
    else
    { aVars->wav_tb=a_wavy_tb; aVars->wav_lr=a_wavy_lr; aVars->wav_z=a_wavy_z; }
    aVars->wav_tb *= DEG2RAD/(sqrt(8*log(2)));   /* Convert from deg FWHM to rad Gaussian sigma */
    aVars->wav_lr *= DEG2RAD/(sqrt(8*log(2)));
    aVars->wav_z  *= DEG2RAD/(sqrt(8*log(2)));
    /* handle chamfers init */
    if (a_chamfers && (!a_chamfers_z && !a_chamfers_lr && !a_chamfers_tb))
    { aVars->chamfer_z=aVars->chamfer_lr=aVars->chamfer_tb=a_chamfers; }
    else
    {
      aVars->chamfer_z=a_chamfers_z;
      aVars->chamfer_lr=a_chamfers_lr;
      aVars->chamfer_tb=a_chamfers_tb;
    }

    aVars->fc_freq  = a_nu;
    aVars->fc_phase = a_phase;
  }

  int Gravity_guide_Trace(double *dt,
        Gravity_guide_Vars_type *aVars,
        double cx, double cy, double cz,
        double cvx, double cvy, double cvz,
        double cxnum, double cxk, double cynum, double cyk,
        double *cnx, double *cny,double *cnz)
  {
    double B, C;
    int    ret=0;
    int    side=0;
    double n1;
    double dt0, dt_min=0;
    int    i;
    double loc_num, loc_nslit;
    int    i_slope=3;

    /* look if there is a previous intersection with guide sides */
    /* A = 0.5 n.g; B = n.v; C = n.(r-W); */
    /* 5=+Z side: n=(0, 0, -l) ; W = (0, 0, l) (at z=l, guide exit)*/
    B = aVars->nz[5]*cvz; C = aVars->nz[5]*(cz - aVars->wz[5]);
    ret = solve_2nd_order(&dt0, NULL, aVars->A[5], B, C);
    if (ret && dt0>1e-10) { dt_min = dt0; side=5; }

    loc_num = cynum; loc_nslit = cyk;
    for (i=4; i>0; i--)
    {
      if (i == 2) { i_slope=1; loc_num = cxnum; loc_nslit = cxk; }

      if (aVars->nzC[i_slope] != 0) {
        n1 = loc_nslit - 2*(loc_num);  /* slope of l/r/u/d sides depends on the channel ! */
        loc_num++; /* use partial computations to alter nz and A */
        aVars->nz[i]= aVars->nzC[i]*n1;
        aVars->A[i] = aVars->Axy[i] + aVars->nz[i]*aVars->gz/2;
      }
      if (i < 3)
      {      B = aVars->nx[i]*cvx + aVars->nz[i]*cvz; C = aVars->nx[i]*(cx-aVars->wx[i]) + aVars->nz[i]*cz; }
      else { B = aVars->ny[i]*cvy + aVars->nz[i]*cvz; C = aVars->ny[i]*(cy-aVars->wy[i]) + aVars->nz[i]*cz; }
      ret = solve_2nd_order(&dt0, NULL, aVars->A[i], B, C);
      if (ret && dt0>1e-10 && (dt0<dt_min || !dt_min))
      { dt_min = dt0; side=i;
        if (aVars->nzC[i] != 0)
        { aVars->norm_n2[i] = aVars->norm_n2xy[i] + aVars->nz[i]*aVars->nz[i];
          aVars->norm_n[i]  = sqrt(aVars->norm_n2[i]); }
      }
     }

    *dt = dt_min;
    /* handles waviness: rotate n vector */
    if (side > 0 && side < 5 && (aVars->wav_z || aVars->wav_lr || aVars->wav_tb))
    {
      double nt_x, nt_y, nt_z;  /* transverse vector */
      double nn_x, nn_y, nn_z;  /* normal vector (tmp) */
      double phi;
      /* normal vector n_z = [ 0,0,1], n_t = n x n_z; */
      vec_prod(nt_x,nt_y,nt_z, aVars->nx[side],aVars->ny[side],aVars->nz[side], 0,0,1);
      /* rotate n with angle wavy_z around n_t -> nn */
      if (aVars->wav_z) {
        phi = aVars->wav_z;
        rotate(nn_x,nn_y,nn_z, aVars->nx[side],aVars->ny[side],aVars->nz[side], aVars->wav_z*randnorm(), nt_x,nt_y,nt_z);
      } else { nn_x=aVars->nx[side]; nn_y=aVars->ny[side]; nn_z=aVars->nz[side]; }
      /* rotate n with angle wavy_{x|y} around n_z -> nt */
      phi = (side <=2) ? aVars->wav_lr : aVars->wav_tb;
      if (phi) {
        rotate(nt_x,nt_y,nt_z, nn_x,nn_y,nn_z, phi*randnorm(), 0,0,1);
      } else { nt_x=nn_x; nt_y=nn_y; nt_z=nn_z; }
      *cnx=nt_x; *cny=nt_y; *cnz=nt_z;
    } else
    { *cnx=aVars->nx[side]; *cny=aVars->ny[side]; *cnz=aVars->nz[side]; }
    return (side);
  }



#endif
#line 7988 "./ISIS_SANS2d.c"

/* Instrument parameters. */
MCNUM mcipL1;
MCNUM mcipA1w;
MCNUM mcipA1h;
MCNUM mcipS6;
MCNUM mcipA2;
MCNUM mcipLmin;
MCNUM mcipLmax;

#define mcNUMIPAR 7
int mcnumipar = 7;
struct mcinputtable_struct mcinputtable[mcNUMIPAR+1] = {
  "L1", &mcipL1, instr_type_double, "3.926", 
  "A1w", &mcipA1w, instr_type_double, "0.030", 
  "A1h", &mcipA1h, instr_type_double, "0.02", 
  "S6", &mcipS6, instr_type_double, "0.006", 
  "A2", &mcipA2, instr_type_double, "0.006", 
  "Lmin", &mcipLmin, instr_type_double, "1", 
  "Lmax", &mcipLmax, instr_type_double, "14", 
  NULL, NULL, instr_type_double, ""
};

/* User declarations from instrument definition. */
#define mccompcurname  ISIS_SANS2d
#define mccompcurtype  INSTRUMENT
#define mccompcurindex 0
#define mcposaISIS_SANS2d coords_set(0,0,0)
#define L1 mcipL1
#define A1w mcipA1w
#define A1h mcipA1h
#define S6 mcipS6
#define A2 mcipA2
#define Lmin mcipLmin
#define Lmax mcipLmax
#undef Lmax
#undef Lmin
#undef A2
#undef S6
#undef A1h
#undef A1w
#undef L1
#undef mcposaISIS_SANS2d
#undef mccompcurindex
#undef mccompcurtype
#undef mccompcurname

/* neutron state table at each component input (local coords) */
/* [x, y, z, vx, vy, vz, t, sx, sy, sz, p] */
MCNUM mccomp_storein[11*30];
/* Components position table (absolute and relative coords) */
Coords mccomp_posa[30];
Coords mccomp_posr[30];
/* Counter for each comp to check for inactive ones */
MCNUM  mcNCounter[30];
MCNUM  mcPCounter[30];
MCNUM  mcP2Counter[30];
#define mcNUMCOMP 29 /* number of components */
/* Counter for PROP ABSORB */
MCNUM  mcAbsorbProp[30];
/* Flag true when previous component acted on the neutron (SCATTER) */
MCNUM mcScattered=0;
/* Flag true when neutron should be restored (RESTORE) */
MCNUM mcRestore=0;
/* Declarations of component definition and setting parameters. */

/* Definition parameters for component 'isis_source' [2]. */
#define mccisis_source_Face "E2" /* declared as a string. May produce warnings at compile */
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

/* Definition parameters for component 'lmon1' [3]. */
#define mcclmon1_nL 140
/* Setting parameters for component 'lmon1' [3]. */
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

/* Setting parameters for component 'psd1' [4]. */
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

/* Setting parameters for component 'bender1' [5]. */
MCNUM mccbender1_w1;
MCNUM mccbender1_h1;
MCNUM mccbender1_w2;
MCNUM mccbender1_h2;
MCNUM mccbender1_l;
MCNUM mccbender1_R0;
MCNUM mccbender1_Qc;
MCNUM mccbender1_alpha;
MCNUM mccbender1_m;
MCNUM mccbender1_W;
MCNUM mccbender1_nslit;
MCNUM mccbender1_d;
MCNUM mccbender1_mleft;
MCNUM mccbender1_mright;
MCNUM mccbender1_mtop;
MCNUM mccbender1_mbottom;
MCNUM mccbender1_nhslit;
MCNUM mccbender1_G;
MCNUM mccbender1_aleft;
MCNUM mccbender1_aright;
MCNUM mccbender1_atop;
MCNUM mccbender1_abottom;
MCNUM mccbender1_wavy;
MCNUM mccbender1_wavy_z;
MCNUM mccbender1_wavy_tb;
MCNUM mccbender1_wavy_lr;
MCNUM mccbender1_chamfers;
MCNUM mccbender1_chamfers_z;
MCNUM mccbender1_chamfers_lr;
MCNUM mccbender1_chamfers_tb;
MCNUM mccbender1_nelements;
MCNUM mccbender1_nu;
MCNUM mccbender1_phase;
char mccbender1_reflect[16384];

/* Setting parameters for component 'bender2' [6]. */
MCNUM mccbender2_w1;
MCNUM mccbender2_h1;
MCNUM mccbender2_w2;
MCNUM mccbender2_h2;
MCNUM mccbender2_l;
MCNUM mccbender2_R0;
MCNUM mccbender2_Qc;
MCNUM mccbender2_alpha;
MCNUM mccbender2_m;
MCNUM mccbender2_W;
MCNUM mccbender2_nslit;
MCNUM mccbender2_d;
MCNUM mccbender2_mleft;
MCNUM mccbender2_mright;
MCNUM mccbender2_mtop;
MCNUM mccbender2_mbottom;
MCNUM mccbender2_nhslit;
MCNUM mccbender2_G;
MCNUM mccbender2_aleft;
MCNUM mccbender2_aright;
MCNUM mccbender2_atop;
MCNUM mccbender2_abottom;
MCNUM mccbender2_wavy;
MCNUM mccbender2_wavy_z;
MCNUM mccbender2_wavy_tb;
MCNUM mccbender2_wavy_lr;
MCNUM mccbender2_chamfers;
MCNUM mccbender2_chamfers_z;
MCNUM mccbender2_chamfers_lr;
MCNUM mccbender2_chamfers_tb;
MCNUM mccbender2_nelements;
MCNUM mccbender2_nu;
MCNUM mccbender2_phase;
char mccbender2_reflect[16384];

/* Setting parameters for component 'bender3' [7]. */
MCNUM mccbender3_w1;
MCNUM mccbender3_h1;
MCNUM mccbender3_w2;
MCNUM mccbender3_h2;
MCNUM mccbender3_l;
MCNUM mccbender3_R0;
MCNUM mccbender3_Qc;
MCNUM mccbender3_alpha;
MCNUM mccbender3_m;
MCNUM mccbender3_W;
MCNUM mccbender3_nslit;
MCNUM mccbender3_d;
MCNUM mccbender3_mleft;
MCNUM mccbender3_mright;
MCNUM mccbender3_mtop;
MCNUM mccbender3_mbottom;
MCNUM mccbender3_nhslit;
MCNUM mccbender3_G;
MCNUM mccbender3_aleft;
MCNUM mccbender3_aright;
MCNUM mccbender3_atop;
MCNUM mccbender3_abottom;
MCNUM mccbender3_wavy;
MCNUM mccbender3_wavy_z;
MCNUM mccbender3_wavy_tb;
MCNUM mccbender3_wavy_lr;
MCNUM mccbender3_chamfers;
MCNUM mccbender3_chamfers_z;
MCNUM mccbender3_chamfers_lr;
MCNUM mccbender3_chamfers_tb;
MCNUM mccbender3_nelements;
MCNUM mccbender3_nu;
MCNUM mccbender3_phase;
char mccbender3_reflect[16384];

/* Setting parameters for component 'bender4' [8]. */
MCNUM mccbender4_w1;
MCNUM mccbender4_h1;
MCNUM mccbender4_w2;
MCNUM mccbender4_h2;
MCNUM mccbender4_l;
MCNUM mccbender4_R0;
MCNUM mccbender4_Qc;
MCNUM mccbender4_alpha;
MCNUM mccbender4_m;
MCNUM mccbender4_W;
MCNUM mccbender4_nslit;
MCNUM mccbender4_d;
MCNUM mccbender4_mleft;
MCNUM mccbender4_mright;
MCNUM mccbender4_mtop;
MCNUM mccbender4_mbottom;
MCNUM mccbender4_nhslit;
MCNUM mccbender4_G;
MCNUM mccbender4_aleft;
MCNUM mccbender4_aright;
MCNUM mccbender4_atop;
MCNUM mccbender4_abottom;
MCNUM mccbender4_wavy;
MCNUM mccbender4_wavy_z;
MCNUM mccbender4_wavy_tb;
MCNUM mccbender4_wavy_lr;
MCNUM mccbender4_chamfers;
MCNUM mccbender4_chamfers_z;
MCNUM mccbender4_chamfers_lr;
MCNUM mccbender4_chamfers_tb;
MCNUM mccbender4_nelements;
MCNUM mccbender4_nu;
MCNUM mccbender4_phase;
char mccbender4_reflect[16384];

/* Setting parameters for component 'bender5' [9]. */
MCNUM mccbender5_w1;
MCNUM mccbender5_h1;
MCNUM mccbender5_w2;
MCNUM mccbender5_h2;
MCNUM mccbender5_l;
MCNUM mccbender5_R0;
MCNUM mccbender5_Qc;
MCNUM mccbender5_alpha;
MCNUM mccbender5_m;
MCNUM mccbender5_W;
MCNUM mccbender5_nslit;
MCNUM mccbender5_d;
MCNUM mccbender5_mleft;
MCNUM mccbender5_mright;
MCNUM mccbender5_mtop;
MCNUM mccbender5_mbottom;
MCNUM mccbender5_nhslit;
MCNUM mccbender5_G;
MCNUM mccbender5_aleft;
MCNUM mccbender5_aright;
MCNUM mccbender5_atop;
MCNUM mccbender5_abottom;
MCNUM mccbender5_wavy;
MCNUM mccbender5_wavy_z;
MCNUM mccbender5_wavy_tb;
MCNUM mccbender5_wavy_lr;
MCNUM mccbender5_chamfers;
MCNUM mccbender5_chamfers_z;
MCNUM mccbender5_chamfers_lr;
MCNUM mccbender5_chamfers_tb;
MCNUM mccbender5_nelements;
MCNUM mccbender5_nu;
MCNUM mccbender5_phase;
char mccbender5_reflect[16384];

/* Setting parameters for component 'bender6' [10]. */
MCNUM mccbender6_w1;
MCNUM mccbender6_h1;
MCNUM mccbender6_w2;
MCNUM mccbender6_h2;
MCNUM mccbender6_l;
MCNUM mccbender6_R0;
MCNUM mccbender6_Qc;
MCNUM mccbender6_alpha;
MCNUM mccbender6_m;
MCNUM mccbender6_W;
MCNUM mccbender6_nslit;
MCNUM mccbender6_d;
MCNUM mccbender6_mleft;
MCNUM mccbender6_mright;
MCNUM mccbender6_mtop;
MCNUM mccbender6_mbottom;
MCNUM mccbender6_nhslit;
MCNUM mccbender6_G;
MCNUM mccbender6_aleft;
MCNUM mccbender6_aright;
MCNUM mccbender6_atop;
MCNUM mccbender6_abottom;
MCNUM mccbender6_wavy;
MCNUM mccbender6_wavy_z;
MCNUM mccbender6_wavy_tb;
MCNUM mccbender6_wavy_lr;
MCNUM mccbender6_chamfers;
MCNUM mccbender6_chamfers_z;
MCNUM mccbender6_chamfers_lr;
MCNUM mccbender6_chamfers_tb;
MCNUM mccbender6_nelements;
MCNUM mccbender6_nu;
MCNUM mccbender6_phase;
char mccbender6_reflect[16384];

/* Setting parameters for component 'bender7' [11]. */
MCNUM mccbender7_w1;
MCNUM mccbender7_h1;
MCNUM mccbender7_w2;
MCNUM mccbender7_h2;
MCNUM mccbender7_l;
MCNUM mccbender7_R0;
MCNUM mccbender7_Qc;
MCNUM mccbender7_alpha;
MCNUM mccbender7_m;
MCNUM mccbender7_W;
MCNUM mccbender7_nslit;
MCNUM mccbender7_d;
MCNUM mccbender7_mleft;
MCNUM mccbender7_mright;
MCNUM mccbender7_mtop;
MCNUM mccbender7_mbottom;
MCNUM mccbender7_nhslit;
MCNUM mccbender7_G;
MCNUM mccbender7_aleft;
MCNUM mccbender7_aright;
MCNUM mccbender7_atop;
MCNUM mccbender7_abottom;
MCNUM mccbender7_wavy;
MCNUM mccbender7_wavy_z;
MCNUM mccbender7_wavy_tb;
MCNUM mccbender7_wavy_lr;
MCNUM mccbender7_chamfers;
MCNUM mccbender7_chamfers_z;
MCNUM mccbender7_chamfers_lr;
MCNUM mccbender7_chamfers_tb;
MCNUM mccbender7_nelements;
MCNUM mccbender7_nu;
MCNUM mccbender7_phase;
char mccbender7_reflect[16384];

/* Setting parameters for component 'bender8' [12]. */
MCNUM mccbender8_w1;
MCNUM mccbender8_h1;
MCNUM mccbender8_w2;
MCNUM mccbender8_h2;
MCNUM mccbender8_l;
MCNUM mccbender8_R0;
MCNUM mccbender8_Qc;
MCNUM mccbender8_alpha;
MCNUM mccbender8_m;
MCNUM mccbender8_W;
MCNUM mccbender8_nslit;
MCNUM mccbender8_d;
MCNUM mccbender8_mleft;
MCNUM mccbender8_mright;
MCNUM mccbender8_mtop;
MCNUM mccbender8_mbottom;
MCNUM mccbender8_nhslit;
MCNUM mccbender8_G;
MCNUM mccbender8_aleft;
MCNUM mccbender8_aright;
MCNUM mccbender8_atop;
MCNUM mccbender8_abottom;
MCNUM mccbender8_wavy;
MCNUM mccbender8_wavy_z;
MCNUM mccbender8_wavy_tb;
MCNUM mccbender8_wavy_lr;
MCNUM mccbender8_chamfers;
MCNUM mccbender8_chamfers_z;
MCNUM mccbender8_chamfers_lr;
MCNUM mccbender8_chamfers_tb;
MCNUM mccbender8_nelements;
MCNUM mccbender8_nu;
MCNUM mccbender8_phase;
char mccbender8_reflect[16384];

/* Setting parameters for component 'bender9' [13]. */
MCNUM mccbender9_w1;
MCNUM mccbender9_h1;
MCNUM mccbender9_w2;
MCNUM mccbender9_h2;
MCNUM mccbender9_l;
MCNUM mccbender9_R0;
MCNUM mccbender9_Qc;
MCNUM mccbender9_alpha;
MCNUM mccbender9_m;
MCNUM mccbender9_W;
MCNUM mccbender9_nslit;
MCNUM mccbender9_d;
MCNUM mccbender9_mleft;
MCNUM mccbender9_mright;
MCNUM mccbender9_mtop;
MCNUM mccbender9_mbottom;
MCNUM mccbender9_nhslit;
MCNUM mccbender9_G;
MCNUM mccbender9_aleft;
MCNUM mccbender9_aright;
MCNUM mccbender9_atop;
MCNUM mccbender9_abottom;
MCNUM mccbender9_wavy;
MCNUM mccbender9_wavy_z;
MCNUM mccbender9_wavy_tb;
MCNUM mccbender9_wavy_lr;
MCNUM mccbender9_chamfers;
MCNUM mccbender9_chamfers_z;
MCNUM mccbender9_chamfers_lr;
MCNUM mccbender9_chamfers_tb;
MCNUM mccbender9_nelements;
MCNUM mccbender9_nu;
MCNUM mccbender9_phase;
char mccbender9_reflect[16384];

/* Setting parameters for component 'bender10' [14]. */
MCNUM mccbender10_w1;
MCNUM mccbender10_h1;
MCNUM mccbender10_w2;
MCNUM mccbender10_h2;
MCNUM mccbender10_l;
MCNUM mccbender10_R0;
MCNUM mccbender10_Qc;
MCNUM mccbender10_alpha;
MCNUM mccbender10_m;
MCNUM mccbender10_W;
MCNUM mccbender10_nslit;
MCNUM mccbender10_d;
MCNUM mccbender10_mleft;
MCNUM mccbender10_mright;
MCNUM mccbender10_mtop;
MCNUM mccbender10_mbottom;
MCNUM mccbender10_nhslit;
MCNUM mccbender10_G;
MCNUM mccbender10_aleft;
MCNUM mccbender10_aright;
MCNUM mccbender10_atop;
MCNUM mccbender10_abottom;
MCNUM mccbender10_wavy;
MCNUM mccbender10_wavy_z;
MCNUM mccbender10_wavy_tb;
MCNUM mccbender10_wavy_lr;
MCNUM mccbender10_chamfers;
MCNUM mccbender10_chamfers_z;
MCNUM mccbender10_chamfers_lr;
MCNUM mccbender10_chamfers_tb;
MCNUM mccbender10_nelements;
MCNUM mccbender10_nu;
MCNUM mccbender10_phase;
char mccbender10_reflect[16384];

/* Definition parameters for component 'lmonb' [15]. */
#define mcclmonb_nL 140
/* Setting parameters for component 'lmonb' [15]. */
char mcclmonb_filename[16384];
MCNUM mcclmonb_xmin;
MCNUM mcclmonb_xmax;
MCNUM mcclmonb_ymin;
MCNUM mcclmonb_ymax;
MCNUM mcclmonb_xwidth;
MCNUM mcclmonb_yheight;
MCNUM mcclmonb_Lmin;
MCNUM mcclmonb_Lmax;
MCNUM mcclmonb_restore_neutron;
int mcclmonb_nowritefile;

/* Setting parameters for component 'psd2' [16]. */
int mccpsd2_nx;
int mccpsd2_ny;
char mccpsd2_filename[16384];
MCNUM mccpsd2_xmin;
MCNUM mccpsd2_xmax;
MCNUM mccpsd2_ymin;
MCNUM mccpsd2_ymax;
MCNUM mccpsd2_xwidth;
MCNUM mccpsd2_yheight;
MCNUM mccpsd2_restore_neutron;

/* Setting parameters for component 'guide_in' [17]. */
MCNUM mccguide_in_xmin;
MCNUM mccguide_in_xmax;
MCNUM mccguide_in_ymin;
MCNUM mccguide_in_ymax;
MCNUM mccguide_in_radius;
MCNUM mccguide_in_xwidth;
MCNUM mccguide_in_yheight;

/* Setting parameters for component 'guide_straight1' [18]. */
MCNUM mccguide_straight1_w1;
MCNUM mccguide_straight1_h1;
MCNUM mccguide_straight1_w2;
MCNUM mccguide_straight1_h2;
MCNUM mccguide_straight1_l;
MCNUM mccguide_straight1_R0;
MCNUM mccguide_straight1_Qc;
MCNUM mccguide_straight1_alpha;
MCNUM mccguide_straight1_m;
MCNUM mccguide_straight1_W;
MCNUM mccguide_straight1_nslit;
MCNUM mccguide_straight1_d;
MCNUM mccguide_straight1_mleft;
MCNUM mccguide_straight1_mright;
MCNUM mccguide_straight1_mtop;
MCNUM mccguide_straight1_mbottom;
MCNUM mccguide_straight1_nhslit;
MCNUM mccguide_straight1_G;
MCNUM mccguide_straight1_aleft;
MCNUM mccguide_straight1_aright;
MCNUM mccguide_straight1_atop;
MCNUM mccguide_straight1_abottom;
MCNUM mccguide_straight1_wavy;
MCNUM mccguide_straight1_wavy_z;
MCNUM mccguide_straight1_wavy_tb;
MCNUM mccguide_straight1_wavy_lr;
MCNUM mccguide_straight1_chamfers;
MCNUM mccguide_straight1_chamfers_z;
MCNUM mccguide_straight1_chamfers_lr;
MCNUM mccguide_straight1_chamfers_tb;
MCNUM mccguide_straight1_nelements;
MCNUM mccguide_straight1_nu;
MCNUM mccguide_straight1_phase;
char mccguide_straight1_reflect[16384];

/* Setting parameters for component 'guide_straight2' [19]. */
MCNUM mccguide_straight2_w1;
MCNUM mccguide_straight2_h1;
MCNUM mccguide_straight2_w2;
MCNUM mccguide_straight2_h2;
MCNUM mccguide_straight2_l;
MCNUM mccguide_straight2_R0;
MCNUM mccguide_straight2_Qc;
MCNUM mccguide_straight2_alpha;
MCNUM mccguide_straight2_m;
MCNUM mccguide_straight2_W;
MCNUM mccguide_straight2_nslit;
MCNUM mccguide_straight2_d;
MCNUM mccguide_straight2_mleft;
MCNUM mccguide_straight2_mright;
MCNUM mccguide_straight2_mtop;
MCNUM mccguide_straight2_mbottom;
MCNUM mccguide_straight2_nhslit;
MCNUM mccguide_straight2_G;
MCNUM mccguide_straight2_aleft;
MCNUM mccguide_straight2_aright;
MCNUM mccguide_straight2_atop;
MCNUM mccguide_straight2_abottom;
MCNUM mccguide_straight2_wavy;
MCNUM mccguide_straight2_wavy_z;
MCNUM mccguide_straight2_wavy_tb;
MCNUM mccguide_straight2_wavy_lr;
MCNUM mccguide_straight2_chamfers;
MCNUM mccguide_straight2_chamfers_z;
MCNUM mccguide_straight2_chamfers_lr;
MCNUM mccguide_straight2_chamfers_tb;
MCNUM mccguide_straight2_nelements;
MCNUM mccguide_straight2_nu;
MCNUM mccguide_straight2_phase;
char mccguide_straight2_reflect[16384];

/* Setting parameters for component 'guide_straight3' [20]. */
MCNUM mccguide_straight3_w1;
MCNUM mccguide_straight3_h1;
MCNUM mccguide_straight3_w2;
MCNUM mccguide_straight3_h2;
MCNUM mccguide_straight3_l;
MCNUM mccguide_straight3_R0;
MCNUM mccguide_straight3_Qc;
MCNUM mccguide_straight3_alpha;
MCNUM mccguide_straight3_m;
MCNUM mccguide_straight3_W;
MCNUM mccguide_straight3_nslit;
MCNUM mccguide_straight3_d;
MCNUM mccguide_straight3_mleft;
MCNUM mccguide_straight3_mright;
MCNUM mccguide_straight3_mtop;
MCNUM mccguide_straight3_mbottom;
MCNUM mccguide_straight3_nhslit;
MCNUM mccguide_straight3_G;
MCNUM mccguide_straight3_aleft;
MCNUM mccguide_straight3_aright;
MCNUM mccguide_straight3_atop;
MCNUM mccguide_straight3_abottom;
MCNUM mccguide_straight3_wavy;
MCNUM mccguide_straight3_wavy_z;
MCNUM mccguide_straight3_wavy_tb;
MCNUM mccguide_straight3_wavy_lr;
MCNUM mccguide_straight3_chamfers;
MCNUM mccguide_straight3_chamfers_z;
MCNUM mccguide_straight3_chamfers_lr;
MCNUM mccguide_straight3_chamfers_tb;
MCNUM mccguide_straight3_nelements;
MCNUM mccguide_straight3_nu;
MCNUM mccguide_straight3_phase;
char mccguide_straight3_reflect[16384];

/* Setting parameters for component 'guide_straight4' [21]. */
MCNUM mccguide_straight4_w1;
MCNUM mccguide_straight4_h1;
MCNUM mccguide_straight4_w2;
MCNUM mccguide_straight4_h2;
MCNUM mccguide_straight4_l;
MCNUM mccguide_straight4_R0;
MCNUM mccguide_straight4_Qc;
MCNUM mccguide_straight4_alpha;
MCNUM mccguide_straight4_m;
MCNUM mccguide_straight4_W;
MCNUM mccguide_straight4_nslit;
MCNUM mccguide_straight4_d;
MCNUM mccguide_straight4_mleft;
MCNUM mccguide_straight4_mright;
MCNUM mccguide_straight4_mtop;
MCNUM mccguide_straight4_mbottom;
MCNUM mccguide_straight4_nhslit;
MCNUM mccguide_straight4_G;
MCNUM mccguide_straight4_aleft;
MCNUM mccguide_straight4_aright;
MCNUM mccguide_straight4_atop;
MCNUM mccguide_straight4_abottom;
MCNUM mccguide_straight4_wavy;
MCNUM mccguide_straight4_wavy_z;
MCNUM mccguide_straight4_wavy_tb;
MCNUM mccguide_straight4_wavy_lr;
MCNUM mccguide_straight4_chamfers;
MCNUM mccguide_straight4_chamfers_z;
MCNUM mccguide_straight4_chamfers_lr;
MCNUM mccguide_straight4_chamfers_tb;
MCNUM mccguide_straight4_nelements;
MCNUM mccguide_straight4_nu;
MCNUM mccguide_straight4_phase;
char mccguide_straight4_reflect[16384];

/* Setting parameters for component 'psd3' [22]. */
int mccpsd3_nx;
int mccpsd3_ny;
char mccpsd3_filename[16384];
MCNUM mccpsd3_xmin;
MCNUM mccpsd3_xmax;
MCNUM mccpsd3_ymin;
MCNUM mccpsd3_ymax;
MCNUM mccpsd3_xwidth;
MCNUM mccpsd3_yheight;
MCNUM mccpsd3_restore_neutron;

/* Setting parameters for component 'aperture1' [23]. */
MCNUM mccaperture1_xmin;
MCNUM mccaperture1_xmax;
MCNUM mccaperture1_ymin;
MCNUM mccaperture1_ymax;
MCNUM mccaperture1_radius;
MCNUM mccaperture1_xwidth;
MCNUM mccaperture1_yheight;

/* Definition parameters for component 'lmonitor2' [24]. */
#define mcclmonitor2_nL 140
/* Setting parameters for component 'lmonitor2' [24]. */
char mcclmonitor2_filename[16384];
MCNUM mcclmonitor2_xmin;
MCNUM mcclmonitor2_xmax;
MCNUM mcclmonitor2_ymin;
MCNUM mcclmonitor2_ymax;
MCNUM mcclmonitor2_xwidth;
MCNUM mcclmonitor2_yheight;
MCNUM mcclmonitor2_Lmin;
MCNUM mcclmonitor2_Lmax;
MCNUM mcclmonitor2_restore_neutron;
int mcclmonitor2_nowritefile;

/* Setting parameters for component 'S6' [25]. */
MCNUM mccS6_xmin;
MCNUM mccS6_xmax;
MCNUM mccS6_ymin;
MCNUM mccS6_ymax;
MCNUM mccS6_radius;
MCNUM mccS6_xwidth;
MCNUM mccS6_yheight;

/* Setting parameters for component 'APERTURE2' [26]. */
MCNUM mccAPERTURE2_xmin;
MCNUM mccAPERTURE2_xmax;
MCNUM mccAPERTURE2_ymin;
MCNUM mccAPERTURE2_ymax;
MCNUM mccAPERTURE2_radius;
MCNUM mccAPERTURE2_xwidth;
MCNUM mccAPERTURE2_yheight;

/* Definition parameters for component 'lmon2' [27]. */
#define mcclmon2_nL 140
/* Setting parameters for component 'lmon2' [27]. */
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

/* Setting parameters for component 'psd4' [28]. */
int mccpsd4_nx;
int mccpsd4_ny;
char mccpsd4_filename[16384];
MCNUM mccpsd4_xmin;
MCNUM mccpsd4_xmax;
MCNUM mccpsd4_ymin;
MCNUM mccpsd4_ymax;
MCNUM mccpsd4_xwidth;
MCNUM mccpsd4_yheight;
MCNUM mccpsd4_restore_neutron;

/* User component declarations. */

/* User declarations for component 'Origin' [1]. */
#define mccompcurname  Origin
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
#line 815 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../contrib/ISIS_moderator.comp"
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
#line 8780 "./ISIS_SANS2d.c"
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

/* User declarations for component 'lmon1' [3]. */
#define mccompcurname  lmon1
#define mccompcurtype  L_monitor
#define mccompcurindex 3
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
#line 57 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/L_monitor.comp"
double L_N[nL];
double L_p[nL], L_p2[nL];
#line 8834 "./ISIS_SANS2d.c"
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

/* User declarations for component 'psd1' [4]. */
#define mccompcurname  psd1
#define mccompcurtype  PSD_monitor
#define mccompcurindex 4
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
#line 8875 "./ISIS_SANS2d.c"
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

/* User declarations for component 'bender1' [5]. */
#define mccompcurname  bender1
#define mccompcurtype  Guide_gravity
#define mccompcurindex 5
#define GVars mccbender1_GVars
#define pTable mccbender1_pTable
#define w1 mccbender1_w1
#define h1 mccbender1_h1
#define w2 mccbender1_w2
#define h2 mccbender1_h2
#define l mccbender1_l
#define R0 mccbender1_R0
#define Qc mccbender1_Qc
#define alpha mccbender1_alpha
#define m mccbender1_m
#define W mccbender1_W
#define nslit mccbender1_nslit
#define d mccbender1_d
#define mleft mccbender1_mleft
#define mright mccbender1_mright
#define mtop mccbender1_mtop
#define mbottom mccbender1_mbottom
#define nhslit mccbender1_nhslit
#define G mccbender1_G
#define aleft mccbender1_aleft
#define aright mccbender1_aright
#define atop mccbender1_atop
#define abottom mccbender1_abottom
#define wavy mccbender1_wavy
#define wavy_z mccbender1_wavy_z
#define wavy_tb mccbender1_wavy_tb
#define wavy_lr mccbender1_wavy_lr
#define chamfers mccbender1_chamfers
#define chamfers_z mccbender1_chamfers_z
#define chamfers_lr mccbender1_chamfers_lr
#define chamfers_tb mccbender1_chamfers_tb
#define nelements mccbender1_nelements
#define nu mccbender1_nu
#define phase mccbender1_phase
#define reflect mccbender1_reflect
#line 334 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide_gravity.comp"
  Gravity_guide_Vars_type GVars;
  t_Table pTable;
#line 8936 "./ISIS_SANS2d.c"
#undef reflect
#undef phase
#undef nu
#undef nelements
#undef chamfers_tb
#undef chamfers_lr
#undef chamfers_z
#undef chamfers
#undef wavy_lr
#undef wavy_tb
#undef wavy_z
#undef wavy
#undef abottom
#undef atop
#undef aright
#undef aleft
#undef G
#undef nhslit
#undef mbottom
#undef mtop
#undef mright
#undef mleft
#undef d
#undef nslit
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
#undef pTable
#undef GVars
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'bender2' [6]. */
#define mccompcurname  bender2
#define mccompcurtype  Guide_gravity
#define mccompcurindex 6
#define GVars mccbender2_GVars
#define pTable mccbender2_pTable
#define w1 mccbender2_w1
#define h1 mccbender2_h1
#define w2 mccbender2_w2
#define h2 mccbender2_h2
#define l mccbender2_l
#define R0 mccbender2_R0
#define Qc mccbender2_Qc
#define alpha mccbender2_alpha
#define m mccbender2_m
#define W mccbender2_W
#define nslit mccbender2_nslit
#define d mccbender2_d
#define mleft mccbender2_mleft
#define mright mccbender2_mright
#define mtop mccbender2_mtop
#define mbottom mccbender2_mbottom
#define nhslit mccbender2_nhslit
#define G mccbender2_G
#define aleft mccbender2_aleft
#define aright mccbender2_aright
#define atop mccbender2_atop
#define abottom mccbender2_abottom
#define wavy mccbender2_wavy
#define wavy_z mccbender2_wavy_z
#define wavy_tb mccbender2_wavy_tb
#define wavy_lr mccbender2_wavy_lr
#define chamfers mccbender2_chamfers
#define chamfers_z mccbender2_chamfers_z
#define chamfers_lr mccbender2_chamfers_lr
#define chamfers_tb mccbender2_chamfers_tb
#define nelements mccbender2_nelements
#define nu mccbender2_nu
#define phase mccbender2_phase
#define reflect mccbender2_reflect
#line 334 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide_gravity.comp"
  Gravity_guide_Vars_type GVars;
  t_Table pTable;
#line 9020 "./ISIS_SANS2d.c"
#undef reflect
#undef phase
#undef nu
#undef nelements
#undef chamfers_tb
#undef chamfers_lr
#undef chamfers_z
#undef chamfers
#undef wavy_lr
#undef wavy_tb
#undef wavy_z
#undef wavy
#undef abottom
#undef atop
#undef aright
#undef aleft
#undef G
#undef nhslit
#undef mbottom
#undef mtop
#undef mright
#undef mleft
#undef d
#undef nslit
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
#undef pTable
#undef GVars
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'bender3' [7]. */
#define mccompcurname  bender3
#define mccompcurtype  Guide_gravity
#define mccompcurindex 7
#define GVars mccbender3_GVars
#define pTable mccbender3_pTable
#define w1 mccbender3_w1
#define h1 mccbender3_h1
#define w2 mccbender3_w2
#define h2 mccbender3_h2
#define l mccbender3_l
#define R0 mccbender3_R0
#define Qc mccbender3_Qc
#define alpha mccbender3_alpha
#define m mccbender3_m
#define W mccbender3_W
#define nslit mccbender3_nslit
#define d mccbender3_d
#define mleft mccbender3_mleft
#define mright mccbender3_mright
#define mtop mccbender3_mtop
#define mbottom mccbender3_mbottom
#define nhslit mccbender3_nhslit
#define G mccbender3_G
#define aleft mccbender3_aleft
#define aright mccbender3_aright
#define atop mccbender3_atop
#define abottom mccbender3_abottom
#define wavy mccbender3_wavy
#define wavy_z mccbender3_wavy_z
#define wavy_tb mccbender3_wavy_tb
#define wavy_lr mccbender3_wavy_lr
#define chamfers mccbender3_chamfers
#define chamfers_z mccbender3_chamfers_z
#define chamfers_lr mccbender3_chamfers_lr
#define chamfers_tb mccbender3_chamfers_tb
#define nelements mccbender3_nelements
#define nu mccbender3_nu
#define phase mccbender3_phase
#define reflect mccbender3_reflect
#line 334 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide_gravity.comp"
  Gravity_guide_Vars_type GVars;
  t_Table pTable;
#line 9104 "./ISIS_SANS2d.c"
#undef reflect
#undef phase
#undef nu
#undef nelements
#undef chamfers_tb
#undef chamfers_lr
#undef chamfers_z
#undef chamfers
#undef wavy_lr
#undef wavy_tb
#undef wavy_z
#undef wavy
#undef abottom
#undef atop
#undef aright
#undef aleft
#undef G
#undef nhslit
#undef mbottom
#undef mtop
#undef mright
#undef mleft
#undef d
#undef nslit
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
#undef pTable
#undef GVars
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'bender4' [8]. */
#define mccompcurname  bender4
#define mccompcurtype  Guide_gravity
#define mccompcurindex 8
#define GVars mccbender4_GVars
#define pTable mccbender4_pTable
#define w1 mccbender4_w1
#define h1 mccbender4_h1
#define w2 mccbender4_w2
#define h2 mccbender4_h2
#define l mccbender4_l
#define R0 mccbender4_R0
#define Qc mccbender4_Qc
#define alpha mccbender4_alpha
#define m mccbender4_m
#define W mccbender4_W
#define nslit mccbender4_nslit
#define d mccbender4_d
#define mleft mccbender4_mleft
#define mright mccbender4_mright
#define mtop mccbender4_mtop
#define mbottom mccbender4_mbottom
#define nhslit mccbender4_nhslit
#define G mccbender4_G
#define aleft mccbender4_aleft
#define aright mccbender4_aright
#define atop mccbender4_atop
#define abottom mccbender4_abottom
#define wavy mccbender4_wavy
#define wavy_z mccbender4_wavy_z
#define wavy_tb mccbender4_wavy_tb
#define wavy_lr mccbender4_wavy_lr
#define chamfers mccbender4_chamfers
#define chamfers_z mccbender4_chamfers_z
#define chamfers_lr mccbender4_chamfers_lr
#define chamfers_tb mccbender4_chamfers_tb
#define nelements mccbender4_nelements
#define nu mccbender4_nu
#define phase mccbender4_phase
#define reflect mccbender4_reflect
#line 334 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide_gravity.comp"
  Gravity_guide_Vars_type GVars;
  t_Table pTable;
#line 9188 "./ISIS_SANS2d.c"
#undef reflect
#undef phase
#undef nu
#undef nelements
#undef chamfers_tb
#undef chamfers_lr
#undef chamfers_z
#undef chamfers
#undef wavy_lr
#undef wavy_tb
#undef wavy_z
#undef wavy
#undef abottom
#undef atop
#undef aright
#undef aleft
#undef G
#undef nhslit
#undef mbottom
#undef mtop
#undef mright
#undef mleft
#undef d
#undef nslit
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
#undef pTable
#undef GVars
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'bender5' [9]. */
#define mccompcurname  bender5
#define mccompcurtype  Guide_gravity
#define mccompcurindex 9
#define GVars mccbender5_GVars
#define pTable mccbender5_pTable
#define w1 mccbender5_w1
#define h1 mccbender5_h1
#define w2 mccbender5_w2
#define h2 mccbender5_h2
#define l mccbender5_l
#define R0 mccbender5_R0
#define Qc mccbender5_Qc
#define alpha mccbender5_alpha
#define m mccbender5_m
#define W mccbender5_W
#define nslit mccbender5_nslit
#define d mccbender5_d
#define mleft mccbender5_mleft
#define mright mccbender5_mright
#define mtop mccbender5_mtop
#define mbottom mccbender5_mbottom
#define nhslit mccbender5_nhslit
#define G mccbender5_G
#define aleft mccbender5_aleft
#define aright mccbender5_aright
#define atop mccbender5_atop
#define abottom mccbender5_abottom
#define wavy mccbender5_wavy
#define wavy_z mccbender5_wavy_z
#define wavy_tb mccbender5_wavy_tb
#define wavy_lr mccbender5_wavy_lr
#define chamfers mccbender5_chamfers
#define chamfers_z mccbender5_chamfers_z
#define chamfers_lr mccbender5_chamfers_lr
#define chamfers_tb mccbender5_chamfers_tb
#define nelements mccbender5_nelements
#define nu mccbender5_nu
#define phase mccbender5_phase
#define reflect mccbender5_reflect
#line 334 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide_gravity.comp"
  Gravity_guide_Vars_type GVars;
  t_Table pTable;
#line 9272 "./ISIS_SANS2d.c"
#undef reflect
#undef phase
#undef nu
#undef nelements
#undef chamfers_tb
#undef chamfers_lr
#undef chamfers_z
#undef chamfers
#undef wavy_lr
#undef wavy_tb
#undef wavy_z
#undef wavy
#undef abottom
#undef atop
#undef aright
#undef aleft
#undef G
#undef nhslit
#undef mbottom
#undef mtop
#undef mright
#undef mleft
#undef d
#undef nslit
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
#undef pTable
#undef GVars
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'bender6' [10]. */
#define mccompcurname  bender6
#define mccompcurtype  Guide_gravity
#define mccompcurindex 10
#define GVars mccbender6_GVars
#define pTable mccbender6_pTable
#define w1 mccbender6_w1
#define h1 mccbender6_h1
#define w2 mccbender6_w2
#define h2 mccbender6_h2
#define l mccbender6_l
#define R0 mccbender6_R0
#define Qc mccbender6_Qc
#define alpha mccbender6_alpha
#define m mccbender6_m
#define W mccbender6_W
#define nslit mccbender6_nslit
#define d mccbender6_d
#define mleft mccbender6_mleft
#define mright mccbender6_mright
#define mtop mccbender6_mtop
#define mbottom mccbender6_mbottom
#define nhslit mccbender6_nhslit
#define G mccbender6_G
#define aleft mccbender6_aleft
#define aright mccbender6_aright
#define atop mccbender6_atop
#define abottom mccbender6_abottom
#define wavy mccbender6_wavy
#define wavy_z mccbender6_wavy_z
#define wavy_tb mccbender6_wavy_tb
#define wavy_lr mccbender6_wavy_lr
#define chamfers mccbender6_chamfers
#define chamfers_z mccbender6_chamfers_z
#define chamfers_lr mccbender6_chamfers_lr
#define chamfers_tb mccbender6_chamfers_tb
#define nelements mccbender6_nelements
#define nu mccbender6_nu
#define phase mccbender6_phase
#define reflect mccbender6_reflect
#line 334 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide_gravity.comp"
  Gravity_guide_Vars_type GVars;
  t_Table pTable;
#line 9356 "./ISIS_SANS2d.c"
#undef reflect
#undef phase
#undef nu
#undef nelements
#undef chamfers_tb
#undef chamfers_lr
#undef chamfers_z
#undef chamfers
#undef wavy_lr
#undef wavy_tb
#undef wavy_z
#undef wavy
#undef abottom
#undef atop
#undef aright
#undef aleft
#undef G
#undef nhslit
#undef mbottom
#undef mtop
#undef mright
#undef mleft
#undef d
#undef nslit
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
#undef pTable
#undef GVars
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'bender7' [11]. */
#define mccompcurname  bender7
#define mccompcurtype  Guide_gravity
#define mccompcurindex 11
#define GVars mccbender7_GVars
#define pTable mccbender7_pTable
#define w1 mccbender7_w1
#define h1 mccbender7_h1
#define w2 mccbender7_w2
#define h2 mccbender7_h2
#define l mccbender7_l
#define R0 mccbender7_R0
#define Qc mccbender7_Qc
#define alpha mccbender7_alpha
#define m mccbender7_m
#define W mccbender7_W
#define nslit mccbender7_nslit
#define d mccbender7_d
#define mleft mccbender7_mleft
#define mright mccbender7_mright
#define mtop mccbender7_mtop
#define mbottom mccbender7_mbottom
#define nhslit mccbender7_nhslit
#define G mccbender7_G
#define aleft mccbender7_aleft
#define aright mccbender7_aright
#define atop mccbender7_atop
#define abottom mccbender7_abottom
#define wavy mccbender7_wavy
#define wavy_z mccbender7_wavy_z
#define wavy_tb mccbender7_wavy_tb
#define wavy_lr mccbender7_wavy_lr
#define chamfers mccbender7_chamfers
#define chamfers_z mccbender7_chamfers_z
#define chamfers_lr mccbender7_chamfers_lr
#define chamfers_tb mccbender7_chamfers_tb
#define nelements mccbender7_nelements
#define nu mccbender7_nu
#define phase mccbender7_phase
#define reflect mccbender7_reflect
#line 334 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide_gravity.comp"
  Gravity_guide_Vars_type GVars;
  t_Table pTable;
#line 9440 "./ISIS_SANS2d.c"
#undef reflect
#undef phase
#undef nu
#undef nelements
#undef chamfers_tb
#undef chamfers_lr
#undef chamfers_z
#undef chamfers
#undef wavy_lr
#undef wavy_tb
#undef wavy_z
#undef wavy
#undef abottom
#undef atop
#undef aright
#undef aleft
#undef G
#undef nhslit
#undef mbottom
#undef mtop
#undef mright
#undef mleft
#undef d
#undef nslit
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
#undef pTable
#undef GVars
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'bender8' [12]. */
#define mccompcurname  bender8
#define mccompcurtype  Guide_gravity
#define mccompcurindex 12
#define GVars mccbender8_GVars
#define pTable mccbender8_pTable
#define w1 mccbender8_w1
#define h1 mccbender8_h1
#define w2 mccbender8_w2
#define h2 mccbender8_h2
#define l mccbender8_l
#define R0 mccbender8_R0
#define Qc mccbender8_Qc
#define alpha mccbender8_alpha
#define m mccbender8_m
#define W mccbender8_W
#define nslit mccbender8_nslit
#define d mccbender8_d
#define mleft mccbender8_mleft
#define mright mccbender8_mright
#define mtop mccbender8_mtop
#define mbottom mccbender8_mbottom
#define nhslit mccbender8_nhslit
#define G mccbender8_G
#define aleft mccbender8_aleft
#define aright mccbender8_aright
#define atop mccbender8_atop
#define abottom mccbender8_abottom
#define wavy mccbender8_wavy
#define wavy_z mccbender8_wavy_z
#define wavy_tb mccbender8_wavy_tb
#define wavy_lr mccbender8_wavy_lr
#define chamfers mccbender8_chamfers
#define chamfers_z mccbender8_chamfers_z
#define chamfers_lr mccbender8_chamfers_lr
#define chamfers_tb mccbender8_chamfers_tb
#define nelements mccbender8_nelements
#define nu mccbender8_nu
#define phase mccbender8_phase
#define reflect mccbender8_reflect
#line 334 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide_gravity.comp"
  Gravity_guide_Vars_type GVars;
  t_Table pTable;
#line 9524 "./ISIS_SANS2d.c"
#undef reflect
#undef phase
#undef nu
#undef nelements
#undef chamfers_tb
#undef chamfers_lr
#undef chamfers_z
#undef chamfers
#undef wavy_lr
#undef wavy_tb
#undef wavy_z
#undef wavy
#undef abottom
#undef atop
#undef aright
#undef aleft
#undef G
#undef nhslit
#undef mbottom
#undef mtop
#undef mright
#undef mleft
#undef d
#undef nslit
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
#undef pTable
#undef GVars
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'bender9' [13]. */
#define mccompcurname  bender9
#define mccompcurtype  Guide_gravity
#define mccompcurindex 13
#define GVars mccbender9_GVars
#define pTable mccbender9_pTable
#define w1 mccbender9_w1
#define h1 mccbender9_h1
#define w2 mccbender9_w2
#define h2 mccbender9_h2
#define l mccbender9_l
#define R0 mccbender9_R0
#define Qc mccbender9_Qc
#define alpha mccbender9_alpha
#define m mccbender9_m
#define W mccbender9_W
#define nslit mccbender9_nslit
#define d mccbender9_d
#define mleft mccbender9_mleft
#define mright mccbender9_mright
#define mtop mccbender9_mtop
#define mbottom mccbender9_mbottom
#define nhslit mccbender9_nhslit
#define G mccbender9_G
#define aleft mccbender9_aleft
#define aright mccbender9_aright
#define atop mccbender9_atop
#define abottom mccbender9_abottom
#define wavy mccbender9_wavy
#define wavy_z mccbender9_wavy_z
#define wavy_tb mccbender9_wavy_tb
#define wavy_lr mccbender9_wavy_lr
#define chamfers mccbender9_chamfers
#define chamfers_z mccbender9_chamfers_z
#define chamfers_lr mccbender9_chamfers_lr
#define chamfers_tb mccbender9_chamfers_tb
#define nelements mccbender9_nelements
#define nu mccbender9_nu
#define phase mccbender9_phase
#define reflect mccbender9_reflect
#line 334 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide_gravity.comp"
  Gravity_guide_Vars_type GVars;
  t_Table pTable;
#line 9608 "./ISIS_SANS2d.c"
#undef reflect
#undef phase
#undef nu
#undef nelements
#undef chamfers_tb
#undef chamfers_lr
#undef chamfers_z
#undef chamfers
#undef wavy_lr
#undef wavy_tb
#undef wavy_z
#undef wavy
#undef abottom
#undef atop
#undef aright
#undef aleft
#undef G
#undef nhslit
#undef mbottom
#undef mtop
#undef mright
#undef mleft
#undef d
#undef nslit
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
#undef pTable
#undef GVars
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'bender10' [14]. */
#define mccompcurname  bender10
#define mccompcurtype  Guide_gravity
#define mccompcurindex 14
#define GVars mccbender10_GVars
#define pTable mccbender10_pTable
#define w1 mccbender10_w1
#define h1 mccbender10_h1
#define w2 mccbender10_w2
#define h2 mccbender10_h2
#define l mccbender10_l
#define R0 mccbender10_R0
#define Qc mccbender10_Qc
#define alpha mccbender10_alpha
#define m mccbender10_m
#define W mccbender10_W
#define nslit mccbender10_nslit
#define d mccbender10_d
#define mleft mccbender10_mleft
#define mright mccbender10_mright
#define mtop mccbender10_mtop
#define mbottom mccbender10_mbottom
#define nhslit mccbender10_nhslit
#define G mccbender10_G
#define aleft mccbender10_aleft
#define aright mccbender10_aright
#define atop mccbender10_atop
#define abottom mccbender10_abottom
#define wavy mccbender10_wavy
#define wavy_z mccbender10_wavy_z
#define wavy_tb mccbender10_wavy_tb
#define wavy_lr mccbender10_wavy_lr
#define chamfers mccbender10_chamfers
#define chamfers_z mccbender10_chamfers_z
#define chamfers_lr mccbender10_chamfers_lr
#define chamfers_tb mccbender10_chamfers_tb
#define nelements mccbender10_nelements
#define nu mccbender10_nu
#define phase mccbender10_phase
#define reflect mccbender10_reflect
#line 334 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide_gravity.comp"
  Gravity_guide_Vars_type GVars;
  t_Table pTable;
#line 9692 "./ISIS_SANS2d.c"
#undef reflect
#undef phase
#undef nu
#undef nelements
#undef chamfers_tb
#undef chamfers_lr
#undef chamfers_z
#undef chamfers
#undef wavy_lr
#undef wavy_tb
#undef wavy_z
#undef wavy
#undef abottom
#undef atop
#undef aright
#undef aleft
#undef G
#undef nhslit
#undef mbottom
#undef mtop
#undef mright
#undef mleft
#undef d
#undef nslit
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
#undef pTable
#undef GVars
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'lmonb' [15]. */
#define mccompcurname  lmonb
#define mccompcurtype  L_monitor
#define mccompcurindex 15
#define nL mcclmonb_nL
#define L_N mcclmonb_L_N
#define L_p mcclmonb_L_p
#define L_p2 mcclmonb_L_p2
#define filename mcclmonb_filename
#define xmin mcclmonb_xmin
#define xmax mcclmonb_xmax
#define ymin mcclmonb_ymin
#define ymax mcclmonb_ymax
#define xwidth mcclmonb_xwidth
#define yheight mcclmonb_yheight
#define Lmin mcclmonb_Lmin
#define Lmax mcclmonb_Lmax
#define restore_neutron mcclmonb_restore_neutron
#define nowritefile mcclmonb_nowritefile
#line 57 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/L_monitor.comp"
double L_N[nL];
double L_p[nL], L_p2[nL];
#line 9755 "./ISIS_SANS2d.c"
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

/* User declarations for component 'psd2' [16]. */
#define mccompcurname  psd2
#define mccompcurtype  PSD_monitor
#define mccompcurindex 16
#define PSD_N mccpsd2_PSD_N
#define PSD_p mccpsd2_PSD_p
#define PSD_p2 mccpsd2_PSD_p2
#define nx mccpsd2_nx
#define ny mccpsd2_ny
#define filename mccpsd2_filename
#define xmin mccpsd2_xmin
#define xmax mccpsd2_xmax
#define ymin mccpsd2_ymin
#define ymax mccpsd2_ymax
#define xwidth mccpsd2_xwidth
#define yheight mccpsd2_yheight
#define restore_neutron mccpsd2_restore_neutron
#line 62 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
  DArray2d PSD_N;
  DArray2d PSD_p;
  DArray2d PSD_p2;
#line 9796 "./ISIS_SANS2d.c"
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

/* User declarations for component 'guide_in' [17]. */
#define mccompcurname  guide_in
#define mccompcurtype  Slit
#define mccompcurindex 17
#define xmin mccguide_in_xmin
#define xmax mccguide_in_xmax
#define ymin mccguide_in_ymin
#define ymax mccguide_in_ymax
#define radius mccguide_in_radius
#define xwidth mccguide_in_xwidth
#define yheight mccguide_in_yheight
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

/* User declarations for component 'guide_straight1' [18]. */
#define mccompcurname  guide_straight1
#define mccompcurtype  Guide_gravity
#define mccompcurindex 18
#define GVars mccguide_straight1_GVars
#define pTable mccguide_straight1_pTable
#define w1 mccguide_straight1_w1
#define h1 mccguide_straight1_h1
#define w2 mccguide_straight1_w2
#define h2 mccguide_straight1_h2
#define l mccguide_straight1_l
#define R0 mccguide_straight1_R0
#define Qc mccguide_straight1_Qc
#define alpha mccguide_straight1_alpha
#define m mccguide_straight1_m
#define W mccguide_straight1_W
#define nslit mccguide_straight1_nslit
#define d mccguide_straight1_d
#define mleft mccguide_straight1_mleft
#define mright mccguide_straight1_mright
#define mtop mccguide_straight1_mtop
#define mbottom mccguide_straight1_mbottom
#define nhslit mccguide_straight1_nhslit
#define G mccguide_straight1_G
#define aleft mccguide_straight1_aleft
#define aright mccguide_straight1_aright
#define atop mccguide_straight1_atop
#define abottom mccguide_straight1_abottom
#define wavy mccguide_straight1_wavy
#define wavy_z mccguide_straight1_wavy_z
#define wavy_tb mccguide_straight1_wavy_tb
#define wavy_lr mccguide_straight1_wavy_lr
#define chamfers mccguide_straight1_chamfers
#define chamfers_z mccguide_straight1_chamfers_z
#define chamfers_lr mccguide_straight1_chamfers_lr
#define chamfers_tb mccguide_straight1_chamfers_tb
#define nelements mccguide_straight1_nelements
#define nu mccguide_straight1_nu
#define phase mccguide_straight1_phase
#define reflect mccguide_straight1_reflect
#line 334 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide_gravity.comp"
  Gravity_guide_Vars_type GVars;
  t_Table pTable;
#line 9879 "./ISIS_SANS2d.c"
#undef reflect
#undef phase
#undef nu
#undef nelements
#undef chamfers_tb
#undef chamfers_lr
#undef chamfers_z
#undef chamfers
#undef wavy_lr
#undef wavy_tb
#undef wavy_z
#undef wavy
#undef abottom
#undef atop
#undef aright
#undef aleft
#undef G
#undef nhslit
#undef mbottom
#undef mtop
#undef mright
#undef mleft
#undef d
#undef nslit
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
#undef pTable
#undef GVars
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'guide_straight2' [19]. */
#define mccompcurname  guide_straight2
#define mccompcurtype  Guide_gravity
#define mccompcurindex 19
#define GVars mccguide_straight2_GVars
#define pTable mccguide_straight2_pTable
#define w1 mccguide_straight2_w1
#define h1 mccguide_straight2_h1
#define w2 mccguide_straight2_w2
#define h2 mccguide_straight2_h2
#define l mccguide_straight2_l
#define R0 mccguide_straight2_R0
#define Qc mccguide_straight2_Qc
#define alpha mccguide_straight2_alpha
#define m mccguide_straight2_m
#define W mccguide_straight2_W
#define nslit mccguide_straight2_nslit
#define d mccguide_straight2_d
#define mleft mccguide_straight2_mleft
#define mright mccguide_straight2_mright
#define mtop mccguide_straight2_mtop
#define mbottom mccguide_straight2_mbottom
#define nhslit mccguide_straight2_nhslit
#define G mccguide_straight2_G
#define aleft mccguide_straight2_aleft
#define aright mccguide_straight2_aright
#define atop mccguide_straight2_atop
#define abottom mccguide_straight2_abottom
#define wavy mccguide_straight2_wavy
#define wavy_z mccguide_straight2_wavy_z
#define wavy_tb mccguide_straight2_wavy_tb
#define wavy_lr mccguide_straight2_wavy_lr
#define chamfers mccguide_straight2_chamfers
#define chamfers_z mccguide_straight2_chamfers_z
#define chamfers_lr mccguide_straight2_chamfers_lr
#define chamfers_tb mccguide_straight2_chamfers_tb
#define nelements mccguide_straight2_nelements
#define nu mccguide_straight2_nu
#define phase mccguide_straight2_phase
#define reflect mccguide_straight2_reflect
#line 334 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide_gravity.comp"
  Gravity_guide_Vars_type GVars;
  t_Table pTable;
#line 9963 "./ISIS_SANS2d.c"
#undef reflect
#undef phase
#undef nu
#undef nelements
#undef chamfers_tb
#undef chamfers_lr
#undef chamfers_z
#undef chamfers
#undef wavy_lr
#undef wavy_tb
#undef wavy_z
#undef wavy
#undef abottom
#undef atop
#undef aright
#undef aleft
#undef G
#undef nhslit
#undef mbottom
#undef mtop
#undef mright
#undef mleft
#undef d
#undef nslit
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
#undef pTable
#undef GVars
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'guide_straight3' [20]. */
#define mccompcurname  guide_straight3
#define mccompcurtype  Guide_gravity
#define mccompcurindex 20
#define GVars mccguide_straight3_GVars
#define pTable mccguide_straight3_pTable
#define w1 mccguide_straight3_w1
#define h1 mccguide_straight3_h1
#define w2 mccguide_straight3_w2
#define h2 mccguide_straight3_h2
#define l mccguide_straight3_l
#define R0 mccguide_straight3_R0
#define Qc mccguide_straight3_Qc
#define alpha mccguide_straight3_alpha
#define m mccguide_straight3_m
#define W mccguide_straight3_W
#define nslit mccguide_straight3_nslit
#define d mccguide_straight3_d
#define mleft mccguide_straight3_mleft
#define mright mccguide_straight3_mright
#define mtop mccguide_straight3_mtop
#define mbottom mccguide_straight3_mbottom
#define nhslit mccguide_straight3_nhslit
#define G mccguide_straight3_G
#define aleft mccguide_straight3_aleft
#define aright mccguide_straight3_aright
#define atop mccguide_straight3_atop
#define abottom mccguide_straight3_abottom
#define wavy mccguide_straight3_wavy
#define wavy_z mccguide_straight3_wavy_z
#define wavy_tb mccguide_straight3_wavy_tb
#define wavy_lr mccguide_straight3_wavy_lr
#define chamfers mccguide_straight3_chamfers
#define chamfers_z mccguide_straight3_chamfers_z
#define chamfers_lr mccguide_straight3_chamfers_lr
#define chamfers_tb mccguide_straight3_chamfers_tb
#define nelements mccguide_straight3_nelements
#define nu mccguide_straight3_nu
#define phase mccguide_straight3_phase
#define reflect mccguide_straight3_reflect
#line 334 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide_gravity.comp"
  Gravity_guide_Vars_type GVars;
  t_Table pTable;
#line 10047 "./ISIS_SANS2d.c"
#undef reflect
#undef phase
#undef nu
#undef nelements
#undef chamfers_tb
#undef chamfers_lr
#undef chamfers_z
#undef chamfers
#undef wavy_lr
#undef wavy_tb
#undef wavy_z
#undef wavy
#undef abottom
#undef atop
#undef aright
#undef aleft
#undef G
#undef nhslit
#undef mbottom
#undef mtop
#undef mright
#undef mleft
#undef d
#undef nslit
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
#undef pTable
#undef GVars
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'guide_straight4' [21]. */
#define mccompcurname  guide_straight4
#define mccompcurtype  Guide_gravity
#define mccompcurindex 21
#define GVars mccguide_straight4_GVars
#define pTable mccguide_straight4_pTable
#define w1 mccguide_straight4_w1
#define h1 mccguide_straight4_h1
#define w2 mccguide_straight4_w2
#define h2 mccguide_straight4_h2
#define l mccguide_straight4_l
#define R0 mccguide_straight4_R0
#define Qc mccguide_straight4_Qc
#define alpha mccguide_straight4_alpha
#define m mccguide_straight4_m
#define W mccguide_straight4_W
#define nslit mccguide_straight4_nslit
#define d mccguide_straight4_d
#define mleft mccguide_straight4_mleft
#define mright mccguide_straight4_mright
#define mtop mccguide_straight4_mtop
#define mbottom mccguide_straight4_mbottom
#define nhslit mccguide_straight4_nhslit
#define G mccguide_straight4_G
#define aleft mccguide_straight4_aleft
#define aright mccguide_straight4_aright
#define atop mccguide_straight4_atop
#define abottom mccguide_straight4_abottom
#define wavy mccguide_straight4_wavy
#define wavy_z mccguide_straight4_wavy_z
#define wavy_tb mccguide_straight4_wavy_tb
#define wavy_lr mccguide_straight4_wavy_lr
#define chamfers mccguide_straight4_chamfers
#define chamfers_z mccguide_straight4_chamfers_z
#define chamfers_lr mccguide_straight4_chamfers_lr
#define chamfers_tb mccguide_straight4_chamfers_tb
#define nelements mccguide_straight4_nelements
#define nu mccguide_straight4_nu
#define phase mccguide_straight4_phase
#define reflect mccguide_straight4_reflect
#line 334 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide_gravity.comp"
  Gravity_guide_Vars_type GVars;
  t_Table pTable;
#line 10131 "./ISIS_SANS2d.c"
#undef reflect
#undef phase
#undef nu
#undef nelements
#undef chamfers_tb
#undef chamfers_lr
#undef chamfers_z
#undef chamfers
#undef wavy_lr
#undef wavy_tb
#undef wavy_z
#undef wavy
#undef abottom
#undef atop
#undef aright
#undef aleft
#undef G
#undef nhslit
#undef mbottom
#undef mtop
#undef mright
#undef mleft
#undef d
#undef nslit
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
#undef pTable
#undef GVars
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'psd3' [22]. */
#define mccompcurname  psd3
#define mccompcurtype  PSD_monitor
#define mccompcurindex 22
#define PSD_N mccpsd3_PSD_N
#define PSD_p mccpsd3_PSD_p
#define PSD_p2 mccpsd3_PSD_p2
#define nx mccpsd3_nx
#define ny mccpsd3_ny
#define filename mccpsd3_filename
#define xmin mccpsd3_xmin
#define xmax mccpsd3_xmax
#define ymin mccpsd3_ymin
#define ymax mccpsd3_ymax
#define xwidth mccpsd3_xwidth
#define yheight mccpsd3_yheight
#define restore_neutron mccpsd3_restore_neutron
#line 62 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
  DArray2d PSD_N;
  DArray2d PSD_p;
  DArray2d PSD_p2;
#line 10193 "./ISIS_SANS2d.c"
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

/* User declarations for component 'aperture1' [23]. */
#define mccompcurname  aperture1
#define mccompcurtype  Slit
#define mccompcurindex 23
#define xmin mccaperture1_xmin
#define xmax mccaperture1_xmax
#define ymin mccaperture1_ymin
#define ymax mccaperture1_ymax
#define radius mccaperture1_radius
#define xwidth mccaperture1_xwidth
#define yheight mccaperture1_yheight
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

/* User declarations for component 'lmonitor2' [24]. */
#define mccompcurname  lmonitor2
#define mccompcurtype  L_monitor
#define mccompcurindex 24
#define nL mcclmonitor2_nL
#define L_N mcclmonitor2_L_N
#define L_p mcclmonitor2_L_p
#define L_p2 mcclmonitor2_L_p2
#define filename mcclmonitor2_filename
#define xmin mcclmonitor2_xmin
#define xmax mcclmonitor2_xmax
#define ymin mcclmonitor2_ymin
#define ymax mcclmonitor2_ymax
#define xwidth mcclmonitor2_xwidth
#define yheight mcclmonitor2_yheight
#define Lmin mcclmonitor2_Lmin
#define Lmax mcclmonitor2_Lmax
#define restore_neutron mcclmonitor2_restore_neutron
#define nowritefile mcclmonitor2_nowritefile
#line 57 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/L_monitor.comp"
double L_N[nL];
double L_p[nL], L_p2[nL];
#line 10255 "./ISIS_SANS2d.c"
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

/* User declarations for component 'S6' [25]. */
#define mccompcurname  S6
#define mccompcurtype  Slit
#define mccompcurindex 25
#define xmin mccS6_xmin
#define xmax mccS6_xmax
#define ymin mccS6_ymin
#define ymax mccS6_ymax
#define radius mccS6_radius
#define xwidth mccS6_xwidth
#define yheight mccS6_yheight
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

/* User declarations for component 'APERTURE2' [26]. */
#define mccompcurname  APERTURE2
#define mccompcurtype  Slit
#define mccompcurindex 26
#define xmin mccAPERTURE2_xmin
#define xmax mccAPERTURE2_xmax
#define ymin mccAPERTURE2_ymin
#define ymax mccAPERTURE2_ymax
#define radius mccAPERTURE2_radius
#define xwidth mccAPERTURE2_xwidth
#define yheight mccAPERTURE2_yheight
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

/* User declarations for component 'lmon2' [27]. */
#define mccompcurname  lmon2
#define mccompcurtype  L_monitor
#define mccompcurindex 27
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
#line 57 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/L_monitor.comp"
double L_N[nL];
double L_p[nL], L_p2[nL];
#line 10341 "./ISIS_SANS2d.c"
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

/* User declarations for component 'psd4' [28]. */
#define mccompcurname  psd4
#define mccompcurtype  PSD_monitor
#define mccompcurindex 28
#define PSD_N mccpsd4_PSD_N
#define PSD_p mccpsd4_PSD_p
#define PSD_p2 mccpsd4_PSD_p2
#define nx mccpsd4_nx
#define ny mccpsd4_ny
#define filename mccpsd4_filename
#define xmin mccpsd4_xmin
#define xmax mccpsd4_xmax
#define ymin mccpsd4_ymin
#define ymax mccpsd4_ymax
#define xwidth mccpsd4_xwidth
#define yheight mccpsd4_yheight
#define restore_neutron mccpsd4_restore_neutron
#line 62 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
  DArray2d PSD_N;
  DArray2d PSD_p;
  DArray2d PSD_p2;
#line 10382 "./ISIS_SANS2d.c"
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
Coords mcposaisis_source, mcposrisis_source;
Rotation mcrotaisis_source, mcrotrisis_source;
Coords mcposalmon1, mcposrlmon1;
Rotation mcrotalmon1, mcrotrlmon1;
Coords mcposapsd1, mcposrpsd1;
Rotation mcrotapsd1, mcrotrpsd1;
Coords mcposabender1, mcposrbender1;
Rotation mcrotabender1, mcrotrbender1;
Coords mcposabender2, mcposrbender2;
Rotation mcrotabender2, mcrotrbender2;
Coords mcposabender3, mcposrbender3;
Rotation mcrotabender3, mcrotrbender3;
Coords mcposabender4, mcposrbender4;
Rotation mcrotabender4, mcrotrbender4;
Coords mcposabender5, mcposrbender5;
Rotation mcrotabender5, mcrotrbender5;
Coords mcposabender6, mcposrbender6;
Rotation mcrotabender6, mcrotrbender6;
Coords mcposabender7, mcposrbender7;
Rotation mcrotabender7, mcrotrbender7;
Coords mcposabender8, mcposrbender8;
Rotation mcrotabender8, mcrotrbender8;
Coords mcposabender9, mcposrbender9;
Rotation mcrotabender9, mcrotrbender9;
Coords mcposabender10, mcposrbender10;
Rotation mcrotabender10, mcrotrbender10;
Coords mcposalmonb, mcposrlmonb;
Rotation mcrotalmonb, mcrotrlmonb;
Coords mcposapsd2, mcposrpsd2;
Rotation mcrotapsd2, mcrotrpsd2;
Coords mcposaguide_in, mcposrguide_in;
Rotation mcrotaguide_in, mcrotrguide_in;
Coords mcposaguide_straight1, mcposrguide_straight1;
Rotation mcrotaguide_straight1, mcrotrguide_straight1;
Coords mcposaguide_straight2, mcposrguide_straight2;
Rotation mcrotaguide_straight2, mcrotrguide_straight2;
Coords mcposaguide_straight3, mcposrguide_straight3;
Rotation mcrotaguide_straight3, mcrotrguide_straight3;
Coords mcposaguide_straight4, mcposrguide_straight4;
Rotation mcrotaguide_straight4, mcrotrguide_straight4;
Coords mcposapsd3, mcposrpsd3;
Rotation mcrotapsd3, mcrotrpsd3;
Coords mcposaaperture1, mcposraperture1;
Rotation mcrotaaperture1, mcrotraperture1;
Coords mcposalmonitor2, mcposrlmonitor2;
Rotation mcrotalmonitor2, mcrotrlmonitor2;
Coords mcposaS6, mcposrS6;
Rotation mcrotaS6, mcrotrS6;
Coords mcposaAPERTURE2, mcposrAPERTURE2;
Rotation mcrotaAPERTURE2, mcrotrAPERTURE2;
Coords mcposalmon2, mcposrlmon2;
Rotation mcrotalmon2, mcrotrlmon2;
Coords mcposapsd4, mcposrpsd4;
Rotation mcrotapsd4, mcrotrpsd4;

MCNUM mcnx, mcny, mcnz, mcnvx, mcnvy, mcnvz, mcnt, mcnsx, mcnsy, mcnsz, mcnp;

/* end declare */

void mcinit(void) {
#define mccompcurname  ISIS_SANS2d
#define mccompcurtype  INSTRUMENT
#define mccompcurindex 0
#define mcposaISIS_SANS2d coords_set(0,0,0)
#define L1 mcipL1
#define A1w mcipA1w
#define A1h mcipA1h
#define S6 mcipS6
#define A2 mcipA2
#define Lmin mcipLmin
#define Lmax mcipLmax
#undef Lmax
#undef Lmin
#undef A2
#undef S6
#undef A1h
#undef A1w
#undef L1
#undef mcposaISIS_SANS2d
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

  SIG_MESSAGE("Origin (Init:Place/Rotate)");
  rot_set_rotation(mcrotaOrigin,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 10503 "./ISIS_SANS2d.c"
  rot_copy(mcrotrOrigin, mcrotaOrigin);
  mcposaOrigin = coords_set(
#line 59 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    0,
#line 59 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    0,
#line 59 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    0);
#line 10512 "./ISIS_SANS2d.c"
  mctc1 = coords_neg(mcposaOrigin);
  mcposrOrigin = rot_apply(mcrotaOrigin, mctc1);
  mcDEBUG_COMPONENT("Origin", mcposaOrigin, mcrotaOrigin)
  mccomp_posa[1] = mcposaOrigin;
  mccomp_posr[1] = mcposrOrigin;
  mcNCounter[1]  = mcPCounter[1] = mcP2Counter[1] = 0;
  mcAbsorbProp[1]= 0;
    /* Component isis_source. */
  /* Setting parameters for component isis_source. */
  SIG_MESSAGE("isis_source (Init:SetPar)");
#line 63 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccisis_source_Emin = - mcipLmax;
#line 63 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccisis_source_Emax = - mcipLmin;
#line 63 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccisis_source_dist = 3.68;
#line 64 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccisis_source_focus_xw = 0.0365;
#line 64 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccisis_source_focus_yh = 0.021;
#line 64 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccisis_source_xwidth = -1;
#line 64 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccisis_source_yheight = -1;
#line 65 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccisis_source_CAngle = 0.0;
#line 65 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccisis_source_SAC = 1;
#line 58 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccisis_source_Lmin = 0;
#line 58 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccisis_source_Lmax = 0;
#line 58 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccisis_source_target_index = + 1;
#line 58 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccisis_source_verbose = 0;
#line 10549 "./ISIS_SANS2d.c"

  SIG_MESSAGE("isis_source (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 67 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    (0.0)*DEG2RAD,
#line 67 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    (0.0)*DEG2RAD,
#line 67 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    (0.0)*DEG2RAD);
#line 10559 "./ISIS_SANS2d.c"
  rot_mul(mctr1, mcrotaOrigin, mcrotaisis_source);
  rot_transpose(mcrotaOrigin, mctr1);
  rot_mul(mcrotaisis_source, mctr1, mcrotrisis_source);
  mctc1 = coords_set(
#line 66 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    0.0,
#line 66 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    0.0,
#line 66 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    0.00001);
#line 10570 "./ISIS_SANS2d.c"
  rot_transpose(mcrotaOrigin, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaisis_source = coords_add(mcposaOrigin, mctc2);
  mctc1 = coords_sub(mcposaOrigin, mcposaisis_source);
  mcposrisis_source = rot_apply(mcrotaisis_source, mctc1);
  mcDEBUG_COMPONENT("isis_source", mcposaisis_source, mcrotaisis_source)
  mccomp_posa[2] = mcposaisis_source;
  mccomp_posr[2] = mcposrisis_source;
  mcNCounter[2]  = mcPCounter[2] = mcP2Counter[2] = 0;
  mcAbsorbProp[2]= 0;
    /* Component lmon1. */
  /* Setting parameters for component lmon1. */
  SIG_MESSAGE("lmon1 (Init:SetPar)");
#line 70 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  if("lmon1.dat") strncpy(mcclmon1_filename, "lmon1.dat" ? "lmon1.dat" : "", 16384); else mcclmon1_filename[0]='\0';
#line 70 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mcclmon1_xmin = -0.04;
#line 71 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mcclmon1_xmax = 0.04;
#line 71 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mcclmon1_ymin = -0.03;
#line 71 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mcclmon1_ymax = 0.03;
#line 51 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mcclmon1_xwidth = 0;
#line 51 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mcclmon1_yheight = 0;
#line 71 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mcclmon1_Lmin = 0.0;
#line 72 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mcclmon1_Lmax = 17.0;
#line 51 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mcclmon1_restore_neutron = 0;
#line 51 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mcclmon1_nowritefile = 0;
#line 10606 "./ISIS_SANS2d.c"

  SIG_MESSAGE("lmon1 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 10613 "./ISIS_SANS2d.c"
  rot_mul(mctr1, mcrotaisis_source, mcrotalmon1);
  rot_transpose(mcrotaisis_source, mctr1);
  rot_mul(mcrotalmon1, mctr1, mcrotrlmon1);
  mctc1 = coords_set(
#line 73 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    0.0,
#line 73 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    0.0,
#line 73 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    3.698);
#line 10624 "./ISIS_SANS2d.c"
  rot_transpose(mcrotaisis_source, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposalmon1 = coords_add(mcposaisis_source, mctc2);
  mctc1 = coords_sub(mcposaisis_source, mcposalmon1);
  mcposrlmon1 = rot_apply(mcrotalmon1, mctc1);
  mcDEBUG_COMPONENT("lmon1", mcposalmon1, mcrotalmon1)
  mccomp_posa[3] = mcposalmon1;
  mccomp_posr[3] = mcposrlmon1;
  mcNCounter[3]  = mcPCounter[3] = mcP2Counter[3] = 0;
  mcAbsorbProp[3]= 0;
    /* Component psd1. */
  /* Setting parameters for component psd1. */
  SIG_MESSAGE("psd1 (Init:SetPar)");
#line 76 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccpsd1_nx = 100;
#line 76 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccpsd1_ny = 100;
#line 76 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  if("psd1.dat") strncpy(mccpsd1_filename, "psd1.dat" ? "psd1.dat" : "", 16384); else mccpsd1_filename[0]='\0';
#line 76 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccpsd1_xmin = -0.05;
#line 77 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccpsd1_xmax = 0.05;
#line 77 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccpsd1_ymin = -0.05;
#line 77 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccpsd1_ymax = 0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccpsd1_xwidth = 0;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccpsd1_yheight = 0;
#line 51 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccpsd1_restore_neutron = 0;
#line 10658 "./ISIS_SANS2d.c"

  SIG_MESSAGE("psd1 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 10665 "./ISIS_SANS2d.c"
  rot_mul(mctr1, mcrotaisis_source, mcrotapsd1);
  rot_transpose(mcrotalmon1, mctr1);
  rot_mul(mcrotapsd1, mctr1, mcrotrpsd1);
  mctc1 = coords_set(
#line 78 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    0.0,
#line 78 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    0.0,
#line 78 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    3.699);
#line 10676 "./ISIS_SANS2d.c"
  rot_transpose(mcrotaisis_source, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposapsd1 = coords_add(mcposaisis_source, mctc2);
  mctc1 = coords_sub(mcposalmon1, mcposapsd1);
  mcposrpsd1 = rot_apply(mcrotapsd1, mctc1);
  mcDEBUG_COMPONENT("psd1", mcposapsd1, mcrotapsd1)
  mccomp_posa[4] = mcposapsd1;
  mccomp_posr[4] = mcposrpsd1;
  mcNCounter[4]  = mcPCounter[4] = mcP2Counter[4] = 0;
  mcAbsorbProp[4]= 0;
    /* Component bender1. */
  /* Setting parameters for component bender1. */
  SIG_MESSAGE("bender1 (Init:SetPar)");
#line 81 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender1_w1 = .0355;
#line 81 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender1_h1 = .020;
#line 81 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender1_w2 = .0355;
#line 81 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender1_h2 = .020;
#line 82 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender1_l = 0.3245;
#line 114 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender1_R0 = 0.995;
#line 114 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender1_Qc = 0.0218;
#line 114 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender1_alpha = 4.38;
#line 114 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender1_m = 1.0;
#line 114 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender1_W = 0.003;
#line 81 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender1_nslit = 9;
#line 81 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender1_d = .0005;
#line 82 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender1_mleft = 1;
#line 82 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender1_mright = 3;
#line 82 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender1_mtop = 1;
#line 82 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender1_mbottom = 1;
#line 115 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender1_nhslit = 1;
#line 115 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender1_G = 0;
#line 116 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender1_aleft = -1;
#line 116 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender1_aright = -1;
#line 116 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender1_atop = -1;
#line 116 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender1_abottom = -1;
#line 82 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender1_wavy = 0;
#line 117 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender1_wavy_z = 0;
#line 117 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender1_wavy_tb = 0;
#line 117 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender1_wavy_lr = 0;
#line 118 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender1_chamfers = 0;
#line 118 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender1_chamfers_z = 0;
#line 118 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender1_chamfers_lr = 0;
#line 118 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender1_chamfers_tb = 0;
#line 118 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender1_nelements = 1;
#line 119 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender1_nu = 0;
#line 119 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender1_phase = 0;
#line 119 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  if("NULL") strncpy(mccbender1_reflect, "NULL" ? "NULL" : "", 16384); else mccbender1_reflect[0]='\0';
#line 10758 "./ISIS_SANS2d.c"

  SIG_MESSAGE("bender1 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 84 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    (0.0)*DEG2RAD,
#line 84 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    (0.137099)*DEG2RAD,
#line 84 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    (0.0)*DEG2RAD);
#line 10768 "./ISIS_SANS2d.c"
  rot_mul(mctr1, mcrotaisis_source, mcrotabender1);
  rot_transpose(mcrotapsd1, mctr1);
  rot_mul(mcrotabender1, mctr1, mcrotrbender1);
  mctc1 = coords_set(
#line 83 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    0,
#line 83 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    0,
#line 83 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    3.7);
#line 10779 "./ISIS_SANS2d.c"
  rot_transpose(mcrotaisis_source, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposabender1 = coords_add(mcposaisis_source, mctc2);
  mctc1 = coords_sub(mcposapsd1, mcposabender1);
  mcposrbender1 = rot_apply(mcrotabender1, mctc1);
  mcDEBUG_COMPONENT("bender1", mcposabender1, mcrotabender1)
  mccomp_posa[5] = mcposabender1;
  mccomp_posr[5] = mcposrbender1;
  mcNCounter[5]  = mcPCounter[5] = mcP2Counter[5] = 0;
  mcAbsorbProp[5]= 0;
    /* Component bender2. */
  /* Setting parameters for component bender2. */
  SIG_MESSAGE("bender2 (Init:SetPar)");
#line 87 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender2_w1 = .0355;
#line 87 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender2_h1 = .020;
#line 87 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender2_w2 = .0355;
#line 87 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender2_h2 = .020;
#line 88 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender2_l = 0.3245;
#line 114 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender2_R0 = 0.995;
#line 114 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender2_Qc = 0.0218;
#line 114 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender2_alpha = 4.38;
#line 114 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender2_m = 1.0;
#line 114 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender2_W = 0.003;
#line 87 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender2_nslit = 9;
#line 87 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender2_d = .0005;
#line 88 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender2_mleft = 1;
#line 88 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender2_mright = 3;
#line 88 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender2_mtop = 1;
#line 88 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender2_mbottom = 1;
#line 115 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender2_nhslit = 1;
#line 115 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender2_G = 0;
#line 116 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender2_aleft = -1;
#line 116 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender2_aright = -1;
#line 116 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender2_atop = -1;
#line 116 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender2_abottom = -1;
#line 88 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender2_wavy = 0;
#line 117 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender2_wavy_z = 0;
#line 117 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender2_wavy_tb = 0;
#line 117 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender2_wavy_lr = 0;
#line 118 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender2_chamfers = 0;
#line 118 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender2_chamfers_z = 0;
#line 118 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender2_chamfers_lr = 0;
#line 118 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender2_chamfers_tb = 0;
#line 118 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender2_nelements = 1;
#line 119 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender2_nu = 0;
#line 119 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender2_phase = 0;
#line 119 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  if("NULL") strncpy(mccbender2_reflect, "NULL" ? "NULL" : "", 16384); else mccbender2_reflect[0]='\0';
#line 10861 "./ISIS_SANS2d.c"

  SIG_MESSAGE("bender2 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 90 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    (0.0)*DEG2RAD,
#line 90 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    (0.1375099)*DEG2RAD,
#line 90 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    (0.0)*DEG2RAD);
#line 10871 "./ISIS_SANS2d.c"
  rot_mul(mctr1, mcrotabender1, mcrotabender2);
  rot_transpose(mcrotabender1, mctr1);
  rot_mul(mcrotabender2, mctr1, mcrotrbender2);
  mctc1 = coords_set(
#line 89 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    0,
#line 89 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    0,
#line 89 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    0.325);
#line 10882 "./ISIS_SANS2d.c"
  rot_transpose(mcrotabender1, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposabender2 = coords_add(mcposabender1, mctc2);
  mctc1 = coords_sub(mcposabender1, mcposabender2);
  mcposrbender2 = rot_apply(mcrotabender2, mctc1);
  mcDEBUG_COMPONENT("bender2", mcposabender2, mcrotabender2)
  mccomp_posa[6] = mcposabender2;
  mccomp_posr[6] = mcposrbender2;
  mcNCounter[6]  = mcPCounter[6] = mcP2Counter[6] = 0;
  mcAbsorbProp[6]= 0;
    /* Component bender3. */
  /* Setting parameters for component bender3. */
  SIG_MESSAGE("bender3 (Init:SetPar)");
#line 93 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender3_w1 = .0355;
#line 93 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender3_h1 = .020;
#line 93 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender3_w2 = .0355;
#line 93 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender3_h2 = .020;
#line 94 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender3_l = 0.3245;
#line 114 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender3_R0 = 0.995;
#line 114 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender3_Qc = 0.0218;
#line 114 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender3_alpha = 4.38;
#line 114 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender3_m = 1.0;
#line 114 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender3_W = 0.003;
#line 93 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender3_nslit = 9;
#line 93 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender3_d = .0005;
#line 94 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender3_mleft = 1;
#line 94 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender3_mright = 3;
#line 94 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender3_mtop = 1;
#line 94 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender3_mbottom = 1;
#line 115 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender3_nhslit = 1;
#line 115 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender3_G = 0;
#line 116 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender3_aleft = -1;
#line 116 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender3_aright = -1;
#line 116 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender3_atop = -1;
#line 116 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender3_abottom = -1;
#line 94 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender3_wavy = 0;
#line 117 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender3_wavy_z = 0;
#line 117 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender3_wavy_tb = 0;
#line 117 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender3_wavy_lr = 0;
#line 118 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender3_chamfers = 0;
#line 118 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender3_chamfers_z = 0;
#line 118 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender3_chamfers_lr = 0;
#line 118 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender3_chamfers_tb = 0;
#line 118 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender3_nelements = 1;
#line 119 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender3_nu = 0;
#line 119 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender3_phase = 0;
#line 119 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  if("NULL") strncpy(mccbender3_reflect, "NULL" ? "NULL" : "", 16384); else mccbender3_reflect[0]='\0';
#line 10964 "./ISIS_SANS2d.c"

  SIG_MESSAGE("bender3 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 96 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    (0.0)*DEG2RAD,
#line 96 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    (0.1375099)*DEG2RAD,
#line 96 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    (0.0)*DEG2RAD);
#line 10974 "./ISIS_SANS2d.c"
  rot_mul(mctr1, mcrotabender2, mcrotabender3);
  rot_transpose(mcrotabender2, mctr1);
  rot_mul(mcrotabender3, mctr1, mcrotrbender3);
  mctc1 = coords_set(
#line 95 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    0,
#line 95 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    0,
#line 95 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    0.325);
#line 10985 "./ISIS_SANS2d.c"
  rot_transpose(mcrotabender2, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposabender3 = coords_add(mcposabender2, mctc2);
  mctc1 = coords_sub(mcposabender2, mcposabender3);
  mcposrbender3 = rot_apply(mcrotabender3, mctc1);
  mcDEBUG_COMPONENT("bender3", mcposabender3, mcrotabender3)
  mccomp_posa[7] = mcposabender3;
  mccomp_posr[7] = mcposrbender3;
  mcNCounter[7]  = mcPCounter[7] = mcP2Counter[7] = 0;
  mcAbsorbProp[7]= 0;
    /* Component bender4. */
  /* Setting parameters for component bender4. */
  SIG_MESSAGE("bender4 (Init:SetPar)");
#line 99 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender4_w1 = .0355;
#line 99 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender4_h1 = .020;
#line 99 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender4_w2 = .0355;
#line 99 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender4_h2 = .020;
#line 100 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender4_l = 0.3245;
#line 114 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender4_R0 = 0.995;
#line 114 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender4_Qc = 0.0218;
#line 114 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender4_alpha = 4.38;
#line 114 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender4_m = 1.0;
#line 114 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender4_W = 0.003;
#line 99 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender4_nslit = 9;
#line 99 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender4_d = .0005;
#line 100 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender4_mleft = 1;
#line 100 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender4_mright = 3;
#line 100 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender4_mtop = 1;
#line 100 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender4_mbottom = 1;
#line 115 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender4_nhslit = 1;
#line 115 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender4_G = 0;
#line 116 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender4_aleft = -1;
#line 116 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender4_aright = -1;
#line 116 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender4_atop = -1;
#line 116 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender4_abottom = -1;
#line 100 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender4_wavy = 0;
#line 117 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender4_wavy_z = 0;
#line 117 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender4_wavy_tb = 0;
#line 117 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender4_wavy_lr = 0;
#line 118 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender4_chamfers = 0;
#line 118 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender4_chamfers_z = 0;
#line 118 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender4_chamfers_lr = 0;
#line 118 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender4_chamfers_tb = 0;
#line 118 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender4_nelements = 1;
#line 119 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender4_nu = 0;
#line 119 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender4_phase = 0;
#line 119 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  if("NULL") strncpy(mccbender4_reflect, "NULL" ? "NULL" : "", 16384); else mccbender4_reflect[0]='\0';
#line 11067 "./ISIS_SANS2d.c"

  SIG_MESSAGE("bender4 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 102 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    (0.0)*DEG2RAD,
#line 102 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    (0.1375099)*DEG2RAD,
#line 102 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    (0.0)*DEG2RAD);
#line 11077 "./ISIS_SANS2d.c"
  rot_mul(mctr1, mcrotabender3, mcrotabender4);
  rot_transpose(mcrotabender3, mctr1);
  rot_mul(mcrotabender4, mctr1, mcrotrbender4);
  mctc1 = coords_set(
#line 101 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    0,
#line 101 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    0,
#line 101 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    0.325);
#line 11088 "./ISIS_SANS2d.c"
  rot_transpose(mcrotabender3, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposabender4 = coords_add(mcposabender3, mctc2);
  mctc1 = coords_sub(mcposabender3, mcposabender4);
  mcposrbender4 = rot_apply(mcrotabender4, mctc1);
  mcDEBUG_COMPONENT("bender4", mcposabender4, mcrotabender4)
  mccomp_posa[8] = mcposabender4;
  mccomp_posr[8] = mcposrbender4;
  mcNCounter[8]  = mcPCounter[8] = mcP2Counter[8] = 0;
  mcAbsorbProp[8]= 0;
    /* Component bender5. */
  /* Setting parameters for component bender5. */
  SIG_MESSAGE("bender5 (Init:SetPar)");
#line 105 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender5_w1 = .0355;
#line 105 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender5_h1 = .020;
#line 105 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender5_w2 = .0355;
#line 105 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender5_h2 = .020;
#line 106 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender5_l = 0.3245;
#line 114 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender5_R0 = 0.995;
#line 114 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender5_Qc = 0.0218;
#line 114 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender5_alpha = 4.38;
#line 114 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender5_m = 1.0;
#line 114 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender5_W = 0.003;
#line 105 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender5_nslit = 9;
#line 105 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender5_d = .0005;
#line 106 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender5_mleft = 1;
#line 106 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender5_mright = 3;
#line 106 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender5_mtop = 1;
#line 106 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender5_mbottom = 1;
#line 115 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender5_nhslit = 1;
#line 115 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender5_G = 0;
#line 116 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender5_aleft = -1;
#line 116 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender5_aright = -1;
#line 116 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender5_atop = -1;
#line 116 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender5_abottom = -1;
#line 106 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender5_wavy = 0;
#line 117 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender5_wavy_z = 0;
#line 117 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender5_wavy_tb = 0;
#line 117 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender5_wavy_lr = 0;
#line 118 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender5_chamfers = 0;
#line 118 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender5_chamfers_z = 0;
#line 118 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender5_chamfers_lr = 0;
#line 118 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender5_chamfers_tb = 0;
#line 118 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender5_nelements = 1;
#line 119 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender5_nu = 0;
#line 119 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender5_phase = 0;
#line 119 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  if("NULL") strncpy(mccbender5_reflect, "NULL" ? "NULL" : "", 16384); else mccbender5_reflect[0]='\0';
#line 11170 "./ISIS_SANS2d.c"

  SIG_MESSAGE("bender5 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 108 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    (0.0)*DEG2RAD,
#line 108 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    (0.1375099)*DEG2RAD,
#line 108 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    (0.0)*DEG2RAD);
#line 11180 "./ISIS_SANS2d.c"
  rot_mul(mctr1, mcrotabender4, mcrotabender5);
  rot_transpose(mcrotabender4, mctr1);
  rot_mul(mcrotabender5, mctr1, mcrotrbender5);
  mctc1 = coords_set(
#line 107 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    0,
#line 107 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    0,
#line 107 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    0.325);
#line 11191 "./ISIS_SANS2d.c"
  rot_transpose(mcrotabender4, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposabender5 = coords_add(mcposabender4, mctc2);
  mctc1 = coords_sub(mcposabender4, mcposabender5);
  mcposrbender5 = rot_apply(mcrotabender5, mctc1);
  mcDEBUG_COMPONENT("bender5", mcposabender5, mcrotabender5)
  mccomp_posa[9] = mcposabender5;
  mccomp_posr[9] = mcposrbender5;
  mcNCounter[9]  = mcPCounter[9] = mcP2Counter[9] = 0;
  mcAbsorbProp[9]= 0;
    /* Component bender6. */
  /* Setting parameters for component bender6. */
  SIG_MESSAGE("bender6 (Init:SetPar)");
#line 111 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender6_w1 = .0355;
#line 111 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender6_h1 = .020;
#line 111 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender6_w2 = .0355;
#line 111 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender6_h2 = .020;
#line 112 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender6_l = 0.3245;
#line 114 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender6_R0 = 0.995;
#line 114 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender6_Qc = 0.0218;
#line 114 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender6_alpha = 4.38;
#line 114 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender6_m = 1.0;
#line 114 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender6_W = 0.003;
#line 111 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender6_nslit = 9;
#line 111 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender6_d = .0005;
#line 112 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender6_mleft = 1;
#line 112 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender6_mright = 3;
#line 112 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender6_mtop = 1;
#line 112 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender6_mbottom = 1;
#line 115 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender6_nhslit = 1;
#line 115 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender6_G = 0;
#line 116 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender6_aleft = -1;
#line 116 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender6_aright = -1;
#line 116 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender6_atop = -1;
#line 116 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender6_abottom = -1;
#line 112 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender6_wavy = 0;
#line 117 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender6_wavy_z = 0;
#line 117 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender6_wavy_tb = 0;
#line 117 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender6_wavy_lr = 0;
#line 118 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender6_chamfers = 0;
#line 118 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender6_chamfers_z = 0;
#line 118 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender6_chamfers_lr = 0;
#line 118 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender6_chamfers_tb = 0;
#line 118 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender6_nelements = 1;
#line 119 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender6_nu = 0;
#line 119 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender6_phase = 0;
#line 119 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  if("NULL") strncpy(mccbender6_reflect, "NULL" ? "NULL" : "", 16384); else mccbender6_reflect[0]='\0';
#line 11273 "./ISIS_SANS2d.c"

  SIG_MESSAGE("bender6 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 114 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    (0.0)*DEG2RAD,
#line 114 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    (0.1375099)*DEG2RAD,
#line 114 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    (0.0)*DEG2RAD);
#line 11283 "./ISIS_SANS2d.c"
  rot_mul(mctr1, mcrotabender5, mcrotabender6);
  rot_transpose(mcrotabender5, mctr1);
  rot_mul(mcrotabender6, mctr1, mcrotrbender6);
  mctc1 = coords_set(
#line 113 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    0,
#line 113 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    0,
#line 113 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    0.325);
#line 11294 "./ISIS_SANS2d.c"
  rot_transpose(mcrotabender5, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposabender6 = coords_add(mcposabender5, mctc2);
  mctc1 = coords_sub(mcposabender5, mcposabender6);
  mcposrbender6 = rot_apply(mcrotabender6, mctc1);
  mcDEBUG_COMPONENT("bender6", mcposabender6, mcrotabender6)
  mccomp_posa[10] = mcposabender6;
  mccomp_posr[10] = mcposrbender6;
  mcNCounter[10]  = mcPCounter[10] = mcP2Counter[10] = 0;
  mcAbsorbProp[10]= 0;
    /* Component bender7. */
  /* Setting parameters for component bender7. */
  SIG_MESSAGE("bender7 (Init:SetPar)");
#line 117 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender7_w1 = .0355;
#line 117 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender7_h1 = .020;
#line 117 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender7_w2 = .0355;
#line 117 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender7_h2 = .020;
#line 118 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender7_l = 0.3245;
#line 114 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender7_R0 = 0.995;
#line 114 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender7_Qc = 0.0218;
#line 114 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender7_alpha = 4.38;
#line 114 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender7_m = 1.0;
#line 114 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender7_W = 0.003;
#line 117 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender7_nslit = 9;
#line 117 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender7_d = .0005;
#line 118 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender7_mleft = 1;
#line 118 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender7_mright = 3;
#line 118 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender7_mtop = 1;
#line 118 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender7_mbottom = 1;
#line 115 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender7_nhslit = 1;
#line 115 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender7_G = 0;
#line 116 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender7_aleft = -1;
#line 116 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender7_aright = -1;
#line 116 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender7_atop = -1;
#line 116 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender7_abottom = -1;
#line 118 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender7_wavy = 0;
#line 117 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender7_wavy_z = 0;
#line 117 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender7_wavy_tb = 0;
#line 117 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender7_wavy_lr = 0;
#line 118 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender7_chamfers = 0;
#line 118 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender7_chamfers_z = 0;
#line 118 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender7_chamfers_lr = 0;
#line 118 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender7_chamfers_tb = 0;
#line 118 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender7_nelements = 1;
#line 119 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender7_nu = 0;
#line 119 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender7_phase = 0;
#line 119 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  if("NULL") strncpy(mccbender7_reflect, "NULL" ? "NULL" : "", 16384); else mccbender7_reflect[0]='\0';
#line 11376 "./ISIS_SANS2d.c"

  SIG_MESSAGE("bender7 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 120 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    (0.0)*DEG2RAD,
#line 120 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    (0.1375099)*DEG2RAD,
#line 120 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    (0.0)*DEG2RAD);
#line 11386 "./ISIS_SANS2d.c"
  rot_mul(mctr1, mcrotabender6, mcrotabender7);
  rot_transpose(mcrotabender6, mctr1);
  rot_mul(mcrotabender7, mctr1, mcrotrbender7);
  mctc1 = coords_set(
#line 119 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    0,
#line 119 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    0,
#line 119 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    0.325);
#line 11397 "./ISIS_SANS2d.c"
  rot_transpose(mcrotabender6, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposabender7 = coords_add(mcposabender6, mctc2);
  mctc1 = coords_sub(mcposabender6, mcposabender7);
  mcposrbender7 = rot_apply(mcrotabender7, mctc1);
  mcDEBUG_COMPONENT("bender7", mcposabender7, mcrotabender7)
  mccomp_posa[11] = mcposabender7;
  mccomp_posr[11] = mcposrbender7;
  mcNCounter[11]  = mcPCounter[11] = mcP2Counter[11] = 0;
  mcAbsorbProp[11]= 0;
    /* Component bender8. */
  /* Setting parameters for component bender8. */
  SIG_MESSAGE("bender8 (Init:SetPar)");
#line 123 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender8_w1 = .0355;
#line 123 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender8_h1 = .020;
#line 123 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender8_w2 = .0355;
#line 123 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender8_h2 = .020;
#line 124 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender8_l = 0.3245;
#line 114 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender8_R0 = 0.995;
#line 114 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender8_Qc = 0.0218;
#line 114 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender8_alpha = 4.38;
#line 114 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender8_m = 1.0;
#line 114 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender8_W = 0.003;
#line 123 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender8_nslit = 9;
#line 123 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender8_d = .0005;
#line 124 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender8_mleft = 1;
#line 124 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender8_mright = 3;
#line 124 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender8_mtop = 1;
#line 124 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender8_mbottom = 1;
#line 115 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender8_nhslit = 1;
#line 115 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender8_G = 0;
#line 116 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender8_aleft = -1;
#line 116 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender8_aright = -1;
#line 116 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender8_atop = -1;
#line 116 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender8_abottom = -1;
#line 124 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender8_wavy = 0;
#line 117 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender8_wavy_z = 0;
#line 117 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender8_wavy_tb = 0;
#line 117 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender8_wavy_lr = 0;
#line 118 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender8_chamfers = 0;
#line 118 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender8_chamfers_z = 0;
#line 118 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender8_chamfers_lr = 0;
#line 118 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender8_chamfers_tb = 0;
#line 118 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender8_nelements = 1;
#line 119 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender8_nu = 0;
#line 119 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender8_phase = 0;
#line 119 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  if("NULL") strncpy(mccbender8_reflect, "NULL" ? "NULL" : "", 16384); else mccbender8_reflect[0]='\0';
#line 11479 "./ISIS_SANS2d.c"

  SIG_MESSAGE("bender8 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 126 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    (0.0)*DEG2RAD,
#line 126 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    (0.1375099)*DEG2RAD,
#line 126 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    (0.0)*DEG2RAD);
#line 11489 "./ISIS_SANS2d.c"
  rot_mul(mctr1, mcrotabender7, mcrotabender8);
  rot_transpose(mcrotabender7, mctr1);
  rot_mul(mcrotabender8, mctr1, mcrotrbender8);
  mctc1 = coords_set(
#line 125 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    0,
#line 125 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    0,
#line 125 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    0.325);
#line 11500 "./ISIS_SANS2d.c"
  rot_transpose(mcrotabender7, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposabender8 = coords_add(mcposabender7, mctc2);
  mctc1 = coords_sub(mcposabender7, mcposabender8);
  mcposrbender8 = rot_apply(mcrotabender8, mctc1);
  mcDEBUG_COMPONENT("bender8", mcposabender8, mcrotabender8)
  mccomp_posa[12] = mcposabender8;
  mccomp_posr[12] = mcposrbender8;
  mcNCounter[12]  = mcPCounter[12] = mcP2Counter[12] = 0;
  mcAbsorbProp[12]= 0;
    /* Component bender9. */
  /* Setting parameters for component bender9. */
  SIG_MESSAGE("bender9 (Init:SetPar)");
#line 129 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender9_w1 = .0355;
#line 129 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender9_h1 = .020;
#line 129 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender9_w2 = .0355;
#line 129 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender9_h2 = .020;
#line 130 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender9_l = 0.3245;
#line 114 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender9_R0 = 0.995;
#line 114 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender9_Qc = 0.0218;
#line 114 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender9_alpha = 4.38;
#line 114 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender9_m = 1.0;
#line 114 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender9_W = 0.003;
#line 129 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender9_nslit = 9;
#line 129 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender9_d = .0005;
#line 130 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender9_mleft = 1;
#line 130 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender9_mright = 3;
#line 130 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender9_mtop = 1;
#line 130 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender9_mbottom = 1;
#line 115 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender9_nhslit = 1;
#line 115 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender9_G = 0;
#line 116 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender9_aleft = -1;
#line 116 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender9_aright = -1;
#line 116 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender9_atop = -1;
#line 116 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender9_abottom = -1;
#line 130 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender9_wavy = 0;
#line 117 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender9_wavy_z = 0;
#line 117 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender9_wavy_tb = 0;
#line 117 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender9_wavy_lr = 0;
#line 118 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender9_chamfers = 0;
#line 118 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender9_chamfers_z = 0;
#line 118 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender9_chamfers_lr = 0;
#line 118 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender9_chamfers_tb = 0;
#line 118 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender9_nelements = 1;
#line 119 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender9_nu = 0;
#line 119 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender9_phase = 0;
#line 119 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  if("NULL") strncpy(mccbender9_reflect, "NULL" ? "NULL" : "", 16384); else mccbender9_reflect[0]='\0';
#line 11582 "./ISIS_SANS2d.c"

  SIG_MESSAGE("bender9 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 132 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    (0.0)*DEG2RAD,
#line 132 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    (0.1375099)*DEG2RAD,
#line 132 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    (0.0)*DEG2RAD);
#line 11592 "./ISIS_SANS2d.c"
  rot_mul(mctr1, mcrotabender8, mcrotabender9);
  rot_transpose(mcrotabender8, mctr1);
  rot_mul(mcrotabender9, mctr1, mcrotrbender9);
  mctc1 = coords_set(
#line 131 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    0,
#line 131 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    0,
#line 131 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    0.325);
#line 11603 "./ISIS_SANS2d.c"
  rot_transpose(mcrotabender8, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposabender9 = coords_add(mcposabender8, mctc2);
  mctc1 = coords_sub(mcposabender8, mcposabender9);
  mcposrbender9 = rot_apply(mcrotabender9, mctc1);
  mcDEBUG_COMPONENT("bender9", mcposabender9, mcrotabender9)
  mccomp_posa[13] = mcposabender9;
  mccomp_posr[13] = mcposrbender9;
  mcNCounter[13]  = mcPCounter[13] = mcP2Counter[13] = 0;
  mcAbsorbProp[13]= 0;
    /* Component bender10. */
  /* Setting parameters for component bender10. */
  SIG_MESSAGE("bender10 (Init:SetPar)");
#line 135 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender10_w1 = .0355;
#line 135 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender10_h1 = .020;
#line 135 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender10_w2 = .0355;
#line 135 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender10_h2 = .020;
#line 136 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender10_l = 0.3245;
#line 114 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender10_R0 = 0.995;
#line 114 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender10_Qc = 0.0218;
#line 114 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender10_alpha = 4.38;
#line 114 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender10_m = 1.0;
#line 114 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender10_W = 0.003;
#line 135 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender10_nslit = 9;
#line 135 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender10_d = .0005;
#line 136 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender10_mleft = 1;
#line 136 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender10_mright = 3;
#line 136 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender10_mtop = 1;
#line 136 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender10_mbottom = 1;
#line 115 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender10_nhslit = 1;
#line 115 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender10_G = 0;
#line 116 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender10_aleft = -1;
#line 116 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender10_aright = -1;
#line 116 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender10_atop = -1;
#line 116 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender10_abottom = -1;
#line 136 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender10_wavy = 0;
#line 117 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender10_wavy_z = 0;
#line 117 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender10_wavy_tb = 0;
#line 117 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender10_wavy_lr = 0;
#line 118 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender10_chamfers = 0;
#line 118 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender10_chamfers_z = 0;
#line 118 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender10_chamfers_lr = 0;
#line 118 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender10_chamfers_tb = 0;
#line 118 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender10_nelements = 1;
#line 119 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender10_nu = 0;
#line 119 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccbender10_phase = 0;
#line 119 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  if("NULL") strncpy(mccbender10_reflect, "NULL" ? "NULL" : "", 16384); else mccbender10_reflect[0]='\0';
#line 11685 "./ISIS_SANS2d.c"

  SIG_MESSAGE("bender10 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 138 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    (0.0)*DEG2RAD,
#line 138 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    (0.1375099)*DEG2RAD,
#line 138 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    (0.0)*DEG2RAD);
#line 11695 "./ISIS_SANS2d.c"
  rot_mul(mctr1, mcrotabender9, mcrotabender10);
  rot_transpose(mcrotabender9, mctr1);
  rot_mul(mcrotabender10, mctr1, mcrotrbender10);
  mctc1 = coords_set(
#line 137 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    0,
#line 137 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    0,
#line 137 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    0.325);
#line 11706 "./ISIS_SANS2d.c"
  rot_transpose(mcrotabender9, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposabender10 = coords_add(mcposabender9, mctc2);
  mctc1 = coords_sub(mcposabender9, mcposabender10);
  mcposrbender10 = rot_apply(mcrotabender10, mctc1);
  mcDEBUG_COMPONENT("bender10", mcposabender10, mcrotabender10)
  mccomp_posa[14] = mcposabender10;
  mccomp_posr[14] = mcposrbender10;
  mcNCounter[14]  = mcPCounter[14] = mcP2Counter[14] = 0;
  mcAbsorbProp[14]= 0;
    /* Component lmonb. */
  /* Setting parameters for component lmonb. */
  SIG_MESSAGE("lmonb (Init:SetPar)");
#line 141 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  if("lmonB.dat") strncpy(mcclmonb_filename, "lmonB.dat" ? "lmonB.dat" : "", 16384); else mcclmonb_filename[0]='\0';
#line 141 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mcclmonb_xmin = -0.018;
#line 142 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mcclmonb_xmax = 0.018;
#line 142 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mcclmonb_ymin = -0.018;
#line 142 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mcclmonb_ymax = 0.018;
#line 51 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mcclmonb_xwidth = 0;
#line 51 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mcclmonb_yheight = 0;
#line 142 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mcclmonb_Lmin = 0.0;
#line 143 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mcclmonb_Lmax = 17.0;
#line 51 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mcclmonb_restore_neutron = 0;
#line 51 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mcclmonb_nowritefile = 0;
#line 11742 "./ISIS_SANS2d.c"

  SIG_MESSAGE("lmonb (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 11749 "./ISIS_SANS2d.c"
  rot_mul(mctr1, mcrotabender10, mcrotalmonb);
  rot_transpose(mcrotabender10, mctr1);
  rot_mul(mcrotalmonb, mctr1, mcrotrlmonb);
  mctc1 = coords_set(
#line 144 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    0.0,
#line 144 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    0.0,
#line 144 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    0.326);
#line 11760 "./ISIS_SANS2d.c"
  rot_transpose(mcrotabender10, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposalmonb = coords_add(mcposabender10, mctc2);
  mctc1 = coords_sub(mcposabender10, mcposalmonb);
  mcposrlmonb = rot_apply(mcrotalmonb, mctc1);
  mcDEBUG_COMPONENT("lmonb", mcposalmonb, mcrotalmonb)
  mccomp_posa[15] = mcposalmonb;
  mccomp_posr[15] = mcposrlmonb;
  mcNCounter[15]  = mcPCounter[15] = mcP2Counter[15] = 0;
  mcAbsorbProp[15]= 0;
    /* Component psd2. */
  /* Setting parameters for component psd2. */
  SIG_MESSAGE("psd2 (Init:SetPar)");
#line 147 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccpsd2_nx = 100;
#line 147 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccpsd2_ny = 100;
#line 147 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  if("psd2.dat") strncpy(mccpsd2_filename, "psd2.dat" ? "psd2.dat" : "", 16384); else mccpsd2_filename[0]='\0';
#line 147 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccpsd2_xmin = -0.025;
#line 148 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccpsd2_xmax = 0.025;
#line 148 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccpsd2_ymin = -0.025;
#line 148 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccpsd2_ymax = 0.025;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccpsd2_xwidth = 0;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccpsd2_yheight = 0;
#line 51 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccpsd2_restore_neutron = 0;
#line 11794 "./ISIS_SANS2d.c"

  SIG_MESSAGE("psd2 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 11801 "./ISIS_SANS2d.c"
  rot_mul(mctr1, mcrotalmonb, mcrotapsd2);
  rot_transpose(mcrotalmonb, mctr1);
  rot_mul(mcrotapsd2, mctr1, mcrotrpsd2);
  mctc1 = coords_set(
#line 149 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    0.0,
#line 149 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    0.0,
#line 149 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    0.001);
#line 11812 "./ISIS_SANS2d.c"
  rot_transpose(mcrotalmonb, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposapsd2 = coords_add(mcposalmonb, mctc2);
  mctc1 = coords_sub(mcposalmonb, mcposapsd2);
  mcposrpsd2 = rot_apply(mcrotapsd2, mctc1);
  mcDEBUG_COMPONENT("psd2", mcposapsd2, mcrotapsd2)
  mccomp_posa[16] = mcposapsd2;
  mccomp_posr[16] = mcposrpsd2;
  mcNCounter[16]  = mcPCounter[16] = mcP2Counter[16] = 0;
  mcAbsorbProp[16]= 0;
    /* Component guide_in. */
  /* Setting parameters for component guide_in. */
  SIG_MESSAGE("guide_in (Init:SetPar)");
#line 153 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_in_xmin = -0.015;
#line 153 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_in_xmax = 0.015;
#line 153 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_in_ymin = -.01;
#line 153 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_in_ymax = + .01;
#line 46 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_in_radius = 0;
#line 46 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_in_xwidth = 0;
#line 46 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_in_yheight = 0;
#line 11840 "./ISIS_SANS2d.c"

  SIG_MESSAGE("guide_in (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 11847 "./ISIS_SANS2d.c"
  rot_mul(mctr1, mcrotapsd2, mcrotaguide_in);
  rot_transpose(mcrotapsd2, mctr1);
  rot_mul(mcrotaguide_in, mctr1, mcrotrguide_in);
  mctc1 = coords_set(
#line 154 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    0,
#line 154 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    0,
#line 154 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    0.2845);
#line 11858 "./ISIS_SANS2d.c"
  rot_transpose(mcrotapsd2, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaguide_in = coords_add(mcposapsd2, mctc2);
  mctc1 = coords_sub(mcposapsd2, mcposaguide_in);
  mcposrguide_in = rot_apply(mcrotaguide_in, mctc1);
  mcDEBUG_COMPONENT("guide_in", mcposaguide_in, mcrotaguide_in)
  mccomp_posa[17] = mcposaguide_in;
  mccomp_posr[17] = mcposrguide_in;
  mcNCounter[17]  = mcPCounter[17] = mcP2Counter[17] = 0;
  mcAbsorbProp[17]= 0;
    /* Component guide_straight1. */
  /* Setting parameters for component guide_straight1. */
  SIG_MESSAGE("guide_straight1 (Init:SetPar)");
#line 160 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight1_w1 = .030;
#line 160 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight1_h1 = .020;
#line 160 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight1_w2 = .030;
#line 160 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight1_h2 = .020;
#line 160 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight1_l = 1.985;
#line 114 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight1_R0 = 0.995;
#line 114 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight1_Qc = 0.0218;
#line 114 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight1_alpha = 4.38;
#line 114 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight1_m = 1.0;
#line 114 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight1_W = 0.003;
#line 114 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight1_nslit = 1;
#line 114 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight1_d = 0.0005;
#line 161 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight1_mleft = 1;
#line 161 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight1_mright = 1;
#line 161 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight1_mtop = 1;
#line 161 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight1_mbottom = 1;
#line 115 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight1_nhslit = 1;
#line 115 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight1_G = 0;
#line 116 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight1_aleft = -1;
#line 116 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight1_aright = -1;
#line 116 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight1_atop = -1;
#line 116 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight1_abottom = -1;
#line 161 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight1_wavy = 0;
#line 117 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight1_wavy_z = 0;
#line 117 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight1_wavy_tb = 0;
#line 117 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight1_wavy_lr = 0;
#line 118 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight1_chamfers = 0;
#line 118 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight1_chamfers_z = 0;
#line 118 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight1_chamfers_lr = 0;
#line 118 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight1_chamfers_tb = 0;
#line 118 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight1_nelements = 1;
#line 119 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight1_nu = 0;
#line 119 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight1_phase = 0;
#line 119 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  if("NULL") strncpy(mccguide_straight1_reflect, "NULL" ? "NULL" : "", 16384); else mccguide_straight1_reflect[0]='\0';
#line 11940 "./ISIS_SANS2d.c"

  SIG_MESSAGE("guide_straight1 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 11947 "./ISIS_SANS2d.c"
  rot_mul(mctr1, mcrotaguide_in, mcrotaguide_straight1);
  rot_transpose(mcrotaguide_in, mctr1);
  rot_mul(mcrotaguide_straight1, mctr1, mcrotrguide_straight1);
  mctc1 = coords_set(
#line 162 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    0,
#line 162 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    0,
#line 162 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    0.0075);
#line 11958 "./ISIS_SANS2d.c"
  rot_transpose(mcrotaguide_in, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaguide_straight1 = coords_add(mcposaguide_in, mctc2);
  mctc1 = coords_sub(mcposaguide_in, mcposaguide_straight1);
  mcposrguide_straight1 = rot_apply(mcrotaguide_straight1, mctc1);
  mcDEBUG_COMPONENT("guide_straight1", mcposaguide_straight1, mcrotaguide_straight1)
  mccomp_posa[18] = mcposaguide_straight1;
  mccomp_posr[18] = mcposrguide_straight1;
  mcNCounter[18]  = mcPCounter[18] = mcP2Counter[18] = 0;
  mcAbsorbProp[18]= 0;
    /* Component guide_straight2. */
  /* Setting parameters for component guide_straight2. */
  SIG_MESSAGE("guide_straight2 (Init:SetPar)");
#line 165 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight2_w1 = .030;
#line 165 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight2_h1 = .020;
#line 165 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight2_w2 = .030;
#line 165 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight2_h2 = .020;
#line 165 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight2_l = 1.985;
#line 114 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight2_R0 = 0.995;
#line 114 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight2_Qc = 0.0218;
#line 114 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight2_alpha = 4.38;
#line 114 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight2_m = 1.0;
#line 114 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight2_W = 0.003;
#line 114 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight2_nslit = 1;
#line 114 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight2_d = 0.0005;
#line 166 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight2_mleft = 1;
#line 166 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight2_mright = 1;
#line 166 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight2_mtop = 1;
#line 166 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight2_mbottom = 1;
#line 115 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight2_nhslit = 1;
#line 115 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight2_G = 0;
#line 116 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight2_aleft = -1;
#line 116 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight2_aright = -1;
#line 116 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight2_atop = -1;
#line 116 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight2_abottom = -1;
#line 166 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight2_wavy = 0;
#line 117 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight2_wavy_z = 0;
#line 117 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight2_wavy_tb = 0;
#line 117 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight2_wavy_lr = 0;
#line 118 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight2_chamfers = 0;
#line 118 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight2_chamfers_z = 0;
#line 118 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight2_chamfers_lr = 0;
#line 118 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight2_chamfers_tb = 0;
#line 118 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight2_nelements = 1;
#line 119 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight2_nu = 0;
#line 119 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight2_phase = 0;
#line 119 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  if("NULL") strncpy(mccguide_straight2_reflect, "NULL" ? "NULL" : "", 16384); else mccguide_straight2_reflect[0]='\0';
#line 12040 "./ISIS_SANS2d.c"

  SIG_MESSAGE("guide_straight2 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 12047 "./ISIS_SANS2d.c"
  rot_mul(mctr1, mcrotaguide_straight1, mcrotaguide_straight2);
  rot_transpose(mcrotaguide_straight1, mctr1);
  rot_mul(mcrotaguide_straight2, mctr1, mcrotrguide_straight2);
  mctc1 = coords_set(
#line 167 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    0,
#line 167 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    0,
#line 167 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    2.000);
#line 12058 "./ISIS_SANS2d.c"
  rot_transpose(mcrotaguide_straight1, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaguide_straight2 = coords_add(mcposaguide_straight1, mctc2);
  mctc1 = coords_sub(mcposaguide_straight1, mcposaguide_straight2);
  mcposrguide_straight2 = rot_apply(mcrotaguide_straight2, mctc1);
  mcDEBUG_COMPONENT("guide_straight2", mcposaguide_straight2, mcrotaguide_straight2)
  mccomp_posa[19] = mcposaguide_straight2;
  mccomp_posr[19] = mcposrguide_straight2;
  mcNCounter[19]  = mcPCounter[19] = mcP2Counter[19] = 0;
  mcAbsorbProp[19]= 0;
    /* Component guide_straight3. */
  /* Setting parameters for component guide_straight3. */
  SIG_MESSAGE("guide_straight3 (Init:SetPar)");
#line 170 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight3_w1 = .030;
#line 170 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight3_h1 = .020;
#line 170 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight3_w2 = .030;
#line 170 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight3_h2 = .020;
#line 170 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight3_l = 1.985;
#line 114 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight3_R0 = 0.995;
#line 114 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight3_Qc = 0.0218;
#line 114 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight3_alpha = 4.38;
#line 114 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight3_m = 1.0;
#line 114 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight3_W = 0.003;
#line 114 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight3_nslit = 1;
#line 114 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight3_d = 0.0005;
#line 171 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight3_mleft = 1;
#line 171 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight3_mright = 1;
#line 171 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight3_mtop = 1;
#line 171 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight3_mbottom = 1;
#line 115 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight3_nhslit = 1;
#line 115 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight3_G = 0;
#line 116 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight3_aleft = -1;
#line 116 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight3_aright = -1;
#line 116 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight3_atop = -1;
#line 116 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight3_abottom = -1;
#line 171 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight3_wavy = 0;
#line 117 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight3_wavy_z = 0;
#line 117 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight3_wavy_tb = 0;
#line 117 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight3_wavy_lr = 0;
#line 118 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight3_chamfers = 0;
#line 118 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight3_chamfers_z = 0;
#line 118 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight3_chamfers_lr = 0;
#line 118 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight3_chamfers_tb = 0;
#line 118 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight3_nelements = 1;
#line 119 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight3_nu = 0;
#line 119 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight3_phase = 0;
#line 119 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  if("NULL") strncpy(mccguide_straight3_reflect, "NULL" ? "NULL" : "", 16384); else mccguide_straight3_reflect[0]='\0';
#line 12140 "./ISIS_SANS2d.c"

  SIG_MESSAGE("guide_straight3 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 12147 "./ISIS_SANS2d.c"
  rot_mul(mctr1, mcrotaguide_straight2, mcrotaguide_straight3);
  rot_transpose(mcrotaguide_straight2, mctr1);
  rot_mul(mcrotaguide_straight3, mctr1, mcrotrguide_straight3);
  mctc1 = coords_set(
#line 172 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    0,
#line 172 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    0,
#line 172 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    2.000);
#line 12158 "./ISIS_SANS2d.c"
  rot_transpose(mcrotaguide_straight2, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaguide_straight3 = coords_add(mcposaguide_straight2, mctc2);
  mctc1 = coords_sub(mcposaguide_straight2, mcposaguide_straight3);
  mcposrguide_straight3 = rot_apply(mcrotaguide_straight3, mctc1);
  mcDEBUG_COMPONENT("guide_straight3", mcposaguide_straight3, mcrotaguide_straight3)
  mccomp_posa[20] = mcposaguide_straight3;
  mccomp_posr[20] = mcposrguide_straight3;
  mcNCounter[20]  = mcPCounter[20] = mcP2Counter[20] = 0;
  mcAbsorbProp[20]= 0;
    /* Component guide_straight4. */
  /* Setting parameters for component guide_straight4. */
  SIG_MESSAGE("guide_straight4 (Init:SetPar)");
#line 175 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight4_w1 = .030;
#line 175 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight4_h1 = .020;
#line 175 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight4_w2 = .030;
#line 175 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight4_h2 = .020;
#line 175 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight4_l = 1.985;
#line 114 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight4_R0 = 0.995;
#line 114 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight4_Qc = 0.0218;
#line 114 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight4_alpha = 4.38;
#line 114 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight4_m = 1.0;
#line 114 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight4_W = 0.003;
#line 114 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight4_nslit = 1;
#line 114 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight4_d = 0.0005;
#line 176 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight4_mleft = 1;
#line 176 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight4_mright = 1;
#line 176 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight4_mtop = 1;
#line 176 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight4_mbottom = 1;
#line 115 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight4_nhslit = 1;
#line 115 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight4_G = 0;
#line 116 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight4_aleft = -1;
#line 116 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight4_aright = -1;
#line 116 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight4_atop = -1;
#line 116 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight4_abottom = -1;
#line 176 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight4_wavy = 0;
#line 117 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight4_wavy_z = 0;
#line 117 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight4_wavy_tb = 0;
#line 117 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight4_wavy_lr = 0;
#line 118 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight4_chamfers = 0;
#line 118 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight4_chamfers_z = 0;
#line 118 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight4_chamfers_lr = 0;
#line 118 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight4_chamfers_tb = 0;
#line 118 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight4_nelements = 1;
#line 119 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight4_nu = 0;
#line 119 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccguide_straight4_phase = 0;
#line 119 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  if("NULL") strncpy(mccguide_straight4_reflect, "NULL" ? "NULL" : "", 16384); else mccguide_straight4_reflect[0]='\0';
#line 12240 "./ISIS_SANS2d.c"

  SIG_MESSAGE("guide_straight4 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 12247 "./ISIS_SANS2d.c"
  rot_mul(mctr1, mcrotaguide_straight3, mcrotaguide_straight4);
  rot_transpose(mcrotaguide_straight3, mctr1);
  rot_mul(mcrotaguide_straight4, mctr1, mcrotrguide_straight4);
  mctc1 = coords_set(
#line 177 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    0,
#line 177 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    0,
#line 177 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    2.000);
#line 12258 "./ISIS_SANS2d.c"
  rot_transpose(mcrotaguide_straight3, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaguide_straight4 = coords_add(mcposaguide_straight3, mctc2);
  mctc1 = coords_sub(mcposaguide_straight3, mcposaguide_straight4);
  mcposrguide_straight4 = rot_apply(mcrotaguide_straight4, mctc1);
  mcDEBUG_COMPONENT("guide_straight4", mcposaguide_straight4, mcrotaguide_straight4)
  mccomp_posa[21] = mcposaguide_straight4;
  mccomp_posr[21] = mcposrguide_straight4;
  mcNCounter[21]  = mcPCounter[21] = mcP2Counter[21] = 0;
  mcAbsorbProp[21]= 0;
    /* Component psd3. */
  /* Setting parameters for component psd3. */
  SIG_MESSAGE("psd3 (Init:SetPar)");
#line 180 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccpsd3_nx = 100;
#line 180 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccpsd3_ny = 100;
#line 180 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  if("psd3.dat") strncpy(mccpsd3_filename, "psd3.dat" ? "psd3.dat" : "", 16384); else mccpsd3_filename[0]='\0';
#line 180 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccpsd3_xmin = -0.030;
#line 181 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccpsd3_xmax = 0.030;
#line 181 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccpsd3_ymin = -0.030;
#line 181 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccpsd3_ymax = 0.030;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccpsd3_xwidth = 0;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccpsd3_yheight = 0;
#line 51 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccpsd3_restore_neutron = 0;
#line 12292 "./ISIS_SANS2d.c"

  SIG_MESSAGE("psd3 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 12299 "./ISIS_SANS2d.c"
  rot_mul(mctr1, mcrotaguide_in, mcrotapsd3);
  rot_transpose(mcrotaguide_straight4, mctr1);
  rot_mul(mcrotapsd3, mctr1, mcrotrpsd3);
  mctc1 = coords_set(
#line 182 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    0.0,
#line 182 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    0.0,
#line 182 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    7.999);
#line 12310 "./ISIS_SANS2d.c"
  rot_transpose(mcrotaguide_in, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposapsd3 = coords_add(mcposaguide_in, mctc2);
  mctc1 = coords_sub(mcposaguide_straight4, mcposapsd3);
  mcposrpsd3 = rot_apply(mcrotapsd3, mctc1);
  mcDEBUG_COMPONENT("psd3", mcposapsd3, mcrotapsd3)
  mccomp_posa[22] = mcposapsd3;
  mccomp_posr[22] = mcposrpsd3;
  mcNCounter[22]  = mcPCounter[22] = mcP2Counter[22] = 0;
  mcAbsorbProp[22]= 0;
    /* Component aperture1. */
  /* Setting parameters for component aperture1. */
  SIG_MESSAGE("aperture1 (Init:SetPar)");
#line 46 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccaperture1_xmin = 0;
#line 46 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccaperture1_xmax = 0;
#line 46 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccaperture1_ymin = 0;
#line 46 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccaperture1_ymax = 0;
#line 46 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccaperture1_radius = 0;
#line 186 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccaperture1_xwidth = mcipA1w;
#line 186 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccaperture1_yheight = mcipA1h;
#line 12338 "./ISIS_SANS2d.c"

  SIG_MESSAGE("aperture1 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 12345 "./ISIS_SANS2d.c"
  rot_mul(mctr1, mcrotaguide_in, mcrotaaperture1);
  rot_transpose(mcrotapsd3, mctr1);
  rot_mul(mcrotaaperture1, mctr1, mcrotraperture1);
  mctc1 = coords_set(
#line 187 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    0,
#line 187 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    0,
#line 187 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    8.000);
#line 12356 "./ISIS_SANS2d.c"
  rot_transpose(mcrotaguide_in, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaaperture1 = coords_add(mcposaguide_in, mctc2);
  mctc1 = coords_sub(mcposapsd3, mcposaaperture1);
  mcposraperture1 = rot_apply(mcrotaaperture1, mctc1);
  mcDEBUG_COMPONENT("aperture1", mcposaaperture1, mcrotaaperture1)
  mccomp_posa[23] = mcposaaperture1;
  mccomp_posr[23] = mcposraperture1;
  mcNCounter[23]  = mcPCounter[23] = mcP2Counter[23] = 0;
  mcAbsorbProp[23]= 0;
    /* Component lmonitor2. */
  /* Setting parameters for component lmonitor2. */
  SIG_MESSAGE("lmonitor2 (Init:SetPar)");
#line 191 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  if("lmonitor2.dat") strncpy(mcclmonitor2_filename, "lmonitor2.dat" ? "lmonitor2.dat" : "", 16384); else mcclmonitor2_filename[0]='\0';
#line 191 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mcclmonitor2_xmin = -0.0155;
#line 192 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mcclmonitor2_xmax = 0.0155;
#line 192 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mcclmonitor2_ymin = -0.0105;
#line 192 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mcclmonitor2_ymax = 0.0105;
#line 51 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mcclmonitor2_xwidth = 0;
#line 51 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mcclmonitor2_yheight = 0;
#line 192 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mcclmonitor2_Lmin = 0.0;
#line 193 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mcclmonitor2_Lmax = 17.0;
#line 51 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mcclmonitor2_restore_neutron = 0;
#line 51 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mcclmonitor2_nowritefile = 0;
#line 12392 "./ISIS_SANS2d.c"

  SIG_MESSAGE("lmonitor2 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 12399 "./ISIS_SANS2d.c"
  rot_mul(mctr1, mcrotaaperture1, mcrotalmonitor2);
  rot_transpose(mcrotaaperture1, mctr1);
  rot_mul(mcrotalmonitor2, mctr1, mcrotrlmonitor2);
  mctc1 = coords_set(
#line 194 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    0.0,
#line 194 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    0.0,
#line 194 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    2.651);
#line 12410 "./ISIS_SANS2d.c"
  rot_transpose(mcrotaaperture1, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposalmonitor2 = coords_add(mcposaaperture1, mctc2);
  mctc1 = coords_sub(mcposaaperture1, mcposalmonitor2);
  mcposrlmonitor2 = rot_apply(mcrotalmonitor2, mctc1);
  mcDEBUG_COMPONENT("lmonitor2", mcposalmonitor2, mcrotalmonitor2)
  mccomp_posa[24] = mcposalmonitor2;
  mccomp_posr[24] = mcposrlmonitor2;
  mcNCounter[24]  = mcPCounter[24] = mcP2Counter[24] = 0;
  mcAbsorbProp[24]= 0;
    /* Component S6. */
  /* Setting parameters for component S6. */
  SIG_MESSAGE("S6 (Init:SetPar)");
#line 46 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccS6_xmin = 0;
#line 46 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccS6_xmax = 0;
#line 46 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccS6_ymin = 0;
#line 46 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccS6_ymax = 0;
#line 198 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccS6_radius = mcipS6;
#line 46 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccS6_xwidth = 0;
#line 46 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccS6_yheight = 0;
#line 12438 "./ISIS_SANS2d.c"

  SIG_MESSAGE("S6 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 12445 "./ISIS_SANS2d.c"
  rot_mul(mctr1, mcrotaaperture1, mcrotaS6);
  rot_transpose(mcrotalmonitor2, mctr1);
  rot_mul(mcrotaS6, mctr1, mcrotrS6);
  mctc1 = coords_set(
#line 199 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    0,
#line 199 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    0,
#line 199 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    2.800);
#line 12456 "./ISIS_SANS2d.c"
  rot_transpose(mcrotaaperture1, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaS6 = coords_add(mcposaaperture1, mctc2);
  mctc1 = coords_sub(mcposalmonitor2, mcposaS6);
  mcposrS6 = rot_apply(mcrotaS6, mctc1);
  mcDEBUG_COMPONENT("S6", mcposaS6, mcrotaS6)
  mccomp_posa[25] = mcposaS6;
  mccomp_posr[25] = mcposrS6;
  mcNCounter[25]  = mcPCounter[25] = mcP2Counter[25] = 0;
  mcAbsorbProp[25]= 0;
    /* Component APERTURE2. */
  /* Setting parameters for component APERTURE2. */
  SIG_MESSAGE("APERTURE2 (Init:SetPar)");
#line 46 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccAPERTURE2_xmin = 0;
#line 46 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccAPERTURE2_xmax = 0;
#line 46 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccAPERTURE2_ymin = 0;
#line 46 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccAPERTURE2_ymax = 0;
#line 204 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccAPERTURE2_radius = mcipA2;
#line 46 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccAPERTURE2_xwidth = 0;
#line 46 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccAPERTURE2_yheight = 0;
#line 12484 "./ISIS_SANS2d.c"

  SIG_MESSAGE("APERTURE2 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 12491 "./ISIS_SANS2d.c"
  rot_mul(mctr1, mcrotaaperture1, mcrotaAPERTURE2);
  rot_transpose(mcrotaS6, mctr1);
  rot_mul(mcrotaAPERTURE2, mctr1, mcrotrAPERTURE2);
  mctc1 = coords_set(
#line 205 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    0,
#line 205 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    0,
#line 205 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    mcipL1);
#line 12502 "./ISIS_SANS2d.c"
  rot_transpose(mcrotaaperture1, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaAPERTURE2 = coords_add(mcposaaperture1, mctc2);
  mctc1 = coords_sub(mcposaS6, mcposaAPERTURE2);
  mcposrAPERTURE2 = rot_apply(mcrotaAPERTURE2, mctc1);
  mcDEBUG_COMPONENT("APERTURE2", mcposaAPERTURE2, mcrotaAPERTURE2)
  mccomp_posa[26] = mcposaAPERTURE2;
  mccomp_posr[26] = mcposrAPERTURE2;
  mcNCounter[26]  = mcPCounter[26] = mcP2Counter[26] = 0;
  mcAbsorbProp[26]= 0;
    /* Component lmon2. */
  /* Setting parameters for component lmon2. */
  SIG_MESSAGE("lmon2 (Init:SetPar)");
#line 208 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  if("lmonitor3.dat") strncpy(mcclmon2_filename, "lmonitor3.dat" ? "lmonitor3.dat" : "", 16384); else mcclmon2_filename[0]='\0';
#line 208 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mcclmon2_xmin = -0.0075;
#line 209 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mcclmon2_xmax = 0.0075;
#line 209 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mcclmon2_ymin = -0.0075;
#line 209 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mcclmon2_ymax = 0.0075;
#line 51 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mcclmon2_xwidth = 0;
#line 51 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mcclmon2_yheight = 0;
#line 209 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mcclmon2_Lmin = 0.0;
#line 210 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mcclmon2_Lmax = 17.0;
#line 51 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mcclmon2_restore_neutron = 0;
#line 51 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mcclmon2_nowritefile = 0;
#line 12538 "./ISIS_SANS2d.c"

  SIG_MESSAGE("lmon2 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 12545 "./ISIS_SANS2d.c"
  rot_mul(mctr1, mcrotaAPERTURE2, mcrotalmon2);
  rot_transpose(mcrotaAPERTURE2, mctr1);
  rot_mul(mcrotalmon2, mctr1, mcrotrlmon2);
  mctc1 = coords_set(
#line 211 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    0.0,
#line 211 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    0.0,
#line 211 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    0.285);
#line 12556 "./ISIS_SANS2d.c"
  rot_transpose(mcrotaAPERTURE2, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposalmon2 = coords_add(mcposaAPERTURE2, mctc2);
  mctc1 = coords_sub(mcposaAPERTURE2, mcposalmon2);
  mcposrlmon2 = rot_apply(mcrotalmon2, mctc1);
  mcDEBUG_COMPONENT("lmon2", mcposalmon2, mcrotalmon2)
  mccomp_posa[27] = mcposalmon2;
  mccomp_posr[27] = mcposrlmon2;
  mcNCounter[27]  = mcPCounter[27] = mcP2Counter[27] = 0;
  mcAbsorbProp[27]= 0;
    /* Component psd4. */
  /* Setting parameters for component psd4. */
  SIG_MESSAGE("psd4 (Init:SetPar)");
#line 214 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccpsd4_nx = 100;
#line 214 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccpsd4_ny = 100;
#line 214 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  if("psd4.dat") strncpy(mccpsd4_filename, "psd4.dat" ? "psd4.dat" : "", 16384); else mccpsd4_filename[0]='\0';
#line 214 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccpsd4_xmin = -0.0075;
#line 215 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccpsd4_xmax = 0.0075;
#line 215 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccpsd4_ymin = -0.0075;
#line 215 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccpsd4_ymax = 0.0075;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccpsd4_xwidth = 0;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccpsd4_yheight = 0;
#line 51 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
  mccpsd4_restore_neutron = 0;
#line 12590 "./ISIS_SANS2d.c"

  SIG_MESSAGE("psd4 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 12597 "./ISIS_SANS2d.c"
  rot_mul(mctr1, mcrotaAPERTURE2, mcrotapsd4);
  rot_transpose(mcrotalmon2, mctr1);
  rot_mul(mcrotapsd4, mctr1, mcrotrpsd4);
  mctc1 = coords_set(
#line 216 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    0.0,
#line 216 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    0.0,
#line 216 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/ISIS_SANS2d/ISIS_SANS2d.instr"
    0.286);
#line 12608 "./ISIS_SANS2d.c"
  rot_transpose(mcrotaAPERTURE2, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposapsd4 = coords_add(mcposaAPERTURE2, mctc2);
  mctc1 = coords_sub(mcposalmon2, mcposapsd4);
  mcposrpsd4 = rot_apply(mcrotapsd4, mctc1);
  mcDEBUG_COMPONENT("psd4", mcposapsd4, mcrotapsd4)
  mccomp_posa[28] = mcposapsd4;
  mccomp_posr[28] = mcposrpsd4;
  mcNCounter[28]  = mcPCounter[28] = mcP2Counter[28] = 0;
  mcAbsorbProp[28]= 0;
  /* Component initializations. */
  /* Initializations for component Origin. */
  SIG_MESSAGE("Origin (Init)");

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
#line 836 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../contrib/ISIS_moderator.comp"
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
#line 12878 "./ISIS_SANS2d.c"
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

  /* Initializations for component lmon1. */
  SIG_MESSAGE("lmon1 (Init)");
#define mccompcurname  lmon1
#define mccompcurtype  L_monitor
#define mccompcurindex 3
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
#line 12951 "./ISIS_SANS2d.c"
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
#define mccompcurindex 4
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
#line 13014 "./ISIS_SANS2d.c"
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

  /* Initializations for component bender1. */
  SIG_MESSAGE("bender1 (Init)");
#define mccompcurname  bender1
#define mccompcurtype  Guide_gravity
#define mccompcurindex 5
#define GVars mccbender1_GVars
#define pTable mccbender1_pTable
#define w1 mccbender1_w1
#define h1 mccbender1_h1
#define w2 mccbender1_w2
#define h2 mccbender1_h2
#define l mccbender1_l
#define R0 mccbender1_R0
#define Qc mccbender1_Qc
#define alpha mccbender1_alpha
#define m mccbender1_m
#define W mccbender1_W
#define nslit mccbender1_nslit
#define d mccbender1_d
#define mleft mccbender1_mleft
#define mright mccbender1_mright
#define mtop mccbender1_mtop
#define mbottom mccbender1_mbottom
#define nhslit mccbender1_nhslit
#define G mccbender1_G
#define aleft mccbender1_aleft
#define aright mccbender1_aright
#define atop mccbender1_atop
#define abottom mccbender1_abottom
#define wavy mccbender1_wavy
#define wavy_z mccbender1_wavy_z
#define wavy_tb mccbender1_wavy_tb
#define wavy_lr mccbender1_wavy_lr
#define chamfers mccbender1_chamfers
#define chamfers_z mccbender1_chamfers_z
#define chamfers_lr mccbender1_chamfers_lr
#define chamfers_tb mccbender1_chamfers_tb
#define nelements mccbender1_nelements
#define nu mccbender1_nu
#define phase mccbender1_phase
#define reflect mccbender1_reflect
#line 339 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide_gravity.comp"
{
  double Gx=0, Gy=-GRAVITY, Gz=0;
  Coords mcLocG;
  int i;

  if (reflect && strlen(reflect) && strcmp(reflect,"NULL") && strcmp(reflect,"0")) {
    if (Table_Read(&pTable, reflect, 1) <= 0) /* read 1st block data from file into pTable */
      exit(fprintf(stderr,"Guide_gravity: %s: can not read file %s\n", NAME_CURRENT_COMP, reflect));
  } else {
    if (W < 0 || R0 < 0 || Qc < 0)
    { fprintf(stderr,"Guide_gravity: %s: W R0 Qc must be >0.\n", NAME_CURRENT_COMP);
      exit(-1); }
  }

  if (nslit <= 0 || nhslit <= 0)
  { fprintf(stderr,"Guide_gravity: %s: nslit nhslit must be >0.\n", NAME_CURRENT_COMP);
    exit(-1); }

  if (!w1 || !h1)
  { fprintf(stderr,"Guide_gravity: %s: input window is closed (w1=h1=0).\n", NAME_CURRENT_COMP);
    exit(-1); }

  if (d*nslit > w1) exit(fprintf(stderr, "Guide_gravity: %s: absorbing walls fill input window. No space left for transmission (d*nslit > w1).\n", NAME_CURRENT_COMP));

  if (!w2) w2=w1;
  if (!h2) h2=h1;

  if (mcgravitation) G=-GRAVITY;
  mcLocG = rot_apply(ROT_A_CURRENT_COMP, coords_set(0,G,0));
  coords_get(mcLocG, &Gx, &Gy, &Gz);

  strcpy(GVars.compcurname, NAME_CURRENT_COMP);

  if (l > 0 && nelements > 0) {

    Gravity_guide_Init(&GVars,
      w1, h1, w2, h2, l, R0,
      Qc, alpha, m, W, nslit, d,
      Gx, Gy, Gz, mleft, mright, mtop,
      mbottom, nhslit, wavy_lr, wavy_tb, wavy_z, wavy,
      chamfers_z, chamfers_lr, chamfers_tb, chamfers,nu,phase,aleft,aright,atop,abottom);
    if (!G) for (i=0; i<5; GVars.A[i++] = 0);
    if (GVars.fc_freq != 0 || GVars.fc_phase != 0) {
      if (w1 != w2 || h1 != h2)
      exit(fprintf(stderr,"Guide_gravity: %s: rotating slit pack must be straight (w1=w2 and h1=h2).\n", NAME_CURRENT_COMP));
      printf("Guide_gravity: %s: Fermi Chopper mode: frequency=%g [Hz] phase=%g [deg]\n",
        NAME_CURRENT_COMP, GVars.fc_freq, GVars.fc_phase);
    }
  } else printf("Guide_gravity: %s: unactivated (l=0 or nelements=0)\n", NAME_CURRENT_COMP);

}
#line 13125 "./ISIS_SANS2d.c"
#undef reflect
#undef phase
#undef nu
#undef nelements
#undef chamfers_tb
#undef chamfers_lr
#undef chamfers_z
#undef chamfers
#undef wavy_lr
#undef wavy_tb
#undef wavy_z
#undef wavy
#undef abottom
#undef atop
#undef aright
#undef aleft
#undef G
#undef nhslit
#undef mbottom
#undef mtop
#undef mright
#undef mleft
#undef d
#undef nslit
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
#undef pTable
#undef GVars
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component bender2. */
  SIG_MESSAGE("bender2 (Init)");
#define mccompcurname  bender2
#define mccompcurtype  Guide_gravity
#define mccompcurindex 6
#define GVars mccbender2_GVars
#define pTable mccbender2_pTable
#define w1 mccbender2_w1
#define h1 mccbender2_h1
#define w2 mccbender2_w2
#define h2 mccbender2_h2
#define l mccbender2_l
#define R0 mccbender2_R0
#define Qc mccbender2_Qc
#define alpha mccbender2_alpha
#define m mccbender2_m
#define W mccbender2_W
#define nslit mccbender2_nslit
#define d mccbender2_d
#define mleft mccbender2_mleft
#define mright mccbender2_mright
#define mtop mccbender2_mtop
#define mbottom mccbender2_mbottom
#define nhslit mccbender2_nhslit
#define G mccbender2_G
#define aleft mccbender2_aleft
#define aright mccbender2_aright
#define atop mccbender2_atop
#define abottom mccbender2_abottom
#define wavy mccbender2_wavy
#define wavy_z mccbender2_wavy_z
#define wavy_tb mccbender2_wavy_tb
#define wavy_lr mccbender2_wavy_lr
#define chamfers mccbender2_chamfers
#define chamfers_z mccbender2_chamfers_z
#define chamfers_lr mccbender2_chamfers_lr
#define chamfers_tb mccbender2_chamfers_tb
#define nelements mccbender2_nelements
#define nu mccbender2_nu
#define phase mccbender2_phase
#define reflect mccbender2_reflect
#line 339 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide_gravity.comp"
{
  double Gx=0, Gy=-GRAVITY, Gz=0;
  Coords mcLocG;
  int i;

  if (reflect && strlen(reflect) && strcmp(reflect,"NULL") && strcmp(reflect,"0")) {
    if (Table_Read(&pTable, reflect, 1) <= 0) /* read 1st block data from file into pTable */
      exit(fprintf(stderr,"Guide_gravity: %s: can not read file %s\n", NAME_CURRENT_COMP, reflect));
  } else {
    if (W < 0 || R0 < 0 || Qc < 0)
    { fprintf(stderr,"Guide_gravity: %s: W R0 Qc must be >0.\n", NAME_CURRENT_COMP);
      exit(-1); }
  }

  if (nslit <= 0 || nhslit <= 0)
  { fprintf(stderr,"Guide_gravity: %s: nslit nhslit must be >0.\n", NAME_CURRENT_COMP);
    exit(-1); }

  if (!w1 || !h1)
  { fprintf(stderr,"Guide_gravity: %s: input window is closed (w1=h1=0).\n", NAME_CURRENT_COMP);
    exit(-1); }

  if (d*nslit > w1) exit(fprintf(stderr, "Guide_gravity: %s: absorbing walls fill input window. No space left for transmission (d*nslit > w1).\n", NAME_CURRENT_COMP));

  if (!w2) w2=w1;
  if (!h2) h2=h1;

  if (mcgravitation) G=-GRAVITY;
  mcLocG = rot_apply(ROT_A_CURRENT_COMP, coords_set(0,G,0));
  coords_get(mcLocG, &Gx, &Gy, &Gz);

  strcpy(GVars.compcurname, NAME_CURRENT_COMP);

  if (l > 0 && nelements > 0) {

    Gravity_guide_Init(&GVars,
      w1, h1, w2, h2, l, R0,
      Qc, alpha, m, W, nslit, d,
      Gx, Gy, Gz, mleft, mright, mtop,
      mbottom, nhslit, wavy_lr, wavy_tb, wavy_z, wavy,
      chamfers_z, chamfers_lr, chamfers_tb, chamfers,nu,phase,aleft,aright,atop,abottom);
    if (!G) for (i=0; i<5; GVars.A[i++] = 0);
    if (GVars.fc_freq != 0 || GVars.fc_phase != 0) {
      if (w1 != w2 || h1 != h2)
      exit(fprintf(stderr,"Guide_gravity: %s: rotating slit pack must be straight (w1=w2 and h1=h2).\n", NAME_CURRENT_COMP));
      printf("Guide_gravity: %s: Fermi Chopper mode: frequency=%g [Hz] phase=%g [deg]\n",
        NAME_CURRENT_COMP, GVars.fc_freq, GVars.fc_phase);
    }
  } else printf("Guide_gravity: %s: unactivated (l=0 or nelements=0)\n", NAME_CURRENT_COMP);

}
#line 13259 "./ISIS_SANS2d.c"
#undef reflect
#undef phase
#undef nu
#undef nelements
#undef chamfers_tb
#undef chamfers_lr
#undef chamfers_z
#undef chamfers
#undef wavy_lr
#undef wavy_tb
#undef wavy_z
#undef wavy
#undef abottom
#undef atop
#undef aright
#undef aleft
#undef G
#undef nhslit
#undef mbottom
#undef mtop
#undef mright
#undef mleft
#undef d
#undef nslit
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
#undef pTable
#undef GVars
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component bender3. */
  SIG_MESSAGE("bender3 (Init)");
#define mccompcurname  bender3
#define mccompcurtype  Guide_gravity
#define mccompcurindex 7
#define GVars mccbender3_GVars
#define pTable mccbender3_pTable
#define w1 mccbender3_w1
#define h1 mccbender3_h1
#define w2 mccbender3_w2
#define h2 mccbender3_h2
#define l mccbender3_l
#define R0 mccbender3_R0
#define Qc mccbender3_Qc
#define alpha mccbender3_alpha
#define m mccbender3_m
#define W mccbender3_W
#define nslit mccbender3_nslit
#define d mccbender3_d
#define mleft mccbender3_mleft
#define mright mccbender3_mright
#define mtop mccbender3_mtop
#define mbottom mccbender3_mbottom
#define nhslit mccbender3_nhslit
#define G mccbender3_G
#define aleft mccbender3_aleft
#define aright mccbender3_aright
#define atop mccbender3_atop
#define abottom mccbender3_abottom
#define wavy mccbender3_wavy
#define wavy_z mccbender3_wavy_z
#define wavy_tb mccbender3_wavy_tb
#define wavy_lr mccbender3_wavy_lr
#define chamfers mccbender3_chamfers
#define chamfers_z mccbender3_chamfers_z
#define chamfers_lr mccbender3_chamfers_lr
#define chamfers_tb mccbender3_chamfers_tb
#define nelements mccbender3_nelements
#define nu mccbender3_nu
#define phase mccbender3_phase
#define reflect mccbender3_reflect
#line 339 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide_gravity.comp"
{
  double Gx=0, Gy=-GRAVITY, Gz=0;
  Coords mcLocG;
  int i;

  if (reflect && strlen(reflect) && strcmp(reflect,"NULL") && strcmp(reflect,"0")) {
    if (Table_Read(&pTable, reflect, 1) <= 0) /* read 1st block data from file into pTable */
      exit(fprintf(stderr,"Guide_gravity: %s: can not read file %s\n", NAME_CURRENT_COMP, reflect));
  } else {
    if (W < 0 || R0 < 0 || Qc < 0)
    { fprintf(stderr,"Guide_gravity: %s: W R0 Qc must be >0.\n", NAME_CURRENT_COMP);
      exit(-1); }
  }

  if (nslit <= 0 || nhslit <= 0)
  { fprintf(stderr,"Guide_gravity: %s: nslit nhslit must be >0.\n", NAME_CURRENT_COMP);
    exit(-1); }

  if (!w1 || !h1)
  { fprintf(stderr,"Guide_gravity: %s: input window is closed (w1=h1=0).\n", NAME_CURRENT_COMP);
    exit(-1); }

  if (d*nslit > w1) exit(fprintf(stderr, "Guide_gravity: %s: absorbing walls fill input window. No space left for transmission (d*nslit > w1).\n", NAME_CURRENT_COMP));

  if (!w2) w2=w1;
  if (!h2) h2=h1;

  if (mcgravitation) G=-GRAVITY;
  mcLocG = rot_apply(ROT_A_CURRENT_COMP, coords_set(0,G,0));
  coords_get(mcLocG, &Gx, &Gy, &Gz);

  strcpy(GVars.compcurname, NAME_CURRENT_COMP);

  if (l > 0 && nelements > 0) {

    Gravity_guide_Init(&GVars,
      w1, h1, w2, h2, l, R0,
      Qc, alpha, m, W, nslit, d,
      Gx, Gy, Gz, mleft, mright, mtop,
      mbottom, nhslit, wavy_lr, wavy_tb, wavy_z, wavy,
      chamfers_z, chamfers_lr, chamfers_tb, chamfers,nu,phase,aleft,aright,atop,abottom);
    if (!G) for (i=0; i<5; GVars.A[i++] = 0);
    if (GVars.fc_freq != 0 || GVars.fc_phase != 0) {
      if (w1 != w2 || h1 != h2)
      exit(fprintf(stderr,"Guide_gravity: %s: rotating slit pack must be straight (w1=w2 and h1=h2).\n", NAME_CURRENT_COMP));
      printf("Guide_gravity: %s: Fermi Chopper mode: frequency=%g [Hz] phase=%g [deg]\n",
        NAME_CURRENT_COMP, GVars.fc_freq, GVars.fc_phase);
    }
  } else printf("Guide_gravity: %s: unactivated (l=0 or nelements=0)\n", NAME_CURRENT_COMP);

}
#line 13393 "./ISIS_SANS2d.c"
#undef reflect
#undef phase
#undef nu
#undef nelements
#undef chamfers_tb
#undef chamfers_lr
#undef chamfers_z
#undef chamfers
#undef wavy_lr
#undef wavy_tb
#undef wavy_z
#undef wavy
#undef abottom
#undef atop
#undef aright
#undef aleft
#undef G
#undef nhslit
#undef mbottom
#undef mtop
#undef mright
#undef mleft
#undef d
#undef nslit
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
#undef pTable
#undef GVars
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component bender4. */
  SIG_MESSAGE("bender4 (Init)");
#define mccompcurname  bender4
#define mccompcurtype  Guide_gravity
#define mccompcurindex 8
#define GVars mccbender4_GVars
#define pTable mccbender4_pTable
#define w1 mccbender4_w1
#define h1 mccbender4_h1
#define w2 mccbender4_w2
#define h2 mccbender4_h2
#define l mccbender4_l
#define R0 mccbender4_R0
#define Qc mccbender4_Qc
#define alpha mccbender4_alpha
#define m mccbender4_m
#define W mccbender4_W
#define nslit mccbender4_nslit
#define d mccbender4_d
#define mleft mccbender4_mleft
#define mright mccbender4_mright
#define mtop mccbender4_mtop
#define mbottom mccbender4_mbottom
#define nhslit mccbender4_nhslit
#define G mccbender4_G
#define aleft mccbender4_aleft
#define aright mccbender4_aright
#define atop mccbender4_atop
#define abottom mccbender4_abottom
#define wavy mccbender4_wavy
#define wavy_z mccbender4_wavy_z
#define wavy_tb mccbender4_wavy_tb
#define wavy_lr mccbender4_wavy_lr
#define chamfers mccbender4_chamfers
#define chamfers_z mccbender4_chamfers_z
#define chamfers_lr mccbender4_chamfers_lr
#define chamfers_tb mccbender4_chamfers_tb
#define nelements mccbender4_nelements
#define nu mccbender4_nu
#define phase mccbender4_phase
#define reflect mccbender4_reflect
#line 339 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide_gravity.comp"
{
  double Gx=0, Gy=-GRAVITY, Gz=0;
  Coords mcLocG;
  int i;

  if (reflect && strlen(reflect) && strcmp(reflect,"NULL") && strcmp(reflect,"0")) {
    if (Table_Read(&pTable, reflect, 1) <= 0) /* read 1st block data from file into pTable */
      exit(fprintf(stderr,"Guide_gravity: %s: can not read file %s\n", NAME_CURRENT_COMP, reflect));
  } else {
    if (W < 0 || R0 < 0 || Qc < 0)
    { fprintf(stderr,"Guide_gravity: %s: W R0 Qc must be >0.\n", NAME_CURRENT_COMP);
      exit(-1); }
  }

  if (nslit <= 0 || nhslit <= 0)
  { fprintf(stderr,"Guide_gravity: %s: nslit nhslit must be >0.\n", NAME_CURRENT_COMP);
    exit(-1); }

  if (!w1 || !h1)
  { fprintf(stderr,"Guide_gravity: %s: input window is closed (w1=h1=0).\n", NAME_CURRENT_COMP);
    exit(-1); }

  if (d*nslit > w1) exit(fprintf(stderr, "Guide_gravity: %s: absorbing walls fill input window. No space left for transmission (d*nslit > w1).\n", NAME_CURRENT_COMP));

  if (!w2) w2=w1;
  if (!h2) h2=h1;

  if (mcgravitation) G=-GRAVITY;
  mcLocG = rot_apply(ROT_A_CURRENT_COMP, coords_set(0,G,0));
  coords_get(mcLocG, &Gx, &Gy, &Gz);

  strcpy(GVars.compcurname, NAME_CURRENT_COMP);

  if (l > 0 && nelements > 0) {

    Gravity_guide_Init(&GVars,
      w1, h1, w2, h2, l, R0,
      Qc, alpha, m, W, nslit, d,
      Gx, Gy, Gz, mleft, mright, mtop,
      mbottom, nhslit, wavy_lr, wavy_tb, wavy_z, wavy,
      chamfers_z, chamfers_lr, chamfers_tb, chamfers,nu,phase,aleft,aright,atop,abottom);
    if (!G) for (i=0; i<5; GVars.A[i++] = 0);
    if (GVars.fc_freq != 0 || GVars.fc_phase != 0) {
      if (w1 != w2 || h1 != h2)
      exit(fprintf(stderr,"Guide_gravity: %s: rotating slit pack must be straight (w1=w2 and h1=h2).\n", NAME_CURRENT_COMP));
      printf("Guide_gravity: %s: Fermi Chopper mode: frequency=%g [Hz] phase=%g [deg]\n",
        NAME_CURRENT_COMP, GVars.fc_freq, GVars.fc_phase);
    }
  } else printf("Guide_gravity: %s: unactivated (l=0 or nelements=0)\n", NAME_CURRENT_COMP);

}
#line 13527 "./ISIS_SANS2d.c"
#undef reflect
#undef phase
#undef nu
#undef nelements
#undef chamfers_tb
#undef chamfers_lr
#undef chamfers_z
#undef chamfers
#undef wavy_lr
#undef wavy_tb
#undef wavy_z
#undef wavy
#undef abottom
#undef atop
#undef aright
#undef aleft
#undef G
#undef nhslit
#undef mbottom
#undef mtop
#undef mright
#undef mleft
#undef d
#undef nslit
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
#undef pTable
#undef GVars
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component bender5. */
  SIG_MESSAGE("bender5 (Init)");
#define mccompcurname  bender5
#define mccompcurtype  Guide_gravity
#define mccompcurindex 9
#define GVars mccbender5_GVars
#define pTable mccbender5_pTable
#define w1 mccbender5_w1
#define h1 mccbender5_h1
#define w2 mccbender5_w2
#define h2 mccbender5_h2
#define l mccbender5_l
#define R0 mccbender5_R0
#define Qc mccbender5_Qc
#define alpha mccbender5_alpha
#define m mccbender5_m
#define W mccbender5_W
#define nslit mccbender5_nslit
#define d mccbender5_d
#define mleft mccbender5_mleft
#define mright mccbender5_mright
#define mtop mccbender5_mtop
#define mbottom mccbender5_mbottom
#define nhslit mccbender5_nhslit
#define G mccbender5_G
#define aleft mccbender5_aleft
#define aright mccbender5_aright
#define atop mccbender5_atop
#define abottom mccbender5_abottom
#define wavy mccbender5_wavy
#define wavy_z mccbender5_wavy_z
#define wavy_tb mccbender5_wavy_tb
#define wavy_lr mccbender5_wavy_lr
#define chamfers mccbender5_chamfers
#define chamfers_z mccbender5_chamfers_z
#define chamfers_lr mccbender5_chamfers_lr
#define chamfers_tb mccbender5_chamfers_tb
#define nelements mccbender5_nelements
#define nu mccbender5_nu
#define phase mccbender5_phase
#define reflect mccbender5_reflect
#line 339 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide_gravity.comp"
{
  double Gx=0, Gy=-GRAVITY, Gz=0;
  Coords mcLocG;
  int i;

  if (reflect && strlen(reflect) && strcmp(reflect,"NULL") && strcmp(reflect,"0")) {
    if (Table_Read(&pTable, reflect, 1) <= 0) /* read 1st block data from file into pTable */
      exit(fprintf(stderr,"Guide_gravity: %s: can not read file %s\n", NAME_CURRENT_COMP, reflect));
  } else {
    if (W < 0 || R0 < 0 || Qc < 0)
    { fprintf(stderr,"Guide_gravity: %s: W R0 Qc must be >0.\n", NAME_CURRENT_COMP);
      exit(-1); }
  }

  if (nslit <= 0 || nhslit <= 0)
  { fprintf(stderr,"Guide_gravity: %s: nslit nhslit must be >0.\n", NAME_CURRENT_COMP);
    exit(-1); }

  if (!w1 || !h1)
  { fprintf(stderr,"Guide_gravity: %s: input window is closed (w1=h1=0).\n", NAME_CURRENT_COMP);
    exit(-1); }

  if (d*nslit > w1) exit(fprintf(stderr, "Guide_gravity: %s: absorbing walls fill input window. No space left for transmission (d*nslit > w1).\n", NAME_CURRENT_COMP));

  if (!w2) w2=w1;
  if (!h2) h2=h1;

  if (mcgravitation) G=-GRAVITY;
  mcLocG = rot_apply(ROT_A_CURRENT_COMP, coords_set(0,G,0));
  coords_get(mcLocG, &Gx, &Gy, &Gz);

  strcpy(GVars.compcurname, NAME_CURRENT_COMP);

  if (l > 0 && nelements > 0) {

    Gravity_guide_Init(&GVars,
      w1, h1, w2, h2, l, R0,
      Qc, alpha, m, W, nslit, d,
      Gx, Gy, Gz, mleft, mright, mtop,
      mbottom, nhslit, wavy_lr, wavy_tb, wavy_z, wavy,
      chamfers_z, chamfers_lr, chamfers_tb, chamfers,nu,phase,aleft,aright,atop,abottom);
    if (!G) for (i=0; i<5; GVars.A[i++] = 0);
    if (GVars.fc_freq != 0 || GVars.fc_phase != 0) {
      if (w1 != w2 || h1 != h2)
      exit(fprintf(stderr,"Guide_gravity: %s: rotating slit pack must be straight (w1=w2 and h1=h2).\n", NAME_CURRENT_COMP));
      printf("Guide_gravity: %s: Fermi Chopper mode: frequency=%g [Hz] phase=%g [deg]\n",
        NAME_CURRENT_COMP, GVars.fc_freq, GVars.fc_phase);
    }
  } else printf("Guide_gravity: %s: unactivated (l=0 or nelements=0)\n", NAME_CURRENT_COMP);

}
#line 13661 "./ISIS_SANS2d.c"
#undef reflect
#undef phase
#undef nu
#undef nelements
#undef chamfers_tb
#undef chamfers_lr
#undef chamfers_z
#undef chamfers
#undef wavy_lr
#undef wavy_tb
#undef wavy_z
#undef wavy
#undef abottom
#undef atop
#undef aright
#undef aleft
#undef G
#undef nhslit
#undef mbottom
#undef mtop
#undef mright
#undef mleft
#undef d
#undef nslit
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
#undef pTable
#undef GVars
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component bender6. */
  SIG_MESSAGE("bender6 (Init)");
#define mccompcurname  bender6
#define mccompcurtype  Guide_gravity
#define mccompcurindex 10
#define GVars mccbender6_GVars
#define pTable mccbender6_pTable
#define w1 mccbender6_w1
#define h1 mccbender6_h1
#define w2 mccbender6_w2
#define h2 mccbender6_h2
#define l mccbender6_l
#define R0 mccbender6_R0
#define Qc mccbender6_Qc
#define alpha mccbender6_alpha
#define m mccbender6_m
#define W mccbender6_W
#define nslit mccbender6_nslit
#define d mccbender6_d
#define mleft mccbender6_mleft
#define mright mccbender6_mright
#define mtop mccbender6_mtop
#define mbottom mccbender6_mbottom
#define nhslit mccbender6_nhslit
#define G mccbender6_G
#define aleft mccbender6_aleft
#define aright mccbender6_aright
#define atop mccbender6_atop
#define abottom mccbender6_abottom
#define wavy mccbender6_wavy
#define wavy_z mccbender6_wavy_z
#define wavy_tb mccbender6_wavy_tb
#define wavy_lr mccbender6_wavy_lr
#define chamfers mccbender6_chamfers
#define chamfers_z mccbender6_chamfers_z
#define chamfers_lr mccbender6_chamfers_lr
#define chamfers_tb mccbender6_chamfers_tb
#define nelements mccbender6_nelements
#define nu mccbender6_nu
#define phase mccbender6_phase
#define reflect mccbender6_reflect
#line 339 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide_gravity.comp"
{
  double Gx=0, Gy=-GRAVITY, Gz=0;
  Coords mcLocG;
  int i;

  if (reflect && strlen(reflect) && strcmp(reflect,"NULL") && strcmp(reflect,"0")) {
    if (Table_Read(&pTable, reflect, 1) <= 0) /* read 1st block data from file into pTable */
      exit(fprintf(stderr,"Guide_gravity: %s: can not read file %s\n", NAME_CURRENT_COMP, reflect));
  } else {
    if (W < 0 || R0 < 0 || Qc < 0)
    { fprintf(stderr,"Guide_gravity: %s: W R0 Qc must be >0.\n", NAME_CURRENT_COMP);
      exit(-1); }
  }

  if (nslit <= 0 || nhslit <= 0)
  { fprintf(stderr,"Guide_gravity: %s: nslit nhslit must be >0.\n", NAME_CURRENT_COMP);
    exit(-1); }

  if (!w1 || !h1)
  { fprintf(stderr,"Guide_gravity: %s: input window is closed (w1=h1=0).\n", NAME_CURRENT_COMP);
    exit(-1); }

  if (d*nslit > w1) exit(fprintf(stderr, "Guide_gravity: %s: absorbing walls fill input window. No space left for transmission (d*nslit > w1).\n", NAME_CURRENT_COMP));

  if (!w2) w2=w1;
  if (!h2) h2=h1;

  if (mcgravitation) G=-GRAVITY;
  mcLocG = rot_apply(ROT_A_CURRENT_COMP, coords_set(0,G,0));
  coords_get(mcLocG, &Gx, &Gy, &Gz);

  strcpy(GVars.compcurname, NAME_CURRENT_COMP);

  if (l > 0 && nelements > 0) {

    Gravity_guide_Init(&GVars,
      w1, h1, w2, h2, l, R0,
      Qc, alpha, m, W, nslit, d,
      Gx, Gy, Gz, mleft, mright, mtop,
      mbottom, nhslit, wavy_lr, wavy_tb, wavy_z, wavy,
      chamfers_z, chamfers_lr, chamfers_tb, chamfers,nu,phase,aleft,aright,atop,abottom);
    if (!G) for (i=0; i<5; GVars.A[i++] = 0);
    if (GVars.fc_freq != 0 || GVars.fc_phase != 0) {
      if (w1 != w2 || h1 != h2)
      exit(fprintf(stderr,"Guide_gravity: %s: rotating slit pack must be straight (w1=w2 and h1=h2).\n", NAME_CURRENT_COMP));
      printf("Guide_gravity: %s: Fermi Chopper mode: frequency=%g [Hz] phase=%g [deg]\n",
        NAME_CURRENT_COMP, GVars.fc_freq, GVars.fc_phase);
    }
  } else printf("Guide_gravity: %s: unactivated (l=0 or nelements=0)\n", NAME_CURRENT_COMP);

}
#line 13795 "./ISIS_SANS2d.c"
#undef reflect
#undef phase
#undef nu
#undef nelements
#undef chamfers_tb
#undef chamfers_lr
#undef chamfers_z
#undef chamfers
#undef wavy_lr
#undef wavy_tb
#undef wavy_z
#undef wavy
#undef abottom
#undef atop
#undef aright
#undef aleft
#undef G
#undef nhslit
#undef mbottom
#undef mtop
#undef mright
#undef mleft
#undef d
#undef nslit
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
#undef pTable
#undef GVars
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component bender7. */
  SIG_MESSAGE("bender7 (Init)");
#define mccompcurname  bender7
#define mccompcurtype  Guide_gravity
#define mccompcurindex 11
#define GVars mccbender7_GVars
#define pTable mccbender7_pTable
#define w1 mccbender7_w1
#define h1 mccbender7_h1
#define w2 mccbender7_w2
#define h2 mccbender7_h2
#define l mccbender7_l
#define R0 mccbender7_R0
#define Qc mccbender7_Qc
#define alpha mccbender7_alpha
#define m mccbender7_m
#define W mccbender7_W
#define nslit mccbender7_nslit
#define d mccbender7_d
#define mleft mccbender7_mleft
#define mright mccbender7_mright
#define mtop mccbender7_mtop
#define mbottom mccbender7_mbottom
#define nhslit mccbender7_nhslit
#define G mccbender7_G
#define aleft mccbender7_aleft
#define aright mccbender7_aright
#define atop mccbender7_atop
#define abottom mccbender7_abottom
#define wavy mccbender7_wavy
#define wavy_z mccbender7_wavy_z
#define wavy_tb mccbender7_wavy_tb
#define wavy_lr mccbender7_wavy_lr
#define chamfers mccbender7_chamfers
#define chamfers_z mccbender7_chamfers_z
#define chamfers_lr mccbender7_chamfers_lr
#define chamfers_tb mccbender7_chamfers_tb
#define nelements mccbender7_nelements
#define nu mccbender7_nu
#define phase mccbender7_phase
#define reflect mccbender7_reflect
#line 339 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide_gravity.comp"
{
  double Gx=0, Gy=-GRAVITY, Gz=0;
  Coords mcLocG;
  int i;

  if (reflect && strlen(reflect) && strcmp(reflect,"NULL") && strcmp(reflect,"0")) {
    if (Table_Read(&pTable, reflect, 1) <= 0) /* read 1st block data from file into pTable */
      exit(fprintf(stderr,"Guide_gravity: %s: can not read file %s\n", NAME_CURRENT_COMP, reflect));
  } else {
    if (W < 0 || R0 < 0 || Qc < 0)
    { fprintf(stderr,"Guide_gravity: %s: W R0 Qc must be >0.\n", NAME_CURRENT_COMP);
      exit(-1); }
  }

  if (nslit <= 0 || nhslit <= 0)
  { fprintf(stderr,"Guide_gravity: %s: nslit nhslit must be >0.\n", NAME_CURRENT_COMP);
    exit(-1); }

  if (!w1 || !h1)
  { fprintf(stderr,"Guide_gravity: %s: input window is closed (w1=h1=0).\n", NAME_CURRENT_COMP);
    exit(-1); }

  if (d*nslit > w1) exit(fprintf(stderr, "Guide_gravity: %s: absorbing walls fill input window. No space left for transmission (d*nslit > w1).\n", NAME_CURRENT_COMP));

  if (!w2) w2=w1;
  if (!h2) h2=h1;

  if (mcgravitation) G=-GRAVITY;
  mcLocG = rot_apply(ROT_A_CURRENT_COMP, coords_set(0,G,0));
  coords_get(mcLocG, &Gx, &Gy, &Gz);

  strcpy(GVars.compcurname, NAME_CURRENT_COMP);

  if (l > 0 && nelements > 0) {

    Gravity_guide_Init(&GVars,
      w1, h1, w2, h2, l, R0,
      Qc, alpha, m, W, nslit, d,
      Gx, Gy, Gz, mleft, mright, mtop,
      mbottom, nhslit, wavy_lr, wavy_tb, wavy_z, wavy,
      chamfers_z, chamfers_lr, chamfers_tb, chamfers,nu,phase,aleft,aright,atop,abottom);
    if (!G) for (i=0; i<5; GVars.A[i++] = 0);
    if (GVars.fc_freq != 0 || GVars.fc_phase != 0) {
      if (w1 != w2 || h1 != h2)
      exit(fprintf(stderr,"Guide_gravity: %s: rotating slit pack must be straight (w1=w2 and h1=h2).\n", NAME_CURRENT_COMP));
      printf("Guide_gravity: %s: Fermi Chopper mode: frequency=%g [Hz] phase=%g [deg]\n",
        NAME_CURRENT_COMP, GVars.fc_freq, GVars.fc_phase);
    }
  } else printf("Guide_gravity: %s: unactivated (l=0 or nelements=0)\n", NAME_CURRENT_COMP);

}
#line 13929 "./ISIS_SANS2d.c"
#undef reflect
#undef phase
#undef nu
#undef nelements
#undef chamfers_tb
#undef chamfers_lr
#undef chamfers_z
#undef chamfers
#undef wavy_lr
#undef wavy_tb
#undef wavy_z
#undef wavy
#undef abottom
#undef atop
#undef aright
#undef aleft
#undef G
#undef nhslit
#undef mbottom
#undef mtop
#undef mright
#undef mleft
#undef d
#undef nslit
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
#undef pTable
#undef GVars
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component bender8. */
  SIG_MESSAGE("bender8 (Init)");
#define mccompcurname  bender8
#define mccompcurtype  Guide_gravity
#define mccompcurindex 12
#define GVars mccbender8_GVars
#define pTable mccbender8_pTable
#define w1 mccbender8_w1
#define h1 mccbender8_h1
#define w2 mccbender8_w2
#define h2 mccbender8_h2
#define l mccbender8_l
#define R0 mccbender8_R0
#define Qc mccbender8_Qc
#define alpha mccbender8_alpha
#define m mccbender8_m
#define W mccbender8_W
#define nslit mccbender8_nslit
#define d mccbender8_d
#define mleft mccbender8_mleft
#define mright mccbender8_mright
#define mtop mccbender8_mtop
#define mbottom mccbender8_mbottom
#define nhslit mccbender8_nhslit
#define G mccbender8_G
#define aleft mccbender8_aleft
#define aright mccbender8_aright
#define atop mccbender8_atop
#define abottom mccbender8_abottom
#define wavy mccbender8_wavy
#define wavy_z mccbender8_wavy_z
#define wavy_tb mccbender8_wavy_tb
#define wavy_lr mccbender8_wavy_lr
#define chamfers mccbender8_chamfers
#define chamfers_z mccbender8_chamfers_z
#define chamfers_lr mccbender8_chamfers_lr
#define chamfers_tb mccbender8_chamfers_tb
#define nelements mccbender8_nelements
#define nu mccbender8_nu
#define phase mccbender8_phase
#define reflect mccbender8_reflect
#line 339 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide_gravity.comp"
{
  double Gx=0, Gy=-GRAVITY, Gz=0;
  Coords mcLocG;
  int i;

  if (reflect && strlen(reflect) && strcmp(reflect,"NULL") && strcmp(reflect,"0")) {
    if (Table_Read(&pTable, reflect, 1) <= 0) /* read 1st block data from file into pTable */
      exit(fprintf(stderr,"Guide_gravity: %s: can not read file %s\n", NAME_CURRENT_COMP, reflect));
  } else {
    if (W < 0 || R0 < 0 || Qc < 0)
    { fprintf(stderr,"Guide_gravity: %s: W R0 Qc must be >0.\n", NAME_CURRENT_COMP);
      exit(-1); }
  }

  if (nslit <= 0 || nhslit <= 0)
  { fprintf(stderr,"Guide_gravity: %s: nslit nhslit must be >0.\n", NAME_CURRENT_COMP);
    exit(-1); }

  if (!w1 || !h1)
  { fprintf(stderr,"Guide_gravity: %s: input window is closed (w1=h1=0).\n", NAME_CURRENT_COMP);
    exit(-1); }

  if (d*nslit > w1) exit(fprintf(stderr, "Guide_gravity: %s: absorbing walls fill input window. No space left for transmission (d*nslit > w1).\n", NAME_CURRENT_COMP));

  if (!w2) w2=w1;
  if (!h2) h2=h1;

  if (mcgravitation) G=-GRAVITY;
  mcLocG = rot_apply(ROT_A_CURRENT_COMP, coords_set(0,G,0));
  coords_get(mcLocG, &Gx, &Gy, &Gz);

  strcpy(GVars.compcurname, NAME_CURRENT_COMP);

  if (l > 0 && nelements > 0) {

    Gravity_guide_Init(&GVars,
      w1, h1, w2, h2, l, R0,
      Qc, alpha, m, W, nslit, d,
      Gx, Gy, Gz, mleft, mright, mtop,
      mbottom, nhslit, wavy_lr, wavy_tb, wavy_z, wavy,
      chamfers_z, chamfers_lr, chamfers_tb, chamfers,nu,phase,aleft,aright,atop,abottom);
    if (!G) for (i=0; i<5; GVars.A[i++] = 0);
    if (GVars.fc_freq != 0 || GVars.fc_phase != 0) {
      if (w1 != w2 || h1 != h2)
      exit(fprintf(stderr,"Guide_gravity: %s: rotating slit pack must be straight (w1=w2 and h1=h2).\n", NAME_CURRENT_COMP));
      printf("Guide_gravity: %s: Fermi Chopper mode: frequency=%g [Hz] phase=%g [deg]\n",
        NAME_CURRENT_COMP, GVars.fc_freq, GVars.fc_phase);
    }
  } else printf("Guide_gravity: %s: unactivated (l=0 or nelements=0)\n", NAME_CURRENT_COMP);

}
#line 14063 "./ISIS_SANS2d.c"
#undef reflect
#undef phase
#undef nu
#undef nelements
#undef chamfers_tb
#undef chamfers_lr
#undef chamfers_z
#undef chamfers
#undef wavy_lr
#undef wavy_tb
#undef wavy_z
#undef wavy
#undef abottom
#undef atop
#undef aright
#undef aleft
#undef G
#undef nhslit
#undef mbottom
#undef mtop
#undef mright
#undef mleft
#undef d
#undef nslit
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
#undef pTable
#undef GVars
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component bender9. */
  SIG_MESSAGE("bender9 (Init)");
#define mccompcurname  bender9
#define mccompcurtype  Guide_gravity
#define mccompcurindex 13
#define GVars mccbender9_GVars
#define pTable mccbender9_pTable
#define w1 mccbender9_w1
#define h1 mccbender9_h1
#define w2 mccbender9_w2
#define h2 mccbender9_h2
#define l mccbender9_l
#define R0 mccbender9_R0
#define Qc mccbender9_Qc
#define alpha mccbender9_alpha
#define m mccbender9_m
#define W mccbender9_W
#define nslit mccbender9_nslit
#define d mccbender9_d
#define mleft mccbender9_mleft
#define mright mccbender9_mright
#define mtop mccbender9_mtop
#define mbottom mccbender9_mbottom
#define nhslit mccbender9_nhslit
#define G mccbender9_G
#define aleft mccbender9_aleft
#define aright mccbender9_aright
#define atop mccbender9_atop
#define abottom mccbender9_abottom
#define wavy mccbender9_wavy
#define wavy_z mccbender9_wavy_z
#define wavy_tb mccbender9_wavy_tb
#define wavy_lr mccbender9_wavy_lr
#define chamfers mccbender9_chamfers
#define chamfers_z mccbender9_chamfers_z
#define chamfers_lr mccbender9_chamfers_lr
#define chamfers_tb mccbender9_chamfers_tb
#define nelements mccbender9_nelements
#define nu mccbender9_nu
#define phase mccbender9_phase
#define reflect mccbender9_reflect
#line 339 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide_gravity.comp"
{
  double Gx=0, Gy=-GRAVITY, Gz=0;
  Coords mcLocG;
  int i;

  if (reflect && strlen(reflect) && strcmp(reflect,"NULL") && strcmp(reflect,"0")) {
    if (Table_Read(&pTable, reflect, 1) <= 0) /* read 1st block data from file into pTable */
      exit(fprintf(stderr,"Guide_gravity: %s: can not read file %s\n", NAME_CURRENT_COMP, reflect));
  } else {
    if (W < 0 || R0 < 0 || Qc < 0)
    { fprintf(stderr,"Guide_gravity: %s: W R0 Qc must be >0.\n", NAME_CURRENT_COMP);
      exit(-1); }
  }

  if (nslit <= 0 || nhslit <= 0)
  { fprintf(stderr,"Guide_gravity: %s: nslit nhslit must be >0.\n", NAME_CURRENT_COMP);
    exit(-1); }

  if (!w1 || !h1)
  { fprintf(stderr,"Guide_gravity: %s: input window is closed (w1=h1=0).\n", NAME_CURRENT_COMP);
    exit(-1); }

  if (d*nslit > w1) exit(fprintf(stderr, "Guide_gravity: %s: absorbing walls fill input window. No space left for transmission (d*nslit > w1).\n", NAME_CURRENT_COMP));

  if (!w2) w2=w1;
  if (!h2) h2=h1;

  if (mcgravitation) G=-GRAVITY;
  mcLocG = rot_apply(ROT_A_CURRENT_COMP, coords_set(0,G,0));
  coords_get(mcLocG, &Gx, &Gy, &Gz);

  strcpy(GVars.compcurname, NAME_CURRENT_COMP);

  if (l > 0 && nelements > 0) {

    Gravity_guide_Init(&GVars,
      w1, h1, w2, h2, l, R0,
      Qc, alpha, m, W, nslit, d,
      Gx, Gy, Gz, mleft, mright, mtop,
      mbottom, nhslit, wavy_lr, wavy_tb, wavy_z, wavy,
      chamfers_z, chamfers_lr, chamfers_tb, chamfers,nu,phase,aleft,aright,atop,abottom);
    if (!G) for (i=0; i<5; GVars.A[i++] = 0);
    if (GVars.fc_freq != 0 || GVars.fc_phase != 0) {
      if (w1 != w2 || h1 != h2)
      exit(fprintf(stderr,"Guide_gravity: %s: rotating slit pack must be straight (w1=w2 and h1=h2).\n", NAME_CURRENT_COMP));
      printf("Guide_gravity: %s: Fermi Chopper mode: frequency=%g [Hz] phase=%g [deg]\n",
        NAME_CURRENT_COMP, GVars.fc_freq, GVars.fc_phase);
    }
  } else printf("Guide_gravity: %s: unactivated (l=0 or nelements=0)\n", NAME_CURRENT_COMP);

}
#line 14197 "./ISIS_SANS2d.c"
#undef reflect
#undef phase
#undef nu
#undef nelements
#undef chamfers_tb
#undef chamfers_lr
#undef chamfers_z
#undef chamfers
#undef wavy_lr
#undef wavy_tb
#undef wavy_z
#undef wavy
#undef abottom
#undef atop
#undef aright
#undef aleft
#undef G
#undef nhslit
#undef mbottom
#undef mtop
#undef mright
#undef mleft
#undef d
#undef nslit
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
#undef pTable
#undef GVars
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component bender10. */
  SIG_MESSAGE("bender10 (Init)");
#define mccompcurname  bender10
#define mccompcurtype  Guide_gravity
#define mccompcurindex 14
#define GVars mccbender10_GVars
#define pTable mccbender10_pTable
#define w1 mccbender10_w1
#define h1 mccbender10_h1
#define w2 mccbender10_w2
#define h2 mccbender10_h2
#define l mccbender10_l
#define R0 mccbender10_R0
#define Qc mccbender10_Qc
#define alpha mccbender10_alpha
#define m mccbender10_m
#define W mccbender10_W
#define nslit mccbender10_nslit
#define d mccbender10_d
#define mleft mccbender10_mleft
#define mright mccbender10_mright
#define mtop mccbender10_mtop
#define mbottom mccbender10_mbottom
#define nhslit mccbender10_nhslit
#define G mccbender10_G
#define aleft mccbender10_aleft
#define aright mccbender10_aright
#define atop mccbender10_atop
#define abottom mccbender10_abottom
#define wavy mccbender10_wavy
#define wavy_z mccbender10_wavy_z
#define wavy_tb mccbender10_wavy_tb
#define wavy_lr mccbender10_wavy_lr
#define chamfers mccbender10_chamfers
#define chamfers_z mccbender10_chamfers_z
#define chamfers_lr mccbender10_chamfers_lr
#define chamfers_tb mccbender10_chamfers_tb
#define nelements mccbender10_nelements
#define nu mccbender10_nu
#define phase mccbender10_phase
#define reflect mccbender10_reflect
#line 339 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide_gravity.comp"
{
  double Gx=0, Gy=-GRAVITY, Gz=0;
  Coords mcLocG;
  int i;

  if (reflect && strlen(reflect) && strcmp(reflect,"NULL") && strcmp(reflect,"0")) {
    if (Table_Read(&pTable, reflect, 1) <= 0) /* read 1st block data from file into pTable */
      exit(fprintf(stderr,"Guide_gravity: %s: can not read file %s\n", NAME_CURRENT_COMP, reflect));
  } else {
    if (W < 0 || R0 < 0 || Qc < 0)
    { fprintf(stderr,"Guide_gravity: %s: W R0 Qc must be >0.\n", NAME_CURRENT_COMP);
      exit(-1); }
  }

  if (nslit <= 0 || nhslit <= 0)
  { fprintf(stderr,"Guide_gravity: %s: nslit nhslit must be >0.\n", NAME_CURRENT_COMP);
    exit(-1); }

  if (!w1 || !h1)
  { fprintf(stderr,"Guide_gravity: %s: input window is closed (w1=h1=0).\n", NAME_CURRENT_COMP);
    exit(-1); }

  if (d*nslit > w1) exit(fprintf(stderr, "Guide_gravity: %s: absorbing walls fill input window. No space left for transmission (d*nslit > w1).\n", NAME_CURRENT_COMP));

  if (!w2) w2=w1;
  if (!h2) h2=h1;

  if (mcgravitation) G=-GRAVITY;
  mcLocG = rot_apply(ROT_A_CURRENT_COMP, coords_set(0,G,0));
  coords_get(mcLocG, &Gx, &Gy, &Gz);

  strcpy(GVars.compcurname, NAME_CURRENT_COMP);

  if (l > 0 && nelements > 0) {

    Gravity_guide_Init(&GVars,
      w1, h1, w2, h2, l, R0,
      Qc, alpha, m, W, nslit, d,
      Gx, Gy, Gz, mleft, mright, mtop,
      mbottom, nhslit, wavy_lr, wavy_tb, wavy_z, wavy,
      chamfers_z, chamfers_lr, chamfers_tb, chamfers,nu,phase,aleft,aright,atop,abottom);
    if (!G) for (i=0; i<5; GVars.A[i++] = 0);
    if (GVars.fc_freq != 0 || GVars.fc_phase != 0) {
      if (w1 != w2 || h1 != h2)
      exit(fprintf(stderr,"Guide_gravity: %s: rotating slit pack must be straight (w1=w2 and h1=h2).\n", NAME_CURRENT_COMP));
      printf("Guide_gravity: %s: Fermi Chopper mode: frequency=%g [Hz] phase=%g [deg]\n",
        NAME_CURRENT_COMP, GVars.fc_freq, GVars.fc_phase);
    }
  } else printf("Guide_gravity: %s: unactivated (l=0 or nelements=0)\n", NAME_CURRENT_COMP);

}
#line 14331 "./ISIS_SANS2d.c"
#undef reflect
#undef phase
#undef nu
#undef nelements
#undef chamfers_tb
#undef chamfers_lr
#undef chamfers_z
#undef chamfers
#undef wavy_lr
#undef wavy_tb
#undef wavy_z
#undef wavy
#undef abottom
#undef atop
#undef aright
#undef aleft
#undef G
#undef nhslit
#undef mbottom
#undef mtop
#undef mright
#undef mleft
#undef d
#undef nslit
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
#undef pTable
#undef GVars
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component lmonb. */
  SIG_MESSAGE("lmonb (Init)");
#define mccompcurname  lmonb
#define mccompcurtype  L_monitor
#define mccompcurindex 15
#define nL mcclmonb_nL
#define L_N mcclmonb_L_N
#define L_p mcclmonb_L_p
#define L_p2 mcclmonb_L_p2
#define filename mcclmonb_filename
#define xmin mcclmonb_xmin
#define xmax mcclmonb_xmax
#define ymin mcclmonb_ymin
#define ymax mcclmonb_ymax
#define xwidth mcclmonb_xwidth
#define yheight mcclmonb_yheight
#define Lmin mcclmonb_Lmin
#define Lmax mcclmonb_Lmax
#define restore_neutron mcclmonb_restore_neutron
#define nowritefile mcclmonb_nowritefile
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
#line 14413 "./ISIS_SANS2d.c"
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

  /* Initializations for component psd2. */
  SIG_MESSAGE("psd2 (Init)");
#define mccompcurname  psd2
#define mccompcurtype  PSD_monitor
#define mccompcurindex 16
#define PSD_N mccpsd2_PSD_N
#define PSD_p mccpsd2_PSD_p
#define PSD_p2 mccpsd2_PSD_p2
#define nx mccpsd2_nx
#define ny mccpsd2_ny
#define filename mccpsd2_filename
#define xmin mccpsd2_xmin
#define xmax mccpsd2_xmax
#define ymin mccpsd2_ymin
#define ymax mccpsd2_ymax
#define xwidth mccpsd2_xwidth
#define yheight mccpsd2_yheight
#define restore_neutron mccpsd2_restore_neutron
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
#line 14476 "./ISIS_SANS2d.c"
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

  /* Initializations for component guide_in. */
  SIG_MESSAGE("guide_in (Init)");
#define mccompcurname  guide_in
#define mccompcurtype  Slit
#define mccompcurindex 17
#define xmin mccguide_in_xmin
#define xmax mccguide_in_xmax
#define ymin mccguide_in_ymin
#define ymax mccguide_in_ymax
#define radius mccguide_in_radius
#define xwidth mccguide_in_xwidth
#define yheight mccguide_in_yheight
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
#line 14526 "./ISIS_SANS2d.c"
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

  /* Initializations for component guide_straight1. */
  SIG_MESSAGE("guide_straight1 (Init)");
#define mccompcurname  guide_straight1
#define mccompcurtype  Guide_gravity
#define mccompcurindex 18
#define GVars mccguide_straight1_GVars
#define pTable mccguide_straight1_pTable
#define w1 mccguide_straight1_w1
#define h1 mccguide_straight1_h1
#define w2 mccguide_straight1_w2
#define h2 mccguide_straight1_h2
#define l mccguide_straight1_l
#define R0 mccguide_straight1_R0
#define Qc mccguide_straight1_Qc
#define alpha mccguide_straight1_alpha
#define m mccguide_straight1_m
#define W mccguide_straight1_W
#define nslit mccguide_straight1_nslit
#define d mccguide_straight1_d
#define mleft mccguide_straight1_mleft
#define mright mccguide_straight1_mright
#define mtop mccguide_straight1_mtop
#define mbottom mccguide_straight1_mbottom
#define nhslit mccguide_straight1_nhslit
#define G mccguide_straight1_G
#define aleft mccguide_straight1_aleft
#define aright mccguide_straight1_aright
#define atop mccguide_straight1_atop
#define abottom mccguide_straight1_abottom
#define wavy mccguide_straight1_wavy
#define wavy_z mccguide_straight1_wavy_z
#define wavy_tb mccguide_straight1_wavy_tb
#define wavy_lr mccguide_straight1_wavy_lr
#define chamfers mccguide_straight1_chamfers
#define chamfers_z mccguide_straight1_chamfers_z
#define chamfers_lr mccguide_straight1_chamfers_lr
#define chamfers_tb mccguide_straight1_chamfers_tb
#define nelements mccguide_straight1_nelements
#define nu mccguide_straight1_nu
#define phase mccguide_straight1_phase
#define reflect mccguide_straight1_reflect
#line 339 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide_gravity.comp"
{
  double Gx=0, Gy=-GRAVITY, Gz=0;
  Coords mcLocG;
  int i;

  if (reflect && strlen(reflect) && strcmp(reflect,"NULL") && strcmp(reflect,"0")) {
    if (Table_Read(&pTable, reflect, 1) <= 0) /* read 1st block data from file into pTable */
      exit(fprintf(stderr,"Guide_gravity: %s: can not read file %s\n", NAME_CURRENT_COMP, reflect));
  } else {
    if (W < 0 || R0 < 0 || Qc < 0)
    { fprintf(stderr,"Guide_gravity: %s: W R0 Qc must be >0.\n", NAME_CURRENT_COMP);
      exit(-1); }
  }

  if (nslit <= 0 || nhslit <= 0)
  { fprintf(stderr,"Guide_gravity: %s: nslit nhslit must be >0.\n", NAME_CURRENT_COMP);
    exit(-1); }

  if (!w1 || !h1)
  { fprintf(stderr,"Guide_gravity: %s: input window is closed (w1=h1=0).\n", NAME_CURRENT_COMP);
    exit(-1); }

  if (d*nslit > w1) exit(fprintf(stderr, "Guide_gravity: %s: absorbing walls fill input window. No space left for transmission (d*nslit > w1).\n", NAME_CURRENT_COMP));

  if (!w2) w2=w1;
  if (!h2) h2=h1;

  if (mcgravitation) G=-GRAVITY;
  mcLocG = rot_apply(ROT_A_CURRENT_COMP, coords_set(0,G,0));
  coords_get(mcLocG, &Gx, &Gy, &Gz);

  strcpy(GVars.compcurname, NAME_CURRENT_COMP);

  if (l > 0 && nelements > 0) {

    Gravity_guide_Init(&GVars,
      w1, h1, w2, h2, l, R0,
      Qc, alpha, m, W, nslit, d,
      Gx, Gy, Gz, mleft, mright, mtop,
      mbottom, nhslit, wavy_lr, wavy_tb, wavy_z, wavy,
      chamfers_z, chamfers_lr, chamfers_tb, chamfers,nu,phase,aleft,aright,atop,abottom);
    if (!G) for (i=0; i<5; GVars.A[i++] = 0);
    if (GVars.fc_freq != 0 || GVars.fc_phase != 0) {
      if (w1 != w2 || h1 != h2)
      exit(fprintf(stderr,"Guide_gravity: %s: rotating slit pack must be straight (w1=w2 and h1=h2).\n", NAME_CURRENT_COMP));
      printf("Guide_gravity: %s: Fermi Chopper mode: frequency=%g [Hz] phase=%g [deg]\n",
        NAME_CURRENT_COMP, GVars.fc_freq, GVars.fc_phase);
    }
  } else printf("Guide_gravity: %s: unactivated (l=0 or nelements=0)\n", NAME_CURRENT_COMP);

}
#line 14631 "./ISIS_SANS2d.c"
#undef reflect
#undef phase
#undef nu
#undef nelements
#undef chamfers_tb
#undef chamfers_lr
#undef chamfers_z
#undef chamfers
#undef wavy_lr
#undef wavy_tb
#undef wavy_z
#undef wavy
#undef abottom
#undef atop
#undef aright
#undef aleft
#undef G
#undef nhslit
#undef mbottom
#undef mtop
#undef mright
#undef mleft
#undef d
#undef nslit
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
#undef pTable
#undef GVars
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component guide_straight2. */
  SIG_MESSAGE("guide_straight2 (Init)");
#define mccompcurname  guide_straight2
#define mccompcurtype  Guide_gravity
#define mccompcurindex 19
#define GVars mccguide_straight2_GVars
#define pTable mccguide_straight2_pTable
#define w1 mccguide_straight2_w1
#define h1 mccguide_straight2_h1
#define w2 mccguide_straight2_w2
#define h2 mccguide_straight2_h2
#define l mccguide_straight2_l
#define R0 mccguide_straight2_R0
#define Qc mccguide_straight2_Qc
#define alpha mccguide_straight2_alpha
#define m mccguide_straight2_m
#define W mccguide_straight2_W
#define nslit mccguide_straight2_nslit
#define d mccguide_straight2_d
#define mleft mccguide_straight2_mleft
#define mright mccguide_straight2_mright
#define mtop mccguide_straight2_mtop
#define mbottom mccguide_straight2_mbottom
#define nhslit mccguide_straight2_nhslit
#define G mccguide_straight2_G
#define aleft mccguide_straight2_aleft
#define aright mccguide_straight2_aright
#define atop mccguide_straight2_atop
#define abottom mccguide_straight2_abottom
#define wavy mccguide_straight2_wavy
#define wavy_z mccguide_straight2_wavy_z
#define wavy_tb mccguide_straight2_wavy_tb
#define wavy_lr mccguide_straight2_wavy_lr
#define chamfers mccguide_straight2_chamfers
#define chamfers_z mccguide_straight2_chamfers_z
#define chamfers_lr mccguide_straight2_chamfers_lr
#define chamfers_tb mccguide_straight2_chamfers_tb
#define nelements mccguide_straight2_nelements
#define nu mccguide_straight2_nu
#define phase mccguide_straight2_phase
#define reflect mccguide_straight2_reflect
#line 339 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide_gravity.comp"
{
  double Gx=0, Gy=-GRAVITY, Gz=0;
  Coords mcLocG;
  int i;

  if (reflect && strlen(reflect) && strcmp(reflect,"NULL") && strcmp(reflect,"0")) {
    if (Table_Read(&pTable, reflect, 1) <= 0) /* read 1st block data from file into pTable */
      exit(fprintf(stderr,"Guide_gravity: %s: can not read file %s\n", NAME_CURRENT_COMP, reflect));
  } else {
    if (W < 0 || R0 < 0 || Qc < 0)
    { fprintf(stderr,"Guide_gravity: %s: W R0 Qc must be >0.\n", NAME_CURRENT_COMP);
      exit(-1); }
  }

  if (nslit <= 0 || nhslit <= 0)
  { fprintf(stderr,"Guide_gravity: %s: nslit nhslit must be >0.\n", NAME_CURRENT_COMP);
    exit(-1); }

  if (!w1 || !h1)
  { fprintf(stderr,"Guide_gravity: %s: input window is closed (w1=h1=0).\n", NAME_CURRENT_COMP);
    exit(-1); }

  if (d*nslit > w1) exit(fprintf(stderr, "Guide_gravity: %s: absorbing walls fill input window. No space left for transmission (d*nslit > w1).\n", NAME_CURRENT_COMP));

  if (!w2) w2=w1;
  if (!h2) h2=h1;

  if (mcgravitation) G=-GRAVITY;
  mcLocG = rot_apply(ROT_A_CURRENT_COMP, coords_set(0,G,0));
  coords_get(mcLocG, &Gx, &Gy, &Gz);

  strcpy(GVars.compcurname, NAME_CURRENT_COMP);

  if (l > 0 && nelements > 0) {

    Gravity_guide_Init(&GVars,
      w1, h1, w2, h2, l, R0,
      Qc, alpha, m, W, nslit, d,
      Gx, Gy, Gz, mleft, mright, mtop,
      mbottom, nhslit, wavy_lr, wavy_tb, wavy_z, wavy,
      chamfers_z, chamfers_lr, chamfers_tb, chamfers,nu,phase,aleft,aright,atop,abottom);
    if (!G) for (i=0; i<5; GVars.A[i++] = 0);
    if (GVars.fc_freq != 0 || GVars.fc_phase != 0) {
      if (w1 != w2 || h1 != h2)
      exit(fprintf(stderr,"Guide_gravity: %s: rotating slit pack must be straight (w1=w2 and h1=h2).\n", NAME_CURRENT_COMP));
      printf("Guide_gravity: %s: Fermi Chopper mode: frequency=%g [Hz] phase=%g [deg]\n",
        NAME_CURRENT_COMP, GVars.fc_freq, GVars.fc_phase);
    }
  } else printf("Guide_gravity: %s: unactivated (l=0 or nelements=0)\n", NAME_CURRENT_COMP);

}
#line 14765 "./ISIS_SANS2d.c"
#undef reflect
#undef phase
#undef nu
#undef nelements
#undef chamfers_tb
#undef chamfers_lr
#undef chamfers_z
#undef chamfers
#undef wavy_lr
#undef wavy_tb
#undef wavy_z
#undef wavy
#undef abottom
#undef atop
#undef aright
#undef aleft
#undef G
#undef nhslit
#undef mbottom
#undef mtop
#undef mright
#undef mleft
#undef d
#undef nslit
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
#undef pTable
#undef GVars
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component guide_straight3. */
  SIG_MESSAGE("guide_straight3 (Init)");
#define mccompcurname  guide_straight3
#define mccompcurtype  Guide_gravity
#define mccompcurindex 20
#define GVars mccguide_straight3_GVars
#define pTable mccguide_straight3_pTable
#define w1 mccguide_straight3_w1
#define h1 mccguide_straight3_h1
#define w2 mccguide_straight3_w2
#define h2 mccguide_straight3_h2
#define l mccguide_straight3_l
#define R0 mccguide_straight3_R0
#define Qc mccguide_straight3_Qc
#define alpha mccguide_straight3_alpha
#define m mccguide_straight3_m
#define W mccguide_straight3_W
#define nslit mccguide_straight3_nslit
#define d mccguide_straight3_d
#define mleft mccguide_straight3_mleft
#define mright mccguide_straight3_mright
#define mtop mccguide_straight3_mtop
#define mbottom mccguide_straight3_mbottom
#define nhslit mccguide_straight3_nhslit
#define G mccguide_straight3_G
#define aleft mccguide_straight3_aleft
#define aright mccguide_straight3_aright
#define atop mccguide_straight3_atop
#define abottom mccguide_straight3_abottom
#define wavy mccguide_straight3_wavy
#define wavy_z mccguide_straight3_wavy_z
#define wavy_tb mccguide_straight3_wavy_tb
#define wavy_lr mccguide_straight3_wavy_lr
#define chamfers mccguide_straight3_chamfers
#define chamfers_z mccguide_straight3_chamfers_z
#define chamfers_lr mccguide_straight3_chamfers_lr
#define chamfers_tb mccguide_straight3_chamfers_tb
#define nelements mccguide_straight3_nelements
#define nu mccguide_straight3_nu
#define phase mccguide_straight3_phase
#define reflect mccguide_straight3_reflect
#line 339 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide_gravity.comp"
{
  double Gx=0, Gy=-GRAVITY, Gz=0;
  Coords mcLocG;
  int i;

  if (reflect && strlen(reflect) && strcmp(reflect,"NULL") && strcmp(reflect,"0")) {
    if (Table_Read(&pTable, reflect, 1) <= 0) /* read 1st block data from file into pTable */
      exit(fprintf(stderr,"Guide_gravity: %s: can not read file %s\n", NAME_CURRENT_COMP, reflect));
  } else {
    if (W < 0 || R0 < 0 || Qc < 0)
    { fprintf(stderr,"Guide_gravity: %s: W R0 Qc must be >0.\n", NAME_CURRENT_COMP);
      exit(-1); }
  }

  if (nslit <= 0 || nhslit <= 0)
  { fprintf(stderr,"Guide_gravity: %s: nslit nhslit must be >0.\n", NAME_CURRENT_COMP);
    exit(-1); }

  if (!w1 || !h1)
  { fprintf(stderr,"Guide_gravity: %s: input window is closed (w1=h1=0).\n", NAME_CURRENT_COMP);
    exit(-1); }

  if (d*nslit > w1) exit(fprintf(stderr, "Guide_gravity: %s: absorbing walls fill input window. No space left for transmission (d*nslit > w1).\n", NAME_CURRENT_COMP));

  if (!w2) w2=w1;
  if (!h2) h2=h1;

  if (mcgravitation) G=-GRAVITY;
  mcLocG = rot_apply(ROT_A_CURRENT_COMP, coords_set(0,G,0));
  coords_get(mcLocG, &Gx, &Gy, &Gz);

  strcpy(GVars.compcurname, NAME_CURRENT_COMP);

  if (l > 0 && nelements > 0) {

    Gravity_guide_Init(&GVars,
      w1, h1, w2, h2, l, R0,
      Qc, alpha, m, W, nslit, d,
      Gx, Gy, Gz, mleft, mright, mtop,
      mbottom, nhslit, wavy_lr, wavy_tb, wavy_z, wavy,
      chamfers_z, chamfers_lr, chamfers_tb, chamfers,nu,phase,aleft,aright,atop,abottom);
    if (!G) for (i=0; i<5; GVars.A[i++] = 0);
    if (GVars.fc_freq != 0 || GVars.fc_phase != 0) {
      if (w1 != w2 || h1 != h2)
      exit(fprintf(stderr,"Guide_gravity: %s: rotating slit pack must be straight (w1=w2 and h1=h2).\n", NAME_CURRENT_COMP));
      printf("Guide_gravity: %s: Fermi Chopper mode: frequency=%g [Hz] phase=%g [deg]\n",
        NAME_CURRENT_COMP, GVars.fc_freq, GVars.fc_phase);
    }
  } else printf("Guide_gravity: %s: unactivated (l=0 or nelements=0)\n", NAME_CURRENT_COMP);

}
#line 14899 "./ISIS_SANS2d.c"
#undef reflect
#undef phase
#undef nu
#undef nelements
#undef chamfers_tb
#undef chamfers_lr
#undef chamfers_z
#undef chamfers
#undef wavy_lr
#undef wavy_tb
#undef wavy_z
#undef wavy
#undef abottom
#undef atop
#undef aright
#undef aleft
#undef G
#undef nhslit
#undef mbottom
#undef mtop
#undef mright
#undef mleft
#undef d
#undef nslit
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
#undef pTable
#undef GVars
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component guide_straight4. */
  SIG_MESSAGE("guide_straight4 (Init)");
#define mccompcurname  guide_straight4
#define mccompcurtype  Guide_gravity
#define mccompcurindex 21
#define GVars mccguide_straight4_GVars
#define pTable mccguide_straight4_pTable
#define w1 mccguide_straight4_w1
#define h1 mccguide_straight4_h1
#define w2 mccguide_straight4_w2
#define h2 mccguide_straight4_h2
#define l mccguide_straight4_l
#define R0 mccguide_straight4_R0
#define Qc mccguide_straight4_Qc
#define alpha mccguide_straight4_alpha
#define m mccguide_straight4_m
#define W mccguide_straight4_W
#define nslit mccguide_straight4_nslit
#define d mccguide_straight4_d
#define mleft mccguide_straight4_mleft
#define mright mccguide_straight4_mright
#define mtop mccguide_straight4_mtop
#define mbottom mccguide_straight4_mbottom
#define nhslit mccguide_straight4_nhslit
#define G mccguide_straight4_G
#define aleft mccguide_straight4_aleft
#define aright mccguide_straight4_aright
#define atop mccguide_straight4_atop
#define abottom mccguide_straight4_abottom
#define wavy mccguide_straight4_wavy
#define wavy_z mccguide_straight4_wavy_z
#define wavy_tb mccguide_straight4_wavy_tb
#define wavy_lr mccguide_straight4_wavy_lr
#define chamfers mccguide_straight4_chamfers
#define chamfers_z mccguide_straight4_chamfers_z
#define chamfers_lr mccguide_straight4_chamfers_lr
#define chamfers_tb mccguide_straight4_chamfers_tb
#define nelements mccguide_straight4_nelements
#define nu mccguide_straight4_nu
#define phase mccguide_straight4_phase
#define reflect mccguide_straight4_reflect
#line 339 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide_gravity.comp"
{
  double Gx=0, Gy=-GRAVITY, Gz=0;
  Coords mcLocG;
  int i;

  if (reflect && strlen(reflect) && strcmp(reflect,"NULL") && strcmp(reflect,"0")) {
    if (Table_Read(&pTable, reflect, 1) <= 0) /* read 1st block data from file into pTable */
      exit(fprintf(stderr,"Guide_gravity: %s: can not read file %s\n", NAME_CURRENT_COMP, reflect));
  } else {
    if (W < 0 || R0 < 0 || Qc < 0)
    { fprintf(stderr,"Guide_gravity: %s: W R0 Qc must be >0.\n", NAME_CURRENT_COMP);
      exit(-1); }
  }

  if (nslit <= 0 || nhslit <= 0)
  { fprintf(stderr,"Guide_gravity: %s: nslit nhslit must be >0.\n", NAME_CURRENT_COMP);
    exit(-1); }

  if (!w1 || !h1)
  { fprintf(stderr,"Guide_gravity: %s: input window is closed (w1=h1=0).\n", NAME_CURRENT_COMP);
    exit(-1); }

  if (d*nslit > w1) exit(fprintf(stderr, "Guide_gravity: %s: absorbing walls fill input window. No space left for transmission (d*nslit > w1).\n", NAME_CURRENT_COMP));

  if (!w2) w2=w1;
  if (!h2) h2=h1;

  if (mcgravitation) G=-GRAVITY;
  mcLocG = rot_apply(ROT_A_CURRENT_COMP, coords_set(0,G,0));
  coords_get(mcLocG, &Gx, &Gy, &Gz);

  strcpy(GVars.compcurname, NAME_CURRENT_COMP);

  if (l > 0 && nelements > 0) {

    Gravity_guide_Init(&GVars,
      w1, h1, w2, h2, l, R0,
      Qc, alpha, m, W, nslit, d,
      Gx, Gy, Gz, mleft, mright, mtop,
      mbottom, nhslit, wavy_lr, wavy_tb, wavy_z, wavy,
      chamfers_z, chamfers_lr, chamfers_tb, chamfers,nu,phase,aleft,aright,atop,abottom);
    if (!G) for (i=0; i<5; GVars.A[i++] = 0);
    if (GVars.fc_freq != 0 || GVars.fc_phase != 0) {
      if (w1 != w2 || h1 != h2)
      exit(fprintf(stderr,"Guide_gravity: %s: rotating slit pack must be straight (w1=w2 and h1=h2).\n", NAME_CURRENT_COMP));
      printf("Guide_gravity: %s: Fermi Chopper mode: frequency=%g [Hz] phase=%g [deg]\n",
        NAME_CURRENT_COMP, GVars.fc_freq, GVars.fc_phase);
    }
  } else printf("Guide_gravity: %s: unactivated (l=0 or nelements=0)\n", NAME_CURRENT_COMP);

}
#line 15033 "./ISIS_SANS2d.c"
#undef reflect
#undef phase
#undef nu
#undef nelements
#undef chamfers_tb
#undef chamfers_lr
#undef chamfers_z
#undef chamfers
#undef wavy_lr
#undef wavy_tb
#undef wavy_z
#undef wavy
#undef abottom
#undef atop
#undef aright
#undef aleft
#undef G
#undef nhslit
#undef mbottom
#undef mtop
#undef mright
#undef mleft
#undef d
#undef nslit
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
#undef pTable
#undef GVars
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component psd3. */
  SIG_MESSAGE("psd3 (Init)");
#define mccompcurname  psd3
#define mccompcurtype  PSD_monitor
#define mccompcurindex 22
#define PSD_N mccpsd3_PSD_N
#define PSD_p mccpsd3_PSD_p
#define PSD_p2 mccpsd3_PSD_p2
#define nx mccpsd3_nx
#define ny mccpsd3_ny
#define filename mccpsd3_filename
#define xmin mccpsd3_xmin
#define xmax mccpsd3_xmax
#define ymin mccpsd3_ymin
#define ymax mccpsd3_ymax
#define xwidth mccpsd3_xwidth
#define yheight mccpsd3_yheight
#define restore_neutron mccpsd3_restore_neutron
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
#line 15117 "./ISIS_SANS2d.c"
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

  /* Initializations for component aperture1. */
  SIG_MESSAGE("aperture1 (Init)");
#define mccompcurname  aperture1
#define mccompcurtype  Slit
#define mccompcurindex 23
#define xmin mccaperture1_xmin
#define xmax mccaperture1_xmax
#define ymin mccaperture1_ymin
#define ymax mccaperture1_ymax
#define radius mccaperture1_radius
#define xwidth mccaperture1_xwidth
#define yheight mccaperture1_yheight
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
#line 15167 "./ISIS_SANS2d.c"
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

  /* Initializations for component lmonitor2. */
  SIG_MESSAGE("lmonitor2 (Init)");
#define mccompcurname  lmonitor2
#define mccompcurtype  L_monitor
#define mccompcurindex 24
#define nL mcclmonitor2_nL
#define L_N mcclmonitor2_L_N
#define L_p mcclmonitor2_L_p
#define L_p2 mcclmonitor2_L_p2
#define filename mcclmonitor2_filename
#define xmin mcclmonitor2_xmin
#define xmax mcclmonitor2_xmax
#define ymin mcclmonitor2_ymin
#define ymax mcclmonitor2_ymax
#define xwidth mcclmonitor2_xwidth
#define yheight mcclmonitor2_yheight
#define Lmin mcclmonitor2_Lmin
#define Lmax mcclmonitor2_Lmax
#define restore_neutron mcclmonitor2_restore_neutron
#define nowritefile mcclmonitor2_nowritefile
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
#line 15220 "./ISIS_SANS2d.c"
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

  /* Initializations for component S6. */
  SIG_MESSAGE("S6 (Init)");
#define mccompcurname  S6
#define mccompcurtype  Slit
#define mccompcurindex 25
#define xmin mccS6_xmin
#define xmax mccS6_xmax
#define ymin mccS6_ymin
#define ymax mccS6_ymax
#define radius mccS6_radius
#define xwidth mccS6_xwidth
#define yheight mccS6_yheight
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
#line 15272 "./ISIS_SANS2d.c"
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

  /* Initializations for component APERTURE2. */
  SIG_MESSAGE("APERTURE2 (Init)");
#define mccompcurname  APERTURE2
#define mccompcurtype  Slit
#define mccompcurindex 26
#define xmin mccAPERTURE2_xmin
#define xmax mccAPERTURE2_xmax
#define ymin mccAPERTURE2_ymin
#define ymax mccAPERTURE2_ymax
#define radius mccAPERTURE2_radius
#define xwidth mccAPERTURE2_xwidth
#define yheight mccAPERTURE2_yheight
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
#line 15316 "./ISIS_SANS2d.c"
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
#define mccompcurindex 27
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
#line 15369 "./ISIS_SANS2d.c"
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

  /* Initializations for component psd4. */
  SIG_MESSAGE("psd4 (Init)");
#define mccompcurname  psd4
#define mccompcurtype  PSD_monitor
#define mccompcurindex 28
#define PSD_N mccpsd4_PSD_N
#define PSD_p mccpsd4_PSD_p
#define PSD_p2 mccpsd4_PSD_p2
#define nx mccpsd4_nx
#define ny mccpsd4_ny
#define filename mccpsd4_filename
#define xmin mccpsd4_xmin
#define xmax mccpsd4_xmax
#define ymin mccpsd4_ymin
#define ymax mccpsd4_ymax
#define xwidth mccpsd4_xwidth
#define yheight mccpsd4_yheight
#define restore_neutron mccpsd4_restore_neutron
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
#line 15432 "./ISIS_SANS2d.c"
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
#define mccompcurtype  Arm
#define mccompcurindex 1
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
#line 1060 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../contrib/ISIS_moderator.comp"
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
#line 15738 "./ISIS_SANS2d.c"
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

  /* TRACE Component lmon1 [3] */
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
#define mccompcurname  lmon1
#define mccompcurtype  L_monitor
#define mccompcurindex 3
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
#line 15895 "./ISIS_SANS2d.c"
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

  /* TRACE Component psd1 [4] */
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
#define mccompcurname  psd1
#define mccompcurtype  PSD_monitor
#define mccompcurindex 4
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
#line 16033 "./ISIS_SANS2d.c"
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

  /* TRACE Component bender1 [5] */
  mccoordschange(mcposrbender1, mcrotrbender1,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component bender1 (without coords transformations) */
  mcJumpTrace_bender1:
  SIG_MESSAGE("bender1 (Trace)");
  mcDEBUG_COMP("bender1")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbCompbender1
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
#define mccompcurname  bender1
#define mccompcurtype  Guide_gravity
#define mccompcurindex 5
#define GVars mccbender1_GVars
#define pTable mccbender1_pTable
{   /* Declarations of bender1=Guide_gravity() SETTING parameters. */
MCNUM w1 = mccbender1_w1;
MCNUM h1 = mccbender1_h1;
MCNUM w2 = mccbender1_w2;
MCNUM h2 = mccbender1_h2;
MCNUM l = mccbender1_l;
MCNUM R0 = mccbender1_R0;
MCNUM Qc = mccbender1_Qc;
MCNUM alpha = mccbender1_alpha;
MCNUM m = mccbender1_m;
MCNUM W = mccbender1_W;
MCNUM nslit = mccbender1_nslit;
MCNUM d = mccbender1_d;
MCNUM mleft = mccbender1_mleft;
MCNUM mright = mccbender1_mright;
MCNUM mtop = mccbender1_mtop;
MCNUM mbottom = mccbender1_mbottom;
MCNUM nhslit = mccbender1_nhslit;
MCNUM G = mccbender1_G;
MCNUM aleft = mccbender1_aleft;
MCNUM aright = mccbender1_aright;
MCNUM atop = mccbender1_atop;
MCNUM abottom = mccbender1_abottom;
MCNUM wavy = mccbender1_wavy;
MCNUM wavy_z = mccbender1_wavy_z;
MCNUM wavy_tb = mccbender1_wavy_tb;
MCNUM wavy_lr = mccbender1_wavy_lr;
MCNUM chamfers = mccbender1_chamfers;
MCNUM chamfers_z = mccbender1_chamfers_z;
MCNUM chamfers_lr = mccbender1_chamfers_lr;
MCNUM chamfers_tb = mccbender1_chamfers_tb;
MCNUM nelements = mccbender1_nelements;
MCNUM nu = mccbender1_nu;
MCNUM phase = mccbender1_phase;
char* reflect = mccbender1_reflect;
#line 392 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide_gravity.comp"
{
  if (l > 0 && nelements > 0) {
    double B, C, dt;
    int    ret, bounces = 0, i=0;
    double this_width, this_height;
    double angle=0;

    if (GVars.fc_freq != 0 || GVars.fc_phase != 0) { /* rotate neutron w/r to guide element */
      /* approximation of rotating straight Fermi Chopper */
      Coords   X = coords_set(x,y,z-l/2);  /* current coordinates of neutron in centered static frame */
      Rotation R;
      double dt=(-z+l/2)/vz; /* time shift to each center of slit package */
      angle=fmod(360*GVars.fc_freq*(t+dt)+GVars.fc_phase, 360); /* in deg */
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

    for (i=0; i<7; GVars.N_reflection[i++] = 0);

    /* propagate to box input (with gravitation) in comp local coords */
    /* A = 0.5 n.g; B = n.v; C = n.(r-W); */
    /* 0=Z0 side: n=(0, 0, -l) ; W = (0, 0, 0) (at z=0, guide input)*/
    B = -l*vz; C = -l*z;

    ret = solve_2nd_order(&dt, NULL, GVars.A[0], B, C);
    if (ret==0) ABSORB;

    if (dt>0.0) PROP_GRAV_DT(dt, GVars.gx, GVars.gy, GVars.gz); else if (angle) ABSORB;
    GVars.N_reflection[6]++;

    this_width  = w1;
    this_height = h1;

  /* check if we are in the box input, else absorb */
    if (fabs(x) > this_width/2 || fabs(y) > this_height/2)
      ABSORB;
    else
    {
      double w_edge, w_adj; /* Channel displacement on X */
      double h_edge, h_adj; /* Channel displacement on Y */
      double w_chnum,h_chnum; /* channel indexes */

      SCATTER;

      /* X: Shift origin to center of channel hit (absorb if hit dividing walls) */
      x += w1/2.0;
      w_chnum = floor(x/(GVars.w1c+d));  /* 0= right side, nslit+1=left side  */
      w_edge  = w_chnum*(GVars.w1c+d);
      if(x - w_edge > GVars.w1c)
      {
        x -= w1/2.0; /* Re-adjust origin */
        ABSORB;
      }
      w_adj = w_edge + (GVars.w1c)/2.0;
      x -= w_adj; w_adj -=  w1/2.0;

      /* Y: Shift origin to center of channel hit (absorb if hit dividing walls) */
      y += h1/2.0;
      h_chnum = floor(y/(GVars.h1c+d));  /* 0= lower side, nslit+1=upper side  */
      h_edge  = h_chnum*(GVars.h1c+d);
      if(y - h_edge > GVars.h1c)
      {
        y -= h1/2.0; /* Re-adjust origin */
        ABSORB;
      }
      h_adj = h_edge + (GVars.h1c)/2.0;
      y -= h_adj; h_adj -=  h1/2.0;

      /* neutron is now in the input window of the guide */
      /* do loops on reflections in the box */
      for(;;)
      {
        /* get intersections for all box sides */
        double q, nx,ny,nz;
        double this_length;
        int side=0;

        bounces++;
        /* now look for intersection with guide sides and exit */
        side = Gravity_guide_Trace(&dt, &GVars, x, y, z,
            vx, vy, vz, w_chnum, nslit, h_chnum, nhslit,
            &nx, &ny, &nz);

        /* only positive dt are valid */
        /* exit reflection loops if no intersection (neutron is after box) */
        if (side == 0 || dt <= 0)
          { if (GVars.warnings < 100)
              fprintf(stderr,"%s: warning: neutron has entered guide, but can not exit !\n", GVars.compcurname);
            GVars.warnings++;
            x += w_adj; y += h_adj; ABSORB; } /* should never occur */

        /* propagate to dt */
        PROP_GRAV_DT(dt, GVars.gx, GVars.gy, GVars.gz);

        /* do reflection on speed for l/r/u/d sides */
        if (side == 5) /* neutron reaches end of guide: end loop and exit comp */
          { GVars.N_reflection[side]++; x += w_adj; y += h_adj; SCATTER; x -= w_adj; y -= h_adj; break; }
        /* else reflection on a guide wall */
        if(GVars.M[side] == 0 || Qc == 0 || R0 == 0)  /* walls are absorbing */
          { x += w_adj; y += h_adj; ABSORB; }
        /* handle chamfers */
        this_width = w1+(w2-w1)*z/l;
        this_height= h1+(h2-h1)*z/l;
        this_length= fmod(z, l/nelements);
        /* absorb on input/output of element parts */
        if (GVars.chamfer_z && (this_length<GVars.chamfer_z || this_length>l/nelements-GVars.chamfer_z))
        { x += w_adj; y += h_adj; ABSORB; }
        /* absorb on l/r/t/b sides */
        if (GVars.chamfer_lr && (side==1 || side==2) && (fabs(y+h_adj)>this_height/2-GVars.chamfer_lr))
        { x += w_adj; y += h_adj; ABSORB; }
        if (GVars.chamfer_tb && (side==3 || side==4) && (fabs(x+w_adj)>this_width/2- GVars.chamfer_tb))
        { x += w_adj; y += h_adj; ABSORB; }
        /* change/mirror velocity: h_f = v - n.2*n.v/|n|^2 */
        GVars.N_reflection[side]++; /* GVars.norm_n2 > 0 was checked at INIT */
        /* compute n.v using current values */
        B = scalar_prod(vx,vy,vz,nx,ny,nz);
        dt = 2*B/GVars.norm_n2[side]; /* 2*n.v/|n|^2 */
        vx -= nx*dt;
        vy -= ny*dt;
        vz -= nz*dt;

        /* compute q and modify neutron weight */
        /* scattering q=|n_i-n_f| = V2Q*|vf - v| = V2Q*2*n.v/|n| */
        q = 2*V2Q*fabs(B)/GVars.norm_n[side];

        if (reflect && strlen(reflect) && strcmp(reflect,"NULL") && strcmp(reflect,"0"))
          TableReflecFunc(q, &pTable, &B);
        else {
          double par[] = {R0, Qc, GVars.Alpha[side], GVars.M[side], W};
          StdReflecFunc(q, par, &B);
        }
        if (B <= 0) { x += w_adj; y += h_adj; ABSORB; }
        else p *= B;
        x += w_adj; y += h_adj; SCATTER; x -= w_adj; y -= h_adj;
        GVars.N_reflection[0]++;
        /* go to the next reflection */
        if (bounces > 1000) ABSORB;
      } /* end for */
      x += w_adj; y += h_adj; /* Re-adjust origin after SCATTER */
    }

    if (GVars.fc_freq != 0 || GVars.fc_phase != 0) { /* rotate back neutron w/r to guide element */
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

  } /* if l */
}
#line 16347 "./ISIS_SANS2d.c"
}   /* End of bender1=Guide_gravity() SETTING parameter declarations. */
#undef pTable
#undef GVars
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompbender1:
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

  /* TRACE Component bender2 [6] */
  mccoordschange(mcposrbender2, mcrotrbender2,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component bender2 (without coords transformations) */
  mcJumpTrace_bender2:
  SIG_MESSAGE("bender2 (Trace)");
  mcDEBUG_COMP("bender2")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbCompbender2
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
#define mccompcurname  bender2
#define mccompcurtype  Guide_gravity
#define mccompcurindex 6
#define GVars mccbender2_GVars
#define pTable mccbender2_pTable
{   /* Declarations of bender2=Guide_gravity() SETTING parameters. */
MCNUM w1 = mccbender2_w1;
MCNUM h1 = mccbender2_h1;
MCNUM w2 = mccbender2_w2;
MCNUM h2 = mccbender2_h2;
MCNUM l = mccbender2_l;
MCNUM R0 = mccbender2_R0;
MCNUM Qc = mccbender2_Qc;
MCNUM alpha = mccbender2_alpha;
MCNUM m = mccbender2_m;
MCNUM W = mccbender2_W;
MCNUM nslit = mccbender2_nslit;
MCNUM d = mccbender2_d;
MCNUM mleft = mccbender2_mleft;
MCNUM mright = mccbender2_mright;
MCNUM mtop = mccbender2_mtop;
MCNUM mbottom = mccbender2_mbottom;
MCNUM nhslit = mccbender2_nhslit;
MCNUM G = mccbender2_G;
MCNUM aleft = mccbender2_aleft;
MCNUM aright = mccbender2_aright;
MCNUM atop = mccbender2_atop;
MCNUM abottom = mccbender2_abottom;
MCNUM wavy = mccbender2_wavy;
MCNUM wavy_z = mccbender2_wavy_z;
MCNUM wavy_tb = mccbender2_wavy_tb;
MCNUM wavy_lr = mccbender2_wavy_lr;
MCNUM chamfers = mccbender2_chamfers;
MCNUM chamfers_z = mccbender2_chamfers_z;
MCNUM chamfers_lr = mccbender2_chamfers_lr;
MCNUM chamfers_tb = mccbender2_chamfers_tb;
MCNUM nelements = mccbender2_nelements;
MCNUM nu = mccbender2_nu;
MCNUM phase = mccbender2_phase;
char* reflect = mccbender2_reflect;
#line 392 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide_gravity.comp"
{
  if (l > 0 && nelements > 0) {
    double B, C, dt;
    int    ret, bounces = 0, i=0;
    double this_width, this_height;
    double angle=0;

    if (GVars.fc_freq != 0 || GVars.fc_phase != 0) { /* rotate neutron w/r to guide element */
      /* approximation of rotating straight Fermi Chopper */
      Coords   X = coords_set(x,y,z-l/2);  /* current coordinates of neutron in centered static frame */
      Rotation R;
      double dt=(-z+l/2)/vz; /* time shift to each center of slit package */
      angle=fmod(360*GVars.fc_freq*(t+dt)+GVars.fc_phase, 360); /* in deg */
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

    for (i=0; i<7; GVars.N_reflection[i++] = 0);

    /* propagate to box input (with gravitation) in comp local coords */
    /* A = 0.5 n.g; B = n.v; C = n.(r-W); */
    /* 0=Z0 side: n=(0, 0, -l) ; W = (0, 0, 0) (at z=0, guide input)*/
    B = -l*vz; C = -l*z;

    ret = solve_2nd_order(&dt, NULL, GVars.A[0], B, C);
    if (ret==0) ABSORB;

    if (dt>0.0) PROP_GRAV_DT(dt, GVars.gx, GVars.gy, GVars.gz); else if (angle) ABSORB;
    GVars.N_reflection[6]++;

    this_width  = w1;
    this_height = h1;

  /* check if we are in the box input, else absorb */
    if (fabs(x) > this_width/2 || fabs(y) > this_height/2)
      ABSORB;
    else
    {
      double w_edge, w_adj; /* Channel displacement on X */
      double h_edge, h_adj; /* Channel displacement on Y */
      double w_chnum,h_chnum; /* channel indexes */

      SCATTER;

      /* X: Shift origin to center of channel hit (absorb if hit dividing walls) */
      x += w1/2.0;
      w_chnum = floor(x/(GVars.w1c+d));  /* 0= right side, nslit+1=left side  */
      w_edge  = w_chnum*(GVars.w1c+d);
      if(x - w_edge > GVars.w1c)
      {
        x -= w1/2.0; /* Re-adjust origin */
        ABSORB;
      }
      w_adj = w_edge + (GVars.w1c)/2.0;
      x -= w_adj; w_adj -=  w1/2.0;

      /* Y: Shift origin to center of channel hit (absorb if hit dividing walls) */
      y += h1/2.0;
      h_chnum = floor(y/(GVars.h1c+d));  /* 0= lower side, nslit+1=upper side  */
      h_edge  = h_chnum*(GVars.h1c+d);
      if(y - h_edge > GVars.h1c)
      {
        y -= h1/2.0; /* Re-adjust origin */
        ABSORB;
      }
      h_adj = h_edge + (GVars.h1c)/2.0;
      y -= h_adj; h_adj -=  h1/2.0;

      /* neutron is now in the input window of the guide */
      /* do loops on reflections in the box */
      for(;;)
      {
        /* get intersections for all box sides */
        double q, nx,ny,nz;
        double this_length;
        int side=0;

        bounces++;
        /* now look for intersection with guide sides and exit */
        side = Gravity_guide_Trace(&dt, &GVars, x, y, z,
            vx, vy, vz, w_chnum, nslit, h_chnum, nhslit,
            &nx, &ny, &nz);

        /* only positive dt are valid */
        /* exit reflection loops if no intersection (neutron is after box) */
        if (side == 0 || dt <= 0)
          { if (GVars.warnings < 100)
              fprintf(stderr,"%s: warning: neutron has entered guide, but can not exit !\n", GVars.compcurname);
            GVars.warnings++;
            x += w_adj; y += h_adj; ABSORB; } /* should never occur */

        /* propagate to dt */
        PROP_GRAV_DT(dt, GVars.gx, GVars.gy, GVars.gz);

        /* do reflection on speed for l/r/u/d sides */
        if (side == 5) /* neutron reaches end of guide: end loop and exit comp */
          { GVars.N_reflection[side]++; x += w_adj; y += h_adj; SCATTER; x -= w_adj; y -= h_adj; break; }
        /* else reflection on a guide wall */
        if(GVars.M[side] == 0 || Qc == 0 || R0 == 0)  /* walls are absorbing */
          { x += w_adj; y += h_adj; ABSORB; }
        /* handle chamfers */
        this_width = w1+(w2-w1)*z/l;
        this_height= h1+(h2-h1)*z/l;
        this_length= fmod(z, l/nelements);
        /* absorb on input/output of element parts */
        if (GVars.chamfer_z && (this_length<GVars.chamfer_z || this_length>l/nelements-GVars.chamfer_z))
        { x += w_adj; y += h_adj; ABSORB; }
        /* absorb on l/r/t/b sides */
        if (GVars.chamfer_lr && (side==1 || side==2) && (fabs(y+h_adj)>this_height/2-GVars.chamfer_lr))
        { x += w_adj; y += h_adj; ABSORB; }
        if (GVars.chamfer_tb && (side==3 || side==4) && (fabs(x+w_adj)>this_width/2- GVars.chamfer_tb))
        { x += w_adj; y += h_adj; ABSORB; }
        /* change/mirror velocity: h_f = v - n.2*n.v/|n|^2 */
        GVars.N_reflection[side]++; /* GVars.norm_n2 > 0 was checked at INIT */
        /* compute n.v using current values */
        B = scalar_prod(vx,vy,vz,nx,ny,nz);
        dt = 2*B/GVars.norm_n2[side]; /* 2*n.v/|n|^2 */
        vx -= nx*dt;
        vy -= ny*dt;
        vz -= nz*dt;

        /* compute q and modify neutron weight */
        /* scattering q=|n_i-n_f| = V2Q*|vf - v| = V2Q*2*n.v/|n| */
        q = 2*V2Q*fabs(B)/GVars.norm_n[side];

        if (reflect && strlen(reflect) && strcmp(reflect,"NULL") && strcmp(reflect,"0"))
          TableReflecFunc(q, &pTable, &B);
        else {
          double par[] = {R0, Qc, GVars.Alpha[side], GVars.M[side], W};
          StdReflecFunc(q, par, &B);
        }
        if (B <= 0) { x += w_adj; y += h_adj; ABSORB; }
        else p *= B;
        x += w_adj; y += h_adj; SCATTER; x -= w_adj; y -= h_adj;
        GVars.N_reflection[0]++;
        /* go to the next reflection */
        if (bounces > 1000) ABSORB;
      } /* end for */
      x += w_adj; y += h_adj; /* Re-adjust origin after SCATTER */
    }

    if (GVars.fc_freq != 0 || GVars.fc_phase != 0) { /* rotate back neutron w/r to guide element */
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

  } /* if l */
}
#line 16660 "./ISIS_SANS2d.c"
}   /* End of bender2=Guide_gravity() SETTING parameter declarations. */
#undef pTable
#undef GVars
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompbender2:
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

  /* TRACE Component bender3 [7] */
  mccoordschange(mcposrbender3, mcrotrbender3,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component bender3 (without coords transformations) */
  mcJumpTrace_bender3:
  SIG_MESSAGE("bender3 (Trace)");
  mcDEBUG_COMP("bender3")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbCompbender3
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
#define mccompcurname  bender3
#define mccompcurtype  Guide_gravity
#define mccompcurindex 7
#define GVars mccbender3_GVars
#define pTable mccbender3_pTable
{   /* Declarations of bender3=Guide_gravity() SETTING parameters. */
MCNUM w1 = mccbender3_w1;
MCNUM h1 = mccbender3_h1;
MCNUM w2 = mccbender3_w2;
MCNUM h2 = mccbender3_h2;
MCNUM l = mccbender3_l;
MCNUM R0 = mccbender3_R0;
MCNUM Qc = mccbender3_Qc;
MCNUM alpha = mccbender3_alpha;
MCNUM m = mccbender3_m;
MCNUM W = mccbender3_W;
MCNUM nslit = mccbender3_nslit;
MCNUM d = mccbender3_d;
MCNUM mleft = mccbender3_mleft;
MCNUM mright = mccbender3_mright;
MCNUM mtop = mccbender3_mtop;
MCNUM mbottom = mccbender3_mbottom;
MCNUM nhslit = mccbender3_nhslit;
MCNUM G = mccbender3_G;
MCNUM aleft = mccbender3_aleft;
MCNUM aright = mccbender3_aright;
MCNUM atop = mccbender3_atop;
MCNUM abottom = mccbender3_abottom;
MCNUM wavy = mccbender3_wavy;
MCNUM wavy_z = mccbender3_wavy_z;
MCNUM wavy_tb = mccbender3_wavy_tb;
MCNUM wavy_lr = mccbender3_wavy_lr;
MCNUM chamfers = mccbender3_chamfers;
MCNUM chamfers_z = mccbender3_chamfers_z;
MCNUM chamfers_lr = mccbender3_chamfers_lr;
MCNUM chamfers_tb = mccbender3_chamfers_tb;
MCNUM nelements = mccbender3_nelements;
MCNUM nu = mccbender3_nu;
MCNUM phase = mccbender3_phase;
char* reflect = mccbender3_reflect;
#line 392 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide_gravity.comp"
{
  if (l > 0 && nelements > 0) {
    double B, C, dt;
    int    ret, bounces = 0, i=0;
    double this_width, this_height;
    double angle=0;

    if (GVars.fc_freq != 0 || GVars.fc_phase != 0) { /* rotate neutron w/r to guide element */
      /* approximation of rotating straight Fermi Chopper */
      Coords   X = coords_set(x,y,z-l/2);  /* current coordinates of neutron in centered static frame */
      Rotation R;
      double dt=(-z+l/2)/vz; /* time shift to each center of slit package */
      angle=fmod(360*GVars.fc_freq*(t+dt)+GVars.fc_phase, 360); /* in deg */
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

    for (i=0; i<7; GVars.N_reflection[i++] = 0);

    /* propagate to box input (with gravitation) in comp local coords */
    /* A = 0.5 n.g; B = n.v; C = n.(r-W); */
    /* 0=Z0 side: n=(0, 0, -l) ; W = (0, 0, 0) (at z=0, guide input)*/
    B = -l*vz; C = -l*z;

    ret = solve_2nd_order(&dt, NULL, GVars.A[0], B, C);
    if (ret==0) ABSORB;

    if (dt>0.0) PROP_GRAV_DT(dt, GVars.gx, GVars.gy, GVars.gz); else if (angle) ABSORB;
    GVars.N_reflection[6]++;

    this_width  = w1;
    this_height = h1;

  /* check if we are in the box input, else absorb */
    if (fabs(x) > this_width/2 || fabs(y) > this_height/2)
      ABSORB;
    else
    {
      double w_edge, w_adj; /* Channel displacement on X */
      double h_edge, h_adj; /* Channel displacement on Y */
      double w_chnum,h_chnum; /* channel indexes */

      SCATTER;

      /* X: Shift origin to center of channel hit (absorb if hit dividing walls) */
      x += w1/2.0;
      w_chnum = floor(x/(GVars.w1c+d));  /* 0= right side, nslit+1=left side  */
      w_edge  = w_chnum*(GVars.w1c+d);
      if(x - w_edge > GVars.w1c)
      {
        x -= w1/2.0; /* Re-adjust origin */
        ABSORB;
      }
      w_adj = w_edge + (GVars.w1c)/2.0;
      x -= w_adj; w_adj -=  w1/2.0;

      /* Y: Shift origin to center of channel hit (absorb if hit dividing walls) */
      y += h1/2.0;
      h_chnum = floor(y/(GVars.h1c+d));  /* 0= lower side, nslit+1=upper side  */
      h_edge  = h_chnum*(GVars.h1c+d);
      if(y - h_edge > GVars.h1c)
      {
        y -= h1/2.0; /* Re-adjust origin */
        ABSORB;
      }
      h_adj = h_edge + (GVars.h1c)/2.0;
      y -= h_adj; h_adj -=  h1/2.0;

      /* neutron is now in the input window of the guide */
      /* do loops on reflections in the box */
      for(;;)
      {
        /* get intersections for all box sides */
        double q, nx,ny,nz;
        double this_length;
        int side=0;

        bounces++;
        /* now look for intersection with guide sides and exit */
        side = Gravity_guide_Trace(&dt, &GVars, x, y, z,
            vx, vy, vz, w_chnum, nslit, h_chnum, nhslit,
            &nx, &ny, &nz);

        /* only positive dt are valid */
        /* exit reflection loops if no intersection (neutron is after box) */
        if (side == 0 || dt <= 0)
          { if (GVars.warnings < 100)
              fprintf(stderr,"%s: warning: neutron has entered guide, but can not exit !\n", GVars.compcurname);
            GVars.warnings++;
            x += w_adj; y += h_adj; ABSORB; } /* should never occur */

        /* propagate to dt */
        PROP_GRAV_DT(dt, GVars.gx, GVars.gy, GVars.gz);

        /* do reflection on speed for l/r/u/d sides */
        if (side == 5) /* neutron reaches end of guide: end loop and exit comp */
          { GVars.N_reflection[side]++; x += w_adj; y += h_adj; SCATTER; x -= w_adj; y -= h_adj; break; }
        /* else reflection on a guide wall */
        if(GVars.M[side] == 0 || Qc == 0 || R0 == 0)  /* walls are absorbing */
          { x += w_adj; y += h_adj; ABSORB; }
        /* handle chamfers */
        this_width = w1+(w2-w1)*z/l;
        this_height= h1+(h2-h1)*z/l;
        this_length= fmod(z, l/nelements);
        /* absorb on input/output of element parts */
        if (GVars.chamfer_z && (this_length<GVars.chamfer_z || this_length>l/nelements-GVars.chamfer_z))
        { x += w_adj; y += h_adj; ABSORB; }
        /* absorb on l/r/t/b sides */
        if (GVars.chamfer_lr && (side==1 || side==2) && (fabs(y+h_adj)>this_height/2-GVars.chamfer_lr))
        { x += w_adj; y += h_adj; ABSORB; }
        if (GVars.chamfer_tb && (side==3 || side==4) && (fabs(x+w_adj)>this_width/2- GVars.chamfer_tb))
        { x += w_adj; y += h_adj; ABSORB; }
        /* change/mirror velocity: h_f = v - n.2*n.v/|n|^2 */
        GVars.N_reflection[side]++; /* GVars.norm_n2 > 0 was checked at INIT */
        /* compute n.v using current values */
        B = scalar_prod(vx,vy,vz,nx,ny,nz);
        dt = 2*B/GVars.norm_n2[side]; /* 2*n.v/|n|^2 */
        vx -= nx*dt;
        vy -= ny*dt;
        vz -= nz*dt;

        /* compute q and modify neutron weight */
        /* scattering q=|n_i-n_f| = V2Q*|vf - v| = V2Q*2*n.v/|n| */
        q = 2*V2Q*fabs(B)/GVars.norm_n[side];

        if (reflect && strlen(reflect) && strcmp(reflect,"NULL") && strcmp(reflect,"0"))
          TableReflecFunc(q, &pTable, &B);
        else {
          double par[] = {R0, Qc, GVars.Alpha[side], GVars.M[side], W};
          StdReflecFunc(q, par, &B);
        }
        if (B <= 0) { x += w_adj; y += h_adj; ABSORB; }
        else p *= B;
        x += w_adj; y += h_adj; SCATTER; x -= w_adj; y -= h_adj;
        GVars.N_reflection[0]++;
        /* go to the next reflection */
        if (bounces > 1000) ABSORB;
      } /* end for */
      x += w_adj; y += h_adj; /* Re-adjust origin after SCATTER */
    }

    if (GVars.fc_freq != 0 || GVars.fc_phase != 0) { /* rotate back neutron w/r to guide element */
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

  } /* if l */
}
#line 16973 "./ISIS_SANS2d.c"
}   /* End of bender3=Guide_gravity() SETTING parameter declarations. */
#undef pTable
#undef GVars
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompbender3:
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

  /* TRACE Component bender4 [8] */
  mccoordschange(mcposrbender4, mcrotrbender4,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component bender4 (without coords transformations) */
  mcJumpTrace_bender4:
  SIG_MESSAGE("bender4 (Trace)");
  mcDEBUG_COMP("bender4")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbCompbender4
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
#define mccompcurname  bender4
#define mccompcurtype  Guide_gravity
#define mccompcurindex 8
#define GVars mccbender4_GVars
#define pTable mccbender4_pTable
{   /* Declarations of bender4=Guide_gravity() SETTING parameters. */
MCNUM w1 = mccbender4_w1;
MCNUM h1 = mccbender4_h1;
MCNUM w2 = mccbender4_w2;
MCNUM h2 = mccbender4_h2;
MCNUM l = mccbender4_l;
MCNUM R0 = mccbender4_R0;
MCNUM Qc = mccbender4_Qc;
MCNUM alpha = mccbender4_alpha;
MCNUM m = mccbender4_m;
MCNUM W = mccbender4_W;
MCNUM nslit = mccbender4_nslit;
MCNUM d = mccbender4_d;
MCNUM mleft = mccbender4_mleft;
MCNUM mright = mccbender4_mright;
MCNUM mtop = mccbender4_mtop;
MCNUM mbottom = mccbender4_mbottom;
MCNUM nhslit = mccbender4_nhslit;
MCNUM G = mccbender4_G;
MCNUM aleft = mccbender4_aleft;
MCNUM aright = mccbender4_aright;
MCNUM atop = mccbender4_atop;
MCNUM abottom = mccbender4_abottom;
MCNUM wavy = mccbender4_wavy;
MCNUM wavy_z = mccbender4_wavy_z;
MCNUM wavy_tb = mccbender4_wavy_tb;
MCNUM wavy_lr = mccbender4_wavy_lr;
MCNUM chamfers = mccbender4_chamfers;
MCNUM chamfers_z = mccbender4_chamfers_z;
MCNUM chamfers_lr = mccbender4_chamfers_lr;
MCNUM chamfers_tb = mccbender4_chamfers_tb;
MCNUM nelements = mccbender4_nelements;
MCNUM nu = mccbender4_nu;
MCNUM phase = mccbender4_phase;
char* reflect = mccbender4_reflect;
#line 392 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide_gravity.comp"
{
  if (l > 0 && nelements > 0) {
    double B, C, dt;
    int    ret, bounces = 0, i=0;
    double this_width, this_height;
    double angle=0;

    if (GVars.fc_freq != 0 || GVars.fc_phase != 0) { /* rotate neutron w/r to guide element */
      /* approximation of rotating straight Fermi Chopper */
      Coords   X = coords_set(x,y,z-l/2);  /* current coordinates of neutron in centered static frame */
      Rotation R;
      double dt=(-z+l/2)/vz; /* time shift to each center of slit package */
      angle=fmod(360*GVars.fc_freq*(t+dt)+GVars.fc_phase, 360); /* in deg */
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

    for (i=0; i<7; GVars.N_reflection[i++] = 0);

    /* propagate to box input (with gravitation) in comp local coords */
    /* A = 0.5 n.g; B = n.v; C = n.(r-W); */
    /* 0=Z0 side: n=(0, 0, -l) ; W = (0, 0, 0) (at z=0, guide input)*/
    B = -l*vz; C = -l*z;

    ret = solve_2nd_order(&dt, NULL, GVars.A[0], B, C);
    if (ret==0) ABSORB;

    if (dt>0.0) PROP_GRAV_DT(dt, GVars.gx, GVars.gy, GVars.gz); else if (angle) ABSORB;
    GVars.N_reflection[6]++;

    this_width  = w1;
    this_height = h1;

  /* check if we are in the box input, else absorb */
    if (fabs(x) > this_width/2 || fabs(y) > this_height/2)
      ABSORB;
    else
    {
      double w_edge, w_adj; /* Channel displacement on X */
      double h_edge, h_adj; /* Channel displacement on Y */
      double w_chnum,h_chnum; /* channel indexes */

      SCATTER;

      /* X: Shift origin to center of channel hit (absorb if hit dividing walls) */
      x += w1/2.0;
      w_chnum = floor(x/(GVars.w1c+d));  /* 0= right side, nslit+1=left side  */
      w_edge  = w_chnum*(GVars.w1c+d);
      if(x - w_edge > GVars.w1c)
      {
        x -= w1/2.0; /* Re-adjust origin */
        ABSORB;
      }
      w_adj = w_edge + (GVars.w1c)/2.0;
      x -= w_adj; w_adj -=  w1/2.0;

      /* Y: Shift origin to center of channel hit (absorb if hit dividing walls) */
      y += h1/2.0;
      h_chnum = floor(y/(GVars.h1c+d));  /* 0= lower side, nslit+1=upper side  */
      h_edge  = h_chnum*(GVars.h1c+d);
      if(y - h_edge > GVars.h1c)
      {
        y -= h1/2.0; /* Re-adjust origin */
        ABSORB;
      }
      h_adj = h_edge + (GVars.h1c)/2.0;
      y -= h_adj; h_adj -=  h1/2.0;

      /* neutron is now in the input window of the guide */
      /* do loops on reflections in the box */
      for(;;)
      {
        /* get intersections for all box sides */
        double q, nx,ny,nz;
        double this_length;
        int side=0;

        bounces++;
        /* now look for intersection with guide sides and exit */
        side = Gravity_guide_Trace(&dt, &GVars, x, y, z,
            vx, vy, vz, w_chnum, nslit, h_chnum, nhslit,
            &nx, &ny, &nz);

        /* only positive dt are valid */
        /* exit reflection loops if no intersection (neutron is after box) */
        if (side == 0 || dt <= 0)
          { if (GVars.warnings < 100)
              fprintf(stderr,"%s: warning: neutron has entered guide, but can not exit !\n", GVars.compcurname);
            GVars.warnings++;
            x += w_adj; y += h_adj; ABSORB; } /* should never occur */

        /* propagate to dt */
        PROP_GRAV_DT(dt, GVars.gx, GVars.gy, GVars.gz);

        /* do reflection on speed for l/r/u/d sides */
        if (side == 5) /* neutron reaches end of guide: end loop and exit comp */
          { GVars.N_reflection[side]++; x += w_adj; y += h_adj; SCATTER; x -= w_adj; y -= h_adj; break; }
        /* else reflection on a guide wall */
        if(GVars.M[side] == 0 || Qc == 0 || R0 == 0)  /* walls are absorbing */
          { x += w_adj; y += h_adj; ABSORB; }
        /* handle chamfers */
        this_width = w1+(w2-w1)*z/l;
        this_height= h1+(h2-h1)*z/l;
        this_length= fmod(z, l/nelements);
        /* absorb on input/output of element parts */
        if (GVars.chamfer_z && (this_length<GVars.chamfer_z || this_length>l/nelements-GVars.chamfer_z))
        { x += w_adj; y += h_adj; ABSORB; }
        /* absorb on l/r/t/b sides */
        if (GVars.chamfer_lr && (side==1 || side==2) && (fabs(y+h_adj)>this_height/2-GVars.chamfer_lr))
        { x += w_adj; y += h_adj; ABSORB; }
        if (GVars.chamfer_tb && (side==3 || side==4) && (fabs(x+w_adj)>this_width/2- GVars.chamfer_tb))
        { x += w_adj; y += h_adj; ABSORB; }
        /* change/mirror velocity: h_f = v - n.2*n.v/|n|^2 */
        GVars.N_reflection[side]++; /* GVars.norm_n2 > 0 was checked at INIT */
        /* compute n.v using current values */
        B = scalar_prod(vx,vy,vz,nx,ny,nz);
        dt = 2*B/GVars.norm_n2[side]; /* 2*n.v/|n|^2 */
        vx -= nx*dt;
        vy -= ny*dt;
        vz -= nz*dt;

        /* compute q and modify neutron weight */
        /* scattering q=|n_i-n_f| = V2Q*|vf - v| = V2Q*2*n.v/|n| */
        q = 2*V2Q*fabs(B)/GVars.norm_n[side];

        if (reflect && strlen(reflect) && strcmp(reflect,"NULL") && strcmp(reflect,"0"))
          TableReflecFunc(q, &pTable, &B);
        else {
          double par[] = {R0, Qc, GVars.Alpha[side], GVars.M[side], W};
          StdReflecFunc(q, par, &B);
        }
        if (B <= 0) { x += w_adj; y += h_adj; ABSORB; }
        else p *= B;
        x += w_adj; y += h_adj; SCATTER; x -= w_adj; y -= h_adj;
        GVars.N_reflection[0]++;
        /* go to the next reflection */
        if (bounces > 1000) ABSORB;
      } /* end for */
      x += w_adj; y += h_adj; /* Re-adjust origin after SCATTER */
    }

    if (GVars.fc_freq != 0 || GVars.fc_phase != 0) { /* rotate back neutron w/r to guide element */
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

  } /* if l */
}
#line 17286 "./ISIS_SANS2d.c"
}   /* End of bender4=Guide_gravity() SETTING parameter declarations. */
#undef pTable
#undef GVars
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompbender4:
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

  /* TRACE Component bender5 [9] */
  mccoordschange(mcposrbender5, mcrotrbender5,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component bender5 (without coords transformations) */
  mcJumpTrace_bender5:
  SIG_MESSAGE("bender5 (Trace)");
  mcDEBUG_COMP("bender5")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbCompbender5
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
#define mccompcurname  bender5
#define mccompcurtype  Guide_gravity
#define mccompcurindex 9
#define GVars mccbender5_GVars
#define pTable mccbender5_pTable
{   /* Declarations of bender5=Guide_gravity() SETTING parameters. */
MCNUM w1 = mccbender5_w1;
MCNUM h1 = mccbender5_h1;
MCNUM w2 = mccbender5_w2;
MCNUM h2 = mccbender5_h2;
MCNUM l = mccbender5_l;
MCNUM R0 = mccbender5_R0;
MCNUM Qc = mccbender5_Qc;
MCNUM alpha = mccbender5_alpha;
MCNUM m = mccbender5_m;
MCNUM W = mccbender5_W;
MCNUM nslit = mccbender5_nslit;
MCNUM d = mccbender5_d;
MCNUM mleft = mccbender5_mleft;
MCNUM mright = mccbender5_mright;
MCNUM mtop = mccbender5_mtop;
MCNUM mbottom = mccbender5_mbottom;
MCNUM nhslit = mccbender5_nhslit;
MCNUM G = mccbender5_G;
MCNUM aleft = mccbender5_aleft;
MCNUM aright = mccbender5_aright;
MCNUM atop = mccbender5_atop;
MCNUM abottom = mccbender5_abottom;
MCNUM wavy = mccbender5_wavy;
MCNUM wavy_z = mccbender5_wavy_z;
MCNUM wavy_tb = mccbender5_wavy_tb;
MCNUM wavy_lr = mccbender5_wavy_lr;
MCNUM chamfers = mccbender5_chamfers;
MCNUM chamfers_z = mccbender5_chamfers_z;
MCNUM chamfers_lr = mccbender5_chamfers_lr;
MCNUM chamfers_tb = mccbender5_chamfers_tb;
MCNUM nelements = mccbender5_nelements;
MCNUM nu = mccbender5_nu;
MCNUM phase = mccbender5_phase;
char* reflect = mccbender5_reflect;
#line 392 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide_gravity.comp"
{
  if (l > 0 && nelements > 0) {
    double B, C, dt;
    int    ret, bounces = 0, i=0;
    double this_width, this_height;
    double angle=0;

    if (GVars.fc_freq != 0 || GVars.fc_phase != 0) { /* rotate neutron w/r to guide element */
      /* approximation of rotating straight Fermi Chopper */
      Coords   X = coords_set(x,y,z-l/2);  /* current coordinates of neutron in centered static frame */
      Rotation R;
      double dt=(-z+l/2)/vz; /* time shift to each center of slit package */
      angle=fmod(360*GVars.fc_freq*(t+dt)+GVars.fc_phase, 360); /* in deg */
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

    for (i=0; i<7; GVars.N_reflection[i++] = 0);

    /* propagate to box input (with gravitation) in comp local coords */
    /* A = 0.5 n.g; B = n.v; C = n.(r-W); */
    /* 0=Z0 side: n=(0, 0, -l) ; W = (0, 0, 0) (at z=0, guide input)*/
    B = -l*vz; C = -l*z;

    ret = solve_2nd_order(&dt, NULL, GVars.A[0], B, C);
    if (ret==0) ABSORB;

    if (dt>0.0) PROP_GRAV_DT(dt, GVars.gx, GVars.gy, GVars.gz); else if (angle) ABSORB;
    GVars.N_reflection[6]++;

    this_width  = w1;
    this_height = h1;

  /* check if we are in the box input, else absorb */
    if (fabs(x) > this_width/2 || fabs(y) > this_height/2)
      ABSORB;
    else
    {
      double w_edge, w_adj; /* Channel displacement on X */
      double h_edge, h_adj; /* Channel displacement on Y */
      double w_chnum,h_chnum; /* channel indexes */

      SCATTER;

      /* X: Shift origin to center of channel hit (absorb if hit dividing walls) */
      x += w1/2.0;
      w_chnum = floor(x/(GVars.w1c+d));  /* 0= right side, nslit+1=left side  */
      w_edge  = w_chnum*(GVars.w1c+d);
      if(x - w_edge > GVars.w1c)
      {
        x -= w1/2.0; /* Re-adjust origin */
        ABSORB;
      }
      w_adj = w_edge + (GVars.w1c)/2.0;
      x -= w_adj; w_adj -=  w1/2.0;

      /* Y: Shift origin to center of channel hit (absorb if hit dividing walls) */
      y += h1/2.0;
      h_chnum = floor(y/(GVars.h1c+d));  /* 0= lower side, nslit+1=upper side  */
      h_edge  = h_chnum*(GVars.h1c+d);
      if(y - h_edge > GVars.h1c)
      {
        y -= h1/2.0; /* Re-adjust origin */
        ABSORB;
      }
      h_adj = h_edge + (GVars.h1c)/2.0;
      y -= h_adj; h_adj -=  h1/2.0;

      /* neutron is now in the input window of the guide */
      /* do loops on reflections in the box */
      for(;;)
      {
        /* get intersections for all box sides */
        double q, nx,ny,nz;
        double this_length;
        int side=0;

        bounces++;
        /* now look for intersection with guide sides and exit */
        side = Gravity_guide_Trace(&dt, &GVars, x, y, z,
            vx, vy, vz, w_chnum, nslit, h_chnum, nhslit,
            &nx, &ny, &nz);

        /* only positive dt are valid */
        /* exit reflection loops if no intersection (neutron is after box) */
        if (side == 0 || dt <= 0)
          { if (GVars.warnings < 100)
              fprintf(stderr,"%s: warning: neutron has entered guide, but can not exit !\n", GVars.compcurname);
            GVars.warnings++;
            x += w_adj; y += h_adj; ABSORB; } /* should never occur */

        /* propagate to dt */
        PROP_GRAV_DT(dt, GVars.gx, GVars.gy, GVars.gz);

        /* do reflection on speed for l/r/u/d sides */
        if (side == 5) /* neutron reaches end of guide: end loop and exit comp */
          { GVars.N_reflection[side]++; x += w_adj; y += h_adj; SCATTER; x -= w_adj; y -= h_adj; break; }
        /* else reflection on a guide wall */
        if(GVars.M[side] == 0 || Qc == 0 || R0 == 0)  /* walls are absorbing */
          { x += w_adj; y += h_adj; ABSORB; }
        /* handle chamfers */
        this_width = w1+(w2-w1)*z/l;
        this_height= h1+(h2-h1)*z/l;
        this_length= fmod(z, l/nelements);
        /* absorb on input/output of element parts */
        if (GVars.chamfer_z && (this_length<GVars.chamfer_z || this_length>l/nelements-GVars.chamfer_z))
        { x += w_adj; y += h_adj; ABSORB; }
        /* absorb on l/r/t/b sides */
        if (GVars.chamfer_lr && (side==1 || side==2) && (fabs(y+h_adj)>this_height/2-GVars.chamfer_lr))
        { x += w_adj; y += h_adj; ABSORB; }
        if (GVars.chamfer_tb && (side==3 || side==4) && (fabs(x+w_adj)>this_width/2- GVars.chamfer_tb))
        { x += w_adj; y += h_adj; ABSORB; }
        /* change/mirror velocity: h_f = v - n.2*n.v/|n|^2 */
        GVars.N_reflection[side]++; /* GVars.norm_n2 > 0 was checked at INIT */
        /* compute n.v using current values */
        B = scalar_prod(vx,vy,vz,nx,ny,nz);
        dt = 2*B/GVars.norm_n2[side]; /* 2*n.v/|n|^2 */
        vx -= nx*dt;
        vy -= ny*dt;
        vz -= nz*dt;

        /* compute q and modify neutron weight */
        /* scattering q=|n_i-n_f| = V2Q*|vf - v| = V2Q*2*n.v/|n| */
        q = 2*V2Q*fabs(B)/GVars.norm_n[side];

        if (reflect && strlen(reflect) && strcmp(reflect,"NULL") && strcmp(reflect,"0"))
          TableReflecFunc(q, &pTable, &B);
        else {
          double par[] = {R0, Qc, GVars.Alpha[side], GVars.M[side], W};
          StdReflecFunc(q, par, &B);
        }
        if (B <= 0) { x += w_adj; y += h_adj; ABSORB; }
        else p *= B;
        x += w_adj; y += h_adj; SCATTER; x -= w_adj; y -= h_adj;
        GVars.N_reflection[0]++;
        /* go to the next reflection */
        if (bounces > 1000) ABSORB;
      } /* end for */
      x += w_adj; y += h_adj; /* Re-adjust origin after SCATTER */
    }

    if (GVars.fc_freq != 0 || GVars.fc_phase != 0) { /* rotate back neutron w/r to guide element */
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

  } /* if l */
}
#line 17599 "./ISIS_SANS2d.c"
}   /* End of bender5=Guide_gravity() SETTING parameter declarations. */
#undef pTable
#undef GVars
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompbender5:
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

  /* TRACE Component bender6 [10] */
  mccoordschange(mcposrbender6, mcrotrbender6,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component bender6 (without coords transformations) */
  mcJumpTrace_bender6:
  SIG_MESSAGE("bender6 (Trace)");
  mcDEBUG_COMP("bender6")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbCompbender6
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
#define mccompcurname  bender6
#define mccompcurtype  Guide_gravity
#define mccompcurindex 10
#define GVars mccbender6_GVars
#define pTable mccbender6_pTable
{   /* Declarations of bender6=Guide_gravity() SETTING parameters. */
MCNUM w1 = mccbender6_w1;
MCNUM h1 = mccbender6_h1;
MCNUM w2 = mccbender6_w2;
MCNUM h2 = mccbender6_h2;
MCNUM l = mccbender6_l;
MCNUM R0 = mccbender6_R0;
MCNUM Qc = mccbender6_Qc;
MCNUM alpha = mccbender6_alpha;
MCNUM m = mccbender6_m;
MCNUM W = mccbender6_W;
MCNUM nslit = mccbender6_nslit;
MCNUM d = mccbender6_d;
MCNUM mleft = mccbender6_mleft;
MCNUM mright = mccbender6_mright;
MCNUM mtop = mccbender6_mtop;
MCNUM mbottom = mccbender6_mbottom;
MCNUM nhslit = mccbender6_nhslit;
MCNUM G = mccbender6_G;
MCNUM aleft = mccbender6_aleft;
MCNUM aright = mccbender6_aright;
MCNUM atop = mccbender6_atop;
MCNUM abottom = mccbender6_abottom;
MCNUM wavy = mccbender6_wavy;
MCNUM wavy_z = mccbender6_wavy_z;
MCNUM wavy_tb = mccbender6_wavy_tb;
MCNUM wavy_lr = mccbender6_wavy_lr;
MCNUM chamfers = mccbender6_chamfers;
MCNUM chamfers_z = mccbender6_chamfers_z;
MCNUM chamfers_lr = mccbender6_chamfers_lr;
MCNUM chamfers_tb = mccbender6_chamfers_tb;
MCNUM nelements = mccbender6_nelements;
MCNUM nu = mccbender6_nu;
MCNUM phase = mccbender6_phase;
char* reflect = mccbender6_reflect;
#line 392 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide_gravity.comp"
{
  if (l > 0 && nelements > 0) {
    double B, C, dt;
    int    ret, bounces = 0, i=0;
    double this_width, this_height;
    double angle=0;

    if (GVars.fc_freq != 0 || GVars.fc_phase != 0) { /* rotate neutron w/r to guide element */
      /* approximation of rotating straight Fermi Chopper */
      Coords   X = coords_set(x,y,z-l/2);  /* current coordinates of neutron in centered static frame */
      Rotation R;
      double dt=(-z+l/2)/vz; /* time shift to each center of slit package */
      angle=fmod(360*GVars.fc_freq*(t+dt)+GVars.fc_phase, 360); /* in deg */
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

    for (i=0; i<7; GVars.N_reflection[i++] = 0);

    /* propagate to box input (with gravitation) in comp local coords */
    /* A = 0.5 n.g; B = n.v; C = n.(r-W); */
    /* 0=Z0 side: n=(0, 0, -l) ; W = (0, 0, 0) (at z=0, guide input)*/
    B = -l*vz; C = -l*z;

    ret = solve_2nd_order(&dt, NULL, GVars.A[0], B, C);
    if (ret==0) ABSORB;

    if (dt>0.0) PROP_GRAV_DT(dt, GVars.gx, GVars.gy, GVars.gz); else if (angle) ABSORB;
    GVars.N_reflection[6]++;

    this_width  = w1;
    this_height = h1;

  /* check if we are in the box input, else absorb */
    if (fabs(x) > this_width/2 || fabs(y) > this_height/2)
      ABSORB;
    else
    {
      double w_edge, w_adj; /* Channel displacement on X */
      double h_edge, h_adj; /* Channel displacement on Y */
      double w_chnum,h_chnum; /* channel indexes */

      SCATTER;

      /* X: Shift origin to center of channel hit (absorb if hit dividing walls) */
      x += w1/2.0;
      w_chnum = floor(x/(GVars.w1c+d));  /* 0= right side, nslit+1=left side  */
      w_edge  = w_chnum*(GVars.w1c+d);
      if(x - w_edge > GVars.w1c)
      {
        x -= w1/2.0; /* Re-adjust origin */
        ABSORB;
      }
      w_adj = w_edge + (GVars.w1c)/2.0;
      x -= w_adj; w_adj -=  w1/2.0;

      /* Y: Shift origin to center of channel hit (absorb if hit dividing walls) */
      y += h1/2.0;
      h_chnum = floor(y/(GVars.h1c+d));  /* 0= lower side, nslit+1=upper side  */
      h_edge  = h_chnum*(GVars.h1c+d);
      if(y - h_edge > GVars.h1c)
      {
        y -= h1/2.0; /* Re-adjust origin */
        ABSORB;
      }
      h_adj = h_edge + (GVars.h1c)/2.0;
      y -= h_adj; h_adj -=  h1/2.0;

      /* neutron is now in the input window of the guide */
      /* do loops on reflections in the box */
      for(;;)
      {
        /* get intersections for all box sides */
        double q, nx,ny,nz;
        double this_length;
        int side=0;

        bounces++;
        /* now look for intersection with guide sides and exit */
        side = Gravity_guide_Trace(&dt, &GVars, x, y, z,
            vx, vy, vz, w_chnum, nslit, h_chnum, nhslit,
            &nx, &ny, &nz);

        /* only positive dt are valid */
        /* exit reflection loops if no intersection (neutron is after box) */
        if (side == 0 || dt <= 0)
          { if (GVars.warnings < 100)
              fprintf(stderr,"%s: warning: neutron has entered guide, but can not exit !\n", GVars.compcurname);
            GVars.warnings++;
            x += w_adj; y += h_adj; ABSORB; } /* should never occur */

        /* propagate to dt */
        PROP_GRAV_DT(dt, GVars.gx, GVars.gy, GVars.gz);

        /* do reflection on speed for l/r/u/d sides */
        if (side == 5) /* neutron reaches end of guide: end loop and exit comp */
          { GVars.N_reflection[side]++; x += w_adj; y += h_adj; SCATTER; x -= w_adj; y -= h_adj; break; }
        /* else reflection on a guide wall */
        if(GVars.M[side] == 0 || Qc == 0 || R0 == 0)  /* walls are absorbing */
          { x += w_adj; y += h_adj; ABSORB; }
        /* handle chamfers */
        this_width = w1+(w2-w1)*z/l;
        this_height= h1+(h2-h1)*z/l;
        this_length= fmod(z, l/nelements);
        /* absorb on input/output of element parts */
        if (GVars.chamfer_z && (this_length<GVars.chamfer_z || this_length>l/nelements-GVars.chamfer_z))
        { x += w_adj; y += h_adj; ABSORB; }
        /* absorb on l/r/t/b sides */
        if (GVars.chamfer_lr && (side==1 || side==2) && (fabs(y+h_adj)>this_height/2-GVars.chamfer_lr))
        { x += w_adj; y += h_adj; ABSORB; }
        if (GVars.chamfer_tb && (side==3 || side==4) && (fabs(x+w_adj)>this_width/2- GVars.chamfer_tb))
        { x += w_adj; y += h_adj; ABSORB; }
        /* change/mirror velocity: h_f = v - n.2*n.v/|n|^2 */
        GVars.N_reflection[side]++; /* GVars.norm_n2 > 0 was checked at INIT */
        /* compute n.v using current values */
        B = scalar_prod(vx,vy,vz,nx,ny,nz);
        dt = 2*B/GVars.norm_n2[side]; /* 2*n.v/|n|^2 */
        vx -= nx*dt;
        vy -= ny*dt;
        vz -= nz*dt;

        /* compute q and modify neutron weight */
        /* scattering q=|n_i-n_f| = V2Q*|vf - v| = V2Q*2*n.v/|n| */
        q = 2*V2Q*fabs(B)/GVars.norm_n[side];

        if (reflect && strlen(reflect) && strcmp(reflect,"NULL") && strcmp(reflect,"0"))
          TableReflecFunc(q, &pTable, &B);
        else {
          double par[] = {R0, Qc, GVars.Alpha[side], GVars.M[side], W};
          StdReflecFunc(q, par, &B);
        }
        if (B <= 0) { x += w_adj; y += h_adj; ABSORB; }
        else p *= B;
        x += w_adj; y += h_adj; SCATTER; x -= w_adj; y -= h_adj;
        GVars.N_reflection[0]++;
        /* go to the next reflection */
        if (bounces > 1000) ABSORB;
      } /* end for */
      x += w_adj; y += h_adj; /* Re-adjust origin after SCATTER */
    }

    if (GVars.fc_freq != 0 || GVars.fc_phase != 0) { /* rotate back neutron w/r to guide element */
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

  } /* if l */
}
#line 17912 "./ISIS_SANS2d.c"
}   /* End of bender6=Guide_gravity() SETTING parameter declarations. */
#undef pTable
#undef GVars
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompbender6:
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

  /* TRACE Component bender7 [11] */
  mccoordschange(mcposrbender7, mcrotrbender7,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component bender7 (without coords transformations) */
  mcJumpTrace_bender7:
  SIG_MESSAGE("bender7 (Trace)");
  mcDEBUG_COMP("bender7")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbCompbender7
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
#define mccompcurname  bender7
#define mccompcurtype  Guide_gravity
#define mccompcurindex 11
#define GVars mccbender7_GVars
#define pTable mccbender7_pTable
{   /* Declarations of bender7=Guide_gravity() SETTING parameters. */
MCNUM w1 = mccbender7_w1;
MCNUM h1 = mccbender7_h1;
MCNUM w2 = mccbender7_w2;
MCNUM h2 = mccbender7_h2;
MCNUM l = mccbender7_l;
MCNUM R0 = mccbender7_R0;
MCNUM Qc = mccbender7_Qc;
MCNUM alpha = mccbender7_alpha;
MCNUM m = mccbender7_m;
MCNUM W = mccbender7_W;
MCNUM nslit = mccbender7_nslit;
MCNUM d = mccbender7_d;
MCNUM mleft = mccbender7_mleft;
MCNUM mright = mccbender7_mright;
MCNUM mtop = mccbender7_mtop;
MCNUM mbottom = mccbender7_mbottom;
MCNUM nhslit = mccbender7_nhslit;
MCNUM G = mccbender7_G;
MCNUM aleft = mccbender7_aleft;
MCNUM aright = mccbender7_aright;
MCNUM atop = mccbender7_atop;
MCNUM abottom = mccbender7_abottom;
MCNUM wavy = mccbender7_wavy;
MCNUM wavy_z = mccbender7_wavy_z;
MCNUM wavy_tb = mccbender7_wavy_tb;
MCNUM wavy_lr = mccbender7_wavy_lr;
MCNUM chamfers = mccbender7_chamfers;
MCNUM chamfers_z = mccbender7_chamfers_z;
MCNUM chamfers_lr = mccbender7_chamfers_lr;
MCNUM chamfers_tb = mccbender7_chamfers_tb;
MCNUM nelements = mccbender7_nelements;
MCNUM nu = mccbender7_nu;
MCNUM phase = mccbender7_phase;
char* reflect = mccbender7_reflect;
#line 392 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide_gravity.comp"
{
  if (l > 0 && nelements > 0) {
    double B, C, dt;
    int    ret, bounces = 0, i=0;
    double this_width, this_height;
    double angle=0;

    if (GVars.fc_freq != 0 || GVars.fc_phase != 0) { /* rotate neutron w/r to guide element */
      /* approximation of rotating straight Fermi Chopper */
      Coords   X = coords_set(x,y,z-l/2);  /* current coordinates of neutron in centered static frame */
      Rotation R;
      double dt=(-z+l/2)/vz; /* time shift to each center of slit package */
      angle=fmod(360*GVars.fc_freq*(t+dt)+GVars.fc_phase, 360); /* in deg */
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

    for (i=0; i<7; GVars.N_reflection[i++] = 0);

    /* propagate to box input (with gravitation) in comp local coords */
    /* A = 0.5 n.g; B = n.v; C = n.(r-W); */
    /* 0=Z0 side: n=(0, 0, -l) ; W = (0, 0, 0) (at z=0, guide input)*/
    B = -l*vz; C = -l*z;

    ret = solve_2nd_order(&dt, NULL, GVars.A[0], B, C);
    if (ret==0) ABSORB;

    if (dt>0.0) PROP_GRAV_DT(dt, GVars.gx, GVars.gy, GVars.gz); else if (angle) ABSORB;
    GVars.N_reflection[6]++;

    this_width  = w1;
    this_height = h1;

  /* check if we are in the box input, else absorb */
    if (fabs(x) > this_width/2 || fabs(y) > this_height/2)
      ABSORB;
    else
    {
      double w_edge, w_adj; /* Channel displacement on X */
      double h_edge, h_adj; /* Channel displacement on Y */
      double w_chnum,h_chnum; /* channel indexes */

      SCATTER;

      /* X: Shift origin to center of channel hit (absorb if hit dividing walls) */
      x += w1/2.0;
      w_chnum = floor(x/(GVars.w1c+d));  /* 0= right side, nslit+1=left side  */
      w_edge  = w_chnum*(GVars.w1c+d);
      if(x - w_edge > GVars.w1c)
      {
        x -= w1/2.0; /* Re-adjust origin */
        ABSORB;
      }
      w_adj = w_edge + (GVars.w1c)/2.0;
      x -= w_adj; w_adj -=  w1/2.0;

      /* Y: Shift origin to center of channel hit (absorb if hit dividing walls) */
      y += h1/2.0;
      h_chnum = floor(y/(GVars.h1c+d));  /* 0= lower side, nslit+1=upper side  */
      h_edge  = h_chnum*(GVars.h1c+d);
      if(y - h_edge > GVars.h1c)
      {
        y -= h1/2.0; /* Re-adjust origin */
        ABSORB;
      }
      h_adj = h_edge + (GVars.h1c)/2.0;
      y -= h_adj; h_adj -=  h1/2.0;

      /* neutron is now in the input window of the guide */
      /* do loops on reflections in the box */
      for(;;)
      {
        /* get intersections for all box sides */
        double q, nx,ny,nz;
        double this_length;
        int side=0;

        bounces++;
        /* now look for intersection with guide sides and exit */
        side = Gravity_guide_Trace(&dt, &GVars, x, y, z,
            vx, vy, vz, w_chnum, nslit, h_chnum, nhslit,
            &nx, &ny, &nz);

        /* only positive dt are valid */
        /* exit reflection loops if no intersection (neutron is after box) */
        if (side == 0 || dt <= 0)
          { if (GVars.warnings < 100)
              fprintf(stderr,"%s: warning: neutron has entered guide, but can not exit !\n", GVars.compcurname);
            GVars.warnings++;
            x += w_adj; y += h_adj; ABSORB; } /* should never occur */

        /* propagate to dt */
        PROP_GRAV_DT(dt, GVars.gx, GVars.gy, GVars.gz);

        /* do reflection on speed for l/r/u/d sides */
        if (side == 5) /* neutron reaches end of guide: end loop and exit comp */
          { GVars.N_reflection[side]++; x += w_adj; y += h_adj; SCATTER; x -= w_adj; y -= h_adj; break; }
        /* else reflection on a guide wall */
        if(GVars.M[side] == 0 || Qc == 0 || R0 == 0)  /* walls are absorbing */
          { x += w_adj; y += h_adj; ABSORB; }
        /* handle chamfers */
        this_width = w1+(w2-w1)*z/l;
        this_height= h1+(h2-h1)*z/l;
        this_length= fmod(z, l/nelements);
        /* absorb on input/output of element parts */
        if (GVars.chamfer_z && (this_length<GVars.chamfer_z || this_length>l/nelements-GVars.chamfer_z))
        { x += w_adj; y += h_adj; ABSORB; }
        /* absorb on l/r/t/b sides */
        if (GVars.chamfer_lr && (side==1 || side==2) && (fabs(y+h_adj)>this_height/2-GVars.chamfer_lr))
        { x += w_adj; y += h_adj; ABSORB; }
        if (GVars.chamfer_tb && (side==3 || side==4) && (fabs(x+w_adj)>this_width/2- GVars.chamfer_tb))
        { x += w_adj; y += h_adj; ABSORB; }
        /* change/mirror velocity: h_f = v - n.2*n.v/|n|^2 */
        GVars.N_reflection[side]++; /* GVars.norm_n2 > 0 was checked at INIT */
        /* compute n.v using current values */
        B = scalar_prod(vx,vy,vz,nx,ny,nz);
        dt = 2*B/GVars.norm_n2[side]; /* 2*n.v/|n|^2 */
        vx -= nx*dt;
        vy -= ny*dt;
        vz -= nz*dt;

        /* compute q and modify neutron weight */
        /* scattering q=|n_i-n_f| = V2Q*|vf - v| = V2Q*2*n.v/|n| */
        q = 2*V2Q*fabs(B)/GVars.norm_n[side];

        if (reflect && strlen(reflect) && strcmp(reflect,"NULL") && strcmp(reflect,"0"))
          TableReflecFunc(q, &pTable, &B);
        else {
          double par[] = {R0, Qc, GVars.Alpha[side], GVars.M[side], W};
          StdReflecFunc(q, par, &B);
        }
        if (B <= 0) { x += w_adj; y += h_adj; ABSORB; }
        else p *= B;
        x += w_adj; y += h_adj; SCATTER; x -= w_adj; y -= h_adj;
        GVars.N_reflection[0]++;
        /* go to the next reflection */
        if (bounces > 1000) ABSORB;
      } /* end for */
      x += w_adj; y += h_adj; /* Re-adjust origin after SCATTER */
    }

    if (GVars.fc_freq != 0 || GVars.fc_phase != 0) { /* rotate back neutron w/r to guide element */
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

  } /* if l */
}
#line 18225 "./ISIS_SANS2d.c"
}   /* End of bender7=Guide_gravity() SETTING parameter declarations. */
#undef pTable
#undef GVars
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompbender7:
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

  /* TRACE Component bender8 [12] */
  mccoordschange(mcposrbender8, mcrotrbender8,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component bender8 (without coords transformations) */
  mcJumpTrace_bender8:
  SIG_MESSAGE("bender8 (Trace)");
  mcDEBUG_COMP("bender8")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbCompbender8
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
#define mccompcurname  bender8
#define mccompcurtype  Guide_gravity
#define mccompcurindex 12
#define GVars mccbender8_GVars
#define pTable mccbender8_pTable
{   /* Declarations of bender8=Guide_gravity() SETTING parameters. */
MCNUM w1 = mccbender8_w1;
MCNUM h1 = mccbender8_h1;
MCNUM w2 = mccbender8_w2;
MCNUM h2 = mccbender8_h2;
MCNUM l = mccbender8_l;
MCNUM R0 = mccbender8_R0;
MCNUM Qc = mccbender8_Qc;
MCNUM alpha = mccbender8_alpha;
MCNUM m = mccbender8_m;
MCNUM W = mccbender8_W;
MCNUM nslit = mccbender8_nslit;
MCNUM d = mccbender8_d;
MCNUM mleft = mccbender8_mleft;
MCNUM mright = mccbender8_mright;
MCNUM mtop = mccbender8_mtop;
MCNUM mbottom = mccbender8_mbottom;
MCNUM nhslit = mccbender8_nhslit;
MCNUM G = mccbender8_G;
MCNUM aleft = mccbender8_aleft;
MCNUM aright = mccbender8_aright;
MCNUM atop = mccbender8_atop;
MCNUM abottom = mccbender8_abottom;
MCNUM wavy = mccbender8_wavy;
MCNUM wavy_z = mccbender8_wavy_z;
MCNUM wavy_tb = mccbender8_wavy_tb;
MCNUM wavy_lr = mccbender8_wavy_lr;
MCNUM chamfers = mccbender8_chamfers;
MCNUM chamfers_z = mccbender8_chamfers_z;
MCNUM chamfers_lr = mccbender8_chamfers_lr;
MCNUM chamfers_tb = mccbender8_chamfers_tb;
MCNUM nelements = mccbender8_nelements;
MCNUM nu = mccbender8_nu;
MCNUM phase = mccbender8_phase;
char* reflect = mccbender8_reflect;
#line 392 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide_gravity.comp"
{
  if (l > 0 && nelements > 0) {
    double B, C, dt;
    int    ret, bounces = 0, i=0;
    double this_width, this_height;
    double angle=0;

    if (GVars.fc_freq != 0 || GVars.fc_phase != 0) { /* rotate neutron w/r to guide element */
      /* approximation of rotating straight Fermi Chopper */
      Coords   X = coords_set(x,y,z-l/2);  /* current coordinates of neutron in centered static frame */
      Rotation R;
      double dt=(-z+l/2)/vz; /* time shift to each center of slit package */
      angle=fmod(360*GVars.fc_freq*(t+dt)+GVars.fc_phase, 360); /* in deg */
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

    for (i=0; i<7; GVars.N_reflection[i++] = 0);

    /* propagate to box input (with gravitation) in comp local coords */
    /* A = 0.5 n.g; B = n.v; C = n.(r-W); */
    /* 0=Z0 side: n=(0, 0, -l) ; W = (0, 0, 0) (at z=0, guide input)*/
    B = -l*vz; C = -l*z;

    ret = solve_2nd_order(&dt, NULL, GVars.A[0], B, C);
    if (ret==0) ABSORB;

    if (dt>0.0) PROP_GRAV_DT(dt, GVars.gx, GVars.gy, GVars.gz); else if (angle) ABSORB;
    GVars.N_reflection[6]++;

    this_width  = w1;
    this_height = h1;

  /* check if we are in the box input, else absorb */
    if (fabs(x) > this_width/2 || fabs(y) > this_height/2)
      ABSORB;
    else
    {
      double w_edge, w_adj; /* Channel displacement on X */
      double h_edge, h_adj; /* Channel displacement on Y */
      double w_chnum,h_chnum; /* channel indexes */

      SCATTER;

      /* X: Shift origin to center of channel hit (absorb if hit dividing walls) */
      x += w1/2.0;
      w_chnum = floor(x/(GVars.w1c+d));  /* 0= right side, nslit+1=left side  */
      w_edge  = w_chnum*(GVars.w1c+d);
      if(x - w_edge > GVars.w1c)
      {
        x -= w1/2.0; /* Re-adjust origin */
        ABSORB;
      }
      w_adj = w_edge + (GVars.w1c)/2.0;
      x -= w_adj; w_adj -=  w1/2.0;

      /* Y: Shift origin to center of channel hit (absorb if hit dividing walls) */
      y += h1/2.0;
      h_chnum = floor(y/(GVars.h1c+d));  /* 0= lower side, nslit+1=upper side  */
      h_edge  = h_chnum*(GVars.h1c+d);
      if(y - h_edge > GVars.h1c)
      {
        y -= h1/2.0; /* Re-adjust origin */
        ABSORB;
      }
      h_adj = h_edge + (GVars.h1c)/2.0;
      y -= h_adj; h_adj -=  h1/2.0;

      /* neutron is now in the input window of the guide */
      /* do loops on reflections in the box */
      for(;;)
      {
        /* get intersections for all box sides */
        double q, nx,ny,nz;
        double this_length;
        int side=0;

        bounces++;
        /* now look for intersection with guide sides and exit */
        side = Gravity_guide_Trace(&dt, &GVars, x, y, z,
            vx, vy, vz, w_chnum, nslit, h_chnum, nhslit,
            &nx, &ny, &nz);

        /* only positive dt are valid */
        /* exit reflection loops if no intersection (neutron is after box) */
        if (side == 0 || dt <= 0)
          { if (GVars.warnings < 100)
              fprintf(stderr,"%s: warning: neutron has entered guide, but can not exit !\n", GVars.compcurname);
            GVars.warnings++;
            x += w_adj; y += h_adj; ABSORB; } /* should never occur */

        /* propagate to dt */
        PROP_GRAV_DT(dt, GVars.gx, GVars.gy, GVars.gz);

        /* do reflection on speed for l/r/u/d sides */
        if (side == 5) /* neutron reaches end of guide: end loop and exit comp */
          { GVars.N_reflection[side]++; x += w_adj; y += h_adj; SCATTER; x -= w_adj; y -= h_adj; break; }
        /* else reflection on a guide wall */
        if(GVars.M[side] == 0 || Qc == 0 || R0 == 0)  /* walls are absorbing */
          { x += w_adj; y += h_adj; ABSORB; }
        /* handle chamfers */
        this_width = w1+(w2-w1)*z/l;
        this_height= h1+(h2-h1)*z/l;
        this_length= fmod(z, l/nelements);
        /* absorb on input/output of element parts */
        if (GVars.chamfer_z && (this_length<GVars.chamfer_z || this_length>l/nelements-GVars.chamfer_z))
        { x += w_adj; y += h_adj; ABSORB; }
        /* absorb on l/r/t/b sides */
        if (GVars.chamfer_lr && (side==1 || side==2) && (fabs(y+h_adj)>this_height/2-GVars.chamfer_lr))
        { x += w_adj; y += h_adj; ABSORB; }
        if (GVars.chamfer_tb && (side==3 || side==4) && (fabs(x+w_adj)>this_width/2- GVars.chamfer_tb))
        { x += w_adj; y += h_adj; ABSORB; }
        /* change/mirror velocity: h_f = v - n.2*n.v/|n|^2 */
        GVars.N_reflection[side]++; /* GVars.norm_n2 > 0 was checked at INIT */
        /* compute n.v using current values */
        B = scalar_prod(vx,vy,vz,nx,ny,nz);
        dt = 2*B/GVars.norm_n2[side]; /* 2*n.v/|n|^2 */
        vx -= nx*dt;
        vy -= ny*dt;
        vz -= nz*dt;

        /* compute q and modify neutron weight */
        /* scattering q=|n_i-n_f| = V2Q*|vf - v| = V2Q*2*n.v/|n| */
        q = 2*V2Q*fabs(B)/GVars.norm_n[side];

        if (reflect && strlen(reflect) && strcmp(reflect,"NULL") && strcmp(reflect,"0"))
          TableReflecFunc(q, &pTable, &B);
        else {
          double par[] = {R0, Qc, GVars.Alpha[side], GVars.M[side], W};
          StdReflecFunc(q, par, &B);
        }
        if (B <= 0) { x += w_adj; y += h_adj; ABSORB; }
        else p *= B;
        x += w_adj; y += h_adj; SCATTER; x -= w_adj; y -= h_adj;
        GVars.N_reflection[0]++;
        /* go to the next reflection */
        if (bounces > 1000) ABSORB;
      } /* end for */
      x += w_adj; y += h_adj; /* Re-adjust origin after SCATTER */
    }

    if (GVars.fc_freq != 0 || GVars.fc_phase != 0) { /* rotate back neutron w/r to guide element */
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

  } /* if l */
}
#line 18538 "./ISIS_SANS2d.c"
}   /* End of bender8=Guide_gravity() SETTING parameter declarations. */
#undef pTable
#undef GVars
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompbender8:
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

  /* TRACE Component bender9 [13] */
  mccoordschange(mcposrbender9, mcrotrbender9,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component bender9 (without coords transformations) */
  mcJumpTrace_bender9:
  SIG_MESSAGE("bender9 (Trace)");
  mcDEBUG_COMP("bender9")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbCompbender9
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
#define mccompcurname  bender9
#define mccompcurtype  Guide_gravity
#define mccompcurindex 13
#define GVars mccbender9_GVars
#define pTable mccbender9_pTable
{   /* Declarations of bender9=Guide_gravity() SETTING parameters. */
MCNUM w1 = mccbender9_w1;
MCNUM h1 = mccbender9_h1;
MCNUM w2 = mccbender9_w2;
MCNUM h2 = mccbender9_h2;
MCNUM l = mccbender9_l;
MCNUM R0 = mccbender9_R0;
MCNUM Qc = mccbender9_Qc;
MCNUM alpha = mccbender9_alpha;
MCNUM m = mccbender9_m;
MCNUM W = mccbender9_W;
MCNUM nslit = mccbender9_nslit;
MCNUM d = mccbender9_d;
MCNUM mleft = mccbender9_mleft;
MCNUM mright = mccbender9_mright;
MCNUM mtop = mccbender9_mtop;
MCNUM mbottom = mccbender9_mbottom;
MCNUM nhslit = mccbender9_nhslit;
MCNUM G = mccbender9_G;
MCNUM aleft = mccbender9_aleft;
MCNUM aright = mccbender9_aright;
MCNUM atop = mccbender9_atop;
MCNUM abottom = mccbender9_abottom;
MCNUM wavy = mccbender9_wavy;
MCNUM wavy_z = mccbender9_wavy_z;
MCNUM wavy_tb = mccbender9_wavy_tb;
MCNUM wavy_lr = mccbender9_wavy_lr;
MCNUM chamfers = mccbender9_chamfers;
MCNUM chamfers_z = mccbender9_chamfers_z;
MCNUM chamfers_lr = mccbender9_chamfers_lr;
MCNUM chamfers_tb = mccbender9_chamfers_tb;
MCNUM nelements = mccbender9_nelements;
MCNUM nu = mccbender9_nu;
MCNUM phase = mccbender9_phase;
char* reflect = mccbender9_reflect;
#line 392 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide_gravity.comp"
{
  if (l > 0 && nelements > 0) {
    double B, C, dt;
    int    ret, bounces = 0, i=0;
    double this_width, this_height;
    double angle=0;

    if (GVars.fc_freq != 0 || GVars.fc_phase != 0) { /* rotate neutron w/r to guide element */
      /* approximation of rotating straight Fermi Chopper */
      Coords   X = coords_set(x,y,z-l/2);  /* current coordinates of neutron in centered static frame */
      Rotation R;
      double dt=(-z+l/2)/vz; /* time shift to each center of slit package */
      angle=fmod(360*GVars.fc_freq*(t+dt)+GVars.fc_phase, 360); /* in deg */
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

    for (i=0; i<7; GVars.N_reflection[i++] = 0);

    /* propagate to box input (with gravitation) in comp local coords */
    /* A = 0.5 n.g; B = n.v; C = n.(r-W); */
    /* 0=Z0 side: n=(0, 0, -l) ; W = (0, 0, 0) (at z=0, guide input)*/
    B = -l*vz; C = -l*z;

    ret = solve_2nd_order(&dt, NULL, GVars.A[0], B, C);
    if (ret==0) ABSORB;

    if (dt>0.0) PROP_GRAV_DT(dt, GVars.gx, GVars.gy, GVars.gz); else if (angle) ABSORB;
    GVars.N_reflection[6]++;

    this_width  = w1;
    this_height = h1;

  /* check if we are in the box input, else absorb */
    if (fabs(x) > this_width/2 || fabs(y) > this_height/2)
      ABSORB;
    else
    {
      double w_edge, w_adj; /* Channel displacement on X */
      double h_edge, h_adj; /* Channel displacement on Y */
      double w_chnum,h_chnum; /* channel indexes */

      SCATTER;

      /* X: Shift origin to center of channel hit (absorb if hit dividing walls) */
      x += w1/2.0;
      w_chnum = floor(x/(GVars.w1c+d));  /* 0= right side, nslit+1=left side  */
      w_edge  = w_chnum*(GVars.w1c+d);
      if(x - w_edge > GVars.w1c)
      {
        x -= w1/2.0; /* Re-adjust origin */
        ABSORB;
      }
      w_adj = w_edge + (GVars.w1c)/2.0;
      x -= w_adj; w_adj -=  w1/2.0;

      /* Y: Shift origin to center of channel hit (absorb if hit dividing walls) */
      y += h1/2.0;
      h_chnum = floor(y/(GVars.h1c+d));  /* 0= lower side, nslit+1=upper side  */
      h_edge  = h_chnum*(GVars.h1c+d);
      if(y - h_edge > GVars.h1c)
      {
        y -= h1/2.0; /* Re-adjust origin */
        ABSORB;
      }
      h_adj = h_edge + (GVars.h1c)/2.0;
      y -= h_adj; h_adj -=  h1/2.0;

      /* neutron is now in the input window of the guide */
      /* do loops on reflections in the box */
      for(;;)
      {
        /* get intersections for all box sides */
        double q, nx,ny,nz;
        double this_length;
        int side=0;

        bounces++;
        /* now look for intersection with guide sides and exit */
        side = Gravity_guide_Trace(&dt, &GVars, x, y, z,
            vx, vy, vz, w_chnum, nslit, h_chnum, nhslit,
            &nx, &ny, &nz);

        /* only positive dt are valid */
        /* exit reflection loops if no intersection (neutron is after box) */
        if (side == 0 || dt <= 0)
          { if (GVars.warnings < 100)
              fprintf(stderr,"%s: warning: neutron has entered guide, but can not exit !\n", GVars.compcurname);
            GVars.warnings++;
            x += w_adj; y += h_adj; ABSORB; } /* should never occur */

        /* propagate to dt */
        PROP_GRAV_DT(dt, GVars.gx, GVars.gy, GVars.gz);

        /* do reflection on speed for l/r/u/d sides */
        if (side == 5) /* neutron reaches end of guide: end loop and exit comp */
          { GVars.N_reflection[side]++; x += w_adj; y += h_adj; SCATTER; x -= w_adj; y -= h_adj; break; }
        /* else reflection on a guide wall */
        if(GVars.M[side] == 0 || Qc == 0 || R0 == 0)  /* walls are absorbing */
          { x += w_adj; y += h_adj; ABSORB; }
        /* handle chamfers */
        this_width = w1+(w2-w1)*z/l;
        this_height= h1+(h2-h1)*z/l;
        this_length= fmod(z, l/nelements);
        /* absorb on input/output of element parts */
        if (GVars.chamfer_z && (this_length<GVars.chamfer_z || this_length>l/nelements-GVars.chamfer_z))
        { x += w_adj; y += h_adj; ABSORB; }
        /* absorb on l/r/t/b sides */
        if (GVars.chamfer_lr && (side==1 || side==2) && (fabs(y+h_adj)>this_height/2-GVars.chamfer_lr))
        { x += w_adj; y += h_adj; ABSORB; }
        if (GVars.chamfer_tb && (side==3 || side==4) && (fabs(x+w_adj)>this_width/2- GVars.chamfer_tb))
        { x += w_adj; y += h_adj; ABSORB; }
        /* change/mirror velocity: h_f = v - n.2*n.v/|n|^2 */
        GVars.N_reflection[side]++; /* GVars.norm_n2 > 0 was checked at INIT */
        /* compute n.v using current values */
        B = scalar_prod(vx,vy,vz,nx,ny,nz);
        dt = 2*B/GVars.norm_n2[side]; /* 2*n.v/|n|^2 */
        vx -= nx*dt;
        vy -= ny*dt;
        vz -= nz*dt;

        /* compute q and modify neutron weight */
        /* scattering q=|n_i-n_f| = V2Q*|vf - v| = V2Q*2*n.v/|n| */
        q = 2*V2Q*fabs(B)/GVars.norm_n[side];

        if (reflect && strlen(reflect) && strcmp(reflect,"NULL") && strcmp(reflect,"0"))
          TableReflecFunc(q, &pTable, &B);
        else {
          double par[] = {R0, Qc, GVars.Alpha[side], GVars.M[side], W};
          StdReflecFunc(q, par, &B);
        }
        if (B <= 0) { x += w_adj; y += h_adj; ABSORB; }
        else p *= B;
        x += w_adj; y += h_adj; SCATTER; x -= w_adj; y -= h_adj;
        GVars.N_reflection[0]++;
        /* go to the next reflection */
        if (bounces > 1000) ABSORB;
      } /* end for */
      x += w_adj; y += h_adj; /* Re-adjust origin after SCATTER */
    }

    if (GVars.fc_freq != 0 || GVars.fc_phase != 0) { /* rotate back neutron w/r to guide element */
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

  } /* if l */
}
#line 18851 "./ISIS_SANS2d.c"
}   /* End of bender9=Guide_gravity() SETTING parameter declarations. */
#undef pTable
#undef GVars
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompbender9:
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

  /* TRACE Component bender10 [14] */
  mccoordschange(mcposrbender10, mcrotrbender10,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component bender10 (without coords transformations) */
  mcJumpTrace_bender10:
  SIG_MESSAGE("bender10 (Trace)");
  mcDEBUG_COMP("bender10")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbCompbender10
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
#define mccompcurname  bender10
#define mccompcurtype  Guide_gravity
#define mccompcurindex 14
#define GVars mccbender10_GVars
#define pTable mccbender10_pTable
{   /* Declarations of bender10=Guide_gravity() SETTING parameters. */
MCNUM w1 = mccbender10_w1;
MCNUM h1 = mccbender10_h1;
MCNUM w2 = mccbender10_w2;
MCNUM h2 = mccbender10_h2;
MCNUM l = mccbender10_l;
MCNUM R0 = mccbender10_R0;
MCNUM Qc = mccbender10_Qc;
MCNUM alpha = mccbender10_alpha;
MCNUM m = mccbender10_m;
MCNUM W = mccbender10_W;
MCNUM nslit = mccbender10_nslit;
MCNUM d = mccbender10_d;
MCNUM mleft = mccbender10_mleft;
MCNUM mright = mccbender10_mright;
MCNUM mtop = mccbender10_mtop;
MCNUM mbottom = mccbender10_mbottom;
MCNUM nhslit = mccbender10_nhslit;
MCNUM G = mccbender10_G;
MCNUM aleft = mccbender10_aleft;
MCNUM aright = mccbender10_aright;
MCNUM atop = mccbender10_atop;
MCNUM abottom = mccbender10_abottom;
MCNUM wavy = mccbender10_wavy;
MCNUM wavy_z = mccbender10_wavy_z;
MCNUM wavy_tb = mccbender10_wavy_tb;
MCNUM wavy_lr = mccbender10_wavy_lr;
MCNUM chamfers = mccbender10_chamfers;
MCNUM chamfers_z = mccbender10_chamfers_z;
MCNUM chamfers_lr = mccbender10_chamfers_lr;
MCNUM chamfers_tb = mccbender10_chamfers_tb;
MCNUM nelements = mccbender10_nelements;
MCNUM nu = mccbender10_nu;
MCNUM phase = mccbender10_phase;
char* reflect = mccbender10_reflect;
#line 392 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide_gravity.comp"
{
  if (l > 0 && nelements > 0) {
    double B, C, dt;
    int    ret, bounces = 0, i=0;
    double this_width, this_height;
    double angle=0;

    if (GVars.fc_freq != 0 || GVars.fc_phase != 0) { /* rotate neutron w/r to guide element */
      /* approximation of rotating straight Fermi Chopper */
      Coords   X = coords_set(x,y,z-l/2);  /* current coordinates of neutron in centered static frame */
      Rotation R;
      double dt=(-z+l/2)/vz; /* time shift to each center of slit package */
      angle=fmod(360*GVars.fc_freq*(t+dt)+GVars.fc_phase, 360); /* in deg */
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

    for (i=0; i<7; GVars.N_reflection[i++] = 0);

    /* propagate to box input (with gravitation) in comp local coords */
    /* A = 0.5 n.g; B = n.v; C = n.(r-W); */
    /* 0=Z0 side: n=(0, 0, -l) ; W = (0, 0, 0) (at z=0, guide input)*/
    B = -l*vz; C = -l*z;

    ret = solve_2nd_order(&dt, NULL, GVars.A[0], B, C);
    if (ret==0) ABSORB;

    if (dt>0.0) PROP_GRAV_DT(dt, GVars.gx, GVars.gy, GVars.gz); else if (angle) ABSORB;
    GVars.N_reflection[6]++;

    this_width  = w1;
    this_height = h1;

  /* check if we are in the box input, else absorb */
    if (fabs(x) > this_width/2 || fabs(y) > this_height/2)
      ABSORB;
    else
    {
      double w_edge, w_adj; /* Channel displacement on X */
      double h_edge, h_adj; /* Channel displacement on Y */
      double w_chnum,h_chnum; /* channel indexes */

      SCATTER;

      /* X: Shift origin to center of channel hit (absorb if hit dividing walls) */
      x += w1/2.0;
      w_chnum = floor(x/(GVars.w1c+d));  /* 0= right side, nslit+1=left side  */
      w_edge  = w_chnum*(GVars.w1c+d);
      if(x - w_edge > GVars.w1c)
      {
        x -= w1/2.0; /* Re-adjust origin */
        ABSORB;
      }
      w_adj = w_edge + (GVars.w1c)/2.0;
      x -= w_adj; w_adj -=  w1/2.0;

      /* Y: Shift origin to center of channel hit (absorb if hit dividing walls) */
      y += h1/2.0;
      h_chnum = floor(y/(GVars.h1c+d));  /* 0= lower side, nslit+1=upper side  */
      h_edge  = h_chnum*(GVars.h1c+d);
      if(y - h_edge > GVars.h1c)
      {
        y -= h1/2.0; /* Re-adjust origin */
        ABSORB;
      }
      h_adj = h_edge + (GVars.h1c)/2.0;
      y -= h_adj; h_adj -=  h1/2.0;

      /* neutron is now in the input window of the guide */
      /* do loops on reflections in the box */
      for(;;)
      {
        /* get intersections for all box sides */
        double q, nx,ny,nz;
        double this_length;
        int side=0;

        bounces++;
        /* now look for intersection with guide sides and exit */
        side = Gravity_guide_Trace(&dt, &GVars, x, y, z,
            vx, vy, vz, w_chnum, nslit, h_chnum, nhslit,
            &nx, &ny, &nz);

        /* only positive dt are valid */
        /* exit reflection loops if no intersection (neutron is after box) */
        if (side == 0 || dt <= 0)
          { if (GVars.warnings < 100)
              fprintf(stderr,"%s: warning: neutron has entered guide, but can not exit !\n", GVars.compcurname);
            GVars.warnings++;
            x += w_adj; y += h_adj; ABSORB; } /* should never occur */

        /* propagate to dt */
        PROP_GRAV_DT(dt, GVars.gx, GVars.gy, GVars.gz);

        /* do reflection on speed for l/r/u/d sides */
        if (side == 5) /* neutron reaches end of guide: end loop and exit comp */
          { GVars.N_reflection[side]++; x += w_adj; y += h_adj; SCATTER; x -= w_adj; y -= h_adj; break; }
        /* else reflection on a guide wall */
        if(GVars.M[side] == 0 || Qc == 0 || R0 == 0)  /* walls are absorbing */
          { x += w_adj; y += h_adj; ABSORB; }
        /* handle chamfers */
        this_width = w1+(w2-w1)*z/l;
        this_height= h1+(h2-h1)*z/l;
        this_length= fmod(z, l/nelements);
        /* absorb on input/output of element parts */
        if (GVars.chamfer_z && (this_length<GVars.chamfer_z || this_length>l/nelements-GVars.chamfer_z))
        { x += w_adj; y += h_adj; ABSORB; }
        /* absorb on l/r/t/b sides */
        if (GVars.chamfer_lr && (side==1 || side==2) && (fabs(y+h_adj)>this_height/2-GVars.chamfer_lr))
        { x += w_adj; y += h_adj; ABSORB; }
        if (GVars.chamfer_tb && (side==3 || side==4) && (fabs(x+w_adj)>this_width/2- GVars.chamfer_tb))
        { x += w_adj; y += h_adj; ABSORB; }
        /* change/mirror velocity: h_f = v - n.2*n.v/|n|^2 */
        GVars.N_reflection[side]++; /* GVars.norm_n2 > 0 was checked at INIT */
        /* compute n.v using current values */
        B = scalar_prod(vx,vy,vz,nx,ny,nz);
        dt = 2*B/GVars.norm_n2[side]; /* 2*n.v/|n|^2 */
        vx -= nx*dt;
        vy -= ny*dt;
        vz -= nz*dt;

        /* compute q and modify neutron weight */
        /* scattering q=|n_i-n_f| = V2Q*|vf - v| = V2Q*2*n.v/|n| */
        q = 2*V2Q*fabs(B)/GVars.norm_n[side];

        if (reflect && strlen(reflect) && strcmp(reflect,"NULL") && strcmp(reflect,"0"))
          TableReflecFunc(q, &pTable, &B);
        else {
          double par[] = {R0, Qc, GVars.Alpha[side], GVars.M[side], W};
          StdReflecFunc(q, par, &B);
        }
        if (B <= 0) { x += w_adj; y += h_adj; ABSORB; }
        else p *= B;
        x += w_adj; y += h_adj; SCATTER; x -= w_adj; y -= h_adj;
        GVars.N_reflection[0]++;
        /* go to the next reflection */
        if (bounces > 1000) ABSORB;
      } /* end for */
      x += w_adj; y += h_adj; /* Re-adjust origin after SCATTER */
    }

    if (GVars.fc_freq != 0 || GVars.fc_phase != 0) { /* rotate back neutron w/r to guide element */
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

  } /* if l */
}
#line 19164 "./ISIS_SANS2d.c"
}   /* End of bender10=Guide_gravity() SETTING parameter declarations. */
#undef pTable
#undef GVars
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompbender10:
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

  /* TRACE Component lmonb [15] */
  mccoordschange(mcposrlmonb, mcrotrlmonb,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component lmonb (without coords transformations) */
  mcJumpTrace_lmonb:
  SIG_MESSAGE("lmonb (Trace)");
  mcDEBUG_COMP("lmonb")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbComplmonb
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
#define mccompcurname  lmonb
#define mccompcurtype  L_monitor
#define mccompcurindex 15
#define nL mcclmonb_nL
#define L_N mcclmonb_L_N
#define L_p mcclmonb_L_p
#define L_p2 mcclmonb_L_p2
{   /* Declarations of lmonb=L_monitor() SETTING parameters. */
char* filename = mcclmonb_filename;
MCNUM xmin = mcclmonb_xmin;
MCNUM xmax = mcclmonb_xmax;
MCNUM ymin = mcclmonb_ymin;
MCNUM ymax = mcclmonb_ymax;
MCNUM xwidth = mcclmonb_xwidth;
MCNUM yheight = mcclmonb_yheight;
MCNUM Lmin = mcclmonb_Lmin;
MCNUM Lmax = mcclmonb_Lmax;
MCNUM restore_neutron = mcclmonb_restore_neutron;
int nowritefile = mcclmonb_nowritefile;
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
#line 19309 "./ISIS_SANS2d.c"
}   /* End of lmonb=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbComplmonb:
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

  /* TRACE Component psd2 [16] */
  mccoordschange(mcposrpsd2, mcrotrpsd2,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component psd2 (without coords transformations) */
  mcJumpTrace_psd2:
  SIG_MESSAGE("psd2 (Trace)");
  mcDEBUG_COMP("psd2")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbComppsd2
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
#define mccompcurname  psd2
#define mccompcurtype  PSD_monitor
#define mccompcurindex 16
#define PSD_N mccpsd2_PSD_N
#define PSD_p mccpsd2_PSD_p
#define PSD_p2 mccpsd2_PSD_p2
{   /* Declarations of psd2=PSD_monitor() SETTING parameters. */
int nx = mccpsd2_nx;
int ny = mccpsd2_ny;
char* filename = mccpsd2_filename;
MCNUM xmin = mccpsd2_xmin;
MCNUM xmax = mccpsd2_xmax;
MCNUM ymin = mccpsd2_ymin;
MCNUM ymax = mccpsd2_ymax;
MCNUM xwidth = mccpsd2_xwidth;
MCNUM yheight = mccpsd2_yheight;
MCNUM restore_neutron = mccpsd2_restore_neutron;
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
#line 19447 "./ISIS_SANS2d.c"
}   /* End of psd2=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbComppsd2:
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

  /* TRACE Component guide_in [17] */
  mccoordschange(mcposrguide_in, mcrotrguide_in,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component guide_in (without coords transformations) */
  mcJumpTrace_guide_in:
  SIG_MESSAGE("guide_in (Trace)");
  mcDEBUG_COMP("guide_in")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbCompguide_in
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
#define mccompcurname  guide_in
#define mccompcurtype  Slit
#define mccompcurindex 17
{   /* Declarations of guide_in=Slit() SETTING parameters. */
MCNUM xmin = mccguide_in_xmin;
MCNUM xmax = mccguide_in_xmax;
MCNUM ymin = mccguide_in_ymin;
MCNUM ymax = mccguide_in_ymax;
MCNUM radius = mccguide_in_radius;
MCNUM xwidth = mccguide_in_xwidth;
MCNUM yheight = mccguide_in_yheight;
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
#line 19574 "./ISIS_SANS2d.c"
}   /* End of guide_in=Slit() SETTING parameter declarations. */
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompguide_in:
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

  /* TRACE Component guide_straight1 [18] */
  mccoordschange(mcposrguide_straight1, mcrotrguide_straight1,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component guide_straight1 (without coords transformations) */
  mcJumpTrace_guide_straight1:
  SIG_MESSAGE("guide_straight1 (Trace)");
  mcDEBUG_COMP("guide_straight1")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbCompguide_straight1
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
#define mccompcurname  guide_straight1
#define mccompcurtype  Guide_gravity
#define mccompcurindex 18
#define GVars mccguide_straight1_GVars
#define pTable mccguide_straight1_pTable
{   /* Declarations of guide_straight1=Guide_gravity() SETTING parameters. */
MCNUM w1 = mccguide_straight1_w1;
MCNUM h1 = mccguide_straight1_h1;
MCNUM w2 = mccguide_straight1_w2;
MCNUM h2 = mccguide_straight1_h2;
MCNUM l = mccguide_straight1_l;
MCNUM R0 = mccguide_straight1_R0;
MCNUM Qc = mccguide_straight1_Qc;
MCNUM alpha = mccguide_straight1_alpha;
MCNUM m = mccguide_straight1_m;
MCNUM W = mccguide_straight1_W;
MCNUM nslit = mccguide_straight1_nslit;
MCNUM d = mccguide_straight1_d;
MCNUM mleft = mccguide_straight1_mleft;
MCNUM mright = mccguide_straight1_mright;
MCNUM mtop = mccguide_straight1_mtop;
MCNUM mbottom = mccguide_straight1_mbottom;
MCNUM nhslit = mccguide_straight1_nhslit;
MCNUM G = mccguide_straight1_G;
MCNUM aleft = mccguide_straight1_aleft;
MCNUM aright = mccguide_straight1_aright;
MCNUM atop = mccguide_straight1_atop;
MCNUM abottom = mccguide_straight1_abottom;
MCNUM wavy = mccguide_straight1_wavy;
MCNUM wavy_z = mccguide_straight1_wavy_z;
MCNUM wavy_tb = mccguide_straight1_wavy_tb;
MCNUM wavy_lr = mccguide_straight1_wavy_lr;
MCNUM chamfers = mccguide_straight1_chamfers;
MCNUM chamfers_z = mccguide_straight1_chamfers_z;
MCNUM chamfers_lr = mccguide_straight1_chamfers_lr;
MCNUM chamfers_tb = mccguide_straight1_chamfers_tb;
MCNUM nelements = mccguide_straight1_nelements;
MCNUM nu = mccguide_straight1_nu;
MCNUM phase = mccguide_straight1_phase;
char* reflect = mccguide_straight1_reflect;
#line 392 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide_gravity.comp"
{
  if (l > 0 && nelements > 0) {
    double B, C, dt;
    int    ret, bounces = 0, i=0;
    double this_width, this_height;
    double angle=0;

    if (GVars.fc_freq != 0 || GVars.fc_phase != 0) { /* rotate neutron w/r to guide element */
      /* approximation of rotating straight Fermi Chopper */
      Coords   X = coords_set(x,y,z-l/2);  /* current coordinates of neutron in centered static frame */
      Rotation R;
      double dt=(-z+l/2)/vz; /* time shift to each center of slit package */
      angle=fmod(360*GVars.fc_freq*(t+dt)+GVars.fc_phase, 360); /* in deg */
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

    for (i=0; i<7; GVars.N_reflection[i++] = 0);

    /* propagate to box input (with gravitation) in comp local coords */
    /* A = 0.5 n.g; B = n.v; C = n.(r-W); */
    /* 0=Z0 side: n=(0, 0, -l) ; W = (0, 0, 0) (at z=0, guide input)*/
    B = -l*vz; C = -l*z;

    ret = solve_2nd_order(&dt, NULL, GVars.A[0], B, C);
    if (ret==0) ABSORB;

    if (dt>0.0) PROP_GRAV_DT(dt, GVars.gx, GVars.gy, GVars.gz); else if (angle) ABSORB;
    GVars.N_reflection[6]++;

    this_width  = w1;
    this_height = h1;

  /* check if we are in the box input, else absorb */
    if (fabs(x) > this_width/2 || fabs(y) > this_height/2)
      ABSORB;
    else
    {
      double w_edge, w_adj; /* Channel displacement on X */
      double h_edge, h_adj; /* Channel displacement on Y */
      double w_chnum,h_chnum; /* channel indexes */

      SCATTER;

      /* X: Shift origin to center of channel hit (absorb if hit dividing walls) */
      x += w1/2.0;
      w_chnum = floor(x/(GVars.w1c+d));  /* 0= right side, nslit+1=left side  */
      w_edge  = w_chnum*(GVars.w1c+d);
      if(x - w_edge > GVars.w1c)
      {
        x -= w1/2.0; /* Re-adjust origin */
        ABSORB;
      }
      w_adj = w_edge + (GVars.w1c)/2.0;
      x -= w_adj; w_adj -=  w1/2.0;

      /* Y: Shift origin to center of channel hit (absorb if hit dividing walls) */
      y += h1/2.0;
      h_chnum = floor(y/(GVars.h1c+d));  /* 0= lower side, nslit+1=upper side  */
      h_edge  = h_chnum*(GVars.h1c+d);
      if(y - h_edge > GVars.h1c)
      {
        y -= h1/2.0; /* Re-adjust origin */
        ABSORB;
      }
      h_adj = h_edge + (GVars.h1c)/2.0;
      y -= h_adj; h_adj -=  h1/2.0;

      /* neutron is now in the input window of the guide */
      /* do loops on reflections in the box */
      for(;;)
      {
        /* get intersections for all box sides */
        double q, nx,ny,nz;
        double this_length;
        int side=0;

        bounces++;
        /* now look for intersection with guide sides and exit */
        side = Gravity_guide_Trace(&dt, &GVars, x, y, z,
            vx, vy, vz, w_chnum, nslit, h_chnum, nhslit,
            &nx, &ny, &nz);

        /* only positive dt are valid */
        /* exit reflection loops if no intersection (neutron is after box) */
        if (side == 0 || dt <= 0)
          { if (GVars.warnings < 100)
              fprintf(stderr,"%s: warning: neutron has entered guide, but can not exit !\n", GVars.compcurname);
            GVars.warnings++;
            x += w_adj; y += h_adj; ABSORB; } /* should never occur */

        /* propagate to dt */
        PROP_GRAV_DT(dt, GVars.gx, GVars.gy, GVars.gz);

        /* do reflection on speed for l/r/u/d sides */
        if (side == 5) /* neutron reaches end of guide: end loop and exit comp */
          { GVars.N_reflection[side]++; x += w_adj; y += h_adj; SCATTER; x -= w_adj; y -= h_adj; break; }
        /* else reflection on a guide wall */
        if(GVars.M[side] == 0 || Qc == 0 || R0 == 0)  /* walls are absorbing */
          { x += w_adj; y += h_adj; ABSORB; }
        /* handle chamfers */
        this_width = w1+(w2-w1)*z/l;
        this_height= h1+(h2-h1)*z/l;
        this_length= fmod(z, l/nelements);
        /* absorb on input/output of element parts */
        if (GVars.chamfer_z && (this_length<GVars.chamfer_z || this_length>l/nelements-GVars.chamfer_z))
        { x += w_adj; y += h_adj; ABSORB; }
        /* absorb on l/r/t/b sides */
        if (GVars.chamfer_lr && (side==1 || side==2) && (fabs(y+h_adj)>this_height/2-GVars.chamfer_lr))
        { x += w_adj; y += h_adj; ABSORB; }
        if (GVars.chamfer_tb && (side==3 || side==4) && (fabs(x+w_adj)>this_width/2- GVars.chamfer_tb))
        { x += w_adj; y += h_adj; ABSORB; }
        /* change/mirror velocity: h_f = v - n.2*n.v/|n|^2 */
        GVars.N_reflection[side]++; /* GVars.norm_n2 > 0 was checked at INIT */
        /* compute n.v using current values */
        B = scalar_prod(vx,vy,vz,nx,ny,nz);
        dt = 2*B/GVars.norm_n2[side]; /* 2*n.v/|n|^2 */
        vx -= nx*dt;
        vy -= ny*dt;
        vz -= nz*dt;

        /* compute q and modify neutron weight */
        /* scattering q=|n_i-n_f| = V2Q*|vf - v| = V2Q*2*n.v/|n| */
        q = 2*V2Q*fabs(B)/GVars.norm_n[side];

        if (reflect && strlen(reflect) && strcmp(reflect,"NULL") && strcmp(reflect,"0"))
          TableReflecFunc(q, &pTable, &B);
        else {
          double par[] = {R0, Qc, GVars.Alpha[side], GVars.M[side], W};
          StdReflecFunc(q, par, &B);
        }
        if (B <= 0) { x += w_adj; y += h_adj; ABSORB; }
        else p *= B;
        x += w_adj; y += h_adj; SCATTER; x -= w_adj; y -= h_adj;
        GVars.N_reflection[0]++;
        /* go to the next reflection */
        if (bounces > 1000) ABSORB;
      } /* end for */
      x += w_adj; y += h_adj; /* Re-adjust origin after SCATTER */
    }

    if (GVars.fc_freq != 0 || GVars.fc_phase != 0) { /* rotate back neutron w/r to guide element */
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

  } /* if l */
}
#line 19885 "./ISIS_SANS2d.c"
}   /* End of guide_straight1=Guide_gravity() SETTING parameter declarations. */
#undef pTable
#undef GVars
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompguide_straight1:
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

  /* TRACE Component guide_straight2 [19] */
  mccoordschange(mcposrguide_straight2, mcrotrguide_straight2,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component guide_straight2 (without coords transformations) */
  mcJumpTrace_guide_straight2:
  SIG_MESSAGE("guide_straight2 (Trace)");
  mcDEBUG_COMP("guide_straight2")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbCompguide_straight2
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
#define mccompcurname  guide_straight2
#define mccompcurtype  Guide_gravity
#define mccompcurindex 19
#define GVars mccguide_straight2_GVars
#define pTable mccguide_straight2_pTable
{   /* Declarations of guide_straight2=Guide_gravity() SETTING parameters. */
MCNUM w1 = mccguide_straight2_w1;
MCNUM h1 = mccguide_straight2_h1;
MCNUM w2 = mccguide_straight2_w2;
MCNUM h2 = mccguide_straight2_h2;
MCNUM l = mccguide_straight2_l;
MCNUM R0 = mccguide_straight2_R0;
MCNUM Qc = mccguide_straight2_Qc;
MCNUM alpha = mccguide_straight2_alpha;
MCNUM m = mccguide_straight2_m;
MCNUM W = mccguide_straight2_W;
MCNUM nslit = mccguide_straight2_nslit;
MCNUM d = mccguide_straight2_d;
MCNUM mleft = mccguide_straight2_mleft;
MCNUM mright = mccguide_straight2_mright;
MCNUM mtop = mccguide_straight2_mtop;
MCNUM mbottom = mccguide_straight2_mbottom;
MCNUM nhslit = mccguide_straight2_nhslit;
MCNUM G = mccguide_straight2_G;
MCNUM aleft = mccguide_straight2_aleft;
MCNUM aright = mccguide_straight2_aright;
MCNUM atop = mccguide_straight2_atop;
MCNUM abottom = mccguide_straight2_abottom;
MCNUM wavy = mccguide_straight2_wavy;
MCNUM wavy_z = mccguide_straight2_wavy_z;
MCNUM wavy_tb = mccguide_straight2_wavy_tb;
MCNUM wavy_lr = mccguide_straight2_wavy_lr;
MCNUM chamfers = mccguide_straight2_chamfers;
MCNUM chamfers_z = mccguide_straight2_chamfers_z;
MCNUM chamfers_lr = mccguide_straight2_chamfers_lr;
MCNUM chamfers_tb = mccguide_straight2_chamfers_tb;
MCNUM nelements = mccguide_straight2_nelements;
MCNUM nu = mccguide_straight2_nu;
MCNUM phase = mccguide_straight2_phase;
char* reflect = mccguide_straight2_reflect;
#line 392 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide_gravity.comp"
{
  if (l > 0 && nelements > 0) {
    double B, C, dt;
    int    ret, bounces = 0, i=0;
    double this_width, this_height;
    double angle=0;

    if (GVars.fc_freq != 0 || GVars.fc_phase != 0) { /* rotate neutron w/r to guide element */
      /* approximation of rotating straight Fermi Chopper */
      Coords   X = coords_set(x,y,z-l/2);  /* current coordinates of neutron in centered static frame */
      Rotation R;
      double dt=(-z+l/2)/vz; /* time shift to each center of slit package */
      angle=fmod(360*GVars.fc_freq*(t+dt)+GVars.fc_phase, 360); /* in deg */
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

    for (i=0; i<7; GVars.N_reflection[i++] = 0);

    /* propagate to box input (with gravitation) in comp local coords */
    /* A = 0.5 n.g; B = n.v; C = n.(r-W); */
    /* 0=Z0 side: n=(0, 0, -l) ; W = (0, 0, 0) (at z=0, guide input)*/
    B = -l*vz; C = -l*z;

    ret = solve_2nd_order(&dt, NULL, GVars.A[0], B, C);
    if (ret==0) ABSORB;

    if (dt>0.0) PROP_GRAV_DT(dt, GVars.gx, GVars.gy, GVars.gz); else if (angle) ABSORB;
    GVars.N_reflection[6]++;

    this_width  = w1;
    this_height = h1;

  /* check if we are in the box input, else absorb */
    if (fabs(x) > this_width/2 || fabs(y) > this_height/2)
      ABSORB;
    else
    {
      double w_edge, w_adj; /* Channel displacement on X */
      double h_edge, h_adj; /* Channel displacement on Y */
      double w_chnum,h_chnum; /* channel indexes */

      SCATTER;

      /* X: Shift origin to center of channel hit (absorb if hit dividing walls) */
      x += w1/2.0;
      w_chnum = floor(x/(GVars.w1c+d));  /* 0= right side, nslit+1=left side  */
      w_edge  = w_chnum*(GVars.w1c+d);
      if(x - w_edge > GVars.w1c)
      {
        x -= w1/2.0; /* Re-adjust origin */
        ABSORB;
      }
      w_adj = w_edge + (GVars.w1c)/2.0;
      x -= w_adj; w_adj -=  w1/2.0;

      /* Y: Shift origin to center of channel hit (absorb if hit dividing walls) */
      y += h1/2.0;
      h_chnum = floor(y/(GVars.h1c+d));  /* 0= lower side, nslit+1=upper side  */
      h_edge  = h_chnum*(GVars.h1c+d);
      if(y - h_edge > GVars.h1c)
      {
        y -= h1/2.0; /* Re-adjust origin */
        ABSORB;
      }
      h_adj = h_edge + (GVars.h1c)/2.0;
      y -= h_adj; h_adj -=  h1/2.0;

      /* neutron is now in the input window of the guide */
      /* do loops on reflections in the box */
      for(;;)
      {
        /* get intersections for all box sides */
        double q, nx,ny,nz;
        double this_length;
        int side=0;

        bounces++;
        /* now look for intersection with guide sides and exit */
        side = Gravity_guide_Trace(&dt, &GVars, x, y, z,
            vx, vy, vz, w_chnum, nslit, h_chnum, nhslit,
            &nx, &ny, &nz);

        /* only positive dt are valid */
        /* exit reflection loops if no intersection (neutron is after box) */
        if (side == 0 || dt <= 0)
          { if (GVars.warnings < 100)
              fprintf(stderr,"%s: warning: neutron has entered guide, but can not exit !\n", GVars.compcurname);
            GVars.warnings++;
            x += w_adj; y += h_adj; ABSORB; } /* should never occur */

        /* propagate to dt */
        PROP_GRAV_DT(dt, GVars.gx, GVars.gy, GVars.gz);

        /* do reflection on speed for l/r/u/d sides */
        if (side == 5) /* neutron reaches end of guide: end loop and exit comp */
          { GVars.N_reflection[side]++; x += w_adj; y += h_adj; SCATTER; x -= w_adj; y -= h_adj; break; }
        /* else reflection on a guide wall */
        if(GVars.M[side] == 0 || Qc == 0 || R0 == 0)  /* walls are absorbing */
          { x += w_adj; y += h_adj; ABSORB; }
        /* handle chamfers */
        this_width = w1+(w2-w1)*z/l;
        this_height= h1+(h2-h1)*z/l;
        this_length= fmod(z, l/nelements);
        /* absorb on input/output of element parts */
        if (GVars.chamfer_z && (this_length<GVars.chamfer_z || this_length>l/nelements-GVars.chamfer_z))
        { x += w_adj; y += h_adj; ABSORB; }
        /* absorb on l/r/t/b sides */
        if (GVars.chamfer_lr && (side==1 || side==2) && (fabs(y+h_adj)>this_height/2-GVars.chamfer_lr))
        { x += w_adj; y += h_adj; ABSORB; }
        if (GVars.chamfer_tb && (side==3 || side==4) && (fabs(x+w_adj)>this_width/2- GVars.chamfer_tb))
        { x += w_adj; y += h_adj; ABSORB; }
        /* change/mirror velocity: h_f = v - n.2*n.v/|n|^2 */
        GVars.N_reflection[side]++; /* GVars.norm_n2 > 0 was checked at INIT */
        /* compute n.v using current values */
        B = scalar_prod(vx,vy,vz,nx,ny,nz);
        dt = 2*B/GVars.norm_n2[side]; /* 2*n.v/|n|^2 */
        vx -= nx*dt;
        vy -= ny*dt;
        vz -= nz*dt;

        /* compute q and modify neutron weight */
        /* scattering q=|n_i-n_f| = V2Q*|vf - v| = V2Q*2*n.v/|n| */
        q = 2*V2Q*fabs(B)/GVars.norm_n[side];

        if (reflect && strlen(reflect) && strcmp(reflect,"NULL") && strcmp(reflect,"0"))
          TableReflecFunc(q, &pTable, &B);
        else {
          double par[] = {R0, Qc, GVars.Alpha[side], GVars.M[side], W};
          StdReflecFunc(q, par, &B);
        }
        if (B <= 0) { x += w_adj; y += h_adj; ABSORB; }
        else p *= B;
        x += w_adj; y += h_adj; SCATTER; x -= w_adj; y -= h_adj;
        GVars.N_reflection[0]++;
        /* go to the next reflection */
        if (bounces > 1000) ABSORB;
      } /* end for */
      x += w_adj; y += h_adj; /* Re-adjust origin after SCATTER */
    }

    if (GVars.fc_freq != 0 || GVars.fc_phase != 0) { /* rotate back neutron w/r to guide element */
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

  } /* if l */
}
#line 20198 "./ISIS_SANS2d.c"
}   /* End of guide_straight2=Guide_gravity() SETTING parameter declarations. */
#undef pTable
#undef GVars
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompguide_straight2:
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

  /* TRACE Component guide_straight3 [20] */
  mccoordschange(mcposrguide_straight3, mcrotrguide_straight3,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component guide_straight3 (without coords transformations) */
  mcJumpTrace_guide_straight3:
  SIG_MESSAGE("guide_straight3 (Trace)");
  mcDEBUG_COMP("guide_straight3")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbCompguide_straight3
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
#define mccompcurname  guide_straight3
#define mccompcurtype  Guide_gravity
#define mccompcurindex 20
#define GVars mccguide_straight3_GVars
#define pTable mccguide_straight3_pTable
{   /* Declarations of guide_straight3=Guide_gravity() SETTING parameters. */
MCNUM w1 = mccguide_straight3_w1;
MCNUM h1 = mccguide_straight3_h1;
MCNUM w2 = mccguide_straight3_w2;
MCNUM h2 = mccguide_straight3_h2;
MCNUM l = mccguide_straight3_l;
MCNUM R0 = mccguide_straight3_R0;
MCNUM Qc = mccguide_straight3_Qc;
MCNUM alpha = mccguide_straight3_alpha;
MCNUM m = mccguide_straight3_m;
MCNUM W = mccguide_straight3_W;
MCNUM nslit = mccguide_straight3_nslit;
MCNUM d = mccguide_straight3_d;
MCNUM mleft = mccguide_straight3_mleft;
MCNUM mright = mccguide_straight3_mright;
MCNUM mtop = mccguide_straight3_mtop;
MCNUM mbottom = mccguide_straight3_mbottom;
MCNUM nhslit = mccguide_straight3_nhslit;
MCNUM G = mccguide_straight3_G;
MCNUM aleft = mccguide_straight3_aleft;
MCNUM aright = mccguide_straight3_aright;
MCNUM atop = mccguide_straight3_atop;
MCNUM abottom = mccguide_straight3_abottom;
MCNUM wavy = mccguide_straight3_wavy;
MCNUM wavy_z = mccguide_straight3_wavy_z;
MCNUM wavy_tb = mccguide_straight3_wavy_tb;
MCNUM wavy_lr = mccguide_straight3_wavy_lr;
MCNUM chamfers = mccguide_straight3_chamfers;
MCNUM chamfers_z = mccguide_straight3_chamfers_z;
MCNUM chamfers_lr = mccguide_straight3_chamfers_lr;
MCNUM chamfers_tb = mccguide_straight3_chamfers_tb;
MCNUM nelements = mccguide_straight3_nelements;
MCNUM nu = mccguide_straight3_nu;
MCNUM phase = mccguide_straight3_phase;
char* reflect = mccguide_straight3_reflect;
#line 392 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide_gravity.comp"
{
  if (l > 0 && nelements > 0) {
    double B, C, dt;
    int    ret, bounces = 0, i=0;
    double this_width, this_height;
    double angle=0;

    if (GVars.fc_freq != 0 || GVars.fc_phase != 0) { /* rotate neutron w/r to guide element */
      /* approximation of rotating straight Fermi Chopper */
      Coords   X = coords_set(x,y,z-l/2);  /* current coordinates of neutron in centered static frame */
      Rotation R;
      double dt=(-z+l/2)/vz; /* time shift to each center of slit package */
      angle=fmod(360*GVars.fc_freq*(t+dt)+GVars.fc_phase, 360); /* in deg */
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

    for (i=0; i<7; GVars.N_reflection[i++] = 0);

    /* propagate to box input (with gravitation) in comp local coords */
    /* A = 0.5 n.g; B = n.v; C = n.(r-W); */
    /* 0=Z0 side: n=(0, 0, -l) ; W = (0, 0, 0) (at z=0, guide input)*/
    B = -l*vz; C = -l*z;

    ret = solve_2nd_order(&dt, NULL, GVars.A[0], B, C);
    if (ret==0) ABSORB;

    if (dt>0.0) PROP_GRAV_DT(dt, GVars.gx, GVars.gy, GVars.gz); else if (angle) ABSORB;
    GVars.N_reflection[6]++;

    this_width  = w1;
    this_height = h1;

  /* check if we are in the box input, else absorb */
    if (fabs(x) > this_width/2 || fabs(y) > this_height/2)
      ABSORB;
    else
    {
      double w_edge, w_adj; /* Channel displacement on X */
      double h_edge, h_adj; /* Channel displacement on Y */
      double w_chnum,h_chnum; /* channel indexes */

      SCATTER;

      /* X: Shift origin to center of channel hit (absorb if hit dividing walls) */
      x += w1/2.0;
      w_chnum = floor(x/(GVars.w1c+d));  /* 0= right side, nslit+1=left side  */
      w_edge  = w_chnum*(GVars.w1c+d);
      if(x - w_edge > GVars.w1c)
      {
        x -= w1/2.0; /* Re-adjust origin */
        ABSORB;
      }
      w_adj = w_edge + (GVars.w1c)/2.0;
      x -= w_adj; w_adj -=  w1/2.0;

      /* Y: Shift origin to center of channel hit (absorb if hit dividing walls) */
      y += h1/2.0;
      h_chnum = floor(y/(GVars.h1c+d));  /* 0= lower side, nslit+1=upper side  */
      h_edge  = h_chnum*(GVars.h1c+d);
      if(y - h_edge > GVars.h1c)
      {
        y -= h1/2.0; /* Re-adjust origin */
        ABSORB;
      }
      h_adj = h_edge + (GVars.h1c)/2.0;
      y -= h_adj; h_adj -=  h1/2.0;

      /* neutron is now in the input window of the guide */
      /* do loops on reflections in the box */
      for(;;)
      {
        /* get intersections for all box sides */
        double q, nx,ny,nz;
        double this_length;
        int side=0;

        bounces++;
        /* now look for intersection with guide sides and exit */
        side = Gravity_guide_Trace(&dt, &GVars, x, y, z,
            vx, vy, vz, w_chnum, nslit, h_chnum, nhslit,
            &nx, &ny, &nz);

        /* only positive dt are valid */
        /* exit reflection loops if no intersection (neutron is after box) */
        if (side == 0 || dt <= 0)
          { if (GVars.warnings < 100)
              fprintf(stderr,"%s: warning: neutron has entered guide, but can not exit !\n", GVars.compcurname);
            GVars.warnings++;
            x += w_adj; y += h_adj; ABSORB; } /* should never occur */

        /* propagate to dt */
        PROP_GRAV_DT(dt, GVars.gx, GVars.gy, GVars.gz);

        /* do reflection on speed for l/r/u/d sides */
        if (side == 5) /* neutron reaches end of guide: end loop and exit comp */
          { GVars.N_reflection[side]++; x += w_adj; y += h_adj; SCATTER; x -= w_adj; y -= h_adj; break; }
        /* else reflection on a guide wall */
        if(GVars.M[side] == 0 || Qc == 0 || R0 == 0)  /* walls are absorbing */
          { x += w_adj; y += h_adj; ABSORB; }
        /* handle chamfers */
        this_width = w1+(w2-w1)*z/l;
        this_height= h1+(h2-h1)*z/l;
        this_length= fmod(z, l/nelements);
        /* absorb on input/output of element parts */
        if (GVars.chamfer_z && (this_length<GVars.chamfer_z || this_length>l/nelements-GVars.chamfer_z))
        { x += w_adj; y += h_adj; ABSORB; }
        /* absorb on l/r/t/b sides */
        if (GVars.chamfer_lr && (side==1 || side==2) && (fabs(y+h_adj)>this_height/2-GVars.chamfer_lr))
        { x += w_adj; y += h_adj; ABSORB; }
        if (GVars.chamfer_tb && (side==3 || side==4) && (fabs(x+w_adj)>this_width/2- GVars.chamfer_tb))
        { x += w_adj; y += h_adj; ABSORB; }
        /* change/mirror velocity: h_f = v - n.2*n.v/|n|^2 */
        GVars.N_reflection[side]++; /* GVars.norm_n2 > 0 was checked at INIT */
        /* compute n.v using current values */
        B = scalar_prod(vx,vy,vz,nx,ny,nz);
        dt = 2*B/GVars.norm_n2[side]; /* 2*n.v/|n|^2 */
        vx -= nx*dt;
        vy -= ny*dt;
        vz -= nz*dt;

        /* compute q and modify neutron weight */
        /* scattering q=|n_i-n_f| = V2Q*|vf - v| = V2Q*2*n.v/|n| */
        q = 2*V2Q*fabs(B)/GVars.norm_n[side];

        if (reflect && strlen(reflect) && strcmp(reflect,"NULL") && strcmp(reflect,"0"))
          TableReflecFunc(q, &pTable, &B);
        else {
          double par[] = {R0, Qc, GVars.Alpha[side], GVars.M[side], W};
          StdReflecFunc(q, par, &B);
        }
        if (B <= 0) { x += w_adj; y += h_adj; ABSORB; }
        else p *= B;
        x += w_adj; y += h_adj; SCATTER; x -= w_adj; y -= h_adj;
        GVars.N_reflection[0]++;
        /* go to the next reflection */
        if (bounces > 1000) ABSORB;
      } /* end for */
      x += w_adj; y += h_adj; /* Re-adjust origin after SCATTER */
    }

    if (GVars.fc_freq != 0 || GVars.fc_phase != 0) { /* rotate back neutron w/r to guide element */
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

  } /* if l */
}
#line 20511 "./ISIS_SANS2d.c"
}   /* End of guide_straight3=Guide_gravity() SETTING parameter declarations. */
#undef pTable
#undef GVars
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompguide_straight3:
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

  /* TRACE Component guide_straight4 [21] */
  mccoordschange(mcposrguide_straight4, mcrotrguide_straight4,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component guide_straight4 (without coords transformations) */
  mcJumpTrace_guide_straight4:
  SIG_MESSAGE("guide_straight4 (Trace)");
  mcDEBUG_COMP("guide_straight4")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbCompguide_straight4
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
#define mccompcurname  guide_straight4
#define mccompcurtype  Guide_gravity
#define mccompcurindex 21
#define GVars mccguide_straight4_GVars
#define pTable mccguide_straight4_pTable
{   /* Declarations of guide_straight4=Guide_gravity() SETTING parameters. */
MCNUM w1 = mccguide_straight4_w1;
MCNUM h1 = mccguide_straight4_h1;
MCNUM w2 = mccguide_straight4_w2;
MCNUM h2 = mccguide_straight4_h2;
MCNUM l = mccguide_straight4_l;
MCNUM R0 = mccguide_straight4_R0;
MCNUM Qc = mccguide_straight4_Qc;
MCNUM alpha = mccguide_straight4_alpha;
MCNUM m = mccguide_straight4_m;
MCNUM W = mccguide_straight4_W;
MCNUM nslit = mccguide_straight4_nslit;
MCNUM d = mccguide_straight4_d;
MCNUM mleft = mccguide_straight4_mleft;
MCNUM mright = mccguide_straight4_mright;
MCNUM mtop = mccguide_straight4_mtop;
MCNUM mbottom = mccguide_straight4_mbottom;
MCNUM nhslit = mccguide_straight4_nhslit;
MCNUM G = mccguide_straight4_G;
MCNUM aleft = mccguide_straight4_aleft;
MCNUM aright = mccguide_straight4_aright;
MCNUM atop = mccguide_straight4_atop;
MCNUM abottom = mccguide_straight4_abottom;
MCNUM wavy = mccguide_straight4_wavy;
MCNUM wavy_z = mccguide_straight4_wavy_z;
MCNUM wavy_tb = mccguide_straight4_wavy_tb;
MCNUM wavy_lr = mccguide_straight4_wavy_lr;
MCNUM chamfers = mccguide_straight4_chamfers;
MCNUM chamfers_z = mccguide_straight4_chamfers_z;
MCNUM chamfers_lr = mccguide_straight4_chamfers_lr;
MCNUM chamfers_tb = mccguide_straight4_chamfers_tb;
MCNUM nelements = mccguide_straight4_nelements;
MCNUM nu = mccguide_straight4_nu;
MCNUM phase = mccguide_straight4_phase;
char* reflect = mccguide_straight4_reflect;
#line 392 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide_gravity.comp"
{
  if (l > 0 && nelements > 0) {
    double B, C, dt;
    int    ret, bounces = 0, i=0;
    double this_width, this_height;
    double angle=0;

    if (GVars.fc_freq != 0 || GVars.fc_phase != 0) { /* rotate neutron w/r to guide element */
      /* approximation of rotating straight Fermi Chopper */
      Coords   X = coords_set(x,y,z-l/2);  /* current coordinates of neutron in centered static frame */
      Rotation R;
      double dt=(-z+l/2)/vz; /* time shift to each center of slit package */
      angle=fmod(360*GVars.fc_freq*(t+dt)+GVars.fc_phase, 360); /* in deg */
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

    for (i=0; i<7; GVars.N_reflection[i++] = 0);

    /* propagate to box input (with gravitation) in comp local coords */
    /* A = 0.5 n.g; B = n.v; C = n.(r-W); */
    /* 0=Z0 side: n=(0, 0, -l) ; W = (0, 0, 0) (at z=0, guide input)*/
    B = -l*vz; C = -l*z;

    ret = solve_2nd_order(&dt, NULL, GVars.A[0], B, C);
    if (ret==0) ABSORB;

    if (dt>0.0) PROP_GRAV_DT(dt, GVars.gx, GVars.gy, GVars.gz); else if (angle) ABSORB;
    GVars.N_reflection[6]++;

    this_width  = w1;
    this_height = h1;

  /* check if we are in the box input, else absorb */
    if (fabs(x) > this_width/2 || fabs(y) > this_height/2)
      ABSORB;
    else
    {
      double w_edge, w_adj; /* Channel displacement on X */
      double h_edge, h_adj; /* Channel displacement on Y */
      double w_chnum,h_chnum; /* channel indexes */

      SCATTER;

      /* X: Shift origin to center of channel hit (absorb if hit dividing walls) */
      x += w1/2.0;
      w_chnum = floor(x/(GVars.w1c+d));  /* 0= right side, nslit+1=left side  */
      w_edge  = w_chnum*(GVars.w1c+d);
      if(x - w_edge > GVars.w1c)
      {
        x -= w1/2.0; /* Re-adjust origin */
        ABSORB;
      }
      w_adj = w_edge + (GVars.w1c)/2.0;
      x -= w_adj; w_adj -=  w1/2.0;

      /* Y: Shift origin to center of channel hit (absorb if hit dividing walls) */
      y += h1/2.0;
      h_chnum = floor(y/(GVars.h1c+d));  /* 0= lower side, nslit+1=upper side  */
      h_edge  = h_chnum*(GVars.h1c+d);
      if(y - h_edge > GVars.h1c)
      {
        y -= h1/2.0; /* Re-adjust origin */
        ABSORB;
      }
      h_adj = h_edge + (GVars.h1c)/2.0;
      y -= h_adj; h_adj -=  h1/2.0;

      /* neutron is now in the input window of the guide */
      /* do loops on reflections in the box */
      for(;;)
      {
        /* get intersections for all box sides */
        double q, nx,ny,nz;
        double this_length;
        int side=0;

        bounces++;
        /* now look for intersection with guide sides and exit */
        side = Gravity_guide_Trace(&dt, &GVars, x, y, z,
            vx, vy, vz, w_chnum, nslit, h_chnum, nhslit,
            &nx, &ny, &nz);

        /* only positive dt are valid */
        /* exit reflection loops if no intersection (neutron is after box) */
        if (side == 0 || dt <= 0)
          { if (GVars.warnings < 100)
              fprintf(stderr,"%s: warning: neutron has entered guide, but can not exit !\n", GVars.compcurname);
            GVars.warnings++;
            x += w_adj; y += h_adj; ABSORB; } /* should never occur */

        /* propagate to dt */
        PROP_GRAV_DT(dt, GVars.gx, GVars.gy, GVars.gz);

        /* do reflection on speed for l/r/u/d sides */
        if (side == 5) /* neutron reaches end of guide: end loop and exit comp */
          { GVars.N_reflection[side]++; x += w_adj; y += h_adj; SCATTER; x -= w_adj; y -= h_adj; break; }
        /* else reflection on a guide wall */
        if(GVars.M[side] == 0 || Qc == 0 || R0 == 0)  /* walls are absorbing */
          { x += w_adj; y += h_adj; ABSORB; }
        /* handle chamfers */
        this_width = w1+(w2-w1)*z/l;
        this_height= h1+(h2-h1)*z/l;
        this_length= fmod(z, l/nelements);
        /* absorb on input/output of element parts */
        if (GVars.chamfer_z && (this_length<GVars.chamfer_z || this_length>l/nelements-GVars.chamfer_z))
        { x += w_adj; y += h_adj; ABSORB; }
        /* absorb on l/r/t/b sides */
        if (GVars.chamfer_lr && (side==1 || side==2) && (fabs(y+h_adj)>this_height/2-GVars.chamfer_lr))
        { x += w_adj; y += h_adj; ABSORB; }
        if (GVars.chamfer_tb && (side==3 || side==4) && (fabs(x+w_adj)>this_width/2- GVars.chamfer_tb))
        { x += w_adj; y += h_adj; ABSORB; }
        /* change/mirror velocity: h_f = v - n.2*n.v/|n|^2 */
        GVars.N_reflection[side]++; /* GVars.norm_n2 > 0 was checked at INIT */
        /* compute n.v using current values */
        B = scalar_prod(vx,vy,vz,nx,ny,nz);
        dt = 2*B/GVars.norm_n2[side]; /* 2*n.v/|n|^2 */
        vx -= nx*dt;
        vy -= ny*dt;
        vz -= nz*dt;

        /* compute q and modify neutron weight */
        /* scattering q=|n_i-n_f| = V2Q*|vf - v| = V2Q*2*n.v/|n| */
        q = 2*V2Q*fabs(B)/GVars.norm_n[side];

        if (reflect && strlen(reflect) && strcmp(reflect,"NULL") && strcmp(reflect,"0"))
          TableReflecFunc(q, &pTable, &B);
        else {
          double par[] = {R0, Qc, GVars.Alpha[side], GVars.M[side], W};
          StdReflecFunc(q, par, &B);
        }
        if (B <= 0) { x += w_adj; y += h_adj; ABSORB; }
        else p *= B;
        x += w_adj; y += h_adj; SCATTER; x -= w_adj; y -= h_adj;
        GVars.N_reflection[0]++;
        /* go to the next reflection */
        if (bounces > 1000) ABSORB;
      } /* end for */
      x += w_adj; y += h_adj; /* Re-adjust origin after SCATTER */
    }

    if (GVars.fc_freq != 0 || GVars.fc_phase != 0) { /* rotate back neutron w/r to guide element */
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

  } /* if l */
}
#line 20824 "./ISIS_SANS2d.c"
}   /* End of guide_straight4=Guide_gravity() SETTING parameter declarations. */
#undef pTable
#undef GVars
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompguide_straight4:
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

  /* TRACE Component psd3 [22] */
  mccoordschange(mcposrpsd3, mcrotrpsd3,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component psd3 (without coords transformations) */
  mcJumpTrace_psd3:
  SIG_MESSAGE("psd3 (Trace)");
  mcDEBUG_COMP("psd3")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbComppsd3
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
#define mccompcurname  psd3
#define mccompcurtype  PSD_monitor
#define mccompcurindex 22
#define PSD_N mccpsd3_PSD_N
#define PSD_p mccpsd3_PSD_p
#define PSD_p2 mccpsd3_PSD_p2
{   /* Declarations of psd3=PSD_monitor() SETTING parameters. */
int nx = mccpsd3_nx;
int ny = mccpsd3_ny;
char* filename = mccpsd3_filename;
MCNUM xmin = mccpsd3_xmin;
MCNUM xmax = mccpsd3_xmax;
MCNUM ymin = mccpsd3_ymin;
MCNUM ymax = mccpsd3_ymax;
MCNUM xwidth = mccpsd3_xwidth;
MCNUM yheight = mccpsd3_yheight;
MCNUM restore_neutron = mccpsd3_restore_neutron;
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
#line 20960 "./ISIS_SANS2d.c"
}   /* End of psd3=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbComppsd3:
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

  /* TRACE Component aperture1 [23] */
  mccoordschange(mcposraperture1, mcrotraperture1,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component aperture1 (without coords transformations) */
  mcJumpTrace_aperture1:
  SIG_MESSAGE("aperture1 (Trace)");
  mcDEBUG_COMP("aperture1")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbCompaperture1
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
#define mccompcurname  aperture1
#define mccompcurtype  Slit
#define mccompcurindex 23
{   /* Declarations of aperture1=Slit() SETTING parameters. */
MCNUM xmin = mccaperture1_xmin;
MCNUM xmax = mccaperture1_xmax;
MCNUM ymin = mccaperture1_ymin;
MCNUM ymax = mccaperture1_ymax;
MCNUM radius = mccaperture1_radius;
MCNUM xwidth = mccaperture1_xwidth;
MCNUM yheight = mccaperture1_yheight;
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
#line 21087 "./ISIS_SANS2d.c"
}   /* End of aperture1=Slit() SETTING parameter declarations. */
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompaperture1:
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

  /* TRACE Component lmonitor2 [24] */
  mccoordschange(mcposrlmonitor2, mcrotrlmonitor2,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component lmonitor2 (without coords transformations) */
  mcJumpTrace_lmonitor2:
  SIG_MESSAGE("lmonitor2 (Trace)");
  mcDEBUG_COMP("lmonitor2")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbComplmonitor2
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
#define mccompcurname  lmonitor2
#define mccompcurtype  L_monitor
#define mccompcurindex 24
#define nL mcclmonitor2_nL
#define L_N mcclmonitor2_L_N
#define L_p mcclmonitor2_L_p
#define L_p2 mcclmonitor2_L_p2
{   /* Declarations of lmonitor2=L_monitor() SETTING parameters. */
char* filename = mcclmonitor2_filename;
MCNUM xmin = mcclmonitor2_xmin;
MCNUM xmax = mcclmonitor2_xmax;
MCNUM ymin = mcclmonitor2_ymin;
MCNUM ymax = mcclmonitor2_ymax;
MCNUM xwidth = mcclmonitor2_xwidth;
MCNUM yheight = mcclmonitor2_yheight;
MCNUM Lmin = mcclmonitor2_Lmin;
MCNUM Lmax = mcclmonitor2_Lmax;
MCNUM restore_neutron = mcclmonitor2_restore_neutron;
int nowritefile = mcclmonitor2_nowritefile;
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
#line 21230 "./ISIS_SANS2d.c"
}   /* End of lmonitor2=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbComplmonitor2:
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

  /* TRACE Component S6 [25] */
  mccoordschange(mcposrS6, mcrotrS6,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component S6 (without coords transformations) */
  mcJumpTrace_S6:
  SIG_MESSAGE("S6 (Trace)");
  mcDEBUG_COMP("S6")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbCompS6
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
#define mccompcurname  S6
#define mccompcurtype  Slit
#define mccompcurindex 25
{   /* Declarations of S6=Slit() SETTING parameters. */
MCNUM xmin = mccS6_xmin;
MCNUM xmax = mccS6_xmax;
MCNUM ymin = mccS6_ymin;
MCNUM ymax = mccS6_ymax;
MCNUM radius = mccS6_radius;
MCNUM xwidth = mccS6_xwidth;
MCNUM yheight = mccS6_yheight;
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
#line 21358 "./ISIS_SANS2d.c"
}   /* End of S6=Slit() SETTING parameter declarations. */
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompS6:
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

  /* TRACE Component APERTURE2 [26] */
  mccoordschange(mcposrAPERTURE2, mcrotrAPERTURE2,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component APERTURE2 (without coords transformations) */
  mcJumpTrace_APERTURE2:
  SIG_MESSAGE("APERTURE2 (Trace)");
  mcDEBUG_COMP("APERTURE2")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbCompAPERTURE2
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
#define mccompcurname  APERTURE2
#define mccompcurtype  Slit
#define mccompcurindex 26
{   /* Declarations of APERTURE2=Slit() SETTING parameters. */
MCNUM xmin = mccAPERTURE2_xmin;
MCNUM xmax = mccAPERTURE2_xmax;
MCNUM ymin = mccAPERTURE2_ymin;
MCNUM ymax = mccAPERTURE2_ymax;
MCNUM radius = mccAPERTURE2_radius;
MCNUM xwidth = mccAPERTURE2_xwidth;
MCNUM yheight = mccAPERTURE2_yheight;
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
#line 21482 "./ISIS_SANS2d.c"
}   /* End of APERTURE2=Slit() SETTING parameter declarations. */
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompAPERTURE2:
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

  /* TRACE Component lmon2 [27] */
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
#define mccompcurname  lmon2
#define mccompcurtype  L_monitor
#define mccompcurindex 27
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
#line 21625 "./ISIS_SANS2d.c"
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

  /* TRACE Component psd4 [28] */
  mccoordschange(mcposrpsd4, mcrotrpsd4,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component psd4 (without coords transformations) */
  mcJumpTrace_psd4:
  SIG_MESSAGE("psd4 (Trace)");
  mcDEBUG_COMP("psd4")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbComppsd4
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
#define mccompcurname  psd4
#define mccompcurtype  PSD_monitor
#define mccompcurindex 28
#define PSD_N mccpsd4_PSD_N
#define PSD_p mccpsd4_PSD_p
#define PSD_p2 mccpsd4_PSD_p2
{   /* Declarations of psd4=PSD_monitor() SETTING parameters. */
int nx = mccpsd4_nx;
int ny = mccpsd4_ny;
char* filename = mccpsd4_filename;
MCNUM xmin = mccpsd4_xmin;
MCNUM xmax = mccpsd4_xmax;
MCNUM ymin = mccpsd4_ymin;
MCNUM ymax = mccpsd4_ymax;
MCNUM xwidth = mccpsd4_xwidth;
MCNUM yheight = mccpsd4_yheight;
MCNUM restore_neutron = mccpsd4_restore_neutron;
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
#line 21763 "./ISIS_SANS2d.c"
}   /* End of psd4=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbComppsd4:
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
#define mccompcurindex 3
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
#line 21877 "./ISIS_SANS2d.c"
}   /* End of lmon1=L_monitor() SETTING parameter declarations. */
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
#define mccompcurindex 4
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
#line 21917 "./ISIS_SANS2d.c"
}   /* End of psd1=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'lmonb'. */
  SIG_MESSAGE("lmonb (Save)");
#define mccompcurname  lmonb
#define mccompcurtype  L_monitor
#define mccompcurindex 15
#define nL mcclmonb_nL
#define L_N mcclmonb_L_N
#define L_p mcclmonb_L_p
#define L_p2 mcclmonb_L_p2
{   /* Declarations of lmonb=L_monitor() SETTING parameters. */
char* filename = mcclmonb_filename;
MCNUM xmin = mcclmonb_xmin;
MCNUM xmax = mcclmonb_xmax;
MCNUM ymin = mcclmonb_ymin;
MCNUM ymax = mcclmonb_ymax;
MCNUM xwidth = mcclmonb_xwidth;
MCNUM yheight = mcclmonb_yheight;
MCNUM Lmin = mcclmonb_Lmin;
MCNUM Lmax = mcclmonb_Lmax;
MCNUM restore_neutron = mcclmonb_restore_neutron;
int nowritefile = mcclmonb_nowritefile;
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
#line 21959 "./ISIS_SANS2d.c"
}   /* End of lmonb=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'psd2'. */
  SIG_MESSAGE("psd2 (Save)");
#define mccompcurname  psd2
#define mccompcurtype  PSD_monitor
#define mccompcurindex 16
#define PSD_N mccpsd2_PSD_N
#define PSD_p mccpsd2_PSD_p
#define PSD_p2 mccpsd2_PSD_p2
{   /* Declarations of psd2=PSD_monitor() SETTING parameters. */
int nx = mccpsd2_nx;
int ny = mccpsd2_ny;
char* filename = mccpsd2_filename;
MCNUM xmin = mccpsd2_xmin;
MCNUM xmax = mccpsd2_xmax;
MCNUM ymin = mccpsd2_ymin;
MCNUM ymax = mccpsd2_ymax;
MCNUM xwidth = mccpsd2_xwidth;
MCNUM yheight = mccpsd2_yheight;
MCNUM restore_neutron = mccpsd2_restore_neutron;
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
#line 21999 "./ISIS_SANS2d.c"
}   /* End of psd2=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'psd3'. */
  SIG_MESSAGE("psd3 (Save)");
#define mccompcurname  psd3
#define mccompcurtype  PSD_monitor
#define mccompcurindex 22
#define PSD_N mccpsd3_PSD_N
#define PSD_p mccpsd3_PSD_p
#define PSD_p2 mccpsd3_PSD_p2
{   /* Declarations of psd3=PSD_monitor() SETTING parameters. */
int nx = mccpsd3_nx;
int ny = mccpsd3_ny;
char* filename = mccpsd3_filename;
MCNUM xmin = mccpsd3_xmin;
MCNUM xmax = mccpsd3_xmax;
MCNUM ymin = mccpsd3_ymin;
MCNUM ymax = mccpsd3_ymax;
MCNUM xwidth = mccpsd3_xwidth;
MCNUM yheight = mccpsd3_yheight;
MCNUM restore_neutron = mccpsd3_restore_neutron;
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
#line 22038 "./ISIS_SANS2d.c"
}   /* End of psd3=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'lmonitor2'. */
  SIG_MESSAGE("lmonitor2 (Save)");
#define mccompcurname  lmonitor2
#define mccompcurtype  L_monitor
#define mccompcurindex 24
#define nL mcclmonitor2_nL
#define L_N mcclmonitor2_L_N
#define L_p mcclmonitor2_L_p
#define L_p2 mcclmonitor2_L_p2
{   /* Declarations of lmonitor2=L_monitor() SETTING parameters. */
char* filename = mcclmonitor2_filename;
MCNUM xmin = mcclmonitor2_xmin;
MCNUM xmax = mcclmonitor2_xmax;
MCNUM ymin = mcclmonitor2_ymin;
MCNUM ymax = mcclmonitor2_ymax;
MCNUM xwidth = mcclmonitor2_xwidth;
MCNUM yheight = mcclmonitor2_yheight;
MCNUM Lmin = mcclmonitor2_Lmin;
MCNUM Lmax = mcclmonitor2_Lmax;
MCNUM restore_neutron = mcclmonitor2_restore_neutron;
int nowritefile = mcclmonitor2_nowritefile;
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
#line 22080 "./ISIS_SANS2d.c"
}   /* End of lmonitor2=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'lmon2'. */
  SIG_MESSAGE("lmon2 (Save)");
#define mccompcurname  lmon2
#define mccompcurtype  L_monitor
#define mccompcurindex 27
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
#line 22123 "./ISIS_SANS2d.c"
}   /* End of lmon2=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'psd4'. */
  SIG_MESSAGE("psd4 (Save)");
#define mccompcurname  psd4
#define mccompcurtype  PSD_monitor
#define mccompcurindex 28
#define PSD_N mccpsd4_PSD_N
#define PSD_p mccpsd4_PSD_p
#define PSD_p2 mccpsd4_PSD_p2
{   /* Declarations of psd4=PSD_monitor() SETTING parameters. */
int nx = mccpsd4_nx;
int ny = mccpsd4_ny;
char* filename = mccpsd4_filename;
MCNUM xmin = mccpsd4_xmin;
MCNUM xmax = mccpsd4_xmax;
MCNUM ymin = mccpsd4_ymin;
MCNUM ymax = mccpsd4_ymax;
MCNUM xwidth = mccpsd4_xwidth;
MCNUM yheight = mccpsd4_yheight;
MCNUM restore_neutron = mccpsd4_restore_neutron;
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
#line 22163 "./ISIS_SANS2d.c"
}   /* End of psd4=PSD_monitor() SETTING parameter declarations. */
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

    if (!mcNCounter[1]) fprintf(stderr, "Warning: No neutron could reach Component[1] Origin\n");
    if (mcAbsorbProp[1]) fprintf(stderr, "Warning: %g events were removed in Component[1] Origin=Arm()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[1]);
    if (!mcNCounter[2]) fprintf(stderr, "Warning: No neutron could reach Component[2] isis_source\n");
    if (mcAbsorbProp[2]) fprintf(stderr, "Warning: %g events were removed in Component[2] isis_source=ISIS_moderator()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[2]);
    if (!mcNCounter[3]) fprintf(stderr, "Warning: No neutron could reach Component[3] lmon1\n");
    if (mcAbsorbProp[3]) fprintf(stderr, "Warning: %g events were removed in Component[3] lmon1=L_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[3]);
  /* User FINALLY code for component 'psd1'. */
  SIG_MESSAGE("psd1 (Finally)");
#define mccompcurname  psd1
#define mccompcurtype  PSD_monitor
#define mccompcurindex 4
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
#line 22210 "./ISIS_SANS2d.c"
}   /* End of psd1=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[4]) fprintf(stderr, "Warning: No neutron could reach Component[4] psd1\n");
    if (mcAbsorbProp[4]) fprintf(stderr, "Warning: %g events were removed in Component[4] psd1=PSD_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[4]);
  /* User FINALLY code for component 'bender1'. */
  SIG_MESSAGE("bender1 (Finally)");
#define mccompcurname  bender1
#define mccompcurtype  Guide_gravity
#define mccompcurindex 5
#define GVars mccbender1_GVars
#define pTable mccbender1_pTable
{   /* Declarations of bender1=Guide_gravity() SETTING parameters. */
MCNUM w1 = mccbender1_w1;
MCNUM h1 = mccbender1_h1;
MCNUM w2 = mccbender1_w2;
MCNUM h2 = mccbender1_h2;
MCNUM l = mccbender1_l;
MCNUM R0 = mccbender1_R0;
MCNUM Qc = mccbender1_Qc;
MCNUM alpha = mccbender1_alpha;
MCNUM m = mccbender1_m;
MCNUM W = mccbender1_W;
MCNUM nslit = mccbender1_nslit;
MCNUM d = mccbender1_d;
MCNUM mleft = mccbender1_mleft;
MCNUM mright = mccbender1_mright;
MCNUM mtop = mccbender1_mtop;
MCNUM mbottom = mccbender1_mbottom;
MCNUM nhslit = mccbender1_nhslit;
MCNUM G = mccbender1_G;
MCNUM aleft = mccbender1_aleft;
MCNUM aright = mccbender1_aright;
MCNUM atop = mccbender1_atop;
MCNUM abottom = mccbender1_abottom;
MCNUM wavy = mccbender1_wavy;
MCNUM wavy_z = mccbender1_wavy_z;
MCNUM wavy_tb = mccbender1_wavy_tb;
MCNUM wavy_lr = mccbender1_wavy_lr;
MCNUM chamfers = mccbender1_chamfers;
MCNUM chamfers_z = mccbender1_chamfers_z;
MCNUM chamfers_lr = mccbender1_chamfers_lr;
MCNUM chamfers_tb = mccbender1_chamfers_tb;
MCNUM nelements = mccbender1_nelements;
MCNUM nu = mccbender1_nu;
MCNUM phase = mccbender1_phase;
char* reflect = mccbender1_reflect;
#line 562 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide_gravity.comp"
{
if (GVars.warnings > 100) {
  fprintf(stderr,"%s: warning: neutron has entered guide, but can not exit !\n", GVars.compcurname);
  fprintf(stderr,"%s: warning: This message has been repeated %g times\n", GVars.compcurname, GVars.warnings);
}
}
#line 22270 "./ISIS_SANS2d.c"
}   /* End of bender1=Guide_gravity() SETTING parameter declarations. */
#undef pTable
#undef GVars
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[5]) fprintf(stderr, "Warning: No neutron could reach Component[5] bender1\n");
    if (mcAbsorbProp[5]) fprintf(stderr, "Warning: %g events were removed in Component[5] bender1=Guide_gravity()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[5]);
  /* User FINALLY code for component 'bender2'. */
  SIG_MESSAGE("bender2 (Finally)");
#define mccompcurname  bender2
#define mccompcurtype  Guide_gravity
#define mccompcurindex 6
#define GVars mccbender2_GVars
#define pTable mccbender2_pTable
{   /* Declarations of bender2=Guide_gravity() SETTING parameters. */
MCNUM w1 = mccbender2_w1;
MCNUM h1 = mccbender2_h1;
MCNUM w2 = mccbender2_w2;
MCNUM h2 = mccbender2_h2;
MCNUM l = mccbender2_l;
MCNUM R0 = mccbender2_R0;
MCNUM Qc = mccbender2_Qc;
MCNUM alpha = mccbender2_alpha;
MCNUM m = mccbender2_m;
MCNUM W = mccbender2_W;
MCNUM nslit = mccbender2_nslit;
MCNUM d = mccbender2_d;
MCNUM mleft = mccbender2_mleft;
MCNUM mright = mccbender2_mright;
MCNUM mtop = mccbender2_mtop;
MCNUM mbottom = mccbender2_mbottom;
MCNUM nhslit = mccbender2_nhslit;
MCNUM G = mccbender2_G;
MCNUM aleft = mccbender2_aleft;
MCNUM aright = mccbender2_aright;
MCNUM atop = mccbender2_atop;
MCNUM abottom = mccbender2_abottom;
MCNUM wavy = mccbender2_wavy;
MCNUM wavy_z = mccbender2_wavy_z;
MCNUM wavy_tb = mccbender2_wavy_tb;
MCNUM wavy_lr = mccbender2_wavy_lr;
MCNUM chamfers = mccbender2_chamfers;
MCNUM chamfers_z = mccbender2_chamfers_z;
MCNUM chamfers_lr = mccbender2_chamfers_lr;
MCNUM chamfers_tb = mccbender2_chamfers_tb;
MCNUM nelements = mccbender2_nelements;
MCNUM nu = mccbender2_nu;
MCNUM phase = mccbender2_phase;
char* reflect = mccbender2_reflect;
#line 562 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide_gravity.comp"
{
if (GVars.warnings > 100) {
  fprintf(stderr,"%s: warning: neutron has entered guide, but can not exit !\n", GVars.compcurname);
  fprintf(stderr,"%s: warning: This message has been repeated %g times\n", GVars.compcurname, GVars.warnings);
}
}
#line 22329 "./ISIS_SANS2d.c"
}   /* End of bender2=Guide_gravity() SETTING parameter declarations. */
#undef pTable
#undef GVars
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[6]) fprintf(stderr, "Warning: No neutron could reach Component[6] bender2\n");
    if (mcAbsorbProp[6]) fprintf(stderr, "Warning: %g events were removed in Component[6] bender2=Guide_gravity()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[6]);
  /* User FINALLY code for component 'bender3'. */
  SIG_MESSAGE("bender3 (Finally)");
#define mccompcurname  bender3
#define mccompcurtype  Guide_gravity
#define mccompcurindex 7
#define GVars mccbender3_GVars
#define pTable mccbender3_pTable
{   /* Declarations of bender3=Guide_gravity() SETTING parameters. */
MCNUM w1 = mccbender3_w1;
MCNUM h1 = mccbender3_h1;
MCNUM w2 = mccbender3_w2;
MCNUM h2 = mccbender3_h2;
MCNUM l = mccbender3_l;
MCNUM R0 = mccbender3_R0;
MCNUM Qc = mccbender3_Qc;
MCNUM alpha = mccbender3_alpha;
MCNUM m = mccbender3_m;
MCNUM W = mccbender3_W;
MCNUM nslit = mccbender3_nslit;
MCNUM d = mccbender3_d;
MCNUM mleft = mccbender3_mleft;
MCNUM mright = mccbender3_mright;
MCNUM mtop = mccbender3_mtop;
MCNUM mbottom = mccbender3_mbottom;
MCNUM nhslit = mccbender3_nhslit;
MCNUM G = mccbender3_G;
MCNUM aleft = mccbender3_aleft;
MCNUM aright = mccbender3_aright;
MCNUM atop = mccbender3_atop;
MCNUM abottom = mccbender3_abottom;
MCNUM wavy = mccbender3_wavy;
MCNUM wavy_z = mccbender3_wavy_z;
MCNUM wavy_tb = mccbender3_wavy_tb;
MCNUM wavy_lr = mccbender3_wavy_lr;
MCNUM chamfers = mccbender3_chamfers;
MCNUM chamfers_z = mccbender3_chamfers_z;
MCNUM chamfers_lr = mccbender3_chamfers_lr;
MCNUM chamfers_tb = mccbender3_chamfers_tb;
MCNUM nelements = mccbender3_nelements;
MCNUM nu = mccbender3_nu;
MCNUM phase = mccbender3_phase;
char* reflect = mccbender3_reflect;
#line 562 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide_gravity.comp"
{
if (GVars.warnings > 100) {
  fprintf(stderr,"%s: warning: neutron has entered guide, but can not exit !\n", GVars.compcurname);
  fprintf(stderr,"%s: warning: This message has been repeated %g times\n", GVars.compcurname, GVars.warnings);
}
}
#line 22388 "./ISIS_SANS2d.c"
}   /* End of bender3=Guide_gravity() SETTING parameter declarations. */
#undef pTable
#undef GVars
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[7]) fprintf(stderr, "Warning: No neutron could reach Component[7] bender3\n");
    if (mcAbsorbProp[7]) fprintf(stderr, "Warning: %g events were removed in Component[7] bender3=Guide_gravity()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[7]);
  /* User FINALLY code for component 'bender4'. */
  SIG_MESSAGE("bender4 (Finally)");
#define mccompcurname  bender4
#define mccompcurtype  Guide_gravity
#define mccompcurindex 8
#define GVars mccbender4_GVars
#define pTable mccbender4_pTable
{   /* Declarations of bender4=Guide_gravity() SETTING parameters. */
MCNUM w1 = mccbender4_w1;
MCNUM h1 = mccbender4_h1;
MCNUM w2 = mccbender4_w2;
MCNUM h2 = mccbender4_h2;
MCNUM l = mccbender4_l;
MCNUM R0 = mccbender4_R0;
MCNUM Qc = mccbender4_Qc;
MCNUM alpha = mccbender4_alpha;
MCNUM m = mccbender4_m;
MCNUM W = mccbender4_W;
MCNUM nslit = mccbender4_nslit;
MCNUM d = mccbender4_d;
MCNUM mleft = mccbender4_mleft;
MCNUM mright = mccbender4_mright;
MCNUM mtop = mccbender4_mtop;
MCNUM mbottom = mccbender4_mbottom;
MCNUM nhslit = mccbender4_nhslit;
MCNUM G = mccbender4_G;
MCNUM aleft = mccbender4_aleft;
MCNUM aright = mccbender4_aright;
MCNUM atop = mccbender4_atop;
MCNUM abottom = mccbender4_abottom;
MCNUM wavy = mccbender4_wavy;
MCNUM wavy_z = mccbender4_wavy_z;
MCNUM wavy_tb = mccbender4_wavy_tb;
MCNUM wavy_lr = mccbender4_wavy_lr;
MCNUM chamfers = mccbender4_chamfers;
MCNUM chamfers_z = mccbender4_chamfers_z;
MCNUM chamfers_lr = mccbender4_chamfers_lr;
MCNUM chamfers_tb = mccbender4_chamfers_tb;
MCNUM nelements = mccbender4_nelements;
MCNUM nu = mccbender4_nu;
MCNUM phase = mccbender4_phase;
char* reflect = mccbender4_reflect;
#line 562 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide_gravity.comp"
{
if (GVars.warnings > 100) {
  fprintf(stderr,"%s: warning: neutron has entered guide, but can not exit !\n", GVars.compcurname);
  fprintf(stderr,"%s: warning: This message has been repeated %g times\n", GVars.compcurname, GVars.warnings);
}
}
#line 22447 "./ISIS_SANS2d.c"
}   /* End of bender4=Guide_gravity() SETTING parameter declarations. */
#undef pTable
#undef GVars
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[8]) fprintf(stderr, "Warning: No neutron could reach Component[8] bender4\n");
    if (mcAbsorbProp[8]) fprintf(stderr, "Warning: %g events were removed in Component[8] bender4=Guide_gravity()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[8]);
  /* User FINALLY code for component 'bender5'. */
  SIG_MESSAGE("bender5 (Finally)");
#define mccompcurname  bender5
#define mccompcurtype  Guide_gravity
#define mccompcurindex 9
#define GVars mccbender5_GVars
#define pTable mccbender5_pTable
{   /* Declarations of bender5=Guide_gravity() SETTING parameters. */
MCNUM w1 = mccbender5_w1;
MCNUM h1 = mccbender5_h1;
MCNUM w2 = mccbender5_w2;
MCNUM h2 = mccbender5_h2;
MCNUM l = mccbender5_l;
MCNUM R0 = mccbender5_R0;
MCNUM Qc = mccbender5_Qc;
MCNUM alpha = mccbender5_alpha;
MCNUM m = mccbender5_m;
MCNUM W = mccbender5_W;
MCNUM nslit = mccbender5_nslit;
MCNUM d = mccbender5_d;
MCNUM mleft = mccbender5_mleft;
MCNUM mright = mccbender5_mright;
MCNUM mtop = mccbender5_mtop;
MCNUM mbottom = mccbender5_mbottom;
MCNUM nhslit = mccbender5_nhslit;
MCNUM G = mccbender5_G;
MCNUM aleft = mccbender5_aleft;
MCNUM aright = mccbender5_aright;
MCNUM atop = mccbender5_atop;
MCNUM abottom = mccbender5_abottom;
MCNUM wavy = mccbender5_wavy;
MCNUM wavy_z = mccbender5_wavy_z;
MCNUM wavy_tb = mccbender5_wavy_tb;
MCNUM wavy_lr = mccbender5_wavy_lr;
MCNUM chamfers = mccbender5_chamfers;
MCNUM chamfers_z = mccbender5_chamfers_z;
MCNUM chamfers_lr = mccbender5_chamfers_lr;
MCNUM chamfers_tb = mccbender5_chamfers_tb;
MCNUM nelements = mccbender5_nelements;
MCNUM nu = mccbender5_nu;
MCNUM phase = mccbender5_phase;
char* reflect = mccbender5_reflect;
#line 562 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide_gravity.comp"
{
if (GVars.warnings > 100) {
  fprintf(stderr,"%s: warning: neutron has entered guide, but can not exit !\n", GVars.compcurname);
  fprintf(stderr,"%s: warning: This message has been repeated %g times\n", GVars.compcurname, GVars.warnings);
}
}
#line 22506 "./ISIS_SANS2d.c"
}   /* End of bender5=Guide_gravity() SETTING parameter declarations. */
#undef pTable
#undef GVars
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[9]) fprintf(stderr, "Warning: No neutron could reach Component[9] bender5\n");
    if (mcAbsorbProp[9]) fprintf(stderr, "Warning: %g events were removed in Component[9] bender5=Guide_gravity()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[9]);
  /* User FINALLY code for component 'bender6'. */
  SIG_MESSAGE("bender6 (Finally)");
#define mccompcurname  bender6
#define mccompcurtype  Guide_gravity
#define mccompcurindex 10
#define GVars mccbender6_GVars
#define pTable mccbender6_pTable
{   /* Declarations of bender6=Guide_gravity() SETTING parameters. */
MCNUM w1 = mccbender6_w1;
MCNUM h1 = mccbender6_h1;
MCNUM w2 = mccbender6_w2;
MCNUM h2 = mccbender6_h2;
MCNUM l = mccbender6_l;
MCNUM R0 = mccbender6_R0;
MCNUM Qc = mccbender6_Qc;
MCNUM alpha = mccbender6_alpha;
MCNUM m = mccbender6_m;
MCNUM W = mccbender6_W;
MCNUM nslit = mccbender6_nslit;
MCNUM d = mccbender6_d;
MCNUM mleft = mccbender6_mleft;
MCNUM mright = mccbender6_mright;
MCNUM mtop = mccbender6_mtop;
MCNUM mbottom = mccbender6_mbottom;
MCNUM nhslit = mccbender6_nhslit;
MCNUM G = mccbender6_G;
MCNUM aleft = mccbender6_aleft;
MCNUM aright = mccbender6_aright;
MCNUM atop = mccbender6_atop;
MCNUM abottom = mccbender6_abottom;
MCNUM wavy = mccbender6_wavy;
MCNUM wavy_z = mccbender6_wavy_z;
MCNUM wavy_tb = mccbender6_wavy_tb;
MCNUM wavy_lr = mccbender6_wavy_lr;
MCNUM chamfers = mccbender6_chamfers;
MCNUM chamfers_z = mccbender6_chamfers_z;
MCNUM chamfers_lr = mccbender6_chamfers_lr;
MCNUM chamfers_tb = mccbender6_chamfers_tb;
MCNUM nelements = mccbender6_nelements;
MCNUM nu = mccbender6_nu;
MCNUM phase = mccbender6_phase;
char* reflect = mccbender6_reflect;
#line 562 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide_gravity.comp"
{
if (GVars.warnings > 100) {
  fprintf(stderr,"%s: warning: neutron has entered guide, but can not exit !\n", GVars.compcurname);
  fprintf(stderr,"%s: warning: This message has been repeated %g times\n", GVars.compcurname, GVars.warnings);
}
}
#line 22565 "./ISIS_SANS2d.c"
}   /* End of bender6=Guide_gravity() SETTING parameter declarations. */
#undef pTable
#undef GVars
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[10]) fprintf(stderr, "Warning: No neutron could reach Component[10] bender6\n");
    if (mcAbsorbProp[10]) fprintf(stderr, "Warning: %g events were removed in Component[10] bender6=Guide_gravity()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[10]);
  /* User FINALLY code for component 'bender7'. */
  SIG_MESSAGE("bender7 (Finally)");
#define mccompcurname  bender7
#define mccompcurtype  Guide_gravity
#define mccompcurindex 11
#define GVars mccbender7_GVars
#define pTable mccbender7_pTable
{   /* Declarations of bender7=Guide_gravity() SETTING parameters. */
MCNUM w1 = mccbender7_w1;
MCNUM h1 = mccbender7_h1;
MCNUM w2 = mccbender7_w2;
MCNUM h2 = mccbender7_h2;
MCNUM l = mccbender7_l;
MCNUM R0 = mccbender7_R0;
MCNUM Qc = mccbender7_Qc;
MCNUM alpha = mccbender7_alpha;
MCNUM m = mccbender7_m;
MCNUM W = mccbender7_W;
MCNUM nslit = mccbender7_nslit;
MCNUM d = mccbender7_d;
MCNUM mleft = mccbender7_mleft;
MCNUM mright = mccbender7_mright;
MCNUM mtop = mccbender7_mtop;
MCNUM mbottom = mccbender7_mbottom;
MCNUM nhslit = mccbender7_nhslit;
MCNUM G = mccbender7_G;
MCNUM aleft = mccbender7_aleft;
MCNUM aright = mccbender7_aright;
MCNUM atop = mccbender7_atop;
MCNUM abottom = mccbender7_abottom;
MCNUM wavy = mccbender7_wavy;
MCNUM wavy_z = mccbender7_wavy_z;
MCNUM wavy_tb = mccbender7_wavy_tb;
MCNUM wavy_lr = mccbender7_wavy_lr;
MCNUM chamfers = mccbender7_chamfers;
MCNUM chamfers_z = mccbender7_chamfers_z;
MCNUM chamfers_lr = mccbender7_chamfers_lr;
MCNUM chamfers_tb = mccbender7_chamfers_tb;
MCNUM nelements = mccbender7_nelements;
MCNUM nu = mccbender7_nu;
MCNUM phase = mccbender7_phase;
char* reflect = mccbender7_reflect;
#line 562 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide_gravity.comp"
{
if (GVars.warnings > 100) {
  fprintf(stderr,"%s: warning: neutron has entered guide, but can not exit !\n", GVars.compcurname);
  fprintf(stderr,"%s: warning: This message has been repeated %g times\n", GVars.compcurname, GVars.warnings);
}
}
#line 22624 "./ISIS_SANS2d.c"
}   /* End of bender7=Guide_gravity() SETTING parameter declarations. */
#undef pTable
#undef GVars
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[11]) fprintf(stderr, "Warning: No neutron could reach Component[11] bender7\n");
    if (mcAbsorbProp[11]) fprintf(stderr, "Warning: %g events were removed in Component[11] bender7=Guide_gravity()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[11]);
  /* User FINALLY code for component 'bender8'. */
  SIG_MESSAGE("bender8 (Finally)");
#define mccompcurname  bender8
#define mccompcurtype  Guide_gravity
#define mccompcurindex 12
#define GVars mccbender8_GVars
#define pTable mccbender8_pTable
{   /* Declarations of bender8=Guide_gravity() SETTING parameters. */
MCNUM w1 = mccbender8_w1;
MCNUM h1 = mccbender8_h1;
MCNUM w2 = mccbender8_w2;
MCNUM h2 = mccbender8_h2;
MCNUM l = mccbender8_l;
MCNUM R0 = mccbender8_R0;
MCNUM Qc = mccbender8_Qc;
MCNUM alpha = mccbender8_alpha;
MCNUM m = mccbender8_m;
MCNUM W = mccbender8_W;
MCNUM nslit = mccbender8_nslit;
MCNUM d = mccbender8_d;
MCNUM mleft = mccbender8_mleft;
MCNUM mright = mccbender8_mright;
MCNUM mtop = mccbender8_mtop;
MCNUM mbottom = mccbender8_mbottom;
MCNUM nhslit = mccbender8_nhslit;
MCNUM G = mccbender8_G;
MCNUM aleft = mccbender8_aleft;
MCNUM aright = mccbender8_aright;
MCNUM atop = mccbender8_atop;
MCNUM abottom = mccbender8_abottom;
MCNUM wavy = mccbender8_wavy;
MCNUM wavy_z = mccbender8_wavy_z;
MCNUM wavy_tb = mccbender8_wavy_tb;
MCNUM wavy_lr = mccbender8_wavy_lr;
MCNUM chamfers = mccbender8_chamfers;
MCNUM chamfers_z = mccbender8_chamfers_z;
MCNUM chamfers_lr = mccbender8_chamfers_lr;
MCNUM chamfers_tb = mccbender8_chamfers_tb;
MCNUM nelements = mccbender8_nelements;
MCNUM nu = mccbender8_nu;
MCNUM phase = mccbender8_phase;
char* reflect = mccbender8_reflect;
#line 562 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide_gravity.comp"
{
if (GVars.warnings > 100) {
  fprintf(stderr,"%s: warning: neutron has entered guide, but can not exit !\n", GVars.compcurname);
  fprintf(stderr,"%s: warning: This message has been repeated %g times\n", GVars.compcurname, GVars.warnings);
}
}
#line 22683 "./ISIS_SANS2d.c"
}   /* End of bender8=Guide_gravity() SETTING parameter declarations. */
#undef pTable
#undef GVars
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[12]) fprintf(stderr, "Warning: No neutron could reach Component[12] bender8\n");
    if (mcAbsorbProp[12]) fprintf(stderr, "Warning: %g events were removed in Component[12] bender8=Guide_gravity()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[12]);
  /* User FINALLY code for component 'bender9'. */
  SIG_MESSAGE("bender9 (Finally)");
#define mccompcurname  bender9
#define mccompcurtype  Guide_gravity
#define mccompcurindex 13
#define GVars mccbender9_GVars
#define pTable mccbender9_pTable
{   /* Declarations of bender9=Guide_gravity() SETTING parameters. */
MCNUM w1 = mccbender9_w1;
MCNUM h1 = mccbender9_h1;
MCNUM w2 = mccbender9_w2;
MCNUM h2 = mccbender9_h2;
MCNUM l = mccbender9_l;
MCNUM R0 = mccbender9_R0;
MCNUM Qc = mccbender9_Qc;
MCNUM alpha = mccbender9_alpha;
MCNUM m = mccbender9_m;
MCNUM W = mccbender9_W;
MCNUM nslit = mccbender9_nslit;
MCNUM d = mccbender9_d;
MCNUM mleft = mccbender9_mleft;
MCNUM mright = mccbender9_mright;
MCNUM mtop = mccbender9_mtop;
MCNUM mbottom = mccbender9_mbottom;
MCNUM nhslit = mccbender9_nhslit;
MCNUM G = mccbender9_G;
MCNUM aleft = mccbender9_aleft;
MCNUM aright = mccbender9_aright;
MCNUM atop = mccbender9_atop;
MCNUM abottom = mccbender9_abottom;
MCNUM wavy = mccbender9_wavy;
MCNUM wavy_z = mccbender9_wavy_z;
MCNUM wavy_tb = mccbender9_wavy_tb;
MCNUM wavy_lr = mccbender9_wavy_lr;
MCNUM chamfers = mccbender9_chamfers;
MCNUM chamfers_z = mccbender9_chamfers_z;
MCNUM chamfers_lr = mccbender9_chamfers_lr;
MCNUM chamfers_tb = mccbender9_chamfers_tb;
MCNUM nelements = mccbender9_nelements;
MCNUM nu = mccbender9_nu;
MCNUM phase = mccbender9_phase;
char* reflect = mccbender9_reflect;
#line 562 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide_gravity.comp"
{
if (GVars.warnings > 100) {
  fprintf(stderr,"%s: warning: neutron has entered guide, but can not exit !\n", GVars.compcurname);
  fprintf(stderr,"%s: warning: This message has been repeated %g times\n", GVars.compcurname, GVars.warnings);
}
}
#line 22742 "./ISIS_SANS2d.c"
}   /* End of bender9=Guide_gravity() SETTING parameter declarations. */
#undef pTable
#undef GVars
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[13]) fprintf(stderr, "Warning: No neutron could reach Component[13] bender9\n");
    if (mcAbsorbProp[13]) fprintf(stderr, "Warning: %g events were removed in Component[13] bender9=Guide_gravity()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[13]);
  /* User FINALLY code for component 'bender10'. */
  SIG_MESSAGE("bender10 (Finally)");
#define mccompcurname  bender10
#define mccompcurtype  Guide_gravity
#define mccompcurindex 14
#define GVars mccbender10_GVars
#define pTable mccbender10_pTable
{   /* Declarations of bender10=Guide_gravity() SETTING parameters. */
MCNUM w1 = mccbender10_w1;
MCNUM h1 = mccbender10_h1;
MCNUM w2 = mccbender10_w2;
MCNUM h2 = mccbender10_h2;
MCNUM l = mccbender10_l;
MCNUM R0 = mccbender10_R0;
MCNUM Qc = mccbender10_Qc;
MCNUM alpha = mccbender10_alpha;
MCNUM m = mccbender10_m;
MCNUM W = mccbender10_W;
MCNUM nslit = mccbender10_nslit;
MCNUM d = mccbender10_d;
MCNUM mleft = mccbender10_mleft;
MCNUM mright = mccbender10_mright;
MCNUM mtop = mccbender10_mtop;
MCNUM mbottom = mccbender10_mbottom;
MCNUM nhslit = mccbender10_nhslit;
MCNUM G = mccbender10_G;
MCNUM aleft = mccbender10_aleft;
MCNUM aright = mccbender10_aright;
MCNUM atop = mccbender10_atop;
MCNUM abottom = mccbender10_abottom;
MCNUM wavy = mccbender10_wavy;
MCNUM wavy_z = mccbender10_wavy_z;
MCNUM wavy_tb = mccbender10_wavy_tb;
MCNUM wavy_lr = mccbender10_wavy_lr;
MCNUM chamfers = mccbender10_chamfers;
MCNUM chamfers_z = mccbender10_chamfers_z;
MCNUM chamfers_lr = mccbender10_chamfers_lr;
MCNUM chamfers_tb = mccbender10_chamfers_tb;
MCNUM nelements = mccbender10_nelements;
MCNUM nu = mccbender10_nu;
MCNUM phase = mccbender10_phase;
char* reflect = mccbender10_reflect;
#line 562 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide_gravity.comp"
{
if (GVars.warnings > 100) {
  fprintf(stderr,"%s: warning: neutron has entered guide, but can not exit !\n", GVars.compcurname);
  fprintf(stderr,"%s: warning: This message has been repeated %g times\n", GVars.compcurname, GVars.warnings);
}
}
#line 22801 "./ISIS_SANS2d.c"
}   /* End of bender10=Guide_gravity() SETTING parameter declarations. */
#undef pTable
#undef GVars
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[14]) fprintf(stderr, "Warning: No neutron could reach Component[14] bender10\n");
    if (mcAbsorbProp[14]) fprintf(stderr, "Warning: %g events were removed in Component[14] bender10=Guide_gravity()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[14]);
    if (!mcNCounter[15]) fprintf(stderr, "Warning: No neutron could reach Component[15] lmonb\n");
    if (mcAbsorbProp[15]) fprintf(stderr, "Warning: %g events were removed in Component[15] lmonb=L_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[15]);
  /* User FINALLY code for component 'psd2'. */
  SIG_MESSAGE("psd2 (Finally)");
#define mccompcurname  psd2
#define mccompcurtype  PSD_monitor
#define mccompcurindex 16
#define PSD_N mccpsd2_PSD_N
#define PSD_p mccpsd2_PSD_p
#define PSD_p2 mccpsd2_PSD_p2
{   /* Declarations of psd2=PSD_monitor() SETTING parameters. */
int nx = mccpsd2_nx;
int ny = mccpsd2_ny;
char* filename = mccpsd2_filename;
MCNUM xmin = mccpsd2_xmin;
MCNUM xmax = mccpsd2_xmax;
MCNUM ymin = mccpsd2_ymin;
MCNUM ymax = mccpsd2_ymax;
MCNUM xwidth = mccpsd2_xwidth;
MCNUM yheight = mccpsd2_yheight;
MCNUM restore_neutron = mccpsd2_restore_neutron;
#line 122 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
{
  destroy_darr2d(PSD_N);
  destroy_darr2d(PSD_p);
  destroy_darr2d(PSD_p2);
}
#line 22838 "./ISIS_SANS2d.c"
}   /* End of psd2=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[16]) fprintf(stderr, "Warning: No neutron could reach Component[16] psd2\n");
    if (mcAbsorbProp[16]) fprintf(stderr, "Warning: %g events were removed in Component[16] psd2=PSD_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[16]);
    if (!mcNCounter[17]) fprintf(stderr, "Warning: No neutron could reach Component[17] guide_in\n");
    if (mcAbsorbProp[17]) fprintf(stderr, "Warning: %g events were removed in Component[17] guide_in=Slit()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[17]);
  /* User FINALLY code for component 'guide_straight1'. */
  SIG_MESSAGE("guide_straight1 (Finally)");
#define mccompcurname  guide_straight1
#define mccompcurtype  Guide_gravity
#define mccompcurindex 18
#define GVars mccguide_straight1_GVars
#define pTable mccguide_straight1_pTable
{   /* Declarations of guide_straight1=Guide_gravity() SETTING parameters. */
MCNUM w1 = mccguide_straight1_w1;
MCNUM h1 = mccguide_straight1_h1;
MCNUM w2 = mccguide_straight1_w2;
MCNUM h2 = mccguide_straight1_h2;
MCNUM l = mccguide_straight1_l;
MCNUM R0 = mccguide_straight1_R0;
MCNUM Qc = mccguide_straight1_Qc;
MCNUM alpha = mccguide_straight1_alpha;
MCNUM m = mccguide_straight1_m;
MCNUM W = mccguide_straight1_W;
MCNUM nslit = mccguide_straight1_nslit;
MCNUM d = mccguide_straight1_d;
MCNUM mleft = mccguide_straight1_mleft;
MCNUM mright = mccguide_straight1_mright;
MCNUM mtop = mccguide_straight1_mtop;
MCNUM mbottom = mccguide_straight1_mbottom;
MCNUM nhslit = mccguide_straight1_nhslit;
MCNUM G = mccguide_straight1_G;
MCNUM aleft = mccguide_straight1_aleft;
MCNUM aright = mccguide_straight1_aright;
MCNUM atop = mccguide_straight1_atop;
MCNUM abottom = mccguide_straight1_abottom;
MCNUM wavy = mccguide_straight1_wavy;
MCNUM wavy_z = mccguide_straight1_wavy_z;
MCNUM wavy_tb = mccguide_straight1_wavy_tb;
MCNUM wavy_lr = mccguide_straight1_wavy_lr;
MCNUM chamfers = mccguide_straight1_chamfers;
MCNUM chamfers_z = mccguide_straight1_chamfers_z;
MCNUM chamfers_lr = mccguide_straight1_chamfers_lr;
MCNUM chamfers_tb = mccguide_straight1_chamfers_tb;
MCNUM nelements = mccguide_straight1_nelements;
MCNUM nu = mccguide_straight1_nu;
MCNUM phase = mccguide_straight1_phase;
char* reflect = mccguide_straight1_reflect;
#line 562 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide_gravity.comp"
{
if (GVars.warnings > 100) {
  fprintf(stderr,"%s: warning: neutron has entered guide, but can not exit !\n", GVars.compcurname);
  fprintf(stderr,"%s: warning: This message has been repeated %g times\n", GVars.compcurname, GVars.warnings);
}
}
#line 22900 "./ISIS_SANS2d.c"
}   /* End of guide_straight1=Guide_gravity() SETTING parameter declarations. */
#undef pTable
#undef GVars
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[18]) fprintf(stderr, "Warning: No neutron could reach Component[18] guide_straight1\n");
    if (mcAbsorbProp[18]) fprintf(stderr, "Warning: %g events were removed in Component[18] guide_straight1=Guide_gravity()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[18]);
  /* User FINALLY code for component 'guide_straight2'. */
  SIG_MESSAGE("guide_straight2 (Finally)");
#define mccompcurname  guide_straight2
#define mccompcurtype  Guide_gravity
#define mccompcurindex 19
#define GVars mccguide_straight2_GVars
#define pTable mccguide_straight2_pTable
{   /* Declarations of guide_straight2=Guide_gravity() SETTING parameters. */
MCNUM w1 = mccguide_straight2_w1;
MCNUM h1 = mccguide_straight2_h1;
MCNUM w2 = mccguide_straight2_w2;
MCNUM h2 = mccguide_straight2_h2;
MCNUM l = mccguide_straight2_l;
MCNUM R0 = mccguide_straight2_R0;
MCNUM Qc = mccguide_straight2_Qc;
MCNUM alpha = mccguide_straight2_alpha;
MCNUM m = mccguide_straight2_m;
MCNUM W = mccguide_straight2_W;
MCNUM nslit = mccguide_straight2_nslit;
MCNUM d = mccguide_straight2_d;
MCNUM mleft = mccguide_straight2_mleft;
MCNUM mright = mccguide_straight2_mright;
MCNUM mtop = mccguide_straight2_mtop;
MCNUM mbottom = mccguide_straight2_mbottom;
MCNUM nhslit = mccguide_straight2_nhslit;
MCNUM G = mccguide_straight2_G;
MCNUM aleft = mccguide_straight2_aleft;
MCNUM aright = mccguide_straight2_aright;
MCNUM atop = mccguide_straight2_atop;
MCNUM abottom = mccguide_straight2_abottom;
MCNUM wavy = mccguide_straight2_wavy;
MCNUM wavy_z = mccguide_straight2_wavy_z;
MCNUM wavy_tb = mccguide_straight2_wavy_tb;
MCNUM wavy_lr = mccguide_straight2_wavy_lr;
MCNUM chamfers = mccguide_straight2_chamfers;
MCNUM chamfers_z = mccguide_straight2_chamfers_z;
MCNUM chamfers_lr = mccguide_straight2_chamfers_lr;
MCNUM chamfers_tb = mccguide_straight2_chamfers_tb;
MCNUM nelements = mccguide_straight2_nelements;
MCNUM nu = mccguide_straight2_nu;
MCNUM phase = mccguide_straight2_phase;
char* reflect = mccguide_straight2_reflect;
#line 562 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide_gravity.comp"
{
if (GVars.warnings > 100) {
  fprintf(stderr,"%s: warning: neutron has entered guide, but can not exit !\n", GVars.compcurname);
  fprintf(stderr,"%s: warning: This message has been repeated %g times\n", GVars.compcurname, GVars.warnings);
}
}
#line 22959 "./ISIS_SANS2d.c"
}   /* End of guide_straight2=Guide_gravity() SETTING parameter declarations. */
#undef pTable
#undef GVars
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[19]) fprintf(stderr, "Warning: No neutron could reach Component[19] guide_straight2\n");
    if (mcAbsorbProp[19]) fprintf(stderr, "Warning: %g events were removed in Component[19] guide_straight2=Guide_gravity()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[19]);
  /* User FINALLY code for component 'guide_straight3'. */
  SIG_MESSAGE("guide_straight3 (Finally)");
#define mccompcurname  guide_straight3
#define mccompcurtype  Guide_gravity
#define mccompcurindex 20
#define GVars mccguide_straight3_GVars
#define pTable mccguide_straight3_pTable
{   /* Declarations of guide_straight3=Guide_gravity() SETTING parameters. */
MCNUM w1 = mccguide_straight3_w1;
MCNUM h1 = mccguide_straight3_h1;
MCNUM w2 = mccguide_straight3_w2;
MCNUM h2 = mccguide_straight3_h2;
MCNUM l = mccguide_straight3_l;
MCNUM R0 = mccguide_straight3_R0;
MCNUM Qc = mccguide_straight3_Qc;
MCNUM alpha = mccguide_straight3_alpha;
MCNUM m = mccguide_straight3_m;
MCNUM W = mccguide_straight3_W;
MCNUM nslit = mccguide_straight3_nslit;
MCNUM d = mccguide_straight3_d;
MCNUM mleft = mccguide_straight3_mleft;
MCNUM mright = mccguide_straight3_mright;
MCNUM mtop = mccguide_straight3_mtop;
MCNUM mbottom = mccguide_straight3_mbottom;
MCNUM nhslit = mccguide_straight3_nhslit;
MCNUM G = mccguide_straight3_G;
MCNUM aleft = mccguide_straight3_aleft;
MCNUM aright = mccguide_straight3_aright;
MCNUM atop = mccguide_straight3_atop;
MCNUM abottom = mccguide_straight3_abottom;
MCNUM wavy = mccguide_straight3_wavy;
MCNUM wavy_z = mccguide_straight3_wavy_z;
MCNUM wavy_tb = mccguide_straight3_wavy_tb;
MCNUM wavy_lr = mccguide_straight3_wavy_lr;
MCNUM chamfers = mccguide_straight3_chamfers;
MCNUM chamfers_z = mccguide_straight3_chamfers_z;
MCNUM chamfers_lr = mccguide_straight3_chamfers_lr;
MCNUM chamfers_tb = mccguide_straight3_chamfers_tb;
MCNUM nelements = mccguide_straight3_nelements;
MCNUM nu = mccguide_straight3_nu;
MCNUM phase = mccguide_straight3_phase;
char* reflect = mccguide_straight3_reflect;
#line 562 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide_gravity.comp"
{
if (GVars.warnings > 100) {
  fprintf(stderr,"%s: warning: neutron has entered guide, but can not exit !\n", GVars.compcurname);
  fprintf(stderr,"%s: warning: This message has been repeated %g times\n", GVars.compcurname, GVars.warnings);
}
}
#line 23018 "./ISIS_SANS2d.c"
}   /* End of guide_straight3=Guide_gravity() SETTING parameter declarations. */
#undef pTable
#undef GVars
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[20]) fprintf(stderr, "Warning: No neutron could reach Component[20] guide_straight3\n");
    if (mcAbsorbProp[20]) fprintf(stderr, "Warning: %g events were removed in Component[20] guide_straight3=Guide_gravity()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[20]);
  /* User FINALLY code for component 'guide_straight4'. */
  SIG_MESSAGE("guide_straight4 (Finally)");
#define mccompcurname  guide_straight4
#define mccompcurtype  Guide_gravity
#define mccompcurindex 21
#define GVars mccguide_straight4_GVars
#define pTable mccguide_straight4_pTable
{   /* Declarations of guide_straight4=Guide_gravity() SETTING parameters. */
MCNUM w1 = mccguide_straight4_w1;
MCNUM h1 = mccguide_straight4_h1;
MCNUM w2 = mccguide_straight4_w2;
MCNUM h2 = mccguide_straight4_h2;
MCNUM l = mccguide_straight4_l;
MCNUM R0 = mccguide_straight4_R0;
MCNUM Qc = mccguide_straight4_Qc;
MCNUM alpha = mccguide_straight4_alpha;
MCNUM m = mccguide_straight4_m;
MCNUM W = mccguide_straight4_W;
MCNUM nslit = mccguide_straight4_nslit;
MCNUM d = mccguide_straight4_d;
MCNUM mleft = mccguide_straight4_mleft;
MCNUM mright = mccguide_straight4_mright;
MCNUM mtop = mccguide_straight4_mtop;
MCNUM mbottom = mccguide_straight4_mbottom;
MCNUM nhslit = mccguide_straight4_nhslit;
MCNUM G = mccguide_straight4_G;
MCNUM aleft = mccguide_straight4_aleft;
MCNUM aright = mccguide_straight4_aright;
MCNUM atop = mccguide_straight4_atop;
MCNUM abottom = mccguide_straight4_abottom;
MCNUM wavy = mccguide_straight4_wavy;
MCNUM wavy_z = mccguide_straight4_wavy_z;
MCNUM wavy_tb = mccguide_straight4_wavy_tb;
MCNUM wavy_lr = mccguide_straight4_wavy_lr;
MCNUM chamfers = mccguide_straight4_chamfers;
MCNUM chamfers_z = mccguide_straight4_chamfers_z;
MCNUM chamfers_lr = mccguide_straight4_chamfers_lr;
MCNUM chamfers_tb = mccguide_straight4_chamfers_tb;
MCNUM nelements = mccguide_straight4_nelements;
MCNUM nu = mccguide_straight4_nu;
MCNUM phase = mccguide_straight4_phase;
char* reflect = mccguide_straight4_reflect;
#line 562 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide_gravity.comp"
{
if (GVars.warnings > 100) {
  fprintf(stderr,"%s: warning: neutron has entered guide, but can not exit !\n", GVars.compcurname);
  fprintf(stderr,"%s: warning: This message has been repeated %g times\n", GVars.compcurname, GVars.warnings);
}
}
#line 23077 "./ISIS_SANS2d.c"
}   /* End of guide_straight4=Guide_gravity() SETTING parameter declarations. */
#undef pTable
#undef GVars
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[21]) fprintf(stderr, "Warning: No neutron could reach Component[21] guide_straight4\n");
    if (mcAbsorbProp[21]) fprintf(stderr, "Warning: %g events were removed in Component[21] guide_straight4=Guide_gravity()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[21]);
  /* User FINALLY code for component 'psd3'. */
  SIG_MESSAGE("psd3 (Finally)");
#define mccompcurname  psd3
#define mccompcurtype  PSD_monitor
#define mccompcurindex 22
#define PSD_N mccpsd3_PSD_N
#define PSD_p mccpsd3_PSD_p
#define PSD_p2 mccpsd3_PSD_p2
{   /* Declarations of psd3=PSD_monitor() SETTING parameters. */
int nx = mccpsd3_nx;
int ny = mccpsd3_ny;
char* filename = mccpsd3_filename;
MCNUM xmin = mccpsd3_xmin;
MCNUM xmax = mccpsd3_xmax;
MCNUM ymin = mccpsd3_ymin;
MCNUM ymax = mccpsd3_ymax;
MCNUM xwidth = mccpsd3_xwidth;
MCNUM yheight = mccpsd3_yheight;
MCNUM restore_neutron = mccpsd3_restore_neutron;
#line 122 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
{
  destroy_darr2d(PSD_N);
  destroy_darr2d(PSD_p);
  destroy_darr2d(PSD_p2);
}
#line 23112 "./ISIS_SANS2d.c"
}   /* End of psd3=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[22]) fprintf(stderr, "Warning: No neutron could reach Component[22] psd3\n");
    if (mcAbsorbProp[22]) fprintf(stderr, "Warning: %g events were removed in Component[22] psd3=PSD_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[22]);
    if (!mcNCounter[23]) fprintf(stderr, "Warning: No neutron could reach Component[23] aperture1\n");
    if (mcAbsorbProp[23]) fprintf(stderr, "Warning: %g events were removed in Component[23] aperture1=Slit()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[23]);
    if (!mcNCounter[24]) fprintf(stderr, "Warning: No neutron could reach Component[24] lmonitor2\n");
    if (mcAbsorbProp[24]) fprintf(stderr, "Warning: %g events were removed in Component[24] lmonitor2=L_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[24]);
    if (!mcNCounter[25]) fprintf(stderr, "Warning: No neutron could reach Component[25] S6\n");
    if (mcAbsorbProp[25]) fprintf(stderr, "Warning: %g events were removed in Component[25] S6=Slit()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[25]);
    if (!mcNCounter[26]) fprintf(stderr, "Warning: No neutron could reach Component[26] APERTURE2\n");
    if (mcAbsorbProp[26]) fprintf(stderr, "Warning: %g events were removed in Component[26] APERTURE2=Slit()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[26]);
    if (!mcNCounter[27]) fprintf(stderr, "Warning: No neutron could reach Component[27] lmon2\n");
    if (mcAbsorbProp[27]) fprintf(stderr, "Warning: %g events were removed in Component[27] lmon2=L_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[27]);
  /* User FINALLY code for component 'psd4'. */
  SIG_MESSAGE("psd4 (Finally)");
#define mccompcurname  psd4
#define mccompcurtype  PSD_monitor
#define mccompcurindex 28
#define PSD_N mccpsd4_PSD_N
#define PSD_p mccpsd4_PSD_p
#define PSD_p2 mccpsd4_PSD_p2
{   /* Declarations of psd4=PSD_monitor() SETTING parameters. */
int nx = mccpsd4_nx;
int ny = mccpsd4_ny;
char* filename = mccpsd4_filename;
MCNUM xmin = mccpsd4_xmin;
MCNUM xmax = mccpsd4_xmax;
MCNUM ymin = mccpsd4_ymin;
MCNUM ymax = mccpsd4_ymax;
MCNUM xwidth = mccpsd4_xwidth;
MCNUM yheight = mccpsd4_yheight;
MCNUM restore_neutron = mccpsd4_restore_neutron;
#line 122 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
{
  destroy_darr2d(PSD_N);
  destroy_darr2d(PSD_p);
  destroy_darr2d(PSD_p2);
}
#line 23158 "./ISIS_SANS2d.c"
}   /* End of psd4=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[28]) fprintf(stderr, "Warning: No neutron could reach Component[28] psd4\n");
    if (mcAbsorbProp[28]) fprintf(stderr, "Warning: %g events were removed in Component[28] psd4=PSD_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[28]);
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
#define mccompcurtype  Arm
#define mccompcurindex 1
#line 40 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Arm.comp"
{
  /* A bit ugly; hard-coded dimensions. */
  
  line(0,0,0,0.2,0,0);
  line(0,0,0,0,0.2,0);
  line(0,0,0,0,0,0.2);
}
#line 23198 "./ISIS_SANS2d.c"
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
#line 1119 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../contrib/ISIS_moderator.comp"
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
#line 23264 "./ISIS_SANS2d.c"
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

  /* MCDISPLAY code for component 'lmon1'. */
  SIG_MESSAGE("lmon1 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "lmon1");
#define mccompcurname  lmon1
#define mccompcurtype  L_monitor
#define mccompcurindex 3
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
#line 120 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/L_monitor.comp"
{
  
  multiline(5, (double)xmin, (double)ymin, 0.0,
               (double)xmax, (double)ymin, 0.0,
               (double)xmax, (double)ymax, 0.0,
               (double)xmin, (double)ymax, 0.0,
               (double)xmin, (double)ymin, 0.0);
}
#line 23315 "./ISIS_SANS2d.c"
}   /* End of lmon1=L_monitor() SETTING parameter declarations. */
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
#define mccompcurindex 4
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
#line 23354 "./ISIS_SANS2d.c"
}   /* End of psd1=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'bender1'. */
  SIG_MESSAGE("bender1 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "bender1");
#define mccompcurname  bender1
#define mccompcurtype  Guide_gravity
#define mccompcurindex 5
#define GVars mccbender1_GVars
#define pTable mccbender1_pTable
{   /* Declarations of bender1=Guide_gravity() SETTING parameters. */
MCNUM w1 = mccbender1_w1;
MCNUM h1 = mccbender1_h1;
MCNUM w2 = mccbender1_w2;
MCNUM h2 = mccbender1_h2;
MCNUM l = mccbender1_l;
MCNUM R0 = mccbender1_R0;
MCNUM Qc = mccbender1_Qc;
MCNUM alpha = mccbender1_alpha;
MCNUM m = mccbender1_m;
MCNUM W = mccbender1_W;
MCNUM nslit = mccbender1_nslit;
MCNUM d = mccbender1_d;
MCNUM mleft = mccbender1_mleft;
MCNUM mright = mccbender1_mright;
MCNUM mtop = mccbender1_mtop;
MCNUM mbottom = mccbender1_mbottom;
MCNUM nhslit = mccbender1_nhslit;
MCNUM G = mccbender1_G;
MCNUM aleft = mccbender1_aleft;
MCNUM aright = mccbender1_aright;
MCNUM atop = mccbender1_atop;
MCNUM abottom = mccbender1_abottom;
MCNUM wavy = mccbender1_wavy;
MCNUM wavy_z = mccbender1_wavy_z;
MCNUM wavy_tb = mccbender1_wavy_tb;
MCNUM wavy_lr = mccbender1_wavy_lr;
MCNUM chamfers = mccbender1_chamfers;
MCNUM chamfers_z = mccbender1_chamfers_z;
MCNUM chamfers_lr = mccbender1_chamfers_lr;
MCNUM chamfers_tb = mccbender1_chamfers_tb;
MCNUM nelements = mccbender1_nelements;
MCNUM nu = mccbender1_nu;
MCNUM phase = mccbender1_phase;
char* reflect = mccbender1_reflect;
#line 571 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide_gravity.comp"
{

  if (l > 0 && nelements > 0) {
    int i,j,n;
    double x1,x2,x3,x4;
    double y1,y2,y3,y4;
    double nel = (nelements > 11 ? 11 : nelements);

    
    for (n=0; n<nel; n++)
    {
      double z0, z1;
      z0 =     n*(l/nel);
      z1 = (n+1)*(l/nel);

      for(j = 0; j < nhslit; j++)
      {
        y1 = j*(GVars.h1c+d)         - h1/2.0;
        y2 = j*(GVars.h2c+d)         - h2/2.0;
        y3 = (j+1)*(GVars.h1c+d) - d - h1/2.0;
        y4 = (j+1)*(GVars.h2c+d) - d - h2/2.0;
        for(i = 0; i < nslit; i++)
        {
          x1 = i*(GVars.w1c+d)         - w1/2.0;
          x2 = i*(GVars.w2c+d)         - w2/2.0;
          x3 = (i+1)*(GVars.w1c+d) - d - w1/2.0;
          x4 = (i+1)*(GVars.w2c+d) - d - w2/2.0;
          multiline(5,
                    x1, y1, z0,
                    x2, y2, z1,
                    x2, y4, z1,
                    x1, y3, z0,
                    x1, y1, z0);
          multiline(5,
                    x3, y1, z0,
                    x4, y2, z1,
                    x4, y4, z1,
                    x3, y3, z0,
                    x3, y1, z0);
        }
        line(-w1/2.0, y1, z0, w1/2.0, y1, z0);
        line(-w2/2.0, y2, z1, w2/2.0, y2, z1);
      }
    }

    if (nu || phase) {
      double radius = sqrt(w1*w1+l*l);
      /* cylinder top/center/bottom  */
      circle("xz", 0,-h1/2,l/2,radius);
      circle("xz", 0,0    ,l/2,radius);
      circle("xz", 0, h1/2,l/2,radius);
    }
  }
  else {
    /* A bit ugly; hard-coded dimensions. */
    
    line(0,0,0,0.2,0,0);
    line(0,0,0,0,0.2,0);
    line(0,0,0,0,0,0.2);
  }

}
#line 23469 "./ISIS_SANS2d.c"
}   /* End of bender1=Guide_gravity() SETTING parameter declarations. */
#undef pTable
#undef GVars
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'bender2'. */
  SIG_MESSAGE("bender2 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "bender2");
#define mccompcurname  bender2
#define mccompcurtype  Guide_gravity
#define mccompcurindex 6
#define GVars mccbender2_GVars
#define pTable mccbender2_pTable
{   /* Declarations of bender2=Guide_gravity() SETTING parameters. */
MCNUM w1 = mccbender2_w1;
MCNUM h1 = mccbender2_h1;
MCNUM w2 = mccbender2_w2;
MCNUM h2 = mccbender2_h2;
MCNUM l = mccbender2_l;
MCNUM R0 = mccbender2_R0;
MCNUM Qc = mccbender2_Qc;
MCNUM alpha = mccbender2_alpha;
MCNUM m = mccbender2_m;
MCNUM W = mccbender2_W;
MCNUM nslit = mccbender2_nslit;
MCNUM d = mccbender2_d;
MCNUM mleft = mccbender2_mleft;
MCNUM mright = mccbender2_mright;
MCNUM mtop = mccbender2_mtop;
MCNUM mbottom = mccbender2_mbottom;
MCNUM nhslit = mccbender2_nhslit;
MCNUM G = mccbender2_G;
MCNUM aleft = mccbender2_aleft;
MCNUM aright = mccbender2_aright;
MCNUM atop = mccbender2_atop;
MCNUM abottom = mccbender2_abottom;
MCNUM wavy = mccbender2_wavy;
MCNUM wavy_z = mccbender2_wavy_z;
MCNUM wavy_tb = mccbender2_wavy_tb;
MCNUM wavy_lr = mccbender2_wavy_lr;
MCNUM chamfers = mccbender2_chamfers;
MCNUM chamfers_z = mccbender2_chamfers_z;
MCNUM chamfers_lr = mccbender2_chamfers_lr;
MCNUM chamfers_tb = mccbender2_chamfers_tb;
MCNUM nelements = mccbender2_nelements;
MCNUM nu = mccbender2_nu;
MCNUM phase = mccbender2_phase;
char* reflect = mccbender2_reflect;
#line 571 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide_gravity.comp"
{

  if (l > 0 && nelements > 0) {
    int i,j,n;
    double x1,x2,x3,x4;
    double y1,y2,y3,y4;
    double nel = (nelements > 11 ? 11 : nelements);

    
    for (n=0; n<nel; n++)
    {
      double z0, z1;
      z0 =     n*(l/nel);
      z1 = (n+1)*(l/nel);

      for(j = 0; j < nhslit; j++)
      {
        y1 = j*(GVars.h1c+d)         - h1/2.0;
        y2 = j*(GVars.h2c+d)         - h2/2.0;
        y3 = (j+1)*(GVars.h1c+d) - d - h1/2.0;
        y4 = (j+1)*(GVars.h2c+d) - d - h2/2.0;
        for(i = 0; i < nslit; i++)
        {
          x1 = i*(GVars.w1c+d)         - w1/2.0;
          x2 = i*(GVars.w2c+d)         - w2/2.0;
          x3 = (i+1)*(GVars.w1c+d) - d - w1/2.0;
          x4 = (i+1)*(GVars.w2c+d) - d - w2/2.0;
          multiline(5,
                    x1, y1, z0,
                    x2, y2, z1,
                    x2, y4, z1,
                    x1, y3, z0,
                    x1, y1, z0);
          multiline(5,
                    x3, y1, z0,
                    x4, y2, z1,
                    x4, y4, z1,
                    x3, y3, z0,
                    x3, y1, z0);
        }
        line(-w1/2.0, y1, z0, w1/2.0, y1, z0);
        line(-w2/2.0, y2, z1, w2/2.0, y2, z1);
      }
    }

    if (nu || phase) {
      double radius = sqrt(w1*w1+l*l);
      /* cylinder top/center/bottom  */
      circle("xz", 0,-h1/2,l/2,radius);
      circle("xz", 0,0    ,l/2,radius);
      circle("xz", 0, h1/2,l/2,radius);
    }
  }
  else {
    /* A bit ugly; hard-coded dimensions. */
    
    line(0,0,0,0.2,0,0);
    line(0,0,0,0,0.2,0);
    line(0,0,0,0,0,0.2);
  }

}
#line 23583 "./ISIS_SANS2d.c"
}   /* End of bender2=Guide_gravity() SETTING parameter declarations. */
#undef pTable
#undef GVars
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'bender3'. */
  SIG_MESSAGE("bender3 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "bender3");
#define mccompcurname  bender3
#define mccompcurtype  Guide_gravity
#define mccompcurindex 7
#define GVars mccbender3_GVars
#define pTable mccbender3_pTable
{   /* Declarations of bender3=Guide_gravity() SETTING parameters. */
MCNUM w1 = mccbender3_w1;
MCNUM h1 = mccbender3_h1;
MCNUM w2 = mccbender3_w2;
MCNUM h2 = mccbender3_h2;
MCNUM l = mccbender3_l;
MCNUM R0 = mccbender3_R0;
MCNUM Qc = mccbender3_Qc;
MCNUM alpha = mccbender3_alpha;
MCNUM m = mccbender3_m;
MCNUM W = mccbender3_W;
MCNUM nslit = mccbender3_nslit;
MCNUM d = mccbender3_d;
MCNUM mleft = mccbender3_mleft;
MCNUM mright = mccbender3_mright;
MCNUM mtop = mccbender3_mtop;
MCNUM mbottom = mccbender3_mbottom;
MCNUM nhslit = mccbender3_nhslit;
MCNUM G = mccbender3_G;
MCNUM aleft = mccbender3_aleft;
MCNUM aright = mccbender3_aright;
MCNUM atop = mccbender3_atop;
MCNUM abottom = mccbender3_abottom;
MCNUM wavy = mccbender3_wavy;
MCNUM wavy_z = mccbender3_wavy_z;
MCNUM wavy_tb = mccbender3_wavy_tb;
MCNUM wavy_lr = mccbender3_wavy_lr;
MCNUM chamfers = mccbender3_chamfers;
MCNUM chamfers_z = mccbender3_chamfers_z;
MCNUM chamfers_lr = mccbender3_chamfers_lr;
MCNUM chamfers_tb = mccbender3_chamfers_tb;
MCNUM nelements = mccbender3_nelements;
MCNUM nu = mccbender3_nu;
MCNUM phase = mccbender3_phase;
char* reflect = mccbender3_reflect;
#line 571 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide_gravity.comp"
{

  if (l > 0 && nelements > 0) {
    int i,j,n;
    double x1,x2,x3,x4;
    double y1,y2,y3,y4;
    double nel = (nelements > 11 ? 11 : nelements);

    
    for (n=0; n<nel; n++)
    {
      double z0, z1;
      z0 =     n*(l/nel);
      z1 = (n+1)*(l/nel);

      for(j = 0; j < nhslit; j++)
      {
        y1 = j*(GVars.h1c+d)         - h1/2.0;
        y2 = j*(GVars.h2c+d)         - h2/2.0;
        y3 = (j+1)*(GVars.h1c+d) - d - h1/2.0;
        y4 = (j+1)*(GVars.h2c+d) - d - h2/2.0;
        for(i = 0; i < nslit; i++)
        {
          x1 = i*(GVars.w1c+d)         - w1/2.0;
          x2 = i*(GVars.w2c+d)         - w2/2.0;
          x3 = (i+1)*(GVars.w1c+d) - d - w1/2.0;
          x4 = (i+1)*(GVars.w2c+d) - d - w2/2.0;
          multiline(5,
                    x1, y1, z0,
                    x2, y2, z1,
                    x2, y4, z1,
                    x1, y3, z0,
                    x1, y1, z0);
          multiline(5,
                    x3, y1, z0,
                    x4, y2, z1,
                    x4, y4, z1,
                    x3, y3, z0,
                    x3, y1, z0);
        }
        line(-w1/2.0, y1, z0, w1/2.0, y1, z0);
        line(-w2/2.0, y2, z1, w2/2.0, y2, z1);
      }
    }

    if (nu || phase) {
      double radius = sqrt(w1*w1+l*l);
      /* cylinder top/center/bottom  */
      circle("xz", 0,-h1/2,l/2,radius);
      circle("xz", 0,0    ,l/2,radius);
      circle("xz", 0, h1/2,l/2,radius);
    }
  }
  else {
    /* A bit ugly; hard-coded dimensions. */
    
    line(0,0,0,0.2,0,0);
    line(0,0,0,0,0.2,0);
    line(0,0,0,0,0,0.2);
  }

}
#line 23697 "./ISIS_SANS2d.c"
}   /* End of bender3=Guide_gravity() SETTING parameter declarations. */
#undef pTable
#undef GVars
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'bender4'. */
  SIG_MESSAGE("bender4 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "bender4");
#define mccompcurname  bender4
#define mccompcurtype  Guide_gravity
#define mccompcurindex 8
#define GVars mccbender4_GVars
#define pTable mccbender4_pTable
{   /* Declarations of bender4=Guide_gravity() SETTING parameters. */
MCNUM w1 = mccbender4_w1;
MCNUM h1 = mccbender4_h1;
MCNUM w2 = mccbender4_w2;
MCNUM h2 = mccbender4_h2;
MCNUM l = mccbender4_l;
MCNUM R0 = mccbender4_R0;
MCNUM Qc = mccbender4_Qc;
MCNUM alpha = mccbender4_alpha;
MCNUM m = mccbender4_m;
MCNUM W = mccbender4_W;
MCNUM nslit = mccbender4_nslit;
MCNUM d = mccbender4_d;
MCNUM mleft = mccbender4_mleft;
MCNUM mright = mccbender4_mright;
MCNUM mtop = mccbender4_mtop;
MCNUM mbottom = mccbender4_mbottom;
MCNUM nhslit = mccbender4_nhslit;
MCNUM G = mccbender4_G;
MCNUM aleft = mccbender4_aleft;
MCNUM aright = mccbender4_aright;
MCNUM atop = mccbender4_atop;
MCNUM abottom = mccbender4_abottom;
MCNUM wavy = mccbender4_wavy;
MCNUM wavy_z = mccbender4_wavy_z;
MCNUM wavy_tb = mccbender4_wavy_tb;
MCNUM wavy_lr = mccbender4_wavy_lr;
MCNUM chamfers = mccbender4_chamfers;
MCNUM chamfers_z = mccbender4_chamfers_z;
MCNUM chamfers_lr = mccbender4_chamfers_lr;
MCNUM chamfers_tb = mccbender4_chamfers_tb;
MCNUM nelements = mccbender4_nelements;
MCNUM nu = mccbender4_nu;
MCNUM phase = mccbender4_phase;
char* reflect = mccbender4_reflect;
#line 571 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide_gravity.comp"
{

  if (l > 0 && nelements > 0) {
    int i,j,n;
    double x1,x2,x3,x4;
    double y1,y2,y3,y4;
    double nel = (nelements > 11 ? 11 : nelements);

    
    for (n=0; n<nel; n++)
    {
      double z0, z1;
      z0 =     n*(l/nel);
      z1 = (n+1)*(l/nel);

      for(j = 0; j < nhslit; j++)
      {
        y1 = j*(GVars.h1c+d)         - h1/2.0;
        y2 = j*(GVars.h2c+d)         - h2/2.0;
        y3 = (j+1)*(GVars.h1c+d) - d - h1/2.0;
        y4 = (j+1)*(GVars.h2c+d) - d - h2/2.0;
        for(i = 0; i < nslit; i++)
        {
          x1 = i*(GVars.w1c+d)         - w1/2.0;
          x2 = i*(GVars.w2c+d)         - w2/2.0;
          x3 = (i+1)*(GVars.w1c+d) - d - w1/2.0;
          x4 = (i+1)*(GVars.w2c+d) - d - w2/2.0;
          multiline(5,
                    x1, y1, z0,
                    x2, y2, z1,
                    x2, y4, z1,
                    x1, y3, z0,
                    x1, y1, z0);
          multiline(5,
                    x3, y1, z0,
                    x4, y2, z1,
                    x4, y4, z1,
                    x3, y3, z0,
                    x3, y1, z0);
        }
        line(-w1/2.0, y1, z0, w1/2.0, y1, z0);
        line(-w2/2.0, y2, z1, w2/2.0, y2, z1);
      }
    }

    if (nu || phase) {
      double radius = sqrt(w1*w1+l*l);
      /* cylinder top/center/bottom  */
      circle("xz", 0,-h1/2,l/2,radius);
      circle("xz", 0,0    ,l/2,radius);
      circle("xz", 0, h1/2,l/2,radius);
    }
  }
  else {
    /* A bit ugly; hard-coded dimensions. */
    
    line(0,0,0,0.2,0,0);
    line(0,0,0,0,0.2,0);
    line(0,0,0,0,0,0.2);
  }

}
#line 23811 "./ISIS_SANS2d.c"
}   /* End of bender4=Guide_gravity() SETTING parameter declarations. */
#undef pTable
#undef GVars
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'bender5'. */
  SIG_MESSAGE("bender5 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "bender5");
#define mccompcurname  bender5
#define mccompcurtype  Guide_gravity
#define mccompcurindex 9
#define GVars mccbender5_GVars
#define pTable mccbender5_pTable
{   /* Declarations of bender5=Guide_gravity() SETTING parameters. */
MCNUM w1 = mccbender5_w1;
MCNUM h1 = mccbender5_h1;
MCNUM w2 = mccbender5_w2;
MCNUM h2 = mccbender5_h2;
MCNUM l = mccbender5_l;
MCNUM R0 = mccbender5_R0;
MCNUM Qc = mccbender5_Qc;
MCNUM alpha = mccbender5_alpha;
MCNUM m = mccbender5_m;
MCNUM W = mccbender5_W;
MCNUM nslit = mccbender5_nslit;
MCNUM d = mccbender5_d;
MCNUM mleft = mccbender5_mleft;
MCNUM mright = mccbender5_mright;
MCNUM mtop = mccbender5_mtop;
MCNUM mbottom = mccbender5_mbottom;
MCNUM nhslit = mccbender5_nhslit;
MCNUM G = mccbender5_G;
MCNUM aleft = mccbender5_aleft;
MCNUM aright = mccbender5_aright;
MCNUM atop = mccbender5_atop;
MCNUM abottom = mccbender5_abottom;
MCNUM wavy = mccbender5_wavy;
MCNUM wavy_z = mccbender5_wavy_z;
MCNUM wavy_tb = mccbender5_wavy_tb;
MCNUM wavy_lr = mccbender5_wavy_lr;
MCNUM chamfers = mccbender5_chamfers;
MCNUM chamfers_z = mccbender5_chamfers_z;
MCNUM chamfers_lr = mccbender5_chamfers_lr;
MCNUM chamfers_tb = mccbender5_chamfers_tb;
MCNUM nelements = mccbender5_nelements;
MCNUM nu = mccbender5_nu;
MCNUM phase = mccbender5_phase;
char* reflect = mccbender5_reflect;
#line 571 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide_gravity.comp"
{

  if (l > 0 && nelements > 0) {
    int i,j,n;
    double x1,x2,x3,x4;
    double y1,y2,y3,y4;
    double nel = (nelements > 11 ? 11 : nelements);

    
    for (n=0; n<nel; n++)
    {
      double z0, z1;
      z0 =     n*(l/nel);
      z1 = (n+1)*(l/nel);

      for(j = 0; j < nhslit; j++)
      {
        y1 = j*(GVars.h1c+d)         - h1/2.0;
        y2 = j*(GVars.h2c+d)         - h2/2.0;
        y3 = (j+1)*(GVars.h1c+d) - d - h1/2.0;
        y4 = (j+1)*(GVars.h2c+d) - d - h2/2.0;
        for(i = 0; i < nslit; i++)
        {
          x1 = i*(GVars.w1c+d)         - w1/2.0;
          x2 = i*(GVars.w2c+d)         - w2/2.0;
          x3 = (i+1)*(GVars.w1c+d) - d - w1/2.0;
          x4 = (i+1)*(GVars.w2c+d) - d - w2/2.0;
          multiline(5,
                    x1, y1, z0,
                    x2, y2, z1,
                    x2, y4, z1,
                    x1, y3, z0,
                    x1, y1, z0);
          multiline(5,
                    x3, y1, z0,
                    x4, y2, z1,
                    x4, y4, z1,
                    x3, y3, z0,
                    x3, y1, z0);
        }
        line(-w1/2.0, y1, z0, w1/2.0, y1, z0);
        line(-w2/2.0, y2, z1, w2/2.0, y2, z1);
      }
    }

    if (nu || phase) {
      double radius = sqrt(w1*w1+l*l);
      /* cylinder top/center/bottom  */
      circle("xz", 0,-h1/2,l/2,radius);
      circle("xz", 0,0    ,l/2,radius);
      circle("xz", 0, h1/2,l/2,radius);
    }
  }
  else {
    /* A bit ugly; hard-coded dimensions. */
    
    line(0,0,0,0.2,0,0);
    line(0,0,0,0,0.2,0);
    line(0,0,0,0,0,0.2);
  }

}
#line 23925 "./ISIS_SANS2d.c"
}   /* End of bender5=Guide_gravity() SETTING parameter declarations. */
#undef pTable
#undef GVars
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'bender6'. */
  SIG_MESSAGE("bender6 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "bender6");
#define mccompcurname  bender6
#define mccompcurtype  Guide_gravity
#define mccompcurindex 10
#define GVars mccbender6_GVars
#define pTable mccbender6_pTable
{   /* Declarations of bender6=Guide_gravity() SETTING parameters. */
MCNUM w1 = mccbender6_w1;
MCNUM h1 = mccbender6_h1;
MCNUM w2 = mccbender6_w2;
MCNUM h2 = mccbender6_h2;
MCNUM l = mccbender6_l;
MCNUM R0 = mccbender6_R0;
MCNUM Qc = mccbender6_Qc;
MCNUM alpha = mccbender6_alpha;
MCNUM m = mccbender6_m;
MCNUM W = mccbender6_W;
MCNUM nslit = mccbender6_nslit;
MCNUM d = mccbender6_d;
MCNUM mleft = mccbender6_mleft;
MCNUM mright = mccbender6_mright;
MCNUM mtop = mccbender6_mtop;
MCNUM mbottom = mccbender6_mbottom;
MCNUM nhslit = mccbender6_nhslit;
MCNUM G = mccbender6_G;
MCNUM aleft = mccbender6_aleft;
MCNUM aright = mccbender6_aright;
MCNUM atop = mccbender6_atop;
MCNUM abottom = mccbender6_abottom;
MCNUM wavy = mccbender6_wavy;
MCNUM wavy_z = mccbender6_wavy_z;
MCNUM wavy_tb = mccbender6_wavy_tb;
MCNUM wavy_lr = mccbender6_wavy_lr;
MCNUM chamfers = mccbender6_chamfers;
MCNUM chamfers_z = mccbender6_chamfers_z;
MCNUM chamfers_lr = mccbender6_chamfers_lr;
MCNUM chamfers_tb = mccbender6_chamfers_tb;
MCNUM nelements = mccbender6_nelements;
MCNUM nu = mccbender6_nu;
MCNUM phase = mccbender6_phase;
char* reflect = mccbender6_reflect;
#line 571 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide_gravity.comp"
{

  if (l > 0 && nelements > 0) {
    int i,j,n;
    double x1,x2,x3,x4;
    double y1,y2,y3,y4;
    double nel = (nelements > 11 ? 11 : nelements);

    
    for (n=0; n<nel; n++)
    {
      double z0, z1;
      z0 =     n*(l/nel);
      z1 = (n+1)*(l/nel);

      for(j = 0; j < nhslit; j++)
      {
        y1 = j*(GVars.h1c+d)         - h1/2.0;
        y2 = j*(GVars.h2c+d)         - h2/2.0;
        y3 = (j+1)*(GVars.h1c+d) - d - h1/2.0;
        y4 = (j+1)*(GVars.h2c+d) - d - h2/2.0;
        for(i = 0; i < nslit; i++)
        {
          x1 = i*(GVars.w1c+d)         - w1/2.0;
          x2 = i*(GVars.w2c+d)         - w2/2.0;
          x3 = (i+1)*(GVars.w1c+d) - d - w1/2.0;
          x4 = (i+1)*(GVars.w2c+d) - d - w2/2.0;
          multiline(5,
                    x1, y1, z0,
                    x2, y2, z1,
                    x2, y4, z1,
                    x1, y3, z0,
                    x1, y1, z0);
          multiline(5,
                    x3, y1, z0,
                    x4, y2, z1,
                    x4, y4, z1,
                    x3, y3, z0,
                    x3, y1, z0);
        }
        line(-w1/2.0, y1, z0, w1/2.0, y1, z0);
        line(-w2/2.0, y2, z1, w2/2.0, y2, z1);
      }
    }

    if (nu || phase) {
      double radius = sqrt(w1*w1+l*l);
      /* cylinder top/center/bottom  */
      circle("xz", 0,-h1/2,l/2,radius);
      circle("xz", 0,0    ,l/2,radius);
      circle("xz", 0, h1/2,l/2,radius);
    }
  }
  else {
    /* A bit ugly; hard-coded dimensions. */
    
    line(0,0,0,0.2,0,0);
    line(0,0,0,0,0.2,0);
    line(0,0,0,0,0,0.2);
  }

}
#line 24039 "./ISIS_SANS2d.c"
}   /* End of bender6=Guide_gravity() SETTING parameter declarations. */
#undef pTable
#undef GVars
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'bender7'. */
  SIG_MESSAGE("bender7 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "bender7");
#define mccompcurname  bender7
#define mccompcurtype  Guide_gravity
#define mccompcurindex 11
#define GVars mccbender7_GVars
#define pTable mccbender7_pTable
{   /* Declarations of bender7=Guide_gravity() SETTING parameters. */
MCNUM w1 = mccbender7_w1;
MCNUM h1 = mccbender7_h1;
MCNUM w2 = mccbender7_w2;
MCNUM h2 = mccbender7_h2;
MCNUM l = mccbender7_l;
MCNUM R0 = mccbender7_R0;
MCNUM Qc = mccbender7_Qc;
MCNUM alpha = mccbender7_alpha;
MCNUM m = mccbender7_m;
MCNUM W = mccbender7_W;
MCNUM nslit = mccbender7_nslit;
MCNUM d = mccbender7_d;
MCNUM mleft = mccbender7_mleft;
MCNUM mright = mccbender7_mright;
MCNUM mtop = mccbender7_mtop;
MCNUM mbottom = mccbender7_mbottom;
MCNUM nhslit = mccbender7_nhslit;
MCNUM G = mccbender7_G;
MCNUM aleft = mccbender7_aleft;
MCNUM aright = mccbender7_aright;
MCNUM atop = mccbender7_atop;
MCNUM abottom = mccbender7_abottom;
MCNUM wavy = mccbender7_wavy;
MCNUM wavy_z = mccbender7_wavy_z;
MCNUM wavy_tb = mccbender7_wavy_tb;
MCNUM wavy_lr = mccbender7_wavy_lr;
MCNUM chamfers = mccbender7_chamfers;
MCNUM chamfers_z = mccbender7_chamfers_z;
MCNUM chamfers_lr = mccbender7_chamfers_lr;
MCNUM chamfers_tb = mccbender7_chamfers_tb;
MCNUM nelements = mccbender7_nelements;
MCNUM nu = mccbender7_nu;
MCNUM phase = mccbender7_phase;
char* reflect = mccbender7_reflect;
#line 571 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide_gravity.comp"
{

  if (l > 0 && nelements > 0) {
    int i,j,n;
    double x1,x2,x3,x4;
    double y1,y2,y3,y4;
    double nel = (nelements > 11 ? 11 : nelements);

    
    for (n=0; n<nel; n++)
    {
      double z0, z1;
      z0 =     n*(l/nel);
      z1 = (n+1)*(l/nel);

      for(j = 0; j < nhslit; j++)
      {
        y1 = j*(GVars.h1c+d)         - h1/2.0;
        y2 = j*(GVars.h2c+d)         - h2/2.0;
        y3 = (j+1)*(GVars.h1c+d) - d - h1/2.0;
        y4 = (j+1)*(GVars.h2c+d) - d - h2/2.0;
        for(i = 0; i < nslit; i++)
        {
          x1 = i*(GVars.w1c+d)         - w1/2.0;
          x2 = i*(GVars.w2c+d)         - w2/2.0;
          x3 = (i+1)*(GVars.w1c+d) - d - w1/2.0;
          x4 = (i+1)*(GVars.w2c+d) - d - w2/2.0;
          multiline(5,
                    x1, y1, z0,
                    x2, y2, z1,
                    x2, y4, z1,
                    x1, y3, z0,
                    x1, y1, z0);
          multiline(5,
                    x3, y1, z0,
                    x4, y2, z1,
                    x4, y4, z1,
                    x3, y3, z0,
                    x3, y1, z0);
        }
        line(-w1/2.0, y1, z0, w1/2.0, y1, z0);
        line(-w2/2.0, y2, z1, w2/2.0, y2, z1);
      }
    }

    if (nu || phase) {
      double radius = sqrt(w1*w1+l*l);
      /* cylinder top/center/bottom  */
      circle("xz", 0,-h1/2,l/2,radius);
      circle("xz", 0,0    ,l/2,radius);
      circle("xz", 0, h1/2,l/2,radius);
    }
  }
  else {
    /* A bit ugly; hard-coded dimensions. */
    
    line(0,0,0,0.2,0,0);
    line(0,0,0,0,0.2,0);
    line(0,0,0,0,0,0.2);
  }

}
#line 24153 "./ISIS_SANS2d.c"
}   /* End of bender7=Guide_gravity() SETTING parameter declarations. */
#undef pTable
#undef GVars
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'bender8'. */
  SIG_MESSAGE("bender8 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "bender8");
#define mccompcurname  bender8
#define mccompcurtype  Guide_gravity
#define mccompcurindex 12
#define GVars mccbender8_GVars
#define pTable mccbender8_pTable
{   /* Declarations of bender8=Guide_gravity() SETTING parameters. */
MCNUM w1 = mccbender8_w1;
MCNUM h1 = mccbender8_h1;
MCNUM w2 = mccbender8_w2;
MCNUM h2 = mccbender8_h2;
MCNUM l = mccbender8_l;
MCNUM R0 = mccbender8_R0;
MCNUM Qc = mccbender8_Qc;
MCNUM alpha = mccbender8_alpha;
MCNUM m = mccbender8_m;
MCNUM W = mccbender8_W;
MCNUM nslit = mccbender8_nslit;
MCNUM d = mccbender8_d;
MCNUM mleft = mccbender8_mleft;
MCNUM mright = mccbender8_mright;
MCNUM mtop = mccbender8_mtop;
MCNUM mbottom = mccbender8_mbottom;
MCNUM nhslit = mccbender8_nhslit;
MCNUM G = mccbender8_G;
MCNUM aleft = mccbender8_aleft;
MCNUM aright = mccbender8_aright;
MCNUM atop = mccbender8_atop;
MCNUM abottom = mccbender8_abottom;
MCNUM wavy = mccbender8_wavy;
MCNUM wavy_z = mccbender8_wavy_z;
MCNUM wavy_tb = mccbender8_wavy_tb;
MCNUM wavy_lr = mccbender8_wavy_lr;
MCNUM chamfers = mccbender8_chamfers;
MCNUM chamfers_z = mccbender8_chamfers_z;
MCNUM chamfers_lr = mccbender8_chamfers_lr;
MCNUM chamfers_tb = mccbender8_chamfers_tb;
MCNUM nelements = mccbender8_nelements;
MCNUM nu = mccbender8_nu;
MCNUM phase = mccbender8_phase;
char* reflect = mccbender8_reflect;
#line 571 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide_gravity.comp"
{

  if (l > 0 && nelements > 0) {
    int i,j,n;
    double x1,x2,x3,x4;
    double y1,y2,y3,y4;
    double nel = (nelements > 11 ? 11 : nelements);

    
    for (n=0; n<nel; n++)
    {
      double z0, z1;
      z0 =     n*(l/nel);
      z1 = (n+1)*(l/nel);

      for(j = 0; j < nhslit; j++)
      {
        y1 = j*(GVars.h1c+d)         - h1/2.0;
        y2 = j*(GVars.h2c+d)         - h2/2.0;
        y3 = (j+1)*(GVars.h1c+d) - d - h1/2.0;
        y4 = (j+1)*(GVars.h2c+d) - d - h2/2.0;
        for(i = 0; i < nslit; i++)
        {
          x1 = i*(GVars.w1c+d)         - w1/2.0;
          x2 = i*(GVars.w2c+d)         - w2/2.0;
          x3 = (i+1)*(GVars.w1c+d) - d - w1/2.0;
          x4 = (i+1)*(GVars.w2c+d) - d - w2/2.0;
          multiline(5,
                    x1, y1, z0,
                    x2, y2, z1,
                    x2, y4, z1,
                    x1, y3, z0,
                    x1, y1, z0);
          multiline(5,
                    x3, y1, z0,
                    x4, y2, z1,
                    x4, y4, z1,
                    x3, y3, z0,
                    x3, y1, z0);
        }
        line(-w1/2.0, y1, z0, w1/2.0, y1, z0);
        line(-w2/2.0, y2, z1, w2/2.0, y2, z1);
      }
    }

    if (nu || phase) {
      double radius = sqrt(w1*w1+l*l);
      /* cylinder top/center/bottom  */
      circle("xz", 0,-h1/2,l/2,radius);
      circle("xz", 0,0    ,l/2,radius);
      circle("xz", 0, h1/2,l/2,radius);
    }
  }
  else {
    /* A bit ugly; hard-coded dimensions. */
    
    line(0,0,0,0.2,0,0);
    line(0,0,0,0,0.2,0);
    line(0,0,0,0,0,0.2);
  }

}
#line 24267 "./ISIS_SANS2d.c"
}   /* End of bender8=Guide_gravity() SETTING parameter declarations. */
#undef pTable
#undef GVars
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'bender9'. */
  SIG_MESSAGE("bender9 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "bender9");
#define mccompcurname  bender9
#define mccompcurtype  Guide_gravity
#define mccompcurindex 13
#define GVars mccbender9_GVars
#define pTable mccbender9_pTable
{   /* Declarations of bender9=Guide_gravity() SETTING parameters. */
MCNUM w1 = mccbender9_w1;
MCNUM h1 = mccbender9_h1;
MCNUM w2 = mccbender9_w2;
MCNUM h2 = mccbender9_h2;
MCNUM l = mccbender9_l;
MCNUM R0 = mccbender9_R0;
MCNUM Qc = mccbender9_Qc;
MCNUM alpha = mccbender9_alpha;
MCNUM m = mccbender9_m;
MCNUM W = mccbender9_W;
MCNUM nslit = mccbender9_nslit;
MCNUM d = mccbender9_d;
MCNUM mleft = mccbender9_mleft;
MCNUM mright = mccbender9_mright;
MCNUM mtop = mccbender9_mtop;
MCNUM mbottom = mccbender9_mbottom;
MCNUM nhslit = mccbender9_nhslit;
MCNUM G = mccbender9_G;
MCNUM aleft = mccbender9_aleft;
MCNUM aright = mccbender9_aright;
MCNUM atop = mccbender9_atop;
MCNUM abottom = mccbender9_abottom;
MCNUM wavy = mccbender9_wavy;
MCNUM wavy_z = mccbender9_wavy_z;
MCNUM wavy_tb = mccbender9_wavy_tb;
MCNUM wavy_lr = mccbender9_wavy_lr;
MCNUM chamfers = mccbender9_chamfers;
MCNUM chamfers_z = mccbender9_chamfers_z;
MCNUM chamfers_lr = mccbender9_chamfers_lr;
MCNUM chamfers_tb = mccbender9_chamfers_tb;
MCNUM nelements = mccbender9_nelements;
MCNUM nu = mccbender9_nu;
MCNUM phase = mccbender9_phase;
char* reflect = mccbender9_reflect;
#line 571 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide_gravity.comp"
{

  if (l > 0 && nelements > 0) {
    int i,j,n;
    double x1,x2,x3,x4;
    double y1,y2,y3,y4;
    double nel = (nelements > 11 ? 11 : nelements);

    
    for (n=0; n<nel; n++)
    {
      double z0, z1;
      z0 =     n*(l/nel);
      z1 = (n+1)*(l/nel);

      for(j = 0; j < nhslit; j++)
      {
        y1 = j*(GVars.h1c+d)         - h1/2.0;
        y2 = j*(GVars.h2c+d)         - h2/2.0;
        y3 = (j+1)*(GVars.h1c+d) - d - h1/2.0;
        y4 = (j+1)*(GVars.h2c+d) - d - h2/2.0;
        for(i = 0; i < nslit; i++)
        {
          x1 = i*(GVars.w1c+d)         - w1/2.0;
          x2 = i*(GVars.w2c+d)         - w2/2.0;
          x3 = (i+1)*(GVars.w1c+d) - d - w1/2.0;
          x4 = (i+1)*(GVars.w2c+d) - d - w2/2.0;
          multiline(5,
                    x1, y1, z0,
                    x2, y2, z1,
                    x2, y4, z1,
                    x1, y3, z0,
                    x1, y1, z0);
          multiline(5,
                    x3, y1, z0,
                    x4, y2, z1,
                    x4, y4, z1,
                    x3, y3, z0,
                    x3, y1, z0);
        }
        line(-w1/2.0, y1, z0, w1/2.0, y1, z0);
        line(-w2/2.0, y2, z1, w2/2.0, y2, z1);
      }
    }

    if (nu || phase) {
      double radius = sqrt(w1*w1+l*l);
      /* cylinder top/center/bottom  */
      circle("xz", 0,-h1/2,l/2,radius);
      circle("xz", 0,0    ,l/2,radius);
      circle("xz", 0, h1/2,l/2,radius);
    }
  }
  else {
    /* A bit ugly; hard-coded dimensions. */
    
    line(0,0,0,0.2,0,0);
    line(0,0,0,0,0.2,0);
    line(0,0,0,0,0,0.2);
  }

}
#line 24381 "./ISIS_SANS2d.c"
}   /* End of bender9=Guide_gravity() SETTING parameter declarations. */
#undef pTable
#undef GVars
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'bender10'. */
  SIG_MESSAGE("bender10 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "bender10");
#define mccompcurname  bender10
#define mccompcurtype  Guide_gravity
#define mccompcurindex 14
#define GVars mccbender10_GVars
#define pTable mccbender10_pTable
{   /* Declarations of bender10=Guide_gravity() SETTING parameters. */
MCNUM w1 = mccbender10_w1;
MCNUM h1 = mccbender10_h1;
MCNUM w2 = mccbender10_w2;
MCNUM h2 = mccbender10_h2;
MCNUM l = mccbender10_l;
MCNUM R0 = mccbender10_R0;
MCNUM Qc = mccbender10_Qc;
MCNUM alpha = mccbender10_alpha;
MCNUM m = mccbender10_m;
MCNUM W = mccbender10_W;
MCNUM nslit = mccbender10_nslit;
MCNUM d = mccbender10_d;
MCNUM mleft = mccbender10_mleft;
MCNUM mright = mccbender10_mright;
MCNUM mtop = mccbender10_mtop;
MCNUM mbottom = mccbender10_mbottom;
MCNUM nhslit = mccbender10_nhslit;
MCNUM G = mccbender10_G;
MCNUM aleft = mccbender10_aleft;
MCNUM aright = mccbender10_aright;
MCNUM atop = mccbender10_atop;
MCNUM abottom = mccbender10_abottom;
MCNUM wavy = mccbender10_wavy;
MCNUM wavy_z = mccbender10_wavy_z;
MCNUM wavy_tb = mccbender10_wavy_tb;
MCNUM wavy_lr = mccbender10_wavy_lr;
MCNUM chamfers = mccbender10_chamfers;
MCNUM chamfers_z = mccbender10_chamfers_z;
MCNUM chamfers_lr = mccbender10_chamfers_lr;
MCNUM chamfers_tb = mccbender10_chamfers_tb;
MCNUM nelements = mccbender10_nelements;
MCNUM nu = mccbender10_nu;
MCNUM phase = mccbender10_phase;
char* reflect = mccbender10_reflect;
#line 571 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide_gravity.comp"
{

  if (l > 0 && nelements > 0) {
    int i,j,n;
    double x1,x2,x3,x4;
    double y1,y2,y3,y4;
    double nel = (nelements > 11 ? 11 : nelements);

    
    for (n=0; n<nel; n++)
    {
      double z0, z1;
      z0 =     n*(l/nel);
      z1 = (n+1)*(l/nel);

      for(j = 0; j < nhslit; j++)
      {
        y1 = j*(GVars.h1c+d)         - h1/2.0;
        y2 = j*(GVars.h2c+d)         - h2/2.0;
        y3 = (j+1)*(GVars.h1c+d) - d - h1/2.0;
        y4 = (j+1)*(GVars.h2c+d) - d - h2/2.0;
        for(i = 0; i < nslit; i++)
        {
          x1 = i*(GVars.w1c+d)         - w1/2.0;
          x2 = i*(GVars.w2c+d)         - w2/2.0;
          x3 = (i+1)*(GVars.w1c+d) - d - w1/2.0;
          x4 = (i+1)*(GVars.w2c+d) - d - w2/2.0;
          multiline(5,
                    x1, y1, z0,
                    x2, y2, z1,
                    x2, y4, z1,
                    x1, y3, z0,
                    x1, y1, z0);
          multiline(5,
                    x3, y1, z0,
                    x4, y2, z1,
                    x4, y4, z1,
                    x3, y3, z0,
                    x3, y1, z0);
        }
        line(-w1/2.0, y1, z0, w1/2.0, y1, z0);
        line(-w2/2.0, y2, z1, w2/2.0, y2, z1);
      }
    }

    if (nu || phase) {
      double radius = sqrt(w1*w1+l*l);
      /* cylinder top/center/bottom  */
      circle("xz", 0,-h1/2,l/2,radius);
      circle("xz", 0,0    ,l/2,radius);
      circle("xz", 0, h1/2,l/2,radius);
    }
  }
  else {
    /* A bit ugly; hard-coded dimensions. */
    
    line(0,0,0,0.2,0,0);
    line(0,0,0,0,0.2,0);
    line(0,0,0,0,0,0.2);
  }

}
#line 24495 "./ISIS_SANS2d.c"
}   /* End of bender10=Guide_gravity() SETTING parameter declarations. */
#undef pTable
#undef GVars
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'lmonb'. */
  SIG_MESSAGE("lmonb (McDisplay)");
  printf("MCDISPLAY: component %s\n", "lmonb");
#define mccompcurname  lmonb
#define mccompcurtype  L_monitor
#define mccompcurindex 15
#define nL mcclmonb_nL
#define L_N mcclmonb_L_N
#define L_p mcclmonb_L_p
#define L_p2 mcclmonb_L_p2
{   /* Declarations of lmonb=L_monitor() SETTING parameters. */
char* filename = mcclmonb_filename;
MCNUM xmin = mcclmonb_xmin;
MCNUM xmax = mcclmonb_xmax;
MCNUM ymin = mcclmonb_ymin;
MCNUM ymax = mcclmonb_ymax;
MCNUM xwidth = mcclmonb_xwidth;
MCNUM yheight = mcclmonb_yheight;
MCNUM Lmin = mcclmonb_Lmin;
MCNUM Lmax = mcclmonb_Lmax;
MCNUM restore_neutron = mcclmonb_restore_neutron;
int nowritefile = mcclmonb_nowritefile;
#line 120 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/L_monitor.comp"
{
  
  multiline(5, (double)xmin, (double)ymin, 0.0,
               (double)xmax, (double)ymin, 0.0,
               (double)xmax, (double)ymax, 0.0,
               (double)xmin, (double)ymax, 0.0,
               (double)xmin, (double)ymin, 0.0);
}
#line 24534 "./ISIS_SANS2d.c"
}   /* End of lmonb=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'psd2'. */
  SIG_MESSAGE("psd2 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "psd2");
#define mccompcurname  psd2
#define mccompcurtype  PSD_monitor
#define mccompcurindex 16
#define PSD_N mccpsd2_PSD_N
#define PSD_p mccpsd2_PSD_p
#define PSD_p2 mccpsd2_PSD_p2
{   /* Declarations of psd2=PSD_monitor() SETTING parameters. */
int nx = mccpsd2_nx;
int ny = mccpsd2_ny;
char* filename = mccpsd2_filename;
MCNUM xmin = mccpsd2_xmin;
MCNUM xmax = mccpsd2_xmax;
MCNUM ymin = mccpsd2_ymin;
MCNUM ymax = mccpsd2_ymax;
MCNUM xwidth = mccpsd2_xwidth;
MCNUM yheight = mccpsd2_yheight;
MCNUM restore_neutron = mccpsd2_restore_neutron;
#line 129 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
{
  multiline(5,
    (double)xmin, (double)ymin, 0.0,
    (double)xmax, (double)ymin, 0.0,
    (double)xmax, (double)ymax, 0.0,
    (double)xmin, (double)ymax, 0.0,
    (double)xmin, (double)ymin, 0.0);
}
#line 24573 "./ISIS_SANS2d.c"
}   /* End of psd2=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'guide_in'. */
  SIG_MESSAGE("guide_in (McDisplay)");
  printf("MCDISPLAY: component %s\n", "guide_in");
#define mccompcurname  guide_in
#define mccompcurtype  Slit
#define mccompcurindex 17
{   /* Declarations of guide_in=Slit() SETTING parameters. */
MCNUM xmin = mccguide_in_xmin;
MCNUM xmax = mccguide_in_xmax;
MCNUM ymin = mccguide_in_ymin;
MCNUM ymax = mccguide_in_ymax;
MCNUM radius = mccguide_in_radius;
MCNUM xwidth = mccguide_in_xwidth;
MCNUM yheight = mccguide_in_yheight;
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
#line 24619 "./ISIS_SANS2d.c"
}   /* End of guide_in=Slit() SETTING parameter declarations. */
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'guide_straight1'. */
  SIG_MESSAGE("guide_straight1 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "guide_straight1");
#define mccompcurname  guide_straight1
#define mccompcurtype  Guide_gravity
#define mccompcurindex 18
#define GVars mccguide_straight1_GVars
#define pTable mccguide_straight1_pTable
{   /* Declarations of guide_straight1=Guide_gravity() SETTING parameters. */
MCNUM w1 = mccguide_straight1_w1;
MCNUM h1 = mccguide_straight1_h1;
MCNUM w2 = mccguide_straight1_w2;
MCNUM h2 = mccguide_straight1_h2;
MCNUM l = mccguide_straight1_l;
MCNUM R0 = mccguide_straight1_R0;
MCNUM Qc = mccguide_straight1_Qc;
MCNUM alpha = mccguide_straight1_alpha;
MCNUM m = mccguide_straight1_m;
MCNUM W = mccguide_straight1_W;
MCNUM nslit = mccguide_straight1_nslit;
MCNUM d = mccguide_straight1_d;
MCNUM mleft = mccguide_straight1_mleft;
MCNUM mright = mccguide_straight1_mright;
MCNUM mtop = mccguide_straight1_mtop;
MCNUM mbottom = mccguide_straight1_mbottom;
MCNUM nhslit = mccguide_straight1_nhslit;
MCNUM G = mccguide_straight1_G;
MCNUM aleft = mccguide_straight1_aleft;
MCNUM aright = mccguide_straight1_aright;
MCNUM atop = mccguide_straight1_atop;
MCNUM abottom = mccguide_straight1_abottom;
MCNUM wavy = mccguide_straight1_wavy;
MCNUM wavy_z = mccguide_straight1_wavy_z;
MCNUM wavy_tb = mccguide_straight1_wavy_tb;
MCNUM wavy_lr = mccguide_straight1_wavy_lr;
MCNUM chamfers = mccguide_straight1_chamfers;
MCNUM chamfers_z = mccguide_straight1_chamfers_z;
MCNUM chamfers_lr = mccguide_straight1_chamfers_lr;
MCNUM chamfers_tb = mccguide_straight1_chamfers_tb;
MCNUM nelements = mccguide_straight1_nelements;
MCNUM nu = mccguide_straight1_nu;
MCNUM phase = mccguide_straight1_phase;
char* reflect = mccguide_straight1_reflect;
#line 571 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide_gravity.comp"
{

  if (l > 0 && nelements > 0) {
    int i,j,n;
    double x1,x2,x3,x4;
    double y1,y2,y3,y4;
    double nel = (nelements > 11 ? 11 : nelements);

    
    for (n=0; n<nel; n++)
    {
      double z0, z1;
      z0 =     n*(l/nel);
      z1 = (n+1)*(l/nel);

      for(j = 0; j < nhslit; j++)
      {
        y1 = j*(GVars.h1c+d)         - h1/2.0;
        y2 = j*(GVars.h2c+d)         - h2/2.0;
        y3 = (j+1)*(GVars.h1c+d) - d - h1/2.0;
        y4 = (j+1)*(GVars.h2c+d) - d - h2/2.0;
        for(i = 0; i < nslit; i++)
        {
          x1 = i*(GVars.w1c+d)         - w1/2.0;
          x2 = i*(GVars.w2c+d)         - w2/2.0;
          x3 = (i+1)*(GVars.w1c+d) - d - w1/2.0;
          x4 = (i+1)*(GVars.w2c+d) - d - w2/2.0;
          multiline(5,
                    x1, y1, z0,
                    x2, y2, z1,
                    x2, y4, z1,
                    x1, y3, z0,
                    x1, y1, z0);
          multiline(5,
                    x3, y1, z0,
                    x4, y2, z1,
                    x4, y4, z1,
                    x3, y3, z0,
                    x3, y1, z0);
        }
        line(-w1/2.0, y1, z0, w1/2.0, y1, z0);
        line(-w2/2.0, y2, z1, w2/2.0, y2, z1);
      }
    }

    if (nu || phase) {
      double radius = sqrt(w1*w1+l*l);
      /* cylinder top/center/bottom  */
      circle("xz", 0,-h1/2,l/2,radius);
      circle("xz", 0,0    ,l/2,radius);
      circle("xz", 0, h1/2,l/2,radius);
    }
  }
  else {
    /* A bit ugly; hard-coded dimensions. */
    
    line(0,0,0,0.2,0,0);
    line(0,0,0,0,0.2,0);
    line(0,0,0,0,0,0.2);
  }

}
#line 24731 "./ISIS_SANS2d.c"
}   /* End of guide_straight1=Guide_gravity() SETTING parameter declarations. */
#undef pTable
#undef GVars
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'guide_straight2'. */
  SIG_MESSAGE("guide_straight2 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "guide_straight2");
#define mccompcurname  guide_straight2
#define mccompcurtype  Guide_gravity
#define mccompcurindex 19
#define GVars mccguide_straight2_GVars
#define pTable mccguide_straight2_pTable
{   /* Declarations of guide_straight2=Guide_gravity() SETTING parameters. */
MCNUM w1 = mccguide_straight2_w1;
MCNUM h1 = mccguide_straight2_h1;
MCNUM w2 = mccguide_straight2_w2;
MCNUM h2 = mccguide_straight2_h2;
MCNUM l = mccguide_straight2_l;
MCNUM R0 = mccguide_straight2_R0;
MCNUM Qc = mccguide_straight2_Qc;
MCNUM alpha = mccguide_straight2_alpha;
MCNUM m = mccguide_straight2_m;
MCNUM W = mccguide_straight2_W;
MCNUM nslit = mccguide_straight2_nslit;
MCNUM d = mccguide_straight2_d;
MCNUM mleft = mccguide_straight2_mleft;
MCNUM mright = mccguide_straight2_mright;
MCNUM mtop = mccguide_straight2_mtop;
MCNUM mbottom = mccguide_straight2_mbottom;
MCNUM nhslit = mccguide_straight2_nhslit;
MCNUM G = mccguide_straight2_G;
MCNUM aleft = mccguide_straight2_aleft;
MCNUM aright = mccguide_straight2_aright;
MCNUM atop = mccguide_straight2_atop;
MCNUM abottom = mccguide_straight2_abottom;
MCNUM wavy = mccguide_straight2_wavy;
MCNUM wavy_z = mccguide_straight2_wavy_z;
MCNUM wavy_tb = mccguide_straight2_wavy_tb;
MCNUM wavy_lr = mccguide_straight2_wavy_lr;
MCNUM chamfers = mccguide_straight2_chamfers;
MCNUM chamfers_z = mccguide_straight2_chamfers_z;
MCNUM chamfers_lr = mccguide_straight2_chamfers_lr;
MCNUM chamfers_tb = mccguide_straight2_chamfers_tb;
MCNUM nelements = mccguide_straight2_nelements;
MCNUM nu = mccguide_straight2_nu;
MCNUM phase = mccguide_straight2_phase;
char* reflect = mccguide_straight2_reflect;
#line 571 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide_gravity.comp"
{

  if (l > 0 && nelements > 0) {
    int i,j,n;
    double x1,x2,x3,x4;
    double y1,y2,y3,y4;
    double nel = (nelements > 11 ? 11 : nelements);

    
    for (n=0; n<nel; n++)
    {
      double z0, z1;
      z0 =     n*(l/nel);
      z1 = (n+1)*(l/nel);

      for(j = 0; j < nhslit; j++)
      {
        y1 = j*(GVars.h1c+d)         - h1/2.0;
        y2 = j*(GVars.h2c+d)         - h2/2.0;
        y3 = (j+1)*(GVars.h1c+d) - d - h1/2.0;
        y4 = (j+1)*(GVars.h2c+d) - d - h2/2.0;
        for(i = 0; i < nslit; i++)
        {
          x1 = i*(GVars.w1c+d)         - w1/2.0;
          x2 = i*(GVars.w2c+d)         - w2/2.0;
          x3 = (i+1)*(GVars.w1c+d) - d - w1/2.0;
          x4 = (i+1)*(GVars.w2c+d) - d - w2/2.0;
          multiline(5,
                    x1, y1, z0,
                    x2, y2, z1,
                    x2, y4, z1,
                    x1, y3, z0,
                    x1, y1, z0);
          multiline(5,
                    x3, y1, z0,
                    x4, y2, z1,
                    x4, y4, z1,
                    x3, y3, z0,
                    x3, y1, z0);
        }
        line(-w1/2.0, y1, z0, w1/2.0, y1, z0);
        line(-w2/2.0, y2, z1, w2/2.0, y2, z1);
      }
    }

    if (nu || phase) {
      double radius = sqrt(w1*w1+l*l);
      /* cylinder top/center/bottom  */
      circle("xz", 0,-h1/2,l/2,radius);
      circle("xz", 0,0    ,l/2,radius);
      circle("xz", 0, h1/2,l/2,radius);
    }
  }
  else {
    /* A bit ugly; hard-coded dimensions. */
    
    line(0,0,0,0.2,0,0);
    line(0,0,0,0,0.2,0);
    line(0,0,0,0,0,0.2);
  }

}
#line 24845 "./ISIS_SANS2d.c"
}   /* End of guide_straight2=Guide_gravity() SETTING parameter declarations. */
#undef pTable
#undef GVars
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'guide_straight3'. */
  SIG_MESSAGE("guide_straight3 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "guide_straight3");
#define mccompcurname  guide_straight3
#define mccompcurtype  Guide_gravity
#define mccompcurindex 20
#define GVars mccguide_straight3_GVars
#define pTable mccguide_straight3_pTable
{   /* Declarations of guide_straight3=Guide_gravity() SETTING parameters. */
MCNUM w1 = mccguide_straight3_w1;
MCNUM h1 = mccguide_straight3_h1;
MCNUM w2 = mccguide_straight3_w2;
MCNUM h2 = mccguide_straight3_h2;
MCNUM l = mccguide_straight3_l;
MCNUM R0 = mccguide_straight3_R0;
MCNUM Qc = mccguide_straight3_Qc;
MCNUM alpha = mccguide_straight3_alpha;
MCNUM m = mccguide_straight3_m;
MCNUM W = mccguide_straight3_W;
MCNUM nslit = mccguide_straight3_nslit;
MCNUM d = mccguide_straight3_d;
MCNUM mleft = mccguide_straight3_mleft;
MCNUM mright = mccguide_straight3_mright;
MCNUM mtop = mccguide_straight3_mtop;
MCNUM mbottom = mccguide_straight3_mbottom;
MCNUM nhslit = mccguide_straight3_nhslit;
MCNUM G = mccguide_straight3_G;
MCNUM aleft = mccguide_straight3_aleft;
MCNUM aright = mccguide_straight3_aright;
MCNUM atop = mccguide_straight3_atop;
MCNUM abottom = mccguide_straight3_abottom;
MCNUM wavy = mccguide_straight3_wavy;
MCNUM wavy_z = mccguide_straight3_wavy_z;
MCNUM wavy_tb = mccguide_straight3_wavy_tb;
MCNUM wavy_lr = mccguide_straight3_wavy_lr;
MCNUM chamfers = mccguide_straight3_chamfers;
MCNUM chamfers_z = mccguide_straight3_chamfers_z;
MCNUM chamfers_lr = mccguide_straight3_chamfers_lr;
MCNUM chamfers_tb = mccguide_straight3_chamfers_tb;
MCNUM nelements = mccguide_straight3_nelements;
MCNUM nu = mccguide_straight3_nu;
MCNUM phase = mccguide_straight3_phase;
char* reflect = mccguide_straight3_reflect;
#line 571 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide_gravity.comp"
{

  if (l > 0 && nelements > 0) {
    int i,j,n;
    double x1,x2,x3,x4;
    double y1,y2,y3,y4;
    double nel = (nelements > 11 ? 11 : nelements);

    
    for (n=0; n<nel; n++)
    {
      double z0, z1;
      z0 =     n*(l/nel);
      z1 = (n+1)*(l/nel);

      for(j = 0; j < nhslit; j++)
      {
        y1 = j*(GVars.h1c+d)         - h1/2.0;
        y2 = j*(GVars.h2c+d)         - h2/2.0;
        y3 = (j+1)*(GVars.h1c+d) - d - h1/2.0;
        y4 = (j+1)*(GVars.h2c+d) - d - h2/2.0;
        for(i = 0; i < nslit; i++)
        {
          x1 = i*(GVars.w1c+d)         - w1/2.0;
          x2 = i*(GVars.w2c+d)         - w2/2.0;
          x3 = (i+1)*(GVars.w1c+d) - d - w1/2.0;
          x4 = (i+1)*(GVars.w2c+d) - d - w2/2.0;
          multiline(5,
                    x1, y1, z0,
                    x2, y2, z1,
                    x2, y4, z1,
                    x1, y3, z0,
                    x1, y1, z0);
          multiline(5,
                    x3, y1, z0,
                    x4, y2, z1,
                    x4, y4, z1,
                    x3, y3, z0,
                    x3, y1, z0);
        }
        line(-w1/2.0, y1, z0, w1/2.0, y1, z0);
        line(-w2/2.0, y2, z1, w2/2.0, y2, z1);
      }
    }

    if (nu || phase) {
      double radius = sqrt(w1*w1+l*l);
      /* cylinder top/center/bottom  */
      circle("xz", 0,-h1/2,l/2,radius);
      circle("xz", 0,0    ,l/2,radius);
      circle("xz", 0, h1/2,l/2,radius);
    }
  }
  else {
    /* A bit ugly; hard-coded dimensions. */
    
    line(0,0,0,0.2,0,0);
    line(0,0,0,0,0.2,0);
    line(0,0,0,0,0,0.2);
  }

}
#line 24959 "./ISIS_SANS2d.c"
}   /* End of guide_straight3=Guide_gravity() SETTING parameter declarations. */
#undef pTable
#undef GVars
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'guide_straight4'. */
  SIG_MESSAGE("guide_straight4 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "guide_straight4");
#define mccompcurname  guide_straight4
#define mccompcurtype  Guide_gravity
#define mccompcurindex 21
#define GVars mccguide_straight4_GVars
#define pTable mccguide_straight4_pTable
{   /* Declarations of guide_straight4=Guide_gravity() SETTING parameters. */
MCNUM w1 = mccguide_straight4_w1;
MCNUM h1 = mccguide_straight4_h1;
MCNUM w2 = mccguide_straight4_w2;
MCNUM h2 = mccguide_straight4_h2;
MCNUM l = mccguide_straight4_l;
MCNUM R0 = mccguide_straight4_R0;
MCNUM Qc = mccguide_straight4_Qc;
MCNUM alpha = mccguide_straight4_alpha;
MCNUM m = mccguide_straight4_m;
MCNUM W = mccguide_straight4_W;
MCNUM nslit = mccguide_straight4_nslit;
MCNUM d = mccguide_straight4_d;
MCNUM mleft = mccguide_straight4_mleft;
MCNUM mright = mccguide_straight4_mright;
MCNUM mtop = mccguide_straight4_mtop;
MCNUM mbottom = mccguide_straight4_mbottom;
MCNUM nhslit = mccguide_straight4_nhslit;
MCNUM G = mccguide_straight4_G;
MCNUM aleft = mccguide_straight4_aleft;
MCNUM aright = mccguide_straight4_aright;
MCNUM atop = mccguide_straight4_atop;
MCNUM abottom = mccguide_straight4_abottom;
MCNUM wavy = mccguide_straight4_wavy;
MCNUM wavy_z = mccguide_straight4_wavy_z;
MCNUM wavy_tb = mccguide_straight4_wavy_tb;
MCNUM wavy_lr = mccguide_straight4_wavy_lr;
MCNUM chamfers = mccguide_straight4_chamfers;
MCNUM chamfers_z = mccguide_straight4_chamfers_z;
MCNUM chamfers_lr = mccguide_straight4_chamfers_lr;
MCNUM chamfers_tb = mccguide_straight4_chamfers_tb;
MCNUM nelements = mccguide_straight4_nelements;
MCNUM nu = mccguide_straight4_nu;
MCNUM phase = mccguide_straight4_phase;
char* reflect = mccguide_straight4_reflect;
#line 571 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide_gravity.comp"
{

  if (l > 0 && nelements > 0) {
    int i,j,n;
    double x1,x2,x3,x4;
    double y1,y2,y3,y4;
    double nel = (nelements > 11 ? 11 : nelements);

    
    for (n=0; n<nel; n++)
    {
      double z0, z1;
      z0 =     n*(l/nel);
      z1 = (n+1)*(l/nel);

      for(j = 0; j < nhslit; j++)
      {
        y1 = j*(GVars.h1c+d)         - h1/2.0;
        y2 = j*(GVars.h2c+d)         - h2/2.0;
        y3 = (j+1)*(GVars.h1c+d) - d - h1/2.0;
        y4 = (j+1)*(GVars.h2c+d) - d - h2/2.0;
        for(i = 0; i < nslit; i++)
        {
          x1 = i*(GVars.w1c+d)         - w1/2.0;
          x2 = i*(GVars.w2c+d)         - w2/2.0;
          x3 = (i+1)*(GVars.w1c+d) - d - w1/2.0;
          x4 = (i+1)*(GVars.w2c+d) - d - w2/2.0;
          multiline(5,
                    x1, y1, z0,
                    x2, y2, z1,
                    x2, y4, z1,
                    x1, y3, z0,
                    x1, y1, z0);
          multiline(5,
                    x3, y1, z0,
                    x4, y2, z1,
                    x4, y4, z1,
                    x3, y3, z0,
                    x3, y1, z0);
        }
        line(-w1/2.0, y1, z0, w1/2.0, y1, z0);
        line(-w2/2.0, y2, z1, w2/2.0, y2, z1);
      }
    }

    if (nu || phase) {
      double radius = sqrt(w1*w1+l*l);
      /* cylinder top/center/bottom  */
      circle("xz", 0,-h1/2,l/2,radius);
      circle("xz", 0,0    ,l/2,radius);
      circle("xz", 0, h1/2,l/2,radius);
    }
  }
  else {
    /* A bit ugly; hard-coded dimensions. */
    
    line(0,0,0,0.2,0,0);
    line(0,0,0,0,0.2,0);
    line(0,0,0,0,0,0.2);
  }

}
#line 25073 "./ISIS_SANS2d.c"
}   /* End of guide_straight4=Guide_gravity() SETTING parameter declarations. */
#undef pTable
#undef GVars
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'psd3'. */
  SIG_MESSAGE("psd3 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "psd3");
#define mccompcurname  psd3
#define mccompcurtype  PSD_monitor
#define mccompcurindex 22
#define PSD_N mccpsd3_PSD_N
#define PSD_p mccpsd3_PSD_p
#define PSD_p2 mccpsd3_PSD_p2
{   /* Declarations of psd3=PSD_monitor() SETTING parameters. */
int nx = mccpsd3_nx;
int ny = mccpsd3_ny;
char* filename = mccpsd3_filename;
MCNUM xmin = mccpsd3_xmin;
MCNUM xmax = mccpsd3_xmax;
MCNUM ymin = mccpsd3_ymin;
MCNUM ymax = mccpsd3_ymax;
MCNUM xwidth = mccpsd3_xwidth;
MCNUM yheight = mccpsd3_yheight;
MCNUM restore_neutron = mccpsd3_restore_neutron;
#line 129 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
{
  multiline(5,
    (double)xmin, (double)ymin, 0.0,
    (double)xmax, (double)ymin, 0.0,
    (double)xmax, (double)ymax, 0.0,
    (double)xmin, (double)ymax, 0.0,
    (double)xmin, (double)ymin, 0.0);
}
#line 25110 "./ISIS_SANS2d.c"
}   /* End of psd3=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'aperture1'. */
  SIG_MESSAGE("aperture1 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "aperture1");
#define mccompcurname  aperture1
#define mccompcurtype  Slit
#define mccompcurindex 23
{   /* Declarations of aperture1=Slit() SETTING parameters. */
MCNUM xmin = mccaperture1_xmin;
MCNUM xmax = mccaperture1_xmax;
MCNUM ymin = mccaperture1_ymin;
MCNUM ymax = mccaperture1_ymax;
MCNUM radius = mccaperture1_radius;
MCNUM xwidth = mccaperture1_xwidth;
MCNUM yheight = mccaperture1_yheight;
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
#line 25156 "./ISIS_SANS2d.c"
}   /* End of aperture1=Slit() SETTING parameter declarations. */
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'lmonitor2'. */
  SIG_MESSAGE("lmonitor2 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "lmonitor2");
#define mccompcurname  lmonitor2
#define mccompcurtype  L_monitor
#define mccompcurindex 24
#define nL mcclmonitor2_nL
#define L_N mcclmonitor2_L_N
#define L_p mcclmonitor2_L_p
#define L_p2 mcclmonitor2_L_p2
{   /* Declarations of lmonitor2=L_monitor() SETTING parameters. */
char* filename = mcclmonitor2_filename;
MCNUM xmin = mcclmonitor2_xmin;
MCNUM xmax = mcclmonitor2_xmax;
MCNUM ymin = mcclmonitor2_ymin;
MCNUM ymax = mcclmonitor2_ymax;
MCNUM xwidth = mcclmonitor2_xwidth;
MCNUM yheight = mcclmonitor2_yheight;
MCNUM Lmin = mcclmonitor2_Lmin;
MCNUM Lmax = mcclmonitor2_Lmax;
MCNUM restore_neutron = mcclmonitor2_restore_neutron;
int nowritefile = mcclmonitor2_nowritefile;
#line 120 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/L_monitor.comp"
{
  
  multiline(5, (double)xmin, (double)ymin, 0.0,
               (double)xmax, (double)ymin, 0.0,
               (double)xmax, (double)ymax, 0.0,
               (double)xmin, (double)ymax, 0.0,
               (double)xmin, (double)ymin, 0.0);
}
#line 25193 "./ISIS_SANS2d.c"
}   /* End of lmonitor2=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'S6'. */
  SIG_MESSAGE("S6 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "S6");
#define mccompcurname  S6
#define mccompcurtype  Slit
#define mccompcurindex 25
{   /* Declarations of S6=Slit() SETTING parameters. */
MCNUM xmin = mccS6_xmin;
MCNUM xmax = mccS6_xmax;
MCNUM ymin = mccS6_ymin;
MCNUM ymax = mccS6_ymax;
MCNUM radius = mccS6_radius;
MCNUM xwidth = mccS6_xwidth;
MCNUM yheight = mccS6_yheight;
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
#line 25240 "./ISIS_SANS2d.c"
}   /* End of S6=Slit() SETTING parameter declarations. */
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'APERTURE2'. */
  SIG_MESSAGE("APERTURE2 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "APERTURE2");
#define mccompcurname  APERTURE2
#define mccompcurtype  Slit
#define mccompcurindex 26
{   /* Declarations of APERTURE2=Slit() SETTING parameters. */
MCNUM xmin = mccAPERTURE2_xmin;
MCNUM xmax = mccAPERTURE2_xmax;
MCNUM ymin = mccAPERTURE2_ymin;
MCNUM ymax = mccAPERTURE2_ymax;
MCNUM radius = mccAPERTURE2_radius;
MCNUM xwidth = mccAPERTURE2_xwidth;
MCNUM yheight = mccAPERTURE2_yheight;
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
#line 25283 "./ISIS_SANS2d.c"
}   /* End of APERTURE2=Slit() SETTING parameter declarations. */
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'lmon2'. */
  SIG_MESSAGE("lmon2 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "lmon2");
#define mccompcurname  lmon2
#define mccompcurtype  L_monitor
#define mccompcurindex 27
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
#line 120 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/L_monitor.comp"
{
  
  multiline(5, (double)xmin, (double)ymin, 0.0,
               (double)xmax, (double)ymin, 0.0,
               (double)xmax, (double)ymax, 0.0,
               (double)xmin, (double)ymax, 0.0,
               (double)xmin, (double)ymin, 0.0);
}
#line 25320 "./ISIS_SANS2d.c"
}   /* End of lmon2=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'psd4'. */
  SIG_MESSAGE("psd4 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "psd4");
#define mccompcurname  psd4
#define mccompcurtype  PSD_monitor
#define mccompcurindex 28
#define PSD_N mccpsd4_PSD_N
#define PSD_p mccpsd4_PSD_p
#define PSD_p2 mccpsd4_PSD_p2
{   /* Declarations of psd4=PSD_monitor() SETTING parameters. */
int nx = mccpsd4_nx;
int ny = mccpsd4_ny;
char* filename = mccpsd4_filename;
MCNUM xmin = mccpsd4_xmin;
MCNUM xmax = mccpsd4_xmax;
MCNUM ymin = mccpsd4_ymin;
MCNUM ymax = mccpsd4_ymax;
MCNUM xwidth = mccpsd4_xwidth;
MCNUM yheight = mccpsd4_yheight;
MCNUM restore_neutron = mccpsd4_restore_neutron;
#line 129 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
{
  multiline(5,
    (double)xmin, (double)ymin, 0.0,
    (double)xmax, (double)ymin, 0.0,
    (double)xmax, (double)ymax, 0.0,
    (double)xmin, (double)ymax, 0.0,
    (double)xmin, (double)ymin, 0.0);
}
#line 25359 "./ISIS_SANS2d.c"
}   /* End of psd4=PSD_monitor() SETTING parameter declarations. */
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
/* end of generated C code ./ISIS_SANS2d.c */
