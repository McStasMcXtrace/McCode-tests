/* Automatically generated file. Do not edit. 
 * Format:     ANSI C source code
 * Creator:    McStas <http://www.mcstas.org>
 * Instrument: /zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr (SE_example2)
 * Date:       Wed Feb 26 19:15:54 2020
 * File:       ./SE_example2.c
 * Compile:    cc -o SE_example2.out ./SE_example2.c 
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

#line 712 "./SE_example2.c"

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

#line 945 "./SE_example2.c"

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

#line 4977 "./SE_example2.c"

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

#line 5337 "./SE_example2.c"
#ifdef MC_TRACE_ENABLED
int mctraceenabled = 1;
#else
int mctraceenabled = 0;
#endif
#define MCSTAS "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../"
int mcdefaultmain = 1;
char mcinstrument_name[] = "SE_example2";
char mcinstrument_source[] = "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr";
char *mcinstrument_exe=NULL; /* will be set to argv[0] in main */
int main(int argc, char *argv[]){return mccode_main(argc, argv);}
void mcinit(void);
void mcraytrace(void);
void mcsave(FILE *);
void mcfinally(void);
void mcdisplay(void);

/* Shared user declarations for all components 'Pol_mirror'. */
#line 89 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Pol_mirror.comp"
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

#line 8339 "./SE_example2.c"

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

#line 8425 "./SE_example2.c"

/* Shared user declarations for all components 'Pol_FieldBox'. */
#line 44 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Pol_FieldBox.comp"

#line 8430 "./SE_example2.c"

/* Shared user declarations for all components 'V_sample'. */
#line 100 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../obsolete/V_sample.comp"
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
#line 8447 "./SE_example2.c"

/* Instrument parameters. */
MCNUM mcipPOL_ANGLE;
MCNUM mcipSAMPLE;
MCNUM mcipBguide;
MCNUM mcipBflip;
MCNUM mcipdBz;
MCNUM mcipLam;
MCNUM mcipdLam;

#define mcNUMIPAR 7
int mcnumipar = 7;
struct mcinputtable_struct mcinputtable[mcNUMIPAR+1] = {
  "POL_ANGLE", &mcipPOL_ANGLE, instr_type_double, "2", 
  "SAMPLE", &mcipSAMPLE, instr_type_double, "0", 
  "Bguide", &mcipBguide, instr_type_double, "0.1", 
  "Bflip", &mcipBflip, instr_type_double, "3e-4", 
  "dBz", &mcipdBz, instr_type_double, "0", 
  "Lam", &mcipLam, instr_type_double, "8", 
  "dLam", &mcipdLam, instr_type_double, "0.8", 
  NULL, NULL, instr_type_double, ""
};

/* User declarations from instrument definition. */
#define mccompcurname  SE_example2
#define mccompcurtype  INSTRUMENT
#define mccompcurindex 0
#define mcposaSE_example2 coords_set(0,0,0)
#define POL_ANGLE mcipPOL_ANGLE
#define SAMPLE mcipSAMPLE
#define Bguide mcipBguide
#define Bflip mcipBflip
#define dBz mcipdBz
#define Lam mcipLam
#define dLam mcipdLam
#line 42 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  double pol_radius;
#line 8485 "./SE_example2.c"
#undef dLam
#undef Lam
#undef dBz
#undef Bflip
#undef Bguide
#undef SAMPLE
#undef POL_ANGLE
#undef mcposaSE_example2
#undef mccompcurindex
#undef mccompcurtype
#undef mccompcurname

/* neutron state table at each component input (local coords) */
/* [x, y, z, vx, vy, vz, t, sx, sy, sz, p] */
MCNUM mccomp_storein[11*33];
/* Components position table (absolute and relative coords) */
Coords mccomp_posa[33];
Coords mccomp_posr[33];
/* Counter for each comp to check for inactive ones */
MCNUM  mcNCounter[33];
MCNUM  mcPCounter[33];
MCNUM  mcP2Counter[33];
#define mcNUMCOMP 32 /* number of components */
/* Counter for PROP ABSORB */
MCNUM  mcAbsorbProp[33];
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

/* Definition parameters for component 'lam_mon' [3]. */
#define mcclam_mon_nL 41
/* Setting parameters for component 'lam_mon' [3]. */
char mcclam_mon_filename[16384];
MCNUM mcclam_mon_xmin;
MCNUM mcclam_mon_xmax;
MCNUM mcclam_mon_ymin;
MCNUM mcclam_mon_ymax;
MCNUM mcclam_mon_xwidth;
MCNUM mcclam_mon_yheight;
MCNUM mcclam_mon_Lmin;
MCNUM mcclam_mon_Lmax;
MCNUM mcclam_mon_restore_neutron;
int mcclam_mon_nowritefile;

/* Definition parameters for component 'hdiv_mon' [5]. */
#define mcchdiv_mon_nh 21
/* Setting parameters for component 'hdiv_mon' [5]. */
char mcchdiv_mon_filename[16384];
MCNUM mcchdiv_mon_xmin;
MCNUM mcchdiv_mon_xmax;
MCNUM mcchdiv_mon_ymin;
MCNUM mcchdiv_mon_ymax;
MCNUM mcchdiv_mon_xwidth;
MCNUM mcchdiv_mon_yheight;
MCNUM mcchdiv_mon_h_maxdiv;
MCNUM mcchdiv_mon_restore_neutron;
int mcchdiv_mon_nowritefile;

/* Definition parameters for component 'polariser2' [6]. */
#define mccpolariser2_rUpFunc StdReflecFunc
#define mccpolariser2_rDownFunc rUpFunc
#define mccpolariser2_rUpPar { 0.99 , 0.0485 , 0.9 , 1 , 0.0002 }
#define mccpolariser2_rDownPar { 0.99 , 0.11 , 0.9 , 2 , 0.0185 }
#define mccpolariser2_useTables 0
#define mccpolariser2_zwidth 0.1
#define mccpolariser2_yheight 0.3
/* Setting parameters for component 'polariser2' [6]. */
MCNUM mccpolariser2_p_reflect;

/* Setting parameters for component 'psd_ap' [7]. */
int mccpsd_ap_nx;
int mccpsd_ap_ny;
char mccpsd_ap_filename[16384];
MCNUM mccpsd_ap_xmin;
MCNUM mccpsd_ap_xmax;
MCNUM mccpsd_ap_ymin;
MCNUM mccpsd_ap_ymax;
MCNUM mccpsd_ap_xwidth;
MCNUM mccpsd_ap_yheight;
MCNUM mccpsd_ap_restore_neutron;

/* Definition parameters for component 'pl_mon0' [8]. */
#define mccpl_mon0_xwidth 0.1
#define mccpl_mon0_yheight 0.01
#define mccpl_mon0_nL 20
#define mccpl_mon0_restore_neutron 1
#define mccpl_mon0_nowritefile 0
/* Setting parameters for component 'pl_mon0' [8]. */
char mccpl_mon0_filename[16384];
MCNUM mccpl_mon0_mx;
MCNUM mccpl_mon0_my;
MCNUM mccpl_mon0_mz;
MCNUM mccpl_mon0_Lmin;
MCNUM mccpl_mon0_Lmax;

/* Definition parameters for component 'flipper1' [9]. */
#define mccflipper1_fieldFunction const_magnetic_field
/* Setting parameters for component 'flipper1' [9]. */
MCNUM mccflipper1_xwidth;
MCNUM mccflipper1_yheight;
MCNUM mccflipper1_zdepth;
MCNUM mccflipper1_Bx;
MCNUM mccflipper1_By;
MCNUM mccflipper1_Bz;
char mccflipper1_filename[16384];

/* Definition parameters for component 'pl_mon1' [10]. */
#define mccpl_mon1_xwidth 0.1
#define mccpl_mon1_yheight 0.01
#define mccpl_mon1_nL 20
#define mccpl_mon1_restore_neutron 1
#define mccpl_mon1_nowritefile 0
/* Setting parameters for component 'pl_mon1' [10]. */
char mccpl_mon1_filename[16384];
MCNUM mccpl_mon1_mx;
MCNUM mccpl_mon1_my;
MCNUM mccpl_mon1_mz;
MCNUM mccpl_mon1_Lmin;
MCNUM mccpl_mon1_Lmax;

/* Definition parameters for component 'guide_field_1' [11]. */
#define mccguide_field_1_fieldFunction const_magnetic_field
/* Setting parameters for component 'guide_field_1' [11]. */
MCNUM mccguide_field_1_xwidth;
MCNUM mccguide_field_1_yheight;
MCNUM mccguide_field_1_zdepth;
MCNUM mccguide_field_1_Bx;
MCNUM mccguide_field_1_By;
MCNUM mccguide_field_1_Bz;
char mccguide_field_1_filename[16384];

/* Setting parameters for component 'guide_field_1_end' [12]. */
MCNUM mccguide_field_1_end_xmin;
MCNUM mccguide_field_1_end_xmax;
MCNUM mccguide_field_1_end_ymin;
MCNUM mccguide_field_1_end_ymax;
MCNUM mccguide_field_1_end_radius;
MCNUM mccguide_field_1_end_xwidth;
MCNUM mccguide_field_1_end_yheight;

/* Definition parameters for component 'pl_mon2' [13]. */
#define mccpl_mon2_xwidth 0.1
#define mccpl_mon2_yheight 0.1
#define mccpl_mon2_nL 20
#define mccpl_mon2_restore_neutron 1
#define mccpl_mon2_nowritefile 0
/* Setting parameters for component 'pl_mon2' [13]. */
char mccpl_mon2_filename[16384];
MCNUM mccpl_mon2_mx;
MCNUM mccpl_mon2_my;
MCNUM mccpl_mon2_mz;
MCNUM mccpl_mon2_Lmin;
MCNUM mccpl_mon2_Lmax;

/* Definition parameters for component 'flipper2_1' [14]. */
#define mccflipper2_1_fieldFunction const_magnetic_field
/* Setting parameters for component 'flipper2_1' [14]. */
MCNUM mccflipper2_1_xwidth;
MCNUM mccflipper2_1_yheight;
MCNUM mccflipper2_1_zdepth;
MCNUM mccflipper2_1_Bx;
MCNUM mccflipper2_1_By;
MCNUM mccflipper2_1_Bz;
char mccflipper2_1_filename[16384];

/* Definition parameters for component 'flipper2_2' [15]. */
#define mccflipper2_2_fieldFunction const_magnetic_field
/* Setting parameters for component 'flipper2_2' [15]. */
MCNUM mccflipper2_2_xwidth;
MCNUM mccflipper2_2_yheight;
MCNUM mccflipper2_2_zdepth;
MCNUM mccflipper2_2_Bx;
MCNUM mccflipper2_2_By;
MCNUM mccflipper2_2_Bz;
char mccflipper2_2_filename[16384];

/* Setting parameters for component 'flipper_2_2_end' [16]. */
MCNUM mccflipper_2_2_end_xmin;
MCNUM mccflipper_2_2_end_xmax;
MCNUM mccflipper_2_2_end_ymin;
MCNUM mccflipper_2_2_end_ymax;
MCNUM mccflipper_2_2_end_radius;
MCNUM mccflipper_2_2_end_xwidth;
MCNUM mccflipper_2_2_end_yheight;

/* Definition parameters for component 'pl_mon3' [18]. */
#define mccpl_mon3_xwidth 0.1
#define mccpl_mon3_yheight 0.1
#define mccpl_mon3_nL 20
#define mccpl_mon3_restore_neutron 1
#define mccpl_mon3_nowritefile 0
/* Setting parameters for component 'pl_mon3' [18]. */
char mccpl_mon3_filename[16384];
MCNUM mccpl_mon3_mx;
MCNUM mccpl_mon3_my;
MCNUM mccpl_mon3_mz;
MCNUM mccpl_mon3_Lmin;
MCNUM mccpl_mon3_Lmax;

/* Setting parameters for component 'sample0' [19]. */
MCNUM mccsample0_radius;
MCNUM mccsample0_thickness;
MCNUM mccsample0_zdepth;
MCNUM mccsample0_Vc;
MCNUM mccsample0_sigma_abs;
MCNUM mccsample0_sigma_inc;
MCNUM mccsample0_radius_i;
MCNUM mccsample0_radius_o;
MCNUM mccsample0_h;
MCNUM mccsample0_focus_r;
MCNUM mccsample0_pack;
MCNUM mccsample0_frac;
MCNUM mccsample0_f_QE;
MCNUM mccsample0_gamma;
MCNUM mccsample0_target_x;
MCNUM mccsample0_target_y;
MCNUM mccsample0_target_z;
MCNUM mccsample0_focus_xw;
MCNUM mccsample0_focus_yh;
MCNUM mccsample0_focus_aw;
MCNUM mccsample0_focus_ah;
MCNUM mccsample0_xwidth;
MCNUM mccsample0_yheight;
MCNUM mccsample0_zthick;
MCNUM mccsample0_rad_sphere;
MCNUM mccsample0_sig_a;
MCNUM mccsample0_sig_i;
MCNUM mccsample0_V0;
int mccsample0_target_index;
MCNUM mccsample0_multiples;

/* Definition parameters for component 'guide_field_2' [20]. */
#define mccguide_field_2_fieldFunction const_magnetic_field
/* Setting parameters for component 'guide_field_2' [20]. */
MCNUM mccguide_field_2_xwidth;
MCNUM mccguide_field_2_yheight;
MCNUM mccguide_field_2_zdepth;
MCNUM mccguide_field_2_Bx;
MCNUM mccguide_field_2_By;
MCNUM mccguide_field_2_Bz;
char mccguide_field_2_filename[16384];

/* Setting parameters for component 'guide_field_2_end' [21]. */
MCNUM mccguide_field_2_end_xmin;
MCNUM mccguide_field_2_end_xmax;
MCNUM mccguide_field_2_end_ymin;
MCNUM mccguide_field_2_end_ymax;
MCNUM mccguide_field_2_end_radius;
MCNUM mccguide_field_2_end_xwidth;
MCNUM mccguide_field_2_end_yheight;

/* Definition parameters for component 'pl_mon40' [22]. */
#define mccpl_mon40_xwidth 0.1
#define mccpl_mon40_yheight 0.1
#define mccpl_mon40_nL 20
#define mccpl_mon40_restore_neutron 0
#define mccpl_mon40_nowritefile 0
/* Setting parameters for component 'pl_mon40' [22]. */
char mccpl_mon40_filename[16384];
MCNUM mccpl_mon40_mx;
MCNUM mccpl_mon40_my;
MCNUM mccpl_mon40_mz;
MCNUM mccpl_mon40_Lmin;
MCNUM mccpl_mon40_Lmax;

/* Definition parameters for component 'flipper3' [23]. */
#define mccflipper3_fieldFunction const_magnetic_field
/* Setting parameters for component 'flipper3' [23]. */
MCNUM mccflipper3_xwidth;
MCNUM mccflipper3_yheight;
MCNUM mccflipper3_zdepth;
MCNUM mccflipper3_Bx;
MCNUM mccflipper3_By;
MCNUM mccflipper3_Bz;
char mccflipper3_filename[16384];

/* Definition parameters for component 'pl_mon4' [24]. */
#define mccpl_mon4_xwidth 0.1
#define mccpl_mon4_yheight 0.1
#define mccpl_mon4_nL 20
#define mccpl_mon4_restore_neutron 1
#define mccpl_mon4_nowritefile 0
/* Setting parameters for component 'pl_mon4' [24]. */
char mccpl_mon4_filename[16384];
MCNUM mccpl_mon4_mx;
MCNUM mccpl_mon4_my;
MCNUM mccpl_mon4_mz;
MCNUM mccpl_mon4_Lmin;
MCNUM mccpl_mon4_Lmax;

/* Setting parameters for component 'psd_ba' [26]. */
int mccpsd_ba_nx;
int mccpsd_ba_ny;
char mccpsd_ba_filename[16384];
MCNUM mccpsd_ba_xmin;
MCNUM mccpsd_ba_xmax;
MCNUM mccpsd_ba_ymin;
MCNUM mccpsd_ba_ymax;
MCNUM mccpsd_ba_xwidth;
MCNUM mccpsd_ba_yheight;
MCNUM mccpsd_ba_restore_neutron;

/* Definition parameters for component 'pl_mon_ba' [27]. */
#define mccpl_mon_ba_xwidth 0.1
#define mccpl_mon_ba_yheight 0.01
#define mccpl_mon_ba_nL 20
#define mccpl_mon_ba_restore_neutron 1
#define mccpl_mon_ba_nowritefile 0
/* Setting parameters for component 'pl_mon_ba' [27]. */
char mccpl_mon_ba_filename[16384];
MCNUM mccpl_mon_ba_mx;
MCNUM mccpl_mon_ba_my;
MCNUM mccpl_mon_ba_mz;
MCNUM mccpl_mon_ba_Lmin;
MCNUM mccpl_mon_ba_Lmax;

/* Definition parameters for component 'analyzer' [28]. */
#define mccanalyzer_rUpFunc StdReflecFunc
#define mccanalyzer_rDownFunc rUpFunc
#define mccanalyzer_rUpPar { 0.99 , 0.0485 , 0.9 , 1 , 0.0002 }
#define mccanalyzer_rDownPar { 0.99 , 0.11 , 0.9 , 2 , 0.0185 }
#define mccanalyzer_useTables 0
#define mccanalyzer_zwidth 0.1
#define mccanalyzer_yheight 0.3
/* Setting parameters for component 'analyzer' [28]. */
MCNUM mccanalyzer_p_reflect;

/* Setting parameters for component 'detector_slit' [29]. */
MCNUM mccdetector_slit_xmin;
MCNUM mccdetector_slit_xmax;
MCNUM mccdetector_slit_ymin;
MCNUM mccdetector_slit_ymax;
MCNUM mccdetector_slit_radius;
MCNUM mccdetector_slit_xwidth;
MCNUM mccdetector_slit_yheight;

/* Definition parameters for component 'pl_mon5' [30]. */
#define mccpl_mon5_xwidth 0.1
#define mccpl_mon5_yheight 0.1
#define mccpl_mon5_nL 20
#define mccpl_mon5_restore_neutron 1
#define mccpl_mon5_nowritefile 0
/* Setting parameters for component 'pl_mon5' [30]. */
char mccpl_mon5_filename[16384];
MCNUM mccpl_mon5_mx;
MCNUM mccpl_mon5_my;
MCNUM mccpl_mon5_mz;
MCNUM mccpl_mon5_Lmin;
MCNUM mccpl_mon5_Lmax;

/* Setting parameters for component 'detector' [31]. */
int mccdetector_nx;
int mccdetector_ny;
char mccdetector_filename[16384];
MCNUM mccdetector_xmin;
MCNUM mccdetector_xmax;
MCNUM mccdetector_ymin;
MCNUM mccdetector_ymax;
MCNUM mccdetector_xwidth;
MCNUM mccdetector_yheight;
MCNUM mccdetector_restore_neutron;

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
#line 8896 "./SE_example2.c"
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
#line 8933 "./SE_example2.c"
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

/* User declarations for component 'lam_mon' [3]. */
#define mccompcurname  lam_mon
#define mccompcurtype  L_monitor
#define mccompcurindex 3
#define nL mcclam_mon_nL
#define L_N mcclam_mon_L_N
#define L_p mcclam_mon_L_p
#define L_p2 mcclam_mon_L_p2
#define filename mcclam_mon_filename
#define xmin mcclam_mon_xmin
#define xmax mcclam_mon_xmax
#define ymin mcclam_mon_ymin
#define ymax mcclam_mon_ymax
#define xwidth mcclam_mon_xwidth
#define yheight mcclam_mon_yheight
#define Lmin mcclam_mon_Lmin
#define Lmax mcclam_mon_Lmax
#define restore_neutron mcclam_mon_restore_neutron
#define nowritefile mcclam_mon_nowritefile
#line 57 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/L_monitor.comp"
double L_N[nL];
double L_p[nL], L_p2[nL];
#line 8976 "./SE_example2.c"
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

/* User declarations for component 'polariser_pt' [4]. */
#define mccompcurname  polariser_pt
#define mccompcurtype  Arm
#define mccompcurindex 4
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'hdiv_mon' [5]. */
#define mccompcurname  hdiv_mon
#define mccompcurtype  Hdiv_monitor
#define mccompcurindex 5
#define nh mcchdiv_mon_nh
#define Div_N mcchdiv_mon_Div_N
#define Div_p mcchdiv_mon_Div_p
#define Div_p2 mcchdiv_mon_Div_p2
#define filename mcchdiv_mon_filename
#define xmin mcchdiv_mon_xmin
#define xmax mcchdiv_mon_xmax
#define ymin mcchdiv_mon_ymin
#define ymax mcchdiv_mon_ymax
#define xwidth mcchdiv_mon_xwidth
#define yheight mcchdiv_mon_yheight
#define h_maxdiv mcchdiv_mon_h_maxdiv
#define restore_neutron mcchdiv_mon_restore_neutron
#define nowritefile mcchdiv_mon_nowritefile
#line 57 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/Hdiv_monitor.comp"
double Div_N[nh];
double Div_p[nh];
double Div_p2[nh];
#line 9026 "./SE_example2.c"
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

/* User declarations for component 'polariser2' [6]. */
#define mccompcurname  polariser2
#define mccompcurtype  Pol_mirror
#define mccompcurindex 6
#define rUpFunc mccpolariser2_rUpFunc
#define rDownFunc mccpolariser2_rDownFunc
#define rUpPar mccpolariser2_rUpPar
#define rDownPar mccpolariser2_rDownPar
#define useTables mccpolariser2_useTables
#define zwidth mccpolariser2_zwidth
#define yheight mccpolariser2_yheight
#define rUpParPtr mccpolariser2_rUpParPtr
#define rDownParPtr mccpolariser2_rDownParPtr
#define p_reflect mccpolariser2_p_reflect
#line 95 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Pol_mirror.comp"
#if (useTables)
  t_Table *rUpParPtr   = 0;
  t_Table *rDownParPtr = 0;
#else
  double rUpParPtr[]   = (double[])rUpPar;
  double rDownParPtr[] = (double[])rDownPar;
#endif
#line 9067 "./SE_example2.c"
#undef p_reflect
#undef rDownParPtr
#undef rUpParPtr
#undef yheight
#undef zwidth
#undef useTables
#undef rDownPar
#undef rUpPar
#undef rDownFunc
#undef rUpFunc
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'psd_ap' [7]. */
#define mccompcurname  psd_ap
#define mccompcurtype  PSD_monitor
#define mccompcurindex 7
#define PSD_N mccpsd_ap_PSD_N
#define PSD_p mccpsd_ap_PSD_p
#define PSD_p2 mccpsd_ap_PSD_p2
#define nx mccpsd_ap_nx
#define ny mccpsd_ap_ny
#define filename mccpsd_ap_filename
#define xmin mccpsd_ap_xmin
#define xmax mccpsd_ap_xmax
#define ymin mccpsd_ap_ymin
#define ymax mccpsd_ap_ymax
#define xwidth mccpsd_ap_xwidth
#define yheight mccpsd_ap_yheight
#define restore_neutron mccpsd_ap_restore_neutron
#line 62 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
  DArray2d PSD_N;
  DArray2d PSD_p;
  DArray2d PSD_p2;
#line 9103 "./SE_example2.c"
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

/* User declarations for component 'pl_mon0' [8]. */
#define mccompcurname  pl_mon0
#define mccompcurtype  MeanPolLambda_monitor
#define mccompcurindex 8
#define xwidth mccpl_mon0_xwidth
#define yheight mccpl_mon0_yheight
#define nL mccpl_mon0_nL
#define restore_neutron mccpl_mon0_restore_neutron
#define nowritefile mccpl_mon0_nowritefile
#define PolL_N mccpl_mon0_PolL_N
#define PolL_p mccpl_mon0_PolL_p
#define PolL_p2 mccpl_mon0_PolL_p2
#define HelpArray mccpl_mon0_HelpArray
#define filename mccpl_mon0_filename
#define mx mccpl_mon0_mx
#define my mccpl_mon0_my
#define mz mccpl_mon0_mz
#define Lmin mccpl_mon0_Lmin
#define Lmax mccpl_mon0_Lmax
#line 59 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/MeanPolLambda_monitor.comp"
double PolL_N[nL];
double PolL_p[nL];
double PolL_p2[nL];
double HelpArray[nL];
#line 9145 "./SE_example2.c"
#undef Lmax
#undef Lmin
#undef mz
#undef my
#undef mx
#undef filename
#undef HelpArray
#undef PolL_p2
#undef PolL_p
#undef PolL_N
#undef nowritefile
#undef restore_neutron
#undef nL
#undef yheight
#undef xwidth
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'flipper1' [9]. */
#define mccompcurname  flipper1
#define mccompcurtype  Pol_FieldBox
#define mccompcurindex 9
#define fieldFunction mccflipper1_fieldFunction
#define field_parameters mccflipper1_field_parameters
#define const_field_parameters mccflipper1_const_field_parameters
#define xwidth mccflipper1_xwidth
#define yheight mccflipper1_yheight
#define zdepth mccflipper1_zdepth
#define Bx mccflipper1_Bx
#define By mccflipper1_By
#define Bz mccflipper1_Bz
#define filename mccflipper1_filename
#line 50 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Pol_FieldBox.comp"
  double const_field_parameters[3];
  void *field_parameters;
#line 9182 "./SE_example2.c"
#undef filename
#undef Bz
#undef By
#undef Bx
#undef zdepth
#undef yheight
#undef xwidth
#undef const_field_parameters
#undef field_parameters
#undef fieldFunction
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'pl_mon1' [10]. */
#define mccompcurname  pl_mon1
#define mccompcurtype  MeanPolLambda_monitor
#define mccompcurindex 10
#define xwidth mccpl_mon1_xwidth
#define yheight mccpl_mon1_yheight
#define nL mccpl_mon1_nL
#define restore_neutron mccpl_mon1_restore_neutron
#define nowritefile mccpl_mon1_nowritefile
#define PolL_N mccpl_mon1_PolL_N
#define PolL_p mccpl_mon1_PolL_p
#define PolL_p2 mccpl_mon1_PolL_p2
#define HelpArray mccpl_mon1_HelpArray
#define filename mccpl_mon1_filename
#define mx mccpl_mon1_mx
#define my mccpl_mon1_my
#define mz mccpl_mon1_mz
#define Lmin mccpl_mon1_Lmin
#define Lmax mccpl_mon1_Lmax
#line 59 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/MeanPolLambda_monitor.comp"
double PolL_N[nL];
double PolL_p[nL];
double PolL_p2[nL];
double HelpArray[nL];
#line 9221 "./SE_example2.c"
#undef Lmax
#undef Lmin
#undef mz
#undef my
#undef mx
#undef filename
#undef HelpArray
#undef PolL_p2
#undef PolL_p
#undef PolL_N
#undef nowritefile
#undef restore_neutron
#undef nL
#undef yheight
#undef xwidth
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'guide_field_1' [11]. */
#define mccompcurname  guide_field_1
#define mccompcurtype  Pol_FieldBox
#define mccompcurindex 11
#define fieldFunction mccguide_field_1_fieldFunction
#define field_parameters mccguide_field_1_field_parameters
#define const_field_parameters mccguide_field_1_const_field_parameters
#define xwidth mccguide_field_1_xwidth
#define yheight mccguide_field_1_yheight
#define zdepth mccguide_field_1_zdepth
#define Bx mccguide_field_1_Bx
#define By mccguide_field_1_By
#define Bz mccguide_field_1_Bz
#define filename mccguide_field_1_filename
#line 50 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Pol_FieldBox.comp"
  double const_field_parameters[3];
  void *field_parameters;
#line 9258 "./SE_example2.c"
#undef filename
#undef Bz
#undef By
#undef Bx
#undef zdepth
#undef yheight
#undef xwidth
#undef const_field_parameters
#undef field_parameters
#undef fieldFunction
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'guide_field_1_end' [12]. */
#define mccompcurname  guide_field_1_end
#define mccompcurtype  Slit
#define mccompcurindex 12
#define xmin mccguide_field_1_end_xmin
#define xmax mccguide_field_1_end_xmax
#define ymin mccguide_field_1_end_ymin
#define ymax mccguide_field_1_end_ymax
#define radius mccguide_field_1_end_radius
#define xwidth mccguide_field_1_end_xwidth
#define yheight mccguide_field_1_end_yheight
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

/* User declarations for component 'pl_mon2' [13]. */
#define mccompcurname  pl_mon2
#define mccompcurtype  MeanPolLambda_monitor
#define mccompcurindex 13
#define xwidth mccpl_mon2_xwidth
#define yheight mccpl_mon2_yheight
#define nL mccpl_mon2_nL
#define restore_neutron mccpl_mon2_restore_neutron
#define nowritefile mccpl_mon2_nowritefile
#define PolL_N mccpl_mon2_PolL_N
#define PolL_p mccpl_mon2_PolL_p
#define PolL_p2 mccpl_mon2_PolL_p2
#define HelpArray mccpl_mon2_HelpArray
#define filename mccpl_mon2_filename
#define mx mccpl_mon2_mx
#define my mccpl_mon2_my
#define mz mccpl_mon2_mz
#define Lmin mccpl_mon2_Lmin
#define Lmax mccpl_mon2_Lmax
#line 59 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/MeanPolLambda_monitor.comp"
double PolL_N[nL];
double PolL_p[nL];
double PolL_p2[nL];
double HelpArray[nL];
#line 9319 "./SE_example2.c"
#undef Lmax
#undef Lmin
#undef mz
#undef my
#undef mx
#undef filename
#undef HelpArray
#undef PolL_p2
#undef PolL_p
#undef PolL_N
#undef nowritefile
#undef restore_neutron
#undef nL
#undef yheight
#undef xwidth
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'flipper2_1' [14]. */
#define mccompcurname  flipper2_1
#define mccompcurtype  Pol_FieldBox
#define mccompcurindex 14
#define fieldFunction mccflipper2_1_fieldFunction
#define field_parameters mccflipper2_1_field_parameters
#define const_field_parameters mccflipper2_1_const_field_parameters
#define xwidth mccflipper2_1_xwidth
#define yheight mccflipper2_1_yheight
#define zdepth mccflipper2_1_zdepth
#define Bx mccflipper2_1_Bx
#define By mccflipper2_1_By
#define Bz mccflipper2_1_Bz
#define filename mccflipper2_1_filename
#line 50 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Pol_FieldBox.comp"
  double const_field_parameters[3];
  void *field_parameters;
#line 9356 "./SE_example2.c"
#undef filename
#undef Bz
#undef By
#undef Bx
#undef zdepth
#undef yheight
#undef xwidth
#undef const_field_parameters
#undef field_parameters
#undef fieldFunction
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'flipper2_2' [15]. */
#define mccompcurname  flipper2_2
#define mccompcurtype  Pol_FieldBox
#define mccompcurindex 15
#define fieldFunction mccflipper2_2_fieldFunction
#define field_parameters mccflipper2_2_field_parameters
#define const_field_parameters mccflipper2_2_const_field_parameters
#define xwidth mccflipper2_2_xwidth
#define yheight mccflipper2_2_yheight
#define zdepth mccflipper2_2_zdepth
#define Bx mccflipper2_2_Bx
#define By mccflipper2_2_By
#define Bz mccflipper2_2_Bz
#define filename mccflipper2_2_filename
#line 50 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Pol_FieldBox.comp"
  double const_field_parameters[3];
  void *field_parameters;
#line 9388 "./SE_example2.c"
#undef filename
#undef Bz
#undef By
#undef Bx
#undef zdepth
#undef yheight
#undef xwidth
#undef const_field_parameters
#undef field_parameters
#undef fieldFunction
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'flipper_2_2_end' [16]. */
#define mccompcurname  flipper_2_2_end
#define mccompcurtype  Slit
#define mccompcurindex 16
#define xmin mccflipper_2_2_end_xmin
#define xmax mccflipper_2_2_end_xmax
#define ymin mccflipper_2_2_end_ymin
#define ymax mccflipper_2_2_end_ymax
#define radius mccflipper_2_2_end_radius
#define xwidth mccflipper_2_2_end_xwidth
#define yheight mccflipper_2_2_end_yheight
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

/* User declarations for component 'sample_mnt' [17]. */
#define mccompcurname  sample_mnt
#define mccompcurtype  Arm
#define mccompcurindex 17
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'pl_mon3' [18]. */
#define mccompcurname  pl_mon3
#define mccompcurtype  MeanPolLambda_monitor
#define mccompcurindex 18
#define xwidth mccpl_mon3_xwidth
#define yheight mccpl_mon3_yheight
#define nL mccpl_mon3_nL
#define restore_neutron mccpl_mon3_restore_neutron
#define nowritefile mccpl_mon3_nowritefile
#define PolL_N mccpl_mon3_PolL_N
#define PolL_p mccpl_mon3_PolL_p
#define PolL_p2 mccpl_mon3_PolL_p2
#define HelpArray mccpl_mon3_HelpArray
#define filename mccpl_mon3_filename
#define mx mccpl_mon3_mx
#define my mccpl_mon3_my
#define mz mccpl_mon3_mz
#define Lmin mccpl_mon3_Lmin
#define Lmax mccpl_mon3_Lmax
#line 59 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/MeanPolLambda_monitor.comp"
double PolL_N[nL];
double PolL_p[nL];
double PolL_p2[nL];
double HelpArray[nL];
#line 9457 "./SE_example2.c"
#undef Lmax
#undef Lmin
#undef mz
#undef my
#undef mx
#undef filename
#undef HelpArray
#undef PolL_p2
#undef PolL_p
#undef PolL_N
#undef nowritefile
#undef restore_neutron
#undef nL
#undef yheight
#undef xwidth
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'sample0' [19]. */
#define mccompcurname  sample0
#define mccompcurtype  V_sample
#define mccompcurindex 19
#define VarsV mccsample0_VarsV
#define radius mccsample0_radius
#define thickness mccsample0_thickness
#define zdepth mccsample0_zdepth
#define Vc mccsample0_Vc
#define sigma_abs mccsample0_sigma_abs
#define sigma_inc mccsample0_sigma_inc
#define radius_i mccsample0_radius_i
#define radius_o mccsample0_radius_o
#define h mccsample0_h
#define focus_r mccsample0_focus_r
#define pack mccsample0_pack
#define frac mccsample0_frac
#define f_QE mccsample0_f_QE
#define gamma mccsample0_gamma
#define target_x mccsample0_target_x
#define target_y mccsample0_target_y
#define target_z mccsample0_target_z
#define focus_xw mccsample0_focus_xw
#define focus_yh mccsample0_focus_yh
#define focus_aw mccsample0_focus_aw
#define focus_ah mccsample0_focus_ah
#define xwidth mccsample0_xwidth
#define yheight mccsample0_yheight
#define zthick mccsample0_zthick
#define rad_sphere mccsample0_rad_sphere
#define sig_a mccsample0_sig_a
#define sig_i mccsample0_sig_i
#define V0 mccsample0_V0
#define target_index mccsample0_target_index
#define multiples mccsample0_multiples
#line 117 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../obsolete/V_sample.comp"
  struct StructVarsV VarsV;
#line 9514 "./SE_example2.c"
#undef multiples
#undef target_index
#undef V0
#undef sig_i
#undef sig_a
#undef rad_sphere
#undef zthick
#undef yheight
#undef xwidth
#undef focus_ah
#undef focus_aw
#undef focus_yh
#undef focus_xw
#undef target_z
#undef target_y
#undef target_x
#undef gamma
#undef f_QE
#undef frac
#undef pack
#undef focus_r
#undef h
#undef radius_o
#undef radius_i
#undef sigma_inc
#undef sigma_abs
#undef Vc
#undef zdepth
#undef thickness
#undef radius
#undef VarsV
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'guide_field_2' [20]. */
#define mccompcurname  guide_field_2
#define mccompcurtype  Pol_FieldBox
#define mccompcurindex 20
#define fieldFunction mccguide_field_2_fieldFunction
#define field_parameters mccguide_field_2_field_parameters
#define const_field_parameters mccguide_field_2_const_field_parameters
#define xwidth mccguide_field_2_xwidth
#define yheight mccguide_field_2_yheight
#define zdepth mccguide_field_2_zdepth
#define Bx mccguide_field_2_Bx
#define By mccguide_field_2_By
#define Bz mccguide_field_2_Bz
#define filename mccguide_field_2_filename
#line 50 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Pol_FieldBox.comp"
  double const_field_parameters[3];
  void *field_parameters;
#line 9567 "./SE_example2.c"
#undef filename
#undef Bz
#undef By
#undef Bx
#undef zdepth
#undef yheight
#undef xwidth
#undef const_field_parameters
#undef field_parameters
#undef fieldFunction
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'guide_field_2_end' [21]. */
#define mccompcurname  guide_field_2_end
#define mccompcurtype  Slit
#define mccompcurindex 21
#define xmin mccguide_field_2_end_xmin
#define xmax mccguide_field_2_end_xmax
#define ymin mccguide_field_2_end_ymin
#define ymax mccguide_field_2_end_ymax
#define radius mccguide_field_2_end_radius
#define xwidth mccguide_field_2_end_xwidth
#define yheight mccguide_field_2_end_yheight
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

/* User declarations for component 'pl_mon40' [22]. */
#define mccompcurname  pl_mon40
#define mccompcurtype  MeanPolLambda_monitor
#define mccompcurindex 22
#define xwidth mccpl_mon40_xwidth
#define yheight mccpl_mon40_yheight
#define nL mccpl_mon40_nL
#define restore_neutron mccpl_mon40_restore_neutron
#define nowritefile mccpl_mon40_nowritefile
#define PolL_N mccpl_mon40_PolL_N
#define PolL_p mccpl_mon40_PolL_p
#define PolL_p2 mccpl_mon40_PolL_p2
#define HelpArray mccpl_mon40_HelpArray
#define filename mccpl_mon40_filename
#define mx mccpl_mon40_mx
#define my mccpl_mon40_my
#define mz mccpl_mon40_mz
#define Lmin mccpl_mon40_Lmin
#define Lmax mccpl_mon40_Lmax
#line 59 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/MeanPolLambda_monitor.comp"
double PolL_N[nL];
double PolL_p[nL];
double PolL_p2[nL];
double HelpArray[nL];
#line 9628 "./SE_example2.c"
#undef Lmax
#undef Lmin
#undef mz
#undef my
#undef mx
#undef filename
#undef HelpArray
#undef PolL_p2
#undef PolL_p
#undef PolL_N
#undef nowritefile
#undef restore_neutron
#undef nL
#undef yheight
#undef xwidth
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'flipper3' [23]. */
#define mccompcurname  flipper3
#define mccompcurtype  Pol_FieldBox
#define mccompcurindex 23
#define fieldFunction mccflipper3_fieldFunction
#define field_parameters mccflipper3_field_parameters
#define const_field_parameters mccflipper3_const_field_parameters
#define xwidth mccflipper3_xwidth
#define yheight mccflipper3_yheight
#define zdepth mccflipper3_zdepth
#define Bx mccflipper3_Bx
#define By mccflipper3_By
#define Bz mccflipper3_Bz
#define filename mccflipper3_filename
#line 50 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Pol_FieldBox.comp"
  double const_field_parameters[3];
  void *field_parameters;
#line 9665 "./SE_example2.c"
#undef filename
#undef Bz
#undef By
#undef Bx
#undef zdepth
#undef yheight
#undef xwidth
#undef const_field_parameters
#undef field_parameters
#undef fieldFunction
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'pl_mon4' [24]. */
#define mccompcurname  pl_mon4
#define mccompcurtype  MeanPolLambda_monitor
#define mccompcurindex 24
#define xwidth mccpl_mon4_xwidth
#define yheight mccpl_mon4_yheight
#define nL mccpl_mon4_nL
#define restore_neutron mccpl_mon4_restore_neutron
#define nowritefile mccpl_mon4_nowritefile
#define PolL_N mccpl_mon4_PolL_N
#define PolL_p mccpl_mon4_PolL_p
#define PolL_p2 mccpl_mon4_PolL_p2
#define HelpArray mccpl_mon4_HelpArray
#define filename mccpl_mon4_filename
#define mx mccpl_mon4_mx
#define my mccpl_mon4_my
#define mz mccpl_mon4_mz
#define Lmin mccpl_mon4_Lmin
#define Lmax mccpl_mon4_Lmax
#line 59 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/MeanPolLambda_monitor.comp"
double PolL_N[nL];
double PolL_p[nL];
double PolL_p2[nL];
double HelpArray[nL];
#line 9704 "./SE_example2.c"
#undef Lmax
#undef Lmin
#undef mz
#undef my
#undef mx
#undef filename
#undef HelpArray
#undef PolL_p2
#undef PolL_p
#undef PolL_N
#undef nowritefile
#undef restore_neutron
#undef nL
#undef yheight
#undef xwidth
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'analyzer_arm' [25]. */
#define mccompcurname  analyzer_arm
#define mccompcurtype  Arm
#define mccompcurindex 25
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'psd_ba' [26]. */
#define mccompcurname  psd_ba
#define mccompcurtype  PSD_monitor
#define mccompcurindex 26
#define PSD_N mccpsd_ba_PSD_N
#define PSD_p mccpsd_ba_PSD_p
#define PSD_p2 mccpsd_ba_PSD_p2
#define nx mccpsd_ba_nx
#define ny mccpsd_ba_ny
#define filename mccpsd_ba_filename
#define xmin mccpsd_ba_xmin
#define xmax mccpsd_ba_xmax
#define ymin mccpsd_ba_ymin
#define ymax mccpsd_ba_ymax
#define xwidth mccpsd_ba_xwidth
#define yheight mccpsd_ba_yheight
#define restore_neutron mccpsd_ba_restore_neutron
#line 62 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
  DArray2d PSD_N;
  DArray2d PSD_p;
  DArray2d PSD_p2;
#line 9753 "./SE_example2.c"
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

/* User declarations for component 'pl_mon_ba' [27]. */
#define mccompcurname  pl_mon_ba
#define mccompcurtype  MeanPolLambda_monitor
#define mccompcurindex 27
#define xwidth mccpl_mon_ba_xwidth
#define yheight mccpl_mon_ba_yheight
#define nL mccpl_mon_ba_nL
#define restore_neutron mccpl_mon_ba_restore_neutron
#define nowritefile mccpl_mon_ba_nowritefile
#define PolL_N mccpl_mon_ba_PolL_N
#define PolL_p mccpl_mon_ba_PolL_p
#define PolL_p2 mccpl_mon_ba_PolL_p2
#define HelpArray mccpl_mon_ba_HelpArray
#define filename mccpl_mon_ba_filename
#define mx mccpl_mon_ba_mx
#define my mccpl_mon_ba_my
#define mz mccpl_mon_ba_mz
#define Lmin mccpl_mon_ba_Lmin
#define Lmax mccpl_mon_ba_Lmax
#line 59 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/MeanPolLambda_monitor.comp"
double PolL_N[nL];
double PolL_p[nL];
double PolL_p2[nL];
double HelpArray[nL];
#line 9795 "./SE_example2.c"
#undef Lmax
#undef Lmin
#undef mz
#undef my
#undef mx
#undef filename
#undef HelpArray
#undef PolL_p2
#undef PolL_p
#undef PolL_N
#undef nowritefile
#undef restore_neutron
#undef nL
#undef yheight
#undef xwidth
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'analyzer' [28]. */
#define mccompcurname  analyzer
#define mccompcurtype  Pol_mirror
#define mccompcurindex 28
#define rUpFunc mccanalyzer_rUpFunc
#define rDownFunc mccanalyzer_rDownFunc
#define rUpPar mccanalyzer_rUpPar
#define rDownPar mccanalyzer_rDownPar
#define useTables mccanalyzer_useTables
#define zwidth mccanalyzer_zwidth
#define yheight mccanalyzer_yheight
#define rUpParPtr mccanalyzer_rUpParPtr
#define rDownParPtr mccanalyzer_rDownParPtr
#define p_reflect mccanalyzer_p_reflect
#line 95 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Pol_mirror.comp"
#if (useTables)
  t_Table *rUpParPtr   = 0;
  t_Table *rDownParPtr = 0;
#else
  double rUpParPtr[]   = (double[])rUpPar;
  double rDownParPtr[] = (double[])rDownPar;
#endif
#line 9837 "./SE_example2.c"
#undef p_reflect
#undef rDownParPtr
#undef rUpParPtr
#undef yheight
#undef zwidth
#undef useTables
#undef rDownPar
#undef rUpPar
#undef rDownFunc
#undef rUpFunc
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'detector_slit' [29]. */
#define mccompcurname  detector_slit
#define mccompcurtype  Slit
#define mccompcurindex 29
#define xmin mccdetector_slit_xmin
#define xmax mccdetector_slit_xmax
#define ymin mccdetector_slit_ymin
#define ymax mccdetector_slit_ymax
#define radius mccdetector_slit_radius
#define xwidth mccdetector_slit_xwidth
#define yheight mccdetector_slit_yheight
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

/* User declarations for component 'pl_mon5' [30]. */
#define mccompcurname  pl_mon5
#define mccompcurtype  MeanPolLambda_monitor
#define mccompcurindex 30
#define xwidth mccpl_mon5_xwidth
#define yheight mccpl_mon5_yheight
#define nL mccpl_mon5_nL
#define restore_neutron mccpl_mon5_restore_neutron
#define nowritefile mccpl_mon5_nowritefile
#define PolL_N mccpl_mon5_PolL_N
#define PolL_p mccpl_mon5_PolL_p
#define PolL_p2 mccpl_mon5_PolL_p2
#define HelpArray mccpl_mon5_HelpArray
#define filename mccpl_mon5_filename
#define mx mccpl_mon5_mx
#define my mccpl_mon5_my
#define mz mccpl_mon5_mz
#define Lmin mccpl_mon5_Lmin
#define Lmax mccpl_mon5_Lmax
#line 59 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/MeanPolLambda_monitor.comp"
double PolL_N[nL];
double PolL_p[nL];
double PolL_p2[nL];
double HelpArray[nL];
#line 9898 "./SE_example2.c"
#undef Lmax
#undef Lmin
#undef mz
#undef my
#undef mx
#undef filename
#undef HelpArray
#undef PolL_p2
#undef PolL_p
#undef PolL_N
#undef nowritefile
#undef restore_neutron
#undef nL
#undef yheight
#undef xwidth
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'detector' [31]. */
#define mccompcurname  detector
#define mccompcurtype  PSD_monitor
#define mccompcurindex 31
#define PSD_N mccdetector_PSD_N
#define PSD_p mccdetector_PSD_p
#define PSD_p2 mccdetector_PSD_p2
#define nx mccdetector_nx
#define ny mccdetector_ny
#define filename mccdetector_filename
#define xmin mccdetector_xmin
#define xmax mccdetector_xmax
#define ymin mccdetector_ymin
#define ymax mccdetector_ymax
#define xwidth mccdetector_xwidth
#define yheight mccdetector_yheight
#define restore_neutron mccdetector_restore_neutron
#line 62 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
  DArray2d PSD_N;
  DArray2d PSD_p;
  DArray2d PSD_p2;
#line 9939 "./SE_example2.c"
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
Coords mcposalam_mon, mcposrlam_mon;
Rotation mcrotalam_mon, mcrotrlam_mon;
Coords mcposapolariser_pt, mcposrpolariser_pt;
Rotation mcrotapolariser_pt, mcrotrpolariser_pt;
Coords mcposahdiv_mon, mcposrhdiv_mon;
Rotation mcrotahdiv_mon, mcrotrhdiv_mon;
Coords mcposapolariser2, mcposrpolariser2;
Rotation mcrotapolariser2, mcrotrpolariser2;
Coords mcposapsd_ap, mcposrpsd_ap;
Rotation mcrotapsd_ap, mcrotrpsd_ap;
Coords mcposapl_mon0, mcposrpl_mon0;
Rotation mcrotapl_mon0, mcrotrpl_mon0;
Coords mcposaflipper1, mcposrflipper1;
Rotation mcrotaflipper1, mcrotrflipper1;
Coords mcposapl_mon1, mcposrpl_mon1;
Rotation mcrotapl_mon1, mcrotrpl_mon1;
Coords mcposaguide_field_1, mcposrguide_field_1;
Rotation mcrotaguide_field_1, mcrotrguide_field_1;
Coords mcposaguide_field_1_end, mcposrguide_field_1_end;
Rotation mcrotaguide_field_1_end, mcrotrguide_field_1_end;
Coords mcposapl_mon2, mcposrpl_mon2;
Rotation mcrotapl_mon2, mcrotrpl_mon2;
Coords mcposaflipper2_1, mcposrflipper2_1;
Rotation mcrotaflipper2_1, mcrotrflipper2_1;
Coords mcposaflipper2_2, mcposrflipper2_2;
Rotation mcrotaflipper2_2, mcrotrflipper2_2;
Coords mcposaflipper_2_2_end, mcposrflipper_2_2_end;
Rotation mcrotaflipper_2_2_end, mcrotrflipper_2_2_end;
Coords mcposasample_mnt, mcposrsample_mnt;
Rotation mcrotasample_mnt, mcrotrsample_mnt;
Coords mcposapl_mon3, mcposrpl_mon3;
Rotation mcrotapl_mon3, mcrotrpl_mon3;
Coords mcposasample0, mcposrsample0;
Rotation mcrotasample0, mcrotrsample0;
Coords mcposaguide_field_2, mcposrguide_field_2;
Rotation mcrotaguide_field_2, mcrotrguide_field_2;
Coords mcposaguide_field_2_end, mcposrguide_field_2_end;
Rotation mcrotaguide_field_2_end, mcrotrguide_field_2_end;
Coords mcposapl_mon40, mcposrpl_mon40;
Rotation mcrotapl_mon40, mcrotrpl_mon40;
Coords mcposaflipper3, mcposrflipper3;
Rotation mcrotaflipper3, mcrotrflipper3;
Coords mcposapl_mon4, mcposrpl_mon4;
Rotation mcrotapl_mon4, mcrotrpl_mon4;
Coords mcposaanalyzer_arm, mcposranalyzer_arm;
Rotation mcrotaanalyzer_arm, mcrotranalyzer_arm;
Coords mcposapsd_ba, mcposrpsd_ba;
Rotation mcrotapsd_ba, mcrotrpsd_ba;
Coords mcposapl_mon_ba, mcposrpl_mon_ba;
Rotation mcrotapl_mon_ba, mcrotrpl_mon_ba;
Coords mcposaanalyzer, mcposranalyzer;
Rotation mcrotaanalyzer, mcrotranalyzer;
Coords mcposadetector_slit, mcposrdetector_slit;
Rotation mcrotadetector_slit, mcrotrdetector_slit;
Coords mcposapl_mon5, mcposrpl_mon5;
Rotation mcrotapl_mon5, mcrotrpl_mon5;
Coords mcposadetector, mcposrdetector;
Rotation mcrotadetector, mcrotrdetector;

MCNUM mcnx, mcny, mcnz, mcnvx, mcnvy, mcnvz, mcnt, mcnsx, mcnsy, mcnsz, mcnp;

/* end declare */

void mcinit(void) {
#define mccompcurname  SE_example2
#define mccompcurtype  INSTRUMENT
#define mccompcurindex 0
#define mcposaSE_example2 coords_set(0,0,0)
#define POL_ANGLE mcipPOL_ANGLE
#define SAMPLE mcipSAMPLE
#define Bguide mcipBguide
#define Bflip mcipBflip
#define dBz mcipdBz
#define Lam mcipLam
#define dLam mcipdLam
#line 48 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
{
  pol_radius=0.3/(POL_ANGLE*DEG2RAD);
}
#line 10040 "./SE_example2.c"
#undef dLam
#undef Lam
#undef dBz
#undef Bflip
#undef Bguide
#undef SAMPLE
#undef POL_ANGLE
#undef mcposaSE_example2
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
#line 39 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  if("NULL") strncpy(mccOrigin_profile, "NULL" ? "NULL" : "", 16384); else mccOrigin_profile[0]='\0';
#line 39 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccOrigin_percent = 10;
#line 39 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccOrigin_flag_save = 0;
#line 39 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccOrigin_minutes = 0;
#line 10073 "./SE_example2.c"

  SIG_MESSAGE("Origin (Init:Place/Rotate)");
  rot_set_rotation(mcrotaOrigin,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 10080 "./SE_example2.c"
  rot_copy(mcrotrOrigin, mcrotaOrigin);
  mcposaOrigin = coords_set(
#line 64 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    0,
#line 64 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    0,
#line 64 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    0);
#line 10089 "./SE_example2.c"
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
#line 67 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccSource_radius = .012;
#line 67 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccSource_yheight = 0.1 * sin ( mcipPOL_ANGLE * DEG2RAD );
#line 67 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccSource_xwidth = 0.02;
#line 67 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccSource_dist = 1;
#line 53 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccSource_focus_xw = .045;
#line 53 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccSource_focus_yh = .12;
#line 54 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccSource_E0 = 0;
#line 54 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccSource_dE = 0;
#line 67 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccSource_lambda0 = mcipLam;
#line 68 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccSource_dlambda = mcipdLam;
#line 68 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccSource_flux = 1e7;
#line 55 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccSource_gauss = 0;
#line 55 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccSource_target_index = + 1;
#line 10126 "./SE_example2.c"

  SIG_MESSAGE("Source (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 10133 "./SE_example2.c"
  rot_mul(mctr1, mcrotaOrigin, mcrotaSource);
  rot_transpose(mcrotaOrigin, mctr1);
  rot_mul(mcrotaSource, mctr1, mcrotrSource);
  mctc1 = coords_set(
#line 69 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    0,
#line 69 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    0,
#line 69 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    0);
#line 10144 "./SE_example2.c"
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
    /* Component lam_mon. */
  /* Setting parameters for component lam_mon. */
  SIG_MESSAGE("lam_mon (Init:SetPar)");
#line 72 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  if("lam_mon") strncpy(mcclam_mon_filename, "lam_mon" ? "lam_mon" : "", 16384); else mcclam_mon_filename[0]='\0';
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mcclam_mon_xmin = -0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mcclam_mon_xmax = 0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mcclam_mon_ymin = -0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mcclam_mon_ymax = 0.05;
#line 73 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mcclam_mon_xwidth = 0.05;
#line 73 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mcclam_mon_yheight = 0.05;
#line 73 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mcclam_mon_Lmin = 1;
#line 73 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mcclam_mon_Lmax = 10;
#line 72 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mcclam_mon_restore_neutron = 1;
#line 51 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mcclam_mon_nowritefile = 0;
#line 10180 "./SE_example2.c"

  SIG_MESSAGE("lam_mon (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 10187 "./SE_example2.c"
  rot_mul(mctr1, mcrotaSource, mcrotalam_mon);
  rot_transpose(mcrotaSource, mctr1);
  rot_mul(mcrotalam_mon, mctr1, mcrotrlam_mon);
  mctc1 = coords_set(
#line 74 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    0,
#line 74 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    0,
#line 74 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    1);
#line 10198 "./SE_example2.c"
  rot_transpose(mcrotaSource, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposalam_mon = coords_add(mcposaSource, mctc2);
  mctc1 = coords_sub(mcposaSource, mcposalam_mon);
  mcposrlam_mon = rot_apply(mcrotalam_mon, mctc1);
  mcDEBUG_COMPONENT("lam_mon", mcposalam_mon, mcrotalam_mon)
  mccomp_posa[3] = mcposalam_mon;
  mccomp_posr[3] = mcposrlam_mon;
  mcNCounter[3]  = mcPCounter[3] = mcP2Counter[3] = 0;
  mcAbsorbProp[3]= 0;
    /* Component polariser_pt. */
  /* Setting parameters for component polariser_pt. */
  SIG_MESSAGE("polariser_pt (Init:SetPar)");

  SIG_MESSAGE("polariser_pt (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 78 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    (0)*DEG2RAD,
#line 78 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    (0)*DEG2RAD,
#line 78 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    (0)*DEG2RAD);
#line 10221 "./SE_example2.c"
  rot_mul(mctr1, mcrotaSource, mcrotapolariser_pt);
  rot_transpose(mcrotalam_mon, mctr1);
  rot_mul(mcrotapolariser_pt, mctr1, mcrotrpolariser_pt);
  mctc1 = coords_set(
#line 77 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    0,
#line 77 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    0,
#line 77 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    1);
#line 10232 "./SE_example2.c"
  rot_transpose(mcrotaSource, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposapolariser_pt = coords_add(mcposaSource, mctc2);
  mctc1 = coords_sub(mcposalam_mon, mcposapolariser_pt);
  mcposrpolariser_pt = rot_apply(mcrotapolariser_pt, mctc1);
  mcDEBUG_COMPONENT("polariser_pt", mcposapolariser_pt, mcrotapolariser_pt)
  mccomp_posa[4] = mcposapolariser_pt;
  mccomp_posr[4] = mcposrpolariser_pt;
  mcNCounter[4]  = mcPCounter[4] = mcP2Counter[4] = 0;
  mcAbsorbProp[4]= 0;
    /* Component hdiv_mon. */
  /* Setting parameters for component hdiv_mon. */
  SIG_MESSAGE("hdiv_mon (Init:SetPar)");
#line 81 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  if("hdiv_mon") strncpy(mcchdiv_mon_filename, "hdiv_mon" ? "hdiv_mon" : "", 16384); else mcchdiv_mon_filename[0]='\0';
#line 51 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mcchdiv_mon_xmin = -0.05;
#line 51 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mcchdiv_mon_xmax = 0.05;
#line 51 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mcchdiv_mon_ymin = -0.05;
#line 51 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mcchdiv_mon_ymax = 0.05;
#line 82 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mcchdiv_mon_xwidth = 0.05;
#line 82 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mcchdiv_mon_yheight = 0.05;
#line 82 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mcchdiv_mon_h_maxdiv = 1.0;
#line 81 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mcchdiv_mon_restore_neutron = 1;
#line 52 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mcchdiv_mon_nowritefile = 0;
#line 10266 "./SE_example2.c"

  SIG_MESSAGE("hdiv_mon (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 84 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    (0)*DEG2RAD,
#line 84 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    (0)*DEG2RAD,
#line 84 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    (90)*DEG2RAD);
#line 10276 "./SE_example2.c"
  rot_mul(mctr1, mcrotapolariser_pt, mcrotahdiv_mon);
  rot_transpose(mcrotapolariser_pt, mctr1);
  rot_mul(mcrotahdiv_mon, mctr1, mcrotrhdiv_mon);
  mctc1 = coords_set(
#line 83 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    0,
#line 83 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    0,
#line 83 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    0);
#line 10287 "./SE_example2.c"
  rot_transpose(mcrotapolariser_pt, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposahdiv_mon = coords_add(mcposapolariser_pt, mctc2);
  mctc1 = coords_sub(mcposapolariser_pt, mcposahdiv_mon);
  mcposrhdiv_mon = rot_apply(mcrotahdiv_mon, mctc1);
  mcDEBUG_COMPONENT("hdiv_mon", mcposahdiv_mon, mcrotahdiv_mon)
  mccomp_posa[5] = mcposahdiv_mon;
  mccomp_posr[5] = mcposrhdiv_mon;
  mcNCounter[5]  = mcPCounter[5] = mcP2Counter[5] = 0;
  mcAbsorbProp[5]= 0;
    /* Component polariser2. */
  /* Setting parameters for component polariser2. */
  SIG_MESSAGE("polariser2 (Init:SetPar)");
#line 89 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccpolariser2_p_reflect = 0;
#line 10303 "./SE_example2.c"

  SIG_MESSAGE("polariser2 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 91 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    (90 - mcipPOL_ANGLE)*DEG2RAD,
#line 91 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    (90)*DEG2RAD,
#line 91 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    (0)*DEG2RAD);
#line 10313 "./SE_example2.c"
  rot_mul(mctr1, mcrotapolariser_pt, mcrotapolariser2);
  rot_transpose(mcrotahdiv_mon, mctr1);
  rot_mul(mcrotapolariser2, mctr1, mcrotrpolariser2);
  mctc1 = coords_set(
#line 90 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    0,
#line 90 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    0,
#line 90 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    0);
#line 10324 "./SE_example2.c"
  rot_transpose(mcrotapolariser_pt, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposapolariser2 = coords_add(mcposapolariser_pt, mctc2);
  mctc1 = coords_sub(mcposahdiv_mon, mcposapolariser2);
  mcposrpolariser2 = rot_apply(mcrotapolariser2, mctc1);
  mcDEBUG_COMPONENT("polariser2", mcposapolariser2, mcrotapolariser2)
  mccomp_posa[6] = mcposapolariser2;
  mccomp_posr[6] = mcposrpolariser2;
  mcNCounter[6]  = mcPCounter[6] = mcP2Counter[6] = 0;
  mcAbsorbProp[6]= 0;
    /* Component psd_ap. */
  /* Setting parameters for component psd_ap. */
  SIG_MESSAGE("psd_ap (Init:SetPar)");
#line 49 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccpsd_ap_nx = 90;
#line 49 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccpsd_ap_ny = 90;
#line 98 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  if("psd_ap") strncpy(mccpsd_ap_filename, "psd_ap" ? "psd_ap" : "", 16384); else mccpsd_ap_filename[0]='\0';
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccpsd_ap_xmin = -0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccpsd_ap_xmax = 0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccpsd_ap_ymin = -0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccpsd_ap_ymax = 0.05;
#line 98 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccpsd_ap_xwidth = 0.1;
#line 98 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccpsd_ap_yheight = 0.01;
#line 98 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccpsd_ap_restore_neutron = 1;
#line 10358 "./SE_example2.c"

  SIG_MESSAGE("psd_ap (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 100 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    (0)*DEG2RAD,
#line 100 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    (0)*DEG2RAD,
#line 100 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    (0)*DEG2RAD);
#line 10368 "./SE_example2.c"
  rot_mul(mctr1, mcrotaSource, mcrotapsd_ap);
  rot_transpose(mcrotapolariser2, mctr1);
  rot_mul(mcrotapsd_ap, mctr1, mcrotrpsd_ap);
  mctc1 = coords_set(
#line 99 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    0,
#line 99 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    0,
#line 99 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    0.32);
#line 10379 "./SE_example2.c"
  rot_transpose(mcrotapolariser_pt, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposapsd_ap = coords_add(mcposapolariser_pt, mctc2);
  mctc1 = coords_sub(mcposapolariser2, mcposapsd_ap);
  mcposrpsd_ap = rot_apply(mcrotapsd_ap, mctc1);
  mcDEBUG_COMPONENT("psd_ap", mcposapsd_ap, mcrotapsd_ap)
  mccomp_posa[7] = mcposapsd_ap;
  mccomp_posr[7] = mcposrpsd_ap;
  mcNCounter[7]  = mcPCounter[7] = mcP2Counter[7] = 0;
  mcAbsorbProp[7]= 0;
    /* Component pl_mon0. */
  /* Setting parameters for component pl_mon0. */
  SIG_MESSAGE("pl_mon0 (Init:SetPar)");
#line 103 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  if("pol_lam0") strncpy(mccpl_mon0_filename, "pol_lam0" ? "pol_lam0" : "", 16384); else mccpl_mon0_filename[0]='\0';
#line 104 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccpl_mon0_mx = 0;
#line 104 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccpl_mon0_my = 0;
#line 104 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccpl_mon0_mz = 1;
#line 104 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccpl_mon0_Lmin = 1;
#line 105 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccpl_mon0_Lmax = 10;
#line 10405 "./SE_example2.c"

  SIG_MESSAGE("pl_mon0 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 107 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    (0)*DEG2RAD,
#line 107 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    (0)*DEG2RAD,
#line 107 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    (0)*DEG2RAD);
#line 10415 "./SE_example2.c"
  rot_mul(mctr1, mcrotaSource, mcrotapl_mon0);
  rot_transpose(mcrotapsd_ap, mctr1);
  rot_mul(mcrotapl_mon0, mctr1, mcrotrpl_mon0);
  mctc1 = coords_set(
#line 106 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    0,
#line 106 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    0,
#line 106 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    0.32);
#line 10426 "./SE_example2.c"
  rot_transpose(mcrotapolariser_pt, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposapl_mon0 = coords_add(mcposapolariser_pt, mctc2);
  mctc1 = coords_sub(mcposapsd_ap, mcposapl_mon0);
  mcposrpl_mon0 = rot_apply(mcrotapl_mon0, mctc1);
  mcDEBUG_COMPONENT("pl_mon0", mcposapl_mon0, mcrotapl_mon0)
  mccomp_posa[8] = mcposapl_mon0;
  mccomp_posr[8] = mcposrpl_mon0;
  mcNCounter[8]  = mcPCounter[8] = mcP2Counter[8] = 0;
  mcAbsorbProp[8]= 0;
    /* Component flipper1. */
  /* Setting parameters for component flipper1. */
  SIG_MESSAGE("flipper1 (Init:SetPar)");
#line 112 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccflipper1_xwidth = 0.1;
#line 112 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccflipper1_yheight = 0.1;
#line 112 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccflipper1_zdepth = 0.02;
#line 113 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccflipper1_Bx = 0;
#line 113 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccflipper1_By = mcipBflip;
#line 113 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccflipper1_Bz = mcipBflip;
#line 38 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  if("") strncpy(mccflipper1_filename, "" ? "" : "", 16384); else mccflipper1_filename[0]='\0';
#line 10454 "./SE_example2.c"

  SIG_MESSAGE("flipper1 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 115 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    (0)*DEG2RAD,
#line 115 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    (0)*DEG2RAD,
#line 115 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    (0)*DEG2RAD);
#line 10464 "./SE_example2.c"
  rot_mul(mctr1, mcrotaSource, mcrotaflipper1);
  rot_transpose(mcrotapl_mon0, mctr1);
  rot_mul(mcrotaflipper1, mctr1, mcrotrflipper1);
  mctc1 = coords_set(
#line 114 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    0,
#line 114 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    0,
#line 114 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    0.32);
#line 10475 "./SE_example2.c"
  rot_transpose(mcrotapolariser_pt, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaflipper1 = coords_add(mcposapolariser_pt, mctc2);
  mctc1 = coords_sub(mcposapl_mon0, mcposaflipper1);
  mcposrflipper1 = rot_apply(mcrotaflipper1, mctc1);
  mcDEBUG_COMPONENT("flipper1", mcposaflipper1, mcrotaflipper1)
  mccomp_posa[9] = mcposaflipper1;
  mccomp_posr[9] = mcposrflipper1;
  mcNCounter[9]  = mcPCounter[9] = mcP2Counter[9] = 0;
  mcAbsorbProp[9]= 0;
    /* Component pl_mon1. */
  /* Setting parameters for component pl_mon1. */
  SIG_MESSAGE("pl_mon1 (Init:SetPar)");
#line 118 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  if("pol_lam1") strncpy(mccpl_mon1_filename, "pol_lam1" ? "pol_lam1" : "", 16384); else mccpl_mon1_filename[0]='\0';
#line 119 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccpl_mon1_mx = 0;
#line 119 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccpl_mon1_my = 1;
#line 119 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccpl_mon1_mz = 0;
#line 119 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccpl_mon1_Lmin = 1;
#line 120 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccpl_mon1_Lmax = 10;
#line 10501 "./SE_example2.c"

  SIG_MESSAGE("pl_mon1 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 122 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    (0)*DEG2RAD,
#line 122 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    (0)*DEG2RAD,
#line 122 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    (0)*DEG2RAD);
#line 10511 "./SE_example2.c"
  rot_mul(mctr1, mcrotaSource, mcrotapl_mon1);
  rot_transpose(mcrotaflipper1, mctr1);
  rot_mul(mcrotapl_mon1, mctr1, mcrotrpl_mon1);
  mctc1 = coords_set(
#line 121 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    0,
#line 121 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    0,
#line 121 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    0.021);
#line 10522 "./SE_example2.c"
  rot_transpose(mcrotaflipper1, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposapl_mon1 = coords_add(mcposaflipper1, mctc2);
  mctc1 = coords_sub(mcposaflipper1, mcposapl_mon1);
  mcposrpl_mon1 = rot_apply(mcrotapl_mon1, mctc1);
  mcDEBUG_COMPONENT("pl_mon1", mcposapl_mon1, mcrotapl_mon1)
  mccomp_posa[10] = mcposapl_mon1;
  mccomp_posr[10] = mcposrpl_mon1;
  mcNCounter[10]  = mcPCounter[10] = mcP2Counter[10] = 0;
  mcAbsorbProp[10]= 0;
    /* Component guide_field_1. */
  /* Setting parameters for component guide_field_1. */
  SIG_MESSAGE("guide_field_1 (Init:SetPar)");
#line 127 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccguide_field_1_xwidth = 0.1;
#line 127 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccguide_field_1_yheight = 0.1;
#line 127 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccguide_field_1_zdepth = 1;
#line 127 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccguide_field_1_Bx = 0;
#line 127 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccguide_field_1_By = 0;
#line 127 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccguide_field_1_Bz = mcipBguide + mcipdBz;
#line 38 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  if("") strncpy(mccguide_field_1_filename, "" ? "" : "", 16384); else mccguide_field_1_filename[0]='\0';
#line 10550 "./SE_example2.c"

  SIG_MESSAGE("guide_field_1 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 10557 "./SE_example2.c"
  rot_mul(mctr1, mcrotaflipper1, mcrotaguide_field_1);
  rot_transpose(mcrotapl_mon1, mctr1);
  rot_mul(mcrotaguide_field_1, mctr1, mcrotrguide_field_1);
  mctc1 = coords_set(
#line 128 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    0,
#line 128 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    0,
#line 128 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    0.021);
#line 10568 "./SE_example2.c"
  rot_transpose(mcrotaflipper1, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaguide_field_1 = coords_add(mcposaflipper1, mctc2);
  mctc1 = coords_sub(mcposapl_mon1, mcposaguide_field_1);
  mcposrguide_field_1 = rot_apply(mcrotaguide_field_1, mctc1);
  mcDEBUG_COMPONENT("guide_field_1", mcposaguide_field_1, mcrotaguide_field_1)
  mccomp_posa[11] = mcposaguide_field_1;
  mccomp_posr[11] = mcposrguide_field_1;
  mcNCounter[11]  = mcPCounter[11] = mcP2Counter[11] = 0;
  mcAbsorbProp[11]= 0;
    /* Component guide_field_1_end. */
  /* Setting parameters for component guide_field_1_end. */
  SIG_MESSAGE("guide_field_1_end (Init:SetPar)");
#line 46 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccguide_field_1_end_xmin = 0;
#line 46 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccguide_field_1_end_xmax = 0;
#line 46 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccguide_field_1_end_ymin = 0;
#line 46 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccguide_field_1_end_ymax = 0;
#line 46 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccguide_field_1_end_radius = 0;
#line 130 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccguide_field_1_end_xwidth = 0.1;
#line 130 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccguide_field_1_end_yheight = 0.01;
#line 10596 "./SE_example2.c"

  SIG_MESSAGE("guide_field_1_end (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 132 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    (0)*DEG2RAD,
#line 132 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    (0)*DEG2RAD,
#line 132 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    (0)*DEG2RAD);
#line 10606 "./SE_example2.c"
  rot_mul(mctr1, mcrotaSource, mcrotaguide_field_1_end);
  rot_transpose(mcrotaguide_field_1, mctr1);
  rot_mul(mcrotaguide_field_1_end, mctr1, mcrotrguide_field_1_end);
  mctc1 = coords_set(
#line 131 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    0,
#line 131 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    0,
#line 131 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    1.001);
#line 10617 "./SE_example2.c"
  rot_transpose(mcrotaguide_field_1, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaguide_field_1_end = coords_add(mcposaguide_field_1, mctc2);
  mctc1 = coords_sub(mcposaguide_field_1, mcposaguide_field_1_end);
  mcposrguide_field_1_end = rot_apply(mcrotaguide_field_1_end, mctc1);
  mcDEBUG_COMPONENT("guide_field_1_end", mcposaguide_field_1_end, mcrotaguide_field_1_end)
  mccomp_posa[12] = mcposaguide_field_1_end;
  mccomp_posr[12] = mcposrguide_field_1_end;
  mcNCounter[12]  = mcPCounter[12] = mcP2Counter[12] = 0;
  mcAbsorbProp[12]= 0;
    /* Component pl_mon2. */
  /* Setting parameters for component pl_mon2. */
  SIG_MESSAGE("pl_mon2 (Init:SetPar)");
#line 135 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  if("pol_lam2") strncpy(mccpl_mon2_filename, "pol_lam2" ? "pol_lam2" : "", 16384); else mccpl_mon2_filename[0]='\0';
#line 136 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccpl_mon2_mx = 0;
#line 136 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccpl_mon2_my = 1;
#line 136 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccpl_mon2_mz = 0;
#line 136 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccpl_mon2_Lmin = 1;
#line 137 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccpl_mon2_Lmax = 10;
#line 10643 "./SE_example2.c"

  SIG_MESSAGE("pl_mon2 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 139 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    (0)*DEG2RAD,
#line 139 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    (0)*DEG2RAD,
#line 139 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    (0)*DEG2RAD);
#line 10653 "./SE_example2.c"
  rot_mul(mctr1, mcrotaSource, mcrotapl_mon2);
  rot_transpose(mcrotaguide_field_1_end, mctr1);
  rot_mul(mcrotapl_mon2, mctr1, mcrotrpl_mon2);
  mctc1 = coords_set(
#line 138 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    0,
#line 138 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    0,
#line 138 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    1.001);
#line 10664 "./SE_example2.c"
  rot_transpose(mcrotaguide_field_1, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposapl_mon2 = coords_add(mcposaguide_field_1, mctc2);
  mctc1 = coords_sub(mcposaguide_field_1_end, mcposapl_mon2);
  mcposrpl_mon2 = rot_apply(mcrotapl_mon2, mctc1);
  mcDEBUG_COMPONENT("pl_mon2", mcposapl_mon2, mcrotapl_mon2)
  mccomp_posa[13] = mcposapl_mon2;
  mccomp_posr[13] = mcposrpl_mon2;
  mcNCounter[13]  = mcPCounter[13] = mcP2Counter[13] = 0;
  mcAbsorbProp[13]= 0;
    /* Component flipper2_1. */
  /* Setting parameters for component flipper2_1. */
  SIG_MESSAGE("flipper2_1 (Init:SetPar)");
#line 143 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccflipper2_1_xwidth = 0.1;
#line 143 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccflipper2_1_yheight = 0.1;
#line 143 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccflipper2_1_zdepth = 0.02;
#line 144 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccflipper2_1_Bx = 0;
#line 144 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccflipper2_1_By = mcipBflip;
#line 144 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccflipper2_1_Bz = mcipBflip;
#line 38 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  if("") strncpy(mccflipper2_1_filename, "" ? "" : "", 16384); else mccflipper2_1_filename[0]='\0';
#line 10692 "./SE_example2.c"

  SIG_MESSAGE("flipper2_1 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 10699 "./SE_example2.c"
  rot_mul(mctr1, mcrotaguide_field_1, mcrotaflipper2_1);
  rot_transpose(mcrotapl_mon2, mctr1);
  rot_mul(mcrotaflipper2_1, mctr1, mcrotrflipper2_1);
  mctc1 = coords_set(
#line 145 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    0,
#line 145 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    0,
#line 145 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    1.001 + 2e-3);
#line 10710 "./SE_example2.c"
  rot_transpose(mcrotaguide_field_1, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaflipper2_1 = coords_add(mcposaguide_field_1, mctc2);
  mctc1 = coords_sub(mcposapl_mon2, mcposaflipper2_1);
  mcposrflipper2_1 = rot_apply(mcrotaflipper2_1, mctc1);
  mcDEBUG_COMPONENT("flipper2_1", mcposaflipper2_1, mcrotaflipper2_1)
  mccomp_posa[14] = mcposaflipper2_1;
  mccomp_posr[14] = mcposrflipper2_1;
  mcNCounter[14]  = mcPCounter[14] = mcP2Counter[14] = 0;
  mcAbsorbProp[14]= 0;
    /* Component flipper2_2. */
  /* Setting parameters for component flipper2_2. */
  SIG_MESSAGE("flipper2_2 (Init:SetPar)");
#line 148 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccflipper2_2_xwidth = 0.1;
#line 148 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccflipper2_2_yheight = 0.1;
#line 148 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccflipper2_2_zdepth = 0.02;
#line 149 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccflipper2_2_Bx = 0;
#line 149 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccflipper2_2_By = - mcipBflip;
#line 149 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccflipper2_2_Bz = mcipBflip;
#line 38 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  if("") strncpy(mccflipper2_2_filename, "" ? "" : "", 16384); else mccflipper2_2_filename[0]='\0';
#line 10738 "./SE_example2.c"

  SIG_MESSAGE("flipper2_2 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 10745 "./SE_example2.c"
  rot_mul(mctr1, mcrotaflipper2_1, mcrotaflipper2_2);
  rot_transpose(mcrotaflipper2_1, mctr1);
  rot_mul(mcrotaflipper2_2, mctr1, mcrotrflipper2_2);
  mctc1 = coords_set(
#line 150 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    0,
#line 150 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    0,
#line 150 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    0.021);
#line 10756 "./SE_example2.c"
  rot_transpose(mcrotaflipper2_1, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaflipper2_2 = coords_add(mcposaflipper2_1, mctc2);
  mctc1 = coords_sub(mcposaflipper2_1, mcposaflipper2_2);
  mcposrflipper2_2 = rot_apply(mcrotaflipper2_2, mctc1);
  mcDEBUG_COMPONENT("flipper2_2", mcposaflipper2_2, mcrotaflipper2_2)
  mccomp_posa[15] = mcposaflipper2_2;
  mccomp_posr[15] = mcposrflipper2_2;
  mcNCounter[15]  = mcPCounter[15] = mcP2Counter[15] = 0;
  mcAbsorbProp[15]= 0;
    /* Component flipper_2_2_end. */
  /* Setting parameters for component flipper_2_2_end. */
  SIG_MESSAGE("flipper_2_2_end (Init:SetPar)");
#line 46 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccflipper_2_2_end_xmin = 0;
#line 46 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccflipper_2_2_end_xmax = 0;
#line 46 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccflipper_2_2_end_ymin = 0;
#line 46 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccflipper_2_2_end_ymax = 0;
#line 46 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccflipper_2_2_end_radius = 0;
#line 152 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccflipper_2_2_end_xwidth = 0.1;
#line 152 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccflipper_2_2_end_yheight = 0.01;
#line 10784 "./SE_example2.c"

  SIG_MESSAGE("flipper_2_2_end (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 10791 "./SE_example2.c"
  rot_mul(mctr1, mcrotaflipper2_2, mcrotaflipper_2_2_end);
  rot_transpose(mcrotaflipper2_2, mctr1);
  rot_mul(mcrotaflipper_2_2_end, mctr1, mcrotrflipper_2_2_end);
  mctc1 = coords_set(
#line 153 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    0,
#line 153 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    0,
#line 153 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    0.0201);
#line 10802 "./SE_example2.c"
  rot_transpose(mcrotaflipper2_2, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaflipper_2_2_end = coords_add(mcposaflipper2_2, mctc2);
  mctc1 = coords_sub(mcposaflipper2_2, mcposaflipper_2_2_end);
  mcposrflipper_2_2_end = rot_apply(mcrotaflipper_2_2_end, mctc1);
  mcDEBUG_COMPONENT("flipper_2_2_end", mcposaflipper_2_2_end, mcrotaflipper_2_2_end)
  mccomp_posa[16] = mcposaflipper_2_2_end;
  mccomp_posr[16] = mcposrflipper_2_2_end;
  mcNCounter[16]  = mcPCounter[16] = mcP2Counter[16] = 0;
  mcAbsorbProp[16]= 0;
    /* Component sample_mnt. */
  /* Setting parameters for component sample_mnt. */
  SIG_MESSAGE("sample_mnt (Init:SetPar)");

  SIG_MESSAGE("sample_mnt (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 10822 "./SE_example2.c"
  rot_mul(mctr1, mcrotaflipper2_2, mcrotasample_mnt);
  rot_transpose(mcrotaflipper_2_2_end, mctr1);
  rot_mul(mcrotasample_mnt, mctr1, mcrotrsample_mnt);
  mctc1 = coords_set(
#line 157 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    0,
#line 157 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    0,
#line 157 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    0.02 + 0.021);
#line 10833 "./SE_example2.c"
  rot_transpose(mcrotaflipper2_2, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposasample_mnt = coords_add(mcposaflipper2_2, mctc2);
  mctc1 = coords_sub(mcposaflipper_2_2_end, mcposasample_mnt);
  mcposrsample_mnt = rot_apply(mcrotasample_mnt, mctc1);
  mcDEBUG_COMPONENT("sample_mnt", mcposasample_mnt, mcrotasample_mnt)
  mccomp_posa[17] = mcposasample_mnt;
  mccomp_posr[17] = mcposrsample_mnt;
  mcNCounter[17]  = mcPCounter[17] = mcP2Counter[17] = 0;
  mcAbsorbProp[17]= 0;
    /* Component pl_mon3. */
  /* Setting parameters for component pl_mon3. */
  SIG_MESSAGE("pl_mon3 (Init:SetPar)");
#line 160 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  if("pol_lam3") strncpy(mccpl_mon3_filename, "pol_lam3" ? "pol_lam3" : "", 16384); else mccpl_mon3_filename[0]='\0';
#line 161 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccpl_mon3_mx = 0;
#line 161 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccpl_mon3_my = 1;
#line 161 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccpl_mon3_mz = 0;
#line 161 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccpl_mon3_Lmin = 1;
#line 162 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccpl_mon3_Lmax = 10;
#line 10859 "./SE_example2.c"

  SIG_MESSAGE("pl_mon3 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 10866 "./SE_example2.c"
  rot_mul(mctr1, mcrotasample_mnt, mcrotapl_mon3);
  rot_transpose(mcrotasample_mnt, mctr1);
  rot_mul(mcrotapl_mon3, mctr1, mcrotrpl_mon3);
  mctc1 = coords_set(
#line 163 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    0,
#line 163 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    0,
#line 163 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    0);
#line 10877 "./SE_example2.c"
  rot_transpose(mcrotasample_mnt, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposapl_mon3 = coords_add(mcposasample_mnt, mctc2);
  mctc1 = coords_sub(mcposasample_mnt, mcposapl_mon3);
  mcposrpl_mon3 = rot_apply(mcrotapl_mon3, mctc1);
  mcDEBUG_COMPONENT("pl_mon3", mcposapl_mon3, mcrotapl_mon3)
  mccomp_posa[18] = mcposapl_mon3;
  mccomp_posr[18] = mcposrpl_mon3;
  mcNCounter[18]  = mcPCounter[18] = mcP2Counter[18] = 0;
  mcAbsorbProp[18]= 0;
    /* Component sample0. */
  /* Setting parameters for component sample0. */
  SIG_MESSAGE("sample0 (Init:SetPar)");
#line 166 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccsample0_radius = 0.0099;
#line 91 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccsample0_thickness = 0;
#line 91 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccsample0_zdepth = 0;
#line 91 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccsample0_Vc = 13.827;
#line 91 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccsample0_sigma_abs = 5.08;
#line 91 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccsample0_sigma_inc = 5.08;
#line 92 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccsample0_radius_i = 0;
#line 92 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccsample0_radius_o = 0;
#line 166 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccsample0_h = 0.05;
#line 92 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccsample0_focus_r = 0;
#line 92 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccsample0_pack = 1;
#line 92 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccsample0_frac = 1;
#line 92 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccsample0_f_QE = 0;
#line 92 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccsample0_gamma = 0;
#line 93 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccsample0_target_x = 0;
#line 93 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccsample0_target_y = 0;
#line 93 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccsample0_target_z = 0;
#line 166 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccsample0_focus_xw = 0.1;
#line 166 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccsample0_focus_yh = 0.01;
#line 94 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccsample0_focus_aw = 0;
#line 94 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccsample0_focus_ah = 0;
#line 94 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccsample0_xwidth = 0;
#line 94 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccsample0_yheight = 0;
#line 94 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccsample0_zthick = 0;
#line 94 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccsample0_rad_sphere = 0;
#line 94 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccsample0_sig_a = 0;
#line 94 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccsample0_sig_i = 0;
#line 94 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccsample0_V0 = 0;
#line 166 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccsample0_target_index = 4;
#line 94 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccsample0_multiples = 1;
#line 10951 "./SE_example2.c"

  SIG_MESSAGE("sample0 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 10958 "./SE_example2.c"
  rot_mul(mctr1, mcrotasample_mnt, mcrotasample0);
  rot_transpose(mcrotapl_mon3, mctr1);
  rot_mul(mcrotasample0, mctr1, mcrotrsample0);
  mctc1 = coords_set(
#line 167 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    0,
#line 167 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    0,
#line 167 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    0);
#line 10969 "./SE_example2.c"
  rot_transpose(mcrotasample_mnt, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposasample0 = coords_add(mcposasample_mnt, mctc2);
  mctc1 = coords_sub(mcposapl_mon3, mcposasample0);
  mcposrsample0 = rot_apply(mcrotasample0, mctc1);
  mcDEBUG_COMPONENT("sample0", mcposasample0, mcrotasample0)
  mccomp_posa[19] = mcposasample0;
  mccomp_posr[19] = mcposrsample0;
  mcNCounter[19]  = mcPCounter[19] = mcP2Counter[19] = 0;
  mcAbsorbProp[19]= 0;
    /* Component guide_field_2. */
  /* Setting parameters for component guide_field_2. */
  SIG_MESSAGE("guide_field_2 (Init:SetPar)");
#line 172 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccguide_field_2_xwidth = 0.1;
#line 172 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccguide_field_2_yheight = 0.1;
#line 172 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccguide_field_2_zdepth = 1;
#line 172 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccguide_field_2_Bx = 0;
#line 172 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccguide_field_2_By = 0;
#line 172 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccguide_field_2_Bz = mcipBguide;
#line 38 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  if("") strncpy(mccguide_field_2_filename, "" ? "" : "", 16384); else mccguide_field_2_filename[0]='\0';
#line 10997 "./SE_example2.c"

  SIG_MESSAGE("guide_field_2 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 11004 "./SE_example2.c"
  rot_mul(mctr1, mcrotasample_mnt, mcrotaguide_field_2);
  rot_transpose(mcrotasample0, mctr1);
  rot_mul(mcrotaguide_field_2, mctr1, mcrotrguide_field_2);
  mctc1 = coords_set(
#line 173 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    0,
#line 173 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    0,
#line 173 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    0.021);
#line 11015 "./SE_example2.c"
  rot_transpose(mcrotasample_mnt, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaguide_field_2 = coords_add(mcposasample_mnt, mctc2);
  mctc1 = coords_sub(mcposasample0, mcposaguide_field_2);
  mcposrguide_field_2 = rot_apply(mcrotaguide_field_2, mctc1);
  mcDEBUG_COMPONENT("guide_field_2", mcposaguide_field_2, mcrotaguide_field_2)
  mccomp_posa[20] = mcposaguide_field_2;
  mccomp_posr[20] = mcposrguide_field_2;
  mcNCounter[20]  = mcPCounter[20] = mcP2Counter[20] = 0;
  mcAbsorbProp[20]= 0;
    /* Component guide_field_2_end. */
  /* Setting parameters for component guide_field_2_end. */
  SIG_MESSAGE("guide_field_2_end (Init:SetPar)");
#line 46 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccguide_field_2_end_xmin = 0;
#line 46 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccguide_field_2_end_xmax = 0;
#line 46 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccguide_field_2_end_ymin = 0;
#line 46 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccguide_field_2_end_ymax = 0;
#line 46 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccguide_field_2_end_radius = 0;
#line 175 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccguide_field_2_end_xwidth = 0.1;
#line 175 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccguide_field_2_end_yheight = 0.01;
#line 11043 "./SE_example2.c"

  SIG_MESSAGE("guide_field_2_end (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 177 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    (0)*DEG2RAD,
#line 177 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    (0)*DEG2RAD,
#line 177 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    (0)*DEG2RAD);
#line 11053 "./SE_example2.c"
  rot_mul(mctr1, mcrotaSource, mcrotaguide_field_2_end);
  rot_transpose(mcrotaguide_field_2, mctr1);
  rot_mul(mcrotaguide_field_2_end, mctr1, mcrotrguide_field_2_end);
  mctc1 = coords_set(
#line 176 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    0,
#line 176 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    0,
#line 176 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    1.001);
#line 11064 "./SE_example2.c"
  rot_transpose(mcrotaguide_field_2, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaguide_field_2_end = coords_add(mcposaguide_field_2, mctc2);
  mctc1 = coords_sub(mcposaguide_field_2, mcposaguide_field_2_end);
  mcposrguide_field_2_end = rot_apply(mcrotaguide_field_2_end, mctc1);
  mcDEBUG_COMPONENT("guide_field_2_end", mcposaguide_field_2_end, mcrotaguide_field_2_end)
  mccomp_posa[21] = mcposaguide_field_2_end;
  mccomp_posr[21] = mcposrguide_field_2_end;
  mcNCounter[21]  = mcPCounter[21] = mcP2Counter[21] = 0;
  mcAbsorbProp[21]= 0;
    /* Component pl_mon40. */
  /* Setting parameters for component pl_mon40. */
  SIG_MESSAGE("pl_mon40 (Init:SetPar)");
#line 181 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  if("pol_lam40") strncpy(mccpl_mon40_filename, "pol_lam40" ? "pol_lam40" : "", 16384); else mccpl_mon40_filename[0]='\0';
#line 182 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccpl_mon40_mx = 0;
#line 182 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccpl_mon40_my = 0;
#line 182 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccpl_mon40_mz = 1;
#line 182 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccpl_mon40_Lmin = 1;
#line 183 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccpl_mon40_Lmax = 10;
#line 11090 "./SE_example2.c"

  SIG_MESSAGE("pl_mon40 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 185 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    (0)*DEG2RAD,
#line 185 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    (0)*DEG2RAD,
#line 185 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    (0)*DEG2RAD);
#line 11100 "./SE_example2.c"
  rot_mul(mctr1, mcrotaSource, mcrotapl_mon40);
  rot_transpose(mcrotaguide_field_2_end, mctr1);
  rot_mul(mcrotapl_mon40, mctr1, mcrotrpl_mon40);
  mctc1 = coords_set(
#line 184 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    0,
#line 184 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    0,
#line 184 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    1.001);
#line 11111 "./SE_example2.c"
  rot_transpose(mcrotaguide_field_2, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposapl_mon40 = coords_add(mcposaguide_field_2, mctc2);
  mctc1 = coords_sub(mcposaguide_field_2_end, mcposapl_mon40);
  mcposrpl_mon40 = rot_apply(mcrotapl_mon40, mctc1);
  mcDEBUG_COMPONENT("pl_mon40", mcposapl_mon40, mcrotapl_mon40)
  mccomp_posa[22] = mcposapl_mon40;
  mccomp_posr[22] = mcposrpl_mon40;
  mcNCounter[22]  = mcPCounter[22] = mcP2Counter[22] = 0;
  mcAbsorbProp[22]= 0;
    /* Component flipper3. */
  /* Setting parameters for component flipper3. */
  SIG_MESSAGE("flipper3 (Init:SetPar)");
#line 188 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccflipper3_xwidth = 0.1;
#line 188 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccflipper3_yheight = 0.1;
#line 188 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccflipper3_zdepth = 0.02;
#line 189 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccflipper3_Bx = 0;
#line 189 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccflipper3_By = - mcipBflip;
#line 189 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccflipper3_Bz = mcipBflip;
#line 38 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  if("") strncpy(mccflipper3_filename, "" ? "" : "", 16384); else mccflipper3_filename[0]='\0';
#line 11139 "./SE_example2.c"

  SIG_MESSAGE("flipper3 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 191 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    (0)*DEG2RAD,
#line 191 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    (0)*DEG2RAD,
#line 191 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    (0)*DEG2RAD);
#line 11149 "./SE_example2.c"
  rot_mul(mctr1, mcrotaSource, mcrotaflipper3);
  rot_transpose(mcrotapl_mon40, mctr1);
  rot_mul(mcrotaflipper3, mctr1, mcrotrflipper3);
  mctc1 = coords_set(
#line 190 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    0,
#line 190 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    0,
#line 190 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    1.001);
#line 11160 "./SE_example2.c"
  rot_transpose(mcrotaguide_field_2, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaflipper3 = coords_add(mcposaguide_field_2, mctc2);
  mctc1 = coords_sub(mcposapl_mon40, mcposaflipper3);
  mcposrflipper3 = rot_apply(mcrotaflipper3, mctc1);
  mcDEBUG_COMPONENT("flipper3", mcposaflipper3, mcrotaflipper3)
  mccomp_posa[23] = mcposaflipper3;
  mccomp_posr[23] = mcposrflipper3;
  mcNCounter[23]  = mcPCounter[23] = mcP2Counter[23] = 0;
  mcAbsorbProp[23]= 0;
    /* Component pl_mon4. */
  /* Setting parameters for component pl_mon4. */
  SIG_MESSAGE("pl_mon4 (Init:SetPar)");
#line 194 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  if("pol_lam4") strncpy(mccpl_mon4_filename, "pol_lam4" ? "pol_lam4" : "", 16384); else mccpl_mon4_filename[0]='\0';
#line 195 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccpl_mon4_mx = 0;
#line 195 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccpl_mon4_my = 0;
#line 195 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccpl_mon4_mz = 1;
#line 195 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccpl_mon4_Lmin = 1;
#line 196 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccpl_mon4_Lmax = 10;
#line 11186 "./SE_example2.c"

  SIG_MESSAGE("pl_mon4 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 198 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    (0)*DEG2RAD,
#line 198 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    (0)*DEG2RAD,
#line 198 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    (0)*DEG2RAD);
#line 11196 "./SE_example2.c"
  rot_mul(mctr1, mcrotaSource, mcrotapl_mon4);
  rot_transpose(mcrotaflipper3, mctr1);
  rot_mul(mcrotapl_mon4, mctr1, mcrotrpl_mon4);
  mctc1 = coords_set(
#line 197 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    0,
#line 197 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    0,
#line 197 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    0.021);
#line 11207 "./SE_example2.c"
  rot_transpose(mcrotaflipper3, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposapl_mon4 = coords_add(mcposaflipper3, mctc2);
  mctc1 = coords_sub(mcposaflipper3, mcposapl_mon4);
  mcposrpl_mon4 = rot_apply(mcrotapl_mon4, mctc1);
  mcDEBUG_COMPONENT("pl_mon4", mcposapl_mon4, mcrotapl_mon4)
  mccomp_posa[24] = mcposapl_mon4;
  mccomp_posr[24] = mcposrpl_mon4;
  mcNCounter[24]  = mcPCounter[24] = mcP2Counter[24] = 0;
  mcAbsorbProp[24]= 0;
    /* Component analyzer_arm. */
  /* Setting parameters for component analyzer_arm. */
  SIG_MESSAGE("analyzer_arm (Init:SetPar)");

  SIG_MESSAGE("analyzer_arm (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 204 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    (0)*DEG2RAD,
#line 204 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    (0)*DEG2RAD,
#line 204 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    (0)*DEG2RAD);
#line 11230 "./SE_example2.c"
  rot_mul(mctr1, mcrotaflipper3, mcrotaanalyzer_arm);
  rot_transpose(mcrotapl_mon4, mctr1);
  rot_mul(mcrotaanalyzer_arm, mctr1, mcrotranalyzer_arm);
  mctc1 = coords_set(
#line 203 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    0,
#line 203 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    0,
#line 203 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    0.32 + 0.021);
#line 11241 "./SE_example2.c"
  rot_transpose(mcrotaflipper3, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaanalyzer_arm = coords_add(mcposaflipper3, mctc2);
  mctc1 = coords_sub(mcposapl_mon4, mcposaanalyzer_arm);
  mcposranalyzer_arm = rot_apply(mcrotaanalyzer_arm, mctc1);
  mcDEBUG_COMPONENT("analyzer_arm", mcposaanalyzer_arm, mcrotaanalyzer_arm)
  mccomp_posa[25] = mcposaanalyzer_arm;
  mccomp_posr[25] = mcposranalyzer_arm;
  mcNCounter[25]  = mcPCounter[25] = mcP2Counter[25] = 0;
  mcAbsorbProp[25]= 0;
    /* Component psd_ba. */
  /* Setting parameters for component psd_ba. */
  SIG_MESSAGE("psd_ba (Init:SetPar)");
#line 49 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccpsd_ba_nx = 90;
#line 49 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccpsd_ba_ny = 90;
#line 207 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  if("psd_ap") strncpy(mccpsd_ba_filename, "psd_ap" ? "psd_ap" : "", 16384); else mccpsd_ba_filename[0]='\0';
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccpsd_ba_xmin = -0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccpsd_ba_xmax = 0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccpsd_ba_ymin = -0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccpsd_ba_ymax = 0.05;
#line 207 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccpsd_ba_xwidth = 0.1;
#line 207 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccpsd_ba_yheight = 0.01;
#line 207 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccpsd_ba_restore_neutron = 1;
#line 11275 "./SE_example2.c"

  SIG_MESSAGE("psd_ba (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 209 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    (0)*DEG2RAD,
#line 209 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    (0)*DEG2RAD,
#line 209 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    (0)*DEG2RAD);
#line 11285 "./SE_example2.c"
  rot_mul(mctr1, mcrotaSource, mcrotapsd_ba);
  rot_transpose(mcrotaanalyzer_arm, mctr1);
  rot_mul(mcrotapsd_ba, mctr1, mcrotrpsd_ba);
  mctc1 = coords_set(
#line 208 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    0,
#line 208 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    0,
#line 208 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    0);
#line 11296 "./SE_example2.c"
  rot_transpose(mcrotaanalyzer_arm, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposapsd_ba = coords_add(mcposaanalyzer_arm, mctc2);
  mctc1 = coords_sub(mcposaanalyzer_arm, mcposapsd_ba);
  mcposrpsd_ba = rot_apply(mcrotapsd_ba, mctc1);
  mcDEBUG_COMPONENT("psd_ba", mcposapsd_ba, mcrotapsd_ba)
  mccomp_posa[26] = mcposapsd_ba;
  mccomp_posr[26] = mcposrpsd_ba;
  mcNCounter[26]  = mcPCounter[26] = mcP2Counter[26] = 0;
  mcAbsorbProp[26]= 0;
    /* Component pl_mon_ba. */
  /* Setting parameters for component pl_mon_ba. */
  SIG_MESSAGE("pl_mon_ba (Init:SetPar)");
#line 212 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  if("lam_pol_ba") strncpy(mccpl_mon_ba_filename, "lam_pol_ba" ? "lam_pol_ba" : "", 16384); else mccpl_mon_ba_filename[0]='\0';
#line 213 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccpl_mon_ba_mx = 1;
#line 213 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccpl_mon_ba_my = 0;
#line 213 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccpl_mon_ba_mz = 0;
#line 213 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccpl_mon_ba_Lmin = 1;
#line 214 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccpl_mon_ba_Lmax = 10;
#line 11322 "./SE_example2.c"

  SIG_MESSAGE("pl_mon_ba (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 216 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    (0)*DEG2RAD,
#line 216 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    (0)*DEG2RAD,
#line 216 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    (0)*DEG2RAD);
#line 11332 "./SE_example2.c"
  rot_mul(mctr1, mcrotaSource, mcrotapl_mon_ba);
  rot_transpose(mcrotapsd_ba, mctr1);
  rot_mul(mcrotapl_mon_ba, mctr1, mcrotrpl_mon_ba);
  mctc1 = coords_set(
#line 215 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    0,
#line 215 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    0,
#line 215 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    0);
#line 11343 "./SE_example2.c"
  rot_transpose(mcrotaanalyzer_arm, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposapl_mon_ba = coords_add(mcposaanalyzer_arm, mctc2);
  mctc1 = coords_sub(mcposapsd_ba, mcposapl_mon_ba);
  mcposrpl_mon_ba = rot_apply(mcrotapl_mon_ba, mctc1);
  mcDEBUG_COMPONENT("pl_mon_ba", mcposapl_mon_ba, mcrotapl_mon_ba)
  mccomp_posa[27] = mcposapl_mon_ba;
  mccomp_posr[27] = mcposrpl_mon_ba;
  mcNCounter[27]  = mcPCounter[27] = mcP2Counter[27] = 0;
  mcAbsorbProp[27]= 0;
    /* Component analyzer. */
  /* Setting parameters for component analyzer. */
  SIG_MESSAGE("analyzer (Init:SetPar)");
#line 221 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccanalyzer_p_reflect = 0;
#line 11359 "./SE_example2.c"

  SIG_MESSAGE("analyzer (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 223 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    (90 - mcipPOL_ANGLE)*DEG2RAD,
#line 223 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    (90)*DEG2RAD,
#line 223 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    (0)*DEG2RAD);
#line 11369 "./SE_example2.c"
  rot_mul(mctr1, mcrotaanalyzer_arm, mcrotaanalyzer);
  rot_transpose(mcrotapl_mon_ba, mctr1);
  rot_mul(mcrotaanalyzer, mctr1, mcrotranalyzer);
  mctc1 = coords_set(
#line 222 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    0,
#line 222 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    0,
#line 222 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    0);
#line 11380 "./SE_example2.c"
  rot_transpose(mcrotaanalyzer_arm, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaanalyzer = coords_add(mcposaanalyzer_arm, mctc2);
  mctc1 = coords_sub(mcposapl_mon_ba, mcposaanalyzer);
  mcposranalyzer = rot_apply(mcrotaanalyzer, mctc1);
  mcDEBUG_COMPONENT("analyzer", mcposaanalyzer, mcrotaanalyzer)
  mccomp_posa[28] = mcposaanalyzer;
  mccomp_posr[28] = mcposranalyzer;
  mcNCounter[28]  = mcPCounter[28] = mcP2Counter[28] = 0;
  mcAbsorbProp[28]= 0;
    /* Component detector_slit. */
  /* Setting parameters for component detector_slit. */
  SIG_MESSAGE("detector_slit (Init:SetPar)");
#line 46 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccdetector_slit_xmin = 0;
#line 46 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccdetector_slit_xmax = 0;
#line 46 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccdetector_slit_ymin = 0;
#line 46 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccdetector_slit_ymax = 0;
#line 46 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccdetector_slit_radius = 0;
#line 228 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccdetector_slit_xwidth = 0.1;
#line 228 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccdetector_slit_yheight = 0.01;
#line 11408 "./SE_example2.c"

  SIG_MESSAGE("detector_slit (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 230 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    (0)*DEG2RAD,
#line 230 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    (0)*DEG2RAD,
#line 230 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    (0)*DEG2RAD);
#line 11418 "./SE_example2.c"
  rot_mul(mctr1, mcrotaSource, mcrotadetector_slit);
  rot_transpose(mcrotaanalyzer, mctr1);
  rot_mul(mcrotadetector_slit, mctr1, mcrotrdetector_slit);
  mctc1 = coords_set(
#line 229 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    0,
#line 229 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    0,
#line 229 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    0.64);
#line 11429 "./SE_example2.c"
  rot_transpose(mcrotaanalyzer_arm, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposadetector_slit = coords_add(mcposaanalyzer_arm, mctc2);
  mctc1 = coords_sub(mcposaanalyzer, mcposadetector_slit);
  mcposrdetector_slit = rot_apply(mcrotadetector_slit, mctc1);
  mcDEBUG_COMPONENT("detector_slit", mcposadetector_slit, mcrotadetector_slit)
  mccomp_posa[29] = mcposadetector_slit;
  mccomp_posr[29] = mcposrdetector_slit;
  mcNCounter[29]  = mcPCounter[29] = mcP2Counter[29] = 0;
  mcAbsorbProp[29]= 0;
    /* Component pl_mon5. */
  /* Setting parameters for component pl_mon5. */
  SIG_MESSAGE("pl_mon5 (Init:SetPar)");
#line 233 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  if("pol_lam5") strncpy(mccpl_mon5_filename, "pol_lam5" ? "pol_lam5" : "", 16384); else mccpl_mon5_filename[0]='\0';
#line 234 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccpl_mon5_mx = 0;
#line 234 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccpl_mon5_my = 0;
#line 234 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccpl_mon5_mz = 1;
#line 234 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccpl_mon5_Lmin = 1;
#line 235 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccpl_mon5_Lmax = 10;
#line 11455 "./SE_example2.c"

  SIG_MESSAGE("pl_mon5 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 237 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    (0)*DEG2RAD,
#line 237 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    (0)*DEG2RAD,
#line 237 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    (0)*DEG2RAD);
#line 11465 "./SE_example2.c"
  rot_mul(mctr1, mcrotaSource, mcrotapl_mon5);
  rot_transpose(mcrotadetector_slit, mctr1);
  rot_mul(mcrotapl_mon5, mctr1, mcrotrpl_mon5);
  mctc1 = coords_set(
#line 236 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    0,
#line 236 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    0,
#line 236 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    0.641);
#line 11476 "./SE_example2.c"
  rot_transpose(mcrotaanalyzer_arm, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposapl_mon5 = coords_add(mcposaanalyzer_arm, mctc2);
  mctc1 = coords_sub(mcposadetector_slit, mcposapl_mon5);
  mcposrpl_mon5 = rot_apply(mcrotapl_mon5, mctc1);
  mcDEBUG_COMPONENT("pl_mon5", mcposapl_mon5, mcrotapl_mon5)
  mccomp_posa[30] = mcposapl_mon5;
  mccomp_posr[30] = mcposrpl_mon5;
  mcNCounter[30]  = mcPCounter[30] = mcP2Counter[30] = 0;
  mcAbsorbProp[30]= 0;
    /* Component detector. */
  /* Setting parameters for component detector. */
  SIG_MESSAGE("detector (Init:SetPar)");
#line 49 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccdetector_nx = 90;
#line 49 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccdetector_ny = 90;
#line 240 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  if("detector") strncpy(mccdetector_filename, "detector" ? "detector" : "", 16384); else mccdetector_filename[0]='\0';
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccdetector_xmin = -0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccdetector_xmax = 0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccdetector_ymin = -0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccdetector_ymax = 0.05;
#line 240 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccdetector_xwidth = 0.1;
#line 240 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccdetector_yheight = 0.1;
#line 51 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  mccdetector_restore_neutron = 0;
#line 11510 "./SE_example2.c"

  SIG_MESSAGE("detector (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 242 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    (0)*DEG2RAD,
#line 242 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    (0)*DEG2RAD,
#line 242 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    (0)*DEG2RAD);
#line 11520 "./SE_example2.c"
  rot_mul(mctr1, mcrotaSource, mcrotadetector);
  rot_transpose(mcrotapl_mon5, mctr1);
  rot_mul(mcrotadetector, mctr1, mcrotrdetector);
  mctc1 = coords_set(
#line 241 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    0,
#line 241 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    0,
#line 241 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
    0.642);
#line 11531 "./SE_example2.c"
  rot_transpose(mcrotaanalyzer_arm, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposadetector = coords_add(mcposaanalyzer_arm, mctc2);
  mctc1 = coords_sub(mcposapl_mon5, mcposadetector);
  mcposrdetector = rot_apply(mcrotadetector, mctc1);
  mcDEBUG_COMPONENT("detector", mcposadetector, mcrotadetector)
  mccomp_posa[31] = mcposadetector;
  mccomp_posr[31] = mcposrdetector;
  mcNCounter[31]  = mcPCounter[31] = mcP2Counter[31] = 0;
  mcAbsorbProp[31]= 0;
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
#line 11568 "./SE_example2.c"
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
#line 11662 "./SE_example2.c"
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

  /* Initializations for component lam_mon. */
  SIG_MESSAGE("lam_mon (Init)");
#define mccompcurname  lam_mon
#define mccompcurtype  L_monitor
#define mccompcurindex 3
#define nL mcclam_mon_nL
#define L_N mcclam_mon_L_N
#define L_p mcclam_mon_L_p
#define L_p2 mcclam_mon_L_p2
#define filename mcclam_mon_filename
#define xmin mcclam_mon_xmin
#define xmax mcclam_mon_xmax
#define ymin mcclam_mon_ymin
#define ymax mcclam_mon_ymax
#define xwidth mcclam_mon_xwidth
#define yheight mcclam_mon_yheight
#define Lmin mcclam_mon_Lmin
#define Lmax mcclam_mon_Lmax
#define restore_neutron mcclam_mon_restore_neutron
#define nowritefile mcclam_mon_nowritefile
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
#line 11724 "./SE_example2.c"
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

  /* Initializations for component polariser_pt. */
  SIG_MESSAGE("polariser_pt (Init)");

  /* Initializations for component hdiv_mon. */
  SIG_MESSAGE("hdiv_mon (Init)");
#define mccompcurname  hdiv_mon
#define mccompcurtype  Hdiv_monitor
#define mccompcurindex 5
#define nh mcchdiv_mon_nh
#define Div_N mcchdiv_mon_Div_N
#define Div_p mcchdiv_mon_Div_p
#define Div_p2 mcchdiv_mon_Div_p2
#define filename mcchdiv_mon_filename
#define xmin mcchdiv_mon_xmin
#define xmax mcchdiv_mon_xmax
#define ymin mcchdiv_mon_ymin
#define ymax mcchdiv_mon_ymax
#define xwidth mcchdiv_mon_xwidth
#define yheight mcchdiv_mon_yheight
#define h_maxdiv mcchdiv_mon_h_maxdiv
#define restore_neutron mcchdiv_mon_restore_neutron
#define nowritefile mcchdiv_mon_nowritefile
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
#line 11789 "./SE_example2.c"
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

  /* Initializations for component polariser2. */
  SIG_MESSAGE("polariser2 (Init)");
#define mccompcurname  polariser2
#define mccompcurtype  Pol_mirror
#define mccompcurindex 6
#define rUpFunc mccpolariser2_rUpFunc
#define rDownFunc mccpolariser2_rDownFunc
#define rUpPar mccpolariser2_rUpPar
#define rDownPar mccpolariser2_rDownPar
#define useTables mccpolariser2_useTables
#define zwidth mccpolariser2_zwidth
#define yheight mccpolariser2_yheight
#define rUpParPtr mccpolariser2_rUpParPtr
#define rDownParPtr mccpolariser2_rDownParPtr
#define p_reflect mccpolariser2_p_reflect
#line 105 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Pol_mirror.comp"
{
#if (useTables)
  rUpParPtr   = (t_Table*) malloc(sizeof(t_Table));
  rDownParPtr = (t_Table*) malloc(sizeof(t_Table));
  if (Table_Read(rUpParPtr, rUpPar, 1) <= 0) {
    fprintf(stderr,"Pol_mirror: %s: can not read file %s\n",
            NAME_CURRENT_COMP, rUpPar);
    exit(1);
  }
  if (Table_Read(rDownParPtr, rDownPar, 1) <= 0) {
    fprintf(stderr,"Pol_mirror: %s: can not read file %s\n",
            NAME_CURRENT_COMP, rDownPar);
    exit(1);
  }
#endif

  if ((zwidth<=0) || (yheight <= 0)) {
    fprintf(stderr, "Pol_mirror: %s: NULL or negative length scale!\n"
            "ERROR      (zwidth=%g,yheight=%g). Exiting\n",
            NAME_CURRENT_COMP, zwidth, yheight);
    exit(1);
  }

}
#line 11848 "./SE_example2.c"
#undef p_reflect
#undef rDownParPtr
#undef rUpParPtr
#undef yheight
#undef zwidth
#undef useTables
#undef rDownPar
#undef rUpPar
#undef rDownFunc
#undef rUpFunc
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component psd_ap. */
  SIG_MESSAGE("psd_ap (Init)");
#define mccompcurname  psd_ap
#define mccompcurtype  PSD_monitor
#define mccompcurindex 7
#define PSD_N mccpsd_ap_PSD_N
#define PSD_p mccpsd_ap_PSD_p
#define PSD_p2 mccpsd_ap_PSD_p2
#define nx mccpsd_ap_nx
#define ny mccpsd_ap_ny
#define filename mccpsd_ap_filename
#define xmin mccpsd_ap_xmin
#define xmax mccpsd_ap_xmax
#define ymin mccpsd_ap_ymin
#define ymax mccpsd_ap_ymax
#define xwidth mccpsd_ap_xwidth
#define yheight mccpsd_ap_yheight
#define restore_neutron mccpsd_ap_restore_neutron
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
#line 11906 "./SE_example2.c"
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

  /* Initializations for component pl_mon0. */
  SIG_MESSAGE("pl_mon0 (Init)");
#define mccompcurname  pl_mon0
#define mccompcurtype  MeanPolLambda_monitor
#define mccompcurindex 8
#define xwidth mccpl_mon0_xwidth
#define yheight mccpl_mon0_yheight
#define nL mccpl_mon0_nL
#define restore_neutron mccpl_mon0_restore_neutron
#define nowritefile mccpl_mon0_nowritefile
#define PolL_N mccpl_mon0_PolL_N
#define PolL_p mccpl_mon0_PolL_p
#define PolL_p2 mccpl_mon0_PolL_p2
#define HelpArray mccpl_mon0_HelpArray
#define filename mccpl_mon0_filename
#define mx mccpl_mon0_mx
#define my mccpl_mon0_my
#define mz mccpl_mon0_mz
#define Lmin mccpl_mon0_Lmin
#define Lmax mccpl_mon0_Lmax
#line 66 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/MeanPolLambda_monitor.comp"
{
int i;

// Check that input parameteters makes sense

if (Lmax<=Lmin) {
    fprintf(stderr, "Pol_monitor: %s: l1 <= l0!\n"
	   "ERROR. Exiting",
           NAME_CURRENT_COMP);
    exit(1);
  }

  if (mx==0 && my==0 && mz==0) {
    fprintf(stderr, "Pol_monitor: %s: NULL vector defined!\n"
	   "ERROR      (mx, my, mz). Exiting",
           NAME_CURRENT_COMP);
    exit(1);
  }

  if ((xwidth<=0) || (yheight <= 0)) {
    fprintf(stderr, "Pol_monitor: %s: Null detection area !\n"
	   "ERROR      (xwidth,yheight). Exiting",
           NAME_CURRENT_COMP);
    exit(1);
  }

  // Initialize variables

  NORM(mx, my, mz);

  for (i=0; i<nL; i++) {

      PolL_N[i] = 0;
      PolL_p[i] = 0;
      PolL_p2[i] = 0;
      HelpArray[i] = 0;
  }
}
#line 11983 "./SE_example2.c"
#undef Lmax
#undef Lmin
#undef mz
#undef my
#undef mx
#undef filename
#undef HelpArray
#undef PolL_p2
#undef PolL_p
#undef PolL_N
#undef nowritefile
#undef restore_neutron
#undef nL
#undef yheight
#undef xwidth
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component flipper1. */
  SIG_MESSAGE("flipper1 (Init)");
#define mccompcurname  flipper1
#define mccompcurtype  Pol_FieldBox
#define mccompcurindex 9
#define fieldFunction mccflipper1_fieldFunction
#define field_parameters mccflipper1_field_parameters
#define const_field_parameters mccflipper1_const_field_parameters
#define xwidth mccflipper1_xwidth
#define yheight mccflipper1_yheight
#define zdepth mccflipper1_zdepth
#define Bx mccflipper1_Bx
#define By mccflipper1_By
#define Bz mccflipper1_Bz
#define filename mccflipper1_filename
#line 55 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Pol_FieldBox.comp"
{
  if (fieldFunction==table_magnetic_field){
      /*initialize the interpolation vertex structure*/
      struct interpolator_struct *interpolator = interpolator_load(filename, 3, 3, NULL);
      interpolator_info(interpolator);
      field_parameters=(void *)interpolator;
  }  else {
      field_parameters = (void *) const_field_parameters;
  }

}
#line 12030 "./SE_example2.c"
#undef filename
#undef Bz
#undef By
#undef Bx
#undef zdepth
#undef yheight
#undef xwidth
#undef const_field_parameters
#undef field_parameters
#undef fieldFunction
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component pl_mon1. */
  SIG_MESSAGE("pl_mon1 (Init)");
#define mccompcurname  pl_mon1
#define mccompcurtype  MeanPolLambda_monitor
#define mccompcurindex 10
#define xwidth mccpl_mon1_xwidth
#define yheight mccpl_mon1_yheight
#define nL mccpl_mon1_nL
#define restore_neutron mccpl_mon1_restore_neutron
#define nowritefile mccpl_mon1_nowritefile
#define PolL_N mccpl_mon1_PolL_N
#define PolL_p mccpl_mon1_PolL_p
#define PolL_p2 mccpl_mon1_PolL_p2
#define HelpArray mccpl_mon1_HelpArray
#define filename mccpl_mon1_filename
#define mx mccpl_mon1_mx
#define my mccpl_mon1_my
#define mz mccpl_mon1_mz
#define Lmin mccpl_mon1_Lmin
#define Lmax mccpl_mon1_Lmax
#line 66 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/MeanPolLambda_monitor.comp"
{
int i;

// Check that input parameteters makes sense

if (Lmax<=Lmin) {
    fprintf(stderr, "Pol_monitor: %s: l1 <= l0!\n"
	   "ERROR. Exiting",
           NAME_CURRENT_COMP);
    exit(1);
  }

  if (mx==0 && my==0 && mz==0) {
    fprintf(stderr, "Pol_monitor: %s: NULL vector defined!\n"
	   "ERROR      (mx, my, mz). Exiting",
           NAME_CURRENT_COMP);
    exit(1);
  }

  if ((xwidth<=0) || (yheight <= 0)) {
    fprintf(stderr, "Pol_monitor: %s: Null detection area !\n"
	   "ERROR      (xwidth,yheight). Exiting",
           NAME_CURRENT_COMP);
    exit(1);
  }

  // Initialize variables

  NORM(mx, my, mz);

  for (i=0; i<nL; i++) {

      PolL_N[i] = 0;
      PolL_p[i] = 0;
      PolL_p2[i] = 0;
      HelpArray[i] = 0;
  }
}
#line 12104 "./SE_example2.c"
#undef Lmax
#undef Lmin
#undef mz
#undef my
#undef mx
#undef filename
#undef HelpArray
#undef PolL_p2
#undef PolL_p
#undef PolL_N
#undef nowritefile
#undef restore_neutron
#undef nL
#undef yheight
#undef xwidth
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component guide_field_1. */
  SIG_MESSAGE("guide_field_1 (Init)");
#define mccompcurname  guide_field_1
#define mccompcurtype  Pol_FieldBox
#define mccompcurindex 11
#define fieldFunction mccguide_field_1_fieldFunction
#define field_parameters mccguide_field_1_field_parameters
#define const_field_parameters mccguide_field_1_const_field_parameters
#define xwidth mccguide_field_1_xwidth
#define yheight mccguide_field_1_yheight
#define zdepth mccguide_field_1_zdepth
#define Bx mccguide_field_1_Bx
#define By mccguide_field_1_By
#define Bz mccguide_field_1_Bz
#define filename mccguide_field_1_filename
#line 55 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Pol_FieldBox.comp"
{
  if (fieldFunction==table_magnetic_field){
      /*initialize the interpolation vertex structure*/
      struct interpolator_struct *interpolator = interpolator_load(filename, 3, 3, NULL);
      interpolator_info(interpolator);
      field_parameters=(void *)interpolator;
  }  else {
      field_parameters = (void *) const_field_parameters;
  }

}
#line 12151 "./SE_example2.c"
#undef filename
#undef Bz
#undef By
#undef Bx
#undef zdepth
#undef yheight
#undef xwidth
#undef const_field_parameters
#undef field_parameters
#undef fieldFunction
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component guide_field_1_end. */
  SIG_MESSAGE("guide_field_1_end (Init)");
#define mccompcurname  guide_field_1_end
#define mccompcurtype  Slit
#define mccompcurindex 12
#define xmin mccguide_field_1_end_xmin
#define xmax mccguide_field_1_end_xmax
#define ymin mccguide_field_1_end_ymin
#define ymax mccguide_field_1_end_ymax
#define radius mccguide_field_1_end_radius
#define xwidth mccguide_field_1_end_xwidth
#define yheight mccguide_field_1_end_yheight
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
#line 12198 "./SE_example2.c"
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

  /* Initializations for component pl_mon2. */
  SIG_MESSAGE("pl_mon2 (Init)");
#define mccompcurname  pl_mon2
#define mccompcurtype  MeanPolLambda_monitor
#define mccompcurindex 13
#define xwidth mccpl_mon2_xwidth
#define yheight mccpl_mon2_yheight
#define nL mccpl_mon2_nL
#define restore_neutron mccpl_mon2_restore_neutron
#define nowritefile mccpl_mon2_nowritefile
#define PolL_N mccpl_mon2_PolL_N
#define PolL_p mccpl_mon2_PolL_p
#define PolL_p2 mccpl_mon2_PolL_p2
#define HelpArray mccpl_mon2_HelpArray
#define filename mccpl_mon2_filename
#define mx mccpl_mon2_mx
#define my mccpl_mon2_my
#define mz mccpl_mon2_mz
#define Lmin mccpl_mon2_Lmin
#define Lmax mccpl_mon2_Lmax
#line 66 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/MeanPolLambda_monitor.comp"
{
int i;

// Check that input parameteters makes sense

if (Lmax<=Lmin) {
    fprintf(stderr, "Pol_monitor: %s: l1 <= l0!\n"
	   "ERROR. Exiting",
           NAME_CURRENT_COMP);
    exit(1);
  }

  if (mx==0 && my==0 && mz==0) {
    fprintf(stderr, "Pol_monitor: %s: NULL vector defined!\n"
	   "ERROR      (mx, my, mz). Exiting",
           NAME_CURRENT_COMP);
    exit(1);
  }

  if ((xwidth<=0) || (yheight <= 0)) {
    fprintf(stderr, "Pol_monitor: %s: Null detection area !\n"
	   "ERROR      (xwidth,yheight). Exiting",
           NAME_CURRENT_COMP);
    exit(1);
  }

  // Initialize variables

  NORM(mx, my, mz);

  for (i=0; i<nL; i++) {

      PolL_N[i] = 0;
      PolL_p[i] = 0;
      PolL_p2[i] = 0;
      HelpArray[i] = 0;
  }
}
#line 12269 "./SE_example2.c"
#undef Lmax
#undef Lmin
#undef mz
#undef my
#undef mx
#undef filename
#undef HelpArray
#undef PolL_p2
#undef PolL_p
#undef PolL_N
#undef nowritefile
#undef restore_neutron
#undef nL
#undef yheight
#undef xwidth
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component flipper2_1. */
  SIG_MESSAGE("flipper2_1 (Init)");
#define mccompcurname  flipper2_1
#define mccompcurtype  Pol_FieldBox
#define mccompcurindex 14
#define fieldFunction mccflipper2_1_fieldFunction
#define field_parameters mccflipper2_1_field_parameters
#define const_field_parameters mccflipper2_1_const_field_parameters
#define xwidth mccflipper2_1_xwidth
#define yheight mccflipper2_1_yheight
#define zdepth mccflipper2_1_zdepth
#define Bx mccflipper2_1_Bx
#define By mccflipper2_1_By
#define Bz mccflipper2_1_Bz
#define filename mccflipper2_1_filename
#line 55 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Pol_FieldBox.comp"
{
  if (fieldFunction==table_magnetic_field){
      /*initialize the interpolation vertex structure*/
      struct interpolator_struct *interpolator = interpolator_load(filename, 3, 3, NULL);
      interpolator_info(interpolator);
      field_parameters=(void *)interpolator;
  }  else {
      field_parameters = (void *) const_field_parameters;
  }

}
#line 12316 "./SE_example2.c"
#undef filename
#undef Bz
#undef By
#undef Bx
#undef zdepth
#undef yheight
#undef xwidth
#undef const_field_parameters
#undef field_parameters
#undef fieldFunction
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component flipper2_2. */
  SIG_MESSAGE("flipper2_2 (Init)");
#define mccompcurname  flipper2_2
#define mccompcurtype  Pol_FieldBox
#define mccompcurindex 15
#define fieldFunction mccflipper2_2_fieldFunction
#define field_parameters mccflipper2_2_field_parameters
#define const_field_parameters mccflipper2_2_const_field_parameters
#define xwidth mccflipper2_2_xwidth
#define yheight mccflipper2_2_yheight
#define zdepth mccflipper2_2_zdepth
#define Bx mccflipper2_2_Bx
#define By mccflipper2_2_By
#define Bz mccflipper2_2_Bz
#define filename mccflipper2_2_filename
#line 55 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Pol_FieldBox.comp"
{
  if (fieldFunction==table_magnetic_field){
      /*initialize the interpolation vertex structure*/
      struct interpolator_struct *interpolator = interpolator_load(filename, 3, 3, NULL);
      interpolator_info(interpolator);
      field_parameters=(void *)interpolator;
  }  else {
      field_parameters = (void *) const_field_parameters;
  }

}
#line 12358 "./SE_example2.c"
#undef filename
#undef Bz
#undef By
#undef Bx
#undef zdepth
#undef yheight
#undef xwidth
#undef const_field_parameters
#undef field_parameters
#undef fieldFunction
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component flipper_2_2_end. */
  SIG_MESSAGE("flipper_2_2_end (Init)");
#define mccompcurname  flipper_2_2_end
#define mccompcurtype  Slit
#define mccompcurindex 16
#define xmin mccflipper_2_2_end_xmin
#define xmax mccflipper_2_2_end_xmax
#define ymin mccflipper_2_2_end_ymin
#define ymax mccflipper_2_2_end_ymax
#define radius mccflipper_2_2_end_radius
#define xwidth mccflipper_2_2_end_xwidth
#define yheight mccflipper_2_2_end_yheight
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
#line 12405 "./SE_example2.c"
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

  /* Initializations for component sample_mnt. */
  SIG_MESSAGE("sample_mnt (Init)");

  /* Initializations for component pl_mon3. */
  SIG_MESSAGE("pl_mon3 (Init)");
#define mccompcurname  pl_mon3
#define mccompcurtype  MeanPolLambda_monitor
#define mccompcurindex 18
#define xwidth mccpl_mon3_xwidth
#define yheight mccpl_mon3_yheight
#define nL mccpl_mon3_nL
#define restore_neutron mccpl_mon3_restore_neutron
#define nowritefile mccpl_mon3_nowritefile
#define PolL_N mccpl_mon3_PolL_N
#define PolL_p mccpl_mon3_PolL_p
#define PolL_p2 mccpl_mon3_PolL_p2
#define HelpArray mccpl_mon3_HelpArray
#define filename mccpl_mon3_filename
#define mx mccpl_mon3_mx
#define my mccpl_mon3_my
#define mz mccpl_mon3_mz
#define Lmin mccpl_mon3_Lmin
#define Lmax mccpl_mon3_Lmax
#line 66 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/MeanPolLambda_monitor.comp"
{
int i;

// Check that input parameteters makes sense

if (Lmax<=Lmin) {
    fprintf(stderr, "Pol_monitor: %s: l1 <= l0!\n"
	   "ERROR. Exiting",
           NAME_CURRENT_COMP);
    exit(1);
  }

  if (mx==0 && my==0 && mz==0) {
    fprintf(stderr, "Pol_monitor: %s: NULL vector defined!\n"
	   "ERROR      (mx, my, mz). Exiting",
           NAME_CURRENT_COMP);
    exit(1);
  }

  if ((xwidth<=0) || (yheight <= 0)) {
    fprintf(stderr, "Pol_monitor: %s: Null detection area !\n"
	   "ERROR      (xwidth,yheight). Exiting",
           NAME_CURRENT_COMP);
    exit(1);
  }

  // Initialize variables

  NORM(mx, my, mz);

  for (i=0; i<nL; i++) {

      PolL_N[i] = 0;
      PolL_p[i] = 0;
      PolL_p2[i] = 0;
      HelpArray[i] = 0;
  }
}
#line 12479 "./SE_example2.c"
#undef Lmax
#undef Lmin
#undef mz
#undef my
#undef mx
#undef filename
#undef HelpArray
#undef PolL_p2
#undef PolL_p
#undef PolL_N
#undef nowritefile
#undef restore_neutron
#undef nL
#undef yheight
#undef xwidth
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component sample0. */
  SIG_MESSAGE("sample0 (Init)");
#define mccompcurname  sample0
#define mccompcurtype  V_sample
#define mccompcurindex 19
#define VarsV mccsample0_VarsV
#define radius mccsample0_radius
#define thickness mccsample0_thickness
#define zdepth mccsample0_zdepth
#define Vc mccsample0_Vc
#define sigma_abs mccsample0_sigma_abs
#define sigma_inc mccsample0_sigma_inc
#define radius_i mccsample0_radius_i
#define radius_o mccsample0_radius_o
#define h mccsample0_h
#define focus_r mccsample0_focus_r
#define pack mccsample0_pack
#define frac mccsample0_frac
#define f_QE mccsample0_f_QE
#define gamma mccsample0_gamma
#define target_x mccsample0_target_x
#define target_y mccsample0_target_y
#define target_z mccsample0_target_z
#define focus_xw mccsample0_focus_xw
#define focus_yh mccsample0_focus_yh
#define focus_aw mccsample0_focus_aw
#define focus_ah mccsample0_focus_ah
#define xwidth mccsample0_xwidth
#define yheight mccsample0_yheight
#define zthick mccsample0_zthick
#define rad_sphere mccsample0_rad_sphere
#define sig_a mccsample0_sig_a
#define sig_i mccsample0_sig_i
#define V0 mccsample0_V0
#define target_index mccsample0_target_index
#define multiples mccsample0_multiples
#line 121 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../obsolete/V_sample.comp"
{
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
}
#line 12593 "./SE_example2.c"
#undef multiples
#undef target_index
#undef V0
#undef sig_i
#undef sig_a
#undef rad_sphere
#undef zthick
#undef yheight
#undef xwidth
#undef focus_ah
#undef focus_aw
#undef focus_yh
#undef focus_xw
#undef target_z
#undef target_y
#undef target_x
#undef gamma
#undef f_QE
#undef frac
#undef pack
#undef focus_r
#undef h
#undef radius_o
#undef radius_i
#undef sigma_inc
#undef sigma_abs
#undef Vc
#undef zdepth
#undef thickness
#undef radius
#undef VarsV
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component guide_field_2. */
  SIG_MESSAGE("guide_field_2 (Init)");
#define mccompcurname  guide_field_2
#define mccompcurtype  Pol_FieldBox
#define mccompcurindex 20
#define fieldFunction mccguide_field_2_fieldFunction
#define field_parameters mccguide_field_2_field_parameters
#define const_field_parameters mccguide_field_2_const_field_parameters
#define xwidth mccguide_field_2_xwidth
#define yheight mccguide_field_2_yheight
#define zdepth mccguide_field_2_zdepth
#define Bx mccguide_field_2_Bx
#define By mccguide_field_2_By
#define Bz mccguide_field_2_Bz
#define filename mccguide_field_2_filename
#line 55 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Pol_FieldBox.comp"
{
  if (fieldFunction==table_magnetic_field){
      /*initialize the interpolation vertex structure*/
      struct interpolator_struct *interpolator = interpolator_load(filename, 3, 3, NULL);
      interpolator_info(interpolator);
      field_parameters=(void *)interpolator;
  }  else {
      field_parameters = (void *) const_field_parameters;
  }

}
#line 12656 "./SE_example2.c"
#undef filename
#undef Bz
#undef By
#undef Bx
#undef zdepth
#undef yheight
#undef xwidth
#undef const_field_parameters
#undef field_parameters
#undef fieldFunction
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component guide_field_2_end. */
  SIG_MESSAGE("guide_field_2_end (Init)");
#define mccompcurname  guide_field_2_end
#define mccompcurtype  Slit
#define mccompcurindex 21
#define xmin mccguide_field_2_end_xmin
#define xmax mccguide_field_2_end_xmax
#define ymin mccguide_field_2_end_ymin
#define ymax mccguide_field_2_end_ymax
#define radius mccguide_field_2_end_radius
#define xwidth mccguide_field_2_end_xwidth
#define yheight mccguide_field_2_end_yheight
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
#line 12703 "./SE_example2.c"
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

  /* Initializations for component pl_mon40. */
  SIG_MESSAGE("pl_mon40 (Init)");
#define mccompcurname  pl_mon40
#define mccompcurtype  MeanPolLambda_monitor
#define mccompcurindex 22
#define xwidth mccpl_mon40_xwidth
#define yheight mccpl_mon40_yheight
#define nL mccpl_mon40_nL
#define restore_neutron mccpl_mon40_restore_neutron
#define nowritefile mccpl_mon40_nowritefile
#define PolL_N mccpl_mon40_PolL_N
#define PolL_p mccpl_mon40_PolL_p
#define PolL_p2 mccpl_mon40_PolL_p2
#define HelpArray mccpl_mon40_HelpArray
#define filename mccpl_mon40_filename
#define mx mccpl_mon40_mx
#define my mccpl_mon40_my
#define mz mccpl_mon40_mz
#define Lmin mccpl_mon40_Lmin
#define Lmax mccpl_mon40_Lmax
#line 66 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/MeanPolLambda_monitor.comp"
{
int i;

// Check that input parameteters makes sense

if (Lmax<=Lmin) {
    fprintf(stderr, "Pol_monitor: %s: l1 <= l0!\n"
	   "ERROR. Exiting",
           NAME_CURRENT_COMP);
    exit(1);
  }

  if (mx==0 && my==0 && mz==0) {
    fprintf(stderr, "Pol_monitor: %s: NULL vector defined!\n"
	   "ERROR      (mx, my, mz). Exiting",
           NAME_CURRENT_COMP);
    exit(1);
  }

  if ((xwidth<=0) || (yheight <= 0)) {
    fprintf(stderr, "Pol_monitor: %s: Null detection area !\n"
	   "ERROR      (xwidth,yheight). Exiting",
           NAME_CURRENT_COMP);
    exit(1);
  }

  // Initialize variables

  NORM(mx, my, mz);

  for (i=0; i<nL; i++) {

      PolL_N[i] = 0;
      PolL_p[i] = 0;
      PolL_p2[i] = 0;
      HelpArray[i] = 0;
  }
}
#line 12774 "./SE_example2.c"
#undef Lmax
#undef Lmin
#undef mz
#undef my
#undef mx
#undef filename
#undef HelpArray
#undef PolL_p2
#undef PolL_p
#undef PolL_N
#undef nowritefile
#undef restore_neutron
#undef nL
#undef yheight
#undef xwidth
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component flipper3. */
  SIG_MESSAGE("flipper3 (Init)");
#define mccompcurname  flipper3
#define mccompcurtype  Pol_FieldBox
#define mccompcurindex 23
#define fieldFunction mccflipper3_fieldFunction
#define field_parameters mccflipper3_field_parameters
#define const_field_parameters mccflipper3_const_field_parameters
#define xwidth mccflipper3_xwidth
#define yheight mccflipper3_yheight
#define zdepth mccflipper3_zdepth
#define Bx mccflipper3_Bx
#define By mccflipper3_By
#define Bz mccflipper3_Bz
#define filename mccflipper3_filename
#line 55 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Pol_FieldBox.comp"
{
  if (fieldFunction==table_magnetic_field){
      /*initialize the interpolation vertex structure*/
      struct interpolator_struct *interpolator = interpolator_load(filename, 3, 3, NULL);
      interpolator_info(interpolator);
      field_parameters=(void *)interpolator;
  }  else {
      field_parameters = (void *) const_field_parameters;
  }

}
#line 12821 "./SE_example2.c"
#undef filename
#undef Bz
#undef By
#undef Bx
#undef zdepth
#undef yheight
#undef xwidth
#undef const_field_parameters
#undef field_parameters
#undef fieldFunction
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component pl_mon4. */
  SIG_MESSAGE("pl_mon4 (Init)");
#define mccompcurname  pl_mon4
#define mccompcurtype  MeanPolLambda_monitor
#define mccompcurindex 24
#define xwidth mccpl_mon4_xwidth
#define yheight mccpl_mon4_yheight
#define nL mccpl_mon4_nL
#define restore_neutron mccpl_mon4_restore_neutron
#define nowritefile mccpl_mon4_nowritefile
#define PolL_N mccpl_mon4_PolL_N
#define PolL_p mccpl_mon4_PolL_p
#define PolL_p2 mccpl_mon4_PolL_p2
#define HelpArray mccpl_mon4_HelpArray
#define filename mccpl_mon4_filename
#define mx mccpl_mon4_mx
#define my mccpl_mon4_my
#define mz mccpl_mon4_mz
#define Lmin mccpl_mon4_Lmin
#define Lmax mccpl_mon4_Lmax
#line 66 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/MeanPolLambda_monitor.comp"
{
int i;

// Check that input parameteters makes sense

if (Lmax<=Lmin) {
    fprintf(stderr, "Pol_monitor: %s: l1 <= l0!\n"
	   "ERROR. Exiting",
           NAME_CURRENT_COMP);
    exit(1);
  }

  if (mx==0 && my==0 && mz==0) {
    fprintf(stderr, "Pol_monitor: %s: NULL vector defined!\n"
	   "ERROR      (mx, my, mz). Exiting",
           NAME_CURRENT_COMP);
    exit(1);
  }

  if ((xwidth<=0) || (yheight <= 0)) {
    fprintf(stderr, "Pol_monitor: %s: Null detection area !\n"
	   "ERROR      (xwidth,yheight). Exiting",
           NAME_CURRENT_COMP);
    exit(1);
  }

  // Initialize variables

  NORM(mx, my, mz);

  for (i=0; i<nL; i++) {

      PolL_N[i] = 0;
      PolL_p[i] = 0;
      PolL_p2[i] = 0;
      HelpArray[i] = 0;
  }
}
#line 12895 "./SE_example2.c"
#undef Lmax
#undef Lmin
#undef mz
#undef my
#undef mx
#undef filename
#undef HelpArray
#undef PolL_p2
#undef PolL_p
#undef PolL_N
#undef nowritefile
#undef restore_neutron
#undef nL
#undef yheight
#undef xwidth
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component analyzer_arm. */
  SIG_MESSAGE("analyzer_arm (Init)");

  /* Initializations for component psd_ba. */
  SIG_MESSAGE("psd_ba (Init)");
#define mccompcurname  psd_ba
#define mccompcurtype  PSD_monitor
#define mccompcurindex 26
#define PSD_N mccpsd_ba_PSD_N
#define PSD_p mccpsd_ba_PSD_p
#define PSD_p2 mccpsd_ba_PSD_p2
#define nx mccpsd_ba_nx
#define ny mccpsd_ba_ny
#define filename mccpsd_ba_filename
#define xmin mccpsd_ba_xmin
#define xmax mccpsd_ba_xmax
#define ymin mccpsd_ba_ymin
#define ymax mccpsd_ba_ymax
#define xwidth mccpsd_ba_xwidth
#define yheight mccpsd_ba_yheight
#define restore_neutron mccpsd_ba_restore_neutron
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
#line 12961 "./SE_example2.c"
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

  /* Initializations for component pl_mon_ba. */
  SIG_MESSAGE("pl_mon_ba (Init)");
#define mccompcurname  pl_mon_ba
#define mccompcurtype  MeanPolLambda_monitor
#define mccompcurindex 27
#define xwidth mccpl_mon_ba_xwidth
#define yheight mccpl_mon_ba_yheight
#define nL mccpl_mon_ba_nL
#define restore_neutron mccpl_mon_ba_restore_neutron
#define nowritefile mccpl_mon_ba_nowritefile
#define PolL_N mccpl_mon_ba_PolL_N
#define PolL_p mccpl_mon_ba_PolL_p
#define PolL_p2 mccpl_mon_ba_PolL_p2
#define HelpArray mccpl_mon_ba_HelpArray
#define filename mccpl_mon_ba_filename
#define mx mccpl_mon_ba_mx
#define my mccpl_mon_ba_my
#define mz mccpl_mon_ba_mz
#define Lmin mccpl_mon_ba_Lmin
#define Lmax mccpl_mon_ba_Lmax
#line 66 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/MeanPolLambda_monitor.comp"
{
int i;

// Check that input parameteters makes sense

if (Lmax<=Lmin) {
    fprintf(stderr, "Pol_monitor: %s: l1 <= l0!\n"
	   "ERROR. Exiting",
           NAME_CURRENT_COMP);
    exit(1);
  }

  if (mx==0 && my==0 && mz==0) {
    fprintf(stderr, "Pol_monitor: %s: NULL vector defined!\n"
	   "ERROR      (mx, my, mz). Exiting",
           NAME_CURRENT_COMP);
    exit(1);
  }

  if ((xwidth<=0) || (yheight <= 0)) {
    fprintf(stderr, "Pol_monitor: %s: Null detection area !\n"
	   "ERROR      (xwidth,yheight). Exiting",
           NAME_CURRENT_COMP);
    exit(1);
  }

  // Initialize variables

  NORM(mx, my, mz);

  for (i=0; i<nL; i++) {

      PolL_N[i] = 0;
      PolL_p[i] = 0;
      PolL_p2[i] = 0;
      HelpArray[i] = 0;
  }
}
#line 13038 "./SE_example2.c"
#undef Lmax
#undef Lmin
#undef mz
#undef my
#undef mx
#undef filename
#undef HelpArray
#undef PolL_p2
#undef PolL_p
#undef PolL_N
#undef nowritefile
#undef restore_neutron
#undef nL
#undef yheight
#undef xwidth
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component analyzer. */
  SIG_MESSAGE("analyzer (Init)");
#define mccompcurname  analyzer
#define mccompcurtype  Pol_mirror
#define mccompcurindex 28
#define rUpFunc mccanalyzer_rUpFunc
#define rDownFunc mccanalyzer_rDownFunc
#define rUpPar mccanalyzer_rUpPar
#define rDownPar mccanalyzer_rDownPar
#define useTables mccanalyzer_useTables
#define zwidth mccanalyzer_zwidth
#define yheight mccanalyzer_yheight
#define rUpParPtr mccanalyzer_rUpParPtr
#define rDownParPtr mccanalyzer_rDownParPtr
#define p_reflect mccanalyzer_p_reflect
#line 105 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Pol_mirror.comp"
{
#if (useTables)
  rUpParPtr   = (t_Table*) malloc(sizeof(t_Table));
  rDownParPtr = (t_Table*) malloc(sizeof(t_Table));
  if (Table_Read(rUpParPtr, rUpPar, 1) <= 0) {
    fprintf(stderr,"Pol_mirror: %s: can not read file %s\n",
            NAME_CURRENT_COMP, rUpPar);
    exit(1);
  }
  if (Table_Read(rDownParPtr, rDownPar, 1) <= 0) {
    fprintf(stderr,"Pol_mirror: %s: can not read file %s\n",
            NAME_CURRENT_COMP, rDownPar);
    exit(1);
  }
#endif

  if ((zwidth<=0) || (yheight <= 0)) {
    fprintf(stderr, "Pol_mirror: %s: NULL or negative length scale!\n"
            "ERROR      (zwidth=%g,yheight=%g). Exiting\n",
            NAME_CURRENT_COMP, zwidth, yheight);
    exit(1);
  }

}
#line 13098 "./SE_example2.c"
#undef p_reflect
#undef rDownParPtr
#undef rUpParPtr
#undef yheight
#undef zwidth
#undef useTables
#undef rDownPar
#undef rUpPar
#undef rDownFunc
#undef rUpFunc
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component detector_slit. */
  SIG_MESSAGE("detector_slit (Init)");
#define mccompcurname  detector_slit
#define mccompcurtype  Slit
#define mccompcurindex 29
#define xmin mccdetector_slit_xmin
#define xmax mccdetector_slit_xmax
#define ymin mccdetector_slit_ymin
#define ymax mccdetector_slit_ymax
#define radius mccdetector_slit_radius
#define xwidth mccdetector_slit_xwidth
#define yheight mccdetector_slit_yheight
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
#line 13145 "./SE_example2.c"
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

  /* Initializations for component pl_mon5. */
  SIG_MESSAGE("pl_mon5 (Init)");
#define mccompcurname  pl_mon5
#define mccompcurtype  MeanPolLambda_monitor
#define mccompcurindex 30
#define xwidth mccpl_mon5_xwidth
#define yheight mccpl_mon5_yheight
#define nL mccpl_mon5_nL
#define restore_neutron mccpl_mon5_restore_neutron
#define nowritefile mccpl_mon5_nowritefile
#define PolL_N mccpl_mon5_PolL_N
#define PolL_p mccpl_mon5_PolL_p
#define PolL_p2 mccpl_mon5_PolL_p2
#define HelpArray mccpl_mon5_HelpArray
#define filename mccpl_mon5_filename
#define mx mccpl_mon5_mx
#define my mccpl_mon5_my
#define mz mccpl_mon5_mz
#define Lmin mccpl_mon5_Lmin
#define Lmax mccpl_mon5_Lmax
#line 66 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/MeanPolLambda_monitor.comp"
{
int i;

// Check that input parameteters makes sense

if (Lmax<=Lmin) {
    fprintf(stderr, "Pol_monitor: %s: l1 <= l0!\n"
	   "ERROR. Exiting",
           NAME_CURRENT_COMP);
    exit(1);
  }

  if (mx==0 && my==0 && mz==0) {
    fprintf(stderr, "Pol_monitor: %s: NULL vector defined!\n"
	   "ERROR      (mx, my, mz). Exiting",
           NAME_CURRENT_COMP);
    exit(1);
  }

  if ((xwidth<=0) || (yheight <= 0)) {
    fprintf(stderr, "Pol_monitor: %s: Null detection area !\n"
	   "ERROR      (xwidth,yheight). Exiting",
           NAME_CURRENT_COMP);
    exit(1);
  }

  // Initialize variables

  NORM(mx, my, mz);

  for (i=0; i<nL; i++) {

      PolL_N[i] = 0;
      PolL_p[i] = 0;
      PolL_p2[i] = 0;
      HelpArray[i] = 0;
  }
}
#line 13216 "./SE_example2.c"
#undef Lmax
#undef Lmin
#undef mz
#undef my
#undef mx
#undef filename
#undef HelpArray
#undef PolL_p2
#undef PolL_p
#undef PolL_N
#undef nowritefile
#undef restore_neutron
#undef nL
#undef yheight
#undef xwidth
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component detector. */
  SIG_MESSAGE("detector (Init)");
#define mccompcurname  detector
#define mccompcurtype  PSD_monitor
#define mccompcurindex 31
#define PSD_N mccdetector_PSD_N
#define PSD_p mccdetector_PSD_p
#define PSD_p2 mccdetector_PSD_p2
#define nx mccdetector_nx
#define ny mccdetector_ny
#define filename mccdetector_filename
#define xmin mccdetector_xmin
#define xmax mccdetector_xmax
#define ymin mccdetector_ymin
#define ymax mccdetector_ymax
#define xwidth mccdetector_xwidth
#define yheight mccdetector_yheight
#define restore_neutron mccdetector_restore_neutron
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
#line 13279 "./SE_example2.c"
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
#line 13450 "./SE_example2.c"
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
#line 13621 "./SE_example2.c"
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

  /* TRACE Component lam_mon [3] */
  mccoordschange(mcposrlam_mon, mcrotrlam_mon,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component lam_mon (without coords transformations) */
  mcJumpTrace_lam_mon:
  SIG_MESSAGE("lam_mon (Trace)");
  mcDEBUG_COMP("lam_mon")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbComplam_mon
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
#define mccompcurname  lam_mon
#define mccompcurtype  L_monitor
#define mccompcurindex 3
#define nL mcclam_mon_nL
#define L_N mcclam_mon_L_N
#define L_p mcclam_mon_L_p
#define L_p2 mcclam_mon_L_p2
{   /* Declarations of lam_mon=L_monitor() SETTING parameters. */
char* filename = mcclam_mon_filename;
MCNUM xmin = mcclam_mon_xmin;
MCNUM xmax = mcclam_mon_xmax;
MCNUM ymin = mcclam_mon_ymin;
MCNUM ymax = mcclam_mon_ymax;
MCNUM xwidth = mcclam_mon_xwidth;
MCNUM yheight = mcclam_mon_yheight;
MCNUM Lmin = mcclam_mon_Lmin;
MCNUM Lmax = mcclam_mon_Lmax;
MCNUM restore_neutron = mcclam_mon_restore_neutron;
int nowritefile = mcclam_mon_nowritefile;
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
#line 13767 "./SE_example2.c"
}   /* End of lam_mon=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbComplam_mon:
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

  /* TRACE Component polariser_pt [4] */
  mccoordschange(mcposrpolariser_pt, mcrotrpolariser_pt,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component polariser_pt (without coords transformations) */
  mcJumpTrace_polariser_pt:
  SIG_MESSAGE("polariser_pt (Trace)");
  mcDEBUG_COMP("polariser_pt")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbComppolariser_pt
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
#define mccompcurname  polariser_pt
#define mccompcurtype  Arm
#define mccompcurindex 4
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbComppolariser_pt:
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

  /* TRACE Component hdiv_mon [5] */
  mccoordschange(mcposrhdiv_mon, mcrotrhdiv_mon,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component hdiv_mon (without coords transformations) */
  mcJumpTrace_hdiv_mon:
  SIG_MESSAGE("hdiv_mon (Trace)");
  mcDEBUG_COMP("hdiv_mon")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbComphdiv_mon
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
#define mccompcurname  hdiv_mon
#define mccompcurtype  Hdiv_monitor
#define mccompcurindex 5
#define nh mcchdiv_mon_nh
#define Div_N mcchdiv_mon_Div_N
#define Div_p mcchdiv_mon_Div_p
#define Div_p2 mcchdiv_mon_Div_p2
{   /* Declarations of hdiv_mon=Hdiv_monitor() SETTING parameters. */
char* filename = mcchdiv_mon_filename;
MCNUM xmin = mcchdiv_mon_xmin;
MCNUM xmax = mcchdiv_mon_xmax;
MCNUM ymin = mcchdiv_mon_ymin;
MCNUM ymax = mcchdiv_mon_ymax;
MCNUM xwidth = mcchdiv_mon_xwidth;
MCNUM yheight = mcchdiv_mon_yheight;
MCNUM h_maxdiv = mcchdiv_mon_h_maxdiv;
MCNUM restore_neutron = mcchdiv_mon_restore_neutron;
int nowritefile = mcchdiv_mon_nowritefile;
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
#line 14016 "./SE_example2.c"
}   /* End of hdiv_mon=Hdiv_monitor() SETTING parameter declarations. */
#undef Div_p2
#undef Div_p
#undef Div_N
#undef nh
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbComphdiv_mon:
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

  /* TRACE Component polariser2 [6] */
  mccoordschange(mcposrpolariser2, mcrotrpolariser2,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component polariser2 (without coords transformations) */
  mcJumpTrace_polariser2:
  SIG_MESSAGE("polariser2 (Trace)");
  mcDEBUG_COMP("polariser2")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbComppolariser2
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
#define mccompcurname  polariser2
#define mccompcurtype  Pol_mirror
#define mccompcurindex 6
#define rUpFunc mccpolariser2_rUpFunc
#define rDownFunc mccpolariser2_rDownFunc
#define rUpPar mccpolariser2_rUpPar
#define rDownPar mccpolariser2_rDownPar
#define useTables mccpolariser2_useTables
#define zwidth mccpolariser2_zwidth
#define yheight mccpolariser2_yheight
#define rUpParPtr mccpolariser2_rUpParPtr
#define rDownParPtr mccpolariser2_rDownParPtr
{   /* Declarations of polariser2=Pol_mirror() SETTING parameters. */
MCNUM p_reflect = mccpolariser2_p_reflect;
#line 131 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Pol_mirror.comp"
{
  double Q, Rup, Rdown, FN, FM, refWeight;
  int reflect = -1;
  int isPolarising = 0;

  // propagate to mirror plane
  PROP_X0;

  if (inside_rectangle(z, y, zwidth, yheight)) {/* Intersect the crystal? */

    // calculate scattering vector magnitude
    Q = fabs(2*vx*V2K);

    // calculate reflection probability
    rUpFunc(Q, rUpParPtr, &Rup);
    rDownFunc(Q, rDownParPtr, &Rdown);

    if(Rup != Rdown) {

      isPolarising = 1;
      GetMonoPolFNFM(Rup, Rdown, &FN, &FM);
      GetMonoPolRefProb(FN, FM, sy, &refWeight);
    } else
      refWeight = Rup;

    // check that refWeight is meaningfull
    if (refWeight <  0) ABSORB;
    if (refWeight >  1) refWeight =1 ;

    // find out if neutrons is reflected or transmitted
    if (p_reflect < 0 || p_reflect > 1) { // reflect OR transmit using mirror reflectivity

      if (rand01()<refWeight) //reflect
        reflect = 1;
      else
        reflect = 0;

    } else { // reflect OR transmit using a specified weighting

      if (rand01()<p_reflect) {
        reflect = 1;          // reflect
        p *= refWeight/p_reflect;
      } else {
        reflect = 0;          // transmit
        p *= (1-refWeight)/(1-p_reflect);
      }
    }

    // set outgoing velocity and polarisation
    if (reflect==1) { // reflect: SCATTERED==1 for reflection

      vx = -vx;
      if(isPolarising)
        SetMonoPolRefOut(FN, FM, refWeight, &sx, &sy, &sz);

    } else {         // transmit: SCATTERED==2 for transmission
      if(isPolarising)
        SetMonoPolTransOut(FN, FM, refWeight, &sx, &sy, &sz);
      SCATTER;
    }

    if(isPolarising)
      if(sx*sx+sy*sy+sz*sz>1.000001)
        fprintf(stderr,"Pol_mirror: %s: Warning: polarisation |s| = %g > 1\n",
              NAME_CURRENT_COMP, sx*sx+sy*sy+sz*sz); // check that polarisation is meaningfull

    SCATTER;

  } /* End intersect the mirror */
  else
  {
    /* neutron will miss the mirror, so don't propagate i.e restore it */
    RESTORE_NEUTRON(INDEX_CURRENT_COMP, x, y, z, vx, vy, vz, t, sx, sy, sz, p);
  }
}
#line 14212 "./SE_example2.c"
/* 'polariser2=Pol_mirror()' component instance extend code */
    SIG_MESSAGE("polariser2 (Trace:Extend)");
#line 93 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  if(!SCATTERED) ABSORB;
#line 14217 "./SE_example2.c"
}   /* End of polariser2=Pol_mirror() SETTING parameter declarations. */
#undef rDownParPtr
#undef rUpParPtr
#undef yheight
#undef zwidth
#undef useTables
#undef rDownPar
#undef rUpPar
#undef rDownFunc
#undef rUpFunc
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbComppolariser2:
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

  /* TRACE Component psd_ap [7] */
  mccoordschange(mcposrpsd_ap, mcrotrpsd_ap,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component psd_ap (without coords transformations) */
  mcJumpTrace_psd_ap:
  SIG_MESSAGE("psd_ap (Trace)");
  mcDEBUG_COMP("psd_ap")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbComppsd_ap
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
#define mccompcurname  psd_ap
#define mccompcurtype  PSD_monitor
#define mccompcurindex 7
#define PSD_N mccpsd_ap_PSD_N
#define PSD_p mccpsd_ap_PSD_p
#define PSD_p2 mccpsd_ap_PSD_p2
{   /* Declarations of psd_ap=PSD_monitor() SETTING parameters. */
int nx = mccpsd_ap_nx;
int ny = mccpsd_ap_ny;
char* filename = mccpsd_ap_filename;
MCNUM xmin = mccpsd_ap_xmin;
MCNUM xmax = mccpsd_ap_xmax;
MCNUM ymin = mccpsd_ap_ymin;
MCNUM ymax = mccpsd_ap_ymax;
MCNUM xwidth = mccpsd_ap_xwidth;
MCNUM yheight = mccpsd_ap_yheight;
MCNUM restore_neutron = mccpsd_ap_restore_neutron;
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
#line 14360 "./SE_example2.c"
}   /* End of psd_ap=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbComppsd_ap:
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

  /* TRACE Component pl_mon0 [8] */
  mccoordschange(mcposrpl_mon0, mcrotrpl_mon0,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component pl_mon0 (without coords transformations) */
  mcJumpTrace_pl_mon0:
  SIG_MESSAGE("pl_mon0 (Trace)");
  mcDEBUG_COMP("pl_mon0")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbComppl_mon0
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
#define mccompcurname  pl_mon0
#define mccompcurtype  MeanPolLambda_monitor
#define mccompcurindex 8
#define xwidth mccpl_mon0_xwidth
#define yheight mccpl_mon0_yheight
#define nL mccpl_mon0_nL
#define restore_neutron mccpl_mon0_restore_neutron
#define nowritefile mccpl_mon0_nowritefile
#define PolL_N mccpl_mon0_PolL_N
#define PolL_p mccpl_mon0_PolL_p
#define PolL_p2 mccpl_mon0_PolL_p2
#define HelpArray mccpl_mon0_HelpArray
{   /* Declarations of pl_mon0=MeanPolLambda_monitor() SETTING parameters. */
char* filename = mccpl_mon0_filename;
MCNUM mx = mccpl_mon0_mx;
MCNUM my = mccpl_mon0_my;
MCNUM mz = mccpl_mon0_mz;
MCNUM Lmin = mccpl_mon0_Lmin;
MCNUM Lmax = mccpl_mon0_Lmax;
#line 106 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/MeanPolLambda_monitor.comp"
{
  int i;
  double pol_proj;
  double lambda;

  PROP_Z0;
  lambda = (2*PI/V2K)/sqrt(vx*vx + vy*vy + vz*vz);

  if (inside_rectangle(x, y, xwidth, yheight) &&
      lambda > Lmin && lambda < Lmax) {

    pol_proj = scalar_prod(mx, my, mz, sx, sy, sz);

    i= floor((lambda - Lmin)*nL/(Lmax - Lmin));

    PolL_N[i]++;
    HelpArray[i] += p;
    PolL_p[i]    += pol_proj*p;
    PolL_p2[i]   += pol_proj*pol_proj*p;

    SCATTER;
  }
  if (restore_neutron) {
    RESTORE_NEUTRON(INDEX_CURRENT_COMP, x, y, z, vx, vy, vz, t, sx, sy, sz, p);
  }
}
#line 14511 "./SE_example2.c"
}   /* End of pl_mon0=MeanPolLambda_monitor() SETTING parameter declarations. */
#undef HelpArray
#undef PolL_p2
#undef PolL_p
#undef PolL_N
#undef nowritefile
#undef restore_neutron
#undef nL
#undef yheight
#undef xwidth
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbComppl_mon0:
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

  /* TRACE Component flipper1 [9] */
  mccoordschange(mcposrflipper1, mcrotrflipper1,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component flipper1 (without coords transformations) */
  mcJumpTrace_flipper1:
  SIG_MESSAGE("flipper1 (Trace)");
  mcDEBUG_COMP("flipper1")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbCompflipper1
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
#define mccompcurname  flipper1
#define mccompcurtype  Pol_FieldBox
#define mccompcurindex 9
#define fieldFunction mccflipper1_fieldFunction
#define field_parameters mccflipper1_field_parameters
#define const_field_parameters mccflipper1_const_field_parameters
{   /* Declarations of flipper1=Pol_FieldBox() SETTING parameters. */
MCNUM xwidth = mccflipper1_xwidth;
MCNUM yheight = mccflipper1_yheight;
MCNUM zdepth = mccflipper1_zdepth;
MCNUM Bx = mccflipper1_Bx;
MCNUM By = mccflipper1_By;
MCNUM Bz = mccflipper1_Bz;
char* filename = mccflipper1_filename;
#line 70 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Pol_FieldBox.comp"
{
    int hit;
    double t0,t1;
    if (hit=box_intersect(&t0,&t1,x,y,z,vx,vy,vz,xwidth,yheight,zdepth)){
        if(t0>0) PROP_DT(t0);
        ((double *)field_parameters)[0]=Bx;
        ((double *)field_parameters)[1]=By;
        ((double *)field_parameters)[2]=Bz;

        mcmagnet_push(fieldFunction,&(ROT_A_CURRENT_COMP),&(POS_A_CURRENT_COMP),0,field_parameters);
        if(t1-t0>0) PROP_DT(t1-t0);

        mcmagnet_pop();
    }

}
#line 14653 "./SE_example2.c"
}   /* End of flipper1=Pol_FieldBox() SETTING parameter declarations. */
#undef const_field_parameters
#undef field_parameters
#undef fieldFunction
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompflipper1:
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

  /* TRACE Component pl_mon1 [10] */
  mccoordschange(mcposrpl_mon1, mcrotrpl_mon1,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component pl_mon1 (without coords transformations) */
  mcJumpTrace_pl_mon1:
  SIG_MESSAGE("pl_mon1 (Trace)");
  mcDEBUG_COMP("pl_mon1")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbComppl_mon1
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
#define mccompcurname  pl_mon1
#define mccompcurtype  MeanPolLambda_monitor
#define mccompcurindex 10
#define xwidth mccpl_mon1_xwidth
#define yheight mccpl_mon1_yheight
#define nL mccpl_mon1_nL
#define restore_neutron mccpl_mon1_restore_neutron
#define nowritefile mccpl_mon1_nowritefile
#define PolL_N mccpl_mon1_PolL_N
#define PolL_p mccpl_mon1_PolL_p
#define PolL_p2 mccpl_mon1_PolL_p2
#define HelpArray mccpl_mon1_HelpArray
{   /* Declarations of pl_mon1=MeanPolLambda_monitor() SETTING parameters. */
char* filename = mccpl_mon1_filename;
MCNUM mx = mccpl_mon1_mx;
MCNUM my = mccpl_mon1_my;
MCNUM mz = mccpl_mon1_mz;
MCNUM Lmin = mccpl_mon1_Lmin;
MCNUM Lmax = mccpl_mon1_Lmax;
#line 106 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/MeanPolLambda_monitor.comp"
{
  int i;
  double pol_proj;
  double lambda;

  PROP_Z0;
  lambda = (2*PI/V2K)/sqrt(vx*vx + vy*vy + vz*vz);

  if (inside_rectangle(x, y, xwidth, yheight) &&
      lambda > Lmin && lambda < Lmax) {

    pol_proj = scalar_prod(mx, my, mz, sx, sy, sz);

    i= floor((lambda - Lmin)*nL/(Lmax - Lmin));

    PolL_N[i]++;
    HelpArray[i] += p;
    PolL_p[i]    += pol_proj*p;
    PolL_p2[i]   += pol_proj*pol_proj*p;

    SCATTER;
  }
  if (restore_neutron) {
    RESTORE_NEUTRON(INDEX_CURRENT_COMP, x, y, z, vx, vy, vz, t, sx, sy, sz, p);
  }
}
#line 14804 "./SE_example2.c"
}   /* End of pl_mon1=MeanPolLambda_monitor() SETTING parameter declarations. */
#undef HelpArray
#undef PolL_p2
#undef PolL_p
#undef PolL_N
#undef nowritefile
#undef restore_neutron
#undef nL
#undef yheight
#undef xwidth
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbComppl_mon1:
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

  /* TRACE Component guide_field_1 [11] */
  mccoordschange(mcposrguide_field_1, mcrotrguide_field_1,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component guide_field_1 (without coords transformations) */
  mcJumpTrace_guide_field_1:
  SIG_MESSAGE("guide_field_1 (Trace)");
  mcDEBUG_COMP("guide_field_1")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbCompguide_field_1
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
#define mccompcurname  guide_field_1
#define mccompcurtype  Pol_FieldBox
#define mccompcurindex 11
#define fieldFunction mccguide_field_1_fieldFunction
#define field_parameters mccguide_field_1_field_parameters
#define const_field_parameters mccguide_field_1_const_field_parameters
{   /* Declarations of guide_field_1=Pol_FieldBox() SETTING parameters. */
MCNUM xwidth = mccguide_field_1_xwidth;
MCNUM yheight = mccguide_field_1_yheight;
MCNUM zdepth = mccguide_field_1_zdepth;
MCNUM Bx = mccguide_field_1_Bx;
MCNUM By = mccguide_field_1_By;
MCNUM Bz = mccguide_field_1_Bz;
char* filename = mccguide_field_1_filename;
#line 70 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Pol_FieldBox.comp"
{
    int hit;
    double t0,t1;
    if (hit=box_intersect(&t0,&t1,x,y,z,vx,vy,vz,xwidth,yheight,zdepth)){
        if(t0>0) PROP_DT(t0);
        ((double *)field_parameters)[0]=Bx;
        ((double *)field_parameters)[1]=By;
        ((double *)field_parameters)[2]=Bz;

        mcmagnet_push(fieldFunction,&(ROT_A_CURRENT_COMP),&(POS_A_CURRENT_COMP),0,field_parameters);
        if(t1-t0>0) PROP_DT(t1-t0);

        mcmagnet_pop();
    }

}
#line 14946 "./SE_example2.c"
}   /* End of guide_field_1=Pol_FieldBox() SETTING parameter declarations. */
#undef const_field_parameters
#undef field_parameters
#undef fieldFunction
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompguide_field_1:
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

  /* TRACE Component guide_field_1_end [12] */
  mccoordschange(mcposrguide_field_1_end, mcrotrguide_field_1_end,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component guide_field_1_end (without coords transformations) */
  mcJumpTrace_guide_field_1_end:
  SIG_MESSAGE("guide_field_1_end (Trace)");
  mcDEBUG_COMP("guide_field_1_end")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbCompguide_field_1_end
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
#define mccompcurname  guide_field_1_end
#define mccompcurtype  Slit
#define mccompcurindex 12
{   /* Declarations of guide_field_1_end=Slit() SETTING parameters. */
MCNUM xmin = mccguide_field_1_end_xmin;
MCNUM xmax = mccguide_field_1_end_xmax;
MCNUM ymin = mccguide_field_1_end_ymin;
MCNUM ymax = mccguide_field_1_end_ymax;
MCNUM radius = mccguide_field_1_end_radius;
MCNUM xwidth = mccguide_field_1_end_xwidth;
MCNUM yheight = mccguide_field_1_end_yheight;
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
#line 15073 "./SE_example2.c"
}   /* End of guide_field_1_end=Slit() SETTING parameter declarations. */
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompguide_field_1_end:
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

  /* TRACE Component pl_mon2 [13] */
  mccoordschange(mcposrpl_mon2, mcrotrpl_mon2,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component pl_mon2 (without coords transformations) */
  mcJumpTrace_pl_mon2:
  SIG_MESSAGE("pl_mon2 (Trace)");
  mcDEBUG_COMP("pl_mon2")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbComppl_mon2
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
#define mccompcurname  pl_mon2
#define mccompcurtype  MeanPolLambda_monitor
#define mccompcurindex 13
#define xwidth mccpl_mon2_xwidth
#define yheight mccpl_mon2_yheight
#define nL mccpl_mon2_nL
#define restore_neutron mccpl_mon2_restore_neutron
#define nowritefile mccpl_mon2_nowritefile
#define PolL_N mccpl_mon2_PolL_N
#define PolL_p mccpl_mon2_PolL_p
#define PolL_p2 mccpl_mon2_PolL_p2
#define HelpArray mccpl_mon2_HelpArray
{   /* Declarations of pl_mon2=MeanPolLambda_monitor() SETTING parameters. */
char* filename = mccpl_mon2_filename;
MCNUM mx = mccpl_mon2_mx;
MCNUM my = mccpl_mon2_my;
MCNUM mz = mccpl_mon2_mz;
MCNUM Lmin = mccpl_mon2_Lmin;
MCNUM Lmax = mccpl_mon2_Lmax;
#line 106 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/MeanPolLambda_monitor.comp"
{
  int i;
  double pol_proj;
  double lambda;

  PROP_Z0;
  lambda = (2*PI/V2K)/sqrt(vx*vx + vy*vy + vz*vz);

  if (inside_rectangle(x, y, xwidth, yheight) &&
      lambda > Lmin && lambda < Lmax) {

    pol_proj = scalar_prod(mx, my, mz, sx, sy, sz);

    i= floor((lambda - Lmin)*nL/(Lmax - Lmin));

    PolL_N[i]++;
    HelpArray[i] += p;
    PolL_p[i]    += pol_proj*p;
    PolL_p2[i]   += pol_proj*pol_proj*p;

    SCATTER;
  }
  if (restore_neutron) {
    RESTORE_NEUTRON(INDEX_CURRENT_COMP, x, y, z, vx, vy, vz, t, sx, sy, sz, p);
  }
}
#line 15221 "./SE_example2.c"
}   /* End of pl_mon2=MeanPolLambda_monitor() SETTING parameter declarations. */
#undef HelpArray
#undef PolL_p2
#undef PolL_p
#undef PolL_N
#undef nowritefile
#undef restore_neutron
#undef nL
#undef yheight
#undef xwidth
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbComppl_mon2:
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

  /* TRACE Component flipper2_1 [14] */
  mccoordschange(mcposrflipper2_1, mcrotrflipper2_1,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component flipper2_1 (without coords transformations) */
  mcJumpTrace_flipper2_1:
  SIG_MESSAGE("flipper2_1 (Trace)");
  mcDEBUG_COMP("flipper2_1")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbCompflipper2_1
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
#define mccompcurname  flipper2_1
#define mccompcurtype  Pol_FieldBox
#define mccompcurindex 14
#define fieldFunction mccflipper2_1_fieldFunction
#define field_parameters mccflipper2_1_field_parameters
#define const_field_parameters mccflipper2_1_const_field_parameters
{   /* Declarations of flipper2_1=Pol_FieldBox() SETTING parameters. */
MCNUM xwidth = mccflipper2_1_xwidth;
MCNUM yheight = mccflipper2_1_yheight;
MCNUM zdepth = mccflipper2_1_zdepth;
MCNUM Bx = mccflipper2_1_Bx;
MCNUM By = mccflipper2_1_By;
MCNUM Bz = mccflipper2_1_Bz;
char* filename = mccflipper2_1_filename;
#line 70 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Pol_FieldBox.comp"
{
    int hit;
    double t0,t1;
    if (hit=box_intersect(&t0,&t1,x,y,z,vx,vy,vz,xwidth,yheight,zdepth)){
        if(t0>0) PROP_DT(t0);
        ((double *)field_parameters)[0]=Bx;
        ((double *)field_parameters)[1]=By;
        ((double *)field_parameters)[2]=Bz;

        mcmagnet_push(fieldFunction,&(ROT_A_CURRENT_COMP),&(POS_A_CURRENT_COMP),0,field_parameters);
        if(t1-t0>0) PROP_DT(t1-t0);

        mcmagnet_pop();
    }

}
#line 15363 "./SE_example2.c"
}   /* End of flipper2_1=Pol_FieldBox() SETTING parameter declarations. */
#undef const_field_parameters
#undef field_parameters
#undef fieldFunction
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompflipper2_1:
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

  /* TRACE Component flipper2_2 [15] */
  mccoordschange(mcposrflipper2_2, mcrotrflipper2_2,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component flipper2_2 (without coords transformations) */
  mcJumpTrace_flipper2_2:
  SIG_MESSAGE("flipper2_2 (Trace)");
  mcDEBUG_COMP("flipper2_2")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbCompflipper2_2
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
#define mccompcurname  flipper2_2
#define mccompcurtype  Pol_FieldBox
#define mccompcurindex 15
#define fieldFunction mccflipper2_2_fieldFunction
#define field_parameters mccflipper2_2_field_parameters
#define const_field_parameters mccflipper2_2_const_field_parameters
{   /* Declarations of flipper2_2=Pol_FieldBox() SETTING parameters. */
MCNUM xwidth = mccflipper2_2_xwidth;
MCNUM yheight = mccflipper2_2_yheight;
MCNUM zdepth = mccflipper2_2_zdepth;
MCNUM Bx = mccflipper2_2_Bx;
MCNUM By = mccflipper2_2_By;
MCNUM Bz = mccflipper2_2_Bz;
char* filename = mccflipper2_2_filename;
#line 70 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Pol_FieldBox.comp"
{
    int hit;
    double t0,t1;
    if (hit=box_intersect(&t0,&t1,x,y,z,vx,vy,vz,xwidth,yheight,zdepth)){
        if(t0>0) PROP_DT(t0);
        ((double *)field_parameters)[0]=Bx;
        ((double *)field_parameters)[1]=By;
        ((double *)field_parameters)[2]=Bz;

        mcmagnet_push(fieldFunction,&(ROT_A_CURRENT_COMP),&(POS_A_CURRENT_COMP),0,field_parameters);
        if(t1-t0>0) PROP_DT(t1-t0);

        mcmagnet_pop();
    }

}
#line 15499 "./SE_example2.c"
}   /* End of flipper2_2=Pol_FieldBox() SETTING parameter declarations. */
#undef const_field_parameters
#undef field_parameters
#undef fieldFunction
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompflipper2_2:
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

  /* TRACE Component flipper_2_2_end [16] */
  mccoordschange(mcposrflipper_2_2_end, mcrotrflipper_2_2_end,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component flipper_2_2_end (without coords transformations) */
  mcJumpTrace_flipper_2_2_end:
  SIG_MESSAGE("flipper_2_2_end (Trace)");
  mcDEBUG_COMP("flipper_2_2_end")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbCompflipper_2_2_end
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
#define mccompcurname  flipper_2_2_end
#define mccompcurtype  Slit
#define mccompcurindex 16
{   /* Declarations of flipper_2_2_end=Slit() SETTING parameters. */
MCNUM xmin = mccflipper_2_2_end_xmin;
MCNUM xmax = mccflipper_2_2_end_xmax;
MCNUM ymin = mccflipper_2_2_end_ymin;
MCNUM ymax = mccflipper_2_2_end_ymax;
MCNUM radius = mccflipper_2_2_end_radius;
MCNUM xwidth = mccflipper_2_2_end_xwidth;
MCNUM yheight = mccflipper_2_2_end_yheight;
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
#line 15626 "./SE_example2.c"
}   /* End of flipper_2_2_end=Slit() SETTING parameter declarations. */
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompflipper_2_2_end:
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

  /* TRACE Component sample_mnt [17] */
  mccoordschange(mcposrsample_mnt, mcrotrsample_mnt,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component sample_mnt (without coords transformations) */
  mcJumpTrace_sample_mnt:
  SIG_MESSAGE("sample_mnt (Trace)");
  mcDEBUG_COMP("sample_mnt")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbCompsample_mnt
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
#define mccompcurname  sample_mnt
#define mccompcurtype  Arm
#define mccompcurindex 17
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompsample_mnt:
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

  /* TRACE Component pl_mon3 [18] */
  mccoordschange(mcposrpl_mon3, mcrotrpl_mon3,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component pl_mon3 (without coords transformations) */
  mcJumpTrace_pl_mon3:
  SIG_MESSAGE("pl_mon3 (Trace)");
  mcDEBUG_COMP("pl_mon3")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbComppl_mon3
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
#define mccompcurname  pl_mon3
#define mccompcurtype  MeanPolLambda_monitor
#define mccompcurindex 18
#define xwidth mccpl_mon3_xwidth
#define yheight mccpl_mon3_yheight
#define nL mccpl_mon3_nL
#define restore_neutron mccpl_mon3_restore_neutron
#define nowritefile mccpl_mon3_nowritefile
#define PolL_N mccpl_mon3_PolL_N
#define PolL_p mccpl_mon3_PolL_p
#define PolL_p2 mccpl_mon3_PolL_p2
#define HelpArray mccpl_mon3_HelpArray
{   /* Declarations of pl_mon3=MeanPolLambda_monitor() SETTING parameters. */
char* filename = mccpl_mon3_filename;
MCNUM mx = mccpl_mon3_mx;
MCNUM my = mccpl_mon3_my;
MCNUM mz = mccpl_mon3_mz;
MCNUM Lmin = mccpl_mon3_Lmin;
MCNUM Lmax = mccpl_mon3_Lmax;
#line 106 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/MeanPolLambda_monitor.comp"
{
  int i;
  double pol_proj;
  double lambda;

  PROP_Z0;
  lambda = (2*PI/V2K)/sqrt(vx*vx + vy*vy + vz*vz);

  if (inside_rectangle(x, y, xwidth, yheight) &&
      lambda > Lmin && lambda < Lmax) {

    pol_proj = scalar_prod(mx, my, mz, sx, sy, sz);

    i= floor((lambda - Lmin)*nL/(Lmax - Lmin));

    PolL_N[i]++;
    HelpArray[i] += p;
    PolL_p[i]    += pol_proj*p;
    PolL_p2[i]   += pol_proj*pol_proj*p;

    SCATTER;
  }
  if (restore_neutron) {
    RESTORE_NEUTRON(INDEX_CURRENT_COMP, x, y, z, vx, vy, vz, t, sx, sy, sz, p);
  }
}
#line 15877 "./SE_example2.c"
}   /* End of pl_mon3=MeanPolLambda_monitor() SETTING parameter declarations. */
#undef HelpArray
#undef PolL_p2
#undef PolL_p
#undef PolL_N
#undef nowritefile
#undef restore_neutron
#undef nL
#undef yheight
#undef xwidth
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbComppl_mon3:
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

  /* TRACE Component sample0 [19] */
  mccoordschange(mcposrsample0, mcrotrsample0,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component sample0 (without coords transformations) */
  mcJumpTrace_sample0:
  SIG_MESSAGE("sample0 (Trace)");
  mcDEBUG_COMP("sample0")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbCompsample0
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
#define mccompcurname  sample0
#define mccompcurtype  V_sample
#define mccompcurindex 19
#define VarsV mccsample0_VarsV
{   /* Declarations of sample0=V_sample() SETTING parameters. */
MCNUM radius = mccsample0_radius;
MCNUM thickness = mccsample0_thickness;
MCNUM zdepth = mccsample0_zdepth;
MCNUM Vc = mccsample0_Vc;
MCNUM sigma_abs = mccsample0_sigma_abs;
MCNUM sigma_inc = mccsample0_sigma_inc;
MCNUM radius_i = mccsample0_radius_i;
MCNUM radius_o = mccsample0_radius_o;
MCNUM h = mccsample0_h;
MCNUM focus_r = mccsample0_focus_r;
MCNUM pack = mccsample0_pack;
MCNUM frac = mccsample0_frac;
MCNUM f_QE = mccsample0_f_QE;
MCNUM gamma = mccsample0_gamma;
MCNUM target_x = mccsample0_target_x;
MCNUM target_y = mccsample0_target_y;
MCNUM target_z = mccsample0_target_z;
MCNUM focus_xw = mccsample0_focus_xw;
MCNUM focus_yh = mccsample0_focus_yh;
MCNUM focus_aw = mccsample0_focus_aw;
MCNUM focus_ah = mccsample0_focus_ah;
MCNUM xwidth = mccsample0_xwidth;
MCNUM yheight = mccsample0_yheight;
MCNUM zthick = mccsample0_zthick;
MCNUM rad_sphere = mccsample0_rad_sphere;
MCNUM sig_a = mccsample0_sig_a;
MCNUM sig_i = mccsample0_sig_i;
MCNUM V0 = mccsample0_V0;
int target_index = mccsample0_target_index;
MCNUM multiples = mccsample0_multiples;
/* 'sample0=V_sample()' component instance has conditional execution */
if (( mcipSAMPLE == 0 ))

#line 180 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../obsolete/V_sample.comp"
{
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
}
#line 16164 "./SE_example2.c"
}   /* End of sample0=V_sample() SETTING parameter declarations. */
#undef VarsV
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompsample0:
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

  /* TRACE Component guide_field_2 [20] */
  mccoordschange(mcposrguide_field_2, mcrotrguide_field_2,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component guide_field_2 (without coords transformations) */
  mcJumpTrace_guide_field_2:
  SIG_MESSAGE("guide_field_2 (Trace)");
  mcDEBUG_COMP("guide_field_2")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbCompguide_field_2
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
#define mccompcurname  guide_field_2
#define mccompcurtype  Pol_FieldBox
#define mccompcurindex 20
#define fieldFunction mccguide_field_2_fieldFunction
#define field_parameters mccguide_field_2_field_parameters
#define const_field_parameters mccguide_field_2_const_field_parameters
{   /* Declarations of guide_field_2=Pol_FieldBox() SETTING parameters. */
MCNUM xwidth = mccguide_field_2_xwidth;
MCNUM yheight = mccguide_field_2_yheight;
MCNUM zdepth = mccguide_field_2_zdepth;
MCNUM Bx = mccguide_field_2_Bx;
MCNUM By = mccguide_field_2_By;
MCNUM Bz = mccguide_field_2_Bz;
char* filename = mccguide_field_2_filename;
#line 70 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Pol_FieldBox.comp"
{
    int hit;
    double t0,t1;
    if (hit=box_intersect(&t0,&t1,x,y,z,vx,vy,vz,xwidth,yheight,zdepth)){
        if(t0>0) PROP_DT(t0);
        ((double *)field_parameters)[0]=Bx;
        ((double *)field_parameters)[1]=By;
        ((double *)field_parameters)[2]=Bz;

        mcmagnet_push(fieldFunction,&(ROT_A_CURRENT_COMP),&(POS_A_CURRENT_COMP),0,field_parameters);
        if(t1-t0>0) PROP_DT(t1-t0);

        mcmagnet_pop();
    }

}
#line 16298 "./SE_example2.c"
}   /* End of guide_field_2=Pol_FieldBox() SETTING parameter declarations. */
#undef const_field_parameters
#undef field_parameters
#undef fieldFunction
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompguide_field_2:
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

  /* TRACE Component guide_field_2_end [21] */
  mccoordschange(mcposrguide_field_2_end, mcrotrguide_field_2_end,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component guide_field_2_end (without coords transformations) */
  mcJumpTrace_guide_field_2_end:
  SIG_MESSAGE("guide_field_2_end (Trace)");
  mcDEBUG_COMP("guide_field_2_end")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbCompguide_field_2_end
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
#define mccompcurname  guide_field_2_end
#define mccompcurtype  Slit
#define mccompcurindex 21
{   /* Declarations of guide_field_2_end=Slit() SETTING parameters. */
MCNUM xmin = mccguide_field_2_end_xmin;
MCNUM xmax = mccguide_field_2_end_xmax;
MCNUM ymin = mccguide_field_2_end_ymin;
MCNUM ymax = mccguide_field_2_end_ymax;
MCNUM radius = mccguide_field_2_end_radius;
MCNUM xwidth = mccguide_field_2_end_xwidth;
MCNUM yheight = mccguide_field_2_end_yheight;
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
#line 16425 "./SE_example2.c"
}   /* End of guide_field_2_end=Slit() SETTING parameter declarations. */
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompguide_field_2_end:
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

  /* TRACE Component pl_mon40 [22] */
  mccoordschange(mcposrpl_mon40, mcrotrpl_mon40,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component pl_mon40 (without coords transformations) */
  mcJumpTrace_pl_mon40:
  SIG_MESSAGE("pl_mon40 (Trace)");
  mcDEBUG_COMP("pl_mon40")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbComppl_mon40
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
#define mccompcurname  pl_mon40
#define mccompcurtype  MeanPolLambda_monitor
#define mccompcurindex 22
#define xwidth mccpl_mon40_xwidth
#define yheight mccpl_mon40_yheight
#define nL mccpl_mon40_nL
#define restore_neutron mccpl_mon40_restore_neutron
#define nowritefile mccpl_mon40_nowritefile
#define PolL_N mccpl_mon40_PolL_N
#define PolL_p mccpl_mon40_PolL_p
#define PolL_p2 mccpl_mon40_PolL_p2
#define HelpArray mccpl_mon40_HelpArray
{   /* Declarations of pl_mon40=MeanPolLambda_monitor() SETTING parameters. */
char* filename = mccpl_mon40_filename;
MCNUM mx = mccpl_mon40_mx;
MCNUM my = mccpl_mon40_my;
MCNUM mz = mccpl_mon40_mz;
MCNUM Lmin = mccpl_mon40_Lmin;
MCNUM Lmax = mccpl_mon40_Lmax;
#line 106 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/MeanPolLambda_monitor.comp"
{
  int i;
  double pol_proj;
  double lambda;

  PROP_Z0;
  lambda = (2*PI/V2K)/sqrt(vx*vx + vy*vy + vz*vz);

  if (inside_rectangle(x, y, xwidth, yheight) &&
      lambda > Lmin && lambda < Lmax) {

    pol_proj = scalar_prod(mx, my, mz, sx, sy, sz);

    i= floor((lambda - Lmin)*nL/(Lmax - Lmin));

    PolL_N[i]++;
    HelpArray[i] += p;
    PolL_p[i]    += pol_proj*p;
    PolL_p2[i]   += pol_proj*pol_proj*p;

    SCATTER;
  }
  if (restore_neutron) {
    RESTORE_NEUTRON(INDEX_CURRENT_COMP, x, y, z, vx, vy, vz, t, sx, sy, sz, p);
  }
}
#line 16573 "./SE_example2.c"
}   /* End of pl_mon40=MeanPolLambda_monitor() SETTING parameter declarations. */
#undef HelpArray
#undef PolL_p2
#undef PolL_p
#undef PolL_N
#undef nowritefile
#undef restore_neutron
#undef nL
#undef yheight
#undef xwidth
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbComppl_mon40:
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

  /* TRACE Component flipper3 [23] */
  mccoordschange(mcposrflipper3, mcrotrflipper3,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component flipper3 (without coords transformations) */
  mcJumpTrace_flipper3:
  SIG_MESSAGE("flipper3 (Trace)");
  mcDEBUG_COMP("flipper3")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbCompflipper3
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
#define mccompcurname  flipper3
#define mccompcurtype  Pol_FieldBox
#define mccompcurindex 23
#define fieldFunction mccflipper3_fieldFunction
#define field_parameters mccflipper3_field_parameters
#define const_field_parameters mccflipper3_const_field_parameters
{   /* Declarations of flipper3=Pol_FieldBox() SETTING parameters. */
MCNUM xwidth = mccflipper3_xwidth;
MCNUM yheight = mccflipper3_yheight;
MCNUM zdepth = mccflipper3_zdepth;
MCNUM Bx = mccflipper3_Bx;
MCNUM By = mccflipper3_By;
MCNUM Bz = mccflipper3_Bz;
char* filename = mccflipper3_filename;
#line 70 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Pol_FieldBox.comp"
{
    int hit;
    double t0,t1;
    if (hit=box_intersect(&t0,&t1,x,y,z,vx,vy,vz,xwidth,yheight,zdepth)){
        if(t0>0) PROP_DT(t0);
        ((double *)field_parameters)[0]=Bx;
        ((double *)field_parameters)[1]=By;
        ((double *)field_parameters)[2]=Bz;

        mcmagnet_push(fieldFunction,&(ROT_A_CURRENT_COMP),&(POS_A_CURRENT_COMP),0,field_parameters);
        if(t1-t0>0) PROP_DT(t1-t0);

        mcmagnet_pop();
    }

}
#line 16715 "./SE_example2.c"
}   /* End of flipper3=Pol_FieldBox() SETTING parameter declarations. */
#undef const_field_parameters
#undef field_parameters
#undef fieldFunction
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompflipper3:
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

  /* TRACE Component pl_mon4 [24] */
  mccoordschange(mcposrpl_mon4, mcrotrpl_mon4,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component pl_mon4 (without coords transformations) */
  mcJumpTrace_pl_mon4:
  SIG_MESSAGE("pl_mon4 (Trace)");
  mcDEBUG_COMP("pl_mon4")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbComppl_mon4
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
#define mccompcurname  pl_mon4
#define mccompcurtype  MeanPolLambda_monitor
#define mccompcurindex 24
#define xwidth mccpl_mon4_xwidth
#define yheight mccpl_mon4_yheight
#define nL mccpl_mon4_nL
#define restore_neutron mccpl_mon4_restore_neutron
#define nowritefile mccpl_mon4_nowritefile
#define PolL_N mccpl_mon4_PolL_N
#define PolL_p mccpl_mon4_PolL_p
#define PolL_p2 mccpl_mon4_PolL_p2
#define HelpArray mccpl_mon4_HelpArray
{   /* Declarations of pl_mon4=MeanPolLambda_monitor() SETTING parameters. */
char* filename = mccpl_mon4_filename;
MCNUM mx = mccpl_mon4_mx;
MCNUM my = mccpl_mon4_my;
MCNUM mz = mccpl_mon4_mz;
MCNUM Lmin = mccpl_mon4_Lmin;
MCNUM Lmax = mccpl_mon4_Lmax;
#line 106 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/MeanPolLambda_monitor.comp"
{
  int i;
  double pol_proj;
  double lambda;

  PROP_Z0;
  lambda = (2*PI/V2K)/sqrt(vx*vx + vy*vy + vz*vz);

  if (inside_rectangle(x, y, xwidth, yheight) &&
      lambda > Lmin && lambda < Lmax) {

    pol_proj = scalar_prod(mx, my, mz, sx, sy, sz);

    i= floor((lambda - Lmin)*nL/(Lmax - Lmin));

    PolL_N[i]++;
    HelpArray[i] += p;
    PolL_p[i]    += pol_proj*p;
    PolL_p2[i]   += pol_proj*pol_proj*p;

    SCATTER;
  }
  if (restore_neutron) {
    RESTORE_NEUTRON(INDEX_CURRENT_COMP, x, y, z, vx, vy, vz, t, sx, sy, sz, p);
  }
}
#line 16866 "./SE_example2.c"
}   /* End of pl_mon4=MeanPolLambda_monitor() SETTING parameter declarations. */
#undef HelpArray
#undef PolL_p2
#undef PolL_p
#undef PolL_N
#undef nowritefile
#undef restore_neutron
#undef nL
#undef yheight
#undef xwidth
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbComppl_mon4:
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

  /* TRACE Component analyzer_arm [25] */
  mccoordschange(mcposranalyzer_arm, mcrotranalyzer_arm,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component analyzer_arm (without coords transformations) */
  mcJumpTrace_analyzer_arm:
  SIG_MESSAGE("analyzer_arm (Trace)");
  mcDEBUG_COMP("analyzer_arm")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbCompanalyzer_arm
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
#define mccompcurname  analyzer_arm
#define mccompcurtype  Arm
#define mccompcurindex 25
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompanalyzer_arm:
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

  /* TRACE Component psd_ba [26] */
  mccoordschange(mcposrpsd_ba, mcrotrpsd_ba,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component psd_ba (without coords transformations) */
  mcJumpTrace_psd_ba:
  SIG_MESSAGE("psd_ba (Trace)");
  mcDEBUG_COMP("psd_ba")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbComppsd_ba
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
#define mccompcurname  psd_ba
#define mccompcurtype  PSD_monitor
#define mccompcurindex 26
#define PSD_N mccpsd_ba_PSD_N
#define PSD_p mccpsd_ba_PSD_p
#define PSD_p2 mccpsd_ba_PSD_p2
{   /* Declarations of psd_ba=PSD_monitor() SETTING parameters. */
int nx = mccpsd_ba_nx;
int ny = mccpsd_ba_ny;
char* filename = mccpsd_ba_filename;
MCNUM xmin = mccpsd_ba_xmin;
MCNUM xmax = mccpsd_ba_xmax;
MCNUM ymin = mccpsd_ba_ymin;
MCNUM ymax = mccpsd_ba_ymax;
MCNUM xwidth = mccpsd_ba_xwidth;
MCNUM yheight = mccpsd_ba_yheight;
MCNUM restore_neutron = mccpsd_ba_restore_neutron;
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
#line 17112 "./SE_example2.c"
}   /* End of psd_ba=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbComppsd_ba:
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

  /* TRACE Component pl_mon_ba [27] */
  mccoordschange(mcposrpl_mon_ba, mcrotrpl_mon_ba,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component pl_mon_ba (without coords transformations) */
  mcJumpTrace_pl_mon_ba:
  SIG_MESSAGE("pl_mon_ba (Trace)");
  mcDEBUG_COMP("pl_mon_ba")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbComppl_mon_ba
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
#define mccompcurname  pl_mon_ba
#define mccompcurtype  MeanPolLambda_monitor
#define mccompcurindex 27
#define xwidth mccpl_mon_ba_xwidth
#define yheight mccpl_mon_ba_yheight
#define nL mccpl_mon_ba_nL
#define restore_neutron mccpl_mon_ba_restore_neutron
#define nowritefile mccpl_mon_ba_nowritefile
#define PolL_N mccpl_mon_ba_PolL_N
#define PolL_p mccpl_mon_ba_PolL_p
#define PolL_p2 mccpl_mon_ba_PolL_p2
#define HelpArray mccpl_mon_ba_HelpArray
{   /* Declarations of pl_mon_ba=MeanPolLambda_monitor() SETTING parameters. */
char* filename = mccpl_mon_ba_filename;
MCNUM mx = mccpl_mon_ba_mx;
MCNUM my = mccpl_mon_ba_my;
MCNUM mz = mccpl_mon_ba_mz;
MCNUM Lmin = mccpl_mon_ba_Lmin;
MCNUM Lmax = mccpl_mon_ba_Lmax;
#line 106 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/MeanPolLambda_monitor.comp"
{
  int i;
  double pol_proj;
  double lambda;

  PROP_Z0;
  lambda = (2*PI/V2K)/sqrt(vx*vx + vy*vy + vz*vz);

  if (inside_rectangle(x, y, xwidth, yheight) &&
      lambda > Lmin && lambda < Lmax) {

    pol_proj = scalar_prod(mx, my, mz, sx, sy, sz);

    i= floor((lambda - Lmin)*nL/(Lmax - Lmin));

    PolL_N[i]++;
    HelpArray[i] += p;
    PolL_p[i]    += pol_proj*p;
    PolL_p2[i]   += pol_proj*pol_proj*p;

    SCATTER;
  }
  if (restore_neutron) {
    RESTORE_NEUTRON(INDEX_CURRENT_COMP, x, y, z, vx, vy, vz, t, sx, sy, sz, p);
  }
}
#line 17263 "./SE_example2.c"
}   /* End of pl_mon_ba=MeanPolLambda_monitor() SETTING parameter declarations. */
#undef HelpArray
#undef PolL_p2
#undef PolL_p
#undef PolL_N
#undef nowritefile
#undef restore_neutron
#undef nL
#undef yheight
#undef xwidth
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbComppl_mon_ba:
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

  /* TRACE Component analyzer [28] */
  mccoordschange(mcposranalyzer, mcrotranalyzer,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component analyzer (without coords transformations) */
  mcJumpTrace_analyzer:
  SIG_MESSAGE("analyzer (Trace)");
  mcDEBUG_COMP("analyzer")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbCompanalyzer
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
#define mccompcurname  analyzer
#define mccompcurtype  Pol_mirror
#define mccompcurindex 28
#define rUpFunc mccanalyzer_rUpFunc
#define rDownFunc mccanalyzer_rDownFunc
#define rUpPar mccanalyzer_rUpPar
#define rDownPar mccanalyzer_rDownPar
#define useTables mccanalyzer_useTables
#define zwidth mccanalyzer_zwidth
#define yheight mccanalyzer_yheight
#define rUpParPtr mccanalyzer_rUpParPtr
#define rDownParPtr mccanalyzer_rDownParPtr
{   /* Declarations of analyzer=Pol_mirror() SETTING parameters. */
MCNUM p_reflect = mccanalyzer_p_reflect;
#line 131 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Pol_mirror.comp"
{
  double Q, Rup, Rdown, FN, FM, refWeight;
  int reflect = -1;
  int isPolarising = 0;

  // propagate to mirror plane
  PROP_X0;

  if (inside_rectangle(z, y, zwidth, yheight)) {/* Intersect the crystal? */

    // calculate scattering vector magnitude
    Q = fabs(2*vx*V2K);

    // calculate reflection probability
    rUpFunc(Q, rUpParPtr, &Rup);
    rDownFunc(Q, rDownParPtr, &Rdown);

    if(Rup != Rdown) {

      isPolarising = 1;
      GetMonoPolFNFM(Rup, Rdown, &FN, &FM);
      GetMonoPolRefProb(FN, FM, sy, &refWeight);
    } else
      refWeight = Rup;

    // check that refWeight is meaningfull
    if (refWeight <  0) ABSORB;
    if (refWeight >  1) refWeight =1 ;

    // find out if neutrons is reflected or transmitted
    if (p_reflect < 0 || p_reflect > 1) { // reflect OR transmit using mirror reflectivity

      if (rand01()<refWeight) //reflect
        reflect = 1;
      else
        reflect = 0;

    } else { // reflect OR transmit using a specified weighting

      if (rand01()<p_reflect) {
        reflect = 1;          // reflect
        p *= refWeight/p_reflect;
      } else {
        reflect = 0;          // transmit
        p *= (1-refWeight)/(1-p_reflect);
      }
    }

    // set outgoing velocity and polarisation
    if (reflect==1) { // reflect: SCATTERED==1 for reflection

      vx = -vx;
      if(isPolarising)
        SetMonoPolRefOut(FN, FM, refWeight, &sx, &sy, &sz);

    } else {         // transmit: SCATTERED==2 for transmission
      if(isPolarising)
        SetMonoPolTransOut(FN, FM, refWeight, &sx, &sy, &sz);
      SCATTER;
    }

    if(isPolarising)
      if(sx*sx+sy*sy+sz*sz>1.000001)
        fprintf(stderr,"Pol_mirror: %s: Warning: polarisation |s| = %g > 1\n",
              NAME_CURRENT_COMP, sx*sx+sy*sy+sz*sz); // check that polarisation is meaningfull

    SCATTER;

  } /* End intersect the mirror */
  else
  {
    /* neutron will miss the mirror, so don't propagate i.e restore it */
    RESTORE_NEUTRON(INDEX_CURRENT_COMP, x, y, z, vx, vy, vz, t, sx, sy, sz, p);
  }
}
#line 17464 "./SE_example2.c"
/* 'analyzer=Pol_mirror()' component instance extend code */
    SIG_MESSAGE("analyzer (Trace:Extend)");
#line 225 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/SE_example2/SE_example2.instr"
  if(!SCATTERED) ABSORB;
#line 17469 "./SE_example2.c"
}   /* End of analyzer=Pol_mirror() SETTING parameter declarations. */
#undef rDownParPtr
#undef rUpParPtr
#undef yheight
#undef zwidth
#undef useTables
#undef rDownPar
#undef rUpPar
#undef rDownFunc
#undef rUpFunc
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompanalyzer:
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

  /* TRACE Component detector_slit [29] */
  mccoordschange(mcposrdetector_slit, mcrotrdetector_slit,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component detector_slit (without coords transformations) */
  mcJumpTrace_detector_slit:
  SIG_MESSAGE("detector_slit (Trace)");
  mcDEBUG_COMP("detector_slit")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbCompdetector_slit
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
#define mccompcurname  detector_slit
#define mccompcurtype  Slit
#define mccompcurindex 29
{   /* Declarations of detector_slit=Slit() SETTING parameters. */
MCNUM xmin = mccdetector_slit_xmin;
MCNUM xmax = mccdetector_slit_xmax;
MCNUM ymin = mccdetector_slit_ymin;
MCNUM ymax = mccdetector_slit_ymax;
MCNUM radius = mccdetector_slit_radius;
MCNUM xwidth = mccdetector_slit_xwidth;
MCNUM yheight = mccdetector_slit_yheight;
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
#line 17602 "./SE_example2.c"
}   /* End of detector_slit=Slit() SETTING parameter declarations. */
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompdetector_slit:
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

  /* TRACE Component pl_mon5 [30] */
  mccoordschange(mcposrpl_mon5, mcrotrpl_mon5,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component pl_mon5 (without coords transformations) */
  mcJumpTrace_pl_mon5:
  SIG_MESSAGE("pl_mon5 (Trace)");
  mcDEBUG_COMP("pl_mon5")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbComppl_mon5
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
#define mccompcurname  pl_mon5
#define mccompcurtype  MeanPolLambda_monitor
#define mccompcurindex 30
#define xwidth mccpl_mon5_xwidth
#define yheight mccpl_mon5_yheight
#define nL mccpl_mon5_nL
#define restore_neutron mccpl_mon5_restore_neutron
#define nowritefile mccpl_mon5_nowritefile
#define PolL_N mccpl_mon5_PolL_N
#define PolL_p mccpl_mon5_PolL_p
#define PolL_p2 mccpl_mon5_PolL_p2
#define HelpArray mccpl_mon5_HelpArray
{   /* Declarations of pl_mon5=MeanPolLambda_monitor() SETTING parameters. */
char* filename = mccpl_mon5_filename;
MCNUM mx = mccpl_mon5_mx;
MCNUM my = mccpl_mon5_my;
MCNUM mz = mccpl_mon5_mz;
MCNUM Lmin = mccpl_mon5_Lmin;
MCNUM Lmax = mccpl_mon5_Lmax;
#line 106 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/MeanPolLambda_monitor.comp"
{
  int i;
  double pol_proj;
  double lambda;

  PROP_Z0;
  lambda = (2*PI/V2K)/sqrt(vx*vx + vy*vy + vz*vz);

  if (inside_rectangle(x, y, xwidth, yheight) &&
      lambda > Lmin && lambda < Lmax) {

    pol_proj = scalar_prod(mx, my, mz, sx, sy, sz);

    i= floor((lambda - Lmin)*nL/(Lmax - Lmin));

    PolL_N[i]++;
    HelpArray[i] += p;
    PolL_p[i]    += pol_proj*p;
    PolL_p2[i]   += pol_proj*pol_proj*p;

    SCATTER;
  }
  if (restore_neutron) {
    RESTORE_NEUTRON(INDEX_CURRENT_COMP, x, y, z, vx, vy, vz, t, sx, sy, sz, p);
  }
}
#line 17750 "./SE_example2.c"
}   /* End of pl_mon5=MeanPolLambda_monitor() SETTING parameter declarations. */
#undef HelpArray
#undef PolL_p2
#undef PolL_p
#undef PolL_N
#undef nowritefile
#undef restore_neutron
#undef nL
#undef yheight
#undef xwidth
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbComppl_mon5:
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

  /* TRACE Component detector [31] */
  mccoordschange(mcposrdetector, mcrotrdetector,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component detector (without coords transformations) */
  mcJumpTrace_detector:
  SIG_MESSAGE("detector (Trace)");
  mcDEBUG_COMP("detector")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbCompdetector
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
#define mccompcurname  detector
#define mccompcurtype  PSD_monitor
#define mccompcurindex 31
#define PSD_N mccdetector_PSD_N
#define PSD_p mccdetector_PSD_p
#define PSD_p2 mccdetector_PSD_p2
{   /* Declarations of detector=PSD_monitor() SETTING parameters. */
int nx = mccdetector_nx;
int ny = mccdetector_ny;
char* filename = mccdetector_filename;
MCNUM xmin = mccdetector_xmin;
MCNUM xmax = mccdetector_xmax;
MCNUM ymin = mccdetector_ymin;
MCNUM ymax = mccdetector_ymax;
MCNUM xwidth = mccdetector_xwidth;
MCNUM yheight = mccdetector_yheight;
MCNUM restore_neutron = mccdetector_restore_neutron;
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
#line 17893 "./SE_example2.c"
}   /* End of detector=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompdetector:
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
#line 18005 "./SE_example2.c"
}   /* End of Origin=Progress_bar() SETTING parameter declarations. */
#undef CurrentTime
#undef EndTime
#undef StartTime
#undef IntermediateCnts
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'lam_mon'. */
  SIG_MESSAGE("lam_mon (Save)");
#define mccompcurname  lam_mon
#define mccompcurtype  L_monitor
#define mccompcurindex 3
#define nL mcclam_mon_nL
#define L_N mcclam_mon_L_N
#define L_p mcclam_mon_L_p
#define L_p2 mcclam_mon_L_p2
{   /* Declarations of lam_mon=L_monitor() SETTING parameters. */
char* filename = mcclam_mon_filename;
MCNUM xmin = mcclam_mon_xmin;
MCNUM xmax = mcclam_mon_xmax;
MCNUM ymin = mcclam_mon_ymin;
MCNUM ymax = mcclam_mon_ymax;
MCNUM xwidth = mcclam_mon_xwidth;
MCNUM yheight = mcclam_mon_yheight;
MCNUM Lmin = mcclam_mon_Lmin;
MCNUM Lmax = mcclam_mon_Lmax;
MCNUM restore_neutron = mcclam_mon_restore_neutron;
int nowritefile = mcclam_mon_nowritefile;
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
#line 18048 "./SE_example2.c"
}   /* End of lam_mon=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'hdiv_mon'. */
  SIG_MESSAGE("hdiv_mon (Save)");
#define mccompcurname  hdiv_mon
#define mccompcurtype  Hdiv_monitor
#define mccompcurindex 5
#define nh mcchdiv_mon_nh
#define Div_N mcchdiv_mon_Div_N
#define Div_p mcchdiv_mon_Div_p
#define Div_p2 mcchdiv_mon_Div_p2
{   /* Declarations of hdiv_mon=Hdiv_monitor() SETTING parameters. */
char* filename = mcchdiv_mon_filename;
MCNUM xmin = mcchdiv_mon_xmin;
MCNUM xmax = mcchdiv_mon_xmax;
MCNUM ymin = mcchdiv_mon_ymin;
MCNUM ymax = mcchdiv_mon_ymax;
MCNUM xwidth = mcchdiv_mon_xwidth;
MCNUM yheight = mcchdiv_mon_yheight;
MCNUM h_maxdiv = mcchdiv_mon_h_maxdiv;
MCNUM restore_neutron = mcchdiv_mon_restore_neutron;
int nowritefile = mcchdiv_mon_nowritefile;
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
#line 18090 "./SE_example2.c"
}   /* End of hdiv_mon=Hdiv_monitor() SETTING parameter declarations. */
#undef Div_p2
#undef Div_p
#undef Div_N
#undef nh
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'psd_ap'. */
  SIG_MESSAGE("psd_ap (Save)");
#define mccompcurname  psd_ap
#define mccompcurtype  PSD_monitor
#define mccompcurindex 7
#define PSD_N mccpsd_ap_PSD_N
#define PSD_p mccpsd_ap_PSD_p
#define PSD_p2 mccpsd_ap_PSD_p2
{   /* Declarations of psd_ap=PSD_monitor() SETTING parameters. */
int nx = mccpsd_ap_nx;
int ny = mccpsd_ap_ny;
char* filename = mccpsd_ap_filename;
MCNUM xmin = mccpsd_ap_xmin;
MCNUM xmax = mccpsd_ap_xmax;
MCNUM ymin = mccpsd_ap_ymin;
MCNUM ymax = mccpsd_ap_ymax;
MCNUM xwidth = mccpsd_ap_xwidth;
MCNUM yheight = mccpsd_ap_yheight;
MCNUM restore_neutron = mccpsd_ap_restore_neutron;
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
#line 18130 "./SE_example2.c"
}   /* End of psd_ap=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'pl_mon0'. */
  SIG_MESSAGE("pl_mon0 (Save)");
#define mccompcurname  pl_mon0
#define mccompcurtype  MeanPolLambda_monitor
#define mccompcurindex 8
#define xwidth mccpl_mon0_xwidth
#define yheight mccpl_mon0_yheight
#define nL mccpl_mon0_nL
#define restore_neutron mccpl_mon0_restore_neutron
#define nowritefile mccpl_mon0_nowritefile
#define PolL_N mccpl_mon0_PolL_N
#define PolL_p mccpl_mon0_PolL_p
#define PolL_p2 mccpl_mon0_PolL_p2
#define HelpArray mccpl_mon0_HelpArray
{   /* Declarations of pl_mon0=MeanPolLambda_monitor() SETTING parameters. */
char* filename = mccpl_mon0_filename;
MCNUM mx = mccpl_mon0_mx;
MCNUM my = mccpl_mon0_my;
MCNUM mz = mccpl_mon0_mz;
MCNUM Lmin = mccpl_mon0_Lmin;
MCNUM Lmax = mccpl_mon0_Lmax;
#line 134 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/MeanPolLambda_monitor.comp"
{
    if (!nowritefile) {
  int i;

  // Average output
  for (i=0; i<nL; i++) {

    if(HelpArray[i] == 0)
      continue;

    PolL_p[i]  /= HelpArray[i];
    // Mcstas uses the error sigma**2=sum p_i**2 for intensities
    // But here we have the mean so the error**2 should be VAR/N:
    // sigma**2 = [sum s_i*s_i*p_i/sum p_i - (sum s_i*p_i/sum p_i)**2]/N
    PolL_p2[i] /= HelpArray[i];
    PolL_p2[i] -= PolL_p[i]*PolL_p[i];
    PolL_p2[i] /= PolL_N[i];
  }

  DETECTOR_OUT_1D("Pol-wavelength monitor",
		  "Wavelength [AA]",
		  "Mean Polarisation",
		  "Wavelength", Lmin, Lmax, nL,
		  &PolL_N[0],&PolL_p[0],&PolL_p2[0],
		  filename);
    }
}
#line 18188 "./SE_example2.c"
}   /* End of pl_mon0=MeanPolLambda_monitor() SETTING parameter declarations. */
#undef HelpArray
#undef PolL_p2
#undef PolL_p
#undef PolL_N
#undef nowritefile
#undef restore_neutron
#undef nL
#undef yheight
#undef xwidth
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'pl_mon1'. */
  SIG_MESSAGE("pl_mon1 (Save)");
#define mccompcurname  pl_mon1
#define mccompcurtype  MeanPolLambda_monitor
#define mccompcurindex 10
#define xwidth mccpl_mon1_xwidth
#define yheight mccpl_mon1_yheight
#define nL mccpl_mon1_nL
#define restore_neutron mccpl_mon1_restore_neutron
#define nowritefile mccpl_mon1_nowritefile
#define PolL_N mccpl_mon1_PolL_N
#define PolL_p mccpl_mon1_PolL_p
#define PolL_p2 mccpl_mon1_PolL_p2
#define HelpArray mccpl_mon1_HelpArray
{   /* Declarations of pl_mon1=MeanPolLambda_monitor() SETTING parameters. */
char* filename = mccpl_mon1_filename;
MCNUM mx = mccpl_mon1_mx;
MCNUM my = mccpl_mon1_my;
MCNUM mz = mccpl_mon1_mz;
MCNUM Lmin = mccpl_mon1_Lmin;
MCNUM Lmax = mccpl_mon1_Lmax;
#line 134 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/MeanPolLambda_monitor.comp"
{
    if (!nowritefile) {
  int i;

  // Average output
  for (i=0; i<nL; i++) {

    if(HelpArray[i] == 0)
      continue;

    PolL_p[i]  /= HelpArray[i];
    // Mcstas uses the error sigma**2=sum p_i**2 for intensities
    // But here we have the mean so the error**2 should be VAR/N:
    // sigma**2 = [sum s_i*s_i*p_i/sum p_i - (sum s_i*p_i/sum p_i)**2]/N
    PolL_p2[i] /= HelpArray[i];
    PolL_p2[i] -= PolL_p[i]*PolL_p[i];
    PolL_p2[i] /= PolL_N[i];
  }

  DETECTOR_OUT_1D("Pol-wavelength monitor",
		  "Wavelength [AA]",
		  "Mean Polarisation",
		  "Wavelength", Lmin, Lmax, nL,
		  &PolL_N[0],&PolL_p[0],&PolL_p2[0],
		  filename);
    }
}
#line 18252 "./SE_example2.c"
}   /* End of pl_mon1=MeanPolLambda_monitor() SETTING parameter declarations. */
#undef HelpArray
#undef PolL_p2
#undef PolL_p
#undef PolL_N
#undef nowritefile
#undef restore_neutron
#undef nL
#undef yheight
#undef xwidth
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'pl_mon2'. */
  SIG_MESSAGE("pl_mon2 (Save)");
#define mccompcurname  pl_mon2
#define mccompcurtype  MeanPolLambda_monitor
#define mccompcurindex 13
#define xwidth mccpl_mon2_xwidth
#define yheight mccpl_mon2_yheight
#define nL mccpl_mon2_nL
#define restore_neutron mccpl_mon2_restore_neutron
#define nowritefile mccpl_mon2_nowritefile
#define PolL_N mccpl_mon2_PolL_N
#define PolL_p mccpl_mon2_PolL_p
#define PolL_p2 mccpl_mon2_PolL_p2
#define HelpArray mccpl_mon2_HelpArray
{   /* Declarations of pl_mon2=MeanPolLambda_monitor() SETTING parameters. */
char* filename = mccpl_mon2_filename;
MCNUM mx = mccpl_mon2_mx;
MCNUM my = mccpl_mon2_my;
MCNUM mz = mccpl_mon2_mz;
MCNUM Lmin = mccpl_mon2_Lmin;
MCNUM Lmax = mccpl_mon2_Lmax;
#line 134 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/MeanPolLambda_monitor.comp"
{
    if (!nowritefile) {
  int i;

  // Average output
  for (i=0; i<nL; i++) {

    if(HelpArray[i] == 0)
      continue;

    PolL_p[i]  /= HelpArray[i];
    // Mcstas uses the error sigma**2=sum p_i**2 for intensities
    // But here we have the mean so the error**2 should be VAR/N:
    // sigma**2 = [sum s_i*s_i*p_i/sum p_i - (sum s_i*p_i/sum p_i)**2]/N
    PolL_p2[i] /= HelpArray[i];
    PolL_p2[i] -= PolL_p[i]*PolL_p[i];
    PolL_p2[i] /= PolL_N[i];
  }

  DETECTOR_OUT_1D("Pol-wavelength monitor",
		  "Wavelength [AA]",
		  "Mean Polarisation",
		  "Wavelength", Lmin, Lmax, nL,
		  &PolL_N[0],&PolL_p[0],&PolL_p2[0],
		  filename);
    }
}
#line 18316 "./SE_example2.c"
}   /* End of pl_mon2=MeanPolLambda_monitor() SETTING parameter declarations. */
#undef HelpArray
#undef PolL_p2
#undef PolL_p
#undef PolL_N
#undef nowritefile
#undef restore_neutron
#undef nL
#undef yheight
#undef xwidth
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'pl_mon3'. */
  SIG_MESSAGE("pl_mon3 (Save)");
#define mccompcurname  pl_mon3
#define mccompcurtype  MeanPolLambda_monitor
#define mccompcurindex 18
#define xwidth mccpl_mon3_xwidth
#define yheight mccpl_mon3_yheight
#define nL mccpl_mon3_nL
#define restore_neutron mccpl_mon3_restore_neutron
#define nowritefile mccpl_mon3_nowritefile
#define PolL_N mccpl_mon3_PolL_N
#define PolL_p mccpl_mon3_PolL_p
#define PolL_p2 mccpl_mon3_PolL_p2
#define HelpArray mccpl_mon3_HelpArray
{   /* Declarations of pl_mon3=MeanPolLambda_monitor() SETTING parameters. */
char* filename = mccpl_mon3_filename;
MCNUM mx = mccpl_mon3_mx;
MCNUM my = mccpl_mon3_my;
MCNUM mz = mccpl_mon3_mz;
MCNUM Lmin = mccpl_mon3_Lmin;
MCNUM Lmax = mccpl_mon3_Lmax;
#line 134 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/MeanPolLambda_monitor.comp"
{
    if (!nowritefile) {
  int i;

  // Average output
  for (i=0; i<nL; i++) {

    if(HelpArray[i] == 0)
      continue;

    PolL_p[i]  /= HelpArray[i];
    // Mcstas uses the error sigma**2=sum p_i**2 for intensities
    // But here we have the mean so the error**2 should be VAR/N:
    // sigma**2 = [sum s_i*s_i*p_i/sum p_i - (sum s_i*p_i/sum p_i)**2]/N
    PolL_p2[i] /= HelpArray[i];
    PolL_p2[i] -= PolL_p[i]*PolL_p[i];
    PolL_p2[i] /= PolL_N[i];
  }

  DETECTOR_OUT_1D("Pol-wavelength monitor",
		  "Wavelength [AA]",
		  "Mean Polarisation",
		  "Wavelength", Lmin, Lmax, nL,
		  &PolL_N[0],&PolL_p[0],&PolL_p2[0],
		  filename);
    }
}
#line 18380 "./SE_example2.c"
}   /* End of pl_mon3=MeanPolLambda_monitor() SETTING parameter declarations. */
#undef HelpArray
#undef PolL_p2
#undef PolL_p
#undef PolL_N
#undef nowritefile
#undef restore_neutron
#undef nL
#undef yheight
#undef xwidth
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'pl_mon40'. */
  SIG_MESSAGE("pl_mon40 (Save)");
#define mccompcurname  pl_mon40
#define mccompcurtype  MeanPolLambda_monitor
#define mccompcurindex 22
#define xwidth mccpl_mon40_xwidth
#define yheight mccpl_mon40_yheight
#define nL mccpl_mon40_nL
#define restore_neutron mccpl_mon40_restore_neutron
#define nowritefile mccpl_mon40_nowritefile
#define PolL_N mccpl_mon40_PolL_N
#define PolL_p mccpl_mon40_PolL_p
#define PolL_p2 mccpl_mon40_PolL_p2
#define HelpArray mccpl_mon40_HelpArray
{   /* Declarations of pl_mon40=MeanPolLambda_monitor() SETTING parameters. */
char* filename = mccpl_mon40_filename;
MCNUM mx = mccpl_mon40_mx;
MCNUM my = mccpl_mon40_my;
MCNUM mz = mccpl_mon40_mz;
MCNUM Lmin = mccpl_mon40_Lmin;
MCNUM Lmax = mccpl_mon40_Lmax;
#line 134 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/MeanPolLambda_monitor.comp"
{
    if (!nowritefile) {
  int i;

  // Average output
  for (i=0; i<nL; i++) {

    if(HelpArray[i] == 0)
      continue;

    PolL_p[i]  /= HelpArray[i];
    // Mcstas uses the error sigma**2=sum p_i**2 for intensities
    // But here we have the mean so the error**2 should be VAR/N:
    // sigma**2 = [sum s_i*s_i*p_i/sum p_i - (sum s_i*p_i/sum p_i)**2]/N
    PolL_p2[i] /= HelpArray[i];
    PolL_p2[i] -= PolL_p[i]*PolL_p[i];
    PolL_p2[i] /= PolL_N[i];
  }

  DETECTOR_OUT_1D("Pol-wavelength monitor",
		  "Wavelength [AA]",
		  "Mean Polarisation",
		  "Wavelength", Lmin, Lmax, nL,
		  &PolL_N[0],&PolL_p[0],&PolL_p2[0],
		  filename);
    }
}
#line 18444 "./SE_example2.c"
}   /* End of pl_mon40=MeanPolLambda_monitor() SETTING parameter declarations. */
#undef HelpArray
#undef PolL_p2
#undef PolL_p
#undef PolL_N
#undef nowritefile
#undef restore_neutron
#undef nL
#undef yheight
#undef xwidth
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'pl_mon4'. */
  SIG_MESSAGE("pl_mon4 (Save)");
#define mccompcurname  pl_mon4
#define mccompcurtype  MeanPolLambda_monitor
#define mccompcurindex 24
#define xwidth mccpl_mon4_xwidth
#define yheight mccpl_mon4_yheight
#define nL mccpl_mon4_nL
#define restore_neutron mccpl_mon4_restore_neutron
#define nowritefile mccpl_mon4_nowritefile
#define PolL_N mccpl_mon4_PolL_N
#define PolL_p mccpl_mon4_PolL_p
#define PolL_p2 mccpl_mon4_PolL_p2
#define HelpArray mccpl_mon4_HelpArray
{   /* Declarations of pl_mon4=MeanPolLambda_monitor() SETTING parameters. */
char* filename = mccpl_mon4_filename;
MCNUM mx = mccpl_mon4_mx;
MCNUM my = mccpl_mon4_my;
MCNUM mz = mccpl_mon4_mz;
MCNUM Lmin = mccpl_mon4_Lmin;
MCNUM Lmax = mccpl_mon4_Lmax;
#line 134 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/MeanPolLambda_monitor.comp"
{
    if (!nowritefile) {
  int i;

  // Average output
  for (i=0; i<nL; i++) {

    if(HelpArray[i] == 0)
      continue;

    PolL_p[i]  /= HelpArray[i];
    // Mcstas uses the error sigma**2=sum p_i**2 for intensities
    // But here we have the mean so the error**2 should be VAR/N:
    // sigma**2 = [sum s_i*s_i*p_i/sum p_i - (sum s_i*p_i/sum p_i)**2]/N
    PolL_p2[i] /= HelpArray[i];
    PolL_p2[i] -= PolL_p[i]*PolL_p[i];
    PolL_p2[i] /= PolL_N[i];
  }

  DETECTOR_OUT_1D("Pol-wavelength monitor",
		  "Wavelength [AA]",
		  "Mean Polarisation",
		  "Wavelength", Lmin, Lmax, nL,
		  &PolL_N[0],&PolL_p[0],&PolL_p2[0],
		  filename);
    }
}
#line 18508 "./SE_example2.c"
}   /* End of pl_mon4=MeanPolLambda_monitor() SETTING parameter declarations. */
#undef HelpArray
#undef PolL_p2
#undef PolL_p
#undef PolL_N
#undef nowritefile
#undef restore_neutron
#undef nL
#undef yheight
#undef xwidth
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'psd_ba'. */
  SIG_MESSAGE("psd_ba (Save)");
#define mccompcurname  psd_ba
#define mccompcurtype  PSD_monitor
#define mccompcurindex 26
#define PSD_N mccpsd_ba_PSD_N
#define PSD_p mccpsd_ba_PSD_p
#define PSD_p2 mccpsd_ba_PSD_p2
{   /* Declarations of psd_ba=PSD_monitor() SETTING parameters. */
int nx = mccpsd_ba_nx;
int ny = mccpsd_ba_ny;
char* filename = mccpsd_ba_filename;
MCNUM xmin = mccpsd_ba_xmin;
MCNUM xmax = mccpsd_ba_xmax;
MCNUM ymin = mccpsd_ba_ymin;
MCNUM ymax = mccpsd_ba_ymax;
MCNUM xwidth = mccpsd_ba_xwidth;
MCNUM yheight = mccpsd_ba_yheight;
MCNUM restore_neutron = mccpsd_ba_restore_neutron;
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
#line 18553 "./SE_example2.c"
}   /* End of psd_ba=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'pl_mon_ba'. */
  SIG_MESSAGE("pl_mon_ba (Save)");
#define mccompcurname  pl_mon_ba
#define mccompcurtype  MeanPolLambda_monitor
#define mccompcurindex 27
#define xwidth mccpl_mon_ba_xwidth
#define yheight mccpl_mon_ba_yheight
#define nL mccpl_mon_ba_nL
#define restore_neutron mccpl_mon_ba_restore_neutron
#define nowritefile mccpl_mon_ba_nowritefile
#define PolL_N mccpl_mon_ba_PolL_N
#define PolL_p mccpl_mon_ba_PolL_p
#define PolL_p2 mccpl_mon_ba_PolL_p2
#define HelpArray mccpl_mon_ba_HelpArray
{   /* Declarations of pl_mon_ba=MeanPolLambda_monitor() SETTING parameters. */
char* filename = mccpl_mon_ba_filename;
MCNUM mx = mccpl_mon_ba_mx;
MCNUM my = mccpl_mon_ba_my;
MCNUM mz = mccpl_mon_ba_mz;
MCNUM Lmin = mccpl_mon_ba_Lmin;
MCNUM Lmax = mccpl_mon_ba_Lmax;
#line 134 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/MeanPolLambda_monitor.comp"
{
    if (!nowritefile) {
  int i;

  // Average output
  for (i=0; i<nL; i++) {

    if(HelpArray[i] == 0)
      continue;

    PolL_p[i]  /= HelpArray[i];
    // Mcstas uses the error sigma**2=sum p_i**2 for intensities
    // But here we have the mean so the error**2 should be VAR/N:
    // sigma**2 = [sum s_i*s_i*p_i/sum p_i - (sum s_i*p_i/sum p_i)**2]/N
    PolL_p2[i] /= HelpArray[i];
    PolL_p2[i] -= PolL_p[i]*PolL_p[i];
    PolL_p2[i] /= PolL_N[i];
  }

  DETECTOR_OUT_1D("Pol-wavelength monitor",
		  "Wavelength [AA]",
		  "Mean Polarisation",
		  "Wavelength", Lmin, Lmax, nL,
		  &PolL_N[0],&PolL_p[0],&PolL_p2[0],
		  filename);
    }
}
#line 18611 "./SE_example2.c"
}   /* End of pl_mon_ba=MeanPolLambda_monitor() SETTING parameter declarations. */
#undef HelpArray
#undef PolL_p2
#undef PolL_p
#undef PolL_N
#undef nowritefile
#undef restore_neutron
#undef nL
#undef yheight
#undef xwidth
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'pl_mon5'. */
  SIG_MESSAGE("pl_mon5 (Save)");
#define mccompcurname  pl_mon5
#define mccompcurtype  MeanPolLambda_monitor
#define mccompcurindex 30
#define xwidth mccpl_mon5_xwidth
#define yheight mccpl_mon5_yheight
#define nL mccpl_mon5_nL
#define restore_neutron mccpl_mon5_restore_neutron
#define nowritefile mccpl_mon5_nowritefile
#define PolL_N mccpl_mon5_PolL_N
#define PolL_p mccpl_mon5_PolL_p
#define PolL_p2 mccpl_mon5_PolL_p2
#define HelpArray mccpl_mon5_HelpArray
{   /* Declarations of pl_mon5=MeanPolLambda_monitor() SETTING parameters. */
char* filename = mccpl_mon5_filename;
MCNUM mx = mccpl_mon5_mx;
MCNUM my = mccpl_mon5_my;
MCNUM mz = mccpl_mon5_mz;
MCNUM Lmin = mccpl_mon5_Lmin;
MCNUM Lmax = mccpl_mon5_Lmax;
#line 134 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/MeanPolLambda_monitor.comp"
{
    if (!nowritefile) {
  int i;

  // Average output
  for (i=0; i<nL; i++) {

    if(HelpArray[i] == 0)
      continue;

    PolL_p[i]  /= HelpArray[i];
    // Mcstas uses the error sigma**2=sum p_i**2 for intensities
    // But here we have the mean so the error**2 should be VAR/N:
    // sigma**2 = [sum s_i*s_i*p_i/sum p_i - (sum s_i*p_i/sum p_i)**2]/N
    PolL_p2[i] /= HelpArray[i];
    PolL_p2[i] -= PolL_p[i]*PolL_p[i];
    PolL_p2[i] /= PolL_N[i];
  }

  DETECTOR_OUT_1D("Pol-wavelength monitor",
		  "Wavelength [AA]",
		  "Mean Polarisation",
		  "Wavelength", Lmin, Lmax, nL,
		  &PolL_N[0],&PolL_p[0],&PolL_p2[0],
		  filename);
    }
}
#line 18675 "./SE_example2.c"
}   /* End of pl_mon5=MeanPolLambda_monitor() SETTING parameter declarations. */
#undef HelpArray
#undef PolL_p2
#undef PolL_p
#undef PolL_N
#undef nowritefile
#undef restore_neutron
#undef nL
#undef yheight
#undef xwidth
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'detector'. */
  SIG_MESSAGE("detector (Save)");
#define mccompcurname  detector
#define mccompcurtype  PSD_monitor
#define mccompcurindex 31
#define PSD_N mccdetector_PSD_N
#define PSD_p mccdetector_PSD_p
#define PSD_p2 mccdetector_PSD_p2
{   /* Declarations of detector=PSD_monitor() SETTING parameters. */
int nx = mccdetector_nx;
int ny = mccdetector_ny;
char* filename = mccdetector_filename;
MCNUM xmin = mccdetector_xmin;
MCNUM xmax = mccdetector_xmax;
MCNUM ymin = mccdetector_ymin;
MCNUM ymax = mccdetector_ymax;
MCNUM xwidth = mccdetector_xwidth;
MCNUM yheight = mccdetector_yheight;
MCNUM restore_neutron = mccdetector_restore_neutron;
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
#line 18720 "./SE_example2.c"
}   /* End of detector=PSD_monitor() SETTING parameter declarations. */
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
#line 18763 "./SE_example2.c"
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
    if (!mcNCounter[3]) fprintf(stderr, "Warning: No neutron could reach Component[3] lam_mon\n");
    if (mcAbsorbProp[3]) fprintf(stderr, "Warning: %g events were removed in Component[3] lam_mon=L_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[3]);
    if (!mcNCounter[4]) fprintf(stderr, "Warning: No neutron could reach Component[4] polariser_pt\n");
    if (mcAbsorbProp[4]) fprintf(stderr, "Warning: %g events were removed in Component[4] polariser_pt=Arm()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[4]);
    if (!mcNCounter[5]) fprintf(stderr, "Warning: No neutron could reach Component[5] hdiv_mon\n");
    if (mcAbsorbProp[5]) fprintf(stderr, "Warning: %g events were removed in Component[5] hdiv_mon=Hdiv_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[5]);
    if (!mcNCounter[6]) fprintf(stderr, "Warning: No neutron could reach Component[6] polariser2\n");
    if (mcAbsorbProp[6]) fprintf(stderr, "Warning: %g events were removed in Component[6] polariser2=Pol_mirror()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[6]);
  /* User FINALLY code for component 'psd_ap'. */
  SIG_MESSAGE("psd_ap (Finally)");
#define mccompcurname  psd_ap
#define mccompcurtype  PSD_monitor
#define mccompcurindex 7
#define PSD_N mccpsd_ap_PSD_N
#define PSD_p mccpsd_ap_PSD_p
#define PSD_p2 mccpsd_ap_PSD_p2
{   /* Declarations of psd_ap=PSD_monitor() SETTING parameters. */
int nx = mccpsd_ap_nx;
int ny = mccpsd_ap_ny;
char* filename = mccpsd_ap_filename;
MCNUM xmin = mccpsd_ap_xmin;
MCNUM xmax = mccpsd_ap_xmax;
MCNUM ymin = mccpsd_ap_ymin;
MCNUM ymax = mccpsd_ap_ymax;
MCNUM xwidth = mccpsd_ap_xwidth;
MCNUM yheight = mccpsd_ap_yheight;
MCNUM restore_neutron = mccpsd_ap_restore_neutron;
#line 122 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
{
  destroy_darr2d(PSD_N);
  destroy_darr2d(PSD_p);
  destroy_darr2d(PSD_p2);
}
#line 18810 "./SE_example2.c"
}   /* End of psd_ap=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[7]) fprintf(stderr, "Warning: No neutron could reach Component[7] psd_ap\n");
    if (mcAbsorbProp[7]) fprintf(stderr, "Warning: %g events were removed in Component[7] psd_ap=PSD_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[7]);
    if (!mcNCounter[8]) fprintf(stderr, "Warning: No neutron could reach Component[8] pl_mon0\n");
    if (mcAbsorbProp[8]) fprintf(stderr, "Warning: %g events were removed in Component[8] pl_mon0=MeanPolLambda_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[8]);
    if (!mcNCounter[9]) fprintf(stderr, "Warning: No neutron could reach Component[9] flipper1\n");
    if (mcAbsorbProp[9]) fprintf(stderr, "Warning: %g events were removed in Component[9] flipper1=Pol_FieldBox()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[9]);
    if (!mcNCounter[10]) fprintf(stderr, "Warning: No neutron could reach Component[10] pl_mon1\n");
    if (mcAbsorbProp[10]) fprintf(stderr, "Warning: %g events were removed in Component[10] pl_mon1=MeanPolLambda_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[10]);
    if (!mcNCounter[11]) fprintf(stderr, "Warning: No neutron could reach Component[11] guide_field_1\n");
    if (mcAbsorbProp[11]) fprintf(stderr, "Warning: %g events were removed in Component[11] guide_field_1=Pol_FieldBox()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[11]);
    if (!mcNCounter[12]) fprintf(stderr, "Warning: No neutron could reach Component[12] guide_field_1_end\n");
    if (mcAbsorbProp[12]) fprintf(stderr, "Warning: %g events were removed in Component[12] guide_field_1_end=Slit()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[12]);
    if (!mcNCounter[13]) fprintf(stderr, "Warning: No neutron could reach Component[13] pl_mon2\n");
    if (mcAbsorbProp[13]) fprintf(stderr, "Warning: %g events were removed in Component[13] pl_mon2=MeanPolLambda_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[13]);
    if (!mcNCounter[14]) fprintf(stderr, "Warning: No neutron could reach Component[14] flipper2_1\n");
    if (mcAbsorbProp[14]) fprintf(stderr, "Warning: %g events were removed in Component[14] flipper2_1=Pol_FieldBox()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[14]);
    if (!mcNCounter[15]) fprintf(stderr, "Warning: No neutron could reach Component[15] flipper2_2\n");
    if (mcAbsorbProp[15]) fprintf(stderr, "Warning: %g events were removed in Component[15] flipper2_2=Pol_FieldBox()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[15]);
    if (!mcNCounter[16]) fprintf(stderr, "Warning: No neutron could reach Component[16] flipper_2_2_end\n");
    if (mcAbsorbProp[16]) fprintf(stderr, "Warning: %g events were removed in Component[16] flipper_2_2_end=Slit()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[16]);
    if (!mcNCounter[17]) fprintf(stderr, "Warning: No neutron could reach Component[17] sample_mnt\n");
    if (mcAbsorbProp[17]) fprintf(stderr, "Warning: %g events were removed in Component[17] sample_mnt=Arm()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[17]);
    if (!mcNCounter[18]) fprintf(stderr, "Warning: No neutron could reach Component[18] pl_mon3\n");
    if (mcAbsorbProp[18]) fprintf(stderr, "Warning: %g events were removed in Component[18] pl_mon3=MeanPolLambda_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[18]);
    if (!mcNCounter[19]) fprintf(stderr, "Warning: No neutron could reach Component[19] sample0\n");
    if (mcAbsorbProp[19]) fprintf(stderr, "Warning: %g events were removed in Component[19] sample0=V_sample()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[19]);
    if (!mcNCounter[20]) fprintf(stderr, "Warning: No neutron could reach Component[20] guide_field_2\n");
    if (mcAbsorbProp[20]) fprintf(stderr, "Warning: %g events were removed in Component[20] guide_field_2=Pol_FieldBox()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[20]);
    if (!mcNCounter[21]) fprintf(stderr, "Warning: No neutron could reach Component[21] guide_field_2_end\n");
    if (mcAbsorbProp[21]) fprintf(stderr, "Warning: %g events were removed in Component[21] guide_field_2_end=Slit()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[21]);
    if (!mcNCounter[22]) fprintf(stderr, "Warning: No neutron could reach Component[22] pl_mon40\n");
    if (mcAbsorbProp[22]) fprintf(stderr, "Warning: %g events were removed in Component[22] pl_mon40=MeanPolLambda_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[22]);
    if (!mcNCounter[23]) fprintf(stderr, "Warning: No neutron could reach Component[23] flipper3\n");
    if (mcAbsorbProp[23]) fprintf(stderr, "Warning: %g events were removed in Component[23] flipper3=Pol_FieldBox()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[23]);
    if (!mcNCounter[24]) fprintf(stderr, "Warning: No neutron could reach Component[24] pl_mon4\n");
    if (mcAbsorbProp[24]) fprintf(stderr, "Warning: %g events were removed in Component[24] pl_mon4=MeanPolLambda_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[24]);
    if (!mcNCounter[25]) fprintf(stderr, "Warning: No neutron could reach Component[25] analyzer_arm\n");
    if (mcAbsorbProp[25]) fprintf(stderr, "Warning: %g events were removed in Component[25] analyzer_arm=Arm()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[25]);
  /* User FINALLY code for component 'psd_ba'. */
  SIG_MESSAGE("psd_ba (Finally)");
#define mccompcurname  psd_ba
#define mccompcurtype  PSD_monitor
#define mccompcurindex 26
#define PSD_N mccpsd_ba_PSD_N
#define PSD_p mccpsd_ba_PSD_p
#define PSD_p2 mccpsd_ba_PSD_p2
{   /* Declarations of psd_ba=PSD_monitor() SETTING parameters. */
int nx = mccpsd_ba_nx;
int ny = mccpsd_ba_ny;
char* filename = mccpsd_ba_filename;
MCNUM xmin = mccpsd_ba_xmin;
MCNUM xmax = mccpsd_ba_xmax;
MCNUM ymin = mccpsd_ba_ymin;
MCNUM ymax = mccpsd_ba_ymax;
MCNUM xwidth = mccpsd_ba_xwidth;
MCNUM yheight = mccpsd_ba_yheight;
MCNUM restore_neutron = mccpsd_ba_restore_neutron;
#line 122 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
{
  destroy_darr2d(PSD_N);
  destroy_darr2d(PSD_p);
  destroy_darr2d(PSD_p2);
}
#line 18882 "./SE_example2.c"
}   /* End of psd_ba=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[26]) fprintf(stderr, "Warning: No neutron could reach Component[26] psd_ba\n");
    if (mcAbsorbProp[26]) fprintf(stderr, "Warning: %g events were removed in Component[26] psd_ba=PSD_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[26]);
    if (!mcNCounter[27]) fprintf(stderr, "Warning: No neutron could reach Component[27] pl_mon_ba\n");
    if (mcAbsorbProp[27]) fprintf(stderr, "Warning: %g events were removed in Component[27] pl_mon_ba=MeanPolLambda_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[27]);
    if (!mcNCounter[28]) fprintf(stderr, "Warning: No neutron could reach Component[28] analyzer\n");
    if (mcAbsorbProp[28]) fprintf(stderr, "Warning: %g events were removed in Component[28] analyzer=Pol_mirror()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[28]);
    if (!mcNCounter[29]) fprintf(stderr, "Warning: No neutron could reach Component[29] detector_slit\n");
    if (mcAbsorbProp[29]) fprintf(stderr, "Warning: %g events were removed in Component[29] detector_slit=Slit()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[29]);
    if (!mcNCounter[30]) fprintf(stderr, "Warning: No neutron could reach Component[30] pl_mon5\n");
    if (mcAbsorbProp[30]) fprintf(stderr, "Warning: %g events were removed in Component[30] pl_mon5=MeanPolLambda_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[30]);
  /* User FINALLY code for component 'detector'. */
  SIG_MESSAGE("detector (Finally)");
#define mccompcurname  detector
#define mccompcurtype  PSD_monitor
#define mccompcurindex 31
#define PSD_N mccdetector_PSD_N
#define PSD_p mccdetector_PSD_p
#define PSD_p2 mccdetector_PSD_p2
{   /* Declarations of detector=PSD_monitor() SETTING parameters. */
int nx = mccdetector_nx;
int ny = mccdetector_ny;
char* filename = mccdetector_filename;
MCNUM xmin = mccdetector_xmin;
MCNUM xmax = mccdetector_xmax;
MCNUM ymin = mccdetector_ymin;
MCNUM ymax = mccdetector_ymax;
MCNUM xwidth = mccdetector_xwidth;
MCNUM yheight = mccdetector_yheight;
MCNUM restore_neutron = mccdetector_restore_neutron;
#line 122 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
{
  destroy_darr2d(PSD_N);
  destroy_darr2d(PSD_p);
  destroy_darr2d(PSD_p2);
}
#line 18926 "./SE_example2.c"
}   /* End of detector=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[31]) fprintf(stderr, "Warning: No neutron could reach Component[31] detector\n");
    if (mcAbsorbProp[31]) fprintf(stderr, "Warning: %g events were removed in Component[31] detector=PSD_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[31]);
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
#line 18971 "./SE_example2.c"
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
#line 19020 "./SE_example2.c"
}   /* End of Source=Source_simple() SETTING parameter declarations. */
#undef srcArea
#undef square
#undef pmul
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'lam_mon'. */
  SIG_MESSAGE("lam_mon (McDisplay)");
  printf("MCDISPLAY: component %s\n", "lam_mon");
#define mccompcurname  lam_mon
#define mccompcurtype  L_monitor
#define mccompcurindex 3
#define nL mcclam_mon_nL
#define L_N mcclam_mon_L_N
#define L_p mcclam_mon_L_p
#define L_p2 mcclam_mon_L_p2
{   /* Declarations of lam_mon=L_monitor() SETTING parameters. */
char* filename = mcclam_mon_filename;
MCNUM xmin = mcclam_mon_xmin;
MCNUM xmax = mcclam_mon_xmax;
MCNUM ymin = mcclam_mon_ymin;
MCNUM ymax = mcclam_mon_ymax;
MCNUM xwidth = mcclam_mon_xwidth;
MCNUM yheight = mcclam_mon_yheight;
MCNUM Lmin = mcclam_mon_Lmin;
MCNUM Lmax = mcclam_mon_Lmax;
MCNUM restore_neutron = mcclam_mon_restore_neutron;
int nowritefile = mcclam_mon_nowritefile;
#line 120 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/L_monitor.comp"
{
  
  multiline(5, (double)xmin, (double)ymin, 0.0,
               (double)xmax, (double)ymin, 0.0,
               (double)xmax, (double)ymax, 0.0,
               (double)xmin, (double)ymax, 0.0,
               (double)xmin, (double)ymin, 0.0);
}
#line 19060 "./SE_example2.c"
}   /* End of lam_mon=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'polariser_pt'. */
  SIG_MESSAGE("polariser_pt (McDisplay)");
  printf("MCDISPLAY: component %s\n", "polariser_pt");
#define mccompcurname  polariser_pt
#define mccompcurtype  Arm
#define mccompcurindex 4
#line 40 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Arm.comp"
{
  /* A bit ugly; hard-coded dimensions. */
  
  line(0,0,0,0.2,0,0);
  line(0,0,0,0,0.2,0);
  line(0,0,0,0,0,0.2);
}
#line 19084 "./SE_example2.c"
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'hdiv_mon'. */
  SIG_MESSAGE("hdiv_mon (McDisplay)");
  printf("MCDISPLAY: component %s\n", "hdiv_mon");
#define mccompcurname  hdiv_mon
#define mccompcurtype  Hdiv_monitor
#define mccompcurindex 5
#define nh mcchdiv_mon_nh
#define Div_N mcchdiv_mon_Div_N
#define Div_p mcchdiv_mon_Div_p
#define Div_p2 mcchdiv_mon_Div_p2
{   /* Declarations of hdiv_mon=Hdiv_monitor() SETTING parameters. */
char* filename = mcchdiv_mon_filename;
MCNUM xmin = mcchdiv_mon_xmin;
MCNUM xmax = mcchdiv_mon_xmax;
MCNUM ymin = mcchdiv_mon_ymin;
MCNUM ymax = mcchdiv_mon_ymax;
MCNUM xwidth = mcchdiv_mon_xwidth;
MCNUM yheight = mcchdiv_mon_yheight;
MCNUM h_maxdiv = mcchdiv_mon_h_maxdiv;
MCNUM restore_neutron = mcchdiv_mon_restore_neutron;
int nowritefile = mcchdiv_mon_nowritefile;
#line 120 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/Hdiv_monitor.comp"
{
    
    multiline(5, (double)xmin, (double)ymin, 0.0,
                 (double)xmax, (double)ymin, 0.0,
                 (double)xmax, (double)ymax, 0.0,
                 (double)xmin, (double)ymax, 0.0,
                 (double)xmin, (double)ymin, 0.0);
}
#line 19119 "./SE_example2.c"
}   /* End of hdiv_mon=Hdiv_monitor() SETTING parameter declarations. */
#undef Div_p2
#undef Div_p
#undef Div_N
#undef nh
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'polariser2'. */
  SIG_MESSAGE("polariser2 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "polariser2");
#define mccompcurname  polariser2
#define mccompcurtype  Pol_mirror
#define mccompcurindex 6
#define rUpFunc mccpolariser2_rUpFunc
#define rDownFunc mccpolariser2_rDownFunc
#define rUpPar mccpolariser2_rUpPar
#define rDownPar mccpolariser2_rDownPar
#define useTables mccpolariser2_useTables
#define zwidth mccpolariser2_zwidth
#define yheight mccpolariser2_yheight
#define rUpParPtr mccpolariser2_rUpParPtr
#define rDownParPtr mccpolariser2_rDownParPtr
{   /* Declarations of polariser2=Pol_mirror() SETTING parameters. */
MCNUM p_reflect = mccpolariser2_p_reflect;
#line 208 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Pol_mirror.comp"
{
  
  rectangle("yz", 0, 0, 0, zwidth, yheight);
}
#line 19151 "./SE_example2.c"
}   /* End of polariser2=Pol_mirror() SETTING parameter declarations. */
#undef rDownParPtr
#undef rUpParPtr
#undef yheight
#undef zwidth
#undef useTables
#undef rDownPar
#undef rUpPar
#undef rDownFunc
#undef rUpFunc
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'psd_ap'. */
  SIG_MESSAGE("psd_ap (McDisplay)");
  printf("MCDISPLAY: component %s\n", "psd_ap");
#define mccompcurname  psd_ap
#define mccompcurtype  PSD_monitor
#define mccompcurindex 7
#define PSD_N mccpsd_ap_PSD_N
#define PSD_p mccpsd_ap_PSD_p
#define PSD_p2 mccpsd_ap_PSD_p2
{   /* Declarations of psd_ap=PSD_monitor() SETTING parameters. */
int nx = mccpsd_ap_nx;
int ny = mccpsd_ap_ny;
char* filename = mccpsd_ap_filename;
MCNUM xmin = mccpsd_ap_xmin;
MCNUM xmax = mccpsd_ap_xmax;
MCNUM ymin = mccpsd_ap_ymin;
MCNUM ymax = mccpsd_ap_ymax;
MCNUM xwidth = mccpsd_ap_xwidth;
MCNUM yheight = mccpsd_ap_yheight;
MCNUM restore_neutron = mccpsd_ap_restore_neutron;
#line 129 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
{
  multiline(5,
    (double)xmin, (double)ymin, 0.0,
    (double)xmax, (double)ymin, 0.0,
    (double)xmax, (double)ymax, 0.0,
    (double)xmin, (double)ymax, 0.0,
    (double)xmin, (double)ymin, 0.0);
}
#line 19195 "./SE_example2.c"
}   /* End of psd_ap=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'pl_mon0'. */
  SIG_MESSAGE("pl_mon0 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "pl_mon0");
#define mccompcurname  pl_mon0
#define mccompcurtype  MeanPolLambda_monitor
#define mccompcurindex 8
#define xwidth mccpl_mon0_xwidth
#define yheight mccpl_mon0_yheight
#define nL mccpl_mon0_nL
#define restore_neutron mccpl_mon0_restore_neutron
#define nowritefile mccpl_mon0_nowritefile
#define PolL_N mccpl_mon0_PolL_N
#define PolL_p mccpl_mon0_PolL_p
#define PolL_p2 mccpl_mon0_PolL_p2
#define HelpArray mccpl_mon0_HelpArray
{   /* Declarations of pl_mon0=MeanPolLambda_monitor() SETTING parameters. */
char* filename = mccpl_mon0_filename;
MCNUM mx = mccpl_mon0_mx;
MCNUM my = mccpl_mon0_my;
MCNUM mz = mccpl_mon0_mz;
MCNUM Lmin = mccpl_mon0_Lmin;
MCNUM Lmax = mccpl_mon0_Lmax;
#line 163 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/MeanPolLambda_monitor.comp"
{
  
  rectangle("xy", 0, 0, 0, xwidth, yheight);
}
#line 19231 "./SE_example2.c"
}   /* End of pl_mon0=MeanPolLambda_monitor() SETTING parameter declarations. */
#undef HelpArray
#undef PolL_p2
#undef PolL_p
#undef PolL_N
#undef nowritefile
#undef restore_neutron
#undef nL
#undef yheight
#undef xwidth
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'flipper1'. */
  SIG_MESSAGE("flipper1 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "flipper1");
#define mccompcurname  flipper1
#define mccompcurtype  Pol_FieldBox
#define mccompcurindex 9
#define fieldFunction mccflipper1_fieldFunction
#define field_parameters mccflipper1_field_parameters
#define const_field_parameters mccflipper1_const_field_parameters
{   /* Declarations of flipper1=Pol_FieldBox() SETTING parameters. */
MCNUM xwidth = mccflipper1_xwidth;
MCNUM yheight = mccflipper1_yheight;
MCNUM zdepth = mccflipper1_zdepth;
MCNUM Bx = mccflipper1_Bx;
MCNUM By = mccflipper1_By;
MCNUM Bz = mccflipper1_Bz;
char* filename = mccflipper1_filename;
#line 88 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Pol_FieldBox.comp"
{
  const int nDash = 10;
  double xw_2,yh_2,zd_2;
  xw_2=xwidth/2.0;yh_2=yheight/2.0;zd_2=zdepth/2.0;
  /*entrance*/
  dashed_line(-xw_2, -yh_2, -zd_2,  xw_2, -yh_2, -zd_2, nDash);
  dashed_line(-xw_2, -yh_2, -zd_2, -xw_2,  yh_2, -zd_2, nDash);
  dashed_line( xw_2,  yh_2, -zd_2, -xw_2,  yh_2, -zd_2, nDash);
  dashed_line( xw_2,  yh_2, -zd_2,  xw_2, -yh_2, -zd_2, nDash);

  /*exit*/
  dashed_line(-xw_2, -yh_2, +zd_2,  xw_2, -yh_2, +zd_2, nDash);
  dashed_line(-xw_2, -yh_2, +zd_2, -xw_2,  yh_2, +zd_2, nDash);
  dashed_line( xw_2,  yh_2, +zd_2, -xw_2,  yh_2, +zd_2, nDash);
  dashed_line( xw_2,  yh_2, +zd_2,  xw_2, -yh_2, +zd_2, nDash);

  /*4 lines to make a box*/
  dashed_line(-xw_2, -yh_2, -zd_2, -xw_2, -yh_2, +zd_2, nDash);
  dashed_line(-xw_2,  yh_2, -zd_2, -xw_2,  yh_2, +zd_2, nDash);
  dashed_line( xw_2, -yh_2, -zd_2,  xw_2, -yh_2, +zd_2, nDash);
  dashed_line( xw_2,  yh_2, -zd_2,  xw_2,  yh_2, +zd_2, nDash);
}
#line 19286 "./SE_example2.c"
}   /* End of flipper1=Pol_FieldBox() SETTING parameter declarations. */
#undef const_field_parameters
#undef field_parameters
#undef fieldFunction
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'pl_mon1'. */
  SIG_MESSAGE("pl_mon1 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "pl_mon1");
#define mccompcurname  pl_mon1
#define mccompcurtype  MeanPolLambda_monitor
#define mccompcurindex 10
#define xwidth mccpl_mon1_xwidth
#define yheight mccpl_mon1_yheight
#define nL mccpl_mon1_nL
#define restore_neutron mccpl_mon1_restore_neutron
#define nowritefile mccpl_mon1_nowritefile
#define PolL_N mccpl_mon1_PolL_N
#define PolL_p mccpl_mon1_PolL_p
#define PolL_p2 mccpl_mon1_PolL_p2
#define HelpArray mccpl_mon1_HelpArray
{   /* Declarations of pl_mon1=MeanPolLambda_monitor() SETTING parameters. */
char* filename = mccpl_mon1_filename;
MCNUM mx = mccpl_mon1_mx;
MCNUM my = mccpl_mon1_my;
MCNUM mz = mccpl_mon1_mz;
MCNUM Lmin = mccpl_mon1_Lmin;
MCNUM Lmax = mccpl_mon1_Lmax;
#line 163 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/MeanPolLambda_monitor.comp"
{
  
  rectangle("xy", 0, 0, 0, xwidth, yheight);
}
#line 19322 "./SE_example2.c"
}   /* End of pl_mon1=MeanPolLambda_monitor() SETTING parameter declarations. */
#undef HelpArray
#undef PolL_p2
#undef PolL_p
#undef PolL_N
#undef nowritefile
#undef restore_neutron
#undef nL
#undef yheight
#undef xwidth
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'guide_field_1'. */
  SIG_MESSAGE("guide_field_1 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "guide_field_1");
#define mccompcurname  guide_field_1
#define mccompcurtype  Pol_FieldBox
#define mccompcurindex 11
#define fieldFunction mccguide_field_1_fieldFunction
#define field_parameters mccguide_field_1_field_parameters
#define const_field_parameters mccguide_field_1_const_field_parameters
{   /* Declarations of guide_field_1=Pol_FieldBox() SETTING parameters. */
MCNUM xwidth = mccguide_field_1_xwidth;
MCNUM yheight = mccguide_field_1_yheight;
MCNUM zdepth = mccguide_field_1_zdepth;
MCNUM Bx = mccguide_field_1_Bx;
MCNUM By = mccguide_field_1_By;
MCNUM Bz = mccguide_field_1_Bz;
char* filename = mccguide_field_1_filename;
#line 88 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Pol_FieldBox.comp"
{
  const int nDash = 10;
  double xw_2,yh_2,zd_2;
  xw_2=xwidth/2.0;yh_2=yheight/2.0;zd_2=zdepth/2.0;
  /*entrance*/
  dashed_line(-xw_2, -yh_2, -zd_2,  xw_2, -yh_2, -zd_2, nDash);
  dashed_line(-xw_2, -yh_2, -zd_2, -xw_2,  yh_2, -zd_2, nDash);
  dashed_line( xw_2,  yh_2, -zd_2, -xw_2,  yh_2, -zd_2, nDash);
  dashed_line( xw_2,  yh_2, -zd_2,  xw_2, -yh_2, -zd_2, nDash);

  /*exit*/
  dashed_line(-xw_2, -yh_2, +zd_2,  xw_2, -yh_2, +zd_2, nDash);
  dashed_line(-xw_2, -yh_2, +zd_2, -xw_2,  yh_2, +zd_2, nDash);
  dashed_line( xw_2,  yh_2, +zd_2, -xw_2,  yh_2, +zd_2, nDash);
  dashed_line( xw_2,  yh_2, +zd_2,  xw_2, -yh_2, +zd_2, nDash);

  /*4 lines to make a box*/
  dashed_line(-xw_2, -yh_2, -zd_2, -xw_2, -yh_2, +zd_2, nDash);
  dashed_line(-xw_2,  yh_2, -zd_2, -xw_2,  yh_2, +zd_2, nDash);
  dashed_line( xw_2, -yh_2, -zd_2,  xw_2, -yh_2, +zd_2, nDash);
  dashed_line( xw_2,  yh_2, -zd_2,  xw_2,  yh_2, +zd_2, nDash);
}
#line 19377 "./SE_example2.c"
}   /* End of guide_field_1=Pol_FieldBox() SETTING parameter declarations. */
#undef const_field_parameters
#undef field_parameters
#undef fieldFunction
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'guide_field_1_end'. */
  SIG_MESSAGE("guide_field_1_end (McDisplay)");
  printf("MCDISPLAY: component %s\n", "guide_field_1_end");
#define mccompcurname  guide_field_1_end
#define mccompcurtype  Slit
#define mccompcurindex 12
{   /* Declarations of guide_field_1_end=Slit() SETTING parameters. */
MCNUM xmin = mccguide_field_1_end_xmin;
MCNUM xmax = mccguide_field_1_end_xmax;
MCNUM ymin = mccguide_field_1_end_ymin;
MCNUM ymax = mccguide_field_1_end_ymax;
MCNUM radius = mccguide_field_1_end_radius;
MCNUM xwidth = mccguide_field_1_end_xwidth;
MCNUM yheight = mccguide_field_1_end_yheight;
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
#line 19423 "./SE_example2.c"
}   /* End of guide_field_1_end=Slit() SETTING parameter declarations. */
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'pl_mon2'. */
  SIG_MESSAGE("pl_mon2 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "pl_mon2");
#define mccompcurname  pl_mon2
#define mccompcurtype  MeanPolLambda_monitor
#define mccompcurindex 13
#define xwidth mccpl_mon2_xwidth
#define yheight mccpl_mon2_yheight
#define nL mccpl_mon2_nL
#define restore_neutron mccpl_mon2_restore_neutron
#define nowritefile mccpl_mon2_nowritefile
#define PolL_N mccpl_mon2_PolL_N
#define PolL_p mccpl_mon2_PolL_p
#define PolL_p2 mccpl_mon2_PolL_p2
#define HelpArray mccpl_mon2_HelpArray
{   /* Declarations of pl_mon2=MeanPolLambda_monitor() SETTING parameters. */
char* filename = mccpl_mon2_filename;
MCNUM mx = mccpl_mon2_mx;
MCNUM my = mccpl_mon2_my;
MCNUM mz = mccpl_mon2_mz;
MCNUM Lmin = mccpl_mon2_Lmin;
MCNUM Lmax = mccpl_mon2_Lmax;
#line 163 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/MeanPolLambda_monitor.comp"
{
  
  rectangle("xy", 0, 0, 0, xwidth, yheight);
}
#line 19456 "./SE_example2.c"
}   /* End of pl_mon2=MeanPolLambda_monitor() SETTING parameter declarations. */
#undef HelpArray
#undef PolL_p2
#undef PolL_p
#undef PolL_N
#undef nowritefile
#undef restore_neutron
#undef nL
#undef yheight
#undef xwidth
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'flipper2_1'. */
  SIG_MESSAGE("flipper2_1 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "flipper2_1");
#define mccompcurname  flipper2_1
#define mccompcurtype  Pol_FieldBox
#define mccompcurindex 14
#define fieldFunction mccflipper2_1_fieldFunction
#define field_parameters mccflipper2_1_field_parameters
#define const_field_parameters mccflipper2_1_const_field_parameters
{   /* Declarations of flipper2_1=Pol_FieldBox() SETTING parameters. */
MCNUM xwidth = mccflipper2_1_xwidth;
MCNUM yheight = mccflipper2_1_yheight;
MCNUM zdepth = mccflipper2_1_zdepth;
MCNUM Bx = mccflipper2_1_Bx;
MCNUM By = mccflipper2_1_By;
MCNUM Bz = mccflipper2_1_Bz;
char* filename = mccflipper2_1_filename;
#line 88 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Pol_FieldBox.comp"
{
  const int nDash = 10;
  double xw_2,yh_2,zd_2;
  xw_2=xwidth/2.0;yh_2=yheight/2.0;zd_2=zdepth/2.0;
  /*entrance*/
  dashed_line(-xw_2, -yh_2, -zd_2,  xw_2, -yh_2, -zd_2, nDash);
  dashed_line(-xw_2, -yh_2, -zd_2, -xw_2,  yh_2, -zd_2, nDash);
  dashed_line( xw_2,  yh_2, -zd_2, -xw_2,  yh_2, -zd_2, nDash);
  dashed_line( xw_2,  yh_2, -zd_2,  xw_2, -yh_2, -zd_2, nDash);

  /*exit*/
  dashed_line(-xw_2, -yh_2, +zd_2,  xw_2, -yh_2, +zd_2, nDash);
  dashed_line(-xw_2, -yh_2, +zd_2, -xw_2,  yh_2, +zd_2, nDash);
  dashed_line( xw_2,  yh_2, +zd_2, -xw_2,  yh_2, +zd_2, nDash);
  dashed_line( xw_2,  yh_2, +zd_2,  xw_2, -yh_2, +zd_2, nDash);

  /*4 lines to make a box*/
  dashed_line(-xw_2, -yh_2, -zd_2, -xw_2, -yh_2, +zd_2, nDash);
  dashed_line(-xw_2,  yh_2, -zd_2, -xw_2,  yh_2, +zd_2, nDash);
  dashed_line( xw_2, -yh_2, -zd_2,  xw_2, -yh_2, +zd_2, nDash);
  dashed_line( xw_2,  yh_2, -zd_2,  xw_2,  yh_2, +zd_2, nDash);
}
#line 19511 "./SE_example2.c"
}   /* End of flipper2_1=Pol_FieldBox() SETTING parameter declarations. */
#undef const_field_parameters
#undef field_parameters
#undef fieldFunction
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'flipper2_2'. */
  SIG_MESSAGE("flipper2_2 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "flipper2_2");
#define mccompcurname  flipper2_2
#define mccompcurtype  Pol_FieldBox
#define mccompcurindex 15
#define fieldFunction mccflipper2_2_fieldFunction
#define field_parameters mccflipper2_2_field_parameters
#define const_field_parameters mccflipper2_2_const_field_parameters
{   /* Declarations of flipper2_2=Pol_FieldBox() SETTING parameters. */
MCNUM xwidth = mccflipper2_2_xwidth;
MCNUM yheight = mccflipper2_2_yheight;
MCNUM zdepth = mccflipper2_2_zdepth;
MCNUM Bx = mccflipper2_2_Bx;
MCNUM By = mccflipper2_2_By;
MCNUM Bz = mccflipper2_2_Bz;
char* filename = mccflipper2_2_filename;
#line 88 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Pol_FieldBox.comp"
{
  const int nDash = 10;
  double xw_2,yh_2,zd_2;
  xw_2=xwidth/2.0;yh_2=yheight/2.0;zd_2=zdepth/2.0;
  /*entrance*/
  dashed_line(-xw_2, -yh_2, -zd_2,  xw_2, -yh_2, -zd_2, nDash);
  dashed_line(-xw_2, -yh_2, -zd_2, -xw_2,  yh_2, -zd_2, nDash);
  dashed_line( xw_2,  yh_2, -zd_2, -xw_2,  yh_2, -zd_2, nDash);
  dashed_line( xw_2,  yh_2, -zd_2,  xw_2, -yh_2, -zd_2, nDash);

  /*exit*/
  dashed_line(-xw_2, -yh_2, +zd_2,  xw_2, -yh_2, +zd_2, nDash);
  dashed_line(-xw_2, -yh_2, +zd_2, -xw_2,  yh_2, +zd_2, nDash);
  dashed_line( xw_2,  yh_2, +zd_2, -xw_2,  yh_2, +zd_2, nDash);
  dashed_line( xw_2,  yh_2, +zd_2,  xw_2, -yh_2, +zd_2, nDash);

  /*4 lines to make a box*/
  dashed_line(-xw_2, -yh_2, -zd_2, -xw_2, -yh_2, +zd_2, nDash);
  dashed_line(-xw_2,  yh_2, -zd_2, -xw_2,  yh_2, +zd_2, nDash);
  dashed_line( xw_2, -yh_2, -zd_2,  xw_2, -yh_2, +zd_2, nDash);
  dashed_line( xw_2,  yh_2, -zd_2,  xw_2,  yh_2, +zd_2, nDash);
}
#line 19560 "./SE_example2.c"
}   /* End of flipper2_2=Pol_FieldBox() SETTING parameter declarations. */
#undef const_field_parameters
#undef field_parameters
#undef fieldFunction
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'flipper_2_2_end'. */
  SIG_MESSAGE("flipper_2_2_end (McDisplay)");
  printf("MCDISPLAY: component %s\n", "flipper_2_2_end");
#define mccompcurname  flipper_2_2_end
#define mccompcurtype  Slit
#define mccompcurindex 16
{   /* Declarations of flipper_2_2_end=Slit() SETTING parameters. */
MCNUM xmin = mccflipper_2_2_end_xmin;
MCNUM xmax = mccflipper_2_2_end_xmax;
MCNUM ymin = mccflipper_2_2_end_ymin;
MCNUM ymax = mccflipper_2_2_end_ymax;
MCNUM radius = mccflipper_2_2_end_radius;
MCNUM xwidth = mccflipper_2_2_end_xwidth;
MCNUM yheight = mccflipper_2_2_end_yheight;
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
#line 19606 "./SE_example2.c"
}   /* End of flipper_2_2_end=Slit() SETTING parameter declarations. */
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'sample_mnt'. */
  SIG_MESSAGE("sample_mnt (McDisplay)");
  printf("MCDISPLAY: component %s\n", "sample_mnt");
#define mccompcurname  sample_mnt
#define mccompcurtype  Arm
#define mccompcurindex 17
#line 40 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Arm.comp"
{
  /* A bit ugly; hard-coded dimensions. */
  
  line(0,0,0,0.2,0,0);
  line(0,0,0,0,0.2,0);
  line(0,0,0,0,0,0.2);
}
#line 19626 "./SE_example2.c"
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'pl_mon3'. */
  SIG_MESSAGE("pl_mon3 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "pl_mon3");
#define mccompcurname  pl_mon3
#define mccompcurtype  MeanPolLambda_monitor
#define mccompcurindex 18
#define xwidth mccpl_mon3_xwidth
#define yheight mccpl_mon3_yheight
#define nL mccpl_mon3_nL
#define restore_neutron mccpl_mon3_restore_neutron
#define nowritefile mccpl_mon3_nowritefile
#define PolL_N mccpl_mon3_PolL_N
#define PolL_p mccpl_mon3_PolL_p
#define PolL_p2 mccpl_mon3_PolL_p2
#define HelpArray mccpl_mon3_HelpArray
{   /* Declarations of pl_mon3=MeanPolLambda_monitor() SETTING parameters. */
char* filename = mccpl_mon3_filename;
MCNUM mx = mccpl_mon3_mx;
MCNUM my = mccpl_mon3_my;
MCNUM mz = mccpl_mon3_mz;
MCNUM Lmin = mccpl_mon3_Lmin;
MCNUM Lmax = mccpl_mon3_Lmax;
#line 163 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/MeanPolLambda_monitor.comp"
{
  
  rectangle("xy", 0, 0, 0, xwidth, yheight);
}
#line 19658 "./SE_example2.c"
}   /* End of pl_mon3=MeanPolLambda_monitor() SETTING parameter declarations. */
#undef HelpArray
#undef PolL_p2
#undef PolL_p
#undef PolL_N
#undef nowritefile
#undef restore_neutron
#undef nL
#undef yheight
#undef xwidth
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'sample0'. */
  SIG_MESSAGE("sample0 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "sample0");
#define mccompcurname  sample0
#define mccompcurtype  V_sample
#define mccompcurindex 19
#define VarsV mccsample0_VarsV
{   /* Declarations of sample0=V_sample() SETTING parameters. */
MCNUM radius = mccsample0_radius;
MCNUM thickness = mccsample0_thickness;
MCNUM zdepth = mccsample0_zdepth;
MCNUM Vc = mccsample0_Vc;
MCNUM sigma_abs = mccsample0_sigma_abs;
MCNUM sigma_inc = mccsample0_sigma_inc;
MCNUM radius_i = mccsample0_radius_i;
MCNUM radius_o = mccsample0_radius_o;
MCNUM h = mccsample0_h;
MCNUM focus_r = mccsample0_focus_r;
MCNUM pack = mccsample0_pack;
MCNUM frac = mccsample0_frac;
MCNUM f_QE = mccsample0_f_QE;
MCNUM gamma = mccsample0_gamma;
MCNUM target_x = mccsample0_target_x;
MCNUM target_y = mccsample0_target_y;
MCNUM target_z = mccsample0_target_z;
MCNUM focus_xw = mccsample0_focus_xw;
MCNUM focus_yh = mccsample0_focus_yh;
MCNUM focus_aw = mccsample0_focus_aw;
MCNUM focus_ah = mccsample0_focus_ah;
MCNUM xwidth = mccsample0_xwidth;
MCNUM yheight = mccsample0_yheight;
MCNUM zthick = mccsample0_zthick;
MCNUM rad_sphere = mccsample0_rad_sphere;
MCNUM sig_a = mccsample0_sig_a;
MCNUM sig_i = mccsample0_sig_i;
MCNUM V0 = mccsample0_V0;
int target_index = mccsample0_target_index;
MCNUM multiples = mccsample0_multiples;
#line 320 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../obsolete/V_sample.comp"
{
  
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
}
#line 19758 "./SE_example2.c"
}   /* End of sample0=V_sample() SETTING parameter declarations. */
#undef VarsV
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'guide_field_2'. */
  SIG_MESSAGE("guide_field_2 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "guide_field_2");
#define mccompcurname  guide_field_2
#define mccompcurtype  Pol_FieldBox
#define mccompcurindex 20
#define fieldFunction mccguide_field_2_fieldFunction
#define field_parameters mccguide_field_2_field_parameters
#define const_field_parameters mccguide_field_2_const_field_parameters
{   /* Declarations of guide_field_2=Pol_FieldBox() SETTING parameters. */
MCNUM xwidth = mccguide_field_2_xwidth;
MCNUM yheight = mccguide_field_2_yheight;
MCNUM zdepth = mccguide_field_2_zdepth;
MCNUM Bx = mccguide_field_2_Bx;
MCNUM By = mccguide_field_2_By;
MCNUM Bz = mccguide_field_2_Bz;
char* filename = mccguide_field_2_filename;
#line 88 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Pol_FieldBox.comp"
{
  const int nDash = 10;
  double xw_2,yh_2,zd_2;
  xw_2=xwidth/2.0;yh_2=yheight/2.0;zd_2=zdepth/2.0;
  /*entrance*/
  dashed_line(-xw_2, -yh_2, -zd_2,  xw_2, -yh_2, -zd_2, nDash);
  dashed_line(-xw_2, -yh_2, -zd_2, -xw_2,  yh_2, -zd_2, nDash);
  dashed_line( xw_2,  yh_2, -zd_2, -xw_2,  yh_2, -zd_2, nDash);
  dashed_line( xw_2,  yh_2, -zd_2,  xw_2, -yh_2, -zd_2, nDash);

  /*exit*/
  dashed_line(-xw_2, -yh_2, +zd_2,  xw_2, -yh_2, +zd_2, nDash);
  dashed_line(-xw_2, -yh_2, +zd_2, -xw_2,  yh_2, +zd_2, nDash);
  dashed_line( xw_2,  yh_2, +zd_2, -xw_2,  yh_2, +zd_2, nDash);
  dashed_line( xw_2,  yh_2, +zd_2,  xw_2, -yh_2, +zd_2, nDash);

  /*4 lines to make a box*/
  dashed_line(-xw_2, -yh_2, -zd_2, -xw_2, -yh_2, +zd_2, nDash);
  dashed_line(-xw_2,  yh_2, -zd_2, -xw_2,  yh_2, +zd_2, nDash);
  dashed_line( xw_2, -yh_2, -zd_2,  xw_2, -yh_2, +zd_2, nDash);
  dashed_line( xw_2,  yh_2, -zd_2,  xw_2,  yh_2, +zd_2, nDash);
}
#line 19805 "./SE_example2.c"
}   /* End of guide_field_2=Pol_FieldBox() SETTING parameter declarations. */
#undef const_field_parameters
#undef field_parameters
#undef fieldFunction
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'guide_field_2_end'. */
  SIG_MESSAGE("guide_field_2_end (McDisplay)");
  printf("MCDISPLAY: component %s\n", "guide_field_2_end");
#define mccompcurname  guide_field_2_end
#define mccompcurtype  Slit
#define mccompcurindex 21
{   /* Declarations of guide_field_2_end=Slit() SETTING parameters. */
MCNUM xmin = mccguide_field_2_end_xmin;
MCNUM xmax = mccguide_field_2_end_xmax;
MCNUM ymin = mccguide_field_2_end_ymin;
MCNUM ymax = mccguide_field_2_end_ymax;
MCNUM radius = mccguide_field_2_end_radius;
MCNUM xwidth = mccguide_field_2_end_xwidth;
MCNUM yheight = mccguide_field_2_end_yheight;
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
#line 19851 "./SE_example2.c"
}   /* End of guide_field_2_end=Slit() SETTING parameter declarations. */
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'pl_mon40'. */
  SIG_MESSAGE("pl_mon40 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "pl_mon40");
#define mccompcurname  pl_mon40
#define mccompcurtype  MeanPolLambda_monitor
#define mccompcurindex 22
#define xwidth mccpl_mon40_xwidth
#define yheight mccpl_mon40_yheight
#define nL mccpl_mon40_nL
#define restore_neutron mccpl_mon40_restore_neutron
#define nowritefile mccpl_mon40_nowritefile
#define PolL_N mccpl_mon40_PolL_N
#define PolL_p mccpl_mon40_PolL_p
#define PolL_p2 mccpl_mon40_PolL_p2
#define HelpArray mccpl_mon40_HelpArray
{   /* Declarations of pl_mon40=MeanPolLambda_monitor() SETTING parameters. */
char* filename = mccpl_mon40_filename;
MCNUM mx = mccpl_mon40_mx;
MCNUM my = mccpl_mon40_my;
MCNUM mz = mccpl_mon40_mz;
MCNUM Lmin = mccpl_mon40_Lmin;
MCNUM Lmax = mccpl_mon40_Lmax;
#line 163 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/MeanPolLambda_monitor.comp"
{
  
  rectangle("xy", 0, 0, 0, xwidth, yheight);
}
#line 19884 "./SE_example2.c"
}   /* End of pl_mon40=MeanPolLambda_monitor() SETTING parameter declarations. */
#undef HelpArray
#undef PolL_p2
#undef PolL_p
#undef PolL_N
#undef nowritefile
#undef restore_neutron
#undef nL
#undef yheight
#undef xwidth
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'flipper3'. */
  SIG_MESSAGE("flipper3 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "flipper3");
#define mccompcurname  flipper3
#define mccompcurtype  Pol_FieldBox
#define mccompcurindex 23
#define fieldFunction mccflipper3_fieldFunction
#define field_parameters mccflipper3_field_parameters
#define const_field_parameters mccflipper3_const_field_parameters
{   /* Declarations of flipper3=Pol_FieldBox() SETTING parameters. */
MCNUM xwidth = mccflipper3_xwidth;
MCNUM yheight = mccflipper3_yheight;
MCNUM zdepth = mccflipper3_zdepth;
MCNUM Bx = mccflipper3_Bx;
MCNUM By = mccflipper3_By;
MCNUM Bz = mccflipper3_Bz;
char* filename = mccflipper3_filename;
#line 88 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Pol_FieldBox.comp"
{
  const int nDash = 10;
  double xw_2,yh_2,zd_2;
  xw_2=xwidth/2.0;yh_2=yheight/2.0;zd_2=zdepth/2.0;
  /*entrance*/
  dashed_line(-xw_2, -yh_2, -zd_2,  xw_2, -yh_2, -zd_2, nDash);
  dashed_line(-xw_2, -yh_2, -zd_2, -xw_2,  yh_2, -zd_2, nDash);
  dashed_line( xw_2,  yh_2, -zd_2, -xw_2,  yh_2, -zd_2, nDash);
  dashed_line( xw_2,  yh_2, -zd_2,  xw_2, -yh_2, -zd_2, nDash);

  /*exit*/
  dashed_line(-xw_2, -yh_2, +zd_2,  xw_2, -yh_2, +zd_2, nDash);
  dashed_line(-xw_2, -yh_2, +zd_2, -xw_2,  yh_2, +zd_2, nDash);
  dashed_line( xw_2,  yh_2, +zd_2, -xw_2,  yh_2, +zd_2, nDash);
  dashed_line( xw_2,  yh_2, +zd_2,  xw_2, -yh_2, +zd_2, nDash);

  /*4 lines to make a box*/
  dashed_line(-xw_2, -yh_2, -zd_2, -xw_2, -yh_2, +zd_2, nDash);
  dashed_line(-xw_2,  yh_2, -zd_2, -xw_2,  yh_2, +zd_2, nDash);
  dashed_line( xw_2, -yh_2, -zd_2,  xw_2, -yh_2, +zd_2, nDash);
  dashed_line( xw_2,  yh_2, -zd_2,  xw_2,  yh_2, +zd_2, nDash);
}
#line 19939 "./SE_example2.c"
}   /* End of flipper3=Pol_FieldBox() SETTING parameter declarations. */
#undef const_field_parameters
#undef field_parameters
#undef fieldFunction
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'pl_mon4'. */
  SIG_MESSAGE("pl_mon4 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "pl_mon4");
#define mccompcurname  pl_mon4
#define mccompcurtype  MeanPolLambda_monitor
#define mccompcurindex 24
#define xwidth mccpl_mon4_xwidth
#define yheight mccpl_mon4_yheight
#define nL mccpl_mon4_nL
#define restore_neutron mccpl_mon4_restore_neutron
#define nowritefile mccpl_mon4_nowritefile
#define PolL_N mccpl_mon4_PolL_N
#define PolL_p mccpl_mon4_PolL_p
#define PolL_p2 mccpl_mon4_PolL_p2
#define HelpArray mccpl_mon4_HelpArray
{   /* Declarations of pl_mon4=MeanPolLambda_monitor() SETTING parameters. */
char* filename = mccpl_mon4_filename;
MCNUM mx = mccpl_mon4_mx;
MCNUM my = mccpl_mon4_my;
MCNUM mz = mccpl_mon4_mz;
MCNUM Lmin = mccpl_mon4_Lmin;
MCNUM Lmax = mccpl_mon4_Lmax;
#line 163 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/MeanPolLambda_monitor.comp"
{
  
  rectangle("xy", 0, 0, 0, xwidth, yheight);
}
#line 19975 "./SE_example2.c"
}   /* End of pl_mon4=MeanPolLambda_monitor() SETTING parameter declarations. */
#undef HelpArray
#undef PolL_p2
#undef PolL_p
#undef PolL_N
#undef nowritefile
#undef restore_neutron
#undef nL
#undef yheight
#undef xwidth
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'analyzer_arm'. */
  SIG_MESSAGE("analyzer_arm (McDisplay)");
  printf("MCDISPLAY: component %s\n", "analyzer_arm");
#define mccompcurname  analyzer_arm
#define mccompcurtype  Arm
#define mccompcurindex 25
#line 40 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Arm.comp"
{
  /* A bit ugly; hard-coded dimensions. */
  
  line(0,0,0,0.2,0,0);
  line(0,0,0,0,0.2,0);
  line(0,0,0,0,0,0.2);
}
#line 20004 "./SE_example2.c"
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'psd_ba'. */
  SIG_MESSAGE("psd_ba (McDisplay)");
  printf("MCDISPLAY: component %s\n", "psd_ba");
#define mccompcurname  psd_ba
#define mccompcurtype  PSD_monitor
#define mccompcurindex 26
#define PSD_N mccpsd_ba_PSD_N
#define PSD_p mccpsd_ba_PSD_p
#define PSD_p2 mccpsd_ba_PSD_p2
{   /* Declarations of psd_ba=PSD_monitor() SETTING parameters. */
int nx = mccpsd_ba_nx;
int ny = mccpsd_ba_ny;
char* filename = mccpsd_ba_filename;
MCNUM xmin = mccpsd_ba_xmin;
MCNUM xmax = mccpsd_ba_xmax;
MCNUM ymin = mccpsd_ba_ymin;
MCNUM ymax = mccpsd_ba_ymax;
MCNUM xwidth = mccpsd_ba_xwidth;
MCNUM yheight = mccpsd_ba_yheight;
MCNUM restore_neutron = mccpsd_ba_restore_neutron;
#line 129 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
{
  multiline(5,
    (double)xmin, (double)ymin, 0.0,
    (double)xmax, (double)ymin, 0.0,
    (double)xmax, (double)ymax, 0.0,
    (double)xmin, (double)ymax, 0.0,
    (double)xmin, (double)ymin, 0.0);
}
#line 20038 "./SE_example2.c"
}   /* End of psd_ba=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'pl_mon_ba'. */
  SIG_MESSAGE("pl_mon_ba (McDisplay)");
  printf("MCDISPLAY: component %s\n", "pl_mon_ba");
#define mccompcurname  pl_mon_ba
#define mccompcurtype  MeanPolLambda_monitor
#define mccompcurindex 27
#define xwidth mccpl_mon_ba_xwidth
#define yheight mccpl_mon_ba_yheight
#define nL mccpl_mon_ba_nL
#define restore_neutron mccpl_mon_ba_restore_neutron
#define nowritefile mccpl_mon_ba_nowritefile
#define PolL_N mccpl_mon_ba_PolL_N
#define PolL_p mccpl_mon_ba_PolL_p
#define PolL_p2 mccpl_mon_ba_PolL_p2
#define HelpArray mccpl_mon_ba_HelpArray
{   /* Declarations of pl_mon_ba=MeanPolLambda_monitor() SETTING parameters. */
char* filename = mccpl_mon_ba_filename;
MCNUM mx = mccpl_mon_ba_mx;
MCNUM my = mccpl_mon_ba_my;
MCNUM mz = mccpl_mon_ba_mz;
MCNUM Lmin = mccpl_mon_ba_Lmin;
MCNUM Lmax = mccpl_mon_ba_Lmax;
#line 163 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/MeanPolLambda_monitor.comp"
{
  
  rectangle("xy", 0, 0, 0, xwidth, yheight);
}
#line 20074 "./SE_example2.c"
}   /* End of pl_mon_ba=MeanPolLambda_monitor() SETTING parameter declarations. */
#undef HelpArray
#undef PolL_p2
#undef PolL_p
#undef PolL_N
#undef nowritefile
#undef restore_neutron
#undef nL
#undef yheight
#undef xwidth
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'analyzer'. */
  SIG_MESSAGE("analyzer (McDisplay)");
  printf("MCDISPLAY: component %s\n", "analyzer");
#define mccompcurname  analyzer
#define mccompcurtype  Pol_mirror
#define mccompcurindex 28
#define rUpFunc mccanalyzer_rUpFunc
#define rDownFunc mccanalyzer_rDownFunc
#define rUpPar mccanalyzer_rUpPar
#define rDownPar mccanalyzer_rDownPar
#define useTables mccanalyzer_useTables
#define zwidth mccanalyzer_zwidth
#define yheight mccanalyzer_yheight
#define rUpParPtr mccanalyzer_rUpParPtr
#define rDownParPtr mccanalyzer_rDownParPtr
{   /* Declarations of analyzer=Pol_mirror() SETTING parameters. */
MCNUM p_reflect = mccanalyzer_p_reflect;
#line 208 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Pol_mirror.comp"
{
  
  rectangle("yz", 0, 0, 0, zwidth, yheight);
}
#line 20111 "./SE_example2.c"
}   /* End of analyzer=Pol_mirror() SETTING parameter declarations. */
#undef rDownParPtr
#undef rUpParPtr
#undef yheight
#undef zwidth
#undef useTables
#undef rDownPar
#undef rUpPar
#undef rDownFunc
#undef rUpFunc
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'detector_slit'. */
  SIG_MESSAGE("detector_slit (McDisplay)");
  printf("MCDISPLAY: component %s\n", "detector_slit");
#define mccompcurname  detector_slit
#define mccompcurtype  Slit
#define mccompcurindex 29
{   /* Declarations of detector_slit=Slit() SETTING parameters. */
MCNUM xmin = mccdetector_slit_xmin;
MCNUM xmax = mccdetector_slit_xmax;
MCNUM ymin = mccdetector_slit_ymin;
MCNUM ymax = mccdetector_slit_ymax;
MCNUM radius = mccdetector_slit_radius;
MCNUM xwidth = mccdetector_slit_xwidth;
MCNUM yheight = mccdetector_slit_yheight;
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
#line 20163 "./SE_example2.c"
}   /* End of detector_slit=Slit() SETTING parameter declarations. */
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'pl_mon5'. */
  SIG_MESSAGE("pl_mon5 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "pl_mon5");
#define mccompcurname  pl_mon5
#define mccompcurtype  MeanPolLambda_monitor
#define mccompcurindex 30
#define xwidth mccpl_mon5_xwidth
#define yheight mccpl_mon5_yheight
#define nL mccpl_mon5_nL
#define restore_neutron mccpl_mon5_restore_neutron
#define nowritefile mccpl_mon5_nowritefile
#define PolL_N mccpl_mon5_PolL_N
#define PolL_p mccpl_mon5_PolL_p
#define PolL_p2 mccpl_mon5_PolL_p2
#define HelpArray mccpl_mon5_HelpArray
{   /* Declarations of pl_mon5=MeanPolLambda_monitor() SETTING parameters. */
char* filename = mccpl_mon5_filename;
MCNUM mx = mccpl_mon5_mx;
MCNUM my = mccpl_mon5_my;
MCNUM mz = mccpl_mon5_mz;
MCNUM Lmin = mccpl_mon5_Lmin;
MCNUM Lmax = mccpl_mon5_Lmax;
#line 163 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/MeanPolLambda_monitor.comp"
{
  
  rectangle("xy", 0, 0, 0, xwidth, yheight);
}
#line 20196 "./SE_example2.c"
}   /* End of pl_mon5=MeanPolLambda_monitor() SETTING parameter declarations. */
#undef HelpArray
#undef PolL_p2
#undef PolL_p
#undef PolL_N
#undef nowritefile
#undef restore_neutron
#undef nL
#undef yheight
#undef xwidth
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'detector'. */
  SIG_MESSAGE("detector (McDisplay)");
  printf("MCDISPLAY: component %s\n", "detector");
#define mccompcurname  detector
#define mccompcurtype  PSD_monitor
#define mccompcurindex 31
#define PSD_N mccdetector_PSD_N
#define PSD_p mccdetector_PSD_p
#define PSD_p2 mccdetector_PSD_p2
{   /* Declarations of detector=PSD_monitor() SETTING parameters. */
int nx = mccdetector_nx;
int ny = mccdetector_ny;
char* filename = mccdetector_filename;
MCNUM xmin = mccdetector_xmin;
MCNUM xmax = mccdetector_xmax;
MCNUM ymin = mccdetector_ymin;
MCNUM ymax = mccdetector_ymax;
MCNUM xwidth = mccdetector_xwidth;
MCNUM yheight = mccdetector_yheight;
MCNUM restore_neutron = mccdetector_restore_neutron;
#line 129 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
{
  multiline(5,
    (double)xmin, (double)ymin, 0.0,
    (double)xmax, (double)ymin, 0.0,
    (double)xmax, (double)ymax, 0.0,
    (double)xmin, (double)ymax, 0.0,
    (double)xmin, (double)ymin, 0.0);
}
#line 20240 "./SE_example2.c"
}   /* End of detector=PSD_monitor() SETTING parameter declarations. */
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
/* end of generated C code ./SE_example2.c */
