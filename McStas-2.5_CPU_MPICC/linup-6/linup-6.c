/* Automatically generated file. Do not edit. 
 * Format:     ANSI C source code
 * Creator:    McStas <http://www.mcstas.org>
 * Instrument: /zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr (TAS1_Vana)
 * Date:       Wed Nov 20 00:57:24 2019
 * File:       ./linup-6.c
 * Compile:    cc -o TAS1_Vana.out ./linup-6.c 
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

#line 712 "./linup-6.c"

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

#line 945 "./linup-6.c"

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

#line 4977 "./linup-6.c"

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

#line 5337 "./linup-6.c"
#ifdef MC_TRACE_ENABLED
int mctraceenabled = 1;
#else
int mctraceenabled = 0;
#endif
#define MCSTAS "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../"
int mcdefaultmain = 1;
char mcinstrument_name[] = "TAS1_Vana";
char mcinstrument_source[] = "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr";
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
#line 5377 "./linup-6.c"

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
#line 5394 "./linup-6.c"

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

#line 5480 "./linup-6.c"

/* Instrument parameters. */
MCNUM mcipPHM;
MCNUM mcipTTM;
MCNUM mcipTT;
MCNUM mcipTTA;
MCNUM mcipC1;
MCNUM mcipOMC1;
MCNUM mcipC2;
MCNUM mcipC3;
MCNUM mcipOMA;

#define mcNUMIPAR 9
int mcnumipar = 9;
struct mcinputtable_struct mcinputtable[mcNUMIPAR+1] = {
  "PHM", &mcipPHM, instr_type_double, "-37.077", 
  "TTM", &mcipTTM, instr_type_double, "-74", 
  "TT", &mcipTT, instr_type_double, "33.52", 
  "TTA", &mcipTTA, instr_type_double, "0", 
  "C1", &mcipC1, instr_type_double, "30", 
  "OMC1", &mcipOMC1, instr_type_double, "5.5", 
  "C2", &mcipC2, instr_type_double, "28", 
  "C3", &mcipC3, instr_type_double, "67", 
  "OMA", &mcipOMA, instr_type_double, "-17.45", 
  NULL, NULL, instr_type_double, ""
};

/* User declarations from instrument definition. */
#define mccompcurname  TAS1_Vana
#define mccompcurtype  INSTRUMENT
#define mccompcurindex 0
#define mcposaTAS1_Vana coords_set(0,0,0)
#define PHM mcipPHM
#define TTM mcipTTM
#define TT mcipTT
#define TTA mcipTTA
#define C1 mcipC1
#define OMC1 mcipOMC1
#define C2 mcipC2
#define C3 mcipC3
#define OMA mcipOMA
#line 48 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
/* Mosaicity used on monochromator and analysator */
double tas1_mono_mosaic = 45; /* Measurements indicate its really 45' */
double tas1_ana_mosaic = 45;  /* Measurements indicate its really 45' */
/* Q vector for bragg scattering with monochromator and analysator */
double tas1_mono_q = 2*1.87325; /* Fake 2nd order scattering for 20meV */
double tas1_mono_r0 = 0.6;
double tas1_ana_q = 1.87325;  /* 20meV */
double tas1_ana_r0 = 0.6;

double OMC1_d;
double alu_focus_x;

double mpos0, mpos1, mpos2, mpos3, mpos4, mpos5, mpos6, mpos7;
double mrot0, mrot1, mrot2, mrot3, mrot4, mrot5, mrot6, mrot7;
#line 5537 "./linup-6.c"
#undef OMA
#undef C3
#undef C2
#undef OMC1
#undef C1
#undef TTA
#undef TT
#undef TTM
#undef PHM
#undef mcposaTAS1_Vana
#undef mccompcurindex
#undef mccompcurtype
#undef mccompcurname

/* neutron state table at each component input (local coords) */
/* [x, y, z, vx, vy, vz, t, sx, sy, sz, p] */
MCNUM mccomp_storein[11*34];
/* Components position table (absolute and relative coords) */
Coords mccomp_posa[34];
Coords mccomp_posr[34];
/* Counter for each comp to check for inactive ones */
MCNUM  mcNCounter[34];
MCNUM  mcPCounter[34];
MCNUM  mcP2Counter[34];
#define mcNUMCOMP 33 /* number of components */
/* Counter for PROP ABSORB */
MCNUM  mcAbsorbProp[34];
/* Flag true when previous component acted on the neutron (SCATTER) */
MCNUM mcScattered=0;
/* Flag true when neutron should be restored (RESTORE) */
MCNUM mcRestore=0;
/* Declarations of component definition and setting parameters. */

/* Setting parameters for component 'a1' [1]. */
char mcca1_profile[16384];
MCNUM mcca1_percent;
MCNUM mcca1_flag_save;
MCNUM mcca1_minutes;

/* Setting parameters for component 'source' [2]. */
MCNUM mccsource_radius;
MCNUM mccsource_yheight;
MCNUM mccsource_xwidth;
MCNUM mccsource_dist;
MCNUM mccsource_focus_xw;
MCNUM mccsource_focus_yh;
MCNUM mccsource_E0;
MCNUM mccsource_dE;
MCNUM mccsource_lambda0;
MCNUM mccsource_dlambda;
MCNUM mccsource_flux;
MCNUM mccsource_gauss;
int mccsource_target_index;

/* Setting parameters for component 'slit1' [3]. */
MCNUM mccslit1_xmin;
MCNUM mccslit1_xmax;
MCNUM mccslit1_ymin;
MCNUM mccslit1_ymax;
MCNUM mccslit1_radius;
MCNUM mccslit1_xwidth;
MCNUM mccslit1_yheight;

/* Setting parameters for component 'slit2' [4]. */
MCNUM mccslit2_xmin;
MCNUM mccslit2_xmax;
MCNUM mccslit2_ymin;
MCNUM mccslit2_ymax;
MCNUM mccslit2_radius;
MCNUM mccslit2_xwidth;
MCNUM mccslit2_yheight;

/* Setting parameters for component 'slit3' [5]. */
MCNUM mccslit3_xmin;
MCNUM mccslit3_xmax;
MCNUM mccslit3_ymin;
MCNUM mccslit3_ymax;
MCNUM mccslit3_radius;
MCNUM mccslit3_xwidth;
MCNUM mccslit3_yheight;

/* Setting parameters for component 'm0' [7]. */
MCNUM mccm0_zmin;
MCNUM mccm0_zmax;
MCNUM mccm0_ymin;
MCNUM mccm0_ymax;
MCNUM mccm0_zwidth;
MCNUM mccm0_yheight;
MCNUM mccm0_mosaich;
MCNUM mccm0_mosaicv;
MCNUM mccm0_r0;
MCNUM mccm0_Q;
MCNUM mccm0_DM;

/* Setting parameters for component 'm1' [8]. */
MCNUM mccm1_zmin;
MCNUM mccm1_zmax;
MCNUM mccm1_ymin;
MCNUM mccm1_ymax;
MCNUM mccm1_zwidth;
MCNUM mccm1_yheight;
MCNUM mccm1_mosaich;
MCNUM mccm1_mosaicv;
MCNUM mccm1_r0;
MCNUM mccm1_Q;
MCNUM mccm1_DM;

/* Setting parameters for component 'm2' [9]. */
MCNUM mccm2_zmin;
MCNUM mccm2_zmax;
MCNUM mccm2_ymin;
MCNUM mccm2_ymax;
MCNUM mccm2_zwidth;
MCNUM mccm2_yheight;
MCNUM mccm2_mosaich;
MCNUM mccm2_mosaicv;
MCNUM mccm2_r0;
MCNUM mccm2_Q;
MCNUM mccm2_DM;

/* Setting parameters for component 'm3' [10]. */
MCNUM mccm3_zmin;
MCNUM mccm3_zmax;
MCNUM mccm3_ymin;
MCNUM mccm3_ymax;
MCNUM mccm3_zwidth;
MCNUM mccm3_yheight;
MCNUM mccm3_mosaich;
MCNUM mccm3_mosaicv;
MCNUM mccm3_r0;
MCNUM mccm3_Q;
MCNUM mccm3_DM;

/* Setting parameters for component 'm4' [11]. */
MCNUM mccm4_zmin;
MCNUM mccm4_zmax;
MCNUM mccm4_ymin;
MCNUM mccm4_ymax;
MCNUM mccm4_zwidth;
MCNUM mccm4_yheight;
MCNUM mccm4_mosaich;
MCNUM mccm4_mosaicv;
MCNUM mccm4_r0;
MCNUM mccm4_Q;
MCNUM mccm4_DM;

/* Setting parameters for component 'm5' [12]. */
MCNUM mccm5_zmin;
MCNUM mccm5_zmax;
MCNUM mccm5_ymin;
MCNUM mccm5_ymax;
MCNUM mccm5_zwidth;
MCNUM mccm5_yheight;
MCNUM mccm5_mosaich;
MCNUM mccm5_mosaicv;
MCNUM mccm5_r0;
MCNUM mccm5_Q;
MCNUM mccm5_DM;

/* Setting parameters for component 'm6' [13]. */
MCNUM mccm6_zmin;
MCNUM mccm6_zmax;
MCNUM mccm6_ymin;
MCNUM mccm6_ymax;
MCNUM mccm6_zwidth;
MCNUM mccm6_yheight;
MCNUM mccm6_mosaich;
MCNUM mccm6_mosaicv;
MCNUM mccm6_r0;
MCNUM mccm6_Q;
MCNUM mccm6_DM;

/* Setting parameters for component 'm7' [14]. */
MCNUM mccm7_zmin;
MCNUM mccm7_zmax;
MCNUM mccm7_ymin;
MCNUM mccm7_ymax;
MCNUM mccm7_zwidth;
MCNUM mccm7_yheight;
MCNUM mccm7_mosaich;
MCNUM mccm7_mosaicv;
MCNUM mccm7_r0;
MCNUM mccm7_Q;
MCNUM mccm7_DM;

/* Setting parameters for component 'slitMS1' [16]. */
MCNUM mccslitMS1_xmin;
MCNUM mccslitMS1_xmax;
MCNUM mccslitMS1_ymin;
MCNUM mccslitMS1_ymax;
MCNUM mccslitMS1_radius;
MCNUM mccslitMS1_xwidth;
MCNUM mccslitMS1_yheight;

/* Setting parameters for component 'slitMS2' [17]. */
MCNUM mccslitMS2_xmin;
MCNUM mccslitMS2_xmax;
MCNUM mccslitMS2_ymin;
MCNUM mccslitMS2_ymax;
MCNUM mccslitMS2_radius;
MCNUM mccslitMS2_xwidth;
MCNUM mccslitMS2_yheight;

/* Setting parameters for component 'c1' [18]. */
MCNUM mccc1_xmin;
MCNUM mccc1_xmax;
MCNUM mccc1_ymin;
MCNUM mccc1_ymax;
MCNUM mccc1_xwidth;
MCNUM mccc1_yheight;
MCNUM mccc1_length;
MCNUM mccc1_divergence;
MCNUM mccc1_transmission;
MCNUM mccc1_divergenceV;

/* Setting parameters for component 'slitMS3' [19]. */
MCNUM mccslitMS3_xmin;
MCNUM mccslitMS3_xmax;
MCNUM mccslitMS3_ymin;
MCNUM mccslitMS3_ymax;
MCNUM mccslitMS3_radius;
MCNUM mccslitMS3_xwidth;
MCNUM mccslitMS3_yheight;

/* Setting parameters for component 'slitMS4' [20]. */
MCNUM mccslitMS4_xmin;
MCNUM mccslitMS4_xmax;
MCNUM mccslitMS4_ymin;
MCNUM mccslitMS4_ymax;
MCNUM mccslitMS4_radius;
MCNUM mccslitMS4_xwidth;
MCNUM mccslitMS4_yheight;

/* Setting parameters for component 'slitMS5' [21]. */
MCNUM mccslitMS5_xmin;
MCNUM mccslitMS5_xmax;
MCNUM mccslitMS5_ymin;
MCNUM mccslitMS5_ymax;
MCNUM mccslitMS5_radius;
MCNUM mccslitMS5_xwidth;
MCNUM mccslitMS5_yheight;

/* Setting parameters for component 'mon' [22]. */
MCNUM mccmon_xmin;
MCNUM mccmon_xmax;
MCNUM mccmon_ymin;
MCNUM mccmon_ymax;
MCNUM mccmon_xwidth;
MCNUM mccmon_yheight;
MCNUM mccmon_restore_neutron;

/* Definition parameters for component 'emon1' [23]. */
#define mccemon1_nE 35
/* Setting parameters for component 'emon1' [23]. */
char mccemon1_filename[16384];
MCNUM mccemon1_xmin;
MCNUM mccemon1_xmax;
MCNUM mccemon1_ymin;
MCNUM mccemon1_ymax;
MCNUM mccemon1_xwidth;
MCNUM mccemon1_yheight;
MCNUM mccemon1_Emin;
MCNUM mccemon1_Emax;
MCNUM mccemon1_restore_neutron;
int mccemon1_nowritefile;

/* Setting parameters for component 'sample' [24]. */
MCNUM mccsample_radius;
MCNUM mccsample_thickness;
MCNUM mccsample_zdepth;
MCNUM mccsample_Vc;
MCNUM mccsample_sigma_abs;
MCNUM mccsample_sigma_inc;
MCNUM mccsample_radius_i;
MCNUM mccsample_radius_o;
MCNUM mccsample_h;
MCNUM mccsample_focus_r;
MCNUM mccsample_pack;
MCNUM mccsample_frac;
MCNUM mccsample_f_QE;
MCNUM mccsample_gamma;
MCNUM mccsample_target_x;
MCNUM mccsample_target_y;
MCNUM mccsample_target_z;
MCNUM mccsample_focus_xw;
MCNUM mccsample_focus_yh;
MCNUM mccsample_focus_aw;
MCNUM mccsample_focus_ah;
MCNUM mccsample_xwidth;
MCNUM mccsample_yheight;
MCNUM mccsample_zthick;
MCNUM mccsample_rad_sphere;
MCNUM mccsample_sig_a;
MCNUM mccsample_sig_i;
MCNUM mccsample_V0;
int mccsample_target_index;
MCNUM mccsample_multiples;

/* Setting parameters for component 'focus_check' [26]. */
int mccfocus_check_nx;
int mccfocus_check_ny;
char mccfocus_check_filename[16384];
MCNUM mccfocus_check_xmin;
MCNUM mccfocus_check_xmax;
MCNUM mccfocus_check_ymin;
MCNUM mccfocus_check_ymax;
MCNUM mccfocus_check_xwidth;
MCNUM mccfocus_check_yheight;
MCNUM mccfocus_check_restore_neutron;

/* Setting parameters for component 'c2' [27]. */
MCNUM mccc2_xmin;
MCNUM mccc2_xmax;
MCNUM mccc2_ymin;
MCNUM mccc2_ymax;
MCNUM mccc2_xwidth;
MCNUM mccc2_yheight;
MCNUM mccc2_length;
MCNUM mccc2_divergence;
MCNUM mccc2_transmission;
MCNUM mccc2_divergenceV;

/* Setting parameters for component 'ana' [28]. */
MCNUM mccana_zmin;
MCNUM mccana_zmax;
MCNUM mccana_ymin;
MCNUM mccana_ymax;
MCNUM mccana_zwidth;
MCNUM mccana_yheight;
MCNUM mccana_mosaich;
MCNUM mccana_mosaicv;
MCNUM mccana_r0;
MCNUM mccana_Q;
MCNUM mccana_DM;

/* Setting parameters for component 'c3' [30]. */
MCNUM mccc3_xmin;
MCNUM mccc3_xmax;
MCNUM mccc3_ymin;
MCNUM mccc3_ymax;
MCNUM mccc3_xwidth;
MCNUM mccc3_yheight;
MCNUM mccc3_length;
MCNUM mccc3_divergence;
MCNUM mccc3_transmission;
MCNUM mccc3_divergenceV;

/* Setting parameters for component 'sng' [31]. */
MCNUM mccsng_xmin;
MCNUM mccsng_xmax;
MCNUM mccsng_ymin;
MCNUM mccsng_ymax;
MCNUM mccsng_xwidth;
MCNUM mccsng_yheight;
MCNUM mccsng_restore_neutron;

/* Definition parameters for component 'emon2' [32]. */
#define mccemon2_nE 35
/* Setting parameters for component 'emon2' [32]. */
char mccemon2_filename[16384];
MCNUM mccemon2_xmin;
MCNUM mccemon2_xmax;
MCNUM mccemon2_ymin;
MCNUM mccemon2_ymax;
MCNUM mccemon2_xwidth;
MCNUM mccemon2_yheight;
MCNUM mccemon2_Emin;
MCNUM mccemon2_Emax;
MCNUM mccemon2_restore_neutron;
int mccemon2_nowritefile;

/* User component declarations. */

/* User declarations for component 'a1' [1]. */
#define mccompcurname  a1
#define mccompcurtype  Progress_bar
#define mccompcurindex 1
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
#line 5934 "./linup-6.c"
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

/* User declarations for component 'source' [2]. */
#define mccompcurname  source
#define mccompcurtype  Source_simple
#define mccompcurindex 2
#define pmul mccsource_pmul
#define square mccsource_square
#define srcArea mccsource_srcArea
#define radius mccsource_radius
#define yheight mccsource_yheight
#define xwidth mccsource_xwidth
#define dist mccsource_dist
#define focus_xw mccsource_focus_xw
#define focus_yh mccsource_focus_yh
#define E0 mccsource_E0
#define dE mccsource_dE
#define lambda0 mccsource_lambda0
#define dlambda mccsource_dlambda
#define flux mccsource_flux
#define gauss mccsource_gauss
#define target_index mccsource_target_index
#line 60 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../sources/Source_simple.comp"
double pmul, srcArea;
int square;
double tx,ty,tz;
#line 5971 "./linup-6.c"
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

/* User declarations for component 'slit1' [3]. */
#define mccompcurname  slit1
#define mccompcurtype  Slit
#define mccompcurindex 3
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

/* User declarations for component 'slit2' [4]. */
#define mccompcurname  slit2
#define mccompcurtype  Slit
#define mccompcurindex 4
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

/* User declarations for component 'slit3' [5]. */
#define mccompcurname  slit3
#define mccompcurtype  Slit
#define mccompcurindex 5
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

/* User declarations for component 'focus_mono' [6]. */
#define mccompcurname  focus_mono
#define mccompcurtype  Arm
#define mccompcurindex 6
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'm0' [7]. */
#define mccompcurname  m0
#define mccompcurtype  Monochromator_flat
#define mccompcurindex 7
#define mos_rms_y mccm0_mos_rms_y
#define mos_rms_z mccm0_mos_rms_z
#define mos_rms_max mccm0_mos_rms_max
#define mono_Q mccm0_mono_Q
#define zmin mccm0_zmin
#define zmax mccm0_zmax
#define ymin mccm0_ymin
#define ymax mccm0_ymax
#define zwidth mccm0_zwidth
#define yheight mccm0_yheight
#define mosaich mccm0_mosaich
#define mosaicv mccm0_mosaicv
#define r0 mccm0_r0
#define Q mccm0_Q
#define DM mccm0_DM
#line 95 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Monochromator_flat.comp"
  double mos_rms_y; /* root-mean-square of mosaic, in radians */
  double mos_rms_z;
  double mos_rms_max;
  double mono_Q;
#line 6090 "./linup-6.c"
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

/* User declarations for component 'm1' [8]. */
#define mccompcurname  m1
#define mccompcurtype  Monochromator_flat
#define mccompcurindex 8
#define mos_rms_y mccm1_mos_rms_y
#define mos_rms_z mccm1_mos_rms_z
#define mos_rms_max mccm1_mos_rms_max
#define mono_Q mccm1_mono_Q
#define zmin mccm1_zmin
#define zmax mccm1_zmax
#define ymin mccm1_ymin
#define ymax mccm1_ymax
#define zwidth mccm1_zwidth
#define yheight mccm1_yheight
#define mosaich mccm1_mosaich
#define mosaicv mccm1_mosaicv
#define r0 mccm1_r0
#define Q mccm1_Q
#define DM mccm1_DM
#line 95 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Monochromator_flat.comp"
  double mos_rms_y; /* root-mean-square of mosaic, in radians */
  double mos_rms_z;
  double mos_rms_max;
  double mono_Q;
#line 6134 "./linup-6.c"
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

/* User declarations for component 'm2' [9]. */
#define mccompcurname  m2
#define mccompcurtype  Monochromator_flat
#define mccompcurindex 9
#define mos_rms_y mccm2_mos_rms_y
#define mos_rms_z mccm2_mos_rms_z
#define mos_rms_max mccm2_mos_rms_max
#define mono_Q mccm2_mono_Q
#define zmin mccm2_zmin
#define zmax mccm2_zmax
#define ymin mccm2_ymin
#define ymax mccm2_ymax
#define zwidth mccm2_zwidth
#define yheight mccm2_yheight
#define mosaich mccm2_mosaich
#define mosaicv mccm2_mosaicv
#define r0 mccm2_r0
#define Q mccm2_Q
#define DM mccm2_DM
#line 95 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Monochromator_flat.comp"
  double mos_rms_y; /* root-mean-square of mosaic, in radians */
  double mos_rms_z;
  double mos_rms_max;
  double mono_Q;
#line 6178 "./linup-6.c"
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

/* User declarations for component 'm3' [10]. */
#define mccompcurname  m3
#define mccompcurtype  Monochromator_flat
#define mccompcurindex 10
#define mos_rms_y mccm3_mos_rms_y
#define mos_rms_z mccm3_mos_rms_z
#define mos_rms_max mccm3_mos_rms_max
#define mono_Q mccm3_mono_Q
#define zmin mccm3_zmin
#define zmax mccm3_zmax
#define ymin mccm3_ymin
#define ymax mccm3_ymax
#define zwidth mccm3_zwidth
#define yheight mccm3_yheight
#define mosaich mccm3_mosaich
#define mosaicv mccm3_mosaicv
#define r0 mccm3_r0
#define Q mccm3_Q
#define DM mccm3_DM
#line 95 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Monochromator_flat.comp"
  double mos_rms_y; /* root-mean-square of mosaic, in radians */
  double mos_rms_z;
  double mos_rms_max;
  double mono_Q;
#line 6222 "./linup-6.c"
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

/* User declarations for component 'm4' [11]. */
#define mccompcurname  m4
#define mccompcurtype  Monochromator_flat
#define mccompcurindex 11
#define mos_rms_y mccm4_mos_rms_y
#define mos_rms_z mccm4_mos_rms_z
#define mos_rms_max mccm4_mos_rms_max
#define mono_Q mccm4_mono_Q
#define zmin mccm4_zmin
#define zmax mccm4_zmax
#define ymin mccm4_ymin
#define ymax mccm4_ymax
#define zwidth mccm4_zwidth
#define yheight mccm4_yheight
#define mosaich mccm4_mosaich
#define mosaicv mccm4_mosaicv
#define r0 mccm4_r0
#define Q mccm4_Q
#define DM mccm4_DM
#line 95 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Monochromator_flat.comp"
  double mos_rms_y; /* root-mean-square of mosaic, in radians */
  double mos_rms_z;
  double mos_rms_max;
  double mono_Q;
#line 6266 "./linup-6.c"
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

/* User declarations for component 'm5' [12]. */
#define mccompcurname  m5
#define mccompcurtype  Monochromator_flat
#define mccompcurindex 12
#define mos_rms_y mccm5_mos_rms_y
#define mos_rms_z mccm5_mos_rms_z
#define mos_rms_max mccm5_mos_rms_max
#define mono_Q mccm5_mono_Q
#define zmin mccm5_zmin
#define zmax mccm5_zmax
#define ymin mccm5_ymin
#define ymax mccm5_ymax
#define zwidth mccm5_zwidth
#define yheight mccm5_yheight
#define mosaich mccm5_mosaich
#define mosaicv mccm5_mosaicv
#define r0 mccm5_r0
#define Q mccm5_Q
#define DM mccm5_DM
#line 95 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Monochromator_flat.comp"
  double mos_rms_y; /* root-mean-square of mosaic, in radians */
  double mos_rms_z;
  double mos_rms_max;
  double mono_Q;
#line 6310 "./linup-6.c"
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

/* User declarations for component 'm6' [13]. */
#define mccompcurname  m6
#define mccompcurtype  Monochromator_flat
#define mccompcurindex 13
#define mos_rms_y mccm6_mos_rms_y
#define mos_rms_z mccm6_mos_rms_z
#define mos_rms_max mccm6_mos_rms_max
#define mono_Q mccm6_mono_Q
#define zmin mccm6_zmin
#define zmax mccm6_zmax
#define ymin mccm6_ymin
#define ymax mccm6_ymax
#define zwidth mccm6_zwidth
#define yheight mccm6_yheight
#define mosaich mccm6_mosaich
#define mosaicv mccm6_mosaicv
#define r0 mccm6_r0
#define Q mccm6_Q
#define DM mccm6_DM
#line 95 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Monochromator_flat.comp"
  double mos_rms_y; /* root-mean-square of mosaic, in radians */
  double mos_rms_z;
  double mos_rms_max;
  double mono_Q;
#line 6354 "./linup-6.c"
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

/* User declarations for component 'm7' [14]. */
#define mccompcurname  m7
#define mccompcurtype  Monochromator_flat
#define mccompcurindex 14
#define mos_rms_y mccm7_mos_rms_y
#define mos_rms_z mccm7_mos_rms_z
#define mos_rms_max mccm7_mos_rms_max
#define mono_Q mccm7_mono_Q
#define zmin mccm7_zmin
#define zmax mccm7_zmax
#define ymin mccm7_ymin
#define ymax mccm7_ymax
#define zwidth mccm7_zwidth
#define yheight mccm7_yheight
#define mosaich mccm7_mosaich
#define mosaicv mccm7_mosaicv
#define r0 mccm7_r0
#define Q mccm7_Q
#define DM mccm7_DM
#line 95 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Monochromator_flat.comp"
  double mos_rms_y; /* root-mean-square of mosaic, in radians */
  double mos_rms_z;
  double mos_rms_max;
  double mono_Q;
#line 6398 "./linup-6.c"
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

/* User declarations for component 'a2' [15]. */
#define mccompcurname  a2
#define mccompcurtype  Arm
#define mccompcurindex 15
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'slitMS1' [16]. */
#define mccompcurname  slitMS1
#define mccompcurtype  Slit
#define mccompcurindex 16
#define xmin mccslitMS1_xmin
#define xmax mccslitMS1_xmax
#define ymin mccslitMS1_ymin
#define ymax mccslitMS1_ymax
#define radius mccslitMS1_radius
#define xwidth mccslitMS1_xwidth
#define yheight mccslitMS1_yheight
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

/* User declarations for component 'slitMS2' [17]. */
#define mccompcurname  slitMS2
#define mccompcurtype  Slit
#define mccompcurindex 17
#define xmin mccslitMS2_xmin
#define xmax mccslitMS2_xmax
#define ymin mccslitMS2_ymin
#define ymax mccslitMS2_ymax
#define radius mccslitMS2_radius
#define xwidth mccslitMS2_xwidth
#define yheight mccslitMS2_yheight
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

/* User declarations for component 'c1' [18]. */
#define mccompcurname  c1
#define mccompcurtype  Collimator_linear
#define mccompcurindex 18
#define slope mccc1_slope
#define slopeV mccc1_slopeV
#define xmin mccc1_xmin
#define xmax mccc1_xmax
#define ymin mccc1_ymin
#define ymax mccc1_ymax
#define xwidth mccc1_xwidth
#define yheight mccc1_yheight
#define length mccc1_length
#define divergence mccc1_divergence
#define transmission mccc1_transmission
#define divergenceV mccc1_divergenceV
#line 51 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Collimator_linear.comp"
double slope, slopeV;
#line 6488 "./linup-6.c"
#undef divergenceV
#undef transmission
#undef divergence
#undef length
#undef yheight
#undef xwidth
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef slopeV
#undef slope
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'slitMS3' [19]. */
#define mccompcurname  slitMS3
#define mccompcurtype  Slit
#define mccompcurindex 19
#define xmin mccslitMS3_xmin
#define xmax mccslitMS3_xmax
#define ymin mccslitMS3_ymin
#define ymax mccslitMS3_ymax
#define radius mccslitMS3_radius
#define xwidth mccslitMS3_xwidth
#define yheight mccslitMS3_yheight
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

/* User declarations for component 'slitMS4' [20]. */
#define mccompcurname  slitMS4
#define mccompcurtype  Slit
#define mccompcurindex 20
#define xmin mccslitMS4_xmin
#define xmax mccslitMS4_xmax
#define ymin mccslitMS4_ymin
#define ymax mccslitMS4_ymax
#define radius mccslitMS4_radius
#define xwidth mccslitMS4_xwidth
#define yheight mccslitMS4_yheight
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

/* User declarations for component 'slitMS5' [21]. */
#define mccompcurname  slitMS5
#define mccompcurtype  Slit
#define mccompcurindex 21
#define xmin mccslitMS5_xmin
#define xmax mccslitMS5_xmax
#define ymin mccslitMS5_ymin
#define ymax mccslitMS5_ymax
#define radius mccslitMS5_radius
#define xwidth mccslitMS5_xwidth
#define yheight mccslitMS5_yheight
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

/* User declarations for component 'mon' [22]. */
#define mccompcurname  mon
#define mccompcurtype  Monitor
#define mccompcurindex 22
#define Nsum mccmon_Nsum
#define psum mccmon_psum
#define p2sum mccmon_p2sum
#define xmin mccmon_xmin
#define xmax mccmon_xmax
#define ymin mccmon_ymin
#define ymax mccmon_ymax
#define xwidth mccmon_xwidth
#define yheight mccmon_yheight
#define restore_neutron mccmon_restore_neutron
#line 52 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/Monitor.comp"
double Nsum;
double psum, p2sum;
#line 6588 "./linup-6.c"
#undef restore_neutron
#undef yheight
#undef xwidth
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef p2sum
#undef psum
#undef Nsum
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'emon1' [23]. */
#define mccompcurname  emon1
#define mccompcurtype  E_monitor
#define mccompcurindex 23
#define nE mccemon1_nE
#define E_N mccemon1_E_N
#define E_p mccemon1_E_p
#define E_p2 mccemon1_E_p2
#define S_p mccemon1_S_p
#define S_pE mccemon1_S_pE
#define S_pE2 mccemon1_S_pE2
#define filename mccemon1_filename
#define xmin mccemon1_xmin
#define xmax mccemon1_xmax
#define ymin mccemon1_ymin
#define ymax mccemon1_ymax
#define xwidth mccemon1_xwidth
#define yheight mccemon1_yheight
#define Emin mccemon1_Emin
#define Emax mccemon1_Emax
#define restore_neutron mccemon1_restore_neutron
#define nowritefile mccemon1_nowritefile
#line 60 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/E_monitor.comp"
double E_N[nE];
double E_p[nE], E_p2[nE];
double S_p, S_pE, S_pE2;
#line 6629 "./linup-6.c"
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

/* User declarations for component 'sample' [24]. */
#define mccompcurname  sample
#define mccompcurtype  V_sample
#define mccompcurindex 24
#define VarsV mccsample_VarsV
#define radius mccsample_radius
#define thickness mccsample_thickness
#define zdepth mccsample_zdepth
#define Vc mccsample_Vc
#define sigma_abs mccsample_sigma_abs
#define sigma_inc mccsample_sigma_inc
#define radius_i mccsample_radius_i
#define radius_o mccsample_radius_o
#define h mccsample_h
#define focus_r mccsample_focus_r
#define pack mccsample_pack
#define frac mccsample_frac
#define f_QE mccsample_f_QE
#define gamma mccsample_gamma
#define target_x mccsample_target_x
#define target_y mccsample_target_y
#define target_z mccsample_target_z
#define focus_xw mccsample_focus_xw
#define focus_yh mccsample_focus_yh
#define focus_aw mccsample_focus_aw
#define focus_ah mccsample_focus_ah
#define xwidth mccsample_xwidth
#define yheight mccsample_yheight
#define zthick mccsample_zthick
#define rad_sphere mccsample_rad_sphere
#define sig_a mccsample_sig_a
#define sig_i mccsample_sig_i
#define V0 mccsample_V0
#define target_index mccsample_target_index
#define multiples mccsample_multiples
#line 117 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../obsolete/V_sample.comp"
  struct StructVarsV VarsV;
#line 6689 "./linup-6.c"
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

/* User declarations for component 'a3' [25]. */
#define mccompcurname  a3
#define mccompcurtype  Arm
#define mccompcurindex 25
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'focus_check' [26]. */
#define mccompcurname  focus_check
#define mccompcurtype  PSD_monitor
#define mccompcurindex 26
#define PSD_N mccfocus_check_PSD_N
#define PSD_p mccfocus_check_PSD_p
#define PSD_p2 mccfocus_check_PSD_p2
#define nx mccfocus_check_nx
#define ny mccfocus_check_ny
#define filename mccfocus_check_filename
#define xmin mccfocus_check_xmin
#define xmax mccfocus_check_xmax
#define ymin mccfocus_check_ymin
#define ymax mccfocus_check_ymax
#define xwidth mccfocus_check_xwidth
#define yheight mccfocus_check_yheight
#define restore_neutron mccfocus_check_restore_neutron
#line 62 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
  DArray2d PSD_N;
  DArray2d PSD_p;
  DArray2d PSD_p2;
#line 6754 "./linup-6.c"
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

/* User declarations for component 'c2' [27]. */
#define mccompcurname  c2
#define mccompcurtype  Collimator_linear
#define mccompcurindex 27
#define slope mccc2_slope
#define slopeV mccc2_slopeV
#define xmin mccc2_xmin
#define xmax mccc2_xmax
#define ymin mccc2_ymin
#define ymax mccc2_ymax
#define xwidth mccc2_xwidth
#define yheight mccc2_yheight
#define length mccc2_length
#define divergence mccc2_divergence
#define transmission mccc2_transmission
#define divergenceV mccc2_divergenceV
#line 51 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Collimator_linear.comp"
double slope, slopeV;
#line 6790 "./linup-6.c"
#undef divergenceV
#undef transmission
#undef divergence
#undef length
#undef yheight
#undef xwidth
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef slopeV
#undef slope
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'ana' [28]. */
#define mccompcurname  ana
#define mccompcurtype  Monochromator_flat
#define mccompcurindex 28
#define mos_rms_y mccana_mos_rms_y
#define mos_rms_z mccana_mos_rms_z
#define mos_rms_max mccana_mos_rms_max
#define mono_Q mccana_mono_Q
#define zmin mccana_zmin
#define zmax mccana_zmax
#define ymin mccana_ymin
#define ymax mccana_ymax
#define zwidth mccana_zwidth
#define yheight mccana_yheight
#define mosaich mccana_mosaich
#define mosaicv mccana_mosaicv
#define r0 mccana_r0
#define Q mccana_Q
#define DM mccana_DM
#line 95 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Monochromator_flat.comp"
  double mos_rms_y; /* root-mean-square of mosaic, in radians */
  double mos_rms_z;
  double mos_rms_max;
  double mono_Q;
#line 6831 "./linup-6.c"
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

/* User declarations for component 'a4' [29]. */
#define mccompcurname  a4
#define mccompcurtype  Arm
#define mccompcurindex 29
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'c3' [30]. */
#define mccompcurname  c3
#define mccompcurtype  Collimator_linear
#define mccompcurindex 30
#define slope mccc3_slope
#define slopeV mccc3_slopeV
#define xmin mccc3_xmin
#define xmax mccc3_xmax
#define ymin mccc3_ymin
#define ymax mccc3_ymax
#define xwidth mccc3_xwidth
#define yheight mccc3_yheight
#define length mccc3_length
#define divergence mccc3_divergence
#define transmission mccc3_transmission
#define divergenceV mccc3_divergenceV
#line 51 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Collimator_linear.comp"
double slope, slopeV;
#line 6877 "./linup-6.c"
#undef divergenceV
#undef transmission
#undef divergence
#undef length
#undef yheight
#undef xwidth
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef slopeV
#undef slope
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'sng' [31]. */
#define mccompcurname  sng
#define mccompcurtype  Monitor
#define mccompcurindex 31
#define Nsum mccsng_Nsum
#define psum mccsng_psum
#define p2sum mccsng_p2sum
#define xmin mccsng_xmin
#define xmax mccsng_xmax
#define ymin mccsng_ymin
#define ymax mccsng_ymax
#define xwidth mccsng_xwidth
#define yheight mccsng_yheight
#define restore_neutron mccsng_restore_neutron
#line 52 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/Monitor.comp"
double Nsum;
double psum, p2sum;
#line 6911 "./linup-6.c"
#undef restore_neutron
#undef yheight
#undef xwidth
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef p2sum
#undef psum
#undef Nsum
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'emon2' [32]. */
#define mccompcurname  emon2
#define mccompcurtype  E_monitor
#define mccompcurindex 32
#define nE mccemon2_nE
#define E_N mccemon2_E_N
#define E_p mccemon2_E_p
#define E_p2 mccemon2_E_p2
#define S_p mccemon2_S_p
#define S_pE mccemon2_S_pE
#define S_pE2 mccemon2_S_pE2
#define filename mccemon2_filename
#define xmin mccemon2_xmin
#define xmax mccemon2_xmax
#define ymin mccemon2_ymin
#define ymax mccemon2_ymax
#define xwidth mccemon2_xwidth
#define yheight mccemon2_yheight
#define Emin mccemon2_Emin
#define Emax mccemon2_Emax
#define restore_neutron mccemon2_restore_neutron
#define nowritefile mccemon2_nowritefile
#line 60 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/E_monitor.comp"
double E_N[nE];
double E_p[nE], E_p2[nE];
double S_p, S_pE, S_pE2;
#line 6952 "./linup-6.c"
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

Coords mcposaa1, mcposra1;
Rotation mcrotaa1, mcrotra1;
Coords mcposasource, mcposrsource;
Rotation mcrotasource, mcrotrsource;
Coords mcposaslit1, mcposrslit1;
Rotation mcrotaslit1, mcrotrslit1;
Coords mcposaslit2, mcposrslit2;
Rotation mcrotaslit2, mcrotrslit2;
Coords mcposaslit3, mcposrslit3;
Rotation mcrotaslit3, mcrotrslit3;
Coords mcposafocus_mono, mcposrfocus_mono;
Rotation mcrotafocus_mono, mcrotrfocus_mono;
Coords mcposam0, mcposrm0;
Rotation mcrotam0, mcrotrm0;
Coords mcposam1, mcposrm1;
Rotation mcrotam1, mcrotrm1;
Coords mcposam2, mcposrm2;
Rotation mcrotam2, mcrotrm2;
Coords mcposam3, mcposrm3;
Rotation mcrotam3, mcrotrm3;
Coords mcposam4, mcposrm4;
Rotation mcrotam4, mcrotrm4;
Coords mcposam5, mcposrm5;
Rotation mcrotam5, mcrotrm5;
Coords mcposam6, mcposrm6;
Rotation mcrotam6, mcrotrm6;
Coords mcposam7, mcposrm7;
Rotation mcrotam7, mcrotrm7;
Coords mcposaa2, mcposra2;
Rotation mcrotaa2, mcrotra2;
Coords mcposaslitMS1, mcposrslitMS1;
Rotation mcrotaslitMS1, mcrotrslitMS1;
Coords mcposaslitMS2, mcposrslitMS2;
Rotation mcrotaslitMS2, mcrotrslitMS2;
Coords mcposac1, mcposrc1;
Rotation mcrotac1, mcrotrc1;
Coords mcposaslitMS3, mcposrslitMS3;
Rotation mcrotaslitMS3, mcrotrslitMS3;
Coords mcposaslitMS4, mcposrslitMS4;
Rotation mcrotaslitMS4, mcrotrslitMS4;
Coords mcposaslitMS5, mcposrslitMS5;
Rotation mcrotaslitMS5, mcrotrslitMS5;
Coords mcposamon, mcposrmon;
Rotation mcrotamon, mcrotrmon;
Coords mcposaemon1, mcposremon1;
Rotation mcrotaemon1, mcrotremon1;
Coords mcposasample, mcposrsample;
Rotation mcrotasample, mcrotrsample;
Coords mcposaa3, mcposra3;
Rotation mcrotaa3, mcrotra3;
Coords mcposafocus_check, mcposrfocus_check;
Rotation mcrotafocus_check, mcrotrfocus_check;
Coords mcposac2, mcposrc2;
Rotation mcrotac2, mcrotrc2;
Coords mcposaana, mcposrana;
Rotation mcrotaana, mcrotrana;
Coords mcposaa4, mcposra4;
Rotation mcrotaa4, mcrotra4;
Coords mcposac3, mcposrc3;
Rotation mcrotac3, mcrotrc3;
Coords mcposasng, mcposrsng;
Rotation mcrotasng, mcrotrsng;
Coords mcposaemon2, mcposremon2;
Rotation mcrotaemon2, mcrotremon2;

MCNUM mcnx, mcny, mcnz, mcnvx, mcnvy, mcnvz, mcnt, mcnsx, mcnsy, mcnsz, mcnp;

/* end declare */

void mcinit(void) {
#define mccompcurname  TAS1_Vana
#define mccompcurtype  INSTRUMENT
#define mccompcurindex 0
#define mcposaTAS1_Vana coords_set(0,0,0)
#define PHM mcipPHM
#define TTM mcipTTM
#define TT mcipTT
#define TTA mcipTTA
#define C1 mcipC1
#define OMC1 mcipOMC1
#define C2 mcipC2
#define C3 mcipC3
#define OMA mcipOMA
#line 65 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
{
  double d = 0.0125;    /* 12.5 mm between slab centers. */
  double phi = 0.5443;    /* Rotation between adjacent slabs. */
  mpos0 = -3.5*d; mrot0 = -3.5*phi;
  mpos1 = -2.5*d; mrot1 = -2.5*phi;
  mpos2 = -1.5*d; mrot2 = -1.5*phi;
  mpos3 = -0.5*d; mrot3 = -0.5*phi;
  mpos4 =  0.5*d; mrot4 =  0.5*phi;
  mpos5 =  1.5*d; mrot5 =  1.5*phi;
  mpos6 =  2.5*d; mrot6 =  2.5*phi;
  mpos7 =  3.5*d; mrot7 =  3.5*phi;

  OMC1_d = OMC1/60.0;
  alu_focus_x = TT >= 0 ? 1000 : -1000;
}
#line 7074 "./linup-6.c"
#undef OMA
#undef C3
#undef C2
#undef OMC1
#undef C1
#undef TTA
#undef TT
#undef TTM
#undef PHM
#undef mcposaTAS1_Vana
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
#line 39 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  if("NULL") strncpy(mcca1_profile, "NULL" ? "NULL" : "", 16384); else mcca1_profile[0]='\0';
#line 39 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mcca1_percent = 10;
#line 39 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mcca1_flag_save = 0;
#line 39 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mcca1_minutes = 0;
#line 7109 "./linup-6.c"

  SIG_MESSAGE("a1 (Init:Place/Rotate)");
  rot_set_rotation(mcrotaa1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 7116 "./linup-6.c"
  rot_copy(mcrotra1, mcrotaa1);
  mcposaa1 = coords_set(
#line 84 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    0,
#line 84 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    0,
#line 84 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    0);
#line 7125 "./linup-6.c"
  mctc1 = coords_neg(mcposaa1);
  mcposra1 = rot_apply(mcrotaa1, mctc1);
  mcDEBUG_COMPONENT("a1", mcposaa1, mcrotaa1)
  mccomp_posa[1] = mcposaa1;
  mccomp_posr[1] = mcposra1;
  mcNCounter[1]  = mcPCounter[1] = mcP2Counter[1] = 0;
  mcAbsorbProp[1]= 0;
    /* Component source. */
  /* Setting parameters for component source. */
  SIG_MESSAGE("source (Init:SetPar)");
#line 87 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccsource_radius = 0.060;
#line 52 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccsource_yheight = 0;
#line 52 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccsource_xwidth = 0;
#line 88 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccsource_dist = 3.288;
#line 89 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccsource_focus_xw = 0.042;
#line 89 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccsource_focus_yh = 0.082;
#line 90 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccsource_E0 = 20;
#line 91 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccsource_dE = 0.82;
#line 54 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccsource_lambda0 = 0;
#line 54 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccsource_dlambda = 0;
#line 55 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccsource_flux = 1;
#line 55 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccsource_gauss = 0;
#line 55 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccsource_target_index = + 1;
#line 7162 "./linup-6.c"

  SIG_MESSAGE("source (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 93 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    (0)*DEG2RAD,
#line 93 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    (0)*DEG2RAD,
#line 93 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    (0)*DEG2RAD);
#line 7172 "./linup-6.c"
  rot_mul(mctr1, mcrotaa1, mcrotasource);
  rot_transpose(mcrotaa1, mctr1);
  rot_mul(mcrotasource, mctr1, mcrotrsource);
  mctc1 = coords_set(
#line 93 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    0,
#line 93 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    0,
#line 93 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    0);
#line 7183 "./linup-6.c"
  rot_transpose(mcrotaa1, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposasource = coords_add(mcposaa1, mctc2);
  mctc1 = coords_sub(mcposaa1, mcposasource);
  mcposrsource = rot_apply(mcrotasource, mctc1);
  mcDEBUG_COMPONENT("source", mcposasource, mcrotasource)
  mccomp_posa[2] = mcposasource;
  mccomp_posr[2] = mcposrsource;
  mcNCounter[2]  = mcPCounter[2] = mcP2Counter[2] = 0;
  mcAbsorbProp[2]= 0;
    /* Component slit1. */
  /* Setting parameters for component slit1. */
  SIG_MESSAGE("slit1 (Init:SetPar)");
#line 96 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccslit1_xmin = -0.020;
#line 96 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccslit1_xmax = 0.065;
#line 97 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccslit1_ymin = -0.075;
#line 97 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccslit1_ymax = 0.075;
#line 46 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccslit1_radius = 0;
#line 46 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccslit1_xwidth = 0;
#line 46 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccslit1_yheight = 0;
#line 7211 "./linup-6.c"

  SIG_MESSAGE("slit1 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 98 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    (0)*DEG2RAD,
#line 98 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    (0)*DEG2RAD,
#line 98 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    (0)*DEG2RAD);
#line 7221 "./linup-6.c"
  rot_mul(mctr1, mcrotaa1, mcrotaslit1);
  rot_transpose(mcrotasource, mctr1);
  rot_mul(mcrotaslit1, mctr1, mcrotrslit1);
  mctc1 = coords_set(
#line 98 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    0,
#line 98 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    0,
#line 98 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    1.1215);
#line 7232 "./linup-6.c"
  rot_transpose(mcrotaa1, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaslit1 = coords_add(mcposaa1, mctc2);
  mctc1 = coords_sub(mcposasource, mcposaslit1);
  mcposrslit1 = rot_apply(mcrotaslit1, mctc1);
  mcDEBUG_COMPONENT("slit1", mcposaslit1, mcrotaslit1)
  mccomp_posa[3] = mcposaslit1;
  mccomp_posr[3] = mcposrslit1;
  mcNCounter[3]  = mcPCounter[3] = mcP2Counter[3] = 0;
  mcAbsorbProp[3]= 0;
    /* Component slit2. */
  /* Setting parameters for component slit2. */
  SIG_MESSAGE("slit2 (Init:SetPar)");
#line 101 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccslit2_xmin = -0.020;
#line 101 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccslit2_xmax = 0.020;
#line 102 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccslit2_ymin = -0.040;
#line 102 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccslit2_ymax = 0.040;
#line 46 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccslit2_radius = 0;
#line 46 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccslit2_xwidth = 0;
#line 46 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccslit2_yheight = 0;
#line 7260 "./linup-6.c"

  SIG_MESSAGE("slit2 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 103 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    (0)*DEG2RAD,
#line 103 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    (0)*DEG2RAD,
#line 103 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    (0)*DEG2RAD);
#line 7270 "./linup-6.c"
  rot_mul(mctr1, mcrotaa1, mcrotaslit2);
  rot_transpose(mcrotaslit1, mctr1);
  rot_mul(mcrotaslit2, mctr1, mcrotrslit2);
  mctc1 = coords_set(
#line 103 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    0,
#line 103 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    0,
#line 103 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    1.900);
#line 7281 "./linup-6.c"
  rot_transpose(mcrotaa1, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaslit2 = coords_add(mcposaa1, mctc2);
  mctc1 = coords_sub(mcposaslit1, mcposaslit2);
  mcposrslit2 = rot_apply(mcrotaslit2, mctc1);
  mcDEBUG_COMPONENT("slit2", mcposaslit2, mcrotaslit2)
  mccomp_posa[4] = mcposaslit2;
  mccomp_posr[4] = mcposrslit2;
  mcNCounter[4]  = mcPCounter[4] = mcP2Counter[4] = 0;
  mcAbsorbProp[4]= 0;
    /* Component slit3. */
  /* Setting parameters for component slit3. */
  SIG_MESSAGE("slit3 (Init:SetPar)");
#line 106 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccslit3_xmin = -0.021;
#line 106 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccslit3_xmax = 0.021;
#line 107 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccslit3_ymin = -0.041;
#line 107 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccslit3_ymax = 0.041;
#line 46 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccslit3_radius = 0;
#line 46 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccslit3_xwidth = 0;
#line 46 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccslit3_yheight = 0;
#line 7309 "./linup-6.c"

  SIG_MESSAGE("slit3 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 108 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    (0)*DEG2RAD,
#line 108 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    (0)*DEG2RAD,
#line 108 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    (0)*DEG2RAD);
#line 7319 "./linup-6.c"
  rot_mul(mctr1, mcrotaa1, mcrotaslit3);
  rot_transpose(mcrotaslit2, mctr1);
  rot_mul(mcrotaslit3, mctr1, mcrotrslit3);
  mctc1 = coords_set(
#line 108 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    0,
#line 108 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    0,
#line 108 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    3.288);
#line 7330 "./linup-6.c"
  rot_transpose(mcrotaa1, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaslit3 = coords_add(mcposaa1, mctc2);
  mctc1 = coords_sub(mcposaslit2, mcposaslit3);
  mcposrslit3 = rot_apply(mcrotaslit3, mctc1);
  mcDEBUG_COMPONENT("slit3", mcposaslit3, mcrotaslit3)
  mccomp_posa[5] = mcposaslit3;
  mccomp_posr[5] = mcposrslit3;
  mcNCounter[5]  = mcPCounter[5] = mcP2Counter[5] = 0;
  mcAbsorbProp[5]= 0;
    /* Component focus_mono. */
  /* Setting parameters for component focus_mono. */
  SIG_MESSAGE("focus_mono (Init:SetPar)");

  SIG_MESSAGE("focus_mono (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 111 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    (0)*DEG2RAD,
#line 111 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    (mcipPHM)*DEG2RAD,
#line 111 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    (0)*DEG2RAD);
#line 7353 "./linup-6.c"
  rot_mul(mctr1, mcrotaa1, mcrotafocus_mono);
  rot_transpose(mcrotaslit3, mctr1);
  rot_mul(mcrotafocus_mono, mctr1, mcrotrfocus_mono);
  mctc1 = coords_set(
#line 111 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    0,
#line 111 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    0,
#line 111 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    3.56);
#line 7364 "./linup-6.c"
  rot_transpose(mcrotaa1, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposafocus_mono = coords_add(mcposaa1, mctc2);
  mctc1 = coords_sub(mcposaslit3, mcposafocus_mono);
  mcposrfocus_mono = rot_apply(mcrotafocus_mono, mctc1);
  mcDEBUG_COMPONENT("focus_mono", mcposafocus_mono, mcrotafocus_mono)
  mccomp_posa[6] = mcposafocus_mono;
  mccomp_posr[6] = mcposrfocus_mono;
  mcNCounter[6]  = mcPCounter[6] = mcP2Counter[6] = 0;
  mcAbsorbProp[6]= 0;
    /* Component m0. */
  /* Setting parameters for component m0. */
  SIG_MESSAGE("m0 (Init:SetPar)");
#line 114 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccm0_zmin = -0.0375;
#line 114 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccm0_zmax = 0.0375;
#line 114 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccm0_ymin = -0.006;
#line 114 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccm0_ymax = 0.006;
#line 65 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccm0_zwidth = 0;
#line 65 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccm0_yheight = 0;
#line 115 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccm0_mosaich = tas1_mono_mosaic;
#line 115 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccm0_mosaicv = tas1_mono_mosaic;
#line 116 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccm0_r0 = tas1_mono_r0;
#line 116 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccm0_Q = tas1_mono_q;
#line 66 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccm0_DM = 0;
#line 7400 "./linup-6.c"

  SIG_MESSAGE("m0 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 118 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    (0)*DEG2RAD,
#line 118 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    (0)*DEG2RAD,
#line 118 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    (mrot0)*DEG2RAD);
#line 7410 "./linup-6.c"
  rot_mul(mctr1, mcrotafocus_mono, mcrotam0);
  rot_transpose(mcrotafocus_mono, mctr1);
  rot_mul(mcrotam0, mctr1, mcrotrm0);
  mctc1 = coords_set(
#line 117 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    0,
#line 117 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    mpos0,
#line 117 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    0);
#line 7421 "./linup-6.c"
  rot_transpose(mcrotafocus_mono, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposam0 = coords_add(mcposafocus_mono, mctc2);
  mctc1 = coords_sub(mcposafocus_mono, mcposam0);
  mcposrm0 = rot_apply(mcrotam0, mctc1);
  mcDEBUG_COMPONENT("m0", mcposam0, mcrotam0)
  mccomp_posa[7] = mcposam0;
  mccomp_posr[7] = mcposrm0;
  mcNCounter[7]  = mcPCounter[7] = mcP2Counter[7] = 0;
  mcAbsorbProp[7]= 0;
    /* Component m1. */
  /* Setting parameters for component m1. */
  SIG_MESSAGE("m1 (Init:SetPar)");
#line 121 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccm1_zmin = -0.0375;
#line 121 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccm1_zmax = 0.0375;
#line 121 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccm1_ymin = -0.006;
#line 121 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccm1_ymax = 0.006;
#line 65 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccm1_zwidth = 0;
#line 65 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccm1_yheight = 0;
#line 122 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccm1_mosaich = tas1_mono_mosaic;
#line 122 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccm1_mosaicv = tas1_mono_mosaic;
#line 123 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccm1_r0 = tas1_mono_r0;
#line 123 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccm1_Q = tas1_mono_q;
#line 66 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccm1_DM = 0;
#line 7457 "./linup-6.c"

  SIG_MESSAGE("m1 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 125 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    (0)*DEG2RAD,
#line 125 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    (0)*DEG2RAD,
#line 125 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    (mrot1)*DEG2RAD);
#line 7467 "./linup-6.c"
  rot_mul(mctr1, mcrotafocus_mono, mcrotam1);
  rot_transpose(mcrotam0, mctr1);
  rot_mul(mcrotam1, mctr1, mcrotrm1);
  mctc1 = coords_set(
#line 124 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    0,
#line 124 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    mpos1,
#line 124 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    0);
#line 7478 "./linup-6.c"
  rot_transpose(mcrotafocus_mono, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposam1 = coords_add(mcposafocus_mono, mctc2);
  mctc1 = coords_sub(mcposam0, mcposam1);
  mcposrm1 = rot_apply(mcrotam1, mctc1);
  mcDEBUG_COMPONENT("m1", mcposam1, mcrotam1)
  mccomp_posa[8] = mcposam1;
  mccomp_posr[8] = mcposrm1;
  mcNCounter[8]  = mcPCounter[8] = mcP2Counter[8] = 0;
  mcAbsorbProp[8]= 0;
    /* Component m2. */
  /* Setting parameters for component m2. */
  SIG_MESSAGE("m2 (Init:SetPar)");
#line 128 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccm2_zmin = -0.0375;
#line 128 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccm2_zmax = 0.0375;
#line 128 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccm2_ymin = -0.006;
#line 128 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccm2_ymax = 0.006;
#line 65 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccm2_zwidth = 0;
#line 65 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccm2_yheight = 0;
#line 129 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccm2_mosaich = tas1_mono_mosaic;
#line 129 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccm2_mosaicv = tas1_mono_mosaic;
#line 130 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccm2_r0 = tas1_mono_r0;
#line 130 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccm2_Q = tas1_mono_q;
#line 66 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccm2_DM = 0;
#line 7514 "./linup-6.c"

  SIG_MESSAGE("m2 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 132 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    (0)*DEG2RAD,
#line 132 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    (0)*DEG2RAD,
#line 132 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    (mrot2)*DEG2RAD);
#line 7524 "./linup-6.c"
  rot_mul(mctr1, mcrotafocus_mono, mcrotam2);
  rot_transpose(mcrotam1, mctr1);
  rot_mul(mcrotam2, mctr1, mcrotrm2);
  mctc1 = coords_set(
#line 131 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    0,
#line 131 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    mpos2,
#line 131 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    0);
#line 7535 "./linup-6.c"
  rot_transpose(mcrotafocus_mono, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposam2 = coords_add(mcposafocus_mono, mctc2);
  mctc1 = coords_sub(mcposam1, mcposam2);
  mcposrm2 = rot_apply(mcrotam2, mctc1);
  mcDEBUG_COMPONENT("m2", mcposam2, mcrotam2)
  mccomp_posa[9] = mcposam2;
  mccomp_posr[9] = mcposrm2;
  mcNCounter[9]  = mcPCounter[9] = mcP2Counter[9] = 0;
  mcAbsorbProp[9]= 0;
    /* Component m3. */
  /* Setting parameters for component m3. */
  SIG_MESSAGE("m3 (Init:SetPar)");
#line 135 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccm3_zmin = -0.0375;
#line 135 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccm3_zmax = 0.0375;
#line 135 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccm3_ymin = -0.006;
#line 135 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccm3_ymax = 0.006;
#line 65 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccm3_zwidth = 0;
#line 65 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccm3_yheight = 0;
#line 136 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccm3_mosaich = tas1_mono_mosaic;
#line 136 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccm3_mosaicv = tas1_mono_mosaic;
#line 137 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccm3_r0 = tas1_mono_r0;
#line 137 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccm3_Q = tas1_mono_q;
#line 66 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccm3_DM = 0;
#line 7571 "./linup-6.c"

  SIG_MESSAGE("m3 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 139 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    (0)*DEG2RAD,
#line 139 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    (0)*DEG2RAD,
#line 139 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    (mrot3)*DEG2RAD);
#line 7581 "./linup-6.c"
  rot_mul(mctr1, mcrotafocus_mono, mcrotam3);
  rot_transpose(mcrotam2, mctr1);
  rot_mul(mcrotam3, mctr1, mcrotrm3);
  mctc1 = coords_set(
#line 138 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    0,
#line 138 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    mpos3,
#line 138 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    0);
#line 7592 "./linup-6.c"
  rot_transpose(mcrotafocus_mono, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposam3 = coords_add(mcposafocus_mono, mctc2);
  mctc1 = coords_sub(mcposam2, mcposam3);
  mcposrm3 = rot_apply(mcrotam3, mctc1);
  mcDEBUG_COMPONENT("m3", mcposam3, mcrotam3)
  mccomp_posa[10] = mcposam3;
  mccomp_posr[10] = mcposrm3;
  mcNCounter[10]  = mcPCounter[10] = mcP2Counter[10] = 0;
  mcAbsorbProp[10]= 0;
    /* Component m4. */
  /* Setting parameters for component m4. */
  SIG_MESSAGE("m4 (Init:SetPar)");
#line 142 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccm4_zmin = -0.0375;
#line 142 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccm4_zmax = 0.0375;
#line 142 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccm4_ymin = -0.006;
#line 142 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccm4_ymax = 0.006;
#line 65 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccm4_zwidth = 0;
#line 65 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccm4_yheight = 0;
#line 143 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccm4_mosaich = tas1_mono_mosaic;
#line 143 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccm4_mosaicv = tas1_mono_mosaic;
#line 144 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccm4_r0 = tas1_mono_r0;
#line 144 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccm4_Q = tas1_mono_q;
#line 66 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccm4_DM = 0;
#line 7628 "./linup-6.c"

  SIG_MESSAGE("m4 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 146 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    (0)*DEG2RAD,
#line 146 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    (0)*DEG2RAD,
#line 146 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    (mrot4)*DEG2RAD);
#line 7638 "./linup-6.c"
  rot_mul(mctr1, mcrotafocus_mono, mcrotam4);
  rot_transpose(mcrotam3, mctr1);
  rot_mul(mcrotam4, mctr1, mcrotrm4);
  mctc1 = coords_set(
#line 145 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    0,
#line 145 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    mpos4,
#line 145 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    0);
#line 7649 "./linup-6.c"
  rot_transpose(mcrotafocus_mono, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposam4 = coords_add(mcposafocus_mono, mctc2);
  mctc1 = coords_sub(mcposam3, mcposam4);
  mcposrm4 = rot_apply(mcrotam4, mctc1);
  mcDEBUG_COMPONENT("m4", mcposam4, mcrotam4)
  mccomp_posa[11] = mcposam4;
  mccomp_posr[11] = mcposrm4;
  mcNCounter[11]  = mcPCounter[11] = mcP2Counter[11] = 0;
  mcAbsorbProp[11]= 0;
    /* Component m5. */
  /* Setting parameters for component m5. */
  SIG_MESSAGE("m5 (Init:SetPar)");
#line 149 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccm5_zmin = -0.0375;
#line 149 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccm5_zmax = 0.0375;
#line 149 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccm5_ymin = -0.006;
#line 149 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccm5_ymax = 0.006;
#line 65 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccm5_zwidth = 0;
#line 65 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccm5_yheight = 0;
#line 150 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccm5_mosaich = tas1_mono_mosaic;
#line 150 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccm5_mosaicv = tas1_mono_mosaic;
#line 151 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccm5_r0 = tas1_mono_r0;
#line 151 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccm5_Q = tas1_mono_q;
#line 66 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccm5_DM = 0;
#line 7685 "./linup-6.c"

  SIG_MESSAGE("m5 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 153 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    (0)*DEG2RAD,
#line 153 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    (0)*DEG2RAD,
#line 153 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    (mrot5)*DEG2RAD);
#line 7695 "./linup-6.c"
  rot_mul(mctr1, mcrotafocus_mono, mcrotam5);
  rot_transpose(mcrotam4, mctr1);
  rot_mul(mcrotam5, mctr1, mcrotrm5);
  mctc1 = coords_set(
#line 152 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    0,
#line 152 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    mpos5,
#line 152 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    0);
#line 7706 "./linup-6.c"
  rot_transpose(mcrotafocus_mono, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposam5 = coords_add(mcposafocus_mono, mctc2);
  mctc1 = coords_sub(mcposam4, mcposam5);
  mcposrm5 = rot_apply(mcrotam5, mctc1);
  mcDEBUG_COMPONENT("m5", mcposam5, mcrotam5)
  mccomp_posa[12] = mcposam5;
  mccomp_posr[12] = mcposrm5;
  mcNCounter[12]  = mcPCounter[12] = mcP2Counter[12] = 0;
  mcAbsorbProp[12]= 0;
    /* Component m6. */
  /* Setting parameters for component m6. */
  SIG_MESSAGE("m6 (Init:SetPar)");
#line 156 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccm6_zmin = -0.0375;
#line 156 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccm6_zmax = 0.0375;
#line 156 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccm6_ymin = -0.006;
#line 156 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccm6_ymax = 0.006;
#line 65 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccm6_zwidth = 0;
#line 65 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccm6_yheight = 0;
#line 157 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccm6_mosaich = tas1_mono_mosaic;
#line 157 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccm6_mosaicv = tas1_mono_mosaic;
#line 158 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccm6_r0 = tas1_mono_r0;
#line 158 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccm6_Q = tas1_mono_q;
#line 66 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccm6_DM = 0;
#line 7742 "./linup-6.c"

  SIG_MESSAGE("m6 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 160 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    (0)*DEG2RAD,
#line 160 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    (0)*DEG2RAD,
#line 160 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    (mrot6)*DEG2RAD);
#line 7752 "./linup-6.c"
  rot_mul(mctr1, mcrotafocus_mono, mcrotam6);
  rot_transpose(mcrotam5, mctr1);
  rot_mul(mcrotam6, mctr1, mcrotrm6);
  mctc1 = coords_set(
#line 159 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    0,
#line 159 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    mpos6,
#line 159 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    0);
#line 7763 "./linup-6.c"
  rot_transpose(mcrotafocus_mono, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposam6 = coords_add(mcposafocus_mono, mctc2);
  mctc1 = coords_sub(mcposam5, mcposam6);
  mcposrm6 = rot_apply(mcrotam6, mctc1);
  mcDEBUG_COMPONENT("m6", mcposam6, mcrotam6)
  mccomp_posa[13] = mcposam6;
  mccomp_posr[13] = mcposrm6;
  mcNCounter[13]  = mcPCounter[13] = mcP2Counter[13] = 0;
  mcAbsorbProp[13]= 0;
    /* Component m7. */
  /* Setting parameters for component m7. */
  SIG_MESSAGE("m7 (Init:SetPar)");
#line 163 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccm7_zmin = -0.0375;
#line 163 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccm7_zmax = 0.0375;
#line 163 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccm7_ymin = -0.006;
#line 163 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccm7_ymax = 0.006;
#line 65 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccm7_zwidth = 0;
#line 65 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccm7_yheight = 0;
#line 164 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccm7_mosaich = tas1_mono_mosaic;
#line 164 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccm7_mosaicv = tas1_mono_mosaic;
#line 165 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccm7_r0 = tas1_mono_r0;
#line 165 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccm7_Q = tas1_mono_q;
#line 66 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccm7_DM = 0;
#line 7799 "./linup-6.c"

  SIG_MESSAGE("m7 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 167 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    (0)*DEG2RAD,
#line 167 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    (0)*DEG2RAD,
#line 167 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    (mrot7)*DEG2RAD);
#line 7809 "./linup-6.c"
  rot_mul(mctr1, mcrotafocus_mono, mcrotam7);
  rot_transpose(mcrotam6, mctr1);
  rot_mul(mcrotam7, mctr1, mcrotrm7);
  mctc1 = coords_set(
#line 166 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    0,
#line 166 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    mpos7,
#line 166 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    0);
#line 7820 "./linup-6.c"
  rot_transpose(mcrotafocus_mono, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposam7 = coords_add(mcposafocus_mono, mctc2);
  mctc1 = coords_sub(mcposam6, mcposam7);
  mcposrm7 = rot_apply(mcrotam7, mctc1);
  mcDEBUG_COMPONENT("m7", mcposam7, mcrotam7)
  mccomp_posa[14] = mcposam7;
  mccomp_posr[14] = mcposrm7;
  mcNCounter[14]  = mcPCounter[14] = mcP2Counter[14] = 0;
  mcAbsorbProp[14]= 0;
    /* Component a2. */
  /* Setting parameters for component a2. */
  SIG_MESSAGE("a2 (Init:SetPar)");

  SIG_MESSAGE("a2 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 170 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    (0)*DEG2RAD,
#line 170 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    (mcipTTM)*DEG2RAD,
#line 170 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    (0)*DEG2RAD);
#line 7843 "./linup-6.c"
  rot_mul(mctr1, mcrotaa1, mcrotaa2);
  rot_transpose(mcrotam7, mctr1);
  rot_mul(mcrotaa2, mctr1, mcrotra2);
  mctc1 = coords_set(
#line 170 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    0,
#line 170 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    0,
#line 170 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    0);
#line 7854 "./linup-6.c"
  rot_transpose(mcrotafocus_mono, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaa2 = coords_add(mcposafocus_mono, mctc2);
  mctc1 = coords_sub(mcposam7, mcposaa2);
  mcposra2 = rot_apply(mcrotaa2, mctc1);
  mcDEBUG_COMPONENT("a2", mcposaa2, mcrotaa2)
  mccomp_posa[15] = mcposaa2;
  mccomp_posr[15] = mcposra2;
  mcNCounter[15]  = mcPCounter[15] = mcP2Counter[15] = 0;
  mcAbsorbProp[15]= 0;
    /* Component slitMS1. */
  /* Setting parameters for component slitMS1. */
  SIG_MESSAGE("slitMS1 (Init:SetPar)");
#line 173 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccslitMS1_xmin = -0.0105;
#line 173 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccslitMS1_xmax = 0.0105;
#line 173 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccslitMS1_ymin = -0.035;
#line 173 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccslitMS1_ymax = 0.035;
#line 46 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccslitMS1_radius = 0;
#line 46 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccslitMS1_xwidth = 0;
#line 46 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccslitMS1_yheight = 0;
#line 7882 "./linup-6.c"

  SIG_MESSAGE("slitMS1 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 174 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    (0)*DEG2RAD,
#line 174 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    (0)*DEG2RAD,
#line 174 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    (0)*DEG2RAD);
#line 7892 "./linup-6.c"
  rot_mul(mctr1, mcrotaa2, mcrotaslitMS1);
  rot_transpose(mcrotaa2, mctr1);
  rot_mul(mcrotaslitMS1, mctr1, mcrotrslitMS1);
  mctc1 = coords_set(
#line 174 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    0,
#line 174 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    0,
#line 174 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    0.565);
#line 7903 "./linup-6.c"
  rot_transpose(mcrotaa2, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaslitMS1 = coords_add(mcposaa2, mctc2);
  mctc1 = coords_sub(mcposaa2, mcposaslitMS1);
  mcposrslitMS1 = rot_apply(mcrotaslitMS1, mctc1);
  mcDEBUG_COMPONENT("slitMS1", mcposaslitMS1, mcrotaslitMS1)
  mccomp_posa[16] = mcposaslitMS1;
  mccomp_posr[16] = mcposrslitMS1;
  mcNCounter[16]  = mcPCounter[16] = mcP2Counter[16] = 0;
  mcAbsorbProp[16]= 0;
    /* Component slitMS2. */
  /* Setting parameters for component slitMS2. */
  SIG_MESSAGE("slitMS2 (Init:SetPar)");
#line 177 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccslitMS2_xmin = -0.0105;
#line 177 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccslitMS2_xmax = 0.0105;
#line 177 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccslitMS2_ymin = -0.035;
#line 177 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccslitMS2_ymax = 0.035;
#line 46 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccslitMS2_radius = 0;
#line 46 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccslitMS2_xwidth = 0;
#line 46 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccslitMS2_yheight = 0;
#line 7931 "./linup-6.c"

  SIG_MESSAGE("slitMS2 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 178 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    (0)*DEG2RAD,
#line 178 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    (0)*DEG2RAD,
#line 178 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    (0)*DEG2RAD);
#line 7941 "./linup-6.c"
  rot_mul(mctr1, mcrotaa2, mcrotaslitMS2);
  rot_transpose(mcrotaslitMS1, mctr1);
  rot_mul(mcrotaslitMS2, mctr1, mcrotrslitMS2);
  mctc1 = coords_set(
#line 178 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    0,
#line 178 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    0,
#line 178 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    0.855);
#line 7952 "./linup-6.c"
  rot_transpose(mcrotaa2, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaslitMS2 = coords_add(mcposaa2, mctc2);
  mctc1 = coords_sub(mcposaslitMS1, mcposaslitMS2);
  mcposrslitMS2 = rot_apply(mcrotaslitMS2, mctc1);
  mcDEBUG_COMPONENT("slitMS2", mcposaslitMS2, mcrotaslitMS2)
  mccomp_posa[17] = mcposaslitMS2;
  mccomp_posr[17] = mcposrslitMS2;
  mcNCounter[17]  = mcPCounter[17] = mcP2Counter[17] = 0;
  mcAbsorbProp[17]= 0;
    /* Component c1. */
  /* Setting parameters for component c1. */
  SIG_MESSAGE("c1 (Init:SetPar)");
#line 181 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccc1_xmin = -0.02;
#line 181 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccc1_xmax = 0.02;
#line 181 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccc1_ymin = -0.0375;
#line 181 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccc1_ymax = 0.0375;
#line 46 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccc1_xwidth = 0;
#line 46 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccc1_yheight = 0;
#line 182 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccc1_length = 0.250;
#line 182 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccc1_divergence = mcipC1;
#line 46 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccc1_transmission = 1;
#line 46 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccc1_divergenceV = 0;
#line 7986 "./linup-6.c"

  SIG_MESSAGE("c1 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 183 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    (0)*DEG2RAD,
#line 183 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    (OMC1_d)*DEG2RAD,
#line 183 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    (0)*DEG2RAD);
#line 7996 "./linup-6.c"
  rot_mul(mctr1, mcrotaa2, mcrotac1);
  rot_transpose(mcrotaslitMS2, mctr1);
  rot_mul(mcrotac1, mctr1, mcrotrc1);
  mctc1 = coords_set(
#line 183 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    0,
#line 183 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    0,
#line 183 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    0.87);
#line 8007 "./linup-6.c"
  rot_transpose(mcrotaa2, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposac1 = coords_add(mcposaa2, mctc2);
  mctc1 = coords_sub(mcposaslitMS2, mcposac1);
  mcposrc1 = rot_apply(mcrotac1, mctc1);
  mcDEBUG_COMPONENT("c1", mcposac1, mcrotac1)
  mccomp_posa[18] = mcposac1;
  mccomp_posr[18] = mcposrc1;
  mcNCounter[18]  = mcPCounter[18] = mcP2Counter[18] = 0;
  mcAbsorbProp[18]= 0;
    /* Component slitMS3. */
  /* Setting parameters for component slitMS3. */
  SIG_MESSAGE("slitMS3 (Init:SetPar)");
#line 46 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccslitMS3_xmin = 0;
#line 46 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccslitMS3_xmax = 0;
#line 46 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccslitMS3_ymin = 0;
#line 46 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccslitMS3_ymax = 0;
#line 185 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccslitMS3_radius = 0.025;
#line 46 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccslitMS3_xwidth = 0;
#line 46 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccslitMS3_yheight = 0;
#line 8035 "./linup-6.c"

  SIG_MESSAGE("slitMS3 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 186 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    (0)*DEG2RAD,
#line 186 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    (0)*DEG2RAD,
#line 186 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    (0)*DEG2RAD);
#line 8045 "./linup-6.c"
  rot_mul(mctr1, mcrotaa2, mcrotaslitMS3);
  rot_transpose(mcrotac1, mctr1);
  rot_mul(mcrotaslitMS3, mctr1, mcrotrslitMS3);
  mctc1 = coords_set(
#line 186 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    0,
#line 186 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    0,
#line 186 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    1.130);
#line 8056 "./linup-6.c"
  rot_transpose(mcrotaa2, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaslitMS3 = coords_add(mcposaa2, mctc2);
  mctc1 = coords_sub(mcposac1, mcposaslitMS3);
  mcposrslitMS3 = rot_apply(mcrotaslitMS3, mctc1);
  mcDEBUG_COMPONENT("slitMS3", mcposaslitMS3, mcrotaslitMS3)
  mccomp_posa[19] = mcposaslitMS3;
  mccomp_posr[19] = mcposrslitMS3;
  mcNCounter[19]  = mcPCounter[19] = mcP2Counter[19] = 0;
  mcAbsorbProp[19]= 0;
    /* Component slitMS4. */
  /* Setting parameters for component slitMS4. */
  SIG_MESSAGE("slitMS4 (Init:SetPar)");
#line 46 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccslitMS4_xmin = 0;
#line 46 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccslitMS4_xmax = 0;
#line 46 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccslitMS4_ymin = 0;
#line 46 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccslitMS4_ymax = 0;
#line 188 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccslitMS4_radius = 0.025;
#line 46 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccslitMS4_xwidth = 0;
#line 46 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccslitMS4_yheight = 0;
#line 8084 "./linup-6.c"

  SIG_MESSAGE("slitMS4 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 189 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    (0)*DEG2RAD,
#line 189 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    (0)*DEG2RAD,
#line 189 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    (0)*DEG2RAD);
#line 8094 "./linup-6.c"
  rot_mul(mctr1, mcrotaa2, mcrotaslitMS4);
  rot_transpose(mcrotaslitMS3, mctr1);
  rot_mul(mcrotaslitMS4, mctr1, mcrotrslitMS4);
  mctc1 = coords_set(
#line 189 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    0,
#line 189 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    0,
#line 189 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    1.180);
#line 8105 "./linup-6.c"
  rot_transpose(mcrotaa2, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaslitMS4 = coords_add(mcposaa2, mctc2);
  mctc1 = coords_sub(mcposaslitMS3, mcposaslitMS4);
  mcposrslitMS4 = rot_apply(mcrotaslitMS4, mctc1);
  mcDEBUG_COMPONENT("slitMS4", mcposaslitMS4, mcrotaslitMS4)
  mccomp_posa[20] = mcposaslitMS4;
  mccomp_posr[20] = mcposrslitMS4;
  mcNCounter[20]  = mcPCounter[20] = mcP2Counter[20] = 0;
  mcAbsorbProp[20]= 0;
    /* Component slitMS5. */
  /* Setting parameters for component slitMS5. */
  SIG_MESSAGE("slitMS5 (Init:SetPar)");
#line 46 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccslitMS5_xmin = 0;
#line 46 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccslitMS5_xmax = 0;
#line 46 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccslitMS5_ymin = 0;
#line 46 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccslitMS5_ymax = 0;
#line 191 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccslitMS5_radius = 0.0275;
#line 46 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccslitMS5_xwidth = 0;
#line 46 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccslitMS5_yheight = 0;
#line 8133 "./linup-6.c"

  SIG_MESSAGE("slitMS5 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 192 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    (0)*DEG2RAD,
#line 192 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    (0)*DEG2RAD,
#line 192 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    (0)*DEG2RAD);
#line 8143 "./linup-6.c"
  rot_mul(mctr1, mcrotaa2, mcrotaslitMS5);
  rot_transpose(mcrotaslitMS4, mctr1);
  rot_mul(mcrotaslitMS5, mctr1, mcrotrslitMS5);
  mctc1 = coords_set(
#line 192 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    0,
#line 192 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    0,
#line 192 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    1.230);
#line 8154 "./linup-6.c"
  rot_transpose(mcrotaa2, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaslitMS5 = coords_add(mcposaa2, mctc2);
  mctc1 = coords_sub(mcposaslitMS4, mcposaslitMS5);
  mcposrslitMS5 = rot_apply(mcrotaslitMS5, mctc1);
  mcDEBUG_COMPONENT("slitMS5", mcposaslitMS5, mcrotaslitMS5)
  mccomp_posa[21] = mcposaslitMS5;
  mccomp_posr[21] = mcposrslitMS5;
  mcNCounter[21]  = mcPCounter[21] = mcP2Counter[21] = 0;
  mcAbsorbProp[21]= 0;
    /* Component mon. */
  /* Setting parameters for component mon. */
  SIG_MESSAGE("mon (Init:SetPar)");
#line 195 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccmon_xmin = -0.025;
#line 195 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccmon_xmax = 0.025;
#line 195 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccmon_ymin = -0.0375;
#line 195 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccmon_ymax = 0.0375;
#line 46 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccmon_xwidth = 0;
#line 46 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccmon_yheight = 0;
#line 46 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccmon_restore_neutron = 0;
#line 8182 "./linup-6.c"

  SIG_MESSAGE("mon (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 196 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    (0)*DEG2RAD,
#line 196 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    (0)*DEG2RAD,
#line 196 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    (0)*DEG2RAD);
#line 8192 "./linup-6.c"
  rot_mul(mctr1, mcrotaa2, mcrotamon);
  rot_transpose(mcrotaslitMS5, mctr1);
  rot_mul(mcrotamon, mctr1, mcrotrmon);
  mctc1 = coords_set(
#line 196 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    0,
#line 196 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    0,
#line 196 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    1.280);
#line 8203 "./linup-6.c"
  rot_transpose(mcrotaa2, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposamon = coords_add(mcposaa2, mctc2);
  mctc1 = coords_sub(mcposaslitMS5, mcposamon);
  mcposrmon = rot_apply(mcrotamon, mctc1);
  mcDEBUG_COMPONENT("mon", mcposamon, mcrotamon)
  mccomp_posa[22] = mcposamon;
  mccomp_posr[22] = mcposrmon;
  mcNCounter[22]  = mcPCounter[22] = mcP2Counter[22] = 0;
  mcAbsorbProp[22]= 0;
    /* Component emon1. */
  /* Setting parameters for component emon1. */
  SIG_MESSAGE("emon1 (Init:SetPar)");
#line 201 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  if("linup_6_1.vmon") strncpy(mccemon1_filename, "linup_6_1.vmon" ? "linup_6_1.vmon" : "", 16384); else mccemon1_filename[0]='\0';
#line 199 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccemon1_xmin = -0.01;
#line 199 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccemon1_xmax = 0.01;
#line 199 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccemon1_ymin = -0.1;
#line 199 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccemon1_ymax = 0.1;
#line 54 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccemon1_xwidth = 0;
#line 54 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccemon1_yheight = 0;
#line 200 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccemon1_Emin = 19.25;
#line 200 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccemon1_Emax = 20.75;
#line 54 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccemon1_restore_neutron = 0;
#line 54 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccemon1_nowritefile = 0;
#line 8239 "./linup-6.c"

  SIG_MESSAGE("emon1 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 202 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    (0)*DEG2RAD,
#line 202 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    (0)*DEG2RAD,
#line 202 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    (0)*DEG2RAD);
#line 8249 "./linup-6.c"
  rot_mul(mctr1, mcrotaa2, mcrotaemon1);
  rot_transpose(mcrotamon, mctr1);
  rot_mul(mcrotaemon1, mctr1, mcrotremon1);
  mctc1 = coords_set(
#line 202 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    0,
#line 202 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    0,
#line 202 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    1.5);
#line 8260 "./linup-6.c"
  rot_transpose(mcrotaa2, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaemon1 = coords_add(mcposaa2, mctc2);
  mctc1 = coords_sub(mcposamon, mcposaemon1);
  mcposremon1 = rot_apply(mcrotaemon1, mctc1);
  mcDEBUG_COMPONENT("emon1", mcposaemon1, mcrotaemon1)
  mccomp_posa[23] = mcposaemon1;
  mccomp_posr[23] = mcposremon1;
  mcNCounter[23]  = mcPCounter[23] = mcP2Counter[23] = 0;
  mcAbsorbProp[23]= 0;
    /* Component sample. */
  /* Setting parameters for component sample. */
  SIG_MESSAGE("sample (Init:SetPar)");
#line 205 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccsample_radius = 0.010;
#line 205 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccsample_thickness = 0.0025;
#line 91 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccsample_zdepth = 0;
#line 91 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccsample_Vc = 13.827;
#line 91 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccsample_sigma_abs = 5.08;
#line 91 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccsample_sigma_inc = 5.08;
#line 92 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccsample_radius_i = 0;
#line 92 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccsample_radius_o = 0;
#line 92 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccsample_h = 0;
#line 92 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccsample_focus_r = 0;
#line 208 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccsample_pack = 1;
#line 92 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccsample_frac = 1;
#line 92 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccsample_f_QE = 0;
#line 92 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccsample_gamma = 0;
#line 93 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccsample_target_x = 0;
#line 93 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccsample_target_y = 0;
#line 93 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccsample_target_z = 0;
#line 207 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccsample_focus_xw = 0.04;
#line 207 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccsample_focus_yh = 0.063;
#line 94 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccsample_focus_aw = 0;
#line 94 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccsample_focus_ah = 0;
#line 94 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccsample_xwidth = 0;
#line 206 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccsample_yheight = 0.019;
#line 94 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccsample_zthick = 0;
#line 94 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccsample_rad_sphere = 0;
#line 94 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccsample_sig_a = 0;
#line 94 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccsample_sig_i = 0;
#line 94 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccsample_V0 = 0;
#line 209 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccsample_target_index = + 2;
#line 94 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccsample_multiples = 1;
#line 8334 "./linup-6.c"

  SIG_MESSAGE("sample (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 210 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    (0)*DEG2RAD,
#line 210 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    (0)*DEG2RAD,
#line 210 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    (0)*DEG2RAD);
#line 8344 "./linup-6.c"
  rot_mul(mctr1, mcrotaa2, mcrotasample);
  rot_transpose(mcrotaemon1, mctr1);
  rot_mul(mcrotasample, mctr1, mcrotrsample);
  mctc1 = coords_set(
#line 210 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    -0.0015,
#line 210 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    0,
#line 210 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    1.565);
#line 8355 "./linup-6.c"
  rot_transpose(mcrotaa2, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposasample = coords_add(mcposaa2, mctc2);
  mctc1 = coords_sub(mcposaemon1, mcposasample);
  mcposrsample = rot_apply(mcrotasample, mctc1);
  mcDEBUG_COMPONENT("sample", mcposasample, mcrotasample)
  mccomp_posa[24] = mcposasample;
  mccomp_posr[24] = mcposrsample;
  mcNCounter[24]  = mcPCounter[24] = mcP2Counter[24] = 0;
  mcAbsorbProp[24]= 0;
    /* Component a3. */
  /* Setting parameters for component a3. */
  SIG_MESSAGE("a3 (Init:SetPar)");

  SIG_MESSAGE("a3 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 213 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    (0)*DEG2RAD,
#line 213 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    (mcipTT)*DEG2RAD,
#line 213 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    (0)*DEG2RAD);
#line 8378 "./linup-6.c"
  rot_mul(mctr1, mcrotaa2, mcrotaa3);
  rot_transpose(mcrotasample, mctr1);
  rot_mul(mcrotaa3, mctr1, mcrotra3);
  mctc1 = coords_set(
#line 213 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    0,
#line 213 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    0,
#line 213 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    0);
#line 8389 "./linup-6.c"
  rot_transpose(mcrotasample, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaa3 = coords_add(mcposasample, mctc2);
  mctc1 = coords_sub(mcposasample, mcposaa3);
  mcposra3 = rot_apply(mcrotaa3, mctc1);
  mcDEBUG_COMPONENT("a3", mcposaa3, mcrotaa3)
  mccomp_posa[25] = mcposaa3;
  mccomp_posr[25] = mcposra3;
  mcNCounter[25]  = mcPCounter[25] = mcP2Counter[25] = 0;
  mcAbsorbProp[25]= 0;
    /* Component focus_check. */
  /* Setting parameters for component focus_check. */
  SIG_MESSAGE("focus_check (Init:SetPar)");
#line 217 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccfocus_check_nx = 20;
#line 217 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccfocus_check_ny = 20;
#line 218 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  if("linup_6.psd") strncpy(mccfocus_check_filename, "linup_6.psd" ? "linup_6.psd" : "", 16384); else mccfocus_check_filename[0]='\0';
#line 216 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccfocus_check_xmin = -0.02;
#line 216 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccfocus_check_xmax = 0.02;
#line 217 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccfocus_check_ymin = -0.0315;
#line 217 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccfocus_check_ymax = 0.0315;
#line 50 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccfocus_check_xwidth = 0;
#line 50 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccfocus_check_yheight = 0;
#line 51 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccfocus_check_restore_neutron = 0;
#line 8423 "./linup-6.c"

  SIG_MESSAGE("focus_check (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 219 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    (0)*DEG2RAD,
#line 219 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    (0)*DEG2RAD,
#line 219 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    (0)*DEG2RAD);
#line 8433 "./linup-6.c"
  rot_mul(mctr1, mcrotaa3, mcrotafocus_check);
  rot_transpose(mcrotaa3, mctr1);
  rot_mul(mcrotafocus_check, mctr1, mcrotrfocus_check);
  mctc1 = coords_set(
#line 219 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    0,
#line 219 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    0,
#line 219 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    0.369999);
#line 8444 "./linup-6.c"
  rot_transpose(mcrotaa3, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposafocus_check = coords_add(mcposaa3, mctc2);
  mctc1 = coords_sub(mcposaa3, mcposafocus_check);
  mcposrfocus_check = rot_apply(mcrotafocus_check, mctc1);
  mcDEBUG_COMPONENT("focus_check", mcposafocus_check, mcrotafocus_check)
  mccomp_posa[26] = mcposafocus_check;
  mccomp_posr[26] = mcposrfocus_check;
  mcNCounter[26]  = mcPCounter[26] = mcP2Counter[26] = 0;
  mcAbsorbProp[26]= 0;
    /* Component c2. */
  /* Setting parameters for component c2. */
  SIG_MESSAGE("c2 (Init:SetPar)");
#line 222 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccc2_xmin = -0.02;
#line 222 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccc2_xmax = 0.02;
#line 222 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccc2_ymin = -0.0315;
#line 222 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccc2_ymax = 0.0315;
#line 46 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccc2_xwidth = 0;
#line 46 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccc2_yheight = 0;
#line 223 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccc2_length = 0.300;
#line 223 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccc2_divergence = mcipC2;
#line 46 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccc2_transmission = 1;
#line 46 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccc2_divergenceV = 0;
#line 8478 "./linup-6.c"

  SIG_MESSAGE("c2 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 224 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    (0)*DEG2RAD,
#line 224 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    (0)*DEG2RAD,
#line 224 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    (0)*DEG2RAD);
#line 8488 "./linup-6.c"
  rot_mul(mctr1, mcrotaa3, mcrotac2);
  rot_transpose(mcrotafocus_check, mctr1);
  rot_mul(mcrotac2, mctr1, mcrotrc2);
  mctc1 = coords_set(
#line 224 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    0,
#line 224 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    0,
#line 224 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    0.370);
#line 8499 "./linup-6.c"
  rot_transpose(mcrotaa3, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposac2 = coords_add(mcposaa3, mctc2);
  mctc1 = coords_sub(mcposafocus_check, mcposac2);
  mcposrc2 = rot_apply(mcrotac2, mctc1);
  mcDEBUG_COMPONENT("c2", mcposac2, mcrotac2)
  mccomp_posa[27] = mcposac2;
  mccomp_posr[27] = mcposrc2;
  mcNCounter[27]  = mcPCounter[27] = mcP2Counter[27] = 0;
  mcAbsorbProp[27]= 0;
    /* Component ana. */
  /* Setting parameters for component ana. */
  SIG_MESSAGE("ana (Init:SetPar)");
#line 227 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccana_zmin = -0.0375;
#line 227 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccana_zmax = 0.0375;
#line 228 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccana_ymin = -0.024;
#line 228 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccana_ymax = 0.024;
#line 65 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccana_zwidth = 0;
#line 65 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccana_yheight = 0;
#line 229 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccana_mosaich = tas1_ana_mosaic;
#line 229 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccana_mosaicv = tas1_ana_mosaic;
#line 230 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccana_r0 = tas1_ana_r0;
#line 230 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccana_Q = tas1_ana_q;
#line 66 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccana_DM = 0;
#line 8535 "./linup-6.c"

  SIG_MESSAGE("ana (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 231 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    (0)*DEG2RAD,
#line 231 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    (mcipOMA)*DEG2RAD,
#line 231 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    (0)*DEG2RAD);
#line 8545 "./linup-6.c"
  rot_mul(mctr1, mcrotaa3, mcrotaana);
  rot_transpose(mcrotac2, mctr1);
  rot_mul(mcrotaana, mctr1, mcrotrana);
  mctc1 = coords_set(
#line 231 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    0,
#line 231 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    0,
#line 231 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    0.770);
#line 8556 "./linup-6.c"
  rot_transpose(mcrotaa3, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaana = coords_add(mcposaa3, mctc2);
  mctc1 = coords_sub(mcposac2, mcposaana);
  mcposrana = rot_apply(mcrotaana, mctc1);
  mcDEBUG_COMPONENT("ana", mcposaana, mcrotaana)
  mccomp_posa[28] = mcposaana;
  mccomp_posr[28] = mcposrana;
  mcNCounter[28]  = mcPCounter[28] = mcP2Counter[28] = 0;
  mcAbsorbProp[28]= 0;
    /* Component a4. */
  /* Setting parameters for component a4. */
  SIG_MESSAGE("a4 (Init:SetPar)");

  SIG_MESSAGE("a4 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 234 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    (0)*DEG2RAD,
#line 234 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    (mcipTTA)*DEG2RAD,
#line 234 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    (0)*DEG2RAD);
#line 8579 "./linup-6.c"
  rot_mul(mctr1, mcrotaa3, mcrotaa4);
  rot_transpose(mcrotaana, mctr1);
  rot_mul(mcrotaa4, mctr1, mcrotra4);
  mctc1 = coords_set(
#line 234 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    0,
#line 234 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    0,
#line 234 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    0);
#line 8590 "./linup-6.c"
  rot_transpose(mcrotaana, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaa4 = coords_add(mcposaana, mctc2);
  mctc1 = coords_sub(mcposaana, mcposaa4);
  mcposra4 = rot_apply(mcrotaa4, mctc1);
  mcDEBUG_COMPONENT("a4", mcposaa4, mcrotaa4)
  mccomp_posa[29] = mcposaa4;
  mccomp_posr[29] = mcposra4;
  mcNCounter[29]  = mcPCounter[29] = mcP2Counter[29] = 0;
  mcAbsorbProp[29]= 0;
    /* Component c3. */
  /* Setting parameters for component c3. */
  SIG_MESSAGE("c3 (Init:SetPar)");
#line 237 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccc3_xmin = -0.02;
#line 237 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccc3_xmax = 0.02;
#line 237 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccc3_ymin = -0.05;
#line 237 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccc3_ymax = 0.05;
#line 46 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccc3_xwidth = 0;
#line 46 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccc3_yheight = 0;
#line 238 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccc3_length = 0.270;
#line 238 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccc3_divergence = mcipC3;
#line 46 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccc3_transmission = 1;
#line 46 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccc3_divergenceV = 0;
#line 8624 "./linup-6.c"

  SIG_MESSAGE("c3 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 239 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    (0)*DEG2RAD,
#line 239 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    (0)*DEG2RAD,
#line 239 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    (0)*DEG2RAD);
#line 8634 "./linup-6.c"
  rot_mul(mctr1, mcrotaa4, mcrotac3);
  rot_transpose(mcrotaa4, mctr1);
  rot_mul(mcrotac3, mctr1, mcrotrc3);
  mctc1 = coords_set(
#line 239 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    0,
#line 239 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    0,
#line 239 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    0.104);
#line 8645 "./linup-6.c"
  rot_transpose(mcrotaa4, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposac3 = coords_add(mcposaa4, mctc2);
  mctc1 = coords_sub(mcposaa4, mcposac3);
  mcposrc3 = rot_apply(mcrotac3, mctc1);
  mcDEBUG_COMPONENT("c3", mcposac3, mcrotac3)
  mccomp_posa[30] = mcposac3;
  mccomp_posr[30] = mcposrc3;
  mcNCounter[30]  = mcPCounter[30] = mcP2Counter[30] = 0;
  mcAbsorbProp[30]= 0;
    /* Component sng. */
  /* Setting parameters for component sng. */
  SIG_MESSAGE("sng (Init:SetPar)");
#line 242 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccsng_xmin = -0.01;
#line 242 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccsng_xmax = 0.01;
#line 242 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccsng_ymin = -0.045;
#line 242 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccsng_ymax = 0.045;
#line 46 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccsng_xwidth = 0;
#line 46 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccsng_yheight = 0;
#line 46 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccsng_restore_neutron = 0;
#line 8673 "./linup-6.c"

  SIG_MESSAGE("sng (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 243 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    (0)*DEG2RAD,
#line 243 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    (0)*DEG2RAD,
#line 243 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    (0)*DEG2RAD);
#line 8683 "./linup-6.c"
  rot_mul(mctr1, mcrotaa4, mcrotasng);
  rot_transpose(mcrotac3, mctr1);
  rot_mul(mcrotasng, mctr1, mcrotrsng);
  mctc1 = coords_set(
#line 243 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    0,
#line 243 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    0,
#line 243 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    0.43);
#line 8694 "./linup-6.c"
  rot_transpose(mcrotaa4, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposasng = coords_add(mcposaa4, mctc2);
  mctc1 = coords_sub(mcposac3, mcposasng);
  mcposrsng = rot_apply(mcrotasng, mctc1);
  mcDEBUG_COMPONENT("sng", mcposasng, mcrotasng)
  mccomp_posa[31] = mcposasng;
  mccomp_posr[31] = mcposrsng;
  mcNCounter[31]  = mcPCounter[31] = mcP2Counter[31] = 0;
  mcAbsorbProp[31]= 0;
    /* Component emon2. */
  /* Setting parameters for component emon2. */
  SIG_MESSAGE("emon2 (Init:SetPar)");
#line 248 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  if("linup_6_2.vmon") strncpy(mccemon2_filename, "linup_6_2.vmon" ? "linup_6_2.vmon" : "", 16384); else mccemon2_filename[0]='\0';
#line 246 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccemon2_xmin = -0.0125;
#line 246 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccemon2_xmax = 0.0125;
#line 246 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccemon2_ymin = -0.05;
#line 246 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccemon2_ymax = 0.05;
#line 54 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccemon2_xwidth = 0;
#line 54 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccemon2_yheight = 0;
#line 247 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccemon2_Emin = 19.25;
#line 247 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccemon2_Emax = 20.75;
#line 54 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccemon2_restore_neutron = 0;
#line 54 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
  mccemon2_nowritefile = 0;
#line 8730 "./linup-6.c"

  SIG_MESSAGE("emon2 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 249 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    (0)*DEG2RAD,
#line 249 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    (0)*DEG2RAD,
#line 249 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    (0)*DEG2RAD);
#line 8740 "./linup-6.c"
  rot_mul(mctr1, mcrotaa4, mcrotaemon2);
  rot_transpose(mcrotasng, mctr1);
  rot_mul(mcrotaemon2, mctr1, mcrotremon2);
  mctc1 = coords_set(
#line 249 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    0,
#line 249 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    0,
#line 249 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/linup-6/linup-6.instr"
    0.430001);
#line 8751 "./linup-6.c"
  rot_transpose(mcrotaa4, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaemon2 = coords_add(mcposaa4, mctc2);
  mctc1 = coords_sub(mcposasng, mcposaemon2);
  mcposremon2 = rot_apply(mcrotaemon2, mctc1);
  mcDEBUG_COMPONENT("emon2", mcposaemon2, mcrotaemon2)
  mccomp_posa[32] = mcposaemon2;
  mccomp_posr[32] = mcposremon2;
  mcNCounter[32]  = mcPCounter[32] = mcP2Counter[32] = 0;
  mcAbsorbProp[32]= 0;
  /* Component initializations. */
  /* Initializations for component a1. */
  SIG_MESSAGE("a1 (Init)");
#define mccompcurname  a1
#define mccompcurtype  Progress_bar
#define mccompcurindex 1
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
#line 8788 "./linup-6.c"
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
#define mccompcurtype  Source_simple
#define mccompcurindex 2
#define pmul mccsource_pmul
#define square mccsource_square
#define srcArea mccsource_srcArea
#define radius mccsource_radius
#define yheight mccsource_yheight
#define xwidth mccsource_xwidth
#define dist mccsource_dist
#define focus_xw mccsource_focus_xw
#define focus_yh mccsource_focus_yh
#define E0 mccsource_E0
#define dE mccsource_dE
#define lambda0 mccsource_lambda0
#define dlambda mccsource_dlambda
#define flux mccsource_flux
#define gauss mccsource_gauss
#define target_index mccsource_target_index
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
#line 8882 "./linup-6.c"
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

  /* Initializations for component slit1. */
  SIG_MESSAGE("slit1 (Init)");
#define mccompcurname  slit1
#define mccompcurtype  Slit
#define mccompcurindex 3
#define xmin mccslit1_xmin
#define xmax mccslit1_xmax
#define ymin mccslit1_ymin
#define ymax mccslit1_ymax
#define radius mccslit1_radius
#define xwidth mccslit1_xwidth
#define yheight mccslit1_yheight
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
#line 8935 "./linup-6.c"
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

  /* Initializations for component slit2. */
  SIG_MESSAGE("slit2 (Init)");
#define mccompcurname  slit2
#define mccompcurtype  Slit
#define mccompcurindex 4
#define xmin mccslit2_xmin
#define xmax mccslit2_xmax
#define ymin mccslit2_ymin
#define ymax mccslit2_ymax
#define radius mccslit2_radius
#define xwidth mccslit2_xwidth
#define yheight mccslit2_yheight
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
#line 8979 "./linup-6.c"
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

  /* Initializations for component slit3. */
  SIG_MESSAGE("slit3 (Init)");
#define mccompcurname  slit3
#define mccompcurtype  Slit
#define mccompcurindex 5
#define xmin mccslit3_xmin
#define xmax mccslit3_xmax
#define ymin mccslit3_ymin
#define ymax mccslit3_ymax
#define radius mccslit3_radius
#define xwidth mccslit3_xwidth
#define yheight mccslit3_yheight
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
#line 9023 "./linup-6.c"
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

  /* Initializations for component focus_mono. */
  SIG_MESSAGE("focus_mono (Init)");

  /* Initializations for component m0. */
  SIG_MESSAGE("m0 (Init)");
#define mccompcurname  m0
#define mccompcurtype  Monochromator_flat
#define mccompcurindex 7
#define mos_rms_y mccm0_mos_rms_y
#define mos_rms_z mccm0_mos_rms_z
#define mos_rms_max mccm0_mos_rms_max
#define mono_Q mccm0_mono_Q
#define zmin mccm0_zmin
#define zmax mccm0_zmax
#define ymin mccm0_ymin
#define ymax mccm0_ymax
#define zwidth mccm0_zwidth
#define yheight mccm0_yheight
#define mosaich mccm0_mosaich
#define mosaicv mccm0_mosaicv
#define r0 mccm0_r0
#define Q mccm0_Q
#define DM mccm0_DM
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
#line 9073 "./linup-6.c"
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

  /* Initializations for component m1. */
  SIG_MESSAGE("m1 (Init)");
#define mccompcurname  m1
#define mccompcurtype  Monochromator_flat
#define mccompcurindex 8
#define mos_rms_y mccm1_mos_rms_y
#define mos_rms_z mccm1_mos_rms_z
#define mos_rms_max mccm1_mos_rms_max
#define mono_Q mccm1_mono_Q
#define zmin mccm1_zmin
#define zmax mccm1_zmax
#define ymin mccm1_ymin
#define ymax mccm1_ymax
#define zwidth mccm1_zwidth
#define yheight mccm1_yheight
#define mosaich mccm1_mosaich
#define mosaicv mccm1_mosaicv
#define r0 mccm1_r0
#define Q mccm1_Q
#define DM mccm1_DM
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
#line 9128 "./linup-6.c"
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

  /* Initializations for component m2. */
  SIG_MESSAGE("m2 (Init)");
#define mccompcurname  m2
#define mccompcurtype  Monochromator_flat
#define mccompcurindex 9
#define mos_rms_y mccm2_mos_rms_y
#define mos_rms_z mccm2_mos_rms_z
#define mos_rms_max mccm2_mos_rms_max
#define mono_Q mccm2_mono_Q
#define zmin mccm2_zmin
#define zmax mccm2_zmax
#define ymin mccm2_ymin
#define ymax mccm2_ymax
#define zwidth mccm2_zwidth
#define yheight mccm2_yheight
#define mosaich mccm2_mosaich
#define mosaicv mccm2_mosaicv
#define r0 mccm2_r0
#define Q mccm2_Q
#define DM mccm2_DM
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
#line 9183 "./linup-6.c"
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

  /* Initializations for component m3. */
  SIG_MESSAGE("m3 (Init)");
#define mccompcurname  m3
#define mccompcurtype  Monochromator_flat
#define mccompcurindex 10
#define mos_rms_y mccm3_mos_rms_y
#define mos_rms_z mccm3_mos_rms_z
#define mos_rms_max mccm3_mos_rms_max
#define mono_Q mccm3_mono_Q
#define zmin mccm3_zmin
#define zmax mccm3_zmax
#define ymin mccm3_ymin
#define ymax mccm3_ymax
#define zwidth mccm3_zwidth
#define yheight mccm3_yheight
#define mosaich mccm3_mosaich
#define mosaicv mccm3_mosaicv
#define r0 mccm3_r0
#define Q mccm3_Q
#define DM mccm3_DM
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
#line 9238 "./linup-6.c"
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

  /* Initializations for component m4. */
  SIG_MESSAGE("m4 (Init)");
#define mccompcurname  m4
#define mccompcurtype  Monochromator_flat
#define mccompcurindex 11
#define mos_rms_y mccm4_mos_rms_y
#define mos_rms_z mccm4_mos_rms_z
#define mos_rms_max mccm4_mos_rms_max
#define mono_Q mccm4_mono_Q
#define zmin mccm4_zmin
#define zmax mccm4_zmax
#define ymin mccm4_ymin
#define ymax mccm4_ymax
#define zwidth mccm4_zwidth
#define yheight mccm4_yheight
#define mosaich mccm4_mosaich
#define mosaicv mccm4_mosaicv
#define r0 mccm4_r0
#define Q mccm4_Q
#define DM mccm4_DM
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
#line 9293 "./linup-6.c"
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

  /* Initializations for component m5. */
  SIG_MESSAGE("m5 (Init)");
#define mccompcurname  m5
#define mccompcurtype  Monochromator_flat
#define mccompcurindex 12
#define mos_rms_y mccm5_mos_rms_y
#define mos_rms_z mccm5_mos_rms_z
#define mos_rms_max mccm5_mos_rms_max
#define mono_Q mccm5_mono_Q
#define zmin mccm5_zmin
#define zmax mccm5_zmax
#define ymin mccm5_ymin
#define ymax mccm5_ymax
#define zwidth mccm5_zwidth
#define yheight mccm5_yheight
#define mosaich mccm5_mosaich
#define mosaicv mccm5_mosaicv
#define r0 mccm5_r0
#define Q mccm5_Q
#define DM mccm5_DM
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
#line 9348 "./linup-6.c"
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

  /* Initializations for component m6. */
  SIG_MESSAGE("m6 (Init)");
#define mccompcurname  m6
#define mccompcurtype  Monochromator_flat
#define mccompcurindex 13
#define mos_rms_y mccm6_mos_rms_y
#define mos_rms_z mccm6_mos_rms_z
#define mos_rms_max mccm6_mos_rms_max
#define mono_Q mccm6_mono_Q
#define zmin mccm6_zmin
#define zmax mccm6_zmax
#define ymin mccm6_ymin
#define ymax mccm6_ymax
#define zwidth mccm6_zwidth
#define yheight mccm6_yheight
#define mosaich mccm6_mosaich
#define mosaicv mccm6_mosaicv
#define r0 mccm6_r0
#define Q mccm6_Q
#define DM mccm6_DM
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
#line 9403 "./linup-6.c"
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

  /* Initializations for component m7. */
  SIG_MESSAGE("m7 (Init)");
#define mccompcurname  m7
#define mccompcurtype  Monochromator_flat
#define mccompcurindex 14
#define mos_rms_y mccm7_mos_rms_y
#define mos_rms_z mccm7_mos_rms_z
#define mos_rms_max mccm7_mos_rms_max
#define mono_Q mccm7_mono_Q
#define zmin mccm7_zmin
#define zmax mccm7_zmax
#define ymin mccm7_ymin
#define ymax mccm7_ymax
#define zwidth mccm7_zwidth
#define yheight mccm7_yheight
#define mosaich mccm7_mosaich
#define mosaicv mccm7_mosaicv
#define r0 mccm7_r0
#define Q mccm7_Q
#define DM mccm7_DM
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
#line 9458 "./linup-6.c"
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

  /* Initializations for component a2. */
  SIG_MESSAGE("a2 (Init)");

  /* Initializations for component slitMS1. */
  SIG_MESSAGE("slitMS1 (Init)");
#define mccompcurname  slitMS1
#define mccompcurtype  Slit
#define mccompcurindex 16
#define xmin mccslitMS1_xmin
#define xmax mccslitMS1_xmax
#define ymin mccslitMS1_ymin
#define ymax mccslitMS1_ymax
#define radius mccslitMS1_radius
#define xwidth mccslitMS1_xwidth
#define yheight mccslitMS1_yheight
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
#line 9513 "./linup-6.c"
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

  /* Initializations for component slitMS2. */
  SIG_MESSAGE("slitMS2 (Init)");
#define mccompcurname  slitMS2
#define mccompcurtype  Slit
#define mccompcurindex 17
#define xmin mccslitMS2_xmin
#define xmax mccslitMS2_xmax
#define ymin mccslitMS2_ymin
#define ymax mccslitMS2_ymax
#define radius mccslitMS2_radius
#define xwidth mccslitMS2_xwidth
#define yheight mccslitMS2_yheight
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
#line 9557 "./linup-6.c"
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

  /* Initializations for component c1. */
  SIG_MESSAGE("c1 (Init)");
#define mccompcurname  c1
#define mccompcurtype  Collimator_linear
#define mccompcurindex 18
#define slope mccc1_slope
#define slopeV mccc1_slopeV
#define xmin mccc1_xmin
#define xmax mccc1_xmax
#define ymin mccc1_ymin
#define ymax mccc1_ymax
#define xwidth mccc1_xwidth
#define yheight mccc1_yheight
#define length mccc1_length
#define divergence mccc1_divergence
#define transmission mccc1_transmission
#define divergenceV mccc1_divergenceV
#line 54 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Collimator_linear.comp"
{
slope = tan(MIN2RAD*divergence);
  slopeV= tan(MIN2RAD*divergenceV);
  if (xwidth  > 0) { xmax = xwidth/2;  xmin = -xmax; }
  if (yheight > 0) { ymax = yheight/2; ymin = -ymax; }

  if ((xmin >= xmax) || (ymin >= ymax)) {
    printf("Collimator_linear: %s: Null slit opening area !\n"
	         "ERROR              (xwidth,yheight,xmin,xmax,ymin,ymax). Exiting",
           NAME_CURRENT_COMP);
    exit(0);
  }

}
#line 9601 "./linup-6.c"
#undef divergenceV
#undef transmission
#undef divergence
#undef length
#undef yheight
#undef xwidth
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef slopeV
#undef slope
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component slitMS3. */
  SIG_MESSAGE("slitMS3 (Init)");
#define mccompcurname  slitMS3
#define mccompcurtype  Slit
#define mccompcurindex 19
#define xmin mccslitMS3_xmin
#define xmax mccslitMS3_xmax
#define ymin mccslitMS3_ymin
#define ymax mccslitMS3_ymax
#define radius mccslitMS3_radius
#define xwidth mccslitMS3_xwidth
#define yheight mccslitMS3_yheight
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
#line 9650 "./linup-6.c"
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

  /* Initializations for component slitMS4. */
  SIG_MESSAGE("slitMS4 (Init)");
#define mccompcurname  slitMS4
#define mccompcurtype  Slit
#define mccompcurindex 20
#define xmin mccslitMS4_xmin
#define xmax mccslitMS4_xmax
#define ymin mccslitMS4_ymin
#define ymax mccslitMS4_ymax
#define radius mccslitMS4_radius
#define xwidth mccslitMS4_xwidth
#define yheight mccslitMS4_yheight
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
#line 9694 "./linup-6.c"
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

  /* Initializations for component slitMS5. */
  SIG_MESSAGE("slitMS5 (Init)");
#define mccompcurname  slitMS5
#define mccompcurtype  Slit
#define mccompcurindex 21
#define xmin mccslitMS5_xmin
#define xmax mccslitMS5_xmax
#define ymin mccslitMS5_ymin
#define ymax mccslitMS5_ymax
#define radius mccslitMS5_radius
#define xwidth mccslitMS5_xwidth
#define yheight mccslitMS5_yheight
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
#line 9738 "./linup-6.c"
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

  /* Initializations for component mon. */
  SIG_MESSAGE("mon (Init)");
#define mccompcurname  mon
#define mccompcurtype  Monitor
#define mccompcurindex 22
#define Nsum mccmon_Nsum
#define psum mccmon_psum
#define p2sum mccmon_p2sum
#define xmin mccmon_xmin
#define xmax mccmon_xmax
#define ymin mccmon_ymin
#define ymax mccmon_ymax
#define xwidth mccmon_xwidth
#define yheight mccmon_yheight
#define restore_neutron mccmon_restore_neutron
#line 57 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/Monitor.comp"
{
psum = 0;
p2sum = 0;
Nsum = 0;

if (xwidth  > 0) { xmax = xwidth/2;  xmin = -xmax; }
    if (yheight > 0) { ymax = yheight/2; ymin = -ymax; }

    if ((xmin >= xmax) || (ymin >= ymax)) {
            printf("Monitor: %s: Null detection area !\n"
                   "ERROR    (xwidth,yheight,xmin,xmax,ymin,ymax). Exiting",
           NAME_CURRENT_COMP);
      exit(0);
    }
}
#line 9781 "./linup-6.c"
#undef restore_neutron
#undef yheight
#undef xwidth
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef p2sum
#undef psum
#undef Nsum
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component emon1. */
  SIG_MESSAGE("emon1 (Init)");
#define mccompcurname  emon1
#define mccompcurtype  E_monitor
#define mccompcurindex 23
#define nE mccemon1_nE
#define E_N mccemon1_E_N
#define E_p mccemon1_E_p
#define E_p2 mccemon1_E_p2
#define S_p mccemon1_S_p
#define S_pE mccemon1_S_pE
#define S_pE2 mccemon1_S_pE2
#define filename mccemon1_filename
#define xmin mccemon1_xmin
#define xmax mccemon1_xmax
#define ymin mccemon1_ymin
#define ymax mccemon1_ymax
#define xwidth mccemon1_xwidth
#define yheight mccemon1_yheight
#define Emin mccemon1_Emin
#define Emax mccemon1_Emax
#define restore_neutron mccemon1_restore_neutron
#define nowritefile mccemon1_nowritefile
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
#line 9841 "./linup-6.c"
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

  /* Initializations for component sample. */
  SIG_MESSAGE("sample (Init)");
#define mccompcurname  sample
#define mccompcurtype  V_sample
#define mccompcurindex 24
#define VarsV mccsample_VarsV
#define radius mccsample_radius
#define thickness mccsample_thickness
#define zdepth mccsample_zdepth
#define Vc mccsample_Vc
#define sigma_abs mccsample_sigma_abs
#define sigma_inc mccsample_sigma_inc
#define radius_i mccsample_radius_i
#define radius_o mccsample_radius_o
#define h mccsample_h
#define focus_r mccsample_focus_r
#define pack mccsample_pack
#define frac mccsample_frac
#define f_QE mccsample_f_QE
#define gamma mccsample_gamma
#define target_x mccsample_target_x
#define target_y mccsample_target_y
#define target_z mccsample_target_z
#define focus_xw mccsample_focus_xw
#define focus_yh mccsample_focus_yh
#define focus_aw mccsample_focus_aw
#define focus_ah mccsample_focus_ah
#define xwidth mccsample_xwidth
#define yheight mccsample_yheight
#define zthick mccsample_zthick
#define rad_sphere mccsample_rad_sphere
#define sig_a mccsample_sig_a
#define sig_i mccsample_sig_i
#define V0 mccsample_V0
#define target_index mccsample_target_index
#define multiples mccsample_multiples
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
#line 9958 "./linup-6.c"
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

  /* Initializations for component a3. */
  SIG_MESSAGE("a3 (Init)");

  /* Initializations for component focus_check. */
  SIG_MESSAGE("focus_check (Init)");
#define mccompcurname  focus_check
#define mccompcurtype  PSD_monitor
#define mccompcurindex 26
#define PSD_N mccfocus_check_PSD_N
#define PSD_p mccfocus_check_PSD_p
#define PSD_p2 mccfocus_check_PSD_p2
#define nx mccfocus_check_nx
#define ny mccfocus_check_ny
#define filename mccfocus_check_filename
#define xmin mccfocus_check_xmin
#define xmax mccfocus_check_xmax
#define ymin mccfocus_check_ymin
#define ymax mccfocus_check_ymax
#define xwidth mccfocus_check_xwidth
#define yheight mccfocus_check_yheight
#define restore_neutron mccfocus_check_restore_neutron
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
#line 10040 "./linup-6.c"
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

  /* Initializations for component c2. */
  SIG_MESSAGE("c2 (Init)");
#define mccompcurname  c2
#define mccompcurtype  Collimator_linear
#define mccompcurindex 27
#define slope mccc2_slope
#define slopeV mccc2_slopeV
#define xmin mccc2_xmin
#define xmax mccc2_xmax
#define ymin mccc2_ymin
#define ymax mccc2_ymax
#define xwidth mccc2_xwidth
#define yheight mccc2_yheight
#define length mccc2_length
#define divergence mccc2_divergence
#define transmission mccc2_transmission
#define divergenceV mccc2_divergenceV
#line 54 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Collimator_linear.comp"
{
slope = tan(MIN2RAD*divergence);
  slopeV= tan(MIN2RAD*divergenceV);
  if (xwidth  > 0) { xmax = xwidth/2;  xmin = -xmax; }
  if (yheight > 0) { ymax = yheight/2; ymin = -ymax; }

  if ((xmin >= xmax) || (ymin >= ymax)) {
    printf("Collimator_linear: %s: Null slit opening area !\n"
	         "ERROR              (xwidth,yheight,xmin,xmax,ymin,ymax). Exiting",
           NAME_CURRENT_COMP);
    exit(0);
  }

}
#line 10090 "./linup-6.c"
#undef divergenceV
#undef transmission
#undef divergence
#undef length
#undef yheight
#undef xwidth
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef slopeV
#undef slope
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component ana. */
  SIG_MESSAGE("ana (Init)");
#define mccompcurname  ana
#define mccompcurtype  Monochromator_flat
#define mccompcurindex 28
#define mos_rms_y mccana_mos_rms_y
#define mos_rms_z mccana_mos_rms_z
#define mos_rms_max mccana_mos_rms_max
#define mono_Q mccana_mono_Q
#define zmin mccana_zmin
#define zmax mccana_zmax
#define ymin mccana_ymin
#define ymax mccana_ymax
#define zwidth mccana_zwidth
#define yheight mccana_yheight
#define mosaich mccana_mosaich
#define mosaicv mccana_mosaicv
#define r0 mccana_r0
#define Q mccana_Q
#define DM mccana_DM
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
#line 10142 "./linup-6.c"
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

  /* Initializations for component a4. */
  SIG_MESSAGE("a4 (Init)");

  /* Initializations for component c3. */
  SIG_MESSAGE("c3 (Init)");
#define mccompcurname  c3
#define mccompcurtype  Collimator_linear
#define mccompcurindex 30
#define slope mccc3_slope
#define slopeV mccc3_slopeV
#define xmin mccc3_xmin
#define xmax mccc3_xmax
#define ymin mccc3_ymin
#define ymax mccc3_ymax
#define xwidth mccc3_xwidth
#define yheight mccc3_yheight
#define length mccc3_length
#define divergence mccc3_divergence
#define transmission mccc3_transmission
#define divergenceV mccc3_divergenceV
#line 54 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Collimator_linear.comp"
{
slope = tan(MIN2RAD*divergence);
  slopeV= tan(MIN2RAD*divergenceV);
  if (xwidth  > 0) { xmax = xwidth/2;  xmin = -xmax; }
  if (yheight > 0) { ymax = yheight/2; ymin = -ymax; }

  if ((xmin >= xmax) || (ymin >= ymax)) {
    printf("Collimator_linear: %s: Null slit opening area !\n"
	         "ERROR              (xwidth,yheight,xmin,xmax,ymin,ymax). Exiting",
           NAME_CURRENT_COMP);
    exit(0);
  }

}
#line 10197 "./linup-6.c"
#undef divergenceV
#undef transmission
#undef divergence
#undef length
#undef yheight
#undef xwidth
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef slopeV
#undef slope
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component sng. */
  SIG_MESSAGE("sng (Init)");
#define mccompcurname  sng
#define mccompcurtype  Monitor
#define mccompcurindex 31
#define Nsum mccsng_Nsum
#define psum mccsng_psum
#define p2sum mccsng_p2sum
#define xmin mccsng_xmin
#define xmax mccsng_xmax
#define ymin mccsng_ymin
#define ymax mccsng_ymax
#define xwidth mccsng_xwidth
#define yheight mccsng_yheight
#define restore_neutron mccsng_restore_neutron
#line 57 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/Monitor.comp"
{
psum = 0;
p2sum = 0;
Nsum = 0;

if (xwidth  > 0) { xmax = xwidth/2;  xmin = -xmax; }
    if (yheight > 0) { ymax = yheight/2; ymin = -ymax; }

    if ((xmin >= xmax) || (ymin >= ymax)) {
            printf("Monitor: %s: Null detection area !\n"
                   "ERROR    (xwidth,yheight,xmin,xmax,ymin,ymax). Exiting",
           NAME_CURRENT_COMP);
      exit(0);
    }
}
#line 10245 "./linup-6.c"
#undef restore_neutron
#undef yheight
#undef xwidth
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef p2sum
#undef psum
#undef Nsum
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component emon2. */
  SIG_MESSAGE("emon2 (Init)");
#define mccompcurname  emon2
#define mccompcurtype  E_monitor
#define mccompcurindex 32
#define nE mccemon2_nE
#define E_N mccemon2_E_N
#define E_p mccemon2_E_p
#define E_p2 mccemon2_E_p2
#define S_p mccemon2_S_p
#define S_pE mccemon2_S_pE
#define S_pE2 mccemon2_S_pE2
#define filename mccemon2_filename
#define xmin mccemon2_xmin
#define xmax mccemon2_xmax
#define ymin mccemon2_ymin
#define ymax mccemon2_ymax
#define xwidth mccemon2_xwidth
#define yheight mccemon2_yheight
#define Emin mccemon2_Emin
#define Emax mccemon2_Emax
#define restore_neutron mccemon2_restore_neutron
#define nowritefile mccemon2_nowritefile
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
#line 10305 "./linup-6.c"
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
  /* SPLIT counter for component sample */
  int mcSplit_sample=0;
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
#define mccompcurtype  Progress_bar
#define mccompcurindex 1
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
#line 10483 "./linup-6.c"
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

  /* TRACE Component source [2] */
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
#define mccompcurname  source
#define mccompcurtype  Source_simple
#define mccompcurindex 2
#define pmul mccsource_pmul
#define square mccsource_square
#define srcArea mccsource_srcArea
{   /* Declarations of source=Source_simple() SETTING parameters. */
MCNUM radius = mccsource_radius;
MCNUM yheight = mccsource_yheight;
MCNUM xwidth = mccsource_xwidth;
MCNUM dist = mccsource_dist;
MCNUM focus_xw = mccsource_focus_xw;
MCNUM focus_yh = mccsource_focus_yh;
MCNUM E0 = mccsource_E0;
MCNUM dE = mccsource_dE;
MCNUM lambda0 = mccsource_lambda0;
MCNUM dlambda = mccsource_dlambda;
MCNUM flux = mccsource_flux;
MCNUM gauss = mccsource_gauss;
int target_index = mccsource_target_index;
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
#line 10654 "./linup-6.c"
}   /* End of source=Source_simple() SETTING parameter declarations. */
#undef srcArea
#undef square
#undef pmul
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompsource:
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

  /* TRACE Component slit1 [3] */
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
#define mccompcurname  slit1
#define mccompcurtype  Slit
#define mccompcurindex 3
{   /* Declarations of slit1=Slit() SETTING parameters. */
MCNUM xmin = mccslit1_xmin;
MCNUM xmax = mccslit1_xmax;
MCNUM ymin = mccslit1_ymin;
MCNUM ymax = mccslit1_ymax;
MCNUM radius = mccslit1_radius;
MCNUM xwidth = mccslit1_xwidth;
MCNUM yheight = mccslit1_yheight;
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
#line 10781 "./linup-6.c"
}   /* End of slit1=Slit() SETTING parameter declarations. */
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompslit1:
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

  /* TRACE Component slit2 [4] */
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
#define mccompcurname  slit2
#define mccompcurtype  Slit
#define mccompcurindex 4
{   /* Declarations of slit2=Slit() SETTING parameters. */
MCNUM xmin = mccslit2_xmin;
MCNUM xmax = mccslit2_xmax;
MCNUM ymin = mccslit2_ymin;
MCNUM ymax = mccslit2_ymax;
MCNUM radius = mccslit2_radius;
MCNUM xwidth = mccslit2_xwidth;
MCNUM yheight = mccslit2_yheight;
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
#line 10905 "./linup-6.c"
}   /* End of slit2=Slit() SETTING parameter declarations. */
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompslit2:
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

  /* TRACE Component slit3 [5] */
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
#define mccompcurname  slit3
#define mccompcurtype  Slit
#define mccompcurindex 5
{   /* Declarations of slit3=Slit() SETTING parameters. */
MCNUM xmin = mccslit3_xmin;
MCNUM xmax = mccslit3_xmax;
MCNUM ymin = mccslit3_ymin;
MCNUM ymax = mccslit3_ymax;
MCNUM radius = mccslit3_radius;
MCNUM xwidth = mccslit3_xwidth;
MCNUM yheight = mccslit3_yheight;
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
#line 11029 "./linup-6.c"
}   /* End of slit3=Slit() SETTING parameter declarations. */
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompslit3:
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

  /* TRACE Component focus_mono [6] */
  mccoordschange(mcposrfocus_mono, mcrotrfocus_mono,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component focus_mono (without coords transformations) */
  mcJumpTrace_focus_mono:
  SIG_MESSAGE("focus_mono (Trace)");
  mcDEBUG_COMP("focus_mono")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbCompfocus_mono
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
#define mccompcurname  focus_mono
#define mccompcurtype  Arm
#define mccompcurindex 6
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompfocus_mono:
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

  /* TRACE Component m0 [7] */
  mccoordschange(mcposrm0, mcrotrm0,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component m0 (without coords transformations) */
  mcJumpTrace_m0:
  SIG_MESSAGE("m0 (Trace)");
  mcDEBUG_COMP("m0")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbCompm0
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
#define mccompcurname  m0
#define mccompcurtype  Monochromator_flat
#define mccompcurindex 7
#define mos_rms_y mccm0_mos_rms_y
#define mos_rms_z mccm0_mos_rms_z
#define mos_rms_max mccm0_mos_rms_max
#define mono_Q mccm0_mono_Q
{   /* Declarations of m0=Monochromator_flat() SETTING parameters. */
MCNUM zmin = mccm0_zmin;
MCNUM zmax = mccm0_zmax;
MCNUM ymin = mccm0_ymin;
MCNUM ymax = mccm0_ymax;
MCNUM zwidth = mccm0_zwidth;
MCNUM yheight = mccm0_yheight;
MCNUM mosaich = mccm0_mosaich;
MCNUM mosaicv = mccm0_mosaicv;
MCNUM r0 = mccm0_r0;
MCNUM Q = mccm0_Q;
MCNUM DM = mccm0_DM;
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
#line 11388 "./linup-6.c"
}   /* End of m0=Monochromator_flat() SETTING parameter declarations. */
#undef mono_Q
#undef mos_rms_max
#undef mos_rms_z
#undef mos_rms_y
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompm0:
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

  /* TRACE Component m1 [8] */
  mccoordschange(mcposrm1, mcrotrm1,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component m1 (without coords transformations) */
  mcJumpTrace_m1:
  SIG_MESSAGE("m1 (Trace)");
  mcDEBUG_COMP("m1")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbCompm1
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
#define mccompcurname  m1
#define mccompcurtype  Monochromator_flat
#define mccompcurindex 8
#define mos_rms_y mccm1_mos_rms_y
#define mos_rms_z mccm1_mos_rms_z
#define mos_rms_max mccm1_mos_rms_max
#define mono_Q mccm1_mono_Q
{   /* Declarations of m1=Monochromator_flat() SETTING parameters. */
MCNUM zmin = mccm1_zmin;
MCNUM zmax = mccm1_zmax;
MCNUM ymin = mccm1_ymin;
MCNUM ymax = mccm1_ymax;
MCNUM zwidth = mccm1_zwidth;
MCNUM yheight = mccm1_yheight;
MCNUM mosaich = mccm1_mosaich;
MCNUM mosaicv = mccm1_mosaicv;
MCNUM r0 = mccm1_r0;
MCNUM Q = mccm1_Q;
MCNUM DM = mccm1_DM;
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
#line 11648 "./linup-6.c"
}   /* End of m1=Monochromator_flat() SETTING parameter declarations. */
#undef mono_Q
#undef mos_rms_max
#undef mos_rms_z
#undef mos_rms_y
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompm1:
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

  /* TRACE Component m2 [9] */
  mccoordschange(mcposrm2, mcrotrm2,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component m2 (without coords transformations) */
  mcJumpTrace_m2:
  SIG_MESSAGE("m2 (Trace)");
  mcDEBUG_COMP("m2")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbCompm2
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
#define mccompcurname  m2
#define mccompcurtype  Monochromator_flat
#define mccompcurindex 9
#define mos_rms_y mccm2_mos_rms_y
#define mos_rms_z mccm2_mos_rms_z
#define mos_rms_max mccm2_mos_rms_max
#define mono_Q mccm2_mono_Q
{   /* Declarations of m2=Monochromator_flat() SETTING parameters. */
MCNUM zmin = mccm2_zmin;
MCNUM zmax = mccm2_zmax;
MCNUM ymin = mccm2_ymin;
MCNUM ymax = mccm2_ymax;
MCNUM zwidth = mccm2_zwidth;
MCNUM yheight = mccm2_yheight;
MCNUM mosaich = mccm2_mosaich;
MCNUM mosaicv = mccm2_mosaicv;
MCNUM r0 = mccm2_r0;
MCNUM Q = mccm2_Q;
MCNUM DM = mccm2_DM;
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
#line 11908 "./linup-6.c"
}   /* End of m2=Monochromator_flat() SETTING parameter declarations. */
#undef mono_Q
#undef mos_rms_max
#undef mos_rms_z
#undef mos_rms_y
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompm2:
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

  /* TRACE Component m3 [10] */
  mccoordschange(mcposrm3, mcrotrm3,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component m3 (without coords transformations) */
  mcJumpTrace_m3:
  SIG_MESSAGE("m3 (Trace)");
  mcDEBUG_COMP("m3")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbCompm3
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
#define mccompcurname  m3
#define mccompcurtype  Monochromator_flat
#define mccompcurindex 10
#define mos_rms_y mccm3_mos_rms_y
#define mos_rms_z mccm3_mos_rms_z
#define mos_rms_max mccm3_mos_rms_max
#define mono_Q mccm3_mono_Q
{   /* Declarations of m3=Monochromator_flat() SETTING parameters. */
MCNUM zmin = mccm3_zmin;
MCNUM zmax = mccm3_zmax;
MCNUM ymin = mccm3_ymin;
MCNUM ymax = mccm3_ymax;
MCNUM zwidth = mccm3_zwidth;
MCNUM yheight = mccm3_yheight;
MCNUM mosaich = mccm3_mosaich;
MCNUM mosaicv = mccm3_mosaicv;
MCNUM r0 = mccm3_r0;
MCNUM Q = mccm3_Q;
MCNUM DM = mccm3_DM;
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
#line 12168 "./linup-6.c"
}   /* End of m3=Monochromator_flat() SETTING parameter declarations. */
#undef mono_Q
#undef mos_rms_max
#undef mos_rms_z
#undef mos_rms_y
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompm3:
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

  /* TRACE Component m4 [11] */
  mccoordschange(mcposrm4, mcrotrm4,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component m4 (without coords transformations) */
  mcJumpTrace_m4:
  SIG_MESSAGE("m4 (Trace)");
  mcDEBUG_COMP("m4")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbCompm4
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
#define mccompcurname  m4
#define mccompcurtype  Monochromator_flat
#define mccompcurindex 11
#define mos_rms_y mccm4_mos_rms_y
#define mos_rms_z mccm4_mos_rms_z
#define mos_rms_max mccm4_mos_rms_max
#define mono_Q mccm4_mono_Q
{   /* Declarations of m4=Monochromator_flat() SETTING parameters. */
MCNUM zmin = mccm4_zmin;
MCNUM zmax = mccm4_zmax;
MCNUM ymin = mccm4_ymin;
MCNUM ymax = mccm4_ymax;
MCNUM zwidth = mccm4_zwidth;
MCNUM yheight = mccm4_yheight;
MCNUM mosaich = mccm4_mosaich;
MCNUM mosaicv = mccm4_mosaicv;
MCNUM r0 = mccm4_r0;
MCNUM Q = mccm4_Q;
MCNUM DM = mccm4_DM;
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
#line 12428 "./linup-6.c"
}   /* End of m4=Monochromator_flat() SETTING parameter declarations. */
#undef mono_Q
#undef mos_rms_max
#undef mos_rms_z
#undef mos_rms_y
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompm4:
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

  /* TRACE Component m5 [12] */
  mccoordschange(mcposrm5, mcrotrm5,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component m5 (without coords transformations) */
  mcJumpTrace_m5:
  SIG_MESSAGE("m5 (Trace)");
  mcDEBUG_COMP("m5")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbCompm5
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
#define mccompcurname  m5
#define mccompcurtype  Monochromator_flat
#define mccompcurindex 12
#define mos_rms_y mccm5_mos_rms_y
#define mos_rms_z mccm5_mos_rms_z
#define mos_rms_max mccm5_mos_rms_max
#define mono_Q mccm5_mono_Q
{   /* Declarations of m5=Monochromator_flat() SETTING parameters. */
MCNUM zmin = mccm5_zmin;
MCNUM zmax = mccm5_zmax;
MCNUM ymin = mccm5_ymin;
MCNUM ymax = mccm5_ymax;
MCNUM zwidth = mccm5_zwidth;
MCNUM yheight = mccm5_yheight;
MCNUM mosaich = mccm5_mosaich;
MCNUM mosaicv = mccm5_mosaicv;
MCNUM r0 = mccm5_r0;
MCNUM Q = mccm5_Q;
MCNUM DM = mccm5_DM;
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
#line 12688 "./linup-6.c"
}   /* End of m5=Monochromator_flat() SETTING parameter declarations. */
#undef mono_Q
#undef mos_rms_max
#undef mos_rms_z
#undef mos_rms_y
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompm5:
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

  /* TRACE Component m6 [13] */
  mccoordschange(mcposrm6, mcrotrm6,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component m6 (without coords transformations) */
  mcJumpTrace_m6:
  SIG_MESSAGE("m6 (Trace)");
  mcDEBUG_COMP("m6")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbCompm6
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
#define mccompcurname  m6
#define mccompcurtype  Monochromator_flat
#define mccompcurindex 13
#define mos_rms_y mccm6_mos_rms_y
#define mos_rms_z mccm6_mos_rms_z
#define mos_rms_max mccm6_mos_rms_max
#define mono_Q mccm6_mono_Q
{   /* Declarations of m6=Monochromator_flat() SETTING parameters. */
MCNUM zmin = mccm6_zmin;
MCNUM zmax = mccm6_zmax;
MCNUM ymin = mccm6_ymin;
MCNUM ymax = mccm6_ymax;
MCNUM zwidth = mccm6_zwidth;
MCNUM yheight = mccm6_yheight;
MCNUM mosaich = mccm6_mosaich;
MCNUM mosaicv = mccm6_mosaicv;
MCNUM r0 = mccm6_r0;
MCNUM Q = mccm6_Q;
MCNUM DM = mccm6_DM;
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
#line 12948 "./linup-6.c"
}   /* End of m6=Monochromator_flat() SETTING parameter declarations. */
#undef mono_Q
#undef mos_rms_max
#undef mos_rms_z
#undef mos_rms_y
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompm6:
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

  /* TRACE Component m7 [14] */
  mccoordschange(mcposrm7, mcrotrm7,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component m7 (without coords transformations) */
  mcJumpTrace_m7:
  SIG_MESSAGE("m7 (Trace)");
  mcDEBUG_COMP("m7")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbCompm7
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
#define mccompcurname  m7
#define mccompcurtype  Monochromator_flat
#define mccompcurindex 14
#define mos_rms_y mccm7_mos_rms_y
#define mos_rms_z mccm7_mos_rms_z
#define mos_rms_max mccm7_mos_rms_max
#define mono_Q mccm7_mono_Q
{   /* Declarations of m7=Monochromator_flat() SETTING parameters. */
MCNUM zmin = mccm7_zmin;
MCNUM zmax = mccm7_zmax;
MCNUM ymin = mccm7_ymin;
MCNUM ymax = mccm7_ymax;
MCNUM zwidth = mccm7_zwidth;
MCNUM yheight = mccm7_yheight;
MCNUM mosaich = mccm7_mosaich;
MCNUM mosaicv = mccm7_mosaicv;
MCNUM r0 = mccm7_r0;
MCNUM Q = mccm7_Q;
MCNUM DM = mccm7_DM;
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
#line 13208 "./linup-6.c"
}   /* End of m7=Monochromator_flat() SETTING parameter declarations. */
#undef mono_Q
#undef mos_rms_max
#undef mos_rms_z
#undef mos_rms_y
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompm7:
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

  /* TRACE Component a2 [15] */
  mccoordschange(mcposra2, mcrotra2,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component a2 (without coords transformations) */
  mcJumpTrace_a2:
  SIG_MESSAGE("a2 (Trace)");
  mcDEBUG_COMP("a2")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbCompa2
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
#define mccompcurname  a2
#define mccompcurtype  Arm
#define mccompcurindex 15
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompa2:
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

  /* TRACE Component slitMS1 [16] */
  mccoordschange(mcposrslitMS1, mcrotrslitMS1,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component slitMS1 (without coords transformations) */
  mcJumpTrace_slitMS1:
  SIG_MESSAGE("slitMS1 (Trace)");
  mcDEBUG_COMP("slitMS1")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbCompslitMS1
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
#define mccompcurname  slitMS1
#define mccompcurtype  Slit
#define mccompcurindex 16
{   /* Declarations of slitMS1=Slit() SETTING parameters. */
MCNUM xmin = mccslitMS1_xmin;
MCNUM xmax = mccslitMS1_xmax;
MCNUM ymin = mccslitMS1_ymin;
MCNUM ymax = mccslitMS1_ymax;
MCNUM radius = mccslitMS1_radius;
MCNUM xwidth = mccslitMS1_xwidth;
MCNUM yheight = mccslitMS1_yheight;
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
#line 13439 "./linup-6.c"
}   /* End of slitMS1=Slit() SETTING parameter declarations. */
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompslitMS1:
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

  /* TRACE Component slitMS2 [17] */
  mccoordschange(mcposrslitMS2, mcrotrslitMS2,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component slitMS2 (without coords transformations) */
  mcJumpTrace_slitMS2:
  SIG_MESSAGE("slitMS2 (Trace)");
  mcDEBUG_COMP("slitMS2")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbCompslitMS2
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
#define mccompcurname  slitMS2
#define mccompcurtype  Slit
#define mccompcurindex 17
{   /* Declarations of slitMS2=Slit() SETTING parameters. */
MCNUM xmin = mccslitMS2_xmin;
MCNUM xmax = mccslitMS2_xmax;
MCNUM ymin = mccslitMS2_ymin;
MCNUM ymax = mccslitMS2_ymax;
MCNUM radius = mccslitMS2_radius;
MCNUM xwidth = mccslitMS2_xwidth;
MCNUM yheight = mccslitMS2_yheight;
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
#line 13563 "./linup-6.c"
}   /* End of slitMS2=Slit() SETTING parameter declarations. */
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompslitMS2:
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

  /* TRACE Component c1 [18] */
  mccoordschange(mcposrc1, mcrotrc1,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component c1 (without coords transformations) */
  mcJumpTrace_c1:
  SIG_MESSAGE("c1 (Trace)");
  mcDEBUG_COMP("c1")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbCompc1
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
#define mccompcurname  c1
#define mccompcurtype  Collimator_linear
#define mccompcurindex 18
#define slope mccc1_slope
#define slopeV mccc1_slopeV
{   /* Declarations of c1=Collimator_linear() SETTING parameters. */
MCNUM xmin = mccc1_xmin;
MCNUM xmax = mccc1_xmax;
MCNUM ymin = mccc1_ymin;
MCNUM ymax = mccc1_ymax;
MCNUM xwidth = mccc1_xwidth;
MCNUM yheight = mccc1_yheight;
MCNUM length = mccc1_length;
MCNUM divergence = mccc1_divergence;
MCNUM transmission = mccc1_transmission;
MCNUM divergenceV = mccc1_divergenceV;
#line 69 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Collimator_linear.comp"
{
    double phi, dt;

    PROP_Z0;
    if (x<xmin || x>xmax || y<ymin || y>ymax)
      ABSORB;
    dt = length/vz;
    PROP_DT(dt);
    if (x<xmin || x>xmax || y<ymin || y>ymax)
      ABSORB;

    if(slope > 0.0)
    {
      phi = fabs(vx/vz);
      if (phi > slope)
        ABSORB;
      else
        p *= transmission*(1.0 - phi/slope);
      SCATTER;
    }
    if (slopeV > 0) {
      phi = fabs(vy/vz);
      if (phi > slopeV)
        ABSORB;
      else
        p *= transmission*(1.0 - phi/slopeV);
      SCATTER;
    }
}
#line 13711 "./linup-6.c"
}   /* End of c1=Collimator_linear() SETTING parameter declarations. */
#undef slopeV
#undef slope
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompc1:
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

  /* TRACE Component slitMS3 [19] */
  mccoordschange(mcposrslitMS3, mcrotrslitMS3,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component slitMS3 (without coords transformations) */
  mcJumpTrace_slitMS3:
  SIG_MESSAGE("slitMS3 (Trace)");
  mcDEBUG_COMP("slitMS3")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbCompslitMS3
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
#define mccompcurname  slitMS3
#define mccompcurtype  Slit
#define mccompcurindex 19
{   /* Declarations of slitMS3=Slit() SETTING parameters. */
MCNUM xmin = mccslitMS3_xmin;
MCNUM xmax = mccslitMS3_xmax;
MCNUM ymin = mccslitMS3_ymin;
MCNUM ymax = mccslitMS3_ymax;
MCNUM radius = mccslitMS3_radius;
MCNUM xwidth = mccslitMS3_xwidth;
MCNUM yheight = mccslitMS3_yheight;
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
#line 13837 "./linup-6.c"
}   /* End of slitMS3=Slit() SETTING parameter declarations. */
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompslitMS3:
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

  /* TRACE Component slitMS4 [20] */
  mccoordschange(mcposrslitMS4, mcrotrslitMS4,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component slitMS4 (without coords transformations) */
  mcJumpTrace_slitMS4:
  SIG_MESSAGE("slitMS4 (Trace)");
  mcDEBUG_COMP("slitMS4")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbCompslitMS4
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
#define mccompcurname  slitMS4
#define mccompcurtype  Slit
#define mccompcurindex 20
{   /* Declarations of slitMS4=Slit() SETTING parameters. */
MCNUM xmin = mccslitMS4_xmin;
MCNUM xmax = mccslitMS4_xmax;
MCNUM ymin = mccslitMS4_ymin;
MCNUM ymax = mccslitMS4_ymax;
MCNUM radius = mccslitMS4_radius;
MCNUM xwidth = mccslitMS4_xwidth;
MCNUM yheight = mccslitMS4_yheight;
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
#line 13961 "./linup-6.c"
}   /* End of slitMS4=Slit() SETTING parameter declarations. */
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompslitMS4:
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

  /* TRACE Component slitMS5 [21] */
  mccoordschange(mcposrslitMS5, mcrotrslitMS5,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component slitMS5 (without coords transformations) */
  mcJumpTrace_slitMS5:
  SIG_MESSAGE("slitMS5 (Trace)");
  mcDEBUG_COMP("slitMS5")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbCompslitMS5
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
#define mccompcurname  slitMS5
#define mccompcurtype  Slit
#define mccompcurindex 21
{   /* Declarations of slitMS5=Slit() SETTING parameters. */
MCNUM xmin = mccslitMS5_xmin;
MCNUM xmax = mccslitMS5_xmax;
MCNUM ymin = mccslitMS5_ymin;
MCNUM ymax = mccslitMS5_ymax;
MCNUM radius = mccslitMS5_radius;
MCNUM xwidth = mccslitMS5_xwidth;
MCNUM yheight = mccslitMS5_yheight;
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
#line 14085 "./linup-6.c"
}   /* End of slitMS5=Slit() SETTING parameter declarations. */
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompslitMS5:
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

  /* TRACE Component mon [22] */
  mccoordschange(mcposrmon, mcrotrmon,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component mon (without coords transformations) */
  mcJumpTrace_mon:
  SIG_MESSAGE("mon (Trace)");
  mcDEBUG_COMP("mon")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbCompmon
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
#define mccompcurname  mon
#define mccompcurtype  Monitor
#define mccompcurindex 22
#define Nsum mccmon_Nsum
#define psum mccmon_psum
#define p2sum mccmon_p2sum
{   /* Declarations of mon=Monitor() SETTING parameters. */
MCNUM xmin = mccmon_xmin;
MCNUM xmax = mccmon_xmax;
MCNUM ymin = mccmon_ymin;
MCNUM ymax = mccmon_ymax;
MCNUM xwidth = mccmon_xwidth;
MCNUM yheight = mccmon_yheight;
MCNUM restore_neutron = mccmon_restore_neutron;
#line 74 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/Monitor.comp"
{
    PROP_Z0;
    if (x>xmin && x<xmax && y>ymin && y<ymax)
    {
      Nsum++;
      psum += p;
      p2sum += p*p;
      SCATTER;
    }
    if (restore_neutron) {
      RESTORE_NEUTRON(INDEX_CURRENT_COMP, x, y, z, vx, vy, vz, t, sx, sy, sz, p);
    }
}
#line 14215 "./linup-6.c"
}   /* End of mon=Monitor() SETTING parameter declarations. */
#undef p2sum
#undef psum
#undef Nsum
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompmon:
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

  /* TRACE Component emon1 [23] */
  mccoordschange(mcposremon1, mcrotremon1,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component emon1 (without coords transformations) */
  mcJumpTrace_emon1:
  SIG_MESSAGE("emon1 (Trace)");
  mcDEBUG_COMP("emon1")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbCompemon1
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
#define mccompcurname  emon1
#define mccompcurtype  E_monitor
#define mccompcurindex 23
#define nE mccemon1_nE
#define E_N mccemon1_E_N
#define E_p mccemon1_E_p
#define E_p2 mccemon1_E_p2
#define S_p mccemon1_S_p
#define S_pE mccemon1_S_pE
#define S_pE2 mccemon1_S_pE2
{   /* Declarations of emon1=E_monitor() SETTING parameters. */
char* filename = mccemon1_filename;
MCNUM xmin = mccemon1_xmin;
MCNUM xmax = mccemon1_xmax;
MCNUM ymin = mccemon1_ymin;
MCNUM ymax = mccemon1_ymax;
MCNUM xwidth = mccemon1_xwidth;
MCNUM yheight = mccemon1_yheight;
MCNUM Emin = mccemon1_Emin;
MCNUM Emax = mccemon1_Emax;
MCNUM restore_neutron = mccemon1_restore_neutron;
int nowritefile = mccemon1_nowritefile;
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
#line 14369 "./linup-6.c"
}   /* End of emon1=E_monitor() SETTING parameter declarations. */
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
  mcabsorbCompemon1:
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

  /* TRACE Component sample [24] */
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
  if (!mcSplit_sample) {                   /* STORE only the first time */
    if (floor(10) > 1) p /= floor(10); /* adapt weight for SPLITed neutron */
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
  } else {
    RESTORE_NEUTRON(24,
      mcnlx,
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
  }
  mcSplit_sample++; /* SPLIT number */
  mcScattered=0;
  mcRestore=0;
  mcNCounter[24]++;
  mcPCounter[24] += p;
  mcP2Counter[24] += p*p;
#define mccompcurname  sample
#define mccompcurtype  V_sample
#define mccompcurindex 24
#define VarsV mccsample_VarsV
{   /* Declarations of sample=V_sample() SETTING parameters. */
MCNUM radius = mccsample_radius;
MCNUM thickness = mccsample_thickness;
MCNUM zdepth = mccsample_zdepth;
MCNUM Vc = mccsample_Vc;
MCNUM sigma_abs = mccsample_sigma_abs;
MCNUM sigma_inc = mccsample_sigma_inc;
MCNUM radius_i = mccsample_radius_i;
MCNUM radius_o = mccsample_radius_o;
MCNUM h = mccsample_h;
MCNUM focus_r = mccsample_focus_r;
MCNUM pack = mccsample_pack;
MCNUM frac = mccsample_frac;
MCNUM f_QE = mccsample_f_QE;
MCNUM gamma = mccsample_gamma;
MCNUM target_x = mccsample_target_x;
MCNUM target_y = mccsample_target_y;
MCNUM target_z = mccsample_target_z;
MCNUM focus_xw = mccsample_focus_xw;
MCNUM focus_yh = mccsample_focus_yh;
MCNUM focus_aw = mccsample_focus_aw;
MCNUM focus_ah = mccsample_focus_ah;
MCNUM xwidth = mccsample_xwidth;
MCNUM yheight = mccsample_yheight;
MCNUM zthick = mccsample_zthick;
MCNUM rad_sphere = mccsample_rad_sphere;
MCNUM sig_a = mccsample_sig_a;
MCNUM sig_i = mccsample_sig_i;
MCNUM V0 = mccsample_V0;
int target_index = mccsample_target_index;
MCNUM multiples = mccsample_multiples;
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
#line 14669 "./linup-6.c"
}   /* End of sample=V_sample() SETTING parameter declarations. */
#undef VarsV
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompsample:
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

  /* TRACE Component a3 [25] */
  mccoordschange(mcposra3, mcrotra3,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component a3 (without coords transformations) */
  mcJumpTrace_a3:
  SIG_MESSAGE("a3 (Trace)");
  mcDEBUG_COMP("a3")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbCompa3
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
#define mccompcurname  a3
#define mccompcurtype  Arm
#define mccompcurindex 25
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompa3:
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

  /* TRACE Component focus_check [26] */
  mccoordschange(mcposrfocus_check, mcrotrfocus_check,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component focus_check (without coords transformations) */
  mcJumpTrace_focus_check:
  SIG_MESSAGE("focus_check (Trace)");
  mcDEBUG_COMP("focus_check")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbCompfocus_check
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
#define mccompcurname  focus_check
#define mccompcurtype  PSD_monitor
#define mccompcurindex 26
#define PSD_N mccfocus_check_PSD_N
#define PSD_p mccfocus_check_PSD_p
#define PSD_p2 mccfocus_check_PSD_p2
{   /* Declarations of focus_check=PSD_monitor() SETTING parameters. */
int nx = mccfocus_check_nx;
int ny = mccfocus_check_ny;
char* filename = mccfocus_check_filename;
MCNUM xmin = mccfocus_check_xmin;
MCNUM xmax = mccfocus_check_xmax;
MCNUM ymin = mccfocus_check_ymin;
MCNUM ymax = mccfocus_check_ymax;
MCNUM xwidth = mccfocus_check_xwidth;
MCNUM yheight = mccfocus_check_yheight;
MCNUM restore_neutron = mccfocus_check_restore_neutron;
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
#line 14907 "./linup-6.c"
}   /* End of focus_check=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompfocus_check:
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

  /* TRACE Component c2 [27] */
  mccoordschange(mcposrc2, mcrotrc2,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component c2 (without coords transformations) */
  mcJumpTrace_c2:
  SIG_MESSAGE("c2 (Trace)");
  mcDEBUG_COMP("c2")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbCompc2
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
#define mccompcurname  c2
#define mccompcurtype  Collimator_linear
#define mccompcurindex 27
#define slope mccc2_slope
#define slopeV mccc2_slopeV
{   /* Declarations of c2=Collimator_linear() SETTING parameters. */
MCNUM xmin = mccc2_xmin;
MCNUM xmax = mccc2_xmax;
MCNUM ymin = mccc2_ymin;
MCNUM ymax = mccc2_ymax;
MCNUM xwidth = mccc2_xwidth;
MCNUM yheight = mccc2_yheight;
MCNUM length = mccc2_length;
MCNUM divergence = mccc2_divergence;
MCNUM transmission = mccc2_transmission;
MCNUM divergenceV = mccc2_divergenceV;
#line 69 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Collimator_linear.comp"
{
    double phi, dt;

    PROP_Z0;
    if (x<xmin || x>xmax || y<ymin || y>ymax)
      ABSORB;
    dt = length/vz;
    PROP_DT(dt);
    if (x<xmin || x>xmax || y<ymin || y>ymax)
      ABSORB;

    if(slope > 0.0)
    {
      phi = fabs(vx/vz);
      if (phi > slope)
        ABSORB;
      else
        p *= transmission*(1.0 - phi/slope);
      SCATTER;
    }
    if (slopeV > 0) {
      phi = fabs(vy/vz);
      if (phi > slopeV)
        ABSORB;
      else
        p *= transmission*(1.0 - phi/slopeV);
      SCATTER;
    }
}
#line 15058 "./linup-6.c"
}   /* End of c2=Collimator_linear() SETTING parameter declarations. */
#undef slopeV
#undef slope
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompc2:
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

  /* TRACE Component ana [28] */
  mccoordschange(mcposrana, mcrotrana,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component ana (without coords transformations) */
  mcJumpTrace_ana:
  SIG_MESSAGE("ana (Trace)");
  mcDEBUG_COMP("ana")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbCompana
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
#define mccompcurname  ana
#define mccompcurtype  Monochromator_flat
#define mccompcurindex 28
#define mos_rms_y mccana_mos_rms_y
#define mos_rms_z mccana_mos_rms_z
#define mos_rms_max mccana_mos_rms_max
#define mono_Q mccana_mono_Q
{   /* Declarations of ana=Monochromator_flat() SETTING parameters. */
MCNUM zmin = mccana_zmin;
MCNUM zmax = mccana_zmax;
MCNUM ymin = mccana_ymin;
MCNUM ymax = mccana_ymax;
MCNUM zwidth = mccana_zwidth;
MCNUM yheight = mccana_yheight;
MCNUM mosaich = mccana_mosaich;
MCNUM mosaicv = mccana_mosaicv;
MCNUM r0 = mccana_r0;
MCNUM Q = mccana_Q;
MCNUM DM = mccana_DM;
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
#line 15316 "./linup-6.c"
}   /* End of ana=Monochromator_flat() SETTING parameter declarations. */
#undef mono_Q
#undef mos_rms_max
#undef mos_rms_z
#undef mos_rms_y
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompana:
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

  /* TRACE Component a4 [29] */
  mccoordschange(mcposra4, mcrotra4,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component a4 (without coords transformations) */
  mcJumpTrace_a4:
  SIG_MESSAGE("a4 (Trace)");
  mcDEBUG_COMP("a4")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbCompa4
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
#define mccompcurname  a4
#define mccompcurtype  Arm
#define mccompcurindex 29
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompa4:
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

  /* TRACE Component c3 [30] */
  mccoordschange(mcposrc3, mcrotrc3,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component c3 (without coords transformations) */
  mcJumpTrace_c3:
  SIG_MESSAGE("c3 (Trace)");
  mcDEBUG_COMP("c3")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbCompc3
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
#define mccompcurname  c3
#define mccompcurtype  Collimator_linear
#define mccompcurindex 30
#define slope mccc3_slope
#define slopeV mccc3_slopeV
{   /* Declarations of c3=Collimator_linear() SETTING parameters. */
MCNUM xmin = mccc3_xmin;
MCNUM xmax = mccc3_xmax;
MCNUM ymin = mccc3_ymin;
MCNUM ymax = mccc3_ymax;
MCNUM xwidth = mccc3_xwidth;
MCNUM yheight = mccc3_yheight;
MCNUM length = mccc3_length;
MCNUM divergence = mccc3_divergence;
MCNUM transmission = mccc3_transmission;
MCNUM divergenceV = mccc3_divergenceV;
#line 69 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Collimator_linear.comp"
{
    double phi, dt;

    PROP_Z0;
    if (x<xmin || x>xmax || y<ymin || y>ymax)
      ABSORB;
    dt = length/vz;
    PROP_DT(dt);
    if (x<xmin || x>xmax || y<ymin || y>ymax)
      ABSORB;

    if(slope > 0.0)
    {
      phi = fabs(vx/vz);
      if (phi > slope)
        ABSORB;
      else
        p *= transmission*(1.0 - phi/slope);
      SCATTER;
    }
    if (slopeV > 0) {
      phi = fabs(vy/vz);
      if (phi > slopeV)
        ABSORB;
      else
        p *= transmission*(1.0 - phi/slopeV);
      SCATTER;
    }
}
#line 15571 "./linup-6.c"
}   /* End of c3=Collimator_linear() SETTING parameter declarations. */
#undef slopeV
#undef slope
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompc3:
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

  /* TRACE Component sng [31] */
  mccoordschange(mcposrsng, mcrotrsng,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component sng (without coords transformations) */
  mcJumpTrace_sng:
  SIG_MESSAGE("sng (Trace)");
  mcDEBUG_COMP("sng")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbCompsng
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
#define mccompcurname  sng
#define mccompcurtype  Monitor
#define mccompcurindex 31
#define Nsum mccsng_Nsum
#define psum mccsng_psum
#define p2sum mccsng_p2sum
{   /* Declarations of sng=Monitor() SETTING parameters. */
MCNUM xmin = mccsng_xmin;
MCNUM xmax = mccsng_xmax;
MCNUM ymin = mccsng_ymin;
MCNUM ymax = mccsng_ymax;
MCNUM xwidth = mccsng_xwidth;
MCNUM yheight = mccsng_yheight;
MCNUM restore_neutron = mccsng_restore_neutron;
#line 74 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/Monitor.comp"
{
    PROP_Z0;
    if (x>xmin && x<xmax && y>ymin && y<ymax)
    {
      Nsum++;
      psum += p;
      p2sum += p*p;
      SCATTER;
    }
    if (restore_neutron) {
      RESTORE_NEUTRON(INDEX_CURRENT_COMP, x, y, z, vx, vy, vz, t, sx, sy, sz, p);
    }
}
#line 15703 "./linup-6.c"
}   /* End of sng=Monitor() SETTING parameter declarations. */
#undef p2sum
#undef psum
#undef Nsum
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompsng:
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

  /* TRACE Component emon2 [32] */
  mccoordschange(mcposremon2, mcrotremon2,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component emon2 (without coords transformations) */
  mcJumpTrace_emon2:
  SIG_MESSAGE("emon2 (Trace)");
  mcDEBUG_COMP("emon2")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
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

#define mcabsorbComp mcabsorbCompemon2
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
#define mccompcurname  emon2
#define mccompcurtype  E_monitor
#define mccompcurindex 32
#define nE mccemon2_nE
#define E_N mccemon2_E_N
#define E_p mccemon2_E_p
#define E_p2 mccemon2_E_p2
#define S_p mccemon2_S_p
#define S_pE mccemon2_S_pE
#define S_pE2 mccemon2_S_pE2
{   /* Declarations of emon2=E_monitor() SETTING parameters. */
char* filename = mccemon2_filename;
MCNUM xmin = mccemon2_xmin;
MCNUM xmax = mccemon2_xmax;
MCNUM ymin = mccemon2_ymin;
MCNUM ymax = mccemon2_ymax;
MCNUM xwidth = mccemon2_xwidth;
MCNUM yheight = mccemon2_yheight;
MCNUM Emin = mccemon2_Emin;
MCNUM Emax = mccemon2_Emax;
MCNUM restore_neutron = mccemon2_restore_neutron;
int nowritefile = mccemon2_nowritefile;
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
#line 15857 "./linup-6.c"
}   /* End of emon2=E_monitor() SETTING parameter declarations. */
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
  mcabsorbCompemon2:
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

  mcabsorbAll:
  /* SPLIT loops in reverse order */
  if (mcSplit_sample && mcSplit_sample < (10)) {
    goto mcJumpTrace_sample;
  }
    else mcSplit_sample=0;

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
#define mccompcurindex 1
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
#line 15974 "./linup-6.c"
}   /* End of a1=Progress_bar() SETTING parameter declarations. */
#undef CurrentTime
#undef EndTime
#undef StartTime
#undef IntermediateCnts
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'mon'. */
  SIG_MESSAGE("mon (Save)");
#define mccompcurname  mon
#define mccompcurtype  Monitor
#define mccompcurindex 22
#define Nsum mccmon_Nsum
#define psum mccmon_psum
#define p2sum mccmon_p2sum
{   /* Declarations of mon=Monitor() SETTING parameters. */
MCNUM xmin = mccmon_xmin;
MCNUM xmax = mccmon_xmax;
MCNUM ymin = mccmon_ymin;
MCNUM ymax = mccmon_ymax;
MCNUM xwidth = mccmon_xwidth;
MCNUM yheight = mccmon_yheight;
MCNUM restore_neutron = mccmon_restore_neutron;
#line 89 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/Monitor.comp"
{
    char title[1024];
    sprintf(title, "Single monitor %s", NAME_CURRENT_COMP);
    DETECTOR_OUT_0D(title, Nsum, psum, p2sum);
}
#line 16006 "./linup-6.c"
}   /* End of mon=Monitor() SETTING parameter declarations. */
#undef p2sum
#undef psum
#undef Nsum
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'emon1'. */
  SIG_MESSAGE("emon1 (Save)");
#define mccompcurname  emon1
#define mccompcurtype  E_monitor
#define mccompcurindex 23
#define nE mccemon1_nE
#define E_N mccemon1_E_N
#define E_p mccemon1_E_p
#define E_p2 mccemon1_E_p2
#define S_p mccemon1_S_p
#define S_pE mccemon1_S_pE
#define S_pE2 mccemon1_S_pE2
{   /* Declarations of emon1=E_monitor() SETTING parameters. */
char* filename = mccemon1_filename;
MCNUM xmin = mccemon1_xmin;
MCNUM xmax = mccemon1_xmax;
MCNUM ymin = mccemon1_ymin;
MCNUM ymax = mccemon1_ymax;
MCNUM xwidth = mccemon1_xwidth;
MCNUM yheight = mccemon1_yheight;
MCNUM Emin = mccemon1_Emin;
MCNUM Emax = mccemon1_Emax;
MCNUM restore_neutron = mccemon1_restore_neutron;
int nowritefile = mccemon1_nowritefile;
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
#line 16053 "./linup-6.c"
}   /* End of emon1=E_monitor() SETTING parameter declarations. */
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

  /* User SAVE code for component 'focus_check'. */
  SIG_MESSAGE("focus_check (Save)");
#define mccompcurname  focus_check
#define mccompcurtype  PSD_monitor
#define mccompcurindex 26
#define PSD_N mccfocus_check_PSD_N
#define PSD_p mccfocus_check_PSD_p
#define PSD_p2 mccfocus_check_PSD_p2
{   /* Declarations of focus_check=PSD_monitor() SETTING parameters. */
int nx = mccfocus_check_nx;
int ny = mccfocus_check_ny;
char* filename = mccfocus_check_filename;
MCNUM xmin = mccfocus_check_xmin;
MCNUM xmax = mccfocus_check_xmax;
MCNUM ymin = mccfocus_check_ymin;
MCNUM ymax = mccfocus_check_ymax;
MCNUM xwidth = mccfocus_check_xwidth;
MCNUM yheight = mccfocus_check_yheight;
MCNUM restore_neutron = mccfocus_check_restore_neutron;
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
#line 16096 "./linup-6.c"
}   /* End of focus_check=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'sng'. */
  SIG_MESSAGE("sng (Save)");
#define mccompcurname  sng
#define mccompcurtype  Monitor
#define mccompcurindex 31
#define Nsum mccsng_Nsum
#define psum mccsng_psum
#define p2sum mccsng_p2sum
{   /* Declarations of sng=Monitor() SETTING parameters. */
MCNUM xmin = mccsng_xmin;
MCNUM xmax = mccsng_xmax;
MCNUM ymin = mccsng_ymin;
MCNUM ymax = mccsng_ymax;
MCNUM xwidth = mccsng_xwidth;
MCNUM yheight = mccsng_yheight;
MCNUM restore_neutron = mccsng_restore_neutron;
#line 89 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/Monitor.comp"
{
    char title[1024];
    sprintf(title, "Single monitor %s", NAME_CURRENT_COMP);
    DETECTOR_OUT_0D(title, Nsum, psum, p2sum);
}
#line 16127 "./linup-6.c"
}   /* End of sng=Monitor() SETTING parameter declarations. */
#undef p2sum
#undef psum
#undef Nsum
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'emon2'. */
  SIG_MESSAGE("emon2 (Save)");
#define mccompcurname  emon2
#define mccompcurtype  E_monitor
#define mccompcurindex 32
#define nE mccemon2_nE
#define E_N mccemon2_E_N
#define E_p mccemon2_E_p
#define E_p2 mccemon2_E_p2
#define S_p mccemon2_S_p
#define S_pE mccemon2_S_pE
#define S_pE2 mccemon2_S_pE2
{   /* Declarations of emon2=E_monitor() SETTING parameters. */
char* filename = mccemon2_filename;
MCNUM xmin = mccemon2_xmin;
MCNUM xmax = mccemon2_xmax;
MCNUM ymin = mccemon2_ymin;
MCNUM ymax = mccemon2_ymax;
MCNUM xwidth = mccemon2_xwidth;
MCNUM yheight = mccemon2_yheight;
MCNUM Emin = mccemon2_Emin;
MCNUM Emax = mccemon2_Emax;
MCNUM restore_neutron = mccemon2_restore_neutron;
int nowritefile = mccemon2_nowritefile;
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
#line 16174 "./linup-6.c"
}   /* End of emon2=E_monitor() SETTING parameter declarations. */
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

  /* User FINALLY code for component 'a1'. */
  SIG_MESSAGE("a1 (Finally)");
#define mccompcurname  a1
#define mccompcurtype  Progress_bar
#define mccompcurindex 1
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
#line 16221 "./linup-6.c"
}   /* End of a1=Progress_bar() SETTING parameter declarations. */
#undef CurrentTime
#undef EndTime
#undef StartTime
#undef IntermediateCnts
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[1]) fprintf(stderr, "Warning: No neutron could reach Component[1] a1\n");
    if (mcAbsorbProp[1]) fprintf(stderr, "Warning: %g events were removed in Component[1] a1=Progress_bar()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[1]);
    if (!mcNCounter[2]) fprintf(stderr, "Warning: No neutron could reach Component[2] source\n");
    if (mcAbsorbProp[2]) fprintf(stderr, "Warning: %g events were removed in Component[2] source=Source_simple()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[2]);
    if (!mcNCounter[3]) fprintf(stderr, "Warning: No neutron could reach Component[3] slit1\n");
    if (mcAbsorbProp[3]) fprintf(stderr, "Warning: %g events were removed in Component[3] slit1=Slit()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[3]);
    if (!mcNCounter[4]) fprintf(stderr, "Warning: No neutron could reach Component[4] slit2\n");
    if (mcAbsorbProp[4]) fprintf(stderr, "Warning: %g events were removed in Component[4] slit2=Slit()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[4]);
    if (!mcNCounter[5]) fprintf(stderr, "Warning: No neutron could reach Component[5] slit3\n");
    if (mcAbsorbProp[5]) fprintf(stderr, "Warning: %g events were removed in Component[5] slit3=Slit()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[5]);
    if (!mcNCounter[6]) fprintf(stderr, "Warning: No neutron could reach Component[6] focus_mono\n");
    if (mcAbsorbProp[6]) fprintf(stderr, "Warning: %g events were removed in Component[6] focus_mono=Arm()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[6]);
    if (!mcNCounter[7]) fprintf(stderr, "Warning: No neutron could reach Component[7] m0\n");
    if (mcAbsorbProp[7]) fprintf(stderr, "Warning: %g events were removed in Component[7] m0=Monochromator_flat()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[7]);
    if (!mcNCounter[8]) fprintf(stderr, "Warning: No neutron could reach Component[8] m1\n");
    if (mcAbsorbProp[8]) fprintf(stderr, "Warning: %g events were removed in Component[8] m1=Monochromator_flat()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[8]);
    if (!mcNCounter[9]) fprintf(stderr, "Warning: No neutron could reach Component[9] m2\n");
    if (mcAbsorbProp[9]) fprintf(stderr, "Warning: %g events were removed in Component[9] m2=Monochromator_flat()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[9]);
    if (!mcNCounter[10]) fprintf(stderr, "Warning: No neutron could reach Component[10] m3\n");
    if (mcAbsorbProp[10]) fprintf(stderr, "Warning: %g events were removed in Component[10] m3=Monochromator_flat()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[10]);
    if (!mcNCounter[11]) fprintf(stderr, "Warning: No neutron could reach Component[11] m4\n");
    if (mcAbsorbProp[11]) fprintf(stderr, "Warning: %g events were removed in Component[11] m4=Monochromator_flat()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[11]);
    if (!mcNCounter[12]) fprintf(stderr, "Warning: No neutron could reach Component[12] m5\n");
    if (mcAbsorbProp[12]) fprintf(stderr, "Warning: %g events were removed in Component[12] m5=Monochromator_flat()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[12]);
    if (!mcNCounter[13]) fprintf(stderr, "Warning: No neutron could reach Component[13] m6\n");
    if (mcAbsorbProp[13]) fprintf(stderr, "Warning: %g events were removed in Component[13] m6=Monochromator_flat()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[13]);
    if (!mcNCounter[14]) fprintf(stderr, "Warning: No neutron could reach Component[14] m7\n");
    if (mcAbsorbProp[14]) fprintf(stderr, "Warning: %g events were removed in Component[14] m7=Monochromator_flat()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[14]);
    if (!mcNCounter[15]) fprintf(stderr, "Warning: No neutron could reach Component[15] a2\n");
    if (mcAbsorbProp[15]) fprintf(stderr, "Warning: %g events were removed in Component[15] a2=Arm()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[15]);
    if (!mcNCounter[16]) fprintf(stderr, "Warning: No neutron could reach Component[16] slitMS1\n");
    if (mcAbsorbProp[16]) fprintf(stderr, "Warning: %g events were removed in Component[16] slitMS1=Slit()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[16]);
    if (!mcNCounter[17]) fprintf(stderr, "Warning: No neutron could reach Component[17] slitMS2\n");
    if (mcAbsorbProp[17]) fprintf(stderr, "Warning: %g events were removed in Component[17] slitMS2=Slit()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[17]);
    if (!mcNCounter[18]) fprintf(stderr, "Warning: No neutron could reach Component[18] c1\n");
    if (mcAbsorbProp[18]) fprintf(stderr, "Warning: %g events were removed in Component[18] c1=Collimator_linear()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[18]);
    if (!mcNCounter[19]) fprintf(stderr, "Warning: No neutron could reach Component[19] slitMS3\n");
    if (mcAbsorbProp[19]) fprintf(stderr, "Warning: %g events were removed in Component[19] slitMS3=Slit()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[19]);
    if (!mcNCounter[20]) fprintf(stderr, "Warning: No neutron could reach Component[20] slitMS4\n");
    if (mcAbsorbProp[20]) fprintf(stderr, "Warning: %g events were removed in Component[20] slitMS4=Slit()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[20]);
    if (!mcNCounter[21]) fprintf(stderr, "Warning: No neutron could reach Component[21] slitMS5\n");
    if (mcAbsorbProp[21]) fprintf(stderr, "Warning: %g events were removed in Component[21] slitMS5=Slit()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[21]);
    if (!mcNCounter[22]) fprintf(stderr, "Warning: No neutron could reach Component[22] mon\n");
    if (mcAbsorbProp[22]) fprintf(stderr, "Warning: %g events were removed in Component[22] mon=Monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[22]);
    if (!mcNCounter[23]) fprintf(stderr, "Warning: No neutron could reach Component[23] emon1\n");
    if (mcAbsorbProp[23]) fprintf(stderr, "Warning: %g events were removed in Component[23] emon1=E_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[23]);
    if (!mcNCounter[24]) fprintf(stderr, "Warning: No neutron could reach Component[24] sample\n");
    if (mcNCounter[24] < 1000*(10)) fprintf(stderr, 
"Warning: Number of events %g reaching SPLIT position Component[24] sample=V_sample()\n"
"         is probably too low. Increase Ncount.\n", mcNCounter[24]);

    if (mcAbsorbProp[24]) fprintf(stderr, "Warning: %g events were removed in Component[24] sample=V_sample()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[24]);
    if (!mcNCounter[25]) fprintf(stderr, "Warning: No neutron could reach Component[25] a3\n");
    if (mcAbsorbProp[25]) fprintf(stderr, "Warning: %g events were removed in Component[25] a3=Arm()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[25]);
  /* User FINALLY code for component 'focus_check'. */
  SIG_MESSAGE("focus_check (Finally)");
#define mccompcurname  focus_check
#define mccompcurtype  PSD_monitor
#define mccompcurindex 26
#define PSD_N mccfocus_check_PSD_N
#define PSD_p mccfocus_check_PSD_p
#define PSD_p2 mccfocus_check_PSD_p2
{   /* Declarations of focus_check=PSD_monitor() SETTING parameters. */
int nx = mccfocus_check_nx;
int ny = mccfocus_check_ny;
char* filename = mccfocus_check_filename;
MCNUM xmin = mccfocus_check_xmin;
MCNUM xmax = mccfocus_check_xmax;
MCNUM ymin = mccfocus_check_ymin;
MCNUM ymax = mccfocus_check_ymax;
MCNUM xwidth = mccfocus_check_xwidth;
MCNUM yheight = mccfocus_check_yheight;
MCNUM restore_neutron = mccfocus_check_restore_neutron;
#line 122 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
{
  destroy_darr2d(PSD_N);
  destroy_darr2d(PSD_p);
  destroy_darr2d(PSD_p2);
}
#line 16307 "./linup-6.c"
}   /* End of focus_check=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[26]) fprintf(stderr, "Warning: No neutron could reach Component[26] focus_check\n");
    if (mcAbsorbProp[26]) fprintf(stderr, "Warning: %g events were removed in Component[26] focus_check=PSD_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[26]);
    if (!mcNCounter[27]) fprintf(stderr, "Warning: No neutron could reach Component[27] c2\n");
    if (mcAbsorbProp[27]) fprintf(stderr, "Warning: %g events were removed in Component[27] c2=Collimator_linear()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[27]);
    if (!mcNCounter[28]) fprintf(stderr, "Warning: No neutron could reach Component[28] ana\n");
    if (mcAbsorbProp[28]) fprintf(stderr, "Warning: %g events were removed in Component[28] ana=Monochromator_flat()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[28]);
    if (!mcNCounter[29]) fprintf(stderr, "Warning: No neutron could reach Component[29] a4\n");
    if (mcAbsorbProp[29]) fprintf(stderr, "Warning: %g events were removed in Component[29] a4=Arm()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[29]);
    if (!mcNCounter[30]) fprintf(stderr, "Warning: No neutron could reach Component[30] c3\n");
    if (mcAbsorbProp[30]) fprintf(stderr, "Warning: %g events were removed in Component[30] c3=Collimator_linear()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[30]);
    if (!mcNCounter[31]) fprintf(stderr, "Warning: No neutron could reach Component[31] sng\n");
    if (mcAbsorbProp[31]) fprintf(stderr, "Warning: %g events were removed in Component[31] sng=Monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[31]);
    if (!mcNCounter[32]) fprintf(stderr, "Warning: No neutron could reach Component[32] emon2\n");
    if (mcAbsorbProp[32]) fprintf(stderr, "Warning: %g events were removed in Component[32] emon2=E_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[32]);
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
#define mccompcurindex 1
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
#line 16364 "./linup-6.c"
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
#define mccompcurtype  Source_simple
#define mccompcurindex 2
#define pmul mccsource_pmul
#define square mccsource_square
#define srcArea mccsource_srcArea
{   /* Declarations of source=Source_simple() SETTING parameters. */
MCNUM radius = mccsource_radius;
MCNUM yheight = mccsource_yheight;
MCNUM xwidth = mccsource_xwidth;
MCNUM dist = mccsource_dist;
MCNUM focus_xw = mccsource_focus_xw;
MCNUM focus_yh = mccsource_focus_yh;
MCNUM E0 = mccsource_E0;
MCNUM dE = mccsource_dE;
MCNUM lambda0 = mccsource_lambda0;
MCNUM dlambda = mccsource_dlambda;
MCNUM flux = mccsource_flux;
MCNUM gauss = mccsource_gauss;
int target_index = mccsource_target_index;
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
#line 16413 "./linup-6.c"
}   /* End of source=Source_simple() SETTING parameter declarations. */
#undef srcArea
#undef square
#undef pmul
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'slit1'. */
  SIG_MESSAGE("slit1 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "slit1");
#define mccompcurname  slit1
#define mccompcurtype  Slit
#define mccompcurindex 3
{   /* Declarations of slit1=Slit() SETTING parameters. */
MCNUM xmin = mccslit1_xmin;
MCNUM xmax = mccslit1_xmax;
MCNUM ymin = mccslit1_ymin;
MCNUM ymax = mccslit1_ymax;
MCNUM radius = mccslit1_radius;
MCNUM xwidth = mccslit1_xwidth;
MCNUM yheight = mccslit1_yheight;
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
#line 16459 "./linup-6.c"
}   /* End of slit1=Slit() SETTING parameter declarations. */
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'slit2'. */
  SIG_MESSAGE("slit2 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "slit2");
#define mccompcurname  slit2
#define mccompcurtype  Slit
#define mccompcurindex 4
{   /* Declarations of slit2=Slit() SETTING parameters. */
MCNUM xmin = mccslit2_xmin;
MCNUM xmax = mccslit2_xmax;
MCNUM ymin = mccslit2_ymin;
MCNUM ymax = mccslit2_ymax;
MCNUM radius = mccslit2_radius;
MCNUM xwidth = mccslit2_xwidth;
MCNUM yheight = mccslit2_yheight;
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
#line 16502 "./linup-6.c"
}   /* End of slit2=Slit() SETTING parameter declarations. */
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'slit3'. */
  SIG_MESSAGE("slit3 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "slit3");
#define mccompcurname  slit3
#define mccompcurtype  Slit
#define mccompcurindex 5
{   /* Declarations of slit3=Slit() SETTING parameters. */
MCNUM xmin = mccslit3_xmin;
MCNUM xmax = mccslit3_xmax;
MCNUM ymin = mccslit3_ymin;
MCNUM ymax = mccslit3_ymax;
MCNUM radius = mccslit3_radius;
MCNUM xwidth = mccslit3_xwidth;
MCNUM yheight = mccslit3_yheight;
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
#line 16545 "./linup-6.c"
}   /* End of slit3=Slit() SETTING parameter declarations. */
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'focus_mono'. */
  SIG_MESSAGE("focus_mono (McDisplay)");
  printf("MCDISPLAY: component %s\n", "focus_mono");
#define mccompcurname  focus_mono
#define mccompcurtype  Arm
#define mccompcurindex 6
#line 40 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Arm.comp"
{
  /* A bit ugly; hard-coded dimensions. */
  
  line(0,0,0,0.2,0,0);
  line(0,0,0,0,0.2,0);
  line(0,0,0,0,0,0.2);
}
#line 16565 "./linup-6.c"
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'm0'. */
  SIG_MESSAGE("m0 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "m0");
#define mccompcurname  m0
#define mccompcurtype  Monochromator_flat
#define mccompcurindex 7
#define mos_rms_y mccm0_mos_rms_y
#define mos_rms_z mccm0_mos_rms_z
#define mos_rms_max mccm0_mos_rms_max
#define mono_Q mccm0_mono_Q
{   /* Declarations of m0=Monochromator_flat() SETTING parameters. */
MCNUM zmin = mccm0_zmin;
MCNUM zmax = mccm0_zmax;
MCNUM ymin = mccm0_ymin;
MCNUM ymax = mccm0_ymax;
MCNUM zwidth = mccm0_zwidth;
MCNUM yheight = mccm0_yheight;
MCNUM mosaich = mccm0_mosaich;
MCNUM mosaicv = mccm0_mosaicv;
MCNUM r0 = mccm0_r0;
MCNUM Q = mccm0_Q;
MCNUM DM = mccm0_DM;
#line 254 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Monochromator_flat.comp"
{
  
  multiline(5, 0.0, (double)ymin, (double)zmin,
               0.0, (double)ymax, (double)zmin,
               0.0, (double)ymax, (double)zmax,
               0.0, (double)ymin, (double)zmax,
               0.0, (double)ymin, (double)zmin);
}
#line 16601 "./linup-6.c"
}   /* End of m0=Monochromator_flat() SETTING parameter declarations. */
#undef mono_Q
#undef mos_rms_max
#undef mos_rms_z
#undef mos_rms_y
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'm1'. */
  SIG_MESSAGE("m1 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "m1");
#define mccompcurname  m1
#define mccompcurtype  Monochromator_flat
#define mccompcurindex 8
#define mos_rms_y mccm1_mos_rms_y
#define mos_rms_z mccm1_mos_rms_z
#define mos_rms_max mccm1_mos_rms_max
#define mono_Q mccm1_mono_Q
{   /* Declarations of m1=Monochromator_flat() SETTING parameters. */
MCNUM zmin = mccm1_zmin;
MCNUM zmax = mccm1_zmax;
MCNUM ymin = mccm1_ymin;
MCNUM ymax = mccm1_ymax;
MCNUM zwidth = mccm1_zwidth;
MCNUM yheight = mccm1_yheight;
MCNUM mosaich = mccm1_mosaich;
MCNUM mosaicv = mccm1_mosaicv;
MCNUM r0 = mccm1_r0;
MCNUM Q = mccm1_Q;
MCNUM DM = mccm1_DM;
#line 254 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Monochromator_flat.comp"
{
  
  multiline(5, 0.0, (double)ymin, (double)zmin,
               0.0, (double)ymax, (double)zmin,
               0.0, (double)ymax, (double)zmax,
               0.0, (double)ymin, (double)zmax,
               0.0, (double)ymin, (double)zmin);
}
#line 16642 "./linup-6.c"
}   /* End of m1=Monochromator_flat() SETTING parameter declarations. */
#undef mono_Q
#undef mos_rms_max
#undef mos_rms_z
#undef mos_rms_y
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'm2'. */
  SIG_MESSAGE("m2 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "m2");
#define mccompcurname  m2
#define mccompcurtype  Monochromator_flat
#define mccompcurindex 9
#define mos_rms_y mccm2_mos_rms_y
#define mos_rms_z mccm2_mos_rms_z
#define mos_rms_max mccm2_mos_rms_max
#define mono_Q mccm2_mono_Q
{   /* Declarations of m2=Monochromator_flat() SETTING parameters. */
MCNUM zmin = mccm2_zmin;
MCNUM zmax = mccm2_zmax;
MCNUM ymin = mccm2_ymin;
MCNUM ymax = mccm2_ymax;
MCNUM zwidth = mccm2_zwidth;
MCNUM yheight = mccm2_yheight;
MCNUM mosaich = mccm2_mosaich;
MCNUM mosaicv = mccm2_mosaicv;
MCNUM r0 = mccm2_r0;
MCNUM Q = mccm2_Q;
MCNUM DM = mccm2_DM;
#line 254 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Monochromator_flat.comp"
{
  
  multiline(5, 0.0, (double)ymin, (double)zmin,
               0.0, (double)ymax, (double)zmin,
               0.0, (double)ymax, (double)zmax,
               0.0, (double)ymin, (double)zmax,
               0.0, (double)ymin, (double)zmin);
}
#line 16683 "./linup-6.c"
}   /* End of m2=Monochromator_flat() SETTING parameter declarations. */
#undef mono_Q
#undef mos_rms_max
#undef mos_rms_z
#undef mos_rms_y
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'm3'. */
  SIG_MESSAGE("m3 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "m3");
#define mccompcurname  m3
#define mccompcurtype  Monochromator_flat
#define mccompcurindex 10
#define mos_rms_y mccm3_mos_rms_y
#define mos_rms_z mccm3_mos_rms_z
#define mos_rms_max mccm3_mos_rms_max
#define mono_Q mccm3_mono_Q
{   /* Declarations of m3=Monochromator_flat() SETTING parameters. */
MCNUM zmin = mccm3_zmin;
MCNUM zmax = mccm3_zmax;
MCNUM ymin = mccm3_ymin;
MCNUM ymax = mccm3_ymax;
MCNUM zwidth = mccm3_zwidth;
MCNUM yheight = mccm3_yheight;
MCNUM mosaich = mccm3_mosaich;
MCNUM mosaicv = mccm3_mosaicv;
MCNUM r0 = mccm3_r0;
MCNUM Q = mccm3_Q;
MCNUM DM = mccm3_DM;
#line 254 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Monochromator_flat.comp"
{
  
  multiline(5, 0.0, (double)ymin, (double)zmin,
               0.0, (double)ymax, (double)zmin,
               0.0, (double)ymax, (double)zmax,
               0.0, (double)ymin, (double)zmax,
               0.0, (double)ymin, (double)zmin);
}
#line 16724 "./linup-6.c"
}   /* End of m3=Monochromator_flat() SETTING parameter declarations. */
#undef mono_Q
#undef mos_rms_max
#undef mos_rms_z
#undef mos_rms_y
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'm4'. */
  SIG_MESSAGE("m4 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "m4");
#define mccompcurname  m4
#define mccompcurtype  Monochromator_flat
#define mccompcurindex 11
#define mos_rms_y mccm4_mos_rms_y
#define mos_rms_z mccm4_mos_rms_z
#define mos_rms_max mccm4_mos_rms_max
#define mono_Q mccm4_mono_Q
{   /* Declarations of m4=Monochromator_flat() SETTING parameters. */
MCNUM zmin = mccm4_zmin;
MCNUM zmax = mccm4_zmax;
MCNUM ymin = mccm4_ymin;
MCNUM ymax = mccm4_ymax;
MCNUM zwidth = mccm4_zwidth;
MCNUM yheight = mccm4_yheight;
MCNUM mosaich = mccm4_mosaich;
MCNUM mosaicv = mccm4_mosaicv;
MCNUM r0 = mccm4_r0;
MCNUM Q = mccm4_Q;
MCNUM DM = mccm4_DM;
#line 254 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Monochromator_flat.comp"
{
  
  multiline(5, 0.0, (double)ymin, (double)zmin,
               0.0, (double)ymax, (double)zmin,
               0.0, (double)ymax, (double)zmax,
               0.0, (double)ymin, (double)zmax,
               0.0, (double)ymin, (double)zmin);
}
#line 16765 "./linup-6.c"
}   /* End of m4=Monochromator_flat() SETTING parameter declarations. */
#undef mono_Q
#undef mos_rms_max
#undef mos_rms_z
#undef mos_rms_y
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'm5'. */
  SIG_MESSAGE("m5 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "m5");
#define mccompcurname  m5
#define mccompcurtype  Monochromator_flat
#define mccompcurindex 12
#define mos_rms_y mccm5_mos_rms_y
#define mos_rms_z mccm5_mos_rms_z
#define mos_rms_max mccm5_mos_rms_max
#define mono_Q mccm5_mono_Q
{   /* Declarations of m5=Monochromator_flat() SETTING parameters. */
MCNUM zmin = mccm5_zmin;
MCNUM zmax = mccm5_zmax;
MCNUM ymin = mccm5_ymin;
MCNUM ymax = mccm5_ymax;
MCNUM zwidth = mccm5_zwidth;
MCNUM yheight = mccm5_yheight;
MCNUM mosaich = mccm5_mosaich;
MCNUM mosaicv = mccm5_mosaicv;
MCNUM r0 = mccm5_r0;
MCNUM Q = mccm5_Q;
MCNUM DM = mccm5_DM;
#line 254 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Monochromator_flat.comp"
{
  
  multiline(5, 0.0, (double)ymin, (double)zmin,
               0.0, (double)ymax, (double)zmin,
               0.0, (double)ymax, (double)zmax,
               0.0, (double)ymin, (double)zmax,
               0.0, (double)ymin, (double)zmin);
}
#line 16806 "./linup-6.c"
}   /* End of m5=Monochromator_flat() SETTING parameter declarations. */
#undef mono_Q
#undef mos_rms_max
#undef mos_rms_z
#undef mos_rms_y
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'm6'. */
  SIG_MESSAGE("m6 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "m6");
#define mccompcurname  m6
#define mccompcurtype  Monochromator_flat
#define mccompcurindex 13
#define mos_rms_y mccm6_mos_rms_y
#define mos_rms_z mccm6_mos_rms_z
#define mos_rms_max mccm6_mos_rms_max
#define mono_Q mccm6_mono_Q
{   /* Declarations of m6=Monochromator_flat() SETTING parameters. */
MCNUM zmin = mccm6_zmin;
MCNUM zmax = mccm6_zmax;
MCNUM ymin = mccm6_ymin;
MCNUM ymax = mccm6_ymax;
MCNUM zwidth = mccm6_zwidth;
MCNUM yheight = mccm6_yheight;
MCNUM mosaich = mccm6_mosaich;
MCNUM mosaicv = mccm6_mosaicv;
MCNUM r0 = mccm6_r0;
MCNUM Q = mccm6_Q;
MCNUM DM = mccm6_DM;
#line 254 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Monochromator_flat.comp"
{
  
  multiline(5, 0.0, (double)ymin, (double)zmin,
               0.0, (double)ymax, (double)zmin,
               0.0, (double)ymax, (double)zmax,
               0.0, (double)ymin, (double)zmax,
               0.0, (double)ymin, (double)zmin);
}
#line 16847 "./linup-6.c"
}   /* End of m6=Monochromator_flat() SETTING parameter declarations. */
#undef mono_Q
#undef mos_rms_max
#undef mos_rms_z
#undef mos_rms_y
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'm7'. */
  SIG_MESSAGE("m7 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "m7");
#define mccompcurname  m7
#define mccompcurtype  Monochromator_flat
#define mccompcurindex 14
#define mos_rms_y mccm7_mos_rms_y
#define mos_rms_z mccm7_mos_rms_z
#define mos_rms_max mccm7_mos_rms_max
#define mono_Q mccm7_mono_Q
{   /* Declarations of m7=Monochromator_flat() SETTING parameters. */
MCNUM zmin = mccm7_zmin;
MCNUM zmax = mccm7_zmax;
MCNUM ymin = mccm7_ymin;
MCNUM ymax = mccm7_ymax;
MCNUM zwidth = mccm7_zwidth;
MCNUM yheight = mccm7_yheight;
MCNUM mosaich = mccm7_mosaich;
MCNUM mosaicv = mccm7_mosaicv;
MCNUM r0 = mccm7_r0;
MCNUM Q = mccm7_Q;
MCNUM DM = mccm7_DM;
#line 254 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Monochromator_flat.comp"
{
  
  multiline(5, 0.0, (double)ymin, (double)zmin,
               0.0, (double)ymax, (double)zmin,
               0.0, (double)ymax, (double)zmax,
               0.0, (double)ymin, (double)zmax,
               0.0, (double)ymin, (double)zmin);
}
#line 16888 "./linup-6.c"
}   /* End of m7=Monochromator_flat() SETTING parameter declarations. */
#undef mono_Q
#undef mos_rms_max
#undef mos_rms_z
#undef mos_rms_y
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'a2'. */
  SIG_MESSAGE("a2 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "a2");
#define mccompcurname  a2
#define mccompcurtype  Arm
#define mccompcurindex 15
#line 40 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Arm.comp"
{
  /* A bit ugly; hard-coded dimensions. */
  
  line(0,0,0,0.2,0,0);
  line(0,0,0,0,0.2,0);
  line(0,0,0,0,0,0.2);
}
#line 16912 "./linup-6.c"
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'slitMS1'. */
  SIG_MESSAGE("slitMS1 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "slitMS1");
#define mccompcurname  slitMS1
#define mccompcurtype  Slit
#define mccompcurindex 16
{   /* Declarations of slitMS1=Slit() SETTING parameters. */
MCNUM xmin = mccslitMS1_xmin;
MCNUM xmax = mccslitMS1_xmax;
MCNUM ymin = mccslitMS1_ymin;
MCNUM ymax = mccslitMS1_ymax;
MCNUM radius = mccslitMS1_radius;
MCNUM xwidth = mccslitMS1_xwidth;
MCNUM yheight = mccslitMS1_yheight;
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
#line 16954 "./linup-6.c"
}   /* End of slitMS1=Slit() SETTING parameter declarations. */
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'slitMS2'. */
  SIG_MESSAGE("slitMS2 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "slitMS2");
#define mccompcurname  slitMS2
#define mccompcurtype  Slit
#define mccompcurindex 17
{   /* Declarations of slitMS2=Slit() SETTING parameters. */
MCNUM xmin = mccslitMS2_xmin;
MCNUM xmax = mccslitMS2_xmax;
MCNUM ymin = mccslitMS2_ymin;
MCNUM ymax = mccslitMS2_ymax;
MCNUM radius = mccslitMS2_radius;
MCNUM xwidth = mccslitMS2_xwidth;
MCNUM yheight = mccslitMS2_yheight;
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
#line 16997 "./linup-6.c"
}   /* End of slitMS2=Slit() SETTING parameter declarations. */
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'c1'. */
  SIG_MESSAGE("c1 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "c1");
#define mccompcurname  c1
#define mccompcurtype  Collimator_linear
#define mccompcurindex 18
#define slope mccc1_slope
#define slopeV mccc1_slopeV
{   /* Declarations of c1=Collimator_linear() SETTING parameters. */
MCNUM xmin = mccc1_xmin;
MCNUM xmax = mccc1_xmax;
MCNUM ymin = mccc1_ymin;
MCNUM ymax = mccc1_ymax;
MCNUM xwidth = mccc1_xwidth;
MCNUM yheight = mccc1_yheight;
MCNUM length = mccc1_length;
MCNUM divergence = mccc1_divergence;
MCNUM transmission = mccc1_transmission;
MCNUM divergenceV = mccc1_divergenceV;
#line 100 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Collimator_linear.comp"
{
  double x;
  int i;

  
  for(x = xmin, i = 0; i <= 3; i++, x += (xmax - xmin)/3.0)
    multiline(5, x, (double)ymin, 0.0, x, (double)ymax, 0.0,
              x, (double)ymax, (double)length, x, (double)ymin, (double)length,
              x, (double)ymin, 0.0);
  line(xmin, ymin, 0,   xmax, ymin, 0);
  line(xmin, ymax, 0,   xmax, ymax, 0);
  line(xmin, ymin, length, xmax, ymin, length);
  line(xmin, ymax, length, xmax, ymax, length);
}
#line 17037 "./linup-6.c"
}   /* End of c1=Collimator_linear() SETTING parameter declarations. */
#undef slopeV
#undef slope
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'slitMS3'. */
  SIG_MESSAGE("slitMS3 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "slitMS3");
#define mccompcurname  slitMS3
#define mccompcurtype  Slit
#define mccompcurindex 19
{   /* Declarations of slitMS3=Slit() SETTING parameters. */
MCNUM xmin = mccslitMS3_xmin;
MCNUM xmax = mccslitMS3_xmax;
MCNUM ymin = mccslitMS3_ymin;
MCNUM ymax = mccslitMS3_ymax;
MCNUM radius = mccslitMS3_radius;
MCNUM xwidth = mccslitMS3_xwidth;
MCNUM yheight = mccslitMS3_yheight;
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
#line 17082 "./linup-6.c"
}   /* End of slitMS3=Slit() SETTING parameter declarations. */
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'slitMS4'. */
  SIG_MESSAGE("slitMS4 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "slitMS4");
#define mccompcurname  slitMS4
#define mccompcurtype  Slit
#define mccompcurindex 20
{   /* Declarations of slitMS4=Slit() SETTING parameters. */
MCNUM xmin = mccslitMS4_xmin;
MCNUM xmax = mccslitMS4_xmax;
MCNUM ymin = mccslitMS4_ymin;
MCNUM ymax = mccslitMS4_ymax;
MCNUM radius = mccslitMS4_radius;
MCNUM xwidth = mccslitMS4_xwidth;
MCNUM yheight = mccslitMS4_yheight;
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
#line 17125 "./linup-6.c"
}   /* End of slitMS4=Slit() SETTING parameter declarations. */
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'slitMS5'. */
  SIG_MESSAGE("slitMS5 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "slitMS5");
#define mccompcurname  slitMS5
#define mccompcurtype  Slit
#define mccompcurindex 21
{   /* Declarations of slitMS5=Slit() SETTING parameters. */
MCNUM xmin = mccslitMS5_xmin;
MCNUM xmax = mccslitMS5_xmax;
MCNUM ymin = mccslitMS5_ymin;
MCNUM ymax = mccslitMS5_ymax;
MCNUM radius = mccslitMS5_radius;
MCNUM xwidth = mccslitMS5_xwidth;
MCNUM yheight = mccslitMS5_yheight;
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
#line 17168 "./linup-6.c"
}   /* End of slitMS5=Slit() SETTING parameter declarations. */
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'mon'. */
  SIG_MESSAGE("mon (McDisplay)");
  printf("MCDISPLAY: component %s\n", "mon");
#define mccompcurname  mon
#define mccompcurtype  Monitor
#define mccompcurindex 22
#define Nsum mccmon_Nsum
#define psum mccmon_psum
#define p2sum mccmon_p2sum
{   /* Declarations of mon=Monitor() SETTING parameters. */
MCNUM xmin = mccmon_xmin;
MCNUM xmax = mccmon_xmax;
MCNUM ymin = mccmon_ymin;
MCNUM ymax = mccmon_ymax;
MCNUM xwidth = mccmon_xwidth;
MCNUM yheight = mccmon_yheight;
MCNUM restore_neutron = mccmon_restore_neutron;
#line 96 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/Monitor.comp"
{
  
  multiline(5, (double)xmin, (double)ymin, 0.0,
               (double)xmax, (double)ymin, 0.0,
               (double)xmax, (double)ymax, 0.0,
               (double)xmin, (double)ymax, 0.0,
               (double)xmin, (double)ymin, 0.0);
}
#line 17200 "./linup-6.c"
}   /* End of mon=Monitor() SETTING parameter declarations. */
#undef p2sum
#undef psum
#undef Nsum
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'emon1'. */
  SIG_MESSAGE("emon1 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "emon1");
#define mccompcurname  emon1
#define mccompcurtype  E_monitor
#define mccompcurindex 23
#define nE mccemon1_nE
#define E_N mccemon1_E_N
#define E_p mccemon1_E_p
#define E_p2 mccemon1_E_p2
#define S_p mccemon1_S_p
#define S_pE mccemon1_S_pE
#define S_pE2 mccemon1_S_pE2
{   /* Declarations of emon1=E_monitor() SETTING parameters. */
char* filename = mccemon1_filename;
MCNUM xmin = mccemon1_xmin;
MCNUM xmax = mccemon1_xmax;
MCNUM ymin = mccemon1_ymin;
MCNUM ymax = mccemon1_ymax;
MCNUM xwidth = mccemon1_xwidth;
MCNUM yheight = mccemon1_yheight;
MCNUM Emin = mccemon1_Emin;
MCNUM Emax = mccemon1_Emax;
MCNUM restore_neutron = mccemon1_restore_neutron;
int nowritefile = mccemon1_nowritefile;
#line 132 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/E_monitor.comp"
{
  
  multiline(5, (double)xmin, (double)ymin, 0.0,
               (double)xmax, (double)ymin, 0.0,
               (double)xmax, (double)ymax, 0.0,
               (double)xmin, (double)ymax, 0.0,
               (double)xmin, (double)ymin, 0.0);
}
#line 17243 "./linup-6.c"
}   /* End of emon1=E_monitor() SETTING parameter declarations. */
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

  /* MCDISPLAY code for component 'sample'. */
  SIG_MESSAGE("sample (McDisplay)");
  printf("MCDISPLAY: component %s\n", "sample");
#define mccompcurname  sample
#define mccompcurtype  V_sample
#define mccompcurindex 24
#define VarsV mccsample_VarsV
{   /* Declarations of sample=V_sample() SETTING parameters. */
MCNUM radius = mccsample_radius;
MCNUM thickness = mccsample_thickness;
MCNUM zdepth = mccsample_zdepth;
MCNUM Vc = mccsample_Vc;
MCNUM sigma_abs = mccsample_sigma_abs;
MCNUM sigma_inc = mccsample_sigma_inc;
MCNUM radius_i = mccsample_radius_i;
MCNUM radius_o = mccsample_radius_o;
MCNUM h = mccsample_h;
MCNUM focus_r = mccsample_focus_r;
MCNUM pack = mccsample_pack;
MCNUM frac = mccsample_frac;
MCNUM f_QE = mccsample_f_QE;
MCNUM gamma = mccsample_gamma;
MCNUM target_x = mccsample_target_x;
MCNUM target_y = mccsample_target_y;
MCNUM target_z = mccsample_target_z;
MCNUM focus_xw = mccsample_focus_xw;
MCNUM focus_yh = mccsample_focus_yh;
MCNUM focus_aw = mccsample_focus_aw;
MCNUM focus_ah = mccsample_focus_ah;
MCNUM xwidth = mccsample_xwidth;
MCNUM yheight = mccsample_yheight;
MCNUM zthick = mccsample_zthick;
MCNUM rad_sphere = mccsample_rad_sphere;
MCNUM sig_a = mccsample_sig_a;
MCNUM sig_i = mccsample_sig_i;
MCNUM V0 = mccsample_V0;
int target_index = mccsample_target_index;
MCNUM multiples = mccsample_multiples;
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
#line 17341 "./linup-6.c"
}   /* End of sample=V_sample() SETTING parameter declarations. */
#undef VarsV
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'a3'. */
  SIG_MESSAGE("a3 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "a3");
#define mccompcurname  a3
#define mccompcurtype  Arm
#define mccompcurindex 25
#line 40 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Arm.comp"
{
  /* A bit ugly; hard-coded dimensions. */
  
  line(0,0,0,0.2,0,0);
  line(0,0,0,0,0.2,0);
  line(0,0,0,0,0,0.2);
}
#line 17362 "./linup-6.c"
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'focus_check'. */
  SIG_MESSAGE("focus_check (McDisplay)");
  printf("MCDISPLAY: component %s\n", "focus_check");
#define mccompcurname  focus_check
#define mccompcurtype  PSD_monitor
#define mccompcurindex 26
#define PSD_N mccfocus_check_PSD_N
#define PSD_p mccfocus_check_PSD_p
#define PSD_p2 mccfocus_check_PSD_p2
{   /* Declarations of focus_check=PSD_monitor() SETTING parameters. */
int nx = mccfocus_check_nx;
int ny = mccfocus_check_ny;
char* filename = mccfocus_check_filename;
MCNUM xmin = mccfocus_check_xmin;
MCNUM xmax = mccfocus_check_xmax;
MCNUM ymin = mccfocus_check_ymin;
MCNUM ymax = mccfocus_check_ymax;
MCNUM xwidth = mccfocus_check_xwidth;
MCNUM yheight = mccfocus_check_yheight;
MCNUM restore_neutron = mccfocus_check_restore_neutron;
#line 129 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
{
  multiline(5,
    (double)xmin, (double)ymin, 0.0,
    (double)xmax, (double)ymin, 0.0,
    (double)xmax, (double)ymax, 0.0,
    (double)xmin, (double)ymax, 0.0,
    (double)xmin, (double)ymin, 0.0);
}
#line 17396 "./linup-6.c"
}   /* End of focus_check=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'c2'. */
  SIG_MESSAGE("c2 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "c2");
#define mccompcurname  c2
#define mccompcurtype  Collimator_linear
#define mccompcurindex 27
#define slope mccc2_slope
#define slopeV mccc2_slopeV
{   /* Declarations of c2=Collimator_linear() SETTING parameters. */
MCNUM xmin = mccc2_xmin;
MCNUM xmax = mccc2_xmax;
MCNUM ymin = mccc2_ymin;
MCNUM ymax = mccc2_ymax;
MCNUM xwidth = mccc2_xwidth;
MCNUM yheight = mccc2_yheight;
MCNUM length = mccc2_length;
MCNUM divergence = mccc2_divergence;
MCNUM transmission = mccc2_transmission;
MCNUM divergenceV = mccc2_divergenceV;
#line 100 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Collimator_linear.comp"
{
  double x;
  int i;

  
  for(x = xmin, i = 0; i <= 3; i++, x += (xmax - xmin)/3.0)
    multiline(5, x, (double)ymin, 0.0, x, (double)ymax, 0.0,
              x, (double)ymax, (double)length, x, (double)ymin, (double)length,
              x, (double)ymin, 0.0);
  line(xmin, ymin, 0,   xmax, ymin, 0);
  line(xmin, ymax, 0,   xmax, ymax, 0);
  line(xmin, ymin, length, xmax, ymin, length);
  line(xmin, ymax, length, xmax, ymax, length);
}
#line 17439 "./linup-6.c"
}   /* End of c2=Collimator_linear() SETTING parameter declarations. */
#undef slopeV
#undef slope
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'ana'. */
  SIG_MESSAGE("ana (McDisplay)");
  printf("MCDISPLAY: component %s\n", "ana");
#define mccompcurname  ana
#define mccompcurtype  Monochromator_flat
#define mccompcurindex 28
#define mos_rms_y mccana_mos_rms_y
#define mos_rms_z mccana_mos_rms_z
#define mos_rms_max mccana_mos_rms_max
#define mono_Q mccana_mono_Q
{   /* Declarations of ana=Monochromator_flat() SETTING parameters. */
MCNUM zmin = mccana_zmin;
MCNUM zmax = mccana_zmax;
MCNUM ymin = mccana_ymin;
MCNUM ymax = mccana_ymax;
MCNUM zwidth = mccana_zwidth;
MCNUM yheight = mccana_yheight;
MCNUM mosaich = mccana_mosaich;
MCNUM mosaicv = mccana_mosaicv;
MCNUM r0 = mccana_r0;
MCNUM Q = mccana_Q;
MCNUM DM = mccana_DM;
#line 254 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Monochromator_flat.comp"
{
  
  multiline(5, 0.0, (double)ymin, (double)zmin,
               0.0, (double)ymax, (double)zmin,
               0.0, (double)ymax, (double)zmax,
               0.0, (double)ymin, (double)zmax,
               0.0, (double)ymin, (double)zmin);
}
#line 17478 "./linup-6.c"
}   /* End of ana=Monochromator_flat() SETTING parameter declarations. */
#undef mono_Q
#undef mos_rms_max
#undef mos_rms_z
#undef mos_rms_y
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'a4'. */
  SIG_MESSAGE("a4 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "a4");
#define mccompcurname  a4
#define mccompcurtype  Arm
#define mccompcurindex 29
#line 40 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Arm.comp"
{
  /* A bit ugly; hard-coded dimensions. */
  
  line(0,0,0,0.2,0,0);
  line(0,0,0,0,0.2,0);
  line(0,0,0,0,0,0.2);
}
#line 17502 "./linup-6.c"
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'c3'. */
  SIG_MESSAGE("c3 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "c3");
#define mccompcurname  c3
#define mccompcurtype  Collimator_linear
#define mccompcurindex 30
#define slope mccc3_slope
#define slopeV mccc3_slopeV
{   /* Declarations of c3=Collimator_linear() SETTING parameters. */
MCNUM xmin = mccc3_xmin;
MCNUM xmax = mccc3_xmax;
MCNUM ymin = mccc3_ymin;
MCNUM ymax = mccc3_ymax;
MCNUM xwidth = mccc3_xwidth;
MCNUM yheight = mccc3_yheight;
MCNUM length = mccc3_length;
MCNUM divergence = mccc3_divergence;
MCNUM transmission = mccc3_transmission;
MCNUM divergenceV = mccc3_divergenceV;
#line 100 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Collimator_linear.comp"
{
  double x;
  int i;

  
  for(x = xmin, i = 0; i <= 3; i++, x += (xmax - xmin)/3.0)
    multiline(5, x, (double)ymin, 0.0, x, (double)ymax, 0.0,
              x, (double)ymax, (double)length, x, (double)ymin, (double)length,
              x, (double)ymin, 0.0);
  line(xmin, ymin, 0,   xmax, ymin, 0);
  line(xmin, ymax, 0,   xmax, ymax, 0);
  line(xmin, ymin, length, xmax, ymin, length);
  line(xmin, ymax, length, xmax, ymax, length);
}
#line 17541 "./linup-6.c"
}   /* End of c3=Collimator_linear() SETTING parameter declarations. */
#undef slopeV
#undef slope
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'sng'. */
  SIG_MESSAGE("sng (McDisplay)");
  printf("MCDISPLAY: component %s\n", "sng");
#define mccompcurname  sng
#define mccompcurtype  Monitor
#define mccompcurindex 31
#define Nsum mccsng_Nsum
#define psum mccsng_psum
#define p2sum mccsng_p2sum
{   /* Declarations of sng=Monitor() SETTING parameters. */
MCNUM xmin = mccsng_xmin;
MCNUM xmax = mccsng_xmax;
MCNUM ymin = mccsng_ymin;
MCNUM ymax = mccsng_ymax;
MCNUM xwidth = mccsng_xwidth;
MCNUM yheight = mccsng_yheight;
MCNUM restore_neutron = mccsng_restore_neutron;
#line 96 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/Monitor.comp"
{
  
  multiline(5, (double)xmin, (double)ymin, 0.0,
               (double)xmax, (double)ymin, 0.0,
               (double)xmax, (double)ymax, 0.0,
               (double)xmin, (double)ymax, 0.0,
               (double)xmin, (double)ymin, 0.0);
}
#line 17575 "./linup-6.c"
}   /* End of sng=Monitor() SETTING parameter declarations. */
#undef p2sum
#undef psum
#undef Nsum
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'emon2'. */
  SIG_MESSAGE("emon2 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "emon2");
#define mccompcurname  emon2
#define mccompcurtype  E_monitor
#define mccompcurindex 32
#define nE mccemon2_nE
#define E_N mccemon2_E_N
#define E_p mccemon2_E_p
#define E_p2 mccemon2_E_p2
#define S_p mccemon2_S_p
#define S_pE mccemon2_S_pE
#define S_pE2 mccemon2_S_pE2
{   /* Declarations of emon2=E_monitor() SETTING parameters. */
char* filename = mccemon2_filename;
MCNUM xmin = mccemon2_xmin;
MCNUM xmax = mccemon2_xmax;
MCNUM ymin = mccemon2_ymin;
MCNUM ymax = mccemon2_ymax;
MCNUM xwidth = mccemon2_xwidth;
MCNUM yheight = mccemon2_yheight;
MCNUM Emin = mccemon2_Emin;
MCNUM Emax = mccemon2_Emax;
MCNUM restore_neutron = mccemon2_restore_neutron;
int nowritefile = mccemon2_nowritefile;
#line 132 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/E_monitor.comp"
{
  
  multiline(5, (double)xmin, (double)ymin, 0.0,
               (double)xmax, (double)ymin, 0.0,
               (double)xmax, (double)ymax, 0.0,
               (double)xmin, (double)ymax, 0.0,
               (double)xmin, (double)ymin, 0.0);
}
#line 17618 "./linup-6.c"
}   /* End of emon2=E_monitor() SETTING parameter declarations. */
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
/* end of generated C code ./linup-6.c */
