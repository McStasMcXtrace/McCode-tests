/* Automatically generated file. Do not edit. 
 * Format:     ANSI C source code
 * Creator:    McStas <http://www.mcstas.org>
 * Instrument: /zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr (PSI_DMC)
 * Date:       Wed Feb 26 19:13:17 2020
 * File:       ./PSI_DMC.c
 * Compile:    cc -o PSI_DMC.out ./PSI_DMC.c 
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

#line 712 "./PSI_DMC.c"

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

#line 945 "./PSI_DMC.c"

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

#line 4977 "./PSI_DMC.c"

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

#line 5337 "./PSI_DMC.c"
#ifdef MC_TRACE_ENABLED
int mctraceenabled = 1;
#else
int mctraceenabled = 0;
#endif
#define MCSTAS "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../"
int mcdefaultmain = 1;
char mcinstrument_name[] = "PSI_DMC";
char mcinstrument_source[] = "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr";
char *mcinstrument_exe=NULL; /* will be set to argv[0] in main */
int main(int argc, char *argv[]){return mccode_main(argc, argv);}
void mcinit(void);
void mcraytrace(void);
void mcsave(FILE *);
void mcfinally(void);
void mcdisplay(void);

/* Shared user declarations for all components 'Source_Maxwell_3'. */
#line 70 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../sources/Source_Maxwell_3.comp"
/* A normalised Maxwellian distribution : Integral over all l = 1 */
double SM3_Maxwell(double l, double temp)
  {
    double a=949.0/temp;
    return 2*a*a*exp(-a/(l*l))/(l*l*l*l*l);
  }
#line 5363 "./PSI_DMC.c"

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

#line 5449 "./PSI_DMC.c"

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

#line 7044 "./PSI_DMC.c"

/* Shared user declarations for all components 'Bender'. */
#line 102 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Bender.comp"

#line 7049 "./PSI_DMC.c"

/* Shared user declarations for all components 'Al_window'. */
#line 47 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../contrib/Al_window.comp"
/* ToDo: Should be component local names. */
#ifndef AL_WINDOW
#define avogadro 6.022 /* 10E23 Atoms per mole (mol-1) */
#define Al_sigma_a .231 /* Absorption cross section per atom (barns) at 2200m/s */
#define Al_sigma_i .0082 /* Incoherent scattering cross section per atom (barns) */
#define Al_rho 2.7 /* density (gcm-3) */
#define Al_mmol 27 /* molar mass Al (gmol-1) */
#define Al_my_s (Al_rho / Al_mmol * Al_sigma_i * avogadro * 10) /* inc. XS (barn) */
#define Al_my_a_v (Al_rho / Al_mmol * Al_sigma_a * avogadro * 10 * 2200 )
/* Define Constants for Polynomial Fit of
  sigma_tot(lambda)=A+B1*X+B2*X^2+B3*X^3+B4*X^4+... */
#define Al_pf_A 1.34722
#define Al_pf_B1 .12409
#define Al_pf_B2 .01078
#define Al_pf_B3 -3.25895e-5
#define Al_pf_B4 3.74731e-6
#define AL_WINDOW
#endif
#line 7071 "./PSI_DMC.c"

/* Shared user declarations for all components 'Monochromator_2foc'. */
#line 83 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../contrib/Monochromator_2foc.comp"

#line 7076 "./PSI_DMC.c"

/* Shared user declarations for all components 'PowderN'. */
#line 217 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../samples/PowderN.comp"
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
#ifndef POWDERN_DECL
#define POWDERN_DECL
/* format definitions in the order {j d F2 DW Dd inv2d q F strain} */
#ifndef Crystallographica
#define Crystallographica { 4,5,7,0,0,0,0,0,0 }
#define Fullprof          { 4,0,8,0,0,5,0,0,0 }
#define Lazy              {17,6,0,0,0,0,0,13,0 }
#define Undefined         { 0,0,0,0,0,0,0,0,0 }
#endif
  
struct line_data
{
double F2;                  /* Value of structure factor */
double q;                   /* Qvector */
int j;                      /* Multiplicity */
double DWfactor;            /* Debye-Waller factor */
double w;                   /* Intrinsic line width */
double Epsilon;             /* Strain=delta_d_d/d shift in ppm */
};

struct line_info_struct
{
struct line_data *list;     /* Reflection array */
int  count;                  /* Number of reflections */
double Dd;
double DWfactor;
double V_0;
double rho;
double at_weight;
double at_nb;
double sigma_a;
double sigma_i;
char   compname[256];
double flag_barns;
int    shape; /* 0 cylinder, 1 box, 2 sphere, 3 OFF file */
int    column_order[9]; /* column signification */
int    flag_warning;
char   type;  /* interaction type of event t=Transmit, i=Incoherent, c=Coherent */
int    itype; /* int representation of type */ 
double dq;    /* wavevector transfer [Angs-1] */
double Epsilon; /* global strain in ppm */
double XsectionFactor;
double my_s_v2_sum;
double my_a_v;
double my_inc;
double lfree; // store mean free path for the last event;
double *w_v,*q_v, *my_s_v2;
double radius_i,xwidth_i,yheight_i,zdepth_i;
double v; /* last velocity (cached) */
double Nq;
int    nb_reuses, nb_refl, nb_refl_count;
double v_min, v_max;
double xs_Nq[CHAR_BUF_LENGTH];
double xs_sum[CHAR_BUF_LENGTH];
double neutron_passed;
long   xs_compute, xs_reuse, xs_calls;
};

  off_struct offdata;

  // PN_list_compare *****************************************************************

  int PN_list_compare (void const *a, void const *b)
  {
     struct line_data const *pa = a;
     struct line_data const *pb = b;
     double s = pa->q - pb->q;

     if (!s) return 0;
     else    return (s < 0 ? -1 : 1);
  } /* PN_list_compare */

  int read_line_data(char *SC_file, struct line_info_struct *info)
  {
    struct line_data *list = NULL;
    int    size = 0;
    t_Table sTable; /* sample data table structure from SC_file */
    int    i=0;
    int    mult_count  =0;
    char   flag=0;
    double q_count=0, j_count=0, F2_count=0;
    char **parsing;
    int    list_count=0;

    if (!SC_file || !strlen(SC_file) || !strcmp(SC_file, "NULL")) {
      MPI_MASTER(
      printf("PowderN: %s: Using incoherent elastic scattering only.\n",
          info->compname);
      );
      info->count = 0;
      return(0);
    }
    Table_Read(&sTable, SC_file, 1); /* read 1st block data from SC_file into sTable*/

    /* parsing of header */
    parsing = Table_ParseHeader(sTable.header,
      "Vc","V_0",
      "sigma_abs","sigma_a ",
      "sigma_inc","sigma_i ",
      "column_j",
      "column_d",
      "column_F2",
      "column_DW",
      "column_Dd",
      "column_inv2d", "column_1/2d", "column_sintheta/lambda",
      "column_q", /* 14 */
      "DW", "Debye_Waller",
      "delta_d_d/d",
      "column_F ",
      "V_rho",
      "density",
      "weight",
      "nb_atoms","multiplicity", /* 23 */
      "column_ppm","column_strain",
      NULL);

    if (parsing) {
      if (parsing[0] && !info->V_0)     info->V_0    =atof(parsing[0]);
      if (parsing[1] && !info->V_0)     info->V_0    =atof(parsing[1]);
      if (parsing[2] && !info->sigma_a) info->sigma_a=atof(parsing[2]);
      if (parsing[3] && !info->sigma_a) info->sigma_a=atof(parsing[3]);
      if (parsing[4] && !info->sigma_i) info->sigma_i=atof(parsing[4]);
      if (parsing[5] && !info->sigma_i) info->sigma_i=atof(parsing[5]);
      if (parsing[6])                   info->column_order[0]=atoi(parsing[6]);
      if (parsing[7])                   info->column_order[1]=atoi(parsing[7]);
      if (parsing[8])                   info->column_order[2]=atoi(parsing[8]);
      if (parsing[9])                   info->column_order[3]=atoi(parsing[9]);
      if (parsing[10])                  info->column_order[4]=atoi(parsing[10]);
      if (parsing[11])                  info->column_order[5]=atoi(parsing[11]);
      if (parsing[12])                  info->column_order[5]=atoi(parsing[12]);
      if (parsing[13])                  info->column_order[5]=atoi(parsing[13]);
      if (parsing[14])                  info->column_order[6]=atoi(parsing[14]);
      if (parsing[15] && info->DWfactor<=0)    info->DWfactor=atof(parsing[15]);
      if (parsing[16] && info->DWfactor<=0)    info->DWfactor=atof(parsing[16]);
      if (parsing[17] && info->Dd <0)          info->Dd      =atof(parsing[17]);
      if (parsing[18])                  info->column_order[7]=atoi(parsing[18]);
      if (parsing[19] && !info->V_0)    info->V_0    =1/atof(parsing[19]);
      if (parsing[20] && !info->rho)    info->rho    =atof(parsing[20]);
      if (parsing[21] && !info->at_weight)     info->at_weight    =atof(parsing[21]);
      if (parsing[22] && info->at_nb <= 1)  info->at_nb    =atof(parsing[22]);
      if (parsing[23] && info->at_nb <= 1)  info->at_nb    =atof(parsing[23]);
      if (parsing[24])                  info->column_order[8]=atoi(parsing[24]);
      if (parsing[25])                  info->column_order[8]=atoi(parsing[25]);
      for (i=0; i<=25; i++) if (parsing[i]) free(parsing[i]);
      free(parsing);
    }

    if (!sTable.rows)
      exit(fprintf(stderr, "PowderN: %s: Error: The number of rows in %s "
       "should be at least %d\n", info->compname, SC_file, 1));
    else
      size = sTable.rows;

    MPI_MASTER(
    Table_Info(sTable);
    printf("PowderN: %s: Reading %d rows from %s\n",
          info->compname, size, SC_file);
    );

    if (info->column_order[0] == 4 && info->flag_barns !=0)
    MPI_MASTER(
      printf("PowderN: %s: Powder file probably of type Crystallographica/Fullprof (lau)\n"
           "WARNING: but F2 unit is set to barns=1 (barns). Intensity might be 100 times too high.\n",
           info->compname);
    );
    if (info->column_order[0] == 17 && info->flag_barns == 0)
    MPI_MASTER(
      printf("PowderN: %s: Powder file probably of type Lazy Pulver (laz)\n"
           "WARNING: but F2 unit is set to barns=0 (fm^2). Intensity might be 100 times too low.\n",
           info->compname);
    );
    /* allocate line_data array */
    list = (struct line_data*)malloc(size*sizeof(struct line_data));

    for (i=0; i<size; i++)
    {
      /*      printf("Reading in line %i\n",i);*/
      double j=0, d=0, w=0, q=0, DWfactor=0, F2=0, Epsilon=0;
      int index;

      if (info->Dd >= 0)      w         = info->Dd;
      if (info->DWfactor > 0) DWfactor  = info->DWfactor;
      if (info->Epsilon)      Epsilon   = info->Epsilon*1e-6;

      /* get data from table using columns {j d F2 DW Dd inv2d q F} */
      /* column indexes start at 1, thus need to substract 1 */
      if (info->column_order[0] >0)
        j = Table_Index(sTable, i, info->column_order[0]-1);
      if (info->column_order[1] >0)
        d = Table_Index(sTable, i, info->column_order[1]-1);
      if (info->column_order[2] >0)
        F2 = Table_Index(sTable, i, info->column_order[2]-1);
      if (info->column_order[3] >0)
        DWfactor = Table_Index(sTable, i, info->column_order[3]-1);
      if (info->column_order[4] >0)
        w = Table_Index(sTable, i, info->column_order[4]-1);
      if (info->column_order[5] >0)
        { d = Table_Index(sTable, i, info->column_order[5]-1);
          d = (d > 0? 1/d/2 : 0); }
      if (info->column_order[6] >0)
        { q = Table_Index(sTable, i, info->column_order[6]-1);
          d = (q > 0 ? 2*PI/q : 0); }
      if (info->column_order[7] >0  && !F2)
        { F2 = Table_Index(sTable, i, info->column_order[7]-1); F2 *= F2; }
      if (info->column_order[8] >0  && !Epsilon)
        { Epsilon = Table_Index(sTable, i, info->column_order[8]-1)*1e-6; }

      /* assign and check values */
      j        = (j > 0 ? j : 0);
      q        = (d > 0 ? 2*PI/d : 0); /* this is q */
      if (Epsilon && fabs(Epsilon) < 1e6) {
        q     -= Epsilon*q; /* dq/q = -delta_d_d/d = -Epsilon */
      }
      DWfactor = (DWfactor > 0 ? DWfactor : 1);
      w        = (w>0 ? w : 0); /* this is q and d relative spreading */
      F2       = (F2 >= 0 ? F2 : 0);
      if (j == 0 || q == 0) {
        MPI_MASTER(
        printf("PowderN: %s: line %i has invalid definition\n"
               "         (mult=0 or q=0 or d=0)\n", info->compname, i);
        );
        continue;
      }
      list[list_count].j = j;
      list[list_count].q = q;
      list[list_count].DWfactor = DWfactor;
      list[list_count].w = w;
      list[list_count].F2= F2;
      list[list_count].Epsilon = Epsilon;

      /* adjust multiplicity if j-column + multiple d-spacing lines */
      /* if  d = previous d, increase line duplication index */
      if (!q_count)      q_count  = q;
      if (!j_count)      j_count  = j;
      if (!F2_count)     F2_count = F2;
      if (fabs(q_count-q) < 0.0001*fabs(q)
       && fabs(F2_count-F2) < 0.0001*fabs(F2) && j_count == j) {
       mult_count++; flag=0; }
      else flag=1;
      if (i == size-1) flag=1;
      /* else if d != previous d : just passed equivalent lines */
      if (flag) {
        if (i == size-1) list_count++;
      /*   if duplication index == previous multiplicity */
      /*      set back multiplicity of previous lines to 1 */
        if ((mult_count && list_count>0)
            && (mult_count == list[list_count-1].j
                || ((list_count < size) && (i == size - 1)
                    && (mult_count == list[list_count].j))) ) {
          MPI_MASTER(
          printf("PowderN: %s: Set multiplicity to 1 for lines [%i:%i]\n"
                  "         (d-spacing %g is duplicated %i times)\n",
            info->compname, list_count-mult_count, list_count-1, list[list_count-1].q, mult_count);
          );
          for (index=list_count-mult_count; index<list_count; list[index++].j = 1);
          mult_count = 1;
          q_count   = q;
          j_count   = j;
          F2_count  = F2;
        }
        if (i == size-1) list_count--;
        flag=0;
      }
      list_count++;
    } /* end for */

    Table_Free(&sTable);

    /* sort the list with increasing q */
    qsort(list, list_count, sizeof(struct line_data),  PN_list_compare);

    MPI_MASTER(
    printf("PowderN: %s: Read %i reflections from file '%s'\n",
      info->compname, list_count, SC_file);
    );

    info->list  = list;
    info->count = list_count;

    return(list_count);
  } /* read_line_data */


/* computes the number of possible reflections (return value), and the total xsection 'sum' */
/* this routine looks for a pre-computed value in the Nq and sum cache tables               */
/* when found, the earch starts from the corresponding lower element in the table           */
int calc_xsect(double v, double *qv, double *my_sv2, int count, double *sum,
  struct line_info_struct *line_info) {
  int    Nq = 0, line=0, line0=0;
  *sum=0;

  /* check if a line_info element has been recorded already */
  if (v >= line_info->v_min && v <= line_info->v_max && line_info->neutron_passed >= CHAR_BUF_LENGTH) {
    line = (int)floor(v - line_info->v_min)*CHAR_BUF_LENGTH/(line_info->v_max - line_info->v_min);
    Nq    = line_info->xs_Nq[line];
    *sum  = line_info->xs_sum[line];
    if (!Nq && *sum == 0) {
      /* not yet set: we compute the sum up to the corresponding speed in the table cache */
      double line_v = line_info->v_min + line*(line_info->v_max - line_info->v_min)/CHAR_BUF_LENGTH;
      for(line0=0; line0<count; line0++) {
        if (qv[line0] <= 2*line_v) { /* q < 2*kf: restrict structural range */
          *sum += my_sv2[line0];
          if (Nq < line0+1) Nq=line0+1; /* determine maximum line index which can scatter */
        } else break;
      }
      line_info->xs_Nq[line] = Nq;
      line_info->xs_sum[line]= *sum;
      line_info->xs_compute++;
    } else line_info->xs_reuse++;
    line0 = Nq;
  }

  line_info->xs_calls++;

  for(line=line0; line<count; line++) {
    if (qv[line] <= 2*v) { /* q < 2*kf: restrict structural range */
      *sum += my_sv2[line];
      if (Nq < line+1) Nq=line+1; /* determine maximum line index which can scatter */
    } else break;
  }

  return(Nq);
} /* calc_xsect */

#endif /* !POWDERN_DECL */

#line 8401 "./PSI_DMC.c"

/* Shared user declarations for all components 'Monitor_nD'. */
#line 216 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/Monitor_nD.comp"
/*******************************************************************************
*
* McStas, neutron ray-tracing package
*         Copyright 1997-2002, All rights reserved
*         Risoe National Laboratory, Roskilde, Denmark
*         Institut Laue Langevin, Grenoble, France
*
* Library: share/monitor_nd-lib.h
*
* %Identification
* Written by: EF
* Date: Aug 28, 2002
* Origin: ILL
* Release: McStas 1.6
* Version: $Revision$
*
* This file is to be imported by the monitor_nd related components
* It handles some shared functions.
*
* Usage: within SHARE
* %include "monitor_nd-lib"
*
*******************************************************************************/

#ifndef MONITOR_ND_LIB_H

#define MONITOR_ND_LIB_H "$Revision$"
#define MONnD_COORD_NMAX  30  /* max number of variables to record */

  typedef struct MonitornD_Defines
  {
    int COORD_NONE  ;
    int COORD_X     ;
    int COORD_Y     ;
    int COORD_Z     ;
    int COORD_RADIUS; 
    int COORD_VX    ;
    int COORD_VY    ;
    int COORD_VZ    ;
    int COORD_V     ;
    int COORD_T     ;
    int COORD_P     ;
    int COORD_SX    ;
    int COORD_SY    ;
    int COORD_SZ    ;
    int COORD_KX    ;
    int COORD_KY    ;
    int COORD_KZ    ;
    int COORD_K     ;
    int COORD_ENERGY;
    int COORD_LAMBDA;
    int COORD_KXY   ;
    int COORD_KYZ   ;
    int COORD_KXZ   ;
    int COORD_VXY   ;
    int COORD_VYZ   ;
    int COORD_VXZ   ;
    int COORD_HDIV  ;
    int COORD_VDIV  ;
    int COORD_ANGLE ;
    int COORD_NCOUNT;
    int COORD_THETA ;
    int COORD_PHI   ;
    int COORD_USER1 ;
    int COORD_USER2 ;
    int COORD_USER3 ;
    int COORD_XY    ;
    int COORD_XZ    ;
    int COORD_YZ    ;
    int COORD_PIXELID;

    /* token modifiers */
    int COORD_VAR   ; /* next token should be a variable or normal option */
    int COORD_MIN   ; /* next token is a min value */
    int COORD_MAX   ; /* next token is a max value */
    int COORD_DIM   ; /* next token is a bin value */
    int COORD_FIL   ; /* next token is a filename */
    int COORD_EVNT  ; /* next token is a buffer size value */
    int COORD_3HE   ; /* next token is a 3He pressure value */
    int COORD_LOG   ; /* next variable will be in log scale */
    int COORD_ABS   ; /* next variable will be in abs scale */
    int COORD_SIGNAL; /* next variable will be the signal var */
    int COORD_AUTO  ; /* set auto limits */

    char TOKEN_DEL[32]; /* token separators */

    char SHAPE_SQUARE; /* shape of the monitor */
    char SHAPE_DISK  ;
    char SHAPE_SPHERE;
    char SHAPE_CYLIND;
    char SHAPE_BANANA; /* cylinder without top/bottom, on restricted angular area */
    char SHAPE_BOX   ;
    char SHAPE_PREVIOUS;
    char SHAPE_OFF;

  } MonitornD_Defines_type;

  typedef struct MonitornD_Variables
  {
    double area;
    double Sphere_Radius     ;
    double Cylinder_Height   ;
    char   Flag_With_Borders ;   /* 2 means xy borders too */
    char   Flag_List         ;   /* 1 store 1 buffer, 2 is list all, 3 list all+append */
    char   Flag_Multiple     ;   /* 1 when n1D, 0 for 2D */
    char   Flag_Verbose      ;
    int    Flag_Shape        ;
    char   Flag_Auto_Limits  ;   /* get limits from first Buffer */
    char   Flag_Absorb       ;   /* monitor is also a slit */
    char   Flag_per_cm2      ;   /* flux is per cm2 */
    char   Flag_log          ;   /* log10 of the flux */
    char   Flag_parallel     ;   /* set neutron state back after detection (parallel components) */
    char   Flag_Binary_List  ;
    char   Flag_capture      ;   /* lambda monitor with lambda/lambda(2200m/s = 1.7985 Angs) weightening */
    int    Flag_signal       ;   /* 0:monitor p, else monitor a mean value */
    int    Flag_mantid       ;   /* 0:normal monitor, else do mantid-event specifics */
    int    Flag_OFF          ;   /* Flag to indicate external geometry from OFF file */
    unsigned long OFF_polyidx;   /* When intersection is done externally by off_intersect, this gives the 
				    polygon number, i.e. pixel index */

    unsigned long Coord_Number      ;   /* total number of variables to monitor, plus intensity (0) */
    unsigned long Coord_NumberNoPixel;  /* same but without counting PixelID */
    unsigned long Buffer_Block      ;   /* Buffer size for list or auto limits */
    unsigned long Neutron_Counter   ;   /* event counter, simulation total counts is mcget_ncount() */
    unsigned long Buffer_Counter    ;   /* index in Buffer size (for realloc) */
    unsigned long Buffer_Size       ;
    int    Coord_Type[MONnD_COORD_NMAX];      /* type of variable */
    char   Coord_Label[MONnD_COORD_NMAX][30]; /* label of variable */
    char   Coord_Var[MONnD_COORD_NMAX][30];   /* short id of variable */
    long   Coord_Bin[MONnD_COORD_NMAX];       /* bins of variable array */
    long   Coord_BinProd[MONnD_COORD_NMAX];   /* product of bins of variable array */
    double Coord_Min[MONnD_COORD_NMAX];
    double Coord_Max[MONnD_COORD_NMAX];
    char   Monitor_Label[MONnD_COORD_NMAX*30];/* Label for monitor */
    char   Mon_File[128];                     /* output file name */

    double cx,cy,cz;
    double cvx, cvy, cvz;
    double ckx, cky, ckz;
    double csx, csy, csz;
    double cEx, cEy, cEz;
    double cs1, cs2, ct, cphi, cp;
    double He3_pressure;
    char   Flag_UsePreMonitor    ;   /* use a previously stored neutron parameter set */
    char   UserName1[128];
    char   UserName2[128];
    char   UserName3[128];
    double UserVariable1;
    double UserVariable2;
    double UserVariable3;
    char   option[CHAR_BUF_LENGTH];

    long long int Nsum;
    double psum, p2sum;
    double **Mon2D_N;
    double **Mon2D_p;
    double **Mon2D_p2;
    double *Mon2D_Buffer;
    unsigned long PixelID;

    double mxmin,mxmax,mymin,mymax,mzmin,mzmax;
    double mean_dx, mean_dy, min_x, min_y, max_x, max_y, mean_p;

    char   compcurname[128];
    Coords compcurpos;

  } MonitornD_Variables_type;

/* monitor_nd-lib function prototypes */
/* ========================================================================= */

void Monitor_nD_Init(MonitornD_Defines_type *, MonitornD_Variables_type *, MCNUM, MCNUM, MCNUM, MCNUM, MCNUM, MCNUM, MCNUM, MCNUM, MCNUM, int);
int Monitor_nD_Trace(MonitornD_Defines_type *, MonitornD_Variables_type *);
MCDETECTOR Monitor_nD_Save(MonitornD_Defines_type *, MonitornD_Variables_type *);
void Monitor_nD_Finally(MonitornD_Defines_type *, MonitornD_Variables_type *);
void Monitor_nD_McDisplay(MonitornD_Defines_type *,
 MonitornD_Variables_type *);

#define MONND_DECLARE(monname) \
  struct MonitornD_Variables *mcmonnd ## monname;
#define MONND_USER_TITLE(monname, num, title) \
  { mcmonnd ## monname = &(MC_GETPAR(monname, Vars)); \
    strcpy(mcmonnd ## monname->UserName ## num, title); }
#define MONND_USER_VALUE(monname, num, value) \
  { mcmonnd ## monname = &(MC_GETPAR(monname, Vars)); \
    mcmonnd ## monname->UserVariable ## num = (value); }

#endif

/* end of monitor_nd-lib.h */
/*******************************************************************************
*
* McStas, neutron ray-tracing package
*         Copyright 1997-2002, All rights reserved
*         Risoe National Laboratory, Roskilde, Denmark
*         Institut Laue Langevin, Grenoble, France
*
* Library: share/monitor_nd-lib.c
*
* %Identification
* Written by: EF
* Date: Aug 28, 2002
* Origin: ILL
* Release: McStas 1.6
* Version: $Revision$
*
* This file is to be imported by the monitor_nd related components
* It handles some shared functions. Embedded within instrument in runtime mode.
*
* Usage: within SHARE
* %include "monitor_nd-lib"
*
*******************************************************************************/

#ifndef MONITOR_ND_LIB_H
#error McStas : please import this library with %include "monitor_nd-lib"
#endif

/* ========================================================================= */
/* Monitor_nD_Init: this routine is used to parse options                    */
/* ========================================================================= */

void Monitor_nD_Init(MonitornD_Defines_type *DEFS,
  MonitornD_Variables_type *Vars,
  MCNUM xwidth,
  MCNUM yheight,
  MCNUM zdepth,
  MCNUM xmin,
  MCNUM xmax,
  MCNUM ymin,
  MCNUM ymax,
  MCNUM zmin,
  MCNUM zmax,
  int offflag)
  {
    long carg = 1;
    char *option_copy, *token;
    char Flag_New_token = 1;
    char Flag_End       = 1;
    char Flag_All       = 0;
    char Flag_No        = 0;
    char Flag_abs       = 0;
    int  Flag_auto      = 0;  /* -1: all, 1: the current variable */
    int  Set_Vars_Coord_Type;
    char Set_Vars_Coord_Label[64];
    char Set_Vars_Coord_Var[64];
    char Short_Label[MONnD_COORD_NMAX][64];
    int  Set_Coord_Mode;
    long i=0, j=0;
    double lmin, lmax, XY=0;
    long t;


    t = (long)time(NULL);

/* initialize DEFS */
/* Variables to monitor */
    DEFS->COORD_NONE   =0;
    DEFS->COORD_X      =1;
    DEFS->COORD_Y      =2;
    DEFS->COORD_Z      =3;
    DEFS->COORD_RADIUS =19;
    DEFS->COORD_VX     =4;
    DEFS->COORD_VY     =5;
    DEFS->COORD_VZ     =6;
    DEFS->COORD_V      =16;
    DEFS->COORD_T      =7;
    DEFS->COORD_P      =8;
    DEFS->COORD_SX     =9;
    DEFS->COORD_SY     =10;
    DEFS->COORD_SZ     =11;
    DEFS->COORD_KX     =12;
    DEFS->COORD_KY     =13;
    DEFS->COORD_KZ     =14;
    DEFS->COORD_K      =15;
    DEFS->COORD_ENERGY =17;
    DEFS->COORD_LAMBDA =18;
    DEFS->COORD_HDIV   =20;
    DEFS->COORD_VDIV   =21;
    DEFS->COORD_ANGLE  =22;
    DEFS->COORD_NCOUNT =23;
    DEFS->COORD_THETA  =24;
    DEFS->COORD_PHI    =25;
    DEFS->COORD_USER1  =26;
    DEFS->COORD_USER2  =27;
    DEFS->COORD_USER3  =28;
    DEFS->COORD_XY     =37;
    DEFS->COORD_YZ     =31;
    DEFS->COORD_XZ     =32;
    DEFS->COORD_VXY    =30;
    DEFS->COORD_VYZ    =34;
    DEFS->COORD_VXZ    =36;
    DEFS->COORD_KXY    =29;
    DEFS->COORD_KYZ    =33;
    DEFS->COORD_KXZ    =35;
    DEFS->COORD_PIXELID=38;

/* token modifiers */
    DEFS->COORD_VAR    =0;    /* next token should be a variable or normal option */
    DEFS->COORD_MIN    =1;    /* next token is a min value */
    DEFS->COORD_MAX    =2;    /* next token is a max value */
    DEFS->COORD_DIM    =3;    /* next token is a bin value */
    DEFS->COORD_FIL    =4;    /* next token is a filename */
    DEFS->COORD_EVNT   =5;    /* next token is a buffer size value */
    DEFS->COORD_3HE    =6;    /* next token is a 3He pressure value */
    DEFS->COORD_LOG    =64;   /* next variable will be in log scale */
    DEFS->COORD_ABS    =128;  /* next variable will be in abs scale */
    DEFS->COORD_SIGNAL =256;  /* next variable will be the signal var */
    DEFS->COORD_AUTO   =512;  /* set auto limits */

    strcpy(DEFS->TOKEN_DEL, " =,;[](){}:");  /* token separators */

    DEFS->SHAPE_SQUARE =0;    /* shape of the monitor */
    DEFS->SHAPE_DISK   =1;
    DEFS->SHAPE_SPHERE =2;
    DEFS->SHAPE_CYLIND =3;
    DEFS->SHAPE_BANANA =4;
    DEFS->SHAPE_BOX    =5;
    DEFS->SHAPE_PREVIOUS=6;
    DEFS->SHAPE_OFF=7;

    Vars->Sphere_Radius     = 0;
    Vars->Cylinder_Height   = 0;
    Vars->Flag_With_Borders = 0;   /* 2 means xy borders too */
    Vars->Flag_List         = 0;   /* 1=store 1 buffer, 2=list all, 3=re-use buffer */
    Vars->Flag_Multiple     = 0;   /* 1 when n1D, 0 for 2D */
    Vars->Flag_Verbose      = 0;
    Vars->Flag_Shape        = DEFS->SHAPE_SQUARE;
    Vars->Flag_Auto_Limits  = 0;   /* get limits from first Buffer */
    Vars->Flag_Absorb       = 0;   /* monitor is also a slit */
    Vars->Flag_per_cm2      = 0;   /* flux is per cm2 */
    Vars->Flag_log          = 0;   /* log10 of the flux */
    Vars->Flag_parallel     = 0;   /* set neutron state back after detection (parallel components) */
    Vars->Flag_Binary_List  = 0;   /* save list as a binary file (smaller) */
    Vars->Coord_Number      = 0;   /* total number of variables to monitor, plus intensity (0) */
    Vars->Coord_NumberNoPixel=0;   /* same but without counting PixelID */

/* Allow to specify size of Monitor_nD buffer via a define*/
#ifndef MONND_BUFSIZ
    Vars->Buffer_Block      = 100000;     /* Buffer size for list or auto limits */
#else
	Vars->Buffer_Block      = MONND_BUFSIZ;     /* Buffer size for list or auto limits */	
#endif
    Vars->Neutron_Counter   = 0;   /* event counter, simulation total counts is mcget_ncount() */
    Vars->Buffer_Counter    = 0;   /* index in Buffer size (for realloc) */
    Vars->Buffer_Size       = 0;
    Vars->UserVariable1     = 0;
    Vars->UserVariable2     = 0;
    Vars->He3_pressure      = 0;
    Vars->Flag_capture      = 0;
    Vars->Flag_signal       = DEFS->COORD_P;
    Vars->Flag_mantid       = 0;
    Vars->Flag_OFF          = offflag;
    Vars->OFF_polyidx       = -1;
    Vars->mean_dx=Vars->mean_dy=0;
    Vars->min_x = Vars->max_x  =0;
    Vars->min_y = Vars->max_y  =0;

    Set_Vars_Coord_Type = DEFS->COORD_NONE;
    Set_Coord_Mode = DEFS->COORD_VAR;

    /* handle size parameters */
    /* normal use is with xwidth, yheight, zdepth */
    /* if xmin,xmax,ymin,ymax,zmin,zmax are non 0, use them */
    if (fabs(xmin-xmax) == 0)
      { Vars->mxmin = -fabs(xwidth)/2; Vars->mxmax = fabs(xwidth)/2; }
    else
      { if (xmin < xmax) {Vars->mxmin = xmin; Vars->mxmax = xmax;}
        else {Vars->mxmin = xmax; Vars->mxmax = xmin;}
      }
    if (fabs(ymin-ymax) == 0)
      { Vars->mymin = -fabs(yheight)/2; Vars->mymax = fabs(yheight)/2; }
    else
      { if (ymin < ymax) {Vars->mymin = ymin; Vars->mymax = ymax;}
        else {Vars->mymin = ymax; Vars->mymax = ymin;}
      }
    if (fabs(zmin-zmax) == 0)
      { Vars->mzmin = -fabs(zdepth)/2; Vars->mzmax = fabs(zdepth)/2; }
    else
      { if (zmin < zmax) {Vars->mzmin = zmin; Vars->mzmax = zmax; }
        else {Vars->mzmin = zmax; Vars->mzmax = zmin; }
      }

    if (fabs(Vars->mzmax-Vars->mzmin) == 0)
      Vars->Flag_Shape        = DEFS->SHAPE_SQUARE;
    else
      Vars->Flag_Shape        = DEFS->SHAPE_BOX;

    if (Vars->Flag_OFF) {
      Vars->Flag_Shape        = DEFS->SHAPE_OFF;
    }
    
    /* parse option string */

    option_copy = (char*)malloc(strlen(Vars->option)+1);
    if (option_copy == NULL)
    {
      fprintf(stderr,"Monitor_nD: %s cannot allocate 'options' copy (%li). Fatal.\n", Vars->compcurname, (long)strlen(Vars->option));
      exit(-1);
    }

    if (strlen(Vars->option))
    {
      Flag_End = 0;
      strcpy(option_copy, Vars->option);
    }

    if (strstr(Vars->option, "cm2") || strstr(Vars->option, "cm^2")) Vars->Flag_per_cm2 = 1;

    if (strstr(Vars->option, "binary") || strstr(Vars->option, "float"))
      Vars->Flag_Binary_List  = 1;
    if (strstr(Vars->option, "double"))
      Vars->Flag_Binary_List  = 2;

    strcpy(Vars->Coord_Label[0],"Intensity");
    strncpy(Vars->Coord_Var[0],"p",30);
    Vars->Coord_Type[0] = DEFS->COORD_P;
    Vars->Coord_Bin[0] = 1;
    Vars->Coord_Min[0] = 0;
    Vars->Coord_Max[0] = FLT_MAX;

    /* default file name is comp_name+dateID */
    sprintf(Vars->Mon_File, "%s_%li", Vars->compcurname, t);

    carg = 1;
    while((Flag_End == 0) && (carg < 128))
    {

      if (Flag_New_token) /* retain previous token or get a new one */
      {
        if (carg == 1) token=(char *)strtok(option_copy,DEFS->TOKEN_DEL);
        else token=(char *)strtok(NULL,DEFS->TOKEN_DEL);
        if (token == NULL) Flag_End=1;
      }
      Flag_New_token = 1;
      if ((token != NULL) && (strlen(token) != 0))
      {
        char iskeyword=0; /* left at 0 when variables are processed, 1 for modifiers */
        int  old_Mode;
        /* change token to lower case */
        for (i=0; i<strlen(token); i++) token[i]=tolower(token[i]);
        /* first handle option values from preceeding keyword token detected */
        old_Mode = Set_Coord_Mode;
        if (Set_Coord_Mode == DEFS->COORD_MAX)  /* max=%i */
        {
          if (!Flag_All)
            Vars->Coord_Max[Vars->Coord_Number] = atof(token);
          else
            for (i = 0; i <= Vars->Coord_Number; Vars->Coord_Max[i++] = atof(token));
          Set_Coord_Mode = DEFS->COORD_VAR; Flag_All = 0;
        }
        if (Set_Coord_Mode == DEFS->COORD_MIN)  /* min=%i */
        {
          if (!Flag_All)
            Vars->Coord_Min[Vars->Coord_Number] = atof(token);
          else
            for (i = 0; i <= Vars->Coord_Number; Vars->Coord_Min[i++] = atof(token));
          Set_Coord_Mode = DEFS->COORD_MAX;
        }
        if (Set_Coord_Mode == DEFS->COORD_DIM)  /* bins=%i */
        {
          if (!Flag_All)
            Vars->Coord_Bin[Vars->Coord_Number] = atoi(token);
          else
            for (i = 0; i <= Vars->Coord_Number; Vars->Coord_Bin[i++] = atoi(token));
          Set_Coord_Mode = DEFS->COORD_VAR; Flag_All = 0;
        }
        if (Set_Coord_Mode == DEFS->COORD_FIL)  /* file=%s */
        {
          if (!Flag_No) strncpy(Vars->Mon_File,token,128);
          else { strcpy(Vars->Mon_File,""); Vars->Coord_Number = 0; Flag_End = 1;}
          Set_Coord_Mode = DEFS->COORD_VAR;
        }
        if (Set_Coord_Mode == DEFS->COORD_EVNT) /* list=%i */
        {
          if (!strcmp(token, "all") || Flag_All) Vars->Flag_List = 2;
          else { i = (long)ceil(atof(token)); if (i) Vars->Buffer_Block = i;
            Vars->Flag_List = 1; }
          Set_Coord_Mode = DEFS->COORD_VAR; Flag_All = 0;
        }
        if (Set_Coord_Mode == DEFS->COORD_3HE)  /* pressure=%g */
        {
            Vars->He3_pressure = atof(token);
            Set_Coord_Mode = DEFS->COORD_VAR; Flag_All = 0;
        }

        /* now look for general option keywords */
        if (!strcmp(token, "borders"))  {Vars->Flag_With_Borders = 1; iskeyword=1; }
        if (!strcmp(token, "verbose"))  {Vars->Flag_Verbose      = 1; iskeyword=1; }
        if (!strcmp(token, "log"))      {Vars->Flag_log          = 1; iskeyword=1; }
        if (!strcmp(token, "abs"))      {Flag_abs                = 1; iskeyword=1; }
        if (!strcmp(token, "multiple")) {Vars->Flag_Multiple     = 1; iskeyword=1; }
        if (!strcmp(token, "list") || !strcmp(token, "events")) {
          Vars->Flag_List = 1; Set_Coord_Mode = DEFS->COORD_EVNT;  }
        if (!strcmp(token, "limits") || !strcmp(token, "min"))
          Set_Coord_Mode = DEFS->COORD_MIN;
        if (!strcmp(token, "slit") || !strcmp(token, "absorb")) {
          Vars->Flag_Absorb = 1;  iskeyword=1; }
        if (!strcmp(token, "max"))  Set_Coord_Mode = DEFS->COORD_MAX;
        if (!strcmp(token, "bins") || !strcmp(token, "dim")) Set_Coord_Mode = DEFS->COORD_DIM;
        if (!strcmp(token, "file") || !strcmp(token, "filename")) {
          Set_Coord_Mode = DEFS->COORD_FIL;
          if (Flag_No) { strcpy(Vars->Mon_File,""); Vars->Coord_Number = 0; Flag_End = 1; }
        }
        if (!strcmp(token, "unactivate")) {
          Flag_End = 1; Vars->Coord_Number = 0; iskeyword=1; }
        if (!strcmp(token, "all"))    { Flag_All = 1;  iskeyword=1; }
        if (!strcmp(token, "sphere")) { Vars->Flag_Shape = DEFS->SHAPE_SPHERE; iskeyword=1; }
        if (!strcmp(token, "cylinder")) { Vars->Flag_Shape = DEFS->SHAPE_CYLIND; iskeyword=1; }
        if (!strcmp(token, "banana")) { Vars->Flag_Shape = DEFS->SHAPE_BANANA; iskeyword=1; }
        if (!strcmp(token, "square")) { Vars->Flag_Shape = DEFS->SHAPE_SQUARE; iskeyword=1; }
        if (!strcmp(token, "disk"))   { Vars->Flag_Shape = DEFS->SHAPE_DISK; iskeyword=1; }
        if (!strcmp(token, "box"))     { Vars->Flag_Shape = DEFS->SHAPE_BOX; iskeyword=1; }
        if (!strcmp(token, "previous")) { Vars->Flag_Shape = DEFS->SHAPE_PREVIOUS; iskeyword=1; }
        if (!strcmp(token, "parallel")){ Vars->Flag_parallel = 1; iskeyword=1; }
        if (!strcmp(token, "capture")) { Vars->Flag_capture = 1; iskeyword=1; }
        if (!strcmp(token, "auto") && (Flag_auto != -1)) {
          Vars->Flag_Auto_Limits = 1;
          if (Flag_All) Flag_auto = -1;
          else          Flag_auto = 1;
          iskeyword=1; Flag_All=0; }
        if (!strcmp(token, "premonitor")) {
          Vars->Flag_UsePreMonitor = 1; iskeyword=1; }
        if (!strcmp(token, "3He_pressure") || !strcmp(token, "pressure")) {
          Vars->He3_pressure = 3; iskeyword=1; }
        if (!strcmp(token, "no") || !strcmp(token, "not")) { Flag_No = 1;  iskeyword=1; }
        if (!strcmp(token, "signal")) Set_Coord_Mode = DEFS->COORD_SIGNAL;
        if (!strcmp(token, "mantid")) { Vars->Flag_mantid = 1; iskeyword=1; }

        /* Mode has changed: this was a keyword or value  ? */
        if (Set_Coord_Mode != old_Mode) iskeyword=1;

        /* now look for variable names to monitor */
        Set_Vars_Coord_Type = DEFS->COORD_NONE; lmin = 0; lmax = 0;

        if (!strcmp(token, "x"))
          { Set_Vars_Coord_Type = DEFS->COORD_X; strcpy(Set_Vars_Coord_Label,"x [m]"); strcpy(Set_Vars_Coord_Var,"x");
          lmin = Vars->mxmin; lmax = Vars->mxmax;
          Vars->Coord_Min[Vars->Coord_Number+1] = Vars->mxmin;
          Vars->Coord_Max[Vars->Coord_Number+1] = Vars->mxmax;}
        if (!strcmp(token, "y"))
          { Set_Vars_Coord_Type = DEFS->COORD_Y; strcpy(Set_Vars_Coord_Label,"y [m]"); strcpy(Set_Vars_Coord_Var,"y");
          lmin = Vars->mymin; lmax = Vars->mymax;
          Vars->Coord_Min[Vars->Coord_Number+1] = Vars->mymin;
          Vars->Coord_Max[Vars->Coord_Number+1] = Vars->mymax;}
        if (!strcmp(token, "z"))
          { Set_Vars_Coord_Type = DEFS->COORD_Z; strcpy(Set_Vars_Coord_Label,"z [m]"); strcpy(Set_Vars_Coord_Var,"z"); lmin = Vars->mzmin; lmax = Vars->mzmax; }
        if (!strcmp(token, "k") || !strcmp(token, "wavevector"))
          { Set_Vars_Coord_Type = DEFS->COORD_K; strcpy(Set_Vars_Coord_Label,"|k| [Angs-1]"); strcpy(Set_Vars_Coord_Var,"k"); lmin = 0; lmax = 10; }
        if (!strcmp(token, "v"))
          { Set_Vars_Coord_Type = DEFS->COORD_V; strcpy(Set_Vars_Coord_Label,"Velocity [m/s]"); strcpy(Set_Vars_Coord_Var,"v"); lmin = 0; lmax = 1000000; }
        if (!strcmp(token, "t") || !strcmp(token, "time") || !strcmp(token, "tof"))
          { Set_Vars_Coord_Type = DEFS->COORD_T; strcpy(Set_Vars_Coord_Label,"TOF [s]"); strcpy(Set_Vars_Coord_Var,"t"); lmin = 0; lmax = .1; }
        if ((!strcmp(token, "p") || !strcmp(token, "i") || !strcmp(token, "intensity") || !strcmp(token, "flux")))
          { Set_Vars_Coord_Type = DEFS->COORD_P;
            strcpy(Set_Vars_Coord_Label,"Intensity");
            strncat(Set_Vars_Coord_Label, " [n/s", 30);
            if (Vars->Flag_per_cm2) strncat(Set_Vars_Coord_Label, "/cm2", 30);
            if (XY > 1 && Vars->Coord_Number)
              strncat(Set_Vars_Coord_Label, "/bin", 30);
            strncat(Set_Vars_Coord_Label, "]", 30);
            strcpy(Set_Vars_Coord_Var,"I");
            lmin = 0; lmax = FLT_MAX;
            if (Flag_auto>0) Flag_auto=0;
          }

        if (!strcmp(token, "vx"))
          { Set_Vars_Coord_Type = DEFS->COORD_VX; strcpy(Set_Vars_Coord_Label,"vx [m/s]"); strcpy(Set_Vars_Coord_Var,"vx"); lmin = -1000; lmax = 1000; }
        if (!strcmp(token, "vy"))
          { Set_Vars_Coord_Type = DEFS->COORD_VY; strcpy(Set_Vars_Coord_Label,"vy [m/s]"); strcpy(Set_Vars_Coord_Var,"vy"); lmin = -1000; lmax = 1000; }
        if (!strcmp(token, "vz"))
          { Set_Vars_Coord_Type = DEFS->COORD_VZ; strcpy(Set_Vars_Coord_Label,"vz [m/s]"); strcpy(Set_Vars_Coord_Var,"vz"); lmin = -10000; lmax = 10000; }
        if (!strcmp(token, "kx"))
          { Set_Vars_Coord_Type = DEFS->COORD_KX; strcpy(Set_Vars_Coord_Label,"kx [Angs-1]"); strcpy(Set_Vars_Coord_Var,"kx"); lmin = -1; lmax = 1; }
        if (!strcmp(token, "ky"))
          { Set_Vars_Coord_Type = DEFS->COORD_KY; strcpy(Set_Vars_Coord_Label,"ky [Angs-1]"); strcpy(Set_Vars_Coord_Var,"ky"); lmin = -1; lmax = 1; }
        if (!strcmp(token, "kz"))
          { Set_Vars_Coord_Type = DEFS->COORD_KZ; strcpy(Set_Vars_Coord_Label,"kz [Angs-1]"); strcpy(Set_Vars_Coord_Var,"kz"); lmin = -10; lmax = 10; }
        if (!strcmp(token, "sx"))
          { Set_Vars_Coord_Type = DEFS->COORD_SX; strcpy(Set_Vars_Coord_Label,"sx [1]"); strcpy(Set_Vars_Coord_Var,"sx"); lmin = -1; lmax = 1; }
        if (!strcmp(token, "sy"))
          { Set_Vars_Coord_Type = DEFS->COORD_SY; strcpy(Set_Vars_Coord_Label,"sy [1]"); strcpy(Set_Vars_Coord_Var,"sy"); lmin = -1; lmax = 1; }
        if (!strcmp(token, "sz"))
          { Set_Vars_Coord_Type = DEFS->COORD_SZ; strcpy(Set_Vars_Coord_Label,"sz [1]"); strcpy(Set_Vars_Coord_Var,"sz"); lmin = -1; lmax = 1; }

        if (!strcmp(token, "energy") || !strcmp(token, "omega") || !strcmp(token, "e"))
          { Set_Vars_Coord_Type = DEFS->COORD_ENERGY; strcpy(Set_Vars_Coord_Label,"Energy [meV]"); strcpy(Set_Vars_Coord_Var,"E"); lmin = 0; lmax = 100; }
        if (!strcmp(token, "lambda") || !strcmp(token, "wavelength") || !strcmp(token, "l"))
          { Set_Vars_Coord_Type = DEFS->COORD_LAMBDA; strcpy(Set_Vars_Coord_Label,"Wavelength [Angs]"); strcpy(Set_Vars_Coord_Var,"L"); lmin = 0; lmax = 100; }
        if (!strcmp(token, "radius") || !strcmp(token, "r"))
          { Set_Vars_Coord_Type = DEFS->COORD_RADIUS; strcpy(Set_Vars_Coord_Label,"Radius [m]"); strcpy(Set_Vars_Coord_Var,"xy"); lmin = 0; lmax = xmax; }
        if (!strcmp(token, "xy"))
          { Set_Vars_Coord_Type = DEFS->COORD_XY; strcpy(Set_Vars_Coord_Label,"Radius (xy) [m]"); strcpy(Set_Vars_Coord_Var,"xy"); lmin = 0; lmax = xmax; }
        if (!strcmp(token, "yz"))
          { Set_Vars_Coord_Type = DEFS->COORD_YZ; strcpy(Set_Vars_Coord_Label,"Radius (yz) [m]"); strcpy(Set_Vars_Coord_Var,"yz"); lmin = 0; lmax = xmax; }
        if (!strcmp(token, "xz"))
          { Set_Vars_Coord_Type = DEFS->COORD_XZ; strcpy(Set_Vars_Coord_Label,"Radius (xz) [m]"); strcpy(Set_Vars_Coord_Var,"xz"); lmin = 0; lmax = xmax; }
        if (!strcmp(token, "vxy"))
          { Set_Vars_Coord_Type = DEFS->COORD_VXY; strcpy(Set_Vars_Coord_Label,"Radial Velocity (xy) [m]"); strcpy(Set_Vars_Coord_Var,"Vxy"); lmin = 0; lmax = 2000; }
        if (!strcmp(token, "kxy"))
          { Set_Vars_Coord_Type = DEFS->COORD_KXY; strcpy(Set_Vars_Coord_Label,"Radial Wavevector (xy) [Angs-1]"); strcpy(Set_Vars_Coord_Var,"Kxy"); lmin = 0; lmax = 2; }
        if (!strcmp(token, "vyz"))
          { Set_Vars_Coord_Type = DEFS->COORD_VYZ; strcpy(Set_Vars_Coord_Label,"Radial Velocity (yz) [m]"); strcpy(Set_Vars_Coord_Var,"Vyz"); lmin = 0; lmax = 2000; }
        if (!strcmp(token, "kyz"))
          { Set_Vars_Coord_Type = DEFS->COORD_KYZ; strcpy(Set_Vars_Coord_Label,"Radial Wavevector (yz) [Angs-1]"); strcpy(Set_Vars_Coord_Var,"Kyz"); lmin = 0; lmax = 2; }
        if (!strcmp(token, "vxz"))
          { Set_Vars_Coord_Type = DEFS->COORD_VXZ; strcpy(Set_Vars_Coord_Label,"Radial Velocity (xz) [m]"); strcpy(Set_Vars_Coord_Var,"Vxz"); lmin = 0; lmax = 2000; }
        if (!strcmp(token, "kxz"))
          { Set_Vars_Coord_Type = DEFS->COORD_KXZ; strcpy(Set_Vars_Coord_Label,"Radial Wavevector (xz) [Angs-1]"); strcpy(Set_Vars_Coord_Var,"Kxz"); lmin = 0; lmax = 2; }
        if (!strcmp(token, "angle") || !strcmp(token, "a"))
          { Set_Vars_Coord_Type = DEFS->COORD_ANGLE; strcpy(Set_Vars_Coord_Label,"Angle [deg]"); strcpy(Set_Vars_Coord_Var,"A"); lmin = -50; lmax = 50; }
        if (!strcmp(token, "hdiv")|| !strcmp(token, "divergence") || !strcmp(token, "xdiv") || !strcmp(token, "hd") || !strcmp(token, "dx"))
          { Set_Vars_Coord_Type = DEFS->COORD_HDIV; strcpy(Set_Vars_Coord_Label,"Hor. Divergence [deg]"); strcpy(Set_Vars_Coord_Var,"hd"); lmin = -5; lmax = 5; }
        if (!strcmp(token, "vdiv") || !strcmp(token, "ydiv") || !strcmp(token, "vd") || !strcmp(token, "dy"))
          { Set_Vars_Coord_Type = DEFS->COORD_VDIV; strcpy(Set_Vars_Coord_Label,"Vert. Divergence [deg]"); strcpy(Set_Vars_Coord_Var,"vd"); lmin = -5; lmax = 5; }
        if (!strcmp(token, "theta") || !strcmp(token, "longitude") || !strcmp(token, "th"))
          { Set_Vars_Coord_Type = DEFS->COORD_THETA; strcpy(Set_Vars_Coord_Label,"Longitude [deg]"); strcpy(Set_Vars_Coord_Var,"th"); lmin = -180; lmax = 180; }
        if (!strcmp(token, "phi") || !strcmp(token, "lattitude") || !strcmp(token, "ph"))
          { Set_Vars_Coord_Type = DEFS->COORD_PHI; strcpy(Set_Vars_Coord_Label,"Lattitude [deg]"); strcpy(Set_Vars_Coord_Var,"ph"); lmin = -180; lmax = 180; }
        if (!strcmp(token, "ncounts") || !strcmp(token, "n") || !strcmp(token, "neutron"))
          { Set_Vars_Coord_Type = DEFS->COORD_NCOUNT; strcpy(Set_Vars_Coord_Label,"Neutron ID [1]"); strcpy(Set_Vars_Coord_Var,"n"); lmin = 0; lmax = mcget_ncount(); if (Flag_auto>0) Flag_auto=0; }
        if (!strcmp(token, "id") || !strcmp(token, "pixel"))
          { Set_Vars_Coord_Type = DEFS->COORD_PIXELID; 
            strcpy(Set_Vars_Coord_Label,"Pixel ID [1]"); 
            strcpy(Set_Vars_Coord_Var,"id"); lmin = 0; lmax = FLT_MAX; 
            if (Flag_auto>0) Flag_auto=0;
            Vars->Flag_List = 1; }
        if (!strcmp(token, "user") || !strcmp(token, "user1") || !strcmp(token, "u1"))
          { Set_Vars_Coord_Type = DEFS->COORD_USER1; strncpy(Set_Vars_Coord_Label,Vars->UserName1,30); strcpy(Set_Vars_Coord_Var,"U1"); lmin = -1e10; lmax = 1e10; }
        if (!strcmp(token, "user2") || !strcmp(token, "u2"))
          { Set_Vars_Coord_Type = DEFS->COORD_USER2; strncpy(Set_Vars_Coord_Label,Vars->UserName2,30); strcpy(Set_Vars_Coord_Var,"U2"); lmin = -1e10; lmax = 1e10; }
        if (!strcmp(token, "user3") || !strcmp(token, "u3"))
          { Set_Vars_Coord_Type = DEFS->COORD_USER3; strncpy(Set_Vars_Coord_Label,Vars->UserName3,30); strcpy(Set_Vars_Coord_Var,"U3"); lmin = -1e10; lmax = 1e10; }

        /* now stores variable keywords detected, if any */
        if (Set_Vars_Coord_Type != DEFS->COORD_NONE)
        {
          int Coord_Number = Vars->Coord_Number;
          if (Vars->Flag_log) { Set_Vars_Coord_Type |= DEFS->COORD_LOG; Vars->Flag_log = 0; }
          if (Flag_abs) { Set_Vars_Coord_Type |= DEFS->COORD_ABS; Flag_abs = 0; }
          if (Flag_auto != 0) { Set_Vars_Coord_Type |= DEFS->COORD_AUTO; 
            if (Flag_auto > 0) Flag_auto = 0; }
          if (Set_Coord_Mode == DEFS->COORD_SIGNAL)
          {
            Coord_Number = 0;
            Vars->Flag_signal = Set_Vars_Coord_Type;
          }
          else
          {
            if (Coord_Number < MONnD_COORD_NMAX)
            { Coord_Number++;
              Vars->Coord_Number = Coord_Number; 
              if (Set_Vars_Coord_Type != DEFS->COORD_PIXELID)
                Vars->Coord_NumberNoPixel++;
            }
            else if (Vars->Flag_Verbose) printf("Monitor_nD: %s reached max number of variables (%i).\n", Vars->compcurname, MONnD_COORD_NMAX);
          }
          Vars->Coord_Type[Coord_Number] = Set_Vars_Coord_Type;
          strncpy(Vars->Coord_Label[Coord_Number], Set_Vars_Coord_Label,30);
          strncpy(Vars->Coord_Var[Coord_Number], Set_Vars_Coord_Var,30);
          if (lmin > lmax) { XY = lmin; lmin=lmax; lmax = XY; }
          Vars->Coord_Min[Coord_Number] = lmin;
          Vars->Coord_Max[Coord_Number] = lmax;
          if (Set_Vars_Coord_Type == DEFS->COORD_NCOUNT || Set_Vars_Coord_Type == DEFS->COORD_PIXELID || Set_Vars_Coord_Type == DEFS->COORD_SIGNAL)
            Vars->Coord_Bin[Coord_Number] = 1;
          else
            Vars->Coord_Bin[Coord_Number] = 20;
          Set_Coord_Mode = DEFS->COORD_VAR;
          Flag_All = 0;
          Flag_No  = 0;
        } else {
          /* no variable name could be read from options */
          if (!iskeyword) {
            if (strcmp(token, "cm2") && strcmp(token, "incoming")
             && strcmp(token, "outgoing") && strcmp(token, "cm2")
             && strcmp(token, "cm^2") && strcmp(token, "float")
             && strcmp(token, "double") && strcmp(token, "binary")
             && strcmp(token, "steradian") && Vars->Flag_Verbose)
              printf("Monitor_nD: %s: unknown '%s' keyword in 'options'. Ignoring.\n", Vars->compcurname, token);
          }
        }
      carg++;
      } /* end if token */
    } /* end while carg */
    free(option_copy);
    if (carg == 128) printf("Monitor_nD: %s reached max number of tokens (%i). Skipping.\n", Vars->compcurname, 128);

    if ((Vars->Flag_Shape == DEFS->SHAPE_BOX) && (fabs(Vars->mzmax - Vars->mzmin) == 0)) Vars->Flag_Shape = DEFS->SHAPE_SQUARE;

    if (Vars->Flag_log == 1) Vars->Coord_Type[0] |= DEFS->COORD_LOG;
    if (Vars->Coord_Number == 0)
    { Vars->Flag_Auto_Limits=0; Vars->Flag_Multiple=0; Vars->Flag_List=0; }

    /* now setting Monitor Name from variable labels */
    strcpy(Vars->Monitor_Label,"");
    XY = 1; /* will contain total bin number */
    for (i = 0; i <= Vars->Coord_Number; i++)
    {
      if (Flag_auto != 0) Vars->Coord_Type[i] |= DEFS->COORD_AUTO;
      Set_Vars_Coord_Type = (Vars->Coord_Type[i] & (DEFS->COORD_LOG-1));
      if ((Set_Vars_Coord_Type == DEFS->COORD_X)
       || (Set_Vars_Coord_Type == DEFS->COORD_Y)
       || (Set_Vars_Coord_Type == DEFS->COORD_Z))
       strcpy(Short_Label[i],"Position");
      else
      if ((Set_Vars_Coord_Type == DEFS->COORD_THETA)
       || (Set_Vars_Coord_Type == DEFS->COORD_PHI)
       || (Set_Vars_Coord_Type == DEFS->COORD_ANGLE))
       strcpy(Short_Label[i],"Angle");
      else
      if ((Set_Vars_Coord_Type == DEFS->COORD_XY)
       || (Set_Vars_Coord_Type == DEFS->COORD_XZ)
       || (Set_Vars_Coord_Type == DEFS->COORD_YZ)
       || (Set_Vars_Coord_Type == DEFS->COORD_RADIUS))
       strcpy(Short_Label[i],"Radius");
      else
      if ((Set_Vars_Coord_Type == DEFS->COORD_VX)
       || (Set_Vars_Coord_Type == DEFS->COORD_VY)
       || (Set_Vars_Coord_Type == DEFS->COORD_VZ)
       || (Set_Vars_Coord_Type == DEFS->COORD_V)
       || (Set_Vars_Coord_Type == DEFS->COORD_VXY)
       || (Set_Vars_Coord_Type == DEFS->COORD_VYZ)
       || (Set_Vars_Coord_Type == DEFS->COORD_VXZ))
       strcpy(Short_Label[i],"Velocity");
      else
      if ((Set_Vars_Coord_Type == DEFS->COORD_KX)
       || (Set_Vars_Coord_Type == DEFS->COORD_KY)
       || (Set_Vars_Coord_Type == DEFS->COORD_KZ)
       || (Set_Vars_Coord_Type == DEFS->COORD_KXY)
       || (Set_Vars_Coord_Type == DEFS->COORD_KYZ)
       || (Set_Vars_Coord_Type == DEFS->COORD_KXZ)
       || (Set_Vars_Coord_Type == DEFS->COORD_K))
       strcpy(Short_Label[i],"Wavevector");
      else
      if ((Set_Vars_Coord_Type == DEFS->COORD_SX)
       || (Set_Vars_Coord_Type == DEFS->COORD_SY)
       || (Set_Vars_Coord_Type == DEFS->COORD_SZ))
       strcpy(Short_Label[i],"Spin");
      else
      if ((Set_Vars_Coord_Type == DEFS->COORD_HDIV)
       || (Set_Vars_Coord_Type == DEFS->COORD_VDIV))
       strcpy(Short_Label[i],"Divergence");
      else
      if (Set_Vars_Coord_Type == DEFS->COORD_ENERGY)
       strcpy(Short_Label[i],"Energy");
      else
      if (Set_Vars_Coord_Type == DEFS->COORD_LAMBDA)
       strcpy(Short_Label[i],"Wavelength");
      else
      if (Set_Vars_Coord_Type == DEFS->COORD_NCOUNT)
       strcpy(Short_Label[i],"Neutron_ID");
      else
      if (Set_Vars_Coord_Type == DEFS->COORD_PIXELID)
       strcpy(Short_Label[i],"Pixel_ID");
      else
      if (Set_Vars_Coord_Type == DEFS->COORD_T)
          strcpy(Short_Label[i],"Time_Of_Flight");
      else
      if (Set_Vars_Coord_Type == DEFS->COORD_P)
          strcpy(Short_Label[i],"Intensity");
      else
      if (Set_Vars_Coord_Type == DEFS->COORD_USER1)
          strncpy(Short_Label[i],Vars->UserName1,30);
      else
      if (Set_Vars_Coord_Type == DEFS->COORD_USER2)
          strncpy(Short_Label[i],Vars->UserName2,30);
      else
      if (Set_Vars_Coord_Type == DEFS->COORD_USER3)
          strncpy(Short_Label[i],Vars->UserName3,30);
      else
          strcpy(Short_Label[i],"Unknown");

      if (Vars->Coord_Type[i] & DEFS->COORD_ABS)
      { strcat(Vars->Coord_Label[i]," (abs)"); }

      if (Vars->Coord_Type[i] & DEFS->COORD_LOG)
      { strcat(Vars->Coord_Label[i]," (log)"); }

      strcat(Vars->Monitor_Label, " ");
      strcat(Vars->Monitor_Label, Short_Label[i]);
      XY *= Vars->Coord_Bin[i];

    } /* end for Short_Label */

    if ((Vars->Coord_Type[0] & (DEFS->COORD_LOG-1)) == DEFS->COORD_P) {
      strncat(Vars->Coord_Label[0], " [n/s", 30);
      if (Vars->Flag_per_cm2) strncat(Vars->Coord_Label[0], "/cm2", 30);

      if (XY > 1 && Vars->Coord_Number)
        strncat(Vars->Coord_Label[0], "/bin", 30);
      strncat(Vars->Coord_Label[0], "]", 30);
    }

    /* update label 'signal per bin' if more than 1 bin */
    if (XY > 1 && Vars->Coord_Number) {
      if (Vars->Flag_capture)
        printf("Monitor_nD: %s: Using capture flux weightening on %ld bins.\n"
               "WARNING     Use binned data with caution, and prefer monitor integral value (I,Ierr).\n", Vars->compcurname, (long)XY);
    }

    strcat(Vars->Monitor_Label, " Monitor");
    if (Vars->Flag_Shape == DEFS->SHAPE_SQUARE) strcat(Vars->Monitor_Label, " (Square)");
    if (Vars->Flag_Shape == DEFS->SHAPE_DISK)   strcat(Vars->Monitor_Label, " (Disk)");
    if (Vars->Flag_Shape == DEFS->SHAPE_SPHERE) strcat(Vars->Monitor_Label, " (Sphere)");
    if (Vars->Flag_Shape == DEFS->SHAPE_CYLIND) strcat(Vars->Monitor_Label, " (Cylinder)");
    if (Vars->Flag_Shape == DEFS->SHAPE_BANANA) strcat(Vars->Monitor_Label, " (Banana)");
    if (Vars->Flag_Shape == DEFS->SHAPE_BOX)    strcat(Vars->Monitor_Label, " (Box)");
    if (Vars->Flag_Shape == DEFS->SHAPE_PREVIOUS) strcat(Vars->Monitor_Label, " (on PREVIOUS)");
    if (Vars->Flag_Shape == DEFS->SHAPE_OFF) strcat(Vars->Monitor_Label, " (OFF geometry)");
    if ((Vars->Flag_Shape == DEFS->SHAPE_CYLIND) || (Vars->Flag_Shape == DEFS->SHAPE_BANANA) || (Vars->Flag_Shape == DEFS->SHAPE_SPHERE) || (Vars->Flag_Shape == DEFS->SHAPE_BOX))
    {
      if (strstr(Vars->option, "incoming"))
      {
        Vars->Flag_Shape = abs(Vars->Flag_Shape);
        strcat(Vars->Monitor_Label, " [in]");
      }
      else /* if strstr(Vars->option, "outgoing")) */
      {
        Vars->Flag_Shape = -abs(Vars->Flag_Shape);
        strcat(Vars->Monitor_Label, " [out]");
      }
    }
    if (Vars->Flag_UsePreMonitor == 1)
    {
        strcat(Vars->Monitor_Label, " at ");
        strncat(Vars->Monitor_Label, Vars->UserName1,30);
    }
    if (Vars->Flag_log == 1) strcat(Vars->Monitor_Label, " [log] ");

    /* now allocate memory to store variables in TRACE */

    /* Vars->Coord_Number  0   : intensity or signal
     * Vars->Coord_Number  1:n : detector variables */

    if ((Vars->Coord_NumberNoPixel != 2) && !Vars->Flag_Multiple && !Vars->Flag_List)
    { Vars->Flag_Multiple = 1; /* default is n1D */
      if (Vars->Coord_Number != Vars->Coord_NumberNoPixel) Vars->Flag_List = 1; }

    /* list and auto limits case : Vars->Flag_List or Vars->Flag_Auto_Limits
     * -> Buffer to flush and suppress after Vars->Flag_Auto_Limits
     */
    if ((Vars->Flag_Auto_Limits || Vars->Flag_List) && Vars->Coord_Number)
    { /* Dim : (Vars->Coord_Number+1)*Vars->Buffer_Block matrix (for p, dp) */
      Vars->Mon2D_Buffer = (double *)malloc((Vars->Coord_Number+1)*Vars->Buffer_Block*sizeof(double));
      if (Vars->Mon2D_Buffer == NULL)
      { printf("Monitor_nD: %s cannot allocate Vars->Mon2D_Buffer (%li). No list and auto limits.\n", Vars->compcurname, Vars->Buffer_Block*(Vars->Coord_Number+1)*sizeof(double)); Vars->Flag_List = 0; Vars->Flag_Auto_Limits = 0; }
      else
      {
        for (i=0; i < (Vars->Coord_Number+1)*Vars->Buffer_Block; Vars->Mon2D_Buffer[i++] = (double)0);
      }
      Vars->Buffer_Size = Vars->Buffer_Block;
    }

    /* 1D and n1D case : Vars->Flag_Multiple */
    if (Vars->Flag_Multiple && Vars->Coord_NumberNoPixel)
    { /* Dim : Vars->Coord_Number*Vars->Coord_Bin[i] vectors */
      Vars->Mon2D_N  = (double **)malloc((Vars->Coord_Number)*sizeof(double *));
      Vars->Mon2D_p  = (double **)malloc((Vars->Coord_Number)*sizeof(double *));
      Vars->Mon2D_p2 = (double **)malloc((Vars->Coord_Number)*sizeof(double *));
      if ((Vars->Mon2D_N == NULL) || (Vars->Mon2D_p == NULL) || (Vars->Mon2D_p2 == NULL))
      { fprintf(stderr,"Monitor_nD: %s n1D cannot allocate Vars->Mon2D_N/p/p2 (%li). Fatal.\n", Vars->compcurname, (Vars->Coord_Number)*sizeof(double *)); exit(-1); }
      for (i= 1; i <= Vars->Coord_Number; i++)
      {
        Vars->Mon2D_N[i-1]  = (double *)malloc(Vars->Coord_Bin[i]*sizeof(double));
        Vars->Mon2D_p[i-1]  = (double *)malloc(Vars->Coord_Bin[i]*sizeof(double));
        Vars->Mon2D_p2[i-1] = (double *)malloc(Vars->Coord_Bin[i]*sizeof(double));
        if ((Vars->Mon2D_N == NULL) || (Vars->Mon2D_p == NULL) || (Vars->Mon2D_p2 == NULL))
        { fprintf(stderr,"Monitor_nD: %s n1D cannot allocate %s Vars->Mon2D_N/p/p2[%li] (%li). Fatal.\n", Vars->compcurname, Vars->Coord_Var[i], i, (Vars->Coord_Bin[i])*sizeof(double *)); exit(-1); }
        else
        {
          for (j=0; j < Vars->Coord_Bin[i]; j++ )
          { Vars->Mon2D_N[i-1][j] = (double)0; Vars->Mon2D_p[i-1][j] = (double)0; Vars->Mon2D_p2[i-1][j] = (double)0; }
        }
      }
    }
    else /* 2D case : Vars->Coord_Number==2 and !Vars->Flag_Multiple and !Vars->Flag_List */
    if ((Vars->Coord_NumberNoPixel == 2) && !Vars->Flag_Multiple)
    { /* Dim : Vars->Coord_Bin[1]*Vars->Coord_Bin[2] matrix */
      Vars->Mon2D_N  = (double **)malloc((Vars->Coord_Bin[1])*sizeof(double *));
      Vars->Mon2D_p  = (double **)malloc((Vars->Coord_Bin[1])*sizeof(double *));
      Vars->Mon2D_p2 = (double **)malloc((Vars->Coord_Bin[1])*sizeof(double *));
      if ((Vars->Mon2D_N == NULL) || (Vars->Mon2D_p == NULL) || (Vars->Mon2D_p2 == NULL))
      { fprintf(stderr,"Monitor_nD: %s 2D cannot allocate %s Vars->Mon2D_N/p/p2 (%li). Fatal.\n", Vars->compcurname, Vars->Coord_Var[1], (Vars->Coord_Bin[1])*sizeof(double *)); exit(-1); }
      for (i= 0; i < Vars->Coord_Bin[1]; i++)
      {
        Vars->Mon2D_N[i]  = (double *)malloc(Vars->Coord_Bin[2]*sizeof(double));
        Vars->Mon2D_p[i]  = (double *)malloc(Vars->Coord_Bin[2]*sizeof(double));
        Vars->Mon2D_p2[i] = (double *)malloc(Vars->Coord_Bin[2]*sizeof(double));
        if ((Vars->Mon2D_N == NULL) || (Vars->Mon2D_p == NULL) || (Vars->Mon2D_p2 == NULL))
        { fprintf(stderr,"Monitor_nD: %s 2D cannot allocate %s Vars->Mon2D_N/p/p2[%li] (%li). Fatal.\n", Vars->compcurname, Vars->Coord_Var[1], i, (Vars->Coord_Bin[2])*sizeof(double *)); exit(-1); }
        else
        {
          for (j=0; j < Vars->Coord_Bin[2]; j++ )
          { Vars->Mon2D_N[i][j] = (double)0; Vars->Mon2D_p[i][j] = (double)0; Vars->Mon2D_p2[i][j] = (double)0; }
        }
      }
    }
    else {
      Vars->Mon2D_N = Vars->Mon2D_p = Vars->Mon2D_p2 = NULL;
    }
      /* no Mon2D allocated for
       * (Vars->Coord_Number != 2) && !Vars->Flag_Multiple && Vars->Flag_List */

    Vars->psum  = 0;
    Vars->p2sum = 0;
    Vars->Nsum  = 0;

    Vars->area  = fabs(Vars->mxmax - Vars->mxmin)*fabs(Vars->mymax - Vars->mymin)*1E4; /* in cm**2 for square and box shapes */
    Vars->Sphere_Radius = fabs(Vars->mxmax - Vars->mxmin)/2;
    if ((abs(Vars->Flag_Shape) == DEFS->SHAPE_DISK) || (abs(Vars->Flag_Shape) == DEFS->SHAPE_SPHERE))
    {
      Vars->area = PI*Vars->Sphere_Radius*Vars->Sphere_Radius*1E4; /* disk shapes */
    }


    if (Vars->area == 0 && abs(Vars->Flag_Shape) != DEFS->SHAPE_PREVIOUS ) {
      if (abs(Vars->Flag_Shape) != DEFS->SHAPE_OFF) {  
	Vars->Coord_Number = 0;
      }
    }
    if (Vars->Coord_Number == 0 && Vars->Flag_Verbose)
      printf("Monitor_nD: %s is unactivated (0D)\n", Vars->compcurname);
    Vars->Cylinder_Height = fabs(Vars->mymax - Vars->mymin);

    if (Vars->Flag_Verbose)
    {
      printf("Monitor_nD: %s is a %s.\n", Vars->compcurname, Vars->Monitor_Label);
      printf("Monitor_nD: version %s with options=%s\n", MONITOR_ND_LIB_H, Vars->option);
    }
    
    /* compute the product of bin dimensions for PixelID */
    Vars->Coord_BinProd[0]=1;
    for (i = 1; i <= Vars->Coord_Number; i++)
      Vars->Coord_BinProd[i]=Vars->Coord_Bin[i]*Vars->Coord_BinProd[i-1];
  } /* end Monitor_nD_Init */

/* ========================================================================= */
/* Monitor_nD_Trace: this routine is used to monitor one propagating neutron */
/* return values: 0=neutron was absorbed, -1=neutron was outside bounds, 1=neutron was measured*/
/* ========================================================================= */

int Monitor_nD_Trace(MonitornD_Defines_type *DEFS, MonitornD_Variables_type *Vars)
{

  double  XY=0, pp=0;
  int     retval;
  long    i =0, j =0;
  double  Coord[MONnD_COORD_NMAX];
  long    Coord_Index[MONnD_COORD_NMAX];
  char    While_End   =0;
  long    While_Buffer=0;
  char    Set_Vars_Coord_Type = DEFS->COORD_NONE;
  
  /* the logic below depends mainly on:
       Flag_List:        1=store 1 buffer, 2=list all, 3=re-use buffer 
       Flag_Auto_Limits: 0 (no auto limits/list), 1 (store events into Buffer), 2 (re-emit store events)
   */

  /* Vars->Flag_Auto_Limits=1: buffer full, we read the Buffer, and determine min and max bounds */
  if ((Vars->Buffer_Counter >= Vars->Buffer_Block) && (Vars->Flag_Auto_Limits == 1) && (Vars->Coord_Number > 0))
  {
    /* auto limits case : get limits in Buffer for each variable */
          /* Dim : (Vars->Coord_Number+1)*Vars->Buffer_Block matrix (for p, dp) */
    if (Vars->Flag_Verbose) printf("Monitor_nD: %s getting %li Auto Limits from List (%li events) in TRACE.\n", Vars->compcurname, Vars->Coord_Number, Vars->Buffer_Counter);
    for (i = 1; i <= Vars->Coord_Number; i++)
    {
      if (Vars->Coord_Type[i] & DEFS->COORD_AUTO)
      {
        Vars->Coord_Min[i] =  FLT_MAX;
        Vars->Coord_Max[i] = -FLT_MAX;
        for (j = 0; j < Vars->Buffer_Counter; j++)
        {
          XY = Vars->Mon2D_Buffer[i+j*(Vars->Coord_Number+1)];  /* scanning variables in Buffer */
          if (XY < Vars->Coord_Min[i]) Vars->Coord_Min[i] = XY;
          if (XY > Vars->Coord_Max[i]) Vars->Coord_Max[i] = XY;
        }
        if  (Vars->Flag_Verbose)  
          printf("  %s: min=%g max=%g\n", Vars->Coord_Var[i], Vars->Coord_Min[i], Vars->Coord_Max[i]);
      }
    }
    Vars->Flag_Auto_Limits = 2;  /* pass to 2nd auto limits step (read Buffer and generate new events to store in histograms) */
  } /* end if Flag_Auto_Limits == 1 */

  /* manage realloc for 'list all' if Buffer size exceeded: flush Buffer to file */
  if ((Vars->Buffer_Counter >= Vars->Buffer_Block) && (Vars->Flag_List >= 2))
  {
    if (Vars->Buffer_Size >= 1000000 || Vars->Flag_List == 3)
    { /* save current (possibly append) and re-use Buffer */
      Monitor_nD_Save(DEFS, Vars);
      Vars->Flag_List = 3;
      Vars->Buffer_Block = Vars->Buffer_Size;
      Vars->Buffer_Counter  = 0;
      Vars->Neutron_Counter = 0;
    }
    else
    {
      Vars->Mon2D_Buffer  = (double *)realloc(Vars->Mon2D_Buffer, (Vars->Coord_Number+1)*(Vars->Neutron_Counter+Vars->Buffer_Block)*sizeof(double));
      if (Vars->Mon2D_Buffer == NULL)
            { printf("Monitor_nD: %s cannot reallocate Vars->Mon2D_Buffer[%li] (%li). Skipping.\n", Vars->compcurname, i, (Vars->Neutron_Counter+Vars->Buffer_Block)*sizeof(double)); Vars->Flag_List = 1; }
      else { Vars->Buffer_Counter = 0; Vars->Buffer_Size = Vars->Neutron_Counter+Vars->Buffer_Block; }
    }
  } /* end if Buffer realloc */

  char    outsidebounds=0;
  while (!While_End)
  { /* we generate Coord[] and Coord_index[] from Buffer (auto limits) or passing neutron */
    if ((Vars->Flag_Auto_Limits == 2) && (Vars->Coord_Number > 0))
    { /* Vars->Flag_Auto_Limits == 2: read back from Buffer (Buffer is filled or auto limits have been computed) */
      if (While_Buffer < Vars->Buffer_Block)
      {
        /* first while loop (While_Buffer) */
        /* auto limits case : scan Buffer within limits and store in Mon2D */
        Coord[0] = pp = Vars->Mon2D_Buffer[While_Buffer*(Vars->Coord_Number+1)];

        for (i = 1; i <= Vars->Coord_Number; i++)
        {
          /* scanning variables in Buffer */
          if (Vars->Coord_Bin[i] <= 1) continue;
          XY = (Vars->Coord_Max[i]-Vars->Coord_Min[i]);

          Coord[i] = Vars->Mon2D_Buffer[i+While_Buffer*(Vars->Coord_Number+1)];
          if (XY > 0) Coord_Index[i] = floor((Coord[i]-Vars->Coord_Min[i])*Vars->Coord_Bin[i]/XY);
          else        Coord_Index[i] = 0;
          if (Vars->Flag_With_Borders)
          {
            if (Coord_Index[i] < 0)                   Coord_Index[i] = 0;
            if (Coord_Index[i] >= Vars->Coord_Bin[i]) Coord_Index[i] = Vars->Coord_Bin[i] - 1;
          }
        } /* end for */
        
        /* update the PixelID, we compute it from the previous variables index */
        if (Vars->Coord_NumberNoPixel < Vars->Coord_Number) /* there is a Pixel variable */
        for (i = 1; i <= Vars->Coord_Number; i++) {
          char Set_Vars_Coord_Type = (Vars->Coord_Type[i] & (DEFS->COORD_LOG-1));
          if (Set_Vars_Coord_Type == DEFS->COORD_PIXELID) {
            char flag_outside=0;
            Coord_Index[i] = Coord[i] = 0;
            for (j= 1; j < i; j++) {
              /* not for 1D variables with Bin=1 such as PixelID, NCOUNT, Intensity */
              if (Vars->Coord_Bin[j] == 1) continue; 
              if (0 > Coord_Index[j] || Coord_Index[j] >= Vars->Coord_Bin[j]) {
                flag_outside=1;
                Coord[i] = 0;
                break;
              }
              Coord[i] += Coord_Index[j]*Vars->Coord_BinProd[j-1];
            }
            if (!flag_outside) {
              Vars->Mon2D_Buffer[i+While_Buffer*(Vars->Coord_Number+1)] = Coord[i];
            }
          } /* end if PixelID */
        }
        While_Buffer++;
      } /* end if in Buffer */
      else /* (While_Buffer >= Vars->Buffer_Block) && (Vars->Flag_Auto_Limits == 2) */
      {
        Vars->Flag_Auto_Limits = 0;
        if (!Vars->Flag_List) /* free Buffer not needed anymore (no list to output) */
        { /* Dim : (Vars->Coord_Number+1)*Vars->Buffer_Block matrix (for p, p2) */
          free(Vars->Mon2D_Buffer); Vars->Mon2D_Buffer = NULL;
        }
        if (Vars->Flag_Verbose) printf("Monitor_nD: %s flushed %li Auto Limits from List (%li) in TRACE.\n", Vars->compcurname, Vars->Coord_Number, Vars->Buffer_Counter);
      }
    } /* if Vars->Flag_Auto_Limits == 2 */
    
    if (Vars->Flag_Auto_Limits != 2 || !Vars->Coord_Number) /* Vars->Flag_Auto_Limits == 0 (no auto limits/list) or 1 (store events into Buffer) */
    {
      /* automatically compute area and steradian solid angle when in AUTO mode */
      /* compute the steradian solid angle incoming on the monitor */
      double v;
      v=sqrt(Vars->cvx*Vars->cvx
            +Vars->cvy*Vars->cvy
            +Vars->cvz*Vars->cvz);
      if (Vars->min_x > Vars->cx) Vars->min_x = Vars->cx;
      if (Vars->max_x < Vars->cx) Vars->max_x = Vars->cx;
      if (Vars->min_y > Vars->cy) Vars->min_y = Vars->cy;
      if (Vars->max_y < Vars->cy) Vars->max_y = Vars->cy;
      Vars->mean_p  += Vars->cp;
      if (v) {
        Vars->mean_dx += Vars->cp*fabs(Vars->cvx/v);
        Vars->mean_dy += Vars->cp*fabs(Vars->cvy/v);
      }

      for (i = 0; i <= Vars->Coord_Number; i++)
      { /* handle current neutron : last while */
        XY = 0;
        Set_Vars_Coord_Type = (Vars->Coord_Type[i] & (DEFS->COORD_LOG-1));
        /* get values for variables to monitor */
        if (Set_Vars_Coord_Type == DEFS->COORD_X) XY = Vars->cx;
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_Y) XY = Vars->cy;
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_Z) XY = Vars->cz;
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_VX) XY = Vars->cvx;
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_VY) XY = Vars->cvy;
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_VZ) XY = Vars->cvz;
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_KX) XY = V2K*Vars->cvx;
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_KY) XY = V2K*Vars->cvy;
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_KZ) XY = V2K*Vars->cvz;
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_SX) XY = Vars->csx;
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_SY) XY = Vars->csy;
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_SZ) XY = Vars->csz;
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_T) XY = Vars->ct;
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_P) XY = Vars->cp;
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_HDIV) XY = RAD2DEG*atan2(Vars->cvx,Vars->cvz);
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_VDIV) XY = RAD2DEG*atan2(Vars->cvy,Vars->cvz);
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_V) XY = sqrt(Vars->cvx*Vars->cvx+Vars->cvy*Vars->cvy+Vars->cvz*Vars->cvz);
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_RADIUS)
          XY = sqrt(Vars->cx*Vars->cx+Vars->cy*Vars->cy+Vars->cz*Vars->cz);
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_XY)
          XY = sqrt(Vars->cx*Vars->cx+Vars->cy*Vars->cy)*(Vars->cx > 0 ? 1 : -1);
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_YZ) XY = sqrt(Vars->cy*Vars->cy+Vars->cz*Vars->cz);
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_XZ)
          XY = sqrt(Vars->cx*Vars->cx+Vars->cz*Vars->cz);
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_VXY) XY = sqrt(Vars->cvx*Vars->cvx+Vars->cvy*Vars->cvy);
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_VXZ) XY = sqrt(Vars->cvx*Vars->cvx+Vars->cvz*Vars->cvz);
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_VYZ) XY = sqrt(Vars->cvy*Vars->cvy+Vars->cvz*Vars->cvz);
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_K) { XY = sqrt(Vars->cvx*Vars->cvx+Vars->cvy*Vars->cvy+Vars->cvz*Vars->cvz);  XY *= V2K; }
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_KXY) { XY = sqrt(Vars->cvx*Vars->cvx+Vars->cvy*Vars->cvy);  XY *= V2K; }
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_KXZ) { XY = sqrt(Vars->cvx*Vars->cvx+Vars->cvz*Vars->cvz);  XY *= V2K; }
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_KYZ) { XY = sqrt(Vars->cvy*Vars->cvy+Vars->cvz*Vars->cvz);  XY *= V2K; }
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_ENERGY) { XY = Vars->cvx*Vars->cvx+Vars->cvy*Vars->cvy+Vars->cvz*Vars->cvz;  XY *= VS2E; }
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_LAMBDA) { XY = sqrt(Vars->cvx*Vars->cvx+Vars->cvy*Vars->cvy+Vars->cvz*Vars->cvz);  XY *= V2K; if (XY != 0) XY = 2*PI/XY; }
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_NCOUNT) XY = Vars->Neutron_Counter;
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_ANGLE)
        {  XY = sqrt(Vars->cvx*Vars->cvx+Vars->cvy*Vars->cvy);
           if (Vars->cvz != 0)
                XY = RAD2DEG*atan2(XY,Vars->cvz)*(Vars->cx > 0 ? 1 : -1);
           else XY = 0;
        }
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_THETA)  { if (Vars->cz != 0) XY = RAD2DEG*atan2(Vars->cx,Vars->cz); }
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_PHI) { if (Vars->cz != 0) XY = RAD2DEG*asin(Vars->cy/Vars->cz); }
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_USER1) XY = Vars->UserVariable1;
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_USER2) XY = Vars->UserVariable2;
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_USER3) XY = Vars->UserVariable3;
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_PIXELID && !Vars->Flag_Auto_Limits) {
          /* compute the PixelID from previous coordinates 
             the PixelID is the product of Coord_Index[i] in the detector geometry 
             pixelID = sum( Coord_Index[j]*prod(Vars->Coord_Bin[1:(j-1)]) )
             
             this does not apply when we store events in the buffer as Coord_Index
             is not set. Then the pixelID will be re-computed during SAVE.
          */
          char flag_outside=0;
          for (j= 1; j < i; j++) {
            /* not for 1D variables with Bin=1 such as PixelID, NCOUNT, Intensity */
            if (Vars->Coord_Bin[j] <= 1) continue; 
            if (0 > Coord_Index[j] || Coord_Index[j] >= Vars->Coord_Bin[j]) { 
              flag_outside=1; XY=0; break;
            }
            XY += Coord_Index[j]*Vars->Coord_BinProd[j-1];
          }
	  if (Vars->Flag_mantid && Vars->Flag_OFF && Vars->OFF_polyidx >=0) XY=Vars->OFF_polyidx;
          if (!flag_outside) XY += Vars->Coord_Min[i];
        }
        
        /* handle 'abs' and 'log' keywords */
        if (Vars->Coord_Type[i] & DEFS->COORD_ABS) XY=fabs(XY);

        if (Vars->Coord_Type[i] & DEFS->COORD_LOG) /* compute log of variable if requested */
        {  if (XY > 0) XY = log(XY)/log(10);
           else        XY = -100; }

        Coord[i] = XY; Coord_Index[i] = 0;
        if (i == 0) { pp = XY; Coord_Index[i] = 0; }
        else {
        /* check bounds for variables which have no automatic limits */
          if ((!Vars->Flag_Auto_Limits || !(Vars->Coord_Type[i] & DEFS->COORD_AUTO)) && Vars->Coord_Bin[i]>1)
          { /* compute index in histograms for each variable to monitor */
            XY = (Vars->Coord_Max[i]-Vars->Coord_Min[i]);
            if (XY > 0) Coord_Index[i] = floor((Coord[i]-Vars->Coord_Min[i])*Vars->Coord_Bin[i]/XY);
            if (Vars->Flag_With_Borders)
            {
              if (Coord_Index[i] >= Vars->Coord_Bin[i]) Coord_Index[i] = Vars->Coord_Bin[i] - 1;
              if (Coord_Index[i] < 0) Coord_Index[i] = 0;
            }
            //if (0 > Coord_Index[i] || Coord_Index[i] >= Vars->Coord_Bin[i])
            //  outsidebounds=1;
          } /* else will get Index later from Buffer when Flag_Auto_Limits == 2 */
        }
        
      } /* end for i */
      While_End = 1;
    }/* end else if Vars->Flag_Auto_Limits == 2 */
    
    /* ====================================================================== */
    /* store n1d/2d neutron from Buffer (Auto_Limits == 2) or current neutron in while */
    if (Vars->Flag_Auto_Limits != 1) /* not when storing auto limits Buffer */
    {
      /* apply per cm2 */
      if (Vars->Flag_per_cm2 && Vars->area != 0)
        pp /= Vars->area;

      /* 2D case : Vars->Coord_Number==2 and !Vars->Flag_Multiple and !Vars->Flag_List */
      if ( Vars->Coord_NumberNoPixel == 2 && !Vars->Flag_Multiple)
      { /* Dim : Vars->Coord_Bin[1]*Vars->Coord_Bin[2] matrix */
        
        i = Coord_Index[1];
        j = Coord_Index[2];
        if (i >= 0 && i < Vars->Coord_Bin[1] && j >= 0 && j < Vars->Coord_Bin[2])
        {
          if (Vars->Mon2D_N) { 
            Vars->Mon2D_N[i][j]++;
            Vars->Mon2D_p[i][j] += pp;
            Vars->Mon2D_p2[i][j] += pp*pp;
          }
        } else {
          outsidebounds=1; 
        }
      } else {
        /* 1D and n1D case : Vars->Flag_Multiple */
        /* Dim : Vars->Coord_Number*Vars->Coord_Bin[i] vectors (intensity is not included) */
          
        for (i= 1; i <= Vars->Coord_Number; i++) {
          j = Coord_Index[i];
          if (j >= 0 && j < Vars->Coord_Bin[i]) {
            if  (Vars->Flag_Multiple && Vars->Mon2D_N) {
              Vars->Mon2D_N[i-1][j]++;
              Vars->Mon2D_p[i-1][j]  += pp;
              Vars->Mon2D_p2[i-1][j] += pp*pp;
            }
          } else { 
            outsidebounds=1;
            break;
          }
        }
      }
    } /* end (Vars->Flag_Auto_Limits != 1) */
    
    if (Vars->Flag_Auto_Limits != 2 && !outsidebounds) /* not when reading auto limits Buffer */
    { /* now store Coord into Buffer (no index needed) if necessary (list or auto limits) */
      if ((Vars->Buffer_Counter < Vars->Buffer_Block) && ((Vars->Flag_List) || (Vars->Flag_Auto_Limits == 1)))
      {
          
        for (i = 0; i <= Vars->Coord_Number; i++)
        {
          Vars->Mon2D_Buffer[i + Vars->Neutron_Counter*(Vars->Coord_Number+1)] = Coord[i];
        }
        Vars->Buffer_Counter++;
        if (Vars->Flag_Verbose && (Vars->Buffer_Counter >= Vars->Buffer_Block) && (Vars->Flag_List == 1)) 
          printf("Monitor_nD: %s %li neutrons stored in List.\n", Vars->compcurname, Vars->Buffer_Counter);
      }
      Vars->Neutron_Counter++;
    } /* end (Vars->Flag_Auto_Limits != 2) */
    
  } /* end while */
  Vars->Nsum++;
  Vars->psum  += pp;
  Vars->p2sum += pp*pp;

  /*determine return value: 1:neutron was in bounds and measured, -1: outside bounds, 0: outside bounds, should be absorbed.*/
  if(outsidebounds){
      if(Vars->Flag_Absorb){
          return 0;
      }else{
          return -1;
      }
  }
  return 1;
} /* end Monitor_nD_Trace */

/* ========================================================================= */
/* Monitor_nD_Save: this routine is used to save data files                  */
/* ========================================================================= */

MCDETECTOR Monitor_nD_Save(MonitornD_Defines_type *DEFS, MonitornD_Variables_type *Vars)
  {
    char   *fname;
    long    i,j;
    double *p0m = NULL;
    double *p1m = NULL;
    double *p2m = NULL;
    char    Coord_X_Label[CHAR_BUF_LENGTH];
    double  min1d, max1d;
    double  min2d, max2d;
    long    bin1d, bin2d;
    char    While_End = 0;
    long    While_Buffer = 0;
    double  XY=0, pp=0;
    double  Coord[MONnD_COORD_NMAX];
    long    Coord_Index[MONnD_COORD_NMAX];
    char    label[CHAR_BUF_LENGTH];
    double  ratio;

    MCDETECTOR detector;

    ratio = 100.0*mcget_run_num()/mcget_ncount();
    if (Vars->Flag_Verbose && Vars->Flag_per_cm2) {
      printf("Monitor_nD: %s: active flat detector area is %g [cm^2], total area is %g [cm^2]\n",
        Vars->compcurname, (Vars->max_x-Vars->min_x)
                          *(Vars->max_y-Vars->min_y)*1E4, Vars->area);
      printf("Monitor_nD: %s: beam solid angle is %g [st] (%g x %g [deg^2])\n",
        Vars->compcurname,
        2*fabs(2*atan(Vars->mean_dx/Vars->mean_p)
         *sin(2*atan(Vars->mean_dy/Vars->mean_p)/2)),
        atan(Vars->mean_dx/Vars->mean_p)*RAD2DEG,
        atan(Vars->mean_dy/Vars->mean_p)*RAD2DEG);
    }

    /* check Buffer flush when end of simulation reached */
    if ((Vars->Buffer_Counter <= Vars->Buffer_Block) && Vars->Flag_Auto_Limits && Vars->Mon2D_Buffer && Vars->Buffer_Counter)
    {
      /* Get Auto Limits */
      if (Vars->Flag_Verbose) printf("Monitor_nD: %s getting %li Auto Limits from List (%li events).\n", Vars->compcurname, Vars->Coord_Number, Vars->Buffer_Counter);

      for (i = 1; i <= Vars->Coord_Number; i++)
      {
        if ((Vars->Coord_Type[i] & DEFS->COORD_AUTO) && Vars->Coord_Bin[i] > 1)
        {
          Vars->Coord_Min[i] = FLT_MAX;
          Vars->Coord_Max[i] = -FLT_MAX;
          for (j = 0; j < Vars->Buffer_Counter; j++)
          {
            XY = Vars->Mon2D_Buffer[i+j*(Vars->Coord_Number+1)];  /* scanning variables in Buffer */
            if (XY < Vars->Coord_Min[i]) Vars->Coord_Min[i] = XY;
            if (XY > Vars->Coord_Max[i]) Vars->Coord_Max[i] = XY;
          }
          if  (Vars->Flag_Verbose)  
            printf("  %s: min=%g max=%g in %li bins\n", Vars->Coord_Var[i], Vars->Coord_Min[i], Vars->Coord_Max[i], Vars->Coord_Bin[i]);
        }
      }
      Vars->Flag_Auto_Limits = 2;  /* pass to 2nd auto limits step */
      Vars->Buffer_Block = Vars->Buffer_Counter;

      while (!While_End)
      { /* we generate Coord[] and Coord_index[] from Buffer (auto limits) */
        /* simulation ended before Buffer was filled. Limits have to be computed, and stored events must be sent into histograms */
        
        if (While_Buffer < Vars->Buffer_Block)
        {
          /* first while loops (While_Buffer) */
          Coord[0] = Vars->Mon2D_Buffer[While_Buffer*(Vars->Coord_Number+1)];

          /* auto limits case : scan Buffer within limits and store in Mon2D */
          for (i = 1; i <= Vars->Coord_Number; i++)
          {
            /* scanning variables in Buffer */
            if (Vars->Coord_Bin[i] <= 1) Coord_Index[i] = 0;
            else {
              XY = (Vars->Coord_Max[i]-Vars->Coord_Min[i]);
              Coord[i] = Vars->Mon2D_Buffer[i+While_Buffer*(Vars->Coord_Number+1)];
              if (XY > 0) Coord_Index[i] = floor((Coord[i]-Vars->Coord_Min[i])*Vars->Coord_Bin[i]/XY);
              else Coord_Index[i] = 0;
              if (Vars->Flag_With_Borders)
              {
                if (Coord_Index[i] < 0) Coord_Index[i] = 0;
                if (Coord_Index[i] >= Vars->Coord_Bin[i]) Coord_Index[i] = Vars->Coord_Bin[i] - 1;
              }
            }
          } /* end for */

          /* update the PixelID, we compute it from the previous variables index */
          for (i = 1; i <= Vars->Coord_Number; i++) {
            char Set_Vars_Coord_Type = (Vars->Coord_Type[i] & (DEFS->COORD_LOG-1));
            if (Set_Vars_Coord_Type == DEFS->COORD_PIXELID) {
              char outsidebounds=0;
              Coord_Index[i] = Coord[i] = 0;
              for (j= 1; j < i; j++) {
                /* not for 1D variables with Bin=1 such as PixelID, NCOUNT, Intensity */
                if (Vars->Coord_Bin[j] == 1) continue; 
                if (0 > Coord_Index[j] || Coord_Index[j] >= Vars->Coord_Bin[j]) {
                  outsidebounds=1;
                  Coord[i] = 0;
                  break;
                }
                Coord[i] += Coord_Index[j]*Vars->Coord_BinProd[j-1];
              }
              if (!outsidebounds) {
                Vars->Mon2D_Buffer[i+While_Buffer*(Vars->Coord_Number+1)] = Coord[i];
              }
            } /* end if PixelID */
          }
          While_Buffer++;
        } /* end if in Buffer */
        else /* (While_Buffer >= Vars->Buffer_Block) && (Vars->Flag_Auto_Limits == 2) */
        {
          Vars->Flag_Auto_Limits = 0;
          While_End = 1;
          if (Vars->Flag_Verbose) printf("Monitor_nD: %s flushed %li Auto Limits from List (%li).\n", Vars->compcurname, Vars->Coord_Number, Vars->Buffer_Counter);
        }

        /* store n1d/2d section from Buffer */

        pp = Coord[0];
        /* apply per cm2 or per st */
        if (Vars->Flag_per_cm2 && Vars->area      != 0)
          pp /= Vars->area;
        
        /* 2D case : Vars->Coord_Number==2 and !Vars->Flag_Multiple and !Vars->Flag_List */
        if (!Vars->Flag_Multiple && Vars->Coord_NumberNoPixel == 2)
        { /* Dim : Vars->Coord_Bin[1]*Vars->Coord_Bin[2] matrix */
          i = Coord_Index[1];
          j = Coord_Index[2];
          if (i >= 0 && i < Vars->Coord_Bin[1] && j >= 0 && j < Vars->Coord_Bin[2])
          {
            if (Vars->Mon2D_N) {
              Vars->Mon2D_N[i][j]++;
              Vars->Mon2D_p[i][j] += pp;
              Vars->Mon2D_p2[i][j] += pp*pp;
            }
          } else if (Vars->Flag_Absorb) pp=0;
        }
        else
        /* 1D and n1D case : Vars->Flag_Multiple */
        { /* Dim : Vars->Coord_Number*Vars->Coord_Bin[i] vectors (intensity is not included) */
          for (i= 1; i <= Vars->Coord_Number; i++)
          {
            j = Coord_Index[i];
            if (j >= 0 && j < Vars->Coord_Bin[i])
            {
              if (Vars->Flag_Multiple && Vars->Mon2D_N) {
                Vars->Mon2D_N[i-1][j]++;
                Vars->Mon2D_p[i-1][j] += pp;
                Vars->Mon2D_p2[i-1][j] += pp*pp;
              }
            } else if (Vars->Flag_Absorb) {
              pp=0; break;
            }
          }
        } /* end store 2D/1D */
        
      } /* end while */
    } /* end Force Get Limits */

    /* write output files (sent to file as p[i*n + j] vectors) */
    if (Vars->Coord_Number == 0)
    {
      double Nsum;
      double psum, p2sum;
      Nsum = Vars->Nsum;
      psum = Vars->psum;
      p2sum= Vars->p2sum;
      if (Vars->Flag_signal != DEFS->COORD_P && Nsum > 0)
      { psum /=Nsum; p2sum /= Nsum*Nsum; }
      /* DETECTOR_OUT_0D(Vars->Monitor_Label, Vars->Nsum, Vars->psum, Vars->p2sum); */
      detector = mcdetector_out_0D(Vars->Monitor_Label, Nsum, psum, p2sum, Vars->compcurname, Vars->compcurpos);
    }
    else
    if (strlen(Vars->Mon_File) > 0)
    {
      fname = (char*)malloc(strlen(Vars->Mon_File)+10*Vars->Coord_Number);
      if (Vars->Flag_List && Vars->Mon2D_Buffer) /* List: DETECTOR_OUT_2D */
      {
       
        if (Vars->Flag_List >= 2) Vars->Buffer_Size = Vars->Neutron_Counter;
        if (Vars->Buffer_Size >= Vars->Neutron_Counter)
          Vars->Buffer_Size = Vars->Neutron_Counter;
        strcpy(fname,Vars->Mon_File);
        if (strchr(Vars->Mon_File,'.') == NULL) strcat(fname, "_list");

        strcpy(Coord_X_Label,"");
        for (i= 0; i <= Vars->Coord_Number; i++)
        {
          strcat(Coord_X_Label, Vars->Coord_Var[i]);
          strcat(Coord_X_Label, " ");
          if (strchr(Vars->Mon_File,'.') == NULL)
          { strcat(fname, "."); strcat(fname, Vars->Coord_Var[i]); }
        }
        if (Vars->Flag_Verbose) printf("Monitor_nD: %s write monitor file %s List (%lix%li).\n", Vars->compcurname, fname,bin2d,bin1d);

        /* handle the type of list output */
        strcpy(label, Vars->Monitor_Label);
        
        detector = mcdetector_out_list(
              label, "List of neutron events", Coord_X_Label,
              -Vars->Buffer_Size, Vars->Coord_Number+1,
              Vars->Mon2D_Buffer,
              fname, Vars->compcurname, Vars->compcurpos);
      }
      if (Vars->Flag_Multiple) /* n1D: DETECTOR_OUT_1D */
      {
        for (i= 0; i < Vars->Coord_Number; i++)
        {

          strcpy(fname,Vars->Mon_File);
          if (strchr(Vars->Mon_File,'.') == NULL)
          { strcat(fname, "."); strcat(fname, Vars->Coord_Var[i+1]); }
          sprintf(Coord_X_Label, "%s monitor", Vars->Coord_Label[i+1]);
          strcpy(label, Coord_X_Label);
          if (Vars->Coord_Bin[i+1] > 0) { /* 1D monitor */
            if (Vars->Flag_Verbose) printf("Monitor_nD: %s write monitor file %s 1D (%li).\n", Vars->compcurname, fname, Vars->Coord_Bin[i+1]);
            min1d = Vars->Coord_Min[i+1];
            max1d = Vars->Coord_Max[i+1];
            if (min1d == max1d) max1d = min1d+1e-6;
            p1m = (double *)malloc(Vars->Coord_Bin[i+1]*sizeof(double));
            p2m = (double *)malloc(Vars->Coord_Bin[i+1]*sizeof(double));
            if (p2m == NULL) /* use Raw Buffer line output */
            {
              if (Vars->Flag_Verbose) printf("Monitor_nD: %s cannot allocate memory for output. Using raw data.\n", Vars->compcurname);
              if (p1m != NULL) free(p1m);
              detector = mcdetector_out_1D(
              label,
              Vars->Coord_Label[i+1],
              Vars->Coord_Label[0],
              Vars->Coord_Var[i+1],
              min1d, max1d,
              Vars->Coord_Bin[i+1],
              Vars->Mon2D_N[i],Vars->Mon2D_p[i],Vars->Mon2D_p2[i],
              fname, Vars->compcurname, Vars->compcurpos);
            } /* if (p2m == NULL) */
            else
            {
              if (Vars->Flag_log != 0)
              {
                XY = FLT_MAX;
                for (j=0; j < Vars->Coord_Bin[i+1]; j++) /* search min of signal */
                  if ((XY > Vars->Mon2D_p[i][j]) && (Vars->Mon2D_p[i][j] > 0)) XY = Vars->Mon2D_p[i][j];
                if (XY <= 0) XY = -log(FLT_MAX)/log(10); else XY = log(XY)/log(10)-1;
              } /* if */

              for (j=0; j < Vars->Coord_Bin[i+1]; j++)
              {
                p1m[j] = Vars->Mon2D_p[i][j];
                p2m[j] = Vars->Mon2D_p2[i][j];
                if (Vars->Flag_signal != DEFS->COORD_P && Vars->Mon2D_N[i][j] > 0)
                { /* normalize mean signal to the number of events */
                  p1m[j] /= Vars->Mon2D_N[i][j];
                  p2m[j] /= Vars->Mon2D_N[i][j]*Vars->Mon2D_N[i][j];
                }
                if (Vars->Flag_log != 0)
                {
                  if ((p1m[j] > 0) && (p2m[j] > 0))
                  {
                    p2m[j] /= p1m[j]*p1m[j];
                    p1m[j] = log(p1m[j])/log(10);
                  }
                  else
                  {
                    p1m[j] = XY;
                    p2m[j] = 0;
                  }
                }
              } /* for */
              detector = mcdetector_out_1D(
                label,
                Vars->Coord_Label[i+1],
                Vars->Coord_Label[0],
                Vars->Coord_Var[i+1],
                min1d, max1d,
                Vars->Coord_Bin[i+1],
                Vars->Mon2D_N[i],p1m,p2m,
                fname, Vars->compcurname, Vars->compcurpos);

            } /* else */
            /* comment out 'free memory' lines to avoid loosing arrays if
               'detector' structure is used by other instrument parts
            if (p1m != NULL) free(p1m); p1m=NULL;
            if (p2m != NULL) free(p2m); p2m=NULL;
            */
          } else { /* 0d monitor */
            detector = mcdetector_out_0D(label, Vars->Mon2D_p[i][0], Vars->Mon2D_p2[i][0], Vars->Mon2D_N[i][0], Vars->compcurname, Vars->compcurpos);
          }


        } /* for */
      } /* if 1D */
      else
      if (Vars->Coord_NumberNoPixel == 2)  /* 2D: DETECTOR_OUT_2D */
      {
        strcpy(fname,Vars->Mon_File);

        p0m = (double *)malloc(Vars->Coord_Bin[1]*Vars->Coord_Bin[2]*sizeof(double));
        p1m = (double *)malloc(Vars->Coord_Bin[1]*Vars->Coord_Bin[2]*sizeof(double));
        p2m = (double *)malloc(Vars->Coord_Bin[1]*Vars->Coord_Bin[2]*sizeof(double));
        if (p2m == NULL)
        {
          if (Vars->Flag_Verbose) printf("Monitor_nD: %s cannot allocate memory for 2D array (%li). Skipping.\n", Vars->compcurname, 3*Vars->Coord_Bin[1]*Vars->Coord_Bin[2]*sizeof(double));
          /* comment out 'free memory' lines to avoid loosing arrays if
               'detector' structure is used by other instrument parts
          if (p0m != NULL) free(p0m);
          if (p1m != NULL) free(p1m);
          */
        }
        else
        {
          if (Vars->Flag_log != 0)
          {
            XY = FLT_MAX;
            for (i= 0; i < Vars->Coord_Bin[1]; i++)
              for (j= 0; j < Vars->Coord_Bin[2]; j++) /* search min of signal */
                if ((XY > Vars->Mon2D_p[i][j]) && (Vars->Mon2D_p[i][j]>0)) XY = Vars->Mon2D_p[i][j];
            if (XY <= 0) XY = -log(FLT_MAX)/log(10); else XY = log(XY)/log(10)-1;
          }
          for (i= 0; i < Vars->Coord_Bin[1]; i++)
          {
            for (j= 0; j < Vars->Coord_Bin[2]; j++)
            {
              long index;
              index = j + i*Vars->Coord_Bin[2];
              p0m[index] = Vars->Mon2D_N[i][j];
              p1m[index] = Vars->Mon2D_p[i][j];
              p2m[index] = Vars->Mon2D_p2[i][j];
              if (Vars->Flag_signal != DEFS->COORD_P && p0m[index] > 0)
              {
                  p1m[index] /= p0m[index];
                  p2m[index] /= p0m[index]*p0m[index];
              }

              if (Vars->Flag_log != 0)
              {
                if ((p1m[index] > 0) && (p2m[index] > 0))
                {
                  p2m[index] /= (p1m[index]*p1m[index]);
                  p1m[index] = log(p1m[index])/log(10);

                }
                else
                {
                  p1m[index] = XY;
                  p2m[index] = 0;
                }
              }
            }
          }
          if (strchr(Vars->Mon_File,'.') == NULL)
          { strcat(fname, "."); strcat(fname, Vars->Coord_Var[1]);
              strcat(fname, "_"); strcat(fname, Vars->Coord_Var[2]); }
          if (Vars->Flag_Verbose) printf("Monitor_nD: %s write monitor file %s 2D (%lix%li).\n", Vars->compcurname, fname, Vars->Coord_Bin[1], Vars->Coord_Bin[2]);

          min1d = Vars->Coord_Min[1];
          max1d = Vars->Coord_Max[1];
          if (min1d == max1d) max1d = min1d+1e-6;
          min2d = Vars->Coord_Min[2];
          max2d = Vars->Coord_Max[2];
          if (min2d == max2d) max2d = min2d+1e-6;
          strcpy(label, Vars->Monitor_Label);
          if (Vars->Coord_Bin[1]*Vars->Coord_Bin[2] > 1
           && Vars->Flag_signal == DEFS->COORD_P)
            strcat(label, " per bin");

          detector = mcdetector_out_2D(
            label,
            Vars->Coord_Label[1],
            Vars->Coord_Label[2],
            min1d, max1d,
            min2d, max2d,
            Vars->Coord_Bin[1],
            Vars->Coord_Bin[2],
            p0m,p1m,p2m,
            fname, Vars->compcurname, Vars->compcurpos);

          /* comment out 'free memory' lines to avoid loosing arrays if
               'detector' structure is used by other instrument parts
          if (p0m != NULL) free(p0m);
          if (p1m != NULL) free(p1m);
          if (p2m != NULL) free(p2m);
          */
        }
      }
      free(fname);
    }
    return(detector);
  } /* end Monitor_nD_Save */

/* ========================================================================= */
/* Monitor_nD_Finally: this routine is used to free memory                   */
/* ========================================================================= */

void Monitor_nD_Finally(MonitornD_Defines_type *DEFS,
  MonitornD_Variables_type *Vars)
  {
    int i;

    /* Now Free memory Mon2D.. */
    if ((Vars->Flag_Auto_Limits || Vars->Flag_List) && Vars->Coord_Number)
    { /* Dim : (Vars->Coord_Number+1)*Vars->Buffer_Block matrix (for p, dp) */
      if (Vars->Mon2D_Buffer != NULL) free(Vars->Mon2D_Buffer);
    }

    /* 1D and n1D case : Vars->Flag_Multiple */
    if (Vars->Flag_Multiple && Vars->Coord_Number)
    { /* Dim : Vars->Coord_Number*Vars->Coord_Bin[i] vectors */
      for (i= 0; i < Vars->Coord_Number; i++)
      {
        free(Vars->Mon2D_N[i]);
        free(Vars->Mon2D_p[i]);
        free(Vars->Mon2D_p2[i]);
      }
      free(Vars->Mon2D_N);
      free(Vars->Mon2D_p);
      free(Vars->Mon2D_p2);
    }


    /* 2D case : Vars->Coord_Number==2 and !Vars->Flag_Multiple and !Vars->Flag_List */
    if ((Vars->Coord_NumberNoPixel == 2) && !Vars->Flag_Multiple)
    { /* Dim : Vars->Coord_Bin[1]*Vars->Coord_Bin[2] matrix */
      for (i= 0; i < Vars->Coord_Bin[1]; i++)
      {
        free(Vars->Mon2D_N[i]);
        free(Vars->Mon2D_p[i]);
        free(Vars->Mon2D_p2[i]);
      }
      free(Vars->Mon2D_N);
      free(Vars->Mon2D_p);
      free(Vars->Mon2D_p2);
    }
  } /* end Monitor_nD_Finally */

/* ========================================================================= */
/* Monitor_nD_McDisplay: this routine is used to display component           */
/* ========================================================================= */

void Monitor_nD_McDisplay(MonitornD_Defines_type *DEFS,
  MonitornD_Variables_type *Vars)
  {
    double radius, h;
    double xmin;
    double xmax;
    double ymin;
    double ymax;
    double zmin;
    double zmax;
    int    i;
    double hdiv_min=-180, hdiv_max=180, vdiv_min=-90, vdiv_max=90;
    char   restricted = 0;

    radius = Vars->Sphere_Radius;
    h = Vars->Cylinder_Height;
    xmin = Vars->mxmin;
    xmax = Vars->mxmax;
    ymin = Vars->mymin;
    ymax = Vars->mymax;
    zmin = Vars->mzmin;
    zmax = Vars->mzmax;

    /* determine if there are angular limits set at start (no auto) in coord_types
     * cylinder/banana: look for hdiv
     * sphere: look for angle, radius (->atan2(val,radius)), hdiv, vdiv
     * this activates a 'restricted' flag, to draw a region as blades on cylinder/sphere
     */
    for (i= 0; i <= Vars->Coord_Number; i++)
    {
      int Set_Vars_Coord_Type;
      Set_Vars_Coord_Type = (Vars->Coord_Type[i] & (DEFS->COORD_LOG-1));
      if (Set_Vars_Coord_Type == DEFS->COORD_HDIV || Set_Vars_Coord_Type == DEFS->COORD_THETA)
      { hdiv_min = Vars->Coord_Min[i]; hdiv_max = Vars->Coord_Max[i]; restricted = 1; }
      else if (Set_Vars_Coord_Type == DEFS->COORD_VDIV || Set_Vars_Coord_Type == DEFS->COORD_PHI)
      { vdiv_min = Vars->Coord_Min[i]; vdiv_max = Vars->Coord_Max[i];restricted = 1;  }
      else if (Set_Vars_Coord_Type == DEFS->COORD_ANGLE)
      { hdiv_min = vdiv_min = Vars->Coord_Min[i];
        hdiv_max = vdiv_max = Vars->Coord_Max[i];
        restricted = 1; }
      else if (Set_Vars_Coord_Type == DEFS->COORD_RADIUS)
      { double angle;
        angle = RAD2DEG*atan2(Vars->Coord_Max[i], radius);
        hdiv_min = vdiv_min = angle;
        hdiv_max = vdiv_max = angle;
        restricted = 1; }
      else if (Set_Vars_Coord_Type == DEFS->COORD_Y && abs(Vars->Flag_Shape) == DEFS->SHAPE_SPHERE)
      {
        vdiv_min = atan2(ymin,radius)*RAD2DEG;
        vdiv_max = atan2(ymax,radius)*RAD2DEG;
        restricted = 1;
      }
    }
    /* full sphere */
    if ((!restricted && (abs(Vars->Flag_Shape) == DEFS->SHAPE_SPHERE))
    || abs(Vars->Flag_Shape) == DEFS->SHAPE_PREVIOUS)
    {
      mcdis_magnify("");
      mcdis_circle("xy",0,0,0,radius);
      mcdis_circle("xz",0,0,0,radius);
      mcdis_circle("yz",0,0,0,radius);
    }
    /* banana/cylinder/sphere portion */
    else
    if (restricted && ((abs(Vars->Flag_Shape) == DEFS->SHAPE_CYLIND)
                    || (abs(Vars->Flag_Shape) == DEFS->SHAPE_BANANA)
                    || (abs(Vars->Flag_Shape) == DEFS->SHAPE_SPHERE)))
    {
      int NH=24, NV=24;
      int ih, iv;
      double width, height;
      int issphere;
      issphere = (abs(Vars->Flag_Shape) == DEFS->SHAPE_SPHERE);
      width = (hdiv_max-hdiv_min)/NH;
      if (!issphere) NV=1; /* cylinder has vertical axis */
      else height= (vdiv_max-vdiv_min)/NV;
      
      /* check width and height of elements (sphere) to make sure the nb
         of plates remains limited */
      if (width < 10  && NH > 1) { width = 10;  NH=(hdiv_max-hdiv_min)/width; width=(hdiv_max-hdiv_min)/NH; }
      if (height < 10 && NV > 1) { height = 10; NV=(vdiv_max-vdiv_min)/height; height= (vdiv_max-vdiv_min)/NV; }
      
      mcdis_magnify("xyz");
      for(ih = 0; ih < NH; ih++)
        for(iv = 0; iv < NV; iv++)
        {
          double theta0, phi0, theta1, phi1;          /* angles in spherical coordinates */
          double x0,y0,z0,x1,y1,z1,x2,y2,z2,x3,y3,z3; /* vertices at plate edges */
          phi0 = (hdiv_min+ width*ih-90)*DEG2RAD;        /* in xz plane */
          phi1 = (hdiv_min+ width*(ih+1)-90)*DEG2RAD;
          if (issphere)
          {
            theta0= (vdiv_min+height* iv + 90)   *DEG2RAD; /* in vertical plane */
            theta1= (vdiv_min+height*(iv+1) + 90)*DEG2RAD;
            if (y0 < ymin) y0=ymin; 
            if (y0 > ymax) y0=ymax;
            if (y1 < ymin) y1=ymin; 
            if (y1 > ymax) y1=ymax;
            
            y0 = -radius*cos(theta0);            /* z with Z vertical */
            y1 = -radius*cos(theta1);
            if (y0 < ymin) y0=ymin;
            if (y0 > ymax) y0=ymax;
            if (y1 < ymin) y1=ymin;
            if (y1 > ymax) y1=ymax;
          } else {
            y0 = ymin;
            y1 = ymax;
            theta0=theta1=90*DEG2RAD;
          }

          x0 = radius*sin(theta0)*cos(phi0); /* x with Z vertical */
          z0 =-radius*sin(theta0)*sin(phi0); /* y with Z vertical */
          x1 = radius*sin(theta1)*cos(phi0); 
          z1 =-radius*sin(theta1)*sin(phi0);
          x2 = radius*sin(theta1)*cos(phi1); 
          z2 =-radius*sin(theta1)*sin(phi1);
          x3 = radius*sin(theta0)*cos(phi1); 
          z3 =-radius*sin(theta0)*sin(phi1);
          y2 = y1; y3 = y0;

          mcdis_multiline(5,
            x0,y0,z0,
            x1,y1,z1,
            x2,y2,z2,
            x3,y3,z3,
            x0,y0,z0);
        }
      if (Vars->Flag_mantid) {
	/* First define the base pixel type */
	double dt, dy;
	dt = (Vars->Coord_Max[1]-Vars->Coord_Min[1])/Vars->Coord_Bin[1];
	dy = (Vars->Coord_Max[2]-Vars->Coord_Min[2])/Vars->Coord_Bin[2];
	printf("MANTID_BANANA_DET:  %g, %g, %g, %g, %g, %li, %li, %g\n", radius, 
	       Vars->Coord_Min[1],Vars->Coord_Max[1], Vars->Coord_Min[2],Vars->Coord_Max[2], Vars->Coord_Bin[1], Vars->Coord_Bin[2], Vars->Coord_Min[4]); 
      }
    }
    /* disk (circle) */
    else
    if (abs(Vars->Flag_Shape) == DEFS->SHAPE_DISK)
    {
      mcdis_magnify("");
      mcdis_circle("xy",0,0,0,radius);
    }
    /* rectangle (square) */
    else
    if (abs(Vars->Flag_Shape) == DEFS->SHAPE_SQUARE)
    {
      mcdis_magnify("xy");
      mcdis_multiline(5, (double)xmin, (double)ymin, 0.0,
             (double)xmax, (double)ymin, 0.0,
             (double)xmax, (double)ymax, 0.0,
             (double)xmin, (double)ymax, 0.0,
             (double)xmin, (double)ymin, 0.0);
      
      if (Vars->Flag_mantid) {
	/* First define the base pixel type */
	double dx, dy;
	dx = (Vars->Coord_Max[1]-Vars->Coord_Min[1])/Vars->Coord_Bin[1];
	dy = (Vars->Coord_Max[2]-Vars->Coord_Min[2])/Vars->Coord_Bin[2];
	printf("MANTID_RECTANGULAR_DET:  %g, %g, %g, %g, %li, %li, %g\n", 
	       Vars->Coord_Min[1],Vars->Coord_Max[1], Vars->Coord_Min[2],Vars->Coord_Max[2], Vars->Coord_Bin[1], Vars->Coord_Bin[2], Vars->Coord_Min[4]);
      }
    }
    /* full cylinder/banana */
    else
    if (!restricted && ((abs(Vars->Flag_Shape) == DEFS->SHAPE_CYLIND) || (abs(Vars->Flag_Shape) == DEFS->SHAPE_BANANA)))
    {
      mcdis_magnify("xyz");
      mcdis_circle("xz", 0,  h/2.0, 0, radius);
      mcdis_circle("xz", 0, -h/2.0, 0, radius);
      mcdis_line(-radius, -h/2.0, 0, -radius, +h/2.0, 0);
      mcdis_line(+radius, -h/2.0, 0, +radius, +h/2.0, 0);
      mcdis_line(0, -h/2.0, -radius, 0, +h/2.0, -radius);
      mcdis_line(0, -h/2.0, +radius, 0, +h/2.0, +radius);
    }
    else
    /* box */
    if (abs(Vars->Flag_Shape) == DEFS->SHAPE_BOX)
    {
      mcdis_magnify("xyz");
      mcdis_multiline(5, xmin, ymin, zmin,
                   xmax, ymin, zmin,
                   xmax, ymax, zmin,
                   xmin, ymax, zmin,
                   xmin, ymin, zmin);
      mcdis_multiline(5, xmin, ymin, zmax,
                   xmax, ymin, zmax,
                   xmax, ymax, zmax,
                   xmin, ymax, zmax,
                   xmin, ymin, zmax);
      mcdis_line(xmin, ymin, zmin, xmin, ymin, zmax);
      mcdis_line(xmax, ymin, zmin, xmax, ymin, zmax);
      mcdis_line(xmin, ymax, zmin, xmin, ymax, zmax);
      mcdis_line(xmax, ymax, zmin, xmax, ymax, zmax);
    }
  } /* end Monitor_nD_McDisplay */

/* end of monitor_nd-lib.c */



#line 10345 "./PSI_DMC.c"

/* Instrument parameters. */
MCNUM mciplambda;
MCNUM mcipR;
MCNUM mcipR_curve;
char* mcipfilename;
MCNUM mcipD_PHI;
MCNUM mcipSHIFT;
MCNUM mcipPACK;
MCNUM mcipDw;
MCNUM mcipBARNS;

#define mcNUMIPAR 9
int mcnumipar = 9;
struct mcinputtable_struct mcinputtable[mcNUMIPAR+1] = {
  "lambda", &mciplambda, instr_type_double, "2.5666", 
  "R", &mcipR, instr_type_double, "0.87", 
  "R_curve", &mcipR_curve, instr_type_double, "0.87", 
  "filename", &mcipfilename, instr_type_string, "Na2Ca3Al2F14.laz", 
  "D_PHI", &mcipD_PHI, instr_type_double, "6", 
  "SHIFT", &mcipSHIFT, instr_type_double, "0", 
  "PACK", &mcipPACK, instr_type_double, "0.7", 
  "Dw", &mcipDw, instr_type_double, "0.8", 
  "BARNS", &mcipBARNS, instr_type_double, "1", 
  NULL, NULL, instr_type_double, ""
};

/* User declarations from instrument definition. */
#define mccompcurname  PSI_DMC
#define mccompcurtype  INSTRUMENT
#define mccompcurindex 0
#define mcposaPSI_DMC coords_set(0,0,0)
#define lambda mciplambda
#define R mcipR
#define R_curve mcipR_curve
#define filename mcipfilename
#define D_PHI mcipD_PHI
#define SHIFT mcipSHIFT
#define PACK mcipPACK
#define Dw mcipDw
#define BARNS mcipBARNS
#line 46 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  double mono_q = 1.8734;
  double OMA;  
  double RV;
  double y_mono = 0.025;
  double NV = 5;
  double d_phi_0;
  double TTM;
  double sample_radius = 0.008/2;
  double sample_height = 0.03;
  double can_radius = 0.0083/2;
  double can_height = 0.0303;
  double can_thick = 0.00015;
  
  /******Mirrorvalues*****/
  
  double alpha;
  double Qc=0.0217;
  double R0=0.995;
  double Mvalue=1.9;
  double W=1.0/250.0;
  
  double alpha_curve;
  double Qc_curve=0.0217;
  double R0_curve= 0.995;
  double Mvalue_curve=2.1;
  double W_curve=1.0/250.0;
  
  double ldiff=0.05;
  /* Curved guide element angle*/
  double angleGuideCurved;

#line 10419 "./PSI_DMC.c"
#undef BARNS
#undef Dw
#undef PACK
#undef SHIFT
#undef D_PHI
#undef filename
#undef R_curve
#undef R
#undef lambda
#undef mcposaPSI_DMC
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

/* Setting parameters for component 'source_arm' [1]. */
char mccsource_arm_profile[16384];
MCNUM mccsource_arm_percent;
MCNUM mccsource_arm_flag_save;
MCNUM mccsource_arm_minutes;

/* Setting parameters for component 'source' [2]. */
MCNUM mccsource_size;
MCNUM mccsource_yheight;
MCNUM mccsource_xwidth;
MCNUM mccsource_Lmin;
MCNUM mccsource_Lmax;
MCNUM mccsource_dist;
MCNUM mccsource_focus_xw;
MCNUM mccsource_focus_yh;
MCNUM mccsource_T1;
MCNUM mccsource_T2;
MCNUM mccsource_T3;
MCNUM mccsource_I1;
MCNUM mccsource_I2;
MCNUM mccsource_I3;
int mccsource_target_index;
MCNUM mccsource_lambda0;
MCNUM mccsource_dlambda;

/* Setting parameters for component 'PSDbefore_guides' [3]. */
int mccPSDbefore_guides_nx;
int mccPSDbefore_guides_ny;
char mccPSDbefore_guides_filename[16384];
MCNUM mccPSDbefore_guides_xmin;
MCNUM mccPSDbefore_guides_xmax;
MCNUM mccPSDbefore_guides_ymin;
MCNUM mccPSDbefore_guides_ymax;
MCNUM mccPSDbefore_guides_xwidth;
MCNUM mccPSDbefore_guides_yheight;
MCNUM mccPSDbefore_guides_restore_neutron;

/* Definition parameters for component 'l_mon_source' [4]. */
#define mccl_mon_source_nL 101
/* Setting parameters for component 'l_mon_source' [4]. */
char mccl_mon_source_filename[16384];
MCNUM mccl_mon_source_xmin;
MCNUM mccl_mon_source_xmax;
MCNUM mccl_mon_source_ymin;
MCNUM mccl_mon_source_ymax;
MCNUM mccl_mon_source_xwidth;
MCNUM mccl_mon_source_yheight;
MCNUM mccl_mon_source_Lmin;
MCNUM mccl_mon_source_Lmax;
MCNUM mccl_mon_source_restore_neutron;
int mccl_mon_source_nowritefile;

/* Setting parameters for component 'guide1' [5]. */
char mccguide1_reflect[16384];
MCNUM mccguide1_w1;
MCNUM mccguide1_h1;
MCNUM mccguide1_w2;
MCNUM mccguide1_h2;
MCNUM mccguide1_l;
MCNUM mccguide1_R0;
MCNUM mccguide1_Qc;
MCNUM mccguide1_alpha;
MCNUM mccguide1_m;
MCNUM mccguide1_W;

/* Setting parameters for component 'PSDbefore_curve' [6]. */
int mccPSDbefore_curve_nx;
int mccPSDbefore_curve_ny;
char mccPSDbefore_curve_filename[16384];
MCNUM mccPSDbefore_curve_xmin;
MCNUM mccPSDbefore_curve_xmax;
MCNUM mccPSDbefore_curve_ymin;
MCNUM mccPSDbefore_curve_ymax;
MCNUM mccPSDbefore_curve_xwidth;
MCNUM mccPSDbefore_curve_yheight;
MCNUM mccPSDbefore_curve_restore_neutron;

/* Setting parameters for component 'guide2' [7]. */
MCNUM mccguide2_w;
MCNUM mccguide2_h;
MCNUM mccguide2_r;
MCNUM mccguide2_Win;
MCNUM mccguide2_k;
MCNUM mccguide2_d;
MCNUM mccguide2_l;
MCNUM mccguide2_R0a;
MCNUM mccguide2_Qca;
MCNUM mccguide2_alphaa;
MCNUM mccguide2_ma;
MCNUM mccguide2_Wa;
MCNUM mccguide2_R0i;
MCNUM mccguide2_Qci;
MCNUM mccguide2_alphai;
MCNUM mccguide2_mi;
MCNUM mccguide2_Wi;
MCNUM mccguide2_R0s;
MCNUM mccguide2_Qcs;
MCNUM mccguide2_alphas;
MCNUM mccguide2_ms;
MCNUM mccguide2_Ws;

/* Setting parameters for component 'PSDafter_curve' [8]. */
int mccPSDafter_curve_nx;
int mccPSDafter_curve_ny;
char mccPSDafter_curve_filename[16384];
MCNUM mccPSDafter_curve_xmin;
MCNUM mccPSDafter_curve_xmax;
MCNUM mccPSDafter_curve_ymin;
MCNUM mccPSDafter_curve_ymax;
MCNUM mccPSDafter_curve_xwidth;
MCNUM mccPSDafter_curve_yheight;
MCNUM mccPSDafter_curve_restore_neutron;

/* Setting parameters for component 'bunker' [9]. */
char mccbunker_reflect[16384];
MCNUM mccbunker_w1;
MCNUM mccbunker_h1;
MCNUM mccbunker_w2;
MCNUM mccbunker_h2;
MCNUM mccbunker_l;
MCNUM mccbunker_R0;
MCNUM mccbunker_Qc;
MCNUM mccbunker_alpha;
MCNUM mccbunker_m;
MCNUM mccbunker_W;

/* Setting parameters for component 'guide3' [10]. */
char mccguide3_reflect[16384];
MCNUM mccguide3_w1;
MCNUM mccguide3_h1;
MCNUM mccguide3_w2;
MCNUM mccguide3_h2;
MCNUM mccguide3_l;
MCNUM mccguide3_R0;
MCNUM mccguide3_Qc;
MCNUM mccguide3_alpha;
MCNUM mccguide3_m;
MCNUM mccguide3_W;

/* Setting parameters for component 'guide4' [11]. */
char mccguide4_reflect[16384];
MCNUM mccguide4_w1;
MCNUM mccguide4_h1;
MCNUM mccguide4_w2;
MCNUM mccguide4_h2;
MCNUM mccguide4_l;
MCNUM mccguide4_R0;
MCNUM mccguide4_Qc;
MCNUM mccguide4_alpha;
MCNUM mccguide4_m;
MCNUM mccguide4_W;

/* Setting parameters for component 'window1' [12]. */
MCNUM mccwindow1_thickness;

/* Definition parameters for component 'ydist_fluxpos' [13]. */
#define mccydist_fluxpos_nx 11
/* Setting parameters for component 'ydist_fluxpos' [13]. */
char mccydist_fluxpos_filename[16384];
MCNUM mccydist_fluxpos_xmin;
MCNUM mccydist_fluxpos_xmax;
MCNUM mccydist_fluxpos_ymin;
MCNUM mccydist_fluxpos_ymax;
MCNUM mccydist_fluxpos_xwidth;
MCNUM mccydist_fluxpos_yheight;
MCNUM mccydist_fluxpos_restore_neutron;
int mccydist_fluxpos_nowritefile;

/* Setting parameters for component 'PSD_fluxpos' [14]. */
int mccPSD_fluxpos_nx;
int mccPSD_fluxpos_ny;
char mccPSD_fluxpos_filename[16384];
MCNUM mccPSD_fluxpos_xmin;
MCNUM mccPSD_fluxpos_xmax;
MCNUM mccPSD_fluxpos_ymin;
MCNUM mccPSD_fluxpos_ymax;
MCNUM mccPSD_fluxpos_xwidth;
MCNUM mccPSD_fluxpos_yheight;
MCNUM mccPSD_fluxpos_restore_neutron;

/* Definition parameters for component 'xdist_flux_pos' [15]. */
#define mccxdist_flux_pos_nx 11
/* Setting parameters for component 'xdist_flux_pos' [15]. */
char mccxdist_flux_pos_filename[16384];
MCNUM mccxdist_flux_pos_xmin;
MCNUM mccxdist_flux_pos_xmax;
MCNUM mccxdist_flux_pos_ymin;
MCNUM mccxdist_flux_pos_ymax;
MCNUM mccxdist_flux_pos_xwidth;
MCNUM mccxdist_flux_pos_yheight;
MCNUM mccxdist_flux_pos_restore_neutron;
int mccxdist_flux_pos_nowritefile;

/* Setting parameters for component 'PSD_fluxposB' [16]. */
int mccPSD_fluxposB_nx;
int mccPSD_fluxposB_ny;
char mccPSD_fluxposB_filename[16384];
MCNUM mccPSD_fluxposB_xmin;
MCNUM mccPSD_fluxposB_xmax;
MCNUM mccPSD_fluxposB_ymin;
MCNUM mccPSD_fluxposB_ymax;
MCNUM mccPSD_fluxposB_xwidth;
MCNUM mccPSD_fluxposB_yheight;
MCNUM mccPSD_fluxposB_restore_neutron;

/* Setting parameters for component 'window2' [17]. */
MCNUM mccwindow2_thickness;

/* Setting parameters for component 'in_slit' [18]. */
MCNUM mccin_slit_xmin;
MCNUM mccin_slit_xmax;
MCNUM mccin_slit_ymin;
MCNUM mccin_slit_ymax;
MCNUM mccin_slit_radius;
MCNUM mccin_slit_xwidth;
MCNUM mccin_slit_yheight;

/* Definition parameters for component 'lambda_in' [19]. */
#define mcclambda_in_nL 128
/* Setting parameters for component 'lambda_in' [19]. */
char mcclambda_in_filename[16384];
MCNUM mcclambda_in_xmin;
MCNUM mcclambda_in_xmax;
MCNUM mcclambda_in_ymin;
MCNUM mcclambda_in_ymax;
MCNUM mcclambda_in_xwidth;
MCNUM mcclambda_in_yheight;
MCNUM mcclambda_in_Lmin;
MCNUM mcclambda_in_Lmax;
MCNUM mcclambda_in_restore_neutron;
int mcclambda_in_nowritefile;

/* Setting parameters for component 'foc_mono' [21]. */
char mccfoc_mono_reflect[16384];
MCNUM mccfoc_mono_zwidth;
MCNUM mccfoc_mono_yheight;
MCNUM mccfoc_mono_gap;
MCNUM mccfoc_mono_NH;
MCNUM mccfoc_mono_NV;
MCNUM mccfoc_mono_mosaich;
MCNUM mccfoc_mono_mosaicv;
MCNUM mccfoc_mono_r0;
MCNUM mccfoc_mono_Q;
MCNUM mccfoc_mono_RV;
MCNUM mccfoc_mono_RH;
MCNUM mccfoc_mono_DM;
MCNUM mccfoc_mono_mosaic;
MCNUM mccfoc_mono_width;
MCNUM mccfoc_mono_height;
MCNUM mccfoc_mono_verbose;

/* Setting parameters for component 'out1_slit' [23]. */
MCNUM mccout1_slit_xmin;
MCNUM mccout1_slit_xmax;
MCNUM mccout1_slit_ymin;
MCNUM mccout1_slit_ymax;
MCNUM mccout1_slit_radius;
MCNUM mccout1_slit_xwidth;
MCNUM mccout1_slit_yheight;

/* Setting parameters for component 'Amoin_slit' [24]. */
MCNUM mccAmoin_slit_xmin;
MCNUM mccAmoin_slit_xmax;
MCNUM mccAmoin_slit_ymin;
MCNUM mccAmoin_slit_ymax;
MCNUM mccAmoin_slit_radius;
MCNUM mccAmoin_slit_xwidth;
MCNUM mccAmoin_slit_yheight;

/* Setting parameters for component 'Bmoin_slit' [25]. */
MCNUM mccBmoin_slit_xmin;
MCNUM mccBmoin_slit_xmax;
MCNUM mccBmoin_slit_ymin;
MCNUM mccBmoin_slit_ymax;
MCNUM mccBmoin_slit_radius;
MCNUM mccBmoin_slit_xwidth;
MCNUM mccBmoin_slit_yheight;

/* Setting parameters for component 'out2_slit' [26]. */
MCNUM mccout2_slit_xmin;
MCNUM mccout2_slit_xmax;
MCNUM mccout2_slit_ymin;
MCNUM mccout2_slit_ymax;
MCNUM mccout2_slit_radius;
MCNUM mccout2_slit_xwidth;
MCNUM mccout2_slit_yheight;

/* Setting parameters for component 'PSD_sample' [27]. */
int mccPSD_sample_nx;
int mccPSD_sample_ny;
char mccPSD_sample_filename[16384];
MCNUM mccPSD_sample_xmin;
MCNUM mccPSD_sample_xmax;
MCNUM mccPSD_sample_ymin;
MCNUM mccPSD_sample_ymax;
MCNUM mccPSD_sample_xwidth;
MCNUM mccPSD_sample_yheight;
MCNUM mccPSD_sample_restore_neutron;

/* Definition parameters for component 'lambda_sample' [28]. */
#define mcclambda_sample_nL 128
/* Setting parameters for component 'lambda_sample' [28]. */
char mcclambda_sample_filename[16384];
MCNUM mcclambda_sample_xmin;
MCNUM mcclambda_sample_xmax;
MCNUM mcclambda_sample_ymin;
MCNUM mcclambda_sample_ymax;
MCNUM mcclambda_sample_xwidth;
MCNUM mcclambda_sample_yheight;
MCNUM mcclambda_sample_Lmin;
MCNUM mcclambda_sample_Lmax;
MCNUM mcclambda_sample_restore_neutron;
int mcclambda_sample_nowritefile;

/* Definition parameters for component 'sample' [30]. */
#define mccsample_format Undefined
/* Setting parameters for component 'sample' [30]. */
char mccsample_reflections[16384];
char mccsample_geometry[16384];
MCNUM mccsample_radius;
MCNUM mccsample_yheight;
MCNUM mccsample_xwidth;
MCNUM mccsample_zdepth;
MCNUM mccsample_thickness;
MCNUM mccsample_pack;
MCNUM mccsample_Vc;
MCNUM mccsample_sigma_abs;
MCNUM mccsample_sigma_inc;
MCNUM mccsample_delta_d_d;
MCNUM mccsample_p_inc;
MCNUM mccsample_p_transmit;
MCNUM mccsample_DW;
MCNUM mccsample_nb_atoms;
MCNUM mccsample_d_omega;
MCNUM mccsample_d_phi;
MCNUM mccsample_tth_sign;
MCNUM mccsample_p_interact;
MCNUM mccsample_concentric;
MCNUM mccsample_density;
MCNUM mccsample_weight;
MCNUM mccsample_barns;
MCNUM mccsample_Strain;
MCNUM mccsample_focus_flip;
int mccsample_target_index;

/* Setting parameters for component 'STOP' [31]. */
MCNUM mccSTOP_xmin;
MCNUM mccSTOP_xmax;
MCNUM mccSTOP_ymin;
MCNUM mccSTOP_ymax;
MCNUM mccSTOP_xwidth;
MCNUM mccSTOP_yheight;
MCNUM mccSTOP_radius;

/* Definition parameters for component 'Detector' [32]. */
#define mccDetector_user1 FLT_MAX
#define mccDetector_user2 FLT_MAX
#define mccDetector_user3 FLT_MAX
/* Setting parameters for component 'Detector' [32]. */
MCNUM mccDetector_xwidth;
MCNUM mccDetector_yheight;
MCNUM mccDetector_zdepth;
MCNUM mccDetector_xmin;
MCNUM mccDetector_xmax;
MCNUM mccDetector_ymin;
MCNUM mccDetector_ymax;
MCNUM mccDetector_zmin;
MCNUM mccDetector_zmax;
MCNUM mccDetector_bins;
MCNUM mccDetector_min;
MCNUM mccDetector_max;
MCNUM mccDetector_restore_neutron;
MCNUM mccDetector_radius;
char mccDetector_options[16384];
char mccDetector_filename[16384];
char mccDetector_geometry[16384];
char mccDetector_username1[16384];
char mccDetector_username2[16384];
char mccDetector_username3[16384];
int mccDetector_nowritefile;

/* User component declarations. */

/* User declarations for component 'source_arm' [1]. */
#define mccompcurname  source_arm
#define mccompcurtype  Progress_bar
#define mccompcurindex 1
#define IntermediateCnts mccsource_arm_IntermediateCnts
#define StartTime mccsource_arm_StartTime
#define EndTime mccsource_arm_EndTime
#define CurrentTime mccsource_arm_CurrentTime
#define profile mccsource_arm_profile
#define percent mccsource_arm_percent
#define flag_save mccsource_arm_flag_save
#define minutes mccsource_arm_minutes
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
#line 10859 "./PSI_DMC.c"
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
#define mccompcurtype  Source_Maxwell_3
#define mccompcurindex 2
#define M mccsource_M
#define l_range mccsource_l_range
#define w_mult mccsource_w_mult
#define w_source mccsource_w_source
#define h_source mccsource_h_source
#define size mccsource_size
#define yheight mccsource_yheight
#define xwidth mccsource_xwidth
#define Lmin mccsource_Lmin
#define Lmax mccsource_Lmax
#define dist mccsource_dist
#define focus_xw mccsource_focus_xw
#define focus_yh mccsource_focus_yh
#define T1 mccsource_T1
#define T2 mccsource_T2
#define T3 mccsource_T3
#define I1 mccsource_I1
#define I2 mccsource_I2
#define I3 mccsource_I3
#define target_index mccsource_target_index
#define lambda0 mccsource_lambda0
#define dlambda mccsource_dlambda
#line 80 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../sources/Source_Maxwell_3.comp"
  double l_range, w_mult;
  double w_source, h_source;
#line 10901 "./PSI_DMC.c"
#undef dlambda
#undef lambda0
#undef target_index
#undef I3
#undef I2
#undef I1
#undef T3
#undef T2
#undef T1
#undef focus_yh
#undef focus_xw
#undef dist
#undef Lmax
#undef Lmin
#undef xwidth
#undef yheight
#undef size
#undef h_source
#undef w_source
#undef w_mult
#undef l_range
#undef M
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'PSDbefore_guides' [3]. */
#define mccompcurname  PSDbefore_guides
#define mccompcurtype  PSD_monitor
#define mccompcurindex 3
#define PSD_N mccPSDbefore_guides_PSD_N
#define PSD_p mccPSDbefore_guides_PSD_p
#define PSD_p2 mccPSDbefore_guides_PSD_p2
#define nx mccPSDbefore_guides_nx
#define ny mccPSDbefore_guides_ny
#define filename mccPSDbefore_guides_filename
#define xmin mccPSDbefore_guides_xmin
#define xmax mccPSDbefore_guides_xmax
#define ymin mccPSDbefore_guides_ymin
#define ymax mccPSDbefore_guides_ymax
#define xwidth mccPSDbefore_guides_xwidth
#define yheight mccPSDbefore_guides_yheight
#define restore_neutron mccPSDbefore_guides_restore_neutron
#line 62 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
  DArray2d PSD_N;
  DArray2d PSD_p;
  DArray2d PSD_p2;
#line 10949 "./PSI_DMC.c"
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

/* User declarations for component 'l_mon_source' [4]. */
#define mccompcurname  l_mon_source
#define mccompcurtype  L_monitor
#define mccompcurindex 4
#define nL mccl_mon_source_nL
#define L_N mccl_mon_source_L_N
#define L_p mccl_mon_source_L_p
#define L_p2 mccl_mon_source_L_p2
#define filename mccl_mon_source_filename
#define xmin mccl_mon_source_xmin
#define xmax mccl_mon_source_xmax
#define ymin mccl_mon_source_ymin
#define ymax mccl_mon_source_ymax
#define xwidth mccl_mon_source_xwidth
#define yheight mccl_mon_source_yheight
#define Lmin mccl_mon_source_Lmin
#define Lmax mccl_mon_source_Lmax
#define restore_neutron mccl_mon_source_restore_neutron
#define nowritefile mccl_mon_source_nowritefile
#line 57 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/L_monitor.comp"
double L_N[nL];
double L_p[nL], L_p2[nL];
#line 10989 "./PSI_DMC.c"
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

/* User declarations for component 'guide1' [5]. */
#define mccompcurname  guide1
#define mccompcurtype  Guide
#define mccompcurindex 5
#define pTable mccguide1_pTable
#define reflect mccguide1_reflect
#define w1 mccguide1_w1
#define h1 mccguide1_h1
#define w2 mccguide1_w2
#define h2 mccguide1_h2
#define l mccguide1_l
#define R0 mccguide1_R0
#define Qc mccguide1_Qc
#define alpha mccguide1_alpha
#define m mccguide1_m
#define W mccguide1_W
#line 69 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide.comp"
t_Table pTable;
#line 11027 "./PSI_DMC.c"
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

/* User declarations for component 'PSDbefore_curve' [6]. */
#define mccompcurname  PSDbefore_curve
#define mccompcurtype  PSD_monitor
#define mccompcurindex 6
#define PSD_N mccPSDbefore_curve_PSD_N
#define PSD_p mccPSDbefore_curve_PSD_p
#define PSD_p2 mccPSDbefore_curve_PSD_p2
#define nx mccPSDbefore_curve_nx
#define ny mccPSDbefore_curve_ny
#define filename mccPSDbefore_curve_filename
#define xmin mccPSDbefore_curve_xmin
#define xmax mccPSDbefore_curve_xmax
#define ymin mccPSDbefore_curve_ymin
#define ymax mccPSDbefore_curve_ymax
#define xwidth mccPSDbefore_curve_xwidth
#define yheight mccPSDbefore_curve_yheight
#define restore_neutron mccPSDbefore_curve_restore_neutron
#line 62 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
  DArray2d PSD_N;
  DArray2d PSD_p;
  DArray2d PSD_p2;
#line 11065 "./PSI_DMC.c"
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

/* User declarations for component 'guide2' [7]. */
#define mccompcurname  guide2
#define mccompcurtype  Bender
#define mccompcurindex 7
#define bk mccguide2_bk
#define mWin mccguide2_mWin
#define w mccguide2_w
#define h mccguide2_h
#define r mccguide2_r
#define Win mccguide2_Win
#define k mccguide2_k
#define d mccguide2_d
#define l mccguide2_l
#define R0a mccguide2_R0a
#define Qca mccguide2_Qca
#define alphaa mccguide2_alphaa
#define ma mccguide2_ma
#define Wa mccguide2_Wa
#define R0i mccguide2_R0i
#define Qci mccguide2_Qci
#define alphai mccguide2_alphai
#define mi mccguide2_mi
#define Wi mccguide2_Wi
#define R0s mccguide2_R0s
#define Qcs mccguide2_Qcs
#define alphas mccguide2_alphas
#define ms mccguide2_ms
#define Ws mccguide2_Ws
#line 108 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Bender.comp"
double bk, mWin;
#line 11113 "./PSI_DMC.c"
#undef Ws
#undef ms
#undef alphas
#undef Qcs
#undef R0s
#undef Wi
#undef mi
#undef alphai
#undef Qci
#undef R0i
#undef Wa
#undef ma
#undef alphaa
#undef Qca
#undef R0a
#undef l
#undef d
#undef k
#undef Win
#undef r
#undef h
#undef w
#undef mWin
#undef bk
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'PSDafter_curve' [8]. */
#define mccompcurname  PSDafter_curve
#define mccompcurtype  PSD_monitor
#define mccompcurindex 8
#define PSD_N mccPSDafter_curve_PSD_N
#define PSD_p mccPSDafter_curve_PSD_p
#define PSD_p2 mccPSDafter_curve_PSD_p2
#define nx mccPSDafter_curve_nx
#define ny mccPSDafter_curve_ny
#define filename mccPSDafter_curve_filename
#define xmin mccPSDafter_curve_xmin
#define xmax mccPSDafter_curve_xmax
#define ymin mccPSDafter_curve_ymin
#define ymax mccPSDafter_curve_ymax
#define xwidth mccPSDafter_curve_xwidth
#define yheight mccPSDafter_curve_yheight
#define restore_neutron mccPSDafter_curve_restore_neutron
#line 62 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
  DArray2d PSD_N;
  DArray2d PSD_p;
  DArray2d PSD_p2;
#line 11163 "./PSI_DMC.c"
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

/* User declarations for component 'bunker' [9]. */
#define mccompcurname  bunker
#define mccompcurtype  Guide
#define mccompcurindex 9
#define pTable mccbunker_pTable
#define reflect mccbunker_reflect
#define w1 mccbunker_w1
#define h1 mccbunker_h1
#define w2 mccbunker_w2
#define h2 mccbunker_h2
#define l mccbunker_l
#define R0 mccbunker_R0
#define Qc mccbunker_Qc
#define alpha mccbunker_alpha
#define m mccbunker_m
#define W mccbunker_W
#line 69 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide.comp"
t_Table pTable;
#line 11199 "./PSI_DMC.c"
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

/* User declarations for component 'guide3' [10]. */
#define mccompcurname  guide3
#define mccompcurtype  Guide
#define mccompcurindex 10
#define pTable mccguide3_pTable
#define reflect mccguide3_reflect
#define w1 mccguide3_w1
#define h1 mccguide3_h1
#define w2 mccguide3_w2
#define h2 mccguide3_h2
#define l mccguide3_l
#define R0 mccguide3_R0
#define Qc mccguide3_Qc
#define alpha mccguide3_alpha
#define m mccguide3_m
#define W mccguide3_W
#line 69 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide.comp"
t_Table pTable;
#line 11234 "./PSI_DMC.c"
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

/* User declarations for component 'guide4' [11]. */
#define mccompcurname  guide4
#define mccompcurtype  Guide
#define mccompcurindex 11
#define pTable mccguide4_pTable
#define reflect mccguide4_reflect
#define w1 mccguide4_w1
#define h1 mccguide4_h1
#define w2 mccguide4_w2
#define h2 mccguide4_h2
#define l mccguide4_l
#define R0 mccguide4_R0
#define Qc mccguide4_Qc
#define alpha mccguide4_alpha
#define m mccguide4_m
#define W mccguide4_W
#line 69 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide.comp"
t_Table pTable;
#line 11269 "./PSI_DMC.c"
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

/* User declarations for component 'ydist_fluxpos' [13]. */
#define mccompcurname  ydist_fluxpos
#define mccompcurtype  PSDlin_monitor
#define mccompcurindex 13
#define nx mccydist_fluxpos_nx
#define PSDlin_N mccydist_fluxpos_PSDlin_N
#define PSDlin_p mccydist_fluxpos_PSDlin_p
#define PSDlin_p2 mccydist_fluxpos_PSDlin_p2
#define filename mccydist_fluxpos_filename
#define xmin mccydist_fluxpos_xmin
#define xmax mccydist_fluxpos_xmax
#define ymin mccydist_fluxpos_ymin
#define ymax mccydist_fluxpos_ymax
#define xwidth mccydist_fluxpos_xwidth
#define yheight mccydist_fluxpos_yheight
#define restore_neutron mccydist_fluxpos_restore_neutron
#define nowritefile mccydist_fluxpos_nowritefile
#line 53 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/PSDlin_monitor.comp"
    double PSDlin_N[nx];
    double PSDlin_p[nx];
    double PSDlin_p2[nx];
#line 11307 "./PSI_DMC.c"
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

/* User declarations for component 'PSD_fluxpos' [14]. */
#define mccompcurname  PSD_fluxpos
#define mccompcurtype  PSD_monitor
#define mccompcurindex 14
#define PSD_N mccPSD_fluxpos_PSD_N
#define PSD_p mccPSD_fluxpos_PSD_p
#define PSD_p2 mccPSD_fluxpos_PSD_p2
#define nx mccPSD_fluxpos_nx
#define ny mccPSD_fluxpos_ny
#define filename mccPSD_fluxpos_filename
#define xmin mccPSD_fluxpos_xmin
#define xmax mccPSD_fluxpos_xmax
#define ymin mccPSD_fluxpos_ymin
#define ymax mccPSD_fluxpos_ymax
#define xwidth mccPSD_fluxpos_xwidth
#define yheight mccPSD_fluxpos_yheight
#define restore_neutron mccPSD_fluxpos_restore_neutron
#line 62 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
  DArray2d PSD_N;
  DArray2d PSD_p;
  DArray2d PSD_p2;
#line 11346 "./PSI_DMC.c"
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

/* User declarations for component 'xdist_flux_pos' [15]. */
#define mccompcurname  xdist_flux_pos
#define mccompcurtype  PSDlin_monitor
#define mccompcurindex 15
#define nx mccxdist_flux_pos_nx
#define PSDlin_N mccxdist_flux_pos_PSDlin_N
#define PSDlin_p mccxdist_flux_pos_PSDlin_p
#define PSDlin_p2 mccxdist_flux_pos_PSDlin_p2
#define filename mccxdist_flux_pos_filename
#define xmin mccxdist_flux_pos_xmin
#define xmax mccxdist_flux_pos_xmax
#define ymin mccxdist_flux_pos_ymin
#define ymax mccxdist_flux_pos_ymax
#define xwidth mccxdist_flux_pos_xwidth
#define yheight mccxdist_flux_pos_yheight
#define restore_neutron mccxdist_flux_pos_restore_neutron
#define nowritefile mccxdist_flux_pos_nowritefile
#line 53 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/PSDlin_monitor.comp"
    double PSDlin_N[nx];
    double PSDlin_p[nx];
    double PSDlin_p2[nx];
#line 11385 "./PSI_DMC.c"
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

/* User declarations for component 'PSD_fluxposB' [16]. */
#define mccompcurname  PSD_fluxposB
#define mccompcurtype  PSD_monitor
#define mccompcurindex 16
#define PSD_N mccPSD_fluxposB_PSD_N
#define PSD_p mccPSD_fluxposB_PSD_p
#define PSD_p2 mccPSD_fluxposB_PSD_p2
#define nx mccPSD_fluxposB_nx
#define ny mccPSD_fluxposB_ny
#define filename mccPSD_fluxposB_filename
#define xmin mccPSD_fluxposB_xmin
#define xmax mccPSD_fluxposB_xmax
#define ymin mccPSD_fluxposB_ymin
#define ymax mccPSD_fluxposB_ymax
#define xwidth mccPSD_fluxposB_xwidth
#define yheight mccPSD_fluxposB_yheight
#define restore_neutron mccPSD_fluxposB_restore_neutron
#line 62 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
  DArray2d PSD_N;
  DArray2d PSD_p;
  DArray2d PSD_p2;
#line 11424 "./PSI_DMC.c"
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

/* User declarations for component 'in_slit' [18]. */
#define mccompcurname  in_slit
#define mccompcurtype  Slit
#define mccompcurindex 18
#define xmin mccin_slit_xmin
#define xmax mccin_slit_xmax
#define ymin mccin_slit_ymin
#define ymax mccin_slit_ymax
#define radius mccin_slit_radius
#define xwidth mccin_slit_xwidth
#define yheight mccin_slit_yheight
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

/* User declarations for component 'lambda_in' [19]. */
#define mccompcurname  lambda_in
#define mccompcurtype  L_monitor
#define mccompcurindex 19
#define nL mcclambda_in_nL
#define L_N mcclambda_in_L_N
#define L_p mcclambda_in_L_p
#define L_p2 mcclambda_in_L_p2
#define filename mcclambda_in_filename
#define xmin mcclambda_in_xmin
#define xmax mcclambda_in_xmax
#define ymin mcclambda_in_ymin
#define ymax mcclambda_in_ymax
#define xwidth mcclambda_in_xwidth
#define yheight mcclambda_in_yheight
#define Lmin mcclambda_in_Lmin
#define Lmax mcclambda_in_Lmax
#define restore_neutron mcclambda_in_restore_neutron
#define nowritefile mcclambda_in_nowritefile
#line 57 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/L_monitor.comp"
double L_N[nL];
double L_p[nL], L_p2[nL];
#line 11486 "./PSI_DMC.c"
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

/* User declarations for component 'sma' [20]. */
#define mccompcurname  sma
#define mccompcurtype  Arm
#define mccompcurindex 20
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'foc_mono' [21]. */
#define mccompcurname  foc_mono
#define mccompcurtype  Monochromator_2foc
#define mccompcurindex 21
#define mos_y mccfoc_mono_mos_y
#define mos_z mccfoc_mono_mos_z
#define mono_Q mccfoc_mono_mono_Q
#define SlabWidth mccfoc_mono_SlabWidth
#define SlabHeight mccfoc_mono_SlabHeight
#define rTable mccfoc_mono_rTable
#define reflect mccfoc_mono_reflect
#define zwidth mccfoc_mono_zwidth
#define yheight mccfoc_mono_yheight
#define gap mccfoc_mono_gap
#define NH mccfoc_mono_NH
#define NV mccfoc_mono_NV
#define mosaich mccfoc_mono_mosaich
#define mosaicv mccfoc_mono_mosaicv
#define r0 mccfoc_mono_r0
#define Q mccfoc_mono_Q
#define RV mccfoc_mono_RV
#define RH mccfoc_mono_RH
#define DM mccfoc_mono_DM
#define mosaic mccfoc_mono_mosaic
#define width mccfoc_mono_width
#define height mccfoc_mono_height
#define verbose mccfoc_mono_verbose
#line 88 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../contrib/Monochromator_2foc.comp"
#ifndef DIV_CUTOFF
#define DIV_CUTOFF 2            /* ~ 10^-5 cutoff. */
#endif
double mos_y; /* mosaic, in radians */
double mos_z;
double mono_Q;
double SlabWidth, SlabHeight;
t_Table rTable;
#line 11550 "./PSI_DMC.c"
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

/* User declarations for component 'msa' [22]. */
#define mccompcurname  msa
#define mccompcurtype  Arm
#define mccompcurindex 22
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'out1_slit' [23]. */
#define mccompcurname  out1_slit
#define mccompcurtype  Slit
#define mccompcurindex 23
#define xmin mccout1_slit_xmin
#define xmax mccout1_slit_xmax
#define ymin mccout1_slit_ymin
#define ymax mccout1_slit_ymax
#define radius mccout1_slit_radius
#define xwidth mccout1_slit_xwidth
#define yheight mccout1_slit_yheight
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

/* User declarations for component 'Amoin_slit' [24]. */
#define mccompcurname  Amoin_slit
#define mccompcurtype  Slit
#define mccompcurindex 24
#define xmin mccAmoin_slit_xmin
#define xmax mccAmoin_slit_xmax
#define ymin mccAmoin_slit_ymin
#define ymax mccAmoin_slit_ymax
#define radius mccAmoin_slit_radius
#define xwidth mccAmoin_slit_xwidth
#define yheight mccAmoin_slit_yheight
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

/* User declarations for component 'Bmoin_slit' [25]. */
#define mccompcurname  Bmoin_slit
#define mccompcurtype  Slit
#define mccompcurindex 25
#define xmin mccBmoin_slit_xmin
#define xmax mccBmoin_slit_xmax
#define ymin mccBmoin_slit_ymin
#define ymax mccBmoin_slit_ymax
#define radius mccBmoin_slit_radius
#define xwidth mccBmoin_slit_xwidth
#define yheight mccBmoin_slit_yheight
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

/* User declarations for component 'out2_slit' [26]. */
#define mccompcurname  out2_slit
#define mccompcurtype  Slit
#define mccompcurindex 26
#define xmin mccout2_slit_xmin
#define xmax mccout2_slit_xmax
#define ymin mccout2_slit_ymin
#define ymax mccout2_slit_ymax
#define radius mccout2_slit_radius
#define xwidth mccout2_slit_xwidth
#define yheight mccout2_slit_yheight
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

/* User declarations for component 'PSD_sample' [27]. */
#define mccompcurname  PSD_sample
#define mccompcurtype  PSD_monitor
#define mccompcurindex 27
#define PSD_N mccPSD_sample_PSD_N
#define PSD_p mccPSD_sample_PSD_p
#define PSD_p2 mccPSD_sample_PSD_p2
#define nx mccPSD_sample_nx
#define ny mccPSD_sample_ny
#define filename mccPSD_sample_filename
#define xmin mccPSD_sample_xmin
#define xmax mccPSD_sample_xmax
#define ymin mccPSD_sample_ymin
#define ymax mccPSD_sample_ymax
#define xwidth mccPSD_sample_xwidth
#define yheight mccPSD_sample_yheight
#define restore_neutron mccPSD_sample_restore_neutron
#line 62 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
  DArray2d PSD_N;
  DArray2d PSD_p;
  DArray2d PSD_p2;
#line 11695 "./PSI_DMC.c"
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

/* User declarations for component 'lambda_sample' [28]. */
#define mccompcurname  lambda_sample
#define mccompcurtype  L_monitor
#define mccompcurindex 28
#define nL mcclambda_sample_nL
#define L_N mcclambda_sample_L_N
#define L_p mcclambda_sample_L_p
#define L_p2 mcclambda_sample_L_p2
#define filename mcclambda_sample_filename
#define xmin mcclambda_sample_xmin
#define xmax mcclambda_sample_xmax
#define ymin mcclambda_sample_ymin
#define ymax mcclambda_sample_ymax
#define xwidth mcclambda_sample_xwidth
#define yheight mcclambda_sample_yheight
#define Lmin mcclambda_sample_Lmin
#define Lmax mcclambda_sample_Lmax
#define restore_neutron mcclambda_sample_restore_neutron
#define nowritefile mcclambda_sample_nowritefile
#line 57 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/L_monitor.comp"
double L_N[nL];
double L_p[nL], L_p2[nL];
#line 11735 "./PSI_DMC.c"
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

/* User declarations for component 'sa_arm' [29]. */
#define mccompcurname  sa_arm
#define mccompcurtype  Arm
#define mccompcurindex 29
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'sample' [30]. */
#define mccompcurname  sample
#define mccompcurtype  PowderN
#define mccompcurindex 30
#define format mccsample_format
#define line_info mccsample_line_info
#define columns mccsample_columns
#define offdata mccsample_offdata
#define reflections mccsample_reflections
#define geometry mccsample_geometry
#define radius mccsample_radius
#define yheight mccsample_yheight
#define xwidth mccsample_xwidth
#define zdepth mccsample_zdepth
#define thickness mccsample_thickness
#define pack mccsample_pack
#define Vc mccsample_Vc
#define sigma_abs mccsample_sigma_abs
#define sigma_inc mccsample_sigma_inc
#define delta_d_d mccsample_delta_d_d
#define p_inc mccsample_p_inc
#define p_transmit mccsample_p_transmit
#define DW mccsample_DW
#define nb_atoms mccsample_nb_atoms
#define d_omega mccsample_d_omega
#define d_phi mccsample_d_phi
#define tth_sign mccsample_tth_sign
#define p_interact mccsample_p_interact
#define concentric mccsample_concentric
#define density mccsample_density
#define weight mccsample_weight
#define barns mccsample_barns
#define Strain mccsample_Strain
#define focus_flip mccsample_focus_flip
#define target_index mccsample_target_index
#line 552 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../samples/PowderN.comp"
  struct line_info_struct line_info;
  int *columns;
  off_struct offdata;
  double tgt_x, tgt_y, tgt_z;

#line 11804 "./PSI_DMC.c"
#undef target_index
#undef focus_flip
#undef Strain
#undef barns
#undef weight
#undef density
#undef concentric
#undef p_interact
#undef tth_sign
#undef d_phi
#undef d_omega
#undef nb_atoms
#undef DW
#undef p_transmit
#undef p_inc
#undef delta_d_d
#undef sigma_inc
#undef sigma_abs
#undef Vc
#undef pack
#undef thickness
#undef zdepth
#undef xwidth
#undef yheight
#undef radius
#undef geometry
#undef reflections
#undef offdata
#undef columns
#undef line_info
#undef format
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'STOP' [31]. */
#define mccompcurname  STOP
#define mccompcurtype  Beamstop
#define mccompcurindex 31
#define xmin mccSTOP_xmin
#define xmax mccSTOP_xmax
#define ymin mccSTOP_ymin
#define ymax mccSTOP_ymax
#define xwidth mccSTOP_xwidth
#define yheight mccSTOP_yheight
#define radius mccSTOP_radius
#undef radius
#undef yheight
#undef xwidth
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'Detector' [32]. */
#define mccompcurname  Detector
#define mccompcurtype  Monitor_nD
#define mccompcurindex 32
#define user1 mccDetector_user1
#define user2 mccDetector_user2
#define user3 mccDetector_user3
#define DEFS mccDetector_DEFS
#define Vars mccDetector_Vars
#define detector mccDetector_detector
#define offdata mccDetector_offdata
#define xwidth mccDetector_xwidth
#define yheight mccDetector_yheight
#define zdepth mccDetector_zdepth
#define xmin mccDetector_xmin
#define xmax mccDetector_xmax
#define ymin mccDetector_ymin
#define ymax mccDetector_ymax
#define zmin mccDetector_zmin
#define zmax mccDetector_zmax
#define bins mccDetector_bins
#define min mccDetector_min
#define max mccDetector_max
#define restore_neutron mccDetector_restore_neutron
#define radius mccDetector_radius
#define options mccDetector_options
#define filename mccDetector_filename
#define geometry mccDetector_geometry
#define username1 mccDetector_username1
#define username2 mccDetector_username2
#define username3 mccDetector_username3
#define nowritefile mccDetector_nowritefile
#line 225 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/Monitor_nD.comp"
  MonitornD_Defines_type DEFS;
  MonitornD_Variables_type Vars;
  MCDETECTOR detector;
  off_struct offdata;
#line 11899 "./PSI_DMC.c"
#undef nowritefile
#undef username3
#undef username2
#undef username1
#undef geometry
#undef filename
#undef options
#undef radius
#undef restore_neutron
#undef max
#undef min
#undef bins
#undef zmax
#undef zmin
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef zdepth
#undef yheight
#undef xwidth
#undef offdata
#undef detector
#undef Vars
#undef DEFS
#undef user3
#undef user2
#undef user1
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

Coords mcposasource_arm, mcposrsource_arm;
Rotation mcrotasource_arm, mcrotrsource_arm;
Coords mcposasource, mcposrsource;
Rotation mcrotasource, mcrotrsource;
Coords mcposaPSDbefore_guides, mcposrPSDbefore_guides;
Rotation mcrotaPSDbefore_guides, mcrotrPSDbefore_guides;
Coords mcposal_mon_source, mcposrl_mon_source;
Rotation mcrotal_mon_source, mcrotrl_mon_source;
Coords mcposaguide1, mcposrguide1;
Rotation mcrotaguide1, mcrotrguide1;
Coords mcposaPSDbefore_curve, mcposrPSDbefore_curve;
Rotation mcrotaPSDbefore_curve, mcrotrPSDbefore_curve;
Coords mcposaguide2, mcposrguide2;
Rotation mcrotaguide2, mcrotrguide2;
Coords mcposaPSDafter_curve, mcposrPSDafter_curve;
Rotation mcrotaPSDafter_curve, mcrotrPSDafter_curve;
Coords mcposabunker, mcposrbunker;
Rotation mcrotabunker, mcrotrbunker;
Coords mcposaguide3, mcposrguide3;
Rotation mcrotaguide3, mcrotrguide3;
Coords mcposaguide4, mcposrguide4;
Rotation mcrotaguide4, mcrotrguide4;
Coords mcposawindow1, mcposrwindow1;
Rotation mcrotawindow1, mcrotrwindow1;
Coords mcposaydist_fluxpos, mcposrydist_fluxpos;
Rotation mcrotaydist_fluxpos, mcrotrydist_fluxpos;
Coords mcposaPSD_fluxpos, mcposrPSD_fluxpos;
Rotation mcrotaPSD_fluxpos, mcrotrPSD_fluxpos;
Coords mcposaxdist_flux_pos, mcposrxdist_flux_pos;
Rotation mcrotaxdist_flux_pos, mcrotrxdist_flux_pos;
Coords mcposaPSD_fluxposB, mcposrPSD_fluxposB;
Rotation mcrotaPSD_fluxposB, mcrotrPSD_fluxposB;
Coords mcposawindow2, mcposrwindow2;
Rotation mcrotawindow2, mcrotrwindow2;
Coords mcposain_slit, mcposrin_slit;
Rotation mcrotain_slit, mcrotrin_slit;
Coords mcposalambda_in, mcposrlambda_in;
Rotation mcrotalambda_in, mcrotrlambda_in;
Coords mcposasma, mcposrsma;
Rotation mcrotasma, mcrotrsma;
Coords mcposafoc_mono, mcposrfoc_mono;
Rotation mcrotafoc_mono, mcrotrfoc_mono;
Coords mcposamsa, mcposrmsa;
Rotation mcrotamsa, mcrotrmsa;
Coords mcposaout1_slit, mcposrout1_slit;
Rotation mcrotaout1_slit, mcrotrout1_slit;
Coords mcposaAmoin_slit, mcposrAmoin_slit;
Rotation mcrotaAmoin_slit, mcrotrAmoin_slit;
Coords mcposaBmoin_slit, mcposrBmoin_slit;
Rotation mcrotaBmoin_slit, mcrotrBmoin_slit;
Coords mcposaout2_slit, mcposrout2_slit;
Rotation mcrotaout2_slit, mcrotrout2_slit;
Coords mcposaPSD_sample, mcposrPSD_sample;
Rotation mcrotaPSD_sample, mcrotrPSD_sample;
Coords mcposalambda_sample, mcposrlambda_sample;
Rotation mcrotalambda_sample, mcrotrlambda_sample;
Coords mcposasa_arm, mcposrsa_arm;
Rotation mcrotasa_arm, mcrotrsa_arm;
Coords mcposasample, mcposrsample;
Rotation mcrotasample, mcrotrsample;
Coords mcposaSTOP, mcposrSTOP;
Rotation mcrotaSTOP, mcrotrSTOP;
Coords mcposaDetector, mcposrDetector;
Rotation mcrotaDetector, mcrotrDetector;

MCNUM mcnx, mcny, mcnz, mcnvx, mcnvy, mcnvz, mcnt, mcnsx, mcnsy, mcnsz, mcnp;

/* end declare */

void mcinit(void) {
#define mccompcurname  PSI_DMC
#define mccompcurtype  INSTRUMENT
#define mccompcurindex 0
#define mcposaPSI_DMC coords_set(0,0,0)
#define lambda mciplambda
#define R mcipR
#define R_curve mcipR_curve
#define filename mcipfilename
#define D_PHI mcipD_PHI
#define SHIFT mcipSHIFT
#define PACK mcipPACK
#define Dw mcipDw
#define BARNS mcipBARNS
#line 80 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
{
  TTM = 2*asin(mono_q*lambda/(4*PI))*RAD2DEG;
  OMA = TTM/2;
  RV = fabs(2*2.82*sin(DEG2RAD*OMA));
  
  angleGuideCurved=-2.0*asin(0.4995 /2.0/3612)/PI*180;
  alpha=(R0-R)/Qc/(Mvalue-1);
  alpha_curve=(R0_curve-R_curve)/Qc_curve/(Mvalue_curve-1);
  
}
#line 12026 "./PSI_DMC.c"
#undef BARNS
#undef Dw
#undef PACK
#undef SHIFT
#undef D_PHI
#undef filename
#undef R_curve
#undef R
#undef lambda
#undef mcposaPSI_DMC
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
    /* Component source_arm. */
  /* Setting parameters for component source_arm. */
  SIG_MESSAGE("source_arm (Init:SetPar)");
#line 39 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  if("NULL") strncpy(mccsource_arm_profile, "NULL" ? "NULL" : "", 16384); else mccsource_arm_profile[0]='\0';
#line 39 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccsource_arm_percent = 10;
#line 39 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccsource_arm_flag_save = 0;
#line 39 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccsource_arm_minutes = 0;
#line 12061 "./PSI_DMC.c"

  SIG_MESSAGE("source_arm (Init:Place/Rotate)");
  rot_set_rotation(mcrotasource_arm,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 12068 "./PSI_DMC.c"
  rot_copy(mcrotrsource_arm, mcrotasource_arm);
  mcposasource_arm = coords_set(
#line 94 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    0,
#line 94 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    0,
#line 94 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    0);
#line 12077 "./PSI_DMC.c"
  mctc1 = coords_neg(mcposasource_arm);
  mcposrsource_arm = rot_apply(mcrotasource_arm, mctc1);
  mcDEBUG_COMPONENT("source_arm", mcposasource_arm, mcrotasource_arm)
  mccomp_posa[1] = mcposasource_arm;
  mccomp_posr[1] = mcposrsource_arm;
  mcNCounter[1]  = mcPCounter[1] = mcP2Counter[1] = 0;
  mcAbsorbProp[1]= 0;
    /* Component source. */
  /* Setting parameters for component source. */
  SIG_MESSAGE("source (Init:SetPar)");
#line 60 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccsource_size = 0;
#line 97 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccsource_yheight = 0.156;
#line 97 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccsource_xwidth = 0.126;
#line 98 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccsource_Lmin = mciplambda - ldiff / 2;
#line 98 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccsource_Lmax = mciplambda + ldiff / 2;
#line 99 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccsource_dist = 1.5;
#line 99 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccsource_focus_xw = 0.02;
#line 99 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccsource_focus_yh = 0.12;
#line 100 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccsource_T1 = 296.16;
#line 101 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccsource_T2 = 40.68;
#line 62 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccsource_T3 = 300;
#line 100 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccsource_I1 = 8.5E11;
#line 101 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccsource_I2 = 5.2E11;
#line 62 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccsource_I3 = 0;
#line 63 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccsource_target_index = + 1;
#line 63 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccsource_lambda0 = 0;
#line 63 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccsource_dlambda = 0;
#line 12122 "./PSI_DMC.c"

  SIG_MESSAGE("source (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 102 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    (0)*DEG2RAD,
#line 102 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    (0)*DEG2RAD,
#line 102 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    (0)*DEG2RAD);
#line 12132 "./PSI_DMC.c"
  rot_mul(mctr1, mcrotasource_arm, mcrotasource);
  rot_transpose(mcrotasource_arm, mctr1);
  rot_mul(mcrotasource, mctr1, mcrotrsource);
  mctc1 = coords_set(
#line 102 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    0,
#line 102 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    0,
#line 102 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    0);
#line 12143 "./PSI_DMC.c"
  rot_transpose(mcrotasource_arm, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposasource = coords_add(mcposasource_arm, mctc2);
  mctc1 = coords_sub(mcposasource_arm, mcposasource);
  mcposrsource = rot_apply(mcrotasource, mctc1);
  mcDEBUG_COMPONENT("source", mcposasource, mcrotasource)
  mccomp_posa[2] = mcposasource;
  mccomp_posr[2] = mcposrsource;
  mcNCounter[2]  = mcPCounter[2] = mcP2Counter[2] = 0;
  mcAbsorbProp[2]= 0;
    /* Component PSDbefore_guides. */
  /* Setting parameters for component PSDbefore_guides. */
  SIG_MESSAGE("PSDbefore_guides (Init:SetPar)");
#line 105 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccPSDbefore_guides_nx = 128;
#line 105 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccPSDbefore_guides_ny = 128;
#line 105 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  if("PSDbefore_guides") strncpy(mccPSDbefore_guides_filename, "PSDbefore_guides" ? "PSDbefore_guides" : "", 16384); else mccPSDbefore_guides_filename[0]='\0';
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccPSDbefore_guides_xmin = -0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccPSDbefore_guides_xmax = 0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccPSDbefore_guides_ymin = -0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccPSDbefore_guides_ymax = 0.05;
#line 106 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccPSDbefore_guides_xwidth = 0.02;
#line 106 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccPSDbefore_guides_yheight = 0.12;
#line 51 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccPSDbefore_guides_restore_neutron = 0;
#line 12177 "./PSI_DMC.c"

  SIG_MESSAGE("PSDbefore_guides (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 12184 "./PSI_DMC.c"
  rot_mul(mctr1, mcrotasource_arm, mcrotaPSDbefore_guides);
  rot_transpose(mcrotasource, mctr1);
  rot_mul(mcrotaPSDbefore_guides, mctr1, mcrotrPSDbefore_guides);
  mctc1 = coords_set(
#line 107 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    0,
#line 107 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    0,
#line 107 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    1.49999);
#line 12195 "./PSI_DMC.c"
  rot_transpose(mcrotasource_arm, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaPSDbefore_guides = coords_add(mcposasource_arm, mctc2);
  mctc1 = coords_sub(mcposasource, mcposaPSDbefore_guides);
  mcposrPSDbefore_guides = rot_apply(mcrotaPSDbefore_guides, mctc1);
  mcDEBUG_COMPONENT("PSDbefore_guides", mcposaPSDbefore_guides, mcrotaPSDbefore_guides)
  mccomp_posa[3] = mcposaPSDbefore_guides;
  mccomp_posr[3] = mcposrPSDbefore_guides;
  mcNCounter[3]  = mcPCounter[3] = mcP2Counter[3] = 0;
  mcAbsorbProp[3]= 0;
    /* Component l_mon_source. */
  /* Setting parameters for component l_mon_source. */
  SIG_MESSAGE("l_mon_source (Init:SetPar)");
#line 110 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  if("lmonsource.dat") strncpy(mccl_mon_source_filename, "lmonsource.dat" ? "lmonsource.dat" : "", 16384); else mccl_mon_source_filename[0]='\0';
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccl_mon_source_xmin = -0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccl_mon_source_xmax = 0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccl_mon_source_ymin = -0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccl_mon_source_ymax = 0.05;
#line 110 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccl_mon_source_xwidth = 0.02;
#line 111 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccl_mon_source_yheight = 0.12;
#line 111 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccl_mon_source_Lmin = 0;
#line 111 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccl_mon_source_Lmax = 20;
#line 51 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccl_mon_source_restore_neutron = 0;
#line 51 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccl_mon_source_nowritefile = 0;
#line 12231 "./PSI_DMC.c"

  SIG_MESSAGE("l_mon_source (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 12238 "./PSI_DMC.c"
  rot_mul(mctr1, mcrotaPSDbefore_guides, mcrotal_mon_source);
  rot_transpose(mcrotaPSDbefore_guides, mctr1);
  rot_mul(mcrotal_mon_source, mctr1, mcrotrl_mon_source);
  mctc1 = coords_set(
#line 112 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    0,
#line 112 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    0,
#line 112 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    1e-9);
#line 12249 "./PSI_DMC.c"
  rot_transpose(mcrotaPSDbefore_guides, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposal_mon_source = coords_add(mcposaPSDbefore_guides, mctc2);
  mctc1 = coords_sub(mcposaPSDbefore_guides, mcposal_mon_source);
  mcposrl_mon_source = rot_apply(mcrotal_mon_source, mctc1);
  mcDEBUG_COMPONENT("l_mon_source", mcposal_mon_source, mcrotal_mon_source)
  mccomp_posa[4] = mcposal_mon_source;
  mccomp_posr[4] = mcposrl_mon_source;
  mcNCounter[4]  = mcPCounter[4] = mcP2Counter[4] = 0;
  mcAbsorbProp[4]= 0;
    /* Component guide1. */
  /* Setting parameters for component guide1. */
  SIG_MESSAGE("guide1 (Init:SetPar)");
#line 58 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  if(0) strncpy(mccguide1_reflect, 0 ? 0 : "", 16384); else mccguide1_reflect[0]='\0';
#line 116 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccguide1_w1 = 0.02;
#line 116 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccguide1_h1 = 0.12;
#line 116 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccguide1_w2 = 0.02;
#line 116 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccguide1_h2 = 0.12;
#line 117 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccguide1_l = 4.66;
#line 117 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccguide1_R0 = R0;
#line 117 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccguide1_Qc = Qc;
#line 117 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccguide1_alpha = alpha;
#line 118 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccguide1_m = 1.8;
#line 118 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccguide1_W = W;
#line 12285 "./PSI_DMC.c"

  SIG_MESSAGE("guide1 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 119 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    (0)*DEG2RAD,
#line 119 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    (0)*DEG2RAD,
#line 119 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    (0)*DEG2RAD);
#line 12295 "./PSI_DMC.c"
  rot_mul(mctr1, mcrotasource_arm, mcrotaguide1);
  rot_transpose(mcrotal_mon_source, mctr1);
  rot_mul(mcrotaguide1, mctr1, mcrotrguide1);
  mctc1 = coords_set(
#line 119 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    0,
#line 119 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    0,
#line 119 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    1.50);
#line 12306 "./PSI_DMC.c"
  rot_transpose(mcrotasource_arm, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaguide1 = coords_add(mcposasource_arm, mctc2);
  mctc1 = coords_sub(mcposal_mon_source, mcposaguide1);
  mcposrguide1 = rot_apply(mcrotaguide1, mctc1);
  mcDEBUG_COMPONENT("guide1", mcposaguide1, mcrotaguide1)
  mccomp_posa[5] = mcposaguide1;
  mccomp_posr[5] = mcposrguide1;
  mcNCounter[5]  = mcPCounter[5] = mcP2Counter[5] = 0;
  mcAbsorbProp[5]= 0;
    /* Component PSDbefore_curve. */
  /* Setting parameters for component PSDbefore_curve. */
  SIG_MESSAGE("PSDbefore_curve (Init:SetPar)");
#line 122 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccPSDbefore_curve_nx = 128;
#line 122 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccPSDbefore_curve_ny = 128;
#line 122 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  if("PSDbefore_curve") strncpy(mccPSDbefore_curve_filename, "PSDbefore_curve" ? "PSDbefore_curve" : "", 16384); else mccPSDbefore_curve_filename[0]='\0';
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccPSDbefore_curve_xmin = -0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccPSDbefore_curve_xmax = 0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccPSDbefore_curve_ymin = -0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccPSDbefore_curve_ymax = 0.05;
#line 123 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccPSDbefore_curve_xwidth = 0.02;
#line 123 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccPSDbefore_curve_yheight = 0.12;
#line 51 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccPSDbefore_curve_restore_neutron = 0;
#line 12340 "./PSI_DMC.c"

  SIG_MESSAGE("PSDbefore_curve (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 12347 "./PSI_DMC.c"
  rot_mul(mctr1, mcrotaguide1, mcrotaPSDbefore_curve);
  rot_transpose(mcrotaguide1, mctr1);
  rot_mul(mcrotaPSDbefore_curve, mctr1, mcrotrPSDbefore_curve);
  mctc1 = coords_set(
#line 124 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    0,
#line 124 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    0,
#line 124 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    4.664);
#line 12358 "./PSI_DMC.c"
  rot_transpose(mcrotaguide1, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaPSDbefore_curve = coords_add(mcposaguide1, mctc2);
  mctc1 = coords_sub(mcposaguide1, mcposaPSDbefore_curve);
  mcposrPSDbefore_curve = rot_apply(mcrotaPSDbefore_curve, mctc1);
  mcDEBUG_COMPONENT("PSDbefore_curve", mcposaPSDbefore_curve, mcrotaPSDbefore_curve)
  mccomp_posa[6] = mcposaPSDbefore_curve;
  mccomp_posr[6] = mcposrPSDbefore_curve;
  mcNCounter[6]  = mcPCounter[6] = mcP2Counter[6] = 0;
  mcAbsorbProp[6]= 0;
    /* Component guide2. */
  /* Setting parameters for component guide2. */
  SIG_MESSAGE("guide2 (Init:SetPar)");
#line 127 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccguide2_w = 0.02;
#line 127 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccguide2_h = 0.12;
#line 127 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccguide2_r = 3612;
#line 94 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccguide2_Win = 0.04;
#line 94 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccguide2_k = 1;
#line 94 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccguide2_d = 0.001;
#line 130 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccguide2_l = 20;
#line 127 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccguide2_R0a = R0_curve;
#line 127 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccguide2_Qca = Qc_curve;
#line 128 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccguide2_alphaa = alpha_curve;
#line 128 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccguide2_ma = Mvalue_curve;
#line 128 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccguide2_Wa = W_curve;
#line 128 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccguide2_R0i = R0_curve;
#line 128 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccguide2_Qci = Qc_curve;
#line 129 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccguide2_alphai = alpha_curve;
#line 129 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccguide2_mi = 1;
#line 129 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccguide2_Wi = W_curve;
#line 129 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccguide2_R0s = R0_curve;
#line 129 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccguide2_Qcs = Qc_curve;
#line 130 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccguide2_alphas = alpha_curve;
#line 130 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccguide2_ms = Mvalue_curve;
#line 130 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccguide2_Ws = W_curve;
#line 12416 "./PSI_DMC.c"

  SIG_MESSAGE("guide2 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 12423 "./PSI_DMC.c"
  rot_mul(mctr1, mcrotaguide1, mcrotaguide2);
  rot_transpose(mcrotaPSDbefore_curve, mctr1);
  rot_mul(mcrotaguide2, mctr1, mcrotrguide2);
  mctc1 = coords_set(
#line 131 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    0,
#line 131 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    0,
#line 131 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    4.69);
#line 12434 "./PSI_DMC.c"
  rot_transpose(mcrotaguide1, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaguide2 = coords_add(mcposaguide1, mctc2);
  mctc1 = coords_sub(mcposaPSDbefore_curve, mcposaguide2);
  mcposrguide2 = rot_apply(mcrotaguide2, mctc1);
  mcDEBUG_COMPONENT("guide2", mcposaguide2, mcrotaguide2)
  mccomp_posa[7] = mcposaguide2;
  mccomp_posr[7] = mcposrguide2;
  mcNCounter[7]  = mcPCounter[7] = mcP2Counter[7] = 0;
  mcAbsorbProp[7]= 0;
    /* Component PSDafter_curve. */
  /* Setting parameters for component PSDafter_curve. */
  SIG_MESSAGE("PSDafter_curve (Init:SetPar)");
#line 134 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccPSDafter_curve_nx = 128;
#line 134 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccPSDafter_curve_ny = 128;
#line 134 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  if("PSDafter_curve") strncpy(mccPSDafter_curve_filename, "PSDafter_curve" ? "PSDafter_curve" : "", 16384); else mccPSDafter_curve_filename[0]='\0';
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccPSDafter_curve_xmin = -0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccPSDafter_curve_xmax = 0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccPSDafter_curve_ymin = -0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccPSDafter_curve_ymax = 0.05;
#line 135 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccPSDafter_curve_xwidth = 0.02;
#line 135 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccPSDafter_curve_yheight = 0.12;
#line 51 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccPSDafter_curve_restore_neutron = 0;
#line 12468 "./PSI_DMC.c"

  SIG_MESSAGE("PSDafter_curve (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 12475 "./PSI_DMC.c"
  rot_mul(mctr1, mcrotaguide2, mcrotaPSDafter_curve);
  rot_transpose(mcrotaguide2, mctr1);
  rot_mul(mcrotaPSDafter_curve, mctr1, mcrotrPSDafter_curve);
  mctc1 = coords_set(
#line 136 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    0,
#line 136 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    0,
#line 136 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    20.0001);
#line 12486 "./PSI_DMC.c"
  rot_transpose(mcrotaguide2, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaPSDafter_curve = coords_add(mcposaguide2, mctc2);
  mctc1 = coords_sub(mcposaguide2, mcposaPSDafter_curve);
  mcposrPSDafter_curve = rot_apply(mcrotaPSDafter_curve, mctc1);
  mcDEBUG_COMPONENT("PSDafter_curve", mcposaPSDafter_curve, mcrotaPSDafter_curve)
  mccomp_posa[8] = mcposaPSDafter_curve;
  mccomp_posr[8] = mcposrPSDafter_curve;
  mcNCounter[8]  = mcPCounter[8] = mcP2Counter[8] = 0;
  mcAbsorbProp[8]= 0;
    /* Component bunker. */
  /* Setting parameters for component bunker. */
  SIG_MESSAGE("bunker (Init:SetPar)");
#line 58 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  if(0) strncpy(mccbunker_reflect, 0 ? 0 : "", 16384); else mccbunker_reflect[0]='\0';
#line 140 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccbunker_w1 = 0.02;
#line 140 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccbunker_h1 = .12;
#line 140 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccbunker_w2 = 0.02;
#line 140 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccbunker_h2 = .12;
#line 141 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccbunker_l = 3.43;
#line 141 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccbunker_R0 = R0;
#line 141 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccbunker_Qc = Qc;
#line 141 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccbunker_alpha = alpha;
#line 142 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccbunker_m = 1.6;
#line 142 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccbunker_W = W;
#line 12522 "./PSI_DMC.c"

  SIG_MESSAGE("bunker (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 143 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    (0)*DEG2RAD,
#line 143 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    (0)*DEG2RAD,
#line 143 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    (0)*DEG2RAD);
#line 12532 "./PSI_DMC.c"
  rot_mul(mctr1, mcrotaguide2, mcrotabunker);
  rot_transpose(mcrotaPSDafter_curve, mctr1);
  rot_mul(mcrotabunker, mctr1, mcrotrbunker);
  mctc1 = coords_set(
#line 143 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    0,
#line 143 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    0,
#line 143 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    20.1502);
#line 12543 "./PSI_DMC.c"
  rot_transpose(mcrotaguide2, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposabunker = coords_add(mcposaguide2, mctc2);
  mctc1 = coords_sub(mcposaPSDafter_curve, mcposabunker);
  mcposrbunker = rot_apply(mcrotabunker, mctc1);
  mcDEBUG_COMPONENT("bunker", mcposabunker, mcrotabunker)
  mccomp_posa[9] = mcposabunker;
  mccomp_posr[9] = mcposrbunker;
  mcNCounter[9]  = mcPCounter[9] = mcP2Counter[9] = 0;
  mcAbsorbProp[9]= 0;
    /* Component guide3. */
  /* Setting parameters for component guide3. */
  SIG_MESSAGE("guide3 (Init:SetPar)");
#line 58 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  if(0) strncpy(mccguide3_reflect, 0 ? 0 : "", 16384); else mccguide3_reflect[0]='\0';
#line 148 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccguide3_w1 = 0.02;
#line 148 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccguide3_h1 = .12;
#line 148 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccguide3_w2 = 0.02;
#line 148 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccguide3_h2 = .12;
#line 149 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccguide3_l = 12.275;
#line 149 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccguide3_R0 = R0;
#line 149 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccguide3_Qc = Qc;
#line 149 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccguide3_alpha = alpha;
#line 150 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccguide3_m = 1.6;
#line 150 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccguide3_W = W;
#line 12579 "./PSI_DMC.c"

  SIG_MESSAGE("guide3 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 151 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    (0)*DEG2RAD,
#line 151 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    (0)*DEG2RAD,
#line 151 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    (0)*DEG2RAD);
#line 12589 "./PSI_DMC.c"
  rot_mul(mctr1, mcrotabunker, mcrotaguide3);
  rot_transpose(mcrotabunker, mctr1);
  rot_mul(mcrotaguide3, mctr1, mcrotrguide3);
  mctc1 = coords_set(
#line 151 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    0,
#line 151 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    0,
#line 151 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    3.56);
#line 12600 "./PSI_DMC.c"
  rot_transpose(mcrotabunker, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaguide3 = coords_add(mcposabunker, mctc2);
  mctc1 = coords_sub(mcposabunker, mcposaguide3);
  mcposrguide3 = rot_apply(mcrotaguide3, mctc1);
  mcDEBUG_COMPONENT("guide3", mcposaguide3, mcrotaguide3)
  mccomp_posa[10] = mcposaguide3;
  mccomp_posr[10] = mcposrguide3;
  mcNCounter[10]  = mcPCounter[10] = mcP2Counter[10] = 0;
  mcAbsorbProp[10]= 0;
    /* Component guide4. */
  /* Setting parameters for component guide4. */
  SIG_MESSAGE("guide4 (Init:SetPar)");
#line 58 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  if(0) strncpy(mccguide4_reflect, 0 ? 0 : "", 16384); else mccguide4_reflect[0]='\0';
#line 154 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccguide4_w1 = 0.02;
#line 154 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccguide4_h1 = .12;
#line 154 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccguide4_w2 = 0.02;
#line 154 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccguide4_h2 = .12;
#line 155 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccguide4_l = 5.66;
#line 155 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccguide4_R0 = R0;
#line 155 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccguide4_Qc = Qc;
#line 155 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccguide4_alpha = alpha;
#line 156 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccguide4_m = 1.6;
#line 156 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccguide4_W = W;
#line 12636 "./PSI_DMC.c"

  SIG_MESSAGE("guide4 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 157 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    (0)*DEG2RAD,
#line 157 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    (0)*DEG2RAD,
#line 157 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    (0)*DEG2RAD);
#line 12646 "./PSI_DMC.c"
  rot_mul(mctr1, mcrotaguide3, mcrotaguide4);
  rot_transpose(mcrotaguide3, mctr1);
  rot_mul(mcrotaguide4, mctr1, mcrotrguide4);
  mctc1 = coords_set(
#line 157 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    0,
#line 157 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    0,
#line 157 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    15.8555);
#line 12657 "./PSI_DMC.c"
  rot_transpose(mcrotabunker, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaguide4 = coords_add(mcposabunker, mctc2);
  mctc1 = coords_sub(mcposaguide3, mcposaguide4);
  mcposrguide4 = rot_apply(mcrotaguide4, mctc1);
  mcDEBUG_COMPONENT("guide4", mcposaguide4, mcrotaguide4)
  mccomp_posa[11] = mcposaguide4;
  mccomp_posr[11] = mcposrguide4;
  mcNCounter[11]  = mcPCounter[11] = mcP2Counter[11] = 0;
  mcAbsorbProp[11]= 0;
    /* Component window1. */
  /* Setting parameters for component window1. */
  SIG_MESSAGE("window1 (Init:SetPar)");
#line 160 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccwindow1_thickness = 0.002;
#line 12673 "./PSI_DMC.c"

  SIG_MESSAGE("window1 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 12680 "./PSI_DMC.c"
  rot_mul(mctr1, mcrotaguide4, mcrotawindow1);
  rot_transpose(mcrotaguide4, mctr1);
  rot_mul(mcrotawindow1, mctr1, mcrotrwindow1);
  mctc1 = coords_set(
#line 161 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    0,
#line 161 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    0,
#line 161 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    5.66 + 1e-9);
#line 12691 "./PSI_DMC.c"
  rot_transpose(mcrotaguide4, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposawindow1 = coords_add(mcposaguide4, mctc2);
  mctc1 = coords_sub(mcposaguide4, mcposawindow1);
  mcposrwindow1 = rot_apply(mcrotawindow1, mctc1);
  mcDEBUG_COMPONENT("window1", mcposawindow1, mcrotawindow1)
  mccomp_posa[12] = mcposawindow1;
  mccomp_posr[12] = mcposrwindow1;
  mcNCounter[12]  = mcPCounter[12] = mcP2Counter[12] = 0;
  mcAbsorbProp[12]= 0;
    /* Component ydist_fluxpos. */
  /* Setting parameters for component ydist_fluxpos. */
  SIG_MESSAGE("ydist_fluxpos (Init:SetPar)");
#line 166 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  if("ydist_fluxpos.dat") strncpy(mccydist_fluxpos_filename, "ydist_fluxpos.dat" ? "ydist_fluxpos.dat" : "", 16384); else mccydist_fluxpos_filename[0]='\0';
#line 46 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccydist_fluxpos_xmin = -0.05;
#line 46 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccydist_fluxpos_xmax = 0.05;
#line 46 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccydist_fluxpos_ymin = -0.05;
#line 46 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccydist_fluxpos_ymax = 0.05;
#line 166 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccydist_fluxpos_xwidth = 0.120;
#line 167 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccydist_fluxpos_yheight = 0.02;
#line 47 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccydist_fluxpos_restore_neutron = 0;
#line 47 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccydist_fluxpos_nowritefile = 0;
#line 12723 "./PSI_DMC.c"

  SIG_MESSAGE("ydist_fluxpos (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 168 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    (0)*DEG2RAD,
#line 168 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    (0)*DEG2RAD,
#line 168 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    (90)*DEG2RAD);
#line 12733 "./PSI_DMC.c"
  rot_mul(mctr1, mcrotawindow1, mcrotaydist_fluxpos);
  rot_transpose(mcrotawindow1, mctr1);
  rot_mul(mcrotaydist_fluxpos, mctr1, mcrotrydist_fluxpos);
  mctc1 = coords_set(
#line 168 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    0,
#line 168 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    0,
#line 168 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    5.66 + 1e-8 + 0.01);
#line 12744 "./PSI_DMC.c"
  rot_transpose(mcrotaguide4, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaydist_fluxpos = coords_add(mcposaguide4, mctc2);
  mctc1 = coords_sub(mcposawindow1, mcposaydist_fluxpos);
  mcposrydist_fluxpos = rot_apply(mcrotaydist_fluxpos, mctc1);
  mcDEBUG_COMPONENT("ydist_fluxpos", mcposaydist_fluxpos, mcrotaydist_fluxpos)
  mccomp_posa[13] = mcposaydist_fluxpos;
  mccomp_posr[13] = mcposrydist_fluxpos;
  mcNCounter[13]  = mcPCounter[13] = mcP2Counter[13] = 0;
  mcAbsorbProp[13]= 0;
    /* Component PSD_fluxpos. */
  /* Setting parameters for component PSD_fluxpos. */
  SIG_MESSAGE("PSD_fluxpos (Init:SetPar)");
#line 171 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccPSD_fluxpos_nx = 100;
#line 171 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccPSD_fluxpos_ny = 100;
#line 171 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  if("xdist_fluxposy.dat") strncpy(mccPSD_fluxpos_filename, "xdist_fluxposy.dat" ? "xdist_fluxposy.dat" : "", 16384); else mccPSD_fluxpos_filename[0]='\0';
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccPSD_fluxpos_xmin = -0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccPSD_fluxpos_xmax = 0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccPSD_fluxpos_ymin = -0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccPSD_fluxpos_ymax = 0.05;
#line 172 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccPSD_fluxpos_xwidth = 0.02;
#line 172 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccPSD_fluxpos_yheight = 0.12;
#line 51 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccPSD_fluxpos_restore_neutron = 0;
#line 12778 "./PSI_DMC.c"

  SIG_MESSAGE("PSD_fluxpos (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 12785 "./PSI_DMC.c"
  rot_mul(mctr1, mcrotaguide4, mcrotaPSD_fluxpos);
  rot_transpose(mcrotaydist_fluxpos, mctr1);
  rot_mul(mcrotaPSD_fluxpos, mctr1, mcrotrPSD_fluxpos);
  mctc1 = coords_set(
#line 173 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    0,
#line 173 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    0,
#line 173 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    5.66 + 1e-7 + 0.01);
#line 12796 "./PSI_DMC.c"
  rot_transpose(mcrotaguide4, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaPSD_fluxpos = coords_add(mcposaguide4, mctc2);
  mctc1 = coords_sub(mcposaydist_fluxpos, mcposaPSD_fluxpos);
  mcposrPSD_fluxpos = rot_apply(mcrotaPSD_fluxpos, mctc1);
  mcDEBUG_COMPONENT("PSD_fluxpos", mcposaPSD_fluxpos, mcrotaPSD_fluxpos)
  mccomp_posa[14] = mcposaPSD_fluxpos;
  mccomp_posr[14] = mcposrPSD_fluxpos;
  mcNCounter[14]  = mcPCounter[14] = mcP2Counter[14] = 0;
  mcAbsorbProp[14]= 0;
    /* Component xdist_flux_pos. */
  /* Setting parameters for component xdist_flux_pos. */
  SIG_MESSAGE("xdist_flux_pos (Init:SetPar)");
#line 176 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  if("xdist_fluxpos.dat") strncpy(mccxdist_flux_pos_filename, "xdist_fluxpos.dat" ? "xdist_fluxpos.dat" : "", 16384); else mccxdist_flux_pos_filename[0]='\0';
#line 46 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccxdist_flux_pos_xmin = -0.05;
#line 46 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccxdist_flux_pos_xmax = 0.05;
#line 46 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccxdist_flux_pos_ymin = -0.05;
#line 46 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccxdist_flux_pos_ymax = 0.05;
#line 176 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccxdist_flux_pos_xwidth = 0.020;
#line 177 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccxdist_flux_pos_yheight = 0.12;
#line 47 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccxdist_flux_pos_restore_neutron = 0;
#line 47 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccxdist_flux_pos_nowritefile = 0;
#line 12828 "./PSI_DMC.c"

  SIG_MESSAGE("xdist_flux_pos (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 12835 "./PSI_DMC.c"
  rot_mul(mctr1, mcrotaPSD_fluxpos, mcrotaxdist_flux_pos);
  rot_transpose(mcrotaPSD_fluxpos, mctr1);
  rot_mul(mcrotaxdist_flux_pos, mctr1, mcrotrxdist_flux_pos);
  mctc1 = coords_set(
#line 178 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    0,
#line 178 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    0,
#line 178 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    1e-9);
#line 12846 "./PSI_DMC.c"
  rot_transpose(mcrotaPSD_fluxpos, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaxdist_flux_pos = coords_add(mcposaPSD_fluxpos, mctc2);
  mctc1 = coords_sub(mcposaPSD_fluxpos, mcposaxdist_flux_pos);
  mcposrxdist_flux_pos = rot_apply(mcrotaxdist_flux_pos, mctc1);
  mcDEBUG_COMPONENT("xdist_flux_pos", mcposaxdist_flux_pos, mcrotaxdist_flux_pos)
  mccomp_posa[15] = mcposaxdist_flux_pos;
  mccomp_posr[15] = mcposrxdist_flux_pos;
  mcNCounter[15]  = mcPCounter[15] = mcP2Counter[15] = 0;
  mcAbsorbProp[15]= 0;
    /* Component PSD_fluxposB. */
  /* Setting parameters for component PSD_fluxposB. */
  SIG_MESSAGE("PSD_fluxposB (Init:SetPar)");
#line 181 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccPSD_fluxposB_nx = 100;
#line 181 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccPSD_fluxposB_ny = 100;
#line 181 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  if("PSD_fluxposB.dat") strncpy(mccPSD_fluxposB_filename, "PSD_fluxposB.dat" ? "PSD_fluxposB.dat" : "", 16384); else mccPSD_fluxposB_filename[0]='\0';
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccPSD_fluxposB_xmin = -0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccPSD_fluxposB_xmax = 0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccPSD_fluxposB_ymin = -0.05;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccPSD_fluxposB_ymax = 0.05;
#line 182 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccPSD_fluxposB_xwidth = 0.02;
#line 182 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccPSD_fluxposB_yheight = 0.12;
#line 51 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccPSD_fluxposB_restore_neutron = 0;
#line 12880 "./PSI_DMC.c"

  SIG_MESSAGE("PSD_fluxposB (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 12887 "./PSI_DMC.c"
  rot_mul(mctr1, mcrotaguide4, mcrotaPSD_fluxposB);
  rot_transpose(mcrotaxdist_flux_pos, mctr1);
  rot_mul(mcrotaPSD_fluxposB, mctr1, mcrotrPSD_fluxposB);
  mctc1 = coords_set(
#line 183 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    0,
#line 183 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    0,
#line 183 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    6.24 -1e-7 -0.01);
#line 12898 "./PSI_DMC.c"
  rot_transpose(mcrotaguide4, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaPSD_fluxposB = coords_add(mcposaguide4, mctc2);
  mctc1 = coords_sub(mcposaxdist_flux_pos, mcposaPSD_fluxposB);
  mcposrPSD_fluxposB = rot_apply(mcrotaPSD_fluxposB, mctc1);
  mcDEBUG_COMPONENT("PSD_fluxposB", mcposaPSD_fluxposB, mcrotaPSD_fluxposB)
  mccomp_posa[16] = mcposaPSD_fluxposB;
  mccomp_posr[16] = mcposrPSD_fluxposB;
  mcNCounter[16]  = mcPCounter[16] = mcP2Counter[16] = 0;
  mcAbsorbProp[16]= 0;
    /* Component window2. */
  /* Setting parameters for component window2. */
  SIG_MESSAGE("window2 (Init:SetPar)");
#line 186 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccwindow2_thickness = 0.002;
#line 12914 "./PSI_DMC.c"

  SIG_MESSAGE("window2 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 12921 "./PSI_DMC.c"
  rot_mul(mctr1, mcrotaPSD_fluxposB, mcrotawindow2);
  rot_transpose(mcrotaPSD_fluxposB, mctr1);
  rot_mul(mcrotawindow2, mctr1, mcrotrwindow2);
  mctc1 = coords_set(
#line 187 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    0,
#line 187 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    0,
#line 187 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    1e-9);
#line 12932 "./PSI_DMC.c"
  rot_transpose(mcrotaPSD_fluxposB, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposawindow2 = coords_add(mcposaPSD_fluxposB, mctc2);
  mctc1 = coords_sub(mcposaPSD_fluxposB, mcposawindow2);
  mcposrwindow2 = rot_apply(mcrotawindow2, mctc1);
  mcDEBUG_COMPONENT("window2", mcposawindow2, mcrotawindow2)
  mccomp_posa[17] = mcposawindow2;
  mccomp_posr[17] = mcposrwindow2;
  mcNCounter[17]  = mcPCounter[17] = mcP2Counter[17] = 0;
  mcAbsorbProp[17]= 0;
    /* Component in_slit. */
  /* Setting parameters for component in_slit. */
  SIG_MESSAGE("in_slit (Init:SetPar)");
#line 192 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccin_slit_xmin = -0.01;
#line 192 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccin_slit_xmax = 0.01;
#line 192 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccin_slit_ymin = -0.06;
#line 192 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccin_slit_ymax = 0.06;
#line 46 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccin_slit_radius = 0;
#line 46 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccin_slit_xwidth = 0;
#line 46 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccin_slit_yheight = 0;
#line 12960 "./PSI_DMC.c"

  SIG_MESSAGE("in_slit (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 12967 "./PSI_DMC.c"
  rot_mul(mctr1, mcrotawindow2, mcrotain_slit);
  rot_transpose(mcrotawindow2, mctr1);
  rot_mul(mcrotain_slit, mctr1, mcrotrin_slit);
  mctc1 = coords_set(
#line 193 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    0,
#line 193 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    0,
#line 193 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    0.0021);
#line 12978 "./PSI_DMC.c"
  rot_transpose(mcrotawindow2, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposain_slit = coords_add(mcposawindow2, mctc2);
  mctc1 = coords_sub(mcposawindow2, mcposain_slit);
  mcposrin_slit = rot_apply(mcrotain_slit, mctc1);
  mcDEBUG_COMPONENT("in_slit", mcposain_slit, mcrotain_slit)
  mccomp_posa[18] = mcposain_slit;
  mccomp_posr[18] = mcposrin_slit;
  mcNCounter[18]  = mcPCounter[18] = mcP2Counter[18] = 0;
  mcAbsorbProp[18]= 0;
    /* Component lambda_in. */
  /* Setting parameters for component lambda_in. */
  SIG_MESSAGE("lambda_in (Init:SetPar)");
#line 196 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  if("L_in.dat") strncpy(mcclambda_in_filename, "L_in.dat" ? "L_in.dat" : "", 16384); else mcclambda_in_filename[0]='\0';
#line 196 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mcclambda_in_xmin = -0.011;
#line 196 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mcclambda_in_xmax = 0.011;
#line 196 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mcclambda_in_ymin = -0.061;
#line 196 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mcclambda_in_ymax = 0.061;
#line 51 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mcclambda_in_xwidth = 0;
#line 51 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mcclambda_in_yheight = 0;
#line 196 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mcclambda_in_Lmin = 0;
#line 196 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mcclambda_in_Lmax = 2 * mciplambda;
#line 51 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mcclambda_in_restore_neutron = 0;
#line 51 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mcclambda_in_nowritefile = 0;
#line 13014 "./PSI_DMC.c"

  SIG_MESSAGE("lambda_in (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 13021 "./PSI_DMC.c"
  rot_mul(mctr1, mcrotain_slit, mcrotalambda_in);
  rot_transpose(mcrotain_slit, mctr1);
  rot_mul(mcrotalambda_in, mctr1, mcrotrlambda_in);
  mctc1 = coords_set(
#line 197 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    0,
#line 197 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    0,
#line 197 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    0.001);
#line 13032 "./PSI_DMC.c"
  rot_transpose(mcrotain_slit, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposalambda_in = coords_add(mcposain_slit, mctc2);
  mctc1 = coords_sub(mcposain_slit, mcposalambda_in);
  mcposrlambda_in = rot_apply(mcrotalambda_in, mctc1);
  mcDEBUG_COMPONENT("lambda_in", mcposalambda_in, mcrotalambda_in)
  mccomp_posa[19] = mcposalambda_in;
  mccomp_posr[19] = mcposrlambda_in;
  mcNCounter[19]  = mcPCounter[19] = mcP2Counter[19] = 0;
  mcAbsorbProp[19]= 0;
    /* Component sma. */
  /* Setting parameters for component sma. */
  SIG_MESSAGE("sma (Init:SetPar)");

  SIG_MESSAGE("sma (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 202 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    (0)*DEG2RAD,
#line 202 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    (OMA)*DEG2RAD,
#line 202 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    (0)*DEG2RAD);
#line 13055 "./PSI_DMC.c"
  rot_mul(mctr1, mcrotain_slit, mcrotasma);
  rot_transpose(mcrotalambda_in, mctr1);
  rot_mul(mcrotasma, mctr1, mcrotrsma);
  mctc1 = coords_set(
#line 202 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    0,
#line 202 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    0,
#line 202 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    0.65);
#line 13066 "./PSI_DMC.c"
  rot_transpose(mcrotain_slit, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposasma = coords_add(mcposain_slit, mctc2);
  mctc1 = coords_sub(mcposalambda_in, mcposasma);
  mcposrsma = rot_apply(mcrotasma, mctc1);
  mcDEBUG_COMPONENT("sma", mcposasma, mcrotasma)
  mccomp_posa[20] = mcposasma;
  mccomp_posr[20] = mcposrsma;
  mcNCounter[20]  = mcPCounter[20] = mcP2Counter[20] = 0;
  mcAbsorbProp[20]= 0;
    /* Component foc_mono. */
  /* Setting parameters for component foc_mono. */
  SIG_MESSAGE("foc_mono (Init:SetPar)");
#line 78 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  if(0) strncpy(mccfoc_mono_reflect, 0 ? 0 : "", 16384); else mccfoc_mono_reflect[0]='\0';
#line 205 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccfoc_mono_zwidth = 0.05;
#line 205 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccfoc_mono_yheight = 0.025;
#line 205 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccfoc_mono_gap = 0.0005;
#line 205 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccfoc_mono_NH = 1;
#line 205 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccfoc_mono_NV = 5;
#line 206 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccfoc_mono_mosaich = 38;
#line 206 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccfoc_mono_mosaicv = 38;
#line 206 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccfoc_mono_r0 = 0.7;
#line 206 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccfoc_mono_Q = mono_q;
#line 206 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccfoc_mono_RV = RV;
#line 207 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccfoc_mono_RH = 0;
#line 78 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccfoc_mono_DM = 0;
#line 78 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccfoc_mono_mosaic = 0;
#line 78 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccfoc_mono_width = 0;
#line 78 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccfoc_mono_height = 0;
#line 78 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccfoc_mono_verbose = 0;
#line 13114 "./PSI_DMC.c"

  SIG_MESSAGE("foc_mono (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 13121 "./PSI_DMC.c"
  rot_mul(mctr1, mcrotasma, mcrotafoc_mono);
  rot_transpose(mcrotasma, mctr1);
  rot_mul(mcrotafoc_mono, mctr1, mcrotrfoc_mono);
  mctc1 = coords_set(
#line 208 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    0,
#line 208 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    0,
#line 208 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    0);
#line 13132 "./PSI_DMC.c"
  rot_transpose(mcrotasma, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposafoc_mono = coords_add(mcposasma, mctc2);
  mctc1 = coords_sub(mcposasma, mcposafoc_mono);
  mcposrfoc_mono = rot_apply(mcrotafoc_mono, mctc1);
  mcDEBUG_COMPONENT("foc_mono", mcposafoc_mono, mcrotafoc_mono)
  mccomp_posa[21] = mcposafoc_mono;
  mccomp_posr[21] = mcposrfoc_mono;
  mcNCounter[21]  = mcPCounter[21] = mcP2Counter[21] = 0;
  mcAbsorbProp[21]= 0;
    /* Component msa. */
  /* Setting parameters for component msa. */
  SIG_MESSAGE("msa (Init:SetPar)");

  SIG_MESSAGE("msa (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 211 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    (0)*DEG2RAD,
#line 211 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    (TTM)*DEG2RAD,
#line 211 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    (0)*DEG2RAD);
#line 13155 "./PSI_DMC.c"
  rot_mul(mctr1, mcrotain_slit, mcrotamsa);
  rot_transpose(mcrotafoc_mono, mctr1);
  rot_mul(mcrotamsa, mctr1, mcrotrmsa);
  mctc1 = coords_set(
#line 211 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    0,
#line 211 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    0,
#line 211 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    0);
#line 13166 "./PSI_DMC.c"
  rot_transpose(mcrotasma, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposamsa = coords_add(mcposasma, mctc2);
  mctc1 = coords_sub(mcposafoc_mono, mcposamsa);
  mcposrmsa = rot_apply(mcrotamsa, mctc1);
  mcDEBUG_COMPONENT("msa", mcposamsa, mcrotamsa)
  mccomp_posa[22] = mcposamsa;
  mccomp_posr[22] = mcposrmsa;
  mcNCounter[22]  = mcPCounter[22] = mcP2Counter[22] = 0;
  mcAbsorbProp[22]= 0;
    /* Component out1_slit. */
  /* Setting parameters for component out1_slit. */
  SIG_MESSAGE("out1_slit (Init:SetPar)");
#line 214 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccout1_slit_xmin = -0.01;
#line 214 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccout1_slit_xmax = 0.01;
#line 214 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccout1_slit_ymin = -0.06;
#line 214 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccout1_slit_ymax = 0.06;
#line 46 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccout1_slit_radius = 0;
#line 46 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccout1_slit_xwidth = 0;
#line 46 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccout1_slit_yheight = 0;
#line 13194 "./PSI_DMC.c"

  SIG_MESSAGE("out1_slit (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 215 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    (0)*DEG2RAD,
#line 215 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    (0)*DEG2RAD,
#line 215 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    (0)*DEG2RAD);
#line 13204 "./PSI_DMC.c"
  rot_mul(mctr1, mcrotamsa, mcrotaout1_slit);
  rot_transpose(mcrotamsa, mctr1);
  rot_mul(mcrotaout1_slit, mctr1, mcrotrout1_slit);
  mctc1 = coords_set(
#line 215 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    0,
#line 215 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    0,
#line 215 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    0.2);
#line 13215 "./PSI_DMC.c"
  rot_transpose(mcrotamsa, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaout1_slit = coords_add(mcposamsa, mctc2);
  mctc1 = coords_sub(mcposamsa, mcposaout1_slit);
  mcposrout1_slit = rot_apply(mcrotaout1_slit, mctc1);
  mcDEBUG_COMPONENT("out1_slit", mcposaout1_slit, mcrotaout1_slit)
  mccomp_posa[23] = mcposaout1_slit;
  mccomp_posr[23] = mcposrout1_slit;
  mcNCounter[23]  = mcPCounter[23] = mcP2Counter[23] = 0;
  mcAbsorbProp[23]= 0;
    /* Component Amoin_slit. */
  /* Setting parameters for component Amoin_slit. */
  SIG_MESSAGE("Amoin_slit (Init:SetPar)");
#line 218 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccAmoin_slit_xmin = -0.01;
#line 218 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccAmoin_slit_xmax = 0.01;
#line 218 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccAmoin_slit_ymin = -0.06;
#line 218 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccAmoin_slit_ymax = 0.06;
#line 46 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccAmoin_slit_radius = 0;
#line 46 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccAmoin_slit_xwidth = 0;
#line 46 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccAmoin_slit_yheight = 0;
#line 13243 "./PSI_DMC.c"

  SIG_MESSAGE("Amoin_slit (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 219 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    (0)*DEG2RAD,
#line 219 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    (0)*DEG2RAD,
#line 219 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    (0)*DEG2RAD);
#line 13253 "./PSI_DMC.c"
  rot_mul(mctr1, mcrotamsa, mcrotaAmoin_slit);
  rot_transpose(mcrotaout1_slit, mctr1);
  rot_mul(mcrotaAmoin_slit, mctr1, mcrotrAmoin_slit);
  mctc1 = coords_set(
#line 219 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    0,
#line 219 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    0,
#line 219 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    0.325);
#line 13264 "./PSI_DMC.c"
  rot_transpose(mcrotamsa, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaAmoin_slit = coords_add(mcposamsa, mctc2);
  mctc1 = coords_sub(mcposaout1_slit, mcposaAmoin_slit);
  mcposrAmoin_slit = rot_apply(mcrotaAmoin_slit, mctc1);
  mcDEBUG_COMPONENT("Amoin_slit", mcposaAmoin_slit, mcrotaAmoin_slit)
  mccomp_posa[24] = mcposaAmoin_slit;
  mccomp_posr[24] = mcposrAmoin_slit;
  mcNCounter[24]  = mcPCounter[24] = mcP2Counter[24] = 0;
  mcAbsorbProp[24]= 0;
    /* Component Bmoin_slit. */
  /* Setting parameters for component Bmoin_slit. */
  SIG_MESSAGE("Bmoin_slit (Init:SetPar)");
#line 222 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccBmoin_slit_xmin = -0.01;
#line 222 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccBmoin_slit_xmax = 0.01;
#line 222 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccBmoin_slit_ymin = -0.06;
#line 222 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccBmoin_slit_ymax = 0.06;
#line 46 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccBmoin_slit_radius = 0;
#line 46 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccBmoin_slit_xwidth = 0;
#line 46 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccBmoin_slit_yheight = 0;
#line 13292 "./PSI_DMC.c"

  SIG_MESSAGE("Bmoin_slit (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 223 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    (0)*DEG2RAD,
#line 223 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    (0)*DEG2RAD,
#line 223 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    (0)*DEG2RAD);
#line 13302 "./PSI_DMC.c"
  rot_mul(mctr1, mcrotamsa, mcrotaBmoin_slit);
  rot_transpose(mcrotaAmoin_slit, mctr1);
  rot_mul(mcrotaBmoin_slit, mctr1, mcrotrBmoin_slit);
  mctc1 = coords_set(
#line 223 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    0,
#line 223 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    0,
#line 223 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    0.525);
#line 13313 "./PSI_DMC.c"
  rot_transpose(mcrotamsa, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaBmoin_slit = coords_add(mcposamsa, mctc2);
  mctc1 = coords_sub(mcposaAmoin_slit, mcposaBmoin_slit);
  mcposrBmoin_slit = rot_apply(mcrotaBmoin_slit, mctc1);
  mcDEBUG_COMPONENT("Bmoin_slit", mcposaBmoin_slit, mcrotaBmoin_slit)
  mccomp_posa[25] = mcposaBmoin_slit;
  mccomp_posr[25] = mcposrBmoin_slit;
  mcNCounter[25]  = mcPCounter[25] = mcP2Counter[25] = 0;
  mcAbsorbProp[25]= 0;
    /* Component out2_slit. */
  /* Setting parameters for component out2_slit. */
  SIG_MESSAGE("out2_slit (Init:SetPar)");
#line 226 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccout2_slit_xmin = -0.01;
#line 226 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccout2_slit_xmax = 0.01;
#line 226 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccout2_slit_ymin = -0.06;
#line 226 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccout2_slit_ymax = 0.06;
#line 46 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccout2_slit_radius = 0;
#line 46 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccout2_slit_xwidth = 0;
#line 46 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccout2_slit_yheight = 0;
#line 13341 "./PSI_DMC.c"

  SIG_MESSAGE("out2_slit (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 227 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    (0)*DEG2RAD,
#line 227 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    (0)*DEG2RAD,
#line 227 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    (0)*DEG2RAD);
#line 13351 "./PSI_DMC.c"
  rot_mul(mctr1, mcrotamsa, mcrotaout2_slit);
  rot_transpose(mcrotaBmoin_slit, mctr1);
  rot_mul(mcrotaout2_slit, mctr1, mcrotrout2_slit);
  mctc1 = coords_set(
#line 227 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    0,
#line 227 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    0,
#line 227 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    0.65);
#line 13362 "./PSI_DMC.c"
  rot_transpose(mcrotamsa, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaout2_slit = coords_add(mcposamsa, mctc2);
  mctc1 = coords_sub(mcposaBmoin_slit, mcposaout2_slit);
  mcposrout2_slit = rot_apply(mcrotaout2_slit, mctc1);
  mcDEBUG_COMPONENT("out2_slit", mcposaout2_slit, mcrotaout2_slit)
  mccomp_posa[26] = mcposaout2_slit;
  mccomp_posr[26] = mcposrout2_slit;
  mcNCounter[26]  = mcPCounter[26] = mcP2Counter[26] = 0;
  mcAbsorbProp[26]= 0;
    /* Component PSD_sample. */
  /* Setting parameters for component PSD_sample. */
  SIG_MESSAGE("PSD_sample (Init:SetPar)");
#line 231 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccPSD_sample_nx = 80;
#line 231 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccPSD_sample_ny = 80;
#line 231 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  if("PSD_sample.dat") strncpy(mccPSD_sample_filename, "PSD_sample.dat" ? "PSD_sample.dat" : "", 16384); else mccPSD_sample_filename[0]='\0';
#line 230 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccPSD_sample_xmin = -0.05;
#line 230 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccPSD_sample_xmax = 0.05;
#line 230 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccPSD_sample_ymin = -0.07;
#line 230 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccPSD_sample_ymax = 0.07;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccPSD_sample_xwidth = 0;
#line 50 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccPSD_sample_yheight = 0;
#line 51 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccPSD_sample_restore_neutron = 0;
#line 13396 "./PSI_DMC.c"

  SIG_MESSAGE("PSD_sample (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 13403 "./PSI_DMC.c"
  rot_mul(mctr1, mcrotamsa, mcrotaPSD_sample);
  rot_transpose(mcrotaout2_slit, mctr1);
  rot_mul(mcrotaPSD_sample, mctr1, mcrotrPSD_sample);
  mctc1 = coords_set(
#line 232 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    0,
#line 232 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    0,
#line 232 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    2.77);
#line 13414 "./PSI_DMC.c"
  rot_transpose(mcrotamsa, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaPSD_sample = coords_add(mcposamsa, mctc2);
  mctc1 = coords_sub(mcposaout2_slit, mcposaPSD_sample);
  mcposrPSD_sample = rot_apply(mcrotaPSD_sample, mctc1);
  mcDEBUG_COMPONENT("PSD_sample", mcposaPSD_sample, mcrotaPSD_sample)
  mccomp_posa[27] = mcposaPSD_sample;
  mccomp_posr[27] = mcposrPSD_sample;
  mcNCounter[27]  = mcPCounter[27] = mcP2Counter[27] = 0;
  mcAbsorbProp[27]= 0;
    /* Component lambda_sample. */
  /* Setting parameters for component lambda_sample. */
  SIG_MESSAGE("lambda_sample (Init:SetPar)");
#line 235 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  if("L_sample.dat") strncpy(mcclambda_sample_filename, "L_sample.dat" ? "L_sample.dat" : "", 16384); else mcclambda_sample_filename[0]='\0';
#line 235 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mcclambda_sample_xmin = - sample_radius;
#line 235 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mcclambda_sample_xmax = sample_radius;
#line 235 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mcclambda_sample_ymin = - sample_height / 2;
#line 235 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mcclambda_sample_ymax = sample_height / 2;
#line 51 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mcclambda_sample_xwidth = 0;
#line 51 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mcclambda_sample_yheight = 0;
#line 235 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mcclambda_sample_Lmin = mciplambda -0.2;
#line 235 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mcclambda_sample_Lmax = mciplambda + 0.2;
#line 51 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mcclambda_sample_restore_neutron = 0;
#line 51 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mcclambda_sample_nowritefile = 0;
#line 13450 "./PSI_DMC.c"

  SIG_MESSAGE("lambda_sample (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 13457 "./PSI_DMC.c"
  rot_mul(mctr1, mcrotamsa, mcrotalambda_sample);
  rot_transpose(mcrotaPSD_sample, mctr1);
  rot_mul(mcrotalambda_sample, mctr1, mcrotrlambda_sample);
  mctc1 = coords_set(
#line 236 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    0,
#line 236 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    0,
#line 236 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    2.81);
#line 13468 "./PSI_DMC.c"
  rot_transpose(mcrotamsa, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposalambda_sample = coords_add(mcposamsa, mctc2);
  mctc1 = coords_sub(mcposaPSD_sample, mcposalambda_sample);
  mcposrlambda_sample = rot_apply(mcrotalambda_sample, mctc1);
  mcDEBUG_COMPONENT("lambda_sample", mcposalambda_sample, mcrotalambda_sample)
  mccomp_posa[28] = mcposalambda_sample;
  mccomp_posr[28] = mcposrlambda_sample;
  mcNCounter[28]  = mcPCounter[28] = mcP2Counter[28] = 0;
  mcAbsorbProp[28]= 0;
    /* Component sa_arm. */
  /* Setting parameters for component sa_arm. */
  SIG_MESSAGE("sa_arm (Init:SetPar)");

  SIG_MESSAGE("sa_arm (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 240 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    (0)*DEG2RAD,
#line 240 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    (0)*DEG2RAD,
#line 240 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    (0)*DEG2RAD);
#line 13491 "./PSI_DMC.c"
  rot_mul(mctr1, mcrotamsa, mcrotasa_arm);
  rot_transpose(mcrotalambda_sample, mctr1);
  rot_mul(mcrotasa_arm, mctr1, mcrotrsa_arm);
  mctc1 = coords_set(
#line 239 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    0,
#line 239 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    0,
#line 239 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    2.82);
#line 13502 "./PSI_DMC.c"
  rot_transpose(mcrotamsa, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposasa_arm = coords_add(mcposamsa, mctc2);
  mctc1 = coords_sub(mcposalambda_sample, mcposasa_arm);
  mcposrsa_arm = rot_apply(mcrotasa_arm, mctc1);
  mcDEBUG_COMPONENT("sa_arm", mcposasa_arm, mcrotasa_arm)
  mccomp_posa[29] = mcposasa_arm;
  mccomp_posr[29] = mcposrsa_arm;
  mcNCounter[29]  = mcPCounter[29] = mcP2Counter[29] = 0;
  mcAbsorbProp[29]= 0;
    /* Component sample. */
  /* Setting parameters for component sample. */
  SIG_MESSAGE("sample (Init:SetPar)");
#line 244 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  if(mcipfilename) strncpy(mccsample_reflections, mcipfilename ? mcipfilename : "", 16384); else mccsample_reflections[0]='\0';
#line 207 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  if("NULL") strncpy(mccsample_geometry, "NULL" ? "NULL" : "", 16384); else mccsample_geometry[0]='\0';
#line 243 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccsample_radius = sample_radius;
#line 243 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccsample_yheight = sample_height;
#line 208 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccsample_xwidth = 0;
#line 208 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccsample_zdepth = 0;
#line 208 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccsample_thickness = 0;
#line 244 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccsample_pack = mcipPACK;
#line 209 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccsample_Vc = 0;
#line 209 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccsample_sigma_abs = 0;
#line 209 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccsample_sigma_inc = 0;
#line 209 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccsample_delta_d_d = 0;
#line 244 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccsample_p_inc = 0;
#line 244 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccsample_p_transmit = 0;
#line 243 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccsample_DW = mcipDw;
#line 210 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccsample_nb_atoms = 1;
#line 210 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccsample_d_omega = 0;
#line 243 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccsample_d_phi = mcipD_PHI;
#line 210 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccsample_tth_sign = 0;
#line 210 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccsample_p_interact = 0;
#line 211 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccsample_concentric = 0;
#line 211 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccsample_density = 0;
#line 211 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccsample_weight = 0;
#line 244 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccsample_barns = mcipBARNS;
#line 211 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccsample_Strain = 0;
#line 211 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccsample_focus_flip = 0;
#line 211 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccsample_target_index = 0;
#line 13570 "./PSI_DMC.c"

  SIG_MESSAGE("sample (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 13577 "./PSI_DMC.c"
  rot_mul(mctr1, mcrotasa_arm, mcrotasample);
  rot_transpose(mcrotasa_arm, mctr1);
  rot_mul(mcrotasample, mctr1, mcrotrsample);
  mctc1 = coords_set(
#line 245 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    0,
#line 245 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    0,
#line 245 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    0);
#line 13588 "./PSI_DMC.c"
  rot_transpose(mcrotasa_arm, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposasample = coords_add(mcposasa_arm, mctc2);
  mctc1 = coords_sub(mcposasa_arm, mcposasample);
  mcposrsample = rot_apply(mcrotasample, mctc1);
  mcDEBUG_COMPONENT("sample", mcposasample, mcrotasample)
  mccomp_posa[30] = mcposasample;
  mccomp_posr[30] = mcposrsample;
  mcNCounter[30]  = mcPCounter[30] = mcP2Counter[30] = 0;
  mcAbsorbProp[30]= 0;
    /* Component STOP. */
  /* Setting parameters for component STOP. */
  SIG_MESSAGE("STOP (Init:SetPar)");
#line 44 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccSTOP_xmin = -0.05;
#line 44 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccSTOP_xmax = 0.05;
#line 44 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccSTOP_ymin = -0.05;
#line 44 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccSTOP_ymax = 0.05;
#line 45 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccSTOP_xwidth = 0;
#line 45 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccSTOP_yheight = 0;
#line 247 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccSTOP_radius = 0.3;
#line 13616 "./PSI_DMC.c"

  SIG_MESSAGE("STOP (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 249 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    (0)*DEG2RAD,
#line 249 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    (0)*DEG2RAD,
#line 249 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    (0)*DEG2RAD);
#line 13626 "./PSI_DMC.c"
  rot_mul(mctr1, mcrotasa_arm, mcrotaSTOP);
  rot_transpose(mcrotasample, mctr1);
  rot_mul(mcrotaSTOP, mctr1, mcrotrSTOP);
  mctc1 = coords_set(
#line 248 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    0,
#line 248 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    0,
#line 248 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    1.4);
#line 13637 "./PSI_DMC.c"
  rot_transpose(mcrotasa_arm, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaSTOP = coords_add(mcposasa_arm, mctc2);
  mctc1 = coords_sub(mcposasample, mcposaSTOP);
  mcposrSTOP = rot_apply(mcrotaSTOP, mctc1);
  mcDEBUG_COMPONENT("STOP", mcposaSTOP, mcrotaSTOP)
  mccomp_posa[31] = mcposaSTOP;
  mccomp_posr[31] = mcposrSTOP;
  mcNCounter[31]  = mcPCounter[31] = mcP2Counter[31] = 0;
  mcAbsorbProp[31]= 0;
    /* Component Detector. */
  /* Setting parameters for component Detector. */
  SIG_MESSAGE("Detector (Init:SetPar)");
#line 253 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccDetector_xwidth = 3.0;
#line 253 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccDetector_yheight = 0.09;
#line 203 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccDetector_zdepth = 0;
#line 204 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccDetector_xmin = 0;
#line 204 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccDetector_xmax = 0;
#line 204 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccDetector_ymin = 0;
#line 204 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccDetector_ymax = 0;
#line 204 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccDetector_zmin = 0;
#line 204 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccDetector_zmax = 0;
#line 253 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccDetector_bins = 400;
#line 253 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccDetector_min = 19.9 + mcipSHIFT;
#line 253 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccDetector_max = 99.9 + mcipSHIFT;
#line 205 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccDetector_restore_neutron = 0;
#line 205 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccDetector_radius = 0;
#line 254 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  if("banana, theta") strncpy(mccDetector_options, "banana, theta" ? "banana, theta" : "", 16384); else mccDetector_options[0]='\0';
#line 253 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  if("detector.dat") strncpy(mccDetector_filename, "detector.dat" ? "detector.dat" : "", 16384); else mccDetector_filename[0]='\0';
#line 206 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  if("NULL") strncpy(mccDetector_geometry, "NULL" ? "NULL" : "", 16384); else mccDetector_geometry[0]='\0';
#line 207 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  if("NULL") strncpy(mccDetector_username1, "NULL" ? "NULL" : "", 16384); else mccDetector_username1[0]='\0';
#line 207 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  if("NULL") strncpy(mccDetector_username2, "NULL" ? "NULL" : "", 16384); else mccDetector_username2[0]='\0';
#line 207 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  if("NULL") strncpy(mccDetector_username3, "NULL" ? "NULL" : "", 16384); else mccDetector_username3[0]='\0';
#line 208 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
  mccDetector_nowritefile = 0;
#line 13693 "./PSI_DMC.c"

  SIG_MESSAGE("Detector (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 256 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    (0)*DEG2RAD,
#line 256 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    (0)*DEG2RAD,
#line 256 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    (180)*DEG2RAD);
#line 13703 "./PSI_DMC.c"
  rot_mul(mctr1, mcrotasa_arm, mcrotaDetector);
  rot_transpose(mcrotaSTOP, mctr1);
  rot_mul(mcrotaDetector, mctr1, mcrotrDetector);
  mctc1 = coords_set(
#line 255 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    0,
#line 255 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    0,
#line 255 "/zhome/89/0/38697/once/McStas-2.5_CPU_MPICC/PSI_DMC/PSI_DMC.instr"
    0);
#line 13714 "./PSI_DMC.c"
  rot_transpose(mcrotasa_arm, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaDetector = coords_add(mcposasa_arm, mctc2);
  mctc1 = coords_sub(mcposaSTOP, mcposaDetector);
  mcposrDetector = rot_apply(mcrotaDetector, mctc1);
  mcDEBUG_COMPONENT("Detector", mcposaDetector, mcrotaDetector)
  mccomp_posa[32] = mcposaDetector;
  mccomp_posr[32] = mcposrDetector;
  mcNCounter[32]  = mcPCounter[32] = mcP2Counter[32] = 0;
  mcAbsorbProp[32]= 0;
  /* Component initializations. */
  /* Initializations for component source_arm. */
  SIG_MESSAGE("source_arm (Init)");
#define mccompcurname  source_arm
#define mccompcurtype  Progress_bar
#define mccompcurindex 1
#define IntermediateCnts mccsource_arm_IntermediateCnts
#define StartTime mccsource_arm_StartTime
#define EndTime mccsource_arm_EndTime
#define CurrentTime mccsource_arm_CurrentTime
#define profile mccsource_arm_profile
#define percent mccsource_arm_percent
#define flag_save mccsource_arm_flag_save
#define minutes mccsource_arm_minutes
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
#line 13751 "./PSI_DMC.c"
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
#define mccompcurtype  Source_Maxwell_3
#define mccompcurindex 2
#define M mccsource_M
#define l_range mccsource_l_range
#define w_mult mccsource_w_mult
#define w_source mccsource_w_source
#define h_source mccsource_h_source
#define size mccsource_size
#define yheight mccsource_yheight
#define xwidth mccsource_xwidth
#define Lmin mccsource_Lmin
#define Lmax mccsource_Lmax
#define dist mccsource_dist
#define focus_xw mccsource_focus_xw
#define focus_yh mccsource_focus_yh
#define T1 mccsource_T1
#define T2 mccsource_T2
#define T3 mccsource_T3
#define I1 mccsource_I1
#define I2 mccsource_I2
#define I3 mccsource_I3
#define target_index mccsource_target_index
#define lambda0 mccsource_lambda0
#define dlambda mccsource_dlambda
#line 85 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../sources/Source_Maxwell_3.comp"
{
  if (target_index && !dist)
  {
    Coords ToTarget;
    double tx,ty,tz;
    ToTarget = coords_sub(POS_A_COMP_INDEX(INDEX_CURRENT_COMP+target_index),POS_A_CURRENT_COMP);
    ToTarget = rot_apply(ROT_A_CURRENT_COMP, ToTarget);
    coords_get(ToTarget, &tx, &ty, &tz);
    dist=sqrt(tx*tx+ty*ty+tz*tz);
  }

  if (size>0) {
    w_source = h_source = size;
  } else {
    w_source = xwidth;
    h_source = yheight;
  }
  if (lambda0) {
    Lmin=lambda0-dlambda;
    Lmax=lambda0+dlambda;
  }
  l_range = Lmax-Lmin;
  w_mult = w_source*h_source*1.0e4;     /* source area correction */
  w_mult *= l_range;            /* wavelength range correction */
  w_mult *= 1.0/mcget_ncount();   /* correct for # neutron rays */

  if (w_source <0 || h_source < 0 || Lmin <= 0 || Lmax <= 0 || dist <= 0 || T1 <= 0 || T2 <= 0|| T3 <= 0 || Lmax<=Lmin) {
      printf("Source_Maxwell_3: %s: Error in input parameter values!\n"
             "ERROR          Exiting\n",
           NAME_CURRENT_COMP);
      exit(0);
  }

}
#line 13826 "./PSI_DMC.c"
#undef dlambda
#undef lambda0
#undef target_index
#undef I3
#undef I2
#undef I1
#undef T3
#undef T2
#undef T1
#undef focus_yh
#undef focus_xw
#undef dist
#undef Lmax
#undef Lmin
#undef xwidth
#undef yheight
#undef size
#undef h_source
#undef w_source
#undef w_mult
#undef l_range
#undef M
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component PSDbefore_guides. */
  SIG_MESSAGE("PSDbefore_guides (Init)");
#define mccompcurname  PSDbefore_guides
#define mccompcurtype  PSD_monitor
#define mccompcurindex 3
#define PSD_N mccPSDbefore_guides_PSD_N
#define PSD_p mccPSDbefore_guides_PSD_p
#define PSD_p2 mccPSDbefore_guides_PSD_p2
#define nx mccPSDbefore_guides_nx
#define ny mccPSDbefore_guides_ny
#define filename mccPSDbefore_guides_filename
#define xmin mccPSDbefore_guides_xmin
#define xmax mccPSDbefore_guides_xmax
#define ymin mccPSDbefore_guides_ymin
#define ymax mccPSDbefore_guides_ymax
#define xwidth mccPSDbefore_guides_xwidth
#define yheight mccPSDbefore_guides_yheight
#define restore_neutron mccPSDbefore_guides_restore_neutron
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
#line 13896 "./PSI_DMC.c"
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

  /* Initializations for component l_mon_source. */
  SIG_MESSAGE("l_mon_source (Init)");
#define mccompcurname  l_mon_source
#define mccompcurtype  L_monitor
#define mccompcurindex 4
#define nL mccl_mon_source_nL
#define L_N mccl_mon_source_L_N
#define L_p mccl_mon_source_L_p
#define L_p2 mccl_mon_source_L_p2
#define filename mccl_mon_source_filename
#define xmin mccl_mon_source_xmin
#define xmax mccl_mon_source_xmax
#define ymin mccl_mon_source_ymin
#define ymax mccl_mon_source_ymax
#define xwidth mccl_mon_source_xwidth
#define yheight mccl_mon_source_yheight
#define Lmin mccl_mon_source_Lmin
#define Lmax mccl_mon_source_Lmax
#define restore_neutron mccl_mon_source_restore_neutron
#define nowritefile mccl_mon_source_nowritefile
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
#line 13955 "./PSI_DMC.c"
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

  /* Initializations for component guide1. */
  SIG_MESSAGE("guide1 (Init)");
#define mccompcurname  guide1
#define mccompcurtype  Guide
#define mccompcurindex 5
#define pTable mccguide1_pTable
#define reflect mccguide1_reflect
#define w1 mccguide1_w1
#define h1 mccguide1_h1
#define w2 mccguide1_w2
#define h2 mccguide1_h2
#define l mccguide1_l
#define R0 mccguide1_R0
#define Qc mccguide1_Qc
#define alpha mccguide1_alpha
#define m mccguide1_m
#define W mccguide1_W
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
#line 14011 "./PSI_DMC.c"
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

  /* Initializations for component PSDbefore_curve. */
  SIG_MESSAGE("PSDbefore_curve (Init)");
#define mccompcurname  PSDbefore_curve
#define mccompcurtype  PSD_monitor
#define mccompcurindex 6
#define PSD_N mccPSDbefore_curve_PSD_N
#define PSD_p mccPSDbefore_curve_PSD_p
#define PSD_p2 mccPSDbefore_curve_PSD_p2
#define nx mccPSDbefore_curve_nx
#define ny mccPSDbefore_curve_ny
#define filename mccPSDbefore_curve_filename
#define xmin mccPSDbefore_curve_xmin
#define xmax mccPSDbefore_curve_xmax
#define ymin mccPSDbefore_curve_ymin
#define ymax mccPSDbefore_curve_ymax
#define xwidth mccPSDbefore_curve_xwidth
#define yheight mccPSDbefore_curve_yheight
#define restore_neutron mccPSDbefore_curve_restore_neutron
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
#line 14071 "./PSI_DMC.c"
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

  /* Initializations for component guide2. */
  SIG_MESSAGE("guide2 (Init)");
#define mccompcurname  guide2
#define mccompcurtype  Bender
#define mccompcurindex 7
#define bk mccguide2_bk
#define mWin mccguide2_mWin
#define w mccguide2_w
#define h mccguide2_h
#define r mccguide2_r
#define Win mccguide2_Win
#define k mccguide2_k
#define d mccguide2_d
#define l mccguide2_l
#define R0a mccguide2_R0a
#define Qca mccguide2_Qca
#define alphaa mccguide2_alphaa
#define ma mccguide2_ma
#define Wa mccguide2_Wa
#define R0i mccguide2_R0i
#define Qci mccguide2_Qci
#define alphai mccguide2_alphai
#define mi mccguide2_mi
#define Wi mccguide2_Wi
#define R0s mccguide2_R0s
#define Qcs mccguide2_Qcs
#define alphas mccguide2_alphas
#define ms mccguide2_ms
#define Ws mccguide2_Ws
#line 112 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Bender.comp"
{
if (r <0)
      { fprintf(stderr,"Bender: error: %s: to bend in the other direction\n", NAME_CURRENT_COMP);
        fprintf(stderr,"        rotate comp on z-axis by 180 deg.\n"); exit(-1); }

      if (k*d > w)
      { fprintf(stderr,"Bender: error: %s has (k*d > w).\n", NAME_CURRENT_COMP);
        exit(-1); }
      if (w*h*r*Win*k == 0)
      { fprintf(stderr,"Bender: error: %s has one of w,h,r,Win,k null.\n", NAME_CURRENT_COMP);
        exit(-1); }
      /* width of one channel + thickness d of partition */
      mWin = Win;
      if (l!= 0 && r != 0) mWin = (double)l/(double)r;
      bk=(w+d)/k;
      if (mcgravitation) fprintf(stderr,"WARNING: Bender: %s: "
        "This component produces wrong results with gravitation !\n",
        NAME_CURRENT_COMP);
}
#line 14138 "./PSI_DMC.c"
#undef Ws
#undef ms
#undef alphas
#undef Qcs
#undef R0s
#undef Wi
#undef mi
#undef alphai
#undef Qci
#undef R0i
#undef Wa
#undef ma
#undef alphaa
#undef Qca
#undef R0a
#undef l
#undef d
#undef k
#undef Win
#undef r
#undef h
#undef w
#undef mWin
#undef bk
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component PSDafter_curve. */
  SIG_MESSAGE("PSDafter_curve (Init)");
#define mccompcurname  PSDafter_curve
#define mccompcurtype  PSD_monitor
#define mccompcurindex 8
#define PSD_N mccPSDafter_curve_PSD_N
#define PSD_p mccPSDafter_curve_PSD_p
#define PSD_p2 mccPSDafter_curve_PSD_p2
#define nx mccPSDafter_curve_nx
#define ny mccPSDafter_curve_ny
#define filename mccPSDafter_curve_filename
#define xmin mccPSDafter_curve_xmin
#define xmax mccPSDafter_curve_xmax
#define ymin mccPSDafter_curve_ymin
#define ymax mccPSDafter_curve_ymax
#define xwidth mccPSDafter_curve_xwidth
#define yheight mccPSDafter_curve_yheight
#define restore_neutron mccPSDafter_curve_restore_neutron
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
#line 14210 "./PSI_DMC.c"
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

  /* Initializations for component bunker. */
  SIG_MESSAGE("bunker (Init)");
#define mccompcurname  bunker
#define mccompcurtype  Guide
#define mccompcurindex 9
#define pTable mccbunker_pTable
#define reflect mccbunker_reflect
#define w1 mccbunker_w1
#define h1 mccbunker_h1
#define w2 mccbunker_w2
#define h2 mccbunker_h2
#define l mccbunker_l
#define R0 mccbunker_R0
#define Qc mccbunker_Qc
#define alpha mccbunker_alpha
#define m mccbunker_m
#define W mccbunker_W
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
#line 14264 "./PSI_DMC.c"
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

  /* Initializations for component guide3. */
  SIG_MESSAGE("guide3 (Init)");
#define mccompcurname  guide3
#define mccompcurtype  Guide
#define mccompcurindex 10
#define pTable mccguide3_pTable
#define reflect mccguide3_reflect
#define w1 mccguide3_w1
#define h1 mccguide3_h1
#define w2 mccguide3_w2
#define h2 mccguide3_h2
#define l mccguide3_l
#define R0 mccguide3_R0
#define Qc mccguide3_Qc
#define alpha mccguide3_alpha
#define m mccguide3_m
#define W mccguide3_W
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
#line 14317 "./PSI_DMC.c"
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

  /* Initializations for component guide4. */
  SIG_MESSAGE("guide4 (Init)");
#define mccompcurname  guide4
#define mccompcurtype  Guide
#define mccompcurindex 11
#define pTable mccguide4_pTable
#define reflect mccguide4_reflect
#define w1 mccguide4_w1
#define h1 mccguide4_h1
#define w2 mccguide4_w2
#define h2 mccguide4_h2
#define l mccguide4_l
#define R0 mccguide4_R0
#define Qc mccguide4_Qc
#define alpha mccguide4_alpha
#define m mccguide4_m
#define W mccguide4_W
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
#line 14370 "./PSI_DMC.c"
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

  /* Initializations for component window1. */
  SIG_MESSAGE("window1 (Init)");

  /* Initializations for component ydist_fluxpos. */
  SIG_MESSAGE("ydist_fluxpos (Init)");
#define mccompcurname  ydist_fluxpos
#define mccompcurtype  PSDlin_monitor
#define mccompcurindex 13
#define nx mccydist_fluxpos_nx
#define PSDlin_N mccydist_fluxpos_PSDlin_N
#define PSDlin_p mccydist_fluxpos_PSDlin_p
#define PSDlin_p2 mccydist_fluxpos_PSDlin_p2
#define filename mccydist_fluxpos_filename
#define xmin mccydist_fluxpos_xmin
#define xmax mccydist_fluxpos_xmax
#define ymin mccydist_fluxpos_ymin
#define ymax mccydist_fluxpos_ymax
#define xwidth mccydist_fluxpos_xwidth
#define yheight mccydist_fluxpos_yheight
#define restore_neutron mccydist_fluxpos_restore_neutron
#define nowritefile mccydist_fluxpos_nowritefile
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
#line 14429 "./PSI_DMC.c"
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

  /* Initializations for component PSD_fluxpos. */
  SIG_MESSAGE("PSD_fluxpos (Init)");
#define mccompcurname  PSD_fluxpos
#define mccompcurtype  PSD_monitor
#define mccompcurindex 14
#define PSD_N mccPSD_fluxpos_PSD_N
#define PSD_p mccPSD_fluxpos_PSD_p
#define PSD_p2 mccPSD_fluxpos_PSD_p2
#define nx mccPSD_fluxpos_nx
#define ny mccPSD_fluxpos_ny
#define filename mccPSD_fluxpos_filename
#define xmin mccPSD_fluxpos_xmin
#define xmax mccPSD_fluxpos_xmax
#define ymin mccPSD_fluxpos_ymin
#define ymax mccPSD_fluxpos_ymax
#define xwidth mccPSD_fluxpos_xwidth
#define yheight mccPSD_fluxpos_yheight
#define restore_neutron mccPSD_fluxpos_restore_neutron
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
#line 14490 "./PSI_DMC.c"
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

  /* Initializations for component xdist_flux_pos. */
  SIG_MESSAGE("xdist_flux_pos (Init)");
#define mccompcurname  xdist_flux_pos
#define mccompcurtype  PSDlin_monitor
#define mccompcurindex 15
#define nx mccxdist_flux_pos_nx
#define PSDlin_N mccxdist_flux_pos_PSDlin_N
#define PSDlin_p mccxdist_flux_pos_PSDlin_p
#define PSDlin_p2 mccxdist_flux_pos_PSDlin_p2
#define filename mccxdist_flux_pos_filename
#define xmin mccxdist_flux_pos_xmin
#define xmax mccxdist_flux_pos_xmax
#define ymin mccxdist_flux_pos_ymin
#define ymax mccxdist_flux_pos_ymax
#define xwidth mccxdist_flux_pos_xwidth
#define yheight mccxdist_flux_pos_yheight
#define restore_neutron mccxdist_flux_pos_restore_neutron
#define nowritefile mccxdist_flux_pos_nowritefile
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
#line 14547 "./PSI_DMC.c"
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

  /* Initializations for component PSD_fluxposB. */
  SIG_MESSAGE("PSD_fluxposB (Init)");
#define mccompcurname  PSD_fluxposB
#define mccompcurtype  PSD_monitor
#define mccompcurindex 16
#define PSD_N mccPSD_fluxposB_PSD_N
#define PSD_p mccPSD_fluxposB_PSD_p
#define PSD_p2 mccPSD_fluxposB_PSD_p2
#define nx mccPSD_fluxposB_nx
#define ny mccPSD_fluxposB_ny
#define filename mccPSD_fluxposB_filename
#define xmin mccPSD_fluxposB_xmin
#define xmax mccPSD_fluxposB_xmax
#define ymin mccPSD_fluxposB_ymin
#define ymax mccPSD_fluxposB_ymax
#define xwidth mccPSD_fluxposB_xwidth
#define yheight mccPSD_fluxposB_yheight
#define restore_neutron mccPSD_fluxposB_restore_neutron
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
#line 14608 "./PSI_DMC.c"
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

  /* Initializations for component window2. */
  SIG_MESSAGE("window2 (Init)");

  /* Initializations for component in_slit. */
  SIG_MESSAGE("in_slit (Init)");
#define mccompcurname  in_slit
#define mccompcurtype  Slit
#define mccompcurindex 18
#define xmin mccin_slit_xmin
#define xmax mccin_slit_xmax
#define ymin mccin_slit_ymin
#define ymax mccin_slit_ymax
#define radius mccin_slit_radius
#define xwidth mccin_slit_xwidth
#define yheight mccin_slit_yheight
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
#line 14661 "./PSI_DMC.c"
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

  /* Initializations for component lambda_in. */
  SIG_MESSAGE("lambda_in (Init)");
#define mccompcurname  lambda_in
#define mccompcurtype  L_monitor
#define mccompcurindex 19
#define nL mcclambda_in_nL
#define L_N mcclambda_in_L_N
#define L_p mcclambda_in_L_p
#define L_p2 mcclambda_in_L_p2
#define filename mcclambda_in_filename
#define xmin mcclambda_in_xmin
#define xmax mcclambda_in_xmax
#define ymin mcclambda_in_ymin
#define ymax mcclambda_in_ymax
#define xwidth mcclambda_in_xwidth
#define yheight mcclambda_in_yheight
#define Lmin mcclambda_in_Lmin
#define Lmax mcclambda_in_Lmax
#define restore_neutron mcclambda_in_restore_neutron
#define nowritefile mcclambda_in_nowritefile
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
#line 14714 "./PSI_DMC.c"
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

  /* Initializations for component sma. */
  SIG_MESSAGE("sma (Init)");

  /* Initializations for component foc_mono. */
  SIG_MESSAGE("foc_mono (Init)");
#define mccompcurname  foc_mono
#define mccompcurtype  Monochromator_2foc
#define mccompcurindex 21
#define mos_y mccfoc_mono_mos_y
#define mos_z mccfoc_mono_mos_z
#define mono_Q mccfoc_mono_mono_Q
#define SlabWidth mccfoc_mono_SlabWidth
#define SlabHeight mccfoc_mono_SlabHeight
#define rTable mccfoc_mono_rTable
#define reflect mccfoc_mono_reflect
#define zwidth mccfoc_mono_zwidth
#define yheight mccfoc_mono_yheight
#define gap mccfoc_mono_gap
#define NH mccfoc_mono_NH
#define NV mccfoc_mono_NV
#define mosaich mccfoc_mono_mosaich
#define mosaicv mccfoc_mono_mosaicv
#define r0 mccfoc_mono_r0
#define Q mccfoc_mono_Q
#define RV mccfoc_mono_RV
#define RH mccfoc_mono_RH
#define DM mccfoc_mono_DM
#define mosaic mccfoc_mono_mosaic
#define width mccfoc_mono_width
#define height mccfoc_mono_height
#define verbose mccfoc_mono_verbose
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
#line 14813 "./PSI_DMC.c"
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

  /* Initializations for component msa. */
  SIG_MESSAGE("msa (Init)");

  /* Initializations for component out1_slit. */
  SIG_MESSAGE("out1_slit (Init)");
#define mccompcurname  out1_slit
#define mccompcurtype  Slit
#define mccompcurindex 23
#define xmin mccout1_slit_xmin
#define xmax mccout1_slit_xmax
#define ymin mccout1_slit_ymin
#define ymax mccout1_slit_ymax
#define radius mccout1_slit_radius
#define xwidth mccout1_slit_xwidth
#define yheight mccout1_slit_yheight
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
#line 14876 "./PSI_DMC.c"
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

  /* Initializations for component Amoin_slit. */
  SIG_MESSAGE("Amoin_slit (Init)");
#define mccompcurname  Amoin_slit
#define mccompcurtype  Slit
#define mccompcurindex 24
#define xmin mccAmoin_slit_xmin
#define xmax mccAmoin_slit_xmax
#define ymin mccAmoin_slit_ymin
#define ymax mccAmoin_slit_ymax
#define radius mccAmoin_slit_radius
#define xwidth mccAmoin_slit_xwidth
#define yheight mccAmoin_slit_yheight
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
#line 14920 "./PSI_DMC.c"
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

  /* Initializations for component Bmoin_slit. */
  SIG_MESSAGE("Bmoin_slit (Init)");
#define mccompcurname  Bmoin_slit
#define mccompcurtype  Slit
#define mccompcurindex 25
#define xmin mccBmoin_slit_xmin
#define xmax mccBmoin_slit_xmax
#define ymin mccBmoin_slit_ymin
#define ymax mccBmoin_slit_ymax
#define radius mccBmoin_slit_radius
#define xwidth mccBmoin_slit_xwidth
#define yheight mccBmoin_slit_yheight
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
#line 14964 "./PSI_DMC.c"
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

  /* Initializations for component out2_slit. */
  SIG_MESSAGE("out2_slit (Init)");
#define mccompcurname  out2_slit
#define mccompcurtype  Slit
#define mccompcurindex 26
#define xmin mccout2_slit_xmin
#define xmax mccout2_slit_xmax
#define ymin mccout2_slit_ymin
#define ymax mccout2_slit_ymax
#define radius mccout2_slit_radius
#define xwidth mccout2_slit_xwidth
#define yheight mccout2_slit_yheight
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
#line 15008 "./PSI_DMC.c"
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

  /* Initializations for component PSD_sample. */
  SIG_MESSAGE("PSD_sample (Init)");
#define mccompcurname  PSD_sample
#define mccompcurtype  PSD_monitor
#define mccompcurindex 27
#define PSD_N mccPSD_sample_PSD_N
#define PSD_p mccPSD_sample_PSD_p
#define PSD_p2 mccPSD_sample_PSD_p2
#define nx mccPSD_sample_nx
#define ny mccPSD_sample_ny
#define filename mccPSD_sample_filename
#define xmin mccPSD_sample_xmin
#define xmax mccPSD_sample_xmax
#define ymin mccPSD_sample_ymin
#define ymax mccPSD_sample_ymax
#define xwidth mccPSD_sample_xwidth
#define yheight mccPSD_sample_yheight
#define restore_neutron mccPSD_sample_restore_neutron
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
#line 15063 "./PSI_DMC.c"
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

  /* Initializations for component lambda_sample. */
  SIG_MESSAGE("lambda_sample (Init)");
#define mccompcurname  lambda_sample
#define mccompcurtype  L_monitor
#define mccompcurindex 28
#define nL mcclambda_sample_nL
#define L_N mcclambda_sample_L_N
#define L_p mcclambda_sample_L_p
#define L_p2 mcclambda_sample_L_p2
#define filename mcclambda_sample_filename
#define xmin mcclambda_sample_xmin
#define xmax mcclambda_sample_xmax
#define ymin mcclambda_sample_ymin
#define ymax mcclambda_sample_ymax
#define xwidth mcclambda_sample_xwidth
#define yheight mcclambda_sample_yheight
#define Lmin mcclambda_sample_Lmin
#define Lmax mcclambda_sample_Lmax
#define restore_neutron mcclambda_sample_restore_neutron
#define nowritefile mcclambda_sample_nowritefile
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
#line 15122 "./PSI_DMC.c"
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

  /* Initializations for component sa_arm. */
  SIG_MESSAGE("sa_arm (Init)");

  /* Initializations for component sample. */
  SIG_MESSAGE("sample (Init)");
#define mccompcurname  sample
#define mccompcurtype  PowderN
#define mccompcurindex 30
#define format mccsample_format
#define line_info mccsample_line_info
#define columns mccsample_columns
#define offdata mccsample_offdata
#define reflections mccsample_reflections
#define geometry mccsample_geometry
#define radius mccsample_radius
#define yheight mccsample_yheight
#define xwidth mccsample_xwidth
#define zdepth mccsample_zdepth
#define thickness mccsample_thickness
#define pack mccsample_pack
#define Vc mccsample_Vc
#define sigma_abs mccsample_sigma_abs
#define sigma_inc mccsample_sigma_inc
#define delta_d_d mccsample_delta_d_d
#define p_inc mccsample_p_inc
#define p_transmit mccsample_p_transmit
#define DW mccsample_DW
#define nb_atoms mccsample_nb_atoms
#define d_omega mccsample_d_omega
#define d_phi mccsample_d_phi
#define tth_sign mccsample_tth_sign
#define p_interact mccsample_p_interact
#define concentric mccsample_concentric
#define density mccsample_density
#define weight mccsample_weight
#define barns mccsample_barns
#define Strain mccsample_Strain
#define focus_flip mccsample_focus_flip
#define target_index mccsample_target_index
#line 560 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../samples/PowderN.comp"
{
  columns = (int[])format;

  int i=0;
  struct line_data *L;
  line_info.Dd       = delta_d_d;
  line_info.DWfactor = DW;
  line_info.V_0      = Vc;
  line_info.rho      = density;
  line_info.at_weight= weight;
  line_info.at_nb    = nb_atoms;
  line_info.sigma_a  = sigma_abs;
  line_info.sigma_i  = sigma_inc;
  line_info.flag_barns=barns;
  line_info.shape    = 0;
  line_info.flag_warning=0;
  line_info.Epsilon  = Strain;
  line_info.radius_i =line_info.xwidth_i=line_info.yheight_i=line_info.zdepth_i=0;
  line_info.v  = 0;
  line_info.Nq = 0;
  line_info.v_min = FLT_MAX; line_info.v_max = 0;
  line_info.neutron_passed=0;
  line_info.nb_reuses = line_info.nb_refl = line_info.nb_refl_count = 0;
  line_info.xs_compute= line_info.xs_reuse= line_info.xs_calls =0;
  for (i=0; i< 9; i++) line_info.column_order[i] = columns[i];
  strncpy(line_info.compname, NAME_CURRENT_COMP, 256);

  line_info.shape=-1; /* -1:no shape, 0:cyl, 1:box, 2:sphere, 3:any-shape  */
  if (geometry && strlen(geometry) && strcmp(geometry, "NULL") && strcmp(geometry, "0")) {
	  if (off_init(geometry, xwidth, yheight, zdepth, 0, &offdata)) {
      line_info.shape=3; thickness=0; concentric=0;
    }
  }
  else if (xwidth && yheight && zdepth)  line_info.shape=1; /* box */
  else if (radius > 0 && yheight)        line_info.shape=0; /* cylinder */
  else if (radius > 0 && !yheight)       line_info.shape=2; /* sphere */

  if (line_info.shape < 0)
    exit(fprintf(stderr,"PowderN: %s: sample has invalid dimensions.\n"
                        "ERROR    Please check parameter values (xwidth, yheight, zdepth, radius).\n", NAME_CURRENT_COMP));
  if (thickness) {
    if (radius && (radius < fabs(thickness))) {
      MPI_MASTER(
      printf("PowderN: %s: hollow sample thickness is larger than its volume (sphere/cylinder).\n"
                     "WARNING  Please check parameter values. Using bulk sample (thickness=0).\n", NAME_CURRENT_COMP);
      );
      thickness=0;
    }
    else if (!radius && (xwidth < 2*fabs(thickness) || yheight < 2*fabs(thickness) || zdepth < 2*fabs(thickness))) {
      MPI_MASTER(
      printf("PowderN: %s: hollow sample thickness is larger than its volume (box).\n"
                     "WARNING  Please check parameter values.\n", NAME_CURRENT_COMP);
      );
    }
  }

  if (concentric && thickness==0) {
    MPI_MASTER(
    printf("PowderN: %s:Can not use concentric mode\n"
           "WARNING     on non hollow shape. Ignoring.\n",
           NAME_CURRENT_COMP);
    );
    concentric=0;
  }

  if (thickness>0) {
    if (radius>thickness) {
      line_info.radius_i=radius-thickness;
    } else {
      if (xwidth>2*thickness)  line_info.xwidth_i =xwidth -2*thickness;
      if (yheight>2*thickness) line_info.yheight_i=yheight-2*thickness;
      if (zdepth>2*thickness)  line_info.zdepth_i =zdepth -2*thickness;
    }
  } else if (thickness<0) {
    thickness = fabs(thickness);
    if (radius) {
      line_info.radius_i=radius;
      radius=line_info.radius_i+thickness;
    } else {
      line_info.xwidth_i =xwidth;
      line_info.yheight_i=yheight;
      line_info.zdepth_i =zdepth;
      xwidth   =xwidth +2*thickness;
      yheight  =yheight+2*thickness;
      zdepth   =zdepth +2*thickness;
    }
  }

  if (!line_info.yheight_i) {
    line_info.yheight_i = yheight;
  }
  if (p_interact) {
    if (p_interact < p_inc) { double tmp=p_interact; p_interact=p_inc; p_inc=tmp; }
    p_transmit = 1-p_interact-p_inc;
  }

  if (p_inc + p_transmit > 1) {
    MPI_MASTER(
    printf("PowderN: %s: You have requested an unmeaningful choice of the 'p_inc' and 'p_transmit' parameters (sum is %g, exeeding 1). Fixing.\n",
                 NAME_CURRENT_COMP, p_inc+p_transmit);
    );
    if (p_inc > p_transmit) p_transmit=1-2*p_inc;
    else p_transmit=1-2*p_inc;
  } else if (p_inc + p_transmit == 1) {
    MPI_MASTER(
    printf("PowderN: %s: You have requested all neutrons be attenuated\n"
           "WARNING  or incoherently scattered!\n", NAME_CURRENT_COMP);
    );
  }

  if (concentric) {
    MPI_MASTER(
    printf("PowderN: %s: Concentric mode - remember to include the 'opposite' copy of this component !\n"
           "WARNING  The equivalent, 'opposite' comp should have concentric=0\n", NAME_CURRENT_COMP);
    );
    if (p_transmit == 0) {
      MPI_MASTER(
      printf("PowderN: %s: Concentric mode and p_transmit==0 !?\n"
             "WARNING  Don't you want any transmitted neutrons?\n", NAME_CURRENT_COMP);
      );
    }
  }

  if (reflections && strlen(reflections) && strcmp(reflections, "NULL") && strcmp(reflections, "0")) {
    i = read_line_data(reflections, &line_info);
    if (i == 0)
      exit(fprintf(stderr,"PowderN: %s: reflection file %s is not valid.\n"
                          "ERROR    Please check file format (laz or lau).\n", NAME_CURRENT_COMP, reflections));
  }

  /* compute the scattering unit density from material weight and density */
  /* the weight of the scattering element is the chemical formula molecular weight
   * times the nb of chemical formulae in the scattering element (nb_atoms) */
  if (!line_info.V_0 && line_info.at_nb > 0
    && line_info.at_weight > 0 && line_info.rho > 0) {
    /* molar volume [cm^3/mol] = weight [g/mol] / density [g/cm^3] */
    /* atom density per Angs^3 = [mol/cm^3] * N_Avogadro *(1e-8)^3 */
    line_info.V_0 = line_info.at_nb
      /(line_info.rho/line_info.at_weight/1e24*6.02214199e23);
  }

  /* the scattering unit cross sections are the chemical formula onces
   * times the nb of chemical formulae in the scattering element */
  if (line_info.at_nb > 0) {
    line_info.sigma_a *= line_info.at_nb; line_info.sigma_i *= line_info.at_nb;
  }

  if (line_info.sigma_a<0) line_info.sigma_a=0;
  if (line_info.sigma_i<0) line_info.sigma_i=0;

  if (line_info.V_0 <= 0)
  MPI_MASTER(
    printf("PowderN: %s: density/unit cell volume is NULL (Vc). Unactivating component.\n", NAME_CURRENT_COMP);
  );

  if (line_info.V_0 > 0 && p_inc && !line_info.sigma_i) {
  MPI_MASTER(
    printf("PowderN: %s: WARNING: You have requested statistics for incoherent scattering but not defined sigma_inc!\n", NAME_CURRENT_COMP);
  );
  }

  if (line_info.flag_barns) { /* Factor 100 to convert from barns to fm^2 */
    line_info.XsectionFactor = 100;
  } else {
    line_info.XsectionFactor = 1;
  }

  if (line_info.V_0 > 0 && i) {
    L = line_info.list;

    line_info.q_v = malloc(line_info.count*sizeof(double));
    line_info.w_v = malloc(line_info.count*sizeof(double));
    line_info.my_s_v2 = malloc(line_info.count*sizeof(double));
    if (!line_info.q_v || !line_info.w_v || !line_info.my_s_v2)
      exit(fprintf(stderr,"PowderN: %s: ERROR allocating memory (init)\n", NAME_CURRENT_COMP));
    for(i=0; i<line_info.count; i++)
    {
      line_info.my_s_v2[i] = 4*PI*PI*PI*pack*(L[i].DWfactor ? L[i].DWfactor : 1)
                 /(line_info.V_0*line_info.V_0*V2K*V2K)
                 *(L[i].j * L[i].F2 / L[i].q)*line_info.XsectionFactor;
      /* Is not yet divided by v^2 */
      /* Squires [3.103] */
      line_info.q_v[i] = L[i].q*K2V;
      line_info.w_v[i] = L[i].w;
    }
  }
  if (line_info.V_0 > 0) {
    /* Is not yet divided by v */
    line_info.my_a_v = pack*line_info.sigma_a/line_info.V_0*2200*100;   // Factor 100 to convert from barns to fm^2
    line_info.my_inc = pack*line_info.sigma_i/line_info.V_0*100;   // Factor 100 to convert from barns to fm^2
    MPI_MASTER(
    printf("PowderN: %s: Vc=%g [Angs] sigma_abs=%g [barn] sigma_inc=%g [barn] reflections=%s\n",
      NAME_CURRENT_COMP, line_info.V_0, line_info.sigma_a, line_info.sigma_i, reflections && strlen(reflections) ? reflections : "NULL");
    );
  }
  
  /* update JS, 1/7/2017
    Get target coordinates relative to the local reference frame.
  */
    if (target_index) {
		Coords ToTarget;
		ToTarget = coords_sub(POS_A_COMP_INDEX(INDEX_CURRENT_COMP+target_index),POS_A_CURRENT_COMP);
		ToTarget = rot_apply(ROT_A_CURRENT_COMP, ToTarget);
		coords_get(ToTarget, &tgt_x, &tgt_y, &tgt_z);
		NORM(tgt_x, tgt_y, tgt_z);
		printf("PowderN: Target direction = (%g %g %g)\n",tgt_x, tgt_y, tgt_z);
	} else {
		tgt_x=0.0;
		tgt_y=0.0;
		tgt_z=1.0;
	}	   

}
#line 15395 "./PSI_DMC.c"
#undef target_index
#undef focus_flip
#undef Strain
#undef barns
#undef weight
#undef density
#undef concentric
#undef p_interact
#undef tth_sign
#undef d_phi
#undef d_omega
#undef nb_atoms
#undef DW
#undef p_transmit
#undef p_inc
#undef delta_d_d
#undef sigma_inc
#undef sigma_abs
#undef Vc
#undef pack
#undef thickness
#undef zdepth
#undef xwidth
#undef yheight
#undef radius
#undef geometry
#undef reflections
#undef offdata
#undef columns
#undef line_info
#undef format
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component STOP. */
  SIG_MESSAGE("STOP (Init)");
#define mccompcurname  STOP
#define mccompcurtype  Beamstop
#define mccompcurindex 31
#define xmin mccSTOP_xmin
#define xmax mccSTOP_xmax
#define ymin mccSTOP_ymin
#define ymax mccSTOP_ymax
#define xwidth mccSTOP_xwidth
#define yheight mccSTOP_yheight
#define radius mccSTOP_radius
#line 50 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Beamstop.comp"
{
if (xwidth  > 0) { xmax = xwidth/2;  xmin = -xmax; }
  if (yheight > 0) { ymax = yheight/2; ymin = -ymax; }

  if (xmin == 0 && xmax == 0 && ymin == 0 & ymax == 0 && radius == 0)
  { fprintf(stderr,"Beamstop: %s: Error: give geometry\n", NAME_CURRENT_COMP); exit(-1); }
}
#line 15451 "./PSI_DMC.c"
#undef radius
#undef yheight
#undef xwidth
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component Detector. */
  SIG_MESSAGE("Detector (Init)");
#define mccompcurname  Detector
#define mccompcurtype  Monitor_nD
#define mccompcurindex 32
#define user1 mccDetector_user1
#define user2 mccDetector_user2
#define user3 mccDetector_user3
#define DEFS mccDetector_DEFS
#define Vars mccDetector_Vars
#define detector mccDetector_detector
#define offdata mccDetector_offdata
#define xwidth mccDetector_xwidth
#define yheight mccDetector_yheight
#define zdepth mccDetector_zdepth
#define xmin mccDetector_xmin
#define xmax mccDetector_xmax
#define ymin mccDetector_ymin
#define ymax mccDetector_ymax
#define zmin mccDetector_zmin
#define zmax mccDetector_zmax
#define bins mccDetector_bins
#define min mccDetector_min
#define max mccDetector_max
#define restore_neutron mccDetector_restore_neutron
#define radius mccDetector_radius
#define options mccDetector_options
#define filename mccDetector_filename
#define geometry mccDetector_geometry
#define username1 mccDetector_username1
#define username2 mccDetector_username2
#define username3 mccDetector_username3
#define nowritefile mccDetector_nowritefile
#line 232 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/Monitor_nD.comp"
{
  char tmp[CHAR_BUF_LENGTH];
  strcpy(Vars.compcurname, NAME_CURRENT_COMP);
  if (options != NULL)
    strncpy(Vars.option, options, CHAR_BUF_LENGTH);
  else {
    strcpy(Vars.option, "x y");
    printf("Monitor_nD: %s has no option specified. Setting to PSD ('x y') monitor.\n", NAME_CURRENT_COMP);
  }
  Vars.compcurpos = POS_A_CURRENT_COMP;

  if (strstr(Vars.option, "source"))
    strcat(Vars.option, " list, x y z vx vy vz t sx sy sz ");

  if (bins) { sprintf(tmp, " all bins=%ld ", (long)bins); strcat(Vars.option, tmp); }
  if (min > -FLT_MAX && max < FLT_MAX) { sprintf(tmp, " all limits=[%g %g]", min, max); strcat(Vars.option, tmp); }
  else if (min > -FLT_MAX) { sprintf(tmp, " all min=%g", min); strcat(Vars.option, tmp); }
  else if (max <  FLT_MAX) { sprintf(tmp, " all max=%g", max); strcat(Vars.option, tmp); }

  strncpy(Vars.UserName1,
    username1 && strlen(username1) && strcmp(username1, "0") && strcmp(username1, "NULL") ?
    username1 : "", 128);
  strncpy(Vars.UserName2,
    username2 && strlen(username2) && strcmp(username2, "0") && strcmp(username2, "NULL") ?
    username2 : "", 128);
  strncpy(Vars.UserName3,
    username3 && strlen(username3) && strcmp(username3, "0") && strcmp(username3, "NULL") ?
    username3 : "", 128);
  if (radius) {
    xwidth = zdepth = 2*radius;
    if (yheight && !strstr(Vars.option, "cylinder") && !strstr(Vars.option, "banana") && !strstr(Vars.option, "sphere"))
      strcat(Vars.option, " banana");
    else if (!yheight && !strstr(Vars.option ,"sphere")) {
      strcat(Vars.option, " sphere");
      yheight=2*radius;
    }
  }
  int offflag=0;
  if (geometry && strlen(geometry) && strcmp(geometry,"0") && strcmp(geometry, "NULL")) {
    if (!off_init(  geometry, xwidth, yheight, zdepth, 1, &offdata )) {
      printf("Monitor_nD: %s could not initiate the OFF geometry %s. \n"
             "            Defaulting to normal Monitor dimensions.\n",
             NAME_CURRENT_COMP, geometry);
      strcpy(geometry, "");
    } else {
      offflag=1;
    }
  }

  if (!radius && !xwidth && !yheight && !zdepth && !xmin && !xmax && !ymin && !ymax &&
    !strstr(Vars.option, "previous") && (!geometry || !strlen(geometry)))
    exit(printf("Monitor_nD: %s has no dimension specified. Aborting (radius, xwidth, yheight, zdepth, previous, geometry).\n", NAME_CURRENT_COMP));

  Monitor_nD_Init(&DEFS, &Vars, xwidth, yheight, zdepth, xmin,xmax,ymin,ymax,zmin,zmax,offflag);

  if (Vars.Flag_OFF) {
    offdata.mantidflag=Vars.Flag_mantid;
    offdata.mantidoffset=Vars.Coord_Min[Vars.Coord_Number-1];
  }


  if (filename && strlen(filename) && strcmp(filename,"NULL") && strcmp(filename,"0"))
    strncpy(Vars.Mon_File, filename, 128);

  /* check if user given filename with ext will be used more than once */
  if ( ((Vars.Flag_Multiple && Vars.Coord_Number > 1) || Vars.Flag_List) && strchr(Vars.Mon_File,'.') )
  { char *XY; XY = strrchr(Vars.Mon_File,'.'); *XY='_'; }

  if (restore_neutron) Vars.Flag_parallel=1;
  detector.m = 0;

#ifdef USE_MPI
MPI_MASTER(
  if (strstr(Vars.option, "auto") && mpi_node_count > 1)
    printf("Monitor_nD: %s is using automatic limits option 'auto' together with MPI.\n"
           "WARNING     this may create incorrect distributions (but integrated flux will be right).\n", NAME_CURRENT_COMP);
);
#endif
}
#line 15576 "./PSI_DMC.c"
#undef nowritefile
#undef username3
#undef username2
#undef username1
#undef geometry
#undef filename
#undef options
#undef radius
#undef restore_neutron
#undef max
#undef min
#undef bins
#undef zmax
#undef zmin
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef zdepth
#undef yheight
#undef xwidth
#undef offdata
#undef detector
#undef Vars
#undef DEFS
#undef user3
#undef user2
#undef user1
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
  /* SPLIT counter for component foc_mono */
  int mcSplit_foc_mono=0;
  /* SPLIT counter for component sample */
  int mcSplit_sample=0;
  /* TRACE Component source_arm [1] */
  mccoordschange(mcposrsource_arm, mcrotrsource_arm,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component source_arm (without coords transformations) */
  mcJumpTrace_source_arm:
  SIG_MESSAGE("source_arm (Trace)");
  mcDEBUG_COMP("source_arm")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompsource_arm
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
#define mccompcurname  source_arm
#define mccompcurtype  Progress_bar
#define mccompcurindex 1
#define IntermediateCnts mccsource_arm_IntermediateCnts
#define StartTime mccsource_arm_StartTime
#define EndTime mccsource_arm_EndTime
#define CurrentTime mccsource_arm_CurrentTime
{   /* Declarations of source_arm=Progress_bar() SETTING parameters. */
char* profile = mccsource_arm_profile;
MCNUM percent = mccsource_arm_percent;
MCNUM flag_save = mccsource_arm_flag_save;
MCNUM minutes = mccsource_arm_minutes;
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
#line 15766 "./PSI_DMC.c"
}   /* End of source_arm=Progress_bar() SETTING parameter declarations. */
#undef CurrentTime
#undef EndTime
#undef StartTime
#undef IntermediateCnts
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompsource_arm:
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
#define mccompcurtype  Source_Maxwell_3
#define mccompcurindex 2
#define M mccsource_M
#define l_range mccsource_l_range
#define w_mult mccsource_w_mult
#define w_source mccsource_w_source
#define h_source mccsource_h_source
{   /* Declarations of source=Source_Maxwell_3() SETTING parameters. */
MCNUM size = mccsource_size;
MCNUM yheight = mccsource_yheight;
MCNUM xwidth = mccsource_xwidth;
MCNUM Lmin = mccsource_Lmin;
MCNUM Lmax = mccsource_Lmax;
MCNUM dist = mccsource_dist;
MCNUM focus_xw = mccsource_focus_xw;
MCNUM focus_yh = mccsource_focus_yh;
MCNUM T1 = mccsource_T1;
MCNUM T2 = mccsource_T2;
MCNUM T3 = mccsource_T3;
MCNUM I1 = mccsource_I1;
MCNUM I2 = mccsource_I2;
MCNUM I3 = mccsource_I3;
int target_index = mccsource_target_index;
MCNUM lambda0 = mccsource_lambda0;
MCNUM dlambda = mccsource_dlambda;
#line 121 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../sources/Source_Maxwell_3.comp"
{
  double v,tau_l,E,lambda,k,r,xf,yf,dx,dy,w_focus;

  t=0;
  z=0;
  x = 0.5*w_source*randpm1();
  y = 0.5*h_source*randpm1();         /* Choose initial position */

  randvec_target_rect_real(&xf, &yf, &r, &w_focus,
		      0, 0, dist, focus_xw, focus_yh, ROT_A_CURRENT_COMP, x, y, z, 2);

  dx = xf-x;
  dy = yf-y;
  r = sqrt(dx*dx+dy*dy+dist*dist);

  lambda = Lmin+l_range*rand01();    /* Choose from uniform distribution */
  k = 2*PI/lambda;
  v = K2V*k;

  vz = v*dist/r;
  vy = v*dy/r;
  vx = v*dx/r;


/*  printf("pos0 (%g %g %g), pos1 (%g %g %g), r: %g, v (%g %g %g), v %g\n",
  x,y,z,xf,yf,dist,r,vx,vy,vz, v);
  printf("l %g, w_focus %g \n", lambda, w_focus);  */

  p *= w_mult*w_focus;                /* Correct for target focusing etc */
  p *= I1*SM3_Maxwell(lambda,T1)+I2*SM3_Maxwell(lambda,T2)+I3*SM3_Maxwell(lambda,T3);
                                        /* Calculate true intensity */
}
#line 15931 "./PSI_DMC.c"
}   /* End of source=Source_Maxwell_3() SETTING parameter declarations. */
#undef h_source
#undef w_source
#undef w_mult
#undef l_range
#undef M
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

  /* TRACE Component PSDbefore_guides [3] */
  mccoordschange(mcposrPSDbefore_guides, mcrotrPSDbefore_guides,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component PSDbefore_guides (without coords transformations) */
  mcJumpTrace_PSDbefore_guides:
  SIG_MESSAGE("PSDbefore_guides (Trace)");
  mcDEBUG_COMP("PSDbefore_guides")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompPSDbefore_guides
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
#define mccompcurname  PSDbefore_guides
#define mccompcurtype  PSD_monitor
#define mccompcurindex 3
#define PSD_N mccPSDbefore_guides_PSD_N
#define PSD_p mccPSDbefore_guides_PSD_p
#define PSD_p2 mccPSDbefore_guides_PSD_p2
{   /* Declarations of PSDbefore_guides=PSD_monitor() SETTING parameters. */
int nx = mccPSDbefore_guides_nx;
int ny = mccPSDbefore_guides_ny;
char* filename = mccPSDbefore_guides_filename;
MCNUM xmin = mccPSDbefore_guides_xmin;
MCNUM xmax = mccPSDbefore_guides_xmax;
MCNUM ymin = mccPSDbefore_guides_ymin;
MCNUM ymax = mccPSDbefore_guides_ymax;
MCNUM xwidth = mccPSDbefore_guides_xwidth;
MCNUM yheight = mccPSDbefore_guides_yheight;
MCNUM restore_neutron = mccPSDbefore_guides_restore_neutron;
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
#line 16070 "./PSI_DMC.c"
}   /* End of PSDbefore_guides=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompPSDbefore_guides:
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

  /* TRACE Component l_mon_source [4] */
  mccoordschange(mcposrl_mon_source, mcrotrl_mon_source,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component l_mon_source (without coords transformations) */
  mcJumpTrace_l_mon_source:
  SIG_MESSAGE("l_mon_source (Trace)");
  mcDEBUG_COMP("l_mon_source")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompl_mon_source
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
#define mccompcurname  l_mon_source
#define mccompcurtype  L_monitor
#define mccompcurindex 4
#define nL mccl_mon_source_nL
#define L_N mccl_mon_source_L_N
#define L_p mccl_mon_source_L_p
#define L_p2 mccl_mon_source_L_p2
{   /* Declarations of l_mon_source=L_monitor() SETTING parameters. */
char* filename = mccl_mon_source_filename;
MCNUM xmin = mccl_mon_source_xmin;
MCNUM xmax = mccl_mon_source_xmax;
MCNUM ymin = mccl_mon_source_ymin;
MCNUM ymax = mccl_mon_source_ymax;
MCNUM xwidth = mccl_mon_source_xwidth;
MCNUM yheight = mccl_mon_source_yheight;
MCNUM Lmin = mccl_mon_source_Lmin;
MCNUM Lmax = mccl_mon_source_Lmax;
MCNUM restore_neutron = mccl_mon_source_restore_neutron;
int nowritefile = mccl_mon_source_nowritefile;
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
#line 16216 "./PSI_DMC.c"
}   /* End of l_mon_source=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompl_mon_source:
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

  /* TRACE Component guide1 [5] */
  mccoordschange(mcposrguide1, mcrotrguide1,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component guide1 (without coords transformations) */
  mcJumpTrace_guide1:
  SIG_MESSAGE("guide1 (Trace)");
  mcDEBUG_COMP("guide1")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompguide1
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
#define mccompcurname  guide1
#define mccompcurtype  Guide
#define mccompcurindex 5
#define pTable mccguide1_pTable
{   /* Declarations of guide1=Guide() SETTING parameters. */
char* reflect = mccguide1_reflect;
MCNUM w1 = mccguide1_w1;
MCNUM h1 = mccguide1_h1;
MCNUM w2 = mccguide1_w2;
MCNUM h2 = mccguide1_h2;
MCNUM l = mccguide1_l;
MCNUM R0 = mccguide1_R0;
MCNUM Qc = mccguide1_Qc;
MCNUM alpha = mccguide1_alpha;
MCNUM m = mccguide1_m;
MCNUM W = mccguide1_W;
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
#line 16445 "./PSI_DMC.c"
}   /* End of guide1=Guide() SETTING parameter declarations. */
#undef pTable
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompguide1:
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

  /* TRACE Component PSDbefore_curve [6] */
  mccoordschange(mcposrPSDbefore_curve, mcrotrPSDbefore_curve,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component PSDbefore_curve (without coords transformations) */
  mcJumpTrace_PSDbefore_curve:
  SIG_MESSAGE("PSDbefore_curve (Trace)");
  mcDEBUG_COMP("PSDbefore_curve")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompPSDbefore_curve
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
#define mccompcurname  PSDbefore_curve
#define mccompcurtype  PSD_monitor
#define mccompcurindex 6
#define PSD_N mccPSDbefore_curve_PSD_N
#define PSD_p mccPSDbefore_curve_PSD_p
#define PSD_p2 mccPSDbefore_curve_PSD_p2
{   /* Declarations of PSDbefore_curve=PSD_monitor() SETTING parameters. */
int nx = mccPSDbefore_curve_nx;
int ny = mccPSDbefore_curve_ny;
char* filename = mccPSDbefore_curve_filename;
MCNUM xmin = mccPSDbefore_curve_xmin;
MCNUM xmax = mccPSDbefore_curve_xmax;
MCNUM ymin = mccPSDbefore_curve_ymin;
MCNUM ymax = mccPSDbefore_curve_ymax;
MCNUM xwidth = mccPSDbefore_curve_xwidth;
MCNUM yheight = mccPSDbefore_curve_yheight;
MCNUM restore_neutron = mccPSDbefore_curve_restore_neutron;
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
#line 16580 "./PSI_DMC.c"
}   /* End of PSDbefore_curve=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompPSDbefore_curve:
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

  /* TRACE Component guide2 [7] */
  mccoordschange(mcposrguide2, mcrotrguide2,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component guide2 (without coords transformations) */
  mcJumpTrace_guide2:
  SIG_MESSAGE("guide2 (Trace)");
  mcDEBUG_COMP("guide2")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompguide2
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
#define mccompcurname  guide2
#define mccompcurtype  Bender
#define mccompcurindex 7
#define bk mccguide2_bk
#define mWin mccguide2_mWin
{   /* Declarations of guide2=Bender() SETTING parameters. */
MCNUM w = mccguide2_w;
MCNUM h = mccguide2_h;
MCNUM r = mccguide2_r;
MCNUM Win = mccguide2_Win;
MCNUM k = mccguide2_k;
MCNUM d = mccguide2_d;
MCNUM l = mccguide2_l;
MCNUM R0a = mccguide2_R0a;
MCNUM Qca = mccguide2_Qca;
MCNUM alphaa = mccguide2_alphaa;
MCNUM ma = mccguide2_ma;
MCNUM Wa = mccguide2_Wa;
MCNUM R0i = mccguide2_R0i;
MCNUM Qci = mccguide2_Qci;
MCNUM alphai = mccguide2_alphai;
MCNUM mi = mccguide2_mi;
MCNUM Wi = mccguide2_Wi;
MCNUM R0s = mccguide2_R0s;
MCNUM Qcs = mccguide2_Qcs;
MCNUM alphas = mccguide2_alphas;
MCNUM ms = mccguide2_ms;
MCNUM Ws = mccguide2_Ws;
#line 133 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Bender.comp"
{
    int i,num,numa,numi;
    double dru,ab,dab,R,Q,Ta,vpl;
    double einmWin,ausmWin,zykmWin,aeumWin,innmWin,ref,innref,aeuref;
    double einzei,auszei,zykzei;

    /* does the neutron hit the bender at the entrance? */
    PROP_Z0;
    if ((fabs(x)<w/2) && (fabs(y)<h/2))
    {
      /*** reflections in the XZ-plane ***/

      /* distance between neutron and concave side of the channel at the entrance */
      dru=floor((w/2-x)/bk)*bk;
      ab=w/2.0-x-dru;

      /* radius of the channel */
      R=r-dru;

      /* does the neutron hit the partition at the entrance? */
      if (ab<bk-d)
      {
        double aeu[] = {R0a, Qca, alphaa, ma, Wa};
        /* velocity in the XZ-plane */
        vpl=sqrt(vx*vx+vz*vz);

        /* divergence of the neutron at the entrance */
        einmWin=atan(vx/vz);

        /* maximal distance between neutron and concave side of the channel */
        dab=R-cos(einmWin)*(R-ab);

        /* reflection angle at the concave side */
        aeumWin=acos((R-dab)/R);

        /* reflection coefficient at the concave side */
        Q=2.0*V2K*vpl*sin(aeumWin);
        StdReflecFunc(Q, aeu, &aeuref);

        /* does the neutron hit the convex side of the channel? */
        innmWin=0.0;
        innref=1.0;
        if (dab>bk-d)
        {
           double inn[] = {R0i, Qci, alphai, mi, Wi};
           /* reflection coefficient at the convex side */
           innmWin=acos((R-dab)/(R-bk+d));
           Q=2.0*V2K*vpl*sin(innmWin);
           StdReflecFunc(Q, inn, &innref);
        }

        /* divergence of the neutron at the exit */
        zykmWin=2.0*(aeumWin-innmWin);
        ausmWin=fmod(mWin+einmWin+aeumWin-innmWin
          *(1.0+SIGN(einmWin)),zykmWin)-zykmWin/2.0;
        ausmWin+=innmWin*SIGN(ausmWin);

        /* number of reflections at the concave side */
        numa=(mWin+einmWin+aeumWin-innmWin*(1.0+SIGN(einmWin)))/zykmWin;

        /* number of reflections at the convex side */
        numi=numa;
        if (ausmWin*einmWin<0)
        {
           if (ausmWin-einmWin>0)
              numi++;
           else
              numi--;
        }

        /* is the reflection coefficient too small? */
        if (((numa>0) && (aeuref<=0)) || ((numi>0) && (innref<=0)))
           ABSORB;

        /* calculation of the neutron probability weight p */
        for (i=1;i<=numa;i++)
            p*=aeuref;
        for (i=1;i<=numi;i++)
            p*=innref;

        /* time to cross the bender */
        Ta=(2*numa*(tan(aeumWin)-tan(innmWin))
          +tan(ausmWin)-tan(einmWin)
          -tan(innmWin)*(SIGN(ausmWin)-SIGN(einmWin)))
          *(R-dab)/vpl;
        t+=Ta;

        /* distance between neutron and concave side of channel at the exit */
        ab=R-(R-dab)/cos(ausmWin);

        /* calculation of the exit coordinates in the XZ-plane */
        x=w/2.0-ab-dru;
        z=r*mWin;
        vx=sin(ausmWin)*vpl;
        vz=cos(ausmWin)*vpl;

        /*** reflections at top and bottom side (Y axis) ***/

        if (vy!=0.0)
        {
          double s[] = {R0s, Qcs, alphas, ms, Ws};
          /* reflection coefficent at the top and bottom side */
          Q=2.0*V2K*fabs(vy);
          StdReflecFunc(Q, s, &ref);

          /* number of reflections at top and bottom */
          einzei=h/2.0/fabs(vy)+y/vy;
          zykzei=h/fabs(vy);
          num=(Ta+einzei)/zykzei;

          /* time between the last reflection and the exit */
          auszei=fmod(Ta+einzei,zykzei);

          /* is the reflection coefficient too small? */
          if ((num>0) && (ref<=0))
             ABSORB;

          /* calculation of the probability weight p */
          for (i=1;i<=num;i++) {
               p*=ref;
               vy*=-1.0; }

          /* calculation of the exit coordinate */
          y=auszei*vy-vy*h/fabs(vy)/2.0;
        } /* if (vy!=0.0) */
        SCATTER;
      } /* if (dab>bk-d)  */
      else
        ABSORB; /* hit separating walls */
    }
    else /* if ((fabs(x)<w/2) && (fabs(y)<h/2))   */
      ABSORB; /* miss entry window */

}
#line 16848 "./PSI_DMC.c"
}   /* End of guide2=Bender() SETTING parameter declarations. */
#undef mWin
#undef bk
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompguide2:
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

  /* TRACE Component PSDafter_curve [8] */
  mccoordschange(mcposrPSDafter_curve, mcrotrPSDafter_curve,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component PSDafter_curve (without coords transformations) */
  mcJumpTrace_PSDafter_curve:
  SIG_MESSAGE("PSDafter_curve (Trace)");
  mcDEBUG_COMP("PSDafter_curve")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompPSDafter_curve
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
#define mccompcurname  PSDafter_curve
#define mccompcurtype  PSD_monitor
#define mccompcurindex 8
#define PSD_N mccPSDafter_curve_PSD_N
#define PSD_p mccPSDafter_curve_PSD_p
#define PSD_p2 mccPSDafter_curve_PSD_p2
{   /* Declarations of PSDafter_curve=PSD_monitor() SETTING parameters. */
int nx = mccPSDafter_curve_nx;
int ny = mccPSDafter_curve_ny;
char* filename = mccPSDafter_curve_filename;
MCNUM xmin = mccPSDafter_curve_xmin;
MCNUM xmax = mccPSDafter_curve_xmax;
MCNUM ymin = mccPSDafter_curve_ymin;
MCNUM ymax = mccPSDafter_curve_ymax;
MCNUM xwidth = mccPSDafter_curve_xwidth;
MCNUM yheight = mccPSDafter_curve_yheight;
MCNUM restore_neutron = mccPSDafter_curve_restore_neutron;
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
#line 16984 "./PSI_DMC.c"
}   /* End of PSDafter_curve=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompPSDafter_curve:
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

  /* TRACE Component bunker [9] */
  mccoordschange(mcposrbunker, mcrotrbunker,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component bunker (without coords transformations) */
  mcJumpTrace_bunker:
  SIG_MESSAGE("bunker (Trace)");
  mcDEBUG_COMP("bunker")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompbunker
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
#define mccompcurname  bunker
#define mccompcurtype  Guide
#define mccompcurindex 9
#define pTable mccbunker_pTable
{   /* Declarations of bunker=Guide() SETTING parameters. */
char* reflect = mccbunker_reflect;
MCNUM w1 = mccbunker_w1;
MCNUM h1 = mccbunker_h1;
MCNUM w2 = mccbunker_w2;
MCNUM h2 = mccbunker_h2;
MCNUM l = mccbunker_l;
MCNUM R0 = mccbunker_R0;
MCNUM Qc = mccbunker_Qc;
MCNUM alpha = mccbunker_alpha;
MCNUM m = mccbunker_m;
MCNUM W = mccbunker_W;
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
#line 17212 "./PSI_DMC.c"
}   /* End of bunker=Guide() SETTING parameter declarations. */
#undef pTable
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompbunker:
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

  /* TRACE Component guide3 [10] */
  mccoordschange(mcposrguide3, mcrotrguide3,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component guide3 (without coords transformations) */
  mcJumpTrace_guide3:
  SIG_MESSAGE("guide3 (Trace)");
  mcDEBUG_COMP("guide3")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompguide3
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
#define mccompcurname  guide3
#define mccompcurtype  Guide
#define mccompcurindex 10
#define pTable mccguide3_pTable
{   /* Declarations of guide3=Guide() SETTING parameters. */
char* reflect = mccguide3_reflect;
MCNUM w1 = mccguide3_w1;
MCNUM h1 = mccguide3_h1;
MCNUM w2 = mccguide3_w2;
MCNUM h2 = mccguide3_h2;
MCNUM l = mccguide3_l;
MCNUM R0 = mccguide3_R0;
MCNUM Qc = mccguide3_Qc;
MCNUM alpha = mccguide3_alpha;
MCNUM m = mccguide3_m;
MCNUM W = mccguide3_W;
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
#line 17438 "./PSI_DMC.c"
}   /* End of guide3=Guide() SETTING parameter declarations. */
#undef pTable
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompguide3:
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

  /* TRACE Component guide4 [11] */
  mccoordschange(mcposrguide4, mcrotrguide4,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component guide4 (without coords transformations) */
  mcJumpTrace_guide4:
  SIG_MESSAGE("guide4 (Trace)");
  mcDEBUG_COMP("guide4")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompguide4
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
#define mccompcurname  guide4
#define mccompcurtype  Guide
#define mccompcurindex 11
#define pTable mccguide4_pTable
{   /* Declarations of guide4=Guide() SETTING parameters. */
char* reflect = mccguide4_reflect;
MCNUM w1 = mccguide4_w1;
MCNUM h1 = mccguide4_h1;
MCNUM w2 = mccguide4_w2;
MCNUM h2 = mccguide4_h2;
MCNUM l = mccguide4_l;
MCNUM R0 = mccguide4_R0;
MCNUM Qc = mccguide4_Qc;
MCNUM alpha = mccguide4_alpha;
MCNUM m = mccguide4_m;
MCNUM W = mccguide4_W;
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
#line 17664 "./PSI_DMC.c"
}   /* End of guide4=Guide() SETTING parameter declarations. */
#undef pTable
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompguide4:
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

  /* TRACE Component window1 [12] */
  mccoordschange(mcposrwindow1, mcrotrwindow1,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component window1 (without coords transformations) */
  mcJumpTrace_window1:
  SIG_MESSAGE("window1 (Trace)");
  mcDEBUG_COMP("window1")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompwindow1
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
#define mccompcurname  window1
#define mccompcurtype  Al_window
#define mccompcurindex 12
{   /* Declarations of window1=Al_window() SETTING parameters. */
MCNUM thickness = mccwindow1_thickness;
#line 68 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../contrib/Al_window.comp"
{
  double v;                     /* Neutron velocity */
  double dt0;     /* Flight times through sample */
  double dist;
  double Al_s_tot_lambda,Al_my_tot,Al_my_a ; /* total XS (barn), total scattering length (m-1), absorption scat. length */
  double lambda; /* neutrons wavelength */

  PROP_Z0;

  dt0=thickness/vz;
  v=sqrt(vx*vx+vy*vy+vz*vz);
  PROP_DT(dt0);
  dist=v*dt0;

  lambda=sqrt(81.81/(VS2E*v*v));
  Al_s_tot_lambda= Al_pf_A+Al_pf_B1*lambda+ Al_pf_B2*lambda*lambda+ Al_pf_B3*lambda*lambda*lambda;
  Al_s_tot_lambda+=Al_pf_B4*lambda*lambda*lambda*lambda;
  Al_my_tot=Al_rho / Al_mmol * Al_s_tot_lambda * avogadro * 10;
  Al_my_a = Al_my_a_v/v;

  p *=exp(-Al_my_a*dist);/* neutron passes window without any interaction */

  /* TODO: scatter in Debye-Scherrer cone */

}
#line 17798 "./PSI_DMC.c"
}   /* End of window1=Al_window() SETTING parameter declarations. */
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompwindow1:
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

  /* TRACE Component ydist_fluxpos [13] */
  mccoordschange(mcposrydist_fluxpos, mcrotrydist_fluxpos,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component ydist_fluxpos (without coords transformations) */
  mcJumpTrace_ydist_fluxpos:
  SIG_MESSAGE("ydist_fluxpos (Trace)");
  mcDEBUG_COMP("ydist_fluxpos")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompydist_fluxpos
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
#define mccompcurname  ydist_fluxpos
#define mccompcurtype  PSDlin_monitor
#define mccompcurindex 13
#define nx mccydist_fluxpos_nx
#define PSDlin_N mccydist_fluxpos_PSDlin_N
#define PSDlin_p mccydist_fluxpos_PSDlin_p
#define PSDlin_p2 mccydist_fluxpos_PSDlin_p2
{   /* Declarations of ydist_fluxpos=PSDlin_monitor() SETTING parameters. */
char* filename = mccydist_fluxpos_filename;
MCNUM xmin = mccydist_fluxpos_xmin;
MCNUM xmax = mccydist_fluxpos_xmax;
MCNUM ymin = mccydist_fluxpos_ymin;
MCNUM ymax = mccydist_fluxpos_ymax;
MCNUM xwidth = mccydist_fluxpos_xwidth;
MCNUM yheight = mccydist_fluxpos_yheight;
MCNUM restore_neutron = mccydist_fluxpos_restore_neutron;
int nowritefile = mccydist_fluxpos_nowritefile;
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
#line 17938 "./PSI_DMC.c"
}   /* End of ydist_fluxpos=PSDlin_monitor() SETTING parameter declarations. */
#undef PSDlin_p2
#undef PSDlin_p
#undef PSDlin_N
#undef nx
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompydist_fluxpos:
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

  /* TRACE Component PSD_fluxpos [14] */
  mccoordschange(mcposrPSD_fluxpos, mcrotrPSD_fluxpos,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component PSD_fluxpos (without coords transformations) */
  mcJumpTrace_PSD_fluxpos:
  SIG_MESSAGE("PSD_fluxpos (Trace)");
  mcDEBUG_COMP("PSD_fluxpos")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompPSD_fluxpos
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
#define mccompcurname  PSD_fluxpos
#define mccompcurtype  PSD_monitor
#define mccompcurindex 14
#define PSD_N mccPSD_fluxpos_PSD_N
#define PSD_p mccPSD_fluxpos_PSD_p
#define PSD_p2 mccPSD_fluxpos_PSD_p2
{   /* Declarations of PSD_fluxpos=PSD_monitor() SETTING parameters. */
int nx = mccPSD_fluxpos_nx;
int ny = mccPSD_fluxpos_ny;
char* filename = mccPSD_fluxpos_filename;
MCNUM xmin = mccPSD_fluxpos_xmin;
MCNUM xmax = mccPSD_fluxpos_xmax;
MCNUM ymin = mccPSD_fluxpos_ymin;
MCNUM ymax = mccPSD_fluxpos_ymax;
MCNUM xwidth = mccPSD_fluxpos_xwidth;
MCNUM yheight = mccPSD_fluxpos_yheight;
MCNUM restore_neutron = mccPSD_fluxpos_restore_neutron;
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
#line 18076 "./PSI_DMC.c"
}   /* End of PSD_fluxpos=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompPSD_fluxpos:
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

  /* TRACE Component xdist_flux_pos [15] */
  mccoordschange(mcposrxdist_flux_pos, mcrotrxdist_flux_pos,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component xdist_flux_pos (without coords transformations) */
  mcJumpTrace_xdist_flux_pos:
  SIG_MESSAGE("xdist_flux_pos (Trace)");
  mcDEBUG_COMP("xdist_flux_pos")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompxdist_flux_pos
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
#define mccompcurname  xdist_flux_pos
#define mccompcurtype  PSDlin_monitor
#define mccompcurindex 15
#define nx mccxdist_flux_pos_nx
#define PSDlin_N mccxdist_flux_pos_PSDlin_N
#define PSDlin_p mccxdist_flux_pos_PSDlin_p
#define PSDlin_p2 mccxdist_flux_pos_PSDlin_p2
{   /* Declarations of xdist_flux_pos=PSDlin_monitor() SETTING parameters. */
char* filename = mccxdist_flux_pos_filename;
MCNUM xmin = mccxdist_flux_pos_xmin;
MCNUM xmax = mccxdist_flux_pos_xmax;
MCNUM ymin = mccxdist_flux_pos_ymin;
MCNUM ymax = mccxdist_flux_pos_ymax;
MCNUM xwidth = mccxdist_flux_pos_xwidth;
MCNUM yheight = mccxdist_flux_pos_yheight;
MCNUM restore_neutron = mccxdist_flux_pos_restore_neutron;
int nowritefile = mccxdist_flux_pos_nowritefile;
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
#line 18219 "./PSI_DMC.c"
}   /* End of xdist_flux_pos=PSDlin_monitor() SETTING parameter declarations. */
#undef PSDlin_p2
#undef PSDlin_p
#undef PSDlin_N
#undef nx
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompxdist_flux_pos:
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

  /* TRACE Component PSD_fluxposB [16] */
  mccoordschange(mcposrPSD_fluxposB, mcrotrPSD_fluxposB,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component PSD_fluxposB (without coords transformations) */
  mcJumpTrace_PSD_fluxposB:
  SIG_MESSAGE("PSD_fluxposB (Trace)");
  mcDEBUG_COMP("PSD_fluxposB")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompPSD_fluxposB
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
#define mccompcurname  PSD_fluxposB
#define mccompcurtype  PSD_monitor
#define mccompcurindex 16
#define PSD_N mccPSD_fluxposB_PSD_N
#define PSD_p mccPSD_fluxposB_PSD_p
#define PSD_p2 mccPSD_fluxposB_PSD_p2
{   /* Declarations of PSD_fluxposB=PSD_monitor() SETTING parameters. */
int nx = mccPSD_fluxposB_nx;
int ny = mccPSD_fluxposB_ny;
char* filename = mccPSD_fluxposB_filename;
MCNUM xmin = mccPSD_fluxposB_xmin;
MCNUM xmax = mccPSD_fluxposB_xmax;
MCNUM ymin = mccPSD_fluxposB_ymin;
MCNUM ymax = mccPSD_fluxposB_ymax;
MCNUM xwidth = mccPSD_fluxposB_xwidth;
MCNUM yheight = mccPSD_fluxposB_yheight;
MCNUM restore_neutron = mccPSD_fluxposB_restore_neutron;
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
#line 18357 "./PSI_DMC.c"
}   /* End of PSD_fluxposB=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompPSD_fluxposB:
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

  /* TRACE Component window2 [17] */
  mccoordschange(mcposrwindow2, mcrotrwindow2,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component window2 (without coords transformations) */
  mcJumpTrace_window2:
  SIG_MESSAGE("window2 (Trace)");
  mcDEBUG_COMP("window2")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompwindow2
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
#define mccompcurname  window2
#define mccompcurtype  Al_window
#define mccompcurindex 17
{   /* Declarations of window2=Al_window() SETTING parameters. */
MCNUM thickness = mccwindow2_thickness;
#line 68 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../contrib/Al_window.comp"
{
  double v;                     /* Neutron velocity */
  double dt0;     /* Flight times through sample */
  double dist;
  double Al_s_tot_lambda,Al_my_tot,Al_my_a ; /* total XS (barn), total scattering length (m-1), absorption scat. length */
  double lambda; /* neutrons wavelength */

  PROP_Z0;

  dt0=thickness/vz;
  v=sqrt(vx*vx+vy*vy+vz*vz);
  PROP_DT(dt0);
  dist=v*dt0;

  lambda=sqrt(81.81/(VS2E*v*v));
  Al_s_tot_lambda= Al_pf_A+Al_pf_B1*lambda+ Al_pf_B2*lambda*lambda+ Al_pf_B3*lambda*lambda*lambda;
  Al_s_tot_lambda+=Al_pf_B4*lambda*lambda*lambda*lambda;
  Al_my_tot=Al_rho / Al_mmol * Al_s_tot_lambda * avogadro * 10;
  Al_my_a = Al_my_a_v/v;

  p *=exp(-Al_my_a*dist);/* neutron passes window without any interaction */

  /* TODO: scatter in Debye-Scherrer cone */

}
#line 18493 "./PSI_DMC.c"
}   /* End of window2=Al_window() SETTING parameter declarations. */
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompwindow2:
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

  /* TRACE Component in_slit [18] */
  mccoordschange(mcposrin_slit, mcrotrin_slit,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component in_slit (without coords transformations) */
  mcJumpTrace_in_slit:
  SIG_MESSAGE("in_slit (Trace)");
  mcDEBUG_COMP("in_slit")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompin_slit
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
#define mccompcurname  in_slit
#define mccompcurtype  Slit
#define mccompcurindex 18
{   /* Declarations of in_slit=Slit() SETTING parameters. */
MCNUM xmin = mccin_slit_xmin;
MCNUM xmax = mccin_slit_xmax;
MCNUM ymin = mccin_slit_ymin;
MCNUM ymax = mccin_slit_ymax;
MCNUM radius = mccin_slit_radius;
MCNUM xwidth = mccin_slit_xwidth;
MCNUM yheight = mccin_slit_yheight;
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
#line 18617 "./PSI_DMC.c"
}   /* End of in_slit=Slit() SETTING parameter declarations. */
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompin_slit:
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

  /* TRACE Component lambda_in [19] */
  mccoordschange(mcposrlambda_in, mcrotrlambda_in,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component lambda_in (without coords transformations) */
  mcJumpTrace_lambda_in:
  SIG_MESSAGE("lambda_in (Trace)");
  mcDEBUG_COMP("lambda_in")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbComplambda_in
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
#define mccompcurname  lambda_in
#define mccompcurtype  L_monitor
#define mccompcurindex 19
#define nL mcclambda_in_nL
#define L_N mcclambda_in_L_N
#define L_p mcclambda_in_L_p
#define L_p2 mcclambda_in_L_p2
{   /* Declarations of lambda_in=L_monitor() SETTING parameters. */
char* filename = mcclambda_in_filename;
MCNUM xmin = mcclambda_in_xmin;
MCNUM xmax = mcclambda_in_xmax;
MCNUM ymin = mcclambda_in_ymin;
MCNUM ymax = mcclambda_in_ymax;
MCNUM xwidth = mcclambda_in_xwidth;
MCNUM yheight = mcclambda_in_yheight;
MCNUM Lmin = mcclambda_in_Lmin;
MCNUM Lmax = mcclambda_in_Lmax;
MCNUM restore_neutron = mcclambda_in_restore_neutron;
int nowritefile = mcclambda_in_nowritefile;
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
#line 18760 "./PSI_DMC.c"
}   /* End of lambda_in=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbComplambda_in:
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

  /* TRACE Component sma [20] */
  mccoordschange(mcposrsma, mcrotrsma,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component sma (without coords transformations) */
  mcJumpTrace_sma:
  SIG_MESSAGE("sma (Trace)");
  mcDEBUG_COMP("sma")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompsma
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
#define mccompcurname  sma
#define mccompcurtype  Arm
#define mccompcurindex 20
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompsma:
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

  /* TRACE Component foc_mono [21] */
  mccoordschange(mcposrfoc_mono, mcrotrfoc_mono,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component foc_mono (without coords transformations) */
  mcJumpTrace_foc_mono:
  SIG_MESSAGE("foc_mono (Trace)");
  mcDEBUG_COMP("foc_mono")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompfoc_mono
  if (!mcSplit_foc_mono) {                   /* STORE only the first time */
    if (floor(10) > 1) p /= floor(10); /* adapt weight for SPLITed neutron */
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
  } else {
    RESTORE_NEUTRON(21,
      mcnlx,
      mcnly,
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
  mcSplit_foc_mono++; /* SPLIT number */
  mcScattered=0;
  mcRestore=0;
  mcNCounter[21]++;
  mcPCounter[21] += p;
  mcP2Counter[21] += p*p;
#define mccompcurname  foc_mono
#define mccompcurtype  Monochromator_2foc
#define mccompcurindex 21
#define mos_y mccfoc_mono_mos_y
#define mos_z mccfoc_mono_mos_z
#define mono_Q mccfoc_mono_mono_Q
#define SlabWidth mccfoc_mono_SlabWidth
#define SlabHeight mccfoc_mono_SlabHeight
#define rTable mccfoc_mono_rTable
{   /* Declarations of foc_mono=Monochromator_2foc() SETTING parameters. */
char* reflect = mccfoc_mono_reflect;
MCNUM zwidth = mccfoc_mono_zwidth;
MCNUM yheight = mccfoc_mono_yheight;
MCNUM gap = mccfoc_mono_gap;
MCNUM NH = mccfoc_mono_NH;
MCNUM NV = mccfoc_mono_NV;
MCNUM mosaich = mccfoc_mono_mosaich;
MCNUM mosaicv = mccfoc_mono_mosaicv;
MCNUM r0 = mccfoc_mono_r0;
MCNUM Q = mccfoc_mono_Q;
MCNUM RV = mccfoc_mono_RV;
MCNUM RH = mccfoc_mono_RH;
MCNUM DM = mccfoc_mono_DM;
MCNUM mosaic = mccfoc_mono_mosaic;
MCNUM width = mccfoc_mono_width;
MCNUM height = mccfoc_mono_height;
MCNUM verbose = mccfoc_mono_verbose;
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
#line 19134 "./PSI_DMC.c"
}   /* End of foc_mono=Monochromator_2foc() SETTING parameter declarations. */
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
  mcabsorbCompfoc_mono:
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

  /* TRACE Component msa [22] */
  mccoordschange(mcposrmsa, mcrotrmsa,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component msa (without coords transformations) */
  mcJumpTrace_msa:
  SIG_MESSAGE("msa (Trace)");
  mcDEBUG_COMP("msa")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompmsa
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
#define mccompcurname  msa
#define mccompcurtype  Arm
#define mccompcurindex 22
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompmsa:
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

  /* TRACE Component out1_slit [23] */
  mccoordschange(mcposrout1_slit, mcrotrout1_slit,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component out1_slit (without coords transformations) */
  mcJumpTrace_out1_slit:
  SIG_MESSAGE("out1_slit (Trace)");
  mcDEBUG_COMP("out1_slit")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompout1_slit
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
#define mccompcurname  out1_slit
#define mccompcurtype  Slit
#define mccompcurindex 23
{   /* Declarations of out1_slit=Slit() SETTING parameters. */
MCNUM xmin = mccout1_slit_xmin;
MCNUM xmax = mccout1_slit_xmax;
MCNUM ymin = mccout1_slit_ymin;
MCNUM ymax = mccout1_slit_ymax;
MCNUM radius = mccout1_slit_radius;
MCNUM xwidth = mccout1_slit_xwidth;
MCNUM yheight = mccout1_slit_yheight;
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
#line 19367 "./PSI_DMC.c"
}   /* End of out1_slit=Slit() SETTING parameter declarations. */
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompout1_slit:
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

  /* TRACE Component Amoin_slit [24] */
  mccoordschange(mcposrAmoin_slit, mcrotrAmoin_slit,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component Amoin_slit (without coords transformations) */
  mcJumpTrace_Amoin_slit:
  SIG_MESSAGE("Amoin_slit (Trace)");
  mcDEBUG_COMP("Amoin_slit")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompAmoin_slit
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
#define mccompcurname  Amoin_slit
#define mccompcurtype  Slit
#define mccompcurindex 24
{   /* Declarations of Amoin_slit=Slit() SETTING parameters. */
MCNUM xmin = mccAmoin_slit_xmin;
MCNUM xmax = mccAmoin_slit_xmax;
MCNUM ymin = mccAmoin_slit_ymin;
MCNUM ymax = mccAmoin_slit_ymax;
MCNUM radius = mccAmoin_slit_radius;
MCNUM xwidth = mccAmoin_slit_xwidth;
MCNUM yheight = mccAmoin_slit_yheight;
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
#line 19491 "./PSI_DMC.c"
}   /* End of Amoin_slit=Slit() SETTING parameter declarations. */
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompAmoin_slit:
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

  /* TRACE Component Bmoin_slit [25] */
  mccoordschange(mcposrBmoin_slit, mcrotrBmoin_slit,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component Bmoin_slit (without coords transformations) */
  mcJumpTrace_Bmoin_slit:
  SIG_MESSAGE("Bmoin_slit (Trace)");
  mcDEBUG_COMP("Bmoin_slit")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompBmoin_slit
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
#define mccompcurname  Bmoin_slit
#define mccompcurtype  Slit
#define mccompcurindex 25
{   /* Declarations of Bmoin_slit=Slit() SETTING parameters. */
MCNUM xmin = mccBmoin_slit_xmin;
MCNUM xmax = mccBmoin_slit_xmax;
MCNUM ymin = mccBmoin_slit_ymin;
MCNUM ymax = mccBmoin_slit_ymax;
MCNUM radius = mccBmoin_slit_radius;
MCNUM xwidth = mccBmoin_slit_xwidth;
MCNUM yheight = mccBmoin_slit_yheight;
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
#line 19615 "./PSI_DMC.c"
}   /* End of Bmoin_slit=Slit() SETTING parameter declarations. */
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompBmoin_slit:
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

  /* TRACE Component out2_slit [26] */
  mccoordschange(mcposrout2_slit, mcrotrout2_slit,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component out2_slit (without coords transformations) */
  mcJumpTrace_out2_slit:
  SIG_MESSAGE("out2_slit (Trace)");
  mcDEBUG_COMP("out2_slit")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompout2_slit
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
#define mccompcurname  out2_slit
#define mccompcurtype  Slit
#define mccompcurindex 26
{   /* Declarations of out2_slit=Slit() SETTING parameters. */
MCNUM xmin = mccout2_slit_xmin;
MCNUM xmax = mccout2_slit_xmax;
MCNUM ymin = mccout2_slit_ymin;
MCNUM ymax = mccout2_slit_ymax;
MCNUM radius = mccout2_slit_radius;
MCNUM xwidth = mccout2_slit_xwidth;
MCNUM yheight = mccout2_slit_yheight;
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
#line 19739 "./PSI_DMC.c"
}   /* End of out2_slit=Slit() SETTING parameter declarations. */
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompout2_slit:
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

  /* TRACE Component PSD_sample [27] */
  mccoordschange(mcposrPSD_sample, mcrotrPSD_sample,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component PSD_sample (without coords transformations) */
  mcJumpTrace_PSD_sample:
  SIG_MESSAGE("PSD_sample (Trace)");
  mcDEBUG_COMP("PSD_sample")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompPSD_sample
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
#define mccompcurname  PSD_sample
#define mccompcurtype  PSD_monitor
#define mccompcurindex 27
#define PSD_N mccPSD_sample_PSD_N
#define PSD_p mccPSD_sample_PSD_p
#define PSD_p2 mccPSD_sample_PSD_p2
{   /* Declarations of PSD_sample=PSD_monitor() SETTING parameters. */
int nx = mccPSD_sample_nx;
int ny = mccPSD_sample_ny;
char* filename = mccPSD_sample_filename;
MCNUM xmin = mccPSD_sample_xmin;
MCNUM xmax = mccPSD_sample_xmax;
MCNUM ymin = mccPSD_sample_ymin;
MCNUM ymax = mccPSD_sample_ymax;
MCNUM xwidth = mccPSD_sample_xwidth;
MCNUM yheight = mccPSD_sample_yheight;
MCNUM restore_neutron = mccPSD_sample_restore_neutron;
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
#line 19873 "./PSI_DMC.c"
}   /* End of PSD_sample=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompPSD_sample:
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

  /* TRACE Component lambda_sample [28] */
  mccoordschange(mcposrlambda_sample, mcrotrlambda_sample,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component lambda_sample (without coords transformations) */
  mcJumpTrace_lambda_sample:
  SIG_MESSAGE("lambda_sample (Trace)");
  mcDEBUG_COMP("lambda_sample")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbComplambda_sample
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
#define mccompcurname  lambda_sample
#define mccompcurtype  L_monitor
#define mccompcurindex 28
#define nL mcclambda_sample_nL
#define L_N mcclambda_sample_L_N
#define L_p mcclambda_sample_L_p
#define L_p2 mcclambda_sample_L_p2
{   /* Declarations of lambda_sample=L_monitor() SETTING parameters. */
char* filename = mcclambda_sample_filename;
MCNUM xmin = mcclambda_sample_xmin;
MCNUM xmax = mcclambda_sample_xmax;
MCNUM ymin = mcclambda_sample_ymin;
MCNUM ymax = mcclambda_sample_ymax;
MCNUM xwidth = mcclambda_sample_xwidth;
MCNUM yheight = mcclambda_sample_yheight;
MCNUM Lmin = mcclambda_sample_Lmin;
MCNUM Lmax = mcclambda_sample_Lmax;
MCNUM restore_neutron = mcclambda_sample_restore_neutron;
int nowritefile = mcclambda_sample_nowritefile;
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
#line 20019 "./PSI_DMC.c"
}   /* End of lambda_sample=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbComplambda_sample:
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

  /* TRACE Component sa_arm [29] */
  mccoordschange(mcposrsa_arm, mcrotrsa_arm,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component sa_arm (without coords transformations) */
  mcJumpTrace_sa_arm:
  SIG_MESSAGE("sa_arm (Trace)");
  mcDEBUG_COMP("sa_arm")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompsa_arm
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
#define mccompcurname  sa_arm
#define mccompcurtype  Arm
#define mccompcurindex 29
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompsa_arm:
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

  /* TRACE Component sample [30] */
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
  } else {
    RESTORE_NEUTRON(30,
      mcnlx,
      mcnly,
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
  mcNCounter[30]++;
  mcPCounter[30] += p;
  mcP2Counter[30] += p*p;
#define mccompcurname  sample
#define mccompcurtype  PowderN
#define mccompcurindex 30
#define format mccsample_format
#define line_info mccsample_line_info
#define columns mccsample_columns
#define offdata mccsample_offdata
{   /* Declarations of sample=PowderN() SETTING parameters. */
char* reflections = mccsample_reflections;
char* geometry = mccsample_geometry;
MCNUM radius = mccsample_radius;
MCNUM yheight = mccsample_yheight;
MCNUM xwidth = mccsample_xwidth;
MCNUM zdepth = mccsample_zdepth;
MCNUM thickness = mccsample_thickness;
MCNUM pack = mccsample_pack;
MCNUM Vc = mccsample_Vc;
MCNUM sigma_abs = mccsample_sigma_abs;
MCNUM sigma_inc = mccsample_sigma_inc;
MCNUM delta_d_d = mccsample_delta_d_d;
MCNUM p_inc = mccsample_p_inc;
MCNUM p_transmit = mccsample_p_transmit;
MCNUM DW = mccsample_DW;
MCNUM nb_atoms = mccsample_nb_atoms;
MCNUM d_omega = mccsample_d_omega;
MCNUM d_phi = mccsample_d_phi;
MCNUM tth_sign = mccsample_tth_sign;
MCNUM p_interact = mccsample_p_interact;
MCNUM concentric = mccsample_concentric;
MCNUM density = mccsample_density;
MCNUM weight = mccsample_weight;
MCNUM barns = mccsample_barns;
MCNUM Strain = mccsample_Strain;
MCNUM focus_flip = mccsample_focus_flip;
int target_index = mccsample_target_index;
#line 774 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../samples/PowderN.comp"
{
  double t0, t1, t2, t3, v, v1,l_full, l, l_1, dt, alpha0, alpha, theta, my_s, my_s_n, sg;
  double solid_angle, neutrontype;
  double arg, tmp_vx, tmp_vy, tmp_vz, vout_x, vout_y, vout_z, nx, ny, nz, pmul=1;
  int    line;
  char   intersect=0;
  char   intersecti=0;
  line_info.type = '\0';
  line_info.itype = 0;

  if (line_info.V_0 > 0 && (line_info.count || line_info.my_inc)) {
    if (line_info.shape == 1) {
      intersect  = box_intersect(&t0, &t3, x, y, z, vx, vy, vz, xwidth, yheight, zdepth);
      intersecti = box_intersect(&t1, &t2, x, y, z, vx, vy, vz, line_info.xwidth_i, line_info.yheight_i, line_info.zdepth_i);
    } else if (line_info.shape == 0) {
      intersect  = cylinder_intersect(&t0, &t3, x, y, z, vx, vy, vz, radius, yheight);
      intersecti = cylinder_intersect(&t1, &t2, x, y, z, vx, vy, vz, line_info.radius_i, line_info.yheight_i);
    } else if (line_info.shape == 2) {
      intersect  = sphere_intersect  (&t0, &t3, x,y,z, vx,vy,vz, radius);
      intersecti = sphere_intersect  (&t1, &t2, x,y,z, vx,vy,vz, line_info.radius_i);
    } else if (line_info.shape == 3) {
      intersect  = off_intersect  (&t0, &t3, NULL, NULL, x,y,z, vx,vy,vz, offdata);
      intersecti = 0;
    }
  }

  if(intersect && t3 >0) {

    if (concentric) {
      /* Set up for concentric case */
      /* 'Remove' the backside of this comp */
      if (!intersecti) {
        t1 = (t3 + t0) /2;
      }
      t2 = t1;
      t3 = t1;
      dt = -1.0*rand01(); /* In case of scattering we will scatter on 'forward' part of sample */
    } else {
      if (!intersecti) {
        t1 = (t3 + t0) /2;
        t2 = t1;
      }
      dt = randpm1(); /* Possibility to scatter at all points in line of sight */
    }

    /* Neutron enters at t=t0. */
    if(t0 < 0) t0=0; /* already in sample */
    if(t1 < 0) t1=0; /* already in inner hollow */
    if(t2 < 0) t2=0; /* already past inner hollow */
    v = sqrt(vx*vx + vy*vy + vz*vz);
    l_full = v * (t3 - t2 + t1 - t0);

    if (line_info.neutron_passed < CHAR_BUF_LENGTH) {
      if (v < line_info.v_min) line_info.v_min = v;
      if (v > line_info.v_max) line_info.v_max = v;
      line_info.neutron_passed++;
    }

    /* Calculate total scattering cross section at relevant velocity */
    if ( fabs(v - line_info.v) < 1e-6) {
        line_info.nb_reuses++;
      } else {
        line_info.Nq = calc_xsect(v, line_info.q_v, line_info.my_s_v2, line_info.count, &line_info.my_s_v2_sum, &line_info);
        line_info.v = v;
        line_info.nb_refl += line_info.Nq;
        line_info.nb_refl_count++;
      }

    if (t3 < 0) {
      t3=0; /* Already past sample?! */
      if (line_info.flag_warning < 100)
      printf("PowderN: %s: Warning: Neutron has already passed us? (Skipped).\n"
             "         In concentric geometry, this may be caused by a missing concentric=0 option in 2nd enclosing instance.\n", NAME_CURRENT_COMP);
      line_info.flag_warning++;
    } else {
      if (dt<0) { /* Calculate scattering point position */
        dt = fabs(dt)*(t1 - t0); /* 'Forward' part */
      } else {
        dt = dt * (t3 - t2) + (t2-t0) ; /* Possibly also 'backside' part */
      }

      my_s = line_info.my_s_v2_sum/(v*v)+line_info.my_inc;
      /* Total attenuation from scattering */
	  line_info.lfree=0;
      neutrontype = rand01();
      /* How to handle this one? Transmit (1) / Incoherent (2) / Coherent (3) ? */
      if (neutrontype < p_transmit) {
        neutrontype = 1;
        l = l_full; /* Passing through, full length */
        PROP_DT(t3);
      } else if (neutrontype >= p_transmit && neutrontype < (p_transmit + p_inc)) {
        neutrontype = 2;
        l = v*dt;       /* Penetration in sample */
        PROP_DT(dt+t0); /* Point of scattering */
        SCATTER;
      } else if (neutrontype >= p_transmit + p_inc) {
        neutrontype = 3;
        l = v*dt;       /* Penetration in sample */
        PROP_DT(dt+t0); /* Point of scattering */
        SCATTER;
      } else {
        exit(fprintf(stderr,"PowderN %s: DEAD - this shouldn't happen!\n", NAME_CURRENT_COMP));
      }

      if (neutrontype == 3) { /* Make coherent scattering event */
        if (line_info.count > 0) {
          /* choose line */
			if (line_info.Nq > 1) line=floor(line_info.Nq*rand01());  /* Select between Nq powder lines */
			else line = 0;
			if (line_info.w_v[line])
				arg = line_info.q_v[line]*(1+line_info.w_v[line]*randnorm())/(2.0*v);
			else
				arg = line_info.q_v[line]/(2.0*v);
			my_s_n = line_info.my_s_v2[line]/(v*v);
			if(fabs(arg) > 1)
				ABSORB;                   /* No bragg scattering possible*/
			if (tth_sign == 0) {
				sg = randpm1();
				if (sg > 0) sg = 1; else sg=-1;
			}
			else
				sg = tth_sign/fabs(tth_sign);
			theta = asin(arg);          /* Bragg scattering law */
			  /* Choose point on Debye-Scherrer cone */
			if (d_phi)
			{ /* relate height of detector to the height on DS cone */
				arg = sin(d_phi*DEG2RAD/2)/sin(2*theta);
				/* If full Debye-Scherrer cone is within d_phi, don't focus */
				if (arg < -1 || arg > 1) d_phi = 0;
				/* Otherwise, determine alpha to rotate from scattering plane
				   into d_phi focusing area*/
				else alpha = 2*asin(arg);
			}
			if (d_phi) {
				/* Focusing */
				alpha = fabs(alpha);
				alpha0 = 0.5*randpm1()*alpha;
				if(focus_flip){
					alpha0+=M_PI_2;
				}
			}
			else
				alpha0 = PI*randpm1();

          /* now find a nearly vertical rotation axis:
           * Either
           *  (v along Z) x (X axis) -> nearly Y axis
           * Or
           *  (v along X) x (Z axis) -> nearly Y axis
           */
		   
  /* update JS, 1/7/2017
    If a target is defined, try to define vertical axis as a normal to the plane 
	defined by the incident neutron velocity and target position. 
	Check that v is not ~ parallel to the target direction.
  */
			double vnorm=0.0;
			if (target_index) {
				vec_prod(tmp_vx, tmp_vy, tmp_vz, vx,vy,vz, tgt_x, tgt_y, tgt_z);
				vnorm = sqrt(tmp_vx*tmp_vx+tmp_vy*tmp_vy+tmp_vz*tmp_vz)/v;
			}
			// no target or direction is nearly parallel to v:
			if (vnorm<0.01) {
				if (fabs(vx/v) < fabs(vz/v)) {
					nx = 1; ny = 0; nz = 0;
				} else {
					nx = 0; ny = 0; nz = 1;
				}
				vec_prod(tmp_vx,tmp_vy,tmp_vz, vx,vy,vz, nx,ny,nz);
			}

			  /* v_out = rotate 'v' by 2*theta around tmp_v: Bragg angle */
			rotate(vout_x,vout_y,vout_z, vx,vy,vz, 2*sg*theta, tmp_vx,tmp_vy,tmp_vz);
			
			/* tmp_v = rotate v_out by alpha0 around 'v' (Debye-Scherrer cone) */
			rotate(tmp_vx,tmp_vy,tmp_vz, vout_x,vout_y,vout_z, alpha0, vx, vy, vz);
			vx = tmp_vx;
			vy = tmp_vy;
			vz = tmp_vz;
			  
			  /* Since now scattered and new direction given, calculate path to exit */
			if (line_info.shape == 1) {
				intersect  = box_intersect(&t0, &t3, x, y, z, vx, vy, vz, xwidth, yheight, zdepth);
				intersecti = box_intersect(&t1, &t2, x, y, z, vx, vy, vz, line_info.xwidth_i, line_info.yheight_i, line_info.zdepth_i);
			} else if (line_info.shape == 0) {
				intersect  = cylinder_intersect(&t0, &t3, x, y, z, vx, vy, vz, radius, yheight);
				intersecti = cylinder_intersect(&t1, &t2, x, y, z, vx, vy, vz, line_info.radius_i, line_info.yheight_i);
			} else if (line_info.shape == 2) {
				intersect  = sphere_intersect  (&t0, &t3, x,y,z, vx,vy,vz, radius);
				intersecti = sphere_intersect  (&t1, &t2, x,y,z, vx,vy,vz, line_info.radius_i);
			} else if (line_info.shape == 3) {
				intersect  = off_intersect  (&t0, &t3, NULL, NULL, x,y,z, vx,vy,vz, offdata);
				intersecti = 0;
			}

			if (!intersect) {
				/* Strange error: did not hit cylinder */
				if (line_info.flag_warning < 100)
				  printf("PowderN: %s: WARNING: Did not hit sample from inside (coh). ABSORB.\n", NAME_CURRENT_COMP);
				line_info.flag_warning++;
				ABSORB;
			}

			if (!intersecti) {
				t1 = (t3 + t0) /2;
				t2 = t1;
			}

			if (concentric && intersecti) {
				/* In case of concentricity, 'remove' backward wall of sample */
				t2 = t1;
				t3 = t1;
			}

			if(t0 < 0) t0=0; /* already in sample */
			if(t1 < 0) t1=0; /* already in inner hollow */
			if(t2 < 0) t2=0; /* already past inner hollow */


			l_1 = v*(t3 - t2 + t1 - t0); /* Length to exit */

			pmul  = line_info.Nq*l_full*my_s_n*exp(-(line_info.my_a_v/v+my_s)*(l+l_1))
									  /(1-(p_inc+p_transmit));
			  /* Correction in case of d_phi focusing - BUT only when d_phi != 0 */
			if (d_phi) pmul *= alpha/PI;

			line_info.type = 'c';
			line_info.itype = 1;
			line_info.dq = line_info.q_v[line]*V2K;
			line_info.lfree=1/(line_info.my_a_v/v+my_s);
        } /* else transmit <-- No powder lines in file */
      }  /* Coherent scattering event */
      else if (neutrontype == 2) {  /* Make incoherent scattering event */
		if (d_omega && d_phi) {
			randvec_target_rect_angular(&vx, &vy, &vz, &solid_angle,
			  tgt_x, tgt_y, tgt_z, d_omega*DEG2RAD, d_phi*DEG2RAD, ROT_A_CURRENT_COMP);
		} else if (d_phi) {
			randvec_target_rect_angular(&vx, &vy, &vz, &solid_angle,
                                      tgt_x, tgt_y, tgt_z,
                                      2*PI, d_phi*DEG2RAD, ROT_A_CURRENT_COMP);
        } else {
          randvec_target_circle(&vx, &vy, &vz,
                                &solid_angle, 0, 0, 1, 0);
        }
        v1 = sqrt(vx*vx+vy*vy+vz*vz);
        vx *= v/v1;
        vy *= v/v1;
        vz *= v/v1;

        /* Since now scattered and new direction given, calculate path to exit */
        if (line_info.shape == 1) {
          intersect  = box_intersect(&t0, &t3, x, y, z, vx, vy, vz, xwidth, yheight, zdepth);
          intersecti = box_intersect(&t1, &t2, x, y, z, vx, vy, vz, line_info.xwidth_i, line_info.yheight_i, line_info.zdepth_i);
        } else if (line_info.shape == 0) {
          intersect  = cylinder_intersect(&t0, &t3, x, y, z, vx, vy, vz, radius, yheight);
          intersecti = cylinder_intersect(&t1, &t2, x, y, z, vx, vy, vz, line_info.radius_i, line_info.yheight_i);
        } else if (line_info.shape == 2) {
          intersect  = sphere_intersect  (&t0, &t3, x,y,z, vx,vy,vz, radius);
          intersecti = sphere_intersect  (&t1, &t2, x,y,z, vx,vy,vz, line_info.radius_i);
        } else if (line_info.shape == 3) {
          intersect  = off_intersect  (&t0, &t3, NULL, NULL, x,y,z, vx,vy,vz, offdata);
          intersecti = 0;
        }

        if (!intersect) {
          /* Strange error: did not hit cylinder */
          if (line_info.flag_warning < 100)
            printf("PowderN: %s: WARNING: Did not hit sample from inside (inc). ABSORB.\n", NAME_CURRENT_COMP);
          line_info.flag_warning++;
          ABSORB;
        }

        if (!intersecti) {
          t1 = (t3 + t0) /2;
          t2 = t1;
        }

        if (concentric && intersecti) {
          /* In case of concentricity, 'remove' backward wall of sample */
          t2 = t1;
          t3 = t1;
        }

        if(t0 < 0) t0=0; /* already in sample */
        if(t1 < 0) t1=0; /* already in inner hollow */
        if(t2 < 0) t2=0; /* already past inner hollow */


        l_1 = v*(t3 - t2 + t1 - t0); /* Length to exit */

        pmul = l_full*line_info.my_inc*exp(-(line_info.my_a_v/v+my_s)*(l+l_1))/(p_inc);
        pmul *= solid_angle/(4*PI);
		line_info.lfree=1/(line_info.my_a_v/v+my_s);
        line_info.type = 'i';
		line_info.itype = 2;

      }  /* Incoherent scattering event */
      else if (neutrontype == 1) {
        /* Make transmitted (absorption-corrected) event */
        /* No coordinate changes here, simply change neutron weight */
        pmul = exp(-(line_info.my_a_v/v+my_s)*(l))/(p_transmit);
        line_info.lfree=1/(line_info.my_a_v/v+my_s);
        line_info.type = 't';
		line_info.itype = 3;
      }
      p *= pmul;
    } /* Neutron leaving since it has passed already */
  } /* else transmit non interacting neutrons */

}
#line 20591 "./PSI_DMC.c"
}   /* End of sample=PowderN() SETTING parameter declarations. */
#undef offdata
#undef columns
#undef line_info
#undef format
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompsample:
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

  /* TRACE Component STOP [31] */
  mccoordschange(mcposrSTOP, mcrotrSTOP,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component STOP (without coords transformations) */
  mcJumpTrace_STOP:
  SIG_MESSAGE("STOP (Trace)");
  mcDEBUG_COMP("STOP")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompSTOP
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
#define mccompcurname  STOP
#define mccompcurtype  Beamstop
#define mccompcurindex 31
{   /* Declarations of STOP=Beamstop() SETTING parameters. */
MCNUM xmin = mccSTOP_xmin;
MCNUM xmax = mccSTOP_xmax;
MCNUM ymin = mccSTOP_ymin;
MCNUM ymax = mccSTOP_ymax;
MCNUM xwidth = mccSTOP_xwidth;
MCNUM yheight = mccSTOP_yheight;
MCNUM radius = mccSTOP_radius;
#line 59 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Beamstop.comp"
{
    double Time = t;
    ALLOW_BACKPROP;
    PROP_Z0;
    Time = t - Time;
    if (((Time>=0) && ((radius!=0) && (x*x + y*y <= radius*radius)))
    || ((Time>=0) && (radius==0) && (x>xmin && x<xmax && y>ymin && y<ymax)))
      ABSORB;
    else
      RESTORE_NEUTRON(INDEX_CURRENT_COMP, x, y, z, vx, vy, vz, t, sx, sy, sz, p);
}
#line 20720 "./PSI_DMC.c"
}   /* End of STOP=Beamstop() SETTING parameter declarations. */
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompSTOP:
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

  /* TRACE Component Detector [32] */
  mccoordschange(mcposrDetector, mcrotrDetector,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component Detector (without coords transformations) */
  mcJumpTrace_Detector:
  SIG_MESSAGE("Detector (Trace)");
  mcDEBUG_COMP("Detector")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompDetector
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
#define mccompcurname  Detector
#define mccompcurtype  Monitor_nD
#define mccompcurindex 32
#define user1 mccDetector_user1
#define user2 mccDetector_user2
#define user3 mccDetector_user3
#define DEFS mccDetector_DEFS
#define Vars mccDetector_Vars
#define detector mccDetector_detector
#define offdata mccDetector_offdata
{   /* Declarations of Detector=Monitor_nD() SETTING parameters. */
MCNUM xwidth = mccDetector_xwidth;
MCNUM yheight = mccDetector_yheight;
MCNUM zdepth = mccDetector_zdepth;
MCNUM xmin = mccDetector_xmin;
MCNUM xmax = mccDetector_xmax;
MCNUM ymin = mccDetector_ymin;
MCNUM ymax = mccDetector_ymax;
MCNUM zmin = mccDetector_zmin;
MCNUM zmax = mccDetector_zmax;
MCNUM bins = mccDetector_bins;
MCNUM min = mccDetector_min;
MCNUM max = mccDetector_max;
MCNUM restore_neutron = mccDetector_restore_neutron;
MCNUM radius = mccDetector_radius;
char* options = mccDetector_options;
char* filename = mccDetector_filename;
char* geometry = mccDetector_geometry;
char* username1 = mccDetector_username1;
char* username2 = mccDetector_username2;
char* username3 = mccDetector_username3;
int nowritefile = mccDetector_nowritefile;
#line 313 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/Monitor_nD.comp"
{
  double  XY=0;
  double  t0 = 0;
  double  t1 = 0;
  double  pp;
  int     intersect   = 0;
  char    Flag_Restore = 0;

  if (user1 != FLT_MAX) Vars.UserVariable1 = user1;
  if (user2 != FLT_MAX) Vars.UserVariable2 = user2;
  if (user3 != FLT_MAX) Vars.UserVariable3 = user3;

  /* this is done automatically
    STORE_NEUTRON(INDEX_CURRENT_COMP, x, y, z, vx, vy, vz, t, sx, sy, sz, p);
  */

  if (geometry && strlen(geometry) && strcmp(geometry,"0") && strcmp(geometry, "NULL"))
  {
    /* determine intersections with object */
    intersect = off_intersect_all(&t0, &t1, NULL, NULL,
       x,y,z, vx, vy, vz, &offdata );
    if (Vars.Flag_mantid) {
      if(intersect) {
        Vars.OFF_polyidx=(offdata.intersects[offdata.nextintersect]).index;
      } else {
        Vars.OFF_polyidx=-1;
      }
    }
  }
  else if ( (abs(Vars.Flag_Shape) == DEFS.SHAPE_SQUARE)
            || (abs(Vars.Flag_Shape) == DEFS.SHAPE_DISK) ) /* square xy or disk xy */
  {
    // propagate to xy plane and find intersection
    // make sure the event is recoverable afterwards
    t0 = t;
    ALLOW_BACKPROP;
    PROP_Z0;
    if ( (t>=t0) && (z==0.0) ) // forward propagation to xy plane was successful
    {
      if (abs(Vars.Flag_Shape) == DEFS.SHAPE_SQUARE)
      {
        // square xy
        intersect = (x>=Vars.mxmin && x<=Vars.mxmax && y>=Vars.mymin && y<=Vars.mymax);
      }
      else
      {
        // disk xy
        intersect = (SQR(x) + SQR(y)) <= SQR(Vars.Sphere_Radius);
      }
    }
    else
    {
      intersect=0;
    }
  }
  else if (abs(Vars.Flag_Shape) == DEFS.SHAPE_SPHERE) /* sphere */
  {
    intersect = sphere_intersect(&t0, &t1, x, y, z, vx, vy, vz, Vars.Sphere_Radius);
  /*      intersect = (intersect && t0 > 0); */
  }
  else if ((abs(Vars.Flag_Shape) == DEFS.SHAPE_CYLIND) || (abs(Vars.Flag_Shape) == DEFS.SHAPE_BANANA)) /* cylinder */
  {
    intersect = cylinder_intersect(&t0, &t1, x, y, z, vx, vy, vz, Vars.Sphere_Radius, Vars.Cylinder_Height);
  }
  else if (abs(Vars.Flag_Shape) == DEFS.SHAPE_BOX) /* box */
  {
    intersect = box_intersect(&t0, &t1, x, y, z, vx, vy, vz,
                              fabs(Vars.mxmax-Vars.mxmin), fabs(Vars.mymax-Vars.mymin), fabs(Vars.mzmax-Vars.mzmin));
  }
  else if (abs(Vars.Flag_Shape) == DEFS.SHAPE_PREVIOUS) /* previous comp */
  { intersect = 1; }

  if (intersect)
  {
    if ((abs(Vars.Flag_Shape) == DEFS.SHAPE_SPHERE) || (abs(Vars.Flag_Shape) == DEFS.SHAPE_CYLIND)
     || (abs(Vars.Flag_Shape) == DEFS.SHAPE_BOX) || (abs(Vars.Flag_Shape) == DEFS.SHAPE_BANANA)
     || (geometry && strlen(geometry) && strcmp(geometry,"0") && strcmp(geometry, "NULL")) )
    {
      /* check if we have to remove the top/bottom with BANANA shape */
      if ((abs(Vars.Flag_Shape) == DEFS.SHAPE_BANANA) && (intersect != 1)) {
        double y0,y1;
        /* propagate to intersection point as temporary variable to check top/bottom */
        y0 = y+t0*vy;
        y1 = y+t1*vy;
        if (fabs(y0) >= Vars.Cylinder_Height/2*0.99) t0 = t1;
        if (fabs(y1) >= Vars.Cylinder_Height/2*0.99) t1 = t0;
      }
      if (t0 < 0 && t1 > 0)
        t0 = t;  /* neutron was already inside ! */
      if (t1 < 0 && t0 > 0) /* neutron exit before entering !! */
        t1 = t;
      /* t0 is now time of incoming intersection with the detection area */
      if ((Vars.Flag_Shape < 0) && (t1 > 0))
        PROP_DT(t1); /* t1 outgoing beam */
      else
        PROP_DT(t0); /* t0 incoming beam */
      /* Final test if we are on lid / bottom of banana/sphere */
      if (abs(Vars.Flag_Shape) == DEFS.SHAPE_BANANA || abs(Vars.Flag_Shape) == DEFS.SHAPE_SPHERE) {
        if (fabs(y) >= Vars.Cylinder_Height/2*0.99) {
          intersect=0;
          Flag_Restore=1;
        }
      }
    }
  }

  if (intersect)
  {
    /* Now get the data to monitor: current or keep from PreMonitor */
    if (Vars.Flag_UsePreMonitor != 1)
    {
      Vars.cp  = p;
      Vars.cx  = x;
      Vars.cvx = vx;
      Vars.csx = sx;
      Vars.cy  = y;
      Vars.cvy = vy;
      Vars.csy = sy;
      Vars.cz  = z;
      Vars.cvz = vz;
      Vars.csz = sz;
      Vars.ct  = t;
    }

    if ((Vars.He3_pressure > 0) && (t1 != t0) && ((abs(Vars.Flag_Shape) == DEFS.SHAPE_SPHERE) || (abs(Vars.Flag_Shape) == DEFS.SHAPE_CYLIND) || (abs(Vars.Flag_Shape) == DEFS.SHAPE_BOX)))
    {
      XY = exp(-7.417*Vars.He3_pressure*fabs(t1-t0)*2*PI*K2V);
      /* will monitor the absorbed part */
      Vars.cp *= 1-XY;
      /* and modify the neutron weight after monitor, only remains 1-p_detect */
      p *= XY;
    }

    if (Vars.Flag_capture)
    {
      XY = sqrt(Vars.cvx*Vars.cvx+Vars.cvy*Vars.cvy+Vars.cvz*Vars.cvz);
      XY *= V2K;
      if (XY != 0) XY = 2*PI/XY; /* lambda. lambda(2200 m/2) = 1.7985 Angs  */
      Vars.cp *= XY/1.7985;
    }

    pp = Monitor_nD_Trace(&DEFS, &Vars);
    if (pp==0.0)
    { ABSORB;
    }
    else if(pp==1)
    {
      SCATTER;
    }

    if (Vars.Flag_parallel) /* back to neutron state before detection */
      Flag_Restore = 1;
  } /* end if intersection */
  else {
    if (Vars.Flag_Absorb && !Vars.Flag_parallel)
    {
      // restore neutron ray before absorbing for correct mcdisplay
      RESTORE_NEUTRON(INDEX_CURRENT_COMP, x, y, z, vx, vy, vz, t, sx, sy, sz, p);
      ABSORB;
    }
    else Flag_Restore = 1;  /* no intersection, back to previous state */
  }

  if (Flag_Restore)
  {
    RESTORE_NEUTRON(INDEX_CURRENT_COMP, x, y, z, vx, vy, vz, t, sx, sy, sz, p);
  }
}
#line 21023 "./PSI_DMC.c"
}   /* End of Detector=Monitor_nD() SETTING parameter declarations. */
#undef offdata
#undef detector
#undef Vars
#undef DEFS
#undef user3
#undef user2
#undef user1
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompDetector:
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
  if (mcSplit_foc_mono && mcSplit_foc_mono < (10)) {
    goto mcJumpTrace_foc_mono;
  }
    else mcSplit_foc_mono=0;

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

  /* User SAVE code for component 'source_arm'. */
  SIG_MESSAGE("source_arm (Save)");
#define mccompcurname  source_arm
#define mccompcurtype  Progress_bar
#define mccompcurindex 1
#define IntermediateCnts mccsource_arm_IntermediateCnts
#define StartTime mccsource_arm_StartTime
#define EndTime mccsource_arm_EndTime
#define CurrentTime mccsource_arm_CurrentTime
{   /* Declarations of source_arm=Progress_bar() SETTING parameters. */
char* profile = mccsource_arm_profile;
MCNUM percent = mccsource_arm_percent;
MCNUM flag_save = mccsource_arm_flag_save;
MCNUM minutes = mccsource_arm_minutes;
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
#line 21140 "./PSI_DMC.c"
}   /* End of source_arm=Progress_bar() SETTING parameter declarations. */
#undef CurrentTime
#undef EndTime
#undef StartTime
#undef IntermediateCnts
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'PSDbefore_guides'. */
  SIG_MESSAGE("PSDbefore_guides (Save)");
#define mccompcurname  PSDbefore_guides
#define mccompcurtype  PSD_monitor
#define mccompcurindex 3
#define PSD_N mccPSDbefore_guides_PSD_N
#define PSD_p mccPSDbefore_guides_PSD_p
#define PSD_p2 mccPSDbefore_guides_PSD_p2
{   /* Declarations of PSDbefore_guides=PSD_monitor() SETTING parameters. */
int nx = mccPSDbefore_guides_nx;
int ny = mccPSDbefore_guides_ny;
char* filename = mccPSDbefore_guides_filename;
MCNUM xmin = mccPSDbefore_guides_xmin;
MCNUM xmax = mccPSDbefore_guides_xmax;
MCNUM ymin = mccPSDbefore_guides_ymin;
MCNUM ymax = mccPSDbefore_guides_ymax;
MCNUM xwidth = mccPSDbefore_guides_xwidth;
MCNUM yheight = mccPSDbefore_guides_yheight;
MCNUM restore_neutron = mccPSDbefore_guides_restore_neutron;
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
#line 21180 "./PSI_DMC.c"
}   /* End of PSDbefore_guides=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'l_mon_source'. */
  SIG_MESSAGE("l_mon_source (Save)");
#define mccompcurname  l_mon_source
#define mccompcurtype  L_monitor
#define mccompcurindex 4
#define nL mccl_mon_source_nL
#define L_N mccl_mon_source_L_N
#define L_p mccl_mon_source_L_p
#define L_p2 mccl_mon_source_L_p2
{   /* Declarations of l_mon_source=L_monitor() SETTING parameters. */
char* filename = mccl_mon_source_filename;
MCNUM xmin = mccl_mon_source_xmin;
MCNUM xmax = mccl_mon_source_xmax;
MCNUM ymin = mccl_mon_source_ymin;
MCNUM ymax = mccl_mon_source_ymax;
MCNUM xwidth = mccl_mon_source_xwidth;
MCNUM yheight = mccl_mon_source_yheight;
MCNUM Lmin = mccl_mon_source_Lmin;
MCNUM Lmax = mccl_mon_source_Lmax;
MCNUM restore_neutron = mccl_mon_source_restore_neutron;
int nowritefile = mccl_mon_source_nowritefile;
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
#line 21222 "./PSI_DMC.c"
}   /* End of l_mon_source=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'PSDbefore_curve'. */
  SIG_MESSAGE("PSDbefore_curve (Save)");
#define mccompcurname  PSDbefore_curve
#define mccompcurtype  PSD_monitor
#define mccompcurindex 6
#define PSD_N mccPSDbefore_curve_PSD_N
#define PSD_p mccPSDbefore_curve_PSD_p
#define PSD_p2 mccPSDbefore_curve_PSD_p2
{   /* Declarations of PSDbefore_curve=PSD_monitor() SETTING parameters. */
int nx = mccPSDbefore_curve_nx;
int ny = mccPSDbefore_curve_ny;
char* filename = mccPSDbefore_curve_filename;
MCNUM xmin = mccPSDbefore_curve_xmin;
MCNUM xmax = mccPSDbefore_curve_xmax;
MCNUM ymin = mccPSDbefore_curve_ymin;
MCNUM ymax = mccPSDbefore_curve_ymax;
MCNUM xwidth = mccPSDbefore_curve_xwidth;
MCNUM yheight = mccPSDbefore_curve_yheight;
MCNUM restore_neutron = mccPSDbefore_curve_restore_neutron;
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
#line 21262 "./PSI_DMC.c"
}   /* End of PSDbefore_curve=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'PSDafter_curve'. */
  SIG_MESSAGE("PSDafter_curve (Save)");
#define mccompcurname  PSDafter_curve
#define mccompcurtype  PSD_monitor
#define mccompcurindex 8
#define PSD_N mccPSDafter_curve_PSD_N
#define PSD_p mccPSDafter_curve_PSD_p
#define PSD_p2 mccPSDafter_curve_PSD_p2
{   /* Declarations of PSDafter_curve=PSD_monitor() SETTING parameters. */
int nx = mccPSDafter_curve_nx;
int ny = mccPSDafter_curve_ny;
char* filename = mccPSDafter_curve_filename;
MCNUM xmin = mccPSDafter_curve_xmin;
MCNUM xmax = mccPSDafter_curve_xmax;
MCNUM ymin = mccPSDafter_curve_ymin;
MCNUM ymax = mccPSDafter_curve_ymax;
MCNUM xwidth = mccPSDafter_curve_xwidth;
MCNUM yheight = mccPSDafter_curve_yheight;
MCNUM restore_neutron = mccPSDafter_curve_restore_neutron;
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
#line 21301 "./PSI_DMC.c"
}   /* End of PSDafter_curve=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'ydist_fluxpos'. */
  SIG_MESSAGE("ydist_fluxpos (Save)");
#define mccompcurname  ydist_fluxpos
#define mccompcurtype  PSDlin_monitor
#define mccompcurindex 13
#define nx mccydist_fluxpos_nx
#define PSDlin_N mccydist_fluxpos_PSDlin_N
#define PSDlin_p mccydist_fluxpos_PSDlin_p
#define PSDlin_p2 mccydist_fluxpos_PSDlin_p2
{   /* Declarations of ydist_fluxpos=PSDlin_monitor() SETTING parameters. */
char* filename = mccydist_fluxpos_filename;
MCNUM xmin = mccydist_fluxpos_xmin;
MCNUM xmax = mccydist_fluxpos_xmax;
MCNUM ymin = mccydist_fluxpos_ymin;
MCNUM ymax = mccydist_fluxpos_ymax;
MCNUM xwidth = mccydist_fluxpos_xwidth;
MCNUM yheight = mccydist_fluxpos_yheight;
MCNUM restore_neutron = mccydist_fluxpos_restore_neutron;
int nowritefile = mccydist_fluxpos_nowritefile;
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
#line 21341 "./PSI_DMC.c"
}   /* End of ydist_fluxpos=PSDlin_monitor() SETTING parameter declarations. */
#undef PSDlin_p2
#undef PSDlin_p
#undef PSDlin_N
#undef nx
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'PSD_fluxpos'. */
  SIG_MESSAGE("PSD_fluxpos (Save)");
#define mccompcurname  PSD_fluxpos
#define mccompcurtype  PSD_monitor
#define mccompcurindex 14
#define PSD_N mccPSD_fluxpos_PSD_N
#define PSD_p mccPSD_fluxpos_PSD_p
#define PSD_p2 mccPSD_fluxpos_PSD_p2
{   /* Declarations of PSD_fluxpos=PSD_monitor() SETTING parameters. */
int nx = mccPSD_fluxpos_nx;
int ny = mccPSD_fluxpos_ny;
char* filename = mccPSD_fluxpos_filename;
MCNUM xmin = mccPSD_fluxpos_xmin;
MCNUM xmax = mccPSD_fluxpos_xmax;
MCNUM ymin = mccPSD_fluxpos_ymin;
MCNUM ymax = mccPSD_fluxpos_ymax;
MCNUM xwidth = mccPSD_fluxpos_xwidth;
MCNUM yheight = mccPSD_fluxpos_yheight;
MCNUM restore_neutron = mccPSD_fluxpos_restore_neutron;
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
#line 21381 "./PSI_DMC.c"
}   /* End of PSD_fluxpos=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'xdist_flux_pos'. */
  SIG_MESSAGE("xdist_flux_pos (Save)");
#define mccompcurname  xdist_flux_pos
#define mccompcurtype  PSDlin_monitor
#define mccompcurindex 15
#define nx mccxdist_flux_pos_nx
#define PSDlin_N mccxdist_flux_pos_PSDlin_N
#define PSDlin_p mccxdist_flux_pos_PSDlin_p
#define PSDlin_p2 mccxdist_flux_pos_PSDlin_p2
{   /* Declarations of xdist_flux_pos=PSDlin_monitor() SETTING parameters. */
char* filename = mccxdist_flux_pos_filename;
MCNUM xmin = mccxdist_flux_pos_xmin;
MCNUM xmax = mccxdist_flux_pos_xmax;
MCNUM ymin = mccxdist_flux_pos_ymin;
MCNUM ymax = mccxdist_flux_pos_ymax;
MCNUM xwidth = mccxdist_flux_pos_xwidth;
MCNUM yheight = mccxdist_flux_pos_yheight;
MCNUM restore_neutron = mccxdist_flux_pos_restore_neutron;
int nowritefile = mccxdist_flux_pos_nowritefile;
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
#line 21421 "./PSI_DMC.c"
}   /* End of xdist_flux_pos=PSDlin_monitor() SETTING parameter declarations. */
#undef PSDlin_p2
#undef PSDlin_p
#undef PSDlin_N
#undef nx
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'PSD_fluxposB'. */
  SIG_MESSAGE("PSD_fluxposB (Save)");
#define mccompcurname  PSD_fluxposB
#define mccompcurtype  PSD_monitor
#define mccompcurindex 16
#define PSD_N mccPSD_fluxposB_PSD_N
#define PSD_p mccPSD_fluxposB_PSD_p
#define PSD_p2 mccPSD_fluxposB_PSD_p2
{   /* Declarations of PSD_fluxposB=PSD_monitor() SETTING parameters. */
int nx = mccPSD_fluxposB_nx;
int ny = mccPSD_fluxposB_ny;
char* filename = mccPSD_fluxposB_filename;
MCNUM xmin = mccPSD_fluxposB_xmin;
MCNUM xmax = mccPSD_fluxposB_xmax;
MCNUM ymin = mccPSD_fluxposB_ymin;
MCNUM ymax = mccPSD_fluxposB_ymax;
MCNUM xwidth = mccPSD_fluxposB_xwidth;
MCNUM yheight = mccPSD_fluxposB_yheight;
MCNUM restore_neutron = mccPSD_fluxposB_restore_neutron;
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
#line 21461 "./PSI_DMC.c"
}   /* End of PSD_fluxposB=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'lambda_in'. */
  SIG_MESSAGE("lambda_in (Save)");
#define mccompcurname  lambda_in
#define mccompcurtype  L_monitor
#define mccompcurindex 19
#define nL mcclambda_in_nL
#define L_N mcclambda_in_L_N
#define L_p mcclambda_in_L_p
#define L_p2 mcclambda_in_L_p2
{   /* Declarations of lambda_in=L_monitor() SETTING parameters. */
char* filename = mcclambda_in_filename;
MCNUM xmin = mcclambda_in_xmin;
MCNUM xmax = mcclambda_in_xmax;
MCNUM ymin = mcclambda_in_ymin;
MCNUM ymax = mcclambda_in_ymax;
MCNUM xwidth = mcclambda_in_xwidth;
MCNUM yheight = mcclambda_in_yheight;
MCNUM Lmin = mcclambda_in_Lmin;
MCNUM Lmax = mcclambda_in_Lmax;
MCNUM restore_neutron = mcclambda_in_restore_neutron;
int nowritefile = mcclambda_in_nowritefile;
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
#line 21503 "./PSI_DMC.c"
}   /* End of lambda_in=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'PSD_sample'. */
  SIG_MESSAGE("PSD_sample (Save)");
#define mccompcurname  PSD_sample
#define mccompcurtype  PSD_monitor
#define mccompcurindex 27
#define PSD_N mccPSD_sample_PSD_N
#define PSD_p mccPSD_sample_PSD_p
#define PSD_p2 mccPSD_sample_PSD_p2
{   /* Declarations of PSD_sample=PSD_monitor() SETTING parameters. */
int nx = mccPSD_sample_nx;
int ny = mccPSD_sample_ny;
char* filename = mccPSD_sample_filename;
MCNUM xmin = mccPSD_sample_xmin;
MCNUM xmax = mccPSD_sample_xmax;
MCNUM ymin = mccPSD_sample_ymin;
MCNUM ymax = mccPSD_sample_ymax;
MCNUM xwidth = mccPSD_sample_xwidth;
MCNUM yheight = mccPSD_sample_yheight;
MCNUM restore_neutron = mccPSD_sample_restore_neutron;
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
#line 21543 "./PSI_DMC.c"
}   /* End of PSD_sample=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'lambda_sample'. */
  SIG_MESSAGE("lambda_sample (Save)");
#define mccompcurname  lambda_sample
#define mccompcurtype  L_monitor
#define mccompcurindex 28
#define nL mcclambda_sample_nL
#define L_N mcclambda_sample_L_N
#define L_p mcclambda_sample_L_p
#define L_p2 mcclambda_sample_L_p2
{   /* Declarations of lambda_sample=L_monitor() SETTING parameters. */
char* filename = mcclambda_sample_filename;
MCNUM xmin = mcclambda_sample_xmin;
MCNUM xmax = mcclambda_sample_xmax;
MCNUM ymin = mcclambda_sample_ymin;
MCNUM ymax = mcclambda_sample_ymax;
MCNUM xwidth = mcclambda_sample_xwidth;
MCNUM yheight = mcclambda_sample_yheight;
MCNUM Lmin = mcclambda_sample_Lmin;
MCNUM Lmax = mcclambda_sample_Lmax;
MCNUM restore_neutron = mcclambda_sample_restore_neutron;
int nowritefile = mcclambda_sample_nowritefile;
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
#line 21585 "./PSI_DMC.c"
}   /* End of lambda_sample=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'Detector'. */
  SIG_MESSAGE("Detector (Save)");
#define mccompcurname  Detector
#define mccompcurtype  Monitor_nD
#define mccompcurindex 32
#define user1 mccDetector_user1
#define user2 mccDetector_user2
#define user3 mccDetector_user3
#define DEFS mccDetector_DEFS
#define Vars mccDetector_Vars
#define detector mccDetector_detector
#define offdata mccDetector_offdata
{   /* Declarations of Detector=Monitor_nD() SETTING parameters. */
MCNUM xwidth = mccDetector_xwidth;
MCNUM yheight = mccDetector_yheight;
MCNUM zdepth = mccDetector_zdepth;
MCNUM xmin = mccDetector_xmin;
MCNUM xmax = mccDetector_xmax;
MCNUM ymin = mccDetector_ymin;
MCNUM ymax = mccDetector_ymax;
MCNUM zmin = mccDetector_zmin;
MCNUM zmax = mccDetector_zmax;
MCNUM bins = mccDetector_bins;
MCNUM min = mccDetector_min;
MCNUM max = mccDetector_max;
MCNUM restore_neutron = mccDetector_restore_neutron;
MCNUM radius = mccDetector_radius;
char* options = mccDetector_options;
char* filename = mccDetector_filename;
char* geometry = mccDetector_geometry;
char* username1 = mccDetector_username1;
char* username2 = mccDetector_username2;
char* username3 = mccDetector_username3;
int nowritefile = mccDetector_nowritefile;
#line 483 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/Monitor_nD.comp"
{
  /* save results, but do not free pointers */
  if (!nowritefile) {
    detector = Monitor_nD_Save(&DEFS, &Vars);
  }
}
#line 21636 "./PSI_DMC.c"
}   /* End of Detector=Monitor_nD() SETTING parameter declarations. */
#undef offdata
#undef detector
#undef Vars
#undef DEFS
#undef user3
#undef user2
#undef user1
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  if (!handle) mcsiminfo_close(); 
} /* end save */
void mcfinally(void) {
  /* User component FINALLY code. */
  mcsiminfo_init(NULL);
  mcsave(mcsiminfo_file); /* save data when simulation ends */

  /* User FINALLY code for component 'source_arm'. */
  SIG_MESSAGE("source_arm (Finally)");
#define mccompcurname  source_arm
#define mccompcurtype  Progress_bar
#define mccompcurindex 1
#define IntermediateCnts mccsource_arm_IntermediateCnts
#define StartTime mccsource_arm_StartTime
#define EndTime mccsource_arm_EndTime
#define CurrentTime mccsource_arm_CurrentTime
{   /* Declarations of source_arm=Progress_bar() SETTING parameters. */
char* profile = mccsource_arm_profile;
MCNUM percent = mccsource_arm_percent;
MCNUM flag_save = mccsource_arm_flag_save;
MCNUM minutes = mccsource_arm_minutes;
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
#line 21683 "./PSI_DMC.c"
}   /* End of source_arm=Progress_bar() SETTING parameter declarations. */
#undef CurrentTime
#undef EndTime
#undef StartTime
#undef IntermediateCnts
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[1]) fprintf(stderr, "Warning: No neutron could reach Component[1] source_arm\n");
    if (mcAbsorbProp[1]) fprintf(stderr, "Warning: %g events were removed in Component[1] source_arm=Progress_bar()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[1]);
    if (!mcNCounter[2]) fprintf(stderr, "Warning: No neutron could reach Component[2] source\n");
    if (mcAbsorbProp[2]) fprintf(stderr, "Warning: %g events were removed in Component[2] source=Source_Maxwell_3()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[2]);
  /* User FINALLY code for component 'PSDbefore_guides'. */
  SIG_MESSAGE("PSDbefore_guides (Finally)");
#define mccompcurname  PSDbefore_guides
#define mccompcurtype  PSD_monitor
#define mccompcurindex 3
#define PSD_N mccPSDbefore_guides_PSD_N
#define PSD_p mccPSDbefore_guides_PSD_p
#define PSD_p2 mccPSDbefore_guides_PSD_p2
{   /* Declarations of PSDbefore_guides=PSD_monitor() SETTING parameters. */
int nx = mccPSDbefore_guides_nx;
int ny = mccPSDbefore_guides_ny;
char* filename = mccPSDbefore_guides_filename;
MCNUM xmin = mccPSDbefore_guides_xmin;
MCNUM xmax = mccPSDbefore_guides_xmax;
MCNUM ymin = mccPSDbefore_guides_ymin;
MCNUM ymax = mccPSDbefore_guides_ymax;
MCNUM xwidth = mccPSDbefore_guides_xwidth;
MCNUM yheight = mccPSDbefore_guides_yheight;
MCNUM restore_neutron = mccPSDbefore_guides_restore_neutron;
#line 122 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
{
  destroy_darr2d(PSD_N);
  destroy_darr2d(PSD_p);
  destroy_darr2d(PSD_p2);
}
#line 21722 "./PSI_DMC.c"
}   /* End of PSDbefore_guides=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[3]) fprintf(stderr, "Warning: No neutron could reach Component[3] PSDbefore_guides\n");
    if (mcAbsorbProp[3]) fprintf(stderr, "Warning: %g events were removed in Component[3] PSDbefore_guides=PSD_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[3]);
    if (!mcNCounter[4]) fprintf(stderr, "Warning: No neutron could reach Component[4] l_mon_source\n");
    if (mcAbsorbProp[4]) fprintf(stderr, "Warning: %g events were removed in Component[4] l_mon_source=L_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[4]);
    if (!mcNCounter[5]) fprintf(stderr, "Warning: No neutron could reach Component[5] guide1\n");
    if (mcAbsorbProp[5]) fprintf(stderr, "Warning: %g events were removed in Component[5] guide1=Guide()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[5]);
  /* User FINALLY code for component 'PSDbefore_curve'. */
  SIG_MESSAGE("PSDbefore_curve (Finally)");
#define mccompcurname  PSDbefore_curve
#define mccompcurtype  PSD_monitor
#define mccompcurindex 6
#define PSD_N mccPSDbefore_curve_PSD_N
#define PSD_p mccPSDbefore_curve_PSD_p
#define PSD_p2 mccPSDbefore_curve_PSD_p2
{   /* Declarations of PSDbefore_curve=PSD_monitor() SETTING parameters. */
int nx = mccPSDbefore_curve_nx;
int ny = mccPSDbefore_curve_ny;
char* filename = mccPSDbefore_curve_filename;
MCNUM xmin = mccPSDbefore_curve_xmin;
MCNUM xmax = mccPSDbefore_curve_xmax;
MCNUM ymin = mccPSDbefore_curve_ymin;
MCNUM ymax = mccPSDbefore_curve_ymax;
MCNUM xwidth = mccPSDbefore_curve_xwidth;
MCNUM yheight = mccPSDbefore_curve_yheight;
MCNUM restore_neutron = mccPSDbefore_curve_restore_neutron;
#line 122 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
{
  destroy_darr2d(PSD_N);
  destroy_darr2d(PSD_p);
  destroy_darr2d(PSD_p2);
}
#line 21762 "./PSI_DMC.c"
}   /* End of PSDbefore_curve=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[6]) fprintf(stderr, "Warning: No neutron could reach Component[6] PSDbefore_curve\n");
    if (mcAbsorbProp[6]) fprintf(stderr, "Warning: %g events were removed in Component[6] PSDbefore_curve=PSD_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[6]);
    if (!mcNCounter[7]) fprintf(stderr, "Warning: No neutron could reach Component[7] guide2\n");
    if (mcAbsorbProp[7]) fprintf(stderr, "Warning: %g events were removed in Component[7] guide2=Bender()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[7]);
  /* User FINALLY code for component 'PSDafter_curve'. */
  SIG_MESSAGE("PSDafter_curve (Finally)");
#define mccompcurname  PSDafter_curve
#define mccompcurtype  PSD_monitor
#define mccompcurindex 8
#define PSD_N mccPSDafter_curve_PSD_N
#define PSD_p mccPSDafter_curve_PSD_p
#define PSD_p2 mccPSDafter_curve_PSD_p2
{   /* Declarations of PSDafter_curve=PSD_monitor() SETTING parameters. */
int nx = mccPSDafter_curve_nx;
int ny = mccPSDafter_curve_ny;
char* filename = mccPSDafter_curve_filename;
MCNUM xmin = mccPSDafter_curve_xmin;
MCNUM xmax = mccPSDafter_curve_xmax;
MCNUM ymin = mccPSDafter_curve_ymin;
MCNUM ymax = mccPSDafter_curve_ymax;
MCNUM xwidth = mccPSDafter_curve_xwidth;
MCNUM yheight = mccPSDafter_curve_yheight;
MCNUM restore_neutron = mccPSDafter_curve_restore_neutron;
#line 122 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
{
  destroy_darr2d(PSD_N);
  destroy_darr2d(PSD_p);
  destroy_darr2d(PSD_p2);
}
#line 21800 "./PSI_DMC.c"
}   /* End of PSDafter_curve=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[8]) fprintf(stderr, "Warning: No neutron could reach Component[8] PSDafter_curve\n");
    if (mcAbsorbProp[8]) fprintf(stderr, "Warning: %g events were removed in Component[8] PSDafter_curve=PSD_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[8]);
    if (!mcNCounter[9]) fprintf(stderr, "Warning: No neutron could reach Component[9] bunker\n");
    if (mcAbsorbProp[9]) fprintf(stderr, "Warning: %g events were removed in Component[9] bunker=Guide()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[9]);
    if (!mcNCounter[10]) fprintf(stderr, "Warning: No neutron could reach Component[10] guide3\n");
    if (mcAbsorbProp[10]) fprintf(stderr, "Warning: %g events were removed in Component[10] guide3=Guide()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[10]);
    if (!mcNCounter[11]) fprintf(stderr, "Warning: No neutron could reach Component[11] guide4\n");
    if (mcAbsorbProp[11]) fprintf(stderr, "Warning: %g events were removed in Component[11] guide4=Guide()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[11]);
    if (!mcNCounter[12]) fprintf(stderr, "Warning: No neutron could reach Component[12] window1\n");
    if (mcAbsorbProp[12]) fprintf(stderr, "Warning: %g events were removed in Component[12] window1=Al_window()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[12]);
    if (!mcNCounter[13]) fprintf(stderr, "Warning: No neutron could reach Component[13] ydist_fluxpos\n");
    if (mcAbsorbProp[13]) fprintf(stderr, "Warning: %g events were removed in Component[13] ydist_fluxpos=PSDlin_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[13]);
  /* User FINALLY code for component 'PSD_fluxpos'. */
  SIG_MESSAGE("PSD_fluxpos (Finally)");
#define mccompcurname  PSD_fluxpos
#define mccompcurtype  PSD_monitor
#define mccompcurindex 14
#define PSD_N mccPSD_fluxpos_PSD_N
#define PSD_p mccPSD_fluxpos_PSD_p
#define PSD_p2 mccPSD_fluxpos_PSD_p2
{   /* Declarations of PSD_fluxpos=PSD_monitor() SETTING parameters. */
int nx = mccPSD_fluxpos_nx;
int ny = mccPSD_fluxpos_ny;
char* filename = mccPSD_fluxpos_filename;
MCNUM xmin = mccPSD_fluxpos_xmin;
MCNUM xmax = mccPSD_fluxpos_xmax;
MCNUM ymin = mccPSD_fluxpos_ymin;
MCNUM ymax = mccPSD_fluxpos_ymax;
MCNUM xwidth = mccPSD_fluxpos_xwidth;
MCNUM yheight = mccPSD_fluxpos_yheight;
MCNUM restore_neutron = mccPSD_fluxpos_restore_neutron;
#line 122 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
{
  destroy_darr2d(PSD_N);
  destroy_darr2d(PSD_p);
  destroy_darr2d(PSD_p2);
}
#line 21846 "./PSI_DMC.c"
}   /* End of PSD_fluxpos=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[14]) fprintf(stderr, "Warning: No neutron could reach Component[14] PSD_fluxpos\n");
    if (mcAbsorbProp[14]) fprintf(stderr, "Warning: %g events were removed in Component[14] PSD_fluxpos=PSD_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[14]);
    if (!mcNCounter[15]) fprintf(stderr, "Warning: No neutron could reach Component[15] xdist_flux_pos\n");
    if (mcAbsorbProp[15]) fprintf(stderr, "Warning: %g events were removed in Component[15] xdist_flux_pos=PSDlin_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[15]);
  /* User FINALLY code for component 'PSD_fluxposB'. */
  SIG_MESSAGE("PSD_fluxposB (Finally)");
#define mccompcurname  PSD_fluxposB
#define mccompcurtype  PSD_monitor
#define mccompcurindex 16
#define PSD_N mccPSD_fluxposB_PSD_N
#define PSD_p mccPSD_fluxposB_PSD_p
#define PSD_p2 mccPSD_fluxposB_PSD_p2
{   /* Declarations of PSD_fluxposB=PSD_monitor() SETTING parameters. */
int nx = mccPSD_fluxposB_nx;
int ny = mccPSD_fluxposB_ny;
char* filename = mccPSD_fluxposB_filename;
MCNUM xmin = mccPSD_fluxposB_xmin;
MCNUM xmax = mccPSD_fluxposB_xmax;
MCNUM ymin = mccPSD_fluxposB_ymin;
MCNUM ymax = mccPSD_fluxposB_ymax;
MCNUM xwidth = mccPSD_fluxposB_xwidth;
MCNUM yheight = mccPSD_fluxposB_yheight;
MCNUM restore_neutron = mccPSD_fluxposB_restore_neutron;
#line 122 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
{
  destroy_darr2d(PSD_N);
  destroy_darr2d(PSD_p);
  destroy_darr2d(PSD_p2);
}
#line 21884 "./PSI_DMC.c"
}   /* End of PSD_fluxposB=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[16]) fprintf(stderr, "Warning: No neutron could reach Component[16] PSD_fluxposB\n");
    if (mcAbsorbProp[16]) fprintf(stderr, "Warning: %g events were removed in Component[16] PSD_fluxposB=PSD_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[16]);
    if (!mcNCounter[17]) fprintf(stderr, "Warning: No neutron could reach Component[17] window2\n");
    if (mcAbsorbProp[17]) fprintf(stderr, "Warning: %g events were removed in Component[17] window2=Al_window()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[17]);
    if (!mcNCounter[18]) fprintf(stderr, "Warning: No neutron could reach Component[18] in_slit\n");
    if (mcAbsorbProp[18]) fprintf(stderr, "Warning: %g events were removed in Component[18] in_slit=Slit()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[18]);
    if (!mcNCounter[19]) fprintf(stderr, "Warning: No neutron could reach Component[19] lambda_in\n");
    if (mcAbsorbProp[19]) fprintf(stderr, "Warning: %g events were removed in Component[19] lambda_in=L_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[19]);
    if (!mcNCounter[20]) fprintf(stderr, "Warning: No neutron could reach Component[20] sma\n");
    if (mcAbsorbProp[20]) fprintf(stderr, "Warning: %g events were removed in Component[20] sma=Arm()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[20]);
    if (!mcNCounter[21]) fprintf(stderr, "Warning: No neutron could reach Component[21] foc_mono\n");
    if (mcNCounter[21] < 1000*(10)) fprintf(stderr, 
"Warning: Number of events %g reaching SPLIT position Component[21] foc_mono=Monochromator_2foc()\n"
"         is probably too low. Increase Ncount.\n", mcNCounter[21]);

    if (mcAbsorbProp[21]) fprintf(stderr, "Warning: %g events were removed in Component[21] foc_mono=Monochromator_2foc()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[21]);
    if (!mcNCounter[22]) fprintf(stderr, "Warning: No neutron could reach Component[22] msa\n");
    if (mcAbsorbProp[22]) fprintf(stderr, "Warning: %g events were removed in Component[22] msa=Arm()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[22]);
    if (!mcNCounter[23]) fprintf(stderr, "Warning: No neutron could reach Component[23] out1_slit\n");
    if (mcAbsorbProp[23]) fprintf(stderr, "Warning: %g events were removed in Component[23] out1_slit=Slit()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[23]);
    if (!mcNCounter[24]) fprintf(stderr, "Warning: No neutron could reach Component[24] Amoin_slit\n");
    if (mcAbsorbProp[24]) fprintf(stderr, "Warning: %g events were removed in Component[24] Amoin_slit=Slit()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[24]);
    if (!mcNCounter[25]) fprintf(stderr, "Warning: No neutron could reach Component[25] Bmoin_slit\n");
    if (mcAbsorbProp[25]) fprintf(stderr, "Warning: %g events were removed in Component[25] Bmoin_slit=Slit()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[25]);
    if (!mcNCounter[26]) fprintf(stderr, "Warning: No neutron could reach Component[26] out2_slit\n");
    if (mcAbsorbProp[26]) fprintf(stderr, "Warning: %g events were removed in Component[26] out2_slit=Slit()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[26]);
  /* User FINALLY code for component 'PSD_sample'. */
  SIG_MESSAGE("PSD_sample (Finally)");
#define mccompcurname  PSD_sample
#define mccompcurtype  PSD_monitor
#define mccompcurindex 27
#define PSD_N mccPSD_sample_PSD_N
#define PSD_p mccPSD_sample_PSD_p
#define PSD_p2 mccPSD_sample_PSD_p2
{   /* Declarations of PSD_sample=PSD_monitor() SETTING parameters. */
int nx = mccPSD_sample_nx;
int ny = mccPSD_sample_ny;
char* filename = mccPSD_sample_filename;
MCNUM xmin = mccPSD_sample_xmin;
MCNUM xmax = mccPSD_sample_xmax;
MCNUM ymin = mccPSD_sample_ymin;
MCNUM ymax = mccPSD_sample_ymax;
MCNUM xwidth = mccPSD_sample_xwidth;
MCNUM yheight = mccPSD_sample_yheight;
MCNUM restore_neutron = mccPSD_sample_restore_neutron;
#line 122 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
{
  destroy_darr2d(PSD_N);
  destroy_darr2d(PSD_p);
  destroy_darr2d(PSD_p2);
}
#line 21941 "./PSI_DMC.c"
}   /* End of PSD_sample=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[27]) fprintf(stderr, "Warning: No neutron could reach Component[27] PSD_sample\n");
    if (mcAbsorbProp[27]) fprintf(stderr, "Warning: %g events were removed in Component[27] PSD_sample=PSD_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[27]);
    if (!mcNCounter[28]) fprintf(stderr, "Warning: No neutron could reach Component[28] lambda_sample\n");
    if (mcAbsorbProp[28]) fprintf(stderr, "Warning: %g events were removed in Component[28] lambda_sample=L_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[28]);
    if (!mcNCounter[29]) fprintf(stderr, "Warning: No neutron could reach Component[29] sa_arm\n");
    if (mcAbsorbProp[29]) fprintf(stderr, "Warning: %g events were removed in Component[29] sa_arm=Arm()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[29]);
  /* User FINALLY code for component 'sample'. */
  SIG_MESSAGE("sample (Finally)");
#define mccompcurname  sample
#define mccompcurtype  PowderN
#define mccompcurindex 30
#define format mccsample_format
#define line_info mccsample_line_info
#define columns mccsample_columns
#define offdata mccsample_offdata
{   /* Declarations of sample=PowderN() SETTING parameters. */
char* reflections = mccsample_reflections;
char* geometry = mccsample_geometry;
MCNUM radius = mccsample_radius;
MCNUM yheight = mccsample_yheight;
MCNUM xwidth = mccsample_xwidth;
MCNUM zdepth = mccsample_zdepth;
MCNUM thickness = mccsample_thickness;
MCNUM pack = mccsample_pack;
MCNUM Vc = mccsample_Vc;
MCNUM sigma_abs = mccsample_sigma_abs;
MCNUM sigma_inc = mccsample_sigma_inc;
MCNUM delta_d_d = mccsample_delta_d_d;
MCNUM p_inc = mccsample_p_inc;
MCNUM p_transmit = mccsample_p_transmit;
MCNUM DW = mccsample_DW;
MCNUM nb_atoms = mccsample_nb_atoms;
MCNUM d_omega = mccsample_d_omega;
MCNUM d_phi = mccsample_d_phi;
MCNUM tth_sign = mccsample_tth_sign;
MCNUM p_interact = mccsample_p_interact;
MCNUM concentric = mccsample_concentric;
MCNUM density = mccsample_density;
MCNUM weight = mccsample_weight;
MCNUM barns = mccsample_barns;
MCNUM Strain = mccsample_Strain;
MCNUM focus_flip = mccsample_focus_flip;
int target_index = mccsample_target_index;
#line 1086 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../samples/PowderN.comp"
{
  free(line_info.list);
  free(line_info.q_v);
  free(line_info.w_v);
  free(line_info.my_s_v2);
  MPI_MASTER(
  if (line_info.flag_warning)
    printf("PowderN: %s: Error messages were repeated %i times with absorbed neutrons.\n",
      NAME_CURRENT_COMP, line_info.flag_warning);

  /* in case this instance is used in a SPLIT, we can recommend the
     optimal iteration value */
  if (line_info.nb_refl_count) {
    double split_iterations = (double)line_info.nb_reuses/line_info.nb_refl_count + 1;
    double split_optimal    = (double)line_info.nb_refl/line_info.nb_refl_count;
    if (split_optimal > split_iterations + 5)
      printf("PowderN: %s: Info: you may highly improve the computation efficiency by using\n"
        "    SPLIT %i COMPONENT %s=PowderN(...)\n"
        "  in the instrument description %s.\n",
        NAME_CURRENT_COMP, (int)split_optimal, NAME_CURRENT_COMP, mcinstrument_source);
  }
  );

}
#line 22018 "./PSI_DMC.c"
}   /* End of sample=PowderN() SETTING parameter declarations. */
#undef offdata
#undef columns
#undef line_info
#undef format
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[30]) fprintf(stderr, "Warning: No neutron could reach Component[30] sample\n");
    if (mcNCounter[30] < 1000*(10)) fprintf(stderr, 
"Warning: Number of events %g reaching SPLIT position Component[30] sample=PowderN()\n"
"         is probably too low. Increase Ncount.\n", mcNCounter[30]);

    if (mcAbsorbProp[30]) fprintf(stderr, "Warning: %g events were removed in Component[30] sample=PowderN()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[30]);
    if (!mcNCounter[31]) fprintf(stderr, "Warning: No neutron could reach Component[31] STOP\n");
    if (mcAbsorbProp[31]) fprintf(stderr, "Warning: %g events were removed in Component[31] STOP=Beamstop()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[31]);
  /* User FINALLY code for component 'Detector'. */
  SIG_MESSAGE("Detector (Finally)");
#define mccompcurname  Detector
#define mccompcurtype  Monitor_nD
#define mccompcurindex 32
#define user1 mccDetector_user1
#define user2 mccDetector_user2
#define user3 mccDetector_user3
#define DEFS mccDetector_DEFS
#define Vars mccDetector_Vars
#define detector mccDetector_detector
#define offdata mccDetector_offdata
{   /* Declarations of Detector=Monitor_nD() SETTING parameters. */
MCNUM xwidth = mccDetector_xwidth;
MCNUM yheight = mccDetector_yheight;
MCNUM zdepth = mccDetector_zdepth;
MCNUM xmin = mccDetector_xmin;
MCNUM xmax = mccDetector_xmax;
MCNUM ymin = mccDetector_ymin;
MCNUM ymax = mccDetector_ymax;
MCNUM zmin = mccDetector_zmin;
MCNUM zmax = mccDetector_zmax;
MCNUM bins = mccDetector_bins;
MCNUM min = mccDetector_min;
MCNUM max = mccDetector_max;
MCNUM restore_neutron = mccDetector_restore_neutron;
MCNUM radius = mccDetector_radius;
char* options = mccDetector_options;
char* filename = mccDetector_filename;
char* geometry = mccDetector_geometry;
char* username1 = mccDetector_username1;
char* username2 = mccDetector_username2;
char* username3 = mccDetector_username3;
int nowritefile = mccDetector_nowritefile;
#line 491 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/Monitor_nD.comp"
{
  /* free pointers */
  Monitor_nD_Finally(&DEFS, &Vars);
}
#line 22072 "./PSI_DMC.c"
}   /* End of Detector=Monitor_nD() SETTING parameter declarations. */
#undef offdata
#undef detector
#undef Vars
#undef DEFS
#undef user3
#undef user2
#undef user1
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[32]) fprintf(stderr, "Warning: No neutron could reach Component[32] Detector\n");
    if (mcAbsorbProp[32]) fprintf(stderr, "Warning: %g events were removed in Component[32] Detector=Monitor_nD()\n"
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

  /* MCDISPLAY code for component 'source_arm'. */
  SIG_MESSAGE("source_arm (McDisplay)");
  printf("MCDISPLAY: component %s\n", "source_arm");
#define mccompcurname  source_arm
#define mccompcurtype  Progress_bar
#define mccompcurindex 1
#define IntermediateCnts mccsource_arm_IntermediateCnts
#define StartTime mccsource_arm_StartTime
#define EndTime mccsource_arm_EndTime
#define CurrentTime mccsource_arm_CurrentTime
{   /* Declarations of source_arm=Progress_bar() SETTING parameters. */
char* profile = mccsource_arm_profile;
MCNUM percent = mccsource_arm_percent;
MCNUM flag_save = mccsource_arm_flag_save;
MCNUM minutes = mccsource_arm_minutes;
#line 147 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../misc/Progress_bar.comp"
{
  
}
#line 22121 "./PSI_DMC.c"
}   /* End of source_arm=Progress_bar() SETTING parameter declarations. */
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
#define mccompcurtype  Source_Maxwell_3
#define mccompcurindex 2
#define M mccsource_M
#define l_range mccsource_l_range
#define w_mult mccsource_w_mult
#define w_source mccsource_w_source
#define h_source mccsource_h_source
{   /* Declarations of source=Source_Maxwell_3() SETTING parameters. */
MCNUM size = mccsource_size;
MCNUM yheight = mccsource_yheight;
MCNUM xwidth = mccsource_xwidth;
MCNUM Lmin = mccsource_Lmin;
MCNUM Lmax = mccsource_Lmax;
MCNUM dist = mccsource_dist;
MCNUM focus_xw = mccsource_focus_xw;
MCNUM focus_yh = mccsource_focus_yh;
MCNUM T1 = mccsource_T1;
MCNUM T2 = mccsource_T2;
MCNUM T3 = mccsource_T3;
MCNUM I1 = mccsource_I1;
MCNUM I2 = mccsource_I2;
MCNUM I3 = mccsource_I3;
int target_index = mccsource_target_index;
MCNUM lambda0 = mccsource_lambda0;
MCNUM dlambda = mccsource_dlambda;
#line 155 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../sources/Source_Maxwell_3.comp"
{
  
  multiline(5, -(double)focus_xw/2.0, -(double)focus_yh/2.0, 0.0,
                (double)focus_xw/2.0, -(double)focus_yh/2.0, 0.0,
                (double)focus_xw/2.0,  (double)focus_yh/2.0, 0.0,
               -(double)focus_xw/2.0,  (double)focus_yh/2.0, 0.0,
               -(double)focus_xw/2.0, -(double)focus_yh/2.0, 0.0);
  if (dist) {
    dashed_line(0,0,0, -focus_xw/2,-focus_yh/2,dist, 4);
    dashed_line(0,0,0,  focus_xw/2,-focus_yh/2,dist, 4);
    dashed_line(0,0,0,  focus_xw/2, focus_yh/2,dist, 4);
    dashed_line(0,0,0, -focus_xw/2, focus_yh/2,dist, 4);
  }
}
#line 22175 "./PSI_DMC.c"
}   /* End of source=Source_Maxwell_3() SETTING parameter declarations. */
#undef h_source
#undef w_source
#undef w_mult
#undef l_range
#undef M
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'PSDbefore_guides'. */
  SIG_MESSAGE("PSDbefore_guides (McDisplay)");
  printf("MCDISPLAY: component %s\n", "PSDbefore_guides");
#define mccompcurname  PSDbefore_guides
#define mccompcurtype  PSD_monitor
#define mccompcurindex 3
#define PSD_N mccPSDbefore_guides_PSD_N
#define PSD_p mccPSDbefore_guides_PSD_p
#define PSD_p2 mccPSDbefore_guides_PSD_p2
{   /* Declarations of PSDbefore_guides=PSD_monitor() SETTING parameters. */
int nx = mccPSDbefore_guides_nx;
int ny = mccPSDbefore_guides_ny;
char* filename = mccPSDbefore_guides_filename;
MCNUM xmin = mccPSDbefore_guides_xmin;
MCNUM xmax = mccPSDbefore_guides_xmax;
MCNUM ymin = mccPSDbefore_guides_ymin;
MCNUM ymax = mccPSDbefore_guides_ymax;
MCNUM xwidth = mccPSDbefore_guides_xwidth;
MCNUM yheight = mccPSDbefore_guides_yheight;
MCNUM restore_neutron = mccPSDbefore_guides_restore_neutron;
#line 129 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
{
  multiline(5,
    (double)xmin, (double)ymin, 0.0,
    (double)xmax, (double)ymin, 0.0,
    (double)xmax, (double)ymax, 0.0,
    (double)xmin, (double)ymax, 0.0,
    (double)xmin, (double)ymin, 0.0);
}
#line 22215 "./PSI_DMC.c"
}   /* End of PSDbefore_guides=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'l_mon_source'. */
  SIG_MESSAGE("l_mon_source (McDisplay)");
  printf("MCDISPLAY: component %s\n", "l_mon_source");
#define mccompcurname  l_mon_source
#define mccompcurtype  L_monitor
#define mccompcurindex 4
#define nL mccl_mon_source_nL
#define L_N mccl_mon_source_L_N
#define L_p mccl_mon_source_L_p
#define L_p2 mccl_mon_source_L_p2
{   /* Declarations of l_mon_source=L_monitor() SETTING parameters. */
char* filename = mccl_mon_source_filename;
MCNUM xmin = mccl_mon_source_xmin;
MCNUM xmax = mccl_mon_source_xmax;
MCNUM ymin = mccl_mon_source_ymin;
MCNUM ymax = mccl_mon_source_ymax;
MCNUM xwidth = mccl_mon_source_xwidth;
MCNUM yheight = mccl_mon_source_yheight;
MCNUM Lmin = mccl_mon_source_Lmin;
MCNUM Lmax = mccl_mon_source_Lmax;
MCNUM restore_neutron = mccl_mon_source_restore_neutron;
int nowritefile = mccl_mon_source_nowritefile;
#line 120 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/L_monitor.comp"
{
  
  multiline(5, (double)xmin, (double)ymin, 0.0,
               (double)xmax, (double)ymin, 0.0,
               (double)xmax, (double)ymax, 0.0,
               (double)xmin, (double)ymax, 0.0,
               (double)xmin, (double)ymin, 0.0);
}
#line 22255 "./PSI_DMC.c"
}   /* End of l_mon_source=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'guide1'. */
  SIG_MESSAGE("guide1 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "guide1");
#define mccompcurname  guide1
#define mccompcurtype  Guide
#define mccompcurindex 5
#define pTable mccguide1_pTable
{   /* Declarations of guide1=Guide() SETTING parameters. */
char* reflect = mccguide1_reflect;
MCNUM w1 = mccguide1_w1;
MCNUM h1 = mccguide1_h1;
MCNUM w2 = mccguide1_w2;
MCNUM h2 = mccguide1_h2;
MCNUM l = mccguide1_l;
MCNUM R0 = mccguide1_R0;
MCNUM Qc = mccguide1_Qc;
MCNUM alpha = mccguide1_alpha;
MCNUM m = mccguide1_m;
MCNUM W = mccguide1_W;
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
#line 22304 "./PSI_DMC.c"
}   /* End of guide1=Guide() SETTING parameter declarations. */
#undef pTable
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'PSDbefore_curve'. */
  SIG_MESSAGE("PSDbefore_curve (McDisplay)");
  printf("MCDISPLAY: component %s\n", "PSDbefore_curve");
#define mccompcurname  PSDbefore_curve
#define mccompcurtype  PSD_monitor
#define mccompcurindex 6
#define PSD_N mccPSDbefore_curve_PSD_N
#define PSD_p mccPSDbefore_curve_PSD_p
#define PSD_p2 mccPSDbefore_curve_PSD_p2
{   /* Declarations of PSDbefore_curve=PSD_monitor() SETTING parameters. */
int nx = mccPSDbefore_curve_nx;
int ny = mccPSDbefore_curve_ny;
char* filename = mccPSDbefore_curve_filename;
MCNUM xmin = mccPSDbefore_curve_xmin;
MCNUM xmax = mccPSDbefore_curve_xmax;
MCNUM ymin = mccPSDbefore_curve_ymin;
MCNUM ymax = mccPSDbefore_curve_ymax;
MCNUM xwidth = mccPSDbefore_curve_xwidth;
MCNUM yheight = mccPSDbefore_curve_yheight;
MCNUM restore_neutron = mccPSDbefore_curve_restore_neutron;
#line 129 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
{
  multiline(5,
    (double)xmin, (double)ymin, 0.0,
    (double)xmax, (double)ymin, 0.0,
    (double)xmax, (double)ymax, 0.0,
    (double)xmin, (double)ymax, 0.0,
    (double)xmin, (double)ymin, 0.0);
}
#line 22340 "./PSI_DMC.c"
}   /* End of PSDbefore_curve=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'guide2'. */
  SIG_MESSAGE("guide2 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "guide2");
#define mccompcurname  guide2
#define mccompcurtype  Bender
#define mccompcurindex 7
#define bk mccguide2_bk
#define mWin mccguide2_mWin
{   /* Declarations of guide2=Bender() SETTING parameters. */
MCNUM w = mccguide2_w;
MCNUM h = mccguide2_h;
MCNUM r = mccguide2_r;
MCNUM Win = mccguide2_Win;
MCNUM k = mccguide2_k;
MCNUM d = mccguide2_d;
MCNUM l = mccguide2_l;
MCNUM R0a = mccguide2_R0a;
MCNUM Qca = mccguide2_Qca;
MCNUM alphaa = mccguide2_alphaa;
MCNUM ma = mccguide2_ma;
MCNUM Wa = mccguide2_Wa;
MCNUM R0i = mccguide2_R0i;
MCNUM Qci = mccguide2_Qci;
MCNUM alphai = mccguide2_alphai;
MCNUM mi = mccguide2_mi;
MCNUM Wi = mccguide2_Wi;
MCNUM R0s = mccguide2_R0s;
MCNUM Qcs = mccguide2_Qcs;
MCNUM alphas = mccguide2_alphas;
MCNUM ms = mccguide2_ms;
MCNUM Ws = mccguide2_Ws;
#line 269 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Bender.comp"
{
  int i;
  double w1c, w2c, h1, h2, L, w1, w2;

  w1c = (w + d)/(double)k;
  w2c = w1c; h1 = h; h2 = h;
  L = r*mWin; w1 = w; w2 = w;

  
  for(i = 0; i < k; i++)
  {
    multiline(5,
              i*w1c - w1/2.0, -h1/2.0, 0.0,
              i*w2c - w2/2.0, -h2/2.0, (double)L,
              i*w2c - w2/2.0,  h2/2.0, (double)L,
              i*w1c - w1/2.0,  h1/2.0, 0.0,
              i*w1c - w1/2.0, -h1/2.0, 0.0);
    multiline(5,
              (i+1)*w1c - d - w1/2.0, -h1/2.0, 0.0,
              (i+1)*w2c - d - w2/2.0, -h2/2.0, (double)L,
              (i+1)*w2c - d - w2/2.0,  h2/2.0, (double)L,
              (i+1)*w1c - d - w1/2.0,  h1/2.0, 0.0,
              (i+1)*w1c - d - w1/2.0, -h1/2.0, 0.0);
  }
  line(-w1/2.0, -h1/2.0, 0.0, w1/2.0, -h1/2.0, 0.0);
  line(-w2/2.0, -h2/2.0, (double)L, w2/2.0, -h2/2.0, (double)L);
}
#line 22408 "./PSI_DMC.c"
}   /* End of guide2=Bender() SETTING parameter declarations. */
#undef mWin
#undef bk
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'PSDafter_curve'. */
  SIG_MESSAGE("PSDafter_curve (McDisplay)");
  printf("MCDISPLAY: component %s\n", "PSDafter_curve");
#define mccompcurname  PSDafter_curve
#define mccompcurtype  PSD_monitor
#define mccompcurindex 8
#define PSD_N mccPSDafter_curve_PSD_N
#define PSD_p mccPSDafter_curve_PSD_p
#define PSD_p2 mccPSDafter_curve_PSD_p2
{   /* Declarations of PSDafter_curve=PSD_monitor() SETTING parameters. */
int nx = mccPSDafter_curve_nx;
int ny = mccPSDafter_curve_ny;
char* filename = mccPSDafter_curve_filename;
MCNUM xmin = mccPSDafter_curve_xmin;
MCNUM xmax = mccPSDafter_curve_xmax;
MCNUM ymin = mccPSDafter_curve_ymin;
MCNUM ymax = mccPSDafter_curve_ymax;
MCNUM xwidth = mccPSDafter_curve_xwidth;
MCNUM yheight = mccPSDafter_curve_yheight;
MCNUM restore_neutron = mccPSDafter_curve_restore_neutron;
#line 129 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
{
  multiline(5,
    (double)xmin, (double)ymin, 0.0,
    (double)xmax, (double)ymin, 0.0,
    (double)xmax, (double)ymax, 0.0,
    (double)xmin, (double)ymax, 0.0,
    (double)xmin, (double)ymin, 0.0);
}
#line 22445 "./PSI_DMC.c"
}   /* End of PSDafter_curve=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'bunker'. */
  SIG_MESSAGE("bunker (McDisplay)");
  printf("MCDISPLAY: component %s\n", "bunker");
#define mccompcurname  bunker
#define mccompcurtype  Guide
#define mccompcurindex 9
#define pTable mccbunker_pTable
{   /* Declarations of bunker=Guide() SETTING parameters. */
char* reflect = mccbunker_reflect;
MCNUM w1 = mccbunker_w1;
MCNUM h1 = mccbunker_h1;
MCNUM w2 = mccbunker_w2;
MCNUM h2 = mccbunker_h2;
MCNUM l = mccbunker_l;
MCNUM R0 = mccbunker_R0;
MCNUM Qc = mccbunker_Qc;
MCNUM alpha = mccbunker_alpha;
MCNUM m = mccbunker_m;
MCNUM W = mccbunker_W;
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
#line 22493 "./PSI_DMC.c"
}   /* End of bunker=Guide() SETTING parameter declarations. */
#undef pTable
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'guide3'. */
  SIG_MESSAGE("guide3 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "guide3");
#define mccompcurname  guide3
#define mccompcurtype  Guide
#define mccompcurindex 10
#define pTable mccguide3_pTable
{   /* Declarations of guide3=Guide() SETTING parameters. */
char* reflect = mccguide3_reflect;
MCNUM w1 = mccguide3_w1;
MCNUM h1 = mccguide3_h1;
MCNUM w2 = mccguide3_w2;
MCNUM h2 = mccguide3_h2;
MCNUM l = mccguide3_l;
MCNUM R0 = mccguide3_R0;
MCNUM Qc = mccguide3_Qc;
MCNUM alpha = mccguide3_alpha;
MCNUM m = mccguide3_m;
MCNUM W = mccguide3_W;
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
#line 22539 "./PSI_DMC.c"
}   /* End of guide3=Guide() SETTING parameter declarations. */
#undef pTable
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'guide4'. */
  SIG_MESSAGE("guide4 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "guide4");
#define mccompcurname  guide4
#define mccompcurtype  Guide
#define mccompcurindex 11
#define pTable mccguide4_pTable
{   /* Declarations of guide4=Guide() SETTING parameters. */
char* reflect = mccguide4_reflect;
MCNUM w1 = mccguide4_w1;
MCNUM h1 = mccguide4_h1;
MCNUM w2 = mccguide4_w2;
MCNUM h2 = mccguide4_h2;
MCNUM l = mccguide4_l;
MCNUM R0 = mccguide4_R0;
MCNUM Qc = mccguide4_Qc;
MCNUM alpha = mccguide4_alpha;
MCNUM m = mccguide4_m;
MCNUM W = mccguide4_W;
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
#line 22585 "./PSI_DMC.c"
}   /* End of guide4=Guide() SETTING parameter declarations. */
#undef pTable
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'window1'. */
  SIG_MESSAGE("window1 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "window1");
#define mccompcurname  window1
#define mccompcurtype  Al_window
#define mccompcurindex 12
{   /* Declarations of window1=Al_window() SETTING parameters. */
MCNUM thickness = mccwindow1_thickness;
#line 95 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../contrib/Al_window.comp"
{
  /* A bit ugly; hard-coded dimensions. */
  
  line(0,0,0,0.2,0,0);
  line(0,0,0,0,0.2,0);
  line(0,0,0,0,0,0.2);
}
#line 22608 "./PSI_DMC.c"
}   /* End of window1=Al_window() SETTING parameter declarations. */
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'ydist_fluxpos'. */
  SIG_MESSAGE("ydist_fluxpos (McDisplay)");
  printf("MCDISPLAY: component %s\n", "ydist_fluxpos");
#define mccompcurname  ydist_fluxpos
#define mccompcurtype  PSDlin_monitor
#define mccompcurindex 13
#define nx mccydist_fluxpos_nx
#define PSDlin_N mccydist_fluxpos_PSDlin_N
#define PSDlin_p mccydist_fluxpos_PSDlin_p
#define PSDlin_p2 mccydist_fluxpos_PSDlin_p2
{   /* Declarations of ydist_fluxpos=PSDlin_monitor() SETTING parameters. */
char* filename = mccydist_fluxpos_filename;
MCNUM xmin = mccydist_fluxpos_xmin;
MCNUM xmax = mccydist_fluxpos_xmax;
MCNUM ymin = mccydist_fluxpos_ymin;
MCNUM ymax = mccydist_fluxpos_ymax;
MCNUM xwidth = mccydist_fluxpos_xwidth;
MCNUM yheight = mccydist_fluxpos_yheight;
MCNUM restore_neutron = mccydist_fluxpos_restore_neutron;
int nowritefile = mccydist_fluxpos_nowritefile;
#line 116 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/PSDlin_monitor.comp"
{
  
  multiline(5, (double)xmin, (double)ymin, 0.0,
               (double)xmax, (double)ymin, 0.0,
               (double)xmax, (double)ymax, 0.0,
               (double)xmin, (double)ymax, 0.0,
               (double)xmin, (double)ymin, 0.0);
}
#line 22643 "./PSI_DMC.c"
}   /* End of ydist_fluxpos=PSDlin_monitor() SETTING parameter declarations. */
#undef PSDlin_p2
#undef PSDlin_p
#undef PSDlin_N
#undef nx
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'PSD_fluxpos'. */
  SIG_MESSAGE("PSD_fluxpos (McDisplay)");
  printf("MCDISPLAY: component %s\n", "PSD_fluxpos");
#define mccompcurname  PSD_fluxpos
#define mccompcurtype  PSD_monitor
#define mccompcurindex 14
#define PSD_N mccPSD_fluxpos_PSD_N
#define PSD_p mccPSD_fluxpos_PSD_p
#define PSD_p2 mccPSD_fluxpos_PSD_p2
{   /* Declarations of PSD_fluxpos=PSD_monitor() SETTING parameters. */
int nx = mccPSD_fluxpos_nx;
int ny = mccPSD_fluxpos_ny;
char* filename = mccPSD_fluxpos_filename;
MCNUM xmin = mccPSD_fluxpos_xmin;
MCNUM xmax = mccPSD_fluxpos_xmax;
MCNUM ymin = mccPSD_fluxpos_ymin;
MCNUM ymax = mccPSD_fluxpos_ymax;
MCNUM xwidth = mccPSD_fluxpos_xwidth;
MCNUM yheight = mccPSD_fluxpos_yheight;
MCNUM restore_neutron = mccPSD_fluxpos_restore_neutron;
#line 129 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
{
  multiline(5,
    (double)xmin, (double)ymin, 0.0,
    (double)xmax, (double)ymin, 0.0,
    (double)xmax, (double)ymax, 0.0,
    (double)xmin, (double)ymax, 0.0,
    (double)xmin, (double)ymin, 0.0);
}
#line 22682 "./PSI_DMC.c"
}   /* End of PSD_fluxpos=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'xdist_flux_pos'. */
  SIG_MESSAGE("xdist_flux_pos (McDisplay)");
  printf("MCDISPLAY: component %s\n", "xdist_flux_pos");
#define mccompcurname  xdist_flux_pos
#define mccompcurtype  PSDlin_monitor
#define mccompcurindex 15
#define nx mccxdist_flux_pos_nx
#define PSDlin_N mccxdist_flux_pos_PSDlin_N
#define PSDlin_p mccxdist_flux_pos_PSDlin_p
#define PSDlin_p2 mccxdist_flux_pos_PSDlin_p2
{   /* Declarations of xdist_flux_pos=PSDlin_monitor() SETTING parameters. */
char* filename = mccxdist_flux_pos_filename;
MCNUM xmin = mccxdist_flux_pos_xmin;
MCNUM xmax = mccxdist_flux_pos_xmax;
MCNUM ymin = mccxdist_flux_pos_ymin;
MCNUM ymax = mccxdist_flux_pos_ymax;
MCNUM xwidth = mccxdist_flux_pos_xwidth;
MCNUM yheight = mccxdist_flux_pos_yheight;
MCNUM restore_neutron = mccxdist_flux_pos_restore_neutron;
int nowritefile = mccxdist_flux_pos_nowritefile;
#line 116 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/PSDlin_monitor.comp"
{
  
  multiline(5, (double)xmin, (double)ymin, 0.0,
               (double)xmax, (double)ymin, 0.0,
               (double)xmax, (double)ymax, 0.0,
               (double)xmin, (double)ymax, 0.0,
               (double)xmin, (double)ymin, 0.0);
}
#line 22720 "./PSI_DMC.c"
}   /* End of xdist_flux_pos=PSDlin_monitor() SETTING parameter declarations. */
#undef PSDlin_p2
#undef PSDlin_p
#undef PSDlin_N
#undef nx
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'PSD_fluxposB'. */
  SIG_MESSAGE("PSD_fluxposB (McDisplay)");
  printf("MCDISPLAY: component %s\n", "PSD_fluxposB");
#define mccompcurname  PSD_fluxposB
#define mccompcurtype  PSD_monitor
#define mccompcurindex 16
#define PSD_N mccPSD_fluxposB_PSD_N
#define PSD_p mccPSD_fluxposB_PSD_p
#define PSD_p2 mccPSD_fluxposB_PSD_p2
{   /* Declarations of PSD_fluxposB=PSD_monitor() SETTING parameters. */
int nx = mccPSD_fluxposB_nx;
int ny = mccPSD_fluxposB_ny;
char* filename = mccPSD_fluxposB_filename;
MCNUM xmin = mccPSD_fluxposB_xmin;
MCNUM xmax = mccPSD_fluxposB_xmax;
MCNUM ymin = mccPSD_fluxposB_ymin;
MCNUM ymax = mccPSD_fluxposB_ymax;
MCNUM xwidth = mccPSD_fluxposB_xwidth;
MCNUM yheight = mccPSD_fluxposB_yheight;
MCNUM restore_neutron = mccPSD_fluxposB_restore_neutron;
#line 129 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
{
  multiline(5,
    (double)xmin, (double)ymin, 0.0,
    (double)xmax, (double)ymin, 0.0,
    (double)xmax, (double)ymax, 0.0,
    (double)xmin, (double)ymax, 0.0,
    (double)xmin, (double)ymin, 0.0);
}
#line 22759 "./PSI_DMC.c"
}   /* End of PSD_fluxposB=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'window2'. */
  SIG_MESSAGE("window2 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "window2");
#define mccompcurname  window2
#define mccompcurtype  Al_window
#define mccompcurindex 17
{   /* Declarations of window2=Al_window() SETTING parameters. */
MCNUM thickness = mccwindow2_thickness;
#line 95 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../contrib/Al_window.comp"
{
  /* A bit ugly; hard-coded dimensions. */
  
  line(0,0,0,0.2,0,0);
  line(0,0,0,0,0.2,0);
  line(0,0,0,0,0,0.2);
}
#line 22784 "./PSI_DMC.c"
}   /* End of window2=Al_window() SETTING parameter declarations. */
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'in_slit'. */
  SIG_MESSAGE("in_slit (McDisplay)");
  printf("MCDISPLAY: component %s\n", "in_slit");
#define mccompcurname  in_slit
#define mccompcurtype  Slit
#define mccompcurindex 18
{   /* Declarations of in_slit=Slit() SETTING parameters. */
MCNUM xmin = mccin_slit_xmin;
MCNUM xmax = mccin_slit_xmax;
MCNUM ymin = mccin_slit_ymin;
MCNUM ymax = mccin_slit_ymax;
MCNUM radius = mccin_slit_radius;
MCNUM xwidth = mccin_slit_xwidth;
MCNUM yheight = mccin_slit_yheight;
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
#line 22827 "./PSI_DMC.c"
}   /* End of in_slit=Slit() SETTING parameter declarations. */
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'lambda_in'. */
  SIG_MESSAGE("lambda_in (McDisplay)");
  printf("MCDISPLAY: component %s\n", "lambda_in");
#define mccompcurname  lambda_in
#define mccompcurtype  L_monitor
#define mccompcurindex 19
#define nL mcclambda_in_nL
#define L_N mcclambda_in_L_N
#define L_p mcclambda_in_L_p
#define L_p2 mcclambda_in_L_p2
{   /* Declarations of lambda_in=L_monitor() SETTING parameters. */
char* filename = mcclambda_in_filename;
MCNUM xmin = mcclambda_in_xmin;
MCNUM xmax = mcclambda_in_xmax;
MCNUM ymin = mcclambda_in_ymin;
MCNUM ymax = mcclambda_in_ymax;
MCNUM xwidth = mcclambda_in_xwidth;
MCNUM yheight = mcclambda_in_yheight;
MCNUM Lmin = mcclambda_in_Lmin;
MCNUM Lmax = mcclambda_in_Lmax;
MCNUM restore_neutron = mcclambda_in_restore_neutron;
int nowritefile = mcclambda_in_nowritefile;
#line 120 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/L_monitor.comp"
{
  
  multiline(5, (double)xmin, (double)ymin, 0.0,
               (double)xmax, (double)ymin, 0.0,
               (double)xmax, (double)ymax, 0.0,
               (double)xmin, (double)ymax, 0.0,
               (double)xmin, (double)ymin, 0.0);
}
#line 22864 "./PSI_DMC.c"
}   /* End of lambda_in=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'sma'. */
  SIG_MESSAGE("sma (McDisplay)");
  printf("MCDISPLAY: component %s\n", "sma");
#define mccompcurname  sma
#define mccompcurtype  Arm
#define mccompcurindex 20
#line 40 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Arm.comp"
{
  /* A bit ugly; hard-coded dimensions. */
  
  line(0,0,0,0.2,0,0);
  line(0,0,0,0,0.2,0);
  line(0,0,0,0,0,0.2);
}
#line 22888 "./PSI_DMC.c"
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'foc_mono'. */
  SIG_MESSAGE("foc_mono (McDisplay)");
  printf("MCDISPLAY: component %s\n", "foc_mono");
#define mccompcurname  foc_mono
#define mccompcurtype  Monochromator_2foc
#define mccompcurindex 21
#define mos_y mccfoc_mono_mos_y
#define mos_z mccfoc_mono_mos_z
#define mono_Q mccfoc_mono_mono_Q
#define SlabWidth mccfoc_mono_SlabWidth
#define SlabHeight mccfoc_mono_SlabHeight
#define rTable mccfoc_mono_rTable
{   /* Declarations of foc_mono=Monochromator_2foc() SETTING parameters. */
char* reflect = mccfoc_mono_reflect;
MCNUM zwidth = mccfoc_mono_zwidth;
MCNUM yheight = mccfoc_mono_yheight;
MCNUM gap = mccfoc_mono_gap;
MCNUM NH = mccfoc_mono_NH;
MCNUM NV = mccfoc_mono_NV;
MCNUM mosaich = mccfoc_mono_mosaich;
MCNUM mosaicv = mccfoc_mono_mosaicv;
MCNUM r0 = mccfoc_mono_r0;
MCNUM Q = mccfoc_mono_Q;
MCNUM RV = mccfoc_mono_RV;
MCNUM RH = mccfoc_mono_RH;
MCNUM DM = mccfoc_mono_DM;
MCNUM mosaic = mccfoc_mono_mosaic;
MCNUM width = mccfoc_mono_width;
MCNUM height = mccfoc_mono_height;
MCNUM verbose = mccfoc_mono_verbose;
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
#line 22958 "./PSI_DMC.c"
}   /* End of foc_mono=Monochromator_2foc() SETTING parameter declarations. */
#undef rTable
#undef SlabHeight
#undef SlabWidth
#undef mono_Q
#undef mos_z
#undef mos_y
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'msa'. */
  SIG_MESSAGE("msa (McDisplay)");
  printf("MCDISPLAY: component %s\n", "msa");
#define mccompcurname  msa
#define mccompcurtype  Arm
#define mccompcurindex 22
#line 40 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Arm.comp"
{
  /* A bit ugly; hard-coded dimensions. */
  
  line(0,0,0,0.2,0,0);
  line(0,0,0,0,0.2,0);
  line(0,0,0,0,0,0.2);
}
#line 22984 "./PSI_DMC.c"
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'out1_slit'. */
  SIG_MESSAGE("out1_slit (McDisplay)");
  printf("MCDISPLAY: component %s\n", "out1_slit");
#define mccompcurname  out1_slit
#define mccompcurtype  Slit
#define mccompcurindex 23
{   /* Declarations of out1_slit=Slit() SETTING parameters. */
MCNUM xmin = mccout1_slit_xmin;
MCNUM xmax = mccout1_slit_xmax;
MCNUM ymin = mccout1_slit_ymin;
MCNUM ymax = mccout1_slit_ymax;
MCNUM radius = mccout1_slit_radius;
MCNUM xwidth = mccout1_slit_xwidth;
MCNUM yheight = mccout1_slit_yheight;
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
#line 23026 "./PSI_DMC.c"
}   /* End of out1_slit=Slit() SETTING parameter declarations. */
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'Amoin_slit'. */
  SIG_MESSAGE("Amoin_slit (McDisplay)");
  printf("MCDISPLAY: component %s\n", "Amoin_slit");
#define mccompcurname  Amoin_slit
#define mccompcurtype  Slit
#define mccompcurindex 24
{   /* Declarations of Amoin_slit=Slit() SETTING parameters. */
MCNUM xmin = mccAmoin_slit_xmin;
MCNUM xmax = mccAmoin_slit_xmax;
MCNUM ymin = mccAmoin_slit_ymin;
MCNUM ymax = mccAmoin_slit_ymax;
MCNUM radius = mccAmoin_slit_radius;
MCNUM xwidth = mccAmoin_slit_xwidth;
MCNUM yheight = mccAmoin_slit_yheight;
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
#line 23069 "./PSI_DMC.c"
}   /* End of Amoin_slit=Slit() SETTING parameter declarations. */
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'Bmoin_slit'. */
  SIG_MESSAGE("Bmoin_slit (McDisplay)");
  printf("MCDISPLAY: component %s\n", "Bmoin_slit");
#define mccompcurname  Bmoin_slit
#define mccompcurtype  Slit
#define mccompcurindex 25
{   /* Declarations of Bmoin_slit=Slit() SETTING parameters. */
MCNUM xmin = mccBmoin_slit_xmin;
MCNUM xmax = mccBmoin_slit_xmax;
MCNUM ymin = mccBmoin_slit_ymin;
MCNUM ymax = mccBmoin_slit_ymax;
MCNUM radius = mccBmoin_slit_radius;
MCNUM xwidth = mccBmoin_slit_xwidth;
MCNUM yheight = mccBmoin_slit_yheight;
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
#line 23112 "./PSI_DMC.c"
}   /* End of Bmoin_slit=Slit() SETTING parameter declarations. */
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'out2_slit'. */
  SIG_MESSAGE("out2_slit (McDisplay)");
  printf("MCDISPLAY: component %s\n", "out2_slit");
#define mccompcurname  out2_slit
#define mccompcurtype  Slit
#define mccompcurindex 26
{   /* Declarations of out2_slit=Slit() SETTING parameters. */
MCNUM xmin = mccout2_slit_xmin;
MCNUM xmax = mccout2_slit_xmax;
MCNUM ymin = mccout2_slit_ymin;
MCNUM ymax = mccout2_slit_ymax;
MCNUM radius = mccout2_slit_radius;
MCNUM xwidth = mccout2_slit_xwidth;
MCNUM yheight = mccout2_slit_yheight;
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
#line 23155 "./PSI_DMC.c"
}   /* End of out2_slit=Slit() SETTING parameter declarations. */
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'PSD_sample'. */
  SIG_MESSAGE("PSD_sample (McDisplay)");
  printf("MCDISPLAY: component %s\n", "PSD_sample");
#define mccompcurname  PSD_sample
#define mccompcurtype  PSD_monitor
#define mccompcurindex 27
#define PSD_N mccPSD_sample_PSD_N
#define PSD_p mccPSD_sample_PSD_p
#define PSD_p2 mccPSD_sample_PSD_p2
{   /* Declarations of PSD_sample=PSD_monitor() SETTING parameters. */
int nx = mccPSD_sample_nx;
int ny = mccPSD_sample_ny;
char* filename = mccPSD_sample_filename;
MCNUM xmin = mccPSD_sample_xmin;
MCNUM xmax = mccPSD_sample_xmax;
MCNUM ymin = mccPSD_sample_ymin;
MCNUM ymax = mccPSD_sample_ymax;
MCNUM xwidth = mccPSD_sample_xwidth;
MCNUM yheight = mccPSD_sample_yheight;
MCNUM restore_neutron = mccPSD_sample_restore_neutron;
#line 129 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
{
  multiline(5,
    (double)xmin, (double)ymin, 0.0,
    (double)xmax, (double)ymin, 0.0,
    (double)xmax, (double)ymax, 0.0,
    (double)xmin, (double)ymax, 0.0,
    (double)xmin, (double)ymin, 0.0);
}
#line 23190 "./PSI_DMC.c"
}   /* End of PSD_sample=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'lambda_sample'. */
  SIG_MESSAGE("lambda_sample (McDisplay)");
  printf("MCDISPLAY: component %s\n", "lambda_sample");
#define mccompcurname  lambda_sample
#define mccompcurtype  L_monitor
#define mccompcurindex 28
#define nL mcclambda_sample_nL
#define L_N mcclambda_sample_L_N
#define L_p mcclambda_sample_L_p
#define L_p2 mcclambda_sample_L_p2
{   /* Declarations of lambda_sample=L_monitor() SETTING parameters. */
char* filename = mcclambda_sample_filename;
MCNUM xmin = mcclambda_sample_xmin;
MCNUM xmax = mcclambda_sample_xmax;
MCNUM ymin = mcclambda_sample_ymin;
MCNUM ymax = mcclambda_sample_ymax;
MCNUM xwidth = mcclambda_sample_xwidth;
MCNUM yheight = mcclambda_sample_yheight;
MCNUM Lmin = mcclambda_sample_Lmin;
MCNUM Lmax = mcclambda_sample_Lmax;
MCNUM restore_neutron = mcclambda_sample_restore_neutron;
int nowritefile = mcclambda_sample_nowritefile;
#line 120 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/L_monitor.comp"
{
  
  multiline(5, (double)xmin, (double)ymin, 0.0,
               (double)xmax, (double)ymin, 0.0,
               (double)xmax, (double)ymax, 0.0,
               (double)xmin, (double)ymax, 0.0,
               (double)xmin, (double)ymin, 0.0);
}
#line 23230 "./PSI_DMC.c"
}   /* End of lambda_sample=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'sa_arm'. */
  SIG_MESSAGE("sa_arm (McDisplay)");
  printf("MCDISPLAY: component %s\n", "sa_arm");
#define mccompcurname  sa_arm
#define mccompcurtype  Arm
#define mccompcurindex 29
#line 40 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Arm.comp"
{
  /* A bit ugly; hard-coded dimensions. */
  
  line(0,0,0,0.2,0,0);
  line(0,0,0,0,0.2,0);
  line(0,0,0,0,0,0.2);
}
#line 23254 "./PSI_DMC.c"
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'sample'. */
  SIG_MESSAGE("sample (McDisplay)");
  printf("MCDISPLAY: component %s\n", "sample");
#define mccompcurname  sample
#define mccompcurtype  PowderN
#define mccompcurindex 30
#define format mccsample_format
#define line_info mccsample_line_info
#define columns mccsample_columns
#define offdata mccsample_offdata
{   /* Declarations of sample=PowderN() SETTING parameters. */
char* reflections = mccsample_reflections;
char* geometry = mccsample_geometry;
MCNUM radius = mccsample_radius;
MCNUM yheight = mccsample_yheight;
MCNUM xwidth = mccsample_xwidth;
MCNUM zdepth = mccsample_zdepth;
MCNUM thickness = mccsample_thickness;
MCNUM pack = mccsample_pack;
MCNUM Vc = mccsample_Vc;
MCNUM sigma_abs = mccsample_sigma_abs;
MCNUM sigma_inc = mccsample_sigma_inc;
MCNUM delta_d_d = mccsample_delta_d_d;
MCNUM p_inc = mccsample_p_inc;
MCNUM p_transmit = mccsample_p_transmit;
MCNUM DW = mccsample_DW;
MCNUM nb_atoms = mccsample_nb_atoms;
MCNUM d_omega = mccsample_d_omega;
MCNUM d_phi = mccsample_d_phi;
MCNUM tth_sign = mccsample_tth_sign;
MCNUM p_interact = mccsample_p_interact;
MCNUM concentric = mccsample_concentric;
MCNUM density = mccsample_density;
MCNUM weight = mccsample_weight;
MCNUM barns = mccsample_barns;
MCNUM Strain = mccsample_Strain;
MCNUM focus_flip = mccsample_focus_flip;
int target_index = mccsample_target_index;
#line 1112 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../samples/PowderN.comp"
{
  if (line_info.V_0) {
    magnify("xyz");
    if (line_info.shape == 0) { /* cyl */
      circle("xz", 0,  yheight/2.0, 0, radius);
      circle("xz", 0, -yheight/2.0, 0, radius);
      line(-radius, -yheight/2.0, 0, -radius, +yheight/2.0, 0);
      line(+radius, -yheight/2.0, 0, +radius, +yheight/2.0, 0);
      line(0, -yheight/2.0, -radius, 0, +yheight/2.0, -radius);
      line(0, -yheight/2.0, +radius, 0, +yheight/2.0, +radius);
      if (thickness) {
        double radius_i=radius-thickness;
        circle("xz", 0,  yheight/2.0, 0, radius_i);
        circle("xz", 0, -yheight/2.0, 0, radius_i);
        line(-radius_i, -yheight/2.0, 0, -radius_i, +yheight/2.0, 0);
        line(+radius_i, -yheight/2.0, 0, +radius_i, +yheight/2.0, 0);
        line(0, -yheight/2.0, -radius_i, 0, +yheight/2.0, -radius_i);
        line(0, -yheight/2.0, +radius_i, 0, +yheight/2.0, +radius_i);
      }
    } else if (line_info.shape == 1) {  /* box */
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
      if (line_info.zdepth_i) {
        xmin = -0.5*line_info.xwidth_i;
        xmax =  0.5*line_info.xwidth_i;
        ymin = -0.5*line_info.yheight_i;
        ymax =  0.5*line_info.yheight_i;
        zmin = -0.5*line_info.zdepth_i;
        zmax =  0.5*line_info.zdepth_i;
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
    } if (line_info.shape == 2) { /* sphere */
      if (line_info.radius_i) {
        circle("xy",0,0,0,line_info.radius_i);
        circle("xz",0,0,0,line_info.radius_i);
        circle("yz",0,0,0,line_info.radius_i);
      }
      circle("xy",0,0,0,radius);
      circle("xz",0,0,0,radius);
      circle("yz",0,0,0,radius);
    } else if (line_info.shape == 3) {	/* OFF file */
      off_display(offdata);
    }
  }
}
#line 23374 "./PSI_DMC.c"
}   /* End of sample=PowderN() SETTING parameter declarations. */
#undef offdata
#undef columns
#undef line_info
#undef format
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'STOP'. */
  SIG_MESSAGE("STOP (McDisplay)");
  printf("MCDISPLAY: component %s\n", "STOP");
#define mccompcurname  STOP
#define mccompcurtype  Beamstop
#define mccompcurindex 31
{   /* Declarations of STOP=Beamstop() SETTING parameters. */
MCNUM xmin = mccSTOP_xmin;
MCNUM xmax = mccSTOP_xmax;
MCNUM ymin = mccSTOP_ymin;
MCNUM ymax = mccSTOP_ymax;
MCNUM xwidth = mccSTOP_xwidth;
MCNUM yheight = mccSTOP_yheight;
MCNUM radius = mccSTOP_radius;
#line 72 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Beamstop.comp"
{
  
  if (radius != 0)
    circle("xy", 0, 0, 0, radius);
  else
    multiline(5, (double)xmin, (double)ymin, 0.0,
               (double)xmax, (double)ymin, 0.0,
               (double)xmax, (double)ymax, 0.0,
               (double)xmin, (double)ymax, 0.0,
               (double)xmin, (double)ymin, 0.0);
}
#line 23410 "./PSI_DMC.c"
}   /* End of STOP=Beamstop() SETTING parameter declarations. */
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'Detector'. */
  SIG_MESSAGE("Detector (McDisplay)");
  printf("MCDISPLAY: component %s\n", "Detector");
#define mccompcurname  Detector
#define mccompcurtype  Monitor_nD
#define mccompcurindex 32
#define user1 mccDetector_user1
#define user2 mccDetector_user2
#define user3 mccDetector_user3
#define DEFS mccDetector_DEFS
#define Vars mccDetector_Vars
#define detector mccDetector_detector
#define offdata mccDetector_offdata
{   /* Declarations of Detector=Monitor_nD() SETTING parameters. */
MCNUM xwidth = mccDetector_xwidth;
MCNUM yheight = mccDetector_yheight;
MCNUM zdepth = mccDetector_zdepth;
MCNUM xmin = mccDetector_xmin;
MCNUM xmax = mccDetector_xmax;
MCNUM ymin = mccDetector_ymin;
MCNUM ymax = mccDetector_ymax;
MCNUM zmin = mccDetector_zmin;
MCNUM zmax = mccDetector_zmax;
MCNUM bins = mccDetector_bins;
MCNUM min = mccDetector_min;
MCNUM max = mccDetector_max;
MCNUM restore_neutron = mccDetector_restore_neutron;
MCNUM radius = mccDetector_radius;
char* options = mccDetector_options;
char* filename = mccDetector_filename;
char* geometry = mccDetector_geometry;
char* username1 = mccDetector_username1;
char* username2 = mccDetector_username2;
char* username3 = mccDetector_username3;
int nowritefile = mccDetector_nowritefile;
#line 497 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/Monitor_nD.comp"
{
  if (geometry && strlen(geometry) && strcmp(geometry,"0") && strcmp(geometry, "NULL"))
  {
    off_display(offdata);
  } else {
    Monitor_nD_McDisplay(&DEFS, &Vars);
  }
}
#line 23460 "./PSI_DMC.c"
}   /* End of Detector=Monitor_nD() SETTING parameter declarations. */
#undef offdata
#undef detector
#undef Vars
#undef DEFS
#undef user3
#undef user2
#undef user1
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
/* end of generated C code ./PSI_DMC.c */
