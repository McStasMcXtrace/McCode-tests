/* Automatically generated file. Do not edit. 
 * Format:     ANSI C source code
 * Creator:    McStas <http://www.mcstas.org>
 * Instrument: /zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr (Test_Fermi)
 * Date:       Wed Nov 20 00:49:34 2019
 * File:       ./Test_Fermi.c
 * Compile:    cc -o Test_Fermi.out ./Test_Fermi.c 
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

#line 712 "./Test_Fermi.c"

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

#line 945 "./Test_Fermi.c"

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

#line 4977 "./Test_Fermi.c"

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

#line 5337 "./Test_Fermi.c"
#ifdef MC_TRACE_ENABLED
int mctraceenabled = 1;
#else
int mctraceenabled = 0;
#endif
#define MCSTAS "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../"
int mcdefaultmain = 1;
char mcinstrument_name[] = "Test_Fermi";
char mcinstrument_source[] = "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr";
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

#line 6839 "./Test_Fermi.c"

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

#line 9774 "./Test_Fermi.c"

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
#line 10130 "./Test_Fermi.c"

/* Shared user declarations for all components 'Guide_channeled'. */
#line 76 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide_channeled.comp"

#line 10135 "./Test_Fermi.c"

/* Shared user declarations for all components 'FermiChopper'. */
#line 97 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/FermiChopper.comp"

#ifndef FermiChopper_TimeAccuracy
#define FermiChopper_TimeAccuracy 1e-9
#define FermiChopper_MAXITER      100

  /* Definition of internal variable structure: all counters */
  struct FermiChopper_struct {

    double omega, ph0, t0;  /* chopper rotation */
    double C_slit;  /* C_slit radius in [m] */
    double L_slit;
    double sum_t;
    double sum_v;
    double sum_N;
    double sum_N_pass;
    /* events */
    long absorb_alreadyinside;
    long absorb_topbottom;
    long absorb_cylentrance;
    long absorb_sideentrance;
    long absorb_notreachentrance;
    long absorb_packentrance;
    long absorb_slitcoating;
    long warn_notreachslitwall;
    long absorb_exitslitpack;
    long absorb_maxiterations;
    long absorb_wrongdirection;
    long absorb_nocontrol;
    long absorb_cylexit;
    long warn_notreachslitoutput;
    char compcurname[256];
  };

/*****************************************************************************
 * FC_zrot: returns Z' in rotating frame, from X,Z and t,omega,ph0
 ****************************************************************************/
double FC_zrot(double X, double Z, double T, struct FermiChopper_struct FCs){
  double omega =FCs.omega;
  double ph0   =FCs.ph0;

  return( Z*cos(omega*T+ph0)-X*sin(omega*T+ph0) );
}


/*****************************************************************************
 * FC_xrot: returns X' in rotating frame, from X,Z and omega,t,ph0
 *          additional coordinate shift in case of curved slits
 ****************************************************************************/
double FC_xrot(double X, double Z, double T, struct FermiChopper_struct FCs){
  double omega =FCs.omega;
  double ph0   =FCs.ph0;
  double C_slit=FCs.C_slit;
  double ret, tmp;

  ret = X*cos(omega*T+ph0)+Z*sin(omega*T+ph0);

  if (C_slit) {
    tmp  = fabs(FC_zrot(X, Z, T, FCs));
    if (tmp < FCs.L_slit/2) {
      tmp  = (FCs.L_slit/2 - tmp)*C_slit;
      ret += (1-sqrt(1-tmp*tmp))/C_slit;
    }
  }
  return( ret );
}

/*****************************************************************************
 * FC_xzrot_dt(x,z,vx,vz, t,dt, type='x' or 'z', FCs)
 *   returns X' or Z' in rotating frame, from X,Z and t,omega,ph0
 *              taking into account propagation with velocity during time dt
 ****************************************************************************/
double FC_xzrot_dt(double x, double z, double vx, double vz,
                   double t, double dt, char type, struct FermiChopper_struct FCs) {
  if (dt) /* with propagation */
    return( (type == 'x' ? FC_xrot(x+vx*dt, z+vz*dt, t+dt, FCs)
                         : FC_zrot(x+vx*dt, z+vz*dt, t+dt, FCs)) );
  else    /* without propagation */
    return( (type == 'x' ? FC_xrot(x,z,t,FCs)
                         : FC_zrot(x,z,t,FCs)) );
}

/*****************************************************************************
 * FC_xzbrent(x,z,vx,vz, t,dt, type='x' or 'z', d, FCs)
 *   solves X'=d and Z'=d with Brent algorithm in time interval [0, dt].
 *           Returns time within [0,dt], from NumRecip in C, chap 9, p360 (zbrent)
 *           ERRORS: return -1 should never occur
 *                          -2 if exceed MAX iteration
 *                          -3 no sign change in range
 ****************************************************************************/
double FC_xzbrent(double x, double z, double vx, double vz,
                  double t, double dt,
                  char type, double d, struct FermiChopper_struct FCs) {

  int iter;
  double a=0,b=dt;
  double c,dd,e,min1,min2;
  double tol=FermiChopper_TimeAccuracy;
  double EPS=FermiChopper_TimeAccuracy;
  double fa=FC_xzrot_dt(x,z,vx,vz, t,a, type, FCs) - d;
  double fb=FC_xzrot_dt(x,z,vx,vz, t,b, type, FCs) - d;
  double fc,p,q,r,s,tol1,xm;

  if (fb*fa > 0.0) return -3;
  fc=fb;
  for (iter=1;iter<=FermiChopper_MAXITER;iter++) {
    if (fb*fc > 0.0) {
      c=a;
      fc=fa;
      e=dd=b-a;
    }
    if (fabs(fc) < fabs(fb)) {
      a=b;
      b=c;
      c=a;
      fa=fb;
      fb=fc;
      fc=fa;
    }
    tol1=2.0*EPS*fabs(b)+0.5*tol;
    xm=0.5*(c-b);
    if (fabs(xm) <= tol1 || fb == 0.0) return b;
    if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
      s=fb/fa;
      if (a == c) {
        p=2.0*xm*s;
        q=1.0-s;
      } else {
        q=fa/fc;
        r=fb/fc;
        p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
        q=(q-1.0)*(r-1.0)*(s-1.0);
      }
      if (p > 0.0)  q = -q;
      p=fabs(p);
      min1=3.0*xm*q-fabs(tol1*q);
      min2=fabs(e*q);
      if (2.0*p < (min1 < min2 ? min1 : min2)) {
        e=dd;
        dd=p/q;
      } else {
        dd=xm;
        e=dd;
      }
    } else {
      dd=xm;
      e=dd;
    }
    a=b;
    fa=fb;
    if (fabs(dd) > tol1)
      b += dd;
    else
      b += (xm > 0.0 ? fabs(tol1) : -fabs(tol1));
    fb=FC_xzrot_dt(x,z,vx,vz, t,b, type, FCs) - d;
  }
  return -2;
} /* FC_xzbrent */

/*****************************************************************************
 * Wrappers to intersection algorithms 
 ****************************************************************************/

double FC_xintersect(double x, double z, double vx, double vz,
                   double t, double dt,
                   double d, struct FermiChopper_struct FCs) {
  return(FC_xzbrent(x, z, vx, vz, t, dt, 'x', d, FCs));
}

double FC_zintersect(double x, double z, double vx, double vz,
                   double t, double dt,
                   double d, struct FermiChopper_struct FCs) {
  return(FC_xzbrent(x, z, vx, vz, t, dt, 'z', d, FCs));
}

#endif

#line 10315 "./Test_Fermi.c"

/* Shared user declarations for all components 'FermiChopper_ILL'. */
#line 103 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../contrib/FermiChopper_ILL.comp"
#ifndef FCILL_TimeAccuracy
#define FCILL_TimeAccuracy 1e-8
#define FCILL_MAXITER      10
/* Definition of internal variable structure: all counters */
struct FermiChopper_ILL_struct {

/** other variables ********************************/
double omega, t0;  /* chopper rotation */
};


/**************** DECLARING FUNCTIONS ***************************************/

/*********** ORTHOGONAL TRANSFORMATION INTO ROTATING FRAME ******************/

/************ X - component ********************/
double xstrich(double X, double Z, double T, double omega, double t0){
  return( X*cos(omega*(T-t0))+Z*sin(omega*(T-t0)) );
}

/************ Z - component ********************/
double zstrich(double X, double Z, double T, double omega, double t0){
  return( Z*cos(omega*(T-t0))-X*sin(omega*(T-t0)) );
}

/*************************NUMERICAL METHODS *********************************/

/*************************** SECANT METHOD FOR... ***************************/

/****************************...X-component *********************************/
double xsecant(double x, double z, double vx, double vz,
               double t, double dt, double d, double omega, double t0){

  double dt1     = 1;
  double counter = 0;
  double t1      = 0;
  double t2      = dt;
  double xr1     = xstrich(x,z,t, omega, t0)-d;
  double xr2     = xstrich(x+vx*t2,z+vz*t2,t+t2, omega, t0)-d;
  double sign;

  while ((fabs(dt1) > FCILL_TimeAccuracy) && (counter < FCILL_MAXITER) && (xr2-xr1)){
    counter++;
    dt1 = (t2-t1)*xr2/(xr2-xr1);
    t2  = t1;
    xr1 = xr2;
    t1 += dt1;
    xr2 = xstrich(x+vx*t1,z+vz*t1,t+t1, omega, t0)-d;
  }

  if(counter >= FCILL_MAXITER) t1 = -2;

  return(t1);
}


/****************************...Z-component *********************************/
double zsecant(double x, double z, double vx, double vz,
               double t, double dt, double d, double omega, double t0) {

  double t1      = 0;
  double t2      = dt;
  double dt1     = 1;
  double counter = 0;
  double zr1     =  zstrich(x,z,t, omega, t0)-d;
  double zr2     =  zstrich(x+vx*t2,z+vz*t2,t+t2, omega, t0)-d;

  while ((fabs(dt1) > FCILL_TimeAccuracy) && (counter < FCILL_MAXITER) && (zr2-zr1)){
    counter++;
    dt1 = (t2-t1)*zr2/(zr2-zr1);
    t2  = t1;
    zr1 = zr2;
    t1 += dt1;
    zr2 = zstrich(x+vx*t1,z+vz*t1,t+t1, omega, t0)-d;
  }

  if(counter >= FCILL_MAXITER) t1=-1;

  return(t1);
}


/*************************** INTERPOLATION METHOD FOR... ********************/

/****************************...X-component *********************************/
double xinterpolation(double x, double z, double vx, double vz,
                      double t, double dt, double d, double omega, double t0){

  double sign;
  double xr3=1, t3=0, t1=0, t2=dt, dt1=dt;
  double counter = 0;
  double xr1      =  xstrich(x,z,t, omega, t0)-d;
  double xr2      =  xstrich(x+vx*dt,z+vz*dt,t+dt, omega, t0)-d;

  while ((fabs(xr3) > FCILL_TimeAccuracy)&&(counter < FCILL_MAXITER)){
    counter++;
    t3 = (t1+t2)*0.5;
    xr3 = xstrich(x+(vx*(t3)),z+(vz*(t3)),t+t3, omega, t0)-d;
    xr2 = xstrich(x+(vx*(t2)),z+(vz*(t2)),t+t2, omega, t0)-d;
    if(xr2*xr3<0) t1=t3;
    else          t2=t3;
  }

  if(counter >= FCILL_MAXITER) t3=-1;

  return(t3);
}


/****************************...Z-component *********************************/
double zinterpolation(double x, double z, double vx, double vz,
                      double t, double dt, double d, double omega, double t0){

  double counter = 0;
  double zr3=1,zr2=0,t3=0,t1=0,t2=dt;

  while ((fabs(zr3)>FCILL_TimeAccuracy)&&(counter<FCILL_MAXITER)) {
    counter++;
    t3 = (t1+t2)*0.5;
    zr3 = zstrich(x+(vx*(t3)),z+(vz*(t3)),t+t3, omega, t0)-d;
    zr2 = zstrich(x+(vx*(t2)),z+(vz*(t2)),t+t2, omega, t0)-d;
    if(zr2*zr3 < 0) t1=t3;
    else            t2=t3;
  }

  if(counter >= FCILL_MAXITER) t3=-1;

  return(t3);
}
#endif
#line 10449 "./Test_Fermi.c"

/* Shared user declarations for all components 'Vitess_ChopperFermi'. */
#line 130 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Vitess_ChopperFermi.comp"
/********************************************************************************************/
/*  VITESS module 'general.c'                                                               */
/*    Elementary functions for all VITESS modules                                           */
/*                                                                                          */
/* The free non-commercial use of these routines is granted providing due credit is given   */
/* to the authors:                                                                          */
/* Friedrich Streffer, Géza Zsigmond, Dietmar Wechsler,                                     */
/* Michael Fromme, Klaus Lieutenant, Sergey Manoshin                                        */ 
/*                                                                                          */
/* Change: K.L.  2002 JAN, reorganized routines                                             */
/* Change: G.Zs. 2002 JUL, new routines                                                     */
/* Change: K.L.  2003 JAN, new functions 'ReadLine', 'StrgLShift', and 'StrgCopy'           */
/* Change: K.L.  2003 FEB, definitions of 'idum' and 'LogFilePtr' from init to general      */
/* Change: K.L.  2003 MAR, new function 'StrgScanLF', additional parameter in 'ReadLine'    */

#ifdef VITESS
 #include "general.h"
#else
#ifndef GENERAL_H
#define GENERAL_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/******************************/
/** Definitions              **/
/******************************/

#ifdef _MSC_VER
#define M_PI            3.14159265358979323846  /* pi */ 
#define M_PI_2          1.57079632679489661923  /* pi/2 */
#endif

#define MN          1.6749284E-27
#define G           9.80665 
#define K           1.380662E-23
#define NA          6.022137E23
#define H           6.626076E-34
#define L_2_E       81805.048
#define E_C         1.6021773E-19

#define TRUE 		1
#define FALSE 		0

#define UP          1
#define DOWN        0

#define ON          1
#define OFF         0

#define GUIDEFLIGHT 1

#define MAX_COLLISIONS      100
#define MAX_CHOPPER_WINDOWS  10
#define LAMBDA_MIN            0.001
#define LAMBDA_MAX          100.0

#define BUFFER_SIZE       10000
#ifndef CHAR_BUF_LENGTH
#define CHAR_BUF_LENGTH    1024
#endif
#define CHAR_BUF_LARGE     5120
#define CHAR_BUF_SMALL      256

#ifdef VITESS
 #define ABSORB continue
#endif

typedef enum 
{	VT_CUBE   = 1,
	VT_CYL    = 2,
	VT_SPHERE = 3
}
SampleGeom;


typedef enum
{
	VT_SOURCE      =   1,
	VT_GUIDE       =  11,
	VT_BENDER      =  12,
	VT_COLLIMATOR  =  13,
	VT_SM_ENSEMBLE =  15,
	VT_SPACE       =  20,
	VT_WINDOW      =  21,
	VT_WND_MULT    =  22,
	VT_GRID        =  23,
	VT_CHOP_DISC   =  31,
	VT_CHOP_FERMI  =  32,
	VT_VEL_SELECT  =  41,
	VT_MONOC_ANALY =  45,
	VT_POL_HE3     =  51,
	VT_POL_SM      =  52,
	VT_POL_MIRROR  =  53,
	VT_FLIP_COIL   =  55,
	VT_FLIP_GRAD   =  56,
	VT_RES_DRABKIN =  59,
	VT_PREC_FIELD  =  60,
	VT_ROT_FIELD   =  61,
	VT_DETECTOR    =  71,
	VT_SMPL_EL_ISO =  81,
	VT_SMPL_INELAST=  83,
	VT_SMPL_SING_X =  84,
	VT_SMPL_POWDER =  85,
	VT_SMPL_S_Q    =  86,
	VT_SMPL_SANS   =  87,
	VT_SMPL_REFL   =  89,
	VT_MONITOR_1   = 101,
	VT_MONITOR_2   = 102,
	VT_MON_POL_1   = 103,
	VT_MON_POL_POS = 104,
	VT_EVAL_ELAST  = 111,
	VT_EVAL_INELAST= 112,
	VT_VISUAL      = 121,
	VT_FRAME       = 131,
	VT_WRITEOUT    = 141,
	VT_TOOL        = 999
}
VtModID;


typedef double VectorType[3];
typedef double DoublePair[2];


/******************************/
/** Structures               **/
/******************************/

typedef struct
{
	double X,Y,Z;
}
CartesianPoint;


typedef struct
{
	double	A, B, C, D;
}
Plane;


typedef struct
{
	double  A, B, C, D, E, F, W, P, Q, R;
}
SurfaceSecond;


typedef struct
{
	char           IDGrp[2];
	unsigned long  IDNo;
}
TotalID;


typedef struct 
{
	TotalID        ID;
	char           Debug;
	short          Color;
	double         Time;
	double         Wavelength;
	double         Probability;
	VectorType     Position;
	VectorType     Vector;
	VectorType     Spin;
}
Neutron;


typedef struct 
{
      double height, r;
} 
CylinderType;


typedef struct 
{
      double height, width, thickness;
} 
CubeType;


typedef struct 
{
      double r;
} 
BallType;


typedef union 
{
    CylinderType Cyl;
    CubeType     Cube;
    BallType     Ball;
} 
SampleGeomType;


typedef struct 
{
  SampleGeom Type;
  VectorType Position;
  VectorType Direction;
  SampleGeomType SG;
} 
SampleType;

typedef struct
{
	VtModID  eModule;
	double   dWPar;    /* width, ...             */
	double   dHPar;    /* height, end width, ... */
	double   dRPar;    /* radius, ...            */
	long     nNumber;  /* number of ....         */
	short    eType;    /* shape, mon. par., ...  */
	char*    pDescr;   /* material, ...          */
}
ModProp;



/******************************/
/** Prototypes               **/
/******************************/

double ENERGY_FROM_LAMBDA(double x);
double LAMBDA_FROM_ENERGY(double x);
double ENERGY_FROM_V   (double x);    
double V_FROM_LAMBDA   (double x);
double LAMBDA_FROM_V(double x);

double ran3       (long * i);
double MonteCarlo (double x, double y);
double vector3rand(double *, double *, double *);

double sq   (double Value);                        /* = Value*Value*/
double atan0(double a, double b);
void   Exchange(double* pValue1, double* pValue2);
double Min(double value1, double value2);
double Max(double value1, double value2);
long   mini(long value1, long value2);
long   maxi(long value1, long value2);

double SolidAngle   (const double dHorAngle, const double dVertAngle);

void   CopyVector   (const VectorType Src, VectorType Dest);
long   MAXV         (const VectorType Vector);
double LengthVector (const VectorType Vector);
double DistVector   (const VectorType Vec1, const VectorType Vec2);
double ScalarProduct(const VectorType Vec1, const VectorType Vec2);
double AngleVectors (VectorType v1, VectorType v2);
double Area(VectorType v1, VectorType v2);
short  NormVector      (VectorType Vector);
void   AddVector       (VectorType Value,  const VectorType Add);
void   SubVector       (VectorType Value,  const VectorType Sub);
void   MultiplyByScalar(VectorType Vector, const double Scalar);
void   RotVector       (double RotMatrix[3][3], VectorType Vector);
void   RotBackVector   (double RotMatrix[3][3], VectorType Vector);
void   FillRMatrixZY   (double RotMatrix[3][3], double roty, double rotz);

FILE * fileOpen(const char *name, char *mode);
void   Error(const char *text);
void   Warning(const char *text);
void   Abort();

int    ReadLine(FILE* pFile, char* pLine, int nStrLen);
void   ReadParString(FILE *fpt, char *stringvar);
double ReadParF(FILE *fpt);
int    ReadParI(FILE *fpt);
void   ReadParComment(FILE *fpt);

void   StrgCopy  (char* sCopy, const char* sOrigin, int nLen);
void   StrgLShift(char* sStr, int kWidth);
long   StrgScanLF(const char* sStr, double* pTable, const int nMax, const int nStart);

#endif


#endif
#include "ctype.h"


#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/(double)MBIG)

long  idum;              /* parameter for the random number generation  */
FILE* LogFilePtr;        /* pointer to the log file stream              */


/****************************************************************************************/
/*  Coversion between physical properties                                               */
/****************************************************************************************/

double ENERGY_FROM_LAMBDA(double x) 
{	
	return(81805.048 / x / x); /*[ueV]*/
}

double LAMBDA_FROM_ENERGY(double x)
{
	return(sqrt(81805.048 / x));  /*[ueV]*/
}

double ENERGY_FROM_V(double x) 
{
	return(0.5227032667573 * x * x); /*[ueV]*/
}

double V_FROM_LAMBDA(double x)
{
	return(395.60346 / x); /* cm/ms 395.6034613488 */
}
double LAMBDA_FROM_V(double x)
{
	return(395.60346 / x); /* cm/ms 395.6034613488 */
}


/****************************************************************************************/
/*  Random Functions                                                                    */
/****************************************************************************************/

/* 'ran3' computes next random variable in [0,1[         */
/* (C) Copr. 1986-92 Numerical Recipes Software #>,1'59. */

double ran3(long *idum)
{
	static int inext,inextp;
	static long ma[56];
	static int iff=0;
	long mj,mk;
	int i,ii,k;

	if (*idum < 0 || iff == 0) {
		iff=1;
		mj=MSEED-(*idum < 0 ? -*idum : *idum);
		mj %= MBIG;
		ma[55]=mj;
		mk=1;
		for (i=1;i<=54;i++) {
			ii=(21*i) % 55;
			ma[ii]=mk;
			mk=mj-mk;
			while (mk < MZ) mk += MBIG;
			mj=ma[ii];
		}
		for (k=1;k<=4;k++)
			for (i=1;i<=55;i++) {
				ma[i] -= ma[1+(i+30) % 55];
				 while (ma[i] < MZ) ma[i] += MBIG;
			}
		inext=0;
		inextp=31;
		*idum=1;
	}
	if (++inext == 56) inext=1;
	if (++inextp == 56) inextp=1;
	mj=ma[inext]-ma[inextp];
	while (mj < MZ) mj += MBIG;  /*modification */
	ma[inext]=mj;
	return FAC*(double)mj;
}
#undef MBIG
#undef MSEED
#undef MZ
#undef FAC

double MonteCarlo(double x, double y)
{
   return((y-x)*ran3(&idum) + x);
}

/* 'vector3rand' simulates randomly orientied vector with unit lenght */
/*  (Manoshin Sergey 28.02.01)                                        */

double vector3rand(double *vx, double *vy, double *vz)
{
  double len;
  len = 0.0;

  while ((len == 0.0)||(len > 1.0))
    {
      *vx = 2.0*ran3(&idum) - 1.0;
      *vy = 2.0*ran3(&idum) - 1.0;
      *vz = 2.0*ran3(&idum) - 1.0;
      len = (*vx)*(*vx) + (*vy)*(*vy) + (*vz)*(*vz);
    }

  /*     Normalising  */
  len = sqrt(len);

  *vx = *vx/len;
  *vy = *vy/len;
  *vz = *vz/len;

  return(len);
}


/****************************************************************************************/
/*  General Functions                                                                   */
/****************************************************************************************/

/* computes square of a real value */

double sq(double Value)
{
	return Value * Value ;
}



/* calculates atan2 in the range (0, 2*M_PI) */

double atan0(double a, double b)
{
	if (b > 0.)		return (double) atan(a / b) ;

	if (b == 0.)	return M_PI_2 ;

	else			return (double) atan(a / b) + M_PI ;
}



/* minimum and maximum of two double or long values */

long mini(long value1, long value2)
{
	if(value1 < value2) return value1 ;
	else                return value2 ;
}

long maxi(long value1, long value2)
{
	if(value1 > value2) return value1 ;
	else                return value2 ;
}

double Min(double value1, double value2)
{
	if(value1 < value2) return value1 ;
	else                return value2 ;
}

double Max(double value1, double value2)
{
	if(value1 > value2) return value1 ;
	else                return value2 ;
}



/* 'Exchange' of two values */

void Exchange(double* pValue1, double* pValue2)
{
	double dHelp;

	dHelp    = *pValue1 ;
	*pValue1 = *pValue2;
	*pValue2 = dHelp;

	return;
}



/****************************************************************************************/
/*  Vector Functions                                                                    */
/****************************************************************************************/

/* 'Copy' copies the contents of Vector 'Src' to vector 'Dest'  */
/*                                                    */
void CopyVector(const VectorType Src, VectorType Dest)
{
	int i;

	for(i=0;i<3;i++)
	{	Dest[i] = Src[i];
	}
}


/* 'MAXV' returns the number of the largest component of 'Vector': 0, 1 or 2  */
/*                                                                            */
long MAXV(const VectorType Vector)
{
	if( (fabs(Vector[0]) > fabs(Vector[1])) && (fabs(Vector[0]) > fabs(Vector[2])))
		return 0;
	else
		if(fabs(Vector[1]) > fabs(Vector[2]))
			return 1;
		else
			return 2;
}


/* 'LengthVector' returns the length of vector 'Vec'  */
/*                                                    */
double LengthVector(const VectorType Vec)
{
	return sqrt(ScalarProduct(Vec,Vec));
}


/* 'NormVector' changes the vector length to 1  */
/*                                              */
short NormVector(VectorType Vector)
{
	long   i;
	double dLen = LengthVector(Vector);	

	if (dLen==0.0)
	{	return FALSE;
	}
	else
	{	for(i=0;i<3;i++)
		{	Vector[i] /= dLen;
		}
		return TRUE;
	}
}


/* 'DistVector' calculates the distance between the points described by Vec1 and Vec2  */
/*                                                                                     */
double DistVector(const VectorType Vec1, const VectorType Vec2)
{
	VectorType Vhlp;

	CopyVector(Vec1, Vhlp) ;
	SubVector (Vhlp, Vec2);
	return LengthVector(Vhlp);
}


/* 'AddVector' adds 'Add' to 'Value' and returns 'Value'  */
/*                                                        */
void AddVector(VectorType Value, const VectorType Add)
{
	int i ;
	VectorType Result ;
	for(i=0;i<3;i++)
		Result[i]=Value[i]+ Add[i] ;
	CopyVector(Result, Value) ;
}


/* 'SubVector' Substracts 'Sub' to 'Value' and returns 'Value' */
/*                                                             */
void SubVector(VectorType Value, const VectorType Sub)
{
	int i ;
	VectorType Result ;
	for(i=0;i<3;i++)
	{
		Result[i]=Value[i]- Sub[i] ;
		Value[i] = Result[i] ;
	}
}


/* 'MultiplyByScalar' multiplies a vector by a scalar */
/*                                                    */
void MultiplyByScalar(VectorType Vector, const double Scalar)
{
	int i ;
	VectorType Result ;
	for(i=0;i<3;i++)
	{
		Result[i] = Scalar * Vector[i] ;
		Vector[i] = Result[i] ;
	}
}


/* 'ScalarProduct' calculates the scalar product of two vectors 'v1' and 'v2' */
/*                                                                            */
double ScalarProduct(const VectorType v1, const VectorType v2)
{
	int		j;
	double	result;

	result = 0.;
	for(j=0;j<3;j++) result += v1[j]*v2[j] ;

	return result ;
}

/* angle between two vectors in degs */

double AngleVectors(VectorType v1, VectorType v2)
{
double theta ;
 
	theta= ScalarProduct(v1, v2) / (double)sqrt(ScalarProduct(v1, v1)) / (double)sqrt(ScalarProduct(v2, v2)) ;
	return 180./M_PI * (double) acos(theta) ;
}

/* area of triangle from two vectors, G.Zs */


double Area(VectorType v1, VectorType v2)
{
return LengthVector(v1) * LengthVector(v2) * 
		
		fabs(sin(acos( ScalarProduct(v1, v2)/(LengthVector(v1) * LengthVector(v2)))) /2.);

}


/* 'RotVector' does essentially a Vector times matrix multiplication   */
/* in order to rotate the Vector. The rotation Matrix may be supplied  */
/* by e.g. 'RotMatrixX'                                                */
/* Author: F. Streffer.                                                */
void RotVector(double RotMatrix[3][3], VectorType Vector)
{
	VectorType TempVec;
	int        i;

	for(i=0;i<3;i++)
		TempVec[i]=ScalarProduct(RotMatrix[i],Vector);
	CopyVector(TempVec, Vector);
}

/* 'RotBackVector' rotates a vector, by multiplication of the Vector  */
/* with the invers of RotMatrix. E.g if RotMatrix is the same as in   */
/* 'RotVector' and is applied to the result of RotVector the original */
/* vector is restored.   (Remember det(RotMatrix)=1)                  */
/* Author: F. Streffer.                                               */
void RotBackVector(double RotMatrix[3][3], VectorType Vector)
{
	VectorType TempVec;
	int        i;

	for(i=0;i<3;i++)
		TempVec[i]=RotMatrix[0][i]*Vector[0]+RotMatrix[1][i]*Vector[1]+RotMatrix[2][i]*Vector[2];
	CopyVector(TempVec, Vector);
}


/* 'FillRMatrixZY' calculates a rotation matrix, which rotates a frame   */
/* at first about the z-axis by 'rotz' and then about the y-axis by 'roty' */
/*  Author: F. Streffer.                                                   */
/*  Change: G. Zs. 16 JUL 2002  rotation convention                        */
void FillRMatrixZY(double RotMatrix[3][3], double roty, double rotz)
{
  double sz, cz, sy, cy;
  long   i,j;

  sy= (double) sin(roty);
  cy= (double) cos(roty);
  sz= (double) sin(rotz);
  cz= (double) cos(rotz);

  /* now, fill the matrix */
  RotMatrix[0][0] =  cy*cz;
  RotMatrix[0][1] =  cy*sz;
  RotMatrix[0][2] =  sy;
  RotMatrix[1][0] = -sz;
  RotMatrix[1][1] =  cz;
  RotMatrix[1][2] =  0.0;
  RotMatrix[2][0] = -sy*cz;
  RotMatrix[2][1] = -sy*sz;
  RotMatrix[2][2] =  cy;

  /* cutoff very small matrix elements */
  for(i=0; i<3; i++)
    for(j=0; j<3; j++)
      if(fabs(RotMatrix[i][j]) < 1e-12) RotMatrix[i][j] = 0.0;
}


/****************************************************************************************/
/*  General I/O Functions                                                               */
/****************************************************************************************/

/* fileOpen open file 'name' and gives pointer back 
   in case of an opening error, a message is written to the LogFile */

FILE * fileOpen(const char *name, char *mode)
{	
	FILE *f=NULL;
	
	f = fopen(name, mode);
	if (f==NULL) 
	{	fprintf(LogFilePtr, "ERROR: Can't open %s!\n", name);
		exit(-1);
	}
	return f;
}


void Error(const char *text)
{
	fprintf(LogFilePtr,"ERROR: %s!\n", text);
	exit(-1);
}


void Warning(const char *text)
{
	fprintf(LogFilePtr,"Warning: %s!\n", text);
}


void Abort()
{
	exit(-1);
}



/****************************************************************************************/
/*  Functions for Reading of Input Data                                                 */
/****************************************************************************************/

/* ReadLine reads next line from file 'pFile' into string 'pLine' that is 
     not empty and not a comment line (beginning with #)
	  returning TRUE if line is found and FALSE otherwise
   it strips comments at the end, leading and succeeding blanks, line feeds, tabs anc cr 
   the maximal number of characters in the string must be given in 'nStrLen'
*/
int
ReadLine(FILE* pFile, char* pLine, int nStrLen)
{
	char *pComment;
	short k, kmax;

	strcpy(pLine, "");
	if (pFile!=NULL)
	{
		while(strlen(pLine)==0  && !feof(pFile))	
		{	
			fgets (pLine, nStrLen, pFile);

			/* delete line feeds, tabs and carriage returns */
			kmax = (short) strlen(pLine);
			for (k=0; k < kmax; k++)
			{	if (pLine[k]=='\n' || pLine[k]=='\t' || pLine[k]=='\r')
					pLine[k]=' ';
			}
			/* strip the comments and leading and succeeding blanks */
			pComment = strchr(pLine, '#');
			if (pComment != NULL)
				*pComment = '\0';
			while (pLine[0]==' ')
			{	StrgLShift(pLine,1);
			}
			while (pLine[strlen(pLine)-1]==' ')
			{	pLine[strlen(pLine)-1]='\0';
			}
		}
	}
	if (strlen(pLine) > 0)
		return TRUE;
	else
		return FALSE;
}


/*  ReadParString(FILE *fpt) reads one string value from parameter file */

void ReadParString(FILE *fpt, char *stringvar)
{
	fscanf(fpt,"%s", stringvar ) ;

	return ;
}


/*  ReadParF(FILE *fpt) reads one double value from parameter file */

double ReadParF(FILE *fpt)
{
	double value ;
	value=0. ;
	fscanf(fpt,"%lf", &value ) ;
	return value;
}


/*  ReadParI(FILE *fpt) reads one integer value from parameter file */

int ReadParI(FILE *fpt)
{
	int value ;
	value=0 ;
	fscanf(fpt,"%d", &value ) ;
	return value;
}


/* ReadParComment(FILE *fpt) reads comment line */

void ReadParComment(FILE *fpt)
{
	char comment[100], *c ;
	c=fgets(comment, 100, fpt) ;
}


/**********************************************************/
/*  String Operations                                     */
/**********************************************************/

/* Copy 'nLen' bytes of 'sOrigin' into the new string 'sCopy' */
void
StrgCopy(char* sCopy, const char* sOrigin, int nLen)
{
	strncpy(sCopy, sOrigin, nLen);
	sCopy[nLen]='\0';
}


/* Shift string 'sStr' 'kWidth' bytes to the left */
void
StrgLShift(char* sStr, int kWidth)
{
	int k, ke;

	ke = strlen(sStr) - kWidth;

	for (k=0; k <= ke; k++)
		sStr[k] = sStr[k+kWidth];
}


/* Scan string 'sStr' and copy all values (but maximally 'nMax') 
   to list 'pTab' of double values,  beginning with value number 'nStart'*/
long
StrgScanLF(const char* sStr, double* pTab, const int nMax, const int nStart)
{
	int    k, n=0;
	char   *pStr, sNumber[31];

	pStr = (char*) sStr;
	n   -= nStart;
	do
	{	/* search of beginning and end of 1st number of (remaining) string */
		k=0;
		/* step forward until first number or control character */
		while (isdigit(pStr[k])==0 && iscntrl(pStr[k])==0) 
			k++; 
		/* step forward until space-like or control character */
		while (isspace(pStr[k])==0 && iscntrl(pStr[k])==0) 
			k++;  

		/* separating first number and adding it to the list */
		if (k > 0)
		{	
			StrgCopy(sNumber, pStr, k);
			if (n >= 0)
				pTab[n] = atof(sNumber);
			n++;	
			pStr += k;
		}
	}
	while (n < nMax && k > 0);

	return(n);
}


/********************************************************************************************/
/* Intersection.c (merger of YTSDefs.c and vectan.c)                                        */
/* Functions that calculate intersection points with various surfaces                       */
/*                                                                                          */
/* The free non-commercial use of these routines is granted providing due credit is given   */
/* to the authors:                                                                          */
/* Friedrich Streffer, Géza Zsigmond, Dietmar Wechsler,                                     */
/* Michael Fromme, Klaus Lieutenant, Sergey Manoshin                                        */ 
/*                                                                                          */
/* Change: K. L. 2002 JAN, reorganized routines                                                 */

#ifdef VT_GRAPH
# include "cpgplot.h"
extern int do_visualise;
#endif

#include <string.h>

#ifdef VITESS
 #include "intersection.h"
#else
#ifndef INTERSECTION_H
#define INTERSECTION_H

/*********************************************************/
/* intersection.h                                        */
/* Functions that calculate intersection points with     */
/* various surfaces                                      */
/*********************************************************/

#ifdef VITESS
 #include "general.h"
#else
#ifndef GENERAL_H
#define GENERAL_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/******************************/
/** Definitions              **/
/******************************/

#ifdef _MSC_VER
#define M_PI            3.14159265358979323846  /* pi */ 
#define M_PI_2          1.57079632679489661923  /* pi/2 */
#endif

#define MN          1.6749284E-27
#define G           9.80665 
#define K           1.380662E-23
#define NA          6.022137E23
#define H           6.626076E-34
#define L_2_E       81805.048
#define E_C         1.6021773E-19

#define TRUE 		1
#define FALSE 		0

#define UP          1
#define DOWN        0

#define ON          1
#define OFF         0

#define GUIDEFLIGHT 1

#define MAX_COLLISIONS      100
#define MAX_CHOPPER_WINDOWS  10
#define LAMBDA_MIN            0.001
#define LAMBDA_MAX          100.0

#define BUFFER_SIZE       10000
#ifndef CHAR_BUF_LENGTH
#define CHAR_BUF_LENGTH    1024
#endif
#define CHAR_BUF_LARGE     5120
#define CHAR_BUF_SMALL      256

#ifdef VITESS
 #define ABSORB continue
#endif

typedef enum 
{	VT_CUBE   = 1,
	VT_CYL    = 2,
	VT_SPHERE = 3
}
SampleGeom;


typedef enum
{
	VT_SOURCE      =   1,
	VT_GUIDE       =  11,
	VT_BENDER      =  12,
	VT_COLLIMATOR  =  13,
	VT_SM_ENSEMBLE =  15,
	VT_SPACE       =  20,
	VT_WINDOW      =  21,
	VT_WND_MULT    =  22,
	VT_GRID        =  23,
	VT_CHOP_DISC   =  31,
	VT_CHOP_FERMI  =  32,
	VT_VEL_SELECT  =  41,
	VT_MONOC_ANALY =  45,
	VT_POL_HE3     =  51,
	VT_POL_SM      =  52,
	VT_POL_MIRROR  =  53,
	VT_FLIP_COIL   =  55,
	VT_FLIP_GRAD   =  56,
	VT_RES_DRABKIN =  59,
	VT_PREC_FIELD  =  60,
	VT_ROT_FIELD   =  61,
	VT_DETECTOR    =  71,
	VT_SMPL_EL_ISO =  81,
	VT_SMPL_INELAST=  83,
	VT_SMPL_SING_X =  84,
	VT_SMPL_POWDER =  85,
	VT_SMPL_S_Q    =  86,
	VT_SMPL_SANS   =  87,
	VT_SMPL_REFL   =  89,
	VT_MONITOR_1   = 101,
	VT_MONITOR_2   = 102,
	VT_MON_POL_1   = 103,
	VT_MON_POL_POS = 104,
	VT_EVAL_ELAST  = 111,
	VT_EVAL_INELAST= 112,
	VT_VISUAL      = 121,
	VT_FRAME       = 131,
	VT_WRITEOUT    = 141,
	VT_TOOL        = 999
}
VtModID;


typedef double VectorType[3];
typedef double DoublePair[2];


/******************************/
/** Structures               **/
/******************************/

typedef struct
{
	double X,Y,Z;
}
CartesianPoint;


typedef struct
{
	double	A, B, C, D;
}
Plane;


typedef struct
{
	double  A, B, C, D, E, F, W, P, Q, R;
}
SurfaceSecond;


typedef struct
{
	char           IDGrp[2];
	unsigned long  IDNo;
}
TotalID;


typedef struct 
{
	TotalID        ID;
	char           Debug;
	short          Color;
	double         Time;
	double         Wavelength;
	double         Probability;
	VectorType     Position;
	VectorType     Vector;
	VectorType     Spin;
}
Neutron;


typedef struct 
{
      double height, r;
} 
CylinderType;


typedef struct 
{
      double height, width, thickness;
} 
CubeType;


typedef struct 
{
      double r;
} 
BallType;


typedef union 
{
    CylinderType Cyl;
    CubeType     Cube;
    BallType     Ball;
} 
SampleGeomType;


typedef struct 
{
  SampleGeom Type;
  VectorType Position;
  VectorType Direction;
  SampleGeomType SG;
} 
SampleType;

typedef struct
{
	VtModID  eModule;
	double   dWPar;    /* width, ...             */
	double   dHPar;    /* height, end width, ... */
	double   dRPar;    /* radius, ...            */
	long     nNumber;  /* number of ....         */
	short    eType;    /* shape, mon. par., ...  */
	char*    pDescr;   /* material, ...          */
}
ModProp;



/******************************/
/** Prototypes               **/
/******************************/

double ENERGY_FROM_LAMBDA(double x);
double LAMBDA_FROM_ENERGY(double x);
double ENERGY_FROM_V   (double x);    
double V_FROM_LAMBDA   (double x);
double LAMBDA_FROM_V(double x);

double ran3       (long * i);
double MonteCarlo (double x, double y);
double vector3rand(double *, double *, double *);

double sq   (double Value);                        /* = Value*Value*/
double atan0(double a, double b);
void   Exchange(double* pValue1, double* pValue2);
double Min(double value1, double value2);
double Max(double value1, double value2);
long   mini(long value1, long value2);
long   maxi(long value1, long value2);

double SolidAngle   (const double dHorAngle, const double dVertAngle);

void   CopyVector   (const VectorType Src, VectorType Dest);
long   MAXV         (const VectorType Vector);
double LengthVector (const VectorType Vector);
double DistVector   (const VectorType Vec1, const VectorType Vec2);
double ScalarProduct(const VectorType Vec1, const VectorType Vec2);
double AngleVectors (VectorType v1, VectorType v2);
double Area(VectorType v1, VectorType v2);
short  NormVector      (VectorType Vector);
void   AddVector       (VectorType Value,  const VectorType Add);
void   SubVector       (VectorType Value,  const VectorType Sub);
void   MultiplyByScalar(VectorType Vector, const double Scalar);
void   RotVector       (double RotMatrix[3][3], VectorType Vector);
void   RotBackVector   (double RotMatrix[3][3], VectorType Vector);
void   FillRMatrixZY   (double RotMatrix[3][3], double roty, double rotz);

FILE * fileOpen(const char *name, char *mode);
void   Error(const char *text);
void   Warning(const char *text);
void   Abort();

int    ReadLine(FILE* pFile, char* pLine, int nStrLen);
void   ReadParString(FILE *fpt, char *stringvar);
double ReadParF(FILE *fpt);
int    ReadParI(FILE *fpt);
void   ReadParComment(FILE *fpt);

void   StrgCopy  (char* sCopy, const char* sOrigin, int nLen);
void   StrgLShift(char* sStr, int kWidth);
long   StrgScanLF(const char* sStr, double* pTable, const int nMax, const int nStart);

#endif


#endif

#define X(x) ISP[x][0]
#define Y(x) ISP[x][1]
#define Z(x) ISP[x][2]


/* for function 'PathThroughBenderGravOrder2' in module bender */
/* ----------------------------------------------------------- */
double NeutronPlaneAngle2             (Neutron *, double, double, double);
double NeutronSurfaceSecIntersectionGr(Neutron *, SurfaceSecond, long);

/* for several modules */
/* ------------------- */
double NeutronPlaneIntersectionGrav(Neutron *, Plane);
double NeutronPlaneIntersection1   (Neutron *, Plane);

/* for modules 'sample_sans', 'sample_powder', 'monochr_analyser' etc. */
/* ------------------------------------------------------------------- */
int PlaneLineIntersect (VectorType LineOffset, VectorType LineDir, VectorType PlaneNormalVector, double PlaneDistane, VectorType Result);
int PlaneLineIntersect2(VectorType LineOffset, VectorType LineDir, VectorType PlaneNormalVector, double PlaneDistane, VectorType Result);
long IntersectionWithHorizontalPlane(double Z , VectorType PosVect , VectorType Dir, VectorType Result);
int  OrderPositions(VectorType Dir, VectorType Pos1, VectorType Pos2);

long IntersectionWithRectangular(VectorType DimSample, VectorType Pos, VectorType Dir, VectorType Pos1, VectorType Pos2);

long LineIntersectsCube    (VectorType Offset, VectorType Direction, CubeType   *Cube,  double t[2]);
long LineIntersectsCylinder(VectorType Offset, VectorType Direction, CylinderType *Cyl, double t[2]);
long IntersectionWithInfiniteCylinder(double DiameterCyl, VectorType Pos, VectorType Dir, VectorType Pos1, VectorType Pos2);
long IntersectionWithCylinder(VectorType DimSample, VectorType Pos, VectorType Dir, VectorType Pos1, VectorType Pos2);
long LineIntersectsSphere  (VectorType Offset, VectorType Direction, BallType  *Sphere, double t[2]);
long IntersectionWithSphere(VectorType DimSample, VectorType Pos, VectorType Dir, VectorType Pos1, VectorType Pos2);

#endif

#endif


/***********************************************************************************/
/* global variables and constants                                                  */
/***********************************************************************************/

VectorType PVector[6] = {{1.0, 0.0, 0.0}, {-1.0, 0.0, 0.0}, {0.0, 1.0, 0.0},
                         {0.0,-1.0, 0.0}, { 0.0, 0.0, 1.0}, {0.0, 0.0,-1.0}};


/***********************************************************************************/
/* prototypes of local functions                                                   */
/***********************************************************************************/

double SolveQuadraticEq(double, double, double);



/***********************************************************************************/
/* global functions                                                                */
/***********************************************************************************/

/***************************************************************************************************/
/* This Function similar NeutronPlaneIntersection, but he is calculate new position of neutron
   in the plane (update) and RETURN time of flight of neutron
   Author: Manoshin Sergey, manoshin@hmi.de  20.02.01 */
/**************************************************************************************************/
double NeutronPlaneIntersection1(Neutron *ThisNeutron, Plane ThisPlane)
{
	double  Time, BB, CC, VelocityReal;
		
	/*	Velocity cm/ms, Time ms, Position cm	*/

	/*	Calculating real velocity  */
		
	VelocityReal = (double)(V_FROM_LAMBDA(ThisNeutron->Wavelength));


	BB = ThisPlane.A*ThisNeutron->Vector[0]
	  +ThisPlane.B*ThisNeutron->Vector[1]
	  +ThisPlane.C*ThisNeutron->Vector[2];
		
	BB = BB*VelocityReal;

	CC = ThisPlane.A*ThisNeutron->Position[0]
	  +ThisPlane.B*ThisNeutron->Position[1]
	  +ThisPlane.C*ThisNeutron->Position[2]
	  +ThisPlane.D;

	/***********************************************************************************/
	/* Now we must to decide equation for find */
	/* Time BB*Time + CC = 0 */
	/***********************************************************************************/
	if(BB != 0.0)
	{
		Time = -CC/BB;
	}
	else
	{
		Time = -1.0e4;
	}

	if (fabs(Time) < 1e-12)
		Time = 0.0;
		
	/***********************************************************************************/
	/* This part make move neutron in the time include gravity effect  */
	/***********************************************************************************/

	ThisNeutron->Position[0] = ThisNeutron->Position[0] + VelocityReal*Time*(ThisNeutron->Vector[0]);
	ThisNeutron->Position[1] = ThisNeutron->Position[1] + VelocityReal*Time*(ThisNeutron->Vector[1]);
	ThisNeutron->Position[2] = ThisNeutron->Position[2] + VelocityReal*Time*(ThisNeutron->Vector[2]);

	/* return Time; ms */

	return Time;
}


/***************************************************************************************************/
/* This Function similar NeutronPlaneIntersection, but he is include GRAVITY effect
   and calculate new position of neutron in the plane (update) and time of flight of neutron (return)
   Author: Manoshin Sergey, manoshin@hmi.de  12.02.01 */
/**************************************************************************************************/
double NeutronPlaneIntersectionGrav(Neutron *ThisNeutron, Plane ThisPlane)
{
	/***********************************************************************************/
	/* This part calculates the time of flight of neutron INCLUDE gravity with   */
	/* plane.									   */
	/***********************************************************************************/

	double  Time, AA, BB, CC, VelocityReal;
		
	/*	Make the koefficients of quadratic equation    */
	/*	G = 9.8 m/c**2, we need to cm/ms**2 (100/1000/1000)  */
	/*	Velocity cm/ms, Time ms, Position cm	*/

	/*	Calculating real velocity  */
		
	VelocityReal = (double)(V_FROM_LAMBDA(ThisNeutron->Wavelength));

	AA = -0.5*(G*1.0e-4)*ThisPlane.C;

	BB = ThisPlane.A*ThisNeutron->Vector[0]
	  +ThisPlane.B*ThisNeutron->Vector[1]
	  +ThisPlane.C*ThisNeutron->Vector[2];
		
	BB = BB*VelocityReal;

	CC = ThisPlane.A*ThisNeutron->Position[0]
	  +ThisPlane.B*ThisNeutron->Position[1]
	  +ThisPlane.C*ThisNeutron->Position[2]
	  +ThisPlane.D;

	/***********************************************************************************/
	/* Now we must to decide quadratic equation for find */
	/* Time AA*Time*Time + BB*Time + CC = 0 */
	/***********************************************************************************/

	Time = SolveQuadraticEq(AA,BB,CC);
		
	/***********************************************************************************/
	/* This part make move neutron in the time include gravity effect  */
	/***********************************************************************************/

	ThisNeutron->Position[0] = ThisNeutron->Position[0] + VelocityReal*Time*(ThisNeutron->Vector[0]);
	ThisNeutron->Position[1] = ThisNeutron->Position[1] + VelocityReal*Time*(ThisNeutron->Vector[1]);
	ThisNeutron->Position[2] = ThisNeutron->Position[2] + VelocityReal*Time*(ThisNeutron->Vector[2]);

	/* Include gravity */

	ThisNeutron->Position[2] = ThisNeutron->Position[2] - 0.5*(G*1.0e-4)*Time*Time;
	ThisNeutron->Vector[2] = ThisNeutron->Vector[2] - ((G*1.0e-4)*Time/VelocityReal);	
		
	/* return Time; ms */

	return Time;
}



/***********************************************************************************/
/* This subroutine calculates the angle between a directed line and a plane. The   */
/* formulae used can be found in any good mathematical handbook.                   */
/* (Modified by Manoshin Sergey for module of velocity not equal 1)                */	
/***********************************************************************************/
double	NeutronPlaneAngle2(Neutron *ThisNeutron, double AP, double BP, double CP)
{
	/* Revised calling parameter */

	double SinGamma=0.0;
	double Velmod;
	double ABC;
		
		
	Velmod = sqrt(ThisNeutron->Vector[0]*ThisNeutron->Vector[0]+
			ThisNeutron->Vector[1]*ThisNeutron->Vector[1] +
			ThisNeutron->Vector[2]*ThisNeutron->Vector[2]);
		
	ABC = sqrt(AP*AP + BP*BP + CP*CP);

	if (ABC != 0.0)
	{
		SinGamma = (AP*ThisNeutron->Vector[0] + BP*ThisNeutron->Vector[1] + CP*ThisNeutron->Vector[2])/ABC;
	}
		
					
	if (Velmod != 0.0)
	{				
		SinGamma = SinGamma/Velmod;
	}

	return(asin(SinGamma));
}



/***************************************************************************************************/
/* This function move neutron from current position to the surface , furthemore he is include
   GRAVITY effect  and time of flight of neutron (return)
   Author: Manoshin Sergey, manoshin@hmi.de  25.03.01 */
/**************************************************************************************************/
double NeutronSurfaceSecIntersectionGr(Neutron *ThisNeutron, SurfaceSecond ThisSurfaceSecond, long keygrav)
{
	/***********************************************************************************/
	/* This part calculates the time of flight of neutron INCLUDE gravity with   */
	/* plane.	This functions is core for transporting neutrons!									   */
	/***********************************************************************************/

	double  Time, AA, BB, CC, VelocityReal;
	double  VX, VY, VZ, X, Y, Z;
		
	/*	Make the koefficients of quadratic equation    */
	/*	G = 9.8 m/c**2, we need to cm/ms**2 (100/1000/1000)  */
	/*	Velocity cm/ms, Time ms, Position cm	*/

	/*	Calculating real velocity, cm/ms  */
		
	VelocityReal = (double)(V_FROM_LAMBDA(ThisNeutron->Wavelength));
		
	/*	Components of velocity, projections, and position coordinats */
	/*	Local copy */
		
	X = ThisNeutron->Position[0];
	Y = ThisNeutron->Position[1];
	Z = ThisNeutron->Position[2];
		
	VX = VelocityReal*ThisNeutron->Vector[0];
	VY = VelocityReal*ThisNeutron->Vector[1];
	VZ = VelocityReal*ThisNeutron->Vector[2];
		
	/*	Find the coefficients of the quadratic equation */	
		
	if (keygrav == 1)
	{
	
/* Corrected: Marz 03 */	

		AA = -0.5*(G*1.0e-4)*ThisSurfaceSecond.F;
	}
	else
	{
		AA = ThisSurfaceSecond.A*VX*VX +
		ThisSurfaceSecond.C*VY*VY +
		ThisSurfaceSecond.E*VZ*VZ +
		ThisSurfaceSecond.P*VX*VY +
		ThisSurfaceSecond.Q*VY*VZ +
		ThisSurfaceSecond.R*VX*VZ;
	}

  BB = 2.0*(ThisSurfaceSecond.A*VX*X + ThisSurfaceSecond.C*VY*Y + ThisSurfaceSecond.E*VZ*Z) +
    ThisSurfaceSecond.B*VX + ThisSurfaceSecond.D*VY + ThisSurfaceSecond.F*VZ +
    ThisSurfaceSecond.P*(X*VY+Y*VX) +
    ThisSurfaceSecond.Q*(Y*VZ+Z*VY) +
    ThisSurfaceSecond.R*(X*VZ+Z*VX);
	

	CC = ThisSurfaceSecond.A*X*X + ThisSurfaceSecond.B*X+
	  ThisSurfaceSecond.C*Y*Y + ThisSurfaceSecond.D*Y+
	  ThisSurfaceSecond.E*Z*Z + ThisSurfaceSecond.F*Z + ThisSurfaceSecond.W +
	  ThisSurfaceSecond.P*X*Y + ThisSurfaceSecond.Q*Y*Z + ThisSurfaceSecond.R*X*Z;



	/***********************************************************************************/
	/* Now we must to decide quadratic equation for find */
	/* Time AA*Time*Time + BB*Time + CC = 0 */
	/***********************************************************************************/

	Time = SolveQuadraticEq(AA,BB,CC);
		
	/***********************************************************************************/
	/* This part make move neutron in the time include gravity effect  */
	/***********************************************************************************/

	ThisNeutron->Position[0] = ThisNeutron->Position[0] + VelocityReal*Time*(ThisNeutron->Vector[0]);
	ThisNeutron->Position[1] = ThisNeutron->Position[1] + VelocityReal*Time*(ThisNeutron->Vector[1]);
	ThisNeutron->Position[2] = ThisNeutron->Position[2] + VelocityReal*Time*(ThisNeutron->Vector[2]);

	/* Include gravity */
	if (keygrav == 1)
	{	
		ThisNeutron->Position[2] = ThisNeutron->Position[2] - 0.5*(G*1.0e-4)*Time*Time;
		ThisNeutron->Vector[2] = ThisNeutron->Vector[2] - ((G*1.0e-4)*Time/VelocityReal);	
	}
	/* return Time; ms */

	return Time;
}

/* ordering points on a trajectory so that Dir shows from Pos1 to Pos2 */

int OrderPositions(VectorType Dir, VectorType Pos1, VectorType Pos2) 
{
VectorType propag, V;

CopyVector(Pos2, propag); SubVector(propag, Pos1);

	if(ScalarProduct(propag, Dir) == 0.) return 0;
	if(ScalarProduct(propag, Dir) < 0.)
	{
		CopyVector(Pos1, V) ;
		CopyVector(Pos2, Pos1) ;
		CopyVector(V, Pos2) ;
	}
	return 1;
}

/******************************************************************************/
/* 'IntersectionWithRectangular'                                              */
/* 'LineIntersectsCube'                                                       */
/* 'LineIntersectsCylinder'                                                   */
/* 'LineIntersectsSphere'                                                     */
/* These function test whether a line defined by the vectors 'Offset'         */
/* and 'Direction' intersect a solid of certain geometry                      */
/* If the line intersects, the functions returns TRUE                         */
/* For 'LineIntersectsCube', 'LineIntersectsCylinder', 'LineIntersectsSphere' */
/* the parameters t[0] and t[1] are calculated                                */
/* (t[i] are proprotional to the time-of-flight until the intersection)       */
/******************************************************************************/

/******************************************************************************/
/* 'IntersectionWithRectangular'      intersection function for a rectangular */
/* (Author: G. Zsigmond)                                                                     */
long IntersectionWithRectangular(VectorType DimSample, VectorType Pos, VectorType Dir, VectorType Pos1, VectorType Pos2)
{
	VectorType	n, pos0, pos1, pos2, pos3, pos4, pos5 ;
	int			k ;

	for(k=0;k<3;k++) Pos1[k] = Pos2[k] = 0. ;

	n[0] = 1. ; n[1] = n[2] = 0. ;

	if(PlaneLineIntersect2(Pos, Dir, n, - DimSample[0]/2, pos0) == 1)
	{
		if( /*(fabs(pos0[0]) <= DimSample[0]/2.) && */(fabs(pos0[1]) <= DimSample[1]/2.) && (fabs(pos0[2]) <= DimSample[2]/2.) ) CopyVector(pos0, Pos1) ;
	}
	if(PlaneLineIntersect2(Pos, Dir, n, + DimSample[0]/2, pos1) == 1)
	{
		if( /*(fabs(pos1[0]) <= DimSample[0]/2) && */(fabs(pos1[1]) <= DimSample[1]/2) && (fabs(pos1[2]) <= DimSample[2]/2) )
		{
			if(LengthVector(Pos1) == 0.) CopyVector(pos1, Pos1) ;
			else CopyVector(pos1, Pos2) ;
		}
	}

	n[1] = 1. ; n[2] = n[0] = 0. ;

	if(PlaneLineIntersect2(Pos, Dir, n, - DimSample[1]/2, pos2) == 1)
	{
		if( (fabs(pos2[0]) <= DimSample[0]/2) && /*(fabs(pos2[1]) <= DimSample[1]/2) &&*/ (fabs(pos2[2]) <= DimSample[2]/2) )
		{
			if(LengthVector(Pos1) == 0.) CopyVector(pos2, Pos1) ;
			else CopyVector(pos2, Pos2) ;
		}
	}
	if(PlaneLineIntersect2(Pos, Dir, n, + DimSample[1]/2, pos3) == 1)
	{
		if( (fabs(pos3[0]) <= DimSample[0]/2) && /*(fabs(pos3[1]) <= DimSample[1]/2) &&*/ (fabs(pos3[2]) <= DimSample[2]/2) )
		{
			if(LengthVector(Pos1) == 0.) CopyVector(pos3, Pos1) ;
			else CopyVector(pos3, Pos2) ;
		}
	}

	n[2] = 1. ; n[0] = n[1] = 0. ;

	if(PlaneLineIntersect2(Pos, Dir, n, - DimSample[2]/2, pos4) == 1)
	{
		if( (fabs(pos4[0]) <= DimSample[0]/2) && (fabs(pos4[1]) <= DimSample[1]/2) /*&& (fabs(pos4[2]) <= DimSample[2]/2)*/ )
		{
			if(LengthVector(Pos1) == 0.) CopyVector(pos4, Pos1) ;
			else CopyVector(pos4, Pos2) ;
		}
	}
	if(PlaneLineIntersect2(Pos, Dir, n, + DimSample[2]/2, pos5) == 1)
	{
		if( (fabs(pos5[0]) <= DimSample[0]/2) && (fabs(pos5[1]) <= DimSample[1]/2) /*&& (fabs(pos5[2]) <= DimSample[2]/2)*/ )
		{
			if(LengthVector(Pos1) == 0.) CopyVector(pos5, Pos1) ;
			else CopyVector(pos5, Pos2) ;
		}
	}
	/*pi(0);pv(pos0);pv(pos1);pv(pos2);pv(pos3);pv(pos4);pv(pos5);*/

	if((LengthVector(Pos1) == 0.) || (LengthVector(Pos2) == 0.)) return 0 ;

	/*ordering intersection positions */
	
	if(OrderPositions(Dir, Pos1, Pos2)==1);

	return 1 ;

}/* End IntersectionWithRectangular() */



/*************************************************************************/
/* 'LineIntersectsCube'    intersection function for a cube              */
/* The function assumes that the corners of the cube are parallel        */
/* to the x-,y- and z-axis and the center of the cube resides at (0,0,0) */
/* (Author: F. Streffer)                                                                 */
long LineIntersectsCube(VectorType Offset, VectorType Direction, CubeType *Cube, double t[2])
{
	long i,j,k;
	double     PDist[6];
	VectorType ISP[2],
				  TestV; 

	PDist[0] = 0.5*Cube->thickness;
	PDist[1] = 0.5*Cube->thickness;
	PDist[2] = 0.5*Cube->width;    
	PDist[3] = 0.5*Cube->width;
	PDist[4] = 0.5*Cube->height,    
	PDist[5] = 0.5*Cube->height;

	TestV[0] = 0.500000000001*Cube->thickness;   /* not equal 0.5*... to prevent rounding errors */
	TestV[1] = 0.500000000001*Cube->width;
	TestV[2] = 0.500000000001*Cube->height;

	k=0;
	i=MAXV(Direction);

	for(j=0;k<2 && j<6; j++) 
	{
		PlaneLineIntersect(Offset, Direction, PVector[j], PDist[j], ISP[k]);
		if((fabs(ISP[k][0])<TestV[0]) &&
		   (fabs(ISP[k][1])<TestV[1]) &&
		   (fabs(ISP[k][2])<TestV[2])) 
		{
			t[k]=(ISP[k][i]-Offset[i])/Direction[i];
			k++;
		}
	}
	/* Ordering: t1 > t0 */
	if(k>1) 
	{	if (t[0] > t[1])
			Exchange(&t[0], &t[1]);
		return TRUE;
	}
	else
	{	return FALSE;
	}
}



/*************************************************************************/
/* 'LineIntersectsCylinders'    intersection function for a cylinder     */
/* the cylinder axis points a long the x-axis                            */
/* (Author: F. Streffer)                                                                 */
long LineIntersectsCylinder(VectorType Offset, VectorType Direction,
                              CylinderType *Cyl, double t[2])
{
	VectorType ISP[2];
	double p=0.0,q,
	       sp, spq, spr,
	       h;           /* half of the cylinder height*/
	long   result, i, iDir;

	/* determine the intersection of the line with an infinite cylinder */
	h   = Cyl->height/2.0;
	spq = Direction[1]*Direction[1] + Direction[2]*Direction[2];
	if(spq > 0.0) 
	{
		p = (Offset[1]*Direction[1] + Offset[2]*Direction[2]) / spq;
		q = (Offset[1]*Offset[1] + Offset[2]*Offset[2] - Cyl->r*Cyl->r) / spq;
		sp = p*p-q;
	} 
	else 
	{
		sp = 0.0;
	}

	if(sp > 0.0) 
	{
		/*Ok, the neutron intersects the infinite cylinder */
		spr=sqrt(sp);
		t[0]=-p-spr;
		t[1]=-p+spr;

		/* Give the t's an order so its easier later on */
		if(t[0] > t[1]) 
		  Exchange(&t[0], &t[1]);

		/* checks, where the line enters and leaves the cylinder */
		/* X[i] = ISP[i][0] (see intersection.h) */
		X(0) = Offset[0]+t[0]*Direction[0]; 
		X(1) = Offset[0]+t[1]*Direction[0];
		if (Direction[0] > 0.0) iDir = 1;
		else                    iDir = 0;
		i=MAXV(Direction);

		/* A) line enters and leaves through the cylinder wall */
		if     (fabs(X(0)) < h  &&  fabs(X(1)) < h) 
		{
			result=TRUE;
		}
		/* B) line enters through the cylinder wall and leaves through one of the cylinder disks */
		else if(fabs(X(0)) < h  &&  fabs(X(1)) >= h) 
		{
			PlaneLineIntersect(Offset, Direction, PVector[1-iDir], h, ISP[1]);
			t[1] = (ISP[1][i]-Offset[i])/Direction[i];
			result=TRUE;
		} 
		/* C) lines enters through a cylinder disk and leaves through the cylinder wall */
		else if(fabs(X(0)) >= h  &&  fabs(X(1)) < h) 
		{
			PlaneLineIntersect(Offset, Direction, PVector[iDir], h, ISP[0]);
			t[0] = (ISP[0][i]-Offset[i])/Direction[i];
			result=TRUE;
		} 
		/* D) line enters through one of the cylinder disks and leaves through the other one */
		else if( (X(0) < -h  &&  X(1) > h) || (X(0) > h  &&  X(1) < -h) ) 
		{
			PlaneLineIntersect(Offset, Direction, PVector[1-iDir], h, ISP[1]);
			PlaneLineIntersect(Offset, Direction, PVector[iDir],   h, ISP[0]);
			t[0] = (ISP[0][i]-Offset[i])/Direction[i];
			t[1] = (ISP[1][i]-Offset[i])/Direction[i];
			result=TRUE;
		}
		else 
		{
			result=FALSE;
		}
	} 
	else 
	{	/* Hmm, just missed, but wait one chance remains */
		result=FALSE;
		if(spq==0.0 && (Offset[1]*Offset[1] + Offset[2]*Offset[2] < Cyl->r*Cyl->r))
		{
			/* ok, really not very propable but ... */
			/* the neutron travels paralles to x-axis inside the cylinder   */
			t[0] = ( h-Offset[0])/Direction[0];
			t[1] = (-h-Offset[0])/Direction[0];
			result=TRUE;
		}
	}

	/* Ordering: t1 > t0 */
	if (result==TRUE  &&  t[0] > t[1])
		Exchange(&t[0], &t[1]);

	return result;
}

/************************************************************************/
/* IntersectionWithInfiniteCylinder: gives coordinates of intersections  */
/* (Author: G. Zsigmond)                                                                */
long	IntersectionWithInfiniteCylinder(double DiameterCyl, VectorType Pos, VectorType Dir, VectorType Pos1, VectorType Pos2)
{
double b, c, delta;

		b = (Pos[0] * Dir[0] + Pos[1] * Dir[1]) / (sq(Dir[0]) + sq(Dir[1])) ;
		c = (sq(Pos[0]) + sq(Pos[1]) - sq(DiameterCyl/2)) / (sq(Dir[0]) + sq(Dir[1])) ;

		if((delta = sq(b) - c) < 0.) return 0 ;

		CopyVector(Dir, Pos1) ;
		MultiplyByScalar(Pos1, - b + (double) sqrt(delta)) ;
		AddVector(Pos1, Pos) ;

		CopyVector(Dir, Pos2) ;
		MultiplyByScalar(Pos2, - b - (double) sqrt(delta)) ;
		AddVector(Pos2, Pos) ;

		return 1 ;
}

/*****************************************************************/
/* Intersection with cylinder: gives coordinates of intersections */
/* (Author: G. Zsigmond)                                                           */
long IntersectionWithCylinder(VectorType DimSample, VectorType Pos, VectorType Dir, VectorType Pos1, VectorType Pos2)
{
	if(IntersectionWithInfiniteCylinder(DimSample[0], Pos, Dir, Pos1, Pos2) == 0) return 0 ;

	/* ordering intersection positions */
	
		if(OrderPositions(Dir, Pos1, Pos2)==1);	

	    if((Pos1[2] > DimSample[2]/2.) && (Pos2[2] > DimSample[2]/2.)) return 0 ;
        if((Pos1[2] < - DimSample[2]/2.) && (Pos2[2] < - DimSample[2]/2.)) return 0 ;
        if((fabs(Pos1[2]) <= DimSample[2]/2.) && (fabs(Pos2[2]) <= DimSample[2]/2.)) return 1 ;/* both through cylinder walls */
		
		if( (Pos2[2] <= DimSample[2]/2.) && (Pos2[2] >= - DimSample[2]/2.) ) /* exit cylinder wall */
		{
			if(Pos1[2] >= DimSample[2]/2.)   /* entrance top */
			{
				if(IntersectionWithHorizontalPlane(DimSample[2]/2., Pos2, Dir, Pos1) == 0) return 0 ;
				return 1 ;
			}
			if(Pos1[2] <= - DimSample[2]/2.) /* entrance bottom */
			{
				if(IntersectionWithHorizontalPlane(- DimSample[2]/2., Pos2, Dir, Pos1) == 0) return 0 ;
				return 1;
			}
		else return 0 ;
		}

        if((Pos1[2] >= DimSample[2]/2.) && (Pos2[2] <= - DimSample[2]/2.)) /* entrance top exit bottom */
		{
			if(IntersectionWithHorizontalPlane(DimSample[2]/2., Pos, Dir, Pos1) == 0) return 0 ;
			if(IntersectionWithHorizontalPlane(- DimSample[2]/2., Pos, Dir, Pos2) == 0) return 0 ;
	        return 1;
        }

        if((Pos1[2] <= - DimSample[2]/2.) && (Pos2[2] >= DimSample[2]/2.)) /* entrance bottom exit top */
		{
			if(IntersectionWithHorizontalPlane(- DimSample[2]/2., Pos, Dir, Pos1) == 0) return 0 ;
			if(IntersectionWithHorizontalPlane(DimSample[2]/2., Pos, Dir, Pos2) == 0) return 0 ;
			return 1;
        }

		if( (Pos1[2] >= - DimSample[2]/2.) && (Pos1[2] <= DimSample[2]/2.) ) /* entrance cylinder wall */
		{
			if(Pos2[2] >= DimSample[2]/2.) /* exit top */
			{
				if(IntersectionWithHorizontalPlane(DimSample[2]/2., Pos1, Dir, Pos2) == 0) return 0 ;
				return 1;
			}
			if(Pos2[2] <= - DimSample[2]/2.) /* exit bottom */
			{
				if(IntersectionWithHorizontalPlane(- DimSample[2]/2., Pos1, Dir, Pos2) == 0) return 0 ;
				return 1;
			}
		else return 0 ;
        }
		else return 0 ;
}


/*****************************************************************************/
/* 'LineIntersectsSphere'    intersection function for a sphere              */
/* (Author: F. Streffer)                                                                 */
long LineIntersectsSphere(VectorType Offset, VectorType Direction,
                            BallType *Sphere, double t[2])
/* similar to Line_Intersects_Cube but for a sphere */
{
	/* lets solve x+y+z-r=0 with x=x_0+t*x_t */
	/* x_t+y_t+z_t=1 since Nin->Vector is normalized */
	/* p equals p/2 of the formual solving quadratic eqs. */
	double p,spq;

	p = ScalarProduct(Direction,Offset);
	spq= p*p-(ScalarProduct(Offset,Offset)-Sphere->r*Sphere->r);
	if(spq > 0.0) 
	{
		spq=sqrt(spq);
		t[0]=-p-spq;
		t[1]=-p+spq;

		/* Ordering: t1 > t0 */
		if (t[0] > t[1])
			Exchange(&t[0], &t[1]);
		return TRUE;
	} 
	else
	{	/* Also if spq is 0.0, i.e. the entrance and exit is the same point */
		/* the neutron is considered to have missed the sample */
		return FALSE;
	}
}


/*****************************************************************************/
/* 'IntersectionWithSphere'  gives coordinates of intersections with a sphere */
/* (Author: G. Zsigmond)                                                                     */
long IntersectionWithSphere(VectorType DimSample, VectorType Pos, VectorType Dir, VectorType Pos1, VectorType Pos2)
{
double b, c, delta, khi1, khi2 ;

	b = ScalarProduct(Pos, Dir) ;
	c = ScalarProduct(Pos, Pos) - sq(DimSample[0]/2) ;
	delta = sq(b) - c ;

	khi1 = - b + (double) sqrt(delta) ;
	khi2 = - b - (double) sqrt(delta) ;

	if(delta < 0.) return 0 ;
	else
	{
		CopyVector(Dir, Pos1) ;
		MultiplyByScalar(Pos1, khi1) ;
		AddVector(Pos1, Pos) ;

		CopyVector(Dir, Pos2) ;
		MultiplyByScalar(Pos2, khi2) ;
		AddVector(Pos2, Pos) ;

		/*ordering intersection positions */
	
		if(OrderPositions(Dir, Pos1, Pos2)==1);	

		return 1 ;
	}
}


/****************************************************************/
/* 'PlaneLineIntersect'                                         */
/* 'PlaneLineIntersect2'                                        */
/* These functions compute the intersect of a line give by some */
/* offset and a Direction, and a plane in the Hessian form.     */
/****************************************************************/
/****************************************************************/
/*  PlaneLineIntersect                                          */
/* function returns TRUE if the line intersects and             */
/*                       'Result' will contain the point        */
/*                  FALSE if the line is parallel to plane      */
/* (Author: F. Streffer)                                                                 */
int PlaneLineIntersect(VectorType LineOffset, VectorType LineDir,
                       VectorType PlaneNormalVector, double PlaneDistane,
                       VectorType Result)
{
	double     help,t;
	int        i;
	double help1;

	help=ScalarProduct(LineDir,PlaneNormalVector);
	if(help==0.0) 
	{
		/* Ok, there is somehow a problem, the line given is parallel to   */
		/* the plane and will never intersect.                             */
		return FALSE;
	} 
	else 
	{
		help1= ScalarProduct(LineOffset,PlaneNormalVector);
		t=(PlaneDistane-help1)/help;
		for(i=0;i<3;i++)
			Result[i]=LineOffset[i]+t*LineDir[i];
	}
	return TRUE;
}

/***************************************/
/* IntersectionWithHorizontalPlane		*/
/* Z:		position of the plane		*/
/* PosVect: point on line				*/
/* Dir:		flight direction			*/
/* (Author: G. Zsigmond)                                                           */

long IntersectionWithHorizontalPlane(double Z , VectorType PosVect , VectorType Dir, VectorType Result)
{
	if(Dir[2] ==0.) return 0 ;

	CopyVector(Dir, Result) ;
	MultiplyByScalar(Result, (Z - PosVect[2])/Dir[2]) ;
	AddVector(Result, PosVect) ;

	return 1;
}


/*********************************************************************/
/*  PlaneLineIntersect2                                              */
/*                                                                   */
/* function returns TRUE, 'Result' contains the point                */
/* if line is parallel to plane, 'Result' contains very high values  */
/* (Author: G. Zsigmond)                                                           */
int PlaneLineIntersect2(VectorType LineOffset, VectorType LineDir,
                       VectorType PlaneNormalVector, double PlaneDistane,
                       VectorType Result)
{
	double     help,t;
	int        i;
	double help1;

	help=ScalarProduct(LineDir,PlaneNormalVector);
	if(help==0.0)
	{
		/* Ok, there is somehow a problem, the line given is parallel to   */
		/* the plane and will never intersect.                             */
		for(i=0;i<3;i++) Result[i]=1.e20; return 1;
	} 
	else 
	{
		help1= ScalarProduct(LineOffset,PlaneNormalVector);
		t=(PlaneDistane-help1)/help;
		for(i=0;i<3;i++)
			Result[i]=LineOffset[i]+t*LineDir[i];
	}
	return TRUE;
}


/***********************************************************************************/
/** local functions                                                               **/
/***********************************************************************************/

/*******************************************************************/
/* THIS FUNCTION SOLVES THE QUADRATIC EQUATION
   AA*X*X+BB*X+CC=0
   Version from 29.01.01, Author Manoshin Sergey  manoshin@hmi.de  */
/*******************************************************************/
double SolveQuadraticEq(double AA, double BB, double CC)
{	
	double X1=0.0, X2=0.0;
	double DD, TEMP;
	double Time=0.0;

	/*	ignore a - very small or eq zero , to solve bx+c=0 */
		
	if ((fabs(AA)*1.0e12) <= fabs(BB)||(AA == 0.0))
	{
		if (BB != 0.0)
		{
			/*	AA=0 BB # 0	  */
			X1 = -CC/BB;
			X2 = X1;
			goto ok15;
		}
		else
		{
			if (CC == 0.0)
			{
				/*	AA=0 BB=0 CC=0	    */
				X1 = 0.0;
				X2 = 0.0;
				goto ok15;
			}
			/*	AA=0 BB=0 CC # 0	  */
			X1 = -1.0e4;
			X2 = X1;
			goto ok15;
		}
	}

	/*  later AA not eq zero */

	/*	square eq. Axx+Bx=0 AA # 0 C=0	*/
	if (CC == 0)
	{
		X1 = -BB/AA;
		X2 = 0.0;
		goto ok15;
	}

	/*	Axx+C=0	  */
	if (BB == 0.0)
	{
		TEMP = -CC/AA;
		if (TEMP < 0.0)
		{
			/*	SQRT(TEMP)   */
			X1 = -1.0e4;
			X2 = -1.0e4;
			goto ok15;
		}
		X1 = sqrt(TEMP);
		X2 = -X1;
		goto ok15;
	}


	/*	eq Axx+Bx+C=0	*/
	DD = BB*BB - 4.0*AA*CC;
	if (DD < 0.0)
	{
		/*	D<0   */
		X1 = -1.0e4;
		X2 = -1.0e4;
		goto ok15;
	}
		
	DD = sqrt(DD);
	if (DD == 0.0)
	{
		X1=-BB/(2.0*AA);
		X2=X1;
		goto ok15;
	}
		
		
	if (BB > 0.0)
	{
		X1 = (-BB-DD)/(2.0*AA);
		X2 = (2.0*CC)/(-BB-DD);
		goto ok15;
	}
		
	if (BB < 0.0)
	{
		X2 = (-BB+DD)/(2.0*AA);
		X1 = (2.0*CC)/(-BB+DD);
		goto ok15;
	}
		
	ok15:
	/*	printf("squares= %e  %e \n",X1,X2);	*/
	/* Now find smallest root and copy this in Time  */

	if (fabs(X1) < 4e-10) 
		X1=0.0;
	if (fabs(X2) < 4e-10) 
		X2=0.0;

	if ((X1 > 0.0) && (X2 > 0.0))
	{
		Time = Min(X1, X2);
	}	
	else
	{
		Time = Max(X1, X2);
	}
		
	return Time;
}	
 	

/*******************************************************************************
*
* McStas, neutron ray-tracing package
*         Copyright 1997-2002, All rights reserved
*         Risoe National Laboratory, Roskilde, Denmark
*         Institut Laue Langevin, Grenoble, France
*
* Library: share/vitess-lib.h
*
* %Identification
* Written by: KN, EF
* Date:   Aug 28, 2002
* Origin: Risoe
* Release: McStas 1.6
* Version: $Revision$
*
* This file is to be imported by the mcstas2vitess perl script
* It handles the way Vitess parses parameters.
* Functions are used by the Virtual_input and Virtual_output
* components.
*
* Usage: within SHARE
* %include "vitess-lib"
*
*******************************************************************************/

#ifndef VITESS_LIB_H
#define VITESS_LIB_H "$Revision$"

#ifndef GENERAL_H
 #include <math.h>
 #include <stdlib.h>
 #include <stdio.h>

 /* The Neutron structure, taken from VITESS 2.3 source code "general.h" */
 typedef double VectorType[3];
 typedef struct
 {
  char           IDGrp[2];
  unsigned long  IDNo;
 }
 TotalID;
 typedef struct
 {
  TotalID        ID;
  char           Debug;
  short          Color;
  double         Time;
  double         Wavelength;
  double         Probability;
  VectorType     Position;
  VectorType     Vector;
  VectorType     Spin;
 }
 Neutron;
#endif

extern char *vitess_infile; /* Neutron input file name, or NULL. */
extern char *vitess_outfile;  /* Neutron output file name, or NULL. */
extern int vitess_tracepoints;  /* If true, use dots as progress-indicator */
extern int vitess_repcnt; /* Number of times to repeat this neutron */
extern int vitess_bufsize;  /* The buffer size for neutron read/write */

/* vitess-lib function prototypes */
/* ========================================================================= */
Neutron mcstas2vitess(double x, double y, double z,
                      double vx, double vy, double vz,
                      double t,
                      double sx, double sy, double sz,
                      double p);
void vitess2mcstas(Neutron neu,
                   double *x, double *y, double *z,
                   double *vx, double *vy, double *vz,
                   double *sx, double *sy, double *sz,
                   double *t, double *p);
void vitess_option_error(char *opt);
void vitess_parseopt(int argc, char *argv[],
         double *dptr[], char dchr[], char **sptr[], char schr[]);

void McInitVt();
void McCleanupVt();
void setParDirectory (char *a);
char* FullParName(char* filename);

#endif

/* end of vitess-lib.h */
/*******************************************************************************
*
* McStas, neutron ray-tracing package
*         Copyright 1997-2002, All rights reserved
*         Risoe National Laboratory, Roskilde, Denmark
*         Institut Laue Langevin, Grenoble, France
*
* Library: share/vitess-lib.c
*
* %Identification
* Written by: KN, EF
* Date:   Aug 28, 2002
* Origin: Risoe
* Release: McStas 1.6
* Version: $Revision$
*
* This file is to be imported by the mcstas2vitess perl script
* It handles the way Vitess parses parameters.
* Functions are imported in the Virtual_imput and Virtual_output
* components. Embedded within instrument if MC_EMBEDDED_RUNTIME is defined.
*
* Usage: within SHARE
* %include "vitess-lib"
*
*******************************************************************************/

#ifndef VITESS_LIB_H
#error McStas : please import this library with %include "vitess-lib"
#endif

/********************************************************************************************/
#ifdef WIN32
# include <fcntl.h>
# include <io.h>
# define cSlash '\\'
#else
# define cSlash '/'
#endif

double dTimeMeas   =  0.0,  /* time of measurement           [s]   */
       dLmbdWant   =  0.0,  /* desired wavelength            [Ang] */
       dFreq       =  0.0;  /* frequency of the soure        [Hz]  */

long     BufferSize;       /* size of the neutron input and output buffer */
Neutron* InputNeutrons;    /* input neutron Buffer */
Neutron* OutputNeutrons;   /* output neutron buffer */

double   wei_min=0.0;      /* Minimal weight for tracing neutron */
long     keygrav=1;
short    bTrace=FALSE,     /* criterion: write trace files */
         bOldFrame=FALSE,  /* criterion: co-ordinate system of prev. module used for current module */
         bSepRate=FALSE;   /* criterion: write separate count rates */
char*    ParDirectory;     /* parameter directory */
int      ParDirectoryLength;

/* Convert McStas state parameters to VITESS Neutron structure. In
   VITESS, the neutron velocity is represented by a wavelength in
   AAngstroem and a unit direction vector, time is in msec and
   positions are in cm.*/
Neutron mcstas2vitess(double x, double y, double z,
          double vx, double vy, double vz,
          double t,
          double sx, double sy, double sz,
          double p)
{
  static unsigned long  i=0;
  double v,s;     /* Neutron speed */
  Neutron neu;      /* Vitess Neutron structure */

  neu.Position[0] = z*100;  /* Convert position from m to cm */
  neu.Position[1] = x*100;
  neu.Position[2] = y*100;
  v = sqrt(vx*vx + vy*vy + vz*vz);
  if(v == 0.0)
  {
    fprintf(stderr, "Error: zero velocity! (mcstas2vitess)\n");
    exit(1);
  }
  neu.Wavelength = 3956.0346/v; /* Convert speed to wavelength */
  neu.Vector[0] = vz/v;   /* Convert velocity to unit direction vector */
  neu.Vector[1] = vx/v;
  neu.Vector[2] = vy/v;
  s = sqrt(sx*sx+sy*sy+sz*sz);
  if(s != 0.0)
  {
    neu.Spin[0] = sz/s;
    neu.Spin[1] = sx/s;
    neu.Spin[2] = sy/s;
  }

  neu.Time = t*1000;    /* Convert time from sec to msec */
  neu.Probability = p;    /* Neutron weight */
  neu.Color = 0;
  neu.Debug = 'N';
  neu.ID.IDGrp[0] = 'A';
  neu.ID.IDGrp[1] = 'A';
  neu.ID.IDNo     = i++;
  return neu;
}

/* Convert VITESS neutron structure to McStas state parameters. In
   VITESS, the neutron velocity is represented by a wavelength in
   AAngstroem and a unit direction vector, time is in msec and
   positions are in cm. */
void vitess2mcstas(Neutron neu,
       double *x, double *y, double *z,
       double *vx, double *vy, double *vz,
       double *t,
       double *sx, double *sy, double *sz,
       double *p)
{
  double v;     /* Neutron speed */

  *x = 0.01*neu.Position[1];  /* Convert position from cm to m */
  *y = 0.01*neu.Position[2];
  *z = 0.01*neu.Position[0];
  if(neu.Wavelength == 0.0)
  {
    fprintf(stderr, "Error: zero wavelength! (mcstas2vitess: )\n");
    exit(1);
  }
  v = 3956.0346/neu.Wavelength; /* Convert wavelength to speed */
  *vx = v*neu.Vector[1];  /* Convert unit direction vector to velocity */
  *vy = v*neu.Vector[2];
  *vz = v*neu.Vector[0];
  *sx = neu.Spin[1];
  *sy = neu.Spin[2];
  *sz = neu.Spin[0];
  *t = 0.001*neu.Time;    /* Convert msec to sec */
  *p = neu.Probability;   /* Neutron weight */
}

/* Standard VITESS option parsing. */
char *vitess_infile;    /* Neutron input file name, or NULL. */
char *vitess_outfile;   /* Neutron output file name, or NULL. */
int vitess_tracepoints;   /* If true, use dots as progress-indicator */
int vitess_repcnt;    /* Number of times to repeat this neutron */
int vitess_bufsize;   /* The buffer size for neutron read/write */

void
vitess_option_error(char *opt)
{
  fprintf(stderr, "Error: Invalid VITESS option '%s'\n", opt);
  exit(1);
}

void
vitess_parseopt(int argc, char *argv[],
    double *dptr[], char dchr[], char **sptr[], char schr[])
{
  int i, j;

  /* Initialize variables to defaults. */
  vitess_infile = NULL;
  vitess_outfile = NULL;
  vitess_tracepoints = 0;
  vitess_repcnt = 1;
  vitess_bufsize = 10000;
  for(i = 0; dptr[i]; i++)
    *dptr[i] = 0;   /* Set all double parameters to zero */
  for(i = 0; sptr[i]; i++)
    *sptr[i] = NULL;    /* Set all string parameters to NULL */

  /* Now loop over all option arguments. */
  for(i = 1; i < argc; i++)
  {
    if(argv[i][0] != '-')
      vitess_option_error(argv[i]);

    if (argv[i][0] == '-' && (argv[i][1] == '-'))
    { switch(argv[i][2])
      {
        case 'f':
          vitess_infile = &argv[i][3];
          break;
        case 'F':
          vitess_outfile = &argv[i][3];
          break;
        case 'J':
          vitess_tracepoints = 1;
          break;
        case 'L':
          if(!freopen(&argv[i][3], "wt", stderr))
          {
            fprintf(stderr, "Can't open %s for output!\n", &argv[i][2]);
            exit(1);
          }
          LogFilePtr=stderr;
          break;
        case 'Z':
          srandom(atol(&(argv[i][3])));
          break;
        case 'G':
          mcgravitation = atol(&(argv[i][3]));
          keygrav = mcgravitation;
          break;
        case 'B':
          vitess_bufsize = atol(&(argv[i][3]));
          break;
        case 'P':
          setParDirectory(&argv[i][3]);
          break;
        case 'U':
          wei_min = atof(&(argv[i][3]));    /* minimal weight for tracing neutron */
      }
    }
    else 
    { /* First look for a matching double parameter. */
      for(j = 0; dchr[j]; j++)
      {
        if(argv[i][1] == dchr[j])
        {
          *dptr[j] = atof(&(argv[i][2]));
          goto end_loop;
        }
      }
      if(!dchr[j])
      {
        /* Then look for a matching string parameter. */
        for(j = 0; schr[j]; j++)
        {
          if(argv[i][1] == schr[j])
          {
            *sptr[j] = &(argv[i][2]);
            goto end_loop;
          }
        }
        if(!schr[j])
          vitess_option_error(argv[i]);
      }
    }
end_loop:
  continue;
  }
}

void WriteInstrData(long nModuleNo, VectorType Pos, double dLength, double dRotZ, double dRotY)
{
  FILE*  pFile=NULL;
  char   *pBuffer, sBuffer[CHAR_BUF_SMALL+1];
  int    m, i;

  /* source module writes header lines */
  if (nModuleNo==0)
  { pFile = fopen( FullParName("instrument.inf"), "w");
    fprintf(pFile, "# No ID    module           len [m]  x [m]   y [m]   z [m]    hor. [deg] ver.      W-Par.       H-Par.       R-Par       number  type Description\n");
    fprintf(pFile, "# ------------------------------------------------------------------------------------------------------------------------------------------------\n");
  }
  /* first module of 2nd, 3rd ... part re-writes file up to end of previous part */
  else if (vitess_infile!=NULL)
  { i=-1;
    pFile = fopen(FullParName("instrument.inf"), "r");
    pBuffer=malloc(CHAR_BUF_SMALL*(nModuleNo+3+3));
    for (m=-2; m<nModuleNo; m++)
    { fgets (sBuffer, sizeof(sBuffer)-1, pFile);
      strcpy(&pBuffer[++i*CHAR_BUF_SMALL], sBuffer);
      if (memcmp(sBuffer, "EOP", 3)==0)
      { fgets (sBuffer, sizeof(sBuffer)-1, pFile);
        strcpy(&pBuffer[++i*CHAR_BUF_SMALL], sBuffer);
      }
    }
    fclose(pFile);
    pFile = fopen( FullParName("instrument.inf"), "w");
    for (m=0; m<=i; m++)
      fprintf(pFile, "%s", &pBuffer[CHAR_BUF_SMALL*m]);
    fprintf(pFile, "EOP\n");
    free(pBuffer);
  }
  /* each other module appends a line */
  else
  { pFile = fopen(FullParName("instrument.inf"), "a");
  }

  if (pFile) {
    char cNF=' ';
    if (bOldFrame) cNF='F';
    fprintf(pFile, "%3ld %3d %-18.18s %7.3f %7.3f %7.3f %7.3f  %8.3f %8.3f  %12.4e %12.4e %12.4e %c %5i %5d\n",
                   nModuleNo, 500, mcinstrument_name, dLength, Pos[0], Pos[1], Pos[2],
                   180.0/M_PI*dRotZ, 180.0/M_PI*dRotY, 0.0, 0.0, 0.0, cNF, 0, 0);
    /* mark end of actual part */
    if (vitess_outfile!=NULL && nModuleNo > 0)
      fprintf(pFile, "EOP\n");
    fclose(pFile);
  }
}

void ReadInstrData(long* pModuleNo, VectorType Pos, double* pLength, double* pRotZ, double* pRotY)
{
  FILE*  pFile=NULL;
  int    nModuleID;
  long   No=0, nDum;
  char   sBuffer[CHAR_BUF_LENGTH]="", sLine[CHAR_BUF_LENGTH]="", sLineH[CHAR_BUF_LENGTH]="";

  *pModuleNo = 0;
  Pos[0]   = Pos[1] = Pos[2] = 0.0;
  *pLength = 0.0;
  *pRotY   = 0.0;
  *pRotZ   = 0.0;

  pFile = fopen(FullParName("instrument.inf"), "r");
  if (pFile)
  {
    if (vitess_infile==NULL)
    { /* Read last line and copy content, except:
		   lines containing F at pos 116-118, they have not a new frame) */
      while (ReadLine(pFile, sBuffer, sizeof(sBuffer)-1))
      { sscanf(sBuffer, "%ld", pModuleNo);
		  // ndig = short(floor(lg10(*pModuleNo));
        if (sBuffer[116]!='F' && sBuffer[117]!='F' && sBuffer[118]!='F') strcpy(sLine, sBuffer);
		}
    }
    else
    {  /* read until end of previous part, if input file is used */
      while (ReadLine(pFile, sBuffer, sizeof(sBuffer)-1))
      { if (memcmp(sBuffer, "EOP", 3)==0)
        {  *pModuleNo = No;
           strcpy(sLine, sLineH);
        }
        else
        { sscanf(sBuffer, "%ld", &No);
          if (sBuffer[116]!='F' && sBuffer[117]!='F' && sBuffer[118]!='F') strcpy(sLineH, sBuffer);
		  }
      }
      if (strlen(sLine)==0) {*pModuleNo = No; strcpy(sLine, sLineH);}
    }
    sscanf(sLine, "%ld %3d %18c %lf %lf %lf %lf %lf %lf",
                  &nDum, &nModuleID, sBuffer, pLength, &Pos[0], &Pos[1], &Pos[2], pRotZ, pRotY);
    *pRotZ *= M_PI/180.0;
    *pRotY *= M_PI/180.0;
    fclose(pFile);
  }
}


vitess_write(double NumNeutRead, double NumNeutWritten, double CntRate, double CntRateSqr,
             double dShiftX,     double dShiftY,        double dShiftZ, 
             double dHorizAngle, double dVertAngle)
{
  double dRotMatrix[3][3], dRotY, dRotZ, dLength,
         CntRateErr;
  long   nModuleNo=0;
  int    k;
  VectorType Shift,  /* Shift of end position        [m] */
             EndPos; /* end position of prev. module [m] */

  fprintf(LogFilePtr,"\n\nVITESS 2.6 / McStas 1.9  module %s\n", mcinstrument_name);

  /* update 'instrument.inf' */
  ReadInstrData(&nModuleNo, EndPos, &dLength, &dRotZ, &dRotY);
  nModuleNo++;
  Shift[0] = dShiftZ;
  Shift[1] = dShiftX;
  Shift[2] = dShiftY;
  FillRMatrixZY(dRotMatrix, dRotY, dRotZ);
  RotBackVector(dRotMatrix, Shift);
  for (k=0; k<3; k++)
    EndPos[k] += Shift[k];
  dLength += LengthVector(Shift);
  dRotZ   += dHorizAngle;
  dRotY   += dVertAngle;
  WriteInstrData(nModuleNo, EndPos, dLength, dRotZ, dRotY);

  CntRateErr = sqrt( sq(CntRate)/NumNeutWritten
                   + (NumNeutWritten*CntRateSqr-sq(CntRate)) / (NumNeutWritten-1) );
  fprintf(stderr, "%2ld number of trajectories read         : %11.0f\n", nModuleNo, NumNeutRead);
  fprintf(stderr, "   number of trajectories written      : %11.0f\n", NumNeutWritten);
  fprintf(stderr, "(time averaged) neutron count rate     : %11.4e +/- %10.3e n/s \n", CntRate, CntRateErr);

}


/**************************************************************/
/* Init does a general program initialization, which is ok    */
/* for all modules of the VITESS program package.             */
/**************************************************************/

void McInitVt()
{
  /* Set some default values */
  BufferSize  = 2;
  LogFilePtr  = stderr;

  if (mcdirname==NULL)
    setParDirectory(getenv("PWD") ? getenv("PWD") : ".");
  else
    setParDirectory(mcdirname);
  idum    = -mcseed;
  keygrav = (long) mcgravitation;      /* key for gravity 1 -yes (default), 0 - no  */

  /* allocte memory for the neutron buffers */
  if((InputNeutrons=(Neutron *)calloc(BufferSize, sizeof(Neutron)))==NULL) {
    fprintf(LogFilePtr, "Couldn't allocate memory for input buffer\n");
    exit(-1);
  }
  if((OutputNeutrons=(Neutron *)calloc(BufferSize, sizeof(Neutron)))==NULL) {
    fprintf(LogFilePtr, "Couldn't allocate memory for output buffer\n");
    exit(-1);
  }

  /* initalize the random number generator */
  ran3(&idum);
}


/*****************************************************************/
/* Cleanup() does last things before the VITESS module is closed */
/* e.g. buffers are flushed and files are closed etc.            */
/* you should also write your OwnCleanup() for your module       */
/*****************************************************************/
void McCleanupVt()
{
  /* release the buffer memory */
  free(InputNeutrons);
  free(OutputNeutrons);
}


void setParDirectory (char *a) {
  int len;
  if ((len = strlen(a))) {
    /* last character should be a slash */
    if (a[len-1] == cSlash) {
      memcpy ((ParDirectory = (char *) malloc(len+1)), a, len);
    } else {
      memcpy ((ParDirectory = (char *) malloc(len+2)), a, len);
      ParDirectory[len++] = cSlash;
      ParDirectory[len] = 0;
    }
    ParDirectoryLength = len;
  }
}

/* Adding path of parameter directory to file name */
char* FullParName(char* filename)
{
  int sel=0;
  char *res, *a=NULL;
  int alen, blen;

  if (filename == 0)
     return 0;

  /* Do not change an absolute path. */
 #ifdef WIN32
  /* we consider a filename with : as absolute */
  if (strstr(filename, ":")) sel = -1;
 #else
  if (filename[0] == '/') sel = -1;
 #endif
  if (sel == -1) {
    alen = 0;
  } else if (sel == 0) {
    a = ParDirectory;
    alen = ParDirectoryLength;
  }
  blen = strlen(filename);
  if ((res = (char *) malloc(alen+blen+1)))
  { if (alen)
      strcpy(res, a);
    else
      strcpy(res,"");
    strcat(res, filename);
  }
  return res;
}

/* end of vitess-lib.c */

#line 13108 "./Test_Fermi.c"

/* Instrument parameters. */
int mcipFermi;
MCNUM mciplambda;
MCNUM mcipwidth_FC;
MCNUM mcipheight_FC;
MCNUM mciplength_FC;
MCNUM mcipFC_Hz;
MCNUM mcipNslit_FC;
MCNUM mcipd_SF;
MCNUM mcipd_FD;
MCNUM mcipphase;

#define mcNUMIPAR 10
int mcnumipar = 10;
struct mcinputtable_struct mcinputtable[mcNUMIPAR+1] = {
  "Fermi", &mcipFermi, instr_type_int, "1", 
  "lambda", &mciplambda, instr_type_double, "3.39", 
  "width_FC", &mcipwidth_FC, instr_type_double, "0.044", 
  "height_FC", &mcipheight_FC, instr_type_double, "0.064", 
  "length_FC", &mciplength_FC, instr_type_double, "0.012", 
  "FC_Hz", &mcipFC_Hz, instr_type_double, "100", 
  "Nslit_FC", &mcipNslit_FC, instr_type_double, "120", 
  "d_SF", &mcipd_SF, instr_type_double, "3", 
  "d_FD", &mcipd_FD, instr_type_double, "3", 
  "phase", &mcipphase, instr_type_double, "271.92", 
  NULL, NULL, instr_type_double, ""
};

/* User declarations from instrument definition. */
#define mccompcurname  Test_Fermi
#define mccompcurtype  INSTRUMENT
#define mccompcurindex 0
#define mcposaTest_Fermi coords_set(0,0,0)
#define Fermi mcipFermi
#define lambda mciplambda
#define width_FC mcipwidth_FC
#define height_FC mcipheight_FC
#define length_FC mciplength_FC
#define FC_Hz mcipFC_Hz
#define Nslit_FC mcipNslit_FC
#define d_SF mcipd_SF
#define d_FD mcipd_FD
#define phase mcipphase
#line 44 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  double time_to_arrival;
  double time_window_width;
#line 13156 "./Test_Fermi.c"
#undef phase
#undef d_FD
#undef d_SF
#undef Nslit_FC
#undef FC_Hz
#undef length_FC
#undef height_FC
#undef width_FC
#undef lambda
#undef Fermi
#undef mcposaTest_Fermi
#undef mccompcurindex
#undef mccompcurtype
#undef mccompcurname

/* neutron state table at each component input (local coords) */
/* [x, y, z, vx, vy, vz, t, sx, sy, sz, p] */
MCNUM mccomp_storein[11*13];
/* Components position table (absolute and relative coords) */
Coords mccomp_posa[13];
Coords mccomp_posr[13];
/* Counter for each comp to check for inactive ones */
MCNUM  mcNCounter[13];
MCNUM  mcPCounter[13];
MCNUM  mcP2Counter[13];
#define mcNUMCOMP 12 /* number of components */
/* Counter for PROP ABSORB */
MCNUM  mcAbsorbProp[13];
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
char mccSource_flux_file[16384];
char mccSource_xdiv_file[16384];
char mccSource_ydiv_file[16384];
MCNUM mccSource_radius;
MCNUM mccSource_dist;
MCNUM mccSource_focus_xw;
MCNUM mccSource_focus_yh;
MCNUM mccSource_focus_aw;
MCNUM mccSource_focus_ah;
MCNUM mccSource_E0;
MCNUM mccSource_dE;
MCNUM mccSource_lambda0;
MCNUM mccSource_dlambda;
MCNUM mccSource_I1;
MCNUM mccSource_yheight;
MCNUM mccSource_xwidth;
MCNUM mccSource_verbose;
MCNUM mccSource_T1;
MCNUM mccSource_flux_file_perAA;
MCNUM mccSource_flux_file_log;
MCNUM mccSource_Lmin;
MCNUM mccSource_Lmax;
MCNUM mccSource_Emin;
MCNUM mccSource_Emax;
MCNUM mccSource_T2;
MCNUM mccSource_I2;
MCNUM mccSource_T3;
MCNUM mccSource_I3;
MCNUM mccSource_zdepth;
int mccSource_target_index;

/* Definition parameters for component 'Monitor1_xt' [3]. */
#define mccMonitor1_xt_user1 FLT_MAX
#define mccMonitor1_xt_user2 FLT_MAX
#define mccMonitor1_xt_user3 FLT_MAX
/* Setting parameters for component 'Monitor1_xt' [3]. */
MCNUM mccMonitor1_xt_xwidth;
MCNUM mccMonitor1_xt_yheight;
MCNUM mccMonitor1_xt_zdepth;
MCNUM mccMonitor1_xt_xmin;
MCNUM mccMonitor1_xt_xmax;
MCNUM mccMonitor1_xt_ymin;
MCNUM mccMonitor1_xt_ymax;
MCNUM mccMonitor1_xt_zmin;
MCNUM mccMonitor1_xt_zmax;
MCNUM mccMonitor1_xt_bins;
MCNUM mccMonitor1_xt_min;
MCNUM mccMonitor1_xt_max;
MCNUM mccMonitor1_xt_restore_neutron;
MCNUM mccMonitor1_xt_radius;
char mccMonitor1_xt_options[16384];
char mccMonitor1_xt_filename[16384];
char mccMonitor1_xt_geometry[16384];
char mccMonitor1_xt_username1[16384];
char mccMonitor1_xt_username2[16384];
char mccMonitor1_xt_username3[16384];
int mccMonitor1_xt_nowritefile;

/* Setting parameters for component 'FC_GuideG' [5]. */
MCNUM mccFC_GuideG_w1;
MCNUM mccFC_GuideG_h1;
MCNUM mccFC_GuideG_w2;
MCNUM mccFC_GuideG_h2;
MCNUM mccFC_GuideG_l;
MCNUM mccFC_GuideG_R0;
MCNUM mccFC_GuideG_Qc;
MCNUM mccFC_GuideG_alpha;
MCNUM mccFC_GuideG_m;
MCNUM mccFC_GuideG_W;
MCNUM mccFC_GuideG_nslit;
MCNUM mccFC_GuideG_d;
MCNUM mccFC_GuideG_mleft;
MCNUM mccFC_GuideG_mright;
MCNUM mccFC_GuideG_mtop;
MCNUM mccFC_GuideG_mbottom;
MCNUM mccFC_GuideG_nhslit;
MCNUM mccFC_GuideG_G;
MCNUM mccFC_GuideG_aleft;
MCNUM mccFC_GuideG_aright;
MCNUM mccFC_GuideG_atop;
MCNUM mccFC_GuideG_abottom;
MCNUM mccFC_GuideG_wavy;
MCNUM mccFC_GuideG_wavy_z;
MCNUM mccFC_GuideG_wavy_tb;
MCNUM mccFC_GuideG_wavy_lr;
MCNUM mccFC_GuideG_chamfers;
MCNUM mccFC_GuideG_chamfers_z;
MCNUM mccFC_GuideG_chamfers_lr;
MCNUM mccFC_GuideG_chamfers_tb;
MCNUM mccFC_GuideG_nelements;
MCNUM mccFC_GuideG_nu;
MCNUM mccFC_GuideG_phase;
char mccFC_GuideG_reflect[16384];

/* Setting parameters for component 'FC_GuideC' [6]. */
MCNUM mccFC_GuideC_w1;
MCNUM mccFC_GuideC_h1;
MCNUM mccFC_GuideC_w2;
MCNUM mccFC_GuideC_h2;
MCNUM mccFC_GuideC_l;
MCNUM mccFC_GuideC_R0;
MCNUM mccFC_GuideC_Qc;
MCNUM mccFC_GuideC_alpha;
MCNUM mccFC_GuideC_m;
MCNUM mccFC_GuideC_nslit;
MCNUM mccFC_GuideC_d;
MCNUM mccFC_GuideC_Qcx;
MCNUM mccFC_GuideC_Qcy;
MCNUM mccFC_GuideC_alphax;
MCNUM mccFC_GuideC_alphay;
MCNUM mccFC_GuideC_W;
MCNUM mccFC_GuideC_mx;
MCNUM mccFC_GuideC_my;
MCNUM mccFC_GuideC_nu;
MCNUM mccFC_GuideC_phase;

/* Setting parameters for component 'FC_McStas' [7]. */
MCNUM mccFC_McStas_phase;
MCNUM mccFC_McStas_radius;
MCNUM mccFC_McStas_nu;
MCNUM mccFC_McStas_w;
MCNUM mccFC_McStas_nslit;
MCNUM mccFC_McStas_R0;
MCNUM mccFC_McStas_Qc;
MCNUM mccFC_McStas_alpha;
MCNUM mccFC_McStas_m;
MCNUM mccFC_McStas_W;
MCNUM mccFC_McStas_length;
MCNUM mccFC_McStas_eff;
MCNUM mccFC_McStas_zero_time;
MCNUM mccFC_McStas_xwidth;
MCNUM mccFC_McStas_verbose;
MCNUM mccFC_McStas_yheight;
MCNUM mccFC_McStas_curvature;
MCNUM mccFC_McStas_delay;

/* Setting parameters for component 'FC_ILL' [8]. */
MCNUM mccFC_ILL_phase;
MCNUM mccFC_ILL_radius;
MCNUM mccFC_ILL_nu;
MCNUM mccFC_ILL_yheight;
MCNUM mccFC_ILL_w;
MCNUM mccFC_ILL_nslit;
MCNUM mccFC_ILL_R0;
MCNUM mccFC_ILL_Qc;
MCNUM mccFC_ILL_alpha;
MCNUM mccFC_ILL_m;
MCNUM mccFC_ILL_W;
MCNUM mccFC_ILL_length;
MCNUM mccFC_ILL_eff;
MCNUM mccFC_ILL_zero_time;
MCNUM mccFC_ILL_xwidth;
MCNUM mccFC_ILL_verbose;

/* Setting parameters for component 'FC_Vitess' [10]. */
char mccFC_Vitess_sGeomFileName[16384];
int mccFC_Vitess_GeomOption;
int mccFC_Vitess_zerotime;
int mccFC_Vitess_Nchannels;
int mccFC_Vitess_Ngates;
MCNUM mccFC_Vitess_freq;
MCNUM mccFC_Vitess_height;
MCNUM mccFC_Vitess_width;
MCNUM mccFC_Vitess_depth;
MCNUM mccFC_Vitess_r_curv;
MCNUM mccFC_Vitess_diameter;
MCNUM mccFC_Vitess_Phase;
MCNUM mccFC_Vitess_wallwidth;

/* Definition parameters for component 'Monitor2_xt' [11]. */
#define mccMonitor2_xt_user1 FLT_MAX
#define mccMonitor2_xt_user2 FLT_MAX
#define mccMonitor2_xt_user3 FLT_MAX
/* Setting parameters for component 'Monitor2_xt' [11]. */
MCNUM mccMonitor2_xt_xwidth;
MCNUM mccMonitor2_xt_yheight;
MCNUM mccMonitor2_xt_zdepth;
MCNUM mccMonitor2_xt_xmin;
MCNUM mccMonitor2_xt_xmax;
MCNUM mccMonitor2_xt_ymin;
MCNUM mccMonitor2_xt_ymax;
MCNUM mccMonitor2_xt_zmin;
MCNUM mccMonitor2_xt_zmax;
MCNUM mccMonitor2_xt_bins;
MCNUM mccMonitor2_xt_min;
MCNUM mccMonitor2_xt_max;
MCNUM mccMonitor2_xt_restore_neutron;
MCNUM mccMonitor2_xt_radius;
char mccMonitor2_xt_options[16384];
char mccMonitor2_xt_filename[16384];
char mccMonitor2_xt_geometry[16384];
char mccMonitor2_xt_username1[16384];
char mccMonitor2_xt_username2[16384];
char mccMonitor2_xt_username3[16384];
int mccMonitor2_xt_nowritefile;

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
#line 13419 "./Test_Fermi.c"
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
#define mccompcurtype  Source_gen
#define mccompcurindex 2
#define p_in mccSource_p_in
#define lambda1 mccSource_lambda1
#define lambda2 mccSource_lambda2
#define lambda3 mccSource_lambda3
#define pTable mccSource_pTable
#define pTable_x mccSource_pTable_x
#define pTable_y mccSource_pTable_y
#define pTable_xmin mccSource_pTable_xmin
#define pTable_xmax mccSource_pTable_xmax
#define pTable_xsum mccSource_pTable_xsum
#define pTable_ymin mccSource_pTable_ymin
#define pTable_ymax mccSource_pTable_ymax
#define pTable_ysum mccSource_pTable_ysum
#define pTable_dxmin mccSource_pTable_dxmin
#define pTable_dxmax mccSource_pTable_dxmax
#define pTable_dymin mccSource_pTable_dymin
#define pTable_dymax mccSource_pTable_dymax
#define flux_file mccSource_flux_file
#define xdiv_file mccSource_xdiv_file
#define ydiv_file mccSource_ydiv_file
#define radius mccSource_radius
#define dist mccSource_dist
#define focus_xw mccSource_focus_xw
#define focus_yh mccSource_focus_yh
#define focus_aw mccSource_focus_aw
#define focus_ah mccSource_focus_ah
#define E0 mccSource_E0
#define dE mccSource_dE
#define lambda0 mccSource_lambda0
#define dlambda mccSource_dlambda
#define I1 mccSource_I1
#define yheight mccSource_yheight
#define xwidth mccSource_xwidth
#define verbose mccSource_verbose
#define T1 mccSource_T1
#define flux_file_perAA mccSource_flux_file_perAA
#define flux_file_log mccSource_flux_file_log
#define Lmin mccSource_Lmin
#define Lmax mccSource_Lmax
#define Emin mccSource_Emin
#define Emax mccSource_Emax
#define T2 mccSource_T2
#define I2 mccSource_I2
#define T3 mccSource_T3
#define I3 mccSource_I3
#define zdepth mccSource_zdepth
#define target_index mccSource_target_index
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

#line 13503 "./Test_Fermi.c"
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

/* User declarations for component 'Monitor1_xt' [3]. */
#define mccompcurname  Monitor1_xt
#define mccompcurtype  Monitor_nD
#define mccompcurindex 3
#define user1 mccMonitor1_xt_user1
#define user2 mccMonitor1_xt_user2
#define user3 mccMonitor1_xt_user3
#define DEFS mccMonitor1_xt_DEFS
#define Vars mccMonitor1_xt_Vars
#define detector mccMonitor1_xt_detector
#define offdata mccMonitor1_xt_offdata
#define xwidth mccMonitor1_xt_xwidth
#define yheight mccMonitor1_xt_yheight
#define zdepth mccMonitor1_xt_zdepth
#define xmin mccMonitor1_xt_xmin
#define xmax mccMonitor1_xt_xmax
#define ymin mccMonitor1_xt_ymin
#define ymax mccMonitor1_xt_ymax
#define zmin mccMonitor1_xt_zmin
#define zmax mccMonitor1_xt_zmax
#define bins mccMonitor1_xt_bins
#define min mccMonitor1_xt_min
#define max mccMonitor1_xt_max
#define restore_neutron mccMonitor1_xt_restore_neutron
#define radius mccMonitor1_xt_radius
#define options mccMonitor1_xt_options
#define filename mccMonitor1_xt_filename
#define geometry mccMonitor1_xt_geometry
#define username1 mccMonitor1_xt_username1
#define username2 mccMonitor1_xt_username2
#define username3 mccMonitor1_xt_username3
#define nowritefile mccMonitor1_xt_nowritefile
#line 224 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/Monitor_nD.comp"
  MonitornD_Defines_type DEFS;
  MonitornD_Variables_type Vars;
  MCDETECTOR detector;
  off_struct offdata;
#line 13592 "./Test_Fermi.c"
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

/* User declarations for component 'FC_Position' [4]. */
#define mccompcurname  FC_Position
#define mccompcurtype  Arm
#define mccompcurindex 4
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'FC_GuideG' [5]. */
#define mccompcurname  FC_GuideG
#define mccompcurtype  Guide_gravity
#define mccompcurindex 5
#define GVars mccFC_GuideG_GVars
#define pTable mccFC_GuideG_pTable
#define w1 mccFC_GuideG_w1
#define h1 mccFC_GuideG_h1
#define w2 mccFC_GuideG_w2
#define h2 mccFC_GuideG_h2
#define l mccFC_GuideG_l
#define R0 mccFC_GuideG_R0
#define Qc mccFC_GuideG_Qc
#define alpha mccFC_GuideG_alpha
#define m mccFC_GuideG_m
#define W mccFC_GuideG_W
#define nslit mccFC_GuideG_nslit
#define d mccFC_GuideG_d
#define mleft mccFC_GuideG_mleft
#define mright mccFC_GuideG_mright
#define mtop mccFC_GuideG_mtop
#define mbottom mccFC_GuideG_mbottom
#define nhslit mccFC_GuideG_nhslit
#define G mccFC_GuideG_G
#define aleft mccFC_GuideG_aleft
#define aright mccFC_GuideG_aright
#define atop mccFC_GuideG_atop
#define abottom mccFC_GuideG_abottom
#define wavy mccFC_GuideG_wavy
#define wavy_z mccFC_GuideG_wavy_z
#define wavy_tb mccFC_GuideG_wavy_tb
#define wavy_lr mccFC_GuideG_wavy_lr
#define chamfers mccFC_GuideG_chamfers
#define chamfers_z mccFC_GuideG_chamfers_z
#define chamfers_lr mccFC_GuideG_chamfers_lr
#define chamfers_tb mccFC_GuideG_chamfers_tb
#define nelements mccFC_GuideG_nelements
#define nu mccFC_GuideG_nu
#define phase mccFC_GuideG_phase
#define reflect mccFC_GuideG_reflect
#line 334 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide_gravity.comp"
  Gravity_guide_Vars_type GVars;
  t_Table pTable;
#line 13676 "./Test_Fermi.c"
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

/* User declarations for component 'FC_GuideC' [6]. */
#define mccompcurname  FC_GuideC
#define mccompcurtype  Guide_channeled
#define mccompcurindex 6
#define w1c mccFC_GuideC_w1c
#define w2c mccFC_GuideC_w2c
#define ww mccFC_GuideC_ww
#define hh mccFC_GuideC_hh
#define whalf mccFC_GuideC_whalf
#define hhalf mccFC_GuideC_hhalf
#define lwhalf mccFC_GuideC_lwhalf
#define lhhalf mccFC_GuideC_lhhalf
#define w1 mccFC_GuideC_w1
#define h1 mccFC_GuideC_h1
#define w2 mccFC_GuideC_w2
#define h2 mccFC_GuideC_h2
#define l mccFC_GuideC_l
#define R0 mccFC_GuideC_R0
#define Qc mccFC_GuideC_Qc
#define alpha mccFC_GuideC_alpha
#define m mccFC_GuideC_m
#define nslit mccFC_GuideC_nslit
#define d mccFC_GuideC_d
#define Qcx mccFC_GuideC_Qcx
#define Qcy mccFC_GuideC_Qcy
#define alphax mccFC_GuideC_alphax
#define alphay mccFC_GuideC_alphay
#define W mccFC_GuideC_W
#define mx mccFC_GuideC_mx
#define my mccFC_GuideC_my
#define nu mccFC_GuideC_nu
#define phase mccFC_GuideC_phase
#line 81 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide_channeled.comp"
double w1c;
double w2c;
double ww, hh;
double whalf, hhalf;
double lwhalf, lhhalf;
#line 13755 "./Test_Fermi.c"
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

/* User declarations for component 'FC_McStas' [7]. */
#define mccompcurname  FC_McStas
#define mccompcurtype  FermiChopper
#define mccompcurindex 7
#define FCVars mccFC_McStas_FCVars
#define phase mccFC_McStas_phase
#define radius mccFC_McStas_radius
#define nu mccFC_McStas_nu
#define w mccFC_McStas_w
#define nslit mccFC_McStas_nslit
#define R0 mccFC_McStas_R0
#define Qc mccFC_McStas_Qc
#define alpha mccFC_McStas_alpha
#define m mccFC_McStas_m
#define W mccFC_McStas_W
#define length mccFC_McStas_length
#define eff mccFC_McStas_eff
#define zero_time mccFC_McStas_zero_time
#define xwidth mccFC_McStas_xwidth
#define verbose mccFC_McStas_verbose
#define yheight mccFC_McStas_yheight
#define curvature mccFC_McStas_curvature
#define delay mccFC_McStas_delay
#line 278 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/FermiChopper.comp"
        struct FermiChopper_struct FCVars;
#line 13813 "./Test_Fermi.c"
#undef delay
#undef curvature
#undef yheight
#undef verbose
#undef xwidth
#undef zero_time
#undef eff
#undef length
#undef W
#undef m
#undef alpha
#undef Qc
#undef R0
#undef nslit
#undef w
#undef nu
#undef radius
#undef phase
#undef FCVars
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'FC_ILL' [8]. */
#define mccompcurname  FC_ILL
#define mccompcurtype  FermiChopper_ILL
#define mccompcurindex 8
#define FCVars mccFC_ILL_FCVars
#define phase mccFC_ILL_phase
#define radius mccFC_ILL_radius
#define nu mccFC_ILL_nu
#define yheight mccFC_ILL_yheight
#define w mccFC_ILL_w
#define nslit mccFC_ILL_nslit
#define R0 mccFC_ILL_R0
#define Qc mccFC_ILL_Qc
#define alpha mccFC_ILL_alpha
#define m mccFC_ILL_m
#define W mccFC_ILL_W
#define length mccFC_ILL_length
#define eff mccFC_ILL_eff
#define zero_time mccFC_ILL_zero_time
#define xwidth mccFC_ILL_xwidth
#define verbose mccFC_ILL_verbose
#line 237 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../contrib/FermiChopper_ILL.comp"
  struct FermiChopper_ILL_struct FCVars;
#line 13860 "./Test_Fermi.c"
#undef verbose
#undef xwidth
#undef zero_time
#undef eff
#undef length
#undef W
#undef m
#undef alpha
#undef Qc
#undef R0
#undef nslit
#undef w
#undef yheight
#undef nu
#undef radius
#undef phase
#undef FCVars
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'Fake_Origin' [9]. */
#define mccompcurname  Fake_Origin
#define mccompcurtype  Arm
#define mccompcurindex 9
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'FC_Vitess' [10]. */
#define mccompcurname  FC_Vitess
#define mccompcurtype  Vitess_ChopperFermi
#define mccompcurindex 10
#define Option mccFC_Vitess_Option
#define CurvGeomOption mccFC_Vitess_CurvGeomOption
#define TOF mccFC_Vitess_TOF
#define WL mccFC_Vitess_WL
#define radius_of_curv mccFC_Vitess_radius_of_curv
#define main_depth mccFC_Vitess_main_depth
#define shift_y mccFC_Vitess_shift_y
#define angle_channel mccFC_Vitess_angle_channel
#define phase0 mccFC_Vitess_phase0
#define y_ch mccFC_Vitess_y_ch
#define x_ch mccFC_Vitess_x_ch
#define coef_pi mccFC_Vitess_coef_pi
#define XFILEName mccFC_Vitess_XFILEName
#define GeomFilePtr mccFC_Vitess_GeomFilePtr
#define Pos mccFC_Vitess_Pos
#define Dir mccFC_Vitess_Dir
#define Neutrons mccFC_Vitess_Neutrons
#define pos_ch mccFC_Vitess_pos_ch
#define omega mccFC_Vitess_omega
#define optimal_wl mccFC_Vitess_optimal_wl
#define sGeomFileName mccFC_Vitess_sGeomFileName
#define GeomOption mccFC_Vitess_GeomOption
#define zerotime mccFC_Vitess_zerotime
#define Nchannels mccFC_Vitess_Nchannels
#define Ngates mccFC_Vitess_Ngates
#define freq mccFC_Vitess_freq
#define height mccFC_Vitess_height
#define width mccFC_Vitess_width
#define depth mccFC_Vitess_depth
#define r_curv mccFC_Vitess_r_curv
#define diameter mccFC_Vitess_diameter
#define Phase mccFC_Vitess_Phase
#define wallwidth mccFC_Vitess_wallwidth
#line 137 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Vitess_ChopperFermi.comp"

#define	STRING_BUFFER 100

/* other global parameters for chopper_fermi */
int        Option,            /* 1: straight FC,  2: curved FC                       */
           CurvGeomOption;    /* 1: ideal shape (nearly parabolic)  2: circular      */
double     TOF,               /* TOF of neutron under consideration                  */
           WL,                /* wavelength of neutron                               */
           radius_of_curv,    /* radius of curvature (curved FC)                     */
           main_depth,        /* max. channel length due to diameter and total_width */
           shift_y=0.,        /* shift to channel actually written to geometry file  */ 
           angle_channel,     /* half of the curvature of a curved Fermi chopper     */
           phase0,            /* chopper phase at TOF of neutron to chopper centre   */
           y_ch[10][2000],    /* position of gates perpendicular to flight direction */
           x_ch[10][2000],    /* position of gates along flight direction            */
           coef_pi;           /* number of half-rotation to reach identical state    */
char*      GeomFileName=sGeomFileName; /* pointer to geometry file name */
FILE*      GeomFilePtr=NULL;  /* pointer to geometry file        */
VectorType Pos,               /* position of neutron             */
           Dir;               /* flight direction of neutron     */
Neutron    Neutrons;          /* parameter set at end of chopper */
   /* instance variables */

double omega;      /* rotation frequency */
double optimal_wl; /* optimal wavelength */
VectorType pos_ch = {0.0, 0.0, 0.0};; /* centre pos. of the FC in the frame of the exit of the prev. comp. [cm] */

#define MCSTAS_SHARE
/********************************************************************************************************************************************************/
/*  VITESS module 'chopper_fermi'                                                            
/*                                                                                           
/* The free non-commercial use of these routines is granted providing due credit is given to 
/* the authors.                                                                              
/*                                                                                           
/* 1.00  Jul 2002  G. Zsigmond	initial version                                             
/* 1.01  Aug 2002  G. Zsigmond	forward all coordinates                                     
/* 1.02  Sep 2002  G. Zsigmond	included more channel windows                               
/* 1.03  Apr 2003  G. Zsigmond	sign correction                                             
/* 1.04  May 2003  G. Zsigmond	info changed                                                
/* 1.05  Jun 2003  G. Zsigmond	modulo function included for safety                         
/* 1.06  Jul 2003  G. Zsigmond	generalised for optional number of pulses; warnings included
/* 1.07  Oct 2003  G. Zsigmond	superfluous modulo function cancelled                       
/* 1.08  Nov 2003  G. Zsigmond	put 2 more windows representing channels, now 6 windows    
/*                               in the big IF loop ">=" changed to ">"                      
/* 1.09  Jan 2004  K. Lieutenant changes for 'instrument.dat'                                
/* 1.10  Jan 2004  G. Zsigmond   back to 4 windows representing channels
/* 1.11  Apr 2004  G. Zsigmond   negative time of flight defined
/* 1.12  Apr 2004  G. Zsigmond   negative time of flight - corrections, set zero time
/* 1.13  May 2004  G. Zsigmond   circular geom option and channel length included in curved fc
/* 1.14  JUL 2004  G. Zsigmond   small change in Init to adapt to new GUI
/* 1.15  OCT 2004  G. Zsigmond   changed to use both even or odd number of channels
/* 1.16  MAY 2005  G. Zsigmond   output changed to give trajectory coordinates at a plane crossing the center of the chopper (to be compatible with zero time option)
/*	                              zero time option fixed to get one peak 
/*                               shadowing cylinder opening activated 
/* 1.17  MAY 2005  G. Zsigmond  new option choice of 4, 6(better,slower) or 8(much better, very slow) gates, 4 gates option adjusted
/* 1.18  MAY 2005  K. Lieutenant changes to use this module in McStas as well               
/********************************************************************************************************************************************************/

#if defined VITESS

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "general.h"
#include "init.h"
#include "intersection.h"
#include "softabort.h"

/* START HEADER STORY */

/* input parameters */
int        Ngates=4,    /* Number of gates forming the channel: 4 (default), 6 or 8 */
           zerotime=0;     /* option: set time (close to) zero                    -z */
long       Nchannels;      /* number of channels of the Fermi chopper             -l */
double     omega,          /* frequency of rotation                       [1/s]   -n */
           height,         /* height of the Fermi chopper                  [cm]   -a */
           width,          /* width of the Fermi chopper                   [cm]   -b */
           depth,          /* channel length of the Fermi chopper          [cm]   -c */
           optimal_wl,     /* wavelength of highest transmission          [Ang]   -L */
           diameter,       /* diameter of the shadowing cylinder           [cm]   -r */
           Phase,          /* dephasing angle at zero time                [deg]   -l */
           wallwidth;      /* thickness of walls separating the channels          -m */
char       sGeomFileName[512];/* name of output file for geometry information        -G */
VectorType pos_ch;         /* centre position of the Fermi chopper        [cm] X,Y,V */

/* other global parameters */
#include "chopper_fermi.h"

/* prototypes and small functions */
void		OutputTransformations(double *tof, double *wl, double *prob, VectorType Pos, VectorType Dir, VectorType SpinVector);
void		ReadParameterFile();
void		ChopperFermiInit(int argc, char *argv[]);
void		ChopperFermiCleanup();
double	asinplus(double val);
double	asinminus(double val);

double		asin2PI(double val)
{ double result;
  if (val>=0) result = (double)asin(val);
  else  result =  2*M_PI + (double)asin(val);
  return result;
}

/* FINISH HEADER STORY */


int main(int argc, char **argv)
{
  long	 i;                 /* index of trajectory      */

  /* Initialize the program according to the parameters given  */

  Init(argc, argv, VT_CHOP_FERMI);
  print_module_name("Fermi-Chopper 1.17m");
  ChopperFermiInit(argc, argv);


  /* Get the neutrons from the file */
  DECLARE_ABORT;

  while((ReadNeutrons())!= 0)
    {
      CHECK;	/* here is what happens to the neutron */

      for(i=0;i<NumNeutGot;i++)

	{
	  CHECK;
#endif

#if defined VITESS || defined MCSTAS_TRACE
	  TOF = InputNeutrons[i].Time;

	  WL = InputNeutrons[i].Wavelength;

	  CopyVector(InputNeutrons[i].Position, Pos); 

	  CopyVector(InputNeutrons[i].Vector, Dir); 

	  Dir[0]	= (double)sqrt(1 - sq(Dir[1]) - sq(Dir[2])); 


	  /* shift to center of Fermi-Chopper */
#ifdef VITESS
	  SubVector(Pos, pos_ch);
#endif

		/*trajectories which do not intersect the entrance and exit window */

	  {double pos[3], n[3];

		  n[1]=n[2]=0.; n[0]=1.;
		  if((PlaneLineIntersect(Pos, Dir, n, - diameter/2., pos))==1)
			{
				  if((pos[2]>=height/2.)||(pos[2]<= - height/2.)||(pos[1]>=diameter/2.)||(pos[1]<= - diameter/2.)) ABSORB;
			}
		  else ABSORB;

		  if((PlaneLineIntersect(Pos, Dir, n, diameter/2., pos))==1)
			{
				  if((pos[2]>=height/2.)||(pos[2]<= - height/2.)||(pos[1]>=diameter/2.)||(pos[1]<= - diameter/2.)) ABSORB;
			}
		  else ABSORB;
	

		/* translates neutron variables for X'= - diameter/2.  */

		{
			VectorType Path;

			TOF = TOF + (- diameter/2. - Pos[0]) / fabs(Dir[0]) / V_FROM_LAMBDA(WL); 
			
			if((TOF<0)&&(Nchannels==1)){fprintf(LogFilePtr,"\nERROR: Single-slit Fermi chopper needs positive flight time at the chopper position! \n"); exit(-1); }
				
			CopyVector(Dir, Path);

			MultiplyByScalar(Path, (- diameter/2. - Pos[0])/ Dir[0] );

			AddVector(Pos, Path);  
		}							/*	 Path = displacement vector */

	
	/* calculate time entering-edge and exiting-edge of 4 windows along the channels */

	  {	long j, k, m;
	  double sq_D_0_1, sq_term, omega_fact, dirpos, vz_pos, phase[10][2000];
 
	  sq_D_0_1 = sq(Dir[0]) + sq(Dir[1]);
	  dirpos = Dir[0]*Pos[1] - Dir[1]*Pos[0];
	  sq_term  = sq(dirpos) / sq_D_0_1;
	  omega_fact = omega / (V_FROM_LAMBDA(WL) * Dir[0]);
	  vz_pos = Pos[1] > 0.0 ? 1.0 : -1.0;

	  phase0 = fmod(Phase + omega*TOF, coef_pi*M_PI); 
	
	  for(k=0; k<Ngates; k++) 
	  {
		for(j=0; j < 2*Nchannels+2; j++) {

		  double x_ch_k_j, y_ch_k_j, sq_x_ch_k_j, Denom_k, Arg_k, arg_k, pha_k_j, y_ch_new_k_j;

		  x_ch_k_j    = x_ch[k][j];
		  y_ch_k_j    = y_ch[k][j];
		  sq_x_ch_k_j = sq(x_ch_k_j);
		  Denom_k     = sqrt( sq_D_0_1 * (sq_x_ch_k_j + sq(y_ch_k_j)) );

		  Arg_k = dirpos / Denom_k;

		  if (fabs(Arg_k) > 1.) {
			Arg_k = vz_pos;
			y_ch_new_k_j = Arg_k * sqrt(sq_term - sq_x_ch_k_j);
		  } else
			y_ch_new_k_j = y_ch_k_j; /* no intersection with trajectory */

		  Denom_k = sqrt( sq_D_0_1 * (sq_x_ch_k_j + sq(y_ch_new_k_j)) );

		  arg_k = (Dir[0]*y_ch_new_k_j - Dir[1]*x_ch_k_j) / Denom_k;

		  if (fabs(arg_k) > 1.) {
			phase[k][j] = 777;} 
			  
		  else {
			pha_k_j = asin(Arg_k) - asin(arg_k); 

			if(x_ch_k_j < 0.) pha_k_j = - pha_k_j; 
								
			phase[k][j] =  pha_k_j - omega_fact * (x_ch_k_j * cos(pha_k_j) - y_ch_new_k_j * sin(pha_k_j) - Pos[0]);
		  }
		}
	  }

	  if(Ngates==4){
		  m = -1; 

		  for(j=0;j<2*Nchannels+1; j++) 
		  {
					if((m == 1)&&(phase[0][j] < phase0 )&&(phase0 < phase[0][j+1])
							   &&(phase[1][j] < phase0 )&&(phase0 < phase[1][j+1])
							   &&(phase[2][j] > phase0 )&&(phase0 > phase[2][j+1])
							   &&(phase[3][j] > phase0 )&&(phase0 > phase[3][j+1]))
					{/*fprintf(LogFilePtr, "j  %d   phases %f   %f    %f\n", j, 57.296*phase[0][j], 57.296*phase[0][j+1], 57.296*phase0);*/ 
											 goto happyend;}

					m = m * (-1); 
		  }
		  
		  
		  /* also tries one turn earlier  */
		  
		  if((phase0 > 0)&&(omega > 0)) phase0 +=  - coef_pi*M_PI;
		  if((phase0 < 0)&&(omega < 0)) phase0 +=  coef_pi*M_PI;

		  m = -1; 

		  for(j=0;j<2*Nchannels+1; j++) 
		  {
					if((m == 1)&&(phase[0][j] < phase0 )&&(phase0 < phase[0][j+1])
							   &&(phase[1][j] < phase0 )&&(phase0 < phase[1][j+1])
							   &&(phase[2][j] > phase0 )&&(phase0 > phase[2][j+1])
							   &&(phase[3][j] > phase0 )&&(phase0 > phase[3][j+1]))
					{/*fprintf(LogFilePtr, "j  %d   phases %f   %f    %f\n", j, 57.296*phase[0][j], 57.296*phase[0][j+1], 57.296*phase0);*/ 
											 goto happyend;}

					m = m * (-1); 
		  }
	  }

	  if(Ngates==6){
		  m = -1; 

		  for(j=0;j<2*Nchannels+1; j++) 
		  {
					if((m == 1)&&(phase[0][j] < phase0 )&&(phase0 < phase[0][j+1])
							   &&(phase[1][j] < phase0 )&&(phase0 < phase[1][j+1])
							   &&(phase[2][j] < phase0 )&&(phase0 < phase[2][j+1])
							   &&(phase[3][j] > phase0 )&&(phase0 > phase[3][j+1])
							   &&(phase[4][j] > phase0 )&&(phase0 > phase[4][j+1])
							   &&(phase[5][j] > phase0 )&&(phase0 > phase[5][j+1]))
					{goto happyend;}
											 
					m = m * (-1); 
		  }
		  
		  /* also tries one turn earlier  */
		  
		  if((phase0 > 0)&&(omega > 0)) phase0 +=  - coef_pi*M_PI;
		  if((phase0 < 0)&&(omega < 0)) phase0 +=  coef_pi*M_PI;

		  m = -1; 

		  for(j=0;j<2*Nchannels+1; j++) 
		  {
					if((m == 1)&&(phase[0][j] < phase0 )&&(phase0 < phase[0][j+1])
							   &&(phase[1][j] < phase0 )&&(phase0 < phase[1][j+1])
							   &&(phase[2][j] < phase0 )&&(phase0 < phase[2][j+1])
							   &&(phase[3][j] > phase0 )&&(phase0 > phase[3][j+1])
							   &&(phase[4][j] > phase0 )&&(phase0 > phase[4][j+1])
							   &&(phase[5][j] > phase0 )&&(phase0 > phase[5][j+1]))
					{goto happyend;}
											 
					m = m * (-1); 
		  }
	  }

	  if(Ngates==8){
		  m = -1; 

		  for(j=0;j<2*Nchannels+1; j++) 
		  {
					if((m == 1)&&(phase[0][j] < phase0 )&&(phase0 < phase[0][j+1])
							   &&(phase[1][j] < phase0 )&&(phase0 < phase[1][j+1])
							   &&(phase[2][j] < phase0 )&&(phase0 < phase[2][j+1])
							   &&(phase[3][j] < phase0 )&&(phase0 < phase[3][j+1])
							   &&(phase[4][j] > phase0 )&&(phase0 > phase[4][j+1])
							   &&(phase[5][j] > phase0 )&&(phase0 > phase[5][j+1])
							   &&(phase[6][j] > phase0 )&&(phase0 > phase[6][j+1])
							   &&(phase[7][j] > phase0 )&&(phase0 > phase[7][j+1]))
					{goto happyend;}
											 
					m = m * (-1); 
		  }
		  
		  /* also tries one turn earlier  */
		  
		  if((phase0 > 0)&&(omega > 0)) phase0 +=  - coef_pi*M_PI;
		  if((phase0 < 0)&&(omega < 0)) phase0 +=  coef_pi*M_PI;

		  m = -1; 

		  for(j=0;j<2*Nchannels+1; j++) 
		  {
					if((m == 1)&&(phase[0][j] < phase0 )&&(phase0 < phase[0][j+1])
							   &&(phase[1][j] < phase0 )&&(phase0 < phase[1][j+1])
							   &&(phase[2][j] < phase0 )&&(phase0 < phase[2][j+1])
							   &&(phase[3][j] < phase0 )&&(phase0 < phase[3][j+1])
							   &&(phase[4][j] > phase0 )&&(phase0 > phase[4][j+1])
							   &&(phase[5][j] > phase0 )&&(phase0 > phase[5][j+1])
							   &&(phase[6][j] > phase0 )&&(phase0 > phase[6][j+1])
							   &&(phase[7][j] > phase0 )&&(phase0 > phase[7][j+1]))
					{goto happyend;}
											 
					m = m * (-1); 
		  }
	  }

					ABSORB;

		  }
	happyend:;

	  /* Output matters */

	/* transmit coordinates which were not changed, the rest overwrite below */
	Neutrons = InputNeutrons[i]; 



	/* translates neutron variables for output - X'= 0. . */

	{
		VectorType Path;

		Neutrons.Time = TOF + (- Pos[0]) / Dir[0] / V_FROM_LAMBDA(WL); 
			
		if(zerotime==1)
		{ 
			Neutrons.Time = fabs(fmod(Neutrons.Time + Phase/omega + coef_pi*M_PI/omega/2., coef_pi*M_PI/omega)) - coef_pi*M_PI/2./omega ;
		}

		CopyVector(Dir, Path);

		MultiplyByScalar(Path, (- Pos[0])/ Dir[0] );

		AddVector(Pos, Path);  /*Path = displacement vector */
		
		CopyVector(Pos, Neutrons.Position);
	}												 

	}
#endif

#if defined VITESS 

	  /* writes output binary file */
	  WriteNeutron(&Neutrons);

	}
}

  /* Do the general cleanup */

my_exit:;

  ChopperFermiCleanup();

  fprintf(LogFilePtr," \n");

  Cleanup(pos_ch[0],pos_ch[1],pos_ch[2], 0.0,0.0);	

  return 0;
}
#endif



#if defined VITESS || defined MCSTAS_SHARE

/* own initialization of the module */

void ChopperFermiInit(int argc, char *argv[])
{
  /*    INPUT  */

 #ifdef VITESS
  GeomFileName="chopper_fermi_g.dat";

  CurvGeomOption = 1;
 #endif

  while(argc>1)
    {
      char *arg = &argv[1][2];

      switch(argv[1][1])
	{
	case 'X':
	  sscanf(arg, "%lf", &pos_ch[0]);
	  break;
	
	case 'Y':
	  sscanf(arg, "%lf", &pos_ch[1]);
	  break;

	case 'V':
	  sscanf(arg, "%lf", &pos_ch[2]);
	  break;

	case 'a':
	  sscanf(arg, "%lf", &height);
	  break;

	case 'b':
	  sscanf(arg, "%lf", &width);
	  break;

	case 'c':
	  sscanf(arg, "%lf", &depth); 
	  break;

	case 'L':
	  sscanf(arg, "%lf", &optimal_wl);
	  break;

	case 'l':
	  sscanf(arg, "%ld", &Nchannels);
	  break;

	case 'm':
	  sscanf(arg, "%lf", &wallwidth);
	  break;

	case 'n':
	  sscanf(arg, "%lf", &omega);
	  if(omega == 0) omega = 1.E-6;
	  omega = omega * 2. * M_PI / 1000.;
	  break;

	case 'q':
	  sscanf(arg, "%lf", &Phase);
	  Phase = Phase * M_PI / 180.;
	  break;

	case 'O':
	  sscanf(arg, "%d", &Option);
	  if((Option!=1)&&(Option!=2)){fprintf(LogFilePtr,"\nERROR: Wrong option! Good options: 1-straight, 2-curved \n"); exit(-1); }
	  break;

	case 'g':
	  sscanf(arg, "%d", &CurvGeomOption);
	  break;

	case 'r':
	  sscanf(arg, "%lf", &diameter);
	  break;

	case 'G':
	  GeomFileName=arg;
	  break;

	case 'z':
	  sscanf(arg, "%d", &zerotime);
	  break;


	}
      argc--;
      argv++;
    }

	if(pos_ch[0] < diameter/2.) {fprintf(LogFilePtr,"\nERROR: Minimum position is diameter/2 \n"); exit(-1); }

	if(Nchannels==1) wallwidth=0.;

	
  /* calculate edge positions */ 
  
  if(Option==1)
    {

      fprintf(LogFilePtr,"\nStraight Fermi chopper option activated\n");

      /* diameter matter */

      main_depth = 2. * sqrt(sq(diameter/2.) - sq(width/2.));
      if(depth > main_depth) {
	fprintf(LogFilePtr,"\nERROR: Diameter too small - not compatible with 'channel length' and 'width'!\nTake min %f cm\n", 2. * sqrt(sq(depth/2.) + sq(width/2.)));

	exit(-1);
      }

	{
	  double add, w_ch;

	  long j, k, m;

	  if((     w_ch = ( width - (Nchannels + 1) * wallwidth ) / Nchannels    ) <= 0.)
	    {fprintf(LogFilePtr,"\nERROR: Channel width =< 0 !\n"); exit(-1);}

		fprintf(LogFilePtr,"Channel width: %f cm\n",w_ch);

	  for(k=0;k<4; k++) y_ch[k][0] =  - width/2.;

	  m = 1;

	  for(j=1; j< 2*Nchannels+2; j++)
	    {

	      if(m == 1) add = wallwidth; else add = w_ch;

	      for(k=0;k<4; k++) 
		  {
				x_ch[k][j] = - depth/2. * (3.- k)/3. + k/3. * depth/2.;

				y_ch[k][j] = y_ch[k][j-1] + add;
		  }

	      m = m * (-1);
	    }

	  for(k=0;k<4; k++) 
	  {		  
		x_ch[k][0] = x_ch[k][1] = x_ch[k][2*Nchannels] = x_ch[k][2*Nchannels+1] = - depth/2. *(3. - k)/3. + k/3. * depth/2.;
	  }

	  /* activating shadowing cylinder */

		x_ch[0][0] = x_ch[0][1] = x_ch[0][2*Nchannels] = x_ch[0][2*Nchannels+1] = - main_depth/2.;
	
		x_ch[3][0] = x_ch[3][1] = x_ch[3][2*Nchannels] = x_ch[3][2*Nchannels+1] = main_depth/2.;
	
	}

    }

  if(Option==2)
    {	
		  if(CurvGeomOption==1)
		  {
						  fprintf(LogFilePtr,"\nCurved Fermi chopper activated \n");

						  fprintf(LogFilePtr,"\nGeometry option: ideally shaped (close to parabolic) long channels ('channel length' inactive parameter) \n");

						  fprintf(LogFilePtr,"\nRadius of curvature (parabolic approximation):\n" 
											 "	%f cm at center\n" 
											 "	%f cm at circumference\n", 
											 V_FROM_LAMBDA(optimal_wl)/2./omega, V_FROM_LAMBDA(optimal_wl)/2./omega*pow((1 + sq(omega*diameter/V_FROM_LAMBDA(optimal_wl))), 1.5));


				  GeomFilePtr = fopen(GeomFileName,"w");
				  {
						double add, w_ch, xx[500], yy[500], tt;

						long j, m, k;

						if((     w_ch = ( width - (Nchannels + 1) * wallwidth ) / Nchannels    ) <= 0.)
						  {fprintf(LogFilePtr,"\nChannel width =< 0 !\n"); exit(-1);}

							fprintf(LogFilePtr,"Channel width: %f cm\n",w_ch);

							y_ch[0][0] =  - width/2.;
							x_ch[0][0] = sqrt(sq(diameter/2.) - sq(y_ch[0][0]));
							x_ch[3][0] = diameter/2. * sin(omega*diameter/V_FROM_LAMBDA(optimal_wl) + acos(2.*y_ch[0][0]/diameter));
							y_ch[3][0] = diameter/2. * cos(omega*diameter/V_FROM_LAMBDA(optimal_wl) + acos(2.*y_ch[0][0]/diameter));

						for(k=0; k<500; k++)
						  { tt = k /500. * 2 * M_PI; xx[k] = diameter/2. * cos(tt); yy[k] = diameter/2. * sin(tt); /* the circle */
						  fprintf(GeomFilePtr, " %ld %lf  %lf\n", 0, xx[k], yy[k]);
						  }

						m = 1;

						for(j=1; j< 2*Nchannels+2; j++)
						{

							if(m == 1) add = wallwidth; else add = w_ch;
							
							y_ch[0][j] = y_ch[0][j-1] + add;

							x_ch[0][j] = - sqrt(sq(diameter/2.) - sq(y_ch[0][j]));

							
							for(k=0; k<500; k++)
							  {
								tt= k /500. * 2. * fabs(x_ch[0][j])/V_FROM_LAMBDA(optimal_wl);
								xx[k] = y_ch[0][j] * sin(omega*tt) + (V_FROM_LAMBDA(optimal_wl)*tt - fabs(x_ch[0][j])) * cos(omega*tt);
								yy[k] = y_ch[0][j] * cos(omega*tt) - (V_FROM_LAMBDA(optimal_wl)*tt - fabs(x_ch[0][j])) * sin(omega*tt);

								fprintf(GeomFilePtr, "%ld   %lf  %lf\n", j, xx[k], yy[k]);		
							  }

								x_ch[1][j] = xx[167];
								y_ch[1][j] = yy[167];
								x_ch[2][j] = xx[333];
								y_ch[2][j] = yy[333];
								x_ch[3][j] = diameter/2. * sin(2.* omega*fabs(x_ch[0][j])/V_FROM_LAMBDA(optimal_wl) + acos(2.*y_ch[0][j]/diameter));
								y_ch[3][j] = diameter/2. * cos(2.* omega*fabs(x_ch[0][j])/V_FROM_LAMBDA(optimal_wl) + acos(2.*y_ch[0][j]/diameter));

								m = m * (-1);
						}
				  }
		  }
		if(CurvGeomOption==2)
		{

		  
					/* diameter matter */

				  main_depth = 2. * sqrt(sq(diameter/2.) - sq(width/2.));
				  if(depth > main_depth) {
					fprintf(LogFilePtr,"\nERROR: Diameter too small - not compatible with 'channel length' and 'width'!\nTake min %f cm\n", 2. * sqrt(sq(depth/2.) + sq(width/2.))); exit(-1);}


		  
				  fprintf(LogFilePtr,"\nCurved Fermi chopper activated \n");

				  fprintf(LogFilePtr,"\nGeometry option: circular shaped channels with fixed length (via 'channel length') \n");

				  radius_of_curv = V_FROM_LAMBDA(optimal_wl)/2./omega;

				  angle_channel = atan(depth/2./radius_of_curv);

				  fprintf(LogFilePtr,"\nRadius of curvature (parabolic approximation):\n" 
									 "optimal_velocity/2./omega = %f cm \n" 
									 "Angle channel: %f deg \n", 
									 radius_of_curv, 180./M_PI*angle_channel);


			  GeomFilePtr = fopen(GeomFileName,"w");
			  {
					double add, w_ch, xx[500], yy[500], tt;

					long j, m, k;

					if((     w_ch = ( width - (Nchannels + 1) * wallwidth ) / Nchannels    ) <= 0.)
					  {fprintf(LogFilePtr,"\nChannel width =< 0 !\n"); exit(-1);}

					fprintf(LogFilePtr,"Channel width: %f cm\n",w_ch);


					for(k=0; k<500; k++)
					  { tt = k /500. * 2 * M_PI; xx[k] = diameter/2. * cos(tt); yy[k] = diameter/2. * sin(tt); /* the circle */
					  fprintf(GeomFilePtr, " %ld %lf  %lf\n", -1, xx[k], yy[k]);
					  }

					m = 1;

					for(j=0; j< 2*Nchannels+2; j++)
					{

						if(m == 1) add = wallwidth; else add = w_ch;

								
						for(k=0; k<500; k++)
						  {
							xx[k]= (2.*k /499.  - 1.) * depth/2.;
							yy[k] = - width/2. - sqrt(sq(radius_of_curv) - sq(depth/2.)) + sqrt(sq(radius_of_curv) - sq(xx[k])) + shift_y;

							fprintf(GeomFilePtr, "%ld   %lf  %lf\n", j, xx[k], yy[k]);		
						  }

							x_ch[0][j] = xx[0];
							y_ch[0][j] = yy[0];
							x_ch[1][j] = xx[167];
							y_ch[1][j] = yy[167];
							x_ch[2][j] = xx[333];
							y_ch[2][j] = yy[333];
							x_ch[3][j] = xx[499];
							y_ch[3][j] = yy[499];

							shift_y += add;
							m = m * (-1);
					}
				  /* activating shadowing cylinder */

					x_ch[0][0] = x_ch[0][1] = x_ch[0][2*Nchannels] = x_ch[0][2*Nchannels+1] = - main_depth/2.;
				
					x_ch[3][0] = x_ch[3][1] = x_ch[3][2*Nchannels] = x_ch[3][2*Nchannels+1] = main_depth/2.;
	

			  }
		}
    }
		
  if(Option==1) coef_pi=1.; else coef_pi=2.;  fprintf(LogFilePtr,"Phase set is %f°.\n", 180./M_PI*fmod(Phase , coef_pi*M_PI));   

}/* End OwnInit */


/* own cleanup of the monochromator/analyser module */

void ChopperFermiCleanup()
{
  if (GeomFilePtr) fclose(GeomFilePtr);
}/* End OwnCleanup */

#endif
 /* include functions */
#undef MCSTAS_SHARE
#line 14657 "./Test_Fermi.c"
#undef wallwidth
#undef Phase
#undef diameter
#undef r_curv
#undef depth
#undef width
#undef height
#undef freq
#undef Ngates
#undef Nchannels
#undef zerotime
#undef GeomOption
#undef sGeomFileName
#undef optimal_wl
#undef omega
#undef pos_ch
#undef Neutrons
#undef Dir
#undef Pos
#undef GeomFilePtr
#undef XFILEName
#undef coef_pi
#undef x_ch
#undef y_ch
#undef phase0
#undef angle_channel
#undef shift_y
#undef main_depth
#undef radius_of_curv
#undef WL
#undef TOF
#undef CurvGeomOption
#undef Option
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'Monitor2_xt' [11]. */
#define mccompcurname  Monitor2_xt
#define mccompcurtype  Monitor_nD
#define mccompcurindex 11
#define user1 mccMonitor2_xt_user1
#define user2 mccMonitor2_xt_user2
#define user3 mccMonitor2_xt_user3
#define DEFS mccMonitor2_xt_DEFS
#define Vars mccMonitor2_xt_Vars
#define detector mccMonitor2_xt_detector
#define offdata mccMonitor2_xt_offdata
#define xwidth mccMonitor2_xt_xwidth
#define yheight mccMonitor2_xt_yheight
#define zdepth mccMonitor2_xt_zdepth
#define xmin mccMonitor2_xt_xmin
#define xmax mccMonitor2_xt_xmax
#define ymin mccMonitor2_xt_ymin
#define ymax mccMonitor2_xt_ymax
#define zmin mccMonitor2_xt_zmin
#define zmax mccMonitor2_xt_zmax
#define bins mccMonitor2_xt_bins
#define min mccMonitor2_xt_min
#define max mccMonitor2_xt_max
#define restore_neutron mccMonitor2_xt_restore_neutron
#define radius mccMonitor2_xt_radius
#define options mccMonitor2_xt_options
#define filename mccMonitor2_xt_filename
#define geometry mccMonitor2_xt_geometry
#define username1 mccMonitor2_xt_username1
#define username2 mccMonitor2_xt_username2
#define username3 mccMonitor2_xt_username3
#define nowritefile mccMonitor2_xt_nowritefile
#line 224 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/Monitor_nD.comp"
  MonitornD_Defines_type DEFS;
  MonitornD_Variables_type Vars;
  MCDETECTOR detector;
  off_struct offdata;
#line 14732 "./Test_Fermi.c"
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

Coords mcposaOrigin, mcposrOrigin;
Rotation mcrotaOrigin, mcrotrOrigin;
Coords mcposaSource, mcposrSource;
Rotation mcrotaSource, mcrotrSource;
Coords mcposaMonitor1_xt, mcposrMonitor1_xt;
Rotation mcrotaMonitor1_xt, mcrotrMonitor1_xt;
Coords mcposaFC_Position, mcposrFC_Position;
Rotation mcrotaFC_Position, mcrotrFC_Position;
Coords mcposaFC_GuideG, mcposrFC_GuideG;
Rotation mcrotaFC_GuideG, mcrotrFC_GuideG;
Coords mcposaFC_GuideC, mcposrFC_GuideC;
Rotation mcrotaFC_GuideC, mcrotrFC_GuideC;
Coords mcposaFC_McStas, mcposrFC_McStas;
Rotation mcrotaFC_McStas, mcrotrFC_McStas;
Coords mcposaFC_ILL, mcposrFC_ILL;
Rotation mcrotaFC_ILL, mcrotrFC_ILL;
Coords mcposaFake_Origin, mcposrFake_Origin;
Rotation mcrotaFake_Origin, mcrotrFake_Origin;
Coords mcposaFC_Vitess, mcposrFC_Vitess;
Rotation mcrotaFC_Vitess, mcrotrFC_Vitess;
Coords mcposaMonitor2_xt, mcposrMonitor2_xt;
Rotation mcrotaMonitor2_xt, mcrotrMonitor2_xt;

MCNUM mcnx, mcny, mcnz, mcnvx, mcnvy, mcnvz, mcnt, mcnsx, mcnsy, mcnsz, mcnp;

/* end declare */

void mcinit(void) {
#define mccompcurname  Test_Fermi
#define mccompcurtype  INSTRUMENT
#define mccompcurindex 0
#define mcposaTest_Fermi coords_set(0,0,0)
#define Fermi mcipFermi
#define lambda mciplambda
#define width_FC mcipwidth_FC
#define height_FC mcipheight_FC
#define length_FC mciplength_FC
#define FC_Hz mcipFC_Hz
#define Nslit_FC mcipNslit_FC
#define d_SF mcipd_SF
#define d_FD mcipd_FD
#define phase mcipphase
#line 49 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
{
  printf("\n%s: ", NAME_CURRENT_COMP);
  switch (Fermi) {
  case 1:
    printf("Using FermiChopper\n"); break;
  case 2:
    printf("Using Vitess_ChopperFermi\n"); break;
  case 3:
    printf("Using FermiChopper_ILL\n"); break;
  case 4:
    printf("Using rotating Guide_gravity\n"); break;
  case 5:
    printf("Using rotating Guide_channeled\n"); break;
  }

  double w = width_FC/Nslit_FC;
  double v = 3956/lambda;

  printf("\nTheor: Lambda=%g [Angs]\n",
    3956*w/2/PI/FC_Hz/length_FC/length_FC/2);

  printf("Theor: Time from source  t=%g [s]\n", d_SF/v);
  printf("Theor: Time to detection t=%g [s]\n", d_FD/v);
  printf("       Time period      dt=%g [s]\n", 1/FC_Hz);
  printf("       Slit pack div       %g [deg] (full width)\n", 2*atan2(w,length_FC)/PI*180);
  printf("       Time window width  =%g [s] (pulse width)\n",
    atan2(w,length_FC)/PI/FC_Hz);
  printf("       Phase           phi=%g [deg]\n", (d_SF/v)/(1/FC_Hz)*360);

  time_to_arrival  = (d_SF/v);
  time_window_width= atan2(w,length_FC)/PI/FC_Hz;
  if (phase == -0) phase=(d_SF/v)/(1/FC_Hz)*360; /* assumes time at source is centered on 0 */
}
#line 14841 "./Test_Fermi.c"
#undef phase
#undef d_FD
#undef d_SF
#undef Nslit_FC
#undef FC_Hz
#undef length_FC
#undef height_FC
#undef width_FC
#undef lambda
#undef Fermi
#undef mcposaTest_Fermi
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
#line 39 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  if("NULL") strncpy(mccOrigin_profile, "NULL" ? "NULL" : "", 16384); else mccOrigin_profile[0]='\0';
#line 39 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccOrigin_percent = 10;
#line 39 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccOrigin_flag_save = 0;
#line 39 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccOrigin_minutes = 0;
#line 14877 "./Test_Fermi.c"

  SIG_MESSAGE("Origin (Init:Place/Rotate)");
  rot_set_rotation(mcrotaOrigin,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 14884 "./Test_Fermi.c"
  rot_copy(mcrotrOrigin, mcrotaOrigin);
  mcposaOrigin = coords_set(
#line 86 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
    0,
#line 86 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
    0,
#line 86 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
    0);
#line 14893 "./Test_Fermi.c"
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
#line 129 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  if("NULL") strncpy(mccSource_flux_file, "NULL" ? "NULL" : "", 16384); else mccSource_flux_file[0]='\0';
#line 129 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  if("NULL") strncpy(mccSource_xdiv_file, "NULL" ? "NULL" : "", 16384); else mccSource_xdiv_file[0]='\0';
#line 129 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  if("NULL") strncpy(mccSource_ydiv_file, "NULL" ? "NULL" : "", 16384); else mccSource_ydiv_file[0]='\0';
#line 130 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccSource_radius = 0.0;
#line 130 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccSource_dist = 0;
#line 90 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccSource_focus_xw = mcipwidth_FC;
#line 90 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccSource_focus_yh = mcipheight_FC;
#line 130 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccSource_focus_aw = 0;
#line 130 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccSource_focus_ah = 0;
#line 131 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccSource_E0 = 0;
#line 131 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccSource_dE = 0;
#line 90 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccSource_lambda0 = mciplambda;
#line 90 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccSource_dlambda = 0.3;
#line 131 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccSource_I1 = 1;
#line 91 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccSource_yheight = mcipheight_FC;
#line 91 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccSource_xwidth = mcipwidth_FC;
#line 132 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccSource_verbose = 0;
#line 132 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccSource_T1 = 0;
#line 133 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccSource_flux_file_perAA = 0;
#line 133 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccSource_flux_file_log = 0;
#line 134 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccSource_Lmin = 0;
#line 134 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccSource_Lmax = 0;
#line 134 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccSource_Emin = 0;
#line 134 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccSource_Emax = 0;
#line 134 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccSource_T2 = 0;
#line 134 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccSource_I2 = 0;
#line 134 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccSource_T3 = 0;
#line 134 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccSource_I3 = 0;
#line 134 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccSource_zdepth = 0;
#line 134 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccSource_target_index = + 1;
#line 14964 "./Test_Fermi.c"

  SIG_MESSAGE("Source (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 14971 "./Test_Fermi.c"
  rot_mul(mctr1, mcrotaOrigin, mcrotaSource);
  rot_transpose(mcrotaOrigin, mctr1);
  rot_mul(mcrotaSource, mctr1, mcrotrSource);
  mctc1 = coords_set(
#line 92 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
    0,
#line 92 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
    0,
#line 92 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
    0);
#line 14982 "./Test_Fermi.c"
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
    /* Component Monitor1_xt. */
  /* Setting parameters for component Monitor1_xt. */
  SIG_MESSAGE("Monitor1_xt (Init:SetPar)");
#line 102 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccMonitor1_xt_xwidth = mcipwidth_FC;
#line 102 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccMonitor1_xt_yheight = mcipheight_FC;
#line 203 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccMonitor1_xt_zdepth = 0;
#line 204 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccMonitor1_xt_xmin = 0;
#line 204 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccMonitor1_xt_xmax = 0;
#line 204 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccMonitor1_xt_ymin = 0;
#line 204 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccMonitor1_xt_ymax = 0;
#line 204 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccMonitor1_xt_zmin = 0;
#line 204 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccMonitor1_xt_zmax = 0;
#line 205 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccMonitor1_xt_bins = 0;
#line 205 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccMonitor1_xt_min = -1e40;
#line 205 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccMonitor1_xt_max = 1e40;
#line 205 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccMonitor1_xt_restore_neutron = 0;
#line 205 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccMonitor1_xt_radius = 0;
#line 102 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  if("multiple x y time, all auto") strncpy(mccMonitor1_xt_options, "multiple x y time, all auto" ? "multiple x y time, all auto" : "", 16384); else mccMonitor1_xt_options[0]='\0';
#line 206 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  if("NULL") strncpy(mccMonitor1_xt_filename, "NULL" ? "NULL" : "", 16384); else mccMonitor1_xt_filename[0]='\0';
#line 206 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  if("NULL") strncpy(mccMonitor1_xt_geometry, "NULL" ? "NULL" : "", 16384); else mccMonitor1_xt_geometry[0]='\0';
#line 207 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  if("NULL") strncpy(mccMonitor1_xt_username1, "NULL" ? "NULL" : "", 16384); else mccMonitor1_xt_username1[0]='\0';
#line 207 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  if("NULL") strncpy(mccMonitor1_xt_username2, "NULL" ? "NULL" : "", 16384); else mccMonitor1_xt_username2[0]='\0';
#line 207 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  if("NULL") strncpy(mccMonitor1_xt_username3, "NULL" ? "NULL" : "", 16384); else mccMonitor1_xt_username3[0]='\0';
#line 208 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccMonitor1_xt_nowritefile = 0;
#line 15038 "./Test_Fermi.c"

  SIG_MESSAGE("Monitor1_xt (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 15045 "./Test_Fermi.c"
  rot_mul(mctr1, mcrotaSource, mcrotaMonitor1_xt);
  rot_transpose(mcrotaSource, mctr1);
  rot_mul(mcrotaMonitor1_xt, mctr1, mcrotrMonitor1_xt);
  mctc1 = coords_set(
#line 103 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
    0,
#line 103 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
    0,
#line 103 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
    mcipd_SF -0.1);
#line 15056 "./Test_Fermi.c"
  rot_transpose(mcrotaSource, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaMonitor1_xt = coords_add(mcposaSource, mctc2);
  mctc1 = coords_sub(mcposaSource, mcposaMonitor1_xt);
  mcposrMonitor1_xt = rot_apply(mcrotaMonitor1_xt, mctc1);
  mcDEBUG_COMPONENT("Monitor1_xt", mcposaMonitor1_xt, mcrotaMonitor1_xt)
  mccomp_posa[3] = mcposaMonitor1_xt;
  mccomp_posr[3] = mcposrMonitor1_xt;
  mcNCounter[3]  = mcPCounter[3] = mcP2Counter[3] = 0;
  mcAbsorbProp[3]= 0;
    /* Component FC_Position. */
  /* Setting parameters for component FC_Position. */
  SIG_MESSAGE("FC_Position (Init:SetPar)");

  SIG_MESSAGE("FC_Position (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 15076 "./Test_Fermi.c"
  rot_mul(mctr1, mcrotaMonitor1_xt, mcrotaFC_Position);
  rot_transpose(mcrotaMonitor1_xt, mctr1);
  rot_mul(mcrotaFC_Position, mctr1, mcrotrFC_Position);
  mctc1 = coords_set(
#line 106 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
    0,
#line 106 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
    0,
#line 106 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
    0.1);
#line 15087 "./Test_Fermi.c"
  rot_transpose(mcrotaMonitor1_xt, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaFC_Position = coords_add(mcposaMonitor1_xt, mctc2);
  mctc1 = coords_sub(mcposaMonitor1_xt, mcposaFC_Position);
  mcposrFC_Position = rot_apply(mcrotaFC_Position, mctc1);
  mcDEBUG_COMPONENT("FC_Position", mcposaFC_Position, mcrotaFC_Position)
  mccomp_posa[4] = mcposaFC_Position;
  mccomp_posr[4] = mcposrFC_Position;
  mcNCounter[4]  = mcPCounter[4] = mcP2Counter[4] = 0;
  mcAbsorbProp[4]= 0;
    /* Component FC_GuideG. */
  /* Setting parameters for component FC_GuideG. */
  SIG_MESSAGE("FC_GuideG (Init:SetPar)");
#line 109 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccFC_GuideG_w1 = mcipwidth_FC;
#line 109 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccFC_GuideG_h1 = mcipheight_FC;
#line 113 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccFC_GuideG_w2 = 0;
#line 113 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccFC_GuideG_h2 = 0;
#line 109 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccFC_GuideG_l = mciplength_FC;
#line 109 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccFC_GuideG_R0 = 0.0;
#line 114 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccFC_GuideG_Qc = 0.0218;
#line 114 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccFC_GuideG_alpha = 4.38;
#line 114 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccFC_GuideG_m = 1.0;
#line 114 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccFC_GuideG_W = 0.003;
#line 109 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccFC_GuideG_nslit = mcipNslit_FC;
#line 109 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccFC_GuideG_d = 0;
#line 115 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccFC_GuideG_mleft = -1;
#line 115 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccFC_GuideG_mright = -1;
#line 115 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccFC_GuideG_mtop = -1;
#line 115 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccFC_GuideG_mbottom = -1;
#line 115 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccFC_GuideG_nhslit = 1;
#line 115 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccFC_GuideG_G = 0;
#line 116 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccFC_GuideG_aleft = -1;
#line 116 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccFC_GuideG_aright = -1;
#line 116 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccFC_GuideG_atop = -1;
#line 116 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccFC_GuideG_abottom = -1;
#line 117 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccFC_GuideG_wavy = 0;
#line 117 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccFC_GuideG_wavy_z = 0;
#line 117 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccFC_GuideG_wavy_tb = 0;
#line 117 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccFC_GuideG_wavy_lr = 0;
#line 118 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccFC_GuideG_chamfers = 0;
#line 118 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccFC_GuideG_chamfers_z = 0;
#line 118 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccFC_GuideG_chamfers_lr = 0;
#line 118 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccFC_GuideG_chamfers_tb = 0;
#line 118 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccFC_GuideG_nelements = 1;
#line 109 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccFC_GuideG_nu = mcipFC_Hz;
#line 109 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccFC_GuideG_phase = mcipphase;
#line 119 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  if("NULL") strncpy(mccFC_GuideG_reflect, "NULL" ? "NULL" : "", 16384); else mccFC_GuideG_reflect[0]='\0';
#line 15169 "./Test_Fermi.c"

  SIG_MESSAGE("FC_GuideG (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 15176 "./Test_Fermi.c"
  rot_mul(mctr1, mcrotaFC_Position, mcrotaFC_GuideG);
  rot_transpose(mcrotaFC_Position, mctr1);
  rot_mul(mcrotaFC_GuideG, mctr1, mcrotrFC_GuideG);
  mctc1 = coords_set(
#line 111 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
    0,
#line 111 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
    0,
#line 111 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
    - mciplength_FC / 2.0);
#line 15187 "./Test_Fermi.c"
  rot_transpose(mcrotaFC_Position, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaFC_GuideG = coords_add(mcposaFC_Position, mctc2);
  mctc1 = coords_sub(mcposaFC_Position, mcposaFC_GuideG);
  mcposrFC_GuideG = rot_apply(mcrotaFC_GuideG, mctc1);
  mcDEBUG_COMPONENT("FC_GuideG", mcposaFC_GuideG, mcrotaFC_GuideG)
  mccomp_posa[5] = mcposaFC_GuideG;
  mccomp_posr[5] = mcposrFC_GuideG;
  mcNCounter[5]  = mcPCounter[5] = mcP2Counter[5] = 0;
  mcAbsorbProp[5]= 0;
    /* Component FC_GuideC. */
  /* Setting parameters for component FC_GuideC. */
  SIG_MESSAGE("FC_GuideC (Init:SetPar)");
#line 117 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccFC_GuideC_w1 = mcipwidth_FC;
#line 117 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccFC_GuideC_h1 = mcipheight_FC;
#line 70 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccFC_GuideC_w2 = 0;
#line 70 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccFC_GuideC_h2 = 0;
#line 117 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccFC_GuideC_l = mciplength_FC;
#line 117 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccFC_GuideC_R0 = 0.0;
#line 71 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccFC_GuideC_Qc = 0;
#line 71 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccFC_GuideC_alpha = 0;
#line 71 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccFC_GuideC_m = 0;
#line 117 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccFC_GuideC_nslit = mcipNslit_FC;
#line 117 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccFC_GuideC_d = 0;
#line 72 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccFC_GuideC_Qcx = 0.0218;
#line 72 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccFC_GuideC_Qcy = 0.0218;
#line 72 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccFC_GuideC_alphax = 4.38;
#line 72 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccFC_GuideC_alphay = 4.38;
#line 72 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccFC_GuideC_W = 0.003;
#line 72 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccFC_GuideC_mx = 1;
#line 72 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccFC_GuideC_my = 1;
#line 117 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccFC_GuideC_nu = mcipFC_Hz;
#line 117 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccFC_GuideC_phase = mcipphase;
#line 15241 "./Test_Fermi.c"

  SIG_MESSAGE("FC_GuideC (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 15248 "./Test_Fermi.c"
  rot_mul(mctr1, mcrotaFC_Position, mcrotaFC_GuideC);
  rot_transpose(mcrotaFC_GuideG, mctr1);
  rot_mul(mcrotaFC_GuideC, mctr1, mcrotrFC_GuideC);
  mctc1 = coords_set(
#line 119 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
    0,
#line 119 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
    0,
#line 119 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
    - mciplength_FC / 2.0);
#line 15259 "./Test_Fermi.c"
  rot_transpose(mcrotaFC_Position, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaFC_GuideC = coords_add(mcposaFC_Position, mctc2);
  mctc1 = coords_sub(mcposaFC_GuideG, mcposaFC_GuideC);
  mcposrFC_GuideC = rot_apply(mcrotaFC_GuideC, mctc1);
  mcDEBUG_COMPONENT("FC_GuideC", mcposaFC_GuideC, mcrotaFC_GuideC)
  mccomp_posa[6] = mcposaFC_GuideC;
  mccomp_posr[6] = mcposrFC_GuideC;
  mcNCounter[6]  = mcPCounter[6] = mcP2Counter[6] = 0;
  mcAbsorbProp[6]= 0;
    /* Component FC_McStas. */
  /* Setting parameters for component FC_McStas. */
  SIG_MESSAGE("FC_McStas (Init:SetPar)");
#line 126 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccFC_McStas_phase = mcipphase;
#line 125 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccFC_McStas_radius = 0.1;
#line 125 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccFC_McStas_nu = mcipFC_Hz;
#line 89 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccFC_McStas_w = 0.00022475;
#line 125 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccFC_McStas_nslit = mcipNslit_FC;
#line 89 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccFC_McStas_R0 = 0.0;
#line 90 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccFC_McStas_Qc = 0.02176;
#line 90 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccFC_McStas_alpha = 2.33;
#line 90 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccFC_McStas_m = 0.0;
#line 90 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccFC_McStas_W = 2e-3;
#line 125 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccFC_McStas_length = mciplength_FC;
#line 90 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccFC_McStas_eff = 0.95;
#line 91 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccFC_McStas_zero_time = 0;
#line 125 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccFC_McStas_xwidth = mcipwidth_FC;
#line 126 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccFC_McStas_verbose = 1;
#line 125 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccFC_McStas_yheight = mcipheight_FC;
#line 92 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccFC_McStas_curvature = 0;
#line 92 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccFC_McStas_delay = 0;
#line 15309 "./Test_Fermi.c"

  SIG_MESSAGE("FC_McStas (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 15316 "./Test_Fermi.c"
  rot_mul(mctr1, mcrotaFC_Position, mcrotaFC_McStas);
  rot_transpose(mcrotaFC_GuideC, mctr1);
  rot_mul(mcrotaFC_McStas, mctr1, mcrotrFC_McStas);
  mctc1 = coords_set(
#line 128 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
    0,
#line 128 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
    0,
#line 128 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
    0);
#line 15327 "./Test_Fermi.c"
  rot_transpose(mcrotaFC_Position, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaFC_McStas = coords_add(mcposaFC_Position, mctc2);
  mctc1 = coords_sub(mcposaFC_GuideC, mcposaFC_McStas);
  mcposrFC_McStas = rot_apply(mcrotaFC_McStas, mctc1);
  mcDEBUG_COMPONENT("FC_McStas", mcposaFC_McStas, mcrotaFC_McStas)
  mccomp_posa[7] = mcposaFC_McStas;
  mccomp_posr[7] = mcposrFC_McStas;
  mcNCounter[7]  = mcPCounter[7] = mcP2Counter[7] = 0;
  mcAbsorbProp[7]= 0;
    /* Component FC_ILL. */
  /* Setting parameters for component FC_ILL. */
  SIG_MESSAGE("FC_ILL (Init:SetPar)");
#line 132 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccFC_ILL_phase = mcipphase;
#line 131 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccFC_ILL_radius = 0.1;
#line 131 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccFC_ILL_nu = mcipFC_Hz;
#line 131 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccFC_ILL_yheight = mcipheight_FC;
#line 96 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccFC_ILL_w = 0.00022475;
#line 132 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccFC_ILL_nslit = mcipNslit_FC;
#line 96 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccFC_ILL_R0 = 0.0;
#line 97 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccFC_ILL_Qc = 0.02176;
#line 97 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccFC_ILL_alpha = 2.33;
#line 97 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccFC_ILL_m = 0.0;
#line 97 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccFC_ILL_W = 2e-3;
#line 132 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccFC_ILL_length = mciplength_FC;
#line 97 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccFC_ILL_eff = 0.95;
#line 98 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccFC_ILL_zero_time = 0;
#line 132 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccFC_ILL_xwidth = mcipwidth_FC;
#line 98 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccFC_ILL_verbose = 0;
#line 15373 "./Test_Fermi.c"

  SIG_MESSAGE("FC_ILL (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 15380 "./Test_Fermi.c"
  rot_mul(mctr1, mcrotaFC_Position, mcrotaFC_ILL);
  rot_transpose(mcrotaFC_McStas, mctr1);
  rot_mul(mcrotaFC_ILL, mctr1, mcrotrFC_ILL);
  mctc1 = coords_set(
#line 134 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
    0,
#line 134 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
    0,
#line 134 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
    0);
#line 15391 "./Test_Fermi.c"
  rot_transpose(mcrotaFC_Position, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaFC_ILL = coords_add(mcposaFC_Position, mctc2);
  mctc1 = coords_sub(mcposaFC_McStas, mcposaFC_ILL);
  mcposrFC_ILL = rot_apply(mcrotaFC_ILL, mctc1);
  mcDEBUG_COMPONENT("FC_ILL", mcposaFC_ILL, mcrotaFC_ILL)
  mccomp_posa[8] = mcposaFC_ILL;
  mccomp_posr[8] = mcposrFC_ILL;
  mcNCounter[8]  = mcPCounter[8] = mcP2Counter[8] = 0;
  mcAbsorbProp[8]= 0;
    /* Component Fake_Origin. */
  /* Setting parameters for component Fake_Origin. */
  SIG_MESSAGE("Fake_Origin (Init:SetPar)");

  SIG_MESSAGE("Fake_Origin (Init:Place/Rotate)");
  rot_set_rotation(mcrotaFake_Origin,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 15411 "./Test_Fermi.c"
  rot_transpose(mcrotaFC_ILL, mctr1);
  rot_mul(mcrotaFake_Origin, mctr1, mcrotrFake_Origin);
  mcposaFake_Origin = coords_set(
#line 137 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
    0,
#line 137 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
    0,
#line 137 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
    0);
#line 15421 "./Test_Fermi.c"
  mctc1 = coords_sub(mcposaFC_ILL, mcposaFake_Origin);
  mcposrFake_Origin = rot_apply(mcrotaFake_Origin, mctc1);
  mcDEBUG_COMPONENT("Fake_Origin", mcposaFake_Origin, mcrotaFake_Origin)
  mccomp_posa[9] = mcposaFake_Origin;
  mccomp_posr[9] = mcposrFake_Origin;
  mcNCounter[9]  = mcPCounter[9] = mcP2Counter[9] = 0;
  mcAbsorbProp[9]= 0;
    /* Component FC_Vitess. */
  /* Setting parameters for component FC_Vitess. */
  SIG_MESSAGE("FC_Vitess (Init:SetPar)");
#line 140 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  if("FC_geom_str.dat") strncpy(mccFC_Vitess_sGeomFileName, "FC_geom_str.dat" ? "FC_geom_str.dat" : "", 16384); else mccFC_Vitess_sGeomFileName[0]='\0';
#line 119 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccFC_Vitess_GeomOption = 0;
#line 119 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccFC_Vitess_zerotime = 0;
#line 140 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccFC_Vitess_Nchannels = mcipNslit_FC;
#line 119 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccFC_Vitess_Ngates = 4;
#line 141 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccFC_Vitess_freq = mcipFC_Hz;
#line 141 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccFC_Vitess_height = mcipheight_FC;
#line 141 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccFC_Vitess_width = mcipwidth_FC;
#line 141 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccFC_Vitess_depth = mciplength_FC;
#line 121 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccFC_Vitess_r_curv = 0.5;
#line 142 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccFC_Vitess_diameter = 0.1;
#line 142 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccFC_Vitess_Phase = mcipphase;
#line 140 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccFC_Vitess_wallwidth = 0;
#line 15458 "./Test_Fermi.c"

  SIG_MESSAGE("FC_Vitess (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 15465 "./Test_Fermi.c"
  rot_mul(mctr1, mcrotaFake_Origin, mcrotaFC_Vitess);
  rot_transpose(mcrotaFake_Origin, mctr1);
  rot_mul(mcrotaFC_Vitess, mctr1, mcrotrFC_Vitess);
  mctc1 = coords_set(
#line 144 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
    0,
#line 144 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
    0,
#line 144 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
    mcipd_SF);
#line 15476 "./Test_Fermi.c"
  rot_transpose(mcrotaFake_Origin, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaFC_Vitess = coords_add(mcposaFake_Origin, mctc2);
  mctc1 = coords_sub(mcposaFake_Origin, mcposaFC_Vitess);
  mcposrFC_Vitess = rot_apply(mcrotaFC_Vitess, mctc1);
  mcDEBUG_COMPONENT("FC_Vitess", mcposaFC_Vitess, mcrotaFC_Vitess)
  mccomp_posa[10] = mcposaFC_Vitess;
  mccomp_posr[10] = mcposrFC_Vitess;
  mcNCounter[10]  = mcPCounter[10] = mcP2Counter[10] = 0;
  mcAbsorbProp[10]= 0;
    /* Component Monitor2_xt. */
  /* Setting parameters for component Monitor2_xt. */
  SIG_MESSAGE("Monitor2_xt (Init:SetPar)");
#line 149 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccMonitor2_xt_xwidth = mcipwidth_FC;
#line 149 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccMonitor2_xt_yheight = mcipheight_FC;
#line 203 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccMonitor2_xt_zdepth = 0;
#line 204 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccMonitor2_xt_xmin = 0;
#line 204 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccMonitor2_xt_xmax = 0;
#line 204 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccMonitor2_xt_ymin = 0;
#line 204 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccMonitor2_xt_ymax = 0;
#line 204 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccMonitor2_xt_zmin = 0;
#line 204 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccMonitor2_xt_zmax = 0;
#line 205 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccMonitor2_xt_bins = 0;
#line 205 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccMonitor2_xt_min = -1e40;
#line 205 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccMonitor2_xt_max = 1e40;
#line 148 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccMonitor2_xt_restore_neutron = 1;
#line 205 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccMonitor2_xt_radius = 0;
#line 147 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  if("multiple x y time, all auto") strncpy(mccMonitor2_xt_options, "multiple x y time, all auto" ? "multiple x y time, all auto" : "", 16384); else mccMonitor2_xt_options[0]='\0';
#line 206 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  if("NULL") strncpy(mccMonitor2_xt_filename, "NULL" ? "NULL" : "", 16384); else mccMonitor2_xt_filename[0]='\0';
#line 206 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  if("NULL") strncpy(mccMonitor2_xt_geometry, "NULL" ? "NULL" : "", 16384); else mccMonitor2_xt_geometry[0]='\0';
#line 207 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  if("NULL") strncpy(mccMonitor2_xt_username1, "NULL" ? "NULL" : "", 16384); else mccMonitor2_xt_username1[0]='\0';
#line 207 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  if("NULL") strncpy(mccMonitor2_xt_username2, "NULL" ? "NULL" : "", 16384); else mccMonitor2_xt_username2[0]='\0';
#line 207 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  if("NULL") strncpy(mccMonitor2_xt_username3, "NULL" ? "NULL" : "", 16384); else mccMonitor2_xt_username3[0]='\0';
#line 208 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  mccMonitor2_xt_nowritefile = 0;
#line 15532 "./Test_Fermi.c"

  SIG_MESSAGE("Monitor2_xt (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 15539 "./Test_Fermi.c"
  rot_mul(mctr1, mcrotaFC_Position, mcrotaMonitor2_xt);
  rot_transpose(mcrotaFC_Vitess, mctr1);
  rot_mul(mcrotaMonitor2_xt, mctr1, mcrotrMonitor2_xt);
  mctc1 = coords_set(
#line 150 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
    0,
#line 150 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
    0,
#line 150 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
    mcipd_FD);
#line 15550 "./Test_Fermi.c"
  rot_transpose(mcrotaFC_Position, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaMonitor2_xt = coords_add(mcposaFC_Position, mctc2);
  mctc1 = coords_sub(mcposaFC_Vitess, mcposaMonitor2_xt);
  mcposrMonitor2_xt = rot_apply(mcrotaMonitor2_xt, mctc1);
  mcDEBUG_COMPONENT("Monitor2_xt", mcposaMonitor2_xt, mcrotaMonitor2_xt)
  mccomp_posa[11] = mcposaMonitor2_xt;
  mccomp_posr[11] = mcposrMonitor2_xt;
  mcNCounter[11]  = mcPCounter[11] = mcP2Counter[11] = 0;
  mcAbsorbProp[11]= 0;
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
#line 15587 "./Test_Fermi.c"
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
#define mccompcurtype  Source_gen
#define mccompcurindex 2
#define p_in mccSource_p_in
#define lambda1 mccSource_lambda1
#define lambda2 mccSource_lambda2
#define lambda3 mccSource_lambda3
#define pTable mccSource_pTable
#define pTable_x mccSource_pTable_x
#define pTable_y mccSource_pTable_y
#define pTable_xmin mccSource_pTable_xmin
#define pTable_xmax mccSource_pTable_xmax
#define pTable_xsum mccSource_pTable_xsum
#define pTable_ymin mccSource_pTable_ymin
#define pTable_ymax mccSource_pTable_ymax
#define pTable_ysum mccSource_pTable_ysum
#define pTable_dxmin mccSource_pTable_dxmin
#define pTable_dxmax mccSource_pTable_dxmax
#define pTable_dymin mccSource_pTable_dymin
#define pTable_dymax mccSource_pTable_dymax
#define flux_file mccSource_flux_file
#define xdiv_file mccSource_xdiv_file
#define ydiv_file mccSource_ydiv_file
#define radius mccSource_radius
#define dist mccSource_dist
#define focus_xw mccSource_focus_xw
#define focus_yh mccSource_focus_yh
#define focus_aw mccSource_focus_aw
#define focus_ah mccSource_focus_ah
#define E0 mccSource_E0
#define dE mccSource_dE
#define lambda0 mccSource_lambda0
#define dlambda mccSource_dlambda
#define I1 mccSource_I1
#define yheight mccSource_yheight
#define xwidth mccSource_xwidth
#define verbose mccSource_verbose
#define T1 mccSource_T1
#define flux_file_perAA mccSource_flux_file_perAA
#define flux_file_log mccSource_flux_file_log
#define Lmin mccSource_Lmin
#define Lmax mccSource_Lmax
#define Emin mccSource_Emin
#define Emax mccSource_Emax
#define T2 mccSource_T2
#define I2 mccSource_I2
#define T3 mccSource_T3
#define I3 mccSource_I3
#define zdepth mccSource_zdepth
#define target_index mccSource_target_index
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
#line 15924 "./Test_Fermi.c"
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

  /* Initializations for component Monitor1_xt. */
  SIG_MESSAGE("Monitor1_xt (Init)");
#define mccompcurname  Monitor1_xt
#define mccompcurtype  Monitor_nD
#define mccompcurindex 3
#define user1 mccMonitor1_xt_user1
#define user2 mccMonitor1_xt_user2
#define user3 mccMonitor1_xt_user3
#define DEFS mccMonitor1_xt_DEFS
#define Vars mccMonitor1_xt_Vars
#define detector mccMonitor1_xt_detector
#define offdata mccMonitor1_xt_offdata
#define xwidth mccMonitor1_xt_xwidth
#define yheight mccMonitor1_xt_yheight
#define zdepth mccMonitor1_xt_zdepth
#define xmin mccMonitor1_xt_xmin
#define xmax mccMonitor1_xt_xmax
#define ymin mccMonitor1_xt_ymin
#define ymax mccMonitor1_xt_ymax
#define zmin mccMonitor1_xt_zmin
#define zmax mccMonitor1_xt_zmax
#define bins mccMonitor1_xt_bins
#define min mccMonitor1_xt_min
#define max mccMonitor1_xt_max
#define restore_neutron mccMonitor1_xt_restore_neutron
#define radius mccMonitor1_xt_radius
#define options mccMonitor1_xt_options
#define filename mccMonitor1_xt_filename
#define geometry mccMonitor1_xt_geometry
#define username1 mccMonitor1_xt_username1
#define username2 mccMonitor1_xt_username2
#define username3 mccMonitor1_xt_username3
#define nowritefile mccMonitor1_xt_nowritefile
#line 231 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/Monitor_nD.comp"
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
#line 16089 "./Test_Fermi.c"
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

  /* Initializations for component FC_Position. */
  SIG_MESSAGE("FC_Position (Init)");

  /* Initializations for component FC_GuideG. */
  SIG_MESSAGE("FC_GuideG (Init)");
#define mccompcurname  FC_GuideG
#define mccompcurtype  Guide_gravity
#define mccompcurindex 5
#define GVars mccFC_GuideG_GVars
#define pTable mccFC_GuideG_pTable
#define w1 mccFC_GuideG_w1
#define h1 mccFC_GuideG_h1
#define w2 mccFC_GuideG_w2
#define h2 mccFC_GuideG_h2
#define l mccFC_GuideG_l
#define R0 mccFC_GuideG_R0
#define Qc mccFC_GuideG_Qc
#define alpha mccFC_GuideG_alpha
#define m mccFC_GuideG_m
#define W mccFC_GuideG_W
#define nslit mccFC_GuideG_nslit
#define d mccFC_GuideG_d
#define mleft mccFC_GuideG_mleft
#define mright mccFC_GuideG_mright
#define mtop mccFC_GuideG_mtop
#define mbottom mccFC_GuideG_mbottom
#define nhslit mccFC_GuideG_nhslit
#define G mccFC_GuideG_G
#define aleft mccFC_GuideG_aleft
#define aright mccFC_GuideG_aright
#define atop mccFC_GuideG_atop
#define abottom mccFC_GuideG_abottom
#define wavy mccFC_GuideG_wavy
#define wavy_z mccFC_GuideG_wavy_z
#define wavy_tb mccFC_GuideG_wavy_tb
#define wavy_lr mccFC_GuideG_wavy_lr
#define chamfers mccFC_GuideG_chamfers
#define chamfers_z mccFC_GuideG_chamfers_z
#define chamfers_lr mccFC_GuideG_chamfers_lr
#define chamfers_tb mccFC_GuideG_chamfers_tb
#define nelements mccFC_GuideG_nelements
#define nu mccFC_GuideG_nu
#define phase mccFC_GuideG_phase
#define reflect mccFC_GuideG_reflect
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
#line 16218 "./Test_Fermi.c"
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

  /* Initializations for component FC_GuideC. */
  SIG_MESSAGE("FC_GuideC (Init)");
#define mccompcurname  FC_GuideC
#define mccompcurtype  Guide_channeled
#define mccompcurindex 6
#define w1c mccFC_GuideC_w1c
#define w2c mccFC_GuideC_w2c
#define ww mccFC_GuideC_ww
#define hh mccFC_GuideC_hh
#define whalf mccFC_GuideC_whalf
#define hhalf mccFC_GuideC_hhalf
#define lwhalf mccFC_GuideC_lwhalf
#define lhhalf mccFC_GuideC_lhhalf
#define w1 mccFC_GuideC_w1
#define h1 mccFC_GuideC_h1
#define w2 mccFC_GuideC_w2
#define h2 mccFC_GuideC_h2
#define l mccFC_GuideC_l
#define R0 mccFC_GuideC_R0
#define Qc mccFC_GuideC_Qc
#define alpha mccFC_GuideC_alpha
#define m mccFC_GuideC_m
#define nslit mccFC_GuideC_nslit
#define d mccFC_GuideC_d
#define Qcx mccFC_GuideC_Qcx
#define Qcy mccFC_GuideC_Qcy
#define alphax mccFC_GuideC_alphax
#define alphay mccFC_GuideC_alphay
#define W mccFC_GuideC_W
#define mx mccFC_GuideC_mx
#define my mccFC_GuideC_my
#define nu mccFC_GuideC_nu
#define phase mccFC_GuideC_phase
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
#line 16333 "./Test_Fermi.c"
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

  /* Initializations for component FC_McStas. */
  SIG_MESSAGE("FC_McStas (Init)");
#define mccompcurname  FC_McStas
#define mccompcurtype  FermiChopper
#define mccompcurindex 7
#define FCVars mccFC_McStas_FCVars
#define phase mccFC_McStas_phase
#define radius mccFC_McStas_radius
#define nu mccFC_McStas_nu
#define w mccFC_McStas_w
#define nslit mccFC_McStas_nslit
#define R0 mccFC_McStas_R0
#define Qc mccFC_McStas_Qc
#define alpha mccFC_McStas_alpha
#define m mccFC_McStas_m
#define W mccFC_McStas_W
#define length mccFC_McStas_length
#define eff mccFC_McStas_eff
#define zero_time mccFC_McStas_zero_time
#define xwidth mccFC_McStas_xwidth
#define verbose mccFC_McStas_verbose
#define yheight mccFC_McStas_yheight
#define curvature mccFC_McStas_curvature
#define delay mccFC_McStas_delay
#line 282 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/FermiChopper.comp"
{

/************************ CALCULATION CONSTANTS *****************************/
  strcpy(FCVars.compcurname, NAME_CURRENT_COMP);
  
  FCVars.omega    = 2*PI*nu;
  if (!phase && delay) {
     FCVars.ph0= fmod(-delay*nu*360,360)*DEG2RAD;
  } else FCVars.ph0      = phase*DEG2RAD;
  FCVars.sum_t=FCVars.sum_v=FCVars.sum_N=FCVars.sum_N_pass=0;

  /* check of input parameters */
  if (nslit < 1) nslit=1;
  if (yheight <= 0) exit(printf("FermiChopper: %s: FATAL: unrealistic cylinder yheight =%g [m]\n", NAME_CURRENT_COMP, yheight));

  if (m <= 0) { m=0; R0=0; }
  if (radius <= 0) {
    printf("FermiChopper: %s: FATAL: Unrealistic cylinder radius radius=%g [m]\n", NAME_CURRENT_COMP, radius);
    exit(-1);
  }
  if (xwidth > 0 && xwidth < radius*2 && nslit > 0) {
    w = xwidth/nslit;
  }
  if (w <= 0) {
    printf("FermiChopper: %s: FATAL: Slits in the package have unrealistic width w=%g [m]\n", NAME_CURRENT_COMP, w);
    exit(-1);
  }
  if (nslit*w > radius*2) {
    nslit = floor(radius/w);
    printf("FermiChopper: %s: Too many slits to fit in the cylinder\n"
           "    Adjusting nslit=%f\n", NAME_CURRENT_COMP, nslit);
  }
  if (length > radius*2) {
    length = 2*sqrt(radius*radius - nslit*w*nslit*w/4);
    printf("FermiChopper: %s: Slit package is longer than the whole\n"
           "    chopper cylinder. Adjusting length=%g [m]\n", NAME_CURRENT_COMP, length);
  }

  if (eff <= 0 || eff > 1) {
    eff = 0.95;
    printf("FermiChopper: %s: Efficiency is unrealistic\n"
           "    Adjusting eff=%f\n", NAME_CURRENT_COMP, eff);
  }
  if (Qc <= 0) { Qc = 0.02176; m = 0; R0 = 0; }
  if (W <= 0) W=1e-6;

  if (curvature) {
    FCVars.C_slit = curvature;
    if (1 < fabs(radius*curvature))
      exit(printf("FermiChopper: %s: Slit curvature is unrealistic\n",
           NAME_CURRENT_COMP));
  }
  FCVars.L_slit = length;
  if (verbose && nu)
    printf("FermiChopper: %s: Frequency nu=%g [Hz] %g [rpm], time frame=%g [s] phase=%g [deg]\n"
      , NAME_CURRENT_COMP, nu, nu*60, 2/nu, FCVars.ph0*RAD2DEG);

  FCVars.absorb_alreadyinside    = 0;
  FCVars.absorb_topbottom        = 0;
  FCVars.absorb_cylentrance      = 0;
  FCVars.absorb_sideentrance     = 0;
  FCVars.absorb_notreachentrance = 0;
  FCVars.absorb_packentrance     = 0;
  FCVars.absorb_slitcoating      = 0;
  FCVars.warn_notreachslitwall   = 0;
  FCVars.absorb_exitslitpack     = 0;
  FCVars.absorb_maxiterations    = 0;
  FCVars.absorb_wrongdirection   = 0;
  FCVars.absorb_nocontrol        = 0;
  FCVars.absorb_cylexit          = 0;
  FCVars.warn_notreachslitoutput = 0;
  
  /* fix for the wrong coordinate frame orientation to come back to McStas XYZ system */
  FCVars.omega *= -1;
  FCVars.ph0   *= -1;
  FCVars.t0     = -FCVars.ph0/FCVars.omega;
}
#line 16468 "./Test_Fermi.c"
#undef delay
#undef curvature
#undef yheight
#undef verbose
#undef xwidth
#undef zero_time
#undef eff
#undef length
#undef W
#undef m
#undef alpha
#undef Qc
#undef R0
#undef nslit
#undef w
#undef nu
#undef radius
#undef phase
#undef FCVars
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component FC_ILL. */
  SIG_MESSAGE("FC_ILL (Init)");
#define mccompcurname  FC_ILL
#define mccompcurtype  FermiChopper_ILL
#define mccompcurindex 8
#define FCVars mccFC_ILL_FCVars
#define phase mccFC_ILL_phase
#define radius mccFC_ILL_radius
#define nu mccFC_ILL_nu
#define yheight mccFC_ILL_yheight
#define w mccFC_ILL_w
#define nslit mccFC_ILL_nslit
#define R0 mccFC_ILL_R0
#define Qc mccFC_ILL_Qc
#define alpha mccFC_ILL_alpha
#define m mccFC_ILL_m
#define W mccFC_ILL_W
#define length mccFC_ILL_length
#define eff mccFC_ILL_eff
#define zero_time mccFC_ILL_zero_time
#define xwidth mccFC_ILL_xwidth
#define verbose mccFC_ILL_verbose
#line 241 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../contrib/FermiChopper_ILL.comp"
{

/************************* INITIALIZE COUNTERS ******************************/

  int i;

/************************ CALCULATION CONSTANTS *****************************/
  FCVars.omega    = 2*PI*nu;
  if (nu && phase) FCVars.t0 = -phase/360.0/nu;

  /* check of input parameters */
  if (m < 0) m == 0;
  if (radius <= 0) {
    printf("FermiChopper_ILL: %s: FATAL: unrealistic cylinder radius radius=%g [m]\n", NAME_CURRENT_COMP, radius);
    exit(-1);
  }
  if (yheight <= 0) 
  	exit(printf("FermiChopper_ILL: %s: FATAL: unrealistic cylinder yheight =%g [m]\n", NAME_CURRENT_COMP, yheight));
  if (xwidth > 0 && xwidth < radius*2 && nslit > 0) {
    w = xwidth/nslit;
  }
  if (w <= 0) {
    printf("FermiChopper_ILL: %s: FATAL: Slits in the package have unrealistic width w=%g [m]\n", NAME_CURRENT_COMP, w);
    exit(-1);
  }
  if (nslit*w > radius*2) {
    nslit = floor(radius/w);
    printf("FermiChopper_ILL: %s: Too many slits to fit in the cylinder\n"
           "Adjusting nslit=%f\n", NAME_CURRENT_COMP, nslit);
  }
  if (length > radius*2) {
    length = sqrt(radius*radius - nslit*w*nslit*w/4);
    printf("FermiChopper_ILL: %s: Slit package is longer than the whole\n"
           "chopper cylinder. Adjusting length=%g [m]\n", NAME_CURRENT_COMP, length);
  }

  if (eff <= 0 || eff > 1) {
    eff = 0.95;
    printf("FermiChopper_ILL: %s: Efficiency is unrealistic\n"
           "Adjusting eff=%f\n", NAME_CURRENT_COMP, eff);
  }
  if (Qc <= 0) { Qc = 0.02176; m = 0; R0=0; }
  if (W <= 0) W=1e-6;
  
  if (verbose && nu)
    printf("FermiChopper_ILL: %s: frequency nu=%g [Hz] %g [rpm], time frame=%g [s] phase=%g [deg]\n"
      , NAME_CURRENT_COMP, nu, nu*60, 2/nu, -FCVars.t0*360*nu);
  
  /* fix for the wrong coordinate frame orientation to come back to McStas XYZ system */
  FCVars.omega *= -1;
}
#line 16566 "./Test_Fermi.c"
#undef verbose
#undef xwidth
#undef zero_time
#undef eff
#undef length
#undef W
#undef m
#undef alpha
#undef Qc
#undef R0
#undef nslit
#undef w
#undef yheight
#undef nu
#undef radius
#undef phase
#undef FCVars
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component Fake_Origin. */
  SIG_MESSAGE("Fake_Origin (Init)");

  /* Initializations for component FC_Vitess. */
  SIG_MESSAGE("FC_Vitess (Init)");
#define mccompcurname  FC_Vitess
#define mccompcurtype  Vitess_ChopperFermi
#define mccompcurindex 10
#define Option mccFC_Vitess_Option
#define CurvGeomOption mccFC_Vitess_CurvGeomOption
#define TOF mccFC_Vitess_TOF
#define WL mccFC_Vitess_WL
#define radius_of_curv mccFC_Vitess_radius_of_curv
#define main_depth mccFC_Vitess_main_depth
#define shift_y mccFC_Vitess_shift_y
#define angle_channel mccFC_Vitess_angle_channel
#define phase0 mccFC_Vitess_phase0
#define y_ch mccFC_Vitess_y_ch
#define x_ch mccFC_Vitess_x_ch
#define coef_pi mccFC_Vitess_coef_pi
#define XFILEName mccFC_Vitess_XFILEName
#define GeomFilePtr mccFC_Vitess_GeomFilePtr
#define Pos mccFC_Vitess_Pos
#define Dir mccFC_Vitess_Dir
#define Neutrons mccFC_Vitess_Neutrons
#define pos_ch mccFC_Vitess_pos_ch
#define omega mccFC_Vitess_omega
#define optimal_wl mccFC_Vitess_optimal_wl
#define sGeomFileName mccFC_Vitess_sGeomFileName
#define GeomOption mccFC_Vitess_GeomOption
#define zerotime mccFC_Vitess_zerotime
#define Nchannels mccFC_Vitess_Nchannels
#define Ngates mccFC_Vitess_Ngates
#define freq mccFC_Vitess_freq
#define height mccFC_Vitess_height
#define width mccFC_Vitess_width
#define depth mccFC_Vitess_depth
#define r_curv mccFC_Vitess_r_curv
#define diameter mccFC_Vitess_diameter
#define Phase mccFC_Vitess_Phase
#define wallwidth mccFC_Vitess_wallwidth
#line 149 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Vitess_ChopperFermi.comp"
{
omega = 0.0;
optimal_wl = 0.0;

double x,y,z;
coords_get(POS_R_COMP_INDEX(INDEX_CURRENT_COMP), &x, &y, &z);
  pos_ch[0] = -100.0 * z;
  pos_ch[1] = -100.0 * x;
  pos_ch[2] = -100.0 * y;

  McInitVt();

  /* transformation of units */
  height    *= 100.0; /* m -> cm */
  width     *= 100.0;
  depth     *= 100.0;
  diameter  *= 100.0;
  r_curv    *= 100.0;
  wallwidth *= 100.0;
  omega      = freq*2*PI/1000.0;  /* 1/s -> 2pi/ms */
  Phase     *= DEG2RAD;

  /* checks and completion of input data */
  CurvGeomOption = (int) GeomOption;
  if (GeomOption > 0 && omega*r_curv==0)
  { 
    printf("Error: 'omega*r_curv' must not be zero for curved Fermi chopper"); exit(-1);
  }
  switch(GeomOption)
  {
    case 0: Option=1; optimal_wl=0.0;    break;
    case 1: Option=2; optimal_wl=LAMBDA_FROM_V(2.0*omega*r_curv); break;
    case 2: Option=2; optimal_wl=LAMBDA_FROM_V(2.0*omega*r_curv); break;
    default: printf("Wrong option! Good options: 0-straight, 1-parabolic, 2-circular");
  }

  ChopperFermiInit(0, NULL);
}
#line 16668 "./Test_Fermi.c"
#undef wallwidth
#undef Phase
#undef diameter
#undef r_curv
#undef depth
#undef width
#undef height
#undef freq
#undef Ngates
#undef Nchannels
#undef zerotime
#undef GeomOption
#undef sGeomFileName
#undef optimal_wl
#undef omega
#undef pos_ch
#undef Neutrons
#undef Dir
#undef Pos
#undef GeomFilePtr
#undef XFILEName
#undef coef_pi
#undef x_ch
#undef y_ch
#undef phase0
#undef angle_channel
#undef shift_y
#undef main_depth
#undef radius_of_curv
#undef WL
#undef TOF
#undef CurvGeomOption
#undef Option
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component Monitor2_xt. */
  SIG_MESSAGE("Monitor2_xt (Init)");
#define mccompcurname  Monitor2_xt
#define mccompcurtype  Monitor_nD
#define mccompcurindex 11
#define user1 mccMonitor2_xt_user1
#define user2 mccMonitor2_xt_user2
#define user3 mccMonitor2_xt_user3
#define DEFS mccMonitor2_xt_DEFS
#define Vars mccMonitor2_xt_Vars
#define detector mccMonitor2_xt_detector
#define offdata mccMonitor2_xt_offdata
#define xwidth mccMonitor2_xt_xwidth
#define yheight mccMonitor2_xt_yheight
#define zdepth mccMonitor2_xt_zdepth
#define xmin mccMonitor2_xt_xmin
#define xmax mccMonitor2_xt_xmax
#define ymin mccMonitor2_xt_ymin
#define ymax mccMonitor2_xt_ymax
#define zmin mccMonitor2_xt_zmin
#define zmax mccMonitor2_xt_zmax
#define bins mccMonitor2_xt_bins
#define min mccMonitor2_xt_min
#define max mccMonitor2_xt_max
#define restore_neutron mccMonitor2_xt_restore_neutron
#define radius mccMonitor2_xt_radius
#define options mccMonitor2_xt_options
#define filename mccMonitor2_xt_filename
#define geometry mccMonitor2_xt_geometry
#define username1 mccMonitor2_xt_username1
#define username2 mccMonitor2_xt_username2
#define username3 mccMonitor2_xt_username3
#define nowritefile mccMonitor2_xt_nowritefile
#line 231 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/Monitor_nD.comp"
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
#line 16819 "./Test_Fermi.c"
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
#line 17005 "./Test_Fermi.c"
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
#define mccompcurtype  Source_gen
#define mccompcurindex 2
#define p_in mccSource_p_in
#define lambda1 mccSource_lambda1
#define lambda2 mccSource_lambda2
#define lambda3 mccSource_lambda3
#define pTable mccSource_pTable
#define pTable_x mccSource_pTable_x
#define pTable_y mccSource_pTable_y
#define pTable_xmin mccSource_pTable_xmin
#define pTable_xmax mccSource_pTable_xmax
#define pTable_xsum mccSource_pTable_xsum
#define pTable_ymin mccSource_pTable_ymin
#define pTable_ymax mccSource_pTable_ymax
#define pTable_ysum mccSource_pTable_ysum
#define pTable_dxmin mccSource_pTable_dxmin
#define pTable_dxmax mccSource_pTable_dxmax
#define pTable_dymin mccSource_pTable_dymin
#define pTable_dymax mccSource_pTable_dymax
{   /* Declarations of Source=Source_gen() SETTING parameters. */
char* flux_file = mccSource_flux_file;
char* xdiv_file = mccSource_xdiv_file;
char* ydiv_file = mccSource_ydiv_file;
MCNUM radius = mccSource_radius;
MCNUM dist = mccSource_dist;
MCNUM focus_xw = mccSource_focus_xw;
MCNUM focus_yh = mccSource_focus_yh;
MCNUM focus_aw = mccSource_focus_aw;
MCNUM focus_ah = mccSource_focus_ah;
MCNUM E0 = mccSource_E0;
MCNUM dE = mccSource_dE;
MCNUM lambda0 = mccSource_lambda0;
MCNUM dlambda = mccSource_dlambda;
MCNUM I1 = mccSource_I1;
MCNUM yheight = mccSource_yheight;
MCNUM xwidth = mccSource_xwidth;
MCNUM verbose = mccSource_verbose;
MCNUM T1 = mccSource_T1;
MCNUM flux_file_perAA = mccSource_flux_file_perAA;
MCNUM flux_file_log = mccSource_flux_file_log;
MCNUM Lmin = mccSource_Lmin;
MCNUM Lmax = mccSource_Lmax;
MCNUM Emin = mccSource_Emin;
MCNUM Emax = mccSource_Emax;
MCNUM T2 = mccSource_T2;
MCNUM I2 = mccSource_I2;
MCNUM T3 = mccSource_T3;
MCNUM I3 = mccSource_I3;
MCNUM zdepth = mccSource_zdepth;
int target_index = mccSource_target_index;
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
#line 17253 "./Test_Fermi.c"
/* 'Source=Source_gen()' component instance extend code */
    SIG_MESSAGE("Source (Trace:Extend)");
#line 94 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  /* 1 ms triangle time window */
  // t = randtriangle()*1e-3;  /* trianglular distribution */
  // t = randpm1()*time_window_width-time_to_arrival;          /* rectangular distribution */
  t = randpm1()*time_window_width;
  vx=vy=0;
#line 17262 "./Test_Fermi.c"
}   /* End of Source=Source_gen() SETTING parameter declarations. */
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

  /* TRACE Component Monitor1_xt [3] */
  mccoordschange(mcposrMonitor1_xt, mcrotrMonitor1_xt,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component Monitor1_xt (without coords transformations) */
  mcJumpTrace_Monitor1_xt:
  SIG_MESSAGE("Monitor1_xt (Trace)");
  mcDEBUG_COMP("Monitor1_xt")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompMonitor1_xt
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
#define mccompcurname  Monitor1_xt
#define mccompcurtype  Monitor_nD
#define mccompcurindex 3
#define user1 mccMonitor1_xt_user1
#define user2 mccMonitor1_xt_user2
#define user3 mccMonitor1_xt_user3
#define DEFS mccMonitor1_xt_DEFS
#define Vars mccMonitor1_xt_Vars
#define detector mccMonitor1_xt_detector
#define offdata mccMonitor1_xt_offdata
{   /* Declarations of Monitor1_xt=Monitor_nD() SETTING parameters. */
MCNUM xwidth = mccMonitor1_xt_xwidth;
MCNUM yheight = mccMonitor1_xt_yheight;
MCNUM zdepth = mccMonitor1_xt_zdepth;
MCNUM xmin = mccMonitor1_xt_xmin;
MCNUM xmax = mccMonitor1_xt_xmax;
MCNUM ymin = mccMonitor1_xt_ymin;
MCNUM ymax = mccMonitor1_xt_ymax;
MCNUM zmin = mccMonitor1_xt_zmin;
MCNUM zmax = mccMonitor1_xt_zmax;
MCNUM bins = mccMonitor1_xt_bins;
MCNUM min = mccMonitor1_xt_min;
MCNUM max = mccMonitor1_xt_max;
MCNUM restore_neutron = mccMonitor1_xt_restore_neutron;
MCNUM radius = mccMonitor1_xt_radius;
char* options = mccMonitor1_xt_options;
char* filename = mccMonitor1_xt_filename;
char* geometry = mccMonitor1_xt_geometry;
char* username1 = mccMonitor1_xt_username1;
char* username2 = mccMonitor1_xt_username2;
char* username3 = mccMonitor1_xt_username3;
int nowritefile = mccMonitor1_xt_nowritefile;
#line 312 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/Monitor_nD.comp"
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
#line 17582 "./Test_Fermi.c"
}   /* End of Monitor1_xt=Monitor_nD() SETTING parameter declarations. */
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
  mcabsorbCompMonitor1_xt:
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

  /* TRACE Component FC_Position [4] */
  mccoordschange(mcposrFC_Position, mcrotrFC_Position,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component FC_Position (without coords transformations) */
  mcJumpTrace_FC_Position:
  SIG_MESSAGE("FC_Position (Trace)");
  mcDEBUG_COMP("FC_Position")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompFC_Position
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
#define mccompcurname  FC_Position
#define mccompcurtype  Arm
#define mccompcurindex 4
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompFC_Position:
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

  /* TRACE Component FC_GuideG [5] */
  mccoordschange(mcposrFC_GuideG, mcrotrFC_GuideG,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component FC_GuideG (without coords transformations) */
  mcJumpTrace_FC_GuideG:
  SIG_MESSAGE("FC_GuideG (Trace)");
  mcDEBUG_COMP("FC_GuideG")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompFC_GuideG
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
#define mccompcurname  FC_GuideG
#define mccompcurtype  Guide_gravity
#define mccompcurindex 5
#define GVars mccFC_GuideG_GVars
#define pTable mccFC_GuideG_pTable
{   /* Declarations of FC_GuideG=Guide_gravity() SETTING parameters. */
MCNUM w1 = mccFC_GuideG_w1;
MCNUM h1 = mccFC_GuideG_h1;
MCNUM w2 = mccFC_GuideG_w2;
MCNUM h2 = mccFC_GuideG_h2;
MCNUM l = mccFC_GuideG_l;
MCNUM R0 = mccFC_GuideG_R0;
MCNUM Qc = mccFC_GuideG_Qc;
MCNUM alpha = mccFC_GuideG_alpha;
MCNUM m = mccFC_GuideG_m;
MCNUM W = mccFC_GuideG_W;
MCNUM nslit = mccFC_GuideG_nslit;
MCNUM d = mccFC_GuideG_d;
MCNUM mleft = mccFC_GuideG_mleft;
MCNUM mright = mccFC_GuideG_mright;
MCNUM mtop = mccFC_GuideG_mtop;
MCNUM mbottom = mccFC_GuideG_mbottom;
MCNUM nhslit = mccFC_GuideG_nhslit;
MCNUM G = mccFC_GuideG_G;
MCNUM aleft = mccFC_GuideG_aleft;
MCNUM aright = mccFC_GuideG_aright;
MCNUM atop = mccFC_GuideG_atop;
MCNUM abottom = mccFC_GuideG_abottom;
MCNUM wavy = mccFC_GuideG_wavy;
MCNUM wavy_z = mccFC_GuideG_wavy_z;
MCNUM wavy_tb = mccFC_GuideG_wavy_tb;
MCNUM wavy_lr = mccFC_GuideG_wavy_lr;
MCNUM chamfers = mccFC_GuideG_chamfers;
MCNUM chamfers_z = mccFC_GuideG_chamfers_z;
MCNUM chamfers_lr = mccFC_GuideG_chamfers_lr;
MCNUM chamfers_tb = mccFC_GuideG_chamfers_tb;
MCNUM nelements = mccFC_GuideG_nelements;
MCNUM nu = mccFC_GuideG_nu;
MCNUM phase = mccFC_GuideG_phase;
char* reflect = mccFC_GuideG_reflect;
/* 'FC_GuideG=Guide_gravity()' component instance has conditional execution */
if (( mcipFermi == 4 ))

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
#line 18005 "./Test_Fermi.c"
/* 'FC_GuideG=Guide_gravity()' component instance extend code */
    SIG_MESSAGE("FC_GuideG (Trace:Extend)");
if (( mcipFermi == 4 )) {

#line 113 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  if (!SCATTERED) ABSORB;
#line 18011 "./Test_Fermi.c"
}

}   /* End of FC_GuideG=Guide_gravity() SETTING parameter declarations. */
#undef pTable
#undef GVars
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompFC_GuideG:
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

  /* TRACE Component FC_GuideC [6] */
  mccoordschange(mcposrFC_GuideC, mcrotrFC_GuideC,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component FC_GuideC (without coords transformations) */
  mcJumpTrace_FC_GuideC:
  SIG_MESSAGE("FC_GuideC (Trace)");
  mcDEBUG_COMP("FC_GuideC")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompFC_GuideC
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
#define mccompcurname  FC_GuideC
#define mccompcurtype  Guide_channeled
#define mccompcurindex 6
#define w1c mccFC_GuideC_w1c
#define w2c mccFC_GuideC_w2c
#define ww mccFC_GuideC_ww
#define hh mccFC_GuideC_hh
#define whalf mccFC_GuideC_whalf
#define hhalf mccFC_GuideC_hhalf
#define lwhalf mccFC_GuideC_lwhalf
#define lhhalf mccFC_GuideC_lhhalf
{   /* Declarations of FC_GuideC=Guide_channeled() SETTING parameters. */
MCNUM w1 = mccFC_GuideC_w1;
MCNUM h1 = mccFC_GuideC_h1;
MCNUM w2 = mccFC_GuideC_w2;
MCNUM h2 = mccFC_GuideC_h2;
MCNUM l = mccFC_GuideC_l;
MCNUM R0 = mccFC_GuideC_R0;
MCNUM Qc = mccFC_GuideC_Qc;
MCNUM alpha = mccFC_GuideC_alpha;
MCNUM m = mccFC_GuideC_m;
MCNUM nslit = mccFC_GuideC_nslit;
MCNUM d = mccFC_GuideC_d;
MCNUM Qcx = mccFC_GuideC_Qcx;
MCNUM Qcy = mccFC_GuideC_Qcy;
MCNUM alphax = mccFC_GuideC_alphax;
MCNUM alphay = mccFC_GuideC_alphay;
MCNUM W = mccFC_GuideC_W;
MCNUM mx = mccFC_GuideC_mx;
MCNUM my = mccFC_GuideC_my;
MCNUM nu = mccFC_GuideC_nu;
MCNUM phase = mccFC_GuideC_phase;
/* 'FC_GuideC=Guide_channeled()' component instance has conditional execution */
if (( mcipFermi == 5 ))

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
#line 18315 "./Test_Fermi.c"
/* 'FC_GuideC=Guide_channeled()' component instance extend code */
    SIG_MESSAGE("FC_GuideC (Trace:Extend)");
if (( mcipFermi == 5 )) {

#line 121 "/zhome/89/0/38697/TESTS/2019-11-20/McStas-2.5_CPU_MPICC/Test_Fermi/Test_Fermi.instr"
  if (!SCATTERED) ABSORB;
#line 18321 "./Test_Fermi.c"
}

}   /* End of FC_GuideC=Guide_channeled() SETTING parameter declarations. */
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
  mcabsorbCompFC_GuideC:
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

  /* TRACE Component FC_McStas [7] */
  mccoordschange(mcposrFC_McStas, mcrotrFC_McStas,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component FC_McStas (without coords transformations) */
  mcJumpTrace_FC_McStas:
  SIG_MESSAGE("FC_McStas (Trace)");
  mcDEBUG_COMP("FC_McStas")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompFC_McStas
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
#define mccompcurname  FC_McStas
#define mccompcurtype  FermiChopper
#define mccompcurindex 7
#define FCVars mccFC_McStas_FCVars
{   /* Declarations of FC_McStas=FermiChopper() SETTING parameters. */
MCNUM phase = mccFC_McStas_phase;
MCNUM radius = mccFC_McStas_radius;
MCNUM nu = mccFC_McStas_nu;
MCNUM w = mccFC_McStas_w;
MCNUM nslit = mccFC_McStas_nslit;
MCNUM R0 = mccFC_McStas_R0;
MCNUM Qc = mccFC_McStas_Qc;
MCNUM alpha = mccFC_McStas_alpha;
MCNUM m = mccFC_McStas_m;
MCNUM W = mccFC_McStas_W;
MCNUM length = mccFC_McStas_length;
MCNUM eff = mccFC_McStas_eff;
MCNUM zero_time = mccFC_McStas_zero_time;
MCNUM xwidth = mccFC_McStas_xwidth;
MCNUM verbose = mccFC_McStas_verbose;
MCNUM yheight = mccFC_McStas_yheight;
MCNUM curvature = mccFC_McStas_curvature;
MCNUM delay = mccFC_McStas_delay;
/* 'FC_McStas=FermiChopper()' component instance has conditional execution */
if (( mcipFermi == 1 ))

#line 361 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/FermiChopper.comp"
{

  /** local CALCULATION VARIABLES**************************************/

  /** Interaction with slit package ***************************/
  double slit_input; /* length of the slits */

  /** Variables for calculating interaction with blades ***************/
  double xp1, zp1, xp2, zp2, vxp1, vxp2, xp3, vzp1;

  /**  Reflections ***********************************************/
  double n1,n2,n3;

  /** variables used for calculating new velocities after reflection **/
  double q;

  /**  Multiple Reflections  ******************************/
  int loopcounter=0;   /* How many reflections happen? */

  /** Time variables *********************************/
  double t3;      /* interaction time at n3 position (slit wall) */
  double t1,t2;   /* cylinder intersection time (at entry and exit of slit pack - n1 n2 - or cylinder) */
  double dt;      /* interaction intervals (e,g, time for exiting slit pack) */

/************************ TIME OF FLIGHT RESET ************************/
  if (zero_time == 1 && nu)
    t -= floor( (t+1/(4*nu))*(2*nu) )/(2*nu);

/************** test, if the neutron interacts with the cylinder ***/
  if (cylinder_intersect (&t1, &t2, x, y, z, vx, vy, vz, radius, yheight))
  {
    if (t1 <= 0) { /* Neutron started inside the cylinder */
      if (verbose > 2 && FCVars.absorb_alreadyinside<FermiChopper_MAXITER) printf("FermiChopper: %s: ABSORB Neutron started inside the cylinder, t1=%8.3g (enter).\n", 
        NAME_CURRENT_COMP, t1);
      FCVars.absorb_alreadyinside++; 
      ABSORB; 
    }    

    dt=t2-t1;     /* total time of flight inside the cylinder  */
    PROP_DT(t1);  /* Propagates neutron to entrance of the cylinder */
    SCATTER;
    
    if (verbose > 2)
      printf("FermiChopper: %s: PROP_DT t1=%8.3g t2=%8.3g xyz=[%8.3g %8.3g %8.3g] v=[%8.3g %8.3g %8.3g] t=%8.3g (IN cyl).\n",
           NAME_CURRENT_COMP, t1, t2, x,y,z,vx,vy,vz,t);

    /* neutron must not enter or leave from top or bottom of cylinder. */
    if (fabs(y) >= yheight/2.0 || fabs(y+vy*dt) >= yheight/2.0) { 
      if (verbose > 2) printf("FermiChopper: %s: ABSORB Neutron hits top/bottom of cylinder, y=%8.3g (enter).\n",
        NAME_CURRENT_COMP, y);
      FCVars.absorb_topbottom++; ABSORB; 
    }

    vxp1 = sqrt(vx*vx+vy*vy+vz*vz);
    FCVars.sum_v += p*vxp1;
    FCVars.sum_t += p*t;
    FCVars.sum_N += p;

    if (zero_time > 1 && FCVars.sum_N) { /* automatic phase adjustment */
      double mean_t, mean_phase;
      mean_t     = FCVars.sum_t/FCVars.sum_N;
      mean_phase = fmod(mean_t*nu*2*PI, 2*PI);
      /* now we shift the phase so that the neutron pulse is centered on the slit pack */
      mean_phase+=radius/vxp1*2*PI*nu;
      FCVars.ph0 = mean_phase;
    }

    /* neutron must enter the cylinder opening: |X'| < full package width*/
    xp1 = FC_xrot(x,z, t,FCVars); /* X'(t) */
    if (fabs(xp1) >= nslit*w/2) { 
      if (verbose > 2) printf("FermiChopper: %s: ABSORB Neutron X is outside cylinder aperture, xp1=%8.3g (enter).\n", 
        NAME_CURRENT_COMP, xp1);
      FCVars.absorb_cylentrance++; ABSORB; 
    }

/*********************** PROPAGATE TO SLIT PACKAGE **************************/

    /* zp1 = Z' at entrance of cylinder Z'(t) */
    zp1  = FC_zrot(x,z, t, FCVars);

    /* Checking on which side of the Chopper the Neutron enters: sign(Z') */
    slit_input = (zp1 > 0 ? length/2 : -length/2);

    /* time shift to reach slit package in [0,time to exit cylinder]: Z'=slit_input */
    /* t3 is used here as a tmp variable, will be redefined in for loop  */
    t3 = FC_zintersect(x,z,vx,vz, t,dt, slit_input, FCVars);

    if( (t3 < 0)||(t3 > dt) ) {
      if (verbose > 1 && FCVars.absorb_notreachentrance < FermiChopper_MAXITER) {
        printf("FermiChopper: %s: Can not reach entrance of slits. dt=%8.3g t3=%8.3g.\n",
        NAME_CURRENT_COMP, dt, t3);
        if (t3 == -3)
            printf("          No sign change to determine intersection (intersection:1)\n");
        else if (t3 == -2)
            printf("          Max iterations reached (intersection:1)\n");
        else if (t3 < 0)
            printf("          Error when solving intersection (intersection:1)\n");
      }
      FCVars.absorb_notreachentrance++;
      ABSORB; /* neutron can not reach slit entrance */
    }

    /* Propagating to the slit package entrance */
    PROP_DT(t3);
    dt -= t3; /* remaining time from slit pack entry to exit of cylinder */
    xp1 = FC_xrot(x,z, t, FCVars); /* should be slit_input */
    zp1 = FC_zrot(x,z, t, FCVars);
    if (mcdotrace) {
      /* indicate position of neutron in mcdisplay */
      xp2 = x; zp2 = z; x = xp1; z=zp1; SCATTER; x=xp2; z=zp2;
    } else SCATTER;
    
    if (verbose > 2)
      printf("FermiChopper: %s: PROP_DT t3=%8.3g dt=%8.3g xyz=[%8.3g %8.3g %8.3g] length=%g (slit enter) z'=%g.\n",
           NAME_CURRENT_COMP, t3, dt, x,y,z, slit_input, zp1);

    /* must have X'< slit package width at package Z'=slit_input */
    if (fabs(xp1) >= nslit*w/2) { 
      if (verbose > 2) printf("FermiChopper: %s: ABSORB Neutron X is outside slit package, xp1=%8.3g (enter).\n", 
        NAME_CURRENT_COMP, xp1);
      FCVars.absorb_packentrance++; ABSORB; 
    }
    
    /* solve Z'=-slit_input for time of exit of slit package */
    /* t3 is used here as a tmp variable, will be redefined in for loop  */
    t3 = FC_zintersect(x,z,vx,vz, t,dt, -slit_input, FCVars);

    if((t3 < FermiChopper_TimeAccuracy)||(t3 > dt)) {
      if (verbose > 1 && FCVars.warn_notreachslitoutput < FermiChopper_MAXITER) {
        printf("FermiChopper: %s: Can not reach exit of slits. dt=%8.3g t3=%8.3g.\n",
          NAME_CURRENT_COMP, dt, t3);
        if (t3 == -3)
          printf("              No sign change to determine intersection (intersection:2)\n");
        else if (t3 == -2)
          printf("              Max iterations reached (intersection:2)\n");
        else if (t3 < 0)
          printf("              Error when solving intersection (intersection:2)\n");
      }
      FCVars.warn_notreachslitoutput++;
      ABSORB;
      /* neutron can not reach slit output. */
    } else dt = t3; /* reduce time interval to [0, time of slit exit] */
    
    /* dt= time shift to go to exit of slit package (or exit of cylinder in case of error) */
    /*
      |---------|
      | /       |
      |o        * (dt)
      |         |
      |---------|
     */

/*********************PROPAGATION INSIDE THE SLIT PACKAGE *******************/

    /* index of slit hit at entrance n1=-N/2 (xp1=-) ... N/2 (xp1=+) */
    n1 = floor(xp1/w);

/******************* BEGIN LOOP FOR MULTIPLE REFLECTIONS ********************/

    for (loopcounter=0; loopcounter<=FermiChopper_MAXITER; loopcounter++) {
      double dt_to_tangent=0; /* time shift to go to tangent intersection */

      /* compute trajectory tangents: m1=Vz'+w.X'(t), m2=Vz'+w.X'(t+dt) */
      xp1 = FC_xrot    (x,z,          t,   FCVars);          /* X'(t)    current position */
      xp2 = FC_xzrot_dt(x,z,vx,vz,    t,   dt, 'x', FCVars); /* X'(t+dt) slit exit */

      /* slit index at the end of the slit: */
      n2 = floor(xp2/w);

      vxp1= FC_xrot    (vx+z*FCVars.omega,vz-x*FCVars.omega,
                                      t,   FCVars);          /* dX'(t)/dt slope at current position*/

      vxp2= FC_xrot    (vx+(z+vz*dt)*FCVars.omega,vz-(x+vx*dt)*FCVars.omega,
                                      t+dt,FCVars);          /* dX'(t+dt)/dt slope at slit exit */

      /* absolute time at tangent intersection, changed to time shift below */
      dt_to_tangent = (vxp1 - vxp2 ? ((xp2 - xp1) + (t*vxp1 - (t+dt)*vxp2))/(vxp1 - vxp2) - t
                                   : -1);

      /* If algorithm fails, take the middle of the interval*/
      if((dt_to_tangent < 0)||(dt_to_tangent > dt)) dt_to_tangent=dt*0.5; 

      /*
           *(dt_to_tangent)
      |---------|
      | /     \ |
      |o       \|(dt)
      |         |
      |---------|
     */

      /* point coordinates at tangent intersection/middle point */
      xp3 = FC_xzrot_dt(x,z,vx,vz, t, dt_to_tangent, 'x', FCVars); /* X'(t+dt_to_tangent) */

      /* slit index at the tangent intersection/middle point */
      n3 = floor(xp3/w);
      
      if (verbose > 2)
        printf("FermiChopper: %s: t3=%8.3g n=[%g %g %g] (time at tangent intersection).\n",
           NAME_CURRENT_COMP, dt_to_tangent, n1, n2, n3);

      /* change slit means there is a reflection/intersection inside */
      if ( (n2!=n1) || (n3!=n1) ) {

        double t3a, t3b, distance_Wa, distance_Wb;
        if (m == 0 || R0 == 0) { 
          if (verbose > 2) printf("FermiChopper: %s: ABSORB Neutron hits absorbing coating (change slit).\n", 
            NAME_CURRENT_COMP);
          FCVars.absorb_slitcoating++; ABSORB; 
        }
        
        /* Choosing the first time it isn't in the slit anymore */
        if(n3!=n1){
          n2 = n3;
        }
        
        /* get position of slit wall towards which neutron is propagating */
        if (n2 > n1) {     /* X' positive side of slit is in principle the first intersection to test*/
          distance_Wa = n1*w+w;
          distance_Wb = n1*w;
        } else {            /* X' negative side of slit */
          distance_Wb = n1*w+w;
          distance_Wa = n1*w;
        }

        /* time shift to reach slit wall point in [0,dt_to_tangent]: X'=distance_W slit_wall */
        t3a = FC_xintersect(x,z,vx,vz, t,dt_to_tangent, distance_Wa, FCVars);
        t3b = FC_xintersect(x,z,vx,vz, t,dt_to_tangent, distance_Wb, FCVars);
        if (t3b < 0) t3 = t3a;
        else if (t3a < 0 && t3b >= 0) t3 = t3b;
        else t3 = (t3a < t3b ? t3a : t3b);

        /* handle case where intersection search fails */
        if((t3 < FermiChopper_TimeAccuracy)||(t3 >= dt_to_tangent)) {
          if (verbose > 1 && FCVars.warn_notreachslitwall < FermiChopper_MAXITER) {
            printf("FermiChopper: %s: Can not reach slit wall (iteration %i). dt=%8.3g t3=%8.3g.\n",
              NAME_CURRENT_COMP, loopcounter, dt_to_tangent, t3);
            if (t3 == -3)
              printf("        No sign change to determine intersection (intersection:3)\n");
            else if (t3 == -2)
              printf("        Max iterations reached (intersection:3)\n");
            else if (t3 < 0 || t3 >= dt_to_tangent)
              printf("        Error when solving intersection (intersection:3)\n");
          }
          /* we try to find a solution in a larger time range */
          dt_to_tangent = dt;
          t3a = FC_xintersect(x,z,vx,vz, t,dt_to_tangent, distance_Wa, FCVars);
          t3b = FC_xintersect(x,z,vx,vz, t,dt_to_tangent, distance_Wb, FCVars);
          if (t3b < 0) t3 = t3a;
          else if (t3a < 0 && t3b >= 0) t3 = t3b;
          else t3 = (t3a < t3b ? t3a : t3b);
        }
        if((t3 < FermiChopper_TimeAccuracy)||(t3 >= dt_to_tangent)) {
          if (verbose > 1 && FCVars.warn_notreachslitwall < FermiChopper_MAXITER) {
            printf("FermiChopper: %s: Can not reach slit wall (iteration %i). dt=%8.3g t3=%8.3g.\n",
              NAME_CURRENT_COMP, loopcounter, dt_to_tangent, t3);
            if (t3 == -3)
              printf("        No sign change to determine intersection (intersection:4)\n");
            else if (t3 == -2)
              printf("        Max iterations reached (intersection:4)\n");
            else if (t3 < 0 || t3 >= dt_to_tangent)
              printf("        Error when solving intersection (intersection:4)\n");
          }
          /* neutron can not reach slit wall. */
          FCVars.warn_notreachslitwall++;
          ABSORB;
        }

        /* Propagate to slit wall point (t+t3) on slit n3 wall */
        PROP_DT(t3); dt -= t3; /* dt: time remaining to slit exit after propagation */
        zp1 = FC_zrot(x,z, t, FCVars); /* Z'(t+t3) */
        
        if (verbose > 2)
          printf("FermiChopper: %s: PROP_DT t3=%8.3g dt=%8.3g xyz=[%8.3g %8.3g %8.3g] (on wall). z'=%g\n",
           NAME_CURRENT_COMP, t3, dt, x,y,z, zp1);
        
        /* check if intersection point is still in the slit package, else exit loop */
        if (fabs(zp1) > length/2) {
          if (verbose > 2) 
            printf("FermiChopper: %s: Neutron is outside slit pack (on slit wall).\n", 
                NAME_CURRENT_COMP);
          break;
        }
        
    /*
      |---o-----|
      | /   \   |
      |/     \  |
      |       \ |(dt)
      |---------|
     */

        /* get velocity in rotating frame, on slit wall */
        vxp1 = FC_xrot(vx,vz, t, FCVars);
        vzp1 = FC_zrot(vx,vz, t, FCVars);

        q    = 2*V2Q*(fabs(vxp1));

        {
          double ref = 0;
          double par[] = {R0, Qc, alpha, m, W};
          StdReflecFunc(q, par, &ref);
          if (ref>0) p *= ref;
          else { 
            if (verbose > 2) printf("FermiChopper: %s: ABSORB Neutron hits absorbing coating (on slit wall).\n", 
              NAME_CURRENT_COMP);
            FCVars.absorb_slitcoating++; ABSORB; 
          } /* Cutoff ~ 1E-10 */
        }
        if (mcdotrace) {
          xp1 = FC_xrot(x,z, t,FCVars);
          zp1 = FC_zrot(x,z, t,FCVars);
          /* indicate position of neutron in mcdisplay */
          xp2 = x; zp2 = z; x = xp1; z=zp1; SCATTER; x=xp2; z=zp2;
        } else SCATTER;

        /* reflect perpendicular velocity and compute new velocity in static frame */
        vxp1 *= -1;
        /* apply transposed Transformation matrix */
        vx = FC_xrot( vxp1,-vzp1, t,FCVars);
        vz = FC_zrot(-vxp1, vzp1, t,FCVars);
        
        /* recompute time to slit exit */
        /* solve Z'=-slit_input for time of exit of slit package */
        t3 = FC_zintersect(x,z,vx,vz, t,dt, -slit_input, FCVars);

        if((t3 < 0)||(t3 > dt)) {
          if (verbose > 1 && FCVars.warn_notreachslitoutput < FermiChopper_MAXITER) {
            printf("FermiChopper: %s: Can not reach exit of slits. dt=%8.3g t3=%8.3g.\n",
              NAME_CURRENT_COMP, dt, t3);
            if (t3 == -3)
              printf("              No sign change to determine intersection (intersection:2)\n");
            else if (t3 == -2)
              printf("              Max iterations reached (intersection:2)\n");
            else if (t3 < 0)
              printf("              Error when solving intersection (intersection:2)\n");
          }
          FCVars.warn_notreachslitoutput++;
          ABSORB;
          /* neutron can not reach slit output. */
        } else dt = t3; /* reduce time interval to [0, time of slit exit] */

      } /* end if changed slit index */
      else { /* neutron remains in the same slit: direct flight */
        break;
      }
    } /* end for */

    if (loopcounter >= FermiChopper_MAXITER) {
      if (verbose > 1 && FCVars.absorb_maxiterations < FermiChopper_MAXITER)
      printf("FermiChopper: %s: Max iterations %i reached inside slit. Absorb.\n",
        NAME_CURRENT_COMP, FermiChopper_MAXITER);
      FCVars.absorb_maxiterations++;
      ABSORB;
    }


    /********************* EXIT SLIT PACKAGE ********************************/

    /* propagation times to cylinder. Should use t2 to exit */
    if (!cylinder_intersect (&t1, &t2, x, y, z, vx, vy, vz, radius, yheight)) { 
      /* must be inside the cylinder */
      if (verbose > 2) printf("FermiChopper: %s: ABSORB Neutron has unexpectidely exited cylinder ! (exiting)\n", 
        NAME_CURRENT_COMP);
      FCVars.absorb_exitslitpack++;  
      ABSORB; } 

    if (t1 > 0) {
      if (verbose > 1 && FCVars.absorb_wrongdirection < FermiChopper_MAXITER)
      printf("FermiChopper: %s: Neutrons are leaving chopper\n"
             "              in the wrong direction. Absorb.\n", NAME_CURRENT_COMP);
      FCVars.absorb_wrongdirection++;
      ABSORB;
    }

    if (t2 <= 0 && FCVars.absorb_nocontrol < FermiChopper_MAXITER) {
      if (verbose > 1)
      printf("FermiChopper: %s: Neutrons are leaving chopper\n"
             "              without any control. Absorb.\n", NAME_CURRENT_COMP);
      FCVars.absorb_nocontrol++;
      ABSORB;
    }

    /* propagate to cylinder surface (exit) */
    PROP_DT(t2);
    SCATTER;
    
    xp1 = FC_xrot(x,z, t,FCVars);
    zp1 = FC_zrot(x,z, t,FCVars);
    
    if (verbose > 2)
      printf("FermiChopper: %s: t1=%8.3g PROP_DT t2=%8.3g xyz=[%8.3g %8.3g %8.3g] (OUT cyl). z'=%g\n",
           NAME_CURRENT_COMP, t1, t2, x,y,z, zp1);

    /* Check if the neutron left the cylinder by its top or bottom */
    if (fabs(y) >= yheight/2) { 
      if (verbose > 2) printf("FermiChopper: %s: ABSORB Neutron hits top/bottom of cylinder, y=%8.3g (exiting)\n", 
        NAME_CURRENT_COMP, y);
      FCVars.absorb_topbottom++; ABSORB; 
    }


    /* must have X'< slit package width at package Z'=cylinder output */
    if (fabs(xp1) >= nslit*w/2) { 
      if (verbose > 2) printf("FermiChopper: %s: ABSORB Neutron X is outside slit package cylinder, xp1=%8.3g (exiting).\n", 
        NAME_CURRENT_COMP, xp1);
      FCVars.absorb_cylexit++; ABSORB; 
    }

    /* Transmission coefficent */
    p = p*eff;          //finite cross section + transmission

    FCVars.sum_N_pass += p;

  } /* end if cylinder_intersect */
  if (!SCATTERED) {
    if (verbose > 2 && 0) printf("FermiChopper: %s: ABSORB Neutron has not interacted with FC\n", 
        NAME_CURRENT_COMP);
    ABSORB;
  }

}
#line 18880 "./Test_Fermi.c"
}   /* End of FC_McStas=FermiChopper() SETTING parameter declarations. */
#undef FCVars
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompFC_McStas:
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

  /* TRACE Component FC_ILL [8] */
  mccoordschange(mcposrFC_ILL, mcrotrFC_ILL,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component FC_ILL (without coords transformations) */
  mcJumpTrace_FC_ILL:
  SIG_MESSAGE("FC_ILL (Trace)");
  mcDEBUG_COMP("FC_ILL")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompFC_ILL
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
#define mccompcurname  FC_ILL
#define mccompcurtype  FermiChopper_ILL
#define mccompcurindex 8
#define FCVars mccFC_ILL_FCVars
{   /* Declarations of FC_ILL=FermiChopper_ILL() SETTING parameters. */
MCNUM phase = mccFC_ILL_phase;
MCNUM radius = mccFC_ILL_radius;
MCNUM nu = mccFC_ILL_nu;
MCNUM yheight = mccFC_ILL_yheight;
MCNUM w = mccFC_ILL_w;
MCNUM nslit = mccFC_ILL_nslit;
MCNUM R0 = mccFC_ILL_R0;
MCNUM Qc = mccFC_ILL_Qc;
MCNUM alpha = mccFC_ILL_alpha;
MCNUM m = mccFC_ILL_m;
MCNUM W = mccFC_ILL_W;
MCNUM length = mccFC_ILL_length;
MCNUM eff = mccFC_ILL_eff;
MCNUM zero_time = mccFC_ILL_zero_time;
MCNUM xwidth = mccFC_ILL_xwidth;
MCNUM verbose = mccFC_ILL_verbose;
/* 'FC_ILL=FermiChopper_ILL()' component instance has conditional execution */
if (( mcipFermi == 3 ))

#line 294 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../contrib/FermiChopper_ILL.comp"
{

  /** local CALCULATION VARIABLES**************************************/

  /** Interaction with slitpacket ***************************/
  double slit_input;   /* length of the slits */
  double zr1,zr2;       /* distance to slitpacket entrance/exit in rotating frame */
  double xr1,xr2;       /* X entrance/exit position in rotating frame */

  /** Variables for calculating interaction with blades ***************/
  double m1,m2;    /* slope of the tangents */
  double b1,b2;    /* y-intersection of tangent */

  /**  Reflections ***********************************************/
  double t3a, t3b, distance_Wa, distance_Wb;
  double n1,n2,n3,n4;

  /** variables used for calculating new velocities after reflection **/
  double q;
  double vper, vpar;
  double arg;

  /**  Multiple Reflections  ******************************/
  int loopcounter=0;   /* How many reflections happen? */

  /** Time variables *********************************/
  double t3;      /* interaction time */
  double dt;      /* interaction intervals */
  double t1,t2;   /* cylinder intersection time */


/************** test, if the neutron interacts with the cylinder ***/
  if (cylinder_intersect (&t1, &t2, x, y, z, vx, vy, vz, radius, yheight)) {

    if (t1 <= 0) {
      if (verbose > 2) 
      	printf("FermiChopper_ILL: %s: ABSORB Neutron started inside the cylinder, t1=%g (enter)\n", 
          NAME_CURRENT_COMP, t1);
      ABSORB;    /* Neutron started inside the cylinder */
    }

    dt=t2-t1;     /* total time of flight inside the cylinder  */
    PROP_DT(t1);  /* Propagates neutron to entrance of the cylinder */
    
    if (verbose > 2)
      printf("FermiChopper_ILL: %s: PROP_DT t1=%8.3g t2=%8.3g xyz=[%8.3g %8.3g %8.3g] v=[%8.3g %8.3g %8.3g] t=%8.3g (IN cyl).\n",
             NAME_CURRENT_COMP, t1, t2, x,y,z,vx,vy,vz,t);

    if(dt > fabs(0.5/FCVars.omega*2*PI) && verbose) {
      printf("FermiChopper_ILL: %s: Frequency too low. Method will fail.\n"
             "              Absorbing neutron\n", NAME_CURRENT_COMP);
      ABSORB;
    }

  /* Checks if neutron enters or leaves from top or bottom of cylinder. */
    if ( fabs(y) > yheight/2 ||
        fabs(y+vy*dt) > yheight/2 ) {
      if (verbose > 2) 
      	printf("FermiChopper_ILL: %s: ABSORB Neutron hits top/bottom of cylinder, y=%8.3g (enter)\n", 
          NAME_CURRENT_COMP, y);
      ABSORB;
    }

  /*  checking wether the neutron can enter the chopper (slit channel) */
    xr1 = xstrich(x,z,t, FCVars.omega, FCVars.t0);
    if(fabs(xr1)>=nslit*w/2) {
      if (verbose > 2) 
      	printf("FermiChopper_ILL: %s: ABSORB Neutron X is outside cylinder aperture, xp1=%8.3g (enter)\n", 
          NAME_CURRENT_COMP, xr1);
      ABSORB;
    }

  /*********************** PROPAGATE TO SLIT PACKAGE **************************/


    /* Checking on which side of the Chopper the Neutron enters******/
    slit_input = 0.5*length;
    zr1 = zstrich(x,z,t, FCVars.omega, FCVars.t0);
    zr2 = zstrich(x+vx*dt,z+vz*dt,t+dt, FCVars.omega, FCVars.t0);
    if(zr1 < 0) slit_input *= -1;

    /*  Checking if the Neutron will hit the slits (Z) */
    zr1  -=  slit_input;
    zr2  -=  slit_input;

    if (zr2*zr1>0) {
      if (verbose > 2) 
      	printf("FermiChopper_ILL: %s: ABSORB Neutron Z does not change sign, zr1=%8.3g zr2=%8.3g (enter)\n", 
          NAME_CURRENT_COMP, zr1,zr2);
      ABSORB;
    }

    /* Calculating where/when Neutron hits the slits (Z) */
    t3 = zsecant(x,z,vx,vz,t,dt,slit_input, FCVars.omega, FCVars.t0);

    if((t3 < 0)||(t3 > dt)) {
      t3 = zinterpolation(x,z,vx,vz,t,dt,slit_input, FCVars.omega, FCVars.t0);
    }
    if((t3 < 0)||(t3 > dt)) {
      if (verbose) 
      	printf("FermiChopper_ILL: %s: Can not reach entrance of slits. dt=%g t3=%g\n", NAME_CURRENT_COMP, dt, t3);
      ABSORB;
    }

    /* Propagating whole system to that point */
    PROP_DT(t3);
    dt -= t3;
    SCATTER;
    
    if (verbose > 2)
      printf("FermiChopper_ILL: %s: PROP_DT t3=%8.3g dt=%8.3g xyz=[%8.3g %8.3g %8.3g] length=%g (slit enter).\n",
             NAME_CURRENT_COMP, t3, dt, x,y,z, slit_input);

    /* Checking if neutron hits the slits entrance window (X) */
    xr1 = xstrich(x,z,t, FCVars.omega, FCVars.t0);
    if(fabs(xr1) >= nslit*w/2) {
      if (verbose > 2) 
      	printf("FermiChopper_ILL: %s: ABSORB Neutron X is outside slit package, xp1=%8.3g (enter)\n", 
        NAME_CURRENT_COMP, xr1);
      ABSORB;
    }

    /* Calculating where/when Neutron leaves the slits (Z) */
    t3 = zsecant(x,z,vx,vz,t,dt,-slit_input, FCVars.omega, FCVars.t0);
    if((t3 < 0) || (t3 > dt)){
      t3 = zinterpolation(x,z,vx,vz,t,dt,-slit_input, FCVars.omega, FCVars.t0);
    }
    if((t3 <= 0) || (t3 > dt)){
      if (verbose) 
      	printf("FermiChopper_ILL: %s: Can not reach exit of slits. dt=%8.3g t3=%8.3g\n", NAME_CURRENT_COMP, dt, t3);
      ABSORB;
    } else dt=t3;

  /********************* PROPAGATION INSIDE THE SLIT PACKET *******************/

    /* Which slit was hit ? */
    n1 = floor(xr1/w);


  /******************* BEGIN LOOP FOR MULTIPLE REFLECTIONS ********************/

    for(loopcounter; loopcounter<=FCILL_MAXITER;loopcounter++){
    	double dt_to_tangent=0; /* time shift to go to tangent intersection */

  /* Calculate most probable time for interaction with blades by using tangents */
      m1 = xstrich(vx,vz,t, FCVars.omega, FCVars.t0)
         + FCVars.omega * zstrich(x,z,t, FCVars.omega, FCVars.t0);
      m2 = xstrich(vx,vz,t+dt, FCVars.omega, FCVars.t0)
         + FCVars.omega * zstrich(x+vx*dt,z+vz*dt,t+dt,FCVars.omega,FCVars.t0);

      b1 = xstrich(x,z,t, FCVars.omega, FCVars.t0) - m1*t;
      b2 = xstrich(x+vx*dt,z+vz*dt,t+dt, FCVars.omega, FCVars.t0) - m2*(t+dt);

      if (m1-m2) dt_to_tangent = ((b2-b1)/(m1-m2))-t;
      else       dt_to_tangent = -1;

      /* If method with tangents doesn't succeed, just take the middle of the interval */
      if((dt_to_tangent < 0)||(dt_to_tangent > dt)) dt_to_tangent=dt*0.5;

      /* Calculate different positions for the neutron to determine interaction. */

      /*...at the end of the slit: */
      n2 = floor(xstrich(x+(vx*dt),z+(vz*dt),t+dt, FCVars.omega, FCVars.t0)/w);

      /*...at the before calculated t3: tangent intersection point */
      n3 = floor(xstrich(x+(vx*dt_to_tangent),z+(vz*dt_to_tangent),t+dt_to_tangent, FCVars.omega, FCVars.t0)/w);

      if (verbose > 2)
        printf("FermiChopper_ILL: %s: t3=%8.3g n=[%g %g %g] (time at tangent intersection).\n",
             NAME_CURRENT_COMP, dt_to_tangent, n1, n2, n3);

      /* Does the neutron stay in the same slit ? */
      if((n2!=n1)||(n3!=n1)){

        /* Choosing the first time it isn't in the slit anymore */
        if(n3!=n1){
          n2 = n3;
        }

      /************ABSORB to save calculation time ******************/
        if (m == 0 || R0 == 0) {
          if (verbose > 2) 
          	printf("FermiChopper_ILL: %s: ABSORB Neutron hits absorbing coating (change slit).\n", 
              NAME_CURRENT_COMP);
          ABSORB;
        }


  /********************** WHEN DOES IT HIT THE BLADE? *************************/

        /*********** SECANT METHOD ****************************/
        /* get position of slit wall towards which neutron is propagating */
        if (n2 > n1) {     /* X' positive side of slit is in principle the first intersection to test*/
          distance_Wa = n1*w+w;
          distance_Wb = n1*w;
        } else {           /* X' negative side of slit */
        	distance_Wb = n1*w+w;
          distance_Wa = n1*w;
        }

        /* time shift to reach slit wall point in [0,dt_to_tangent]: X'=distance_W slit_wall */
        t3a = xsecant(x,z,vx,vz,t,dt,distance_Wa, FCVars.omega, FCVars.t0);
        t3b = xsecant(x,z,vx,vz,t,dt,distance_Wb, FCVars.omega, FCVars.t0);
        if (t3b < 0) t3 = t3a;
        else if (t3a < 0 && t3b > 0) t3 = t3b;
        else t3 = (t3a < t3b ? t3a : t3b);

        /***** INTERPOLATION USED WHEN SECANT METHOD FAILS ****/
        /* try second intersection method in case of failure */
        if ((t3 < 0) || (t3 > dt)) {
        	t3a = xinterpolation(x,z,vx,vz,t,dt,distance_Wa, FCVars.omega, FCVars.t0);
        	t3b = xinterpolation(x,z,vx,vz,t,dt,distance_Wb, FCVars.omega, FCVars.t0);
        	if (t3b < 0) t3 = t3a;
		      else if (t3a < 0 && t3b > 0) t3 = t3b;
		      else t3 = (t3a < t3b ? t3a : t3b);
        }

        /* Check for errors in calculation*******/
        if ((t3 < 0) || (t3 > dt)) {
          if (verbose) 
          	printf("FermiChopper_ILL: %s: Reflecting interpolation Problem. dt=%8.3g t3=%8.3g\n", 
          	NAME_CURRENT_COMP, dt, t3);
          ABSORB;
        }
        /* Propagate whole system to that point */
        PROP_DT(t3); dt -= t3;
        if (verbose > 2)
          printf("FermiChopper_ILL: %s: PROP_DT t3=%8.3g dt=%8.3g xyz=[%8.3g %8.3g %8.3g] (on wall).\n",
             NAME_CURRENT_COMP, t3, dt, x,y,z);

        /* Check if this point is inside the slit packet */
        if(fabs(zstrich(x,z,t, FCVars.omega, FCVars.t0)) > fabs(slit_input)){
          if (verbose > 2) 
          	printf("FermiChopper_ILL: %s: Neutron is outside slit pack (on slit wall).\n", 
                NAME_CURRENT_COMP);
          break;
        }

  /******************** REFLECTION ALGORITHM ********************************/
        vper    = xstrich(vx,vz,t, FCVars.omega, FCVars.t0); /* perpendicular velocity (to blade) */
        vpar    = zstrich(vx,vz,t, FCVars.omega, FCVars.t0);  /* parallel velocity (to blade) */
        q       = 2*MS2AA*(fabs(vper));

        if (q > Qc && W){
          arg = (q-m*Qc)/W;
          if (arg < 10.0) p *= 0.5*(1-tanh(arg))*(1-alpha*(q-Qc));
          else {
            if (verbose > 2) 
            	printf("FermiChopper_ILL: %s: ABSORB Neutron hits absorbing coating (on slit wall).\n", 
                NAME_CURRENT_COMP);
            ABSORB;
          }
        }

        if (R0 != 0.0){
          p *= R0;

          vper *= (-1);   /* Mirroring perpendicular velocity */

          /**************SET NEW VELOCITIES***********/
          vx =  vper*cos(FCVars.omega*(t-FCVars.t0))
             -  vpar*sin(FCVars.omega*(t-FCVars.t0));
          vz =  vper*sin(FCVars.omega*(t-FCVars.t0))
             +  vpar*cos(FCVars.omega*(t-FCVars.t0));
          SCATTER;
        } else {
          if (verbose > 2) 
          	printf("FermiChopper_ILL: %s: ABSORB Neutron hits absorbing coating (R0=0).\n", 
              NAME_CURRENT_COMP);
         ABSORB;
        }


        /* Recalculating when Neutron will leave the slitpacket */
        t3 = zsecant(x,z,vx,vz,t,dt,-slit_input,FCVars.omega,FCVars.t0);
        if((t3 < 0) || (t3 > dt)) {
          t3=zinterpolation(x,z,vx,vz,t,dt,-slit_input,
                             FCVars.omega,FCVars.t0);
        }
        /* Check for errors in calculation*******/
        if ((t3 < 0) || (t3 > dt)) {
          if (verbose) 
          	printf("FermiChopper_ILL: %s: Reflecting interpolation Problem. dt=%8.3g t3=%8.3g\n", 
          	NAME_CURRENT_COMP, dt, t3);
          ABSORB;
        } else  dt=t3;
      } /* end if n2 != n2 != n3 */
      else break;
    } /* end for */
  /********************* END OF THE FOR LOOP **********************************/

    /****New time of cylinder intersection will be calculated**********/
    if (!cylinder_intersect (&t1, &t2, x, y, z, vx, vy, vz, radius, yheight)) {
    	if (verbose > 2) 
    		printf("FermiChopper_ILL: %s: ABSORB Neutron has unexpectidely exited cylinder ! (exiting)\n", 
    			NAME_CURRENT_COMP);
    	ABSORB; 
    }

    if (t1 > 0 && verbose) {
      printf("FermiChopper_ILL: %s: Neutrons are leaving chopper in the wrong direction! \n", NAME_CURRENT_COMP);
    }

    if (t2 <= 0 && verbose) {
      printf("FermiChopper_ILL: %s: Neutrons are leaving chopper without any control\n", NAME_CURRENT_COMP);
    }

  /*********** PROPAGATE TO CYLINDER SURFACE ***********************************/
    PROP_DT(t2);
    SCATTER;
    
    if (verbose > 2)
      printf("FermiChopper_ILL: %s: t1=%8.3g PROP_DT t2=%8.3g xyz=[%8.3g %8.3g %8.3g] (OUT cyl).\n",
             NAME_CURRENT_COMP, t1, t2, x,y,z);

    /*****Checking if the neutron left the cylinder by his top or bottom **/
    if  ( fabs(y) > yheight/2 ){
      if (verbose > 2) 
      	printf("FermiChopper_ILL: %s: ABSORB Neutron hits top/bottom of cylinder, y=%8.3g (exiting)\n", 
          NAME_CURRENT_COMP, y);
      ABSORB;
    }


    /*****Checking if neutron hits chopper exit ***/
    if(fabs(xstrich(x,z,t,FCVars.omega,FCVars.t0))>=nslit*w/2){
      if (verbose > 2) 
      	printf("FermiChopper_ILL: %s: ABSORB Neutron X is outside slit package cylinder, xp1=%8.3g (exiting)\n", 
          NAME_CURRENT_COMP, xstrich(x,z,t,FCVars.omega,FCVars.t0));
      ABSORB;
    }

    /**** Transmission coefficent******/
    p = p*eff;          //finite cross section + transmission

  } /* end if cylinder_intersect */
  else {
    if (verbose > 2 && 0) 
    	printf("FermiChopper_ILL: %s: ABSORB Neutron has not interacted with FC\n", 
        NAME_CURRENT_COMP);
    ABSORB;
  }

/************************ TIME OF FLIGHT RESET ************************/
  if (zero_time && nu)
    t -= (((int)((t+1/(4*nu))/(1/(2*nu))))*(1/(2*nu)));
}
#line 19354 "./Test_Fermi.c"
}   /* End of FC_ILL=FermiChopper_ILL() SETTING parameter declarations. */
#undef FCVars
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompFC_ILL:
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

  /* TRACE Component Fake_Origin [9] */
  mccoordschange(mcposrFake_Origin, mcrotrFake_Origin,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component Fake_Origin (without coords transformations) */
  mcJumpTrace_Fake_Origin:
  SIG_MESSAGE("Fake_Origin (Trace)");
  mcDEBUG_COMP("Fake_Origin")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompFake_Origin
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
#define mccompcurname  Fake_Origin
#define mccompcurtype  Arm
#define mccompcurindex 9
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompFake_Origin:
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

  /* TRACE Component FC_Vitess [10] */
  mccoordschange(mcposrFC_Vitess, mcrotrFC_Vitess,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component FC_Vitess (without coords transformations) */
  mcJumpTrace_FC_Vitess:
  SIG_MESSAGE("FC_Vitess (Trace)");
  mcDEBUG_COMP("FC_Vitess")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompFC_Vitess
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
#define mccompcurname  FC_Vitess
#define mccompcurtype  Vitess_ChopperFermi
#define mccompcurindex 10
#define Option mccFC_Vitess_Option
#define CurvGeomOption mccFC_Vitess_CurvGeomOption
#define TOF mccFC_Vitess_TOF
#define WL mccFC_Vitess_WL
#define radius_of_curv mccFC_Vitess_radius_of_curv
#define main_depth mccFC_Vitess_main_depth
#define shift_y mccFC_Vitess_shift_y
#define angle_channel mccFC_Vitess_angle_channel
#define phase0 mccFC_Vitess_phase0
#define y_ch mccFC_Vitess_y_ch
#define x_ch mccFC_Vitess_x_ch
#define coef_pi mccFC_Vitess_coef_pi
#define XFILEName mccFC_Vitess_XFILEName
#define GeomFilePtr mccFC_Vitess_GeomFilePtr
#define Pos mccFC_Vitess_Pos
#define Dir mccFC_Vitess_Dir
#define Neutrons mccFC_Vitess_Neutrons
#define pos_ch mccFC_Vitess_pos_ch
#define omega mccFC_Vitess_omega
#define optimal_wl mccFC_Vitess_optimal_wl
{   /* Declarations of FC_Vitess=Vitess_ChopperFermi() SETTING parameters. */
char* sGeomFileName = mccFC_Vitess_sGeomFileName;
int GeomOption = mccFC_Vitess_GeomOption;
int zerotime = mccFC_Vitess_zerotime;
int Nchannels = mccFC_Vitess_Nchannels;
int Ngates = mccFC_Vitess_Ngates;
MCNUM freq = mccFC_Vitess_freq;
MCNUM height = mccFC_Vitess_height;
MCNUM width = mccFC_Vitess_width;
MCNUM depth = mccFC_Vitess_depth;
MCNUM r_curv = mccFC_Vitess_r_curv;
MCNUM diameter = mccFC_Vitess_diameter;
MCNUM Phase = mccFC_Vitess_Phase;
MCNUM wallwidth = mccFC_Vitess_wallwidth;
/* 'FC_Vitess=Vitess_ChopperFermi()' component instance has conditional execution */
if (( mcipFermi == 2 ))

#line 189 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Vitess_ChopperFermi.comp"
{
 int i=0;
 InputNeutrons[i] = mcstas2vitess(x, y, z, vx, vy, vz, t, sx, sy, sz, p);
 #define  MCSTAS_TRACE
/********************************************************************************************************************************************************/
/*  VITESS module 'chopper_fermi'                                                            
/*                                                                                           
/* The free non-commercial use of these routines is granted providing due credit is given to 
/* the authors.                                                                              
/*                                                                                           
/* 1.00  Jul 2002  G. Zsigmond	initial version                                             
/* 1.01  Aug 2002  G. Zsigmond	forward all coordinates                                     
/* 1.02  Sep 2002  G. Zsigmond	included more channel windows                               
/* 1.03  Apr 2003  G. Zsigmond	sign correction                                             
/* 1.04  May 2003  G. Zsigmond	info changed                                                
/* 1.05  Jun 2003  G. Zsigmond	modulo function included for safety                         
/* 1.06  Jul 2003  G. Zsigmond	generalised for optional number of pulses; warnings included
/* 1.07  Oct 2003  G. Zsigmond	superfluous modulo function cancelled                       
/* 1.08  Nov 2003  G. Zsigmond	put 2 more windows representing channels, now 6 windows    
/*                               in the big IF loop ">=" changed to ">"                      
/* 1.09  Jan 2004  K. Lieutenant changes for 'instrument.dat'                                
/* 1.10  Jan 2004  G. Zsigmond   back to 4 windows representing channels
/* 1.11  Apr 2004  G. Zsigmond   negative time of flight defined
/* 1.12  Apr 2004  G. Zsigmond   negative time of flight - corrections, set zero time
/* 1.13  May 2004  G. Zsigmond   circular geom option and channel length included in curved fc
/* 1.14  JUL 2004  G. Zsigmond   small change in Init to adapt to new GUI
/* 1.15  OCT 2004  G. Zsigmond   changed to use both even or odd number of channels
/* 1.16  MAY 2005  G. Zsigmond   output changed to give trajectory coordinates at a plane crossing the center of the chopper (to be compatible with zero time option)
/*	                              zero time option fixed to get one peak 
/*                               shadowing cylinder opening activated 
/* 1.17  MAY 2005  G. Zsigmond  new option choice of 4, 6(better,slower) or 8(much better, very slow) gates, 4 gates option adjusted
/* 1.18  MAY 2005  K. Lieutenant changes to use this module in McStas as well               
/********************************************************************************************************************************************************/

#if defined VITESS

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "general.h"
#include "init.h"
#include "intersection.h"
#include "softabort.h"

/* START HEADER STORY */

/* input parameters */
int        Ngates=4,    /* Number of gates forming the channel: 4 (default), 6 or 8 */
           zerotime=0;     /* option: set time (close to) zero                    -z */
long       Nchannels;      /* number of channels of the Fermi chopper             -l */
double     omega,          /* frequency of rotation                       [1/s]   -n */
           height,         /* height of the Fermi chopper                  [cm]   -a */
           width,          /* width of the Fermi chopper                   [cm]   -b */
           depth,          /* channel length of the Fermi chopper          [cm]   -c */
           optimal_wl,     /* wavelength of highest transmission          [Ang]   -L */
           diameter,       /* diameter of the shadowing cylinder           [cm]   -r */
           Phase,          /* dephasing angle at zero time                [deg]   -l */
           wallwidth;      /* thickness of walls separating the channels          -m */
char       sGeomFileName[512];/* name of output file for geometry information        -G */
VectorType pos_ch;         /* centre position of the Fermi chopper        [cm] X,Y,V */

/* other global parameters */
#include "chopper_fermi.h"

/* prototypes and small functions */
void		OutputTransformations(double *tof, double *wl, double *prob, VectorType Pos, VectorType Dir, VectorType SpinVector);
void		ReadParameterFile();
void		ChopperFermiInit(int argc, char *argv[]);
void		ChopperFermiCleanup();
double	asinplus(double val);
double	asinminus(double val);

double		asin2PI(double val)
{ double result;
  if (val>=0) result = (double)asin(val);
  else  result =  2*M_PI + (double)asin(val);
  return result;
}

/* FINISH HEADER STORY */


int main(int argc, char **argv)
{
  long	 i;                 /* index of trajectory      */

  /* Initialize the program according to the parameters given  */

  Init(argc, argv, VT_CHOP_FERMI);
  print_module_name("Fermi-Chopper 1.17m");
  ChopperFermiInit(argc, argv);


  /* Get the neutrons from the file */
  DECLARE_ABORT;

  while((ReadNeutrons())!= 0)
    {
      CHECK;	/* here is what happens to the neutron */

      for(i=0;i<NumNeutGot;i++)

	{
	  CHECK;
#endif

#if defined VITESS || defined MCSTAS_TRACE
	  TOF = InputNeutrons[i].Time;

	  WL = InputNeutrons[i].Wavelength;

	  CopyVector(InputNeutrons[i].Position, Pos); 

	  CopyVector(InputNeutrons[i].Vector, Dir); 

	  Dir[0]	= (double)sqrt(1 - sq(Dir[1]) - sq(Dir[2])); 


	  /* shift to center of Fermi-Chopper */
#ifdef VITESS
	  SubVector(Pos, pos_ch);
#endif

		/*trajectories which do not intersect the entrance and exit window */

	  {double pos[3], n[3];

		  n[1]=n[2]=0.; n[0]=1.;
		  if((PlaneLineIntersect(Pos, Dir, n, - diameter/2., pos))==1)
			{
				  if((pos[2]>=height/2.)||(pos[2]<= - height/2.)||(pos[1]>=diameter/2.)||(pos[1]<= - diameter/2.)) ABSORB;
			}
		  else ABSORB;

		  if((PlaneLineIntersect(Pos, Dir, n, diameter/2., pos))==1)
			{
				  if((pos[2]>=height/2.)||(pos[2]<= - height/2.)||(pos[1]>=diameter/2.)||(pos[1]<= - diameter/2.)) ABSORB;
			}
		  else ABSORB;
	

		/* translates neutron variables for X'= - diameter/2.  */

		{
			VectorType Path;

			TOF = TOF + (- diameter/2. - Pos[0]) / fabs(Dir[0]) / V_FROM_LAMBDA(WL); 
			
			if((TOF<0)&&(Nchannels==1)){fprintf(LogFilePtr,"\nERROR: Single-slit Fermi chopper needs positive flight time at the chopper position! \n"); exit(-1); }
				
			CopyVector(Dir, Path);

			MultiplyByScalar(Path, (- diameter/2. - Pos[0])/ Dir[0] );

			AddVector(Pos, Path);  
		}							/*	 Path = displacement vector */

	
	/* calculate time entering-edge and exiting-edge of 4 windows along the channels */

	  {	long j, k, m;
	  double sq_D_0_1, sq_term, omega_fact, dirpos, vz_pos, phase[10][2000];
 
	  sq_D_0_1 = sq(Dir[0]) + sq(Dir[1]);
	  dirpos = Dir[0]*Pos[1] - Dir[1]*Pos[0];
	  sq_term  = sq(dirpos) / sq_D_0_1;
	  omega_fact = omega / (V_FROM_LAMBDA(WL) * Dir[0]);
	  vz_pos = Pos[1] > 0.0 ? 1.0 : -1.0;

	  phase0 = fmod(Phase + omega*TOF, coef_pi*M_PI); 
	
	  for(k=0; k<Ngates; k++) 
	  {
		for(j=0; j < 2*Nchannels+2; j++) {

		  double x_ch_k_j, y_ch_k_j, sq_x_ch_k_j, Denom_k, Arg_k, arg_k, pha_k_j, y_ch_new_k_j;

		  x_ch_k_j    = x_ch[k][j];
		  y_ch_k_j    = y_ch[k][j];
		  sq_x_ch_k_j = sq(x_ch_k_j);
		  Denom_k     = sqrt( sq_D_0_1 * (sq_x_ch_k_j + sq(y_ch_k_j)) );

		  Arg_k = dirpos / Denom_k;

		  if (fabs(Arg_k) > 1.) {
			Arg_k = vz_pos;
			y_ch_new_k_j = Arg_k * sqrt(sq_term - sq_x_ch_k_j);
		  } else
			y_ch_new_k_j = y_ch_k_j; /* no intersection with trajectory */

		  Denom_k = sqrt( sq_D_0_1 * (sq_x_ch_k_j + sq(y_ch_new_k_j)) );

		  arg_k = (Dir[0]*y_ch_new_k_j - Dir[1]*x_ch_k_j) / Denom_k;

		  if (fabs(arg_k) > 1.) {
			phase[k][j] = 777;} 
			  
		  else {
			pha_k_j = asin(Arg_k) - asin(arg_k); 

			if(x_ch_k_j < 0.) pha_k_j = - pha_k_j; 
								
			phase[k][j] =  pha_k_j - omega_fact * (x_ch_k_j * cos(pha_k_j) - y_ch_new_k_j * sin(pha_k_j) - Pos[0]);
		  }
		}
	  }

	  if(Ngates==4){
		  m = -1; 

		  for(j=0;j<2*Nchannels+1; j++) 
		  {
					if((m == 1)&&(phase[0][j] < phase0 )&&(phase0 < phase[0][j+1])
							   &&(phase[1][j] < phase0 )&&(phase0 < phase[1][j+1])
							   &&(phase[2][j] > phase0 )&&(phase0 > phase[2][j+1])
							   &&(phase[3][j] > phase0 )&&(phase0 > phase[3][j+1]))
					{/*fprintf(LogFilePtr, "j  %d   phases %f   %f    %f\n", j, 57.296*phase[0][j], 57.296*phase[0][j+1], 57.296*phase0);*/ 
											 goto happyend;}

					m = m * (-1); 
		  }
		  
		  
		  /* also tries one turn earlier  */
		  
		  if((phase0 > 0)&&(omega > 0)) phase0 +=  - coef_pi*M_PI;
		  if((phase0 < 0)&&(omega < 0)) phase0 +=  coef_pi*M_PI;

		  m = -1; 

		  for(j=0;j<2*Nchannels+1; j++) 
		  {
					if((m == 1)&&(phase[0][j] < phase0 )&&(phase0 < phase[0][j+1])
							   &&(phase[1][j] < phase0 )&&(phase0 < phase[1][j+1])
							   &&(phase[2][j] > phase0 )&&(phase0 > phase[2][j+1])
							   &&(phase[3][j] > phase0 )&&(phase0 > phase[3][j+1]))
					{/*fprintf(LogFilePtr, "j  %d   phases %f   %f    %f\n", j, 57.296*phase[0][j], 57.296*phase[0][j+1], 57.296*phase0);*/ 
											 goto happyend;}

					m = m * (-1); 
		  }
	  }

	  if(Ngates==6){
		  m = -1; 

		  for(j=0;j<2*Nchannels+1; j++) 
		  {
					if((m == 1)&&(phase[0][j] < phase0 )&&(phase0 < phase[0][j+1])
							   &&(phase[1][j] < phase0 )&&(phase0 < phase[1][j+1])
							   &&(phase[2][j] < phase0 )&&(phase0 < phase[2][j+1])
							   &&(phase[3][j] > phase0 )&&(phase0 > phase[3][j+1])
							   &&(phase[4][j] > phase0 )&&(phase0 > phase[4][j+1])
							   &&(phase[5][j] > phase0 )&&(phase0 > phase[5][j+1]))
					{goto happyend;}
											 
					m = m * (-1); 
		  }
		  
		  /* also tries one turn earlier  */
		  
		  if((phase0 > 0)&&(omega > 0)) phase0 +=  - coef_pi*M_PI;
		  if((phase0 < 0)&&(omega < 0)) phase0 +=  coef_pi*M_PI;

		  m = -1; 

		  for(j=0;j<2*Nchannels+1; j++) 
		  {
					if((m == 1)&&(phase[0][j] < phase0 )&&(phase0 < phase[0][j+1])
							   &&(phase[1][j] < phase0 )&&(phase0 < phase[1][j+1])
							   &&(phase[2][j] < phase0 )&&(phase0 < phase[2][j+1])
							   &&(phase[3][j] > phase0 )&&(phase0 > phase[3][j+1])
							   &&(phase[4][j] > phase0 )&&(phase0 > phase[4][j+1])
							   &&(phase[5][j] > phase0 )&&(phase0 > phase[5][j+1]))
					{goto happyend;}
											 
					m = m * (-1); 
		  }
	  }

	  if(Ngates==8){
		  m = -1; 

		  for(j=0;j<2*Nchannels+1; j++) 
		  {
					if((m == 1)&&(phase[0][j] < phase0 )&&(phase0 < phase[0][j+1])
							   &&(phase[1][j] < phase0 )&&(phase0 < phase[1][j+1])
							   &&(phase[2][j] < phase0 )&&(phase0 < phase[2][j+1])
							   &&(phase[3][j] < phase0 )&&(phase0 < phase[3][j+1])
							   &&(phase[4][j] > phase0 )&&(phase0 > phase[4][j+1])
							   &&(phase[5][j] > phase0 )&&(phase0 > phase[5][j+1])
							   &&(phase[6][j] > phase0 )&&(phase0 > phase[6][j+1])
							   &&(phase[7][j] > phase0 )&&(phase0 > phase[7][j+1]))
					{goto happyend;}
											 
					m = m * (-1); 
		  }
		  
		  /* also tries one turn earlier  */
		  
		  if((phase0 > 0)&&(omega > 0)) phase0 +=  - coef_pi*M_PI;
		  if((phase0 < 0)&&(omega < 0)) phase0 +=  coef_pi*M_PI;

		  m = -1; 

		  for(j=0;j<2*Nchannels+1; j++) 
		  {
					if((m == 1)&&(phase[0][j] < phase0 )&&(phase0 < phase[0][j+1])
							   &&(phase[1][j] < phase0 )&&(phase0 < phase[1][j+1])
							   &&(phase[2][j] < phase0 )&&(phase0 < phase[2][j+1])
							   &&(phase[3][j] < phase0 )&&(phase0 < phase[3][j+1])
							   &&(phase[4][j] > phase0 )&&(phase0 > phase[4][j+1])
							   &&(phase[5][j] > phase0 )&&(phase0 > phase[5][j+1])
							   &&(phase[6][j] > phase0 )&&(phase0 > phase[6][j+1])
							   &&(phase[7][j] > phase0 )&&(phase0 > phase[7][j+1]))
					{goto happyend;}
											 
					m = m * (-1); 
		  }
	  }

					ABSORB;

		  }
	happyend:;

	  /* Output matters */

	/* transmit coordinates which were not changed, the rest overwrite below */
	Neutrons = InputNeutrons[i]; 



	/* translates neutron variables for output - X'= 0. . */

	{
		VectorType Path;

		Neutrons.Time = TOF + (- Pos[0]) / Dir[0] / V_FROM_LAMBDA(WL); 
			
		if(zerotime==1)
		{ 
			Neutrons.Time = fabs(fmod(Neutrons.Time + Phase/omega + coef_pi*M_PI/omega/2., coef_pi*M_PI/omega)) - coef_pi*M_PI/2./omega ;
		}

		CopyVector(Dir, Path);

		MultiplyByScalar(Path, (- Pos[0])/ Dir[0] );

		AddVector(Pos, Path);  /*Path = displacement vector */
		
		CopyVector(Pos, Neutrons.Position);
	}												 

	}
#endif

#if defined VITESS 

	  /* writes output binary file */
	  WriteNeutron(&Neutrons);

	}
}

  /* Do the general cleanup */

my_exit:;

  ChopperFermiCleanup();

  fprintf(LogFilePtr," \n");

  Cleanup(pos_ch[0],pos_ch[1],pos_ch[2], 0.0,0.0);	

  return 0;
}
#endif



#if defined VITESS || defined MCSTAS_SHARE

/* own initialization of the module */

void ChopperFermiInit(int argc, char *argv[])
{
  /*    INPUT  */

 #ifdef VITESS
  GeomFileName="chopper_fermi_g.dat";

  CurvGeomOption = 1;
 #endif

  while(argc>1)
    {
      char *arg = &argv[1][2];

      switch(argv[1][1])
	{
	case 'X':
	  sscanf(arg, "%lf", &pos_ch[0]);
	  break;
	
	case 'Y':
	  sscanf(arg, "%lf", &pos_ch[1]);
	  break;

	case 'V':
	  sscanf(arg, "%lf", &pos_ch[2]);
	  break;

	case 'a':
	  sscanf(arg, "%lf", &height);
	  break;

	case 'b':
	  sscanf(arg, "%lf", &width);
	  break;

	case 'c':
	  sscanf(arg, "%lf", &depth); 
	  break;

	case 'L':
	  sscanf(arg, "%lf", &optimal_wl);
	  break;

	case 'l':
	  sscanf(arg, "%ld", &Nchannels);
	  break;

	case 'm':
	  sscanf(arg, "%lf", &wallwidth);
	  break;

	case 'n':
	  sscanf(arg, "%lf", &omega);
	  if(omega == 0) omega = 1.E-6;
	  omega = omega * 2. * M_PI / 1000.;
	  break;

	case 'q':
	  sscanf(arg, "%lf", &Phase);
	  Phase = Phase * M_PI / 180.;
	  break;

	case 'O':
	  sscanf(arg, "%d", &Option);
	  if((Option!=1)&&(Option!=2)){fprintf(LogFilePtr,"\nERROR: Wrong option! Good options: 1-straight, 2-curved \n"); exit(-1); }
	  break;

	case 'g':
	  sscanf(arg, "%d", &CurvGeomOption);
	  break;

	case 'r':
	  sscanf(arg, "%lf", &diameter);
	  break;

	case 'G':
	  GeomFileName=arg;
	  break;

	case 'z':
	  sscanf(arg, "%d", &zerotime);
	  break;


	}
      argc--;
      argv++;
    }

	if(pos_ch[0] < diameter/2.) {fprintf(LogFilePtr,"\nERROR: Minimum position is diameter/2 \n"); exit(-1); }

	if(Nchannels==1) wallwidth=0.;

	
  /* calculate edge positions */ 
  
  if(Option==1)
    {

      fprintf(LogFilePtr,"\nStraight Fermi chopper option activated\n");

      /* diameter matter */

      main_depth = 2. * sqrt(sq(diameter/2.) - sq(width/2.));
      if(depth > main_depth) {
	fprintf(LogFilePtr,"\nERROR: Diameter too small - not compatible with 'channel length' and 'width'!\nTake min %f cm\n", 2. * sqrt(sq(depth/2.) + sq(width/2.)));

	exit(-1);
      }

	{
	  double add, w_ch;

	  long j, k, m;

	  if((     w_ch = ( width - (Nchannels + 1) * wallwidth ) / Nchannels    ) <= 0.)
	    {fprintf(LogFilePtr,"\nERROR: Channel width =< 0 !\n"); exit(-1);}

		fprintf(LogFilePtr,"Channel width: %f cm\n",w_ch);

	  for(k=0;k<4; k++) y_ch[k][0] =  - width/2.;

	  m = 1;

	  for(j=1; j< 2*Nchannels+2; j++)
	    {

	      if(m == 1) add = wallwidth; else add = w_ch;

	      for(k=0;k<4; k++) 
		  {
				x_ch[k][j] = - depth/2. * (3.- k)/3. + k/3. * depth/2.;

				y_ch[k][j] = y_ch[k][j-1] + add;
		  }

	      m = m * (-1);
	    }

	  for(k=0;k<4; k++) 
	  {		  
		x_ch[k][0] = x_ch[k][1] = x_ch[k][2*Nchannels] = x_ch[k][2*Nchannels+1] = - depth/2. *(3. - k)/3. + k/3. * depth/2.;
	  }

	  /* activating shadowing cylinder */

		x_ch[0][0] = x_ch[0][1] = x_ch[0][2*Nchannels] = x_ch[0][2*Nchannels+1] = - main_depth/2.;
	
		x_ch[3][0] = x_ch[3][1] = x_ch[3][2*Nchannels] = x_ch[3][2*Nchannels+1] = main_depth/2.;
	
	}

    }

  if(Option==2)
    {	
		  if(CurvGeomOption==1)
		  {
						  fprintf(LogFilePtr,"\nCurved Fermi chopper activated \n");

						  fprintf(LogFilePtr,"\nGeometry option: ideally shaped (close to parabolic) long channels ('channel length' inactive parameter) \n");

						  fprintf(LogFilePtr,"\nRadius of curvature (parabolic approximation):\n" 
											 "	%f cm at center\n" 
											 "	%f cm at circumference\n", 
											 V_FROM_LAMBDA(optimal_wl)/2./omega, V_FROM_LAMBDA(optimal_wl)/2./omega*pow((1 + sq(omega*diameter/V_FROM_LAMBDA(optimal_wl))), 1.5));


				  GeomFilePtr = fopen(GeomFileName,"w");
				  {
						double add, w_ch, xx[500], yy[500], tt;

						long j, m, k;

						if((     w_ch = ( width - (Nchannels + 1) * wallwidth ) / Nchannels    ) <= 0.)
						  {fprintf(LogFilePtr,"\nChannel width =< 0 !\n"); exit(-1);}

							fprintf(LogFilePtr,"Channel width: %f cm\n",w_ch);

							y_ch[0][0] =  - width/2.;
							x_ch[0][0] = sqrt(sq(diameter/2.) - sq(y_ch[0][0]));
							x_ch[3][0] = diameter/2. * sin(omega*diameter/V_FROM_LAMBDA(optimal_wl) + acos(2.*y_ch[0][0]/diameter));
							y_ch[3][0] = diameter/2. * cos(omega*diameter/V_FROM_LAMBDA(optimal_wl) + acos(2.*y_ch[0][0]/diameter));

						for(k=0; k<500; k++)
						  { tt = k /500. * 2 * M_PI; xx[k] = diameter/2. * cos(tt); yy[k] = diameter/2. * sin(tt); /* the circle */
						  fprintf(GeomFilePtr, " %ld %lf  %lf\n", 0, xx[k], yy[k]);
						  }

						m = 1;

						for(j=1; j< 2*Nchannels+2; j++)
						{

							if(m == 1) add = wallwidth; else add = w_ch;
							
							y_ch[0][j] = y_ch[0][j-1] + add;

							x_ch[0][j] = - sqrt(sq(diameter/2.) - sq(y_ch[0][j]));

							
							for(k=0; k<500; k++)
							  {
								tt= k /500. * 2. * fabs(x_ch[0][j])/V_FROM_LAMBDA(optimal_wl);
								xx[k] = y_ch[0][j] * sin(omega*tt) + (V_FROM_LAMBDA(optimal_wl)*tt - fabs(x_ch[0][j])) * cos(omega*tt);
								yy[k] = y_ch[0][j] * cos(omega*tt) - (V_FROM_LAMBDA(optimal_wl)*tt - fabs(x_ch[0][j])) * sin(omega*tt);

								fprintf(GeomFilePtr, "%ld   %lf  %lf\n", j, xx[k], yy[k]);		
							  }

								x_ch[1][j] = xx[167];
								y_ch[1][j] = yy[167];
								x_ch[2][j] = xx[333];
								y_ch[2][j] = yy[333];
								x_ch[3][j] = diameter/2. * sin(2.* omega*fabs(x_ch[0][j])/V_FROM_LAMBDA(optimal_wl) + acos(2.*y_ch[0][j]/diameter));
								y_ch[3][j] = diameter/2. * cos(2.* omega*fabs(x_ch[0][j])/V_FROM_LAMBDA(optimal_wl) + acos(2.*y_ch[0][j]/diameter));

								m = m * (-1);
						}
				  }
		  }
		if(CurvGeomOption==2)
		{

		  
					/* diameter matter */

				  main_depth = 2. * sqrt(sq(diameter/2.) - sq(width/2.));
				  if(depth > main_depth) {
					fprintf(LogFilePtr,"\nERROR: Diameter too small - not compatible with 'channel length' and 'width'!\nTake min %f cm\n", 2. * sqrt(sq(depth/2.) + sq(width/2.))); exit(-1);}


		  
				  fprintf(LogFilePtr,"\nCurved Fermi chopper activated \n");

				  fprintf(LogFilePtr,"\nGeometry option: circular shaped channels with fixed length (via 'channel length') \n");

				  radius_of_curv = V_FROM_LAMBDA(optimal_wl)/2./omega;

				  angle_channel = atan(depth/2./radius_of_curv);

				  fprintf(LogFilePtr,"\nRadius of curvature (parabolic approximation):\n" 
									 "optimal_velocity/2./omega = %f cm \n" 
									 "Angle channel: %f deg \n", 
									 radius_of_curv, 180./M_PI*angle_channel);


			  GeomFilePtr = fopen(GeomFileName,"w");
			  {
					double add, w_ch, xx[500], yy[500], tt;

					long j, m, k;

					if((     w_ch = ( width - (Nchannels + 1) * wallwidth ) / Nchannels    ) <= 0.)
					  {fprintf(LogFilePtr,"\nChannel width =< 0 !\n"); exit(-1);}

					fprintf(LogFilePtr,"Channel width: %f cm\n",w_ch);


					for(k=0; k<500; k++)
					  { tt = k /500. * 2 * M_PI; xx[k] = diameter/2. * cos(tt); yy[k] = diameter/2. * sin(tt); /* the circle */
					  fprintf(GeomFilePtr, " %ld %lf  %lf\n", -1, xx[k], yy[k]);
					  }

					m = 1;

					for(j=0; j< 2*Nchannels+2; j++)
					{

						if(m == 1) add = wallwidth; else add = w_ch;

								
						for(k=0; k<500; k++)
						  {
							xx[k]= (2.*k /499.  - 1.) * depth/2.;
							yy[k] = - width/2. - sqrt(sq(radius_of_curv) - sq(depth/2.)) + sqrt(sq(radius_of_curv) - sq(xx[k])) + shift_y;

							fprintf(GeomFilePtr, "%ld   %lf  %lf\n", j, xx[k], yy[k]);		
						  }

							x_ch[0][j] = xx[0];
							y_ch[0][j] = yy[0];
							x_ch[1][j] = xx[167];
							y_ch[1][j] = yy[167];
							x_ch[2][j] = xx[333];
							y_ch[2][j] = yy[333];
							x_ch[3][j] = xx[499];
							y_ch[3][j] = yy[499];

							shift_y += add;
							m = m * (-1);
					}
				  /* activating shadowing cylinder */

					x_ch[0][0] = x_ch[0][1] = x_ch[0][2*Nchannels] = x_ch[0][2*Nchannels+1] = - main_depth/2.;
				
					x_ch[3][0] = x_ch[3][1] = x_ch[3][2*Nchannels] = x_ch[3][2*Nchannels+1] = main_depth/2.;
	

			  }
		}
    }
		
  if(Option==1) coef_pi=1.; else coef_pi=2.;  fprintf(LogFilePtr,"Phase set is %f°.\n", 180./M_PI*fmod(Phase , coef_pi*M_PI));   

}/* End OwnInit */


/* own cleanup of the monochromator/analyser module */

void ChopperFermiCleanup()
{
  if (GeomFilePtr) fclose(GeomFilePtr);
}/* End OwnCleanup */

#endif

 #undef   MCSTAS_TRACE
 vitess2mcstas(Neutrons, &x, &y, &z, &vx, &vy, &vz, &t, &sx, &sy, &sz, &p);
}
#line 20307 "./Test_Fermi.c"
}   /* End of FC_Vitess=Vitess_ChopperFermi() SETTING parameter declarations. */
#undef optimal_wl
#undef omega
#undef pos_ch
#undef Neutrons
#undef Dir
#undef Pos
#undef GeomFilePtr
#undef XFILEName
#undef coef_pi
#undef x_ch
#undef y_ch
#undef phase0
#undef angle_channel
#undef shift_y
#undef main_depth
#undef radius_of_curv
#undef WL
#undef TOF
#undef CurvGeomOption
#undef Option
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompFC_Vitess:
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

  /* TRACE Component Monitor2_xt [11] */
  mccoordschange(mcposrMonitor2_xt, mcrotrMonitor2_xt,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component Monitor2_xt (without coords transformations) */
  mcJumpTrace_Monitor2_xt:
  SIG_MESSAGE("Monitor2_xt (Trace)");
  mcDEBUG_COMP("Monitor2_xt")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompMonitor2_xt
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
#define mccompcurname  Monitor2_xt
#define mccompcurtype  Monitor_nD
#define mccompcurindex 11
#define user1 mccMonitor2_xt_user1
#define user2 mccMonitor2_xt_user2
#define user3 mccMonitor2_xt_user3
#define DEFS mccMonitor2_xt_DEFS
#define Vars mccMonitor2_xt_Vars
#define detector mccMonitor2_xt_detector
#define offdata mccMonitor2_xt_offdata
{   /* Declarations of Monitor2_xt=Monitor_nD() SETTING parameters. */
MCNUM xwidth = mccMonitor2_xt_xwidth;
MCNUM yheight = mccMonitor2_xt_yheight;
MCNUM zdepth = mccMonitor2_xt_zdepth;
MCNUM xmin = mccMonitor2_xt_xmin;
MCNUM xmax = mccMonitor2_xt_xmax;
MCNUM ymin = mccMonitor2_xt_ymin;
MCNUM ymax = mccMonitor2_xt_ymax;
MCNUM zmin = mccMonitor2_xt_zmin;
MCNUM zmax = mccMonitor2_xt_zmax;
MCNUM bins = mccMonitor2_xt_bins;
MCNUM min = mccMonitor2_xt_min;
MCNUM max = mccMonitor2_xt_max;
MCNUM restore_neutron = mccMonitor2_xt_restore_neutron;
MCNUM radius = mccMonitor2_xt_radius;
char* options = mccMonitor2_xt_options;
char* filename = mccMonitor2_xt_filename;
char* geometry = mccMonitor2_xt_geometry;
char* username1 = mccMonitor2_xt_username1;
char* username2 = mccMonitor2_xt_username2;
char* username3 = mccMonitor2_xt_username3;
int nowritefile = mccMonitor2_xt_nowritefile;
#line 312 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/Monitor_nD.comp"
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
#line 20630 "./Test_Fermi.c"
}   /* End of Monitor2_xt=Monitor_nD() SETTING parameter declarations. */
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
  mcabsorbCompMonitor2_xt:
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
#line 20746 "./Test_Fermi.c"
}   /* End of Origin=Progress_bar() SETTING parameter declarations. */
#undef CurrentTime
#undef EndTime
#undef StartTime
#undef IntermediateCnts
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'Monitor1_xt'. */
  SIG_MESSAGE("Monitor1_xt (Save)");
#define mccompcurname  Monitor1_xt
#define mccompcurtype  Monitor_nD
#define mccompcurindex 3
#define user1 mccMonitor1_xt_user1
#define user2 mccMonitor1_xt_user2
#define user3 mccMonitor1_xt_user3
#define DEFS mccMonitor1_xt_DEFS
#define Vars mccMonitor1_xt_Vars
#define detector mccMonitor1_xt_detector
#define offdata mccMonitor1_xt_offdata
{   /* Declarations of Monitor1_xt=Monitor_nD() SETTING parameters. */
MCNUM xwidth = mccMonitor1_xt_xwidth;
MCNUM yheight = mccMonitor1_xt_yheight;
MCNUM zdepth = mccMonitor1_xt_zdepth;
MCNUM xmin = mccMonitor1_xt_xmin;
MCNUM xmax = mccMonitor1_xt_xmax;
MCNUM ymin = mccMonitor1_xt_ymin;
MCNUM ymax = mccMonitor1_xt_ymax;
MCNUM zmin = mccMonitor1_xt_zmin;
MCNUM zmax = mccMonitor1_xt_zmax;
MCNUM bins = mccMonitor1_xt_bins;
MCNUM min = mccMonitor1_xt_min;
MCNUM max = mccMonitor1_xt_max;
MCNUM restore_neutron = mccMonitor1_xt_restore_neutron;
MCNUM radius = mccMonitor1_xt_radius;
char* options = mccMonitor1_xt_options;
char* filename = mccMonitor1_xt_filename;
char* geometry = mccMonitor1_xt_geometry;
char* username1 = mccMonitor1_xt_username1;
char* username2 = mccMonitor1_xt_username2;
char* username3 = mccMonitor1_xt_username3;
int nowritefile = mccMonitor1_xt_nowritefile;
#line 482 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/Monitor_nD.comp"
{
  /* save results, but do not free pointers */
  if (!nowritefile) {
    detector = Monitor_nD_Save(&DEFS, &Vars);
  }
}
#line 20797 "./Test_Fermi.c"
}   /* End of Monitor1_xt=Monitor_nD() SETTING parameter declarations. */
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

  /* User SAVE code for component 'FC_McStas'. */
  SIG_MESSAGE("FC_McStas (Save)");
#define mccompcurname  FC_McStas
#define mccompcurtype  FermiChopper
#define mccompcurindex 7
#define FCVars mccFC_McStas_FCVars
{   /* Declarations of FC_McStas=FermiChopper() SETTING parameters. */
MCNUM phase = mccFC_McStas_phase;
MCNUM radius = mccFC_McStas_radius;
MCNUM nu = mccFC_McStas_nu;
MCNUM w = mccFC_McStas_w;
MCNUM nslit = mccFC_McStas_nslit;
MCNUM R0 = mccFC_McStas_R0;
MCNUM Qc = mccFC_McStas_Qc;
MCNUM alpha = mccFC_McStas_alpha;
MCNUM m = mccFC_McStas_m;
MCNUM W = mccFC_McStas_W;
MCNUM length = mccFC_McStas_length;
MCNUM eff = mccFC_McStas_eff;
MCNUM zero_time = mccFC_McStas_zero_time;
MCNUM xwidth = mccFC_McStas_xwidth;
MCNUM verbose = mccFC_McStas_verbose;
MCNUM yheight = mccFC_McStas_yheight;
MCNUM curvature = mccFC_McStas_curvature;
MCNUM delay = mccFC_McStas_delay;
#line 785 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/FermiChopper.comp"
{
  double mean_k, mean_v, mean_t, mean_w=0, mean_L=0;

  if (FCVars.sum_N) {
    mean_v = FCVars.sum_v/FCVars.sum_N;
    mean_t = FCVars.sum_t/FCVars.sum_N;
    mean_k = V2K*mean_v;
    if (mean_k) mean_L = 2*PI/mean_k;
    mean_w = VS2E*mean_v*mean_v;
    /* base opening time */
    double div, mean_phase;
    div        = atan2(w,length)/PI*180;
    mean_phase = fmod(mean_t*nu*360, 360);
    mean_phase+=radius/mean_v*360*nu;
    if (mean_phase > 180) mean_phase -= 360;

    if (!FCVars.sum_N_pass)
      printf("FermiChopper: %s: No neutron can pass the chopper.\n", NAME_CURRENT_COMP);
    if (!FCVars.sum_N_pass || verbose)
      printf("FermiChopper: %s\n"
           "              Mean velocity     v     = %g [m/s]\n"
           "              Mean wavelength   lambda= %g [Angs]\n"
           "              Mean energy       omega = %g [meV]\n"
           "              Mean arrival time t     = %g [s]\n"
           "              Mean phase              = %g [deg] (%s)\n"
           "              Slit pack divergence    = %g [deg] (full width)\n"
           "              Opening time      dt    = %g [s]\n"
           "              Intensity reaching FC   = %g [n/s]\n"
           "              Intensity passing  FC   = %g [n/s]\n"
           , NAME_CURRENT_COMP,
           mean_v, mean_L, mean_w, mean_t, mean_phase,
           (zero_time > 1 ? "set automatically" : "use phase=-(this value) to optimize"),
           2*div,
           (nu ? fabs(div/PI/nu) : 1),
           FCVars.sum_N,
           FCVars.sum_N_pass);
    if (!FCVars.sum_N_pass || verbose) {
      printf("FermiChopper: %s: Lost events anaylsis\n"
             "              Already inside:            %li\n"
             "              By Top/Bottom of cylinder: %li\n"
             "              At cylinder entrance:      %li\n"
             "              Hit cyl. entrance sides:   %li\n"
             "              Can't prop. to slit pack:  %li\n"
             "              At slit pack entrance:     %li\n"
             "              On absorbing slit coating: %li\n"
             "              Exiting slit pack:         %li\n"
             "              Too many iterations:       %li\n"
             "              Prop. in wrong direction : %li\n"
             "              Mad neutron (no control):  %li\n"
             "              At cylinder exit:          %li\n"
      , NAME_CURRENT_COMP,
      FCVars.absorb_alreadyinside,
      FCVars.absorb_topbottom,
      FCVars.absorb_cylentrance,
      FCVars.absorb_sideentrance,
      FCVars.absorb_notreachentrance,
      FCVars.absorb_packentrance,
      FCVars.absorb_slitcoating,
      
      FCVars.absorb_exitslitpack,
      FCVars.absorb_maxiterations,
      FCVars.absorb_wrongdirection,
      FCVars.absorb_nocontrol,
      FCVars.absorb_cylexit);
      
      if (FCVars.warn_notreachslitwall || FCVars.warn_notreachslitoutput)
        printf("Warning:      Can not reach slit wall:   %li\n"
               "Warning:      Can not reach slit output: %li\n",
        FCVars.warn_notreachslitwall,
        FCVars.warn_notreachslitoutput);
    }

  } else {
    printf("FermiChopper: %s: No neutron can reach the chopper.\n", NAME_CURRENT_COMP);
  }
}
#line 20912 "./Test_Fermi.c"
}   /* End of FC_McStas=FermiChopper() SETTING parameter declarations. */
#undef FCVars
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'Monitor2_xt'. */
  SIG_MESSAGE("Monitor2_xt (Save)");
#define mccompcurname  Monitor2_xt
#define mccompcurtype  Monitor_nD
#define mccompcurindex 11
#define user1 mccMonitor2_xt_user1
#define user2 mccMonitor2_xt_user2
#define user3 mccMonitor2_xt_user3
#define DEFS mccMonitor2_xt_DEFS
#define Vars mccMonitor2_xt_Vars
#define detector mccMonitor2_xt_detector
#define offdata mccMonitor2_xt_offdata
{   /* Declarations of Monitor2_xt=Monitor_nD() SETTING parameters. */
MCNUM xwidth = mccMonitor2_xt_xwidth;
MCNUM yheight = mccMonitor2_xt_yheight;
MCNUM zdepth = mccMonitor2_xt_zdepth;
MCNUM xmin = mccMonitor2_xt_xmin;
MCNUM xmax = mccMonitor2_xt_xmax;
MCNUM ymin = mccMonitor2_xt_ymin;
MCNUM ymax = mccMonitor2_xt_ymax;
MCNUM zmin = mccMonitor2_xt_zmin;
MCNUM zmax = mccMonitor2_xt_zmax;
MCNUM bins = mccMonitor2_xt_bins;
MCNUM min = mccMonitor2_xt_min;
MCNUM max = mccMonitor2_xt_max;
MCNUM restore_neutron = mccMonitor2_xt_restore_neutron;
MCNUM radius = mccMonitor2_xt_radius;
char* options = mccMonitor2_xt_options;
char* filename = mccMonitor2_xt_filename;
char* geometry = mccMonitor2_xt_geometry;
char* username1 = mccMonitor2_xt_username1;
char* username2 = mccMonitor2_xt_username2;
char* username3 = mccMonitor2_xt_username3;
int nowritefile = mccMonitor2_xt_nowritefile;
#line 482 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/Monitor_nD.comp"
{
  /* save results, but do not free pointers */
  if (!nowritefile) {
    detector = Monitor_nD_Save(&DEFS, &Vars);
  }
}
#line 20960 "./Test_Fermi.c"
}   /* End of Monitor2_xt=Monitor_nD() SETTING parameter declarations. */
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
#line 21007 "./Test_Fermi.c"
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
  /* User FINALLY code for component 'Source'. */
  SIG_MESSAGE("Source (Finally)");
#define mccompcurname  Source
#define mccompcurtype  Source_gen
#define mccompcurindex 2
#define p_in mccSource_p_in
#define lambda1 mccSource_lambda1
#define lambda2 mccSource_lambda2
#define lambda3 mccSource_lambda3
#define pTable mccSource_pTable
#define pTable_x mccSource_pTable_x
#define pTable_y mccSource_pTable_y
#define pTable_xmin mccSource_pTable_xmin
#define pTable_xmax mccSource_pTable_xmax
#define pTable_xsum mccSource_pTable_xsum
#define pTable_ymin mccSource_pTable_ymin
#define pTable_ymax mccSource_pTable_ymax
#define pTable_ysum mccSource_pTable_ysum
#define pTable_dxmin mccSource_pTable_dxmin
#define pTable_dxmax mccSource_pTable_dxmax
#define pTable_dymin mccSource_pTable_dymin
#define pTable_dymax mccSource_pTable_dymax
{   /* Declarations of Source=Source_gen() SETTING parameters. */
char* flux_file = mccSource_flux_file;
char* xdiv_file = mccSource_xdiv_file;
char* ydiv_file = mccSource_ydiv_file;
MCNUM radius = mccSource_radius;
MCNUM dist = mccSource_dist;
MCNUM focus_xw = mccSource_focus_xw;
MCNUM focus_yh = mccSource_focus_yh;
MCNUM focus_aw = mccSource_focus_aw;
MCNUM focus_ah = mccSource_focus_ah;
MCNUM E0 = mccSource_E0;
MCNUM dE = mccSource_dE;
MCNUM lambda0 = mccSource_lambda0;
MCNUM dlambda = mccSource_dlambda;
MCNUM I1 = mccSource_I1;
MCNUM yheight = mccSource_yheight;
MCNUM xwidth = mccSource_xwidth;
MCNUM verbose = mccSource_verbose;
MCNUM T1 = mccSource_T1;
MCNUM flux_file_perAA = mccSource_flux_file_perAA;
MCNUM flux_file_log = mccSource_flux_file_log;
MCNUM Lmin = mccSource_Lmin;
MCNUM Lmax = mccSource_Lmax;
MCNUM Emin = mccSource_Emin;
MCNUM Emax = mccSource_Emax;
MCNUM T2 = mccSource_T2;
MCNUM I2 = mccSource_I2;
MCNUM T3 = mccSource_T3;
MCNUM I3 = mccSource_I3;
MCNUM zdepth = mccSource_zdepth;
int target_index = mccSource_target_index;
#line 571 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../sources/Source_gen.comp"
{
  Table_Free(&pTable);
  Table_Free(&pTable_x);
  Table_Free(&pTable_y);
}
#line 21078 "./Test_Fermi.c"
}   /* End of Source=Source_gen() SETTING parameter declarations. */
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

    if (!mcNCounter[2]) fprintf(stderr, "Warning: No neutron could reach Component[2] Source\n");
    if (mcAbsorbProp[2]) fprintf(stderr, "Warning: %g events were removed in Component[2] Source=Source_gen()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[2]);
  /* User FINALLY code for component 'Monitor1_xt'. */
  SIG_MESSAGE("Monitor1_xt (Finally)");
#define mccompcurname  Monitor1_xt
#define mccompcurtype  Monitor_nD
#define mccompcurindex 3
#define user1 mccMonitor1_xt_user1
#define user2 mccMonitor1_xt_user2
#define user3 mccMonitor1_xt_user3
#define DEFS mccMonitor1_xt_DEFS
#define Vars mccMonitor1_xt_Vars
#define detector mccMonitor1_xt_detector
#define offdata mccMonitor1_xt_offdata
{   /* Declarations of Monitor1_xt=Monitor_nD() SETTING parameters. */
MCNUM xwidth = mccMonitor1_xt_xwidth;
MCNUM yheight = mccMonitor1_xt_yheight;
MCNUM zdepth = mccMonitor1_xt_zdepth;
MCNUM xmin = mccMonitor1_xt_xmin;
MCNUM xmax = mccMonitor1_xt_xmax;
MCNUM ymin = mccMonitor1_xt_ymin;
MCNUM ymax = mccMonitor1_xt_ymax;
MCNUM zmin = mccMonitor1_xt_zmin;
MCNUM zmax = mccMonitor1_xt_zmax;
MCNUM bins = mccMonitor1_xt_bins;
MCNUM min = mccMonitor1_xt_min;
MCNUM max = mccMonitor1_xt_max;
MCNUM restore_neutron = mccMonitor1_xt_restore_neutron;
MCNUM radius = mccMonitor1_xt_radius;
char* options = mccMonitor1_xt_options;
char* filename = mccMonitor1_xt_filename;
char* geometry = mccMonitor1_xt_geometry;
char* username1 = mccMonitor1_xt_username1;
char* username2 = mccMonitor1_xt_username2;
char* username3 = mccMonitor1_xt_username3;
int nowritefile = mccMonitor1_xt_nowritefile;
#line 490 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/Monitor_nD.comp"
{
  /* free pointers */
  Monitor_nD_Finally(&DEFS, &Vars);
}
#line 21142 "./Test_Fermi.c"
}   /* End of Monitor1_xt=Monitor_nD() SETTING parameter declarations. */
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

    if (!mcNCounter[3]) fprintf(stderr, "Warning: No neutron could reach Component[3] Monitor1_xt\n");
    if (mcAbsorbProp[3]) fprintf(stderr, "Warning: %g events were removed in Component[3] Monitor1_xt=Monitor_nD()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[3]);
    if (!mcNCounter[4]) fprintf(stderr, "Warning: No neutron could reach Component[4] FC_Position\n");
    if (mcAbsorbProp[4]) fprintf(stderr, "Warning: %g events were removed in Component[4] FC_Position=Arm()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[4]);
  /* User FINALLY code for component 'FC_GuideG'. */
  SIG_MESSAGE("FC_GuideG (Finally)");
#define mccompcurname  FC_GuideG
#define mccompcurtype  Guide_gravity
#define mccompcurindex 5
#define GVars mccFC_GuideG_GVars
#define pTable mccFC_GuideG_pTable
{   /* Declarations of FC_GuideG=Guide_gravity() SETTING parameters. */
MCNUM w1 = mccFC_GuideG_w1;
MCNUM h1 = mccFC_GuideG_h1;
MCNUM w2 = mccFC_GuideG_w2;
MCNUM h2 = mccFC_GuideG_h2;
MCNUM l = mccFC_GuideG_l;
MCNUM R0 = mccFC_GuideG_R0;
MCNUM Qc = mccFC_GuideG_Qc;
MCNUM alpha = mccFC_GuideG_alpha;
MCNUM m = mccFC_GuideG_m;
MCNUM W = mccFC_GuideG_W;
MCNUM nslit = mccFC_GuideG_nslit;
MCNUM d = mccFC_GuideG_d;
MCNUM mleft = mccFC_GuideG_mleft;
MCNUM mright = mccFC_GuideG_mright;
MCNUM mtop = mccFC_GuideG_mtop;
MCNUM mbottom = mccFC_GuideG_mbottom;
MCNUM nhslit = mccFC_GuideG_nhslit;
MCNUM G = mccFC_GuideG_G;
MCNUM aleft = mccFC_GuideG_aleft;
MCNUM aright = mccFC_GuideG_aright;
MCNUM atop = mccFC_GuideG_atop;
MCNUM abottom = mccFC_GuideG_abottom;
MCNUM wavy = mccFC_GuideG_wavy;
MCNUM wavy_z = mccFC_GuideG_wavy_z;
MCNUM wavy_tb = mccFC_GuideG_wavy_tb;
MCNUM wavy_lr = mccFC_GuideG_wavy_lr;
MCNUM chamfers = mccFC_GuideG_chamfers;
MCNUM chamfers_z = mccFC_GuideG_chamfers_z;
MCNUM chamfers_lr = mccFC_GuideG_chamfers_lr;
MCNUM chamfers_tb = mccFC_GuideG_chamfers_tb;
MCNUM nelements = mccFC_GuideG_nelements;
MCNUM nu = mccFC_GuideG_nu;
MCNUM phase = mccFC_GuideG_phase;
char* reflect = mccFC_GuideG_reflect;
#line 562 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Guide_gravity.comp"
{
if (GVars.warnings > 100) {
  fprintf(stderr,"%s: warning: neutron has entered guide, but can not exit !\n", GVars.compcurname);
  fprintf(stderr,"%s: warning: This message has been repeated %g times\n", GVars.compcurname, GVars.warnings);
}
}
#line 21208 "./Test_Fermi.c"
}   /* End of FC_GuideG=Guide_gravity() SETTING parameter declarations. */
#undef pTable
#undef GVars
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[5]) fprintf(stderr, "Warning: No neutron could reach Component[5] FC_GuideG\n");
    if (mcAbsorbProp[5]) fprintf(stderr, "Warning: %g events were removed in Component[5] FC_GuideG=Guide_gravity()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[5]);
    if (!mcNCounter[6]) fprintf(stderr, "Warning: No neutron could reach Component[6] FC_GuideC\n");
    if (mcAbsorbProp[6]) fprintf(stderr, "Warning: %g events were removed in Component[6] FC_GuideC=Guide_channeled()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[6]);
    if (!mcNCounter[7]) fprintf(stderr, "Warning: No neutron could reach Component[7] FC_McStas\n");
    if (mcAbsorbProp[7]) fprintf(stderr, "Warning: %g events were removed in Component[7] FC_McStas=FermiChopper()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[7]);
    if (!mcNCounter[8]) fprintf(stderr, "Warning: No neutron could reach Component[8] FC_ILL\n");
    if (mcAbsorbProp[8]) fprintf(stderr, "Warning: %g events were removed in Component[8] FC_ILL=FermiChopper_ILL()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[8]);
    if (!mcNCounter[9]) fprintf(stderr, "Warning: No neutron could reach Component[9] Fake_Origin\n");
    if (mcAbsorbProp[9]) fprintf(stderr, "Warning: %g events were removed in Component[9] Fake_Origin=Arm()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[9]);
  /* User FINALLY code for component 'FC_Vitess'. */
  SIG_MESSAGE("FC_Vitess (Finally)");
#define mccompcurname  FC_Vitess
#define mccompcurtype  Vitess_ChopperFermi
#define mccompcurindex 10
#define Option mccFC_Vitess_Option
#define CurvGeomOption mccFC_Vitess_CurvGeomOption
#define TOF mccFC_Vitess_TOF
#define WL mccFC_Vitess_WL
#define radius_of_curv mccFC_Vitess_radius_of_curv
#define main_depth mccFC_Vitess_main_depth
#define shift_y mccFC_Vitess_shift_y
#define angle_channel mccFC_Vitess_angle_channel
#define phase0 mccFC_Vitess_phase0
#define y_ch mccFC_Vitess_y_ch
#define x_ch mccFC_Vitess_x_ch
#define coef_pi mccFC_Vitess_coef_pi
#define XFILEName mccFC_Vitess_XFILEName
#define GeomFilePtr mccFC_Vitess_GeomFilePtr
#define Pos mccFC_Vitess_Pos
#define Dir mccFC_Vitess_Dir
#define Neutrons mccFC_Vitess_Neutrons
#define pos_ch mccFC_Vitess_pos_ch
#define omega mccFC_Vitess_omega
#define optimal_wl mccFC_Vitess_optimal_wl
{   /* Declarations of FC_Vitess=Vitess_ChopperFermi() SETTING parameters. */
char* sGeomFileName = mccFC_Vitess_sGeomFileName;
int GeomOption = mccFC_Vitess_GeomOption;
int zerotime = mccFC_Vitess_zerotime;
int Nchannels = mccFC_Vitess_Nchannels;
int Ngates = mccFC_Vitess_Ngates;
MCNUM freq = mccFC_Vitess_freq;
MCNUM height = mccFC_Vitess_height;
MCNUM width = mccFC_Vitess_width;
MCNUM depth = mccFC_Vitess_depth;
MCNUM r_curv = mccFC_Vitess_r_curv;
MCNUM diameter = mccFC_Vitess_diameter;
MCNUM Phase = mccFC_Vitess_Phase;
MCNUM wallwidth = mccFC_Vitess_wallwidth;
#line 199 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Vitess_ChopperFermi.comp"
{
 ChopperFermiCleanup();
 McCleanupVt();
}
#line 21270 "./Test_Fermi.c"
}   /* End of FC_Vitess=Vitess_ChopperFermi() SETTING parameter declarations. */
#undef optimal_wl
#undef omega
#undef pos_ch
#undef Neutrons
#undef Dir
#undef Pos
#undef GeomFilePtr
#undef XFILEName
#undef coef_pi
#undef x_ch
#undef y_ch
#undef phase0
#undef angle_channel
#undef shift_y
#undef main_depth
#undef radius_of_curv
#undef WL
#undef TOF
#undef CurvGeomOption
#undef Option
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[10]) fprintf(stderr, "Warning: No neutron could reach Component[10] FC_Vitess\n");
    if (mcAbsorbProp[10]) fprintf(stderr, "Warning: %g events were removed in Component[10] FC_Vitess=Vitess_ChopperFermi()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[10]);
  /* User FINALLY code for component 'Monitor2_xt'. */
  SIG_MESSAGE("Monitor2_xt (Finally)");
#define mccompcurname  Monitor2_xt
#define mccompcurtype  Monitor_nD
#define mccompcurindex 11
#define user1 mccMonitor2_xt_user1
#define user2 mccMonitor2_xt_user2
#define user3 mccMonitor2_xt_user3
#define DEFS mccMonitor2_xt_DEFS
#define Vars mccMonitor2_xt_Vars
#define detector mccMonitor2_xt_detector
#define offdata mccMonitor2_xt_offdata
{   /* Declarations of Monitor2_xt=Monitor_nD() SETTING parameters. */
MCNUM xwidth = mccMonitor2_xt_xwidth;
MCNUM yheight = mccMonitor2_xt_yheight;
MCNUM zdepth = mccMonitor2_xt_zdepth;
MCNUM xmin = mccMonitor2_xt_xmin;
MCNUM xmax = mccMonitor2_xt_xmax;
MCNUM ymin = mccMonitor2_xt_ymin;
MCNUM ymax = mccMonitor2_xt_ymax;
MCNUM zmin = mccMonitor2_xt_zmin;
MCNUM zmax = mccMonitor2_xt_zmax;
MCNUM bins = mccMonitor2_xt_bins;
MCNUM min = mccMonitor2_xt_min;
MCNUM max = mccMonitor2_xt_max;
MCNUM restore_neutron = mccMonitor2_xt_restore_neutron;
MCNUM radius = mccMonitor2_xt_radius;
char* options = mccMonitor2_xt_options;
char* filename = mccMonitor2_xt_filename;
char* geometry = mccMonitor2_xt_geometry;
char* username1 = mccMonitor2_xt_username1;
char* username2 = mccMonitor2_xt_username2;
char* username3 = mccMonitor2_xt_username3;
int nowritefile = mccMonitor2_xt_nowritefile;
#line 490 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/Monitor_nD.comp"
{
  /* free pointers */
  Monitor_nD_Finally(&DEFS, &Vars);
}
#line 21337 "./Test_Fermi.c"
}   /* End of Monitor2_xt=Monitor_nD() SETTING parameter declarations. */
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

    if (!mcNCounter[11]) fprintf(stderr, "Warning: No neutron could reach Component[11] Monitor2_xt\n");
    if (mcAbsorbProp[11]) fprintf(stderr, "Warning: %g events were removed in Component[11] Monitor2_xt=Monitor_nD()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[11]);
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
#line 21386 "./Test_Fermi.c"
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
#define mccompcurtype  Source_gen
#define mccompcurindex 2
#define p_in mccSource_p_in
#define lambda1 mccSource_lambda1
#define lambda2 mccSource_lambda2
#define lambda3 mccSource_lambda3
#define pTable mccSource_pTable
#define pTable_x mccSource_pTable_x
#define pTable_y mccSource_pTable_y
#define pTable_xmin mccSource_pTable_xmin
#define pTable_xmax mccSource_pTable_xmax
#define pTable_xsum mccSource_pTable_xsum
#define pTable_ymin mccSource_pTable_ymin
#define pTable_ymax mccSource_pTable_ymax
#define pTable_ysum mccSource_pTable_ysum
#define pTable_dxmin mccSource_pTable_dxmin
#define pTable_dxmax mccSource_pTable_dxmax
#define pTable_dymin mccSource_pTable_dymin
#define pTable_dymax mccSource_pTable_dymax
{   /* Declarations of Source=Source_gen() SETTING parameters. */
char* flux_file = mccSource_flux_file;
char* xdiv_file = mccSource_xdiv_file;
char* ydiv_file = mccSource_ydiv_file;
MCNUM radius = mccSource_radius;
MCNUM dist = mccSource_dist;
MCNUM focus_xw = mccSource_focus_xw;
MCNUM focus_yh = mccSource_focus_yh;
MCNUM focus_aw = mccSource_focus_aw;
MCNUM focus_ah = mccSource_focus_ah;
MCNUM E0 = mccSource_E0;
MCNUM dE = mccSource_dE;
MCNUM lambda0 = mccSource_lambda0;
MCNUM dlambda = mccSource_dlambda;
MCNUM I1 = mccSource_I1;
MCNUM yheight = mccSource_yheight;
MCNUM xwidth = mccSource_xwidth;
MCNUM verbose = mccSource_verbose;
MCNUM T1 = mccSource_T1;
MCNUM flux_file_perAA = mccSource_flux_file_perAA;
MCNUM flux_file_log = mccSource_flux_file_log;
MCNUM Lmin = mccSource_Lmin;
MCNUM Lmax = mccSource_Lmax;
MCNUM Emin = mccSource_Emin;
MCNUM Emax = mccSource_Emax;
MCNUM T2 = mccSource_T2;
MCNUM I2 = mccSource_I2;
MCNUM T3 = mccSource_T3;
MCNUM I3 = mccSource_I3;
MCNUM zdepth = mccSource_zdepth;
int target_index = mccSource_target_index;
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
#line 21499 "./Test_Fermi.c"
}   /* End of Source=Source_gen() SETTING parameter declarations. */
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

  /* MCDISPLAY code for component 'Monitor1_xt'. */
  SIG_MESSAGE("Monitor1_xt (McDisplay)");
  printf("MCDISPLAY: component %s\n", "Monitor1_xt");
#define mccompcurname  Monitor1_xt
#define mccompcurtype  Monitor_nD
#define mccompcurindex 3
#define user1 mccMonitor1_xt_user1
#define user2 mccMonitor1_xt_user2
#define user3 mccMonitor1_xt_user3
#define DEFS mccMonitor1_xt_DEFS
#define Vars mccMonitor1_xt_Vars
#define detector mccMonitor1_xt_detector
#define offdata mccMonitor1_xt_offdata
{   /* Declarations of Monitor1_xt=Monitor_nD() SETTING parameters. */
MCNUM xwidth = mccMonitor1_xt_xwidth;
MCNUM yheight = mccMonitor1_xt_yheight;
MCNUM zdepth = mccMonitor1_xt_zdepth;
MCNUM xmin = mccMonitor1_xt_xmin;
MCNUM xmax = mccMonitor1_xt_xmax;
MCNUM ymin = mccMonitor1_xt_ymin;
MCNUM ymax = mccMonitor1_xt_ymax;
MCNUM zmin = mccMonitor1_xt_zmin;
MCNUM zmax = mccMonitor1_xt_zmax;
MCNUM bins = mccMonitor1_xt_bins;
MCNUM min = mccMonitor1_xt_min;
MCNUM max = mccMonitor1_xt_max;
MCNUM restore_neutron = mccMonitor1_xt_restore_neutron;
MCNUM radius = mccMonitor1_xt_radius;
char* options = mccMonitor1_xt_options;
char* filename = mccMonitor1_xt_filename;
char* geometry = mccMonitor1_xt_geometry;
char* username1 = mccMonitor1_xt_username1;
char* username2 = mccMonitor1_xt_username2;
char* username3 = mccMonitor1_xt_username3;
int nowritefile = mccMonitor1_xt_nowritefile;
#line 496 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/Monitor_nD.comp"
{
  if (geometry && strlen(geometry) && strcmp(geometry,"0") && strcmp(geometry, "NULL"))
  {
    off_display(offdata);
  } else {
    Monitor_nD_McDisplay(&DEFS, &Vars);
  }
}
#line 21566 "./Test_Fermi.c"
}   /* End of Monitor1_xt=Monitor_nD() SETTING parameter declarations. */
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

  /* MCDISPLAY code for component 'FC_Position'. */
  SIG_MESSAGE("FC_Position (McDisplay)");
  printf("MCDISPLAY: component %s\n", "FC_Position");
#define mccompcurname  FC_Position
#define mccompcurtype  Arm
#define mccompcurindex 4
#line 40 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Arm.comp"
{
  /* A bit ugly; hard-coded dimensions. */
  
  line(0,0,0,0.2,0,0);
  line(0,0,0,0,0.2,0);
  line(0,0,0,0,0,0.2);
}
#line 21593 "./Test_Fermi.c"
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'FC_GuideG'. */
  SIG_MESSAGE("FC_GuideG (McDisplay)");
  printf("MCDISPLAY: component %s\n", "FC_GuideG");
#define mccompcurname  FC_GuideG
#define mccompcurtype  Guide_gravity
#define mccompcurindex 5
#define GVars mccFC_GuideG_GVars
#define pTable mccFC_GuideG_pTable
{   /* Declarations of FC_GuideG=Guide_gravity() SETTING parameters. */
MCNUM w1 = mccFC_GuideG_w1;
MCNUM h1 = mccFC_GuideG_h1;
MCNUM w2 = mccFC_GuideG_w2;
MCNUM h2 = mccFC_GuideG_h2;
MCNUM l = mccFC_GuideG_l;
MCNUM R0 = mccFC_GuideG_R0;
MCNUM Qc = mccFC_GuideG_Qc;
MCNUM alpha = mccFC_GuideG_alpha;
MCNUM m = mccFC_GuideG_m;
MCNUM W = mccFC_GuideG_W;
MCNUM nslit = mccFC_GuideG_nslit;
MCNUM d = mccFC_GuideG_d;
MCNUM mleft = mccFC_GuideG_mleft;
MCNUM mright = mccFC_GuideG_mright;
MCNUM mtop = mccFC_GuideG_mtop;
MCNUM mbottom = mccFC_GuideG_mbottom;
MCNUM nhslit = mccFC_GuideG_nhslit;
MCNUM G = mccFC_GuideG_G;
MCNUM aleft = mccFC_GuideG_aleft;
MCNUM aright = mccFC_GuideG_aright;
MCNUM atop = mccFC_GuideG_atop;
MCNUM abottom = mccFC_GuideG_abottom;
MCNUM wavy = mccFC_GuideG_wavy;
MCNUM wavy_z = mccFC_GuideG_wavy_z;
MCNUM wavy_tb = mccFC_GuideG_wavy_tb;
MCNUM wavy_lr = mccFC_GuideG_wavy_lr;
MCNUM chamfers = mccFC_GuideG_chamfers;
MCNUM chamfers_z = mccFC_GuideG_chamfers_z;
MCNUM chamfers_lr = mccFC_GuideG_chamfers_lr;
MCNUM chamfers_tb = mccFC_GuideG_chamfers_tb;
MCNUM nelements = mccFC_GuideG_nelements;
MCNUM nu = mccFC_GuideG_nu;
MCNUM phase = mccFC_GuideG_phase;
char* reflect = mccFC_GuideG_reflect;
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
#line 21704 "./Test_Fermi.c"
}   /* End of FC_GuideG=Guide_gravity() SETTING parameter declarations. */
#undef pTable
#undef GVars
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'FC_GuideC'. */
  SIG_MESSAGE("FC_GuideC (McDisplay)");
  printf("MCDISPLAY: component %s\n", "FC_GuideC");
#define mccompcurname  FC_GuideC
#define mccompcurtype  Guide_channeled
#define mccompcurindex 6
#define w1c mccFC_GuideC_w1c
#define w2c mccFC_GuideC_w2c
#define ww mccFC_GuideC_ww
#define hh mccFC_GuideC_hh
#define whalf mccFC_GuideC_whalf
#define hhalf mccFC_GuideC_hhalf
#define lwhalf mccFC_GuideC_lwhalf
#define lhhalf mccFC_GuideC_lhhalf
{   /* Declarations of FC_GuideC=Guide_channeled() SETTING parameters. */
MCNUM w1 = mccFC_GuideC_w1;
MCNUM h1 = mccFC_GuideC_h1;
MCNUM w2 = mccFC_GuideC_w2;
MCNUM h2 = mccFC_GuideC_h2;
MCNUM l = mccFC_GuideC_l;
MCNUM R0 = mccFC_GuideC_R0;
MCNUM Qc = mccFC_GuideC_Qc;
MCNUM alpha = mccFC_GuideC_alpha;
MCNUM m = mccFC_GuideC_m;
MCNUM nslit = mccFC_GuideC_nslit;
MCNUM d = mccFC_GuideC_d;
MCNUM Qcx = mccFC_GuideC_Qcx;
MCNUM Qcy = mccFC_GuideC_Qcy;
MCNUM alphax = mccFC_GuideC_alphax;
MCNUM alphay = mccFC_GuideC_alphay;
MCNUM W = mccFC_GuideC_W;
MCNUM mx = mccFC_GuideC_mx;
MCNUM my = mccFC_GuideC_my;
MCNUM nu = mccFC_GuideC_nu;
MCNUM phase = mccFC_GuideC_phase;
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
#line 21778 "./Test_Fermi.c"
}   /* End of FC_GuideC=Guide_channeled() SETTING parameter declarations. */
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

  /* MCDISPLAY code for component 'FC_McStas'. */
  SIG_MESSAGE("FC_McStas (McDisplay)");
  printf("MCDISPLAY: component %s\n", "FC_McStas");
#define mccompcurname  FC_McStas
#define mccompcurtype  FermiChopper
#define mccompcurindex 7
#define FCVars mccFC_McStas_FCVars
{   /* Declarations of FC_McStas=FermiChopper() SETTING parameters. */
MCNUM phase = mccFC_McStas_phase;
MCNUM radius = mccFC_McStas_radius;
MCNUM nu = mccFC_McStas_nu;
MCNUM w = mccFC_McStas_w;
MCNUM nslit = mccFC_McStas_nslit;
MCNUM R0 = mccFC_McStas_R0;
MCNUM Qc = mccFC_McStas_Qc;
MCNUM alpha = mccFC_McStas_alpha;
MCNUM m = mccFC_McStas_m;
MCNUM W = mccFC_McStas_W;
MCNUM length = mccFC_McStas_length;
MCNUM eff = mccFC_McStas_eff;
MCNUM zero_time = mccFC_McStas_zero_time;
MCNUM xwidth = mccFC_McStas_xwidth;
MCNUM verbose = mccFC_McStas_verbose;
MCNUM yheight = mccFC_McStas_yheight;
MCNUM curvature = mccFC_McStas_curvature;
MCNUM delay = mccFC_McStas_delay;
#line 863 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/FermiChopper.comp"
{
  double index_x=0;
  double index_z=0;
  double xpos, zpos;
  double Nz,Nx;
  double ymin=-yheight/2.0;
  double ymax= yheight/2.0;

  double omega=FCVars.omega;
  FCVars.omega=0;
  // FCVars.ph0  =0;
  Nz = (FCVars.C_slit ?  4 : 1);
  Nx = (nslit > 11    ? 11 : nslit);
  FCVars.C_slit *= -1;
  magnify("xz");
  /* cylinder top/center/bottom  */
  circle("xz", 0,ymax,0,radius);
  circle("xz", 0,0   ,0,radius);
  circle("xz", 0,ymin,0,radius);
  /* vertical lines to make a kind of volume */
  line( 0  ,ymin,-radius, 0  ,ymax,-radius);
  line( 0  ,ymin, radius, 0  ,ymax, radius);
  line(-radius,ymin, 0  ,-radius,ymax, 0  );
  line( radius,ymin, 0  , radius,ymax, 0  );
  /* slit package */
  for (index_x = -Nx/2; index_x < Nx/2; index_x++) {
    for (index_z = -Nz/2; index_z < Nz/2; index_z++) {
      double xs1, zs1, zs2;
      double xp1, xp2, zp1, zp2;
      zs1 = length*index_z/Nz;
      zs2 = length*(index_z+1)/Nz;
      xs1 = w*nslit*index_x/Nx;
      xp1 = FC_xrot(xs1, zs1, 0, FCVars);
      xp2 = FC_xrot(xs1, zs2, 0, FCVars);
      zp1 = FC_zrot(xs1, zs1, 0, FCVars);
      zp2 = FC_zrot(xs1, zs2, 0, FCVars);
      multiline(5, xp1, ymin, zp1,
                   xp1, ymax, zp1,
                   xp2, ymax, zp2,
                   xp2, ymin, zp2,
                   xp1, ymin, zp1);
    }
  }
  /* cylinder inner sides containing slit package */
  double xp1, xp2, zp1, zp2;
  xpos = nslit*w/2;
  zpos = sqrt(radius*radius - xpos*xpos);
  xp1 = FC_xrot(xpos, -zpos, 0, FCVars);
  xp2 = FC_xrot(xpos, +zpos, 0, FCVars);
  zp1 = FC_zrot(xpos, -zpos, 0, FCVars);
  zp2 = FC_zrot(xpos, +zpos, 0, FCVars);
  multiline(5,  xp1, ymin, zp1,
                xp1, ymax, zp1,
                xp2, ymax, zp2,
                xp2, ymin, zp2,
                xp1, ymin, zp1);
  xpos *= -1;
  xp1 = FC_xrot(xpos, -zpos, 0, FCVars);
  xp2 = FC_xrot(xpos, +zpos, 0, FCVars);
  zp1 = FC_zrot(xpos, -zpos, 0, FCVars);
  zp2 = FC_zrot(xpos, +zpos, 0, FCVars);
  multiline(5,  xp1, ymin, zp1,
                xp1, ymax, zp1,
                xp2, ymax, zp2,
                xp2, ymin, zp2,
                xp1, ymin, zp1);
  FCVars.omega=omega;
}
#line 21887 "./Test_Fermi.c"
}   /* End of FC_McStas=FermiChopper() SETTING parameter declarations. */
#undef FCVars
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'FC_ILL'. */
  SIG_MESSAGE("FC_ILL (McDisplay)");
  printf("MCDISPLAY: component %s\n", "FC_ILL");
#define mccompcurname  FC_ILL
#define mccompcurtype  FermiChopper_ILL
#define mccompcurindex 8
#define FCVars mccFC_ILL_FCVars
{   /* Declarations of FC_ILL=FermiChopper_ILL() SETTING parameters. */
MCNUM phase = mccFC_ILL_phase;
MCNUM radius = mccFC_ILL_radius;
MCNUM nu = mccFC_ILL_nu;
MCNUM yheight = mccFC_ILL_yheight;
MCNUM w = mccFC_ILL_w;
MCNUM nslit = mccFC_ILL_nslit;
MCNUM R0 = mccFC_ILL_R0;
MCNUM Qc = mccFC_ILL_Qc;
MCNUM alpha = mccFC_ILL_alpha;
MCNUM m = mccFC_ILL_m;
MCNUM W = mccFC_ILL_W;
MCNUM length = mccFC_ILL_length;
MCNUM eff = mccFC_ILL_eff;
MCNUM zero_time = mccFC_ILL_zero_time;
MCNUM xwidth = mccFC_ILL_xwidth;
MCNUM verbose = mccFC_ILL_verbose;
#line 643 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../contrib/FermiChopper_ILL.comp"
{
  double index=0;
  double xpos, zpos;
  double ymax = yheight/2; 
  double ymin = -ymax;
  
  /* cylinder top/center/bottom  */
  circle("xz", 0,ymax,0,radius);
  circle("xz", 0,0   ,0,radius);
  circle("xz", 0,ymin,0,radius);
  /* vertical lines to make a kind of volume */
  line( 0  ,ymin,-radius, 0  ,ymax,-radius);
  line( 0  ,ymin, radius, 0  ,ymax, radius);
  line(-radius,ymin, 0  ,-radius,ymax, 0  );
  line( radius,ymin, 0  , radius,ymax, 0  );
  /* slit package */
  index = -nslit/2;
  zpos  = length/2;
  for (index = -nslit/2; index < nslit/2; index++) {
    xpos = index*w;
    multiline(5, xpos, ymin, -zpos,
                 xpos, ymax, -zpos,
                 xpos, ymax, +zpos,
                 xpos, ymin, +zpos,
                 xpos, ymin, -zpos);
  }
  /* cylinder inner sides containing slit package */
  xpos = nslit*w/2;
  zpos = sqrt(radius*radius - xpos*xpos);
  multiline(5,   xpos, ymin, -zpos,
                 xpos, ymax, -zpos,
                 xpos, ymax, +zpos,
                 xpos, ymin, +zpos,
                 xpos, ymin, -zpos);
  xpos *= -1;
  multiline(5,   xpos, ymin, -zpos,
                 xpos, ymax, -zpos,
                 xpos, ymax, +zpos,
                 xpos, ymin, +zpos,
                 xpos, ymin, -zpos);
}
#line 21960 "./Test_Fermi.c"
}   /* End of FC_ILL=FermiChopper_ILL() SETTING parameter declarations. */
#undef FCVars
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'Fake_Origin'. */
  SIG_MESSAGE("Fake_Origin (McDisplay)");
  printf("MCDISPLAY: component %s\n", "Fake_Origin");
#define mccompcurname  Fake_Origin
#define mccompcurtype  Arm
#define mccompcurindex 9
#line 40 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Arm.comp"
{
  /* A bit ugly; hard-coded dimensions. */
  
  line(0,0,0,0.2,0,0);
  line(0,0,0,0,0.2,0);
  line(0,0,0,0,0,0.2);
}
#line 21981 "./Test_Fermi.c"
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'FC_Vitess'. */
  SIG_MESSAGE("FC_Vitess (McDisplay)");
  printf("MCDISPLAY: component %s\n", "FC_Vitess");
#define mccompcurname  FC_Vitess
#define mccompcurtype  Vitess_ChopperFermi
#define mccompcurindex 10
#define Option mccFC_Vitess_Option
#define CurvGeomOption mccFC_Vitess_CurvGeomOption
#define TOF mccFC_Vitess_TOF
#define WL mccFC_Vitess_WL
#define radius_of_curv mccFC_Vitess_radius_of_curv
#define main_depth mccFC_Vitess_main_depth
#define shift_y mccFC_Vitess_shift_y
#define angle_channel mccFC_Vitess_angle_channel
#define phase0 mccFC_Vitess_phase0
#define y_ch mccFC_Vitess_y_ch
#define x_ch mccFC_Vitess_x_ch
#define coef_pi mccFC_Vitess_coef_pi
#define XFILEName mccFC_Vitess_XFILEName
#define GeomFilePtr mccFC_Vitess_GeomFilePtr
#define Pos mccFC_Vitess_Pos
#define Dir mccFC_Vitess_Dir
#define Neutrons mccFC_Vitess_Neutrons
#define pos_ch mccFC_Vitess_pos_ch
#define omega mccFC_Vitess_omega
#define optimal_wl mccFC_Vitess_optimal_wl
{   /* Declarations of FC_Vitess=Vitess_ChopperFermi() SETTING parameters. */
char* sGeomFileName = mccFC_Vitess_sGeomFileName;
int GeomOption = mccFC_Vitess_GeomOption;
int zerotime = mccFC_Vitess_zerotime;
int Nchannels = mccFC_Vitess_Nchannels;
int Ngates = mccFC_Vitess_Ngates;
MCNUM freq = mccFC_Vitess_freq;
MCNUM height = mccFC_Vitess_height;
MCNUM width = mccFC_Vitess_width;
MCNUM depth = mccFC_Vitess_depth;
MCNUM r_curv = mccFC_Vitess_r_curv;
MCNUM diameter = mccFC_Vitess_diameter;
MCNUM Phase = mccFC_Vitess_Phase;
MCNUM wallwidth = mccFC_Vitess_wallwidth;
#line 205 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../optics/Vitess_ChopperFermi.comp"
{
  double index=0;
  double xpos, zpos, ymin, ymax, rad, w_ch;
  Nchannels = (Nchannels > 11 ? 11 : Nchannels);
  w_ch = (width - (Nchannels+1) * wallwidth) / Nchannels/100;
  rad  = diameter/2.0/100;
  ymin = -height/2.0/100;
  ymax =  height/2.0/100;
  
  /* cylinder top/center/bottom  */
  circle("xz", 0,ymax,0,rad);
  circle("xz", 0,0   ,0,rad);
  circle("xz", 0,ymin,0,rad);
  /* vertical lines to make a kind of volume */
  line( 0  ,ymin,-rad, 0  ,ymax,-rad);
  line( 0  ,ymin, rad, 0  ,ymax, rad);
  line(-rad,ymin, 0  ,-rad,ymax, 0  );
  line( rad,ymin, 0  , rad,ymax, 0  );
  /* slit package */
  index = -Nchannels/2;
  zpos  = depth/2/100;
  for (index = -Nchannels/2; index < Nchannels/2; index++) {
    xpos = index*w_ch;
    multiline(5, xpos, ymin, -zpos,
                 xpos, ymax, -zpos,
                 xpos, ymax, +zpos,
                 xpos, ymin, +zpos,
                 xpos, ymin, -zpos);
  }
  /* cylinder inner sides containing slit package */
  xpos = Nchannels*w_ch/2;
  zpos = sqrt(rad*rad - xpos*xpos);
  multiline(5,   xpos, ymin, -zpos,
                 xpos, ymax, -zpos,
                 xpos, ymax, +zpos,
                 xpos, ymin, +zpos,
                 xpos, ymin, -zpos);
  xpos *= -1;
  multiline(5,   xpos, ymin, -zpos,
                 xpos, ymax, -zpos,
                 xpos, ymax, +zpos,
                 xpos, ymin, +zpos,
                 xpos, ymin, -zpos);
}
#line 22071 "./Test_Fermi.c"
}   /* End of FC_Vitess=Vitess_ChopperFermi() SETTING parameter declarations. */
#undef optimal_wl
#undef omega
#undef pos_ch
#undef Neutrons
#undef Dir
#undef Pos
#undef GeomFilePtr
#undef XFILEName
#undef coef_pi
#undef x_ch
#undef y_ch
#undef phase0
#undef angle_channel
#undef shift_y
#undef main_depth
#undef radius_of_curv
#undef WL
#undef TOF
#undef CurvGeomOption
#undef Option
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'Monitor2_xt'. */
  SIG_MESSAGE("Monitor2_xt (McDisplay)");
  printf("MCDISPLAY: component %s\n", "Monitor2_xt");
#define mccompcurname  Monitor2_xt
#define mccompcurtype  Monitor_nD
#define mccompcurindex 11
#define user1 mccMonitor2_xt_user1
#define user2 mccMonitor2_xt_user2
#define user3 mccMonitor2_xt_user3
#define DEFS mccMonitor2_xt_DEFS
#define Vars mccMonitor2_xt_Vars
#define detector mccMonitor2_xt_detector
#define offdata mccMonitor2_xt_offdata
{   /* Declarations of Monitor2_xt=Monitor_nD() SETTING parameters. */
MCNUM xwidth = mccMonitor2_xt_xwidth;
MCNUM yheight = mccMonitor2_xt_yheight;
MCNUM zdepth = mccMonitor2_xt_zdepth;
MCNUM xmin = mccMonitor2_xt_xmin;
MCNUM xmax = mccMonitor2_xt_xmax;
MCNUM ymin = mccMonitor2_xt_ymin;
MCNUM ymax = mccMonitor2_xt_ymax;
MCNUM zmin = mccMonitor2_xt_zmin;
MCNUM zmax = mccMonitor2_xt_zmax;
MCNUM bins = mccMonitor2_xt_bins;
MCNUM min = mccMonitor2_xt_min;
MCNUM max = mccMonitor2_xt_max;
MCNUM restore_neutron = mccMonitor2_xt_restore_neutron;
MCNUM radius = mccMonitor2_xt_radius;
char* options = mccMonitor2_xt_options;
char* filename = mccMonitor2_xt_filename;
char* geometry = mccMonitor2_xt_geometry;
char* username1 = mccMonitor2_xt_username1;
char* username2 = mccMonitor2_xt_username2;
char* username3 = mccMonitor2_xt_username3;
int nowritefile = mccMonitor2_xt_nowritefile;
#line 496 "/zhome/89/0/38697/McStas/mcstas/2.5/tools/Python/mcrun/../mccodelib/../../../monitors/Monitor_nD.comp"
{
  if (geometry && strlen(geometry) && strcmp(geometry,"0") && strcmp(geometry, "NULL"))
  {
    off_display(offdata);
  } else {
    Monitor_nD_McDisplay(&DEFS, &Vars);
  }
}
#line 22141 "./Test_Fermi.c"
}   /* End of Monitor2_xt=Monitor_nD() SETTING parameter declarations. */
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
/* end of generated C code ./Test_Fermi.c */
