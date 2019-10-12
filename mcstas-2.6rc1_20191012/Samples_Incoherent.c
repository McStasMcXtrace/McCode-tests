/* Automatically generated file. Do not edit. 
 * Format:     ANSI C source code
 * Creator:    McStas <http://www.mcstas.org>
 * Instrument: Samples_Incoherent.instr (Samples_Incoherent)
 * Date:       Sat Oct 12 09:53:48 2019
 * File:       Samples_Incoherent.c
 * Compile:    cc -o Samples_Incoherent.out Samples_Incoherent.c 
 * CFLAGS=
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

#line 712 "Samples_Incoherent.c"

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

#line 945 "Samples_Incoherent.c"

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

#line 4977 "Samples_Incoherent.c"

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

#line 5337 "Samples_Incoherent.c"
#ifdef MC_TRACE_ENABLED
int mctraceenabled = 1;
#else
int mctraceenabled = 0;
#endif
#define MCSTAS "/usr/share/mcstas/2.6rc1/"
int mcdefaultmain = 1;
char mcinstrument_name[] = "Samples_Incoherent";
char mcinstrument_source[] = "Samples_Incoherent.instr";
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
#line 6104 "Samples_Incoherent.c"

/* Shared user declarations for all components 'V_sample'. */
#line 100 "/usr/share/mcstas/2.6rc1/obsolete/V_sample.comp"
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
#line 6121 "Samples_Incoherent.c"

/* Shared user declarations for all components 'PowderN'. */
#line 217 "/usr/share/mcstas/2.6rc1/samples/PowderN.comp"
/* used for reading data table from file */
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

#line 8888 "Samples_Incoherent.c"

/* Shared user declarations for all components 'Single_crystal'. */
#line 296 "/usr/share/mcstas/2.6rc1/samples/Single_crystal.comp"
/* used for reading data table from file */


/* Declare structures and functions only once in each instrument. */
#ifndef SINGLE_CRYSTAL_DECL
#define SINGLE_CRYSTAL_DECL

#ifndef Mosaic_AB_Undefined
#define Mosaic_AB_Undefined {0,0, 0,0,0, 0,0,0}
#endif

struct hkl_data
{
int h,k,l;                  /* Indices for this reflection */
double F2;                  /* Value of structure factor */
double tau_x, tau_y, tau_z; /* Coordinates in reciprocal space */
double tau;                 /* Length of (tau_x, tau_y, tau_z) */
      double u1x, u1y, u1z;       /* First axis of local coordinate system */
      double u2x, u2y, u2z;       /* Second axis of local coordinate system */
      double u3x, u3y, u3z;       /* Third axis of local coordinate system */
      double sig123;              /* The product sig1*sig2*sig3 = volume of spot */
      double m1, m2, m3;          /* Diagonal matrix representation of Gauss */
      double cutoff;              /* Cutoff value for Gaussian tails */
    };

  struct tau_data
    {
      int index;                  /* Index into reflection table */
      double refl;
      double xsect;
      /* The following vectors are in local koordinates. */
      double rho_x, rho_y, rho_z; /* The vector ki - tau */
      double rho;                 /* Length of rho vector */
      double ox, oy, oz;          /* Origin of Ewald sphere tangent plane */
      double b1x, b1y, b1z;       /* Spanning vectors of Ewald sphere tangent */
      double b2x, b2y, b2z;
      double l11, l12, l22;       /* Cholesky decomposition L of 2D Gauss */
      double y0x, y0y;            /* 2D Gauss center in tangent plane */
    };

  struct hkl_info_struct
    {
      struct hkl_data *list;      /* Reflection array */
      int count;                  /* Number of reflections */
      struct tau_data *tau_list;  /* Reflections close to Ewald Sphere */
      double m_delta_d_d;         /* Delta-d/d FWHM */
      double m_ax,m_ay,m_az;      /* First unit cell axis (direct space, AA) */
      double m_bx,m_by,m_bz;      /* Second unit cell axis */
      double m_cx,m_cy,m_cz;      /* Third unit cell axis */
      double asx,asy,asz;         /* First reciprocal lattice axis (1/AA) */
      double bsx,bsy,bsz;         /* Second reciprocal lattice axis */
      double csx,csy,csz;         /* Third reciprocal lattice axis */
      double m_a, m_b, m_c;       /* length of lattice parameter lengths */
      double m_aa, m_bb, m_cc;    /* lattice angles */
      double sigma_a, sigma_i;    /* abs and inc X sect */
      double rho;                 /* density */
      double at_weight;           /* atomic weight */
      double at_nb;               /* nb of atoms in a cell */
      double V0;                  /* Unit cell volume (AA**3) */
      int    column_order[5];     /* column signification [h,k,l,F,F2] */
      int    recip;               /* Flag to indicate if recip or direct cell axes given */
      int    shape;               /* 0:cylinder, 1:box, 2:sphere 3:any shape*/
      int    flag_warning;        /* number of warnings */
      char   type;                /* type of last event: t=transmit,c=coherent or i=incoherent */
      int    h,k,l;               /* last coherent scattering momentum transfer indices */
      int    tau_count;           /* Number of reflections within cutoff */
      double coh_refl, coh_xsect; /* cross section computed with last tau_list */
      double kix, kiy, kiz;       /* last incoming neutron ki */
      int    nb_reuses, nb_refl, nb_refl_count;
    };

  int SX_list_compare (void const *a, void const *b)
  {
     struct hkl_data const *pa = a;
     struct hkl_data const *pb = b;
     double s = pa->tau - pb->tau;

     if (!s) return 0;
     else    return (s < 0 ? -1 : 1);
  } /* PN_list_compare */

  /* ------------------------------------------------------------------------ */
  int
  read_hkl_data(char *SC_file, struct hkl_info_struct *info,
      double SC_mosaic, double SC_mosaic_a, double SC_mosaic_b, double SC_mosaic_c, double *SC_mosaic_AB)
  {
    struct hkl_data *list = NULL;
    int size = 0;
    t_Table sTable; /* sample data table structure from SC_file */
    int i=0;
    double tmp_x, tmp_y, tmp_z;
    char **parsing;
    char flag=0;
    double nb_atoms=1;

    if (!SC_file || !strlen(SC_file) || !strcmp(SC_file,"NULL") || !strcmp(SC_file,"0")) {
      info->count = 0;
      flag=1;
    }
    if (!flag) {
      Table_Read(&sTable, SC_file, 1); /* read 1st block data from SC_file into sTable*/
      if (sTable.columns < 4) {
        fprintf(stderr, "Single_crystal: Error: The number of columns in %s should be at least %d for [h,k,l,F2]\n", SC_file, 4);
        return(0);
      }
      if (!sTable.rows) {
        fprintf(stderr, "Single_crystal: Error: The number of rows in %s should be at least %d\n", SC_file, 1);
        return(0);
      } else size = sTable.rows;

      /* parsing of header */
      parsing = Table_ParseHeader(sTable.header,
        "sigma_abs","sigma_a ",
        "sigma_inc","sigma_i ",
        "column_h",
        "column_k",
        "column_l",
        "column_F ",
        "column_F2",
        "Delta_d/d",
        "lattice_a ",
        "lattice_b ",
        "lattice_c ",
        "lattice_aa",
        "lattice_bb",
        "lattice_cc",
        "nb_atoms","multiplicity",
        NULL);

      if (parsing) {
        if (parsing[0] && !info->sigma_a) info->sigma_a=atof(parsing[0]);
        if (parsing[1] && !info->sigma_a) info->sigma_a=atof(parsing[1]);
        if (parsing[2] && !info->sigma_i) info->sigma_i=atof(parsing[2]);
        if (parsing[3] && !info->sigma_i) info->sigma_i=atof(parsing[3]);
        if (parsing[4])                   info->column_order[0]=atoi(parsing[4]);
        if (parsing[5])                   info->column_order[1]=atoi(parsing[5]);
        if (parsing[6])                   info->column_order[2]=atoi(parsing[6]);
        if (parsing[7])                   info->column_order[3]=atoi(parsing[7]);
        if (parsing[8])                   info->column_order[4]=atoi(parsing[8]);
        if (parsing[9] && info->m_delta_d_d <0) info->m_delta_d_d=atof(parsing[9]);
        if (parsing[10] && !info->m_a)    info->m_a =atof(parsing[10]);
        if (parsing[11] && !info->m_b)    info->m_b =atof(parsing[11]);
        if (parsing[12] && !info->m_c)    info->m_c =atof(parsing[12]);
        if (parsing[13] && !info->m_aa)   info->m_aa=atof(parsing[13]);
        if (parsing[14] && !info->m_bb)   info->m_bb=atof(parsing[14]);
        if (parsing[15] && !info->m_cc)   info->m_cc=atof(parsing[15]);
        if (parsing[16])   nb_atoms=atof(parsing[16]);
        if (parsing[17])   nb_atoms=atof(parsing[17]);
        for (i=0; i<=17; i++) if (parsing[i]) free(parsing[i]);
        free(parsing);
      }
    }

    if (nb_atoms > 1) { info->sigma_a *= nb_atoms; info->sigma_i *= nb_atoms; }

    /* special cases for the structure definition */
    if (info->m_ax || info->m_ay || info->m_az) info->m_a=0; /* means we specify by hand the vectors */
    if (info->m_bx || info->m_by || info->m_bz) info->m_b=0;
    if (info->m_cx || info->m_cy || info->m_cz) info->m_c=0;

    /* compute the norm from vector a if missing */
    if (info->m_ax || info->m_ay || info->m_az) {
      double as=sqrt(info->m_ax*info->m_ax+info->m_ay*info->m_ay+info->m_az*info->m_az);
      if (!info->m_bx && !info->m_by && !info->m_bz) info->m_a=info->m_b=as;
      if (!info->m_cx && !info->m_cy && !info->m_cz) info->m_a=info->m_c=as;
    }
    if (info->m_a && !info->m_b) info->m_b=info->m_a;
    if (info->m_b && !info->m_c) info->m_c=info->m_b;

    /* compute the lattive angles if not set from data file. Not used when in vector mode. */
    if (info->m_a && !info->m_aa) info->m_aa=90;
    if (info->m_aa && !info->m_bb) info->m_bb=info->m_aa;
    if (info->m_bb && !info->m_cc) info->m_cc=info->m_bb;

    /* parameters consistency checks */
    if (!info->m_ax && !info->m_ay && !info->m_az && !info->m_a) {
      fprintf(stderr,
              "Single_crystal: Error: Wrong a lattice vector definition\n");
      return(0);
    }
    if (!info->m_bx && !info->m_by && !info->m_bz && !info->m_b) {
      fprintf(stderr,
              "Single_crystal: Error: Wrong b lattice vector definition\n");
      return(0);
    }
    if (!info->m_cx && !info->m_cy && !info->m_cz && !info->m_c) {
      fprintf(stderr,
              "Single_crystal: Error: Wrong c lattice vector definition\n");
      return(0);
    }
    if (info->m_aa && info->m_bb && info->m_cc && info->recip) {
      fprintf(stderr,
              "Single_crystal: Error: Selecting reciprocal cell and angles is unmeaningful\n");
      return(0);
    }

    /* when lengths a,b,c + angles are given (instead of vectors a,b,c) */
    if (info->m_aa && info->m_bb && info->m_cc)
    {
      double as,bs,cs;
      if (info->m_a) as = info->m_a;
      else as = sqrt(info->m_ax*info->m_ax+info->m_ay*info->m_ay+info->m_az*info->m_az);
      if (info->m_b) bs = info->m_b;
      else bs = sqrt(info->m_bx*info->m_bx+info->m_by*info->m_by+info->m_bz*info->m_bz);
      if (info->m_c) cs = info->m_c;
      else cs =  sqrt(info->m_cx*info->m_cx+info->m_cy*info->m_cy+info->m_cz*info->m_cz);

      info->m_bz = as; info->m_by = 0; info->m_bx = 0;
      info->m_az = bs*cos(info->m_cc*DEG2RAD);
      info->m_ay = bs*sin(info->m_cc*DEG2RAD);
      info->m_ax = 0;
      info->m_cz = cs*cos(info->m_bb*DEG2RAD);
      info->m_cy = cs*(cos(info->m_aa*DEG2RAD)-cos(info->m_cc*DEG2RAD)*cos(info->m_bb*DEG2RAD))
                     /sin(info->m_cc*DEG2RAD);
      info->m_cx = sqrt(cs*cs - info->m_cz*info->m_cz - info->m_cy*info->m_cy);

      printf("Single_crystal: %s structure a=%g b=%g c=%g aa=%g bb=%g cc=%g ",
        (flag ? "INC" : SC_file), as, bs, cs, info->m_aa, info->m_bb, info->m_cc);
    } else {
      if (!info->recip) {
        printf("Single_crystal: %s structure a=[%g,%g,%g] b=[%g,%g,%g] c=[%g,%g,%g] ",
               (flag ? "INC" : SC_file), info->m_ax ,info->m_ay ,info->m_az,
               info->m_bx ,info->m_by ,info->m_bz,
               info->m_cx ,info->m_cy ,info->m_cz);
      } else {
        printf("Single_crystal: %s structure a*=[%g,%g,%g] b*=[%g,%g,%g] c*=[%g,%g,%g] ",
               (flag ? "INC" : SC_file), info->m_ax ,info->m_ay ,info->m_az,
               info->m_bx ,info->m_by ,info->m_bz,
               info->m_cx ,info->m_cy ,info->m_cz);
      }
    }
    /* Compute reciprocal or direct lattice vectors. */
    if (!info->recip) {
      vec_prod(tmp_x, tmp_y, tmp_z,
               info->m_bx, info->m_by, info->m_bz,
               info->m_cx, info->m_cy, info->m_cz);
      info->V0 = fabs(scalar_prod(info->m_ax, info->m_ay, info->m_az, tmp_x, tmp_y, tmp_z));
      printf("V0=%g\n", info->V0);

      info->asx = 2*PI/info->V0*tmp_x;
      info->asy = 2*PI/info->V0*tmp_y;
      info->asz = 2*PI/info->V0*tmp_z;
      vec_prod(tmp_x, tmp_y, tmp_z, info->m_cx, info->m_cy, info->m_cz, info->m_ax, info->m_ay, info->m_az);
      info->bsx = 2*PI/info->V0*tmp_x;
      info->bsy = 2*PI/info->V0*tmp_y;
      info->bsz = 2*PI/info->V0*tmp_z;
      vec_prod(tmp_x, tmp_y, tmp_z, info->m_ax, info->m_ay, info->m_az, info->m_bx, info->m_by, info->m_bz);
      info->csx = 2*PI/info->V0*tmp_x;
      info->csy = 2*PI/info->V0*tmp_y;
      info->csz = 2*PI/info->V0*tmp_z;
    } else {
      info->asx = info->m_ax;
      info->asy = info->m_ay;
      info->asz = info->m_az;
      info->bsx = info->m_bx;
      info->bsy = info->m_by;
      info->bsz = info->m_bz;
      info->csx = info->m_cx;
      info->csy = info->m_cy;
      info->csz = info->m_cz;

      vec_prod(tmp_x, tmp_y, tmp_z,
               info->bsx/(2*PI), info->bsy/(2*PI), info->bsz/(2*PI),
               info->csx/(2*PI), info->csy/(2*PI), info->csz/(2*PI));
      info->V0 = 1/fabs(scalar_prod(info->asx/(2*PI), info->asy/(2*PI), info->asz/(2*PI), tmp_x, tmp_y, tmp_z));
      printf("V0=%g\n", info->V0);

      /*compute the direct cell parameters, ofr completeness*/
      info->m_ax = tmp_x*info->V0;
      info->m_ay = tmp_y*info->V0;
      info->m_az = tmp_z*info->V0;
      vec_prod(tmp_x, tmp_y, tmp_z,info->csx/(2*PI), info->csy/(2*PI), info->csz/(2*PI),info->asx/(2*PI), info->asy/(2*PI), info->asz/(2*PI));
      info->m_bx = tmp_x*info->V0;
      info->m_by = tmp_y*info->V0;
      info->m_bz = tmp_z*info->V0;
      vec_prod(tmp_x, tmp_y, tmp_z,info->asx/(2*PI), info->asy/(2*PI), info->asz/(2*PI),info->bsx/(2*PI), info->bsy/(2*PI), info->bsz/(2*PI));
      info->m_cx = tmp_x*info->V0;
      info->m_cy = tmp_y*info->V0;
      info->m_cz = tmp_z*info->V0;
    }

    if (flag) return(-1);

    if (!info->column_order[0] || !info->column_order[1] || !info->column_order[2]) {
      fprintf(stderr,
              "Single_crystal: Error: Wrong h,k,l column definition\n");
      return(0);
    }
    if (!info->column_order[3] && !info->column_order[4]) {
      fprintf(stderr,
              "Single_crystal: Error: Wrong F,F2 column definition\n");
      return(0);
    }

    /* allocate hkl_data array */
    list = (struct hkl_data*)malloc(size*sizeof(struct hkl_data));

    for (i=0; i<size; i++)
    {
      double h=0, k=0, l=0, F2=0;
      double b1[3], b2[3];
      double sig1, sig2, sig3;

      /* get data from table */
      h = Table_Index(sTable, i, info->column_order[0]-1);
      k = Table_Index(sTable, i, info->column_order[1]-1);
      l = Table_Index(sTable, i, info->column_order[2]-1);
      if (info->column_order[3])
      { F2= Table_Index(sTable, i, info->column_order[3]-1); F2 *= F2; }
      else if (info->column_order[4])
        F2= Table_Index(sTable, i, info->column_order[4]-1);

      list[i].h = h;
      list[i].k = k;
      list[i].l = l;
      list[i].F2 = F2;

      /* Precompute some values */
      list[i].tau_x = h*info->asx + k*info->bsx + l*info->csx;
      list[i].tau_y = h*info->asy + k*info->bsy + l*info->csy;
      list[i].tau_z = h*info->asz + k*info->bsz + l*info->csz;
      list[i].tau = sqrt(list[i].tau_x*list[i].tau_x +
                         list[i].tau_y*list[i].tau_y +
                         list[i].tau_z*list[i].tau_z);
      list[i].u1x = list[i].tau_x/list[i].tau;
      list[i].u1y = list[i].tau_y/list[i].tau;
      list[i].u1z = list[i].tau_z/list[i].tau;
      sig1 = FWHM2RMS*info->m_delta_d_d*list[i].tau;

      /* Find two arbitrary axes perpendicular to tau and each other. */
      normal_vec(b1[0], b1[1], b1[2],
                 list[i].u1x, list[i].u1y, list[i].u1z);
      vec_prod(b2[0], b2[1], b2[2],
               list[i].u1x, list[i].u1y, list[i].u1z,
               b1[0], b1[1], b1[2]);

      /* Find the two mosaic axes perpendicular to tau. */
      if(SC_mosaic > 0) {
        /* Use isotropic mosaic. */
        list[i].u2x = b1[0];
        list[i].u2y = b1[1];
        list[i].u2z = b1[2];
        sig2 = FWHM2RMS*list[i].tau*MIN2RAD*SC_mosaic;
        list[i].u3x = b2[0];
        list[i].u3y = b2[1];
        list[i].u3z = b2[2];
        sig3 = FWHM2RMS*list[i].tau*MIN2RAD*SC_mosaic;
      } else if(SC_mosaic_a > 0 && SC_mosaic_b > 0 && SC_mosaic_c > 0) {
        /* Use anisotropic mosaic. */
        fprintf(stderr,"Single_crystal: Warning: you are using an experimental feature:\n"
          "  anistropic mosaicity. Please examine your data carefully.\n");
        /* compute the jacobian of (tau_v,tau_n) from rotations around the unit cell vectors. */
        struct hkl_data *l =&(list[i]);
        double xia_x,xia_y,xia_z,xib_x,xib_y,xib_z,xic_x,xic_y,xic_z;
        /*input parameters are in arc minutes*/
        double sig_fi_a=SC_mosaic_a*MIN2RAD;
        double sig_fi_b=SC_mosaic_b*MIN2RAD;
        double sig_fi_c=SC_mosaic_c*MIN2RAD;
        if(info->m_a==0) info->m_a=sqrt(scalar_prod( info->m_ax,info->m_ay,info->m_az,info->m_ax,info->m_ay,info->m_az));
        if(info->m_b==0) info->m_b=sqrt(scalar_prod( info->m_bx,info->m_by,info->m_bz,info->m_bx,info->m_by,info->m_bz));
        if(info->m_c==0) info->m_c=sqrt(scalar_prod( info->m_cx,info->m_cy,info->m_cz,info->m_cx,info->m_cy,info->m_cz));

        l->u2x = b1[0];
        l->u2y = b1[1];
        l->u2z = b1[2];
        l->u3x = b2[0];
        l->u3y = b2[1];
        l->u3z = b2[2];

        xia_x=l->tau_x-(M_2_PI*h/info->m_a)*info->asx;
        xia_y=l->tau_y-(M_2_PI*h/info->m_a)*info->asy;
        xia_z=l->tau_z-(M_2_PI*h/info->m_a)*info->asz;
        xib_x=l->tau_x-(M_2_PI*h/info->m_b)*info->bsx;
        xib_y=l->tau_y-(M_2_PI*h/info->m_b)*info->bsy;
        xib_z=l->tau_z-(M_2_PI*h/info->m_b)*info->bsz;
        xic_x=l->tau_x-(M_2_PI*h/info->m_c)*info->csx;
        xic_y=l->tau_y-(M_2_PI*h/info->m_c)*info->csy;
        xic_z=l->tau_z-(M_2_PI*h/info->m_c)*info->csz;

        double xia=sqrt(xia_x*xia_x + xia_y*xia_y + xia_z*xia_z);
        double xib=sqrt(xib_x*xib_x + xib_y*xib_y + xib_z*xib_z);
        double xic=sqrt(xic_x*xic_x + xic_y*xic_y + xic_z*xic_z);

        vec_prod(tmp_x,tmp_y,tmp_z,l->tau_x,l->tau_y,l->tau_z, l->u2x,l->u2y,l->u2z);
        double J_n_fia= xia/info->m_a/l->tau*scalar_prod(info->asx,info->asy,info->asz,tmp_x,tmp_y,tmp_z);
        vec_prod(tmp_x,tmp_y,tmp_z,l->tau_x,l->tau_y,l->tau_z, l->u2x,l->u2y,l->u2z);
        double J_n_fib= xib/info->m_b/l->tau*scalar_prod(info->bsx,info->bsy,info->bsz,tmp_x,tmp_y,tmp_z);
        vec_prod(tmp_x,tmp_y,tmp_z,l->tau_x,l->tau_y,l->tau_z, l->u2x,l->u2y,l->u2z);
        double J_n_fic= xic/info->m_c/l->tau*scalar_prod(info->csx,info->csy,info->csz,tmp_x,tmp_y,tmp_z);

        vec_prod(tmp_x,tmp_y,tmp_z,l->tau_x,l->tau_y,l->tau_z, l->u3x,l->u3y,l->u3z);
        double J_v_fia= xia/info->m_a/l->tau*scalar_prod(info->asx,info->asy,info->asz,tmp_x,tmp_y,tmp_z);
        vec_prod(tmp_x,tmp_y,tmp_z,l->tau_x,l->tau_y,l->tau_z, l->u3x,l->u3y,l->u3z);
        double J_v_fib= xib/info->m_b/l->tau*scalar_prod(info->bsx,info->bsy,info->bsz,tmp_x,tmp_y,tmp_z);
        vec_prod(tmp_x,tmp_y,tmp_z,l->tau_x,l->tau_y,l->tau_z, l->u3x,l->u3y,l->u3z);
        double J_v_fic= xic/info->m_c/l->tau*scalar_prod(info->csx,info->csy,info->csz,tmp_x,tmp_y,tmp_z);

        /*with the jacobian we can compute the sigmas in terms of the orthogonal vectors u2 and u3*/
        sig2=sig_fi_a*fabs(J_v_fia) + sig_fi_b*fabs(J_v_fib) + sig_fi_c*fabs(J_v_fic);
        sig3=sig_fi_a*fabs(J_n_fia) + sig_fi_b*fabs(J_n_fib) + sig_fi_c*fabs(J_n_fic);
      } else if (SC_mosaic_AB[0]!=0 && SC_mosaic_AB[1]!=0){
        if ( (SC_mosaic_AB[2]==0 && SC_mosaic_AB[3]==0 && SC_mosaic_AB[4]==0) || (SC_mosaic_AB[5]==0 && SC_mosaic_AB[6]==0 && SC_mosaic_AB[7]==0) ){
          fprintf(stderr,"Single_crystal: Error: in-plane mosaics are specified but one (or both)\n"
              "  in-plane reciprocal vector is the zero vector\n");
          return(0);
        }
        fprintf(stderr,"Single_crystal: Warning: you are using an experimental feature: \n"
              "  \"in-plane\" anistropic mosaicity. Please examine your data carefully.\n");

        /*for given reflection in list - compute linear comb of tau_a and tau_b*/
        /*check for not in plane - f.i. check if (tau_a X tau_b).tau_i)==0*/
        struct hkl_data *l =&(list[i]);
        double det,c1,c2,sig_tau_c;
        double em_x,em_y,em_z, tmp_x,tmp_y,tmp_z;
        double tau_a[3],tau_b[3];
        /*convert Miller indices to taus*/
        if(info->m_a==0) info->m_a=sqrt(scalar_prod( info->m_ax,info->m_ay,info->m_az,info->m_ax,info->m_ay,info->m_az));
        if(info->m_b==0) info->m_b=sqrt(scalar_prod( info->m_bx,info->m_by,info->m_bz,info->m_bx,info->m_by,info->m_bz));
        if(info->m_c==0) info->m_c=sqrt(scalar_prod( info->m_cx,info->m_cy,info->m_cz,info->m_cx,info->m_cy,info->m_cz));
        tau_a[0]=M_2_PI*( (SC_mosaic_AB[2]/info->m_a)*info->asx + (SC_mosaic_AB[3]/info->m_b)*info->bsx + (SC_mosaic_AB[4]/info->m_c)*info->csx );
        tau_a[1]=M_2_PI*( (SC_mosaic_AB[2]/info->m_a)*info->asy + (SC_mosaic_AB[3]/info->m_b)*info->bsy + (SC_mosaic_AB[4]/info->m_c)*info->csy );
        tau_a[2]=M_2_PI*( (SC_mosaic_AB[2]/info->m_a)*info->asz + (SC_mosaic_AB[3]/info->m_b)*info->bsz + (SC_mosaic_AB[4]/info->m_c)*info->csz );
        tau_b[0]=M_2_PI*( (SC_mosaic_AB[5]/info->m_a)*info->asx + (SC_mosaic_AB[6]/info->m_b)*info->bsx + (SC_mosaic_AB[7]/info->m_c)*info->csx );
        tau_b[1]=M_2_PI*( (SC_mosaic_AB[5]/info->m_a)*info->asy + (SC_mosaic_AB[6]/info->m_b)*info->bsy + (SC_mosaic_AB[7]/info->m_c)*info->csy );
        tau_b[2]=M_2_PI*( (SC_mosaic_AB[5]/info->m_a)*info->asz + (SC_mosaic_AB[6]/info->m_b)*info->bsz + (SC_mosaic_AB[7]/info->m_c)*info->csz );

        /*check determinants to see how we should compute the linear combination of a and b (to match c)*/
        if ((det=tau_a[0]*tau_b[1]-tau_a[1]*tau_b[0])!=0){
          c1= (l->tau_x*tau_b[1] - l->tau_y*tau_b[0])/det;
          c2= (tau_a[0]*l->tau_y - tau_a[1]*l->tau_x)/det;
        }else if ((det=tau_a[1]*tau_b[2]-tau_a[2]*tau_b[1])!=0){
          c1= (l->tau_y*tau_b[2] - l->tau_z*tau_b[1])/det;
          c2= (tau_a[1]*l->tau_z - tau_a[2]*l->tau_y)/det;
        }else if ((det=tau_a[0]*tau_b[2]-tau_a[2]*tau_b[0])!=0){
          c1= (l->tau_x*tau_b[2] - l->tau_z*tau_b[0])/det;
          c2= (tau_a[0]*l->tau_z - tau_a[2]*l->tau_x)/det;
        }
        if ((c1==0) && (c2==0)){
          fprintf(stderr,"Single_crystal: Warning: reflection tau[%i]=(%g %g %g) "
          "has no component in defined mosaic plane\n",
          i, l->tau_x,l->tau_y,l->tau_z);
        }
        /*compute linear combination => sig_tau_i = | c1*sig_tau_a + c2*sig_tau_b |  - also add in the minute to radian scaling factor*/;
        sig_tau_c = MIN2RAD*sqrt(c1*SC_mosaic_AB[0]*c1*SC_mosaic_AB[0] + c2*SC_mosaic_AB[1]*c2*SC_mosaic_AB[1]);
        l->u2x = b1[0]; l->u2y = b1[1]; l->u2z = b1[2];
        l->u3x = b2[0]; l->u3y = b2[1]; l->u3z = b2[2];

        /*so now let's compute the rotation around planenormal tau_a X tau_b*/
        /*g_bar (unit normal of rotation plane) = tau_a X tau_b / norm(tau_a X tau_b)*/
        vec_prod(tmp_x,tmp_y,tmp_z, tau_a[0],tau_a[1],tau_a[2],tau_b[0],tau_b[1],tau_b[2]);
        vec_prod(em_x,em_y,em_z, l->tau_x, l->tau_y, l->tau_z, tmp_x,tmp_y,tmp_z);
        NORM(em_x,em_y,em_z);
        sig2 = l->tau*sig_tau_c*fabs(scalar_prod(em_x,em_y,em_z, l->u2x,l->u2y,l->u2z));
        sig3 = l->tau*sig_tau_c*fabs(scalar_prod(em_x,em_y,em_z, l->u3x,l->u3y,l->u3z));
        /*protect against collapsing gaussians. These seem to be sensible values.*/
        if (sig2<1e-5) sig2=1e-5;
        if (sig3<1e-5) sig3=1e-5;
      }
      else {
        fprintf(stderr,
                "Single_crystal: Error: EITHER mosaic OR (mosaic_a, mosaic_b, mosaic_c)\n"
                "  must be given and be >0.\n");
        return(0);
      }
      list[i].sig123 = sig1*sig2*sig3;
      list[i].m1 = 1/(2*sig1*sig1);
      list[i].m2 = 1/(2*sig2*sig2);
      list[i].m3 = 1/(2*sig3*sig3);
      /* Set Gauss cutoff to 5 times the maximal sigma. */
      if(sig1 > sig2)
        if(sig1 > sig3)
          list[i].cutoff = 5*sig1;
        else
          list[i].cutoff = 5*sig3;
      else
        if(sig2 > sig3)
          list[i].cutoff = 5*sig2;
        else
          list[i].cutoff = 5*sig3;
    }
    Table_Free(&sTable);

    /* sort the list with increasing tau */
    qsort(list, i, sizeof(struct hkl_data),  SX_list_compare);

    info->list = list;
    info->count = i;
    info->tau_list = malloc(i*sizeof(*info->tau_list));
    if(!info->tau_list)
    {
      fprintf(stderr, "Single_crystal: Error: Out of memory!\n");
      return(0);
    }
    return(info->count);
  } /* read_hkl_data */

  /* ------------------------------------------------------------------------ */
  /* hkl_search
    search the HKL reflections which are on the Ewald sphere
    input:
      L,T,count,V0: constants for all calls
      kix,kiy,kiz: may be different for each call
    this function returns:
      tau_count (return), coh_refl, coh_xsect, T (updated elements in the array up to [j])
   */
  int hkl_search(struct hkl_data *L, struct tau_data *T, int count, double V0,
    double kix, double kiy, double kiz, double tau_max,
    double *coh_refl, double *coh_xsect)
  {
    double rho, rho_x, rho_y, rho_z;
    double diff;
    int    i,j;
    double ox,oy,oz;
    double b1x,b1y,b1z, b2x,b2y,b2z, kx, ky, kz, nx, ny, nz;
    double n11, n22, n12, det_N, inv_n11, inv_n22, inv_n12, l11, l22, l12,  det_L;
    double Bt_D_O_x, Bt_D_O_y, y0x, y0y, alpha;

    double ki = sqrt(kix*kix+kiy*kiy+kiz*kiz);

    /* Common factor in coherent cross-section */
    double xsect_factor = pow(2*PI, 5.0/2.0)/(V0*ki*ki);

    for(i = j = 0; i < count; i++)
      {
    /* Assuming reflections are sorted, stop search when max tau exceeded. */
        if(L[i].tau > tau_max)
          break;
        /* Check if this reciprocal lattice point is close enough to the
           Ewald sphere to make scattering possible. */
        rho_x = kix - L[i].tau_x;
        rho_y = kiy - L[i].tau_y;
        rho_z = kiz - L[i].tau_z;
        rho = sqrt(rho_x*rho_x + rho_y*rho_y + rho_z*rho_z);
        diff = fabs(rho - ki);

        /* Check if scattering is possible (cutoff of Gaussian tails). */
        if(diff <= L[i].cutoff)
        {
          /* Store reflection. */
          T[j].index = i;
          /* Get ki vector in local coordinates. */
          kx = kix*L[i].u1x + kiy*L[i].u1y + kiz*L[i].u1z;
          ky = kix*L[i].u2x + kiy*L[i].u2y + kiz*L[i].u2z;
          kz = kix*L[i].u3x + kiy*L[i].u3y + kiz*L[i].u3z;
          T[j].rho_x = kx - L[i].tau;
          T[j].rho_y = ky;
          T[j].rho_z = kz;
          T[j].rho = rho;
          /* Compute the tangent plane of the Ewald sphere. */
          nx = T[j].rho_x/T[j].rho;
          ny = T[j].rho_y/T[j].rho;
          nz = T[j].rho_z/T[j].rho;
          ox = (ki - T[j].rho)*nx;
          oy = (ki - T[j].rho)*ny;
          oz = (ki - T[j].rho)*nz;
          T[j].ox = ox;
          T[j].oy = oy;
          T[j].oz = oz;
          /* Compute unit vectors b1 and b2 that span the tangent plane. */
          normal_vec(b1x, b1y, b1z, nx, ny, nz);
          vec_prod(b2x, b2y, b2z, nx, ny, nz, b1x, b1y, b1z);
          T[j].b1x = b1x;
          T[j].b1y = b1y;
          T[j].b1z = b1z;
          T[j].b2x = b2x;
          T[j].b2y = b2y;
          T[j].b2z = b2z;
          /* Compute the 2D projection of the 3D Gauss of the reflection. */
          /* The symmetric 2x2 matrix N describing the 2D gauss. */
          n11 = L[i].m1*b1x*b1x + L[i].m2*b1y*b1y + L[i].m3*b1z*b1z;
          n12 = L[i].m1*b1x*b2x + L[i].m2*b1y*b2y + L[i].m3*b1z*b2z;
          n22 = L[i].m1*b2x*b2x + L[i].m2*b2y*b2y + L[i].m3*b2z*b2z;
          /* The (symmetric) inverse matrix of N. */
          det_N = n11*n22 - n12*n12;
          inv_n11 = n22/det_N;
          inv_n12 = -n12/det_N;
          inv_n22 = n11/det_N;
          /* The Cholesky decomposition of 1/2*inv_n (lower triangular L). */
          l11 = sqrt(inv_n11/2);
          l12 = inv_n12/(2*l11);
          l22 = sqrt(inv_n22/2 - l12*l12);
          T[j].l11 = l11;
          T[j].l12 = l12;
          T[j].l22 = l22;
          det_L = l11*l22;
          /* The product B^T D o. */
          Bt_D_O_x = b1x*L[i].m1*ox + b1y*L[i].m2*oy + b1z*L[i].m3*oz;
          Bt_D_O_y = b2x*L[i].m1*ox + b2y*L[i].m2*oy + b2z*L[i].m3*oz;
          /* Center of 2D Gauss in plane coordinates. */
          y0x = -(Bt_D_O_x*inv_n11 + Bt_D_O_y*inv_n12);
          y0y = -(Bt_D_O_x*inv_n12 + Bt_D_O_y*inv_n22);
          T[j].y0x = y0x;
          T[j].y0y = y0y;
          /* Factor alpha for the distance of the 2D Gauss from the origin. */
          alpha = L[i].m1*ox*ox + L[i].m2*oy*oy + L[i].m3*oz*oz -
                       (y0x*y0x*n11 + y0y*y0y*n22 + 2*y0x*y0y*n12);
          T[j].refl = xsect_factor*det_L*exp(-alpha)/L[i].sig123; /* intensity of that Bragg */
          *coh_refl += T[j].refl;                                  /* total scatterable intensity */
          T[j].xsect = T[j].refl*L[i].F2;
          *coh_xsect += T[j].xsect;
          j++;
        }

      } /* end for */
        return (j); // this is 'tau_count', i.e. number of reachable reflections
    } /* end hkl_search */

    int hkl_select(struct tau_data *T, int tau_count, double coh_refl, double *sum) {
      int j;
      double r = rand0max(coh_refl);
      *sum = 0;
      for(j = 0; j < tau_count; j++)
      {
        *sum += T[j].refl;
        if(*sum > r) break;
      }
      return j;
    }

    /* Functions for "reorientation", powder and PG modes */
    /* Powder, forward */
    void randrotate(double *nx, double *ny, double *nz, double a, double b, double c) {
      double x1, y1, z1, x2, y2, z2;
      rotate(x1, y1, z1, *nx,*ny,*nz, a, 1, 0, 0); /* <1> = rot(<n>,a) */
      rotate(x2, y2, z2,  x1, y1, z1, b, 0, 1, 0); /* <2> = rot(<1>,b) */
      rotate(*nx,*ny,*nz, x2, y2, z2, c, 0, 0, 1); /* <n> = rot(<2>,c) */
    }
    /* Powder, back */
    void randderotate(double *nx, double *ny, double *nz, double a, double b, double c) {
      double x1, y1, z1, x2, y2, z2;
      rotate(x1, y1, z1, *nx,*ny,*nz, -c, 0,0,1);
      rotate(x2, y2, z2,  x1, y1, z1, -b, 0,1,0);
      rotate(*nx,*ny,*nz, x2, y2, z2, -a, 1,0,0);
    }
    /* PG, forward */
    void PGrotate(double *nx, double *ny, double *nz, double a, double csx, double csy, double csz) {
      /* Currently assumes c-axis along 'x', ought to be generalized... */
      double nvx, nvy, nvz;
      rotate(nvx,nvy,nvz, *nx, *ny, *nz, a, csx, csy, csz);
      *nx = nvx; *ny = nvy; *nz = nvz;
    }
    /* PG, back */
    void PGderotate(double *nx, double *ny, double *nz, double a, double csx, double csy, double csz) {
      /* Currently assumes c-axis along 'x', ought to be generalized... */
      double nvx, nvy, nvz;
      rotate(nvx,nvy,nvz, *nx, *ny, *nz, -a, csx, csy, csz);
      *nx = nvx; *ny = nvy; *nz = nvz;
    }



    /* rotate vector counterclockwise */
    void vec_rotate_2d(double* x, double* y, double angle) {
        double c, s;
        double newx, newy;

        c = cos(angle);
        s = sin(angle);

        newx = *x*c - *y*s;
        newy = *x*s + *y*c;

        *x = newx;
        *y = newy;
    }


#endif /* !SINGLE_CRYSTAL_DECL */
#line 9560 "Samples_Incoherent.c"

/* Shared user declarations for all components 'Isotropic_Sqw'. */
#line 295 "/usr/share/mcstas/2.6rc1/samples/Isotropic_Sqw.comp"

#ifndef ISOTROPIC_SQW
#define ISOTROPIC_SQW $Revision$

/* {j d F2 DW Dd inv2d q F} + { Sq if j == -1}*/
#ifndef Crystallographica
#define Crystallographica { 4,5,7,0,0,0,0, 0,0 }
#define Fullprof          { 4,0,8,0,0,5,0, 0,0 }
#define Undefined         { 0,0,0,0,0,0,0, 0,0 }
#define Lazy              {17,6,0,0,0,0,0,13,0 }
#endif
/* special case for [q,Sq] table */
#define qSq               {-1,0,0,0,0,0,1, 0,0 }




/* For the density of states S(w) */
struct Sqw_W_struct
{
  double omega;        /* omega value for the data block */
  double cumul_proba;  /* cumulated intensity (between 0 and 1) */
};

/* For the S(q|w) probabilities */
struct Sqw_Q_struct
{
   double Q;           /* omega value for the data block */
   double cumul_proba; /* normalized cumulated probability */
};

struct Sqw_Data_struct /* contains normalized Sqw data for probabilities, coh and inc */
{
  struct Sqw_W_struct  *SW;     /* P(w)  ~ density of states */
  struct Sqw_Q_struct **SQW;    /* P(Q|w)= probability of each Q with w */

  long  *SW_lookup;
  long **QW_lookup;
  t_Table Sqw; /* S(q,w) rebin from file in range -w_max:w_max and 0:q_max, with exp(-hw/kT) weight */
  t_Table iqSq;         /* sigma(Ei) = sigma/2/Ki^2 * \int q S(q,w) dq dw up to 2*Ki_max */
  long   q_bins;
  long   w_bins;        /* length of q and w vectors/axes from file */
  double q_max, q_step; /* min=0      */
  double w_max, w_step; /* min=-w_max */
  long   lookup_length;
  char   filename[80];
  double intensity;
  double Ei_max;        /* max neutron incoming energy for Sigma=iqSq table */
  long   iqSq_length;
  char   type;
  double q_min_file;
};

struct Sqw_sample_struct { /* global parameters gathered as a structure */
  char   compname[256];

  struct Sqw_Data_struct Data_inc;
  struct Sqw_Data_struct Data_coh;

  double s_abs, s_coh, s_inc; /* material constants */
  double my_s;
  double my_a_v;
  double mat_rho;
  double mat_weight;
  double mat_density;
  double Temperature; /* temperature from the data file */
  int    shape;       /* 0:cylinder, 1:box, 2:sphere 3:any shape*/

  double sqw_threshold;       /* options to tune S(q,w) */
  double sqw_classical;
  double sqw_norm;

  double barns;               /* for powders */
  double Dd, DWfactor;

  double T2E;                 /* constants */
  char   Q_correction[256];
  double sqSE2K;

  int    maxloop;             /* flags to monitor caught warnings */
  int    minevents;
  long   neutron_removed;
  long   neutron_enter;
  long   neutron_pmult;
  long   neutron_exit;
  char   verbose_output;
  int    column_order[9];     /* column signification */
  long   lookup_length;

  double dq, dw; /* q/w transfer */
  char   type;   /* interaction type: c(coherent),             i(incoherent),
                                      V(isotropic incoherent), t(transmitted) */
  /* store information from the last event */
  double ki_x,ki_y,ki_z,kf_x,kf_y,kf_z;
  double ti, tf;
  double vi, vf;
  double ki, kf;
  double theta;

  double mean_scatt;      /* stat to show at the end */
  double mean_abs;
  double psum_scatt;
  double single_coh;
  double single_inc;
  double multi;

  double rw, rq;
};

#include <stdio.h>
#include <math.h>

/* sets a Data S(q,w) to 'NULL' */
void Sqw_Data_init(struct Sqw_Data_struct *Sqw_Data)
{
  Sqw_Data->q_bins       =0;
  Sqw_Data->w_bins       =0;
  Sqw_Data->q_max        =0;
  Sqw_Data->q_step       =1;
  Sqw_Data->w_max        =0;
  Sqw_Data->w_step       =1;
  Sqw_Data->Ei_max       = 0;
  Sqw_Data->lookup_length=100; /* length of lookup tables */
  Sqw_Data->intensity    =0;
  strcpy(Sqw_Data->filename, "");
  Sqw_Data->SW           =NULL;
  Sqw_Data->SQW          =NULL;
  Sqw_Data->SW_lookup    =NULL;
  Sqw_Data->QW_lookup    =NULL;
  Sqw_Data->iqSq_length  =100;
  Sqw_Data->type         = ' ';
  Sqw_Data->q_min_file   = 0;
}

off_struct offdata;

/* gaussian distribution to appply around Bragg peaks in a powder */
double Sqw_powder_gauss(double x, double mean, double rms) {
  return (exp(-(x-mean)*(x-mean)/(2*rms*rms))/(sqrt(2*PI)*rms));
}

/* Sqw_quantum_correction
*
* Return the 'quantum correction factor Q so that:
*
*   S(q, w) = Q(w) S*(q,w) 
*   S(q,-w) = exp(-hw/kT) S(q,w)
*   S(q, w) = exp( hw/kT) S(q,-w)
*
* with S*=classical limit and Q(w) defined below. For omega > 0, S(q,w) > S(q,-w)
*
* input:
*   w: energy      [meV]
*   T: temperature [K]
*   type: 'Schofield' or 'Boltzmann'        Q = exp(hw/kT/2)
*         'harmonic'  or 'Bader'            Q = hw/kT./(1-exp(-hw/kT))
*         'standard'  or 'Frommhold'        Q = 2./(1+exp(-hw/kT)) [recommended]
*
* References:
*  B. Hehr, http://www.lib.ncsu.edu/resolver/1840.16/7422 PhD manuscript (2010).
*  S. A. Egorov, K. F. Everitt and J. L. Skinner. J. Phys. Chem., 103, 9494 (1999).
*  P. Schofield. Phys. Rev. Lett., 4, 239 (1960).
*  J. S. Bader and B. J. Berne. J. Chem. Phys., 100, 8359 (1994).
*  T. D. Hone and G. A. Voth. J. Chem. Phys., 121, 6412 (2004).
*  L. Frommhold. Collision-induced absorption in gases, 1 st ed., Cambridge
*    Monographs on Atomic, Molecular, and Chemical Physics, Vol. 2,
*    Cambridge Univ. Press: London (1993).
        
 */
double Sqw_quantum_correction(double hw, double T, char *type) {
  double Q   = 1;
  double kT  = T/11.605;  /* [K] -> [meV = 1000*KB/e] */
  if (!hw || !T) return 1;
  if (type == NULL || !strcmp(type, "standard") 
                   || !strcmp(type, "Frommhold") || !strcmp(type, "default"))
    Q = 2/(1+exp(-hw/kT));
  if (!strcmp(type, "Schofield") || !strcmp(type, "Boltzmann"))
    Q = exp(hw/kT/2);
  if (!strcmp(type, "harmonic") || !strcmp(type, "Bader"))
    Q = hw/kT/(1-exp(-hw/kT));
    
  return Q;
}

/*****************************************************************************
* Sqw_read_PowderN: Read PowderN data files
*   Returns t_Table array or NULL in case of error
* Used in : Sqw_readfile (1)
*****************************************************************************/
t_Table *Sqw_read_PowderN(struct Sqw_sample_struct *Sqw, t_Table sqwTable)
{
  struct line_data
  {
    double F2;                  /* Value of structure factor */
    double q;                   /* Q vector */
    int j;                      /* Multiplicity */
    double DWfactor;            /* Debye-Waller factor */
    double w;                   /* Intrinsic line width */
  };
  struct line_data *list = NULL;
  double q_count=0, j_count=0, F2_count=0;
  int    mult_count  =0;
  double q_step      =FLT_MAX;
  long   size        =0;
  int    i, index;
  double q_min=0, q_max=0;
  char   flag=0;
  int    list_count=0;
  double q_step_cur;
  char   flag_qSq = 0;

  t_Table *retTable;

  flag_qSq = (Sqw->column_order[8]>0 && Sqw->column_order[6]>0);

  MPI_MASTER(
  if (Sqw->column_order[0] == 4 && Sqw->barns !=0)
    printf("Isotropic_sqw: %s: Powder file probably of type Crystallographica/Fullprof (lau)\n"
           "WARNING:       but F2 unit is set to powder_barns=1 (barns). Intensity might be 100 times too high.\n",
           Sqw->compname);
  if (Sqw->column_order[0] == 17 && Sqw->barns == 0)
    printf("Isotropic_sqw: %s: Powder file probably of type Lazy Pulver (laz)\n"
           "WARNING:       but F2 unit is set to powder_barns=0 (fm^2). Intensity might be 100 times too low.\n",
           Sqw->compname);
  );
  size = sqwTable.rows;
  MPI_MASTER(
  if (Sqw->verbose_output > 0) {
    printf("Isotropic_sqw: Converting %ld %s from %s into S(q,w) data\n",
        size, flag_qSq ? "S(q)" : "powder lines", sqwTable.filename);
  }
  );
  /* allocate line_data array */
  list = (struct line_data*)malloc(size*sizeof(struct line_data));

  for (i=0; i<size; i++)
    {
      double j=0, d=0, w=0, DWfactor=0, F2=0, Sq=-1, q=0;
      int index;

      if (Sqw->Dd >= 0)      w         = Sqw->Dd;
      if (Sqw->DWfactor > 0) DWfactor  = Sqw->DWfactor;

      /* get data from table using columns {j d F2 DW Dd inv2d q} + { Sq }*/
      /* column indexes start at 1, thus need to substract 1 */
      if (Sqw->column_order[0]>0)
        j = Table_Index(sqwTable, i, Sqw->column_order[0]-1);
      if (Sqw->column_order[1]>0)
        d = Table_Index(sqwTable, i, Sqw->column_order[1]-1);
      if (Sqw->column_order[2]>0)
        F2 = Table_Index(sqwTable, i, Sqw->column_order[2]-1);
      if (Sqw->column_order[3]>0)
        DWfactor = Table_Index(sqwTable, i, Sqw->column_order[3]-1);
      if (Sqw->column_order[4]>0)
        w = Table_Index(sqwTable, i, Sqw->column_order[4]-1);
      if (Sqw->column_order[5]>0)  {
        d = Table_Index(sqwTable, i, Sqw->column_order[5]-1); if (d) d = 1/d/2; }
      if (Sqw->column_order[6]>0)
        q = Table_Index(sqwTable, i, Sqw->column_order[6]-1);
      if (Sqw->column_order[7]>0 && !F2)
        {F2= Table_Index(sqwTable, i, Sqw->column_order[7]-1); F2 *= F2;}

      if (Sqw->column_order[8]>0)
        Sq= Table_Index(sqwTable, i, Sqw->column_order[8]-1);

      if (q > 0 && Sq >= 0) F2 = Sq;
      if (d > 0 && q <= 0)  q = 2*PI/d;

      /* assign and check values */
      j = (j > 0 ? j : 0);
      if (flag_qSq) j=1;
      DWfactor = (DWfactor > 0 ? DWfactor : 1);
      w = (w>0 ? w : 0);
      F2 = (F2 >= 0 ? F2 : 0);
      d = (q > 0 ? 2*PI/d : 0);
      if (j == 0 || d == 0 || q == 0) {
        MPI_MASTER(
        printf("Isotropic_sqw: %s: Warning: line %i has invalid definition\n"
               "         (mult=0 or q=0 or d=0)\n", Sqw->compname, i);
        );
        continue;
      }
      list[list_count].j = j;
      list[list_count].q = q;
      list[list_count].DWfactor = DWfactor;
      list[list_count].w = w;
      list[list_count].F2= F2; /* or S(q) if flag_qSq */

      if (q_max < d) q_max = q;
      if (q_min > d) q_min = q;
      if (list_count > 1) {
        q_step_cur = fabs(list[list_count].q - list[list_count-1].q);
        if (q_step_cur > 1e-5 && (!q_step || q_step_cur < q_step))
         q_step = q_step_cur;
      }

      /* adjust multiplicity if j-column + multiple d-spacing lines */
      /* if  d = previous d, increase line duplication index */
      if (!q_count)     q_count = q;
      if (!j_count)     j_count = j;
      if (!F2_count)    F2_count= F2;
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
        if (Sqw->verbose_output > 2 && (mult_count == list[list_count-1].j
        || (mult_count == list[list_count].j && i == size-1))) {
          MPI_MASTER(
          printf("Isotropic_Sqw: %s: Setting multiplicity to 1 for lines [%i:%i]\n"
                  "         (d-spacing %g is duplicated %i times)\n",
            Sqw->compname, list_count-mult_count, list_count-1, list[list_count-1].q, mult_count);
          );
          for (index=list_count-mult_count; index<list_count; list[index++].j = 1);
          mult_count   = 1;
          q_count = q;
          j_count = j;
          F2_count= F2;
        }
        if (i == size-1) list_count--;
        flag=0;
      }
      list_count++;
    } /* end for */

  /* now builds new Table_Array to continue with Sqw_readfile */
  if (q_max == q_min || !q_step) return(NULL);
  if (!flag_qSq)
    size = 3*q_max/q_step; /* set a default of 3 q values per line */
  else size = list_count;
  /* update the value of q_step */
  q_step = q_max/size;
  MPI_MASTER(
  if (Sqw->verbose_output > 0)
    printf("Isotropic_sqw: q range [%g:%g], creating %li elements vector\n",
        q_min, q_max, size);
  );

  retTable = (t_Table*)calloc(4, sizeof(t_Table));
  if (!retTable) printf("Isotropic_Sqw: ERROR: Cannot allocate PowderN->Sqw table.\n");
  else {
    char *header;
    if (!Table_Init(&retTable[0], size, 1))
      { printf("Isotropic_Sqw: ERROR Cannot allocate q-axis [%li] from Powder lines.\n", size); return(NULL); }
    if (!Table_Init(&retTable[1], 1, 1))
      { printf("Isotropic_Sqw: ERROR Cannot allocate w-axis from Powder lines.\n"); return(NULL); }
    if (!Table_Init(&retTable[2], size, 1))
      { printf("Isotropic_Sqw: ERROR Cannot allocate Sqw [%li] from Powder lines.\n", size); return(NULL); }
    Table_Init(&retTable[3], 0,0);

    header = malloc(64); if (header)
    { retTable[0].header = header; strcpy(retTable[0].header, "q"); }
    header = malloc(64); if (header)
    { retTable[1].header = header; strcpy(retTable[1].header, "w"); }
    header = malloc(64); if (header)
    { retTable[2].header = header; strcpy(retTable[2].header, "Sqw"); }
    for (i=0; i < 4; i++) {
      retTable[i].array_length = 3;
      retTable[i].block_number = i+1;
    }
    if (!flag_qSq)
      for (i=0; i<size; i++)
        retTable[0].data[i]  = i*q_max/size;
    for (i=0; i<list_count; i++) { /* loop on each Bragg peak */
      double peak_qmin, peak_qmax,factor,q;
      if (list[i].w > 0 && !flag_qSq) {
        peak_qmin = list[i].q*(1 - list[i].w*3);
        peak_qmax = list[i].q*(1 + list[i].w*3);
      } else { /* Dirac peak, no width */
        peak_qmin = peak_qmax = list[i].q;
      }
      /* S(q) intensity is here */
      factor = list[i].j*(list[i].DWfactor ? list[i].DWfactor : 1)
               *Sqw->mat_rho*PI/2
               /(Sqw->type == 'c' ? Sqw->s_coh : Sqw->s_inc)*list[i].F2/list[i].q/list[i].q;
      if (Sqw->barns) factor *= 100;
      for (q=peak_qmin; q <= peak_qmax; q += q_step) {
        index = (long)floor(size*q/q_max);
        if (index < 0) index=0;
        else if (index >= size) index = size-1;
        if (flag_qSq) {
          retTable[2].data[index] += list[i].F2;
          retTable[0].data[index]  = list[i].q;
        } else {
          if (list[i].w <=0 || list[i].w*q < q_step) /* step function */
            retTable[2].data[index] += factor/q_step;
          else /* gaussian */
            retTable[2].data[index] += factor
                  * Sqw_powder_gauss(q, list[i].q, list[i].w*list[i].q);
        }
      }
    } /* end for i */
    Table_Stat(&retTable[0]); Table_Stat(&retTable[1]); Table_Stat(&retTable[2]);
    Sqw->sqw_norm = 0; /* F2 are normalized already */
  }

  return(retTable);
} /* Sqw_read_PowderN */

/*****************************************************************************
*  Sqw_search_SW: For a given random number 'randnum', search for the bin
*   containing  the corresponding Sqw->SW
*  Choose an energy in the projected S(w) distribution
* Used in : TRACE (1)
*****************************************************************************/
int Sqw_search_SW(struct Sqw_Data_struct Sqw, double randnum)
{
  int index_w=0;

  if (randnum <0) randnum=0;
  if (randnum >1) randnum=1;

  if (Sqw.w_bins == 1) return(0);
  /* benefit from fast lookup table if exists */
  if (Sqw.SW_lookup) {
    index_w = Sqw.SW_lookup[(long)floor(randnum*Sqw.lookup_length)]-1;
    if (index_w<0) index_w=0;
  }

  while (index_w < Sqw.w_bins && (&(Sqw.SW[index_w]) != NULL) && (randnum > Sqw.SW[index_w].cumul_proba))
      index_w++;
  if (index_w >= Sqw.w_bins) index_w = Sqw.w_bins-1;

  if (&(Sqw.SW[index_w]) == NULL)
  {
      printf("Isotropic_Sqw: Warning: No corresponding value in the SW. randnum too big.\n");
      printf("  index_w=%i ; randnum=%f ; Sqw.SW[index_w-1].cumul_proba=%f (Sqw_search_SW)\n",
            index_w, randnum, Sqw.SW[index_w-1].cumul_proba);
      return index_w-1;
  }
  else
      return (index_w);
}

/*****************************************************************************
*  Sqw_search_Q_proba_per_w: For a given random number randnum, search for
*   the bin containing the corresponding Sqw.SW in the Q probablility grid
*  Choose a momentum in the S(q|w) distribution
*  index is given by Sqw_search_SW
* Used in : TRACE (1)
*****************************************************************************/
int Sqw_search_Q_proba_per_w(struct Sqw_Data_struct Sqw, double randnum, int index_w)
{
  int index_q=0;

  if (randnum <0) randnum=0;
  if (randnum >1) randnum=1;

  /* benefit from fast lookup table if exists */
  if (Sqw.QW_lookup && Sqw.QW_lookup[index_w]) {
    index_q = Sqw.QW_lookup[index_w][(long)floor(randnum*Sqw.lookup_length)]-1;
    if (index_q<0) index_q=0;
  }

  while (index_q < Sqw.q_bins && (&(Sqw.SQW[index_w][index_q]) != NULL)
    && (randnum > Sqw.SQW[index_w][index_q].cumul_proba)) {
      index_q++;
  }
  if (index_q >= Sqw.q_bins) index_q = Sqw.q_bins-1;

  if (&(Sqw.SQW[index_w][index_q]) == NULL)
    return -1;
  else
    return (index_q);
}

/*****************************************************************************
* compute the effective total cross section \int q S(q,w) dw dq
* for incoming neutron energy 0 < Ei < 2*w_max, and
* integration range w=-w_max:Ei and q=Q0:Q1 with
*   Q0 = SE2Q*(sqrt(Ei)-sqrt(Ei-w))=|Ki-Kf|
*   Q1 = SE2Q*(sqrt(Ei)+sqrt(Ei-w))=|Ki+Kf|
* The data to use is Sqw_Data->Sqw, and the limits are Sqw_Data->w_max Sqw_Data->q_max
*   Returns the integral value
* Used in: Sqw_readfile (1)
*****************************************************************************/
double Sqw_integrate_iqSq(struct Sqw_Data_struct *Sqw_Data, double Ei)
{
  long   index_w;
  double iqSq = 0;
  /* w=Ei-Ef q=ki-kf w>0 neutron looses energy, Stokes, Ef = Ei-w > 0, Kf =|Ki-q| > 0 */
  for (index_w=0; index_w < Sqw_Data->w_bins; index_w++) {
    long   index_q;
    double w = -Sqw_Data->w_max + index_w * Sqw_Data->w_step; /* in the Sqw table */
    if (w <= Ei) {       /* integration range w=-w_max:Ei, Ef = Ei-w > 0 */
      double sq=0, Q0=0, Q1=0;
      sq = sqrt(Ei-w);  /* always real as test was true before */
      Q0 = SE2V*V2K*(sqrt(Ei)-sq);
      Q1 = SE2V*V2K*(sqrt(Ei)+sq);

      for (index_q=0; index_q < Sqw_Data->q_bins; index_q++) {
        double q=(double)index_q * Sqw_Data->q_step;
        /* add 'pixel' = q S(q,w) */
        if (Q0 <= q && q <= Q1) iqSq += q*Table_Index(Sqw_Data->Sqw, index_q, index_w);
      }
    }
  }
  /* multiply by 'pixel' size = dq dw */
  return(iqSq * Sqw_Data->q_step * Sqw_Data->w_step);
} /* Sqw_integrate_iqSq */


/*****************************************************************************
* Sqw_diagnosis: Computes Sqw_classical, moments and physical quantities
*                make consistency checks, and output some data files
*   Return: output files and information displayed
* Used in: Sqw_init (2) only by MASTER node with MPI
*****************************************************************************/
void Sqw_diagnosis(struct Sqw_sample_struct *Sqw, struct Sqw_Data_struct *Sqw_Data)
{

  t_Table  Sqw_cl;             /* the Sqw symmetric/classical version (T-> Inf) */
  t_Table  Gqw;                /* the generalized density of states as of Carpenter and Price, J Non Cryst Sol 92 (1987) 153 */
  t_Table  Sqw_moments[7];     /* M0=S(q) M1=E_r M3 w_c w_l M0_cl=S_cl(q) G(w) */
  t_Table  w_c, w_l;
  long     index_q, index_w;
  char     c[CHAR_BUF_LENGTH]; /* temporary variable */
  long     q_min_index = 0;

  char     do_coh=0, do_inc=0;
  double   q_min =0;
  double   u2    =0, S0=1;
  long     u2_count=0;

  if (!Sqw_Data || !Sqw_Data->intensity) return; /* nothing to do with empty S(q,w) */

  if (Sqw_Data->type=='c') do_coh = 1;
  if (Sqw_Data->type=='i') do_inc = 1;

  q_min = Sqw_Data->q_min_file;
  if (q_min <= 0) q_min = Sqw_Data->q_step;

  /* test if there is only one S(q,w) available */
  if (!((Sqw->Data_inc).intensity) || !((Sqw->Data_coh).intensity))
    do_coh = do_inc = 1; /* do both if only one file given */

  if (Sqw->Temperature > 0) {
    if (!Table_Init(&Sqw_cl, Sqw_Data->q_bins, Sqw_Data->w_bins)) {
      printf("Isotropic_Sqw: %s: Cannot allocate S_cl(q,w) Table (%lix%i).\n"
             "WARNING          Skipping S(q,w) diagnosis.\n",
             Sqw->compname, Sqw_Data->q_bins, 1);
      return;
    }
    sprintf(Sqw_cl.filename,
      "S(q,w)_cl from %s (dynamic structure factor, classical)",
      Sqw_Data->filename);
    Sqw_cl.block_number = 1;
    Sqw_cl.min_x        = 0;
    Sqw_cl.max_x        = Sqw_Data->q_max;
    Sqw_cl.step_x       = Sqw_Data->q_step;
  }

  /* initialize moments and 1D stuff */
  for (index_q=0; index_q < 6; index_q++) {
    if (!Table_Init(&Sqw_moments[index_q], Sqw_Data->q_bins, 1)) {
      printf("Isotropic_Sqw: %s: Cannot allocate S(q,w) moment %ld Table (%lix%i).\n"
           "WARNING          Skipping S(q,w) diagnosis.\n",
           Sqw->compname, index_q, Sqw_Data->q_bins, 1);
      Table_Free(&Sqw_cl);
      return;
    }
    Sqw_moments[index_q].block_number = 1;
    Sqw_moments[index_q].min_x  = 0;
    Sqw_moments[index_q].max_x  = Sqw_Data->q_max;
    Sqw_moments[index_q].step_x = Sqw_Data->q_step;
  }
  index_q=6;
  Table_Init(&Sqw_moments[index_q], Sqw_Data->w_bins, 1);
  Sqw_moments[index_q].block_number = 1;
  Sqw_moments[index_q].min_x  = -Sqw_Data->w_max;
  Sqw_moments[index_q].max_x  =  Sqw_Data->w_max;
  Sqw_moments[index_q].step_x =  Sqw_Data->w_step;

  /* set Table titles */
  sprintf(Sqw_moments[0].filename,
    "S(q)=M0(q) from %s [int S(q,w) dw]",
    Sqw_Data->filename);
  sprintf(Sqw_moments[1].filename,
    "M1(q) 1-st moment from %s [int w S(q,w) dw] = HBAR^2*q^2/2/m (f-sum rule, recoil, Lovesey T1 Eq 3.63 p72, Egelstaff p196)",
    Sqw_Data->filename);
  sprintf(Sqw_moments[2].filename,
    "M3(q) 3-rd moment from %s [int w^3 S(q,w) dw] = M1(q)*w_l^2(q)",
    Sqw_Data->filename);
  sprintf(Sqw_moments[3].filename,
    "w_c(q) = sqrt(M1(q)/M0(q)*2kT) collective excitation from %s (Lovesey T1 Eq 5.38 p180, p211 Eq 5.204). Gaussian half-width of the S(q,w) classical",
    Sqw_Data->filename);
  sprintf(Sqw_moments[4].filename,
    "w_l(q) = sqrt(M3(q)/M1(q)) harmonic frequency from %s (Lovesey T1 5.39 p 180)",
    Sqw_Data->filename);
  sprintf(Sqw_moments[5].filename,
    "S_cl(q)=M0_cl(q) from %s [int S_cl(q,w) dw]",
    Sqw_Data->filename);
  sprintf(Sqw_moments[6].filename,
    "G(w) generalized effective density of states from %s (Carpenter J Non Cryst Sol 92 (1987) 153)",
    Sqw_Data->filename);

  for   (index_q=0; index_q < Sqw_Data->q_bins; index_q++) {
    double q           = index_q*Sqw_Data->q_step; /* q value in Sqw_full ; q_min = 0 */
    double sq          = 0;              /* S(q) = w0 = 0-th moment */
    double w1          = 0;              /* first  moment      \int w     Sqw dw */
    double w3          = 0;              /* third  moment      \int w^3   Sqw dw */
    double sq_cl       = 0;              /* S(q) = M0 = 0-th moment classical */
    double w_c         = 0;
    double w_l         = 0;

    for (index_w=0; index_w < Sqw_Data->w_bins; index_w++) {

      double w = -Sqw_Data->w_max + index_w*Sqw_Data->w_step; /* w value in Sqw_full */
      double sqw_cl      =0;
      double sqw_full    =0;

      sqw_full = Table_Index(Sqw_Data->Sqw, index_q, index_w);

      /* Sqw moments */
      if (w && Sqw_Data->w_bins) {
        double tmp;
        tmp  = sqw_full*Sqw_Data->w_step;
        tmp *= w;   w1 += tmp;
        tmp *= w*w; w3 += tmp;
      }

      /* compute classical Sqw and S(q)_cl */
      if (Sqw->Temperature > 0) {
        double n;
        sqw_cl = sqw_full * Sqw_quantum_correction(-w,Sqw->Temperature,Sqw->Q_correction);
        if (!Table_SetElement(&Sqw_cl, index_q, index_w, sqw_cl))
          printf("Isotropic_Sqw: %s: "
                 "Error when setting Sqw_cl[%li q=%g,%li w=%g]=%g from file %s\n",
                 Sqw->compname, index_q, q, index_w, w, sqw_cl, Sqw_Data->filename);
        sq_cl += sqw_cl;
      }
      sq    += sqw_full;
    } /* for index_w */

    sq    *= Sqw_Data->w_step;         /* S(q) = \int S(q,w) dw = structure factor */
    sq_cl *= Sqw_Data->w_step;
    /* find minimal reliable q value (not interpolated) */
    if (q >= q_min && !q_min_index && sq) {
      q_min_index = index_q;
      q_min       = q;
      if (0.9 < sq)
        S0          = sq; /* minimum reliable S(q) */
      else S0 = 1;
    }
    /* compute <u^2> = <3 * ln(S(q)) / q^2> */
    if (q_min_index && q && S0 && sq) {
      u2      += 3 * log(sq/S0) /q/q;
      u2_count++;
    }

    /* store moment values (q) as M0=S(q) M1=E_r M3 w_c w_l M0_cl=S_cl(q) */
    Table_SetElement(&Sqw_moments[0],    index_q, 0, sq);
    Table_SetElement(&Sqw_moments[1],    index_q, 0, w1);
    Table_SetElement(&Sqw_moments[2],    index_q, 0, w3);
    if (w1 > 0 && sq && Sqw->Temperature > 0) {
      double w_c = sqrt(w1/sq*2*Sqw->Temperature*Sqw->T2E);  /* HBAR^2 q^2 kT /m/ S(q) */
      Table_SetElement(&Sqw_moments[3],    index_q, 0, w_c); /* collective dispersion */
    }
    if (w1 && w3*w1 > 0) {
      double w_l = sqrt(w3/w1);
      Table_SetElement(&Sqw_moments[4],    index_q, 0, w_l); /* harmonic dispersion */
    }
    if (Sqw->Temperature > 0)
      Table_SetElement(&Sqw_moments[5],    index_q, 0, sq_cl);

  } /* for index_q */



  /* display some usefull information, only once in MPI mode (MASTER) */
  if (Sqw->Temperature > 0) {
    double Da         = 1.660538921e-27;  /* [kg] unified atomic mass unit = Dalton = 1 g/mol */
  #ifndef KB
    double KB         = 1.3806503e-23;    /* [J/K] */
    /* HBAR   = 1.05457168e-34 */
  #endif
    double CELE       = 1.602176487e-19;  /* [C] Elementary charge CODATA 2006 'e' */
    double meV2Hz     = CELE/HBAR/1000/2/PI; /* 1 meV = 241.80e9 GHz */
    double gqw_sum    = 0;

    /* classical Sqw */
    sprintf(c, "%s_%s_cl.sqw", Sqw->compname, Sqw_Data->type=='c' ? "coh" : "inc");
    Table_Write(Sqw_cl, c, "Momentum [Angs-1]", "'S(q,w)*exp(hw/2kT) classical limit' Energy [meV]",
        0,Sqw_Data->q_max,-Sqw_Data->w_max,Sqw_Data->w_max);
    Table_Free(&Sqw_cl);

    if (u2_count) u2 /= u2_count;

    MPI_MASTER(
    if (do_coh || do_inc)
      printf("Isotropic_Sqw: %s: "
             "Physical constants from the S(q,w) %s for T=%g [K]. Values are estimates.\n",
             Sqw->compname, Sqw_Data->filename, Sqw->Temperature);
    if (do_coh) {
      if (Sqw->mat_weight) {
        double LAMBDA     = HBAR*2*PI/sqrt(2*PI*Sqw->mat_weight*Da*KB*Sqw->Temperature)*1e10;   /* in [Angs] */
        double z          = Sqw->mat_rho * LAMBDA*LAMBDA*LAMBDA;  /* fugacity , rho=N/V in [Angs-3]*/
        double mu         = KB*Sqw->Temperature*log(z);       /* perfect gas chemical potential */
        printf("# De Broglie wavelength     LAMBDA=%g [Angs]\n", LAMBDA);
        printf("# Fugacity                       z=%g (from Egelstaff p32 Eq 2.31)\n", z);
        printf("# Chemical potential            mu=%g [eV] (eq. perfect gas)\n", mu/CELE);
      }

      /* compute isothermal sound velocity and compressibility */
      /* get the S(q_min) value and the corresponding w_c */

      if (q_min_index > 0 && q_min && q_min < 0.6) {
        double w_c = Table_Index(Sqw_moments[3], q_min_index, 0); /* meV */
        /* HBAR = [J*s] */
        double c_T = 2*PI*w_c*meV2Hz/q_min/1e10;                  /* meV*Angs -> m/s */
        double ChiT= S0/(KB*Sqw->Temperature*Sqw->mat_rho*1e30);
        printf("# Isothermal compressibility Chi_T=%g [Pa-1] (Egelstaff  p201 Eq 10.21) at q=%g [Angs-1]\n",
          ChiT, q_min);
        printf("# Isothermal sound velocity    c_T=%g [m/s]  (Lovesey T1 p210 Eq 5.197) at q=%g [Angs-1]\n",
          c_T, q_min);

        /* Computation if C11 is rather tricky as it is obtained from w_l, which is usually quite noisy
         * This means that the obtained values are not reliable from C = rho c_l^2 (Egelstaff Eq 14.10b p284)
         * C44 = rho c_c^2 ~ C11/3
         */
        double w_l = Table_Index(Sqw_moments[4], q_min_index, 0); /* meV */
        double c_l = 2*PI*w_l*meV2Hz/q_min/1e10;                  /* meV*Angs -> m/s */
        double C11 = (Sqw->mat_weight*Da)*(Sqw->mat_rho*1e30)*c_l*c_l;
        printf("# Elastic modulus              C11=%g [GPa]  (Egelstaff Eq 14.10b p284) [rough estimate] at q=%g [Angs-1]\n",
            C11/1e9, q_min);
      }
    }
    if (do_inc) {
      /* display the mean square displacement from S(q) = exp(-<u^2>q^2/3)
           <u^2>= <3 * ln(S(q)) / q^2>
       */
      if (u2_count && u2) {
        printf("# Mean square displacement   <u^2>=%g [Angs^2] (<3 * ln(S(q)) / q^2>)\n", u2);
      }

      /* compute the mean diffusion coefficient D=w_c/q^2 */
      /* FWHM of gaussian is Gamma*RMS2FWHM, only in diffusive regime (Q < 0.2 Angs-1) */
      if (q_min_index > 0 && q_min && q_min < 0.6) {
        double w_c = Table_Index(Sqw_moments[3], q_min_index, 0);
        double D   = 2*PI*w_c*meV2Hz/q_min/q_min/1e14*RMS2FWHM/2; /* meV*Angs^2 -> mm^2/s */
        printf("# Diffusion coefficient          D=%g [mm^2/s] (Egelstaff p220)\n", D);
        if (u2_count && u2 && D)
          printf("# Jump relaxation time         tau=%g [ns] (Egelstaff Eq 11.8 p220)\n", u2*1e-2/6/D);
      }
    }
    ); /* MPI_MASTER */

    /* density of states (generalized) */
    if (!Table_Init(&Gqw, Sqw_Data->q_bins, Sqw_Data->w_bins)) {
      printf("Isotropic_Sqw: %s: Cannot allocate G(q,w) Table (%lix%i).\n"
             "WARNING          Skipping S(q,w) diagnosis.\n",
             Sqw->compname, Sqw_Data->q_bins, 1);
        return;
    }
    sprintf(Gqw.filename,
      "G(q,w) from %s (generalized density of states, Carpenter J Non Cryst Sol 92 (1987) 153)",
      Sqw_Data->filename);
    Gqw.block_number = 1;
    Gqw.min_x        = 0;
    Gqw.max_x        = Sqw_Data->q_max;
    Gqw.step_x       = Sqw_Data->q_step;

    for (index_w=0; index_w < Sqw_Data->w_bins; index_w++) {
      double w        = -Sqw_Data->w_max + index_w*Sqw_Data->w_step; /* w value in Sqw_full */
      double gw       = 0;
      for   (index_q=0; index_q < Sqw_Data->q_bins; index_q++) {
        double q        = index_q*Sqw_Data->q_step; /* q value in Sqw_full ; q_min = 0 */
        double sqw_full = Table_Index(Sqw_Data->Sqw, index_q, index_w);
        double n        = 1/(exp(w/(Sqw->Temperature*Sqw->T2E))-1); /* Bose factor */
        double DW       = q && u2 ? exp(2*u2*q*q/6) : 1;            /* Debye-Weller factor */
        double gqw      = q && n+1 ? sqw_full*DW*2*(Sqw->mat_weight*Da)*w/(n+1)/q/q : 0;
        if (!Table_SetElement(&Gqw, index_q, index_w, gqw))
          printf("Isotropic_Sqw: %s: "
                 "Error when setting Gqw[%li q=%g,%li w=%g]=%g from file %s\n",
                 Sqw->compname, index_q, q, index_w, w, gqw, Sqw_Data->filename);
        gw      += gqw;
        gqw_sum += gqw;
      }
      Table_SetElement(&Sqw_moments[6],    index_w, 0, gw);
    }

    /* normalize the density of states */
    for (index_w=0; index_w < Sqw_Data->w_bins; index_w++) {
      double gw = Table_Index(Sqw_moments[6], index_w, 0);
      Table_SetElement(&Sqw_moments[6], index_w, 0, gw / gqw_sum);
      for   (index_q=0; index_q < Sqw_Data->q_bins; index_q++) {
        double gqw = Table_Index(Gqw, index_q, index_w);
        Table_SetElement(&Gqw, index_q, index_w, gqw / gqw_sum);
      }
    }

    /* write Gqw and free memory */
    if (Sqw_Data->w_bins > 1) {
      sprintf(c, "%s_%s.gqw", Sqw->compname, Sqw_Data->type=='c' ? "coh" : "inc");
        Table_Write(Gqw, c, "Momentum [Angs-1]", "'Generalized density of states' Energy [meV]",
          0,Sqw_Data->q_max,-Sqw_Data->w_max,Sqw_Data->w_max);
      Table_Free(&Gqw);
    }
  } /* if T>0 */

  /* write all tables to disk M0=S(q) M1=E_r M3 w_c w_l M0_cl=S_cl(q) */
  if (Sqw_Data->w_bins > 1) {
    sprintf(c, "%s_%s.m1",  Sqw->compname, Sqw_Data->type=='c' ? "coh" : "inc");
    Table_Write(Sqw_moments[1], c, "Momentum [Angs-1]", "int w S(q,w) dw (recoil) q^2/2m [meV]",
      0,Sqw_Data->q_max,0,0);
    sprintf(c, "%s_%s.w_l", Sqw->compname, Sqw_Data->type=='c' ? "coh" : "inc");
    Table_Write(Sqw_moments[4], c, "Momentum [Angs-1]", "w_l(q) harmonic frequency [meV]",
      0,Sqw_Data->q_max,0,0);
    sprintf(c, "%s_%s.sqw", Sqw->compname, Sqw_Data->type=='c' ? "coh" : "inc");
    Table_Write(Sqw_Data->Sqw, c, "Momentum [Angs-1]", "'S(q,w) dynamical structure factor [meV-1]' Energy [meV]",
      0,Sqw_Data->q_max,-Sqw_Data->w_max,Sqw_Data->w_max);

    if (Sqw->Temperature > 0) {
      sprintf(c, "%s_%s.w_c",    Sqw->compname, Sqw_Data->type=='c' ? "coh" : "inc");
      Table_Write(Sqw_moments[3], c, "Momentum [Angs-1]", "w_c(q) collective excitation [meV]", 0,Sqw_Data->q_max,0,0);
      sprintf(c, "%s_%s_cl.sq",  Sqw->compname, Sqw_Data->type=='c' ? "coh" : "inc");
      Table_Write(Sqw_moments[5], c, "Momentum [Angs-1]", "int S_cl(q,w) dw",
        0,Sqw_Data->q_max,0,0);
      sprintf(c, "%s_%s.gw",  Sqw->compname, Sqw_Data->type=='c' ? "coh" : "inc");
      Table_Write(Sqw_moments[6], c, "Energy [meV]", "'Generalized effective density of states' Energy [meV]",
        -Sqw_Data->w_max,Sqw_Data->w_max,0,0);

    }
  }
  sprintf(c, "%s_%s.sq",    Sqw->compname, Sqw_Data->type=='c' ? "coh" : "inc");
  Table_Write(Sqw_moments[0], c, "Momentum [Angs-1]","S(q) = int S(q,w) dw", 0,Sqw_Data->q_max,0,0);
  sprintf(c, "%s_%s.sigma", Sqw->compname, Sqw_Data->type=='c' ? "coh" : "inc");
  Table_Write(Sqw_Data->iqSq, c, "Energy [meV]", "sigma kf/ki int q S(q,w) dw scattering cross section [barns]", 0,0,0,0);

  /* free Tables */
  for (index_q=0; index_q < 7; Table_Free(&Sqw_moments[index_q++]));

} /* Sqw_diagnosis */

/*****************************************************************************
* Sqw_readfile: Read Sqw data files
*   Returns Sqw_Data_struct or NULL in case of error
* Used in : Sqw_init (2)
*****************************************************************************/
struct Sqw_Data_struct *Sqw_readfile(
  struct Sqw_sample_struct *Sqw, char *file, struct Sqw_Data_struct *Sqw_Data)
{

  t_Table *Table_Array= NULL;
  long     nblocks    = 0;
  char     flag       = 0;

  t_Table  Sqw_full, iqSq; /* the Sqw (non symmetric) and total scattering X section */

  double   sum=0;
  double   mat_at_nb=1;
  double   iq2Sq=0;
  long    *SW_lookup=NULL;
  long   **QW_lookup=NULL;
  char   **parsing  =NULL;

  long   index_q, index_w;
  double q_min_file, q_max_file, q_step_file;
  long   q_bins_file;
  double w_min_file, w_max_file, w_step_file;
  long   w_bins_file;
  double q_max, q_step;
  long   q_bins;
  double w_max, w_step;
  long   w_bins;

  double alpha=0;

  double M1          = 0;
  double M1_cl       = 0;
  double T           = 0;
  double T_file      = 0;
  long   T_count     = 0;
  long   M1_count    = 0;
  long   M1_cl_count = 0;

  /* setup default */
  Sqw_Data_init(Sqw_Data);

  if (!file || !strlen(file) || !strcmp(file, "NULL") || !strcmp(file, "0")) return(Sqw_Data);
  /* read the Sqw file */
  Table_Array = Table_Read_Array(file, &nblocks);
  strncpy(Sqw_Data->filename, file, 80);
  if (!Table_Array) return(NULL);

  /* (1) parsing of header ================================================== */
  parsing = Table_ParseHeader(Table_Array[0].header,
    "Vc","V_0",
    "sigma_abs","sigma_a ",
    "sigma_inc","sigma_i ",
    "column_j", /* 6 */
    "column_d",
    "column_F2",
    "column_DW",
    "column_Dd",
    "column_inv2d", "column_1/2d", "column_sintheta_lambda",
    "column_q", /* 14 */
    "sigma_coh","sigma_c ",
    "Temperature",
    "column_Sq",
    "column_F ", /* 19 */
    "V_rho",
    "density",
    "weight",
    "nb_atoms","multiplicity",
    "classical",
    NULL);
  if (parsing) {
    int i;
    if (parsing[0] && !Sqw->mat_rho)      Sqw->mat_rho    =1/atof(parsing[0]);
    if (parsing[1] && !Sqw->mat_rho)      Sqw->mat_rho    =1/atof(parsing[1]);
    if (parsing[2] && !Sqw->s_abs)    Sqw->s_abs  =  atof(parsing[2]);
    if (parsing[3] && !Sqw->s_abs)    Sqw->s_abs  =  atof(parsing[3]);
    if (parsing[4] && !Sqw->s_inc)    Sqw->s_inc  =  atof(parsing[4]);
    if (parsing[5] && !Sqw->s_inc)    Sqw->s_inc  =  atof(parsing[5]);
    if (parsing[6])                   Sqw->column_order[0]=atoi(parsing[6]);
    if (parsing[7])                   Sqw->column_order[1]=atoi(parsing[7]);
    if (parsing[8])                   Sqw->column_order[2]=atoi(parsing[8]);
    if (parsing[9])                   Sqw->column_order[3]=atoi(parsing[9]);
    if (parsing[10])                  Sqw->column_order[4]=atoi(parsing[10]);
    if (parsing[11])                  Sqw->column_order[5]=atoi(parsing[11]);
    if (parsing[12])                  Sqw->column_order[5]=atoi(parsing[12]);
    if (parsing[13])                  Sqw->column_order[5]=atoi(parsing[13]);
    if (parsing[14])                  Sqw->column_order[6]=atoi(parsing[14]);
    if (parsing[15] && !Sqw->s_coh)   Sqw->s_coh=atof(parsing[15]);
    if (parsing[16] && !Sqw->s_coh)   Sqw->s_coh=atof(parsing[16]);
    if (parsing[17] && !Sqw->Temperature) Sqw->Temperature=atof(parsing[17]); /* from user or file */
    if (parsing[17] )                 T_file=atof(parsing[17]); /* from file */
    if (parsing[18])                  Sqw->column_order[8]=atoi(parsing[18]);
    if (parsing[19])                  Sqw->column_order[7]=atoi(parsing[19]);
    if (parsing[20] && !Sqw->mat_rho)     Sqw->mat_rho    =atof(parsing[20]);
    if (parsing[21] && !Sqw->mat_density) Sqw->mat_density=atof(parsing[21]);
    if (parsing[22] && !Sqw->mat_weight)  Sqw->mat_weight =atof(parsing[22]);
    if (parsing[23] )                 mat_at_nb   =atof(parsing[23]);
    if (parsing[24] )                 mat_at_nb   =atof(parsing[24]);
    if (parsing[25] )                 { /* classical is found in the header */
      char *endptr;
      double value = strtod(parsing[25], &endptr);
      if (*endptr == *parsing[25]) {
        if (Sqw->sqw_classical < 0) Sqw->sqw_classical = 1;
      } else                        Sqw->sqw_classical = value;
    }
    for (i=0; i<=25; i++) if (parsing[i]) free(parsing[i]);
    free(parsing);
  }

  /* compute the scattering unit density from material weight and density */
  /* the weight of the scattering element is the chemical formula molecular weight
   * times the nb of chemical formulae in the scattering element (nb_atoms) */
  if (!Sqw->mat_rho && Sqw->mat_density > 0 && Sqw->mat_weight > 0 && mat_at_nb > 0) {
    /* molar volume [cm^3/mol] = weight [g/mol] / density [g/cm^3] */
    /* atom density per Angs^3 = [mol/cm^3] * N_Avogadro *(1e-8)^3 */
    Sqw->mat_rho = Sqw->mat_density/(Sqw->mat_weight*mat_at_nb)/1e24*NA;
    MPI_MASTER(
    if (Sqw->verbose_output > 0)
      printf("Isotropic_Sqw: %s: Computing scattering unit density V_rho=%g [AA^-3] from density=%g [g/cm^3] weight=%g [g/mol].\n",
        Sqw->compname, Sqw->mat_rho, Sqw->mat_density, Sqw->mat_weight);
    );
  }

  /* the scattering unit cross sections are the chemical formula ones
   * times the nb of chemical formulae in the scattering element */
  if (mat_at_nb > 0) {
    Sqw->s_abs *= mat_at_nb; Sqw->s_inc *= mat_at_nb; Sqw->s_coh *= mat_at_nb;
  }

  if (nblocks) {
    if (nblocks == 1) {
      /* import Powder file */
      t_Table *newTable   = NULL;
      newTable = Sqw_read_PowderN(Sqw, Table_Array[0]);
      if (!newTable) {
        MPI_MASTER(
        printf("Isotropic_Sqw: %s: ERROR importing powder line file %s.\n"
               "               Check format definition.\n",
              Sqw->compname, file);
        );
        exit(-1);
      } else flag=0;
      Table_Free_Array(Table_Array);
      Table_Array = newTable;
    } else if (nblocks != 3) {
      MPI_MASTER(
      printf("Isotropic_Sqw: %s: ERROR "
             "File %s contains %li block%s instead of 3.\n",
              Sqw->compname, file, nblocks, (nblocks == 1 ? "" : "s"));
      );
    } else { flag=0; Sqw->barns=0; /* Sqw files do not use powder_barns */ }
  }

  /* print some info about Sqw files */
  if (flag) Sqw->verbose_output = 2;

  if (flag) {
    MPI_MASTER(
    if (nblocks) printf("ERROR          Wrong file format.\n"
      "               Disabling contribution.\n"
      "               File must contain 3 blocks for [q,w,sqw] or Powder file (1 block, laz,lau).\n");
    );
    return(Sqw_Data);
  }

  sprintf(Table_Array[0].filename, "%s#q",   file);
  sprintf(Table_Array[1].filename, "%s#w",   file);
  sprintf(Table_Array[2].filename, "%s#sqw", file);

  MPI_MASTER(
  if (nblocks && Sqw->verbose_output > 2) {
    printf("Isotropic_Sqw: %s file read, analysing...\n", file);
    Table_Info_Array(Table_Array);
  }
  );

  /* (2) compute range for full +/- w and allocate S(q,w) =================== */

  /* get the q,w extend of the table from the file */
  q_bins_file = Table_Array[0].rows*Table_Array[0].columns;
  w_bins_file = Table_Array[1].rows*Table_Array[1].columns;

  /* is there enough qw data in file to proceed ? */
  if (q_bins_file <= 1 || w_bins_file <= 0) {
    MPI_MASTER(
    printf("Isotropic_Sqw: %s: Data file %s has incomplete q or omega information (%lix%li).\n"
           "ERROR          Exiting.\n",
      Sqw->compname, file, q_bins_file, w_bins_file);
    );
    return(Sqw_Data);
  }

  q_min_file  = Table_Array[0].min_x; q_max_file = Table_Array[0].max_x;
  q_step_file = Table_Array[0].step_x ? Table_Array[0].step_x : (q_max_file - q_min_file)/(Table_Array[0].rows*Table_Array[0].columns);
  w_min_file  = Table_Array[1].min_x; w_max_file = Table_Array[1].max_x;
  w_step_file = Table_Array[1].step_x;

  /* create a regular extended q,w and Sqw tables applying the exp(-hw/kT) factor */
  q_max  = q_max_file;
  q_bins = (q_step_file ?   q_max/q_step_file : q_bins_file)+1;
  q_step = q_bins-1 > 0 ?   q_max/(q_bins-1) : 1;
  w_max  = fabs(w_max_file);
  if (fabs(w_min_file) > fabs(w_max_file)) w_max = fabs(w_min_file);
  /* w_min =-w_max */
  w_bins = (w_step_file ? (long)(2*w_max/w_step_file) : 0)+1; /* twice the initial w range */
  w_step = w_bins-1 > 0 ? 2*w_max/(w_bins-1) : 1;             /* that is +/- w_max         */

  /* create the Sqw table in full range */
  if (!Table_Init(&Sqw_full, q_bins, w_bins)) {
    printf("Isotropic_Sqw: %s: Cannot allocate Sqw_full Table (%lix%li).\n"
           "ERROR          Exiting.\n",
      Sqw->compname, q_bins, w_bins);
    return(NULL);
  }
  sprintf(Sqw_full.filename, "S(q,w) from %s (dynamic structure factor)", file);
  Sqw_full.block_number = 1;

  Sqw_Data->q_bins = q_bins; Sqw_Data->q_max = q_max; Sqw_Data->q_step= q_step;
  Sqw_Data->w_bins = w_bins; Sqw_Data->w_max = w_max; Sqw_Data->w_step= w_step;
  Sqw_Data->q_min_file = q_min_file;

  /* build an energy symmetric Sqw data set with detailed balance there-in, so
   * that we can both compute effective scattering Xsection, probability distributions
   * that is S(q) and \int q S(q).
   * We scan the new Sqw table elements with regular qw binning and search for their
   * equivalent element in the Sqw file data set. This is slower than doing the opposite.
   * We could be scanning all file elements, and fill the new table, but in the
   * process some empty spaces may appear when the initial file binning is not regular
   * in qw, leading to gaps in the new table.
   */

  /* (3) we build q and w lookup table for conversion file -> sqw_full ====== */
  MPI_MASTER(
  if (Sqw->verbose_output > 2)
    printf("Isotropic_Sqw: %s: Creating Sqw_full... (%s, %s)\n",
      Sqw->compname, file, Sqw->type=='c' ? "coh" : "inc");
  );

  double w_file2full[w_bins]; /* lookup table for fast file -> Sqw_full allocation */

  for (index_w=0; index_w < w_bins; w_file2full[index_w++]=0);

  for (index_w=0; index_w < w_bins; index_w++) {

    double w = -w_max + index_w*w_step; /* w value in Sqw_full */
    double index_w_file=0;              /* w index in Sqw file */
    char   found=0;
    for (index_w_file=0; index_w_file < w_bins_file; index_w_file++) {
      double w0=Table_Index(Table_Array[1], (long)index_w_file,  0);
      double w1=Table_Index(Table_Array[1], (long)index_w_file+1,0);
      /* test if we are in Stokes */
      if (w0 > w1) {
        double tmp=w0; w0=w1; w1=tmp;
      }
      if (w0 <= w && w < w1) {
        /* w ~ w_file exists in file, usually on w > 0 side Stokes, neutron looses energy */
        index_w_file += w1-w0 ? (w-w0)/(w1-w0) : 0; /* may correspond with a position in-betwwen two w elements */
        found=1;
        break;
      }
    }
    /* test if we are in anti-Stokes */
    if (!found)
    for (index_w_file=0; index_w_file < w_bins_file; index_w_file++) {
      double w0=Table_Index(Table_Array[1], (long)index_w_file,  0);
      double w1=Table_Index(Table_Array[1], (long)index_w_file+1,0);
      /* test if we are in anti-Stokes */
      if (w0 > w1) {
        double tmp=w0; w0=w1; w1=tmp;
      }
      if (w0 <= -w && -w < w1) {            /* w value is mirrored from the opposite side in file */
        index_w_file += w1-w0 ? (-w-w0)/(w1-w0) : 0;
        index_w_file = -index_w_file;               /* in this case, index value is set to negative */
        break;
      }
    }
    w_file2full[index_w] = index_w_file;
  }

  double q_file2full[q_bins];
  for   (index_q=0; index_q < q_bins; q_file2full[index_q++]=0);

  for (index_q=0; index_q < q_bins; index_q++) {

    double q           = index_q*q_step; /* q value in Sqw_full ; q_min = 0 */
    double index_q_file= 0;              /* q index in Sqw file */

    /* search for q value in the initial file data set */
    if      (q <= q_min_file) index_q_file=0;
    else if (q >= q_max_file) index_q_file=q_bins_file-1;
    else
    for (index_q_file=0; index_q_file < q_bins_file; index_q_file++) {
      double q0=Table_Index(Table_Array[0], (long)index_q_file,  0);
      double q1=Table_Index(Table_Array[0], (long)index_q_file+1,0);
      if (q0 <= q && q <= q1) {
        index_q_file += q1-q0 ? (q-q0)/(q1-q0) : 0; /* may correspond with a position in-betwwen two q elements */
        break;
      }
    }
    q_file2full[index_q] = index_q_file;
  }

  /* (4) now we build Sqw on full Q,W ranges, using the Q,W table lookup above -> Sqw_full */
  for (index_q=0; index_q < q_bins; index_q++) {

    double q           = index_q*q_step; /* q value in Sqw_full ; q_min = 0 */
    double index_q_file= 0;              /* q index in Sqw file */

    /* get q value in the initial file data set */
    index_q_file = q_file2full[index_q];

    /* now scan energy elements in Sqw full, and search these in file data */
    for (index_w=0; index_w < w_bins; index_w++) {
      double w = -w_max + index_w*w_step; /* w value in Sqw_full */
      double index_w_file=0;              /* w index in Sqw file */
      double sqw_file    =0;              /* Sqw(index_q, index_w) value interpolated from file */

      /* search for w value in the file data set, negative when mirrored */
      index_w_file = w_file2full[index_w];
      /* get Sqw_file element, with bi-linear interpolation from file */
      /* when the initial file does not contain the energy, the opposite element (-w) is used */
      sqw_file     = Table_Value2d(Table_Array[2], index_q_file, fabs(index_w_file));
      /* apply the minimum threshold to remove noisy background in S(q,w) */
      if (sqw_file < Sqw->sqw_threshold) sqw_file = 0;
      else if (index_w_file < 0)         sqw_file = -sqw_file; /* negative == mirrored from other side */

      if (!Table_SetElement(&Sqw_full, index_q, index_w, sqw_file))
        printf("Isotropic_Sqw: %s: "
               "Error when setting Sqw[%li q=%g,%li w=%g]=%g from file %s\n",
               Sqw->compname, index_q, q, index_w, w, fabs(sqw_file), file);
    } /* for index_w */
  } /* for index_q */

  /* free memory and store limits for new full Sqw table */
  Table_Free_Array(Table_Array);

  /* if only one S(q,w) side is given, it is symmetrised by mirroring, then M1=0 */

  /* (5) test if the Sqw_full is classical or not by computing the 1st moment (=0 for classical) */
  /* also compute temperature (quantum case) from file if not set */
  for (index_q=0; index_q < q_bins; index_q++) {

    double q           = index_q*q_step; /* q value in Sqw_full ; q_min = 0 */

    for (index_w=0; index_w < w_bins; index_w++) {
      double w        = -w_max + index_w*w_step; /* w value in Sqw_full */
      double sqw_full = Table_Index(Sqw_full, index_q, index_w);
      long   index_mw = w_bins-1-index_w;        /* opposite w index in S(q,w) */
      double sqw_opp  = Table_Index(Sqw_full, index_q, index_mw);
      double T_defined= T_file ? T_file : Sqw->Temperature; /* T better from file, else from user */

      /* the analysis must be done only on values which exist on both sides */
      /* as integrals must be symmetric, and Bose factor requires both sides as well */
      if (sqw_full > 0 && sqw_opp > 0) {
        /* compute temperature from Bose factor */
        if (sqw_opp != sqw_full) {
          T      += fabs(w/log(sqw_opp/sqw_full)/Sqw->T2E);
          T_count++;
        }
        /* we first assume Sqw is quantum. M1_cl should be 0, M1 should be recoil */
        M1      += w*sqw_full*w_step;
        M1_count++;
        /* we assume it is quantum (non symmetric) and check that its symmetrized version has M1_cl=0 */
        if (T_defined > 0) {
          sqw_opp = sqw_full * Sqw_quantum_correction(-w, T_defined,Sqw->Q_correction); /* Sqw_cl */
          M1_cl      += w*sqw_opp*w_step;
          M1_cl_count++;
        } else if (Sqw->mat_weight) {
          /* T=0 ? would compute the M1_cl = M1 - recoil energy, assuming we have a quantum S(q,w) in file */
          /* the M1(quantum) = (MNEUTRON/m)*2.0725*q^2 recoil energy */
          double Da = 1.660538921e-27; /* atomic mass unit */
          double Er = (MNEUTRON/Sqw->mat_weight/Da)*2.0725*q*q; /* recoil for one scattering unit in the cell [meV] Schober JDN16 p239 */
          M1_cl      += M1 - Er;
          M1_cl_count++;
        }
      } /* both side from file */
    } /*index_w */
  } /*index_q */

  if (T_count)     T     /= T_count;     /* mean temperature from Bose ratio */
  if (M1_count)    M1    /= M1_count;
  if (M1_cl_count) M1_cl /= M1_cl_count; /* mean energy value along q range */

  /* determine if we use a classical or quantum S(q,w) */
  if (Sqw->sqw_classical < 0) {
    if (fabs(M1) < 2*w_step) {
      Sqw->sqw_classical = 1; /* the initial Sqw from file seems to be centered, thus classical */
    } else if (fabs(M1_cl) < fabs(M1)) {
      /* M1 for classical is closer to 0 than for quantum one */
      Sqw->sqw_classical = 0; /* initial data from file seems to be quantum (non classical) */
    } else { /* M1_cl > M1 > 2*w_step */
      MPI_MASTER(
      printf("Isotropic_Sqw: %s: I do not know if S(q,w) data is classical or quantum.\n"
             "WARNING        First moment M1=%g M1_cl=%g for file %s. Defaulting to classical case.\n",
                 Sqw->compname, M1, M1_cl, file);
      );
    }
  }
  if (Sqw->sqw_classical < 0) Sqw->sqw_classical=1; /* default when quantum/classical choice is not set */
  
  /* (5b) set temperature. Temperatures defined are:
  *   T_file:           temperature read from the .sqw file
  *   T:                temperature computed from the quantum Sqw using detailed balance
  *   Sqw->Temperature: temperature defined by user, or read from file when not set
  */
 
  
  /* display some warnings about the computed temperature */
  if (T > 3000) T=0; /* unrealistic */
  if (!T_file && T) {
    MPI_MASTER(
    if (Sqw->verbose_output > 0) {
      printf(  "Isotropic_Sqw: %s: Temperature computed from S(q,w) data from %s is T=%g [K].\n",
        Sqw->compname, file, T);
    );
    }
  }
  
  if (Sqw->Temperature == 0) {
    Sqw->Temperature = T_file ? T_file : T; /* 0:  not set: we use file value, else computed */
  } else if (Sqw->Temperature ==-1) {
    Sqw->Temperature = 0;                   /* -1: no detailed balance. Display message at end of INIT */
  } else if (Sqw->Temperature ==-2) {
    Sqw->Temperature = T ? T : T_file;      /* -2: use guessed value when available */
  } /* else let value as it is (e.g. >0) */
  
  if (Sqw->verbose_output > 0 && Sqw->Temperature) {
    MPI_MASTER(
    printf(  "Isotropic_Sqw: %s: Temperature set to T=%g [K]\n", Sqw->compname, Sqw->Temperature);
    );
  }

  MPI_MASTER(
  if (Sqw->verbose_output > 0 && w_bins > 1)
    printf("Isotropic_Sqw: %s: S(q,w) data from %s (%s) assumed to be %s.\n",
      Sqw->compname, file, Sqw->type=='c' ? "coh" : "inc",
      Sqw->sqw_classical ? "classical (symmetrised in energy)" : "non-classical (includes Bose factor, non symmetric in energy)");
  );

  /* (6) apply detailed balance on Sqw_full for classical or quantum S(q,w) */
  /* compute the \int q^2 S(q) for normalisation */
  
  MPI_MASTER(
  if (Sqw->sqw_classical && Sqw->verbose_output > 0 && Sqw->Temperature > 0)
    printf("Isotropic_Sqw: %s: Applying exp(hw/2kT) factor with T=%g [K] on %s file (classical/symmetric) using '%s' quantum correction\n",
      Sqw->compname, Sqw->Temperature, file, Sqw->Q_correction);
  );
  for   (index_q=0; index_q < q_bins; index_q++) {
    double sq          = 0;
    double q           = index_q*q_step;  /* q value in Sqw_full ; q_min = 0 */
    for (index_w=0; index_w < w_bins; index_w++) {
      double w = -w_max + index_w*w_step; /* w value in Sqw_full */
      double balance   = 1;               /* detailed balance factor, default is 1 */
      double sqw_full  = Table_Index(Sqw_full, index_q, index_w);

      /* do we use a symmetric S(q,w) from real G(r,t) from e.g. MD ? */
      
      if (Sqw->sqw_classical && Sqw->Temperature > 0) /* data is symmetric, we apply Bose factor */
        /* un-symmetrize Sqw(file) * exp(hw/kT/2) on both sides */
        balance = Sqw_quantum_correction(w, Sqw->Temperature, Sqw->Q_correction);
      else if (!Sqw->sqw_classical) {  /* data is quantum (contains Bose) */
        if (sqw_full < 0) { /* quantum but mirrored/symmetric value (was missing in file) */
          if (T)
            /* prefer to use T computed from file for mirroring */
            balance *= exp(w/(T*Sqw->T2E));                /* apply Bose where missing */
          else if (Sqw->Temperature)
            balance *= exp(w/(Sqw->Temperature*Sqw->T2E)); /* apply Bose where missing */
        }
        /* test if T computed from file matches requested T, else apply correction */
        if (T && Sqw->Temperature > 0 && Sqw->Temperature != T) {
          balance *= exp(-w/(T*Sqw->T2E)/2);                /* make it a classical data set: remove computed/read T from quantum data file */
          balance *= exp( w/(Sqw->Temperature*Sqw->T2E)/2); /* then apply Bose to requested temperature */
        }
      }
      
      /* update Sqw value */
      sqw_full = fabs(sqw_full)*balance;
      Table_SetElement(&Sqw_full, index_q, index_w, sqw_full);
      /* sum up the S(q) (0-th moment) = w0 */
      sq       += sqw_full;
    } /* index_w */
    sq    *= w_step;         /* S(q)  = \int S(q,w) dw    = structure factor */
    iq2Sq += q*q*sq*q_step;  /* iq2Sq = \int q^2 S(q) dq  = used for g-sum rule (normalisation) */
    sum   += sq*q_step;      /* |S|   = \int S(q,w) dq dw = total integral value in file */
  } /* index_q */

  if (!sum) {
    MPI_MASTER(
    printf("Isotropic_Sqw: %s: No valid data in the selected (Q,w) range: sum(S)=0\n"
           "ERROR          Available Sqw data is\n",
      Sqw->compname);
    printf("                 q=[%g:%g] w=[%g:%g]\n",
           q_min_file, q_max_file,
           w_min_file, w_max_file);
    );
    Table_Free(&Sqw_full);
    return(NULL);
  }

  /* norm S(q ,w) = sum(S)*q_range/q_bins on total q,w range from file */
  sum *= (q_max_file - q_min_file)/q_bins_file;

  /* (7) renormalization of S(q,w) ========================================== */

  if      ( Sqw->sqw_norm >0) alpha=Sqw->sqw_norm;
  else if (!Sqw->sqw_norm)    alpha=1;

  if (!alpha && iq2Sq) { /* compute theoretical |S| norm */
    /* Eq (2.44) from H.E. Fischer et al, Rep. Prog. Phys., 69 (2006) 233 */
    alpha =
      (q_max*q_max*q_max/3 - (Sqw->type == 'c' ? 2*PI*PI*Sqw->mat_rho : 0))
      /iq2Sq;
  }

  if (alpha < 0) {
    MPI_MASTER(
    printf("Isotropic_Sqw: %s: normalisation factor is negative. rho=%g [Angs^-3] may be too high.\n"
           "WARNING        Disabling renormalization i.e. keeping initial S(q,w).\n",
      Sqw->compname, Sqw->mat_rho);
    );
    alpha=0;
  }

  /* apply normalization on S(q,w) */
  if (alpha && alpha != 1) {
    sum *= alpha;
    for (index_q=0; index_q < q_bins ; index_q++) {
      for (index_w=0; index_w < w_bins; index_w++)
        Table_SetElement(&Sqw_full, index_q, index_w,
          Table_Index(Sqw_full, index_q, index_w) * alpha);
    }
  }

  Sqw_Data->intensity       = sum;

  Table_Stat(&Sqw_full);
  Sqw_full.min_x        = 0;
  Sqw_full.max_x        = q_max;
  Sqw_full.step_x       = q_step;

  MPI_MASTER(
  if (Sqw->verbose_output > 0) {
    printf("Isotropic_Sqw: %s: Generated %s %scoherent Sqw\n"
           "                   q=[%g:%g Angs-1] w=[%g:%g meV] |S|=%g size=[%lix%li] sigma=%g [barns]\n",
           Sqw->compname, file, (Sqw->type == 'i' ? "in" : ""),
           q_min_file, q_max_file,
           w_min_file, w_max_file, Sqw_Data->intensity,
           q_bins, Sqw_Data->w_bins,
           (Sqw->type == 'i' ? Sqw->s_inc : Sqw->s_coh));
    if (w_max < 1e-2)
      printf("               Mainly elastic scattering.\n");
    if (Sqw->sqw_norm >0 && Sqw->sqw_norm != 1)
      printf("                   normalization factor S(q,w)*%g (user)\n", alpha);
    else if (Sqw->sqw_norm<0)
      printf("                   normalization factor S(q,w)*%g (auto) \\int q^2 S(q) dq=%g\n", alpha, iq2Sq);
  }
  );

  /* (8) Compute total cross section ======================================== */

  /* now compute the effective total cross section  (Sqw_integrate_iqSq)
        sigma(Ei) = sigma/2/Ki^2 * \int q S(q,w) dw dq
   * for each incoming neutron energy 0 < Ei < 2*w_max, and
   * integration range w=-Ei:w_max and q=Q0:Q1 with
   *   Q0 = SE2Q*(sqrt(E)-sqrt(E+w))
   *   Q1 = SE2Q*(sqrt(E)+sqrt(E+w))
   */

  Sqw_Data->lookup_length = Sqw->lookup_length;
  Sqw_Data->iqSq_length   = Sqw->lookup_length;
  /* this length should be greater when w_max=0 for e.g. elastic only */
  if (w_bins <= 1) Sqw_Data->iqSq_length = q_bins;

  if (!Table_Init(&iqSq, Sqw_Data->iqSq_length, 1)) {
    /* from 0 to 2*Ki_max */
    printf("Isotropic_Sqw: %s: Cannot allocate [int q S(q,w) dq dw] array (%li bytes).\n"
           "ERROR          Exiting.\n",
      Sqw->compname, Sqw->lookup_length*sizeof(double));
    Table_Free(&Sqw_full);
    return(NULL);
  }

  /* compute the maximum incoming energy that can be handled */
  Sqw_Data->Ei_max = 2*w_max;

  /* Checked in different ways in Powder and "proper inelastic" case */
  if (w_step==1) {
    /* Powder */
    double Ei_max_q = (q_max*K2V)*(q_max*K2V)*VS2E/2;
    if (Ei_max_q > Sqw_Data->Ei_max) Sqw_Data->Ei_max = Ei_max_q;
  } else {
    /* Proper inelastic */
    /* check if the q-range will limit the integration */
    if ((q_max*K2V)*(q_max*K2V)*VS2E/2 > Sqw_Data->Ei_max) {
      /* then scan Ei until we pass q_max */
      for (index_w=0; index_w < Sqw_Data->iqSq_length; index_w++) {
	double Ei = index_w*2*w_max/Sqw_Data->iqSq_length;
	if ( (Ei > w_max && sqrt(Ei)+sqrt(Ei-w_max) >= q_max/(SE2V*V2K))
	     || sqrt(Ei)+sqrt(Ei+w_max) >= q_max/(SE2V*V2K))
	  if (Ei < Sqw_Data->Ei_max) {
	    Sqw_Data->Ei_max = Ei;
	    break;
	  }
      }
    }
  }

  MPI_MASTER(
  if (Sqw->verbose_output >= 2)
    printf("Isotropic_Sqw: %s: Creating Sigma(Ei=0:%g [meV]) with %li entries...(%s %s)\n",
      Sqw->compname, Sqw_Data->Ei_max, Sqw_Data->iqSq_length, file, Sqw->type=='c' ? "coh" : "inc");
  );
  Sqw_Data->Sqw  = Sqw_full; /* store the S(q,w) Table (matrix) for Sqw_integrate_iqSq */

  /* this loop takes time: use MPI when available */

  for (index_w=0; index_w < Sqw_Data->iqSq_length; index_w++) {

    /* Ei = energy of incoming neutron, typically 0:w_max or 0:2*q_max */
    double Ei; /* neutron energy value in Sqw_full, up to 2*w_max */
    double Ki, Vi;
    double Sigma=0;
    Ei = index_w*Sqw_Data->Ei_max/Sqw_Data->iqSq_length;
    Vi = sqrt(Ei/VS2E);
    Ki = V2K*Vi;
    /* sigma(Ei) = sigma/2/Ki^2 * \int q S(q,w) dq dw */
    /* Eq (6) from E. Farhi et al. J. Comp. Phys. 228 (2009) 5251 */
    Sigma = Ki <= 0 ? 0 : (Sqw->type=='c' ? Sqw->s_coh : Sqw->s_inc)
                          /2/Ki/Ki * Sqw_integrate_iqSq(Sqw_Data, Ei);
    Table_SetElement(&iqSq, index_w, 0, Sigma );
  }

  sprintf(iqSq.filename, "[sigma/2Ki^2 int q S(q,w) dq dw] from %s", file);
  iqSq.min_x  = 0;
  iqSq.max_x  = Sqw_Data->Ei_max;
  iqSq.step_x = Sqw_Data->Ei_max/Sqw_Data->iqSq_length;
  iqSq.block_number = 1;

  Sqw_Data->iqSq = iqSq;     /* store the sigma(Ei) = \int q S(q,w) dq dw Table (vector) */

  /* (9) Compute P(w) probability =========================================== */

  /* set up 'density of states' */
  /* uses: Sqw_full and w axes */
  Sqw_Data->SW =
    (struct Sqw_W_struct*)calloc(w_bins, sizeof(struct Sqw_W_struct));

  if (!Sqw_Data->SW) {
    printf("Isotropic_Sqw: %s: Cannot allocate SW (%li bytes).\n"
           "ERROR          Exiting.\n",
      Sqw->compname, w_bins*sizeof(struct Sqw_W_struct));
    Table_Free(&Sqw_full);
    Table_Free(&iqSq);
    return(NULL);
  }
  sum = 0;
  for (index_w=0; index_w < w_bins ; index_w++) {
    double local_val = 0;
    double w         = -w_max + index_w * w_step;
    for (index_q=0; index_q < q_bins ; index_q++) { /* integrate on all q values */
      local_val += Table_Index(Sqw_full, index_q, index_w)*q_step*index_q*q_step; /* q*S(q,w) */
    }
    Sqw_Data->SW[index_w].omega = w;
    sum                  += local_val; /* S(w)=\int S(q,w) dq */
    /* compute cumulated probability */
    Sqw_Data->SW[index_w].cumul_proba = local_val;
    if (index_w) Sqw_Data->SW[index_w].cumul_proba += Sqw_Data->SW[index_w-1].cumul_proba;
    else         Sqw_Data->SW[index_w].cumul_proba = 0;
  }

  if (!sum) {
    MPI_MASTER(
    printf("Isotropic_Sqw: %s: Total S(q,w) intensity is NULL.\n"
           "ERROR          Exiting.\n", Sqw->compname);
    );
    Table_Free(&Sqw_full);
    Table_Free(&iqSq);
    return(NULL);
  }

  /* normalize Pw distribution to 0:1 range */
  for (index_w=0; index_w < w_bins ; index_w++) {
    Sqw_Data->SW[index_w].cumul_proba /= Sqw_Data->SW[w_bins-1].cumul_proba;
  }

  if (Sqw->verbose_output > 2) {
    MPI_MASTER(
    printf("Isotropic_Sqw: %s: Generated normalized SW[%li] in range [0:%g]\n",
      Sqw->compname, w_bins, Sqw_Data->SW[w_bins-1].cumul_proba);
    );
  }

  /* (10) Compute P(Q|w) probability ======================================== */

  /* set up Q probability table per w bin */
  /* uses:  Sqw_full */
  Sqw_Data->SQW =
    (struct Sqw_Q_struct**)calloc(w_bins, sizeof(struct Sqw_Q_struct*));

  if (!Sqw_Data->SQW) {
    printf("Isotropic_Sqw: %s: Cannot allocate P(Q|w) array (%li bytes).\n"
           "ERROR          Exiting.\n",
      Sqw->compname, w_bins*sizeof(struct Sqw_Q_struct*));
    Table_Free(&Sqw_full);
    Table_Free(&iqSq);
    return(NULL);
  }
  for (index_w=0; index_w < w_bins ; index_w++) {
    Sqw_Data->SQW[index_w]=
        (struct Sqw_Q_struct*)calloc(q_bins, sizeof(struct Sqw_Q_struct));

    if (!Sqw_Data->SQW[index_w]) {
      printf("Isotropic_Sqw: %s: Cannot allocate P(Q|w)[%li] (%li bytes).\n"
             "ERROR          Exiting.\n",
        Sqw->compname, index_w, q_bins*sizeof(struct Sqw_Q_struct));
      Table_Free(&Sqw_full);
      Table_Free(&iqSq);
      return(NULL);
    }
    /* set P(Q|W) and compute total intensity */
    for (index_q=0; index_q < q_bins ; index_q++) {
      double q  = index_q * q_step;
      Sqw_Data->SQW[index_w][index_q].Q     = q;

      /* compute cumulated probability and take into account Jacobian with additional 'q' factor */
      Sqw_Data->SQW[index_w][index_q].cumul_proba = q*Table_Index(Sqw_full, index_q, index_w); /* q*S(q,w) */
      if (index_q) Sqw_Data->SQW[index_w][index_q].cumul_proba += Sqw_Data->SQW[index_w][index_q-1].cumul_proba;
      else Sqw_Data->SQW[index_w][index_q].cumul_proba = 0;
    }
    /* normalize P(q|w) distribution to 0:1 range */
    for (index_q=0; index_q < q_bins ;
    	Sqw_Data->SQW[index_w][index_q++].cumul_proba /= Sqw_Data->SQW[index_w][q_bins-1].cumul_proba
    );

  }
  if (Sqw->verbose_output > 2) {
    MPI_MASTER(
    printf("Isotropic_Sqw: %s: Generated P(Q|w)\n",
      Sqw->compname);
    );
  }

  /* (11) generate quick lookup tables for SW and SQW ======================= */

  SW_lookup = (long*)calloc(Sqw->lookup_length, sizeof(long));

  if (!SW_lookup) {
    printf("Isotropic_Sqw: %s: Cannot allocate SW_lookup (%li bytes).\n"
           "Warning        Will be slower.\n",
      Sqw->compname, Sqw->lookup_length*sizeof(long));
  } else {
    int i;
    for (i=0; i < Sqw->lookup_length; i++) {
      double w = (double)i/(double)Sqw->lookup_length; /* a random number tabulated value */
      SW_lookup[i] = Sqw_search_SW(*Sqw_Data, w);
    }
    Sqw_Data->SW_lookup = SW_lookup;
  }
  QW_lookup = (long**)calloc(w_bins, sizeof(long*));

  if (!QW_lookup) {
    printf("Isotropic_Sqw: %s: Cannot allocate QW_lookup (%li bytes).\n"
           "Warning        Will be slower.\n",
      Sqw->compname, w_bins*sizeof(long*));
  } else {
    for (index_w=0; index_w < w_bins ; index_w++) {
      QW_lookup[index_w] =
        (long*)calloc(Sqw->lookup_length, sizeof(long));
      if (!QW_lookup[index_w]) {
        printf("Isotropic_Sqw: %s: Cannot allocate QW_lookup[%li] (%li bytes).\n"
               "Warning        Will be slower.\n",
        Sqw->compname, index_w, Sqw->lookup_length*sizeof(long));
        free(QW_lookup); Sqw_Data->QW_lookup = QW_lookup = NULL; break;
      } else {
        int i;
        for (i=0; i < Sqw->lookup_length; i++) {
          double w = (double)i/(double)Sqw->lookup_length; /* a random number tabulated value */
          QW_lookup[index_w][i] = Sqw_search_Q_proba_per_w(*Sqw_Data, w, index_w);
        }
      }
    }
    Sqw_Data->QW_lookup = QW_lookup;
  }
  if ((Sqw_Data->QW_lookup || Sqw_Data->SW_lookup) && Sqw->verbose_output > 2) {
    MPI_MASTER(
    printf("Isotropic_Sqw: %s: Generated lookup tables with %li entries\n",
      Sqw->compname, Sqw->lookup_length);
    );
  }

  return(Sqw_Data);
} /* end Sqw_readfile */

/*****************************************************************************
* Sqw_init_read: Read coherent/incoherent Sqw data files
*   Returns Sqw total intensity, or 0 (error)
* Used in : INITIALIZE (1)
*****************************************************************************/
double Sqw_init(struct Sqw_sample_struct *Sqw, char *file_coh, char *file_inc)
{
  double ret=0;

  /* read files */
  struct Sqw_Data_struct *d_inc, *d_coh;
  Sqw->type = 'i';
  d_inc = Sqw_readfile(Sqw, file_inc, &(Sqw->Data_inc));
  Sqw->type = 'c';
  d_coh = Sqw_readfile(Sqw, file_coh, &(Sqw->Data_coh));

  if (d_inc && !d_inc->intensity && Sqw->s_inc>0) {
    MPI_MASTER(
    if (Sqw->verbose_output > 0)
      printf("Isotropic_Sqw: %s: Using Isotropic elastic incoherent scattering (sigma=%g [barns])\n", Sqw->compname, Sqw->s_inc);
    );
    ret=1;
  }

  if (!d_inc || !d_coh) return(0);

  d_coh->type = 'c';
  Sqw->Data_inc.type = 'i';
  MPI_MASTER(
  if (d_coh && !d_coh->intensity && Sqw->s_coh)
    printf("Isotropic_Sqw: %s: Coherent scattering Sqw intensity is null.\n"
           "Warning        Disabling coherent scattering.\n", Sqw->compname);
  );
  if (d_inc && d_coh && d_inc->intensity && d_coh->intensity) {
    char msg[80];
    strcpy(msg, "");
    /* check dimensions/limits for Q, Energy in coh and inc Tables */
    if (d_inc->q_bins  != d_coh->q_bins)
      strcpy(msg, "Q axis size");
    if (d_inc->w_bins  != d_coh->w_bins)
      strcpy(msg, "Energy axis size");
    if (d_inc->q_max != d_coh->q_max)
      strcpy(msg, "Q axis limits");
    if (d_inc->w_max != d_coh->w_max)
      strcpy(msg, "Energy axis limits");
    MPI_MASTER(
    if (strlen(msg)) {
      printf("Isotropic_Sqw: %s: Sqw data from files %s and %s do not match\n"
             "WARNING        wrong %s\n",
             Sqw->compname, file_coh, file_inc, msg);
    }
    );
  }

  if (!ret) ret=d_inc->intensity+d_coh->intensity;
  return(ret);
} /* Sqw_init */

#endif /* definied ISOTROPIC_SQW */
#line 11258 "Samples_Incoherent.c"

/* Shared user declarations for all components 'Incoherent'. */
#line 128 "/usr/share/mcstas/2.6rc1/samples/Incoherent.comp"



struct StructVarsInc
{
double  sigma_a; /* Absorption cross section per atom (barns) */
    double  sigma_i; /* Incoherent scattering cross section per atom (barns) */
    double  rho;     /* Density of atoms (AA-3) */
    double  my_s;
    double  my_a_v;
    int     shape;       /* 0 cylinder, 1 box, 2 sphere, 3 OFF file */
    double  aw,ah;       /* rectangular angular dimensions */
    double  xw,yh;       /* rectangular metrical dimensions */
    double  tx,ty,tz;    /* target coords */
    };
#line 11277 "Samples_Incoherent.c"

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

#line 11363 "Samples_Incoherent.c"

/* Shared user declarations for all components 'Monitor_nD'. */
#line 216 "/usr/share/mcstas/2.6rc1/monitors/Monitor_nD.comp"
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



#line 13307 "Samples_Incoherent.c"

/* Instrument parameters. */
MCNUM mcipL_min;
MCNUM mcipL_max;
MCNUM mcipSAMPLE;
MCNUM mcipSTOP;
MCNUM mcipV_multiples;
MCNUM mcipOrder;
MCNUM mcipINC;
MCNUM mcipABS;
MCNUM mcipDB;
char* mcipISISFACE;

#define mcNUMIPAR 10
int mcnumipar = 10;
struct mcinputtable_struct mcinputtable[mcNUMIPAR+1] = {
  "L_min", &mcipL_min, instr_type_double, "0.5", 
  "L_max", &mcipL_max, instr_type_double, "7", 
  "SAMPLE", &mcipSAMPLE, instr_type_double, "", 
  "STOP", &mcipSTOP, instr_type_double, "0", 
  "V_multiples", &mcipV_multiples, instr_type_double, "1", 
  "Order", &mcipOrder, instr_type_double, "0", 
  "INC", &mcipINC, instr_type_double, "5.08", 
  "ABS", &mcipABS, instr_type_double, "5.08", 
  "DB", &mcipDB, instr_type_double, "0", 
  "ISISFACE", &mcipISISFACE, instr_type_string, "hydrogen", 
  NULL, NULL, instr_type_double, ""
};

/* User declarations from instrument definition. */
#define mccompcurname  Samples_Incoherent
#define mccompcurtype  INSTRUMENT
#define mccompcurindex 0
#define mcposaSamples_Incoherent coords_set(0,0,0)
#define L_min mcipL_min
#define L_max mcipL_max
#define SAMPLE mcipSAMPLE
#define STOP mcipSTOP
#define V_multiples mcipV_multiples
#define Order mcipOrder
#define INC mcipINC
#define ABS mcipABS
#define DB mcipDB
#define ISISFACE mcipISISFACE
#line 85 "Samples_Incoherent.instr"
  double e_min,e_max;
#line 13354 "Samples_Incoherent.c"
#undef ISISFACE
#undef DB
#undef ABS
#undef INC
#undef Order
#undef V_multiples
#undef STOP
#undef SAMPLE
#undef L_max
#undef L_min
#undef mcposaSamples_Incoherent
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

/* Definition parameters for component 'Source' [1]. */
#define mccSource_Face mcipISISFACE /* declared as a string. May produce warnings at compile */
/* Setting parameters for component 'Source' [1]. */
MCNUM mccSource_Emin;
MCNUM mccSource_Emax;
MCNUM mccSource_dist;
MCNUM mccSource_focus_xw;
MCNUM mccSource_focus_yh;
MCNUM mccSource_xwidth;
MCNUM mccSource_yheight;
MCNUM mccSource_CAngle;
MCNUM mccSource_SAC;
MCNUM mccSource_Lmin;
MCNUM mccSource_Lmax;
int mccSource_target_index;
MCNUM mccSource_verbose;

/* Setting parameters for component 'Coord' [2]. */
char mccCoord_profile[16384];
MCNUM mccCoord_percent;
MCNUM mccCoord_flag_save;
MCNUM mccCoord_minutes;

/* Setting parameters for component 'Vsample' [3]. */
MCNUM mccVsample_radius;
MCNUM mccVsample_thickness;
MCNUM mccVsample_zdepth;
MCNUM mccVsample_Vc;
MCNUM mccVsample_sigma_abs;
MCNUM mccVsample_sigma_inc;
MCNUM mccVsample_radius_i;
MCNUM mccVsample_radius_o;
MCNUM mccVsample_h;
MCNUM mccVsample_focus_r;
MCNUM mccVsample_pack;
MCNUM mccVsample_frac;
MCNUM mccVsample_f_QE;
MCNUM mccVsample_gamma;
MCNUM mccVsample_target_x;
MCNUM mccVsample_target_y;
MCNUM mccVsample_target_z;
MCNUM mccVsample_focus_xw;
MCNUM mccVsample_focus_yh;
MCNUM mccVsample_focus_aw;
MCNUM mccVsample_focus_ah;
MCNUM mccVsample_xwidth;
MCNUM mccVsample_yheight;
MCNUM mccVsample_zthick;
MCNUM mccVsample_rad_sphere;
MCNUM mccVsample_sig_a;
MCNUM mccVsample_sig_i;
MCNUM mccVsample_V0;
int mccVsample_target_index;
MCNUM mccVsample_multiples;

/* Definition parameters for component 'PowN' [4]. */
#define mccPowN_format Undefined
/* Setting parameters for component 'PowN' [4]. */
char mccPowN_reflections[16384];
char mccPowN_geometry[16384];
MCNUM mccPowN_radius;
MCNUM mccPowN_yheight;
MCNUM mccPowN_xwidth;
MCNUM mccPowN_zdepth;
MCNUM mccPowN_thickness;
MCNUM mccPowN_pack;
MCNUM mccPowN_Vc;
MCNUM mccPowN_sigma_abs;
MCNUM mccPowN_sigma_inc;
MCNUM mccPowN_delta_d_d;
MCNUM mccPowN_p_inc;
MCNUM mccPowN_p_transmit;
MCNUM mccPowN_DW;
MCNUM mccPowN_nb_atoms;
MCNUM mccPowN_d_omega;
MCNUM mccPowN_d_phi;
MCNUM mccPowN_tth_sign;
MCNUM mccPowN_p_interact;
MCNUM mccPowN_concentric;
MCNUM mccPowN_density;
MCNUM mccPowN_weight;
MCNUM mccPowN_barns;
MCNUM mccPowN_Strain;
MCNUM mccPowN_focus_flip;
int mccPowN_target_index;

/* Definition parameters for component 'SX' [5]. */
#define mccSX_mosaic_AB Mosaic_AB_Undefined
/* Setting parameters for component 'SX' [5]. */
char mccSX_reflections[16384];
char mccSX_geometry[16384];
MCNUM mccSX_xwidth;
MCNUM mccSX_yheight;
MCNUM mccSX_zdepth;
MCNUM mccSX_radius;
MCNUM mccSX_delta_d_d;
MCNUM mccSX_mosaic;
MCNUM mccSX_mosaic_a;
MCNUM mccSX_mosaic_b;
MCNUM mccSX_mosaic_c;
MCNUM mccSX_recip_cell;
MCNUM mccSX_barns;
MCNUM mccSX_ax;
MCNUM mccSX_ay;
MCNUM mccSX_az;
MCNUM mccSX_bx;
MCNUM mccSX_by;
MCNUM mccSX_bz;
MCNUM mccSX_cx;
MCNUM mccSX_cy;
MCNUM mccSX_cz;
MCNUM mccSX_p_transmit;
MCNUM mccSX_sigma_abs;
MCNUM mccSX_sigma_inc;
MCNUM mccSX_aa;
MCNUM mccSX_bb;
MCNUM mccSX_cc;
MCNUM mccSX_order;
MCNUM mccSX_RX;
MCNUM mccSX_RY;
MCNUM mccSX_powder;
MCNUM mccSX_PG;
MCNUM mccSX_deltak;

/* Definition parameters for component 'Sqw' [6]. */
#define mccSqw_powder_format Undefined
/* Setting parameters for component 'Sqw' [6]. */
char mccSqw_Sqw_coh[16384];
char mccSqw_Sqw_inc[16384];
char mccSqw_geometry[16384];
MCNUM mccSqw_radius;
MCNUM mccSqw_thickness;
MCNUM mccSqw_xwidth;
MCNUM mccSqw_yheight;
MCNUM mccSqw_zdepth;
MCNUM mccSqw_threshold;
int mccSqw_order;
MCNUM mccSqw_T;
MCNUM mccSqw_verbose;
MCNUM mccSqw_d_phi;
int mccSqw_concentric;
MCNUM mccSqw_rho;
MCNUM mccSqw_sigma_abs;
MCNUM mccSqw_sigma_coh;
MCNUM mccSqw_sigma_inc;
MCNUM mccSqw_classical;
MCNUM mccSqw_powder_Dd;
MCNUM mccSqw_powder_DW;
MCNUM mccSqw_powder_Vc;
MCNUM mccSqw_density;
MCNUM mccSqw_weight;
MCNUM mccSqw_p_interact;
MCNUM mccSqw_norm;
MCNUM mccSqw_powder_barns;
char mccSqw_quantum_correction[16384];

/* Setting parameters for component 'IncS' [7]. */
char mccIncS_geometry[16384];
MCNUM mccIncS_radius;
MCNUM mccIncS_xwidth;
MCNUM mccIncS_yheight;
MCNUM mccIncS_zdepth;
MCNUM mccIncS_thickness;
MCNUM mccIncS_target_x;
MCNUM mccIncS_target_y;
MCNUM mccIncS_target_z;
MCNUM mccIncS_focus_r;
MCNUM mccIncS_focus_xw;
MCNUM mccIncS_focus_yh;
MCNUM mccIncS_focus_aw;
MCNUM mccIncS_focus_ah;
int mccIncS_target_index;
MCNUM mccIncS_pack;
MCNUM mccIncS_p_interact;
MCNUM mccIncS_f_QE;
MCNUM mccIncS_gamma;
MCNUM mccIncS_Etrans;
MCNUM mccIncS_deltaE;
MCNUM mccIncS_sigma_abs;
MCNUM mccIncS_sigma_inc;
MCNUM mccIncS_Vc;
MCNUM mccIncS_concentric;
MCNUM mccIncS_order;

/* Definition parameters for component 'PowNg' [8]. */
#define mccPowNg_format Undefined
/* Setting parameters for component 'PowNg' [8]. */
char mccPowNg_reflections[16384];
char mccPowNg_geometry[16384];
MCNUM mccPowNg_radius;
MCNUM mccPowNg_yheight;
MCNUM mccPowNg_xwidth;
MCNUM mccPowNg_zdepth;
MCNUM mccPowNg_thickness;
MCNUM mccPowNg_pack;
MCNUM mccPowNg_Vc;
MCNUM mccPowNg_sigma_abs;
MCNUM mccPowNg_sigma_inc;
MCNUM mccPowNg_delta_d_d;
MCNUM mccPowNg_p_inc;
MCNUM mccPowNg_p_transmit;
MCNUM mccPowNg_DW;
MCNUM mccPowNg_nb_atoms;
MCNUM mccPowNg_d_omega;
MCNUM mccPowNg_d_phi;
MCNUM mccPowNg_tth_sign;
MCNUM mccPowNg_p_interact;
MCNUM mccPowNg_concentric;
MCNUM mccPowNg_density;
MCNUM mccPowNg_weight;
MCNUM mccPowNg_barns;
MCNUM mccPowNg_Strain;
MCNUM mccPowNg_focus_flip;
int mccPowNg_target_index;

/* Definition parameters for component 'SXg' [9]. */
#define mccSXg_mosaic_AB Mosaic_AB_Undefined
/* Setting parameters for component 'SXg' [9]. */
char mccSXg_reflections[16384];
char mccSXg_geometry[16384];
MCNUM mccSXg_xwidth;
MCNUM mccSXg_yheight;
MCNUM mccSXg_zdepth;
MCNUM mccSXg_radius;
MCNUM mccSXg_delta_d_d;
MCNUM mccSXg_mosaic;
MCNUM mccSXg_mosaic_a;
MCNUM mccSXg_mosaic_b;
MCNUM mccSXg_mosaic_c;
MCNUM mccSXg_recip_cell;
MCNUM mccSXg_barns;
MCNUM mccSXg_ax;
MCNUM mccSXg_ay;
MCNUM mccSXg_az;
MCNUM mccSXg_bx;
MCNUM mccSXg_by;
MCNUM mccSXg_bz;
MCNUM mccSXg_cx;
MCNUM mccSXg_cy;
MCNUM mccSXg_cz;
MCNUM mccSXg_p_transmit;
MCNUM mccSXg_sigma_abs;
MCNUM mccSXg_sigma_inc;
MCNUM mccSXg_aa;
MCNUM mccSXg_bb;
MCNUM mccSXg_cc;
MCNUM mccSXg_order;
MCNUM mccSXg_RX;
MCNUM mccSXg_RY;
MCNUM mccSXg_powder;
MCNUM mccSXg_PG;
MCNUM mccSXg_deltak;

/* Definition parameters for component 'Sqwg' [10]. */
#define mccSqwg_powder_format Undefined
/* Setting parameters for component 'Sqwg' [10]. */
char mccSqwg_Sqw_coh[16384];
char mccSqwg_Sqw_inc[16384];
char mccSqwg_geometry[16384];
MCNUM mccSqwg_radius;
MCNUM mccSqwg_thickness;
MCNUM mccSqwg_xwidth;
MCNUM mccSqwg_yheight;
MCNUM mccSqwg_zdepth;
MCNUM mccSqwg_threshold;
int mccSqwg_order;
MCNUM mccSqwg_T;
MCNUM mccSqwg_verbose;
MCNUM mccSqwg_d_phi;
int mccSqwg_concentric;
MCNUM mccSqwg_rho;
MCNUM mccSqwg_sigma_abs;
MCNUM mccSqwg_sigma_coh;
MCNUM mccSqwg_sigma_inc;
MCNUM mccSqwg_classical;
MCNUM mccSqwg_powder_Dd;
MCNUM mccSqwg_powder_DW;
MCNUM mccSqwg_powder_Vc;
MCNUM mccSqwg_density;
MCNUM mccSqwg_weight;
MCNUM mccSqwg_p_interact;
MCNUM mccSqwg_norm;
MCNUM mccSqwg_powder_barns;
char mccSqwg_quantum_correction[16384];

/* Setting parameters for component 'IncSg' [11]. */
char mccIncSg_geometry[16384];
MCNUM mccIncSg_radius;
MCNUM mccIncSg_xwidth;
MCNUM mccIncSg_yheight;
MCNUM mccIncSg_zdepth;
MCNUM mccIncSg_thickness;
MCNUM mccIncSg_target_x;
MCNUM mccIncSg_target_y;
MCNUM mccIncSg_target_z;
MCNUM mccIncSg_focus_r;
MCNUM mccIncSg_focus_xw;
MCNUM mccIncSg_focus_yh;
MCNUM mccIncSg_focus_aw;
MCNUM mccIncSg_focus_ah;
int mccIncSg_target_index;
MCNUM mccIncSg_pack;
MCNUM mccIncSg_p_interact;
MCNUM mccIncSg_f_QE;
MCNUM mccIncSg_gamma;
MCNUM mccIncSg_Etrans;
MCNUM mccIncSg_deltaE;
MCNUM mccIncSg_sigma_abs;
MCNUM mccIncSg_sigma_inc;
MCNUM mccIncSg_Vc;
MCNUM mccIncSg_concentric;
MCNUM mccIncSg_order;

/* Definition parameters for component 'SqwV' [12]. */
#define mccSqwV_powder_format Undefined
/* Setting parameters for component 'SqwV' [12]. */
char mccSqwV_Sqw_coh[16384];
char mccSqwV_Sqw_inc[16384];
char mccSqwV_geometry[16384];
MCNUM mccSqwV_radius;
MCNUM mccSqwV_thickness;
MCNUM mccSqwV_xwidth;
MCNUM mccSqwV_yheight;
MCNUM mccSqwV_zdepth;
MCNUM mccSqwV_threshold;
int mccSqwV_order;
MCNUM mccSqwV_T;
MCNUM mccSqwV_verbose;
MCNUM mccSqwV_d_phi;
int mccSqwV_concentric;
MCNUM mccSqwV_rho;
MCNUM mccSqwV_sigma_abs;
MCNUM mccSqwV_sigma_coh;
MCNUM mccSqwV_sigma_inc;
MCNUM mccSqwV_classical;
MCNUM mccSqwV_powder_Dd;
MCNUM mccSqwV_powder_DW;
MCNUM mccSqwV_powder_Vc;
MCNUM mccSqwV_density;
MCNUM mccSqwV_weight;
MCNUM mccSqwV_p_interact;
MCNUM mccSqwV_norm;
MCNUM mccSqwV_powder_barns;
char mccSqwV_quantum_correction[16384];

/* Setting parameters for component 'Dirbeam' [13]. */
int mccDirbeam_nx;
int mccDirbeam_ny;
char mccDirbeam_filename[16384];
MCNUM mccDirbeam_xmin;
MCNUM mccDirbeam_xmax;
MCNUM mccDirbeam_ymin;
MCNUM mccDirbeam_ymax;
MCNUM mccDirbeam_xwidth;
MCNUM mccDirbeam_yheight;
MCNUM mccDirbeam_restore_neutron;

/* Setting parameters for component 'Stop' [14]. */
MCNUM mccStop_xmin;
MCNUM mccStop_xmax;
MCNUM mccStop_ymin;
MCNUM mccStop_ymax;
MCNUM mccStop_xwidth;
MCNUM mccStop_yheight;
MCNUM mccStop_radius;

/* Definition parameters for component 'Tof_Sphere_4pi' [15]. */
#define mccTof_Sphere_4pi_user1 FLT_MAX
#define mccTof_Sphere_4pi_user2 FLT_MAX
#define mccTof_Sphere_4pi_user3 FLT_MAX
/* Setting parameters for component 'Tof_Sphere_4pi' [15]. */
MCNUM mccTof_Sphere_4pi_xwidth;
MCNUM mccTof_Sphere_4pi_yheight;
MCNUM mccTof_Sphere_4pi_zdepth;
MCNUM mccTof_Sphere_4pi_xmin;
MCNUM mccTof_Sphere_4pi_xmax;
MCNUM mccTof_Sphere_4pi_ymin;
MCNUM mccTof_Sphere_4pi_ymax;
MCNUM mccTof_Sphere_4pi_zmin;
MCNUM mccTof_Sphere_4pi_zmax;
MCNUM mccTof_Sphere_4pi_bins;
MCNUM mccTof_Sphere_4pi_min;
MCNUM mccTof_Sphere_4pi_max;
MCNUM mccTof_Sphere_4pi_restore_neutron;
MCNUM mccTof_Sphere_4pi_radius;
char mccTof_Sphere_4pi_options[16384];
char mccTof_Sphere_4pi_filename[16384];
char mccTof_Sphere_4pi_geometry[16384];
char mccTof_Sphere_4pi_username1[16384];
char mccTof_Sphere_4pi_username2[16384];
char mccTof_Sphere_4pi_username3[16384];
int mccTof_Sphere_4pi_nowritefile;

/* Definition parameters for component 'PSD_Sphere_4pi' [16]. */
#define mccPSD_Sphere_4pi_nx 180
#define mccPSD_Sphere_4pi_ny 180
/* Setting parameters for component 'PSD_Sphere_4pi' [16]. */
char mccPSD_Sphere_4pi_filename[16384];
MCNUM mccPSD_Sphere_4pi_radius;
MCNUM mccPSD_Sphere_4pi_restore_neutron;
int mccPSD_Sphere_4pi_nowritefile;

/* Definition parameters for component 'Sphere_4pi' [18]. */
#define mccSphere_4pi_nx 180
#define mccSphere_4pi_ny 180
/* Setting parameters for component 'Sphere_4pi' [18]. */
char mccSphere_4pi_filename[16384];
MCNUM mccSphere_4pi_radius;
MCNUM mccSphere_4pi_restore_neutron;
int mccSphere_4pi_nowritefile;

/* User component declarations. */

/* User declarations for component 'Source' [1]. */
#define mccompcurname  Source
#define mccompcurtype  ISIS_moderator
#define mccompcurindex 1
#define Face mccSource_Face
#define p_in mccSource_p_in
#define Tnpts mccSource_Tnpts
#define scaleSize mccSource_scaleSize
#define angleArea mccSource_angleArea
#define Nsim mccSource_Nsim
#define Ncount mccSource_Ncount
#define TS mccSource_TS
#define rtE0 mccSource_rtE0
#define rtE1 mccSource_rtE1
#define rtmodX mccSource_rtmodX
#define rtmodY mccSource_rtmodY
#define TargetStation mccSource_TargetStation
#define CurrentWeight mccSource_CurrentWeight
#define Emin mccSource_Emin
#define Emax mccSource_Emax
#define dist mccSource_dist
#define focus_xw mccSource_focus_xw
#define focus_yh mccSource_focus_yh
#define xwidth mccSource_xwidth
#define yheight mccSource_yheight
#define CAngle mccSource_CAngle
#define SAC mccSource_SAC
#define Lmin mccSource_Lmin
#define Lmax mccSource_Lmax
#define target_index mccSource_target_index
#define verbose mccSource_verbose
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
#line 13852 "Samples_Incoherent.c"
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

/* User declarations for component 'Coord' [2]. */
#define mccompcurname  Coord
#define mccompcurtype  Progress_bar
#define mccompcurindex 2
#define IntermediateCnts mccCoord_IntermediateCnts
#define StartTime mccCoord_StartTime
#define EndTime mccCoord_EndTime
#define CurrentTime mccCoord_CurrentTime
#define profile mccCoord_profile
#define percent mccCoord_percent
#define flag_save mccCoord_flag_save
#define minutes mccCoord_minutes
#line 44 "/usr/share/mcstas/2.6rc1/misc/Progress_bar.comp"
#ifndef PROGRESS_BAR
#define PROGRESS_BAR
#else
#error Only one Progress_bar component may be used in an instrument definition.
#endif

double IntermediateCnts;
time_t StartTime;
time_t EndTime;
time_t CurrentTime;
#line 13907 "Samples_Incoherent.c"
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

/* User declarations for component 'Vsample' [3]. */
#define mccompcurname  Vsample
#define mccompcurtype  V_sample
#define mccompcurindex 3
#define VarsV mccVsample_VarsV
#define radius mccVsample_radius
#define thickness mccVsample_thickness
#define zdepth mccVsample_zdepth
#define Vc mccVsample_Vc
#define sigma_abs mccVsample_sigma_abs
#define sigma_inc mccVsample_sigma_inc
#define radius_i mccVsample_radius_i
#define radius_o mccVsample_radius_o
#define h mccVsample_h
#define focus_r mccVsample_focus_r
#define pack mccVsample_pack
#define frac mccVsample_frac
#define f_QE mccVsample_f_QE
#define gamma mccVsample_gamma
#define target_x mccVsample_target_x
#define target_y mccVsample_target_y
#define target_z mccVsample_target_z
#define focus_xw mccVsample_focus_xw
#define focus_yh mccVsample_focus_yh
#define focus_aw mccVsample_focus_aw
#define focus_ah mccVsample_focus_ah
#define xwidth mccVsample_xwidth
#define yheight mccVsample_yheight
#define zthick mccVsample_zthick
#define rad_sphere mccVsample_rad_sphere
#define sig_a mccVsample_sig_a
#define sig_i mccVsample_sig_i
#define V0 mccVsample_V0
#define target_index mccVsample_target_index
#define multiples mccVsample_multiples
#line 117 "/usr/share/mcstas/2.6rc1/obsolete/V_sample.comp"
  struct StructVarsV VarsV;
#line 13957 "Samples_Incoherent.c"
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

/* User declarations for component 'PowN' [4]. */
#define mccompcurname  PowN
#define mccompcurtype  PowderN
#define mccompcurindex 4
#define format mccPowN_format
#define line_info mccPowN_line_info
#define columns mccPowN_columns
#define offdata mccPowN_offdata
#define reflections mccPowN_reflections
#define geometry mccPowN_geometry
#define radius mccPowN_radius
#define yheight mccPowN_yheight
#define xwidth mccPowN_xwidth
#define zdepth mccPowN_zdepth
#define thickness mccPowN_thickness
#define pack mccPowN_pack
#define Vc mccPowN_Vc
#define sigma_abs mccPowN_sigma_abs
#define sigma_inc mccPowN_sigma_inc
#define delta_d_d mccPowN_delta_d_d
#define p_inc mccPowN_p_inc
#define p_transmit mccPowN_p_transmit
#define DW mccPowN_DW
#define nb_atoms mccPowN_nb_atoms
#define d_omega mccPowN_d_omega
#define d_phi mccPowN_d_phi
#define tth_sign mccPowN_tth_sign
#define p_interact mccPowN_p_interact
#define concentric mccPowN_concentric
#define density mccPowN_density
#define weight mccPowN_weight
#define barns mccPowN_barns
#define Strain mccPowN_Strain
#define focus_flip mccPowN_focus_flip
#define target_index mccPowN_target_index
#line 551 "/usr/share/mcstas/2.6rc1/samples/PowderN.comp"
  struct line_info_struct line_info;
  int *columns;
  off_struct offdata;
  double tgt_x, tgt_y, tgt_z;

#line 14034 "Samples_Incoherent.c"
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

/* User declarations for component 'SX' [5]. */
#define mccompcurname  SX
#define mccompcurtype  Single_crystal
#define mccompcurindex 5
#define mosaic_AB mccSX_mosaic_AB
#define hkl_info mccSX_hkl_info
#define offdata mccSX_offdata
#define reflections mccSX_reflections
#define geometry mccSX_geometry
#define xwidth mccSX_xwidth
#define yheight mccSX_yheight
#define zdepth mccSX_zdepth
#define radius mccSX_radius
#define delta_d_d mccSX_delta_d_d
#define mosaic mccSX_mosaic
#define mosaic_a mccSX_mosaic_a
#define mosaic_b mccSX_mosaic_b
#define mosaic_c mccSX_mosaic_c
#define recip_cell mccSX_recip_cell
#define barns mccSX_barns
#define ax mccSX_ax
#define ay mccSX_ay
#define az mccSX_az
#define bx mccSX_bx
#define by mccSX_by
#define bz mccSX_bz
#define cx mccSX_cx
#define cy mccSX_cy
#define cz mccSX_cz
#define p_transmit mccSX_p_transmit
#define sigma_abs mccSX_sigma_abs
#define sigma_inc mccSX_sigma_inc
#define aa mccSX_aa
#define bb mccSX_bb
#define cc mccSX_cc
#define order mccSX_order
#define RX mccSX_RX
#define RY mccSX_RY
#define powder mccSX_powder
#define PG mccSX_PG
#define deltak mccSX_deltak
#line 971 "/usr/share/mcstas/2.6rc1/samples/Single_crystal.comp"
  struct hkl_info_struct hkl_info;
  off_struct             offdata;
#line 14114 "Samples_Incoherent.c"
#undef deltak
#undef PG
#undef powder
#undef RY
#undef RX
#undef order
#undef cc
#undef bb
#undef aa
#undef sigma_inc
#undef sigma_abs
#undef p_transmit
#undef cz
#undef cy
#undef cx
#undef bz
#undef by
#undef bx
#undef az
#undef ay
#undef ax
#undef barns
#undef recip_cell
#undef mosaic_c
#undef mosaic_b
#undef mosaic_a
#undef mosaic
#undef delta_d_d
#undef radius
#undef zdepth
#undef yheight
#undef xwidth
#undef geometry
#undef reflections
#undef offdata
#undef hkl_info
#undef mosaic_AB
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'Sqw' [6]. */
#define mccompcurname  Sqw
#define mccompcurtype  Isotropic_Sqw
#define mccompcurindex 6
#define powder_format mccSqw_powder_format
#define VarSqw mccSqw_VarSqw
#define columns mccSqw_columns
#define offdata mccSqw_offdata
#define Sqw_coh mccSqw_Sqw_coh
#define Sqw_inc mccSqw_Sqw_inc
#define geometry mccSqw_geometry
#define radius mccSqw_radius
#define thickness mccSqw_thickness
#define xwidth mccSqw_xwidth
#define yheight mccSqw_yheight
#define zdepth mccSqw_zdepth
#define threshold mccSqw_threshold
#define order mccSqw_order
#define T mccSqw_T
#define verbose mccSqw_verbose
#define d_phi mccSqw_d_phi
#define concentric mccSqw_concentric
#define rho mccSqw_rho
#define sigma_abs mccSqw_sigma_abs
#define sigma_coh mccSqw_sigma_coh
#define sigma_inc mccSqw_sigma_inc
#define classical mccSqw_classical
#define powder_Dd mccSqw_powder_Dd
#define powder_DW mccSqw_powder_DW
#define powder_Vc mccSqw_powder_Vc
#define density mccSqw_density
#define weight mccSqw_weight
#define p_interact mccSqw_p_interact
#define norm mccSqw_norm
#define powder_barns mccSqw_powder_barns
#define quantum_correction mccSqw_quantum_correction
#line 1999 "/usr/share/mcstas/2.6rc1/samples/Isotropic_Sqw.comp"
  struct Sqw_sample_struct VarSqw;
  int *columns;
  off_struct offdata;
#line 14196 "Samples_Incoherent.c"
#undef quantum_correction
#undef powder_barns
#undef norm
#undef p_interact
#undef weight
#undef density
#undef powder_Vc
#undef powder_DW
#undef powder_Dd
#undef classical
#undef sigma_inc
#undef sigma_coh
#undef sigma_abs
#undef rho
#undef concentric
#undef d_phi
#undef verbose
#undef T
#undef order
#undef threshold
#undef zdepth
#undef yheight
#undef xwidth
#undef thickness
#undef radius
#undef geometry
#undef Sqw_inc
#undef Sqw_coh
#undef offdata
#undef columns
#undef VarSqw
#undef powder_format
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'IncS' [7]. */
#define mccompcurname  IncS
#define mccompcurtype  Incoherent
#define mccompcurindex 7
#define VarsInc mccIncS_VarsInc
#define offdata mccIncS_offdata
#define geometry mccIncS_geometry
#define radius mccIncS_radius
#define xwidth mccIncS_xwidth
#define yheight mccIncS_yheight
#define zdepth mccIncS_zdepth
#define thickness mccIncS_thickness
#define target_x mccIncS_target_x
#define target_y mccIncS_target_y
#define target_z mccIncS_target_z
#define focus_r mccIncS_focus_r
#define focus_xw mccIncS_focus_xw
#define focus_yh mccIncS_focus_yh
#define focus_aw mccIncS_focus_aw
#define focus_ah mccIncS_focus_ah
#define target_index mccIncS_target_index
#define pack mccIncS_pack
#define p_interact mccIncS_p_interact
#define f_QE mccIncS_f_QE
#define gamma mccIncS_gamma
#define Etrans mccIncS_Etrans
#define deltaE mccIncS_deltaE
#define sigma_abs mccIncS_sigma_abs
#define sigma_inc mccIncS_sigma_inc
#define Vc mccIncS_Vc
#define concentric mccIncS_concentric
#define order mccIncS_order
#line 149 "/usr/share/mcstas/2.6rc1/samples/Incoherent.comp"

  struct StructVarsInc VarsInc;
  off_struct offdata;

#line 14270 "Samples_Incoherent.c"
#undef order
#undef concentric
#undef Vc
#undef sigma_inc
#undef sigma_abs
#undef deltaE
#undef Etrans
#undef gamma
#undef f_QE
#undef p_interact
#undef pack
#undef target_index
#undef focus_ah
#undef focus_aw
#undef focus_yh
#undef focus_xw
#undef focus_r
#undef target_z
#undef target_y
#undef target_x
#undef thickness
#undef zdepth
#undef yheight
#undef xwidth
#undef radius
#undef geometry
#undef offdata
#undef VarsInc
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'PowNg' [8]. */
#define mccompcurname  PowNg
#define mccompcurtype  PowderN
#define mccompcurindex 8
#define format mccPowNg_format
#define line_info mccPowNg_line_info
#define columns mccPowNg_columns
#define offdata mccPowNg_offdata
#define reflections mccPowNg_reflections
#define geometry mccPowNg_geometry
#define radius mccPowNg_radius
#define yheight mccPowNg_yheight
#define xwidth mccPowNg_xwidth
#define zdepth mccPowNg_zdepth
#define thickness mccPowNg_thickness
#define pack mccPowNg_pack
#define Vc mccPowNg_Vc
#define sigma_abs mccPowNg_sigma_abs
#define sigma_inc mccPowNg_sigma_inc
#define delta_d_d mccPowNg_delta_d_d
#define p_inc mccPowNg_p_inc
#define p_transmit mccPowNg_p_transmit
#define DW mccPowNg_DW
#define nb_atoms mccPowNg_nb_atoms
#define d_omega mccPowNg_d_omega
#define d_phi mccPowNg_d_phi
#define tth_sign mccPowNg_tth_sign
#define p_interact mccPowNg_p_interact
#define concentric mccPowNg_concentric
#define density mccPowNg_density
#define weight mccPowNg_weight
#define barns mccPowNg_barns
#define Strain mccPowNg_Strain
#define focus_flip mccPowNg_focus_flip
#define target_index mccPowNg_target_index
#line 551 "/usr/share/mcstas/2.6rc1/samples/PowderN.comp"
  struct line_info_struct line_info;
  int *columns;
  off_struct offdata;
  double tgt_x, tgt_y, tgt_z;

#line 14344 "Samples_Incoherent.c"
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

/* User declarations for component 'SXg' [9]. */
#define mccompcurname  SXg
#define mccompcurtype  Single_crystal
#define mccompcurindex 9
#define mosaic_AB mccSXg_mosaic_AB
#define hkl_info mccSXg_hkl_info
#define offdata mccSXg_offdata
#define reflections mccSXg_reflections
#define geometry mccSXg_geometry
#define xwidth mccSXg_xwidth
#define yheight mccSXg_yheight
#define zdepth mccSXg_zdepth
#define radius mccSXg_radius
#define delta_d_d mccSXg_delta_d_d
#define mosaic mccSXg_mosaic
#define mosaic_a mccSXg_mosaic_a
#define mosaic_b mccSXg_mosaic_b
#define mosaic_c mccSXg_mosaic_c
#define recip_cell mccSXg_recip_cell
#define barns mccSXg_barns
#define ax mccSXg_ax
#define ay mccSXg_ay
#define az mccSXg_az
#define bx mccSXg_bx
#define by mccSXg_by
#define bz mccSXg_bz
#define cx mccSXg_cx
#define cy mccSXg_cy
#define cz mccSXg_cz
#define p_transmit mccSXg_p_transmit
#define sigma_abs mccSXg_sigma_abs
#define sigma_inc mccSXg_sigma_inc
#define aa mccSXg_aa
#define bb mccSXg_bb
#define cc mccSXg_cc
#define order mccSXg_order
#define RX mccSXg_RX
#define RY mccSXg_RY
#define powder mccSXg_powder
#define PG mccSXg_PG
#define deltak mccSXg_deltak
#line 971 "/usr/share/mcstas/2.6rc1/samples/Single_crystal.comp"
  struct hkl_info_struct hkl_info;
  off_struct             offdata;
#line 14424 "Samples_Incoherent.c"
#undef deltak
#undef PG
#undef powder
#undef RY
#undef RX
#undef order
#undef cc
#undef bb
#undef aa
#undef sigma_inc
#undef sigma_abs
#undef p_transmit
#undef cz
#undef cy
#undef cx
#undef bz
#undef by
#undef bx
#undef az
#undef ay
#undef ax
#undef barns
#undef recip_cell
#undef mosaic_c
#undef mosaic_b
#undef mosaic_a
#undef mosaic
#undef delta_d_d
#undef radius
#undef zdepth
#undef yheight
#undef xwidth
#undef geometry
#undef reflections
#undef offdata
#undef hkl_info
#undef mosaic_AB
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'Sqwg' [10]. */
#define mccompcurname  Sqwg
#define mccompcurtype  Isotropic_Sqw
#define mccompcurindex 10
#define powder_format mccSqwg_powder_format
#define VarSqw mccSqwg_VarSqw
#define columns mccSqwg_columns
#define offdata mccSqwg_offdata
#define Sqw_coh mccSqwg_Sqw_coh
#define Sqw_inc mccSqwg_Sqw_inc
#define geometry mccSqwg_geometry
#define radius mccSqwg_radius
#define thickness mccSqwg_thickness
#define xwidth mccSqwg_xwidth
#define yheight mccSqwg_yheight
#define zdepth mccSqwg_zdepth
#define threshold mccSqwg_threshold
#define order mccSqwg_order
#define T mccSqwg_T
#define verbose mccSqwg_verbose
#define d_phi mccSqwg_d_phi
#define concentric mccSqwg_concentric
#define rho mccSqwg_rho
#define sigma_abs mccSqwg_sigma_abs
#define sigma_coh mccSqwg_sigma_coh
#define sigma_inc mccSqwg_sigma_inc
#define classical mccSqwg_classical
#define powder_Dd mccSqwg_powder_Dd
#define powder_DW mccSqwg_powder_DW
#define powder_Vc mccSqwg_powder_Vc
#define density mccSqwg_density
#define weight mccSqwg_weight
#define p_interact mccSqwg_p_interact
#define norm mccSqwg_norm
#define powder_barns mccSqwg_powder_barns
#define quantum_correction mccSqwg_quantum_correction
#line 1999 "/usr/share/mcstas/2.6rc1/samples/Isotropic_Sqw.comp"
  struct Sqw_sample_struct VarSqw;
  int *columns;
  off_struct offdata;
#line 14506 "Samples_Incoherent.c"
#undef quantum_correction
#undef powder_barns
#undef norm
#undef p_interact
#undef weight
#undef density
#undef powder_Vc
#undef powder_DW
#undef powder_Dd
#undef classical
#undef sigma_inc
#undef sigma_coh
#undef sigma_abs
#undef rho
#undef concentric
#undef d_phi
#undef verbose
#undef T
#undef order
#undef threshold
#undef zdepth
#undef yheight
#undef xwidth
#undef thickness
#undef radius
#undef geometry
#undef Sqw_inc
#undef Sqw_coh
#undef offdata
#undef columns
#undef VarSqw
#undef powder_format
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'IncSg' [11]. */
#define mccompcurname  IncSg
#define mccompcurtype  Incoherent
#define mccompcurindex 11
#define VarsInc mccIncSg_VarsInc
#define offdata mccIncSg_offdata
#define geometry mccIncSg_geometry
#define radius mccIncSg_radius
#define xwidth mccIncSg_xwidth
#define yheight mccIncSg_yheight
#define zdepth mccIncSg_zdepth
#define thickness mccIncSg_thickness
#define target_x mccIncSg_target_x
#define target_y mccIncSg_target_y
#define target_z mccIncSg_target_z
#define focus_r mccIncSg_focus_r
#define focus_xw mccIncSg_focus_xw
#define focus_yh mccIncSg_focus_yh
#define focus_aw mccIncSg_focus_aw
#define focus_ah mccIncSg_focus_ah
#define target_index mccIncSg_target_index
#define pack mccIncSg_pack
#define p_interact mccIncSg_p_interact
#define f_QE mccIncSg_f_QE
#define gamma mccIncSg_gamma
#define Etrans mccIncSg_Etrans
#define deltaE mccIncSg_deltaE
#define sigma_abs mccIncSg_sigma_abs
#define sigma_inc mccIncSg_sigma_inc
#define Vc mccIncSg_Vc
#define concentric mccIncSg_concentric
#define order mccIncSg_order
#line 149 "/usr/share/mcstas/2.6rc1/samples/Incoherent.comp"

  struct StructVarsInc VarsInc;
  off_struct offdata;

#line 14580 "Samples_Incoherent.c"
#undef order
#undef concentric
#undef Vc
#undef sigma_inc
#undef sigma_abs
#undef deltaE
#undef Etrans
#undef gamma
#undef f_QE
#undef p_interact
#undef pack
#undef target_index
#undef focus_ah
#undef focus_aw
#undef focus_yh
#undef focus_xw
#undef focus_r
#undef target_z
#undef target_y
#undef target_x
#undef thickness
#undef zdepth
#undef yheight
#undef xwidth
#undef radius
#undef geometry
#undef offdata
#undef VarsInc
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'SqwV' [12]. */
#define mccompcurname  SqwV
#define mccompcurtype  Isotropic_Sqw
#define mccompcurindex 12
#define powder_format mccSqwV_powder_format
#define VarSqw mccSqwV_VarSqw
#define columns mccSqwV_columns
#define offdata mccSqwV_offdata
#define Sqw_coh mccSqwV_Sqw_coh
#define Sqw_inc mccSqwV_Sqw_inc
#define geometry mccSqwV_geometry
#define radius mccSqwV_radius
#define thickness mccSqwV_thickness
#define xwidth mccSqwV_xwidth
#define yheight mccSqwV_yheight
#define zdepth mccSqwV_zdepth
#define threshold mccSqwV_threshold
#define order mccSqwV_order
#define T mccSqwV_T
#define verbose mccSqwV_verbose
#define d_phi mccSqwV_d_phi
#define concentric mccSqwV_concentric
#define rho mccSqwV_rho
#define sigma_abs mccSqwV_sigma_abs
#define sigma_coh mccSqwV_sigma_coh
#define sigma_inc mccSqwV_sigma_inc
#define classical mccSqwV_classical
#define powder_Dd mccSqwV_powder_Dd
#define powder_DW mccSqwV_powder_DW
#define powder_Vc mccSqwV_powder_Vc
#define density mccSqwV_density
#define weight mccSqwV_weight
#define p_interact mccSqwV_p_interact
#define norm mccSqwV_norm
#define powder_barns mccSqwV_powder_barns
#define quantum_correction mccSqwV_quantum_correction
#line 1999 "/usr/share/mcstas/2.6rc1/samples/Isotropic_Sqw.comp"
  struct Sqw_sample_struct VarSqw;
  int *columns;
  off_struct offdata;
#line 14653 "Samples_Incoherent.c"
#undef quantum_correction
#undef powder_barns
#undef norm
#undef p_interact
#undef weight
#undef density
#undef powder_Vc
#undef powder_DW
#undef powder_Dd
#undef classical
#undef sigma_inc
#undef sigma_coh
#undef sigma_abs
#undef rho
#undef concentric
#undef d_phi
#undef verbose
#undef T
#undef order
#undef threshold
#undef zdepth
#undef yheight
#undef xwidth
#undef thickness
#undef radius
#undef geometry
#undef Sqw_inc
#undef Sqw_coh
#undef offdata
#undef columns
#undef VarSqw
#undef powder_format
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'Dirbeam' [13]. */
#define mccompcurname  Dirbeam
#define mccompcurtype  PSD_monitor
#define mccompcurindex 13
#define PSD_N mccDirbeam_PSD_N
#define PSD_p mccDirbeam_PSD_p
#define PSD_p2 mccDirbeam_PSD_p2
#define nx mccDirbeam_nx
#define ny mccDirbeam_ny
#define filename mccDirbeam_filename
#define xmin mccDirbeam_xmin
#define xmax mccDirbeam_xmax
#define ymin mccDirbeam_ymin
#define ymax mccDirbeam_ymax
#define xwidth mccDirbeam_xwidth
#define yheight mccDirbeam_yheight
#define restore_neutron mccDirbeam_restore_neutron
#line 62 "/usr/share/mcstas/2.6rc1/monitors/PSD_monitor.comp"
  DArray2d PSD_N;
  DArray2d PSD_p;
  DArray2d PSD_p2;
#line 14711 "Samples_Incoherent.c"
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

/* User declarations for component 'Stop' [14]. */
#define mccompcurname  Stop
#define mccompcurtype  Beamstop
#define mccompcurindex 14
#define xmin mccStop_xmin
#define xmax mccStop_xmax
#define ymin mccStop_ymin
#define ymax mccStop_ymax
#define xwidth mccStop_xwidth
#define yheight mccStop_yheight
#define radius mccStop_radius
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

/* User declarations for component 'Tof_Sphere_4pi' [15]. */
#define mccompcurname  Tof_Sphere_4pi
#define mccompcurtype  Monitor_nD
#define mccompcurindex 15
#define user1 mccTof_Sphere_4pi_user1
#define user2 mccTof_Sphere_4pi_user2
#define user3 mccTof_Sphere_4pi_user3
#define DEFS mccTof_Sphere_4pi_DEFS
#define Vars mccTof_Sphere_4pi_Vars
#define detector mccTof_Sphere_4pi_detector
#define offdata mccTof_Sphere_4pi_offdata
#define xwidth mccTof_Sphere_4pi_xwidth
#define yheight mccTof_Sphere_4pi_yheight
#define zdepth mccTof_Sphere_4pi_zdepth
#define xmin mccTof_Sphere_4pi_xmin
#define xmax mccTof_Sphere_4pi_xmax
#define ymin mccTof_Sphere_4pi_ymin
#define ymax mccTof_Sphere_4pi_ymax
#define zmin mccTof_Sphere_4pi_zmin
#define zmax mccTof_Sphere_4pi_zmax
#define bins mccTof_Sphere_4pi_bins
#define min mccTof_Sphere_4pi_min
#define max mccTof_Sphere_4pi_max
#define restore_neutron mccTof_Sphere_4pi_restore_neutron
#define radius mccTof_Sphere_4pi_radius
#define options mccTof_Sphere_4pi_options
#define filename mccTof_Sphere_4pi_filename
#define geometry mccTof_Sphere_4pi_geometry
#define username1 mccTof_Sphere_4pi_username1
#define username2 mccTof_Sphere_4pi_username2
#define username3 mccTof_Sphere_4pi_username3
#define nowritefile mccTof_Sphere_4pi_nowritefile
#line 225 "/usr/share/mcstas/2.6rc1/monitors/Monitor_nD.comp"
  MonitornD_Defines_type DEFS;
  MonitornD_Variables_type Vars;
  MCDETECTOR detector;
  off_struct offdata;
#line 14788 "Samples_Incoherent.c"
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

/* User declarations for component 'PSD_Sphere_4pi' [16]. */
#define mccompcurname  PSD_Sphere_4pi
#define mccompcurtype  PSD_monitor_4PI
#define mccompcurindex 16
#define nx mccPSD_Sphere_4pi_nx
#define ny mccPSD_Sphere_4pi_ny
#define PSD_N mccPSD_Sphere_4pi_PSD_N
#define PSD_p mccPSD_Sphere_4pi_PSD_p
#define PSD_p2 mccPSD_Sphere_4pi_PSD_p2
#define filename mccPSD_Sphere_4pi_filename
#define radius mccPSD_Sphere_4pi_radius
#define restore_neutron mccPSD_Sphere_4pi_restore_neutron
#define nowritefile mccPSD_Sphere_4pi_nowritefile
#line 55 "/usr/share/mcstas/2.6rc1/monitors/PSD_monitor_4PI.comp"
double PSD_N[nx][ny];
double PSD_p[nx][ny];
double PSD_p2[nx][ny];
#line 14838 "Samples_Incoherent.c"
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

/* User declarations for component 'Coord3' [17]. */
#define mccompcurname  Coord3
#define mccompcurtype  Arm
#define mccompcurindex 17
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'Sphere_4pi' [18]. */
#define mccompcurname  Sphere_4pi
#define mccompcurtype  PSD_monitor_4PI
#define mccompcurindex 18
#define nx mccSphere_4pi_nx
#define ny mccSphere_4pi_ny
#define PSD_N mccSphere_4pi_PSD_N
#define PSD_p mccSphere_4pi_PSD_p
#define PSD_p2 mccSphere_4pi_PSD_p2
#define filename mccSphere_4pi_filename
#define radius mccSphere_4pi_radius
#define restore_neutron mccSphere_4pi_restore_neutron
#define nowritefile mccSphere_4pi_nowritefile
#line 55 "/usr/share/mcstas/2.6rc1/monitors/PSD_monitor_4PI.comp"
double PSD_N[nx][ny];
double PSD_p[nx][ny];
double PSD_p2[nx][ny];
#line 14877 "Samples_Incoherent.c"
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

Coords mcposaSource, mcposrSource;
Rotation mcrotaSource, mcrotrSource;
Coords mcposaCoord, mcposrCoord;
Rotation mcrotaCoord, mcrotrCoord;
Coords mcposaVsample, mcposrVsample;
Rotation mcrotaVsample, mcrotrVsample;
Coords mcposaPowN, mcposrPowN;
Rotation mcrotaPowN, mcrotrPowN;
Coords mcposaSX, mcposrSX;
Rotation mcrotaSX, mcrotrSX;
Coords mcposaSqw, mcposrSqw;
Rotation mcrotaSqw, mcrotrSqw;
Coords mcposaIncS, mcposrIncS;
Rotation mcrotaIncS, mcrotrIncS;
Coords mcposaPowNg, mcposrPowNg;
Rotation mcrotaPowNg, mcrotrPowNg;
Coords mcposaSXg, mcposrSXg;
Rotation mcrotaSXg, mcrotrSXg;
Coords mcposaSqwg, mcposrSqwg;
Rotation mcrotaSqwg, mcrotrSqwg;
Coords mcposaIncSg, mcposrIncSg;
Rotation mcrotaIncSg, mcrotrIncSg;
Coords mcposaSqwV, mcposrSqwV;
Rotation mcrotaSqwV, mcrotrSqwV;
Coords mcposaDirbeam, mcposrDirbeam;
Rotation mcrotaDirbeam, mcrotrDirbeam;
Coords mcposaStop, mcposrStop;
Rotation mcrotaStop, mcrotrStop;
Coords mcposaTof_Sphere_4pi, mcposrTof_Sphere_4pi;
Rotation mcrotaTof_Sphere_4pi, mcrotrTof_Sphere_4pi;
Coords mcposaPSD_Sphere_4pi, mcposrPSD_Sphere_4pi;
Rotation mcrotaPSD_Sphere_4pi, mcrotrPSD_Sphere_4pi;
Coords mcposaCoord3, mcposrCoord3;
Rotation mcrotaCoord3, mcrotrCoord3;
Coords mcposaSphere_4pi, mcposrSphere_4pi;
Rotation mcrotaSphere_4pi, mcrotrSphere_4pi;

MCNUM mcnx, mcny, mcnz, mcnvx, mcnvy, mcnvz, mcnt, mcnsx, mcnsy, mcnsz, mcnp;

/* end declare */

void mcinit(void) {
#define mccompcurname  Samples_Incoherent
#define mccompcurtype  INSTRUMENT
#define mccompcurindex 0
#define mcposaSamples_Incoherent coords_set(0,0,0)
#define L_min mcipL_min
#define L_max mcipL_max
#define SAMPLE mcipSAMPLE
#define STOP mcipSTOP
#define V_multiples mcipV_multiples
#define Order mcipOrder
#define INC mcipINC
#define ABS mcipABS
#define DB mcipDB
#define ISISFACE mcipISISFACE
#line 88 "Samples_Incoherent.instr"
{
   e_min=81.799/L_max/L_max;
   e_max=81.799/L_min/L_min;
}
#line 14952 "Samples_Incoherent.c"
#undef ISISFACE
#undef DB
#undef ABS
#undef INC
#undef Order
#undef V_multiples
#undef STOP
#undef SAMPLE
#undef L_max
#undef L_min
#undef mcposaSamples_Incoherent
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
    /* Component Source. */
  /* Setting parameters for component Source. */
  SIG_MESSAGE("Source (Init:SetPar)");
#line 96 "Samples_Incoherent.instr"
  mccSource_Emin = e_min;
#line 96 "Samples_Incoherent.instr"
  mccSource_Emax = e_max;
#line 96 "Samples_Incoherent.instr"
  mccSource_dist = 8.299;
#line 96 "Samples_Incoherent.instr"
  mccSource_focus_xw = 0.005;
#line 97 "Samples_Incoherent.instr"
  mccSource_focus_yh = 0.005;
#line 97 "Samples_Incoherent.instr"
  mccSource_xwidth = .12;
#line 97 "Samples_Incoherent.instr"
  mccSource_yheight = .115;
#line 97 "Samples_Incoherent.instr"
  mccSource_CAngle = 0;
#line 98 "Samples_Incoherent.instr"
  mccSource_SAC = 1;
#line 58 "Samples_Incoherent.instr"
  mccSource_Lmin = 0;
#line 58 "Samples_Incoherent.instr"
  mccSource_Lmax = 0;
#line 58 "Samples_Incoherent.instr"
  mccSource_target_index = + 1;
#line 58 "Samples_Incoherent.instr"
  mccSource_verbose = 0;
#line 15006 "Samples_Incoherent.c"

  SIG_MESSAGE("Source (Init:Place/Rotate)");
  rot_set_rotation(mcrotaSource,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 15013 "Samples_Incoherent.c"
  rot_copy(mcrotrSource, mcrotaSource);
  mcposaSource = coords_set(
#line 99 "Samples_Incoherent.instr"
    0,
#line 99 "Samples_Incoherent.instr"
    0,
#line 99 "Samples_Incoherent.instr"
    0);
#line 15022 "Samples_Incoherent.c"
  mctc1 = coords_neg(mcposaSource);
  mcposrSource = rot_apply(mcrotaSource, mctc1);
  mcDEBUG_COMPONENT("Source", mcposaSource, mcrotaSource)
  mccomp_posa[1] = mcposaSource;
  mccomp_posr[1] = mcposrSource;
  mcNCounter[1]  = mcPCounter[1] = mcP2Counter[1] = 0;
  mcAbsorbProp[1]= 0;
    /* Component Coord. */
  /* Setting parameters for component Coord. */
  SIG_MESSAGE("Coord (Init:SetPar)");
#line 39 "Samples_Incoherent.instr"
  if("NULL") strncpy(mccCoord_profile, "NULL" ? "NULL" : "", 16384); else mccCoord_profile[0]='\0';
#line 39 "Samples_Incoherent.instr"
  mccCoord_percent = 10;
#line 39 "Samples_Incoherent.instr"
  mccCoord_flag_save = 0;
#line 39 "Samples_Incoherent.instr"
  mccCoord_minutes = 0;
#line 15041 "Samples_Incoherent.c"

  SIG_MESSAGE("Coord (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 15048 "Samples_Incoherent.c"
  rot_mul(mctr1, mcrotaSource, mcrotaCoord);
  rot_transpose(mcrotaSource, mctr1);
  rot_mul(mcrotaCoord, mctr1, mcrotrCoord);
  mctc1 = coords_set(
#line 102 "Samples_Incoherent.instr"
    0.0,
#line 102 "Samples_Incoherent.instr"
    0.0,
#line 102 "Samples_Incoherent.instr"
    8.3000);
#line 15059 "Samples_Incoherent.c"
  rot_transpose(mcrotaSource, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaCoord = coords_add(mcposaSource, mctc2);
  mctc1 = coords_sub(mcposaSource, mcposaCoord);
  mcposrCoord = rot_apply(mcrotaCoord, mctc1);
  mcDEBUG_COMPONENT("Coord", mcposaCoord, mcrotaCoord)
  mccomp_posa[2] = mcposaCoord;
  mccomp_posr[2] = mcposrCoord;
  mcNCounter[2]  = mcPCounter[2] = mcP2Counter[2] = 0;
  mcAbsorbProp[2]= 0;
    /* Component Vsample. */
  /* Setting parameters for component Vsample. */
  SIG_MESSAGE("Vsample (Init:SetPar)");
#line 91 "Samples_Incoherent.instr"
  mccVsample_radius = 0;
#line 91 "Samples_Incoherent.instr"
  mccVsample_thickness = 0;
#line 108 "Samples_Incoherent.instr"
  mccVsample_zdepth = 0.005;
#line 91 "Samples_Incoherent.instr"
  mccVsample_Vc = 13.827;
#line 107 "Samples_Incoherent.instr"
  mccVsample_sigma_abs = mcipABS;
#line 107 "Samples_Incoherent.instr"
  mccVsample_sigma_inc = mcipINC;
#line 92 "Samples_Incoherent.instr"
  mccVsample_radius_i = 0;
#line 92 "Samples_Incoherent.instr"
  mccVsample_radius_o = 0;
#line 92 "Samples_Incoherent.instr"
  mccVsample_h = 0;
#line 92 "Samples_Incoherent.instr"
  mccVsample_focus_r = 0;
#line 92 "Samples_Incoherent.instr"
  mccVsample_pack = 1;
#line 92 "Samples_Incoherent.instr"
  mccVsample_frac = 1;
#line 92 "Samples_Incoherent.instr"
  mccVsample_f_QE = 0;
#line 92 "Samples_Incoherent.instr"
  mccVsample_gamma = 0;
#line 93 "Samples_Incoherent.instr"
  mccVsample_target_x = 0;
#line 93 "Samples_Incoherent.instr"
  mccVsample_target_y = 0;
#line 93 "Samples_Incoherent.instr"
  mccVsample_target_z = 0;
#line 93 "Samples_Incoherent.instr"
  mccVsample_focus_xw = 0;
#line 93 "Samples_Incoherent.instr"
  mccVsample_focus_yh = 0;
#line 94 "Samples_Incoherent.instr"
  mccVsample_focus_aw = 0;
#line 94 "Samples_Incoherent.instr"
  mccVsample_focus_ah = 0;
#line 108 "Samples_Incoherent.instr"
  mccVsample_xwidth = 0.005;
#line 108 "Samples_Incoherent.instr"
  mccVsample_yheight = 0.005;
#line 94 "Samples_Incoherent.instr"
  mccVsample_zthick = 0;
#line 94 "Samples_Incoherent.instr"
  mccVsample_rad_sphere = 0;
#line 94 "Samples_Incoherent.instr"
  mccVsample_sig_a = 0;
#line 94 "Samples_Incoherent.instr"
  mccVsample_sig_i = 0;
#line 94 "Samples_Incoherent.instr"
  mccVsample_V0 = 0;
#line 94 "Samples_Incoherent.instr"
  mccVsample_target_index = 0;
#line 106 "Samples_Incoherent.instr"
  mccVsample_multiples = mcipV_multiples;
#line 15133 "Samples_Incoherent.c"

  SIG_MESSAGE("Vsample (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 15140 "Samples_Incoherent.c"
  rot_mul(mctr1, mcrotaCoord, mcrotaVsample);
  rot_transpose(mcrotaCoord, mctr1);
  rot_mul(mcrotaVsample, mctr1, mcrotrVsample);
  mctc1 = coords_set(
#line 109 "Samples_Incoherent.instr"
    0,
#line 109 "Samples_Incoherent.instr"
    0,
#line 109 "Samples_Incoherent.instr"
    0);
#line 15151 "Samples_Incoherent.c"
  rot_transpose(mcrotaCoord, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaVsample = coords_add(mcposaCoord, mctc2);
  mctc1 = coords_sub(mcposaCoord, mcposaVsample);
  mcposrVsample = rot_apply(mcrotaVsample, mctc1);
  mcDEBUG_COMPONENT("Vsample", mcposaVsample, mcrotaVsample)
  mccomp_posa[3] = mcposaVsample;
  mccomp_posr[3] = mcposrVsample;
  mcNCounter[3]  = mcPCounter[3] = mcP2Counter[3] = 0;
  mcAbsorbProp[3]= 0;
    /* Component PowN. */
  /* Setting parameters for component PowN. */
  SIG_MESSAGE("PowN (Init:SetPar)");
#line 111 "Samples_Incoherent.instr"
  if("") strncpy(mccPowN_reflections, "" ? "" : "", 16384); else mccPowN_reflections[0]='\0';
#line 207 "Samples_Incoherent.instr"
  if("NULL") strncpy(mccPowN_geometry, "NULL" ? "NULL" : "", 16384); else mccPowN_geometry[0]='\0';
#line 208 "Samples_Incoherent.instr"
  mccPowN_radius = 0;
#line 113 "Samples_Incoherent.instr"
  mccPowN_yheight = 0.005;
#line 113 "Samples_Incoherent.instr"
  mccPowN_xwidth = 0.005;
#line 113 "Samples_Incoherent.instr"
  mccPowN_zdepth = 0.005;
#line 208 "Samples_Incoherent.instr"
  mccPowN_thickness = 0;
#line 209 "Samples_Incoherent.instr"
  mccPowN_pack = 1;
#line 111 "Samples_Incoherent.instr"
  mccPowN_Vc = 13.827;
#line 112 "Samples_Incoherent.instr"
  mccPowN_sigma_abs = mcipABS;
#line 112 "Samples_Incoherent.instr"
  mccPowN_sigma_inc = mcipINC;
#line 209 "Samples_Incoherent.instr"
  mccPowN_delta_d_d = 0;
#line 111 "Samples_Incoherent.instr"
  mccPowN_p_inc = 0.5;
#line 111 "Samples_Incoherent.instr"
  mccPowN_p_transmit = 0.5;
#line 210 "Samples_Incoherent.instr"
  mccPowN_DW = 0;
#line 210 "Samples_Incoherent.instr"
  mccPowN_nb_atoms = 1;
#line 210 "Samples_Incoherent.instr"
  mccPowN_d_omega = 0;
#line 111 "Samples_Incoherent.instr"
  mccPowN_d_phi = 0;
#line 210 "Samples_Incoherent.instr"
  mccPowN_tth_sign = 0;
#line 210 "Samples_Incoherent.instr"
  mccPowN_p_interact = 0;
#line 211 "Samples_Incoherent.instr"
  mccPowN_concentric = 0;
#line 211 "Samples_Incoherent.instr"
  mccPowN_density = 0;
#line 211 "Samples_Incoherent.instr"
  mccPowN_weight = 0;
#line 211 "Samples_Incoherent.instr"
  mccPowN_barns = 1;
#line 211 "Samples_Incoherent.instr"
  mccPowN_Strain = 0;
#line 211 "Samples_Incoherent.instr"
  mccPowN_focus_flip = 0;
#line 211 "Samples_Incoherent.instr"
  mccPowN_target_index = 0;
#line 15219 "Samples_Incoherent.c"

  SIG_MESSAGE("PowN (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 15226 "Samples_Incoherent.c"
  rot_mul(mctr1, mcrotaCoord, mcrotaPowN);
  rot_transpose(mcrotaVsample, mctr1);
  rot_mul(mcrotaPowN, mctr1, mcrotrPowN);
  mctc1 = coords_set(
#line 114 "Samples_Incoherent.instr"
    0,
#line 114 "Samples_Incoherent.instr"
    0,
#line 114 "Samples_Incoherent.instr"
    0);
#line 15237 "Samples_Incoherent.c"
  rot_transpose(mcrotaCoord, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaPowN = coords_add(mcposaCoord, mctc2);
  mctc1 = coords_sub(mcposaVsample, mcposaPowN);
  mcposrPowN = rot_apply(mcrotaPowN, mctc1);
  mcDEBUG_COMPONENT("PowN", mcposaPowN, mcrotaPowN)
  mccomp_posa[4] = mcposaPowN;
  mccomp_posr[4] = mcposrPowN;
  mcNCounter[4]  = mcPCounter[4] = mcP2Counter[4] = 0;
  mcAbsorbProp[4]= 0;
    /* Component SX. */
  /* Setting parameters for component SX. */
  SIG_MESSAGE("SX (Init:SetPar)");
#line 116 "Samples_Incoherent.instr"
  if("") strncpy(mccSX_reflections, "" ? "" : "", 16384); else mccSX_reflections[0]='\0';
#line 282 "Samples_Incoherent.instr"
  if(0) strncpy(mccSX_geometry, 0 ? 0 : "", 16384); else mccSX_geometry[0]='\0';
#line 118 "Samples_Incoherent.instr"
  mccSX_xwidth = 0.005;
#line 118 "Samples_Incoherent.instr"
  mccSX_yheight = 0.005;
#line 118 "Samples_Incoherent.instr"
  mccSX_zdepth = 0.005;
#line 283 "Samples_Incoherent.instr"
  mccSX_radius = 0;
#line 283 "Samples_Incoherent.instr"
  mccSX_delta_d_d = 1e-4;
#line 284 "Samples_Incoherent.instr"
  mccSX_mosaic = -1;
#line 284 "Samples_Incoherent.instr"
  mccSX_mosaic_a = -1;
#line 284 "Samples_Incoherent.instr"
  mccSX_mosaic_b = -1;
#line 284 "Samples_Incoherent.instr"
  mccSX_mosaic_c = -1;
#line 285 "Samples_Incoherent.instr"
  mccSX_recip_cell = 0;
#line 117 "Samples_Incoherent.instr"
  mccSX_barns = 1;
#line 116 "Samples_Incoherent.instr"
  mccSX_ax = 3.0282;
#line 286 "Samples_Incoherent.instr"
  mccSX_ay = 0;
#line 286 "Samples_Incoherent.instr"
  mccSX_az = 0;
#line 287 "Samples_Incoherent.instr"
  mccSX_bx = 0;
#line 116 "Samples_Incoherent.instr"
  mccSX_by = 3.0282;
#line 287 "Samples_Incoherent.instr"
  mccSX_bz = 0;
#line 288 "Samples_Incoherent.instr"
  mccSX_cx = 0;
#line 288 "Samples_Incoherent.instr"
  mccSX_cy = 0;
#line 116 "Samples_Incoherent.instr"
  mccSX_cz = 3.0282 / 2;
#line 289 "Samples_Incoherent.instr"
  mccSX_p_transmit = 0.001;
#line 117 "Samples_Incoherent.instr"
  mccSX_sigma_abs = mcipABS;
#line 117 "Samples_Incoherent.instr"
  mccSX_sigma_inc = mcipINC;
#line 290 "Samples_Incoherent.instr"
  mccSX_aa = 0;
#line 290 "Samples_Incoherent.instr"
  mccSX_bb = 0;
#line 290 "Samples_Incoherent.instr"
  mccSX_cc = 0;
#line 116 "Samples_Incoherent.instr"
  mccSX_order = mcipOrder;
#line 290 "Samples_Incoherent.instr"
  mccSX_RX = 0;
#line 290 "Samples_Incoherent.instr"
  mccSX_RY = 0;
#line 290 "Samples_Incoherent.instr"
  mccSX_powder = 0;
#line 290 "Samples_Incoherent.instr"
  mccSX_PG = 0;
#line 291 "Samples_Incoherent.instr"
  mccSX_deltak = 1e-6;
#line 15319 "Samples_Incoherent.c"

  SIG_MESSAGE("SX (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 15326 "Samples_Incoherent.c"
  rot_mul(mctr1, mcrotaCoord, mcrotaSX);
  rot_transpose(mcrotaPowN, mctr1);
  rot_mul(mcrotaSX, mctr1, mcrotrSX);
  mctc1 = coords_set(
#line 119 "Samples_Incoherent.instr"
    0,
#line 119 "Samples_Incoherent.instr"
    0,
#line 119 "Samples_Incoherent.instr"
    0);
#line 15337 "Samples_Incoherent.c"
  rot_transpose(mcrotaCoord, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaSX = coords_add(mcposaCoord, mctc2);
  mctc1 = coords_sub(mcposaPowN, mcposaSX);
  mcposrSX = rot_apply(mcrotaSX, mctc1);
  mcDEBUG_COMPONENT("SX", mcposaSX, mcrotaSX)
  mccomp_posa[5] = mcposaSX;
  mccomp_posr[5] = mcposrSX;
  mcNCounter[5]  = mcPCounter[5] = mcP2Counter[5] = 0;
  mcAbsorbProp[5]= 0;
    /* Component Sqw. */
  /* Setting parameters for component Sqw. */
  SIG_MESSAGE("Sqw (Init:SetPar)");
#line 270 "Samples_Incoherent.instr"
  if(0) strncpy(mccSqw_Sqw_coh, 0 ? 0 : "", 16384); else mccSqw_Sqw_coh[0]='\0';
#line 270 "Samples_Incoherent.instr"
  if(0) strncpy(mccSqw_Sqw_inc, 0 ? 0 : "", 16384); else mccSqw_Sqw_inc[0]='\0';
#line 270 "Samples_Incoherent.instr"
  if(0) strncpy(mccSqw_geometry, 0 ? 0 : "", 16384); else mccSqw_geometry[0]='\0';
#line 271 "Samples_Incoherent.instr"
  mccSqw_radius = 0;
#line 271 "Samples_Incoherent.instr"
  mccSqw_thickness = 0;
#line 123 "Samples_Incoherent.instr"
  mccSqw_xwidth = 0.005;
#line 123 "Samples_Incoherent.instr"
  mccSqw_yheight = 0.005;
#line 123 "Samples_Incoherent.instr"
  mccSqw_zdepth = 0.005;
#line 273 "Samples_Incoherent.instr"
  mccSqw_threshold = 1e-20;
#line 121 "Samples_Incoherent.instr"
  mccSqw_order = mcipOrder;
#line 273 "Samples_Incoherent.instr"
  mccSqw_T = 0;
#line 273 "Samples_Incoherent.instr"
  mccSqw_verbose = 1;
#line 273 "Samples_Incoherent.instr"
  mccSqw_d_phi = 0;
#line 273 "Samples_Incoherent.instr"
  mccSqw_concentric = 0;
#line 121 "Samples_Incoherent.instr"
  mccSqw_rho = 1 / 13.827;
#line 122 "Samples_Incoherent.instr"
  mccSqw_sigma_abs = mcipABS;
#line 122 "Samples_Incoherent.instr"
  mccSqw_sigma_coh = 0;
#line 122 "Samples_Incoherent.instr"
  mccSqw_sigma_inc = mcipINC;
#line 274 "Samples_Incoherent.instr"
  mccSqw_classical = -1;
#line 275 "Samples_Incoherent.instr"
  mccSqw_powder_Dd = 0;
#line 275 "Samples_Incoherent.instr"
  mccSqw_powder_DW = 0;
#line 275 "Samples_Incoherent.instr"
  mccSqw_powder_Vc = 0;
#line 275 "Samples_Incoherent.instr"
  mccSqw_density = 0;
#line 275 "Samples_Incoherent.instr"
  mccSqw_weight = 0;
#line 276 "Samples_Incoherent.instr"
  mccSqw_p_interact = -1;
#line 276 "Samples_Incoherent.instr"
  mccSqw_norm = -1;
#line 276 "Samples_Incoherent.instr"
  mccSqw_powder_barns = 1;
#line 276 "Samples_Incoherent.instr"
  if("Frommhold") strncpy(mccSqw_quantum_correction, "Frommhold" ? "Frommhold" : "", 16384); else mccSqw_quantum_correction[0]='\0';
#line 15407 "Samples_Incoherent.c"

  SIG_MESSAGE("Sqw (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 15414 "Samples_Incoherent.c"
  rot_mul(mctr1, mcrotaCoord, mcrotaSqw);
  rot_transpose(mcrotaSX, mctr1);
  rot_mul(mcrotaSqw, mctr1, mcrotrSqw);
  mctc1 = coords_set(
#line 124 "Samples_Incoherent.instr"
    0,
#line 124 "Samples_Incoherent.instr"
    0,
#line 124 "Samples_Incoherent.instr"
    0);
#line 15425 "Samples_Incoherent.c"
  rot_transpose(mcrotaCoord, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaSqw = coords_add(mcposaCoord, mctc2);
  mctc1 = coords_sub(mcposaSX, mcposaSqw);
  mcposrSqw = rot_apply(mcrotaSqw, mctc1);
  mcDEBUG_COMPONENT("Sqw", mcposaSqw, mcrotaSqw)
  mccomp_posa[6] = mcposaSqw;
  mccomp_posr[6] = mcposrSqw;
  mcNCounter[6]  = mcPCounter[6] = mcP2Counter[6] = 0;
  mcAbsorbProp[6]= 0;
    /* Component IncS. */
  /* Setting parameters for component IncS. */
  SIG_MESSAGE("IncS (Init:SetPar)");
#line 118 "Samples_Incoherent.instr"
  if(0) strncpy(mccIncS_geometry, 0 ? 0 : "", 16384); else mccIncS_geometry[0]='\0';
#line 118 "Samples_Incoherent.instr"
  mccIncS_radius = 0;
#line 128 "Samples_Incoherent.instr"
  mccIncS_xwidth = 0.005;
#line 128 "Samples_Incoherent.instr"
  mccIncS_yheight = 0.005;
#line 128 "Samples_Incoherent.instr"
  mccIncS_zdepth = 0.005;
#line 118 "Samples_Incoherent.instr"
  mccIncS_thickness = 0;
#line 119 "Samples_Incoherent.instr"
  mccIncS_target_x = 0;
#line 119 "Samples_Incoherent.instr"
  mccIncS_target_y = 0;
#line 119 "Samples_Incoherent.instr"
  mccIncS_target_z = 0;
#line 119 "Samples_Incoherent.instr"
  mccIncS_focus_r = 0;
#line 120 "Samples_Incoherent.instr"
  mccIncS_focus_xw = 0;
#line 120 "Samples_Incoherent.instr"
  mccIncS_focus_yh = 0;
#line 120 "Samples_Incoherent.instr"
  mccIncS_focus_aw = 0;
#line 120 "Samples_Incoherent.instr"
  mccIncS_focus_ah = 0;
#line 120 "Samples_Incoherent.instr"
  mccIncS_target_index = 0;
#line 121 "Samples_Incoherent.instr"
  mccIncS_pack = 1;
#line 121 "Samples_Incoherent.instr"
  mccIncS_p_interact = 1;
#line 121 "Samples_Incoherent.instr"
  mccIncS_f_QE = 0;
#line 121 "Samples_Incoherent.instr"
  mccIncS_gamma = 0;
#line 121 "Samples_Incoherent.instr"
  mccIncS_Etrans = 0;
#line 121 "Samples_Incoherent.instr"
  mccIncS_deltaE = 0;
#line 127 "Samples_Incoherent.instr"
  mccIncS_sigma_abs = mcipABS;
#line 127 "Samples_Incoherent.instr"
  mccIncS_sigma_inc = mcipINC;
#line 122 "Samples_Incoherent.instr"
  mccIncS_Vc = 13.827;
#line 122 "Samples_Incoherent.instr"
  mccIncS_concentric = 0;
#line 126 "Samples_Incoherent.instr"
  mccIncS_order = mcipOrder;
#line 15491 "Samples_Incoherent.c"

  SIG_MESSAGE("IncS (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 15498 "Samples_Incoherent.c"
  rot_mul(mctr1, mcrotaCoord, mcrotaIncS);
  rot_transpose(mcrotaSqw, mctr1);
  rot_mul(mcrotaIncS, mctr1, mcrotrIncS);
  mctc1 = coords_set(
#line 129 "Samples_Incoherent.instr"
    0,
#line 129 "Samples_Incoherent.instr"
    0,
#line 129 "Samples_Incoherent.instr"
    0);
#line 15509 "Samples_Incoherent.c"
  rot_transpose(mcrotaCoord, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaIncS = coords_add(mcposaCoord, mctc2);
  mctc1 = coords_sub(mcposaSqw, mcposaIncS);
  mcposrIncS = rot_apply(mcrotaIncS, mctc1);
  mcDEBUG_COMPONENT("IncS", mcposaIncS, mcrotaIncS)
  mccomp_posa[7] = mcposaIncS;
  mccomp_posr[7] = mcposrIncS;
  mcNCounter[7]  = mcPCounter[7] = mcP2Counter[7] = 0;
  mcAbsorbProp[7]= 0;
    /* Component PowNg. */
  /* Setting parameters for component PowNg. */
  SIG_MESSAGE("PowNg (Init:SetPar)");
#line 133 "Samples_Incoherent.instr"
  if("") strncpy(mccPowNg_reflections, "" ? "" : "", 16384); else mccPowNg_reflections[0]='\0';
#line 135 "Samples_Incoherent.instr"
  if("cube.off") strncpy(mccPowNg_geometry, "cube.off" ? "cube.off" : "", 16384); else mccPowNg_geometry[0]='\0';
#line 208 "Samples_Incoherent.instr"
  mccPowNg_radius = 0;
#line 208 "Samples_Incoherent.instr"
  mccPowNg_yheight = 0;
#line 135 "Samples_Incoherent.instr"
  mccPowNg_xwidth = 0.005;
#line 208 "Samples_Incoherent.instr"
  mccPowNg_zdepth = 0;
#line 208 "Samples_Incoherent.instr"
  mccPowNg_thickness = 0;
#line 209 "Samples_Incoherent.instr"
  mccPowNg_pack = 1;
#line 133 "Samples_Incoherent.instr"
  mccPowNg_Vc = 13.827;
#line 134 "Samples_Incoherent.instr"
  mccPowNg_sigma_abs = mcipABS;
#line 134 "Samples_Incoherent.instr"
  mccPowNg_sigma_inc = mcipINC;
#line 209 "Samples_Incoherent.instr"
  mccPowNg_delta_d_d = 0;
#line 133 "Samples_Incoherent.instr"
  mccPowNg_p_inc = 0.5;
#line 133 "Samples_Incoherent.instr"
  mccPowNg_p_transmit = 0.5;
#line 210 "Samples_Incoherent.instr"
  mccPowNg_DW = 0;
#line 210 "Samples_Incoherent.instr"
  mccPowNg_nb_atoms = 1;
#line 210 "Samples_Incoherent.instr"
  mccPowNg_d_omega = 0;
#line 133 "Samples_Incoherent.instr"
  mccPowNg_d_phi = 0;
#line 210 "Samples_Incoherent.instr"
  mccPowNg_tth_sign = 0;
#line 210 "Samples_Incoherent.instr"
  mccPowNg_p_interact = 0;
#line 211 "Samples_Incoherent.instr"
  mccPowNg_concentric = 0;
#line 211 "Samples_Incoherent.instr"
  mccPowNg_density = 0;
#line 211 "Samples_Incoherent.instr"
  mccPowNg_weight = 0;
#line 211 "Samples_Incoherent.instr"
  mccPowNg_barns = 1;
#line 211 "Samples_Incoherent.instr"
  mccPowNg_Strain = 0;
#line 211 "Samples_Incoherent.instr"
  mccPowNg_focus_flip = 0;
#line 211 "Samples_Incoherent.instr"
  mccPowNg_target_index = 0;
#line 15577 "Samples_Incoherent.c"

  SIG_MESSAGE("PowNg (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 15584 "Samples_Incoherent.c"
  rot_mul(mctr1, mcrotaCoord, mcrotaPowNg);
  rot_transpose(mcrotaIncS, mctr1);
  rot_mul(mcrotaPowNg, mctr1, mcrotrPowNg);
  mctc1 = coords_set(
#line 136 "Samples_Incoherent.instr"
    0,
#line 136 "Samples_Incoherent.instr"
    0,
#line 136 "Samples_Incoherent.instr"
    0);
#line 15595 "Samples_Incoherent.c"
  rot_transpose(mcrotaCoord, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaPowNg = coords_add(mcposaCoord, mctc2);
  mctc1 = coords_sub(mcposaIncS, mcposaPowNg);
  mcposrPowNg = rot_apply(mcrotaPowNg, mctc1);
  mcDEBUG_COMPONENT("PowNg", mcposaPowNg, mcrotaPowNg)
  mccomp_posa[8] = mcposaPowNg;
  mccomp_posr[8] = mcposrPowNg;
  mcNCounter[8]  = mcPCounter[8] = mcP2Counter[8] = 0;
  mcAbsorbProp[8]= 0;
    /* Component SXg. */
  /* Setting parameters for component SXg. */
  SIG_MESSAGE("SXg (Init:SetPar)");
#line 138 "Samples_Incoherent.instr"
  if("") strncpy(mccSXg_reflections, "" ? "" : "", 16384); else mccSXg_reflections[0]='\0';
#line 140 "Samples_Incoherent.instr"
  if("cube.off") strncpy(mccSXg_geometry, "cube.off" ? "cube.off" : "", 16384); else mccSXg_geometry[0]='\0';
#line 140 "Samples_Incoherent.instr"
  mccSXg_xwidth = 0.005;
#line 283 "Samples_Incoherent.instr"
  mccSXg_yheight = 0;
#line 283 "Samples_Incoherent.instr"
  mccSXg_zdepth = 0;
#line 283 "Samples_Incoherent.instr"
  mccSXg_radius = 0;
#line 283 "Samples_Incoherent.instr"
  mccSXg_delta_d_d = 1e-4;
#line 284 "Samples_Incoherent.instr"
  mccSXg_mosaic = -1;
#line 284 "Samples_Incoherent.instr"
  mccSXg_mosaic_a = -1;
#line 284 "Samples_Incoherent.instr"
  mccSXg_mosaic_b = -1;
#line 284 "Samples_Incoherent.instr"
  mccSXg_mosaic_c = -1;
#line 285 "Samples_Incoherent.instr"
  mccSXg_recip_cell = 0;
#line 139 "Samples_Incoherent.instr"
  mccSXg_barns = 1;
#line 138 "Samples_Incoherent.instr"
  mccSXg_ax = 3.0282;
#line 286 "Samples_Incoherent.instr"
  mccSXg_ay = 0;
#line 286 "Samples_Incoherent.instr"
  mccSXg_az = 0;
#line 287 "Samples_Incoherent.instr"
  mccSXg_bx = 0;
#line 138 "Samples_Incoherent.instr"
  mccSXg_by = 3.0282;
#line 287 "Samples_Incoherent.instr"
  mccSXg_bz = 0;
#line 288 "Samples_Incoherent.instr"
  mccSXg_cx = 0;
#line 288 "Samples_Incoherent.instr"
  mccSXg_cy = 0;
#line 138 "Samples_Incoherent.instr"
  mccSXg_cz = 3.0282 / 2;
#line 289 "Samples_Incoherent.instr"
  mccSXg_p_transmit = 0.001;
#line 139 "Samples_Incoherent.instr"
  mccSXg_sigma_abs = mcipABS;
#line 139 "Samples_Incoherent.instr"
  mccSXg_sigma_inc = mcipINC;
#line 290 "Samples_Incoherent.instr"
  mccSXg_aa = 0;
#line 290 "Samples_Incoherent.instr"
  mccSXg_bb = 0;
#line 290 "Samples_Incoherent.instr"
  mccSXg_cc = 0;
#line 138 "Samples_Incoherent.instr"
  mccSXg_order = mcipOrder;
#line 290 "Samples_Incoherent.instr"
  mccSXg_RX = 0;
#line 290 "Samples_Incoherent.instr"
  mccSXg_RY = 0;
#line 290 "Samples_Incoherent.instr"
  mccSXg_powder = 0;
#line 290 "Samples_Incoherent.instr"
  mccSXg_PG = 0;
#line 291 "Samples_Incoherent.instr"
  mccSXg_deltak = 1e-6;
#line 15677 "Samples_Incoherent.c"

  SIG_MESSAGE("SXg (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 15684 "Samples_Incoherent.c"
  rot_mul(mctr1, mcrotaCoord, mcrotaSXg);
  rot_transpose(mcrotaPowNg, mctr1);
  rot_mul(mcrotaSXg, mctr1, mcrotrSXg);
  mctc1 = coords_set(
#line 141 "Samples_Incoherent.instr"
    0,
#line 141 "Samples_Incoherent.instr"
    0,
#line 141 "Samples_Incoherent.instr"
    0);
#line 15695 "Samples_Incoherent.c"
  rot_transpose(mcrotaCoord, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaSXg = coords_add(mcposaCoord, mctc2);
  mctc1 = coords_sub(mcposaPowNg, mcposaSXg);
  mcposrSXg = rot_apply(mcrotaSXg, mctc1);
  mcDEBUG_COMPONENT("SXg", mcposaSXg, mcrotaSXg)
  mccomp_posa[9] = mcposaSXg;
  mccomp_posr[9] = mcposrSXg;
  mcNCounter[9]  = mcPCounter[9] = mcP2Counter[9] = 0;
  mcAbsorbProp[9]= 0;
    /* Component Sqwg. */
  /* Setting parameters for component Sqwg. */
  SIG_MESSAGE("Sqwg (Init:SetPar)");
#line 270 "Samples_Incoherent.instr"
  if(0) strncpy(mccSqwg_Sqw_coh, 0 ? 0 : "", 16384); else mccSqwg_Sqw_coh[0]='\0';
#line 270 "Samples_Incoherent.instr"
  if(0) strncpy(mccSqwg_Sqw_inc, 0 ? 0 : "", 16384); else mccSqwg_Sqw_inc[0]='\0';
#line 145 "Samples_Incoherent.instr"
  if("cube.off") strncpy(mccSqwg_geometry, "cube.off" ? "cube.off" : "", 16384); else mccSqwg_geometry[0]='\0';
#line 271 "Samples_Incoherent.instr"
  mccSqwg_radius = 0;
#line 271 "Samples_Incoherent.instr"
  mccSqwg_thickness = 0;
#line 145 "Samples_Incoherent.instr"
  mccSqwg_xwidth = 0.005;
#line 272 "Samples_Incoherent.instr"
  mccSqwg_yheight = 0;
#line 272 "Samples_Incoherent.instr"
  mccSqwg_zdepth = 0;
#line 273 "Samples_Incoherent.instr"
  mccSqwg_threshold = 1e-20;
#line 143 "Samples_Incoherent.instr"
  mccSqwg_order = mcipOrder;
#line 273 "Samples_Incoherent.instr"
  mccSqwg_T = 0;
#line 273 "Samples_Incoherent.instr"
  mccSqwg_verbose = 1;
#line 273 "Samples_Incoherent.instr"
  mccSqwg_d_phi = 0;
#line 273 "Samples_Incoherent.instr"
  mccSqwg_concentric = 0;
#line 143 "Samples_Incoherent.instr"
  mccSqwg_rho = 1 / 13.827;
#line 144 "Samples_Incoherent.instr"
  mccSqwg_sigma_abs = mcipABS;
#line 144 "Samples_Incoherent.instr"
  mccSqwg_sigma_coh = 0;
#line 144 "Samples_Incoherent.instr"
  mccSqwg_sigma_inc = mcipINC;
#line 274 "Samples_Incoherent.instr"
  mccSqwg_classical = -1;
#line 275 "Samples_Incoherent.instr"
  mccSqwg_powder_Dd = 0;
#line 275 "Samples_Incoherent.instr"
  mccSqwg_powder_DW = 0;
#line 275 "Samples_Incoherent.instr"
  mccSqwg_powder_Vc = 0;
#line 275 "Samples_Incoherent.instr"
  mccSqwg_density = 0;
#line 275 "Samples_Incoherent.instr"
  mccSqwg_weight = 0;
#line 276 "Samples_Incoherent.instr"
  mccSqwg_p_interact = -1;
#line 276 "Samples_Incoherent.instr"
  mccSqwg_norm = -1;
#line 276 "Samples_Incoherent.instr"
  mccSqwg_powder_barns = 1;
#line 276 "Samples_Incoherent.instr"
  if("Frommhold") strncpy(mccSqwg_quantum_correction, "Frommhold" ? "Frommhold" : "", 16384); else mccSqwg_quantum_correction[0]='\0';
#line 15765 "Samples_Incoherent.c"

  SIG_MESSAGE("Sqwg (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 15772 "Samples_Incoherent.c"
  rot_mul(mctr1, mcrotaCoord, mcrotaSqwg);
  rot_transpose(mcrotaSXg, mctr1);
  rot_mul(mcrotaSqwg, mctr1, mcrotrSqwg);
  mctc1 = coords_set(
#line 146 "Samples_Incoherent.instr"
    0,
#line 146 "Samples_Incoherent.instr"
    0,
#line 146 "Samples_Incoherent.instr"
    0);
#line 15783 "Samples_Incoherent.c"
  rot_transpose(mcrotaCoord, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaSqwg = coords_add(mcposaCoord, mctc2);
  mctc1 = coords_sub(mcposaSXg, mcposaSqwg);
  mcposrSqwg = rot_apply(mcrotaSqwg, mctc1);
  mcDEBUG_COMPONENT("Sqwg", mcposaSqwg, mcrotaSqwg)
  mccomp_posa[10] = mcposaSqwg;
  mccomp_posr[10] = mcposrSqwg;
  mcNCounter[10]  = mcPCounter[10] = mcP2Counter[10] = 0;
  mcAbsorbProp[10]= 0;
    /* Component IncSg. */
  /* Setting parameters for component IncSg. */
  SIG_MESSAGE("IncSg (Init:SetPar)");
#line 150 "Samples_Incoherent.instr"
  if("cube.off") strncpy(mccIncSg_geometry, "cube.off" ? "cube.off" : "", 16384); else mccIncSg_geometry[0]='\0';
#line 118 "Samples_Incoherent.instr"
  mccIncSg_radius = 0;
#line 150 "Samples_Incoherent.instr"
  mccIncSg_xwidth = 0.005;
#line 118 "Samples_Incoherent.instr"
  mccIncSg_yheight = 0;
#line 118 "Samples_Incoherent.instr"
  mccIncSg_zdepth = 0;
#line 118 "Samples_Incoherent.instr"
  mccIncSg_thickness = 0;
#line 119 "Samples_Incoherent.instr"
  mccIncSg_target_x = 0;
#line 119 "Samples_Incoherent.instr"
  mccIncSg_target_y = 0;
#line 119 "Samples_Incoherent.instr"
  mccIncSg_target_z = 0;
#line 119 "Samples_Incoherent.instr"
  mccIncSg_focus_r = 0;
#line 120 "Samples_Incoherent.instr"
  mccIncSg_focus_xw = 0;
#line 120 "Samples_Incoherent.instr"
  mccIncSg_focus_yh = 0;
#line 120 "Samples_Incoherent.instr"
  mccIncSg_focus_aw = 0;
#line 120 "Samples_Incoherent.instr"
  mccIncSg_focus_ah = 0;
#line 120 "Samples_Incoherent.instr"
  mccIncSg_target_index = 0;
#line 121 "Samples_Incoherent.instr"
  mccIncSg_pack = 1;
#line 121 "Samples_Incoherent.instr"
  mccIncSg_p_interact = 1;
#line 121 "Samples_Incoherent.instr"
  mccIncSg_f_QE = 0;
#line 121 "Samples_Incoherent.instr"
  mccIncSg_gamma = 0;
#line 121 "Samples_Incoherent.instr"
  mccIncSg_Etrans = 0;
#line 121 "Samples_Incoherent.instr"
  mccIncSg_deltaE = 0;
#line 149 "Samples_Incoherent.instr"
  mccIncSg_sigma_abs = mcipABS;
#line 149 "Samples_Incoherent.instr"
  mccIncSg_sigma_inc = mcipINC;
#line 122 "Samples_Incoherent.instr"
  mccIncSg_Vc = 13.827;
#line 122 "Samples_Incoherent.instr"
  mccIncSg_concentric = 0;
#line 148 "Samples_Incoherent.instr"
  mccIncSg_order = mcipOrder;
#line 15849 "Samples_Incoherent.c"

  SIG_MESSAGE("IncSg (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 15856 "Samples_Incoherent.c"
  rot_mul(mctr1, mcrotaCoord, mcrotaIncSg);
  rot_transpose(mcrotaSqwg, mctr1);
  rot_mul(mcrotaIncSg, mctr1, mcrotrIncSg);
  mctc1 = coords_set(
#line 151 "Samples_Incoherent.instr"
    0,
#line 151 "Samples_Incoherent.instr"
    0,
#line 151 "Samples_Incoherent.instr"
    0);
#line 15867 "Samples_Incoherent.c"
  rot_transpose(mcrotaCoord, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaIncSg = coords_add(mcposaCoord, mctc2);
  mctc1 = coords_sub(mcposaSqwg, mcposaIncSg);
  mcposrIncSg = rot_apply(mcrotaIncSg, mctc1);
  mcDEBUG_COMPONENT("IncSg", mcposaIncSg, mcrotaIncSg)
  mccomp_posa[11] = mcposaIncSg;
  mccomp_posr[11] = mcposrIncSg;
  mcNCounter[11]  = mcPCounter[11] = mcP2Counter[11] = 0;
  mcAbsorbProp[11]= 0;
    /* Component SqwV. */
  /* Setting parameters for component SqwV. */
  SIG_MESSAGE("SqwV (Init:SetPar)");
#line 155 "Samples_Incoherent.instr"
  if("V.laz") strncpy(mccSqwV_Sqw_coh, "V.laz" ? "V.laz" : "", 16384); else mccSqwV_Sqw_coh[0]='\0';
#line 270 "Samples_Incoherent.instr"
  if(0) strncpy(mccSqwV_Sqw_inc, 0 ? 0 : "", 16384); else mccSqwV_Sqw_inc[0]='\0';
#line 270 "Samples_Incoherent.instr"
  if(0) strncpy(mccSqwV_geometry, 0 ? 0 : "", 16384); else mccSqwV_geometry[0]='\0';
#line 271 "Samples_Incoherent.instr"
  mccSqwV_radius = 0;
#line 271 "Samples_Incoherent.instr"
  mccSqwV_thickness = 0;
#line 156 "Samples_Incoherent.instr"
  mccSqwV_xwidth = 0.005;
#line 156 "Samples_Incoherent.instr"
  mccSqwV_yheight = 0.005;
#line 156 "Samples_Incoherent.instr"
  mccSqwV_zdepth = 0.005;
#line 273 "Samples_Incoherent.instr"
  mccSqwV_threshold = 1e-20;
#line 273 "Samples_Incoherent.instr"
  mccSqwV_order = 0;
#line 273 "Samples_Incoherent.instr"
  mccSqwV_T = 0;
#line 273 "Samples_Incoherent.instr"
  mccSqwV_verbose = 1;
#line 273 "Samples_Incoherent.instr"
  mccSqwV_d_phi = 0;
#line 273 "Samples_Incoherent.instr"
  mccSqwV_concentric = 0;
#line 274 "Samples_Incoherent.instr"
  mccSqwV_rho = 0;
#line 274 "Samples_Incoherent.instr"
  mccSqwV_sigma_abs = 0;
#line 274 "Samples_Incoherent.instr"
  mccSqwV_sigma_coh = 0;
#line 274 "Samples_Incoherent.instr"
  mccSqwV_sigma_inc = 0;
#line 274 "Samples_Incoherent.instr"
  mccSqwV_classical = -1;
#line 275 "Samples_Incoherent.instr"
  mccSqwV_powder_Dd = 0;
#line 275 "Samples_Incoherent.instr"
  mccSqwV_powder_DW = 0;
#line 275 "Samples_Incoherent.instr"
  mccSqwV_powder_Vc = 0;
#line 275 "Samples_Incoherent.instr"
  mccSqwV_density = 0;
#line 275 "Samples_Incoherent.instr"
  mccSqwV_weight = 0;
#line 276 "Samples_Incoherent.instr"
  mccSqwV_p_interact = -1;
#line 276 "Samples_Incoherent.instr"
  mccSqwV_norm = -1;
#line 276 "Samples_Incoherent.instr"
  mccSqwV_powder_barns = 1;
#line 276 "Samples_Incoherent.instr"
  if("Frommhold") strncpy(mccSqwV_quantum_correction, "Frommhold" ? "Frommhold" : "", 16384); else mccSqwV_quantum_correction[0]='\0';
#line 15937 "Samples_Incoherent.c"

  SIG_MESSAGE("SqwV (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 15944 "Samples_Incoherent.c"
  rot_mul(mctr1, mcrotaCoord, mcrotaSqwV);
  rot_transpose(mcrotaIncSg, mctr1);
  rot_mul(mcrotaSqwV, mctr1, mcrotrSqwV);
  mctc1 = coords_set(
#line 157 "Samples_Incoherent.instr"
    0,
#line 157 "Samples_Incoherent.instr"
    0,
#line 157 "Samples_Incoherent.instr"
    0);
#line 15955 "Samples_Incoherent.c"
  rot_transpose(mcrotaCoord, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaSqwV = coords_add(mcposaCoord, mctc2);
  mctc1 = coords_sub(mcposaIncSg, mcposaSqwV);
  mcposrSqwV = rot_apply(mcrotaSqwV, mctc1);
  mcDEBUG_COMPONENT("SqwV", mcposaSqwV, mcrotaSqwV)
  mccomp_posa[12] = mcposaSqwV;
  mccomp_posr[12] = mcposrSqwV;
  mcNCounter[12]  = mcPCounter[12] = mcP2Counter[12] = 0;
  mcAbsorbProp[12]= 0;
    /* Component Dirbeam. */
  /* Setting parameters for component Dirbeam. */
  SIG_MESSAGE("Dirbeam (Init:SetPar)");
#line 161 "Samples_Incoherent.instr"
  mccDirbeam_nx = 128;
#line 161 "Samples_Incoherent.instr"
  mccDirbeam_ny = 128;
#line 161 "Samples_Incoherent.instr"
  if("beamstopPSD") strncpy(mccDirbeam_filename, "beamstopPSD" ? "beamstopPSD" : "", 16384); else mccDirbeam_filename[0]='\0';
#line 50 "Samples_Incoherent.instr"
  mccDirbeam_xmin = -0.05;
#line 50 "Samples_Incoherent.instr"
  mccDirbeam_xmax = 0.05;
#line 50 "Samples_Incoherent.instr"
  mccDirbeam_ymin = -0.05;
#line 50 "Samples_Incoherent.instr"
  mccDirbeam_ymax = 0.05;
#line 161 "Samples_Incoherent.instr"
  mccDirbeam_xwidth = 0.01;
#line 161 "Samples_Incoherent.instr"
  mccDirbeam_yheight = 0.01;
#line 51 "Samples_Incoherent.instr"
  mccDirbeam_restore_neutron = 0;
#line 15989 "Samples_Incoherent.c"

  SIG_MESSAGE("Dirbeam (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 15996 "Samples_Incoherent.c"
  rot_mul(mctr1, mcrotaSqw, mcrotaDirbeam);
  rot_transpose(mcrotaSqwV, mctr1);
  rot_mul(mcrotaDirbeam, mctr1, mcrotrDirbeam);
  mctc1 = coords_set(
#line 162 "Samples_Incoherent.instr"
    0,
#line 162 "Samples_Incoherent.instr"
    0,
#line 162 "Samples_Incoherent.instr"
    0.039);
#line 16007 "Samples_Incoherent.c"
  rot_transpose(mcrotaSqw, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaDirbeam = coords_add(mcposaSqw, mctc2);
  mctc1 = coords_sub(mcposaSqwV, mcposaDirbeam);
  mcposrDirbeam = rot_apply(mcrotaDirbeam, mctc1);
  mcDEBUG_COMPONENT("Dirbeam", mcposaDirbeam, mcrotaDirbeam)
  mccomp_posa[13] = mcposaDirbeam;
  mccomp_posr[13] = mcposrDirbeam;
  mcNCounter[13]  = mcPCounter[13] = mcP2Counter[13] = 0;
  mcAbsorbProp[13]= 0;
    /* Component Stop. */
  /* Setting parameters for component Stop. */
  SIG_MESSAGE("Stop (Init:SetPar)");
#line 44 "Samples_Incoherent.instr"
  mccStop_xmin = -0.05;
#line 44 "Samples_Incoherent.instr"
  mccStop_xmax = 0.05;
#line 44 "Samples_Incoherent.instr"
  mccStop_ymin = -0.05;
#line 44 "Samples_Incoherent.instr"
  mccStop_ymax = 0.05;
#line 45 "Samples_Incoherent.instr"
  mccStop_xwidth = 0;
#line 45 "Samples_Incoherent.instr"
  mccStop_yheight = 0;
#line 164 "Samples_Incoherent.instr"
  mccStop_radius = 0.005;
#line 16035 "Samples_Incoherent.c"

  SIG_MESSAGE("Stop (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 16042 "Samples_Incoherent.c"
  rot_mul(mctr1, mcrotaSqw, mcrotaStop);
  rot_transpose(mcrotaDirbeam, mctr1);
  rot_mul(mcrotaStop, mctr1, mcrotrStop);
  mctc1 = coords_set(
#line 165 "Samples_Incoherent.instr"
    0,
#line 165 "Samples_Incoherent.instr"
    0,
#line 165 "Samples_Incoherent.instr"
    0.04);
#line 16053 "Samples_Incoherent.c"
  rot_transpose(mcrotaSqw, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaStop = coords_add(mcposaSqw, mctc2);
  mctc1 = coords_sub(mcposaDirbeam, mcposaStop);
  mcposrStop = rot_apply(mcrotaStop, mctc1);
  mcDEBUG_COMPONENT("Stop", mcposaStop, mcrotaStop)
  mccomp_posa[14] = mcposaStop;
  mccomp_posr[14] = mcposrStop;
  mcNCounter[14]  = mcPCounter[14] = mcP2Counter[14] = 0;
  mcAbsorbProp[14]= 0;
    /* Component Tof_Sphere_4pi. */
  /* Setting parameters for component Tof_Sphere_4pi. */
  SIG_MESSAGE("Tof_Sphere_4pi (Init:SetPar)");
#line 168 "Samples_Incoherent.instr"
  mccTof_Sphere_4pi_xwidth = 0.1;
#line 203 "Samples_Incoherent.instr"
  mccTof_Sphere_4pi_yheight = 0;
#line 203 "Samples_Incoherent.instr"
  mccTof_Sphere_4pi_zdepth = 0;
#line 204 "Samples_Incoherent.instr"
  mccTof_Sphere_4pi_xmin = 0;
#line 204 "Samples_Incoherent.instr"
  mccTof_Sphere_4pi_xmax = 0;
#line 204 "Samples_Incoherent.instr"
  mccTof_Sphere_4pi_ymin = 0;
#line 204 "Samples_Incoherent.instr"
  mccTof_Sphere_4pi_ymax = 0;
#line 204 "Samples_Incoherent.instr"
  mccTof_Sphere_4pi_zmin = 0;
#line 204 "Samples_Incoherent.instr"
  mccTof_Sphere_4pi_zmax = 0;
#line 205 "Samples_Incoherent.instr"
  mccTof_Sphere_4pi_bins = 0;
#line 205 "Samples_Incoherent.instr"
  mccTof_Sphere_4pi_min = -1e40;
#line 205 "Samples_Incoherent.instr"
  mccTof_Sphere_4pi_max = 1e40;
#line 205 "Samples_Incoherent.instr"
  mccTof_Sphere_4pi_restore_neutron = 0;
#line 205 "Samples_Incoherent.instr"
  mccTof_Sphere_4pi_radius = 0;
#line 168 "Samples_Incoherent.instr"
  if("sphere abs lambda, limits=[0.38,7.13], bins=2000, file=TOF") strncpy(mccTof_Sphere_4pi_options, "sphere abs lambda, limits=[0.38,7.13], bins=2000, file=TOF" ? "sphere abs lambda, limits=[0.38,7.13], bins=2000, file=TOF" : "", 16384); else mccTof_Sphere_4pi_options[0]='\0';
#line 206 "Samples_Incoherent.instr"
  if("NULL") strncpy(mccTof_Sphere_4pi_filename, "NULL" ? "NULL" : "", 16384); else mccTof_Sphere_4pi_filename[0]='\0';
#line 206 "Samples_Incoherent.instr"
  if("NULL") strncpy(mccTof_Sphere_4pi_geometry, "NULL" ? "NULL" : "", 16384); else mccTof_Sphere_4pi_geometry[0]='\0';
#line 207 "Samples_Incoherent.instr"
  if("NULL") strncpy(mccTof_Sphere_4pi_username1, "NULL" ? "NULL" : "", 16384); else mccTof_Sphere_4pi_username1[0]='\0';
#line 207 "Samples_Incoherent.instr"
  if("NULL") strncpy(mccTof_Sphere_4pi_username2, "NULL" ? "NULL" : "", 16384); else mccTof_Sphere_4pi_username2[0]='\0';
#line 207 "Samples_Incoherent.instr"
  if("NULL") strncpy(mccTof_Sphere_4pi_username3, "NULL" ? "NULL" : "", 16384); else mccTof_Sphere_4pi_username3[0]='\0';
#line 208 "Samples_Incoherent.instr"
  mccTof_Sphere_4pi_nowritefile = 0;
#line 16109 "Samples_Incoherent.c"

  SIG_MESSAGE("Tof_Sphere_4pi (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 16116 "Samples_Incoherent.c"
  rot_mul(mctr1, mcrotaCoord, mcrotaTof_Sphere_4pi);
  rot_transpose(mcrotaStop, mctr1);
  rot_mul(mcrotaTof_Sphere_4pi, mctr1, mcrotrTof_Sphere_4pi);
  mctc1 = coords_set(
#line 169 "Samples_Incoherent.instr"
    0,
#line 169 "Samples_Incoherent.instr"
    0,
#line 169 "Samples_Incoherent.instr"
    0);
#line 16127 "Samples_Incoherent.c"
  rot_transpose(mcrotaCoord, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaTof_Sphere_4pi = coords_add(mcposaCoord, mctc2);
  mctc1 = coords_sub(mcposaStop, mcposaTof_Sphere_4pi);
  mcposrTof_Sphere_4pi = rot_apply(mcrotaTof_Sphere_4pi, mctc1);
  mcDEBUG_COMPONENT("Tof_Sphere_4pi", mcposaTof_Sphere_4pi, mcrotaTof_Sphere_4pi)
  mccomp_posa[15] = mcposaTof_Sphere_4pi;
  mccomp_posr[15] = mcposrTof_Sphere_4pi;
  mcNCounter[15]  = mcPCounter[15] = mcP2Counter[15] = 0;
  mcAbsorbProp[15]= 0;
    /* Component PSD_Sphere_4pi. */
  /* Setting parameters for component PSD_Sphere_4pi. */
  SIG_MESSAGE("PSD_Sphere_4pi (Init:SetPar)");
#line 171 "Samples_Incoherent.instr"
  if("PSD") strncpy(mccPSD_Sphere_4pi_filename, "PSD" ? "PSD" : "", 16384); else mccPSD_Sphere_4pi_filename[0]='\0';
#line 171 "Samples_Incoherent.instr"
  mccPSD_Sphere_4pi_radius = 0.051;
#line 49 "Samples_Incoherent.instr"
  mccPSD_Sphere_4pi_restore_neutron = 0;
#line 49 "Samples_Incoherent.instr"
  mccPSD_Sphere_4pi_nowritefile = 0;
#line 16149 "Samples_Incoherent.c"

  SIG_MESSAGE("PSD_Sphere_4pi (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 16156 "Samples_Incoherent.c"
  rot_mul(mctr1, mcrotaCoord, mcrotaPSD_Sphere_4pi);
  rot_transpose(mcrotaTof_Sphere_4pi, mctr1);
  rot_mul(mcrotaPSD_Sphere_4pi, mctr1, mcrotrPSD_Sphere_4pi);
  mctc1 = coords_set(
#line 172 "Samples_Incoherent.instr"
    0,
#line 172 "Samples_Incoherent.instr"
    0,
#line 172 "Samples_Incoherent.instr"
    0);
#line 16167 "Samples_Incoherent.c"
  rot_transpose(mcrotaCoord, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaPSD_Sphere_4pi = coords_add(mcposaCoord, mctc2);
  mctc1 = coords_sub(mcposaTof_Sphere_4pi, mcposaPSD_Sphere_4pi);
  mcposrPSD_Sphere_4pi = rot_apply(mcrotaPSD_Sphere_4pi, mctc1);
  mcDEBUG_COMPONENT("PSD_Sphere_4pi", mcposaPSD_Sphere_4pi, mcrotaPSD_Sphere_4pi)
  mccomp_posa[16] = mcposaPSD_Sphere_4pi;
  mccomp_posr[16] = mcposrPSD_Sphere_4pi;
  mcNCounter[16]  = mcPCounter[16] = mcP2Counter[16] = 0;
  mcAbsorbProp[16]= 0;
    /* Component Coord3. */
  /* Setting parameters for component Coord3. */
  SIG_MESSAGE("Coord3 (Init:SetPar)");

  SIG_MESSAGE("Coord3 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 16187 "Samples_Incoherent.c"
  rot_mul(mctr1, mcrotaSource, mcrotaCoord3);
  rot_transpose(mcrotaPSD_Sphere_4pi, mctr1);
  rot_mul(mcrotaCoord3, mctr1, mcrotrCoord3);
  mctc1 = coords_set(
#line 175 "Samples_Incoherent.instr"
    0.0,
#line 175 "Samples_Incoherent.instr"
    0.0,
#line 175 "Samples_Incoherent.instr"
    8.3000);
#line 16198 "Samples_Incoherent.c"
  rot_transpose(mcrotaSource, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaCoord3 = coords_add(mcposaSource, mctc2);
  mctc1 = coords_sub(mcposaPSD_Sphere_4pi, mcposaCoord3);
  mcposrCoord3 = rot_apply(mcrotaCoord3, mctc1);
  mcDEBUG_COMPONENT("Coord3", mcposaCoord3, mcrotaCoord3)
  mccomp_posa[17] = mcposaCoord3;
  mccomp_posr[17] = mcposrCoord3;
  mcNCounter[17]  = mcPCounter[17] = mcP2Counter[17] = 0;
  mcAbsorbProp[17]= 0;
    /* Component Sphere_4pi. */
  /* Setting parameters for component Sphere_4pi. */
  SIG_MESSAGE("Sphere_4pi (Init:SetPar)");
#line 180 "Samples_Incoherent.instr"
  if("Events") strncpy(mccSphere_4pi_filename, "Events" ? "Events" : "", 16384); else mccSphere_4pi_filename[0]='\0';
#line 180 "Samples_Incoherent.instr"
  mccSphere_4pi_radius = 0.052;
#line 49 "Samples_Incoherent.instr"
  mccSphere_4pi_restore_neutron = 0;
#line 49 "Samples_Incoherent.instr"
  mccSphere_4pi_nowritefile = 0;
#line 16220 "Samples_Incoherent.c"

  SIG_MESSAGE("Sphere_4pi (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 16227 "Samples_Incoherent.c"
  rot_mul(mctr1, mcrotaCoord, mcrotaSphere_4pi);
  rot_transpose(mcrotaCoord3, mctr1);
  rot_mul(mcrotaSphere_4pi, mctr1, mcrotrSphere_4pi);
  mctc1 = coords_set(
#line 181 "Samples_Incoherent.instr"
    0,
#line 181 "Samples_Incoherent.instr"
    0,
#line 181 "Samples_Incoherent.instr"
    0);
#line 16238 "Samples_Incoherent.c"
  rot_transpose(mcrotaCoord, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaSphere_4pi = coords_add(mcposaCoord, mctc2);
  mctc1 = coords_sub(mcposaCoord3, mcposaSphere_4pi);
  mcposrSphere_4pi = rot_apply(mcrotaSphere_4pi, mctc1);
  mcDEBUG_COMPONENT("Sphere_4pi", mcposaSphere_4pi, mcrotaSphere_4pi)
  mccomp_posa[18] = mcposaSphere_4pi;
  mccomp_posr[18] = mcposrSphere_4pi;
  mcNCounter[18]  = mcPCounter[18] = mcP2Counter[18] = 0;
  mcAbsorbProp[18]= 0;
  /* Component initializations. */
  /* Initializations for component Source. */
  SIG_MESSAGE("Source (Init)");
#define mccompcurname  Source
#define mccompcurtype  ISIS_moderator
#define mccompcurindex 1
#define Face mccSource_Face
#define p_in mccSource_p_in
#define Tnpts mccSource_Tnpts
#define scaleSize mccSource_scaleSize
#define angleArea mccSource_angleArea
#define Nsim mccSource_Nsim
#define Ncount mccSource_Ncount
#define TS mccSource_TS
#define rtE0 mccSource_rtE0
#define rtE1 mccSource_rtE1
#define rtmodX mccSource_rtmodX
#define rtmodY mccSource_rtmodY
#define TargetStation mccSource_TargetStation
#define CurrentWeight mccSource_CurrentWeight
#define Emin mccSource_Emin
#define Emax mccSource_Emax
#define dist mccSource_dist
#define focus_xw mccSource_focus_xw
#define focus_yh mccSource_focus_yh
#define xwidth mccSource_xwidth
#define yheight mccSource_yheight
#define CAngle mccSource_CAngle
#define SAC mccSource_SAC
#define Lmin mccSource_Lmin
#define Lmax mccSource_Lmax
#define target_index mccSource_target_index
#define verbose mccSource_verbose
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
#line 16505 "Samples_Incoherent.c"
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

  /* Initializations for component Coord. */
  SIG_MESSAGE("Coord (Init)");
#define mccompcurname  Coord
#define mccompcurtype  Progress_bar
#define mccompcurindex 2
#define IntermediateCnts mccCoord_IntermediateCnts
#define StartTime mccCoord_StartTime
#define EndTime mccCoord_EndTime
#define CurrentTime mccCoord_CurrentTime
#define profile mccCoord_profile
#define percent mccCoord_percent
#define flag_save mccCoord_flag_save
#define minutes mccCoord_minutes
#line 57 "/usr/share/mcstas/2.6rc1/misc/Progress_bar.comp"
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
#line 16562 "Samples_Incoherent.c"
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

  /* Initializations for component Vsample. */
  SIG_MESSAGE("Vsample (Init)");
#define mccompcurname  Vsample
#define mccompcurtype  V_sample
#define mccompcurindex 3
#define VarsV mccVsample_VarsV
#define radius mccVsample_radius
#define thickness mccVsample_thickness
#define zdepth mccVsample_zdepth
#define Vc mccVsample_Vc
#define sigma_abs mccVsample_sigma_abs
#define sigma_inc mccVsample_sigma_inc
#define radius_i mccVsample_radius_i
#define radius_o mccVsample_radius_o
#define h mccVsample_h
#define focus_r mccVsample_focus_r
#define pack mccVsample_pack
#define frac mccVsample_frac
#define f_QE mccVsample_f_QE
#define gamma mccVsample_gamma
#define target_x mccVsample_target_x
#define target_y mccVsample_target_y
#define target_z mccVsample_target_z
#define focus_xw mccVsample_focus_xw
#define focus_yh mccVsample_focus_yh
#define focus_aw mccVsample_focus_aw
#define focus_ah mccVsample_focus_ah
#define xwidth mccVsample_xwidth
#define yheight mccVsample_yheight
#define zthick mccVsample_zthick
#define rad_sphere mccVsample_rad_sphere
#define sig_a mccVsample_sig_a
#define sig_i mccVsample_sig_i
#define V0 mccVsample_V0
#define target_index mccVsample_target_index
#define multiples mccVsample_multiples
#line 121 "/usr/share/mcstas/2.6rc1/obsolete/V_sample.comp"
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
#line 16669 "Samples_Incoherent.c"
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

  /* Initializations for component PowN. */
  SIG_MESSAGE("PowN (Init)");
#define mccompcurname  PowN
#define mccompcurtype  PowderN
#define mccompcurindex 4
#define format mccPowN_format
#define line_info mccPowN_line_info
#define columns mccPowN_columns
#define offdata mccPowN_offdata
#define reflections mccPowN_reflections
#define geometry mccPowN_geometry
#define radius mccPowN_radius
#define yheight mccPowN_yheight
#define xwidth mccPowN_xwidth
#define zdepth mccPowN_zdepth
#define thickness mccPowN_thickness
#define pack mccPowN_pack
#define Vc mccPowN_Vc
#define sigma_abs mccPowN_sigma_abs
#define sigma_inc mccPowN_sigma_inc
#define delta_d_d mccPowN_delta_d_d
#define p_inc mccPowN_p_inc
#define p_transmit mccPowN_p_transmit
#define DW mccPowN_DW
#define nb_atoms mccPowN_nb_atoms
#define d_omega mccPowN_d_omega
#define d_phi mccPowN_d_phi
#define tth_sign mccPowN_tth_sign
#define p_interact mccPowN_p_interact
#define concentric mccPowN_concentric
#define density mccPowN_density
#define weight mccPowN_weight
#define barns mccPowN_barns
#define Strain mccPowN_Strain
#define focus_flip mccPowN_focus_flip
#define target_index mccPowN_target_index
#line 559 "/usr/share/mcstas/2.6rc1/samples/PowderN.comp"
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
#line 16955 "Samples_Incoherent.c"
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

  /* Initializations for component SX. */
  SIG_MESSAGE("SX (Init)");
#define mccompcurname  SX
#define mccompcurtype  Single_crystal
#define mccompcurindex 5
#define mosaic_AB mccSX_mosaic_AB
#define hkl_info mccSX_hkl_info
#define offdata mccSX_offdata
#define reflections mccSX_reflections
#define geometry mccSX_geometry
#define xwidth mccSX_xwidth
#define yheight mccSX_yheight
#define zdepth mccSX_zdepth
#define radius mccSX_radius
#define delta_d_d mccSX_delta_d_d
#define mosaic mccSX_mosaic
#define mosaic_a mccSX_mosaic_a
#define mosaic_b mccSX_mosaic_b
#define mosaic_c mccSX_mosaic_c
#define recip_cell mccSX_recip_cell
#define barns mccSX_barns
#define ax mccSX_ax
#define ay mccSX_ay
#define az mccSX_az
#define bx mccSX_bx
#define by mccSX_by
#define bz mccSX_bz
#define cx mccSX_cx
#define cy mccSX_cy
#define cz mccSX_cz
#define p_transmit mccSX_p_transmit
#define sigma_abs mccSX_sigma_abs
#define sigma_inc mccSX_sigma_inc
#define aa mccSX_aa
#define bb mccSX_bb
#define cc mccSX_cc
#define order mccSX_order
#define RX mccSX_RX
#define RY mccSX_RY
#define powder mccSX_powder
#define PG mccSX_PG
#define deltak mccSX_deltak
#line 976 "/usr/share/mcstas/2.6rc1/samples/Single_crystal.comp"
{
  double as, bs, cs;
  int i=0;

  /* transfer input parameters */
  hkl_info.m_delta_d_d = delta_d_d;
  hkl_info.m_a  = 0;
  hkl_info.m_b  = 0;
  hkl_info.m_c  = 0;
  hkl_info.m_aa = aa;
  hkl_info.m_bb = bb;
  hkl_info.m_cc = cc;
  hkl_info.m_ax = ax;
  hkl_info.m_ay = ay;
  hkl_info.m_az = az;
  hkl_info.m_bx = bx;
  hkl_info.m_by = by;
  hkl_info.m_bz = bz;
  hkl_info.m_cx = cx;
  hkl_info.m_cy = cy;
  hkl_info.m_cz = cz;
  hkl_info.sigma_a = sigma_abs;
  hkl_info.sigma_i = sigma_inc;
  hkl_info.recip   = recip_cell;

  /* default format h,k,l,F,F2  */
  hkl_info.column_order[0]=1;
  hkl_info.column_order[1]=2;
  hkl_info.column_order[2]=3;
  hkl_info.column_order[3]=0;
  hkl_info.column_order[4]=7;
  hkl_info.kix = hkl_info.kiy = hkl_info.kiz = 0;
  hkl_info.nb_reuses = hkl_info.nb_refl = hkl_info.nb_refl_count = 0;
  hkl_info.tau_count = 0;

  /*this is necessary to allow a numerical array to be passed through as a DEFINITION parameter*/
  double mosaic_ABin[]=mosaic_AB;
  /* Read in structure factors, and do some pre-calculations. */
  if (!read_hkl_data(reflections, &hkl_info, mosaic, mosaic_a, mosaic_b, mosaic_c, mosaic_ABin)) {
    printf("Single_crystal: %s: Error: Aborting.\n", NAME_CURRENT_COMP);
    exit(0);
  }

  if (hkl_info.sigma_a<0) hkl_info.sigma_a=0;
  if (hkl_info.sigma_i<0) hkl_info.sigma_i=0;

  if (hkl_info.count)
    printf("Single_crystal: %s: Read %d reflections from file '%s'\n",
      NAME_CURRENT_COMP, hkl_info.count, reflections);
  else printf("Single_crystal: %s: Using incoherent elastic scattering only sigma=%g.\n",
      NAME_CURRENT_COMP, hkl_info.sigma_i);

  hkl_info.shape=-1; /* -1:no shape, 0:cyl, 1:box, 2:sphere, 3:any-shape  */
  if (geometry && strlen(geometry) && strcmp(geometry, "NULL") && strcmp(geometry, "0")) {
          if (off_init(geometry, xwidth, yheight, zdepth, 0, &offdata)) {
      hkl_info.shape=3;
    }
  }
  else if (xwidth && yheight && zdepth)  hkl_info.shape=1; /* box */
  else if (radius > 0 && yheight)        hkl_info.shape=0; /* cylinder */
  else if (radius > 0 && !yheight)       hkl_info.shape=2; /* sphere */

  if (hkl_info.shape < 0)
    exit(fprintf(stderr,"Single_crystal: %s: sample has invalid dimensions.\n"
                        "ERROR           Please check parameter values (xwidth, yheight, zdepth, radius).\n", NAME_CURRENT_COMP));

  printf("Single_crystal: %s: Vc=%g [Angs] sigma_abs=%g [barn] sigma_inc=%g [barn] reflections=%s\n",
      NAME_CURRENT_COMP, hkl_info.V0, hkl_info.sigma_a, hkl_info.sigma_i,
      reflections && strlen(reflections) ? reflections : "NULL");

  if (powder && PG)
    exit(fprintf(stderr,"Single_crystal: %s: powder and PG modes can not be used together!\n"
             "ERROR           Please use EITHER powder or PG mode.\n", NAME_CURRENT_COMP));

  if (powder && !(order==1)) {
    fprintf(stderr,"Single_crystal: %s: powder mode means implicit choice of no multiple scattering!\n"
            "WARNING setting order=1\n", NAME_CURRENT_COMP);
    order=1;
  }

  if (PG && !(order==1)) {
    fprintf(stderr,"Single_crystal: %s: PG mode means implicit choice of no multiple scattering!\n"
            "WARNING setting order=1\n", NAME_CURRENT_COMP);
    order=1;
  }


}
#line 17122 "Samples_Incoherent.c"
#undef deltak
#undef PG
#undef powder
#undef RY
#undef RX
#undef order
#undef cc
#undef bb
#undef aa
#undef sigma_inc
#undef sigma_abs
#undef p_transmit
#undef cz
#undef cy
#undef cx
#undef bz
#undef by
#undef bx
#undef az
#undef ay
#undef ax
#undef barns
#undef recip_cell
#undef mosaic_c
#undef mosaic_b
#undef mosaic_a
#undef mosaic
#undef delta_d_d
#undef radius
#undef zdepth
#undef yheight
#undef xwidth
#undef geometry
#undef reflections
#undef offdata
#undef hkl_info
#undef mosaic_AB
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component Sqw. */
  SIG_MESSAGE("Sqw (Init)");
#define mccompcurname  Sqw
#define mccompcurtype  Isotropic_Sqw
#define mccompcurindex 6
#define powder_format mccSqw_powder_format
#define VarSqw mccSqw_VarSqw
#define columns mccSqw_columns
#define offdata mccSqw_offdata
#define Sqw_coh mccSqw_Sqw_coh
#define Sqw_inc mccSqw_Sqw_inc
#define geometry mccSqw_geometry
#define radius mccSqw_radius
#define thickness mccSqw_thickness
#define xwidth mccSqw_xwidth
#define yheight mccSqw_yheight
#define zdepth mccSqw_zdepth
#define threshold mccSqw_threshold
#define order mccSqw_order
#define T mccSqw_T
#define verbose mccSqw_verbose
#define d_phi mccSqw_d_phi
#define concentric mccSqw_concentric
#define rho mccSqw_rho
#define sigma_abs mccSqw_sigma_abs
#define sigma_coh mccSqw_sigma_coh
#define sigma_inc mccSqw_sigma_inc
#define classical mccSqw_classical
#define powder_Dd mccSqw_powder_Dd
#define powder_DW mccSqw_powder_DW
#define powder_Vc mccSqw_powder_Vc
#define density mccSqw_density
#define weight mccSqw_weight
#define p_interact mccSqw_p_interact
#define norm mccSqw_norm
#define powder_barns mccSqw_powder_barns
#define quantum_correction mccSqw_quantum_correction
#line 2008 "/usr/share/mcstas/2.6rc1/samples/Isotropic_Sqw.comp"
{
  int i;
  /* check for parameters */
  columns = (int[])powder_format;

  VarSqw.verbose_output= verbose;
  VarSqw.shape = -1; /* -1:no shape, 0:cyl, 1:box, 2:sphere, 3:any-shape  */
  if (geometry && strlen(geometry) && strcmp(geometry, "NULL") && strcmp(geometry, "0")) {
	  if (off_init(geometry, xwidth, yheight, zdepth, 0, &offdata)) {
      VarSqw.shape=3; thickness=0; concentric=0;
    }
  }
  else if (xwidth && yheight && zdepth)  VarSqw.shape=1; /* box */
  else if (radius > 0 && yheight)        VarSqw.shape=0; /* cylinder */
  else if (radius > 0 && !yheight)       VarSqw.shape=2; /* sphere */

  if (VarSqw.shape < 0)
    exit(fprintf(stderr,"Isotropic_Sqw: %s: sample has invalid dimensions.\n"
                        "ERROR          Please check parameter values (xwidth, yheight, zdepth, radius).\n", NAME_CURRENT_COMP));



  if (thickness) {
    if (radius && (radius < fabs(thickness) )) {
      MPI_MASTER(
      fprintf(stderr,"Isotropic_Sqw: %s: hollow sample thickness is larger than its volume (sphere/cylinder).\n"
                     "WARNING        Please check parameter values. Using bulk sample (thickness=0).\n", NAME_CURRENT_COMP);
      );
      thickness=0;
    }
    else if (!radius && (xwidth < 2*fabs(thickness) || yheight < 2*fabs(thickness) || zdepth < 2*fabs(thickness))) {
      MPI_MASTER(
      fprintf(stderr,"Isotropic_Sqw: %s: hollow sample thickness is larger than its volume (box).\n"
                     "WARNING        Please check parameter values.\n", NAME_CURRENT_COMP);
      );
    }
  }
  MPI_MASTER(
  if (VarSqw.verbose_output) {
    switch (VarSqw.shape) {
      case 0: printf("Isotropic_Sqw: %s: is a %scylinder: radius=%f thickness=%f height=%f [J Comp Phys 228 (2009) 5251]\n",
              NAME_CURRENT_COMP, (thickness ? "hollow " : ""),
              radius,fabs(thickness),yheight);
              break;
      case 1: printf("Isotropic_Sqw: %s: is a %sbox: width=%f height=%f depth=%f \n",
              NAME_CURRENT_COMP, (thickness ? "hollow " : ""), xwidth,yheight,zdepth);
              break;
      case 2: printf("Isotropic_Sqw: %s: is a %ssphere: radius=%f thickness=%f\n",
              NAME_CURRENT_COMP, (thickness ? "hollow " : ""),
              radius,fabs(thickness));
              break;
      case 3: printf("Isotropic_Sqw: %s: is a volume defined from file %s\n",
              NAME_CURRENT_COMP, geometry);
    }
  }
  );

  if (concentric && !thickness) {
    MPI_MASTER(
    printf("Isotropic_Sqw: %s:Can not use concentric mode\n"
           "WARNING        on non hollow shape. Ignoring.\n",
           NAME_CURRENT_COMP);
    );
    concentric=0;
  }

  strncpy(VarSqw.compname, NAME_CURRENT_COMP, 256);
  VarSqw.T2E       =(1/11.605);   /* Kelvin to meV = 1000*KB/e */
  VarSqw.sqSE2K    = (V2K*SE2V)*(V2K*SE2V);
  VarSqw.sqw_threshold = (threshold > 0 ? threshold : 0);
  VarSqw.s_abs     = sigma_abs;
  VarSqw.s_coh     = sigma_coh;
  VarSqw.s_inc     = sigma_inc; /* s_scatt member initialized in Sqw_init */
  VarSqw.maxloop   = 100;       /* atempts to close triangle */
  VarSqw.minevents = 100;       /* minimal # of events required to get dynamical range */
  VarSqw.neutron_removed = 0;
  VarSqw.neutron_enter   = 0;
  VarSqw.neutron_pmult   = 0;
  VarSqw.neutron_exit    = 0;
  VarSqw.mat_rho       = rho;
  VarSqw.sqw_norm  = norm;
  VarSqw.mean_scatt= 0;
  VarSqw.mean_abs  = 0;
  VarSqw.psum_scatt= 0;
  VarSqw.single_coh= 0;
  VarSqw.single_inc= 0;
  VarSqw.multi     = 0;
  VarSqw.barns     = powder_barns;
  VarSqw.sqw_classical = classical;
  VarSqw.lookup_length=100;
  VarSqw.mat_weight    = weight;
  VarSqw.mat_density   = density;
  if (quantum_correction && strlen(quantum_correction)) 
    strncpy(VarSqw.Q_correction, quantum_correction, 256);
  else
    strncpy(VarSqw.Q_correction, "default", 256);

  /* PowderN compatibility members */
  VarSqw.Dd        = powder_Dd;
  VarSqw.DWfactor  = powder_DW;
  VarSqw.Temperature= T;
  for (i=0; i< 9; i++) VarSqw.column_order[i] = columns[i];
  VarSqw.column_order[8] = (VarSqw.column_order[0] >= 0 ? 0 : 2);

  /* optional ways to define rho */
  if (!VarSqw.mat_rho && powder_Vc > 0)
    VarSqw.mat_rho = 1/powder_Vc;
  /* import the data files ================================================== */
  if (!Sqw_init(&VarSqw, Sqw_coh, Sqw_inc)) {
    MPI_MASTER(
    printf("Isotropic_Sqw: %s: ERROR importing data files (Sqw_init coh=%s inc=%s).\n",NAME_CURRENT_COMP, Sqw_coh, Sqw_inc);
    );
  }
  if ( VarSqw.s_coh < 0) VarSqw.s_coh=0;
  if ( VarSqw.s_inc < 0) VarSqw.s_inc=0;
  if ( VarSqw.s_abs < 0) VarSqw.s_abs=0;
  if ((VarSqw.s_coh > 0 || VarSqw.s_inc > 0) && VarSqw.mat_rho <= 0) {
    MPI_MASTER(
    printf("Isotropic_Sqw: %s: WARNING: Null density (V_rho). Unactivating component.\n",NAME_CURRENT_COMP);
    );
    VarSqw.s_coh=VarSqw.s_inc=0;
  }
  /* 100: convert from barns to fm^2 */
  VarSqw.my_a_v  =(VarSqw.mat_rho*100*VarSqw.s_abs*2200);
  VarSqw.my_s    =(VarSqw.mat_rho*100*(VarSqw.s_coh>0 ? VarSqw.s_coh : 0
                                     +VarSqw.s_inc>0 ? VarSqw.s_inc : 0));
  MPI_MASTER(
  if ((VarSqw.s_coh > 0 || VarSqw.s_inc > 0) && !VarSqw.Temperature
   && (VarSqw.Data_coh.intensity || VarSqw.Data_inc.intensity)
   && VarSqw.verbose_output)
    printf("Isotropic_Sqw: %s: Sample temperature not defined (T=0).\n"
           "Warning        Disabling detailed balance.\n", NAME_CURRENT_COMP);
  if (VarSqw.s_coh<=0 && VarSqw.s_inc<=0) {
    printf("Isotropic_Sqw: %s: Scattering cross section is zero\n"
           "ERROR          (sigma_coh, sigma_inc).\n",NAME_CURRENT_COMP);
  }
  );
  if (d_phi) d_phi = fabs(d_phi)*DEG2RAD;

  if (d_phi > PI) d_phi = 0; /* V_scatt on 4*PI */

  if (d_phi && order != 1) {
    MPI_MASTER(
    printf("Isotropic_Sqw: %s: Focusing can only apply for single\n"
           "               scattering. Setting to order=1.\n",
           NAME_CURRENT_COMP);
    );
    order = 1;
  }

  /* request statistics */
  if (VarSqw.verbose_output > 1) {
    Sqw_diagnosis(&VarSqw, &VarSqw.Data_coh);
    Sqw_diagnosis(&VarSqw, &VarSqw.Data_inc);
  }

  for (i=0; i < 2; i++) {
    struct Sqw_Data_struct Data_sqw;
    Data_sqw =  (i == 0 ? VarSqw.Data_coh : VarSqw.Data_inc);
    Table_Free(&(Data_sqw.Sqw));
  }

/* end INITIALIZE */
}
#line 17366 "Samples_Incoherent.c"
#undef quantum_correction
#undef powder_barns
#undef norm
#undef p_interact
#undef weight
#undef density
#undef powder_Vc
#undef powder_DW
#undef powder_Dd
#undef classical
#undef sigma_inc
#undef sigma_coh
#undef sigma_abs
#undef rho
#undef concentric
#undef d_phi
#undef verbose
#undef T
#undef order
#undef threshold
#undef zdepth
#undef yheight
#undef xwidth
#undef thickness
#undef radius
#undef geometry
#undef Sqw_inc
#undef Sqw_coh
#undef offdata
#undef columns
#undef VarSqw
#undef powder_format
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component IncS. */
  SIG_MESSAGE("IncS (Init)");
#define mccompcurname  IncS
#define mccompcurtype  Incoherent
#define mccompcurindex 7
#define VarsInc mccIncS_VarsInc
#define offdata mccIncS_offdata
#define geometry mccIncS_geometry
#define radius mccIncS_radius
#define xwidth mccIncS_xwidth
#define yheight mccIncS_yheight
#define zdepth mccIncS_zdepth
#define thickness mccIncS_thickness
#define target_x mccIncS_target_x
#define target_y mccIncS_target_y
#define target_z mccIncS_target_z
#define focus_r mccIncS_focus_r
#define focus_xw mccIncS_focus_xw
#define focus_yh mccIncS_focus_yh
#define focus_aw mccIncS_focus_aw
#define focus_ah mccIncS_focus_ah
#define target_index mccIncS_target_index
#define pack mccIncS_pack
#define p_interact mccIncS_p_interact
#define f_QE mccIncS_f_QE
#define gamma mccIncS_gamma
#define Etrans mccIncS_Etrans
#define deltaE mccIncS_deltaE
#define sigma_abs mccIncS_sigma_abs
#define sigma_inc mccIncS_sigma_inc
#define Vc mccIncS_Vc
#define concentric mccIncS_concentric
#define order mccIncS_order
#line 156 "/usr/share/mcstas/2.6rc1/samples/Incoherent.comp"
{
  VarsInc.shape=-1; /* -1:no shape, 0:cyl, 1:box, 2:sphere, 3:any-shape  */
  if (geometry && strlen(geometry) && strcmp(geometry, "NULL") && strcmp(geometry, "0")) {
	  if (off_init(geometry, xwidth, yheight, zdepth, 0, &offdata)) {
      VarsInc.shape=3; thickness=0; concentric=0;
    }
  }
  else if (xwidth && yheight && zdepth)  VarsInc.shape=1; /* box */
  else if (radius > 0 &&  yheight)       VarsInc.shape=0; /* cylinder */
  else if (radius > 0 && !yheight)       VarsInc.shape=2; /* sphere */

  if (VarsInc.shape < 0)
    exit(fprintf(stderr,"Incoherent: %s: sample has invalid dimensions.\n"
                        "ERROR       Please check parameter values (xwidth, yheight, zdepth, radius).\n", NAME_CURRENT_COMP));
  if (thickness) {
    if (radius && (radius < thickness || ( yheight && (yheight < 2*thickness)))) {
      fprintf(stderr,"Incoherent: %s: hollow sample thickness is larger than its volume (sphere/cylinder).\n"
                     "WARNING     Please check parameter values. Using bulk sample (thickness=0).\n", NAME_CURRENT_COMP);
      thickness=0;
    }
    else if (!radius && (xwidth < 2*thickness || yheight < 2*thickness || zdepth < 2*thickness)) {
      fprintf(stderr,"Incoherent: %s: hollow sample thickness is larger than its volume (box).\n"
                     "WARNING     Please check parameter values. Using bulk sample (thickness=0).\n", NAME_CURRENT_COMP);
      thickness=0;
    }
  }

  if (concentric && thickness<=0) {
    printf("Incoherent: %s:Can not use concentric mode\n"
           "WARNING     on non hollow shape. Ignoring.\n",
           NAME_CURRENT_COMP);
    concentric=0;
  }

  VarsInc.sigma_a= sigma_abs;
  VarsInc.sigma_i= sigma_inc;
  VarsInc.rho    = (pack/Vc);
  VarsInc.my_s   = (VarsInc.rho * 100 * VarsInc.sigma_i);
  VarsInc.my_a_v = (VarsInc.rho * 100 * VarsInc.sigma_a);

  /* now compute target coords if a component index is supplied */
  VarsInc.tx= VarsInc.ty=VarsInc.tz=0;
  if (!target_index && !target_x && !target_y && !target_z) target_index=1;
  if (target_index)
  {
    Coords ToTarget;
    ToTarget = coords_sub(POS_A_COMP_INDEX(INDEX_CURRENT_COMP+target_index),POS_A_CURRENT_COMP);
    ToTarget = rot_apply(ROT_A_CURRENT_COMP, ToTarget);
    coords_get(ToTarget, &VarsInc.tx, &VarsInc.ty, &VarsInc.tz);
  }
  else
  { VarsInc.tx = target_x; VarsInc.ty = target_y; VarsInc.tz = target_z; }

  if (!(VarsInc.tx || VarsInc.ty || VarsInc.tz)) {
    MPI_MASTER(
    printf("Incoherent: %s: The target is not defined. Using direct beam (Z-axis).\n",
      NAME_CURRENT_COMP);
    );
    VarsInc.tz=1;
  }

  /* different ways of setting rectangular area */
  VarsInc.aw  = VarsInc.ah = 0;
  if (focus_xw) { VarsInc.xw = focus_xw; }
  if (focus_yh) { VarsInc.yh = focus_yh; }
  if (focus_aw) { VarsInc.aw = DEG2RAD*focus_aw; }
  if (focus_ah) { VarsInc.ah = DEG2RAD*focus_ah; }

  MPI_MASTER(
  printf("Incoherent: %s: Vc=%g [Angs] sigma_abs=%g [barn] sigma_inc=%g [barn]\n",
      NAME_CURRENT_COMP, Vc, VarsInc.sigma_a, VarsInc.sigma_i);
  );

}
#line 17511 "Samples_Incoherent.c"
#undef order
#undef concentric
#undef Vc
#undef sigma_inc
#undef sigma_abs
#undef deltaE
#undef Etrans
#undef gamma
#undef f_QE
#undef p_interact
#undef pack
#undef target_index
#undef focus_ah
#undef focus_aw
#undef focus_yh
#undef focus_xw
#undef focus_r
#undef target_z
#undef target_y
#undef target_x
#undef thickness
#undef zdepth
#undef yheight
#undef xwidth
#undef radius
#undef geometry
#undef offdata
#undef VarsInc
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component PowNg. */
  SIG_MESSAGE("PowNg (Init)");
#define mccompcurname  PowNg
#define mccompcurtype  PowderN
#define mccompcurindex 8
#define format mccPowNg_format
#define line_info mccPowNg_line_info
#define columns mccPowNg_columns
#define offdata mccPowNg_offdata
#define reflections mccPowNg_reflections
#define geometry mccPowNg_geometry
#define radius mccPowNg_radius
#define yheight mccPowNg_yheight
#define xwidth mccPowNg_xwidth
#define zdepth mccPowNg_zdepth
#define thickness mccPowNg_thickness
#define pack mccPowNg_pack
#define Vc mccPowNg_Vc
#define sigma_abs mccPowNg_sigma_abs
#define sigma_inc mccPowNg_sigma_inc
#define delta_d_d mccPowNg_delta_d_d
#define p_inc mccPowNg_p_inc
#define p_transmit mccPowNg_p_transmit
#define DW mccPowNg_DW
#define nb_atoms mccPowNg_nb_atoms
#define d_omega mccPowNg_d_omega
#define d_phi mccPowNg_d_phi
#define tth_sign mccPowNg_tth_sign
#define p_interact mccPowNg_p_interact
#define concentric mccPowNg_concentric
#define density mccPowNg_density
#define weight mccPowNg_weight
#define barns mccPowNg_barns
#define Strain mccPowNg_Strain
#define focus_flip mccPowNg_focus_flip
#define target_index mccPowNg_target_index
#line 559 "/usr/share/mcstas/2.6rc1/samples/PowderN.comp"
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
#line 17794 "Samples_Incoherent.c"
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

  /* Initializations for component SXg. */
  SIG_MESSAGE("SXg (Init)");
#define mccompcurname  SXg
#define mccompcurtype  Single_crystal
#define mccompcurindex 9
#define mosaic_AB mccSXg_mosaic_AB
#define hkl_info mccSXg_hkl_info
#define offdata mccSXg_offdata
#define reflections mccSXg_reflections
#define geometry mccSXg_geometry
#define xwidth mccSXg_xwidth
#define yheight mccSXg_yheight
#define zdepth mccSXg_zdepth
#define radius mccSXg_radius
#define delta_d_d mccSXg_delta_d_d
#define mosaic mccSXg_mosaic
#define mosaic_a mccSXg_mosaic_a
#define mosaic_b mccSXg_mosaic_b
#define mosaic_c mccSXg_mosaic_c
#define recip_cell mccSXg_recip_cell
#define barns mccSXg_barns
#define ax mccSXg_ax
#define ay mccSXg_ay
#define az mccSXg_az
#define bx mccSXg_bx
#define by mccSXg_by
#define bz mccSXg_bz
#define cx mccSXg_cx
#define cy mccSXg_cy
#define cz mccSXg_cz
#define p_transmit mccSXg_p_transmit
#define sigma_abs mccSXg_sigma_abs
#define sigma_inc mccSXg_sigma_inc
#define aa mccSXg_aa
#define bb mccSXg_bb
#define cc mccSXg_cc
#define order mccSXg_order
#define RX mccSXg_RX
#define RY mccSXg_RY
#define powder mccSXg_powder
#define PG mccSXg_PG
#define deltak mccSXg_deltak
#line 976 "/usr/share/mcstas/2.6rc1/samples/Single_crystal.comp"
{
  double as, bs, cs;
  int i=0;

  /* transfer input parameters */
  hkl_info.m_delta_d_d = delta_d_d;
  hkl_info.m_a  = 0;
  hkl_info.m_b  = 0;
  hkl_info.m_c  = 0;
  hkl_info.m_aa = aa;
  hkl_info.m_bb = bb;
  hkl_info.m_cc = cc;
  hkl_info.m_ax = ax;
  hkl_info.m_ay = ay;
  hkl_info.m_az = az;
  hkl_info.m_bx = bx;
  hkl_info.m_by = by;
  hkl_info.m_bz = bz;
  hkl_info.m_cx = cx;
  hkl_info.m_cy = cy;
  hkl_info.m_cz = cz;
  hkl_info.sigma_a = sigma_abs;
  hkl_info.sigma_i = sigma_inc;
  hkl_info.recip   = recip_cell;

  /* default format h,k,l,F,F2  */
  hkl_info.column_order[0]=1;
  hkl_info.column_order[1]=2;
  hkl_info.column_order[2]=3;
  hkl_info.column_order[3]=0;
  hkl_info.column_order[4]=7;
  hkl_info.kix = hkl_info.kiy = hkl_info.kiz = 0;
  hkl_info.nb_reuses = hkl_info.nb_refl = hkl_info.nb_refl_count = 0;
  hkl_info.tau_count = 0;

  /*this is necessary to allow a numerical array to be passed through as a DEFINITION parameter*/
  double mosaic_ABin[]=mosaic_AB;
  /* Read in structure factors, and do some pre-calculations. */
  if (!read_hkl_data(reflections, &hkl_info, mosaic, mosaic_a, mosaic_b, mosaic_c, mosaic_ABin)) {
    printf("Single_crystal: %s: Error: Aborting.\n", NAME_CURRENT_COMP);
    exit(0);
  }

  if (hkl_info.sigma_a<0) hkl_info.sigma_a=0;
  if (hkl_info.sigma_i<0) hkl_info.sigma_i=0;

  if (hkl_info.count)
    printf("Single_crystal: %s: Read %d reflections from file '%s'\n",
      NAME_CURRENT_COMP, hkl_info.count, reflections);
  else printf("Single_crystal: %s: Using incoherent elastic scattering only sigma=%g.\n",
      NAME_CURRENT_COMP, hkl_info.sigma_i);

  hkl_info.shape=-1; /* -1:no shape, 0:cyl, 1:box, 2:sphere, 3:any-shape  */
  if (geometry && strlen(geometry) && strcmp(geometry, "NULL") && strcmp(geometry, "0")) {
          if (off_init(geometry, xwidth, yheight, zdepth, 0, &offdata)) {
      hkl_info.shape=3;
    }
  }
  else if (xwidth && yheight && zdepth)  hkl_info.shape=1; /* box */
  else if (radius > 0 && yheight)        hkl_info.shape=0; /* cylinder */
  else if (radius > 0 && !yheight)       hkl_info.shape=2; /* sphere */

  if (hkl_info.shape < 0)
    exit(fprintf(stderr,"Single_crystal: %s: sample has invalid dimensions.\n"
                        "ERROR           Please check parameter values (xwidth, yheight, zdepth, radius).\n", NAME_CURRENT_COMP));

  printf("Single_crystal: %s: Vc=%g [Angs] sigma_abs=%g [barn] sigma_inc=%g [barn] reflections=%s\n",
      NAME_CURRENT_COMP, hkl_info.V0, hkl_info.sigma_a, hkl_info.sigma_i,
      reflections && strlen(reflections) ? reflections : "NULL");

  if (powder && PG)
    exit(fprintf(stderr,"Single_crystal: %s: powder and PG modes can not be used together!\n"
             "ERROR           Please use EITHER powder or PG mode.\n", NAME_CURRENT_COMP));

  if (powder && !(order==1)) {
    fprintf(stderr,"Single_crystal: %s: powder mode means implicit choice of no multiple scattering!\n"
            "WARNING setting order=1\n", NAME_CURRENT_COMP);
    order=1;
  }

  if (PG && !(order==1)) {
    fprintf(stderr,"Single_crystal: %s: PG mode means implicit choice of no multiple scattering!\n"
            "WARNING setting order=1\n", NAME_CURRENT_COMP);
    order=1;
  }


}
#line 17961 "Samples_Incoherent.c"
#undef deltak
#undef PG
#undef powder
#undef RY
#undef RX
#undef order
#undef cc
#undef bb
#undef aa
#undef sigma_inc
#undef sigma_abs
#undef p_transmit
#undef cz
#undef cy
#undef cx
#undef bz
#undef by
#undef bx
#undef az
#undef ay
#undef ax
#undef barns
#undef recip_cell
#undef mosaic_c
#undef mosaic_b
#undef mosaic_a
#undef mosaic
#undef delta_d_d
#undef radius
#undef zdepth
#undef yheight
#undef xwidth
#undef geometry
#undef reflections
#undef offdata
#undef hkl_info
#undef mosaic_AB
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component Sqwg. */
  SIG_MESSAGE("Sqwg (Init)");
#define mccompcurname  Sqwg
#define mccompcurtype  Isotropic_Sqw
#define mccompcurindex 10
#define powder_format mccSqwg_powder_format
#define VarSqw mccSqwg_VarSqw
#define columns mccSqwg_columns
#define offdata mccSqwg_offdata
#define Sqw_coh mccSqwg_Sqw_coh
#define Sqw_inc mccSqwg_Sqw_inc
#define geometry mccSqwg_geometry
#define radius mccSqwg_radius
#define thickness mccSqwg_thickness
#define xwidth mccSqwg_xwidth
#define yheight mccSqwg_yheight
#define zdepth mccSqwg_zdepth
#define threshold mccSqwg_threshold
#define order mccSqwg_order
#define T mccSqwg_T
#define verbose mccSqwg_verbose
#define d_phi mccSqwg_d_phi
#define concentric mccSqwg_concentric
#define rho mccSqwg_rho
#define sigma_abs mccSqwg_sigma_abs
#define sigma_coh mccSqwg_sigma_coh
#define sigma_inc mccSqwg_sigma_inc
#define classical mccSqwg_classical
#define powder_Dd mccSqwg_powder_Dd
#define powder_DW mccSqwg_powder_DW
#define powder_Vc mccSqwg_powder_Vc
#define density mccSqwg_density
#define weight mccSqwg_weight
#define p_interact mccSqwg_p_interact
#define norm mccSqwg_norm
#define powder_barns mccSqwg_powder_barns
#define quantum_correction mccSqwg_quantum_correction
#line 2008 "/usr/share/mcstas/2.6rc1/samples/Isotropic_Sqw.comp"
{
  int i;
  /* check for parameters */
  columns = (int[])powder_format;

  VarSqw.verbose_output= verbose;
  VarSqw.shape = -1; /* -1:no shape, 0:cyl, 1:box, 2:sphere, 3:any-shape  */
  if (geometry && strlen(geometry) && strcmp(geometry, "NULL") && strcmp(geometry, "0")) {
	  if (off_init(geometry, xwidth, yheight, zdepth, 0, &offdata)) {
      VarSqw.shape=3; thickness=0; concentric=0;
    }
  }
  else if (xwidth && yheight && zdepth)  VarSqw.shape=1; /* box */
  else if (radius > 0 && yheight)        VarSqw.shape=0; /* cylinder */
  else if (radius > 0 && !yheight)       VarSqw.shape=2; /* sphere */

  if (VarSqw.shape < 0)
    exit(fprintf(stderr,"Isotropic_Sqw: %s: sample has invalid dimensions.\n"
                        "ERROR          Please check parameter values (xwidth, yheight, zdepth, radius).\n", NAME_CURRENT_COMP));



  if (thickness) {
    if (radius && (radius < fabs(thickness) )) {
      MPI_MASTER(
      fprintf(stderr,"Isotropic_Sqw: %s: hollow sample thickness is larger than its volume (sphere/cylinder).\n"
                     "WARNING        Please check parameter values. Using bulk sample (thickness=0).\n", NAME_CURRENT_COMP);
      );
      thickness=0;
    }
    else if (!radius && (xwidth < 2*fabs(thickness) || yheight < 2*fabs(thickness) || zdepth < 2*fabs(thickness))) {
      MPI_MASTER(
      fprintf(stderr,"Isotropic_Sqw: %s: hollow sample thickness is larger than its volume (box).\n"
                     "WARNING        Please check parameter values.\n", NAME_CURRENT_COMP);
      );
    }
  }
  MPI_MASTER(
  if (VarSqw.verbose_output) {
    switch (VarSqw.shape) {
      case 0: printf("Isotropic_Sqw: %s: is a %scylinder: radius=%f thickness=%f height=%f [J Comp Phys 228 (2009) 5251]\n",
              NAME_CURRENT_COMP, (thickness ? "hollow " : ""),
              radius,fabs(thickness),yheight);
              break;
      case 1: printf("Isotropic_Sqw: %s: is a %sbox: width=%f height=%f depth=%f \n",
              NAME_CURRENT_COMP, (thickness ? "hollow " : ""), xwidth,yheight,zdepth);
              break;
      case 2: printf("Isotropic_Sqw: %s: is a %ssphere: radius=%f thickness=%f\n",
              NAME_CURRENT_COMP, (thickness ? "hollow " : ""),
              radius,fabs(thickness));
              break;
      case 3: printf("Isotropic_Sqw: %s: is a volume defined from file %s\n",
              NAME_CURRENT_COMP, geometry);
    }
  }
  );

  if (concentric && !thickness) {
    MPI_MASTER(
    printf("Isotropic_Sqw: %s:Can not use concentric mode\n"
           "WARNING        on non hollow shape. Ignoring.\n",
           NAME_CURRENT_COMP);
    );
    concentric=0;
  }

  strncpy(VarSqw.compname, NAME_CURRENT_COMP, 256);
  VarSqw.T2E       =(1/11.605);   /* Kelvin to meV = 1000*KB/e */
  VarSqw.sqSE2K    = (V2K*SE2V)*(V2K*SE2V);
  VarSqw.sqw_threshold = (threshold > 0 ? threshold : 0);
  VarSqw.s_abs     = sigma_abs;
  VarSqw.s_coh     = sigma_coh;
  VarSqw.s_inc     = sigma_inc; /* s_scatt member initialized in Sqw_init */
  VarSqw.maxloop   = 100;       /* atempts to close triangle */
  VarSqw.minevents = 100;       /* minimal # of events required to get dynamical range */
  VarSqw.neutron_removed = 0;
  VarSqw.neutron_enter   = 0;
  VarSqw.neutron_pmult   = 0;
  VarSqw.neutron_exit    = 0;
  VarSqw.mat_rho       = rho;
  VarSqw.sqw_norm  = norm;
  VarSqw.mean_scatt= 0;
  VarSqw.mean_abs  = 0;
  VarSqw.psum_scatt= 0;
  VarSqw.single_coh= 0;
  VarSqw.single_inc= 0;
  VarSqw.multi     = 0;
  VarSqw.barns     = powder_barns;
  VarSqw.sqw_classical = classical;
  VarSqw.lookup_length=100;
  VarSqw.mat_weight    = weight;
  VarSqw.mat_density   = density;
  if (quantum_correction && strlen(quantum_correction)) 
    strncpy(VarSqw.Q_correction, quantum_correction, 256);
  else
    strncpy(VarSqw.Q_correction, "default", 256);

  /* PowderN compatibility members */
  VarSqw.Dd        = powder_Dd;
  VarSqw.DWfactor  = powder_DW;
  VarSqw.Temperature= T;
  for (i=0; i< 9; i++) VarSqw.column_order[i] = columns[i];
  VarSqw.column_order[8] = (VarSqw.column_order[0] >= 0 ? 0 : 2);

  /* optional ways to define rho */
  if (!VarSqw.mat_rho && powder_Vc > 0)
    VarSqw.mat_rho = 1/powder_Vc;
  /* import the data files ================================================== */
  if (!Sqw_init(&VarSqw, Sqw_coh, Sqw_inc)) {
    MPI_MASTER(
    printf("Isotropic_Sqw: %s: ERROR importing data files (Sqw_init coh=%s inc=%s).\n",NAME_CURRENT_COMP, Sqw_coh, Sqw_inc);
    );
  }
  if ( VarSqw.s_coh < 0) VarSqw.s_coh=0;
  if ( VarSqw.s_inc < 0) VarSqw.s_inc=0;
  if ( VarSqw.s_abs < 0) VarSqw.s_abs=0;
  if ((VarSqw.s_coh > 0 || VarSqw.s_inc > 0) && VarSqw.mat_rho <= 0) {
    MPI_MASTER(
    printf("Isotropic_Sqw: %s: WARNING: Null density (V_rho). Unactivating component.\n",NAME_CURRENT_COMP);
    );
    VarSqw.s_coh=VarSqw.s_inc=0;
  }
  /* 100: convert from barns to fm^2 */
  VarSqw.my_a_v  =(VarSqw.mat_rho*100*VarSqw.s_abs*2200);
  VarSqw.my_s    =(VarSqw.mat_rho*100*(VarSqw.s_coh>0 ? VarSqw.s_coh : 0
                                     +VarSqw.s_inc>0 ? VarSqw.s_inc : 0));
  MPI_MASTER(
  if ((VarSqw.s_coh > 0 || VarSqw.s_inc > 0) && !VarSqw.Temperature
   && (VarSqw.Data_coh.intensity || VarSqw.Data_inc.intensity)
   && VarSqw.verbose_output)
    printf("Isotropic_Sqw: %s: Sample temperature not defined (T=0).\n"
           "Warning        Disabling detailed balance.\n", NAME_CURRENT_COMP);
  if (VarSqw.s_coh<=0 && VarSqw.s_inc<=0) {
    printf("Isotropic_Sqw: %s: Scattering cross section is zero\n"
           "ERROR          (sigma_coh, sigma_inc).\n",NAME_CURRENT_COMP);
  }
  );
  if (d_phi) d_phi = fabs(d_phi)*DEG2RAD;

  if (d_phi > PI) d_phi = 0; /* V_scatt on 4*PI */

  if (d_phi && order != 1) {
    MPI_MASTER(
    printf("Isotropic_Sqw: %s: Focusing can only apply for single\n"
           "               scattering. Setting to order=1.\n",
           NAME_CURRENT_COMP);
    );
    order = 1;
  }

  /* request statistics */
  if (VarSqw.verbose_output > 1) {
    Sqw_diagnosis(&VarSqw, &VarSqw.Data_coh);
    Sqw_diagnosis(&VarSqw, &VarSqw.Data_inc);
  }

  for (i=0; i < 2; i++) {
    struct Sqw_Data_struct Data_sqw;
    Data_sqw =  (i == 0 ? VarSqw.Data_coh : VarSqw.Data_inc);
    Table_Free(&(Data_sqw.Sqw));
  }

/* end INITIALIZE */
}
#line 18205 "Samples_Incoherent.c"
#undef quantum_correction
#undef powder_barns
#undef norm
#undef p_interact
#undef weight
#undef density
#undef powder_Vc
#undef powder_DW
#undef powder_Dd
#undef classical
#undef sigma_inc
#undef sigma_coh
#undef sigma_abs
#undef rho
#undef concentric
#undef d_phi
#undef verbose
#undef T
#undef order
#undef threshold
#undef zdepth
#undef yheight
#undef xwidth
#undef thickness
#undef radius
#undef geometry
#undef Sqw_inc
#undef Sqw_coh
#undef offdata
#undef columns
#undef VarSqw
#undef powder_format
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component IncSg. */
  SIG_MESSAGE("IncSg (Init)");
#define mccompcurname  IncSg
#define mccompcurtype  Incoherent
#define mccompcurindex 11
#define VarsInc mccIncSg_VarsInc
#define offdata mccIncSg_offdata
#define geometry mccIncSg_geometry
#define radius mccIncSg_radius
#define xwidth mccIncSg_xwidth
#define yheight mccIncSg_yheight
#define zdepth mccIncSg_zdepth
#define thickness mccIncSg_thickness
#define target_x mccIncSg_target_x
#define target_y mccIncSg_target_y
#define target_z mccIncSg_target_z
#define focus_r mccIncSg_focus_r
#define focus_xw mccIncSg_focus_xw
#define focus_yh mccIncSg_focus_yh
#define focus_aw mccIncSg_focus_aw
#define focus_ah mccIncSg_focus_ah
#define target_index mccIncSg_target_index
#define pack mccIncSg_pack
#define p_interact mccIncSg_p_interact
#define f_QE mccIncSg_f_QE
#define gamma mccIncSg_gamma
#define Etrans mccIncSg_Etrans
#define deltaE mccIncSg_deltaE
#define sigma_abs mccIncSg_sigma_abs
#define sigma_inc mccIncSg_sigma_inc
#define Vc mccIncSg_Vc
#define concentric mccIncSg_concentric
#define order mccIncSg_order
#line 156 "/usr/share/mcstas/2.6rc1/samples/Incoherent.comp"
{
  VarsInc.shape=-1; /* -1:no shape, 0:cyl, 1:box, 2:sphere, 3:any-shape  */
  if (geometry && strlen(geometry) && strcmp(geometry, "NULL") && strcmp(geometry, "0")) {
	  if (off_init(geometry, xwidth, yheight, zdepth, 0, &offdata)) {
      VarsInc.shape=3; thickness=0; concentric=0;
    }
  }
  else if (xwidth && yheight && zdepth)  VarsInc.shape=1; /* box */
  else if (radius > 0 &&  yheight)       VarsInc.shape=0; /* cylinder */
  else if (radius > 0 && !yheight)       VarsInc.shape=2; /* sphere */

  if (VarsInc.shape < 0)
    exit(fprintf(stderr,"Incoherent: %s: sample has invalid dimensions.\n"
                        "ERROR       Please check parameter values (xwidth, yheight, zdepth, radius).\n", NAME_CURRENT_COMP));
  if (thickness) {
    if (radius && (radius < thickness || ( yheight && (yheight < 2*thickness)))) {
      fprintf(stderr,"Incoherent: %s: hollow sample thickness is larger than its volume (sphere/cylinder).\n"
                     "WARNING     Please check parameter values. Using bulk sample (thickness=0).\n", NAME_CURRENT_COMP);
      thickness=0;
    }
    else if (!radius && (xwidth < 2*thickness || yheight < 2*thickness || zdepth < 2*thickness)) {
      fprintf(stderr,"Incoherent: %s: hollow sample thickness is larger than its volume (box).\n"
                     "WARNING     Please check parameter values. Using bulk sample (thickness=0).\n", NAME_CURRENT_COMP);
      thickness=0;
    }
  }

  if (concentric && thickness<=0) {
    printf("Incoherent: %s:Can not use concentric mode\n"
           "WARNING     on non hollow shape. Ignoring.\n",
           NAME_CURRENT_COMP);
    concentric=0;
  }

  VarsInc.sigma_a= sigma_abs;
  VarsInc.sigma_i= sigma_inc;
  VarsInc.rho    = (pack/Vc);
  VarsInc.my_s   = (VarsInc.rho * 100 * VarsInc.sigma_i);
  VarsInc.my_a_v = (VarsInc.rho * 100 * VarsInc.sigma_a);

  /* now compute target coords if a component index is supplied */
  VarsInc.tx= VarsInc.ty=VarsInc.tz=0;
  if (!target_index && !target_x && !target_y && !target_z) target_index=1;
  if (target_index)
  {
    Coords ToTarget;
    ToTarget = coords_sub(POS_A_COMP_INDEX(INDEX_CURRENT_COMP+target_index),POS_A_CURRENT_COMP);
    ToTarget = rot_apply(ROT_A_CURRENT_COMP, ToTarget);
    coords_get(ToTarget, &VarsInc.tx, &VarsInc.ty, &VarsInc.tz);
  }
  else
  { VarsInc.tx = target_x; VarsInc.ty = target_y; VarsInc.tz = target_z; }

  if (!(VarsInc.tx || VarsInc.ty || VarsInc.tz)) {
    MPI_MASTER(
    printf("Incoherent: %s: The target is not defined. Using direct beam (Z-axis).\n",
      NAME_CURRENT_COMP);
    );
    VarsInc.tz=1;
  }

  /* different ways of setting rectangular area */
  VarsInc.aw  = VarsInc.ah = 0;
  if (focus_xw) { VarsInc.xw = focus_xw; }
  if (focus_yh) { VarsInc.yh = focus_yh; }
  if (focus_aw) { VarsInc.aw = DEG2RAD*focus_aw; }
  if (focus_ah) { VarsInc.ah = DEG2RAD*focus_ah; }

  MPI_MASTER(
  printf("Incoherent: %s: Vc=%g [Angs] sigma_abs=%g [barn] sigma_inc=%g [barn]\n",
      NAME_CURRENT_COMP, Vc, VarsInc.sigma_a, VarsInc.sigma_i);
  );

}
#line 18350 "Samples_Incoherent.c"
#undef order
#undef concentric
#undef Vc
#undef sigma_inc
#undef sigma_abs
#undef deltaE
#undef Etrans
#undef gamma
#undef f_QE
#undef p_interact
#undef pack
#undef target_index
#undef focus_ah
#undef focus_aw
#undef focus_yh
#undef focus_xw
#undef focus_r
#undef target_z
#undef target_y
#undef target_x
#undef thickness
#undef zdepth
#undef yheight
#undef xwidth
#undef radius
#undef geometry
#undef offdata
#undef VarsInc
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component SqwV. */
  SIG_MESSAGE("SqwV (Init)");
#define mccompcurname  SqwV
#define mccompcurtype  Isotropic_Sqw
#define mccompcurindex 12
#define powder_format mccSqwV_powder_format
#define VarSqw mccSqwV_VarSqw
#define columns mccSqwV_columns
#define offdata mccSqwV_offdata
#define Sqw_coh mccSqwV_Sqw_coh
#define Sqw_inc mccSqwV_Sqw_inc
#define geometry mccSqwV_geometry
#define radius mccSqwV_radius
#define thickness mccSqwV_thickness
#define xwidth mccSqwV_xwidth
#define yheight mccSqwV_yheight
#define zdepth mccSqwV_zdepth
#define threshold mccSqwV_threshold
#define order mccSqwV_order
#define T mccSqwV_T
#define verbose mccSqwV_verbose
#define d_phi mccSqwV_d_phi
#define concentric mccSqwV_concentric
#define rho mccSqwV_rho
#define sigma_abs mccSqwV_sigma_abs
#define sigma_coh mccSqwV_sigma_coh
#define sigma_inc mccSqwV_sigma_inc
#define classical mccSqwV_classical
#define powder_Dd mccSqwV_powder_Dd
#define powder_DW mccSqwV_powder_DW
#define powder_Vc mccSqwV_powder_Vc
#define density mccSqwV_density
#define weight mccSqwV_weight
#define p_interact mccSqwV_p_interact
#define norm mccSqwV_norm
#define powder_barns mccSqwV_powder_barns
#define quantum_correction mccSqwV_quantum_correction
#line 2008 "/usr/share/mcstas/2.6rc1/samples/Isotropic_Sqw.comp"
{
  int i;
  /* check for parameters */
  columns = (int[])powder_format;

  VarSqw.verbose_output= verbose;
  VarSqw.shape = -1; /* -1:no shape, 0:cyl, 1:box, 2:sphere, 3:any-shape  */
  if (geometry && strlen(geometry) && strcmp(geometry, "NULL") && strcmp(geometry, "0")) {
	  if (off_init(geometry, xwidth, yheight, zdepth, 0, &offdata)) {
      VarSqw.shape=3; thickness=0; concentric=0;
    }
  }
  else if (xwidth && yheight && zdepth)  VarSqw.shape=1; /* box */
  else if (radius > 0 && yheight)        VarSqw.shape=0; /* cylinder */
  else if (radius > 0 && !yheight)       VarSqw.shape=2; /* sphere */

  if (VarSqw.shape < 0)
    exit(fprintf(stderr,"Isotropic_Sqw: %s: sample has invalid dimensions.\n"
                        "ERROR          Please check parameter values (xwidth, yheight, zdepth, radius).\n", NAME_CURRENT_COMP));



  if (thickness) {
    if (radius && (radius < fabs(thickness) )) {
      MPI_MASTER(
      fprintf(stderr,"Isotropic_Sqw: %s: hollow sample thickness is larger than its volume (sphere/cylinder).\n"
                     "WARNING        Please check parameter values. Using bulk sample (thickness=0).\n", NAME_CURRENT_COMP);
      );
      thickness=0;
    }
    else if (!radius && (xwidth < 2*fabs(thickness) || yheight < 2*fabs(thickness) || zdepth < 2*fabs(thickness))) {
      MPI_MASTER(
      fprintf(stderr,"Isotropic_Sqw: %s: hollow sample thickness is larger than its volume (box).\n"
                     "WARNING        Please check parameter values.\n", NAME_CURRENT_COMP);
      );
    }
  }
  MPI_MASTER(
  if (VarSqw.verbose_output) {
    switch (VarSqw.shape) {
      case 0: printf("Isotropic_Sqw: %s: is a %scylinder: radius=%f thickness=%f height=%f [J Comp Phys 228 (2009) 5251]\n",
              NAME_CURRENT_COMP, (thickness ? "hollow " : ""),
              radius,fabs(thickness),yheight);
              break;
      case 1: printf("Isotropic_Sqw: %s: is a %sbox: width=%f height=%f depth=%f \n",
              NAME_CURRENT_COMP, (thickness ? "hollow " : ""), xwidth,yheight,zdepth);
              break;
      case 2: printf("Isotropic_Sqw: %s: is a %ssphere: radius=%f thickness=%f\n",
              NAME_CURRENT_COMP, (thickness ? "hollow " : ""),
              radius,fabs(thickness));
              break;
      case 3: printf("Isotropic_Sqw: %s: is a volume defined from file %s\n",
              NAME_CURRENT_COMP, geometry);
    }
  }
  );

  if (concentric && !thickness) {
    MPI_MASTER(
    printf("Isotropic_Sqw: %s:Can not use concentric mode\n"
           "WARNING        on non hollow shape. Ignoring.\n",
           NAME_CURRENT_COMP);
    );
    concentric=0;
  }

  strncpy(VarSqw.compname, NAME_CURRENT_COMP, 256);
  VarSqw.T2E       =(1/11.605);   /* Kelvin to meV = 1000*KB/e */
  VarSqw.sqSE2K    = (V2K*SE2V)*(V2K*SE2V);
  VarSqw.sqw_threshold = (threshold > 0 ? threshold : 0);
  VarSqw.s_abs     = sigma_abs;
  VarSqw.s_coh     = sigma_coh;
  VarSqw.s_inc     = sigma_inc; /* s_scatt member initialized in Sqw_init */
  VarSqw.maxloop   = 100;       /* atempts to close triangle */
  VarSqw.minevents = 100;       /* minimal # of events required to get dynamical range */
  VarSqw.neutron_removed = 0;
  VarSqw.neutron_enter   = 0;
  VarSqw.neutron_pmult   = 0;
  VarSqw.neutron_exit    = 0;
  VarSqw.mat_rho       = rho;
  VarSqw.sqw_norm  = norm;
  VarSqw.mean_scatt= 0;
  VarSqw.mean_abs  = 0;
  VarSqw.psum_scatt= 0;
  VarSqw.single_coh= 0;
  VarSqw.single_inc= 0;
  VarSqw.multi     = 0;
  VarSqw.barns     = powder_barns;
  VarSqw.sqw_classical = classical;
  VarSqw.lookup_length=100;
  VarSqw.mat_weight    = weight;
  VarSqw.mat_density   = density;
  if (quantum_correction && strlen(quantum_correction)) 
    strncpy(VarSqw.Q_correction, quantum_correction, 256);
  else
    strncpy(VarSqw.Q_correction, "default", 256);

  /* PowderN compatibility members */
  VarSqw.Dd        = powder_Dd;
  VarSqw.DWfactor  = powder_DW;
  VarSqw.Temperature= T;
  for (i=0; i< 9; i++) VarSqw.column_order[i] = columns[i];
  VarSqw.column_order[8] = (VarSqw.column_order[0] >= 0 ? 0 : 2);

  /* optional ways to define rho */
  if (!VarSqw.mat_rho && powder_Vc > 0)
    VarSqw.mat_rho = 1/powder_Vc;
  /* import the data files ================================================== */
  if (!Sqw_init(&VarSqw, Sqw_coh, Sqw_inc)) {
    MPI_MASTER(
    printf("Isotropic_Sqw: %s: ERROR importing data files (Sqw_init coh=%s inc=%s).\n",NAME_CURRENT_COMP, Sqw_coh, Sqw_inc);
    );
  }
  if ( VarSqw.s_coh < 0) VarSqw.s_coh=0;
  if ( VarSqw.s_inc < 0) VarSqw.s_inc=0;
  if ( VarSqw.s_abs < 0) VarSqw.s_abs=0;
  if ((VarSqw.s_coh > 0 || VarSqw.s_inc > 0) && VarSqw.mat_rho <= 0) {
    MPI_MASTER(
    printf("Isotropic_Sqw: %s: WARNING: Null density (V_rho). Unactivating component.\n",NAME_CURRENT_COMP);
    );
    VarSqw.s_coh=VarSqw.s_inc=0;
  }
  /* 100: convert from barns to fm^2 */
  VarSqw.my_a_v  =(VarSqw.mat_rho*100*VarSqw.s_abs*2200);
  VarSqw.my_s    =(VarSqw.mat_rho*100*(VarSqw.s_coh>0 ? VarSqw.s_coh : 0
                                     +VarSqw.s_inc>0 ? VarSqw.s_inc : 0));
  MPI_MASTER(
  if ((VarSqw.s_coh > 0 || VarSqw.s_inc > 0) && !VarSqw.Temperature
   && (VarSqw.Data_coh.intensity || VarSqw.Data_inc.intensity)
   && VarSqw.verbose_output)
    printf("Isotropic_Sqw: %s: Sample temperature not defined (T=0).\n"
           "Warning        Disabling detailed balance.\n", NAME_CURRENT_COMP);
  if (VarSqw.s_coh<=0 && VarSqw.s_inc<=0) {
    printf("Isotropic_Sqw: %s: Scattering cross section is zero\n"
           "ERROR          (sigma_coh, sigma_inc).\n",NAME_CURRENT_COMP);
  }
  );
  if (d_phi) d_phi = fabs(d_phi)*DEG2RAD;

  if (d_phi > PI) d_phi = 0; /* V_scatt on 4*PI */

  if (d_phi && order != 1) {
    MPI_MASTER(
    printf("Isotropic_Sqw: %s: Focusing can only apply for single\n"
           "               scattering. Setting to order=1.\n",
           NAME_CURRENT_COMP);
    );
    order = 1;
  }

  /* request statistics */
  if (VarSqw.verbose_output > 1) {
    Sqw_diagnosis(&VarSqw, &VarSqw.Data_coh);
    Sqw_diagnosis(&VarSqw, &VarSqw.Data_inc);
  }

  for (i=0; i < 2; i++) {
    struct Sqw_Data_struct Data_sqw;
    Data_sqw =  (i == 0 ? VarSqw.Data_coh : VarSqw.Data_inc);
    Table_Free(&(Data_sqw.Sqw));
  }

/* end INITIALIZE */
}
#line 18585 "Samples_Incoherent.c"
#undef quantum_correction
#undef powder_barns
#undef norm
#undef p_interact
#undef weight
#undef density
#undef powder_Vc
#undef powder_DW
#undef powder_Dd
#undef classical
#undef sigma_inc
#undef sigma_coh
#undef sigma_abs
#undef rho
#undef concentric
#undef d_phi
#undef verbose
#undef T
#undef order
#undef threshold
#undef zdepth
#undef yheight
#undef xwidth
#undef thickness
#undef radius
#undef geometry
#undef Sqw_inc
#undef Sqw_coh
#undef offdata
#undef columns
#undef VarSqw
#undef powder_format
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component Dirbeam. */
  SIG_MESSAGE("Dirbeam (Init)");
#define mccompcurname  Dirbeam
#define mccompcurtype  PSD_monitor
#define mccompcurindex 13
#define PSD_N mccDirbeam_PSD_N
#define PSD_p mccDirbeam_PSD_p
#define PSD_p2 mccDirbeam_PSD_p2
#define nx mccDirbeam_nx
#define ny mccDirbeam_ny
#define filename mccDirbeam_filename
#define xmin mccDirbeam_xmin
#define xmax mccDirbeam_xmax
#define ymin mccDirbeam_ymin
#define ymax mccDirbeam_ymax
#define xwidth mccDirbeam_xwidth
#define yheight mccDirbeam_yheight
#define restore_neutron mccDirbeam_restore_neutron
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
#line 18665 "Samples_Incoherent.c"
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

  /* Initializations for component Stop. */
  SIG_MESSAGE("Stop (Init)");
#define mccompcurname  Stop
#define mccompcurtype  Beamstop
#define mccompcurindex 14
#define xmin mccStop_xmin
#define xmax mccStop_xmax
#define ymin mccStop_ymin
#define ymax mccStop_ymax
#define xwidth mccStop_xwidth
#define yheight mccStop_yheight
#define radius mccStop_radius
#line 50 "/usr/share/mcstas/2.6rc1/optics/Beamstop.comp"
{
if (xwidth  > 0) { xmax = xwidth/2;  xmin = -xmax; }
  if (yheight > 0) { ymax = yheight/2; ymin = -ymax; }

  if (xmin == 0 && xmax == 0 && ymin == 0 & ymax == 0 && radius == 0)
  { fprintf(stderr,"Beamstop: %s: Error: give geometry\n", NAME_CURRENT_COMP); exit(-1); }
}
#line 18703 "Samples_Incoherent.c"
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

  /* Initializations for component Tof_Sphere_4pi. */
  SIG_MESSAGE("Tof_Sphere_4pi (Init)");
#define mccompcurname  Tof_Sphere_4pi
#define mccompcurtype  Monitor_nD
#define mccompcurindex 15
#define user1 mccTof_Sphere_4pi_user1
#define user2 mccTof_Sphere_4pi_user2
#define user3 mccTof_Sphere_4pi_user3
#define DEFS mccTof_Sphere_4pi_DEFS
#define Vars mccTof_Sphere_4pi_Vars
#define detector mccTof_Sphere_4pi_detector
#define offdata mccTof_Sphere_4pi_offdata
#define xwidth mccTof_Sphere_4pi_xwidth
#define yheight mccTof_Sphere_4pi_yheight
#define zdepth mccTof_Sphere_4pi_zdepth
#define xmin mccTof_Sphere_4pi_xmin
#define xmax mccTof_Sphere_4pi_xmax
#define ymin mccTof_Sphere_4pi_ymin
#define ymax mccTof_Sphere_4pi_ymax
#define zmin mccTof_Sphere_4pi_zmin
#define zmax mccTof_Sphere_4pi_zmax
#define bins mccTof_Sphere_4pi_bins
#define min mccTof_Sphere_4pi_min
#define max mccTof_Sphere_4pi_max
#define restore_neutron mccTof_Sphere_4pi_restore_neutron
#define radius mccTof_Sphere_4pi_radius
#define options mccTof_Sphere_4pi_options
#define filename mccTof_Sphere_4pi_filename
#define geometry mccTof_Sphere_4pi_geometry
#define username1 mccTof_Sphere_4pi_username1
#define username2 mccTof_Sphere_4pi_username2
#define username3 mccTof_Sphere_4pi_username3
#define nowritefile mccTof_Sphere_4pi_nowritefile
#line 232 "/usr/share/mcstas/2.6rc1/monitors/Monitor_nD.comp"
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
#line 18828 "Samples_Incoherent.c"
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

  /* Initializations for component PSD_Sphere_4pi. */
  SIG_MESSAGE("PSD_Sphere_4pi (Init)");
#define mccompcurname  PSD_Sphere_4pi
#define mccompcurtype  PSD_monitor_4PI
#define mccompcurindex 16
#define nx mccPSD_Sphere_4pi_nx
#define ny mccPSD_Sphere_4pi_ny
#define PSD_N mccPSD_Sphere_4pi_PSD_N
#define PSD_p mccPSD_Sphere_4pi_PSD_p
#define PSD_p2 mccPSD_Sphere_4pi_PSD_p2
#define filename mccPSD_Sphere_4pi_filename
#define radius mccPSD_Sphere_4pi_radius
#define restore_neutron mccPSD_Sphere_4pi_restore_neutron
#define nowritefile mccPSD_Sphere_4pi_nowritefile
#line 61 "/usr/share/mcstas/2.6rc1/monitors/PSD_monitor_4PI.comp"
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
#line 18887 "Samples_Incoherent.c"
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

  /* Initializations for component Coord3. */
  SIG_MESSAGE("Coord3 (Init)");

  /* Initializations for component Sphere_4pi. */
  SIG_MESSAGE("Sphere_4pi (Init)");
#define mccompcurname  Sphere_4pi
#define mccompcurtype  PSD_monitor_4PI
#define mccompcurindex 18
#define nx mccSphere_4pi_nx
#define ny mccSphere_4pi_ny
#define PSD_N mccSphere_4pi_PSD_N
#define PSD_p mccSphere_4pi_PSD_p
#define PSD_p2 mccSphere_4pi_PSD_p2
#define filename mccSphere_4pi_filename
#define radius mccSphere_4pi_radius
#define restore_neutron mccSphere_4pi_restore_neutron
#define nowritefile mccSphere_4pi_nowritefile
#line 61 "/usr/share/mcstas/2.6rc1/monitors/PSD_monitor_4PI.comp"
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
#line 18930 "Samples_Incoherent.c"
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
  /* TRACE Component Source [1] */
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
#define mccompcurname  Source
#define mccompcurtype  ISIS_moderator
#define mccompcurindex 1
#define Face mccSource_Face
#define p_in mccSource_p_in
#define Tnpts mccSource_Tnpts
#define scaleSize mccSource_scaleSize
#define angleArea mccSource_angleArea
#define Nsim mccSource_Nsim
#define Ncount mccSource_Ncount
#define TS mccSource_TS
#define rtE0 mccSource_rtE0
#define rtE1 mccSource_rtE1
#define rtmodX mccSource_rtmodX
#define rtmodY mccSource_rtmodY
#define TargetStation mccSource_TargetStation
#define CurrentWeight mccSource_CurrentWeight
{   /* Declarations of Source=ISIS_moderator() SETTING parameters. */
MCNUM Emin = mccSource_Emin;
MCNUM Emax = mccSource_Emax;
MCNUM dist = mccSource_dist;
MCNUM focus_xw = mccSource_focus_xw;
MCNUM focus_yh = mccSource_focus_yh;
MCNUM xwidth = mccSource_xwidth;
MCNUM yheight = mccSource_yheight;
MCNUM CAngle = mccSource_CAngle;
MCNUM SAC = mccSource_SAC;
MCNUM Lmin = mccSource_Lmin;
MCNUM Lmax = mccSource_Lmax;
int target_index = mccSource_target_index;
MCNUM verbose = mccSource_verbose;
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
#line 19129 "Samples_Incoherent.c"
}   /* End of Source=ISIS_moderator() SETTING parameter declarations. */
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
  mcabsorbCompSource:
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

  /* TRACE Component Coord [2] */
  mccoordschange(mcposrCoord, mcrotrCoord,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component Coord (without coords transformations) */
  mcJumpTrace_Coord:
  SIG_MESSAGE("Coord (Trace)");
  mcDEBUG_COMP("Coord")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompCoord
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
#define mccompcurname  Coord
#define mccompcurtype  Progress_bar
#define mccompcurindex 2
#define IntermediateCnts mccCoord_IntermediateCnts
#define StartTime mccCoord_StartTime
#define EndTime mccCoord_EndTime
#define CurrentTime mccCoord_CurrentTime
{   /* Declarations of Coord=Progress_bar() SETTING parameters. */
char* profile = mccCoord_profile;
MCNUM percent = mccCoord_percent;
MCNUM flag_save = mccCoord_flag_save;
MCNUM minutes = mccCoord_minutes;
#line 70 "/usr/share/mcstas/2.6rc1/misc/Progress_bar.comp"
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
#line 19301 "Samples_Incoherent.c"
}   /* End of Coord=Progress_bar() SETTING parameter declarations. */
#undef CurrentTime
#undef EndTime
#undef StartTime
#undef IntermediateCnts
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompCoord:
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

  /* TRACE Component Vsample [3] */
  mccoordschange(mcposrVsample, mcrotrVsample,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component Vsample (without coords transformations) */
  mcJumpTrace_Vsample:
  SIG_MESSAGE("Vsample (Trace)");
  mcDEBUG_COMP("Vsample")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompVsample
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
#define mccompcurname  Vsample
#define mccompcurtype  V_sample
#define mccompcurindex 3
#define VarsV mccVsample_VarsV
{   /* Declarations of Vsample=V_sample() SETTING parameters. */
MCNUM radius = mccVsample_radius;
MCNUM thickness = mccVsample_thickness;
MCNUM zdepth = mccVsample_zdepth;
MCNUM Vc = mccVsample_Vc;
MCNUM sigma_abs = mccVsample_sigma_abs;
MCNUM sigma_inc = mccVsample_sigma_inc;
MCNUM radius_i = mccVsample_radius_i;
MCNUM radius_o = mccVsample_radius_o;
MCNUM h = mccVsample_h;
MCNUM focus_r = mccVsample_focus_r;
MCNUM pack = mccVsample_pack;
MCNUM frac = mccVsample_frac;
MCNUM f_QE = mccVsample_f_QE;
MCNUM gamma = mccVsample_gamma;
MCNUM target_x = mccVsample_target_x;
MCNUM target_y = mccVsample_target_y;
MCNUM target_z = mccVsample_target_z;
MCNUM focus_xw = mccVsample_focus_xw;
MCNUM focus_yh = mccVsample_focus_yh;
MCNUM focus_aw = mccVsample_focus_aw;
MCNUM focus_ah = mccVsample_focus_ah;
MCNUM xwidth = mccVsample_xwidth;
MCNUM yheight = mccVsample_yheight;
MCNUM zthick = mccVsample_zthick;
MCNUM rad_sphere = mccVsample_rad_sphere;
MCNUM sig_a = mccVsample_sig_a;
MCNUM sig_i = mccVsample_sig_i;
MCNUM V0 = mccVsample_V0;
int target_index = mccVsample_target_index;
MCNUM multiples = mccVsample_multiples;
/* 'Vsample=V_sample()' component instance has conditional execution */
if (( mcipSAMPLE == 1 ))

#line 180 "/usr/share/mcstas/2.6rc1/obsolete/V_sample.comp"
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
#line 19583 "Samples_Incoherent.c"
}   /* End of Vsample=V_sample() SETTING parameter declarations. */
#undef VarsV
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompVsample:
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

  /* TRACE Component PowN [4] */
  mccoordschange(mcposrPowN, mcrotrPowN,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component PowN (without coords transformations) */
  mcJumpTrace_PowN:
  SIG_MESSAGE("PowN (Trace)");
  mcDEBUG_COMP("PowN")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompPowN
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
#define mccompcurname  PowN
#define mccompcurtype  PowderN
#define mccompcurindex 4
#define format mccPowN_format
#define line_info mccPowN_line_info
#define columns mccPowN_columns
#define offdata mccPowN_offdata
{   /* Declarations of PowN=PowderN() SETTING parameters. */
char* reflections = mccPowN_reflections;
char* geometry = mccPowN_geometry;
MCNUM radius = mccPowN_radius;
MCNUM yheight = mccPowN_yheight;
MCNUM xwidth = mccPowN_xwidth;
MCNUM zdepth = mccPowN_zdepth;
MCNUM thickness = mccPowN_thickness;
MCNUM pack = mccPowN_pack;
MCNUM Vc = mccPowN_Vc;
MCNUM sigma_abs = mccPowN_sigma_abs;
MCNUM sigma_inc = mccPowN_sigma_inc;
MCNUM delta_d_d = mccPowN_delta_d_d;
MCNUM p_inc = mccPowN_p_inc;
MCNUM p_transmit = mccPowN_p_transmit;
MCNUM DW = mccPowN_DW;
MCNUM nb_atoms = mccPowN_nb_atoms;
MCNUM d_omega = mccPowN_d_omega;
MCNUM d_phi = mccPowN_d_phi;
MCNUM tth_sign = mccPowN_tth_sign;
MCNUM p_interact = mccPowN_p_interact;
MCNUM concentric = mccPowN_concentric;
MCNUM density = mccPowN_density;
MCNUM weight = mccPowN_weight;
MCNUM barns = mccPowN_barns;
MCNUM Strain = mccPowN_Strain;
MCNUM focus_flip = mccPowN_focus_flip;
int target_index = mccPowN_target_index;
/* 'PowN=PowderN()' component instance has conditional execution */
if (( mcipSAMPLE == 2 ))

#line 773 "/usr/share/mcstas/2.6rc1/samples/PowderN.comp"
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
#line 20034 "Samples_Incoherent.c"
}   /* End of PowN=PowderN() SETTING parameter declarations. */
#undef offdata
#undef columns
#undef line_info
#undef format
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompPowN:
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

  /* TRACE Component SX [5] */
  mccoordschange(mcposrSX, mcrotrSX,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component SX (without coords transformations) */
  mcJumpTrace_SX:
  SIG_MESSAGE("SX (Trace)");
  mcDEBUG_COMP("SX")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompSX
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
#define mccompcurname  SX
#define mccompcurtype  Single_crystal
#define mccompcurindex 5
#define mosaic_AB mccSX_mosaic_AB
#define hkl_info mccSX_hkl_info
#define offdata mccSX_offdata
{   /* Declarations of SX=Single_crystal() SETTING parameters. */
char* reflections = mccSX_reflections;
char* geometry = mccSX_geometry;
MCNUM xwidth = mccSX_xwidth;
MCNUM yheight = mccSX_yheight;
MCNUM zdepth = mccSX_zdepth;
MCNUM radius = mccSX_radius;
MCNUM delta_d_d = mccSX_delta_d_d;
MCNUM mosaic = mccSX_mosaic;
MCNUM mosaic_a = mccSX_mosaic_a;
MCNUM mosaic_b = mccSX_mosaic_b;
MCNUM mosaic_c = mccSX_mosaic_c;
MCNUM recip_cell = mccSX_recip_cell;
MCNUM barns = mccSX_barns;
MCNUM ax = mccSX_ax;
MCNUM ay = mccSX_ay;
MCNUM az = mccSX_az;
MCNUM bx = mccSX_bx;
MCNUM by = mccSX_by;
MCNUM bz = mccSX_bz;
MCNUM cx = mccSX_cx;
MCNUM cy = mccSX_cy;
MCNUM cz = mccSX_cz;
MCNUM p_transmit = mccSX_p_transmit;
MCNUM sigma_abs = mccSX_sigma_abs;
MCNUM sigma_inc = mccSX_sigma_inc;
MCNUM aa = mccSX_aa;
MCNUM bb = mccSX_bb;
MCNUM cc = mccSX_cc;
MCNUM order = mccSX_order;
MCNUM RX = mccSX_RX;
MCNUM RY = mccSX_RY;
MCNUM powder = mccSX_powder;
MCNUM PG = mccSX_PG;
MCNUM deltak = mccSX_deltak;
/* 'SX=Single_crystal()' component instance has conditional execution */
if (( mcipSAMPLE == 3 ))

#line 1066 "/usr/share/mcstas/2.6rc1/samples/Single_crystal.comp"
{
  double t1, t2=0;                /* Entry and exit times in sample */
  struct hkl_data *L;           /* Structure factor list */
  int i;                        /* Index into structure factor list */
  struct tau_data *T;           /* List of reflections close to Ewald sphere */
  int j;                        /* Index into reflection list */
  int event_counter;            /* scattering event counter */
  double kix, kiy, kiz, ki;     /* Initial wave vector [1/AA] */
  double kfx, kfy, kfz;         /* Final wave vector */
  double v;                     /* Neutron velocity */
  double rho_x, rho_y, rho_z;   /* the vector ki - tau */
  double rho;
  double diff;                  /* Deviation from Bragg condition */
  double ox, oy, oz;            /* Origin of Ewald sphere tangent plane */
  double b1x, b1y, b1z;         /* First vector spanning tangent plane */
  double b2x, b2y, b2z;         /* Second vector spanning tangent plane */
  double n11, n12, n22;         /* 2D Gauss description matrix N */
  double det_N;                 /* Determinant of N */
  double inv_n11, inv_n12, inv_n22; /* Inverse of N */
  double l11, l12, l22;         /* Cholesky decomposition L of 1/2*inv(N) */
  double det_L;                 /* Determinant of L */
  double Bt_D_O_x, Bt_D_O_y;    /* Temporaries */
  double y0x, y0y;              /* Center of 2D Gauss in plane coordinates */
  double alpha;                 /* Offset of 2D Gauss center from 3D center */
  double V0;                    /* Volume of unit cell */
  double l_full;                /* Neutron path length for transmission */
  double l;                     /* Path length to scattering event */
  double abs_xsect, abs_xlen;   /* Absorption cross section and length */
  double inc_xsect, inc_xlen;   /* Incoherent scattering cross section and length */
  double coh_xlen;              /* Coherent cross section and length */
  double tot_xsect, tot_xlen;   /* Total cross section and length */
  double z1, z2, y1, y2;        /* Temporaries to choose kf from 2D Gauss */
  double adjust, r, sum;        /* Temporaries */

  double p_trans;               /* Transmission probability */
  double mc_trans, mc_interact; /* Transmission, interaction MC choices */
  int    intersect=0;
  double theta, phi;            /* rotation angles for curved lattice option */

  double curv_xangle;
  double curv_yangle;

  double _vx;
  double _vy;
  double _vz;


  /* Intersection neutron trajectory / sample (sample surface) */
  if (hkl_info.shape == 0)
    intersect = cylinder_intersect(&t1, &t2, x, y, z, vx, vy, vz, radius, yheight);
  else if (hkl_info.shape == 1)
    intersect = box_intersect(&t1, &t2, x, y, z, vx, vy, vz, xwidth, yheight, zdepth);
  else if (hkl_info.shape == 2)
    intersect = sphere_intersect(&t1, &t2, x, y, z, vx, vy, vz, radius);
  else if (hkl_info.shape == 3)
    intersect = off_intersect(&t1, &t2, NULL, NULL, x, y, z, vx, vy, vz, offdata );

  if (t2 < 0) intersect=0;  /* we passed sample volume already */

  if(intersect)
  {                         /* Neutron intersects crystal */
    if(t1 > 0)
      PROP_DT(t1);          /* Move to crystal surface if not inside */
    v  = sqrt(vx*vx + vy*vy + vz*vz);
    ki = V2K*v;
    event_counter = 0;
    abs_xsect = hkl_info.sigma_a*2200/v;
    inc_xsect = hkl_info.sigma_i;
    V0= hkl_info.V0;
    abs_xlen  = abs_xsect/V0;
    inc_xlen  = inc_xsect/V0;

    /* Scalar cross sections for inc/abs are given in barns, so we need a scaling factor of 100
       to get scattering lengths in m, since V0 is assumed to be in AA*/
    abs_xlen *= 100; inc_xlen *= 100;

    L = hkl_info.list;
    T = hkl_info.tau_list;
    hkl_info.type = '\0';

    do {  /* Loop over multiple scattering events */
      /* Angles for powder randomization */
      double Alpha, Beta, Gamma;

      if (hkl_info.shape == 0)
        intersect = cylinder_intersect(&t1, &t2, x, y, z, vx, vy, vz, radius, yheight);
      else if (hkl_info.shape == 1)
        intersect = box_intersect(&t1, &t2, x, y, z, vx, vy, vz, xwidth, yheight, zdepth);
      else if (hkl_info.shape == 2)
        intersect = sphere_intersect(&t1, &t2, x, y, z, vx, vy, vz, radius);
      else if (hkl_info.shape == 3)
        intersect = off_intersect(&t1, &t2, NULL, NULL, x, y, z, vx, vy, vz, offdata );
      if(!intersect || t2*v < -1e-9 || t1*v > 1e-9)
      {
        /* neutron is leaving the sample */
        if (hkl_info.flag_warning < 100)
          fprintf(stderr,
                "Single_crystal: %s: Warning: neutron has unexpectedly left the crystal!\n"
                "                t1=%g t2=%g x=%g y=%g z=%g vx=%g vy=%g vz=%g\n",
                NAME_CURRENT_COMP, t1, t2, x, y, z, vx, vy, vz);
        hkl_info.flag_warning++;
        break;
      }

      l_full = t2*v;

      /* (1). Compute incoming wave vector ki */
      if (powder) { /* orientation of crystallite is random */
        Alpha = randpm1()*PI*powder;
        Beta  = randpm1()*PI/2;
        Gamma = randpm1()*PI;
        randrotate(&vx, &vy, &vz, Alpha, Beta, Gamma);
      }
      if (PG) { /* orientation of crystallite is random along <c> axis */
        Alpha = rand01()*2*PI*PG;
        PGrotate(&vx, &vy, &vz, Alpha, hkl_info.csx, hkl_info.csy, hkl_info.csz);
      }



      /* ------------------------------------------------------------------------- */
      /* lattice curvature option: rotate neutron velocity */
      /* WARNING: cannot be used together with the PG c-rotation! */
      curv_xangle = 0;
      curv_yangle = 0;

      _vx = vx;
      _vy = vy;
      _vz = vz;

      if(RY) { /* rotate v around x axis based on y pos, for vertical focus */
          curv_yangle = atan2(y, RY);
          vec_rotate_2d(&vy,&vz, curv_yangle);
          vec_rotate_2d(&sy,&sz, curv_yangle);

          /*changing y,z actually curves the crystal, not only the planes*/
          /*comment out if only curvature of the lattice planes is needed*/
          vec_rotate_2d(&y,&z, curv_yangle);
      }
      if(RX) { /* rotate v around y axis based on x pos, for horizontal focus */
          curv_xangle = atan2(x, RX);
          vec_rotate_2d(&vx,&vz, curv_xangle);
          vec_rotate_2d(&sx,&sz, curv_xangle);

          /*changing x,z actually curves the crystal, not only the planes*/
          /*comment out if only curvature of the lattice planes is needed*/
          vec_rotate_2d(&x,&z, curv_xangle);
      }

      kix = V2K*vx;
      kiy = V2K*vy;
      kiz = V2K*vz;
      vx = _vx;
      vy = _vy;
      vz = _vz;
      /* ------------------------------------------------------------------------- */



      /* (2). Intersection of Ewald sphere with reciprocal lattice points */

      /* in case we use 'SPLIT' then consecutive neutrons can be identical when entering here
         and we may skip the hkl_search call */
      if ( fabs(kix - hkl_info.kix) < deltak
        && fabs(kiy - hkl_info.kiy) < deltak
        && fabs(kiz - hkl_info.kiz) < deltak) {
        hkl_info.nb_reuses++;
      } else {
        /* Max possible tau for this ki with 5*sigma delta-d/d cutoff. */
        double tau_max   = 2*ki/(1 - 5*hkl_info.m_delta_d_d);
        double coh_xsect = 0, coh_refl = 0;

        /* call hkl_search */
        hkl_info.tau_count = hkl_search(L, T, hkl_info.count, hkl_info.V0,
              kix, kiy, kiz, tau_max,
              &coh_refl, &coh_xsect); /* CPU consuming */

        /* store ki so that we can check for further SPLIT iterations */
        if (event_counter == 0 ) { /* only for incoming neutron */
          hkl_info.kix = kix;
          hkl_info.kiy = kiy;
          hkl_info.kiz = kiz;
        }

        hkl_info.coh_refl  = coh_refl;
        hkl_info.coh_xsect = coh_xsect;
        hkl_info.nb_refl += hkl_info.tau_count;
        hkl_info.nb_refl_count++;
      }

      /* (3). Probabilities of the different possible interactions. */
      /* Cross-sections are in barns = 10**-28 m**2, and unit cell volumes are
         in AA**3 = 10**-30 m**2. Hence a factor of 100 is used to convert
         scattering lengths to m**-1 */
      coh_xlen = hkl_info.coh_xsect/V0;
      if (barns) {
        coh_xlen *= 100;
      } /* else assume fm^2 */
      tot_xlen = abs_xlen + inc_xlen + coh_xlen;

      /* (5). Transmission */
      p_trans = exp(-tot_xlen*l_full);
      if(!event_counter && p_transmit >= 0 && p_transmit <= 1) {
        mc_trans = p_transmit; /* first event */
      } else {
        mc_trans = p_trans;
      }
      hkl_info.type = 't';
      mc_interact = 1 - mc_trans;
      if(mc_trans > 0 && (mc_trans >= 1 || rand01() < mc_trans))  /* Transmit */
      {
        p *= p_trans/mc_trans;
        intersect=0;
        if (powder) { /* orientation of crystallite is longer random */
          randderotate(&vx, &vy, &vz, Alpha, Beta, Gamma);
        }
        if (PG) { /* orientation of crystallite is longer random */
          PGderotate(&vx, &vy, &vz, Alpha, hkl_info.m_cx, hkl_info.m_cz, hkl_info.m_cz);
        }
        break;
      }
      if(tot_xlen <= 0){
        ABSORB;
      }
      if(mc_interact <= 0)        /* Protect against rounding errors */
        { intersect=0;
          if (powder) { /* orientation of crystallite is no longer random */
            randderotate(&vx, &vy, &vz, Alpha, Beta, Gamma);
          }
          if (PG) { /* orientation of crystallite is no longer random, rotation around <c> */
            PGderotate(&vx, &vy, &vz, Alpha, hkl_info.csx, hkl_info.csy, hkl_info.csz);
          }
          break;
        }
      if (!event_counter) p *= fabs(1 - p_trans)/mc_interact;
      /* Select a point at which to scatter the neutron, taking
         secondary extinction into account. */
      /* dP(l) = exp(-tot_xlen*l)dl
         P(l<l_0) = [-1/tot_xlen*exp(-tot_xlen*l)]_0^l_0
                  = (1 - exp(-tot_xlen*l0))/tot_xlen
         l = -log(1 - tot_xlen*rand0max(P(l<l_full)))/tot_xlen
       */
      if(tot_xlen*l_full < 1e-6)
        /* For very weak scattering, use simple uniform sampling of scattering
           point to avoid rounding errors. */
        l = rand0max(l_full);
      else
        l = -log(1 - rand0max((1 - exp(-tot_xlen*l_full))))/tot_xlen;
      PROP_DT(l/v);
      event_counter++;

      /* (4). Account for the probability of sigma_abs */
      p *= (coh_xlen + inc_xlen)/tot_xlen;
      /* Choose between coherent and incoherent scattering */
      if(coh_xlen == 0 || rand0max(coh_xlen + inc_xlen) <= inc_xlen)
      {
        /* (6). Incoherent scattering */
        randvec_target_circle(&kix, &kiy, &kiz, NULL, vx, vy, vz, 0);
        vx = kix; /* ki vector is used as tmp var with norm v */
        vy = kiy;
        vz = kiz; /* Go for next scattering event */
        hkl_info.type = 'i';
      } else {
        /* 7. Coherent scattering. Select reciprocal lattice point. */
        if(hkl_info.coh_refl <= 0){
          ABSORB;
        }
        sum = 0;
        j = hkl_select(T, hkl_info.tau_count, hkl_info.coh_refl, &sum);

        if(j >= hkl_info.tau_count)
        {
          if (hkl_info.flag_warning < 100)
            fprintf(stderr, "Single_crystal: Error: Illegal tau search "
              "(r=%g, sum=%g, j=%i, tau_count=%i).\n", r, sum, j , hkl_info.tau_count);
          hkl_info.flag_warning++;
          j = hkl_info.tau_count - 1;
        }
        i = T[j].index;
        /* (8). Pick scattered wavevector kf from 2D Gauss distribution. */
        z1 = randnorm();
        z2 = randnorm();
        y1 = T[j].l11*z1 + T[j].y0x;
        y2 = T[j].l12*z1 + T[j].l22*z2 + T[j].y0y;
        kfx = T[j].rho_x + T[j].ox + T[j].b1x*y1 + T[j].b2x*y2;
        kfy = T[j].rho_y + T[j].oy + T[j].b1y*y1 + T[j].b2y*y2;
        kfz = T[j].rho_z + T[j].oz + T[j].b1z*y1 + T[j].b2z*y2;
        /* Normalize kf to length of ki, to account for planer
          approximation of the Ewald sphere. */
        adjust = ki/sqrt(kfx*kfx + kfy*kfy + kfz*kfz);
        kfx *= adjust;
        kfy *= adjust;
        kfz *= adjust;
        /* Adjust neutron weight (see manual for explanation). */
        p *= T[j].xsect*hkl_info.coh_refl/(hkl_info.coh_xsect*T[j].refl);

        vx = K2V*(L[i].u1x*kfx + L[i].u2x*kfy + L[i].u3x*kfz);
        vy = K2V*(L[i].u1y*kfx + L[i].u2y*kfy + L[i].u3y*kfz);
        vz = K2V*(L[i].u1z*kfx + L[i].u2z*kfy + L[i].u3z*kfz);
        hkl_info.type = 'c';
        hkl_info.h    = L[i].h;
        hkl_info.k    = L[i].k;
        hkl_info.l    = L[i].l;
      }



      /* ------------------------------------------------------------------------- */
      /* lattice curvature option: rotate back neutron velocity */
      if(RX) {
          vec_rotate_2d(&vx,&vz, -curv_xangle);
          vec_rotate_2d(&sx,&sz, -curv_xangle);

          /*changing x,z actually curves the crystal, not only the planes*/
          /*comment out if only curvature of the lattice planes is needed*/
          vec_rotate_2d(&x,&z, -curv_xangle);
      }
      if(RY) {
          vec_rotate_2d(&vy,&vz, -curv_yangle);
          vec_rotate_2d(&sy,&sz, -curv_yangle);

          /*changing y,z actually curves the crystal, not only the planes*/
          /*comment out if only curvature of the lattice planes is needed*/
          vec_rotate_2d(&y,&z, -curv_yangle);
      }
      /* ------------------------------------------------------------------------- */



      SCATTER;
      if (powder) { /* orientation of crystallite is no longer random */
        randderotate(&vx, &vy, &vz, Alpha, Beta, Gamma);
      }
      if (PG) { /* orientation of crystallite is longer random */
        PGderotate(&vx, &vy, &vz, Alpha, hkl_info.csx, hkl_info.csy, hkl_info.csz);
      }
      /* exit if multiple scattering order has been reached */
      if (order && event_counter >= order) { intersect=0; break; }
      /* Repeat loop for next scattering event. */
    } while (intersect); /* end do (intersect) (multiple scattering loop) */
  } /* if intersect */
}
#line 20526 "Samples_Incoherent.c"
}   /* End of SX=Single_crystal() SETTING parameter declarations. */
#undef offdata
#undef hkl_info
#undef mosaic_AB
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompSX:
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

  /* TRACE Component Sqw [6] */
  mccoordschange(mcposrSqw, mcrotrSqw,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component Sqw (without coords transformations) */
  mcJumpTrace_Sqw:
  SIG_MESSAGE("Sqw (Trace)");
  mcDEBUG_COMP("Sqw")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompSqw
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
#define mccompcurname  Sqw
#define mccompcurtype  Isotropic_Sqw
#define mccompcurindex 6
#define powder_format mccSqw_powder_format
#define VarSqw mccSqw_VarSqw
#define columns mccSqw_columns
#define offdata mccSqw_offdata
{   /* Declarations of Sqw=Isotropic_Sqw() SETTING parameters. */
char* Sqw_coh = mccSqw_Sqw_coh;
char* Sqw_inc = mccSqw_Sqw_inc;
char* geometry = mccSqw_geometry;
MCNUM radius = mccSqw_radius;
MCNUM thickness = mccSqw_thickness;
MCNUM xwidth = mccSqw_xwidth;
MCNUM yheight = mccSqw_yheight;
MCNUM zdepth = mccSqw_zdepth;
MCNUM threshold = mccSqw_threshold;
int order = mccSqw_order;
MCNUM T = mccSqw_T;
MCNUM verbose = mccSqw_verbose;
MCNUM d_phi = mccSqw_d_phi;
int concentric = mccSqw_concentric;
MCNUM rho = mccSqw_rho;
MCNUM sigma_abs = mccSqw_sigma_abs;
MCNUM sigma_coh = mccSqw_sigma_coh;
MCNUM sigma_inc = mccSqw_sigma_inc;
MCNUM classical = mccSqw_classical;
MCNUM powder_Dd = mccSqw_powder_Dd;
MCNUM powder_DW = mccSqw_powder_DW;
MCNUM powder_Vc = mccSqw_powder_Vc;
MCNUM density = mccSqw_density;
MCNUM weight = mccSqw_weight;
MCNUM p_interact = mccSqw_p_interact;
MCNUM norm = mccSqw_norm;
MCNUM powder_barns = mccSqw_powder_barns;
char* quantum_correction = mccSqw_quantum_correction;
/* 'Sqw=Isotropic_Sqw()' component instance has conditional execution */
if (( mcipSAMPLE == 4 ))

#line 2176 "/usr/share/mcstas/2.6rc1/samples/Isotropic_Sqw.comp"
{

int    intersect=0;     /* flag to continue/stop */
double t0,  t1,  t2,  t3; /* times for intersections */
double dt0, dt1, dt2, dt; /* time intervals */
double k=0, Ei=0;
double v=0, vf=0;
double d_path;        /* total path length for straight trajectory */
double my_a;          /* absorption cross-section scaled to velocity (2200) */
double ws, p_scatt;   /* probability for scattering/absorption and for */
                      /* interaction along d_path */
double tmp_rand;      /* temporary var */
double ratio_w=0, ratio_q=0; /* variables for bilinear interpolation */
double q11, q21, q22, q12;
double omega=0;       /* energy transfer */
double q=0;           /* wavevector transfer */
long   index_w;       /* energy index for table look-up SW */
long   index_q;       /* Q index for table look-up P(Q|w) */
double theta=0, costheta=0; /* for the choice of kf direction */
double u1x,u1y,u1z;
double u2x,u2y,u2z;
double u0x,u0y,u0z;
int    index_counter;
int    flag=0;
int    flag_concentric=0;
int    flag_ishollow=0;
double solid_angle=0;
double my_t=0;
double p_mult=1;
double mc_trans, p_trans, mc_scatt;
double coh=0, inc=0;
struct Sqw_Data_struct Data_sqw;


/* Store Initial neutron state */

VarSqw.ki_x = V2K*vx;
VarSqw.ki_y = V2K*vy;
VarSqw.ki_z = V2K*vz;
VarSqw.ti   = t;
VarSqw.vi   = 0;
VarSqw.ki   = 0;
VarSqw.type = '\0';

do { /* Main interaction loop. Ends with intersect=0 */

  /* Intersection neutron trajectory / sample (sample surface) */
  if (VarSqw.s_coh > 0 || VarSqw.s_inc > 0) {
    if (thickness >= 0) {
      if (VarSqw.shape==0)
        intersect=cylinder_intersect(&t0,&t3, x,y,z,vx,vy,vz, radius,yheight);
      else if (VarSqw.shape==1)
        intersect=box_intersect     (&t0,&t3, x,y,z,vx,vy,vz, xwidth,yheight,zdepth);
      else if (VarSqw.shape==2)
        intersect=sphere_intersect  (&t0,&t3, x,y,z,vx,vy,vz, radius);
      else if (VarSqw.shape == 3)
        intersect=off_intersect(&t0, &t3, NULL, NULL, x, y, z, vx, vy, vz, offdata );
    } else {
      if (VarSqw.shape==0)
        intersect=cylinder_intersect(&t0,&t3, x,y,z,vx,vy,vz, radius-thickness,
          yheight-2*thickness > 0 ? yheight-2*thickness : yheight);
      else if (VarSqw.shape==1)
        intersect=box_intersect     (&t0,&t3, x,y,z,vx,vy,vz,
          xwidth-2*thickness > 0 ?  xwidth-2*thickness : xwidth,
          yheight-2*thickness > 0 ? yheight-2*thickness : yheight,
          zdepth-2*thickness > 0 ?  zdepth-2*thickness : zdepth);
      else if (VarSqw.shape==2)
        intersect=sphere_intersect  (&t0,&t3, x,y,z,vx,vy,vz, radius-thickness);
      else if (VarSqw.shape == 3)
        intersect=off_intersect(&t0, &t3, NULL, NULL, x, y, z, vx, vy, vz, offdata );
    }
  } else intersect=0;

  /* Computing the intermediate times */
  if (intersect) {
    flag_ishollow = 0;
    if (thickness > 0) {
      if (VarSqw.shape==0 && cylinder_intersect(&t1,&t2, x,y,z,vx,vy,vz, radius-thickness,
        yheight-2*thickness > 0 ? yheight-2*thickness : yheight))
        flag_ishollow=1;
      else if (VarSqw.shape==2 && sphere_intersect   (&t1,&t2, x,y,z,vx,vy,vz, radius-thickness))
        flag_ishollow=1;
      else if (VarSqw.shape==1 && box_intersect(&t1,&t2, x,y,z,vx,vy,vz,
        xwidth-2*thickness > 0 ? xwidth-2*thickness : xwidth,
        yheight-2*thickness > 0 ? yheight-2*thickness : yheight,
        zdepth-2*thickness > 0 ? zdepth-2*thickness : zdepth))
        flag_ishollow=1;
    } else if (thickness<0) {
      if (VarSqw.shape==0 && cylinder_intersect(&t1,&t2, x,y,z,vx,vy,vz, radius,yheight))
        flag_ishollow=1;
      else if (VarSqw.shape==2 && sphere_intersect   (&t1,&t2, x,y,z,vx,vy,vz, radius))
        flag_ishollow=1;
      else if (VarSqw.shape==1 && box_intersect(&t1,&t2, x,y,z,vx,vy,vz, xwidth, yheight, zdepth))
        flag_ishollow=1;
    }
    if (!flag_ishollow) t1 = t2 = t3; /* no empty space inside */
  } else break; /* neutron does not hit sample: transmitted  */

  if (intersect) { /* the neutron hits the sample */

    if (t0 > 0) {  /* we are before the sample */
      PROP_DT(t0); /* propagates neutron to the entry of the sample */
    } else if (t1 > 0 && t1 > t0) { /* we are inside first part of the sample */
      /* no propagation, stay inside */
    } else if (t2 > 0 && t2 > t1) { /* we are in the hole */
      PROP_DT(t2); /* propagate to inner surface of 2nd part of sample */
    } else if (t3 > 0 && t3 > t2) { /* we are in the 2nd part of sample */
      /* no propagation, stay inside */
    }

    dt0=t1-(t0 > 0 ? t0 : 0); /* Time in first part of hollow/cylinder/box */
    dt1=t2-(t1 > 0 ? t1 : 0); /* Time in hole */
    dt2=t3-(t2 > 0 ? t2 : 0); /* Time in 2nd part of hollow cylinder */

    if (dt0 < 0) dt0 = 0;
    if (dt1 < 0) dt1 = 0;
    if (dt2 < 0) dt2 = 0;

    /* initialize concentric mode */
    if (concentric && !flag_concentric && t0 >= 0
     && VarSqw.shape==0 && thickness) {
      flag_concentric=1;
    }

    if (flag_concentric == 1) {
      dt1=dt2=0; /* force exit when reaching hole/2nd part */
    }

    if (!dt0 && !dt2) {
      intersect = 0; /* the sample was passed entirely */
      break;
    }

    VarSqw.neutron_enter++;
    p_mult = 1;
    if (!v) {
      v  = vx*vx+vy*vy+vz*vz;
      v = sqrt(v);
    }
    k  = V2K*v;
    Ei = VS2E*v*v;

    if (!VarSqw.vi) VarSqw.vi = v;
    if (!VarSqw.ki) VarSqw.ki = k;

    if (v <= 0) {
      printf("Isotropic_Sqw: %s: ERROR: Null velocity !\n",NAME_CURRENT_COMP);
      VarSqw.neutron_removed++;
      ABSORB; /* should never occur */
    }

    /* check for scattering event */
    my_a   = VarSqw.my_a_v / v; /* absorption 'mu' */
    /* compute total scattering X section */
    /* \int q S(q) dq /2 /ki^2 sigma  OR  bare Xsection*/
    /* contains the 4*PI*kf/ki factor */
    coh = VarSqw.s_coh;
    inc = VarSqw.s_inc;
    if (k && VarSqw.s_coh>0 && VarSqw.Data_coh.intensity) {
      double Ei       = VS2E*v*v;
      double index_Ei = Ei / (VarSqw.Data_coh.Ei_max/VarSqw.Data_coh.iqSq_length);
      coh = Table_Value2d(VarSqw.Data_coh.iqSq, index_Ei, 0);
    }
    if (k && VarSqw.s_inc>0 && VarSqw.Data_inc.intensity) {
      double Ei       = VS2E*v*v;
      double index_Ei = Ei / (VarSqw.Data_inc.Ei_max/VarSqw.Data_inc.iqSq_length);
      inc = Table_Value2d(VarSqw.Data_inc.iqSq, index_Ei, 0);
    }
    if (coh<0) coh=0;
    if (inc<0) inc=0;
    VarSqw.my_s    =(VarSqw.mat_rho*100*(coh + inc));

    my_t = my_a + VarSqw.my_s;  /* total scattering Xsect */
    if (my_t <= 0) {
      if (VarSqw.neutron_removed<VarSqw.maxloop) printf("Isotropic_Sqw: %s: ERROR: Null total cross section %g. Removing event.\n",
        NAME_CURRENT_COMP, my_t);
      VarSqw.neutron_removed++;
      ABSORB; /* should never occur */
    } else if (VarSqw.my_s <= 0) {
      if (VarSqw.verbose_output > 1 && VarSqw.neutron_removed<VarSqw.maxloop)
        printf("Isotropic_Sqw: %s: Warning: Null scattering cross section %g. Ignoring.\n",
          NAME_CURRENT_COMP, VarSqw.my_s);
      VarSqw.my_s = 0;
    }

    /* Proba of scattering vs absorption (integrating along the whole trajectory) */
    ws = VarSqw.my_s/my_t;  /* (inc+coh)/(inc+coh+abs) */
    d_path = v*( dt0 +dt2 );    /* total path lenght in sample */
    /* Proba of transmission/interaction along length d_path */
    p_trans = exp(-my_t*d_path);
    p_scatt = 1 - p_trans; /* portion of beam which scatters */

    flag = 0; /* flag used for propagation to exit point before ending */

    /* are we next to the exit ? probably no scattering (avoid rounding errors) */
    if (VarSqw.my_s*d_path <= 4e-7) {
      flag = 1;           /* No interaction before the exit */
    }
    /* force a given fraction of the beam to scatter */
    if (p_interact>0 && p_interact<=1) {
      /* we force a portion of the beam to interact */
      /* This is used to improve statistics on single scattering (and multiple) */
      if (!SCATTERED) mc_trans = 1-p_interact;
      else            mc_trans = 1-p_interact/(4*SCATTERED+1); /* reduce effect on multi scatt */
    } else {
      mc_trans = p_trans; /* 1 - p_scatt */
    }
    mc_scatt = 1 - mc_trans; /* portion of beam to scatter (or force to) */
    if (mc_scatt <= 0 || mc_scatt>1) flag=1;
    /* MC choice: Interaction or transmission ? */
    if (!flag && mc_scatt > 0 && (mc_scatt >= 1 || rand01() < mc_scatt)) { /* Interaction neutron/sample */
      p_mult *= ws; /* Update weight ; account for absorption and retain scattered fraction */
      /* we have chosen portion mc_scatt of beam instead of p_scatt, so we compensate */
      if (!mc_scatt) ABSORB;
      p_mult *= fabs(p_scatt/mc_scatt); /* lower than 1 */
    } else {
      flag = 1; /* Transmission : no interaction neutron/sample */
      if (!VarSqw.type) VarSqw.type = 't';
      if (!mc_trans) ABSORB;
      p_mult *= fabs(p_trans/mc_trans);  /* attenuate beam by portion which is scattered (and left along) */
    }

    if (flag) { /* propagate to exit of sample and finish */
      intersect = 0;
      p *= p_mult; /* apply absorption correction */
      PROP_DT(dt0+dt2);
      break; /* exit main multi scatt while loop */
    }
  } /* end if intersect the neutron hits the sample */
  else break;

  if (intersect) { /* scattering event */
    double kf=0, kf1, kf2;
    /* mean scattering probability and absorption fraction */
    VarSqw.mean_scatt += (1-exp(-VarSqw.my_s*d_path))*p;
    VarSqw.mean_abs   += (1-ws)*p;
    VarSqw.psum_scatt += p;

    /* Decaying exponential distribution of the path length before scattering */
    /* Select a point at which to scatter the neutron, taking
         secondary extinction into account. */
    if (my_t*d_path < 1e-6)
    /* For very weak scattering, use simple uniform sampling of scattering
       point to avoid rounding errors. */
      dt = rand0max(d_path); /* length */
    else
      dt = -log(1 - rand0max((1 - exp(-my_t*d_path)))) / my_t; /* length */
    dt /= v; /* Time from present position to scattering point */

    /* If t0 is in hole, propagate to next part of the hollow cylinder */
    if (dt1 > 0 && dt0 > 0 && dt > dt0) dt += dt1;

    /* Neutron propagation to the scattering point */
    PROP_DT(dt);

    /* choice between coherent/incoherent scattering */
    tmp_rand = rand01();
    /* local description at the scattering point (scat probability for atom) */
    tmp_rand *= (coh+inc);

    flag=0;
    if (VarSqw.s_inc>0 && tmp_rand < inc) {
      /* CASE 1: incoherent case */
      if (!VarSqw.Data_inc.intensity) {
        /* CASE 1a: no incoherent Sqw from file, use isotropic V-like */
        if (d_phi && order == 1) {
          randvec_target_rect_angular(&u1x, &u1y, &u1z, &solid_angle,
              vx, vy, vz, 2*PI, d_phi, ROT_A_CURRENT_COMP);
          p_mult *= solid_angle/4/PI; /* weighted by focused range to total range */
        } else
          randvec_target_circle(&u1x, &u1y, &u1z, NULL, vx, vy, vz, 0);

        vx = u1x; vy = u1y; vz = u1z;
        vf = v; kf = k;
        if (!VarSqw.type) VarSqw.type = 'v';
        SCATTER;
      } else {
        /* CASE 1b: incoherent Sqw from file */
        if (VarSqw.Data_inc.intensity) {
          Data_sqw = VarSqw.Data_inc;
          if (!VarSqw.type) VarSqw.type = 'i';
          flag = 1;
        }
      }
    } else if (VarSqw.s_coh>0 && tmp_rand > VarSqw.s_inc) {
      if (VarSqw.Data_coh.intensity) {
        /* CASE2: coherent case */
        Data_sqw = VarSqw.Data_coh;
        if (!VarSqw.type) VarSqw.type = 'c';
        flag = 1;
      }
    }

    if (flag) { /* true when S(q,w) table exists (Data_sqw) */

      double alpha=0, alpha0;
      /* give us a limited number of tries for scattering: choose W then Q */
      for (index_counter=VarSqw.maxloop; index_counter > 0 ; index_counter--) {

        /* MC choice: energy transfer w=Ei-Ef in the S(w) = SW */
        omega = 0;
        tmp_rand = rand01();
        /* energy index for rand > cumul SW */
        index_w  = Sqw_search_SW(Data_sqw, tmp_rand);
        VarSqw.rw = (double)index_w;
        if (index_w >= 0 && &(Data_sqw.SW[index_w]) != NULL) {
          if (Data_sqw.w_bins > 1) {
            double w1, w2;
            if (index_w > 0) { /* interpolate linearly energy */
              ratio_w = (tmp_rand                         - Data_sqw.SW[index_w-1].cumul_proba)
                       /(Data_sqw.SW[index_w].cumul_proba - Data_sqw.SW[index_w-1].cumul_proba);
              /* ratio_w=0 omega[index_w-1], ratio=1 omega[index] */
              w1 = Data_sqw.SW[index_w-1].omega; w2 = Data_sqw.SW[index_w].omega;
            } else { /* index_w = 0 interpolate to 0 energy */
              /* ratio_w=0 omega=0, ratio=1 omega[index] */
              w1 = Data_sqw.SW[index_w].omega; w2= Data_sqw.SW[index_w+1].omega;
              if (!w2 && index_w+1 < Data_sqw.w_bins)
                w2= Data_sqw.SW[index_w+1].omega;
              if (Data_sqw.w_bins && Data_sqw.SW[index_w].cumul_proba) {
                ratio_w = tmp_rand/Data_sqw.SW[index_w].cumul_proba;
              } else ratio_w=0;
            }
            if (ratio_w<0) ratio_w=0; else if (ratio_w>1) ratio_w=1;
            omega = (1-ratio_w)*w1 + ratio_w*w2;
          } else {
            ratio_w = 0;
            omega = Data_sqw.SW[index_w].omega;
          }
        } else {
          if (VarSqw.verbose_output >= 3 && VarSqw.neutron_removed<VarSqw.maxloop)
            printf("Isotropic_Sqw: %s: Warning: No suitable w transfer for index_w=%li.\n",
              NAME_CURRENT_COMP, index_w);
          continue; /* no W value: try again with an other energy transfer */
        }

        /* MC choice: momentum transfer Q in P(Q|w) */
        tmp_rand = rand01();

        /* momentum index for rand > cumul SQ|W */
        index_q  = Sqw_search_Q_proba_per_w(Data_sqw, tmp_rand, index_w);
        VarSqw.rq = (double)index_q;

        if (index_q >= 0 && &(Data_sqw.SQW[index_w]) != NULL) {
          if (Data_sqw.q_bins > 1 && index_q > 0) {
            if (index_w > 0 && Data_sqw.w_bins > 1) {
              /* bilinear interpolation on - side: index_w > 0, index_q > 0 */
              ratio_q = (tmp_rand - Data_sqw.SQW[index_w][index_q-1].cumul_proba)
                       /(Data_sqw.SQW[index_w][index_q].cumul_proba
                       - Data_sqw.SQW[index_w][index_q-1].cumul_proba);
              q22 = Data_sqw.SQW[index_w]  [index_q].Q;
              q11 = Data_sqw.SQW[index_w-1][index_q-1].Q;
              q21 = Data_sqw.SQW[index_w]  [index_q-1].Q;
              q12 = Data_sqw.SQW[index_w-1][index_q].Q;
              if (ratio_q<0) ratio_q=0; else if (ratio_q>1) ratio_q=1;
              q = (1-ratio_w)*(1-ratio_q)*q11+ratio_w*(1-ratio_q)*q21
                + ratio_w*ratio_q*q22        +(1-ratio_w)*ratio_q*q12;
            } else { /* bilinear interpolation on + side: index_w=0, index_q > 0 */
              ratio_q = (tmp_rand - Data_sqw.SQW[index_w][index_q-1].cumul_proba)
                       /(Data_sqw.SQW[index_w][index_q].cumul_proba
                       - Data_sqw.SQW[index_w][index_q-1].cumul_proba);
              q11 = Data_sqw.SQW[index_w]  [index_q-1].Q;
              q12 = Data_sqw.SQW[index_w]  [index_q].Q;
              if (ratio_q<0) ratio_q=0; else if (ratio_q>1) ratio_q=1;
              if (index_w < Data_sqw.w_bins-1 && Data_sqw.w_bins > 1) {
                q22 = Data_sqw.SQW[index_w+1][index_q].Q;
                q21 = Data_sqw.SQW[index_w+1][index_q-1].Q;
                q = (1-ratio_w)*(1-ratio_q)*q11+ratio_w*(1-ratio_q)*q21
                  + ratio_w*ratio_q*q22        +(1-ratio_w)*ratio_q*q12;
              } else {
                q    = (1-ratio_q)*q11  + ratio_q*q12;
              }
            }
          } else {
            q    = Data_sqw.SQW[index_w][index_q].Q;
          }
        } else {
          if (VarSqw.verbose_output >= 3 && VarSqw.neutron_removed<VarSqw.maxloop)
            printf("Isotropic_Sqw: %s: Warning: No suitable q transfer for w=%g.\n",
              NAME_CURRENT_COMP, omega);
          VarSqw.neutron_removed++;
          continue; /* no Q value for this w choice */
        }

        /* Search for length of final wave vector kf */
        /* kf is such that : hbar*w = hbar*hbar/2/m*(k*k - kf*kf) */
        /* acceptable values for kf are kf1 and kf2 */
        if (!solve_2nd_order(&kf1, &kf2, 1, 0, -k*k + VarSqw.sqSE2K*omega)) {
          if (VarSqw.verbose_output >= 3 && VarSqw.neutron_removed<VarSqw.maxloop)
            printf("Isotropic_Sqw: %s: Warning: imaginary root for w=%g q=%g Ei=%g (triangle can not close)\n",
            NAME_CURRENT_COMP, omega, q, Ei);
          VarSqw.neutron_removed++;
          continue; /* all roots are imaginary */
        }

        /* kf1 and kf2 are opposite */
        kf = fabs(kf1);
        vf = K2V*kf;

        /* Search of the direction of kf such that : q = ki - kf */
        /* cos theta = (ki2+kf2-q2)/(2ki kf) */

        costheta= (k*k+kf*kf-q*q)/(2*kf*k); /* this is cos(theta) */

        if (-1 < costheta && costheta < 1) {
          break; /* satisfies q momentum conservation */
        }
/*      else continue; */

        /* exit for loop on success */
      } /* end for index_counter */

      if (!index_counter) { /* for loop ended: failure for scattering */
        intersect=0; /* Could not scatter: finish multiple scattering loop */
        if (VarSqw.verbose_output >= 2 && VarSqw.neutron_removed<VarSqw.maxloop)
          printf("Isotropic_Sqw: %s: Warning: No scattering [q,w] conditions\n"
               "               last try (%i): type=%c w=%g q=%g cos(theta)=%g k=%g\n",
          NAME_CURRENT_COMP, VarSqw.maxloop, (VarSqw.type ? VarSqw.type : '-'), omega, q, costheta, k);
        VarSqw.neutron_removed++;
        if (order && SCATTERED != order) ABSORB;
        break;       /* finish multiple scattering loop */
      }

      /* scattering angle from ki to DS cone */
      theta = acos(costheta);

      /* Choose point on Debye-Scherrer cone */
      if (order == 1 && d_phi)
      { /* relate height of detector to the height on DS cone */
        double cone_focus;
        cone_focus = sin(d_phi/2)/sin(theta);
        /* If full Debye-Scherrer cone is within d_phi, don't focus */
        if (cone_focus < -1 || cone_focus > 1) d_phi = 0;
        /* Otherwise, determine alpha to rotate from scattering plane
            into d_phi focusing area*/
        else alpha = 2*asin(cone_focus);
        if (d_phi) p_mult *= alpha/PI;
      }
      if (d_phi) {
        /* Focusing */
        alpha = fabs(alpha);
        /* Trick to get scattering for pos/neg theta's */
        alpha0= 2*rand01()*alpha;
        if (alpha0 > alpha) {
          alpha0=PI+(alpha0-1.5*alpha);
        } else {
          alpha0=alpha0-0.5*alpha;
        }
      }
      else
        alpha0 = PI*randpm1();

      /* now find a nearly vertical rotation axis (u1) :
	       * Either
	       *  (v along Z) x (X axis) -> nearly Y axis
	       * Or
	       *  (v along X) x (Z axis) -> nearly Y axis
	       */
	    if (fabs(scalar_prod(1,0,0,vx/v,vy/v,vz/v)) < fabs(scalar_prod(0,0,1,vx/v,vy/v,vz/v))) {
        u1x = 1; u1y = u1z = 0;
	    } else {
        u1x = u1y = 0; u1z = 1;
	    }
	    vec_prod(u2x,u2y,u2z, vx,vy,vz, u1x,u1y,u1z);

      /* handle case where v and aim are parallel */
      if (!u2x && !u2y && !u2z) { u2x=u2z=0; u2y=1; }

      /* u1 = rotate 'v' by theta around u2: DS scattering angle, nearly in horz plane */
      rotate(u1x,u1y,u1z, vx,vy,vz, theta, u2x,u2y,u2z);

      /* u0 = rotate u1 by alpha0 around 'v' (Debye-Scherrer cone) */
      rotate(u0x,u0y,u0z, u1x,u1y,u1z, alpha0, vx, vy, vz);
      NORM(u0x,u0y,u0z);
      vx = u0x*vf;
      vy = u0y*vf;
      vz = u0z*vf;

      SCATTER;

      v = vf; k = kf; /* for next iteration */

    } /* end if (flag) */

    VarSqw.neutron_exit++;
    p *= p_mult;
    if (p_mult > 1) VarSqw.neutron_pmult++;

    /* test for a given multiple order */
    if (order && SCATTERED >= order) {
      intersect=0; /* reached required number of SCATTERing */
      break;       /* finish multiple scattering loop */
    }

  } /* end if (intersect) scattering event  */

} while (intersect); /* end do (intersect) (multiple scattering loop) */

/* Store Final neutron state */
VarSqw.kf_x = V2K*vx;
VarSqw.kf_y = V2K*vy;
VarSqw.kf_z = V2K*vz;
VarSqw.tf   = t;
VarSqw.vf   = v;
VarSqw.kf   = k;
VarSqw.theta= theta;

if (SCATTERED) {



  if (SCATTERED == 1) {
    if (VarSqw.type == 'c') VarSqw.single_coh += p;
    else                    VarSqw.single_inc += p;
    VarSqw.dq = sqrt((VarSqw.kf_x-VarSqw.ki_x)*(VarSqw.kf_x-VarSqw.ki_x)
                  +(VarSqw.kf_y-VarSqw.ki_y)*(VarSqw.kf_y-VarSqw.ki_y)
                  +(VarSqw.kf_z-VarSqw.ki_z)*(VarSqw.kf_z-VarSqw.ki_z));
    VarSqw.dw = VS2E*(VarSqw.vf*VarSqw.vf - VarSqw.vi*VarSqw.vi);
  } else VarSqw.multi += p;

} else VarSqw.dq=VarSqw.dw=0;

/* end TRACE */
}
#line 21193 "Samples_Incoherent.c"
}   /* End of Sqw=Isotropic_Sqw() SETTING parameter declarations. */
#undef offdata
#undef columns
#undef VarSqw
#undef powder_format
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompSqw:
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

  /* TRACE Component IncS [7] */
  mccoordschange(mcposrIncS, mcrotrIncS,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component IncS (without coords transformations) */
  mcJumpTrace_IncS:
  SIG_MESSAGE("IncS (Trace)");
  mcDEBUG_COMP("IncS")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompIncS
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
#define mccompcurname  IncS
#define mccompcurtype  Incoherent
#define mccompcurindex 7
#define VarsInc mccIncS_VarsInc
#define offdata mccIncS_offdata
{   /* Declarations of IncS=Incoherent() SETTING parameters. */
char* geometry = mccIncS_geometry;
MCNUM radius = mccIncS_radius;
MCNUM xwidth = mccIncS_xwidth;
MCNUM yheight = mccIncS_yheight;
MCNUM zdepth = mccIncS_zdepth;
MCNUM thickness = mccIncS_thickness;
MCNUM target_x = mccIncS_target_x;
MCNUM target_y = mccIncS_target_y;
MCNUM target_z = mccIncS_target_z;
MCNUM focus_r = mccIncS_focus_r;
MCNUM focus_xw = mccIncS_focus_xw;
MCNUM focus_yh = mccIncS_focus_yh;
MCNUM focus_aw = mccIncS_focus_aw;
MCNUM focus_ah = mccIncS_focus_ah;
int target_index = mccIncS_target_index;
MCNUM pack = mccIncS_pack;
MCNUM p_interact = mccIncS_p_interact;
MCNUM f_QE = mccIncS_f_QE;
MCNUM gamma = mccIncS_gamma;
MCNUM Etrans = mccIncS_Etrans;
MCNUM deltaE = mccIncS_deltaE;
MCNUM sigma_abs = mccIncS_sigma_abs;
MCNUM sigma_inc = mccIncS_sigma_inc;
MCNUM Vc = mccIncS_Vc;
MCNUM concentric = mccIncS_concentric;
MCNUM order = mccIncS_order;
/* 'IncS=Incoherent()' component instance has conditional execution */
if (( mcipSAMPLE == 5 ))

#line 232 "/usr/share/mcstas/2.6rc1/samples/Incoherent.comp"
{
  double t0, t3;                /* Entry/exit time for outer surface */
  double t1, t2;                /* Entry/exit time for inner surface */
  double dt0, dt1, dt2, dt;     /* Flight times through sample */
  double v=0;                   /* Neutron velocity */
  double d_path;                /* Flight path length for non-scattered neutron */
  double l_i, l_o=0;            /* Flight path lenght in/out for scattered neutron */
  double my_a=0,my_t=0;         /* Velocity-dependent attenuation factor and total Xsec */
  double solid_angle=0;         /* Solid angle of target as seen from scattering point */
  double aim_x=0, aim_y=0, aim_z=1;   /* Position of target relative to scattering point */
  double v_i, v_f, E_i, E_f; /* initial and final energies and velocities */
  double dE;                 /* Energy transfer */
  int    intersect=0;
  int    flag_concentric=0;
  int    flag=0;
  double mc_trans, p_trans, mc_scatt, p_scatt, ws;
  double p_mult=1;

  do { /* Main interaction loop. Ends with intersect=0 */

    /* Intersection neutron trajectory / sample (sample surface) */
    if (VarsInc.shape == 0)
      intersect = cylinder_intersect(&t0, &t3, x, y, z, vx, vy, vz, radius, yheight);
    else if (VarsInc.shape == 1)
      intersect = box_intersect(&t0, &t3, x, y, z, vx, vy, vz, xwidth, yheight, zdepth);
    else if (VarsInc.shape == 2)
      intersect = sphere_intersect(&t0, &t3, x, y, z, vx, vy, vz, radius);
    else if (VarsInc.shape == 3)
      intersect = off_intersect(&t0, &t3, NULL, NULL, x, y, z, vx, vy, vz, offdata );
    if (intersect) {
      int flag_ishollow = 0;
      if (thickness>0) {
        if (VarsInc.shape==0 && cylinder_intersect(&t1,&t2, x,y,z,vx,vy,vz, radius-thickness,yheight-2*thickness))
          flag_ishollow=1;
        else if (VarsInc.shape==2 && sphere_intersect   (&t1,&t2, x,y,z,vx,vy,vz, radius-thickness))
          flag_ishollow=1;
        else if (VarsInc.shape==1 && box_intersect(&t1,&t2, x,y,z,vx,vy,vz, xwidth-2*thickness, yheight-2*thickness, zdepth-2*thickness))
          flag_ishollow = 1;
      }
      if (!flag_ishollow) t1 = t2 = t3; /* no empty space inside */

      dt0 = t1-t0;                /* Time in sample, ingoing */
      dt1 = t2-t1;                /* Time in hole */
      dt2 = t3-t2;                /* Time in sample, outgoing */

      if (t0 > 0) {  /* we are before the sample */
        PROP_DT(t0); /* propagates neutron to the entry of the sample */
      } else if (t1 > 0 && t1 > t0) { /* we are inside first part of the sample */
        /* no propagation, stay inside */
      } else if (t2 > 0 && t2 > t1) { /* we are in the hole */
        PROP_DT(t2); /* propagate to inner surface of 2nd part of sample */
      } else if (t3 > 0 && t3 > t2) { /* we are in the 2nd part of sample */
        /* no propagation, stay inside */
      }

      dt0=t1-(t0 > 0 ? t0 : 0); /* Time in first part of hollow/cylinder/box */
      dt1=t2-(t1 > 0 ? t1 : 0); /* Time in hole */
      dt2=t3-(t2 > 0 ? t2 : 0); /* Time in 2nd part of hollow cylinder */

      if (dt0 < 0) dt0 = 0;
      if (dt1 < 0) dt1 = 0;
      if (dt2 < 0) dt2 = 0;

      /* initialize concentric mode */
      if (concentric && !flag_concentric && t0 >= 0
       && VarsInc.shape==0 && thickness>0) {
        flag_concentric=1;
      }

      if (flag_concentric == 1) {
        dt1=dt2=0; /* force exit when reaching hole/2nd part */
      }

      if (!dt0 && !dt2) {
        intersect = 0; /* the sample was passed entirely */
        break;
      }

      p_mult = 1;
      if (!v) v = sqrt(vx*vx + vy*vy + vz*vz);
      if (v) my_a = VarsInc.my_a_v*(2200/v);
      else {
        printf("Incoherent: %s: ERROR: Null velocity\n",NAME_CURRENT_COMP);
        ABSORB; /* should never occur */
      }

      my_t = my_a + VarsInc.my_s;  /* total scattering Xsect (tmp var) */
      if (my_t <= 0) {
        printf("Incoherent: %s: ERROR: Null total cross section %g. Removing event.\n",
          NAME_CURRENT_COMP, my_t);
        ABSORB; /* should never occur */
      }
      d_path = v * (dt0 + dt2);   /* Length of full path through sample */
      /* Proba of scattering vs absorption (integrating along the whole trajectory) */
      ws = VarsInc.my_s/my_t;  /* (inc+coh)/(inc+coh+abs) */
      /* Proba of transmission along length d_path */
      p_trans = exp(-my_t*d_path);
      p_scatt = 1 - p_trans; /* portion of beam which scatters */
      flag = 0; /* flag used for propagation to exit point before ending */
      /* are we next to the exit ? probably no scattering (avoid rounding errors) */
      if (VarsInc.my_s*d_path <= 4e-7) {
        flag = 1;           /* No interaction before the exit */
      }
      /* force a given fraction of the beam to scatter */
      if (p_interact>0 && p_interact<=1) {
        /* we force a portion of the beam to interact */
        /* This is used to improve statistics on single scattering (and multiple) */
        if (!SCATTERED) mc_trans = 1-p_interact;
        else            mc_trans = 1-p_interact/(4*SCATTERED+1); /* reduce effect on multi scatt */
      } else {
        mc_trans = p_trans; /* 1 - p_scatt */
      }
      mc_scatt = 1 - mc_trans; /* portion of beam to scatter (or force to) */
      if (mc_scatt <= 0 || mc_scatt>1) flag=1;
      /* MC choice: Interaction or transmission ? */
      if (!flag && mc_scatt > 0 && (mc_scatt >= 1 || (rand01()) < mc_scatt)) { /* Interaction neutron/sample */
        p_mult *= ws; /* Update weight ; account for absorption and retain scattered fraction */
        if (!mc_scatt) ABSORB;
        /* we have chosen portion mc_scatt of beam instead of p_scatt, so we compensate */
        p_mult *= fabs(p_scatt/mc_scatt); /* lower than 1 */
      } else {
        flag = 1; /* Transmission : no interaction neutron/sample */
        if (!mc_trans) ABSORB;
        p_mult *= fabs(p_trans/mc_trans);  /* attenuate beam by portion which is scattered (and left along) */
      }

      if (flag) { /* propagate to exit of sample and finish */
        intersect = 0;
        p *= p_mult; /* apply absorption correction */
        PROP_DT(dt0+dt2);
        break; /* exit main multi scatt while loop */
      }
      if (my_t*d_path < 1e-6)
      /* For very weak scattering, use simple uniform sampling of scattering
         point to avoid rounding errors. */
        dt = rand0max(d_path); /* length */
      else
        dt = -log(1 - rand0max((1 - exp(-my_t*d_path)))) / my_t; /* length */
      l_i = dt;/* Penetration in sample: scattering+abs */
      dt /= v; /* Time from present position to scattering point */

      /* If t0 is in hole, propagate to next part of the hollow cylinder */
      if (dt1 > 0 && dt0 > 0 && dt > dt0) dt += dt1;

      PROP_DT(dt);                /* Point of scattering */

      if ((VarsInc.tx || VarsInc.ty || VarsInc.tz)) {
        aim_x = VarsInc.tx-x;       /* Vector pointing at target (anal./det.) */
        aim_y = VarsInc.ty-y;
        aim_z = VarsInc.tz-z;
      }
      if(VarsInc.aw && VarsInc.ah) {
        randvec_target_rect_angular(&vx, &vy, &vz, &solid_angle,
          aim_x, aim_y, aim_z, VarsInc.aw, VarsInc.ah, ROT_A_CURRENT_COMP);
      } else if(VarsInc.xw && VarsInc.yh) {
        randvec_target_rect(&vx, &vy, &vz, &solid_angle,
          aim_x, aim_y, aim_z, VarsInc.xw, VarsInc.yh, ROT_A_CURRENT_COMP);
      } else {
        randvec_target_circle(&vx, &vy, &vz, &solid_angle, aim_x, aim_y, aim_z, focus_r);
      }
      NORM(vx, vy, vz);

      v_i = v;          /* Store initial velocity in case of quasielastic */
      E_i = VS2E*v_i*v_i;
      if (deltaE==0) { 
	if (rand01()<f_QE)	/* Quasielastic contribution */
	  {
	    dE = gamma*tan(PI/2*randpm1());
	    E_f = E_i + dE;
	    if (E_f <= 0)
	      ABSORB;
	    v_f = SE2V*sqrt(E_f);
	    v = v_f;
	    /*          printf("vi: %g Ei: %g dE: %g Ef %g vf: %g v: %g \n",
			v_i,E_i,dE,E_f,v_f,v); */
	  }
      } else {
	E_f = E_i - Etrans + deltaE*randpm1(); // E_f=E0;
	v_f = SE2V*sqrt(E_f);
	v = v_f;
      }

      vx *= v;
      vy *= v;
      vz *= v;

      /* We do not consider scattering from 2nd part (outgoing) */
      p_mult *= solid_angle/4/PI;
      p *= p_mult;

      /* Polarisation part (1/3 NSF, 2/3 SF) */
      sx *= -1.0/3.0;
      sy *= -1.0/3.0;
      sz *= -1.0/3.0;

      SCATTER;

      /* test for a given multiple order */
      if (order && SCATTERED >= order) {
        intersect=0; /* reached required number of SCATTERing */
        break;       /* finish multiple scattering loop */
      }
    } /* end if intersect */
  } while (intersect); /* end do (intersect) (multiple scattering loop) */
}
#line 21539 "Samples_Incoherent.c"
}   /* End of IncS=Incoherent() SETTING parameter declarations. */
#undef offdata
#undef VarsInc
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompIncS:
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

  /* TRACE Component PowNg [8] */
  mccoordschange(mcposrPowNg, mcrotrPowNg,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component PowNg (without coords transformations) */
  mcJumpTrace_PowNg:
  SIG_MESSAGE("PowNg (Trace)");
  mcDEBUG_COMP("PowNg")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompPowNg
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
#define mccompcurname  PowNg
#define mccompcurtype  PowderN
#define mccompcurindex 8
#define format mccPowNg_format
#define line_info mccPowNg_line_info
#define columns mccPowNg_columns
#define offdata mccPowNg_offdata
{   /* Declarations of PowNg=PowderN() SETTING parameters. */
char* reflections = mccPowNg_reflections;
char* geometry = mccPowNg_geometry;
MCNUM radius = mccPowNg_radius;
MCNUM yheight = mccPowNg_yheight;
MCNUM xwidth = mccPowNg_xwidth;
MCNUM zdepth = mccPowNg_zdepth;
MCNUM thickness = mccPowNg_thickness;
MCNUM pack = mccPowNg_pack;
MCNUM Vc = mccPowNg_Vc;
MCNUM sigma_abs = mccPowNg_sigma_abs;
MCNUM sigma_inc = mccPowNg_sigma_inc;
MCNUM delta_d_d = mccPowNg_delta_d_d;
MCNUM p_inc = mccPowNg_p_inc;
MCNUM p_transmit = mccPowNg_p_transmit;
MCNUM DW = mccPowNg_DW;
MCNUM nb_atoms = mccPowNg_nb_atoms;
MCNUM d_omega = mccPowNg_d_omega;
MCNUM d_phi = mccPowNg_d_phi;
MCNUM tth_sign = mccPowNg_tth_sign;
MCNUM p_interact = mccPowNg_p_interact;
MCNUM concentric = mccPowNg_concentric;
MCNUM density = mccPowNg_density;
MCNUM weight = mccPowNg_weight;
MCNUM barns = mccPowNg_barns;
MCNUM Strain = mccPowNg_Strain;
MCNUM focus_flip = mccPowNg_focus_flip;
int target_index = mccPowNg_target_index;
/* 'PowNg=PowderN()' component instance has conditional execution */
if (( mcipSAMPLE == 6 ))

#line 773 "/usr/share/mcstas/2.6rc1/samples/PowderN.comp"
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
#line 21991 "Samples_Incoherent.c"
}   /* End of PowNg=PowderN() SETTING parameter declarations. */
#undef offdata
#undef columns
#undef line_info
#undef format
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompPowNg:
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

  /* TRACE Component SXg [9] */
  mccoordschange(mcposrSXg, mcrotrSXg,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component SXg (without coords transformations) */
  mcJumpTrace_SXg:
  SIG_MESSAGE("SXg (Trace)");
  mcDEBUG_COMP("SXg")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompSXg
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
#define mccompcurname  SXg
#define mccompcurtype  Single_crystal
#define mccompcurindex 9
#define mosaic_AB mccSXg_mosaic_AB
#define hkl_info mccSXg_hkl_info
#define offdata mccSXg_offdata
{   /* Declarations of SXg=Single_crystal() SETTING parameters. */
char* reflections = mccSXg_reflections;
char* geometry = mccSXg_geometry;
MCNUM xwidth = mccSXg_xwidth;
MCNUM yheight = mccSXg_yheight;
MCNUM zdepth = mccSXg_zdepth;
MCNUM radius = mccSXg_radius;
MCNUM delta_d_d = mccSXg_delta_d_d;
MCNUM mosaic = mccSXg_mosaic;
MCNUM mosaic_a = mccSXg_mosaic_a;
MCNUM mosaic_b = mccSXg_mosaic_b;
MCNUM mosaic_c = mccSXg_mosaic_c;
MCNUM recip_cell = mccSXg_recip_cell;
MCNUM barns = mccSXg_barns;
MCNUM ax = mccSXg_ax;
MCNUM ay = mccSXg_ay;
MCNUM az = mccSXg_az;
MCNUM bx = mccSXg_bx;
MCNUM by = mccSXg_by;
MCNUM bz = mccSXg_bz;
MCNUM cx = mccSXg_cx;
MCNUM cy = mccSXg_cy;
MCNUM cz = mccSXg_cz;
MCNUM p_transmit = mccSXg_p_transmit;
MCNUM sigma_abs = mccSXg_sigma_abs;
MCNUM sigma_inc = mccSXg_sigma_inc;
MCNUM aa = mccSXg_aa;
MCNUM bb = mccSXg_bb;
MCNUM cc = mccSXg_cc;
MCNUM order = mccSXg_order;
MCNUM RX = mccSXg_RX;
MCNUM RY = mccSXg_RY;
MCNUM powder = mccSXg_powder;
MCNUM PG = mccSXg_PG;
MCNUM deltak = mccSXg_deltak;
/* 'SXg=Single_crystal()' component instance has conditional execution */
if (( mcipSAMPLE == 7 ))

#line 1066 "/usr/share/mcstas/2.6rc1/samples/Single_crystal.comp"
{
  double t1, t2=0;                /* Entry and exit times in sample */
  struct hkl_data *L;           /* Structure factor list */
  int i;                        /* Index into structure factor list */
  struct tau_data *T;           /* List of reflections close to Ewald sphere */
  int j;                        /* Index into reflection list */
  int event_counter;            /* scattering event counter */
  double kix, kiy, kiz, ki;     /* Initial wave vector [1/AA] */
  double kfx, kfy, kfz;         /* Final wave vector */
  double v;                     /* Neutron velocity */
  double rho_x, rho_y, rho_z;   /* the vector ki - tau */
  double rho;
  double diff;                  /* Deviation from Bragg condition */
  double ox, oy, oz;            /* Origin of Ewald sphere tangent plane */
  double b1x, b1y, b1z;         /* First vector spanning tangent plane */
  double b2x, b2y, b2z;         /* Second vector spanning tangent plane */
  double n11, n12, n22;         /* 2D Gauss description matrix N */
  double det_N;                 /* Determinant of N */
  double inv_n11, inv_n12, inv_n22; /* Inverse of N */
  double l11, l12, l22;         /* Cholesky decomposition L of 1/2*inv(N) */
  double det_L;                 /* Determinant of L */
  double Bt_D_O_x, Bt_D_O_y;    /* Temporaries */
  double y0x, y0y;              /* Center of 2D Gauss in plane coordinates */
  double alpha;                 /* Offset of 2D Gauss center from 3D center */
  double V0;                    /* Volume of unit cell */
  double l_full;                /* Neutron path length for transmission */
  double l;                     /* Path length to scattering event */
  double abs_xsect, abs_xlen;   /* Absorption cross section and length */
  double inc_xsect, inc_xlen;   /* Incoherent scattering cross section and length */
  double coh_xlen;              /* Coherent cross section and length */
  double tot_xsect, tot_xlen;   /* Total cross section and length */
  double z1, z2, y1, y2;        /* Temporaries to choose kf from 2D Gauss */
  double adjust, r, sum;        /* Temporaries */

  double p_trans;               /* Transmission probability */
  double mc_trans, mc_interact; /* Transmission, interaction MC choices */
  int    intersect=0;
  double theta, phi;            /* rotation angles for curved lattice option */

  double curv_xangle;
  double curv_yangle;

  double _vx;
  double _vy;
  double _vz;


  /* Intersection neutron trajectory / sample (sample surface) */
  if (hkl_info.shape == 0)
    intersect = cylinder_intersect(&t1, &t2, x, y, z, vx, vy, vz, radius, yheight);
  else if (hkl_info.shape == 1)
    intersect = box_intersect(&t1, &t2, x, y, z, vx, vy, vz, xwidth, yheight, zdepth);
  else if (hkl_info.shape == 2)
    intersect = sphere_intersect(&t1, &t2, x, y, z, vx, vy, vz, radius);
  else if (hkl_info.shape == 3)
    intersect = off_intersect(&t1, &t2, NULL, NULL, x, y, z, vx, vy, vz, offdata );

  if (t2 < 0) intersect=0;  /* we passed sample volume already */

  if(intersect)
  {                         /* Neutron intersects crystal */
    if(t1 > 0)
      PROP_DT(t1);          /* Move to crystal surface if not inside */
    v  = sqrt(vx*vx + vy*vy + vz*vz);
    ki = V2K*v;
    event_counter = 0;
    abs_xsect = hkl_info.sigma_a*2200/v;
    inc_xsect = hkl_info.sigma_i;
    V0= hkl_info.V0;
    abs_xlen  = abs_xsect/V0;
    inc_xlen  = inc_xsect/V0;

    /* Scalar cross sections for inc/abs are given in barns, so we need a scaling factor of 100
       to get scattering lengths in m, since V0 is assumed to be in AA*/
    abs_xlen *= 100; inc_xlen *= 100;

    L = hkl_info.list;
    T = hkl_info.tau_list;
    hkl_info.type = '\0';

    do {  /* Loop over multiple scattering events */
      /* Angles for powder randomization */
      double Alpha, Beta, Gamma;

      if (hkl_info.shape == 0)
        intersect = cylinder_intersect(&t1, &t2, x, y, z, vx, vy, vz, radius, yheight);
      else if (hkl_info.shape == 1)
        intersect = box_intersect(&t1, &t2, x, y, z, vx, vy, vz, xwidth, yheight, zdepth);
      else if (hkl_info.shape == 2)
        intersect = sphere_intersect(&t1, &t2, x, y, z, vx, vy, vz, radius);
      else if (hkl_info.shape == 3)
        intersect = off_intersect(&t1, &t2, NULL, NULL, x, y, z, vx, vy, vz, offdata );
      if(!intersect || t2*v < -1e-9 || t1*v > 1e-9)
      {
        /* neutron is leaving the sample */
        if (hkl_info.flag_warning < 100)
          fprintf(stderr,
                "Single_crystal: %s: Warning: neutron has unexpectedly left the crystal!\n"
                "                t1=%g t2=%g x=%g y=%g z=%g vx=%g vy=%g vz=%g\n",
                NAME_CURRENT_COMP, t1, t2, x, y, z, vx, vy, vz);
        hkl_info.flag_warning++;
        break;
      }

      l_full = t2*v;

      /* (1). Compute incoming wave vector ki */
      if (powder) { /* orientation of crystallite is random */
        Alpha = randpm1()*PI*powder;
        Beta  = randpm1()*PI/2;
        Gamma = randpm1()*PI;
        randrotate(&vx, &vy, &vz, Alpha, Beta, Gamma);
      }
      if (PG) { /* orientation of crystallite is random along <c> axis */
        Alpha = rand01()*2*PI*PG;
        PGrotate(&vx, &vy, &vz, Alpha, hkl_info.csx, hkl_info.csy, hkl_info.csz);
      }



      /* ------------------------------------------------------------------------- */
      /* lattice curvature option: rotate neutron velocity */
      /* WARNING: cannot be used together with the PG c-rotation! */
      curv_xangle = 0;
      curv_yangle = 0;

      _vx = vx;
      _vy = vy;
      _vz = vz;

      if(RY) { /* rotate v around x axis based on y pos, for vertical focus */
          curv_yangle = atan2(y, RY);
          vec_rotate_2d(&vy,&vz, curv_yangle);
          vec_rotate_2d(&sy,&sz, curv_yangle);

          /*changing y,z actually curves the crystal, not only the planes*/
          /*comment out if only curvature of the lattice planes is needed*/
          vec_rotate_2d(&y,&z, curv_yangle);
      }
      if(RX) { /* rotate v around y axis based on x pos, for horizontal focus */
          curv_xangle = atan2(x, RX);
          vec_rotate_2d(&vx,&vz, curv_xangle);
          vec_rotate_2d(&sx,&sz, curv_xangle);

          /*changing x,z actually curves the crystal, not only the planes*/
          /*comment out if only curvature of the lattice planes is needed*/
          vec_rotate_2d(&x,&z, curv_xangle);
      }

      kix = V2K*vx;
      kiy = V2K*vy;
      kiz = V2K*vz;
      vx = _vx;
      vy = _vy;
      vz = _vz;
      /* ------------------------------------------------------------------------- */



      /* (2). Intersection of Ewald sphere with reciprocal lattice points */

      /* in case we use 'SPLIT' then consecutive neutrons can be identical when entering here
         and we may skip the hkl_search call */
      if ( fabs(kix - hkl_info.kix) < deltak
        && fabs(kiy - hkl_info.kiy) < deltak
        && fabs(kiz - hkl_info.kiz) < deltak) {
        hkl_info.nb_reuses++;
      } else {
        /* Max possible tau for this ki with 5*sigma delta-d/d cutoff. */
        double tau_max   = 2*ki/(1 - 5*hkl_info.m_delta_d_d);
        double coh_xsect = 0, coh_refl = 0;

        /* call hkl_search */
        hkl_info.tau_count = hkl_search(L, T, hkl_info.count, hkl_info.V0,
              kix, kiy, kiz, tau_max,
              &coh_refl, &coh_xsect); /* CPU consuming */

        /* store ki so that we can check for further SPLIT iterations */
        if (event_counter == 0 ) { /* only for incoming neutron */
          hkl_info.kix = kix;
          hkl_info.kiy = kiy;
          hkl_info.kiz = kiz;
        }

        hkl_info.coh_refl  = coh_refl;
        hkl_info.coh_xsect = coh_xsect;
        hkl_info.nb_refl += hkl_info.tau_count;
        hkl_info.nb_refl_count++;
      }

      /* (3). Probabilities of the different possible interactions. */
      /* Cross-sections are in barns = 10**-28 m**2, and unit cell volumes are
         in AA**3 = 10**-30 m**2. Hence a factor of 100 is used to convert
         scattering lengths to m**-1 */
      coh_xlen = hkl_info.coh_xsect/V0;
      if (barns) {
        coh_xlen *= 100;
      } /* else assume fm^2 */
      tot_xlen = abs_xlen + inc_xlen + coh_xlen;

      /* (5). Transmission */
      p_trans = exp(-tot_xlen*l_full);
      if(!event_counter && p_transmit >= 0 && p_transmit <= 1) {
        mc_trans = p_transmit; /* first event */
      } else {
        mc_trans = p_trans;
      }
      hkl_info.type = 't';
      mc_interact = 1 - mc_trans;
      if(mc_trans > 0 && (mc_trans >= 1 || rand01() < mc_trans))  /* Transmit */
      {
        p *= p_trans/mc_trans;
        intersect=0;
        if (powder) { /* orientation of crystallite is longer random */
          randderotate(&vx, &vy, &vz, Alpha, Beta, Gamma);
        }
        if (PG) { /* orientation of crystallite is longer random */
          PGderotate(&vx, &vy, &vz, Alpha, hkl_info.m_cx, hkl_info.m_cz, hkl_info.m_cz);
        }
        break;
      }
      if(tot_xlen <= 0){
        ABSORB;
      }
      if(mc_interact <= 0)        /* Protect against rounding errors */
        { intersect=0;
          if (powder) { /* orientation of crystallite is no longer random */
            randderotate(&vx, &vy, &vz, Alpha, Beta, Gamma);
          }
          if (PG) { /* orientation of crystallite is no longer random, rotation around <c> */
            PGderotate(&vx, &vy, &vz, Alpha, hkl_info.csx, hkl_info.csy, hkl_info.csz);
          }
          break;
        }
      if (!event_counter) p *= fabs(1 - p_trans)/mc_interact;
      /* Select a point at which to scatter the neutron, taking
         secondary extinction into account. */
      /* dP(l) = exp(-tot_xlen*l)dl
         P(l<l_0) = [-1/tot_xlen*exp(-tot_xlen*l)]_0^l_0
                  = (1 - exp(-tot_xlen*l0))/tot_xlen
         l = -log(1 - tot_xlen*rand0max(P(l<l_full)))/tot_xlen
       */
      if(tot_xlen*l_full < 1e-6)
        /* For very weak scattering, use simple uniform sampling of scattering
           point to avoid rounding errors. */
        l = rand0max(l_full);
      else
        l = -log(1 - rand0max((1 - exp(-tot_xlen*l_full))))/tot_xlen;
      PROP_DT(l/v);
      event_counter++;

      /* (4). Account for the probability of sigma_abs */
      p *= (coh_xlen + inc_xlen)/tot_xlen;
      /* Choose between coherent and incoherent scattering */
      if(coh_xlen == 0 || rand0max(coh_xlen + inc_xlen) <= inc_xlen)
      {
        /* (6). Incoherent scattering */
        randvec_target_circle(&kix, &kiy, &kiz, NULL, vx, vy, vz, 0);
        vx = kix; /* ki vector is used as tmp var with norm v */
        vy = kiy;
        vz = kiz; /* Go for next scattering event */
        hkl_info.type = 'i';
      } else {
        /* 7. Coherent scattering. Select reciprocal lattice point. */
        if(hkl_info.coh_refl <= 0){
          ABSORB;
        }
        sum = 0;
        j = hkl_select(T, hkl_info.tau_count, hkl_info.coh_refl, &sum);

        if(j >= hkl_info.tau_count)
        {
          if (hkl_info.flag_warning < 100)
            fprintf(stderr, "Single_crystal: Error: Illegal tau search "
              "(r=%g, sum=%g, j=%i, tau_count=%i).\n", r, sum, j , hkl_info.tau_count);
          hkl_info.flag_warning++;
          j = hkl_info.tau_count - 1;
        }
        i = T[j].index;
        /* (8). Pick scattered wavevector kf from 2D Gauss distribution. */
        z1 = randnorm();
        z2 = randnorm();
        y1 = T[j].l11*z1 + T[j].y0x;
        y2 = T[j].l12*z1 + T[j].l22*z2 + T[j].y0y;
        kfx = T[j].rho_x + T[j].ox + T[j].b1x*y1 + T[j].b2x*y2;
        kfy = T[j].rho_y + T[j].oy + T[j].b1y*y1 + T[j].b2y*y2;
        kfz = T[j].rho_z + T[j].oz + T[j].b1z*y1 + T[j].b2z*y2;
        /* Normalize kf to length of ki, to account for planer
          approximation of the Ewald sphere. */
        adjust = ki/sqrt(kfx*kfx + kfy*kfy + kfz*kfz);
        kfx *= adjust;
        kfy *= adjust;
        kfz *= adjust;
        /* Adjust neutron weight (see manual for explanation). */
        p *= T[j].xsect*hkl_info.coh_refl/(hkl_info.coh_xsect*T[j].refl);

        vx = K2V*(L[i].u1x*kfx + L[i].u2x*kfy + L[i].u3x*kfz);
        vy = K2V*(L[i].u1y*kfx + L[i].u2y*kfy + L[i].u3y*kfz);
        vz = K2V*(L[i].u1z*kfx + L[i].u2z*kfy + L[i].u3z*kfz);
        hkl_info.type = 'c';
        hkl_info.h    = L[i].h;
        hkl_info.k    = L[i].k;
        hkl_info.l    = L[i].l;
      }



      /* ------------------------------------------------------------------------- */
      /* lattice curvature option: rotate back neutron velocity */
      if(RX) {
          vec_rotate_2d(&vx,&vz, -curv_xangle);
          vec_rotate_2d(&sx,&sz, -curv_xangle);

          /*changing x,z actually curves the crystal, not only the planes*/
          /*comment out if only curvature of the lattice planes is needed*/
          vec_rotate_2d(&x,&z, -curv_xangle);
      }
      if(RY) {
          vec_rotate_2d(&vy,&vz, -curv_yangle);
          vec_rotate_2d(&sy,&sz, -curv_yangle);

          /*changing y,z actually curves the crystal, not only the planes*/
          /*comment out if only curvature of the lattice planes is needed*/
          vec_rotate_2d(&y,&z, -curv_yangle);
      }
      /* ------------------------------------------------------------------------- */



      SCATTER;
      if (powder) { /* orientation of crystallite is no longer random */
        randderotate(&vx, &vy, &vz, Alpha, Beta, Gamma);
      }
      if (PG) { /* orientation of crystallite is longer random */
        PGderotate(&vx, &vy, &vz, Alpha, hkl_info.csx, hkl_info.csy, hkl_info.csz);
      }
      /* exit if multiple scattering order has been reached */
      if (order && event_counter >= order) { intersect=0; break; }
      /* Repeat loop for next scattering event. */
    } while (intersect); /* end do (intersect) (multiple scattering loop) */
  } /* if intersect */
}
#line 22483 "Samples_Incoherent.c"
}   /* End of SXg=Single_crystal() SETTING parameter declarations. */
#undef offdata
#undef hkl_info
#undef mosaic_AB
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompSXg:
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

  /* TRACE Component Sqwg [10] */
  mccoordschange(mcposrSqwg, mcrotrSqwg,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component Sqwg (without coords transformations) */
  mcJumpTrace_Sqwg:
  SIG_MESSAGE("Sqwg (Trace)");
  mcDEBUG_COMP("Sqwg")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompSqwg
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
#define mccompcurname  Sqwg
#define mccompcurtype  Isotropic_Sqw
#define mccompcurindex 10
#define powder_format mccSqwg_powder_format
#define VarSqw mccSqwg_VarSqw
#define columns mccSqwg_columns
#define offdata mccSqwg_offdata
{   /* Declarations of Sqwg=Isotropic_Sqw() SETTING parameters. */
char* Sqw_coh = mccSqwg_Sqw_coh;
char* Sqw_inc = mccSqwg_Sqw_inc;
char* geometry = mccSqwg_geometry;
MCNUM radius = mccSqwg_radius;
MCNUM thickness = mccSqwg_thickness;
MCNUM xwidth = mccSqwg_xwidth;
MCNUM yheight = mccSqwg_yheight;
MCNUM zdepth = mccSqwg_zdepth;
MCNUM threshold = mccSqwg_threshold;
int order = mccSqwg_order;
MCNUM T = mccSqwg_T;
MCNUM verbose = mccSqwg_verbose;
MCNUM d_phi = mccSqwg_d_phi;
int concentric = mccSqwg_concentric;
MCNUM rho = mccSqwg_rho;
MCNUM sigma_abs = mccSqwg_sigma_abs;
MCNUM sigma_coh = mccSqwg_sigma_coh;
MCNUM sigma_inc = mccSqwg_sigma_inc;
MCNUM classical = mccSqwg_classical;
MCNUM powder_Dd = mccSqwg_powder_Dd;
MCNUM powder_DW = mccSqwg_powder_DW;
MCNUM powder_Vc = mccSqwg_powder_Vc;
MCNUM density = mccSqwg_density;
MCNUM weight = mccSqwg_weight;
MCNUM p_interact = mccSqwg_p_interact;
MCNUM norm = mccSqwg_norm;
MCNUM powder_barns = mccSqwg_powder_barns;
char* quantum_correction = mccSqwg_quantum_correction;
/* 'Sqwg=Isotropic_Sqw()' component instance has conditional execution */
if (( mcipSAMPLE == 8 ))

#line 2176 "/usr/share/mcstas/2.6rc1/samples/Isotropic_Sqw.comp"
{

int    intersect=0;     /* flag to continue/stop */
double t0,  t1,  t2,  t3; /* times for intersections */
double dt0, dt1, dt2, dt; /* time intervals */
double k=0, Ei=0;
double v=0, vf=0;
double d_path;        /* total path length for straight trajectory */
double my_a;          /* absorption cross-section scaled to velocity (2200) */
double ws, p_scatt;   /* probability for scattering/absorption and for */
                      /* interaction along d_path */
double tmp_rand;      /* temporary var */
double ratio_w=0, ratio_q=0; /* variables for bilinear interpolation */
double q11, q21, q22, q12;
double omega=0;       /* energy transfer */
double q=0;           /* wavevector transfer */
long   index_w;       /* energy index for table look-up SW */
long   index_q;       /* Q index for table look-up P(Q|w) */
double theta=0, costheta=0; /* for the choice of kf direction */
double u1x,u1y,u1z;
double u2x,u2y,u2z;
double u0x,u0y,u0z;
int    index_counter;
int    flag=0;
int    flag_concentric=0;
int    flag_ishollow=0;
double solid_angle=0;
double my_t=0;
double p_mult=1;
double mc_trans, p_trans, mc_scatt;
double coh=0, inc=0;
struct Sqw_Data_struct Data_sqw;


/* Store Initial neutron state */

VarSqw.ki_x = V2K*vx;
VarSqw.ki_y = V2K*vy;
VarSqw.ki_z = V2K*vz;
VarSqw.ti   = t;
VarSqw.vi   = 0;
VarSqw.ki   = 0;
VarSqw.type = '\0';

do { /* Main interaction loop. Ends with intersect=0 */

  /* Intersection neutron trajectory / sample (sample surface) */
  if (VarSqw.s_coh > 0 || VarSqw.s_inc > 0) {
    if (thickness >= 0) {
      if (VarSqw.shape==0)
        intersect=cylinder_intersect(&t0,&t3, x,y,z,vx,vy,vz, radius,yheight);
      else if (VarSqw.shape==1)
        intersect=box_intersect     (&t0,&t3, x,y,z,vx,vy,vz, xwidth,yheight,zdepth);
      else if (VarSqw.shape==2)
        intersect=sphere_intersect  (&t0,&t3, x,y,z,vx,vy,vz, radius);
      else if (VarSqw.shape == 3)
        intersect=off_intersect(&t0, &t3, NULL, NULL, x, y, z, vx, vy, vz, offdata );
    } else {
      if (VarSqw.shape==0)
        intersect=cylinder_intersect(&t0,&t3, x,y,z,vx,vy,vz, radius-thickness,
          yheight-2*thickness > 0 ? yheight-2*thickness : yheight);
      else if (VarSqw.shape==1)
        intersect=box_intersect     (&t0,&t3, x,y,z,vx,vy,vz,
          xwidth-2*thickness > 0 ?  xwidth-2*thickness : xwidth,
          yheight-2*thickness > 0 ? yheight-2*thickness : yheight,
          zdepth-2*thickness > 0 ?  zdepth-2*thickness : zdepth);
      else if (VarSqw.shape==2)
        intersect=sphere_intersect  (&t0,&t3, x,y,z,vx,vy,vz, radius-thickness);
      else if (VarSqw.shape == 3)
        intersect=off_intersect(&t0, &t3, NULL, NULL, x, y, z, vx, vy, vz, offdata );
    }
  } else intersect=0;

  /* Computing the intermediate times */
  if (intersect) {
    flag_ishollow = 0;
    if (thickness > 0) {
      if (VarSqw.shape==0 && cylinder_intersect(&t1,&t2, x,y,z,vx,vy,vz, radius-thickness,
        yheight-2*thickness > 0 ? yheight-2*thickness : yheight))
        flag_ishollow=1;
      else if (VarSqw.shape==2 && sphere_intersect   (&t1,&t2, x,y,z,vx,vy,vz, radius-thickness))
        flag_ishollow=1;
      else if (VarSqw.shape==1 && box_intersect(&t1,&t2, x,y,z,vx,vy,vz,
        xwidth-2*thickness > 0 ? xwidth-2*thickness : xwidth,
        yheight-2*thickness > 0 ? yheight-2*thickness : yheight,
        zdepth-2*thickness > 0 ? zdepth-2*thickness : zdepth))
        flag_ishollow=1;
    } else if (thickness<0) {
      if (VarSqw.shape==0 && cylinder_intersect(&t1,&t2, x,y,z,vx,vy,vz, radius,yheight))
        flag_ishollow=1;
      else if (VarSqw.shape==2 && sphere_intersect   (&t1,&t2, x,y,z,vx,vy,vz, radius))
        flag_ishollow=1;
      else if (VarSqw.shape==1 && box_intersect(&t1,&t2, x,y,z,vx,vy,vz, xwidth, yheight, zdepth))
        flag_ishollow=1;
    }
    if (!flag_ishollow) t1 = t2 = t3; /* no empty space inside */
  } else break; /* neutron does not hit sample: transmitted  */

  if (intersect) { /* the neutron hits the sample */

    if (t0 > 0) {  /* we are before the sample */
      PROP_DT(t0); /* propagates neutron to the entry of the sample */
    } else if (t1 > 0 && t1 > t0) { /* we are inside first part of the sample */
      /* no propagation, stay inside */
    } else if (t2 > 0 && t2 > t1) { /* we are in the hole */
      PROP_DT(t2); /* propagate to inner surface of 2nd part of sample */
    } else if (t3 > 0 && t3 > t2) { /* we are in the 2nd part of sample */
      /* no propagation, stay inside */
    }

    dt0=t1-(t0 > 0 ? t0 : 0); /* Time in first part of hollow/cylinder/box */
    dt1=t2-(t1 > 0 ? t1 : 0); /* Time in hole */
    dt2=t3-(t2 > 0 ? t2 : 0); /* Time in 2nd part of hollow cylinder */

    if (dt0 < 0) dt0 = 0;
    if (dt1 < 0) dt1 = 0;
    if (dt2 < 0) dt2 = 0;

    /* initialize concentric mode */
    if (concentric && !flag_concentric && t0 >= 0
     && VarSqw.shape==0 && thickness) {
      flag_concentric=1;
    }

    if (flag_concentric == 1) {
      dt1=dt2=0; /* force exit when reaching hole/2nd part */
    }

    if (!dt0 && !dt2) {
      intersect = 0; /* the sample was passed entirely */
      break;
    }

    VarSqw.neutron_enter++;
    p_mult = 1;
    if (!v) {
      v  = vx*vx+vy*vy+vz*vz;
      v = sqrt(v);
    }
    k  = V2K*v;
    Ei = VS2E*v*v;

    if (!VarSqw.vi) VarSqw.vi = v;
    if (!VarSqw.ki) VarSqw.ki = k;

    if (v <= 0) {
      printf("Isotropic_Sqw: %s: ERROR: Null velocity !\n",NAME_CURRENT_COMP);
      VarSqw.neutron_removed++;
      ABSORB; /* should never occur */
    }

    /* check for scattering event */
    my_a   = VarSqw.my_a_v / v; /* absorption 'mu' */
    /* compute total scattering X section */
    /* \int q S(q) dq /2 /ki^2 sigma  OR  bare Xsection*/
    /* contains the 4*PI*kf/ki factor */
    coh = VarSqw.s_coh;
    inc = VarSqw.s_inc;
    if (k && VarSqw.s_coh>0 && VarSqw.Data_coh.intensity) {
      double Ei       = VS2E*v*v;
      double index_Ei = Ei / (VarSqw.Data_coh.Ei_max/VarSqw.Data_coh.iqSq_length);
      coh = Table_Value2d(VarSqw.Data_coh.iqSq, index_Ei, 0);
    }
    if (k && VarSqw.s_inc>0 && VarSqw.Data_inc.intensity) {
      double Ei       = VS2E*v*v;
      double index_Ei = Ei / (VarSqw.Data_inc.Ei_max/VarSqw.Data_inc.iqSq_length);
      inc = Table_Value2d(VarSqw.Data_inc.iqSq, index_Ei, 0);
    }
    if (coh<0) coh=0;
    if (inc<0) inc=0;
    VarSqw.my_s    =(VarSqw.mat_rho*100*(coh + inc));

    my_t = my_a + VarSqw.my_s;  /* total scattering Xsect */
    if (my_t <= 0) {
      if (VarSqw.neutron_removed<VarSqw.maxloop) printf("Isotropic_Sqw: %s: ERROR: Null total cross section %g. Removing event.\n",
        NAME_CURRENT_COMP, my_t);
      VarSqw.neutron_removed++;
      ABSORB; /* should never occur */
    } else if (VarSqw.my_s <= 0) {
      if (VarSqw.verbose_output > 1 && VarSqw.neutron_removed<VarSqw.maxloop)
        printf("Isotropic_Sqw: %s: Warning: Null scattering cross section %g. Ignoring.\n",
          NAME_CURRENT_COMP, VarSqw.my_s);
      VarSqw.my_s = 0;
    }

    /* Proba of scattering vs absorption (integrating along the whole trajectory) */
    ws = VarSqw.my_s/my_t;  /* (inc+coh)/(inc+coh+abs) */
    d_path = v*( dt0 +dt2 );    /* total path lenght in sample */
    /* Proba of transmission/interaction along length d_path */
    p_trans = exp(-my_t*d_path);
    p_scatt = 1 - p_trans; /* portion of beam which scatters */

    flag = 0; /* flag used for propagation to exit point before ending */

    /* are we next to the exit ? probably no scattering (avoid rounding errors) */
    if (VarSqw.my_s*d_path <= 4e-7) {
      flag = 1;           /* No interaction before the exit */
    }
    /* force a given fraction of the beam to scatter */
    if (p_interact>0 && p_interact<=1) {
      /* we force a portion of the beam to interact */
      /* This is used to improve statistics on single scattering (and multiple) */
      if (!SCATTERED) mc_trans = 1-p_interact;
      else            mc_trans = 1-p_interact/(4*SCATTERED+1); /* reduce effect on multi scatt */
    } else {
      mc_trans = p_trans; /* 1 - p_scatt */
    }
    mc_scatt = 1 - mc_trans; /* portion of beam to scatter (or force to) */
    if (mc_scatt <= 0 || mc_scatt>1) flag=1;
    /* MC choice: Interaction or transmission ? */
    if (!flag && mc_scatt > 0 && (mc_scatt >= 1 || rand01() < mc_scatt)) { /* Interaction neutron/sample */
      p_mult *= ws; /* Update weight ; account for absorption and retain scattered fraction */
      /* we have chosen portion mc_scatt of beam instead of p_scatt, so we compensate */
      if (!mc_scatt) ABSORB;
      p_mult *= fabs(p_scatt/mc_scatt); /* lower than 1 */
    } else {
      flag = 1; /* Transmission : no interaction neutron/sample */
      if (!VarSqw.type) VarSqw.type = 't';
      if (!mc_trans) ABSORB;
      p_mult *= fabs(p_trans/mc_trans);  /* attenuate beam by portion which is scattered (and left along) */
    }

    if (flag) { /* propagate to exit of sample and finish */
      intersect = 0;
      p *= p_mult; /* apply absorption correction */
      PROP_DT(dt0+dt2);
      break; /* exit main multi scatt while loop */
    }
  } /* end if intersect the neutron hits the sample */
  else break;

  if (intersect) { /* scattering event */
    double kf=0, kf1, kf2;
    /* mean scattering probability and absorption fraction */
    VarSqw.mean_scatt += (1-exp(-VarSqw.my_s*d_path))*p;
    VarSqw.mean_abs   += (1-ws)*p;
    VarSqw.psum_scatt += p;

    /* Decaying exponential distribution of the path length before scattering */
    /* Select a point at which to scatter the neutron, taking
         secondary extinction into account. */
    if (my_t*d_path < 1e-6)
    /* For very weak scattering, use simple uniform sampling of scattering
       point to avoid rounding errors. */
      dt = rand0max(d_path); /* length */
    else
      dt = -log(1 - rand0max((1 - exp(-my_t*d_path)))) / my_t; /* length */
    dt /= v; /* Time from present position to scattering point */

    /* If t0 is in hole, propagate to next part of the hollow cylinder */
    if (dt1 > 0 && dt0 > 0 && dt > dt0) dt += dt1;

    /* Neutron propagation to the scattering point */
    PROP_DT(dt);

    /* choice between coherent/incoherent scattering */
    tmp_rand = rand01();
    /* local description at the scattering point (scat probability for atom) */
    tmp_rand *= (coh+inc);

    flag=0;
    if (VarSqw.s_inc>0 && tmp_rand < inc) {
      /* CASE 1: incoherent case */
      if (!VarSqw.Data_inc.intensity) {
        /* CASE 1a: no incoherent Sqw from file, use isotropic V-like */
        if (d_phi && order == 1) {
          randvec_target_rect_angular(&u1x, &u1y, &u1z, &solid_angle,
              vx, vy, vz, 2*PI, d_phi, ROT_A_CURRENT_COMP);
          p_mult *= solid_angle/4/PI; /* weighted by focused range to total range */
        } else
          randvec_target_circle(&u1x, &u1y, &u1z, NULL, vx, vy, vz, 0);

        vx = u1x; vy = u1y; vz = u1z;
        vf = v; kf = k;
        if (!VarSqw.type) VarSqw.type = 'v';
        SCATTER;
      } else {
        /* CASE 1b: incoherent Sqw from file */
        if (VarSqw.Data_inc.intensity) {
          Data_sqw = VarSqw.Data_inc;
          if (!VarSqw.type) VarSqw.type = 'i';
          flag = 1;
        }
      }
    } else if (VarSqw.s_coh>0 && tmp_rand > VarSqw.s_inc) {
      if (VarSqw.Data_coh.intensity) {
        /* CASE2: coherent case */
        Data_sqw = VarSqw.Data_coh;
        if (!VarSqw.type) VarSqw.type = 'c';
        flag = 1;
      }
    }

    if (flag) { /* true when S(q,w) table exists (Data_sqw) */

      double alpha=0, alpha0;
      /* give us a limited number of tries for scattering: choose W then Q */
      for (index_counter=VarSqw.maxloop; index_counter > 0 ; index_counter--) {

        /* MC choice: energy transfer w=Ei-Ef in the S(w) = SW */
        omega = 0;
        tmp_rand = rand01();
        /* energy index for rand > cumul SW */
        index_w  = Sqw_search_SW(Data_sqw, tmp_rand);
        VarSqw.rw = (double)index_w;
        if (index_w >= 0 && &(Data_sqw.SW[index_w]) != NULL) {
          if (Data_sqw.w_bins > 1) {
            double w1, w2;
            if (index_w > 0) { /* interpolate linearly energy */
              ratio_w = (tmp_rand                         - Data_sqw.SW[index_w-1].cumul_proba)
                       /(Data_sqw.SW[index_w].cumul_proba - Data_sqw.SW[index_w-1].cumul_proba);
              /* ratio_w=0 omega[index_w-1], ratio=1 omega[index] */
              w1 = Data_sqw.SW[index_w-1].omega; w2 = Data_sqw.SW[index_w].omega;
            } else { /* index_w = 0 interpolate to 0 energy */
              /* ratio_w=0 omega=0, ratio=1 omega[index] */
              w1 = Data_sqw.SW[index_w].omega; w2= Data_sqw.SW[index_w+1].omega;
              if (!w2 && index_w+1 < Data_sqw.w_bins)
                w2= Data_sqw.SW[index_w+1].omega;
              if (Data_sqw.w_bins && Data_sqw.SW[index_w].cumul_proba) {
                ratio_w = tmp_rand/Data_sqw.SW[index_w].cumul_proba;
              } else ratio_w=0;
            }
            if (ratio_w<0) ratio_w=0; else if (ratio_w>1) ratio_w=1;
            omega = (1-ratio_w)*w1 + ratio_w*w2;
          } else {
            ratio_w = 0;
            omega = Data_sqw.SW[index_w].omega;
          }
        } else {
          if (VarSqw.verbose_output >= 3 && VarSqw.neutron_removed<VarSqw.maxloop)
            printf("Isotropic_Sqw: %s: Warning: No suitable w transfer for index_w=%li.\n",
              NAME_CURRENT_COMP, index_w);
          continue; /* no W value: try again with an other energy transfer */
        }

        /* MC choice: momentum transfer Q in P(Q|w) */
        tmp_rand = rand01();

        /* momentum index for rand > cumul SQ|W */
        index_q  = Sqw_search_Q_proba_per_w(Data_sqw, tmp_rand, index_w);
        VarSqw.rq = (double)index_q;

        if (index_q >= 0 && &(Data_sqw.SQW[index_w]) != NULL) {
          if (Data_sqw.q_bins > 1 && index_q > 0) {
            if (index_w > 0 && Data_sqw.w_bins > 1) {
              /* bilinear interpolation on - side: index_w > 0, index_q > 0 */
              ratio_q = (tmp_rand - Data_sqw.SQW[index_w][index_q-1].cumul_proba)
                       /(Data_sqw.SQW[index_w][index_q].cumul_proba
                       - Data_sqw.SQW[index_w][index_q-1].cumul_proba);
              q22 = Data_sqw.SQW[index_w]  [index_q].Q;
              q11 = Data_sqw.SQW[index_w-1][index_q-1].Q;
              q21 = Data_sqw.SQW[index_w]  [index_q-1].Q;
              q12 = Data_sqw.SQW[index_w-1][index_q].Q;
              if (ratio_q<0) ratio_q=0; else if (ratio_q>1) ratio_q=1;
              q = (1-ratio_w)*(1-ratio_q)*q11+ratio_w*(1-ratio_q)*q21
                + ratio_w*ratio_q*q22        +(1-ratio_w)*ratio_q*q12;
            } else { /* bilinear interpolation on + side: index_w=0, index_q > 0 */
              ratio_q = (tmp_rand - Data_sqw.SQW[index_w][index_q-1].cumul_proba)
                       /(Data_sqw.SQW[index_w][index_q].cumul_proba
                       - Data_sqw.SQW[index_w][index_q-1].cumul_proba);
              q11 = Data_sqw.SQW[index_w]  [index_q-1].Q;
              q12 = Data_sqw.SQW[index_w]  [index_q].Q;
              if (ratio_q<0) ratio_q=0; else if (ratio_q>1) ratio_q=1;
              if (index_w < Data_sqw.w_bins-1 && Data_sqw.w_bins > 1) {
                q22 = Data_sqw.SQW[index_w+1][index_q].Q;
                q21 = Data_sqw.SQW[index_w+1][index_q-1].Q;
                q = (1-ratio_w)*(1-ratio_q)*q11+ratio_w*(1-ratio_q)*q21
                  + ratio_w*ratio_q*q22        +(1-ratio_w)*ratio_q*q12;
              } else {
                q    = (1-ratio_q)*q11  + ratio_q*q12;
              }
            }
          } else {
            q    = Data_sqw.SQW[index_w][index_q].Q;
          }
        } else {
          if (VarSqw.verbose_output >= 3 && VarSqw.neutron_removed<VarSqw.maxloop)
            printf("Isotropic_Sqw: %s: Warning: No suitable q transfer for w=%g.\n",
              NAME_CURRENT_COMP, omega);
          VarSqw.neutron_removed++;
          continue; /* no Q value for this w choice */
        }

        /* Search for length of final wave vector kf */
        /* kf is such that : hbar*w = hbar*hbar/2/m*(k*k - kf*kf) */
        /* acceptable values for kf are kf1 and kf2 */
        if (!solve_2nd_order(&kf1, &kf2, 1, 0, -k*k + VarSqw.sqSE2K*omega)) {
          if (VarSqw.verbose_output >= 3 && VarSqw.neutron_removed<VarSqw.maxloop)
            printf("Isotropic_Sqw: %s: Warning: imaginary root for w=%g q=%g Ei=%g (triangle can not close)\n",
            NAME_CURRENT_COMP, omega, q, Ei);
          VarSqw.neutron_removed++;
          continue; /* all roots are imaginary */
        }

        /* kf1 and kf2 are opposite */
        kf = fabs(kf1);
        vf = K2V*kf;

        /* Search of the direction of kf such that : q = ki - kf */
        /* cos theta = (ki2+kf2-q2)/(2ki kf) */

        costheta= (k*k+kf*kf-q*q)/(2*kf*k); /* this is cos(theta) */

        if (-1 < costheta && costheta < 1) {
          break; /* satisfies q momentum conservation */
        }
/*      else continue; */

        /* exit for loop on success */
      } /* end for index_counter */

      if (!index_counter) { /* for loop ended: failure for scattering */
        intersect=0; /* Could not scatter: finish multiple scattering loop */
        if (VarSqw.verbose_output >= 2 && VarSqw.neutron_removed<VarSqw.maxloop)
          printf("Isotropic_Sqw: %s: Warning: No scattering [q,w] conditions\n"
               "               last try (%i): type=%c w=%g q=%g cos(theta)=%g k=%g\n",
          NAME_CURRENT_COMP, VarSqw.maxloop, (VarSqw.type ? VarSqw.type : '-'), omega, q, costheta, k);
        VarSqw.neutron_removed++;
        if (order && SCATTERED != order) ABSORB;
        break;       /* finish multiple scattering loop */
      }

      /* scattering angle from ki to DS cone */
      theta = acos(costheta);

      /* Choose point on Debye-Scherrer cone */
      if (order == 1 && d_phi)
      { /* relate height of detector to the height on DS cone */
        double cone_focus;
        cone_focus = sin(d_phi/2)/sin(theta);
        /* If full Debye-Scherrer cone is within d_phi, don't focus */
        if (cone_focus < -1 || cone_focus > 1) d_phi = 0;
        /* Otherwise, determine alpha to rotate from scattering plane
            into d_phi focusing area*/
        else alpha = 2*asin(cone_focus);
        if (d_phi) p_mult *= alpha/PI;
      }
      if (d_phi) {
        /* Focusing */
        alpha = fabs(alpha);
        /* Trick to get scattering for pos/neg theta's */
        alpha0= 2*rand01()*alpha;
        if (alpha0 > alpha) {
          alpha0=PI+(alpha0-1.5*alpha);
        } else {
          alpha0=alpha0-0.5*alpha;
        }
      }
      else
        alpha0 = PI*randpm1();

      /* now find a nearly vertical rotation axis (u1) :
	       * Either
	       *  (v along Z) x (X axis) -> nearly Y axis
	       * Or
	       *  (v along X) x (Z axis) -> nearly Y axis
	       */
	    if (fabs(scalar_prod(1,0,0,vx/v,vy/v,vz/v)) < fabs(scalar_prod(0,0,1,vx/v,vy/v,vz/v))) {
        u1x = 1; u1y = u1z = 0;
	    } else {
        u1x = u1y = 0; u1z = 1;
	    }
	    vec_prod(u2x,u2y,u2z, vx,vy,vz, u1x,u1y,u1z);

      /* handle case where v and aim are parallel */
      if (!u2x && !u2y && !u2z) { u2x=u2z=0; u2y=1; }

      /* u1 = rotate 'v' by theta around u2: DS scattering angle, nearly in horz plane */
      rotate(u1x,u1y,u1z, vx,vy,vz, theta, u2x,u2y,u2z);

      /* u0 = rotate u1 by alpha0 around 'v' (Debye-Scherrer cone) */
      rotate(u0x,u0y,u0z, u1x,u1y,u1z, alpha0, vx, vy, vz);
      NORM(u0x,u0y,u0z);
      vx = u0x*vf;
      vy = u0y*vf;
      vz = u0z*vf;

      SCATTER;

      v = vf; k = kf; /* for next iteration */

    } /* end if (flag) */

    VarSqw.neutron_exit++;
    p *= p_mult;
    if (p_mult > 1) VarSqw.neutron_pmult++;

    /* test for a given multiple order */
    if (order && SCATTERED >= order) {
      intersect=0; /* reached required number of SCATTERing */
      break;       /* finish multiple scattering loop */
    }

  } /* end if (intersect) scattering event  */

} while (intersect); /* end do (intersect) (multiple scattering loop) */

/* Store Final neutron state */
VarSqw.kf_x = V2K*vx;
VarSqw.kf_y = V2K*vy;
VarSqw.kf_z = V2K*vz;
VarSqw.tf   = t;
VarSqw.vf   = v;
VarSqw.kf   = k;
VarSqw.theta= theta;

if (SCATTERED) {



  if (SCATTERED == 1) {
    if (VarSqw.type == 'c') VarSqw.single_coh += p;
    else                    VarSqw.single_inc += p;
    VarSqw.dq = sqrt((VarSqw.kf_x-VarSqw.ki_x)*(VarSqw.kf_x-VarSqw.ki_x)
                  +(VarSqw.kf_y-VarSqw.ki_y)*(VarSqw.kf_y-VarSqw.ki_y)
                  +(VarSqw.kf_z-VarSqw.ki_z)*(VarSqw.kf_z-VarSqw.ki_z));
    VarSqw.dw = VS2E*(VarSqw.vf*VarSqw.vf - VarSqw.vi*VarSqw.vi);
  } else VarSqw.multi += p;

} else VarSqw.dq=VarSqw.dw=0;

/* end TRACE */
}
#line 23150 "Samples_Incoherent.c"
}   /* End of Sqwg=Isotropic_Sqw() SETTING parameter declarations. */
#undef offdata
#undef columns
#undef VarSqw
#undef powder_format
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompSqwg:
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

  /* TRACE Component IncSg [11] */
  mccoordschange(mcposrIncSg, mcrotrIncSg,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component IncSg (without coords transformations) */
  mcJumpTrace_IncSg:
  SIG_MESSAGE("IncSg (Trace)");
  mcDEBUG_COMP("IncSg")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompIncSg
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
#define mccompcurname  IncSg
#define mccompcurtype  Incoherent
#define mccompcurindex 11
#define VarsInc mccIncSg_VarsInc
#define offdata mccIncSg_offdata
{   /* Declarations of IncSg=Incoherent() SETTING parameters. */
char* geometry = mccIncSg_geometry;
MCNUM radius = mccIncSg_radius;
MCNUM xwidth = mccIncSg_xwidth;
MCNUM yheight = mccIncSg_yheight;
MCNUM zdepth = mccIncSg_zdepth;
MCNUM thickness = mccIncSg_thickness;
MCNUM target_x = mccIncSg_target_x;
MCNUM target_y = mccIncSg_target_y;
MCNUM target_z = mccIncSg_target_z;
MCNUM focus_r = mccIncSg_focus_r;
MCNUM focus_xw = mccIncSg_focus_xw;
MCNUM focus_yh = mccIncSg_focus_yh;
MCNUM focus_aw = mccIncSg_focus_aw;
MCNUM focus_ah = mccIncSg_focus_ah;
int target_index = mccIncSg_target_index;
MCNUM pack = mccIncSg_pack;
MCNUM p_interact = mccIncSg_p_interact;
MCNUM f_QE = mccIncSg_f_QE;
MCNUM gamma = mccIncSg_gamma;
MCNUM Etrans = mccIncSg_Etrans;
MCNUM deltaE = mccIncSg_deltaE;
MCNUM sigma_abs = mccIncSg_sigma_abs;
MCNUM sigma_inc = mccIncSg_sigma_inc;
MCNUM Vc = mccIncSg_Vc;
MCNUM concentric = mccIncSg_concentric;
MCNUM order = mccIncSg_order;
/* 'IncSg=Incoherent()' component instance has conditional execution */
if (( mcipSAMPLE == 9 ))

#line 232 "/usr/share/mcstas/2.6rc1/samples/Incoherent.comp"
{
  double t0, t3;                /* Entry/exit time for outer surface */
  double t1, t2;                /* Entry/exit time for inner surface */
  double dt0, dt1, dt2, dt;     /* Flight times through sample */
  double v=0;                   /* Neutron velocity */
  double d_path;                /* Flight path length for non-scattered neutron */
  double l_i, l_o=0;            /* Flight path lenght in/out for scattered neutron */
  double my_a=0,my_t=0;         /* Velocity-dependent attenuation factor and total Xsec */
  double solid_angle=0;         /* Solid angle of target as seen from scattering point */
  double aim_x=0, aim_y=0, aim_z=1;   /* Position of target relative to scattering point */
  double v_i, v_f, E_i, E_f; /* initial and final energies and velocities */
  double dE;                 /* Energy transfer */
  int    intersect=0;
  int    flag_concentric=0;
  int    flag=0;
  double mc_trans, p_trans, mc_scatt, p_scatt, ws;
  double p_mult=1;

  do { /* Main interaction loop. Ends with intersect=0 */

    /* Intersection neutron trajectory / sample (sample surface) */
    if (VarsInc.shape == 0)
      intersect = cylinder_intersect(&t0, &t3, x, y, z, vx, vy, vz, radius, yheight);
    else if (VarsInc.shape == 1)
      intersect = box_intersect(&t0, &t3, x, y, z, vx, vy, vz, xwidth, yheight, zdepth);
    else if (VarsInc.shape == 2)
      intersect = sphere_intersect(&t0, &t3, x, y, z, vx, vy, vz, radius);
    else if (VarsInc.shape == 3)
      intersect = off_intersect(&t0, &t3, NULL, NULL, x, y, z, vx, vy, vz, offdata );
    if (intersect) {
      int flag_ishollow = 0;
      if (thickness>0) {
        if (VarsInc.shape==0 && cylinder_intersect(&t1,&t2, x,y,z,vx,vy,vz, radius-thickness,yheight-2*thickness))
          flag_ishollow=1;
        else if (VarsInc.shape==2 && sphere_intersect   (&t1,&t2, x,y,z,vx,vy,vz, radius-thickness))
          flag_ishollow=1;
        else if (VarsInc.shape==1 && box_intersect(&t1,&t2, x,y,z,vx,vy,vz, xwidth-2*thickness, yheight-2*thickness, zdepth-2*thickness))
          flag_ishollow = 1;
      }
      if (!flag_ishollow) t1 = t2 = t3; /* no empty space inside */

      dt0 = t1-t0;                /* Time in sample, ingoing */
      dt1 = t2-t1;                /* Time in hole */
      dt2 = t3-t2;                /* Time in sample, outgoing */

      if (t0 > 0) {  /* we are before the sample */
        PROP_DT(t0); /* propagates neutron to the entry of the sample */
      } else if (t1 > 0 && t1 > t0) { /* we are inside first part of the sample */
        /* no propagation, stay inside */
      } else if (t2 > 0 && t2 > t1) { /* we are in the hole */
        PROP_DT(t2); /* propagate to inner surface of 2nd part of sample */
      } else if (t3 > 0 && t3 > t2) { /* we are in the 2nd part of sample */
        /* no propagation, stay inside */
      }

      dt0=t1-(t0 > 0 ? t0 : 0); /* Time in first part of hollow/cylinder/box */
      dt1=t2-(t1 > 0 ? t1 : 0); /* Time in hole */
      dt2=t3-(t2 > 0 ? t2 : 0); /* Time in 2nd part of hollow cylinder */

      if (dt0 < 0) dt0 = 0;
      if (dt1 < 0) dt1 = 0;
      if (dt2 < 0) dt2 = 0;

      /* initialize concentric mode */
      if (concentric && !flag_concentric && t0 >= 0
       && VarsInc.shape==0 && thickness>0) {
        flag_concentric=1;
      }

      if (flag_concentric == 1) {
        dt1=dt2=0; /* force exit when reaching hole/2nd part */
      }

      if (!dt0 && !dt2) {
        intersect = 0; /* the sample was passed entirely */
        break;
      }

      p_mult = 1;
      if (!v) v = sqrt(vx*vx + vy*vy + vz*vz);
      if (v) my_a = VarsInc.my_a_v*(2200/v);
      else {
        printf("Incoherent: %s: ERROR: Null velocity\n",NAME_CURRENT_COMP);
        ABSORB; /* should never occur */
      }

      my_t = my_a + VarsInc.my_s;  /* total scattering Xsect (tmp var) */
      if (my_t <= 0) {
        printf("Incoherent: %s: ERROR: Null total cross section %g. Removing event.\n",
          NAME_CURRENT_COMP, my_t);
        ABSORB; /* should never occur */
      }
      d_path = v * (dt0 + dt2);   /* Length of full path through sample */
      /* Proba of scattering vs absorption (integrating along the whole trajectory) */
      ws = VarsInc.my_s/my_t;  /* (inc+coh)/(inc+coh+abs) */
      /* Proba of transmission along length d_path */
      p_trans = exp(-my_t*d_path);
      p_scatt = 1 - p_trans; /* portion of beam which scatters */
      flag = 0; /* flag used for propagation to exit point before ending */
      /* are we next to the exit ? probably no scattering (avoid rounding errors) */
      if (VarsInc.my_s*d_path <= 4e-7) {
        flag = 1;           /* No interaction before the exit */
      }
      /* force a given fraction of the beam to scatter */
      if (p_interact>0 && p_interact<=1) {
        /* we force a portion of the beam to interact */
        /* This is used to improve statistics on single scattering (and multiple) */
        if (!SCATTERED) mc_trans = 1-p_interact;
        else            mc_trans = 1-p_interact/(4*SCATTERED+1); /* reduce effect on multi scatt */
      } else {
        mc_trans = p_trans; /* 1 - p_scatt */
      }
      mc_scatt = 1 - mc_trans; /* portion of beam to scatter (or force to) */
      if (mc_scatt <= 0 || mc_scatt>1) flag=1;
      /* MC choice: Interaction or transmission ? */
      if (!flag && mc_scatt > 0 && (mc_scatt >= 1 || (rand01()) < mc_scatt)) { /* Interaction neutron/sample */
        p_mult *= ws; /* Update weight ; account for absorption and retain scattered fraction */
        if (!mc_scatt) ABSORB;
        /* we have chosen portion mc_scatt of beam instead of p_scatt, so we compensate */
        p_mult *= fabs(p_scatt/mc_scatt); /* lower than 1 */
      } else {
        flag = 1; /* Transmission : no interaction neutron/sample */
        if (!mc_trans) ABSORB;
        p_mult *= fabs(p_trans/mc_trans);  /* attenuate beam by portion which is scattered (and left along) */
      }

      if (flag) { /* propagate to exit of sample and finish */
        intersect = 0;
        p *= p_mult; /* apply absorption correction */
        PROP_DT(dt0+dt2);
        break; /* exit main multi scatt while loop */
      }
      if (my_t*d_path < 1e-6)
      /* For very weak scattering, use simple uniform sampling of scattering
         point to avoid rounding errors. */
        dt = rand0max(d_path); /* length */
      else
        dt = -log(1 - rand0max((1 - exp(-my_t*d_path)))) / my_t; /* length */
      l_i = dt;/* Penetration in sample: scattering+abs */
      dt /= v; /* Time from present position to scattering point */

      /* If t0 is in hole, propagate to next part of the hollow cylinder */
      if (dt1 > 0 && dt0 > 0 && dt > dt0) dt += dt1;

      PROP_DT(dt);                /* Point of scattering */

      if ((VarsInc.tx || VarsInc.ty || VarsInc.tz)) {
        aim_x = VarsInc.tx-x;       /* Vector pointing at target (anal./det.) */
        aim_y = VarsInc.ty-y;
        aim_z = VarsInc.tz-z;
      }
      if(VarsInc.aw && VarsInc.ah) {
        randvec_target_rect_angular(&vx, &vy, &vz, &solid_angle,
          aim_x, aim_y, aim_z, VarsInc.aw, VarsInc.ah, ROT_A_CURRENT_COMP);
      } else if(VarsInc.xw && VarsInc.yh) {
        randvec_target_rect(&vx, &vy, &vz, &solid_angle,
          aim_x, aim_y, aim_z, VarsInc.xw, VarsInc.yh, ROT_A_CURRENT_COMP);
      } else {
        randvec_target_circle(&vx, &vy, &vz, &solid_angle, aim_x, aim_y, aim_z, focus_r);
      }
      NORM(vx, vy, vz);

      v_i = v;          /* Store initial velocity in case of quasielastic */
      E_i = VS2E*v_i*v_i;
      if (deltaE==0) { 
	if (rand01()<f_QE)	/* Quasielastic contribution */
	  {
	    dE = gamma*tan(PI/2*randpm1());
	    E_f = E_i + dE;
	    if (E_f <= 0)
	      ABSORB;
	    v_f = SE2V*sqrt(E_f);
	    v = v_f;
	    /*          printf("vi: %g Ei: %g dE: %g Ef %g vf: %g v: %g \n",
			v_i,E_i,dE,E_f,v_f,v); */
	  }
      } else {
	E_f = E_i - Etrans + deltaE*randpm1(); // E_f=E0;
	v_f = SE2V*sqrt(E_f);
	v = v_f;
      }

      vx *= v;
      vy *= v;
      vz *= v;

      /* We do not consider scattering from 2nd part (outgoing) */
      p_mult *= solid_angle/4/PI;
      p *= p_mult;

      /* Polarisation part (1/3 NSF, 2/3 SF) */
      sx *= -1.0/3.0;
      sy *= -1.0/3.0;
      sz *= -1.0/3.0;

      SCATTER;

      /* test for a given multiple order */
      if (order && SCATTERED >= order) {
        intersect=0; /* reached required number of SCATTERing */
        break;       /* finish multiple scattering loop */
      }
    } /* end if intersect */
  } while (intersect); /* end do (intersect) (multiple scattering loop) */
}
#line 23496 "Samples_Incoherent.c"
}   /* End of IncSg=Incoherent() SETTING parameter declarations. */
#undef offdata
#undef VarsInc
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompIncSg:
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

  /* TRACE Component SqwV [12] */
  mccoordschange(mcposrSqwV, mcrotrSqwV,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component SqwV (without coords transformations) */
  mcJumpTrace_SqwV:
  SIG_MESSAGE("SqwV (Trace)");
  mcDEBUG_COMP("SqwV")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompSqwV
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
#define mccompcurname  SqwV
#define mccompcurtype  Isotropic_Sqw
#define mccompcurindex 12
#define powder_format mccSqwV_powder_format
#define VarSqw mccSqwV_VarSqw
#define columns mccSqwV_columns
#define offdata mccSqwV_offdata
{   /* Declarations of SqwV=Isotropic_Sqw() SETTING parameters. */
char* Sqw_coh = mccSqwV_Sqw_coh;
char* Sqw_inc = mccSqwV_Sqw_inc;
char* geometry = mccSqwV_geometry;
MCNUM radius = mccSqwV_radius;
MCNUM thickness = mccSqwV_thickness;
MCNUM xwidth = mccSqwV_xwidth;
MCNUM yheight = mccSqwV_yheight;
MCNUM zdepth = mccSqwV_zdepth;
MCNUM threshold = mccSqwV_threshold;
int order = mccSqwV_order;
MCNUM T = mccSqwV_T;
MCNUM verbose = mccSqwV_verbose;
MCNUM d_phi = mccSqwV_d_phi;
int concentric = mccSqwV_concentric;
MCNUM rho = mccSqwV_rho;
MCNUM sigma_abs = mccSqwV_sigma_abs;
MCNUM sigma_coh = mccSqwV_sigma_coh;
MCNUM sigma_inc = mccSqwV_sigma_inc;
MCNUM classical = mccSqwV_classical;
MCNUM powder_Dd = mccSqwV_powder_Dd;
MCNUM powder_DW = mccSqwV_powder_DW;
MCNUM powder_Vc = mccSqwV_powder_Vc;
MCNUM density = mccSqwV_density;
MCNUM weight = mccSqwV_weight;
MCNUM p_interact = mccSqwV_p_interact;
MCNUM norm = mccSqwV_norm;
MCNUM powder_barns = mccSqwV_powder_barns;
char* quantum_correction = mccSqwV_quantum_correction;
/* 'SqwV=Isotropic_Sqw()' component instance has conditional execution */
if (( mcipSAMPLE == 10 ))

#line 2176 "/usr/share/mcstas/2.6rc1/samples/Isotropic_Sqw.comp"
{

int    intersect=0;     /* flag to continue/stop */
double t0,  t1,  t2,  t3; /* times for intersections */
double dt0, dt1, dt2, dt; /* time intervals */
double k=0, Ei=0;
double v=0, vf=0;
double d_path;        /* total path length for straight trajectory */
double my_a;          /* absorption cross-section scaled to velocity (2200) */
double ws, p_scatt;   /* probability for scattering/absorption and for */
                      /* interaction along d_path */
double tmp_rand;      /* temporary var */
double ratio_w=0, ratio_q=0; /* variables for bilinear interpolation */
double q11, q21, q22, q12;
double omega=0;       /* energy transfer */
double q=0;           /* wavevector transfer */
long   index_w;       /* energy index for table look-up SW */
long   index_q;       /* Q index for table look-up P(Q|w) */
double theta=0, costheta=0; /* for the choice of kf direction */
double u1x,u1y,u1z;
double u2x,u2y,u2z;
double u0x,u0y,u0z;
int    index_counter;
int    flag=0;
int    flag_concentric=0;
int    flag_ishollow=0;
double solid_angle=0;
double my_t=0;
double p_mult=1;
double mc_trans, p_trans, mc_scatt;
double coh=0, inc=0;
struct Sqw_Data_struct Data_sqw;


/* Store Initial neutron state */

VarSqw.ki_x = V2K*vx;
VarSqw.ki_y = V2K*vy;
VarSqw.ki_z = V2K*vz;
VarSqw.ti   = t;
VarSqw.vi   = 0;
VarSqw.ki   = 0;
VarSqw.type = '\0';

do { /* Main interaction loop. Ends with intersect=0 */

  /* Intersection neutron trajectory / sample (sample surface) */
  if (VarSqw.s_coh > 0 || VarSqw.s_inc > 0) {
    if (thickness >= 0) {
      if (VarSqw.shape==0)
        intersect=cylinder_intersect(&t0,&t3, x,y,z,vx,vy,vz, radius,yheight);
      else if (VarSqw.shape==1)
        intersect=box_intersect     (&t0,&t3, x,y,z,vx,vy,vz, xwidth,yheight,zdepth);
      else if (VarSqw.shape==2)
        intersect=sphere_intersect  (&t0,&t3, x,y,z,vx,vy,vz, radius);
      else if (VarSqw.shape == 3)
        intersect=off_intersect(&t0, &t3, NULL, NULL, x, y, z, vx, vy, vz, offdata );
    } else {
      if (VarSqw.shape==0)
        intersect=cylinder_intersect(&t0,&t3, x,y,z,vx,vy,vz, radius-thickness,
          yheight-2*thickness > 0 ? yheight-2*thickness : yheight);
      else if (VarSqw.shape==1)
        intersect=box_intersect     (&t0,&t3, x,y,z,vx,vy,vz,
          xwidth-2*thickness > 0 ?  xwidth-2*thickness : xwidth,
          yheight-2*thickness > 0 ? yheight-2*thickness : yheight,
          zdepth-2*thickness > 0 ?  zdepth-2*thickness : zdepth);
      else if (VarSqw.shape==2)
        intersect=sphere_intersect  (&t0,&t3, x,y,z,vx,vy,vz, radius-thickness);
      else if (VarSqw.shape == 3)
        intersect=off_intersect(&t0, &t3, NULL, NULL, x, y, z, vx, vy, vz, offdata );
    }
  } else intersect=0;

  /* Computing the intermediate times */
  if (intersect) {
    flag_ishollow = 0;
    if (thickness > 0) {
      if (VarSqw.shape==0 && cylinder_intersect(&t1,&t2, x,y,z,vx,vy,vz, radius-thickness,
        yheight-2*thickness > 0 ? yheight-2*thickness : yheight))
        flag_ishollow=1;
      else if (VarSqw.shape==2 && sphere_intersect   (&t1,&t2, x,y,z,vx,vy,vz, radius-thickness))
        flag_ishollow=1;
      else if (VarSqw.shape==1 && box_intersect(&t1,&t2, x,y,z,vx,vy,vz,
        xwidth-2*thickness > 0 ? xwidth-2*thickness : xwidth,
        yheight-2*thickness > 0 ? yheight-2*thickness : yheight,
        zdepth-2*thickness > 0 ? zdepth-2*thickness : zdepth))
        flag_ishollow=1;
    } else if (thickness<0) {
      if (VarSqw.shape==0 && cylinder_intersect(&t1,&t2, x,y,z,vx,vy,vz, radius,yheight))
        flag_ishollow=1;
      else if (VarSqw.shape==2 && sphere_intersect   (&t1,&t2, x,y,z,vx,vy,vz, radius))
        flag_ishollow=1;
      else if (VarSqw.shape==1 && box_intersect(&t1,&t2, x,y,z,vx,vy,vz, xwidth, yheight, zdepth))
        flag_ishollow=1;
    }
    if (!flag_ishollow) t1 = t2 = t3; /* no empty space inside */
  } else break; /* neutron does not hit sample: transmitted  */

  if (intersect) { /* the neutron hits the sample */

    if (t0 > 0) {  /* we are before the sample */
      PROP_DT(t0); /* propagates neutron to the entry of the sample */
    } else if (t1 > 0 && t1 > t0) { /* we are inside first part of the sample */
      /* no propagation, stay inside */
    } else if (t2 > 0 && t2 > t1) { /* we are in the hole */
      PROP_DT(t2); /* propagate to inner surface of 2nd part of sample */
    } else if (t3 > 0 && t3 > t2) { /* we are in the 2nd part of sample */
      /* no propagation, stay inside */
    }

    dt0=t1-(t0 > 0 ? t0 : 0); /* Time in first part of hollow/cylinder/box */
    dt1=t2-(t1 > 0 ? t1 : 0); /* Time in hole */
    dt2=t3-(t2 > 0 ? t2 : 0); /* Time in 2nd part of hollow cylinder */

    if (dt0 < 0) dt0 = 0;
    if (dt1 < 0) dt1 = 0;
    if (dt2 < 0) dt2 = 0;

    /* initialize concentric mode */
    if (concentric && !flag_concentric && t0 >= 0
     && VarSqw.shape==0 && thickness) {
      flag_concentric=1;
    }

    if (flag_concentric == 1) {
      dt1=dt2=0; /* force exit when reaching hole/2nd part */
    }

    if (!dt0 && !dt2) {
      intersect = 0; /* the sample was passed entirely */
      break;
    }

    VarSqw.neutron_enter++;
    p_mult = 1;
    if (!v) {
      v  = vx*vx+vy*vy+vz*vz;
      v = sqrt(v);
    }
    k  = V2K*v;
    Ei = VS2E*v*v;

    if (!VarSqw.vi) VarSqw.vi = v;
    if (!VarSqw.ki) VarSqw.ki = k;

    if (v <= 0) {
      printf("Isotropic_Sqw: %s: ERROR: Null velocity !\n",NAME_CURRENT_COMP);
      VarSqw.neutron_removed++;
      ABSORB; /* should never occur */
    }

    /* check for scattering event */
    my_a   = VarSqw.my_a_v / v; /* absorption 'mu' */
    /* compute total scattering X section */
    /* \int q S(q) dq /2 /ki^2 sigma  OR  bare Xsection*/
    /* contains the 4*PI*kf/ki factor */
    coh = VarSqw.s_coh;
    inc = VarSqw.s_inc;
    if (k && VarSqw.s_coh>0 && VarSqw.Data_coh.intensity) {
      double Ei       = VS2E*v*v;
      double index_Ei = Ei / (VarSqw.Data_coh.Ei_max/VarSqw.Data_coh.iqSq_length);
      coh = Table_Value2d(VarSqw.Data_coh.iqSq, index_Ei, 0);
    }
    if (k && VarSqw.s_inc>0 && VarSqw.Data_inc.intensity) {
      double Ei       = VS2E*v*v;
      double index_Ei = Ei / (VarSqw.Data_inc.Ei_max/VarSqw.Data_inc.iqSq_length);
      inc = Table_Value2d(VarSqw.Data_inc.iqSq, index_Ei, 0);
    }
    if (coh<0) coh=0;
    if (inc<0) inc=0;
    VarSqw.my_s    =(VarSqw.mat_rho*100*(coh + inc));

    my_t = my_a + VarSqw.my_s;  /* total scattering Xsect */
    if (my_t <= 0) {
      if (VarSqw.neutron_removed<VarSqw.maxloop) printf("Isotropic_Sqw: %s: ERROR: Null total cross section %g. Removing event.\n",
        NAME_CURRENT_COMP, my_t);
      VarSqw.neutron_removed++;
      ABSORB; /* should never occur */
    } else if (VarSqw.my_s <= 0) {
      if (VarSqw.verbose_output > 1 && VarSqw.neutron_removed<VarSqw.maxloop)
        printf("Isotropic_Sqw: %s: Warning: Null scattering cross section %g. Ignoring.\n",
          NAME_CURRENT_COMP, VarSqw.my_s);
      VarSqw.my_s = 0;
    }

    /* Proba of scattering vs absorption (integrating along the whole trajectory) */
    ws = VarSqw.my_s/my_t;  /* (inc+coh)/(inc+coh+abs) */
    d_path = v*( dt0 +dt2 );    /* total path lenght in sample */
    /* Proba of transmission/interaction along length d_path */
    p_trans = exp(-my_t*d_path);
    p_scatt = 1 - p_trans; /* portion of beam which scatters */

    flag = 0; /* flag used for propagation to exit point before ending */

    /* are we next to the exit ? probably no scattering (avoid rounding errors) */
    if (VarSqw.my_s*d_path <= 4e-7) {
      flag = 1;           /* No interaction before the exit */
    }
    /* force a given fraction of the beam to scatter */
    if (p_interact>0 && p_interact<=1) {
      /* we force a portion of the beam to interact */
      /* This is used to improve statistics on single scattering (and multiple) */
      if (!SCATTERED) mc_trans = 1-p_interact;
      else            mc_trans = 1-p_interact/(4*SCATTERED+1); /* reduce effect on multi scatt */
    } else {
      mc_trans = p_trans; /* 1 - p_scatt */
    }
    mc_scatt = 1 - mc_trans; /* portion of beam to scatter (or force to) */
    if (mc_scatt <= 0 || mc_scatt>1) flag=1;
    /* MC choice: Interaction or transmission ? */
    if (!flag && mc_scatt > 0 && (mc_scatt >= 1 || rand01() < mc_scatt)) { /* Interaction neutron/sample */
      p_mult *= ws; /* Update weight ; account for absorption and retain scattered fraction */
      /* we have chosen portion mc_scatt of beam instead of p_scatt, so we compensate */
      if (!mc_scatt) ABSORB;
      p_mult *= fabs(p_scatt/mc_scatt); /* lower than 1 */
    } else {
      flag = 1; /* Transmission : no interaction neutron/sample */
      if (!VarSqw.type) VarSqw.type = 't';
      if (!mc_trans) ABSORB;
      p_mult *= fabs(p_trans/mc_trans);  /* attenuate beam by portion which is scattered (and left along) */
    }

    if (flag) { /* propagate to exit of sample and finish */
      intersect = 0;
      p *= p_mult; /* apply absorption correction */
      PROP_DT(dt0+dt2);
      break; /* exit main multi scatt while loop */
    }
  } /* end if intersect the neutron hits the sample */
  else break;

  if (intersect) { /* scattering event */
    double kf=0, kf1, kf2;
    /* mean scattering probability and absorption fraction */
    VarSqw.mean_scatt += (1-exp(-VarSqw.my_s*d_path))*p;
    VarSqw.mean_abs   += (1-ws)*p;
    VarSqw.psum_scatt += p;

    /* Decaying exponential distribution of the path length before scattering */
    /* Select a point at which to scatter the neutron, taking
         secondary extinction into account. */
    if (my_t*d_path < 1e-6)
    /* For very weak scattering, use simple uniform sampling of scattering
       point to avoid rounding errors. */
      dt = rand0max(d_path); /* length */
    else
      dt = -log(1 - rand0max((1 - exp(-my_t*d_path)))) / my_t; /* length */
    dt /= v; /* Time from present position to scattering point */

    /* If t0 is in hole, propagate to next part of the hollow cylinder */
    if (dt1 > 0 && dt0 > 0 && dt > dt0) dt += dt1;

    /* Neutron propagation to the scattering point */
    PROP_DT(dt);

    /* choice between coherent/incoherent scattering */
    tmp_rand = rand01();
    /* local description at the scattering point (scat probability for atom) */
    tmp_rand *= (coh+inc);

    flag=0;
    if (VarSqw.s_inc>0 && tmp_rand < inc) {
      /* CASE 1: incoherent case */
      if (!VarSqw.Data_inc.intensity) {
        /* CASE 1a: no incoherent Sqw from file, use isotropic V-like */
        if (d_phi && order == 1) {
          randvec_target_rect_angular(&u1x, &u1y, &u1z, &solid_angle,
              vx, vy, vz, 2*PI, d_phi, ROT_A_CURRENT_COMP);
          p_mult *= solid_angle/4/PI; /* weighted by focused range to total range */
        } else
          randvec_target_circle(&u1x, &u1y, &u1z, NULL, vx, vy, vz, 0);

        vx = u1x; vy = u1y; vz = u1z;
        vf = v; kf = k;
        if (!VarSqw.type) VarSqw.type = 'v';
        SCATTER;
      } else {
        /* CASE 1b: incoherent Sqw from file */
        if (VarSqw.Data_inc.intensity) {
          Data_sqw = VarSqw.Data_inc;
          if (!VarSqw.type) VarSqw.type = 'i';
          flag = 1;
        }
      }
    } else if (VarSqw.s_coh>0 && tmp_rand > VarSqw.s_inc) {
      if (VarSqw.Data_coh.intensity) {
        /* CASE2: coherent case */
        Data_sqw = VarSqw.Data_coh;
        if (!VarSqw.type) VarSqw.type = 'c';
        flag = 1;
      }
    }

    if (flag) { /* true when S(q,w) table exists (Data_sqw) */

      double alpha=0, alpha0;
      /* give us a limited number of tries for scattering: choose W then Q */
      for (index_counter=VarSqw.maxloop; index_counter > 0 ; index_counter--) {

        /* MC choice: energy transfer w=Ei-Ef in the S(w) = SW */
        omega = 0;
        tmp_rand = rand01();
        /* energy index for rand > cumul SW */
        index_w  = Sqw_search_SW(Data_sqw, tmp_rand);
        VarSqw.rw = (double)index_w;
        if (index_w >= 0 && &(Data_sqw.SW[index_w]) != NULL) {
          if (Data_sqw.w_bins > 1) {
            double w1, w2;
            if (index_w > 0) { /* interpolate linearly energy */
              ratio_w = (tmp_rand                         - Data_sqw.SW[index_w-1].cumul_proba)
                       /(Data_sqw.SW[index_w].cumul_proba - Data_sqw.SW[index_w-1].cumul_proba);
              /* ratio_w=0 omega[index_w-1], ratio=1 omega[index] */
              w1 = Data_sqw.SW[index_w-1].omega; w2 = Data_sqw.SW[index_w].omega;
            } else { /* index_w = 0 interpolate to 0 energy */
              /* ratio_w=0 omega=0, ratio=1 omega[index] */
              w1 = Data_sqw.SW[index_w].omega; w2= Data_sqw.SW[index_w+1].omega;
              if (!w2 && index_w+1 < Data_sqw.w_bins)
                w2= Data_sqw.SW[index_w+1].omega;
              if (Data_sqw.w_bins && Data_sqw.SW[index_w].cumul_proba) {
                ratio_w = tmp_rand/Data_sqw.SW[index_w].cumul_proba;
              } else ratio_w=0;
            }
            if (ratio_w<0) ratio_w=0; else if (ratio_w>1) ratio_w=1;
            omega = (1-ratio_w)*w1 + ratio_w*w2;
          } else {
            ratio_w = 0;
            omega = Data_sqw.SW[index_w].omega;
          }
        } else {
          if (VarSqw.verbose_output >= 3 && VarSqw.neutron_removed<VarSqw.maxloop)
            printf("Isotropic_Sqw: %s: Warning: No suitable w transfer for index_w=%li.\n",
              NAME_CURRENT_COMP, index_w);
          continue; /* no W value: try again with an other energy transfer */
        }

        /* MC choice: momentum transfer Q in P(Q|w) */
        tmp_rand = rand01();

        /* momentum index for rand > cumul SQ|W */
        index_q  = Sqw_search_Q_proba_per_w(Data_sqw, tmp_rand, index_w);
        VarSqw.rq = (double)index_q;

        if (index_q >= 0 && &(Data_sqw.SQW[index_w]) != NULL) {
          if (Data_sqw.q_bins > 1 && index_q > 0) {
            if (index_w > 0 && Data_sqw.w_bins > 1) {
              /* bilinear interpolation on - side: index_w > 0, index_q > 0 */
              ratio_q = (tmp_rand - Data_sqw.SQW[index_w][index_q-1].cumul_proba)
                       /(Data_sqw.SQW[index_w][index_q].cumul_proba
                       - Data_sqw.SQW[index_w][index_q-1].cumul_proba);
              q22 = Data_sqw.SQW[index_w]  [index_q].Q;
              q11 = Data_sqw.SQW[index_w-1][index_q-1].Q;
              q21 = Data_sqw.SQW[index_w]  [index_q-1].Q;
              q12 = Data_sqw.SQW[index_w-1][index_q].Q;
              if (ratio_q<0) ratio_q=0; else if (ratio_q>1) ratio_q=1;
              q = (1-ratio_w)*(1-ratio_q)*q11+ratio_w*(1-ratio_q)*q21
                + ratio_w*ratio_q*q22        +(1-ratio_w)*ratio_q*q12;
            } else { /* bilinear interpolation on + side: index_w=0, index_q > 0 */
              ratio_q = (tmp_rand - Data_sqw.SQW[index_w][index_q-1].cumul_proba)
                       /(Data_sqw.SQW[index_w][index_q].cumul_proba
                       - Data_sqw.SQW[index_w][index_q-1].cumul_proba);
              q11 = Data_sqw.SQW[index_w]  [index_q-1].Q;
              q12 = Data_sqw.SQW[index_w]  [index_q].Q;
              if (ratio_q<0) ratio_q=0; else if (ratio_q>1) ratio_q=1;
              if (index_w < Data_sqw.w_bins-1 && Data_sqw.w_bins > 1) {
                q22 = Data_sqw.SQW[index_w+1][index_q].Q;
                q21 = Data_sqw.SQW[index_w+1][index_q-1].Q;
                q = (1-ratio_w)*(1-ratio_q)*q11+ratio_w*(1-ratio_q)*q21
                  + ratio_w*ratio_q*q22        +(1-ratio_w)*ratio_q*q12;
              } else {
                q    = (1-ratio_q)*q11  + ratio_q*q12;
              }
            }
          } else {
            q    = Data_sqw.SQW[index_w][index_q].Q;
          }
        } else {
          if (VarSqw.verbose_output >= 3 && VarSqw.neutron_removed<VarSqw.maxloop)
            printf("Isotropic_Sqw: %s: Warning: No suitable q transfer for w=%g.\n",
              NAME_CURRENT_COMP, omega);
          VarSqw.neutron_removed++;
          continue; /* no Q value for this w choice */
        }

        /* Search for length of final wave vector kf */
        /* kf is such that : hbar*w = hbar*hbar/2/m*(k*k - kf*kf) */
        /* acceptable values for kf are kf1 and kf2 */
        if (!solve_2nd_order(&kf1, &kf2, 1, 0, -k*k + VarSqw.sqSE2K*omega)) {
          if (VarSqw.verbose_output >= 3 && VarSqw.neutron_removed<VarSqw.maxloop)
            printf("Isotropic_Sqw: %s: Warning: imaginary root for w=%g q=%g Ei=%g (triangle can not close)\n",
            NAME_CURRENT_COMP, omega, q, Ei);
          VarSqw.neutron_removed++;
          continue; /* all roots are imaginary */
        }

        /* kf1 and kf2 are opposite */
        kf = fabs(kf1);
        vf = K2V*kf;

        /* Search of the direction of kf such that : q = ki - kf */
        /* cos theta = (ki2+kf2-q2)/(2ki kf) */

        costheta= (k*k+kf*kf-q*q)/(2*kf*k); /* this is cos(theta) */

        if (-1 < costheta && costheta < 1) {
          break; /* satisfies q momentum conservation */
        }
/*      else continue; */

        /* exit for loop on success */
      } /* end for index_counter */

      if (!index_counter) { /* for loop ended: failure for scattering */
        intersect=0; /* Could not scatter: finish multiple scattering loop */
        if (VarSqw.verbose_output >= 2 && VarSqw.neutron_removed<VarSqw.maxloop)
          printf("Isotropic_Sqw: %s: Warning: No scattering [q,w] conditions\n"
               "               last try (%i): type=%c w=%g q=%g cos(theta)=%g k=%g\n",
          NAME_CURRENT_COMP, VarSqw.maxloop, (VarSqw.type ? VarSqw.type : '-'), omega, q, costheta, k);
        VarSqw.neutron_removed++;
        if (order && SCATTERED != order) ABSORB;
        break;       /* finish multiple scattering loop */
      }

      /* scattering angle from ki to DS cone */
      theta = acos(costheta);

      /* Choose point on Debye-Scherrer cone */
      if (order == 1 && d_phi)
      { /* relate height of detector to the height on DS cone */
        double cone_focus;
        cone_focus = sin(d_phi/2)/sin(theta);
        /* If full Debye-Scherrer cone is within d_phi, don't focus */
        if (cone_focus < -1 || cone_focus > 1) d_phi = 0;
        /* Otherwise, determine alpha to rotate from scattering plane
            into d_phi focusing area*/
        else alpha = 2*asin(cone_focus);
        if (d_phi) p_mult *= alpha/PI;
      }
      if (d_phi) {
        /* Focusing */
        alpha = fabs(alpha);
        /* Trick to get scattering for pos/neg theta's */
        alpha0= 2*rand01()*alpha;
        if (alpha0 > alpha) {
          alpha0=PI+(alpha0-1.5*alpha);
        } else {
          alpha0=alpha0-0.5*alpha;
        }
      }
      else
        alpha0 = PI*randpm1();

      /* now find a nearly vertical rotation axis (u1) :
	       * Either
	       *  (v along Z) x (X axis) -> nearly Y axis
	       * Or
	       *  (v along X) x (Z axis) -> nearly Y axis
	       */
	    if (fabs(scalar_prod(1,0,0,vx/v,vy/v,vz/v)) < fabs(scalar_prod(0,0,1,vx/v,vy/v,vz/v))) {
        u1x = 1; u1y = u1z = 0;
	    } else {
        u1x = u1y = 0; u1z = 1;
	    }
	    vec_prod(u2x,u2y,u2z, vx,vy,vz, u1x,u1y,u1z);

      /* handle case where v and aim are parallel */
      if (!u2x && !u2y && !u2z) { u2x=u2z=0; u2y=1; }

      /* u1 = rotate 'v' by theta around u2: DS scattering angle, nearly in horz plane */
      rotate(u1x,u1y,u1z, vx,vy,vz, theta, u2x,u2y,u2z);

      /* u0 = rotate u1 by alpha0 around 'v' (Debye-Scherrer cone) */
      rotate(u0x,u0y,u0z, u1x,u1y,u1z, alpha0, vx, vy, vz);
      NORM(u0x,u0y,u0z);
      vx = u0x*vf;
      vy = u0y*vf;
      vz = u0z*vf;

      SCATTER;

      v = vf; k = kf; /* for next iteration */

    } /* end if (flag) */

    VarSqw.neutron_exit++;
    p *= p_mult;
    if (p_mult > 1) VarSqw.neutron_pmult++;

    /* test for a given multiple order */
    if (order && SCATTERED >= order) {
      intersect=0; /* reached required number of SCATTERing */
      break;       /* finish multiple scattering loop */
    }

  } /* end if (intersect) scattering event  */

} while (intersect); /* end do (intersect) (multiple scattering loop) */

/* Store Final neutron state */
VarSqw.kf_x = V2K*vx;
VarSqw.kf_y = V2K*vy;
VarSqw.kf_z = V2K*vz;
VarSqw.tf   = t;
VarSqw.vf   = v;
VarSqw.kf   = k;
VarSqw.theta= theta;

if (SCATTERED) {



  if (SCATTERED == 1) {
    if (VarSqw.type == 'c') VarSqw.single_coh += p;
    else                    VarSqw.single_inc += p;
    VarSqw.dq = sqrt((VarSqw.kf_x-VarSqw.ki_x)*(VarSqw.kf_x-VarSqw.ki_x)
                  +(VarSqw.kf_y-VarSqw.ki_y)*(VarSqw.kf_y-VarSqw.ki_y)
                  +(VarSqw.kf_z-VarSqw.ki_z)*(VarSqw.kf_z-VarSqw.ki_z));
    VarSqw.dw = VS2E*(VarSqw.vf*VarSqw.vf - VarSqw.vi*VarSqw.vi);
  } else VarSqw.multi += p;

} else VarSqw.dq=VarSqw.dw=0;

/* end TRACE */
}
#line 24162 "Samples_Incoherent.c"
}   /* End of SqwV=Isotropic_Sqw() SETTING parameter declarations. */
#undef offdata
#undef columns
#undef VarSqw
#undef powder_format
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompSqwV:
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

  /* TRACE Component Dirbeam [13] */
  mccoordschange(mcposrDirbeam, mcrotrDirbeam,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component Dirbeam (without coords transformations) */
  mcJumpTrace_Dirbeam:
  SIG_MESSAGE("Dirbeam (Trace)");
  mcDEBUG_COMP("Dirbeam")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompDirbeam
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
#define mccompcurname  Dirbeam
#define mccompcurtype  PSD_monitor
#define mccompcurindex 13
#define PSD_N mccDirbeam_PSD_N
#define PSD_p mccDirbeam_PSD_p
#define PSD_p2 mccDirbeam_PSD_p2
{   /* Declarations of Dirbeam=PSD_monitor() SETTING parameters. */
int nx = mccDirbeam_nx;
int ny = mccDirbeam_ny;
char* filename = mccDirbeam_filename;
MCNUM xmin = mccDirbeam_xmin;
MCNUM xmax = mccDirbeam_xmax;
MCNUM ymin = mccDirbeam_ymin;
MCNUM ymax = mccDirbeam_ymax;
MCNUM xwidth = mccDirbeam_xwidth;
MCNUM yheight = mccDirbeam_yheight;
MCNUM restore_neutron = mccDirbeam_restore_neutron;
/* 'Dirbeam=PSD_monitor()' component instance has conditional execution */
if (( mcipDB ))

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
#line 24302 "Samples_Incoherent.c"
}   /* End of Dirbeam=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompDirbeam:
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

  /* TRACE Component Stop [14] */
  mccoordschange(mcposrStop, mcrotrStop,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component Stop (without coords transformations) */
  mcJumpTrace_Stop:
  SIG_MESSAGE("Stop (Trace)");
  mcDEBUG_COMP("Stop")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompStop
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
#define mccompcurname  Stop
#define mccompcurtype  Beamstop
#define mccompcurindex 14
{   /* Declarations of Stop=Beamstop() SETTING parameters. */
MCNUM xmin = mccStop_xmin;
MCNUM xmax = mccStop_xmax;
MCNUM ymin = mccStop_ymin;
MCNUM ymax = mccStop_ymax;
MCNUM xwidth = mccStop_xwidth;
MCNUM yheight = mccStop_yheight;
MCNUM radius = mccStop_radius;
/* 'Stop=Beamstop()' component instance has conditional execution */
if (( mcipSTOP == 1 ))

#line 59 "/usr/share/mcstas/2.6rc1/optics/Beamstop.comp"
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
#line 24432 "Samples_Incoherent.c"
}   /* End of Stop=Beamstop() SETTING parameter declarations. */
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompStop:
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

  /* TRACE Component Tof_Sphere_4pi [15] */
  mccoordschange(mcposrTof_Sphere_4pi, mcrotrTof_Sphere_4pi,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component Tof_Sphere_4pi (without coords transformations) */
  mcJumpTrace_Tof_Sphere_4pi:
  SIG_MESSAGE("Tof_Sphere_4pi (Trace)");
  mcDEBUG_COMP("Tof_Sphere_4pi")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompTof_Sphere_4pi
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
#define mccompcurname  Tof_Sphere_4pi
#define mccompcurtype  Monitor_nD
#define mccompcurindex 15
#define user1 mccTof_Sphere_4pi_user1
#define user2 mccTof_Sphere_4pi_user2
#define user3 mccTof_Sphere_4pi_user3
#define DEFS mccTof_Sphere_4pi_DEFS
#define Vars mccTof_Sphere_4pi_Vars
#define detector mccTof_Sphere_4pi_detector
#define offdata mccTof_Sphere_4pi_offdata
{   /* Declarations of Tof_Sphere_4pi=Monitor_nD() SETTING parameters. */
MCNUM xwidth = mccTof_Sphere_4pi_xwidth;
MCNUM yheight = mccTof_Sphere_4pi_yheight;
MCNUM zdepth = mccTof_Sphere_4pi_zdepth;
MCNUM xmin = mccTof_Sphere_4pi_xmin;
MCNUM xmax = mccTof_Sphere_4pi_xmax;
MCNUM ymin = mccTof_Sphere_4pi_ymin;
MCNUM ymax = mccTof_Sphere_4pi_ymax;
MCNUM zmin = mccTof_Sphere_4pi_zmin;
MCNUM zmax = mccTof_Sphere_4pi_zmax;
MCNUM bins = mccTof_Sphere_4pi_bins;
MCNUM min = mccTof_Sphere_4pi_min;
MCNUM max = mccTof_Sphere_4pi_max;
MCNUM restore_neutron = mccTof_Sphere_4pi_restore_neutron;
MCNUM radius = mccTof_Sphere_4pi_radius;
char* options = mccTof_Sphere_4pi_options;
char* filename = mccTof_Sphere_4pi_filename;
char* geometry = mccTof_Sphere_4pi_geometry;
char* username1 = mccTof_Sphere_4pi_username1;
char* username2 = mccTof_Sphere_4pi_username2;
char* username3 = mccTof_Sphere_4pi_username3;
int nowritefile = mccTof_Sphere_4pi_nowritefile;
#line 313 "/usr/share/mcstas/2.6rc1/monitors/Monitor_nD.comp"
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
#line 24735 "Samples_Incoherent.c"
}   /* End of Tof_Sphere_4pi=Monitor_nD() SETTING parameter declarations. */
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
  mcabsorbCompTof_Sphere_4pi:
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

  /* TRACE Component PSD_Sphere_4pi [16] */
  mccoordschange(mcposrPSD_Sphere_4pi, mcrotrPSD_Sphere_4pi,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component PSD_Sphere_4pi (without coords transformations) */
  mcJumpTrace_PSD_Sphere_4pi:
  SIG_MESSAGE("PSD_Sphere_4pi (Trace)");
  mcDEBUG_COMP("PSD_Sphere_4pi")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompPSD_Sphere_4pi
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
#define mccompcurname  PSD_Sphere_4pi
#define mccompcurtype  PSD_monitor_4PI
#define mccompcurindex 16
#define nx mccPSD_Sphere_4pi_nx
#define ny mccPSD_Sphere_4pi_ny
#define PSD_N mccPSD_Sphere_4pi_PSD_N
#define PSD_p mccPSD_Sphere_4pi_PSD_p
#define PSD_p2 mccPSD_Sphere_4pi_PSD_p2
{   /* Declarations of PSD_Sphere_4pi=PSD_monitor_4PI() SETTING parameters. */
char* filename = mccPSD_Sphere_4pi_filename;
MCNUM radius = mccPSD_Sphere_4pi_radius;
MCNUM restore_neutron = mccPSD_Sphere_4pi_restore_neutron;
int nowritefile = mccPSD_Sphere_4pi_nowritefile;
#line 74 "/usr/share/mcstas/2.6rc1/monitors/PSD_monitor_4PI.comp"
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
#line 24889 "Samples_Incoherent.c"
}   /* End of PSD_Sphere_4pi=PSD_monitor_4PI() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef ny
#undef nx
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompPSD_Sphere_4pi:
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

  /* TRACE Component Coord3 [17] */
  mccoordschange(mcposrCoord3, mcrotrCoord3,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component Coord3 (without coords transformations) */
  mcJumpTrace_Coord3:
  SIG_MESSAGE("Coord3 (Trace)");
  mcDEBUG_COMP("Coord3")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompCoord3
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
#define mccompcurname  Coord3
#define mccompcurtype  Arm
#define mccompcurindex 17
/* 'Coord3=Arm()' component instance extend code */
    SIG_MESSAGE("Coord3 (Trace:Extend)");
#line 177 "Samples_Incoherent.instr"
  p=1;
#line 25003 "Samples_Incoherent.c"
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompCoord3:
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

  /* TRACE Component Sphere_4pi [18] */
  mccoordschange(mcposrSphere_4pi, mcrotrSphere_4pi,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component Sphere_4pi (without coords transformations) */
  mcJumpTrace_Sphere_4pi:
  SIG_MESSAGE("Sphere_4pi (Trace)");
  mcDEBUG_COMP("Sphere_4pi")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompSphere_4pi
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
#define mccompcurname  Sphere_4pi
#define mccompcurtype  PSD_monitor_4PI
#define mccompcurindex 18
#define nx mccSphere_4pi_nx
#define ny mccSphere_4pi_ny
#define PSD_N mccSphere_4pi_PSD_N
#define PSD_p mccSphere_4pi_PSD_p
#define PSD_p2 mccSphere_4pi_PSD_p2
{   /* Declarations of Sphere_4pi=PSD_monitor_4PI() SETTING parameters. */
char* filename = mccSphere_4pi_filename;
MCNUM radius = mccSphere_4pi_radius;
MCNUM restore_neutron = mccSphere_4pi_restore_neutron;
int nowritefile = mccSphere_4pi_nowritefile;
#line 74 "/usr/share/mcstas/2.6rc1/monitors/PSD_monitor_4PI.comp"
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
#line 25149 "Samples_Incoherent.c"
}   /* End of Sphere_4pi=PSD_monitor_4PI() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef ny
#undef nx
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompSphere_4pi:
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

  /* User SAVE code for component 'Coord'. */
  SIG_MESSAGE("Coord (Save)");
#define mccompcurname  Coord
#define mccompcurtype  Progress_bar
#define mccompcurindex 2
#define IntermediateCnts mccCoord_IntermediateCnts
#define StartTime mccCoord_StartTime
#define EndTime mccCoord_EndTime
#define CurrentTime mccCoord_CurrentTime
{   /* Declarations of Coord=Progress_bar() SETTING parameters. */
char* profile = mccCoord_profile;
MCNUM percent = mccCoord_percent;
MCNUM flag_save = mccCoord_flag_save;
MCNUM minutes = mccCoord_minutes;
#line 115 "/usr/share/mcstas/2.6rc1/misc/Progress_bar.comp"
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
#line 25263 "Samples_Incoherent.c"
}   /* End of Coord=Progress_bar() SETTING parameter declarations. */
#undef CurrentTime
#undef EndTime
#undef StartTime
#undef IntermediateCnts
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'Dirbeam'. */
  SIG_MESSAGE("Dirbeam (Save)");
#define mccompcurname  Dirbeam
#define mccompcurtype  PSD_monitor
#define mccompcurindex 13
#define PSD_N mccDirbeam_PSD_N
#define PSD_p mccDirbeam_PSD_p
#define PSD_p2 mccDirbeam_PSD_p2
{   /* Declarations of Dirbeam=PSD_monitor() SETTING parameters. */
int nx = mccDirbeam_nx;
int ny = mccDirbeam_ny;
char* filename = mccDirbeam_filename;
MCNUM xmin = mccDirbeam_xmin;
MCNUM xmax = mccDirbeam_xmax;
MCNUM ymin = mccDirbeam_ymin;
MCNUM ymax = mccDirbeam_ymax;
MCNUM xwidth = mccDirbeam_xwidth;
MCNUM yheight = mccDirbeam_yheight;
MCNUM restore_neutron = mccDirbeam_restore_neutron;
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
#line 25303 "Samples_Incoherent.c"
}   /* End of Dirbeam=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'Tof_Sphere_4pi'. */
  SIG_MESSAGE("Tof_Sphere_4pi (Save)");
#define mccompcurname  Tof_Sphere_4pi
#define mccompcurtype  Monitor_nD
#define mccompcurindex 15
#define user1 mccTof_Sphere_4pi_user1
#define user2 mccTof_Sphere_4pi_user2
#define user3 mccTof_Sphere_4pi_user3
#define DEFS mccTof_Sphere_4pi_DEFS
#define Vars mccTof_Sphere_4pi_Vars
#define detector mccTof_Sphere_4pi_detector
#define offdata mccTof_Sphere_4pi_offdata
{   /* Declarations of Tof_Sphere_4pi=Monitor_nD() SETTING parameters. */
MCNUM xwidth = mccTof_Sphere_4pi_xwidth;
MCNUM yheight = mccTof_Sphere_4pi_yheight;
MCNUM zdepth = mccTof_Sphere_4pi_zdepth;
MCNUM xmin = mccTof_Sphere_4pi_xmin;
MCNUM xmax = mccTof_Sphere_4pi_xmax;
MCNUM ymin = mccTof_Sphere_4pi_ymin;
MCNUM ymax = mccTof_Sphere_4pi_ymax;
MCNUM zmin = mccTof_Sphere_4pi_zmin;
MCNUM zmax = mccTof_Sphere_4pi_zmax;
MCNUM bins = mccTof_Sphere_4pi_bins;
MCNUM min = mccTof_Sphere_4pi_min;
MCNUM max = mccTof_Sphere_4pi_max;
MCNUM restore_neutron = mccTof_Sphere_4pi_restore_neutron;
MCNUM radius = mccTof_Sphere_4pi_radius;
char* options = mccTof_Sphere_4pi_options;
char* filename = mccTof_Sphere_4pi_filename;
char* geometry = mccTof_Sphere_4pi_geometry;
char* username1 = mccTof_Sphere_4pi_username1;
char* username2 = mccTof_Sphere_4pi_username2;
char* username3 = mccTof_Sphere_4pi_username3;
int nowritefile = mccTof_Sphere_4pi_nowritefile;
#line 483 "/usr/share/mcstas/2.6rc1/monitors/Monitor_nD.comp"
{
  /* save results, but do not free pointers */
  if (!nowritefile) {
    detector = Monitor_nD_Save(&DEFS, &Vars);
  }
}
#line 25353 "Samples_Incoherent.c"
}   /* End of Tof_Sphere_4pi=Monitor_nD() SETTING parameter declarations. */
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

  /* User SAVE code for component 'PSD_Sphere_4pi'. */
  SIG_MESSAGE("PSD_Sphere_4pi (Save)");
#define mccompcurname  PSD_Sphere_4pi
#define mccompcurtype  PSD_monitor_4PI
#define mccompcurindex 16
#define nx mccPSD_Sphere_4pi_nx
#define ny mccPSD_Sphere_4pi_ny
#define PSD_N mccPSD_Sphere_4pi_PSD_N
#define PSD_p mccPSD_Sphere_4pi_PSD_p
#define PSD_p2 mccPSD_Sphere_4pi_PSD_p2
{   /* Declarations of PSD_Sphere_4pi=PSD_monitor_4PI() SETTING parameters. */
char* filename = mccPSD_Sphere_4pi_filename;
MCNUM radius = mccPSD_Sphere_4pi_radius;
MCNUM restore_neutron = mccPSD_Sphere_4pi_restore_neutron;
int nowritefile = mccPSD_Sphere_4pi_nowritefile;
#line 107 "/usr/share/mcstas/2.6rc1/monitors/PSD_monitor_4PI.comp"
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
#line 25394 "Samples_Incoherent.c"
}   /* End of PSD_Sphere_4pi=PSD_monitor_4PI() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef ny
#undef nx
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'Sphere_4pi'. */
  SIG_MESSAGE("Sphere_4pi (Save)");
#define mccompcurname  Sphere_4pi
#define mccompcurtype  PSD_monitor_4PI
#define mccompcurindex 18
#define nx mccSphere_4pi_nx
#define ny mccSphere_4pi_ny
#define PSD_N mccSphere_4pi_PSD_N
#define PSD_p mccSphere_4pi_PSD_p
#define PSD_p2 mccSphere_4pi_PSD_p2
{   /* Declarations of Sphere_4pi=PSD_monitor_4PI() SETTING parameters. */
char* filename = mccSphere_4pi_filename;
MCNUM radius = mccSphere_4pi_radius;
MCNUM restore_neutron = mccSphere_4pi_restore_neutron;
int nowritefile = mccSphere_4pi_nowritefile;
#line 107 "/usr/share/mcstas/2.6rc1/monitors/PSD_monitor_4PI.comp"
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
#line 25433 "Samples_Incoherent.c"
}   /* End of Sphere_4pi=PSD_monitor_4PI() SETTING parameter declarations. */
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

    if (!mcNCounter[1]) fprintf(stderr, "Warning: No neutron could reach Component[1] Source\n");
    if (mcAbsorbProp[1]) fprintf(stderr, "Warning: %g events were removed in Component[1] Source=ISIS_moderator()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[1]);
  /* User FINALLY code for component 'Coord'. */
  SIG_MESSAGE("Coord (Finally)");
#define mccompcurname  Coord
#define mccompcurtype  Progress_bar
#define mccompcurindex 2
#define IntermediateCnts mccCoord_IntermediateCnts
#define StartTime mccCoord_StartTime
#define EndTime mccCoord_EndTime
#define CurrentTime mccCoord_CurrentTime
{   /* Declarations of Coord=Progress_bar() SETTING parameters. */
char* profile = mccCoord_profile;
MCNUM percent = mccCoord_percent;
MCNUM flag_save = mccCoord_flag_save;
MCNUM minutes = mccCoord_minutes;
#line 133 "/usr/share/mcstas/2.6rc1/misc/Progress_bar.comp"
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
#line 25480 "Samples_Incoherent.c"
}   /* End of Coord=Progress_bar() SETTING parameter declarations. */
#undef CurrentTime
#undef EndTime
#undef StartTime
#undef IntermediateCnts
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[2]) fprintf(stderr, "Warning: No neutron could reach Component[2] Coord\n");
    if (mcAbsorbProp[2]) fprintf(stderr, "Warning: %g events were removed in Component[2] Coord=Progress_bar()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[2]);
    if (!mcNCounter[3]) fprintf(stderr, "Warning: No neutron could reach Component[3] Vsample\n");
    if (mcAbsorbProp[3]) fprintf(stderr, "Warning: %g events were removed in Component[3] Vsample=V_sample()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[3]);
  /* User FINALLY code for component 'PowN'. */
  SIG_MESSAGE("PowN (Finally)");
#define mccompcurname  PowN
#define mccompcurtype  PowderN
#define mccompcurindex 4
#define format mccPowN_format
#define line_info mccPowN_line_info
#define columns mccPowN_columns
#define offdata mccPowN_offdata
{   /* Declarations of PowN=PowderN() SETTING parameters. */
char* reflections = mccPowN_reflections;
char* geometry = mccPowN_geometry;
MCNUM radius = mccPowN_radius;
MCNUM yheight = mccPowN_yheight;
MCNUM xwidth = mccPowN_xwidth;
MCNUM zdepth = mccPowN_zdepth;
MCNUM thickness = mccPowN_thickness;
MCNUM pack = mccPowN_pack;
MCNUM Vc = mccPowN_Vc;
MCNUM sigma_abs = mccPowN_sigma_abs;
MCNUM sigma_inc = mccPowN_sigma_inc;
MCNUM delta_d_d = mccPowN_delta_d_d;
MCNUM p_inc = mccPowN_p_inc;
MCNUM p_transmit = mccPowN_p_transmit;
MCNUM DW = mccPowN_DW;
MCNUM nb_atoms = mccPowN_nb_atoms;
MCNUM d_omega = mccPowN_d_omega;
MCNUM d_phi = mccPowN_d_phi;
MCNUM tth_sign = mccPowN_tth_sign;
MCNUM p_interact = mccPowN_p_interact;
MCNUM concentric = mccPowN_concentric;
MCNUM density = mccPowN_density;
MCNUM weight = mccPowN_weight;
MCNUM barns = mccPowN_barns;
MCNUM Strain = mccPowN_Strain;
MCNUM focus_flip = mccPowN_focus_flip;
int target_index = mccPowN_target_index;
#line 1085 "/usr/share/mcstas/2.6rc1/samples/PowderN.comp"
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
#line 25556 "Samples_Incoherent.c"
}   /* End of PowN=PowderN() SETTING parameter declarations. */
#undef offdata
#undef columns
#undef line_info
#undef format
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[4]) fprintf(stderr, "Warning: No neutron could reach Component[4] PowN\n");
    if (mcAbsorbProp[4]) fprintf(stderr, "Warning: %g events were removed in Component[4] PowN=PowderN()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[4]);
  /* User FINALLY code for component 'SX'. */
  SIG_MESSAGE("SX (Finally)");
#define mccompcurname  SX
#define mccompcurtype  Single_crystal
#define mccompcurindex 5
#define mosaic_AB mccSX_mosaic_AB
#define hkl_info mccSX_hkl_info
#define offdata mccSX_offdata
{   /* Declarations of SX=Single_crystal() SETTING parameters. */
char* reflections = mccSX_reflections;
char* geometry = mccSX_geometry;
MCNUM xwidth = mccSX_xwidth;
MCNUM yheight = mccSX_yheight;
MCNUM zdepth = mccSX_zdepth;
MCNUM radius = mccSX_radius;
MCNUM delta_d_d = mccSX_delta_d_d;
MCNUM mosaic = mccSX_mosaic;
MCNUM mosaic_a = mccSX_mosaic_a;
MCNUM mosaic_b = mccSX_mosaic_b;
MCNUM mosaic_c = mccSX_mosaic_c;
MCNUM recip_cell = mccSX_recip_cell;
MCNUM barns = mccSX_barns;
MCNUM ax = mccSX_ax;
MCNUM ay = mccSX_ay;
MCNUM az = mccSX_az;
MCNUM bx = mccSX_bx;
MCNUM by = mccSX_by;
MCNUM bz = mccSX_bz;
MCNUM cx = mccSX_cx;
MCNUM cy = mccSX_cy;
MCNUM cz = mccSX_cz;
MCNUM p_transmit = mccSX_p_transmit;
MCNUM sigma_abs = mccSX_sigma_abs;
MCNUM sigma_inc = mccSX_sigma_inc;
MCNUM aa = mccSX_aa;
MCNUM bb = mccSX_bb;
MCNUM cc = mccSX_cc;
MCNUM order = mccSX_order;
MCNUM RX = mccSX_RX;
MCNUM RY = mccSX_RY;
MCNUM powder = mccSX_powder;
MCNUM PG = mccSX_PG;
MCNUM deltak = mccSX_deltak;
#line 1410 "/usr/share/mcstas/2.6rc1/samples/Single_crystal.comp"
{
  MPI_MASTER(
  if (hkl_info.flag_warning)
    fprintf(stderr, "Single_crystal: %s: Error message was repeated %i times with absorbed neutrons.\n",
      NAME_CURRENT_COMP, hkl_info.flag_warning);

  /* in case this instance is used in a SPLIT, we can recommend the
     optimal iteration value */
  if (hkl_info.nb_refl_count) {
    double split_iterations = (double)hkl_info.nb_reuses/hkl_info.nb_refl_count + 1;
    double split_optimal    = (double)hkl_info.nb_refl/hkl_info.nb_refl_count;
    if (split_optimal > split_iterations + 2)
      printf("Single_crystal: %s: Info: you may highly improve the computation efficiency by using\n"
        "    SPLIT %i COMPONENT %s=Single_crystal(order=1, ...)\n"
        "  in the instrument description %s.\n",
        NAME_CURRENT_COMP, (int)split_optimal, NAME_CURRENT_COMP, mcinstrument_source);
  }
  );
}
#line 25631 "Samples_Incoherent.c"
}   /* End of SX=Single_crystal() SETTING parameter declarations. */
#undef offdata
#undef hkl_info
#undef mosaic_AB
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[5]) fprintf(stderr, "Warning: No neutron could reach Component[5] SX\n");
    if (mcAbsorbProp[5]) fprintf(stderr, "Warning: %g events were removed in Component[5] SX=Single_crystal()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[5]);
  /* User FINALLY code for component 'Sqw'. */
  SIG_MESSAGE("Sqw (Finally)");
#define mccompcurname  Sqw
#define mccompcurtype  Isotropic_Sqw
#define mccompcurindex 6
#define powder_format mccSqw_powder_format
#define VarSqw mccSqw_VarSqw
#define columns mccSqw_columns
#define offdata mccSqw_offdata
{   /* Declarations of Sqw=Isotropic_Sqw() SETTING parameters. */
char* Sqw_coh = mccSqw_Sqw_coh;
char* Sqw_inc = mccSqw_Sqw_inc;
char* geometry = mccSqw_geometry;
MCNUM radius = mccSqw_radius;
MCNUM thickness = mccSqw_thickness;
MCNUM xwidth = mccSqw_xwidth;
MCNUM yheight = mccSqw_yheight;
MCNUM zdepth = mccSqw_zdepth;
MCNUM threshold = mccSqw_threshold;
int order = mccSqw_order;
MCNUM T = mccSqw_T;
MCNUM verbose = mccSqw_verbose;
MCNUM d_phi = mccSqw_d_phi;
int concentric = mccSqw_concentric;
MCNUM rho = mccSqw_rho;
MCNUM sigma_abs = mccSqw_sigma_abs;
MCNUM sigma_coh = mccSqw_sigma_coh;
MCNUM sigma_inc = mccSqw_sigma_inc;
MCNUM classical = mccSqw_classical;
MCNUM powder_Dd = mccSqw_powder_Dd;
MCNUM powder_DW = mccSqw_powder_DW;
MCNUM powder_Vc = mccSqw_powder_Vc;
MCNUM density = mccSqw_density;
MCNUM weight = mccSqw_weight;
MCNUM p_interact = mccSqw_p_interact;
MCNUM norm = mccSqw_norm;
MCNUM powder_barns = mccSqw_powder_barns;
char* quantum_correction = mccSqw_quantum_correction;
#line 2701 "/usr/share/mcstas/2.6rc1/samples/Isotropic_Sqw.comp"
{
  int  k;

  if (VarSqw.s_coh > 0 || VarSqw.s_inc > 0)
  for (k=0; k < 2; k++) {
    struct Sqw_Data_struct Data_sqw;

    Data_sqw =  (k == 0 ? VarSqw.Data_coh : VarSqw.Data_inc);
    /* Data_sqw->Sqw has already been freed at end of INIT */
    Table_Free(&(Data_sqw.iqSq));

    if (Data_sqw.SW)           free(Data_sqw.SW);
    if (Data_sqw.SQW)          free(Data_sqw.SQW);
    if (Data_sqw.SW_lookup)    free(Data_sqw.SW_lookup);
    if (Data_sqw.QW_lookup)    free(Data_sqw.QW_lookup);
  } /* end for */

#ifdef USE_MPI
  if (mpi_node_count > 1) {
    double tmp;
    tmp = (double)VarSqw.neutron_removed; mc_MPI_Sum(&tmp, 1); VarSqw.neutron_removed=(long)tmp;
    tmp = (double)VarSqw.neutron_exit;    mc_MPI_Sum(&tmp, 1); VarSqw.neutron_exit=(long)tmp;
    tmp = (double)VarSqw.neutron_pmult;   mc_MPI_Sum(&tmp, 1); VarSqw.neutron_pmult=(long)tmp;
    mc_MPI_Sum(&VarSqw.mean_scatt, 1);
    mc_MPI_Sum(&VarSqw.psum_scatt, 1);
    mc_MPI_Sum(&VarSqw.mean_abs, 1);
    mc_MPI_Sum(&VarSqw.single_coh, 1);
    mc_MPI_Sum(&VarSqw.single_inc, 1);
    mc_MPI_Sum(&VarSqw.multi, 1);
  }
#endif
  MPI_MASTER(
  if (VarSqw.neutron_removed)
    printf("Isotropic_Sqw: %s: %li neutron events (out of %li) that should have\n"
           "               scattered were transmitted because scattering conditions\n"
           "WARNING        could not be satisfied after %i tries.\n",
          NAME_CURRENT_COMP, VarSqw.neutron_removed,
          VarSqw.neutron_exit+VarSqw.neutron_removed, VarSqw.maxloop);
  if (VarSqw.neutron_pmult)
    printf("Isotropic_Sqw: %s: %li neutron events (out of %li) reached\n"
           "WARNING        unrealistic weight. The S(q,w) norm might be too high.\n",
          NAME_CURRENT_COMP, VarSqw.neutron_pmult, VarSqw.neutron_exit);

  if (VarSqw.verbose_output >= 1 && VarSqw.psum_scatt > 0) {
    printf("Isotropic_Sqw: %s: Scattering fraction=%g of incoming intensity\n"
           "               Absorption fraction           =%g\n",
           NAME_CURRENT_COMP,
           VarSqw.mean_scatt/VarSqw.psum_scatt, VarSqw.mean_abs/VarSqw.psum_scatt);
    printf("               Single   scattering intensity =%g (coh=%g inc=%g)\n"
           "               Multiple scattering intensity =%g\n",
           VarSqw.single_coh+VarSqw.single_inc, VarSqw.single_coh, VarSqw.single_inc, VarSqw.multi);
    );
  }

/* end FINALLY */
}
#line 25737 "Samples_Incoherent.c"
}   /* End of Sqw=Isotropic_Sqw() SETTING parameter declarations. */
#undef offdata
#undef columns
#undef VarSqw
#undef powder_format
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[6]) fprintf(stderr, "Warning: No neutron could reach Component[6] Sqw\n");
    if (mcAbsorbProp[6]) fprintf(stderr, "Warning: %g events were removed in Component[6] Sqw=Isotropic_Sqw()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[6]);
    if (!mcNCounter[7]) fprintf(stderr, "Warning: No neutron could reach Component[7] IncS\n");
    if (mcAbsorbProp[7]) fprintf(stderr, "Warning: %g events were removed in Component[7] IncS=Incoherent()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[7]);
  /* User FINALLY code for component 'PowNg'. */
  SIG_MESSAGE("PowNg (Finally)");
#define mccompcurname  PowNg
#define mccompcurtype  PowderN
#define mccompcurindex 8
#define format mccPowNg_format
#define line_info mccPowNg_line_info
#define columns mccPowNg_columns
#define offdata mccPowNg_offdata
{   /* Declarations of PowNg=PowderN() SETTING parameters. */
char* reflections = mccPowNg_reflections;
char* geometry = mccPowNg_geometry;
MCNUM radius = mccPowNg_radius;
MCNUM yheight = mccPowNg_yheight;
MCNUM xwidth = mccPowNg_xwidth;
MCNUM zdepth = mccPowNg_zdepth;
MCNUM thickness = mccPowNg_thickness;
MCNUM pack = mccPowNg_pack;
MCNUM Vc = mccPowNg_Vc;
MCNUM sigma_abs = mccPowNg_sigma_abs;
MCNUM sigma_inc = mccPowNg_sigma_inc;
MCNUM delta_d_d = mccPowNg_delta_d_d;
MCNUM p_inc = mccPowNg_p_inc;
MCNUM p_transmit = mccPowNg_p_transmit;
MCNUM DW = mccPowNg_DW;
MCNUM nb_atoms = mccPowNg_nb_atoms;
MCNUM d_omega = mccPowNg_d_omega;
MCNUM d_phi = mccPowNg_d_phi;
MCNUM tth_sign = mccPowNg_tth_sign;
MCNUM p_interact = mccPowNg_p_interact;
MCNUM concentric = mccPowNg_concentric;
MCNUM density = mccPowNg_density;
MCNUM weight = mccPowNg_weight;
MCNUM barns = mccPowNg_barns;
MCNUM Strain = mccPowNg_Strain;
MCNUM focus_flip = mccPowNg_focus_flip;
int target_index = mccPowNg_target_index;
#line 1085 "/usr/share/mcstas/2.6rc1/samples/PowderN.comp"
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
#line 25813 "Samples_Incoherent.c"
}   /* End of PowNg=PowderN() SETTING parameter declarations. */
#undef offdata
#undef columns
#undef line_info
#undef format
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[8]) fprintf(stderr, "Warning: No neutron could reach Component[8] PowNg\n");
    if (mcAbsorbProp[8]) fprintf(stderr, "Warning: %g events were removed in Component[8] PowNg=PowderN()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[8]);
  /* User FINALLY code for component 'SXg'. */
  SIG_MESSAGE("SXg (Finally)");
#define mccompcurname  SXg
#define mccompcurtype  Single_crystal
#define mccompcurindex 9
#define mosaic_AB mccSXg_mosaic_AB
#define hkl_info mccSXg_hkl_info
#define offdata mccSXg_offdata
{   /* Declarations of SXg=Single_crystal() SETTING parameters. */
char* reflections = mccSXg_reflections;
char* geometry = mccSXg_geometry;
MCNUM xwidth = mccSXg_xwidth;
MCNUM yheight = mccSXg_yheight;
MCNUM zdepth = mccSXg_zdepth;
MCNUM radius = mccSXg_radius;
MCNUM delta_d_d = mccSXg_delta_d_d;
MCNUM mosaic = mccSXg_mosaic;
MCNUM mosaic_a = mccSXg_mosaic_a;
MCNUM mosaic_b = mccSXg_mosaic_b;
MCNUM mosaic_c = mccSXg_mosaic_c;
MCNUM recip_cell = mccSXg_recip_cell;
MCNUM barns = mccSXg_barns;
MCNUM ax = mccSXg_ax;
MCNUM ay = mccSXg_ay;
MCNUM az = mccSXg_az;
MCNUM bx = mccSXg_bx;
MCNUM by = mccSXg_by;
MCNUM bz = mccSXg_bz;
MCNUM cx = mccSXg_cx;
MCNUM cy = mccSXg_cy;
MCNUM cz = mccSXg_cz;
MCNUM p_transmit = mccSXg_p_transmit;
MCNUM sigma_abs = mccSXg_sigma_abs;
MCNUM sigma_inc = mccSXg_sigma_inc;
MCNUM aa = mccSXg_aa;
MCNUM bb = mccSXg_bb;
MCNUM cc = mccSXg_cc;
MCNUM order = mccSXg_order;
MCNUM RX = mccSXg_RX;
MCNUM RY = mccSXg_RY;
MCNUM powder = mccSXg_powder;
MCNUM PG = mccSXg_PG;
MCNUM deltak = mccSXg_deltak;
#line 1410 "/usr/share/mcstas/2.6rc1/samples/Single_crystal.comp"
{
  MPI_MASTER(
  if (hkl_info.flag_warning)
    fprintf(stderr, "Single_crystal: %s: Error message was repeated %i times with absorbed neutrons.\n",
      NAME_CURRENT_COMP, hkl_info.flag_warning);

  /* in case this instance is used in a SPLIT, we can recommend the
     optimal iteration value */
  if (hkl_info.nb_refl_count) {
    double split_iterations = (double)hkl_info.nb_reuses/hkl_info.nb_refl_count + 1;
    double split_optimal    = (double)hkl_info.nb_refl/hkl_info.nb_refl_count;
    if (split_optimal > split_iterations + 2)
      printf("Single_crystal: %s: Info: you may highly improve the computation efficiency by using\n"
        "    SPLIT %i COMPONENT %s=Single_crystal(order=1, ...)\n"
        "  in the instrument description %s.\n",
        NAME_CURRENT_COMP, (int)split_optimal, NAME_CURRENT_COMP, mcinstrument_source);
  }
  );
}
#line 25888 "Samples_Incoherent.c"
}   /* End of SXg=Single_crystal() SETTING parameter declarations. */
#undef offdata
#undef hkl_info
#undef mosaic_AB
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[9]) fprintf(stderr, "Warning: No neutron could reach Component[9] SXg\n");
    if (mcAbsorbProp[9]) fprintf(stderr, "Warning: %g events were removed in Component[9] SXg=Single_crystal()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[9]);
  /* User FINALLY code for component 'Sqwg'. */
  SIG_MESSAGE("Sqwg (Finally)");
#define mccompcurname  Sqwg
#define mccompcurtype  Isotropic_Sqw
#define mccompcurindex 10
#define powder_format mccSqwg_powder_format
#define VarSqw mccSqwg_VarSqw
#define columns mccSqwg_columns
#define offdata mccSqwg_offdata
{   /* Declarations of Sqwg=Isotropic_Sqw() SETTING parameters. */
char* Sqw_coh = mccSqwg_Sqw_coh;
char* Sqw_inc = mccSqwg_Sqw_inc;
char* geometry = mccSqwg_geometry;
MCNUM radius = mccSqwg_radius;
MCNUM thickness = mccSqwg_thickness;
MCNUM xwidth = mccSqwg_xwidth;
MCNUM yheight = mccSqwg_yheight;
MCNUM zdepth = mccSqwg_zdepth;
MCNUM threshold = mccSqwg_threshold;
int order = mccSqwg_order;
MCNUM T = mccSqwg_T;
MCNUM verbose = mccSqwg_verbose;
MCNUM d_phi = mccSqwg_d_phi;
int concentric = mccSqwg_concentric;
MCNUM rho = mccSqwg_rho;
MCNUM sigma_abs = mccSqwg_sigma_abs;
MCNUM sigma_coh = mccSqwg_sigma_coh;
MCNUM sigma_inc = mccSqwg_sigma_inc;
MCNUM classical = mccSqwg_classical;
MCNUM powder_Dd = mccSqwg_powder_Dd;
MCNUM powder_DW = mccSqwg_powder_DW;
MCNUM powder_Vc = mccSqwg_powder_Vc;
MCNUM density = mccSqwg_density;
MCNUM weight = mccSqwg_weight;
MCNUM p_interact = mccSqwg_p_interact;
MCNUM norm = mccSqwg_norm;
MCNUM powder_barns = mccSqwg_powder_barns;
char* quantum_correction = mccSqwg_quantum_correction;
#line 2701 "/usr/share/mcstas/2.6rc1/samples/Isotropic_Sqw.comp"
{
  int  k;

  if (VarSqw.s_coh > 0 || VarSqw.s_inc > 0)
  for (k=0; k < 2; k++) {
    struct Sqw_Data_struct Data_sqw;

    Data_sqw =  (k == 0 ? VarSqw.Data_coh : VarSqw.Data_inc);
    /* Data_sqw->Sqw has already been freed at end of INIT */
    Table_Free(&(Data_sqw.iqSq));

    if (Data_sqw.SW)           free(Data_sqw.SW);
    if (Data_sqw.SQW)          free(Data_sqw.SQW);
    if (Data_sqw.SW_lookup)    free(Data_sqw.SW_lookup);
    if (Data_sqw.QW_lookup)    free(Data_sqw.QW_lookup);
  } /* end for */

#ifdef USE_MPI
  if (mpi_node_count > 1) {
    double tmp;
    tmp = (double)VarSqw.neutron_removed; mc_MPI_Sum(&tmp, 1); VarSqw.neutron_removed=(long)tmp;
    tmp = (double)VarSqw.neutron_exit;    mc_MPI_Sum(&tmp, 1); VarSqw.neutron_exit=(long)tmp;
    tmp = (double)VarSqw.neutron_pmult;   mc_MPI_Sum(&tmp, 1); VarSqw.neutron_pmult=(long)tmp;
    mc_MPI_Sum(&VarSqw.mean_scatt, 1);
    mc_MPI_Sum(&VarSqw.psum_scatt, 1);
    mc_MPI_Sum(&VarSqw.mean_abs, 1);
    mc_MPI_Sum(&VarSqw.single_coh, 1);
    mc_MPI_Sum(&VarSqw.single_inc, 1);
    mc_MPI_Sum(&VarSqw.multi, 1);
  }
#endif
  MPI_MASTER(
  if (VarSqw.neutron_removed)
    printf("Isotropic_Sqw: %s: %li neutron events (out of %li) that should have\n"
           "               scattered were transmitted because scattering conditions\n"
           "WARNING        could not be satisfied after %i tries.\n",
          NAME_CURRENT_COMP, VarSqw.neutron_removed,
          VarSqw.neutron_exit+VarSqw.neutron_removed, VarSqw.maxloop);
  if (VarSqw.neutron_pmult)
    printf("Isotropic_Sqw: %s: %li neutron events (out of %li) reached\n"
           "WARNING        unrealistic weight. The S(q,w) norm might be too high.\n",
          NAME_CURRENT_COMP, VarSqw.neutron_pmult, VarSqw.neutron_exit);

  if (VarSqw.verbose_output >= 1 && VarSqw.psum_scatt > 0) {
    printf("Isotropic_Sqw: %s: Scattering fraction=%g of incoming intensity\n"
           "               Absorption fraction           =%g\n",
           NAME_CURRENT_COMP,
           VarSqw.mean_scatt/VarSqw.psum_scatt, VarSqw.mean_abs/VarSqw.psum_scatt);
    printf("               Single   scattering intensity =%g (coh=%g inc=%g)\n"
           "               Multiple scattering intensity =%g\n",
           VarSqw.single_coh+VarSqw.single_inc, VarSqw.single_coh, VarSqw.single_inc, VarSqw.multi);
    );
  }

/* end FINALLY */
}
#line 25994 "Samples_Incoherent.c"
}   /* End of Sqwg=Isotropic_Sqw() SETTING parameter declarations. */
#undef offdata
#undef columns
#undef VarSqw
#undef powder_format
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[10]) fprintf(stderr, "Warning: No neutron could reach Component[10] Sqwg\n");
    if (mcAbsorbProp[10]) fprintf(stderr, "Warning: %g events were removed in Component[10] Sqwg=Isotropic_Sqw()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[10]);
    if (!mcNCounter[11]) fprintf(stderr, "Warning: No neutron could reach Component[11] IncSg\n");
    if (mcAbsorbProp[11]) fprintf(stderr, "Warning: %g events were removed in Component[11] IncSg=Incoherent()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[11]);
  /* User FINALLY code for component 'SqwV'. */
  SIG_MESSAGE("SqwV (Finally)");
#define mccompcurname  SqwV
#define mccompcurtype  Isotropic_Sqw
#define mccompcurindex 12
#define powder_format mccSqwV_powder_format
#define VarSqw mccSqwV_VarSqw
#define columns mccSqwV_columns
#define offdata mccSqwV_offdata
{   /* Declarations of SqwV=Isotropic_Sqw() SETTING parameters. */
char* Sqw_coh = mccSqwV_Sqw_coh;
char* Sqw_inc = mccSqwV_Sqw_inc;
char* geometry = mccSqwV_geometry;
MCNUM radius = mccSqwV_radius;
MCNUM thickness = mccSqwV_thickness;
MCNUM xwidth = mccSqwV_xwidth;
MCNUM yheight = mccSqwV_yheight;
MCNUM zdepth = mccSqwV_zdepth;
MCNUM threshold = mccSqwV_threshold;
int order = mccSqwV_order;
MCNUM T = mccSqwV_T;
MCNUM verbose = mccSqwV_verbose;
MCNUM d_phi = mccSqwV_d_phi;
int concentric = mccSqwV_concentric;
MCNUM rho = mccSqwV_rho;
MCNUM sigma_abs = mccSqwV_sigma_abs;
MCNUM sigma_coh = mccSqwV_sigma_coh;
MCNUM sigma_inc = mccSqwV_sigma_inc;
MCNUM classical = mccSqwV_classical;
MCNUM powder_Dd = mccSqwV_powder_Dd;
MCNUM powder_DW = mccSqwV_powder_DW;
MCNUM powder_Vc = mccSqwV_powder_Vc;
MCNUM density = mccSqwV_density;
MCNUM weight = mccSqwV_weight;
MCNUM p_interact = mccSqwV_p_interact;
MCNUM norm = mccSqwV_norm;
MCNUM powder_barns = mccSqwV_powder_barns;
char* quantum_correction = mccSqwV_quantum_correction;
#line 2701 "/usr/share/mcstas/2.6rc1/samples/Isotropic_Sqw.comp"
{
  int  k;

  if (VarSqw.s_coh > 0 || VarSqw.s_inc > 0)
  for (k=0; k < 2; k++) {
    struct Sqw_Data_struct Data_sqw;

    Data_sqw =  (k == 0 ? VarSqw.Data_coh : VarSqw.Data_inc);
    /* Data_sqw->Sqw has already been freed at end of INIT */
    Table_Free(&(Data_sqw.iqSq));

    if (Data_sqw.SW)           free(Data_sqw.SW);
    if (Data_sqw.SQW)          free(Data_sqw.SQW);
    if (Data_sqw.SW_lookup)    free(Data_sqw.SW_lookup);
    if (Data_sqw.QW_lookup)    free(Data_sqw.QW_lookup);
  } /* end for */

#ifdef USE_MPI
  if (mpi_node_count > 1) {
    double tmp;
    tmp = (double)VarSqw.neutron_removed; mc_MPI_Sum(&tmp, 1); VarSqw.neutron_removed=(long)tmp;
    tmp = (double)VarSqw.neutron_exit;    mc_MPI_Sum(&tmp, 1); VarSqw.neutron_exit=(long)tmp;
    tmp = (double)VarSqw.neutron_pmult;   mc_MPI_Sum(&tmp, 1); VarSqw.neutron_pmult=(long)tmp;
    mc_MPI_Sum(&VarSqw.mean_scatt, 1);
    mc_MPI_Sum(&VarSqw.psum_scatt, 1);
    mc_MPI_Sum(&VarSqw.mean_abs, 1);
    mc_MPI_Sum(&VarSqw.single_coh, 1);
    mc_MPI_Sum(&VarSqw.single_inc, 1);
    mc_MPI_Sum(&VarSqw.multi, 1);
  }
#endif
  MPI_MASTER(
  if (VarSqw.neutron_removed)
    printf("Isotropic_Sqw: %s: %li neutron events (out of %li) that should have\n"
           "               scattered were transmitted because scattering conditions\n"
           "WARNING        could not be satisfied after %i tries.\n",
          NAME_CURRENT_COMP, VarSqw.neutron_removed,
          VarSqw.neutron_exit+VarSqw.neutron_removed, VarSqw.maxloop);
  if (VarSqw.neutron_pmult)
    printf("Isotropic_Sqw: %s: %li neutron events (out of %li) reached\n"
           "WARNING        unrealistic weight. The S(q,w) norm might be too high.\n",
          NAME_CURRENT_COMP, VarSqw.neutron_pmult, VarSqw.neutron_exit);

  if (VarSqw.verbose_output >= 1 && VarSqw.psum_scatt > 0) {
    printf("Isotropic_Sqw: %s: Scattering fraction=%g of incoming intensity\n"
           "               Absorption fraction           =%g\n",
           NAME_CURRENT_COMP,
           VarSqw.mean_scatt/VarSqw.psum_scatt, VarSqw.mean_abs/VarSqw.psum_scatt);
    printf("               Single   scattering intensity =%g (coh=%g inc=%g)\n"
           "               Multiple scattering intensity =%g\n",
           VarSqw.single_coh+VarSqw.single_inc, VarSqw.single_coh, VarSqw.single_inc, VarSqw.multi);
    );
  }

/* end FINALLY */
}
#line 26103 "Samples_Incoherent.c"
}   /* End of SqwV=Isotropic_Sqw() SETTING parameter declarations. */
#undef offdata
#undef columns
#undef VarSqw
#undef powder_format
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[12]) fprintf(stderr, "Warning: No neutron could reach Component[12] SqwV\n");
    if (mcAbsorbProp[12]) fprintf(stderr, "Warning: %g events were removed in Component[12] SqwV=Isotropic_Sqw()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[12]);
  /* User FINALLY code for component 'Dirbeam'. */
  SIG_MESSAGE("Dirbeam (Finally)");
#define mccompcurname  Dirbeam
#define mccompcurtype  PSD_monitor
#define mccompcurindex 13
#define PSD_N mccDirbeam_PSD_N
#define PSD_p mccDirbeam_PSD_p
#define PSD_p2 mccDirbeam_PSD_p2
{   /* Declarations of Dirbeam=PSD_monitor() SETTING parameters. */
int nx = mccDirbeam_nx;
int ny = mccDirbeam_ny;
char* filename = mccDirbeam_filename;
MCNUM xmin = mccDirbeam_xmin;
MCNUM xmax = mccDirbeam_xmax;
MCNUM ymin = mccDirbeam_ymin;
MCNUM ymax = mccDirbeam_ymax;
MCNUM xwidth = mccDirbeam_xwidth;
MCNUM yheight = mccDirbeam_yheight;
MCNUM restore_neutron = mccDirbeam_restore_neutron;
#line 122 "/usr/share/mcstas/2.6rc1/monitors/PSD_monitor.comp"
{
  destroy_darr2d(PSD_N);
  destroy_darr2d(PSD_p);
  destroy_darr2d(PSD_p2);
}
#line 26140 "Samples_Incoherent.c"
}   /* End of Dirbeam=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[13]) fprintf(stderr, "Warning: No neutron could reach Component[13] Dirbeam\n");
    if (mcAbsorbProp[13]) fprintf(stderr, "Warning: %g events were removed in Component[13] Dirbeam=PSD_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[13]);
    if (!mcNCounter[14]) fprintf(stderr, "Warning: No neutron could reach Component[14] Stop\n");
    if (mcAbsorbProp[14]) fprintf(stderr, "Warning: %g events were removed in Component[14] Stop=Beamstop()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[14]);
  /* User FINALLY code for component 'Tof_Sphere_4pi'. */
  SIG_MESSAGE("Tof_Sphere_4pi (Finally)");
#define mccompcurname  Tof_Sphere_4pi
#define mccompcurtype  Monitor_nD
#define mccompcurindex 15
#define user1 mccTof_Sphere_4pi_user1
#define user2 mccTof_Sphere_4pi_user2
#define user3 mccTof_Sphere_4pi_user3
#define DEFS mccTof_Sphere_4pi_DEFS
#define Vars mccTof_Sphere_4pi_Vars
#define detector mccTof_Sphere_4pi_detector
#define offdata mccTof_Sphere_4pi_offdata
{   /* Declarations of Tof_Sphere_4pi=Monitor_nD() SETTING parameters. */
MCNUM xwidth = mccTof_Sphere_4pi_xwidth;
MCNUM yheight = mccTof_Sphere_4pi_yheight;
MCNUM zdepth = mccTof_Sphere_4pi_zdepth;
MCNUM xmin = mccTof_Sphere_4pi_xmin;
MCNUM xmax = mccTof_Sphere_4pi_xmax;
MCNUM ymin = mccTof_Sphere_4pi_ymin;
MCNUM ymax = mccTof_Sphere_4pi_ymax;
MCNUM zmin = mccTof_Sphere_4pi_zmin;
MCNUM zmax = mccTof_Sphere_4pi_zmax;
MCNUM bins = mccTof_Sphere_4pi_bins;
MCNUM min = mccTof_Sphere_4pi_min;
MCNUM max = mccTof_Sphere_4pi_max;
MCNUM restore_neutron = mccTof_Sphere_4pi_restore_neutron;
MCNUM radius = mccTof_Sphere_4pi_radius;
char* options = mccTof_Sphere_4pi_options;
char* filename = mccTof_Sphere_4pi_filename;
char* geometry = mccTof_Sphere_4pi_geometry;
char* username1 = mccTof_Sphere_4pi_username1;
char* username2 = mccTof_Sphere_4pi_username2;
char* username3 = mccTof_Sphere_4pi_username3;
int nowritefile = mccTof_Sphere_4pi_nowritefile;
#line 491 "/usr/share/mcstas/2.6rc1/monitors/Monitor_nD.comp"
{
  /* free pointers */
  Monitor_nD_Finally(&DEFS, &Vars);
}
#line 26192 "Samples_Incoherent.c"
}   /* End of Tof_Sphere_4pi=Monitor_nD() SETTING parameter declarations. */
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

    if (!mcNCounter[15]) fprintf(stderr, "Warning: No neutron could reach Component[15] Tof_Sphere_4pi\n");
    if (mcAbsorbProp[15]) fprintf(stderr, "Warning: %g events were removed in Component[15] Tof_Sphere_4pi=Monitor_nD()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[15]);
    if (!mcNCounter[16]) fprintf(stderr, "Warning: No neutron could reach Component[16] PSD_Sphere_4pi\n");
    if (mcAbsorbProp[16]) fprintf(stderr, "Warning: %g events were removed in Component[16] PSD_Sphere_4pi=PSD_monitor_4PI()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[16]);
    if (!mcNCounter[17]) fprintf(stderr, "Warning: No neutron could reach Component[17] Coord3\n");
    if (mcAbsorbProp[17]) fprintf(stderr, "Warning: %g events were removed in Component[17] Coord3=Arm()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[17]);
    if (!mcNCounter[18]) fprintf(stderr, "Warning: No neutron could reach Component[18] Sphere_4pi\n");
    if (mcAbsorbProp[18]) fprintf(stderr, "Warning: %g events were removed in Component[18] Sphere_4pi=PSD_monitor_4PI()\n"
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

  /* MCDISPLAY code for component 'Source'. */
  SIG_MESSAGE("Source (McDisplay)");
  printf("MCDISPLAY: component %s\n", "Source");
#define mccompcurname  Source
#define mccompcurtype  ISIS_moderator
#define mccompcurindex 1
#define Face mccSource_Face
#define p_in mccSource_p_in
#define Tnpts mccSource_Tnpts
#define scaleSize mccSource_scaleSize
#define angleArea mccSource_angleArea
#define Nsim mccSource_Nsim
#define Ncount mccSource_Ncount
#define TS mccSource_TS
#define rtE0 mccSource_rtE0
#define rtE1 mccSource_rtE1
#define rtmodX mccSource_rtmodX
#define rtmodY mccSource_rtmodY
#define TargetStation mccSource_TargetStation
#define CurrentWeight mccSource_CurrentWeight
{   /* Declarations of Source=ISIS_moderator() SETTING parameters. */
MCNUM Emin = mccSource_Emin;
MCNUM Emax = mccSource_Emax;
MCNUM dist = mccSource_dist;
MCNUM focus_xw = mccSource_focus_xw;
MCNUM focus_yh = mccSource_focus_yh;
MCNUM xwidth = mccSource_xwidth;
MCNUM yheight = mccSource_yheight;
MCNUM CAngle = mccSource_CAngle;
MCNUM SAC = mccSource_SAC;
MCNUM Lmin = mccSource_Lmin;
MCNUM Lmax = mccSource_Lmax;
int target_index = mccSource_target_index;
MCNUM verbose = mccSource_verbose;
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
#line 26289 "Samples_Incoherent.c"
}   /* End of Source=ISIS_moderator() SETTING parameter declarations. */
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

  /* MCDISPLAY code for component 'Coord'. */
  SIG_MESSAGE("Coord (McDisplay)");
  printf("MCDISPLAY: component %s\n", "Coord");
#define mccompcurname  Coord
#define mccompcurtype  Progress_bar
#define mccompcurindex 2
#define IntermediateCnts mccCoord_IntermediateCnts
#define StartTime mccCoord_StartTime
#define EndTime mccCoord_EndTime
#define CurrentTime mccCoord_CurrentTime
{   /* Declarations of Coord=Progress_bar() SETTING parameters. */
char* profile = mccCoord_profile;
MCNUM percent = mccCoord_percent;
MCNUM flag_save = mccCoord_flag_save;
MCNUM minutes = mccCoord_minutes;
#line 147 "/usr/share/mcstas/2.6rc1/misc/Progress_bar.comp"
{
  
}
#line 26328 "Samples_Incoherent.c"
}   /* End of Coord=Progress_bar() SETTING parameter declarations. */
#undef CurrentTime
#undef EndTime
#undef StartTime
#undef IntermediateCnts
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'Vsample'. */
  SIG_MESSAGE("Vsample (McDisplay)");
  printf("MCDISPLAY: component %s\n", "Vsample");
#define mccompcurname  Vsample
#define mccompcurtype  V_sample
#define mccompcurindex 3
#define VarsV mccVsample_VarsV
{   /* Declarations of Vsample=V_sample() SETTING parameters. */
MCNUM radius = mccVsample_radius;
MCNUM thickness = mccVsample_thickness;
MCNUM zdepth = mccVsample_zdepth;
MCNUM Vc = mccVsample_Vc;
MCNUM sigma_abs = mccVsample_sigma_abs;
MCNUM sigma_inc = mccVsample_sigma_inc;
MCNUM radius_i = mccVsample_radius_i;
MCNUM radius_o = mccVsample_radius_o;
MCNUM h = mccVsample_h;
MCNUM focus_r = mccVsample_focus_r;
MCNUM pack = mccVsample_pack;
MCNUM frac = mccVsample_frac;
MCNUM f_QE = mccVsample_f_QE;
MCNUM gamma = mccVsample_gamma;
MCNUM target_x = mccVsample_target_x;
MCNUM target_y = mccVsample_target_y;
MCNUM target_z = mccVsample_target_z;
MCNUM focus_xw = mccVsample_focus_xw;
MCNUM focus_yh = mccVsample_focus_yh;
MCNUM focus_aw = mccVsample_focus_aw;
MCNUM focus_ah = mccVsample_focus_ah;
MCNUM xwidth = mccVsample_xwidth;
MCNUM yheight = mccVsample_yheight;
MCNUM zthick = mccVsample_zthick;
MCNUM rad_sphere = mccVsample_rad_sphere;
MCNUM sig_a = mccVsample_sig_a;
MCNUM sig_i = mccVsample_sig_i;
MCNUM V0 = mccVsample_V0;
int target_index = mccVsample_target_index;
MCNUM multiples = mccVsample_multiples;
#line 320 "/usr/share/mcstas/2.6rc1/obsolete/V_sample.comp"
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
#line 26423 "Samples_Incoherent.c"
}   /* End of Vsample=V_sample() SETTING parameter declarations. */
#undef VarsV
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'PowN'. */
  SIG_MESSAGE("PowN (McDisplay)");
  printf("MCDISPLAY: component %s\n", "PowN");
#define mccompcurname  PowN
#define mccompcurtype  PowderN
#define mccompcurindex 4
#define format mccPowN_format
#define line_info mccPowN_line_info
#define columns mccPowN_columns
#define offdata mccPowN_offdata
{   /* Declarations of PowN=PowderN() SETTING parameters. */
char* reflections = mccPowN_reflections;
char* geometry = mccPowN_geometry;
MCNUM radius = mccPowN_radius;
MCNUM yheight = mccPowN_yheight;
MCNUM xwidth = mccPowN_xwidth;
MCNUM zdepth = mccPowN_zdepth;
MCNUM thickness = mccPowN_thickness;
MCNUM pack = mccPowN_pack;
MCNUM Vc = mccPowN_Vc;
MCNUM sigma_abs = mccPowN_sigma_abs;
MCNUM sigma_inc = mccPowN_sigma_inc;
MCNUM delta_d_d = mccPowN_delta_d_d;
MCNUM p_inc = mccPowN_p_inc;
MCNUM p_transmit = mccPowN_p_transmit;
MCNUM DW = mccPowN_DW;
MCNUM nb_atoms = mccPowN_nb_atoms;
MCNUM d_omega = mccPowN_d_omega;
MCNUM d_phi = mccPowN_d_phi;
MCNUM tth_sign = mccPowN_tth_sign;
MCNUM p_interact = mccPowN_p_interact;
MCNUM concentric = mccPowN_concentric;
MCNUM density = mccPowN_density;
MCNUM weight = mccPowN_weight;
MCNUM barns = mccPowN_barns;
MCNUM Strain = mccPowN_Strain;
MCNUM focus_flip = mccPowN_focus_flip;
int target_index = mccPowN_target_index;
#line 1111 "/usr/share/mcstas/2.6rc1/samples/PowderN.comp"
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
#line 26545 "Samples_Incoherent.c"
}   /* End of PowN=PowderN() SETTING parameter declarations. */
#undef offdata
#undef columns
#undef line_info
#undef format
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'SX'. */
  SIG_MESSAGE("SX (McDisplay)");
  printf("MCDISPLAY: component %s\n", "SX");
#define mccompcurname  SX
#define mccompcurtype  Single_crystal
#define mccompcurindex 5
#define mosaic_AB mccSX_mosaic_AB
#define hkl_info mccSX_hkl_info
#define offdata mccSX_offdata
{   /* Declarations of SX=Single_crystal() SETTING parameters. */
char* reflections = mccSX_reflections;
char* geometry = mccSX_geometry;
MCNUM xwidth = mccSX_xwidth;
MCNUM yheight = mccSX_yheight;
MCNUM zdepth = mccSX_zdepth;
MCNUM radius = mccSX_radius;
MCNUM delta_d_d = mccSX_delta_d_d;
MCNUM mosaic = mccSX_mosaic;
MCNUM mosaic_a = mccSX_mosaic_a;
MCNUM mosaic_b = mccSX_mosaic_b;
MCNUM mosaic_c = mccSX_mosaic_c;
MCNUM recip_cell = mccSX_recip_cell;
MCNUM barns = mccSX_barns;
MCNUM ax = mccSX_ax;
MCNUM ay = mccSX_ay;
MCNUM az = mccSX_az;
MCNUM bx = mccSX_bx;
MCNUM by = mccSX_by;
MCNUM bz = mccSX_bz;
MCNUM cx = mccSX_cx;
MCNUM cy = mccSX_cy;
MCNUM cz = mccSX_cz;
MCNUM p_transmit = mccSX_p_transmit;
MCNUM sigma_abs = mccSX_sigma_abs;
MCNUM sigma_inc = mccSX_sigma_inc;
MCNUM aa = mccSX_aa;
MCNUM bb = mccSX_bb;
MCNUM cc = mccSX_cc;
MCNUM order = mccSX_order;
MCNUM RX = mccSX_RX;
MCNUM RY = mccSX_RY;
MCNUM powder = mccSX_powder;
MCNUM PG = mccSX_PG;
MCNUM deltak = mccSX_deltak;
#line 1431 "/usr/share/mcstas/2.6rc1/samples/Single_crystal.comp"
{
  
  if (hkl_info.shape == 0) {        /* cylinder */
    circle("xz", 0,  yheight/2.0, 0, radius);
    circle("xz", 0, -yheight/2.0, 0, radius);
    line(-radius, -yheight/2.0, 0, -radius, +yheight/2.0, 0);
    line(+radius, -yheight/2.0, 0, +radius, +yheight/2.0, 0);
    line(0, -yheight/2.0, -radius, 0, +yheight/2.0, -radius);
    line(0, -yheight/2.0, +radius, 0, +yheight/2.0, +radius);
  }
  else if (hkl_info.shape == 1) {         /* box */
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
  else if (hkl_info.shape == 2) {        /* sphere */
    circle("xy", 0,  0.0, 0, radius);
    circle("xz", 0,  0.0, 0, radius);
    circle("yz", 0,  0.0, 0, radius);
  }
  else if (hkl_info.shape == 3) {        /* OFF file */
    off_display(offdata);
  }
}
#line 26641 "Samples_Incoherent.c"
}   /* End of SX=Single_crystal() SETTING parameter declarations. */
#undef offdata
#undef hkl_info
#undef mosaic_AB
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'Sqw'. */
  SIG_MESSAGE("Sqw (McDisplay)");
  printf("MCDISPLAY: component %s\n", "Sqw");
#define mccompcurname  Sqw
#define mccompcurtype  Isotropic_Sqw
#define mccompcurindex 6
#define powder_format mccSqw_powder_format
#define VarSqw mccSqw_VarSqw
#define columns mccSqw_columns
#define offdata mccSqw_offdata
{   /* Declarations of Sqw=Isotropic_Sqw() SETTING parameters. */
char* Sqw_coh = mccSqw_Sqw_coh;
char* Sqw_inc = mccSqw_Sqw_inc;
char* geometry = mccSqw_geometry;
MCNUM radius = mccSqw_radius;
MCNUM thickness = mccSqw_thickness;
MCNUM xwidth = mccSqw_xwidth;
MCNUM yheight = mccSqw_yheight;
MCNUM zdepth = mccSqw_zdepth;
MCNUM threshold = mccSqw_threshold;
int order = mccSqw_order;
MCNUM T = mccSqw_T;
MCNUM verbose = mccSqw_verbose;
MCNUM d_phi = mccSqw_d_phi;
int concentric = mccSqw_concentric;
MCNUM rho = mccSqw_rho;
MCNUM sigma_abs = mccSqw_sigma_abs;
MCNUM sigma_coh = mccSqw_sigma_coh;
MCNUM sigma_inc = mccSqw_sigma_inc;
MCNUM classical = mccSqw_classical;
MCNUM powder_Dd = mccSqw_powder_Dd;
MCNUM powder_DW = mccSqw_powder_DW;
MCNUM powder_Vc = mccSqw_powder_Vc;
MCNUM density = mccSqw_density;
MCNUM weight = mccSqw_weight;
MCNUM p_interact = mccSqw_p_interact;
MCNUM norm = mccSqw_norm;
MCNUM powder_barns = mccSqw_powder_barns;
char* quantum_correction = mccSqw_quantum_correction;
#line 2761 "/usr/share/mcstas/2.6rc1/samples/Isotropic_Sqw.comp"
{
  if (VarSqw.s_coh > 0 || VarSqw.s_inc > 0) {
    
    if(VarSqw.shape==1)
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

      if (thickness) {
        xmin = -0.5*xwidth+thickness;
        xmax = -xmin;
        ymin = -0.5*yheight+thickness;
        ymax = -ymin;
        zmin = -0.5*zdepth+thickness;
        zmax = -zmin;
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
    else if(VarSqw.shape==0)
    {
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
    } else if(VarSqw.shape==2) {
      if (thickness) {
        double radius_i=radius-thickness;
        circle("xy",0,0,0,radius_i);
        circle("xz",0,0,0,radius_i);
        circle("yz",0,0,0,radius_i);
      }
      circle("xy",0,0,0,radius);
      circle("xz",0,0,0,radius);
      circle("yz",0,0,0,radius);
    } else if (VarSqw.shape == 3) {	/* OFF file */
      off_display(offdata);
    }
  }
/* end MCDISPLAY */
}
#line 26772 "Samples_Incoherent.c"
}   /* End of Sqw=Isotropic_Sqw() SETTING parameter declarations. */
#undef offdata
#undef columns
#undef VarSqw
#undef powder_format
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'IncS'. */
  SIG_MESSAGE("IncS (McDisplay)");
  printf("MCDISPLAY: component %s\n", "IncS");
#define mccompcurname  IncS
#define mccompcurtype  Incoherent
#define mccompcurindex 7
#define VarsInc mccIncS_VarsInc
#define offdata mccIncS_offdata
{   /* Declarations of IncS=Incoherent() SETTING parameters. */
char* geometry = mccIncS_geometry;
MCNUM radius = mccIncS_radius;
MCNUM xwidth = mccIncS_xwidth;
MCNUM yheight = mccIncS_yheight;
MCNUM zdepth = mccIncS_zdepth;
MCNUM thickness = mccIncS_thickness;
MCNUM target_x = mccIncS_target_x;
MCNUM target_y = mccIncS_target_y;
MCNUM target_z = mccIncS_target_z;
MCNUM focus_r = mccIncS_focus_r;
MCNUM focus_xw = mccIncS_focus_xw;
MCNUM focus_yh = mccIncS_focus_yh;
MCNUM focus_aw = mccIncS_focus_aw;
MCNUM focus_ah = mccIncS_focus_ah;
int target_index = mccIncS_target_index;
MCNUM pack = mccIncS_pack;
MCNUM p_interact = mccIncS_p_interact;
MCNUM f_QE = mccIncS_f_QE;
MCNUM gamma = mccIncS_gamma;
MCNUM Etrans = mccIncS_Etrans;
MCNUM deltaE = mccIncS_deltaE;
MCNUM sigma_abs = mccIncS_sigma_abs;
MCNUM sigma_inc = mccIncS_sigma_inc;
MCNUM Vc = mccIncS_Vc;
MCNUM concentric = mccIncS_concentric;
MCNUM order = mccIncS_order;
#line 439 "/usr/share/mcstas/2.6rc1/samples/Incoherent.comp"
{
  
  if (geometry && strlen(geometry) && strcmp(geometry, "NULL") && strcmp(geometry, "0")) {	/* OFF file */
    off_display(offdata);
  }
  else if (radius > 0 &&  yheight) {	/* cylinder */
      cylinder(0,0,0,radius,yheight,0, 0,1,0);
      cylinder(0,0,0,radius-thickness,yheight,0, 0,1,0);
  }
  else if (xwidth && yheight) { 	/* box/rectangle XY */
      box(0,0,0,xwidth, yheight, zdepth);
  }
  else if (radius > 0 && !yheight) {	/* sphere */
      sphere(0,0,0,radius,0);
  }
}
#line 26834 "Samples_Incoherent.c"
}   /* End of IncS=Incoherent() SETTING parameter declarations. */
#undef offdata
#undef VarsInc
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'PowNg'. */
  SIG_MESSAGE("PowNg (McDisplay)");
  printf("MCDISPLAY: component %s\n", "PowNg");
#define mccompcurname  PowNg
#define mccompcurtype  PowderN
#define mccompcurindex 8
#define format mccPowNg_format
#define line_info mccPowNg_line_info
#define columns mccPowNg_columns
#define offdata mccPowNg_offdata
{   /* Declarations of PowNg=PowderN() SETTING parameters. */
char* reflections = mccPowNg_reflections;
char* geometry = mccPowNg_geometry;
MCNUM radius = mccPowNg_radius;
MCNUM yheight = mccPowNg_yheight;
MCNUM xwidth = mccPowNg_xwidth;
MCNUM zdepth = mccPowNg_zdepth;
MCNUM thickness = mccPowNg_thickness;
MCNUM pack = mccPowNg_pack;
MCNUM Vc = mccPowNg_Vc;
MCNUM sigma_abs = mccPowNg_sigma_abs;
MCNUM sigma_inc = mccPowNg_sigma_inc;
MCNUM delta_d_d = mccPowNg_delta_d_d;
MCNUM p_inc = mccPowNg_p_inc;
MCNUM p_transmit = mccPowNg_p_transmit;
MCNUM DW = mccPowNg_DW;
MCNUM nb_atoms = mccPowNg_nb_atoms;
MCNUM d_omega = mccPowNg_d_omega;
MCNUM d_phi = mccPowNg_d_phi;
MCNUM tth_sign = mccPowNg_tth_sign;
MCNUM p_interact = mccPowNg_p_interact;
MCNUM concentric = mccPowNg_concentric;
MCNUM density = mccPowNg_density;
MCNUM weight = mccPowNg_weight;
MCNUM barns = mccPowNg_barns;
MCNUM Strain = mccPowNg_Strain;
MCNUM focus_flip = mccPowNg_focus_flip;
int target_index = mccPowNg_target_index;
#line 1111 "/usr/share/mcstas/2.6rc1/samples/PowderN.comp"
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
#line 26957 "Samples_Incoherent.c"
}   /* End of PowNg=PowderN() SETTING parameter declarations. */
#undef offdata
#undef columns
#undef line_info
#undef format
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'SXg'. */
  SIG_MESSAGE("SXg (McDisplay)");
  printf("MCDISPLAY: component %s\n", "SXg");
#define mccompcurname  SXg
#define mccompcurtype  Single_crystal
#define mccompcurindex 9
#define mosaic_AB mccSXg_mosaic_AB
#define hkl_info mccSXg_hkl_info
#define offdata mccSXg_offdata
{   /* Declarations of SXg=Single_crystal() SETTING parameters. */
char* reflections = mccSXg_reflections;
char* geometry = mccSXg_geometry;
MCNUM xwidth = mccSXg_xwidth;
MCNUM yheight = mccSXg_yheight;
MCNUM zdepth = mccSXg_zdepth;
MCNUM radius = mccSXg_radius;
MCNUM delta_d_d = mccSXg_delta_d_d;
MCNUM mosaic = mccSXg_mosaic;
MCNUM mosaic_a = mccSXg_mosaic_a;
MCNUM mosaic_b = mccSXg_mosaic_b;
MCNUM mosaic_c = mccSXg_mosaic_c;
MCNUM recip_cell = mccSXg_recip_cell;
MCNUM barns = mccSXg_barns;
MCNUM ax = mccSXg_ax;
MCNUM ay = mccSXg_ay;
MCNUM az = mccSXg_az;
MCNUM bx = mccSXg_bx;
MCNUM by = mccSXg_by;
MCNUM bz = mccSXg_bz;
MCNUM cx = mccSXg_cx;
MCNUM cy = mccSXg_cy;
MCNUM cz = mccSXg_cz;
MCNUM p_transmit = mccSXg_p_transmit;
MCNUM sigma_abs = mccSXg_sigma_abs;
MCNUM sigma_inc = mccSXg_sigma_inc;
MCNUM aa = mccSXg_aa;
MCNUM bb = mccSXg_bb;
MCNUM cc = mccSXg_cc;
MCNUM order = mccSXg_order;
MCNUM RX = mccSXg_RX;
MCNUM RY = mccSXg_RY;
MCNUM powder = mccSXg_powder;
MCNUM PG = mccSXg_PG;
MCNUM deltak = mccSXg_deltak;
#line 1431 "/usr/share/mcstas/2.6rc1/samples/Single_crystal.comp"
{
  
  if (hkl_info.shape == 0) {        /* cylinder */
    circle("xz", 0,  yheight/2.0, 0, radius);
    circle("xz", 0, -yheight/2.0, 0, radius);
    line(-radius, -yheight/2.0, 0, -radius, +yheight/2.0, 0);
    line(+radius, -yheight/2.0, 0, +radius, +yheight/2.0, 0);
    line(0, -yheight/2.0, -radius, 0, +yheight/2.0, -radius);
    line(0, -yheight/2.0, +radius, 0, +yheight/2.0, +radius);
  }
  else if (hkl_info.shape == 1) {         /* box */
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
  else if (hkl_info.shape == 2) {        /* sphere */
    circle("xy", 0,  0.0, 0, radius);
    circle("xz", 0,  0.0, 0, radius);
    circle("yz", 0,  0.0, 0, radius);
  }
  else if (hkl_info.shape == 3) {        /* OFF file */
    off_display(offdata);
  }
}
#line 27053 "Samples_Incoherent.c"
}   /* End of SXg=Single_crystal() SETTING parameter declarations. */
#undef offdata
#undef hkl_info
#undef mosaic_AB
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'Sqwg'. */
  SIG_MESSAGE("Sqwg (McDisplay)");
  printf("MCDISPLAY: component %s\n", "Sqwg");
#define mccompcurname  Sqwg
#define mccompcurtype  Isotropic_Sqw
#define mccompcurindex 10
#define powder_format mccSqwg_powder_format
#define VarSqw mccSqwg_VarSqw
#define columns mccSqwg_columns
#define offdata mccSqwg_offdata
{   /* Declarations of Sqwg=Isotropic_Sqw() SETTING parameters. */
char* Sqw_coh = mccSqwg_Sqw_coh;
char* Sqw_inc = mccSqwg_Sqw_inc;
char* geometry = mccSqwg_geometry;
MCNUM radius = mccSqwg_radius;
MCNUM thickness = mccSqwg_thickness;
MCNUM xwidth = mccSqwg_xwidth;
MCNUM yheight = mccSqwg_yheight;
MCNUM zdepth = mccSqwg_zdepth;
MCNUM threshold = mccSqwg_threshold;
int order = mccSqwg_order;
MCNUM T = mccSqwg_T;
MCNUM verbose = mccSqwg_verbose;
MCNUM d_phi = mccSqwg_d_phi;
int concentric = mccSqwg_concentric;
MCNUM rho = mccSqwg_rho;
MCNUM sigma_abs = mccSqwg_sigma_abs;
MCNUM sigma_coh = mccSqwg_sigma_coh;
MCNUM sigma_inc = mccSqwg_sigma_inc;
MCNUM classical = mccSqwg_classical;
MCNUM powder_Dd = mccSqwg_powder_Dd;
MCNUM powder_DW = mccSqwg_powder_DW;
MCNUM powder_Vc = mccSqwg_powder_Vc;
MCNUM density = mccSqwg_density;
MCNUM weight = mccSqwg_weight;
MCNUM p_interact = mccSqwg_p_interact;
MCNUM norm = mccSqwg_norm;
MCNUM powder_barns = mccSqwg_powder_barns;
char* quantum_correction = mccSqwg_quantum_correction;
#line 2761 "/usr/share/mcstas/2.6rc1/samples/Isotropic_Sqw.comp"
{
  if (VarSqw.s_coh > 0 || VarSqw.s_inc > 0) {
    
    if(VarSqw.shape==1)
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

      if (thickness) {
        xmin = -0.5*xwidth+thickness;
        xmax = -xmin;
        ymin = -0.5*yheight+thickness;
        ymax = -ymin;
        zmin = -0.5*zdepth+thickness;
        zmax = -zmin;
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
    else if(VarSqw.shape==0)
    {
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
    } else if(VarSqw.shape==2) {
      if (thickness) {
        double radius_i=radius-thickness;
        circle("xy",0,0,0,radius_i);
        circle("xz",0,0,0,radius_i);
        circle("yz",0,0,0,radius_i);
      }
      circle("xy",0,0,0,radius);
      circle("xz",0,0,0,radius);
      circle("yz",0,0,0,radius);
    } else if (VarSqw.shape == 3) {	/* OFF file */
      off_display(offdata);
    }
  }
/* end MCDISPLAY */
}
#line 27184 "Samples_Incoherent.c"
}   /* End of Sqwg=Isotropic_Sqw() SETTING parameter declarations. */
#undef offdata
#undef columns
#undef VarSqw
#undef powder_format
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'IncSg'. */
  SIG_MESSAGE("IncSg (McDisplay)");
  printf("MCDISPLAY: component %s\n", "IncSg");
#define mccompcurname  IncSg
#define mccompcurtype  Incoherent
#define mccompcurindex 11
#define VarsInc mccIncSg_VarsInc
#define offdata mccIncSg_offdata
{   /* Declarations of IncSg=Incoherent() SETTING parameters. */
char* geometry = mccIncSg_geometry;
MCNUM radius = mccIncSg_radius;
MCNUM xwidth = mccIncSg_xwidth;
MCNUM yheight = mccIncSg_yheight;
MCNUM zdepth = mccIncSg_zdepth;
MCNUM thickness = mccIncSg_thickness;
MCNUM target_x = mccIncSg_target_x;
MCNUM target_y = mccIncSg_target_y;
MCNUM target_z = mccIncSg_target_z;
MCNUM focus_r = mccIncSg_focus_r;
MCNUM focus_xw = mccIncSg_focus_xw;
MCNUM focus_yh = mccIncSg_focus_yh;
MCNUM focus_aw = mccIncSg_focus_aw;
MCNUM focus_ah = mccIncSg_focus_ah;
int target_index = mccIncSg_target_index;
MCNUM pack = mccIncSg_pack;
MCNUM p_interact = mccIncSg_p_interact;
MCNUM f_QE = mccIncSg_f_QE;
MCNUM gamma = mccIncSg_gamma;
MCNUM Etrans = mccIncSg_Etrans;
MCNUM deltaE = mccIncSg_deltaE;
MCNUM sigma_abs = mccIncSg_sigma_abs;
MCNUM sigma_inc = mccIncSg_sigma_inc;
MCNUM Vc = mccIncSg_Vc;
MCNUM concentric = mccIncSg_concentric;
MCNUM order = mccIncSg_order;
#line 439 "/usr/share/mcstas/2.6rc1/samples/Incoherent.comp"
{
  
  if (geometry && strlen(geometry) && strcmp(geometry, "NULL") && strcmp(geometry, "0")) {	/* OFF file */
    off_display(offdata);
  }
  else if (radius > 0 &&  yheight) {	/* cylinder */
      cylinder(0,0,0,radius,yheight,0, 0,1,0);
      cylinder(0,0,0,radius-thickness,yheight,0, 0,1,0);
  }
  else if (xwidth && yheight) { 	/* box/rectangle XY */
      box(0,0,0,xwidth, yheight, zdepth);
  }
  else if (radius > 0 && !yheight) {	/* sphere */
      sphere(0,0,0,radius,0);
  }
}
#line 27246 "Samples_Incoherent.c"
}   /* End of IncSg=Incoherent() SETTING parameter declarations. */
#undef offdata
#undef VarsInc
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'SqwV'. */
  SIG_MESSAGE("SqwV (McDisplay)");
  printf("MCDISPLAY: component %s\n", "SqwV");
#define mccompcurname  SqwV
#define mccompcurtype  Isotropic_Sqw
#define mccompcurindex 12
#define powder_format mccSqwV_powder_format
#define VarSqw mccSqwV_VarSqw
#define columns mccSqwV_columns
#define offdata mccSqwV_offdata
{   /* Declarations of SqwV=Isotropic_Sqw() SETTING parameters. */
char* Sqw_coh = mccSqwV_Sqw_coh;
char* Sqw_inc = mccSqwV_Sqw_inc;
char* geometry = mccSqwV_geometry;
MCNUM radius = mccSqwV_radius;
MCNUM thickness = mccSqwV_thickness;
MCNUM xwidth = mccSqwV_xwidth;
MCNUM yheight = mccSqwV_yheight;
MCNUM zdepth = mccSqwV_zdepth;
MCNUM threshold = mccSqwV_threshold;
int order = mccSqwV_order;
MCNUM T = mccSqwV_T;
MCNUM verbose = mccSqwV_verbose;
MCNUM d_phi = mccSqwV_d_phi;
int concentric = mccSqwV_concentric;
MCNUM rho = mccSqwV_rho;
MCNUM sigma_abs = mccSqwV_sigma_abs;
MCNUM sigma_coh = mccSqwV_sigma_coh;
MCNUM sigma_inc = mccSqwV_sigma_inc;
MCNUM classical = mccSqwV_classical;
MCNUM powder_Dd = mccSqwV_powder_Dd;
MCNUM powder_DW = mccSqwV_powder_DW;
MCNUM powder_Vc = mccSqwV_powder_Vc;
MCNUM density = mccSqwV_density;
MCNUM weight = mccSqwV_weight;
MCNUM p_interact = mccSqwV_p_interact;
MCNUM norm = mccSqwV_norm;
MCNUM powder_barns = mccSqwV_powder_barns;
char* quantum_correction = mccSqwV_quantum_correction;
#line 2761 "/usr/share/mcstas/2.6rc1/samples/Isotropic_Sqw.comp"
{
  if (VarSqw.s_coh > 0 || VarSqw.s_inc > 0) {
    
    if(VarSqw.shape==1)
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

      if (thickness) {
        xmin = -0.5*xwidth+thickness;
        xmax = -xmin;
        ymin = -0.5*yheight+thickness;
        ymax = -ymin;
        zmin = -0.5*zdepth+thickness;
        zmax = -zmin;
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
    else if(VarSqw.shape==0)
    {
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
    } else if(VarSqw.shape==2) {
      if (thickness) {
        double radius_i=radius-thickness;
        circle("xy",0,0,0,radius_i);
        circle("xz",0,0,0,radius_i);
        circle("yz",0,0,0,radius_i);
      }
      circle("xy",0,0,0,radius);
      circle("xz",0,0,0,radius);
      circle("yz",0,0,0,radius);
    } else if (VarSqw.shape == 3) {	/* OFF file */
      off_display(offdata);
    }
  }
/* end MCDISPLAY */
}
#line 27376 "Samples_Incoherent.c"
}   /* End of SqwV=Isotropic_Sqw() SETTING parameter declarations. */
#undef offdata
#undef columns
#undef VarSqw
#undef powder_format
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'Dirbeam'. */
  SIG_MESSAGE("Dirbeam (McDisplay)");
  printf("MCDISPLAY: component %s\n", "Dirbeam");
#define mccompcurname  Dirbeam
#define mccompcurtype  PSD_monitor
#define mccompcurindex 13
#define PSD_N mccDirbeam_PSD_N
#define PSD_p mccDirbeam_PSD_p
#define PSD_p2 mccDirbeam_PSD_p2
{   /* Declarations of Dirbeam=PSD_monitor() SETTING parameters. */
int nx = mccDirbeam_nx;
int ny = mccDirbeam_ny;
char* filename = mccDirbeam_filename;
MCNUM xmin = mccDirbeam_xmin;
MCNUM xmax = mccDirbeam_xmax;
MCNUM ymin = mccDirbeam_ymin;
MCNUM ymax = mccDirbeam_ymax;
MCNUM xwidth = mccDirbeam_xwidth;
MCNUM yheight = mccDirbeam_yheight;
MCNUM restore_neutron = mccDirbeam_restore_neutron;
#line 129 "/usr/share/mcstas/2.6rc1/monitors/PSD_monitor.comp"
{
  multiline(5,
    (double)xmin, (double)ymin, 0.0,
    (double)xmax, (double)ymin, 0.0,
    (double)xmax, (double)ymax, 0.0,
    (double)xmin, (double)ymax, 0.0,
    (double)xmin, (double)ymin, 0.0);
}
#line 27415 "Samples_Incoherent.c"
}   /* End of Dirbeam=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'Stop'. */
  SIG_MESSAGE("Stop (McDisplay)");
  printf("MCDISPLAY: component %s\n", "Stop");
#define mccompcurname  Stop
#define mccompcurtype  Beamstop
#define mccompcurindex 14
{   /* Declarations of Stop=Beamstop() SETTING parameters. */
MCNUM xmin = mccStop_xmin;
MCNUM xmax = mccStop_xmax;
MCNUM ymin = mccStop_ymin;
MCNUM ymax = mccStop_ymax;
MCNUM xwidth = mccStop_xwidth;
MCNUM yheight = mccStop_yheight;
MCNUM radius = mccStop_radius;
#line 72 "/usr/share/mcstas/2.6rc1/optics/Beamstop.comp"
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
#line 27450 "Samples_Incoherent.c"
}   /* End of Stop=Beamstop() SETTING parameter declarations. */
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'Tof_Sphere_4pi'. */
  SIG_MESSAGE("Tof_Sphere_4pi (McDisplay)");
  printf("MCDISPLAY: component %s\n", "Tof_Sphere_4pi");
#define mccompcurname  Tof_Sphere_4pi
#define mccompcurtype  Monitor_nD
#define mccompcurindex 15
#define user1 mccTof_Sphere_4pi_user1
#define user2 mccTof_Sphere_4pi_user2
#define user3 mccTof_Sphere_4pi_user3
#define DEFS mccTof_Sphere_4pi_DEFS
#define Vars mccTof_Sphere_4pi_Vars
#define detector mccTof_Sphere_4pi_detector
#define offdata mccTof_Sphere_4pi_offdata
{   /* Declarations of Tof_Sphere_4pi=Monitor_nD() SETTING parameters. */
MCNUM xwidth = mccTof_Sphere_4pi_xwidth;
MCNUM yheight = mccTof_Sphere_4pi_yheight;
MCNUM zdepth = mccTof_Sphere_4pi_zdepth;
MCNUM xmin = mccTof_Sphere_4pi_xmin;
MCNUM xmax = mccTof_Sphere_4pi_xmax;
MCNUM ymin = mccTof_Sphere_4pi_ymin;
MCNUM ymax = mccTof_Sphere_4pi_ymax;
MCNUM zmin = mccTof_Sphere_4pi_zmin;
MCNUM zmax = mccTof_Sphere_4pi_zmax;
MCNUM bins = mccTof_Sphere_4pi_bins;
MCNUM min = mccTof_Sphere_4pi_min;
MCNUM max = mccTof_Sphere_4pi_max;
MCNUM restore_neutron = mccTof_Sphere_4pi_restore_neutron;
MCNUM radius = mccTof_Sphere_4pi_radius;
char* options = mccTof_Sphere_4pi_options;
char* filename = mccTof_Sphere_4pi_filename;
char* geometry = mccTof_Sphere_4pi_geometry;
char* username1 = mccTof_Sphere_4pi_username1;
char* username2 = mccTof_Sphere_4pi_username2;
char* username3 = mccTof_Sphere_4pi_username3;
int nowritefile = mccTof_Sphere_4pi_nowritefile;
#line 497 "/usr/share/mcstas/2.6rc1/monitors/Monitor_nD.comp"
{
  if (geometry && strlen(geometry) && strcmp(geometry,"0") && strcmp(geometry, "NULL"))
  {
    off_display(offdata);
  } else {
    Monitor_nD_McDisplay(&DEFS, &Vars);
  }
}
#line 27500 "Samples_Incoherent.c"
}   /* End of Tof_Sphere_4pi=Monitor_nD() SETTING parameter declarations. */
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

  /* MCDISPLAY code for component 'PSD_Sphere_4pi'. */
  SIG_MESSAGE("PSD_Sphere_4pi (McDisplay)");
  printf("MCDISPLAY: component %s\n", "PSD_Sphere_4pi");
#define mccompcurname  PSD_Sphere_4pi
#define mccompcurtype  PSD_monitor_4PI
#define mccompcurindex 16
#define nx mccPSD_Sphere_4pi_nx
#define ny mccPSD_Sphere_4pi_ny
#define PSD_N mccPSD_Sphere_4pi_PSD_N
#define PSD_p mccPSD_Sphere_4pi_PSD_p
#define PSD_p2 mccPSD_Sphere_4pi_PSD_p2
{   /* Declarations of PSD_Sphere_4pi=PSD_monitor_4PI() SETTING parameters. */
char* filename = mccPSD_Sphere_4pi_filename;
MCNUM radius = mccPSD_Sphere_4pi_radius;
MCNUM restore_neutron = mccPSD_Sphere_4pi_restore_neutron;
int nowritefile = mccPSD_Sphere_4pi_nowritefile;
#line 121 "/usr/share/mcstas/2.6rc1/monitors/PSD_monitor_4PI.comp"
{
  
  circle("xy",0,0,0,radius);
  circle("xz",0,0,0,radius);
  circle("yz",0,0,0,radius);
}
#line 27536 "Samples_Incoherent.c"
}   /* End of PSD_Sphere_4pi=PSD_monitor_4PI() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef ny
#undef nx
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'Coord3'. */
  SIG_MESSAGE("Coord3 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "Coord3");
#define mccompcurname  Coord3
#define mccompcurtype  Arm
#define mccompcurindex 17
#line 40 "/usr/share/mcstas/2.6rc1/optics/Arm.comp"
{
  /* A bit ugly; hard-coded dimensions. */
  
  line(0,0,0,0.2,0,0);
  line(0,0,0,0,0.2,0);
  line(0,0,0,0,0,0.2);
}
#line 27561 "Samples_Incoherent.c"
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'Sphere_4pi'. */
  SIG_MESSAGE("Sphere_4pi (McDisplay)");
  printf("MCDISPLAY: component %s\n", "Sphere_4pi");
#define mccompcurname  Sphere_4pi
#define mccompcurtype  PSD_monitor_4PI
#define mccompcurindex 18
#define nx mccSphere_4pi_nx
#define ny mccSphere_4pi_ny
#define PSD_N mccSphere_4pi_PSD_N
#define PSD_p mccSphere_4pi_PSD_p
#define PSD_p2 mccSphere_4pi_PSD_p2
{   /* Declarations of Sphere_4pi=PSD_monitor_4PI() SETTING parameters. */
char* filename = mccSphere_4pi_filename;
MCNUM radius = mccSphere_4pi_radius;
MCNUM restore_neutron = mccSphere_4pi_restore_neutron;
int nowritefile = mccSphere_4pi_nowritefile;
#line 121 "/usr/share/mcstas/2.6rc1/monitors/PSD_monitor_4PI.comp"
{
  
  circle("xy",0,0,0,radius);
  circle("xz",0,0,0,radius);
  circle("yz",0,0,0,radius);
}
#line 27589 "Samples_Incoherent.c"
}   /* End of Sphere_4pi=PSD_monitor_4PI() SETTING parameter declarations. */
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
/* end of generated C code Samples_Incoherent.c */
