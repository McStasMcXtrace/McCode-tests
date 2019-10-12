/* Automatically generated file. Do not edit. 
 * Format:     ANSI C source code
 * Creator:    McStas <http://www.mcstas.org>
 * Instrument: ISIS_MERLIN.instr (ISIS_MERLIN)
 * Date:       Fri Oct 11 20:26:33 2019
 * File:       ISIS_MERLIN.c
 * Compile:    cc -o ISIS_MERLIN.out ISIS_MERLIN.c 
 * CFLAGS=
 */


#define MCCODE_STRING "McStas 2.5 - Dec. 12, 2018"
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



/* the version string is replaced when building distribution with mkdist */
#ifndef MCCODE_STRING
#define MCCODE_STRING "McStas 2.5 - Dec. 12, 2018"
#endif

#ifndef MCCODE_DATE
#define MCCODE_DATE "Dec. 12, 2018"
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
#  define PI 3.14159265358979323846
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
mcstatic inline void vec_prod_func(double *x, double *y, double *z,
		double x1, double y1, double z1, double x2, double y2, double z2);

mcstatic inline double scalar_prod(
		double x1, double y1, double z1, double x2, double y2, double z2);

#define NORM(x,y,z) \
	norm_func(&x, &y, &z)
mcstatic inline void norm_func(double *x, double *y, double *z) {
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
mcstatic inline void normal_vec_func(double *nx, double *ny, double *nz,
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
mcstatic inline void coords_norm(Coords* c);

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

#line 695 "ISIS_MERLIN.c"

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

#line 928 "ISIS_MERLIN.c"

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

  m=abs(m); n=abs(n); p=abs(p); /* make sure dimensions are positive */
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
* mcruninfo_out: output simulation tags/info (both in SIM and data files)
* Used in: mcsiminfo_init (ascii case), mcdetector_out_xD_ascii
*******************************************************************************/
static void mcruninfo_out(char *pre, FILE *f)
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
      if (mcinputtable[i].par){
	/* Parameters with a default value */
	if(mcinputtable[i].val && strlen(mcinputtable[i].val)){
	  (*mcinputtypes[mcinputtable[i].type].printer)(Parameters, mcinputtable[i].par);
	  fprintf(f, "%sParam: %s=%s\n", pre, mcinputtable[i].name, Parameters);
        /* ... and those without */
	}else{
	  fprintf(f, "%sParam: %s=NULL\n", pre, mcinputtable[i].name);
	}
      } 
  }
} /* mcruninfo_out */

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
  if (xl && strlen(xl)) { strncpy(xvar, xl, CHAR_BUF_LENGTH); xvar[2]='\0'; }
  else strcpy(xvar, "x");
  if (yl && strlen(yl)) { strncpy(yvar, yl, CHAR_BUF_LENGTH); yvar[2]='\0'; }
  else strcpy(yvar, "y");

  MCDETECTOR detector;

  /* import and perform basic detector analysis (and handle MPI_Reduce) */
  if (abs(m) == 1) {/* n>1 on Y, m==1 on X: 1D, no X axis*/
    detector = mcdetector_import(mcformat,
      c, (t ? t : MCCODE_STRING " 1D data"),
      n, 1, 1,
      yl, "", "Signal per bin",
      yvar, "(I,Ierr)", "I",
      y1, y2, x1, x2, 0, 0, f,
      p0, p1, p2, posa); /* write Detector: line */
  } else if (abs(n)==1) {/* m>1 on X, n==1 on Y: 1D, no Y axis*/
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
                  1,abs(m),1,abs(n),
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
  mcruninfo_out("  ", stdout);
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

mcstatic inline void coords_norm(Coords* c) {
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
mcstatic inline void vec_prod_func(double *x, double *y, double *z,
		double x1, double y1, double z1,
		double x2, double y2, double z2) {
    *x = (y1)*(z2) - (y2)*(z1);
    *y = (z1)*(x2) - (z2)*(x1);
    *z = (x1)*(y2) - (x2)*(y1);
}

/**
 * Scalar product: use coords_sp for Coords.
 */
mcstatic inline double scalar_prod(
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
mcstatic inline void normal_vec_func(double *nx, double *ny, double *nz,
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

#line 4947 "ISIS_MERLIN.c"

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

#line 5307 "ISIS_MERLIN.c"
#ifdef MC_TRACE_ENABLED
int mctraceenabled = 1;
#else
int mctraceenabled = 0;
#endif
#define MCSTAS "/usr/share/mcstas/2.5/"
int mcdefaultmain = 1;
char mcinstrument_name[] = "ISIS_MERLIN";
char mcinstrument_source[] = "ISIS_MERLIN.instr";
char *mcinstrument_exe=NULL; /* will be set to argv[0] in main */
int main(int argc, char *argv[]){return mccode_main(argc, argv);}
void mcinit(void);
void mcraytrace(void);
void mcsave(FILE *);
void mcfinally(void);
void mcdisplay(void);

/* Shared user declarations for all components 'Guide_channeled'. */
#line 76 "/usr/share/mcstas/2.5/optics/Guide_channeled.comp"
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
            tr->table_ref=(t_Table *)calloc(1,sizeof(t_Table));
            /*copy the contents of the table handle*/
            *(tr->table_ref)= *((t_Table *) item);
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
    Table_File_List_Handler(STORE,tab,0);
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

      if (!hfile && mcinstrument_source && strlen(mcinstrument_source)) /* search in instrument source location */
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
      if (!hfile && mcinstrument_exe && strlen(mcinstrument_exe)) /* search in PWD instrument executable location */
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
                    malloc_size = count_in_array+CHAR_BUF_LENGTH;
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
        (Table->filename ? Table->filename : ""),
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
      Table.filename ? Table.filename : "", buffer);
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

#line 6911 "ISIS_MERLIN.c"

/* Shared user declarations for all components 'FermiChopper'. */
#line 97 "/usr/share/mcstas/2.5/optics/FermiChopper.comp"

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

#line 7091 "ISIS_MERLIN.c"

/* Shared user declarations for all components 'Incoherent'. */
#line 128 "/usr/share/mcstas/2.5/samples/Incoherent.comp"

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
      fprintf(stderr,"Error: First line in %s is too long (=%d). Possibly the line is not terminated by '\\n'.\n" 
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
      fprintf(stderr, "Error: can not read [xyz] coordinates for vertex %li in file %s (interoff/off_init). Read %i values.\n", 
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
      fprintf(stderr, "Error: can not read polygon %i length in file %s (interoff/off_init)\n", 
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
    sprintf(pixelinfo, "%u,%u,%u,%i,%g,%g,%g,%g,%g,%g", data.mantidoffset+pixel, data.mantidoffset, data.mantidoffset+data.polySize-1, nbVertex, cmx, cmy, cmz, x1-cmx, y1-cmy, z1-cmz);
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
#line 8101 "ISIS_MERLIN.c"

/* Instrument parameters. */
MCNUM mcipE_min;
MCNUM mcipE_max;
MCNUM mcipm;
MCNUM mcipE_foc;
MCNUM mcipnu_hz;
int mciptype;

#define mcNUMIPAR 6
int mcnumipar = 6;
struct mcinputtable_struct mcinputtable[mcNUMIPAR+1] = {
  "E_min", &mcipE_min, instr_type_double, "1.0", 
  "E_max", &mcipE_max, instr_type_double, "2000.0", 
  "m", &mcipm, instr_type_double, "3.0", 
  "E_foc", &mcipE_foc, instr_type_double, "80.0", 
  "nu_hz", &mcipnu_hz, instr_type_double, "200.0", 
  "type", &mciptype, instr_type_int, "3", 
  NULL, NULL, instr_type_double, ""
};

/* User declarations from instrument definition. */
#define mccompcurname  ISIS_MERLIN
#define mccompcurtype  INSTRUMENT
#define mccompcurindex 0
#define mcposaISIS_MERLIN coords_set(0,0,0)
#define E_min mcipE_min
#define E_max mcipE_max
#define m mcipm
#define E_foc mcipE_foc
#define nu_hz mcipnu_hz
#define type mciptype
#line 49 "ISIS_MERLIN.instr"
double slit_curv,num_slits,width,len;
double lam,phase_time,t_det,t_fer;
#line 8137 "ISIS_MERLIN.c"
#undef type
#undef nu_hz
#undef E_foc
#undef m
#undef E_max
#undef E_min
#undef mcposaISIS_MERLIN
#undef mccompcurindex
#undef mccompcurtype
#undef mccompcurname

/* neutron state table at each component input (local coords) */
/* [x, y, z, vx, vy, vz, t, sx, sy, sz, p] */
MCNUM mccomp_storein[11*17];
/* Components position table (absolute and relative coords) */
Coords mccomp_posa[17];
Coords mccomp_posr[17];
/* Counter for each comp to check for inactive ones */
MCNUM  mcNCounter[17];
MCNUM  mcPCounter[17];
MCNUM  mcP2Counter[17];
#define mcNUMCOMP 16 /* number of components */
/* Counter for PROP ABSORB */
MCNUM  mcAbsorbProp[17];
/* Flag true when previous component acted on the neutron (SCATTER) */
MCNUM mcScattered=0;
/* Flag true when neutron should be restored (RESTORE) */
MCNUM mcRestore=0;
/* Declarations of component definition and setting parameters. */

/* Setting parameters for component 'isis_mod' [2]. */
char mccisis_mod_Face[16384];
MCNUM mccisis_mod_E0;
MCNUM mccisis_mod_E1;
MCNUM mccisis_mod_modPosition;
MCNUM mccisis_mod_xw;
MCNUM mccisis_mod_yh;
MCNUM mccisis_mod_focus_xw;
MCNUM mccisis_mod_focus_yh;
MCNUM mccisis_mod_dist;
int mccisis_mod_verbose;
MCNUM mccisis_mod_beamcurrent;

/* Definition parameters for component 'TOF_mod' [3]. */
#define mccTOF_mod_nt 150
/* Setting parameters for component 'TOF_mod' [3]. */
char mccTOF_mod_filename[16384];
MCNUM mccTOF_mod_xmin;
MCNUM mccTOF_mod_xmax;
MCNUM mccTOF_mod_ymin;
MCNUM mccTOF_mod_ymax;
MCNUM mccTOF_mod_xwidth;
MCNUM mccTOF_mod_yheight;
MCNUM mccTOF_mod_tmin;
MCNUM mccTOF_mod_tmax;
MCNUM mccTOF_mod_dt;
MCNUM mccTOF_mod_restore_neutron;
int mccTOF_mod_nowritefile;

/* Setting parameters for component 'shutter_guide' [4]. */
MCNUM mccshutter_guide_w1;
MCNUM mccshutter_guide_h1;
MCNUM mccshutter_guide_w2;
MCNUM mccshutter_guide_h2;
MCNUM mccshutter_guide_l;
MCNUM mccshutter_guide_R0;
MCNUM mccshutter_guide_Qc;
MCNUM mccshutter_guide_alpha;
MCNUM mccshutter_guide_m;
MCNUM mccshutter_guide_nslit;
MCNUM mccshutter_guide_d;
MCNUM mccshutter_guide_Qcx;
MCNUM mccshutter_guide_Qcy;
MCNUM mccshutter_guide_alphax;
MCNUM mccshutter_guide_alphay;
MCNUM mccshutter_guide_W;
MCNUM mccshutter_guide_mx;
MCNUM mccshutter_guide_my;
MCNUM mccshutter_guide_nu;
MCNUM mccshutter_guide_phase;

/* Setting parameters for component 'conv_guide' [5]. */
MCNUM mccconv_guide_w1;
MCNUM mccconv_guide_h1;
MCNUM mccconv_guide_w2;
MCNUM mccconv_guide_h2;
MCNUM mccconv_guide_l;
MCNUM mccconv_guide_R0;
MCNUM mccconv_guide_Qc;
MCNUM mccconv_guide_alpha;
MCNUM mccconv_guide_m;
MCNUM mccconv_guide_nslit;
MCNUM mccconv_guide_d;
MCNUM mccconv_guide_Qcx;
MCNUM mccconv_guide_Qcy;
MCNUM mccconv_guide_alphax;
MCNUM mccconv_guide_alphay;
MCNUM mccconv_guide_W;
MCNUM mccconv_guide_mx;
MCNUM mccconv_guide_my;
MCNUM mccconv_guide_nu;
MCNUM mccconv_guide_phase;

/* Setting parameters for component 'between_guide' [6]. */
MCNUM mccbetween_guide_w1;
MCNUM mccbetween_guide_h1;
MCNUM mccbetween_guide_w2;
MCNUM mccbetween_guide_h2;
MCNUM mccbetween_guide_l;
MCNUM mccbetween_guide_R0;
MCNUM mccbetween_guide_Qc;
MCNUM mccbetween_guide_alpha;
MCNUM mccbetween_guide_m;
MCNUM mccbetween_guide_nslit;
MCNUM mccbetween_guide_d;
MCNUM mccbetween_guide_Qcx;
MCNUM mccbetween_guide_Qcy;
MCNUM mccbetween_guide_alphax;
MCNUM mccbetween_guide_alphay;
MCNUM mccbetween_guide_W;
MCNUM mccbetween_guide_mx;
MCNUM mccbetween_guide_my;
MCNUM mccbetween_guide_nu;
MCNUM mccbetween_guide_phase;

/* Setting parameters for component 'fermi' [7]. */
MCNUM mccfermi_phase;
MCNUM mccfermi_radius;
MCNUM mccfermi_nu;
MCNUM mccfermi_w;
MCNUM mccfermi_nslit;
MCNUM mccfermi_R0;
MCNUM mccfermi_Qc;
MCNUM mccfermi_alpha;
MCNUM mccfermi_m;
MCNUM mccfermi_W;
MCNUM mccfermi_length;
MCNUM mccfermi_eff;
MCNUM mccfermi_zero_time;
MCNUM mccfermi_xwidth;
MCNUM mccfermi_verbose;
MCNUM mccfermi_yheight;
MCNUM mccfermi_curvature;
MCNUM mccfermi_delay;

/* Definition parameters for component 'TOF_fermi' [8]. */
#define mccTOF_fermi_nt 50
/* Setting parameters for component 'TOF_fermi' [8]. */
char mccTOF_fermi_filename[16384];
MCNUM mccTOF_fermi_xmin;
MCNUM mccTOF_fermi_xmax;
MCNUM mccTOF_fermi_ymin;
MCNUM mccTOF_fermi_ymax;
MCNUM mccTOF_fermi_xwidth;
MCNUM mccTOF_fermi_yheight;
MCNUM mccTOF_fermi_tmin;
MCNUM mccTOF_fermi_tmax;
MCNUM mccTOF_fermi_dt;
MCNUM mccTOF_fermi_restore_neutron;
int mccTOF_fermi_nowritefile;

/* Setting parameters for component 'final_Guide' [9]. */
MCNUM mccfinal_Guide_w1;
MCNUM mccfinal_Guide_h1;
MCNUM mccfinal_Guide_w2;
MCNUM mccfinal_Guide_h2;
MCNUM mccfinal_Guide_l;
MCNUM mccfinal_Guide_R0;
MCNUM mccfinal_Guide_Qc;
MCNUM mccfinal_Guide_alpha;
MCNUM mccfinal_Guide_m;
MCNUM mccfinal_Guide_nslit;
MCNUM mccfinal_Guide_d;
MCNUM mccfinal_Guide_Qcx;
MCNUM mccfinal_Guide_Qcy;
MCNUM mccfinal_Guide_alphax;
MCNUM mccfinal_Guide_alphay;
MCNUM mccfinal_Guide_W;
MCNUM mccfinal_Guide_mx;
MCNUM mccfinal_Guide_my;
MCNUM mccfinal_Guide_nu;
MCNUM mccfinal_Guide_phase;

/* Definition parameters for component 'horz_div' [10]. */
#define mcchorz_div_nh 30
/* Setting parameters for component 'horz_div' [10]. */
char mcchorz_div_filename[16384];
MCNUM mcchorz_div_xmin;
MCNUM mcchorz_div_xmax;
MCNUM mcchorz_div_ymin;
MCNUM mcchorz_div_ymax;
MCNUM mcchorz_div_xwidth;
MCNUM mcchorz_div_yheight;
MCNUM mcchorz_div_h_maxdiv;
MCNUM mcchorz_div_restore_neutron;
int mcchorz_div_nowritefile;

/* Definition parameters for component 'PSD' [11]. */
#define mccPSD_nx 20
#define mccPSD_ny 20
/* Setting parameters for component 'PSD' [11]. */
char mccPSD_filename[16384];
MCNUM mccPSD_xmin;
MCNUM mccPSD_xmax;
MCNUM mccPSD_ymin;
MCNUM mccPSD_ymax;
MCNUM mccPSD_xwidth;
MCNUM mccPSD_yheight;
MCNUM mccPSD_restore_neutron;
int mccPSD_nowritefile;

/* Definition parameters for component 'det_lmon' [12]. */
#define mccdet_lmon_nL 100
/* Setting parameters for component 'det_lmon' [12]. */
char mccdet_lmon_filename[16384];
MCNUM mccdet_lmon_xmin;
MCNUM mccdet_lmon_xmax;
MCNUM mccdet_lmon_ymin;
MCNUM mccdet_lmon_ymax;
MCNUM mccdet_lmon_xwidth;
MCNUM mccdet_lmon_yheight;
MCNUM mccdet_lmon_Lmin;
MCNUM mccdet_lmon_Lmax;
MCNUM mccdet_lmon_restore_neutron;
int mccdet_lmon_nowritefile;

/* Setting parameters for component 'sample' [13]. */
char mccsample_geometry[16384];
MCNUM mccsample_radius;
MCNUM mccsample_xwidth;
MCNUM mccsample_yheight;
MCNUM mccsample_zdepth;
MCNUM mccsample_thickness;
MCNUM mccsample_target_x;
MCNUM mccsample_target_y;
MCNUM mccsample_target_z;
MCNUM mccsample_focus_r;
MCNUM mccsample_focus_xw;
MCNUM mccsample_focus_yh;
MCNUM mccsample_focus_aw;
MCNUM mccsample_focus_ah;
int mccsample_target_index;
MCNUM mccsample_pack;
MCNUM mccsample_p_interact;
MCNUM mccsample_f_QE;
MCNUM mccsample_gamma;
MCNUM mccsample_Etrans;
MCNUM mccsample_deltaE;
MCNUM mccsample_sigma_abs;
MCNUM mccsample_sigma_inc;
MCNUM mccsample_Vc;
MCNUM mccsample_concentric;
MCNUM mccsample_order;

/* Definition parameters for component 'detector' [15]. */
#define mccdetector_nt 25
/* Setting parameters for component 'detector' [15]. */
char mccdetector_filename[16384];
MCNUM mccdetector_xmin;
MCNUM mccdetector_xmax;
MCNUM mccdetector_ymin;
MCNUM mccdetector_ymax;
MCNUM mccdetector_xwidth;
MCNUM mccdetector_yheight;
MCNUM mccdetector_tmin;
MCNUM mccdetector_tmax;
MCNUM mccdetector_dt;
MCNUM mccdetector_restore_neutron;
int mccdetector_nowritefile;

/* User component declarations. */

/* User declarations for component 'Origin' [1]. */
#define mccompcurname  Origin
#define mccompcurtype  Arm
#define mccompcurindex 1
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'isis_mod' [2]. */
#define mccompcurname  isis_mod
#define mccompcurtype  ViewModISIS
#define mccompcurindex 2
#define Face mccisis_mod_Face
#define E0 mccisis_mod_E0
#define E1 mccisis_mod_E1
#define modPosition mccisis_mod_modPosition
#define xw mccisis_mod_xw
#define yh mccisis_mod_yh
#define focus_xw mccisis_mod_focus_xw
#define focus_yh mccisis_mod_focus_yh
#define dist mccisis_mod_dist
#define verbose mccisis_mod_verbose
#define beamcurrent mccisis_mod_beamcurrent
#line 79 "/usr/share/mcstas/2.5/contrib/ViewModISIS.comp"

#include <ctype.h>

typedef struct
{
  int nEnergy;        ///< Number of energy bins
  int nTime;          ///< number of time bins

  double XAxis;
  double ZAxis;

  double rdumMid;     ///< tally time Window mid point
  double timeOffset;  ///< Time separation
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
  int findRDUM(char*);
  int findTimeOffset(char*,double*);
  FILE* openFile(char*);
  FILE* openFileTest(char*);
  void processHeader(FILE*);
  void readHtable(FILE*,const double,const double);
  int timeStart(char*);
  int timeEnd(char*);
  int energyBin(char*,double,double,double*,double*);
  int notComment(char*);
  void orderEnergy(double*,double*);
  void cutToNumber(char*);
  double convertEnergy(double);
  double EtoLambda(double);
  double calcRDum(double*);

  double strArea();

  /* global variables */

  double p_in;            /* Polorization term (from McSTAS) */
  int Tnpts;              /* Number of points in parameteriation */

  double scaleSize;        /* correction for the actual area of the moderator viewed */

  double angleArea;       /* Area seen by the window  */

  double Nsim;		/* Total number of neutrons to be simulated */

  int Ncount;          /* Number of neutron simulate so far*/
  Source TS;

  /* runtime variables*/
  double fullAngle;
  double rtE0,rtE1;       /* runtime Energy minima and maxima so we can use angstroms as negative input */


double**
matrix(const int m,const int n)
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


double
polInterp(double* X,double* Y,int Psize,double Aim)
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

int
binSearch(int Npts,double* AR,double V)
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

int
cmdnumberD(char *mc,double* num)
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

int
notComment(char* Line)
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

void
orderEnergy(double* A,double *B)
{
  double tmp;
  if (*A>*B)
    {
      tmp=*A;
      *A=*B;
      *B=tmp;
    }
  return;
}

int
timeStart(char* Line)
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

void
cutToNumber(char* Line)
{
  int len,i;

  len=strlen(Line);
  for(i=0;i<len && !isdigit(Line[i]) && Line[i]!='-';i++)
    Line[i]=' ';
  return;
}

int
findTimeOffset(char* Line,double* TO)
  /*!
    Determine the time offset from the file
    \return 1 on success
   */
{
  int len,i;
  double D;

  len=strlen(Line);
  for(i=0;i<len && Line[i]!='T';i++);
  if (len-i<11) return 0;

  if (strncmp(Line+i,"TimeOffset",11) &&
      cmdnumberD(Line+i+11,&D))
    {
      *TO=D/100.0;
      return 1;
    }
  return 0;
}

int
findRDUM(char* Line)
 /*!
   Search for a word rdum at the start of
   the line.
   \param Line :: Line to search
   \returns 1 on success 0 on failure
 */
{
  int len,i;

  len=strlen(Line);
  for(i=0;i<len && Line[i]!='r';i++);
  if (len-i<4) return 0;
  return (strncmp(Line+i,"rdum",4)) ? 0 : 1;
}

double 
convertEnergy(double E)
  /*!
    Convert the energy from eV [not change] or 
    from negative -ve which is angstrom
  */
{
  return (E>0.0) ? E : 81.793936/(E*E);
}

double 
EtoLambda(double E)
  /*!
    Convert the energy from eV [not change] o
    to lambda [A]
  */
{
  return sqrt(81.793936/E);
}
  

int
timeEnd(char* Line)
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

int
energyBin(char* Line,double Einit,double Eend,double* Ea,double* Eb)
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

double
calcFraction(double EI,double EE,double Ea,double Eb)
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

double
calcRDum(double* RPts)
  /*!
    Caluclate the mean distance from the window centre point
   */
{
  double sum,zMid;
  int i;
  extern Source TS;

  sum=0.0;
  for(i=0;i<4;i++)
    {
      fprintf(stderr,"RDUM %g %g %g\n",RPts[i*3],RPts[i*3+1],RPts[i*3+2]);
      sum+=sqrt(RPts[i*3]*RPts[i*3]+
		RPts[i*3+1]*RPts[i*3+1]);
    }
  sum/=400.0;   // Convert to metres from cm
  return sum;
}

void
processHeader(FILE* TFile)
  /*!
    Determine the window and the time offset
    \param TFile :: file [rewound on exit]
  */
{
  char ss[255];          /* BIG space for line */
  int rdumCnt;
  double Ea,Eb;
  double RPts[12];
  int timeOffsetFlag;
  int i,j;
  double D;
  

  extern Source TS;
  
  rdumCnt=0;
  timeOffsetFlag=0;

  while(fgets(ss,255,TFile))
    {
      if (!timeOffsetFlag)
	timeOffsetFlag=findTimeOffset(ss,&TS.timeOffset);

      if (!rdumCnt)
	rdumCnt=findRDUM(ss);
      if (rdumCnt && rdumCnt<5)
	{
	  cutToNumber(ss);
	  for(i=0;i<3 && cmdnumberD(ss,&D);i++)
	    RPts[(rdumCnt-1)*3+i]=D;
	  rdumCnt++;
	}
      // EXIT CONDITION:
      if (rdumCnt*timeOffsetFlag==5)
	{
	  for(j=0;j<3;j++)
	    {
	      TS.ZAxis+=(RPts[3+j]-RPts[j])*(RPts[3+j]-RPts[j]);
	      TS.XAxis+=(RPts[6+j]-RPts[3+j])*(RPts[6+j]-RPts[3+j]);
	    }	


	  TS.ZAxis=sqrt(TS.ZAxis)/100.0;   // convert to metres from cm
	  TS.XAxis=sqrt(TS.XAxis)/100.0;
              if (!modPosition)
                 {
	          TS.ZAxis=yh;
	          TS.XAxis=xw;
                  }
	  fprintf(stderr,"Time off sec == %g \n",TS.timeOffset);
	  fprintf(stderr,"mod size z == %g \n",TS.ZAxis);
	  // TS.rdumMid=calcRDum(RPts);
	  TS.rdumMid=TS.timeOffset; // Goran
	  return;
	}
    }
  return;
}



void
readHtable(FILE* TFile,const double Einit,const double Eend)
/*!
  Process a general h.o file to create an integrated
  table of results from Einit -> Eend
  \param Einit :: inital Energy
  \param Eend  :: final energy
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

  extern Source TS;

  int i;
  /*!
    Status Flag::
    Ftime=1 :: [time ] Reading Time : Data : Err [Exit on Total]

    Double Read File to determine how many bins and
    memery size
  */
  Ea=0.0;
  Eb=0.0;
  fprintf(stderr,"Energy == %g %g\n",Einit,Eend);
  eIndex= -1;
  Ftime=0;
  tIndex=0;
  TS.nTime=0;
  TS.nEnergy=0;
  
  processHeader(TFile);

  // Read file and get time bins:
  while(fgets(ss,255,TFile) && Eend>Ea)
    {
      if (notComment(ss))
	{
          if (!Ftime)
	    {
	      // find :: energy bin: <Number> to <Number> 
	      if (energyBin(ss,Einit,Eend,&Ea,&Eb))  
		{
		  if (eIndex==0)
		    TS.nTime=tIndex;
		  eIndex++;
		}
	      else if (timeStart(ss))    // find :: <spc> time 
		{
		  Ftime=1;
		  tIndex=0;
		}
	    }
	  else  // In the time data section
	    {
	      if (timeEnd(ss))     // Found "total"
		Ftime=0;
	      else
		{
		  // Need to read the line in the case of first run
		  if (TS.nTime==0)
		    {
		      if (cmdnumberD(ss,&T) &&
			  cmdnumberD(ss,&D))
			tIndex++;
		    }
		}
	    }
	}
    }
  // 
  // Plus 2 since we have a 0 counter and we have missed the last line.
  // 
  TS.nEnergy=eIndex+2;
  if (!TS.nTime && tIndex)  // SMALL TEST if reading first energy bin
    TS.nTime=tIndex;
  // printf("tIndex %d %d %d %d \n",tIndex,eIndex,TS.nEnergy,TS.nTime);

  /* SECOND TIME THROUGH:: */
  rewind(TFile);

  TS.Flux=matrix(TS.nEnergy,TS.nTime);
  TS.EInt=(double*) malloc(TS.nEnergy*sizeof(double));     // Runing integral of energy 
  TS.TimeBin=(double*) malloc(TS.nTime*sizeof(double));
  TS.EnergyBin=(double*) malloc(TS.nEnergy*sizeof(double));
  if (verbose) fprintf(stderr,"NUMBER %d %d \n",TS.nEnergy,TS.nTime);

  Tsum=0.0;
  Ea=0.0;
  Eb=0.0;
  eIndex=-1;
  Ftime=0;
  tIndex=0;
  TS.EInt[0]=0.0;

  // Read file and get time bins
  while(fgets(ss,255,TFile) && Eend>Ea)
    {
      if (notComment(ss))
	{
          if (!Ftime)
	    {
	      if (energyBin(ss,Einit,Eend,&Ea,&Eb))
		{
		  eIndex++;
		  TS.EnergyBin[eIndex]=(Einit>Ea) ? Einit : Ea;
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
		  TS.EInt[eIndex+1]=Tsum;
		  
		  if (verbose) fprintf(stderr,"Energy %g %g \n",EtoLambda(Ea),(TS.EInt[eIndex+1]-TS.EInt[eIndex])*3.744905847e14);
		}
	      else
		{
		  // Need to read the line in the case of first run
		  if (cmdnumberD(ss,&T) &&
		      cmdnumberD(ss,&D))
		    {
		      TS.TimeBin[tIndex]=T/1e8;     // convert Time into second (from shakes)
		      Tsum+=D*Efraction;
		      TS.Flux[eIndex][tIndex]=Tsum;
		      tIndex++;
		    }
		}
	    }
	}
    }
  TS.EnergyBin[eIndex+1]=Eend;
  TS.Total=Tsum;

  // printf("tIndex %d %d %d \n",tIndex,eIndex,TS.nTime);
  // printf("Tsum %g \n",Tsum);
  // fprintf(stderr,"ebin1 ebinN %g %g\n",TS.EnergyBin[0],TS.EnergyBin[TS.nEnergy-1]);

  return;
}

void
getPoint(double* TV,double* EV,double* lim1, double* lim2)
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

  extern Source TS;
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
  Rend=R=TS.Total*R0;
  // This gives Eint[Epnt-1] > R > Eint[Epnt]
  Epnt=binSearch(TS.nEnergy-1,TS.EInt,R);
  Tpnt=binSearch(TS.nTime-1,TS.Flux[Epnt-1],R);

  if(Epnt<0 || R<TS.Flux[Epnt-1][Tpnt-1] || R >TS.Flux[Epnt-1][Tpnt] )
    {
      fprintf(stderr,"outside bin limits Tpnt/Epnt problem  %12.6e %12.6e %12.6e \n",
      TS.Flux[Epnt-1][Tpnt-1],R,TS.Flux[Epnt-1][Tpnt]);
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
      *EV=(Epnt>TS.nEnergy-1) ? TS.EnergyBin[Epnt+1] : TS.EnergyBin[Epnt];
    }

  if (Tpnt==0 || Epnt==0 || Tpnt==TS.nTime)
    {
      fprintf(stderr,"BIG ERROR WITH Tpnt: %d and Epnt: %d\n",Tpnt,Epnt);
      exit(1);
    }
  *TV=TS.TimeBin[Tpnt-1];
  TRange=TS.TimeBin[Tpnt]-TS.TimeBin[Tpnt-1];
  Tspread=TS.Flux[Epnt-1][Tpnt]-TS.Flux[Epnt-1][Tpnt-1];
  R-=TS.Flux[Epnt-1][Tpnt-1];

  R/=Tspread;
  *TV+=TRange*R;


  R1=TS.EInt[Epnt-1]+Espread*rand01();
  iStart=Epnt>lowPoly ? Epnt-lowPoly : 0;                  // max(Epnt-halfPoly,0)
  iEnd=TS.nEnergy>Epnt+highPoly ? Epnt+highPoly : TS.nEnergy-1;  // min(nEnergy-1,Epnt+highPoly

  *EV=polInterp(TS.EInt+iStart,TS.EnergyBin+iStart,1+iEnd-iStart,R1);

  if(*TV < TS.TimeBin[Tpnt-1] || *TV > TS.TimeBin[Tpnt])
    {
      fprintf(stderr,"%d Tpnt %d Tval %g Epnt %d \n",TS.nTime,Tpnt,*TV,Epnt);
      fprintf(stderr,"TBoundary == %12.6e,%g , %12.6e \n\n",TS.TimeBin[Tpnt-1],*TV,TS.TimeBin[Tpnt]);
    }
  

  if(*EV < *lim1 || *EV > *lim2)
    {
      fprintf(stderr,"outside boundaries\n Epnt= %d, Tpnt= %d binlo %g|%g| binhi %g \n",
	      Epnt,Tpnt,TS.EnergyBin[Epnt-1],*EV,TS.EnergyBin[Epnt]);

      fprintf(stderr,"TS == %g %g :: %d %d \n",TS.EInt[Epnt-1],TS.EInt[Epnt],iStart,iEnd);
      fprintf(stderr,"Points (%g) == ",R1);
      for(i=0;i<iEnd-iStart;i++)
	fprintf(stderr,"Time start %g %g",*(TS.EInt+i+iStart),*(TS.EnergyBin+iStart+i));
      fprintf(stderr,"\n");
      exit(1);
    }
  return;
}

int
cmdnumberI(char *mc,int* num,const int len)
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


FILE* 
openFile(char* FileName)
/*
  Tries to open the file:
  (i) In current working directory 
  (ii) In MC_Path directory + ISIS_tables extension
*/
{
  FILE* efile=0;
  int fLen;
  char extFileName[1024];
  char testFileName[2048];
  char mct[2048];
  
  fLen=strlen(FileName);

  if (fLen>512) 
    {
      fprintf(stderr,"Filename excessivley long [EXIT]:\n %s\n",FileName);
      exit(1);
    }
  
  
  strcpy(extFileName,FileName);
  strcpy(extFileName+fLen,".mcstas");
  
  /* Now check for the requested moderator file. In terms of precedence, 
     1) Use MCTABLES location if available
     2) Is a local file available in PWD?
     3) Is there an ISIS_tables in PWD?
     4) Is the file available from the MCSTAS/contrib/ISIS_tables folder?
     - otherwise fail */

  fprintf(stderr,"Searching for %s in multiple locations... -- \n",FileName);

  /* 1) Is MCTABLES set and file located there? */
  if (getenv("MCTABLES"))
    {
      strcpy(mct, getenv("MCTABLES"));
      sprintf(testFileName, "%s%c%s", mct, MC_PATHSEP_C, FileName);
      efile=fopen(testFileName,"r");
      if (efile) {
	fprintf(stderr," - Found in MCTABLES folder %s!\n",mct);
	return efile;
      }
    }

  /* 2) Is the file located in working dir? */
  efile=fopen(FileName,"r");
  if (efile) {
    fprintf(stderr," - Found in current directory!\n");
    return efile;
  }
  
  efile=fopen(extFileName,"r");
  if (efile) return efile; 
  
  /* 3) Is the file in a local 'tablefiles' folder? */
  sprintf(testFileName,"%s%c%s","ISIS_tables",MC_PATHSEP_C,FileName);
  efile=fopen(testFileName,"r");
  if (efile) {
    fprintf(stderr," - Found in local ISIS_tables directory!\n");
    return efile;
  }

  /* 4) Is the file available within the MCSTAS install dir? */
  sprintf(testFileName,"%s%c%s%c%s%c%s",MCSTAS,MC_PATHSEP_C,"contrib",MC_PATHSEP_C,"ISIS_tables",MC_PATHSEP_C,FileName);
  efile=fopen(testFileName,"r");
  if (efile) {
    fprintf(stderr," - Found in MCSTAS system dir: \n %s%c%s%c%s%\n",MCSTAS,MC_PATHSEP_C,"contrib",MC_PATHSEP_C,"ISIS_tables");
    return efile;
  }

  /* If we reach here, the file was not found - raise fatal error */
  fprintf(stderr,"% - Not found! \nPlease check your McStas installation and/or MCTABLES environment variable!\n",FileName);
  exit(1);
}

double strArea()
{
  /*
    Returns the mean Str view of the viewport
    This integrates over each point on the window focus_xw to focus_yh
    View port is symmetric so use only 1/4 of the view
    for the calcuation.
    Control Values TS.XAxis TS.ZAxis focus_xw focus_yh
  */
  
  double A;
  double Vx,Vy;        // view temp points
  double Mx,My;        // moderator x,y
  double D2;           // Distance ^2
  double projArea;     // viewport projection to moderator
  int i,j,aa,bb;       // loop variables
  int NStep;

  NStep=50;
  D2=dist*dist;
  A=0.0;
  fprintf(stderr,"TS axis == %g %g\n",TS.XAxis,TS.ZAxis);
  fprintf(stderr,"AW axis == %g %g\n",focus_xw,focus_yh);
  for(i=0;i<NStep;i++)              // Mod X
    {
      Mx=(i+0.5)*TS.XAxis/(NStep*2.0);
      for(j=0;j<NStep;j++)         // Mod Y
	{
	  My=(j+0.5)*TS.ZAxis/(NStep*2.0);
	  // Position on moderator == (Mx,My)
	  for(aa= -NStep;aa<=NStep;aa++)  //view port
	    for(bb= -NStep;bb<=NStep;bb++)
	      {
		Vx=aa*focus_xw/(2.0*NStep+1.0);
		Vy=bb*focus_yh/(2.0*NStep+1.0);
		A+=1.0/((Mx-Vx)*(Mx-Vx)+(My-Vy)*(My-Vy)+D2);
	      }
	}
    }
  //change to Mx*My 
  A*= (TS.XAxis*TS.ZAxis)/(NStep*NStep*(2.0*NStep+1.0)*(2.0*NStep+1.0));
  // Correct for the area of the viewport. (tables are per cm^2)
  A*=focus_xw*focus_yh*10000;

  //   testing only - Goran              
    // if (!modPosition)
        // {
            // projArea=focus_xw*focus_yh/tally/tally*(tally+dist)*(tally+dist);
	    // A*=TS.XAxis*TS.ZAxis/(TS.XAxis*TS.ZAxis-projArea);
         // }  

  fprintf(stderr,"Viewport == %g %g Moderator size == (%g * %g) m^2 \n",focus_xw,focus_yh,TS.XAxis,TS.ZAxis);
  fprintf(stderr,"Dist == %g (metres) \n",dist);
  fprintf(stderr,"Viewport Solid angle == %g str\n",A/(focus_xw*focus_yh*10000.0));
  fprintf(stderr,"Solid angle used == %g str\n",A);
  return A;
}

#line 9395 "ISIS_MERLIN.c"
#undef beamcurrent
#undef verbose
#undef dist
#undef focus_yh
#undef focus_xw
#undef yh
#undef xw
#undef modPosition
#undef E1
#undef E0
#undef Face
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'TOF_mod' [3]. */
#define mccompcurname  TOF_mod
#define mccompcurtype  TOF_monitor
#define mccompcurindex 3
#define nt mccTOF_mod_nt
#define TOF_N mccTOF_mod_TOF_N
#define TOF_p mccTOF_mod_TOF_p
#define TOF_p2 mccTOF_mod_TOF_p2
#define t_min mccTOF_mod_t_min
#define t_max mccTOF_mod_t_max
#define delta_t mccTOF_mod_delta_t
#define filename mccTOF_mod_filename
#define xmin mccTOF_mod_xmin
#define xmax mccTOF_mod_xmax
#define ymin mccTOF_mod_ymin
#define ymax mccTOF_mod_ymax
#define xwidth mccTOF_mod_xwidth
#define yheight mccTOF_mod_yheight
#define tmin mccTOF_mod_tmin
#define tmax mccTOF_mod_tmax
#define dt mccTOF_mod_dt
#define restore_neutron mccTOF_mod_restore_neutron
#define nowritefile mccTOF_mod_nowritefile
#line 54 "/usr/share/mcstas/2.5/monitors/TOF_monitor.comp"
double TOF_N[nt];
double TOF_p[nt];
double TOF_p2[nt];
double t_min, t_max, delta_t;
#line 9439 "ISIS_MERLIN.c"
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

/* User declarations for component 'shutter_guide' [4]. */
#define mccompcurname  shutter_guide
#define mccompcurtype  Guide_channeled
#define mccompcurindex 4
#define w1c mccshutter_guide_w1c
#define w2c mccshutter_guide_w2c
#define ww mccshutter_guide_ww
#define hh mccshutter_guide_hh
#define whalf mccshutter_guide_whalf
#define hhalf mccshutter_guide_hhalf
#define lwhalf mccshutter_guide_lwhalf
#define lhhalf mccshutter_guide_lhhalf
#define w1 mccshutter_guide_w1
#define h1 mccshutter_guide_h1
#define w2 mccshutter_guide_w2
#define h2 mccshutter_guide_h2
#define l mccshutter_guide_l
#define R0 mccshutter_guide_R0
#define Qc mccshutter_guide_Qc
#define alpha mccshutter_guide_alpha
#define m mccshutter_guide_m
#define nslit mccshutter_guide_nslit
#define d mccshutter_guide_d
#define Qcx mccshutter_guide_Qcx
#define Qcy mccshutter_guide_Qcy
#define alphax mccshutter_guide_alphax
#define alphay mccshutter_guide_alphay
#define W mccshutter_guide_W
#define mx mccshutter_guide_mx
#define my mccshutter_guide_my
#define nu mccshutter_guide_nu
#define phase mccshutter_guide_phase
#line 80 "/usr/share/mcstas/2.5/optics/Guide_channeled.comp"
double w1c;
double w2c;
double ww, hh;
double whalf, hhalf;
double lwhalf, lhhalf;
#line 9501 "ISIS_MERLIN.c"
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

/* User declarations for component 'conv_guide' [5]. */
#define mccompcurname  conv_guide
#define mccompcurtype  Guide_channeled
#define mccompcurindex 5
#define w1c mccconv_guide_w1c
#define w2c mccconv_guide_w2c
#define ww mccconv_guide_ww
#define hh mccconv_guide_hh
#define whalf mccconv_guide_whalf
#define hhalf mccconv_guide_hhalf
#define lwhalf mccconv_guide_lwhalf
#define lhhalf mccconv_guide_lhhalf
#define w1 mccconv_guide_w1
#define h1 mccconv_guide_h1
#define w2 mccconv_guide_w2
#define h2 mccconv_guide_h2
#define l mccconv_guide_l
#define R0 mccconv_guide_R0
#define Qc mccconv_guide_Qc
#define alpha mccconv_guide_alpha
#define m mccconv_guide_m
#define nslit mccconv_guide_nslit
#define d mccconv_guide_d
#define Qcx mccconv_guide_Qcx
#define Qcy mccconv_guide_Qcy
#define alphax mccconv_guide_alphax
#define alphay mccconv_guide_alphay
#define W mccconv_guide_W
#define mx mccconv_guide_mx
#define my mccconv_guide_my
#define nu mccconv_guide_nu
#define phase mccconv_guide_phase
#line 80 "/usr/share/mcstas/2.5/optics/Guide_channeled.comp"
double w1c;
double w2c;
double ww, hh;
double whalf, hhalf;
double lwhalf, lhhalf;
#line 9572 "ISIS_MERLIN.c"
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

/* User declarations for component 'between_guide' [6]. */
#define mccompcurname  between_guide
#define mccompcurtype  Guide_channeled
#define mccompcurindex 6
#define w1c mccbetween_guide_w1c
#define w2c mccbetween_guide_w2c
#define ww mccbetween_guide_ww
#define hh mccbetween_guide_hh
#define whalf mccbetween_guide_whalf
#define hhalf mccbetween_guide_hhalf
#define lwhalf mccbetween_guide_lwhalf
#define lhhalf mccbetween_guide_lhhalf
#define w1 mccbetween_guide_w1
#define h1 mccbetween_guide_h1
#define w2 mccbetween_guide_w2
#define h2 mccbetween_guide_h2
#define l mccbetween_guide_l
#define R0 mccbetween_guide_R0
#define Qc mccbetween_guide_Qc
#define alpha mccbetween_guide_alpha
#define m mccbetween_guide_m
#define nslit mccbetween_guide_nslit
#define d mccbetween_guide_d
#define Qcx mccbetween_guide_Qcx
#define Qcy mccbetween_guide_Qcy
#define alphax mccbetween_guide_alphax
#define alphay mccbetween_guide_alphay
#define W mccbetween_guide_W
#define mx mccbetween_guide_mx
#define my mccbetween_guide_my
#define nu mccbetween_guide_nu
#define phase mccbetween_guide_phase
#line 80 "/usr/share/mcstas/2.5/optics/Guide_channeled.comp"
double w1c;
double w2c;
double ww, hh;
double whalf, hhalf;
double lwhalf, lhhalf;
#line 9643 "ISIS_MERLIN.c"
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

/* User declarations for component 'fermi' [7]. */
#define mccompcurname  fermi
#define mccompcurtype  FermiChopper
#define mccompcurindex 7
#define FCVars mccfermi_FCVars
#define phase mccfermi_phase
#define radius mccfermi_radius
#define nu mccfermi_nu
#define w mccfermi_w
#define nslit mccfermi_nslit
#define R0 mccfermi_R0
#define Qc mccfermi_Qc
#define alpha mccfermi_alpha
#define m mccfermi_m
#define W mccfermi_W
#define length mccfermi_length
#define eff mccfermi_eff
#define zero_time mccfermi_zero_time
#define xwidth mccfermi_xwidth
#define verbose mccfermi_verbose
#define yheight mccfermi_yheight
#define curvature mccfermi_curvature
#define delay mccfermi_delay
#line 278 "/usr/share/mcstas/2.5/optics/FermiChopper.comp"
        struct FermiChopper_struct FCVars;
#line 9701 "ISIS_MERLIN.c"
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

/* User declarations for component 'TOF_fermi' [8]. */
#define mccompcurname  TOF_fermi
#define mccompcurtype  TOF_monitor
#define mccompcurindex 8
#define nt mccTOF_fermi_nt
#define TOF_N mccTOF_fermi_TOF_N
#define TOF_p mccTOF_fermi_TOF_p
#define TOF_p2 mccTOF_fermi_TOF_p2
#define t_min mccTOF_fermi_t_min
#define t_max mccTOF_fermi_t_max
#define delta_t mccTOF_fermi_delta_t
#define filename mccTOF_fermi_filename
#define xmin mccTOF_fermi_xmin
#define xmax mccTOF_fermi_xmax
#define ymin mccTOF_fermi_ymin
#define ymax mccTOF_fermi_ymax
#define xwidth mccTOF_fermi_xwidth
#define yheight mccTOF_fermi_yheight
#define tmin mccTOF_fermi_tmin
#define tmax mccTOF_fermi_tmax
#define dt mccTOF_fermi_dt
#define restore_neutron mccTOF_fermi_restore_neutron
#define nowritefile mccTOF_fermi_nowritefile
#line 54 "/usr/share/mcstas/2.5/monitors/TOF_monitor.comp"
double TOF_N[nt];
double TOF_p[nt];
double TOF_p2[nt];
double t_min, t_max, delta_t;
#line 9753 "ISIS_MERLIN.c"
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

/* User declarations for component 'final_Guide' [9]. */
#define mccompcurname  final_Guide
#define mccompcurtype  Guide_channeled
#define mccompcurindex 9
#define w1c mccfinal_Guide_w1c
#define w2c mccfinal_Guide_w2c
#define ww mccfinal_Guide_ww
#define hh mccfinal_Guide_hh
#define whalf mccfinal_Guide_whalf
#define hhalf mccfinal_Guide_hhalf
#define lwhalf mccfinal_Guide_lwhalf
#define lhhalf mccfinal_Guide_lhhalf
#define w1 mccfinal_Guide_w1
#define h1 mccfinal_Guide_h1
#define w2 mccfinal_Guide_w2
#define h2 mccfinal_Guide_h2
#define l mccfinal_Guide_l
#define R0 mccfinal_Guide_R0
#define Qc mccfinal_Guide_Qc
#define alpha mccfinal_Guide_alpha
#define m mccfinal_Guide_m
#define nslit mccfinal_Guide_nslit
#define d mccfinal_Guide_d
#define Qcx mccfinal_Guide_Qcx
#define Qcy mccfinal_Guide_Qcy
#define alphax mccfinal_Guide_alphax
#define alphay mccfinal_Guide_alphay
#define W mccfinal_Guide_W
#define mx mccfinal_Guide_mx
#define my mccfinal_Guide_my
#define nu mccfinal_Guide_nu
#define phase mccfinal_Guide_phase
#line 80 "/usr/share/mcstas/2.5/optics/Guide_channeled.comp"
double w1c;
double w2c;
double ww, hh;
double whalf, hhalf;
double lwhalf, lhhalf;
#line 9815 "ISIS_MERLIN.c"
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

/* User declarations for component 'horz_div' [10]. */
#define mccompcurname  horz_div
#define mccompcurtype  Hdiv_monitor
#define mccompcurindex 10
#define nh mcchorz_div_nh
#define Div_N mcchorz_div_Div_N
#define Div_p mcchorz_div_Div_p
#define Div_p2 mcchorz_div_Div_p2
#define filename mcchorz_div_filename
#define xmin mcchorz_div_xmin
#define xmax mcchorz_div_xmax
#define ymin mcchorz_div_ymin
#define ymax mcchorz_div_ymax
#define xwidth mcchorz_div_xwidth
#define yheight mcchorz_div_yheight
#define h_maxdiv mcchorz_div_h_maxdiv
#define restore_neutron mcchorz_div_restore_neutron
#define nowritefile mcchorz_div_nowritefile
#line 57 "/usr/share/mcstas/2.5/monitors/Hdiv_monitor.comp"
double Div_N[nh];
double Div_p[nh];
double Div_p2[nh];
#line 9870 "ISIS_MERLIN.c"
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

/* User declarations for component 'PSD' [11]. */
#define mccompcurname  PSD
#define mccompcurtype  PSD_monitor
#define mccompcurindex 11
#define nx mccPSD_nx
#define ny mccPSD_ny
#define PSD_N mccPSD_PSD_N
#define PSD_p mccPSD_PSD_p
#define PSD_p2 mccPSD_PSD_p2
#define filename mccPSD_filename
#define xmin mccPSD_xmin
#define xmax mccPSD_xmax
#define ymin mccPSD_ymin
#define ymax mccPSD_ymax
#define xwidth mccPSD_xwidth
#define yheight mccPSD_yheight
#define restore_neutron mccPSD_restore_neutron
#define nowritefile mccPSD_nowritefile
#line 56 "/usr/share/mcstas/2.5/monitors/PSD_monitor.comp"
double PSD_N[nx][ny];
double PSD_p[nx][ny];
double PSD_p2[nx][ny];
#line 9911 "ISIS_MERLIN.c"
#undef nowritefile
#undef restore_neutron
#undef yheight
#undef xwidth
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef filename
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef ny
#undef nx
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'det_lmon' [12]. */
#define mccompcurname  det_lmon
#define mccompcurtype  L_monitor
#define mccompcurindex 12
#define nL mccdet_lmon_nL
#define L_N mccdet_lmon_L_N
#define L_p mccdet_lmon_L_p
#define L_p2 mccdet_lmon_L_p2
#define filename mccdet_lmon_filename
#define xmin mccdet_lmon_xmin
#define xmax mccdet_lmon_xmax
#define ymin mccdet_lmon_ymin
#define ymax mccdet_lmon_ymax
#define xwidth mccdet_lmon_xwidth
#define yheight mccdet_lmon_yheight
#define Lmin mccdet_lmon_Lmin
#define Lmax mccdet_lmon_Lmax
#define restore_neutron mccdet_lmon_restore_neutron
#define nowritefile mccdet_lmon_nowritefile
#line 57 "/usr/share/mcstas/2.5/monitors/L_monitor.comp"
double L_N[nL];
double L_p[nL], L_p2[nL];
#line 9952 "ISIS_MERLIN.c"
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

/* User declarations for component 'sample' [13]. */
#define mccompcurname  sample
#define mccompcurtype  Incoherent
#define mccompcurindex 13
#define VarsInc mccsample_VarsInc
#define offdata mccsample_offdata
#define geometry mccsample_geometry
#define radius mccsample_radius
#define xwidth mccsample_xwidth
#define yheight mccsample_yheight
#define zdepth mccsample_zdepth
#define thickness mccsample_thickness
#define target_x mccsample_target_x
#define target_y mccsample_target_y
#define target_z mccsample_target_z
#define focus_r mccsample_focus_r
#define focus_xw mccsample_focus_xw
#define focus_yh mccsample_focus_yh
#define focus_aw mccsample_focus_aw
#define focus_ah mccsample_focus_ah
#define target_index mccsample_target_index
#define pack mccsample_pack
#define p_interact mccsample_p_interact
#define f_QE mccsample_f_QE
#define gamma mccsample_gamma
#define Etrans mccsample_Etrans
#define deltaE mccsample_deltaE
#define sigma_abs mccsample_sigma_abs
#define sigma_inc mccsample_sigma_inc
#define Vc mccsample_Vc
#define concentric mccsample_concentric
#define order mccsample_order
#line 148 "/usr/share/mcstas/2.5/samples/Incoherent.comp"

  struct StructVarsInc VarsInc;
  off_struct offdata;

#line 10009 "ISIS_MERLIN.c"
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

/* User declarations for component 'detectorarm' [14]. */
#define mccompcurname  detectorarm
#define mccompcurtype  Arm
#define mccompcurindex 14
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'detector' [15]. */
#define mccompcurname  detector
#define mccompcurtype  TOF_monitor
#define mccompcurindex 15
#define nt mccdetector_nt
#define TOF_N mccdetector_TOF_N
#define TOF_p mccdetector_TOF_p
#define TOF_p2 mccdetector_TOF_p2
#define t_min mccdetector_t_min
#define t_max mccdetector_t_max
#define delta_t mccdetector_delta_t
#define filename mccdetector_filename
#define xmin mccdetector_xmin
#define xmax mccdetector_xmax
#define ymin mccdetector_ymin
#define ymax mccdetector_ymax
#define xwidth mccdetector_xwidth
#define yheight mccdetector_yheight
#define tmin mccdetector_tmin
#define tmax mccdetector_tmax
#define dt mccdetector_dt
#define restore_neutron mccdetector_restore_neutron
#define nowritefile mccdetector_nowritefile
#line 54 "/usr/share/mcstas/2.5/monitors/TOF_monitor.comp"
double TOF_N[nt];
double TOF_p[nt];
double TOF_p2[nt];
double t_min, t_max, delta_t;
#line 10078 "ISIS_MERLIN.c"
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

Coords mcposaOrigin, mcposrOrigin;
Rotation mcrotaOrigin, mcrotrOrigin;
Coords mcposaisis_mod, mcposrisis_mod;
Rotation mcrotaisis_mod, mcrotrisis_mod;
Coords mcposaTOF_mod, mcposrTOF_mod;
Rotation mcrotaTOF_mod, mcrotrTOF_mod;
Coords mcposashutter_guide, mcposrshutter_guide;
Rotation mcrotashutter_guide, mcrotrshutter_guide;
Coords mcposaconv_guide, mcposrconv_guide;
Rotation mcrotaconv_guide, mcrotrconv_guide;
Coords mcposabetween_guide, mcposrbetween_guide;
Rotation mcrotabetween_guide, mcrotrbetween_guide;
Coords mcposafermi, mcposrfermi;
Rotation mcrotafermi, mcrotrfermi;
Coords mcposaTOF_fermi, mcposrTOF_fermi;
Rotation mcrotaTOF_fermi, mcrotrTOF_fermi;
Coords mcposafinal_Guide, mcposrfinal_Guide;
Rotation mcrotafinal_Guide, mcrotrfinal_Guide;
Coords mcposahorz_div, mcposrhorz_div;
Rotation mcrotahorz_div, mcrotrhorz_div;
Coords mcposaPSD, mcposrPSD;
Rotation mcrotaPSD, mcrotrPSD;
Coords mcposadet_lmon, mcposrdet_lmon;
Rotation mcrotadet_lmon, mcrotrdet_lmon;
Coords mcposasample, mcposrsample;
Rotation mcrotasample, mcrotrsample;
Coords mcposadetectorarm, mcposrdetectorarm;
Rotation mcrotadetectorarm, mcrotrdetectorarm;
Coords mcposadetector, mcposrdetector;
Rotation mcrotadetector, mcrotrdetector;

MCNUM mcnx, mcny, mcnz, mcnvx, mcnvy, mcnvz, mcnt, mcnsx, mcnsy, mcnsz, mcnp;

/* end declare */

void mcinit(void) {
#define mccompcurname  ISIS_MERLIN
#define mccompcurtype  INSTRUMENT
#define mccompcurindex 0
#define mcposaISIS_MERLIN coords_set(0,0,0)
#define E_min mcipE_min
#define E_max mcipE_max
#define m mcipm
#define E_foc mcipE_foc
#define nu_hz mcipnu_hz
#define type mciptype
#line 56 "ISIS_MERLIN.instr"
{

if (type==1)
  {slit_curv=1/0.82;
   num_slits=55;
   width=0.00114;
   len=0.099;
   fprintf(stderr,"MERLIN b chopper selected");
  }
if (type==2)
  {slit_curv=1/1.3;
   num_slits=28;
   width=0.00228;
   len=0.099;
   fprintf(stderr,"MERLIN Sloppy selected selected");
  }
if (type==3)
  {slit_curv=0;
   num_slits=350;
   width=0.0002;
   len=0.01;
   fprintf(stderr,"Merlin Gd selected selected");
  }

/* calculate opening time of the fermi chopper at 10m from moderator*/
phase_time = (2.28e-3*10.0/sqrt(E_foc));
printf("\nphase_time: %f\n\n",phase_time);
t_det=1.43*phase_time;
t_fer=1.02*phase_time;
lam = sqrt(81.81/E_foc);
}
#line 10180 "ISIS_MERLIN.c"
#undef type
#undef nu_hz
#undef E_foc
#undef m
#undef E_max
#undef E_min
#undef mcposaISIS_MERLIN
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
#line 10210 "ISIS_MERLIN.c"
  rot_copy(mcrotrOrigin, mcrotaOrigin);
  mcposaOrigin = coords_set(
#line 92 "ISIS_MERLIN.instr"
    0,
#line 92 "ISIS_MERLIN.instr"
    0,
#line 92 "ISIS_MERLIN.instr"
    0);
#line 10219 "ISIS_MERLIN.c"
  mctc1 = coords_neg(mcposaOrigin);
  mcposrOrigin = rot_apply(mcrotaOrigin, mctc1);
  mcDEBUG_COMPONENT("Origin", mcposaOrigin, mcrotaOrigin)
  mccomp_posa[1] = mcposaOrigin;
  mccomp_posr[1] = mcposrOrigin;
  mcNCounter[1]  = mcPCounter[1] = mcP2Counter[1] = 0;
  mcAbsorbProp[1]= 0;
    /* Component isis_mod. */
  /* Setting parameters for component isis_mod. */
  SIG_MESSAGE("isis_mod (Init:SetPar)");
#line 96 "ISIS_MERLIN.instr"
  if("TS1_S04_Merlin.mcstas") strncpy(mccisis_mod_Face, "TS1_S04_Merlin.mcstas" ? "TS1_S04_Merlin.mcstas" : "", 16384); else mccisis_mod_Face[0]='\0';
#line 96 "ISIS_MERLIN.instr"
  mccisis_mod_E0 = mcipE_min;
#line 96 "ISIS_MERLIN.instr"
  mccisis_mod_E1 = mcipE_max;
#line 97 "ISIS_MERLIN.instr"
  mccisis_mod_modPosition = 0;
#line 97 "ISIS_MERLIN.instr"
  mccisis_mod_xw = 0.12;
#line 97 "ISIS_MERLIN.instr"
  mccisis_mod_yh = 0.12;
#line 97 "ISIS_MERLIN.instr"
  mccisis_mod_focus_xw = 0.094;
#line 97 "ISIS_MERLIN.instr"
  mccisis_mod_focus_yh = 0.094;
#line 97 "ISIS_MERLIN.instr"
  mccisis_mod_dist = 1.7;
#line 76 "ISIS_MERLIN.instr"
  mccisis_mod_verbose = 0;
#line 76 "ISIS_MERLIN.instr"
  mccisis_mod_beamcurrent = 1;
#line 10252 "ISIS_MERLIN.c"

  SIG_MESSAGE("isis_mod (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 10259 "ISIS_MERLIN.c"
  rot_mul(mctr1, mcrotaOrigin, mcrotaisis_mod);
  rot_transpose(mcrotaOrigin, mctr1);
  rot_mul(mcrotaisis_mod, mctr1, mcrotrisis_mod);
  mctc1 = coords_set(
#line 98 "ISIS_MERLIN.instr"
    0,
#line 98 "ISIS_MERLIN.instr"
    0,
#line 98 "ISIS_MERLIN.instr"
    0);
#line 10270 "ISIS_MERLIN.c"
  rot_transpose(mcrotaOrigin, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaisis_mod = coords_add(mcposaOrigin, mctc2);
  mctc1 = coords_sub(mcposaOrigin, mcposaisis_mod);
  mcposrisis_mod = rot_apply(mcrotaisis_mod, mctc1);
  mcDEBUG_COMPONENT("isis_mod", mcposaisis_mod, mcrotaisis_mod)
  mccomp_posa[2] = mcposaisis_mod;
  mccomp_posr[2] = mcposrisis_mod;
  mcNCounter[2]  = mcPCounter[2] = mcP2Counter[2] = 0;
  mcAbsorbProp[2]= 0;
    /* Component TOF_mod. */
  /* Setting parameters for component TOF_mod. */
  SIG_MESSAGE("TOF_mod (Init:SetPar)");
#line 102 "ISIS_MERLIN.instr"
  if("mod.m") strncpy(mccTOF_mod_filename, "mod.m" ? "mod.m" : "", 16384); else mccTOF_mod_filename[0]='\0';
#line 102 "ISIS_MERLIN.instr"
  mccTOF_mod_xmin = -0.025;
#line 103 "ISIS_MERLIN.instr"
  mccTOF_mod_xmax = 0.025;
#line 103 "ISIS_MERLIN.instr"
  mccTOF_mod_ymin = -0.025;
#line 103 "ISIS_MERLIN.instr"
  mccTOF_mod_ymax = 0.025;
#line 48 "ISIS_MERLIN.instr"
  mccTOF_mod_xwidth = 0;
#line 48 "ISIS_MERLIN.instr"
  mccTOF_mod_yheight = 0;
#line 103 "ISIS_MERLIN.instr"
  mccTOF_mod_tmin = 0;
#line 103 "ISIS_MERLIN.instr"
  mccTOF_mod_tmax = 150;
#line 48 "ISIS_MERLIN.instr"
  mccTOF_mod_dt = 1.0;
#line 48 "ISIS_MERLIN.instr"
  mccTOF_mod_restore_neutron = 0;
#line 48 "ISIS_MERLIN.instr"
  mccTOF_mod_nowritefile = 0;
#line 10308 "ISIS_MERLIN.c"

  SIG_MESSAGE("TOF_mod (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 10315 "ISIS_MERLIN.c"
  rot_mul(mctr1, mcrotaOrigin, mcrotaTOF_mod);
  rot_transpose(mcrotaisis_mod, mctr1);
  rot_mul(mcrotaTOF_mod, mctr1, mcrotrTOF_mod);
  mctc1 = coords_set(
#line 104 "ISIS_MERLIN.instr"
    0,
#line 104 "ISIS_MERLIN.instr"
    0,
#line 104 "ISIS_MERLIN.instr"
    0.001);
#line 10326 "ISIS_MERLIN.c"
  rot_transpose(mcrotaOrigin, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaTOF_mod = coords_add(mcposaOrigin, mctc2);
  mctc1 = coords_sub(mcposaisis_mod, mcposaTOF_mod);
  mcposrTOF_mod = rot_apply(mcrotaTOF_mod, mctc1);
  mcDEBUG_COMPONENT("TOF_mod", mcposaTOF_mod, mcrotaTOF_mod)
  mccomp_posa[3] = mcposaTOF_mod;
  mccomp_posr[3] = mcposrTOF_mod;
  mcNCounter[3]  = mcPCounter[3] = mcP2Counter[3] = 0;
  mcAbsorbProp[3]= 0;
    /* Component shutter_guide. */
  /* Setting parameters for component shutter_guide. */
  SIG_MESSAGE("shutter_guide (Init:SetPar)");
#line 108 "ISIS_MERLIN.instr"
  mccshutter_guide_w1 = 0.094;
#line 108 "ISIS_MERLIN.instr"
  mccshutter_guide_h1 = 0.094;
#line 108 "ISIS_MERLIN.instr"
  mccshutter_guide_w2 = 0.094;
#line 108 "ISIS_MERLIN.instr"
  mccshutter_guide_h2 = 0.094;
#line 108 "ISIS_MERLIN.instr"
  mccshutter_guide_l = 2.0;
#line 71 "ISIS_MERLIN.instr"
  mccshutter_guide_R0 = 0.995;
#line 71 "ISIS_MERLIN.instr"
  mccshutter_guide_Qc = 0;
#line 71 "ISIS_MERLIN.instr"
  mccshutter_guide_alpha = 0;
#line 71 "ISIS_MERLIN.instr"
  mccshutter_guide_m = 0;
#line 71 "ISIS_MERLIN.instr"
  mccshutter_guide_nslit = 1;
#line 71 "ISIS_MERLIN.instr"
  mccshutter_guide_d = 0.0005;
#line 72 "ISIS_MERLIN.instr"
  mccshutter_guide_Qcx = 0.0218;
#line 72 "ISIS_MERLIN.instr"
  mccshutter_guide_Qcy = 0.0218;
#line 108 "ISIS_MERLIN.instr"
  mccshutter_guide_alphax = 4.38;
#line 108 "ISIS_MERLIN.instr"
  mccshutter_guide_alphay = 4.38;
#line 109 "ISIS_MERLIN.instr"
  mccshutter_guide_W = 3e-3;
#line 109 "ISIS_MERLIN.instr"
  mccshutter_guide_mx = mcipm;
#line 109 "ISIS_MERLIN.instr"
  mccshutter_guide_my = mcipm;
#line 72 "ISIS_MERLIN.instr"
  mccshutter_guide_nu = 0;
#line 72 "ISIS_MERLIN.instr"
  mccshutter_guide_phase = 0;
#line 10380 "ISIS_MERLIN.c"

  SIG_MESSAGE("shutter_guide (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 10387 "ISIS_MERLIN.c"
  rot_mul(mctr1, mcrotaOrigin, mcrotashutter_guide);
  rot_transpose(mcrotaTOF_mod, mctr1);
  rot_mul(mcrotashutter_guide, mctr1, mcrotrshutter_guide);
  mctc1 = coords_set(
#line 110 "ISIS_MERLIN.instr"
    0,
#line 110 "ISIS_MERLIN.instr"
    0,
#line 110 "ISIS_MERLIN.instr"
    1.7);
#line 10398 "ISIS_MERLIN.c"
  rot_transpose(mcrotaOrigin, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposashutter_guide = coords_add(mcposaOrigin, mctc2);
  mctc1 = coords_sub(mcposaTOF_mod, mcposashutter_guide);
  mcposrshutter_guide = rot_apply(mcrotashutter_guide, mctc1);
  mcDEBUG_COMPONENT("shutter_guide", mcposashutter_guide, mcrotashutter_guide)
  mccomp_posa[4] = mcposashutter_guide;
  mccomp_posr[4] = mcposrshutter_guide;
  mcNCounter[4]  = mcPCounter[4] = mcP2Counter[4] = 0;
  mcAbsorbProp[4]= 0;
    /* Component conv_guide. */
  /* Setting parameters for component conv_guide. */
  SIG_MESSAGE("conv_guide (Init:SetPar)");
#line 114 "ISIS_MERLIN.instr"
  mccconv_guide_w1 = 0.094;
#line 114 "ISIS_MERLIN.instr"
  mccconv_guide_h1 = 0.094;
#line 114 "ISIS_MERLIN.instr"
  mccconv_guide_w2 = 0.067;
#line 114 "ISIS_MERLIN.instr"
  mccconv_guide_h2 = 0.067;
#line 114 "ISIS_MERLIN.instr"
  mccconv_guide_l = 4.760;
#line 71 "ISIS_MERLIN.instr"
  mccconv_guide_R0 = 0.995;
#line 71 "ISIS_MERLIN.instr"
  mccconv_guide_Qc = 0;
#line 71 "ISIS_MERLIN.instr"
  mccconv_guide_alpha = 0;
#line 71 "ISIS_MERLIN.instr"
  mccconv_guide_m = 0;
#line 71 "ISIS_MERLIN.instr"
  mccconv_guide_nslit = 1;
#line 71 "ISIS_MERLIN.instr"
  mccconv_guide_d = 0.0005;
#line 72 "ISIS_MERLIN.instr"
  mccconv_guide_Qcx = 0.0218;
#line 72 "ISIS_MERLIN.instr"
  mccconv_guide_Qcy = 0.0218;
#line 114 "ISIS_MERLIN.instr"
  mccconv_guide_alphax = 4.38;
#line 114 "ISIS_MERLIN.instr"
  mccconv_guide_alphay = 4.38;
#line 115 "ISIS_MERLIN.instr"
  mccconv_guide_W = 3e-3;
#line 115 "ISIS_MERLIN.instr"
  mccconv_guide_mx = mcipm;
#line 115 "ISIS_MERLIN.instr"
  mccconv_guide_my = mcipm;
#line 72 "ISIS_MERLIN.instr"
  mccconv_guide_nu = 0;
#line 72 "ISIS_MERLIN.instr"
  mccconv_guide_phase = 0;
#line 10452 "ISIS_MERLIN.c"

  SIG_MESSAGE("conv_guide (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 10459 "ISIS_MERLIN.c"
  rot_mul(mctr1, mcrotaOrigin, mcrotaconv_guide);
  rot_transpose(mcrotashutter_guide, mctr1);
  rot_mul(mcrotaconv_guide, mctr1, mcrotrconv_guide);
  mctc1 = coords_set(
#line 116 "ISIS_MERLIN.instr"
    0,
#line 116 "ISIS_MERLIN.instr"
    0,
#line 116 "ISIS_MERLIN.instr"
    3.704);
#line 10470 "ISIS_MERLIN.c"
  rot_transpose(mcrotaOrigin, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaconv_guide = coords_add(mcposaOrigin, mctc2);
  mctc1 = coords_sub(mcposashutter_guide, mcposaconv_guide);
  mcposrconv_guide = rot_apply(mcrotaconv_guide, mctc1);
  mcDEBUG_COMPONENT("conv_guide", mcposaconv_guide, mcrotaconv_guide)
  mccomp_posa[5] = mcposaconv_guide;
  mccomp_posr[5] = mcposrconv_guide;
  mcNCounter[5]  = mcPCounter[5] = mcP2Counter[5] = 0;
  mcAbsorbProp[5]= 0;
    /* Component between_guide. */
  /* Setting parameters for component between_guide. */
  SIG_MESSAGE("between_guide (Init:SetPar)");
#line 120 "ISIS_MERLIN.instr"
  mccbetween_guide_w1 = 0.0629;
#line 120 "ISIS_MERLIN.instr"
  mccbetween_guide_h1 = 0.0629;
#line 120 "ISIS_MERLIN.instr"
  mccbetween_guide_w2 = 0.0594;
#line 120 "ISIS_MERLIN.instr"
  mccbetween_guide_h2 = 0.0594;
#line 120 "ISIS_MERLIN.instr"
  mccbetween_guide_l = 0.640;
#line 71 "ISIS_MERLIN.instr"
  mccbetween_guide_R0 = 0.995;
#line 71 "ISIS_MERLIN.instr"
  mccbetween_guide_Qc = 0;
#line 71 "ISIS_MERLIN.instr"
  mccbetween_guide_alpha = 0;
#line 71 "ISIS_MERLIN.instr"
  mccbetween_guide_m = 0;
#line 71 "ISIS_MERLIN.instr"
  mccbetween_guide_nslit = 1;
#line 71 "ISIS_MERLIN.instr"
  mccbetween_guide_d = 0.0005;
#line 72 "ISIS_MERLIN.instr"
  mccbetween_guide_Qcx = 0.0218;
#line 72 "ISIS_MERLIN.instr"
  mccbetween_guide_Qcy = 0.0218;
#line 120 "ISIS_MERLIN.instr"
  mccbetween_guide_alphax = 4.38;
#line 120 "ISIS_MERLIN.instr"
  mccbetween_guide_alphay = 4.38;
#line 121 "ISIS_MERLIN.instr"
  mccbetween_guide_W = 3e-3;
#line 121 "ISIS_MERLIN.instr"
  mccbetween_guide_mx = mcipm;
#line 121 "ISIS_MERLIN.instr"
  mccbetween_guide_my = mcipm;
#line 72 "ISIS_MERLIN.instr"
  mccbetween_guide_nu = 0;
#line 72 "ISIS_MERLIN.instr"
  mccbetween_guide_phase = 0;
#line 10524 "ISIS_MERLIN.c"

  SIG_MESSAGE("between_guide (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 10531 "ISIS_MERLIN.c"
  rot_mul(mctr1, mcrotaOrigin, mcrotabetween_guide);
  rot_transpose(mcrotaconv_guide, mctr1);
  rot_mul(mcrotabetween_guide, mctr1, mcrotrbetween_guide);
  mctc1 = coords_set(
#line 122 "ISIS_MERLIN.instr"
    0,
#line 122 "ISIS_MERLIN.instr"
    0,
#line 122 "ISIS_MERLIN.instr"
    9.182);
#line 10542 "ISIS_MERLIN.c"
  rot_transpose(mcrotaOrigin, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposabetween_guide = coords_add(mcposaOrigin, mctc2);
  mctc1 = coords_sub(mcposaconv_guide, mcposabetween_guide);
  mcposrbetween_guide = rot_apply(mcrotabetween_guide, mctc1);
  mcDEBUG_COMPONENT("between_guide", mcposabetween_guide, mcrotabetween_guide)
  mccomp_posa[6] = mcposabetween_guide;
  mccomp_posr[6] = mcposrbetween_guide;
  mcNCounter[6]  = mcPCounter[6] = mcP2Counter[6] = 0;
  mcAbsorbProp[6]= 0;
    /* Component fermi. */
  /* Setting parameters for component fermi. */
  SIG_MESSAGE("fermi (Init:SetPar)");
#line 88 "ISIS_MERLIN.instr"
  mccfermi_phase = 0;
#line 125 "ISIS_MERLIN.instr"
  mccfermi_radius = 0.055;
#line 125 "ISIS_MERLIN.instr"
  mccfermi_nu = - mcipnu_hz;
#line 126 "ISIS_MERLIN.instr"
  mccfermi_w = width;
#line 126 "ISIS_MERLIN.instr"
  mccfermi_nslit = num_slits;
#line 126 "ISIS_MERLIN.instr"
  mccfermi_R0 = 0.0;
#line 127 "ISIS_MERLIN.instr"
  mccfermi_Qc = 0.02176;
#line 127 "ISIS_MERLIN.instr"
  mccfermi_alpha = 2.33;
#line 127 "ISIS_MERLIN.instr"
  mccfermi_m = 0;
#line 90 "ISIS_MERLIN.instr"
  mccfermi_W = 2e-3;
#line 127 "ISIS_MERLIN.instr"
  mccfermi_length = len;
#line 127 "ISIS_MERLIN.instr"
  mccfermi_eff = 0.95;
#line 127 "ISIS_MERLIN.instr"
  mccfermi_zero_time = 0;
#line 91 "ISIS_MERLIN.instr"
  mccfermi_xwidth = 0;
#line 91 "ISIS_MERLIN.instr"
  mccfermi_verbose = 0;
#line 126 "ISIS_MERLIN.instr"
  mccfermi_yheight = 0.07;
#line 127 "ISIS_MERLIN.instr"
  mccfermi_curvature = - slit_curv;
#line 125 "ISIS_MERLIN.instr"
  mccfermi_delay = phase_time;
#line 10592 "ISIS_MERLIN.c"

  SIG_MESSAGE("fermi (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 10599 "ISIS_MERLIN.c"
  rot_mul(mctr1, mcrotaOrigin, mcrotafermi);
  rot_transpose(mcrotabetween_guide, mctr1);
  rot_mul(mcrotafermi, mctr1, mcrotrfermi);
  mctc1 = coords_set(
#line 128 "ISIS_MERLIN.instr"
    0,
#line 128 "ISIS_MERLIN.instr"
    0,
#line 128 "ISIS_MERLIN.instr"
    10);
#line 10610 "ISIS_MERLIN.c"
  rot_transpose(mcrotaOrigin, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposafermi = coords_add(mcposaOrigin, mctc2);
  mctc1 = coords_sub(mcposabetween_guide, mcposafermi);
  mcposrfermi = rot_apply(mcrotafermi, mctc1);
  mcDEBUG_COMPONENT("fermi", mcposafermi, mcrotafermi)
  mccomp_posa[7] = mcposafermi;
  mccomp_posr[7] = mcposrfermi;
  mcNCounter[7]  = mcPCounter[7] = mcP2Counter[7] = 0;
  mcAbsorbProp[7]= 0;
    /* Component TOF_fermi. */
  /* Setting parameters for component TOF_fermi. */
  SIG_MESSAGE("TOF_fermi (Init:SetPar)");
#line 132 "ISIS_MERLIN.instr"
  if("fermi.m") strncpy(mccTOF_fermi_filename, "fermi.m" ? "fermi.m" : "", 16384); else mccTOF_fermi_filename[0]='\0';
#line 132 "ISIS_MERLIN.instr"
  mccTOF_fermi_xmin = -0.025;
#line 133 "ISIS_MERLIN.instr"
  mccTOF_fermi_xmax = 0.025;
#line 133 "ISIS_MERLIN.instr"
  mccTOF_fermi_ymin = -0.025;
#line 133 "ISIS_MERLIN.instr"
  mccTOF_fermi_ymax = 0.025;
#line 48 "ISIS_MERLIN.instr"
  mccTOF_fermi_xwidth = 0;
#line 48 "ISIS_MERLIN.instr"
  mccTOF_fermi_yheight = 0;
#line 133 "ISIS_MERLIN.instr"
  mccTOF_fermi_tmin = ( 1e6 * t_fer ) -12;
#line 133 "ISIS_MERLIN.instr"
  mccTOF_fermi_tmax = ( 1e6 * t_fer ) + 22;
#line 48 "ISIS_MERLIN.instr"
  mccTOF_fermi_dt = 1.0;
#line 48 "ISIS_MERLIN.instr"
  mccTOF_fermi_restore_neutron = 0;
#line 48 "ISIS_MERLIN.instr"
  mccTOF_fermi_nowritefile = 0;
#line 10648 "ISIS_MERLIN.c"

  SIG_MESSAGE("TOF_fermi (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 10655 "ISIS_MERLIN.c"
  rot_mul(mctr1, mcrotaOrigin, mcrotaTOF_fermi);
  rot_transpose(mcrotafermi, mctr1);
  rot_mul(mcrotaTOF_fermi, mctr1, mcrotrTOF_fermi);
  mctc1 = coords_set(
#line 134 "ISIS_MERLIN.instr"
    0,
#line 134 "ISIS_MERLIN.instr"
    0,
#line 134 "ISIS_MERLIN.instr"
    10.2);
#line 10666 "ISIS_MERLIN.c"
  rot_transpose(mcrotaOrigin, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaTOF_fermi = coords_add(mcposaOrigin, mctc2);
  mctc1 = coords_sub(mcposafermi, mcposaTOF_fermi);
  mcposrTOF_fermi = rot_apply(mcrotaTOF_fermi, mctc1);
  mcDEBUG_COMPONENT("TOF_fermi", mcposaTOF_fermi, mcrotaTOF_fermi)
  mccomp_posa[8] = mcposaTOF_fermi;
  mccomp_posr[8] = mcposrTOF_fermi;
  mcNCounter[8]  = mcPCounter[8] = mcP2Counter[8] = 0;
  mcAbsorbProp[8]= 0;
    /* Component final_Guide. */
  /* Setting parameters for component final_Guide. */
  SIG_MESSAGE("final_Guide (Init:SetPar)");
#line 138 "ISIS_MERLIN.instr"
  mccfinal_Guide_w1 = 0.0568;
#line 138 "ISIS_MERLIN.instr"
  mccfinal_Guide_h1 = 0.0568;
#line 138 "ISIS_MERLIN.instr"
  mccfinal_Guide_w2 = 0.0506;
#line 138 "ISIS_MERLIN.instr"
  mccfinal_Guide_h2 = 0.0506;
#line 138 "ISIS_MERLIN.instr"
  mccfinal_Guide_l = 1.10;
#line 71 "ISIS_MERLIN.instr"
  mccfinal_Guide_R0 = 0.995;
#line 71 "ISIS_MERLIN.instr"
  mccfinal_Guide_Qc = 0;
#line 71 "ISIS_MERLIN.instr"
  mccfinal_Guide_alpha = 0;
#line 71 "ISIS_MERLIN.instr"
  mccfinal_Guide_m = 0;
#line 71 "ISIS_MERLIN.instr"
  mccfinal_Guide_nslit = 1;
#line 71 "ISIS_MERLIN.instr"
  mccfinal_Guide_d = 0.0005;
#line 72 "ISIS_MERLIN.instr"
  mccfinal_Guide_Qcx = 0.0218;
#line 72 "ISIS_MERLIN.instr"
  mccfinal_Guide_Qcy = 0.0218;
#line 138 "ISIS_MERLIN.instr"
  mccfinal_Guide_alphax = 4.38;
#line 138 "ISIS_MERLIN.instr"
  mccfinal_Guide_alphay = 4.38;
#line 139 "ISIS_MERLIN.instr"
  mccfinal_Guide_W = 3e-3;
#line 139 "ISIS_MERLIN.instr"
  mccfinal_Guide_mx = mcipm;
#line 139 "ISIS_MERLIN.instr"
  mccfinal_Guide_my = mcipm;
#line 72 "ISIS_MERLIN.instr"
  mccfinal_Guide_nu = 0;
#line 72 "ISIS_MERLIN.instr"
  mccfinal_Guide_phase = 0;
#line 10720 "ISIS_MERLIN.c"

  SIG_MESSAGE("final_Guide (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 10727 "ISIS_MERLIN.c"
  rot_mul(mctr1, mcrotaOrigin, mcrotafinal_Guide);
  rot_transpose(mcrotaTOF_fermi, mctr1);
  rot_mul(mcrotafinal_Guide, mctr1, mcrotrfinal_Guide);
  mctc1 = coords_set(
#line 140 "ISIS_MERLIN.instr"
    0,
#line 140 "ISIS_MERLIN.instr"
    0,
#line 140 "ISIS_MERLIN.instr"
    10.277);
#line 10738 "ISIS_MERLIN.c"
  rot_transpose(mcrotaOrigin, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposafinal_Guide = coords_add(mcposaOrigin, mctc2);
  mctc1 = coords_sub(mcposaTOF_fermi, mcposafinal_Guide);
  mcposrfinal_Guide = rot_apply(mcrotafinal_Guide, mctc1);
  mcDEBUG_COMPONENT("final_Guide", mcposafinal_Guide, mcrotafinal_Guide)
  mccomp_posa[9] = mcposafinal_Guide;
  mccomp_posr[9] = mcposrfinal_Guide;
  mcNCounter[9]  = mcPCounter[9] = mcP2Counter[9] = 0;
  mcAbsorbProp[9]= 0;
    /* Component horz_div. */
  /* Setting parameters for component horz_div. */
  SIG_MESSAGE("horz_div (Init:SetPar)");
#line 143 "ISIS_MERLIN.instr"
  if("horz.div.m") strncpy(mcchorz_div_filename, "horz.div.m" ? "horz.div.m" : "", 16384); else mcchorz_div_filename[0]='\0';
#line 144 "ISIS_MERLIN.instr"
  mcchorz_div_xmin = -0.025;
#line 144 "ISIS_MERLIN.instr"
  mcchorz_div_xmax = 0.025;
#line 144 "ISIS_MERLIN.instr"
  mcchorz_div_ymin = -0.025;
#line 144 "ISIS_MERLIN.instr"
  mcchorz_div_ymax = 0.025;
#line 52 "ISIS_MERLIN.instr"
  mcchorz_div_xwidth = 0;
#line 52 "ISIS_MERLIN.instr"
  mcchorz_div_yheight = 0;
#line 145 "ISIS_MERLIN.instr"
  mcchorz_div_h_maxdiv = 0.3 * lam + .4;
#line 52 "ISIS_MERLIN.instr"
  mcchorz_div_restore_neutron = 0;
#line 52 "ISIS_MERLIN.instr"
  mcchorz_div_nowritefile = 0;
#line 10772 "ISIS_MERLIN.c"

  SIG_MESSAGE("horz_div (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 10779 "ISIS_MERLIN.c"
  rot_mul(mctr1, mcrotaOrigin, mcrotahorz_div);
  rot_transpose(mcrotafinal_Guide, mctr1);
  rot_mul(mcrotahorz_div, mctr1, mcrotrhorz_div);
  mctc1 = coords_set(
#line 146 "ISIS_MERLIN.instr"
    0,
#line 146 "ISIS_MERLIN.instr"
    0,
#line 146 "ISIS_MERLIN.instr"
    11.79);
#line 10790 "ISIS_MERLIN.c"
  rot_transpose(mcrotaOrigin, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposahorz_div = coords_add(mcposaOrigin, mctc2);
  mctc1 = coords_sub(mcposafinal_Guide, mcposahorz_div);
  mcposrhorz_div = rot_apply(mcrotahorz_div, mctc1);
  mcDEBUG_COMPONENT("horz_div", mcposahorz_div, mcrotahorz_div)
  mccomp_posa[10] = mcposahorz_div;
  mccomp_posr[10] = mcposrhorz_div;
  mcNCounter[10]  = mcPCounter[10] = mcP2Counter[10] = 0;
  mcAbsorbProp[10]= 0;
    /* Component PSD. */
  /* Setting parameters for component PSD. */
  SIG_MESSAGE("PSD (Init:SetPar)");
#line 150 "ISIS_MERLIN.instr"
  if("Output.psd") strncpy(mccPSD_filename, "Output.psd" ? "Output.psd" : "", 16384); else mccPSD_filename[0]='\0';
#line 149 "ISIS_MERLIN.instr"
  mccPSD_xmin = -0.04;
#line 149 "ISIS_MERLIN.instr"
  mccPSD_xmax = 0.04;
#line 149 "ISIS_MERLIN.instr"
  mccPSD_ymin = -0.04;
#line 149 "ISIS_MERLIN.instr"
  mccPSD_ymax = 0.04;
#line 50 "ISIS_MERLIN.instr"
  mccPSD_xwidth = 0;
#line 50 "ISIS_MERLIN.instr"
  mccPSD_yheight = 0;
#line 50 "ISIS_MERLIN.instr"
  mccPSD_restore_neutron = 0;
#line 50 "ISIS_MERLIN.instr"
  mccPSD_nowritefile = 0;
#line 10822 "ISIS_MERLIN.c"

  SIG_MESSAGE("PSD (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 10829 "ISIS_MERLIN.c"
  rot_mul(mctr1, mcrotaOrigin, mcrotaPSD);
  rot_transpose(mcrotahorz_div, mctr1);
  rot_mul(mcrotaPSD, mctr1, mcrotrPSD);
  mctc1 = coords_set(
#line 151 "ISIS_MERLIN.instr"
    0,
#line 151 "ISIS_MERLIN.instr"
    0,
#line 151 "ISIS_MERLIN.instr"
    11.795);
#line 10840 "ISIS_MERLIN.c"
  rot_transpose(mcrotaOrigin, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaPSD = coords_add(mcposaOrigin, mctc2);
  mctc1 = coords_sub(mcposahorz_div, mcposaPSD);
  mcposrPSD = rot_apply(mcrotaPSD, mctc1);
  mcDEBUG_COMPONENT("PSD", mcposaPSD, mcrotaPSD)
  mccomp_posa[11] = mcposaPSD;
  mccomp_posr[11] = mcposrPSD;
  mcNCounter[11]  = mcPCounter[11] = mcP2Counter[11] = 0;
  mcAbsorbProp[11]= 0;
    /* Component det_lmon. */
  /* Setting parameters for component det_lmon. */
  SIG_MESSAGE("det_lmon (Init:SetPar)");
#line 154 "ISIS_MERLIN.instr"
  if("det_Lmon.m") strncpy(mccdet_lmon_filename, "det_Lmon.m" ? "det_Lmon.m" : "", 16384); else mccdet_lmon_filename[0]='\0';
#line 153 "ISIS_MERLIN.instr"
  mccdet_lmon_xmin = -0.005;
#line 153 "ISIS_MERLIN.instr"
  mccdet_lmon_xmax = 0.005;
#line 154 "ISIS_MERLIN.instr"
  mccdet_lmon_ymin = -0.005;
#line 154 "ISIS_MERLIN.instr"
  mccdet_lmon_ymax = 0.005;
#line 51 "ISIS_MERLIN.instr"
  mccdet_lmon_xwidth = 0;
#line 51 "ISIS_MERLIN.instr"
  mccdet_lmon_yheight = 0;
#line 154 "ISIS_MERLIN.instr"
  mccdet_lmon_Lmin = 0.0;
#line 154 "ISIS_MERLIN.instr"
  mccdet_lmon_Lmax = 10;
#line 51 "ISIS_MERLIN.instr"
  mccdet_lmon_restore_neutron = 0;
#line 51 "ISIS_MERLIN.instr"
  mccdet_lmon_nowritefile = 0;
#line 10876 "ISIS_MERLIN.c"

  SIG_MESSAGE("det_lmon (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 10883 "ISIS_MERLIN.c"
  rot_mul(mctr1, mcrotaOrigin, mcrotadet_lmon);
  rot_transpose(mcrotaPSD, mctr1);
  rot_mul(mcrotadet_lmon, mctr1, mcrotrdet_lmon);
  mctc1 = coords_set(
#line 155 "ISIS_MERLIN.instr"
    0,
#line 155 "ISIS_MERLIN.instr"
    0,
#line 155 "ISIS_MERLIN.instr"
    11.796);
#line 10894 "ISIS_MERLIN.c"
  rot_transpose(mcrotaOrigin, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposadet_lmon = coords_add(mcposaOrigin, mctc2);
  mctc1 = coords_sub(mcposaPSD, mcposadet_lmon);
  mcposrdet_lmon = rot_apply(mcrotadet_lmon, mctc1);
  mcDEBUG_COMPONENT("det_lmon", mcposadet_lmon, mcrotadet_lmon)
  mccomp_posa[12] = mcposadet_lmon;
  mccomp_posr[12] = mcposrdet_lmon;
  mcNCounter[12]  = mcPCounter[12] = mcP2Counter[12] = 0;
  mcAbsorbProp[12]= 0;
    /* Component sample. */
  /* Setting parameters for component sample. */
  SIG_MESSAGE("sample (Init:SetPar)");
#line 118 "ISIS_MERLIN.instr"
  if(0) strncpy(mccsample_geometry, 0 ? 0 : "", 16384); else mccsample_geometry[0]='\0';
#line 158 "ISIS_MERLIN.instr"
  mccsample_radius = 0.02;
#line 118 "ISIS_MERLIN.instr"
  mccsample_xwidth = 0;
#line 158 "ISIS_MERLIN.instr"
  mccsample_yheight = 0.04;
#line 118 "ISIS_MERLIN.instr"
  mccsample_zdepth = 0;
#line 158 "ISIS_MERLIN.instr"
  mccsample_thickness = 0.001;
#line 119 "ISIS_MERLIN.instr"
  mccsample_target_x = 0;
#line 119 "ISIS_MERLIN.instr"
  mccsample_target_y = 0;
#line 119 "ISIS_MERLIN.instr"
  mccsample_target_z = 0;
#line 119 "ISIS_MERLIN.instr"
  mccsample_focus_r = 0;
#line 159 "ISIS_MERLIN.instr"
  mccsample_focus_xw = 0.025;
#line 159 "ISIS_MERLIN.instr"
  mccsample_focus_yh = 0.1;
#line 120 "ISIS_MERLIN.instr"
  mccsample_focus_aw = 0;
#line 120 "ISIS_MERLIN.instr"
  mccsample_focus_ah = 0;
#line 159 "ISIS_MERLIN.instr"
  mccsample_target_index = 2;
#line 121 "ISIS_MERLIN.instr"
  mccsample_pack = 1;
#line 121 "ISIS_MERLIN.instr"
  mccsample_p_interact = 1;
#line 121 "ISIS_MERLIN.instr"
  mccsample_f_QE = 0;
#line 121 "ISIS_MERLIN.instr"
  mccsample_gamma = 0;
#line 121 "ISIS_MERLIN.instr"
  mccsample_Etrans = 0;
#line 121 "ISIS_MERLIN.instr"
  mccsample_deltaE = 0;
#line 159 "ISIS_MERLIN.instr"
  mccsample_sigma_abs = 0.0;
#line 159 "ISIS_MERLIN.instr"
  mccsample_sigma_inc = 100;
#line 122 "ISIS_MERLIN.instr"
  mccsample_Vc = 13.827;
#line 122 "ISIS_MERLIN.instr"
  mccsample_concentric = 0;
#line 122 "ISIS_MERLIN.instr"
  mccsample_order = 0;
#line 10960 "ISIS_MERLIN.c"

  SIG_MESSAGE("sample (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 10967 "ISIS_MERLIN.c"
  rot_mul(mctr1, mcrotaOrigin, mcrotasample);
  rot_transpose(mcrotadet_lmon, mctr1);
  rot_mul(mcrotasample, mctr1, mcrotrsample);
  mctc1 = coords_set(
#line 160 "ISIS_MERLIN.instr"
    0,
#line 160 "ISIS_MERLIN.instr"
    0,
#line 160 "ISIS_MERLIN.instr"
    11.8);
#line 10978 "ISIS_MERLIN.c"
  rot_transpose(mcrotaOrigin, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposasample = coords_add(mcposaOrigin, mctc2);
  mctc1 = coords_sub(mcposadet_lmon, mcposasample);
  mcposrsample = rot_apply(mcrotasample, mctc1);
  mcDEBUG_COMPONENT("sample", mcposasample, mcrotasample)
  mccomp_posa[13] = mcposasample;
  mccomp_posr[13] = mcposrsample;
  mcNCounter[13]  = mcPCounter[13] = mcP2Counter[13] = 0;
  mcAbsorbProp[13]= 0;
    /* Component detectorarm. */
  /* Setting parameters for component detectorarm. */
  SIG_MESSAGE("detectorarm (Init:SetPar)");

  SIG_MESSAGE("detectorarm (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 165 "ISIS_MERLIN.instr"
    (0)*DEG2RAD,
#line 165 "ISIS_MERLIN.instr"
    (20.0)*DEG2RAD,
#line 165 "ISIS_MERLIN.instr"
    (0)*DEG2RAD);
#line 11001 "ISIS_MERLIN.c"
  rot_mul(mctr1, mcrotasample, mcrotadetectorarm);
  rot_transpose(mcrotasample, mctr1);
  rot_mul(mcrotadetectorarm, mctr1, mcrotrdetectorarm);
  mctc1 = coords_set(
#line 164 "ISIS_MERLIN.instr"
    0,
#line 164 "ISIS_MERLIN.instr"
    0,
#line 164 "ISIS_MERLIN.instr"
    0);
#line 11012 "ISIS_MERLIN.c"
  rot_transpose(mcrotasample, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposadetectorarm = coords_add(mcposasample, mctc2);
  mctc1 = coords_sub(mcposasample, mcposadetectorarm);
  mcposrdetectorarm = rot_apply(mcrotadetectorarm, mctc1);
  mcDEBUG_COMPONENT("detectorarm", mcposadetectorarm, mcrotadetectorarm)
  mccomp_posa[14] = mcposadetectorarm;
  mccomp_posr[14] = mcposrdetectorarm;
  mcNCounter[14]  = mcPCounter[14] = mcP2Counter[14] = 0;
  mcAbsorbProp[14]= 0;
    /* Component detector. */
  /* Setting parameters for component detector. */
  SIG_MESSAGE("detector (Init:SetPar)");
#line 169 "ISIS_MERLIN.instr"
  if("detector.m") strncpy(mccdetector_filename, "detector.m" ? "detector.m" : "", 16384); else mccdetector_filename[0]='\0';
#line 169 "ISIS_MERLIN.instr"
  mccdetector_xmin = -0.0125;
#line 170 "ISIS_MERLIN.instr"
  mccdetector_xmax = 0.0125;
#line 170 "ISIS_MERLIN.instr"
  mccdetector_ymin = -0.05;
#line 170 "ISIS_MERLIN.instr"
  mccdetector_ymax = 0.05;
#line 48 "ISIS_MERLIN.instr"
  mccdetector_xwidth = 0;
#line 48 "ISIS_MERLIN.instr"
  mccdetector_yheight = 0;
#line 170 "ISIS_MERLIN.instr"
  mccdetector_tmin = ( 1e6 * t_det ) -20;
#line 170 "ISIS_MERLIN.instr"
  mccdetector_tmax = ( 1e6 * t_det ) + 30;
#line 48 "ISIS_MERLIN.instr"
  mccdetector_dt = 1.0;
#line 48 "ISIS_MERLIN.instr"
  mccdetector_restore_neutron = 0;
#line 48 "ISIS_MERLIN.instr"
  mccdetector_nowritefile = 0;
#line 11050 "ISIS_MERLIN.c"

  SIG_MESSAGE("detector (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 11057 "ISIS_MERLIN.c"
  rot_mul(mctr1, mcrotadetectorarm, mcrotadetector);
  rot_transpose(mcrotadetectorarm, mctr1);
  rot_mul(mcrotadetector, mctr1, mcrotrdetector);
  mctc1 = coords_set(
#line 171 "ISIS_MERLIN.instr"
    0,
#line 171 "ISIS_MERLIN.instr"
    0,
#line 171 "ISIS_MERLIN.instr"
    2.5);
#line 11068 "ISIS_MERLIN.c"
  rot_transpose(mcrotadetectorarm, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposadetector = coords_add(mcposadetectorarm, mctc2);
  mctc1 = coords_sub(mcposadetectorarm, mcposadetector);
  mcposrdetector = rot_apply(mcrotadetector, mctc1);
  mcDEBUG_COMPONENT("detector", mcposadetector, mcrotadetector)
  mccomp_posa[15] = mcposadetector;
  mccomp_posr[15] = mcposrdetector;
  mcNCounter[15]  = mcPCounter[15] = mcP2Counter[15] = 0;
  mcAbsorbProp[15]= 0;
  /* Component initializations. */
  /* Initializations for component Origin. */
  SIG_MESSAGE("Origin (Init)");

  /* Initializations for component isis_mod. */
  SIG_MESSAGE("isis_mod (Init)");
#define mccompcurname  isis_mod
#define mccompcurtype  ViewModISIS
#define mccompcurindex 2
#define Face mccisis_mod_Face
#define E0 mccisis_mod_E0
#define E1 mccisis_mod_E1
#define modPosition mccisis_mod_modPosition
#define xw mccisis_mod_xw
#define yh mccisis_mod_yh
#define focus_xw mccisis_mod_focus_xw
#define focus_yh mccisis_mod_focus_yh
#define dist mccisis_mod_dist
#define verbose mccisis_mod_verbose
#define beamcurrent mccisis_mod_beamcurrent
#line 1043 "/usr/share/mcstas/2.5/contrib/ViewModISIS.comp"
{
  /* READ IN THE ENERGY FILE */
  FILE* Tfile;
  

  Nsim=mcget_ncount();  // Number of points requested.

  Tfile=openFile(Face);	              // Get open file
  rtE0=convertEnergy(E0);
  rtE1=convertEnergy(E1);
  orderEnergy(&rtE0,&rtE1);

  readHtable(Tfile,rtE0,rtE1);
  fclose(Tfile);

  /**********************************************************************/

  Tnpts=0;
  Ncount=0;

  fprintf(stderr,"Face == %s \n",Face);
  fprintf(stderr,"Number of Energy Points == %d\n",TS.nEnergy);
  if (dist<0)
    {
      dist=TS.rdumMid;
      fprintf(stderr,"Setting distance to moderator surface == %g\n",
	      dist);
    }
  /* Do solid angle correction */
  angleArea= strArea();

  fprintf(stderr,"Totals:: %g %d %d \n",TS.Total,TS.nEnergy,TS.nTime);

}
#line 11134 "ISIS_MERLIN.c"
#undef beamcurrent
#undef verbose
#undef dist
#undef focus_yh
#undef focus_xw
#undef yh
#undef xw
#undef modPosition
#undef E1
#undef E0
#undef Face
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component TOF_mod. */
  SIG_MESSAGE("TOF_mod (Init)");
#define mccompcurname  TOF_mod
#define mccompcurtype  TOF_monitor
#define mccompcurindex 3
#define nt mccTOF_mod_nt
#define TOF_N mccTOF_mod_TOF_N
#define TOF_p mccTOF_mod_TOF_p
#define TOF_p2 mccTOF_mod_TOF_p2
#define t_min mccTOF_mod_t_min
#define t_max mccTOF_mod_t_max
#define delta_t mccTOF_mod_delta_t
#define filename mccTOF_mod_filename
#define xmin mccTOF_mod_xmin
#define xmax mccTOF_mod_xmax
#define ymin mccTOF_mod_ymin
#define ymax mccTOF_mod_ymax
#define xwidth mccTOF_mod_xwidth
#define yheight mccTOF_mod_yheight
#define tmin mccTOF_mod_tmin
#define tmax mccTOF_mod_tmax
#define dt mccTOF_mod_dt
#define restore_neutron mccTOF_mod_restore_neutron
#define nowritefile mccTOF_mod_nowritefile
#line 61 "/usr/share/mcstas/2.5/monitors/TOF_monitor.comp"
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
#line 11207 "ISIS_MERLIN.c"
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

  /* Initializations for component shutter_guide. */
  SIG_MESSAGE("shutter_guide (Init)");
#define mccompcurname  shutter_guide
#define mccompcurtype  Guide_channeled
#define mccompcurindex 4
#define w1c mccshutter_guide_w1c
#define w2c mccshutter_guide_w2c
#define ww mccshutter_guide_ww
#define hh mccshutter_guide_hh
#define whalf mccshutter_guide_whalf
#define hhalf mccshutter_guide_hhalf
#define lwhalf mccshutter_guide_lwhalf
#define lhhalf mccshutter_guide_lhhalf
#define w1 mccshutter_guide_w1
#define h1 mccshutter_guide_h1
#define w2 mccshutter_guide_w2
#define h2 mccshutter_guide_h2
#define l mccshutter_guide_l
#define R0 mccshutter_guide_R0
#define Qc mccshutter_guide_Qc
#define alpha mccshutter_guide_alpha
#define m mccshutter_guide_m
#define nslit mccshutter_guide_nslit
#define d mccshutter_guide_d
#define Qcx mccshutter_guide_Qcx
#define Qcy mccshutter_guide_Qcy
#define alphax mccshutter_guide_alphax
#define alphay mccshutter_guide_alphay
#define W mccshutter_guide_W
#define mx mccshutter_guide_mx
#define my mccshutter_guide_my
#define nu mccshutter_guide_nu
#define phase mccshutter_guide_phase
#line 88 "/usr/share/mcstas/2.5/optics/Guide_channeled.comp"
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
#line 11305 "ISIS_MERLIN.c"
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

  /* Initializations for component conv_guide. */
  SIG_MESSAGE("conv_guide (Init)");
#define mccompcurname  conv_guide
#define mccompcurtype  Guide_channeled
#define mccompcurindex 5
#define w1c mccconv_guide_w1c
#define w2c mccconv_guide_w2c
#define ww mccconv_guide_ww
#define hh mccconv_guide_hh
#define whalf mccconv_guide_whalf
#define hhalf mccconv_guide_hhalf
#define lwhalf mccconv_guide_lwhalf
#define lhhalf mccconv_guide_lhhalf
#define w1 mccconv_guide_w1
#define h1 mccconv_guide_h1
#define w2 mccconv_guide_w2
#define h2 mccconv_guide_h2
#define l mccconv_guide_l
#define R0 mccconv_guide_R0
#define Qc mccconv_guide_Qc
#define alpha mccconv_guide_alpha
#define m mccconv_guide_m
#define nslit mccconv_guide_nslit
#define d mccconv_guide_d
#define Qcx mccconv_guide_Qcx
#define Qcy mccconv_guide_Qcy
#define alphax mccconv_guide_alphax
#define alphay mccconv_guide_alphay
#define W mccconv_guide_W
#define mx mccconv_guide_mx
#define my mccconv_guide_my
#define nu mccconv_guide_nu
#define phase mccconv_guide_phase
#line 88 "/usr/share/mcstas/2.5/optics/Guide_channeled.comp"
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
#line 11412 "ISIS_MERLIN.c"
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

  /* Initializations for component between_guide. */
  SIG_MESSAGE("between_guide (Init)");
#define mccompcurname  between_guide
#define mccompcurtype  Guide_channeled
#define mccompcurindex 6
#define w1c mccbetween_guide_w1c
#define w2c mccbetween_guide_w2c
#define ww mccbetween_guide_ww
#define hh mccbetween_guide_hh
#define whalf mccbetween_guide_whalf
#define hhalf mccbetween_guide_hhalf
#define lwhalf mccbetween_guide_lwhalf
#define lhhalf mccbetween_guide_lhhalf
#define w1 mccbetween_guide_w1
#define h1 mccbetween_guide_h1
#define w2 mccbetween_guide_w2
#define h2 mccbetween_guide_h2
#define l mccbetween_guide_l
#define R0 mccbetween_guide_R0
#define Qc mccbetween_guide_Qc
#define alpha mccbetween_guide_alpha
#define m mccbetween_guide_m
#define nslit mccbetween_guide_nslit
#define d mccbetween_guide_d
#define Qcx mccbetween_guide_Qcx
#define Qcy mccbetween_guide_Qcy
#define alphax mccbetween_guide_alphax
#define alphay mccbetween_guide_alphay
#define W mccbetween_guide_W
#define mx mccbetween_guide_mx
#define my mccbetween_guide_my
#define nu mccbetween_guide_nu
#define phase mccbetween_guide_phase
#line 88 "/usr/share/mcstas/2.5/optics/Guide_channeled.comp"
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
#line 11519 "ISIS_MERLIN.c"
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

  /* Initializations for component fermi. */
  SIG_MESSAGE("fermi (Init)");
#define mccompcurname  fermi
#define mccompcurtype  FermiChopper
#define mccompcurindex 7
#define FCVars mccfermi_FCVars
#define phase mccfermi_phase
#define radius mccfermi_radius
#define nu mccfermi_nu
#define w mccfermi_w
#define nslit mccfermi_nslit
#define R0 mccfermi_R0
#define Qc mccfermi_Qc
#define alpha mccfermi_alpha
#define m mccfermi_m
#define W mccfermi_W
#define length mccfermi_length
#define eff mccfermi_eff
#define zero_time mccfermi_zero_time
#define xwidth mccfermi_xwidth
#define verbose mccfermi_verbose
#define yheight mccfermi_yheight
#define curvature mccfermi_curvature
#define delay mccfermi_delay
#line 282 "/usr/share/mcstas/2.5/optics/FermiChopper.comp"
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
#line 11654 "ISIS_MERLIN.c"
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

  /* Initializations for component TOF_fermi. */
  SIG_MESSAGE("TOF_fermi (Init)");
#define mccompcurname  TOF_fermi
#define mccompcurtype  TOF_monitor
#define mccompcurindex 8
#define nt mccTOF_fermi_nt
#define TOF_N mccTOF_fermi_TOF_N
#define TOF_p mccTOF_fermi_TOF_p
#define TOF_p2 mccTOF_fermi_TOF_p2
#define t_min mccTOF_fermi_t_min
#define t_max mccTOF_fermi_t_max
#define delta_t mccTOF_fermi_delta_t
#define filename mccTOF_fermi_filename
#define xmin mccTOF_fermi_xmin
#define xmax mccTOF_fermi_xmax
#define ymin mccTOF_fermi_ymin
#define ymax mccTOF_fermi_ymax
#define xwidth mccTOF_fermi_xwidth
#define yheight mccTOF_fermi_yheight
#define tmin mccTOF_fermi_tmin
#define tmax mccTOF_fermi_tmax
#define dt mccTOF_fermi_dt
#define restore_neutron mccTOF_fermi_restore_neutron
#define nowritefile mccTOF_fermi_nowritefile
#line 61 "/usr/share/mcstas/2.5/monitors/TOF_monitor.comp"
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
#line 11735 "ISIS_MERLIN.c"
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

  /* Initializations for component final_Guide. */
  SIG_MESSAGE("final_Guide (Init)");
#define mccompcurname  final_Guide
#define mccompcurtype  Guide_channeled
#define mccompcurindex 9
#define w1c mccfinal_Guide_w1c
#define w2c mccfinal_Guide_w2c
#define ww mccfinal_Guide_ww
#define hh mccfinal_Guide_hh
#define whalf mccfinal_Guide_whalf
#define hhalf mccfinal_Guide_hhalf
#define lwhalf mccfinal_Guide_lwhalf
#define lhhalf mccfinal_Guide_lhhalf
#define w1 mccfinal_Guide_w1
#define h1 mccfinal_Guide_h1
#define w2 mccfinal_Guide_w2
#define h2 mccfinal_Guide_h2
#define l mccfinal_Guide_l
#define R0 mccfinal_Guide_R0
#define Qc mccfinal_Guide_Qc
#define alpha mccfinal_Guide_alpha
#define m mccfinal_Guide_m
#define nslit mccfinal_Guide_nslit
#define d mccfinal_Guide_d
#define Qcx mccfinal_Guide_Qcx
#define Qcy mccfinal_Guide_Qcy
#define alphax mccfinal_Guide_alphax
#define alphay mccfinal_Guide_alphay
#define W mccfinal_Guide_W
#define mx mccfinal_Guide_mx
#define my mccfinal_Guide_my
#define nu mccfinal_Guide_nu
#define phase mccfinal_Guide_phase
#line 88 "/usr/share/mcstas/2.5/optics/Guide_channeled.comp"
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
#line 11833 "ISIS_MERLIN.c"
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

  /* Initializations for component horz_div. */
  SIG_MESSAGE("horz_div (Init)");
#define mccompcurname  horz_div
#define mccompcurtype  Hdiv_monitor
#define mccompcurindex 10
#define nh mcchorz_div_nh
#define Div_N mcchorz_div_Div_N
#define Div_p mcchorz_div_Div_p
#define Div_p2 mcchorz_div_Div_p2
#define filename mcchorz_div_filename
#define xmin mcchorz_div_xmin
#define xmax mcchorz_div_xmax
#define ymin mcchorz_div_ymin
#define ymax mcchorz_div_ymax
#define xwidth mcchorz_div_xwidth
#define yheight mcchorz_div_yheight
#define h_maxdiv mcchorz_div_h_maxdiv
#define restore_neutron mcchorz_div_restore_neutron
#define nowritefile mcchorz_div_nowritefile
#line 62 "/usr/share/mcstas/2.5/monitors/Hdiv_monitor.comp"
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
#line 11908 "ISIS_MERLIN.c"
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

  /* Initializations for component PSD. */
  SIG_MESSAGE("PSD (Init)");
#define mccompcurname  PSD
#define mccompcurtype  PSD_monitor
#define mccompcurindex 11
#define nx mccPSD_nx
#define ny mccPSD_ny
#define PSD_N mccPSD_PSD_N
#define PSD_p mccPSD_PSD_p
#define PSD_p2 mccPSD_PSD_p2
#define filename mccPSD_filename
#define xmin mccPSD_xmin
#define xmax mccPSD_xmax
#define ymin mccPSD_ymin
#define ymax mccPSD_ymax
#define xwidth mccPSD_xwidth
#define yheight mccPSD_yheight
#define restore_neutron mccPSD_restore_neutron
#define nowritefile mccPSD_nowritefile
#line 61 "/usr/share/mcstas/2.5/monitors/PSD_monitor.comp"
{
int i,j;

if (xwidth  > 0) { xmax = xwidth/2;  xmin = -xmax; }
    if (yheight > 0) { ymax = yheight/2; ymin = -ymax; }

    if ((xmin >= xmax) || (ymin >= ymax)) {
            printf("PSD_monitor: %s: Null detection area !\n"
                   "ERROR        (xwidth,yheight,xmin,xmax,ymin,ymax). Exiting",
           NAME_CURRENT_COMP);
      exit(0);
    }

    for (i=0; i<nx; i++)
     for (j=0; j<ny; j++)
     {
      PSD_N[i][j] = 0;
      PSD_p[i][j] = 0;
      PSD_p2[i][j] = 0;
     }
}
#line 11968 "ISIS_MERLIN.c"
#undef nowritefile
#undef restore_neutron
#undef yheight
#undef xwidth
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef filename
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef ny
#undef nx
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component det_lmon. */
  SIG_MESSAGE("det_lmon (Init)");
#define mccompcurname  det_lmon
#define mccompcurtype  L_monitor
#define mccompcurindex 12
#define nL mccdet_lmon_nL
#define L_N mccdet_lmon_L_N
#define L_p mccdet_lmon_L_p
#define L_p2 mccdet_lmon_L_p2
#define filename mccdet_lmon_filename
#define xmin mccdet_lmon_xmin
#define xmax mccdet_lmon_xmax
#define ymin mccdet_lmon_ymin
#define ymax mccdet_lmon_ymax
#define xwidth mccdet_lmon_xwidth
#define yheight mccdet_lmon_yheight
#define Lmin mccdet_lmon_Lmin
#define Lmax mccdet_lmon_Lmax
#define restore_neutron mccdet_lmon_restore_neutron
#define nowritefile mccdet_lmon_nowritefile
#line 62 "/usr/share/mcstas/2.5/monitors/L_monitor.comp"
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
#line 12028 "ISIS_MERLIN.c"
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
#define mccompcurtype  Incoherent
#define mccompcurindex 13
#define VarsInc mccsample_VarsInc
#define offdata mccsample_offdata
#define geometry mccsample_geometry
#define radius mccsample_radius
#define xwidth mccsample_xwidth
#define yheight mccsample_yheight
#define zdepth mccsample_zdepth
#define thickness mccsample_thickness
#define target_x mccsample_target_x
#define target_y mccsample_target_y
#define target_z mccsample_target_z
#define focus_r mccsample_focus_r
#define focus_xw mccsample_focus_xw
#define focus_yh mccsample_focus_yh
#define focus_aw mccsample_focus_aw
#define focus_ah mccsample_focus_ah
#define target_index mccsample_target_index
#define pack mccsample_pack
#define p_interact mccsample_p_interact
#define f_QE mccsample_f_QE
#define gamma mccsample_gamma
#define Etrans mccsample_Etrans
#define deltaE mccsample_deltaE
#define sigma_abs mccsample_sigma_abs
#define sigma_inc mccsample_sigma_inc
#define Vc mccsample_Vc
#define concentric mccsample_concentric
#define order mccsample_order
#line 155 "/usr/share/mcstas/2.5/samples/Incoherent.comp"
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
#line 12156 "ISIS_MERLIN.c"
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

  /* Initializations for component detectorarm. */
  SIG_MESSAGE("detectorarm (Init)");

  /* Initializations for component detector. */
  SIG_MESSAGE("detector (Init)");
#define mccompcurname  detector
#define mccompcurtype  TOF_monitor
#define mccompcurindex 15
#define nt mccdetector_nt
#define TOF_N mccdetector_TOF_N
#define TOF_p mccdetector_TOF_p
#define TOF_p2 mccdetector_TOF_p2
#define t_min mccdetector_t_min
#define t_max mccdetector_t_max
#define delta_t mccdetector_delta_t
#define filename mccdetector_filename
#define xmin mccdetector_xmin
#define xmax mccdetector_xmax
#define ymin mccdetector_ymin
#define ymax mccdetector_ymax
#define xwidth mccdetector_xwidth
#define yheight mccdetector_yheight
#define tmin mccdetector_tmin
#define tmax mccdetector_tmax
#define dt mccdetector_dt
#define restore_neutron mccdetector_restore_neutron
#define nowritefile mccdetector_nowritefile
#line 61 "/usr/share/mcstas/2.5/monitors/TOF_monitor.comp"
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
#line 12249 "ISIS_MERLIN.c"
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

  /* TRACE Component isis_mod [2] */
  mccoordschange(mcposrisis_mod, mcrotrisis_mod,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component isis_mod (without coords transformations) */
  mcJumpTrace_isis_mod:
  SIG_MESSAGE("isis_mod (Trace)");
  mcDEBUG_COMP("isis_mod")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompisis_mod
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
#define mccompcurname  isis_mod
#define mccompcurtype  ViewModISIS
#define mccompcurindex 2
{   /* Declarations of isis_mod=ViewModISIS() SETTING parameters. */
char* Face = mccisis_mod_Face;
MCNUM E0 = mccisis_mod_E0;
MCNUM E1 = mccisis_mod_E1;
MCNUM modPosition = mccisis_mod_modPosition;
MCNUM xw = mccisis_mod_xw;
MCNUM yh = mccisis_mod_yh;
MCNUM focus_xw = mccisis_mod_focus_xw;
MCNUM focus_yh = mccisis_mod_focus_yh;
MCNUM dist = mccisis_mod_dist;
int verbose = mccisis_mod_verbose;
MCNUM beamcurrent = mccisis_mod_beamcurrent;
#line 1079 "/usr/share/mcstas/2.5/contrib/ViewModISIS.comp"
{
  double v,r,E;
  double xf,yf,dx,dy;    /* mxp ->max var in param space */
  double Ival,Tval,Eval;

  Ncount++;

  x += TS.XAxis*(0.5-rand01());            /* Get point on shutter enterance */
  y += TS.ZAxis*(0.5-rand01());            /* Get point on shutter enterance */
    
  // if (modPosition) z=TS.rdumMid;    //   Goran          

  xf = focus_xw*(0.5-rand01());          /* Choose focusing position uniformly */
  yf = focus_yh*(0.5-rand01());          /* Choose focusing position uniformly */

  getPoint(&Tval,&Eval,&rtE0,&rtE1);
  
//  Ival=TS.Total*6.2415093e+12; // * 1.1879451;  /* ( of proton in 1uAmp) * (1-cos(30))*2*Pi  */
//  Ival=TS.Total*5.9294338e+12; // Number of proton in 1uAmp * 0.95 (0.95 because of muon target loss)
  Ival=TS.Total*6.2415093e+12; // Number of proton in 1uAmp

  dx = xf-x;
  dy = yf-y;
  r = sqrt(dx*dx+dy*dy+dist*dist);      // Actual distance to point 
  v = SE2V*sqrt(Eval);                  // Calculate the velocity 
  vz = v*fabs(dist)/r;
  vy = v*dy/r;
  vx = v*dx/r;


   t=Tval-(TS.rdumMid-TS.timeOffset)/vz; ///vz;
   
     if (modPosition)
     {
     t+=TS.rdumMid/vz;
     }

  p=beamcurrent*angleArea*Ival/Nsim;

  //   testing only - Goran              
  //  if (!modPosition) p*=TS.timeOffset;


  if (verbose && !(Ncount % 100000))
    {
      fprintf(stderr,"FPos[%d]=> %g %g %g  \n",Ncount,x,y,z);
      fprintf(stderr,"FDir[%d]=> %g %g %g  \n",Ncount,vx,vy,vz);
      fprintf(stderr,"PlaneAxis %g %g \n",TS.XAxis,fullAngle);
      fprintf(stderr,"RMID %g \n",TS.rdumMid);
      fprintf(stderr,"TimeMid[%d]=> %g\n",Ncount,TS.rdumMid);
      fprintf(stderr,"TimeOffset[%d]=> %g\n",Ncount,TS.timeOffset);
      fprintf(stderr,"TimeTval[%d]=> %g\n",Ncount,Tval);
      fprintf(stderr,"TimeZero[%d]=> %g\n",Ncount,t);
    }  

}
#line 12545 "ISIS_MERLIN.c"
}   /* End of isis_mod=ViewModISIS() SETTING parameter declarations. */
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompisis_mod:
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

  /* TRACE Component TOF_mod [3] */
  mccoordschange(mcposrTOF_mod, mcrotrTOF_mod,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component TOF_mod (without coords transformations) */
  mcJumpTrace_TOF_mod:
  SIG_MESSAGE("TOF_mod (Trace)");
  mcDEBUG_COMP("TOF_mod")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompTOF_mod
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
#define mccompcurname  TOF_mod
#define mccompcurtype  TOF_monitor
#define mccompcurindex 3
#define nt mccTOF_mod_nt
#define TOF_N mccTOF_mod_TOF_N
#define TOF_p mccTOF_mod_TOF_p
#define TOF_p2 mccTOF_mod_TOF_p2
#define t_min mccTOF_mod_t_min
#define t_max mccTOF_mod_t_max
#define delta_t mccTOF_mod_delta_t
{   /* Declarations of TOF_mod=TOF_monitor() SETTING parameters. */
char* filename = mccTOF_mod_filename;
MCNUM xmin = mccTOF_mod_xmin;
MCNUM xmax = mccTOF_mod_xmax;
MCNUM ymin = mccTOF_mod_ymin;
MCNUM ymax = mccTOF_mod_ymax;
MCNUM xwidth = mccTOF_mod_xwidth;
MCNUM yheight = mccTOF_mod_yheight;
MCNUM tmin = mccTOF_mod_tmin;
MCNUM tmax = mccTOF_mod_tmax;
MCNUM dt = mccTOF_mod_dt;
MCNUM restore_neutron = mccTOF_mod_restore_neutron;
int nowritefile = mccTOF_mod_nowritefile;
#line 95 "/usr/share/mcstas/2.5/monitors/TOF_monitor.comp"
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
#line 12689 "ISIS_MERLIN.c"
}   /* End of TOF_mod=TOF_monitor() SETTING parameter declarations. */
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
  mcabsorbCompTOF_mod:
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

  /* TRACE Component shutter_guide [4] */
  mccoordschange(mcposrshutter_guide, mcrotrshutter_guide,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component shutter_guide (without coords transformations) */
  mcJumpTrace_shutter_guide:
  SIG_MESSAGE("shutter_guide (Trace)");
  mcDEBUG_COMP("shutter_guide")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompshutter_guide
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
#define mccompcurname  shutter_guide
#define mccompcurtype  Guide_channeled
#define mccompcurindex 4
#define w1c mccshutter_guide_w1c
#define w2c mccshutter_guide_w2c
#define ww mccshutter_guide_ww
#define hh mccshutter_guide_hh
#define whalf mccshutter_guide_whalf
#define hhalf mccshutter_guide_hhalf
#define lwhalf mccshutter_guide_lwhalf
#define lhhalf mccshutter_guide_lhhalf
{   /* Declarations of shutter_guide=Guide_channeled() SETTING parameters. */
MCNUM w1 = mccshutter_guide_w1;
MCNUM h1 = mccshutter_guide_h1;
MCNUM w2 = mccshutter_guide_w2;
MCNUM h2 = mccshutter_guide_h2;
MCNUM l = mccshutter_guide_l;
MCNUM R0 = mccshutter_guide_R0;
MCNUM Qc = mccshutter_guide_Qc;
MCNUM alpha = mccshutter_guide_alpha;
MCNUM m = mccshutter_guide_m;
MCNUM nslit = mccshutter_guide_nslit;
MCNUM d = mccshutter_guide_d;
MCNUM Qcx = mccshutter_guide_Qcx;
MCNUM Qcy = mccshutter_guide_Qcy;
MCNUM alphax = mccshutter_guide_alphax;
MCNUM alphay = mccshutter_guide_alphay;
MCNUM W = mccshutter_guide_W;
MCNUM mx = mccshutter_guide_mx;
MCNUM my = mccshutter_guide_my;
MCNUM nu = mccshutter_guide_nu;
MCNUM phase = mccshutter_guide_phase;
#line 130 "/usr/share/mcstas/2.5/optics/Guide_channeled.comp"
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
#line 12995 "ISIS_MERLIN.c"
}   /* End of shutter_guide=Guide_channeled() SETTING parameter declarations. */
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
  mcabsorbCompshutter_guide:
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

  /* TRACE Component conv_guide [5] */
  mccoordschange(mcposrconv_guide, mcrotrconv_guide,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component conv_guide (without coords transformations) */
  mcJumpTrace_conv_guide:
  SIG_MESSAGE("conv_guide (Trace)");
  mcDEBUG_COMP("conv_guide")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompconv_guide
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
#define mccompcurname  conv_guide
#define mccompcurtype  Guide_channeled
#define mccompcurindex 5
#define w1c mccconv_guide_w1c
#define w2c mccconv_guide_w2c
#define ww mccconv_guide_ww
#define hh mccconv_guide_hh
#define whalf mccconv_guide_whalf
#define hhalf mccconv_guide_hhalf
#define lwhalf mccconv_guide_lwhalf
#define lhhalf mccconv_guide_lhhalf
{   /* Declarations of conv_guide=Guide_channeled() SETTING parameters. */
MCNUM w1 = mccconv_guide_w1;
MCNUM h1 = mccconv_guide_h1;
MCNUM w2 = mccconv_guide_w2;
MCNUM h2 = mccconv_guide_h2;
MCNUM l = mccconv_guide_l;
MCNUM R0 = mccconv_guide_R0;
MCNUM Qc = mccconv_guide_Qc;
MCNUM alpha = mccconv_guide_alpha;
MCNUM m = mccconv_guide_m;
MCNUM nslit = mccconv_guide_nslit;
MCNUM d = mccconv_guide_d;
MCNUM Qcx = mccconv_guide_Qcx;
MCNUM Qcy = mccconv_guide_Qcy;
MCNUM alphax = mccconv_guide_alphax;
MCNUM alphay = mccconv_guide_alphay;
MCNUM W = mccconv_guide_W;
MCNUM mx = mccconv_guide_mx;
MCNUM my = mccconv_guide_my;
MCNUM nu = mccconv_guide_nu;
MCNUM phase = mccconv_guide_phase;
#line 130 "/usr/share/mcstas/2.5/optics/Guide_channeled.comp"
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
#line 13302 "ISIS_MERLIN.c"
}   /* End of conv_guide=Guide_channeled() SETTING parameter declarations. */
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
  mcabsorbCompconv_guide:
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

  /* TRACE Component between_guide [6] */
  mccoordschange(mcposrbetween_guide, mcrotrbetween_guide,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component between_guide (without coords transformations) */
  mcJumpTrace_between_guide:
  SIG_MESSAGE("between_guide (Trace)");
  mcDEBUG_COMP("between_guide")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompbetween_guide
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
#define mccompcurname  between_guide
#define mccompcurtype  Guide_channeled
#define mccompcurindex 6
#define w1c mccbetween_guide_w1c
#define w2c mccbetween_guide_w2c
#define ww mccbetween_guide_ww
#define hh mccbetween_guide_hh
#define whalf mccbetween_guide_whalf
#define hhalf mccbetween_guide_hhalf
#define lwhalf mccbetween_guide_lwhalf
#define lhhalf mccbetween_guide_lhhalf
{   /* Declarations of between_guide=Guide_channeled() SETTING parameters. */
MCNUM w1 = mccbetween_guide_w1;
MCNUM h1 = mccbetween_guide_h1;
MCNUM w2 = mccbetween_guide_w2;
MCNUM h2 = mccbetween_guide_h2;
MCNUM l = mccbetween_guide_l;
MCNUM R0 = mccbetween_guide_R0;
MCNUM Qc = mccbetween_guide_Qc;
MCNUM alpha = mccbetween_guide_alpha;
MCNUM m = mccbetween_guide_m;
MCNUM nslit = mccbetween_guide_nslit;
MCNUM d = mccbetween_guide_d;
MCNUM Qcx = mccbetween_guide_Qcx;
MCNUM Qcy = mccbetween_guide_Qcy;
MCNUM alphax = mccbetween_guide_alphax;
MCNUM alphay = mccbetween_guide_alphay;
MCNUM W = mccbetween_guide_W;
MCNUM mx = mccbetween_guide_mx;
MCNUM my = mccbetween_guide_my;
MCNUM nu = mccbetween_guide_nu;
MCNUM phase = mccbetween_guide_phase;
#line 130 "/usr/share/mcstas/2.5/optics/Guide_channeled.comp"
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
#line 13609 "ISIS_MERLIN.c"
}   /* End of between_guide=Guide_channeled() SETTING parameter declarations. */
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
  mcabsorbCompbetween_guide:
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

  /* TRACE Component fermi [7] */
  mccoordschange(mcposrfermi, mcrotrfermi,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component fermi (without coords transformations) */
  mcJumpTrace_fermi:
  SIG_MESSAGE("fermi (Trace)");
  mcDEBUG_COMP("fermi")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompfermi
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
#define mccompcurname  fermi
#define mccompcurtype  FermiChopper
#define mccompcurindex 7
#define FCVars mccfermi_FCVars
{   /* Declarations of fermi=FermiChopper() SETTING parameters. */
MCNUM phase = mccfermi_phase;
MCNUM radius = mccfermi_radius;
MCNUM nu = mccfermi_nu;
MCNUM w = mccfermi_w;
MCNUM nslit = mccfermi_nslit;
MCNUM R0 = mccfermi_R0;
MCNUM Qc = mccfermi_Qc;
MCNUM alpha = mccfermi_alpha;
MCNUM m = mccfermi_m;
MCNUM W = mccfermi_W;
MCNUM length = mccfermi_length;
MCNUM eff = mccfermi_eff;
MCNUM zero_time = mccfermi_zero_time;
MCNUM xwidth = mccfermi_xwidth;
MCNUM verbose = mccfermi_verbose;
MCNUM yheight = mccfermi_yheight;
MCNUM curvature = mccfermi_curvature;
MCNUM delay = mccfermi_delay;
#line 361 "/usr/share/mcstas/2.5/optics/FermiChopper.comp"
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
#line 14165 "ISIS_MERLIN.c"
}   /* End of fermi=FermiChopper() SETTING parameter declarations. */
#undef FCVars
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompfermi:
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

  /* TRACE Component TOF_fermi [8] */
  mccoordschange(mcposrTOF_fermi, mcrotrTOF_fermi,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component TOF_fermi (without coords transformations) */
  mcJumpTrace_TOF_fermi:
  SIG_MESSAGE("TOF_fermi (Trace)");
  mcDEBUG_COMP("TOF_fermi")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompTOF_fermi
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
#define mccompcurname  TOF_fermi
#define mccompcurtype  TOF_monitor
#define mccompcurindex 8
#define nt mccTOF_fermi_nt
#define TOF_N mccTOF_fermi_TOF_N
#define TOF_p mccTOF_fermi_TOF_p
#define TOF_p2 mccTOF_fermi_TOF_p2
#define t_min mccTOF_fermi_t_min
#define t_max mccTOF_fermi_t_max
#define delta_t mccTOF_fermi_delta_t
{   /* Declarations of TOF_fermi=TOF_monitor() SETTING parameters. */
char* filename = mccTOF_fermi_filename;
MCNUM xmin = mccTOF_fermi_xmin;
MCNUM xmax = mccTOF_fermi_xmax;
MCNUM ymin = mccTOF_fermi_ymin;
MCNUM ymax = mccTOF_fermi_ymax;
MCNUM xwidth = mccTOF_fermi_xwidth;
MCNUM yheight = mccTOF_fermi_yheight;
MCNUM tmin = mccTOF_fermi_tmin;
MCNUM tmax = mccTOF_fermi_tmax;
MCNUM dt = mccTOF_fermi_dt;
MCNUM restore_neutron = mccTOF_fermi_restore_neutron;
int nowritefile = mccTOF_fermi_nowritefile;
#line 95 "/usr/share/mcstas/2.5/monitors/TOF_monitor.comp"
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
#line 14310 "ISIS_MERLIN.c"
}   /* End of TOF_fermi=TOF_monitor() SETTING parameter declarations. */
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
  mcabsorbCompTOF_fermi:
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

  /* TRACE Component final_Guide [9] */
  mccoordschange(mcposrfinal_Guide, mcrotrfinal_Guide,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component final_Guide (without coords transformations) */
  mcJumpTrace_final_Guide:
  SIG_MESSAGE("final_Guide (Trace)");
  mcDEBUG_COMP("final_Guide")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompfinal_Guide
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
#define mccompcurname  final_Guide
#define mccompcurtype  Guide_channeled
#define mccompcurindex 9
#define w1c mccfinal_Guide_w1c
#define w2c mccfinal_Guide_w2c
#define ww mccfinal_Guide_ww
#define hh mccfinal_Guide_hh
#define whalf mccfinal_Guide_whalf
#define hhalf mccfinal_Guide_hhalf
#define lwhalf mccfinal_Guide_lwhalf
#define lhhalf mccfinal_Guide_lhhalf
{   /* Declarations of final_Guide=Guide_channeled() SETTING parameters. */
MCNUM w1 = mccfinal_Guide_w1;
MCNUM h1 = mccfinal_Guide_h1;
MCNUM w2 = mccfinal_Guide_w2;
MCNUM h2 = mccfinal_Guide_h2;
MCNUM l = mccfinal_Guide_l;
MCNUM R0 = mccfinal_Guide_R0;
MCNUM Qc = mccfinal_Guide_Qc;
MCNUM alpha = mccfinal_Guide_alpha;
MCNUM m = mccfinal_Guide_m;
MCNUM nslit = mccfinal_Guide_nslit;
MCNUM d = mccfinal_Guide_d;
MCNUM Qcx = mccfinal_Guide_Qcx;
MCNUM Qcy = mccfinal_Guide_Qcy;
MCNUM alphax = mccfinal_Guide_alphax;
MCNUM alphay = mccfinal_Guide_alphay;
MCNUM W = mccfinal_Guide_W;
MCNUM mx = mccfinal_Guide_mx;
MCNUM my = mccfinal_Guide_my;
MCNUM nu = mccfinal_Guide_nu;
MCNUM phase = mccfinal_Guide_phase;
#line 130 "/usr/share/mcstas/2.5/optics/Guide_channeled.comp"
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
#line 14616 "ISIS_MERLIN.c"
}   /* End of final_Guide=Guide_channeled() SETTING parameter declarations. */
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
  mcabsorbCompfinal_Guide:
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

  /* TRACE Component horz_div [10] */
  mccoordschange(mcposrhorz_div, mcrotrhorz_div,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component horz_div (without coords transformations) */
  mcJumpTrace_horz_div:
  SIG_MESSAGE("horz_div (Trace)");
  mcDEBUG_COMP("horz_div")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbComphorz_div
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
#define mccompcurname  horz_div
#define mccompcurtype  Hdiv_monitor
#define mccompcurindex 10
#define nh mcchorz_div_nh
#define Div_N mcchorz_div_Div_N
#define Div_p mcchorz_div_Div_p
#define Div_p2 mcchorz_div_Div_p2
{   /* Declarations of horz_div=Hdiv_monitor() SETTING parameters. */
char* filename = mcchorz_div_filename;
MCNUM xmin = mcchorz_div_xmin;
MCNUM xmax = mcchorz_div_xmax;
MCNUM ymin = mcchorz_div_ymin;
MCNUM ymax = mcchorz_div_ymax;
MCNUM xwidth = mcchorz_div_xwidth;
MCNUM yheight = mcchorz_div_yheight;
MCNUM h_maxdiv = mcchorz_div_h_maxdiv;
MCNUM restore_neutron = mcchorz_div_restore_neutron;
int nowritefile = mcchorz_div_nowritefile;
#line 85 "/usr/share/mcstas/2.5/monitors/Hdiv_monitor.comp"
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
#line 14766 "ISIS_MERLIN.c"
}   /* End of horz_div=Hdiv_monitor() SETTING parameter declarations. */
#undef Div_p2
#undef Div_p
#undef Div_N
#undef nh
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbComphorz_div:
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

  /* TRACE Component PSD [11] */
  mccoordschange(mcposrPSD, mcrotrPSD,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component PSD (without coords transformations) */
  mcJumpTrace_PSD:
  SIG_MESSAGE("PSD (Trace)");
  mcDEBUG_COMP("PSD")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompPSD
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
#define mccompcurname  PSD
#define mccompcurtype  PSD_monitor
#define mccompcurindex 11
#define nx mccPSD_nx
#define ny mccPSD_ny
#define PSD_N mccPSD_PSD_N
#define PSD_p mccPSD_PSD_p
#define PSD_p2 mccPSD_PSD_p2
{   /* Declarations of PSD=PSD_monitor() SETTING parameters. */
char* filename = mccPSD_filename;
MCNUM xmin = mccPSD_xmin;
MCNUM xmax = mccPSD_xmax;
MCNUM ymin = mccPSD_ymin;
MCNUM ymax = mccPSD_ymax;
MCNUM xwidth = mccPSD_xwidth;
MCNUM yheight = mccPSD_yheight;
MCNUM restore_neutron = mccPSD_restore_neutron;
int nowritefile = mccPSD_nowritefile;
#line 83 "/usr/share/mcstas/2.5/monitors/PSD_monitor.comp"
{
    int i,j;

    PROP_Z0;
    if (x>xmin && x<xmax && y>ymin && y<ymax)
    {
      i = floor((x - xmin)*nx/(xmax - xmin));
      j = floor((y - ymin)*ny/(ymax - ymin));
      PSD_N[i][j]++;
      PSD_p[i][j] += p;
      PSD_p2[i][j] += p*p;
      SCATTER;
    }
    if (restore_neutron) {
      RESTORE_NEUTRON(INDEX_CURRENT_COMP, x, y, z, vx, vy, vz, t, sx, sy, sz, p);
    }
}
#line 14908 "ISIS_MERLIN.c"
}   /* End of PSD=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef ny
#undef nx
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompPSD:
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

  /* TRACE Component det_lmon [12] */
  mccoordschange(mcposrdet_lmon, mcrotrdet_lmon,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component det_lmon (without coords transformations) */
  mcJumpTrace_det_lmon:
  SIG_MESSAGE("det_lmon (Trace)");
  mcDEBUG_COMP("det_lmon")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompdet_lmon
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
#define mccompcurname  det_lmon
#define mccompcurtype  L_monitor
#define mccompcurindex 12
#define nL mccdet_lmon_nL
#define L_N mccdet_lmon_L_N
#define L_p mccdet_lmon_L_p
#define L_p2 mccdet_lmon_L_p2
{   /* Declarations of det_lmon=L_monitor() SETTING parameters. */
char* filename = mccdet_lmon_filename;
MCNUM xmin = mccdet_lmon_xmin;
MCNUM xmax = mccdet_lmon_xmax;
MCNUM ymin = mccdet_lmon_ymin;
MCNUM ymax = mccdet_lmon_ymax;
MCNUM xwidth = mccdet_lmon_xwidth;
MCNUM yheight = mccdet_lmon_yheight;
MCNUM Lmin = mccdet_lmon_Lmin;
MCNUM Lmax = mccdet_lmon_Lmax;
MCNUM restore_neutron = mccdet_lmon_restore_neutron;
int nowritefile = mccdet_lmon_nowritefile;
#line 84 "/usr/share/mcstas/2.5/monitors/L_monitor.comp"
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
#line 15056 "ISIS_MERLIN.c"
}   /* End of det_lmon=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompdet_lmon:
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

  /* TRACE Component sample [13] */
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
#define mccompcurname  sample
#define mccompcurtype  Incoherent
#define mccompcurindex 13
#define VarsInc mccsample_VarsInc
#define offdata mccsample_offdata
{   /* Declarations of sample=Incoherent() SETTING parameters. */
char* geometry = mccsample_geometry;
MCNUM radius = mccsample_radius;
MCNUM xwidth = mccsample_xwidth;
MCNUM yheight = mccsample_yheight;
MCNUM zdepth = mccsample_zdepth;
MCNUM thickness = mccsample_thickness;
MCNUM target_x = mccsample_target_x;
MCNUM target_y = mccsample_target_y;
MCNUM target_z = mccsample_target_z;
MCNUM focus_r = mccsample_focus_r;
MCNUM focus_xw = mccsample_focus_xw;
MCNUM focus_yh = mccsample_focus_yh;
MCNUM focus_aw = mccsample_focus_aw;
MCNUM focus_ah = mccsample_focus_ah;
int target_index = mccsample_target_index;
MCNUM pack = mccsample_pack;
MCNUM p_interact = mccsample_p_interact;
MCNUM f_QE = mccsample_f_QE;
MCNUM gamma = mccsample_gamma;
MCNUM Etrans = mccsample_Etrans;
MCNUM deltaE = mccsample_deltaE;
MCNUM sigma_abs = mccsample_sigma_abs;
MCNUM sigma_inc = mccsample_sigma_inc;
MCNUM Vc = mccsample_Vc;
MCNUM concentric = mccsample_concentric;
MCNUM order = mccsample_order;
#line 231 "/usr/share/mcstas/2.5/samples/Incoherent.comp"
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
#line 15400 "ISIS_MERLIN.c"
}   /* End of sample=Incoherent() SETTING parameter declarations. */
#undef offdata
#undef VarsInc
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompsample:
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

  /* TRACE Component detectorarm [14] */
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
#define mccompcurname  detectorarm
#define mccompcurtype  Arm
#define mccompcurindex 14
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompdetectorarm:
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

  /* TRACE Component detector [15] */
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
#define mccompcurname  detector
#define mccompcurtype  TOF_monitor
#define mccompcurindex 15
#define nt mccdetector_nt
#define TOF_N mccdetector_TOF_N
#define TOF_p mccdetector_TOF_p
#define TOF_p2 mccdetector_TOF_p2
#define t_min mccdetector_t_min
#define t_max mccdetector_t_max
#define delta_t mccdetector_delta_t
{   /* Declarations of detector=TOF_monitor() SETTING parameters. */
char* filename = mccdetector_filename;
MCNUM xmin = mccdetector_xmin;
MCNUM xmax = mccdetector_xmax;
MCNUM ymin = mccdetector_ymin;
MCNUM ymax = mccdetector_ymax;
MCNUM xwidth = mccdetector_xwidth;
MCNUM yheight = mccdetector_yheight;
MCNUM tmin = mccdetector_tmin;
MCNUM tmax = mccdetector_tmax;
MCNUM dt = mccdetector_dt;
MCNUM restore_neutron = mccdetector_restore_neutron;
int nowritefile = mccdetector_nowritefile;
#line 95 "/usr/share/mcstas/2.5/monitors/TOF_monitor.comp"
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
#line 15649 "ISIS_MERLIN.c"
}   /* End of detector=TOF_monitor() SETTING parameter declarations. */
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
  mcabsorbCompdetector:
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

  /* User SAVE code for component 'TOF_mod'. */
  SIG_MESSAGE("TOF_mod (Save)");
#define mccompcurname  TOF_mod
#define mccompcurtype  TOF_monitor
#define mccompcurindex 3
#define nt mccTOF_mod_nt
#define TOF_N mccTOF_mod_TOF_N
#define TOF_p mccTOF_mod_TOF_p
#define TOF_p2 mccTOF_mod_TOF_p2
#define t_min mccTOF_mod_t_min
#define t_max mccTOF_mod_t_max
#define delta_t mccTOF_mod_delta_t
{   /* Declarations of TOF_mod=TOF_monitor() SETTING parameters. */
char* filename = mccTOF_mod_filename;
MCNUM xmin = mccTOF_mod_xmin;
MCNUM xmax = mccTOF_mod_xmax;
MCNUM ymin = mccTOF_mod_ymin;
MCNUM ymax = mccTOF_mod_ymax;
MCNUM xwidth = mccTOF_mod_xwidth;
MCNUM yheight = mccTOF_mod_yheight;
MCNUM tmin = mccTOF_mod_tmin;
MCNUM tmax = mccTOF_mod_tmax;
MCNUM dt = mccTOF_mod_dt;
MCNUM restore_neutron = mccTOF_mod_restore_neutron;
int nowritefile = mccTOF_mod_nowritefile;
#line 115 "/usr/share/mcstas/2.5/monitors/TOF_monitor.comp"
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
#line 15771 "ISIS_MERLIN.c"
}   /* End of TOF_mod=TOF_monitor() SETTING parameter declarations. */
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

  /* User SAVE code for component 'fermi'. */
  SIG_MESSAGE("fermi (Save)");
#define mccompcurname  fermi
#define mccompcurtype  FermiChopper
#define mccompcurindex 7
#define FCVars mccfermi_FCVars
{   /* Declarations of fermi=FermiChopper() SETTING parameters. */
MCNUM phase = mccfermi_phase;
MCNUM radius = mccfermi_radius;
MCNUM nu = mccfermi_nu;
MCNUM w = mccfermi_w;
MCNUM nslit = mccfermi_nslit;
MCNUM R0 = mccfermi_R0;
MCNUM Qc = mccfermi_Qc;
MCNUM alpha = mccfermi_alpha;
MCNUM m = mccfermi_m;
MCNUM W = mccfermi_W;
MCNUM length = mccfermi_length;
MCNUM eff = mccfermi_eff;
MCNUM zero_time = mccfermi_zero_time;
MCNUM xwidth = mccfermi_xwidth;
MCNUM verbose = mccfermi_verbose;
MCNUM yheight = mccfermi_yheight;
MCNUM curvature = mccfermi_curvature;
MCNUM delay = mccfermi_delay;
#line 785 "/usr/share/mcstas/2.5/optics/FermiChopper.comp"
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
#line 15886 "ISIS_MERLIN.c"
}   /* End of fermi=FermiChopper() SETTING parameter declarations. */
#undef FCVars
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'TOF_fermi'. */
  SIG_MESSAGE("TOF_fermi (Save)");
#define mccompcurname  TOF_fermi
#define mccompcurtype  TOF_monitor
#define mccompcurindex 8
#define nt mccTOF_fermi_nt
#define TOF_N mccTOF_fermi_TOF_N
#define TOF_p mccTOF_fermi_TOF_p
#define TOF_p2 mccTOF_fermi_TOF_p2
#define t_min mccTOF_fermi_t_min
#define t_max mccTOF_fermi_t_max
#define delta_t mccTOF_fermi_delta_t
{   /* Declarations of TOF_fermi=TOF_monitor() SETTING parameters. */
char* filename = mccTOF_fermi_filename;
MCNUM xmin = mccTOF_fermi_xmin;
MCNUM xmax = mccTOF_fermi_xmax;
MCNUM ymin = mccTOF_fermi_ymin;
MCNUM ymax = mccTOF_fermi_ymax;
MCNUM xwidth = mccTOF_fermi_xwidth;
MCNUM yheight = mccTOF_fermi_yheight;
MCNUM tmin = mccTOF_fermi_tmin;
MCNUM tmax = mccTOF_fermi_tmax;
MCNUM dt = mccTOF_fermi_dt;
MCNUM restore_neutron = mccTOF_fermi_restore_neutron;
int nowritefile = mccTOF_fermi_nowritefile;
#line 115 "/usr/share/mcstas/2.5/monitors/TOF_monitor.comp"
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
#line 15930 "ISIS_MERLIN.c"
}   /* End of TOF_fermi=TOF_monitor() SETTING parameter declarations. */
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

  /* User SAVE code for component 'horz_div'. */
  SIG_MESSAGE("horz_div (Save)");
#define mccompcurname  horz_div
#define mccompcurtype  Hdiv_monitor
#define mccompcurindex 10
#define nh mcchorz_div_nh
#define Div_N mcchorz_div_Div_N
#define Div_p mcchorz_div_Div_p
#define Div_p2 mcchorz_div_Div_p2
{   /* Declarations of horz_div=Hdiv_monitor() SETTING parameters. */
char* filename = mcchorz_div_filename;
MCNUM xmin = mcchorz_div_xmin;
MCNUM xmax = mcchorz_div_xmax;
MCNUM ymin = mcchorz_div_ymin;
MCNUM ymax = mcchorz_div_ymax;
MCNUM xwidth = mcchorz_div_xwidth;
MCNUM yheight = mcchorz_div_yheight;
MCNUM h_maxdiv = mcchorz_div_h_maxdiv;
MCNUM restore_neutron = mcchorz_div_restore_neutron;
int nowritefile = mcchorz_div_nowritefile;
#line 107 "/usr/share/mcstas/2.5/monitors/Hdiv_monitor.comp"
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
#line 15975 "ISIS_MERLIN.c"
}   /* End of horz_div=Hdiv_monitor() SETTING parameter declarations. */
#undef Div_p2
#undef Div_p
#undef Div_N
#undef nh
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'PSD'. */
  SIG_MESSAGE("PSD (Save)");
#define mccompcurname  PSD
#define mccompcurtype  PSD_monitor
#define mccompcurindex 11
#define nx mccPSD_nx
#define ny mccPSD_ny
#define PSD_N mccPSD_PSD_N
#define PSD_p mccPSD_PSD_p
#define PSD_p2 mccPSD_PSD_p2
{   /* Declarations of PSD=PSD_monitor() SETTING parameters. */
char* filename = mccPSD_filename;
MCNUM xmin = mccPSD_xmin;
MCNUM xmax = mccPSD_xmax;
MCNUM ymin = mccPSD_ymin;
MCNUM ymax = mccPSD_ymax;
MCNUM xwidth = mccPSD_xwidth;
MCNUM yheight = mccPSD_yheight;
MCNUM restore_neutron = mccPSD_restore_neutron;
int nowritefile = mccPSD_nowritefile;
#line 101 "/usr/share/mcstas/2.5/monitors/PSD_monitor.comp"
{
    if (!nowritefile) {
      DETECTOR_OUT_2D(
          "PSD monitor",
          "X position [cm]",
          "Y position [cm]",
          xmin*100.0, xmax*100.0, ymin*100.0, ymax*100.0,
          nx, ny,
          &PSD_N[0][0],&PSD_p[0][0],&PSD_p2[0][0],
          filename);
    }
}
#line 16018 "ISIS_MERLIN.c"
}   /* End of PSD=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef ny
#undef nx
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'det_lmon'. */
  SIG_MESSAGE("det_lmon (Save)");
#define mccompcurname  det_lmon
#define mccompcurtype  L_monitor
#define mccompcurindex 12
#define nL mccdet_lmon_nL
#define L_N mccdet_lmon_L_N
#define L_p mccdet_lmon_L_p
#define L_p2 mccdet_lmon_L_p2
{   /* Declarations of det_lmon=L_monitor() SETTING parameters. */
char* filename = mccdet_lmon_filename;
MCNUM xmin = mccdet_lmon_xmin;
MCNUM xmax = mccdet_lmon_xmax;
MCNUM ymin = mccdet_lmon_ymin;
MCNUM ymax = mccdet_lmon_ymax;
MCNUM xwidth = mccdet_lmon_xwidth;
MCNUM yheight = mccdet_lmon_yheight;
MCNUM Lmin = mccdet_lmon_Lmin;
MCNUM Lmax = mccdet_lmon_Lmax;
MCNUM restore_neutron = mccdet_lmon_restore_neutron;
int nowritefile = mccdet_lmon_nowritefile;
#line 107 "/usr/share/mcstas/2.5/monitors/L_monitor.comp"
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
#line 16062 "ISIS_MERLIN.c"
}   /* End of det_lmon=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'detector'. */
  SIG_MESSAGE("detector (Save)");
#define mccompcurname  detector
#define mccompcurtype  TOF_monitor
#define mccompcurindex 15
#define nt mccdetector_nt
#define TOF_N mccdetector_TOF_N
#define TOF_p mccdetector_TOF_p
#define TOF_p2 mccdetector_TOF_p2
#define t_min mccdetector_t_min
#define t_max mccdetector_t_max
#define delta_t mccdetector_delta_t
{   /* Declarations of detector=TOF_monitor() SETTING parameters. */
char* filename = mccdetector_filename;
MCNUM xmin = mccdetector_xmin;
MCNUM xmax = mccdetector_xmax;
MCNUM ymin = mccdetector_ymin;
MCNUM ymax = mccdetector_ymax;
MCNUM xwidth = mccdetector_xwidth;
MCNUM yheight = mccdetector_yheight;
MCNUM tmin = mccdetector_tmin;
MCNUM tmax = mccdetector_tmax;
MCNUM dt = mccdetector_dt;
MCNUM restore_neutron = mccdetector_restore_neutron;
int nowritefile = mccdetector_nowritefile;
#line 115 "/usr/share/mcstas/2.5/monitors/TOF_monitor.comp"
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
#line 16109 "ISIS_MERLIN.c"
}   /* End of detector=TOF_monitor() SETTING parameter declarations. */
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

  if (!handle) mcsiminfo_close(); 
} /* end save */
void mcfinally(void) {
  /* User component FINALLY code. */
  mcsiminfo_init(NULL);
  mcsave(mcsiminfo_file); /* save data when simulation ends */

    if (!mcNCounter[1]) fprintf(stderr, "Warning: No neutron could reach Component[1] Origin\n");
    if (mcAbsorbProp[1]) fprintf(stderr, "Warning: %g events were removed in Component[1] Origin=Arm()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[1]);
    if (!mcNCounter[2]) fprintf(stderr, "Warning: No neutron could reach Component[2] isis_mod\n");
    if (mcAbsorbProp[2]) fprintf(stderr, "Warning: %g events were removed in Component[2] isis_mod=ViewModISIS()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[2]);
    if (!mcNCounter[3]) fprintf(stderr, "Warning: No neutron could reach Component[3] TOF_mod\n");
    if (mcAbsorbProp[3]) fprintf(stderr, "Warning: %g events were removed in Component[3] TOF_mod=TOF_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[3]);
    if (!mcNCounter[4]) fprintf(stderr, "Warning: No neutron could reach Component[4] shutter_guide\n");
    if (mcAbsorbProp[4]) fprintf(stderr, "Warning: %g events were removed in Component[4] shutter_guide=Guide_channeled()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[4]);
    if (!mcNCounter[5]) fprintf(stderr, "Warning: No neutron could reach Component[5] conv_guide\n");
    if (mcAbsorbProp[5]) fprintf(stderr, "Warning: %g events were removed in Component[5] conv_guide=Guide_channeled()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[5]);
    if (!mcNCounter[6]) fprintf(stderr, "Warning: No neutron could reach Component[6] between_guide\n");
    if (mcAbsorbProp[6]) fprintf(stderr, "Warning: %g events were removed in Component[6] between_guide=Guide_channeled()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[6]);
    if (!mcNCounter[7]) fprintf(stderr, "Warning: No neutron could reach Component[7] fermi\n");
    if (mcAbsorbProp[7]) fprintf(stderr, "Warning: %g events were removed in Component[7] fermi=FermiChopper()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[7]);
    if (!mcNCounter[8]) fprintf(stderr, "Warning: No neutron could reach Component[8] TOF_fermi\n");
    if (mcAbsorbProp[8]) fprintf(stderr, "Warning: %g events were removed in Component[8] TOF_fermi=TOF_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[8]);
    if (!mcNCounter[9]) fprintf(stderr, "Warning: No neutron could reach Component[9] final_Guide\n");
    if (mcAbsorbProp[9]) fprintf(stderr, "Warning: %g events were removed in Component[9] final_Guide=Guide_channeled()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[9]);
    if (!mcNCounter[10]) fprintf(stderr, "Warning: No neutron could reach Component[10] horz_div\n");
    if (mcAbsorbProp[10]) fprintf(stderr, "Warning: %g events were removed in Component[10] horz_div=Hdiv_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[10]);
    if (!mcNCounter[11]) fprintf(stderr, "Warning: No neutron could reach Component[11] PSD\n");
    if (mcAbsorbProp[11]) fprintf(stderr, "Warning: %g events were removed in Component[11] PSD=PSD_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[11]);
    if (!mcNCounter[12]) fprintf(stderr, "Warning: No neutron could reach Component[12] det_lmon\n");
    if (mcAbsorbProp[12]) fprintf(stderr, "Warning: %g events were removed in Component[12] det_lmon=L_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[12]);
    if (!mcNCounter[13]) fprintf(stderr, "Warning: No neutron could reach Component[13] sample\n");
    if (mcAbsorbProp[13]) fprintf(stderr, "Warning: %g events were removed in Component[13] sample=Incoherent()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[13]);
    if (!mcNCounter[14]) fprintf(stderr, "Warning: No neutron could reach Component[14] detectorarm\n");
    if (mcAbsorbProp[14]) fprintf(stderr, "Warning: %g events were removed in Component[14] detectorarm=Arm()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[14]);
    if (!mcNCounter[15]) fprintf(stderr, "Warning: No neutron could reach Component[15] detector\n");
    if (mcAbsorbProp[15]) fprintf(stderr, "Warning: %g events were removed in Component[15] detector=TOF_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[15]);
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
#line 40 "/usr/share/mcstas/2.5/optics/Arm.comp"
{
  /* A bit ugly; hard-coded dimensions. */
  
  line(0,0,0,0.2,0,0);
  line(0,0,0,0,0.2,0);
  line(0,0,0,0,0,0.2);
}
#line 16188 "ISIS_MERLIN.c"
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'isis_mod'. */
  SIG_MESSAGE("isis_mod (McDisplay)");
  printf("MCDISPLAY: component %s\n", "isis_mod");
#define mccompcurname  isis_mod
#define mccompcurtype  ViewModISIS
#define mccompcurindex 2
{   /* Declarations of isis_mod=ViewModISIS() SETTING parameters. */
char* Face = mccisis_mod_Face;
MCNUM E0 = mccisis_mod_E0;
MCNUM E1 = mccisis_mod_E1;
MCNUM modPosition = mccisis_mod_modPosition;
MCNUM xw = mccisis_mod_xw;
MCNUM yh = mccisis_mod_yh;
MCNUM focus_xw = mccisis_mod_focus_xw;
MCNUM focus_yh = mccisis_mod_focus_yh;
MCNUM dist = mccisis_mod_dist;
int verbose = mccisis_mod_verbose;
MCNUM beamcurrent = mccisis_mod_beamcurrent;
#line 1137 "/usr/share/mcstas/2.5/contrib/ViewModISIS.comp"
{
  double cirp=0.0,cirq=0.3,pi=3.141592654;
  int pp=0; /* circle drawing parameter*/



  magnify("xy");
  multiline(5,-0.5*TS.XAxis,-0.5*TS.ZAxis,0.0,
	    0.5*TS.XAxis,-0.5*TS.ZAxis,0.0,
	    0.5*TS.XAxis,0.5*TS.ZAxis,0.0,
	    -0.5*TS.XAxis,0.5*TS.ZAxis,0.0,
	    -0.5*TS.XAxis,-0.5*TS.ZAxis,0.0);
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
#line 16238 "ISIS_MERLIN.c"
}   /* End of isis_mod=ViewModISIS() SETTING parameter declarations. */
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'TOF_mod'. */
  SIG_MESSAGE("TOF_mod (McDisplay)");
  printf("MCDISPLAY: component %s\n", "TOF_mod");
#define mccompcurname  TOF_mod
#define mccompcurtype  TOF_monitor
#define mccompcurindex 3
#define nt mccTOF_mod_nt
#define TOF_N mccTOF_mod_TOF_N
#define TOF_p mccTOF_mod_TOF_p
#define TOF_p2 mccTOF_mod_TOF_p2
#define t_min mccTOF_mod_t_min
#define t_max mccTOF_mod_t_max
#define delta_t mccTOF_mod_delta_t
{   /* Declarations of TOF_mod=TOF_monitor() SETTING parameters. */
char* filename = mccTOF_mod_filename;
MCNUM xmin = mccTOF_mod_xmin;
MCNUM xmax = mccTOF_mod_xmax;
MCNUM ymin = mccTOF_mod_ymin;
MCNUM ymax = mccTOF_mod_ymax;
MCNUM xwidth = mccTOF_mod_xwidth;
MCNUM yheight = mccTOF_mod_yheight;
MCNUM tmin = mccTOF_mod_tmin;
MCNUM tmax = mccTOF_mod_tmax;
MCNUM dt = mccTOF_mod_dt;
MCNUM restore_neutron = mccTOF_mod_restore_neutron;
int nowritefile = mccTOF_mod_nowritefile;
#line 128 "/usr/share/mcstas/2.5/monitors/TOF_monitor.comp"
{
  
  multiline(5, (double)xmin, (double)ymin, 0.0,
               (double)xmax, (double)ymin, 0.0,
               (double)xmax, (double)ymax, 0.0,
               (double)xmin, (double)ymax, 0.0,
               (double)xmin, (double)ymin, 0.0);
}
#line 16279 "ISIS_MERLIN.c"
}   /* End of TOF_mod=TOF_monitor() SETTING parameter declarations. */
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

  /* MCDISPLAY code for component 'shutter_guide'. */
  SIG_MESSAGE("shutter_guide (McDisplay)");
  printf("MCDISPLAY: component %s\n", "shutter_guide");
#define mccompcurname  shutter_guide
#define mccompcurtype  Guide_channeled
#define mccompcurindex 4
#define w1c mccshutter_guide_w1c
#define w2c mccshutter_guide_w2c
#define ww mccshutter_guide_ww
#define hh mccshutter_guide_hh
#define whalf mccshutter_guide_whalf
#define hhalf mccshutter_guide_hhalf
#define lwhalf mccshutter_guide_lwhalf
#define lhhalf mccshutter_guide_lhhalf
{   /* Declarations of shutter_guide=Guide_channeled() SETTING parameters. */
MCNUM w1 = mccshutter_guide_w1;
MCNUM h1 = mccshutter_guide_h1;
MCNUM w2 = mccshutter_guide_w2;
MCNUM h2 = mccshutter_guide_h2;
MCNUM l = mccshutter_guide_l;
MCNUM R0 = mccshutter_guide_R0;
MCNUM Qc = mccshutter_guide_Qc;
MCNUM alpha = mccshutter_guide_alpha;
MCNUM m = mccshutter_guide_m;
MCNUM nslit = mccshutter_guide_nslit;
MCNUM d = mccshutter_guide_d;
MCNUM Qcx = mccshutter_guide_Qcx;
MCNUM Qcy = mccshutter_guide_Qcy;
MCNUM alphax = mccshutter_guide_alphax;
MCNUM alphay = mccshutter_guide_alphay;
MCNUM W = mccshutter_guide_W;
MCNUM mx = mccshutter_guide_mx;
MCNUM my = mccshutter_guide_my;
MCNUM nu = mccshutter_guide_nu;
MCNUM phase = mccshutter_guide_phase;
#line 296 "/usr/share/mcstas/2.5/optics/Guide_channeled.comp"
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
#line 16358 "ISIS_MERLIN.c"
}   /* End of shutter_guide=Guide_channeled() SETTING parameter declarations. */
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

  /* MCDISPLAY code for component 'conv_guide'. */
  SIG_MESSAGE("conv_guide (McDisplay)");
  printf("MCDISPLAY: component %s\n", "conv_guide");
#define mccompcurname  conv_guide
#define mccompcurtype  Guide_channeled
#define mccompcurindex 5
#define w1c mccconv_guide_w1c
#define w2c mccconv_guide_w2c
#define ww mccconv_guide_ww
#define hh mccconv_guide_hh
#define whalf mccconv_guide_whalf
#define hhalf mccconv_guide_hhalf
#define lwhalf mccconv_guide_lwhalf
#define lhhalf mccconv_guide_lhhalf
{   /* Declarations of conv_guide=Guide_channeled() SETTING parameters. */
MCNUM w1 = mccconv_guide_w1;
MCNUM h1 = mccconv_guide_h1;
MCNUM w2 = mccconv_guide_w2;
MCNUM h2 = mccconv_guide_h2;
MCNUM l = mccconv_guide_l;
MCNUM R0 = mccconv_guide_R0;
MCNUM Qc = mccconv_guide_Qc;
MCNUM alpha = mccconv_guide_alpha;
MCNUM m = mccconv_guide_m;
MCNUM nslit = mccconv_guide_nslit;
MCNUM d = mccconv_guide_d;
MCNUM Qcx = mccconv_guide_Qcx;
MCNUM Qcy = mccconv_guide_Qcy;
MCNUM alphax = mccconv_guide_alphax;
MCNUM alphay = mccconv_guide_alphay;
MCNUM W = mccconv_guide_W;
MCNUM mx = mccconv_guide_mx;
MCNUM my = mccconv_guide_my;
MCNUM nu = mccconv_guide_nu;
MCNUM phase = mccconv_guide_phase;
#line 296 "/usr/share/mcstas/2.5/optics/Guide_channeled.comp"
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
#line 16438 "ISIS_MERLIN.c"
}   /* End of conv_guide=Guide_channeled() SETTING parameter declarations. */
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

  /* MCDISPLAY code for component 'between_guide'. */
  SIG_MESSAGE("between_guide (McDisplay)");
  printf("MCDISPLAY: component %s\n", "between_guide");
#define mccompcurname  between_guide
#define mccompcurtype  Guide_channeled
#define mccompcurindex 6
#define w1c mccbetween_guide_w1c
#define w2c mccbetween_guide_w2c
#define ww mccbetween_guide_ww
#define hh mccbetween_guide_hh
#define whalf mccbetween_guide_whalf
#define hhalf mccbetween_guide_hhalf
#define lwhalf mccbetween_guide_lwhalf
#define lhhalf mccbetween_guide_lhhalf
{   /* Declarations of between_guide=Guide_channeled() SETTING parameters. */
MCNUM w1 = mccbetween_guide_w1;
MCNUM h1 = mccbetween_guide_h1;
MCNUM w2 = mccbetween_guide_w2;
MCNUM h2 = mccbetween_guide_h2;
MCNUM l = mccbetween_guide_l;
MCNUM R0 = mccbetween_guide_R0;
MCNUM Qc = mccbetween_guide_Qc;
MCNUM alpha = mccbetween_guide_alpha;
MCNUM m = mccbetween_guide_m;
MCNUM nslit = mccbetween_guide_nslit;
MCNUM d = mccbetween_guide_d;
MCNUM Qcx = mccbetween_guide_Qcx;
MCNUM Qcy = mccbetween_guide_Qcy;
MCNUM alphax = mccbetween_guide_alphax;
MCNUM alphay = mccbetween_guide_alphay;
MCNUM W = mccbetween_guide_W;
MCNUM mx = mccbetween_guide_mx;
MCNUM my = mccbetween_guide_my;
MCNUM nu = mccbetween_guide_nu;
MCNUM phase = mccbetween_guide_phase;
#line 296 "/usr/share/mcstas/2.5/optics/Guide_channeled.comp"
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
#line 16518 "ISIS_MERLIN.c"
}   /* End of between_guide=Guide_channeled() SETTING parameter declarations. */
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

  /* MCDISPLAY code for component 'fermi'. */
  SIG_MESSAGE("fermi (McDisplay)");
  printf("MCDISPLAY: component %s\n", "fermi");
#define mccompcurname  fermi
#define mccompcurtype  FermiChopper
#define mccompcurindex 7
#define FCVars mccfermi_FCVars
{   /* Declarations of fermi=FermiChopper() SETTING parameters. */
MCNUM phase = mccfermi_phase;
MCNUM radius = mccfermi_radius;
MCNUM nu = mccfermi_nu;
MCNUM w = mccfermi_w;
MCNUM nslit = mccfermi_nslit;
MCNUM R0 = mccfermi_R0;
MCNUM Qc = mccfermi_Qc;
MCNUM alpha = mccfermi_alpha;
MCNUM m = mccfermi_m;
MCNUM W = mccfermi_W;
MCNUM length = mccfermi_length;
MCNUM eff = mccfermi_eff;
MCNUM zero_time = mccfermi_zero_time;
MCNUM xwidth = mccfermi_xwidth;
MCNUM verbose = mccfermi_verbose;
MCNUM yheight = mccfermi_yheight;
MCNUM curvature = mccfermi_curvature;
MCNUM delay = mccfermi_delay;
#line 863 "/usr/share/mcstas/2.5/optics/FermiChopper.comp"
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
#line 16627 "ISIS_MERLIN.c"
}   /* End of fermi=FermiChopper() SETTING parameter declarations. */
#undef FCVars
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'TOF_fermi'. */
  SIG_MESSAGE("TOF_fermi (McDisplay)");
  printf("MCDISPLAY: component %s\n", "TOF_fermi");
#define mccompcurname  TOF_fermi
#define mccompcurtype  TOF_monitor
#define mccompcurindex 8
#define nt mccTOF_fermi_nt
#define TOF_N mccTOF_fermi_TOF_N
#define TOF_p mccTOF_fermi_TOF_p
#define TOF_p2 mccTOF_fermi_TOF_p2
#define t_min mccTOF_fermi_t_min
#define t_max mccTOF_fermi_t_max
#define delta_t mccTOF_fermi_delta_t
{   /* Declarations of TOF_fermi=TOF_monitor() SETTING parameters. */
char* filename = mccTOF_fermi_filename;
MCNUM xmin = mccTOF_fermi_xmin;
MCNUM xmax = mccTOF_fermi_xmax;
MCNUM ymin = mccTOF_fermi_ymin;
MCNUM ymax = mccTOF_fermi_ymax;
MCNUM xwidth = mccTOF_fermi_xwidth;
MCNUM yheight = mccTOF_fermi_yheight;
MCNUM tmin = mccTOF_fermi_tmin;
MCNUM tmax = mccTOF_fermi_tmax;
MCNUM dt = mccTOF_fermi_dt;
MCNUM restore_neutron = mccTOF_fermi_restore_neutron;
int nowritefile = mccTOF_fermi_nowritefile;
#line 128 "/usr/share/mcstas/2.5/monitors/TOF_monitor.comp"
{
  
  multiline(5, (double)xmin, (double)ymin, 0.0,
               (double)xmax, (double)ymin, 0.0,
               (double)xmax, (double)ymax, 0.0,
               (double)xmin, (double)ymax, 0.0,
               (double)xmin, (double)ymin, 0.0);
}
#line 16669 "ISIS_MERLIN.c"
}   /* End of TOF_fermi=TOF_monitor() SETTING parameter declarations. */
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

  /* MCDISPLAY code for component 'final_Guide'. */
  SIG_MESSAGE("final_Guide (McDisplay)");
  printf("MCDISPLAY: component %s\n", "final_Guide");
#define mccompcurname  final_Guide
#define mccompcurtype  Guide_channeled
#define mccompcurindex 9
#define w1c mccfinal_Guide_w1c
#define w2c mccfinal_Guide_w2c
#define ww mccfinal_Guide_ww
#define hh mccfinal_Guide_hh
#define whalf mccfinal_Guide_whalf
#define hhalf mccfinal_Guide_hhalf
#define lwhalf mccfinal_Guide_lwhalf
#define lhhalf mccfinal_Guide_lhhalf
{   /* Declarations of final_Guide=Guide_channeled() SETTING parameters. */
MCNUM w1 = mccfinal_Guide_w1;
MCNUM h1 = mccfinal_Guide_h1;
MCNUM w2 = mccfinal_Guide_w2;
MCNUM h2 = mccfinal_Guide_h2;
MCNUM l = mccfinal_Guide_l;
MCNUM R0 = mccfinal_Guide_R0;
MCNUM Qc = mccfinal_Guide_Qc;
MCNUM alpha = mccfinal_Guide_alpha;
MCNUM m = mccfinal_Guide_m;
MCNUM nslit = mccfinal_Guide_nslit;
MCNUM d = mccfinal_Guide_d;
MCNUM Qcx = mccfinal_Guide_Qcx;
MCNUM Qcy = mccfinal_Guide_Qcy;
MCNUM alphax = mccfinal_Guide_alphax;
MCNUM alphay = mccfinal_Guide_alphay;
MCNUM W = mccfinal_Guide_W;
MCNUM mx = mccfinal_Guide_mx;
MCNUM my = mccfinal_Guide_my;
MCNUM nu = mccfinal_Guide_nu;
MCNUM phase = mccfinal_Guide_phase;
#line 296 "/usr/share/mcstas/2.5/optics/Guide_channeled.comp"
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
#line 16748 "ISIS_MERLIN.c"
}   /* End of final_Guide=Guide_channeled() SETTING parameter declarations. */
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

  /* MCDISPLAY code for component 'horz_div'. */
  SIG_MESSAGE("horz_div (McDisplay)");
  printf("MCDISPLAY: component %s\n", "horz_div");
#define mccompcurname  horz_div
#define mccompcurtype  Hdiv_monitor
#define mccompcurindex 10
#define nh mcchorz_div_nh
#define Div_N mcchorz_div_Div_N
#define Div_p mcchorz_div_Div_p
#define Div_p2 mcchorz_div_Div_p2
{   /* Declarations of horz_div=Hdiv_monitor() SETTING parameters. */
char* filename = mcchorz_div_filename;
MCNUM xmin = mcchorz_div_xmin;
MCNUM xmax = mcchorz_div_xmax;
MCNUM ymin = mcchorz_div_ymin;
MCNUM ymax = mcchorz_div_ymax;
MCNUM xwidth = mcchorz_div_xwidth;
MCNUM yheight = mcchorz_div_yheight;
MCNUM h_maxdiv = mcchorz_div_h_maxdiv;
MCNUM restore_neutron = mcchorz_div_restore_neutron;
int nowritefile = mcchorz_div_nowritefile;
#line 120 "/usr/share/mcstas/2.5/monitors/Hdiv_monitor.comp"
{
    
    multiline(5, (double)xmin, (double)ymin, 0.0,
                 (double)xmax, (double)ymin, 0.0,
                 (double)xmax, (double)ymax, 0.0,
                 (double)xmin, (double)ymax, 0.0,
                 (double)xmin, (double)ymin, 0.0);
}
#line 16792 "ISIS_MERLIN.c"
}   /* End of horz_div=Hdiv_monitor() SETTING parameter declarations. */
#undef Div_p2
#undef Div_p
#undef Div_N
#undef nh
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'PSD'. */
  SIG_MESSAGE("PSD (McDisplay)");
  printf("MCDISPLAY: component %s\n", "PSD");
#define mccompcurname  PSD
#define mccompcurtype  PSD_monitor
#define mccompcurindex 11
#define nx mccPSD_nx
#define ny mccPSD_ny
#define PSD_N mccPSD_PSD_N
#define PSD_p mccPSD_PSD_p
#define PSD_p2 mccPSD_PSD_p2
{   /* Declarations of PSD=PSD_monitor() SETTING parameters. */
char* filename = mccPSD_filename;
MCNUM xmin = mccPSD_xmin;
MCNUM xmax = mccPSD_xmax;
MCNUM ymin = mccPSD_ymin;
MCNUM ymax = mccPSD_ymax;
MCNUM xwidth = mccPSD_xwidth;
MCNUM yheight = mccPSD_yheight;
MCNUM restore_neutron = mccPSD_restore_neutron;
int nowritefile = mccPSD_nowritefile;
#line 115 "/usr/share/mcstas/2.5/monitors/PSD_monitor.comp"
{
  
  multiline(5, (double)xmin, (double)ymin, 0.0,
               (double)xmax, (double)ymin, 0.0,
               (double)xmax, (double)ymax, 0.0,
               (double)xmin, (double)ymax, 0.0,
               (double)xmin, (double)ymin, 0.0);
}
#line 16832 "ISIS_MERLIN.c"
}   /* End of PSD=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef ny
#undef nx
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'det_lmon'. */
  SIG_MESSAGE("det_lmon (McDisplay)");
  printf("MCDISPLAY: component %s\n", "det_lmon");
#define mccompcurname  det_lmon
#define mccompcurtype  L_monitor
#define mccompcurindex 12
#define nL mccdet_lmon_nL
#define L_N mccdet_lmon_L_N
#define L_p mccdet_lmon_L_p
#define L_p2 mccdet_lmon_L_p2
{   /* Declarations of det_lmon=L_monitor() SETTING parameters. */
char* filename = mccdet_lmon_filename;
MCNUM xmin = mccdet_lmon_xmin;
MCNUM xmax = mccdet_lmon_xmax;
MCNUM ymin = mccdet_lmon_ymin;
MCNUM ymax = mccdet_lmon_ymax;
MCNUM xwidth = mccdet_lmon_xwidth;
MCNUM yheight = mccdet_lmon_yheight;
MCNUM Lmin = mccdet_lmon_Lmin;
MCNUM Lmax = mccdet_lmon_Lmax;
MCNUM restore_neutron = mccdet_lmon_restore_neutron;
int nowritefile = mccdet_lmon_nowritefile;
#line 120 "/usr/share/mcstas/2.5/monitors/L_monitor.comp"
{
  
  multiline(5, (double)xmin, (double)ymin, 0.0,
               (double)xmax, (double)ymin, 0.0,
               (double)xmax, (double)ymax, 0.0,
               (double)xmin, (double)ymax, 0.0,
               (double)xmin, (double)ymin, 0.0);
}
#line 16874 "ISIS_MERLIN.c"
}   /* End of det_lmon=L_monitor() SETTING parameter declarations. */
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
#define mccompcurtype  Incoherent
#define mccompcurindex 13
#define VarsInc mccsample_VarsInc
#define offdata mccsample_offdata
{   /* Declarations of sample=Incoherent() SETTING parameters. */
char* geometry = mccsample_geometry;
MCNUM radius = mccsample_radius;
MCNUM xwidth = mccsample_xwidth;
MCNUM yheight = mccsample_yheight;
MCNUM zdepth = mccsample_zdepth;
MCNUM thickness = mccsample_thickness;
MCNUM target_x = mccsample_target_x;
MCNUM target_y = mccsample_target_y;
MCNUM target_z = mccsample_target_z;
MCNUM focus_r = mccsample_focus_r;
MCNUM focus_xw = mccsample_focus_xw;
MCNUM focus_yh = mccsample_focus_yh;
MCNUM focus_aw = mccsample_focus_aw;
MCNUM focus_ah = mccsample_focus_ah;
int target_index = mccsample_target_index;
MCNUM pack = mccsample_pack;
MCNUM p_interact = mccsample_p_interact;
MCNUM f_QE = mccsample_f_QE;
MCNUM gamma = mccsample_gamma;
MCNUM Etrans = mccsample_Etrans;
MCNUM deltaE = mccsample_deltaE;
MCNUM sigma_abs = mccsample_sigma_abs;
MCNUM sigma_inc = mccsample_sigma_inc;
MCNUM Vc = mccsample_Vc;
MCNUM concentric = mccsample_concentric;
MCNUM order = mccsample_order;
#line 438 "/usr/share/mcstas/2.5/samples/Incoherent.comp"
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
#line 16936 "ISIS_MERLIN.c"
}   /* End of sample=Incoherent() SETTING parameter declarations. */
#undef offdata
#undef VarsInc
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'detectorarm'. */
  SIG_MESSAGE("detectorarm (McDisplay)");
  printf("MCDISPLAY: component %s\n", "detectorarm");
#define mccompcurname  detectorarm
#define mccompcurtype  Arm
#define mccompcurindex 14
#line 40 "/usr/share/mcstas/2.5/optics/Arm.comp"
{
  /* A bit ugly; hard-coded dimensions. */
  
  line(0,0,0,0.2,0,0);
  line(0,0,0,0,0.2,0);
  line(0,0,0,0,0,0.2);
}
#line 16958 "ISIS_MERLIN.c"
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'detector'. */
  SIG_MESSAGE("detector (McDisplay)");
  printf("MCDISPLAY: component %s\n", "detector");
#define mccompcurname  detector
#define mccompcurtype  TOF_monitor
#define mccompcurindex 15
#define nt mccdetector_nt
#define TOF_N mccdetector_TOF_N
#define TOF_p mccdetector_TOF_p
#define TOF_p2 mccdetector_TOF_p2
#define t_min mccdetector_t_min
#define t_max mccdetector_t_max
#define delta_t mccdetector_delta_t
{   /* Declarations of detector=TOF_monitor() SETTING parameters. */
char* filename = mccdetector_filename;
MCNUM xmin = mccdetector_xmin;
MCNUM xmax = mccdetector_xmax;
MCNUM ymin = mccdetector_ymin;
MCNUM ymax = mccdetector_ymax;
MCNUM xwidth = mccdetector_xwidth;
MCNUM yheight = mccdetector_yheight;
MCNUM tmin = mccdetector_tmin;
MCNUM tmax = mccdetector_tmax;
MCNUM dt = mccdetector_dt;
MCNUM restore_neutron = mccdetector_restore_neutron;
int nowritefile = mccdetector_nowritefile;
#line 128 "/usr/share/mcstas/2.5/monitors/TOF_monitor.comp"
{
  
  multiline(5, (double)xmin, (double)ymin, 0.0,
               (double)xmax, (double)ymin, 0.0,
               (double)xmax, (double)ymax, 0.0,
               (double)xmin, (double)ymax, 0.0,
               (double)xmin, (double)ymin, 0.0);
}
#line 16998 "ISIS_MERLIN.c"
}   /* End of detector=TOF_monitor() SETTING parameter declarations. */
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
/* end of generated C code ISIS_MERLIN.c */
